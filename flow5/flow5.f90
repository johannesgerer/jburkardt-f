program main

!*****************************************************************************80
!
!! MAIN is the main program for FLOW5.
!
!  Discussion:
!
!    FLOW5 solves a steady state incompressible flow problem using finite elements.
!
!    FLOW5 is a pared down and simplified version of the research code FLOW4.  
!
!    FLOW4 includes parameterization, sensitivity analysis, the use of 
!    isoparametric elements, a curved bump on the bottom of the channel region,
!    more choices for quadrature points and other features which have 
!    been removed from FLOW5.
!
!    FLOW5 and FLOW4 were developed by John Burkardt, based on programs
!    of Janet Peterson and Max Gunzburger.
!
!    FLOW5 writes an output file suitable for use with the TECPLOT
!    graphics program.
!
!
!    The program works on an underlying fluid flow problem, whose
!    behavior is determined by a particular version of the Navier Stokes
!    equations.
!
!    The fluid flow in the region is described by three state functions:
!
!      U(X,Y), the horizontal velocity,
!      V(X,Y), the vertical velocity, 
!      P(X,Y), the pressure.  
!
!    In theory, these functions may be determined once we know the partial 
!    differential equations that govern them within the region, and
!    the value of the functions or certain derivatives of them along the 
!    boundary of the region.
!
!    For our work, we assume that at every point within the flow region, 
!    the functions obey the Navier Stokes equations for stationary, 
!    incompressible, viscous flow:
!
!    - nu * (ddU/dxdx + ddU/dydy) + U * dU/dx + V * dU/dy + dP/dx  = 0
!
!    - nu * (ddV/dxdx + ddV/dydy) + U * dV/dx + V * dV/dy + dP/dy  = 0
!
!    dU/dx + dV/dy = 0
!
!    Here, nu is a physical parameter called the "dynamic viscosity,"
!
!    We prefer the equivalent formulation (when nu is nonzero):
!
!    - (ddU/dxdx + ddU/dydy) + inv_nu * (U*dU/dx + V*dU/dy + dP/dx)  = 0
!
!    - (ddV/dxdx + ddV/dydy) + inv_nu * (U*dV/dx + V*dV/dy + dP/dy)  = 0
!
!    dU/dx + dV/dy = 0
!
!    where inv_nu=1/nu.
!
!
!    To complete the specification of the problem, we specify boundary 
!    conditions for the flow functions.
!
!    There are two problems built in to the code, a channel flow, and
!    a driven cavity.  
!
!    The boundary conditions for the channel flow are:
!
!      The values of U and V are specified along the left boundary;
!      U and V must be zero along the upper and lower walls;
!      dU/dn must be zero, and V must be zero, at the outflow;
!      P must be zero at a single point on the boundary.
!
!    The boundary conditions for the driven cavity are:
!
!      The values of U and V are specified along the top boundary;
!      U and V must be zero along the left, right, and bottom walls;
!      P must be zero at a single point on the boundary.
!
!*****************************************************************************80
!
!  DERIVATION OF FINITE ELEMENT EQUATIONS
!
!  Except for special cases, such as the Poiseuille flow solution 
!  discussed elsewhere, there are no methods of producing the exact 
!  solution functions U, V and P for a general Navier Stokes problem.  
!  In order to get any insight into flow problems, we must replace the 
!  original problem by one that is much weaker.  It's important that the 
!  weaker problem be possible to solve, and that the solutions produced 
!  are in general close to solutions of the original problem, and that 
!  these solutions can be made even closer,if desired.
!
!  A standard method of doing this is to use the method of finite 
!  elements.
!
!  To do so, we assume that instead of being smooth but otherwise 
!  completely arbitrary functions, that U, V and P are representable 
!  as linear combinations of a finite set of basis functions.
!
!  We multiply the first two equations by an arbitrary velocity basis
!  function Wi, and the third equation by an arbitrary pressure basis
!  function Qi, and integrate over the region.  The integrand of the
!  resulting finite element equations is then transformed, using
!  integration by parts, into:
!
!    U-Eqn(I):
!
!      (dU/dx*dWi/dx + dU/dy*dWi/dy) 
!        + Re*(U*dU/dx + V*dU/dy + dP/dx ) * Wi
!
!    V-Eqn(I):
!
!      (dV/dx*dWi/dx + dV/dy*dWi/dy) 
!        + Re*(U*dV/dx + V*dV/dy + dP/dy ) * Wi
!
!    P-Eqn(I):
!
!      (dU/dx + dV/dy) * Qi
!
!
!  These integrands may be rewritten using the program's variable names:
!
!
!    dUdx*dwidx + dUdy*dwidy + inv_nu*(U*dUdx+V*dUdy+dPdx) * wi
!
!    dVdx*dwidx + dVdy*dwidy + inv_nu*(U*dVdx+V*dVdy+dPdy) * wi
!
!    ( dUdx + dVdy ) * qi
!
!
!  This system of nonlinear equations is then solved by Newton's method.
!  That means that we have to differentiate each nonlinear equation
!  with respect to the unknowns, getting the Jacobian matrix, and
!  solving DF(X) * DEL(X) = -F(X).  If we abuse notation, we can
!  consider the linear system DF(X) * DEL(X).
!
!  Here, variables U, V and P in capital letters are to be solved for, 
!  but the same variable names in lowercase represent the current
!  values of those same variables.
!
!
!  d U-Eqn(I) / d U-coefficient * U coefficient:
!
!    dUdx*dwidx + dUdy*dwidy + inv_nu*(U*dudx+u*dUdx+v*dUdy)*wi
!
!  d U-Eqn(I) / d V coefficient * V coefficient:
!
!    inv_nu*V*dudy*wi
!
!  d U-Eqn(I) / d P coefficient * P coefficient:
!
!    inv_nu*dPdx*wi
!
!  d V-Eqn(I) / d U coefficient * U coefficient:
!
!    inv_nu*U*dvdx*wi
!
!  d V-Eqn(I) / d V coefficient * V coefficient:
!
!    dVdx*dwidx + dVdy*dwidy + inv_nu*(u*dVdx+v*dVdy+V*dvdy)*wi
!
!  d V-Eqn(I) / d P coefficient * P coefficient:
!
!    inv_nu*dPdy*wi
!
!  d P-Eqn(I) / d U coefficient * U coefficient:
!
!    dUdx * qi
!
!  d P-Eqn(I) / d V coefficient * V coefficient:
!
!    dVdx * qi
!
!
!*****************************************************************************80
!
!  SENSITIVITIES
!
!
!  This version of the program does not include any sensitivity 
!  calculations, but they would be easy to add.
!
!
!  Suppose that U, V and P depend in some way on a parameter
!  Z, and let us consider differentiating each of the three above
!  equations with respect to Z.  Then we interchange differentiation
!  where desired, and come up with equations for the SENSITIVITIES.
!
!  Now the sensitivities should be written as (dU/dZ, dV/dZ, dP/dZ).
!  In the ensuing equations, we will write them as (U, V, P), but
!  now the lower case letters (u, v, p) represent the current values
!  of the original fluid flow quantities.
!
!  As a simple example, consider the inverse dynamic viscosity to be a 
!  parameter, and regard the original flow variables as dependent on INV_NU.
!  Then we can compute the sensitivity equations by differentiating
!  the Navier Stokes equations with respect to INV_NU:
!
!  d U-Eqn(I) / d inv_nu:
!
!    dUdx*dwidx + dUdy*dwidy
!      + inv_nu*(U*dudx+u*dUdx+V*dudy+v*dUdy+dPdx)*wi
!      +        (u*dudx+v*dudy+dpdx)*wi = 0
!
!  d V-Eqn(I) / d inv_nu:
!
!    dVdx*dwidx + dVdy*dwidy
!      + inv_nu*(U*dvdx+u*dVdx+V*dvdy+v*dVdy+dPdy)*wi
!      +        (u*dvdx+v*dvdy+dpdy)*wi = 0
!
!  d P-Eqn(I) / d inv_nu:
!
!    (dUdx + dVdy) * qi = 0
!
!  Boundary conditions:
!
!    0 everywhere.
!
!  In the case of the INV_NU parameter, we carry the "extra" terms
!  (u*dudx+v*dudy+dpdx)*wi and (u*dvdx+v*dvdy+dpdy)*wi.
!  to the right hand side, and treat them as source terms.  
!
!
!*****************************************************************************80
!
!  POISEUILLE FLOW
!
!
!  Consider a horizontal channel of constant height H, and of length L.
!
!  Suppose a parabolic inflow is specified at the left hand opening,
!  where X=0, of the form
!
!    U(0,Y) = S * Y * (H-Y)
!    V(0,Y) = 0
!    P(0,Y) = 0
!
!  where S is any value.
!
!  Then the following functions (U,V,P) solve the Navier Stokes
!  equations in the region:
!
!    U(X,Y) = S * Y * (H-Y)
!    V(X,Y) = 0
!    P(X,Y) = -2 * S * X * RE
!
!  The standard problem we use has H = 3, L = 10, and chooses 
!  S = (4/9)*Lambda, so that the maximum value of the parabolic inflow, 
!  at Y = H/2, is Lambda.  Then our formula becomes:
!
!    U(X,Y) = (4/9) * Lambda * Y * (3-Y)
!    V(X,Y) = 0
!    P(X,Y) = -2 * (4/9) * Lambda * X * RE
!
!  (In the full application, we include a bump on the lower wall of
!  the channel.  The presence of this bump disrupts the flow, and
!  the formula above is no longer valid.)
!
!*****************************************************************************80
!
!  DRIVEN CAVITY FLOW
!
!  The driven cavity flow is more interesting than Poiseuille flow,
!  but we do not have a closed form solution for it.
!
!  The flow region is, we will say, a square box of unit side, with
!  walls on the sides and bottom, and open on the top.  A tangential
!  force is applied uniformly along the top, which we will enforce
!  with a boundary condition of the form:
!
!    U(X,1) = 1
!    V(X,1) = 0
!
!*****************************************************************************80
!
!  GEOMETRIC INFORMATION
!
!  The following information concerns programming details for this
!  implementation of the finite element method.  It is intended to help
!  the reader understand, test, and modify the code.
!
!
!  1) The finite element nodes
!
!  If the region is rectangular, then FLOW5 places the nodes in such a 
!  way that they are evenly spaced in the X direction, and in the Y 
!  direction, although these two spacings may be different.
!
!  The first node is in the lower left corner.  The second node is the 
!  one immediately above the first, and then numbering proceeds 
!  upwards, and then over to the next column.  For instance:
!
!  Y=3.00       13          26          39          42          65
!  Y=2.75       12          25          38          41          64
!  Y=2.50       11          24          37          50          63
!  Y=2.25       10          23          36          49          62
!  Y=2.00        9          22          35          48          61
!  Y=1.75        8          21          34          47          60
!  Y=1.50        7          20          33          46          59
!  Y=1.25        6          19          32          45          58
!  Y=1.00        5          18          31          44          57
!  Y=0.75        4          17          30          43          56
!  Y=0.50        3          16          29          42          55
!  Y=0.25        2          15          28          41          54
!  Y=0.00        1          14          27          40          53
!
!            X=0.00      X=0.25      X=0.50      X=0.75      X=1.00
!
!
!  2) The basic elements
!
!
!   2--5--3          2
!   |    /          /|
!   |   /          / |
!   4  6          4  5
!   | /          /   |
!   |/          /    |
!   1          1--6--3
!
!
!  3) The quadrature points
!
!
!    For NQUAD = 3:
!
!   .--2--.          .
!   |    /          /|
!   |   /          / |
!   1  3          1  2
!   | /          /   |
!   |/          /    |
!   .          .--3--.
!
!
!
!  4) The elements in the grid
!
!  Here is a schematic of the 24 elements defined by the nodes shown
!  in the earlier diagram:
!
!
!             13--26--39--42--65
!              | 11  / | 23  / |
!              |    /  |    /  |
!             12  25  38  41  64
!              | /     | /     |
!              |/   12 |/   24 |
!             11--24--37--50--63
!              |  9  / | 21  / |
!              |    /  |    /  |
!             10  23  36  49  62
!              | /     | /     |
!              |/   10 |/   22 |
!              9--22--35--48--61
!              |  7  / | 19  / |
!              |    /  |    /  |
!              8  21  34  47  60
!              | /     | /     |
!              |/    8 |/   20 |
!              7--20--33--46--59
!              |  5  / | 17  / |
!              |    /  |    /  |
!              6  19  32  45  58
!              | /     | /     |
!              |/    6 |/   18 |
!              5--18--31--44--57
!              |  3  / | 15  / |
!              |    /  |    /  |
!              4  17  30  43  56
!              | /     | /     |
!              |/    4 |/   16 |
!              3--16--29--42--55
!              |  1  / | 13  / |
!              |    /  |    /  |
!              2  15  28  41  54
!              | /     | /     |
!              |/    2 |/   14 |
!              1--14--27--40--53
!
!
!  5) Numbering for a sample problem.
!
!  Here is how the first 92 unknowns would be numbered, for a channel 
!  problem, with NY=7 and NX=21.  
!
!
!  Y=3.00    U31 V32 P33   U58 V59     U90 V91 P92
!  Y=2.75    U29 V30       U56 V29     U88 V89      
!  Y=2.50    U26 V27 P28   U54 V27     U85 V86 P87  
!  Y=2.25    U24 V25       U52 V25     U83 V84      
!  Y=2.00    U21 V22 P23   U50 V23     U80 V81 P82  
!  Y=1.75    U19 V20       U48 V21     U78 V79      
!  Y=1.50    U16 V17 P18   U46 V19     U75 V76 P77  
!  Y=1.25    U14 V15       U44 V17     U73 V74      
!  Y=1.00    U11 V12 P13   U42 V15     U70 V71 P72  
!  Y=0.75    U09 V10       U40 V41     U68 V69      
!  Y=0.50    U06 V07 P08   U38 V39     U65 V66 P67  
!  Y=0.25    U04 V05       U36 V37     U63 V64      
!  Y=0.00    U01 V02 P03   U34 V35     U60 V61 P62
!
!            X=0.00        X=0.25      X=0.50       ...
!
!
!*****************************************************************************80
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    real ( kind = 8 ) A(MAXROW,MAXEQN), contains a matrix in LINPACK general 
!    band storage mode.  The two dimensional array is of logical dimensions 
!    NROW by NEQN.
!
!    real ( kind = 8 ) AREA(MAXQUAD,MAXELM), contains a common factor 
!    multiplying the term associated with a quadrature point in a given 
!    element, namely,
!      AREA(IQUAD,IELEM) = Ar(IELEM) * WQUAD(IQUAD)
!    or, if the element is isoparametric,
!      AREA(IQUAD,IELEM) = DET * Ar(IELEM) * WQUAD(IQUAD)
!    Here Ar(IELEM) represents the area of element IELEM.
!
!    character ( len = 2 ) EQN(MAXEQN).
!    EQN records the "type" of each equation that will be generated, 
!    and which is associated with an unknown.  Note that most 
!    boundary conditions do not result in an equation.  The current 
!    values are:
!    'U'  The horizontal momentum equation.
!    'UB' The condition U=0 applied at a node on the bump.
!    'UI' The condition U=UInflow(Y,Lambda) at the inflow.
!    'UW' The condition U=0 applied at a node on a fixed wall.
!    'V'  The vertical momentum equation.
!    'VB' The condition V=0 applied at a node on the bump.
!    'VI' The condition V=VInflow(Y,Lambda) at the inflow.
!    'VW' The condition V=0 applied at a node on a fixed wall.
!    'P'  The continuity equation.
!    'PB' The condition P=0 applied at (XMAX,YMAX).
!
!    real ( kind = 8 ) ETAQUAD(MAXQUAD), the "Eta" quadrature coordinates.
!
!    real ( kind = 8 ) G(MAXEQN), G is the current solution vector, in which 
!    are stored the finite element coefficients that define the velocity
!    and pressure functions, U, V and P.
!
!    integer IERROR, an error flag.
!    0, no error occurred in this routine.
!    nonzero, an error occurred.
!
!    integer INDX(3,MAXNP).
!    INDX(I,J) contains, for each node J, the index of U, V and P 
!    at that node, or 0 or a negative value.
!    If K=INDX(I,J) is positive, then the value of the degree
!    of freedom is stored in the solution vector entry G(K),
!    and an equation will be generated to determine its value.
!    If INDX(I,J) is not positive, then no equation is
!    generated to determine for variable I at node J, either 
!    because the variable is specified in some other way, or 
!    because (in the case of pressure), there is no coefficient 
!    associated with that node.
!
!    integer IPIVOT(MAXEQN), pivoting space for the linear system solver.
!
!    integer ISOTRI(MAXELM), 0/1, the element IS NOT/IS isoparametric.
!
!    integer MAXELM, the maximum number of elements.
!
!    integer MAXEQN, the maximum number of equations.
!
!    integer MAXNEW, the maximum number of Newton iterations.
!
!    integer MAXNP, the maximum number of nodes.
!
!    integer MAXQUAD, the maximum number of quadrature points.
!
!    integer MAXROW, the first dimension of the matrix A.
!
!    integer NELEM, the number of elements.
!
!    integer NEQN, the number of finite element equations.
!
!    integer NLBAND, the lower bandwidth of the matrix A.
!
!    integer NODE(6,MAXELM), the nodes that make up each element.
!    The local ordering of the nodes is suggested by this diagram:
!
!          2
!         /|
!        4 5
!       /  |
!      1-6-3
!
!    integer NP, the number of nodes.  NP=(2*NX-1)*(2*NY-1).
!
!    integer NROW, the number of rows needed to store the matrix A.
!
!    integer NX, controls the spacing of nodes and elements in
!    the X direction.  There are 2*NX-1 nodes along various
!    lines in the X direction.
!    Roughly speaking, NX (or 2*NX) is the number of elements along
!    a line in the X direction.
!
!    integer NY, controls the spacing of nodes and elements in
!    the Y direction.  There are 2*NY-1 nodes along various
!    lines in the Y direction.
!    Roughly speaking, NY (or 2*NY) is the number of elements along
!    a line in the Y direction.
!
!    real ( kind = 8 ) PHI(MAXQUAD,6,10,MAXELM).
!    PHI contains the value of a basis function, its derivative,
!    or other information, evaluated at a quadrature point.
!    For a particular element I, quadrature point J, and basis
!    function K, we use the following shorthand for the 
!    entries of PHI:
!      W, dWdX, dWdY
!      Q, dQdX, dQdY
!      dXsidX, dXsidY, dEtadX, dEtadY
!    W is the quadratic basis function associated with velocity,
!    Q the linear basis function associated with pressure,
!    Xsi and Eta the reference coordinates for the point.
!    In particular, PHI(J,K,1,I) is the value of the quadratic
!    basis function associated with local node K in element I,
!    evaluated at quadrature point J.
!    Note that PHI(J,K,4,I)=PHI(J,K,5,I)=PHI(J,K,6,I)=0 for
!    K = 4, 5, or 6, since there are only three linear basis
!    functions.
!
!    real ( kind = 8 ) PNORM, the maximum absolute value of the pressure.
!
!    character ( len = 20 ) REGION, specifies the flow region, 
!    'CHANNEL' or 'CAVITY'.
!
!    real ( kind = 8 ) RES(MAXEQN), the residual.
!
!    real ( kind = 8 ) INV_NU, the value of the inverse dynamic viscosity.
!
!    real ( kind = 8 ) TOLNEW, the convergence tolerance for the Newton
!    iteration.
!
!    real ( kind = 8 ) UVNORM, the maximum velocity magnitude.
!
!    real ( kind = 8 ) WQUAD(MAXQUAD), the quadrature weights.
!
!    real ( kind = 8 ) XC(MAXNP), the X coordinates of the nodes.
!
!    real ( kind = 8 ) XQUAD(MAXQUAD,MAXELM), the X quadrature coordinates.
!
!    real ( kind = 8 ) XRANGE, the width of the region.
!
!    real ( kind = 8 ) XSIQUAD(MAXQUAD), the "Xsi" quadrature coordinates.
!
!    real ( kind = 8 ) YC(MAXNP), the Y coordinates of the nodes.
!
!    real ( kind = 8 ) YQUAD(MAXQUAD,MAXELM), the Y quadrature coordinates.
!
!    real ( kind = 8 ) YRANGE, the height of the region.
!
  implicit none
!
!  NX and NY can't be any larger than MAXNX and MAXNY.  You can
!  increase them here.
!
  integer, parameter :: maxnx = 41
  integer, parameter :: maxny = 13
  integer, parameter :: maxquad = 3
!
!  The following quantities depend on the values of MAXNX and MAXNY.
!
!  The assignment of MAXROW assumes that the nodes are ordered starting
!  at the bottom left corner, proceeding upwards in a column, and then
!  moving back to the bottom of the next column to the right.  For
!  our choice of elements, this allows us to estimate the greatest
!  difference between indices of two variables which occur in the
!  same equation.  
!
!  If our estimate is wrong, then SETBAN will catch it.
!
  integer, parameter :: maxrow = 29*maxny
  integer, parameter :: maxelm = 2*(maxnx-1)*(maxny-1)
  integer, parameter :: maxeqn = 2*(2*maxnx-1)*(2*maxny-1)+maxnx*maxny
  integer, parameter :: maxnp = (2*maxnx-1)*(2*maxny-1)

  real ( kind = 8 ) a(maxrow,maxeqn)
  real ( kind = 8 ) area(maxquad,maxelm)
  character ( len = 2 ) eqn(maxeqn)
  real ( kind = 8 ) etaquad(maxquad)
  real ( kind = 8 ) g(maxeqn)
  integer i
  integer ierror
  integer indx(3,maxnp)
  integer ipivot(maxeqn)
  integer isotri(maxelm)
  integer j
  integer k
  integer l
  integer maxnew
  integer nelem
  integer neqn
  integer nlband
  integer node(6,maxelm)
  integer np
  integer nquad
  integer nrow
  integer nx
  integer ny
  real ( kind = 8 ) p(maxnp)
  real ( kind = 8 ) phi(maxquad,6,10,maxelm)
  real ( kind = 8 ) pnorm
  character ( len = 20 ) region
  real ( kind = 8 ) res(maxeqn)
  real ( kind = 8 ) inv_nu
  real ( kind = 8 ) tolnew
  real ( kind = 8 ) u(maxnp)
  real ( kind = 8 ) v(maxnp)
  real ( kind = 8 ) uvnorm
  real ( kind = 8 ) wquad(maxquad)
  real ( kind = 8 ) xc(maxnp)
  real ( kind = 8 ) xquad(maxquad,maxelm)
  real ( kind = 8 ) xrange
  real ( kind = 8 ) xsiquad(maxquad)
  real ( kind = 8 ) yc(maxnp)
  real ( kind = 8 ) yquad(3,maxelm)
  real ( kind = 8 ) yrange

  call timestamp ( )
!
!  Print the program name, date, and computer name.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FLOW5'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Last modified on 01 June 2000.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Maximum problem grid size :'
  write ( *, '(a,i8)' ) '    MAXNX = ', maxnx
  write ( *, '(a,i8)' ) '    MAXNY = ', maxny
  write ( *, '(a)' ) ' '
!
!  Initialize some of the variables.
!
  a(1:maxrow,1:maxeqn) = 0.0D+00
  area(1:maxquad,1:maxelm) = 0.0D+00
  eqn(1:maxeqn) = '??'  
  etaquad(1:maxquad) = 0.0D+00
  g(1:maxeqn) = 0.0D+00
  ierror = 0
  indx(1:3,1:maxnp) = 0
  ipivot(1:maxeqn) = 0
  isotri(1:maxelm) = 0
  maxnew = 0
  nelem = 0
  neqn = 0
  nlband = 0
  node(1:6,1:maxelm) = 0
  np = 0
  nquad = 0
  nrow = 0
  nx = 0
  ny = 0
  phi(1:maxquad,1:6,1:10,1:maxelm) = 0.0D+00
  pnorm = 0.0D+00
  region = '??'
  res(1:maxeqn) = 0.0D+00
  inv_nu = 0.0D+00
  tolnew = 0.0D+00
  uvnorm = 0.0D+00
  wquad(1:maxquad) = 0.0D+00 
  xc(1:maxnp) = 0.0D+00
  xquad(1:maxquad,1:maxelm) = 0.0D+00
  xrange = 0.0D+00 
  xsiquad(1:maxquad) = 0.0D+00
  yc(1:maxnp) = 0.0D+00
  yquad(1:maxquad,1:maxelm) = 0.0D+00
  yrange = 0.0D+00
!
!  Now we set some values which will define all the others.
!
  maxnew = 10
  nquad = 3

  nx = 11
  if ( maxnx < nx ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FLOW5 - Fatal error!'
    write ( *, '(a)' ) '  MAXNX < NX.'
    write ( *, '(a,i8)' ) '  NX = ', nx
    write ( *, '(a,i8)' ) '  MAXNX = ', maxnx
    stop 
  end if
     
  ny = 11
  if ( maxny < ny ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FLOW5 - Fatal error!'
    write ( *, '(a)' ) '  MAXNY < NY.'
    write ( *, '(a,i8)' ) '  NY = ', ny
    write ( *, '(a,i8)' ) '  MAXNY = ', maxny
    stop 
  end if
  
  region = 'CAVITY'
  inv_nu = 1.0D+00      
  tolnew = 0.0000000001D+00
!
!  Set the total width and height of the region.
!
  if ( region == 'CHANNEL' ) then
    xrange = 10.0D+00
    yrange = 3.0D+00
  else if ( region == 'CAVITY' ) then
    xrange = 1.0D+00
    yrange = 1.0D+00
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FLOW5 - Fatal error!'
    write ( *, '(a)' ) '  An unexpected region was specified.'
    write ( *, '(a)' ) '  REGION = ' // trim ( region )
    stop
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FLOW5:'
  write ( *, '(a)' ) '  The flow region is ' // trim ( region )
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of horizontal elements, NX =', nx
  write ( *, '(a,i8)' ) '  Number of vertical elements, NY =  ', ny
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) &
    '  Maximum number of Newton iterations, MAXNEW =', maxnew
  write ( *, '(a,g14.6)' ) '  Newton iteration tolerance, TOLNEW =', tolnew
!
!  Now compute NELEM, the number of elements and NP, the number of 
!  nodes, based on the user choices for NX and NY.
!
  nelem = 2 * (nx-1) * (ny-1)
  np = (2*nx-1) * (2*ny-1)    

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FLOW5:'
  write ( *, '(a,i8)' ) '  The number of elements, NELEM =', nelem
  write ( *, '(a,i8)' ) '  The number of nodes, NP =', np
!
!  Set up the element information in EQN, INDX, ISOTRI, NEQN, and NODE.
!
  call setnod ( eqn, indx, isotri, maxelm, maxeqn, maxnp,  &
    neqn, node, np, nx, ny, region )  

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FLOW5:'
  write ( *, '(a,i8)' ) '  The number of unknowns is NEQN = ', neqn
!
!  Figure out NLBAND, the lower matrix bandwidth, and NROW,
!  the number of rows we will use in the matrix A.
!
  call setban ( indx, maxelm, maxnp, maxrow, nelem, nlband, node, nrow )
  
  write ( *, '(a)' ) ' ' 
  write ( *, '(a)' ) 'FLOW5:'
  write ( *, '(a,i8)' ) '  Lower bandwidth NLBAND =      ', nlband
  write ( *, '(a,i8)' ) '  Required matrix rows NROW =   ', nrow
  write ( *, '(a,i8)' ) '  Maximum matrix rows, MAXROW = ', maxrow
!
!  Set XC and YC, the coordinates of the nodes.
!
  call setxy ( maxnp, nx, ny, xc, xrange, yc, yrange )
!
!  Set XQUAD, YQUAD, the locations of the quadrature points,
!  WQUAD, the quadrature weights, and AREA, the element areas.
!
  call setquad ( area, etaquad, isotri, maxelm, maxnp, maxquad, &
    nelem, node, nquad, wquad, xc, xquad, xsiquad, yc, yquad )
!
!  Set the value of the basis functions at all quadrature points.
!
  call setbas ( area, etaquad, isotri, maxelm, maxnp, maxquad, &
    nelem, node, nquad, phi, xc, xquad, xsiquad, yc, yquad )
!
!  Solve the nonlinear system.
!
  call newton ( a, area, eqn, g, ierror, indx, ipivot, maxelm, &
    maxeqn, maxnew, maxnp, maxquad, nelem, neqn, nlband, node, &
    np, nquad, nrow, phi, region, res, inv_nu, tolnew, yc, yrange )
 
  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FLOW5 - Fatal error!'
    write ( *, '(a)' ) '  The Newton iteration failed!'
    stop
  else
    call uvpnorm ( g, indx, maxeqn, maxnp, np, pnorm, uvnorm )
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FLOW5:'
    write ( *, '(a)' ) '  NEWTON computed the Navier Stokes solution.'
    write ( *, '(a,g14.6)' ) '  Norm of maximum velocity, UVNORM =', uvnorm
    write ( *, '(a,g14.6)' ) '  Norm of maximum pressure, PNORM = ', pnorm
  end if
!
!  Interpolate values of the pressure at nodes which do not
!  have an associated pressure unknown, and create the P array.
!
  call intprs ( g, indx, maxelm, maxeqn, maxnp, nelem, node, p )
!
!  Copy the U and V arrays out of G.
!
  do i = 1, np
    u(i) = g ( indx(1,i) )
    v(i) = g ( indx(2,i) )
  end do
!
!  Write U, V and P out to a TECPLOT file.
!
  call write_tecplot_file ( maxelm, maxnp, nelem, node, np, p, u, v, xc, yc )
  
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FLOW5:'
  write ( *, '(a)' ) '  Wrote the TECPLOT data file.'
!
!  Writes out the data in a very simple format.
!
  call write_display_file ( ierror, maxelm, maxnp, nelem, node, np, p, &
    u, v, xc, yc )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FLOW5:'
  write ( *, '(a)' ) '  Wrote the DISPLAY data files.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FLOW5:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine bsp ( q, dqdx, dqdy, ielem, iq, maxelm, maxnp, node, xc, xq, yc, yq )

!*****************************************************************************80
!
!! BSP evaluates the linear basis functions for pressure.
!
!  Discussion:
!
!    Here is a typical finite element associated with pressure:
!
!        2
!       /|
!      / |
!     /  |
!    1---3
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) Q, the value of the IQ-th basis at (XQ,YQ).
!
!    Output, real ( kind = 8 ) DQDX, DQDY, the X and Y derivatives of the 
!    IQ-th basis function at (XQ,YQ).
!
!    Input, integer IELEM, the index of the element.
!
!    Input, integer IQ, the index of the basis function.
!
!    Input, integer NELEM, the number of elements.
!
!    Input, integer NODE(6,MAXELM), the nodes that make up each element.
!
!    Input, integer NP, the number of nodes.
!
!    Input, real ( kind = 8 ) XC(MAXNP), the X coordinates of the nodes.
!
!    Input, real ( kind = 8 ) XQ, the X coordinate of the evaluation point.
!
!    Input, real ( kind = 8 ) YC(MAXNP), the Y coordinates of the nodes.
!
!    Input, real ( kind = 8 ) YQ, the Y coordinate of the evaluation point.
!
  implicit none

  integer maxelm
  integer maxnp

  real ( kind = 8 ) q
  real ( kind = 8 ) dqdx
  real ( kind = 8 ) dqdy
  real ( kind = 8 ) d
  integer i1
  integer i2
  integer i3
  integer ielem
  integer iq
  integer iq1
  integer iq2
  integer iq3
  integer node(6,maxelm)
  real ( kind = 8 ) xc(maxnp)
  real ( kind = 8 ) xq
  real ( kind = 8 ) yc(maxnp)
  real ( kind = 8 ) yq

  if ( iq < 1 .or. 6 < iq ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BSP - Fatal error!'
    write ( *, '(a,i8)' ) '  The requested basis function is IQ =', iq
    write ( *, '(a)' ) '  but only values from 1 to 6 are legal.'
    stop
  else if ( 4 <= iq .and. iq <= 6 ) then
    q = 0.0D+00
    dqdx = 0.0D+00
    dqdy = 0.0D+00
    return
  end if
 
  iq1 = iq
  iq2 = mod ( iq, 3 ) + 1
  iq3 = mod ( iq+1, 3 ) + 1
 
  i1 = node(iq1,ielem)
  i2 = node(iq2,ielem)
  i3 = node(iq3,ielem)
 
  d = ( xc(i2) - xc(i1) ) * ( yc(i3) - yc(i1) ) &
     -( xc(i3) - xc(i1) ) * ( yc(i2) - yc(i1) )
 
  dqdx = ( yc(i2) - yc(i3) ) / d
  dqdy = ( xc(i3) - xc(i2) ) / d
 
  q = 1.0D+00 + dqdx * ( xq - xc(i1) ) + dqdy * ( yq - yc(i1) )
 
  return
end
subroutine dgb_fa ( a, lda, n, ml, mu, ipivot, info )

!*****************************************************************************80
!
!! DGB_FA factors a matrix stored in LINPACK general band storage.
!
!  Discussion:
!
!    The matrix is stored in the array using LINPACK general band storage.
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
!    The total number of rows needed in A is 2*ML+MU+1.
!    The ML+MU by ML+MU upper left triangle and the
!    ML by ML lower right triangle are not referenced.
!
!  Modified:
!
!    02 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) A(LDA,N), the matrix in band storage.  The
!    columns of the matrix are stored in the columns of the array,
!    and the diagonals of the matrix are stored in rows ML+1 through
!    2*ML+MU+1.  On return, A has been overwritten by the LU factors.
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least 2*ML+MU+1.
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N-1.
!
!    Output, integer IPIVOT(N), the pivot vector.
!
!    Output, integer INFO, singularity flag.
!    0, no singularity detected.
!    nonzero, the factorization failed on the INFO-th step.
!
  implicit none

  integer lda
  integer n

  real ( kind = 8 ) a(lda,n)
  integer i
  integer i0
  integer ierror
  integer info
  integer ipivot(n)
  integer j
  integer j0
  integer j1
  integer ju
  integer jz
  integer k
  integer l
  integer lm
  integer m
  integer ml
  integer mm
  integer mu
  real t

  m = ml + mu + 1
  info = 0
!
!  Zero out the initial fill-in columns.
!
  j0 = mu + 2
  j1 = min ( n, m ) - 1

  do jz = j0, j1
    i0 = m + 1 - jz
    do i = i0, ml
      a(i,jz) = 0.0D+00
    end do
  end do

  jz = j1
  ju = 0

  do k = 1, n-1
!
!  Zero out the next fill-in column.
!
    jz = jz + 1
    if ( jz <= n ) then
      do i = 1, ml
        a(i,jz) = 0.0D+00
      end do
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

    ipivot(k) = l + k - m
!
!  Zero pivot implies this column already triangularized.
!
    if ( a(l,k) == 0.0D+00 ) then
      info = k
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DGB_FA - Warning!'
      write ( *, '(a,i8)' ) '  Zero pivot on step ', info
      return
    end if
!
!  Interchange if necessary.
!
    if ( l /= m ) then
      t = a(l,k)
      a(l,k) = a(m,k)
      a(m,k) = t
    end if
!
!  Compute multipliers.
!
    do j = m+1, m+lm
      a(j,k) = - a(j,k) / a(m,k)
    end do
!
!  Row elimination with column indexing.
!
    ju = max ( ju, mu+ipivot(k) )
    ju = min ( ju, n )
    mm = m

    do j = k+1, ju

      l = l - 1
      mm = mm - 1
      t = a(l,j)
      if ( l /= mm ) then
        a(l,j) = a(mm,j)
        a(mm,j) = t
      end if

      do i = 1, lm
        a(mm+i,j) = a(mm+i,j) + t * a(m+i,k)
      end do

    end do

  end do

  ipivot(n) = n
  if ( a(m,n) == 0.0D+00 ) then
    info = n
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DGB_FA - Warning!'
    write ( *, '(a,i8)' ) '  Zero pivot on step ', info
  end if

  return
end
subroutine dgb_sl ( a, lda, n, ml, mu, ipivot, b, job )

!*****************************************************************************80
!
!! DGB_SL solves a system factored by DGB_FA.
!
!  Modified:
!
!    02 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A(LDA,N), the LU factors from SGB_FA.
!
!    Input, integer LDA, the leading dimension of the array.
!    LDA must be at least 2*ML+MU+1.
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N-1.
!
!    Input, integer IPIVOT(N), the pivot vector from SGB_FA.
!
!    Input/output, real ( kind = 8 ) B(N).
!    On input, the right hand side vector.
!    On output, the solution.
!
!    Input, integer JOB.
!    0, solve A*X=B.
!    nonzero, solve transpose(A)*X=B.
!
  implicit none

  integer lda
  integer n

  real ( kind = 8 ) a(lda,n)
  real ( kind = 8 ) b(n)
  integer ierror
  integer ipivot(n)
  integer j
  integer job
  integer k
  integer l
  integer la
  integer lb
  integer lm
  integer m
  integer ml
  integer mu
  real t

  m = mu + ml + 1
!
!  Solve A * X = B.
!
  if ( job == 0 ) then
!
!  Solve L * Y = B.
!
    if ( 1 <= ml ) then

      do k = 1, n-1

        lm = min ( ml, n-k )
        l = ipivot(k)
        t = b(l)

        if ( l /= k ) then
          b(l) = b(k)
          b(k) = t
        end if

        do j = 1, lm
          b(k+j) = b(k+j) + t * a(m+j,k)
        end do

      end do
    end if
!
!  Solve U * X = Y.
!
    do k = n, 1, -1

      b(k) = b(k) / a(m,k)
      lm = min ( k, m ) - 1
      la = m - lm
      lb = k - lm
      t = -b(k)

      do j = 0, lm-1
        b(lb+j) = b(lb+j) + t * a(la+j,k)
      end do

    end do
!
!  Solve transpose(A) * X = B.
!
  else
!
!  Solve transpose(U) * Y = B.
!
    do k = 1, n
      lm = min ( k, m ) - 1
      la = m - lm
      lb = k - lm
      t = 0.0D+00
      do j = 0, lm-1
        t = t + a(la+j,k) * b(lb+j)
      end do
      b(k) = ( b(k) - t ) / a(m,k)
    end do
!
!  Solve transpose(L) * X = Y.
!
    if ( 1 <= ml ) then

      do k = n-1, 1, -1

        lm = min ( ml, n-k )

        t = 0.0D+00
        do j = 1, lm
          t = t + a(m+j,k) * b(k+j)
        end do

        b(k) = b(k) + t
        l = ipivot(k)

        if ( l /= k ) then
          t = b(l)
          b(l) = b(k)
          b(k) = t
        end if

      end do

    end if

  end if

  return
end
subroutine fp ( a, area, eqn, g, indx, maxelm, maxeqn, maxnp, maxquad, nelem, &
  neqn, nlband, node, np, nquad, nrow, phi, inv_nu )

!*****************************************************************************80
!
!! FP computes the jacobian matrix A of the Navier Stokes residual. 
!
!  Discussion:
!
!    Essentially, the routine FP is the "derivative" of the routine FX,
!    with respect to the coefficients G.
!
!
!    The differentiated Navier Stokes functions have the form:
!
!
!    d U-Eqn(I) / d U-Coef(J):
!
!      Integral dWj/dx * dWi/dx + dWj/dy * dWi/dy
!      + inv_nu * (Wj*dUold/dx + Uold*dWj/dx+ Vold*dWj/dy) * Wi dx dy
!
!    d U-Eqn(I) / d V-Coef(J):
!
!      Integral inv_nu * Wj * dUold/dy * Wi dx dy
!
!    d U-Eqn(I) / d P-Coef(J):
!
!      Integral inv_nu * dQj/dx * Wi dx dy
!
!    d V-Eqn(I) / d U-Coef(J):
!
!      Integral inv_nu * Wj*dVold/dx * Wi dx dy
!
!    d V-Eqn(I) / d V-Coef(J):
!
!      Integral dWj/dx * dWi/dx + dWj/dy * dWi/dy
!      + inv_nu * (Uold*dWj/dx + Wj*dVold/dy + Vold*dWj/dy) * Wi dx dy
!
!    d V-Eqn(I) / d P-Coef(J):
!
!      Integral inv_nu * dQj/dy * Wi dx dy
!
!    d P-Eqn(I) / d U-Coef(J):
!
!      Integral dWj/dx * Qi dx dy
!
!    d P-Eqn(I) / d V-Coef(J):
!
!      Integral dWj/dy * Qi dx dy
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A(NROW,MAXEQN), the
!    value of D F(I)/D X(J) for each of the NEQN residual
!    functions F(I) with respect to each of the unknown
!    coefficients X(J).
!
!    Input, real ( kind = 8 ) AREA(MAXQUAD,MAXELM).
!
!    AREA contains a common factor multiplying the term associated
!    with a quadrature point in a given element, namely,
!
!      AREA(IQUAD,IELEM) = Ar(IELEM) * WQUAD(IQUAD)
!
!    or, if the element is isoparametric,
!
!      AREA(IQUAD,IELEM) = DET * Ar(IELEM) * WQUAD(IQUAD)
!
!    Here Ar(IELEM) represents the area of element IELEM.
!
!    Input, character ( len = 2 ) EQN(MAXEQN).
!    EQN records the "type" of each equation that will be generated, 
!    and which is associated with an unknown.  Note that most 
!    boundary conditions do not result in an equation.  The current 
!    values are:
!
!    'U'  The horizontal momentum equation.
!    'UB' The condition U=0 applied at a node on the bump.
!    'UI' The condition U=UInflow(Y,Lambda) at the inflow.
!    'UW' The condition U=0 applied at a node on a fixed wall.
!
!    'V'  The vertical momentum equation.
!    'VB' The condition V=0 applied at a node on the bump.
!    'VI' The condition V=VInflow(Y,Lambda) at the inflow.
!    'VW' The condition V=0 applied at a node on a fixed wall.
!
!    'P'  The continuity equation.
!    'PB' The condition P=0 applied at (XMAX,YMAX).
!
!    Input, real ( kind = 8 ) G(MAXEQN, the current solution vector, in 
!    which are stored the finite element coefficients that define the 
!    velocity and pressure functions, U, V and P.
!
!    Input, integer INDX(3,MAXNP).  
!
!    INDX(I,J) contains, for each node J, the index of U, V and P 
!    at that node, or 0 or a negative value.
! 
!    If K=INDX(I,J) is positive, then the value of the degree
!    of freedom is stored in the solution vector entry G(K),
!    and an equation will be generated to determine its value.
!
!    If INDX(I,J) is not positive, then no equation is
!    generated to determine for variable I at node J, either 
!    because the variable is specified in some other way, or 
!    because (in the case of pressure), there is no coefficient 
!    associated with that node.
!
!    Input, integer NELEM, the number of elements.
! 
!    Input, integer NEQN, the number of finite element equations.
! 
!    Input, integer NLBAND, the lower bandwidth of the matrix A.  
! 
!    Input, integer NODE(6,MAXELM), the nodes that make up each element.
! 
!    Input, integer NP, the number of nodes.  NP=(2*NX-1)*(2*NY-1).
! 
!    Input, integer NROW, the number of rows need to store the matrix A.
!
!    Input, real ( kind = 8 ) PHI(MAXQUAD,6,10,MAXELM).  
!
!    PHI contains the value of a basis function, its derivative,
!    or other information, evaluated at a quadrature point.
! 
!    For a particular element I, quadrature point J, and basis
!    function K, we use the following shorthand for the 
!    entries of PHI:
!
!      W, dWdX, dWdY
!      Q, dQdX, dQdY
!      dXsidX, dXsidY, dEtadX, dEtadY
!
!    W is the quadratic basis function associated with velocity,
!    Q the linear basis function associated with pressure,
!    Xsi and Eta the reference coordinates for the point.
!        
!    In particular, PHI(J,K,1,I) is the value of the quadratic 
!    basis function associated with local node K in element I, 
!    evaluated at quadrature point J.
! 
!    Note that PHI(J,K,4,I)=PHI(J,K,5,I)=PHI(J,K,6,I)=0 for 
!    K=4, 5, or 6, since there are only three linear basis functions.
!
!    Input, real ( kind = 8 ) INV_NU, the value of the inverse dynamic viscosity.
!
  implicit none

  integer maxelm
  integer maxeqn
  integer maxnp
  integer maxquad
  integer nrow

  real ( kind = 8 ) a(nrow,maxeqn)
  real ( kind = 8 ) ar
  real ( kind = 8 ) area(maxquad,maxelm)
  real ( kind = 8 ) dpdx
  real ( kind = 8 ) dpdy
  real ( kind = 8 ) dqjdx
  real ( kind = 8 ) dqjdy
  real ( kind = 8 ) dudx
  real ( kind = 8 ) dudy
  real ( kind = 8 ) dvdx
  real ( kind = 8 ) dvdy
  real ( kind = 8 ) dwidx
  real ( kind = 8 ) dwidy
  real ( kind = 8 ) dwjdx
  real ( kind = 8 ) dwjdy
  character ( len = 2 ) eqn(maxeqn)
  real ( kind = 8 ) g(maxeqn)
  integer i
  integer ielem
  integer ihor
  integer indx(3,maxnp)
  integer ip
  integer iprs
  integer iq
  integer iquad
  integer iuse
  integer iver
  integer j
  integer jhor
  integer jp
  integer jprs
  integer jq
  integer jver
  integer nelem
  integer neqn
  integer nlband
  integer node(6,maxelm)
  integer np
  integer nquad
  real ( kind = 8 ) p
  real ( kind = 8 ) phi(maxquad,6,10,maxelm)
  real ( kind = 8 ) qi
  real ( kind = 8 ) inv_nu
  real ( kind = 8 ) term
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) wi
  real ( kind = 8 ) wj
!
!  Zero out the matrix.
!
  do i = 1, nrow
    do j = 1, neqn
      a(i,j) = 0.0D+00
    end do
  end do
!
!  Approximate the integral by summing over all elements.
!
  do ielem = 1, nelem
!
!  Evaluate the integrand at the quadrature points.
!
    do iquad = 1, nquad
 
      ar = area(iquad,ielem)
!
!  For the given quadrature point, evaluate P, U and V.
!
      call uvalq ( dpdx, dpdy, dudx, dudy, dvdx, dvdy, g, ielem, &
        indx, iquad, maxelm, maxeqn, maxnp, maxquad, &
        node, p, phi, u, v )
!
!  Look at all the basis functions in the element IELEM.
!
      do iq = 1, 6
 
        ip = node(iq,ielem)
 
        wi = phi(iquad,iq,1,ielem)
        dwidx = phi(iquad,iq,2,ielem)
        dwidy = phi(iquad,iq,3,ielem)
        qi = phi(iquad,iq,4,ielem)
 
        ihor = indx(1,ip)
        iver = indx(2,ip)
        iprs = indx(3,ip)
!
!  Now compute the derivatives of the functions associated
!  with U, V and P, with respect to the coefficients associated
!  with basis vectors at each node of the element.
!
        do jq = 1, 6
 
          jp = node(jq,ielem)
 
          wj = phi(iquad,jq,1,ielem)
          dwjdx = phi(iquad,jq,2,ielem)
          dwjdy = phi(iquad,jq,3,ielem)
 
          dqjdx = phi(iquad,jq,5,ielem)
          dqjdy = phi(iquad,jq,6,ielem)
 
          jhor = indx(1,jp)
          jver = indx(2,jp)
          jprs = indx(3,jp)
!
!  Contributions of the JHOR horizontal velocity to the U, V, and 
!  P equations.
!
          iuse = ihor - jhor + 2*nlband + 1

          if ( eqn(ihor) == 'U' ) then

            term = ar * ( dwjdx * dwidx + dwjdy * dwidy + &
              inv_nu * ( wj * dudx + u * dwjdx + v * dwjdy ) * wi )

            a(iuse,jhor) = a(iuse,jhor) + term

          end if
 
          if ( eqn(iver) == 'V' ) then
            iuse = iver - jhor + 2*nlband + 1
            term = ar * ( inv_nu * wj * dvdx * wi )
            a(iuse,jhor) = a(iuse,jhor) + term
          end if
 
          if ( 0 < iprs ) then
            if ( eqn(iprs) == 'P' ) then
              iuse = iprs - jhor + 2*nlband + 1
              term = ar * dwjdx * qi
              a(iuse,jhor) = a(iuse,jhor) + term
            end if
          end if
!
!  Contributions of the JVER vertical velocity variable to the
!  U, V and P equations.
!
          if ( eqn(ihor) == 'U' ) then
            iuse = ihor - jver + 2*nlband + 1
            term = ar * inv_nu * wj * dudy * wi
            a(iuse,jver) = a(iuse,jver) + term
          end if

          iuse = iver - jver + 2*nlband + 1
          if ( eqn(iver) == 'V' ) then
            term = ar * ( dwjdx * dwidx + dwjdy * dwidy + inv_nu * &
              ( u * dwjdx + wj * dvdy + v * dwjdy ) * wi )
            a(iuse,jver) = a(iuse,jver) + term
          end if
 
          if ( 0 < iprs ) then
            if ( eqn(iprs) == 'P' ) then
              iuse = iprs - jver + 2*nlband + 1
              term = ar * dwjdy * qi
              a(iuse,jver) = a(iuse,jver) + term
            end if
          end if
!
!  Contributions of the JPRS pressure to the U and V equations.
!
          if ( 0 < jprs ) then
 
            if ( eqn(ihor) == 'U' ) then
              iuse = ihor - jprs + 2*nlband + 1
              term = ar * inv_nu * dqjdx * wi
              a(iuse,jprs) = a(iuse,jprs) + term
            end if
 
            if ( eqn(iver) == 'V' ) then
              iuse = iver - jprs + 2*nlband + 1
              term = ar * inv_nu * dqjdy * wi
              a(iuse,jprs) = a(iuse,jprs) + term
            end if
 
          end if
 
        end do
      end do
    end do
  end do
!
!  Set up the equations that enforce boundary conditions.
!
  do ip = 1, np

    ihor = indx(1,ip)
    iver = indx(2,ip)
    iprs = indx(3,ip)

    if ( eqn(ihor) == 'UB' .or.eqn(ihor) == 'UI' .or. eqn(ihor) == 'UW' ) then

      a(2*nlband+1,ihor) = 1.0D+00

    end if

    if ( eqn(iver) == 'VB' .or. eqn(iver) == 'VI' .or. eqn(iver) == 'VW' ) then

      a(2*nlband+1,iver) = 1.0D+00

    end if

    if ( 0 < iprs ) then
      if ( eqn(iprs) == 'PB' ) then
        a(2*nlband+1,iprs) = 1.0D+00
      end if
    end if

  end do

  return
end
subroutine fx ( area, eqn, g, indx, maxelm, maxeqn, maxnp, maxquad, nelem, &
  neqn, node, nquad, phi, region, res, inv_nu, yc, yrange )

!*****************************************************************************80
!
!! FX computes the residual RES of the Navier Stokes equations.
!
!  Discussion:
!
!    The Navier Stokes variables are 
!
!      U, the horizontal velocity,
!      V, the vertical velocity,
!      P, the fluid pressure.
!
!    The finite element form of the Navier Stokes equations is:
!
!      Integral
!
!        dU/dx * dW/dx + dU/dy * dW/dy
!      + inv_nu * (U*dU/dx + V*dU/dy + dP/dx) * W dx dy = 0
!
!      Integral
!
!        dV/dx * dW/dx + dV/dy * dW/dy
!      + inv_nu * (U*dV/dx + V*dV/dy + dP/dy) * W dx dy = 0
!
!      Integral
!
!        (dU/dx + dV/dy) * Q dx dy = 0
!
!    where
!
!      INV_NU is the inverse of the dynamic viscosity,
!      W is a basis function for U and V, 
!      Q is a basis function for P.
!
!  Modified:
!
!    14 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) AREA(MAXQUAD,MAXELM).
!
!    AREA contains a common factor multiplying the term associated
!    with a quadrature point in a given element, namely,
!
!      AREA(IQUAD,IELEM) = Ar(IELEM) * WQUAD(IQUAD)
!
!    or, if the element is isoparametric,
!
!      AREA(IQUAD,IELEM) = DET * Ar(IELEM) * WQUAD(IQUAD)
!
!    Here Ar(IELEM) represents the area of element IELEM.
!
!    Input, character ( len = 2 ) EQN(MAXEQN).
!    EQN records the "type" of each equation that will be generated, 
!    and which is associated with an unknown.  Note that most 
!    boundary conditions do not result in an equation.  The current 
!    values are:
!
!    'U'  The horizontal momentum equation.
!    'UB' The condition U=0 applied at a node on the bump.
!    'UI' The condition U=UInflow(Y,Lambda) at the inflow.
!    'UW' The condition U=0 applied at a node on a fixed wall.
!
!    'V'  The vertical momentum equation.
!    'VB' The condition V=0 applied at a node on the bump.
!    'VI' The condition V=VInflow(Y,Lambda) at the inflow.
!    'VW' The condition V=0 applied at a node on a fixed wall.
!
!    'P'  The continuity equation.
!    'PB' The condition P=0 applied at (XMAX,YMAX).
!
!    Input, real ( kind = 8 ) G(MAXEQN).
!
!    G is the current solution vector, in which are stored 
!    the finite element coefficients that define the velocity
!    and pressure functions, U, V and P.
!
!    Input, integer INDX(3,MAXNP).  
!
!    INDX(I,J) contains, for each node J, the index of U, V and P 
!    at that node, or 0 or a negative value.
! 
!    If K=INDX(I,J) is positive, then the value of the degree
!    of freedom is stored in the solution vector entry G(K),
!    and an equation will be generated to determine its value.
!
!    If INDX(I,J) is not positive, then no equation is
!    generated to determine for variable I at node J, either 
!    because the variable is specified in some other way, or 
!    because (in the case of pressure), there is no coefficient 
!    associated with that node.
!
!    Input, integer NELEM, the number of elements.
! 
!    Input, integer NEQN, the number of finite element equations.
! 
!    Input, integer NODE(6,MAXELM), the nodes that make up each element.
! 
!    Input, integer NP, the number of nodes.  NP=(2*NX-1)*(2*NY-1).
! 
!    Input, real ( kind = 8 ) PHI(MAXQUAD,6,10,MAXELM).  
!
!    PHI contains the value of a basis function, its derivative,
!    or other information, evaluated at a quadrature point.
! 
!    For a particular element I, quadrature point J, and basis
!    function K, we use the following shorthand for the 
!    entries of PHI:
!
!      W, dWdX, dWdY
!      Q, dQdX, dQdY
!      dXsidX, dXsidY, dEtadX, dEtadY
!
!    W is the quadratic basis function associated with velocity,
!    Q the linear basis function associated with pressure,
!    Xsi and Eta the reference coordinates for the point.
!        
!    In particular, PHI(J,K,1,I) is the value of the quadratic 
!    basis function associated with local node K in element I, 
!    evaluated at quadrature point J.
! 
!    Note that PHI(J,K,4,I)=PHI(J,K,5,I)=PHI(J,K,6,I)=0 for 
!    K=4, 5, or 6, since there are only three linear basis
!    functions.
!
!    Input, character ( len = 20 ) REGION, the flow region, 'CHANNEL' 
!    or 'CAVITY'.
!
!    Output, real ( kind = 8 ) RES(MAXEQN), the residual.
!
!    Input, real ( kind = 8 ) INV_NU, the inverse dynamic viscosity.
!
!    Input, real ( kind = 8 ) YC(MAXNP), the Y coordinates of the nodes.
! 
!    Input, real ( kind = 8 ) YRANGE, the height of the region.
!
  implicit none

  integer maxelm
  integer maxeqn
  integer maxnp
  integer maxquad

  real ( kind = 8 ) ar
  real ( kind = 8 ) area(maxquad,maxelm)
  real ( kind = 8 ) dpdx
  real ( kind = 8 ) dpdy
  real ( kind = 8 ) dudx
  real ( kind = 8 ) dudy
  real ( kind = 8 ) dvdx
  real ( kind = 8 ) dvdy
  real ( kind = 8 ) dwidx
  real ( kind = 8 ) dwidy
  character ( len = 2 ) eqn(maxeqn)
  real ( kind = 8 ) g(maxeqn)
  integer i
  integer ielem
  integer ihor
  integer indx(3,maxnp)
  integer ip
  integer iprs
  integer iq
  integer iquad
  integer iver
  integer nelem
  integer neqn
  integer node(6,maxelm)
  integer nquad
  real ( kind = 8 ) p
  real ( kind = 8 ) phi(maxquad,6,10,maxelm)
  real ( kind = 8 ) qi
  character ( len = 20 ) region
  real ( kind = 8 ) res(maxeqn)
  real ( kind = 8 ) inv_nu
  real ( kind = 8 ) u
  real ( kind = 8 ) ubc
  real ( kind = 8 ) v
  real ( kind = 8 ) vbc
  real ( kind = 8 ) wi
  real ( kind = 8 ) yc(maxnp)
  real ( kind = 8 ) yrange
!
!  Zero out the residual vector.
!
  res(1:neqn) = 0.0D+00
!
!  Consider an element.
!
  do ielem = 1, nelem
!
!  Evaluate the integrand at each of the quadrature points.
!
    do iquad = 1, nquad
 
      ar = area(iquad,ielem)
!
!  Evaluate P, U and V.
!
      call uvalq ( dpdx, dpdy, dudx, dudy, dvdx, dvdy, g, ielem, &
        indx, iquad, maxelm, maxeqn, maxnp, maxquad, &
        node, p, phi, u, v )
!
!  Look at all basis functions in the element IELEM.
!
      do iq = 1, 6
 
        ip = node(iq,ielem)
 
        wi = phi(iquad,iq,1,ielem)
        dwidx = phi(iquad,iq,2,ielem)
        dwidy = phi(iquad,iq,3,ielem)
        qi = phi(iquad,iq,4,ielem)
!
!  The horizontal velocity equations.
!
        ihor = indx(1,ip)

        if ( eqn(ihor) == 'U' ) then

          res(ihor) = res(ihor) + ar * ( dudx * dwidx + dudy * dwidy &
            + inv_nu * ( u * dudx + v * dudy + dpdx ) * wi )

        else if ( eqn(ihor) == 'UB' ) then

          res(ihor) = g(ihor)

        else if ( eqn(ihor) == 'UI' ) then

          if ( region == 'CHANNEL' ) then
            ubc = yc(ip) * ( yrange - yc(ip) )
          else if ( region == 'CAVITY' ) then
            ubc = 1.0D+00
          end if

          res(ihor) = g(ihor) - ubc

        else if ( eqn(ihor) == 'UW' ) then

          res(ihor) = g(ihor)

        end if
!
!  The vertical velocity equations.
!
        iver = indx(2,ip)

        if ( eqn(iver) == 'V' ) then

          res(iver) = res(iver) + ar * ( dvdx * dwidx + dvdy * dwidy &
            + inv_nu * ( u * dvdx + v * dvdy + dpdy ) * wi )

        else if ( eqn(iver) == 'VB' ) then

          res(iver) = g(iver)

        else if ( eqn(iver) == 'VI' ) then

          vbc = 0.0D+00

          res(iver) = g(iver) - vbc

        else if ( eqn(iver) == 'VW' ) then

          res(iver) = g(iver)

        end if
!
!  The pressure equations.
!
        iprs = indx(3,ip)

        if ( 0 < iprs ) then
          if ( eqn(iprs) == 'P' ) then
            res(iprs) = res(iprs) + ar * ( dudx + dvdy ) * qi
          else if ( eqn(iprs) == 'PB' ) then
            res(iprs) = g(iprs)
          end if
        end if
 
      end do
    end do
  end do
 
  return
end
subroutine intprs ( g, indx, maxelm, maxeqn, maxnp, nelem, node, p )

!*****************************************************************************80
!
!! INTPRS interpolates the pressure at the midside nodes.  
!
!  Discussion:
!
!    This is only needed when writing out plot information.
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) G(MAXEQN).
!    G is the current solution estimate for the full problem, 
!    containing pressure and velocity coefficients.  The vector 
!    INDX must be used to index this data.
!
!    Input, integer INDX(3,MAXNP).  
!    INDX(I,J) contains, for each node J, the global index of U, 
!    V and P at that node, or 0 or a negative value.  The global
!    index of U, V, or P is the index of the coefficient vector
!    that contains the value of the finite element coefficient
!    associated with the corresponding basis function at the
!    given node.
!    If K=INDX(I,J) is positive, then the value of the degree
!    of freedom is stored in the solution vector entry GFL(K),
!    and an equation will be generated to determine its value.
!    If INDX(I,J) is not positive, then no equation is
!    generated to determine for variable I at node J, either 
!    because the variable is specified in some other way, or 
!    because (in the case of pressure), there is no coefficient 
!    associated with that node.
!
!    Input, integer NELEM, the number of elements.
! 
!    Input, integer NEQN, the number of finite element equations.
!
!    Input, integer NODE(6,MAXELM), the nodes that make up each element.
! 
!    Input, integer NP, the number of nodes.
!
!    Input, real P(MAXNP), the pressure.
!
  implicit none

  integer maxelm
  integer maxeqn
  integer maxnp

  real ( kind = 8 ) g(maxeqn)
  integer ielem
  integer in1
  integer in2
  integer in3
  integer in4
  integer in5
  integer in6
  integer indx(3,maxnp)
  integer nelem
  integer node(6,maxelm)
  real ( kind = 8 ) p(maxnp)
!
!  For each element,...
!
  do ielem = 1, nelem
!
!  Get the six global node numbers.
!
    in1 = node(1,ielem)
    in2 = node(2,ielem)
    in3 = node(3,ielem)
    in4 = node(4,ielem)
    in5 = node(5,ielem)
    in6 = node(6,ielem)
!
!  Read off the three computed values, and average the other three.
!
    p(in1) = g ( indx(3,in1) )
    p(in2) = g ( indx(3,in2) )
    p(in3) = g ( indx(3,in3) )
    p(in4) = 0.5D+00 * ( p(in1) + p(in2) )
    p(in5) = 0.5D+00 * ( p(in2) + p(in3) )
    p(in6) = 0.5D+00 * ( p(in3) + p(in1) )

  end do

  return
end
subroutine newton ( a, area, eqn, g, ierror, indx, ipivot, maxelm, maxeqn, &
  maxnew, maxnp, maxquad, nelem, neqn, nlband, node, np, nquad, nrow, phi, &
  region, res, inv_nu, tolnew, yc, yrange )

!*****************************************************************************80
!
!! NEWTON applies Newton iteration, seeking a solution of FX(G) = 0.
!
!  Discussion:
!
!    NEWTON is given an initial estimate of the solution of the nonlinear
!    state equations in G, and seeks a better solution.
!
!    The exact solution would have a zero residual, as computed by
!    the routine FX.  NEWTON uses Newton's method to seek a solution
!    whose maximum residual is no more than TOLNEW.  The routine FP
!    is used to compute the Jacobian of the residual functions.
!
!  Modified:
!
!    14 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Workspace, real ( kind = 8 ) A(NROW,MAXEQN), is used to store
!    the jacobian matrix, and then its LU factors.
!
!    Input, real ( kind = 8 ) AREA(MAXQUAD,NELEM).
!    AREA contains a common factor multiplying the term associated
!    with a quadrature point in a given element, namely,
!      AREA(IQUAD,IELEM) = Ar(IELEM) * WQUAD(IQUAD)
!    or, if the element is isoparametric,
!      AREA(IQUAD,IELEM) = DET * Ar(IELEM) * WQUAD(IQUAD)
!    Here Ar(IELEM) represents the area of element IELEM.
!
!    Input, character ( len = 2 ) EQN(MAXEQN).
!    EQN records the "type" of each equation that will be 
!    generated, and which is associated with an unknown.  Note that 
!    most boundary conditions do not result in an equation.  The 
!    current values are:
!
!    'U'  The horizontal momentum equation.
!    'UB' The condition U=0 applied at a node on the bump.
!    'UI' The condition U=UInflow(Y,Lambda) at the inflow.
!    'UW' The condition U=0 applied at a node on a fixed wall.
!
!    'V'  The vertical momentum equation.
!    'VB' The condition V=0 applied at a node on the bump.
!    'VI' The condition V=VInflow(Y,Lambda) at the inflow.
!    'VW' The condition V=0 applied at a node on a fixed wall.
!
!    'P'  The continuity equation.
!    'PB' The condition P=0 applied at (XMAX,YMAX).
!
!    Input/output, real ( kind = 8 ) G(MAXEQN).
!
!    On input, some estimate for the solution.  If no estimate
!    is known, set G to zero.
!
!    On output, G is the improved solution computed by Newton's method.
!
!    Output, integer IERROR, error flag.
!    0, no error occurred.
!    1, an error occurred, and the improved solution could not be computed.
!
!    Input, integer INDX(3,MAXNP).  
!
!    INDX(I,J) contains, for each node J, the index of U, V and P 
!    at that node, or 0 or a negative value.
! 
!    If K=INDX(I,J) is positive, then the value of the degree
!    of freedom is stored in the solution vector entry G(K),
!    and an equation will be generated to determine its value.
!
!    If INDX(I,J) is not positive, then no equation is
!    generated to determine for variable I at node J, either 
!    because the variable is specified in some other way, or 
!    because (in the case of pressure), there is no coefficient 
!    associated with that node.
!
!    Workspace, integer IPIVOT(MAXEQN), pivot space needed by the
!    matrix factorization and solving routines.
!
!    Input, integer MAXNEW, the maximum number of Newton iterations.
!
!    Input, integer NELEM, the number of elements.
! 
!    Input, integer NEQN, the number of finite element equations.
! 
!    Input, integer NLBAND, the lower bandwidth of the matrix A.
! 
!    Input, integer NODE(6,MAXELM), the nodes that make up each element.
! 
!    Input, integer NP, the number of nodes.  NP=(2*NX-1)*(2*NY-1).
! 
!    Input, integer NROW, the number of rows need to store the matrix A.
!
!    Input, real ( kind = 8 ) PHI(MAXQUAD,6,10,MAXELM).  
!
!    PHI contains the value of a basis function, its derivative,
!    or other information, evaluated at a quadrature point.
! 
!    For a particular element I, quadrature point J, and basis
!    function K, we use the following shorthand for the 
!    entries of PHI:
!
!    For a particular element I, quadrature point J, and basis
!    function K, we use the following shorthand for the 
!    entries of PHI:
!
!      W, dWdX, dWdY
!      Q, dQdX, dQdY
!      dXsidX, dXsidY, dEtadX, dEtadY
!
!    W is the quadratic basis function associated with velocity,
!    Q the linear basis function associated with pressure,
!    Xsi and Eta the reference coordinates for the point.
!        
!    In particular, PHI(J,K,1,I) is the value of the quadratic 
!    basis function associated with local node K in element I, 
!    evaluated at quadrature point J.
! 
!    Note that PHI(J,K,4,I)=PHI(J,K,5,I)=PHI(J,K,6,I)=0 for 
!    K=4, 5, or 6, since there are only three linear basis
!    functions.
!
!    Input, character ( len = 20 ) REGION, the flow region, 'CHANNEL' 
!    or 'CAVITY'.
!
!    Workspace, real ( kind = 8 ) RES(MAXEQN), the residual.
!
!    Input, real ( kind = 8 ) INV_NU, the inverse dynamic viscosity.
!
!    Input, real ( kind = 8 ) TOLNEW, the Newton tolerance.
!    NEWTON is asked to find an approximate solution so that
!    the maximum absolute value of all the residuals is no more
!    than TOLNEW.  A value such as 10E-7 is often reasonable,
!    though this depends on the actual equations being solved.
!
!    Input, real ( kind = 8 ) YC(MAXNP), the Y coordinates of the nodes.
!
!    Input, real ( kind = 8 ) YRANGE, the height of the region.
! 
  implicit none

  integer maxelm
  integer maxeqn
  integer maxnp
  integer maxquad
  integer nrow

  real ( kind = 8 ) a(nrow,maxeqn)
  real ( kind = 8 ) area(maxquad,maxelm)
  real ( kind = 8 ) dmax
  character ( len = 2 ) eqn(maxeqn)
  real ( kind = 8 ) g(maxeqn)
  integer i
  integer idmax
  integer ierror
  integer indx(3,maxnp)
  integer info
  integer ipivot(maxeqn)
  integer irmax
  integer iter
  integer ixmax
  integer job
  integer maxnew
  integer nelem
  integer neqn
  integer nlband
  integer node(6,maxelm)
  integer np
  integer nquad
  real ( kind = 8 ) phi(maxquad,6,10,maxelm)
  character ( len = 20 ) region
  real ( kind = 8 ) res(maxeqn)
  real ( kind = 8 ) inv_nu
  real ( kind = 8 ) rmax
  real ( kind = 8 ) rmax0
  real ( kind = 8 ) tolnew
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmax0
  real ( kind = 8 ) yc(maxnp)
  real ( kind = 8 ) yrange

  ierror = 0
  iter = 0
!
!  Compute the max-norm of the initial X value.
!
  call r8vec_amax ( neqn, g, xmax, ixmax )
  xmax0 = xmax
!
!  Evaluate the residual RES of the initial X value.
!
  call fx ( area, eqn, g, indx, maxelm, maxeqn, maxnp, maxquad, nelem, neqn, &
    node, nquad, phi, region, res, inv_nu, yc, yrange )
 
  call r8vec_amax ( neqn, res, rmax, irmax )
  rmax0 = rmax
!
!  Carry out MAXNEW steps of Newton iteration.
!
  do iter = 1, maxnew
!
!  Get the jacobian matrix A.
!
    call fp ( a, area, eqn, g, indx, maxelm, maxeqn, maxnp, maxquad, &
      nelem, neqn, nlband, node, np, nquad, nrow, phi, inv_nu )
!
!  Factor the jacobian matrix.
!
    call dgb_fa ( a, nrow, neqn, nlband, nlband, ipivot, info )

    if ( info /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'Newton - Fatal error!'
      write ( *, '(a)' ) '  The jacobian is singular.'
      write ( *, '(a,i8)' ) '  The factor routine returns INFO=', info
      ierror = 1
      return
    end if
!
!  Solve the linear system A*DX=RES
!
    job = 0
    call dgb_sl ( a, nrow, neqn, nlband, nlband, ipivot, res, job )
 
    call r8vec_amax ( neqn, res, dmax, idmax )
!
!  The next iterate is GNEW = G - DX.
!
    g(1:neqn) = g(1:neqn) - res(1:neqn)

    call r8vec_amax ( neqn, g, xmax, ixmax )
!
!  Evaluate the residual FX of the current estimated solution.
!
    call fx ( area, eqn, g, indx, maxelm, maxeqn, maxnp, maxquad, &
      nelem, neqn, node, nquad, phi, region, res, inv_nu, yc, yrange )
 
    call r8vec_amax ( neqn, res, rmax, irmax )
!
!  Accept the iterate if the residual is small enough.
!
    if ( rmax <= tolnew ) then
      return
    end if
!
!  Reject the iterate if the residual has grown too large.
!
    if ( 10.0D+00 * ( rmax0 + tolnew ) < rmax .and. 1 < iter ) then
      ierror = 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'Newton - Warning!'
      write ( *, '(a,i8)' ) '  Residual too big on step     ', iter
      write ( *, '(a,g14.6)' ) '  The final stepsize was        ', dmax
      write ( *, '(a,g14.6)' ) '  The initial X norm was        ', xmax0
      write ( *, '(a,g14.6)' ) '  The final X norm was          ', xmax
      write ( *, '(a,g14.6)' ) '  Initial residual =            ', rmax0
      write ( *, '(a,g14.6)' ) '  Current residual =            ', rmax
      return
    end if
 
  end do
!
!  The iteration has failed to converge, or may actually
!  have been terminated early.
! 
  ierror = 1
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Newton - Warning!'
  write ( *, '(a,i8,a)' ) '  No convergence after ', maxnew, ' steps.'
  write ( *, '(a,g14.6)' ) '  The final stepsize was        ', dmax
  write ( *, '(a,g14.6)' ) '  The initial X norm was        ', xmax0
  write ( *, '(a,g14.6)' ) '  The final X norm was          ', xmax
  write ( *, '(a,g14.6)' ) '  The initial residual norm was ', rmax0
  write ( *, '(a,g14.6)' ) '  The final residual norm was   ', rmax
 
  return
end
subroutine qbf ( ielem, in, w, dwdx, dwdy, maxelm, maxnp, node, xc, xq, yc, yq )

!*****************************************************************************80
!
!! QBF evaluates a nonisoparametric quadratic basis function.
!
!  Diagram:
!
!      ^
!      |    2
!      |    |\
!   Y  |    5 4
!      |    |  \
!      |    3-6-1
!      |
!      +------------>
!             X
!
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer IELEM, the index of the element.
!
!    Input, integer IN, the index of the basis function.
!    This will be a value between 1 and 6.  Functions
!    1 through 3 are associated with corners, 4 though 6 with sides.
!
!    Output, real ( kind = 8 ) W, DWDX, DWDY, the value of the
!    IN-th basis  function and its X and Y derivatives, at the given point.
!
!    Input, integer NELEM, the number of elements.
!
!    Input, integer NODE(6,MAXELM), the nodes that make up each element.
!
!    Input, integer NP, the number of nodes.
!
!    Input, real ( kind = 8 ) XC(MAXNP), the X coordinates of the nodes.
!
!    Input, real ( kind = 8 ) XQ, the X coordinate of the evaluation point.
!
!    Input, real ( kind = 8 ) YC(MAXNP), the Y coordinates of the nodes.
!
!    Input, real ( kind = 8 ) YQ, the Y coordinate of the evaluation point.
!
  implicit none

  integer maxelm
  integer maxnp

  real ( kind = 8 ) c
  real ( kind = 8 ) d
  real ( kind = 8 ) dwdx
  real ( kind = 8 ) dwdy
  integer i1
  integer i2
  integer i3
  integer ielem
  integer in
  integer in1
  integer in2
  integer in3
  integer node(6,maxelm)
  real ( kind = 8 ) s
  real ( kind = 8 ) t
  real ( kind = 8 ) w
  real ( kind = 8 ) xc(maxnp)
  real ( kind = 8 ) xq
  real ( kind = 8 ) yc(maxnp)
  real ( kind = 8 ) yq
!
!  Case 1: We are inquiring about a basis function associated
!  with a corner.
!
!  Notice that the basis function W is zero exactly if
!  T is 0 or T is 1/2.
!
!  IN1, IN2, and IN3 are the local node numbers of the three
!  corner nodes, and I1, I2 and I3 are the corresponding
!  global node numbers, which are used to look up the X and
!  Y coordinates of the nodes.
!
  if ( 1 <= in .and. in <= 3 ) then
 
    in1 = in
    in2 = mod ( in, 3 ) + 1
    in3 = mod ( in+1, 3 ) + 1
 
    i1 = node(in1,ielem)
    i2 = node(in2,ielem)
    i3 = node(in3,ielem)
 
    d = ( xc(i2) - xc(i1) ) * ( yc(i3) - yc(i1) ) &
      - ( xc(i3) - xc(i1) ) * ( yc(i2) - yc(i1) )
 
    t = 1.0D+00 + ( &
        ( xq     - xc(i1) ) * ( yc(i2) - yc(i3) ) &
      + ( xc(i3) - xc(i2) ) * ( yq     - yc(i1) ) ) / d
 
    w = t * ( 2.0D+00 * t - 1.0D+00 )
 
    dwdx = ( yc(i2) - yc(i3) ) * ( 4.0D+00 * t - 1.0D+00 ) / d
    dwdy = ( xc(i3) - xc(i2) ) * ( 4.0D+00 * t - 1.0D+00 ) / d
!
!  Case 2: We are inquiring about a basis function associated
!  with a midpoint.
!
  else if ( 4 <= in .and. in <= 6 ) then
 
    in1 = in - 3
    in2 = mod ( in-3, 3 ) + 1
    in3 = mod ( in-2, 3 ) + 1
 
    i1 = node(in1,ielem)
    i2 = node(in2,ielem)
    i3 = node(in3,ielem)
 
    d =     ( xc(i2) - xc(i1) ) * ( yc(i3) - yc(i1) ) &
          - ( xc(i3) - xc(i1) ) * ( yc(i2) - yc(i1) )
 
    c =     ( xc(i3) - xc(i2) ) * ( yc(i1) - yc(i2) ) &
          - ( xc(i1) - xc(i2) ) * ( yc(i3) - yc(i2) )
 
    t = 1.0 + ( &
        ( xq     - xc(i1) ) * ( yc(i2) - yc(i3) ) &
      + ( xc(i3) - xc(i2) ) * ( yq     - yc(i1) ) ) / d
 
    s = 1.0 + ( &
        ( xq     - xc(i2) ) * ( yc(i3) - yc(i1) ) &
      + ( xc(i1) - xc(i3) ) * ( yq     - yc(i2) ) ) / c
 
    w = 4.0D+00 * s * t
    dwdx = 4.0D+00 * ( ( yc(i3) - yc(i1) ) * t / c &
                     + ( yc(i2) - yc(i3) ) * s / d )
    dwdy = 4.0D+00 * ( ( xc(i1) - xc(i3) ) * t / c &
                     + ( xc(i3) - xc(i2) ) * s / d )
 
  else
 
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'QBF - Fatal error!'
    write ( *, '(a,i8)' ) '  Request for basis function IN = ', in
    write ( *, '(a)' ) '  but IN must be between 1 and 6.'
    stop
 
  end if
 
  return
end
subroutine r8vec_amax ( n, x, xmax, ixmax )

!*****************************************************************************80
!
!! R8VEC_AMAX returns the maximum absolute value in a real vector.
!
!  Modified:
!
!    22 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the array.
!
!    Input, real ( kind = 8 ) X(N), the array.
!
!    Output, real ( kind = 8 ) XMAX, the largest absolute value of the entries in
!    the array.
!
!    Output, integer IXMAX, the index of the entry of largest absolute value.
!
  implicit none

  integer n

  integer i
  integer ixmax
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xmax

  if ( n <= 0 ) then

    xmax = 0.0D+00
    ixmax = 0

  else

    xmax = abs ( x(1) )
    ixmax = 1

    do i = 2, n

      if ( xmax < abs ( x(i) ) ) then
        xmax = abs ( x(i) )
        ixmax = i
      end if

    end do

  end if

  return
end
subroutine refbsp ( q, dqdx, dqdy, detadx, detady, iq, dxsidx, dxsidy, eta, &
  xsi )

!*****************************************************************************80
!
!! REFBSP evaluates a reference element linear basis function.
!
!  Discussion:
!
!    The basis function and its X and Y derivatives are evaluated at a 
!    particular point (X,Y) in a particular element, by referring to the 
!    corresponding points (XSI,ETA) in the reference triangle.
!
!    It is assumed that we already know the value of the jacobian
!    of the isoparametric transformation between the (XSI, ETA) and
!    (X, Y) spaces.  The four entries of the jacobian are
!    symbolically named DETADX, DETADY, DXSIDX and DXSIDY, and
!    we know that the jacobian gives us the following relation
!    between derivatives with respect to XSI and ETA, and derivatives
!    with respect to X and Y:
!
!      dF/dX = dF/dXsi dXsi/dX + dF/dEta dEta/dX
!      dF/dY = dF/dXsi dXsi/dY + dF/dEta dEta/dY
!
!    Here is a graph of the (XSI, ETA) reference triangle.
!
!          ^
!          |
!        1 +        2
!          |       /|
!    ETA   |      / |
!          |     /  |
!        0 +    1---3
!          |
!          +----+---+--->
!               0   1
!
!                XSI
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) Q, DQDX, DQDY, the value of the basis
!    function, and its derivatives with respect to X and Y, at
!    the point (ETA,XSI).
!
!    Input, real ( kind = 8 ) DETADX, DETADY, the partial derivative
!    d ETA/d X and d ETA/d Y at (ETA,XSI).
!
!    Input, integer IQ, the local node number, between 1 and
!    3, whose basis function is being evaluated.
!
!    Input, real ( kind = 8 ) DXSIDX, DXSIDY, the partial derivative
!    d XSI/d X and d XSI/d Y at (ETA,XSI).
!
!    Input, real ( kind = 8 ) ETA, XSI, the local coordinates of the
!    point at which the basis information is desired.
!
  implicit none

  real ( kind = 8 ) detadx
  real ( kind = 8 ) detady
  real ( kind = 8 ) dqdeta
  real ( kind = 8 ) dqdx
  real ( kind = 8 ) dqdxsi
  real ( kind = 8 ) dqdy
  real ( kind = 8 ) dxsidx
  real ( kind = 8 ) dxsidy
  real ( kind = 8 ) eta
  integer iq
  real ( kind = 8 ) q
  real ( kind = 8 ) xsi

  if ( iq == 1 ) then
    q = 1.0D+00 - xsi
    dqdxsi = -1.0D+00
    dqdeta =  0.0D+00
  else if ( iq == 2 ) then
    q = eta
    dqdxsi = 0.0D+00
    dqdeta = 1.0D+00
  else if ( iq == 3 ) then
    q = xsi - eta
    dqdxsi = 1.0D+00
    dqdeta = -1.0D+00
  else if ( 4 <= iq .and. iq <= 6 ) then
    q = 0.0D+00
    dqdxsi = 0.0D+00
    dqdeta = 0.0D+00
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RefBSP - Fatal error!'
    write ( *, '(a,i8)' ) '  Request for basis function IQ=', iq
    write ( *, '(a)' ) '  but IQ must be between 1 and 6.'
    stop
  end if
 
  dqdx = dqdxsi * dxsidx + dqdeta * detadx
  dqdy = dqdxsi * dxsidy + dqdeta * detady
 
  return
end
subroutine refqbf ( w, dwdx, dwdy, detadx, detady, dxsidx, dxsidy, eta, iq, &
  xsi )

!*****************************************************************************80
!
!! REFQBF evaluates a reference element quadratic basis function.
!
!
!  There are six possible quadratic basis functions.  This routine
!  evaluates just one of them, and its X and Y derivatives, at a 
!  particular point in a particular element, by referring to the 
!  reference triangle.
!
!  The point we are interested in is referred to by its coordinates
!  in the reference triangle.  That is, we are given coordinates
!  (XSI, ETA), even though, physically, we are interested
!  in points in (X, Y) space.
!
!  It is assumed that we already know the value of the jacobian
!  of the isoparametric transformation between the (XSI, ETA) and
!  (X, Y) spaces.  The four entries of the jacobian are
!  symbolically named DETADX, DETADY, DXSIDX and DXSIDY, and
!  we know that the jacobian gives us the following relation
!  between derivatives with respect to XSI and ETA, and derivatives
!  with respect to X and Y:
!
!    d F(X,Y)/dX     (d XSI/dX  d ETA/dX )   ( d F(XSI, ETA)/d XSI )
!    d F(X,Y)/dY  =  (d XSI/dY  d ETA/dY ) * ( d F(XSI, ETA)/d ETA )
!
!  Here is a graph of the (XSI, ETA) reference triangle.
!
!        ^
!        |
!      1 +        2
!        |       /|
!  ETA   |      4 5
!        |     /  |
!      0 +    1-6-3
!        |
!        +----+---+--->
!             0   1
!
!              XSI
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) W, DWDX, DWDY, the value of the basis
!    function, and its derivatives with respect to X and Y, at
!    the point (XSI,ETA).
!
!    Input, real ( kind = 8 ) DETADX, DETADY, the partial derivative
!    d ETA/d X and d ETA/d Y at (XSI,ETA).
!
!    Input, real ( kind = 8 ) DXSIDX, DXSIDY, the partial derivative
!    d XSI/d X and d XSI/d Y at (XSI,ETA).
!
!    Input, real ( kind = 8 ) ETA, the ETA coordinate of the point.
!
!    Input, integer IQ, the local node number, between 1 and
!    6, whose basis function is being evaluated.
!
!    Input, real ( kind = 8 ) XSI, the XSI coordinate of the point.
!
  implicit none
!
  real ( kind = 8 ) detadx
  real ( kind = 8 ) detady
  real ( kind = 8 ) dwdeta
  real ( kind = 8 ) dwdx
  real ( kind = 8 ) dwdxsi
  real ( kind = 8 ) dwdy
  real ( kind = 8 ) dxsidx
  real ( kind = 8 ) dxsidy
  real ( kind = 8 ) eta
  integer iq
  real ( kind = 8 ) w
  real ( kind = 8 ) xsi
!
!  Evaluate W, the quadratic basis function.
!  Evaluate DWDXSI and DWDETA, the partial derivatives d W/d XSI
!  and d W/d ETA.
!
!  Basis 1 is zero if XSI=0.5 or XSI=1.
!
  if ( iq == 1 ) then

    w = ( 2.0D+00 * xsi - 1.0D+00 ) * ( xsi - 1.0D+00 )
    dwdxsi = - 3.0D+00 + 4.0D+00 * xsi
    dwdeta = 0.0D+00
!
!  Basis 2 is zero if ETA=0 or ETA=0.5.
!
  else if ( iq == 2 ) then

    w = eta * ( 2.0D+00 * eta - 1.0D+00 )
    dwdxsi = 0.0D+00
    dwdeta = - 1.0D+00 + 4.0D+00 * eta
!
!  Basis 3 is zero if XSI=ETA, or XSI=ETA+0.5D+00
!
  else if ( iq == 3 ) then

    w = ( xsi - eta ) * ( 2.0D+00 * xsi - 2.0D+00 * eta - 1.0D+00 )
    dwdxsi = - 1.0D+00 + 4.0D+00 * xsi - 4.0D+00 * eta
    dwdeta = 1.0D+00 - 4.0D+00 * xsi + 4.0D+00 * eta
!
!  Basis 4 is zero if ETA=0 or XSI=1.
!
  else if ( iq == 4 ) then

    w = 4.0D+00 * eta * ( 1.0D+00 - xsi )
    dwdxsi = - 4.0D+00 * eta
    dwdeta = 4.0D+00 - 4.0D+00 * xsi
!
!  Basis 5 is zero if ETA=0 or XSI=ETA.
!
  else if ( iq == 5 ) then

    w = 4.0D+00 * eta * ( xsi - eta )
    dwdxsi = 4.0D+00 * eta
    dwdeta = 4.0D+00 * xsi - 8.0D+00 * eta
!
!  Basis 6 is zero if XSI=ETA or XSI=1.
!
  else if ( iq == 6 ) then

    w = 4.0D+00 * ( xsi - eta ) * ( 1.0D+00 - xsi )
    dwdxsi = 4.0D+00 - 8.0D+00 * xsi + 4.0D+00 * eta
    dwdeta = - 4.0D+00 + 4.0D+00 * xsi
!
!  Stop if we were given an unexpected value of IQ.
!
  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RefQBF - Fatal error!'
    write ( *, '(a)' ) '  A basis function index must be between 1 and 6,'
    write ( *, '(a,i8)' ) '  but you input the value IQ=', iq
    stop

  end if
!
!  Convert the d W/d XSI and d W/d ETA derivatives to d W/d X
!  and d W/d Y.
!
  dwdx = dwdxsi * dxsidx + dwdeta * detadx
  dwdy = dwdxsi * dxsidy + dwdeta * detady
 
  return
end
subroutine setban ( indx, maxelm, maxnp, maxrow, nelem, nlband, node, nrow )

!*****************************************************************************80
!
!! SETBAN computes the lower band width of the Jacobian matrix.
!
!  Discussion:
!
!    It also finds NROW, the total number of rows required to store the 
!    matrix in LINPACK general band storage format.
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer INDX(3,MAXNP).  
!
!    INDX(I,J) contains, for each node J, the index of U, V and P 
!    at that node, or 0 or a negative value.
! 
!    If K=INDX(I,J) is positive, then the value of the degree
!    of freedom is stored in the solution vector entry G(K),
!    and an equation will be generated to determine its value.
!
!    If INDX(I,J) is not positive, then no equation is
!    generated to determine for variable I at node J, either 
!    because the variable is specified in some other way, or 
!    because (in the case of pressure), there is no coefficient 
!    associated with that node.
!
!    Input, integer MAXROW, the first dimension of the matrix A.
!
!    Input, integer NELEM, the number of elements.
! 
!    Output, integer NLBAND, the lower bandwidth of the matrix A.
! 
!    Input, integer NODE(6,MAXELM), the nodes that make up each element.
! 
!    Input, integer NP, the number of nodes.  NP=(2*NX-1)*(2*NY-1).
! 
!    Output, integer NROW, the number of rows needed to store the matrix A.
!
  implicit none

  integer maxelm
  integer maxnp

  integer i
  integer ielem
  integer indx(3,maxnp)
  integer ip
  integer ipp
  integer iq
  integer iqq
  integer iuk
  integer iukk
  integer j
  integer maxrow
  integer nelem
  integer nlband
  integer node(6,maxelm)
  integer nrow

  nlband = 0
 
  do ielem = 1, nelem
    do iq = 1, 6
      ip = node(iq,ielem)
      do iuk = 1, 3
        i = indx(iuk,ip)
        if ( 0 < i ) then
          do iqq = 1, 6
            ipp = node(iqq,ielem)
            do iukk = 1, 3
              j = indx(iukk,ipp)
              if ( 0 < j ) then
                if ( nlband < j-i ) then
                  nlband = j - i
                end if
              end if
            end do
          end do
        end if
      end do
    end do
  end do
 
  nrow = 3 * nlband + 1
 
  if ( maxrow < nrow ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SetBan - Fatal error!'
    write ( *, '(a,i8)' ) '  NROW is too large!  NROW=     ', nrow
    write ( *, '(a,i8)' ) '  The maximum allowed is MAXROW=', maxrow
    stop
  end if
 
  return
end
subroutine setbas ( area, etaquad, isotri, maxelm, maxnp, maxquad, &
  nelem, node, nquad, phi, xc, xquad, xsiquad, yc, yquad )

!*****************************************************************************80
!
!! SETBAS evaluates the basis functions at each quadrature point.  
!
!  Discussion:
!
!    The basis functions are computed and saved in this way for efficiency.
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) AREA(MAXQUAD,MAXELM).
!
!    AREA contains the area of each element.  These values are
!    needed when computed the integrals associated with the
!    finite element method.
!
!    For runs in which the region is allowed to change from
!    step to step, AREA must be recalculated at each step.
!
!    Input, real ( kind = 8 ) ETAQUAD(MAXQUAD), the "Eta" quadrature coordinates.
!
!    Input, integer ISOTRI(MAXELM), 0/1, the element IS NOT/IS isoparametric.
!
!    Input, integer NELEM, the number of elements.
!
!    Input, integer NODE(6,MAXELM), the nodes that make up each element.
!
!    Input, integer NP, the number of nodes.
!
!    Output, real ( kind = 8 ) PHI(MAXQUAD,6,10,MAXELM),
!    contains lots of basis function values.  In particular,
!
!    PHI(I,J,K,1) contains the value of the basis function
!    associated with velocity (U or V), in the I-th element,
!    at the J-th quadrature point, associated with the
!    K-th node.
!
!    PHI(I,J,K,2) contains the X derivative, and
!    PHI(I,J,K,3) contains the Y derivative.
!
!    PHI(I,J,K,4) contains the value of the basis function
!    associated with pressure (P) in the I-th element,
!    at the J-th quadrature point, associated with the K-th node.
!
!    Input, real ( kind = 8 ) XC(MAXNP), the X coordinates of the nodes.
!
!    Input, real ( kind = 8 ) XQUAD(MAXQUAD,NELEM), the X quadrature coordinates.
!
!    Input, real ( kind = 8 ) XSIQUAD(MAXQUAD), the "Xsi" quadrature coordinates.
!
!    Input, real ( kind = 8 ) YC(MAXNP), the Y coordinates of the nodes.
!
!    Input, real ( kind = 8 ) YQUAD(MAXQUAD,MAXELM), the Y quadrature coordinates.
!
  implicit none

  integer maxelm
  integer maxnp
  integer maxquad

  real ( kind = 8 ) area(maxquad,maxelm)
  real ( kind = 8 ) det
  real ( kind = 8 ) detadx
  real ( kind = 8 ) detady
  real ( kind = 8 ) dqdx
  real ( kind = 8 ) dqdy
  real ( kind = 8 ) dwdx
  real ( kind = 8 ) dwdy
  real ( kind = 8 ) dxsidx
  real ( kind = 8 ) dxsidy
  real ( kind = 8 ) eta
  real ( kind = 8 ) etaquad(maxquad)
  integer ielem
  integer iq
  integer iquad
  integer isotri(maxelm)
  integer nelem
  integer node(6,maxelm)
  integer nquad
  real ( kind = 8 ) phi(maxquad,6,10,maxelm)
  real ( kind = 8 ) q
  real ( kind = 8 ) w
  real ( kind = 8 ) xc(maxnp)
  real ( kind = 8 ) xquad(maxquad,maxelm)
  real ( kind = 8 ) xq
  real ( kind = 8 ) xsi
  real ( kind = 8 ) xsiquad(maxquad)
  real ( kind = 8 ) yc(maxnp)
  real ( kind = 8 ) yquad(maxquad,maxelm)
  real ( kind = 8 ) yq
!
!  Consider a particular element,
!  and a particular quadrature point (XQ,YQ) in that element.
!
!  Compute, at (XQ,YQ), the local values of the jacobian matrix
!  and its determinant.
!
!  Adjust the AREA array   
!
  do ielem = 1, nelem
 
    do iquad = 1, nquad
 
      xq = xquad(iquad,ielem)
      yq = yquad(iquad,ielem)
 
      if ( isotri(ielem) == 1 ) then

        eta = etaquad(iquad)
        xsi = xsiquad(iquad)

        call trans ( det, detadx, detady, dxsidx, dxsidy, eta, &
          ielem, maxelm, maxnp, node, xc, xsi, yc )

        area(iquad,ielem) = det * area(iquad,ielem)

      end if
!
!  Now consider each of the basis functions associated with a
!  node in the given element.
!
      do iq = 1, 6
!
!  If the element is NOT isoparametric, compute the basis values
!  directly.
!
        if ( isotri(ielem) == 0 ) then
 
          call bsp ( q, dqdx, dqdy, ielem, iq, maxelm, maxnp, &
            node, xc, xq, yc, yq )
 
          call qbf ( ielem, iq, w, dwdx, dwdy, maxelm, maxnp, &
            node, xc, xq, yc, yq )
 
          dxsidx = 1.0D+00
          dxsidy = 0.0D+00
          detadx = 0.0D+00
          detady = 1.0D+00
!
!  For isoparametric elements, use the reference triangle method.
!
        else
 
          call refqbf ( w, dwdx, dwdy, detadx, detady, dxsidx, &
            dxsidy, eta, iq, xsi )
 
          call refbsp ( q, dqdx, dqdy, detadx, detady, iq, dxsidx, &
            dxsidy, eta, xsi )
 
        end if
!
!  Store the values into PHI.
!
        phi(iquad,iq,1,ielem) = w
        phi(iquad,iq,2,ielem) = dwdx
        phi(iquad,iq,3,ielem) = dwdy
        phi(iquad,iq,4,ielem) = q
        phi(iquad,iq,5,ielem) = dqdx
        phi(iquad,iq,6,ielem) = dqdy
 
        phi(iquad,iq,7,ielem) = dxsidx
        phi(iquad,iq,8,ielem) = dxsidy
        phi(iquad,iq,9,ielem) = detadx
        phi(iquad,iq,10,ielem) = detady     

      end do
    end do
  end do

  return
end
subroutine setquad ( area, etaquad, isotri, maxelm, maxnp, maxquad, nelem, &
  node, nquad, wquad, xc, xquad, xsiquad, yc, yquad )

!*****************************************************************************80
!
!! SETQUAD sets the abscissas and weights for quadrature in a triangle.
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) AREA(MAXQUAD,MAXELM).
!
!    AREA contains a common factor multiplying the term associated
!    with a quadrature point in a given element, namely,
!
!      AREA(IQUAD,IELEM) = Ar(IELEM) * WQUAD(IQUAD)
!
!    or, if the element is isoparametric,
!
!      AREA(IQUAD,IELEM) = DET * Ar(IELEM) * WQUAD(IQUAD)
!
!    Here Ar(IELEM) represents the area of element IELEM.
!
!    Output, real ( kind = 8 ) ETAQUAD(MAXQUAD), the "Eta" quadrature
!    coordinates.
!
!    Input, integer ISOTRI(MAXELM), 0/1, the element IS NOT/IS isoparametric.
! 
!    Input, integer NELEM, the number of elements.
!
!    Input, integer NODE(6,MAXELM), the nodes that make up each element.
!
!    Input, integer NP, the number of nodes.
!
!    Output, real ( kind = 8 ) WQUAD(MAXQUAD), the quadrature weights.
!
!    Input, real ( kind = 8 ) XC(MAXNP), the X coordinates of the nodes.
!
!    Output, real ( kind = 8 ) XQUAD(MAXQUAD,MAXELM), the X quadrature coordinates.
!
!    Output, real ( kind = 8 ) XSIQUAD(MAXQUAD), the "Xsi" quadrature coordinates.
!
!    Input, real ( kind = 8 ) YC(MAXNP), the Y coordinates of the nodes.
!
!    Output, real ( kind = 8 ) YQUAD(MAXQUAD,MAXELM), the Y quadrature coordinates. 
!
  implicit none

  integer maxelm
  integer maxnp
  integer maxquad

  real ( kind = 8 ) area(maxquad,maxelm)
  real ( kind = 8 ) eta
  real ( kind = 8 ) etaquad(maxquad)
  integer i
  integer ielem
  integer ip1
  integer ip2
  integer ip3
  integer iquad
  integer isotri(maxelm)
  integer nelem
  integer node(6,maxelm)
  integer nquad
  real ( kind = 8 ) wquad(maxquad)
  real ( kind = 8 ) x
  real ( kind = 8 ) xc(maxnp)
  real ( kind = 8 ) xquad(maxquad,maxelm)
  real ( kind = 8 ) xsi
  real ( kind = 8 ) xsiquad(maxquad)
  real ( kind = 8 ) y
  real ( kind = 8 ) yc(maxnp)
  real ( kind = 8 ) yquad(maxquad,maxelm)
!
!  Set the weights and abscissas, depending on NQUAD.
!
  if ( nquad == 3 ) then

    wquad(1:3) = 1.0D+00 / 6.0D+00

    xsiquad(1) = 0.5D+00
    etaquad(1) = 0.5D+00

    xsiquad(2) = 1.0D+00
    etaquad(2) = 0.5D+00

    xsiquad(3) = 0.5D+00
    etaquad(3) = 0.0D+00

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SETQUAD - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of NQUAD = ', nquad
    write ( *, '(a)' ) '  Legal values are 3.'
    stop

  end if
!
!  Set the X, Y coordinates of quadrature points for each element.
!
  do ielem = 1, nelem
 
    do iquad = 1, nquad

      xsi = xsiquad(iquad)
      eta = etaquad(iquad)

      call xofxsi ( eta, ielem, maxelm, maxnp, node, x, xc, xsi, y, yc )

      xquad(iquad,ielem) = x
      yquad(iquad,ielem) = y

    end do
!
!  We only calculate true areas for nonisoparametric elements.
!
    ip1 = node(1,ielem)
    ip2 = node(2,ielem)
    ip3 = node(3,ielem)

    do iquad = 1, nquad

      if ( isotri(ielem) == 0 ) then
 
        area(iquad,ielem) = wquad(iquad) * abs ( &
               (yc(ip1)+yc(ip2)) * (xc(ip2)-xc(ip1)) &
             + (yc(ip2)+yc(ip3)) * (xc(ip3)-xc(ip2)) &
             + (yc(ip3)+yc(ip1)) * (xc(ip1)-xc(ip3)) )
 
      else
 
        area(iquad,ielem) = wquad(iquad)
 
      end if
  
    end do

  end do
 
  return
end
subroutine setnod ( eqn, indx, isotri, maxelm, maxeqn, maxnp, neqn, node, np, &
  nx, ny, region )

!*****************************************************************************80
!
!! SETNOD assigns numbers to the nodes and elements.
!
!  Discussion:
!
!    It also decides which elements shall be isoparametric, (ISOTRI) and 
!    assigns six nodes to each, via the NODE array.  
!
!    It associates global unknown indices with each node (INDX), and
!    computes NEQN, the number of unknowns and of equations, and
!    compares that to the maximum allowed value, MAXEQN.
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = 2 ) EQN(MAXEQN).
!    EQN records the "type" of each equation that will be generated, 
!    and which is associated with an unknown.  Note that most 
!    boundary conditions do not result in an equation.  The current 
!    values are:
!
!    'U'  The horizontal momentum equation.
!    'UB' The condition U=0 applied at a node on the bump.
!    'UI' The condition U=UInflow(Y,Lambda) at the inflow.
!    'UW' The condition U=0 applied at a node on a fixed wall.
!
!    'V'  The vertical momentum equation.
!    'VB' The condition V=0 applied at a node on the bump.
!    'VI' The condition V=VInflow(Y,Lambda) at the inflow.
!    'VW' The condition V=0 applied at a node on a fixed wall.
!
!    'P'  The continuity equation.
!    'PB' The condition P=0 applied at (XMAX,YMAX).
!
!    Output, integer INDX(3,MAXNP).  
!
!    INDX(I,J) contains, for each node J, the index of U, V and P 
!    at that node, or 0 or a negative value.
! 
!    If K=INDX(I,J) is positive, then the value of the degree
!    of freedom is stored in the solution vector entry G(K),
!    and an equation will be generated to determine its value.
!
!    If INDX(I,J) is not positive, then no equation is
!    generated to determine for variable I at node J, either 
!    because the variable is specified in some other way, or 
!    because (in the case of pressure), there is no coefficient 
!    associated with that node.
!
!    Output, integer ISOTRI(MAXELM), 0/1, the element IS NOT/IS isoparametric.
! 
!    Input, integer MAXEQN, the maximum number of equations.
!
!    Input, integer NELEM, the number of elements.
!
!    Output, integer NEQN, the number of finite element equations.
!
!    Input, integer NODE(6,MAXELM), the nodes that make up each element.
!
!    Input, integer NP, the number of nodes.
!
!    Input, integer NX.
!
!    NX controls the spacing of nodes and elements in
!    the X direction.  There are 2*NX-1 nodes along various
!    lines in the X direction.
!
!    Roughly speaking, NX (or 2*NX) is the number of elements along
!    a line in the X direction.
! 
!    Input, integer NY.
!
!    NY controls the spacing of nodes and elements in
!    the Y direction.  There are 2*NY-1 nodes along various
!    lines in the Y direction.
!
!    Roughly speaking, NY (or 2*NY) is the number of elements along
!    a line in the Y direction.
! 
!    Input, character ( len = 20 ) REGION, the flow region, 
!    'CHANNEL' or 'CAVITY'.
!
  implicit none

  integer maxelm
  integer maxeqn
  integer maxnp

  character ( len = 2 ) eqn(maxeqn)
  integer icol
  integer icol2
  integer ielem
  integer indx(3,maxnp)
  integer ip
  integer irow
  integer irow2
  integer isotri(maxelm)
  integer neqn
  integer node(6,maxelm)
  integer np
  integer nx
  integer ny
  character ( len = 20 ) region
!
!  Compute the global node numbers that will be assigned to the 
!  beginning and ending of the bump.  These numbers are only used to 
!  determine which elements are isoparametric.
!
!
!  Consider each of the NP nodes, which logically lie in an MX by MY
!  rectangular array.  A pair of new elements must be generated every
!  time we reach a node that lies in an odd row and column, (except for
!  the top row, and last column, of course).  At every node, we
!  will have to decide how many equations to generate.
!
  ielem = 0
  neqn = 0
 
  do ip = 1, np
!
!  Determine the row and column of this node, and also whether each
!  of these quantities is odd or even.
!
    icol = ( ( ip - 1 ) / ( 2 * ny - 1 ) ) + 1
    irow = mod ( ( ip - 1 ), 2 * ny - 1 ) + 1
 
    icol2 = mod ( icol, 2 )
    irow2 = mod ( irow, 2 )
!
!  If both the row and the column are odd, and we're not in the last
!  column or top row, then we can define two new triangular elements 
!  based at the node.
!
!  Given the following arrangement of nodes, for instance:
!
!    05 10 15 20 25
!    04 09 14 19 24
!    03 08 13 18 23
!    02 07 12 17 22
!    01 06 11 16 21
!
!  when we arrive at node 13, we will define
!
!    element 7: (25, 13, 15, 19, 14, 20)
!    element 8: (13, 25, 23, 19, 24, 18)
!
    if ( irow2 == 1 .and. icol2 == 1 .and. icol /= 2 * nx - 1 .and. &
         irow /= 2 * ny - 1 ) then
 
      ielem = ielem + 1
 
      node(1,ielem) = ip + 2 * ( 2 * ny - 1 ) + 2
      node(2,ielem) = ip
      node(3,ielem) = ip + 2
      node(4,ielem) = ip + ( 2 * ny - 1 ) + 1
      node(5,ielem) = ip + 1
      node(6,ielem) = ip + ( 2 * ny - 1 ) + 2
 
      if ( region == 'CHANNEL' ) then

        isotri(ielem) = 0

      else if ( region == 'CAVITY' ) then

        isotri(ielem) = 0

      end if
 
      ielem = ielem + 1
 
      node(1,ielem) = ip
      node(2,ielem) = ip + 2 * ( 2 * ny - 1 ) + 2
      node(3,ielem) = ip + 2 * ( 2 * ny - 1 )
      node(4,ielem) = ip + ( 2 * ny - 1 ) + 1
      node(5,ielem) = ip + 2 * ( 2 * ny - 1 ) + 1
      node(6,ielem) = ip + ( 2 * ny - 1 )
 
      if ( region == 'CHANNEL' ) then
 
        isotri(ielem) = 0

      else if ( region == 'CAVITY' ) then

        isotri(ielem) = 0

      end if
 
    end if

    if ( maxeqn < neqn + 2 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SetNod - Fatal error!'
      write ( *, '(a)' ) '  Too many unknowns!'
      write ( *, '(a,i8)' ) '  Processing node IP=', ip
      write ( *, '(a,i8)' ) '  The maximum allowed is MAXEQN=', maxeqn
      write ( *, '(a,i8)' ) '  This problem requires NEQN=', neqn + 2
      stop
    end if
!
!  Now determine what equations to associate with this node.
!
    if ( region == 'CHANNEL' ) then
!
!  The node lies on the left hand inflow boundary.
!  The horizontal and vertical velocities are specified.
!
      if ( icol == 1 .and. 1 < irow .and. irow < 2 * ny-1 ) then
 
        neqn = neqn + 1
        indx(1,ip) = neqn
        eqn(neqn) = 'UI'
 
        neqn = neqn + 1
        indx(2,ip) = neqn
        eqn(neqn) = 'VI'
!
!  The node lies on the right hand boundary.
!  The horizontal velocity is an unknown, the vertical velocity is 
!  zero.
!
      else if ( icol == 2 * nx - 1 .and. 1 < irow .and. irow < 2 * ny - 1 ) then
 
        neqn = neqn + 1
        indx(1,ip) = neqn
        eqn(neqn) = 'U'
 
        neqn = neqn + 1
        indx(2,ip) = neqn
        eqn(neqn) = 'VW'
!
!  The node lies on a fixed wall.
!  The horizontal and vertical velocities are zero.
!
      else if ( icol == 1 .or. icol == 2 * nx - 1 .or. irow == 1 .or. &
        irow == 2 * ny - 1 ) then
 
        neqn = neqn + 1
        indx(1,ip) = neqn
        eqn(neqn) = 'UW'
 
        neqn = neqn + 1
        indx(2,ip) = neqn
        eqn(neqn) = 'VW'
!
!  The node is a normal interior node.
!  The horizontal and vertical velocities are unknown.
!
      else
 
        neqn = neqn + 1
        indx(1,ip) = neqn
        eqn(neqn) = 'U'
 
        neqn = neqn + 1
        indx(2,ip) = neqn
        eqn(neqn) = 'V'
 
      end if
!
!  CAVITY:
!
    else if ( region == 'CAVITY' ) then

      if ( irow == 2 * ny - 1 ) then
 
        neqn = neqn + 1
        indx(1,ip) = neqn
        eqn(neqn) = 'UI'
 
        neqn = neqn + 1
        indx(2,ip) = neqn
        eqn(neqn) = 'VI'

      else if ( icol == 1 .or. icol == 2 * nx - 1 .or. irow == 1 ) then

        neqn = neqn + 1
        indx(1,ip) = neqn
        eqn(neqn) = 'UW'
 
        neqn = neqn + 1
        indx(2,ip) = neqn
        eqn(neqn) = 'VW'

      else

        neqn = neqn + 1
        indx(1,ip) = neqn
        eqn(neqn) = 'U'
 
        neqn = neqn + 1
        indx(2,ip) = neqn
        eqn(neqn) = 'V'

      end if

    end if
!
!  On nodes in an odd row and column, add a pressure equation.
!
    if ( irow2 == 1 .and. icol2 == 1 ) then

      neqn = neqn + 1

      if ( maxeqn < neqn ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SetNod - Fatal error!'
        write ( *, '(a)' ) '  Too many unknowns!'
        write ( *, '(a,i8)' ) '  Processing node IP=', ip
        write ( *, '(a,i8)' ) '  The maximum allowed is MAXEQN=', maxeqn
        write ( *, '(a,i8)' ) '  This problem requires NEQN=', neqn
        stop
      end if

      indx(3,ip) = neqn
      eqn(neqn) = 'P'
    else
      indx(3,ip) = 0
    end if
 
  end do
!
!  The last equation, which is guaranteed to be a pressure equation,
!  is replaced by a pressure boundary condition, associated with
!  an unknown.  (Even though we know this pressure will be zero).
!
  eqn(neqn) = 'PB'
 
  return
end
subroutine setxy ( maxnp, nx, ny, xc, xrange, yc, yrange )

!*****************************************************************************80
!
!! SETXY sets the X and Y coordinates of the nodes.
!
!  Discussion:
!
!    The nodes are numbered from the lower left corner, up and to the
!    right.  For example:
!
!      5  10  15
!      4   9  14
!      3   8  13
!      2   7  12
!      1   6  11
!
!  Modified:
!
!    12 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NP, the number of nodes.
!
!    Input, integer NX.  NX controls the spacing of nodes and elements in
!    the X direction.  There are 2*NX-1 nodes along various
!    lines in the X direction.  Roughly speaking, NX (or 2*NX) is the number 
!    of elements along a line in the X direction.
! 
!    Input, integer NY.  NY controls the spacing of nodes and elements in
!    the Y direction.  There are 2*NY-1 nodes along various
!    lines in the Y direction.  Roughly speaking, NY (or 2*NY) is the number 
!    of elements along a line in the Y direction.
! 
!    Output, real ( kind = 8 ) XC(MAXNP), the X coordinates of the nodes.
!
!    Input, real ( kind = 8 ) XRANGE, the width of the region.
!
!    Output, real ( kind = 8 ) YC(MAXNP), the Y coordinates of the nodes.
!
!    Input, real ( kind = 8 ) YRANGE, the height of the region.
!
  implicit none

  integer maxnp

  integer i
  integer ip
  integer j
  integer nx
  integer ny
  real ( kind = 8 ) x
  real ( kind = 8 ) xc(maxnp)
  real ( kind = 8 ) xrange
  real ( kind = 8 ) y
  real ( kind = 8 ) yc(maxnp)
  real ( kind = 8 ) yrange

  ip = 0
 
  do i = 1, 2 * nx-1

    x = real ( i - 1, kind = 8 ) * xrange / real ( 2 * nx - 2, kind = 8 )

    do j = 1, 2*ny-1
 
      ip = ip + 1

      y = real ( j - 1, kind = 8 ) * yrange / real ( 2 * ny - 2, kind = 8 )
 
      xc(ip) = x
      yc(ip) = y
 
    end do
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
subroutine trans ( det, detadx, detady, dxsidx, dxsidy, eta, ielem, maxelm, &
  maxnp, node, xc, xsi, yc )

!*****************************************************************************80
!
!! TRANS calculates the biquadratic reference element transformation.
!
!  Discussion:
!
!    This transformation maps the reference element in (XSI,ETA) space 
!    into a particular isoparametric element in (X,Y) space.
!
!    We know everything about the isoparametric element once we
!    specify the location of its six nodes.
!
!    The routine computes the entries of the jacobian of the transformation
!    and the determinant of the jacobian.  Essentially, the jacobian
!    records the relationship between derivatives with respect to XSI
!    and ETA and a point in the reference element, and derivatives
!    with respect to X and Y of the same function as defined in the
!    isoparametric element.
!
!    The four entries of the jacobian are symbolically named DETADX,
!    DETADY, DXSIDX and DXSIDY, and we know that the jacobian gives
!    us the following relation between derivatives with respect to
!    XSI and ETA, and derivatives with respect to X and Y:
!
!      d F(X,Y)/dX     (d XSI/dX  d ETA/dX )   ( d F(XSI, ETA)/d XSI )
!      d F(X,Y)/dY  =  (d XSI/dY  d ETA/dY ) * ( d F(XSI, ETA)/d ETA )
!
!    Here is a graph of the (XSI, ETA) reference triangle.
!
!          ^
!          |
!        1 +        2
!          |       /|
!    ETA   |      4 5
!          |     /  |
!        0 +    1-6-3
!          |
!          +----+---+--->
!               0   1 
! 
!              XSI
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) DET, the determinant of the 
!    jacobian of the transformation from the reference to the
!    (physical) isoparametric element.
!
!    Output, real ( kind = 8 ) DETADX, DETADY, the partial 
!    derivative d ETA/d X and d ETA/d Y at (XSI,ETA).
!
!    Output, real ( kind = 8 ) DXSIDX, DXSIDY, the partial 
!    derivative d XSI/d X and d XSI/d Y at (XSI,ETA).
!
!    Input, real ( kind = 8 ) ETA, the ETA coordinate of the point.
!
!    Input, integer IELEM, the index of the element.
!
!    Input, integer NELEM, the number of elements.
!
!    Input, integer NODE(6,MAXELM), the nodes that make up each element.
!
!    Input, integer NP, the number of nodes.
!
!    Input, real ( kind = 8 ) XC(MAXNP), the X coordinates of the nodes.
!
!    Input, real ( kind = 8 ) XSI, the XSI coordinate of the point.
!
!    Input, real ( kind = 8 ) YC(MAXNP), the Y coordinates of the nodes.
!
  implicit none

  integer maxelm
  integer maxnp

  real ( kind = 8 ) a1
  real ( kind = 8 ) a2
  real ( kind = 8 ) b1
  real ( kind = 8 ) b2
  real ( kind = 8 ) c1
  real ( kind = 8 ) c2
  real ( kind = 8 ) d1
  real ( kind = 8 ) d2
  real ( kind = 8 ) det
  real ( kind = 8 ) detadx
  real ( kind = 8 ) detady
  real ( kind = 8 ) dxdeta
  real ( kind = 8 ) dxdxsi
  real ( kind = 8 ) dxsidx
  real ( kind = 8 ) dxsidy
  real ( kind = 8 ) dydeta
  real ( kind = 8 ) dydxsi
  real ( kind = 8 ) e1
  real ( kind = 8 ) e2
  real ( kind = 8 ) eta
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer i
  integer ielem
  integer node(6,maxelm)
  real ( kind = 8 ) x
  real ( kind = 8 ) xn(6)
  real ( kind = 8 ) xc(maxnp)
  real ( kind = 8 ) xsi
  real ( kind = 8 ) y
  real ( kind = 8 ) yn(6)
  real ( kind = 8 ) yc(maxnp)
!
!  Pick off the X, Y coordinates of the nodes and store them
!  in two short lists.
!
  do i = 1, 6
    xn(i) = xc ( node(i,ielem) )
    yn(i) = yc ( node(i,ielem) )
  end do
!
!  Set the coefficients in the transformation
!
!    (XSI,ETA) --> (X,Y).
!
!  The mapping has the form:
!
!    X(XSI,ETA) = A1 * XSI**2 + B1 * XSI*ETA + C1 * ETA**2
!               + D1 * XSI    + E1 * ETA     + F1
!
!    Y(XSI,ETA) = A2 * XSI**2 + B2 * XSI*ETA + C2 * ETA**2
!               + D2 * XSI    + E2 * ETA     + F2
!
  a1 =  2.0D+00 *xn(1)+2.0D+00 *xn(3)-4.0D+00*xn(6)
  b1 = -4.0D+00 *xn(3)-4.0D+00 *xn(4)+4.0D+00*xn(5)+4.0D+00*xn(6)
  c1 =  2.0D+00 *xn(2)+2.0D+00 *xn(3)-4.0D+00*xn(5)
  d1 = -3.0*xn(1)    -xn(3)+4.0D+00*xn(6)
  e1 =     -xn(2)    +xn(3)+4.0D+00*xn(4)-4.0D+00 *xn(6)
  f1 =      xn(1)
 
  a2 =  2.0D+00 *yn(1)+2.0D+00 *yn(3)-4.0D+00*yn(6)
  b2 = -4.0D+00*yn(3)-4.0D+00 *yn(4)+4.0D+00*yn(5)+4.0D+00*yn(6)
  c2 =  2.0D+00 *yn(2)+2.0D+00 *yn(3)-4.0D+00*yn(5)
  d2 = -3.0*yn(1)    -yn(3)+4.0D+00*yn(6)
  e2 =     -yn(2)    +yn(3)+4.0D+00*yn(4)-4.0D+00 *yn(6)
  f2 =      yn(1)
!
!  Compute the partial derivatives at the point (XSI,ETA).
!  This is the jacobian matrix
!
!    J: (XSI,ETA) --> (X,Y).
!
  dxdxsi = 2.0D+00 *a1*xsi +     b1*eta + d1
  dxdeta =     b1*xsi + 2.0D+00 *c1*eta + e1
 
  dydxsi = 2.0D+00 *a2*xsi +     b2*eta + d2
  dydeta =     b2*xsi + 2.0D+00 *c2*eta + e2
!
!  Compute the determinant of the jacobian matrix:
!
!    J: (XSI,ETA) --> (X,Y)
!
  det = dxdxsi * dydeta - dxdeta * dydxsi
!
!  Watch out for a zero determinant.
!
  if ( det == 0.0D+00 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Trans - Fatal error!'
    write ( *, '(a)' ) '  J: (XSI,ETA) --> (X,Y) is singular!'
    write ( *, '(a,i8)' ) '  This occurred for element number ', ielem
    write ( *, '(a,2g14.66)' ) '  Local coordinates XSI,ETA=', xsi, eta

    x = a1*xsi**2 + b1*xsi*eta + c1*eta**2 + d1*xsi + e1*eta + f1
    y = a2*xsi**2 + b2*xsi*eta + c2*eta**2 + d2*xsi + e2*eta + f2

    write ( *, '(a,2g14.6)' ) '  Global coordinates X,Y=', x, y
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The X, Y nodes were:'
    write ( *, '(a)' ) ' '

    do i = 1, 6
      write ( *, '(2g14.6)' ) xn(i), yn(i)
    end do
 
    stop
  end if
!
!  Compute
!
!    d ETA/d X, d ETA/d Y, d XSI/d X, d XSI/d Y
!
!  by inverting the jacobian matrix
!
!    J: (XSI,ETA) --> (X,Y)
!
!  to get the jacobian matrix
!
!    J: (X,Y) --> (XSI,ETA).
!
!  This uses the simple fact that the inverse of
!
!    (a b)
!    (c d)
!
!  is
!
!    1/(ad-bc) * ( d -b)
!                (-c  a)
!
  dxsidx =  dydeta / det
  dxsidy = -dxdeta / det
 
  detadx = -dydxsi / det
  detady =  dxdxsi / det
 
  return
end
subroutine uvpnorm ( g, indx, maxeqn, maxnp, np, pnorm, uvnorm )

!*****************************************************************************80
!
!! UVPNORM returns the "norm" of the computed finite element solution.  
!
!  Definition:
!
!    We define the norm of a solution G = (U,V,P) to be two numbers:
!
!    UVNORM = the maximum velocity magnitude at a node,
!    PNORM = the maximum pressure at a node.
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!  
!    Input, real ( kind = 8 ) G(MAXEQN).
!
!    G is the current solution vector, in which are stored 
!    the finite element coefficients that define the velocity
!    and pressure functions, U, V and P.
!
!    Input, integer INDX(3,MAXNP).  
!
!    INDX(I,J) contains, for each node J, the index of U, V and P 
!    at that node, or 0 or a negative value.
! 
!    If K=INDX(I,J) is positive, then the value of the degree
!    of freedom is stored in the solution vector entry G(K),
!    and an equation will be generated to determine its value.
!
!    If INDX(I,J) is not positive, then no equation is
!    generated to determine for variable I at node J, either 
!    because the variable is specified in some other way, or 
!    because (in the case of pressure), there is no coefficient 
!    associated with that node.
!
!    Input, integer NEQN, the number of finite element equations.
! 
!    Input, integer NP, the number of nodes.
!
!    Output, real ( kind = 8 ) PNORM, the maximum absolute value
!    pressure coefficient.
!
!    Output, real ( kind = 8 ) UVNORM, the maximum velocity magnitude.
!
  implicit none

  integer maxeqn
  integer maxnp

  real ( kind = 8 ) g(maxeqn)
  integer i
  integer indx(3,maxnp)
  integer np
  real ( kind = 8 ) p
  real ( kind = 8 ) pnorm
  real ( kind = 8 ) u
  real ( kind = 8 ) uvnorm
  real ( kind = 8 ) v

  uvnorm = 0.0D+00
  pnorm = 0.0D+00

  do i = 1, np

    u = g ( indx(1,i) )
    v = g ( indx(2,i) )
    uvnorm = max ( uvnorm, sqrt ( u**2 + v**2 ) )

    if ( 0 < indx(3,i) ) then
      p = g ( indx(3,i) )
      pnorm = max ( pnorm, abs ( p ) )
    end if

  end do

  return
end
subroutine uvalq ( dpdx, dpdy, dudx, dudy, dvdx, dvdy, g, ielem, indx, iquad, &
  maxelm, maxeqn, maxnp, maxquad, node, p, phi, u, v )

!*****************************************************************************80
!
!! UVALQ evaluates velocities and pressure at a quadrature point in an element.
!
!  Discussion:
!
!    The X and Y derivatives of the quantities are also evaluated.
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) DPDX, DPDY, the derivatives of the
!    pressure function with respect to X and Y.
!
!    Output, real ( kind = 8 ) DUDX, DUDY, the derivatives of the
!    horizontal velocity function with respect to X and Y.
!
!    Output, real ( kind = 8 ) DVDX, DVDY, the derivatives of the
!    vertical velocity function with respect to X and Y.
!
!    Input, real ( kind = 8 ) G(MAXEQN), the current solution vector, in 
!    which are stored the finite element coefficients that define the 
!    velocity and pressure functions, U, V and P.
!
!    Input, integer IELEM, the index of the element.
!
!    Input, integer INDX(3,MAXNP).  
!
!    INDX(I,J) contains, for each node J, the index of U, V and P 
!    at that node, or 0 or a negative value.
! 
!    If K=INDX(I,J) is positive, then the value of the degree
!    of freedom is stored in the solution vector entry G(K),
!    and an equation will be generated to determine its value.
!
!    If INDX(I,J) is not positive, then no equation is
!    generated to determine for variable I at node J, either 
!    because the variable is specified in some other way, or 
!    because (in the case of pressure), there is no coefficient 
!    associated with that node.
!
!    Input, integer IQUAD, the index of the quadrature point.
!
!    Input, integer NELEM, the number of elements.
! 
!    Input, integer NEQN, the number of finite element equations.
! 
!    Input, integer NODE(6,MAXELM), the nodes that make up each element.
! 
!    Input, integer NP, the number of nodes.  NP=(2*NX-1)*(2*NY-1).
! 
!    Output, real ( kind = 8 ) P, the value of the pressure.
!
!    Input, real ( kind = 8 ) PHI(MAXQUAD,6,10,MAXELM).  
!
!    PHI contains the value of a basis function, its derivative,
!    or other information, evaluated at a quadrature point.
! 
!    For a particular element I, quadrature point J, and basis
!    function K, we use the following shorthand for the 
!    entries of PHI:
!
!      W, dWdX, dWdY
!      Q, dQdX, dQdY
!      dXsidX, dXsidY, dEtadX, dEtadY
!
!    W is the quadratic basis function associated with velocity,
!    Q the linear basis function associated with pressure,
!    Xsi and Eta the reference coordinates for the point.
!        
!    In particular, PHI(J,K,1,I) is the value of the quadratic 
!    basis function associated with local node K in element I, 
!    evaluated at quadrature point J.
! 
!    Note that PHI(J,K,4,I)=PHI(J,K,5,I)=PHI(J,K,6,I)=0 for 
!    K=4, 5, or 6, since there are only three linear basis
!    functions.
!
!    Output, real ( kind = 8 ) U, the value of the horizontal velocity.
!
!    Output, real ( kind = 8 ) V, the value of the vertical velocity.
!
  implicit none

  integer maxelm
  integer maxeqn
  integer maxnp
  integer maxquad

  real ( kind = 8 ) coef
  real ( kind = 8 ) dpdx
  real ( kind = 8 ) dpdy
  real ( kind = 8 ) dqdx
  real ( kind = 8 ) dqdy
  real ( kind = 8 ) dudx
  real ( kind = 8 ) dudy
  real ( kind = 8 ) dvdx
  real ( kind = 8 ) dvdy
  real ( kind = 8 ) dwdx
  real ( kind = 8 ) dwdy
  real ( kind = 8 ) g(maxeqn)
  integer ielem
  integer indx(3,maxnp)
  integer ip
  integer iq
  integer iquad
  integer iun
  integer node(6,maxelm)
  real ( kind = 8 ) p
  real ( kind = 8 ) phi(maxquad,6,10,maxelm)
  real ( kind = 8 ) q
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) w
!
!  Start all the functions at zero.
!
  p = 0.0D+00
  u = 0.0D+00
  v = 0.0D+00
  dpdx = 0.0D+00
  dpdy = 0.0D+00
  dudx = 0.0D+00
  dudy = 0.0D+00
  dvdx = 0.0D+00
  dvdy = 0.0D+00
!
!  Now each of these functions is represented as the sum of
!  coefficients times basis functions.  In this particular
!  element, at this particular quadrature point, we know that
!  exactly 6 basis functions are nonzero.  So if
!  we simply look up the values of the basis functions (and
!  their X and Y derivatives), and multiply by the appropriate
!  coefficients, we can evaluate the functions.
!
!  W, DWDX and DWDY represent the value of a quadratic basis
!  function and its X and Y derivative.
!
!  Q, DQDX and DQDY represent the value of a linear basis
!  function and its X and Y derivatives.
!
  do iq = 1, 6
 
    w = phi(iquad,iq,1,ielem)
    dwdx = phi(iquad,iq,2,ielem)
    dwdy = phi(iquad,iq,3,ielem)
 
    q = phi(iquad,iq,4,ielem)
    dqdx = phi(iquad,iq,5,ielem)
    dqdy = phi(iquad,iq,6,ielem)
!
!  Now that we have the basis function values, we need to look
!  up the coefficient COEF that multiplies the basis function.
!
    ip = node(iq,ielem)

    iun = indx(1,ip)
    coef = g(iun)
    u = u + coef * w
    dudx = dudx + coef * dwdx
    dudy = dudy + coef * dwdy

    iun = indx(2,ip)
    coef = g(iun)
    v = v + coef * w
    dvdx = dvdx + coef * dwdx
    dvdy = dvdy + coef * dwdy

    iun = indx(3,ip)

    if ( 0 < iun ) then
      coef = g(iun)
      p = p + coef * q
      dpdx = dpdx + coef * dqdx
      dpdy = dpdy + coef * dqdy
    end if
 
  end do
 
  return
end
subroutine write_display_file ( ierror, maxelm, maxnp, nelem, node, np, p, &
  u, v, xc, yc )

!*****************************************************************************80
!
!! WRITE_DISPLAY_FILE writes solution information to two files.
!
!  Discussion:
!
!    The data may be used to plot the finite element mesh and solution.
!
!    The first file, "elements.txt", contains the NODE array.
!    The second file, "uvp.txt", contains X, Y, U, V, and P at each node.  
!
!  Modified:
!
!    03 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer IERROR, error flag.
!    0, no error occurred.
!    nonzero, an error occurred.
!
!    Input, integer NELEM, the number of elements.
!
!    Input, integer NODE(6,MAXELM), the nodes that make up each element.
!
!    Input, integer NP, the number of nodes.  NP=(2*NX-1)*(2*NY-1).
!
!    Input, real P(MAXNP), U(MAXNP), V(MAXNP), the pressure, horizontal and
!    vertical velocity at each node.
!
!    Input, real XC(MAXNP), YC(MAXNP), the X and Y coordinates of the nodes.
!
  implicit none

  integer maxelm
  integer maxnp

  integer i
  integer ielem
  integer ierror
  integer ios
  integer nelem
  integer node(6,maxelm)
  integer np
  real ( kind = 8 ) p(maxnp)
  real ( kind = 8 ) u(maxnp)
  real ( kind = 8 ) v(maxnp)
  real ( kind = 8 ) xc(maxnp)
  real ( kind = 8 ) yc(maxnp)

  ierror = 0
!
!  Open the element data file.
!
  open ( unit = 2, file = 'elements.txt', status = 'replace', iostat = ios )

  if ( ios /= 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WRITE_DISPLAY_FILE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the element file!'
    return
  end if

  do ielem = 1, nelem
    write ( 2, '(6i8)' ) node(1:6,ielem)
  end do

  close ( unit = 2 )
!
!  Open the node data file.
!
  open ( unit = 2, file = 'uvp.txt', status = 'replace', iostat = ios )

  if ( ios /= 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WRITE_DISPLAY_FILE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the node file.'
    return
  end if

  do i = 1, np
    write ( 2, '(2f12.6,3g14.6)' ) xc(i), yc(i), u(i), v(i), p(i)
  end do

  close ( unit = 2 )
  
  return
end      
subroutine write_tecplot_file ( maxelm, maxnp, nelem, node, np, p, u, v, xc, yc )

!*****************************************************************************80
!
!! WRITE_TECPLOT_FILE writes out solution information for use with TECPLOT.
!
!  Modified:
!
!    28 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NELEM, the number of elements.
! 
!    Input, integer NODE(6,MAXELM), the nodes that make up each element.
! 
!    Input, integer NP, the number of nodes.
!
!    Input, real P(MAXNP), the pressure.
!
!    Input, real U(MAXNP), the horizontal velocity.
!
!    Input, real V(MAXNP), the vertical velocity.
!
!    Input, real XC(MAXNP), the X coordinates of the nodes.
!
!    Input, real YC(MAXNP), the Y coordinates of the nodes.
!
  implicit none

  integer maxelm
  integer maxnp

  integer i
  integer ielem
  integer nelem
  integer node(6,maxelm)
  integer np
  real ( kind = 8 ) p(maxnp)
  real ( kind = 8 ) u(maxnp)
  real ( kind = 8 ) v(maxnp)
  real ( kind = 8 ) xc(maxnp)
  real ( kind = 8 ) yc(maxnp)

  open ( unit = 10, file = 'flow5.tec', form = 'formatted', &
    access = 'sequential', status = 'replace' )

  write ( 10, * ) 'Title="FLOW5 data"'
  write ( 10, * ) 'Variables="X","Y","P","U","V"'
  write ( 10, * ) 'Zone N=', np, ', E=', 4*nelem, ', F=FEPOINT, ET=TRIANGLE'
!
!  Write out the data at each node.
!
  do i = 1, np
    write ( 10, '(5g15.6)' ) xc(i), yc(i), p(i), u(i), v(i)
  end do
!
!  Write out the data that defines the elements.
!  Each 6 node quadratic element must be described as 4 linear elements.
!
  do ielem = 1, nelem
    write ( 10, '(3i8)' ) node(1,ielem), node(4,ielem), node(6,ielem)
    write ( 10, '(3i8)' ) node(2,ielem), node(5,ielem), node(4,ielem)
    write ( 10, '(3i8)' ) node(3,ielem), node(6,ielem), node(5,ielem)
    write ( 10, '(3i8)' ) node(4,ielem), node(5,ielem), node(6,ielem)
  end do

  close ( unit = 10 )

  return
end
subroutine xofxsi ( eta, ielem, maxelm, maxnp, node, x, xc, xsi, y, yc )

!*****************************************************************************80
!
!! XOFXSI computes X and Y given XSI and ETA coordinates.
!
!  Diagram:
!
!    Here is a graph of the (XSI, ETA) reference triangle.
!
!          ^
!          |
!        1 +        2
!          |       /|
!    ETA   |      4 5
!          |     /  |
!        0 +    1-6-3
!          |
!          +----+---+--->
!               0   1
!
!              XSI
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ETA, the ETA coordinate of the point.
!
!    Input, integer IELEM, the index of the element.
!
!    Input, integer NELEM, the number of elements.
!
!    Input, integer NODE(6,MAXELM), the nodes that make up each element.
!
!    Input, integer NP, the number of nodes.
!
!    Output, real ( kind = 8 ) X, the X coordinate of the point.
!
!    Input, real ( kind = 8 ) XC(MAXNP), the X coordinates of the nodes.
!
!    Input, real ( kind = 8 ) XSI, the XSI coordinate of the point.
!
!    Output, real ( kind = 8 ) Y, the Y coordinate of the point.
!
!    Input, real ( kind = 8 ) YC(MAXNP), the Y coordinates of the nodes.
!
  implicit none

  integer maxelm
  integer maxnp

  real ( kind = 8 ) a1
  real ( kind = 8 ) a2
  real ( kind = 8 ) b1
  real ( kind = 8 ) b2
  real ( kind = 8 ) c1
  real ( kind = 8 ) c2
  real ( kind = 8 ) d1
  real ( kind = 8 ) d2
  real ( kind = 8 ) e1
  real ( kind = 8 ) e2
  real ( kind = 8 ) eta
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer i
  integer ielem
  integer node(6,maxelm)
  real ( kind = 8 ) x
  real ( kind = 8 ) xn(6)
  real ( kind = 8 ) xc(maxnp)
  real ( kind = 8 ) xsi
  real ( kind = 8 ) y
  real ( kind = 8 ) yn(6)
  real ( kind = 8 ) yc(maxnp)
!
!  Pick off the X, Y coordinates of the nodes and store them
!  in two short lists.
!
  do i = 1, 6
    xn(i) = xc ( node(i,ielem) )
    yn(i) = yc ( node(i,ielem) )
  end do
!
!  Set the coefficients in the transformation
!
!    (XSI,ETA) --> (X,Y).
!
!  The mapping has the form:
!
!    X(ETA,XSI) = A1 * XSI**2 + B1 * XSI*ETA + C1 * ETA**2
!               + D1 * XSI    + E1 * ETA     + F1
!
!    Y(ETA,XSI) = A2 * XSI**2 + B2 * XSI*ETA + C2 * ETA**2
!               + D2 * XSI    + E2 * ETA     + F2
!
  a1 =  2.0D+00 *xn(1)+2.0D+00 *xn(3)-4.0D+00*xn(6)
  b1 = -4.0D+00*xn(3)-4.0D+00 *xn(4)+4.0D+00 * xn(5)+4.0D+00 * xn(6)
  c1 =  2.0D+00 *xn(2)+2.0D+00 *xn(3)-4.0D+00*xn(5)
  d1 = -3.0D+00 * xn(1)      -xn(3)+4.0D+00*xn(6)
  e1 =       -xn(2)      +xn(3)+4.0D+00*xn(4)-4.0D+00 *xn(6)
  f1 =        xn(1)
 
  a2 =  2.0D+00 *yn(1)+2.0D+00 *yn(3)-4.0D+00*yn(6)
  b2 = -4.0D+00*yn(3)-4.0D+00 *yn(4)+4.0D+00 * yn(5)+4.0D+00 * yn(6)
  c2 =  2.0D+00 * yn(2) + 2.0D+00 * yn(3) - 4.0D+00 * yn(5)
  d2 = -3.0D+00 * yn(1)      -yn(3)+4.0D+00*yn(6)
  e2 =       -yn(2)      + yn(3) + 4.0D+00 * yn(4)-4.0D+00 *yn(6)
  f2 =        yn(1)
 
  x = a1*xsi**2 + b1*xsi*eta + c1*eta**2 + d1*xsi + e1*eta + f1
 
  y = a2*xsi**2 + b2*xsi*eta + c2*eta**2 + d2*xsi + e2*eta + f2
 
  return
end
