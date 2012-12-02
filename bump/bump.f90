program main

!*****************************************************************************80
!
!! MAIN is the main program for BUMP.
!
!  Discussion:
!
!    BUMP is the main program for the flow-with-bump code.
!
!    The Navier Stokes equations:
!
!    The primitive variable formulation involves horizontal velocity U,
!    vertical velocity V, and pressure P.  The equations are:
!
!      U dUdx + V dUdy + dPdx - mu*(ddU/dxdx + ddU/dydy) = F1
!
!      U dVdx + V dVdy + dPdy - mu*(ddV/dxdx + ddV/dydy) = F2
!
!      dUdx + dVdy = 0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    John Burkardt, Max Gunzburger, Janet Peterson,
!    Discretization of Cost and Sensitivities in Shape Optimization,
!    in Computation and Control IV,
!    edited by Bowers and Lund,
!    Proceedings of the Fourth Bozeman Conference on Computation
!    and Control,
!    Birkhaeuser, 1995.
!
!  Local Parameters:
!
!    Local, real ( kind = 8 ) A(MAXROW,NEQN), contains the matrix
!    in LINPACK general band storage mode.
!
!    Local, real ( kind = 8 ) APROF, the value of the parameter for
!    which the target profile was generated.
!
!    Local, real ( kind = 8 ) AREA(NELEMN), the area of the elements.
!
!    Local, real ( kind = 8 ) DCDA(MY), the sensitivities.
!
!    Local, real ( kind = 8 ) F(NEQN), is used as a right hand side
!    vector in LINSYS, and is overwritten there by the
!    new solution (which is then copied into g).
!
!    Local, real ( kind = 8 ) G(NEQN), the current solution vector,
!    in which are stored pressures and velocities.
!
!    Local, integer ( kind = 4 ) IBUMP, determines where isoparametric elements
!    will be used.
!    0, no isoparametric elements will be used.
!    1, isoparametric elements will be used only for the
!       elements which directly impinge on the bump.
!    2, isoparametric elements will be used for all elements which
!       are above the bump.
!    3, isoparametric elements will be used for all elements.
!
!    Local, integer ( kind = 4 ) INDX(MAXNP,2), contains, for each node I, the
!    index of the U and V velocities at that node, or 0.
!
!    Local, integer ( kind = 4 ) INSC(NP), contains, for each node I, the
!    index of the pressure P at that node, or 0.
!
!    Local, integer ( kind = 4 ) ISOTRI(NELEMN), contains, for each element,
!    a 1 if the element is isoparametric, or 0 otherwise.
!
!    Local, integer ( kind = 4 ) IWRITE, controls the amount of output printed.
!
!    Local, logical LONG.
!    .TRUE. if region is "long and thin", and
!    .FALSE. if region is "tall and skinny".
!
!    Local, integer ( kind = 4 ) MAXELM, the maximum number of elements.
!    This is actually simply the first dimension of the array NODE.
!
!    Local, integer ( kind = 4 ) MAXEQN, the maximum number of equations,
!    and functions.
!
!    Local, integer ( kind = 4 ) MAXNEW, the maximum number of Newton steps to take
!    in one iteration.
!
!    Local, integer ( kind = 4 ) MAXNP, the maximum number of nodes.
!    This is actually simply the first dimension of the array INDX.
!
!    Local, integer ( kind = 4 ) MAXROW, the first dimension of the matrix A.
!
!    Local, integer ( kind = 4 ) MAXSEC, the maximum number of secant steps to take.
!
!    Local, integer ( kind = 4 ) MX, the number of nodes in the X direction.
!
!    Local, integer ( kind = 4 ) MY, the number of nodes in the Y direction.
!
!    Local, integer ( kind = 4 ) NBAND, the bandwidth of the linear system.
!
!    Local, integer ( kind = 4 ) NBLEFT, the column at which the left corner of the
!    bump lies.
!
!    Local, integer ( kind = 4 ) NBRITE, the column at which the right corner of the
!    bump lies.
!
!    Local, integer ( kind = 4 ) NELEMN, the number of elements.
!
!    Local, integer ( kind = 4 ) NEQN, the number of equations, and functions.
!
!    Local, integer ( kind = 4 ) NLBAND, the lower bandwidth of the matrix.
!
!    Local, integer ( kind = 4 ) NNODES, the number of nodes per element, which
!    is 6 for the quadratic triangular elements in use here.
!
!    Local, integer ( kind = 4 ) NODE(MAXELM,NNODES), contains, for each element,
!    the global node indices of the element nodes.
!
!    Local, integer ( kind = 4 ) NODEX0, the node whose Y coordinate is zero,
!    and whose X coordinate is XPROF.  This is the first
!    node along the line where the velocity profile is
!    measured.
!
!    Local, integer ( kind = 4 ) NP, the number of nodes.
!
!    Local, integer ( kind = 4 ) NQUAD, the number of quadrature points per
!    element.
!
!    Local, integer ( kind = 4 ) NROW, the number of rows need to store the
!    matrix A.
!
!    Local, integer ( kind = 4 ) NUMNEW, total number of Newton steps taken.
!
!    Local, integer ( kind = 4 ) NX, controls the spacing of nodes and elements in
!    the X direction.  There are 2*NX+1 nodes in the X
!    direction.
!
!    Local, integer ( kind = 4 ) NY, controls the spacing of nodes and elements in
!    the Y direction.  There are 2*NY+1 nodes in the Y
!    direction.
!
!    Local, real ( kind = 8 ) PHI(NELEMN,NQUAD,NNODES,3).  Each entry of PHI
!    contains the value of a quadratic basis function or its derivative,
!    evaluated at a quadrature point.
!    In particular, PHI(I,J,K,1) is the value of the quadratic basis
!    function associated with local node K in element I, evaluatated
!    at quadrature point J.
!    PHI(I,J,K,2) is the X derivative of that same basis function,
!    PHI(I,J,K,3) is the Y derivative of that same basis function.
!
!    Local, real ( kind = 8 ) PSI(NELEMN,NQUAD,NNODES).  Each entry of PSI
!    contains the value of a linear basis function evaluated at a
!    quadrature point.
!    PSI(I,J,K) is the value of the linear basis function associated
!    with local node K in element I, evaluated at quadrature point J.
!
!    Local, real ( kind = 8 ) RES(NEQN), holds the residual.
!
!    Local, real ( kind = 8 ) REYNLD, the value of the Reynolds
!    number.
!
!    Local, real ( kind = 8 ) SENS(NEQN), the sensitivities.
!
!    Local, real ( kind = 8 ) TOLNEW, convergence tolerance for Newton
!    iteration.
!
!    Local, real ( kind = 8 ) TOLSEC, convergence tolerance for secant
!    iteration.
!
!    Local, real ( kind = 8 ) UPROF(MY), the horizontal velocity
!    along the profile line.
!
!    Local, real ( kind = 8 ) XBLEFT, the left X coordinate of the bump.
!
!    Local, real ( kind = 8 ) XBRITE, the right X coordinate of the bump.
!
!    Local, real ( kind = 8 ) XC(NP), the X coordinates of the nodes.
!
!    Local, real ( kind = 8 ) XLNGTH, the length of the region.
!
!    Local, real ( kind = 8 ) XM(NELEMN,NQUAD), contains the X coordinates
!    of the quadrature points for each element.
!
!    Local, real ( kind = 8 ) XPROF, the X coordinate at which the
!    profile is measured.  This value should be a grid value!
!
!    Local, integer ( kind = 4 ) XY_UNIT, the unit to which XY graphics is written.
!
!    Local, real ( kind = 8 ) YC(NP), the Y coordinates of the nodes.
!
!    Local, real ( kind = 8 ) YLNGTH, the height of the region.
!
!    Local, real ( kind = 8 ) YM(NELEMN,NQUAD), contains the Y coordinates
!    of the quadrature points for each element.
!
  implicit none

  integer ( kind = 4 ), parameter :: maxnew = 4
  integer ( kind = 4 ), parameter :: maxsec = 10
  integer ( kind = 4 ), parameter :: nx = 21
  integer ( kind = 4 ), parameter :: ny = 7
!
!  This assignment should really read (maxrow = 27*min(nx,ny))
!
  integer ( kind = 4 ), parameter :: maxrow = 27*ny
  integer ( kind = 4 ), parameter :: nelemn = 2*(nx-1)*(ny-1)
  integer ( kind = 4 ), parameter :: mx = 2*nx-1
  integer ( kind = 4 ), parameter :: my = 2*ny-1
  integer ( kind = 4 ), parameter :: maxeqn = 2*mx*my+nx*ny
  integer ( kind = 4 ), parameter :: np = mx*my
  integer ( kind = 4 ), parameter :: nnodes = 6
  integer ( kind = 4 ), parameter :: nquad = 3

  real ( kind = 8 ) a(maxrow,maxeqn)
  real ( kind = 8 ) anew
  real ( kind = 8 ) anext
  real ( kind = 8 ) aold
  real ( kind = 8 ) aprof
  real ( kind = 8 ) area(nelemn)
  real ( kind = 8 ) cpu1
  real ( kind = 8 ) cpu2
  real ( kind = 8 ) dcda(my)
  real ( kind = 8 ), external :: ddot
  real ( kind = 8 ) f(maxeqn)
  real ( kind = 8 ) g(maxeqn)
  real ( kind = 8 ) gr(my,my)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibump
  integer ( kind = 4 ) iline(my)
  integer ( kind = 4 ) indx(np,2)
  integer ( kind = 4 ) insc(np)
  integer ( kind = 4 ) isotri(nelemn)
  integer ( kind = 4 ) iter
  integer ( kind = 4 ) itype
  integer ( kind = 4 ) iukk
  integer ( kind = 4 ) iwrite
  integer ( kind = 4 ) j
  logical long
  integer ( kind = 4 ) nband
  integer ( kind = 4 ) neqn
  integer ( kind = 4 ) nlband
  integer ( kind = 4 ) node(nelemn,nnodes)
  integer ( kind = 4 ) npara
  integer ( kind = 4 ) nrow
  integer ( kind = 4 ) numnew
  integer ( kind = 4 ) numsec
  real ( kind = 8 ) para
  real ( kind = 8 ) phi(nelemn,nquad,nnodes,3)
  real ( kind = 8 ) psi(nelemn,nquad,nnodes)
  real ( kind = 8 ) r(my)
  real ( kind = 8 ) res(maxeqn)
  real ( kind = 8 ) reynld
  real ( kind = 8 ) rjpnew
  real ( kind = 8 ) rjpold
  real ( kind = 8 ) rtemp(my)
  real ( kind = 8 ) sens(maxeqn)
  real tarray(2)
  real ( kind = 8 ) temp
  real ( kind = 8 ) test
  real ( kind = 8 ) tolnew
  real ( kind = 8 ) tolsec
  real ( kind = 8 ) uprof(my)
  character ( len = 80 ) uv_file
  integer ( kind = 4 ) uv_unit
  real ( kind = 8 ) xbleft
  real ( kind = 8 ) xbrite
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) xlngth
  real ( kind = 8 ) xm(nelemn,nquad)
  real ( kind = 8 ) xprof
  character ( len = 80 ) xy_file
  integer ( kind = 4 ) xy_unit
  real ( kind = 8 ) yc(np)
  real ( kind = 8 ) ylngth
  real ( kind = 8 ) ym(nelemn,nquad)
  real ( kind = 8 ) ypert

  call cpu_time ( cpu1 )

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BUMP'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Control problem for channel flow over a bump.'
!
!  Set the input data.
!
  anew = 0.0D+00
  anext = 0.3
  aold = 0.0D+00
  aprof = 0.25D+00
  ibump = 2
  iwrite = 1
  numnew = 0
  numsec = 0
  reynld = 1.0D+00
  rjpnew = 0.0D+00
  rjpold = 0.0D+00
  tolnew = 0.0001
  tolsec = 0.0001
  uv_file = 'uv_000.txt'
  xbleft = 1.0D+00
  xbrite = 3.0D+00
  xlngth = 10.0D+00
  xprof = 4.0D+00
  xy_file = 'xy_000.txt'
  ylngth = 3.0D+00

  write ( *, * ) ' '
  write ( *, * ) '  The bump will be generated with a height of ', aprof
  write ( *, * ) ' '
  write ( *, * ) '  NX = ', nx
  write ( *, * ) '  NY = ', ny
  write ( *, * ) '  Number of elements = ', nelemn
  write ( *, * ) '  Reynolds number =  ', reynld
  write ( *, * ) '  Secant tolerance = ', tolsec
  write ( *, * ) '  Newton tolerance = ', tolnew
!
!  SETGRD constructs grid, numbers unknowns, calculates areas,
!  and points for midpoint quadrature rule.
!
  call setgrd ( ibump, indx, insc, isotri, iwrite, long, maxeqn, mx, my, &
    nelemn, neqn, nnodes, node, np, nx, ny, xbleft, xbrite, xlngth )
!
!  Set the X Y coordinates of grid points.
!
  ypert = aprof
  call setxy ( iwrite, long, mx, my, np, nx, ny, xc, xlngth, yc, ylngth, ypert )
!
!  Set the quadrature points.
!
  call setqud ( area, isotri, iwrite, nelemn, nnodes, node, np, nquad, &
    xc, xm, yc, ym )
!
!  Set basis functions values at grid points
!
  call setbas ( isotri, nelemn, nnodes, node, np, nquad, phi, psi, xc, xm, &
    yc, ym )
!
!  Find points on velocity profile sampling line
!
  call setlin ( iline, indx, iwrite, long, mx, my, np, nx, ny, xlngth, xprof )
!
!  Compute the bandwidth.
!
  call setban ( indx, insc, maxrow, nband, nelemn, nlband, nnodes, node, &
    np, nrow )
!
!  Generate the solution of the Navier Stokes equations, using
!  an initial estimate of G = 0.
!
  g(1:neqn) = 0.0D+00

  call nstoke ( a, area, f, g, indx, insc, isotri, maxnew, &
    maxrow, nband, nelemn, neqn, nlband, nnodes, node, np, nquad, &
    nrow, numnew, phi, psi, reynld, tolnew, xc, xm, yc, ym )
!
!  Compute the residual for the solution G.
!
  call resid ( area, g, indx, insc, isotri, iwrite, nelemn, neqn, nnodes, &
    node, np, nquad, phi, psi, res, reynld, xc, xm, yc, ym )
!
!  Copy the flow along the profile line.
!
  call getg ( g, iline, my, neqn, uprof )

  if ( 1 <= iwrite ) then
    write(*,*)' '
    write(*,*)'Velocity profile:'
    write(*,*)' '
    write(*,'(5g14.6)') (uprof(i),i = 1,my)
  end if
!
!  Calculate Gram matrix GR and vector R = line integral of UPROF*phi
!
  call gram ( gr, iline, indx, iwrite, my, nelemn, nnodes, node, &
    np, r, uprof, xc, xprof, yc )
!
!  XY_WRITE creates a TABLE "XY" file
!
  call file_name_inc ( xy_file )

  call get_unit ( xy_unit )

  open ( unit = xy_unit, file = xy_file, status = 'replace' )

  call xy_write ( xy_unit, np, xc, yc )

  close ( unit = xy_unit )

  call file_name_inc ( uv_file )

  call get_unit ( uv_unit )

  open ( unit = uv_unit, file = uv_file, status = 'replace' )

  call uv_write ( f, indx, uv_unit, neqn, np, yc )

  close ( unit = uv_unit )
!
!  Destroy information about true solution before beginning
!  secant iteration.
!
  g(1:neqn) = 0.0D+00
!
!  Secant iteration loop:
!
  do iter = 1, maxsec

    write(*,*) ' '
    write(*,*) 'Secant iteration ',iter

    numsec = numsec + 1
!
!  Update the grid.
!
    ypert = anew
    call setxy ( iwrite, long, mx, my, np, nx, ny, xc, xlngth, yc, &
      ylngth, ypert )
!
!  Set the quadrature points.
!
    call setqud ( area, isotri, iwrite, nelemn, nnodes, node, np, nquad, &
      xc, xm, yc, ym )
!
!  Set basis functions values at grid points.
!
    call setbas(isotri,nelemn,nnodes,node,np,nquad,phi,psi,xc,xm,yc,ym)
!
!  Solve for the flow at the new value of the parameter.
!
    call nstoke(a,area,f,g,indx,insc,isotri,maxnew, &
      maxrow,nband,nelemn,neqn,nlband,nnodes,node,np,nquad, &
      nrow,numnew,phi,psi,reynld,tolnew,xc,xm,yc,ym)
!
!  Get the velocity profile UPROF along sampling line.
!
    call getg ( g, iline, my, neqn, uprof )

    if ( 1 <= iwrite ) then
      write ( *, * ) ' '
      write ( *, * ) 'Velocity profile:'
      write ( *, * ) ' '
      write ( *, '(5g14.6)' ) uprof(1:my)
    end if
!
!  Solve linear system for sensitivities.
!
    itype = -2

    call linsys ( a, area, sens, g, indx, insc, isotri, itype, maxrow, &
      nband, nelemn, neqn, nlband, nnodes, node, np, nquad, nrow, &
      phi, psi, reynld, xc, xm, yc, ym )
!
!  Get the sensitivities DCDA along sampling line.
!
    call getg ( sens, iline, my, neqn, dcda )

    if ( 2 <= iwrite ) then
      write(*,*)' '
      write(*,*)'Sensitivities:'
      write(*,*)' '
      write(*,'(5g14.6)') dcda(1:my)
    end if
!
!  Evaluate J prime at current value of parameter where J is
!  the functional to be minimized.
!
!  JPRIME = 2.0D+00 * DCDA(I) * ( GR(I,J) * UPROF(J) - R(I) )
!
    rjpnew = 0.0D+00
    do i = 1, my
      temp = -r(i)
      do j = 1, my
        temp = temp + gr(i,j) * uprof(j)
      end do
      rjpnew = rjpnew + 2.0D+00 * dcda(i) * temp
    end do
!
!  Write data to graphics files.
!
    call file_name_inc ( xy_file )

    call get_unit ( xy_unit )

    open ( unit = xy_unit, file = xy_file, status = 'replace' )

    call xy_write ( xy_unit, np, xc, yc )

    close ( unit = xy_unit )

    call file_name_inc ( uv_file )

    open ( unit = uv_unit, file = uv_file, status = 'replace' )

    call uv_write ( f, indx, uv_unit, neqn, np, yc )

    close ( unit = uv_unit )

    write ( *, * ) ' '
    write ( *, * ) '  Parameter = ', anew, ' J prime=', rjpnew
!
!  Update the estimate of the parameter using the secant step.
!
    if ( 1 < iter ) then
      anext = aold - rjpold * ( anew - aold ) / ( rjpnew - rjpold )
    end if

    aold = anew
    anew = anext
    rjpold = rjpnew
    if ( anew /= 0.0D+00 ) then
      test = abs ( anew - aold ) / anew
    else
      test = 0.0D+00
    end if

    write ( *, * ) '  New value of parameter = ', anew
    write ( *, * ) '  Convergence test = ', test

    if ( abs ( ( anew - aold ) ) <= abs ( anew ) * tolsec .and. 1 < iter ) then
      write ( *, * ) 'Secant iteration converged.'
      go to 40
    end if

  end do

  write ( *, * ) '  Secant iteration failed to converge.'

 40   continue
!
!  Produce total CPU time used.
!
  call cpu_time ( cpu2 )

  write ( *, '(a)' ) ' '
  write ( *, * ) '  Total execution time = ', cpu2 - cpu1, ' seconds.'
  write ( *, * ) '  Number of secant steps = ', numsec
  write ( *, * ) '  Number of Newton steps = ', numnew
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BUMP:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
function bsp ( it, iq, id, nelemn, nnodes, node, np, xc, xq, yc, yq )

!*****************************************************************************80
!
!! BSP evaluates the linear basis function associated with pressure.
!
!  Discussion:
!
!    The local reference element:
!
!    ^
!    |
!    1  3
!    |  |\
!    |  | \
!    |  |  \
!    |  |   \
!    0  1----2
!    |
!    +--0----1---->
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IT, the element in which the basis function is defined.
!
!    Input, integer ( kind = 4 ) IQ, the index (1, 2 or 3) of the local node at
!    which the basis function is to be evaluated.
!
!    Input, integer ( kind = 4 ) ID, the index (1, 2 or 3) of the local node
!    associated with the basis function.
!
  implicit none

  integer ( kind = 4 ) nelemn
  integer ( kind = 4 ) nnodes
  integer ( kind = 4 ) np

  real ( kind = 8 ) bsp
  real ( kind = 8 ) d
  integer ( kind = 4 ) g1
  integer ( kind = 4 ) g2
  integer ( kind = 4 ) g3
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) id
  integer ( kind = 4 ) iq
  integer ( kind = 4 ) it
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) l2
  integer ( kind = 4 ) l3
  integer ( kind = 4 ) node(nelemn,nnodes)
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) xq
  real ( kind = 8 ) yc(np)
  real ( kind = 8 ) yq
!
!  L1, L2, L3 are the appropriate local node indices.
!
  l1 = iq
  l2 = i4_wrap ( iq + 1, 1, 3 )
  l3 = i4_wrap ( iq + 2, 1, 3 )
!
!  G1, G2, G3 are the global node indices.
!
  g1 = node(it,l1)
  g2 = node(it,l2)
  g3 = node(it,l3)

  d = (xc(g2)-xc(g1))*(yc(g3)-yc(g1))-(xc(g3)-xc(g1))*(yc(g2)-yc(g1))

  if ( id == 1 ) then
    bsp = 1.0D+00 +((yc(g2)-yc(g3))*(xq-xc(g1))+(xc(g3)-xc(g2))*(yq-yc(g1)))/d
  else if ( id == 2 ) then
    bsp = (yc(g2)-yc(g3))/d
  else if ( id == 3 ) then
    bsp = (xc(g3)-xc(g2))/d
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BSP - Fatal error!'
    write ( *, '(a)' ) '  Illegal local index value for linear basis.'
    write ( *, '(a)' ) '  Legal values are 1, 2 or 3.'
    write ( *, '(a,i6)' ) '  The input value was ID = ', id
    stop
  end if

  return
end
subroutine daxpy ( n, da, dx, incx, dy, incy )

!*****************************************************************************80
!
!! DAXPY computes constant times a vector plus a vector.
!
!  Discussion:
!
!    Uses unrolled loops for increments equal to one.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2005
!
!  Author:
!
!    Jack Dongarra
!
!  Reference:
!
!    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, (Society for Industrial and Applied Mathematics),
!    3600 University City Science Center,
!    Philadelphia, PA, 19104-2688.
!    ISBN 0-89871-172-X
!
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    Algorithm 539,
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

  real ( kind = 8 ) da
  real ( kind = 8 ) dx(*)
  real ( kind = 8 ) dy(*)
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

  if ( da  == 0.0D+00 ) then
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

    do i = 1, m
      dy(i) = dy(i) + da * dx(i)
    end do

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
!    This routine uses unrolled loops for increments equal to one.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2005
!
!  Author:
!
!    Jack Dongarra
!
!  Reference:
!
!    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, (Society for Industrial and Applied Mathematics),
!    3600 University City Science Center,
!    Philadelphia, PA, 19104-2688.
!    ISBN 0-89871-172-X
!
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    Algorithm 539,
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

  real ( kind = 8 ) ddot
  real ( kind = 8 ) dtemp
  real ( kind = 8 ) dx(*)
  real ( kind = 8 ) dy(*)
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
subroutine dgbfa ( abd, lda, n, ml, mu, ipvt, info )

!*****************************************************************************80
!
!! DGBFA factors a real band matrix by elimination.
!
!  Discussion:
!
!    DGBFA is usually called by DGBCO, but it can be called
!    directly with a saving in time if RCOND is not needed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2005
!
!  Author:
!
!    FORTRAN90 translation by John Burkardt.
!
!  Reference:
!
!    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, (Society for Industrial and Applied Mathematics),
!    3600 University City Science Center,
!    Philadelphia, PA, 19104-2688.
!    ISBN 0-89871-172-X
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) ABD(LDA,N).  On input, the matrix in band
!    storage.  The columns of the matrix are stored in the columns of ABD
!    and the diagonals of the matrix are stored in rows ML+1 through
!    2*ML+MU+1 of ABD.  On output, an upper triangular matrix in band storage
!    and the multipliers which were used to obtain it.  The factorization
!    can be written A = L*U where L is a product of permutation and unit lower
!    triangular matrices and U is upper triangular.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the array ABD.
!    2*ML + MU + 1 <= LDA is required.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) ML, MU, the number of diagonals below and above the
!    main diagonal.  0 <= ML < N, 0 <= MU < N.
!
!    Output, integer ( kind = 4 ) IPVT(N), the pivot indices.
!
!    Output, integer ( kind = 4 ) INFO, error flag.
!    0, normal value.
!    K, if U(K,K) == 0.0D+00.  This is not an error condition for this
!      subroutine, but it does indicate that DGBSL will divide by zero if
!      called.  Use RCOND in DGBCO for a reliable indication of singularity.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 8 ) abd(lda,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i0
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) idamax
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j0
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) ju
  integer ( kind = 4 ) jz
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lm
  integer ( kind = 4 ) m
  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mm
  integer ( kind = 4 ) mu
  real ( kind = 8 ) t

  m = ml + mu + 1
  info = 0
!
!  Zero initial fill-in columns.
!
  j0 = mu + 2
  j1 = min ( n, m ) - 1

  do jz = j0, j1
    i0 = m + 1 - jz
    do i = i0, ml
      abd(i,jz) = 0.0D+00
    end do
  end do

  jz = j1
  ju = 0
!
!  Gaussian elimination with partial pivoting.
!
  do k = 1, n-1
!
!  Zero out the next fill-in column.
!
    jz = jz + 1
    if ( jz <= n ) then
      abd(1:ml,jz) = 0.0D+00
    end if
!
!  Find L = pivot index.
!
    lm = min ( ml, n-k )
    l = idamax ( lm+1, abd(m,k), 1 ) + m - 1
    ipvt(k) = l + k - m
!
!  Zero pivot implies this column already triangularized.
!
    if ( abd(l,k) == 0.0D+00 ) then

      info = k
!
!  Interchange if necessary.
!
    else

      if ( l /= m ) then
        t = abd(l,k)
        abd(l,k) = abd(m,k)
        abd(m,k) = t
      end if
!
!  Compute multipliers.
!
      t = -1.0D+00 / abd(m,k)
      call dscal ( lm, t, abd(m+1,k), 1 )
!
!  Row elimination with column indexing.
!
      ju = min ( max ( ju, mu+ipvt(k) ), n )
      mm = m

      do j = k+1, ju
        l = l - 1
        mm = mm - 1
        t = abd(l,j)
        if ( l /= mm ) then
          abd(l,j) = abd(mm,j)
          abd(mm,j) = t
        end if
        call daxpy ( lm, t, abd(m+1,k), 1, abd(mm+1,j), 1 )
      end do

    end if

  end do

  ipvt(n) = n

  if ( abd(m,n) == 0.0D+00 ) then
    info = n
  end if

  return
end
subroutine dgbsl ( abd, lda, n, ml, mu, ipvt, b, job )

!*****************************************************************************80
!
!! DGBSL solves a real banded system factored by DGBCO or DGBFA.
!
!  Discussion:
!
!    DGBSL can solve either A * X = B  or  A' * X = B.
!
!    A division by zero will occur if the input factor contains a
!    zero on the diagonal.  Technically this indicates singularity
!    but it is often caused by improper arguments or improper
!    setting of LDA.  It will not occur if the subroutines are
!    called correctly and if DGBCO has set 0.0 < RCOND
!    or DGBFA has set INFO == 0.
!
!    To compute inverse(A) * C  where C is a matrix with P columns:
!
!      call dgbco ( abd, lda, n, ml, mu, ipvt, rcond, z )
!
!      if ( rcond is too small ) then
!        exit
!      end if
!
!      do j = 1, p
!        call dgbsl ( abd, lda, n, ml, mu, ipvt, c(1,j), 0 )
!      end do
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2005
!
!  Author:
!
!    FORTRAN90 translation by John Burkardt.
!
!  Reference:
!
!    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, (Society for Industrial and Applied Mathematics),
!    3600 University City Science Center,
!    Philadelphia, PA, 19104-2688.
!    ISBN 0-89871-172-X
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ABD(LDA,N), the output from DGBCO or DGBFA.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the array ABD.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) ML, MU, the number of diagonals below and above the
!    main diagonal.  0 <= ML < N, 0 <= MU < N.
!
!    Input, integer ( kind = 4 ) IPVT(N), the pivot vector from DGBCO or DGBFA.
!
!    Input/output, real ( kind = 8 ) B(N).  On input, the right hand side.
!    On output, the solution.
!
!    Input, integer ( kind = 4 ) JOB, job choice.
!    0, solve A*X=B.
!    nonzero, solve A'*X=B.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 8 ) abd(lda,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) job
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) la
  integer ( kind = 4 ) lb
  integer ( kind = 4 ) lm
  integer ( kind = 4 ) m
  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
  real ( kind = 8 ) ddot
  real ( kind = 8 ) t

  m = mu + ml + 1
!
!  JOB = 0, Solve A * x = b.
!
!  First solve L * y = b.
!
  if ( job == 0 ) then

    if ( 0 < ml ) then

      do k = 1, n-1
        lm = min ( ml, n-k )
        l = ipvt(k)
        t = b(l)
        if ( l /= k ) then
          b(l) = b(k)
          b(k) = t
        end if
      call daxpy ( lm, t, abd(m+1,k), 1, b(k+1), 1 )
      end do

    end if
!
!  Now solve U * x = y.
!
    do k = n, 1, -1
      b(k) = b(k) / abd(m,k)
      lm = min ( k, m ) - 1
      la = m - lm
      lb = k - lm
      t = -b(k)
      call daxpy ( lm, t, abd(la,k), 1, b(lb), 1 )
    end do
!
!  JOB nonzero, solve A' * x = b.
!
!  First solve U' * y = b.
!
  else

    do k = 1, n
      lm = min ( k, m ) - 1
      la = m - lm
      lb = k - lm
      t = ddot ( lm, abd(la,k), 1, b(lb), 1 )
      b(k) = ( b(k) - t ) / abd(m,k)
    end do
!
!  Now solve L' * x = y.
!
    if ( 0 < ml ) then

      do k = n-1, 1, -1
        lm = min ( ml, n-k )
        b(k) = b(k) + ddot ( lm, abd(m+1,k), 1, b(k+1), 1 )
        l = ipvt(k)
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
subroutine dscal ( n, sa, x, incx )

!*****************************************************************************80
!
!! DSCAL scales a vector by a constant.
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
!    Jack Dongarra
!
!  Reference:
!
!    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, (Society for Industrial and Applied Mathematics),
!    3600 University City Science Center,
!    Philadelphia, PA, 19104-2688.
!    ISBN 0-89871-172-X
!
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
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
!    Input/output, real ( kind = 8 ) X(*), the vector to be scaled.
!
!    Input, integer ( kind = 4 ) INCX, the increment between successive entries of X.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 8 ) sa
  real ( kind = 8 ) x(*)

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
subroutine file_name_inc ( file_name )

!*****************************************************************************80
!
!! FILE_NAME_INC increments a partially numeric filename.
!
!  Discussion:
!
!    It is assumed that the digits in the name, whether scattered or
!    connected, represent a number that is to be increased by 1 on
!    each call.  If this number is all 9's on input, the output number
!    is all 0's.  Non-numeric letters of the name are unaffected.
!
!    If the name is empty, then the routine stops.
!
!    If the name contains no digits, the empty string is returned.
!
!  Example:
!
!      Input            Output
!      -----            ------
!      'a7to11.txt'     'a7to12.txt'
!      'a7to99.txt'     'a8to00.txt'
!      'a9to99.txt'     'a0to00.txt'
!      'cat.txt'        ' '
!      ' '              STOP!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) FILE_NAME.
!    On input, a character string to be incremented.
!    On output, the incremented string.
!
  implicit none

  character c
  integer ( kind = 4 ) change
  integer ( kind = 4 ) digit
  character ( len = * ) file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) lens

  lens = len_trim ( file_name )

  if ( lens <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_NAME_INC - Fatal error!'
    write ( *, '(a)' ) '  The input string is empty.'
    stop
  end if

  change = 0

  do i = lens, 1, -1

    c = file_name(i:i)

    if ( lge ( c, '0' ) .and. lle ( c, '9' ) ) then

      change = change + 1

      digit = ichar ( c ) - 48
      digit = digit + 1

      if ( digit == 10 ) then
        digit = 0
      end if

      c = char ( digit + 48 )

      file_name(i:i) = c

      if ( c /= '0' ) then
        return
      end if

    end if

  end do

  if ( change == 0 ) then
    file_name = ' '
    return
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
!    A "free" FORTRAN unit number is an integer ( kind = 4 ) between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is an integer ( kind = 4 ) between 1 and 99, representing a
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
subroutine getg ( f, iline, my, neqn, u )

!*****************************************************************************80
!
!! GETG extracts the values of a quantity along the profile line.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) F(NEQN), the finite element coefficients.
!
!    Input, integer ( kind = 4 ) ILINE(MY), the index of the finite element 
!    coefficient associated with each node on the profile.
!
!    Input, integer ( kind = 4 ) MY, the number of nodes along the profile line.
!
!    Input, integer ( kind = 4 ) NEQN, the number of finite element coefficients.
!
!    Output, real ( kind = 8 ) U(MY), the finite element coefficients
!    associated with each node of the profile.
!
  implicit none

  integer ( kind = 4 ) my
  integer ( kind = 4 ) neqn

  real ( kind = 8 ) f(neqn)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iline(my)
  integer ( kind = 4 ) j
  real ( kind = 8 ) u(my)

  do i = 1, my
    u(i) = 0.0D+00
    j = iline(i)
    if ( j <= 0 ) then
      u(i) = 0.0D+00
    else
      u(i) = f(j)
    end if
  end do

  return
end
subroutine gram ( gr, iline, indx, iwrite, my, nelemn, nnodes, node, &
  np, r, uprof, xc, xprof, yc )

!*****************************************************************************80
!
!! GRAM computes and stores the Gram matrix.
!
!  Discussion:
!
!    The routine computes the Gram matrix and the vector
!    whose components are the line integral of UI*phi(j).
!
!    A three-point Gauss quadrature rule is used to evaluate the line integral.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) GR(MY,MY), the Gram matrix.
!
!    Input, integer ( kind = 4 ) ILINE(MY), the list of global unknown numbers
!    along the profile line.
!
!    Input, integer ( kind = 4 ) INDX(NP,2), contains, for each node I, the
!    index of the U and V velocities at that node, or 0.
!
!    Input, integer ( kind = 4 ) IWRITE, controls the amount of output printed.
!
!    Input, integer ( kind = 4 ) MY, the number of nodes along the profile line.
!
!    Input, integer ( kind = 4 ) NELEMN, the number of elements.
!
!    Input, integer ( kind = 4 ) NNODES, the number of nodes per element, which
!    is 6 for the quadratic triangular elements in use here.
!
!    Input, integer ( kind = 4 ) NODE(NELEMN,NNODES), contains, for each element,
!    the global node indices of the element nodes.
!
!    Input, integer ( kind = 4 ) NP, the number of nodes.
!
!    Output, real ( kind = 8 ) R(MY), the line integral of UPROF * PHI.
!
!    Input, real ( kind = 8 ) UPROF(MY), the horizontal velocity
!    along the profile line.
!
!    Input, real ( kind = 8 ) XC(NP), the X coordinates of the nodes.
!
!    Input, real ( kind = 8 ) XPROF, the X coordinate at which the
!    profile is measured.  This value should be a grid value!
!
!    Input, real ( kind = 8 ) YC(NP), the Y coordinates of the nodes.
!
  implicit none

  integer ( kind = 4 ) my
  integer ( kind = 4 ) nelemn
  integer ( kind = 4 ) nnodes
  integer ( kind = 4 ) np

  real ( kind = 8 ) ar
  real ( kind = 8 ) bb
  real ( kind = 8 ) bbb
  real ( kind = 8 ) bbx
  real ( kind = 8 ) bby
  real ( kind = 8 ) bma2
  real ( kind = 8 ) bx
  real ( kind = 8 ) by
  real ( kind = 8 ) gr(my,my)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) igetl
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) iline(my)
  integer ( kind = 4 ) indx(np,2)
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) ipp
  integer ( kind = 4 ) iq
  integer ( kind = 4 ) iqq
  integer ( kind = 4 ) iquad
  integer ( kind = 4 ) it
  integer ( kind = 4 ) iun
  integer ( kind = 4 ) iwrite
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kk
  integer ( kind = 4 ) node(nelemn,nnodes)
  real ( kind = 8 ) r(my)
  real ( kind = 8 ) ubc
  real ( kind = 8 ) ubdry
  real ( kind = 8 ) uiqdpt
  real ( kind = 8 ) uprof(my)
  real ( kind = 8 ) wt(3)
  real ( kind = 8 ) x
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) xprof
  real ( kind = 8 ) y
  real ( kind = 8 ) yc(np)
  real ( kind = 8 ) yq(3)
!
!  Values for 3 point Gauss quadrature.
!
  wt(1) = 5.0D+00 / 9.0D+00
  wt(2) = 8.0D+00 / 9.0D+00
  wt(3) = 5.0D+00 / 9.0D+00

  yq(1) = -0.7745966692D+00
  yq(2) =  0.0D+00
  yq(3) =  0.7745966692D+00

  r(1:my) = 0.0D+00
  gr(1:my,1:my) = 0.0D+00
!
!  Compute the line integral by looping over intervals along line
!  using three point Gauss quadrature
!
  do it = 1, nelemn
!
!  Check to see if we are in a triangle with a side along line
!  x = xprof
!
      k = node(it,1)
      kk = node(it,2)

      if ( 1.0D-04 < abs ( xc(k) - xprof ) .or. &
         1.0D-04 < abs ( xc(kk) - xprof ) ) then
        go to 70
      end if

      do iquad = 1, 3

        bma2 = ( yc(kk) - yc(k) ) / 2.0D+00
        ar = bma2 * wt(iquad)
        x = xprof
        y = yc(k) + bma2 * ( yq(iquad) + 1.0D+00 )
!
!  Compute U internal at quadrature points.
!
        uiqdpt = 0.0D+00

        do iq = 1, nnodes
          if ( iq == 1 .or. iq == 2 .or. iq == 4 ) then
            call qbf ( x, y, it, iq, bb, bx, by, nelemn, nnodes, &
              node, np, xc, yc )
            ip = node(it,iq)
            iun = indx(ip,1)
            if ( 0 < iun ) then
              ii = igetl(iun,iline,my)
              uiqdpt = uiqdpt + bb * uprof(ii)
            else if ( iun == -1) then
              ubc = ubdry(1,yc(ip))
              uiqdpt = uiqdpt + bb * ubc
            end if
          end if
        end do
!
!  Only loop over nodes lying on line x = xprof
!
        do iq = 1, nnodes
          if ( iq == 1 .or. iq == 2 .or. iq == 4 ) then
            ip = node(it,iq)
            call qbf(x,y,it,iq,bb,bx,by,nelemn,nnodes,node,np,xc,yc)
            i = indx(ip,1)
            if ( 0 < i ) then
              ii = igetl(i,iline,my)
              r(ii) = r(ii) + bb * uiqdpt * ar
              do iqq = 1, nnodes
                if ( iqq == 1 .or. iqq == 2 .or. iqq == 4 ) then
                  ipp = node(it,iqq)
                  call qbf(x,y,it,iqq,bbb,bbx,bby,nelemn,nnodes,node,np,xc,yc)
                  j = indx(ipp,1)
                  if ( j /= 0 ) then
                    jj = igetl(j,iline,my)
                    gr(ii,jj) = gr(ii,jj)+bb*bbb*ar
                  end if
                end if
              end do
            end if
          end if
        end do

      end do

 70     continue

  end do

  if ( 3 <= iwrite ) then
    write(*,*)' '
    write(*,*)'Gram matrix:'
    write(*,*)' '
    do i = 1, my
      do j = 1, my
        write(*,*)i,j,gr(i,j)
      end do
    end do
    write(*,*)' '
    write(*,*)'R vector:'
    write(*,*)' '
    do i = 1, my
      write(*,*)r(i)
    end do

  end if

  return
end
function i4_modp ( i, j )

!*****************************************************************************80
!
!! I4_MODP returns the nonnegative remainder of integer ( kind = 4 ) division.
!
!  Discussion:
!
!    If
!      NREM = I4_MODP ( I, J )
!      NMULT = ( I - NREM ) / J
!    then
!      I = J * NMULT + NREM
!    where NREM is always nonnegative.
!
!    The MOD function computes a result with the same sign as the
!    quantity being divided.  Thus, suppose you had an angle A,
!    and you wanted to ensure that it was between 0 and 360.
!    Then mod(A,360) would do, if A was positive, but if A
!    was negative, your result would be between -360 and 0.
!
!    On the other hand, I4_MODP(A,360) is between 0 and 360, always.
!
!  Example:
!
!        I     J     MOD I4_MODP    Factorization
!
!      107    50       7       7    107 =  2 *  50 + 7
!      107   -50       7       7    107 = -2 * -50 + 7
!     -107    50      -7      43   -107 = -3 *  50 + 43
!     -107   -50      -7      43   -107 =  3 * -50 + 43
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
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
!    Input, integer ( kind = 4 ) I, the number to be divided.
!
!    Input, integer ( kind = 4 ) J, the number that divides I.
!
!    Output, integer ( kind = 4 ) I4_MODP, the nonnegative remainder when I is
!    divided by J.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ) j
  integer ( kind = 4 ) value

  if ( j == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_MODP - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal divisor J = ', j
    stop
  end if

  value = mod ( i, j )

  if ( value < 0 ) then
    value = value + abs ( j )
  end if

  i4_modp = value

  return
end
function i4_wrap ( ival, ilo, ihi )

!*****************************************************************************80
!
!! I4_WRAP forces an I4 to lie between given limits by wrapping.
!
!  Example:
!
!    ILO = 4, IHI = 8
!
!    I  Value
!
!    -2     8
!    -1     4
!     0     5
!     1     6
!     2     7
!     3     8
!     4     4
!     5     5
!     6     6
!     7     7
!     8     8
!     9     4
!    10     5
!    11     6
!    12     7
!    13     8
!    14     4
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 August 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IVAL, an integer ( kind = 4 ) value.
!
!    Input, integer ( kind = 4 ) ILO, IHI, the desired bounds for the value.
!
!    Output, integer ( kind = 4 ) I4_WRAP, a "wrapped" version of IVAL.
!
  implicit none

  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) value
  integer ( kind = 4 ) wide

  jlo = min ( ilo, ihi )
  jhi = max ( ilo, ihi )

  wide = jhi - jlo + 1

  if ( wide == 1 ) then
    value = jlo
  else
    value = jlo + i4_modp ( ival - jlo, wide )
  end if

  i4_wrap = value

  return
end
function idamax ( n, dx, incx )

!*****************************************************************************80
!
!! IDAMAX finds the index of the vector element of maximum absolute value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 April 1999
!
!  Reference:
!
!    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, (Society for Industrial and Applied Mathematics),
!    3600 University City Science Center,
!    Philadelphia, PA, 19104-2688.
!    ISBN 0-89871-172-X
!
!    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    Algorithm 539,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, real ( kind = 8 ) X(*), the vector to be examined.
!
!    Input, integer ( kind = 4 ) INCX, the increment between successive entries of SX.
!
!    Output, integer ( kind = 4 ) IDAMAX, the index of the element of SX of maximum
!    absolute value.
!
  implicit none

  real ( kind = 8 ) dmax
  real ( kind = 8 ) dx(*)
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
function igetl ( i, iline, my )

!*****************************************************************************80
!
!! IGETL gets the local unknown number along the profile line.
!
!  Discussion:
!
!    For our problem, the profile line is specified by X = XPROF.
!    We are given a global unknown number, and need to determine the
!    local unknown number.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the global unknown number.
!
!    Input, integer ( kind = 4 ) ILINE(MY), the list of global unknown numbers
!    along the profile line.
!
!    Input, integer ( kind = 4 ) MY, the number of nodes on the profile line.
!
!    Output, integer ( kind = 4 ) IGETL, the index of the node in the profile line
!    at which global unknown I occurs, or -1 if there is no such entry.
!
  implicit none

  integer ( kind = 4 ) my

  integer ( kind = 4 ) i
  integer ( kind = 4 ) igetl
  integer ( kind = 4 ) iline(my)
  integer ( kind = 4 ) j

  igetl = -1

  do j = 1, my
    if ( iline(j) == i ) then
      igetl = j
      return
    end if
  end do

  return
end
subroutine linsys ( a, area, f, g, indx, insc, isotri, itype, maxrow, &
  nband, nelemn, neqn, nlband, nnodes, node, np, nquad, nrow, &
  phi, psi, reynld, xc, xm, yc, ym )

!*****************************************************************************80
!
!! LINSYS solves the linearized Navier Stokes equation.
!
!  Discussion:
!
!    ITYPE = -1 for solutions of the Navier Stokes equation.
!    ITYPE = -2 for solutions of the sensitivity equation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) AREA(NELEMN), the area of the elements.
!
!    Input, integer ( kind = 4 ) INDX(MAXNP,2), contains, for each node I, the
!    index of the U and V velocities at that node, or 0.
!
!    Input, integer ( kind = 4 ) NNODES, the number of nodes per element, which
!    is 6 for the quadratic triangular elements in use here.
!
!  Local Parameters:
!
!    Local, integer ( kind = 4 ) IPIVOT(NEQN), pivoting space.
!
  implicit none

  integer ( kind = 4 ) maxrow
  integer ( kind = 4 ) nelemn
  integer ( kind = 4 ) neqn
  integer ( kind = 4 ) nnodes
  integer ( kind = 4 ) np
  integer ( kind = 4 ) nquad

  real ( kind = 8 ) a(maxrow,neqn)
  real ( kind = 8 ) aij
  real ( kind = 8 ) ar
  real ( kind = 8 ) area(nelemn)
  real ( kind = 8 ) bb
  real ( kind = 8 ) bbb
  real ( kind = 8 ) bbbl
  real ( kind = 8 ) bbl
  real ( kind = 8 ) bbx
  real ( kind = 8 ) bby
  real ( kind = 8 ) bx
  real ( kind = 8 ) by
  real ( kind = 8 ) det
  real ( kind = 8 ) etax
  real ( kind = 8 ) etay
  real ( kind = 8 ) f(neqn)
  real ( kind = 8 ) g(neqn)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihor
  integer ( kind = 4 ) indx(np,2)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) insc(np)
  integer ( kind = 4 ) ioff
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) ipivot(neqn)
  integer ( kind = 4 ) ipp
  integer ( kind = 4 ) iprs
  integer ( kind = 4 ) iq
  integer ( kind = 4 ) iqq
  integer ( kind = 4 ) iquad
  integer ( kind = 4 ) isotri(nelemn)
  integer ( kind = 4 ) it
  integer ( kind = 4 ) itype
  integer ( kind = 4 ) iuse
  integer ( kind = 4 ) iver
  integer ( kind = 4 ) j
  integer ( kind = 4 ) job
  integer ( kind = 4 ) jp
  integer ( kind = 4 ) ju
  integer ( kind = 4 ) juse
  integer ( kind = 4 ) jv
  integer ( kind = 4 ) nband
  integer ( kind = 4 ) nlband
  integer ( kind = 4 ) node(nelemn,nnodes)
  integer ( kind = 4 ) nrow
  real ( kind = 8 ) phi(nelemn,nquad,nnodes,3)
  real ( kind = 8 ) psi(nelemn,nquad,nnodes)
  real ( kind = 8 ) reynld
  real ( kind = 8 ) ubc
  real ( kind = 8 ) ubdry
  real ( kind = 8 ) ubump
  real ( kind = 8 ) un(2)
  real ( kind = 8 ) unx(2)
  real ( kind = 8 ) uny(2)
  real ( kind = 8 ) visc
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) xix
  real ( kind = 8 ) xiy
  real ( kind = 8 ) xm(nelemn,nquad)
  real ( kind = 8 ) xq
  real ( kind = 8 ) yc(np)
  real ( kind = 8 ) ym(nelemn,nquad)
  real ( kind = 8 ) yq

  ioff = nlband + nlband + 1

  visc = 1.0D+00 / reynld

  f(1:neqn) = 0.0D+00
  a(1:nrow,1:neqn) = 0.0D+00

  do it = 1, nelemn

    ar = area(it) / 3.0D+00

    do iquad = 1, nquad

      yq = ym(it,iquad)
      xq = xm(it,iquad)

      if ( isotri(it) == 1 ) then
        call trans(det,etax,etay,it,nelemn,nnodes,node,np,xc,xix,xiy,xq,yc,yq)
        ar = det * area(it) / 3.0D+00
      end if

      call uval ( etax, etay, g, indx, isotri, it, nelemn, neqn, &
        nnodes, node, np, un, uny, unx, xc, xix, xiy, xq, yc, yq )
!
!  For each basis function:
!
      do iq = 1, nnodes

        ip = node(it,iq)
        bb = phi(it,iquad,iq,1)
        bx = phi(it,iquad,iq,2)
        by = phi(it,iquad,iq,3)
        bbl = psi(it,iquad,iq)
        ihor = indx(ip,1)
        iver = indx(ip,2)
        iprs = insc(ip)

        if ( 0 < ihor ) then
          f(ihor) = f(ihor) + ar * bb*(un(1)*unx(1)+un(2)*uny(1))
        end if

        if ( 0 < iver ) then
          f(iver) = f(iver) + ar * bb*(un(1)*unx(2)+un(2)*uny(2))
        end if
!
!  For another basis function,
!
        do iqq = 1, nnodes

          ipp = node(it,iqq)
          bbb = phi(it,iquad,iqq,1)
          bbx = phi(it,iquad,iqq,2)
          bby = phi(it,iquad,iqq,3)
          bbbl = psi(it,iquad,iqq)
          ju = indx(ipp,1)
          jv = indx(ipp,2)
          jp = insc(ipp)
!
!  Horizontal velocity variable
!
          if ( 0 < ju ) then

            if ( 0 < ihor ) then
              iuse = ihor-ju+ioff
              a(iuse,ju) = a(iuse,ju)+ar*(visc*(by*bby+bx*bbx) &
                +bb*(bbb*unx(1)+bbx*un(1)+bby*un(2)))
            end if

            if ( 0 < iver ) then
              iuse = iver-ju+ioff
              a(iuse,ju) = a(iuse,ju)+ar*bb*bbb*unx(2)
            end if

            if ( 0 < iprs ) then
              iuse = iprs-ju+ioff
              a(iuse,ju) = a(iuse,ju)+ar*bbx*bbl
            end if

          else if ( ju == itype ) then

            if ( ju == -1 ) then
              ubc = ubdry(1,yc(ipp))
            else if ( ju == -2 ) then
              ubc = ubump(g,indx,ipp,iqq,isotri,it,1,nelemn, &
                neqn,nnodes,node,np,xc,yc)
            end if

            if ( 0 < ihor ) then
              aij = ar*(visc*(by*bby+bx*bbx) &
                +bb*(bbb*unx(1)+bbx*un(1)+bby*un(2)))
              f(ihor) = f(ihor)-ubc*aij
            end if

            if ( 0 < iver ) then
              aij = ar*bb*bbb*unx(2)
              f(iver) = f(iver)-ubc*aij
            end if

            if ( 0 < iprs ) then
              aij = ar*bbx*bbl
              f(iprs) = f(iprs)-ubc*aij
            end if

          end if
!
!  Vertical velocity variable
!
          if ( 0 < jv ) then

            if ( 0 < ihor ) then
              iuse = ihor-jv+ioff
              a(iuse,jv) = a(iuse,jv)+ar*bb*bbb*uny(1)
            end if

            if ( 0 < iver ) then
              iuse = iver-jv+ioff
              a(iuse,jv) = a(iuse,jv)+ar*(visc*(by*bby+bx*bbx) &
                +bb*(bbb*uny(2)+bby*un(2)+bbx*un(1)))
            end if

            if ( 0 < iprs ) then
              iuse = iprs-jv+ioff
              a(iuse,jv) = a(iuse,jv)+ar*bby*bbl
            end if

          else if ( jv == itype ) then

            if ( jv == -1 ) then
              ubc = ubdry(2,yc(ipp))
            else if ( jv == -2 ) then
              ubc = ubump(g,indx,ipp,iqq,isotri,it,2,nelemn, &
                neqn,nnodes,node,np,xc,yc)
            end if

            if ( 0 < ihor ) then
              aij = ar*bb*bbb*uny(1)
              f(ihor) = f(ihor)-ubc*aij
            end if

            if ( 0 < iver ) then
              aij = ar*(visc*(by*bby+bx*bbx) &
                +bb*(bbb*uny(2)+bby*un(2)+bbx*un(1)))
              f(iver) = f(iver)-ubc*aij
            end if

            if ( 0 < iprs ) then
              aij = ar * bby * bbl
              f(iprs) = f(iprs) - ubc * aij
            end if

          end if
!
!  Pressure variable
!
          if ( 0 < jp ) then

            if ( 0 < ihor ) then
              iuse = ihor-jp+ioff
              a(iuse,jp) = a(iuse,jp) - ar * bx * bbbl
            end if

            if ( 0 < iver ) then
              iuse = iver-jp+ioff
              a(iuse,jp) = a(iuse,jp) - ar * by * bbbl
            end if

          end if

        end do
      end do
    end do
  end do
!
!  The last equation is "reset" to require that the last pressure
!  be zero.
!
  f(neqn) = 0.0D+00
  do j = neqn-nlband, neqn-1
    i = neqn-j+ioff
    a(i,j) = 0.0D+00
  end do
  a(nband,neqn) = 1.0D+00
!
!  Factor the matrix
!
  call dgbfa ( a, maxrow, neqn, nlband, nlband, ipivot, info )

  if (info /= 0) then
    write (*,*) ' '
    write (*,*) 'LINSYS - fatal error!'
    write (*,*) 'DGBFA returns INFO = ',info
    stop
  end if
!
!  Solve the linear system
!
  job = 0
  call dgbsl ( a, maxrow, neqn, nlband, nlband, ipivot, f, job )

  return
end
subroutine nstoke ( a, area, f, g, indx, insc, isotri, maxnew, &
  maxrow, nband, nelemn, neqn, nlband, nnodes, node, np, nquad, &
  nrow, numnew, phi, psi, reynld, tolnew, xc, xm, yc, ym )

!*****************************************************************************80
!
!! NSTOKE solves the Navier Stokes equation using Taylor-Hood elements.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) AREA(NELEMN), the area of the elements.
!
!    Input, integer ( kind = 4 ) INDX(MAXNP,2), contains, for each node I, the
!    index of the U and V velocities at that node, or 0.
!
!    Input, integer ( kind = 4 ) NNODES, the number of nodes per element, which
!    is 6 for the quadratic triangular elements in use here.
!
  implicit none

  integer ( kind = 4 ) maxrow
  integer ( kind = 4 ) nelemn
  integer ( kind = 4 ) neqn
  integer ( kind = 4 ) nnodes
  integer ( kind = 4 ) np
  integer ( kind = 4 ) nquad

  real ( kind = 8 ) a(maxrow,neqn)
  real ( kind = 8 ) area(nelemn)
  real ( kind = 8 ) diff
  real ( kind = 8 ) f(neqn)
  real ( kind = 8 ) g(neqn)
  integer ( kind = 4 ) idamax
  integer ( kind = 4 ) indx(np,2)
  integer ( kind = 4 ) insc(np)
  integer ( kind = 4 ) isotri(nelemn)
  integer ( kind = 4 ) iter
  integer ( kind = 4 ) itype
  integer ( kind = 4 ) maxnew
  integer ( kind = 4 ) nband
  integer ( kind = 4 ) nlband
  integer ( kind = 4 ) node(nelemn,nnodes)
  integer ( kind = 4 ) nrow
  integer ( kind = 4 ) numnew
  real ( kind = 8 ) phi(nelemn,nquad,nnodes,3)
  real ( kind = 8 ) psi(nelemn,nquad,nnodes)
  real ( kind = 8 ) reynld
  real ( kind = 8 ) tolnew
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) xm(nelemn,nquad)
  real ( kind = 8 ) yc(np)
  real ( kind = 8 ) ym(nelemn,nquad)
!
!  G contains an initial estimate of the solution.
!
  do iter = 1, maxnew

    numnew = numnew + 1

    itype = -1
    call linsys ( a, area, f, g, indx, insc, isotri, itype, maxrow, &
      nband, nelemn, neqn, nlband, nnodes, node, np, nquad, nrow, &
      phi, psi, reynld, xc, xm, yc, ym )
!
!  Check for convergence.
!
    g(1:neqn) = g(1:neqn) - f(1:neqn)

    diff = abs ( g(idamax(neqn,g,1)) )
    write(*,*) 'NSTOKE: Iteration ', iter, ' MaxNorm(diff) = ', diff

    g(1:neqn) = f(1:neqn)

    if ( diff <= tolnew ) then
      write(*,*) 'NSTOKE converged.'
      exit
    end if

    if ( iter == maxnew ) then
      write(*,*) 'NSTOKE failed!'
      stop
    end if

  end do

  return
end
subroutine pval ( g, insc, long, mx, my, nelemn, neqn, nnodes, node, &
  np, press )

!*****************************************************************************80
!
!! PVAL computes a table of pressures at all the nodes.
!
!  Discussion:
!
!    This table is needed in order to prepare output data for PLOT3D.
!
!    This routine does not handle general isoparametric elements
!    accurately.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) MX, the number of nodes in the X direction.
!
!    Input, integer ( kind = 4 ) MY, the number of nodes in the Y direction.
!
!    Input, integer ( kind = 4 ) NNODES, the number of nodes per element, which
!    is 6 for the quadratic triangular elements in use here.
!
  implicit none

  integer ( kind = 4 ) mx
  integer ( kind = 4 ) my
  integer ( kind = 4 ) nelemn
  integer ( kind = 4 ) neqn
  integer ( kind = 4 ) nnodes
  integer ( kind = 4 ) np

  real ( kind = 8 ) g(neqn)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) insc(np)
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) iq
  integer ( kind = 4 ) it
  integer ( kind = 4 ) ivar
  integer ( kind = 4 ) j
  logical long
  integer ( kind = 4 ) node(nelemn,nnodes)
  real ( kind = 8 ) press(mx,my)

  press(1:mx,1:my) = 0.0D+00
!
!  Read the pressures where they are computed.
!  These are "(odd, odd)" points.
!
  do it = 1, nelemn
    do iq = 1, 3
      ip = node(it,iq)
      ivar = insc(ip)
      if ( long ) then
        i = ((ip-1)/my)+1
        j = mod(ip-1,my)+1
      else
        i = mod(ip-1,mx)+1
        j = ((ip-1)/mx)+1
      end if

      if ( 0 < ivar ) then
        press(i,j) = g(ivar)
      end if

    end do
  end do
!
!  Interpolate the pressures at points (even, odd) and (odd, even).
!
  do i = 2, mx-1, 2
    do j = 1, my, 2
      press(i,j) = 0.5D+00 * (press(i-1,j)+press(i+1,j))
    end do
  end do

  do j = 2, my-1, 2
    do i = 1, mx, 2
      press(i,j) = 0.5D+00 * (press(i,j-1)+press(i,j+1))
    end do
  end do
!
!  Interpolate the pressures at points (even,even).
!
  do j = 2, my-1, 2
    do i = 2, mx-1, 2
      press(i,j) = 0.5D+00 * (press(i-1,j-1)+press(i+1,j+1))
    end do
  end do

  return
end
subroutine qbf ( xq, yq, it, in, bb, bx, by, nelemn, nnodes, node, np, xc, yc )

!*****************************************************************************80
!
!! QBF evaluates the quadratic basis functions.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NNODES, the number of nodes per element, which
!    is 6 for the quadratic triangular elements in use here.
!
  implicit none

  integer ( kind = 4 ) nelemn
  integer ( kind = 4 ) nnodes
  integer ( kind = 4 ) np

  real ( kind = 8 ) bb
  real ( kind = 8 ) bx
  real ( kind = 8 ) by
  real ( kind = 8 ) c
  real ( kind = 8 ) d
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) in
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) in3
  integer ( kind = 4 ) inn
  integer ( kind = 4 ) it
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j3
  integer ( kind = 4 ) node(nelemn,nnodes)
  real ( kind = 8 ) s
  real ( kind = 8 ) t
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) xq
  real ( kind = 8 ) yc(np)
  real ( kind = 8 ) yq

  if ( in <= 3 ) then
    in1 = in
    in2 = mod(in,3)+1
    in3 = mod(in+1,3)+1
    i1 = node(it,in1)
    i2 = node(it,in2)
    i3 = node(it,in3)
    d = (xc(i2)-xc(i1))*(yc(i3)-yc(i1))-(xc(i3)-xc(i1))*(yc(i2)-yc(i1))
    t = 1.0D+00 +((yc(i2)-yc(i3))*(xq-xc(i1))+(xc(i3)-xc(i2))*(yq-yc(i1)))/d
    bb = t*(2.0D+00*t-1.0D+00)
    bx = (yc(i2)-yc(i3))*(4.0D+00*t-1.0D+00)/d
    by = (xc(i3)-xc(i2))*(4.0D+00*t-1.0D+00)/d
  else
    inn = in-3
    in1 = inn
    in2 = mod(inn,3)+1
    in3 = mod(inn+1,3)+1
    i1 = node(it,in1)
    i2 = node(it,in2)
    i3 = node(it,in3)
    j1 = i2
    j2 = i3
    j3 = i1
    d = (xc(i2)-xc(i1))*(yc(i3)-yc(i1))-(xc(i3)-xc(i1))*(yc(i2)-yc(i1))
    c = (xc(j2)-xc(j1))*(yc(j3)-yc(j1))-(xc(j3)-xc(j1))*(yc(j2)-yc(j1))
    t = 1.0D+00+((yc(i2)-yc(i3))*(xq-xc(i1))+(xc(i3)-xc(i2))*(yq-yc(i1)))/d
    s = 1.0D+00+((yc(j2)-yc(j3))*(xq-xc(j1))+(xc(j3)-xc(j2))*(yq-yc(j1)))/c
    bb = 4.0D+00 * s*t
    bx = 4.0D+00 * (t*(yc(j2)-yc(j3))/c+s*(yc(i2)-yc(i3))/d)
    by = 4.0D+00 * (t*(xc(j3)-xc(j2))/c+s*(xc(i3)-xc(i2))/d)
  end if

  return
end
function refbsp ( xq, yq, iq )

!*****************************************************************************80
!
!! REFBSP evaluates the linear basis functions in a reference triangle.
!
!  Discussion:
!
!    The reference triangle given here is not the best.  The nodes
!    are not given in the usual positions, and are listed in
!    clockwise order instead of counterclockwise order!
!
!  Diagram:
!
!        2
!       /|
!      / |
!     /  |
!    1---3
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XQ, YQ, the coordinates of a point in
!    the reference triangle.
!
!    Input, integer ( kind = 4 ) IQ, the index of a basis function in the reference triangle.
!
!    Output, real ( kind = 8 ) REFBSP, the value of the IQ-th basis function
!    at the point (XQ,YQ) in the reference triangle.
!
  implicit none

  integer ( kind = 4 ) iq
  real ( kind = 8 ) refbsp
  real ( kind = 8 ) xq
  real ( kind = 8 ) yq

  if ( iq == 1 ) then
    refbsp = 1.0D+00 - xq
  else if ( iq == 2 ) then
    refbsp = yq
  else if ( iq == 3 ) then
    refbsp = xq - yq
  end if

  return
end
subroutine refqbf ( x, y, in, bb, bx, by, etax, etay, xix, xiy )

!*****************************************************************************80
!
!! REFQBF evaluates the quadratic basis functions on reference triangle
!
!  Diagram:
!
!    3
!    |\
!    6 5
!    |  \
!    1-4-2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision X, Y, the coordinates of a point in
!    the reference triangle.
!
!    Input, integer ( kind = 4 ) IN, the index of a basis function in the reference triangle.
!
!    Output, double precision BB, BX, BY, the value of the basis function
!    and its X and Y derivatives at the point (X,Y) in the reference triangle.
!
!    Input, real ( kind = 8 ) ETAX, ETAY, XIX, XIY, the values of
!    dETA/dX, dETA/dY, dXI/dX and dXI/dY at (X,Y).
!
  implicit none

  real ( kind = 8 ) bb
  real ( kind = 8 ) bx
  real ( kind = 8 ) by
  real ( kind = 8 ) etax
  real ( kind = 8 ) etay
  integer ( kind = 4 ) in
  real ( kind = 8 ) tbx
  real ( kind = 8 ) tby
  real ( kind = 8 ) x
  real ( kind = 8 ) xix
  real ( kind = 8 ) xiy
  real ( kind = 8 ) y

  if ( in == 1) then
    bb = 1.0D+00-3.0D+00*x+2.0D+00*x*x
    tbx = -3.0D+00+4.0D+00*x
    tby = 0.0D+00
  else if ( in == 2) then
    bb = -y+2.0D+00*y*y
    tbx = 0.0D+00
    tby = -1.0D+00+4.0D+00*y
  else if (in == 3) then
    bb = -x+2.0D+00*x*x+y-4.0D+00*x*y+2*y*y
    tbx = -1.0D+00+4.0D+00*x-4.0D+00*y
    tby = 1.0D+00-4.0D+00*x+4.0D+00*y
  else if ( in == 4) then
    bb = 4.0D+00*y-4.0D+00*x*y
    tbx = -4.0D+00*y
    tby = 4.0D+00-4.0D+00*x
  else if ( in == 5) then
    bb = 4.0D+00*x*y-4.0D+00*y*y
    tbx = 4.0D+00*y
    tby = 4.0D+00*x-8.0D+00*y
  else if ( in == 6) then
    bb = 4.0D+00*x-4.0D+00*x*x-4.0D+00*y+4.0D+00*x*y
    tbx = 4.0D+00-8.0D+00*x+4.0D+00*y
    tby = -4.0D+00+4.0D+00*x
  else
    write(*,*)'REFQBF - Illegal value of IN = ',in
    stop
  end if

  bx = tbx * xix + tby * etax
  by = tbx * xiy + tby * etay

  return
end
subroutine resid ( area, g, indx, insc, isotri, iwrite, nelemn, neqn, &
  nnodes, node, np, nquad, phi, psi, res, reynld, xc, xm, yc, ym )

!*****************************************************************************80
!
!! RESID computes the residual.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) AREA(NELEMN), the area of the elements.
!
!    Input, integer ( kind = 4 ) INDX(MAXNP,2), contains, for each node I, the
!    index of the U and V velocities at that node, or 0.
!
!    Input, integer ( kind = 4 ) NNODES, the number of nodes per element, which
!    is 6 for the quadratic triangular elements in use here.
!
  implicit none

  integer ( kind = 4 ) nelemn
  integer ( kind = 4 ) neqn
  integer ( kind = 4 ) nnodes
  integer ( kind = 4 ) np
  integer ( kind = 4 ) nquad

  real ( kind = 8 ) aij
  real ( kind = 8 ) ar
  real ( kind = 8 ) area(nelemn)
  real ( kind = 8 ) bb
  real ( kind = 8 ) bbb
  real ( kind = 8 ) bbbl
  real ( kind = 8 ) bbl
  real ( kind = 8 ) bbx
  real ( kind = 8 ) bby
  real ( kind = 8 ) bx
  real ( kind = 8 ) by
  real ( kind = 8 ) det
  real ( kind = 8 ) etax
  real ( kind = 8 ) etay
  real ( kind = 8 ) g(neqn)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibad
  integer ( kind = 4 ) ihor
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) indx(np,2)
  integer ( kind = 4 ) insc(np)
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) ipp
  integer ( kind = 4 ) iprs
  integer ( kind = 4 ) iq
  integer ( kind = 4 ) iqq
  integer ( kind = 4 ) iquad
  integer ( kind = 4 ) isotri(nelemn)
  integer ( kind = 4 ) it
  integer ( kind = 4 ) itype
  integer ( kind = 4 ) iver
  integer ( kind = 4 ) iwrite
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jp
  integer ( kind = 4 ) ju
  integer ( kind = 4 ) jv
  integer ( kind = 4 ) node(nelemn,nnodes)
  real ( kind = 8 ) phi(nelemn,nquad,nnodes,3)
  real ( kind = 8 ) psi(nelemn,nquad,nnodes)
  real ( kind = 8 ) res(neqn)
  real ( kind = 8 ) reynld
  real ( kind = 8 ) rmax
  real ( kind = 8 ) test
  real ( kind = 8 ) ubc
  real ( kind = 8 ) ubdry
  real ( kind = 8 ) ubump
  real ( kind = 8 ) un(2)
  real ( kind = 8 ) unx(2)
  real ( kind = 8 ) uny(2)
  real ( kind = 8 ) visc
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) xix
  real ( kind = 8 ) xiy
  real ( kind = 8 ) xm(nelemn,nquad)
  real ( kind = 8 ) xq
  real ( kind = 8 ) yc(np)
  real ( kind = 8 ) ym(nelemn,nquad)
  real ( kind = 8 ) yq

  itype = -1
  visc = 1.0D+00 / reynld

  res(1:neqn) = 0.0D+00

  do it = 1, nelemn
    ar = area(it) / 3.0D+00
    do iquad = 1, nquad
      yq = ym(it,iquad)
      xq = xm(it,iquad)
      if ( isotri(it) == 1 ) then
        call trans(det,etax,etay,it,nelemn,nnodes,node,np,xc,xix,xiy,xq,yc,yq)
        ar = det * area(it) / 3.0D+00
      end if

      call uval(etax,etay,g,indx,isotri,it,nelemn,neqn, &
        nnodes,node,np,un,uny,unx,xc,xix,xiy,xq,yc,yq)
!
!  For each basis function:
!
      do iq = 1, nnodes

        ip = node(it,iq)
        bb = phi(it,iquad,iq,1)
        bx = phi(it,iquad,iq,2)
        by = phi(it,iquad,iq,3)
        bbl = psi(it,iquad,iq)
        iprs = insc(ip)
        ihor = indx(ip,1)
        iver = indx(ip,2)

        if ( 0 < ihor ) then
          res(ihor) = res(ihor)-ar*bb*(un(1)*unx(1)+un(2)*uny(1))
        end if

        if ( 0 < iver ) then
          res(iver) = res(iver)-ar*bb*(un(1)*unx(2)+un(2)*uny(2))
        end if
!
!  For another basis function,
!
        do iqq = 1, nnodes

          ipp = node(it,iqq)
          bbb = phi(it,iquad,iqq,1)
          bbx = phi(it,iquad,iqq,2)
          bby = phi(it,iquad,iqq,3)
          bbbl = psi(it,iquad,iqq)
          ju = indx(ipp,1)
          jv = indx(ipp,2)
          jp = insc(ipp)
!
!  Horizontal velocity variable
!
          if ( 0 < ju ) then

            if ( 0 < ihor ) then
              res(ihor) = res(ihor)+ar*(visc*(by*bby+bx*bbx) &
                +bb*(bbb*unx(1)+bbx*un(1)+bby*un(2)))*g(ju)
            end if

            if ( 0 < iver ) then
              res(iver) = res(iver)+ar*bb*bbb*unx(2)*g(ju)
            end if

            if ( 0 < iprs ) then
              res(iprs) = res(iprs)+ar*bbx*bbl*g(ju)
            end if

          else if ( ju == itype ) then

            if ( ju == -2 ) then
              ubc = ubump(g,indx,ipp,iqq,isotri,it,1,nelemn, &
                neqn,nnodes,node,np,xc,yc)
            else if ( ju == -1 ) then
              ubc = ubdry(1,yc(ipp))
            end if

            if ( 0 < ihor ) then
              aij = ar*(visc*(by*bby+bx*bbx) &
                +bb*(bbb*unx(1)+bbx*un(1)+bby*un(2)))
              res(ihor) = res(ihor)+ubc*aij
            end if

            if ( 0 < iver ) then
              aij = ar*bb*bbb*unx(2)
              res(iver) = res(iver)+ubc*aij
            end if

            if ( 0 < iprs ) then
              aij = ar*bbx*bbl
              res(iprs) = res(iprs)+ubc*aij
            end if

          end if
!
!  Vertical velocity variable
!
          if ( 0 < jv ) then

            if ( 0 < ihor ) then
              res(ihor) = res(ihor)+ar*bb*bbb*uny(1)*g(jv)
            end if

            if ( 0 < iver ) then
              res(iver) = res(iver)+ar*(visc*(by*bby+bx*bbx) &
                +bb*(bbb*uny(2)+bby*un(2)+bbx*un(1)))*g(jv)
            end if

            if ( 0 < iprs ) then
              res(iprs) = res(iprs)+ar*bby*bbl*g(jv)
            end if

          else if ( jv == itype ) then

            if ( jv == -2 ) then
              ubc = ubump(g,indx,ipp,iqq,isotri,it,2,nelemn, &
                neqn,nnodes,node,np,xc,yc)
            else if ( jv == -1 ) then
              ubc = ubdry(2,yc(ipp))
            end if

            if ( 0 < ihor ) then
              aij = ar*bb*bbb*uny(1)
              res(ihor) = res(ihor)+ubc*aij
            end if

            if ( 0 < iver ) then
              aij = ar*(visc*(by*bby+bx*bbx) &
                +bb*(bbb*uny(2)+bby*un(2)+bbx*un(1)))
              res(iver) = res(iver)+ubc*aij
            end if

            if ( 0 < iprs ) then
              aij = ar*bby*bbl
              res(iprs) = res(iprs)+ubc*aij
            end if

          end if
!
!  Pressure variable
!
          if ( 0 < jp ) then

            if ( 0 < ihor ) then
              res(ihor) = res(ihor)-ar*bx*bbbl*g(jp)
            end if

            if ( 0 < iver ) then
              res(iver) = res(iver)-ar*by*bbbl*g(jp)
            end if

          end if

        end do
      end do
    end do
  end do
!
!  The last equation is "reset" to require that the last pressure
!  be zero.
!
  res(neqn) = g(neqn)

  rmax = 0.0D+00
  imax = 0
  ibad = 0

  do i = 1, neqn

    test = abs(res(i))

    if ( rmax < test ) then
      rmax = test
      imax = i
    end if

    if ( 1.0D-03 < test ) then
      ibad = ibad+1
    end if

  end do

  if ( 1 <= iwrite ) then
    write(*,*)' '
    write(*,*)'RESIDUAL INFORMATION:'
    write(*,*)' '
    write(*,*)'Worst residual is number ',IMAX
    write(*,*)'of magnitude ',RMAX
    write(*,*)' '
    write(*,*)'Number of "bad" residuals is ',IBAD,' out of ',NEQN
    write(*,*)' '
  end if

  if ( 2 <= iwrite ) then
    write(*,*)'Raw residuals:'
    write(*,*)' '
    i = 0
    do j = 1, np

      if ( 0 < indx(j,1) ) then
        i = i+1
        if ( abs(res(i)) <= 1.0D-03 ) then
          write(*,'(1x,a1,2i5,g14.6)')'U',i,j,res(i)
        else
          write(*,'(a1,a1,2i5,g14.6)')'*','U',i,j,res(i)
        end if
      end if

      if ( 0 < indx(j,2) ) then
        i = i+1
        if ( abs(res(i)) <= 1.0D-03 ) then
          write(*,'(1x,a1,2i5,g14.6)')'V',i,j,res(i)
        else
          write(*,'(a1,a1,2i5,g14.6)')'*','V',i,j,res(i)
        end if
      end if

      if ( 0 < insc(j) ) then
        i = i+1
        if ( abs(res(i)) <= 1.0D-03 ) then
          write(*,'(1x,a1,2i5,g14.6)')'P',i,j,res(i)
        else
          write(*,'(a1,a1,2i5,g14.6)')'*','P',i,j,res(i)
        end if
      end if

    end do

  end if
  return
end
subroutine setban ( indx, insc, maxrow, nband, nelemn, nlband, nnodes, &
  node, np, nrow )

!*****************************************************************************80
!
!! SETBAN computes the half band width.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) INDX(MAXNP,2), contains, for each node I, the
!    index of the U and V velocities at that node, or 0.
!
!    Input, integer ( kind = 4 ) NNODES, the number of nodes per element, which
!    is 6 for the quadratic triangular elements in use here.
!
  implicit none

  integer ( kind = 4 ) nelemn
  integer ( kind = 4 ) nnodes
  integer ( kind = 4 ) np

  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(np,2)
  integer ( kind = 4 ) insc(np)
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) ipp
  integer ( kind = 4 ) iq
  integer ( kind = 4 ) iqq
  integer ( kind = 4 ) it
  integer ( kind = 4 ) iuk
  integer ( kind = 4 ) iukk
  integer ( kind = 4 ) j
  integer ( kind = 4 ) maxrow
  integer ( kind = 4 ) nband
  integer ( kind = 4 ) nlband
  integer ( kind = 4 ) node(nelemn,nnodes)
  integer ( kind = 4 ) nrow

  nlband = 0

  do it = 1, nelemn
    do iq = 1, nnodes
      ip = node(it,iq)
      do iuk = 1, 3
        if (iuk == 3) then
          i = insc(ip)
        else
          i = indx(ip,iuk)
        end if
        if ( 0 < i ) then
          do iqq = 1, nnodes
            ipp = node(it,iqq)
            do iukk = 1, 3
              if (iukk == 3) then
                j = insc(ipp)
              else
                j = indx(ipp,iukk)
              end if
              if ( 0 < j ) then
                nlband = max(nlband,j-i)
              end if
            end do
          end do
        end if
      end do
    end do
  end do

  nband = nlband+nlband+1
  nrow = nlband+nlband+nlband+1

  write(*,*)' '
  write(*,*)'SETBAN:'
  write(*,*)' '
  write(*,*)'  Lower bandwidth = ',nlband
  write(*,*)'  Total bandwidth = ',nband
  write(*,*)'  Required matrix rows = ',nrow

  if ( maxrow < nrow ) then
    write(*,*)'SETBAN - NROW is too large!'
    write(*,*)'The maximum allowed is ',maxrow
    stop
  end if

  return
end
subroutine setbas ( isotri, nelemn, nnodes, node, np, nquad, phi, psi, &
  xc, xm, yc, ym )

!*****************************************************************************80
!
!! SETBAS evaluates the basis functions at each integration point.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NNODES, the number of nodes per element, which
!    is 6 for the quadratic triangular elements in use here.
!
  implicit none

  integer ( kind = 4 ) nelemn
  integer ( kind = 4 ) nnodes
  integer ( kind = 4 ) np
  integer ( kind = 4 ) nquad

  real ( kind = 8 ) bb
  real ( kind = 8 ) bsp
  real ( kind = 8 ) bx
  real ( kind = 8 ) by
  real ( kind = 8 ) det
  real ( kind = 8 ) etax
  real ( kind = 8 ) etay
  integer ( kind = 4 ) iq
  integer ( kind = 4 ) isotri(nelemn)
  integer ( kind = 4 ) it
  integer ( kind = 4 ) j
  integer ( kind = 4 ) node(nelemn,nnodes)
  real ( kind = 8 ) phi(nelemn,nquad,nnodes,3)
  real ( kind = 8 ) psi(nelemn,nquad,nnodes)
  real ( kind = 8 ) refbsp
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) xix
  real ( kind = 8 ) xiy
  real ( kind = 8 ) xm(nelemn,nquad)
  real ( kind = 8 ) xq
  real ( kind = 8 ) yc(np)
  real ( kind = 8 ) ym(nelemn,nquad)
  real ( kind = 8 ) yq

  do it = 1, nelemn
    do j = 1, nquad
      xq = xm(it,j)
      yq = ym(it,j)
      call trans(det,etax,etay,it,nelemn,nnodes,node,np,xc,xix,xiy,xq,yc,yq)
      do iq = 1, nnodes
        if ( isotri(it) == 0 ) then
          psi(it,j,iq) = bsp(it,iq,1,nelemn,nnodes,node,np,xc,xq,yc,yq)
          call qbf(xq,yq,it,iq,bb,bx,by,nelemn,nnodes,node,np,xc,yc)
        else
          call refqbf(xq,yq,iq,bb,bx,by,etax,etay,xix,xiy)
          psi(it,j,iq) = refbsp(xq,yq,iq)
        end if
        phi(it,j,iq,1) = bb
        phi(it,j,iq,2) = bx
        phi(it,j,iq,3) = by
      end do
    end do
  end do

  return
end
subroutine setgrd ( ibump, indx, insc, isotri, iwrite, long, maxeqn, mx, my, &
  nelemn, neqn, nnodes, node, np, nx, ny, xbleft, xbrite, xlngth )

!*****************************************************************************80
!
!! SETGRD sets up the geometric grid.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) INDX(MAXNP,2), contains, for each node I, the
!    index of the U and V velocities at that node, or 0.
!
!    Input, integer ( kind = 4 ) MX, the number of nodes in the X direction.
!
!    Input, integer ( kind = 4 ) MY, the number of nodes in the Y direction.
!
!    Input, integer ( kind = 4 ) NNODES, the number of nodes per element, which
!    is 6 for the quadratic triangular elements in use here.
!
  implicit none

  integer ( kind = 4 ) my
  integer ( kind = 4 ) nelemn
  integer ( kind = 4 ) nnodes
  integer ( kind = 4 ) np

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibump
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) icnt
  integer ( kind = 4 ) ielemn
  integer ( kind = 4 ) indx(np,2)
  integer ( kind = 4 ) insc(np)
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) ip1
  integer ( kind = 4 ) ip2
  integer ( kind = 4 ) isotri(nelemn)
  integer ( kind = 4 ) it
  integer ( kind = 4 ) iwrite
  integer ( kind = 4 ) jc
  integer ( kind = 4 ) jcnt
  logical long
  integer ( kind = 4 ) maxeqn
  integer ( kind = 4 ) mx
  integer ( kind = 4 ) nbleft
  integer ( kind = 4 ) nbrite
  integer ( kind = 4 ) neqn
  integer ( kind = 4 ) node(nelemn,nnodes)
  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny
  real ( kind = 8 ) xbleft
  real ( kind = 8 ) xbrite
  real ( kind = 8 ) xlngth

  write(*,*)' '
  write(*,*)'SETGRD:'
  write(*,*)' '
!
!  Determine whether region is long or skinny.  This will determine
!  how we number the nodes and elements.
!
  if ( ny < nx ) then
    long = .true.
    write(*,*)'Using vertical ordering.'
  else
    long = .false.
    write(*,*)'Using horizontal ordering.'
  end if
!
!  Report on isoparametric strategy:
!
  if ( ibump == 0 ) then
    write(*,*)'No isoparametric elements will be used.'
  else if ( ibump == 1 ) then
    write(*,*)'Isoparametric elements directly on bump.'
  else if ( ibump == 2 ) then
    write(*,*)'All elements above bump are isoparametric.'
  else if ( ibump == 3 ) then
    write(*,*)'All elements are isoparametric.'
  else
    write(*,*)'Unexpected value of IBUMP = ',ibump
    stop
  end if
!
!  Compute node locations of bump corners based on X coordinates
!
  nbleft = nint(xbleft*(mx-1)/xlngth)+1
  nbrite = nint(xbrite*(mx-1)/xlngth)+1
  write(*,*)'Bump extends from ',xbleft,' at node ',nbleft
  write(*,*)'               to ',xbrite,' at node ',nbrite
!
!  Assign nodes to elements.
!
  neqn = 0
  ielemn = 0

  do ip = 1, np

    if ( long ) then
      ic = ((ip-1)/my)+1
      jc = mod((ip-1),my)+1
    else
      ic = mod((ip-1),mx)+1
      jc = ((ip-1)/mx)+1
    end if

    icnt = mod(ic,2)
    jcnt = mod(jc,2)
!
!  If both the row count and the column count are odd,
!  and we're not in the last row or top column,
!  then we can define two new triangular elements based at the node.
!
!  For horizontal ordering,
!  given the following arrangement of nodes, for instance:
!
!    21 22 23 24 25
!    16 17 18 19 20
!    11 12 13 14 15
!    06 07 08 09 10
!    01 02 03 04 05
!
!  when we arrive at node 13, we will define
!
!  element 7: (13, 23, 25, 18, 24, 19)
!  element 8: (13, 25, 15, 19, 20, 14)
!
!
!  For vertical ordering,
!  given the following arrangement of nodes, for instance:
!
!    05 10 15 20 25
!    04 09 14 19 24
!    03 08 13 18 23
!    02 07 12 17 22
!    01 06 11 16 21
!
!  when we arrive at node 13, we will define
!
!  element 7: (13, 25, 23, 19, 24, 18)
!  element 8: (13, 15, 25, 14, 20, 19)
!
    if ( (icnt == 1.and.jcnt == 1).and.(ic /= mx).and.(jc /= my) ) then

      if ( long ) then

        ip1 = ip+my
        ip2 = ip+my+my
        ielemn = ielemn+1
        node(ielemn,1) = ip
        node(ielemn,2) = ip+2
        node(ielemn,3) = ip2+2
        node(ielemn,4) = ip+1
        node(ielemn,5) = ip1+2
        node(ielemn,6) = ip1+1

        if ( ibump == 0 ) then
          isotri(ielemn) = 0
        else if ( ibump == 1 ) then
          isotri(ielemn) = 0
        else if ( ibump == 2 ) then
          if ( nbleft <= ic .and. ic < nbrite ) then
            isotri(ielemn) = 1
          else
            isotri(ielemn) = 0
          end if
        else
          isotri(ielemn) = 1
        end if

        ielemn = ielemn+1
        node(ielemn,1) = ip
        node(ielemn,2) = ip2+2
        node(ielemn,3) = ip2
        node(ielemn,4) = ip1+1
        node(ielemn,5) = ip2+1
        node(ielemn,6) = ip1

        if ( ibump == 0 ) then
          isotri(ielemn) = 0
        else if ( ibump == 1 ) then
          if ( jc == 1 .and. nbleft <= ic .and. ic < nbrite ) then
            isotri(ielemn) = 1
          else
            isotri(ielemn) = 0
          end if
        else if ( ibump == 2 ) then
          if ( nbleft <= ic .and. ic < nbrite ) then
            isotri(ielemn) = 1
          else
            isotri(ielemn) = 0
          end if
        else
          isotri(ielemn) = 1
        end if

      else

        ip1 = ip+mx
        ip2 = ip+mx+mx

        ielemn = ielemn+1
        node(ielemn,1) = ip
        node(ielemn,2) = ip2
        node(ielemn,3) = ip2+2
        node(ielemn,4) = ip1
        node(ielemn,5) = ip2+1
        node(ielemn,6) = ip1+1

        if ( ibump == 0 ) then
          isotri(ielemn) = 0
        else if ( ibump == 1 ) then
          isotri(ielemn) = 0
        else if ( ibump == 2 ) then
          if ( nbleft <= ic .and. ic < nbrite ) then
            isotri(ielemn) = 1
          else
            isotri(ielemn) = 0
          end if
        else
          isotri(ielemn) = 1
        end if

        ielemn = ielemn+1
        node(ielemn,1) = ip
        node(ielemn,2) = ip2+2
        node(ielemn,3) = ip+2
        node(ielemn,4) = ip1+1
        node(ielemn,5) = ip1+2
        node(ielemn,6) = ip+1

        if ( ibump == 0 ) then
          isotri(ielemn) = 0
        else if ( ibump == 1 ) then
          if ( jc == 1 .and. nbleft <= ic .and. ic < nbrite ) then
            isotri(ielemn) = 1
          else
            isotri(ielemn) = 0
          end if
        else if ( ibump == 2 ) then
          if ( nbleft <= ic .and. ic < nbrite ) then
            isotri(ielemn) = 1
          else
            isotri(ielemn) = 0
          end if
        else
          isotri(ielemn) = 1
        end if

      end if
    end if
!
!  Left hand boundary, horizontal and vertical velocities specified.
!
    if ( ic == 1.and.1 < jc .and. jc < my ) then
      indx(ip,1) = -1
      indx(ip,2) = -1
!
!  Right hand boundary, horizontal velocities unknown, vertical
!  velocities specified
!
    else if ( ic == mx.and.1 < jc .and. jc < my) then
      neqn = neqn+1
      indx(ip,1) = neqn
      indx(ip,2) = 0
!
!  Lower boundary, with isoperimetric triangle
!
    else if ( jc == 1 .and. isotri(ielemn) == 1 ) then
      indx(ip,1) = -2
      indx(ip,2) = -2
!
!  Otherwise, just a wall
!
    else if ( ic == 1 .or. ic == mx .or. jc == 1 .or. jc == my ) then
      indx(ip,1) = 0
      indx(ip,2) = 0
!
!  Otherwise, a normal interior node where both velocities are unknown
!
    else
      neqn = neqn+2
      indx(ip,1) = neqn-1
      indx(ip,2) = neqn
    end if

    if ( jcnt == 1 .and. icnt == 1 ) then
      neqn = neqn+1
      insc(ip) = neqn
    else
      insc(ip) = 0
    end if

  end do

  if ( 1 <= iwrite ) then
    write(*,*)' '
    write(*,*)'     I     INDX 1, INDX 2, INSC'
    write(*,*)' '
    do i = 1, np
      write(*,'(4i5)')i,indx(i,1),indx(i,2),insc(i)
    end do
    write(*,*)' '
    write(*,*)'Isoparametric triangles:'
    write(*,*)' '
    do i = 1, nelemn
      if ( isotri(i) == 1)write(*,*)i
    end do
    write(*,*)' '
    write(*,*)'   IT   NODE(IT,*)'
    write(*,*)' '
    do it = 1, nelemn
      write(*,'(7i6)') it,(node(it,i),i = 1,6)
    end do
  end if

  write(*,*)' '
  write(*,*)'SETGRD: Number of unknowns = ',neqn

  if ( maxeqn < neqn ) then
    write(*,*)'SETGRD - Too many unknowns!'
    write(*,*)'The maximum allowed is MAXEQN = ',maxeqn
    write(*,*)'This problem requires NEQN = ',neqn
    stop
  end if

  return
end
subroutine setlin ( iline, indx, iwrite, long, mx, my, np, nx, ny, &
  xlngth, xprof )

!*****************************************************************************80
!
!! SETLIN determines unknown numbers along the profile line.
!
!  Discussion:
!
!    For our problem, the profile line has the equation X = XPROF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) INDX(MAXNP,2), contains, for each node I, the
!    index of the U and V velocities at that node, or 0.
!
!    Input, integer ( kind = 4 ) MX, the number of nodes in the X direction.
!
!    Input, integer ( kind = 4 ) MY, the number of nodes in the Y direction.
!
  implicit none

  integer ( kind = 4 ) my
  integer ( kind = 4 ) np

  integer ( kind = 4 ) i
  integer ( kind = 4 ) iline(my)
  integer ( kind = 4 ) indx(np,2)
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) itemp
  integer ( kind = 4 ) iwrite
  logical long
  integer ( kind = 4 ) mx
  integer ( kind = 4 ) nodex0
  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny
  real ( kind = 8 ) xlngth
  real ( kind = 8 ) xprof
!
!  Determine the number of a node on the profile line.
!
  itemp = nint ( ( 2.0D+00 * real ( nx - 1, kind = 8 ) * xprof ) / xlngth )

  if ( long ) then
    nodex0 = itemp * ( 2 * ny - 1 ) + 1
  else
    nodex0 = itemp + 1
  end if

  write(*,*)' '
  write(*,*)'SETLIN:'
  write(*,*)' '
  write(*,*)'  Profile generated at X = ',xprof
  write(*,*)'  which is above node  = ',nodex0

  do i = 1, my
    if ( long ) then
      ip = nodex0+(i-1)
    else
      ip = nodex0+mx*(i-1)
    end if
    iline(i) = indx(ip,1)
  end do

  if ( 1 <= iwrite ) then
    write(*,*)' '
    write(*,*)'  Indices of unknowns along the profile line:'
    write(*,*)' '
    write(*,'(5i5)') iline(1:my)
  end if

  return
end
subroutine setqud ( area, isotri, iwrite, nelemn, nnodes, node, np, nquad, &
  xc, xm, yc, ym )

!*****************************************************************************80
!
!! SETQUD sets midpoint quadrature rule information.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) AREA(NELEMN), the area of the elements.
!
!    Input, integer ( kind = 4 ) NNODES, the number of nodes per element, which
!    is 6 for the quadratic triangular elements in use here.
!
  implicit none

  integer ( kind = 4 ) nelemn
  integer ( kind = 4 ) nnodes
  integer ( kind = 4 ) np
  integer ( kind = 4 ) nquad

  real ( kind = 8 ) area(nelemn)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ip1
  integer ( kind = 4 ) ip2
  integer ( kind = 4 ) ip3
  integer ( kind = 4 ) isotri(nelemn)
  integer ( kind = 4 ) it
  integer ( kind = 4 ) iwrite
  integer ( kind = 4 ) j
  integer ( kind = 4 ) node(nelemn,nnodes)
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) xm(nelemn,nquad)
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) y3
  real ( kind = 8 ) yc(np)
  real ( kind = 8 ) ym(nelemn,nquad)

  do it = 1, nelemn

    ip1 = node(it,1)
    ip2 = node(it,2)
    ip3 = node(it,3)
    x1 = xc(ip1)
    x2 = xc(ip2)
    x3 = xc(ip3)
    y1 = yc(ip1)
    y2 = yc(ip2)
    y3 = yc(ip3)

    if ( isotri(it) == 0 ) then
      xm(it,1) = 0.5D+00*(x1+x2)
      xm(it,2) = 0.5D+00*(x2+x3)
      xm(it,3) = 0.5D+00*(x3+x1)
      ym(it,1) = 0.5D+00*(y1+y2)
      ym(it,2) = 0.5D+00*(y2+y3)
      ym(it,3) = 0.5D+00*(y3+y1)
      area(it) = 0.5D+00*abs((y1+y2)*(x2-x1)+(y2+y3)*(x3-x2)+(y3+y1)*(x1-x3))
    else
      xm(it,1) = 0.5D+00
      ym(it,1) = 0.5D+00
      xm(it,2) = 1.0D+00
      ym(it,2) = 0.5D+00
      xm(it,3) = 0.5D+00
      ym(it,3) = 0.0D+00
      area(it) = 0.5D+00
    end if

  end do

  if ( 3 <= iwrite ) then
    write(*,*)' '
    write(*,*)'SETQUD: Element Areas and Quadrature points:'
    write(*,*)' '
    do i = 1, nelemn
      write(*,*)i,area(i)
      do j = 1, nquad
        write(*,*)i,j,xm(i,j),ym(i,j)
      end do
    end do
  end if

  return
end
subroutine setxy ( iwrite, long, mx, my, np, nx, ny, xc, xlngth, yc, &
  ylngth, ypert )

!*****************************************************************************80
!
!! SETXY sets the grid coordinates based on the value of the parameter.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) MX, the number of nodes in the X direction.
!
!    Input, integer ( kind = 4 ) MY, the number of nodes in the Y direction.
!
  implicit none

  integer ( kind = 4 ) my
  integer ( kind = 4 ) np

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) iwrite
  integer ( kind = 4 ) jc
  logical long
  integer ( kind = 4 ) mx
  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) xlngth
  real ( kind = 8 ) ybot
  real ( kind = 8 ) yc(np)
  real ( kind = 8 ) ylngth
  real ( kind = 8 ) ylo
  real ( kind = 8 ) ypert

  do ip = 1, np

    if ( long ) then
      ic = ((ip-1)/my)+1
      jc = mod((ip-1),my)+1
    else
      ic = mod((ip-1),mx)+1
      jc = ((ip-1)/mx)+1
    end if

    xc(ip) = (ic-1) * xlngth / (2*nx-2)

    ybot = -ypert * ( xc(ip) - 3.0D+00 ) * ( xc(ip) - 1.0D+00 )
    ylo = max ( 0.0D+00, ybot )

    yc(ip) = ( (my-jc)*ylo + (jc-1)*ylngth ) / (2*ny-2)

  end do

  if ( 2 <= iwrite ) then
    write(*,*)' '
    write(*,*)'SETGRD:'
    write(*,*)' '
    write(*,*)'     I     XC     YC'
    write(*,*)' '
    do i = 1, np
      write(*,'(i5,2f12.5)')i,xc(i),yc(i)
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
!    May 31 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 February 2005
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

  write ( *, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
subroutine trans ( det, etax, etay, it, nelemn, nnodes, node, np, xc, &
  xix, xiy, xq, yc, yq )

!*****************************************************************************80
!
!! TRANS calculates the element transformation mapping.
!
!  Discussion:
!
!    The element transformation mapping maps the reference element
!    to a particular isoparametric element.  The routine also outputs
!    the determinant and partial derivatives of the transformation.
!
!    Diagram of the mapping from reference to isoparametric element:
!
!          2                           2
!    E    /|                          / \
!    t   4 5       ===== >     Y     4   5
!    a  /  |                        /     \
!      1-6-3                       1--6----3
!
!      Xi                             X
!
!    The form of the quadratic mapping is:
!
!      x = a1 * xi^2 + b1 * xi * eta + c1 * eta^2 + d1 * xi + e1 * eta + f1
!      y = a2 * xi^2 + b2 * xi * eta + c2 * eta^2 + d2 * xi + e2 * eta + f2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) DET, the determinant of the transformation,
!    evaluated at (XQ,YQ).
!
!    Output, real ( kind = 8 ) ETAX, ETAY, the Jacobian matrix entries
!    dETA/dX and dETA/dY evaluated at (XQ,YQ).
!
!    Input, integer ( kind = 4 ) IT, the index of the element containing (XQ,YQ).
!
!    Input, integer ( kind = 4 ) NELEMN, the number of elements.
!
!    Input, integer ( kind = 4 ) NNODES, the number of nodes per element, which
!    is 6 for the quadratic triangular elements in use here.
!
!    Input, integer ( kind = 4 ) NODE(NELEMN,NNODES), contains the indices of
!    the nodes that make up each element.
!
!    Input, integer ( kind = 4 ) NP, the number of nodes.
!
!    Input, real ( kind = 8 ) XC(NP), the X coordinates of nodes.
!
!    Output, real ( kind = 8 ) XIX, XIY, the Jacobian matrix entries
!    dXI/dX and dXI/dY evaluated at (XQ,YQ).
!
!    Input, real ( kind = 8 ) XQ, the X coordinate of the point at
!    which the mapping is to be evaluated.
!
!    Input, real ( kind = 8 ) YC(NP), the Y coordinates of nodes.
!
!    Input, real ( kind = 8 ) YQ, the Y coordinate of the point at
!    which the mapping is to be evaluated.
!
  implicit none

  integer ( kind = 4 ) nelemn
  integer ( kind = 4 ) nnodes
  integer ( kind = 4 ) np

  real ( kind = 8 ) a1
  real ( kind = 8 ) a2
  real ( kind = 8 ) b1
  real ( kind = 8 ) b2
  real ( kind = 8 ) c1
  real ( kind = 8 ) c2
  real ( kind = 8 ) d1
  real ( kind = 8 ) d2
  real ( kind = 8 ) det
  real ( kind = 8 ) dxdeta
  real ( kind = 8 ) dxdxi
  real ( kind = 8 ) dydeta
  real ( kind = 8 ) dydxi
  real ( kind = 8 ) e1
  real ( kind = 8 ) e2
  real ( kind = 8 ) etax
  real ( kind = 8 ) etay
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) i4
  integer ( kind = 4 ) i5
  integer ( kind = 4 ) i6
  integer ( kind = 4 ) it
  integer ( kind = 4 ) node(nelemn,nnodes)
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ) x4
  real ( kind = 8 ) x5
  real ( kind = 8 ) x6
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) xix
  real ( kind = 8 ) xiy
  real ( kind = 8 ) xq
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) y3
  real ( kind = 8 ) y4
  real ( kind = 8 ) y5
  real ( kind = 8 ) y6
  real ( kind = 8 ) yc(np)
  real ( kind = 8 ) yq

  i1 = node(it,1)
  i2 = node(it,2)
  i3 = node(it,3)
  i4 = node(it,4)
  i5 = node(it,5)
  i6 = node(it,6)
  x1 = xc(i1)
  y1 = yc(i1)
  x2 = xc(i2)
  y2 = yc(i2)
  x3 = xc(i3)
  y3 = yc(i3)
  x4 = xc(i4)
  y4 = yc(i4)
  x5 = xc(i5)
  y5 = yc(i5)
  x6 = xc(i6)
  y6 = yc(i6)
!
!  Set the coefficients in the transformation:
!
!  X = X(XI,ETA)
!  Y = Y(XI,ETA)
!
  a1 = 2.0D+00*x3-4.0D+00*x6+2.0D+00*x1
  b1 = -4.0D+00*x3-4.0D+00*x4+4.0D+00*x5+4.0D+00*x6
  c1 = 2.0D+00*x2+2.0D+00*x3-4.0D+00*x5
  d1 = -3.0D+00*x1-x3+4.0D+00*x6
  e1 = -x2+x3+4.0D+00*x4-4.0D+00*x6

  a2 = 2.0D+00*y3-4.0D+00*y6+2.0D+00*y1
  b2 = -4.0D+00*y3-4.0D+00*y4+4.0D+00*y5+4.0D+00*y6
  c2 = 2.0D+00*y2+2.0D+00*y3-4.0D+00*y5
  d2 = -3.0D+00*y1-y3+4.0D+00*y6
  e2 = -y2+y3+4.0D+00*y4-4.0D+00*y6
!
!  Compute partial derivatives d x/deta, d x/dxi, d y/deta/ dy/d xi,
!  at point (xq,yq) in reference triangle.  This is the Jacobian:
!
!   ( dx/dxi    dx/deta )
!   ( dy/dxi    dy/deta )
!
  dxdxi = 2.0D+00*a1*xq+b1*yq+d1
  dxdeta = b1*xq+2.0D+00*c1*yq+e1
  dydxi = 2.0D+00*a2*xq+b2*yq+d2
  dydeta = b2*xq+2.0D+00*c2*yq+e2
!
!  Compute the determinant of the transformation.
!
  det = (2.0D+00*a1*b2-2.0D+00*a2*b1)*xq*xq &
    +(4.0D+00*a1*c2-4.0D+00*a2*c1)*xq*yq &
    +(2.0D+00*b1*c2-2.0D+00*b2*c1)*yq*yq &
    +(2.0D+00*a1*e2+b2*d1-b1*d2-2.0D+00*a2*e1)*xq &
    +(2.0D+00*c2*d1+b1*e2-b2*e1-2.0D+00*c1*d2)*yq+d1*e2-d2*e1
!
!  Compute the inverse jacobian 
!
!    ( dxi/dx   dxi/dy  )
!    ( deta/dx  deta/dy )
!
  xix  =   dydeta  / det
  xiy  = - dxdeta / det
  etax = - dydxi / det
  etay =   dxdxi / det

  return
end
function ubdry ( iuk, yy )

!*****************************************************************************80
!
!! UBDRY sets the parabolic inflow in terms of the value of the parameter
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IUK, the index of the unknown.
!    1, the horizontal velocity.
!    2, the vertical velocity.
!
!    Input, real ( kind = 8 ) YY, the Y coordinate of the boundary point.
!
!    Output, real ( kind = 8 ) UBDRY, the value of the prescribed boundary
!    flow component at this coordinate on the boundary.
!
  implicit none

  integer ( kind = 4 ) iuk
  real ( kind = 8 ) ubdry
  real ( kind = 8 ) yy

  if ( iuk == 1 ) then
    ubdry = ( -2.0D+00 * yy + 6.0D+00 ) * yy / 9.0D+00
  else
    ubdry = 0.0D+00
  end if

  return
end
function ubump ( g, indx, ip, iqq, isotri, it, iukk, nelemn, neqn, &
  nnodes, node, np, xc, yc )

!*****************************************************************************80
!
!! UBUMP calculates the sensitivity dU/dA on the bump.
!
!  Discussion:
!
!    This routine sets
!
!      dU/dA = -uy * phi
!      dV/dA = -vy * phi
!
!    where PHI is the shape of the bump.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) INDX(MAXNP,2), contains, for each node I, the
!    index of the U and V velocities at that node, or 0.
!
  implicit none

  integer ( kind = 4 ) nelemn
  integer ( kind = 4 ) neqn
  integer ( kind = 4 ) nnodes
  integer ( kind = 4 ) np

  real ( kind = 8 ) det
  real ( kind = 8 ) etax
  real ( kind = 8 ) etay
  real ( kind = 8 ) g(neqn)
  integer ( kind = 4 ) indx(np,2)
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) iqq
  integer ( kind = 4 ) isotri(nelemn)
  integer ( kind = 4 ) it
  integer ( kind = 4 ) iukk
  integer ( kind = 4 ) node(nelemn,nnodes)
  real ( kind = 8 ) ubump
  real ( kind = 8 ) un(2)
  real ( kind = 8 ) unx(2)
  real ( kind = 8 ) uny(2)
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) xix
  real ( kind = 8 ) xiy
  real ( kind = 8 ) xq
  real ( kind = 8 ) yc(np)
  real ( kind = 8 ) yq

  if ( isotri(it) == 0 ) then
    xq = xc(ip)
    yq = yc(ip)
  else
    if ( iqq == 1 ) then
      xq = 0.0D+00
      yq = 0.0D+00
    else if ( iqq == 2 ) then
      xq = 1.0D+00
      yq = 1.0D+00
    else if ( iqq == 3 ) then
      xq = 1.0D+00
      yq = 0.0D+00
    else if ( iqq == 4 ) then
      xq = 0.5D+00
      yq = 0.5D+00
    else if ( iqq == 5 ) then
      xq = 1.0D+00
      yq = 0.5D+00
    else if ( iqq == 6 ) then
      xq = 0.5D+00
      yq = 0.0D+00
    end if
    call trans(det,etax,etay,it,nelemn,nnodes,node,np,xc,xix,xiy,xq,yc,yq)
  end if
!
!  Calculate value of uy and vy (old solutions) at node
!
  call uval(etax,etay,g,indx,isotri,it,nelemn,neqn, &
    nnodes,node,np,un,uny,unx,xc,xix,xiy,xq,yc,yq)

  if ( iukk == 1 ) then
    ubump = uny(1) * (xc(ip)-1.0D+00) * (xc(ip)-3.0D+00)
  else if ( iukk == 2 ) then
    ubump = uny(2) * (xc(ip)-1.0D+00) * (xc(ip)-3.0D+00)
  else
    write(*,*)'UBUMP called for iukk = ',iukk
    stop
  end if

  return
end
subroutine uval ( etax, etay, g, indx, isotri, it, nelemn, neqn, nnodes, &
  node, np, un, uny, unx, xc, xix, xiy, xq, yc, yq )

!*****************************************************************************80
!
!! UVAL evaluates the velocities at a given quadrature point.
!
!  Discussion:
!
!    This routine also computes the spatial derivatives of the velocities.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) INDX(MAXNP,2), contains, for each node I, the
!    index of the U and V velocities at that node, or 0.
!
!    Input, integer ( kind = 4 ) NNODES, the number of nodes per element, which
!    is 6 for the quadratic triangular elements in use here.
!
  implicit none

  integer ( kind = 4 ) nelemn
  integer ( kind = 4 ) neqn
  integer ( kind = 4 ) nnodes
  integer ( kind = 4 ) np

  real ( kind = 8 ) bb
  real ( kind = 8 ) bx
  real ( kind = 8 ) by
  real ( kind = 8 ) etax
  real ( kind = 8 ) etay
  real ( kind = 8 ) g(neqn)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(np,2)
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) iq
  integer ( kind = 4 ) isotri(nelemn)
  integer ( kind = 4 ) it
  integer ( kind = 4 ) iuk
  integer ( kind = 4 ) iun
  integer ( kind = 4 ) node(nelemn,nnodes)
  real ( kind = 8 ) ubc
  real ( kind = 8 ) ubdry
  real ( kind = 8 ) un(2)
  real ( kind = 8 ) unx(2)
  real ( kind = 8 ) uny(2)
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) xix
  real ( kind = 8 ) xiy
  real ( kind = 8 ) xq
  real ( kind = 8 ) yc(np)
  real ( kind = 8 ) yq

  un(1:2) = 0.0D+00
  unx(1:2) = 0.0D+00
  uny(1:2) = 0.0D+00

  do iq = 1, nnodes

    if ( isotri(it) == 1 ) then
      call refqbf(xq,yq,iq,bb,bx,by,etax,etay,xix,xiy)
    else
      call qbf(xq,yq,it,iq,bb,bx,by,nelemn,nnodes,node,np,xc,yc)
    end if
    ip = node(it,iq)

    do iuk = 1, 2
      iun = indx(ip,iuk)
      if ( 0 < iun ) then
        un(iuk) = un(iuk)+bb*g(iun)
        unx(iuk) = unx(iuk)+bx*g(iun)
        uny(iuk) = uny(iuk)+by*g(iun)
      else if ( iun == -1 ) then
        ubc = ubdry(iuk,yc(ip))
        un(iuk) = un(iuk)+bb*ubc
        unx(iuk) = unx(iuk)+bx*ubc
        uny(iuk) = uny(iuk)+by*ubc
      end if
    end do

  end do

  return
end
subroutine uv_write ( f, indx, uv_unit, neqn, np, yc )

!*****************************************************************************80
!
!! UV_WRITE writes a velocity file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) INDX(MAXNP,2), contains, for each node I, the
!    index of the U and V velocities at that node, or 0.
!
!    Input, integer ( kind = 4 ) UV_UNIT, the FORTRAN unit number associated with
!    the output file.
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations.
!
!    Input, integer ( kind = 4 ) NP, the number of nodes.
!
!    Input, real ( kind = 8 ) YC(NP), the Y coordinates of nodes.
!
  implicit none

  integer ( kind = 4 ) neqn
  integer ( kind = 4 ) np

  real ( kind = 8 ) f(neqn)
  integer ( kind = 4 ) indx(np,2)
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) k
  real ( kind = 8 ) u
  real ( kind = 8 ) ubdry
  integer ( kind = 4 ) uv_unit
  real ( kind = 8 ) v
  real ( kind = 8 ) yc(np)

  do ip = 1, np

    k = indx(ip,1)

    if ( k < 0 ) then
      u = ubdry ( 1, yc(ip) )
    else if ( k == 0 ) then
      u = 0.0D+00
    else
      u = f(k)
    end if

    k = indx(ip,2)

    if ( k < 0 ) then
      v = ubdry ( 2, yc(ip) )
    else if ( k == 0 ) then
      v = 0.0D+00
    else
      v = f(k)
    end if

    write ( uv_unit, '(2x,g14.6,2x,g14.6)' ) u, v

  end do

  return
end
subroutine xy_write ( xy_unit, np, xc, yc )

!*****************************************************************************80
!
!! XY_WRITE creates node coordinate data.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) XY_UNIT, the FORTRAN unit number associated with
!    the output file.
!
!    Input, integer ( kind = 4 ) NP, the number of nodes.
!
!    Input, real ( kind = 8 ) XC(NP), YC(NP), the coordinates of nodes.
!
  implicit none

  integer ( kind = 4 ) np

  integer ( kind = 4 ) ip
  real ( kind = 8 ) xc(np)
  integer ( kind = 4 ) xy_unit
  real ( kind = 8 ) yc(np)

  do ip = 1, np
    write ( xy_unit, '(2x,g14.6,2x,g14.6)' ) xc(ip), yc(ip)
  end do

  return
end
