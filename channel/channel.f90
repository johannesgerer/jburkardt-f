program main

!*****************************************************************************80
!
!! MAIN is the main program for CHANNEL.
!
!  Discussion:
!
!    CHANNEL solves the channel flow problem using the finite element method.
!
!    This program solves a fluid flow problem on the unit square.
!
!    The fluid flow problem is formulated in terms of primitive variables
!    u, v, and p.
!
!    This code tries to match a downstream profile by altering
!    one parameter, the value of the inflow parameter at the mid-height
!    of the channel.
!
!    Piecewise linear functions on triangles approximate
!    the pressure and quadratics on triangles approximate the velocity.
!    This is the "Taylor-Hood" finite element basis.
!
!    The primitive variable formulation of the Navier Stokes equations
!    involves horizontal velocity U,
!    vertical velocity V, and pressure P.  The equations are:
!
!      U dUdx + V dUdy + dPdx - mu*(ddU/dxdx + ddU/dydy) = F1
!      U dVdx + V dVdy + dPdy - mu*(ddV/dxdx + ddV/dydy) = F2
!      dUdx + dVdy = 0
!
!    When reformulated into finite element form, with PHI(i) being the I-th
!    common basis function for U and V, and PSI(i) the I-th basis function
!    for P, these equations become:
!
!    Integral (U dUdx PHI(i) + V dUdy PHI(I) - P dPHI(i)/dx
!              + mu (dUdx dPHI(i)/dx + dUdy dPHI(i)/dy) ) =
!    Integral (F1 * PHI(i))
!
!    Integral (U dVdx PHI(i) + V dVdy PHI(I) - P dPHI(i)/dy
!              + mu (dVdx dPHI(i)/dx + dVdy dPHI(i)/dy) ) =
!    Integral (F2 * PHI(i))
!
!    Integral (dUdx PSI(i) + dVdx PSI(i)) = 0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Journal:
!
!    4 September 1992
!    Added LONG to GDUMP output.
!
!    2 September 1992:
!    Set g to zero before secant iteration.
!    Made sure to update g BEFORE returning from NSTOKE.
!    New code converges in 3 secant steps, 13 Newton steps,
!    and 48 seconds.
!
!    27 August 1992: Corrected a mistake in XYDUMP that only occurred
!    for NX <= NY.
!
!    26 August 1992: broke out SETLIN, SETQUD and SETBAS for
!    similarity with BUMP.
!
!    22 August 1992: Moved bandwidth calculations to SETBAN, mainly for
!    the benefit of the BUMP computation.
!
!    20 August 1992: Moved linear system setup and solve out of NSTOKE,
!    merged it with solvlin, and called it "LINSYS".
!
!    19 August 1992.  After I changed the real ( kind = 8 ) of the matrix
!    A in NSTOKE from (MAXROW,NEQN) to (NROW,NEQN) I found performance
!    was degraded.  The time went up from about 47 seconds to 57 seconds!
!    (I made other, also theoretically harmless changes at the same time).
!    Changing back to MAXROW brought the time down to 50 seconds.  What's
!    going on?  This is definitely a performance problem in LAPACK!
!
!    18 August 1992, the residual calculation has been healed, somehow!
!    That means that I can, if I wish, think about a jacobian.
!
!    13 August 1992, the way that the last pressure was set to 1 was not
!    working properly.  I don't know why, but I changed the code so that
!    the last pressure was set to zero, and the contour plots smoothed out,
!    and the secant iteration converged in three steps, rather than four.
!    So I'm also dropping the entire "PRESET" routine, which forced the
!    average pressure to zero.  What's the point?  Perhaps that's important
!    when you're modifying the grid, though.  So I won't actually delete
!    the code.
!
!    11 August 1992, the original scheme for the last pressure equation
!    does NOT force the average pressure to be zero.  If you don't
!    believe me, repeat the calculation of PMEAN after the adjustment.
!    I've fixed that.
!
!  Variables
!
!    real ( kind = 8 ) A(MAXROW,MAXEQN), the banded matrix
!    used in NSTOKE and SOLVLIN.  It is stored according to
!    LINPACK/LAPACK general band storage mode.
!
!    real ( kind = 8 ) AREA(NELEMN), contains the area of each
!    element.
!
!    real ( kind = 8 ) F(MAXEQN).  After the call to NSTOKE,
!    F contains the solution to the current Navier Stokes problem.
!
!    character*30 FILEG, the filename of the graphics data file
!    to be dumped for the DISPLAY program.
!
!    character*30 FILEU, the filename of the graphics UV data file
!    to be dumped for the PLOT3D program.
!
!    character*30 FILEX, the filename of the graphics XY data file
!    to be dumped for the PLOT3D program.
!
!    real ( kind = 8 ) G(MAXEQN).  After the call to SOLVLIN,
!    G contains the sensitivities.
!
!    integer INDX(NP,2), records, for each node, whether
!    any velocity unknowns are associated with the node.
!    INDX(I,1) records information for horizontal velocities.
!    If it is 0, then no unknown is associated with the
!    node, and a zero horizontal velocity is assumed.
!    If it is -1, then no unknown is associated with the
!    node, but a velocity is specified via the UBDRY routine
!    If it is positive, then the value is the index in the
!    F and G arrays of the unknown coefficient.
!    INDX(I,2) records information for vertical velocities.
!    If it is 0, then no unknown is associated with the
!    node, and a zero vertical velocity is assumed.
!    If it is positive, then the value is the index in the
!    F and G arrays of the unknown coefficient.
!
!    integer INSC(NP), records, for each node, whether an
!    unknown pressure is associated with the node.
!    If INSC(I) is zero, then no unknown is associated with the
!    node, and a pressure of 0 is assumed.
!    If INSC(I) is positive, then the value is the index in the
!    F and G arrays of the unknown coefficient.
!
!    integer IOUNIT, the FORTRAN unit number for the file
!    into which graphics data will be written.
!
!    integer IPIVOT(MAXEQN), pivot information used by the linear
!    solver.
!
!    integer ISOTRI(NELEMN), records whether or not a given
!    element is isometric or not.
!    0, element is not isometric.
!    1, element is isometric.
!
!    Input, integer IVUNIT, the output unit to which the data is to
!    be written.  The data file is unformated.
!
!    integer IWRITE, controls the amount of output produced by the
!    program.
!    0, minimal output.
!    1, normal output, plus graphics data files created by GDUMP, UVDUMP and XYDUMP.
!    2, copious output.
!
!    Input, integer IXUNIT, the output unit to which the data is to
!    be written.  The data file is unformated.
!
!    LOGICAL LONG,
!    .TRUE. if region is "long and thin", and
!    .FALSE. if region is "tall and skinny".
!    This determines how we number nodes, elements, and variables.
!
!    integer MAXEQN, the maximum number of equations allowed.
!
!    integer MAXNEW, the maximum number of Newton steps per iteration.
!
!    integer MAXROW, the maximum row real ( kind = 8 ) of the coefficient
!    matrix that is allowed.
!
!    integer MAXSEC, the maximum number of secant steps allowed.
!
!    integer MX, MX = 2*NX-1, the total number of grid points on
!    the horizontal side of the region.
!
!    integer MY, MY = 2*NY-1, the total number of grid points on
!    the vertical side of the region.
!
!    integer NELEMN, the number of elements used.
!
!    integer NEQN, the number of equations or functions for
!    the full system.
!
!    integer NLBAND, the number of diagonals below the main diagonal
!    of the matrix A which are nonzero.
!
!    integer NNODES, the number of nodes per element, 6.
!
!    integer NODE(NELEMN,6), records the global node numbers of the
!    6 nodes that make up each element.
!
!    integer NODEX0, the lowest numbered node in the column of
!    nodes where the profile is measured.
!
!    integer NP, the number of nodes.
!
!    integer NQUAD, the number of quadrature points, currently
!    set to 3.
!
!    integer NROW, the used row real ( kind = 8 ) of the coefficient
!    matrix.
!
!    integer NUMNEW, total number of Newton iterations taken in
!    the NSTOKE routine during the entire run.
!
!    integer NUMSEC, the number of secant steps taken.
!
!    integer NX, the number of "main" grid points on the horizontal
!    side of the region.
!
!    integer NY, the number of "main" grid points on the vertical
!    side of the region.
!
!    real PHI(NELEMN,NQUAD,NNODES,3).  Each entry of PHI contains
!    the value of a quadratic basis function or its derivative,
!    evaluated at a quadrature point.
!    In particular, PHI(I,J,K,1) is the value of the quadratic basis
!    function associated with local node K in element I, evaluatated
!    at quadrature point J.
!    PHI(I,J,K,2) is the X derivative of that same basis function,
!    PHI(I,J,K,3) is the Y derivative of that same basis function.
!
!    real PSI(NELEMN,NQUAD,NNODES).  Each entry of PSI contains
!    the value of a linear basis function evaluated at a
!    quadrature point.
!    PSI(I,J,K) is the value of the linear basis function associated
!    with local node K in element I, evaluated at quadrature point J.
!
!    real RES(MAXEQN), contains the residuals.
!
!    real REYNLD, the value of the Reynolds number.  In the
!    program's system of units, viscosity = 1 / REYNLD.
!
!    real RJPNEW, the derivative with respect to the
!    parameter A of the functional J.
!
!    real TARRAY(2), an array needed in order to store
!    results of a call to the UNIX CPU timing routine ETIME.
!
!    real TOLNEW, the convergence tolerance for the Newton
!    iteration in NSTOKE.
!
!    real TOLSEC, the convergence tolerance for the secant
!    iteration in the main program.
!
!    real XC(NP), XC(I) is the X coordinate of node I.
!
!    real XLNGTH, the length of the region.
!
!    real XM(NELEMN,NQUAD), XM(IT,I) is the X coordinate of the
!    I-th quadrature point in element IT.
!
!    real YC(NP), YC(I) is the Y coordinate of node I.
!
!    YLNGTH real YLNGTH, the height of the region.
!
!    real YM(NELEMN,NQUAD).  YM(IT,I) is the Y coordinate of the
!    I-th quadrature point in element IT.
!
  implicit none

  integer, parameter :: nx = 21
  integer, parameter :: ny = 7

  integer, parameter :: maxrow = 27*ny
  integer, parameter :: nelemn = 2*(nx-1)*(ny-1)
  integer, parameter :: mx = 2*nx-1
  integer, parameter :: my = 2*ny-1
  integer, parameter :: np = mx*my
  integer, parameter :: maxeqn = 2*mx*my+nx*ny
  integer, parameter :: nnodes = 6
  integer, parameter :: nquad = 3

  real ( kind = 8 ) a(maxrow,maxeqn)
  real ( kind = 8 ) a2
  real ( kind = 8 ) abound
  real ( kind = 8 ) anew
  real ( kind = 8 ) aold
  real ( kind = 8 ) area(nelemn)
  real ( kind = 8 ) dcda(my)
  real ( kind = 8 ) f(maxeqn)
  character fileg*30
  character fileu*30
  character filex*30
  real ( kind = 8 ) g(maxeqn)
  real ( kind = 8 ) gr(my,my)
  integer i
  integer iline(my)
  integer indx(np,2)
  integer insc(np)
  integer iounit
  integer ipivot(maxeqn)
  integer isotri(nelemn)
  integer iter
  integer ivunit
  integer iwrite
  integer ixunit
  integer j
  logical long
  integer maxnew
  integer maxsec
  integer nband
  integer neqn
  integer nlband
  integer node(nelemn,nnodes)
  integer nodex0
  integer npara
  integer nrow
  integer numnew
  integer numsec
  real ( kind = 8 ) para
  real ( kind = 8 ) phi(nelemn,nquad,nnodes,3)
  real ( kind = 8 ) psi(nelemn,nquad,nnodes)
  real ( kind = 8 ) r(my)
  real ( kind = 8 ) res(maxeqn)
  real ( kind = 8 ) reynld
  real ( kind = 8 ) rjpnew
  real ( kind = 8 ) rjpold
  real ( kind = 8 ) rtemp(my)
  real tarray(2)
  real ( kind = 8 ) temp
  real ( kind = 8 ) test
  real ( kind = 8 ) tolnew
  real ( kind = 8 ) tolsec
  real ( kind = 8 ) ui(my)
  real ( kind = 8 ) unew(my)
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) xlngth
  real ( kind = 8 ) xm(nelemn,nquad)
  real ( kind = 8 ) yc(np)
  real ( kind = 8 ) ylngth
  real ( kind = 8 ) ym(nelemn,nquad)

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CHANNEL'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Channel flow control problem'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Last modified:'
  write ( *, '(a)' ) '    4 September 1992.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Flow control problem:'
  write ( *, '(a)' ) '    Inflow controlled by one parameter.'
  write ( *, '(a)' ) '    Velocities measured along vertical line.'
  write ( *, '(a)' ) '    Try to match specified velocity profile.'
!
!  Set input data
!
  fileg = 'display.txt'
  fileu = 'uv.txt'
  filex = 'xy.txt'
  iounit = 2
  ivunit = 4
  iwrite = 2
  ixunit = 3
  maxnew = 10
  maxsec = 8
  npara = 1
  numnew = 0
  numsec = 0
  reynld = 1.0D+00
  rjpnew = 0.0D+00
  tolnew = 1.0D-04
  tolsec = 1.0D-06
  xlngth = 10.0D+00
  ylngth = 3.0D+00

  write (*,*) ' '
  write (*,*) 'NX = ', nx
  write (*,*) 'NY = ', ny
  write (*,*) 'Number of elements = ', nelemn
  write (*,*) 'Reynolds number = ', reynld
  write (*,*) 'Secant tolerance = ', tolsec
  write (*,*) 'Newton tolerance = ', tolnew
  write (*,*) ' '
!
!  SETGRD constructs grid, numbers unknowns, calculates areas,
!  and points for midpoint quadrature rule.
!
  call setgrd (indx, insc, isotri, iwrite, long, maxeqn, mx, my, &
    nelemn, neqn, nnodes, node, np, nx, ny )
!
!  Compute the bandwidth
!
  call setban(indx,insc,maxrow,nband,nelemn,nlband,nnodes, &
    node,np,nrow)
!
!  Record variable numbers along profile sampling line.
!
  call setlin(iline,indx,iwrite,long,mx,my,nodex0,np, &
    nx,ny,xlngth)
!
!  Set the coordinates of grid points.
!
  call setxy(iwrite,long,mx,my,np,nx,ny,xc,xlngth,yc,ylngth)
!
!  Set quadrature points
!
  call setqud(area,nelemn,nnodes,node,np,nquad,xc,xm,yc,ym)
!
!  Evaluate basis functions at quadrature points
!
  call setbas(nelemn,nnodes,node,np,nquad,phi,psi,xc,xm,yc,ym)
!
!  NSTOKE now solves the Navier Stokes problem for an inflow
!  parameter of 1.0.
!
  para = 1.0D+00
  write (*,*) ' '
  write (*,*) 'Solve Navier Stokes problem with parameter = ',para
  write (*,*) 'for profile at x = ', xc(nodex0)
  g(1:neqn) = 1.0D+00

  call nstoke (a,area,f,g,indx,insc,ipivot,iwrite, &
    maxnew,maxrow,nelemn,neqn,nlband,nnodes,node, &
    np,nquad,nrow,numnew,para,phi,psi,reynld,tolnew,yc)
!
!  RESID computes the residual at the given solution
!
  if ( 1 <= iwrite ) then
    call resid (area,f,indx,insc,iwrite,nelemn,neqn, &
      nnodes,node,np,nquad,para,phi,psi,res,reynld,yc)
  end if
!
!  GETG computes the internal velocity profile at X = XC(NODEX0), which will
!  be used to measure the goodness-of-fit of the later solutions.
!
  call getg ( f, iline, my, neqn, ui )

  if ( 1 <= iwrite ) then
    write (*,*) ' '
    write (*,*) 'U profile:'
    write (*,*) ' '
    write (*,'(5g14.6)') ui(1:my)
  end if
!
!  GRAM generates the Gram matrix GR and the vector
!  R = line integral of ui*phi
!
  call gram (gr,iline,indx,iwrite,my,nelemn,nnodes,node, &
    nodex0,np,para,r,ui,xc,yc)
!
!  GDUMP dumps information for graphics display by DISPLAY.
!
  if ( .false. ) then
    write (*,*) 'Writing graphics data to file '//fileg
    call delete(fileg)
    open (unit = iounit,file=fileg,form='formatted',status='new', &
      err = 50)
    rjpnew = 0.0D+00
    call gdump (f,indx,insc,iounit,isotri,long,nelemn,neqn, &
      nnodes,node,np,npara,nx,ny,para,reynld,rjpnew,xc,yc)
  end if
!
!  Write the XY data to a file.
!
  if ( .false. ) then
    call delete(filex)
    open(unit = ixunit,file=filex,form='formatted',status='new')
    call xy_plot3d (ixunit,long,np,nx,ny,xc,yc)
    close(unit = ixunit)
  else
    call delete ( filex )
    open ( unit = ixunit, file = filex, form = 'formatted', &
      status = 'new')
    call xy_table ( ixunit, np, xc, yc )
    close ( unit = ixunit )
  end if
!
!  Write the velocity data to a file.
!
  if ( .false. ) then
    call delete(fileu)
    open(unit = ivunit,file=fileu,form='formatted',status='new')
    call uv_plot3d (f,indx,insc,ivunit,long,mx,my, &
      nelemn,neqn,nnodes,node,np,para,a,reynld,yc)
    close(unit = ivunit)
  else
    call delete(fileu)
    open(unit = ivunit,file=fileu,form='formatted',status='new')
    call uv_table ( f, indx, ivunit, neqn, np, para, yc )
    close(unit = ivunit)
  end if
!
!  Destroy information about true solution
!
  f(1:neqn) = 0.0D+00
  g(1:neqn) = 0.0D+00
!
!  Secant iteration loop
!
  aold = 0.0D+00
  rjpold = 0.0D+00
  anew = 0.1D+00

  do iter = 1, maxsec

    numsec = numsec+1
    write (*,*) ' '
    write (*,*) 'Secant iteration ',iter
!
!  Solve for unew at new value of parameter anew
!
    write (*,*) ' '
    write (*,*) 'Solving Navier Stokes problem for parameter = ',anew
!
!  Use solution F at previous value of parameter for starting point.
!
    call dcopy(neqn,f,1,g,1)
    para = anew

    call nstoke (a,area,f,g,indx,insc,ipivot,iwrite, &
      maxnew,maxrow,nelemn,neqn,nlband,nnodes,node, &
      np,nquad,nrow,numnew,para,phi,psi,reynld,tolnew,yc)
!
!  Get velocity profile
!
    call getg ( f, iline, my, neqn, unew )

    if ( 1 <= iwrite ) then
      write (*,*) ' '
      write (*,*) 'Velocity profile:'
      write (*,*) ' '
      write (*,'(5g14.6)') unew(1:my)
    end if
!
!  Solve linear system for du/da
!
    para = anew
    abound = 1.0D+00
    call linsys (a,area,g,f,indx,insc,ipivot, &
      maxrow,nelemn,neqn,nlband,nnodes,node, &
      np,nquad,nrow,para,abound,phi,psi,reynld,yc)
!
!  Output in DCDA
!
    call getg ( g, iline, my, neqn, dcda )

    if ( 2 <= iwrite ) then
      write (*,*) ' '
      write (*,*) 'Sensitivities:'
      write (*,*) ' '
      write (*,'(5g14.6)') dcda(1:my)
    end if
!
!  Evaluate J prime at current value of parameter where J is
!  functional to be minimized.
!
!  JPRIME = 2.0 * DCDA(I) * (GR(I,J)*UNEW(J)-R(I))
!
    rjpnew = 0.0D+00
    do i = 1, my
      temp = -r(i)
      do j = 1, my
        temp = temp + gr(i,j) * unew(j)
      end do
      rjpnew = rjpnew + 2.0D+00 * dcda(i) * temp
    end do

    write (*,*) ' '
    write (*,*) 'Parameter  = ',anew,' J prime = ',rjpnew
!
!  Dump information for graphics
!
    if ( .false. ) then
      para = anew
      call gdump (f,indx,insc,iounit,isotri,long,nelemn,neqn, &
        nnodes,node,np,npara,nx,ny,para,reynld,rjpnew,xc,yc)
    end if
!
!  Update the estimate of the parameter using the secant step
!
    if (iter == 1) then
      a2 = 0.5D+00
    else
      a2 = aold-rjpold*(anew-aold)/(rjpnew-rjpold)
    end if

    aold = anew
    anew = a2
    rjpold = rjpnew
    test = abs(anew-aold)/abs(anew)

    write (*,*) 'New value of parameter = ',anew
    write(*,*)'Convergence test = ',test

    if (abs(anew-aold) < abs(anew)*tolsec) then
      write (*,*) 'Secant iteration converged.'
      go to 40
    end if

  end do

  write (*,*) 'Secant iteration failed to converge.'
   40 continue

  write (*,*)'Number of secant steps = ', numsec
  write (*,*)'Number of Newton steps = ', numnew
!
!  Close graphics file
!
  if ( .false. ) then
    close ( unit = iounit )
  end if
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CHANNEL:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
!
!  Error opening graphics file
!
   50 continue
  write (*,*) 'CHANNEL could not open the graphics file!'
  stop
end
function bsp ( xq, yq, it, iq, id, nelemn, nnodes, node, np, xc, yc )

!*****************************************************************************80
!
!! BSP evaluates the linear basis functions associated with pressure.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer nelemn
  integer nnodes
  integer np

  real ( kind = 8 ) bsp
  real ( kind = 8 ) d
  integer i1
  integer i2
  integer i3
  integer id
  integer iq
  integer iq1
  integer iq2
  integer iq3
  integer it
  integer node(nelemn,nnodes)
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) xq
  real ( kind = 8 ) yc(np)
  real ( kind = 8 ) yq

  iq1 = iq
  iq2 = mod(iq,3)+1
  iq3 = mod(iq+1,3)+1
  i1 = node(it,iq1)
  i2 = node(it,iq2)
  i3 = node(it,iq3)
  d = (xc(i2)-xc(i1))*(yc(i3)-yc(i1))-(xc(i3)-xc(i1))*(yc(i2)-yc(i1))

  if (id == 1) then

    bsp = 1.0+((yc(i2)-yc(i3))*(xq-xc(i1))+(xc(i3)-xc(i2))*(yq-yc(i1)))/d

  else if (id == 2) then

    bsp = (yc(i2)-yc(i3))/d

  else if (id == 3) then

    bsp = (xc(i3)-xc(i2))/d

  else

    write (*,*) 'BSP - fatal error!'
    write (*,*) 'unknown value of id = ',id
    stop

  end if

  return
end
subroutine daxpy ( n, da, dx, incx, dy, incy )

!*****************************************************************************80
!
!! DAXPY computes constant times a vector plus a vector.
!
!  Modified:
!
!    20 January 2007
!
!  Author:
!
!    Jack Dongarra
!
  real ( kind = 8 ) dx(*),dy(*),da
  integer i,incx,incy,ix,iy,m,n
!
  if ( n <= 0)return
  if (da  ==  0.0d+00 ) return
  if ( incx == 1.and.incy == 1)go to 20
!
!        code for unequal increments or equal increments
!          not equal to 1
!
  ix = 1
  iy = 1
  if ( incx < 0)ix = (-n+1)*incx + 1
  if ( incy < 0)iy = (-n+1)*incy + 1
  do i = 1,n
    dy(iy) = dy(iy) + da*dx(ix)
    ix = ix + incx
    iy = iy + incy
  end do
  return
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
   20 m = mod(n,4)
  if (  m  ==  0 ) go to 40
  do 30 i = 1,m
    dy(i) = dy(i) + da*dx(i)
   30 continue
  if (  n  <  4 ) return
   40 continue

  do i = m+1, n, 4
    dy(i) = dy(i) + da*dx(i)
    dy(i + 1) = dy(i + 1) + da*dx(i + 1)
    dy(i + 2) = dy(i + 2) + da*dx(i + 2)
    dy(i + 3) = dy(i + 3) + da*dx(i + 3)
  end do

  return
end
subroutine dcopy ( n, dx, incx, dy, incy )

!*****************************************************************************80
!
!! DCOPY copies a vector, x, to a vector, y.
!
!  Modified:
!
!    20 January 2007
!
!  Author:
!
!    Jack Dongarra
!
  real ( kind = 8 ) dx(*),dy(*)
  integer i,incx,incy,ix,iy,m,n

  if ( n <= 0)return
  if ( incx == 1.and.incy == 1)go to 20
!
!        code for unequal increments or equal increments
!          not equal to 1
!
  ix = 1
  iy = 1
  if ( incx < 0)ix = (-n+1)*incx + 1
  if ( incy < 0)iy = (-n+1)*incy + 1
  do i = 1,n
    dy(iy) = dx(ix)
    ix = ix + incx
    iy = iy + incy
  end do
  return
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
   20 m = mod(n,7)
  if (  m  ==  0 ) go to 40

  dy(1:m) = dx(1:m)

  if (  n  <  7 ) return
   40 continue
  do i = m+1, n ,7
    dy(i) = dx(i)
    dy(i + 1) = dx(i + 1)
    dy(i + 2) = dx(i + 2)
    dy(i + 3) = dx(i + 3)
    dy(i + 4) = dx(i + 4)
    dy(i + 5) = dx(i + 5)
    dy(i + 6) = dx(i + 6)
  end do

  return
end
function ddot ( n, dx, incx, dy, incy )

!*****************************************************************************80
!
!! DDOT forms the dot product of two vectors.
!
!  Modified:
!
!    20 January 2007
!
!  Author:
!
!    Jack Dongarra
!
  real ( kind = 8 ) ddot
  real ( kind = 8 ) dx(*),dy(*),dtemp
  integer i,incx,incy,ix,iy,m,n

  ddot = 0.0d+00
  dtemp = 0.0d+00
  if ( n <= 0)return
  if ( incx == 1.and.incy == 1)go to 20
!
!        code for unequal increments or equal increments
!          not equal to 1
!
  ix = 1
  iy = 1
  if ( incx < 0)ix = (-n+1)*incx + 1
  if ( incy < 0)iy = (-n+1)*incy + 1
  do i = 1,n
    dtemp = dtemp + dx(ix)*dy(iy)
    ix = ix + incx
    iy = iy + incy
  end do

  ddot = dtemp
  return
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
   20 m = mod(n,5)
  if (  m  ==  0 ) go to 40
  do i = 1,m
    dtemp = dtemp + dx(i)*dy(i)
  end do

  if (  n  <  5 ) go to 60
   40 continue
  do i = m+1, n, 5
    dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) + &
        dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
  end do
   60 ddot = dtemp
  return
end
subroutine delete (filnam)

!*****************************************************************************80
!
!! DELETE deletes a file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 January 2007
!
!  Author:
!
!    John Burkardt
!
  character filnam*(*)
  integer ios

  open (unit = 99,file=filnam,status='old', iostat = ios )

  if ( ios /= 0 ) then
    return
  end if

  close (unit = 99,status='delete', iostat = ios)

  return
end
subroutine dgbfa ( abd, lda, n, ml, mu, ipvt, info )

!*****************************************************************************80
!
!! DGBFA factors a band matrix by elimination.
!
!  Discussion:
!
!     dgbfa is usually called by dgbco, but it can be called
!     directly with a saving in time if  rcond  is not needed.
!
!  Modified:
!
!    20 January 2007
!
!  Author:
!
!    Cleve Moler
!
!  Parameters:
!
!     on entry
!
!        abd     real ( kind = 8 )(lda, n)
!                contains the matrix in band storage.  the columns
!                of the matrix are stored in the columns of  abd  and
!                the diagonals of the matrix are stored in rows
!                ml+1 through 2*ml+mu+1 of  abd .
!                see the comments below for details.
!
!        lda     integer
!                the leading dimension of the array  abd .
!                2*ml + mu + 1 <= LDA.
!
!        n       integer
!                the order of the original matrix.
!
!        ml      integer
!                number of diagonals below the main diagonal.
!                0 <= ml < n .
!
!        mu      integer
!                number of diagonals above the main diagonal.
!                0 <= mu < n .
!                more efficient if  ml <= mu .
!     on return
!
!        abd     an upper triangular matrix in band storage and
!                the multipliers which were used to obtain it.
!                the factorization can be written  a = l*u  where
!                l  is a product of permutation and unit lower
!                triangular matrices and  u  is upper triangular.
!
!        ipvt    integer(n)
!                an integer vector of pivot indices.
!
!        info    integer
!                = 0  normal value.
!                = k  if  u(k,k) == 0.0D+00 .  this is not an error
!                     condition for this subroutine, but it does
!                     indicate that dgbsl will divide by zero if
!                     called.  use  rcond  in dgbco for a reliable
!                     indication of singularity.
!
!     band storage
!
!           if  a  is a band matrix, the following program segment
!           will set up the input.
!
!                   ml = (band width below the diagonal)
!                   mu = (band width above the diagonal)
!                   m = ml + mu + 1
!                   do j = 1, n
!                      i1 = max ( 1, j-mu )
!                      i2 = min ( n, j+ml )
!                      do i = i1, i2
!                         k = i - j + m
!                         abd(k,j) = a(i,j)
!                      end do
!                   end do
!
!           this uses rows  ml+1  through  2*ml+mu+1  of  abd .
!           in addition, the first  ml  rows in  abd  are used for
!           elements generated during the triangularization.
!           the total number of rows needed in  abd  is  2*ml+mu+1 .
!           the  ml+mu by ml+mu  upper left triangle and the
!           ml by ml  lower right triangle are not referenced.
!
!     linpack. this version dated 08/14/78 .
!     cleve moler, university of new mexico, argonne national lab.
!
  integer lda
  integer n

  real ( kind = 8 ) abd(lda,n)
  integer i
  integer i0
  integer info
  integer ipvt(n)
  integer idamax
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
!  Zero next fill-in column.
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
!! DGBSL solves a banded system factored by DGBFA.
!
!  Discussion:
!
!    SGBSL can solve either a * x = b  or  trans(a) * x = b.
!
!  Parameters:
!
!     on entry
!
!        abd     real ( kind = 8 )(lda, n)
!                the output from dgbco or dgbfa.
!
!        lda     integer
!                the leading dimension of the array  abd .
!
!        n       integer
!                the order of the original matrix.
!
!        ml      integer
!                number of diagonals below the main diagonal.
!
!        mu      integer
!                number of diagonals above the main diagonal.
!
!        ipvt    integer(n)
!                the pivot vector from dgbco or dgbfa.
!
!        b       real ( kind = 8 )(n)
!                the right hand side vector.
!
!        job     integer
!                = 0         to solve  a*x = b ,
!                = nonzero   to solve  trans(a)*x = b , where
!                            trans(a)  is the transpose.
!
!     on return
!
!        b       the solution vector  x .
!
!     error condition
!
!        a division by zero will occur if the input factor contains a
!        zero on the diagonal.  technically this indicates singularity
!        but it is often caused by improper arguments or improper
!        setting of lda .  it will not occur if the subroutines are
!        called correctly and if dgbco has set 0.0 < RCOND
!        or dgbfa has set info == 0 .
!
!     to compute  inverse(a) * c  where  c  is a matrix
!     with  p  columns
!           call dgbco ( abd, lda, n, ml, mu, ipvt, rcond, z )
!           if (rcond is too small) go to ...
!           do j = 1, p
!              call dgbsl ( abd, lda, n, ml, mu, ipvt, c(1,j), 0 )
!           end do
!
!     linpack. this version dated 08/14/78 .
!     cleve moler, university of new mexico, argonne national lab.
!
  integer lda
  integer n

  real ( kind = 8 ) abd(lda,n)
  real ( kind = 8 ) b(n)
  integer ipvt(n)
  integer job
  integer k
  integer l
  integer la
  integer lb
  integer lm
  integer m
  integer ml
  integer mu
  real ( kind = 8 ) ddot
  real ( kind = 8 ) t

  m = mu + ml + 1
!
!  JOB = 0, Solve  a * x = b.
!
!  First solve l*y = b.
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
!  Now solve u*x = y.
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
!  JOB nonzero, solve  trans(a) * x = b.
!
!  First solve  trans(u)*y = b.
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
!  Now solve trans(l)*x = y
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
subroutine dscal ( n, da, dx, incx )

!*****************************************************************************80
!
!! DSCAL scales a vector by a constant.
!
!     uses unrolled loops for increment equal to one.
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx  <=  0.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
  real ( kind = 8 ) da,dx(*)
  integer i,incx,m,n,nincx
!
  if (  n <= 0 .or. incx <= 0 )return
  if ( incx == 1)go to 20
!
!        code for increment not equal to 1
!
  nincx = n*incx
  do i = 1,nincx,incx
    dx(i) = da*dx(i)
  end do
  return
!
!        code for increment equal to 1
!
!
!        clean-up loop
!
   20 continue

  m = mod(n,5)
  if (  m  ==  0 ) go to 40
  dx(1:m) = da*dx(1:m)
  if (  n  <  5 ) return

40 continue

  do i = m+1,n,5
    dx(i) = da*dx(i)
    dx(i + 1) = da*dx(i + 1)
    dx(i + 2) = da*dx(i + 2)
    dx(i + 3) = da*dx(i + 3)
    dx(i + 4) = da*dx(i + 4)
  end do

  return
end
subroutine gdump (f,indx,insc,iounit,isotri,long,nelemn,neqn, &
  nnodes,node,np,npara,nx,ny,para,reynld,rjpnew,xc,yc)

!*****************************************************************************80
!
!! GDUMP writes information to a file.
!
!  Discussion:
!
!    The information can be used to create
!    graphics images.  In order to keep things simple, exactly one
!    value, real or integer, is written per record.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NPARA, the number of parameters.  Fixed at 1
!    for now.
!
!    Input, real ( kind = 8 ) PARA(MAXPAR), the parameters.
!
  integer nelemn
  integer neqn
  integer nnodes
  integer np

  real ( kind = 8 ) f(neqn)
  real ( kind = 8 ) fval
  integer i
  integer indx(np,2)
  integer insc(np)
  integer iounit
  integer, save :: iset = 0
  integer isotri(nelemn)
  integer j
  logical long
  integer node(nelemn,nnodes)
  integer npara
  integer nx
  integer ny
  real ( kind = 8 ) para
  real ( kind = 8 ) reynld
  real ( kind = 8 ) rjpnew
  real ( kind = 8 ) ubdry
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) yc(np)

  iset = iset+1

  write(iounit,*)long
  write (iounit,*) nelemn
  write (iounit,*) np
  write (iounit,*) npara
  write (iounit,*) nx
  write (iounit,*) ny
!
!  Pressures
!
  do i = 1, np
    j = insc(i)
    if (j <= 0) then
      fval = 0.0D+00
    else
      fval = f(j)
    end if
    write (iounit,*) fval
  end do
!
!  Horizontal velocities, U
!
  do i = 1, np
    j = indx(i,1)
    if (j == 0) then
      fval = 0.0D+00
    else if (j < 0) then
      fval = ubdry(yc(i),para)
    else
      fval = f(j)
    end if
    write (iounit,*) fval
  end do
!
!  Vertical velocities, V
!
  do i = 1, np
    j = indx(i,2)
    if (j <= 0) then
      fval = 0.0D+00
    else
      fval = f(j)
    end if
    write (iounit,*) fval
  end do

  do i = 1, np
    write (iounit,*) indx(i,1)
    write (iounit,*) indx(i,2)
  end do

  do i = 1, np
    write (iounit,*) insc(i)
  end do

  do i = 1, nelemn
    write (iounit,*) isotri(i)
  end do

  do i = 1, nelemn
    do j = 1, 6
      write (iounit,*) node(i,j)
    end do
  end do

  write (iounit,*) para
  write (iounit,*) reynld
  write (iounit,*) rjpnew

  do i = 1, np
    write (iounit,*) xc(i)
  end do

  do i = 1, np
    write (iounit,*) yc(i)
  end do

  write (*,*) 'GDUMP wrote data set ',iset,' to file.'

  return
end
subroutine getg ( f, iline, my, neqn, u )

!*****************************************************************************80
!
!! GETG outputs field values along the profile line X = XZERO.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer neqn
  integer my

  real ( kind = 8 ) f(neqn)
  integer iline(my)
  integer j
  integer k
  real ( kind = 8 ) u(my)

  do j = 1, my
    k = iline(j)
    if ( 0 < k ) then
      u(j) = f(k)
    else
      u(j) = 0.0D+00
    end if
  end do

  return
end
subroutine gram ( gr, iline, indx, iwrite, my, nelemn, nnodes, node, &
  nodex0, np, para, r, ui, xc, yc )

!*****************************************************************************80
!
!! GRAM computes the Gram matrix, GR(I,J) = INTEGRAL PHI(I)*PHI(J).
!
!  and the vector R(I) = INTEGRAL UI*PHI(I).
!
!  The integrals are computed along the line where the profile is
!  specified.  The three point Gauss quadrature rule is used for the
!  line integral.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer my
  integer nelemn
  integer nnodes
  integer np

  real ( kind = 8 ) ar
  real ( kind = 8 ) bb
  real ( kind = 8 ) bbb
  real ( kind = 8 ) bbx
  real ( kind = 8 ) bby
  real ( kind = 8 ) bma2
  real ( kind = 8 ) bx
  real ( kind = 8 ) by
  real ( kind = 8 ) gr(my,my)
  integer i
  integer igetl
  integer ii
  integer iline(my)
  integer indx(np,2)
  integer ip
  integer ipp
  integer iq
  integer iqq
  integer iquad
  integer it
  integer iun
  integer iwrite
  integer j
  integer jj
  integer k
  integer kk
  integer node(nelemn,nnodes)
  integer nodex0
  real ( kind = 8 ) para
  real ( kind = 8 ) r(my)
  real ( kind = 8 ) ubc
  real ( kind = 8 ) ubdry
  real ( kind = 8 ) ui(my)
  real ( kind = 8 ) uiqdpt
  real ( kind = 8 ) wt(3)
  real ( kind = 8 ) x
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) xzero
  real ( kind = 8 ) y
  real ( kind = 8 ) yq(3)
  real ( kind = 8 ) yc(np)
!
!  Input values for 3 point Gauss quadrature
!
  wt(1) = 5.0D+00 / 9.0D+00
  wt(2) = 8.0D+00 / 9.0D+00
  wt(3) = wt(1)
  yq(1) = -0.7745966692D+00
  yq(2) = 0.0D+00
  yq(3) = -yq(1)
!
!  zero arrays
!
  r(1:my) = 0.0D+00
  gr(1:my,1:my) = 0.0D+00
!
!  Compute line integral by looping over intervals along line
!  using three point Gauss quadrature
!
  xzero = xc(nodex0)
  do 70 it = 1, nelemn
!
!  Check to see if we are in a triangle with a side along line
!  x = xzero.  If not, skip out
!
    k = node(it,1)
    kk = node(it,2)

    if ( 1.0D-04 < abs(xc(k)-xzero) ) then
      cycle
    end if

    if ( 1.0D-04 < abs(xc(kk)-xzero) ) then
      cycle
    end if

    do 60 iquad = 1, 3
      bma2 = (yc(kk)-yc(k))/2.0
      ar = bma2*wt(iquad)
      x = xzero
      y = yc(k)+bma2*(yq(iquad)+1.0D+00 )
!
!  Compute u internal at quadrature points
!
      uiqdpt = 0
      do 30 iq = 1, nnodes
        if ( 4 < iq ) go to 30
        if (iq == 3) go to 30
        call qbf (x,y,it,iq,bb,bx,by,nelemn,nnodes,node,np,xc,yc)
        ip = node(it,iq)
        iun = indx(ip,1)
        if ( 0 < iun ) then
          ii = igetl(iun,iline,my)
          uiqdpt = uiqdpt+bb*ui(ii)
        else if (iun < 0) then
          ubc = ubdry(yc(ip),para)
          uiqdpt = uiqdpt+bb*ubc
        end if
   30     continue
!
!  Only loop over nodes lying on line x = xzero
!
      do 50 iq = 1, nnodes
        if ( iq == 1.or.iq == 2.or.iq == 4 ) then
          ip = node(it,iq)
          call qbf (x,y,it,iq,bb,bx,by,nelemn,nnodes,node,np,xc,yc)
          i = indx(ip,1)
          if (i <= 0) go to 50
          ii = igetl(i,iline,my)
          r(ii) = r(ii)+bb*uiqdpt*ar

          do iqq = 1, nnodes
            if ( iqq == 1.or.iqq == 2.or.iqq == 4 ) then
              ipp = node(it,iqq)
              call qbf (x,y,it,iqq,bbb,bbx,bby,nelemn,nnodes, &
                node,np,xc,yc)
              j = indx(ipp,1)
              if (j /= 0) then
                jj = igetl(j,iline,my)
                gr(ii,jj) = gr(ii,jj)+bb*bbb*ar
              end if
            end if
          end do

        end if
   50     continue
   60   continue
   70 continue

  if ( 2 <= iwrite ) then
    write(*,*)' '
    write(*,*)'Gram matrix:'
    write(*,*)' '
    do i = 1,my
      do j = 1,my
        write(*,*)i,j,gr(i,j)
      end do
    end do
    write(*,*)' '
    write(*,*)'R vector:'
    write(*,*)' '
    do i = 1,my
      write(*,*)i,r(i)
    end do
  end if

  return
end
function idamax ( n, dx, incx )

!*****************************************************************************80
!
!! IDAMAX finds the index of element having max. absolute value.
!
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx  <=  0.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
  integer idamax
  real ( kind = 8 ) dmax
  real ( kind = 8 ) dx(*)
  integer i,incx,ix,n

  idamax = 0
  if (  n < 1 .or. incx <= 0 ) return
  idamax = 1
  if ( n == 1)return
  if ( incx == 1)go to 20
!
!        code for increment not equal to 1
!
  ix = 1
  dmax = dabs(dx(1))
  ix = ix + incx
  do i = 2,n
     if ( dabs(dx(ix)) <= dmax) go to 5
     idamax = i
     dmax = dabs(dx(ix))
    5    ix = ix + incx
  end do

  return
!
!        code for increment equal to 1
!
   20 continue

  dmax = abs ( dx(1) )

  do i = 2,n
    if ( dmax < abs ( dx(i) ) ) then
      idamax = i
      dmax = abs ( dx(i) )
    end if
  end do

  return
end
function igetl ( k, iline, my )

!*****************************************************************************80
!
!! IGETL gets the local unknown number along the profile line.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer my

  integer igetl
  integer iline(my)
  integer j
  integer k

  do j = 1, my
    if (iline(j) == k) then
      igetl = j
      return
    end if
  end do

  write ( *, * ) ' '
  write (*,*) 'IGETL - fatal error!'
  write (*,*) '  Unable to get local unknown number for '
  write (*,*) '  Global variable number ',k
  igetl = 0
  stop
end
subroutine linsys ( a, area, f, g, indx, insc, ipivot, maxrow, nelemn, &
  neqn, nlband, nnodes, node, np, nquad, nrow, para1, para2, phi, psi, &
  reynld, yc )

!*****************************************************************************80
!
!! LINSYS sets up and solves the linear system.
!
!  Discussion:
!
!    The G array contains the previous solution.
!
!    The F array contains the right hand side initially and then the
!    current solution.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 February 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer maxrow
  integer nelemn
  integer neqn
  integer nnodes
  integer np
  integer nquad
  integer nrow

  real ( kind = 8 ) a(maxrow,neqn)
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
  real ( kind = 8 ) f(neqn)
  real ( kind = 8 ) g(neqn)
  integer i
  integer idamax
  integer ihor
  integer indx(np,2)
  integer info
  integer insc(np)
  integer ioff
  integer ip
  integer ipivot(neqn)
  integer ipp
  integer iprs
  integer iq
  integer iqq
  integer iquad
  integer it
  integer iuse
  integer iver
  integer j
  integer job
  integer jp
  integer ju
  integer jv
  integer nlband
  integer node(nelemn,nnodes)
  real ( kind = 8 ) para1
  real ( kind = 8 ) para2
  real ( kind = 8 ) phi(nelemn,nquad,nnodes,3)
  real ( kind = 8 ) psi(nelemn,nquad,nnodes)
  real ( kind = 8 ) reynld
  real ( kind = 8 ) ubdry
  real ( kind = 8 ) un(2)
  real ( kind = 8 ) unx(2)
  real ( kind = 8 ) uny(2)
  real ( kind = 8 ) uu
  real ( kind = 8 ) visc
  real ( kind = 8 ) yc(np)

  ioff = nlband + nlband + 1
  visc = 1.0D+00 / reynld
  f(1:neqn) = 0.0D+00
  a(1:nrow,1:neqn) = 0.0D+00
!
!  For each element,
!
  do it = 1, nelemn

    ar = area(it) / 3.0D+00
!
!  and for each quadrature point in the element,
!
    do iquad = 1, nquad
!
!  Evaluate velocities at quadrature point
!
      call uval (g,indx,iquad,it,nelemn,neqn,nnodes,node, &
        np,nquad,para1,phi,un,unx,uny,yc)
!
!  For each basis function,
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
          f(ihor) = f(ihor)+ar*bb*(un(1)*unx(1)+un(2)*uny(1))
        end if

        if ( 0 < iver ) then
          f(iver) = f(iver)+ar*bb*(un(1)*unx(2)+un(2)*uny(2))
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
                + bb*(bbb*unx(1)+bbx*un(1)+bby*un(2)))
            end if

            if ( 0 < iver ) then
              iuse = iver-ju+ioff
              a(iuse,ju) = a(iuse,ju)+ar*bb*bbb*unx(2)
            end if

            if ( 0 < iprs ) then
              iuse = iprs-ju+ioff
              a(iuse,ju) = a(iuse,ju)+ar*bbx*bbl
            end if

          else if ( ju < 0 ) then

            uu = ubdry(yc(ipp),para2)
            if ( 0 < ihor ) then
              f(ihor) = f(ihor)-ar*uu*(visc*(by*bby+bx*bbx) &
                + bb*(bbb*unx(1)+bbx*un(1)+bby*un(2)))
            end if

            if ( 0 < iver ) then
              f(iver) = f(iver)-ar*uu*bb*bbb*unx(2)
            end if

            if ( 0 < iprs ) then
              f(iprs) = f(iprs)-ar*uu*bbx*bbl
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

          end if
!
!  Pressure variable
!
          if ( 0 < jp ) then

            if ( 0 < ihor ) then
              iuse = ihor-jp+ioff
              a(iuse,jp) = a(iuse,jp)-ar*bx*bbbl
            end if

            if ( 0 < iver ) then
              iuse = iver-jp+ioff
              a(iuse,jp) = a(iuse,jp)-ar*by*bbbl
            end if

          end if

        end do
      end do
    end do
  end do
!
!  To avoid singularity of the pressure system, the last pressure
!  is simply assigned a value of 0.
!
  f(neqn) = 0.0D+00
  do j = neqn-nlband, neqn-1
    i = neqn-j+ioff
    a(i,j) = 0.0D+00
  end do
  a(ioff,neqn) = 1.0D+00
!
!  Factor the matrix
!
  call dgbfa ( a, maxrow, neqn, nlband, nlband, ipivot, info )

  if ( info /= 0 ) then
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
subroutine nstoke (a,area,f,g,indx,insc,ipivot,iwrite,maxnew,maxrow, &
  nelemn,neqn,nlband,nnodes,node,np,nquad,nrow,numnew,para,phi,psi, &
  reynld,tolnew,yc)

!*****************************************************************************80
!
!! NSTOKE solves the Navier Stokes equation using Taylor-Hood elements.
!
!  The G array contains the previous iterate.
!
!  The F array contains the right hand side initially and then the
!  current iterate.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer maxrow
  integer nelemn
  integer neqn
  integer nnodes
  integer np
  integer nquad
  integer nrow

  real ( kind = 8 ) a(maxrow,neqn)
  real ( kind = 8 ) area(nelemn)
  double precision diff
  real ( kind = 8 ) f(neqn)
  real ( kind = 8 ) g(neqn)
  integer idamax
  integer indx(np,2)
  integer insc(np)
  integer ipivot(neqn)
  integer iter
  integer iwrite
  integer maxnew
  integer nlband
  integer node(nelemn,nnodes)
  integer numnew
  real ( kind = 8 ) para
  real ( kind = 8 ) phi(nelemn,nquad,nnodes,3)
  real ( kind = 8 ) psi(nelemn,nquad,nnodes)
  real ( kind = 8 ) reynld
  real ( kind = 8 ) tolnew
  real ( kind = 8 ) yc(np)

  do iter = 1, maxnew

    numnew = numnew + 1

    call linsys (a,area,f,g,indx,insc,ipivot, &
      maxrow,nelemn,neqn,nlband,nnodes,node, &
      np,nquad,nrow,para,para,phi,psi,reynld,yc)
!
!  Check for convergence
!
    g(1:neqn) = g(1:neqn) - f(1:neqn)
    diff = abs ( g(idamax(neqn,g,1)) )

    if ( 1 <= iwrite ) then
      write(*,*)'NSTOKE iteration ',iter,' Mnorm = ',diff
    end if

    g(1:neqn) = f(1:neqn)

    if ( diff <= tolnew ) then
      write (*,*) 'Navier Stokes iteration converged in ', &
        iter,' iterations.'
      return
    end if

  end do

  write (*,*) 'Navier Stokes solution did not converge!'

  return
end
subroutine pval (g,insc,long,mx,my,nelemn,neqn,nnodes,node,np,press)

!*****************************************************************************80
!
!! PVAL computes a table of pressures.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer nelemn
  integer neqn
  integer nnodes
  integer np

  real ( kind = 8 ) g(neqn)
  integer i
  integer insc(np)
  integer ip
  integer iq
  integer it
  integer ivar
  integer j
  logical long
  integer mx
  integer my
  integer node(nelemn,nnodes)
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
      else
        press(i,j) = 0.0D+00
      end if

    end do
  end do
!
!  Interpolate the pressures at points (even, odd) and (odd, even).
!
  do i = 2,mx-1,2
    do j = 1,my,2
      press(i,j) = 0.5D+00*(press(i-1,j)+press(i+1,j))
    end do
  end do

  do j = 2,my-1,2
    do i = 1,mx,2
      press(i,j) = 0.5D+00*(press(i,j-1)+press(i,j+1))
    end do
  end do
!
!  Interpolate the pressures at points (even,even).
!
  do j = 2,my-1,2
    do i = 2,mx-1,2
      press(i,j) = 0.5D+00*(press(i-1,j-1)+press(i+1,j+1))
    end do
  end do

  return
end
subroutine qbf (x,y,it,in,bb,bx,by,nelemn,nnodes,node,np,xc,yc)

!*****************************************************************************80
!
!! QBF evaluates a quadratic basis function in a triangle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer nelemn
  integer nnodes
  integer np

  real ( kind = 8 ) bb
  real ( kind = 8 ) bx
  real ( kind = 8 ) by
  real ( kind = 8 ) c
  real ( kind = 8 ) d
  integer i1
  integer i2
  integer i3
  integer in
  integer in1
  integer in2
  integer in3
  integer inn
  integer it
  integer j1
  integer j2
  integer j3
  integer node(nelemn,nnodes)
  real ( kind = 8 ) s
  real ( kind = 8 ) t
  real ( kind = 8 ) x
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) y
  real ( kind = 8 ) yc(np)

  if (in <= 3) then
    in1 = in
    in2 = mod(in,3)+1
    in3 = mod(in+1,3)+1
    i1 = node(it,in1)
    i2 = node(it,in2)
    i3 = node(it,in3)
    d = (xc(i2)-xc(i1))*(yc(i3)-yc(i1))-(xc(i3)-xc(i1))*(yc(i2)-yc(i1))
    t = 1.0+((yc(i2)-yc(i3))*(x-xc(i1))+(xc(i3)-xc(i2))*(y-yc(i1)))/d
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
    t = 1.0D+00+((yc(i2)-yc(i3))*(x-xc(i1))+(xc(i3)-xc(i2))*(y-yc(i1)))/d
    s = 1.0D+00+((yc(j2)-yc(j3))*(x-xc(j1))+(xc(j3)-xc(j2))*(y-yc(j1)))/c
    bb = 4.0D+00*s*t
    bx = 4.0D+00*(t*(yc(j2)-yc(j3))/c+s*(yc(i2)-yc(i3))/d)
    by = 4.0D+00*(t*(xc(j3)-xc(j2))/c+s*(xc(i3)-xc(i2))/d)
  end if

  return
end
subroutine resid (area,g,indx,insc,iwrite,nelemn,neqn,nnodes, &
  node,np,nquad,para,phi,psi,res,reynld,yc)

!*****************************************************************************80
!
!! RESID computes the residual.
!
!  Discussion:
!
!    The G array contains the current iterate.
!
!    The RES array will contain the value of the residual.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer nelemn
  integer neqn
  integer nnodes
  integer np
  integer nquad

  real ( kind = 8 ) aijpu
  real ( kind = 8 ) aijpv
  real ( kind = 8 ) aijup
  real ( kind = 8 ) aijuu
  real ( kind = 8 ) aijuv
  real ( kind = 8 ) aijvp
  real ( kind = 8 ) aijvu
  real ( kind = 8 ) aijvv
  real ( kind = 8 ) ar
  real ( kind = 8 ) area(nelemn)
  real ( kind = 8 ) bb
  real ( kind = 8 ) bbb
  real ( kind = 8 ) bbbl
  real ( kind = 8 ) bbbx
  real ( kind = 8 ) bbby
  real ( kind = 8 ) bbl
  real ( kind = 8 ) bbx
  real ( kind = 8 ) bby
  real ( kind = 8 ) bx
  real ( kind = 8 ) by
  real ( kind = 8 ) g(neqn)
  integer i
  integer ibad
  integer ihor
  integer imax
  integer indx(np,2)
  integer insc(np)
  integer ip
  integer ipp
  integer iprs
  integer iq
  integer iqq
  integer iquad
  integer it
  integer iver
  integer iwrite
  integer j
  integer jp
  integer ju
  integer jv
  integer node(nelemn,nnodes)
  real ( kind = 8 ) para
  real ( kind = 8 ) phi(nelemn,nquad,nnodes,3)
  real ( kind = 8 ) psi(nelemn,nquad,nnodes)
  real ( kind = 8 ) res(neqn)
  real ( kind = 8 ) reynld
  real ( kind = 8 ) rmax
  real ( kind = 8 ) test
  real ( kind = 8 ) ubdry
  real ( kind = 8 ) un(2)
  real ( kind = 8 ) unx(2)
  real ( kind = 8 ) uny(2)
  real ( kind = 8 ) uu
  real ( kind = 8 ) visc
  real ( kind = 8 ) yc(np)

  visc = 1.0D+00 / reynld

  res(1:neqn) = 0.0D+00
!
!  For each element,
!
  do 90 it = 1, nelemn

    ar = area(it) / 3.0D+00
!
!  and for each quadrature point in the element,
!
    do 80 iquad = 1, nquad
!
!  Evaluate velocities at quadrature point
!
      call uval (g,indx,iquad,it,nelemn,neqn,nnodes,node, &
        np,nquad,para,phi,un,unx,uny,yc)
!
!  For each basis function,
!
      do 70 iq = 1, nnodes
        ip = node(it,iq)
        bb = phi(it,iquad,iq,1)
        bx = phi(it,iquad,iq,2)
        by = phi(it,iquad,iq,3)
        bbl = psi(it,iquad,iq)
        iprs = insc(ip)
        ihor = indx(ip,1)
        iver = indx(ip,2)

        if ( 0 < ihor ) then
          res(ihor) = res(ihor)+(un(1)*unx(1)+un(2)*uny(1))*bb*ar
        end if

        if ( 0 < iver ) then
          res(iver) = res(iver)+(un(1)*unx(2)+un(2)*uny(2))*bb*ar
        end if
!
!  For another basis function,
!
        do 50 iqq = 1, nnodes
          ipp = node(it,iqq)
          bbb = phi(it,iquad,iqq,1)
          bbx = phi(it,iquad,iqq,2)
          bby = phi(it,iquad,iqq,3)
          bbbl = psi(it,iquad,iqq)
          ju = indx(ipp,1)
          jv = indx(ipp,2)
          jp = insc(ipp)

          if ( 0 < ju ) then
            if ( 0 < ihor ) then
              aijuu = visc*(by*bby+bx*bbx) &
                + bb*(bbb*unx(1)+bbx*un(1)+bby*un(2))
              res(ihor) = res(ihor)+aijuu*ar*g(ju)
            end if
            if ( 0 < iver ) then
              aijvu = bb*bbb*unx(2)
              res(iver) = res(iver)+aijvu*ar*g(ju)
            end if
            if ( 0 < iprs ) then
              aijpu = bbx*bbl
              res(iprs) = res(iprs)+aijpu*ar*g(ju)
            end if
          else if ( ju < 0 ) then
            uu = ubdry(yc(ipp),para)
            if ( 0 < ihor ) then
              aijuu = visc*(by*bby+bx*bbx) &
                + bb*(bbb*unx(1)+bbx*un(1)+bby*un(2))
              res(ihor) = res(ihor)+ar*aijuu*uu
            end if
            if ( 0 < iver ) then
              aijvu = bb*bbb*unx(2)
              res(iver) = res(iver)+ar*aijvu*uu
            end if
            if ( 0 < iprs ) then
              aijpu = bbx*bbl
              res(iprs) = res(iprs)+ar*aijpu*uu
            end if
          end if

          if ( 0 < jv ) then
            if ( 0 < ihor ) then
              aijuv = bb*bbb*uny(1)
              res(ihor) = res(ihor)+aijuv*ar*g(jv)
            end if
            if ( 0 < iver ) then
              aijvv = visc*(by*bby+bx*bbx) &
                +bb*(bbb*uny(2)+bby*un(2)+bbx*un(1))
              res(iver) = res(iver)+aijvv*ar*g(jv)
            end if
            if ( 0  < iprs ) then
              aijpv = bby*bbl
              res(iprs) = res(iprs)+aijpv*ar*g(jv)
            end if
          end if

          if ( 0 < jp ) then
            if ( 0 < ihor ) then
              aijup = -bx*bbbl
              res(ihor) = res(ihor)+aijup*ar*g(jp)
            end if
            if ( 0 < iver ) then
              aijvp = -by*bbbl
              res(iver) = res(iver)+aijvp*ar*g(jp)
            end if
          end if

   50       continue
   70     continue
   80   continue
   90 continue

  res(neqn) = g(neqn)

  rmax = 0.0D+00
  imax = 0
  ibad = 0

  do i = 1,neqn

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

    do j = 1,np

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
subroutine setban(indx,insc,maxrow,nband,nelemn,nlband,nnodes, &
  node,np,nrow)

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
!    20 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer nelemn
  integer nnodes
  integer np

  integer i
  integer indx(np,2)
  integer insc(np)
  integer ip
  integer ipp
  integer iq
  integer iqq
  integer it
  integer iuk
  integer iukk
  integer j
  integer maxrow
  integer nband
  integer nlband
  integer node(nelemn,nnodes)
  integer nrow

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
              nlband = max(nlband,j-i)
            end do
          end do
        end if

      end do

    end do
  end do

  nband = nlband+nlband+1
  nrow = nlband+nlband+nlband+1

  write (*,*) 'Lower bandwidth = ',nlband
  write (*,*) 'Total bandwidth = ',nband
  write (*,*) 'NROW  = ',nrow
  if ( maxrow < nrow ) then
    write(*,*)'SETBAN - NROW is too large!'
    write(*,*)'The maximum allowed is ',maxrow
    stop
  end if

  return
end
subroutine setbas(nelemn,nnodes,node,np,nquad,phi,psi,xc,xm,yc,ym)

!*****************************************************************************80
!
!! SETBAS computes the basis functions at each integration point.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer nelemn
  integer nnodes
  integer np
  integer nquad

  real ( kind = 8 ) bb
  real ( kind = 8 ) bsp
  real ( kind = 8 ) bx
  real ( kind = 8 ) by
  integer iq
  integer it
  integer j
  integer node(nelemn,nnodes)
  real ( kind = 8 ) phi(nelemn,nquad,nnodes,3)
  real ( kind = 8 ) psi(nelemn,nquad,nnodes)
  real ( kind = 8 ) x
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) xm(nelemn,nquad)
  real ( kind = 8 ) y
  real ( kind = 8 ) yc(np)
  real ( kind = 8 ) ym(nelemn,nquad)

  do it = 1,nelemn
    do j = 1,nquad
      x = xm(it,j)
      y = ym(it,j)
      do iq = 1,6
        psi(it,j,iq) = bsp(x,y,it,iq,1,nelemn,nnodes,node,np,xc,yc)
        call qbf(x,y,it,iq,bb,bx,by,nelemn,nnodes,node,np,xc,yc)
        phi(it,j,iq,1) = bb
        phi(it,j,iq,2) = bx
        phi(it,j,iq,3) = by
      end do
    end do
  end do
  return
end
subroutine setgrd (indx,insc,isotri,iwrite,long,maxeqn,mx,my, &
  nelemn,neqn,nnodes,node,np,nx,ny)

!*****************************************************************************80
!
!! SETGRD sets up the grid for the problem..
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer my
  integer nelemn
  integer nnodes
  integer np

  integer i
  integer ic
  integer icnt
  integer ielemn
  integer indx(np,2)
  integer insc(np)
  integer ip
  integer ip1
  integer ip2
  integer isotri(nelemn)
  integer it
  integer iwrite
  integer jc
  integer jcnt
  logical long
  integer maxeqn
  integer mx
  integer neqn
  integer node(nelemn,nnodes)
  integer nx
  integer ny
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
!  Set parameters for Taylor Hood element
!
  write (*,*) ' '
  write (*,*) 'SETGRD: Taylor Hood element'
!
!  Construct grid coordinates, elements, and ordering of unknowns
!
  neqn = 0
  ielemn = 0

  do ip = 1,np

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
    if ( (icnt == 1.and.jcnt.eq.1).and. &
       (ic /= mx).and.(jc /= my) ) then

      if ( long ) then
        ip1 = ip+my
        ip2 = ip+my+my
        ielemn = ielemn+1
        node(ielemn,1) = ip
        node(ielemn,2) = ip2+2
        node(ielemn,3) = ip2
        node(ielemn,4) = ip1+1
        node(ielemn,5) = ip2+1
        node(ielemn,6) = ip1
        isotri(ielemn) = 0
        ielemn = ielemn+1
        node(ielemn,1) = ip
        node(ielemn,2) = ip+2
        node(ielemn,3) = ip2+2
        node(ielemn,4) = ip+1
        node(ielemn,5) = ip1+2
        node(ielemn,6) = ip1+1
        isotri(ielemn) = 0

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
        isotri(ielemn) = 0
        ielemn = ielemn+1
        node(ielemn,1) = ip
        node(ielemn,2) = ip2+2
        node(ielemn,3) = ip+2
        node(ielemn,4) = ip1+1
        node(ielemn,5) = ip1+2
        node(ielemn,6) = ip+1
        isotri(ielemn) = 0
      end if
    end if
!
!  Consider whether velocity unknowns should be associated with this node.
!
!  If we are in column 1, horizontal velocities are specified, and
!  vertical velocities are zero.
!
    if ( ic == 1.and.1 < jc.and.jc < my ) then
      indx(ip,1) = -1
      indx(ip,2) = 0
!
!  If we are in column MX, horizontal velocities are unknown, and
!  vertical velocities are zero.
!
    else if ( ic == mx.and.1 < jc.and.jc < my ) then
      neqn = neqn+1
      indx(ip,1) = neqn
      indx(ip,2) = 0
!
!  Otherwise, if we are in row 1 or row MY, both horizontal and
!  vertical velocities are zero.
!
    else if ( jc == 1.or.jc == my ) then
      indx(ip,1) = 0
      indx(ip,2) = 0
!
!  Otherwise, we are at an interior node
!
    else
      neqn = neqn+2
      indx(ip,1) = neqn-1
      indx(ip,2) = neqn
    end if
!
!  Consider whether a pressure unknown should be associated with this node.
!  The answer is yes if both nodes are odd.
!
    if (jcnt == 1.and.icnt == 1) then
      neqn = neqn+1
      insc(ip) = neqn
    else
      insc(ip) = 0
    end if

  end do
!
!  If debugging is requested, print out data.
!
  if ( 2 <= iwrite ) then
    write(*,*)' '
    write(*,*)'    I      INDX 1 & 2, INSC'
    write(*,*)' '
    do i = 1,np
      write (*,'(2xi6,2x,i6,2x,i6,2x,i6)') i,indx(i,1:2),insc(i)
    end do
    write(*,*)' '
    write(*,*)'    IT    NODE(IT,1:6)'
    write(*,*)' '
    do it = 1,nelemn
      write (*,'(7i6)') it, node(it,1:6)
    end do
  end if

  write (*,*) 'Number of unknowns = ',neqn
  if ( maxeqn < neqn ) then
    write(*,*)'SETGRD - Too many unknowns!'
    write(*,*)'The maximum allowed is MAXEQN = ',maxeqn
    write(*,*)'This problem requires NEQN = ',neqn
    stop
  end if

  return
end
subroutine setlin(iline,indx,iwrite,long,mx,my,nodex0,np, &
  nx,ny,xlngth)

!*****************************************************************************80
!
!! SETLIN gets the unknown indices along the profile line.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer my
  integer np

  integer i
  integer iline(my)
  integer indx(np,2)
  integer ip
  integer itemp
  integer iwrite
  logical long
  integer mx
  integer nodex0
  integer nx
  integer ny
  real ( kind = 8 ) xlngth
!
!  Determine the number of a node on the profile line
!
  itemp = nint((2.0D+00*(nx-1)*9.0D+00)/xlngth)
  if ( long ) then
    nodex0 = itemp*(2*ny-1)+1
  else
    nodex0 = itemp+1
  end if

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
    write (*,*) 'SETLIN: unknown numbers along line:'
    write(*,*)' '
    write (*,'(1X,15I5)') (iline(i),i = 1,my)
    write(*,*)' '
  end if

  return
end
subroutine setqud(area,nelemn,nnodes,node,np,nquad,xc,xm,yc,ym)

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
!    20 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer nelemn
  integer nnodes
  integer np
  integer nquad

  real ( kind = 8 ) area(nelemn)
  integer ip1
  integer ip2
  integer ip3
  integer it
  integer node(nelemn,nnodes)
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
    xm(it,1) = 0.5D+00*(x1+x2)
    xm(it,2) = 0.5D+00*(x2+x3)
    xm(it,3) = 0.5D+00*(x3+x1)
    ym(it,1) = 0.5D+00*(y1+y2)
    ym(it,2) = 0.5D+00*(y2+y3)
    ym(it,3) = 0.5D+00*(y3+y1)
    area(it) = 0.5D+00*abs((y1+y2)*(x2-x1)+(y2+y3)*(x3-x2) &
      +(y3+y1)*(x1-x3))
  end do

  return
end
subroutine setxy(iwrite,long,mx,my,np,nx,ny,xc,xlngth,yc,ylngth)

!*****************************************************************************80
!
!! SETXY sets the X, Y coordinates of grid points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer np

  integer i
  integer ic
  integer ip
  integer iwrite
  integer jc
  logical long
  integer mx
  integer my
  integer nx
  integer ny
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) xlngth
  real ( kind = 8 ) yc(np)
  real ( kind = 8 ) ylngth
!
!  Construct grid coordinates
!
  do ip = 1,np

    if ( long ) then
      ic = ((ip-1)/my)+1
      jc = mod((ip-1),my)+1
    else
      ic = mod((ip-1),mx)+1
      jc = ((ip-1)/mx)+1
    end if

    xc(ip) = (ic-1)*xlngth/(2*nx-2)
    yc(ip) = (jc-1)*ylngth/(2*ny-2)

  end do

  if ( 2 <= iwrite ) then
    write(*,*)' '
    write(*,*)'    I      XC           YC'
    write(*,*)' '
    do i = 1,np
      write (*,'(i5,2f12.5)')i,xc(i),yc(i)
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

  character ( len = 8  ) ampm
  integer   ( kind = 4 ) d
  integer   ( kind = 4 ) h
  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) mm
  character ( len = 9  ), parameter, dimension(12) :: month = (/ &
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
function ubdry ( y, para )

!*****************************************************************************80
!
!! UBDRY sets the parabolic inflow in terms of the value of the parameter.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) para
  real ( kind = 8 ) ubdry
  real ( kind = 8 ) y

  ubdry = 4.0D+00 * para * y * ( 3.0D+00 - y ) / 9.0D+00

  return
end
subroutine uval ( g, indx, iquad, it, nelemn, neqn, nnodes, node, &
  np, nquad, para, phi, un, unx, uny, yc )

!*****************************************************************************80
!
!! UVAL evaluates the velocities at a given point in a particular triangle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer nelemn
  integer neqn
  integer nnodes
  integer np
  integer nquad

  real ( kind = 8 ) bb
  real ( kind = 8 ) bx
  real ( kind = 8 ) by
  real ( kind = 8 ) g(neqn)
  integer indx(np,2)
  integer ip
  integer iq
  integer iquad
  integer it
  integer iuk
  integer iun
  integer node(nelemn,nnodes)
  real ( kind = 8 ) para
  real ( kind = 8 ) phi(nelemn,nquad,nnodes,3)
  real ( kind = 8 ) ubc
  real ( kind = 8 ) ubdry
  real ( kind = 8 ) un(2)
  real ( kind = 8 ) unx(2)
  real ( kind = 8 ) uny(2)
  real ( kind = 8 ) yc(np)

  un(1:2) = 0.0D+00
  uny(1:2) = 0.0D+00
  unx(1:2) = 0.0D+00

  do iq = 1, nnodes

    ip = node(it,iq)
    bb = phi(it,iquad,iq,1)
    bx = phi(it,iquad,iq,2)
    by = phi(it,iquad,iq,3)

    do iuk = 1, 2

      iun = indx(ip,iuk)

      if ( 0 < iun ) then
        un(iuk) = un(iuk)+bb*g(iun)
        unx(iuk) = unx(iuk)+bx*g(iun)
        uny(iuk) = uny(iuk)+by*g(iun)
      else if (iun < 0) then
        ip = node(it,iq)
        ubc = ubdry(yc(ip),para)
        un(iuk) = un(iuk)+bb*ubc
        unx(iuk) = unx(iuk)+bx*ubc
        uny(iuk) = uny(iuk)+by*ubc
      end if

    end do

  end do

  return
end
subroutine uv_plot3d (f,indx,insc,ivunit,long,mx,my, &
  nelemn,neqn,nnodes,node,np,para,press,reynld,yc)

!*****************************************************************************80
!
!! UV_PLOT3D creates a velocity file for use by PLOT3D.
!
!  Given the following set of nodes:
!
!    A  B  C
!    D  E  F
!    G  H  I
!
!  the file will have the form:
!
!    D, U(G), V(G), P
!    D, U(H), V(H), P
!    D, U(I), V(I), P
!    D, U(D), V(D), P
!    D, U(E), V(E), P
!    D, U(F), V(F), P
!    D, U(A), V(A), P
!    D, U(B), V(B), P
!    D, U(C), V(C), P
!
!  Here both D and P are set to 1 for now, representing dummy values
!  of density and pressure.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) PARA, the value of the parameter.
!
!    Workspace, real ( kind = 8 ) PRESS(MX,MY), used to hold the
!    computed pressures.
!
  implicit none

  integer nelemn
  integer neqn
  integer nnodes
  integer np

  real ( kind = 8 ) alpha
  real ( kind = 8 ) dval
  real ( kind = 8 ) f(neqn)
  real ( kind = 8 ) fsmach
  integer i
  integer ii
  integer indx(np,2)
  integer insc(np)
  integer ip
  integer, save :: iset = 0
  integer ivunit
  integer j
  integer k
  logical long
  integer mx
  integer my
  integer node(nelemn,nnodes)
  real ( kind = 8 ) para
  real ( kind = 8 ) press(mx,my)
  real ( kind = 8 ) reynld
  real ( kind = 8 ) time
  real ( kind = 8 ) ubdry
  real ( kind = 8 ) uval
  real ( kind = 8 ) vval
  real ( kind = 8 ) yc(np)

  iset = iset+1

  dval = 1.0D+00
  fsmach = 1.0D+00
  alpha = 1.0D+00
  time = 1.0D+00

  call pval (f,insc,long,mx,my,nelemn,neqn,nnodes,node,np,press)
!
!  If NY < NX, then nodes with a constant Y value are numbered consecutively.
!
  if ( long ) then
    write(ivunit,'(2I5)')mx,my
    write(ivunit,'(4G15.5)')fsmach,alpha,reynld,time
    do ii = 1,4
      do j = 1,my
        do i = 1,mx
          ip = (i-1)*my+j
          if ( ii == 1 ) then
            write(ivunit,'(G15.5)')dval
          else if ( ii == 2 ) then
            k = indx(ip,1)
            if (k == 0) then
              uval = 0.0
            else if (k < 0) then
              uval = ubdry(yc(ip),para)
            else
              uval = f(k)
            end if
            write(ivunit,'(G15.5)')uval
          else if ( ii == 3 ) then
            k = indx(ip,2)
            if (k == 0) then
              vval = 0.0D+00
            else
              vval = f(k)
            end if
            write(ivunit,'(G15.5)')vval
          else
            write(ivunit,'(G15.5)')press(i,j)
          end if
        end do
      end do
    end do
!
!  If NX < NY, then nodes with a constant X value are numbered consecutively.
!
  else
    write(ivunit,'(2I5)')mx,my
    write(ivunit,'(4G15.5)')fsmach,alpha,reynld,time
    do ii = 1,4
      do i = 1,mx
        do j = 1,my
          if ( ii == 1 ) then
             write(ivunit,'(G15.5)')dval
          else if ( ii == 2 ) then
            ip = (i-1)*my+j
            k = indx(ip,1)
            if (k == 0) then
              uval = 0.0D+00
            else if (k < 0) then
              uval = ubdry(yc(i),para)
            else
              uval = f(k)
            end if
            write(ivunit,'(G15.5)')uval
          else if ( ii == 3 ) then
            k = indx(ip,2)
            if (k == 0) then
              vval = 0.0D+00
            else
              vval = f(k)
            end if
            write(ivunit,'(G15.5)')vval
          else
            write(ivunit,'(G15.5)')press(i,j)
          end if
        end do
      end do
    end do
  end if

  write (*,*) 'UV_PLOT3D wrote data set ',iset,' to file.'

  return
end
subroutine uv_table ( f, indx, ivunit, neqn, np, para, yc )

!*****************************************************************************80
!
!! UV_TABLE creates a velocity table file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 February 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer neqn
  integer np

  real ( kind = 8 ) f(neqn)
  integer i
  integer indx(np,2)
  integer ip
  integer ivunit
  integer k
  real ( kind = 8 ) para
  real ( kind = 8 ) ubdry
  real ( kind = 8 ) uval
  real ( kind = 8 ) vval
  real ( kind = 8 ) yc(np)

  do ip = 1, np

    k = indx(ip,1)
    if ( k == 0 ) then
      uval = 0.0D+00
    else if ( k < 0 ) then
      uval = ubdry ( yc(ip), para )
    else
      uval = f(k)
    end if

    k = indx(ip,2)
    if ( k == 0 ) then
      vval = 0.0D+00
    else
      vval = f(k)
    end if

    write ( ivunit, '(2x,g14.6,2x,g14.6)' ) uval, vval

  end do

  return
end
subroutine xy_plot3d ( ixunit, long, np, nx, ny, xc, yc )

!*****************************************************************************80
!
!! XY_PLOT3D creates a grid file for use by PLOT3D.
!
!  Given the following set of nodes:
!
!    A  B  C
!    D  E  F
!    G  H  I
!
!  the file will have the form:
!
!    X(G), X(H), X(I), X(D), X(E), X(F), X(A), X(B), X(C),
!    Y(G), Y(H), Y(I), Y(D), Y(E), Y(F), Y(A), Y(B), Y(C).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer np

  integer i
  integer ip
  integer, save :: iset = 0
  integer ixunit
  integer j
  logical long
  integer mx
  integer my
  integer nx
  integer ny
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) yc(np)

  iset = iset+1

  mx = 2*nx-1
  my = 2*ny-1
!
!  If NY < NX, then nodes with a constant Y value are numbered consecutively.
!
  if ( long ) then
    write(ixunit,'(2I15)')mx,my
    do i = 1,my
      do j = 1,mx
        ip = (j-1)*my+i
        write(ixunit,'(G15.5)')xc(ip)
      end do
    end do

    do i = 1,my
      do j = 1,mx
        ip = (j-1)*my+i
        write(ixunit,'(G15.5)')yc(ip)
      end do
    end do
!
!  If NX < NY, then nodes with a constant X value are numbered consecutively.
!
  else

    write(ixunit,'(2I15)')my,mx

    do j = 1,mx
      do i = 1,my
        ip = (j-1)*my+i
        write(ixunit,'(G15.5)')xc(ip)
      end do
    end do

    do j = 1,mx
      do i = 1,my
        ip = (j-1)*my+i
        write(ixunit,'(G15.5)')yc(ip)
      end do
    end do

  end if

  write (*,*) 'XYDUMP wrote data set ',iset,' to file.'

  return
end
subroutine xy_table ( ixunit, np, xc, yc )

!*****************************************************************************80
!
!! XY_TABLE creates an XY table file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 February 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer np

  integer ip
  integer ixunit
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) yc(np)

  do ip = 1, np
    write ( ixunit, '(2x,g14.6,2x,g14.6)' ) xc(ip), yc(ip)
  end do

  return
end
