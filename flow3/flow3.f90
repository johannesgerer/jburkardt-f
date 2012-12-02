program main

!*****************************************************************************80
!
!! MAIN is the main program for FLOW3.
!
!  Discussion:
!
!    FLOW3 controls a package that solves a fluid flow problem.
!
!    MAXNX =    21,    31,    41,    61,     81,     121,      161,
!    MAXNY =     7,    10,    13,    19,     25,      37,       49,
!    H =       1/4,   1/6,   1/8,  1/12,   1/16,    1/24,     1/32,
!           0.25, 0.166, 0.125, 0.083, 0.0625, 0.04166,  0.03125,
!
!    The assignment of MAXROW should really read (maxrow = 28*min(nx,ny)).
!
!  Usage:
!
!    flow3 < input_file > output_file
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Max Gunzburger,
!    Finite Element Methods for Viscous Incompressible Flows,
!    A Guide to Theory, Practice, and Algorithms,
!    Academic Press, 1989,
!    ISBN: 0-12-307350-2,
!    LC: TA357.G86.
!
  implicit none

  integer ( kind = 4 ), parameter :: maxnx = 41
  integer ( kind = 4 ), parameter :: maxny = 13

  integer ( kind = 4 ), parameter :: liv = 60
  integer ( kind = 4 ), parameter :: maxdim = 3
  integer ( kind = 4 ), parameter :: maxparb = 5
  integer ( kind = 4 ), parameter :: maxparf = 5
  integer ( kind = 4 ), parameter :: npe = 6

  integer ( kind = 4 ), parameter :: maxrow = 28 * maxny
  integer ( kind = 4 ), parameter :: maxelm = 2 * ( maxnx - 1 ) * ( maxny - 1 )
  integer ( kind = 4 ), parameter :: maxeqn = 2*(2*maxnx-1)*(2*maxny-1)+maxnx*maxny
  integer ( kind = 4 ), parameter :: maxnp = (2*maxnx-1)*(2*maxny-1)
  integer ( kind = 4 ), parameter :: maxpar = maxparb + maxparf + 1

  integer ( kind = 4 ), parameter :: lv = 78+maxpar*(maxpar+21)/2

  real ( kind = 8 ) a(maxrow,maxeqn)
  real ( kind = 8 ) area(maxelm,3)
  real ( kind = 8 ) base(maxpar,maxdim)
  real ( kind = 8 ) cost
  real ( kind = 8 ) costar
  real ( kind = 8 ) costb
  real ( kind = 8 ) costp
  real ( kind = 8 ) costu
  real ( kind = 8 ) costv
  real ( kind = 8 ) dir(maxpar,maxdim)
  real ( kind = 8 ) disjac
  real ( kind = 8 ) dopt(maxpar)
  real ( kind = 8 ) dpara3(maxpar)
  real ( kind = 8 ) dparsn(maxpar)
  real ( kind = 8 ) dparfd(maxpar)
  real ( kind = 8 ) dparfdc(maxpar)
  real ( kind = 8 ) dpdyn(maxnp)
  real ( kind = 8 ) dudyn(maxnp)
  real ( kind = 8 ) dvdyn(maxnp)
  real ( kind = 8 ) dydpn(maxnp,maxparb)
  character ( len = 2 ) eqn(maxeqn)
  real ( kind = 8 ) etan(6)
  real ( kind = 8 ) etaq(3)
  character ( len = 80 ) plot_file
  character ( len = 80 ) march_file
  real ( kind = 8 ) flarea
  real ( kind = 8 ) g(maxeqn)
  real ( kind = 8 ) gdif(maxeqn,maxpar)
  real ( kind = 8 ) gdifc(maxeqn,maxpar)
  real ( kind = 8 ) gold(maxeqn)
  real ( kind = 8 ) gopt(maxpar)
  real ( kind = 8 ) gradf(maxeqn,maxpar)
  real ( kind = 8 ) gtar(maxeqn)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibc
  integer ( kind = 4 ) ibump
  integer ( kind = 4 ) idfd
  integer ( kind = 4 ) ids
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ifds
  integer ( kind = 4 ) igrad
  integer ( kind = 4 ) igrid
  integer ( kind = 4 ) igunit
  integer ( kind = 4 ) ijac
  integer ( kind = 4 ) indx(maxnp,3)
  integer ( kind = 4 ) iopt(maxpar)
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) ipivot(maxeqn)
  integer ( kind = 4 ) iplot
  integer ( kind = 4 ) ipred
  integer ( kind = 4 ) iprdat
  integer ( kind = 4 ) ishapb
  integer ( kind = 4 ) ishapbt
  integer ( kind = 4 ) ishapf
  integer ( kind = 4 ) ishapft
  integer ( kind = 4 ) ismooth
  integer ( kind = 4 ) isotri(maxelm)
  integer ( kind = 4 ) istep1
  integer ( kind = 4 ) istep2
  integer ( kind = 4 ) itar
  integer ( kind = 4 ) itemp
  integer ( kind = 4 ) itunit
  integer ( kind = 4 ) itype
  integer ( kind = 4 ) iuval
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) ivopt(liv)
  integer ( kind = 4 ) iwrite
  integer ( kind = 4 ) jjac
  integer ( kind = 4 ) jstep1
  integer ( kind = 4 ) jstep2
  logical lmat
  logical lval
  integer ( kind = 4 ) maxnew
  integer ( kind = 4 ) maxstp
!
!  TEMPORARY
!
  integer ( kind = 4 ) nabor(19,maxnp)
  integer ( kind = 4 ) ndim
  integer ( kind = 4 ) nelem
  integer ( kind = 4 ) neqn
  integer ( kind = 4 ) nlband
  integer ( kind = 4 ) node(maxelm,npe)
  integer ( kind = 4 ) nopt
  integer ( kind = 4 ) np
  integer ( kind = 4 ) npar
  integer ( kind = 4 ) nparb
  integer ( kind = 4 ) nparbt
  integer ( kind = 4 ) nparf
  integer ( kind = 4 ) nprof(2*maxny-1)
  integer ( kind = 4 ) nrow
  integer ( kind = 4 ) nstep3
  integer ( kind = 4 ) numel(maxnp)
  integer ( kind = 4 ) numstp
  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny
  real ( kind = 8 ) para(maxpar)
  real ( kind = 8 ) para1(maxpar)
  real ( kind = 8 ) para2(maxpar)
  real ( kind = 8 ) para3(maxpar)
  real ( kind = 8 ) parjac(maxpar)
  real ( kind = 8 ) parnew(maxpar)
  real ( kind = 8 ) parold(maxpar)
  real ( kind = 8 ) partar(maxpar)
  real ( kind = 8 ) phi(maxelm,3,6,10)
  real ( kind = 8 ) res(maxeqn)
  real ( kind = 8 ) sens(maxeqn,maxpar)
  real ( kind = 8 ) splbmp(4,maxparb+2,0:maxparb)
  real ( kind = 8 ) splflo(4,maxparf+2,0:maxparf)
  real ( kind = 8 ) stpmax
  character ( len = 20 ) syseqn
  real ( kind = 8 ) taubmp(maxparb+2)
  real ( kind = 8 ) tauflo(maxparf+2)
  character ( len = 80 ) title
  real ( kind = 8 ) tolnew
  real ( kind = 8 ) tolopt
  real ( kind = 8 ) u1
  real ( kind = 8 ) u1max
  real ( kind = 8 ) u2
  real ( kind = 8 ) u2max
  real ( kind = 8 ) vopt(lv)
  real ( kind = 8 ) wateb
  real ( kind = 8 ) wateb1
  real ( kind = 8 ) wateb2
  real ( kind = 8 ) watep
  real ( kind = 8 ) wateu
  real ( kind = 8 ) watev
  real ( kind = 8 ) wquad(3)
  real ( kind = 8 ) xbl
  real ( kind = 8 ) xbleft
  real ( kind = 8 ) xbltar
  real ( kind = 8 ) xbord(maxnx)
  real ( kind = 8 ) xbr
  real ( kind = 8 ) xbrite
  real ( kind = 8 ) xbrtar
  real ( kind = 8 ) xc(maxnp)
  real ( kind = 8 ) xopt(maxpar)
  real ( kind = 8 ) xquad(maxelm,3)
  real ( kind = 8 ) xprof
  real ( kind = 8 ) xsin(6)
  real ( kind = 8 ) xsiq(3)
  real ( kind = 8 ) y2max
  real ( kind = 8 ) ybl
  real ( kind = 8 ) ybleft
  real ( kind = 8 ) ybltar
  real ( kind = 8 ) ybord(maxny)
  real ( kind = 8 ) ybr
  real ( kind = 8 ) ybrite
  real ( kind = 8 ) ybrtar
  real ( kind = 8 ) yc(maxnp)
  real ( kind = 8 ) yquad(maxelm,3)
  real ( kind = 8 ) yval

  call timestamp ( )

  iprdat = 0
  ival = 0
!
!  Print the program name, date, and computer name.
!
  write ( *, * ) ' '
  write ( *, * ) 'FLOW3'
  write ( *, * ) '  FORTRAN90 version'
  write ( *, * ) '  This is the version of 01 January 2001.'
  write ( *, * ) '  The maximum problem size is ',maxnx,' by ',maxny
!
!  Initialize the variables.
!
  call init ( area,cost,costar,costb,costp,costu,costv,disjac,dopt,dpara3, &
    dparsn,dparfd,dparfdc,dpdyn,dudyn,dvdyn,dydpn,eqn,etan,plot_file, &
    march_file,g,gdif,gold,gopt,gradf,gtar,ibc,ibump,idfd,ids,ierror, &
    ifds,igrad,igrid,igunit,ijac,indx,iopt,iplot,ipred,ishapb,ishapbt, &
    ishapf,ishapft,ismooth,isotri,istep1,istep2,itar,itunit,itype,ivopt, &
    iwrite,jjac,jstep1,jstep2,liv,lv,maxelm,maxeqn,maxnew,maxnp,maxnx,maxny, &
    maxpar,maxparb,maxstp,node,nopt,npar,nparb,nparf,npe,nstep3,nx,ny,para1, &
    para2,para3,parjac,partar,sens,stpmax,syseqn,tolnew,tolopt,vopt,wateb, &
    wateb1,wateb2,watep,wateu,watev,xbleft,xbltar,xbord,xbrite,xbrtar,xprof, &
    xsin,ybleft,ybltar,ybord,ybrite,ybrtar)
!
!  Read the user input.
!
  call input ( disjac,plot_file,march_file,ibc,ibump,idfd,ids,ierror,ifds,igrad, &
    igrid,ijac,iopt,iplot,ipred,ishapb,ishapbt,ishapf,ishapft,ismooth,istep1, &
    istep2,itar,itype,iwrite,jjac,jstep1,jstep2,maxnew,maxnx,maxny,maxpar, &
    maxstp,npar,nparb,nparf,nstep3,nx,ny,para1,para2,para3,partar,stpmax, &
    syseqn,tolnew,tolopt,wateb,wateb1,wateb2,watep,wateu,watev,xbleft, &
    xbltar,xbord,xbrite,xbrtar,xprof,ybleft,ybltar,ybord,ybrite,ybrtar)

  if ( ierror == 0 ) then
    call imemry ( 'get', 'Points_Opt', ival )
    ival = 1
    call imemry ( 'inc', 'Restarts', ival )
    write ( *, * ) ' '
    write ( *, * ) 'Flow - GO signal from user input.'
  else
    write ( *, * ) ' '
    write ( *, * ) 'Flow - STOP signal from user.'
    call pr_work
    call after ( plot_file, march_file, igunit, itunit )
  end if
!
!  Set the number of optimization variables.
!
  nopt = sum ( iopt(1:npar) )
  nelem = 2 * ( nx - 1 ) * ( ny - 1 )
  np = ( 2 * nx - 1 ) * ( 2 * ny - 1 )

  write ( *, * ) ' '
  write ( *, * ) 'FLOW3 - Note.'
  write ( *, * ) '  Number of elements  = ', nelem
  write ( *, * ) '  Number of nodes =   ', np
!
!  Check the data.
!
  call chkdat ( ibump,idfd,ids,ifds,igrad,ijac,iopt,ipred,itype,jjac,maxpar, &
    maxparb,maxparf,nopt,npar,nparb,nparf,nstep3,para1,partar,xbleft,xbrite)
!
!  Print the information that the user can change via the input
!  file.  But only do this once.  The output from this routine
!  is extensive, and tedious to see repeatedly.
!
  if ( iprdat == 0 ) then

    iprdat = iprdat+1

    call pr_dat ( disjac,plot_file,march_file,ibc,ibump,idfd,ids,ifds,igrad, &
      igrid,ijac,iopt,iplot,ipred,ishapb,ishapbt,ishapf,ishapft,ismooth, &
      istep1,istep2,itar,itype,iwrite,jjac,jstep1,jstep2,maxnew,maxpar, &
      maxstp,nopt,npar,nparb,nparf,nstep3,nx,ny,para1,para2,para3,partar, &
      stpmax,syseqn,tolnew,tolopt,wateb,wateb1,wateb2,watep,wateu,watev, &
      xbleft,xbltar,xbord,xbrite,xbrtar,xprof,ybleft,ybltar,ybord,ybrite,ybrtar)

  end if
!
!  Open the plot file.
!
  call plot_file_open ( plot_file, igunit, iplot )
!
!  Open the marching file.
!
  call march_file_open ( march_file, itunit )
!
!  Solve the optimization problem,
!  which involves repeated evaluation of the functional
!  at various flow solutions (U,V,P)+(FLO,BUM,NU_INV).
!
!  Set data specific to target calculations.
!
  lval = .true.
  call lmemry('set','target',lval)

  xbl = xbltar
  xbr = xbrtar
  ybl = ybltar
  ybr = ybrtar
  parnew(1:npar) = partar(1:npar)
!
!  Set up the elements and nodes.
!
  call node_set ( eqn, ibump, indx, isotri, maxeqn, nelem, neqn, node, np, npe, &
    nx, ny, xbl, xbr )
!
!  TEMPORARY
!
!  Set the neighbor array.
!
  call setnab ( nabor, np, nx, ny )
!
!  Find points on velocity profile sampling line.
!
  itemp = nint ( ( 2.0D+00 * real ( nx - 1 ) * xprof ) / 10.0D+00 )

  do i = 1,2*ny-1
    ip = itemp*(2*ny-1)+i
    nprof(i) = ip
  end do
!
!  Get matrix bandwidth.
!
  call setban ( indx, maxrow, nelem, neqn, nlband, node, np, npe, nrow )
!
!  Demand that SETXY set internal node coordinates, even if IGRID = 2.
!
  lval = .true.
  call lmemry ( 'set', 'need_xy', lval )
!
!  Demand that SETBAS set the basis functions.
!
  lval = .true.
  call lmemry ( 'set', 'need_phi', lval )

  if ( itar == 0 ) then

    write ( *, * ) ' '
    write ( *, * ) 'FLOW3'
    write ( *, * ) '  Computing the target solution, GTAR.'

    para(1:npar) = partar(1:npar)

    call flosol ( a,area,disjac,eqn,etan,etaq,flarea,g,ierror,igrid,ijac, &
      indx,ipivot,ishapbt,ishapft,isotri,iwrite,jjac,maxnew,nelem,neqn,nlband, &
      node,np,npar,nparb,nparf,npe,nrow,nx,ny,para,parjac,phi,res,splbmp, &
      splflo, &
      syseqn,taubmp,tauflo,tolnew,wquad,xbl,xbord,xbr,xc,xquad,xsin,xsiq, &
      ybl,ybord,ybr,yc,yquad)

    if ( ierror /= 0 ) then
      write ( *, * ) ' '
      write ( *, * ) 'FLOW3 - Fatal error!'
      write ( *, * ) '  Could not compute the target solution.'
      stop
    end if

    gtar(1:neqn) = g(1:neqn)
!
!  Compute the finite coefficient differences dG/dPara for the target solution.
!
    call getgrd ( a,area,cost,disjac,dpara3,eqn,etan,etaq,flarea,g, &
      gdif,gtar,idfd,ierror,igrid,ijac,indx,iopt,ipivot,ishapbt,ishapft, &
      isotri,iwrite,jjac,maxnew,nelem,neqn,nlband,node,np,npar,nparb,nparf, &
      npe,nprof,nrow,nx,ny,para,parjac,phi,res,splbmp,splflo,syseqn,taubmp, &
      tauflo,tolnew,wateb,watep,wateu,watev,wquad,xbl,xbord,xbr,xc,xquad, &
      xsin,xsiq,ybl,ybord,ybr,yc,yquad)

    if ( ierror /= 0 ) then
      write ( *, * ) ' '
      write ( *, * ) 'FLOW3 - Fatal error!'
      write ( *, * ) '  GETGRD returns IERROR = ',ierror
      stop
    end if

    numstp = 0
!
!  Compute the cost COST associated with the solution GTAR.
!
    call get_cost ( cost, costb, costp, costu, costv, g, gtar, indx, &
      ishapbt, neqn, np, nparb, nprof, ny, splbmp, taubmp, wateb, watep, &
      wateu, watev, xbl, xbr, ybl, ybr, yc )

    if ( iwrite == 1 ) then
      title = 'Cost associated with target solution.'
      call pr_cost1 ( cost, title )
    else if ( iwrite >= 2 ) then
      title = 'Cost associated with target solution.'
      call pr_cost2 ( cost, costb, costp, costu, costv, title, wateb, watep, &
        wateu, watev )
    end if
!
!  Print stuff out.
!
    if ( ifds /= 0 ) then
      call cost_gradient ( dparfd,g,gtar,indx,ishapbt,neqn,np,npar,nparb, &
        nparf,nprof,ny, &
        gdif,splbmp,taubmp,wateb,watep,wateu,watev,xbl,xbr,ybl,ybr,yc)
    end if

    if ( ifds /= 0 ) then
      call cost_gradient ( dparfdc,g,gtar,indx,ishapbt,neqn,np,npar,nparb, &
        nparf,nprof, &
        ny,gdifc,splbmp,taubmp,wateb,watep,wateu,watev,xbl,xbr,ybl,ybr,yc)
    end if

    if ( ids /= 0 ) then
      call cost_gradient ( dparsn,g,gtar,indx,ishapbt,neqn,np,npar,nparb, &
        nparf,nprof,ny, &
        sens,splbmp,taubmp,wateb,watep,wateu,watev,xbl,xbr,ybl,ybr,yc)
    end if

    title = 'Computed target solution'

    call pr_solution ( dpara3, dparfd, dparfdc, dparsn, g, gdif, gdifc, gtar, &
      idfd, ids, ifds, indx, iwrite, neqn, np, npar, nparb, nparf, nprof, &
      numstp, ny, parnew, sens, title, yc )

    if ( iwrite >= 4 ) then
      call xy_print ( maxnp, np, xc, yc )
    end if

    costar = cost

  else if ( itar == 1 ) then

    write ( *, * ) ' '
    write ( *, * ) 'FLOW3'
    write ( *, * ) '  Calling GETTAR2 to get target data.'

    nparbt = 0
    ishapbt = 0

    call xy_set ( igrid, ishapbt, np, nparbt, nx, ny, splbmp, taubmp, &
      xbltar, xbord, xbrtar, xc, ybltar, ybord, ybrtar, yc )

    g(1:neqn) = 0.0D+00

    u1max = 9.0D+00 / 4.0D+00
    y2max = (6.0D+00-2.0D+00*sqrt(3.0D+00))/4.0D+00
    u2max = (3.0D+00-y2max)*(1.5D+00-y2max)*y2max
!
!  WARNING!  The only reason I can get away with referencing
!  INDX in the main program (where it has first dimension MAXNP)
!  is because I reference column 1.  Very lucky!
!
    do i = 1,2*ny-1
      iuval = indx(nprof(i),1)
      yval = real ( i - 1, kind = 8 )*3.0D+00/ real ( 2 * ny - 2, kind = 8 )
      yc(nprof(i)) = yval
      u1 = (3.0D+00-yval)*yval/u1max
      u2 = (3.0D+00-yval)*(1.5D+00-yval)*yval/u2max
      g(iuval) = partar(1)*u1+partar(2)*u2
    end do
!
!  Compute the cost, relative to zero.
!
    gtar(1:neqn) = 0.0D+00

    call disc_cost ( costp, costu, costv, g, gtar, indx, neqn, np, nprof, ny, yc )

    cost = costu
    costar = cost

    call pr_profile ( g, indx, neqn, np, nprof, ny, yc )

    gtar(1:neqn) = g(1:neqn)

  end if
!
!  Write target information to plot file.
!
  if ( iplot /= 0 ) then

    call plot_file_write ( eqn, g, gdif, igunit, indx, isotri, iwrite, nelem, &
      neqn, node, np, npar, npe, nprof, nx, ny, partar, sens, xc, xprof, yc )

  end if
!
!  Turn off target calculation flag.
!
  lval = .false.
  call lmemry ( 'set', 'target', lval )

  call pr_work
!
!  Zero stuff out.
!
  ival = 0
  call imemry('zero','nothing',ival)

  parjac(1:npar) = 0.0D+00

  lmat = .false.
  call lmemry('set','have_fp',lmat)
!
!  Demand that the internal nodes coordinates be reset by SETXY,
!  even if IGRID = 2, for the next calculation.
!
  lval = .true.
  call lmemry('set','need_xy',lval)

  lval = .true.
  call lmemry('set','need_phi',lval)
!
!  Set data specific to the optimization.
!
  xbl = xbleft
  xbr = xbrite
  ybl = ybleft
  ybr = ybrite
!
!  Set the elements and nodes.
!
  if ( xbleft/= xbltar .or. xbrite /= xbrtar .or. ybleft /= ybltar .or. &
    ybrite /= ybrtar ) then

    call node_set ( eqn,ibump,indx,isotri,maxeqn,nelem,neqn,node,np,npe, &
      nx,ny,xbl,xbr)

  end if
!
!  Is this a march?
!
  if ( itype == 1 .or. itype == 2 .or. itype == 4 ) then

    if ( itype == 1 ) then
      ndim = 1
    else if ( itype == 2 ) then
      ndim = 2
    else if ( itype == 4 ) then
      ndim = 3
    end if

    call march ( a,area,base,dir,disjac,dpara3,dparsn,dparfd,dparfdc,dpdyn, &
      dudyn,dvdyn,dydpn,eqn,etan,etaq,g,gdif,gdifc,gold, &
      gradf,gtar,ibc,idfd,ids,ifds,igrid,igunit,ijac,indx,iopt,ipivot,iplot, &
      ipred,ishapb,ishapf,ismooth,isotri,istep1,istep2,itunit,itype,iwrite, &
      jjac,jstep1,jstep2,maxnew,maxstp,ndim,nelem,neqn,nlband,node,np,npar, &
      nparb,nparf,npe,nprof,nrow,nstep3,numel,nx,ny,para,para1,para2,para3, &
      parjac, &
      parnew,parold,phi,res,sens,splbmp,splflo,stpmax,syseqn,taubmp,tauflo, &
      tolnew,wateb,wateb1,wateb2,watep,wateu,watev,wquad,xbl,xbord,xbr,xc, &
      xprof,xquad,xsin,xsiq,ybl,ybord,ybr,yc,yquad)
!
!  ...or a sensitivity optimization?
!
  else if ( itype == 3 ) then

    call qsolve ( a,area,costar,disjac,dopt,dpara3,dparfd,dparfdc,dparsn,dpdyn, &
      dudyn,dvdyn,dydpn,eqn,etan,etaq,g,gdif,gdifc,gold,gopt, &
      gradf,gtar,ibc,idfd,ids,ifds,igrad,igrid,igunit,ijac,indx,iopt,ipivot, &
      iplot,ipred,ishapb,ishapf,ismooth,isotri,itunit,itype,ivopt,iwrite, &
      jjac,liv,lv,maxnew,maxstp,nelem,neqn,nlband,node,nopt,np,npar,nparb, &
      nparf,npe,nprof,nrow,numel,nx,ny,para,para1,parjac,parnew,partar, &
      phi,res,sens,splbmp,splflo,stpmax,syseqn,taubmp,tauflo,tolnew,tolopt, &
      vopt,wateb,watep,wateu,watev,wquad,xbl,xbord,xbr,xc,xopt,xprof,xquad, &
      xsin,xsiq,ybl,ybord,ybr,yc,yquad)
!
!  ...or a one point computation?
!
  else if ( itype == 5 ) then

    call rsolve ( a,area,disjac,dpara3,dparfd,dparfdc,dparsn,dpdyn,dudyn, &
      dvdyn,dydpn,eqn,etan,etaq,g,gdif,gdifc,gradf,gtar, &
      ibc,idfd,ids,ifds,igrid,igunit,ijac,indx,iopt,ipivot,iplot,ipred, &
      ishapb,ishapf,ismooth,isotri,itype,iwrite,jjac,maxnew,nelem,neqn, &
      nlband,node,np,npar,nparb,nparf,npe,nprof,nrow,numel,nx,ny,para,para1, &
      parjac,phi,res,sens,splbmp,splflo,stpmax,syseqn,taubmp, &
      tauflo,tolnew,wateb,watep,wateu,watev,wquad,xbl,xbord,xbr,xc,xprof, &
      xquad,xsin,xsiq,ybl,ybord,ybr,yc,yquad)
!
!  ...or an optimization using function values only?
!
  else if ( itype == 7 ) then

    call osolve ( a,area,disjac,dopt,dpara3,dparfd,dparfdc,dparsn,dpdyn,dudyn, &
      dvdyn,dydpn,eqn,etan,etaq,g,gdif,gdifc,gold,gradf, &
      gtar,ibc,idfd,ids,ifds,igrid,igunit,ijac,indx,iopt,ipivot,iplot,ipred, &
      ishapb,ishapf,ismooth,isotri,itunit,itype,ivopt,iwrite,jjac,liv,lv, &
      maxnew,maxstp,nelem,neqn,nlband,node,nopt,np,npar,nparb,nparf,npe, &
      nprof, &
      nrow,numel,nx,ny,para,para1,parjac,parnew,parold,phi,res,sens,splbmp, &
      splflo,stpmax,syseqn,taubmp,tauflo,tolnew,tolopt,vopt,wateb,watep, &
      wateu,watev,wquad,xbl,xbord,xbr,xc,xopt,xprof,xquad,xsin,xsiq,ybl, &
      ybord,ybr,yc,yquad)

  else

    write ( *, * ) ' '
    write ( *, * ) 'FLOW3 - Fatal error!'
    write ( *, * ) '  Unknown value of ITYPE = ', itype

  end if

  call pr_work

  call after ( plot_file, march_file, igunit, itunit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FLOW3'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine after ( plot_file, march_file, igunit, itunit )

!*****************************************************************************80
!
!! AFTER shuts down files and stops.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
  implicit none

  character ( len = * ) plot_file
  character ( len = * ) march_file
  integer ( kind = 4 ) igunit
  integer ( kind = 4 ) itunit
!
!  Close the marching file.
!
  if ( itunit /= 0 ) then
    close ( unit = itunit )
    write ( *, * ) ' '
    write ( *, * ) 'AFTER - Closing the marching file "' // &
      trim ( march_file ) // '".'
  end if
!
!  Close the graphics file.
!
  if ( igunit /= 0 ) then
    close ( unit = igunit )
    write ( *, * ) ' '
    write ( *, * ) 'AFTER - Closing the graphics file "' // &
      trim ( plot_file ) // '".'
  end if

  return
end
subroutine bsp ( q, dqdx, dqdy, ielem, iq, nelem, node, np, npe, xc, xq, yc, &
  yq )

!*****************************************************************************80
!
!! BSP evaluates the basis functions associated with pressure.
!
!  Discussion:
!
!    Here is a picture of a typical finite element associated with pressure:
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
!    01 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) Q, the value of the IQ-th basis
!    function at the point with global coordinates (XQ,YQ).
!
!    Output, real ( kind = 8 ) DQDX, DQDY, the X and Y
!    derivatives of the IQ-th basis function at the point
!    with global coordinates (XQ,YQ).
!
!    Input, integer ( kind = 4 ) IELEM, the global element number about which
!    we are inquiring.
!
!    Input, integer ( kind = 4 ) IQ, the index of the desired basis
!    function.  This is also the node of the reference
!    triangle which is associated with the basis function.
!
!    Basis function IQ is 1 at node IQ, and zero at the
!    other two nodes.
!
!    Input, integer ( kind = 4 ) NELEM, the number of elements.
!
!    Input, integer ( kind = 4 ) NODE(MAXELM,6).  NODE(I,J) is
!    the global node number of the J-th node in the I-th element.
!
!    Input, integer ( kind = 4 ) NP, the number of nodes.
!
!    Input, real ( kind = 8 ) XC(NP), the global X coordinates
!    of the element nodes.
!
!    Input, real ( kind = 8 ) XQ, the global X coordinate of
!    the point in which we are interested.
!
!    Input, real ( kind = 8 ) YC(NP), the global Y coordinates
!    of the element nodes.
!
!    Input, real ( kind = 8 ) YQ, the global Y coordinate of
!    the point in which we are interested.
!
  implicit none

  integer ( kind = 4 ) nelem
  integer ( kind = 4 ) np
  integer ( kind = 4 ) npe

  real ( kind = 8 ) q
  real ( kind = 8 ) dqdx
  real ( kind = 8 ) dqdy
  real ( kind = 8 ) d
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) ielem
  integer ( kind = 4 ) iq
  integer ( kind = 4 ) iq1
  integer ( kind = 4 ) iq2
  integer ( kind = 4 ) iq3
  integer ( kind = 4 ) node(nelem,npe)
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) xq
  real ( kind = 8 ) yc(np)
  real ( kind = 8 ) yq

  if ( iq < 1 .or. npe < iq ) then
    write ( *, * ) ' '
    write ( *, * ) 'BSP - Fatal error!'
    write ( *, * ) '  The requested basis function is IQ = ',iq
    write ( *, * ) '  but only values from 1 to 6 are legal.'
    stop
  else if ( iq >= 4 .and. iq <= npe ) then
    q = 0.0D+00
    dqdx = 0.0D+00
    dqdy = 0.0D+00
    return
  end if

  iq1 = iq
  iq2 = mod ( iq, 3 ) + 1
  iq3 = mod ( iq+1, 3 ) + 1

  i1 = node(ielem,iq1)
  i2 = node(ielem,iq2)
  i3 = node(ielem,iq3)

  d = ( xc(i2) - xc(i1) ) * ( yc(i3) - yc(i1) ) &
    - ( xc(i3) - xc(i1) ) * ( yc(i2) - yc(i1) )

  dqdx = ( yc(i2) - yc(i3) ) / d
  dqdy = ( xc(i3) - xc(i2) ) / d

  q = 1.0D+00 + dqdx * ( xq - xc(i1) ) + dqdy * ( yq - yc(i1) )

  return
end
subroutine bump_bc ( dudyn,dvdyn,ip,iparb,ishapb,np,nparb,shape,splbmp,taubmp, &
  ubc,vbc,xc)

!*****************************************************************************80
!
!! BUMP_BC computes the boundary conditions for velocity sensitivities.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
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
!    Input, integer ( kind = 4 ) IP, the global node number.
!
!    Input, integer ( kind = 4 ) IPARB, the bump parameter with respect to
!    which the sensitivities are desired.
!
!    Output, real ( kind = 8 ) UBC, the boundary condition for the 
!    horizontal velocity sensitivity.
!
!    Output, real ( kind = 8 ) VBC, the boundary condition for the vertical 
!    velocity sensitivity.
!
  implicit none

  integer ( kind = 4 ) np

  real ( kind = 8 ) dudyn(np)
  real ( kind = 8 ) dvdyn(np)
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) iparb
  integer ( kind = 4 ) ishapb
  integer ( kind = 4 ) jderiv
  integer ( kind = 4 ) nparb
  real ( kind = 8 ) shape
  real ( kind = 8 ) splbmp(4,nparb+2,0:nparb)
  real ( kind = 8 ) taubmp(nparb+2)
  real ( kind = 8 ) ubc
  real ( kind = 8 ) vbc
  real ( kind = 8 ) xc(np)
!
!  Determine the value of the "basis" shape function associated
!  with parameter IPARB, evaluated at node IP.
!
  if ( ishapb == 1 ) then
    call plval1 ( iparb+1, nparb+2, xc(ip), taubmp, shape )
  else if ( ishapb == 2 ) then
    call pqval1 ( iparb+1, nparb+2, xc(ip), taubmp, shape )
  else if ( ishapb == 3 ) then
    jderiv = 0
    call ppvalu ( taubmp, splbmp(1,1,iparb), nparb+1, 4, xc(ip), jderiv, shape )
  end if

  ubc = - dudyn(ip) * shape
  vbc = - dvdyn(ip) * shape

  return
end
subroutine bump_bc1 ( g,indx,ip,iparb,ishapb,neqn,np,nparb,shape,splbmp,taubmp, &
  ubc,vbc,xc,yc)

!*****************************************************************************80
!
!! BUMP_BC1 computes the value of the boundary conditions for the
!  horizontal and vertical velocity sensitivities with respect
!  to a given shape parameter using a two point finite difference
!  formula.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
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
!    Input, real ( kind = 8 ) G(MAXEQN), the finite element coefficients.
!
!    Input, integer ( kind = 4 ) IP, the global node number.
!
!    Input, integer ( kind = 4 ) IPARB, the bump parameter with respect to
!    which the sensitivities are desired.
!
!    Output, real ( kind = 8 ) UBC, the boundary condition for the
!    horizontal velocity sensitivity.
!
!    Output, real ( kind = 8 ) VBC, the boundary condition for the vertical 
!    velocity sensitivity.
!
  implicit none

  integer ( kind = 4 ) neqn
  integer ( kind = 4 ) np

  real ( kind = 8 ) dudy
  real ( kind = 8 ) dvdy
  real ( kind = 8 ) g(neqn)
  integer ( kind = 4 ) ihor
  integer ( kind = 4 ) ihorn
  integer ( kind = 4 ) indx(np,3)
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) iparb
  integer ( kind = 4 ) ishapb
  integer ( kind = 4 ) iver
  integer ( kind = 4 ) ivern
  integer ( kind = 4 ) jderiv
  integer ( kind = 4 ) nparb
  real ( kind = 8 ) shape
  real ( kind = 8 ) splbmp(4,nparb+2,0:nparb)
  real ( kind = 8 ) taubmp(nparb+2)
  real ( kind = 8 ) ubc
  real ( kind = 8 ) vbc
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) yc(np)
!
!  Determine the value of the "basis" shape function associated
!  with parameter IPARB, evaluated at node IP.
!
  if ( ishapb == 1 ) then
    call plval1 ( iparb+1, nparb+2, xc(ip), taubmp, shape )
  else if ( ishapb == 2 ) then
    call pqval1 ( iparb+1, nparb+2, xc(ip), taubmp, shape )
  else if ( ishapb == 3 ) then
    jderiv = 0
    call ppvalu ( taubmp, splbmp(1,1,iparb), nparb+1, 4, xc(ip), jderiv, shape )
  end if
!
!  Estimate dUdY and dVdY at the node.
!
  ihor = indx(ip,1)
  ihorn = indx(ip+1,1)
  dudy = (g(ihorn)-g(ihor))/(yc(ip+1)-yc(ip))

  iver = indx(ip,2)
  ivern = indx(ip+1,2)
  dvdy = (g(ivern)-g(iver))/(yc(ip+1)-yc(ip))
!
!  Set the boundary conditions.
!
  ubc = - dudy * shape
  vbc = - dvdy * shape

  return
end
subroutine bump_bc2 ( g, indx, ip, iparb, ishapb, neqn, np, nparb, shape, &
  splbmp, taubmp, ubc, vbc, xc, yc )

!*****************************************************************************80
!
!! BUMP_BC2 evaluates the velocity sensitivity boundary conditions.
!
!  Discussion:
!
!    The routine evaluates the boundary conditions for the
!    horizontal and vertical velocity sensitivities with respect
!    to a given shape parameter using a three point finite difference
!    formula.
!
!    We derive that formula from the following Taylor series:
!
!      u0 = u0
!      u1 = u0 + h1 u0' + h1**2 u0"/2 + O(h1**3)
!      u2 = u0 + h2 u0' + h2**2 u0"/2 + O(h2**3)
!
!    arriving at:
!
!      u0' = (h1**2 u2 - h2**2 u1 + (h2**2-h1**2) u0) / (h1*h2*(h1-h2))
!        + O(max(h1,h2)**2)
!
!    Note that these Taylor series are really only valid when all three
!    nodes lie in one element, so that U is a smooth polynomial.  This
!    is not true for midside nodes, though to correct this would be
!    painful.
!
!    When all three nodes do lie in one element, and if u is at
!    most a quadratic function, then the formula should be exact.
!    This is the case for the velocities U and V at a node.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) G(MAXEQN), the finite element coefficients.
!
!    Input, integer ( kind = 4 ) IP, the global node number.
!
!    Input, integer ( kind = 4 ) IPARB, the bump parameter with respect to
!    which the sensitivities are desired.
!
!    Output, real ( kind = 8 ) UBC, the boundary condition for the 
!    horizontal velocity sensitivity.
!
!    Output, real ( kind = 8 ) VBC, the boundary condition for the 
!    vertical velocity sensitivity.
!
  implicit none

  integer ( kind = 4 ) neqn
  integer ( kind = 4 ) np

  real ( kind = 8 ) dudy
  real ( kind = 8 ) dvdy
  real ( kind = 8 ) g(neqn)
  real ( kind = 8 ) h1
  real ( kind = 8 ) h2
  integer ( kind = 4 ) indx(np,3)
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) iparb
  integer ( kind = 4 ) ishapb
  integer ( kind = 4 ) jderiv
  integer ( kind = 4 ) nparb
  real ( kind = 8 ) shape
  real ( kind = 8 ) splbmp(4,nparb+2,0:nparb)
  real ( kind = 8 ) taubmp(nparb+2)
  real ( kind = 8 ) u0
  real ( kind = 8 ) u1
  real ( kind = 8 ) u2
  real ( kind = 8 ) ubc
  real ( kind = 8 ) v0
  real ( kind = 8 ) v1
  real ( kind = 8 ) v2
  real ( kind = 8 ) vbc
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) yc(np)
!
!  Determine the value of the "basis" shape function associated
!  with parameter IPARB, evaluated at node IP.
!
  if ( ishapb == 1 ) then
    call plval1 ( iparb+1, nparb+2, xc(ip), taubmp, shape )
  else if ( ishapb == 2 ) then
    call pqval1 ( iparb+1, nparb+2, xc(ip), taubmp, shape )
  else if ( ishapb == 3 ) then
    jderiv = 0
    call ppvalu ( taubmp, splbmp(1,1,iparb), nparb+1, 4, xc(ip), jderiv, shape )
  end if
!
!  Estimate dUdY and dVdY at the node.
!
  h1 = yc(ip+1)-yc(ip)
  h2 = yc(ip+2)-yc(ip)

  u0 = g(indx(ip,1))
  u1 = g(indx(ip+1,1))
  u2 = g(indx(ip+2,1))

  dudy = (h1**2*u2 - h2**2*u1 + (h2**2-h1**2)*u0)/(h1*h2*(h1-h2))

  v0 = g(indx(ip,2))
  v1 = g(indx(ip+1,2))
  v2 = g(indx(ip+2,2))

  dvdy = (h1**2*v2 - h2**2*v1 + (h2**2-h1**2)*v0)/(h1*h2*(h1-h2))
!
!  Set the boundary conditions.
!
  ubc = - dudy * shape
  vbc = - dvdy * shape

  return
end
subroutine bump_cost ( costb,ishapb,nparb,splbmp,taubmp,xbl,xbr,ybl,ybr)

!*****************************************************************************80
!
!! BUMP_COST evaluates the cost of the bump control.
!
!  Discussion:
!
!    The bump connects the points (XBL,YBL) and (XBR,YBR).
!
!    Compute its "cost" by comparing its slope to the slope of the
!    straight line that connects those two points.
!
!      COSTB = Integral (XBL <= X <= XBR) (Bump'(X) - Line'(X))**2 dX
!
!    Here, Bump(X) represents the function describing the shape
!    of the bump, and Line(X) represents the straight line which
!    simply joins the two endpoints, (XBL,YBL) and (XBR,YBR).
!
!    This integral is approximated by numerical integration.
!
!    The interval between XBL and XBR is divided into NPARB+1
!    intervals, over each of which the bump's height is described
!    by a cubic.
!
!    For each such interval, pick NQUAD1 quadrature points,
!    evaluate the derivative of the bump function there, and
!    subtract the slope of the straight line.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) COSTB, the integral of the difference of the
!    derivatives of the straight line joining the two straight line
!    line segments of the bottom, and the bump that is
!    actually drawn there.
!
!    This measures the cost of bump control.
!
  implicit none

  integer ( kind = 4 ) nparb
  integer ( kind = 4 ), parameter :: nquad1 = 5

  real ( kind = 8 ) costb
  real ( kind = 8 ) cprime
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ishapb
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jderiv
  real ( kind = 8 ) slope
  real ( kind = 8 ) splbmp(4,nparb+2,0:nparb)
  real ( kind = 8 ) taubmp(nparb+2)
  real ( kind = 8 ) wquad1(nquad1)
  real ( kind = 8 ) xbl
  real ( kind = 8 ) xbr
  real ( kind = 8 ) xleft
  real ( kind = 8 ) xquad1(nquad1)
  real ( kind = 8 ) xrite
  real ( kind = 8 ) xx
  real ( kind = 8 ) ybl
  real ( kind = 8 ) ybr
  real ( kind = 8 ) yvec(nparb+2)

  costb = 0.0D+00

  if ( nparb == 0 ) then
    return
  end if

  if ( xbl >= xbr ) then
    return
  end if
!
!  Get the Gauss weights and abscissas.
!
  call gquad1 ( nquad1, wquad1, xquad1 )
!
!  Get the slope of the line joining the endpoints of the bump.
!
  slope = (ybr-ybl)/(xbr-xbl)

  do i = 1,nparb+1

    xleft = ( real ( nparb - i + 2, kind = 8 ) * xbl &
            + real (         i - 1, kind = 8 ) * xbr ) &
            / real ( nparb +     1, kind = 8 )

    xrite = ( real ( nparb - i + 1, kind = 8 ) * xbl &
            + real (         i,     kind = 8 ) * xbr ) &
            / real ( nparb     + 1, kind = 8)

    do j = 1,nquad1

      xx = 0.5D+00*((1.0D+00+xquad1(j))*xrite+(1.0D+00-xquad1(j))*xleft)

      if ( ishapb == 1 ) then
        yvec(1:nparb+2) = splbmp(1,1:nparb+2,0)
        call pldx(nparb+2,xx,taubmp,cprime,yvec)
      else if ( ishapb == 2 ) then
        yvec(1:nparb+2) = splbmp(1,1:nparb+2,0)
        call pqdx(nparb+2,xx,taubmp,cprime,yvec)
      else if ( ishapb == 3 ) then
        jderiv = 1
        call ppvalu(taubmp,splbmp,nparb+1,4,xx,jderiv,cprime)
      end if

      costb = costb + 0.5D+00 * wquad1(j) * ( xrite - xleft ) &
        * ( cprime - slope )**2

    end do

  end do

  return
end
subroutine bump_der ( dpara, ishapb, npar, nparb, nparf, splbmp, taubmp, &
  wateb, xbl, xbr, ybl, ybr )

!*****************************************************************************80
!
!! BUMP_DER differentiates the bump control cost with respect to the parameters.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!
  implicit none

  integer ( kind = 4 ) npar
  integer ( kind = 4 ) nparb
  integer ( kind = 4 ), parameter :: nquad1 = 5

  real ( kind = 8 ) cprime
  real ( kind = 8 ) dpara(npar)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ishapb
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jderiv
  integer ( kind = 4 ) k
  integer ( kind = 4 ) nparf
  real ( kind = 8 ) slope
  real ( kind = 8 ) splbmp(4,nparb+2,0:nparb)
  real ( kind = 8 ) taubmp(nparb+2)
  real ( kind = 8 ) temp
  real ( kind = 8 ) wateb
  real ( kind = 8 ) wquad1(nquad1)
  real ( kind = 8 ) xbl
  real ( kind = 8 ) xbr
  real ( kind = 8 ) xleft
  real ( kind = 8 ) xquad1(nquad1)
  real ( kind = 8 ) xrite
  real ( kind = 8 ) xx
  real ( kind = 8 ) ybl
  real ( kind = 8 ) ybr
  real ( kind = 8 ) yvec(nparb+2)

  if ( nparb == 0 ) then
    return
  end if

  if ( xbr  <=  xbl )then
    return
  end if
!
!  Get the Gauss weights and abscissas.
!
  call gquad1(nquad1,wquad1,xquad1)
!
!  Get the slope of the straight line connecting the endpoints.
!
  slope = (ybr-ybl)/(xbr-xbl)

  do i = 1,nparb+1

    xleft = ( real ( nparb - i + 2, kind = 8 ) * xbl &
            + real (         i - 1, kind = 8 ) * xbr ) &
            / real ( nparb     + 1, kind = 8 )

    xrite = ( real ( nparb - i + 1, kind = 8 ) * xbl &
            + real (         i,     kind = 8 ) * xbr ) &
            / real ( nparb     + 1, kind = 8 )

    do j = 1,nquad1

      xx = 0.5D+00*((1.0D+00+xquad1(j))*xrite+(1.0D+00-xquad1(j))*xleft)

      if ( ishapb == 1 ) then
        yvec(1:nparb+2) = splbmp(1,1:nparb+2,0)
        call pldx(nparb+2,xx,taubmp,cprime,yvec)
      else if ( ishapb == 2 ) then
        yvec(1:nparb+2) = splbmp(1,1:nparb+2,0)
        call pqdx(nparb+2,xx,taubmp,cprime,yvec)
      else if ( ishapb == 3 ) then
        jderiv = 1
        call ppvalu(taubmp,splbmp,nparb+1,4,xx,jderiv,cprime)
      end if

      do k = 1,nparb

        if ( ishapb == 1 ) then
          call pldx1(k+1,nparb+2,xx,taubmp,temp)
        else if ( ishapb == 2 ) then
          call pqdx1(k+1,nparb+2,xx,taubmp,temp)
        else if ( ishapb == 3 ) then
          jderiv = 1
          call ppvalu(taubmp,splbmp(1,1,k),nparb+1,4,xx,jderiv,temp)
        end if

        dpara(nparf+k) = dpara(nparf+k)+wateb*wquad1(j)*(xrite-xleft)* &
          (cprime-slope)*temp

      end do

    end do

  end do

  return
end
subroutine bump_fx ( area,dudyn,dvdyn,eqn,g,ibc,indx,ipar,ishapb,nelem,neqn, &
  node,np,npar,nparb,nparf,npe,para,phi,res,sens,splbmp,taubmp,xc,yc)

!*****************************************************************************80
!
!! BUMP_FX measures the residual in the bump sensitivity equations.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) nelem
  integer ( kind = 4 ) neqn
  integer ( kind = 4 ) np
  integer ( kind = 4 ) npar
  integer ( kind = 4 ) nparf
  integer ( kind = 4 ) npe

  real ( kind = 8 ) ar
  real ( kind = 8 ) area(nelem,3)
  real ( kind = 8 ) dpdx
  real ( kind = 8 ) dpdy
  real ( kind = 8 ) dppdx
  real ( kind = 8 ) dppdy
  real ( kind = 8 ) dudx
  real ( kind = 8 ) dudy
  real ( kind = 8 ) dudyn(np)
  real ( kind = 8 ) dupdx
  real ( kind = 8 ) dupdy
  real ( kind = 8 ) dvdx
  real ( kind = 8 ) dvdy
  real ( kind = 8 ) dvdyn(np)
  real ( kind = 8 ) dvpdx
  real ( kind = 8 ) dvpdy
  real ( kind = 8 ) dwidx
  real ( kind = 8 ) dwidy
  character ( len = 2 ) eqn(neqn)
  real ( kind = 8 ) g(neqn)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibc
  integer ( kind = 4 ) ielem
  integer ( kind = 4 ) ihor
  integer ( kind = 4 ) indx(np,3)
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) ipar
  integer ( kind = 4 ) iparb
  integer ( kind = 4 ) iprs
  integer ( kind = 4 ) iq
  integer ( kind = 4 ) iquad
  integer ( kind = 4 ) ishapb
  integer ( kind = 4 ) iver
  integer ( kind = 4 ) node(nelem,npe)
  integer ( kind = 4 ) nparb
  real ( kind = 8 ) p
  real ( kind = 8 ) para(npar)
  real ( kind = 8 ) phi(nelem,3,6,10)
  real ( kind = 8 ) pp
  real ( kind = 8 ) qi
  real ( kind = 8 ) res(neqn)
  real ( kind = 8 ) nu_inv
  real ( kind = 8 ) sens(neqn,npar)
  real ( kind = 8 ) shape
  real ( kind = 8 ) splbmp(4,nparb+2,0:nparb)
  real ( kind = 8 ) taubmp(nparb+2)
  real ( kind = 8 ) u
  real ( kind = 8 ) ubc
  real ( kind = 8 ) up
  real ( kind = 8 ) v
  real ( kind = 8 ) vbc
  real ( kind = 8 ) vp
  real ( kind = 8 ) wi
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) yc(np)

  nu_inv = para(nparf+nparb+1)
  iparb = ipar-nparf

  res(1:neqn) = 0.0D+00
!
!  Approximate the integral by summing over all elements.
!
  do ielem = 1, nelem
!
!  Evaluate the integrand at the quadrature points.
!
    do iquad = 1,3

      ar = area(ielem,iquad)
!
!  For the given quadrature point, evaluate P, U, V.
!
      call uvalq ( dpdx,dpdy,dudx,dudy,dvdx,dvdy,g,ielem,indx,iquad,nelem,neqn, &
        node,np,npe,p,phi,u,v)

      call upvalq ( dppdx,dppdy,dupdx,dupdy,dvpdx,dvpdy,sens(1,ipar),ielem, &
        indx,iquad,nelem,neqn,node,np,npe,phi,pp,up,vp)
!
!  The only equations that have a contribution from this element
!  are those associated with basis functions for the element.
!  These, in turn, are associated with the nodes of the element.
!
!  So now we consider each node in the element.
!
      do iq = 1,6

        ip = node(ielem,iq)

        wi = phi(ielem,iquad,iq,1)
        dwidx = phi(ielem,iquad,iq,2)
        dwidy = phi(ielem,iquad,iq,3)
        qi = phi(ielem,iquad,iq,4)

        ihor = indx(ip,1)
        if ( eqn(ihor) == 'U' ) then
          res(ihor) = res(ihor)+ar*(dupdx*dwidx+dupdy*dwidy &
            +nu_inv*(up*dudx+u*dupdx+vp*dudy+v*dupdy+dppdx)*wi)
        else if ( eqn(ihor) == 'UB' ) then
          if ( ibc == 0 ) then
            call bump_bc(dudyn,dvdyn,ip,iparb,ishapb,np,nparb, &
              shape,splbmp,taubmp,ubc,vbc,xc)
          else if ( ibc == 1 ) then
            call bump_bc1(g,indx,ip,iparb,ishapb,neqn,np,nparb, &
              shape,splbmp,taubmp,ubc,vbc,xc,yc)
          else if ( ibc == 2 ) then
            call bump_bc2(g,indx,ip,iparb,ishapb,neqn,np,nparb, &
              shape,splbmp,taubmp,ubc,vbc,xc,yc)
          else if ( ibc == 3 ) then
            call bump_bc(dudyn,dvdyn,ip,iparb,ishapb,np,nparb, &
              shape,splbmp,taubmp,ubc,vbc,xc)
          end if
          res(ihor) = sens(ihor,ipar)-ubc
        else if ( eqn(ihor) == 'UI' ) then
          res(ihor) = sens(ihor,ipar)
        else if ( eqn(ihor) == 'UW' ) then
          res(ihor) = sens(ihor,ipar)
        end if

        iver = indx(ip,2)
        if ( eqn(iver) == 'V' ) then
          res(iver) = res(iver)+ar*(dvpdx*dwidx+dvpdy*dwidy &
            +nu_inv*(up*dvdx+u*dvpdx+vp*dvdy+v*dvpdy+dppdy)*wi)
        else if ( eqn(iver) == 'VB' ) then
          if ( ibc == 0 ) then
            call bump_bc(dudyn,dvdyn,ip,iparb,ishapb,np,nparb, &
              shape,splbmp,taubmp,ubc,vbc,xc)
          else if ( ibc == 1 ) then
            call bump_bc1(g,indx,ip,iparb,ishapb,neqn,np,nparb, &
              shape,splbmp,taubmp,ubc,vbc,xc,yc)
          else if ( ibc == 2 ) then
            call bump_bc2(g,indx,ip,iparb,ishapb,neqn,np,nparb, &
              shape,splbmp,taubmp,ubc,vbc,xc,yc)
          else if ( ibc == 3 ) then
            call bump_bc(dudyn,dvdyn,ip,iparb,ishapb,np,nparb, &
              shape,splbmp,taubmp,ubc,vbc,xc)
          end if
          res(iver) = sens(iver,ipar)-vbc
        else if ( eqn(iver) == 'VI' ) then
          res(iver) = sens(iver,ipar)
        else if ( eqn(iver) == 'VW' ) then
          res(iver) = sens(iver,ipar)
        end if

        iprs = indx(ip,3)
        if ( iprs>0 ) then
          if ( nu_inv == 0.0D+00 ) then
            res(iprs) = sens(iprs,ipar)
          else if ( eqn(iprs) == 'P' ) then
            res(iprs) = res(iprs)+ar*(dupdx+dvpdy)*qi
          else if ( eqn(iprs) == 'PB' ) then
            res(iprs) = sens(iprs,ipar)
          end if
        end if

      end do
    end do
  end do

  return
end
subroutine bump_sen ( dudyn,dvdyn,eqn,f,g,ibc,indx,ipar,ishapb,neqn,np,nparb, &
  nparf,splbmp,taubmp,xc,yc)

!*****************************************************************************80
!
!! BUMP_SEN sets up the right hand side F associated with the
!  sensitivities of a given flow solution (U,V,P) with respect to the
!  IPAR-th bump parameter.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) neqn
  integer ( kind = 4 ) np
  integer ( kind = 4 ) nparf

  real ( kind = 8 ) dudyn(np)
  real ( kind = 8 ) dvdyn(np)
  character ( len = 2 ) eqn(neqn)
  real ( kind = 8 ) f(neqn)
  real ( kind = 8 ) g(neqn)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibc
  integer ( kind = 4 ) ihor
  integer ( kind = 4 ) indx(np,3)
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) ipar
  integer ( kind = 4 ) iparb
  integer ( kind = 4 ) ishapb
  integer ( kind = 4 ) iver
  integer ( kind = 4 ) nparb
  real ( kind = 8 ) shape
  real ( kind = 8 ) splbmp(4,nparb+2,0:nparb)
  real ( kind = 8 ) taubmp(nparb+2)
  real ( kind = 8 ) ubc
  real ( kind = 8 ) vbc
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) yc(np)

  iparb = ipar-nparf

  f(1:neqn) = 0.0D+00

  do ip = 1,np

    ihor = indx(ip,1)
    if ( eqn(ihor) == 'UB' ) then
      if ( ibc == 0 ) then
        call bump_bc(dudyn,dvdyn,ip,iparb,ishapb,np,nparb,shape,splbmp,taubmp, &
          ubc,vbc,xc)
      else if ( ibc == 1 ) then
        call bump_bc1(g,indx,ip,iparb,ishapb,neqn,np,nparb,shape,splbmp,taubmp, &
          ubc,vbc,xc,yc)
      else if ( ibc == 2 ) then
        call bump_bc2(g,indx,ip,iparb,ishapb,neqn,np,nparb,shape,splbmp,taubmp, &
          ubc,vbc,xc,yc)
      else if ( ibc == 3 ) then
        call bump_bc(dudyn,dvdyn,ip,iparb,ishapb,np,nparb,shape,splbmp,taubmp, &
          ubc,vbc,xc)
      end if
      f(ihor) = ubc
    end if

    iver = indx(ip,2)
    if ( eqn(iver) == 'VB' ) then
      if ( ibc == 0 ) then
        call bump_bc(dudyn,dvdyn,ip,iparb,ishapb,np,nparb,shape,splbmp,taubmp, &
          ubc,vbc,xc)
      else if ( ibc == 1 ) then
        call bump_bc1(g,indx,ip,iparb,ishapb,neqn,np,nparb, &
          shape,splbmp,taubmp,ubc,vbc,xc,yc)
      else if ( ibc == 2 ) then
        call bump_bc2(g,indx,ip,iparb,ishapb,neqn,np,nparb, &
          shape,splbmp,taubmp,ubc,vbc,xc,yc)
      else if ( ibc == 3 ) then
        call bump_bc(dudyn,dvdyn,ip,iparb,ishapb,np,nparb, &
          shape,splbmp,taubmp,ubc,vbc,xc)
      end if
      f(iver) = vbc
    end if

  end do

  return
end
subroutine bump_spl ( ishapb, npar, nparb, nparf, par, splbmp, taubmp, xbl, &
  xbr, ybl, ybr )

!*****************************************************************************80
!
!! BUMP_SPL sets up or updates the spline data that describes the bump.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ISHAPB.
!    1, the bump is modeled by C0 linear splines.
!    2, the bump is modeled by C0 quadratic splines.
!    3, the bump is modeled by C1 cubic splines.
!
  implicit none

  integer ( kind = 4 ) npar
  integer ( kind = 4 ) nparb
  integer ( kind = 4 ) nparf

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibcbeg
  integer ( kind = 4 ) ibcend
  integer ( kind = 4 ) ipar
  integer ( kind = 4 ) ishapb
  real ( kind = 8 ) par(npar)
  real ( kind = 8 ) splbmp(4,nparb+2,0:nparb)
  real ( kind = 8 ) taubmp(nparb+2)
  real ( kind = 8 ) xbl
  real ( kind = 8 ) xbr
  real ( kind = 8 ) ybl
  real ( kind = 8 ) ybr

  if ( nparb <= 0 ) then
    return
  end if
!
!  Set up the bump arrays, including:
!
!    TAUBMP, containing the abscissas, which never change,
!    SPLBMP(1,I,0), the location of the bump at abscissa I,
!    SPLBMP(1,I,IPAR), the I-th coefficient of the partial derivative
!    of the bump shape with respect to the IPAR-th bump parameter.
!
  do i = 1, nparb+2
    taubmp(i) = ( real ( nparb + 2 - i ) * xbl + real ( i - 1 ) * xbr ) &
      / real ( nparb + 1 )
  end do

  do i = 1, nparb+2

    if ( i == 1 ) then
      splbmp(1,i,0) = ybl
    else if ( 2  <= i .and. i <= nparb+1 ) then
      splbmp(1,i,0) = par(nparf+i-1)
    else if ( i == nparb+2 ) then
      splbmp(1,i,0) = ybr
    end if

  end do

  if ( ishapb == 3 ) then

    do i = 1, nparb+2

      do ipar = 1, nparb

        if ( ipar+1 /=  i ) then
          splbmp(1,i,ipar) = 0.0D+00
        else
          splbmp(1,i,ipar) = 1.0D+00
        end if

      end do

    end do

    ibcbeg = 0
    ibcend = 0
    do i = 0, nparb
      call cubspl ( taubmp, splbmp(1,1,i), nparb+2, ibcbeg, ibcend )
    end do

  end if

  return
end
subroutine ch_cap ( c )

!*****************************************************************************80
!
!! CH_CAP capitalizes a single character.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character C, the character to capitalize.
!
  implicit none

  character c
  integer ( kind = 4 ) itemp

  itemp = ichar ( c )

  if ( 97 <= itemp .and. itemp <= 122 ) then
    c = char ( itemp - 32 )
  end if

  return
end
subroutine chkdat(ibump,idfd,ids,ifds,igrad,ijac,iopt,ipred,itype,jjac, &
  maxpar,maxparb,maxparf,nopt,npar,nparb,nparf,nstep3,para1,partar,xbleft, &
  xbrite)

!*****************************************************************************80
!
!! CHKDAT performs some simple checks on the input data.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IGRAD, the cost gradient approximation option.
!    0, No cost gradient approximation is made.
!    1, Chain rule on discretized sensitivities.
!    2, Chain rule on finite coefficient differences.
!    3, Chain rule on corrected finite coefficient differences.
!    4, Direct finite cost differences.
!
  implicit none

  integer ( kind = 4 ) maxpar

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibump
  integer ( kind = 4 ) idfd
  integer ( kind = 4 ) ids
  integer ( kind = 4 ) ifds
  integer ( kind = 4 ) igrad
  integer ( kind = 4 ) ijac
  integer ( kind = 4 ) iopt(maxpar)
  integer ( kind = 4 ) ipred
  integer ( kind = 4 ) isum
  integer ( kind = 4 ) itype
  integer ( kind = 4 ) jjac
  integer ( kind = 4 ) maxparb
  integer ( kind = 4 ) maxparf
  integer ( kind = 4 ) nopt
  integer ( kind = 4 ) npar
  integer ( kind = 4 ) nparb
  integer ( kind = 4 ) nparf
  integer ( kind = 4 ) nstep3
  real ( kind = 8 ) para1(maxpar)
  real ( kind = 8 ) partar(maxpar)
  real ( kind = 8 ) psum
  real ( kind = 8 ) xbleft
  real ( kind = 8 ) xbrite
!
!  IBUMP
!
  if ( ibump < 0 .or. ibump > 3 ) then
    write ( *, * ) ' '
    write ( *, * ) 'CHKDAT - Fatal error!'
    write ( *, * ) '  IBUMP is out of bounds!'
    write ( *, * ) '  Current value of IBUMP = ',ibump
    write ( *, * ) '  Legal values must be between 0 and 3.'
    stop
  end if
!
!  IGRAD
!
  if ( igrad < 0 .or. igrad > 4 ) then
    write ( *, * ) ' '
    write ( *, * ) 'CHKDAT - Fatal error!'
    write ( *, * ) '  IGRAD must be between 0 and 4,'
    write ( *, * ) '  but your value is ',igrad
    stop
  end if

  if ( igrad == 0 ) then
    if ( itype == 3 ) then
      write ( *, * ) ' '
      write ( *, * ) 'CHKDAT - Fatal error!'
      write ( *, * ) '  A cost gradient approximation MUST be made'
      write ( *, * ) '  when ITYPE = ',itype
      write ( *, * ) '  but you have chosen IGRAD = ',igrad
      stop
    end if
  else if ( igrad == 1 ) then
    if ( ids == 0 ) then
      write ( *, * ) ' '
      write ( *, * ) 'CHKDAT - Fatal error!'
      write ( *, * ) '  Your cost gradient choice IGRAD = ',igrad
      write ( *, * ) '  requires a nonzero value of IDS,'
      write ( *, * ) '  but your value is IDS = ',ids
      stop
    end if
  else if ( igrad == 2 ) then
    if ( ifds == 0 ) then
      write ( *, * ) ' '
      write ( *, * ) 'CHKDAT - Fatal error!'
      write ( *, * ) '  Your cost gradient choice IGRAD = ',igrad
      write ( *, * ) '  requires a nonzero value of IFDS,'
      write ( *, * ) '  but your value is IFDS = ',ifds
      stop
    end if
  else if ( igrad == 3 ) then
    if ( ifds == 0 ) then
      write ( *, * ) ' '
      write ( *, * ) 'CHKDAT - Fatal error!'
      write ( *, * ) '  Your cost gradient choice IGRAD = ',igrad
      write ( *, * ) '  requires a nonzero value of IFDS,'
      write ( *, * ) '  but your value is IFDS = ',ifds
      stop
    end if
  else if ( igrad == 4 ) then
    if ( idfd == 0 ) then
      write ( *, * ) ' '
      write ( *, * ) 'CHKDAT - Fatal error!'
      write ( *, * ) '  Your cost gradient choice IGRAD = ',igrad
      write ( *, * ) '  requires a nonzero value of IDFD,'
      write ( *, * ) '  but your value is IDFD = ',idfd
      stop
    end if
  end if
!
!  IJAC
!
  if ( ijac<1 ) then
    write ( *, * ) ' '
    write ( *, * ) 'CHKDAT - Fatal error!'
    write ( *, * ) '  IJAC is out of bounds!'
    write ( *, * ) '  Current value of IJAC = ',ijac
    write ( *, * ) '  Legal values must be 1 or greater.'
    stop
  end if
!
!  IOPT
!
  do i = 1,npar
    if ( iopt(i)/= 0 .and. iopt(i)/=1 ) then
      write ( *, * ) ' '
      write ( *, * ) 'CHKDAT - Fatal error!'
      write ( *, * ) '  IOPT(*) must be 0 or 1, but'
      write ( *, * ) '  IOPT(',I,') = ',iopt(i)
      stop
    end if
  end do
!
!  IPRED
!
  if ( ipred<0 .or. ipred>4 ) then
    write ( *, * ) ' '
    write ( *, * ) 'CHKDAT - Fatal error!'
    write ( *, * ) '  IPRED must be between 0 and 4,'
    write ( *, * ) '  but your value is ',ipred
    stop
  end if

  if ( ipred == 1 .and. ids==0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'CHKDAT - Fatal error!'
    write ( *, * ) '  IPRED = ',ipred
    write ( *, * ) '  but IDS = ',ids
    stop
  else if ( ipred == 2 .and. ifds==0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'CHKDAT - Fatal error!'
    write ( *, * ) '  IPRED = ',ipred
    write ( *, * ) '  but IFDS = ',ifds
    stop
  else if ( ipred == 3 .and. ifds==0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'CHKDAT - Fatal error!'
    write ( *, * ) '  IPRED = ',ipred
    write ( *, * ) '  but IFDS = ',ifds
    stop
  else if ( ipred == 4 .and. ids==0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'CHKDAT - Fatal error!'
    write ( *, * ) '  IPRED = ',ipred
    write ( *, * ) '  but IDS = ',ids
    stop
  end if
!
!  ITYPE
!
  if ( itype<1 .or. itype>7 ) then
    write ( *, * ) ' '
    write ( *, * ) 'CHKDAT - Fatal error!'
    write ( *, * ) '  Illegal value of ITYPE = ',itype
    write ( *, * ) '  Legal values are between 1 and 7.'
    write ( *, * ) '  ChkDat forces a STOP!'
    stop
  end if
!
!  JJAC
!
  if ( jjac<0 .or. jjac>3 ) then
    write ( *, * ) ' '
    write ( *, * ) 'CHKDAT - Fatal error!'
    write ( *, * ) '  Illegal value of JJAC = ',jjac
    write ( *, * ) '  Legal values are between 0 and 3.'
    write ( *, * ) '  ChkDat forces a STOP!'
    stop
  end if
!
!  NOPT
!
  if ( nopt <= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'CHKDAT - Fatal error!'
    write ( *, * ) '  The number of free parameters, NOPT = ',nopt
    write ( *, * ) '  but this value must be positive!'
    stop
  end if
!
!  NPAR
!
  if ( npar>maxpar ) then
    write ( *, * ) ' '
    write ( *, * ) 'CHKDAT - Fatal error!'
    write ( *, * ) '  NPAR is out of bounds!'
    write ( *, * ) '  Current value of NPAR = ',npar
    write ( *, * ) '  Maximum legal value, MAXPARA = ',maxpar
    stop
  end if
!
!  NPARB
!
  if ( nparb<0 .or. nparb>maxparb ) then
    write ( *, * ) ' '
    write ( *, * ) 'CHKDAT - Fatal error!'
    write ( *, * ) '  NPARB is out of bounds!'
    write ( *, * ) '  Input value of NPARB = ',nparb
    write ( *, * ) '  Maximum legal value, MAXPARB = ',maxparb
    stop
  end if
!
!  NPARF
!
  if ( nparf <= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'CHKDAT - Fatal error!'
    write ( *, * ) '  NPARF is out of bounds!'
    write ( *, * ) '  The input value of NPARF is ',nparf
    write ( *, * ) '  But NPARF must be at least 1.'
    stop
  end if

  if ( nparf>maxparf ) then
    write ( *, * ) ' '
    write ( *, * ) 'CHKDAT - Fatal error!'
    write ( *, * ) '  NPARF is too big!'
    write ( *, * ) '  Input value of NPARF = ',nparf
    write ( *, * ) '  Maximum legal value, MAXPARF = ',maxparf
    stop
  end if
!
!  NSTEP3
!
  if ( itype == 4 ) then
    if ( nstep3<1 ) then
      write ( *, * ) ' '
      write ( *, * ) 'CHKDAT - Fatal error!'
      write ( *, * ) ' Nonpositive value for NSTEP3 = ',nstep3
      stop
    else if ( nstep3 == 1 ) then
      write ( *, * ) 'CHKDAT - Warning!'
      write ( *, * ) '  NSTEP3 = 1 is an unusual value!'
    end if
  end if
!
!  PARA1(1:NPARF)
!
  if ( itype == 3 ) then

    psum = sum ( abs ( para1(1:nparf) ) )

    if ( psum == 0.0D+00 ) then

      isum = 0
      do i = 1,nparf
        isum = isum+abs(iopt(i))
      end do

      if ( isum == 0 ) then
        write ( *, * ) ' '
        write ( *, * ) 'CHKDAT - Fatal error!'
        write ( *, * ) '  You are optimizing, ITYPE = ',itype
        write ( *, * ) '  But all inflows are zero,'
        write ( *, * ) '  and no inflow may vary.'
        stop
      end if

    end if
  end if
!
!  PARA1(NPARF+NPARB+1), the value of NU_INV.
!
  if ( para1(nparf+nparb+1)<0.0D+00 ) then
    write ( *, * ) ' '
    write ( *, * ) 'CHKDAT - Warning!'
    write ( *, * ) '  NU_INV entry of PARA1 is negative.'
    write ( *, * ) '  This is being changed to one.'
    para1(nparf+nparb+1) = 1.0D+00
  end if
!
!  PARTAR(NPARF+NPARB+1), the value of NU_INV.
!
  if ( partar(nparf+nparb+1)<0.0D+00 ) then
    write ( *, * ) ' '
    write ( *, * ) 'CHKDAT - Warning!'
    write ( *, * ) '  NU_INV entry of PARTAR is negative.'
    write ( *, * ) '  This is being changed to one.'
    partar(nparf+nparb+1) = 1.0D+00
  end if
!
!  XBLEFT, XBRITE
!
  if ( nparb>0 ) then
    if ( xbleft >= xbrite ) then
      write ( *, * ) ' '
      write ( *, * ) 'CHKDAT - Fatal error!'
      write ( *, * ) '  XBLEFT >= XBRITE.'
      write ( *, * ) '  XBLEFT = ',xbleft
      write ( *, * ) '  XBRITE = ',xbrite
      stop
    end if
  end if

  return
end
subroutine chkopt(cost,costar,dparfd,dparfdc,dparsn,g,gtar,ids,ifds,neqn, &
  npar,para,partar)

!*****************************************************************************80
!
!! CHKOPT is called at the end of an optimization to check the results.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) neqn
  integer ( kind = 4 ) npar

  real ( kind = 8 ) cdist
  real ( kind = 8 ) cost
  real ( kind = 8 ) costar
  real ( kind = 8 ) cpfd
  real ( kind = 8 ) cpfdc
  real ( kind = 8 ) cpsn
  real ( kind = 8 ) dcfd
  real ( kind = 8 ) dcfdc
  real ( kind = 8 ) dcsn
  real ( kind = 8 ) dparfd(npar)
  real ( kind = 8 ) dparfdc(npar)
  real ( kind = 8 ) dparsn(npar)
  real ( kind = 8 ) g(neqn)
  real ( kind = 8 ) gdist
  real ( kind = 8 ) gtar(neqn)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ids
  integer ( kind = 4 ) ifds
  real ( kind = 8 ) para(npar)
  real ( kind = 8 ) partar(npar)
  real ( kind = 8 ) pdist

  cdist = cost-costar

  gdist = 0.0D+00
  do i = 1,neqn
    gdist = gdist+(gtar(i)-g(i))**2
  end do
  gdist = sqrt(gdist)

  pdist = 0.0D+00
  do i = 1,npar
    pdist = pdist+(partar(i)-para(i))**2
  end do
  pdist = sqrt(pdist)

  write ( *, * ) ' '
  write ( *, * ) 'ChkOpt:'
  write ( *, * ) '  L2 Distance from target solution =   ',gdist
  write ( *, * ) '  L2 Distance from target parameters = ',pdist
  write ( *, * ) '  Distance from target cost =          ',cdist
  write ( *, * ) ' '
  write ( *, * ) '  Estimated cost change if we moved to target:'
  write ( *, * ) ' '

  if ( ids/= 0 ) then

    dcsn = 0.0D+00
    do i = 1,npar
      dcsn = dcsn+dparsn(i)*(partar(i)-para(i))
    end do
    write ( *, * ) '  Sensitivities:               ',dcsn
    cpsn = 0.0D+00
    do i = 1,npar
      cpsn = cpsn+dparsn(i)**2
    end do
    cpsn = sqrt(cpsn)
    write ( *, * ) '  L2 norm of disc. sens. cost gradients:  ',cpsn
  end if

  if ( ifds/= 0 ) then

    dcfd = 0.0D+00
    do i = 1,npar
      dcfd = dcfd+dparfd(i)*(partar(i)-para(i))
    end do
    write ( *, * ) '  Finite differences:          ',dcfd

    cpfd = 0.0D+00
    do i = 1,npar
      cpfd = cpfd+dparfd(i)**2
    end do
    cpfd = sqrt(cpfd)
    write ( *, * ) '  L2 norm of fd cost gradients:           ',cpfd

  end if

  dcfdc = 0.0D+00
  do i = 1,npar
    dcfdc = dcfdc+dparfdc(i)*(partar(i)-para(i))
  end do
  write ( *, * ) '  Corrected finite differences:',dcfdc

  cpfdc = 0.0D+00
  do i = 1,npar
    cpfdc = cpfdc+dparfdc(i)**2
  end do
  cpfdc = sqrt(cpfdc)
  write ( *, * ) '  L2 norm of corrected fd cost gradients: ',cpfdc

  return
end
subroutine chrctd(string,dval,ierror,lchar)

!*****************************************************************************80
!
!! CHRCTD accepts a string of characters, and tries to extract a
!  real real number from the initial part of the
!  string.
!
!  CHRCTD will read as many characters as possible until it reaches
!  the end of the string, or encounters a character which cannot be
!  part of the number.
!
!  Legal input is:
!
!     1 blanks,
!     2 '+' or '-' sign,
!     3 integer part,
!     4 decimal point,
!     5 fraction part,
!     6 'E' or 'e' or 'D' or 'd', exponent marker,
!     7 exponent sign,
!     8 exponent integer part,
!     9 exponent decimal point,
!    10 exponent fraction part,
!    11 blanks,
!    12 final comma,
!
!  with most quantities optional.
!
!  Examples:
!
!    STRING            DVAL
!
!    '1'               1.0D+00
!    '     1   '       1.0D+00
!    '1A'              1.0D+00
!    '12,34,56'        12.0D+00
!    '  34 7'          34.0D+00
!    '-1E2ABCD'        -100.0D+00
!    '-1X2ABCD'        -1.0D+00
!    ' 2E-1'           0.2
!    '23.45'           23.45
!    '-4.2E+2'         -420.0D+00
!    '17d2'            1700.0D+00
!    '-14e-2'         -0.14
!    'e2'              100.0D+00
!    '-12.73e-9.23'   -12.73 * 10.0**(-9.23)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!  STRING Input, CHARACTER*(*) STRING, the string containing the
!         data to be read.  Reading will begin at position 1 and
!         terminate at the end of the string, or when no more
!         characters can be read to form a legal real.  Blanks,
!         commas, or other nonnumeric data will, in particular,
!         cause the conversion to halt.
!
!  DVAL   Output, real ( kind = 8 ) DVAL, the value that was read
!         from the string.
!
!  IERROR Output, integer ( kind = 4 ) IERROR, error flag.
!
!         0, no errors occurred.
!
!         1, 2, 6 or 7, the input number was garbled.  The
!         value of IERROR is the last type of input successfully
!         read.  For instance, 1 means initial blanks, 2 means
!         a plus or minus sign, and so on.
!
!  LCHAR  Output, integer ( kind = 4 ) LCHAR, the number of characters read from
!         STRING to form the number, including any terminating
!         characters such as a trailing comma or blanks.
!
  implicit none

  character chrtmp
  real ( kind = 8 ) dval
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ihave
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) iterm
  integer ( kind = 4 ) jbot
  integer ( kind = 4 ) jsgn
  integer ( kind = 4 ) jtop
  integer ( kind = 4 ) lchar
  logical s_eqi
  integer ( kind = 4 ) nchar
  integer ( kind = 4 ) ndig
  real ( kind = 8 ) rbot
  real ( kind = 8 ) rexp
  real ( kind = 8 ) rtop
  character ( len = * ) string

  nchar = len(string)

  ierror = 0
  dval = 0.0D+00
  lchar = -1
  isgn = 1
  rtop = 0.0D+00
  rbot = 1.0D+00
  jsgn = 1
  jtop = 0
  jbot = 1
  ihave = 1
  iterm = 0

10    continue

  lchar = lchar+1
  chrtmp = string(lchar+1:lchar+1)
!
!  Blank character.
!
  if ( chrtmp == ' ' ) then
    if ( ihave == 2 .or. ihave==6.or.ihave==7 ) then
      iterm = 1
    else if ( ihave>1 ) then
      ihave = 11
    end if
!
!  Comma
!
  else if ( chrtmp == ',' ) then
    if ( ihave/= 1 ) then
      iterm = 1
      ihave = 12
      lchar = lchar+1
    end if
!
!  Minus sign.
!
  else if ( chrtmp == '-' ) then
    if ( ihave == 1 ) then
      ihave = 2
      isgn = -1
    else if ( ihave == 6 ) then
      ihave = 7
      jsgn = -1
    else
      iterm = 1
    end if
!
!  Plus sign.
!
  else if ( chrtmp == '+' ) then
    if ( ihave == 1 ) then
      ihave = 2
    else if ( ihave == 6 ) then
      ihave = 7
    else
      iterm = 1
    end if
!
!  Decimal point.
!
  else if ( chrtmp == '.' ) then
    if ( ihave<4 ) then
      ihave = 4
    else if ( ihave >= 6 .and. ihave<=8 ) then
      ihave = 9
    else
      iterm = 1
    end if
!
!  Exponent marker.
!
  else if ( s_eqi(chrtmp,'e') .or. s_eqi(chrtmp,'d')  ) then
    if ( ihave<6 ) then
      ihave = 6
    else
      iterm = 1
    end if
!
!  Digit.
!
  else if ( ihave<11 .and.  lge(chrtmp,'0').and.lle(chrtmp,'9')  ) then

    if ( ihave <= 2 ) then
      ihave = 3
    else if ( ihave == 4 ) then
      ihave = 5
    else if ( ihave == 6 .or. ihave==7 ) then
      ihave = 8
    else if ( ihave == 9 ) then
      ihave = 10
    end if

    read(chrtmp,'(i1)')ndig

    if ( ihave == 3 ) then
      rtop = 10*rtop+ndig
    else if ( ihave == 5 ) then
      rtop = 10*rtop+ndig
      rbot = 10*rbot
    else if ( ihave == 8 ) then
      jtop = 10*jtop+ndig
    else if ( ihave == 10 ) then
      jtop = 10*jtop+ndig
      jbot = 10*jbot
    end if
!
!  Anything else is regarded as a terminator.
!
  else
    iterm = 1
  end if
!
!  If we haven't seen a terminator, and we haven't examined the
!  entire string, go get the next character.
!
  if ( iterm/= 1 .and. lchar+1<nchar)go to 10
!
!  If we haven't seen a terminator, and we have examined the
!  entire string, then we're done, and LCHAR is equal to NCHAR.
!
  if ( iterm/= 1 .and. lchar+1 == nchar)lchar=nchar
!
!  Number seems to have terminated.  Have we got a legal number?
!  Not if we terminated in states 1, 2, 6 or 7!
!
  if ( ihave == 1 .or. ihave==2.or.ihave==6.or.ihave==7 ) then
    ierror = ihave
    write ( *, * ) ' '
    write ( *, * ) 'CHRCTD - Fatal error!'
    write ( *, * ) '  Illegal or nonnumeric input!'
    return
  end if
!
!  Number seems OK.  Form it.
!
  if ( jtop == 0 ) then
    rexp = 1.0D+00
  else
    if ( jbot == 1 ) then
      rexp = 10.0D+00**(jsgn*jtop)
    else
      rexp = real ( jsgn * jtop, kind = 8 )
      rexp = rexp / real ( jbot, kind = 8 )
      rexp = 10.0D+00**rexp
    end if
  end if

  dval = real ( isgn, kind = 8 ) * rexp * rtop / rbot

  return
end
subroutine chrcti(string,intval,ierror,lchar)

!*****************************************************************************80
!
!! CHRCTI accepts a STRING of characters and reads an integer
!  from STRING into INTVAL.  The STRING must begin with an integer
!  but that may be followed by other information.
!
!  CHRCTI will read as many characters as possible until it reaches
!  the end of the STRING, or encounters a character which cannot be
!  part of the number.
!
!  Legal input is
!
!    blanks,
!    initial sign,
!    integer ( kind = 4 ) part,
!    blanks,
!    final comma,
!
!  with most quantities optional.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!  STRING Input, CHARACTER*(*) STRING, the string containing the
!         data to be read.  Reading will begin at position 1 and
!         terminate at the end of the string, or when no more
!         characters can be read to form a legal integer.  Blanks,
!         commas, or other nonnumeric data will, in particular,
!         cause the conversion to halt.
!
!         Sample results:
!
!         STRING            INTVAL
!
!         '1'               1
!         '     1   '       1
!         '1A'              1
!         '12,34,56'        12
!         '  34 7'          34
!         '-1E2ABCD'        -100
!         '-1X2ABCD'        -1
!         ' 2E-1'           0
!         '23.45'           23
!
!  INTVAL Output, integer ( kind = 4 ) INTVAL, the integer read from the string.
!
!  IERROR Output, integer ( kind = 4 ) IERROR, error flag.
!         0 if no errors,
!         Value of IHAVE when error occurred otherwise.
!
!  LCHAR  Output, integer ( kind = 4 ) LCHAR, number of characters read from
!         STRING to form the number.
!
  implicit none

  character chrtmp
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ihave
  integer ( kind = 4 ) intval
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) iterm
  integer ( kind = 4 ) itop
  integer ( kind = 4 ) lchar
  integer ( kind = 4 ) nchar
  integer ( kind = 4 ) ndig
  character ( len = * ) string

  nchar = len(string)

  ierror = 0
  intval = 0
  lchar = -1
  isgn = 1
  itop = 0
  ihave = 1
  iterm = 0

10    continue

  lchar = lchar+1
  chrtmp = string(lchar+1:lchar+1)

  if ( chrtmp == ' ' ) then
    if ( ihave == 2 ) then
      iterm = 1
    else if ( ihave == 3 ) then
      ihave = 11
    end if
  else if ( chrtmp == ',' ) then
    if ( ihave/= 1 ) then
      iterm = 1
      ihave = 12
      lchar = lchar+1
    end if
  else if ( chrtmp == '-' ) then
    if ( ihave == 1 ) then
      ihave = 2
      isgn = -1
    else
      iterm = 1
    end if
  else if ( chrtmp == '+' ) then
    if ( ihave == 1 ) then
      ihave = 2
    else
      iterm = 1
    end if
  else if ( lge(chrtmp,'0') .and. lle(chrtmp,'9').and.ihave<11 ) then
    ihave = 3
    read(chrtmp,'(i1)')ndig
    itop = 10*itop+ndig
  else
    iterm = 1
  end if

  if ( iterm/= 1 .and. lchar+1<nchar)go to 10
  if ( iterm/= 1 .and. lchar+1 == nchar)lchar=nchar
!
!  Number seems to have terminated.  Have we got a legal number?
!
  if ( ihave == 1 .or. ihave==2 ) then
    ierror = ihave
    write ( *, * ) ' '
    write ( *, * ) 'CHRCTI - Fatal error!'
    write ( *, * ) '  IERROR = ',ierror
    write ( *, * ) '  Illegal or nonnumeric input:'
    write(*,'(1x,a)')string
    return
  end if
!
!  Number seems OK.  Form it.
!
  intval = isgn*itop
  return
end
subroutine cost_gradient ( dpara, g, gtar, indx, ishapb, neqn, np, npar, &
  nparb, nparf, nprof, ny, sens, splbmp, taubmp, wateb, watep, wateu, watev, &
  xbl, xbr, ybl, ybr, yc )

!*****************************************************************************80
!
!! COST_GRADIENT returns the gradient of the cost functional.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) DPARA(NPAR), the partial derivative of the cost with respect
!    to each parameter, computed using sensitivities or finite differences.
!
  implicit none

  integer ( kind = 4 ) neqn
  integer ( kind = 4 ) np
  integer ( kind = 4 ) npar
  integer ( kind = 4 ) nparb
  integer ( kind = 4 ) nparf
  integer ( kind = 4 ) ny

  real ( kind = 8 ) dpara(npar)
  real ( kind = 8 ) g(neqn)
  real ( kind = 8 ) gtar(neqn)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(np,3)
  integer ( kind = 4 ) ishapb
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) nprof(2*ny-1)
  real ( kind = 8 ) sens(neqn,npar)
  real ( kind = 8 ) splbmp(4,nparb+2,0:nparb)
  real ( kind = 8 ) taubmp(nparb+2)
  real ( kind = 8 ) wateb
  real ( kind = 8 ) watep
  real ( kind = 8 ) wateu
  real ( kind = 8 ) watev
  real ( kind = 8 ) xbl
  real ( kind = 8 ) xbr
  real ( kind = 8 ) ybl
  real ( kind = 8 ) ybr
  real ( kind = 8 ) yc(np)

  ival = 1
  call imemry('inc','COST_GRADIENT_calls',ival)

  dpara(1:npar) = 0.0D+00

  call bump_der ( dpara,ishapb,npar,nparb,nparf,splbmp,taubmp,wateb, &
    xbl,xbr,ybl,ybr)

  call disc_der ( dpara,g,gtar,indx,neqn,np,npar,nprof,ny,sens,watep, &
    wateu,watev,yc)

  return
end
subroutine cubspl(tau,c,n,ibcbeg,ibcend)

!*****************************************************************************80
!
!! CUBSPL is given data and boundary conditions for a cubic
!  spline, and computes information that defines that cubic
!  spline.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    Original FORTRAN77 version by Carl DeBoor.
!    FORTRAN90 version by John Burkardt
!
!  Parameters:
!
!  TAU    Input, real ( kind = 8 ) TAU(N), the abscissas or X values of
!         the data points.  The entries of TAU are assumed to be
!         strictly increasing.
!
!  N      Input, integer ( kind = 4 ) N, the number of data points.  N is
!         assumed to be at least 2.
!
!  C      Input/output, real ( kind = 8 ) C(4,N).
!
!         On input, if IBCBEG or IBCBEG is 1 or 2, then C(2,1)
!         or C(2,N) should have been set to the desired derivative
!         values, as described further under IBCBEG and IBCEND.
!
!         On output, C contains the polynomial coefficients of
!         the cubic interpolating spline with interior knots
!         TAU(2) through TAU(N-1).
!
!         In the interval interval (TAU(I), TAU(I+1)), the spline
!         F is given by
!
!           F(X) = C(1,I)+H*C(2,I)+(1/2)*H*H*C(3,I)
!                  +(1/6)*H*H*H*C(4,I)
!
!         where H = X - TAU(I).  The routine PPVALU may be used to
!         evaluate F or its derivatives from TAU, C, L = N-1,
!         and K = 4.
!
!  IBCBEG,
!  IBCEND Input, integer ( kind = 4 ) IBCBEG, IBCEND, boundary condition
!         indicators.
!
!         IBCBEG = 0 means no boundary condition at TAU(1) is given.
!         In this case, the "not-a-knot condition" is used.  That
!         is, the jump in the third derivative across TAU(2) is
!         forced to zero.  Thus the first and the second cubic
!         polynomial pieces are made to coincide.
!
!         IBCBEG = 1 means that the slope at TAU(1) is to equal the
!         input value C(2,1).
!
!         IBCBEG = 2 means that the second derivative at TAU(1) is
!         to equal C(2,1).
!
!         IBCEND = 0, 1, or 2 has analogous meaning concerning the
!         boundary condition at TAU(N), with the additional
!         information taken from C(2,N).
!
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c(4,n)
  real ( kind = 8 ) divdf1
  real ( kind = 8 ) divdf3
  real ( kind = 8 ) dtau
  real ( kind = 8 ) g
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibcbeg
  integer ( kind = 4 ) ibcend
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) j
  integer ( kind = 4 ) m
  real ( kind = 8 ) tau(n)

  ival = 1
  call imemry('inc','CubSpl_calls',ival)
!
!  A tridiagonal linear system for the unknown slopes S(I) of
!  F at TAU(I), I = 1,..., N, is generated and then solved by Gauss
!  elimination, with S(I) ending up in C(2,I), for all I.
!
!  C(3,*) and C(4,*) are used initially for temporary storage.
!
!  Store first differences of the TAU sequence in C(3,*).
!
!  Store first divided difference of data in C(4,*).
!
  do m = 2,n
    c(3,m) = tau(m)-tau(m-1)
    c(4,m) = (c(1,m)-c(1,m-1))/c(3,m)
  end do
!
!  Construct the first equation from the boundary condition, of
!  the form:
!
!    C(4,1)*S(1) + C(3,1)*S(2) = C(2,1)
!
  if ( ibcbeg == 1 ) then
    c(4,1) = 1.0D+00
    c(3,1) = 0.0D+00
    go to 60
  end if

  if ( ibcbeg <= 1 ) then
!
!  No condition at left end and N = 2.
!
    if ( n <= 2 ) then
      c(4,1) = 1.0D+00
      c(3,1) = 1.0D+00
      c(2,1) = 2.0D+00 * c(4,2)
      go to 120
    end if
!
!  Not-a-knot condition at left end and N is greater than 2.
!
    c(4,1) = c(3,3)
    c(3,1) = c(3,2)+c(3,3)
    c(2,1) = ((c(3,2)+2.0D+00*c(3,1))*c(4,2)*c(3,3)+c(3,2)**2*c(4,3)) /c(3,1)
    go to 70

  end if
!
!  Second derivative prescribed at left end.
!
  c(4,1) = 2.0D+00
  c(3,1) = 1.0D+00
  c(2,1) = 3.0D+00 * c(4,2)-c(3,2)/2.0D+00 * c(2,1)

60    continue

  if ( n == 2)go to 120
!
!  If there are interior knots, generate the corresponding
!  equations and carry out the forward pass of Gauss elimination,
!  after which the M-th equation reads:
!
!    C(4,M)*S(M) + C(3,M)*S(M+1) = C(2,M).
!
70    continue

  do m = 2,n-1
    g = -c(3,m+1)/c(4,m-1)
    c(2,m) = g*c(2,m-1)+3.0D+00*(c(3,m)*c(4,m+1)+c(3,m+1)*c(4,m))
    c(4,m) = g*c(3,m-1)+2.0D+00*(c(3,m)+c(3,m+1))
  end do
!
!  Construct last equation from the second boundary condition, of
!  the form
!
!    (-G*C(4,N-1))*S(N-1) + C(4,N)*S(N) = C(2,N)
!
!  If slope is prescribed at right end, one can go directly to
!  back-substitution, since the C array happens to be set up just
!  right for it at this point.
!
  if ( ibcend == 1)go to 160
  if ( ibcend>1)go to 110
!
!  Not-a-knot and N >= 3, and either N>3 or also not-a-knot
!  at left end point.
!
  if ( n/= 3 .or. ibcbeg/=0 ) then
    g = c(3,n-1)+c(3,n)
    c(2,n) = ((c(3,n)+2.0D+00*g)*c(4,n)*c(3,n-1)+c(3,n)**2*(c(1,n-1)-c(1,n-2))/c(3,n-1))/g
    g = -g/c(4,n-1)
    c(4,n) = c(3,n-1)
    c(4,n) = g*c(3,n-1)+c(4,n)
    c(2,n) = (g*c(2,n-1)+c(2,n))/c(4,n)
    go to 160
  end if
!
!  Either (N = 3 and not-a-knot also at left) or (N=2 and not not-a-
!  knot at left end point).
!
100   continue

  c(2,n) = 2.0D+00*c(4,n)
  c(4,n) = 1.0D+00
  g = -1.0D+00/c(4,n-1)
  c(4,n) = g*c(3,n-1)+c(4,n)
  c(2,n) = (g*c(2,n-1)+c(2,n))/c(4,n)
  go to 160
!
!  Second derivative prescribed at right endpoint.
!
110   continue

  c(2,n) = 3.0D+00*c(4,n)+c(3,n)/2.0D+00*c(2,n)
  c(4,n) = 2.0D+00
  g = -1.0D+00/c(4,n-1)
  c(4,n) = g*c(3,n-1)+c(4,n)
  c(2,n) = (g*c(2,n-1)+c(2,n))/c(4,n)
  go to 160

120   continue

  if ( ibcend/= 1 ) then

    if ( ibcend>1 ) then
      c(2,n) = 3.0D+00*c(4,n)+c(3,n)/2.0D+00*c(2,n)
      c(4,n) = 2.0D+00
      g = -1.0D+00/c(4,n-1)
      c(4,n) = g*c(3,n-1)+c(4,n)
      c(2,n) = (g*c(2,n-1)+c(2,n))/c(4,n)
      go to 160
    end if

    if ( ibcbeg>0)go to 100
!
!  Not-a-knot at right endpoint and at left endpoint and N = 2.
!
    c(2,n) = c(4,n)

  end if
!
!  Carry out the back substitution
!
160   continue

  do j = n-1,1,-1
    c(2,j) = (c(2,j)-c(3,j)*c(2,j+1))/c(4,j)
  end do
!
!  Generate cubic coefficients in each interval, that is, the
!  derivatives at its left endpoint, from value and slope at its
!  endpoints.
!
  do i = 2,n
    dtau = c(3,i)
    divdf1 = (c(1,i)-c(1,i-1))/dtau
    divdf3 = c(2,i-1)+c(2,i)-2.0D+00*divdf1
    c(3,i-1) = 2.0D+00*(divdf1-c(2,i-1)-divdf3)/dtau
    c(4,i-1) = (divdf3/dtau)*(6.0D+00/dtau)
  end do

  return
end
subroutine disc_cost ( costp, costu, costv, g, gtar, indx, neqn, np, nprof, &
  ny, yc )

!*****************************************************************************80
!
!! DISC_COST computes the discrepancy cost integrals.
!
!  Discussion:
!
!    The discrepancy cost integrals measure the discrepancies between the
!    target and computed pressure, horizontal and vertical velocities along
!    the profile line.
!
!    This integration scheme assumes that the profile line, and
!    the element sides that define it, are straight.  Otherwise,
!    the integration scheme used is not correct.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ), parameter :: nquad1 = 5

  integer ( kind = 4 ) neqn
  integer ( kind = 4 ) np
  integer ( kind = 4 ) ny

  real ( kind = 8 ) bval
  real ( kind = 8 ) costp
  real ( kind = 8 ) costu
  real ( kind = 8 ) costv
  real ( kind = 8 ) g(neqn)
  real ( kind = 8 ) gtar(neqn)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) indx(np,3)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) npol
  integer ( kind = 4 ) nprof(2*ny-1)
  real ( kind = 8 ) pcof(2)
  real ( kind = 8 ) pval
  real ( kind = 8 ) ucof(3)
  real ( kind = 8 ) uval
  real ( kind = 8 ) vcof(3)
  real ( kind = 8 ) vval
  real ( kind = 8 ) wquad1(nquad1)
  real ( kind = 8 ) xquad1(nquad1)
  real ( kind = 8 ) yc(np)
  real ( kind = 8 ) yhi
  real ( kind = 8 ) ylo
  real ( kind = 8 ) ypol(3)
  real ( kind = 8 ) yval
!
!  Get the weights and abscissas to approximate a line integral.
!
  call gquad1 ( nquad1, wquad1, xquad1 )
!
!  Compute the integral of the difference squared between the
!  current velocity and the target values.
!
  costu = 0.0D+00
  costv = 0.0D+00
!
!  The line along which we integrate is broken into NY-1
!  subintervals, over each of which, U and V are represented
!  by quadratic functions.
!
  do i = 1, ny-1
!
!  Get the values of U and V at the beginning, middle, and
!  end of the subinterval.  Use these to compute the quadratic
!  representation of U and V for any point on the subinterval.
!
    ylo = yc(nprof(2*i-1))
    yhi = yc(nprof(2*i+1))

    npol = 3

    do k = 1, npol

      ii = 2 * i - 2 + k
      ypol(k) = yc(nprof(ii))

      j = indx(nprof(ii),1)
      ucof(k) = g(j) - gtar(j)

      j = indx(nprof(ii),2)
      vcof(k) = g(j) - gtar(j)

    end do
!
!  Evaluate the discrepancy at each quadrature point.
!
    do j = 1, nquad1

      call rint_to_rint ( -1.0D+00, +1.0D+00, xquad1(j), ylo, yhi, yval )

      uval = 0.0D+00
      vval = 0.0D+00

      do k = 1, npol
        call lbase ( k, npol, bval, ypol, yval )
        uval = uval + bval * ucof(k)
        vval = vval + bval * vcof(k)
      end do

      costu = costu + 0.5D+00 * wquad1(j) * ( yhi - ylo ) * uval**2
      costv = costv + 0.5D+00 * wquad1(j) * ( yhi - ylo ) * vval**2

    end do
  end do
!
!  Compute the square root of the integral of the difference
!  squared between the current pressure and the target values.
!
  costp = 0.0D+00

  do i = 1, ny-1

    ylo = yc(nprof(2*i-1))
    yhi = yc(nprof(2*i+1))

    npol = 2

    do k = 1, npol

      ii = 2*i-3+2*k

      ypol(k) = yc(nprof(ii))

      j = indx(nprof(ii),3)
      if ( j  <= 0 ) then
        pcof(k) = 0.0D+00
      else
        pcof(k) = g(j)-gtar(j)
      end if

    end do

    do j = 1, nquad1

      call rint_to_rint ( -1.0D+00, +1.0D+00, xquad1(j), ylo, yhi, yval )

      pval = 0.0D+00

      do k = 1, npol
        call lbase ( k, npol, bval, ypol, yval )
        pval = pval + bval * pcof(k)
      end do

      costp = costp + 0.5D+00 * wquad1(j) * ( yhi - ylo ) * pval**2

    end do
  end do

  return
end
subroutine disc_der ( dpara, g, gtar, indx, neqn, np, npar, nprof, ny, sens, &
  watep, wateu, watev, yc )

!*****************************************************************************80
!
!! DISC_DER computes the derivative of the discrepancy cost integrals
!
!  Discussion:
!
!    The discrepancy integrals for pressure, horizontal and vertical
!    velocities are to be differentiated with respect to the
!    various parameters.
!
!    In practice, the dummy parameter DPARA will actually be DPARA1 or
!    DPARA2, and the dummy input parameter SENS will actually be SENS or
!    GDIF, depending on which estimate of the cost derivatives is
!    desired.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) DPARA(NPAR), contains the derivatives of 
!    the cost function with respect to the various parameters.  In particular,
!    DPARA(I) = D cost / D parameter(I).
!
!    Input, real ( kind = 8 ) G(NEQN).
!    G is the current solution vector, in which are stored 
!    pressures and velocities.
!
!    Input, real ( kind = 8 ) GTAR(NEQN),
!    the target solution vector.
!
!    Input, integer ( kind = 4 ) INDX(NP,3).  
!    INDX contains, for each node I, the index of U, V and P at 
!    that node, or 0 or a negative value.
!    If K=INDX(I,J) is positive, then the value of the degree
!    of freedom is stored in the solution vector entry G(K).
!    If INDX(I,J) is positive, then that means that a degree of
!    freedom for variable J (U, V or P) is associated with node
!    I, and an equation will be generated to determine its value.
!    If INDX(I,J) is zero, then that means the the value of variable
!    J (U, V or P) has been specified at node I.  No equation is
!    generated to determine its value.
!
!    Input, integer ( kind = 4 ) NEQN, the number of equations, and functions.
!
!    Input, integer ( kind = 4 ) NP, the number of nodes.  NP=(2*NX-1)*(2*NY-1).
!
!    Input, integer ( kind = 4 ) NPAR.
!    The number of parameters.  NPAR = NPARF + NPARB + 1.
!    The parameters control the shape of the inflow,
!    the shape of the bump obstacle, and the strength of the
!    flow.
!
!    Input, integer ( kind = 4 ) NPROF(2*MAXNY-1), the indices of the nodes along 
!    the profile line.
!
!    Input, integer ( kind = 4 ) NY.
!    NY controls the spacing of nodes and elements in
!    the Y direction.  There are 2*NY-1 nodes along various
!    lines in the Y direction.
!    Roughly speaking, NY is the number of elements along
!    a line in the Y direction.
!
!    Input, real ( kind = 8 ) SENS(MAXEQN,NPAR), the sensitivities.
!    SENS(I,J) contains the sensitivity of the I-th unknown
!    with respect to the J-th parameter.
!
!    Input, real ( kind = 8 ) WATEP, WATEU, WATEV.
!    These are weights used in computing the overall cost 
!    function based on the costs of the flow discrepancy.
!
!    Input, real ( kind = 8 ) YC(NP), the Y coordinates of the nodes
!
  implicit none

  integer ( kind = 4 ) mpar
  integer ( kind = 4 ) nquad1

  parameter (mpar = 10)
  parameter (nquad1 = 5)

  integer ( kind = 4 ) neqn
  integer ( kind = 4 ) np
  integer ( kind = 4 ) npar
  integer ( kind = 4 ) ny

  real ( kind = 8 ) bval
  real ( kind = 8 ) dpara(npar)
  real ( kind = 8 ) dpcof(mpar,mpar)
  real ( kind = 8 ) dpval(mpar)
  real ( kind = 8 ) ducof(mpar,mpar)
  real ( kind = 8 ) duval(mpar)
  real ( kind = 8 ) dvcof(mpar,mpar)
  real ( kind = 8 ) dvval(mpar)
  real ( kind = 8 ) g(neqn)
  real ( kind = 8 ) gtar(neqn)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) indx(np,3)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) npol
  integer ( kind = 4 ) nprof(2*ny-1)
  real ( kind = 8 ) pcof(2)
  real ( kind = 8 ) pval
  real ( kind = 8 ) sens(neqn,npar)
  real ( kind = 8 ) ucof(3)
  real ( kind = 8 ) uval
  real ( kind = 8 ) vcof(3)
  real ( kind = 8 ) vval
  real ( kind = 8 ) watep
  real ( kind = 8 ) wateu
  real ( kind = 8 ) watev
  real ( kind = 8 ) wquad1(nquad1)
  real ( kind = 8 ) xquad1(nquad1)
  real ( kind = 8 ) yc(np)
  real ( kind = 8 ) yhi
  real ( kind = 8 ) ylo
  real ( kind = 8 ) ypol(3)
  real ( kind = 8 ) yval

  if ( npar > mpar ) then
    write ( *, * ) ' '
    write ( *, * ) 'DISC_DER - Fatal error!'
    write ( *, * ) '  The number of parameters is too high.'
    write ( *, * ) '  This routine can handle NPAR = ',mpar
    write ( *, * ) '  Your problem has NPAR = ',npar
    stop
  end if
!
!  Get the Gauss weights and abscissas to approximate a line integral.
!
  call gquad1 ( nquad1, wquad1, xquad1 )
!
!  The line along which we integrate is broken into NY-1
!  subintervals, over each of which, U and V are represented
!  by quadratic functions.
!
  do i = 1, ny-1
!
!  Get the values of U and V at the beginning, middle, and
!  end of the subinterval.  Use these to compute the quadratic
!  representation of U and V for any point on the subinterval.
!
    ylo = yc(nprof(2*i-1))
    yhi = yc(nprof(2*i+1))

    npol = 3

    do k = 1,npol

      ii = 2*i-2+k

      ypol(k) = yc(nprof(ii))

      j = indx(nprof(ii),1)
      ucof(k) = g(j)-gtar(j)

      do l = 1,npar
        ducof(k,l) = sens(j,l)
      end do

      j = indx(nprof(ii),2)
      vcof(k) = g(j)-gtar(j)

      do l = 1,npar
        dvcof(k,l) = sens(j,l)
      end do

    end do
!
!  Evaluate the discrepancy at each quadrature point.
!
    do j = 1,nquad1

      call rint_to_rint ( -1.0D+00, +1.0D+00, xquad1(j), ylo, yhi, yval )

      uval = 0.0D+00
      vval = 0.0D+00
      do l = 1,npar
        duval(l) = 0.0D+00
        dvval(l) = 0.0D+00
      end do

      do k = 1,npol
        call lbase(k,npol,bval,ypol,yval)
        uval = uval+bval*ucof(k)
        vval = vval+bval*vcof(k)
        do l = 1,npar
          duval(l) = duval(l)+bval*ducof(k,l)
          dvval(l) = dvval(l)+bval*dvcof(k,l)
        end do
      end do

      do l = 1,npar
        dpara(l) = dpara(l)+0.5D+00*wquad1(j)*(yhi-ylo) &
          *2.0D+00*(wateu*uval*duval(l)+watev*vval*dvval(l))
      end do

    end do
  end do
!
!  Compute the square root of the integral of the difference
!  squared between the current pressure and the target values.
!
  do i = 1,ny-1

    ylo = yc(nprof(2*i-1))
    yhi = yc(nprof(2*i+1))

    npol = 2

    do k = 1,npol

      ii = 2*i-3+2*k

      ypol(k) = yc(nprof(ii))

      j = indx(nprof(ii),3)
      if ( j <= 0 ) then
        pcof(k) = 0.0D+00
      else
        pcof(k) = g(j)-gtar(j)
      end if

      do l = 1,npar
        if ( j <= 0 ) then
          dpcof(k,l) = 0.0D+00
        else
          dpcof(k,l) = sens(j,l)
        end if
      end do

    end do

    do j = 1,nquad1

      call rint_to_rint ( -1.0D+00, +1.0D+00, xquad1(j), ylo, yhi, yval )

      pval = 0.0D+00
      do l = 1,npar
        dpval(l) = 0.0D+00
      end do

      do k = 1,npol
        call lbase(k,npol,bval,ypol,yval)
        pval = pval+bval*pcof(k)
        do l = 1,npar
          dpval(l) = dpval(l)+bval*dpcof(k,l)
        end do

      end do

      do l = 1,npar
        dpara(l) = dpara(l) + 0.5D+00 * wquad1(j) * ( yhi - ylo ) &
          * 2.0D+00 * watep * pval * dpval(l)
      end do

    end do
  end do

  return
end
subroutine dmemry(action,name,dval)

!*****************************************************************************80
!
!! DMEMRY allows the user to define the name of a real
!  variable, set it, increment it, or get the value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!  ACTION Input, Character*(*) ACTION, desired action.
!
!         'Init', reset all values to zero, wipe out all names.
!         'Name', add a variable of the given name.
!         'Inc',  increment variable NAME by DVAL.
!         'Set',  set variable NAME to DVAL.
!         'Get',  return value of NAME in DVAL.
!         'Zero', reset all values to zero.
!
!  NAME   Input, Character*(*) NAME, the name of the variable.
!
!  DVAL   Input/output, real ( kind = 8 ) DVAL.
!
!         For the 'Inc' and 'Set' commands, DVAL must contain the
!         increment or set value.
!
!         For the 'Get' command, DVAL will contain the value of the
!         named variable on output.
!
  implicit none

  integer ( kind = 4 ) maxnam
  parameter (maxnam = 100)

  character ( len = * ) action
  real ( kind = 8 ) dval
  real ( kind = 8 ) dvals(maxnam)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) lenc
  logical s_eqi
  character ( len = * ) name
  character ( len = 20 ) names(maxnam)
  integer ( kind = 4 ) numnam

  save dvals
  save names
  save numnam
!
!  Initialize everything.
!
  if ( s_eqi(action,'init') ) then

    numnam = 0

    do i = 1,maxnam
      dvals(i) = 0.0D+00
      names(i) = ' '
    end do
!
!  Name something.
!
  else if ( s_eqi(action,'name') ) then

    do i = 1,numnam

      if ( s_eqi(name,names(i)) ) then
        write ( *, * ) ' '
        write ( *, * ) 'DMemry - Warning!'
        write(*,'(''  There is ALREADY a variable '',a)')name
        return
      end if

    end do

    if ( numnam<maxnam ) then
      numnam = numnam+1
      names(numnam) = name
    else
    write ( *, * ) ' '
      write ( *, * ) 'DMemry - Fatal error!'
      write ( *, * ) '  No more name space.'
      stop
    end if
!
!  Increment something.
!
  else if ( s_eqi(action,'inc') ) then

    do i = 1,numnam

      if ( s_eqi(name,names(i)) ) then
        dvals(i) = dvals(i)+dval
        return
      end if

    end do

    write ( *, * ) ' '
    write ( *, * ) 'DMemry - Fatal error!'
    write ( *, * ) '  Attempt to increment unknown variable.'
    write(*,'(''  Variable name is '',a)')name
    stop
!
!  Set something.
!
  else if ( s_eqi(action,'set') ) then

    do i = 1,numnam

      if ( s_eqi(name,names(i)) ) then
        dvals(i) = dval
        return
      end if

    end do

    write ( *, * ) ' '
    write ( *, * ) 'DMemry - Fatal error!'
    write ( *, * ) '  Attempt to set unknown variable.'
    write(*,'(''  Variable name is '',a)')name
    stop
!
!  Get something.
!
  else if ( s_eqi(action,'get') ) then

    do i = 1,numnam

      if ( s_eqi(name,names(i)) ) then
        dval = dvals(i)
        return
      end if

    end do

    write ( *, * ) ' '
    write ( *, * ) 'DMemry - Fatal error!'
    write ( *, * ) '  Attempt to get value of unknown variable.'
    write(*,'(''  Variable name is '',a)')name
    stop
!
!  Initialize everything.
!
  else if ( s_eqi(action,'zero') ) then

    do i = 1,numnam
      dvals(i) = 0.0D+00
    end do

    write ( *, * ) ' '
    write ( *, * ) 'DMemry - Note:'
    write ( *, * ) '  All data has been reset to zero.'
!
!  Unknown action.
!
  else

    write ( *, * ) ' '
    write ( *, * ) 'DMemry - Fatal error!'
    write ( *, * ) '  Unrecognized action requested.'
    write ( *, '(1x,a)')action
    stop

  end if

  return
end
subroutine flo_spl_set ( ishapf, npar, nparf, par, splflo, tauflo )

!*****************************************************************************80
!
!! FLO_SPL_SET sets up or updates the spline data that describes the inflow.
!
!  Discussion:
!
!    It does this for the target parameters and the feasible parameters.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ISHAPF, determines the type of spline used.
!    1, Inflow modeled by C0 piecewise linears.
!    2, Inflow modeled by C0 piecewise quadratics.
!    3, Inflow modeled by C1 cubic splines.
!
  implicit none

  integer ( kind = 4 ) npar
  integer ( kind = 4 ) nparf

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibcbeg
  integer ( kind = 4 ) ibcend
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ipar
  integer ( kind = 4 ) ishapf
  real ( kind = 8 ) par(npar)
  real ( kind = 8 ) splflo(4,nparf+2,0:nparf)
  real ( kind = 8 ) tauflo(nparf+2)
!
!  The abscissas TAUFLO are evenly spaced.
!
  call r4vec_even ( 0.0D+00, 3.0D+00, nparf+2, tauflo )
!
!  The data values SPLFLO are the parameters, padded with a zero at each end.
!
  splflo(1,1,0) = 0.0D+00
  splflo(1,2:nparf+1,0) = par(1:nparf)
  splflo(1,nparf+2,0) = 0.0D+00
!
!  Set up SPLFLO(1,I,IPAR), the I-th coefficient of the partial derivative
!  of the inflow with respect to the IPAR-th inflow parameter.
!
  if ( ishapf == 3 ) then

    do i = 1,nparf+2

      do ipar = 1,nparf
        if ( ipar+1/= i ) then
          splflo(1,i,ipar) = 0.0D+00
        else
          splflo(1,i,ipar) = 1.0D+00
        end if
      end do

    end do

    ibcbeg = 0
    ibcend = 0
    do i = 0,nparf
      call cubspl(tauflo,splflo(1,1,i),nparf+2,ibcbeg,ibcend)
    end do

  end if

  return
end
subroutine flo_spl_val ( ishapf, nparf, splflo, tauflo, ubc, vbc, yy )

!***************************************************************************
!
!! FLO_SPL_VAL computes boundary velocities specified by a parameterized spline.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ISHAPF, determines the type of spline used.
!    1, Inflow modeled by C0 piecewise linears.
!    2, Inflow modeled by C0 piecewise quadratics.
!    3, Inflow modeled by C1 cubic splines.
!
  implicit none

  integer ( kind = 4 ) nparf

  integer ( kind = 4 ) ishapf
  integer ( kind = 4 ) jderiv
  real ( kind = 8 ) splflo(4,nparf+2,0:nparf)
  real ( kind = 8 ) tauflo(nparf+2)
  real ( kind = 8 ) ubc
  real ( kind = 8 ) vbc
  real ( kind = 8 ) yvec(nparf+2)
  real ( kind = 8 ) yy

  if ( ishapf == 1 ) then
    yvec(1:nparf+2) = splflo(1,1:nparf+2,0)
    call plval( nparf+2, yy, tauflo, ubc, yvec )
  else if ( ishapf == 2 ) then
    yvec(1:nparf+2) = splflo(1,1:nparf+2,0)
    call pqval ( nparf+2, yy, tauflo, ubc, yvec )
  else if ( ishapf == 3 ) then
    jderiv = 0
    call ppvalu ( tauflo, splflo(1,1,0), nparf+1, 4, yy, jderiv, ubc )
  else
    write ( *, * ) ' '
    write ( *, * ) 'FLO_SPL_VAL - Fatal error!'
    write ( *, * ) '  Illegal value of ISHAPF = ', ishapf
    stop
  end if

  vbc = 0.0D+00

  return
end
subroutine floduv ( ipar, ishapf, nparf, splflo, tauflo, ubc, vbc, yy )

!***************************************************************************
!
!! FLODUV differentiates the parabolic inflow with respect to parameter.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) UBC, the partial derivative of the horizontal
!    boundary velocity component with respect to a given
!    parameter at the specified point.
!
!    Output, real ( kind = 8 ) VBC, the partial derivative of the vertical
!    boundary velocity component with respect to a given
!    parameter at the specified point.
!
  implicit none

  integer ( kind = 4 ) nparf

  integer ( kind = 4 ) ipar
  integer ( kind = 4 ) ishapf
  integer ( kind = 4 ) jderiv
  real ( kind = 8 ) splflo(4,nparf+2,0:nparf)
  real ( kind = 8 ) tauflo(nparf+2)
  real ( kind = 8 ) ubc
  real ( kind = 8 ) vbc
  real ( kind = 8 ) yy

  ubc = 0.0D+00
  vbc = 0.0D+00
!
!  Handle the case where there is only one inflow parameter, and it
!  is specified as a constant.
!
  if ( nparf <= 0 .and. ipar == 0 ) then
    ubc = 4.0D+00 * ( 3.0D+00 - yy ) * yy / 9.0D+00
    return
  end if

  if ( ipar<1 .or. nparf < ipar ) then
    return
  end if
!
!  Evaluate the basis function which is zero at the other parameters
!  and 1 at parameter IPAR.
!
  if ( ishapf == 1 ) then
    call plval1 ( ipar+1, nparf+2, yy, tauflo, ubc )
  else if ( ishapf == 2 ) then
    call pqval1 ( ipar+1, nparf+2, yy, tauflo, ubc )
  else if ( ishapf == 3 ) then
    jderiv = 0
    call ppvalu ( tauflo, splflo(1,1,ipar), nparf+1, 4, yy, jderiv, ubc )
  else
    write ( *, * ) ' '
    write ( *, * ) 'FLODUV - Fatal error!'
    write ( *, * ) '  Illegal value of ISHAPF = ',ishapf
    stop
  end if

  return
end
subroutine flosen ( eqn,f,indx,ipar,ishapf,neqn,np,nparf,splflo,tauflo,yc)

!*****************************************************************************80
!
!! FLOSEN sets up the right hand side for a sensitivity parameter.
!
!  Discussion:
!
!    This routine sets up the right hand side F associated with the
!    sensitivities of a given flow solution (U,V,P) with
!    respect to the IPAR-th inflow parameter.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) F(NEQN), the right hand
!    side of the sensitivity equations associated with
!    the IPAR-th inflow parameter.
!
  implicit none

  integer ( kind = 4 ) neqn
  integer ( kind = 4 ) np
  integer ( kind = 4 ) nparf

  character ( len = 2 ) eqn(neqn)
  real ( kind = 8 ) f(neqn)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihor
  integer ( kind = 4 ) indx(np,3)
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) ipar
  integer ( kind = 4 ) ishapf
  integer ( kind = 4 ) iver
  real ( kind = 8 ) splflo(4,nparf+2,0:nparf)
  real ( kind = 8 ) tauflo(nparf+2)
  real ( kind = 8 ) ubc
  real ( kind = 8 ) vbc
  real ( kind = 8 ) yc(np)

  f(1:neqn) = 0.0D+00

  do ip = 1,np

    ihor = indx(ip,1)
    if ( eqn(ihor) == 'UI' ) then
      call floduv(ipar,ishapf,nparf,splflo,tauflo,ubc,vbc,yc(ip))
      f(ihor) = ubc
    end if

    iver = indx(ip,2)
    if ( eqn(iver) == 'VI' ) then
      call floduv(ipar,ishapf,nparf,splflo,tauflo,ubc,vbc,yc(ip))
      f(iver) = vbc
    end if

  end do

  return
end
subroutine flosol ( a,area,disjac,eqn,etan,etaq,flarea,g,ierror,igrid,ijac, &
  indx,ipivot,ishapb,ishapf,isotri,iwrite,jjac,maxnew,nelem,neqn,nlband,node, &
  np,npar,nparb,nparf,npe,nrow,nx,ny,para,parjac,phi,res,splbmp,splflo,syseqn, &
  taubmp,tauflo,tolnew,wquad,xbl,xbord,xbr,xc,xquad,xsin,xsiq,ybl,ybord,ybr, &
  yc,yquad)

!*****************************************************************************80
!
!! FLOSOL solves the flow equations for given parameter values.
!
!  Discussion:
!
!    This routine is given a set of flow parameters in PARA, and an
!    approximate solution vector G, and proceeds to set up the
!    constraints associated with PARA, and use Newton iteration
!    to correct G to a solution that satisfies the constraints
!    to within some tolerance.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) G(NEQN).
!    G is the computed solution vector, in which are stored
!    pressures and velocities.
!
!    Output, integer ( kind = 4 ) IERROR.
!    Error flag.  0 means no error, nonzero means an error.
!
  implicit none

  integer ( kind = 4 ) nelem
  integer ( kind = 4 ) neqn
  integer ( kind = 4 ) np
  integer ( kind = 4 ) npar
  integer ( kind = 4 ) nparb
  integer ( kind = 4 ) nparf
  integer ( kind = 4 ) npe
  integer ( kind = 4 ) nrow
  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny

  real ( kind = 8 ) a(nrow,neqn)
  real ( kind = 8 ) area(nelem,3)
  real ( kind = 8 ) disjac
  character ( len = 2 ) eqn(neqn)
  real ( kind = 8 ) etan(6)
  real ( kind = 8 ) etaq(3)
  real ( kind = 8 ) flarea
  real ( kind = 8 ) g(neqn)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) igrid
  integer ( kind = 4 ) ijac
  integer ( kind = 4 ) indx(np,3)
  integer ( kind = 4 ) ipivot(neqn)
  integer ( kind = 4 ) ishapb
  integer ( kind = 4 ) ishapf
  integer ( kind = 4 ) isotri(nelem)
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) iwrite
  integer ( kind = 4 ) jjac
  integer ( kind = 4 ) maxnew
  integer ( kind = 4 ) nlband
  integer ( kind = 4 ) node(nelem,npe)
  real ( kind = 8 ) para(npar)
  real ( kind = 8 ) parjac(npar)
  real ( kind = 8 ) phi(nelem,3,6,10)
  real ( kind = 8 ) res(neqn)
  real ( kind = 8 ) splbmp(4,nparb+2,0:nparb)
  real ( kind = 8 ) splflo(4,nparf+2,0:nparf)
  character ( len = 20 ) syseqn
  real ( kind = 8 ) taubmp(nparb+2)
  real ( kind = 8 ) tauflo(nparf+2)
  real ( kind = 8 ) tolnew
  real ( kind = 8 ) wquad(3)
  real ( kind = 8 ) xbl
  real ( kind = 8 ) xbord(nx)
  real ( kind = 8 ) xbr
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) xquad(nelem,3)
  real ( kind = 8 ) xsin(6)
  real ( kind = 8 ) xsiq(3)
  real ( kind = 8 ) ybl
  real ( kind = 8 ) ybord(ny)
  real ( kind = 8 ) ybr
  real ( kind = 8 ) yc(np)
  real ( kind = 8 ) yquad(nelem,3)

  ival = 1
  call imemry('inc','FloSol_calls',ival)
!
!  Set the spline coefficients for the bump.
!
  call bump_spl ( ishapb,npar,nparb,nparf,para,splbmp,taubmp,xbl,xbr,ybl,ybr)
!
!  Set the spline coefficients for the inflow.
!
  call flo_spl_set ( ishapf, npar, nparf, para, splflo, tauflo )
!
!  Set the X and Y coordinates of the nodes that form the grid.
!
  call xy_set ( igrid, ishapb, np, nparb, nx, ny, splbmp, taubmp, xbl, &
    xbord, xbr, xc, ybl, ybord, ybr, yc )
!
!  Set the quadrature points, which move every step if there
!  are bump parameters.
!
  call setqxy ( area,etan,etaq,isotri,nelem,node,np,npe,wquad,xc,xquad,xsin, &
    xsiq,yc,yquad)
!
!  Set the value of the basis functions at all quadrature points.
!
  call setbas ( area,etaq,flarea,isotri,nelem,node,np,npe,phi,xc,xquad, &
    xsiq,yc,yquad)
!
!  Solve the nonlinear system.
!
  call newton ( a,area,disjac,eqn,g,ierror,ijac,indx,ipivot,ishapf,iwrite, &
    jjac,maxnew,nelem,neqn,nlband,node,np,npar,nparb,nparf,npe,nrow,para, &
    parjac,phi,res,splflo,syseqn,tauflo,tolnew,yc)

  if ( ierror/= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'FLOSOL - Fatal error!'
    write ( *, * ) '  Newton failed!'
    write ( *, * ) '  The parameters at which failure occurred:'
    write ( *, * ) ' '
    call pr_parameter(nparb,nparf,para)
    ierror = 1
    return
  end if

  return
end
subroutine fprime ( a,area,eqn,g,indx,nelem,neqn,nlband,node,np,npar,nparb, &
  nparf,npe,nrow,para,phi,syseqn)

!*****************************************************************************80
!
!! FPRIME computes the jacobian of the Navier Stokes or Stokes functions.
!
!  Discussion:
!
!  The differentiated Navier Stokes functions have the form:
!
!  d U-Eqn/d U-Coef:
!
!    Integral
!
!      dWj/dx * dWi/dx + dWj/dy * dWi/dy
!    + nu_inv * (Wj*dUold/dx + Uold*dWj/dx+ Vold*dWj/dy) * Wi dx dy
!
!  d U-Eqn/d V-Coef:
!
!    Integral
!
!    nu_inv * Wj*dUold/dy * Wi dx dy
!
!  d U-Eqn/d P-Coef:
!
!    Integral
!
!    nu_inv * dQj/dx * Wi dx dy
!
!  d V-Eqn/d U-Coef:
!
!    Integral
!
!    nu_inv * Wj*dVold/dx * Wi dx dy
!
!  d V-Eqn/d V-Coef:
!
!    Integral
!
!      dWj/dx * dWi/dx + dWj/dy * dWi/dy
!    + nu_inv * (Uold*dWj/dx + Wj*dVold/dy + Vold*dWj/dy) * Wi dx dy
!
!  d V-Eqn/d P-Coef:
!
!    Integral
!
!    nu_inv * dQj/dy * Wi dx dy
!
!  d P-Eqn/d U-Coef:
!
!    Integral
!
!      dWj/dx * Qi dx dy
!
!    Integral
!
!      dWj/dy * Qi dx dy
!
!
!  The differentiated Stokes functions have the form:
!
!
!  d U-Eqn/d U-Coef:
!
!    Integral
!
!      dWj/dx * dWi/dx + dWj/dy * dWi/dy dx dy
!
!  d U-Eqn/d P-Coef:
!
!    Integral
!
!    nu_inv * dQj/dx * Wi dx dy
!
!  d V-Eqn/d V-Coef:
!
!    Integral
!
!      dWj/dx * dWi/dx + dWj/dy * dWi/dy dx dy
!
!  d V-Eqn/d P-Coef:
!
!    Integral
!
!    nu_inv * dQj/dy * Wi dx dy
!
!  d P-Eqn/d U-Coef:
!
!    Integral
!
!      dWj/dx * Qi dx dy
!
!    Integral
!
!      dWj/dy * Qi dx dy
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A(NROW,NEQN), the value of D F(I)/D X(J) for each of the
!    NEQN residual functions F(I) with respect to each of the unknown
!    coefficients X(J).
!
  implicit none

  integer ( kind = 4 ) nelem
  integer ( kind = 4 ) neqn
  integer ( kind = 4 ) np
  integer ( kind = 4 ) npar
  integer ( kind = 4 ) nparf
  integer ( kind = 4 ) npe
  integer ( kind = 4 ) nrow

  real ( kind = 8 ) a(nrow,neqn)
  real ( kind = 8 ) ar
  real ( kind = 8 ) area(nelem,3)
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
  character ( len = 2 ) eqn(neqn)
  real ( kind = 8 ) g(neqn)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ielem
  integer ( kind = 4 ) ihor
  integer ( kind = 4 ) indx(np,3)
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) iprs
  integer ( kind = 4 ) iq
  integer ( kind = 4 ) iquad
  integer ( kind = 4 ) iuse
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) iver
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhor
  integer ( kind = 4 ) jp
  integer ( kind = 4 ) jprs
  integer ( kind = 4 ) jq
  integer ( kind = 4 ) jver
  logical s_eqi
  integer ( kind = 4 ) nlband
  integer ( kind = 4 ) node(nelem,npe)
  integer ( kind = 4 ) nparb
  real ( kind = 8 ) pold
  real ( kind = 8 ) para(npar)
  real ( kind = 8 ) phi(nelem,3,6,10)
  real ( kind = 8 ) qi
  real ( kind = 8 ) nu_inv
  character ( len = * ) syseqn
  real ( kind = 8 ) term
  real ( kind = 8 ) uold
  real ( kind = 8 ) vold
  real ( kind = 8 ) wi
  real ( kind = 8 ) wj

  ival = 1
  call imemry('inc','Fprime_calls',ival)

  nu_inv = para(nparf+nparb+1)

  a(1:nrow,1:neqn) = 0.0D+00
!
!  Approximate the integral by summing over all elements.
!
  do ielem = 1,nelem
!
!  Evaluate the integrand at the quadrature points.
!
    do iquad = 1,3

      ar = area(ielem,iquad)
!
!  For the given quadrature point, evaluate P, U and V.
!
      call uvalq ( dpdx,dpdy,dudx,dudy,dvdx,dvdy,g,ielem,indx,iquad,nelem,neqn, &
        node,np,npe,pold,phi,uold,vold)
!
!  Consider each node in the element.
!
      do iq = 1,6

        ip = node(ielem,iq)

        wi = phi(ielem,iquad,iq,1)
        dwidx = phi(ielem,iquad,iq,2)
        dwidy = phi(ielem,iquad,iq,3)
        qi = phi(ielem,iquad,iq,4)

        ihor = indx(ip,1)
        iver = indx(ip,2)
        iprs = indx(ip,3)
!
!  Now compute the derivatives of the functions associated
!  with U, V and P, with respect to the coefficients associated
!  with basis vectors at each node of the element.
!
        do jq = 1,6

          jp = node(ielem,jq)

          wj = phi(ielem,iquad,jq,1)
          dwjdx = phi(ielem,iquad,jq,2)
          dwjdy = phi(ielem,iquad,jq,3)

          dqjdx = phi(ielem,iquad,jq,5)
          dqjdy = phi(ielem,iquad,jq,6)

          jhor = indx(jp,1)
          jver = indx(jp,2)
          jprs = indx(jp,3)
!
!  Contributions of the JHOR horizontal velocity to the U, V, and
!  P equations.
!
          iuse = ihor-jhor+2*nlband+1

          if ( eqn(ihor) == 'U' ) then

            if ( s_eqi(syseqn,'NavierStokes') ) then
              term = ar*(dwjdx*dwidx+dwjdy*dwidy+ &
                nu_inv*(wj*dudx+uold*dwjdx+vold*dwjdy)*wi)
            else if ( s_eqi(syseqn,'Stokes') ) then
              term = ar*(dwjdx*dwidx+dwjdy*dwidy)
            end if

            a(iuse,jhor) = a(iuse,jhor)+term

          end if

          if ( eqn(iver) == 'V' ) then
            if ( s_eqi(syseqn,'NavierStokes') ) then
              iuse = iver-jhor+2*nlband+1
              term = ar*(nu_inv*wj*dvdx*wi)
              a(iuse,jhor) = a(iuse,jhor)+term
            end if
          end if

          if ( iprs>0 ) then
            if ( eqn(iprs) == 'P' ) then
              iuse = iprs-jhor+2*nlband+1
              term = ar*dwjdx*qi
              a(iuse,jhor) = a(iuse,jhor)+term
            end if
          end if
!
!  Contributions of the JVER vertical velocity variable to the
!  U, V and P equations.
!
          if ( eqn(ihor) == 'U' ) then
            if ( s_eqi(syseqn,'NavierStokes') ) then
              iuse = ihor-jver+2*nlband+1
              term = ar*nu_inv*wj*dudy*wi
              a(iuse,jver) = a(iuse,jver)+term
            end if
          end if

          iuse = iver-jver+2*nlband+1
          if ( eqn(iver) == 'V' ) then
            if ( s_eqi(syseqn,'NavierStokes') ) then
              term = ar*(dwjdx*dwidx+dwjdy*dwidy &
                +nu_inv*(uold*dwjdx+wj*dvdy+vold*dwjdy)*wi)
            else if ( s_eqi(syseqn,'Stokes') ) then
              term = ar*(dwjdx*dwidx+dwjdy*dwidy)
            end if
            a(iuse,jver) = a(iuse,jver)+term
          end if

          if ( iprs>0 ) then
            if ( eqn(iprs) == 'P' ) then
              iuse = iprs-jver+2*nlband+1
              term = ar*dwjdy*qi
              a(iuse,jver) = a(iuse,jver)+term
            end if
          end if
!
!  Contributions of the JPRS pressure to the U and V equations.
!
          if ( jprs>0 ) then

            if ( eqn(ihor) == 'U' ) then
              iuse = ihor-jprs+2*nlband+1
              term = ar*nu_inv*dqjdx*wi
              a(iuse,jprs) = a(iuse,jprs)+term
            end if

            if ( eqn(iver) == 'V' ) then
              iuse = iver-jprs+2*nlband+1
              term = ar*nu_inv*dqjdy*wi
              a(iuse,jprs) = a(iuse,jprs)+term
            end if

          end if

        end do
      end do
    end do
  end do
!
!  Set up the equations that enforce boundary conditions.
!
  do ip = 1,np

    ihor = indx(ip,1)
    iver = indx(ip,2)
    iprs = indx(ip,3)

    if ( eqn(ihor) == 'UB' .or. eqn(ihor)=='UI'.or.eqn(ihor)=='UW' ) then
      a(2*nlband+1,ihor) = 1.0D+00
    end if

    if ( eqn(iver) == 'VB' .or. eqn(iver)=='VI'.or.eqn(iver)=='VW' ) then
      a(2*nlband+1,iver) = 1.0D+00
    end if

    if ( iprs>0 ) then
      if ( eqn(iprs) == 'PB' ) then
        a(2*nlband+1,iprs) = 1.0D+00
      end if
    end if

  end do

  return
end
subroutine fprnab ( anab,area,eqn,g,indx,nelem,neqn,node,np,npar,nparb,nparf, &
  npe,ny,para,phi,syseqn)

!*****************************************************************************80
!
!! FPRNAB computes the jacobian and stores it in a neighbor matrix.
!
!  Discussion:
!
!    This routine computes the jacobian of the Navier Stokes or Stokes 
!    residual functions and stores it in the neighbor matrix ANAB.
!
!  The differentiated Navier Stokes functions have the form:
!
!
!  d U-Eqn/d U-Coef:
!
!    Integral
!
!      dWj/dx * dWi/dx + dWj/dy * dWi/dy
!    + nu_inv * (Wj*dUold/dx + Uold*dWj/dx+ Vold*dWj/dy) * Wi dx dy
!
!  d U-Eqn/d V-Coef:
!
!    Integral
!
!    nu_inv * Wj*dUold/dy * Wi dx dy
!
!  d U-Eqn/d P-Coef:
!
!    Integral
!
!    nu_inv * dQj/dx * Wi dx dy
!
!  d V-Eqn/d U-Coef:
!
!    Integral
!
!    nu_inv * Wj*dVold/dx * Wi dx dy
!
!  d V-Eqn/d V-Coef:
!
!    Integral
!
!      dWj/dx * dWi/dx + dWj/dy * dWi/dy
!    + nu_inv * (Uold*dWj/dx + Wj*dVold/dy + Vold*dWj/dy) * Wi dx dy
!
!  d V-Eqn/d P-Coef:
!
!    Integral
!
!    nu_inv * dQj/dy * Wi dx dy
!
!  d P-Eqn/d U-Coef:
!
!    Integral
!
!      dWj/dx * Qi dx dy
!
!    Integral
!
!      dWj/dy * Qi dx dy
!
!
!  The differentiated Stokes functions have the form:
!
!
!  d U-Eqn/d U-Coef:
!
!    Integral
!
!      dWj/dx * dWi/dx + dWj/dy * dWi/dy dx dy
!
!  d U-Eqn/d P-Coef:
!
!    Integral
!
!    nu_inv * dQj/dx * Wi dx dy
!
!  d V-Eqn/d V-Coef:
!
!    Integral
!
!      dWj/dx * dWi/dx + dWj/dy * dWi/dy dx dy
!
!  d V-Eqn/d P-Coef:
!
!    Integral
!
!    nu_inv * dQj/dy * Wi dx dy
!
!  d P-Eqn/d U-Coef:
!
!    Integral
!
!      dWj/dx * Qi dx dy
!
!    Integral
!
!      dWj/dy * Qi dx dy
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) ANAB(3,3,19,NP), contains the value of D F(I)/D X(J) 
!    for each of the NEQN residual functions F(I) with respect to each 
!    of the unknown coefficients X(J), stored in a neighbor format.
!
  implicit none

  integer ( kind = 4 ) nelem
  integer ( kind = 4 ) neqn
  integer ( kind = 4 ) np
  integer ( kind = 4 ) npar
  integer ( kind = 4 ) nparf
  integer ( kind = 4 ) npe

  real ( kind = 8 ) anab(3,3,19,np)
  real ( kind = 8 ) ar
  real ( kind = 8 ) area(nelem,3)
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
  character ( len = 2 ) eqn(neqn)
  real ( kind = 8 ) g(neqn)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ielem
  integer ( kind = 4 ) ihor
  integer ( kind = 4 ) indx(np,3)
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) iprs
  integer ( kind = 4 ) iq
  integer ( kind = 4 ) iquad
  integer ( kind = 4 ) iuse
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) iver
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jcol
  integer ( kind = 4 ) jp
  integer ( kind = 4 ) jq
  integer ( kind = 4 ) jrow
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  logical s_eqi
  integer ( kind = 4 ) my
  integer ( kind = 4 ) node(nelem,npe)
  integer ( kind = 4 ) nparb
  integer ( kind = 4 ) ny
  real ( kind = 8 ) pold
  real ( kind = 8 ) para(npar)
  real ( kind = 8 ) phi(nelem,3,6,10)
  real ( kind = 8 ) qi
  real ( kind = 8 ) nu_inv
  character ( len = * ) syseqn
  real ( kind = 8 ) term
  real ( kind = 8 ) uold
  real ( kind = 8 ) vold
  real ( kind = 8 ) wi
  real ( kind = 8 ) wj

  external s_eqi

  ival = 1
  call imemry('inc','Fprime_calls',ival)

  nu_inv = para(nparf+nparb+1)

  anab(1:3,1:3,1:19,1:np) = 0.0D+00
!
!  Approximate the integral by summing over all elements.
!
  do ielem = 1,nelem
!
!  Evaluate the integrand at the quadrature points.
!
    do iquad = 1,3

      ar = area(ielem,iquad)
!
!  For the given quadrature point, evaluate P, U and V.
!
      call uvalq ( dpdx,dpdy,dudx,dudy,dvdx,dvdy,g,ielem,indx, &
        iquad,nelem,neqn,node,np,npe,pold,phi,uold,vold)
!
!  Consider each node in the element.
!
      do iq = 1,6

        ip = node(ielem,iq)

        wi = phi(ielem,iquad,iq,1)
        dwidx = phi(ielem,iquad,iq,2)
        dwidy = phi(ielem,iquad,iq,3)
        qi = phi(ielem,iquad,iq,4)

        ihor = indx(ip,1)
        iver = indx(ip,2)
        iprs = indx(ip,3)
!
!  Now compute the derivatives of the functions associated
!  with U, V and P, with respect to the coefficients associated
!  with basis vectors at each node of the element.
!
        do jq = 1,6

          jp = node(ielem,jq)

          jcol = ((jp-1)/(2*ny-1))+1
          jrow = jp-(jcol-1)*(2*ny-1)

          wj = phi(ielem,iquad,jq,1)
          dwjdx = phi(ielem,iquad,jq,2)
          dwjdy = phi(ielem,iquad,jq,3)

          dqjdx = phi(ielem,iquad,jq,5)
          dqjdy = phi(ielem,iquad,jq,6)

          call nodnab(ip,jp,my,iuse)
!
!  Contributions of the JHOR variable to various equations.
!
          if ( eqn(ihor) == 'U' ) then

            if ( s_eqi(syseqn,'NavierStokes') ) then
              term = ar*(dwjdx*dwidx+dwjdy*dwidy+ &
                nu_inv*(wj*dudx+uold*dwjdx+vold*dwjdy)*wi)
            else if ( s_eqi(syseqn,'Stokes') ) then
              term = ar*(dwjdx*dwidx+dwjdy*dwidy)
            end if

            anab(1,1,iuse,ip) = anab(1,1,iuse,ip)+term

          end if

          if ( eqn(iver) == 'V' ) then
            if ( s_eqi(syseqn,'NavierStokes') ) then
              term = ar*(nu_inv*wj*dvdx*wi)
              anab(2,1,iuse,ip) = anab(2,1,iuse,ip)+term
            end if
          end if

          if ( iprs>0 ) then
            if ( eqn(iprs) == 'P' ) then
              term = ar*dwjdx*qi
              anab(3,1,iuse,ip) = anab(3,1,iuse,ip)+term
            end if
          end if
!
!  Contributions of the JVER variable to various equations.
!
          if ( eqn(ihor) == 'U' ) then
            if ( s_eqi(syseqn,'NavierStokes') ) then
              term = ar*nu_inv*wj*dudy*wi
              anab(1,2,iuse,ip) = anab(1,2,iuse,ip)+term
            end if
          end if

          if ( eqn(iver) == 'V' ) then
            if ( s_eqi(syseqn,'NavierStokes') ) then
              term = ar*(dwjdx*dwidx+dwjdy*dwidy &
                +nu_inv*(uold*dwjdx+wj*dvdy+vold*dwjdy)*wi)
            else if ( s_eqi(syseqn,'Stokes') ) then
              term = ar*(dwjdx*dwidx+dwjdy*dwidy)
            end if
            anab(2,2,iuse,ip) = anab(2,2,iuse,ip)+term
          end if

          if ( iprs>0 ) then
            if ( eqn(iprs) == 'P' ) then
              term = ar*dwjdy*qi
              anab(3,2,iuse,ip) = anab(3,2,iuse,ip)+term
            end if
          end if
!
!  Contributions of the JPRS variable to various equations.
!
          if ( mod(jrow,2) == 1 .and. mod(jcol,2)==1 ) then

            if ( eqn(ihor) == 'U' ) then
              term = ar*nu_inv*dqjdx*wi
              anab(1,3,iuse,ip) = anab(1,3,iuse,ip)+term
            end if

            if ( eqn(iver) == 'V' ) then
              term = ar*nu_inv*dqjdy*wi
              anab(2,3,iuse,ip) = anab(2,3,iuse,ip)+term
            end if

          end if

        end do
      end do
    end do
  end do
!
!  Set up the equations that enforce boundary conditions.
!
  do ip = 1,np

    ihor = indx(ip,1)
    iver = indx(ip,2)
    iprs = indx(ip,3)

    if ( eqn(ihor) == 'UB' .or. eqn(ihor)=='UI'.or.eqn(ihor)=='UW' ) then
      anab(1,1,10,ip) = 1.0D+00
    end if

    if ( eqn(iver) == 'VB' .or. eqn(iver)=='VI'.or.eqn(iver)=='VW' ) then
      anab(2,2,10,ip) = 1.0D+00
    end if

    if ( iprs>0 ) then
      if ( eqn(iprs) == 'PB' ) then
        anab(3,3,10,ip) = 1.0D+00
      end if
    end if

  end do
  return
end
subroutine fx ( area,eqn,g,indx,ishapf,nelem,neqn,node,np,npar,nparb,nparf, &
  npe,para,phi,res,splflo,syseqn,tauflo,yc)

!*****************************************************************************80
!
!! FX computes the residual of the Navier Stokes or Stokes equations.
!
!  Discussion:
!
!    The Navier Stokes equations have the form:
!
!    Integral
!
!      dU/dx * dW/dx + dU/dy * dW/dy
!    + nu_inv * (U*dU/dx + V*dU/dy + dP/dx) * W dx dy = 0
!
!    Integral
!
!      dV/dx * dW/dx + dV/dy * dW/dy
!    + nu_inv * (U*dV/dx + V*dV/dy + dP/dy) * W dx dy = 0
!
!    Integral
!
!      (dU/dx + dV/dy) * Q dx dy = 0
!
!    The Stokes equations have the form:
!
!    Integral
!
!      dU/dx * dW/dx + dU/dy * dW/dy + nu_inv * dP/dx * W dx dy = 0
!
!    Integral
!
!      dV/dx * dW/dx + dV/dy * dW/dy + nu_inv * dP/dy * W dx dy = 0
!
!    Integral
!
!      (dU/dx + dV/dy) * Q dx dy = 0
!
!    Here W is a basis function for U and V, and Q is a basis
!    function for P.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) RES(NEQN), contains the value
!    of the residual.
!
  implicit none

  integer ( kind = 4 ) nelem
  integer ( kind = 4 ) neqn
  integer ( kind = 4 ) np
  integer ( kind = 4 ) npar
  integer ( kind = 4 ) nparf
  integer ( kind = 4 ) npe

  real ( kind = 8 ) ar
  real ( kind = 8 ) area(nelem,3)
  real ( kind = 8 ) dpdx
  real ( kind = 8 ) dpdy
  real ( kind = 8 ) dudx
  real ( kind = 8 ) dudy
  real ( kind = 8 ) dvdx
  real ( kind = 8 ) dvdy
  real ( kind = 8 ) dwidx
  real ( kind = 8 ) dwidy
  character ( len = 2 ) eqn(neqn)
  real ( kind = 8 ) g(neqn)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ielem
  integer ( kind = 4 ) ihor
  integer ( kind = 4 ) indx(np,3)
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) iprs
  integer ( kind = 4 ) iq
  integer ( kind = 4 ) iquad
  integer ( kind = 4 ) ishapf
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) iver
  logical s_eqi
  integer ( kind = 4 ) node(nelem,npe)
  integer ( kind = 4 ) nparb
  real ( kind = 8 ) p
  real ( kind = 8 ) para(npar)
  real ( kind = 8 ) phi(nelem,3,6,10)
  real ( kind = 8 ) qi
  real ( kind = 8 ) res(neqn)
  real ( kind = 8 ) nu_inv
  real ( kind = 8 ) splflo(4,nparf+2,0:nparf)
  character ( len = 20 ) syseqn
  real ( kind = 8 ) tauflo(nparf+2)
  real ( kind = 8 ) u
  real ( kind = 8 ) ubc
  real ( kind = 8 ) v
  real ( kind = 8 ) vbc
  real ( kind = 8 ) wi
  real ( kind = 8 ) yc(np)

  external s_eqi

  ival = 1
  call imemry('inc','Fx_calls',ival)

  nu_inv = para(nparf+nparb+1)

  res(1:neqn) = 0.0D+00
!
!  Consider an element.
!
  do ielem = 1,nelem
!
!  Evaluate the integrand at the quadrature points.
!
    do iquad = 1,3

      ar = area(ielem,iquad)
!
!  Evaluate P, U and V.
!
      call uvalq ( dpdx,dpdy,dudx,dudy,dvdx,dvdy,g,ielem,indx, &
        iquad,nelem,neqn,node,np,npe,p,phi,u,v)
!
!  Look at nearby basis functions.
!
      do iq = 1,6

        ip = node(ielem,iq)

        wi = phi(ielem,iquad,iq,1)
        dwidx = phi(ielem,iquad,iq,2)
        dwidy = phi(ielem,iquad,iq,3)
        qi = phi(ielem,iquad,iq,4)
!
!  The horizontal velocity equations.
!
        ihor = indx(ip,1)

        if ( eqn(ihor) == 'U' ) then

          if ( s_eqi(syseqn,'NavierStokes') ) then
            res(ihor) = res(ihor)+ar*(dudx*dwidx + dudy*dwidy &
              +nu_inv*(u*dudx+v*dudy+dpdx)*wi )
          else if ( s_eqi(syseqn,'Stokes') ) then
            res(ihor) = res(ihor)+ar*(dudx*dwidx + dudy*dwidy &
              +nu_inv*dpdx*wi )
          end if

        else if ( eqn(ihor) == 'UB' ) then

          res(ihor) = g(ihor)

        else if ( eqn(ihor) == 'UI' ) then

          call flo_spl_val ( ishapf, nparf, splflo, tauflo, ubc, vbc, yc(ip) )
          res(ihor) = g(ihor)-ubc

        else if ( eqn(ihor) == 'UW' ) then

          res(ihor) = g(ihor)

        end if
!
!  The vertical velocity equations.
!
        iver = indx(ip,2)

        if ( eqn(iver) == 'V' ) then

          if ( s_eqi(syseqn,'NavierStokes') ) then
            res(iver) = res(iver)+ar*(dvdx*dwidx + dvdy*dwidy &
              +nu_inv*(u*dvdx+v*dvdy+dpdy)*wi )
          else if ( s_eqi(syseqn,'Stokes') ) then
            res(iver) = res(iver)+ar*(dvdx*dwidx + dvdy*dwidy &
              +nu_inv*dpdy*wi )
          end if

        else if ( eqn(iver) == 'VB' ) then

          res(iver) = g(iver)

        else if ( eqn(iver) == 'VI' ) then

          call flo_spl_val ( ishapf, nparf, splflo, tauflo, ubc, vbc, yc(ip) )
          res(iver) = g(iver)-vbc

        else if ( eqn(iver) == 'VW' ) then

          res(iver) = g(iver)

        end if
!
!  The pressure equations.
!
        iprs = indx(ip,3)
        if ( iprs>0 ) then
          if ( eqn(iprs) == 'P' ) then
            res(iprs) = res(iprs)+ar*(dudx+dvdy)*qi
          else if ( eqn(iprs) == 'PB' ) then
            res(iprs) = g(iprs)
          end if
        end if

      end do
    end do
  end do

  return
end
subroutine get_cost ( cost, costb, costp, costu, costv, g, gtar, indx, &
  ishapb, neqn, np, nparb, nprof, ny, splbmp, taubmp, wateb, watep, wateu, &
  watev, xbl, xbr, ybl, ybr, yc )

!*****************************************************************************80
!
!! GET_COST determines the cost of a solution G given a target solution GTAR.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) neqn
  integer ( kind = 4 ) np
  integer ( kind = 4 ) nparb
  integer ( kind = 4 ) ny

  real ( kind = 8 ) cost
  real ( kind = 8 ) costb
  real ( kind = 8 ) costp
  real ( kind = 8 ) costu
  real ( kind = 8 ) costv
  real ( kind = 8 ) g(neqn)
  real ( kind = 8 ) gtar(neqn)
  integer ( kind = 4 ) indx(np,3)
  integer ( kind = 4 ) ishapb
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) nprof(2*ny-1)
  real ( kind = 8 ) splbmp(4,nparb+2,0:nparb)
  real ( kind = 8 ) taubmp(nparb+2)
  real ( kind = 8 ) wateb
  real ( kind = 8 ) watep
  real ( kind = 8 ) wateu
  real ( kind = 8 ) watev
  real ( kind = 8 ) xbl
  real ( kind = 8 ) xbr
  real ( kind = 8 ) ybl
  real ( kind = 8 ) ybr
  real ( kind = 8 ) yc(np)

  ival = 1
  call imemry ( 'inc', 'GET_COST_calls', ival )

  call bump_cost ( costb, ishapb, nparb, splbmp, taubmp, xbl, xbr, ybl, ybr )

  call disc_cost ( costp, costu, costv, g, gtar, indx, neqn, np, nprof, ny, yc )

  cost = wateb * costb + watep * costp + wateu * costu + watev * costv

  return
end
subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is an integer between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
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
!    Output, integer ( kind = 4 ) IUNIT.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5 and 6).
!
!    Otherwise, IUNIT is an integer between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  logical lopen

  iunit = 0

  do i = 1, 99

    if ( i /=  5 .and. i /= 6 ) then

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
subroutine getdu ( dpdyn,dudyn,dvdyn,etan,g,indx,isotri,nelem,neqn,node,np, &
  npe,numel,xc,xsin,yc)

!*****************************************************************************80
!
!! GETDU estimates spatial derivatives of state variables at nodes.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) nelem
  integer ( kind = 4 ) neqn
  integer ( kind = 4 ) np
  integer ( kind = 4 ) npe

  real ( kind = 8 ) det
  real ( kind = 8 ) detadx
  real ( kind = 8 ) detady
  real ( kind = 8 ) dpdx
  real ( kind = 8 ) dpdy
  real ( kind = 8 ) dpdyn(np)
  real ( kind = 8 ) dudx
  real ( kind = 8 ) dudy
  real ( kind = 8 ) dudyn(np)
  real ( kind = 8 ) dvdx
  real ( kind = 8 ) dvdy
  real ( kind = 8 ) dvdyn(np)
  real ( kind = 8 ) dxsidx
  real ( kind = 8 ) dxsidy
  real ( kind = 8 ) eta
  real ( kind = 8 ) etan(6)
  real ( kind = 8 ) g(neqn)
  integer ( kind = 4 ) ielem
  integer ( kind = 4 ) indx(np,3)
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) iq
  integer ( kind = 4 ) isotri(nelem)
  integer ( kind = 4 ) node(nelem,npe)
  integer ( kind = 4 ) numel(np)
  real ( kind = 8 ) p
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) xq
  real ( kind = 8 ) xsi
  real ( kind = 8 ) xsin(6)
  real ( kind = 8 ) yc(np)
  real ( kind = 8 ) yq
!
!  Estimate dUdY, dVdY, dPdY at each node.
!
!  We have a problem, because these quantities are not continuously
!  defined in all directions at nodes, which lie on the border between
!  two or more elements.
!
!  In order to assign some reasonable value to these quantities,
!  we look at each node, count the number of elements in which it
!  lies, evaluate the quantity within each element, and average the
!  result.
!
  numel(1:np) = 0
  dpdyn(1:np) = 0.0D+00
  dudyn(1:np) = 0.0D+00
  dvdyn(1:np) = 0.0D+00

  do ielem = 1,nelem

    do iq = 1,6

      ip = node(ielem,iq)
      numel(ip) = numel(ip)+1

      xq = xc(ip)
      yq = yc(ip)

      eta = etan(iq)
      xsi = xsin(iq)

      if ( isotri(ielem) == 2 ) then
        call trans ( det,detadx,detady,dxsidx,dxsidy,eta,ielem,nelem,node,np, &
          npe,xc,xsi,yc)
      end if

      call uval ( detadx,detady,dpdx,dpdy,dudx,dudy,dvdx,dvdy, &
        dxsidx,dxsidy,eta,g,ielem,indx,isotri,nelem,neqn, &
        node,np,npe,p,u,v,xc,xq,xsi,yc,yq)

      dpdyn(ip) = dpdyn(ip) + dpdy
      dudyn(ip) = dudyn(ip) + dudy
      dvdyn(ip) = dvdyn(ip) + dvdy
    end do

  end do
!
!  Take the average value of the quantities over all the
!  different elements along whose boundaries they are defined.
!
  do ip = 1, np
    dpdyn(ip) = dpdyn(ip) / real ( numel(ip), kind = 8 )
    dudyn(ip) = dudyn(ip) / real ( numel(ip), kind = 8 )
    dvdyn(ip) = dvdyn(ip) / real ( numel(ip), kind = 8 )
  end do

  return
end
subroutine getdu4(dudyn,np,nx,ny)

!*****************************************************************************80
!
!! GETDU4 uses the Zienkiewicz-Zhou technique to attempt to improve the
!  accuracy of the computed values of dPdY, dUdY and dVdY at corner
!  nodes.
!
!  It differs from GETDU2 in that it also modifies the values of
!  dUdY at the midside nodes.
!
!  Here is a picture of the element patch around an interior node "C":
!
!    .----.----NN---NNE--NNEE
!    |        /|        /|
!    |       / |       / |
!    |      /  |      /  |
!    |     /   |     /   |
!    .    NW   N    NE   NEE
!    |   /     |   /     |
!    |  /      |  /      |
!    | /       | /       |
!    |/        |/        |
!    WW---W----C----E----EE
!    |        /|        /|
!    |       / |       / |
!    |      /  |      /  |
!    |     /   |     /   |
!    SWW  SW   S    SE   .
!    |   /     |   /     |
!    |  /      |  /      |
!    | /       | /       |
!    |/        |/        |
!    SSWW-SSW--SS---.----.
!
!  The midside nodes SWW, SSW, SW, W, NW, S, N, SE, E, NE, NNE, NEE
!  are all quadrature points, at which a higher rater of convergence
!  is expected in the solution.  The nodes labeled SSWW, WW, SS, C,
!  NN, EE and NNEE are corner nodes, whose values of dUdY we wish
!  to improve.  On this step, we will try to compute improved values
!  of the node "C", and the side nodes.  The step will be repeated
!  so that every pressure node gets to be the "C" node of a patch
!  (except for corner nodes).  Nodes labeled "." are not included in
!  the patch, and are shown only for completeness.
!
!  The value of the quantity dUdY at each of the midside nodes is
!  taken as data, to be fitted by a least squares polynomial based
!  on the six basis functions (1, x, y, x*x, x*y, y*y).  For the
!  case where the node C is in the interior of the region, this
!  makes 12 equations in 6 unknowns.
!
!  The above diagram is assigned "logical" coordinates that
!  range from 0  <=  x,y <= 4, with SWW having coordinates (0,1)
!  and C having coordinates (2,2).  (These are really a sort of
!  XSI, ETA coordinate system, rather than a "physical" X, Y
!  system.)
!
!  Using standard discrete least squares methods, we write our
!  original set of equations as
!
!    AC x = b
!
!  where x is the set of 6 unknown polynomial coefficients,
!  row i of matrix AC is the value of (1, x, y, x*x, x*y, y*y)
!  at the i-th of our 12 quadrature points.
!
!  To get a solvable system, we multiply by the transpose of AC,
!  getting the square system known as the normal equations:
!
!    ACT AC x = ACT b
!
!  The matrix ACT AC only has to be set up and factored once, and
!  we can use the factors over and over for all the patches
!  associated with interior nodes.  (If we really like this method,
!  we should stop using the normal equations and use instead
!  a QR routine like the LINPACK SQRDC/SQRSL).
!
!  Once the coefficients of the polynomial are found, the least
!  squares polynomial may be evaluated at any point in the
!  "patch" of six elements.  In this routine, we choose to
!  evaluate the polynomial only at the point C, and use this
!  value in place of the value computed by the finite element
!  method.  The values at the midside nodes are unchanged.
!
!  Special treatment must be made for nodes which lie on the
!  north, east, west or south sides, or in the corners.
!
!  In the case of side nodes, we simply adapt the process to
!  account for 7 data values, using the six basis functions.
!
!  In the case of a corner node, we simply apply the relevant
!  value computed for the patch along the north or south boundary
!  nearest to the node.
!
!  In particular, separate matrices AN, AE, AW and AS must be set
!  up to represent the least squares linear systems that occur for
!  nodes which lie along the north, east, west or south boundaries.
!
!
!  Since the midside nodes can occur in more than one element,
!  we ADD the improved estimate to a running total, and average
!  later.
!
!  Note: this routine explicitly assumes that the nodes are ordered
!  from top to bottom, then from left to right, as in this example:
!
!    5 10 15
!    4  9 14
!    3  8 13
!    2  7 12
!    1  6 11
!
!    It also assumes that EVERY node has an associated value of DUDYN,
!    and that the pressure nodes, at which we will smooth the value of
!    DUDYN, are exactly those nodes in an odd row and odd column.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, REAL DUDYN(NP).
!    On input, DUDYN contains the original value of DUDY at every node.
!    On output, DUDYN contains the smoothed values of dUdY.
!
!    Input, INTEGER NP, the number of nodes.
!
!    Input, INTEGER NX, the number of columns of pressure nodes
!    in the X direction.  The actual number of columns of nodes
!    is 2*NX-1.
!
!    Input, INTEGER NY, the number of rows of pressure nodes
!    in the Y direction.  The actual number of rows of nodes
!    is 2*NY-1.
!
  implicit none

  integer ( kind = 4 ), parameter :: ncoef = 6
  integer ( kind = 4 ), parameter :: npmax = 7889

  integer ( kind = 4 ) np

  real ( kind = 8 ) ac(12,ncoef)
  real ( kind = 8 ) ae(7,ncoef)
  real ( kind = 8 ) an(7,ncoef)
  real ( kind = 8 ) as(7,ncoef)
  real ( kind = 8 ) aw(7,ncoef)
  real ( kind = 8 ) atac(ncoef,ncoef)
  real ( kind = 8 ) atae(ncoef,ncoef)
  real ( kind = 8 ) atan(ncoef,ncoef)
  real ( kind = 8 ) atas(ncoef,ncoef)
  real ( kind = 8 ) ataw(ncoef,ncoef)
  real ( kind = 8 ) b(ncoef)
  real ( kind = 8 ) dat(12)
  real ( kind = 8 ) dudyn(np)
  real ( kind = 8 ) dudyn2(npmax)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) icol
  integer ( kind = 4 ) ie
  integer ( kind = 4 ) iee
  integer ( kind = 4 ) in
  integer ( kind = 4 ) ine
  integer ( kind = 4 ) inee
  integer ( kind = 4 ) info
  integer ( kind = 4 ) inn
  integer ( kind = 4 ) inne
  integer ( kind = 4 ) inw
  integer ( kind = 4 ) ipivc(ncoef)
  integer ( kind = 4 ) ipive(ncoef)
  integer ( kind = 4 ) ipivn(ncoef)
  integer ( kind = 4 ) ipivs(ncoef)
  integer ( kind = 4 ) ipivw(ncoef)
  integer ( kind = 4 ) irow
  integer ( kind = 4 ) is
  integer ( kind = 4 ) ise
  integer ( kind = 4 ) iss
  integer ( kind = 4 ) issw
  integer ( kind = 4 ) isw
  integer ( kind = 4 ) isww
  integer ( kind = 4 ) iw
  integer ( kind = 4 ) iww
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) ndat
  integer ( kind = 4 ) nodat(12)
  integer ( kind = 4 ) nrhs
  integer ( kind = 4 ) nrep(npmax)
  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny
  real ( kind = 8 ) temp
  real ( kind = 8 ) x
  real ( kind = 8 ) xdatc(12)
  real ( kind = 8 ) xdatn(7)
  real ( kind = 8 ) xdate(7)
  real ( kind = 8 ) xdatw(7)
  real ( kind = 8 ) xdats(7)
  real ( kind = 8 ) y
  real ( kind = 8 ) ydatc(12)
  real ( kind = 8 ) ydatn(7)
  real ( kind = 8 ) ydate(7)
  real ( kind = 8 ) ydatw(7)
  real ( kind = 8 ) ydats(7)

  if ( np>npmax ) then
    write ( *, * ) ' '
    write ( *, * ) 'GETDU4 - Fatal error!'
    write ( *, * ) '  NP must be less than NPMAX, but'
    write ( *, * ) '  NP = ',np
    write ( *, * ) '  NPMAX = ',npmax
    write ( *, * ) '  Change NP or NPMAX and try again!'
    stop
  end if

  do i = 1,np
    nrep(i) = 0
    dudyn2(i) = 0.0D+00
  end do
!
!  Here is a picture of the element patch around a node on the
!  northern boundary:
!
!
!    WW---W----C----E----EE
!    |        /|        /|
!    |       / |       / |
!    |      /  |      /  |
!    |     /   |     /   |
!    SWW  SW   S    SE   .
!    |   /     |   /     |
!    |  /      |  /      |
!    | /       | /       |
!    |/        |/        |
!    SSWW-SSW--SS---.----.
!
!  NORTH
!
  ndat = 7

  xdatn(1) = 0.0D+00
  ydatn(1) = 1.0D+00

  xdatn(2) = 1.0D+00
  ydatn(2) = 0.0D+00

  xdatn(3) = 1.0D+00
  ydatn(3) = 1.0D+00

  xdatn(4) = 1.0D+00
  ydatn(4) = 2.0D+00

  xdatn(5) = 2.0D+00
  ydatn(5) = 1.0D+00

  xdatn(6) = 3.0D+00
  ydatn(6) = 1.0D+00

  xdatn(7) = 3.0D+00
  ydatn(7) = 2.0D+00
!
!  Compute matrix of basis polynomials (1,x,y,x**2,xy,y**2) evaluated
!  at the data points.
!
  do i = 1,ndat
    an(i,1) = 1.0D+00
    an(i,2) = xdatn(i)
    an(i,3) = ydatn(i)
    an(i,4) = xdatn(i)**2
    an(i,5) = xdatn(i)*ydatn(i)
    an(i,6) = ydatn(i)**2
  end do
!
!  Compute the normal equations matrix.
!
  do i = 1,ncoef
    do j = 1,ncoef
      atan(i,j) = 0.0D+00
      do k = 1,ndat
        atan(i,j) = atan(i,j)+an(k,i)*an(k,j)
      end do
    end do
  end do
!
!  Factor the normal equations matrix.
!
  call sgetrf(ncoef,ncoef,atan,ncoef,ipivn,info)

  if ( info/= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'GETDU2 - Fatal error!'
    write ( *, * ) '  SGETRF returned INFO = ',info
    write ( *, * ) '  while factoring matrix ATAN.'
    stop
  end if
!
!  EAST
!
!  Here is a picture of the element patch around a node
!  on the eastern boundary:
!
!    .----.----NN
!    |        /|
!    |       / |
!    |      /  |
!    |     /   |
!    .    NW   N
!    |   /     |
!    |  /      |
!    | /       |
!    |/        |
!    WW---W----C
!    |        /|
!    |       / |
!    |      /  |
!    |     /   |
!    SWW  SW   S
!    |   /     |
!    |  /      |
!    | /       |
!    |/        |
!    SSWW-SSW--SS
!
  ndat = 7

  xdate(1) = 0.0D+00
  ydate(1) = 1.0D+00

  xdate(2) = 1.0D+00
  ydate(2) = 0.0D+00

  xdate(3) = 1.0D+00
  ydate(3) = 1.0D+00

  xdate(4) = 1.0D+00
  ydate(4) = 2.0D+00

  xdate(5) = 1.0D+00
  ydate(5) = 3.0D+00

  xdate(6) = 2.0D+00
  ydate(6) = 1.0D+00

  xdate(7) = 2.0D+00
  ydate(7) = 3.0D+00
!
!  Compute matrix of basis polynomials (1,x,y,x**2,xy,y**2) evaluated
!  at the data points.
!
  do i = 1,ndat
    ae(i,1) = 1.0D+00
    ae(i,2) = xdate(i)
    ae(i,3) = ydate(i)
    ae(i,4) = xdate(i)**2
    ae(i,5) = xdate(i)*ydate(i)
    ae(i,6) = ydate(i)**2
  end do
!
!  Compute the normal equations matrix.
!
  do i = 1,ncoef
    do j = 1,ncoef
      atae(i,j) = 0.0D+00
      do k = 1,ndat
        atae(i,j) = atae(i,j)+ae(k,i)*ae(k,j)
      end do
    end do
  end do
!
!  Factor the normal equations matrix.
!
  call sgetrf(ncoef,ncoef,atae,ncoef,ipive,info)

  if ( info/= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'GETDU2 - Fatal error!'
    write ( *, * ) '  SGETRF returned INFO = ',info
    write ( *, * ) '  while factoring matrix ATAE.'
    stop
  end if
!
!  WEST
!
!  Here is a picture of the element patch around a node on the
!  western boundary:
!
!    NN---NNE--NNEE
!    |        /|
!    |       / |
!    |      /  |
!    |     /   |
!    N    NE   NEE
!    |   /     |
!    |  /      |
!    | /       |
!    |/        |
!    C----E----EE
!    |        /|
!    |       / |
!    |      /  |
!    |     /   |
!    S    SE   .
!    |   /     |
!    |  /      |
!    | /       |
!    |/        |
!    SS---.----.
!
  ndat = 7

  xdatw(1) = 2.0D+00
  ydatw(1) = 1.0D+00

  xdatw(2) = 2.0D+00
  ydatw(2) = 3.0D+00

  xdatw(3) = 3.0D+00
  ydatw(3) = 1.0D+00

  xdatw(4) = 3.0D+00
  ydatw(4) = 2.0D+00

  xdatw(5) = 3.0D+00
  ydatw(5) = 3.0D+00

  xdatw(6) = 3.0D+00
  ydatw(6) = 4.0D+00

  xdatw(7) = 4.0D+00
  ydatw(7) = 3.0D+00
!
!  Compute matrix of basis polynomials (1,x,y,x**2,xy,y**2) evaluated
!  at the data points.
!
  do i = 1,ndat
    aw(i,1) = 1.0D+00
    aw(i,2) = xdatw(i)
    aw(i,3) = ydatw(i)
    aw(i,4) = xdatw(i)**2
    aw(i,5) = xdatw(i)*ydatw(i)
    aw(i,6) = ydatw(i)**2
  end do
!
!  Compute the normal equations matrix.
!
  do i = 1,ncoef
    do j = 1,ncoef
      ataw(i,j) = 0.0D+00
      do k = 1,ndat
        ataw(i,j) = ataw(i,j)+aw(k,i)*aw(k,j)
      end do
    end do
  end do
!
!  Factor the normal equations matrix.
!
  call sgetrf ( ncoef, ncoef, ataw, ncoef, ipivw, info )

  if ( info/= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'GETDU2 - Fatal error!'
    write ( *, * ) '  SGETRF returned INFO = ',info
    write ( *, * ) '  while factoring matrix ATAW.'
    stop
  end if
!
!  SOUTH
!
!  Here is a picture of the element patch around a node on the
!  southern boundary:
!
!    .----.----NN---NNE--NNEE
!    |        /|        /|
!    |       / |       / |
!    |      /  |      /  |
!    |     /   |     /   |
!    .    NW   N    NE   NEE
!    |   /     |   /     |
!    |  /      |  /      |
!    | /       | /       |
!    |/        |/        |
!    WW---W----C----E----EE
!
  ndat = 7

  xdats(1) = 1.0D+00
  ydats(1) = 2.0D+00

  xdats(2) = 1.0D+00
  ydats(2) = 3.0D+00

  xdats(3) = 2.0D+00
  ydats(3) = 3.0D+00

  xdats(4) = 3.0D+00
  ydats(4) = 2.0D+00

  xdats(5) = 3.0D+00
  ydats(5) = 3.0D+00

  xdats(6) = 3.0D+00
  ydats(6) = 4.0D+00

  xdats(7) = 4.0D+00
  ydats(7) = 3.0D+00
!
!  Compute matrix of basis polynomials (1,x,y,x**2,xy,y**2) evaluated
!  at the data points.
!
  do i = 1,ndat
    as(i,1) = 1.0D+00
    as(i,2) = xdats(i)
    as(i,3) = ydats(i)
    as(i,4) = xdats(i)**2
    as(i,5) = xdats(i)*ydats(i)
    as(i,6) = ydats(i)**2
  end do
!
!  Compute the normal equations matrix.
!
  do i = 1,ncoef
    do j = 1,ncoef
      atas(i,j) = 0.0D+00
      do k = 1,ndat
        atas(i,j) = atas(i,j)+as(k,i)*as(k,j)
      end do
    end do
  end do
!
!  Factor the normal equations matrix.
!
  call sgetrf(ncoef,ncoef,atas,ncoef,ipivs,info)

  if ( info/= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'GETDU2 - Fatal error!'
    write ( *, * ) '  SGETRF returned INFO = ',info
    write ( *, * ) '  while factoring matrix ATAS.'
    stop
  end if
!
!  CENTRAL.
!
  ndat = 12

  xdatc(1) = 0.0D+00
  ydatc(1) = 1.0D+00

  xdatc(2) = 1.0D+00
  ydatc(2) = 0.0D+00

  xdatc(3) = 1.0D+00
  ydatc(3) = 1.0D+00

  xdatc(4) = 1.0D+00
  ydatc(4) = 2.0D+00

  xdatc(5) = 1.0D+00
  ydatc(5) = 3.0D+00

  xdatc(6) = 2.0D+00
  ydatc(6) = 1.0D+00

  xdatc(7) = 2.0D+00
  ydatc(7) = 3.0D+00

  xdatc(8) = 3.0D+00
  ydatc(8) = 1.0D+00

  xdatc(9) = 3.0D+00
  ydatc(9) = 2.0D+00

  xdatc(10) = 3.0D+00
  ydatc(10) = 3.0D+00

  xdatc(11) = 3.0D+00
  ydatc(11) = 4.0D+00

  xdatc(12) = 4.0D+00
  ydatc(12) = 3.0D+00
!
!  Compute matrix of basis polynomials (1,x,y,x**2,xy,y**2) evaluated
!  at the data points.
!
  do i = 1,ndat
    ac(i,1) = 1.0D+00
    ac(i,2) = xdatc(i)
    ac(i,3) = ydatc(i)
    ac(i,4) = xdatc(i)**2
    ac(i,5) = xdatc(i)*ydatc(i)
    ac(i,6) = ydatc(i)**2
  end do
!
!  Compute the normal equations matrix.
!
  do i = 1,ncoef
    do j = 1,ncoef
      atac(i,j) = 0.0D+00
      do k = 1,ndat
        atac(i,j) = atac(i,j)+ac(k,i)*ac(k,j)
      end do
    end do
  end do
!
!  Factor the normal equations matrix.
!
  call sgetrf(ncoef,ncoef,atac,ncoef,ipivc,info)

  if ( info/= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'GETDU2 - Fatal error!'
    write ( *, * ) '  SGETRF returned INFO = ',info
    write ( *, * ) '  while factoring matrix ATAC.'
    stop
  end if
!
!  Each computation has its "center" at a pressure node, so we
!  can organize our do loops by counting through each presure node.
!
  do icol = 1,nx
    do irow = 1,ny

      ic = 2*(icol-1)*(2*ny-1)+(2*irow-1)

      isww = ic-2*(2*ny-1)-1
      iww = isww+1
      issw = ic-(2*ny-1)-2
      isw = issw+1
      iw = isw+1
      inw = iw+1
      iss = ic-2
      is = iss+1
      in = ic+1
      inn = ic+2
      ise = ic+(2*ny-1)-1
      ie = ise+1
      ine = ie+1
      inne = ine+1
      iee = ic+2*(2*ny-1)
      inee = iee+1
!
!  CASE: CORNER NODES, SKIP EM
!
      if ( icol == 1 .and. irow==1 ) then
      else if ( icol == 1 .and. irow==ny ) then
      else if ( icol == nx .and. irow==1 ) then
      else if ( icol == nx .and. irow==ny ) then
!
!  CASE: THE NORTHERN BOUNDARY
!
      else if ( irow == ny ) then

        nodat(1) = isww
        nodat(2) = issw
        nodat(3) = isw
        nodat(4) = iw
        nodat(5) = is
        nodat(6) = ise
        nodat(7) = ie

        ndat = 7
!
!  Copy out the value of dUdY at the data points.
!
        do i = 1,ndat
          dat(i) = dudyn(nodat(i))
        end do
!
!  Compute right hand side of normal equations.
!
        do i = 1,ncoef
          b(i) = 0.0D+00
          do j = 1,ndat
            b(i) = b(i)+an(j,i)*dat(j)
          end do
        end do
!
!  Solve the normal equations.
!
        nrhs = 1
        call sgetrs('N',ncoef,nrhs,atan,ncoef,ipivn,b,ncoef,info)
!
!  Now evaluate the least squares polynomial at the nodes.
!
        if ( icol == 2 ) then

          x = 0.0D+00
          y = 2.0D+00
          call lspoly(b,ncoef,x,y,temp)
          nrep(iww) = nrep(iww)+1
          dudyn2(iww) = dudyn2(iww)+temp

          x = 0.0D+00
          y = 1.0D+00
          call lspoly(b,ncoef,x,y,temp)
          nrep(isww) = nrep(isww)+1
          dudyn2(isww) = dudyn2(isww)+temp

        end if

        do i = 1,ndat
          x = xdatn(i)
          y = ydatn(i)
          call lspoly(b,ncoef,x,y,temp)
          nrep(nodat(i)) = nrep(nodat(i))+1
          dudyn2(nodat(i)) = dudyn2(nodat(i))+temp
        end do

        x = 2.0D+00
        y = 2.0D+00
        call lspoly(b,ncoef,x,y,temp)
        nrep(ic) = nrep(ic)+1
        dudyn2(ic) = dudyn2(ic)+temp

        if ( icol == nx-1 ) then
          x = 4.0D+00
          y = 2.0D+00
          call lspoly(b,ncoef,x,y,temp)
          nrep(iee) = nrep(iee)+1
          dudyn2(iee) = dudyn2(iee)+temp
        end if
!
!  CASE: THE EASTERN BOUNDARY
!
      else if ( icol == nx ) then

        nodat(1) = isww
        nodat(2) = issw
        nodat(3) = isw
        nodat(4) = iw
        nodat(5) = inw
        nodat(6) = is
        nodat(7) = in

        ndat = 7
!
!  Copy out the value of dUdY at the data points.
!
        do i = 1,ndat
          dat(i) = dudyn(nodat(i))
        end do
!
!  Compute right hand side of normal equations.
!
        do i = 1,ncoef
          b(i) = 0.0D+00
          do j = 1,ndat
            b(i) = b(i)+ae(j,i)*dat(j)
          end do
        end do
!
!  Solve the normal equations.
!
        nrhs = 1
        call sgetrs('N',ncoef,nrhs,atae,ncoef,ipive,b,ncoef,info)
!
!  Now evaluate the least squares polynomial at the nodes.
!
        if ( irow == ny-1 ) then

          x = 2.0D+00
          y = 4.0D+00
          call lspoly(b,ncoef,x,y,temp)
          nrep(inn) = nrep(inn)+1
          dudyn2(inn) = dudyn2(inn)+temp

        end if

        if ( irow == ny-1 ) then

          x = 2.0D+00
          y = 1.0D+00
          call lspoly(b,ncoef,x,y,temp)
          nrep(is) = nrep(is)+1
          dudyn2(is) = dudyn2(is)+temp

          x = 2.0D+00
          y = 0.0D+00
          call lspoly(b,ncoef,x,y,temp)
          nrep(iss) = nrep(iss)+1
          dudyn2(iss) = dudyn2(iss)+temp

        end if

        do i = 1,ndat
          x = xdate(i)
          y = ydate(i)
          call lspoly(b,ncoef,x,y,temp)
          nrep(nodat(i)) = nrep(nodat(i))+1
          dudyn2(nodat(i)) = dudyn2(nodat(i))+temp
        end do

        x = 2.0D+00
        y = 2.0D+00
        call lspoly(b,ncoef,x,y,temp)
        nrep(ic) = nrep(ic)+1
        dudyn2(ic) = dudyn2(ic)+temp
!
!  CASE: THE WESTERN BOUNDARY
!
      else if ( icol == 1 ) then

        nodat(1) = is
        nodat(2) = in
        nodat(3) = ise
        nodat(4) = ie
        nodat(5) = ine
        nodat(6) = inne
        nodat(7) = inee

        ndat = 7
!
!  Copy out the value of dUdY at the data points.
!
        do i = 1,ndat
          dat(i) = dudyn(nodat(i))
        end do
!
!  Compute right hand side of normal equations.
!
        do i = 1,ncoef
          b(i) = 0.0D+00
          do j = 1,ndat
            b(i) = b(i)+aw(j,i)*dat(j)
          end do
        end do
!
!  Solve the normal equations.
!
        nrhs = 1
        call sgetrs('N',ncoef,nrhs,ataw,ncoef,ipivw,b,ncoef,info)
!
!  Now evaluate the least squares polynomial at the nodes.
!
        if ( irow == ny-1 ) then

          x = 2.0D+00
          y = 4.0D+00
          call lspoly(b,ncoef,x,y,temp)
          nrep(inn) = nrep(inn)+1
          dudyn2(inn) = dudyn2(inn)+temp

          x = 3.0D+00
          y = 4.0D+00
          call lspoly(b,ncoef,x,y,temp)
          nrep(inne) = nrep(inne)+1
          dudyn2(inne) = dudyn2(inne)+temp

        end if

        if ( irow == 2 ) then

          x = 2.0D+00
          y = 0.0D+00
          call lspoly(b,ncoef,x,y,temp)
          nrep(iss) = nrep(iss)+1
          dudyn2(iss) = dudyn2(iss)+temp

        end if

        do i = 1,ndat
          x = xdatw(i)
          y = ydatw(i)
          call lspoly(b,ncoef,x,y,temp)
          nrep(nodat(i)) = nrep(nodat(i))+1
          dudyn2(nodat(i)) = dudyn2(nodat(i))+temp
        end do

        x = 2.0D+00
        y = 2.0D+00
        call lspoly(b,ncoef,x,y,temp)
        nrep(ic) = nrep(ic)+1
        dudyn2(ic) = dudyn2(ic)+temp
!
!  CASE: THE SOUTHERN BOUNDARY
!
      else if ( irow == 1 ) then

        nodat(1) = iw
        nodat(2) = inw
        nodat(3) = in
        nodat(4) = ie
        nodat(5) = ine
        nodat(6) = inne
        nodat(7) = inee

        ndat = 7
!
!  Copy out the value of dUdY at the data points.
!
        do i = 1,ndat
          dat(i) = dudyn(nodat(i))
        end do
!
!  Compute right hand side of normal equations.
!
        do i = 1,ncoef
          b(i) = 0.0D+00
          do j = 1,ndat
            b(i) = b(i)+as(j,i)*dat(j)
          end do
        end do
!
!  Solve the normal equations.
!
        nrhs = 1
        call sgetrs('N',ncoef,nrhs,atas,ncoef,ipivs,b,ncoef,info)
!
!  Now evaluate the least squares polynomial at the nodes.
!
        if ( icol == 2 ) then
          x = 0.0D+00
          y = 2.0D+00
          call lspoly(b,ncoef,x,y,temp)
          nrep(iww) = nrep(iww)+1
          dudyn2(iww) = dudyn2(iww)+temp
        end if

        do i = 1,ndat
          x = xdats(i)
          y = ydats(i)
          call lspoly(b,ncoef,x,y,temp)
          nrep(nodat(i)) = nrep(nodat(i))+1
          dudyn2(nodat(i)) = dudyn2(nodat(i))+temp
        end do

        if ( icol == nx-1 ) then

          x = 4.0D+00
          y = 2.0D+00
          call lspoly(b,ncoef,x,y,temp)
          nrep(iee) = nrep(iee)+1
          dudyn2(iee) = dudyn2(iee)+temp

          x = 4.0D+00
          y = 3.0D+00
          call lspoly(b,ncoef,x,y,temp)
          nrep(inee) = nrep(inee)+1
          dudyn2(inee) = dudyn2(inee)+temp

        end if

        x = 2.0D+00
        y = 2.0D+00
        call lspoly(b,ncoef,x,y,temp)
        nrep(ic) = nrep(ic)+1
        dudyn2(ic) = dudyn2(ic)+temp
!
!  CASE: CENTRAL
!
      else

        nodat(1) = isww
        nodat(2) = issw
        nodat(3) = isw
        nodat(4) = iw
        nodat(5) = inw
        nodat(6) = is
        nodat(7) = in
        nodat(8) = ise
        nodat(9) = ie
        nodat(10) = ine
        nodat(11) = inne
        nodat(12) = inee

        ndat = 12
!
!  Copy out the value of dUdY at the data points.
!
        do i = 1,ndat
          dat(i) = dudyn(nodat(i))
        end do
!
!  Compute right hand side of normal equations.
!
        do i = 1,ncoef
          b(i) = 0.0D+00
          do j = 1,ndat
            b(i) = b(i)+ac(j,i)*dat(j)
          end do
        end do
!
!  Solve the normal equations.
!
        nrhs = 1
        call sgetrs('N',ncoef,nrhs,atac,ncoef,ipivc,b,ncoef,info)
!
!  Now evaluate the least squares polynomial at the nodes.
!
        do i = 1,ndat
          x = xdatc(i)
          y = ydatc(i)
          call lspoly(b,ncoef,x,y,temp)
          nrep(nodat(i)) = nrep(nodat(i))+1
          dudyn2(nodat(i)) = dudyn2(nodat(i))+temp
        end do

        x = 2.0D+00
        y = 2.0D+00
        call lspoly(b,ncoef,x,y,temp)
        nrep(ic) = nrep(ic)+1
        dudyn2(ic) = dudyn2(ic)+temp

      end if

    end do
  end do
!
!  Copy the smoothed data.
!
  do i = 1,np

    if ( nrep(i) == 0 ) then
      write ( *, * ) ' '
      write ( *, * ) 'GETDU4 - Fatal error!'
      write ( *, * ) '  NREP  <= 0 for node ',I
      stop
    end if

    dudyn(i) = dudyn2(i) / real ( nrep(i), kind = 8 )
  end do

  return
end
subroutine getfix ( dpdyn, dudyn, dvdyn, dydpn, gradf, indx, iopt, ishapb, &
  neqn, np, npar, nparb, nparf, splbmp, taubmp, xbl, xbr, xc, yc )

!*****************************************************************************80
!
!! GETFIX corrects the finite difference estimate of the sensitivities.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) GRADF(NEQN,NPAR), a correction to the finite difference
!    estimate of the sensitivities.
!
  implicit none

  integer ( kind = 4 ) neqn
  integer ( kind = 4 ) np
  integer ( kind = 4 ) npar
  integer ( kind = 4 ) nparb
  integer ( kind = 4 ) nparf

  real ( kind = 8 ) dpdyn(np)
  real ( kind = 8 ) dudyn(np)
  real ( kind = 8 ) dvdyn(np)
  real ( kind = 8 ) dydpn(np,nparb)
  real ( kind = 8 ) gradf(neqn,npar)
  integer ( kind = 4 ) ihor
  integer ( kind = 4 ) indx(np,3)
  integer ( kind = 4 ) iopt(npar)
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) ipar
  integer ( kind = 4 ) iparb
  integer ( kind = 4 ) iprs
  integer ( kind = 4 ) ishapb
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) iver
  integer ( kind = 4 ) jderiv
  real ( kind = 8 ) shape
  real ( kind = 8 ) splbmp(4,nparb+2,0:nparb)
  real ( kind = 8 ) taubmp(nparb+2)
  real ( kind = 8 ) xbl
  real ( kind = 8 ) xbr
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) xq
  real ( kind = 8 ) ybot
  real ( kind = 8 ) yc(np)
  real ( kind = 8 ) yq
  real ( kind = 8 ) yvec(nparb+2)

  ival = 1
  call imemry('inc','GetFix_calls',ival)
!
!  Compute dYdParameter
!
  do ip = 1,np

    xq = xc(ip)
    yq = yc(ip)

    do iparb = 1,nparb

      dydpn(ip,iparb) = 0.0D+00
!
!  Nodes only move if their X coordinate is (strictly) between XBL
!  and XBR.
!
      if ( xbl < xq .and. xq < xbr ) then
!
!  Compute YBOT, the height of the bump below the current point.
!  and SHAPE, the height at XQ of the basis function
!  associated with parameter IPARB.
!
        if ( ishapb == 1 ) then
          yvec(1:nparb+2) = splbmp(1,1:nparb+2,0)
          call plval ( nparb+2, xq, taubmp, ybot, yvec )
          call plval1 ( iparb+1, nparb+2, xq, taubmp, shape )
        else if ( ishapb == 2 ) then
          yvec(1:nparb+2) = splbmp(1,1:nparb+2,0)
          call pqval ( nparb+2, xq, taubmp, ybot, yvec )
          call pqval1 ( iparb+1, nparb+2, xq, taubmp, shape )
        else if ( ishapb == 3 ) then
          jderiv = 0
          call ppvalu ( taubmp, splbmp(1,1,0), nparb+1, 4, xq, jderiv, ybot )
          jderiv = 0
          call ppvalu ( taubmp, splbmp(1,1,iparb), nparb+1, 4, xq, jderiv, shape )
        end if

        dydpn(ip,iparb) = ( ( 3.0D+00 - yq ) / ( 3.0D+00 - ybot ) ) * shape

      end if

   end do

  end do
!
!  Use these results to compute corrections for those entries of the
!  finite difference vector that correspond to shape parameters, and which
!  are based at nodes which move during the estimation.
!
  do ipar = nparf+1,nparf+nparb

    if ( iopt(ipar) == 1 ) then

      do ip = 1,np

        ihor = indx(ip,1)
        gradf(ihor,ipar) = dudyn(ip)*dydpn(ip,ipar-nparf)

        iver = indx(ip,2)
        gradf(iver,ipar) = dvdyn(ip)*dydpn(ip,ipar-nparf)

        iprs = indx(ip,3)
        if ( 0 < iprs ) then
          gradf(iprs,ipar) = dpdyn(ip)*dydpn(ip,ipar-nparf)
        end if

      end do

    end if

  end do

  return
end
subroutine getgrd ( a,area,cost,disjac,dpara3,eqn,etan,etaq,flarea,g, &
  gdif,gtar,idfd,ierror,igrid,ijac,indx,iopt,ipivot,ishapb,ishapf, &
  isotri,iwrite,jjac,maxnew,nelem,neqn,nlband,node,np,npar,nparb,nparf,npe, &
  nprof, &
  nrow,nx,ny,para,parjac,phi,res,splbmp,splflo,syseqn,taubmp,tauflo,tolnew, &
  wateb,watep,wateu,watev,wquad,xbl,xbord,xbr,xc,xquad,xsin,xsiq,ybl,ybord, &
  ybr,yc,yquad)

!*****************************************************************************80
!
!! GETGRD estimates the derivatives of the variables with respect to parameters.
!
!  Discussion:
!
!    Finite differences are used.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) PARA(NPAR), the current parameter values.
!
!    G the corresponding state variables, and
!
!    COST the resulting cost.
!
  implicit none

  integer ( kind = 4 ) nelem
  integer ( kind = 4 ) neqn
  integer ( kind = 4 ) np
  integer ( kind = 4 ) npar
  integer ( kind = 4 ) nparb
  integer ( kind = 4 ) nparf
  integer ( kind = 4 ) npe
  integer ( kind = 4 ) nrow
  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny

  real ( kind = 8 ) a(nrow,neqn)
  real ( kind = 8 ) area(nelem,3)
  real ( kind = 8 ) cost
  real ( kind = 8 ) cost1
  real ( kind = 8 ) costb
  real ( kind = 8 ) costp
  real ( kind = 8 ) costu
  real ( kind = 8 ) costv
  real ( kind = 8 ) disjac
  real ( kind = 8 ) dpara3(npar)
  character ( len = 2 ) eqn(neqn)
  real ( kind = 8 ) etan(6)
  real ( kind = 8 ) etaq(3)
  real ( kind = 8 ) flarea
  real ( kind = 8 ) g(neqn)
  real ( kind = 8 ) g1(neqn)
  real ( kind = 8 ) gdif(neqn,npar)
  real ( kind = 8 ) gtar(neqn)
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) idfd
  integer ( kind = 4 ) ieqn
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) igrid
  integer ( kind = 4 ) ijac
  integer ( kind = 4 ) indx(np,3)
  integer ( kind = 4 ) iopt(npar)
  integer ( kind = 4 ) ipar
  integer ( kind = 4 ) ipivot(neqn)
  integer ( kind = 4 ) ishapb
  integer ( kind = 4 ) ishapf
  integer ( kind = 4 ) isotri(nelem)
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) iwrite
  integer ( kind = 4 ) jjac
  logical lgdif
  integer ( kind = 4 ) maxnew
  integer ( kind = 4 ) nlband
  integer ( kind = 4 ) node(nelem,npe)
  integer ( kind = 4 ) nprof(2*ny-1)
  real ( kind = 8 ) para(npar)
  real ( kind = 8 ) parjac(npar)
  real ( kind = 8 ) phi(nelem,3,6,10)
  real ( kind = 8 ) psave
  real ( kind = 8 ) res(neqn)
  real ( kind = 8 ) splbmp(4,nparb+2,0:nparb)
  real ( kind = 8 ) splflo(4,nparf+2,0:nparf)
  character ( len = 20 ) syseqn
  real ( kind = 8 ) taubmp(nparb+2)
  real ( kind = 8 ) tauflo(nparf+2)
  real ( kind = 8 ) tolnew
  real ( kind = 8 ) wateb
  real ( kind = 8 ) watep
  real ( kind = 8 ) wateu
  real ( kind = 8 ) watev
  real ( kind = 8 ) wquad(3)
  real ( kind = 8 ) xbl
  real ( kind = 8 ) xbord(nx)
  real ( kind = 8 ) xbr
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) xquad(nelem,3)
  real ( kind = 8 ) xsin(6)
  real ( kind = 8 ) xsiq(3)
  real ( kind = 8 ) ybl
  real ( kind = 8 ) ybord(ny)
  real ( kind = 8 ) ybr
  real ( kind = 8 ) yc(np)
  real ( kind = 8 ) yquad(nelem,3)

  ival = 1
  call imemry ( 'inc', 'GetGrd_calls', ival )
!
!  Set a logical flag, so that all routines know that we are doing
!  a finite difference computation.
!
  lgdif = .true.
  call lmemry ( 'set', 'fd_grad', lgdif )
!
!  Consider each parameter.
!
  do ipar = 1, npar
!
!  If the parameter is allowed to vary,
!
    if ( iopt(ipar) == 1 ) then

      g1(1:neqn) = g(1:neqn)
!
!  ...then perturb the parameter by the amount H.
!
      psave = para(ipar)
      h = sign ( 1.0D+00, psave ) * ( abs ( psave ) + 1.0D+00 ) * sqrt ( epsilon ( h ) )
      para(ipar) = para(ipar) + h
!
!  Solve the flow problem.
!
      ival = 1
      call imemry('inc','Points_FD',ival)

      call flosol ( a,area,disjac,eqn,etan,etaq,flarea,g1,ierror,igrid,ijac, &
        indx,ipivot,ishapb,ishapf,isotri,iwrite,jjac,maxnew,nelem,neqn,nlband, &
        node,np,npar,nparb,nparf,npe,nrow,nx,ny,para,parjac,phi,res,splbmp, &
        splflo, &
        syseqn,taubmp,tauflo,tolnew,wquad,xbl,xbord,xbr,xc,xquad,xsin,xsiq, &
        ybl,ybord,ybr,yc,yquad)

      if ( ierror/= 0 ) then
        write ( *, * ) ' '
        write ( *, * ) 'GETGRD - Fatal error!'
        write ( *, * ) '  FLOSOL returns IERROR = ',ierror
        lgdif = .false.
        call lmemry('set','fd_grad',lgdif)
        return
      end if
!
!  Set the estimate of dG/dP(Ipar)
!
      do ieqn = 1,neqn
        gdif(ieqn,ipar) = ( g1(ieqn) - g(ieqn) ) / h
      end do
!
!  If the direct finite difference estimate of the cost gradient
!  is desired, then compute the COST of the perturbed solution
!  and evaluate the finite difference quotient.
!
      if ( idfd/= 0 ) then

        call get_cost(cost1,costb,costp,costu,costv,g1,gtar,indx,ishapb,neqn, &
          np,nparb,nprof,ny,splbmp,taubmp,wateb,watep,wateu,watev,xbl,xbr, &
          ybl,ybr,yc)

        dpara3(ipar) = (cost1-cost)/h

      end if
!
!  Restore the value of the parameter.
!
      para(ipar) = psave
!
!  If the parameter is not varied, set dG/dP(Ipar) and DPARA3 to zero.
!
    else

      gdif(1:neqn,ipar) = 0.0D+00

      dpara3(ipar) = 0.0D+00

    end if

  end do
!
!  Restore the old data.
!
  call flosol ( a,area,disjac,eqn,etan,etaq,flarea,g,ierror,igrid,ijac, &
    indx,ipivot,ishapb,ishapf,isotri,iwrite,jjac,maxnew,nelem,neqn,nlband, &
    node,np,npar,nparb,nparf,npe,nrow,nx,ny,para,parjac,phi,res,splbmp,splflo, &
    syseqn,taubmp,tauflo,tolnew,wquad,xbl,xbord,xbr,xc,xquad,xsin,xsiq,ybl, &
    ybord,ybr,yc,yquad)

  if ( ierror/= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'GETGRD - Fatal error!'
    write ( *, * ) '  FLOSOL returns IERROR = ',ierror
    lgdif = .false.
    call lmemry('set','fd_grad',lgdif)
    return
  end if
!
!  Turn off the flag that warns other routines we are doing a
!  finite difference computation.
!
  lgdif = .false.
  call lmemry('set','fd_grad',lgdif)

  return
end
subroutine getsen(a,area,dudyn,dvdyn,eqn,g,ibc,indx,iopt,ipivot,ishapb, &
  ishapf,nelem,neqn,nlband,node,np,npar,nparb,nparf,npe,nrow,phi,sens,splbmp, &
  splflo,taubmp,tauflo,xc,yc)

!*****************************************************************************80
!
!! GETSEN computes the sensitivities of the state variables.
!
!  Discussion:
!
!    This routine computes the sensitivities of the state variables
!    U, V and P with respect to the parameters.
!
!    It assumes that the jacobian matrix A has already been factored.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A(NROW,NEQN), contains the
!    value of D F(I)/D X(J) for each of the NEQN residual
!    functions F(I) with respect to each of the unknown
!    coefficients X(J).
!
!    Output, real ( kind = 8 ) SENS(MAXEQN,NPAR).
!    SENS(I,J) contains the sensitivity of the I-th unknown
!    with respect to the J-th parameter.
!
  implicit none

  integer ( kind = 4 ) nelem
  integer ( kind = 4 ) neqn
  integer ( kind = 4 ) np
  integer ( kind = 4 ) npar
  integer ( kind = 4 ) nparb
  integer ( kind = 4 ) nparf
  integer ( kind = 4 ) npe
  integer ( kind = 4 ) nrow

  real ( kind = 8 ) a(nrow,neqn)
  real ( kind = 8 ) area(nelem,3)
  real ( kind = 8 ) dudyn(np)
  real ( kind = 8 ) dvdyn(np)
  character ( len = 2 ) eqn(neqn)
  real ( kind = 8 ) g(neqn)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibc
  integer ( kind = 4 ) indx(np,3)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) iopt(npar)
  integer ( kind = 4 ) ipar
  integer ( kind = 4 ) ipivot(neqn)
  integer ( kind = 4 ) ishapb
  integer ( kind = 4 ) ishapf
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) nlband
  integer ( kind = 4 ) node(nelem,npe)
  real ( kind = 8 ) phi(nelem,3,6,10)
  real ( kind = 8 ) sens(neqn,npar)
  real ( kind = 8 ) splbmp(4,nparb+2,0:nparb)
  real ( kind = 8 ) splflo(4,nparf+2,0:nparf)
  real ( kind = 8 ) taubmp(nparb+2)
  real ( kind = 8 ) tauflo(nparf+2)
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) yc(np)

  ival = 1
  call imemry('inc','GetSen_calls',ival)
!
!  Set up the right hand sides.
!
  do ipar = 1,npar

    if ( iopt(ipar) == 1 ) then

      if ( 1 <= ipar .and. ipar<=nparf ) then

        call flosen ( eqn,sens(1,ipar),indx,ipar,ishapf,neqn,np,nparf,splflo, &
          tauflo,yc)

      else if ( nparf+1 <= ipar .and. ipar<=nparf+nparb ) then

        call bump_sen ( dudyn,dvdyn,eqn,sens(1,ipar),g,ibc,indx,ipar,ishapb,neqn, &
          np,nparb,nparf,splbmp,taubmp,xc,yc)

      else if ( ipar == nparf+nparb+1 ) then

        call nusen ( area,eqn,sens(1,ipar),g,indx,nelem,neqn,node,np,npe,phi)

      end if

    else

      do i = 1,neqn
        sens(i,ipar) = 0.0D+00
      end do

    end if

  end do
!
!  Solve the linear systems.
!
  call sgbtrs('N',neqn,nlband,nlband,npar,a,nrow,ipivot,sens,neqn,info)

  ival = 1
  call imemry('inc','Solve_calls',ival)
  ival = npar
  call imemry('inc','Solve_sys',ival)

  if ( info/= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'GetSen - Warning!'
    write ( *, * ) '  Failure solving for sensitivity ',ipar
    write ( *, * ) '  SGBTRS returns nonzero INFO = ',info
  end if

  return
end
subroutine gquad1 ( nquad1, wquad1, xquad1 )

!*****************************************************************************80
!
!! GQUAD1 returns a one-dimensional Gauss quadrature rule.
!
!  Discussion:
!
!    The integral of a function F(X) over the interval [-1,1]
!
!      Integral (-1 to 1) F(X) DX
!
!    may then be approximated by
!
!      Sum (I = 1 to NQUAD1) WQUAD1(I) * F(XQUAD1(I))
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NQUAD1.
!    The user specifies the rule desired by setting NQUAD1
!    to 3 or 5.  Any other value is illegal, and will cause
!    GQUAD1 to stop.
!
!    Output, real ( kind = 8 ) WQUAD1(NQUAD1).
!    WQUAD1(I) is the weight factor corresponding to the
!    I-th quadrature point.
!
!    Output, real ( kind = 8 ) XQUAD1(NQUAD1).
!    XQUAD1(I) is the I-th quadrature point.
!
  implicit none

  integer ( kind = 4 ) nquad1

  real ( kind = 8 ) wquad1(nquad1)
  real ( kind = 8 ) xquad1(nquad1)

  if ( nquad1 == 3 ) then

    xquad1(1) = -0.7745966692D+00
    xquad1(2) =  0.0000000000D+00
    xquad1(3) =  0.7745966692D+00

    wquad1(1) = 5.0D+00 / 9.0D+00
    wquad1(2) = 8.0D+00 / 9.0D+00
    wquad1(3) = 5.0D+00 / 9.0D+00

  else if ( nquad1 == 5 ) then

    xquad1(1) = -0.906179845938664D+00
    xquad1(2) = -0.538469310105683D+00
    xquad1(3) =  0.000000000000000D+00
    xquad1(4) =  0.538469310105683D+00
    xquad1(5) =  0.906179845938664D+00

    wquad1(1) = 0.236926885056189D+00
    wquad1(2) = 0.478628670499366D+00
    wquad1(3) = 0.568888888888889D+00
    wquad1(4) = 0.478628670499366D+00
    wquad1(5) = 0.236926885056189D+00

  else

    write ( *, * ) ' '
    write ( *, * ) 'GQUAD1 - Fatal error!'
    write ( *, * ) '  An illegal value of NQUAD1 was input.'
    write ( *, * ) '  Only NQUAD1 = 3 or 5 are legal.'
    write ( *, * ) '  The input value was ',nquad1
    write ( *, * ) '  The code is stopping now.'
    stop

  end if

  return
end
function ilaenv(ispec,name,opts,n1,n2,n3,n4)

!*****************************************************************************80
!
!! ILAENV sets LAPACK environmental variables.
!
!  Discussion:
!
!    ILAENV is called from the LAPACK routines to choose problem-dependent
!    parameters for the local environment.  See ISPEC for a description of
!    the parameters.
!
!    This version provides a set of parameters which should give good,
!    but not optimal, performance on many of the currently available
!    computers.  Users are encouraged to modify this subroutine to set
!    the tuning parameters for their particular machine using the option
!    and problem size information in the arguments.
!
!    This routine will not function correctly if it is converted to all
!    lower case.  Converting it to all upper case is allowed.
!
!  Parameters:
!
!  ISPEC   (input) INTEGER
!          Specifies the parameter to be returned as the value of
!          ILAENV.
!          = 1: the optimal blocksize; if this value is 1, an unblocked
!               algorithm will give the best performance.
!          = 2: the minimum block size for which the block routine
!               should be used; if the usable block size is less than
!               this value, an unblocked routine should be used.
!          = 3: the crossover point (in a block routine, for N less
!               than this value, an unblocked routine should be used)
!          = 4: the number of shifts, used in the nonsymmetric
!               eigenvalue routines
!          = 5: the minimum column dimension for blocking to be used;
!               rectangular blocks must have dimension at least k by m,
!               where k is given by ILAENV(2,...) and m by ILAENV(5,...)
!          = 6: the crossover point for the SVD (when reducing an m by n
!               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
!               this value, a QR factorization is used first to reduce
!               the matrix to a triangular form.)
!          = 7: the number of processors
!          = 8: the crossover point for the multishift QR and QZ methods
!               for nonsymmetric eigenvalue problems.
!
!  NAME    (input) CHARACTER*(*)
!          The name of the calling subroutine, in either upper case or
!          lower case.
!
!  OPTS    (input) CHARACTER*(*)
!          The character options to the subroutine NAME, concatenated
!          into a single character string.  For example, UPLO = 'U',
!          TRANS = 'T', and DIAG = 'N' for a triangular routine would
!          be specified as OPTS = 'UTN'.
!
!  N1      (input) INTEGER
!  N2      (input) INTEGER
!  N3      (input) INTEGER
!  N4      (input) INTEGER
!          Problem dimensions for the subroutine NAME; these may not all
!          be required.
!
! (ILAENV) (output) INTEGER
!          >= 0: the value of the parameter specified by ISPEC
!          < 0:  if ILAENV = -k, the k-th argument had an illegal value.
!
!  Further Details
!   == =============
!
!  The following conventions have been used when calling ILAENV from the
!  LAPACK routines:
!  1)  OPTS is a concatenation of all of the character options to
!    subroutine NAME, in the same order that they appear in the
!      argument list for NAME, even if they are not used in determining
!      the value of the parameter specified by ISPEC.
!  2)  The problem dimensions N1, N2, N3, N4 are specified in the order
!      that they appear in the argument list for NAME.  N1 is used
!      first, N2 second, and so on, and unused problem dimensions are
!      passed a value of -1.
!  3)  The parameter value returned by ILAENV is checked for validity in
!      the calling subroutine.  For example, ILAENV is used to retrieve
!      the optimal blocksize for STRTRI as follows:
!
!      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
!      IF( NB.LE.1 ) NB = MAX( 1, N )
!
  implicit none

  character ( len = * )    name, opts
  integer ( kind = 4 ) ilaenv
  integer ( kind = 4 )            ispec, n1, n2, n3, n4
!     ..
  logical            cname, sname
  character        c1
  character ( len = 2 )        c2, c4
  character ( len = 3 )        c3
  character ( len = 6 )        subnam
  integer ( kind = 4 )            i, ic, iz, nb, nbmin, nx
!
!     ..
!     .. Executable Statements ..
!
  go to ( 100, 100, 100, 400, 500, 600, 700, 800 ) ispec
!
!     Invalid value for ISPEC
!
  ilaenv = -1
  return
!
  100 continue
!
!     Convert NAME to upper case if the first character is lower case.
!
  ilaenv = 1
  subnam = name
  ic = ichar( subnam( 1:1 ) )
  iz = ichar( 'z' )
  if (  iz == 90 .or. iz==122 ) then
!
!        ASCII character set
!
     if (  ic >= 97 .and. ic<=122 ) then
        subnam( 1:1 ) = char( ic-32 )
        do 10 i = 2, 6
           ic = ichar( subnam( i:i ) )
           if (  ic >= 97 .and. ic<=122 )subnam( i:i ) = char( ic-32 )
   10       continue
     end if
!
  else if (  iz == 233 .or. iz==169 ) then
!
!        EBCDIC character set
!
     if (  ( ic >= 129 .and. ic<=137 ) .or. &
            ( ic >= 145 .and. ic<=153 ) .or. &
            ( ic >= 162 .and. ic<=169 ) ) then
        subnam( 1:1 ) = char( ic+64 )
        do 20 i = 2, 6
           ic = ichar( subnam( i:i ) )
           if (  ( ic >= 129 .and. ic<=137 ) .or. &
                  ( ic >= 145 .and. ic<=153 ) .or. &
                  ( ic >= 162 .and. ic<=169 ) ) &
                 subnam( i:i ) = char( ic+64 )
   20       continue
     end if
!
  else if (  iz == 218 .or. iz==250 ) then
!
!        Prime machines:  ASCII+128
!
     if (  ic >= 225 .and. ic<=250 ) then
        subnam( 1:1 ) = char( ic-32 )
        do 30 i = 2, 6
           ic = ichar( subnam( i:i ) )
           if (  ic >= 225 .and. ic<=250 )subnam( i:i ) = char( ic-32 )
   30       continue
     end if
  end if

  c1 = subnam( 1:1 )
  sname = c1 == 's' .or. c1=='d'
  cname = c1 == 'c' .or. c1=='z'
  if (  .not.( cname .or. sname ) )return
  c2 = subnam( 2:3 )
  c3 = subnam( 4:6 )
  c4 = c3( 2:3 )

  go to ( 110, 200, 300 ) ispec

  110 continue
!
!     ISPEC = 1:  block size
!
!     In these examples, separate code is provided for setting NB for
!     real and complex.  We assume that NB will take the same value in
!     single or real.
!
  nb = 1
!
  if (  c2 == 'ge' ) then
     if (  c3 == 'trf' ) then
        if (  sname ) then
           nb = 64
        else
           nb = 64
        end if
     else if (  c3 == 'qrf' .or. c3=='rqf' .or. c3=='lqf' .or.c3=='qlf' ) then
        if (  sname ) then
           nb = 32
        else
           nb = 32
        end if
     else if (  c3 == 'hrd' ) then
        if (  sname ) then
           nb = 32
        else
           nb = 32
        end if
     else if (  c3 == 'brd' ) then
        if (  sname ) then
           nb = 32
        else
           nb = 32
        end if
     else if (  c3 == 'tri' ) then
        if (  sname ) then
           nb = 64
        else
           nb = 64
        end if
     end if
  else if (  c2 == 'po' ) then
     if (  c3 == 'trf' ) then
        if (  sname ) then
           nb = 64
        else
           nb = 64
        end if
     end if
  else if (  c2 == 'sy' ) then
     if (  c3 == 'trf' ) then
        if (  sname ) then
           nb = 64
        else
           nb = 64
        end if
     else if (  sname .and. c3 == 'trd' ) then
        nb = 1
     else if (  sname .and. c3 == 'gst' ) then
        nb = 64
     end if
  else if (  cname .and. c2 == 'he' ) then
     if (  c3 == 'trf' ) then
        nb = 64
     else if (  c3 == 'trd' ) then
        nb = 1
     else if (  c3 == 'gst' ) then
        nb = 64
     end if
  else if (  sname .and. c2 == 'or' ) then
     if (  c3( 1:1 ) == 'g' ) then
        if (  c4 == 'qr' .or. c4=='rq' .or. c4=='lq' .or. &
               c4 == 'ql' .or. c4=='hr' .or. c4=='tr' .or. &
               c4 == 'br' ) then
           nb = 32
        end if
     else if (  c3( 1:1 ) == 'm' ) then
        if (  c4 == 'qr' .or. c4=='rq' .or. c4=='lq' .or. &
               c4 == 'ql' .or. c4=='hr' .or. c4=='tr' .or. &
               c4 == 'br' ) then
           nb = 32
        end if
     end if
  else if (  cname .and. c2 == 'un' ) then
     if (  c3( 1:1 ) == 'g' ) then
        if (  c4 == 'qr' .or. c4=='rq' .or. c4=='lq' .or. &
               c4 == 'ql' .or. c4=='hr' .or. c4=='tr' .or. &
              c4 == 'br' ) then
           nb = 32
        end if
     else if (  c3( 1:1 ) == 'm' ) then
        if (  c4 == 'qr' .or. c4=='rq' .or. c4=='lq' .or. &
               c4 == 'ql' .or. c4=='hr' .or. c4=='tr' .or. &
               c4 == 'br' ) then
           nb = 32
        end if
     end if
  else if (  c2 == 'gb' ) then
     if (  c3 == 'trf' ) then
        if (  sname ) then
           if (  n4 <= 64 ) then
              nb = 1
           else
              nb = 32
           end if
        else
           if (  n4 <= 64 ) then
              nb = 1
           else
              nb = 32
           end if
        end if
     end if
  else if (  c2 == 'pb' ) then
     if (  c3 == 'trf' ) then
        if (  sname ) then
           if (  n2 <= 64 ) then
              nb = 1
           else
              nb = 32
           end if
        else
           if (  n2 <= 64 ) then
              nb = 1
           else
              nb = 32
           end if
        end if
     end if
  else if (  c2 == 'tr' ) then
     if (  c3 == 'tri' ) then
        if (  sname ) then
           nb = 64
        else
           nb = 64
        end if
     end if
  else if (  c2 == 'la' ) then
     if (  c3 == 'uum' ) then
        if (  sname ) then
           nb = 64
        else
           nb = 64
        end if
     end if
  else if (  sname .and. c2 == 'st' ) then
     if (  c3 == 'ebz' ) then
        nb = 1
     end if
  end if
  ilaenv = nb
  return
!
  200 continue
!
!     ISPEC = 2:  minimum block size
!
  nbmin = 2
  if (  c2 == 'ge' ) then
     if (  c3 == 'qrf' .or. c3=='rqf' .or. c3=='lqf' .or. &
           c3 == 'qlf' ) then
        if (  sname ) then
           nbmin = 2
        else
           nbmin = 2
        end if
     else if (  c3 == 'hrd' ) then
        if (  sname ) then
           nbmin = 2
        else
           nbmin = 2
        end if
     else if (  c3 == 'brd' ) then
        if (  sname ) then
           nbmin = 2
        else
           nbmin = 2
        end if
     else if (  c3 == 'tri' ) then
        if (  sname ) then
           nbmin = 2
        else
           nbmin = 2
        end if
     end if
  else if (  c2 == 'sy' ) then
     if (  c3 == 'trf' ) then
        if (  sname ) then
           nbmin = 2
        else
           nbmin = 2
        end if
     else if (  sname .and. c3 == 'trd' ) then
        nbmin = 2
     end if
  else if (  cname .and. c2 == 'he' ) then
     if (  c3 == 'trd' ) then
        nbmin = 2
     end if
  else if (  sname .and. c2 == 'or' ) then
     if (  c3( 1:1 ) == 'g' ) then
        if (  c4 == 'qr' .or. c4=='rq' .or. c4=='lq' .or. &
              c4 == 'ql' .or. c4=='hr' .or. c4=='tr' .or. &
               c4 == 'br' ) then
           nbmin = 2
        end if
     else if (  c3( 1:1 ) == 'm' ) then
        if (  c4 == 'qr' .or. c4=='rq' .or. c4=='lq' .or. &
              c4 == 'ql' .or. c4=='hr' .or. c4=='tr' .or. &
              c4 == 'br' ) then
           nbmin = 2
        end if
     end if
  else if (  cname .and. c2 == 'un' ) then
     if (  c3( 1:1 ) == 'g' ) then
        if (  c4 == 'qr' .or. c4=='rq' .or. c4=='lq' .or. &
              c4 == 'ql' .or. c4=='hr' .or. c4=='tr' .or. &
              c4 == 'br' ) then
           nbmin = 2
        end if
     else if (  c3( 1:1 ) == 'm' ) then
        if (  c4 == 'qr' .or. c4=='rq' .or. c4=='lq' .or. &
              c4 == 'ql' .or. c4=='hr' .or. c4=='tr' .or. &
              c4 == 'br' ) then
           nbmin = 2
        end if
     end if
  end if
  ilaenv = nbmin
  return
!
  300 continue
!
!     ISPEC = 3:  crossover point
!
  nx = 0
  if (  c2 == 'ge' ) then
     if (  c3 == 'qrf' .or. c3=='rqf' .or. c3=='lqf' .or.c3=='qlf' ) then
        if (  sname ) then
           nx = 128
        else
           nx = 128
        end if
     else if (  c3 == 'hrd' ) then
        if (  sname ) then
           nx = 128
        else
           nx = 128
        end if
     else if (  c3 == 'brd' ) then
        if (  sname ) then
           nx = 128
        else
           nx = 128
        end if
     end if
  else if (  c2 == 'sy' ) then
     if (  sname .and. c3 == 'trd' ) then
        nx = 1
     end if
  else if (  cname .and. c2 == 'he' ) then
     if (  c3 == 'trd' ) then
        nx = 1
     end if
  else if (  sname .and. c2 == 'or' ) then
     if (  c3( 1:1 ) == 'g' ) then
        if (  c4 == 'qr' .or. c4=='rq' .or. c4=='lq' .or. &
              c4 == 'ql' .or. c4=='hr' .or. c4=='tr' .or. &
              c4 == 'br' ) then
           nx = 128
        end if
     end if
  else if (  cname .and. c2 == 'un' ) then
     if (  c3( 1:1 ) == 'g' ) then
        if (  c4 == 'qr' .or. c4=='rq' .or. c4=='lq' .or. &
              c4 == 'ql' .or. c4=='hr' .or. c4=='tr' .or. &
             c4 == 'br' ) then
           nx = 128
        end if
     end if
  end if
  ilaenv = nx
  return
!
  400 continue
!
!     ISPEC = 4:  number of shifts (used by xHSEQR)
!
  ilaenv = 6
  return
!
  500 continue
!
!     ISPEC = 5:  minimum column dimension (not used)
!
  ilaenv = 2
  return
!
  600 continue
!
!     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD)
!
  ilaenv = int( real ( min ( n1, n2 ), kind = 8 ) * 1.6D+00 )

  return
!
  700 continue
!
!     ISPEC = 7:  number of processors (not used)
!
  ilaenv = 1
  return
!
  800 continue
!
!     ISPEC = 8:  crossover point for multishift (used by xHSEQR)
!
  ilaenv = 50
  return
end
subroutine imemry(action,name,ival)

!*****************************************************************************80
!
!! IMEMRY sets or gets integer values from an internal memory.
!
!  Discussion:
!
!    The routine allows the user to define the name of an integer variable,
!    set it, increment it, or get the value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, Character*(*) ACTION, desired action.
!    'Init', reset all values to zero, wipe out all names.
!    'Name', add a variable of the given name.
!    'Inc',  increment variable NAME by IVAL.
!    'Set',  set variable NAME to IVAL.
!    'Get',  return value of NAME in IVAL.
!    'Zero', reset all values to zero.
!
!    Input, Character*(*) NAME, the name of the variable.
!
!    Input/output, Integer IVAL.
!    For the 'Inc' and 'Set' commands, IVAL must contain the
!    increment or set value.
!    For the 'Get' command, IVAL will contain the value of the
!    named variable on output.
!
  implicit none

  integer ( kind = 4 ) maxnam
  parameter (maxnam = 100)

  character ( len = * ) action
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) ivals(maxnam)
  logical s_eqi
  character ( len = * ) name
  character ( len = 20 ) names(maxnam)
  integer ( kind = 4 ) numnam

  save ivals
  save names
  save numnam
!
!  Initialize everything.
!
  if ( s_eqi(action,'init') ) then

    numnam = 0

    do i = 1,maxnam
      ivals(i) = 0
      names(i) = ' '
    end do
!
!  Name something.
!
  else if ( s_eqi(action,'name') ) then

    do i = 1,numnam

      if ( s_eqi(name,names(i)) ) then
        write ( *, * ) ' '
        write ( *, * ) 'IMemry - Warning!'
        write(*,'(''  There is ALREADY a variable '',a)') trim ( name )
        return
      end if

    end do

    if ( numnam<maxnam ) then
      numnam = numnam+1
      names(numnam) = name
    else
    write ( *, * ) ' '
      write ( *, * ) 'IMemry - Fatal error!'
      write ( *, * ) '  No more name space.'
      stop
    end if
!
!  Increment something.
!
  else if ( s_eqi(action,'inc') ) then

    do i = 1,numnam

      if ( s_eqi(name,names(i)) ) then
        ivals(i) = ivals(i)+ival
        return
      end if

    end do

    write ( *, * ) ' '
    write ( *, * ) 'IMemry - Fatal error!'
    write ( *, * ) '  Attempt to increment unknown variable.'
    write(*,'(''  Variable name is '',a)') trim ( name )
    stop
!
!  Set something.
!
  else if ( s_eqi(action,'set') ) then

    do i = 1,numnam

      if ( s_eqi(name,names(i)) ) then
        ivals(i) = ival
        return
      end if

    end do

    write ( *, * ) ' '
    write ( *, * ) 'IMemry - Fatal error!'
    write ( *, * ) '  Attempt to set unknown variable.'
    write(*,'(''  Variable name is '',a)') trim ( name )
    stop
!
!  Get something.
!
  else if ( s_eqi(action,'get') ) then

    do i = 1,numnam

      if ( s_eqi(name,names(i)) ) then
        ival = ivals(i)
        return
      end if

    end do

    write ( *, * ) ' '
    write ( *, * ) 'IMemry - Fatal error!'
    write ( *, * ) '  Attempt to get value of unknown variable.'
    write(*,'(''  Variable name is '',a)') trim ( name )
    stop
!
!  Initialize everything.
!
  else if ( s_eqi(action,'zero') ) then

    do i = 1,numnam
      ivals(i) = 0
    end do

    write ( *, * ) ' '
    write ( *, * ) 'IMemry - Note:'
    write ( *, * ) '  All data has been reset to zero.'
!
!  Unknown action.
!
  else

    write ( *, * ) ' '
    write ( *, * ) 'IMemry - Fatal error!'
    write ( *, * ) '  Unrecognized action requested.'
    write(*,'(1x,a)')action
    stop

  end if

  return
end
subroutine init(area,cost,costar,costb,costp,costu,costv,disjac,dopt,dpara3, &
  dparsn,dparfd,dparfdc,dpdyn,dudyn,dvdyn,dydpn,eqn,etan,plot_file,march_file, &
  g,gdif,gold,gopt,gradf,gtar,ibc,ibump,idfd,ids,ierror,ifds,igrad, &
  igrid,igunit,ijac,indx,iopt,iplot,ipred,ishapb,ishapbt,ishapf,ishapft, &
  ismooth,isotri,istep1,istep2,itar,itunit,itype,ivopt,iwrite,jjac,jstep1, &
  jstep2,liv,lv,maxelm,maxeqn,maxnew,maxnp,maxnx,maxny,maxpar,maxparb, &
  maxstp,node,nopt,npar,nparb,nparf,npe,nstep3,nx,ny,para1,para2,para3,parjac, &
  partar,sens,stpmax,syseqn,tolnew,tolopt,vopt,wateb,wateb1,wateb2,watep, &
  wateu,watev,xbleft,xbltar,xbord,xbrite,xbrtar,xprof,xsin,ybleft,ybltar, &
  ybord,ybrite,ybrtar)

!*****************************************************************************80
!
!! INIT initializes the program parameters.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) liv
  integer ( kind = 4 ) lv
  integer ( kind = 4 ) maxelm
  integer ( kind = 4 ) maxeqn
  integer ( kind = 4 ) maxnp
  integer ( kind = 4 ) maxnx
  integer ( kind = 4 ) maxny
  integer ( kind = 4 ) maxpar
  integer ( kind = 4 ) maxparb
  integer ( kind = 4 ) npe

  real ( kind = 8 ) area(maxelm,3)
  real ( kind = 8 ) cost
  real ( kind = 8 ) costar
  real ( kind = 8 ) costb
  real ( kind = 8 ) costp
  real ( kind = 8 ) costu
  real ( kind = 8 ) costv
  real ( kind = 8 ) disjac
  real ( kind = 8 ) dopt(maxpar)
  real ( kind = 8 ) dpara3(maxpar)
  real ( kind = 8 ) dparsn(maxpar)
  real ( kind = 8 ) dparfd(maxpar)
  real ( kind = 8 ) dparfdc(maxpar)
  real ( kind = 8 ) dpdyn(maxnp)
  real ( kind = 8 ) dudyn(maxnp)
  real ( kind = 8 ) dval
  real ( kind = 8 ) dvdyn(maxnp)
  real ( kind = 8 ) dydpn(maxnp,maxparb)
  character ( len = 2 ) eqn(maxeqn)
  real ( kind = 8 ) etan(6)
  character ( len = * ) plot_file
  character ( len = * ) march_file
  real ( kind = 8 ) g(maxeqn)
  real ( kind = 8 ) gdif(maxeqn,maxpar)
  real ( kind = 8 ) gold(maxeqn)
  real ( kind = 8 ) gopt(maxpar)
  real ( kind = 8 ) gradf(maxeqn,maxpar)
  real ( kind = 8 ) gtar(maxeqn)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibc
  integer ( kind = 4 ) ibump
  integer ( kind = 4 ) idfd
  integer ( kind = 4 ) ids
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ifds
  integer ( kind = 4 ) igrad
  integer ( kind = 4 ) igrid
  integer ( kind = 4 ) igunit
  integer ( kind = 4 ) ijac
  integer ( kind = 4 ) indx(maxnp,3)
  integer ( kind = 4 ) iopt(maxpar)
  integer ( kind = 4 ) iplot
  integer ( kind = 4 ) ipred
  integer ( kind = 4 ) ishapb
  integer ( kind = 4 ) ishapbt
  integer ( kind = 4 ) ishapf
  integer ( kind = 4 ) ishapft
  integer ( kind = 4 ) ismooth
  integer ( kind = 4 ) isotri(maxelm)
  integer ( kind = 4 ) istep1
  integer ( kind = 4 ) istep2
  integer ( kind = 4 ) itar
  integer ( kind = 4 ) itunit
  integer ( kind = 4 ) itype
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) ivopt(liv)
  integer ( kind = 4 ) iwrite
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jjac
  integer ( kind = 4 ) jstep1
  integer ( kind = 4 ) jstep2
  logical lval
  integer ( kind = 4 ) maxnew
  integer ( kind = 4 ) maxstp
  integer ( kind = 4 ) node(maxelm,npe)
  integer ( kind = 4 ) nopt
  integer ( kind = 4 ) npar
  integer ( kind = 4 ) nparb
  integer ( kind = 4 ) nparf
  integer ( kind = 4 ) nstep3
  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny
  real ( kind = 8 ) para1(maxpar)
  real ( kind = 8 ) para2(maxpar)
  real ( kind = 8 ) para3(maxpar)
  real ( kind = 8 ) parjac(maxpar)
  real ( kind = 8 ) partar(maxpar)
  real ( kind = 8 ) sens(maxeqn,maxpar)
  real ( kind = 8 ) stpmax
  character ( len = 20 ) syseqn
  real ( kind = 8 ) tolnew
  real ( kind = 8 ) tolopt
  real ( kind = 8 ) vopt(lv)
  real ( kind = 8 ) wateb
  real ( kind = 8 ) wateb1
  real ( kind = 8 ) wateb2
  real ( kind = 8 ) watep
  real ( kind = 8 ) wateu
  real ( kind = 8 ) watev
  real ( kind = 8 ) xbleft
  real ( kind = 8 ) xbltar
  real ( kind = 8 ) xbord(maxnx)
  real ( kind = 8 ) xbrite
  real ( kind = 8 ) xbrtar
  real ( kind = 8 ) xprof
  real ( kind = 8 ) xsin(6)
  real ( kind = 8 ) ybleft
  real ( kind = 8 ) ybltar
  real ( kind = 8 ) ybord(maxny)
  real ( kind = 8 ) ybrite
  real ( kind = 8 ) ybrtar
!
!  Set values of some counters saved in DMEMRY, IMEMRY and LMEMRY.
!
  dval = 0.0D+00
  call dmemry('init','nothing',dval)

  ival = 0
  call imemry('init','nothing',ival)

  call imemry('name','CubSpl_calls',ival)
  call imemry('name','Factor_calls',ival)
  call imemry('name','FloSol_calls',ival)
  call imemry('name','Fprime_calls',ival)
  call imemry('name','Fx_calls',ival)
  call imemry('name','GET_COST_calls',ival)
  call imemry('name','COST_GRADIENT_calls',ival)
  call imemry('name','GetFix_calls',ival)
  call imemry('name','GetGrd_calls',ival)
  call imemry('name','GetSen_calls',ival)
  call imemry('name','Newton_calls',ival)
  call imemry('name','Newton_damp',ival)
  call imemry('name','Newton_fail',ival)
  call imemry('name','Newton_falter',ival)
  call imemry('name','Newton_steps',ival)
  call imemry('name','Newton_zero',ival)
  call imemry('name','PltWrt_calls',ival)
  call imemry('name','Points_Con',ival)
  call imemry('name','Points_FD',ival)
  call imemry('name','Points_March',ival)
  call imemry('name','Points_Opt',ival)
  call imemry('name','Points_Tar',ival)
  call imemry('name','Restarts',ival)
  call imemry('name','SetBas_calls',ival)
  call imemry('name','NODE_SET_calls',ival)
  call imemry('name','SetQXY_calls',ival)
  call imemry('name','XY_SET_calls',ival)
  call imemry('name','Snoit_calls',ival)
  call imemry('name','SolCon_calls',ival)
  call imemry('name','SolCon_steps',ival)
  call imemry('name','Solve_calls',ival)
  call imemry('name','Solve_sys',ival)
  call imemry('name','Sumit_calls',ival)
  call imemry('name','Trans_calls',ival)
  call imemry('name','UVal_calls',ival)
  call imemry('name','UValQ_calls',ival)
  call imemry('name','UpValQ_calls',ival)
  call imemry('name','Xofxsi_calls',ival)
!
!  Initialize the logical flags.
!
  lval = .false.
  call lmemry('init','nothing',lval)
!
!  Add the names of individual logical flags.
!
  lval = .false.
  call lmemry('name','bcfile',lval)
  call lmemry('name','fd_grad',lval)
  call lmemry('name','have_fp',lval)
  call lmemry('name','need_phi',lval)
  call lmemry('name','need_xy',lval)
  call lmemry('name','target',lval)
!
!  Assign initial values to the flags.
!
  lval = .false.
  call lmemry('set','bcfile',lval)
  lval = .false.
  call lmemry('set','fd_grad',lval)
  lval = .false.
  call lmemry('set','have_fp',lval)
  lval = .true.
  call lmemry('set','need_phi',lval)
  lval = .true.
  call lmemry('set','need_xy',lval)
  lval = .false.
  call lmemry('set','target',lval)
!
!  Set plain old variables.
!
  area(1:maxelm,1:3) = 0.0D+00
  cost = 0.0D+00
  costar = 0.0D+00
  costb = 0.0D+00
  costp = 0.0D+00
  costu = 0.0D+00
  costv = 0.0D+00
  disjac = 0.0D+00
  dopt(1:maxpar) = 1.0D+00
  dpara3(1:maxpar) = 0.0D+00
  dparfd(1:maxpar) = 0.0D+00
  dparfdc(1:maxpar) = 0.0D+00
  dparsn(1:maxpar) = 0.0D+00
  dpdyn(1:maxnp) = 0.0D+00
  dudyn(1:maxnp) = 0.0D+00
  dvdyn(1:maxnp) = 0.0D+00
  dydpn(1:maxnp,1:maxparb) = 0.0D+00
  eqn(1:maxeqn) = '??'
  etan(1:6) = (/ 0.0D+00, 1.0D+00, 0.0D+00, 0.5D+00, 0.5D+00, 0.0D+00 /)
  plot_file = 'display.dat'
  march_file = 'march.txt'
  g(1:maxeqn) = 0.0D+00
  gold(1:maxeqn) = 0.0D+00
  gopt(1:maxpar) = 0.0D+00
  gdif(1:maxeqn,1:maxpar) = 0.0D+00
  gradf(1:maxeqn,1:maxpar) = 0.0D+00
  gtar(1:maxeqn) = 0.0D+00
  ibc = 0
  ibump = 2
  idfd = 1
  ids = 1
  ierror = 0
  ifds = 1
  igrad = 1
  igrid = 0
  igunit = 0
  ijac = 1
  indx(1:maxnp,1:3) = 0
  iopt(1:maxpar) = 0
  iplot = 0
  ipred = 1
  ishapb = 2
  ishapbt = 2
  ishapf = 2
  ishapft = 2
  ismooth = 0
  isotri(1:maxelm) = 0
  istep1 = 1
  istep2 = 1
  itar = 0
  itunit = 0
  itype = 3
  iunit = 0
  ivopt(1:liv) = 0
  iwrite = 0
  jjac = 1
  jstep1 = 1
  jstep2 = 1
  maxnew = 10
  maxstp = 10
  node(1:maxelm,1:npe) = 0
  nopt = 0
  npar = 1
  nparb = 0
  nparf = 0
  nstep3 = 0
  nx = 11
  ny = 4
  para1(1:maxpar) = 0.0D+00
  para2(1:maxpar) = 0.0D+00
  para3(1:maxpar) = 0.0D+00
  parjac(1:maxpar) = 0.0D+00
  partar(1:maxpar) = 0.0D+00
  sens(1:maxeqn,1:maxpar) = 0.0D+00
  stpmax = 1.0D+00
  syseqn = 'NavierStokes'
  tolnew = sqrt ( epsilon ( tolnew ) )
  tolopt = sqrt ( epsilon ( tolopt ) )
  vopt(1:lv) = 0.0D+00
  wateb = 1.0D+00
  wateb1 = 0.0D+00
  wateb2 = 1.0D+00
  watep = 1.0D+00
  wateu = 1.0D+00
  watev = 1.0D+00
  xbleft = 1.0D+00
  xbltar = 1.0D+00
  xbord(1:maxnx) = 0.0D+00
  xbrite = 3.0D+00
  xbrtar = 3.0D+00
  xprof = 3.0D+00
  xsin(1:6) = (/ 0.0D+00, 1.0D+00, 1.0D+00, 0.5D+00, 1.0D+00, 0.5D+00 /)
  ybleft = 0.0D+00
  ybltar = 0.0D+00
  ybord(1:maxny) = 0.0D+00
  ybrite = 0.0D+00
  ybrtar = 0.0D+00

  return
end
subroutine input(disjac,plot_file,march_file,ibc,ibump,idfd,ids,ierror,ifds, &
  igrad,igrid,ijac,iopt,iplot,ipred,ishapb,ishapbt,ishapf,ishapft,ismooth, &
  istep1,istep2,itar,itype,iwrite,jjac,jstep1,jstep2,maxnew,maxnx,maxny, &
  maxpar,maxstp,npar,nparb,nparf,nstep3,nx,ny,para1,para2,para3,partar, &
  stpmax,syseqn,tolnew,tolopt,wateb,wateb1,wateb2,watep,wateu,watev,xbleft, &
  xbltar,xbord,xbrite,xbrtar,xprof,ybleft,ybltar,ybord,ybrite,ybrtar)

!*****************************************************************************80
!
!! INPUT reads user input of the form "name = value".
!
!  Discussion:
!
!    For information on the meaning and legal values of the variables,
!    please refer to the glossary!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 December 2000
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) maxnx
  integer ( kind = 4 ) maxny
  integer ( kind = 4 ) maxpar

  real ( kind = 8 ) disjac
  character ( len = * ) plot_file
  character ( len = * ) march_file
  integer ( kind = 4 ) ibc
  integer ( kind = 4 ) ibump
  integer ( kind = 4 ) idfd
  integer ( kind = 4 ) ids
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ifds
  integer ( kind = 4 ) igrad
  integer ( kind = 4 ) igrid
  integer ( kind = 4 ) ijac
  integer ( kind = 4 ) iopt(maxpar)
  integer ( kind = 4 ) iplot
  integer ( kind = 4 ) ipred
  integer ( kind = 4 ) ishapb
  integer ( kind = 4 ) ishapbt
  integer ( kind = 4 ) ishapf
  integer ( kind = 4 ) ishapft
  integer ( kind = 4 ) ismooth
  integer ( kind = 4 ) istep1
  integer ( kind = 4 ) istep2
  integer ( kind = 4 ) itar
  integer ( kind = 4 ) itype
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) iwrite
  integer ( kind = 4 ) jjac
  integer ( kind = 4 ) jstep1
  integer ( kind = 4 ) jstep2
  integer ( kind = 4 ) lchar
  integer ( kind = 4 ) lchar2
  logical s_eqi
  integer ( kind = 4 ) maxnew
  integer ( kind = 4 ) maxstp
  character ( len = 80 ) name
  integer ( kind = 4 ) npar
  integer ( kind = 4 ) nparb
  integer ( kind = 4 ) nparf
  integer ( kind = 4 ) nrec
  integer ( kind = 4 ) nstep3
  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny
  real ( kind = 8 ) para1(maxpar)
  real ( kind = 8 ) para2(maxpar)
  real ( kind = 8 ) para3(maxpar)
  real ( kind = 8 ) partar(maxpar)
  character ( len = 80 ) rhs
  real ( kind = 8 ) stpmax
  character ( len = 20 ) syseqn
  real ( kind = 8 ) tolnew
  real ( kind = 8 ) tolopt
  real ( kind = 8 ) value
  real ( kind = 8 ) wateb
  real ( kind = 8 ) wateb1
  real ( kind = 8 ) wateb2
  real ( kind = 8 ) watep
  real ( kind = 8 ) wateu
  real ( kind = 8 ) watev
  real ( kind = 8 ) xbleft
  real ( kind = 8 ) xbltar
  real ( kind = 8 ) xbord(maxnx)
  real ( kind = 8 ) xbrite
  real ( kind = 8 ) xbrtar
  real ( kind = 8 ) xprof
  real ( kind = 8 ) ybleft
  real ( kind = 8 ) ybltar
  real ( kind = 8 ) ybord(maxny)
  real ( kind = 8 ) ybrite
  real ( kind = 8 ) ybrtar

  write ( *, * ) ' '
  write ( *, * ) 'INPUT:'
  write ( *, * ) '  Read control information from the user.'
  write ( *, * ) ' '

  nrec = 0
!
!  Read the next line of input.
!
  do
!
!  Decode the line, assuming it has the form "NAME = VALUE"
!
    call namels ( name, ierror, rhs, value )

    lchar = len_trim ( name )

    if ( ierror == 0 ) then
!
!  Echo the input line.  If the input quantity had an "integer"
!  name, then print it as an integer.
!
      if ( s_eqi ( name(1:10), 'march_file' ) .or.  &
           s_eqi ( name(1:9), 'plot_file' ) ) then

        write ( *, '(a,'' = '',a)' ) trim ( name ), trim ( rhs )

      else if ( s_eqi(name(1:6),'syseqn') ) then

        write ( *, '(a,'' = '',a)' ) trim ( name ), trim ( rhs )

      else if ( (lge(name(1:1),'I') .and. lle(name(1:1),'N')) .or.  &
        (lge(name(1:1),'i') .and. lle(name(1:1),'n')) ) then

        write ( *, '(a,'' = '',i14)' ) trim ( name ), int ( value )

      else

        write ( *, '(a,'' = '',g14.6)' ) trim ( name ), value

      end if

      if ( s_eqi(name,'disjac') ) then
        disjac = value
      else if ( s_eqi(name,'plot_file') ) then
        plot_file = rhs(1:30)
      else if ( s_eqi(name,'march_file') ) then
        march_file = rhs(1:30)
      else if ( s_eqi(name,'ibc') ) then
        ibc = int(value)
      else if ( s_eqi(name,'ibump') ) then
        ibump = int(value)
      else if ( s_eqi(name,'idfd') ) then
        idfd = int(value)
      else if ( s_eqi(name,'ids') ) then
        ids = int(value)
      else if ( s_eqi(name,'ifds') ) then
        ifds = int(value)
      else if ( s_eqi(name,'igrad') ) then
        igrad = int(value)
      else if ( s_eqi(name,'igrid') ) then
        igrid = int(value)
      else if ( s_eqi(name,'ijac') ) then
        ijac = int(value)
      else if ( s_eqi(name(1:5),'iopt(') ) then

        call chrcti(name(6:),ival,ierror,lchar)

        if ( ierror/= 0 ) then
          write ( *, * ) ' '
          write ( *, * ) 'INPUT - Fatal error!'
          write ( *, * ) '  ChrCTI returned nonzero error flag!'
          stop
        end if

        if ( ival<1 .or. ival>maxpar ) then
          write ( *, * ) ' '
          write ( *, * ) 'INPUT - Fatal error!'
          write ( *, * ) '  Index of IOPT is out of bounds!'
          stop
        end if

        iopt(ival) = int(value)

      else if ( s_eqi(name,'iplot') ) then
        iplot = int(value)
      else if ( s_eqi(name,'ipred') ) then
        ipred = int(value)
      else if ( s_eqi(name,'ishapb') ) then
        ishapb = int(value)
      else if ( s_eqi(name,'ishapbt') ) then
        ishapbt = int(value)
      else if ( s_eqi(name,'ishapf') ) then
        ishapf = int(value)
      else if ( s_eqi(name,'ishapft') ) then
        ishapft = int(value)
      else if ( s_eqi(name,'ismooth') ) then
        ismooth = int(value)
      else if ( s_eqi(name,'istep1') ) then
        istep1 = int(value)
      else if ( s_eqi(name,'istep2') ) then
        istep2 = int(value)
      else if ( s_eqi(name,'itar') ) then
        itar = int(value)
      else if ( s_eqi(name,'itype') ) then
        itype = int(value)
      else if ( s_eqi(name,'iwrite') ) then
        iwrite = int(value)
      else if ( s_eqi(name,'jjac') ) then
        jjac = int(value)
      else if ( s_eqi(name,'jstep1') ) then
        jstep1 = int(value)
      else if ( s_eqi(name,'jstep2') ) then
        jstep2 = int(value)
      else if ( s_eqi(name,'maxnew') ) then
        maxnew = int(value)
      else if ( s_eqi(name,'maxstp') ) then
        maxstp = int(value)
      else if ( s_eqi(name,'nparb') ) then
        nparb = int(value)
        npar = nparf+nparb+1
      else if ( s_eqi(name,'nparf') ) then
        nparf = int(value)
        npar = nparf+nparb+1
      else if ( s_eqi(name,'nx') ) then
        nx = int(value)
      else if ( s_eqi(name,'ny') ) then
        ny = int(value)
      else if ( s_eqi(name,'nstep3') ) then
        nstep3 = int(value)
      else if ( s_eqi(name(1:6),'para1(') ) then

        call chrcti(name(7:),ival,ierror,lchar)

        if ( ierror/= 0 ) then
          write ( *, * ) ' '
          write ( *, * ) 'INPUT - Fatal error!'
          write ( *, * ) '  ChrCTI returned nonzero error flag!'
          stop
        end if

        if ( ival<1 .or. ival>maxpar ) then
          write ( *, * ) ' '
          write ( *, * ) 'INPUT - Fatal error!'
          write ( *, * ) '  Index of PARA1 is out of bounds!'
          stop
        end if

        para1(ival) = value

      else if ( s_eqi(name(1:6),'para2(') ) then

        call chrcti(name(7:),ival,ierror,lchar)

        if ( ierror/= 0 ) then
          write ( *, * ) ' '
          write ( *, * ) 'INPUT - Fatal error!'
          write ( *, * ) '  ChrCTI returned nonzero error flag!'
          stop
        end if

        if ( ival<1 .or. ival>maxpar ) then
          write ( *, * ) ' '
          write ( *, * ) 'INPUT - Fatal error!'
          write ( *, * ) '  Index of PARA2 is out of bounds!'
          stop
        end if

        para2(ival) = value

      else if ( s_eqi(name(1:6),'para3(') ) then

        call chrcti(name(7:),ival,ierror,lchar)

        if ( ierror/= 0 ) then
          write ( *, * ) ' '
          write ( *, * ) 'INPUT - Fatal error!'
          write ( *, * ) '  ChrCTI returned nonzero error flag!'
          stop
        end if

        if ( ival<1 .or. ival>maxpar ) then
          write ( *, * ) ' '
          write ( *, * ) 'INPUT - Fatal error!'
          write ( *, * ) '  Index of PARA3 is out of bounds!'
          stop
        end if

        para3(ival) = value

      else if ( s_eqi(name(1:7),'partar(') ) then

        call chrcti(name(8:),ival,ierror,lchar)

        if ( ierror/= 0 ) then
          write ( *, * ) ' '
          write ( *, * ) 'INPUT - Fatal error!'
          write ( *, * ) '  ChrCTI returned nonzero error flag!'
          stop
        end if

        if ( ival<1 .or. ival>maxpar ) then
          write ( *, * ) ' '
          write ( *, * ) 'INPUT - Fatal error!'
          write ( *, * ) '  Index of PARTAR is out of bounds!'
          stop
        end if

        partar(ival) = value

      else if ( s_eqi(name,'stpmax') ) then
        stpmax = value
      else if ( s_eqi(name,'syseqn') ) then
        syseqn = rhs(1:20)
      else if ( s_eqi(name,'tolnew') ) then
        tolnew = value
      else if ( s_eqi(name,'tolopt') ) then
        tolopt = value
      else if ( s_eqi(name,'wateu') ) then
        wateu = value
      else if ( s_eqi(name,'watev') ) then
        watev = value
      else if ( s_eqi(name,'watep') ) then
        watep = value
      else if ( s_eqi(name,'wateb') ) then
        wateb = value
      else if ( s_eqi(name,'wateb1') ) then
        wateb1 = value
      else if ( s_eqi(name,'wateb2') ) then
        wateb2 = value
      else if ( s_eqi(name,'xbleft') ) then
        xbleft = value
      else if ( s_eqi(name,'xbltar') ) then
        xbltar = value
      else if ( s_eqi(name(1:6),'xbord(') ) then

        call chrcti(name(7:),ival,ierror,lchar)

        if ( ierror/= 0 ) then
          write ( *, * ) ' '
          write ( *, * ) 'INPUT - Fatal error!'
          write ( *, * ) '  ChrCTI returned nonzero error flag!'
          stop
        end if

        if ( ival<1 .or. ival>2*maxnx-1 ) then
          write ( *, * ) ' '
          write ( *, * ) 'INPUT - Fatal error!'
          write ( *, * ) '  Index of XBORD is out of bounds!'
          stop
        end if

        xbord(ival) = value

      else if ( s_eqi(name,'xbrite') ) then
        xbrite = value
      else if ( s_eqi(name,'xbrtar') ) then
        xbrtar = value
      else if ( s_eqi(name,'xprof') ) then
        xprof = value
      else if ( s_eqi(name,'ybleft') ) then
        ybleft = value
      else if ( s_eqi(name,'ybltar') ) then
        ybltar = value
      else if ( s_eqi(name(1:6),'ybord(') ) then

        call chrcti(name(7:),ival,ierror,lchar)

        if ( ierror/= 0 ) then
          write ( *, * ) ' '
          write ( *, * ) 'INPUT - Fatal error!'
          write ( *, * ) '  ChrCTI returned nonzero error flag!'
          stop
        end if

        if ( ival<1 .or. ival>2*maxny-1 ) then
          write ( *, * ) ' '
          write ( *, * ) 'INPUT - Fatal error!'
          write ( *, * ) '  Index of YBORD is out of bounds!'
          stop
        end if

        ybord(ival) = value

      else if ( s_eqi(name,'ybrite') ) then
        ybrite = value
      else if ( s_eqi(name,'ybrtar') ) then
        ybrtar = value
!
!  Unknown name.
!
      else

        write ( *, * ) ' '
        write ( *, * ) 'INPUT - Unknown variable!'
        write ( *, * ) '  Variable name = "' // trim ( name ) // '".'
        write ( *, * ) '  Assigned value = ',value
        write ( *, * ) ' '

      end if
!
!  IERROR = 2, possible "STOP" or "GO" statement.
!
    else if ( ierror == 2 ) then

      if ( s_eqi(name,'go') ) then
        write ( *, * ) ' '
        write ( *, * ) 'GO command!'
        ierror = 0
        exit
      else if ( s_eqi(name,'stop') ) then
        ierror = 1
        write ( *, * ) 'STOP command!'
        exit
      else
        write ( *, * ) ' '
        write ( *, * ) 'INPUT - Fatal error!'
        write ( *, * ) '  NAMELS error of type ',ierror
        write ( *, * ) '  Card follows:'
        write(*,'(a)')name
        write ( *, * ) ' '
        stop
      end if
!
!  IERROR = 1, blank line.
!
    else if ( ierror == 1 ) then
      write ( *, * ) ' '
!
!  IERROR = 3, or 4, miscellaneous error.
!
    else if ( ierror == 3 .or. ierror==4 ) then
      write ( *, * ) ' '
      write ( *, * ) 'INPUT - Warning!'
      write ( *, * ) '  NAMELS error of type ',ierror
      write ( *, * ) ' '
!
!  IERROR = 6, comment.
!
    else if ( ierror == 6 ) then
      write(*,'(a)')name(1:lchar)
!
!  IERROR = 5, hard end of input.
!
    else if ( ierror == 5 ) then
      write ( *, * ) ' '
      write ( *, * ) 'INPUT - Warning!'
      write ( *, * ) '  "Hard" end of input.'
      write ( *, * ) '  A total of ',nrec,' records were read.'
      write ( *, * ) ' '
      exit
!
!  IERROR = 7, "soft" end of input.
!
    else if ( ierror == 7 ) then
      write ( *, * ) ' '
      write ( *, * ) 'INPUT - Warning!'
      write ( *, * ) '  "Soft" end of input.'
      write ( *, * ) '  A total of ',nrec,' records were read.'
      exit
!
!  Unrecognized error.
!
    else
      write ( *, * ) ' '
      write ( *, * ) 'INPUT - Fatal error!'
      write ( *, * ) '  Unrecognized error from NAMELS.'
      write ( *, * ) '  IERROR = ',ierror
      write ( *, * ) '  Forcing a STOP!'
      stop
    end if

    nrec = nrec+1

  end do

  return
end
subroutine interv(xt,lxt,x,left,mflag)

!*****************************************************************************80
!
!! INTERV computes LEFT, the maximum value of I so that
!
!    1  <=  I <= LXT
!
!  and
!
!    XT(I)  <=  X.
!
!  The routine is designed to be efficient in the common situation
!  that it is called repeatedly, with X taken from an increasing
!  or decreasing sequence.
!
!  This will happen when a piecewise polynomial is to be graphed.
!  The first guess for LEFT is therefore taken to be the value
!  returned at the previous call and stored in the local variable
!  ILO.
!
!  A first check ascertains that ILO.LT.LXT.  This is necessary
!  since the present call may have nothing to do with the previous
!  call.  Then, if XT(ILO)  <=  X < XT(ILO+1), we set LEFT=ILO
!  and are done after just three comparisons.
!
!  Otherwise, we repeatedly double the difference ISTEP = IHI-ILO
!  while also moving ILO and IHI in the direction of X, until
!    XT(ILO)  <=  X < XT(IHI)
!  after which we use bisection to get, in addition, ILO+1 = IHI.
!  LEFT = ILO is then returned.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    Original FORTRAN77 version by Carl DeBoor.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!  XT     Input, real ( kind = 8 ) XT(LXT), a nondecreasing sequence of
!         values.
!
!  LXT    Input, integer ( kind = 4 ) LXT, the dimension of XT.
!
!  X      Input, real ( kind = 8 ) X, the point whose location with
!         respect to the sequence XT is to be determined.
!
!  LEFT,
!  MFLAG  Output, integer ( kind = 4 ) LEFT, integer MFLAG, whose value is
!
!         1     -1      if               X <  XT(1)
!         I      0      if   XT(I)   <=  X < XT(I+1)
!         LXT    1      if  XT(LXT)  <=  X
!
!        In particular, MFLAG = 0 is the 'usual' case.  MFLAG/=0
!        indicates that X lies outside the half open interval
!        XT(1) <= Y<XT(LXT).  The asymmetric treatment of the
!        interval is due to the decision to make all piecewise
!        polynomials continuous from the right.
!
  implicit none

  integer ( kind = 4 ) lxt

  integer ( kind = 4 ) left
  integer ( kind = 4 ) mflag
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) istep
  integer ( kind = 4 ) middle
  real ( kind = 8 ) x
  real ( kind = 8 ) xt(lxt)

  save ilo

  data ilo / 1 /

  ihi = ilo+1

  if ( ihi >= lxt ) then

    if ( x >= xt(lxt))go to 110

    if ( lxt <= 1 ) then
      mflag = -1
      left = 1
      return
    end if

    ilo = lxt-1
    ihi = lxt

  end if
!
  if (x >= xt(ihi))go to 40

  if ( x >= xt(ilo) ) then
    mflag = 0
    left = ilo
    return
  end if
!
!  Now X < XT(ILO).  Decrease ILO to capture X.
!
  istep = 1

   31 continue

  ihi = ilo
  ilo = ihi - istep

  if ( ilo > 1 ) then
    if (x >= xt(ilo))go to 50
    istep = istep*2
    go to 31
  end if

  ilo = 1

  if ( x<xt(1) ) then
    mflag = -1
    left = 1
    return
  end if

  go to 50
!
!  Now X  = > XT(IHI).  Increase IHI to capture X.
!
   40 continue

  istep = 1

   41 continue

  ilo = ihi
  ihi = ilo + istep

  if ( ihi<lxt ) then
    if ( x<xt(ihi))go to 50
    istep = istep*2
    go to 41
  end if

  if (x >= xt(lxt))go to 110

  ihi = lxt
!
!  Now XT(ILO)  <=  X < XT(IHI).  Narrow the interval.
!
   50 continue

  middle = (ilo + ihi)/2

  if ( middle == ilo ) then
    mflag = 0
    left = ilo
    return
  end if
!
!  It is assumed that MIDDLE = ILO in case IHI = ILO+1.
!
  if ( x >= xt(middle) ) then
    ilo = middle
  else
    ihi = middle
  end if

  go to 50
!
!  Set output and return.
!
  110 continue

  mflag = 1

  if ( x == xt(lxt) ) then
    mflag = 0
  end if

  do left = lxt,1,-1
    if ( xt(left)<xt(lxt))return
  end do

  return
end
function isamax ( n, x, incx )

!*****************************************************************************80
!
!! ISAMAX finds the index of the vector element of maximum absolute value.
!
!  Modified:
!
!    08 April 1999
!
!  Reference:
!
!    Lawson, Hanson, Kincaid, Krogh,
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
!    Output, integer ( kind = 4 ) ISAMAX, the index of the element of SX of maximum
!    absolute value.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) isamax
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) n
  real ( kind = 8 ) samax
  real ( kind = 8 ) x(*)

  if ( n  <= 0 ) then

    isamax = 0

  else if ( n == 1 ) then

    isamax = 1

  else if ( incx == 1 ) then

    isamax = 1
    samax = abs ( x(1) )

    do i = 2, n

      if ( abs ( x(i) ) > samax ) then
        isamax = i
        samax = abs ( x(i) )
      end if

    end do

  else

    if ( incx >= 0 ) then
      ix = 1
    else
      ix = ( - n + 1 ) * incx + 1
    end if

    isamax = 1
    samax = abs ( x(ix) )

    ix = ix + incx

    do i = 2, n
      if ( abs ( x(ix) ) > samax ) then
        isamax = i
        samax = abs ( x(ix) )
      end if
      ix = ix + incx
    end do

  end if

  return
end
subroutine lbase ( ival, npol, pval, xpol, xval )

!*****************************************************************************80
!
!! LBASE evalualates the IVAL-th Lagrange polynomial.
!
!  Discussion:
!
!    The polynomial is defined by a set of NPOL points XPOL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IVAL, the polynomial to evaluate.
!    IVAL should be between 1 and NPOL.
!
!    Input, integer ( kind = 4 ) NPOL, the number of points that define
!    the Lagrange polynomials.
!
!    Output, real ( kind = 8 ) PVAL, the value of the IVAL-th
!    Lagrange polynomial at the point XVAL.
!
!    Input, real ( kind = 8 ) XPOL(NPOL), the abscissas of the
!    Lagrange polynomials.  The entries in XPOL should be
!    distinct.
!
!    Input, real ( kind = 8 ) XVAL, the point at which the
!    IVAL-th Lagrange polynomial is to be evaluated.
!
  implicit none

  integer ( kind = 4 ) npol

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ival
  real ( kind = 8 ) pval
  real ( kind = 8 ) xpol(npol)
  real ( kind = 8 ) xval

  pval = 1.0D+00
  do i = 1,npol
    if ( i/= ival ) then
      pval = pval*(xval-xpol(i))/(xpol(ival)-xpol(i))
    end if
  end do

  return
end
subroutine lmemry ( action, name, lval )

!*****************************************************************************80
!
!! LMEMRY allows the user to define the name of a logical variable,
!  set it, increment it, or get the value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character*(*) ACTION, desired action.
!    'Init', reset all values to FALSE, wipe out all names.
!    'Name', add a variable of the given name.
!    'Not',  replace value of NAME by NOT.NAME.
!    'Set',  set variable NAME to LVAL.
!    'Get',  return value of NAME in LVAL.
!
!    Input, character*(*) NAME, the name of the variable.
!
!    Input/output, logical LVAL.
!    For the 'Set' command, LVAL must contain the set value.
!    For the 'Get' command, LVAL will contain the value of the
!    named variable on output.
!
  implicit none

  integer ( kind = 4 ) maxnam
  parameter (maxnam = 100)

  character ( len = * ) action
  integer ( kind = 4 ) i
  integer ( kind = 4 ) lenc
  logical lval
  logical lvals(maxnam)
  logical s_eqi
  character ( len = * ) name
  character ( len = 20 ) names(maxnam)
  integer ( kind = 4 ) numnam

  save lvals
  save names
  save numnam
!
!  Initialize everything.
!
  if ( s_eqi(action,'init') ) then

    numnam = 0

    do i = 1,maxnam
      lvals(i) = .false.
      names(i) = ' '
    end do
!
!  Name something.
!
  else if ( s_eqi(action,'name') ) then

    do i = 1,numnam

      if ( s_eqi(name,names(i)) ) then
        write ( *, * ) ' '
        write ( *, * ) 'LMemry - Warning!'
        write(*,'(''  There is ALREADY a variable '',a)') trim ( name )
        return
      end if

    end do

    if ( numnam<maxnam ) then
      numnam = numnam+1
      names(numnam) = name
    else
      write ( *, * ) ' '
      write ( *, * ) 'LMemry - Fatal error!'
      write ( *, * ) '  No more name space.'
      stop
    end if
!
!  Switch something.
!
  else if ( s_eqi(action,'not') ) then

    do i = 1,numnam

      if ( s_eqi(name,names(i)) ) then
        lvals(i) = .not.lvals(i)
        return
      end if

    end do

    write ( *, * ) ' '
    write ( *, * ) 'LMemry - Fatal error!'
    write ( *, * ) '  Attempt to NOT unknown variable.'
    write(*,'(''  Variable name is '',a)')name
    stop
!
!  Set something.
!
  else if ( s_eqi(action,'set') ) then

    do i = 1,numnam

      if ( s_eqi(name,names(i)) ) then
        lvals(i) = lval
        return
      end if

    end do

    write ( *, * ) ' '
    write ( *, * ) 'LMemry - Fatal error!'
    write ( *, * ) '  Attempt to set unknown variable.'
    write(*,'(''  Variable name is '',a)')name
    stop
!
!  Get something.
!
  else if ( s_eqi(action,'get') ) then

    do i = 1,numnam

      if ( s_eqi(name,names(i)) ) then
        lval = lvals(i)
        return
      end if

    end do

    write ( *, * ) ' '
    write ( *, * ) 'LMemry - Fatal error!'
    write ( *, * ) '  Attempt to get value of unknown variable.'
    write(*,'(''  Variable name is '',a)')name
    stop

  else
    write ( *, * ) ' '
    write ( *, * ) 'LMemry - Fatal error!'
    write ( *, * ) '  Unrecognized action.'
    write(*,'(1x,a)')action
    stop
  end if

  return
end
function lsame(ca,cb)

!*****************************************************************************80
!
!! LSAME returns .TRUE. if CA is the same letter as CB regardless of
!  case.
!
!  Parameters:
!
!  CA      (input) CHARACTER*1
!
!  CB      (input) CHARACTER*1
!          CA and CB specify the single characters to be compared.
!
!    Output, logical LSAME, is TRUE if the characters are equivalent.
!
  implicit none

  character          ca, cb
  logical lsame
  integer ( kind = 4 )            inta, intb, zcode
!     ..
!     .. Executable Statements ..
!
!     Test if the characters are equal
!
  lsame = ca == cb
  if (  lsame )return
!
!     Now test for equivalence if both characters are alphabetic.
!
  zcode = ichar( 'z' )
!
!     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
!     machines, on which ICHAR returns a value with bit 8 set.
!     ICHAR('A') on Prime machines returns 193 which is the same as
!     ICHAR('A') on an EBCDIC machine.
!
  inta = ichar( ca )
  intb = ichar( cb )
!
  if (  zcode == 90 .or. zcode==122 ) then
!
!        ASCII is assumed - ZCODE is the ASCII code of either lower or
!        upper case 'Z'.
!
     if (  inta >= 97 .and. inta<=122 ) inta = inta - 32
     if (  intb >= 97 .and. intb<=122 ) intb = intb - 32
!
  else if (  zcode == 233 .or. zcode==169 ) then
!
!        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
!        upper case 'Z'.
!
     if (  inta >= 129 .and. inta<=137 .or. &
           inta >= 145 .and. inta<=153 .or. &
           inta >= 162 .and. inta<=169 ) inta = inta + 64
     if (  intb >= 129 .and. intb<=137 .or. &
           intb >= 145 .and. intb<=153 .or. &
           intb >= 162 .and. intb<=169 ) intb = intb + 64

  else if (  zcode == 218 .or. zcode==250 ) then
!
!        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
!        plus 128 of either lower or upper case 'Z'.
!
     if (  inta >= 225 .and. inta<=250 ) inta = inta - 32
     if ( intb >= 225 .and. intb <=250 ) intb = intb - 32
  end if

  lsame = inta == intb

  return
end
subroutine lspoly ( coef, ncoef, x, y, val )

!*****************************************************************************80
!
!! LSPOLY evaluates a polynomial in (1, x, y, x*x, xy, y*y) at (x,y)..
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) COEF(NCOEF), the polynomial coefficients.
!
!    Input, integer ( kind = 4 ) NCOEF, the number of coefficients.
!
!    Input, real ( kind = 8 ) X, Y, the values of X and Y.
!
!    Output, real ( kind = 8 ) VAL, the value of the polynomial.
!
  implicit none

  integer ( kind = 4 ) ncoef

  real ( kind = 8 ) coef(ncoef)
  real ( kind = 8 ) val
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  val = coef(1) &
      + coef(2) * x &
      + coef(3) * y &
      + coef(4) * x * x &
      + coef(5) * x * y &
      + coef(6) * y * y

  return
end
subroutine march ( a,area,base,dir,disjac,dpara3,dparsn,dparfd,dparfdc,dpdyn, &
  dudyn,dvdyn,dydpn,eqn,etan,etaq,g,gdif,gdifc,gold,gradf, &
  gtar,ibc,idfd,ids,ifds,igrid,igunit,ijac,indx,iopt,ipivot,iplot,ipred, &
  ishapb,ishapf,ismooth,isotri,istep1,istep2,itunit,itype,iwrite,jjac,jstep1, &
  jstep2,maxnew,maxstp,ndim,nelem,neqn,nlband,node,np,npar,nparb,nparf,npe, &
  nprof, &
  nrow,nstep3,numel,nx,ny,para,para1,para2,para3,parjac,parnew,parold,phi, &
  res,sens,splbmp,splflo,stpmax,syseqn,taubmp,tauflo,tolnew,wateb,wateb1, &
  wateb2,watep,wateu,watev,wquad,xbl,xbord,xbr,xc,xprof,xquad,xsin,xsiq,ybl, &
  ybord,ybr,yc,yquad)

!*****************************************************************************80
!
!! MARCH carries out a one, two or three dimensional "march".
!
!  In a one dimensional march:
!
!    Step:                      ParNew
!
!    IStep1                     Para1
!    IStep2                     Para2
!
!    Para(IStep) = ( (IStep2-IStep)*Para1 + (IStep-IStep1)*Para2)
!                  / (IStep2-IStep1)
!
!  In a two dimensional march:
!
!    Step:                      ParNew
!
!    IStep1, JStep1             Para1
!    IStep2, JStep2             Para3
!
!    Para(Istep,Jstep) = Para1 + (IStep-IStep1)*Dir1
!    + (J (JStep-JStep1)*Dir2
!
!  The vector Para2 is used to determine the orthogonal directions Dir1
!  and Dir2.  Dir1 is the direction from Para1 to Para2.  Dir2 is the
!  vector which is normal to Dir1 and passes through Para3.
!
!  In a three dimensional march:
!
!    Step:                      ParNew    WateB
!
!    IStep1, JStep1, 1          Para1     WateB1
!    IStep2, JStep2, 1          Para3     WateB1
!
!    IStep1, JStep1, NStep3     Para1     WateB2
!    IStep2, JStep2, NStep3     Para3     WateB2
!
!    Para(IStep,JStep.KStep)  = Para1 + (IStep-IStep1)*Dir1
!                                     + (JStep-JStep1)*Dir2
!    WateB(IStep,JStep,KStep) = ((NStep3-KStep)*WateB1+(KStep-1)*WateB2)
!                               / (NStep3-1)
!
!  The vector Para2 is used to determine the orthogonal directions Dir1
!  and Dir2.  Dir1 is the direction from Para1 to Para2.  Dir2 is the
!  vector which is normal to Dir1 and passes through Para3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ndim
  integer ( kind = 4 ) nelem
  integer ( kind = 4 ) neqn
  integer ( kind = 4 ) np
  integer ( kind = 4 ) npar
  integer ( kind = 4 ) nparb
  integer ( kind = 4 ) nparf
  integer ( kind = 4 ) npe
  integer ( kind = 4 ) nrow
  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny

  real ( kind = 8 ) a(nrow,neqn)
  real ( kind = 8 ) area(nelem,3)
  real ( kind = 8 ) base(npar,ndim)
  real ( kind = 8 ) cost
  real ( kind = 8 ) costb
  real ( kind = 8 ) costp
  real ( kind = 8 ) costu
  real ( kind = 8 ) costv
  real ( kind = 8 ) dir(npar,ndim)
  real ( kind = 8 ) disjac
  real ( kind = 8 ) dist(3)
  real ( kind = 8 ) dpara3(npar)
  real ( kind = 8 ) dparsn(npar)
  real ( kind = 8 ) dparfd(npar)
  real ( kind = 8 ) dparfdc(npar)
  real ( kind = 8 ) dpdyn(np)
  real ( kind = 8 ) dudyn(np)
  real ( kind = 8 ) dvdyn(np)
  real ( kind = 8 ) dydpn(np,nparb)
  character ( len = 2 ) eqn(neqn)
  real ( kind = 8 ) etan(6)
  real ( kind = 8 ) etaq(3)
  real ( kind = 8 ) flarea
  real ( kind = 8 ) g(neqn)
  real ( kind = 8 ) g3(neqn)
  real ( kind = 8 ) gdif(neqn,npar)
  real ( kind = 8 ) gdifc(neqn,npar)
  real ( kind = 8 ) gold(neqn)
  real ( kind = 8 ) gpoptfd(3)
  real ( kind = 8 ) gpoptsn(3)
  real ( kind = 8 ) gpopt3(3)
  real ( kind = 8 ) gradf(neqn,npar)
  real ( kind = 8 ) gtar(neqn)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibc
  integer ( kind = 4 ) idfd
  integer ( kind = 4 ) ids
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ifds
  integer ( kind = 4 ) igrid
  integer ( kind = 4 ) igunit
  integer ( kind = 4 ) ijac
  integer ( kind = 4 ) indx(np,3)
  integer ( kind = 4 ) iopt(npar)
  integer ( kind = 4 ) ipar
  integer ( kind = 4 ) ipivot(neqn)
  integer ( kind = 4 ) iplot
  integer ( kind = 4 ) ipred
  integer ( kind = 4 ) ishapb
  integer ( kind = 4 ) ishapf
  integer ( kind = 4 ) ismooth
  integer ( kind = 4 ) isotri(nelem)
  integer ( kind = 4 ) istep
  integer ( kind = 4 ) istep1
  integer ( kind = 4 ) istep2
  integer ( kind = 4 ) itunit
  integer ( kind = 4 ) itype
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) iwrite
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jjac
  integer ( kind = 4 ) jstep
  integer ( kind = 4 ) jstep1
  integer ( kind = 4 ) jstep2
  integer ( kind = 4 ) kstep
  integer ( kind = 4 ) maxnew
  integer ( kind = 4 ) maxstp
  integer ( kind = 4 ) nlband
  integer ( kind = 4 ) node(nelem,npe)
  integer ( kind = 4 ) nprof(2*ny-1)
  integer ( kind = 4 ) nstep1
  integer ( kind = 4 ) nstep2
  integer ( kind = 4 ) nstep3
  integer ( kind = 4 ) numel(np)
  integer ( kind = 4 ) numstp
  real ( kind = 8 ) para(npar)
  real ( kind = 8 ) para1(npar)
  real ( kind = 8 ) para2(npar)
  real ( kind = 8 ) para3(npar)
  real ( kind = 8 ) parjac(npar)
  real ( kind = 8 ) parnew(npar)
  real ( kind = 8 ) parold(npar)
  real ( kind = 8 ) phi(nelem,3,6,10)
  real ( kind = 8 ) res(neqn)
  real ( kind = 8 ) sens(neqn,npar)
  real ( kind = 8 ) splbmp(4,nparb+2,0:nparb)
  real ( kind = 8 ) splflo(4,nparf+2,0:nparf)
  real ( kind = 8 ) stpmax
  character ( len = 20 ) syseqn
  real ( kind = 8 ) taubmp(nparb+2)
  real ( kind = 8 ) tauflo(nparf+2)
  character ( len = 80 ) title
  real ( kind = 8 ) tolnew
  real ( kind = 8 ) wateb
  real ( kind = 8 ) wateb1
  real ( kind = 8 ) wateb2
  real ( kind = 8 ) watep
  real ( kind = 8 ) wateu
  real ( kind = 8 ) watev
  real ( kind = 8 ) wquad(3)
  real ( kind = 8 ) xbl
  real ( kind = 8 ) xbord(nx)
  real ( kind = 8 ) xbr
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) xprof
  real ( kind = 8 ) xquad(nelem,3)
  real ( kind = 8 ) xsin(6)
  real ( kind = 8 ) xsiq(3)
  real ( kind = 8 ) ybl
  real ( kind = 8 ) ybord(ny)
  real ( kind = 8 ) ybr
  real ( kind = 8 ) yc(np)
  real ( kind = 8 ) yquad(nelem,3)

  if ( maxstp  <= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'MARCH - Warning!'
    write ( *, * ) '  No marching steps were requested!'
    return
  end if

  if ( istep1 == istep2 ) then
    write ( *, * ) ' '
    write ( *, * ) 'MARCH - Fatal error!'
    write ( *, * ) '  ISTEP1 = ISTEP2.'
    return
  end if

  write ( *, * ) ' '
  write ( *, * ) 'MARCH - Begin march!'

  numstp = 0
  gpoptfd(1:3) = 0.0D+00
  gpoptsn(1:3) = 0.0D+00
  gpopt3(1:3) = 0.0D+00
!
!  Set up certain marching indices and variables so that they
!  are correct no matter whether this is a 1, 2 or 3 dimensional march.
!
  if ( ndim == 1 ) then
    jstep1 = 1
    jstep2 = 1
    nstep1 = maxstp
    nstep2 = 1
    nstep3 = 1
    wateb1 = wateb
    wateb2 = wateb
  else if ( ndim == 2 ) then
    nstep1 = maxstp
    nstep2 = maxstp
    nstep3 = 1
  else if ( ndim == 3 ) then
    nstep1 = maxstp
    nstep2 = maxstp
    wateb1 = wateb
    wateb2 = wateb
  end if

  if ( itunit /= 0 ) then
    write(itunit,*)nstep1,nstep2,nstep3
  end if
!
!  Compute the direction vectors for the marching plane.
!  In 3 dimensions, we only compute two direction vectors, since the
!  bump weight is assumed to be the only variable in the third
!  direction.
!
  do ipar = 1, npar

    dir(ipar,1) = para2(ipar) - para1(ipar)

    if ( ndim >= 2 ) then
      dir(ipar,2) = para3(ipar)-para1(ipar)
    end if

    if ( ndim >= 3 ) then
      dir(ipar,3) = 0.0D+00
    end if

  end do
!
!  Set up the vectors that span the marching space.
!
  do i = 1, npar

    base(i,1) = para2(i)-para1(i)

    if ( ndim >= 2 ) then
      base(i,2) = para3(i)-para1(i)
    end if

    if ( ndim >= 3 ) then
      base(i,3) = 0.0D+00
    end if

  end do
!
!  Use orthonormalization to get a marching basis.
!
  call probas(base,npar,ndim)
!
!  Compute DIST, the length of a single marching step in each direction.
!
  if ( ndim == 1 ) then

    dist(1) = 0.0D+00
    do j = 1, npar
      dist(1) = dist(1) + ( para2(j) - para1(j) )**2
    end do
    dist(1) = sqrt ( dist(1) )
    dist(1) = dist(1) / real ( istep2 - istep1, kind = 8 )
    write ( *, * ) ' '
    write ( *, * ) 'MARCH: Distance = ',dist(1)

  else if ( ndim == 2 .or. ndim == 3 ) then

    do j = 1,npar
      dir(j,1) = para3(j)-para1(j)
    end do

    write ( *, * ) ' '
    write ( *, * ) 'MARCH: Para3-Para1:'
    write ( *, * ) ' '
    call pr_parameter(nparb,nparf,dir)

    call projec(base,npar,2,dir,dist)

    dist(1) = dist(1) / real ( istep2 - istep1, kind = 8 )
    dist(2) = dist(2) / real ( jstep2 - jstep1, kind = 8 )
    write ( *, * ) ' '
    write ( *, * ) 'MARCH: Distances = ',dist(1),dist(2)
  end if
!
!  Now do the marching.
!
  do kstep = 1,nstep3

    if ( nstep3>1 ) then
      wateb = ( real ( nstep3 - kstep,     kind = 8 ) * wateb1 &
              + real (          kstep - 1, kind = 8 ) * wateb2 ) &
              / real ( nstep3         - 1, kind = 8 )

      write ( *, * ) ' '
      write ( *, * ) 'MARCH - Value of WATEB is now ',wateb

    end if
!
!  Because of special nature of the third dimension, we
!  wait til now to start G and PARNEW at zero.
!
    g(1:neqn) = 0.0D+00
    parnew(1:npar) = 0.0D+00
!
!   parnew(nparf+nparb+1) = 1.0D+00
!
!  Now compute MAXSTP lines of solutions, each containing
!  MAXSTP points.
!
    do jstep = 1,nstep2

      do istep = 1,nstep1

        ival = 1
        call imemry('inc','Points_Opt',ival)
!
!  Copy the old solution and parameters.
!
        do i = 1,neqn
          gold(i) = g(i)
        end do

        do i = 1,npar
          parold(i) = parnew(i)
        end do
!
!  Set the current parameters, PARNEW.
!
        if ( ndim == 1 ) then

          do j = 1,npar
            parnew(j) = ( real ( istep2 - istep,          kind = 8 ) * para1(j)   &
                        + real (          istep - istep1, kind = 8 ) * para2(j) ) &
                        / real ( istep2         - istep1, kind = 8 )
          end do

        else

          do j = 1,npar
            parnew(j) = para1(j) &
              + real ( istep - istep1, kind = 8 ) * dist(1) * base(j,1) &
              + real ( jstep - jstep1, kind = 8 ) * dist(2) * base(j,2)
          end do

        end if
!
!  For the first step along a line of solutions, retrieve the
!  solution computed on the previous line.
!
        if ( istep == 1 .and. jstep /= 1 ) then
          g(1:neqn) = g3(1:neqn)
        end if
!
!  Compute the flow variables G that correspond to PARNEW.
!
        call solcon(a,area,disjac,dpara3,dpdyn,dudyn,dvdyn,dydpn,eqn, &
          etan,etaq,flarea,g,gdif,gdifc,gold,gradf,gtar,ibc,idfd,ids, &
          ierror,ifds,igrid,ijac,indx,iopt,ipivot,ipred,ishapb,ishapf, &
          ismooth,isotri,itype,iwrite,jjac,maxnew,nelem,neqn,nlband,node,np, &
          npar,nparb,nparf,npe,nprof,nrow,numel,nx,ny,para,parjac,parnew,parold, &
          phi,res,sens,splbmp,splflo,stpmax,syseqn,taubmp,tauflo,tolnew,wateb, &
          watep,wateu,watev,wquad,xbl,xbord,xbr,xc,xquad,xsin,xsiq,ybl,ybord, &
          ybr,yc,yquad)

        if ( ierror/= 0 ) then
          write ( *, * ) ' '
          write ( *, * ) 'MARCH - Fatal error!'
          write ( *, * ) '  This solution point is not valid!'
          write ( *, * ) '  SolCon returned IERROR = ',ierror
          write ( *, * ) '  Step ',' I,J,K = ',istep,jstep,kstep
          call pr_parameter(nparb,nparf,para)
          return
        end if
!
!  For the first step along a line of solutions, save the
!  solution for use on the next line.
!
        if ( istep == 1 ) then
          g3(1:neqn) = g(1:neqn)
        end if
!
!  Compute the cost COST associated with the solution G,
!  determined by the parameters PARNEW.
!
        call get_cost(cost,costb,costp,costu,costv,g,gtar,indx,ishapb,neqn,np, &
          nparb,nprof,ny,splbmp,taubmp,wateb,watep,wateu,watev,xbl,xbr,ybl, &
          ybr,yc)

        if ( iwrite == 1 ) then
          title = 'Cost'
          call pr_cost1 ( cost, title )
        else if ( iwrite >= 2 ) then
          title = 'Cost'
          call pr_cost2(cost,costb,costp,costu,costv,title,wateb,watep,wateu, &
            watev)
        end if
!
!  Print stuff.
!
        numstp = numstp + 1

        if ( ifds/= 0 ) then
          call cost_gradient(dparfd,g,gtar,indx,ishapb,neqn,np,npar, &
            nparb,nparf,nprof,ny, &
            gdif,splbmp,taubmp,wateb,watep,wateu,watev,xbl,xbr,ybl,ybr,yc)
        end if

        if ( ifds/= 0 ) then
          call cost_gradient(dparfdc,g,gtar,indx,ishapb,neqn,np,npar, &
            nparb,nparf,nprof,ny, &
            gdifc,splbmp,taubmp,wateb,watep,wateu,watev,xbl,xbr,ybl,ybr,yc)
        end if

        if ( ids/= 0 ) then
          call cost_gradient(dparsn,g,gtar,indx,ishapb,neqn,np,npar, &
            nparb,nparf,nprof,ny, &
            sens,splbmp,taubmp,wateb,watep,wateu,watev,xbl,xbr,ybl,ybr,yc)
        end if

        title = 'Marching solution point'

        call pr_solution ( dpara3,dparfd,dparfdc,dparsn,g,gdif,gdifc,gtar,idfd,ids, &
          ifds,indx,iwrite,neqn,np,npar,nparb,nparf,nprof,numstp,ny,parnew, &
          sens,title,yc)

        if ( itunit/= 0 ) then
          call projec(base,npar,ndim,dpara3,gpopt3)
          call projec(base,npar,ndim,dparfd,gpoptfd)
          call projec(base,npar,ndim,dparsn,gpoptsn)
          write(itunit,'(1x,6g14.6)')cost,(parnew(i),i = 1,npar)
          write(itunit,'(1x,5g14.6)')(gpoptsn(i),i = 1,ndim)
          write(itunit,'(1x,5g14.6)')(gpoptfd(i),i = 1,ndim)
          write(itunit,'(1x,5g14.6)')(gpopt3(i),i = 1,ndim)
        end if

        para(1:npar) = parnew(1:npar)
!
!  Shall we save DISPLAY graphics information to a file?
!
        if ( iplot < 0 ) then
          if ( (istep == istep1 .and. jstep==jstep1.and.kstep==1) .or.  &
               (istep == istep2 .and. jstep==jstep2.and.kstep==nstep3) ) then

            call plot_file_write ( eqn, g, gdif, igunit, indx, isotri, iwrite, nelem, &
              neqn, node, np, npar, npe, nprof, nx, ny, para, sens, xc, xprof, yc )
          end if
        end if

      end do
    end do
  end do
!
!  The last record in a marching file that includes finite differences
!  should be the lengths of the sides.
!
  if ( itunit /= 0 ) then
    write ( itunit,'(1x,5g14.6)') dist(1:ndim)
  end if

  return
end
subroutine march_file_open ( march_file, itunit )

!*****************************************************************************80
!
!! MARCH_FILE_OPEN opens the marching file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) MARCH_FILE, the file into which the marching
!    information will be stored.  If MARCH_FILE = 'none', or ' ', then
!    no marching file is opened.
!
!    Output, integer ( kind = 4 ) ITUNIT, the FORTRAN unit used for writing data to the
!    marching file.
!
  implicit none

  character ( len = * ) march_file
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) itunit
  logical s_eqi
!
!  Return immediately if no marching file was requested.
!
  if ( s_eqi ( march_file, 'none' ) .or. march_file == ' ' ) then
    itunit = 0
    return
  end if

  if ( itunit == 0 ) then

    write ( *, * ) ' '
    write ( *, * ) 'MARCH_FILE_OPEN - Note:'
    write ( *, * ) '  Opening the marching file "'// trim ( march_file ) // '".'
!
!  Get a free FORTRAN unit number.
!
    call get_unit ( itunit )
!
!  Open the file.
!
    open ( unit = itunit, file = march_file, status = 'replace', &
      form = 'formatted', access = 'sequential', iostat = ios )

    if ( ios /= 0 ) then
      write ( *, * ) ' '
      write ( *, * ) 'MARCH_FILE_OPEN - Warning!'
      write ( *, * ) '  The marching file could not be opened.'
      itunit = 0
    end if
!
!  Marching file is already open.
!
  else
    write ( *, * ) ' '
    write ( *, * ) 'MARCH_FILE_OPEN - Note:'
    write ( *, * ) '  The marching file is already opened.'
    write ( *, * ) '  New information will be appended.'
  end if

  return
end
subroutine namels ( name, ierror, rhs, value )

!*****************************************************************************80
!
!! NAMELS interprets input data in namelist form.
!
!  Discussion:
!
!    The routine reads a line of user input which is similar in form
!    to NAMELIST input, and returns the name of the variable
!    and its value.
!
!    NAMELS is a simple program, and can only handle simple input.
!    In particular, it cannot handle:
!
!    * multiple assignments on one line,
!    * a single assignment extended over multiple lines,
!    * assignments to character or complex variables,
!    * assignments to arrays.
!
!    Typical input would be of the form:
!
!      name = value
!
!    including, for instance:
!
!      a = 1.0D+00
!      n = -17
!      scale = +5.3E-2
!
!    Spaces are ignored, and case is not important.  Integer values
!    will be returned as real, but this is never a
!    problem as long as the integers are "small".
!
!    If a line begins with the character "#", it is assumed to be
!    a comment, and is ignored.  IERROR is returned as 6.
!
!    If a line begins with the characters "end-of-input", it is
!    assumed to be an "end-of-input" marker, and IERROR is returned
!    as 7.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) NAME.
!    NAME contains the left hand side of the assignment
!    statement.
!    Normally, this will be the name of a variable.
!    If the input line was blank, then NAME will equal ' '.
!    If an error occurred while trying to process the
!    input line, NAME will contain the text of the line..
!    If the line began with "#", then NAME will contain the
!    text of the line.
!    If the line equals "end-of-input", then NAME will contain
!    the text of the line.
!
!    Output, integer ( kind = 4 ) IERROR.
!    0, no errors were detected.
!    1, the line was blank.
!    2, the line did not contain an " = " sign.
!    3, the line did not contain a variable name to the
!       left of the " = " sign.
!    4, the right hand side of the assignment did not make sense.
!    5, end of input.
!    6, the line began with "#", signifying a comment.
!       The text of the line is returned in NAME.
!    7, the line began with "end-of-input".
!
!    Output, real ( kind = 8 ) VALUE.
!    VALUE contains the right hand side of the assignment
!    statement.
!    Normally, this will be a real value.
!    But if the input line was blank, or if an error occurred
!    while trying to process the input line, or if input
!    terminated, then VALUE will simply be set to 0.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) lchar
  logical s_eqi
  character ( len = 80 ) line
  character ( len = * ) name
  integer ( kind = 4 ) nchar
  character ( len = 80 ) rhs
  real ( kind = 8 ) value
!
!  Set default values
!
  ierror = 0
  name = ' '
  value = 0.0D+00
!
!  Read a line
!
  read ( *, '(a)', end = 20 ) line
!
!  Empty lines are OK
!
  if ( len_trim ( line ) <= 0 ) then
    ierror = 1
    return
  end if
!
!  Check for a comment.
!
  if ( line(1:1) == '#' ) then
    ierror = 6
    name = line
    return
  end if
!
!  Check for "end-of-line".
!
  if ( s_eqi(line,'end-of-input') ) then
    ierror = 7
    name = line
    return
  end if
!
!  Does the line contain an = sign?
!
  if ( index ( line, '=' ) <= 0 ) then
    ierror = 2
    value = 0
    name = line
    return
  end if
!
!  Find the name of the variable to be assigned.
!
  call s_before_ss_copy ( line, '=', name )

  call s_blank_delete ( name )

  if ( len_trim ( name ) <= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'NAMELS - Warning!'
    write ( *, * ) '  The following input line was ignored, because'
    write ( *, * ) '  there was no variable name on the left hand'
    write ( *, * ) '  side of the assignment statement:'
    write(*,'(a)')line
    write ( *, * ) ' '
    ierror = 3
    return
  end if
!
!  Read the value, as a real number.
!
  nchar = index ( line, '=' )

  rhs = adjustl ( line(nchar+1:) )

  if ( &
    s_eqi ( name, 'march_file' ) .or. &
    s_eqi ( name, 'plot_file' ) .or. &
    s_eqi ( name, 'syseqn' ) ) then
    return
  end if

  call chrctd(line(nchar+1:),value,ierror,lchar)

  if ( ierror/= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'NAMELS - Warning!'
    write ( *, * ) '  The following input line was ignored, because'
    write ( *, * ) '  the right hand side of the assignment statement'
    write ( *, * ) '  did not seem to make sense:'
    write(*,'(a)')line
    write ( *, * ) ' '
    ierror = 4
  end if
  return
!
!  On end of input, return.
!
20    continue
  write ( *, * ) ' '
  write ( *, * ) 'NAMELS - Reached end of input.'
  write ( *, * ) ' '
  ierror = 5

  return
end
subroutine newton ( a,area,disjac,eqn,g,ierror,ijac,indx,ipivot,ishapf, &
  iwrite,jjac,maxnew,nelem,neqn,nlband,node,np,npar,nparb,nparf,npe,nrow,para, &
  parjac,phi,res,splflo,syseqn,tauflo,tolnew,yc)

!*****************************************************************************80
!
!! NEWTON solves the nonlinear system of equations.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) A(NROW,NEQN), used for the value of the jacobian
!    matrix, or for the LU factors of that same matrix.
!    In some cases, NEWTON may use the LU factors stored in
!    A from a previous iteration.  Otherwise, NEWTON will
!    compute the jacobian and factor it.
!
!    Input/output, real ( kind = 8 ) G(NEQN).
!    On input, G contains an initial estimate of the solution.
!    On output, G contains an improved estimate of the
!    solution.
!
!    Workspace, real G2(NEQN).
!    G2 contains the Newton iterate with the lowest residual norm.
!    If we need to restart the Newton iteration, we restart it
!    at G2.
!
!  IERROR Output, integer ( kind = 4 ) IERROR, error flag.
!         0, no error.
!         nonzero, an error occurred.  The calculation was halted
!         without convergence.
!
  implicit none

  integer ( kind = 4 ) nelem
  integer ( kind = 4 ) neqn
  integer ( kind = 4 ) np
  integer ( kind = 4 ) npar
  integer ( kind = 4 ) nparf
  integer ( kind = 4 ) npe
  integer ( kind = 4 ) nrow

  real ( kind = 8 ) a(nrow,neqn)
  real ( kind = 8 ) area(nelem,3)
  real ( kind = 8 ) disjac
  real ( kind = 8 ) dist
  real ( kind = 8 ) dmax
  character ( len = 2 ) eqn(neqn)
  logical falter
  real ( kind = 8 ) g(neqn)
  real ( kind = 8 ) g2(neqn)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) isamax
  integer ( kind = 4 ) idmax
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ijac
  integer ( kind = 4 ) indx(np,3)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipivot(neqn)
  integer ( kind = 4 ) irmax
  integer ( kind = 4 ) ishapf
  integer ( kind = 4 ) iter
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) iwrite
  integer ( kind = 4 ) ixmax
  integer ( kind = 4 ) jjac
  logical leval
  logical lgdif
  logical lmat
  logical ltarg
  logical lzero
  integer ( kind = 4 ) maxnew
  integer ( kind = 4 ) nlband
  integer ( kind = 4 ) node(nelem,npe)
  integer ( kind = 4 ) nparb
  integer ( kind = 4 ) nsys
  real ( kind = 8 ) para(npar)
  real ( kind = 8 ) parjac(npar)
  real ( kind = 8 ) phi(nelem,3,6,10)
  real ( kind = 8 ) res(neqn)
  real ( kind = 8 ) rmax
  real ( kind = 8 ) rmax0
  real ( kind = 8 ) rmax2
  real ( kind = 8 ) splflo(4,nparf+2,0:nparf)
  character ( len = 20 ) syseqn
  real ( kind = 8 ) tauflo(nparf+2)
  real ( kind = 8 ) tolnew
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmax0
  real ( kind = 8 ) xmax2
  real ( kind = 8 ) yc(np)
!
  ival = 1
  call imemry('inc','Newton_calls',ival)
!
!  Initialize the flag that notes whether we have evaluated the
!  jacobian during this process.
!
  leval = .false.
!
!  Find out if this is a finite difference calculation.
!
  lgdif = .false.
  call lmemry('get','fd_grad',lgdif)
!
!  Find out if this is a target calculation.
!
  ltarg = .false.
  call lmemry('get','target',ltarg)
!
!  Find out if we have a valid, factored, jacobian available.
!
  lmat = .false.
  call lmemry('get','have_fp',lmat)

  nsys = 1
  falter = .false.
!
!  If there is no inflow, then the answer is easy.  G = 0.
!  For some reason, the code has trouble figuring this out,
!  so make it easy for it!
!
  lzero = .true.
  do i = 1,nparf
    if ( para(i)/= 0.0D+00 )lzero=.false.
  end do

  if ( lzero ) then
    g(1:neqn) = 0.0D+00
    return
  end if
!
!  Based on the user input quantities JJAC and KJAC, and
!  LGDIF, which tells us whether we are calculating a finite difference,
!  and LMAT, which tells us whether we have a valid jacobian,
!  determine the frequency with which we expect to update
!  the jacobian during this process.
!
!  We may revise this value if the process fails.
!
  if ( jjac == 0 ) then
    lmat = .false.
  else if ( jjac == 1 ) then
    if ( .not.lgdif ) then
      lmat = .false.
    end if
  else if ( jjac == 2 ) then
  else if ( jjac == 3 ) then
  end if
!
!  Start G2 at G.
!  Update G2, if necessary, so that it always contains the point
!  with the lowest residual.
!
  g2(1:neqn) = g(1:neqn)
!
!  If the first Newton iteration failed, you may want to try again
!  by coming back here.
!
10    continue

  ierror = 0
  iter = 0
!
!  Compute the norm of the initial X value.
!
  ixmax = isamax(neqn,g,1)
  xmax = abs(g(ixmax))
  xmax0 = xmax
  xmax2 = xmax

  if ( iwrite >= 3 ) then
    write ( *, * ) ' '
    write ( *, * ) 'NEWTON - Step ', iter, ' Max(X) =  ',xmax,' index=',ixmax
  end if
!
!  Evaluate the residual of the initial X value.
!
  call fx ( area,eqn,g,indx,ishapf,nelem,neqn,node,np,npar,nparb,nparf,npe, &
    para,phi,res,splflo,syseqn,tauflo,yc)

  irmax = isamax(neqn,res,1)
  rmax = abs(res(irmax))

  rmax0 = rmax
  rmax2 = rmax

  if ( iwrite >= 3 ) then
    write ( *, * ) 'NEWTON - Step ', iter, ' Max(FX) = ',rmax,' index=',irmax
  end if
!
!  Accept the point immediately,
!    if FX is small enough, and
!    this is NOT a finite difference point.
!
!  In this case, the jacobian will not be evaluated and factored.
!
  if ( rmax <= tolnew .and. (.not.lgdif) ) then
    ival = 1
    call imemry('inc','Newton_zero',ival)
    if ( iwrite >= 3 ) then
      write ( *, * ) 'NEWTON - Iterate ',iter,' accepted.'
    end if
    return
  end if
!
!  The initial X value is NOT acceptable.  We must carry out
!  Newton iteration, and attempt to improve it.
!
  do iter = 1, maxnew

    ival = 1
    call imemry('inc','Newton_steps',ival)
!
!  If we have a valid, factored jacobian already, then we may
!  reuse it, if it's not too old, and if we're allowed.
!
    if ( lmat ) then
      if ( jjac == 0 ) then
        lmat = .false.
      else if ( jjac == 1 ) then
        if ( (.not.lgdif) .and. mod(iter-1,ijac) == 0 ) then
          lmat = .false.
        end if
      else if ( jjac == 2 ) then
      else if ( jjac == 3 ) then
        dist = 0.0D+00
        do i = 1,npar
          dist = dist+(para(i)-parjac(i))**2
        end do
        dist = sqrt(dist)
        if ( dist>disjac ) then
          lmat = .false.
        end if
      end if
    end if
!
!  If it's time, evaluate and factor the jacobian.
!
    if ( .not.lmat ) then

      call fprime ( a,area,eqn,g,indx,nelem,neqn,nlband,node,np,npar,nparb, &
        nparf,npe,nrow,para,phi,syseqn)

      call sgbtrf(neqn,neqn,nlband,nlband,a,nrow,ipivot,info)

      ival = 1
      call imemry('inc','Factor_calls',ival)

      if ( info/= 0 ) then

        ival = 1
        call imemry('inc','Newton_fail',ival)

        write ( *, * ) ' '
        write ( *, * ) 'NEWTON - Fatal error!'
        write ( *, * ) '  The jacobian is singular.'
        write ( *, * ) '  SGBTRF returns INFO = ',info
        ierror = 1
        return
      else
        leval = .true.
        lmat = .true.
        call lmemry('set','have_fp',lmat)
        do i = 1,npar
          parjac(i) = para(i)
        end do
      end if

    end if
!
!  Solve the linear system A*DX = RES
!
    call sgbtrs('N',neqn,nlband,nlband,nsys,a,nrow,ipivot,res,neqn,info)

    ival = 1
    call imemry('inc','Solve_calls',ival)
    ival = nsys
    call imemry('inc','Solve_sys',ival)

    if ( info/= 0 ) then

      ival = 1
      call imemry('inc','Newton_fail',ival)

      write ( *, * ) ' '
      write ( *, * ) 'NEWTON - Fatal error!'
      write ( *, * ) '  SGBTRS returns nonzero INFO = ',info
      ierror = 1
      return

    end if

    idmax = isamax ( neqn, res, 1 )
    dmax = abs ( res ( idmax ) )

    if ( iwrite >= 3 ) then
      write ( *, * ) ' '
      write ( *, * ) 'NEWTON - Step ', iter, ' Max(DX) = ',dmax,' Index=',idmax
    end if
!
!  Update the estimated solution G.
!
    lzero = .true.
    do i = 1,nparf
      if ( para(i)/= 0.0D+00 ) then
        lzero=.false.
      end if
    end do

    if ( lzero ) then
      g(1:neqn) = 0.0D+00
    else
      g(1:neqn) = g(1:neqn) - res(1:neqn)
    end if
!
!  Compute the norm of the current X value.
!
    ixmax = isamax ( neqn, g, 1 )
    xmax = abs ( g(ixmax) )

    if ( iwrite >= 3 ) then
      write ( *, * ) ' '
      write ( *, * ) 'NEWTON - Step ', iter, ' Max(X) =  ',xmax,' index=',ixmax
    end if
!
!  Evaluate the residual of the current estimated solution.
!
    call fx ( area,eqn,g,indx,ishapf,nelem,neqn,node,np,npar, &
      nparb,nparf,npe,para,phi,res,splflo,syseqn,tauflo,yc)

    irmax = isamax ( neqn, res, 1 )
    rmax = abs ( res(irmax) )

    if ( iwrite >= 3 ) then
      write ( *, * ) 'NEWTON - Step ', iter, ' Max(FX) = ',rmax,' index=',irmax
    end if
!
!  If RMAX is less than RMAX2, copy current G into G2.
!
    if ( rmax<rmax2 ) then
      xmax2 = xmax
      rmax2 = rmax
      g2(1:neqn) = g(1:neqn)
    end if
!
!  Accept the iterate if the residual is small enough.
!
    if ( rmax <= tolnew ) then
      if ( iwrite >= 3 .or. falter ) then
        write ( *, * ) 'NEWTON - Iterate ',iter,' accepted.'
      end if
      return
    end if
!
!  Reject the iterate if the residual has grown too large.
!
    if ( rmax>10.0D+00*(rmax0+tolnew) .and. iter>1 ) then
      ierror = 1
      write ( *, * ) ' '
      write ( *, * ) 'NEWTON - Warning!'
      write ( *, * ) '  Residual too big on step      ',iter
      write ( *, * ) '  The final stepsize was        ',dmax
      write ( *, * ) '  The initial X norm was        ',xmax0
      write ( *, * ) '  The final X norm was          ',xmax
      write ( *, * ) '  Initial residual  =           ',rmax0
      write ( *, * ) '  Current residual  =           ',rmax
      go to 20
    end if

  end do
!
!  The iteration has failed to converge, or may actually
!  have been terminated early.
!
  ierror = 1

  write ( *, * ) ' '
  write ( *, * ) 'NEWTON - Warning!'
  write ( *, * ) '  No Newton convergence after   ',maxnew,' steps.'
  write ( *, * ) '  The final stepsize was        ',dmax
  write ( *, * ) '  The initial X norm was        ',xmax0
  write ( *, * ) '  The final X norm was          ',xmax
  write ( *, * ) '  The initial residual norm was ',rmax0
  write ( *, * ) '  The final residual norm was   ',rmax

20 continue

  write ( *, * ) ' '
  write ( *, * ) 'NEWTON - Current parameters:'
  write ( *, * ) ' '
  call pr_parameter ( nparb, nparf, para )
!
!  If the Newton process did not converge, and we are using the
!  very miserly jacobian method, then set a flag to evaluate the
!  jacobian and try again.
!
!  But we will only do this if we did not already evaluate the
!  jacobian during this particular process.
!
  if ( jjac == 0 ) then

  else if ( jjac == 1 ) then

  else if ( jjac == 2 ) then

    if ( .not.leval ) then

      write ( *, * ) ' '
      write ( *, * ) 'NEWTON - Note:'
      write ( *, * ) '  Retrying Newton process with new jacobian.'
      g(1:neqn) = g2(1:neqn)
      rmax0 = rmax2
      xmax0 = xmax2

      ival = 1
      call imemry('inc','Newton_falter',ival)

      lmat = .false.
      falter = .true.
      go to 10

    end if

  end if

  write ( *, * ) ' '
  write ( *, * ) 'NEWTON - Note:'
  write ( *, * ) '  The failed Newton process could not be retried.'

  if ( lgdif ) then
    write ( *, * ) '  This is a finite difference point calculation.'
  else
    write ( *, * ) '  This is NOT a finite difference point calculation.'
  end if

  if ( ltarg ) then
    write ( *, * ) '  This is a target point calculation.'
  else
    write ( *, * ) '  This is NOT a target point calculation.'
  end if

  ival = 1
  call imemry('inc','Newton_fail',ival)

  return
end
subroutine node_set ( eqn, ibump, indx, isotri, maxeqn, nelem, neqn, node, &
  np, npe, nx, ny, xbl, xbr )

!*****************************************************************************80
!
!! NODE_SET assigns numbers to the nodes.
!
!  Discussion:
!
!    It numbers the elements, decides which elements shall be
!    isoparametric, (ISOTRI) and assigns six nodes to each (NODE).
!
!    It associates global unknown indices with each node (INDX).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = 2 ) EQN(MAXEQN).
!    EQN records the "type" of each equation that will be generated
!    and, which is associated with an unknown.  Note that most
!    conditions do not result in an equation.  The current values are:
!      'U'  The horizontal momentum equation.
!      'V'  The vertical momentum equation.
!      'P'  The continuity equation.
!      'PB' A pressure boundary condition.
!
!    Input, integer ( kind = 4 ) IBUMP.
!    IBUMP determines where isoparametric elements will be used.
!    0, no isoparametric elements will be used.
!    1, isoparametric elements will be used only for the
!       elements which directly impinge on the bump.
!    2, isoparametric elements will be used for all elements above the bump.
!    3, isoparametric elements will be used for all elements.
!
!    Output, integer ( kind = 4 ) INDX(NP,3).
!    INDX contains, for each node I, the index of U, V and P at
!    that node, or 0 or a negative value.
!    If K = INDX(I,J) is positive, then the value of the degree
!    of freedom is stored in the solution vector entry G(K).
!    If INDX(I,J) is positive, then that means that a degree of
!    freedom for variable J (U, V or P) is associated with node
!    I, and an equation will be generated to determine its value.
!    If INDX(I,J) is 0, then that means the the value of variable
!    J (U, V or P) has been specified at node I.  No equation is
!    generated to determine its value.
!
!    Output, integer ( kind = 4 ) ISOTRI(NELEM).
!    0, the element is NOT isoparametric.  The six node triangle has straight sides.
!    1, the element is isoparametric.  The six node triangle
!       has curved sides.  Many computations involving such an
!       element must be computed by using a reference triangle,
!       and evaluating the jacobian of a transformation between
!       that triangle and the element.
!
!    Input, integer ( kind = 4 ) MAXEQN, the maximum number of equations allowed.
!
!    Input, integer ( kind = 4 ) MX.
!    The number of nodes in the X direction.  MX = 2*NX-1.
!
!    Input, integer ( kind = 4 ) MY.
!    The number of nodes in the Y direction.  MY = 2*NY-1.
!
!    Output, integer ( kind = 4 ) NEQN.
!    The number of equations, and functions.
!
!    Output, integer ( kind = 4 ) NODE(MAXELM,6).
!    NODE contains, for each element, the global indices of the element nodes.
!
!    Input, integer ( kind = 4 ) NP.
!    The number of nodes.  NP = MX*MY=(2*NX-1)*(2*NY-1).
!
!    Input, real ( kind = 8 ) XBL.
!    The X coordinate of the left corner of the bump.
!
!    Input, real ( kind = 8 ) XBR.
!    The X coordinate of the right corner of the bump.
!
  implicit none

  integer ( kind = 4 ) maxeqn
  integer ( kind = 4 ) nelem
  integer ( kind = 4 ) np
  integer ( kind = 4 ) npe

  character ( len = 2 ) eqn(maxeqn)
  integer ( kind = 4 ) ibump
  integer ( kind = 4 ) icol
  integer ( kind = 4 ) icol2
  integer ( kind = 4 ) ielem
  integer ( kind = 4 ) ieqn
  integer ( kind = 4 ) indx(np,3)
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) irow
  integer ( kind = 4 ) irow2
  integer ( kind = 4 ) isotri(nelem)
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) nbleft
  integer ( kind = 4 ) nbrite
  integer ( kind = 4 ) neqn
  integer ( kind = 4 ) node(nelem,npe)
  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny
  real ( kind = 8 ) xbl
  real ( kind = 8 ) xbr
!
!  Tell the call counter that NODE_SET has been called again.
!
  ival = 1
  call imemry('inc','NODE_SET_calls',ival)
!
!  Compute the global node numbers that will be assigned to the
!  beginning and ending of the bump.  These numbers are only used to
!  determine which elements are isoparametric.
!
!
!  Here, we explicitly assume the region is 10.0D+00 units long.
!
  nbleft = nint(xbl*(2*nx-2)/10.0D+00)+1
  nbrite = nint(xbr*(2*nx-2)/10.0D+00)+1
!
!  Consider each of the NP nodes, which logically lie in an MX by MY
!  rectangular array.  A pair of new elements must be generated every
!  time we reach a node that lies in an odd row and column, (except for
!  the top row, and last column, of course).  At every node, we
!  will have to decide how many equations to generate.
!
  ielem = 0
  ieqn = 0

  do ip = 1,np
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
!    element 7: (13, 25, 23, 19, 24, 18)
!    element 8: (13, 15, 25, 14, 20, 19)
!
    if ( irow2 == 1 .and. icol2 == 1 .and. &
         icol /= 2*nx-1 .and. irow /= 2*ny-1 ) then

      ielem = ielem+1

      node(ielem,1) = ip
      node(ielem,2) = ip+2
      node(ielem,3) = ip+2*(2*ny-1)+2
      node(ielem,4) = ip+1
      node(ielem,5) = ip+(2*ny-1)+2
      node(ielem,6) = ip+(2*ny-1)+1

      if ( ibump == 0 ) then

        isotri(ielem) = 0

      else if ( ibump == 1 ) then

        if ( icol >= nbleft .and. icol < nbrite ) then
          isotri(ielem) = 1
        else
          isotri(ielem) = 0
        end if

      else if ( ibump == 2 ) then

        if ( icol >= nbleft .and. icol < nbrite ) then
          isotri(ielem) = 2
        else
          isotri(ielem) = 0
        end if

      else

        isotri(ielem) = 2

      end if

      ielem = ielem+1

      node(ielem,1) = ip
      node(ielem,2) = ip+2*(2*ny-1)+2
      node(ielem,3) = ip+2*(2*ny-1)
      node(ielem,4) = ip+(2*ny-1)+1
      node(ielem,5) = ip+2*(2*ny-1)+1
      node(ielem,6) = ip+(2*ny-1)

      if ( ibump == 0 ) then

        isotri(ielem) = 0

      else if ( ibump == 1 ) then

        if ( irow == 1 .and. icol >= nbleft .and. icol < nbrite ) then
          isotri(ielem) = 2
        else if ( icol >= nbleft .and. icol<nbrite ) then
          isotri(ielem) = 1
        else
          isotri(ielem) = 0
        end if

      else if ( ibump == 2 ) then

        if ( icol >= nbleft .and. icol<nbrite ) then
          isotri(ielem) = 2
        else
          isotri(ielem) = 0
        end if

      else

        isotri(ielem) = 2

      end if

    end if
!
!  Now determine what equations to associate with this node.
!
!  The node lies on the left hand inflow boundary.
!  The horizontal and vertical velocities are specified.
!
    if ( icol == 1 .and. 1<irow .and. irow<2*ny-1 ) then

      ieqn = ieqn+1
      indx(ip,1) = ieqn
      eqn(ieqn) = 'UI'

      ieqn = ieqn+1
      indx(ip,2) = ieqn
      eqn(ieqn) = 'VI'
!
!  The node lies on the right hand boundary.
!  The horizontal velocity is an unknown, the vertical velocity is zero.
!
    else if ( icol == 2*nx-1 .and. 1<irow.and.irow<2*ny-1 ) then

      ieqn = ieqn+1
      indx(ip,1) = ieqn
      eqn(ieqn) = 'U'

      ieqn = ieqn+1
      indx(ip,2) = ieqn
      eqn(ieqn) = 'VW'
!
!  The node lies on the moving bump surface.
!  The horizontal and vertical velocities are zero.
!
    else if ( irow == 1 .and. icol>nbleft.and.icol<nbrite ) then

      ieqn = ieqn+1
      indx(ip,1) = ieqn
      eqn(ieqn) = 'UB'

      ieqn = ieqn+1
      indx(ip,2) = ieqn
      eqn(ieqn) = 'VB'
!
!  The node lies on a fixed wall.
!  The horizontal and vertical velocities are zero.
!
    else if ( icol == 1 .or. icol==2*nx-1.or.(irow==1 .and. icol<=nbleft).or. &
      (irow == 1 .and. icol >= nbrite) .or. irow==2*ny-1 ) then

      ieqn = ieqn+1
      indx(ip,1) = ieqn
      eqn(ieqn) = 'UW'

      ieqn = ieqn+1
      indx(ip,2) = ieqn
      eqn(ieqn) = 'VW'
!
!  The node is a normal interior node.
!  The horizontal and vertical velocities are unknown.
!
    else

      ieqn = ieqn+1
      indx(ip,1) = ieqn
      eqn(ieqn) = 'U'

      ieqn = ieqn+1
      indx(ip,2) = ieqn
      eqn(ieqn) = 'V'

    end if
!
!  On nodes in an odd row and column, add a pressure equation.
!
    if ( irow2 == 1 .and. icol2==1 ) then
      ieqn = ieqn+1
      indx(ip,3) = ieqn
      eqn(ieqn) = 'P'
    else
      indx(ip,3) = 0
    end if

  end do
!
!  Set the total number of equations
!
  neqn = ieqn
!
!  The last equation, which is guaranteed to be a pressure equation,
!  is replaced by a pressure boundary condition, associated with
!  an unknown.  (Even though we know this pressure will be zero).
!
  eqn(neqn) = 'PB'
!
!  Make sure we haven't exceeded the maximum number of equations.
!
  if ( neqn > maxeqn ) then
    write ( *, * ) ' '
    write ( *, * ) 'NODE_SET - Fatal error!'
    write ( *, * ) '  Too many unknowns!'
    write ( *, * ) '  The maximum allowed is MAXEQN = ',maxeqn
    write ( *, * ) '  This problem requires NEQN = ',neqn
    stop
  else
    write ( *, * ) ' '
    write ( *, * ) 'NODE_SET - The number of unknowns is NEQN = ',neqn
  end if

  return
end
subroutine nodnab(i,j,my,nab)

!*****************************************************************************80
!
!! NODNAB returns the neighbor index number of node J relative to node I.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) my
  integer ( kind = 4 ) nab

  if ( j == i-2*my-2 ) then
    nab = 1
  else if ( j == i-2*my-1 ) then
    nab = 2
  else if ( j == i-2*my ) then
    nab = 3
  else if ( j == i-my-2 ) then
    nab = 4
  else if ( j == i-my-1 ) then
    nab = 5
  else if ( j == i-my ) then
    nab = 6
  else if ( j == i-my+1 ) then
    nab = 7
  else if ( j == i-2 ) then
    nab = 8
  else if ( j == i-1 ) then
    nab = 9
  else if ( j == i ) then
    nab = 10
  else if ( j == i+1 ) then
    nab = 11
  else if ( j == i+2 ) then
    nab = 12
  else if ( j == i+my-1 ) then
    nab = 13
  else if ( j == i+my ) then
    nab = 14
  else if ( j == i+my+1 ) then
    nab = 15
  else if ( j == i+my+2 ) then
    nab = 16
  else if ( j == i+2*my ) then
    nab = 17
  else if ( j == i+2*my+1 ) then
    nab = 18
  else if ( j == i+2*my+2 ) then
    nab = 19
  else
    nab = 0
  end if

  return
end
subroutine nusen ( area, eqn, f, g, indx, nelem, neqn, node, np, npe, phi )

!*****************************************************************************80
!
!! NUSEN sets the right hand side of the NU_INV sensitivity equation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) F(NEQN), the right hand side of the sensitivity equations
!    associated with the NU_INV parameter.
!
  implicit none

  integer ( kind = 4 ) nelem
  integer ( kind = 4 ) neqn
  integer ( kind = 4 ) np
  integer ( kind = 4 ) npe

  real ( kind = 8 ) ar
  real ( kind = 8 ) area(nelem,3)
  real ( kind = 8 ) dpdx
  real ( kind = 8 ) dpdy
  real ( kind = 8 ) dudx
  real ( kind = 8 ) dudy
  real ( kind = 8 ) dvdx
  real ( kind = 8 ) dvdy
  character ( len = 2 ) eqn(neqn)
  real ( kind = 8 ) f(neqn)
  real ( kind = 8 ) g(neqn)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ielem
  integer ( kind = 4 ) ihor
  integer ( kind = 4 ) indx(np,3)
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) iq
  integer ( kind = 4 ) iquad
  integer ( kind = 4 ) iver
  integer ( kind = 4 ) node(nelem,npe)
  real ( kind = 8 ) p
  real ( kind = 8 ) phi(nelem,3,6,10)
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) wi
!
!  Zero out the right hand side vector.
!
  f(1:neqn) = 0.0D+00
!
!  Approximate the integral by summing over all elements.
!
  do ielem = 1,nelem
!
!  In each element, evaluate the integrand at the quadrature
!  points.
!
    do iquad = 1,3

      ar = area(ielem,iquad)
!
!  Evaluate the state solution at the quadrature point.
!
      call uvalq ( dpdx,dpdy,dudx,dudy,dvdx,dvdy,g,ielem,indx,iquad,nelem, &
        neqn,node,np,npe,p,phi,u,v)
!
!  The only equations that have a contribution from this element
!  are those associated with basis functions for the element.
!  These, in turn, are associated with the nodes of the element.
!
!  So now we consider each node in the element.
!
      do iq = 1,6

        ip = node(ielem,iq)
        wi = phi(ielem,iquad,iq,1)

        ihor = indx(ip,1)
        if ( eqn(ihor) == 'U' ) then
          f(ihor) = f(ihor)-ar*(u*dudx+v*dudy+dpdx)*wi
        end if

        iver = indx(ip,2)
        if ( eqn(iver) == 'V' ) then
          f(iver) = f(iver)-ar*(u*dvdx+v*dvdy+dpdy)*wi
        end if

      end do
    end do
  end do

  return
end
subroutine osolve ( a,area,disjac,dopt,dpara3,dparfd,dparfdc,dparsn,dpdyn, &
  dudyn,dvdyn,dydpn,eqn,etan,etaq,g,gdif,gdifc,gold,gradf,gtar,ibc, &
  idfd,ids,ifds,igrid,igunit,ijac,indx,iopt,ipivot,iplot,ipred,ishapb,ishapf, &
  ismooth,isotri,itunit,itype,ivopt,iwrite,jjac,liv,lv,maxnew,maxstp,nelem, &
  neqn,nlband,node,nopt,np,npar,nparb,nparf,npe,nprof,nrow,numel,nx,ny,para, &
  para1, &
  parjac,parnew,parold,phi,res,sens,splbmp,splflo,stpmax,syseqn,taubmp,tauflo, &
  tolnew,tolopt,vopt,wateb,watep,wateu,watev,wquad,xbl,xbord,xbr,xc,xopt, &
  xprof,xquad,xsin,xsiq,ybl,ybord,ybr,yc,yquad )

!*****************************************************************************80
!
!! OSOLVE carries out the optimization algorithm for a fixed grid,
!  using only functional values for the optimization.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) liv
  integer ( kind = 4 ) lv
  integer ( kind = 4 ) nelem
  integer ( kind = 4 ) neqn
  integer ( kind = 4 ) np
  integer ( kind = 4 ) npar
  integer ( kind = 4 ) nparb
  integer ( kind = 4 ) nparf
  integer ( kind = 4 ) npe
  integer ( kind = 4 ) nrow
  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny

  real ( kind = 8 ) a(nrow,neqn)
  real ( kind = 8 ) area(nelem,3)
  real ( kind = 8 ) cost
  real ( kind = 8 ) costb
  real ( kind = 8 ) costp
  real ( kind = 8 ) costu
  real ( kind = 8 ) costv
  real ( kind = 8 ) disjac
  real ( kind = 8 ) dopt(npar)
  real ( kind = 8 ) dpara3(npar)
  real ( kind = 8 ) dparfd(npar)
  real ( kind = 8 ) dparfdc(npar)
  real ( kind = 8 ) dparsn(npar)
  real ( kind = 8 ) dpdyn(np)
  real ( kind = 8 ) dudyn(np)
  real ( kind = 8 ) dvdyn(np)
  real ( kind = 8 ) dydpn(np,nparb)
  character ( len = 2 ) eqn(neqn)
  real ( kind = 8 ) etan(6)
  real ( kind = 8 ) etaq(3)
  real ( kind = 8 ) flarea
  real ( kind = 8 ) g(neqn)
  real ( kind = 8 ) gdif(neqn,npar)
  real ( kind = 8 ) gdifc(neqn,npar)
  real ( kind = 8 ) gold(neqn)
  real ( kind = 8 ) gradf(neqn,npar)
  real ( kind = 8 ) gtar(neqn)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibc
  integer ( kind = 4 ) idfd
  integer ( kind = 4 ) ids
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ifds
  integer ( kind = 4 ) igrid
  integer ( kind = 4 ) igunit
  integer ( kind = 4 ) ijac
  integer ( kind = 4 ) indx(np,3)
  integer ( kind = 4 ) iopt(npar)
  integer ( kind = 4 ) ipivot(neqn)
  integer ( kind = 4 ) iplot
  integer ( kind = 4 ) ipred
  integer ( kind = 4 ) ishapb
  integer ( kind = 4 ) ishapf
  integer ( kind = 4 ) ismooth
  integer ( kind = 4 ) isotri(nelem)
  integer ( kind = 4 ) itunit
  integer ( kind = 4 ) itype
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) ivopt(liv)
  integer ( kind = 4 ) iwrite
  integer ( kind = 4 ) jjac
  integer ( kind = 4 ) jopt
  integer ( kind = 4 ) maxnew
  integer ( kind = 4 ) maxstp
  integer ( kind = 4 ) nfail
  integer ( kind = 4 ) nlband
  integer ( kind = 4 ) node(nelem,npe)
  integer ( kind = 4 ) nopt
  integer ( kind = 4 ) nprof(2*ny-1)
  integer ( kind = 4 ) numel(np)
  integer ( kind = 4 ) numstp
  real ( kind = 8 ) para(npar)
  real ( kind = 8 ) para1(npar)
  real ( kind = 8 ) parjac(npar)
  real ( kind = 8 ) parnew(npar)
  real ( kind = 8 ) parold(npar)
  real ( kind = 8 ) phi(nelem,3,6,10)
  real ( kind = 8 ) res(neqn)
  real ( kind = 8 ) sens(neqn,npar)
  real ( kind = 8 ) splbmp(4,nparb+2,0:nparb)
  real ( kind = 8 ) splflo(4,nparf+2,0:nparf)
  real ( kind = 8 ) stpmax
  character ( len = 20 ) syseqn
  real ( kind = 8 ) taubmp(nparb+2)
  real ( kind = 8 ) tauflo(nparf+2)
  character ( len = 80 ) title
  real ( kind = 8 ) tolnew
  real ( kind = 8 ) tolopt
  real ( kind = 8 ) vopt(lv)
  real ( kind = 8 ) wateb
  real ( kind = 8 ) watep
  real ( kind = 8 ) wateu
  real ( kind = 8 ) watev
  real ( kind = 8 ) wquad(3)
  real ( kind = 8 ) xbl
  real ( kind = 8 ) xbord(nx)
  real ( kind = 8 ) xbr
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) xopt(npar)
  real ( kind = 8 ) xprof
  real ( kind = 8 ) xquad(nelem,3)
  real ( kind = 8 ) xsin(6)
  real ( kind = 8 ) xsiq(3)
  real ( kind = 8 ) ybl
  real ( kind = 8 ) ybord(ny)
  real ( kind = 8 ) ybr
  real ( kind = 8 ) yc(np)
  real ( kind = 8 ) yquad(nelem,3)
!
!  Initialize some things.
!
  nfail = 0
  numstp = 0
!
!  Set default values for the optimizer.
!
  ival = 2
  call deflt(ival,ivopt,liv,lv,vopt)

  if ( tolopt>0.0D+00 ) then
    vopt(31) = tolopt
    vopt(32) = tolopt
    vopt(33) = tolopt
    vopt(34) = tolopt
    vopt(37) = tolopt
  end if

  ivopt(1) = 12
  ivopt(19) = 0
!
!  Initialize the solution vectors G, GOLD, and the old parameters
!  PAROLD, to zero.
!
  g(1:neqn) = 0.0D+00
  gold(1:neqn) = 0.0D+00
  parnew(1:npar) = para1(1:npar)
  parold(1:npar) = 0.0D+00
!
!     parold(nparf+nparb+1) = 1.0D+00
!
!  Optimization loop
!
  do
!
!  Call the optimizer to get a new set of parameter values, PARNEW.
!
    parold(1:npar) = parnew(1:npar)

    jopt = 0
    do i = 1,npar
      if ( iopt(i) == 1 ) then
        jopt = jopt+1
        xopt(jopt) = parnew(i)
      end if
    end do

    call snoit ( dopt, cost, ivopt, liv, lv, nopt, vopt, xopt )

    ival = 1
    call imemry('inc','Snoit_calls',ival)

    jopt = 0
    do i = 1,npar
      if ( iopt(i) == 1 ) then
        jopt = jopt+1
        parnew(i) = xopt(jopt)
      end if
    end do
!
!  Check return from optimizer.
!
    if ( ivopt(1) >= 3 .and. ivopt(1) <= 8 ) then

      write ( *, * ) ' '
      write ( *, * ) 'OSOLVE - Convergence to a minimizer was achieved!'
      write ( *, * ) '  IVOPT(1) = ',ivopt(1)
      exit

    else if ( ivopt(1)>8 ) then

      write ( *, * ) ' '
      write ( *, * ) 'OSOLVE - Warning!'
      write ( *, * ) '  IVOPT(1) = ',ivopt(1)
      exit

    else if ( ivopt(1) == 1 ) then

      numstp = numstp+1

      call flosol ( a,area,disjac,eqn,etan,etaq,flarea,g,ierror,igrid,ijac, &
        indx,ipivot,ishapb,ishapf,isotri,iwrite,jjac,maxnew,nelem,neqn, &
        nlband,node,np,npar,nparb,nparf,npe,nrow,nx,ny,parnew,parjac,phi,res, &
        splbmp,splflo,syseqn,taubmp,tauflo,tolnew,wquad,xbl,xbord,xbr,xc, &
        xquad,xsin,xsiq,ybl,ybord,ybr,yc,yquad)

!   call solcon(a,area,disjac,dpara3,dpdyn,dudyn,dvdyn,dydpn,eqn,etan, &
!     etaq,flarea,g,gdif,gdifc,gold,gradf,gtar,ibc,idfd,ids,ierror, &
!     ifds,igrid,ijac,indx,iopt,ipivot,ipred,ishapb,ishapf,ismooth,isotri, &
!     itype,iwrite,jjac,maxnew,nelem,neqn,nlband,node,np,npar,nparb,nparf, &
!     npe,nprof,nrow,numel,nx,ny,para,parjac,parnew,parold,phi,res,sens,splbmp, &
!     splflo,stpmax,syseqn,taubmp,tauflo,tolnew,wateb,watep,wateu,watev,wquad, &
!     xbl,xbord,xbr,xc,xquad,xsin,xsiq,ybl,ybord,ybr,yc,yquad)

      if ( ierror /= 0 ) then
        nfail = nfail+1
        if ( nfail <= 1 ) then
          write ( *, * ) ' '
          write ( *, * ) 'OSOLVE - Warning!'
          write ( *, * ) '  SolCon returns IERROR = ',ierror
          write ( *, * ) '  Requesting that SUMIT try a smaller step.'
          ivopt(2) = 1
          cycle
        else
          write ( *, * ) ' '
          write ( *, * ) 'OSOLVE - Fatal error!'
          write ( *, * ) '  SolCon returns IERROR = ',ierror
          exit
        end if
      else
        nfail = 0
      end if
!
!  Compute the cost COST associated with the solution G,
!  determined by the parameters PARNEW.
!
      call get_cost ( cost,costb,costp,costu,costv,g,gtar,indx,ishapb,neqn,np, &
        nparb,nprof,ny,splbmp,taubmp,wateb,watep,wateu,watev,xbl,xbr,ybl,ybr,yc)

      if ( iwrite == 1 ) then
        title = 'Cost'
        call pr_cost1 ( cost, title )
      else if ( iwrite >= 2 ) then
        title = 'Cost'
        call pr_cost2 ( cost, costb, costp,costu,costv,title,wateb,watep, &
          wateu, watev )
      end if
!
!  Print stuff.
!
      if ( ifds /= 0 ) then
        call cost_gradient ( dparfd,g,gtar,indx,ishapb,neqn,np,npar, &
          nparb,nparf,nprof,ny, &
          gdif,splbmp,taubmp,wateb,watep,wateu,watev,xbl,xbr,ybl,ybr,yc)
      end if

      if ( ifds /= 0 ) then
        call cost_gradient ( dparfdc,g,gtar,indx,ishapb,neqn,np,npar, &
          nparb,nparf,nprof,ny, &
          gdifc,splbmp,taubmp,wateb,watep,wateu,watev,xbl,xbr,ybl,ybr,yc)
      end if

      if ( ids /= 0 ) then
        call cost_gradient ( dparsn,g,gtar,indx,ishapb,neqn,np,npar, &
          nparb,nparf,nprof,ny, &
          sens,splbmp,taubmp,wateb,watep,wateu,watev,xbl,xbr,ybl,ybr,yc)
      end if

      title = 'Optimization solution'
      call pr_solution (dpara3,dparfd,dparfdc,dparsn,g,gdif,gdifc,gtar,idfd,ids,ifds, &
        indx,iwrite,neqn,np,npar,nparb,nparf,nprof,numstp,ny,parnew,sens,title,yc)

      if ( itunit /= 0 ) then
        write(itunit,'(1x,6g14.6)')cost,(parnew(i),i = 1,npar)
        write(itunit,'(1x,5g14.6)')(dparsn(i),i = 1,npar)
        write(itunit,'(1x,5g14.6)')(dparfd(i),i = 1,npar)
        write(itunit,'(1x,5g14.6)')(dpara3(i),i = 1,npar)
      end if

      if ( iplot > 0 .and. mod(numstp,iplot) == 0 ) then

        call plot_file_write ( eqn, g, gdif, igunit, indx, isotri, iwrite, nelem, &
          neqn, node, np, npar, npe, nprof, nx, ny, para, sens, xc, xprof, yc )

      end if

      if ( numstp >= maxstp ) then

        if ( iwrite<2 ) then
          write ( *, * ) ' '
          write ( *, * ) 'OSOLVE - Parameters at point ',numstp
          call pr_parameter(nparb,nparf,parnew)

          title = 'Cost'
          call pr_cost1(cost,title)
        end if

        write ( *, * ) ' '
        write ( *, * ) 'OSOLVE - Warning!'
        write ( *, * ) '  Number of steps exceeded.'
        write ( *, * ) ' '
        exit
      end if

    else if ( ivopt(1) == 2 ) then

!       lgdif = .true.
!       call lmemry('set','fd_grad',lgdif)

      call flosol ( a,area,disjac,eqn,etan,etaq,flarea,g,ierror,igrid,ijac, &
        indx,ipivot,ishapb,ishapf,isotri,iwrite,jjac,maxnew,nelem,neqn, &
        nlband,node,np,npar,nparb,nparf,npe,nrow,nx,ny,parnew,parjac,phi,res, &
        splbmp,splflo,syseqn,taubmp,tauflo,tolnew,wquad,xbl,xbord,xbr,xc, &
        xquad,xsin,xsiq,ybl,ybord,ybr,yc,yquad)

!   call solcon(a,area,disjac,dpara3,dpdyn,dudyn,dvdyn,dydpn,eqn,etan, &
!     etaq,flarea,g,gdif,gdifc,gold,gradf,gtar,ibc,idfd,ids,ierror, &
!     ifds,igrid,ijac,indx,iopt,ipivot,ipred,ishapb,ishapf,ismooth,isotri, &
!     itype,iwrite,jjac,maxnew,nelem,neqn,nlband,node,np,npar,nparb,nparf, &
!     npe,nprof,nrow,numel,nx,ny,para,parjac,parnew,parold,phi,res,sens,splbmp, &
!     splflo,stpmax,syseqn,taubmp,tauflo,tolnew,wateb,watep,wateu,watev,wquad, &
!     xbl,xbord,xbr,xc,xquad,xsin,xsiq,ybl,ybord,ybr,yc,yquad)

!       lgdif = .false.
!       call lmemry('set','fd_grad',lgdif)

      if ( ierror /= 0 ) then
        write ( *, * ) ' '
        write ( *, * ) 'OSOLVE - Fatal error!'
        write ( *, * ) '  SolCon returned IERROR = ',ierror
        stop
      end if

      call get_cost ( cost, costb, costp, costu, costv, g, gtar, indx, &
        ishapb, neqn, np, nparb, nprof, ny, splbmp, taubmp, wateb, watep, &
        wateu, watev, xbl, xbr, ybl, ybr, yc )

    else

      write ( *, * ) ' '
      write ( *, * ) 'OSOLVE - Warning!'
      write ( *, * ) '  Unknown value of IVOPT(1) = ',ivopt(1)
      exit

    end if

  end do

  if ( iplot < 0 ) then

    call plot_file_write ( eqn, g, gdif, igunit, indx, isotri, iwrite, nelem, &
      neqn, node, np, npar, npe, nprof, nx, ny, para, sens, xc, xprof, yc )

  end if

  return
end
subroutine pldx ( nvec, xval, xvec, yder, yvec )

!*****************************************************************************80
!
!! PLDX evaluates the derivative of a piecewise linear function.
!
!  Discussion:
!
!    Note that if XVAL falls to the left of XVEC(1), then YDER = 0,
!    and similarly, if XDER is greater than XVEC(NVEC), YVAL = 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NVEC, the number of abscissas and coefficients
!    that define the piecewise linear.  NVEC must be odd, and
!    at least 3.
!
!    Input, real ( kind = 8 ) XVAL, the point at which the
!    derivative with respect to X is to be evaluated.
!
!    Input, real ( kind = 8 ) XVEC(NVEC), the abscissas of the
!    function.  These should be distinct and in ascending order.
!
!    Output, real ( kind = 8 ) YDER, the value of the derivative of
!    the piecewise linear function with respect to X, at the point
!    XVAL.
!
!    Input, real ( kind = 8 ) YVEC(NVEC), the value of the
!    piecewise linear function at each of the abscissas.
!
  implicit none

  integer ( kind = 4 ) nvec

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ival
  real ( kind = 8 ) xval
  real ( kind = 8 ) xvec(nvec)
  real ( kind = 8 ) yder
  real ( kind = 8 ) yvec(nvec)
!
!  Step 1: Check if XVAL lies outside the intervals.
!
  if ( xval <= xvec(1) ) then
    yder = 0
    return
  else if ( xval >= xvec(nvec) ) then
    yder = 0
    return
  end if
!
!  Step 2: Find index I so that XVEC(I)  <=  XVAL < XVEC(I+1)
!
  do i = 1,nvec-1

    if ( xvec(i)  <=  xval .and. xval <= xvec(i+1) ) then
      ival = i
      go to 10
    end if

  end do

  write ( *, * ) ' '
  write ( *, * ) 'PLDX - Fatal error!'
  write ( *, * ) '  Could not bracket XVAL = ',xval
  stop

10    continue
!
!  Step 3: Evaluate the slope of the linear function at XVAL.
!
  i = ival

  yder = (yvec(i+1)-yvec(i))/(xvec(i+1)-xvec(i))

  return
end
subroutine pldx1(ivec,nvec,xval,xvec,yder)

!*****************************************************************************80
!
!! PLDX1 evaluates the X derivative of a piecewise linear basis function.
!
!  Discussion:
!
!    The IVEC-th basis function is 1 at the IVEC-th node and 0 at the others.
!
!    Note that if XVAL falls to the left of XVEC(1), then YDER = 0,
!    and similarly, if XVAL is greater than XVEC(NVEC), YDER = 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IVEC, the coefficient with respect to which
!    the partial derivative is desired.
!
!    Input, integer ( kind = 4 ) NVEC, the number of abscissas and coefficients
!    that define the piecewise linear.  NVEC must be odd, and
!    at least 3.
!
!    Input, real ( kind = 8 ) XVAL, the point at which the function
!    is to be evaluated.
!
!    Input, real ( kind = 8 ) XVEC(NVEC), the abscissas of the
!    function.  These should be distinct and in ascending order.
!
!    Output, real ( kind = 8 ) YDER, the value of the derivative of
!    the piecewise linear function at the point XVAL.
!
  implicit none

  integer ( kind = 4 ) nvec

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) ivec
  real ( kind = 8 ) xval
  real ( kind = 8 ) xvec(nvec)
  real ( kind = 8 ) yder
!
!  Step 1: Check if XVAL lies outside the intervals.
!
  if ( xval <= xvec(1) ) then
    yder = 0.0D+00
    return
  else if ( xval >= xvec(nvec) ) then
    yder = 0.0D+00
    return
  end if
!
!  Step 2: Find index I so that XVEC(I)  <=  XVAL < XVEC(I+1)
!
  do i = 1,nvec-1

    if ( xvec(i) <= xval .and. xval<=xvec(i+1) ) then
      ival = i
      go to 10
    end if

  end do

  write ( *, * ) ' '
  write ( *, * ) 'PLDX1 - Fatal error!'
  write ( *, * ) '  Could not bracket XVAL = ',xval
  stop

10    continue
!
!  Step 3: Evaluate the slope of the IVEC-th linear function at XVAL.
!
  i = ival
  if ( ival == ivec ) then
    yder = (0.0D+00 - 1.0D+00 )/(xvec(ival+1)-xvec(ival))
  else if ( ival+1 == ivec ) then
    yder = (1.0D+00 - 0.0D+00 )/(xvec(ival+1)-xvec(ival))
  else
    yder = 0.0D+00
  end if

  return
end
subroutine plot_file_open ( plot_file, igunit, iplot )

!*****************************************************************************80
!
!! PLOT_FILE_OPEN opens the plotting file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) PLOT_FILE, the name of the file into which
!    the DISPLAY graphics information will be stored.
!
!    Input/output, integer ( kind = 4 ) IGUNIT.
!    On input, if IGUNIT is zero, then the routine believes
!    that the graphics unit has not yet been opened.
!    If the FORTRAN unit has already been opened, then IGUNIT
!    should be nonzero, and the routine will know not to try
!    to open the file, since it is already open.
!    On output, IGUNIT is the FORTRAN unit used for writing data
!    to the plotfile PLOT_FILE.
!
!    Input/output, integer ( kind = 4 ) IPLOT.
!    On input, IPLOT has been set by the user.  On output, IPLOT
!    may have been reset to 0, if the graphics file could not
!    be opened.
!    IPLOT controls whether or not graphics files
!    suitable for use with DISPLAY will be created.
!    IPLOT = 0 means no such file will be created.
!    IPLOT>0 means plot data will be generated for each step
!    which is evenly divisible by IPLOT, or which is less than
!    or equal to IPLOT.
!    IPLOT = -1 means plot data will be generated for the target, the
!    first and the last steps only.
!
  implicit none

  character ( len = * ) plot_file
  integer ( kind = 4 ) igunit
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iplot

  if ( iplot == 0 ) then
    return
  end if

  if ( igunit == 0 ) then

    write ( *, * ) ' '
    write ( *, * ) 'PLOT_FILE_OPEN - Note:'
    write ( *, * ) '  Opening the plot file "'// trim ( plot_file ) // '".'
    write ( *, * ) ' '

    call get_unit ( igunit )

    open ( unit = igunit, file = plot_file, status = 'replace', &
      form = 'formatted', access = 'sequential', iostat = ios )

    if ( ios /= 0 ) then
      write ( *, * ) ' '
      write ( *, * ) 'PLOT_FILE_OPEN - Warning!'
      write ( *, * ) '  The plot file could not be opened.'
      write ( *, * ) '  Resetting IPLOT to 0.'
      iplot = 0
    end if

  else

    write ( *, * ) ' '
    write ( *, * ) 'PLOT_FILE_OPEN - Note'
    write ( *, * ) '  The plot file is already open.'
    write ( *, * ) '  New information will be appended to it.'

  end if

  return
end
subroutine plot_file_write ( eqn, g, gdif, igunit, indx, isotri, iwrite, nelem, &
  neqn, node, np, npar, npe, nprof, nx, ny, para, sens, xc, xprof, yc )

!*****************************************************************************80
!
!! PLOT_FILE_WRITE writes geometry and solution information to a plot file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NPE, the number of nodes per element.
!
  implicit none

  integer ( kind = 4 ) nelem
  integer ( kind = 4 ) neqn
  integer ( kind = 4 ) np
  integer ( kind = 4 ) npar
  integer ( kind = 4 ) npe
  integer ( kind = 4 ) ny

  character ( len = 2 ) ctemp
  character ( len = 2 ) eqn(neqn)
  real ( kind = 8 ) g(neqn)
  real ( kind = 8 ) gdif(neqn,npar)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icheck
  integer ( kind = 4 ) igunit
  integer ( kind = 4 ) ihor
  integer ( kind = 4 ) indx(np,3)
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) ipar
  integer ( kind = 4 ) iprs
  integer ( kind = 4 ), save :: iset = 0
  integer ( kind = 4 ) isotri(nelem)
  integer ( kind = 4 ) iuk
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) iver
  integer ( kind = 4 ) iwrite
  integer ( kind = 4 ) j
  integer ( kind = 4 ) node(nelem,npe)
  integer ( kind = 4 ) nprof(2*ny-1)
  integer ( kind = 4 ) nsen
  integer ( kind = 4 ) nx
  real ( kind = 8 ) para(npar)
  real ( kind = 8 ) rtemp
  real ( kind = 8 ) rtemp2
  real ( kind = 8 ) sens(neqn,npar)
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) xprof
  real ( kind = 8 ) yc(np)
!
!  Is this right?
!
  nsen = npar
!
!  Record number of calls.
!
  ival = 1
  call imemry ( 'inc', 'PltWrt_calls', ival )

  iset = iset+1
!
!  Number of elements, nodes, parameters,
!  elements in the X direction, elements in the Y direction.
!
  write ( igunit, * ) nelem
  write ( igunit, * ) np
  write ( igunit, * ) npar
  write ( igunit, * ) npe
  write ( igunit, * ) nsen
  write ( igunit, * ) nx
  write ( igunit, * ) ny
!
!  Pressures, P.
!
  do ip = 1, np
    iprs = indx(ip,3)
    if ( iprs <= 0 ) then
      rtemp = 0.0D+00
    else
      rtemp = g(iprs)
    end if
    write ( igunit, * ) rtemp
  end do
!
!  Horizontal velocities, U.
!
  do ip = 1, np
    ihor = indx(ip,1)
    write ( igunit, * ) g(ihor)
  end do
!
!  Vertical velocities, V
!
  do ip = 1, np

    iver = indx(ip,2)
    write ( igunit, * ) g(iver)
  end do
!
!  Indicator of element type (isoparametric or not).
!
  do i = 1, nelem
    write ( igunit, * ) isotri(i)
  end do
!
!  Nodes that make up each element.
!
  do j = 1, npe
    do i = 1, nelem
      write ( igunit, * ) node(i,j)
    end do
  end do
!
!  Indices of the nodes along the profile line.
!
  do i = 1, 2*ny-1
    write ( igunit, * ) nprof(i)
  end do
!
!  Parameters.
!
  do i = 1,npar
    write ( igunit, * ) para(i)
  end do
!
!  Pressure sensitivities, dP/dpar
!
  do ipar = 1,npar

    do ip = 1,np
      iuk = 3
      iprs = indx(ip,iuk)
      rtemp = 0.0D+00
      rtemp2 = 0.0D+00
      if ( iprs > 0 ) then
        rtemp = sens(iprs,ipar)
        rtemp2 = gdif(iprs,ipar)
      end if
      write ( igunit, * ) rtemp, rtemp2
    end do
!
!  Horizontal velocity sensitivities, dU/dpar
!
    do ip = 1,np
      ihor = indx(ip,1)
      write ( igunit, * ) sens(ihor,ipar), gdif(ihor,ipar)
    end do
!
!  Vertical velocity sensitivities, dV/dpar
!
    do ip = 1,np
      iver = indx(ip,2)
      write ( igunit, * ) sens(iver,ipar), gdif(iver,ipar)
    end do

  end do
!
!  X coordinates of nodes.
!
  do i = 1,np
    write ( igunit, * ) xc(i)
  end do
!
!  X coordinate of profile line.
!
  write ( igunit, * ) xprof
!
!  Y coordinates of nodes.
!
  do i = 1, np
    write ( igunit, * ) yc(i)
  end do
!
!  Nodal equation types.
!
  do i = 1,np
    ihor = indx(i,1)
    iver = indx(i,2)
    iprs = indx(i,3)
    if ( iprs <= 0 ) then
      ctemp = '  '
    else
      ctemp = eqn(iprs)
    end if
    write ( igunit, '(3a2)' ) eqn(ihor), eqn(iver), ctemp
  end do
!
!  Write a check at the the end.
!
  icheck = 1953
  write ( igunit, * ) icheck

  if ( iwrite >= 1 ) then
    write ( *, * ) ' '
    write ( *, * ) 'PltWrt wrote data set ',iset,' to file.'
  end if

  return
end
subroutine plval ( nvec, xval, xvec, yval, yvec )

!*****************************************************************************80
!
!! PLVAL evaluates a piecewise linear function at a given point.
!
!  Discussion:
!
!    Note that if XVAL falls to the left of XVEC(1), then YVAL = YVEC(1),
!    and similarly, if XVAL is greater than XVEC(NVEC), YVAL = YVEC(NVEC).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NVEC, the number of abscissas and coefficients
!    that define the piecewise linear.  NVEC must be at least 1.
!
!    Input, real ( kind = 8 ) XVAL, the point at which the function
!    is to be evaluated.
!
!  XVEC   Input, real ( kind = 8 ) XVEC(NVEC), the abscissas of the
!       function.  These should be distinct and in ascending order.
!
!  YVAL   Output, real ( kind = 8 ) YVAL, the value of the piecewise
!         linear function at the point XVAL.
!
!  YVEC   Input, real ( kind = 8 ) YVEC(NVEC), the value of the piecewise
!       function at each of the abscissas.
!
  implicit none

  integer ( kind = 4 ) nvec

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ival
  real ( kind = 8 ) xval
  real ( kind = 8 ) xvec(nvec)
  real ( kind = 8 ) yval
  real ( kind = 8 ) yvec(nvec)
!
!  Step 1: Check if XVAL lies outside the intervals.
!
  if ( xval <= xvec(1) ) then
    yval = yvec(1)
    return
  else if ( xval >= xvec(nvec) ) then
    yval = yvec(nvec)
    return
  end if
!
!  Step 2: Find index I so that XVEC(I)  <=  XVAL < XVEC(I+1)
!
  do i = 1,nvec-1

    if ( xvec(i) <= xval .and. xval<=xvec(i+1) ) then
      ival = i
      go to 10
    end if

  end do

  write ( *, * ) ' '
  write ( *, * ) 'PLVAL - Fatal error!'
  write ( *, * ) '  Could not bracket XVAL = ',xval
  stop

10    continue
!
!  Step 3: Evaluate the linear function at XVAL.
!
  i = ival

  if ( xval == xvec(i+1) ) then
    yval = yvec(i+1)
  else if ( xval == xvec(i) ) then
    yval = yvec(i)
  else
    yval = ( yvec(i)*(xvec(i+1)-xval)+yvec(i+1)*(xval-xvec(i)) )  &
      / (xvec(i+1)-xvec(i))
  end if

  return
end
subroutine plval1(ivec,nvec,xval,xvec,yval)

!*****************************************************************************80
!
!! PLVAL1 evaluates the piecewise linear basis function.
!
!  Discussion:
!
!    The IVEC-th basis function is 1 at node IVEC and 0 at the other nodes.
!
!    Note that if XVAL falls to the left of XVEC(1), then YVAL = 0,
!    and similarly, if XVAL is greater than XVEC(NVEC), YVAL = 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IVEC, the coefficient with respect to which
!    the partial derivative is desired.
!
!    Input, integer ( kind = 4 ) NVEC, the number of abscissas and coefficients
!    that define the piecewise linear.  NVEC must be odd, and
!    at least 3.
!
!    Input, real ( kind = 8 ) XVAL, the point at which the function
!    is to be evaluated.
!
!  XVEC   Input, real ( kind = 8 ) XVEC(NVEC), the abscissas of the
!       function.  These should be distinct and in ascending order.
!
!  YDER   Output, real ( kind = 8 ) YDER, the value of the derivative of
!         the piecewise linear function at the point XVAL.
!
  implicit none

  integer ( kind = 4 ) nvec

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) ivec
  real ( kind = 8 ) xval
  real ( kind = 8 ) xvec(nvec)
  real ( kind = 8 ) yval
!
!  Step 1: Check if XVAL lies outside the intervals.
!
  if ( xval <= xvec(1) ) then
    yval = 0.0D+00
    return
  else if ( xval >= xvec(nvec) ) then
    yval = 0.0D+00
    return
  end if
!
!  Step 2: Find index I so that XVEC(I)  <=  XVAL < XVEC(I+1)
!
  do i = 1,nvec-1

    if ( xvec(i) <= xval .and. xval<=xvec(i+1) ) then
      ival = i
      go to 10
    end if

  end do

  write ( *, * ) ' '
  write ( *, * ) 'PLVAL1 - Fatal error!'
  write ( *, * ) '  Could not bracket XVAL = ',xval
  stop

10    continue
!
!  Step 3: Determine the index of the left endpoint of the least and
!  greatest intervals that IVEC can affect.
!
  i = ival
  if ( ival == ivec ) then
    if ( xval == xvec(ival) ) then
      yval = 1.0D+00
    else
      yval = (xvec(ival+1)-xval)/(xvec(ival+1)-xvec(ival))
    end if
  else if ( ival+1 == ivec ) then
    if ( xval == xvec(ival+1) ) then
      yval = 1.0D+00
    else
      yval = (xval-xvec(ival))/(xvec(ival+1)-xvec(ival))
    end if
  else
    yval = 0.0D+00
  end if

  return
end
subroutine ppvalu ( break, coef, l, k, x, jderiv, value )

!*****************************************************************************80
!
!! PPVALU evaluates a piecewise polynomial function or its derivative.
!
!  Discussion:
!
!    PPVALU calculates the value at X of the JDERIV-th derivative of
!    the piecewise polynomial function F from its piecewise
!    polynomial representation.
!
!    The interval index I, appropriate for X, is found through a
!    call to INTERV.  The formula for the JDERIV-th derivative
!    of F is then evaluated by nested multiplication.
!
!    The J-th derivative of F is given by:
!
!      (d**j)f(x) =
!        coef(j+1,i) + h * (
!        coef(j+2,i) + h * (
!        ...
!        coef(k-1,i) + h * (
!        coef(k,i) / (k-j-1) ) / (k-j-2) ... ) / 2 ) / 1
!
!    with
!
!      H=X-BREAK(I)
!
!    and
!
!      i = max( 1 , max( j ,  break(j) <= x , 1 .le. j .le. l ) ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 20009
!
!  Author:
!
!    Original FORTRAN77 version by Carl DeBoor.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer Verlag.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) BREAK(L+1), real COEF(*), integer L, for
!    piecewise polynomial representation of the function F to
!    be evaluated.
!
!    Input, integer ( kind = 4 ) K, the order of the polynomial pieces
!    that make up the function F.  The usual value for
!    K is 4, signifying a piecewise cubic polynomial.
!
!    Input, real ( kind = 8 ) X, the point at which to evaluate F or
!    of its derivatives.
!
!    Input, integer ( kind = 4 ) JDERIV, the order of the derivative to be
!    evaluated.  If JDERIV is 0, then F itself is evaluated,
!    which is actually the most common case.  It is assumed
!    that JDERIV is zero or positive.
!
!    Output, real ( kind = 8 ) VALUE, the value of the JDERIV-th
!    derivative of F at X.
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) l

  real ( kind = 8 ) break(l+1)
  real ( kind = 8 ) coef(k,l)
  real ( kind = 8 ) fmmjdr
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) jderiv
  integer ( kind = 4 ) m
  integer ( kind = 4 ) ndummy
  real ( kind = 8 ) value
  real ( kind = 8 ) x

  value = 0.0D+00

  fmmjdr = k - jderiv
!
!  Derivatives of order K or higher are identically zero.
!
  if ( k <= jderiv ) then
    return
  end if
!
!  Find the index I of the largest breakpoint to the left of X.
!
  call interv ( break, l+1, x, i, ndummy )
!
!  Evaluate the JDERIV-th derivative of the I-th polynomial piece at X.
!
  h = x - break(i)
  m = k

  do

    value = ( value / fmmjdr ) * h + coef(m,i)
    m = m - 1
    fmmjdr = fmmjdr - 1.0D+00

    if ( fmmjdr <= 0.0D+00 ) then
      exit
    end if

  end do

  return
end
subroutine pqdx(nvec,xval,xvec,yder,yvec)

!*****************************************************************************80
!
!! PQDX evaluates the derivative of a piecewise quadratic function.
!
!  Discussion:
!
!    Note that if XDER falls to the left of XVEC(1), then YVAL = 0,
!    and similarly, if XVAL is greater than XVEC(NVEC), YDER = 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NVEC, the number of abscissas and coefficients
!    that define the piecewise quadratic.  NVEC must be odd, and
!    at least 3.
!
!  XVAL   Input, real ( kind = 8 ) XVAL, the point at which the
!         derivative with respect to X is to be evaluated.
!
!  XVEC   Input, real ( kind = 8 ) XVEC(NVEC), the abscissas of the
!       function.  These should be distinct and in ascending order.
!
!  YDER   Output, real ( kind = 8 ) YDER, the value of the derivative
!         of the piecewise  quadratic function with respect to X,
!         at the point XVAL.
!
!  YVEC   Input, real ( kind = 8 ) YVEC(NVEC), the value of the piecewise
!         quadratic function at each of the abscissas.
!
  implicit none

  integer ( kind = 4 ) nvec

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ival
  real ( kind = 8 ) xval
  real ( kind = 8 ) xvec(nvec)
  real ( kind = 8 ) yder
  real ( kind = 8 ) yvec(nvec)
!
!  Step 0: Check data.
!
  if ( nvec<3 ) then
    write ( *, * ) ' '
    write ( *, * ) 'PQDX - Fatal error.'
    write ( *, * ) '  NVEC is ',nvec
    write ( *, * ) '  but NVEC must be at least 3.'
    stop
  end if

  if ( mod(nvec,2)/= 1 ) then
    write ( *, * ) ' '
    write ( *, * ) 'PQDX - Fatal error!'
    write ( *, * ) '  Even value of NVEC = ',nvec
    stop
  end if
!
!  Step 1: Find odd index I so that XVEC(I)  <=  XVAL < XVEC(I+2)
!
  if ( xval <= xvec(1) ) then
    yder = yvec(1)
    return
  else if ( xval >= xvec(nvec) ) then
    yder = yvec(nvec)
    return
  end if

  do i = 1,nvec-2,2

    if ( xvec(i) <= xval .and. xval<=xvec(i+2) ) then
      ival = i
      go to 10
    end if

  end do

  write ( *, * ) ' '
  write ( *, * ) 'PQDX - Fatal error!'
  write ( *, * ) '  Could not bracket XVAL = ',xval
  stop

10    continue
!
!  Step 2: Evaluate the derivative of the quadratic function at XVAL.
!
  i = ival

  yder = yvec(i)*(2*xval-xvec(i+1)-xvec(i+2)) &
         /((xvec(i)-xvec(i+1))*(xvec(i)-xvec(i+2))) &
         +yvec(i+1)*(2*xval-xvec(i)-xvec(i+2)) &
         /((xvec(i+1)-xvec(i))*(xvec(i+1)-xvec(i+2))) &
         +yvec(i+1)*(2*xval-xvec(i)-xvec(i+2)) &
         /((xvec(i+1)-xvec(i))*(xvec(i+1)-xvec(i+2)))

  return
end
subroutine pqdx1(ivec,nvec,xval,xvec,yder)

!*****************************************************************************80
!
!! PQDX1 evaluates the X derivative of the piecewise quadratic basis function.
!
!  Discussion:
!
!    The IVEC-th basis function is 1 at the IVEC-th node and 0 at the others.
!
!    Note that if XVAL falls to the left of XVEC(1), then YDER = 0,
!    and similarly, if XVAL is greater than XVEC(NVEC), YDER = 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IVEC, the coefficient with respect to which
!    the partial derivative is desired.
!
!    Input, integer ( kind = 4 ) NVEC, the number of abscissas and coefficients
!    that define the piecewise quadratic.  NVEC must be odd, and
!    at least 3.
!
!    Input, real ( kind = 8 ) XVAL, the point at which the function
!    be evaluated.
!
!    Input, real ( kind = 8 ) XVEC(NVEC), the abscissas of the
!    function.  These should be distinct and in ascending order.
!
!    Output, real ( kind = 8 ) YDER, the value of the derivative of
!    the piecewise quadratic function at the point XVAL.
!
  implicit none

  integer ( kind = 4 ) nvec

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) ivec
  real ( kind = 8 ) xval
  real ( kind = 8 ) xvec(nvec)
  real ( kind = 8 ) yder
!
!  Step 0: Check data.
!
  if ( nvec<3 ) then
    write ( *, * ) ' '
    write ( *, * ) 'PQDX1 - Fatal error!'
    write ( *, * ) '  NVEC = ',nvec
    write ( *, * ) '  but NVEC must be at least 3.'
    stop
  end if

  if ( mod(nvec,2)/= 1 ) then
    write ( *, * ) ' '
    write ( *, * ) 'PQDX1 - Fatal error!'
    write ( *, * ) '  Even value of NVEC = ',nvec
    stop
  end if
!
!  Step 1: Find odd index I so that XVEC(I)  <=  XVAL < XVEC(I+2)
!
  if ( xval <= xvec(1) ) then
    yder = 0
    return
  else if ( xval >= xvec(nvec) ) then
    yder = 0
    return
  end if

  do i = 1,nvec-2,2

    if ( xvec(i) <= xval .and. xval<=xvec(i+2) ) then
      ival = i
      go to 10
    end if

  end do

  write ( *, * ) ' '
  write ( *, * ) 'PQDX1 - Fatal error!'
  write ( *, * ) '  Could not bracket XVAL = ',xval
  stop

10    continue
!
!  Step 2: Determine the index of the left endpoint of the least and
!  greatest intervals that IVEC can affect.
!
  if ( mod(ivec,2) == 0 ) then
    ilo = ivec-1
    ihi = ivec-1
  else
    ilo = max(ivec-2,1)
    ihi = ivec
  end if
!
!  Step 3: If XVAL is outside of the intervals that IVEC can affect,
!  the derivative is zero.
!
  if ( ival<ilo .or. ival>ihi ) then
    yder = 0
    return
  end if
!
!  Step 3: Evaluate the derivative of the quadratic function at XVAL.
!
  i = ival

  if ( ivec == ival ) then
    yder = (2.0D+00*xval-xvec(i+1)-xvec(i+2)) &
           /((xvec(i)-xvec(i+1))*(xvec(i)-xvec(i+2)))
  else if ( ivec == ival+1 ) then
     yder = (2.0D+00*xval-xvec(i)-xvec(i+2)) &
         /((xvec(i+1)-xvec(i))*(xvec(i+1)-xvec(i+2)))
  else if ( ivec == ival+2 ) then
      yder = (2.0D+00*xval-xvec(i)-xvec(i+1)) &
         /((xvec(i+2)-xvec(i))*(xvec(i+2)-xvec(i+1)))
  else
    write ( *, * ) ' '
    write ( *, * ) 'PQDX1 - Fatal error!'
    write ( *, * ) '  IVEC = ',ivec
    write ( *, * ) '  IVAL = ',ival
  end if

  return
end
subroutine pqval(nvec,xval,xvec,yval,yvec)

!*****************************************************************************80
!
!! PQVAL evaluates a piecewise quadratic function at a given point.
!
!  Discussion:
!
!    The piecewise quadratic is defined by NVEC values, where NVEC
!    is odd, and at least 3.  The function is defined by specifying
!    a list of nodes XVEC(I), and specifying its value YVEC(I) at each
!    node.
!
!    The function will be a quadratic polynomial over each of
!    (NVEC-1)/2 intervals that are made up a set of three consecutive
!    nodes, with the first one odd.  Thus, XVEC(1), XVEC(2) and XVEC(3)
!    lie in the first interval.
!
!    At the odd nodes, the quadratic that defines the function may
!    change, but the function remains continuous there, though not
!    differentiable.
!
!    Note that if XVAL falls to the left of XVEC(1), then YVAL = YVEC(1),
!    and similarly, if XVAL is greater than XVEC(NVEC), YVAL = YVEC(NVEC).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NVEC, the number of abscissas and coefficients
!    that define the piecewise quadratic.
!    NVEC must be odd, and at least 3.
!
!    Input, real ( kind = 8 ) XVAL, the point at which the function
!    is be evaluated.
!
!    Input, real ( kind = 8 ) XVEC(NVEC), the abscissas of the
!    function.  These should be distinct and in ascending order.
!
!    Output, real ( kind = 8 ) YVAL, the value of the piecewise
!    quadratic function at the point XVAL.
!
!    Input, real ( kind = 8 ) YVEC(NVEC), the value of the
!    piecewise quadratic function at each of the abscissas.
!
  implicit none

  integer ( kind = 4 ) nvec

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ival
  real ( kind = 8 ) xval
  real ( kind = 8 ) xvec(nvec)
  real ( kind = 8 ) yval
  real ( kind = 8 ) yvec(nvec)
!
!  Step 0: Check data.
!
  if ( nvec<3 ) then
    write ( *, * ) ' '
    write ( *, * ) 'PQVAL - Fatal error!'
    write ( *, * ) '  Value of NVEC = ',nvec
    write ( *, * ) '  but NVEC must be at least 3.'
    stop
  end if

  if ( mod(nvec,2)/= 1 ) then
    write ( *, * ) ' '
    write ( *, * ) 'PQVAL - Fatal error!'
    write ( *, * ) '  Even value of NVEC = ',nvec
    stop
  end if
!
!  Step 1: Find odd index I so that XVEC(I)  <=  XVAL < XVEC(I+2)
!
  if ( xval <= xvec(1) ) then
    yval = yvec(1)
    return
  else if ( xval >= xvec(nvec) ) then
    yval = yvec(nvec)
    return
  end if

  do i = 1,nvec-2,2

    if ( xvec(i) <= xval .and. xval<=xvec(i+2) ) then
      ival = i
      go to 10
    end if

  end do

  write ( *, * ) ' '
  write ( *, * ) 'PQVal - Fatal error!'
  write ( *, * ) '  Could not bracket XVAL = ',xval
  write ( *, * ) '  There are ',nvec,' nodes.'
  write ( *, * ) '  First node is at ',xvec(1)
  write ( *, * ) '  Last node is at  ',xvec(nvec)

  do i = 1,nvec
    write(*,*)xvec(i)
  end do
  stop

10    continue
!
!  Step 2: Evaluate the quadratic function at XVAL.
!
  i = ival

  yval = yvec(i)*(xval-xvec(i+1)) * (xval-xvec(i+2)) &
         /((xvec(i)-xvec(i+1))*(xvec(i)-xvec(i+2))) &
         +yvec(i+1)*(xval-xvec(i)) * (xval-xvec(i+2)) &
         /((xvec(i+1)-xvec(i))*(xvec(i+1)-xvec(i+2))) &
         +yvec(i+2)*(xval-xvec(i)) * (xval-xvec(i+1)) &
         /((xvec(i+2)-xvec(i))*(xvec(i+2)-xvec(i+1)))

  return
end
subroutine pqval1(ivec,nvec,xval,xvec,yval)

!*****************************************************************************80
!
!! PQVAL1 evaluates the piecewise quadratic basis function.
!
!  Discussion:
!
!    The IVEC-th basis function is 1 at node IVEC and 0 at the other nodes.
!
!    Note that if XVAL falls to the left of XVEC(1), then YVAL = 0,
!    and similarly, if XVAL is greater than XVEC(NVEC), YVAL = 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IVEC, the coefficient with respect to which
!    the partial derivative is desired.
!
!    Input, integer ( kind = 4 ) NVEC, the number of abscissas and coefficients
!    that define the piecewise quadratic.  NVEC must be odd, and
!    at least 3.
!
!    Input, real ( kind = 8 ) XVAL, the point at which the function
!    is to be evaluated.
!
!    Input, real ( kind = 8 ) XVEC(NVEC), the abscissas of the
!    function.  These should be distinct and in ascending order.
!
!    Output, real ( kind = 8 ) YDER, the value of the derivative of
!    the piecewise quadratic function at the point XVAL.
!
  implicit none

  integer ( kind = 4 ) nvec

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) ivec
  real ( kind = 8 ) xval
  real ( kind = 8 ) xvec(nvec)
  real ( kind = 8 ) yval
!
!  Step 0: Check data.
!
  if ( nvec<3 ) then
    write ( *, * ) ' '
    write ( *, * ) 'PQVAL1 - Fatal error!'
    write ( *, * ) '  Value of NVEC is ',nvec
    write ( *, * ) '  but NVEC must be at least 3.'
    stop
  end if

  if ( mod(nvec,2)/= 1 ) then
    write ( *, * ) ' '
    write ( *, * ) 'PQVAL1 - Fatal error!'
    write ( *, * ) '  Even value of NVEC = ',nvec
    stop
  end if
!
!  Step 1: Find odd index I so that XVEC(I)  <=  XVAL < XVEC(I+2)
!
  if ( xval <= xvec(1) ) then
    yval = 0.0D+00
    return
  else if ( xval >= xvec(nvec) ) then
    yval = 0.0D+00
    return
  end if

  do i = 1,nvec-2,2

    if ( xvec(i) <= xval .and. xval<=xvec(i+2) ) then
      ival = i
      go to 10
    end if

  end do

  write ( *, * ) ' '
  write ( *, * ) 'PQVAL1 - Fatal error!'
  write ( *, * ) '  Could not bracket XVAL = ',xval
  stop

10    continue
!
!  Step 2: Determine the index of the left endpoint of the least and
!  greatest intervals that IVEC can affect.
!
  if ( mod(ivec,2) == 0 ) then
    ilo = ivec-1
    ihi = ivec-1
  else
    ilo = max(ivec-2,1)
    ihi = ivec
  end if
!
!  Step 3: If XVAL is outside of the intervals that IVEC can affect,
!  the value is zero.
!
  if ( ival<ilo .or. ival>ihi ) then
    yval = 0
    return
  end if
!
!  Step 3: Evaluate the quadratic function at XVAL.
!
  i = ival

  if ( ivec == ival ) then
    yval = (xval-xvec(i+1)) * (xval-xvec(i+2)) &
           /((xvec(i)-xvec(i+1))*(xvec(i)-xvec(i+2)))
  else if ( ivec == ival+1 ) then
     yval = (xval-xvec(i)) * (xval-xvec(i+2)) &
         /((xvec(i+1)-xvec(i))*(xvec(i+1)-xvec(i+2)))
  else if ( ivec == ival+2 ) then
      yval = (xval-xvec(i)) * (xval-xvec(i+1)) &
         /((xvec(i+2)-xvec(i))*(xvec(i+2)-xvec(i+1)))
  else
    write ( *, * ) ' '
    write ( *, * ) 'PQVAL1 - Fatal error!'
    write ( *, * ) '  IVEC = ',ivec
    write ( *, * ) '  IVAL = ',ival
  end if

  return
end
subroutine pr_bump ( dudyn,dvdyn,eqn,g,ibc,indx,iparb,ishapb,neqn,np,nparb, &
  splbmp,taubmp,xc,yc )

!*****************************************************************************80
!
!! PR_BUMP prints out the boundary conditions for the bump sensitivities.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) neqn
  integer ( kind = 4 ) np

  real ( kind = 8 ) dudyn(np)
  real ( kind = 8 ) dvdyn(np)
  character ( len = 2 ) eqn(neqn)
  real ( kind = 8 ) g(neqn)
  integer ( kind = 4 ) ibc
  integer ( kind = 4 ) ihor
  integer ( kind = 4 ) indx(np,3)
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) iparb
  integer ( kind = 4 ) ishapb
  integer ( kind = 4 ) nparb
  real ( kind = 8 ) shape
  real ( kind = 8 ) splbmp(4,nparb+2,0:nparb)
  real ( kind = 8 ) taubmp(nparb+2)
  real ( kind = 8 ) ubc
  real ( kind = 8 ) vbc
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) yc(np)

  write ( *, * ) ' '
  write ( *, * ) 'PR_BUMP: Boundary conditions on bump:'
  write ( *, * ) ' '

  do ip = 1,np
    ihor = indx(ip,1)
    if ( eqn(ihor) == 'UB' ) then

      write ( *, * ) ' '
      write ( *, * ) 'Node ',ip
      write ( *, * ) 'X, Y = ',xc(ip),yc(ip)

      call bump_bc(dudyn,dvdyn,ip,iparb,ishapb,np,nparb,shape,splbmp,taubmp, &
        ubc,vbc,xc)

      write ( *, * ) 'Shape  = ',shape
      write ( *, * ) 'FE:    ',ubc,vbc

      call bump_bc1(g,indx,ip,iparb,ishapb,neqn,np,nparb,shape,splbmp,taubmp, &
        ubc,vbc,xc,yc)

      write ( *, * ) 'FD2:   ',ubc,vbc

      call bump_bc2(g,indx,ip,iparb,ishapb,neqn,np,nparb,shape,splbmp,taubmp, &
        ubc,vbc,xc,yc)

      write ( *, * ) 'FD3:   ',ubc,vbc

      if ( ibc == 3 ) then
        call bump_bc(dudyn,dvdyn,ip,iparb,ishapb,np,nparb,shape,splbmp,taubmp, &
          ubc,vbc,xc)

        write ( *, * ) 'Regrid:', ubc, vbc
      end if

    end if
  end do

  return
end
subroutine pr_cost1 ( cost, title )

!*****************************************************************************80
!
!! PR_COST1 prints out the current cost function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
  implicit none

  real ( kind = 8 ) cost
  character ( len = * ) title

  if ( len_trim ( title ) > 0  ) then
    write ( *, '(a,a,1x,g11.4)' ) 'PR_COST1: ', trim ( title ), cost
  else
    write ( *, '(a,1x,g11.4)' ) 'PR_COST1: Cost:', cost
  end if

  return
end
subroutine pr_cost2 ( cost, costb, costp, costu, costv, title, wateb, watep, &
  wateu, watev )

!*****************************************************************************80
!
!! PR_COST2 prints out the current cost function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
  implicit none

  real ( kind = 8 ) cost
  real ( kind = 8 ) costb
  real ( kind = 8 ) costp
  real ( kind = 8 ) costu
  real ( kind = 8 ) costv
  character ( len = * ) title
  real ( kind = 8 ) wateb
  real ( kind = 8 ) watep
  real ( kind = 8 ) wateu
  real ( kind = 8 ) watev

  write ( *, * ) ' '
  write ( *, * ) 'PR_COST2:'

  if ( len_trim ( title ) > 0 ) then
    write ( *, '(a)' ) trim ( title )
    write ( *, * ) ' '
  end if

  write ( *, '(a)' ) '                      B         P' // &
    '         U         V'

  write ( *, '(a,g11.4)' ) '  Cost function:   ', cost

  write ( *, '(a,4g10.3)' ) '  Component costs: ', costb, costp, costu, costv

  write ( *, '(a,4g10.3)' ) '  Weighted costs:  ', wateb * costb, &
    watep * costp, wateu * costu, watev * costv

  return
end
subroutine pr_cost_sen ( dpara, npar, title )

!*****************************************************************************80
!
!! PR_COST_SEN prints out the cost sensitivities.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) DPARA(NPAR), an estimate for the derivatives of the cost
!    function with respect to the various parameters.
!
!      DPARA(I) = D cost / D parameter(I).
!
!    Input, integer ( kind = 4 ) NPAR, the number of parameters.  NPAR = NPARF + NPARB + 1.
!
  implicit none

  integer ( kind = 4 ) npar

  real ( kind = 8 ) dpara(npar)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  character ( len = * ) title

  write ( *, * ) ' '
  if ( len_trim ( title ) > 0 ) then
    write ( *, '(a)' ) trim ( title )
    write ( *, * ) ' '
  end if

  do ilo = 1, npar, 5
    ihi = min ( ilo+4, npar )
    write ( *, '(5g13.5)' ) dpara(ilo:ihi)
  end do

  return
end
subroutine pr_dat(disjac,plot_file,march_file,ibc,ibump,idfd,ids,ifds,igrad, &
  igrid,ijac,iopt,iplot,ipred,ishapb,ishapbt,ishapf,ishapft,ismooth,istep1, &
  istep2,itar,itype,iwrite,jjac,jstep1,jstep2,maxnew,maxpar,maxstp,nopt,npar, &
  nparb,nparf,nstep3,nx,ny,para1,para2,para3,partar,stpmax,syseqn,tolnew, &
  tolopt,wateb,wateb1,wateb2,watep,wateu,watev,xbleft,xbltar,xbord,xbrite, &
  xbrtar,xprof,ybleft,ybltar,ybord,ybrite,ybrtar)

!*****************************************************************************80
!
!! PR_DAT prints the user input file data.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) maxpar

  real ( kind = 8 ) disjac
  character ( len = * ) plot_file
  character ( len = * ) march_file
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibc
  integer ( kind = 4 ) ibump
  integer ( kind = 4 ) idfd
  integer ( kind = 4 ) ids
  integer ( kind = 4 ) ifds
  integer ( kind = 4 ) igrad
  integer ( kind = 4 ) igrid
  integer ( kind = 4 ) ijac
  integer ( kind = 4 ) iopt(maxpar)
  integer ( kind = 4 ) iplot
  integer ( kind = 4 ) ipred
  integer ( kind = 4 ) ishapb
  integer ( kind = 4 ) ishapbt
  integer ( kind = 4 ) ishapf
  integer ( kind = 4 ) ishapft
  integer ( kind = 4 ) ismooth
  integer ( kind = 4 ) istep1
  integer ( kind = 4 ) istep2
  integer ( kind = 4 ) itar
  integer ( kind = 4 ) itype
  integer ( kind = 4 ) iwrite
  integer ( kind = 4 ) jjac
  integer ( kind = 4 ) jstep1
  integer ( kind = 4 ) jstep2
  integer ( kind = 4 ) maxnew
  integer ( kind = 4 ) maxstp
  integer ( kind = 4 ) nopt
  integer ( kind = 4 ) npar
  integer ( kind = 4 ) nparb
  integer ( kind = 4 ) nparf
  integer ( kind = 4 ) nstep3
  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny
  real ( kind = 8 ) para1(maxpar)
  real ( kind = 8 ) para2(maxpar)
  real ( kind = 8 ) para3(maxpar)
  real ( kind = 8 ) partar(maxpar)
  real ( kind = 8 ) stpmax
  character ( len = 20 ) syseqn
  real ( kind = 8 ) tolnew
  real ( kind = 8 ) tolopt
  character ( len = 6 ) type
  real ( kind = 8 ) wateb
  real ( kind = 8 ) wateb1
  real ( kind = 8 ) wateb2
  real ( kind = 8 ) watep
  real ( kind = 8 ) wateu
  real ( kind = 8 ) watev
  real ( kind = 8 ) xbleft
  real ( kind = 8 ) xbltar
  real ( kind = 8 ) xbord(nx)
  real ( kind = 8 ) xbrite
  real ( kind = 8 ) xbrtar
  real ( kind = 8 ) xprof
  real ( kind = 8 ) ybleft
  real ( kind = 8 ) ybltar
  real ( kind = 8 ) ybord(ny)
  real ( kind = 8 ) ybrite
  real ( kind = 8 ) ybrtar
  character ( len = 3 ) yesno

  write ( *, * ) ' '
  write ( *, * ) 'PR_DAT:'
  write ( *, * ) '  Values of user-definable variables:'
  write ( *, * ) ' '
  write ( *, * ) '  Number of horizontal elements, NX = ',nx
  write ( *, * ) '  Number of vertical elements, NY =   ',ny
  write ( *, * ) ' '
  write ( *, * ) '  NPARF = ',nparf,' inflow parameters'
  write ( *, * ) '  NPARB = ',nparb,' bump parameters'
  write ( *, * ) '  NPARR = 1 inverse viscosity parameters.'
  write ( *, * ) '  NPAR =  ',npar,' total parameters'
  write ( *, * ) ' '
  write ( *, * ) 'Parameter  Type  Free to Vary?'
  write ( *, * ) ' '
  do i = 1,npar
    if ( i <= nparf ) then
      type = 'Inflow'
    else if ( i <= nparf+nparb ) then
      type = 'Shape'
    else
      type = 'NU_INV'
    end if
    if ( iopt(i) == 0 ) then
      yesno = 'No'
    else
      yesno = 'Yes'
    end if
    write(*,'(1x,i5,2x,a6,2x,a3)')i,type,yesno
  end do

  write ( *, * ) ' '
  write ( *, * ) 'NOPT = ',nopt,' optimization parameters.'
  write ( *, * ) ' '
!
!  Target data.
!
  if ( itar == 0 ) then

    write ( *, * ) ' '
    write ( *, * ) 'ITAR = 0:'
    write ( *, * ) '  The target data is computed from a flow.'
    write ( *, * ) ' '

    write ( *, * ) 'ISHAPBT = ',ishapbt
    if ( ishapbt == 1 ) then
      write ( *, * ) '  Target bump modeled by C0 linear splines.'
    else if ( ishapbt == 2 ) then
      write ( *, * ) '  Target bump modeled by C0 quadratic splines.'
    else if ( ishapbt == 3 ) then
      write ( *, * ) '  Target bump modeled by C1 cubic splines.'
    end if

    write ( *, * ) 'ISHAPFT = ',ishapft
    if ( ishapft == 1 ) then
      write ( *, * ) '  Target inflow modeled by C0 linear splines.'
    else if ( ishapft == 2 ) then
      write ( *, * ) '  Target inflow modeled by C0 quadratic splines.'
    else if ( ishapft == 3 ) then
      write ( *, * ) '  Target inflow modeled by C1 cubic splines.'
    end if

    if ( nparb > 0 ) then
      write ( *, * ) ' '
      write ( *, * ) '  Target bump from x,y =   ',xbltar,ybltar
      write ( *, * ) '                to x,y =   ',xbrtar,ybrtar
    end if

    write ( *, * ) ' '
    write ( *, * ) '  PARTAR, the parameter values at target point:'

    call pr_parameter ( nparb, nparf, partar )

  else if ( itar == 1 ) then

    write ( *, * ) ' '
    write ( *, * ) 'ITAR = 1:'
    write ( *, * ) '  The target data is an arbitrary formula.'

  end if
!
!  Feasible data.
!
  write ( *, * ) ' '
  write ( *, * ) 'ISHAPB = ',ishapb
  if ( ishapb == 1 ) then
    write ( *, * ) '  Bump modeled by C0 piecewise linears.'
  else if ( ishapb == 2 ) then
    write ( *, * ) '  Bump modeled by C0 piecewise quadratics.'
  else if ( ishapb == 3 ) then
    write ( *, * ) '  Bump modeled by C1 cubic splines.'
  end if

  write ( *, * ) ' '
  write ( *, * ) 'ISHAPF = ',ishapf
  if ( ishapf == 1 ) then
    write ( *, * ) '  Inflow modeled by C0 piecewise linears.'
  else if ( ishapf == 2 ) then
    write ( *, * ) '  Inflow modeled by C0 piecewise quadratics.'
  else if ( ishapf == 3 ) then
    write ( *, * ) '  Inflow modeled by C1 cubic splines.'
  end if

  if ( nparb>0 ) then
    write ( *, * ) ' '
    write ( *, * ) '  Feasible bump from x,y =  ',xbleft,ybleft
    write ( *, * ) '                   to x,y = ',xbrite,ybrite
  end if

  write ( *, * ) ' '
  write ( *, * ) '  The flow discrepancy is measured at XPROF = ',xprof
!
!  Cost function.
!
  write ( *, * ) ' '
  write ( *, * ) ' '
  write ( *, * ) '  Cost function is weighted sum of these costs:'
  write ( *, * ) ' '
  write ( *, * ) ' '
  write ( *, * ) '  U, V and P discrepancies at profile line,'
  write ( *, * ) '  SQRT(Integral (InflowSlope)**2 ),'
  write ( *, * ) '  SQRT(Integral (BumpSlope-LineSlope)**2).'
  write ( *, * ) ' '
  write ( *, * ) 'Weight factors:'
  write ( *, * ) ' '
  write ( *, * ) '  Bump control cost,   WATEB =  ',wateb
  write ( *, * ) '  Pressure discrepancy, WATEP = ',watep
  write ( *, * ) '  U discrepancy, WATEU =        ',wateu
  write ( *, * ) '  V discrepancy, WATEV =        ',watev
!
  write ( *, * ) ' '
  write ( *, * ) 'ISMOOTH = ',ismooth
  if ( ismooth == 0 ) then
    write ( *, * ) '  Use averaged dUdY, dVdY, dPdY.'
  else if ( ismooth == 1 ) then
    write ( *, * ) '  This option has been deleted!'
    stop
  else if ( ismooth == 2 ) then
    write ( *, * ) '  This option has been deleted!'
    stop
  else if ( ismooth == 3 ) then
    write ( *, * ) '  Use averaged dUdY, dVdY, dPdY,'
    write ( *, * ) '  then least-squares smooth all nodes.'
  else if ( ismooth == 4 ) then
    write ( *, * ) '  Use averaged dUdY, dVdY, (NO dPdY!),'
    write ( *, * ) '  then least-squares smooth all nodes.'
  end if

  write ( *, * ) ' '
  if ( igrid == 0 ) then
    write ( *, * ) 'IGRID = 0, X and Y nodes are equally spaced.'
  else if ( igrid == 1 ) then
    write ( *, * ) 'IGRID = 1, the user supplies X and Y node spacings:'
    write ( *, * ) ' '
    do i = 1,nx
      write(*,*)i,' XBord(I) = ',xbord(i)
    end do
    do i = 1,ny
      write(*,*)i,' YBord(I) = ',ybord(i)
    end do
  else if ( igrid == 2 ) then
    write ( *, * ) 'IGRID = 2, what does that mean?'
  end if
!
!  Data relating to Newton iteration.
!
  write ( *, * ) ' '
  write ( *, * ) ' '
  write ( *, * ) 'Generation of starting point for Newton:'
  write ( *, * ) '  IPRED = ',ipred
  if ( ipred == 0 ) then
    write ( *, * ) '  G(I) = GOLD(I).'
  else if ( ipred == 1 ) then
    write ( *, * ) '  G(I) = GOLD(I) + SENS(I,J) * (PAR(J)-PAROLD(J))'
  else if ( ipred == 2 ) then
    write ( *, * ) '  G(I) = GOLD(I) + GDIF(I,J) * (PAR(J)-PAROLD(J))'
  else if ( ipred == 3 ) then
    write ( *, * ) '  G(I) = GOLD(I)+(GDIF(I,J)-GRADF(I,J))*(PAR(J)-PAROLD(J))'
  end if

  write ( *, * ) ' '

  if ( ijac == 1 ) then
    write ( *, * ) '  IJAC = 1, Jacobian is evaluated on every step.'
  else if ( ijac>1 ) then
    write ( *, * ) '  IJAC = N, Jacobian is evaluated on steps'
    write ( *, * ) '  0, N, 2*N, 3*N, ...'
  end if

  if ( jjac == 0 ) then
    write ( *, * ) '  JJAC = 0, update the jacobian for regular points'
    write ( *, * ) '  and for finite difference points.  A new point'
    write ( *, * ) '  always gets a new jacobian.'
  else if ( jjac == 1 ) then
    write ( *, * ) '  JJAC = 1, do NOT update the jacobian for'
    write ( *, * ) '  finite difference points.  A new regular point'
    write ( *, * ) '  always gets a new jacobian.'
  else if ( jjac == 2 ) then
    write ( *, * ) '  JJAC = 2, Only update the jacobian when necessary.'
  else if ( jjac == 3 ) then
    write ( *, * ) '  JJAC = 3, evaluate Jacobian when new parameters '
    write ( *, * ) '  differ from last set where Jacobian was'
    write ( *, * ) '  evaluated by at least ',disjac
  end if

  write ( *, * ) ' '
  write ( *, * ) ' '
  write ( *, * ) '  Up to MAXNEW = ',maxnew,' Newton iterations.'
  write ( *, * ) '  TOLNEW, Newton iteration tolerance = ',tolnew
!
!  Data relating to continuation.
!
  write ( *, * ) '  The continuation stepsize is STPMAX = ',stpmax
!
!  Derivative information.
!
  write ( *, * ) ' '
  write ( *, * ) 'Derivative information:'
  write ( *, * ) ' '
  if ( ids == 0 ) then
    write ( *, * ) 'Discretized sensitivities will NOT be computed.'
  else
    write ( *, * ) 'Discretized sensitivities will be computed.'
  end if

  if ( ifds == 0 ) then
    write ( *, * ) 'Finite difference sensitivities will NOT be computed.'
  else
    write ( *, * ) 'Finite difference sensitivities will be computed.'
  end if

  if ( idfd == 0 ) then
    write ( *, * ) 'Direct cost finite differences will NOT be computed.'
  else
    write ( *, * ) 'Direct cost finite differences will be computed.'
  end if
!
!  Cost gradient approximation.
!
  write ( *, * ) ' '
  write ( *, * ) 'Cost gradient approximation option IGRAD = ',igrad
  if ( igrad == 0 ) then
    write ( *, * ) '  No cost gradient approximation is made.'
  else if ( igrad == 1 ) then
    write ( *, * ) '  Chain rule on discretized sensitivities.'
  else if ( igrad == 2 ) then
    write ( *, * ) '  Chain rule on finite coefficient differences.'
  else if ( igrad == 3 ) then
    write ( *, * ) '  Chain rule on corrected finite cost differences.'
  else if ( igrad == 4 ) then
    write ( *, * ) '  Direct finite cost differences.'
  end if
!
!  ITYPE = 1, 1D March
!
  if ( itype == 1 ) then

    write ( *, * ) ' '
    write ( *, * ) 'ITYPE = ', itype
    write ( *, * ) '  This is a 1D march.'
    write ( *, * ) ' '
    write ( *, * ) '  The first parameter set will be associated with'
    write ( *, * ) '  point number          ',istep1
    write ( *, * ) '  The second with point ',istep2
    write ( *, * ) ' '
    write ( *, * ) '  The first parameter set:'
    write ( *, * ) ' '
    call pr_parameter(nparb,nparf,para1)
    write ( *, * ) ' '
    write ( *, * ) '  The second parameter set:'
    write ( *, * ) ' '
    call pr_parameter(nparb,nparf,para2)
!
!  ITYPE = 2, 2D March
!
  else if ( itype == 2 ) then

    write ( *, * ) ' '
    write ( *, * ) 'ITYPE = ', itype
    write ( *, * ) '  This is a 2D march.'
    write ( *, * ) ' '
    write ( *, * ) '  The first parameter set will be associated with'
    write ( *, * ) '  point number          ',istep1,jstep1
    write ( *, * ) '  The third with point  ',istep2,jstep2
    write ( *, * ) ' '
    write ( *, * ) '  The first parameter set:'
    write ( *, * ) ' '
    call pr_parameter(nparb,nparf,para1)
    write ( *, * ) ' '
    write ( *, * ) '  The second parameter set:'
    write ( *, * ) ' '
    call pr_parameter(nparb,nparf,para2)
    write ( *, * ) ' '
    write ( *, * ) '  Third parameter set:'
    write ( *, * ) ' '
    call pr_parameter(nparb,nparf,para3)
!
!  ITYPE = 3, Optimization.
!
  else if ( itype == 3 ) then

    write ( *, * ) ' '
    write ( *, * ) 'ITYPE = ', itype
    write ( *, * ) '  This is an optimization run.'
    write ( *, * ) '  TOLOPT, optimization tolerance = ',tolopt
    write ( *, * ) '  At most MAXSTP = ', maxstp, ' optimization steps will' &
      // ' be used.'
    write ( *, * ) '  The starting parameter set:'
    write ( *, * ) ' '
    call pr_parameter(nparb,nparf,para1)
!
!  ITYPE = 4, 3D March
!
  else if ( itype == 4 ) then

    write ( *, * ) ' '
    write ( *, * ) 'ITYPE = ', itype
    write ( *, * ) '  This is a 3D march.'
    write ( *, * ) ' '
    write ( *, * ) '  The first parameter set will be associated with'
    write ( *, * ) '  point number          ',istep1,jstep1
    write ( *, * ) '  The third with point  ',istep2,jstep2
    write ( *, * ) ' '
    write ( *, * ) '  The first parameter set:'
    write ( *, * ) ' '
    call pr_parameter(nparb,nparf,para1)
    write ( *, * ) ' '
    write ( *, * ) '  The second parameter set:'
    write ( *, * ) ' '
    call pr_parameter(nparb,nparf,para2)
    write ( *, * ) ' '
    write ( *, * ) '  Third parameter set:'
    write ( *, * ) ' '
    call pr_parameter(nparb,nparf,para3)

    write ( *, * ) ' '
    if ( nstep3>1 ) then
      write ( *, * ) '  The penalty parameter WATEB will vary from ', &
        wateb1,' to ',wateb2,' in ',nstep3,' steps.'
    else if ( nstep3 == 1 ) then
      write ( *, * ) '  This is an unusual march, since NSTEP3 = 1.'
      write ( *, * ) '  The penalty parameter WATEB will be set to', &
        wateb1,' and cannot march to ',wateb2
    end if
    write ( *, * ) ' '
    write ( *, * ) '  Functional values will be written to MARCH_FILE:'
    write(*,'(5x,a)') trim ( march_file )
!
!  ITYPE = 5, one step Navier-Stokes solve.
!
  else if ( itype == 5 ) then
    write ( *, * ) ' '
    write ( *, * ) 'ITYPE = ', itype
    write ( *, * ) '  This is a one-step Navier Stokes solution.'
!
!  ITYPE = 7, Optimization.
!
  else if ( itype == 7 ) then

    write ( *, * ) ' '
    write ( *, * ) 'ITYPE = ', itype
    write ( *, * ) '  This is an optimization run.'
    write ( *, * ) '  TOLOPT, optimization tolerance = ',tolopt
    write ( *, * ) '  At most ',maxstp,' optimization steps will'//' be used.'
    write ( *, * ) '  The starting parameters:'
    write ( *, * ) ' '
    call pr_parameter ( nparb, nparf, para1 )
    write ( *, * ) ' '
    write ( *, * ) '  Internal approximations to the gradient are made.'

  end if

  write ( *, * ) ' '
  write ( *, * ) 'IBC = ',ibc
  if ( ibc == 0 ) then
    write ( *, * ) '  Bump sensitivity boundary condition uses'
    write ( *, * ) '  finite element estimate of dUdY, dVdY.'
  else if ( ibc == 1 ) then
    write ( *, * ) '  Bump sensitivity boundary condition uses'
    write ( *, * ) '  two point finite difference estimate of '
    write ( *, * ) '  dUdY, dVdY.'
  else if ( ibc == 2 ) then
    write ( *, * ) '  Bump sensitivity boundary condition uses'
    write ( *, * ) '  three point finite difference estimate of '
    write ( *, * ) '  dUdY, dVdY.'
  else if ( ibc == 3 ) then
    write ( *, * ) '  Bump sensitivity boundary condition uses'
    write ( *, * ) '  data from run on a finer grid.'
  end if

  write ( *, * ) ' '
  if ( ibump == 0 ) then
    write ( *, * ) 'IBUMP = 0:'
    write ( *, * ) '  No isoparametric elements will be used.'
    write ( *, * ) '  Y coordinates of midside nodes above the bump'
    write ( *, * ) '  will be adjusted to preserve straight sides.'
  else if ( ibump == 1 ) then
    write ( *, * ) 'IBUMP = 1:'
    write ( *, * ) '  Isoparametric elements directly on bump.'
    write ( *, * ) '  Y coordinates of midside nodes above the bump'
    write ( *, * ) '  will be adjusted to preserve straight sides.'
  else if ( ibump == 2 ) then
    write ( *, * ) 'IBUMP = 2:'
    write ( *, * ) '  All elements above bump are isoparametric.'
    write ( *, * ) '  Y coordinates of midside nodes above the bump'
    write ( *, * ) '  need not lie on a straight line.'
  else if ( ibump == 3 ) then
    write ( *, * ) 'IBUMP = 3:'
    write ( *, * ) '  All elements are isoparametric.'
    write ( *, * ) '  Y coordinates of midside nodes above the bump'
    write ( *, * ) '  need not lie on a straight line.'
  else
    write ( *, * ) ' '
    write ( *, * ) 'PR_DAT - Fatal error!'
    write ( *, * ) '  Unexpected value of IBUMP = ',ibump
    stop
  end if

  write ( *, * ) ' '
  if ( iplot < 0 ) then
    write ( *, * ) '  DISPLAY graphics written for target and '
    write ( *, * ) '    last points.'
  else if ( iplot == 0 ) then
    write ( *, * ) '  DISPLAY graphics file is not written.'
  else
    write ( *, * ) '  DISPLAY graphics written every ',iplot,' points.'
  end if

  if ( iplot/= 0 ) then
    write ( *, * ) '  Graphics data will be written to the file:'
    write(*,'(4x,a)') trim ( plot_file )
  end if

  write ( *, * ) '  Diagnostic output option, IWRITE = ',iwrite
  write ( *, * ) ' '
  write ( *, * ) '  The flow system to be solved is ' // trim ( syseqn )

  return
end
subroutine pr_disc ( g, gtar, indx, neqn, np, nprof, ny, yc )

!*****************************************************************************80
!
!! PR_DISC prints the discrepancy along the profile line.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) neqn
  integer ( kind = 4 ) np
  integer ( kind = 4 ) ny

  real ( kind = 8 ) g(neqn)
  real ( kind = 8 ) gtar(neqn)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(np,3)
  integer ( kind = 4 ) nprof(2*ny-1)
  real ( kind = 8 ) p
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) yc(np)

  write ( *, * ) ' '
  write ( *, * ) 'PR_DISC:'
  write ( *, * ) '  Discrepancy along the profile line between the'
  write ( *, * ) '  computed solution coefficient vector and the'
  write ( *, * ) '  target solution coefficient vector.'
  write ( *, * ) ' '
  write ( *, * ) 'Node     Y         U         V         P'
  write ( *, * ) ' '
  do i = 1,2*ny-1

    u = g(indx(nprof(i),1)) - gtar(indx(nprof(i),1))
    v = g(indx(nprof(i),2)) - gtar(indx(nprof(i),2))

    if ( indx(nprof(i),3) > 0 ) then
      p = g(indx(nprof(i),3)) - gtar(indx(nprof(i),3))
    else
      p = 0.0D+00
    end if

    if ( indx(nprof(i),3) > 0 ) then
      write(*,'(i5,4g11.3)') nprof(i), yc(nprof(i)), u, v, p
    else
      write(*,'(i5,4g11.3)') nprof(i), yc(nprof(i)), u, v
    end if

  end do

  return
end
subroutine pr_du2(dudyn,g,ic,indx,neqn,np,nx,ny,xc,yc)

!*****************************************************************************80
!
!! PR_DU2 prints a comparison of dUdY and the "improved" version dUdY2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) neqn
  integer ( kind = 4 ) np

  real ( kind = 8 ) dudyn(np)
  real ( kind = 8 ) g(neqn)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) icol(19)
  integer ( kind = 4 ) ie
  integer ( kind = 4 ) iee
  integer ( kind = 4 ) in
  integer ( kind = 4 ) indx(np,3)
  integer ( kind = 4 ) ine
  integer ( kind = 4 ) inee
  integer ( kind = 4 ) inn
  integer ( kind = 4 ) inne
  integer ( kind = 4 ) innee
  integer ( kind = 4 ) inw
  integer ( kind = 4 ) irow(19)
  integer ( kind = 4 ) is
  integer ( kind = 4 ) ise
  integer ( kind = 4 ) iss
  integer ( kind = 4 ) issw
  integer ( kind = 4 ) issww
  integer ( kind = 4 ) isw
  integer ( kind = 4 ) isww
  integer ( kind = 4 ) iw
  integer ( kind = 4 ) iww
  integer ( kind = 4 ) j
  character ( len = 4 ) label(19)
  integer ( kind = 4 ) noddat(19)
  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny
  real ( kind = 8 ) u
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) yc(np)

  irow(10) = mod(ic,2*ny-1)
  icol(10) = 1+(ic-1)/(2*ny-1)

  issww = ic-2*(2*ny-1)-2
  isww = issww+1
  iww = isww+1
  issw = ic-(2*ny-1)-2
  isw = issw+1
  iw = isw+1
  inw = iw+1
  iss = ic-2
  is = iss+1
  in = ic+1
  inn = in+1
  ise = ic+(2*ny-1)-1
  ie = ise+1
  ine = ie+1
  inne = ine+1
  iee = ic+2*(2*ny-1)
  inee = iee+1
  innee = inee+1

  noddat(1) = issww
  noddat(2) = isww
  noddat(3) = iww
  noddat(4) = issw
  noddat(5) = isw
  noddat(6) = iw
  noddat(7) = inw
  noddat(8) = iss
  noddat(9) = is
  noddat(10) = ic
  noddat(11) = in
  noddat(12) = inn
  noddat(13) = ise
  noddat(14) = ie
  noddat(15) = ine
  noddat(16) = inne
  noddat(17) = iee
  noddat(18) = inee
  noddat(19) = innee

  irow(1) = irow(10)-2
  irow(2) = irow(10)-1
  irow(3) = irow(10)
  irow(4) = irow(10)-2
  irow(5) = irow(10)-1
  irow(6) = irow(10)
  irow(7) = irow(10)+1
  irow(8) = irow(10)-2
  irow(9) = irow(10)-1
  irow(11) = irow(10)+1
  irow(12) = irow(10)+2
  irow(13) = irow(10)-1
  irow(14) = irow(10)
  irow(15) = irow(10)+1
  irow(16) = irow(10)+2
  irow(17) = irow(10)
  irow(18) = irow(10)+1
  irow(19) = irow(10)+2

  icol(1) = icol(10)-2
  icol(2) = icol(10)-2
  icol(3) = icol(10)-2
  icol(4) = icol(10)-1
  icol(5) = icol(10)-1
  icol(6) = icol(10)-1
  icol(7) = icol(10)-1
  icol(8) = icol(10)
  icol(9) = icol(10)
  icol(11) = icol(10)
  icol(12) = icol(10)
  icol(13) = icol(10)+1
  icol(14) = icol(10)+1
  icol(15) = icol(10)+1
  icol(16) = icol(10)+1
  icol(17) = icol(10)+2
  icol(18) = icol(10)+2
  icol(19) = icol(10)+2

  label(1) = 'ssww'
  label(2) = 'sww'
  label(3) = 'ww'
  label(4) = 'ssw'
  label(5) = 'sw'
  label(6) = 'w'
  label(7) = 'nw'
  label(8) = 'ss'
  label(9) = 's'
  label(10) = 'c'
  label(11) = 'n'
  label(12) = 'nn'
  label(13) = 'se'
  label(14) = 'e'
  label(15) = 'ne'
  label(16) = 'nne'
  label(17) = 'ee'
  label(18) = 'nee'
  label(19) = 'nnee'

  write ( *, * ) ' '
  write ( *, * ) 'PR_DU2'
  write ( *, * ) '  Values centered at node ',ic
  write ( *, * ) '  in row ',irow(10),' and column ',icol(10)
  write ( *, * ) '  at ',xc(ic),yc(ic)
  write ( *, * ) ' '
  write ( *, * ) '  Local Node, Global Node, Label, dUdY, U, Y'
  write ( *, * ) ' '

  do i = 1,19
    if ( irow(i) >= 1 .and. irow(i)<=2*ny-1.and. &
         icol(i) >= 1 .and. icol(i)<=2*nx-1 ) then
      j = noddat(i)
      u = g(indx(j,1))
      write(*,'(1x,2i8,2x,a4,2x,3g14.6)')i,j,label(i),dudyn(j),u,yc(j)
    end if
  end do

  return
end
subroutine pr_fx3(eqn,indx,neqn,np,res,xc,yc)

!*****************************************************************************80
!
!! PR_FX3 prints out the maximum of the P, U, and V residuals.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) neqn
  integer ( kind = 4 ) np

  character ( len = 2 ) eqn(neqn)
  integer ( kind = 4 ) indx(np,3)
  integer ( kind = 4 ) ipmax
  integer ( kind = 4 ) iumax
  integer ( kind = 4 ) ivmax
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jpmax
  integer ( kind = 4 ) jumax
  integer ( kind = 4 ) jvmax
  integer ( kind = 4 ) k
  real ( kind = 8 ) pmax
  real ( kind = 8 ) res(neqn)
  real ( kind = 8 ) rk
  real ( kind = 8 ) umax
  real ( kind = 8 ) vmax
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) yc(np)

  umax = 0.0D+00
  vmax = 0.0D+00
  pmax = 0.0D+00
  iumax = 1
  ivmax = 1
  ipmax = 1
  jumax = 1
  jvmax = 1
  jpmax = 1

  do j = 1,np

    k = indx(j,1)
    rk = abs(res(k))
    if ( rk >= umax ) then
      umax = rk
      iumax = j
      jumax = k
    end if

    k = indx(j,2)
    rk = abs(res(k))
    if ( rk >= vmax ) then
      vmax = rk
      ivmax = j
      jvmax = k
    end if

    k = indx(j,3)
    if ( k>0 ) then
      rk = abs(res(k))
      if ( rk >= pmax ) then
        pmax = rk
        ipmax = j
        jpmax = k
      end if
    end if

  end do

  write ( *, * ) ' '
  write ( *, * ) 'Variable, EqnType, Eqn #, MaxRes, Node, X, Y'
  write ( *, * ) ' '
  write(*,'(5x,a1,3x,a2,2x,i5,g14.6,i5,2g14.6)') &
    'U',eqn(jumax),jumax,umax,iumax,xc(iumax),yc(iumax)
  write(*,'(5x,a1,3x,a2,2x,i5,g14.6,i5,2g14.6)') &
    'V',eqn(jvmax),jvmax,vmax,ivmax,xc(ivmax),yc(ivmax)
  write(*,'(5x,a1,3x,a2,2x,i5,g14.6,i5,2g14.6)') &
    'P',eqn(jpmax),jpmax,pmax,ipmax,xc(ipmax),yc(ipmax)

  return
end
subroutine pr_gs2 ( flarea,gdif,gdifc,gradf,indx,iopt,neqn,np,npar,ny,sens,xc,yc)

!*****************************************************************************80
!
!! PR_GS2 prints the maximum finite coefficient differences and sensitivities.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) neqn
  integer ( kind = 4 ) np
  integer ( kind = 4 ) npar

  real ( kind = 8 ) flarea
  real ( kind = 8 ) gdif(neqn,npar)
  real ( kind = 8 ) gdifc(neqn,npar)
  real ( kind = 8 ) gradf(neqn,npar)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icol
  integer ( kind = 4 ) ieqn
  integer ( kind = 4 ) indx(np,3)
  integer ( kind = 4 ) iopt(npar)
  integer ( kind = 4 ) ipar
  integer ( kind = 4 ) ipmax
  integer ( kind = 4 ) irow
  integer ( kind = 4 ) iumax
  integer ( kind = 4 ) ivmax
  integer ( kind = 4 ) j
  character ( len = 15 ) label
  integer ( kind = 4 ) npmax
  integer ( kind = 4 ) numax
  integer ( kind = 4 ) nvmax
  integer ( kind = 4 ) ny
  real ( kind = 8 ) pmax
  real ( kind = 8 ) psum
  real ( kind = 8 ) sens(neqn,npar)
  real ( kind = 8 ) umax
  real ( kind = 8 ) usum
  real ( kind = 8 ) val
  real ( kind = 8 ) vmax
  real ( kind = 8 ) vsum
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) yc(np)

  write ( *, * ) ' '
  write ( *, * ) 'PR_GS2:'
  write ( *, * ) '  DS = Discretized Sensitivities;'
  write ( *, * ) '  FCD = Finite Coefficient Differences;'
  write ( *, * ) ' '
  write ( *, * ) '  Compute Max ( DS - FCD )'
  write ( *, * ) ' '
  write ( *, * ) '          Par  Var  L1-Norm  Max-Norm  Index  Node  ' &
    //'X       Y  Row  Column'
  do i = 1,6

    write ( *, * ) ' '
    if ( i == 1 ) then
      label = 'DS'
    else if ( i == 2 ) then
      label = 'FCD'
    else if ( i == 3 ) then
      label = 'FCD-FDFix'
    else if ( i == 4 ) then
      label = 'FDFix'
    else if ( i == 5 ) then
      label = 'FCD-DS'
    else if ( i == 6 ) then
      label = 'FCD-FDFix-DS'
    end if

    do ipar = 1,npar

      if ( iopt(ipar) == 1 ) then

        ieqn = 0

        ipmax = 0
        iumax = 0
        ivmax = 0

        npmax = 0
        numax = 0
        nvmax = 0

        pmax = 0.0D+00
        umax = 0.0D+00
        vmax = 0.0D+00

        psum = 0.0D+00
        usum = 0.0D+00
        vsum = 0.0D+00

        do j = 1,np

          ieqn = ieqn+1

          if ( i == 1 ) then
            val = sens(ieqn,ipar)
          else if ( i == 2 ) then
            val = gdif(ieqn,ipar)
          else if ( i == 3 ) then
            val = gdifc(ieqn,ipar)
          else if ( i == 4 ) then
            val = gradf(ieqn,ipar)
          else if ( i == 5 ) then
            val = sens(ieqn,ipar)-gdif(ieqn,ipar)
          else if ( i == 6 ) then
            val = sens(ieqn,ipar)-gdifc(ieqn,ipar)
          end if

          usum = usum+abs(val)

          if ( abs(val)>umax ) then
            umax = abs(val)
            numax = j
            iumax = ieqn
          end if

          ieqn = ieqn+1

          if ( i == 1 ) then
            val = sens(ieqn,ipar)
          else if ( i == 2 ) then
            val = gdif(ieqn,ipar)
          else if ( i == 3 ) then
            val = gdifc(ieqn,ipar)
          else if ( i == 4 ) then
            val = gradf(ieqn,ipar)
          else if ( i == 5 ) then
            val = sens(ieqn,ipar)-gdif(ieqn,ipar)
          else if ( i == 6 ) then
            val = sens(ieqn,ipar)-gdifc(ieqn,ipar)
          end if

          vsum = vsum+abs(val)

          if ( abs(val) > vmax ) then
            vmax = abs(val)
            nvmax = j
            ivmax = ieqn
          end if

          if ( indx(j,3) > 0 ) then

            ieqn = ieqn+1

            if ( i == 1 ) then
              val = sens(ieqn,ipar)
            else if ( i == 2 ) then
              val = gdif(ieqn,ipar)
            else if ( i == 3 ) then
              val = gdifc(ieqn,ipar)
            else if ( i == 4 ) then
              val = gradf(ieqn,ipar)
            else if ( i == 5 ) then
              val = sens(ieqn,ipar)-gdif(ieqn,ipar)
            else if ( i == 6 ) then
              val = sens(ieqn,ipar)-gdifc(ieqn,ipar)
            end if

            psum = psum+abs(val)

            if ( abs(val) > pmax ) then
              pmax = abs(val)
              npmax = j
              ipmax = ieqn
            end if

          end if

        end do

        numax = max ( numax, 1 )
        nvmax = max ( nvmax, 1 )
        npmax = max ( npmax, 1 )
!
!  Normalize L1 quantities by dividing by current flow area.
!
        usum = usum / flarea
        vsum = vsum / flarea
        psum = psum / flarea

        write ( *, * ) ' '

        icol = ((numax-1)/(2*ny-1))+1
        irow = numax-(icol-1)*(2*ny-1)
        write(*,'(a15,1x,i1,1x,a1,2g11.3,2i8,2f8.3,2i4)') &
          label,ipar,'U',usum,umax,iumax,numax,xc(numax),yc(numax),irow,icol

        icol = ((nvmax-1)/(2*ny-1))+1
        irow = nvmax-(icol-1)*(2*ny-1)
        write(*,'(a15,1x,i1,1x,a1,2g11.3,2i8,2f8.3,2i4)') &
          label,ipar,'V',vsum,vmax,ivmax,nvmax,xc(nvmax),yc(nvmax),irow,icol

        icol = ((npmax-1)/(2*ny-1))+1
        irow = npmax-(icol-1)*(2*ny-1)
        write(*,'(a15,1x,i1,1x,a1,2g11.3,2i8,2f8.3,2i4)') &
          label,ipar,'P',psum,pmax,ipmax,npmax,xc(npmax),yc(npmax),irow,icol

      end if

    end do

  end do

  return
end
subroutine pr_parameter ( nparb, nparf, para )

!*****************************************************************************80
!
!! PR_PARAMETER prints out the current parameters.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) nparb
  integer ( kind = 4 ) nparf

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  real ( kind = 8 ) para(nparf+nparb+1)

  write ( *, * ) ' '

  do ilo = 1, nparf, 5

    ihi = min ( ilo+4, nparf )

    if ( ilo == 1 ) then
      write(*,'(''  Inflow '',5g14.6)') para(ilo:ihi)
    else
      write(*,'(''         '',5g14.6)') para(ilo:ihi)
    end if

  end do

  do ilo = nparf+1,nparf+nparb,5

    ihi = min(nparf+ilo+4,nparf+nparb)

    if ( ilo == nparf+1 ) then
      write(*,'(''  Bump   '',5g14.6)') para(ilo:ihi)
    else
      write(*,'(''         '',5g14.6)') para(ilo:ihi)
    end if

  end do

  write ( *, '(''  NU_INV '',g14.6)' ) para(nparf+nparb+1)

  return
end
subroutine pr_profile ( g, indx, neqn, np, nprof, ny, yc )

!*****************************************************************************80
!
!! PR_PROFILE prints out the solution along the profile line.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) neqn
  integer ( kind = 4 ) np
  integer ( kind = 4 ) ny

  real ( kind = 8 ) g(neqn)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(np,3)
  integer ( kind = 4 ) nprof(2*ny-1)
  real ( kind = 8 ) p
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) yc(np)

  write ( *, * ) ' '
  write ( *, * ) 'PR_PROFILE:'
  write ( *, * ) '  Flow values along the profile line.'
  write ( *, * ) ' '
  write ( *, * ) 'Node     Y         U         V         P'
  write ( *, * ) ' '

  do i = 1, 2*ny-1

    u = g(indx(nprof(i),1))
    v = g(indx(nprof(i),2))

    if ( indx(nprof(i),3) > 0 ) then
      p = g(indx(nprof(i),3))
    else
      p = 0.0D+00
    end if

    if ( indx(nprof(i),3) > 0 ) then
      write ( *, '(i5,4g11.3)' ) nprof(i), yc(nprof(i)), u, v, p
    else
      write ( *, '(i5,4g11.3)' ) nprof(i), yc(nprof(i)), u, v
    end if

  end do

  return
end
subroutine pr_puv ( g, indx, neqn, np, title )

!*****************************************************************************80
!
!! PR_PUV prints the nodal values of the finite differences and sensitivities.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) neqn
  integer ( kind = 4 ) np

  real ( kind = 8 ) g(neqn)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(np,3)
  integer ( kind = 4 ) j
  real ( kind = 8 ) pnorm
  character ( len = * ) title
  real ( kind = 8 ) unorm
  real ( kind = 8 ) vnorm

  pnorm = 0.0D+00
  unorm = 0.0D+00
  vnorm = 0.0D+00

  do i = 1, np

    j = indx(i,1)
    unorm = max ( unorm, abs ( g(j) ) )

    j = indx(i,2)
    vnorm = max ( vnorm, abs ( g(j) ) )

    j = indx(i,3)
    if ( j > 0 ) then
      pnorm = max ( pnorm, abs ( g(j) ) )
    end if

  end do

  write ( *, * ) ' '
  write ( *, * ) 'PR_PUV'
  write ( *, * ) '  Maximum norms of a solution vector:'
  if ( len_trim ( title ) > 0 ) then
    write ( *, '(4x,a)' ) trim ( title )
  end if
  write ( *, * ) ' '
  write ( *, * ) '  U = ', unorm
  write ( *, * ) '  V = ', vnorm
  write ( *, * ) '  P = ', pnorm

  return
end
subroutine pr_solution ( dpara3, dparfd, dparfdc, dparsn, g, gdif, gdifc, &
  gtar, idfd, ids, ifds, indx, iwrite, neqn, np, npar, nparb, nparf, nprof, &
  numstp, ny, parnew, sens, title, yc )

!*****************************************************************************80
!
!! PR_SOLUTION prints out information about a single solution.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) neqn
  integer ( kind = 4 ) np
  integer ( kind = 4 ) npar
  integer ( kind = 4 ) nparb
  integer ( kind = 4 ) nparf
  integer ( kind = 4 ) ny

  real ( kind = 8 ) dpara3(npar)
  real ( kind = 8 ) dparfd(npar)
  real ( kind = 8 ) dparfdc(npar)
  real ( kind = 8 ) dparsn(npar)
  real ( kind = 8 ) g(neqn)
  real ( kind = 8 ) gdif(neqn,npar)
  real ( kind = 8 ) gdifc(neqn,npar)
  real ( kind = 8 ) gtar(neqn)
  integer ( kind = 4 ) idfd
  integer ( kind = 4 ) ids
  integer ( kind = 4 ) ifds
  integer ( kind = 4 ) indx(np,3)
  integer ( kind = 4 ) iwrite
  integer ( kind = 4 ) nprof(2*ny-1)
  integer ( kind = 4 ) numstp
  real ( kind = 8 ) parnew(npar)
  real ( kind = 8 ) sens(neqn,npar)
  character ( len = * ) title
  character ( len = 80 ) title_local
  real ( kind = 8 ) yc(np)
!
!  Print out stuff.
!
  if ( iwrite >= 1 ) then

    write ( *, * ) ' '
    write ( *, * ) 'PR_SOLUTION:'

    if ( len_trim ( title ) > 0 ) then
      write ( *, '(4x,a)' ) trim ( title )
      write ( *, * ) ' '
    end if

    call pr_parameter ( nparb, nparf, parnew )

  end if

  if ( iwrite >= 3 ) then

    title_local = 'Solution coefficients'
    call pr_puv ( g, indx, neqn, np, title_local )

  end if

  if ( iwrite >= 3 ) then
    if ( ids /= 0 ) then

      title_local = 'Discretized sensitivities'
      call pr_puv ( sens(1,2), indx, neqn, np, title_local )
      title_local = 'Chain rule on discretized sensitivities:'
      call pr_cost_sen ( dparsn, npar, title_local )

    end if
  end if

  if ( iwrite >= 3 ) then
    if ( ifds /= 0 ) then

      title_local = 'Finite coefficient differences'
      call pr_puv ( gdif(1,2), indx, neqn, np, title_local )
      title_local = 'Chain rule on finite difference sensitivities:'
      call pr_cost_sen ( dparfd, npar, title_local )

    end if
  end if

  if ( iwrite >= 3 ) then
    if ( ifds /= 0 ) then

      title_local = 'Adjusted finite coefficient differences'
      call pr_puv ( gdifc(1,2), indx, neqn, np, title_local )
      title_local = 'Chain rule on adjusted finite difference sens:'
      call pr_cost_sen ( dparfdc, npar, title_local )

    end if
  end if
!
!  Print finite cost gradient differences.
!
  if ( idfd /= 0 ) then
    if ( iwrite >= 3 ) then
      title_local = 'Finite cost gradient differences:'
      call pr_cost_sen ( dpara3, npar, title_local )
    end if
  end if

  if ( iwrite >= 3 ) then
    call pr_profile ( g, indx, neqn, np, nprof, ny, yc )
  end if

  if ( iwrite >= 3 ) then
    call pr_disc ( g, gtar, indx, neqn, np, nprof, ny, yc )
  end if

  return
end
subroutine pr_spl_data ( nspl, spl, tau )

!*****************************************************************************80
!
!! PR_SPL_DATA prints the raw spline data for the simple cases ISHAPE = 1 or 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) nspl

  integer ( kind = 4 ) i
  real ( kind = 8 ) spl(nspl)
  real ( kind = 8 ) tau(nspl)

  write ( *, * ) ' '
  write ( *, * ) 'PR_SPL_DATA:'
  write ( *, * ) '  Raw spline data.'
  write ( *, * ) ' '
  write ( *, * ) 'I, TAU(I), SPL(I)'

  do i = 1, nspl
    write ( *, * ) i, tau(i), spl(i)
  end do

  return
end
subroutine pr_spln ( jderiv, ipar, ishape, npts, npar, spl, tau, xhi, xlo )

!*****************************************************************************80
!
!! PR_SPLN prints a spline interpolant or its derivatives.
!
!  Discussion:
!
!    It can also print out the value or derivatives of one of the
!    "Lagrangian" spline functions, associated with any of NPAR
!    parameters.
!
!    The IPAR-th Lagrangian spline is 1 at the point associated with
!    the IPAR-th parameter, and 0 at the points associated with the
!    other parameters.  Using the notation S(IPAR)(X) for the value
!    of the IPAR-th Lagrangian spline at the point X, we can write:
!
!      S(X) = Sum (IPAR=1 to NPAR) C(IPAR) * S(IPAR)(X)
!
!    where S(X) is our resulting spline, which has the value C(IPAR)
!    at the IPAR-th point.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) npar

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ipar
  integer ( kind = 4 ) ishape
  integer ( kind = 4 ) jderiv
  integer ( kind = 4 ) npts
  real ( kind = 8 ) spl(4,npar+2,0:npar)
  real ( kind = 8 ) tau(npar+2)
  real ( kind = 8 ) xhi
  real ( kind = 8 ) xlo
  real ( kind = 8 ) xx
  real ( kind = 8 ) yvec(npar+2)
  real ( kind = 8 ) yy

  if ( npts<1 ) then
    write ( *, * ) ' '
    write ( *, * ) 'PR_SPLN - Serious error!'
    write ( *, * ) '  NPTS must be at least 1, but you entered ',npts
    return
  end if

  write ( *, * ) ' '
  write ( *, * ) 'PR_SPLN:'
  write ( *, * ) ' '

  yvec(1:npar+2) = spl(1,1:npar+2,0)

  do i = 1, npts

    if ( npts == 1 ) then
      xx = 0.5D+00 * (xlo+xhi)
    else
      xx = ( real ( npts - i,     kind = 8 ) * xlo &
           + real (        i - 1, kind = 8 ) * xhi ) &
           / real ( npts     - 1, kind = 8 )
    end if

    if ( ishape == 1 ) then

      if ( ipar == 0 ) then

        if ( jderiv == 0 ) then
          call plval ( npar+2, xx, tau, yy, yvec )
        else
          call pldx ( npar+2, xx, tau, yy, yvec )
        end if
      else
        call plval1 ( ipar+1, npar+2, xx, tau, yy )
      end if

    else if ( ishape == 2 ) then

      if ( ipar == 0 ) then
        if ( jderiv == 0 ) then
          call pqval ( npar+2, xx, tau, yy, yvec )
        else
          call pqdx ( npar+2, xx, tau, yy, yvec )
        end if
      else
        call pqval1 ( ipar+1, npar+2, xx, tau, yy )
      end if

    else if ( ishape == 3 ) then

      call ppvalu ( tau, spl(1,1,ipar), npar+1, 4, xx, jderiv, yy )

    end if

    write ( *, '(2g12.4)' ) xx, yy

  end do

  return
end
subroutine pr_work

!*****************************************************************************80
!
!! PR_WORK reports the amount of work carried out.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ival

  ival = 0

  write ( *, * ) ' '
  write ( *, * ) 'PR_WORK:'
  write ( *, * ) ' '
  call imemry('get','Points_Con',ival)
  write(*,*)ival,' continuation auxilliary points calculated.'

  call imemry('get','Points_FD',ival)
  write(*,*)ival,' finite difference auxilliary points calculated.'

  call imemry('get','Points_March',ival)
  write(*,*)ival,' marching points calculated.'

  call imemry('get','Points_Opt',ival)
  write(*,*)ival,' optimization points calculated.'

  call imemry('get','Points_Tar',ival)
  write(*,*)ival,' target points calculated.'

  write ( *, * ) ' '

  call imemry('get','CubSpl_calls',ival)
  write(*,*)ival,' calls to CubSpl.'

  call imemry('get','Factor_calls',ival)
  write(*,*)ival,' calls to SGBTRF, LU factor routine.'

  call imemry('get','Solve_calls',ival)
  write(*,*)ival,' calls to SGBTRS, the LU solver.'

  call imemry('get','FloSol_calls',ival)
  write(*,*)ival,' calls to FloSol.'

  call imemry('get','Fprime_calls',ival)
  write(*,*)ival,' calls to FPrime.'

  call imemry('get','Fx_calls',ival)
  write(*,*)ival,' calls to FX.'

  call imemry('get','GET_COST_calls',ival)
  write(*,*)ival,' calls to GET_COST.'

  call imemry('get','COST_GRADIENT_calls',ival)
  write(*,*)ival,' calls to COST_GRADIENT.'

  call imemry('get','GetFix_calls',ival)
  write(*,*)ival,' calls to GetFix.'

  call imemry('get','GetGrd_calls',ival)
  write(*,*)ival,' calls to GetGrd.'

  call imemry('get','GetSen_calls',ival)
  write(*,*)ival,' calls to GetSen.'

  call imemry('get','Newton_calls',ival)
  write(*,*)ival,' calls to Newton.'

  call imemry('get','PltWrt_calls',ival)
  write(*,*)ival,' calls to PltWrt.'

  call imemry('get','SetBas_calls',ival)
  write(*,*)ival,' calls to SetBas.'

  call imemry('get','NODE_SET_calls',ival)
  write(*,*)ival,' calls to NODE_SET.'

  call imemry('get','SetQXY_calls',ival)
  write(*,*)ival,' calls to SetQXY.'

  call imemry('get','XY_SET_calls',ival)
  write(*,*)ival,' calls to XY_SET.'

  call imemry('get','Snoit_calls',ival)
  write(*,*)ival,' calls to Snoit.'

  call imemry('get','SolCon_calls',ival)
  write(*,*)ival,' calls to SolCon.'

  call imemry('get','Sumit_calls',ival)
  write(*,*)ival,' calls to Sumit.'

  call imemry('get','Trans_calls',ival)
  write(*,*)ival,' calls to Trans.'

  call imemry('get','UVal_calls',ival)
  write(*,*)ival,' calls to UVal.'

  call imemry('get','UValQ_calls',ival)
  write(*,*)ival,' calls to UValQ.'

  call imemry('get','UpValQ_calls',ival)
  write(*,*)ival,' calls to UpValQ.'

  call imemry('get','Xofxsi_calls',ival)
  write(*,*)ival,' calls to Xofxsi.'

  write ( *, * ) ' '

  call imemry('get','SolCon_steps',ival)
  write(*,*)ival,' continuation steps by SolCon.'

  call imemry('get','Newton_damp',ival)
  write(*,*)ival,' attempted dampings of the Newton iteration.'

  call imemry('get','Newton_falter',ival)
  write(*,*)ival,' falterings of the Newton iteration.'

  call imemry('get','Newton_fail',ival)
  write(*,*)ival,' total failures of the Newton iteration.'

  call imemry('get','Newton_steps',ival)
  write(*,*)ival,' Newton iteration steps.'

  call imemry('get','Newton_zero',ival)
  write(*,*)ival,' zero-step Newton processes.'

  call imemry('get','Solve_sys',ival)
  write(*,*)ival,' right hand sides for LU solver.'

  call imemry('get','Restarts',ival)
  write(*,*)ival,' starts and restarts'

  return
end
subroutine probas(base,m,n)

!*****************************************************************************80
!
!! PROBAS orthonormalizes N vectors of length M.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) BASE(M,N).
!    On input, the columns of BASE contain N vectors, each of
!    length M, which span a space of dimension N.
!    On output, the columns of BASE form a basis for the same
!    space, that is, the columns are of unit Euclidean norm, and
!    orthogonal.
!
!    Input, integer ( kind = 4 ) M, the dimension of the higher order space.
!
!    Input, integer ( kind = 4 ) N, the dimension of the lower order space.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) base(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) temp
!
!  For each column J:
!
  do j = 1,n
!
!  ...consider a previous column, I, ...
    do i = 1,j-1
!
!  ...compute the projection of column J onto column I...
!
      temp = 0.0D+00
      do k = 1,m
        temp = temp+base(k,i)*base(k,j)
      end do
!
!  ...subtract off this projection...
!
      do k = 1,m
        base(k,j) = base(k,j)-temp*base(k,i)
      end do
    end do
!
!  Then compute the Euclidean norm of what's left of column J...
!
    temp = 0.0D+00
    do i = 1,m
      temp = temp+base(i,j)**2
    end do
    temp = sqrt(temp)
!
!  ...and normalize the column.
!
    if ( temp>0.0D+00 ) then
      do i = 1,m
        base(i,j) = base(i,j)/temp
      end do
    else
      write ( *, * ) 'PROBAS - Warning!'
      write ( *, * ) '  Column ',j,' of the basis is now 0.'
    end if

  end do

  return
end
subroutine projec(base,m,n,vecm,vecn)

!*****************************************************************************80
!
!! PROJEC projects an M vector into an N dimensional subspace.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) BASE(M,N).  The columns of BASE
!    contain N vectors, each of length M, which form the
!    orthonormal basis for a space of dimension N.
!
!    Input, integer ( kind = 4 ) M, the dimension of the higher order space.
!
!    Input, integer ( kind = 4 ) N, the dimension of the lower order space.
!
!    Input, real ( kind = 8 ) VECM(M), is an M dimensional vector.
!
!    Output, real ( kind = 8 ) VECN(N), the projection of VECM into
!    the lower dimensional space.  These values represent
!    coordinates in the lower order space.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) base(m,n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) vecm(m)
  real ( kind = 8 ) vecn(n)
!
!  The J-th coordinate of the projection of the vector
!  is simply the dot product of the vector with basis vector J.
!
  do j = 1,n
    vecn(j) = 0.0D+00
    do k = 1,m
      vecn(j) = vecn(j)+vecm(k)*base(k,j)
    end do
  end do

  return
end
subroutine qbf ( ielem,in,w,dwdx,dwdy,nelem,node,np,npe,xc,xq,yc,yq)

!*****************************************************************************80
!
!! QBF evaluates a quadratic basis function in a nonisoparametric element.
!
!  Diagram:
!
!      ^
!      |        2
!      |       /|
!   Y  |      4 5
!      |     /  |
!      |    1-6-3
!      |
!      +------------>
!             X
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IELEM, the number of the element we are
!    examining.  This will be a value between 1 and NELEM.
!
!    Input, integer ( kind = 4 ) IN, the number of the basis function we
!    want.  This will be a value between 1 and 6.  Functions
!    1 through 3 are associated with corners, 4 though 6
!    with sides.
!
!    Output, real ( kind = 8 ) W, DWDX, DWDY, the value of the
!    IN-th basisfunction and its X and Y derivatives, at the
!    given point.
!
!    Input, integer ( kind = 4 ) NELEM, the number of elements.
!
!    Input, integer ( kind = 4 ) NODE(MAXELM,6), contains the numbers
!    of the nodes that make up each element.  Element number
!    I is associated with nodes NODE(I,1) through NODE(I,6).
!
!    Input, integer ( kind = 4 ) NP, the number of nodes.
!
!    Input, real ( kind = 8 ) XC(NP), the X coordinates of the
!    nodes.
!
!    Input, real ( kind = 8 ) XQ, the X coordinate of the point
!    where the basis function is to be evaluated.
!
!    Input, real ( kind = 8 ) YC(NP), the Y coordinates of the nodes
!
!    Input, real ( kind = 8 ) YQ, the Y coordinate of the point wher
!    the basis function is to be evaluated.
!
  implicit none

  integer ( kind = 4 ) nelem
  integer ( kind = 4 ) np
  integer ( kind = 4 ) npe

  real ( kind = 8 ) c
  real ( kind = 8 ) d
  real ( kind = 8 ) dwdx
  real ( kind = 8 ) dwdy
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) ielem
  integer ( kind = 4 ) in
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) in3
  integer ( kind = 4 ) node(nelem,npe)
  real ( kind = 8 ) s
  real ( kind = 8 ) t
  real ( kind = 8 ) w
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) xq
  real ( kind = 8 ) yc(np)
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
    in2 = mod(in,3)+1
    in3 = mod(in+1,3)+1

    i1 = node(ielem,in1)
    i2 = node(ielem,in2)
    i3 = node(ielem,in3)

    d = (xc(i2)-xc(i1))*(yc(i3)-yc(i1))-(xc(i3)-xc(i1))*(yc(i2)-yc(i1))

    t = 1.0D+00+( (xq    -xc(i1))*(yc(i2)-yc(i3))+(xc(i3)-xc(i2))*(yq    -yc(i1)) )/d

    w = t*(2.0D+00*t-1.0D+00)

    dwdx = (yc(i2)-yc(i3)) * (4.0D+00*t-1.0D+00)/d
    dwdy = (xc(i3)-xc(i2)) * (4.0D+00*t-1.0D+00)/d
!
!  Case 2: We are inquiring about a basis function associated
!  with a midpoint.
!
  else if ( in >= 4 .and. in<=6 ) then

    in1 = in-3
    in2 = mod(in-3,3)+1
    in3 = mod(in-2,3)+1

    i1 = node(ielem,in1)
    i2 = node(ielem,in2)
    i3 = node(ielem,in3)

    d =  (xc(i2)-xc(i1))*(yc(i3)-yc(i1))-(xc(i3)-xc(i1))*(yc(i2)-yc(i1))

    c =  (xc(i3)-xc(i2))*(yc(i1)-yc(i2))-(xc(i1)-xc(i2))*(yc(i3)-yc(i2))

    t = 1.0D+00+( (xq    -xc(i1))*(yc(i2)-yc(i3))+(xc(i3)-xc(i2))*(yq    -yc(i1)) )/d

    s = 1.0D+00+( (xq    -xc(i2))*(yc(i3)-yc(i1))+(xc(i1)-xc(i3))*(yq    -yc(i2)) )/c

    w = 4.0D+00 * s*t
    dwdx = 4.0D+00 * ((yc(i3)-yc(i1))*t/c + (yc(i2)-yc(i3))*s/d)
    dwdy = 4.0D+00 * ((xc(i1)-xc(i3))*t/c + (xc(i3)-xc(i2))*s/d)

  else

    write ( *, * ) ' '
    write ( *, * ) 'QBF - Fatal error!'
    write ( *, * ) '  Request for basis function IN = ',in
    write ( *, * ) '  but IN must be between 1 and 6.'
    stop

  end if

  return
end
subroutine qsolve(a,area,costar,disjac,dopt,dpara3,dparfd,dparfdc,dparsn, &
  dpdyn,dudyn,dvdyn,dydpn,eqn,etan,etaq,g,gdif,gdifc,gold,gopt, &
  gradf,gtar,ibc,idfd,ids,ifds,igrad,igrid,igunit,ijac,indx,iopt,ipivot,iplot, &
  ipred,ishapb,ishapf,ismooth,isotri,itunit,itype,ivopt,iwrite,jjac,liv,lv, &
  maxnew,maxstp,nelem,neqn,nlband,node,nopt,np,npar,nparb,nparf,npe, &
  nprof,nrow, &
  numel,nx,ny,para,para1,parjac,parnew,partar,phi,res,sens,splbmp, &
  splflo,stpmax,syseqn,taubmp,tauflo,tolnew,tolopt,vopt,wateb,watep,wateu, &
  watev,wquad,xbl,xbord,xbr,xc,xopt,xprof,xquad,xsin,xsiq,ybl,ybord,ybr,yc, &
  yquad)

!*****************************************************************************80
!
!! QSOLVE seeks an optimal set of parameter values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) PARA1(NPAR), the initial estimate of the 
!    optimal parameters.
!
  implicit none

  integer ( kind = 4 ) liv
  integer ( kind = 4 ) lv
  integer ( kind = 4 ) nelem
  integer ( kind = 4 ) neqn
  integer ( kind = 4 ) np
  integer ( kind = 4 ) npar
  integer ( kind = 4 ) nparb
  integer ( kind = 4 ) nparf
  integer ( kind = 4 ) npe
  integer ( kind = 4 ) nrow
  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny

  real ( kind = 8 ) a(nrow,neqn)
  real ( kind = 8 ) area(nelem,3)
  real ( kind = 8 ) cgrad
  real ( kind = 8 ) cost
  real ( kind = 8 ) costb
  real ( kind = 8 ) costp
  real ( kind = 8 ) costu
  real ( kind = 8 ) costv
  real ( kind = 8 ) costar
  real ( kind = 8 ) disjac
  real ( kind = 8 ) dopt(npar)
  real ( kind = 8 ) dpara3(npar)
  real ( kind = 8 ) dparfd(npar)
  real ( kind = 8 ) dparfdc(npar)
  real ( kind = 8 ) dparsn(npar)
  real ( kind = 8 ) dpdyn(np)
  real ( kind = 8 ) dudyn(np)
  real ( kind = 8 ) dvdyn(np)
  real ( kind = 8 ) dydpn(np,nparb)
  character ( len = 2 ) eqn(neqn)
  real ( kind = 8 ) etan(6)
  real ( kind = 8 ) etaq(3)
  real ( kind = 8 ) flarea
  real ( kind = 8 ) g(neqn)
  real ( kind = 8 ) gdif(neqn,npar)
  real ( kind = 8 ) gdifc(neqn,npar)
  real ( kind = 8 ) gold(neqn)
  real ( kind = 8 ) gopt(npar)
  real ( kind = 8 ) gradf(neqn,npar)
  real ( kind = 8 ) gtar(neqn)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibc
  integer ( kind = 4 ) idfd
  integer ( kind = 4 ) ids
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ifds
  integer ( kind = 4 ) igrad
  integer ( kind = 4 ) igrid
  integer ( kind = 4 ) igunit
  integer ( kind = 4 ) ijac
  integer ( kind = 4 ) indx(np,3)
  integer ( kind = 4 ) iopt(npar)
  integer ( kind = 4 ) ipivot(neqn)
  integer ( kind = 4 ) iplot
  integer ( kind = 4 ) ipred
  integer ( kind = 4 ) ishapb
  integer ( kind = 4 ) ishapf
  integer ( kind = 4 ) ismooth
  integer ( kind = 4 ) isotri(nelem)
  integer ( kind = 4 ) itunit
  integer ( kind = 4 ) itype
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) ivopt(liv)
  integer ( kind = 4 ) iwrite
  integer ( kind = 4 ) jjac
  integer ( kind = 4 ) jopt
  integer ( kind = 4 ) maxnew
  integer ( kind = 4 ) maxstp
  integer ( kind = 4 ) nfail
  integer ( kind = 4 ) nlband
  integer ( kind = 4 ) node(nelem,npe)
  integer ( kind = 4 ) nopt
  integer ( kind = 4 ) nprof(2*ny-1)
  integer ( kind = 4 ) numel(np)
  integer ( kind = 4 ) numstp
  real ( kind = 8 ) para(npar)
  real ( kind = 8 ) para1(npar)
  real ( kind = 8 ) parjac(npar)
  real ( kind = 8 ) parnew(npar)
  real ( kind = 8 ) parold(npar)
  real ( kind = 8 ) partar(npar)
  real ( kind = 8 ) phi(nelem,3,6,10)
  real ( kind = 8 ) res(neqn)
  real ( kind = 8 ) sens(neqn,npar)
  real ( kind = 8 ) splbmp(4,nparb+2,0:nparb)
  real ( kind = 8 ) splflo(4,nparf+2,0:nparf)
  real ( kind = 8 ) stpmax
  character ( len = 20 ) syseqn
  real ( kind = 8 ) taubmp(nparb+2)
  real ( kind = 8 ) tauflo(nparf+2)
  character ( len = 80 ) title
  real ( kind = 8 ) tolnew
  real ( kind = 8 ) tolopt
  real ( kind = 8 ) vopt(lv)
  real ( kind = 8 ) wateb
  real ( kind = 8 ) watep
  real ( kind = 8 ) wateu
  real ( kind = 8 ) watev
  real ( kind = 8 ) wquad(3)
  real ( kind = 8 ) xbl
  real ( kind = 8 ) xbord(nx)
  real ( kind = 8 ) xbr
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) xopt(npar)
  real ( kind = 8 ) xprof
  real ( kind = 8 ) xquad(nelem,3)
  real ( kind = 8 ) xsin(6)
  real ( kind = 8 ) xsiq(3)
  real ( kind = 8 ) ybl
  real ( kind = 8 ) ybord(ny)
  real ( kind = 8 ) ybr
  real ( kind = 8 ) yc(np)
  real ( kind = 8 ) yquad(nelem,3)
!
!  Initialize some things.
!
  nfail = 0
  numstp = 0
!
!  Set default values of the optimizer variables.
!
  ival = 2
  call deflt ( ival, ivopt, liv, lv, vopt )
!
!  Adjust the values of some optimizer variables.
!
  if ( tolopt > 0.0D+00 ) then
    vopt(31) = tolopt
    vopt(32) = tolopt
    vopt(33) = tolopt
    vopt(34) = tolopt
    vopt(37) = tolopt
  end if

  ivopt(1) = 12
  ivopt(19) = 0

  dpara3(1:npar) = 0.0D+00
  dparfd(1:npar) = 0.0D+00
  dparfdc(1:npar) = 0.0D+00
  dparsn(1:npar) = 0.0D+00

  g(1:neqn) = 0.0D+00
  parnew(1:npar) = para1(1:npar)

  gold(1:neqn) = 0.0D+00
  parold(1:npar-1) = 0.0D+00
  parold(nparf+nparb+1) = para1(nparf+nparb+1)
!
!  Optimization loop
!
10    continue
!
!  Call the optimizer to get a new set of parameter values, PARA.
!
  write ( *, * ) 'QSOLVE: Top of optimization loop.'

  cgrad = 0.0D+00
  jopt = 0

  do i = 1,npar

    if ( iopt(i) == 1 ) then

      jopt = jopt+1
      xopt(jopt) = parnew(i)

      if ( igrad == 1 ) then
        gopt(jopt) = dparsn(i)
      else if ( igrad == 2 ) then
        gopt(jopt) = dparfd(i)
      else if ( igrad == 3 ) then
        gopt(jopt) = dparfdc(i)
      else if ( igrad == 4 ) then
        gopt(jopt) = dpara3(i)
      end if

      cgrad = max ( cgrad, abs ( gopt(jopt) ) )

    end if

  end do

11    continue

  call sumit ( dopt, cost, gopt, ivopt, liv, lv, nopt, vopt, xopt )

  ival = 1
  call imemry ( 'inc', 'Sumit_calls', ival )

  jopt = 0
  do i = 1,npar
    if ( iopt(i) == 1 ) then
      jopt = jopt+1
      parnew(i) = xopt(jopt)
      para(i) = xopt(jopt)
    end if
  end do
!
!  Did the optimizer converge, so that we're done?
!
  if ( ivopt(1) >= 3 .and. ivopt(1)<=8 ) then

    write ( *, * ) ' '
    write ( *, * ) 'QSOLVE - Convergence to a minimizer was achieved!'
    write ( *, * ) '  IVOPT(1) = ',ivopt(1)

    call chkopt(cost,costar,dparfd,dparfdc,dparsn,g,gtar,ids,ifds,neqn,npar, &
      parnew,partar)

    go to 20
!
!  ...or did the optimizer find a serious problem, so that we're "finished"?
!
  else if ( ivopt(1)>8 ) then

    write ( *, * ) ' '
    write ( *, * ) 'QSOLVE - Warning!'
    write ( *, * ) '  IVOPT(1) = ', ivopt(1)
    go to 20
!
!  The optimizer has not converged, and it hasn't had a failure.
!
!  The optimizer has returned a new, rough, guess for the
!  minimizer, and wants us to evaluate the cost function there.
!
!  To do this, we must find the flow solution (U,V,P) corresponding
!  to the given parameters.
!
  else if ( ivopt(1) == 1 ) then

    numstp = numstp + 1
!
!  Now we are ready to request that SOLCON compute the solution with
!  parameters PARNEW, given the solution GOLD at PAROLD.
!
    write ( *, * ) 'QSOLVE calls FLOSOL'

    call flosol ( a,area,disjac,eqn,etan,etaq,flarea,g,ierror,igrid,ijac, &
      indx,ipivot,ishapb,ishapf,isotri,iwrite,jjac,maxnew,nelem,neqn, &
      nlband,node,np,npar,nparb,nparf,npe,nrow,nx,ny,para,parjac,phi,res, &
      splbmp,splflo,syseqn,taubmp,tauflo,tolnew,wquad,xbl,xbord,xbr,xc, &
      xquad,xsin,xsiq,ybl,ybord,ybr,yc,yquad)

!   call solcon(a,area,disjac,dpara3,dpdyn,dudyn,dvdyn,dydpn,eqn,etan, &
!     etaq,flarea,g,gdif,gdifc,gold,gradf,gtar,ibc,idfd,ids,ierror, &
!     ifds,igrid,ijac,indx,iopt,ipivot,ipred,ishapb,ishapf,ismooth,isotri, &
!     itype,iwrite,jjac,maxnew,nelem,neqn,nlband,node,np,npar,nparb,nparf, &
!     npe,nprof,nrow,numel,nx,ny,para,parjac,parnew,parold,phi,res,sens,splbmp, &
!     splflo,stpmax,syseqn,taubmp,tauflo,tolnew,wateb,watep,wateu,watev, &
!     wquad,xbl,xbord,xbr,xc,xquad,xsin,xsiq,ybl,ybord,ybr,yc,yquad)

    write ( *, * ) 'QSOLVE has returned from FLOSOL'

    if ( ierror /= 0 ) then
      nfail = nfail+1
      if ( nfail <= 1 ) then
        write ( *, * ) ' '
        write ( *, * ) 'QSOLVE - Warning!'
        write ( *, * ) '  SOLCON returns IERROR = ',ierror
        write ( *, * ) '  Requesting that SUMIT try a smaller step.'
        ivopt(2) = 1
        go to 11
      else
        write ( *, * ) ' '
        write ( *, * ) 'QSOLVE - Fatal error!'
        write ( *, * ) '  SOLCON returns IERROR = ',ierror
        return
      end if
    else
      nfail = 0
    end if

    parold(1:npar) = parnew(1:npar)
    gold(1:neqn) = g(1:neqn)
!
!  Compute the cost COST associated with the solution G,
!  determined by the parameters PARNEW.
!
    call get_cost(cost,costb,costp,costu,costv,g,gtar,indx,ishapb,neqn,np, &
      nparb,nprof,ny,splbmp,taubmp,wateb,watep,wateu,watev,xbl,xbr,ybl,ybr,yc)

    if ( iwrite == 1 ) then
      title = 'Cost'
      call pr_cost1(cost,title)
    else if ( iwrite >= 2 ) then
      title = 'Cost'
      call pr_cost2(cost,costb,costp,costu,costv,title,wateb,watep,wateu,watev)
    end if
!
!  Print stuff out.
!
    if ( ifds /= 0 ) then
      call cost_gradient(dparfd,g,gtar,indx,ishapb,neqn,np,npar, &
        nparb,nparf,nprof, &
        ny,gdif,splbmp,taubmp,wateb,watep,wateu,watev,xbl,xbr,ybl,ybr,yc)
    end if

    if ( ifds /= 0 ) then
      call cost_gradient(dparfdc,g,gtar,indx,ishapb,neqn,np,npar,nparb,nparf,nprof, &
        ny,gdifc,splbmp,taubmp,wateb,watep,wateu,watev,xbl,xbr,ybl,ybr,yc)
    end if

    if ( ids /= 0 ) then
      call cost_gradient(dparsn,g,gtar,indx,ishapb,neqn,np,npar,nparb,nparf,nprof, &
      ny,sens,splbmp,taubmp,wateb,watep,wateu,watev,xbl,xbr,ybl,ybr,yc)
    end if

    title = 'Optimization solution'
    call pr_solution ( dpara3,dparfd,dparfdc,dparsn,g,gdif,gdifc,gtar,idfd,ids,ifds, &
      indx,iwrite,neqn,np,npar,nparb,nparf,nprof,numstp,ny,parnew,sens,title,yc)

    if ( iwrite >= 3 ) then
      call pr_gs2 ( flarea,gdif,gdifc,gradf,indx,iopt,neqn,np,npar,ny,sens,xc,yc)
    end if

    if ( itunit /= 0 ) then
      write(itunit,'(1x,6g14.6)')cost,(parnew(i),i = 1,npar)
      write(itunit,'(1x,5g14.6)')(dparsn(i),i = 1,npar)
      write(itunit,'(1x,5g14.6)')(dparfd(i),i = 1,npar)
      write(itunit,'(1x,5g14.6)')(dpara3(i),i = 1,npar)
    end if

    if ( iplot > 0 .and. mod ( numstp, iplot ) == 0 ) then

      call plot_file_write ( eqn, g, gdif, igunit, indx, isotri, iwrite, nelem, &
        neqn, node, np, npar, npe, nprof, nx, ny, para, sens, xc, xprof, yc )

    end if

    if ( numstp >= maxstp ) then

      if ( iwrite < 2 ) then
        write ( *, * ) ' '
        write ( *, * ) 'QSOLVE - Point ',numstp
        call pr_parameter(nparb,nparf,parnew)

        title = 'Cost'
        call pr_cost1(cost,title)
      end if

      write ( *, * ) ' '
      write ( *, * ) 'QSOLVE - Warning!'
      write ( *, * ) '  Number of steps exceeded.'
      write ( *, * ) ' '
      write ( *, * ) 'Last gradient vector was:'
      write ( *, * ) ' '
      do i = 1, npar
        if ( igrad == 1 ) then
          write(*,*)dparsn(i)
        else if ( igrad == 2 ) then
          write(*,*)dparfd(i)
        else if ( igrad == 3 ) then
          write(*,*)dparfdc(i)
        else if ( igrad == 4 ) then
          write(*,*)dpara3(i)
        end if
      end do

      go to 20
    end if

    go to 10
!
!  Did the optimizer ask us to return the derivative values
!  associated with the current set of parameters?
!
  else if ( ivopt(1) == 2 ) then

    go to 10
!
!  We can't figure out what SUMIT is trying to tell us.
!
  else
    write ( *, * ) ' '
    write ( *, * ) 'QSOLVE - Warning!'
    write ( *, * ) '  Unknown value of IVOPT(1) = ',ivopt(1)
    go to 20
  end if

20    continue
!
!  If IPLOT < 0, we want to write out graphics information
!  for the last point computed.
!
  if ( iplot < 0 ) then
    call plot_file_write ( eqn, g, gdif, igunit, indx, isotri, iwrite, nelem, &
      neqn, node, np, npar, npe, nprof, nx, ny, para, sens, xc, xprof, yc )
  end if

  return
end
subroutine refbsp ( q, dqdx, dqdy, detadx, detady, iq, dxsidx, dxsidy, eta, &
  xsi )

!*****************************************************************************80
!
!! REFBSP evaluates a linear basis function in the reference triangle.
!
!  Discussion:
!
!    The routine evaluates one of the three linear basis functions,
!    and its X and Y derivatives, at a particular point (X,Y)
!    in a particular element, by referring to the corresponding
!    points (XSI,ETA) in the reference triangle.
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
!    Here is a graph of the (XSI, ETA) reference triangle we will use.
!
!        ^
!        |
!      1 +        2
!        |       /|
!  ETA   |      / |
!        |     /  |
!      0 +    1---3
!        |
!        +----+---+--->
!             0   1
!
!              XSI
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 January 2001
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
!    Input, integer ( kind = 4 ) IQ, the local node number, between 1 and
!    3, whose basis function is being evaluated.
!
!    Input, real ( kind = 8 ) DXSIDX, DXSIDY, the partial derivative
!    d XSI/d X and d XSI/d Y at (ETA,XSI).
!
!    Input, real ( kind = 8 ) ETA, XSI, the local coordinates of the
!    at which the basis information is desired.
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
  integer ( kind = 4 ) iq
  real ( kind = 8 ) q
  real ( kind = 8 ) xsi

  if ( iq == 1 ) then
    q = 1.0D+00-xsi
    dqdxsi = -1.0D+00
    dqdeta = 0.0D+00
  else if ( iq == 2 ) then
    q = eta
    dqdxsi = 0.0D+00
    dqdeta = 1.0D+00
  else if ( iq == 3 ) then
    q = xsi-eta
    dqdxsi = 1.0D+00
    dqdeta = -1.0D+00
  else if ( iq >= 4 .and. iq<=6 ) then
    q = 0.0D+00
    dqdxsi = 0.0D+00
    dqdeta = 0.0D+00
  else
    write ( *, * ) ' '
    write ( *, * ) 'REFBSP - Fatal error!'
    write ( *, * ) '  Request for basis function IQ = ',iq
    write ( *, * ) '  but IQ must be between 1 and 6.'
    stop
  end if

  dqdx = dqdxsi * dxsidx + dqdeta * detadx
  dqdy = dqdxsi * dxsidy + dqdeta * detady

  return
end
subroutine refqbf ( w, dwdx, dwdy, detadx, detady, dxsidx, dxsidy, eta, &
  iq, xsi )

!*****************************************************************************80
!
!! REFQBF evaluates one of the six quadratic basis functions.
!
!  Discussion:
!
!    This routine evaluates one of the six quadratic basis functions,
!    and its X and Y derivatives, at a particular point in a
!    particular element, by referring to the reference triangle.
!
!    The point we are interested in is referred to by its coordinates
!    in the reference triangle.  That is, we are given coordinates
!    (XSI, ETA), even though, physically, we are interested
!    in points in (X, Y) space.
!
!    It is assumed that we already know the value of the jacobian
!    of the isoparametric transformation between the (XSI, ETA) and
!    (X, Y) spaces.  The four entries of the jacobian are
!    symbolically named DETADX, DETADY, DXSIDX and DXSIDY, and
!    we know that the jacobian gives us the following relation
!    between derivatives with respect to XSI and ETA, and derivatives
!    with respect to X and Y:
!
!      d F(X,Y)/dX     (d XSI/dX  d ETA/dX )   ( d F(XSI, ETA)/d XSI )
!      d F(X,Y)/dY  =  (d XSI/dY  d ETA/dY ) * ( d F(XSI, ETA)/d ETA )
!
!    Here is a graph of the (XSI, ETA) reference triangle we will
!    use.
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
!                XSI
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 January 2001
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
!    Input, integer ( kind = 4 ) IQ, the local node number, between 1 and
!    6, whose basis function is being evaluated.
!
!    Input, real ( kind = 8 ) XSI, the XSI coordinate of the point.
!
  implicit none

  real ( kind = 8 ) detadx
  real ( kind = 8 ) detady
  real ( kind = 8 ) dwdeta
  real ( kind = 8 ) dwdx
  real ( kind = 8 ) dwdxsi
  real ( kind = 8 ) dwdy
  real ( kind = 8 ) dxsidx
  real ( kind = 8 ) dxsidy
  real ( kind = 8 ) eta
  integer ( kind = 4 ) iq
  real ( kind = 8 ) w
  real ( kind = 8 ) xsi
!
!  Evaluate W, the quadratic basis function.
!  Evaluate DWDXSI and DWDETA, the partial derivatives d W/d XSI
!  and d W/d ETA.
!
!  Basis 1 is zero if XSI = 0.5 or XSI=1.
!
  if ( iq == 1 ) then
    w =  (2.0D+00*xsi-1.0D+00) * (xsi-1.0D+00)
    dwdxsi = -3.0D+00+4.0D+00*xsi
    dwdeta = 0.0D+00
!
!  Basis 2 is zero if ETA = 0 or ETA=0.5.
!
  else if ( iq == 2 ) then
    w =  eta * (2.0D+00*eta-1.0D+00)
    dwdxsi = 0.0D+00
    dwdeta = -1.0D+00+4.0D+00*eta
!
!  Basis 3 is zero if XSI = ETA, or XSI=ETA+0.5
!
  else if ( iq == 3 ) then
    w =  (xsi-eta) * (2.0D+00*xsi-2.0D+00*eta-1.0D+00)
    dwdxsi = -1.0D+00+4.0D+00*xsi-4.0D+00*eta
    dwdeta = 1.0D+00-4.0D+00*xsi+4.0D+00*eta
!
!  Basis 4 is zero if ETA = 0 or XSI=1.
!
  else if ( iq == 4 ) then
    w =  4.0D+00 * eta * (1.0D+00-xsi)
    dwdxsi = -4.0D+00*eta
    dwdeta = 4.0D+00-4.0D+00*xsi
!
!  Basis 5 is zero if ETA = 0 or XSI=ETA.
!
  else if ( iq == 5 ) then
    w = 4.0D+00 * eta * (xsi-eta)
    dwdxsi = 4.0D+00*eta
    dwdeta = 4.0D+00*xsi-8.0D+00*eta
!
!  Basis 6 is zero if XSI = ETA or XSI=1.
!
  else if ( iq == 6 ) then
    w = 4.0D+00 * (xsi-eta) * (1.0D+00-xsi)
    dwdxsi = 4.0D+00-8.0D+00*xsi+4.0D+00*eta
    dwdeta = -4.0D+00+4.0D+00*xsi
!
!  Stop if we were given an unexpected value of IQ.
!
  else
    write ( *, * ) ' '
    write ( *, * ) 'REFQBF - Fatal error!'
    write ( *, * ) '  A basis function index must be between 1 and 6,'
    write ( *, * ) '  but you input the value IQ = ',iq
    stop
  end if
!
!  Convert the d W/d XSI and d W/d ETA derivatives to d W/d X
!  and d W/d Y.
!
  dwdx = dwdxsi*dxsidx + dwdeta*detadx
  dwdy = dwdxsi*dxsidy + dwdeta*detady

  return
end
subroutine rint_to_rint ( rmin, rmax, r, r2min, r2max, r2 )

!*****************************************************************************80
!
!! RINT_TO_RINT maps a real interval to another real interval.
!
!  Formula:
!
!    R2 : =  R2MIN + ( R2MAX - R2MIN ) * ( R - RMIN ) / ( RMAX - RMIN )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) RMIN, RMAX, the first real range.
!
!    Input, real ( kind = 8 ) R, the real number to be converted.
!
!    Input, real ( kind = 8 ) R2MAX, R2MIN, the second real range.
!
!    Output, real ( kind = 8 ) R2, the corresponding value in the range [R2MIN,R2MAX].
!
  implicit none

  real ( kind = 8 ) r
  real ( kind = 8 ) rmax
  real ( kind = 8 ) rmin
  real ( kind = 8 ) r2
  real ( kind = 8 ) r2max
  real ( kind = 8 ) r2min

  if ( rmax == rmin ) then

    r2 = ( r2max + r2min ) / 2.0D+00

  else

    r2 = ( ( ( rmax - r ) * r2min + ( r - rmin ) * r2max ) / ( rmax - rmin ) )

  end if

  return
end
subroutine rsolve(a,area,disjac,dpara3,dparfd,dparfdc,dparsn,dpdyn,dudyn, &
  dvdyn,dydpn,eqn,etan,etaq,g,gdif,gdifc,gradf,gtar,ibc, &
  idfd,ids,ifds,igrid,igunit,ijac,indx,iopt,ipivot,iplot,ipred,ishapb,ishapf, &
  ismooth,isotri,itype,iwrite,jjac,maxnew,nelem,neqn,nlband,node,np,npar, &
  nparb,nparf,npe,nprof,nrow,numel,nx,ny,para,para1,parjac,phi,res, &
  sens,splbmp,splflo,stpmax,syseqn,taubmp,tauflo,tolnew,wateb,watep,wateu, &
  watev,wquad,xbl,xbord,xbr,xc,xprof,xquad,xsin,xsiq,ybl,ybord,ybr,yc,yquad)

!*****************************************************************************80
!
!! RSOLVE computes the flow solution for a given set of parameters.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) PARA(NPAR), the parameters associated with the problem.
!
  implicit none

  integer ( kind = 4 ) nelem
  integer ( kind = 4 ) neqn
  integer ( kind = 4 ) np
  integer ( kind = 4 ) npar
  integer ( kind = 4 ) nparb
  integer ( kind = 4 ) nparf
  integer ( kind = 4 ) npe
  integer ( kind = 4 ) nrow
  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny

  real ( kind = 8 ) a(nrow,neqn)
  real ( kind = 8 ) area(nelem,3)
  real ( kind = 8 ) disjac
  real ( kind = 8 ) dpara3(npar)
  real ( kind = 8 ) dparfd(npar)
  real ( kind = 8 ) dparfdc(npar)
  real ( kind = 8 ) dparsn(npar)
  real ( kind = 8 ) dpdyn(np)
  real ( kind = 8 ) dudyn(np)
  real ( kind = 8 ) dvdyn(np)
  real ( kind = 8 ) dydpn(np,nparb)
  character ( len = 2 ) eqn(neqn)
  real ( kind = 8 ) etan(6)
  real ( kind = 8 ) etaq(3)
  real ( kind = 8 ) flarea
  real ( kind = 8 ) g(neqn)
  real ( kind = 8 ) gdif(neqn,npar)
  real ( kind = 8 ) gdifc(neqn,npar)
  real ( kind = 8 ) gradf(neqn,npar)
  real ( kind = 8 ) gtar(neqn)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibc
  integer ( kind = 4 ) idfd
  integer ( kind = 4 ) ids
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ifds
  integer ( kind = 4 ) igrid
  integer ( kind = 4 ) igunit
  integer ( kind = 4 ) ijac
  integer ( kind = 4 ) indx(np,3)
  integer ( kind = 4 ) iopt(npar)
  integer ( kind = 4 ) ipivot(neqn)
  integer ( kind = 4 ) iplot
  integer ( kind = 4 ) ipred
  integer ( kind = 4 ) ishapb
  integer ( kind = 4 ) ishapf
  integer ( kind = 4 ) ismooth
  integer ( kind = 4 ) isotri(nelem)
  integer ( kind = 4 ) itype
  integer ( kind = 4 ) iwrite
  integer ( kind = 4 ) jjac
  integer ( kind = 4 ) maxnew
  integer ( kind = 4 ) nlband
  integer ( kind = 4 ) node(nelem,npe)
  integer ( kind = 4 ) nprof(2*ny-1)
  integer ( kind = 4 ) numel(np)
  integer ( kind = 4 ) numstp
  real ( kind = 8 ) para(npar)
  real ( kind = 8 ) para1(npar)
  real ( kind = 8 ) parjac(npar)
  real ( kind = 8 ) phi(nelem,3,6,10)
  real ( kind = 8 ) res(neqn)
  real ( kind = 8 ) sens(neqn,npar)
  real ( kind = 8 ) splbmp(4,nparb+2,0:nparb)
  real ( kind = 8 ) splflo(4,nparf+2,0:nparf)
  real ( kind = 8 ) stpmax
  character ( len = 20 ) syseqn
  real ( kind = 8 ) taubmp(nparb+2)
  real ( kind = 8 ) tauflo(nparf+2)
  character ( len = 80 ) title
  real ( kind = 8 ) tolnew
  real ( kind = 8 ) wateb
  real ( kind = 8 ) watep
  real ( kind = 8 ) wateu
  real ( kind = 8 ) watev
  real ( kind = 8 ) wquad(3)
  real ( kind = 8 ) xbl
  real ( kind = 8 ) xbord(nx)
  real ( kind = 8 ) xbr
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) xprof
  real ( kind = 8 ) xquad(nelem,3)
  real ( kind = 8 ) xsin(6)
  real ( kind = 8 ) xsiq(3)
  real ( kind = 8 ) ybl
  real ( kind = 8 ) ybord(ny)
  real ( kind = 8 ) ybr
  real ( kind = 8 ) yc(np)
  real ( kind = 8 ) yquad(nelem,3)
!
!  Initialize some things.
!
  g(1:neqn) = 0.0D+00
  dpara3(1:npar) = 0.0D+00
  dparfd(1:npar) = 0.0D+00
  dparfdc(1:npar) = 0.0D+00
  dparsn(1:npar) = 0.0D+00
!
!  I don't know where the parameters are!
!
  para(1:npar) = para1(1:npar)
!
!  Shouldn't we call FLOSOL, rather than SOLCON?
!
  call flosol ( a,area,disjac,eqn,etan,etaq,flarea,g,ierror,igrid,ijac, &
    indx,ipivot,ishapb,ishapf,isotri,iwrite,jjac,maxnew,nelem,neqn, &
    nlband,node,np,npar,nparb,nparf,npe,nrow,nx,ny,para,parjac,phi,res, &
    splbmp,splflo,syseqn,taubmp,tauflo,tolnew,wquad,xbl,xbord,xbr,xc, &
    xquad,xsin,xsiq,ybl,ybord,ybr,yc,yquad)

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'RSOLVE - Fatal error!'
    write ( *, * ) '  FLOSOL returns IERROR = ', ierror
    return
  end if
!
!  Print stuff out.
!
  numstp = 1

  if ( ifds /= 0 ) then
    call cost_gradient(dparfd,g,gtar,indx,ishapb,neqn,np,npar,nparb,nparf,nprof,ny, &
      gdif,splbmp,taubmp,wateb,watep,wateu,watev,xbl,xbr,ybl,ybr,yc)
  end if

  if ( ifds /= 0 ) then
    call cost_gradient(dparfdc,g,gtar,indx,ishapb,neqn,np,npar,nparb,nparf,nprof,ny, &
      gdifc,splbmp,taubmp,wateb,watep,wateu,watev,xbl,xbr,ybl,ybr,yc)
  end if

  if ( ids /= 0 ) then
    call cost_gradient(dparsn,g,gtar,indx,ishapb,neqn,np,npar, &
      nparb,nparf,nprof,ny, &
      sens,splbmp,taubmp,wateb,watep,wateu,watev,xbl,xbr,ybl,ybr,yc)
  end if

  title = 'Solution computed by RSOLVE:'

  call pr_solution ( dpara3, dparfd, dparfdc, dparsn, g, gdif, gdifc, gtar, &
    idfd, ids, ifds, indx, iwrite, neqn, np, npar, nparb, nparf, nprof, numstp, &
    ny, para, sens, title, yc )

  if ( iwrite >= 3 ) then
    call pr_gs2(flarea,gdif,gdifc,gradf,indx,iopt,neqn,np,npar,ny,sens,xc,yc)
  end if

  if ( iplot == 1 .or. iplot == -1 ) then
    call plot_file_write ( eqn, g, gdif, igunit, indx, isotri, iwrite, nelem, &
      neqn, node, np, npar, npe, nprof, nx, ny, para, sens, xc, xprof, yc )
  end if

  return
end
subroutine r4vec_even ( alo, ahi, n, a )

!*****************************************************************************80
!
!! R4VEC_EVEN returns N real values, evenly spaced between ALO and AHI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ALO, AHI, the low and high values.
!
!    Input, integer ( kind = 4 ) N, the number of values.
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
      a(i) = ( real ( n - i ) * alo + real ( i - 1 ) * ahi ) / real ( n - 1 )
    end do

  end if

  return
end
subroutine s_before_ss_copy ( s, ss, s2 )

!*****************************************************************************80
!
!! S_BEFORE_SS_COPY copies a string up to a given substring.
!
!  Discussion:
!
!    S and S2 can be the same object, in which case the string is
!    overwritten by a copy of itself up to the substring, followed
!    by blanks.
!
!  Example:
!
!    Input:
!
!      S = 'ABCDEFGH'
!      SS = 'EF'
!
!    Output:
!
!      S2 = 'ABCD'.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string to be copied.
!
!    Input, character ( len = * ) SS, the substring before which the copy stops.
!
!    Output, character ( len = * ) S2, the copied portion of S.
!
  implicit none

  integer ( kind = 4 ) last
  integer ( kind = 4 ) last_s2
  character ( len = * ) s
  character ( len = * ) s2
  character ( len = * ) ss
!
!  Find the first occurrence of the substring.
!
  last = index ( s, ss )
!
!  If the substring doesn't occur at all, behave as though it begins
!  just after the string terminates.
!
!  Now redefine LAST to point to the last character to copy before
!  the substring begins.
!
  if ( last == 0 ) then
    last = len ( s )
  else
    last = last - 1
  end if
!
!  Now adjust again in case the copy holder is "short".
!
  last_s2 = len ( s2 )

  last = min ( last, last_s2 )
!
!  Copy the beginning of the string.
!  Presumably, compilers now understand that if LAST is 0, we don't
!  copy anything.
!  Clear out the rest of the copy.
!
  s2(1:last) = s(1:last)
  s2(last+1:last_s2) = ' '

  return
end
subroutine s_blank_delete ( s )

!*****************************************************************************80
!
!! S_BLANK_DELETE removes blanks from a string, left justifying the remainder.
!
!  Discussion:
!
!    All TAB characters are also removed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) S, the string to be transformed.
!
  implicit none

  character c
  integer ( kind = 4 ) iget
  integer ( kind = 4 ) iput
  integer ( kind = 4 ) nchar
  character ( len = * ) s
  character, parameter :: TAB = char ( 9 )

  iput = 0
  nchar = len_trim ( s )

  do iget = 1, nchar

    c = s(iget:iget)

    if ( c /= ' ' .and. c /= TAB ) then
      iput = iput + 1
      s(iput:iput) = c
    end if

  end do

  s(iput+1:nchar) = ' '

  return
end
subroutine s_cap ( s )

!*****************************************************************************80
!
!! S_CAP replaces any lowercase letters by uppercase ones in a string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) S, the string to be transformed.
!
  implicit none

  character c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) nchar
  character ( len = * ) s

  nchar = len_trim ( s )

  do i = 1, nchar

    c = s(i:i)
    call ch_cap ( c )
    s(i:i) = c

  end do

  return
end
function s_eqi ( s1, s2 )

!*****************************************************************************80
!
!! S_EQI is a case insensitive comparison of two strings for equality.
!
!  Example:
!
!    S_EQI ( 'Anjana', 'ANJANA' ) is .TRUE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S1, S2, the strings to compare.
!
!    Output, logical S_EQI, the result of the comparison.
!
  implicit none

  character c1
  character c2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) len1
  integer ( kind = 4 ) len2
  integer ( kind = 4 ) lenc
  logical s_eqi
  character ( len = * ) s1
  character ( len = * ) s2

  len1 = len ( s1 )
  len2 = len ( s2 )
  lenc = min ( len1, len2 )

  s_eqi = .false.

  do i = 1, lenc

    c1 = s1(i:i)
    c2 = s2(i:i)
    call ch_cap ( c1 )
    call ch_cap ( c2 )

    if ( c1 /=  c2 ) then
      return
    end if

  end do

  do i = lenc + 1, len1
    if ( s1(i:i) /=  ' ' ) then
      return
    end if
  end do

  do i = lenc + 1, len2
    if ( s2(i:i) /=  ' ' ) then
      return
    end if
  end do

  s_eqi = .true.

  return
end
subroutine scopy(n,dx,incx,dy,incy)

!*****************************************************************************80
!
!! SCOPY copies a vector, x, to a vector, y.
!
  implicit none

  real ( kind = 8 ) dx(*)
  real ( kind = 8 ) dy(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) incy
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) iy
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  if ( n <= 0)return

  if ( incx == 1 .and. incy==1)go to 20

  ix = 1
  iy = 1
  if ( incx<0)ix = (-n+1)*incx + 1
  if ( incy<0)iy = (-n+1)*incy + 1

  do i = 1,n
    dy(iy) = dx(ix)
    ix = ix + incx
    iy = iy + incy
  end do

  return
!
   20 m = mod(n,7)

  do i = 1,m
    dy(i) = dx(i)
  end do

  do i = m+1,n,7
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
function sdot(n,dx,incx,dy,incy)

!*****************************************************************************80
!
!! SDOT forms the dot product of two vectors.
!
  implicit none

  real ( kind = 8 ) sdot
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

  sdot = 0.0D+00
  dtemp = 0.0D+00

  if ( n <= 0)return

  if ( incx == 1 .and. incy==1)go to 20

  if ( incx<0 ) then
    ix = (-n+1)*incx + 1
  else
    ix = 1
  end if

  if ( incy<0 ) then
    iy = (-n+1)*incy + 1
  else
    iy = 1
  end if

  do i = 1,n
    dtemp = dtemp + dx(ix)*dy(iy)
    ix = ix + incx
    iy = iy + incy
  end do

  sdot = dtemp
  return
!
   20 m = mod(n,5)

  do i = 1,m
    dtemp = dtemp + dx(i)*dy(i)
  end do

  do i = m+1,n,5
    dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) + &
      dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
  end do

  sdot = dtemp

  return
end
subroutine setban ( indx,maxrow,nelem,neqn,nlband,node,np,npe,nrow)

!*****************************************************************************80
!
!! SETBAN computes the half band width of the Jacobian matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) nelem
  integer ( kind = 4 ) np
  integer ( kind = 4 ) npe

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ido
  integer ( kind = 4 ) ielem
  integer ( kind = 4 ) indx(np,3)
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) ipp
  integer ( kind = 4 ) iq
  integer ( kind = 4 ) iqq
  integer ( kind = 4 ) iuk
  integer ( kind = 4 ) iukk
  integer ( kind = 4 ) j
  integer ( kind = 4 ) maxrow
  integer ( kind = 4 ) neqn
  integer ( kind = 4 ) nlband
  integer ( kind = 4 ) node(nelem,npe)
  integer ( kind = 4 ) nrow
!
  nlband = 0

  do ielem = 1, nelem
    do iq = 1, npe
      ip = node(ielem,iq)
      do iuk = 1, 3
        i = indx(ip,iuk)
        if ( i>0 ) then
          do iqq = 1, npe
            ipp = node(ielem,iqq)
            do iukk = 1, 3
              j = indx(ipp,iukk)
              if ( j>0 ) then
                if ( j-i>nlband ) then
                  nlband = j-i
                end if
              end if
            end do
          end do
        end if
      end do
    end do
  end do

  nrow = nlband+nlband+nlband+1

  ido = 0
  if ( ido == 1 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SETBAN:'
    write ( *, * ) '  Lower bandwidth NLBAND =       ',nlband
    write ( *, * ) '  Total bandwidth =            ',2*nlband+1
    write ( *, * ) '  Factor work =                ',2*neqn*nlband**2
    write ( *, * ) '  Solve work =                 ',3*neqn*nlband
    write ( *, * ) '  Required matrix rows NROW =    ',nrow
    write ( *, * ) '  Number of equations NEQN =   ',neqn
    write ( *, * ) '  Matrix storage 19 * NEQN =   ',19*neqn
    write ( *, * ) '  Matrix storage NROW * NEQN = ',nrow*neqn
    write ( *, * ) '  Matrix storage NEQN * NEQN = ',neqn*neqn
  end if

  if ( nrow>maxrow ) then
    write ( *, * ) ' '
    write ( *, * ) 'SETBAN - Fatal error!'
    write ( *, * ) '  NROW is too large!  NROW =      ',nrow
    write ( *, * ) '  The maximum allowed is MAXROW = ',maxrow
    stop
  end if

  return
end
subroutine setbas ( area,etaq,flarea,isotri,nelem,node,np,npe,phi,xc,xquad, &
  xsiq,yc,yquad)

!*****************************************************************************80
!
!! SETBAS evaluates the basis functions at each quadrature point.
!
!  Discussion:
!
!    The basis functions are computed and saved in this way for efficiency.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) AREA(NELEM).
!    AREA contains the area of each element.  These values are
!    needed when computed the integrals associated with the
!    finite element method.
!    For runs in which the region is allowed to change from
!    step to step, AREA must be recalculated at each step.
!
!    Input, integer ( kind = 4 ) ISOTRI(NELEM).
!    0, the element is NOT isoparametric.  The six node
!    triangle has straight sides.
!    1, the element is isoparametric.  The six node triangle
!    has curved sides.  Many computations involving such an
!    element must be computed by using a reference triangle,
!    and evaluating the jacobian of a transformation between
!    that triangle and the element.
!
!  NELEM  Input, integer ( kind = 4 ) NELEM, the number of elements.
!
!  NODE   Input, integer ( kind = 4 ) NODE(MAXELM,6), contains the numbers
!         of the nodes that make up each element.  Element number
!         I is associated with nodes NODE(I,1) through NODE(I,6).
!
!  NP     Input, integer ( kind = 4 ) NP, the number of nodes.
!
!  PHI    Output, real ( kind = 8 ) PHI(NELEM,3,6,6),
!         contains lots of basis function values.  In particular,
!
!           PHI(I,J,K,1) contains the value of the basis function
!           associated with velocity (U or V), in the I-th element,
!           at the J-th quadrature point, associated with the
!           K-th node.
!
!           PHI(I,J,K,2) contains the X derivative, and
!           PHI(I,J,K,3) contains the Y derivative.
!
!           PHI(I,J,K,4) contains the value of the basis function
!           associated with pressure (P) in the I-th element,
!           at the J-th quadrature point, associated with the
!           K-th node.
!
!  XC     Input, real ( kind = 8 ) XC(NP), contains the X coordinates
!         of the nodes.
!
!  XQUAD  Input, real ( kind = 8 ) XQUAD(NELEM,3), contains the
!         X coordinates  of the quadrature points in a given element.
!
!  YC     Input, real ( kind = 8 ) YC(NP), contains the Y coordinates
!         of the nodes.
!
!  YQUAD  Input, real ( kind = 8 ) YQUAD(NELEM,3), contains the
!         Y coordinates of the quadrature points in a given element.
!
  implicit none

  integer ( kind = 4 ) nelem
  integer ( kind = 4 ) np
  integer ( kind = 4 ) npe

  real ( kind = 8 ) area(nelem,3)
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
  real ( kind = 8 ) etaq(3)
  real ( kind = 8 ) flarea
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ielem
  integer ( kind = 4 ) iq
  integer ( kind = 4 ) isotri(nelem)
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) j
  logical lval
  integer ( kind = 4 ) node(nelem,npe)
  real ( kind = 8 ) phi(nelem,3,6,10)
  real ( kind = 8 ) q
  real ( kind = 8 ) w
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) xquad(nelem,3)
  real ( kind = 8 ) xq
  real ( kind = 8 ) xsi
  real ( kind = 8 ) xsiq(3)
  real ( kind = 8 ) yc(np)
  real ( kind = 8 ) yquad(nelem,3)
  real ( kind = 8 ) yq
!
  ival = 1
  call imemry('inc','SetBas_calls',ival)

  lval = .false.
  call lmemry('get','need_phi',lval)
!
!  Consider a particular element,
!  and a particular quadrature point (XQ,YQ) in that element.
!
!  Compute, at (XQ,YQ), the local values of the jacobian matrix
!  and its determinant.
!
!  Adjust the AREA array
!
  do ielem = 1,nelem

    if ( isotri(ielem) == 1 .or. isotri(ielem)==2.or.lval ) then

      do j = 1,3

        xq = xquad(ielem,j)
        yq = yquad(ielem,j)

        if ( isotri(ielem) == 2 ) then
          eta = etaq(j)
          xsi = xsiq(j)
          call trans ( det,detadx,detady,dxsidx,dxsidy,eta,ielem,nelem,node, &
            np,npe,xc,xsi,yc)
          area(ielem,j) = det*area(ielem,j)
        end if
!
!  Now consider each of the basis functions associated with a
!  node in the given element.
!
        do iq = 1,6
!
!  If the element is NOT isoparametric, compute the basis values
!  directly.
!
!  For isoparametric elements, use the reference triangle method.
!
          if ( isotri(ielem) == 0 .or. isotri(ielem)==1 ) then

            call bsp ( q,dqdx,dqdy,ielem,iq,nelem,node,np,npe,xc,xq,yc,yq)

            call qbf ( ielem,iq,w,dwdx,dwdy,nelem,node,np,npe,xc,xq,yc,yq)

            dxsidx = 1.0D+00
            dxsidy = 0.0D+00
            detadx = 0.0D+00
            detady = 1.0D+00

          else

            call refqbf ( w,dwdx,dwdy,detadx,detady,dxsidx,dxsidy,eta,iq,xsi)

            call refbsp ( q,dqdx,dqdy,detadx,detady,iq,dxsidx,dxsidy,eta,xsi)

          end if
!
!  Store the values into PHI.
!
          phi(ielem,j,iq,1) = w
          phi(ielem,j,iq,2) = dwdx
          phi(ielem,j,iq,3) = dwdy
          phi(ielem,j,iq,4) = q
          phi(ielem,j,iq,5) = dqdx
          phi(ielem,j,iq,6) = dqdy

          phi(ielem,j,iq,7) = dxsidx
          phi(ielem,j,iq,8) = dxsidy
          phi(ielem,j,iq,9) = detadx
          phi(ielem,j,iq,10) = detady

        end do
      end do

    end if

  end do
!
!  Compute the current flow area.
!
  flarea = 0.0D+00
  do ielem = 1,nelem
    do i = 1,3
      flarea = flarea+area(ielem,i)
    end do
  end do

  lval = .false.
  call lmemry('set','need_phi',lval)

  return
end
subroutine setnab ( nabor, np, nx, ny )

!*****************************************************************************80
!
!! SETNAB sets up the node neighbor array.
!
!  Discussion:
!
!    The neighbor array allows us to quickly find out, for each node,
!    the nodes that are its "neighbors", that is, which appear together
!    in some element.  There can be at most 18 such nodes,
!    but we include the node as its own neighbor, for a maximum of 19.
!
!    The computation is carried out in a straightforward way.  We rely
!    on the (contingent) facts that
!
!    * the grid is laid out in a rectangular fashion, with 2*NY-1 nodes
!      in the vertical direction and 2*NX-1 nodes in the horizontal direction.
!
!    * the rectangular grid is split into triangles by drawing a line
!      from south west to north east.
!
!    This algorithm produces the list of neighbors in order.
!
!    NABOR(K,I) = J means, if J is not zero, that node J is
!    (the K-th) neighbor of node I.
!
!     K     J           Direction
!
!     1     I-2*MY-2    SSWW
!     2     I-2*MY-1    SWW
!     3     I-2*MY      WW
!     4     I-MY-2      SSW
!     5     I-MY-1      SW
!     6     I-MY        W
!     7     I-MY+1      NW
!     8     I-2         SS
!     9     I-1         S
!    10     I           C
!    11     I+1         N
!    12     I+2         NN
!    13     I+MY-1      SE
!    14     I+MY        E
!    15     I+MY+1      NE
!    16     I+MY+2      NNE
!    17     I+2*MY      EE
!    18     I+2*MY+1    NEE
!    19     I+2*MY+2    NNEE
!
!    This pattern makes it possible to compute, from the fact that
!    J and I are neighbors, and the values of J and I alone, the
!    appropriate index of NABOR.  This in turn makes it easy to
!    assemble the matrix A in neighbor form.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) NABOR(19,NP), the node neigbors of each node.
!
!    Input, integer ( kind = 4 ) NP, the number of nodes.
!
!    Input, integer ( kind = 4 ) NX, NY, the number of elements in the X and Y
!    directions.  The number of nodes in a line in each direction
!    is 2*NX-1 and 2*NY-1.
!
  implicit none

  integer ( kind = 4 ) np

  integer ( kind = 4 ) i
  integer ( kind = 4 ) icol
  integer ( kind = 4 ) irow
  integer ( kind = 4 ) j
  integer ( kind = 4 ) mx
  integer ( kind = 4 ) my
  integer ( kind = 4 ) nabor(19,np)
  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny

  my = 2 * ny - 1
  mx = 2 * nx - 1

  do i = 1, np

    icol = ( i - 1 ) / my + 1
    irow = i - ( icol - 1 ) * my

    do j = 1, 19
      nabor(j,i) = 0
    end do
!
!  Example neighborhood of node 201:
!
!             203--213--223
!             /|        /|
!            / |       / |
!           /  |      /  |
!          /   |     /   |
!        191  202  212  222
!        /     |   /     |
!       /      |  /      |
!      /       | /       |
!     /        |/        |
!   181--190--201--211--221
!    |        /|        /
!    |       / |       /
!    |      /  |      /
!    |     /   |     /
!   180  189  200  210
!    |   /     |   /
!    |  /      |  /
!    | /       | /
!    |/        |/
!   179--188--199
!
    if ( mod ( icol, 2 ) == 1 .and. mod ( irow, 2 ) == 1 ) then

      if ( icol > 1 ) then

        if ( irow > 1 ) then
          nabor(1,i) = i - 2 * my - 2
          nabor(2,i) = i - 2 * my - 1
        end if

        nabor(3,i) = i - 2 * my

        if ( irow > 1 ) then
          nabor(4,i) = i - my - 2
          nabor(5,i) = i - my - 1
        end if

        nabor(6,i) = i - my

        if ( irow < my ) then
          nabor(7,i) = i - my + 1
        end if

      end if

      if ( irow > 1 ) then
        nabor(8,i) = i - 2
        nabor(9,i) = i - 1
      end if

      nabor(10,i) = i

      if ( irow < my ) then
        nabor(11,i) = i + 1
        nabor(12,i) = i + 2
      end if

      if ( icol < mx ) then

        if ( irow > 1 ) then
          nabor(13,i) = i + my - 1
        end if

        nabor(14,i) = i + my

        if ( irow < my ) then
          nabor(15,i) = i + my + 1
          nabor(16,i) = i + my + 2
        end if

        nabor(17,i) = i + 2 * my

        if ( irow < my ) then
          nabor(18,i) = i + 2 * my + 1
          nabor(19,i) = i + 2 * my + 2
        end if

      end if
!
!  Example neighborhood of node 202:
!
!             203--213--223
!             /|        /
!            / |       /
!           /  |      /
!          /   |     /
!        191  202  212
!        /     |   /
!       /      |  /
!      /       | /
!     /        |/
!   181--190--201
!
    else if ( mod ( icol, 2 ) == 1 .and. mod ( irow, 2 ) == 0 ) then

      if ( icol > 1 ) then
        nabor(2,i) = i - 2 * my - 1
        nabor(5,i) = i - my - 1
        nabor(6,i) = i - my
      end if

      nabor(9,i)  = i - 1
      nabor(10,i) = i
      nabor(11,i) = i + 1

      if ( icol < mx ) then
        nabor(14,i) = i + my
        nabor(15,i) = i + my + 1
        nabor(18,i) = i + 2 * my + 1
      end if
!
!  Example neighborhood of node 211:
!
!             223
!             /|
!            / |
!           /  |
!          /   |
!        212  222
!        /     |
!       /      |
!      /       |
!     /        |
!   201--211--221
!    |        /
!    |       /
!    |      /
!    |     /
!   200  210
!    |   /
!    |  /
!    | /
!    |/
!   199
!
    else if ( mod ( icol, 2 ) == 0 .and. mod ( irow, 2 ) == 1 ) then

       if ( irow > 1 ) then
         nabor(4,i) = i - my - 2
         nabor(5,i) = i - my - 1
       end if

       nabor(6,i) = i - my

       if ( irow > 1 ) then
         nabor(9,i) = i - 1
       end if

       nabor(10,i) = i

       if ( irow < my ) then
         nabor(11,i) = i + 1
       end if

       nabor(14,i) = i + my

       if ( irow < my ) then
         nabor(15,i) = i + my + 1
         nabor(16,i) = i + my + 2
       end if
!
!  Example neighborhood of node 212:
!
!   203--213--223
!    |        /|
!    |       / |
!    |      /  |
!    |     /   |
!   202  212  222
!    |   /     |
!    |  /      |
!    | /       |
!    |/        |
!   201--211--221
!
    else if ( mod ( icol, 2 ) == 0 .and. mod ( irow, 2 ) == 0 ) then

      nabor(5,i)  = i - my - 1
      nabor(6,i)  = i - my
      nabor(7,i)  = i - my + 1
      nabor(9,i)  = i - 1
      nabor(10,i) = i
      nabor(11,i) = i + 1
      nabor(13,i) = i + my - 1
      nabor(14,i) = i + my
      nabor(15,i) = i + my + 1

    end if

  end do

  return
end
subroutine setqxy ( area,etan,etaq,isotri,nelem,node,np,npe,wquad,xc,xquad,xsin, &
  xsiq,yc,yquad)

!*****************************************************************************80
!
!! SETQXY sets the abscissas and weights for a quadrature rule on a triangle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) AREA2(NELEM,3).
!
!    Output, real ( kind = 8 ) XQUAD(NELEM,3), YQUAD(NELEM,3), the X and Y coordinates
!    of the quadrature points for each element.
!
  implicit none

  integer ( kind = 4 ) nelem
  integer ( kind = 4 ) np
  integer ( kind = 4 ) npe

  real ( kind = 8 ) area(nelem,3)
  real ( kind = 8 ) eta
  real ( kind = 8 ) etan(6)
  real ( kind = 8 ) etaq(3)
  integer ( kind = 4 ) ielem
  integer ( kind = 4 ) ip1
  integer ( kind = 4 ) ip2
  integer ( kind = 4 ) ip3
  integer ( kind = 4 ) iquad
  integer ( kind = 4 ) isotri(nelem)
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) node(nelem,npe)
  real ( kind = 8 ) wquad(3)
  real ( kind = 8 ) x
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) xquad(nelem,3)
  real ( kind = 8 ) xsi
  real ( kind = 8 ) xsin(6)
  real ( kind = 8 ) xsiq(3)
  real ( kind = 8 ) y
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) y3
  real ( kind = 8 ) yc(np)
  real ( kind = 8 ) yquad(nelem,3)

  ival = 1
  call imemry('inc','SetQXY_calls',ival)
!
!  Set the weights.
!
  wquad(1) = 1.0D+00 / 6.0D+00
  wquad(2) = 1.0D+00 / 6.0D+00
  wquad(3) = 1.0D+00 / 6.0D+00
!
!  Set the quadrature points in the reference element.
!
  xsiq(1) = 0.5D+00
  etaq(1) = 0.5D+00
  xsiq(2) = 1.0D+00
  etaq(2) = 0.5D+00
  xsiq(3) = 0.5D+00
  etaq(3) = 0.0D+00
!
!  Set the X, Y coordinates of quadrature points for each element.
!
  do ielem = 1,nelem

    ip1 = node(ielem,1)
    ip2 = node(ielem,2)
    ip3 = node(ielem,3)

    x1 = xc(ip1)
    x2 = xc(ip2)
    x3 = xc(ip3)

    y1 = yc(ip1)
    y2 = yc(ip2)
    y3 = yc(ip3)

    if ( isotri(ielem) == 0 .or. isotri(ielem)==1 ) then

      xquad(ielem,1) = 0.5D+00 * (x1+x2)
      yquad(ielem,1) = 0.5D+00 * (y1+y2)

      xquad(ielem,2) = 0.5D+00 * (x2+x3)
      yquad(ielem,2) = 0.5D+00 * (y2+y3)

      xquad(ielem,3) = 0.5D+00 * (x3+x1)
      yquad(ielem,3) = 0.5D+00 * (y3+y1)

    else

      xsi = 0.5D+00*(xsin(1)+xsin(2))
      eta = 0.5D+00*(etan(1)+etan(2))
      call xofxsi ( eta,ielem,nelem,node,np,npe,x,xc,xsi,y,yc)
      xquad(ielem,1) = x
      yquad(ielem,1) = y

      xsi = 0.5D+00*(xsin(2)+xsin(3))
      eta = 0.5D+00*(etan(2)+etan(3))
      call xofxsi ( eta,ielem,nelem,node,np,npe,x,xc,xsi,y,yc)
      xquad(ielem,2) = x
      yquad(ielem,2) = y

      xsi = 0.5D+00*(xsin(3)+xsin(1))
      eta = 0.5D+00*(etan(3)+etan(1))
      call xofxsi ( eta,ielem,nelem,node,np,npe,x,xc,xsi,y,yc)
      xquad(ielem,3) = x
      yquad(ielem,3) = y

    end if
!
!  We only calculate true areas for nonisoparametric elements.
!
    do iquad = 1,3

      if ( isotri(ielem) == 0 .or. isotri(ielem)==1 ) then

        area(ielem,iquad) = abs((y1+y2)*(x2-x1)+(y2+y3)*(x3-x2) &
          +(y3+y1)*(x1-x3))*wquad(iquad)

      else

        area(ielem,iquad) = wquad(iquad)

      end if

    end do

  end do

  return
end
subroutine sgbscn(n,kl,ku,a,lda,nonzer,npiv,nzer)

!*****************************************************************************80
!
!! SGBSCN scans a matrix stored in LINPACK/LAPACK "general band" mode,
!  and returns the number of nonzero entries, and the number of zero
!  entries.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of columns of A.
!
!  KL     Input, integer ( kind = 4 ) KL, the number of subdiagonals within the
!         band of A.
!
!  KU     Input, integer ( kind = 4 ) KU, the number of superdiagonals within the
!         band of A.  KU >= 0.
!
!  A      Input, real ( kind = 8 ) A(LDA,N),
!         The matrix A in band storage, in rows KL+1 to
!         2*KL+KU+1; rows 1 to KL of the array need not be set.
!         The j-th column of A is stored in the j-th column of the
!         array AB as follows:
!
!           A(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
!
!  LDA    Input, integer ( kind = 4 ) LDA, the leading dimension of the array A.
!
!  NONZER Output, integer ( kind = 4 ) NONZER, the number of nonzero entries in
!         A.
!
!  NPIV   Output, integer ( kind = 4 ) NPIV, the number of pivot entries in A,
!         which is simply KL*N.
!
!  NZER   Output, integer ( kind = 4 ) NZER, the number of zero entries in A,
!         including the pivot entries.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(lda,n)
  integer ( kind = 4 ) irow
  integer ( kind = 4 ) j
  integer ( kind = 4 ) kl
  integer ( kind = 4 ) ku
  integer ( kind = 4 ) nonzer
  integer ( kind = 4 ) npiv
  integer ( kind = 4 ) nzer

  nonzer = 0
  npiv = kl*n
  nzer = kl*n

  do irow = kl+1,2*kl+ku+1
    do j = 1,n
      if ( a(irow,j)/= 0.0D+00 ) then
        nonzer = nonzer+1
      else
        nzer = nzer+1
      end if
    end do
  end do

  return
end
subroutine sgbtf2(m,n,kl,ku,ab,ldab,ipiv,info)

!*****************************************************************************80
!
!! SGBTF2 computes an LU factorization of a real m-by-n band matrix A
!  using partial pivoting with row interchanges.
!
!  This is the unblocked version of the algorithm, calling Level 2 BLAS.
!
!  Parameters:
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  KL      (input) INTEGER
!          The number of subdiagonals within the band of A.  KL >= 0.
!
!  KU      (input) INTEGER
!          The number of superdiagonals within the band of A.  KU >= 0.
!
!  AB      (input/output) real array, dimension (LDAB,N)
!          On entry, the matrix A in band storage, in rows KL+1 to
!          2*KL+KU+1; rows 1 to KL of the array need not be set.
!          The j-th column of A is stored in the j-th column of the
!          array AB as follows:
!          AB(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
!
!          On exit, details of the factorization: U is stored as an
!          upper triangular band matrix with KL+KU superdiagonals in
!          rows 1 to KL+KU+1, and the multipliers used during the
!          factorization are stored in rows KL+KU+2 to 2*KL+KU+1.
!          See below for further details.
!
!  LDAB    (input) INTEGER
!          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
!
!  IPIV    (output) INTEGER array, dimension (min(M,N))
!          The pivot indices; for 1  <=  i <= min(M,N), row i of the
!          matrix was interchanged with row IPIV(i).
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!          > 0: if INFO = +i, U(i,i) is exactly zero. The factorization
!               has been completed, but the factor U is exactly
!               singular, and division by zero will occur if it is used
!               to solve a system of equations.
!
!  Further Details
!
!  The band storage scheme is illustrated by the following example, when
!  M = N = 6, KL = 2, KU = 1:
!
!  On entry:                       On exit:
!
!      *    *    *    +    +    +       *    *    *   u14  u25  u36
!      *    *    +    +    +    +       *    *   u13  u24  u35  u46
!      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
!     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
!     a21  a32  a43  a54  a65   *      m21  m32  m43  m54  m65   *
!     a31  a42  a53  a64   *    *      m31  m42  m53  m64   *    *
!
!  Array elements marked * are not used by the routine; elements marked
!  + need not be set on entry, but are required by the routine to store
!  elements of U, because of fill-in resulting from the row
!  interchanges.
!
  implicit none

  integer ( kind = 4 )            info, kl, ku, ldab, m, n
!     ..
!     .. Array Arguments ..
  integer ( kind = 4 )            ipiv( * )
  real ( kind = 8 )   ab( ldab, * )
!     ..
!
  real ( kind = 8 )   one, zero
  parameter          ( one = 1.0, zero = 0.0D+00 )
!     ..
!     .. Local Scalars ..
  integer ( kind = 4 )            i, j, jp, ju, km, kv
!     ..
!     .. External Functions ..
  integer ( kind = 4 )            isamax
  external           isamax
!     ..
!     .. External Subroutines ..
  external           sger, sscal, sswap, xerbla
!
!     KV is the number of superdiagonals in the factor U, allowing for
!     fill-in.
!
  kv = ku + kl
!
!     Test the input parameters.
!
  info = 0
  if (  m<0 ) then
     info = -1
  else if (  n<0 ) then
     info = -2
  else if (  kl<0 ) then
     info = -3
  else if (  ku<0 ) then
     info = -4
  else if (  ldab<kl+kv+1 ) then
     info = -6
  end if
  if (  info/= 0 ) then
     call xerbla( 'sgbtf2', -info )
     return
  end if
!
!     Quick return if possible
!
  if (  m == 0 .or. n==0 ) return
!
!     Gaussian elimination with partial pivoting
!
!     Set fill-in elements in columns KU+2 to KV to zero.
!
  do j = ku + 2, min( kv, n )
    do i = kv - j + 2, kl
      ab( i, j ) = zero
    end do
  end do
!
!     JU is the index of the last column affected by the current stage
!     of the factorization.
!
  ju = 1
!
  do 40 j = 1, min( m, n )
!
!        Set fill-in elements in column J+KV to zero.
!
     if (  j+kv <= n ) then
        do 30 i = 1, kl
           ab( i, j+kv ) = zero
   30       continue
     end if
!
!        Find pivot and test for singularity. KM is the number of
!        subdiagonal elements in the current column.
!
     km = min( kl, m-j )
     jp = isamax( km+1, ab( kv+1, j ), 1 )
     ipiv( j ) = jp + j - 1
     if (  ab( kv+jp, j )/= zero ) then
        ju = max( ju, min( j+ku+jp-1, n ) )
!
!           Apply interchange to columns J to JU.
!
        if (  jp/= 1 ) call sswap( ju-j+1, ab( kv+jp, j ), ldab-1, &
                          ab( kv+1, j ), ldab-1 )

        if (  km>0 ) then
!
!              Compute multipliers.
!
           call sscal( km, one / ab( kv+1, j ), ab( kv+2, j ), 1 )
!
!              Update trailing submatrix within the band.
!
           if (  ju>j ) call sger( km, ju-j, -one, ab( kv+2, j ), 1, &
                            ab( kv, j+1 ), ldab-1, ab( kv+1, j+1 ),ldab-1 )
        end if
     else
!
!           If pivot is zero, set INFO to the index of the pivot
!           unless a zero pivot has already been found.
!
        if (  info == 0 )info = j
     end if
   40 continue
  return
end
subroutine sgbtrf(m,n,kl,ku,ab,ldab,ipiv,info)

!*****************************************************************************80
!
!! SGBTRF computes an LU factorization of a real m-by-n band matrix A
!  using partial pivoting with row interchanges.
!
!  This is the blocked version of the algorithm, calling Level 3 BLAS.
!
!  Parameters:
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  KL      (input) INTEGER
!          The number of subdiagonals within the band of A.  KL >= 0.
!
!  KU      (input) INTEGER
!          The number of superdiagonals within the band of A.  KU >= 0.
!
!  AB      (input/output) real array, dimension (LDAB,N)
!          On entry, the matrix A in band storage, in rows KL+1 to
!          2*KL+KU+1; rows 1 to KL of the array need not be set.
!          The j-th column of A is stored in the j-th column of the
!          array AB as follows:
!          AB(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
!
!          On exit, details of the factorization: U is stored as an
!          upper triangular band matrix with KL+KU superdiagonals in
!          rows 1 to KL+KU+1, and the multipliers used during the
!          factorization are stored in rows KL+KU+2 to 2*KL+KU+1.
!          See below for further details.
!
!  LDAB    (input) INTEGER
!          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
!
!  IPIV    (output) INTEGER array, dimension (min(M,N))
!          The pivot indices; for 1  <=  i <= min(M,N), row i of the
!          matrix was interchanged with row IPIV(i).
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!          > 0: if INFO = +i, U(i,i) is exactly zero. The factorization
!               has been completed, but the factor U is exactly
!               singular, and division by zero will occur if it is used
!               to solve a system of equations.
!
!  Further Details
!
!  The band storage scheme is illustrated by the following example, when
!  M = N = 6, KL = 2, KU = 1:
!
!  On entry:                       On exit:
!
!      *    *    *    +    +    +       *    *    *   u14  u25  u36
!      *    *    +    +    +    +       *    *   u13  u24  u35  u46
!      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
!     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
!     a21  a32  a43  a54  a65   *      m21  m32  m43  m54  m65   *
!     a31  a42  a53  a64   *    *      m31  m42  m53  m64   *    *
!
!  Array elements marked * are not used by the routine; elements marked
!  + need not be set on entry, but are required by the routine to store
!  elements of U because of fill-in resulting from the row interchanges.
!
  implicit none
  integer ( kind = 4 )            info, kl, ku, ldab, m, n
!     ..
!     .. Array Arguments ..
  integer ( kind = 4 )            ipiv( * )
  real ( kind = 8 )   ab( ldab, * )
!     ..
!
  real ( kind = 8 )   one, zero
  parameter          ( one = 1.0D+00, zero = 0.0D+00 )
  integer ( kind = 4 )            nbmax, ldwork
  parameter          ( nbmax = 64, ldwork = nbmax+1 )
!     ..
!     .. Local Scalars ..
  integer ( kind = 4 )            i, i2, i3, ii, ip, j, j2, j3, jb, jj, jm, jp
  integer ( kind = 4 )                   ju, k2, km, kv, nb, nw
  real ( kind = 8 )   temp
!     ..
!     .. Local Arrays ..
  real ( kind = 8 )   work13( ldwork, nbmax ),work31( ldwork, nbmax )
!     ..
!     .. External Functions ..
  integer ( kind = 4 )            isamax, ilaenv
  external           isamax, ilaenv
!     ..
!     .. External Subroutines ..
  external           scopy, sgbtf2, sgemm, sger, slaswp, sscal,sswap, strsm, xerbla
!
!     KV is the number of superdiagonals in the factor U, allowing for
!     fill-in
!
  kv = ku + kl
!
!     Test the input parameters.
!
  info = 0
  if (  m<0 ) then
     info = -1
  else if (  n<0 ) then
     info = -2
  else if (  kl<0 ) then
     info = -3
  else if (  ku<0 ) then
     info = -4
  else if (  ldab<kl+kv+1 ) then
     info = -6
  end if
  if (  info/= 0 ) then
     call xerbla( 'sgbtrf', -info )
     return
  end if
!
!     Quick return if possible
!
  if (  m == 0 .or. n==0 )return
!
!     Determine the block size for this environment
!
  nb = ilaenv( 1, 'sgbtrf', ' ', m, n, kl, ku )
!
!     The block size must not exceed the limit set by the size of the
!     local arrays WORK13 and WORK31.
!
  nb = min( nb, nbmax )
!
  if (  nb <= 1 .or. nb>kl ) then
!
!        Use unblocked code
!
     call sgbtf2( m, n, kl, ku, ab, ldab, ipiv, info )
  else
!
!        Use blocked code
!
!        Zero the superdiagonal elements of the work array WORK13
!
     do j = 1, nb
        do i = 1, j - 1
           work13( i, j ) = zero
        end do
     end do
!
!        Zero the subdiagonal elements of the work array WORK31
!
     do 40 j = 1, nb
        do 30 i = j + 1, nb
           work31( i, j ) = zero
   30       continue
   40    continue
!
!        Gaussian elimination with partial pivoting
!
!        Set fill-in elements in columns KU+2 to KV to zero
!
     do 60 j = ku + 2, min( kv, n )
        do 50 i = kv - j + 2, kl
           ab( i, j ) = zero
   50       continue
   60    continue
!
!        JU is the index of the last column affected by the current
!        stage of the factorization
!
     ju = 1
!
     do 180 j = 1, min( m, n ), nb
        jb = min( nb, min( m, n )-j+1 )
!
!           The active part of the matrix is partitioned
!
!              A11   A12   A13
!              A21   A22   A23
!              A31   A32   A33
!
!           Here A11, A21 and A31 denote the current block of JB columns
!           which is about to be factorized. The number of rows in the
!           partitioning are JB, I2, I3 respectively, and the numbers
!           of columns are JB, J2, J3. The superdiagonal elements of A13
!           and the subdiagonal elements of A31 lie outside the band.
!
        i2 = min( kl-jb, m-j-jb+1 )
        i3 = min( jb, m-j-kl+1 )
!
!           J2 and J3 are computed after JU has been updated.
!
!           Factorize the current block of JB columns
!
        do 80 jj = j, j + jb - 1
!
!              Set fill-in elements in column JJ+KV to zero
!
           if (  jj+kv <= n ) then
              do 70 i = 1, kl
                 ab( i, jj+kv ) = zero
   70             continue
           end if
!
!              Find pivot and test for singularity. KM is the number of
!              subdiagonal elements in the current column.
!
           km = min( kl, m-jj )
           jp = isamax( km+1, ab( kv+1, jj ), 1 )
           ipiv( jj ) = jp + jj - j
           if (  ab( kv+jp, jj )/= zero ) then
              ju = max( ju, min( jj+ku+jp-1, n ) )
              if (  jp/= 1 ) then
!
!                    Apply interchange to columns J to J+JB-1
!
                 if (  jp+jj-1<j+kl ) then
!
                    call sswap( jb, ab( kv+1+jj-j, j ), ldab-1, &
                                  ab( kv+jp+jj-j, j ), ldab-1 )
                 else
!
!                       The interchange affects columns J to JJ-1 of A31
!                       which are stored in the work array WORK31
!
                    call sswap( jj-j, ab( kv+1+jj-j, j ), ldab-1, &
                                  work31( jp+jj-j-kl, 1 ), ldwork )

                    call sswap( j+jb-jj, ab( kv+1, jj ), ldab-1, &
                                   ab( kv+jp, jj ), ldab-1 )
                 end if
              end if
!
!                 Compute multipliers
!
              call sscal( km, one / ab( kv+1, jj ), ab( kv+2, jj ), 1 )
!
!                 Update trailing submatrix within the band and within
!                 the current block. JM is the index of the last column
!                 which needs to be updated.
!
              jm = min( ju, j+jb-1 )
              if (  jm>jj ) call sger( km, jm-jj, -one, ab( kv+2, jj ), 1, &
                ab( kv, jj+1 ), ldab-1,ab( kv+1, jj+1 ), ldab-1 )
           else
!
!                 If pivot is zero, set INFO to the index of the pivot
!                 unless a zero pivot has already been found.
!
              if (  info == 0 ) info = jj
           end if
!
!              Copy current column of A31 into the work array WORK31
!
           nw = min( jj-j+1, i3 )
           if (  nw>0 ) call scopy( nw, ab( kv+kl+1-jj+j, jj ), 1, &
                             work31( 1, jj-j+1 ), 1 )
   80       continue
        if (  j+jb <= n ) then
!
!              Apply the row interchanges to the other blocks.
!
           j2 = min( ju-j+1, kv ) - jb
           j3 = max( 0, ju-j-kv+1 )
!
!              Use SLASWP to apply the row interchanges to A12, A22, and
!              A32.
!
           call slaswp( j2, ab( kv+1-jb, j+jb ), ldab-1, 1, jb,ipiv( j ), 1 )
!
!              Adjust the pivot indices.
!
           do i = j, j + jb - 1
              ipiv( i ) = ipiv( i ) + j - 1
           end do
!
!              Apply the row interchanges to A13, A23, and A33
!              columnwise.
!
           k2 = j - 1 + jb + j2
           do 110 i = 1, j3
              jj = k2 + i
              do 100 ii = j + i - 1, j + jb - 1
                 ip = ipiv( ii )
                 if (  ip/= ii ) then
                    temp = ab( kv+1+ii-jj, jj )
                    ab( kv+1+ii-jj, jj ) = ab( kv+1+ip-jj, jj )
                    ab( kv+1+ip-jj, jj ) = temp
                 end if
  100             continue
  110          continue
!
!              Update the relevant part of the trailing submatrix
!
           if (  j2>0 ) then
!
!                 Update A12
!
              call strsm( 'left', 'lower', 'no transpose', 'unit', &
                jb, j2, one, ab( kv+1, j ), ldab-1, &
                ab( kv+1-jb, j+jb ), ldab-1 )
!
              if (  i2>0 ) then
!
!                    Update A22
!
                 call sgemm( 'no transpose', 'no transpose', i2, j2, &
                   jb, -one, ab( kv+1+jb, j ), ldab-1, &
                   ab( kv+1-jb, j+jb ), ldab-1, one, &
                   ab( kv+1, j+jb ), ldab-1 )
              end if
!
              if (  i3>0 ) then
!
!                    Update A32
!
                 call sgemm( 'no transpose', 'no transpose', i3, j2, &
                                 jb, -one, work31, ldwork, &
                                 ab( kv+1-jb, j+jb ), ldab-1, one, &
                                 ab( kv+kl+1-jb, j+jb ), ldab-1 )
              end if
           end if
!
           if (  j3>0 ) then
!
!                 Copy the lower triangle of A13 into the work array
!                 WORK13
!
              do 130 jj = 1, j3
                 do 120 ii = jj, jb
                    work13( ii, jj ) = ab( ii-jj+1, jj+j+kv-1 )
  120                continue
  130             continue
!
!                 Update A13 in the work array
!
              call strsm( 'left', 'lower', 'no transpose', 'unit', &
                              jb, j3, one, ab( kv+1, j ), ldab-1, &
                              work13, ldwork )
!
              if (  i2>0 ) then
!
!                    Update A23
!
                 call sgemm( 'no transpose', 'no transpose', i2, j3, &
                                 jb, -one, ab( kv+1+jb, j ), ldab-1, &
                                 work13, ldwork, one, ab( 1+jb, j+kv ), &
                                 ldab-1 )
              end if
!
              if (  i3>0 ) then
!
!                    Update A33
!
                 call sgemm( 'no transpose', 'no transpose', i3, j3, &
                                 jb, -one, work31, ldwork, work13, &
                                 ldwork, one, ab( 1+kl, j+kv ), ldab-1 )
              end if
!
!                 Copy the lower triangle of A13 back into place
!
              do 150 jj = 1, j3
                 do 140 ii = jj, jb
                    ab( ii-jj+1, jj+j+kv-1 ) = work13( ii, jj )
  140                continue
  150             continue
           end if
        else
!
!              Adjust the pivot indices.
!
           do 160 i = j, j + jb - 1
              ipiv( i ) = ipiv( i ) + j - 1
  160          continue
        end if
!
!           Partially undo the interchanges in the current block to
!           restore the upper triangular form of A31 and copy the upper
!           triangle of A31 back into place
!
        do 170 jj = j + jb - 1, j, -1
           jp = ipiv( jj ) - jj + 1
           if (  jp/= 1 ) then
!
!                 Apply interchange to columns J to JJ-1
!
              if (  jp+jj-1<j+kl ) then
!
!                    The interchange does not affect A31
!
                 call sswap( jj-j, ab( kv+1+jj-j, j ), ldab-1, &
                               ab( kv+jp+jj-j, j ), ldab-1 )
              else
!
!                    The interchange does affect A31
!
                 call sswap( jj-j, ab( kv+1+jj-j, j ), ldab-1, &
                                 work31( jp+jj-j-kl, 1 ), ldwork )
              end if
           end if
!
!              Copy the current column of A31 back into place
!
           nw = min( i3, jj-j+1 )
           if (  nw>0 ) call scopy( nw, work31( 1, jj-j+1 ), 1, &
                              ab( kv+kl+1-jj+j, jj ), 1 )
  170       continue
  180    continue
  end if
!
  return
end
subroutine sgbtrs(trans,n,kl,ku,nrhs,ab,ldab,ipiv,b,ldb,info)

!*****************************************************************************80
!
!! SGBTRS solves a system of linear equations
!     A * X = B  or  A' * X = B
!  with a general band matrix A using the LU factorization computed
!  by SGBTRF.
!
!  Parameters:
!
!  TRANS   (input) CHARACTER*1
!          Specifies the form of the system of equations.
!          = 'N':  A * X = B  (No transpose)
!          = 'T':  A'* X = B  (Transpose)
!          = 'C':  A'* X = B  (Conjugate transpose = Transpose)
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  KL      (input) INTEGER
!          The number of subdiagonals within the band of A.  KL >= 0.
!
!  KU      (input) INTEGER
!          The number of superdiagonals within the band of A.  KU >= 0.
!
!  NRHS    (input) INTEGER
!          The number of right hand sides, i.e., the number of columns
!          of the matrix B.  NRHS >= 0.
!
!  AB      (input) real array, dimension (LDAB,N)
!          Details of the LU factorization of the band matrix A, as
!          computed by SGBTRF.  U is stored as an upper triangular band
!          matrix with KL+KU superdiagonals in rows 1 to KL+KU+1, and
!          the multipliers used during the factorization are stored in
!          rows KL+KU+2 to 2*KL+KU+1.
!
!  LDAB    (input) INTEGER
!          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
!
!  IPIV    (input) INTEGER array, dimension (N)
!          The pivot indices; for 1  <=  i <= N, row i of the matrix was
!          interchanged with row IPIV(i).
!
!  B       (input/output) real array, dimension (LDB,NRHS)
!          On entry, the right hand side vectors B for the system of
!          linear equations.
!          On exit, the solution vectors, X.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0: if INFO = -k, the k-th argument had an illegal value
!
  implicit none
  character          trans
  integer ( kind = 4 )            info, kl, ku, ldab, ldb, n, nrhs
!     ..
!     .. Array Arguments ..
  integer ( kind = 4 )            ipiv( * )
  real ( kind = 8 )   ab( ldab, * ), b( ldb, * )

  real ( kind = 8 )   one
  parameter          ( one = 1.0D+00 )
!     ..
!     .. Local Scalars ..
  logical            lnoti, notran
  integer ( kind = 4 )            i, j, kd, l, lm
!     ..
!     .. External Functions ..
  logical            lsame
  external           lsame
!     ..
!     .. External Subroutines ..
  external           sgemv, sger, sswap, stbsv, xerbla
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
  info = 0
  notran = lsame( trans, 'n' )
  if (  .not.notran .and. .not.lsame( trans, 't' ) .and. .not. &
        lsame( trans, 'c' ) ) then
     info = -1
  else if (  n<0 ) then
     info = -2
  else if (  kl<0 ) then
     info = -3
  else if (  ku<0 ) then
     info = -4
  else if (  nrhs<0 ) then
     info = -5
  else if (  ldab<( 2*kl+ku+1 ) ) then
     info = -7
  else if (  ldb<max( 1, n ) ) then
     info = -10
  end if
  if (  info/= 0 ) then
     call xerbla( 'sgbtrs', -info )
     return
  end if
!
!     Quick return if possible
!
  if (  n == 0 .or. nrhs==0 ) return
!
  kd = ku + kl + 1
  lnoti = kl>0
!
  if (  notran ) then
!
!        Solve  A*X = B.
!
!        Solve L*X = B, overwriting B with X.
!
!        L is represented as a product of permutations and unit lower
!        triangular matrices L = P(1) * L(1) * ... * P(n-1) * L(n-1),
!        where each transformation L(i) is a rank-one modification of
!        the identity matrix.
!
     if (  lnoti ) then
        do 10 j = 1, n - 1
           lm = min( kl, n-j )
           l = ipiv( j )
           if (  l/= j ) call sswap( nrhs, b( l, 1 ), ldb, b( j, 1 ), ldb )
           call sger( lm, nrhs, -one, ab( kd+1, j ), 1, b( j, 1 ), &
                         ldb, b( j+1, 1 ), ldb )
   10       continue
     end if
!
     do 20 i = 1, nrhs
!
!           Solve U*X = B, overwriting B with X.
!
        call stbsv( 'upper', 'no transpose', 'non-unit', n, kl+ku, &
                      ab, ldab, b( 1, i ), 1 )
   20    continue
!
  else
!
!        Solve A'*X = B.
!
     do 30 i = 1, nrhs
!
!           Solve U'*X = B, overwriting B with X.
!
        call stbsv( 'upper', 'transpose', 'non-unit', n, kl+ku, ab, &
                      ldab, b( 1, i ), 1 )
   30    continue
!
!        Solve L'*X = B, overwriting B with X.
!
     if (  lnoti ) then
        do 40 j = n - 1, 1, -1
           lm = min( kl, n-j )
           call sgemv( 'transpose', lm, nrhs, -one, b( j+1, 1 ), &
                         ldb, ab( kd+1, j ), 1, one, b( j, 1 ), ldb )
           l = ipiv( j )
           if (  l/= j )call sswap( nrhs, b( l, 1 ), ldb, b( j, 1 ), ldb )
   40       continue
     end if
  end if
  return
end
subroutine sgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)

!*****************************************************************************80
!
!! SGEMM  performs one of the matrix-matrix operations
!
!     C : =  alpha*op( A )*op( B ) + beta*C,
!
!  where  op( X ) is one of
!
!     op( X ) = X   or   op( X ) = X',
!
!  alpha and beta are scalars, and A, B and C are matrices, with op( A )
!  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
!
!  Parameters
!
!  TRANSA - CHARACTER*1.
!           On entry, TRANSA specifies the form of op( A ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSA = 'N' or 'n',  op( A ) = A.
!
!              TRANSA = 'T' or 't',  op( A ) = A'.
!
!              TRANSA = 'C' or 'c',  op( A ) = A'.
!
!           Unchanged on exit.
!
!  TRANSB - CHARACTER*1.
!           On entry, TRANSB specifies the form of op( B ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSB = 'N' or 'n',  op( B ) = B.
!
!              TRANSB = 'T' or 't',  op( B ) = B'.
!
!              TRANSB = 'C' or 'c',  op( B ) = B'.
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry,  M  specifies  the number  of rows  of the  matrix
!           op( A )  and of the  matrix  C.  M  must  be at least  zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry,  N  specifies the number  of columns of the matrix
!           op( B ) and the number of columns of the matrix C. N must be
!           at least zero.
!           Unchanged on exit.
!
!  K      - INTEGER.
!           On entry,  K  specifies  the number of columns of the matrix
!           op( A ) and the number of rows of the matrix op( B ). K must
!           be at least  zero.
!           Unchanged on exit.
!
!  ALPHA  - real.
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - real array of DIMENSION ( LDA, ka ), where ka is
!           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
!           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
!           part of the array  A  must contain the matrix  A,  otherwise
!           the leading  k by m  part of the array  A  must contain  the
!           matrix A.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
!           LDA must be at least  max( 1, m ), otherwise  LDA must be at
!           least  max( 1, k ).
!           Unchanged on exit.
!
!  B      - real array of DIMENSION ( LDB, kb ), where kb is
!           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
!           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
!           part of the array  B  must contain the matrix  B,  otherwise
!           the leading  n by k  part of the array  B  must contain  the
!           matrix B.
!           Unchanged on exit.
!
!  LDB    - INTEGER.
!           On entry, LDB specifies the first dimension of B as declared
!           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
!           LDB must be at least  max( 1, k ), otherwise  LDB must be at
!           least  max( 1, n ).
!           Unchanged on exit.
!
!  BETA   - real.
!           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
!           supplied as zero then C need not be set on input.
!           Unchanged on exit.
!
!  C      - real array of DIMENSION ( LDC, n ).
!           Before entry, the leading  m by n  part of the array  C must
!           contain the matrix  C,  except when  beta  is zero, in which
!           case C need not be set on entry.
!           On exit, the array  C  is overwritten by the  m by n  matrix
!           ( alpha*op( A )*op( B ) + beta*C ).
!
!  LDC    - INTEGER.
!           On entry, LDC specifies the first dimension of C as declared
!           in  the  calling  (sub)  program.   LDC  must  be  at  least
!           max( 1, m ).
!           Unchanged on exit.
!
!
!  Level 3 Blas routine.
!
!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
!
!
  implicit none
  character        transa, transb
  integer ( kind = 4 )            m, n, k, lda, ldb, ldc
  real ( kind = 8 )   alpha, beta
!     .. Array Arguments ..
  real ( kind = 8 )   a( lda, * ), b( ldb, * ), c( ldc, * )
!     ..
!
  logical            lsame
  external           lsame
!     .. External Subroutines ..
  external           xerbla
!
!     .. Local Scalars ..
  logical            nota, notb
  integer ( kind = 4 )            i, info, j, l, nrowa, nrowb
  real ( kind = 8 )   temp
!     .. Parameters ..
  real ( kind = 8 )   one         , zero
  parameter        ( one = 1.0D+00, zero = 0.0D+00 )
!     ..
!     .. Executable Statements ..
!
!     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
!     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
!     and  columns of  A  and the  number of  rows  of  B  respectively.
!
  nota  = lsame( transa, 'n' )
  notb  = lsame( transb, 'n' )
  if (  nota  ) then
     nrowa = m
  else
     nrowa = k
  end if
  if (  notb  ) then
     nrowb = k
  else
     nrowb = n
  end if
!
!     Test the input parameters.
!
  info = 0
  if (       ( .not.nota                 ) .and.  &
             ( .not.lsame( transa, 'c' ) ) .and.  &
              ( .not.lsame( transa, 't' ) )       ) then
     info = 1
  else if (  ( .not.notb                 ) .and.  &
             ( .not.lsame( transb, 'c' ) ) .and.  &
             ( .not.lsame( transb, 't' ) )       ) then
     info = 2
  else if (  m  <0                ) then
     info = 3
  else if (  n  <0                ) then
     info = 4
  else if (  k  <0                ) then
     info = 5
  else if (  lda<max( 1, nrowa )  ) then
     info = 8
  else if (  ldb<max( 1, nrowb )  ) then
     info = 10
  else if (  ldc<max( 1, m     )  ) then
     info = 13
  end if
  if (  info/= 0  ) then
     call xerbla( 'sgemm ', info )
     return
  end if
!
!     Quick return if possible.
!
  if (  ( m == 0 ) .or. ( n==0 ).or. &
        ( ( ( alpha == zero ) .or. ( k==0 ) ) .and. ( beta==one ) ) ) return
!
!     And if  alpha == zero.
!
  if (  alpha == zero  ) then
     if (  beta == zero  ) then
        do 20, j = 1, n
           do 10, i = 1, m
              c( i, j ) = zero
   10          continue
   20       continue
     else
        do 40, j = 1, n
           do 30, i = 1, m
              c( i, j ) = beta*c( i, j )
   30          continue
   40       continue
     end if
     return
  end if
!
!     Start the operations.
!
  if (  notb  ) then
     if (  nota  ) then
!
!           Form  C : =  alpha*A*B + beta*C.
!
        do 90, j = 1, n
           if (  beta == zero  ) then
              do 50, i = 1, m
                 c( i, j ) = zero
   50             continue
           else if (  beta/= one  ) then
              do 60, i = 1, m
                 c( i, j ) = beta*c( i, j )
   60             continue
           end if
           do 80, l = 1, k
              if (  b( l, j )/= zero  ) then
                 temp = alpha*b( l, j )
                 do 70, i = 1, m
                    c( i, j ) = c( i, j ) + temp*a( i, l )
   70                continue
              end if
   80          continue
   90       continue
     else
!
!           Form  C : =  alpha*A'*B + beta*C
!
        do 120, j = 1, n
           do 110, i = 1, m
              temp = zero
              do 100, l = 1, k
                 temp = temp + a( l, i )*b( l, j )
  100             continue
              if (  beta == zero  ) then
                 c( i, j ) = alpha*temp
              else
                 c( i, j ) = alpha*temp + beta*c( i, j )
              end if
  110          continue
  120       continue
     end if
  else
     if (  nota  ) then
!
!           Form  C : =  alpha*A*B' + beta*C
!
        do 170, j = 1, n
           if (  beta == zero  ) then
              do 130, i = 1, m
                 c( i, j ) = zero
  130             continue
           else if (  beta/= one  ) then
              do 140, i = 1, m
                 c( i, j ) = beta*c( i, j )
  140             continue
           end if
           do 160, l = 1, k
              if (  b( j, l )/= zero  ) then
                 temp = alpha*b( j, l )
                 do 150, i = 1, m
                    c( i, j ) = c( i, j ) + temp*a( i, l )
  150                continue
              end if
  160          continue
  170       continue
     else
!
!           Form  C : =  alpha*A'*B' + beta*C
!
        do 200, j = 1, n
           do 190, i = 1, m
              temp = zero
              do 180, l = 1, k
                 temp = temp + a( l, i )*b( j, l )
  180             continue
              if (  beta == zero  ) then
                 c( i, j ) = alpha*temp
              else
                 c( i, j ) = alpha*temp + beta*c( i, j )
              end if
  190          continue
  200       continue
     end if
  end if
!
  return
end
subroutine sgemv(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)

!*****************************************************************************80
!
!! SGEMV  performs one of the matrix-vector operations
!
!     y : =  alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
!
!  where alpha and beta are scalars, x and y are vectors and A is an
!  m by n matrix.
!
!  Parameters
!   == ========
!
!  TRANS  - CHARACTER*1.
!           On entry, TRANS specifies the operation to be performed as
!           follows:
!
!              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
!
!              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
!
!              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry, M specifies the number of rows of the matrix A.
!           M must be at least zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - real.
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - real array of DIMENSION ( LDA, n ).
!           Before entry, the leading m by n part of the array A must
!           contain the matrix of coefficients.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, m ).
!           Unchanged on exit.
!
!  X      - real array of DIMENSION at least
!           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
!           and at least
!           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
!           Before entry, the incremented array X must contain the
!           vector x.
!           Unchanged on exit.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  BETA   - real.
!           On entry, BETA specifies the scalar beta. When BETA is
!           supplied as zero then Y need not be set on input.
!           Unchanged on exit.
!
!  Y      - real array of DIMENSION at least
!           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
!           and at least
!           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
!           Before entry with BETA non-zero, the incremented array Y
!           must contain the vector y. On exit, Y is overwritten by the
!           updated vector y.
!
!  INCY   - INTEGER.
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!           Unchanged on exit.
!
!
!  Level 2 Blas routine.
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
!
  implicit none
  real ( kind = 8 )   alpha, beta
  integer ( kind = 4 )            incx, incy, lda, m, n
  character        trans
!     .. Array Arguments ..
  real ( kind = 8 )   a( lda, * ), x( * ), y( * )

  real ( kind = 8 )   one         , zero
  parameter        ( one = 1.0D+00, zero = 0.0D+00 )
!     .. Local Scalars ..
  real ( kind = 8 )   temp
  integer ( kind = 4 )            i, info, ix, iy, j, jx, jy, kx, ky, lenx, leny
!     .. External Functions ..
  logical            lsame
  external           lsame
!     .. External Subroutines ..
  external           xerbla
!
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
  info = 0
  if     ( .not.lsame( trans, 'n' ) .and.  &
             .not.lsame( trans, 't' ) .and.  &
             .not.lsame( trans, 'c' )       ) then
     info = 1
  else if (  m<0  ) then
     info = 2
  else if (  n<0  ) then
     info = 3
  else if (  lda<max( 1, m )  ) then
     info = 6
  else if (  incx == 0  ) then
     info = 8
  else if (  incy == 0  ) then
     info = 11
  end if
  if (  info/= 0  ) then
     call xerbla( 'sgemv ', info )
     return
  end if
!
!     Quick return if possible.
!
  if (  ( m == 0 ) .or. ( n==0 ).or. &
        ( ( alpha == zero ) .and. ( beta==one ) ) )return
!
!     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
!     up the start points in  X  and  Y.
!
  if (  lsame( trans, 'n' )  ) then
     lenx = n
     leny = m
  else
     lenx = m
     leny = n
  end if
  if (  incx>0  ) then
     kx = 1
  else
     kx = 1 - ( lenx - 1 )*incx
  end if
  if (  incy>0  ) then
     ky = 1
  else
     ky = 1 - ( leny - 1 )*incy
  end if
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
!     First form  y : =  beta*y.
!
  if (  beta/= one  ) then
     if (  incy == 1  ) then
        if (  beta == zero  ) then
           do 10, i = 1, leny
              y( i ) = zero
   10          continue
        else
           do 20, i = 1, leny
              y( i ) = beta*y( i )
   20          continue
        end if
     else
        iy = ky
        if (  beta == zero  ) then
           do 30, i = 1, leny
              y( iy ) = zero
              iy      = iy   + incy
   30          continue
        else
           do 40, i = 1, leny
              y( iy ) = beta*y( iy )
              iy      = iy           + incy
   40          continue
        end if
     end if
  end if
  if (  alpha == zero )return
  if (  lsame( trans, 'n' )  ) then
!
!        Form  y : =  alpha*A*x + y.
!
     jx = kx
     if (  incy == 1  ) then
        do 60, j = 1, n
           if (  x( jx )/= zero  ) then
              temp = alpha*x( jx )
              do 50, i = 1, m
                 y( i ) = y( i ) + temp*a( i, j )
   50             continue
           end if
           jx = jx + incx
   60       continue
     else
        do 80, j = 1, n
           if (  x( jx )/= zero  ) then
              temp = alpha*x( jx )
              iy   = ky
              do 70, i = 1, m
                 y( iy ) = y( iy ) + temp*a( i, j )
                 iy      = iy      + incy
   70             continue
           end if
           jx = jx + incx
   80       continue
     end if
  else
!
!        Form  y : =  alpha*A'*x + y.
!
     jy = ky
     if (  incx == 1  ) then
        do 100, j = 1, n
           temp = zero
           do 90, i = 1, m
              temp = temp + a( i, j )*x( i )
   90          continue
           y( jy ) = y( jy ) + alpha*temp
           jy      = jy      + incy
  100       continue
     else
        do 120, j = 1, n
           temp = zero
           ix   = kx
           do 110, i = 1, m
              temp = temp + a( i, j )*x( ix )
              ix   = ix   + incx
  110          continue
           y( jy ) = y( jy ) + alpha*temp
           jy      = jy      + incy
  120       continue
     end if
  end if
!
  return
end
subroutine sger( m, n, alpha, x, incx, y, incy, a, lda )

!*****************************************************************************80
!
!! SGER   performs the rank 1 operation
!
!     A : =  alpha*x*y' + A,
!
!  where alpha is a scalar, x is an m element vector, y is an n element
!  vector and A is an m by n matrix.
!
!  Parameters
!   == ========
!
!  M      - INTEGER.
!           On entry, M specifies the number of rows of the matrix A.
!           M must be at least zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - real.
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  X      - real array of dimension at least
!           ( 1 + ( m - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the m
!           element vector x.
!           Unchanged on exit.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  Y      - real array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCY ) ).
!           Before entry, the incremented array Y must contain the n
!           element vector y.
!           Unchanged on exit.
!
!  INCY   - INTEGER.
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!           Unchanged on exit.
!
!  A      - real array of DIMENSION ( LDA, n ).
!           Before entry, the leading m by n part of the array A must
!           contain the matrix of coefficients. On exit, A is
!           overwritten by the updated matrix.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, m ).
!           Unchanged on exit.
!
!
!  Level 2 Blas routine.
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
!
  implicit none
  real ( kind = 8 )   alpha
  integer ( kind = 4 )            incx, incy, lda, m, n
!     .. Array Arguments ..
  real ( kind = 8 )   a( lda, * ), x( * ), y( * )

  real ( kind = 8 )   zero
  parameter        ( zero = 0.0D+00 )
!     .. Local Scalars ..
  real ( kind = 8 )   temp
  integer ( kind = 4 )            i, info, ix, j, jy, kx
!     .. External Subroutines ..
  external           xerbla
!
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
  info = 0
  if     ( m<0  ) then
     info = 1
  else if (  n<0  ) then
     info = 2
  else if (  incx == 0  ) then
     info = 5
  else if (  incy == 0  ) then
     info = 7
  else if (  lda<max( 1, m )  ) then
     info = 9
  end if
  if (  info/= 0  ) then
     call xerbla( 'sger  ', info )
     return
  end if
!
!     Quick return if possible.
!
  if (  ( m == 0 ) .or. ( n==0 ).or.( alpha==zero ) ) return
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
  if (  incy>0  ) then
     jy = 1
  else
     jy = 1 - ( n - 1 )*incy
  end if
  if (  incx == 1  ) then
     do 20, j = 1, n
        if (  y( jy )/= zero  ) then
           temp = alpha*y( jy )
           do 10, i = 1, m
              a( i, j ) = a( i, j ) + x( i )*temp
   10          continue
        end if
        jy = jy + incy
   20    continue
  else
     if (  incx>0  ) then
        kx = 1
     else
        kx = 1 - ( m - 1 )*incx
     end if
     do 40, j = 1, n
        if (  y( jy )/= zero  ) then
           temp = alpha*y( jy )
           ix   = kx
           do 30, i = 1, m
              a( i, j ) = a( i, j ) + x( ix )*temp
              ix        = ix        + incx
   30          continue
        end if
        jy = jy + incy
   40    continue
  end if
!
  return
!
!     End of SGER  .
!
end
subroutine sgetf2(m,n,a,lda,ipiv,info)

!*****************************************************************************80
!
!! SGETF2 computes an LU factorization of a general m-by-n matrix A
!  using partial pivoting with row interchanges.
!
!  The factorization has the form
!     A = P * L * U
!  where P is a permutation matrix, L is lower triangular with unit
!  diagonal elements (lower trapezoidal if m > n), and U is upper
!  triangular (upper trapezoidal if m < n).
!
!  This is the right-looking Level 2 BLAS version of the algorithm.
!
!  Arguments
!   == =======
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input/output) real array, dimension (LDA,N)
!          On entry, the m by n matrix to be factored.
!          On exit, the factors L and U from the factorization
!          A = P*L*U; the unit diagonal elements of L are not stored.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  IPIV    (output) INTEGER array, dimension (min(M,N))
!          The pivot indices; for 1  <=  i <= min(M,N), row i of the
!          matrix was interchanged with row IPIV(i).
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -k, the k-th argument had an illegal value
!          > 0: if INFO = k, U(k,k) is exactly zero. The factorization
!               has been completed, but the factor U is exactly
!               singular, and division by zero will occur if it is used
!               to solve a system of equations.
!
  implicit none
  INTEGER            INFO, LDA, M, N
!     ..
!     .. Array Arguments ..
  INTEGER            IPIV( * )
  real ( kind = 8 )   A( LDA, * )
  real ( kind = 8 )   ONE, ZERO
  PARAMETER          ( ONE = 1.0D+00, ZERO = 0.0D+00 )
!     ..
!     .. Local Scalars ..
  INTEGER            J, JP
!     ..
!     .. External Functions ..
  INTEGER            isamax
  EXTERNAL           isamax
!     ..
!     .. External Subroutines ..
  EXTERNAL           SGER, SSCAL, SSWAP, XERBLA
!     ..
!     .. Intrinsic Functions ..
  INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
  INFO = 0
  IF( M.LT.0 ) THEN
     INFO = -1
  ELSE IF( N.LT.0 ) THEN
     INFO = -2
  ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
     INFO = -4
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'SGETF2', -INFO )
     RETURN
  END IF
!
!     Quick return if possible
!
  IF( M.EQ.0 .OR. N.EQ.0 )RETURN
!
  DO 10 J = 1, MIN( M, N )
!
!        Find pivot and test for singularity.
!
     JP = J - 1 + isamax( M-J+1, A( J, J ), 1 )
     IPIV( J ) = JP
     IF( A( JP, J ).NE.ZERO ) THEN
!
!           Apply the interchange to columns 1:N.
!
        IF( JP.NE.J )CALL SSWAP( N, A( J, 1 ), LDA, A( JP, 1 ), LDA )
!
!           Compute elements J+1:M of J-th column.
!
        IF( J.LT.M )CALL SSCAL( M-J, ONE / A( J, J ), A( J+1, J ), 1 )
!
     ELSE IF( INFO.EQ.0 ) THEN
!
        INFO = J
     END IF
!
     IF( J.LT.MIN( M, N ) ) THEN
!
!           Update trailing submatrix.
!
        CALL SGER( M-J, N-J, -ONE, A( J+1, J ), 1, A( J, J+1 ), LDA, &
                    A( J+1, J+1 ), LDA )
     END IF
   10 CONTINUE
  RETURN
!
!     End of SGETF2
!
  END
subroutine sgetrf(m,n,a,lda,ipiv,info)

!*****************************************************************************80
!
!! SGETRF computes an LU factorization of a general M-by-N matrix A
!  using partial pivoting with row interchanges.
!
!  The factorization has the form
!     A = P * L * U
!  where P is a permutation matrix, L is lower triangular with unit
!  diagonal elements (lower trapezoidal if m > n), and U is upper
!  triangular (upper trapezoidal if m < n).
!
!  This is the right-looking Level 3 BLAS version of the algorithm.
!
!  Arguments
!   == =======
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input/output) real array, dimension (LDA,N)
!          On entry, the M-by-N matrix to be factored.
!          On exit, the factors L and U from the factorization
!          A = P*L*U; the unit diagonal elements of L are not stored.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  IPIV    (output) INTEGER array, dimension (min(M,N))
!          The pivot indices; for 1  <=  i <= min(M,N), row i of the
!          matrix was interchanged with row IPIV(i).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
!                has been completed, but the factor U is exactly
!                singular, and division by zero will occur if it is used
!                to solve a system of equations.
!
  implicit none

  INTEGER            INFO, LDA, M, N
!     ..
!     .. Array Arguments ..
  INTEGER            IPIV( * )
  real ( kind = 8 )   A( LDA, * )
  real ( kind = 8 )   ONE
  PARAMETER          ( ONE = 1.0D+00 )
!     ..
!     .. Local Scalars ..
  INTEGER            I, IINFO, J, JB, NB
!     ..
!     .. External Subroutines ..
  EXTERNAL           SGEMM, SGETF2, SLASWP, STRSM, XERBLA
!     ..
!     .. External Functions ..
  INTEGER            ILAENV
  EXTERNAL           ILAENV
!     ..
!     .. Intrinsic Functions ..
  INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
  INFO = 0
  IF( M.LT.0 ) THEN
     INFO = -1
  ELSE IF( N.LT.0 ) THEN
     INFO = -2
  ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
     INFO = -4
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'SGETRF', -INFO )
     RETURN
  END IF
!
!     Quick return if possible
!
  IF( M.EQ.0 .OR. N.EQ.0 )RETURN
!
!     Determine the block size for this environment.
!
  NB = ILAENV( 1, 'SGETRF', ' ', M, N, -1, -1 )
  IF( NB.LE.1 .OR. NB.GE.MIN( M, N ) ) THEN
!
!        Use unblocked code.
!
     CALL SGETF2( M, N, A, LDA, IPIV, INFO )
  ELSE
!
!        Use blocked code.
!
     DO 20 J = 1, MIN( M, N ), NB
        JB = MIN( MIN( M, N )-J+1, NB )
!
!           Factor diagonal and subdiagonal blocks and test for exact
!           singularity.
!
        CALL SGETF2( M-J+1, JB, A( J, J ), LDA, IPIV( J ), IINFO )
!
!           Adjust INFO and the pivot indices.
!
        IF( INFO.EQ.0 .AND. IINFO.GT.0 )INFO = IINFO + J - 1
        DO 10 I = J, MIN( M, J+JB-1 )
           IPIV( I ) = J - 1 + IPIV( I )
   10       CONTINUE
!
!           Apply interchanges to columns 1:J-1.
!
        CALL SLASWP( J-1, A, LDA, J, J+JB-1, IPIV, 1 )
!
        IF( J+JB.LE.N ) THEN
!
!              Apply interchanges to columns J+JB:N.
!
           CALL SLASWP( N-J-JB+1, A( 1, J+JB ), LDA, J, J+JB-1,IPIV, 1 )
!
!              Compute block row of U.
!
           CALL STRSM( 'Left', 'Lower', 'No transpose', 'Unit', JB, &
             N-J-JB+1, ONE, A( J, J ), LDA, A( J, J+JB ),LDA )
           IF( J+JB.LE.M ) THEN
!
!                 Update trailing submatrix.
!
              CALL SGEMM( 'No transpose', 'No transpose', M-J-JB+1, &
                N-J-JB+1, JB, -ONE, A( J+JB, J ), LDA, &
                A( J, J+JB ), LDA, ONE, A( J+JB, J+JB ),LDA )
           END IF
        END IF
   20    CONTINUE
  END IF
  RETURN
END
subroutine sgetrs(trans,n,nrhs,a,lda,ipiv,b,ldb,info)

!*****************************************************************************80
!
!! SGETRS solves a system of linear equations
!     A * X = B  or  A' * X = B
!  with a general N-by-N matrix A using the LU factorization computed
!  by SGETRF.
!
!  Arguments
!   == =======
!
!  TRANS   (input) CHARACTER*1
!          Specifies the form of the system of equations:
!          = 'N':  A * X = B  (No transpose)
!          = 'T':  A'* X = B  (Transpose)
!          = 'C':  A'* X = B  (Conjugate transpose = Transpose)
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  NRHS    (input) INTEGER
!          The number of right hand sides, i.e., the number of columns
!          of the matrix B.  NRHS >= 0.
!
!  A       (input) real array, dimension (LDA,N)
!          The factors L and U from the factorization A = P*L*U
!          as computed by SGETRF.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  IPIV    (input) INTEGER array, dimension (N)
!          The pivot indices from SGETRF; for 1 <= i<=N, row i of the
!          matrix was interchanged with row IPIV(i).
!
!  B       (input/output) real array, dimension (LDB,NRHS)
!          On entry, the right hand side matrix B.
!          On exit, the solution matrix X.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!

  implicit none

  CHARACTER          TRANS
  INTEGER            INFO, LDA, LDB, N, NRHS
!     ..
!     .. Array Arguments ..
  INTEGER            IPIV( * )
  real ( kind = 8 )   A( LDA, * ), B( LDB, * )
  real ( kind = 8 )   ONE
  PARAMETER          ( ONE = 1.0D+00 )
!     ..
!     .. Local Scalars ..
  LOGICAL            NOTRAN
!     ..
!     .. External Functions ..
  LOGICAL            LSAME
  EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
  EXTERNAL           SLASWP, STRSM, XERBLA
!     ..
!     .. Intrinsic Functions ..
  INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
  INFO = 0
  NOTRAN = LSAME( TRANS, 'N' )
  IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT. &
         LSAME( TRANS, 'C' ) ) THEN
     INFO = -1
  ELSE IF( N.LT.0 ) THEN
     INFO = -2
  ELSE IF( NRHS.LT.0 ) THEN
     INFO = -3
  ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
     INFO = -5
  ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
     INFO = -8
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'SGETRS', -INFO )
     RETURN
  END IF
!
!     Quick return if possible
!
  IF( N.EQ.0 .OR. NRHS.EQ.0 ) RETURN
!
  IF( NOTRAN ) THEN
!
!        Solve A * X = B.
!
!        Apply row interchanges to the right hand sides.
!
     CALL SLASWP( NRHS, B, LDB, 1, N, IPIV, 1 )
!
!        Solve L*X = B, overwriting B with X.
!
     CALL STRSM( 'Left', 'Lower', 'No transpose', 'Unit', N, NRHS, &
                   ONE, A, LDA, B, LDB )
!
!        Solve U*X = B, overwriting B with X.
!
     CALL STRSM( 'Left', 'Upper', 'No transpose', 'Non-unit', N, &
                   NRHS, ONE, A, LDA, B, LDB )
  ELSE
!
!        Solve A' * X = B.
!
!        Solve U'*X = B, overwriting B with X.
!
     CALL STRSM( 'Left', 'Upper', 'Transpose', 'Non-unit', N, NRHS, &
                   ONE, A, LDA, B, LDB )
!
!        Solve L'*X = B, overwriting B with X.
!
     CALL STRSM( 'Left', 'Lower', 'Transpose', 'Unit', N, NRHS, ONE, &
                   A, LDA, B, LDB )
!
!        Apply row interchanges to the solution vectors.
!
     CALL SLASWP( NRHS, B, LDB, 1, N, IPIV, -1 )
  END IF

  RETURN
END
subroutine slaswp(n,a,lda,k1,k2,ipiv,incx)

!*****************************************************************************80
!
!! SLASWP performs a series of row interchanges on the matrix A.
!  One row interchange is initiated for each of rows K1 through K2 of A.
!
!  Arguments
!   == =======
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.
!
!  A       (input/output) real array, dimension (LDA,N)
!          On entry, the matrix of column dimension N to which the row
!          interchanges will be applied.
!          On exit, the permuted matrix.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.
!
!  K1      (input) INTEGER
!          The first element of IPIV for which a row interchange will
!          be done.
!
!  K2      (input) INTEGER
!          The last element of IPIV for which a row interchange will
!          be done.
!
!  IPIV    (input) INTEGER array, dimension (M*abs(INCX))
!          The vector of pivot indices.  Only the elements in positions
!          K1 through K2 of IPIV are accessed.
!          IPIV(K) = L implies rows K and L are to be interchanged.
!
!  INCX    (input) INTEGER
!          The increment between successive values of IPIV.  If IPIV
!          is negative, the pivots are applied in reverse order.
!
!
  implicit none

  integer ( kind = 4 )            incx, k1, k2, lda, n
!     ..
!     .. Array Arguments ..
  integer ( kind = 4 )            ipiv( * )
  real ( kind = 8 )   a( lda, * )
  integer ( kind = 4 )            i, ip, ix
!     ..
!     .. External Subroutines ..
  external           sswap
!     ..
!     .. Executable Statements ..
!
!     Interchange row I with row IPIV(I) for each of rows K1 through K2.
!
  if (  incx == 0 )return
  if (  incx>0 ) then
     ix = k1
  else
     ix = 1 + ( 1-k2 )*incx
  end if
  if (  incx == 1 ) then
     do 10 i = k1, k2
        ip = ipiv( i )
        if (  ip/= i )call sswap( n, a( i, 1 ), lda, a( ip, 1 ), lda )
   10    continue
  else if (  incx>1 ) then
     do 20 i = k1, k2
        ip = ipiv( ix )
        if (  ip/= i )call sswap( n, a( i, 1 ), lda, a( ip, 1 ), lda )
        ix = ix + incx
   20    continue
  else if (  incx<0 ) then
     do 30 i = k2, k1, -1
        ip = ipiv( ix )
        if (  ip/= i )call sswap( n, a( i, 1 ), lda, a( ip, 1 ), lda )
        ix = ix + incx
   30    continue
  end if
!
  return
end
subroutine solcon ( a, area, disjac,dpara3,dpdyn,dudyn,dvdyn,dydpn,eqn, &
  etan,etaq,flarea,g,gdif,gdifc,gold,gradf,gtar,ibc,idfd,ids,ierror, &
  ifds,igrid,ijac,indx,iopt,ipivot,ipred,ishapb,ishapf,ismooth,isotri,itype, &
  iwrite,jjac,maxnew,nelem,neqn,nlband,node,np,npar,nparb,nparf,npe,nprof,nrow, &
  numel,nx,ny,para,parjac,parnew,parold,phi,res,sens,splbmp,splflo,stpmax, &
  syseqn,taubmp,tauflo,tolnew,wateb,watep,wateu,watev,wquad,xbl,xbord,xbr,xc, &
  xquad,xsin,xsiq,ybl,ybord,ybr,yc,yquad)

!*****************************************************************************80
!
!! SOLCON computes the flow solution for a given set of parameters.
!
!  Discussion:
!
!    The routine is given a set of old parameter values in PAROLD,
!    and a corresponding flow solution G.
!
!    It is also given a set of parameter values in PARNEW,
!    and its goal is to compute the corresponding flow solution G.
!
!    It does this using continuation from the old solution.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) nelem
  integer ( kind = 4 ) neqn
  integer ( kind = 4 ) np
  integer ( kind = 4 ) npar
  integer ( kind = 4 ) nparb
  integer ( kind = 4 ) nparf
  integer ( kind = 4 ) npe
  integer ( kind = 4 ) nrow
  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny

  real ( kind = 8 ) a(nrow,neqn)
  real ( kind = 8 ) area(nelem,3)
  real ( kind = 8 ) cost
  real ( kind = 8 ) costb
  real ( kind = 8 ) costp
  real ( kind = 8 ) costu
  real ( kind = 8 ) costv
  real ( kind = 8 ) disjac
  real ( kind = 8 ) dist
  real ( kind = 8 ) dpara3(npar)
  real ( kind = 8 ) dpdyn(np)
  real ( kind = 8 ) dudyn(np)
  real ( kind = 8 ) dvdyn(np)
  real ( kind = 8 ) dydpn(np,nparb)
  character ( len = 2 ) eqn(neqn)
  real ( kind = 8 ) etan(6)
  real ( kind = 8 ) etaq(3)
  real ( kind = 8 ) flarea
  real ( kind = 8 ) g(neqn)
  real ( kind = 8 ) gdif(neqn,npar)
  real ( kind = 8 ) gdifc(neqn,npar)
  real ( kind = 8 ) gold(neqn)
  real ( kind = 8 ) gradf(neqn,npar)
  real ( kind = 8 ) gtar(neqn)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibc
  integer ( kind = 4 ) icon
  integer ( kind = 4 ) idfd
  integer ( kind = 4 ) ids
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ifds
  integer ( kind = 4 ) igrid
  integer ( kind = 4 ) ijac
  integer ( kind = 4 ) indx(np,3)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) iopt(npar)
  integer ( kind = 4 ) iopt2(15)
  integer ( kind = 4 ) ipivot(neqn)
  integer ( kind = 4 ) ipred
  integer ( kind = 4 ) ishapb
  integer ( kind = 4 ) ishapf
  integer ( kind = 4 ) ismooth
  integer ( kind = 4 ) isotri(nelem)
  integer ( kind = 4 ) itype
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) iwrite
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jjac
  logical lmat
  logical ltarg
  logical lzero
  integer ( kind = 4 ) maxnew
  integer ( kind = 4 ) ncon
  integer ( kind = 4 ) nlband
  integer ( kind = 4 ) node(nelem,npe)
  integer ( kind = 4 ) nprof(2*ny-1)
  integer ( kind = 4 ) numel(np)
  real ( kind = 8 ) para(npar)
  real ( kind = 8 ) parjac(npar)
  real ( kind = 8 ) parnew(npar)
  real ( kind = 8 ) parold(npar)
  real ( kind = 8 ) phi(nelem,3,6,10)
  real ( kind = 8 ) res(neqn)
  real ( kind = 8 ) sens(neqn,npar)
  real ( kind = 8 ) splbmp(4,nparb+2,0:nparb)
  real ( kind = 8 ) splflo(4,nparf+2,0:nparf)
  real ( kind = 8 ) stpmax
  real ( kind = 8 ) stpsiz
  character ( len = 20 ) syseqn
  real ( kind = 8 ) taubmp(nparb+2)
  real ( kind = 8 ) tauflo(nparf+2)
  real ( kind = 8 ) tolnew
  real ( kind = 8 ) wateb
  real ( kind = 8 ) watep
  real ( kind = 8 ) wateu
  real ( kind = 8 ) watev
  real ( kind = 8 ) wquad(3)
  real ( kind = 8 ) xbl
  real ( kind = 8 ) xbord(nx)
  real ( kind = 8 ) xbr
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) xquad(nelem,3)
  real ( kind = 8 ) xsin(6)
  real ( kind = 8 ) xsiq(3)
  real ( kind = 8 ) ybl
  real ( kind = 8 ) ybord(ny)
  real ( kind = 8 ) ybr
  real ( kind = 8 ) yc(np)
  real ( kind = 8 ) yquad(nelem,3)
!
  ival = 1
  call imemry ( 'inc', 'SolCon_calls', ival )
!
!  Refuse to run this problem if the new or old inverse viscosity
!  is not positive.
!
  if ( parnew(nparf+nparb+1) < 0.0D+00 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SOLCON - Fatal error!'
    write ( *, * ) '  PARNEW(NU_INV) = ', parnew(nparf+nparb+1)
    stop
  end if

  if ( parold(nparf+nparb+1) < 0.0D+00 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SolCon - Warning!'
    write ( *, * ) '  PAROLD(NU_INV) = ', parold(nparf+nparb+1)
  end if
!
!  Compute the L2 distance from PAROLD to PARNEW.
!
  dist = 0.0D+00
  do i = 1, nparf+nparb+1
    dist = dist+(parnew(i)-parold(i))**2
  end do
  dist = sqrt ( dist )
!
!  Compute the number of steps.
!
  lzero = .true.
  do i = 1, nparf
    if ( parnew(i) /= 0.0D+00 ) then
      lzero = .false.
    end if
  end do

  if ( lzero ) then
    ncon = 1
  else
    ncon = int(dist/stpmax)+1
    ncon = max(ncon,1)
  end if

  stpsiz = dist / real ( ncon, kind = 8 )
!
!  This is a clumsy attempt to handle situations where RE is not
!  generally a variable, but does need to be varied at startup,
!  when we move from a dummy solution at RE = 1 to a real solution
!  at RE = 100, for instance.  The critical thing is to ensure that
!  the appropriate RE sensitivities are computed.
!
  do i = 1,npar
    iopt2(i) = iopt(i)
  end do

  do i = 1,npar
    if ( iopt2(i) == 0 ) then
      if ( parnew(i)/= parold(i) ) then
        iopt2(i) = 1
      end if
    end if
  end do

10    continue

  g(1:neqn) = gold(1:neqn)
!
!  Begin continuation loop.
!
  do icon = 1,ncon

    ival = 1
    call imemry('inc','SolCon_steps',ival)
!
!  Compute PARA, the value of the parameters at the next point.
!
    do i = 1,npar
      para(i) = ((ncon-icon)*parold(i)+icon*parnew(i)) / real ( ncon, kind = 8 )
    end do
!
!  Estimate G, the next solution.
!
    if ( icon>0 ) then
!
!  Predict the value of the next point.
!
!  IPRED = 0, G is unchanged.
!  IPRED = 1 or 3, use the sensitivities
!  IPRED = 2 or 4, use the finite differences.
!
      if ( ipred == 0 ) then

      else if ( ipred == 1 ) then

        do i = 1,neqn
          do j = 1,npar
            g(i) = g(i)+sens(i,j)*(parnew(j)-parold(j)) / real ( ncon, kind = 8 )
          end do
        end do

      else if ( ipred == 2 ) then

        do i = 1,neqn
          do j = 1,npar
            g(i) = g(i)+gdif(i,j)*(parnew(j)-parold(j)) / real ( ncon, kind = 8 )
          end do
        end do

      else if ( ipred == 3 ) then

        do i = 1,neqn
          do j = 1,npar
            g(i) = g(i)+(gdif(i,j)-gradf(i,j))*(parnew(j)-parold(j))/ real ( ncon, kind = 8 )
          end do
        end do

      else if ( ipred == 4 ) then

        do i = 1,neqn
          do j = 1,npar
            g(i) = g(i)+(sens(i,j)+gradf(i,j))*(parnew(j)-parold(j))/ real ( ncon, kind = 8 )
          end do
        end do

      else

        write ( *, * ) ' '
        write ( *, * ) 'SOLCON - Fatal error!'
        write ( *, * ) '  Unknown value of IPRED = ',ipred
        stop

      end if

    end if

    ltarg = .false.
    call lmemry('get','target',ltarg)
!
!  Try to get the exact solution G for the current parameters PARA.
!
    ival = 1
    if ( icon<ncon ) then
      call imemry('inc','Points_Con',ival)
    else if ( ltarg ) then
      call imemry('inc','Points_Tar',ival)
    else if ( itype == 3 ) then
      call imemry('inc','Points_Opt',ival)
    else if ( itype == 7 ) then
      call imemry('inc','Points_Opt',ival)
    else
      call imemry('inc','Points_March',ival)
    end if

    lzero = .true.
    do i = 1,nparf
      if ( para(i)/= 0.0D+00 ) then
        lzero = .false.
      end if
    end do

    if ( lzero ) then

      g(1:neqn) = 0.0D+00

    else

      call flosol ( a,area,disjac,eqn,etan,etaq,flarea,g,ierror,igrid,ijac, &
        indx,ipivot,ishapb,ishapf,isotri,iwrite,jjac,maxnew,nelem,neqn, &
        nlband,node,np,npar,nparb,nparf,npe,nrow,nx,ny,para,parjac,phi,res, &
        splbmp,splflo,syseqn,taubmp,tauflo,tolnew,wquad,xbl,xbord,xbr,xc, &
        xquad,xsin,xsiq,ybl,ybord,ybr,yc,yquad)

      if ( ierror/= 0 ) then
        write ( *, * ) ' '
        write ( *, * ) 'SOLCON - Fatal error!'
        write ( *, * ) '  FLOSOL returned nonzero error flag.'
        write ( *, * ) '  Failure on step ',icon,' of ',ncon
        write ( *, * ) '  Stepsize was ',stpsiz
        write ( *, * ) ' '
        write ( *, * ) '  Parameters were:'
        write ( *, * ) ' '
        call pr_parameter(nparb,nparf,para)
        go to 20
      end if

    end if
!
!  Get the finite difference sensitivities and the direct cost
!  finite differences.
!
    lzero = .true.
    do i = 1,nparf
      if ( para(i)/= 0.0D+00 ) then
        lzero = .false.
      end if
    end do

    if ( lzero .or. idfd /= 0 .or. ifds /= 0 ) then

      if ( lzero .and. idfd == 0.and.ifds==0 ) then
        write ( *, * ) 'SOLCON - Calling GETGRD, zero inflow.'
      end if

      if ( idfd/= 0 ) then
        call get_cost(cost,costb,costp,costu,costv,g,gtar,indx,ishapb,neqn, &
          np,nparb,nprof,ny,splbmp,taubmp,wateb,watep,wateu,watev,xbl,xbr, &
          ybl,ybr,yc)
      end if

      call getgrd ( a,area,cost,disjac,dpara3,eqn,etan,etaq,flarea,g, &
        gdif,gtar,idfd,ierror,igrid,ijac,indx,iopt2,ipivot,ishapb, &
        ishapf,isotri,iwrite,jjac,maxnew,nelem,neqn,nlband,node,np,npar, &
        nparb,nparf,npe,nprof,nrow,nx,ny,para,parjac,phi,res,splbmp,splflo, &
        syseqn,taubmp,tauflo,tolnew,wateb,watep,wateu,watev,wquad,xbl,xbord, &
        xbr,xc,xquad,xsin,xsiq,ybl,ybord,ybr,yc,yquad)

      if ( ierror/= 0 ) then
        write ( *, * ) ' '
        write ( *, * ) 'SOLCON - Fatal error!'
        write ( *, * ) '  GETGRD returns IERROR = ',ierror
        return
      end if

    end if
!
!  Get dPdY, dUdY and dVdY
!
    call getdu ( dpdyn,dudyn,dvdyn,etan,g,indx,isotri,nelem,neqn,node,np, &
      npe,numel,xc,xsin,yc)
!
!  If requested, smooth the data.
!
    if ( ismooth == 3 ) then
      call getdu4 ( dudyn, np, nx, ny )
      call getdu4 ( dvdyn, np, nx, ny )
      call getdu4(dpdyn,np,nx,ny)
    else if ( ismooth == 4 ) then
      call getdu4(dudyn,np,nx,ny)
      call getdu4(dvdyn,np,nx,ny)
    end if
!
!  Get the finite difference correction GRADF
!
    call getfix(dpdyn,dudyn,dvdyn,dydpn,gradf,indx,iopt2,ishapb,neqn,np,npar, &
      nparb,nparf,splbmp,taubmp,xbl,xbr,xc,yc)
!
!  Compute GDIFC, the corrected finite difference estimate of the
!  sensitivity.
!
    do i = 1,neqn
      do j = 1,npar
        gdifc(i,j) = gdif(i,j)-gradf(i,j)
      end do
    end do
!
!  Get the sensitivities.
!
!  In the special case where the inflow is entirely zero, use the
!  finite differences instead.  (We suspect the actual jacobian may be
!  singular.)  Of course, if we didn't compute the finite differences,
!  we're in trouble!
!
    if ( ids/= 0 ) then

      if ( lzero ) then

        do i = 1,neqn
          do j = 1,npar
            sens(i,j) = gdif(i,j)
          end do
        end do

      else
!
!  If necessary, evaluate and factor the jacobian.
!
        lmat = .false.
        call lmemry('get','have_fp',lmat)
!
!  TEMPORARY
!  This seems to fix a problem with the computation of sensitivities.
!  I don't know why this problem is happening.
!  The value of PARJAC is equal to PARA, so that's not the problem.
!  Perhaps the value of UOLD, VOLD, POLD?
!
        lmat = .false.
!
        if ( .not.lmat ) then

          call fprime ( a,area,eqn,g,indx,nelem,neqn,nlband,node, &
            np,npar,nparb,nparf,npe,nrow,para,phi,syseqn)

          call sgbtrf(neqn,neqn,nlband,nlband,a,nrow,ipivot,info)

          ival = 1
          call imemry('inc','Factor_calls',ival)

          if ( info/= 0 ) then

            write ( *, * ) ' '
            write ( *, * ) 'SOLCON - Fatal error!'
            write ( *, * ) '  Jacobian factorization failed.'
            write ( *, * ) '  SGBTRF returns nonzero INFO = ',info
            ierror = 1
            return
          else
            lmat = .true.
            call lmemry('set','have_fp',lmat)
            if ( .not.ltarg ) then
              parjac(1:npar) = para(1:npar)
            end if
          end if

        end if
!
!  Solve the sensitivity equations.
!
        call getsen ( a,area,dudyn,dvdyn,eqn,g,ibc,indx,iopt2,ipivot,ishapb, &
          ishapf,nelem,neqn,nlband,node,np,npar,nparb,nparf,npe,nrow, &
          phi,sens,splbmp,splflo,taubmp,tauflo,xc,yc)

      end if

    end if

  end do
!
!  After successful continuation, update GOLD and PAROLD.
!
  gold(1:neqn) = g(1:neqn)
  para(1:npar) = parnew(1:npar)
  parold(1:npar) = parnew(1:npar)

  return
!
!  Continuation failed.  If we only used 1 step, try again.
!
20    continue

  if ( ncon == 1 ) then
    ncon = 5
    stpsiz = dist / real ( ncon, kind = 8 )
    write ( *, * ) ' '
    write ( *, * ) 'SOLCON - Note:'
    write ( *, * ) '  Retrying continuation taking ',ncon,' steps.'
    write ( *, * ) '  Distance is ',dist
    write ( *, * ) '  New stepsize is ',stpsiz
    go to 10
  else
    write ( *, * ) ' '
    write ( *, * ) 'SOLCON - Note:'
    write ( *, * ) '  Continuation was not retried, ncon = ',ncon
  end if

  return
end
subroutine sscal(n,da,dx,incx)

!*****************************************************************************80
!
!! SSCAL scales a vector by a constant.
!
  implicit none

  real ( kind = 8 ) da
  real ( kind = 8 ) dx(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  if ( n <= 0)return

  if ( incx == 1)go to 20

  if ( incx<0 ) then
    ix = (-n+1)*incx + 1
  else
    ix = 1
  end if

  do i = 1,n
    dx(ix) = da*dx(ix)
    ix = ix + incx
  end do

  return

   20 m = mod(n,5)

  do i = 1,m
    dx(i) = da*dx(i)
  end do

  do i = m+1,n,5
    dx(i) = da*dx(i)
    dx(i + 1) = da*dx(i + 1)
    dx(i + 2) = da*dx(i + 2)
    dx(i + 3) = da*dx(i + 3)
    dx(i + 4) = da*dx(i + 4)
  end do

  return
end
subroutine sswap(n,dx,incx,dy,incy)

!*****************************************************************************80
!
!! SSWAP interchanges two vectors.
!
  implicit none

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

  if ( n <= 0)return

  if ( incx == 1 .and. incy==1)go to 20

  ix = 1
  iy = 1
  if ( incx<0)ix = (-n+1)*incx + 1
  if ( incy<0)iy = (-n+1)*incy + 1

  do i = 1,n
    dtemp = dx(ix)
    dx(ix) = dy(iy)
    dy(iy) = dtemp
    ix = ix + incx
    iy = iy + incy
  end do

  return

   20 m = mod(n,3)

  do i = 1,m
    dtemp = dx(i)
    dx(i) = dy(i)
    dy(i) = dtemp
  end do

  do i = m+1,n,3
    dtemp = dx(i)
    dx(i) = dy(i)
    dy(i) = dtemp
    dtemp = dx(i + 1)
    dx(i + 1) = dy(i + 1)
    dy(i + 1) = dtemp
    dtemp = dx(i + 2)
    dx(i + 2) = dy(i + 2)
    dy(i + 2) = dtemp
  end do

  return
end
subroutine stbsv(uplo,trans,diag,n,k,a,lda,x,incx)

!*****************************************************************************80
!
!! STBSV  solves one of the systems of equations
!
!     A*x = b,   or   A'*x = b,
!
!  where b and x are n element vectors and A is an n by n unit, or
!  non-unit, upper or lower triangular band matrix, with ( k + 1 )
!  diagonals.
!
!  No test for singularity or near-singularity is included in this
!  routine. Such tests must be performed before calling this routine.
!
!  Parameters
!   == ========
!
!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the matrix is an upper or
!           lower triangular matrix as follows:
!
!              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!
!              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!
!           Unchanged on exit.
!
!  TRANS  - CHARACTER*1.
!           On entry, TRANS specifies the equations to be solved as
!           follows:
!
!              TRANS = 'N' or 'n'   A*x = b.
!
!              TRANS = 'T' or 't'   A'*x = b.
!
!              TRANS = 'C' or 'c'   A'*x = b.
!
!           Unchanged on exit.
!
!  DIAG   - CHARACTER*1.
!           On entry, DIAG specifies whether or not A is unit
!           triangular as follows:
!
!              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!
!              DIAG = 'N' or 'n'   A is not assumed to be unit
!                                  triangular.
!
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the order of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  K      - INTEGER.
!           On entry with UPLO = 'U' or 'u', K specifies the number of
!           super-diagonals of the matrix A.
!           On entry with UPLO = 'L' or 'l', K specifies the number of
!           sub-diagonals of the matrix A.
!           K must satisfy  0  <=  K.
!           Unchanged on exit.
!
!  A      - real array of DIMENSION ( LDA, n ).
!           Before entry with UPLO = 'U' or 'u', the leading ( k + 1 )
!           by n part of the array A must contain the upper triangular
!           band part of the matrix of coefficients, supplied column by
!           column, with the leading diagonal of the matrix in row
!           ( k + 1 ) of the array, the first super-diagonal starting at
!           position 2 in row k, and so on. The top left k by k triangle
!           of the array A is not referenced.
!           The following program segment will transfer an upper
!           triangular band matrix from conventional full matrix storage
!           to band storage:
!
!                 DO 20, J = 1, N
!                    M = K + 1 - J
!                    DO 10, I = MAX( 1, J - K ), J
!                       A( M + I, J ) = matrix( I, J )
!              10    CONTINUE
!              20 CONTINUE
!
!           Before entry with UPLO = 'L' or 'l', the leading ( k + 1 )
!           by n part of the array A must contain the lower triangular
!           band part of the matrix of coefficients, supplied column by
!           column, with the leading diagonal of the matrix in row 1 of
!           the array, the first sub-diagonal starting at position 1 in
!           row 2, and so on. The bottom right k by k triangle of the
!           array A is not referenced.
!           The following program segment will transfer a lower
!           triangular band matrix from conventional full matrix storage
!           to band storage:
!
!                 DO 20, J = 1, N
!                    M = 1 - J
!                    DO 10, I = J, MIN( N, J + K )
!                       A( M + I, J ) = matrix( I, J )
!              10    CONTINUE
!              20 CONTINUE
!
!           Note that when DIAG = 'U' or 'u' the elements of the array A
!           corresponding to the diagonal elements of the matrix are not
!           referenced, but are assumed to be unity.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           ( k + 1 ).
!           Unchanged on exit.
!
!  X      - real array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the n
!           element right-hand side vector b. On exit, X is overwritten
!           with the solution vector x.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!
!  Level 2 Blas routine.
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
!
  implicit none

  integer ( kind = 4 ) incx, k, lda, n
  character        diag, trans, uplo
  real ( kind = 8 )   a( lda, * ), x( * )

  real ( kind = 8 )   zero
  parameter        ( zero = 0.0D+00 )
!     .. Local Scalars ..
  real ( kind = 8 )   temp
  integer ( kind = 4 )            i, info, ix, j, jx, kplus1, kx, l
  logical            nounit
!     .. External Functions ..
  logical            lsame
  external           lsame
!     .. External Subroutines ..
  external           xerbla
!
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
  info = 0
  if     ( .not.lsame( uplo , 'u' ) .and.  &
             .not.lsame( uplo , 'l' )       ) then
     info = 1
  else if (  .not.lsame( trans, 'n' ) .and.  &
              .not.lsame( trans, 't' ) .and.  &
              .not.lsame( trans, 'c' )       ) then
     info = 2
  else if (  .not.lsame( diag , 'u' ) .and.  &
             .not.lsame( diag , 'n' )       ) then
     info = 3
  else if (  n<0  ) then
     info = 4
  else if (  k<0  ) then
     info = 5
  else if (  lda<( k + 1 )  ) then
     info = 7
  else if (  incx == 0  ) then
     info = 9
  end if
  if (  info/= 0  ) then
     call xerbla( 'stbsv ', info )
     return
  end if
!
!     Quick return if possible.
!
  if (  n == 0 )return
!
  nounit = lsame( diag, 'n' )
!
!     Set up the start point in X if the increment is not unity. This
!     will be  ( N - 1 )*INCX  too small for descending loops.
!
  if (  incx <= 0  ) then
     kx = 1 - ( n - 1 )*incx
  else if (  incx/= 1  ) then
     kx = 1
  end if
!
!     Start the operations. In this version the elements of A are
!     accessed by sequentially with one pass through A.
!
  if (  lsame( trans, 'n' )  ) then
!
!        Form  x : =  inv( A )*x.
!
     if (  lsame( uplo, 'u' )  ) then
        kplus1 = k + 1
        if (  incx == 1  ) then
           do 20, j = n, 1, -1
              if (  x( j )/= zero  ) then
                 l = kplus1 - j
                 if (  nounit )x( j ) = x( j )/a( kplus1, j )
                 temp = x( j )
                 do 10, i = j - 1, max( 1, j - k ), -1
                    x( i ) = x( i ) - temp*a( l + i, j )
   10                continue
              end if
   20          continue
        else
           kx = kx + ( n - 1 )*incx
           jx = kx
           do 40, j = n, 1, -1
              kx = kx - incx
              if (  x( jx )/= zero  ) then
                 ix = kx
                 l  = kplus1 - j
                 if (  nounit )x( jx ) = x( jx )/a( kplus1, j )
                 temp = x( jx )
                 do 30, i = j - 1, max( 1, j - k ), -1
                    x( ix ) = x( ix ) - temp*a( l + i, j )
                    ix      = ix      - incx
   30                continue
              end if
              jx = jx - incx
   40          continue
        end if
     else
        if (  incx == 1  ) then
           do 60, j = 1, n
              if (  x( j )/= zero  ) then
                 l = 1 - j
                 if (  nounit )x( j ) = x( j )/a( 1, j )
                 temp = x( j )
                 do 50, i = j + 1, min( n, j + k )
                    x( i ) = x( i ) - temp*a( l + i, j )
   50                continue
              end if
   60          continue
        else
           jx = kx
           do 80, j = 1, n
              kx = kx + incx
              if (  x( jx )/= zero  ) then
                 ix = kx
                 l  = 1  - j
                 if (  nounit )x( jx ) = x( jx )/a( 1, j )
                 temp = x( jx )
                 do 70, i = j + 1, min( n, j + k )
                    x( ix ) = x( ix ) - temp*a( l + i, j )
                    ix      = ix      + incx
   70                continue
              end if
              jx = jx + incx
   80          continue
        end if
     end if
  else
!
!        Form  x : =  inv( A')*x.
!
     if (  lsame( uplo, 'u' )  ) then
        kplus1 = k + 1
        if (  incx == 1  ) then
           do 100, j = 1, n
              temp = x( j )
              l    = kplus1 - j
              do 90, i = max( 1, j - k ), j - 1
                 temp = temp - a( l + i, j )*x( i )
   90             continue
              if (  nounit )temp = temp/a( kplus1, j )
              x( j ) = temp
  100          continue
        else
           jx = kx
           do 120, j = 1, n
              temp = x( jx )
              ix   = kx
              l    = kplus1  - j
              do 110, i = max( 1, j - k ), j - 1
                 temp = temp - a( l + i, j )*x( ix )
                 ix   = ix   + incx
  110             continue
              if (  nounit ) temp = temp/a( kplus1, j )
              x( jx ) = temp
              jx      = jx   + incx
              if (  j>k )kx = kx + incx
  120          continue
        end if
     else
        if (  incx == 1  ) then
           do 140, j = n, 1, -1
              temp = x( j )
              l    = 1      - j
              do 130, i = min( n, j + k ), j + 1, -1
                 temp = temp - a( l + i, j )*x( i )
  130             continue
              if (  nounit )temp = temp/a( 1, j )
              x( j ) = temp
  140          continue
        else
           kx = kx + ( n - 1 )*incx
           jx = kx
           do 160, j = n, 1, -1
              temp = x( jx )
              ix   = kx
              l    = 1       - j
              do 150, i = min( n, j + k ), j + 1, -1
                 temp = temp - a( l + i, j )*x( ix )
                 ix   = ix   - incx
  150             continue
              if (  nounit )temp = temp/a( 1, j )
              x( jx ) = temp
              jx      = jx   - incx
              if (  ( n - j ) >= k )kx = kx - incx
  160          continue
        end if
     end if
  end if

  return
end
subroutine strsm(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)

!*****************************************************************************80
!
!! STRSM  solves one of the matrix equations
!
!     op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
!
!  where alpha is a scalar, X and B are m by n matrices, A is a unit, or
!  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
!
!     op( A ) = A   or   op( A ) = A'.
!
!  The matrix X is overwritten on B.
!
!  Parameters
!   == ========
!
!  SIDE   - CHARACTER*1.
!           On entry, SIDE specifies whether op( A ) appears on the left
!           or right of X as follows:
!
!              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
!
!              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
!
!           Unchanged on exit.
!
!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the matrix A is an upper or
!           lower triangular matrix as follows:
!
!              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!
!              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!
!           Unchanged on exit.
!
!  TRANSA - CHARACTER*1.
!           On entry, TRANSA specifies the form of op( A ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSA = 'N' or 'n'   op( A ) = A.
!
!              TRANSA = 'T' or 't'   op( A ) = A'.
!
!              TRANSA = 'C' or 'c'   op( A ) = A'.
!
!           Unchanged on exit.
!
!  DIAG   - CHARACTER*1.
!           On entry, DIAG specifies whether or not A is unit triangular
!           as follows:
!
!              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!
!              DIAG = 'N' or 'n'   A is not assumed to be unit
!                                  triangular.
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry, M specifies the number of rows of B. M must be at
!           least zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of B.  N must be
!           at least zero.
!           Unchanged on exit.
!
!  ALPHA  - real.
!           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
!           zero then  A is not referenced and  B need not be set before
!           entry.
!           Unchanged on exit.
!
!  A      - real array of DIMENSION ( LDA, k ), where k is m
!           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
!           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
!           upper triangular part of the array  A must contain the upper
!           triangular matrix  and the strictly lower triangular part of
!           A is not referenced.
!           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
!           lower triangular part of the array  A must contain the lower
!           triangular matrix  and the strictly upper triangular part of
!           A is not referenced.
!           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
!           A  are not referenced either,  but are assumed to be  unity.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
!           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
!           then LDA must be at least max( 1, n ).
!           Unchanged on exit.
!
!  B      - real array of DIMENSION ( LDB, n ).
!           Before entry,  the leading  m by n part of the array  B must
!           contain  the  right-hand  side  matrix  B,  and  on exit  is
!           overwritten by the solution matrix  X.
!
!  LDB    - INTEGER.
!           On entry, LDB specifies the first dimension of B as declared
!           in  the  calling  (sub)  program.   LDB  must  be  at  least
!           max( 1, m ).
!           Unchanged on exit.
!
!
!  Level 3 Blas routine.
!
!
!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
!
!
  implicit none

  character        side, uplo, transa, diag
  integer ( kind = 4 )            m, n, lda, ldb
  real ( kind = 8 )   alpha
  real ( kind = 8 )   a( lda, * ), b( ldb, * )
  logical            lsame
  external           lsame
!     .. External Subroutines ..
  external           xerbla
!
!     .. Local Scalars ..
  logical            lside, nounit, upper
  integer ( kind = 4 )            i, info, j, k, nrowa
  real ( kind = 8 )   temp
!     .. Parameters ..
  real ( kind = 8 )   one         , zero
  parameter        ( one = 1.0D+00, zero = 0.0D+00 )
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
  lside  = lsame( side  , 'l' )
  if (  lside  ) then
     nrowa = m
  else
     nrowa = n
  end if
  nounit = lsame( diag  , 'n' )
  upper  = lsame( uplo  , 'u' )
!
  info   = 0
  if (       ( .not.lside                ) .and.  &
             ( .not.lsame( side  , 'r' ) )       ) then
     info = 1
  else if (  ( .not.upper                ) .and.  &
              ( .not.lsame( uplo  , 'l' ) )       ) then
     info = 2
  else if (  ( .not.lsame( transa, 'n' ) ) .and.  &
              ( .not.lsame( transa, 't' ) ) .and.  &
              ( .not.lsame( transa, 'c' ) )       ) then
     info = 3
  else if (  ( .not.lsame( diag  , 'u' ) ) .and.  &
              ( .not.lsame( diag  , 'n' ) )       ) then
     info = 4
  else if (  m  <0                ) then
     info = 5
  else if (  n  <0                ) then
     info = 6
  else if (  lda<max( 1, nrowa )  ) then
     info = 9
  else if (  ldb<max( 1, m     )  ) then
     info = 11
  end if
  if (  info/= 0  ) then
     call xerbla( 'strsm ', info )
     return
  end if
!
!     Quick return if possible.
!
  if (  n == 0 )return
!
!     And when  alpha == zero.
!
  if (  alpha == zero  ) then
     do 20, j = 1, n
        do 10, i = 1, m
           b( i, j ) = zero
   10       continue
   20    continue
     return
  end if
!
!     Start the operations.
!
  if (  lside  ) then
     if (  lsame( transa, 'n' )  ) then
!
!           Form  B : =  alpha*inv( A )*B.
!
        if (  upper  ) then
           do 60, j = 1, n
              if (  alpha/= one  ) then
                 do 30, i = 1, m
                    b( i, j ) = alpha*b( i, j )
   30                continue
              end if
              do 50, k = m, 1, -1
                 if (  b( k, j )/= zero  ) then
                    if (  nounit )b( k, j ) = b( k, j )/a( k, k )
                    do 40, i = 1, k - 1
                       b( i, j ) = b( i, j ) - b( k, j )*a( i, k )
   40                   continue
                 end if
   50             continue
   60          continue
        else
           do 100, j = 1, n
              if (  alpha/= one  ) then
                 do 70, i = 1, m
                    b( i, j ) = alpha*b( i, j )
   70                continue
              end if
              do 90 k = 1, m
                 if (  b( k, j )/= zero  ) then
                    if (  nounit ) b( k, j ) = b( k, j )/a( k, k )
                    do 80, i = k + 1, m
                       b( i, j ) = b( i, j ) - b( k, j )*a( i, k )
   80                   continue
                 end if
   90             continue
  100          continue
        end if
     else
!
!           Form  B : =  alpha*inv( A' )*B.
!
        if (  upper  ) then
           do 130, j = 1, n
              do 120, i = 1, m
                 temp = alpha*b( i, j )
                 do 110, k = 1, i - 1
                    temp = temp - a( k, i )*b( k, j )
  110                continue
                 if (  nounit )temp = temp/a( i, i )
                 b( i, j ) = temp
  120             continue
  130          continue
        else
           do 160, j = 1, n
              do 150, i = m, 1, -1
                 temp = alpha*b( i, j )
                 do 140, k = i + 1, m
                    temp = temp - a( k, i )*b( k, j )
  140                continue
                 if (  nounit )temp = temp/a( i, i )
                 b( i, j ) = temp
  150             continue
  160          continue
        end if
     end if
  else
     if (  lsame( transa, 'n' )  ) then
!
!           Form  B : =  alpha*B*inv( A ).
!
        if (  upper  ) then
           do 210, j = 1, n
              if (  alpha/= one  ) then
                 do 170, i = 1, m
                    b( i, j ) = alpha*b( i, j )
  170                continue
              end if
              do 190, k = 1, j - 1
                 if (  a( k, j )/= zero  ) then
                    do 180, i = 1, m
                       b( i, j ) = b( i, j ) - a( k, j )*b( i, k )
  180                   continue
                 end if
  190             continue
              if (  nounit  ) then
                 temp = one/a( j, j )
                 do 200, i = 1, m
                    b( i, j ) = temp*b( i, j )
  200                continue
              end if
  210          continue
        else
           do 260, j = n, 1, -1
              if (  alpha/= one  ) then
                 do 220, i = 1, m
                    b( i, j ) = alpha*b( i, j )
  220                continue
              end if
              do 240, k = j + 1, n
                 if (  a( k, j )/= zero  ) then
                    do 230, i = 1, m
                       b( i, j ) = b( i, j ) - a( k, j )*b( i, k )
  230                   continue
                 end if
  240             continue
              if (  nounit  ) then
                 temp = one/a( j, j )
                 do 250, i = 1, m
                   b( i, j ) = temp*b( i, j )
  250                continue
              end if
  260          continue
        end if
     else
!
!           Form  B : =  alpha*B*inv( A' ).
!
        if (  upper  ) then
           do 310, k = n, 1, -1
              if (  nounit  ) then
                 temp = one/a( k, k )
                 do 270, i = 1, m
                    b( i, k ) = temp*b( i, k )
  270                continue
              end if
              do 290, j = 1, k - 1
                 if (  a( j, k )/= zero  ) then
                    temp = a( j, k )
                    do 280, i = 1, m
                       b( i, j ) = b( i, j ) - temp*b( i, k )
  280                   continue
                 end if
  290             continue
              if (  alpha/= one  ) then
                 do 300, i = 1, m
                    b( i, k ) = alpha*b( i, k )
  300                continue
              end if
  310          continue
        else
           do 360, k = 1, n
              if (  nounit  ) then
                 temp = one/a( k, k )
                 do 320, i = 1, m
                    b( i, k ) = temp*b( i, k )
  320                continue
              end if
              do 340, j = k + 1, n
                 if (  a( j, k )/= zero  ) then
                    temp = a( j, k )
                    do 330, i = 1, m
                       b( i, j ) = b( i, j ) - temp*b( i, k )
  330                   continue
                 end if
  340             continue
              if (  alpha/= one  ) then
                 do 350, i = 1, m
                    b( i, k ) = alpha*b( i, k )
  350                continue
              end if
  360          continue
        end if
     end if
  end if
!
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
  integer ( kind = 4 ) d
  character ( len = 8 ) date
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  character ( len = 10 )  time
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y
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
subroutine trans ( det, detadx, detady, dxsidx, dxsidy, eta, ielem, nelem, &
  node, np, npe, xc, xsi, yc )

!*****************************************************************************80
!
!! TRANS calculates the mapping from reference to physical elements.
!
!  Discussion:
!
!    The routine determines the biquadratic transformation which maps the
!    reference element in (XSI,ETA) space into a particular physical
!    isoparametric element in (X,Y) space.
!
!    We know everything about the isoparametric element once we
!    specify the location of its six nodes.
!
!    TRANS computes the entries of the jacobian of the transformation
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
!  Refence Element Diagram:
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
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) DET, the determinant of the jacobian of the transformation 
!    between the reference and isoparametric elements.
!
!    Output, real ( kind = 8 ) DETADX, DETADY, the partial
!    derivative d ETA/d X and d ETA/d Y at (XSI,ETA).
!
!    Output, real ( kind = 8 ) DXSIDX, DXSIDY, the partial
!    derivative d XSI/d X and d XSI/d Y at (XSI,ETA).
!
!    Input, real ( kind = 8 ) ETA, the ETA coordinate of the point.
!
!    Input, integer ( kind = 4 ) IELEM, the number of the isoparametric
!    element we are examining.
!
!    Input, integer ( kind = 4 ) NODE(MAXELM,6), contains the numbers
!    of the nodes that make up each element.  Element number
!    I is associated with nodes NODE(I,1) through NODE(I,6).
!
!    Input, integer ( kind = 4 ) NP, the number of nodes.
!
!    Input, real ( kind = 8 ) XC(NP), the X coordinates of the nodes.
!
!    Input, real ( kind = 8 ) XSI, the XSI coordinate of the point.
!
!    Input, real ( kind = 8 ) YC(NP), the Y coordinates of the nodes.
!
  implicit none

  integer ( kind = 4 ) nelem
  integer ( kind = 4 ) np
  integer ( kind = 4 ) npe

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
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ielem
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) node(nelem,npe)
  real ( kind = 8 ) x
  real ( kind = 8 ) xn(6)
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) xsi
  real ( kind = 8 ) y
  real ( kind = 8 ) yn(6)
  real ( kind = 8 ) yc(np)
!
!  Pick off the X, Y coordinates of the nodes and store them
!  in two short lists.
!
  do i = 1, 6
    xn(i) = xc(node(ielem,i))
    yn(i) = yc(node(ielem,i))
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
  a1 =  2.0D+00*xn(1)+2.0D+00*xn(3)-4.0D+00*xn(6)
  b1 = -4.0D+00*xn(3)-4.0D+00*xn(4)+4.0D+00*xn(5)+4.0D+00*xn(6)
  c1 =  2.0D+00*xn(2)+2.0D+00*xn(3)-4.0D+00*xn(5)
  d1 = -3.0D+00*xn(1)      -xn(3)+4.0D+00*xn(6)
  e1 =       -xn(2)      +xn(3)+4.0D+00*xn(4)-4.0D+00*xn(6)
  f1 =        xn(1)

  a2 =  2.0D+00*yn(1)+2.0D+00*yn(3)-4.0D+00*yn(6)
  b2 = -4.0D+00*yn(3)-4.0D+00*yn(4)+4.0D+00*yn(5)+4.0*yn(6)
  c2 =  2.0D+00*yn(2)+2.0D+00*yn(3)-4.0D+00*yn(5)
  d2 = -3.0D+00*yn(1)      -yn(3)+4.0D+00*yn(6)
  e2 =       -yn(2)      +yn(3)+4.0D+00*yn(4)-4.0D+00*yn(6)
  f2 =        yn(1)
!
!  Compute the partial derivatives at the point (XSI,ETA).
!  This is the jacobian matrix
!
!    J: (XSI,ETA) --> (X,Y).
!
  dxdxsi =  2.0D+00*a1*xsi +       b1*eta + d1
  dxdeta =        b1*xsi + 2.0D+00*c1*eta + e1

  dydxsi =  2.0D+00*a2*xsi +       b2*eta + d2
  dydeta =        b2*xsi + 2.0D+00*c2*eta + e2
!
!  Compute the determinant of the jacobian matrix:
!
!    J: (XSI,ETA) --> (X,Y)
!
  det = dxdxsi*dydeta-dxdeta*dydxsi
!
!  Watch out for a zero determinant.
!
  if ( det == 0.0D+00 ) then
    write ( *, * ) ' '
    write ( *, * ) 'TRANS - Fatal error!'
    write ( *, * ) '  The jacobian J: (XSI,ETA) --> (X,Y) is singular!'
    write ( *, * ) '  This occurred for element number ',ielem
    write ( *, * ) '  Local coordinates XSI,ETA = ',xsi,eta
    x = a1*xsi**2+b1*xsi*eta+c1*eta**2+d1*xsi+e1*eta+f1
    y = a2*xsi**2+b2*xsi*eta+c2*eta**2+d2*xsi+e2*eta+f2
    write ( *, * ) '  Global coordinates X,Y = ',x,y
    write ( *, * ) ' '
    write ( *, * ) '  The X, Y nodes were:'
    write ( *, * ) ' '
    do i = 1,6
      write(*,*)xn(i),yn(i)
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
  dxsidx =  dydeta/det
  dxsidy = -dxdeta/det

  detadx = -dydxsi/det
  detady =  dxdxsi/det

  return
end
subroutine upvalq ( dppdx,dppdy,dupdx,dupdy,dvpdx,dvpdy,gp,ielem,indx,iquad, &
  nelem,neqn,node,np,npe,phi,pp,up,vp)

!*****************************************************************************80
!
!! UPVALQ evaluates sensitivities at a quadrature point in a given element.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IELEM, the element in which the point lies
!    at which the quantities are desired.
!
  implicit none

  integer ( kind = 4 ) nelem
  integer ( kind = 4 ) neqn
  integer ( kind = 4 ) np
  integer ( kind = 4 ) npe

  real ( kind = 8 ) coef
  real ( kind = 8 ) dppdx
  real ( kind = 8 ) dppdy
  real ( kind = 8 ) dqdx
  real ( kind = 8 ) dqdy
  real ( kind = 8 ) dupdx
  real ( kind = 8 ) dupdy
  real ( kind = 8 ) dvpdx
  real ( kind = 8 ) dvpdy
  real ( kind = 8 ) dwdx
  real ( kind = 8 ) dwdy
  real ( kind = 8 ) gp(neqn)
  integer ( kind = 4 ) ielem
  integer ( kind = 4 ) indx(np,3)
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) iq
  integer ( kind = 4 ) iquad
  integer ( kind = 4 ) iun
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) node(nelem,npe)
  real ( kind = 8 ) phi(nelem,3,6,10)
  real ( kind = 8 ) pp
  real ( kind = 8 ) q
  real ( kind = 8 ) up
  real ( kind = 8 ) vp
  real ( kind = 8 ) w

  ival = 1
  call imemry('inc','UpValQ_calls',ival)
!
!  Start all the functions at zero.
!
  pp = 0.0D+00
  up = 0.0D+00
  vp = 0.0D+00
  dppdx = 0.0D+00
  dppdy = 0.0D+00
  dupdx = 0.0D+00
  dupdy = 0.0D+00
  dvpdx = 0.0D+00
  dvpdy = 0.0D+00
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
  do iq = 1,6

    w = phi(ielem,iquad,iq,1)
    dwdx = phi(ielem,iquad,iq,2)
    dwdy = phi(ielem,iquad,iq,3)

    q = phi(ielem,iquad,iq,4)
    dqdx = phi(ielem,iquad,iq,5)
    dqdy = phi(ielem,iquad,iq,6)
!
!  Now that we have the basis function values, we need to look
!  up the coefficient COEF that multiplies the basis function.
!
    ip = node(ielem,iq)

    iun = indx(ip,1)
    coef = gp(iun)
    up = up+coef*w
    dupdx = dupdx+coef*dwdx
    dupdy = dupdy+coef*dwdy

    iun = indx(ip,2)
    coef = gp(iun)
    vp = vp+coef*w
    dvpdx = dvpdx+coef*dwdx
    dvpdy = dvpdy+coef*dwdy

    iun = indx(ip,3)
    if ( iun>0 ) then
      coef = gp(iun)
      pp = pp+coef*q
      dppdx = dppdx+coef*dqdx
      dppdy = dppdy+coef*dqdy
    end if

  end do

  return
end
subroutine uval ( detadx,detady,dpdx,dpdy,dudx,dudy,dvdx,dvdy,dxsidx,dxsidy, &
  eta,g,ielem,indx,isotri,nelem,neqn,node,np,npe,p,u,v,xc,xq,xsi,yc,yq)

!*****************************************************************************80
!
!! UVAL evaluates the velocities and pressure at any point in a given element.
!
!  Discussion:
!
!    UVAL also computes the X and Y derivatives of the quantities.
!
!    If the element is not isoparametric, then UVAL requires the
!    physical X and Y coordinates of the point.
!
!    If the element is isoparametric, UVAL requires the XSI, ETA
!    coordinates of the point.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) DETADX, DETADY, the partial derivative
!    d ETA/d X and d ETA/d Y at (XSI,ETA).
!
!    Input, real ( kind = 8 ) DXSIDX, DXSIDY, the partial derivative
!    d XSI/d X and d XSI/d Y at (XSI,ETA).
!
!    Input, real ( kind = 8 ) ETA, the ETA coordinate of the point,
!    needed only if the element is isoparametric.
!
!    Input, real ( kind = 8 ) G(NEQN), the computed solution vector, in which are stored
!    pressures and velocities.
!
!    Input, integer ( kind = 4 ) IELEM, the element in which the point lies
!    at which the quantities are desired.
!
!    Input, real ( kind = 8 ) XQ, the X coordinate of the point.
!
!    Input, real ( kind = 8 ) XSI, the XSI coordinate of the point,
!    needed only if the element is isoparametric.
!
!    Input, real ( kind = 8 ) YQ, the Y coordinate of the point.
!
  implicit none

  integer ( kind = 4 ) nelem
  integer ( kind = 4 ) neqn
  integer ( kind = 4 ) np
  integer ( kind = 4 ) npe

  real ( kind = 8 ) coef
  real ( kind = 8 ) detadx
  real ( kind = 8 ) detady
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
  real ( kind = 8 ) dxsidx
  real ( kind = 8 ) dxsidy
  real ( kind = 8 ) eta
  real ( kind = 8 ) g(neqn)
  integer ( kind = 4 ) ielem
  integer ( kind = 4 ) indx(np,3)
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) iq
  integer ( kind = 4 ) isotri(nelem)
  integer ( kind = 4 ) iun
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) node(nelem,npe)
  real ( kind = 8 ) p
  real ( kind = 8 ) q
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) w
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) xq
  real ( kind = 8 ) xsi
  real ( kind = 8 ) yc(np)
  real ( kind = 8 ) yq

  ival = 1
  call imemry('inc','UVal_calls',ival)

  p = 0.0D+00
  u = 0.0D+00
  v = 0.0D+00
  dpdx = 0.0D+00
  dpdy = 0.0D+00
  dudx = 0.0D+00
  dudy = 0.0D+00
  dvdx = 0.0D+00
  dvdy = 0.0D+00

  do iq = 1,6
!
!  Evaluate the basis functions W and Q, and their derivatives
!  DQDX, DQDY, DWDX, DWDY at XQ, YQ.
!
    if ( isotri(ielem) == 0 .or. isotri(ielem)==1 ) then

      call qbf ( ielem,iq,w,dwdx,dwdy,nelem,node,np,npe,xc,xq,yc,yq)

      call bsp ( q,dqdx,dqdy,ielem,iq,nelem,node,np,npe,xc,xq,yc,yq)

    else

      call refqbf ( w,dwdx,dwdy,detadx,detady,dxsidx,dxsidy,eta,iq,xsi)

      call refbsp ( q,dqdx,dqdy,detadx,detady,iq,dxsidx,dxsidy,eta,xsi)

    end if
!
!  Compute the coefficients at the node at XP, YP.
!
    ip = node(ielem,iq)

    iun = indx(ip,1)
    coef = g(iun)
    u = u+coef*w
    dudx = dudx+coef*dwdx
    dudy = dudy+coef*dwdy

    iun = indx(ip,2)
    coef = g(iun)
    v = v+coef*w
    dvdx = dvdx+coef*dwdx
    dvdy = dvdy+coef*dwdy

    iun = indx(ip,3)
    if ( iun>0 ) then
      coef = g(iun)
      p = p+coef*q
      dpdx = dpdx+coef*dqdx
      dpdy = dpdy+coef*dqdy
    end if

  end do

  return
end
subroutine uvalq ( dpdx,dpdy,dudx,dudy,dvdx,dvdy,g,ielem,indx,iquad,nelem, &
  neqn,node,np,npe,p,phi,u,v)

!*****************************************************************************80
!
!! UVALQ evaluates the velocities and pressure at a quadrature point.
!
!  Discusion:
!
!    The X and Y derivatives are also evaluated.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 January 2001
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
!    Output, real ( kind = 8 ) P, the value of the pressure.
!
!    Output, real ( kind = 8 ) U, the value of the horizontal
!    velocity.
!
!    Output, real ( kind = 8 ) V, the value of the vertical
!    velocity.
!
  implicit none

  integer ( kind = 4 ) nelem
  integer ( kind = 4 ) neqn
  integer ( kind = 4 ) np
  integer ( kind = 4 ) npe

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
  real ( kind = 8 ) g(neqn)
  integer ( kind = 4 ) ielem
  integer ( kind = 4 ) indx(np,3)
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) iq
  integer ( kind = 4 ) iquad
  integer ( kind = 4 ) iun
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) node(nelem,npe)
  real ( kind = 8 ) p
  real ( kind = 8 ) phi(nelem,3,6,10)
  real ( kind = 8 ) q
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) w

  ival = 1
  call imemry('inc','UValQ_calls',ival)
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
  do iq = 1,6

    w = phi(ielem,iquad,iq,1)
    dwdx = phi(ielem,iquad,iq,2)
    dwdy = phi(ielem,iquad,iq,3)

    q = phi(ielem,iquad,iq,4)
    dqdx = phi(ielem,iquad,iq,5)
    dqdy = phi(ielem,iquad,iq,6)
!
!  Now that we have the basis function values, we need to look
!  up the coefficient COEF that multiplies the basis function.
!
    ip = node(ielem,iq)

    iun = indx(ip,1)
    coef = g(iun)
    u = u+coef*w
    dudx = dudx+coef*dwdx
    dudy = dudy+coef*dwdy

    iun = indx(ip,2)
    coef = g(iun)
    v = v+coef*w
    dvdx = dvdx+coef*dwdx
    dvdy = dvdy+coef*dwdy

    iun = indx(ip,3)
    if ( iun>0 ) then
      coef = g(iun)
      p = p+coef*q
      dpdx = dpdx+coef*dqdx
      dpdy = dpdy+coef*dqdy
    end if

  end do

  return
end
subroutine xerbla(srname,info)

!*****************************************************************************80
!
!! XERBLA  is an error handler for the LAPACK routines.
!
!  Discussion:
!
!    It is called by an LAPACK routine if an input parameter has an
!    invalid value.  A message is printed and execution stops.
!
!  Parameters:
!
!    Input, character*6, SRNAME, the routine which called XERBLA.
!
!    Input, integer ( kind = 4 ) INFO, the position of the invalid parameter in 
!    the parameter list of the calling routine.
!
  implicit none

  character ( len = 6 ) srname
  integer ( kind = 4 ) info

  write( *, fmt = 9999 )srname, info

  stop

 9999 format( ' ** on entry to ', a6, ' parameter number ', i2, ' had ', &
          'an illegal value' )
end
subroutine xofxsi ( eta,ielem,nelem,node,np,npe,x,xc,xsi,y,yc)

!*****************************************************************************80
!
!! XOFXSI maps reference coordinates to physical coordinates.
!
!  Discussion:
!
!    Here is a graph of the (XSI, ETA) reference triangle we will use.
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
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ETA, the ETA coordinate of the point.
!
!    Input, integer ( kind = 4 ) IELEM, the number of the isoparametric
!    element we are examining.
!
!    Input, integer ( kind = 4 ) MAXELM, the maximum number of elements.
!
!    Input, integer ( kind = 4 ) NODE(MAXELM,6), contains the numbers
!    of the nodes that make up each element.  Element number
!    I is associated with nodes NODE(I,1) through NODE(I,6).
!
!    Input, integer ( kind = 4 ) NP, the number of nodes.
!
!    Output, real ( kind = 8 ) X, the X coordinate of the point.
!
!    Input, real ( kind = 8 ) XC(NP), the X coordinates of the
!    nodes.
!
!    Input, real ( kind = 8 ) XSI, the XSI coordinate of the point.
!
!    Output, real ( kind = 8 ) Y, the Y coordinate of the point.
!
!    Input, real ( kind = 8 ) YC(NP), the Y coordinates of the
!    nodes.
!
  implicit none

  integer ( kind = 4 ) nelem
  integer ( kind = 4 ) np
  integer ( kind = 4 ) npe

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
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ielem
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) node(nelem,npe)
  real ( kind = 8 ) x
  real ( kind = 8 ) xn(6)
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) xsi
  real ( kind = 8 ) y
  real ( kind = 8 ) yn(6)
  real ( kind = 8 ) yc(np)
!
  ival = 1
  call imemry('inc','Xofxsi_calls',ival)
!
!  Pick off the X, Y coordinates of the nodes and store them
!  in two short lists.
!
  do i = 1,6
    xn(i) = xc(node(ielem,i))
    yn(i) = yc(node(ielem,i))
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
  a1 =  2.0D+00*xn(1)+2.0D+00*xn(3)-4.0D+00*xn(6)
  b1 = -4.0D+00*xn(3)-4.0D+00*xn(4)+4.0D+00*xn(5)+4.0D+00*xn(6)
  c1 =  2.0D+00*xn(2)+2.0D+00*xn(3)-4.0D+00*xn(5)
  d1 = -3.0D+00*xn(1)      -xn(3)+4.0D+00*xn(6)
  e1 =       -xn(2)      +xn(3)+4.0D+00*xn(4)-4.0D+00*xn(6)
  f1 =        xn(1)

  a2 =  2.0D+00*yn(1)+2.0D+00*yn(3)-4.0D+00*yn(6)
  b2 = -4.0D+00*yn(3)-4.0D+00*yn(4)+4.0D+00*yn(5)+4.0D+00*yn(6)
  c2 =  2.0D+00*yn(2)+2.0D+00*yn(3)-4.0D+00*yn(5)
  d2 = -3.0D+00*yn(1)      -yn(3)+4.0D+00*yn(6)
  e2 =       -yn(2)      +yn(3)+4.0D+00*yn(4)-4.0D+00*yn(6)
  f2 =        yn(1)

  x = a1*xsi**2 + b1*xsi*eta + c1*eta**2 + d1*xsi + e1*eta + f1

  y = a2*xsi**2 + b2*xsi*eta + c2*eta**2 + d2*xsi + e2*eta + f2

  return
end
subroutine xy_print ( maxnp, np, xc, yc )

!*****************************************************************************80
!
!! XY_PRINT prints the X and Y coordinates of the nodes.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) MAXNP, the maximum number of nodes.
!
!    Input, integer ( kind = 4 ) NP, the number of nodes.
!
!    Input, real ( kind = 8 ) XC(MAXNP), the X coordinates of the nodes.
!
!    Input, real ( kind = 8 ) YC(MAXNP), the Y coordinates of the nodes.
!
  implicit none

  integer ( kind = 4 ) maxnp

  integer ( kind = 4 ) i
  integer ( kind = 4 ) np
  real ( kind = 8 ) xc(maxnp)
  real ( kind = 8 ) xold
  real ( kind = 8 ) yc(maxnp)
  real ( kind = 8 ) yold

  write ( *, * ) ' '
  write ( *, * ) 'XY_PRINT'
  write ( *, * ) ' '
  write ( *, * ) '  I, XC(I), YC(I)'
  write ( *, * ) ' '

  xold = 0.0D+00
  yold = 0.0D+00

  do i = 1, np

    if ( xold /= xc(i) .and. yold /= yc(i) ) then
      write ( *, * ) ' '
    end if

    write ( *, '(i8,2g14.6)' ) i, xc(i), yc(i)

    xold = xc(i)
    yold = yc(i)

  end do

  return
end
subroutine xy_set ( igrid, ishapb, np, nparb, nx, ny, splbmp, taubmp, xbl, &
  xbord, xbr, xc, ybl, ybord, ybr, yc )

!*****************************************************************************80
!
!! XY_SET sets the X and Y coordinates of the nodes.
!
!  Discussion:
!
!    This routine assumes that the nodes are numbered
!    in "stacks", starting with the least X and Y coordinates,
!    then fixing X and running through all values of Y, then
!    increasing X to the next value and running through all
!    values of Y, and so on.  For example:
!
!      5  10  15
!      4   9  14
!      3   8  13
!      2   7  12
!      1   6  11
!
!    This allows us to treat the vectors XC and YC as two dimensional
!    arrays, instead of the one dimensional vectors they are declared
!    to be.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) np
  integer ( kind = 4 ) nparb
  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny

  integer ( kind = 4 ) i
  integer ( kind = 4 ) igrid
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) ishapb
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jderiv
  logical lval
  real ( kind = 8 ) splbmp(4,nparb+2,0:nparb)
  real ( kind = 8 ) taubmp(nparb+2)
  real ( kind = 8 ) x
  real ( kind = 8 ) xbl
  real ( kind = 8 ) xbord(nx)
  real ( kind = 8 ) xbr
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) y
  real ( kind = 8 ) ybl
  real ( kind = 8 ) ybot
  real ( kind = 8 ) ybord(ny)
  real ( kind = 8 ) ybr
  real ( kind = 8 ) yc(np)
  real ( kind = 8 ) yvec(nparb+2)
  real ( kind = 8 ) yy
!
  ival = 1
  call imemry ( 'inc', 'XY_SET_calls', ival )

  yvec(1:nparb+2) = splbmp(1,1:nparb+2,0)

  ip = 0

  do i = 1, 2*nx-1
    do j = 1, 2*ny-1

      ip = ip+1

      if ( igrid == 0 ) then
        x = real ( ( i - 1 ) ) * 10.0D+00 / real ( 2 * nx - 2 )
      else if ( mod(i,2) == 1 ) then
        x = xbord((i+1)/2)
      else
        x = 0.5D+00 * ( xbord(i/2) + xbord(1+i/2) )
      end if

      if ( abs(x-xbl)*(2*nx-2) <= 0.5D+00 ) then
        x = xbl
      else if ( abs(x-xbr)*(2*nx-2) <= 0.5D+00 ) then
        x = xbr
      end if

      if ( x <= xbl ) then
        ybot = ybl
      else if ( xbl <= x .and. x <= xbr ) then
        if ( ishapb == 1 ) then
          call plval ( nparb+2, x, taubmp, ybot, yvec )
        else if ( ishapb == 2 ) then
          call pqval ( nparb+2, x, taubmp, ybot, yvec )
        else if ( ishapb == 3 ) then
          jderiv = 0
          call ppvalu ( taubmp, splbmp(1,1,0), nparb+1, 4, x, jderiv, ybot )
        end if
      else
        ybot = ybr
      end if

      if ( igrid == 0 ) then
        y = ( real ( 2 * ny - j - 1 ) * ybot &
            + real (          j - 1, kind = 8 ) * 3.0D+00 ) &
            / real ( 2 * ny     - 2, kind = 8 )
      else 
        if ( mod(j,2) == 1 ) then
          yy = ybord((j+1)/2)
        else
          yy = 0.5D+00*(ybord(j/2)+ybord(1+j/2))
        end if
        y = ( ( 3.0D+00 - yy ) * ybot + yy * 3.0D+00 ) / 3.0D+00
      end if

      xc(ip) = x
      yc(ip) = y

    end do
  end do
!
!  Note that the internal XY coordinates have been set.
!
  lval = .false.
  call lmemry ( 'set', 'need_xy', lval )

  return
end
