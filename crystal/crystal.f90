program main

!*****************************************************************************80
!
!! MAIN is the main program for CRYSTAL.
!
!  Discussion:
!
!    CRYSTAL carries out a simulation of the formation of a silicon 
!    crystal.
!
!  Modified:
!
!    22 March 2002
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: liv = 60
  integer, parameter :: npar = 1

  integer, parameter :: lv = 77+npar*(npar+17)/2

  real ( kind = 8 ) cost
  external cryfun
  real ( kind = 8 ) d(npar)
  external dummy
  integer i
  integer ido
  integer iopt
  integer ipar(1)
  integer iv(liv)
  integer nfun
  real ( kind = 8 ) par(npar)
  real ( kind = 8 ) rpar(1)
  real ( kind = 8 ) tarray(2)
  real ( kind = 8 ) temp
  real ( kind = 8 ) v(lv)

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CRYSTAL:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Version of 25 January 1996'

  iopt = 0

  par(1) = 1713.0D+00
!
!  If IOPT = 0, just call CRYFUN once.
!
  if ( iopt == 0 ) then

    call cryfun ( npar, par, nfun, cost, ipar, rpar, dummy )
!
!  If IOPT = 1, then carry out an optimization.
!
  else

    d(1:npar) = 1.0D+00
!
!  Set the 611 work arrays to default values.
!
    ido = 2
    call deflt ( ido, iv, liv, lv, v )
!
!  Tell SMSNO that we have called DEFLT already.
!
    iv(1) = 12
!
!  Set the maximum number of iterations.
!
    iv(18) = 40

    call smsno ( npar, d, par, cryfun, iv, liv, lv, v, ipar, rpar, dummy )

  end if


  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CRYSTAL'
  write ( *, '(a)' ) '  Normal end of execution.'

  stop
end
subroutine dummy

!*****************************************************************************80
!
!! DUMMY is a dummy subroutine needed as formal input to CRYFUN.
!
!  Modified:
!
!    22 March 2002
!
!  Author:
!
!    John Burkardt
!
  implicit none

  return
end
subroutine cryfun ( npar, par, nfun, cost, ipar, rpar, dummy )

!*****************************************************************************80
!
!! CRYFUN evaluates the cost function for the minimizing software.
!
!
  implicit none

  integer, parameter :: maxbot = 20
  integer, parameter :: ni = 64
  integer, parameter :: nj = 64
  integer, parameter :: nk = 14
  integer npar
  integer, parameter :: ns = 10

  real ( kind = 8 ) ae1(ni,nj)
  real ( kind = 8 ) ae2(ni,nj)
  real ( kind = 8 ) ak1(ni,nj)
  real ( kind = 8 ) ak2(ni,nj)
  real ( kind = 8 ) area(ni,nj)
  real ( kind = 8 ) areal
  real ( kind = 8 ) areas
  real ( kind = 8 ) areat
  real ( kind = 8 ) b
  real ( kind = 8 ) b1jbl(nj)
  real ( kind = 8 ) b2jbl(nj)
  real ( kind = 8 ) birad
  real ( kind = 8 ) bo
  real ( kind = 8 ) cappa
  real ( kind = 8 ) cd
  real ( kind = 8 ) ce1
  real ( kind = 8 ) ce2
  real ( kind = 8 ) cfo
  real ( kind = 8 ) cinc
  real ( kind = 8 ) cinco
  real ( kind = 8 ) cmu
  real ( kind = 8 ) cost
  real ( kind = 8 ) cvn
  real ( kind = 8 ) delt
  real ( kind = 8 ) dtm
  external dummy
  real ( kind = 8 ) epsad
  real ( kind = 8 ) epsil
  real ( kind = 8 ) ewall
  real ( kind = 8 ) f(ni,nj,ns)
  real ( kind = 8 ) fcsl
  real ( kind = 8 ) fjeta(ni,nj)
  real ( kind = 8 ) fjksi(ni,nj)
  real ( kind = 8 ) fks
  real ( kind = 8 ) fksl
  real ( kind = 8 ) fma
  real ( kind = 8 ) fmax(ns)
  real ( kind = 8 ) fn(ni,nj,ns)
  real ( kind = 8 ) fnu
  real ( kind = 8 ) fo(ni,nj,ns)
  real ( kind = 8 ) fr
  real ( kind = 8 ) frsl
  real ( kind = 8 ) gamt(ni,nj)
  real ( kind = 8 ) grash
  real ( kind = 8 ) hamag
  real ( kind = 8 ) heta(ni,nj)
  real ( kind = 8 ) hf
  real ( kind = 8 ) hksi(ni,nj)
  integer i
  integer icost
  integer icrys
  integer inturb
  integer ipar(1)
  integer iplot
  integer ipref
  integer iprint
  integer iter
  integer izone
  integer j
  integer jcrys
  integer jpref
  integer k
  integer l0
  integer l1
  integer last
  integer lastt
  logical lblk(ns)
  logical lconv
  logical lortho
  logical lsolve(ns)
  integer m0
  integer m1
  integer mode
  integer nbot
  integer, save :: ncall = 0
  integer ndt
  integer nf
  integer nfun
  integer np
  integer npc
  integer nsolve(ns)
  integer ntimes(ns)
  real ( kind = 8 ) orth
  real ( kind = 8 ) par(npar)
  real ( kind = 8 ) pr
  real ( kind = 8 ) r(ni,nj)
  real ( kind = 8 ) ra
  real ( kind = 8 ) rdtm
  real ( kind = 8 ) re
  real ( kind = 8 ) recb
  real ( kind = 8 ) rect
  real ( kind = 8 ) relax(nk)
  real ( kind = 8 ) res(ns)
  real ( kind = 8 ) rho(ni,nj)
  real ( kind = 8 ) rhocon
  real ( kind = 8 ) rpar(1)
  real ( kind = 8 ) rpr
  real ( kind = 8 ) rueta(ni,nj)
  real ( kind = 8 ) ruksi(ni,nj)
  real ( kind = 8 ) sige
  real ( kind = 8 ) sigk
  real ( kind = 8 ) sigma
  real ( kind = 8 ) sigt
  real ( kind = 8 ) smax
  real ( kind = 8 ) smooth
  real ( kind = 8 ) ssum
  real ( kind = 8 ) stel
  real ( kind = 8 ) stes
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) t3
  real ( kind = 8 ) tal
  real ( kind = 8 ) tanca
  real ( kind = 8 ) tanca2
  real ( kind = 8 ) tas
  real ( kind = 8 ) tend
  real ( kind = 8 ) tf
  real ( kind = 8 ) tinit
  character ( len = 25 )  title(ns)
  real ( kind = 8 ) tnow
  real ( kind = 8 ) tw
  real ( kind = 8 ) vave
  real ( kind = 8 ) vol(ni,nj)
  real ( kind = 8 ) x(ni,nj)
  real ( kind = 8 ) xbot(maxbot)
  real ( kind = 8 ) xc(ni,nj)
  real ( kind = 8 ) xlen
  real ( kind = 8 ) y(ni,nj)
  real ( kind = 8 ) ybot(maxbot)
  real ( kind = 8 ) yc(ni,nj)
  real ( kind = 8 ) ylen

  ncall = ncall+1
  write ( *, * ) 'CRYFUN: PAR(1) = ',par(1)
!
!  Initialize the data.
!
  call inidat(ae1,ae2,ak1,ak2,area,b,birad,bo,cappa,cd,ce1, &
    ce2,cfo,cinc,cmu,cost,cvn,delt,dtm,epsad,epsil,ewall,f,fcsl,fks, &
    fksl,fma,fn,fnu,fo,fr,frsl,gamt,grash,hamag,heta,hf, &
    hksi,icost,icrys,inturb,iplot,ipref,iprint,jcrys, &
    jpref,l0,l1,last,lastt,lblk,lortho,lsolve,m0,m1,maxbot,mode, &
    nbot,ndt,ni,nj,nk,np,npar,npc,ns,nsolve,ntimes,orth, &
    par,pr,ra,rdtm,re,recb,rect,relax,rho,rhocon,rpr,rueta,ruksi,sige,sigk, &
    sigma,sigt,smooth,stel,stes,tal,tanca,tanca2,tas,tend,tf, &
    tinit,title,tnow,tw,vave,vol,x,xbot,xc,xlen,y,ybot,yc,ylen)
!
!  Print out data.
!
  if ( ncall == 1 ) then
    call prdat(b,birad,bo,cappa,cfo,cvn,delt,dtm,fcsl,fksl,fma, &
      fnu,fr,frsl,grash,icost,icrys,inturb,iprint, &
      jcrys,l0,last,lastt,m0,mode,nbot,ns,nsolve,ntimes,orth, &
      pr,ra,rdtm,recb,rect,rhocon,smooth,stel,stes,tanca,tanca2, &
      tend,tf,tinit,title,tw,vave,xlen,ylen)
  end if
!
!  Generate the initial grid XC, YC.
!
  call inigrd ( icrys, jcrys, l0, m0, nbot, xbot, xc, xlen, ybot, yc, ylen )
!
!  Adapt the grid.
!
  do k = 1, 10
    call adapt ( cvn, epsad, icrys, iprint, jcrys, l1, m1, nbot, &
      orth, smooth, xbot, xc, ybot, yc )
  end do
!
!  Once XC and YC are determined, compute the control volume areas.
!
  call doarea ( area, areal, areas, areat, icrys, jcrys, l0, m0, &
    ni, nj, xc, yc )
!
!  Each time iteration begins at this point.
!
10    continue
!
!  Begin the iterative solution of the state equations in the
!  solid zone.
!
  izone = 1
  l1 = l0
  m1 = jcrys
!
!  Only the temperature (variable 5) needs to be solved for.
!
  lsolve(1) = .false.
  lsolve(2) = .false.
  lsolve(3) = .false.
  lsolve(4) = .false.
  lsolve(5) = .true.
  lsolve(6) = .false.
  lsolve(7) = .false.
  lsolve(8) = .false.
  lsolve(9) = .false.
  lsolve(10) = .false.

  do iter = 0, last

    lconv = .true.
!
!  Set the X and Y coordinates of the primary nodes from XC, YC,
!  the coordinates of the corner nodes.
!
!  Note that this is only done for the left portion of the region!
!
    call setx ( l1, m1, ni, nj, x, xc, y, yc )
!
!  Compute various geometric quantities required for computing
!  derivatives.
!
    call setgeo(ae1,ae2,ak1,ak2,heta,hksi,l1,m1,mode,ni,nj, &
      r,vol,x,xc,y,yc)
!
!  Estimate pressure and momentum at control volume interfaces.
!
    if ( ndt == 0 ) then
      if ( iter == 0 ) then
        call initl(fjeta,fjksi,heta,hksi,l1,m1,ni,nj,f(1,1,3),rho, &
          rueta,ruksi,f(1,1,1),f(1,1,2),x,y)
      end if
    end if

    call setup(ae1,ae2,ak1,ak2,b1jbl, &
      b2jbl,birad,cappa,cd,ce1,ce2,cmu,epsil,ewall,f, &
      fcsl,fjeta,fjksi,fksl,fma,fmax,fn,fo,frsl,gamt,grash,hamag, &
      heta,hksi,inturb,icrys,iter,izone,jcrys,l1,lblk,lconv,lortho, &
      lsolve,m1,mode,nf,np,npc,nsolve,ntimes,pr,r,rdtm,re,recb, &
      rect,res,relax,rho,rhocon,rpr,rueta,ruksi,sige,sigk,sigt, &
      smax,ssum,stel,stes,tal,tas,tf,tw,vol,x,xc,y,yc)

    if ( iprint > 0 ) then
      call output(cfo,iter,izone,res,smax,ssum,f(1,1,5), &
        tnow,f(1,1,1),f(1,1,6))
    end if

    if ( lconv .and. iter >= 5 ) then
      go to 20
    end if

  end do

  if ( iprint > 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'CRYFUN - Warning!'
    write ( *, * ) '  The solid iteration has not converged'
    write ( *, * ) '  after ',iter,' iterations.'
    write ( *, * ) ' '
    write ( *, * ) '  Solid norm and max relative change:'
    write ( *, * ) ' '
    do i = 1, ns
      if ( lsolve(i) ) then
        fmax(i) = maxval ( f(1:l0,1:m0,i) )
        write(*,'(i3,2x,a25,2g14.6)') i, title(i), fmax(i), res(i)
      end if
    end do
  end if
!
!  Zone 2 (LIQUID) calculation
!
20    continue

  f(1:icrys,1:jcrys,5) = fo(1:icrys,1:jcrys,5)
  f(1:icrys,1:jcrys,9) = fo(1:icrys,1:jcrys,9)

  izone = 2
  l1 = icrys
  m1 = m0

  rueta(1:icrys,2) = 0.0D+00
  rueta(1:icrys,m0) = 0.0D+00
  rueta(1,1:m0) = 0.0D+00
  rueta(icrys,1:m0) = 0.0D+00
  ruksi(1:icrys,1) = 0.0D+00
  ruksi(1:icrys,m0) = 0.0D+00
  ruksi(2,1:m0) = 0.0D+00
  ruksi(icrys,1:m0) = 0.0D+00

  if ( inturb == 0 ) then
    lsolve(1) = .true.
    lsolve(2) = .true.
    lsolve(3) = .true.
    lsolve(4) = .true.
    lsolve(5) = .true.
    lsolve(6) = .false.
    lsolve(7) = .false.
    lsolve(8) = .false.
    lsolve(9) = .false.
    lsolve(10) = .false.
  else
    lsolve(1) = .true.
    lsolve(2) = .true.
    lsolve(3) = .true.
    lsolve(4) = .true.
    lsolve(5) = .true.
    lsolve(6) = .false.
    lsolve(7) = .true.
    lsolve(8) = .true.
    lsolve(9) = .false.
    lsolve(10) = .false.
  end if

  do iter = 1, last

    lconv = .true.
!
!  Set the turbulent viscosity.
!
    if ( inturb == 0 ) then
      do i = 2, l1-1
        do j = 2, m1-1
          gamt(i,j) = 0.0D+00
        end do
      end do
    else if ( inturb == 1 ) then
      do i = 2, l1-1
        do j = 2, m1-1
          gamt(i,j) = ( 1.0D+00 - relax(12) ) * gamt(i,j) &
            +relax(12)*cmu*rho(i,j)*f(i,j,7)**2/f(i,j,8)
        end do
      end do
    end if
!
!  Set the X, Y values.
!
    call setx ( l1, m1, ni, nj, x, xc, y, yc )
!
!  Compute various geometric quantities required for computing
!  derivatives.
!
    call setgeo(ae1,ae2,ak1,ak2,heta,hksi,l1,m1,mode,ni,nj, &
      r,vol,x,xc,y,yc)
!
!  Set the right hand side of certain flux boundary conditions.
!
    call setup(ae1,ae2,ak1,ak2,b1jbl, &
         b2jbl,birad,cappa,cd,ce1,ce2,cmu,epsil,ewall,f, &
         fcsl,fjeta,fjksi,fksl,fma,fmax,fn,fo,frsl,gamt,grash,hamag, &
         heta,hksi,inturb,icrys,iter,izone,jcrys,l1,lblk,lconv, &
         lortho, &
         lsolve,m1,mode,nf,np,npc,nsolve,ntimes,pr,r,rdtm,re,recb, &
         rect,res,relax,rho,rhocon,rpr,rueta,ruksi,sige,sigk,sigt, &
         smax,ssum,stel,stes,tal,tas,tf,tw,vol,x,xc,y,yc)
!
!  Compute the stream function PSI.
!
    f(2,2,10) = 0.0D+00

    do i = 2, icrys

      if ( i > 2 ) then
        f(i,2,10) = f(i-1,2,10) - rueta(i-1,2) * ae1(i-1,2) &
          - 0.5D+00 * ( ruksi(i-1,1) + ruksi(i,1) ) * ae2(i-1,2)
      end if

      do j = 3, m1-1
        t1 = (rueta(i,j-1)+rueta(i,j))
        t2 = (rueta(i-1,j-1)+rueta(i-1,j))
        t3 = 0.5D+00*(hksi(i,j-1)*t2+hksi(i-1,j-1)*t1)/(hksi(i,j-1) &
          +hksi(i-1,j-1))
        f(i,j,10) = f(i,j-1,10)+ruksi(i,j-1)*ak1(i,j-1)-t3*ak2(i,j-1)
      end do

    end do

    f(1,1,10) = f(2,2,10)

    do j = 2, m0
      f(1,j,10) = f(2,j,10)
    end do

    do i = 2, l0
      f(i,1,10) = f(i,2,10)
    end do
!
!  Optional printout.
!
    if ( 0 < iprint ) then
      call output(cfo,iter,izone,res,smax,ssum,f(1,1,5), &
        tnow,f(1,1,1),f(1,1,6))
    end if

    if ( lconv ) go to 30

  end do

  if ( iprint > 0 ) then

    write ( *, * ) ' '
    write ( *, * ) 'CRYFUN - Warning!'
    write ( *, * ) '  The liquid iteration has not converged'
    write ( *, * ) '  after ',iter,' iterations.'
    write ( *, * ) ' '
    write ( *, * ) '  Liquid norms and max relative change:'
    write ( *, * ) ' '
    do i = 1, ns
      if ( lsolve(i) ) then
        fmax(i) = maxval ( f(1:l0,1:m0,i) )
        write(*,'(i3,2x,a25,2g14.6)') i, title(i), fmax(i), res(i)
      end if
    end do

  end if
!
!  We have computed the state variables for the current time.
!
30    continue
!
!  Evaluate the cost function integrand at the current time.
!
  cinco = cinc
  call setcst ( area, cinc, icost, icrys, m0, ni, nj, f(1,1,5), tf, &
    f(1,1,1), f(1,1,2), vave )
!
!  Update the total cost, by adding the estimated contribution
!  from the current time interval.
!
!  TEMPORARY
!
  if ( ndt > 0 ) then
!   cost = cost + 0.5 * delt *( cinco + cinc )
    cost = cost + 0.5 * dtm * ( cinco + cinc )
  end if
!
!  If we've reached the end time, then write out final data,
!  possibly save a restart file, and stop.
!
  if ( ndt >= lastt ) then

    if ( iprint > 0 ) then
      call pmod(ipref,jcrys,jpref,l1,m1,f(1,1,3))
    end if

    cost = cost / ( sqrt ( areal ) * ( tend - tinit ) )
!
!  If requested, write out plot data.
!  Sadly, it is necessary to call SETX to generate the X, Y arrays
!  for the entire region.  Previous calls only generate them for
!  a portion of the region.
!
    if ( iplot == 1 ) then

      call setx ( l0, m0, ni, nj, x, xc, y, yc )

      call rswrit ( cost, f(1,1,9), gamt, icrys, jcrys, l0, m0, nbot, &
        f(1,1,3), f(1,1,4), f(1,1,10), &
        rueta, ruksi, f(1,1,5), f(1,1,8), f(1,1,7), tnow, f(1,1,1), &
        f(1,1,2), f(1,1,6), x, xbot, xc, y, yc, ybot )

    end if

    write ( *, * ) 'CRYFUN: COST =        ',cost

    return

  end if
!
!  Update the number of steps, and the current time.
!
  ndt = ndt+1
  tnow = tinit+ndt*dtm
!
!  Save a copy of the current data.
!
  do i = 1, l0
    do j = 1, m0
      do k = 1, ns
        fo(i,j,k) = f(i,j,k)
      end do
    end do
  end do

  call movgrd(b1jbl,b2jbl,bo,delt,fksl,fr,frsl,icrys,iprint,jcrys, &
    l0,m0,f(1,1,3),pr,re,stel,tanca,tanca2,xc,xlen,yc)

  l1 = l0
  m1 = m0
  call adapt(cvn,epsad,icrys,iprint,jcrys,l1,m1,nbot, &
    orth,smooth,xbot,xc,ybot,yc)
!
!  Once XC and YC are determined, compute the control volume areas.
!
  call doarea ( area, areal, areas, areat, icrys, jcrys, l0, m0, &
    ni, nj, xc, yc )

  go to 10
end
subroutine adapt ( cvn, epsad, icrys, iprint, jcrys, l1, m1, nbot, &
  orth, smooth, xbot, xc, ybot, yc )

!*****************************************************************************80
!
!! ADAPT executes MAGG, the multizone adaptive grid generation.
!
  implicit none

  integer nbot
  integer, parameter :: ni = 64
  integer, parameter :: nj = 64

  real ( kind = 8 ) a1
  real ( kind = 8 ) a11
  real ( kind = 8 ) a12
  real ( kind = 8 ) a2
  real ( kind = 8 ) a21
  real ( kind = 8 ) a22
  real ( kind = 8 ) a3
  real ( kind = 8 ) b1
  real ( kind = 8 ) b2
  real ( kind = 8 ) b3
  real ( kind = 8 ) c1
  real ( kind = 8 ) c2
  real ( kind = 8 ) c3
  real ( kind = 8 ) cvn
  real ( kind = 8 ) det
  real ( kind = 8 ) dot
  real ( kind = 8 ) dx
  real ( kind = 8 ) dxc(ni,nj)
  real ( kind = 8 ) dxmax
  real ( kind = 8 ) dy
  real ( kind = 8 ) dyc(ni,nj)
  real ( kind = 8 ) epsad
  real ( kind = 8 ) gn
  real ( kind = 8 ) gt
  real ( kind = 8 ) gx
  real ( kind = 8 ) gy
  integer i
  integer icrys
  integer iprint
  integer j
  integer jcrys
  integer l1
  integer m1
  integer nin
  real ( kind = 8 ) orth
  real ( kind = 8 ) pb(ni,nj)
  real ( kind = 8 ) pc1
  real ( kind = 8 ) pc2
  real ( kind = 8 ) pd(ni,nj)
  real ( kind = 8 ) ratio
  real ( kind = 8 ) res1
  real ( kind = 8 ) res2
  real ( kind = 8 ) rmax
  real ( kind = 8 ) s1
  real ( kind = 8 ) s2
  real ( kind = 8 ) scale
  real ( kind = 8 ) smooth
  real ( kind = 8 ) ssmot
  real ( kind = 8 ) ssweg
  real ( kind = 8 ) vmag2
  real ( kind = 8 ) wn
  real ( kind = 8 ) wt
  real ( kind = 8 ) wx
  real ( kind = 8 ) wy
  real ( kind = 8 ) xbot(nbot)
  real ( kind = 8 ) xc(ni,nj)
  real ( kind = 8 ) xdisp
  real ( kind = 8 ) xij
  real ( kind = 8 ) xjac
  real ( kind = 8 ) xn
  real ( kind = 8 ) xnn
  real ( kind = 8 ) xpp
  real ( kind = 8 ) xt
  real ( kind = 8 ) xtn
  real ( kind = 8 ) xtt
  real ( kind = 8 ) ybot(nbot)
  real ( kind = 8 ) yc(ni,nj)
  real ( kind = 8 ) ydisp
  real ( kind = 8 ) yij
  real ( kind = 8 ) yn
  real ( kind = 8 ) ynn
  real ( kind = 8 ) ypp
  real ( kind = 8 ) yt
  real ( kind = 8 ) ytn
  real ( kind = 8 ) ytt

  dxmax = 0.0D+00

  dxc(1:l1,1:m1) = 0.0D+00
  dyc(1:l1,1:m1) = 0.0D+00

  do i = 2, l1
    do j = 2, m1

      if ( i < icrys ) then
        pc1 = 0.95D+00 * ((i-(2+icrys)/2)/10.0D+00 )**2 + 0.05D+00
      else
        pc1 = 0.8D+00 * ((i-l1)/12.0)**2 + 0.2D+00
      end if

      if ( j <= jcrys ) then
        pc2 = 0.9D+00 * ((j-(2+jcrys)/2)/10.0D+00)**2 + 0.05D+00
      else
        pc2 = 0.9D+00 * ((j-(m1+jcrys)/2)/10.0D+00)**2 + 0.05D+00
      end if

      pb(i,j) = pc1 * pc2

      pd(i,j) = pc1 * pc2

    end do
  end do
!
!  This DO loop iteration is an iterative solution of a linear system.
!
  do nin = 1, 100

    ssmot = 0.0D+00
    ssweg = 0.0D+00
    rmax = 0.0D+00

    if ( iprint > 0 ) then

      if ( nin == 1 ) then
        write ( *, * ) ' '
        write ( *, * ) 'ADAPT: Iter   DXMAX     XC(20,20) YC(20,20)' // &
          ' XC(3,M1)  YC(3,M1)  RMAX'
      end if

      if ( mod(nin,50) == 0 ) then
        write(*,'(6x,i5,6g10.3)') &
          nin,dxmax,xc(20,20),yc(20,20),xc(3,m1),yc(3,m1),rmax
      end if
    end if

    do j = 3, m1-1
      do i = 3, l1-1

        xt = 0.5D+00 * ( xc(i+1,j) - xc(i-1,j) )
        xn = 0.5D+00 * ( xc(i,j+1) - xc(i,j-1) )
        yt = 0.5D+00 * ( yc(i+1,j) - yc(i-1,j) )
        yn = 0.5D+00 * ( yc(i,j+1) - yc(i,j-1) )

        xjac = xt*yn-xn*yt
!
!  Refer to table 2.1, page 14, of Hui Zhang's thesis, for
!  a description of these coefficients.
!
        a1 = -2.0*smooth*(xt*yt+xn*yn)*(xn*xn+yn*yn)/xjac**3
        a2 = 4.0*smooth*(xt*yt+xn*yn)*(xt*xn+yt*yn)/xjac**3
        a3 = -2.0*smooth*(xt*yt+xn*yn)*(xt*xt+yt*yt)/xjac**3

        b1 = 2.0*smooth*(yt*yt+yn*yn)*(xn*xn+yn*yn)/xjac**3
        b2 = -4.0*smooth*(yt*yt+yn*yn)*(xt*xn+yt*yn)/xjac**3
        b3 = 2.0*smooth*(yt*yt+yn*yn)*(xt*xt+yt*yt)/xjac**3

        c1 = 2.0*smooth*(xt*xt+xn*xn)*(xn*xn+yn*yn)/xjac**3
        c2 = -4.0*smooth*(xt*xt+xn*xn)*(xt*xn+yt*yn)/xjac**3
        c3 = 2.0*smooth*(xt*xt+xn*xn)*(xt*xt+yt*yt)/xjac**3

        a1 = a1+orth*pd(i,j)*xn*yn
        a2 = a2+orth*pd(i,j)*(xt*yn+xn*yt)
        a3 = a3+orth*pd(i,j)*xt*yt

        b1 = b1+orth*pd(i,j)*xn*xn
        b2 = b2+2.0*orth*pd(i,j)*(2.0*xt*xn+yt*yn)
        b3 = b3+orth*pd(i,j)*xt*xt

        c1 = c1+orth*pd(i,j)*yn*yn
        c2 = c2+2.0*orth*pd(i,j)*(xt*xn+2.0*yt*yn)
        c3 = c3+orth*pd(i,j)*yt*yt

        a1 = a1-2.0*cvn*pb(i,j)*xn*yn
        a2 = a2+2.0*cvn*pb(i,j)*(xt*yn+xn*yt)
        a3 = a3-2.0*cvn*pb(i,j)*xt*yt

        b1 = b1+2.0*cvn*pb(i,j)*yn*yn
        b2 = b2-4.0*cvn*pb(i,j)*yt*yn
        b3 = b3+2.0*cvn*pb(i,j)*yt*yt

        c1 = c1+2.0*cvn*pb(i,j)*xn*xn
        c2 = c2-4.0*cvn*pb(i,j)*xt*xn
        c3 = c3+2.0*cvn*pb(i,j)*xt*xt
!
!  Now compute terms from the derivatives of the weight functions.
!
        gt = 0.5*(pd(i+1,j)-pd(i-1,j))
        gn = 0.5*(pd(i,j+1)-pd(i,j-1))
        gx = orth*0.5*(gt*yn-gn*yt)*(xt*xn+yt*yn)**2/xjac
        gy = orth*0.5*(gn*xt-gt*xn)*(xt*xn+yt*yn)**2/xjac

        wt = 0.5*(pb(i+1,j)-pb(i-1,j))
        wn = 0.5*(pb(i,j+1)-pb(i,j-1))
        wx = cvn*xjac*(yn*wt-yt*wn)
        wy = cvn*xjac*(xt*wn-xn*wt)
!
!  (RES1, RES2) is the residual that should be driven to zero.
!
        xtt = xc(i+1,j)-2.0*xc(i,j)+xc(i-1,j)
        xnn = xc(i,j+1)-2.0*xc(i,j)+xc(i,j-1)
        ytt = yc(i+1,j)-2.0*yc(i,j)+yc(i-1,j)
        ynn = yc(i,j+1)-2.0*yc(i,j)+yc(i,j-1)
        xtn = 0.25*(xc(i+1,j+1)+xc(i-1,j-1)-xc(i-1,j+1)-xc(i+1,j-1))
        ytn = 0.25*(yc(i+1,j+1)+yc(i-1,j-1)-yc(i-1,j+1)-yc(i+1,j-1))

        res1 = b1*xtt+b2*xtn+b3*xnn+a1*ytt+a2*ytn+a3*ynn+gx+wx
        res2 = a1*xtt+a2*xtn+a3*xnn+c1*ytt+c2*ytn+c3*ynn+gy+wy

        if ( abs(res1)+abs(res2) > rmax ) then
          rmax = abs(res1)+abs(res2)
        end if

        a11 = -2.0D+00 * (b1+b3)
        a12 = -2.0D+00 * (a1+a3)
        a22 = -2.0D+00 * (c1+c3)
        a21 = a12

        det = a11*a22-a12*a21
!
!  (DX, DY) is the pointwise solution of the linear system.
!
        dx = (res2*a21-res1*a22)/det
        dy = (res1*a12-res2*a11)/det
!
!  Now look at how movement of XC(I,J), YC(I,J) affects the neighbors
!
!    (I+1,J-1)     (I+1,J)   (I+1,J+1)
!
!    (I,  J-1)     (I,  J)   (I,  J+1)
!
!    (I-1,J-1)     (I-1,J)   (I-1,J+1)
!
!  to determine the linear factor SCALE that will multiply (DX,DY).
!
        ratio = 0.25D+00
        xij = xc(i,j)
        yij = yc(i,j)

        xdisp = xc(i+1,j)-xij
        ydisp = yc(i+1,j)-yij
        vmag2 = xdisp*xdisp+ydisp*ydisp
        dot = dx*xdisp+dy*ydisp
        ratio = max(dot/vmag2,ratio)

        xdisp = xc(i+1,j+1)-xij
        ydisp = yc(i+1,j+1)-yij
        vmag2 = xdisp*xdisp+ydisp*ydisp
        dot = dx*xdisp+dy*ydisp
        ratio = max(dot/vmag2,ratio)

        xdisp = xc(i,j+1)-xij
        ydisp = yc(i,j+1)-yij
        vmag2 = xdisp*xdisp+ydisp*ydisp
        dot = dx*xdisp+dy*ydisp
        ratio = max(dot/vmag2,ratio)

        xdisp = xc(i-1,j+1)-xij
        ydisp = yc(i-1,j+1)-yij
        vmag2 = xdisp*xdisp+ydisp*ydisp
        dot = dx*xdisp+dy*ydisp
        ratio = max(dot/vmag2,ratio)

        xdisp = xc(i-1,j)-xij
        ydisp = yc(i-1,j)-yij
        vmag2 = xdisp*xdisp+ydisp*ydisp
        dot = dx*xdisp+dy*ydisp
        ratio = max(dot/vmag2,ratio)

        xdisp = xc(i-1,j-1)-xij
        ydisp = yc(i-1,j-1)-yij
        vmag2 = xdisp*xdisp+ydisp*ydisp
        dot = dx*xdisp+dy*ydisp
        ratio = max(dot/vmag2,ratio)

        xdisp = xc(i,j-1)-xij
        ydisp = yc(i,j-1)-yij
        vmag2 = xdisp*xdisp+ydisp*ydisp
        dot = dx*xdisp+dy*ydisp
        ratio = max(dot/vmag2,ratio)

        xdisp = xc(i+1,j-1)-xij
        ydisp = yc(i+1,j-1)-yij
        vmag2 = xdisp*xdisp+ydisp*ydisp
        dot = dx*xdisp+dy*ydisp
        ratio = max(dot/vmag2,ratio)

        scale = min ( 0.5D+00, 0.25D+00 / ratio )

        dxc(i,j) = scale*dx
        dyc(i,j) = scale*dy

        ssmot = ssmot+(xt*xt+xn*xn+yt*yt+yn*yn)/xjac
        ssweg = ssweg+cvn*pb(i,j)*xjac**2

      end do
    end do
!
!  Move the points.
!
    do j = 3, m1-1
      do i = 3, l1-1

        if ( i /= icrys ) then
          xc(i,j) = xc(i,j)+dxc(i,j)
        end if

        if ( j /= jcrys ) then
          yc(i,j) = yc(i,j)+dyc(i,j)
        end if

      end do
    end do

    do i = icrys+1, l1-1
      yc(i,jcrys) = 0.4D+00 &
        + 0.02D+00 * sin(xc(i,jcrys) * 20.0D+00 * 3.14159D+00 )
    end do
!
!  Handle the nodes along the bottom row, from 10 positions to the right
!  of the crystal, to the right hand wall.
!
    do j = jcrys+10, m1-1

      ratio = 0.25D+00
      call findp(j,xc(3,j),yc(3,j),xpp,ypp,xc,yc)
      dx = xpp-xc(2,j)
      dy = ypp-yc(2,j)

      xdisp = xc(2,j+1)-xc(2,j)
      ydisp = yc(2,j+1)-yc(2,j)
      vmag2 = xdisp*xdisp+ydisp*ydisp
      dot = dx*xdisp+dy*ydisp
      ratio = max(dot/vmag2,ratio)
      s1 = vmag2

      xdisp = xc(2,j-1)-xc(2,j)
      ydisp = yc(2,j-1)-yc(2,j)
      vmag2 = xdisp*xdisp+ydisp*ydisp
      dot = dx*xdisp+dy*ydisp
      ratio = max(dot/vmag2,ratio)
      s2 = max(vmag2,s1)

      scale = min ( 0.3D+00, 0.25D+00 / ratio )
      s1 = (scale*dx)**2+(scale*dy)**2

      if ( s1 < s2 ) then
        xc(2,j) = xc(2,j)+scale*dx
        yc(2,j) = yc(2,j)+scale*dy
      end if

      call cubic(nbot,yc(2,j),ybot,xc(2,j),xbot)

      yc(l1,j) = yc(l1-1,j)

    end do
!
!  Handle the nodes along the bottom row, from the left axis
!  of symmetry, to 9 positions beyond the crystal.
!
    do j = 2, jcrys+9

      yc(2,j) = yc(3,j)

      call cubic(nbot,yc(2,j),ybot,xc(2,j),xbot)

      yc(l1,j) = yc(l1-1,j)

    end do

    do i = 3, l1
      xc(i,2) = xc(i,3)
      xc(i,m1) = xc(i,m1-1)
      if ( xc(i,m1-1) <= (xc(2,m1)+0.002D+00 * float(i-2) ) ) then
        xc(i,m1) = xc(2,m1)+0.002D+00 * float(i-2)
      end if
    end do
!
!  Compute the maximum movement of all corner nodes.
!
    dxmax = abs(dxc(3,3))
    do j = 3, m1-1
      do i = 3, l1-1
        dxmax = max(dxmax,abs(dxc(i,j)))
        dxmax = max(dxmax,abs(dyc(i,j)))
      end do
    end do
!
!  Copy values to the dummy corner nodes with I = 1 or J=1.
!
    xc(1,1) = xc(2,2)
    yc(1,1) = yc(2,2)

    do j = 2, m1
      xc(1,j) = xc(2,j)
      yc(1,j) = yc(2,j)
    end do

    do i = 2, l1
      xc(i,1) = xc(i,2)
      yc(i,1) = yc(i,2)
    end do
!
!  Check for convergence.
!
    if ( dxmax <= epsad ) then
      return
    end if

  end do

  return
end
subroutine cubic ( nbot, ynew, ybot, xnew, xbot )

!*****************************************************************************80
!
!! CUBIC constructs a cubic spline through data.
!
!  Discussion:
!
!    CUBIC is given NBOT points (YBOT,XBOT) which lie along the curve
!    defining the bottom of the crucible.
!
!    CUBIC is also given a value YNEW which lies between YBOT(1) and
!    YBOT(NBOT).
!
!    CUBIC constructs a cubic spline through the data, and evaluates it at
!    YNEW, returning the value as XNEW.
!
!  Modified:
!
!    22 March 2002
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer nbot
  integer, parameter :: ni = 64

  real ( kind = 8 ) dome
  integer i
  integer ist
  real ( kind = 8 ) p2
  real ( kind = 8 ) p3
  real ( kind = 8 ) shift
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) tm(ni)
  real ( kind = 8 ) tt(ni)
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) xbot(nbot)
  real ( kind = 8 ) xl(ni)
  real ( kind = 8 ) xnew
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) ybot(nbot)
  real ( kind = 8 ) yl(ni)
  real ( kind = 8 ) ymin
  real ( kind = 8 ) ynew

  if ( ynew < ybot(1) .or. ybot(nbot) < ynew ) then
    write ( *, * ) ' '
    write ( *, * ) 'CUBIC - Fatal error!'
    write ( *, * ) '  YNEW  =  ',ynew
    write ( *, * ) '  outside of range YBOT(1) = ',ybot(1)
    write ( *, * ) '  through YBOT(NBOT) = ',ybot(nbot)
    stop
  end if

  xl(1:ni) = 0.0D+00
  yl(1:ni) = 0.0D+00
  tt(1:ni) = 0.0D+00
  tm(1:ni) = 0.0D+00

  yl(1) = ybot(1)+(ybot(1)-ybot(3))
  yl(2) = ybot(2)+(ybot(1)-ybot(3))

  do i = 3, nbot+2
    yl(i) = ybot(i-2)
    xl(i) = xbot(i-2)
  end do

  yl(nbot+3) = yl(nbot+1)+(yl(nbot+2)-yl(nbot))
  yl(nbot+4) = yl(nbot+2)+(yl(nbot+2)-yl(nbot))

  shift = (xl(3)-xl(4))/(yl(3)-yl(4))-(xl(4)-xl(5))/(yl(4)-yl(5))

  xl(2) = xl(3)+(yl(2)-yl(3))*(shift+(xl(3)-xl(4))/(yl(3)-yl(4)))
  xl(1) = xl(2)+(yl(1)-yl(2))*(shift+(xl(2)-xl(3))/(yl(2)-yl(3)))

  shift = (xl(nbot+2)-xl(nbot+1))/(yl(nbot+2)-yl(nbot+1)) &
    -(xl(nbot+1)-xl(nbot))/(yl(nbot+1)-yl(nbot))

  xl(nbot+3) = xl(nbot+2)+(yl(nbot+3)-yl(nbot+2))* &
    (shift+(xl(nbot+2)-xl(nbot+1))/(yl(nbot+2)-yl(nbot+1)))

  xl(nbot+4) = xl(nbot+3)+(yl(nbot+4)-yl(nbot+3))* &
    (shift+(xl(nbot+3)-xl(nbot+2))/(yl(nbot+3)-yl(nbot+2)))

  do i = 1, nbot+3
    tm(i) = (xl(i+1)-xl(i)) / (yl(i+1)-yl(i))
  end do

  do i = 3, nbot+2
    dome = abs ( tm(i+1) - tm(i) ) + abs ( tm(i-1) + tm(i-2) )
    if ( dome < 1.0E-10) then
      tt(i) = 0.0D+00
    else
      tt(i) = (abs(tm(i+1)-tm(i))*tm(i-1) &
        +abs(tm(i-1)-tm(i-2))*tm(i))/dome
    end if
  end do
!
!  Find the node YL(IST) which is nearest to YNEW.
!
  ymin = abs(yl(3)-ynew)
  ist = 3

  do i = 3, nbot+2
    if ( abs(yl(i)-ynew) < ymin ) then
      ist = i
      ymin = abs(yl(i)-ynew)
    end if
  end do

  if ( (yl(ist)-ynew) > 0.0D+00 ) then
    y1 = yl(ist-1)
    x1 = xl(ist-1)
    y2 = yl(ist)
    x2 = xl(ist)
    t1 = tt(ist-1)
    t2 = tt(ist)
  else
    y1 = yl(ist)
    x1 = xl(ist)
    y2 = yl(ist+1)
    x2 = xl(ist+1)
    t1 = tt(ist)
    t2 = tt(ist+1)
  end if
!
!  Evaluate the spline at YNEW.
!
  p2 = (3.0*(x2-x1)/(y2-y1)-2.0*t1-t2)/(y2-y1)
  p3 = (t1+t2-2.0*(x2-x1)/(y2-y1))/(y2-y1)**2

  xnew = x1+t1*(ynew-y1)+p2*(ynew-y1)**2+p3*(ynew-y1)**3

  return
end
subroutine diflow ( acof, diff, flow )

!*****************************************************************************80
!
!! DIFLOW computes the convection-diffusion coefficient.
!
!  Discussion:
!
!    DIFLOW computes ACOF, given FLOW, the mass velocity RHO*U,
!    and DIFF, the value of GAMMA/DELX.
!
!    Several schemes are available, but currently the power law is used.
!
!    See Patankar, chapter 5.
!
  implicit none

  real ( kind = 8 ) acof
  real ( kind = 8 ) diff
  real ( kind = 8 ) flow

  acof = diff

  if ( diff == 0.0D+00 ) then
    return
  end if
!
!  Power Law.
!
!  Define
!    FE = rho*u
!    DE = gamma/delx
!    Then the coefficient is
!
!    AE  =  DE * Max(0, (1-0.1*Abs(FE)/DE)**5) + Max(0,-FE)
!
  acof = diff * max(1.0d-10,(1.0D+00-0.1D+00*abs(flow/diff))**5)

  if ( flow < 0.0D+00 ) then
    acof = acof - flow
  end if

  return
end
subroutine doarea ( area, areal, areas, areat, icrys, jcrys, l0, m0, &
  ni, nj, xc, yc )

!*****************************************************************************80
!
!! DOAREA computes the area of the control volumes.
!
!  Discussion:
!
!    DOAREA is given (XC,YC), the locations of the "corners" of the
!    control volumes, and calculates the area of each control volume.
!
!    DOAREA uses the fact that the area of a polygon can be computed
!    by
!
!      AREA  =  0.5 * SUM (I=1 to N) X(I) * (Y(I+1)-Y(I-1))
!
!    where Y(0) should be replaced by Y(N), and Y(N+1) by Y(1).
!
!    Using this formula, AREA may come out negative, depending on
!    whether the nodes are given in clockwise or counterclockwise order.
!
  implicit none

  integer ni
  integer nj

  real ( kind = 8 ) area(ni,nj)
  real ( kind = 8 ) areal
  real ( kind = 8 ) areas
  real ( kind = 8 ) areat
  integer i
  integer icrys
  integer j
  integer jcrys
  integer l0
  integer m0
  integer nbad
  real ( kind = 8 ) xc(ni,nj)
  real ( kind = 8 ) yc(ni,nj)

  do i = 1, l0
    do j = 1, m0

      if ( i == 1 .or. i == l0 .or. j == 1 .or. j == m0 ) then
        area(i,j) = 0.0D+00
      else

        area(i,j) = 0.5D+00 * ( &
           xc(i,  j)*  (yc(i+1,j)  -yc(i,  j+1)) &
          +xc(i+1,j)*  (yc(i+1,j+1)-yc(i,  j)) &
          +xc(i+1,j+1)*(yc(i,  j+1)-yc(i+1,j)) &
          +xc(i,  j+1)*(yc(i,  j)  -yc(i+1,j+1)))

      end if

    end do
  end do
!
!  Check for illegal zero length sides.
!
  nbad = 0

  do i = 2, l0-1
    do j = 2, m0-1
      if ( xc(i,j) == xc(i+1,j) .and. yc(i,j) == yc(i+1,j) ) then
        nbad = nbad+1
        write ( *, * ) ' '
        write ( *, * ) 'SETARE - Fatal error!'
        write ( *, * ) '  Zero length side for cell I,J:',i,j
        write ( *, * ) 'XC,YC(I,J) = XC,YC(I+1,J)=',xc(i,j),yc(i,j)
      else if ( xc(i+1,j) == xc(i+1,j+1) .and. yc(i+1,j) == yc(i+1,j+1) ) then
        nbad = nbad+1
        write ( *, * ) ' '
        write ( *, * ) 'SETARE - Fatal error!'
        write ( *, * ) '  Zero length side for cell I,J:',i,j
        write ( *, * ) 'XC,YC(I+1,J) = XC,YC(I+1,J+1)=',xc(i+1,j),yc(i+1,j)
      else if ( xc(i+1,j+1) == xc(i,j+1) .and. yc(i+1,j+1) == yc(i,j+1) ) then
        nbad = nbad+1
        write ( *, * ) ' '
        write ( *, * ) 'SETARE - Fatal error!'
        write ( *, * ) '  Zero length side for cell I,J:',i,j
        write ( *, * ) 'XC,YC(I+1,J+1) = XC,YC(I,J+1)=',xc(i+1,j+1),yc(i+1,j+1)
      else if ( xc(i,j+1) == xc(i,j) .and. yc(i,j+1) == yc(i,j) ) then
        nbad = nbad+1
        write ( *, * ) ' '
        write ( *, * ) 'SETARE - Fatal error!'
        write ( *, * ) '  Zero length side for cell I,J:',i,j
        write ( *, * ) 'XC,YC(I,J+1) = XC,YC(I,J)=',xc(i,j+1),yc(i,j+1)
      end if
    end do
  end do

  if ( nbad > 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SETARE - Fatal error!'
    write ( *, * ) '  A total of ',nbad,' zero length cell sides.'
    stop
  end if
!
!  Check for illegal zero area cells.
!
  nbad = 0
  do i = 2, l0-1
    do j = 2, m0-1
      if ( area(i,j) == 0.0D+00 ) then
        nbad = nbad+1
        write ( *, * ) ' '
        write ( *, * ) 'SETARE - Fatal error!'
        write ( *, * ) '  Zero area for cell I = ',i,' J=',j
      end if
    end do
  end do

  if ( nbad > 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SETARE - Fatal error!'
    write ( *, * ) '  A total of ',nbad,' null cells.'
    stop
  end if
!
!  Compute liquid, solid, and total areas.
!
  areal = 0.0D+00
  do i = 1, icrys-1
    do j = 1, m0
      areal = areal + area(i,j)
    end do
  end do

  areas = 0.0D+00
  do i = icrys, l0
    do j = 1, jcrys
      areas = areas + area(i,j)
    end do
  end do

  areat = sum ( area(1:l0,1:m0) )

  return
end
subroutine findp ( j, xold, yold, xnew, ynew, xc, yc )

!*****************************************************************************80
!
!! FINDP finds the boundary nodes for a Neumann BC.
!
!  Discussion:
!
!    A Neumann boundary condition is applied to the derivative of a quantity.
!
!    FINDP is given a point (XOLD,YOLD).
!
!    FINDP returns a pair of values (XNEW,YNEW), which represent ???
!
  implicit none

  integer, parameter :: ni = 64
  integer, parameter :: nj = 64

  real ( kind = 8 ) a1
  real ( kind = 8 ) ai
  real ( kind = 8 ) aj
  real ( kind = 8 ) alp
  real ( kind = 8 ) b1
  real ( kind = 8 ) bet
  real ( kind = 8 ) c1
  real ( kind = 8 ) dels
  real ( kind = 8 ) denm
  real ( kind = 8 ) dxb1
  real ( kind = 8 ) dxb2
  real ( kind = 8 ) dyb1
  real ( kind = 8 ) dyb2
  integer j
  real ( kind = 8 ) r0
  real ( kind = 8 ) r1
  real ( kind = 8 ) r11
  real ( kind = 8 ) r2
  real ( kind = 8 ) r22
  real ( kind = 8 ) r33
  real ( kind = 8 ) slpi
  real ( kind = 8 ) xc(ni,nj)
  real ( kind = 8 ) xnew
  real ( kind = 8 ) xnew1
  real ( kind = 8 ) xnew2
  real ( kind = 8 ) xold
  real ( kind = 8 ) yc(ni,nj)
  real ( kind = 8 ) yn1
  real ( kind = 8 ) yn2
  real ( kind = 8 ) ynew
  real ( kind = 8 ) yold
!
!  Consider the three points P, Q and R, which have a constant
!  KSI coordinate, and an increasing ETA coordinate:
!
!    P = ( XC(2,J-1), YC(2,J-1) )
!    Q = ( XC(2,J),   YC(2,J) )
!    R = ( XC(2,J+1), YC(2,J+1) )
!
!
!  ^  R  =  [2,J+1]
!  |
!  E
!  T  Q  =  [2,J]
!  A
!  |
!  |  P  =  [2,J-1]
!  |
!  +------KSI------------->
!
!  The slope of the line from P to Q is (YQ-YP)/(XQ-XP), and
!  the slope of the line from Q to R is (YR-YQ)/(XR-XQ).
!
!  If these slopes are close enough, we can assume the three points
!  lie on a straight line.  So look at the size of
!    DENM  =  (YQ-YP)*(XR-XQ)-(YR-YQ)*(XQ-XP).
!
  dxb1 = xc(2,j)-xc(2,j-1)
  dxb2 = xc(2,j+1)-xc(2,j)
  dyb1 = yc(2,j)-yc(2,j-1)
  dyb2 = yc(2,j+1)-yc(2,j)
  denm = dyb2*dxb1-dyb1*dxb2

  if ( abs(denm) < 0.0001D+00 ) then

    if ( abs(dxb1) >= 0.001D+00 ) then
      slpi = -dyb1 / dxb1
      aj = yc(2,j) + slpi*xc(2,j)
      ai = xold - slpi*yold
      ynew = (aj-ai*slpi) / (1.0+slpi*slpi)
      xnew = slpi*ynew+ai
    else
      ynew = yold
      xnew = xc(2,j)
    end if
!
!  Use 3 points to find a circle.
!  (X-ALP)**2+(Y-BET)**2 = R0
!
!
!  Find the cross point between circle and straight line
!  x = ai+slpi*y get a1*y**2-2*b1*y+c1=0.
!
  else

    r11 = xc(2,j-1)**2+yc(2,j-1)**2
    r22 = xc(2,j)**2+yc(2,j)**2
    r33 = xc(2,j+1)**2+yc(2,j+1)**2

    r1 = 0.5D+00 * ( r22 - r11 )
    r2 = 0.5D+00 * ( r33 - r22 )

    alp = (r1*dyb2-r2*dyb1)/denm
    bet = (r2*dxb1-r1*dxb2)/denm

    r0 = (alp-xc(2,j))**2+(bet-yc(2,j))**2

    slpi = (xold-alp)/(yold-bet)
    ai = alp-slpi*bet

    a1 = 1.0 + slpi*slpi
    b1 = slpi * (alp-ai)+bet
    c1 = (alp-ai)**2 + bet*bet-r0

    dels = sqrt ( b1*b1 - a1 * c1 )

    yn1 = (b1+dels) / a1
    yn2 = (b1-dels) / a1

    xnew1 = slpi*yn1+ai
    xnew2 = slpi*yn2+ai

    r1 = (xnew1-xc(2,j))**2+(yn1-yc(2,j))**2
    r2 = (xnew2-xc(2,j))**2+(yn2-yc(2,j))**2

    if ( r1 < r2 ) then
      ynew = yn1
    else
      ynew = yn2
    end if

    xnew = slpi*ynew+ai

  end if

  return
end
subroutine flux ( ae1, ae2, aim, aip, ajm, ajp, ak1, ak2, ap, b1jbl, &
  b2jbl, cappa, cd, cmu, con, epsil, f, fjeta, fjksi, fmax, fn, fo, &
  gam, heta, hksi, icrys, isol, iter, izone, jcrys, l1, lblk, lconv, &
  lortho, m1, mode, nf, nsolve, ntimes, r, relax, res, rueta, ruksi )

!*****************************************************************************80
!
!! FLUX computes the value of the flux of various quantities.
!
!  Discussion:
!
!    It looks like if ISOL is 0, you just go through the big loop once.
!
  implicit none

  integer, parameter :: ni = 64
  integer, parameter :: nj = 64
  integer, parameter :: nk = 14
  integer, parameter :: nmaxij = 64
  integer, parameter :: ns = 10

  real ( kind = 8 ) acof
  real ( kind = 8 ) ae1(ni,nj)
  real ( kind = 8 ) ae2(ni,nj)
  real ( kind = 8 ) aim(ni,nj)
  real ( kind = 8 ) aip(ni,nj)
  real ( kind = 8 ) ajm(ni,nj)
  real ( kind = 8 ) ajp(ni,nj)
  real ( kind = 8 ) ak1(ni,nj)
  real ( kind = 8 ) ak2(ni,nj)
  real ( kind = 8 ) ap(ni,nj)
  real ( kind = 8 ) ap0(ni,nj)
  real ( kind = 8 ) b1jbl(nj)
  real ( kind = 8 ) b2jbl(nj)
  real ( kind = 8 ) b3jbl(ni)
  real ( kind = 8 ) cappa
  real ( kind = 8 ) cd
  real ( kind = 8 ) cmu
  real ( kind = 8 ) con(ni,nj)
  real ( kind = 8 ) con0(ni,nj)
  real ( kind = 8 ) diff
  real ( kind = 8 ) epsil
  real ( kind = 8 ) error
  real ( kind = 8 ) f(ni,nj,ns)
  real ( kind = 8 ) fjbi1(nj,ns)
  real ( kind = 8 ) fjbj1(ni,ns)
  real ( kind = 8 ) fjbl1(nj,ns)
  real ( kind = 8 ) fjbm1(ni,ns)
  real ( kind = 8 ) fjeta(ni,nj)
  real ( kind = 8 ) fjksi(ni,nj)
  real ( kind = 8 ) flow
  real ( kind = 8 ) fmax(ns)
  real ( kind = 8 ) fn(ni,nj,ns)
  real ( kind = 8 ) fo(ni,nj,ns)
  real ( kind = 8 ) gam(ni,nj)
  real ( kind = 8 ) heta(ni,nj)
  real ( kind = 8 ) hksi(ni,nj)
  integer i
  integer ii
  integer icrys
  integer isol
  integer iter
  integer izone
  integer j
  integer jcrys
  integer jj
  integer l1
  logical lblk(ns)
  logical lconv
  logical lortho
  integer m1
  integer mode
  integer nf
  integer nsolve(ns)
  integer nt
  integer ntimes(ns)
  real ( kind = 8 ) pt(nmaxij)
  real ( kind = 8 ) qt(nmaxij)
  real ( kind = 8 ) r(ni,nj)
  real ( kind = 8 ) relax(nk)
  real ( kind = 8 ) res(ns)
  real ( kind = 8 ) rueta(ni,nj)
  real ( kind = 8 ) ruksi(ni,nj)
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) t3
  real ( kind = 8 ) t4
  real ( kind = 8 ) tem1
  real ( kind = 8 ) tem2
  real ( kind = 8 ) xje
  real ( kind = 8 ) xjn
  real ( kind = 8 ) xjw
  real ( kind = 8 ) xjs
  real ( kind = 8 ) xme
  real ( kind = 8 ) xmn
  real ( kind = 8 ) xmw
  real ( kind = 8 ) xms
!
!  Save copies of CON and AP.
!
  con0(1:l1,1:m1) = con(1:l1,1:m1)
  ap0(1:l1,1:m1) = ap(1:l1,1:m1)
!
!  I think ITER = 0 occurs only for Temperature in the Solid zone...
!
  if ( iter == 0 ) then
    if ( isol /= 0 ) then

      call solve1(aim,aip,ajm,ajp,ap,cappa,cd,cmu,con,f,fo, &
        heta,hksi,icrys,izone,jcrys,l1,m1,nf)

      call solve2(aim,aip,ajm,ajp,ap,con,f,l1,lblk,m1,nf,nsolve)

    end if
    return
  end if

  do nt = 1, ntimes(nf)
!
!  Find the flux on the first step only.
!
    if ( nt == 1 ) then

      do j = 1, m1

        diff = gam(2,j)/(0.5*hksi(2,j))
        flow = -ruksi(2,j)
        call diflow(acof,diff,flow)
        qt(2) = ruksi(2,j)*f(2,j,nf)+acof*(f(1,j,nf)-f(2,j,nf))

        do i = 3, l1-1
          jj = j
          if ( j == 1) jj = 2
          if ( j == m1) jj = m1-1
          diff = gam(i,jj)*gam(i-1,jj)/(0.5*hksi(i,jj)*gam(i-1,jj) &
            +0.5*hksi(i-1,jj)*gam(i,jj))
          flow = -ruksi(i,j)
          call diflow(acof,diff,flow)
          qt(i) = ruksi(i,j)*f(i,j,nf)+acof*(f(i-1,j,nf)-f(i,j,nf))
        end do

        diff = gam(l1-1,j)/(0.5D+00*hksi(l1-1,j))
        flow = ruksi(l1,j)
        call diflow(acof,diff,flow)
        qt(l1) = ruksi(l1,j)*f(l1-1,j,nf)+acof*(f(l1-1,j,nf)-f(l1,j,nf))

        do i = 2, l1
          fjksi(i,j) = qt(i)
        end do

      end do
!
!  Along J.
!
      do i = 1, l1

        diff = gam(i,2)/(0.5D+00*heta(i,2))
        flow = -rueta(i,2)
        call diflow(acof,diff,flow)
        qt(2) = rueta(i,2)*f(i,2,nf)+acof*(f(i,1,nf)-f(i,2,nf))

        do j = 3, m1-1

          if ( i == 1 ) then
            ii = 2
          else if ( i == l1 ) then
            ii = l1-1
          else
            ii = i
          end if

          diff = gam(ii,j)*gam(ii,j-1)/(0.5*heta(ii,j)*gam(ii,j-1) &
            +0.5*heta(ii,j-1)*gam(ii,j))
          flow = -rueta(i,j)
          call diflow(acof,diff,flow)
          qt(j) = rueta(i,j)*f(i,j,nf)+acof*(f(i,j-1,nf)-f(i,j,nf))

        end do

        diff = gam(i,m1-1)/(0.5*heta(i,m1-1))
        flow = rueta(i,m1)
        call diflow(acof,diff,flow)
        qt(m1) = rueta(i,m1)*f(i,m1-1,nf)+acof*(f(i,m1-1,nf)-f(i,m1,nf))

        fjeta(i,2:m1) = qt(2:m1)

      end do
!
!  Finish jksi and jeta calculation.
!  Construct boundary flux for j type boundary condition.
!  Boundary condition con(1,j),*,*,* go into jksi and jeta
!  ak2/ak1,* is project on xi - direction from eta - direction
!  variable, in program let ak2 = 0.0, because sym. orth.
!
      do j = 2, m1-1

        if ( gam(1,j) == 0.0D+00 ) then
          fjksi(2,j) = (con(1,j)*heta(1,j)+ 0.5D+00 * ak2(2,j) &
            *(fjeta(1,j)+fjeta(1,j+1)))/ak1(2,j)
        end if

        if ( gam(l1,j) == 0.0D+00 ) then
          fjksi(l1,j) = (con(l1,j)*heta(l1,j)+ 0.5D+00 * ak2(l1,j) &
            *(fjeta(l1,j)+fjeta(l1,j+1)))/ak1(l1,j)
        end if

      end do

      do i = 2, l1-1

        if ( gam(i,1) == 0.0D+00 ) then
          fjeta(i,2) = (con(i,1)*hksi(i,1)+ 0.5D+00 * ae2(i,2) &
            *(fjksi(i,1)+fjksi(i+1,1)))/ae1(i,2)
        end if

        if ( gam(i,m1) == 0.0D+00 ) then
          fjeta(i,m1) = (con(i,m1)*hksi(i,m1)+ 0.5D+00 * ae2(i,m1) &
            *(fjksi(i,m1)+fjksi(i+1,m1)))/ae1(i,m1)
        end if

      end do
!
!  (2-31) and (2-42) find out vector J cdot vector n or J_n in bound.
!
      do i = 2, l1-1

        t1 = 0.5*(fjksi(i,1)+fjksi(i+1,1))
        t2 = 0.5*(fjksi(i,m1)+fjksi(i+1,m1))
        fjbj1(i,nf) = (fjeta(i,2)*ae1(i,2)-t1*ae2(i,2))/hksi(i,1)
        fjbm1(i,nf) = (fjeta(i,m1)*ae1(i,m1)-t2*ae2(i,m1))/hksi(i,m1)

        if ( mode == 1 ) then
          fjbj1(i,nf) = fjbj1(i,nf)/r(i,2)
          fjbm1(i,nf) = fjbm1(i,nf)/r(i,m1)
        end if

      end do

      do j = 2, m1-1

        t1 = 0.5*(fjeta(1,j)+fjeta(1,j+1))
        t2 = 0.5*(fjeta(l1,j)+fjeta(l1,j+1))
        fjbi1(j,nf) = (fjksi(2,j)*ak1(2,j)-t1*ak2(2,j))/heta(1,j)
        fjbl1(j,nf) = (fjksi(l1,j)*ak1(l1,j) &
          -t2*ak2(l1,j))/heta(l1,j)

        if ( mode == 1 ) then
          fjbi1(j,nf) = fjbi1(j,nf)/r(2,j)
          fjbl1(j,nf) = fjbl1(j,nf)/r(l1,j)
        end if

      end do
!
!  This code is only for the temperature variable.
!
      if ( nf == 5 ) then

        do j = 2, jcrys-1

          t1 = gam(icrys-1,j)/hksi(icrys-1,j)*(f(icrys-1,j,5)-f(icrys,j,5))
          b1jbl(j) = t1*ak1(icrys,j)/heta(icrys,j)

          if ( mode == 1 ) then
            b1jbl(j) = b1jbl(j)/r(icrys,j)
          end if

        end do

        b1jbl(jcrys) = b1jbl(jcrys+1)
        b1jbl(m1) = b1jbl(m1-1)

        do j = 2, jcrys-1

          t1 = gam(icrys+1,j)/hksi(icrys+1,j)*(f(icrys,j,5)-f(icrys+1,j,5))
          b2jbl(j) = t1*ak1(icrys,j)/heta(icrys,j)

          if ( mode == 1 ) then
            b2jbl(j) = b2jbl(j)/r(icrys,j)
          end if

        end do

        b2jbl(jcrys) = b2jbl(jcrys+1)
        b2jbl(m1) = b2jbl(m1-1)

        do j = jcrys+1, m1-1

          t1 = gam(icrys-1,j)/hksi(icrys-1,j)*(f(icrys-1,j,5)-f(icrys,j,5))
          t2 = gam(icrys-1,j)/heta(icrys,j)*(f(icrys,j+1,5)-f(icrys,j,5))
          b3jbl(j) = (t1*ak1(icrys,j)-t2*ak2(icrys,j))/heta(icrys,j)

          if ( mode == 1 ) then
            b3jbl(j) = b3jbl(j)/r(icrys,j)
          end if

        end do

        b3jbl(jcrys) = b3jbl(jcrys+1)
        b3jbl(m1) = b3jbl(m1-1)

      end if

    else
!
!  This code is done if this is NOT the first iteration.
!
!  correct j from known phi and j
!  correct jksi
!
      do j = 2, m1-1

        if ( gam(1,j) /= 0.0D+00 ) then
          fjksi(2,j) = fjksi(2,j)+aim(2,j)/ak1(2,j)* &
            (f(1,j,nf)-f(2,j,nf))+ruksi(2,j)*f(2,j,nf)
        end if

        do i = 3, l1

          if ( gam(l1,j) /= 0.0D+00 .or. i.ne.l1 ) then
            fjksi(i,j) = fjksi(i,j)+aip(i-1,j)/ak1(i,j)*(f(i-1,j,nf) &
              -f(i,j,nf))+ruksi(i,j)*f(i-1,j,nf)
          end if

        end do

      end do
!
!  Correct FJETA.
!
      do i = 2, l1-1

        if ( gam(i,1) /= 0.0D+00 ) then
          fjeta(i,2) = fjeta(i,2)+ajm(i,2)/ae1(i,2)* &
            (f(i,1,nf)-f(i,2,nf))+rueta(i,2)*f(i,2,nf)
        end if

        do j = 3, m1

          if ( gam(i,m1) /= 0.0D+00 .or. j.ne.m1 ) then
            fjeta(i,j) = fjeta(i,j)+ajp(i,j-1)/ae1(i,j)* &
              (f(i,j-1,nf)-f(i,j,nf))+rueta(i,j)*f(i,j-1,nf)
          end if

        end do
      end do

    end if
!
!  End of "Is this the first iteration or not" block.
!
!  Restore old values of CON and AP.
!
    do i = 1, l1
      do j = 1, m1
        con(i,j) = con0(i,j)
        ap(i,j) = ap0(i,j)
      end do
    end do

    if ( .not. lortho ) then

      do j = 2, m1-1

        do i = 1, l1
          pt(i) = 0.5*(fjeta(i,j)+fjeta(i,j+1))
          qt(i) = 0.5*(rueta(i,j)+rueta(i,j+1))
        end do

        do i = 2, l1-1
          t1 = hksi(i,j)/(hksi(i+1,j)+hksi(i,j))
          t2 = hksi(i+1,j)/(hksi(i+1,j)+hksi(i,j))
          t3 = hksi(i,j)/(hksi(i-1,j)+hksi(i,j))
          t4 = hksi(i-1,j)/(hksi(i-1,j)+hksi(i,j))
          xje = t2*pt(i)+t1*pt(i+1)
          xme = t2*qt(i)+t1*qt(i+1)
          xjw = t4*pt(i)+t3*pt(i-1)
          xmw = t4*qt(i)+t3*qt(i-1)
          con(i,j) = con(i,j)+(xje-f(i,j,nf)*xme)*ak2(i+1,j) &
            -(xjw-f(i,j,nf)*xmw)*ak2(i,j)
        end do
      end do

      do i = 2, l1-1

        do j = 1, m1
          pt(j) = 0.5D+00 * (fjksi(i,j)+fjksi(i+1,j))
          qt(j) = 0.5D+00 * (ruksi(i,j)+ruksi(i+1,j))
        end do

        do j = 2, m1-1
          tem1 = heta(i,j+1)+heta(i,j)
          tem2 = heta(i,j-1)+heta(i,j)
          t1 = heta(i,j)/tem1
          t2 = heta(i,j+1)/tem1
          t3 = heta(i,j)/tem2
          t4 = heta(i,j-1)/tem2
          xjn = t2*pt(j)+t1*pt(j+1)
          xmn = t2*qt(j)+t1*qt(j+1)
          xjs = t4*pt(j)+t3*pt(j-1)
          xms = t4*qt(j)+t3*qt(j-1)
          con(i,j) = con(i,j)+(xjn-f(i,j,nf)*xmn)*ae2(i,j+1) &
            -(xjs-f(i,j,nf)*xms)*ae2(i,j)
        end do
      end do

    end if
!
!  calculate jhat
!  step 5 use equation (2-97) and (2-102), (2-103), (2-104)
!  calculate boundary data from J and other purpose for b_sp.
!  from here to end. pc  =  jksi hat or jeta hat
!
!  jksi hat
!
    do j = 2, m1-1

      f(2,j,4) = 0.0D+00

      if ( gam(1,j) == 0.0D+00 ) then

        t1 = fjksi(2,j)-ruksi(2,j)*f(2,j,nf)
        diff = gam(2,j)/(0.5D+00 * hksi(2,j))
        flow = -ruksi(2,j)
        call diflow(acof,diff,flow)

        if ( acof /= 0.0D+00 ) then
          f(1,j,nf) = f(2,j,nf)+t1/acof
        else
          f(1,j,nf) = f(2,j,nf)
        end if

        f(2,j,4) = fjksi(2,j)

      end if

      f(l1,j,4) = 0.0D+00

      if ( gam(l1,j) == 0.0D+00 ) then

        t1 = fjksi(l1,j)-ruksi(l1,j)*f(l1-1,j,nf)
        diff = gam(l1-1,j)/(0.5*hksi(l1-1,j))
        flow = ruksi(l1,j)
        call diflow(acof,diff,flow)

        if ( acof /= 0.0D+00 ) then
          f(l1,j,nf) = f(l1-1,j,nf)-t1/acof
        else
          f(l1,j,nf) = f(l1-1,j,nf)
        end if

        f(l1,j,4) = fjksi(l1,j)

      end if

      do i = 3, l1-1
        f(i,j,4) = 0.0D+00
      end do

    end do

    do j = 2, m1-1
      do i = 2, l1
        fjksi(i,j) = f(i,j,4)
      end do
    end do
!
!  jeta hat
!
    do i = 2, l1-1

      f(i,2,4) = 0.0D+00

      if ( gam(i,1) == 0.0D+00 ) then

        t1 = fjeta(i,2)-rueta(i,2)*f(i,2,nf)
        diff = gam(i,2)/( 0.5D+00 * heta(i,2) )
        flow = -rueta(i,2)
        call diflow(acof,diff,flow)

        if ( acof /= 0.0D+00 ) then
          f(i,1,nf) = f(i,2,nf)+t1/acof
        else
          f(i,1,nf) = f(i,2,nf)
        end if

        f(i,2,4) = fjeta(i,2)

      end if

      f(i,m1,4) = 0.0D+00

      if ( gam(i,m1) == 0.0D+00 ) then

        t1 = fjeta(i,m1)-rueta(i,m1)*f(i,m1-1,nf)
        diff = gam(i,m1-1)/(0.5*heta(i,m1-1))
        flow = rueta(i,m1)
        call diflow(acof,diff,flow)

        if ( acof /= 0.0D+00 ) then
          f(i,m1,nf) = f(i,m1-1,nf)-t1/acof
        else
          f(i,m1,nf) = f(i,m1-1,nf)
        end if

        f(i,m1,4) = fjeta(i,m1)

      end if

      do j = 3, m1-1
        f(i,j,4) = 0.0D+00
      end do

    end do

    do i = 2, l1-1
      do j = 2, m1
        fjeta(i,j) = f(i,j,4)
      end do
    end do
!
!  Calculate and solve phy equation
!  (2-106) and (2-107) b  =  b_s + b_no + b_sp, t1=b_sp
!
    do i = 2, l1-1
      do j = 2, m1-1
        con(i,j) = con(i,j)+fjksi(i,j)*ak1(i,j)-fjksi(i+1,j)*ak1(i+1,j) &
          +fjeta(i,j)*ae1(i,j)-fjeta(i,j+1)*ae1(i,j+1)
      end do
    end do

    if ( isol == 0) then
      return
    end if

    do i = 2, l1-1
      do j = 2, m1-1
        ap(i,j) = ap(i,j)/relax(nf)
        con(i,j) = con(i,j)+(1.0-relax(nf))*ap(i,j)*f(i,j,nf)
      end do
    end do

    call solve1(aim,aip,ajm,ajp,ap,cappa,cd,cmu,con,f,fo, &
      heta,hksi,icrys,izone,jcrys,l1,m1,nf)

    call solve2(aim,aip,ajm,ajp,ap,con,f,l1,lblk,m1,nf,nsolve)
!
!  BEGIN NEW.
!
!  Compute the largest solution component.
!
    fmax(nf) = 0.0D+00
    do i = 1, l1
      do j = 1, m1
        fmax(nf) = max(fmax(nf),abs(f(i,j,nf)))
      end do
    end do
!
!  Compare current and previous iterates.
!
    if ( nt > 1 ) then
      res(nf) = 0.0D+00
      do i = 1, l1
        do j = 1, m1

          error = abs(f(i,j,nf)-fn(i,j,nf))

          if ( fmax(nf) > 0.0D+00 ) then
            error = error/fmax(nf)
          end if

          res(nf) = max(res(nf),error)

        end do
      end do
    end if
!
!  Save the new solution in FN.
!
    do i = 1, l1
      do j = 1, m1
        fn(i,j,nf) = f(i,j,nf)
      end do
    end do

    if ( res(nf) < epsil .and. nt > 1 ) then
!
!     write ( *, * ) 'FLUX - Convergence for variable ',nf,' on step ',nt
!     write ( *, * ) '  RES(NF) = ',res(nf),' FMAX(NF)=',fmax(nf)
!
      return
    end if

  end do
!
!  End of big iteration
!
!  Compute the largest solution component.
!
!     fmax(nf) = 0.0D+00
!     do i = 1, l1
!   do j = 1, m1
!     fmax(nf) = max(fmax(nf),abs(f(i,j,nf)))
!   end do
!     end do
!
!  Compute the largest relative change in the solution.
!
!     res(nf) = 0.0D+00
!     do i = 1, l1
!   do j = 1, m1
!
!     error = abs(f(i,j,nf)-fn(i,j,nf))
!
!     if ( fmax(nf) > 0.0D+00 ) then
!       error = error/fmax(nf)
!     end if
!
!     res(nf) = max(res(nf),error)
!
!   end do
!     end do
!
!  Decide whether to declare convergence or not.
!
  if ( res(nf) > epsil ) then
    lconv = .false.
  end if
!
!  Save the new solution in FN.
!
!     do i = 1, l1
!   do j = 1, m1
!     fn(i,j,nf) = f(i,j,nf)
!   end do
!     end do

!     nt = ntimes(nf)
!     write ( *, * ) 'FLUX - No convergence for variable ',nf,' on step ',nt
!     write ( *, * ) '  RES(NF) = ',res(nf),' FMAX(NF)=',fmax(nf)

  return
end
subroutine gamsor ( ap, birad, ce1, ce2, cmu, con, ewall, f, fcsl, fksl, &
  fma, fo, frsl, gam, gamt, grash, hamag, heta, hksi, icrys, inturb, &
  izone, jcrys, l1, lsolve, m1, mode, nf, pr, r, rdtm, re, recb, rect, &
  rho, rhocon, rpr, sige, sigk, sigt, stel, stes, tal, tas, tf, tw, &
  x, xc, y, yc )

!*****************************************************************************80
!
!! GAMSOR sets the coefficients for the transport problems.
!
!  Discussion:
!
!    These coefficients are associated with U, V, T, W, TK, TE and E.
!
!    Note that GAM(I,1) = 0 means the gradient is zero at wall J=1.
!
  implicit none

  integer, parameter :: ni = 64
  integer, parameter :: nj = 64
  integer, parameter :: nmaxij = 64
  integer, parameter :: ns = 10

  real ( kind = 8 ) ap(ni,nj)
  real ( kind = 8 ) ar
  real ( kind = 8 ) birad
  real ( kind = 8 ) ce1
  real ( kind = 8 ) ce2
  real ( kind = 8 ) cm4
  real ( kind = 8 ) cmu
  real ( kind = 8 ) con(ni,nj)
  real ( kind = 8 ) dcdr
  real ( kind = 8 ) depdr
  real ( kind = 8 ) depdx
  real ( kind = 8 ) dudx(ni,nj)
  real ( kind = 8 ) dudy(ni,nj)
  real ( kind = 8 ) dvdx(ni,nj)
  real ( kind = 8 ) dvdy(ni,nj)
  real ( kind = 8 ) dwdx(ni,nj)
  real ( kind = 8 ) dwdy(ni,nj)
  real ( kind = 8 ) ewall
  real ( kind = 8 ) f(ni,nj,ns)
  real ( kind = 8 ) fcsl
  real ( kind = 8 ) fksl
  real ( kind = 8 ) fma
  real ( kind = 8 ) fo(ni,nj,ns)
  real ( kind = 8 ) frsl
  real ( kind = 8 ) gam(ni,nj)
  real ( kind = 8 ) gamt(ni,nj)
  real ( kind = 8 ) gdudx(ni,nj)
  real ( kind = 8 ) gdudy(ni,nj)
  real ( kind = 8 ) gdvdx(ni,nj)
  real ( kind = 8 ) gdvdy(ni,nj)
  real ( kind = 8 ) grash
  real ( kind = 8 ) hamag
  real ( kind = 8 ) heta(ni,nj)
  real ( kind = 8 ) hksi(ni,nj)
  integer i
  integer icrys
  integer inturb
  integer izone
  integer j
  integer jcrys
  integer l1
  logical lsolve(ns)
  integer m1
  integer mode
  integer modelm
  integer nf
  real ( kind = 8 ) p11
  real ( kind = 8 ) p12
  real ( kind = 8 ) p21
  real ( kind = 8 ) p22
  real ( kind = 8 ) p31
  real ( kind = 8 ) p32
  real ( kind = 8 ) p41
  real ( kind = 8 ) p42
  real ( kind = 8 ) pbig
  real ( kind = 8 ) pr
  real ( kind = 8 ) prod(ni,nj)
  real ( kind = 8 ) r(ni,nj)
  real ( kind = 8 ) rdtm
  real ( kind = 8 ) re
  real ( kind = 8 ) recb
  real ( kind = 8 ) rect
  real ( kind = 8 ) rho(ni,nj)
  real ( kind = 8 ) rhocon
  real ( kind = 8 ) rpr
  real ( kind = 8 ) sige
  real ( kind = 8 ) sigk
  real ( kind = 8 ) sigt
  real ( kind = 8 ) stel
  real ( kind = 8 ) stes
  real ( kind = 8 ) tal
  real ( kind = 8 ) tas
  real ( kind = 8 ) tf
  real ( kind = 8 ) tw
  real ( kind = 8 ) x(ni,nj)
  real ( kind = 8 ) xc(ni,nj)
  real ( kind = 8 ) xplus(nmaxij)
  real ( kind = 8 ) y(ni,nj)
  real ( kind = 8 ) yc(ni,nj)
  real ( kind = 8 ) yplus(nmaxij)

  ap(1:l1,1:m1) = - rdtm * rho(1:l1,1:m1)
  con(1:l1,1:m1) = 0.0D+00
  gam(1:l1,1:m1) = gamt(1:l1,1:m1) + rhocon / re
!
!  Horizontal velocity.
!
  if ( nf == 1 ) then

    do i = 2, l1-1
      do j = 2, m1-1
        con(i,j) = fo(i,j,1)*rho(i,j)*rdtm+rho(i,j)*f(i,j,5)*grash/re**2
      end do
    end do

    if ( inturb == 1 ) then

      gdudx(1:l1,1:m1) = 0.0D+00
      gdudy(1:l1,1:m1) = 0.0D+00
      gdvdx(1:l1,1:m1) = 0.0D+00
      gdvdy(1:l1,1:m1) = 0.0D+00

      do i = 2, l1-1
        do j = 2, m1-1
          call gradnt(heta,hksi,i,j,f(1,1,1),dudx(i,j),dudy(i,j),xc,yc)
          call gradnt(heta,hksi,i,j,f(1,1,2),dvdx(i,j),dvdy(i,j),xc,yc)
          gdudx(i,j) = gamt(i,j)*dudx(i,j)
          gdudy(i,j) = gamt(i,j)*dudy(i,j)
          gdvdx(i,j) = gamt(i,j)*dvdx(i,j)
          gdvdy(i,j) = gamt(i,j)*dvdy(i,j)
        end do
      end do

      do i = 2, l1-1
        do j = 2, m1-1
          call gradnt(heta,hksi,i,j,gdudx,p11,p12,xc,yc)
          call gradnt(heta,hksi,i,j,gdudy,p21,p22,xc,yc)
          call gradnt(heta,hksi,i,j,gdvdx,p31,p32,xc,yc)
          call gradnt(heta,hksi,i,j,gdvdy,p41,p42,xc,yc)
          con(i,j) = con(i,j)+p11+p32
        end do
      end do

      cm4 = cmu**0.25

      if ( mode == 0 ) then
        do i = 2, l1-1

          yplus(i) = re*rho(i,2)*sqrt(f(i,2,7))*cm4*0.5*heta(i,2)/rhocon

          if ( yplus(i) <= 11.63 ) then
            gam(i,1) = rhocon/re
          else
            gam(i,1) = yplus(i)/(2.5*log(ewall*yplus(i)))
          end if

        end do
      end if

      do i = 2, l1-1

        yplus(i) = re*rho(i,m1-1)*sqrt(f(i,m1-1,7))*cm4*0.5*heta(i,m1-1)/rhocon

        if ( yplus(i) <= 11.63 ) then
          gam(i,m1) = rhocon/re
        else
          gam(i,m1) = yplus(i)/(2.5*log(ewall*yplus(i)))
        end if

      end do

      do j = 2, m1-1

        xplus(j) = re*rho(2,j)*sqrt(f(2,j,7))*cm4*0.5*hksi(2,j)/rhocon

        if ( xplus(j) <= 11.63 ) then
          gam(1,j) = rhocon/re
        else
          gam(1,j) = xplus(j)/(2.5*log(ewall*xplus(j)))
        end if

      end do

      do j = 2, m1-1

        xplus(j) = re*rho(l1-1,j)*sqrt(f(l1-1,j,7))*cm4*0.5*hksi(l1-1,j)/rhocon

        if ( xplus(j) <= 11.63 ) then
          gam(l1,j) = rhocon/re
        else
          gam(l1,j) = xplus(j)/(2.5*log(ewall*xplus(j)))
        end if

      end do

    end if

    do i = 2, l1-1
      gam(i,1) = 0.0D+00
    end do

    do j = 2, m1-1
      f(1,j,1) = 0.0D+00
      f(l1,j,1) = 0.0D+00
    end do

    f(l1,1,1) = 0.0D+00
    do j = jcrys+1, m1-1
      f(l1,j,1) = fo(l1,j,2)*(xc(l1,j+1)-xc(l1,j))/(yc(l1,j+1)-yc(l1,j))
    end do
    f(l1,jcrys+1,1) = 0.0D+00
!
!  Vertical velocity.
!
  else if ( nf == 2 ) then

    do i = 2, l1-1
      do j = 2, m1-1
        con(i,j) = fo(i,j,2)*rho(i,j)*rdtm
        ap(i,j) = ap(i,j)-rho(i,j)*hamag**2/re
        if ( mode == 1 ) then
          con(i,j) = con(i,j)+rho(i,j)*f(i,j,6)**2/r(i,j)**3
          ap(i,j) = ap(i,j)-rho(i,j)/r(i,j)**2
        end if
      end do
    end do

    if ( inturb == 1 ) then

      gdudx(1:l1,1:m1) = 0.0D+00
      gdudy(1:l1,1:m1) = 0.0D+00
      gdvdx(1:l1,1:m1) = 0.0D+00
      gdvdy(1:l1,1:m1) = 0.0D+00

      do i = 2, l1-1
        do j = 2, m1-1
          call gradnt(heta,hksi,i,j,f(1,1,1),dudx(i,j),dudy(i,j),xc,yc)
          call gradnt(heta,hksi,i,j,f(1,1,2),dvdx(i,j),dvdy(i,j),xc,yc)
          gdudx(i,j) = gamt(i,j)*dudx(i,j)
          gdudy(i,j) = gamt(i,j)*dudy(i,j)
          gdvdx(i,j) = gamt(i,j)*dvdx(i,j)
          gdvdy(i,j) = gamt(i,j)*dvdy(i,j)
        end do
      end do

      do i = 2, l1-1
        do j = 2, m1-1
          call gradnt(heta,hksi,i,j,gdudx,p11,p12,xc,yc)
          call gradnt(heta,hksi,i,j,gdudy,p21,p22,xc,yc)
          call gradnt(heta,hksi,i,j,gdvdx,p31,p32,xc,yc)
          call gradnt(heta,hksi,i,j,gdvdy,p41,p42,xc,yc)
          con(i,j) = con(i,j)+p21+p42
        end do
      end do

      cm4 = cmu**0.25D+00

      if ( mode == 0 ) then
        do i = 2, l1-1

          yplus(i) = re*rho(i,2)*sqrt(f(i,2,7))*cm4*0.5*heta(i,2)/rhocon
          gam(i,1) = rhocon/re
          if ( yplus(i) > 11.63 ) then
            gam(i,1) = yplus(i) / (2.5*log(ewall*yplus(i)))
          end if

        end do
      end if

      do i = 2, l1-1

        yplus(i) = re * rho(i,m1-1) * sqrt ( f(i,m1-1,7) ) * cm4 &
          * 0.5 * heta(i,m1-1) / rhocon
        gam(i,m1) = rhocon / re
        if ( yplus(i) > 11.63 ) then
          gam(i,m1) = yplus(i) / (2.5*log(ewall*yplus(i)))
        end if

      end do

      do j = 2, m1-1
        xplus(j) = re*rho(2,j)*sqrt(f(2,j,7))*cm4*0.5*hksi(2,j)/rhocon
        gam(1,j) = rhocon / re
        if ( xplus(j) > 11.63 ) then
          gam(1,j) = xplus(j) / (2.5*log(ewall*xplus(j)))
        end if
      end do

      do j = 2, m1-1
        xplus(j) = re*rho(l1-1,j)*sqrt(f(l1-1,j,7))*cm4*0.5*hksi(l1-1,j)/rhocon
        gam(l1,j) = rhocon/re
        if ( xplus(j) > 11.63 ) then
          gam(l1,j) = xplus(j)/(2.5*log(ewall*xplus(j)))
        end if
      end do

    end if

    do i = 2, l1-1
      f(i,1,2) = 0.0D+00
    end do

    do j = jcrys+1, m1-1
      dcdr = (f(l1,j+1,5)-f(l1,j,5))/(y(l1,j+1)-y(l1,j))
      f(l1,j,2) = f(l1-1,j,2)+(xc(l1,j)-xc(l1-1,j))*dcdr*(fma/(re*pr))
    end do

    f(l1,jcrys+1,2) = 0.0D+00
    do j = 2, jcrys
      f(1,j,2) = 0.0D+00
      f(l1,j,2) = 0.0D+00
    end do

    f(l1,m1-1,2) = 0.0D+00
    f(l1,m1,2) = 0.0D+00
!
!  Temperature computations in the solid zone.
!
  else if ( nf == 5 .and. izone == 1 ) then

    do i = 1, l1
      do j = 1, m1
        gam(i,j) = rpr*fksl/fcsl/frsl
      end do
    end do

    do i = icrys+1, l1
      gam(i,1) = 0.0D+00
    end do

    do i = 2, l1-1
      do j = 2, m1-1
        con(i,j) = fo(i,j,5)*rho(i,j)*rdtm
      end do
    end do

    do j = 2, jcrys
      f(l1,j,5) = -stes / stel
    end do

    do i = icrys+1, l1
      f(i,m1,5) = f(i,m1-1,5)-birad*((f(i,m1,5)*(tw-tf)+tf)**4-tas**4) &
        *(y(i,m1)-y(i,m1-1))/100000000.0D+00
    end do
!
!  Temperature calculations in the liquid zone.
!
  else if ( nf == 5 .and. izone == 2 ) then

    do i = 1, l1
      do j = 1, m1
        gam(i,j) = rpr + gamt(i,j) / sigt
      end do
    end do

    do i = 2, l1-1
      do j = 2, m1-1
        con(i,j) = fo(i,j,5)*rho(i,j)*rdtm
      end do
    end do

    if ( inturb == 1 ) then

      cm4 = cmu**0.25D+00
      pbig = 9.24D+00*((pr/sigt)**0.75-1.0)*(1.0+0.28*exp( - 0.007D+00 * pr/sigt))

      if ( mode == 0) then
        do i = 2, l1-1
          yplus(i) = re*rho(i,2)*sqrt(f(i,2,7))*cm4 &
            * 0.5D+00 * heta(i,2)/rhocon
          gam(i,1) = rpr
          if ( yplus(i) > 11.63 ) then
            gam(i,1) = max(rhocon/re, &
              yplus(i)/(sigt*(2.5*log(ewall*yplus(i))+pbig)))
          end if
        end do
      end if

      do i = 2, l1-1
        yplus(i) = re*rho(i,m1-1)*sqrt(f(i,m1-1,7))*cm4*0.5*heta(i,m1-1)/rhocon
        gam(i,m1) = rpr
        if ( yplus(i) > 11.63D+00 ) then
          gam(i,m1) = max(rhocon/re,yplus(i)/ &
            (sigt*(2.5*log(ewall*yplus(i))+pbig)))
        end if
      end do

      do j = 2, m1-1
        xplus(j) = re*rho(2,j)*sqrt(f(2,j,7))*cm4*0.5*hksi(2,j)/rhocon
        gam(1,j) = rpr
        if ( xplus(j) > 11.63 ) then
          gam(1,j) = max(rhocon/re, &
            xplus(j)/(sigt*(2.5*log(ewall*xplus(j))+pbig)))
        end if
      end do

      do j = 2, m1-1
        xplus(j) = re*rho(l1-1,j)*sqrt(f(l1-1,j,7)) &
          *cm4*0.5*hksi(l1-1,j) / rhocon
        gam(l1,j) = rpr
        if ( xplus(j) > 11.63 ) then
          gam(l1,j) = max(rhocon/re, &
            xplus(j)/(sigt*(2.5*log(ewall*xplus(j))+pbig)))
        end if
      end do

    end if

    do j = 2, jcrys
      f(l1,j,5) = 0.0D+00
    end do

    do j = jcrys+1, m1
      f(l1,j,5) = f(l1-1,j,5)-0.1*birad*((f(l1,j,5)*(tw-tf)+tf)**4-tal**4) &
        *(x(l1,j)-x(l1-1,j)) / 100000000.0D+00
    end do

    f(1,1:m1,5) = 0.5D+00 + 0.5D+00 * yc(2,1:m1)
    f(1:l1,m1,5) = 1.0D+00
    gam(1:l1,1) = 0.0D+00
!
!  Angular momentum.
!
  else if ( nf == 6 ) then
!
!  modelm = 1. based on rectangular de/dx of Sabhapathy and Salcudean
!  modelm = 2. based on curvilinear de/dx of Sabhapathy and Salcudean.
!  modelm = 3. same as 2 but different numerical treatment for source
!    terms based on strong magnetic. 1,2 is good for weak Ha
!    3 is good for moderate.
!  modelm = 4. based on strong magnetic. rotation is eliminated.
!  modelm = 5. is a test.
!
    if ( hamag == 0.0D+00 ) then
      modelm = 0.0D+00
    else if ( hamag > 0.0D+00 .and. hamag <= 20.0D+00 ) then
      modelm = 2
    else if ( hamag > 20.0D+00 .and. hamag <= 40.0D+00 ) then
      modelm = 5
    else if ( hamag > 40.0D+00 ) then
      modelm = 4
    end if

    do i = 2, l1-1
      do j = 2, m1-1

        con(i,j) = fo(i,j,6)*rho(i,j)*rdtm

        if ( modelm == 1 ) then

          con(i,j) = con(i,j)-rho(i,j)*hamag**2/re/r(i,j)* &
            (f(i+1,j,9)-f(i,j,9))/(x(i+1,j)-x(i,j))

        else if ( modelm == 2 ) then

          call gradnt(heta,hksi,i,j,f(1,1,9),depdx,depdr,xc,yc)

          if ( j >= 3 ) then
            con(i,j) = con(i,j)-rho(i,j)*hamag**2/re/r(i,j)*depdx
          end if

        else if ( modelm == 3 ) then

          call gradnt(heta,hksi,i,j,f(1,1,9),depdx,depdr,xc,yc)

          if ( j >= 3 ) then
            con(i,j) = con(i,j)-rho(i,j)*hamag**2/re/r(i,j)*depdx &
              +rho(i,j)*hamag**2/re*f(i,j,6)
            ap(i,j) = ap(i,j)-rho(i,j)*hamag**2/re
          end if

        else if ( modelm == 4 ) then

          ap(i,j) = ap(i,j)-rho(i,j)*hamag**2/re

        else if ( modelm == 5 ) then

          call gradnt(heta,hksi,i,j,f(1,1,9),depdx,depdr,xc,yc)

          if ( j >= 20 ) then
            con(i,j) = con(i,j)-rho(i,j)*hamag**2/re/r(i,j)*depdx
          end if

          if ( j <= 20 ) then
            ap(i,j) = ap(i,j)-rho(i,j)*hamag**2/re
          end if

        end if
!
!  -2/r*dOmega/dr /re
!
        ar = 2.0/r(i,j)/(y(i,j)-y(i,j-1))/re
        con(i,j) = con(i,j)+rho(i,j)*ar*f(i,j-1,6)
        ap(i,j) = ap(i,j)-rho(i,j)*ar

      end do
    end do

    do i = 1, l1
      f(i,m1,6) = recb
      f(i,1,6) = 0.0D+00
    end do

    do j = jcrys+1, m1
      gam(l1,j) = 0.0D+00
    end do

    do j = 1, m1

      f(1,j,6) = recb*r(2,j)**2

      if ( j > jcrys ) then
        f(l1,j,6) = f(l1-1,j,6)
      else
        f(l1,j,6) = rect*r(l1,j)**2
      end if

    end do
!
!  Turbulent kinetic energy.
!
  else if ( nf == 7 ) then

    do i = 2, l1-1
      do j = 2, m1-1
        con(i,j) = fo(i,j,7)*rdtm
        gam(i,j) = rhocon/re+gamt(i,j)/sigk
      end do
    end do

    do i = 2, l1-1
      do j = 2, m1-1

        call gradnt(heta,hksi,i,j,f(1,1,1),dudx(i,j),dudy(i,j),xc,yc)
        call gradnt(heta,hksi,i,j,f(1,1,2),dvdx(i,j),dvdy(i,j),xc,yc)
        prod(i,j) = gamt(i,j)*(2.0*(dudx(i,j)**2+dvdy(i,j)**2) &
          +(dudy(i,j)+dvdx(i,j))**2)

        if ( mode == 1 ) then
          prod(i,j) = prod(i,j)+gamt(i,j)*2.0*(f(i,j,2)/r(i,j))**2
        end if

        if ( lsolve(6) ) then
          call gradnt(heta,hksi,i,j,f(1,1,6),dwdx(i,j),dwdy(i,j),xc,yc)
          prod(i,j) = prod(i,j)+gamt(i,j)*(dwdx(i,j)**2 &
            +(dwdy(i,j)-2.0*f(i,j,6)/r(i,j))**2)/r(i,j)**2
        end if

        con(i,j) = con(i,j)+prod(i,j)
        ap(i,j) = ap(i,j)-cmu*rho(i,j)**2*abs(f(i,j,7))/(gamt(i,j)+1.e-5)

      end do
    end do
!
!  Turbulent dissipation.
!
  else if ( nf == 8 ) then

    do i = 2, l1-1
      do j = 2, m1-1
        con(i,j) = fo(i,j,8)*rdtm
        gam(i,j) = rhocon/re+gamt(i,j)/sige
      end do
    end do

    do i = 2, l1-1
      do j = 2, m1-1
        con(i,j) = con(i,j)+ce1*cmu*rho(i,j)*abs(f(i,j,7))*prod(i,j) &
          /((gamt(i,j)+1.e-5))
        ap(i,j) = ap(i,j)-ce2*cmu*rho(i,j)**2*abs(f(i,j,7)) &
          /(gamt(i,j)+1.e-5)
      end do
    end do

  end if

  return
end
subroutine gradnt ( heta, hksi, i, j, phi, dphidx, dphidy, xc, yc )

!*****************************************************************************80
!
!! GRADNT calculates gradients at the primary nodes.
!
!  Discussion:
!
!    GRADNT is given the value of a state quantity PHI at the primary nodes,
!    the KSI and ETA spacing between primary nodes, and the coordinates of
!    the corner nodes that define the control volumes.
!
!    GRADNT calculates the gradient (dPHI/dX, dPHI/dY) at the primary nodes.
!
  implicit none

  integer, parameter :: ni = 64
  integer, parameter :: nj = 64

  real ( kind = 8 ) detadx
  real ( kind = 8 ) detady
  real ( kind = 8 ) dpdeta
  real ( kind = 8 ) dpdksi
  real ( kind = 8 ) dphidx
  real ( kind = 8 ) dphidy
  real ( kind = 8 ) dxdeta
  real ( kind = 8 ) dxdksi
  real ( kind = 8 ) dksidx
  real ( kind = 8 ) dksidy
  real ( kind = 8 ) dydeta
  real ( kind = 8 ) dydksi
  real ( kind = 8 ) heta(ni,nj)
  real ( kind = 8 ) hksi(ni,nj)
  integer i
  integer j
  real ( kind = 8 ) phi(ni,nj)
  real ( kind = 8 ) phia
  real ( kind = 8 ) phib
  real ( kind = 8 ) phic
  real ( kind = 8 ) phid
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) volume
  real ( kind = 8 ) xc(ni,nj)
  real ( kind = 8 ) xd
  real ( kind = 8 ) xm
  real ( kind = 8 ) xp
  real ( kind = 8 ) xu
  real ( kind = 8 ) yc(i,j)
  real ( kind = 8 ) yd
  real ( kind = 8 ) ym
  real ( kind = 8 ) yp
  real ( kind = 8 ) yu
!
!  Interpolate the value of PHI at the cell volume interfaces
!  from its value at primary nodes.
!
!  ^               (I,J+1)
!  |
!  |               PHI(B)
!  E
!  T  (I-1,J)  PHI(C)  (I,J)  PHI(A)  (I+1,J)
!  A
!  |               PHI(D)
!  |
!  |               (I,J-1)
!  |
!  +-------------------KSI------------->
!
  t1 = hksi(i,j)/(hksi(i,j)+hksi(i+1,j))
  t2 = hksi(i+1,j)/(hksi(i,j)+hksi(i+1,j))
  phia = t1*phi(i+1,j)+t2*phi(i,j)

  t1 = heta(i,j)/(heta(i,j)+heta(i,j+1))
  t2 = heta(i,j+1)/(heta(i,j)+heta(i,j+1))
  phib = t1*phi(i,j+1)+t2*phi(i,j)

  t1 = hksi(i,j)/(hksi(i,j)+hksi(i-1,j))
  t2 = hksi(i-1,j)/(hksi(i,j)+hksi(i-1,j))
  phic = t1*phi(i-1,j)+t2*phi(i,j)

  t1 = heta(i,j)/(heta(i,j)+heta(i,j-1))
  t2 = heta(i,j-1)/(heta(i,j)+heta(i,j-1))
  phid = t1*phi(i,j-1)+t2*phi(i,j)
!
!  Now subtract opposing values, assuming a nominal spacing of
!  1 for KSI and ETA, to get dPHI/dKSI and dPHI/dETA.
!
  dpdksi = phia-phic
  dpdeta = phib-phid
!
!  Now, using the coordinates of the "corner nodes", locate the midpoints
!  of the control volume interfaces that surround primary node (I,J).
!
!    ^
!    |     [I,J+1]  Up        [I+1,J+1]
!    E
!    T     Minus    (I,J)      Plus
!    A
!    |     [I,J]    Down       [I+1,J]
!    |
!    +-----XSI---->
!
  xp = 0.5 * (xc(i+1,j+1)+xc(i+1,j))
  yp = 0.5 * (yc(i+1,j+1)+yc(i+1,j))

  xm = 0.5 * (xc(i,j+1)+xc(i,j))
  ym = 0.5 * (yc(i,j+1)+yc(i,j))

  xu = 0.5 * (xc(i,j+1)+xc(i+1,j+1))
  yu = 0.5 * (yc(i,j+1)+yc(i+1,j+1))

  xd = 0.5 * (xc(i,j)+xc(i+1,j))
  yd = 0.5 * (yc(i,j)+yc(i+1,j))
!
!  From the coordinates of the control volume wall centers,
!  estimate dX/dKSI, dX/dETA, dY/dKSI and dY/dETA at the primary node,
!  again assuming a spacing of 1 for KSI and ETA between control
!  volume wall centers.
!
  dxdksi = xp-xm
  dxdeta = xu-xd
  dydksi = yp-ym
  dydeta = yu-yd
!
!  Compute the determinant of the Jacobian (XSI,ETA)-->(X,Y).
!
!  J(XSI,ETA)  =  ( dX/dKSI  dX/dETA)
!           ( dY/dKSI  dY/dETA)
!
  volume = dxdksi*dydeta-dxdeta*dydksi
!
!  Now compute the elements of the inverse Jacobian.
!
  dksidx = dydeta/volume
  dksidy = -dxdeta/volume
  detadx = -dydksi/volume
  detady = dxdksi/volume
!
!  Now finally compute dPHI/dX and dPHI/dY using the chain rule.
!
  dphidx = dpdksi*dksidx+dpdeta*detadx
  dphidy = dpdeta*detady+dpdksi*dksidy

  return
end
subroutine inidat ( ae1, ae2, ak1, ak2, area, b, birad, bo, cappa, cd, ce1, &
  ce2, cfo, cinc, cmu, cost, cvn, delt, dtm, epsad, epsil, ewall, f, fcsl, &
  fks,fksl,fma,fn,fnu,fo,fr,frsl,gamt,grash,hamag,heta,hf,hksi, &
  icost,icrys,inturb,iplot,ipref,iprint,jcrys, &
  jpref,l0,l1,last,lastt,lblk,lortho,lsolve,m0,m1,maxbot,mode, &
  nbot,ndt,ni,nj,nk,np,npar,npc,ns,nsolve,ntimes,orth,par,pr,ra,rdtm, &
  re,recb,rect,relax,rho,rhocon,rpr,rueta,ruksi,sige,sigk, &
  sigma,sigt,smooth,stel,stes,tal,tanca,tanca2,tas,tend,tf,tinit, &
  title,tnow,tw,vave,vol,x,xbot,xc,xlen,y,ybot,yc,ylen)

!*****************************************************************************80
!
!! INIDAT sets the initial values of certain data.
!
!  Modified:
!
!    22 March 2002
!
!  Parameters:
!
!    Output, integer IPLOT, controls restart file creation.
!    0, do not create a restart file.
!    1, create a restart file.
!
  implicit none

  integer maxbot
  integer ni
  integer nj
  integer nk
  integer npar
  integer ns

  real ( kind = 8 ) ae1(ni,nj)
  real ( kind = 8 ) ae2(ni,nj)
  real ( kind = 8 ) ak1(ni,nj)
  real ( kind = 8 ) ak2(ni,nj)
  real ( kind = 8 ) area(ni,nj)
  real ( kind = 8 ) b
  real ( kind = 8 ) birad
  real ( kind = 8 ) bo
  real ( kind = 8 ) cappa
  real ( kind = 8 ) cd
  real ( kind = 8 ) ce1
  real ( kind = 8 ) ce2
  real ( kind = 8 ) cfo
  real ( kind = 8 ) cinc
  real ( kind = 8 ) cl
  real ( kind = 8 ) cmu
  real ( kind = 8 ) cost
  real ( kind = 8 ) cs
  real ( kind = 8 ) cvn
  real ( kind = 8 ) delt
  real ( kind = 8 ) dtm
  real ( kind = 8 ) epsad
  real ( kind = 8 ) epsil
  real ( kind = 8 ) ewall
  real ( kind = 8 ) f(ni,nj,ns)
  real ( kind = 8 ) fcsl
  real ( kind = 8 ) fkl
  real ( kind = 8 ) fks
  real ( kind = 8 ) fksl
  real ( kind = 8 ) fma
  real ( kind = 8 ) fn(ni,nj,ns)
  real ( kind = 8 ) fnu
  real ( kind = 8 ) fo(ni,nj,ns)
  real ( kind = 8 ) fr
  real ( kind = 8 ) frsl
  real ( kind = 8 ) gamt(ni,nj)
  real ( kind = 8 ) grash
  real ( kind = 8 ) grav
  real ( kind = 8 ) hamag
  real ( kind = 8 ) heta(ni,nj)
  real ( kind = 8 ) hf
  real ( kind = 8 ) hksi(ni,nj)
  integer i
  integer icost
  integer icrys
  integer inturb
  integer iplot
  integer ipref
  integer iprint
  integer j
  integer jcrys
  integer jpref
  integer k
  integer l0
  integer l1
  integer last
  integer lastt
  logical lblk(ns)
  logical lortho
  logical lsolve(ns)
  integer m0
  integer m1
  integer mode
  integer nbot
  integer ndt
  integer np
  integer npc
  integer nsolve(ns)
  integer ntimes(ns)
  real ( kind = 8 ) orth
  real ( kind = 8 ) par(npar)
  real ( kind = 8 ) pr
  real ( kind = 8 ) ra
  real ( kind = 8 ) rdtm
  real ( kind = 8 ) re
  real ( kind = 8 ) recb
  real ( kind = 8 ) rect
  real ( kind = 8 ) relax(nk)
  real ( kind = 8 ) rho(ni,nj)
  real ( kind = 8 ) rhocon
  real ( kind = 8 ) rhol
  real ( kind = 8 ) rhos
  real ( kind = 8 ) rpr
  real ( kind = 8 ) rueta(ni,nj)
  real ( kind = 8 ) ruksi(ni,nj)
  real ( kind = 8 ) sige
  real ( kind = 8 ) sigk
  real ( kind = 8 ) sigma
  real ( kind = 8 ) sigt
  real ( kind = 8 ) smooth
  real ( kind = 8 ) stel
  real ( kind = 8 ) stes
  real ( kind = 8 ) tal
  real ( kind = 8 ) tanca
  real ( kind = 8 ) tanca2
  real ( kind = 8 ) tas
  real ( kind = 8 ) tend
  real ( kind = 8 ) tf
  real ( kind = 8 ) tin
  real ( kind = 8 ) tinit
  character ( len = 25 )  title(ns)
  real ( kind = 8 ) tnow
  real ( kind = 8 ) tw
  real ( kind = 8 ) vave
  real ( kind = 8 ) vol(ni,nj)
  real ( kind = 8 ) x(ni,nj)
  real ( kind = 8 ) xbot(maxbot)
  real ( kind = 8 ) xc(ni,nj)
  real ( kind = 8 ) xlen
  real ( kind = 8 ) y(ni,nj)
  real ( kind = 8 ) ybot(maxbot)
  real ( kind = 8 ) yc(ni,nj)
  real ( kind = 8 ) ylen
!
!  Set things.
!
  ae1(1:ni,1:nj) = 0.0D+00
  ae2(1:ni,1:nj) = 0.0D+00
  ak1(1:ni,1:nj) = 0.0D+00
  ak2(1:ni,1:nj) = 0.0D+00
  area(1:ni,1:nj) = 0.0D+00
  b = 0.1778D+00
  cappa = 0.4187D+00
  cd = 1.0D+00
  ce1 = 1.44D+00
  ce2 = 1.92D+00
  cinc = 0.0D+00
  cl = 1000.0D+00
  cmu = 0.09D+00
  cost = 0.0D+00
  cs = 1000.0D+00
  cvn = 200000.0D+00
  dtm = 100.0D+00
  epsad = 0.0001D+00
  epsil = 0.0001D+00
  ewall = 9.793D+00
  f(1:ni,1:nj,1:ns) = 0.0D+00
  fkl = 64.0D+00
  fks = 22.0D+00
  fma = -1000.0D+00
  fn(1:ni,1:nj,1:ns) = 0.0D+00
  fnu = 3.0E-7
  fo(1:ni,1:nj,1:ns) = 0.0D+00
  fr = 1.0E-10
  grash = 10000000.0D+00
  grav = 9.81D+00
  hamag = 0.0D+00
  heta(1:ni,1:nj) = 0.0D+00
  hf = 1800000.0D+00
  hksi(1:ni,1:nj) = 0.0D+00
!
!  ICOST = 1, U**2+V**2
!    2, (T-TF)**2
!    3, (VAVE-SQRT(U**2+V**2))
!
  icost = 1
!
!  TEMPORARY
!
  icost = 3

  icrys = 42
  inturb = 0
  iplot = 1
  ipref = 11
!
!  IPRINT = 0, don't print very much.
!     1, print intermediate information.
!
  iprint = 0

  jcrys = 22
  jpref = 30
  l0 = 62
  l1 = 62
!
!  Interstep iteration number.
!  Use LAST = 40 for more accurate results.
!  Use LAST = 10 for quicker results.
!
  last = 40
!
!  Set the number of timesteps.
!
!  For tests, set LASTT = 2.
!  For a real ( kind = 8 ) run, try LASTT = 200.
!
  lastt = 2
  lblk(1:4) = .false.
  lblk(5:10) = .true.
  lortho = .false.
  lsolve(1:ns) = .false.
  m0 = 42
  m1 = 42
  mode = 1
  ndt = 0
  np = 3
  npc = 4
!
!  Number of linear iterations.
!
  nsolve(1) = 3
  nsolve(2) = 3
  nsolve(3) = 1
  nsolve(4) = 1
  nsolve(5) = 3
  nsolve(6) = 3
  nsolve(7) = 3
  nsolve(8) = 3
  nsolve(9) = 1
  nsolve(10) = 1
!
!  Number of nonlinear iterations.
!
  ntimes(1) = 5
  ntimes(2) = 5
  ntimes(3) = 1
  ntimes(4) = 1
  ntimes(5) = 3
  ntimes(6) = 3
  ntimes(7) = 3
  ntimes(8) = 3
  ntimes(9) = 1
  ntimes(10) = 1

  orth = 50000.0D+00
  pr = 0.015D+00
  rdtm = 0.0D+00
  re = 1.0D+00
  recb = 0.0D+00
  rect = 0.0D+00

  relax(1) = 0.3D+00
  relax(2) = 0.3D+00
  relax(3) = 0.8D+00
  relax(4) = 1.0D+00
  relax(5) = 0.7D+00
  relax(6) = 0.2D+00
  relax(7) = 0.5D+00
  relax(8) = 0.5D+00
  relax(9) = 1.0D+00
  relax(10) = 1.0D+00
  relax(11) = 1.0D+00
  relax(12) = 0.6D+00
  relax(13) = 1.0D+00
  relax(14) = 1.0D+00

  rhocon = 1.0D+00
  rhol = 2490.0D+00
  rhos = 2490.0D+00
  rueta(1:ni,1:nj) = 0.0D+00
  ruksi(1:ni,1:nj) = 0.0D+00
  sige = 1.3D+00
  sigk = 1.0D+00
  sigma = 0.72D+00
  sigt = 0.9D+00
  smooth = 0.0001D+00
  tal = 1523.0D+00
  tanca = 1.0D+00
  tanca2 = 1.0D+00
  tas = 1523.0D+00
  tf = 1683.0D+00
  tin = 1523.0D+00
  tinit = 0.0D+00

  title(1) = 'U velocity'
  title(2) = 'V velocity'
  title(3) = 'Pressure'
  title(4) = 'Corrected pressure'
  title(5) = 'Temperature'
  title(6) = 'Rotational velocity'
  title(7) = 'Turbulent dissipation'
  title(8) = 'Turbulent energy'
  title(9) = 'Magnetic stream function'
  title(10) = 'Stream function'
!
!  TEMPORARY
!
  tw = 1713.0D+00
  tw = par(1)

  vave = 1850.0D+00
  vol(1:ni,1:nj) = 0.0D+00
  x(1:ni,1:nj) = 0.0D+00
  xc(1:ni,1:nj) = 0.0D+00
  xlen = 1.2D+00
  y(1:ni,1:nj) = 0.0D+00
  yc(1:ni,1:nj) = 0.0D+00
  ylen = 1.0D+00
!
!  Set things that depend on things.
!
  birad = 0.25D+00 * 5.67D+00 * 1.5D+00 * 0.0254D+00 / fks /( tw - tf ) / re
  bo = (rhol*grav*b**2)/sigma
  cfo = fnu/(b*b)

  if ( inturb == 1 ) then
    f(1:ni,1:nj,7) = 0.005D+00
  end if

  fcsl = cs / cl
  fksl = fks/fkl
  frsl = rhos/rhol
  ra = grash*pr
  rho(1:ni,1:nj) = rhocon
  rpr = rhocon/(pr*re)
  stel = cl*(tw-tf)/hf
  stes = cs*(tf-tin)/hf
  tend = tinit+dtm*lastt
  tnow = tinit
!
!  Set things that depend on things that depend on things.
!
  delt = cfo*dtm

  if ( inturb == 1 ) then
    f(1:ni,1:nj,8) = f(1:ni,1:nj,7)**1.5D+00 / 0.006D+00
    gamt(1:ni,1:nj) = 0.006D+00 * cmu * rho(1:ni,1:nj) * f(1:ni,1:nj,7)**0.5D+00
  end if
!
!  Set the crucible shape.
!
  xbot(1) = 0.0D+00
  ybot(1) = 0.0D+00

  xbot(2) = 0.015D+00
  ybot(2) = 0.3D+00

  xbot(3) = 0.03D+00
  ybot(3) = 0.5D+00

  xbot(4) = 0.046D+00
  ybot(4) = 0.6D+00

  xbot(5) = 0.07D+00
  ybot(5) = 0.7D+00

  xbot(6) = 0.10D+00
  ybot(6) = 0.8D+00

  xbot(7) = 0.14D+00
  ybot(7) = 0.9D+00

  xbot(8) = 0.18D+00
  ybot(8) = 0.95D+00

  xbot(9) = 0.25D+00
  ybot(9) = 1.0D+00

  nbot = 9

  return
end
subroutine inigrd ( icrys, jcrys, l0, m0, nbot, xbot, xc, xlen, ybot, &
  yc, ylen )

!*****************************************************************************80
!
!! INIGRD makes an initial assignment of the grid points XC, YC.
!
  implicit none

  integer nbot
  integer, parameter :: ni = 64
  integer, parameter :: nj = 64

  real ( kind = 8 ) free
  integer i
  integer icrys
  integer j
  integer jcrys
  integer l0
  integer licrys
  integer ll1
  integer llst1
  integer m0
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) wave
  real ( kind = 8 ) xbot(nbot)
  real ( kind = 8 ) xc(ni,nj)
  real ( kind = 8 ) xlen
  real ( kind = 8 ) ybot(nbot)
  real ( kind = 8 ) yc(ni,nj)
  real ( kind = 8 ) ylen

  ll1 = icrys/2+1
  llst1 = jcrys/2+1
  licrys = (jcrys+m0)/2+1

  free = 0.5D+00 * xlen
!
!  Set the Y coordinates
!
  do i = 2, l0

    if ( i < ll1 ) then
      t2 = free * (0.5D+00*(dble(i-2)/dble(ll1-2))**1.5D+00)
    else if ( i < icrys ) then
      t2 = free * (1.0D+00 - 0.5D+00 * (dble(icrys-i)/dble(icrys-ll1))**1.5D+00)
    else
      t2 = free + (xlen-free)*(dble(i-icrys)/dble(l0-icrys))**1.2D+00
    end if

    do j = 2, m0

      wave = 0.4D+00 + 0.02D+00 * sin( t2 * 20.0D+00 * 3.14159D+00 )

      if ( j < llst1 ) then
        t1 = wave/2.0D+00 * (dble(j-2)/dble(llst1-2))**1.5D+00
      else if ( j < jcrys ) then
        t1 = wave-wave / 2.0D+00 * (dble(jcrys-j)/dble(jcrys-llst1))**1.5D+00
      else if ( j < licrys ) then
        t1 = 0.4D+00 + 0.3D+00 * (dble(j-jcrys)/dble(licrys-jcrys))**1.5D+00
      else
        t1 = 1.0D+00 - 0.3D+00 * (dble(m0-j)/dble(m0-licrys))**1.5D+00
      end if

      yc(i,j) = t1*ylen

    end do
  end do
!
!  Set the XC coordinates along the bottom boundary.
!
  xc(2,2) = xbot(1)

  do j = 3, m0-1

    call cubic(nbot,yc(2,j),ybot,xc(2,j),xbot)

  end do

  xc(2,m0) = xbot(nbot)
!
!  Now assign the XC coordinates of the other nodes.
!  Here, we explicitly assume that the free surface and crystal
!  boundaries occur at a coordinate value of 0.6.
!
  do j = 2, m0
    do i = 3, l0

      if ( i < ll1 ) then

        t1 = 0.5*(dble(i-2)/dble(ll1-2))**1.5
        xc(i,j) = (1.0-t1)*xc(2,j)+t1*free

      else if ( i <= icrys ) then

        t1 = 1.0D+00 - 0.5D+00 * (dble(icrys-i)/dble(icrys-ll1))**1.5D+00
        xc(i,j) = (1.0-t1)*xc(2,j)+t1*free

      else

        t1 = (dble(i-icrys)/dble(l0-icrys))**1.2D+00

        if ( j < jcrys ) then
          xc(i,j) = xc(icrys,j)+t1*((xlen-free)+0.5*yc(l0,jcrys)**2 &
            -0.5*yc(l0,j)**2)
        else
          xc(i,j) = xc(icrys,j)+t1*(xlen-free)
        end if

      end if

    end do
  end do
!
!  Copy values to dummy nodes with I = 1 or J=1.
!
  xc(1,1) = xc(2,2)
  yc(1,1) = yc(2,2)

  do j = 2, m0
    xc(1,j) = xc(2,j)
    yc(1,j) = yc(2,j)
  end do

  do i = 2, l0
    xc(i,1) = xc(i,2)
    yc(i,1) = yc(i,2)
  end do

  return
end
subroutine initl ( fjeta, fjksi, heta, hksi, l1, m1, ni, nj, p, rho, &
  rueta, ruksi, u, v, x, y )

!*****************************************************************************80
!
!! INITL estimates quantities at control volume interfaces.
!
!  Discussion:
!
!    INITL estimates the values of P, and of RHO*dU/dKSI and RHO*dU/dETA
!    at the interfaces of the control volumes by averaging values associated
!    with primary nodes.
!
  implicit none

  integer ni
  integer nj

  real ( kind = 8 ) fjeta(ni,nj)
  real ( kind = 8 ) fjksi(ni,nj)
  real ( kind = 8 ) heta(ni,nj)
  real ( kind = 8 ) hksi(ni,nj)
  integer i
  integer j
  integer l1
  integer m1
  real ( kind = 8 ) p(ni,nj)
  real ( kind = 8 ) rho(ni,nj)
  real ( kind = 8 ) rhov
  real ( kind = 8 ) rueta(ni,nj)
  real ( kind = 8 ) ruksi(ni,nj)
  real ( kind = 8 ) tt
  real ( kind = 8 ) u(ni,nj)
  real ( kind = 8 ) ucs
  real ( kind = 8 ) ulen
  real ( kind = 8 ) v(ni,nj)
  real ( kind = 8 ) vcs
  real ( kind = 8 ) x(ni,nj)
  real ( kind = 8 ) xcomp
  real ( kind = 8 ) y(ni,nj)
  real ( kind = 8 ) ycomp

  do i = 2, l1
    do j = 1, m1

      if ( hksi(i,j)+hksi(i-1,j) /= 0.0D+00 ) then
        tt = hksi(i,j)/(hksi(i,j)+hksi(i-1,j))
      else
        write ( *, * ) 'HKSI(I,J),HKSI(I-1,J) = 0, I=',i,' J=',j
        tt = 0.0D+00
      end if

      ucs = tt*u(i-1,j)+(1.0D+00-tt)*u(i,j)
      vcs = tt*v(i-1,j)+(1.0D+00-tt)*v(i,j)
      fjksi(i,j) = tt*p(i-1,j)+(1.0D+00-tt)*p(i,j)
      rhov = tt*rho(i-1,j)+(1.0D+00-tt)*rho(i,j)

      xcomp = x(i,j)-x(i-1,j)
      ycomp = y(i,j)-y(i-1,j)
      ulen = sqrt(xcomp**2+ycomp**2)

      if ( ulen > 0.0D+00 ) then
        ruksi(i,j) = rhov*(ucs*xcomp+vcs*ycomp)/ulen
      else
        ruksi(i,j) = 0.0D+00
      end if

    end do
  end do

  do i = 1, l1
    do j = 2, m1
 
      if ( heta(i,j)+heta(i,j-1) /= 0.0D+00 ) then
        tt = heta(i,j)/(heta(i,j)+heta(i,j-1))
      else
        tt = 0.0D+00
        write ( *, * ) 'HETA(I,J),HETA(I,J-1) = 0, I=',i,' J=',j
      end if

      ucs = tt*u(i,j-1)+(1.0-tt)*u(i,j)
      vcs = tt*v(i,j-1)+(1.0-tt)*v(i,j)
      fjeta(i,j) = tt*p(i,j-1)+(1.0-tt)*p(i,j)
      rhov = tt*rho(i,j-1)+(1.0-tt)*rho(i,j)

      xcomp = x(i,j)-x(i,j-1)
      ycomp = y(i,j)-y(i,j-1)
      ulen = sqrt(xcomp**2+ycomp**2)

      if ( ulen > 0.0D+00 ) then
        rueta(i,j) = rhov*(ucs*xcomp+vcs*ycomp)/ulen
      else
        rueta(i,j) = 0.0D+00
      end if

    end do
  end do

  return
end
subroutine movgrd ( b1jbl, b2jbl, bo, delt, fksl, fr, frsl, icrys, &
  iprint, jcrys, l0, m0, p, pr, re, stel, tanca, tanca2, xc, xlen, yc )

!*****************************************************************************80
!
!! MOVGRD calculates the new position of the interface and free surface.
!
  implicit none

  integer, parameter :: ni = 64
  integer, parameter :: nj = 64

  real ( kind = 8 ) ac
  real ( kind = 8 ) as
  real ( kind = 8 ) b1jbl(nj)
  real ( kind = 8 ) b2jbl(nj)
  real ( kind = 8 ) bo
  real ( kind = 8 ) coeff
  real ( kind = 8 ) coff
  real ( kind = 8 ) d2hdy2(nj)
  real ( kind = 8 ) delt
  real ( kind = 8 ) dhdy(nj)
  real ( kind = 8 ) dhdym(nj)
  real ( kind = 8 ) dhdyp(nj)
  real ( kind = 8 ) dlen
  real ( kind = 8 ) dx1
  real ( kind = 8 ) dxm1
  real ( kind = 8 ) dy1
  real ( kind = 8 ) dym1
  real ( kind = 8 ) fksl
  real ( kind = 8 ) flow(ni)
  real ( kind = 8 ) flomax
  real ( kind = 8 ) flomin
  real ( kind = 8 ) fr
  real ( kind = 8 ) frsl
  real ( kind = 8 ) high(ni)
  real ( kind = 8 ) hilev
  real ( kind = 8 ) himax
  real ( kind = 8 ) himin
  integer i
  integer icrys
  integer interl
  integer iprint
  integer j
  integer jcrys
  integer k
  integer l0
  integer m0
  real ( kind = 8 ) omega
  real ( kind = 8 ) omega2
  real ( kind = 8 ) p(ni,nj)
  real ( kind = 8 ) pcsurf
  real ( kind = 8 ) pr
  real ( kind = 8 ) pres(ni)
  real ( kind = 8 ) pullv
  real ( kind = 8 ) re
  real ( kind = 8 ) rr(nj)
  real ( kind = 8 ) stel
  real ( kind = 8 ) tanca
  real ( kind = 8 ) tanca2
  real ( kind = 8 ) xb(ni)
  real ( kind = 8 ) xc(ni,nj)
  real ( kind = 8 ) xlen
  real ( kind = 8 ) yc(ni,nj)
  real ( kind = 8 ) zb1(ni)

  do j = 2, m0
    zb1(j) = xc(icrys,j)
  end do

  do k = 1, 5

    do j = jcrys+1, m0-1
      dhdy(j) = (xc(icrys,j+1)-xc(icrys,j))/(yc(icrys,j+1)-yc(icrys,j))
      dhdym(j) = (xc(icrys,j)-xc(icrys,j-1))/(yc(icrys,j)-yc(icrys,j-1))
      dhdyp(j) = (xc(icrys,j+1)-xc(icrys,j))/(yc(icrys,j+1)-yc(icrys,j))
    end do

    do j = jcrys+1, m0-1
      d2hdy2(j) = (dhdyp(j)-dhdym(j))/(yc(icrys,j+1)-yc(icrys,j))
      rr(j) = d2hdy2(j)/(1.0+(dhdy(j)*dhdy(j)))**1.5 &
        +dhdy(j)/yc(icrys,j)/(1.0+(dhdy(j)*dhdy(j)))**0.5
    end do

    pres(jcrys+1) = 0.0D+00

    do j = jcrys+1, m0-2
      pres(j+1) = pres(j)+(p(icrys-1,j+1)-p(icrys-1,j))/yc(icrys-1,j+1)
    end do

    himax = 0.0D+00
    himin = 0.0D+00
    interl = (jcrys+m0)/2

    do j = jcrys+1, m0-1
      pcsurf = pres(j)-pres(interl)
      high(j) = fr*pcsurf+rr(j)/bo
    end do

    hilev = high(interl)

    omega = 0.0D+00

    do j = jcrys+1, m0-1
      high(j) = omega*((high(j)-hilev)-(xc(icrys,j)-xc(icrys,interl)))
      himax = max(himax,high(j))
      himin = min(himin,high(j))
      xc(icrys,j) = xc(icrys,j)+high(j)
    end do

    if ( omega /= 0.0D+00 ) then
      xc(icrys,jcrys) = xc(icrys,jcrys+1)+ &
        tanca*(yc(icrys,jcrys+1)-yc(icrys,jcrys))
      xc(icrys,m0) = xc(icrys,m0-1) + &
        tanca2*(yc(icrys,m0)-yc(icrys,m0-1))
    end if

    flomax = 0.0D+00
    flomin = 0.0D+00

    if ( iprint > 0 ) then
      if ( k == 5 ) then
        write ( *, * ) ' '
        write ( *, * ) 'MOVGRD:'
        write ( *, * ) '  Himin = ',himin
        write ( *, * ) '  Himax = ',himax
        write ( *, * ) '  Omega = ',omega
        write ( *, * ) ' '
        write(*,'(7f10.4)')high(jcrys+1),high(jcrys+2),high(m0-10), &
          high(m0-8),high(m0-4),high(m0-2),high(m0-1)
        write(*,'(7f9.4)')xc(icrys,jcrys),xc(icrys,jcrys+1), &
          xc(icrys,jcrys+4),xc(icrys,m0-8),xc(icrys,m0-4), &
          xc(icrys,m0-1),xc(icrys,m0)
        write ( *, * ) ' '
        write ( *, * ) 'MOVGRD:'
        write ( *, * ) '  Free surface movement'
      end if
    end if

  end do

  omega2 = 0.0D+00

  do j = 2, jcrys-1

    dxm1 = xc(icrys,j+1)-xc(icrys,j)
    dym1 = yc(icrys,j+1)-yc(icrys,j)
    dlen = sqrt(dxm1*dxm1+dym1*dym1)
    as = -dxm1/dlen
    ac = dym1/dlen
    dx1 = delt*frsl*(b1jbl(j)-fksl*b2jbl(j))*stel/(re*pr)*ac
    dy1 = delt*frsl*(b1jbl(j)-fksl*b2jbl(j))*stel/(re*pr)*as

    if ( j == jcrys-1 ) then
      flow(j) = dx1
    else
      flow(j) = dx1+dy1**2/dx1
    end if

  end do
!
!  What is PULLV?
!
  pullv = flow(jcrys-1)

  do j = 3, jcrys-1
    flomax = max(flomax,(flow(j)-pullv))
    flomin = min(flomin,(flow(j)-pullv))
  end do

  omega2 = min(omega2, 0.003D+00 / (flomax+1.0e-6))
  omega2 = min(omega2, 0.003D+00 / (-flomin+1.0e-6))

  do j = 3, jcrys-1
    flow(j) = omega2*(flow(j)-pullv)+(xc(icrys,jcrys)-zb1(jcrys))
    xb(j) = xc(icrys,j)+flow(j)
    xc(icrys,j) = xb(j)
  end do

  xc(icrys,2) = xc(icrys,3)

  if ( iprint > 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'MOVGRD:'
    write ( *, * ) '  Flomin  =  ',flomin
    write ( *, * ) '  Flomax  =  ',flomax
    write ( *, * ) '  Omega  =   ',omega2
    write ( *, * ) ' '
    write(*,'(2x,7f9.4)')xb(jcrys-1),xb(jcrys-2),xb(8),xb(6),xb(4),xb(3),xb(2)
    write(*,'(2x,7f9.4)')b1jbl(jcrys-1),-fksl*b2jbl(jcrys-1), &
      b1jbl(jcrys-2),-fksl*b2jbl(jcrys-2),b1jbl(jcrys-5),-fksl*b2jbl(jcrys-5)
    write(*,'(2x,7f9.4)')b1jbl(jcrys-1)-fksl*b2jbl(jcrys-1), &
      b1jbl(jcrys-2)-fksl*b2jbl(jcrys-2),b1jbl(jcrys-5)-fksl*b2jbl(jcrys-5), &
      b1jbl(jcrys-10)-fksl*b2jbl(jcrys-10),b1jbl(jcrys-15)-fksl*b2jbl(jcrys-15), &
      b1jbl(jcrys-17)-fksl*b2jbl(jcrys-17),b1jbl(jcrys-19)-fksl*b2jbl(jcrys-19)
    write(*,'(2x,7f9.4)')flow(21),flow(19),flow(17),flow(15), &
      flow(13),flow(11),flow(10)
    write(*,'(2x,7f9.4)')flow(9),flow(8),flow(7),flow(6),flow(5),flow(4),flow(3)

    write ( *, * ) ' '
    write ( *, * ) 'MOVGRD:'
    write ( *, * ) '  Solid-liquid interface movement'
    write ( *, * ) ' '
  end if

  do j = 2, m0

    coeff = (xc(2,j)-xc(icrys,j))/(xc(2,j)-zb1(j))
    coff = (xlen-xc(icrys,j))/(xlen-zb1(j))

    do i = 2, icrys-1
      xc(i,j) = xc(2,j)-(xc(2,j)-xc(i,j))*coeff
    end do

    do i = icrys+1, l0
      xc(i,j) = xlen-(xlen-xc(i,j))*coff
      if ( j <= jcrys .and. i == l0 ) then
        xc(l0,j) = xlen+0.5D+00*yc(l0,jcrys)**2-0.5D+00*yc(l0,j)**2
      end if
    end do

  end do

  return
end
subroutine output ( cfo, iter, izone, res, smax, ssum, t, tnow, u, w )

!*****************************************************************************80
!
!! OUTPUT prints information about the current solution.
!
  implicit none

  integer, parameter :: ni = 64
  integer, parameter :: nj = 64
  integer, parameter :: ns = 10

  real ( kind = 8 ) cfo
  integer iter
  integer izone
  real ( kind = 8 ) res(ns)
  real ( kind = 8 ) smax
  real ( kind = 8 ) ssum
  real ( kind = 8 ) t(ni,nj)
  real ( kind = 8 ) tnow
  real ( kind = 8 ) u(ni,nj)
  real ( kind = 8 ) w(ni,nj)
!
!  Solid zone output.
!
  if ( izone == 1 ) then

    if ( iter == 0 ) then
      write ( *, * ) ' '
      write ( *, * ) 'OUTPUT:'
      write ( *, * ) '  Current time  =  ',tnow
      write ( *, * ) '   =  ',tnow*cfo,' seconds.'
      write ( *, * ) '  Solid zone results.'
      write ( *, * ) '  Iter    Res(5)      T(50,21)    T(50,22)'
      write ( *, * ) ' '
    end if

    write(*,'(i4,2x,10g12.4)')iter,res(5),t(50,21),t(50,22)
!
!  Liquid zone output.
!
  else if ( izone == 2 ) then

    if ( iter == 1 ) then
      write ( *, * ) ' '
      write ( *, * ) ' '
      write ( *, * ) 'OUTPUT:'
      write ( *, * ) '  Liquid zone results.'
      write ( *, * ) '  SMAX is the maximum local mass imbalance.'
      write ( *, * ) '  SSUM is the total mass imbalance.'
      write ( *, * ) ' '
      write ( *, * ) ' Iter    SMAX        SSUM        U(5,5)      W(5,5)      T(5,5)'
      write ( *, * ) ' '
    end if

    write(*,'(i4,2x,10g12.4)')iter,smax,ssum,u(5,5),w(5,5),t(5,5)

  end if

  return
end
subroutine pmod ( ipref, jcrys, jpref, l1, m1, p )

!*****************************************************************************80
!
!! PMOD extends pressure values from primary nodes to corner nodes.
!
!  Discussion:
!
!    PMOD carries out some simple operations to extend the value of the
!    pressure to primary nodes at the corners, and to normalize the
!    pressure by subtracting off its value at the reference point.
!
  implicit none

  integer, parameter :: ni = 64
  integer, parameter :: nj = 64

  integer i
  integer ipref
  integer j
  integer jcrys
  integer jpref
  integer l1
  integer m1
  real ( kind = 8 ) p(ni,nj)
  real ( kind = 8 ) pref
!
!  Extrapolate to get pressures on the boundary corners.
!
  p(1,1) = p(2,1)+p(1,2)-p(2,2)
  p(l1,1) = p(l1-1,1)+p(l1,2)-p(l1-1,2)
  p(1,jcrys) = p(2,jcrys)+p(1,jcrys-1)-p(2,jcrys-1)
  p(l1,jcrys) = p(l1-1,jcrys)+p(l1,jcrys-1)-p(l1-1,jcrys-1)
!
!  Subtract off the reference pressure.
!
  pref = p(ipref,jpref)

  do j = 1, m1
    do i = 1, l1
      p(i,j) = p(i,j) - pref
    end do
  end do

  return
end
subroutine prdat ( b, birad, bo, cappa, cfo, cvn, delt, dtm, fcsl, fksl, &
  fma, fnu, fr, frsl, grash, icost, icrys, inturb, iprint, &
  jcrys, l0, last, lastt, m0, mode, nbot, ns, nsolve, ntimes, orth, &
  pr, ra, rdtm, recb, rect, rhocon, smooth, stel, stes, tanca, tanca2, tend, &
  tf, tinit, title, tw, vave, xlen, ylen )

!*****************************************************************************80
!
!! PRDAT prints out the initial values of certain data.
!
  implicit none

  integer ns

  real ( kind = 8 ) b
  real ( kind = 8 ) birad
  real ( kind = 8 ) bo
  real ( kind = 8 ) cappa
  real ( kind = 8 ) cfo
  real ( kind = 8 ) cvn
  real ( kind = 8 ) delt
  real ( kind = 8 ) dtm
  real ( kind = 8 ) fcsl
  real ( kind = 8 ) fksl
  real ( kind = 8 ) fma
  real ( kind = 8 ) fnu
  real ( kind = 8 ) fr
  real ( kind = 8 ) frsl
  real ( kind = 8 ) grash
  integer i
  integer icost
  integer icrys
  integer inturb
  integer iprint
  integer jcrys
  integer l0
  integer last
  integer lastt
  integer m0
  integer mode
  integer nbot
  integer nsolve(ns)
  integer ntimes(ns)
  real ( kind = 8 ) orth
  real ( kind = 8 ) pr
  real ( kind = 8 ) ra
  real ( kind = 8 ) rdtm
  real ( kind = 8 ) recb
  real ( kind = 8 ) rect
  real ( kind = 8 ) rhocon
  real ( kind = 8 ) smooth
  real ( kind = 8 ) stel
  real ( kind = 8 ) stes
  real ( kind = 8 ) tanca
  real ( kind = 8 ) tanca2
  real ( kind = 8 ) tend
  real ( kind = 8 ) tf
  real ( kind = 8 ) tinit
  character ( len = 25 ) title(ns)
  real ( kind = 8 ) tw
  real ( kind = 8 ) vave
  real ( kind = 8 ) xlen
  real ( kind = 8 ) ylen

  write ( *, * ) ' '
  write ( *, * ) 'PRDAT:'
  write ( *, * ) ' '
  write ( *, * ) '  B, Crucible radius  =                 ',b
  write ( *, * ) '  BIRAD, thermal coefficient  =         ',birad
  write ( *, * ) '  BO, the Bond number  =                ',bo
  write ( *, * ) '  CAPPA, wall function constant  =      ',cappa
  write ( *, * ) '  CFO, time scale  =                    ',cfo
  write ( *, * ) '  CVN, ADAPT cell volume weight  =      ',cvn
  write ( *, * ) '  DELT, time step in seconds  =         ',delt
  write ( *, * ) '  DTM, dimensionless time step  =       ',dtm
  write ( *, * ) '  FCSL  =  CS/CL =                      ',fcsl
  write ( *, * ) '  FKSL  =  FKS/FKL =                    ',fksl
  write ( *, * ) '  FMA, Marangoni number FMA  =          ',fma
  write ( *, * ) '  FR, Froude number  =                  ',fr
  write ( *, * ) '  FRSL  =  RHOS/RHOL =                  ',frsl
  write ( *, * ) '  GRASH, Grashof number  =              ',grash
  write ( *, * ) '  ICOST, cost functional  =             ',icost
  write ( *, * ) '  ICRYS, maximum I of crystal  =        ',icrys
  write ( *, * ) '  INTURB, turbulence option  =          ',inturb
  write ( *, * ) '  IPRINT, printing option  =            ',iprint
  write ( *, * ) '  JCRYS, maximum J of crystal  =        ',jcrys
  write ( *, * ) '  L0, number of I nodes  =              ',l0
  write ( *, * ) '  LAST, number of zone iterations on  '
  write ( *, * ) '    each time step  =                   ',last
  write ( *, * ) '  LASTT, number of time steps        =  ',lastt
  write ( *, * ) '  M0, number of J nodes  =              ',m0
  write ( *, * ) '  MODE, 0 cartesian, 1 axisymmetric  =  ',mode
  write ( *, * ) '  NBOT, number of boundary points  =    ',nbot
  write ( *, * ) '  ORTH, ADAPT orthogonality weight  =   ',orth
  write ( *, * ) '  PR, Prandtl number PR  =              ',pr
  write ( *, * ) '  RA, Rayleigh number  =                ',ra
  write ( *, * ) '  RBD  =  FMA/RA =                      ',fma/ra
  write ( *, * ) '  RDTM,  =  1/DTM or 0 =                ',rdtm
  write ( *, * ) '  RECB, crucible Reynolds number  =     ',recb
  write ( *, * ) '  RECT, crystal Reynolds number  =      ',rect
  write ( *, * ) '  RHOCON, density constant  =           ',rhocon
  write ( *, * ) '  SMOOTH, ADAPT smoothness weight  =    ',smooth
  write ( *, * ) '  STEL, liquid Stefan number STEL  =    ',stel
  write ( *, * ) '  STES, solid Stefan number STES  =     ',stes
  write ( *, * ) '  TANCA  =                              ',tanca
  write ( *, * ) '  TANCA2  =                             ',tanca2
  write ( *, * ) '  TEND, dimensionless end time  =       ',tend
  write ( *, * ) '  TF, crystal melting temperature  =    ',tf
  write ( *, * ) '  TINIT, dimensionless start time  =    ',tinit
  write ( *, * ) '  TW, the wall temperature  =           ',tw
  write ( *, * ) '  VAVE, desired average velocity  =     ',vave
  write ( *, * ) '  WEBER  =  FR*BO =                     ',fr*bo
  write ( *, * ) '  XLEN  =  problem region length =      ',xlen
  write ( *, * ) '  YLEN  =  problem region height =      ',ylen
  write ( *, * ) '  Characteristic velocity FNU/B  =      ',fnu/b
  write ( *, * ) ' '
  write ( *, * ) 'Variable                     Nonlinear    Linear'
  write ( *, * ) '                             Iterations   Iterations'
  write ( *, * ) ' '
  do i = 1, ns
    write ( *, * ) title(i),ntimes(i),nsolve(i)
  end do

  return
end
subroutine resid ( aim, aip, ajm, ajp, ap, con, f, l1, m1, nf )

!*****************************************************************************80
!
!! RESID computes the linear equation residual at interior primary nodes.
!
  implicit none

  integer, parameter :: ni = 64
  integer, parameter :: nj = 64
  integer, parameter :: ns = 10

  real ( kind = 8 ) aim(ni,nj)
  real ( kind = 8 ) aip(ni,nj)
  real ( kind = 8 ) ajm(ni,nj)
  real ( kind = 8 ) ajp(ni,nj)
  real ( kind = 8 ) ap(ni,nj)
  real ( kind = 8 ) con(ni,nj)
  real ( kind = 8 ) conmax
  real ( kind = 8 ) diamax
  real ( kind = 8 ) f(ni,nj,ns)
  integer i
  integer imax
  integer iprint
  integer j
  integer jmax
  integer l1
  integer m1
  integer nf
  real ( kind = 8 ) offmax
  real ( kind = 8 ) ratio
  real ( kind = 8 ) ratmin
  real ( kind = 8 ) res
  real ( kind = 8 ) resmax
  real ( kind = 8 ) sum
  real ( kind = 8 ) xmax

  iprint = 0
  conmax = 0.0D+00
  diamax = 0.0D+00
  offmax = 0.0D+00
  ratmin = 100000.0D+00
  resmax = 0.0D+00
  xmax = 0.0D+00

  imax = 0
  jmax = 0

  do i = 2, l1-1
    do j = 2, m1-1

      res = con(i,j)-ap(i,j)*f(i,j,nf)+aim(i,j)*f(i-1,j,nf) &
        +aip(i,j)*f(i+1,j,nf)+ajm(i,j)*f(i,j-1,nf)+ajp(i,j)*f(i,j+1,nf)

      if ( abs(res) >= resmax ) then
        resmax = abs(res)
        imax = i
        jmax = j
      end if

      if ( abs(con(i,j)) >= conmax ) then
        conmax = abs(con(i,j))
      end if

      if ( abs(f(i,j,nf)) >= xmax ) then
        xmax = abs(f(i,j,nf))
      end if

      if ( abs(ap(i,j)) >= diamax ) then
        diamax = abs(ap(i,j))
      end if

      sum = abs(aim(i,j))+abs(aip(i,j))+abs(ajm(i,j))+abs(ajp(i,j))
      if ( sum >= offmax ) then
        offmax = sum
      end if

      if ( sum /= 0.0D+00 ) then
        ratio = abs(ap(i,j))/sum
        if ( ratio < ratmin ) then
          ratmin = ratio
        end if
      end if

    end do
  end do

  if ( iprint == 1 ) then
    write ( *, * ) ' '
    write ( *, * ) 'RESID'
    write ( *, * ) '  Maxixum residual for variable ',nf,' is ',resmax
    write ( *, * ) '  at I,J = ',imax,jmax
    i = imax
    j = jmax
    write ( *, * ) ' '
    write ( *, * ) '  A(I,J)*X(I,J) =     ',ap(i,j),f(i,j,nf)
    write ( *, * ) '  A(I-1,J)*X(I-1,J) = ',aim(i,j),f(i-1,j,nf)
    write ( *, * ) '  A(I+1,J)*X(I+1,J) = ',aip(i,j),f(i+1,j,nf)
    write ( *, * ) '  A(I,J-1)*X(I,J-1) = ',ajm(i,j),f(i,j-1,nf)
    write ( *, * ) '  A(I,J+1)*X(I,J+1) = ',ajp(i,j),f(i,j+1,nf)
    write ( *, * ) '  CON(I,J) =          ',con(i,j)
    write ( *, * ) ' '
    write ( *, * ) '  Norm of variable is ',xmax
    write ( *, * ) '  Norm of RHS is      ',conmax
    write ( *, * ) '  Norm of diagonal is ',diamax
    write ( *, * ) '  Norm of off-diag is ',offmax
    write ( *, * ) '  Min diag/off-diag   ',ratmin
  end if

  return
end
subroutine rswrit ( cost, e, gamt, icrys, jcrys, l0, m0, nbot, p, pc, &
  psi, rueta, ruksi, t, te, tk, tnow, u, v, w, x, xbot, xc, y, ybot, yc )

!*****************************************************************************80
!
!! RSWRIT writes out restart information.
!
!  Modified:
!
!    22 March 2002
!
  implicit none

  integer nbot
  integer, parameter :: ni = 64
  integer, parameter :: nj = 64

  real ( kind = 8 ) cost
  real ( kind = 8 ) e(ni,nj)
  character ( len = 80 ) :: file_name = 'crystal_rs.txt'
  real ( kind = 8 ) gamt(ni,nj)
  integer i
  integer icrys
  integer j
  integer jcrys
  integer l0
  integer m0
  real ( kind = 8 ) p(ni,nj)
  real ( kind = 8 ) pc(ni,nj)
  real ( kind = 8 ) psi(ni,nj)
  real ( kind = 8 ) rueta(ni,nj)
  real ( kind = 8 ) ruksi(ni,nj)
  real ( kind = 8 ) t(ni,nj)
  real ( kind = 8 ) te(ni,nj)
  real ( kind = 8 ) tk(ni,nj)
  real ( kind = 8 ) tnow
  real ( kind = 8 ) u(ni,nj)
  real ( kind = 8 ) v(ni,nj)
  real ( kind = 8 ) w(ni,nj)
  real ( kind = 8 ) x(ni,nj)
  real ( kind = 8 ) xbot(nbot)
  real ( kind = 8 ) xc(ni,nj)
  real ( kind = 8 ) y(ni,nj)
  real ( kind = 8 ) ybot(nbot)
  real ( kind = 8 ) yc(ni,nj)

  write ( *, '(a)' )' '
  write ( *, '(a)' )'RSWRIT:'
  write ( *, '(a)' )'  Writing restart information to ' // trim ( file_name )
  write ( *, '(a)' )' '

  open ( unit = 10, file = file_name, status = 'replace' )

  write ( 10, '(e12.5)' ) cost
  write ( 10, '(i6)' ) l0
  write ( 10, '(i6)' ) jcrys
  write ( 10, '(i6)' ) icrys
  write ( 10, '(i6)' ) m0
  write ( 10, '(i6)' ) nbot
  write ( 10, '(e12.5)' ) tnow

  write ( 10, '(6e12.5)' ) ((e(i,j),i = 1,l0),j=1,m0)
  write ( 10, '(6e12.5)' ) ((gamt(i,j),i = 1,l0),j=1,m0)
  write ( 10, '(6e12.5)' ) ((p(i,j),i = 1,l0),j=1,m0)
  write ( 10, '(6e12.5)' ) ((pc(i,j),i = 1,l0),j=1,m0)
  write ( 10, '(6e12.5)' ) ((psi(i,j),i = 1,l0),j=1,m0)
  write ( 10, '(6e12.5)' ) ((rueta(i,j),i = 1,l0),j=1,m0)
  write ( 10, '(6e12.5)' ) ((ruksi(i,j),i = 1,l0),j=1,m0)
  write ( 10, '(6e12.5)' ) ((t(i,j),i = 1,l0),j=1,m0)
  write ( 10, '(6e12.5)' ) ((te(i,j),i = 1,l0),j=1,m0)
  write ( 10, '(6e12.5)' ) ((tk(i,j),i = 1,l0),j=1,m0)
  write ( 10, '(6e12.5)' ) ((u(i,j),i = 1,l0),j=1,m0)
  write ( 10, '(6e12.5)' ) ((v(i,j),i = 1,l0),j=1,m0)
  write ( 10, '(6e12.5)' ) ((w(i,j),i = 1,l0),j=1,m0)
  write ( 10, '(6e12.5)' ) ((x(i,j),i = 1,l0),j=1,m0)
  write ( 10, '(6e12.5)' ) xbot(1:nbot)
  write ( 10, '(6e12.5)' ) ((xc(i,j),i = 1,l0),j=1,m0)
  write ( 10, '(6e12.5)' ) ((y(i,j),i = 1,l0),j=1,m0)
  write ( 10, '(6e12.5)' ) ybot(1:nbot)
  write ( 10, '(6e12.5)' ) ((yc(i,j),i = 1,l0),j=1,m0)

  close ( unit = 10 )

  return
end
subroutine setcst ( area, cinc, icost, icrys, m0, ni, nj, t, tf, u, &
  v, vave )

!*****************************************************************************80
!
!! SETCST computes a portion of the cost functional.
!
!  Discussion:
!
!    SETCST computes a portion on the cost functional, which is currently
!    the integral over time and space of the norm of the fluid velocity.
!
!    The integral to be computed for a particular time is:
!
!      ICOST = 1:
!
!        F(T)  =  SQRT( Integral ( (X,Y) in LIQUID(T))
!               (U(T,X,Y)**2 + V(T,X,Y)**2) dX dY )
!
!      ICOST = 2:
!
!        F(T)  =  SQRT( Integral ( (X,Y) in LIQUID(T))
!               (T(T,X,Y)-TF)**2 dX dY )
!
!      ICOST = 3:
!
!        F(T)  =  SQRT( Integral ( (X,Y) in LIQUID(T))
!               ( Vave - SQRT(U(T,X,Y)**2 + V(T,X,Y)**2)) dX dY )
!
!    and the total cost is
!
!      COST  =  Integral (T=TINIT to TEND) F(T) dT/ (TEND-TINIT)
!
!    Currently, a crude approximation is made to the integrals.
!    For instance, the area of the flow region, LIQUID(T), is computed
!    by summing up the estimated areas of each control volume.  These
!    areas, in turn, are estimated by computing the jacobian at the
!    center of the control volume.
!
!    +-----V-----+
!    |           |
!    |           |
!    U     N     U
!    |           |
!    |           |
!    +-----V-----+
!
!  Modified:
!
!    22 March 2002
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ni
  integer nj

  real ( kind = 8 ) area(ni,nj)
  real ( kind = 8 ) cinc
  integer i
  integer icost
  integer icrys
  integer j
  integer m0
  real ( kind = 8 ) t(ni,nj)
  real ( kind = 8 ) tf
  real ( kind = 8 ) u(ni,nj)
  real ( kind = 8 ) v(ni,nj)
  real ( kind = 8 ) vave

  if ( icost == 1 ) then

    cinc = 0.0D+00
    do i = 1, icrys-1
      do j = 1, m0
        cinc = cinc+(u(i,j)**2+v(i,j)**2)*area(i,j)
      end do
    end do

    cinc = sqrt(cinc)

  else if ( icost == 2 ) then

    cinc = 0.0D+00
    do i = 1, icrys-1
      do j = 1, m0
        cinc = cinc+(t(i,j)-tf)**2*area(i,j)
      end do
    end do

    cinc = sqrt(cinc)

  else if ( icost == 3 ) then

    cinc = 0.0D+00
    do i = 1, icrys-1
      do j = 1, m0
        cinc = cinc+area(i,j)*(vave-sqrt(u(i,j)**2+v(i,j)**2))**2
      end do
    end do

    cinc = sqrt(cinc)

  end if

  return
end
subroutine setgeo ( ae1, ae2, ak1, ak2, heta, hksi, l1, m1, mode, ni, &
  nj, r, vol, x, xc, y, yc )

!*****************************************************************************80
!
!! SETGEO calculates various geometric quantities.
!
!  Discussion:
!
!    SETGEO is given (XC,YC), the locations of the "corners" of the
!    control volumes, and calculates various related geometric parameters,
!    including:
!
!      X, Y the position of the primary nodes,
!      HETA and HKSI, dH/dETA and dH/dKSI,
!      the Jacobian,
!      AK1, AE1, dALPHA/dKSI, dALPHA/dETA,
!      AK2, AE2, dBETA/dKSI, dBETA/dETA,
!      R, a weight factor for cylindrical geometries.
!      VOL, the volume of the control volumes.
!
  implicit none

  integer ni
  integer nj

  real ( kind = 8 ) ae1(ni,nj)
  real ( kind = 8 ) ae2(ni,nj)
  real ( kind = 8 ) ak1(ni,nj)
  real ( kind = 8 ) ak2(ni,nj)
  real ( kind = 8 ) dxdeta
  real ( kind = 8 ) dxdksi
  real ( kind = 8 ) dydeta
  real ( kind = 8 ) dydksi
  real ( kind = 8 ) heta(ni,nj)
  real ( kind = 8 ) hksi(ni,nj)
  integer i
  integer j
  integer l1
  integer m1
  integer mode
  real ( kind = 8 ) r(ni,nj)
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) t3
  real ( kind = 8 ) vol(ni,nj)
  real ( kind = 8 ) x(ni,nj)
  real ( kind = 8 ) xc(ni,nj)
  real ( kind = 8 ) xd
  real ( kind = 8 ) xjacb
  real ( kind = 8 ) xm
  real ( kind = 8 ) xp
  real ( kind = 8 ) xu
  real ( kind = 8 ) y(ni,nj)
  real ( kind = 8 ) yc(ni,nj)
  real ( kind = 8 ) yd
  real ( kind = 8 ) ym
  real ( kind = 8 ) yp
  real ( kind = 8 ) yu
!
!  Calculate the local scale factors HETA and HKSI, and
!  the control volumes VOL.

  do i = 2, l1-1
    do j = 2, m1-1
!
!  Compute the coordinates of the mid-side nodes.
!
!    ^
!    |     [I,J+1]  Up        [I+1,J+1]
!    E
!    T     Minus    (I,J)      Plus
!    A
!    |     [I,J]    Down       [I+1,J]
!    |
!    +-----XSI---->
!
      xp = 0.5D+00 * (xc(i+1,j+1)+xc(i+1,j))
      yp = 0.5D+00 * (yc(i+1,j+1)+yc(i+1,j))

      xm = 0.5D+00 * (xc(i,j+1)+xc(i,j))
      ym = 0.5D+00 * (yc(i,j+1)+yc(i,j))

      xu = 0.5D+00 * (xc(i,j+1)+xc(i+1,j+1))
      yu = 0.5D+00 * (yc(i,j+1)+yc(i+1,j+1))

      xd = 0.5D+00 * (xc(i,j)+xc(i+1,j))
      yd = 0.5D+00 * (yc(i,j)+yc(i+1,j))
!
!  The mid-side nodes differ by Delta KSI or Delta ETA  =  1.
!  Use finite differences to estimate derivatives like dX/dKSI.
!
      dxdksi = xp - xm
      dxdeta = xu - xd
      dydksi = yp - ym
      dydeta = yu - yd
!
!  Now compute the length of the line segments that connect
!  opposing mid-side nodes.
!
      hksi(i,j) = sqrt(dxdksi**2+dydksi**2)
      heta(i,j) = sqrt(dxdeta**2+dydeta**2)
!
!  Now compute the area of the control volume, which is just
!  the determinant of the Jacobian matrix.
!
      vol(i,j) = dxdksi*dydeta-dxdeta*dydksi

    end do
  end do
!
!  Take care of values along the borders J = 1 and J=M1.
!
  do i = 1, l1

    heta(i,1) = 0.0D+00
    heta(i,m1) = 0.0D+00
    vol(i,1) = 0.0D+00
    vol(i,m1) = 0.0D+00

    if ( i /= 1 .and. i /= l1 ) then
      hksi(i,1) = sqrt((xc(i+1,2)-xc(i,2))**2+(yc(i+1,2)-yc(i,2))**2)
      hksi(i,m1) = sqrt((xc(i+1,m1)-xc(i,m1))**2+(yc(i+1,m1)-yc(i,m1))**2)
    else
      hksi(i,1) = 0.0D+00
      hksi(i,m1) = 0.0D+00
    end if

  end do
!
!  Take care of values along the borders I = 1 and I=L1.
!
  do j = 1, m1

    hksi(1,j) = 0.0D+00
    hksi(l1,j) = 0.0D+00
    vol(1,j) = 0.0D+00
    vol(l1,j) = 0.0D+00

    if ( j /= 1 .and. j /= m1 ) then
      heta(1,j) = sqrt((xc(2,j)-xc(2,j+1))**2+(yc(2,j)-yc(2,j+1))**2)
      heta(l1,j) = sqrt((xc(l1,j)-xc(l1,j+1))**2+(yc(l1,j)-yc(l1,j+1))**2)
    else
      heta(1,j) = 0.0D+00
      heta(l1,j) = 0.0D+00
    end if

  end do
!
!  Calculate dALPHA/dKSI, dBETA/dKSI, dALPHA/dETA, dBETA/dETA,
!  the areas on the control-volume faces.
!
  do j = 2, m1-1
    do i = 2, l1

      dxdeta = xc(i,j+1)-xc(i,j)
      dydeta = yc(i,j+1)-yc(i,j)
      dxdksi = x(i,j)-x(i-1,j)
      dydksi = y(i,j)-y(i-1,j)

      if ( i == 2 .or. i == l1 ) then
        dxdksi = dxdksi * 2.0D+00
        dydksi = dydksi * 2.0D+00
      end if

      t1 = dxdeta**2+dydeta**2
      t2 = dxdksi**2+dydksi**2
      t3 = dxdksi*dxdeta+dydeta*dydksi

      xjacb = dxdksi*dydeta-dxdeta*dydksi

      if ( xjacb /= 0.0D+00 ) then
        ak1(i,j) = sqrt(t2)*t1/xjacb
        ak2(i,j) = sqrt(t1)*t3/xjacb
      else
        ak1(i,j) = 0.0D+00
        ak2(i,j) = 0.0D+00
      end if

    end do
  end do

  do i = 2, l1-1
    do j = 2, m1

      dxdeta = x(i,j)-x(i,j-1)
      dydeta = y(i,j)-y(i,j-1)
      dxdksi = xc(i+1,j)-xc(i,j)
      dydksi = yc(i+1,j)-yc(i,j)

      if ( j == 2 .or. j == m1 ) then
        dxdeta = dxdeta*2.0D+00
        dydeta = dydeta*2.0D+00
      end if

      t1 = dxdeta**2+dydeta**2
      t2 = dxdksi**2+dydksi**2
      t3 = dxdksi*dxdeta+dydeta*dydksi

      xjacb = dxdksi*dydeta-dxdeta*dydksi

      if ( xjacb /= 0.0D+00 ) then
        ae1(i,j) = sqrt(t1)*t2/xjacb
        ae2(i,j) = sqrt(t2)*t3/xjacb
      else
        ae1(i,j) = 0.0D+00
        ae2(i,j) = 0.0D+00
      end if

    end do
  end do
!
!  Set the cylindrical geometry weight factor.
!
  if ( mode == 0 ) then

    r(1:l1,1:m1) = 1.0D+00

  else if ( mode == 1 ) then

    do i = 1, l1
      do j = 1, m1

        r(i,j) = y(i,j)

        ak1(i,j) = ak1(i,j)*r(i,j)
        ak2(i,j) = ak2(i,j)*r(i,j)
        ae1(i,j) = ae1(i,j)*r(i,j)
        ae2(i,j) = ae2(i,j)*r(i,j)
        vol(i,j) = vol(i,j)*r(i,j)

      end do
    end do

  end if

  return
end
subroutine setup ( ae1, ae2, ak1, ak2, b1jbl, b2jbl, birad, cappa, &
  cd,ce1,ce2,cmu,epsil,ewall,f,fcsl,fjeta,fjksi,fksl,fma, &
  fmax,fn,fo,frsl,gamt,grash,hamag,heta,hksi,inturb,icrys, &
  iter,izone,jcrys,l1,lblk,lconv,lortho,lsolve,m1, &
  mode,nf,np,npc,nsolve,ntimes,pr,r,rdtm,re,recb,rect,res, &
  relax,rho,rhocon,rpr,rueta,ruksi,sige,sigk,sigt,smax,ssum, &
  stel, stes, tal, tas, tf, tw, vol, x, xc, y, yc )

!*****************************************************************************80
!
!! SETUP calculates the coefficients of the equations, and solves them.
!
  implicit none

  integer, parameter :: ni = 64
  integer, parameter :: nj = 64
  integer, parameter :: nk = 14
  integer, parameter :: nmaxij = 64
  integer, parameter :: ns = 10

  real ( kind = 8 ) a1
  real ( kind = 8 ) a2
  real ( kind = 8 ) acof
  real ( kind = 8 ) ae1(ni,nj)
  real ( kind = 8 ) ae2(ni,nj)
  real ( kind = 8 ) aim(ni,nj)
  real ( kind = 8 ) aip(ni,nj)
  real ( kind = 8 ) ajm(ni,nj)
  real ( kind = 8 ) ajp(ni,nj)
  real ( kind = 8 ) ak1(ni,nj)
  real ( kind = 8 ) ak2(ni,nj)
  real ( kind = 8 ) ap(ni,nj)
  real ( kind = 8 ) apr
  real ( kind = 8 ) b1jbl(nj)
  real ( kind = 8 ) b2jbl(nj)
  real ( kind = 8 ) birad
  real ( kind = 8 ) biu(nj)
  real ( kind = 8 ) biv(nj)
  real ( kind = 8 ) bju(ni)
  real ( kind = 8 ) bjv(ni)
  real ( kind = 8 ) blu(nj)
  real ( kind = 8 ) blv(nj)
  real ( kind = 8 ) bmu(ni)
  real ( kind = 8 ) bmv(ni)
  real ( kind = 8 ) bpi
  real ( kind = 8 ) bpi1
  real ( kind = 8 ) bpj
  real ( kind = 8 ) bpj1
  real ( kind = 8 ) bpm
  real ( kind = 8 ) cappa
  real ( kind = 8 ) cd
  real ( kind = 8 ) ce1
  real ( kind = 8 ) ce2
  real ( kind = 8 ) cmu
  real ( kind = 8 ) cofu(ni,nj,5)
  real ( kind = 8 ) cofv(ni,nj,5)
  real ( kind = 8 ) con(ni,nj)
  real ( kind = 8 ) conu(ni,nj)
  real ( kind = 8 ) conv(ni,nj)
  real ( kind = 8 ) denom
  real ( kind = 8 ) diff
  real ( kind = 8 ) dpeta
  real ( kind = 8 ) dpksi
  real ( kind = 8 ) dxdeta
  real ( kind = 8 ) dxdksi
  real ( kind = 8 ) dydeta
  real ( kind = 8 ) dydksi
  real ( kind = 8 ) em
  real ( kind = 8 ) ep
  real ( kind = 8 ) epsil
  real ( kind = 8 ) ewall
  real ( kind = 8 ) f(ni,nj,ns)
  real ( kind = 8 ) fcsl
  real ( kind = 8 ) fjeta(ni,nj)
  real ( kind = 8 ) fjksi(ni,nj)
  real ( kind = 8 ) fksl
  real ( kind = 8 ) flow
  real ( kind = 8 ) fma
  real ( kind = 8 ) fmax(ns)
  real ( kind = 8 ) fn(ni,nj,ns)
  real ( kind = 8 ) fo(ni,ni,ns)
  real ( kind = 8 ) frc
  real ( kind = 8 ) frsl
  real ( kind = 8 ) gam(ni,nj)
  real ( kind = 8 ) gamt(ni,nj)
  real ( kind = 8 ) grash
  real ( kind = 8 ) hamag
  real ( kind = 8 ) heta(ni,nj)
  real ( kind = 8 ) hetap(ni,nj)
  real ( kind = 8 ) hksi(ni,nj)
  real ( kind = 8 ) hksip(ni,nj)
  real ( kind = 8 ) hutop(ni,nj)
  real ( kind = 8 ) hvtop(ni,nj)
  integer i
  integer icrys
  integer inturb
  integer isol
  integer iter
  integer izone
  integer j
  integer jcrys
  integer l1
  logical lblk(ns)
  logical lconv
  logical lortho
  logical lsolve(ns)
  integer m1
  integer mode
  integer n
  integer nf
  integer np
  integer npc
  integer nsolve(ns)
  integer ntimes(ns)
  real ( kind = 8 ) peta(ni,nj)
  real ( kind = 8 ) pksi(ni,nj)
  real ( kind = 8 ) pr
  real ( kind = 8 ) qt(nmaxij)
  real ( kind = 8 ) r(ni,nj)
  real ( kind = 8 ) rdtm
  real ( kind = 8 ) re
  real ( kind = 8 ) recb
  real ( kind = 8 ) rect
  real ( kind = 8 ) relax(nk)
  real ( kind = 8 ) res(ns)
  real ( kind = 8 ) rho(ni,nj)
  real ( kind = 8 ) rhocon
  real ( kind = 8 ) rpr
  real ( kind = 8 ) rueij
  real ( kind = 8 ) rueij1
  real ( kind = 8 ) ruet
  real ( kind = 8 ) rueta(ni,nj)
  real ( kind = 8 ) ruki1j
  real ( kind = 8 ) rukij
  real ( kind = 8 ) ruksi(ni,nj)
  real ( kind = 8 ) rukt
  real ( kind = 8 ) sige
  real ( kind = 8 ) sigk
  real ( kind = 8 ) sigt
  real ( kind = 8 ) smax
  real ( kind = 8 ) soor
  real ( kind = 8 ) ssum
  real ( kind = 8 ) stel
  real ( kind = 8 ) stes
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) t3
  real ( kind = 8 ) t4
  real ( kind = 8 ) tal
  real ( kind = 8 ) tas
  real ( kind = 8 ) tem1
  real ( kind = 8 ) tem2
  real ( kind = 8 ) temp
  real ( kind = 8 ) tf
  real ( kind = 8 ) tmp(ni,nj)
  real ( kind = 8 ) tw
  real ( kind = 8 ) vol(ni,nj)
  real ( kind = 8 ) x(ni,nj)
  real ( kind = 8 ) xc(ni,nj)
  real ( kind = 8 ) xd
  real ( kind = 8 ) xl
  real ( kind = 8 ) xr
  real ( kind = 8 ) xu
  real ( kind = 8 ) y(ni,nj)
  real ( kind = 8 ) yc(ni,nj)
  real ( kind = 8 ) yd
  real ( kind = 8 ) yl
  real ( kind = 8 ) yr
  real ( kind = 8 ) yu

  do n = 1, ns

    if ( lsolve(n) ) then

      nf = n
!
!  N = 1, U VELOCITY.
!
      if ( n == 1 ) then

        call gamsor(ap,birad,ce1,ce2,cmu,con,ewall,f,fcsl,fksl, &
          fma,fo,frsl,gam,gamt,grash,hamag,heta,hksi,icrys,inturb, &
          izone,jcrys,l1,lsolve,m1,mode,nf,pr,r,rdtm,re,recb,rect, &
          rho,rhocon,rpr,sige,sigk,sigt,stel,stes,tal,tas,tf,tw,x, &
          xc,y,yc)

        do j = 2, m1-1

          diff = 2.0D+00*gam(2,j)/hksi(2,j)
          flow = -ruksi(2,j)
          call diflow(acof,diff,flow)
          aim(2,j) = acof*ak1(2,j)

          diff = 2.0D+00*gam(l1-1,j)/hksi(l1-1,j)
          flow = ruksi(l1,j)
          call diflow(acof,diff,flow)
          aip(l1-1,j) = acof*ak1(l1,j)

          do i = 2, l1-2
            diff = gam(i,j)*gam(i+1,j)/(0.5D+00*hksi(i,j)*gam(i+1,j) &
              +0.5D+00*hksi(i+1,j)*gam(i,j))
            flow = ruksi(i+1,j)
            call diflow(acof,diff,flow)
            aip(i,j) = acof*ak1(i+1,j)
            aim(i+1,j) = (acof+ruksi(i+1,j))*ak1(i+1,j)
          end do

        end do

        do i = 2, l1-1

          diff = 2.0D+00*gam(i,2)/heta(i,2)
          flow = -rueta(i,2)
          call diflow(acof,diff,flow)
          ajm(i,2) = acof*ae1(i,2)

          diff = 2.0D+00*gam(i,m1-1)/heta(i,m1-1)
          flow = rueta(i,m1)
          call diflow(acof,diff,flow)
          ajp(i,m1-1) = acof*ae1(i,m1)

          do j = 2, m1-2
            diff = gam(i,j)*gam(i,j+1)/(0.5D+00*heta(i,j)*gam(i,j+1) &
              +0.5D+00*heta(i,j+1)*gam(i,j))
            flow = rueta(i,j+1)
            call diflow(acof,diff,flow)
            ajp(i,j) = acof*ae1(i,j+1)
            ajm(i,j+1) = (acof+rueta(i,j+1))*ae1(i,j+1)
          end do

        end do

        do j = 1, m1
          biu(j) = gam(1,j)
          blu(j) = gam(l1,j)
        end do

        do i = 1, l1
          bju(i) = gam(i,1)
          bmu(i) = gam(i,m1)
        end do

        do i = 1, l1
          do j = 1, m1
            conu(i,j) = con(i,j)
            cofu(i,j,1) = ap(i,j)
            cofu(i,j,2) = aip(i,j)
            cofu(i,j,3) = aim(i,j)
            cofu(i,j,4) = ajp(i,j)
            cofu(i,j,5) = ajm(i,j)
          end do
        end do

        do i = 2, l1-1
          do j = 2, m1-1
            con(i,j) = (con(i,j)+ap(i,j)*f(i,j,1))*vol(i,j)
          end do
        end do

        isol = 0

        call flux(ae1,ae2,aim,aip,ajm,ajp,ak1,ak2,ap,b1jbl,b2jbl, &
          cappa,cd,cmu,con,epsil,f,fjeta,fjksi,fmax,fn,fo,gam,heta,hksi, &
          icrys,isol,iter,izone,jcrys,l1,lblk,lconv,lortho, &
          m1,mode,nf,nsolve,ntimes,r,relax,res,rueta,ruksi)

        do i = 2, l1-1
          do j = 2, m1-1
            hutop(i,j) = (aip(i,j)*f(i+1,j,1)+aim(i,j)*f(i-1,j,1) &
              +ajp(i,j)*f(i,j+1,1)+ajm(i,j)*f(i,j-1,1)+con(i,j))/vol(i,j)
          end do
        end do
!
!  N = 2, V VELOCITY.
!
      else if ( n == 2 ) then

        call gamsor(ap,birad,ce1,ce2,cmu,con,ewall,f,fcsl,fksl, &
          fma,fo,frsl,gam,gamt,grash,hamag,heta,hksi,icrys,inturb, &
          izone,jcrys,l1,lsolve,m1,mode,nf,pr,r,rdtm,re,recb,rect, &
          rho,rhocon,rpr,sige,sigk,sigt,stel,stes,tal,tas,tf,tw,x,xc,y,yc)

        do j = 2, m1-1

          diff = 2.0D+00*gam(2,j)/hksi(2,j)
          flow = -ruksi(2,j)
          call diflow(acof,diff,flow)
          aim(2,j) = acof*ak1(2,j)

          diff = 2.0D+00*gam(l1-1,j)/hksi(l1-1,j)
          flow = ruksi(l1,j)
          call diflow(acof,diff,flow)
          aip(l1-1,j) = acof*ak1(l1,j)

          do i = 2, l1-2
            diff = gam(i,j)*gam(i+1,j)/(0.5D+00*hksi(i,j)*gam(i+1,j) &
              +0.5D+00*hksi(i+1,j)*gam(i,j))
            flow = ruksi(i+1,j)
            call diflow(acof,diff,flow)
            aip(i,j) = acof*ak1(i+1,j)
            aim(i+1,j) = (acof+ruksi(i+1,j))*ak1(i+1,j)
          end do

        end do

        do i = 2, l1-1

          diff = 2.0D+00*gam(i,2)/heta(i,2)
          flow = -rueta(i,2)
          call diflow(acof,diff,flow)
          ajm(i,2) = acof*ae1(i,2)

          diff = 2.0D+00*gam(i,m1-1)/heta(i,m1-1)
          flow = rueta(i,m1)
          call diflow(acof,diff,flow)
          ajp(i,m1-1) = acof*ae1(i,m1)

          do j = 2, m1-2
            diff = gam(i,j)*gam(i,j+1)/(0.5D+00*heta(i,j)*gam(i,j+1) &
              +0.5D+00*heta(i,j+1)*gam(i,j))
            flow = rueta(i,j+1)
            call diflow(acof,diff,flow)
            ajp(i,j) = acof*ae1(i,j+1)
            ajm(i,j+1) = (acof+rueta(i,j+1))*ae1(i,j+1)
          end do

        end do

        biv(1:m1) = gam(1,1:m1)
        blv(1:m1) = gam(l1,1:m1)

        do i = 1, l1
          bjv(i) = gam(i,1)
          bmv(i) = gam(i,m1)
        end do

        do i = 1, l1
          do j = 1, m1
            conv(i,j) = con(i,j)
            cofv(i,j,1) = ap(i,j)
            cofv(i,j,2) = aip(i,j)
            cofv(i,j,3) = aim(i,j)
            cofv(i,j,4) = ajp(i,j)
            cofv(i,j,5) = ajm(i,j)
          end do
        end do

        do i = 2, l1-1
          do j = 2, m1-1
            con(i,j) = (con(i,j)+ap(i,j)*f(i,j,2))*vol(i,j)
          end do
        end do

        isol = 0

        call flux(ae1,ae2,aim,aip,ajm,ajp,ak1,ak2,ap,b1jbl,b2jbl, &
          cappa,cd,cmu,con,epsil,f,fjeta,fjksi,fmax,fn,fo,gam,heta,hksi, &
          icrys,isol,iter,izone,jcrys,l1,lblk,lconv,lortho, &
          m1,mode,nf,nsolve,ntimes,r,relax,res,rueta,ruksi)

        do i = 2, l1-1
          do j = 2, m1-1
            hvtop(i,j) = (aip(i,j)*f(i+1,j,2)+aim(i,j)*f(i-1,j,2) &
              +ajp(i,j)*f(i,j+1,2)+ajm(i,j)*f(i,j-1,2)+con(i,j))/vol(i,j)
            dxdksi = 0.5*(xc(i+1,j+1)+xc(i+1,j)-xc(i,j+1)-xc(i,j))
            dydksi = 0.5*(yc(i+1,j+1)+yc(i+1,j)-yc(i,j+1)-yc(i,j))
            hksip(i,j) = 0.5*(dxdksi*hutop(i,j)+dydksi*hvtop(i,j))
            dxdeta = 0.5*(xc(i+1,j+1)+xc(i,j+1)-xc(i+1,j)-xc(i,j))
            dydeta = 0.5*(yc(i+1,j+1)+yc(i,j+1)-yc(i+1,j)-yc(i,j))
            hetap(i,j) = 0.5*(dxdeta*hutop(i,j)+dydeta*hvtop(i,j))
          end do
        end do
!
!  N = 3, P PRESSURE.
!
!  HKSIP (stored in CON0), HETAP (stored in AP0)
!  are based on momentum interpolation.
!  T3 is based on the two-dimenisonal correction of MIS
!  T4 is based on the secondary correction of 2-d of MIS
!
      else if ( n == 3 ) then

        do i = 2, l1-1
          do j = 2, m1-1
            con(i,j) = 0.0D+00
            ap(i,j) = (aip(i,j)+aim(i,j)+ajm(i,j)+ajp(i,j))/vol(i,j)/rho(i,j)
            tmp(i,j) = ap(i,j)
          end do
        end do

        if ( .not.lortho ) then

          do j = 2, m1-1

            do i = 1, l1
              qt(i) = 0.5 * (rueta(i,j)+rueta(i,j+1))
            end do

            do i = 2, l1-1
              tem1 = hksi(i+1,j)+hksi(i,j)
              tem2 = hksi(i-1,j)+hksi(i,j)
              t1 = hksi(i,j)/tem1
              t2 = hksi(i+1,j)/tem1
              t3 = hksi(i,j)/tem2
              t4 = hksi(i-1,j)/tem2
              con(i,j) = con(i,j)+(t1*qt(i+1)+t2*qt(i))*ak2(i+1,j) &
                -(t4*qt(i)+t3*qt(i-1))*ak2(i,j)
            end do

          end do

          do i = 2, l1-1

            do j = 1, m1
              qt(j) = 0.5 * (ruksi(i,j)+ruksi(i+1,j))
            end do

            do j = 2, m1-1
              tem1 = heta(i,j+1)+heta(i,j)
              tem2 = heta(i,j-1)+heta(i,j)
              t1 = heta(i,j)/tem1
              t2 = heta(i,j+1)/tem1
              t3 = heta(i,j)/tem2
              t4 = heta(i,j-1)/tem2
              con(i,j) = con(i,j)+(t1*qt(j+1)+t2*qt(j))*ae2(i,j+1) &
                -(t4*qt(j)+t3*qt(j-1))*ae2(i,j)
            end do
          end do

        end if

        do j = 2, m1-1

          temp = 0.0D+00

          do i = 2, l1-2

            a1 = 0.5D+00 * (ak1(i,j)+ak1(i+1,j))
            a2 = 0.5D+00 * (ak1(i+1,j)+ak1(i+2,j))
            denom = 0.5D+00 * (ap(i,j)*hksi(i,j)/a1+ap(i+1,j)*hksi(i+1,j)/a2)
            apr = 2.0D+00 * denom/(hksi(i,j)+hksi(i+1,j))

            dxdksi = 0.5D+00 * (xc(i+1,j+1)+xc(i+1,j)-xc(i,j+1)-xc(i,j))
            dydksi = 0.5D+00 * (yc(i+1,j+1)+yc(i+1,j)-yc(i,j+1)-yc(i,j))
            rukij = (dxdksi*f(i,j,1)+dydksi*f(i,j,2))*rho(i,j)/hksi(i,j)*a1

            dxdksi = 0.5D+00 * (xc(i+2,j+1)+xc(i+2,j)-xc(i+1,j+1)-xc(i+1,j))
            dydksi = 0.5D+00 * (yc(i+2,j+1)+yc(i+2,j)-yc(i+1,j+1)-yc(i+1,j))
            ruki1j = (dxdksi*f(i+1,j,1)+dydksi*f(i+1,j,2))*rho(i+1,j) &
              /hksi(i+1,j)*a2

            t1 = rukij*(apr*hksi(i+1,j)-ap(i,j)*hksi(i,j)/a1)
            t2 = ruki1j*(apr*hksi(i,j)-ap(i+1,j)*hksi(i+1,j)/a2)
            t3 = 0.5 * (t1+t2)

            if ( i == 2 ) then
              rukt = ruksi(2,j)
            else
              rukt = temp
            end if

            temp = ruksi(i+1,j)
            frc = hksi(i+1,j)/(hksi(i,j)+hksi(i+1,j))
            t4 = 0.5*(frc*(ruksi(i+1,j)*ak1(i+1,j)-rukt*ak1(i,j)) &
              +(1.0-frc)*(ruksi(i+1,j)*ak1(i+1,j)-ruksi(i+2,j)*ak1(i+2,j)))

            ruksi(i+1,j) = ((hksip(i,j)+hksip(i+1,j)+t3)/denom+t4)/ak1(i+1,j)

            aip(i,j) = 1.0D+00 / denom
            aim(i+1,j) = aip(i,j)

          end do

          aip(l1-1,j) = 0.0D+00
          aim(2,j) = 0.0D+00

        end do

        do i = 2, l1-1
          do j = 2, m1-2

            a1 = 0.5D+00*(ae1(i,j)+ae1(i,j+1))
            a2 = 0.5D+00*(ae1(i,j+1)+ae1(i,j+2))
            denom = 0.5D+00*(ap(i,j)*heta(i,j)/a1+ap(i,j+1)*heta(i,j+1)/a2)
            apr = 2.0D+00*denom/(heta(i,j)+heta(i,j+1))

            dxdeta = 0.5D+00*(xc(i+1,j+1)+xc(i,j+1)-xc(i+1,j)-xc(i,j))
            dydeta = 0.5D+00*(yc(i+1,j+1)+yc(i,j+1)-yc(i+1,j)-yc(i,j))
            rueij = (dxdeta*f(i,j,1)+dydeta*f(i,j,2))*rho(i,j)/heta(i,j)*a1

            dxdeta = 0.5D+00*(xc(i+1,j+2)+xc(i,j+2)-xc(i+1,j+1)-xc(i,j+1))
            dydeta = 0.5D+00*(yc(i+1,j+2)+yc(i,j+2)-yc(i+1,j+1)-yc(i,j+1))
            rueij1 = (dxdeta*f(i,j+1,1)+dydeta*f(i,j+1,2))*rho(i,j+1) &
              /heta(i,j+1)*a2

            t1 = rueij*(apr*heta(i,j+1)-ap(i,j)*heta(i,j)/a1)
            t2 = rueij1*(apr*heta(i,j)-ap(i,j+1)*heta(i,j+1)/a2)
            t3 = 0.5*(t1+t2)
            ruet = rueta(i,2)
            if ( j /= 2) ruet = temp
            temp = rueta(i,j+1)
            frc = heta(i,j+1)/(heta(i,j)+heta(i,j+1))
            t4 = 0.5*(frc*(rueta(i,j+1)*ae1(i,j+1)-ruet*ae1(i,j)) &
              +(1.0-frc)*(rueta(i,j+1)*ae1(i,j+1)-rueta(i,j+2)*ae1(i,j+2)))

            rueta(i,j+1) = ((hetap(i,j)+hetap(i,j+1)+t3)/denom+t4)/ae1(i,j+1)

            ajp(i,j) = 1.0D+00/denom
            ajm(i,j+1) = ajp(i,j)

          end do

          ajp(i,m1-1) = 0.0D+00
          ajm(i,2) = 0.0D+00

        end do

        do i = 2, l1-1
          do j = 2, m1-1
            con(i,j) = con(i,j)+ak1(i,j)*ruksi(i,j)-ak1(i+1,j)*ruksi(i+1,j) &
              +ae1(i,j)*rueta(i,j)-ae1(i,j+1)*rueta(i,j+1)
            ap(i,j) = aip(i,j)+aim(i,j)+ajp(i,j)+ajm(i,j)
          end do
        end do

        smax = 0.0D+00
        ssum = 0.0D+00
        do i = 2, l1-1
          do j = 2, m1-1
            soor = ap(i,j)*f(i,j,3)-aip(i,j)*f(i+1,j,3) &
              -aim(i,j)*f(i-1,j,3)-ajp(i,j)*f(i,j+1,3) &
              -ajm(i,j)*f(i,j-1,3)-con(i,j)
            ssum = ssum+soor/rho(i,j)
            smax = max(abs(soor)/1.0d4,smax)
          end do
        end do

        do i = 2, l1-1
          do j = 2, m1-1
            ap(i,j) = ap(i,j)/relax(np)
            con(i,j) = con(i,j)+(1.0-relax(np))*ap(i,j)*f(i,j,3)
          end do
        end do

        call solve1(aim,aip,ajm,ajp,ap,cappa,cd,cmu,con,f,fo, &
          heta,hksi,icrys,izone,jcrys,l1,m1,nf)

        call solve2(aim,aip,ajm,ajp,ap,con,f,l1,lblk,m1,nf,nsolve)

        do j = 2, m1-1
          do i = 2, l1-2
            ruksi(i+1,j) = ruksi(i+1,j)+aip(i,j)/ak1(i+1,j) &
              *(f(i,j,nf)-f(i+1,j,nf))
          end do
        end do

        do i = 2, l1-1
          do j = 2, m1-2
            rueta(i,j+1) = rueta(i,j+1)+ajp(i,j)/ae1(i,j+1) &
              *(f(i,j,nf)-f(i,j+1,nf))
          end do
        end do

        do j = 2, m1-1

          do i = 2, l1-2
 
            bpi = f(i,j,3)
            bpi1 = f(i+1,j,3)-hksip(i,j)-hksip(i+1,j)
            em = hksi(i,j)*tmp(i,j)/(ak1(i,j)+ak1(i+1,j))
            ep = hksi(i+1,j)*tmp(i+1,j)/(ak1(i+1,j)+ak1(i+2,j))

            t1 = 0.5*em*ep*(ruksi(i+2,j)*ak1(i+2,j) &
              -ruksi(i,j)*ak1(i,j))/(em+ep)
            bpm = (bpi*ep+bpi1*em)/(em+ep)+t1
            pksi(i+1,j) = bpm+hksip(i,j)

          end do

          pksi(l1,j) = 2.0*f(l1-1,j,nf)-pksi(l1-1,j)
          pksi(2,j) = 2.0*f(2,j,nf)-pksi(3,j)
          f(1,j,nf) = pksi(2,j)
          f(l1,j,nf) = pksi(l1,j)

        end do

        do i = 2, l1-1
          do j = 2, m1-2
            bpj = f(i,j,3)
            bpj1 = f(i,j+1,3)-hetap(i,j)-hetap(i,j+1)
            em = heta(i,j)*tmp(i,j)/(ae1(i,j)+ae1(i,j+1))
            ep = heta(i,j+1)*tmp(i,j+1)/(ae1(i,j+1)+ae1(i,j+2))

            t1 = 0.5*em*ep*(rueta(i,j+2)*ae1(i,j+2) &
              -rueta(i,j)*ae1(i,j))/(em+ep)
            bpm = (bpj*ep+bpj1*em)/(em+ep)+t1
            peta(i,j+1) = bpm+hetap(i,j)

          end do

          peta(i,m1) = 2.0*f(i,m1-1,nf)-peta(i,m1-1)
          peta(i,2) = 2.0*f(i,2,nf)-peta(i,3)
          f(i,1,nf) = peta(i,2)
          f(i,m1,nf) = peta(i,m1)

        end do

        do j = 1, m1
          gam(1,j) = biu(j)
          gam(l1,j) = blu(j)
        end do

        do i = 1, l1
          gam(i,1) = bju(i)
          gam(i,m1) = bmu(i)
        end do

        do i = 1, l1
          do j = 1, m1
            con(i,j) = conu(i,j)
            ap(i,j) = cofu(i,j,1)
            aip(i,j) = cofu(i,j,2)
            aim(i,j) = cofu(i,j,3)
            ajp(i,j) = cofu(i,j,4)
            ajm(i,j) = cofu(i,j,5)
          end do
        end do

        do i = 2, l1-1

          if ( gam(i,1) == 0.0D+00 ) then
            ajm(i,2) = 0.0D+00
          end if

          if ( gam(i,m1) == 0.0D+00 ) then
            ajp(i,m1-1) = 0.0D+00
          end if

        end do

        do j = 2, m1-1

          if ( gam(1,j) == 0.0D+00 ) then
            aim(2,j) = 0.0D+00
          end if

          if ( gam(l1,j) == 0.0D+00 ) then
            aip(l1-1,j) = 0.0D+00
          end if

        end do

        do i = 2, l1-1
          do j = 2, m1-1

            xr = 0.5*(xc(i+1,j+1)+xc(i+1,j))
            yr = 0.5*(yc(i+1,j+1)+yc(i+1,j))
            xl = 0.5*(xc(i,j+1)+xc(i,j))
            yl = 0.5*(yc(i,j+1)+yc(i,j))
            xu = 0.5*(xc(i+1,j+1)+xc(i,j+1))
            yu = 0.5*(yc(i+1,j+1)+yc(i,j+1))
            xd = 0.5*(xc(i+1,j)+xc(i,j))
            yd = 0.5*(yc(i+1,j)+yc(i,j))
            dxdksi = xr-xl
            dydksi = yr-yl
            dxdeta = xu-xd
            dydeta = yu-yd
            dpksi = pksi(i+1,j)-pksi(i,j)
            dpeta = peta(i,j+1)-peta(i,j)
            t1 = -dydeta*dpksi+dydksi*dpeta
            t2 = dxdeta*dpksi-dxdksi*dpeta

            if ( mode == 1 ) then
              t1 = t1*r(i,j)
              t2 = t2*r(i,j)
            end if

            tmp(i,j) = t2
            con(i,j) = con(i,j)*vol(i,j)+t1
            ap(i,j) = -ap(i,j)*vol(i,j)+aip(i,j)+aim(i,j)+ajp(i,j)+ajm(i,j)

          end do
        end do

        nf = 1
        isol = 1

        call flux(ae1,ae2,aim,aip,ajm,ajp,ak1,ak2,ap,b1jbl,b2jbl,cappa, &
          cd,cmu,con,epsil,f,fjeta,fjksi,fmax,fn,fo,gam,heta,hksi, &
          icrys,isol,iter,izone,jcrys,l1,lblk,lconv,lortho, &
          m1,mode,nf,nsolve,ntimes,r,relax,res,rueta,ruksi)

        do j = 1, m1
          gam(1,j) = biv(j)
          gam(l1,j) = blv(j)
        end do

        do i = 1, l1
          gam(i,1) = bjv(i)
          gam(i,m1) = bmv(i)
        end do

        do i = 1, l1
          do j = 1, m1
            con(i,j) = conv(i,j)
            ap(i,j) = cofv(i,j,1)
            aip(i,j) = cofv(i,j,2)
            aim(i,j) = cofv(i,j,3)
            ajp(i,j) = cofv(i,j,4)
            ajm(i,j) = cofv(i,j,5)
          end do
        end do

        nf = 2

        do i = 2, l1-1

          if ( gam(i,1) == 0.0D+00 ) then
            ajm(i,2) = 0.0D+00
          end if

          if ( gam(i,m1) == 0.0D+00 ) then
            ajp(i,m1-1) = 0.0D+00
          end if

        end do

        do j = 2, m1-1

          if ( gam(1,j) == 0.0D+00 ) then
            aim(2,j) = 0.0D+00
          end if

          if ( gam(l1,j) == 0.0D+00 ) then
            aip(l1-1,j) = 0.0D+00
          end if

        end do

        do i = 2, l1-1
          do j = 2, m1-1
            con(i,j) = con(i,j)*vol(i,j)+tmp(i,j)
            ap(i,j) = -ap(i,j)*vol(i,j)+aip(i,j)+aim(i,j)+ajp(i,j)+ajm(i,j)
          end do
        end do

        isol = 1

        call flux(ae1,ae2,aim,aip,ajm,ajp,ak1,ak2,ap,b1jbl,b2jbl, &
          cappa,cd,cmu,con,epsil,f,fjeta,fjksi,fmax,fn,fo,gam,heta,hksi, &
          icrys,isol,iter,izone,jcrys,l1,lblk,lconv,lortho, &
          m1,mode,nf,nsolve,ntimes,r,relax,res,rueta,ruksi)
!
!  N = 4, PC CORRECTED PRESSURE.
!
      else if ( n == 4 ) then

        do j = 1, m1
          gam(1,j) = biu(j)
          gam(l1,j) = blu(j)
        end do

        do i = 1, l1
          gam(i,1) = bju(i)
          gam(i,m1) = bmu(i)
        end do

        do i = 1, l1
          do j = 1, m1
            con(i,j) = conu(i,j)
            ap(i,j) = cofu(i,j,1)
            aip(i,j) = cofu(i,j,2)
            aim(i,j) = cofu(i,j,3)
            ajp(i,j) = cofu(i,j,4)
            ajm(i,j) = cofu(i,j,5)
          end do
        end do

        do i = 2, l1-1
          do j = 2, m1-1
            con(i,j) = (con(i,j)+ap(i,j)*f(i,j,1))*vol(i,j)
          end do
        end do

        nf = 1

        isol = 0

        call flux(ae1,ae2,aim,aip,ajm,ajp,ak1,ak2,ap,b1jbl,b2jbl, &
          cappa,cd,cmu,con,epsil,f,fjeta,fjksi,fmax,fn,fo,gam,heta,hksi, &
          icrys,isol,iter,izone,jcrys,l1,lblk,lconv,lortho, &
          m1,mode,nf,nsolve,ntimes,r,relax,res,rueta,ruksi)

        nf = npc

        do i = 2, l1-1
          do j = 2, m1-1
            hutop(i,j) = (aip(i,j)*f(i+1,j,1)+aim(i,j)*f(i-1,j,1) &
              +ajp(i,j)*f(i,j+1,1)+ajm(i,j)*f(i,j-1,1)+con(i,j))/vol(i,j)
          end do
        end do

        do j = 2, m1-1
          do i = 2, l1-1
            con(i,j) = 0.0D+00
            ap(i,j) = 0.0D+00
          end do
        end do

        do j = 1, m1
          gam(1,j) = biv(j)
          gam(l1,j) = blv(j)
        end do

        do i = 1, l1
          gam(i,1) = bjv(i)
          gam(i,m1) = bmv(i)
        end do

        do i = 1, l1
          do j = 1, m1
            con(i,j) = conv(i,j)
            ap(i,j) = cofv(i,j,1)
            aip(i,j) = cofv(i,j,2)
            aim(i,j) = cofv(i,j,3)
            ajp(i,j) = cofv(i,j,4)
            ajm(i,j) = cofv(i,j,5)
          end do
        end do

        do i = 2, l1-1
          do j = 2, m1-1
            con(i,j) = (con(i,j)+ap(i,j)*f(i,j,2))*vol(i,j)
          end do
        end do

        nf = 2
        isol = 0

        call flux(ae1,ae2,aim,aip,ajm,ajp,ak1,ak2,ap,b1jbl,b2jbl, &
          cappa,cd,cmu,con,epsil,f,fjeta,fjksi,fmax,fn,fo,gam,heta,hksi, &
          icrys,isol,iter,izone,jcrys,l1,lblk,lconv,lortho, &
          m1,mode,nf,nsolve,ntimes,r,relax,res,rueta,ruksi)

        nf = npc

        do i = 2, l1-1
          do j = 2, m1-1
            hvtop(i,j) = (aip(i,j)*f(i+1,j,2)+aim(i,j)*f(i-1,j,2) &
              +ajp(i,j)*f(i,j+1,2)+ajm(i,j)*f(i,j-1,2)+con(i,j))/vol(i,j)
            dxdksi = 0.5*(xc(i+1,j+1)+xc(i+1,j)-xc(i,j+1)-xc(i,j))
            dydksi = 0.5*(yc(i+1,j+1)+yc(i+1,j)-yc(i,j+1)-yc(i,j))
            hksip(i,j) = 0.5*(dxdksi*hutop(i,j)+dydksi*hvtop(i,j))
            dxdeta = 0.5*(xc(i+1,j+1)+xc(i,j+1)-xc(i+1,j)-xc(i,j))
            dydeta = 0.5*(yc(i+1,j+1)+yc(i,j+1)-yc(i+1,j)-yc(i,j))
            hetap(i,j) = 0.5*(dxdeta*hutop(i,j)+dydeta*hvtop(i,j))
          end do
        end do

        do i = 2, l1-1
          do j = 2, m1-1
            con(i,j) = 0.0D+00
            ap(i,j) = (aip(i,j)+aim(i,j)+ajm(i,j)+ajp(i,j))/vol(i,j)/rho(i,j)
            tmp(i,j) = ap(i,j)
          end do
        end do

        if ( .not. lortho ) then

          do j = 2, m1-1

            do i = 1, l1
              qt(i) = 0.5*(rueta(i,j)+rueta(i,j+1))
            end do

            do i = 2, l1-1
              tem1 = hksi(i+1,j)+hksi(i,j)
              tem2 = hksi(i-1,j)+hksi(i,j)
              t1 = hksi(i,j)/tem1
              t2 = hksi(i+1,j)/tem1
              t3 = hksi(i,j)/tem2
              t4 = hksi(i-1,j)/tem2
              con(i,j) = con(i,j)+(t1*qt(i+1)+t2*qt(i))*ak2(i+1,j) &
                -(t4*qt(i)+t3*qt(i-1))*ak2(i,j)
            end do

          end do

          do i = 2, l1-1

            do j = 1, m1
              qt(j) = 0.5*(ruksi(i,j)+ruksi(i+1,j))
            end do

            do j = 2, m1-1
              tem1 = heta(i,j+1)+heta(i,j)
              tem2 = heta(i,j-1)+heta(i,j)
              t1 = heta(i,j)/tem1
              t2 = heta(i,j+1)/tem1
              t3 = heta(i,j)/tem2
              t4 = heta(i,j-1)/tem2
              con(i,j) = con(i,j)+(t1*qt(j+1)+t2*qt(j))*ae2(i,j+1) &
                -(t4*qt(j)+t3*qt(j-1))*ae2(i,j)
            end do
          end do

        end if

        do j = 2, m1-1

          temp = 0.0D+00

          do i = 2, l1-2

            a1 = 0.5*(ak1(i,j)+ak1(i+1,j))
            a2 = 0.5*(ak1(i+1,j)+ak1(i+2,j))
            denom = 0.5*(ap(i,j)*hksi(i,j)/a1+ap(i+1,j)*hksi(i+1,j)/a2)
            apr = 2.0*denom/(hksi(i,j)+hksi(i+1,j))

            dxdksi = 0.5*(xc(i+1,j+1)+xc(i+1,j)-xc(i,j+1)-xc(i,j))
            dydksi = 0.5*(yc(i+1,j+1)+yc(i+1,j)-yc(i,j+1)-yc(i,j))
            rukij = (dxdksi*f(i,j,1)+dydksi*f(i,j,2))*rho(i,j)/hksi(i,j)*a1

            dxdksi = 0.5*(xc(i+2,j+1)+xc(i+2,j)-xc(i+1,j+1)-xc(i+1,j))
            dydksi = 0.5*(yc(i+2,j+1)+yc(i+2,j)-yc(i+1,j+1)-yc(i+1,j))
            ruki1j = (dxdksi*f(i+1,j,1)+dydksi*f(i+1,j,2))*rho(i+1,j) &
              /hksi(i+1,j)*a2

            t1 = rukij*(apr*hksi(i+1,j)-ap(i,j)*hksi(i,j)/a1)
            t2 = ruki1j*(apr*hksi(i,j)-ap(i+1,j)*hksi(i+1,j)/a2)
            t3 = 0.5*(t1+t2)

            if ( i == 2 ) then
              rukt = ruksi(2,j)
            else
              rukt = temp
            end if

            temp = ruksi(i+1,j)
            frc = hksi(i+1,j)/(hksi(i,j)+hksi(i+1,j))
            t4 = 0.5*(frc*(ruksi(i+1,j)*ak1(i+1,j)-rukt*ak1(i,j)) &
              +(1.0-frc)*(ruksi(i+1,j)*ak1(i+1,j) &
              -ruksi(i+2,j)*ak1(i+2,j)))

            ruksi(i+1,j) = ((hksip(i,j)+hksip(i+1,j)+t3)/denom+t4)/ak1(i+1,j)

            ruksi(i+1,j) = ruksi(i+1,j)+(f(i,j,3)-f(i+1,j,3))/denom/ak1(i+1,j)

            aip(i,j) = 1.0D+00/denom
            aim(i+1,j) = aip(i,j)

          end do

          aip(l1-1,j) = 0.0D+00
          aim(2,j) = 0.0D+00

        end do

        do i = 2, l1-1
          do j = 2, m1-2

            a1 = 0.5*(ae1(i,j)+ae1(i,j+1))
            a2 = 0.5*(ae1(i,j+1)+ae1(i,j+2))
            denom = 0.5*(ap(i,j)*heta(i,j)/a1+ap(i,j+1)*heta(i,j+1)/a2)
            apr = 2.0*denom/(heta(i,j)+heta(i,j+1))

            dxdeta = 0.5*(xc(i+1,j+1)+xc(i,j+1)-xc(i+1,j)-xc(i,j))
            dydeta = 0.5*(yc(i+1,j+1)+yc(i,j+1)-yc(i+1,j)-yc(i,j))
            rueij = (dxdeta*f(i,j,1)+dydeta*f(i,j,2))*rho(i,j)/heta(i,j)*a1

            dxdeta = 0.5*(xc(i+1,j+2)+xc(i,j+2)-xc(i+1,j+1)-xc(i,j+1))
            dydeta = 0.5*(yc(i+1,j+2)+yc(i,j+2)-yc(i+1,j+1)-yc(i,j+1))
            rueij1 = (dxdeta*f(i,j+1,1)+dydeta*f(i,j+1,2))*rho(i,j+1) &
              /heta(i,j+1)*a2

            t1 = rueij*(apr*heta(i,j+1)-ap(i,j)*heta(i,j)/a1)
            t2 = rueij1*(apr*heta(i,j)-ap(i,j+1)*heta(i,j+1)/a2)
            t3 = 0.5*(t1+t2)
            ruet = rueta(i,2)
            if ( j /= 2) ruet = temp
            temp = rueta(i,j+1)
            frc = heta(i,j+1)/(heta(i,j)+heta(i,j+1))
            t4 = 0.5*(frc*(rueta(i,j+1)*ae1(i,j+1)-ruet*ae1(i,j)) &
              +(1.0-frc)*(rueta(i,j+1)*ae1(i,j+1)-rueta(i,j+2)*ae1(i,j+2)))

            rueta(i,j+1) = ((hetap(i,j)+hetap(i,j+1)+t3)/denom+t4)/ae1(i,j+1)

            rueta(i,j+1) = rueta(i,j+1) &
              +(f(i,j,3)-f(i,j+1,3))/denom/ae1(i,j+1)

            ajp(i,j) = 1.0D+00/denom
            ajm(i,j+1) = ajp(i,j)

          end do

          ajp(i,m1-1) = 0.0D+00
          ajm(i,2) = 0.0D+00

        end do

        do i = 2, l1-1
          do j = 2, m1-1
            con(i,j) = con(i,j)+ak1(i,j)*ruksi(i,j) &
              -ak1(i+1,j)*ruksi(i+1,j) &
              +ae1(i,j)*rueta(i,j)-ae1(i,j+1)*rueta(i,j+1)
            ap(i,j) = aip(i,j)+aim(i,j)+ajp(i,j)+ajm(i,j)
          end do
        end do

        call solve1(aim,aip,ajm,ajp,ap,cappa,cd,cmu,con,f,fo, &
          heta,hksi,icrys,izone,jcrys,l1,m1,nf)

        call solve2(aim,aip,ajm,ajp,ap,con,f,l1,lblk,m1,nf,nsolve)

        do j = 2, m1-1
          do i = 2, l1-2
            ruksi(i+1,j) = ruksi(i+1,j)+aip(i,j)/ak1(i+1,j) &
              *(f(i,j,nf)-f(i+1,j,nf))
          end do
        end do

        do i = 2, l1-1
          do j = 2, m1-2
            rueta(i,j+1) = rueta(i,j+1)+ajp(i,j)/ae1(i,j+1) &
              *(f(i,j,nf)-f(i,j+1,nf))
          end do
        end do

        do j = 2, m1-1

          do i = 2, l1-2

            bpi = f(i,j,3)
            bpi1 = f(i+1,j,3)-hksip(i,j)-hksip(i+1,j)
            em = hksi(i,j)*tmp(i,j)/(ak1(i,j)+ak1(i+1,j))
            ep = hksi(i+1,j)*tmp(i+1,j)/(ak1(i+1,j)+ak1(i+2,j))
            pksi(i+1,j) = (f(i,j,nf)*ep+f(i+1,j,nf)*em)/(em+ep)

          end do

          pksi(l1,j) = 2.0*f(l1-1,j,nf)-pksi(l1-1,j)
          pksi(2,j) = 2.0*f(2,j,nf)-pksi(3,j)
          f(1,j,nf) = pksi(2,j)
          f(l1,j,nf) = pksi(l1,j)

        end do

        do i = 2, l1-1
          do j = 2, m1-2

            bpj = f(i,j,3)
            bpj1 = f(i,j+1,3)-hetap(i,j)-hetap(i,j+1)
            em = heta(i,j)*tmp(i,j)/(ae1(i,j)+ae1(i,j+1))
            ep = heta(i,j+1)*tmp(i,j+1)/(ae1(i,j+1)+ae1(i,j+2))
            peta(i,j+1) = (f(i,j,nf)*ep+f(i,j+1,nf)*em)/(em+ep)

          end do

          peta(i,m1) = 2.0*f(i,m1-1,nf)-peta(i,m1-1)
          peta(i,2) = 2.0*f(i,2,nf)-peta(i,3)
          f(i,1,nf) = peta(i,2)
          f(i,m1,nf) = peta(i,m1)

        end do

        do i = 2, l1-1
          do j = 2, m1-1

            xr = 0.5*(xc(i+1,j+1)+xc(i+1,j))
            yr = 0.5*(yc(i+1,j+1)+yc(i+1,j))
            xl = 0.5*(xc(i,j+1)+xc(i,j))
            yl = 0.5*(yc(i,j+1)+yc(i,j))
            xu = 0.5*(xc(i+1,j+1)+xc(i,j+1))
            yu = 0.5*(yc(i+1,j+1)+yc(i,j+1))
            xd = 0.5*(xc(i+1,j)+xc(i,j))
            yd = 0.5*(yc(i+1,j)+yc(i,j))
            dxdksi = xr-xl
            dydksi = yr-yl
            dxdeta = xu-xd
            dydeta = yu-yd
            dpksi = pksi(i+1,j)-pksi(i,j)
            dpeta = peta(i,j+1)-peta(i,j)
            t1 = -dydeta*dpksi+dydksi*dpeta
            t2 = dxdeta*dpksi-dxdksi*dpeta

            if ( mode == 1 ) then
              t1 = t1*r(i,j)
              t2 = t2*r(i,j)
            end if

            f(i,j,1) = f(i,j,1)+t1/tmp(i,j)/rho(i,j)/vol(i,j)
            f(i,j,2) = f(i,j,2)+t2/tmp(i,j)/rho(i,j)/vol(i,j)

          end do
        end do
!
!  N > 4, T, W, TK, TE, E, PSI
!
      else

        call gamsor(ap,birad,ce1,ce2,cmu,con,ewall,f,fcsl,fksl, &
          fma,fo,frsl,gam,gamt,grash,hamag,heta,hksi,icrys,inturb, &
          izone,jcrys,l1,lsolve,m1,mode,nf,pr,r,rdtm,re,recb,rect, &
          rho,rhocon,rpr,sige,sigk,sigt,stel,stes,tal,tas,tf,tw,x, &
          xc,y,yc)

        do j = 2, m1-1

          diff = 2.0*gam(1,j)/hksi(2,j)
          flow = -ruksi(2,j)
          call diflow(acof,diff,flow)
          aim(2,j) = acof*ak1(2,j)

          diff = 2.0*gam(l1,j)/hksi(l1-1,j)
          flow = ruksi(l1,j)
          call diflow(acof,diff,flow)
          aip(l1-1,j) = acof*ak1(l1,j)

          do i = 2, l1-2
            diff = gam(i,j)*gam(i+1,j)/(0.5*hksi(i,j)*gam(i+1,j) &
              +0.5*hksi(i+1,j)*gam(i,j))
            flow = ruksi(i+1,j)
            call diflow(acof,diff,flow)
            aip(i,j) = acof*ak1(i+1,j)
            aim(i+1,j) = (acof+ruksi(i+1,j))*ak1(i+1,j)
          end do

        end do

        do i = 2, l1-1

          diff = 2.0*gam(i,1)/heta(i,2)
          flow = -rueta(i,2)
          call diflow(acof,diff,flow)
          ajm(i,2) = acof*ae1(i,2)

          diff = 2.0*gam(i,m1)/heta(i,m1-1)
          flow = rueta(i,m1)
          call diflow(acof,diff,flow)
          ajp(i,m1-1) = acof*ae1(i,m1)

          do j = 2, m1-2
            diff = gam(i,j)*gam(i,j+1)/(0.5*heta(i,j)*gam(i,j+1) &
              +0.5*heta(i,j+1)*gam(i,j))
            flow = rueta(i,j+1)
            call diflow(acof,diff,flow)
            ajp(i,j) = acof*ae1(i,j+1)
            ajm(i,j+1) = (acof+rueta(i,j+1))*ae1(i,j+1)
          end do

        end do

        do i = 2, l1-1
          do j = 2, m1-1
            con(i,j) = con(i,j)*vol(i,j)
            ap(i,j) = -ap(i,j)*vol(i,j)+aip(i,j)+aim(i,j)+ajp(i,j)+ajm(i,j)
          end do
        end do

        isol = 1

        call flux(ae1,ae2,aim,aip,ajm,ajp,ak1,ak2,ap,b1jbl,b2jbl, &
          cappa,cd,cmu,con,epsil,f,fjeta,fjksi,fmax,fn,fo,gam,heta,hksi, &
          icrys,isol,iter,izone,jcrys,l1,lblk,lconv,lortho, &
          m1,mode,nf,nsolve,ntimes,r,relax,res,rueta,ruksi)

      end if

    end if

  end do

  return
end
subroutine setx ( l, m, ni, nj, x, xc, y, yc )

!*****************************************************************************80
!
!! SETX locates the primary nodes from the corner nodes.
!
!  Discussion:
!
!    SETX is given (XC,YC), the locations of the "corners" of the
!    control volumes, and calculates the locations of the primary nodes
!    (X,Y).
!
!    (XC,YC) sits in the center of the control volume.  On the other
!    hand, the point (X,Y) is NOT necessarily the center of the four
!    nearby values of XC and YC.
!
!    SETX may be called to do a portion of the region, or the
!    full region (in which case L = L0, M=M0).
!
  implicit none

  integer ni
  integer nj

  integer i
  integer j
  integer l
  integer m
  real ( kind = 8 ) x(ni,nj)
  real ( kind = 8 ) xc(ni,nj)
  real ( kind = 8 ) y(ni,nj)
  real ( kind = 8 ) yc(ni,nj)
!
!  The central set of values.
!
  do i = 2, l-1
    do j = 2, m-1
      x(i,j) = 0.25*(xc(i,j)+xc(i+1,j)+xc(i,j+1)+xc(i+1,j+1))
      y(i,j) = 0.25*(yc(i,j)+yc(i+1,j)+yc(i,j+1)+yc(i+1,j+1))
    end do
  end do
!
!  The first and last columns, I = 2 to L-1, J=1 and J=M.
!
  do i = 2, l-1
    x(i,1) = 0.5*(xc(i,2)+xc(i+1,2))
    y(i,1) = 0.5*(yc(i,2)+yc(i+1,2))
    x(i,m) = 0.5*(xc(i,m)+xc(i+1,m))
    y(i,m) = 0.5*(yc(i,m)+yc(i+1,m))
  end do
!
!  The first and last rows, I = 1 and I=L, J=2 to M-1.
!
  do j = 2, m-1
    x(1,j) = 0.5*(xc(2,j)+xc(2,j+1))
    y(1,j) = 0.5*(yc(2,j)+yc(2,j+1))
    x(l,j) = 0.5*(xc(l,j)+xc(l,j+1))
    y(l,j) = 0.5*(yc(l,j)+yc(l,j+1))
  end do
!
!  The corners.
!
  x(1,1) = xc(2,2)
  y(1,1) = yc(2,2)

  x(1,m) = xc(2,m)
  y(1,m) = yc(2,m)

  x(l,1) = xc(l,2)
  y(l,1) = yc(l,2)

  x(l,m) = xc(l,m)
  y(l,m) = yc(l,m)

  return
end
subroutine solve1 ( aim, aip, ajm, ajp, ap, cappa, cd, cmu, con, f, fo, &
  heta, hksi, icrys, izone, jcrys, l1, m1, nf )

!*****************************************************************************80
!
!! SOLVE1 sets or modifies the linear system to be handled by SOLVE2.
!
  implicit none

  integer, parameter :: ni = 64
  integer, parameter :: nj = 64
  integer, parameter :: ns = 10

  real ( kind = 8 ) aim(ni,nj)
  real ( kind = 8 ) aip(ni,nj)
  real ( kind = 8 ) ajm(ni,nj)
  real ( kind = 8 ) ajp(ni,nj)
  real ( kind = 8 ) ap(ni,nj)
  real ( kind = 8 ) cappa
  real ( kind = 8 ) cd
  real ( kind = 8 ) cm4
  real ( kind = 8 ) cmu
  real ( kind = 8 ) con(ni,nj)
  real ( kind = 8 ) f(ni,nj,ns)
  real ( kind = 8 ) fo(ni,nj,ns)
  real ( kind = 8 ) heta(ni,nj)
  real ( kind = 8 ) hksi(ni,nj)
  integer i
  integer icrys
  integer izone
  integer j
  integer jcrys
  integer l1
  integer m1
  integer nf
!
!  Temperature calculation in solid zone.
!
  if ( izone == 1 .and. nf == 5 ) then
    do i = 1, icrys
      do j = 1, m1
        ap(i,j) = 1.0D+00
        aip(i,j) = 0.0D+00
        aim(i,j) = 0.0D+00
        ajp(i,j) = 0.0D+00
        ajm(i,j) = 0.0D+00
        con(i,j) = fo(i,j,5)
      end do
    end do
  end if
!
!  Magnetic stream function calculation in solid zone.
!
  if ( izone == 1 .and. nf == 9 ) then
    do i = 1, icrys-1
      do j = 1, m1
        ap(i,j) = 1.0D+00
        aip(i,j) = 0.0D+00
        aim(i,j) = 0.0D+00
        ajp(i,j) = 0.0D+00
        ajm(i,j) = 0.0D+00
        con(i,j) = fo(i,j,9)
      end do
    end do
  end if
!
!  U, V, P, PC, or T calculation in liquid zone.
!
  if ( izone == 2 .and. nf <= 5 ) then
    do j = 2, jcrys
      ap(icrys,j) = 1.0D+00
      aip(icrys,j) = 0.0D+00
      aim(icrys,j) = 0.0D+00
      ajp(icrys,j) = 0.0D+00
      ajm(icrys,j) = 0.0D+00
      con(icrys,j) = 0.0D+00
    end do
  end if
!
!  Turbulent dissipation (TE) calculation.
!
  if ( nf == 8 ) then

    cm4 = cmu**0.25

    do i = 2, l1-1
      con(i,2) = cd*f(i,2,7)**1.5/(cappa*cm4*0.5*heta(i,2))
      ap(i,2) = 1.0
      con(i,m1-1) = cd*f(i,m1-1,7)**1.5/(cappa*cm4*0.5*heta(i,m1-1))
      ap(i,m1-1) = 1.0
    end do

    do j = 2, m1-1
      con(2,j) = cd*f(2,j,7)**1.5/(cappa*cm4*0.5*hksi(2,j))
      ap(2,j) = 1.0
      con(l1-1,j) = cd*f(l1-1,j,7)**1.5/(cappa*cm4*0.5*hksi(l1-1,j))
      ap(l1-1,j) = 1.0
    end do

  end if

  return
end
subroutine solve2 ( aim, aip, ajm, ajp, ap, con, f, l1, lblk, m1, &
  nf, nsolve )

!*****************************************************************************80
!
!! SOLVE2 is the tridiagonal matrix solver.
!
!  Discussion:
!
!    SOLVE2 must solve a set of linear equations defined on a five point
!    stencil, of the form:
!
!    A(I,J)*  U(I,J)
!      - A(I-1,J)*U(I-1,J) - A(I+1,J)*U(I+1,J)
!      - A(I,J-1)*U(I,J-1) - A(I,J+1)*U(I,J+1)  =  CON.
!
!    The equations must be solved for indices 2 < =  I < L1-1
!    and 2 <= J <= M1-1.
!
!    If I = 1 or I=L1 or J=1 or J=M1, then U(I,J) already contains its
!    correct value, and these correct values must be included in the
!    above equations.
!
!    The solution of the equations is stored in F(*,*,NF).
!
  implicit none

  integer, parameter :: ni = 64
  integer, parameter :: nj = 64
  integer, parameter :: nmaxij = 64
  integer, parameter :: ns = 10

  real ( kind = 8 ) aim(ni,nj)
  real ( kind = 8 ) aip(ni,nj)
  real ( kind = 8 ) ajm(ni,nj)
  real ( kind = 8 ) ajp(ni,nj)
  real ( kind = 8 ) ap(ni,nj)
  real ( kind = 8 ) bl
  real ( kind = 8 ) blc
  real ( kind = 8 ) blm
  real ( kind = 8 ) blp
  real ( kind = 8 ) con(ni,nj)
  real ( kind = 8 ) denom
  real ( kind = 8 ) f(ni,nj,ns)
  integer i
  integer j
  integer l1
  logical lblk(ns)
  integer m1
  integer nf
  integer nsolve(ns)
  integer nt
  real ( kind = 8 ) pt(nmaxij)
  real ( kind = 8 ) qt(nmaxij)

  do nt = 1, nsolve(nf)

    if ( lblk(nf) ) then

      pt(1) = 0.0D+00
      qt(1) = 0.0D+00

      do i = 2, l1-1

        bl = 0.0D+00
        blp = 0.0D+00
        blm = 0.0D+00
        blc = 0.0D+00

        do j = 2, m1-1
          bl = bl+ap(i,j)
          if ( j /= m1-1) bl = bl-ajp(i,j)
          if ( j /= 2) bl = bl-ajm(i,j)
          blp = blp+aip(i,j)
          blm = blm+aim(i,j)
          blc = blc+con(i,j)+aip(i,j)*f(i+1,j,nf)+aim(i,j)*f(i-1,j,nf) &
            +ajp(i,j)*f(i,j+1,nf)+ajm(i,j)*f(i,j-1,nf)-ap(i,j)*f(i,j,nf)
        end do

        denom = bl-pt(i-1)*blm

        if ( abs(denom/bl) < 1.0E-10 ) then
          pt(i) = 0.0D+00
          qt(i) = 0.0D+00
        else
          pt(i) = blp/denom
          qt(i) = (blc+blm*qt(i-1))/denom
        end if

      end do

      bl = 0.0D+00
      do i = l1-1, 2, -1

        bl = bl*pt(i)+qt(i)

        do j = 2, m1-1
          f(i,j,nf) = f(i,j,nf)+bl
        end do

      end do

      pt(1) = 0.0D+00
      qt(1) = 0.0D+00

      do j = 2, m1-1

        bl = 0.0D+00
        blp = 0.0D+00
        blm = 0.0D+00
        blc = 0.0D+00

        do i = 2, l1-1
          bl = bl+ap(i,j)
          if ( i /= l1-1) bl = bl-aip(i,j)
          if ( i /= 2) bl = bl-aim(i,j)
          blp = blp+ajp(i,j)
          blm = blm+ajm(i,j)
          blc = blc+con(i,j)+aip(i,j)*f(i+1,j,nf)+aim(i,j)*f(i-1,j,nf) &
            +ajp(i,j)*f(i,j+1,nf)+ajm(i,j)*f(i,j-1,nf)-ap(i,j)*f(i,j,nf)
        end do

        denom = bl-pt(j-1)*blm

        if ( abs(denom/bl) < 1.0E-10 ) then
          pt(j) = 0.0D+00
          qt(j) = 0.0D+00
        else
          pt(j) = blp/denom
          qt(j) = (blc+blm*qt(j-1))/denom
        end if

      end do

      bl = 0.0D+00
      do j = m1-1, 2, -1
        bl = bl*pt(j)+qt(j)
        do i = 2, l1-1
          f(i,j,nf) = f(i,j,nf)+bl
        end do
      end do

    end if
!
!  Sweep 1: For each column J increasing.
!
    do j = 2, m1-1

      pt(1) = 0.0D+00
      qt(1) = f(1,j,nf)

      do i = 2, l1-1
        pt(i) = aip(i,j)/(ap(i,j)-pt(i-1)*aim(i,j))
        qt(i) = (con(i,j)+ajp(i,j)*f(i,j+1,nf)+ajm(i,j)*f(i,j-1,nf) &
          +aim(i,j)*qt(i-1))/(ap(i,j)-pt(i-1)*aim(i,j))
      end do

      do i = l1-1, 2, -1
        f(i,j,nf) = f(i+1,j,nf)*pt(i)+qt(i)
      end do

    end do
!
!  Sweep 2: For each column J decreasing.
!
    do j = m1-2, 2, -1

      pt(1) = 0.0D+00
      qt(1) = f(1,j,nf)

      do i = 2, l1-1
        pt(i) = aip(i,j)/(ap(i,j)-pt(i-1)*aim(i,j))
        qt(i) = (con(i,j)+ajp(i,j)*f(i,j+1,nf)+ajm(i,j)*f(i,j-1,nf) &
          +aim(i,j)*qt(i-1))/(ap(i,j)-pt(i-1)*aim(i,j))
      end do

      do i = l1-1, 2, -1
        f(i,j,nf) = f(i+1,j,nf)*pt(i)+qt(i)
      end do

    end do
!
!  Sweep 3: For each row I increasing.
!
    do i = 2, l1-1

      pt(1) = 0.0D+00
      qt(1) = f(i,1,nf)

      do j = 2, m1-1
        pt(j) = ajp(i,j)/(ap(i,j)-pt(j-1)*ajm(i,j))
        qt(j) = (con(i,j)+aip(i,j)*f(i+1,j,nf)+aim(i,j)*f(i-1,j,nf) &
          +ajm(i,j)*qt(j-1))/(ap(i,j)-pt(j-1)*ajm(i,j))
      end do

      do j = m1-1, 2, -1
        f(i,j,nf) = f(i,j+1,nf)*pt(j)+qt(j)
     end do

    end do
!
!  Sweep 4: For each row I decreasing.
!
    do i = l1-1, 2, -1

      pt(1) = 0.0D+00
      qt(1) = f(i,1,nf)

      do j = 2, m1-1
        pt(j) = ajp(i,j)/(ap(i,j)-pt(j-1)*ajm(i,j))
        qt(j) = (con(i,j)+aip(i,j)*f(i+1,j,nf)+aim(i,j)*f(i-1,j,nf) &
          +ajm(i,j)*qt(j-1))/(ap(i,j)-pt(j-1)*ajm(i,j))
      end do

      do j = m1-1, 2, -1
        f(i,j,nf) = f(i,j+1,nf)*pt(j)+qt(j)
      end do

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
subroutine assst ( iv, liv, lv, v )

!*******************************************************************************
!
!! ASSST assesses a candidate step.
!
!  Discussion:
!
!    This subroutine is called by an unconstrained minimization
!    routine to assess the next candidate step.  it may recommend one
!    of several courses of action, such as accepting the step, recom-
!    puting it using the same or a new quadratic model, or halting due
!    to convergence or false convergence.  See the return code listing
!    below.
!
!  Reference: 
!
!    John Dennis, David Gay, and Roy Welsch,
!    An Adaptive Nonlinear Least-squares Algorithm,
!    ACM Transactions on Mathematical Software, 
!    Volume 7, Number 3, 1981.
!
!    M J D Powell,
!    A Fortran Subroutine for Solving Systems of Nonlinear Algebraic Equations,
!    in Numerical Methods for Nonlinear Algebraic Equations, 
!    edited by Philip Rabinowitz, 
!    Gordon and Breach, London, 1970.
!
!  Parameters: 
!
!  iv (i/o) integer parameter and scratch vector -- see description
!             below of iv values referenced.
!
! liv (in)  length of iv array.
!
!  lv (in)  length of v array.
!
!   v (i/o) real parameter and scratch vector -- see description
!             below of v values referenced.
!
!   iv values referenced 
!
!    iv(irc) (i/o) on input for the first step tried in a new iteration,
!             iv(irc) should be set to 3 or 4 (the value to which it is
!             set when step is definitely to be accepted).  on input
!             after step has been recomputed, iv(irc) should be
!             unchanged since the previous return of assst.
!                on output, iv(irc) is a return code having one of the
!             following values...
!                  1 = switch models or try smaller step.
!                  2 = switch models or accept step.
!                  3 = accept step and determine v(radfac) by gradient
!                       tests.
!                  4 = accept step, v(radfac) has been determined.
!                  5 = recompute step (using the same model).
!                  6 = recompute step with radius = v(lmaxs) but do not
!                       evaulate the objective function.
!                  7 = x-convergence (see v(xctol)).
!                  8 = relative function convergence (see v(rfctol)).
!                  9 = both x- and relative function convergence.
!                 10 = absolute function convergence (see v(afctol)).
!                 11 = singular convergence (see v(lmaxs)).
!                 12 = false convergence (see v(xftol)).
!                 13 = iv(irc) was out of range on input.
!             return code i has precdence over i+1 for i = 9, 10, 11.
! iv(mlstgd) (i/o) saved value of iv(model).
!  iv(model) (i/o) on input, iv(model) should be an integer identifying
!             the current quadratic model of the objective function.
!             if a previous step yielded a better function reduction,
!             then iv(model) will be set to iv(mlstgd) on output.
! iv(nfcall) (in)  invocation count for the objective function.
! iv(nfgcal) (i/o) value of iv(nfcall) at step that gave the biggest
!             function reduction this iteration.  iv(nfgcal) remains
!             unchanged until a function reduction is obtained.
! iv(radinc) (i/o) the number of radius increases (or minus the number
!             of decreases) so far this iteration.
! iv(restor) (out) set to 1 if v(f) has been restored and x should be
!             restored to its initial value, to 2 if x should be saved,
!             to 3 if x should be restored from the saved value, and to
!             0 otherwise.
!  iv(stage) (i/o) count of the number of models tried so far in the
!             current iteration.
! iv(stglim) (in)  maximum number of models to consider.
! iv(switch) (out) set to 0 unless a new model is being tried and it
!             gives a smaller function value than the previous model,
!             in which case assst sets iv(switch) = 1.
! iv(toobig) (in)  is nonzero if step was too big (e.g. if it caused
!             overflow).
!   iv(xirc) (i/o) value that iv(irc) would have in the absence of
!             convergence, false convergence, and oversized steps.
!
!   v values referenced 
!
! v(afctol) (in)  absolute function convergence tolerance.  if the
!             absolute value of the current function value v(f) is less
!             than v(afctol), then assst returns with iv(irc) = 10.
! v(decfac) (in)  factor by which to decrease radius when iv(toobig) is
!             nonzero.
! v(dstnrm) (in)  the 2-norm of d*step.
! v(dstsav) (i/o) value of v(dstnrm) on saved step.
!   v(dst0) (in)  the 2-norm of d times the newton step (when defined,
!             i.e., for v(nreduc) >= 0).
!      v(f) (i/o) on both input and output, v(f) is the objective func-
!             tion value at x.  if x is restored to a previous value,
!             then v(f) is restored to the corresponding value.
!   v(fdif) (out) the function reduction v(f0) - v(f) (for the output
!             value of v(f) if an earlier step gave a bigger function
!             decrease, and for the input value of v(f) otherwise).
! v(flstgd) (i/o) saved value of v(f).
!     v(f0) (in)  objective function value at start of iteration.
! v(gtslst) (i/o) value of v(gtstep) on saved step.
! v(gtstep) (in)  inner product between step and gradient.
! v(incfac) (in)  minimum factor by which to increase radius.
!  v(lmaxs) (in)  maximum reasonable step size (and initial step bound).
!             if the actual function decrease is no more than twice
!             what was predicted, if a return with iv(irc) = 7, 8, 9,
!             or 10 does not occur, if v(dstnrm) > v(lmaxs), and if
!             v(preduc) <= v(sctol) * abs(v(f0)), then assst re-
!             turns with iv(irc) = 11.  if so doing appears worthwhile,
!             then assst repeats this test with v(preduc) computed for
!             a step of length v(lmaxs) (by a return with iv(irc) = 6).
! v(nreduc) (i/o)  function reduction predicted by quadratic model for
!             newton step.  if assst is called with iv(irc) = 6, i.e.,
!             if v(preduc) has been computed with radius = v(lmaxs) for
!             use in the singular convervence test, then v(nreduc) is
!             set to -v(preduc) before the latter is restored.
! v(plstgd) (i/o) value of v(preduc) on saved step.
! v(preduc) (i/o) function reduction predicted by quadratic model for
!             current step.
! v(radfac) (out) factor to be used in determining the new radius,
!             which should be v(radfac)*dst, where  dst  is either the
!             output value of v(dstnrm) or the 2-norm of
!             diag(newd)*step  for the output value of step and the
!             updated version, newd, of the scale vector d.  for
!             iv(irc) = 3, v(radfac) = 1.0D+00 is returned.
! v(rdfcmn) (in)  minimum value for v(radfac) in terms of the input
!             value of v(dstnrm) -- suggested value = 0.1.
! v(rdfcmx) (in)  maximum value for v(radfac) -- suggested value = 4.0.
!  v(reldx) (in) scaled relative change in x caused by step, computed
!             (e.g.) by function  reldst  as
!                 max (d(i)*abs(x(i)-x0(i)), 1 <= i <= p) /
!                    max (d(i)*(abs(x(i))+abs(x0(i))), 1 <= i <= p).
! v(rfctol) (in)  relative function convergence tolerance.  if the
!             actual function reduction is at most twice what was pre-
!             dicted and  v(nreduc) <= v(rfctol)*abs(v(f0)),  then
!             assst returns with iv(irc) = 8 or 9.
! v(stppar) (in)  marquardt parameter -- 0 means full newton step.
! v(tuner1) (in)  tuning constant used to decide if the function
!             reduction was much less than expected.  suggested
!             value = 0.1.
! v(tuner2) (in)  tuning constant used to decide if the function
!             reduction was large enough to accept step.  suggested
!             value = 10**-4.
! v(tuner3) (in)  tuning constant used to decide if the radius
!             should be increased.  suggested value = 0.75.
!  v(xctol) (in)  x-convergence criterion.  if step is a newton step
!             (v(stppar) = 0) having v(reldx) <= v(xctol) and giving
!             at most twice the predicted function decrease, then
!             assst returns iv(irc) = 7 or 9.
!  v(xftol) (in)  false convergence tolerance.  if step gave no or only
!             a small function decrease and v(reldx) <= v(xftol),
!             then assst returns with iv(irc) = 12.
!
!  notes  
!
!   application and usage restrictions 
!
!        this routine is called as part of the nl2sol (nonlinear
!     least-squares) package.  it may be used in any unconstrained
!     minimization solver that uses dogleg, goldfeld-quandt-trotter,
!     or levenberg-marquardt steps.
!
!   algorithm notes 
!
!        see (1) for further discussion of the assessing and model
!     switching strategies.  while nl2sol considers only two models,
!     assst is designed to handle any number of models.
!
!   usage notes 
!
!        on the first call of an iteration, only the i/o variables
!     step, x, iv(irc), iv(model), v(f), v(dstnrm), v(gtstep), and
!     v(preduc) need have been initialized.  between calls, no i/o
!     values execpt step, x, iv(model), v(f) and the stopping toler-
!     ances should be changed.
!        after a return for convergence or false convergence, one can
!     change the stopping tolerances and call assst again, in which
!     case the stopping tests will be repeated.
!
!   history 
!
!        john dennis designed much of this routine, starting with
!     ideas in (2). roy welsch suggested the model switching strategy.
!        david gay and stephen peters cast this subroutine into a more
!     portable form (winter 1977), and david gay cast it into its
!     present form (fall 1978).
!
  integer liv, lv
  integer iv(liv)
  real ( kind = 8 ) v(lv)
  logical goodx
  integer i, nfc
  real ( kind = 8 ) emax, emaxs, gts, rfac1, xmax
  real ( kind = 8 ) half, one, onep2, two
  integer afctol, decfac, dstnrm, dstsav, dst0, f, fdif, flstgd, f0
  integer gtslst, gtstep, incfac, irc, lmaxs, mlstgd, model, nfcall
  integer nfgcal, nreduc, plstgd, preduc, radfac, radinc, rdfcmn
  integer rdfcmx, reldx, restor, rfctol, sctol, stage, stglim
  integer stppar, switch, toobig, tuner1, tuner2, tuner3, xctol
  integer xftol, xirc

  parameter ( half=0.5d+0, one=1.d+0, onep2=1.2d+0, two=2.d+0)
  parameter ( irc=29, mlstgd=32, model=5, nfcall=6, nfgcal=7 )
  parameter ( radinc=8, restor=9, stage=10, stglim=11, switch=12 )
  parameter ( toobig=2, xirc=13)
  parameter (afctol=31, decfac=22, dstnrm=2, dst0=3, dstsav=18 )
  parameter (f=10, fdif=11, flstgd=12, f0=13, gtslst=14, gtstep=4 )
  parameter (incfac=23, lmaxs=36, nreduc=6, plstgd=15, preduc=7 )
  parameter (radfac=16, rdfcmn=24, rdfcmx=25, reldx=17, rfctol=32 )
  parameter (sctol=37, stppar=5, tuner1=26, tuner2=27, tuner3=28 )
  parameter (xctol=33, xftol=34)

  nfc = iv(nfcall)
  iv(switch) = 0
  iv(restor) = 0
  rfac1 = one
  goodx = .true.
  i = iv(irc)

  if (i >= 1 .and. i <= 12) then
        go to (20,30,10,10,40,280,220,220,220,220,220,170), i
  end if

  iv(irc) = 13
  return
!
!  Initialize for new iteration.
!
 10   iv(stage) = 1
  iv(radinc) = 0
  v(flstgd) = v(f0)
  if (iv(toobig) == 0) go to 110
     iv(stage) = -1
     iv(xirc) = i
     go to 60
!
!  Step was recomputed with new model or smaller radius 
!  first decide which 
!
 20   if (iv(model) /= iv(mlstgd)) go to 30
!
!  Old model retained, smaller radius tried 
!  do not consider any more new models this iteration 
!
     iv(stage) = iv(stglim)
     iv(radinc) = -1
     go to 110
!
!  A new model is being tried.  decide whether to keep it. 
!
 30   iv(stage) = iv(stage) + 1
!
!  Now we add the possibiltiy that step was recomputed with 
!  the same model, perhaps because of an oversized step.    
!
 40   if (iv(stage) > 0) go to 50
!
!  Step was recomputed because it was too big. 
!
     if (iv(toobig) /= 0) go to 60
!
!  Restore iv(stage) and pick up where we left off. 
!
     iv(stage) = -iv(stage)
     i = iv(xirc)
     go to (20, 30, 110, 110, 70), i

 50   if (iv(toobig) == 0) go to 70
!
!  Handle oversize step 
!
  if (iv(radinc) > 0) go to 80
     iv(stage) = -iv(stage)
     iv(xirc) = iv(irc)

 60      v(radfac) = v(decfac)
     iv(radinc) = iv(radinc) - 1
     iv(irc) = 5
     iv(restor) = 1
     return

 70   if (v(f) < v(flstgd)) go to 110
!
!  The new step is a loser.  restore old model. 
!
  if (iv(model) == iv(mlstgd)) go to 80
     iv(model) = iv(mlstgd)
     iv(switch) = 1
!
!  Restore step, etc. only if a previous step decreased v(f).
!
 80   if (v(flstgd) >= v(f0)) go to 110
     iv(restor) = 1
     v(f) = v(flstgd)
     v(preduc) = v(plstgd)
     v(gtstep) = v(gtslst)
     if (iv(switch) == 0) rfac1 = v(dstnrm) / v(dstsav)
     v(dstnrm) = v(dstsav)
     nfc = iv(nfgcal)
     goodx = .false.

 110  v(fdif) = v(f0) - v(f)
  if (v(fdif) > v(tuner2) * v(preduc)) go to 140
  if(iv(radinc)>0) go to 140
!
!         no (or only a trivial) function decrease
!         so try new model or smaller radius
!
     if (v(f) < v(f0)) go to 120
          iv(mlstgd) = iv(model)
          v(flstgd) = v(f)
          v(f) = v(f0)
          iv(restor) = 1
          go to 130
 120     iv(nfgcal) = nfc
 130     iv(irc) = 1
     if (iv(stage) < iv(stglim)) go to 160
          iv(irc) = 5
          iv(radinc) = iv(radinc) - 1
          go to 160
!
!  Nontrivial function decrease achieved 
!
 140  iv(nfgcal) = nfc
  rfac1 = 1.0D+00
  v(dstsav) = v(dstnrm)
  if (v(fdif) > v(preduc)*v(tuner1)) go to 190
!
!  Decrease was much less than predicted -- either change models
!  or accept step with decreased radius.
!
  if (iv(stage) >= iv(stglim)) go to 150
!
!  Consider switching models 
!
     iv(irc) = 2
     go to 160
!
!  Accept step with decreased radius 
!
 150  iv(irc) = 4
!
!   set v(radfac) to fletcher*s decrease factor 
!
 160  iv(xirc) = iv(irc)
  emax = v(gtstep) + v(fdif)
  v(radfac) = half * rfac1

  if (emax < v(gtstep)) then
    v(radfac) = rfac1 * max (v(rdfcmn),half * v(gtstep)/emax)
  end if
!
!  Do false convergence test 
!
 170  if (v(reldx) <= v(xftol)) go to 180
     iv(irc) = iv(xirc)
     if (v(f) < v(f0)) go to 200
          go to 230

 180  iv(irc) = 12
  go to 240
!
!  Handle good function decrease 
!
 190  if (v(fdif) < (-v(tuner3) * v(gtstep))) go to 210
!
!  Increasing radius looks worthwhile.  see if we just
!  recomputed step with a decreased radius or restored step
!  after recomputing it with a larger radius.
!
  if (iv(radinc) < 0) go to 210
  if (iv(restor) == 1) go to 210
!
!  We did not.  try a longer step unless this was a newton step.
!
     v(radfac) = v(rdfcmx)
     gts = v(gtstep)
     if (v(fdif) < (half/v(radfac) - 1.0D+00 ) * gts) then
       v(radfac) = max (v(incfac), half*gts/(gts + v(fdif)))
     end if
     iv(irc) = 4
     if (v(stppar) == 0.0D+00 ) go to 230
     if (v(dst0) >= 0.0D+00 .and. (v(dst0) < two*v(dstnrm) &
              .or. v(nreduc) < onep2*v(fdif)))  then
       go to 230
     end if
!
!  Step was not a newton step.  recompute it with a larger radius.
!
          iv(irc) = 5
          iv(radinc) = iv(radinc) + 1
!
!  Save values corresponding to good step 
!
 200  v(flstgd) = v(f)
  iv(mlstgd) = iv(model)
  if (iv(restor) /= 1) iv(restor) = 2
  v(dstsav) = v(dstnrm)
  iv(nfgcal) = nfc
  v(plstgd) = v(preduc)
  v(gtslst) = v(gtstep)
  go to 230
!
!  Accept step with radius unchanged.
!
 210  v(radfac) = 1.0D+00
  iv(irc) = 3
  go to 230
!
!  Come here for a restart after convergence.
!
 220  iv(irc) = iv(xirc)
  if (v(dstsav) >= 0.0D+00 ) go to 240
     iv(irc) = 12
     go to 240
!
!  Perform convergence tests.
!
 230  iv(xirc) = iv(irc)
 240  if (iv(restor) == 1 .and. v(flstgd) < v(f0)) iv(restor) = 3
  if (abs(v(f)) < v(afctol)) iv(irc) = 10

  if (half * v(fdif) > v(preduc)) then
    return
  end if

  emax = v(rfctol) * abs(v(f0))
  emaxs = v(sctol) * abs(v(f0))
  if (v(dstnrm) > v(lmaxs) .and. v(preduc) <= emaxs) then
    iv(irc) = 11
  end if
  if (v(dst0) < 0.0D+00 ) go to 250
  i = 0

  if ((v(nreduc) > 0.0D+00 .and. v(nreduc) <= emax) .or. &
      (v(nreduc) == 0.0D+00 .and. v(preduc) == 0.0D+00 )) then
    i = 2
  end if

  if (v(stppar) == 0.0D+00 .and. v(reldx) <= v(xctol) .and. goodx) then
    i = i + 1
  end if

  if (i > 0) iv(irc) = i + 6
!
!  Consider recomputing step of length v(lmaxs) for singular
!  convergence test.
!
 250  if (iv(irc) > 5 .and. iv(irc) /= 12) then
     return
  end if

  if (v(dstnrm) > v(lmaxs)) go to 260
     if (v(preduc) >= emaxs) then
       return
     end if
          if (v(dst0) <= 0.0D+00 ) go to 270
               if (half * v(dst0) <= v(lmaxs)) then
                 return
               end if
                    go to 270
 260  if (half * v(dstnrm) <= v(lmaxs)) then
        return
      end if
  xmax = v(lmaxs) / v(dstnrm)
  if (xmax * (two - xmax) * v(preduc) >= emaxs) then
    return
  end if
 270  if (v(nreduc) < 0.0D+00 ) go to 290
!
!   recompute v(preduc) for use in singular convergence test 
!
  v(gtslst) = v(gtstep)
  v(dstsav) = v(dstnrm)
  if (iv(irc) == 12) v(dstsav) = -v(dstsav)
  v(plstgd) = v(preduc)
  i = iv(restor)
  iv(restor) = 2
  if (i == 3) iv(restor) = 0
  iv(irc) = 6
  return
!
!  Perform singular convergence test with recomputed v(preduc) 
!
 280  v(gtstep) = v(gtslst)
  v(dstnrm) = abs(v(dstsav))
  iv(irc) = iv(xirc)
  if (v(dstsav) <= 0.0D+00 ) iv(irc) = 12
  v(nreduc) = -v(preduc)
  v(preduc) = v(plstgd)
  iv(restor) = 3

 290  if (-v(nreduc) <= v(rfctol) * abs(v(f0))) iv(irc) = 11

  return
end
subroutine dbdog ( dig, lv, n, nwtstp, step, v )

!*******************************************************************************
!
!! DBDOG: compute a double dogleg step.
!
!  Discussion:
!
!    This subroutine computes a candidate step (for use in an 
!    unconstrained minimization code) by the double dogleg algorithm of
!    dennis and mei (ref. 1), which is a variation on powell*s dogleg
!    scheme (ref. 2, p. 95).
!
!    let  g  and  h  be the current gradient and hessian approxima-
!    tion respectively and let d be the current scale vector.  this
!    routine assumes dig = diag(d)**-2 * g  and  nwtstp = h**-1 * g.
!    the step computed is the same one would get by replacing g and h
!    by  diag(d)**-1 * g  and  diag(d)**-1 * h * diag(d)**-1,
!    computing step, and translating step back to the original
!    variables, i.e., premultiplying it by diag(d)**-1.
!
!  Reference:
!
!    John Dennis, Howell Mei,
!    Two New Unconstrained Optimization Algorithms Which Use 
!    Function and Gradient Values, 
!    Journal of Optimization Theory and Applications,
!    Volume 28, pages 453-482, 1979.
!
!    M J D Powell,
!    A Hybrid Method for Non-linear Equations,
!    in Numerical Methods for Non-linear Equations, 
!    edited by Philip Rabinowitz, 
!    Gordon and Breach, London, 1970.
!
!  Parameters: 
!
!    dig (input) diag(d)**-2 * g -- see algorithm notes.
!      g (input) the current gradient vector.
!     lv (input) length of v.
!      n (input) number of components in  dig, g, nwtstp,  and  step.
! nwtstp (input) negative newton step -- see algorithm notes.
!   step (output) the computed step.
!      v (i/o) values array, the following components of which are
!             used here...
! v(bias)   (input) bias for relaxed newton step, which is v(bias) of
!             the way from the full newton to the fully relaxed newton
!             step.  recommended value = 0.8 .
! v(dgnorm) (input) 2-norm of diag(d)**-1 * g -- see algorithm notes.
! v(dstnrm) (output) 2-norm of diag(d) * step, which is v(radius)
!             unless v(stppar) = 0 -- see algorithm notes.
! v(dst0) (input) 2-norm of diag(d) * nwtstp -- see algorithm notes.
! v(grdfac) (output) the coefficient of  dig  in the step returned --
!             step(i) = v(grdfac)*dig(i) + v(nwtfac)*nwtstp(i).
! v(gthg)   (input) square-root of (dig**t) * (hessian) * dig -- see
!             algorithm notes.
! v(gtstep) (output) inner product between g and step.
! v(nreduc) (output) function reduction predicted for the full newton
!             step.
! v(nwtfac) (output) the coefficient of  nwtstp  in the step returned --
!             see v(grdfac) above.
! v(preduc) (output) function reduction predicted for the step returned.
! v(radius) (input) the trust region radius.  d times the step returned
!             has 2-norm v(radius) unless v(stppar) = 0.
! v(stppar) (output) code telling how step was computed... 0 means a
!             full newton step.  between 0 and 1 means v(stppar) of the
!             way from the newton to the relaxed newton step.  between
!             1 and 2 means a true double dogleg step, v(stppar) - 1 of
!             the way from the relaxed newton to the Cauchy step.
!             greater than 2 means 1 / (v(stppar) - 1) times the Cauchy
!             step.
!
  integer lv
  integer n

  real ( kind = 8 ) dig(n), nwtstp(n), step(n), v(lv)
  external dotprd, v2norm
  real ( kind = 8 ) dotprd, v2norm
  real ( kind = 8 ) cfact, cnorm, ctrnwt, ghinvg, femnsq, gnorm
  real ( kind = 8 ) nwtnrm, relax, rlambd, t, t1, t2
  real ( kind = 8 ) half, two
  integer bias, dgnorm, dstnrm, dst0, grdfac, gthg, gtstep
  integer nreduc, nwtfac, preduc, radius, stppar
  parameter (half=0.5d+0, two=2.d+0)
  parameter (bias=43, dgnorm=1, dstnrm=2, dst0=3, grdfac=45 )
  parameter ( gthg=44, gtstep=4, nreduc=6, nwtfac=46, preduc=7 )
  parameter ( radius=8, stppar=5)

  nwtnrm = v(dst0)
  rlambd = 1.0D+00
  if (nwtnrm > 0.0D+00 ) rlambd = v(radius) / nwtnrm
  gnorm = v(dgnorm)
  ghinvg = two * v(nreduc)
  v(grdfac) = 0.0D+00
  v(nwtfac) = 0.0D+00
  if (rlambd < 1.0D+00 ) go to 30
!
!  The Newton step is inside the trust region.
!
     v(stppar) = 0.0D+00
     v(dstnrm) = nwtnrm
     v(gtstep) = -ghinvg
     v(preduc) = v(nreduc)
     v(nwtfac) = -1.0D+00
     step(1:n) = -nwtstp(1:n)
     return

 30   v(dstnrm) = v(radius)
  cfact = (gnorm / v(gthg))**2
!
!  Cauchy step = -cfact * g.
!
  cnorm = gnorm * cfact
  relax = 1.0D+00 - v(bias) * ( 1.0D+00 - gnorm*cnorm/ghinvg)
  if (rlambd < relax) go to 50
!
!  Step is between relaxed Newton and full Newton steps.
!
     v(stppar) = 1.0D+00 -  (rlambd - relax) / ( 1.0D+00 - relax)
     t = -rlambd
     v(gtstep) = t * ghinvg
     v(preduc) = rlambd * ( 1.0D+00 - half*rlambd) * ghinvg
     v(nwtfac) = t
     step(1:n) = t * nwtstp(1:n)
     return

 50   if (cnorm < v(radius)) go to 70
!
!  The Cauchy step lies outside the trust region --
!  step = scaled Cauchy step.
!
     t = -v(radius) / gnorm
     v(grdfac) = t
     v(stppar) = 1.0D+00  +  cnorm / v(radius)
     v(gtstep) = -v(radius) * gnorm
  v(preduc) = v(radius)*(gnorm - half*v(radius)*(v(gthg)/gnorm)**2)
     step(1:n) = t * dig(1:n)
     return
!
!  Compute dogleg step between Cauchy and relaxed Newton 
!  femur = relaxed newton step minus Cauchy step.
!
 70   ctrnwt = cfact * relax * ghinvg / gnorm
!
!  ctrnwt = inner product of Cauchy and relaxed Newton steps,
!  scaled by gnorm**-1.
!
  t1 = ctrnwt - gnorm*cfact**2
!
!  t1 = inner prod. of femur and Cauchy step, scaled by gnorm**-1.
!
  t2 = v(radius)*(v(radius)/gnorm) - gnorm*cfact**2
  t = relax * nwtnrm
  femnsq = (t/gnorm)*t - ctrnwt - t1
!
!  femnsq = square of 2-norm of femur, scaled by gnorm**-1.
!
  t = t2 / (t1 + sqrt(t1**2 + femnsq*t2))
!
!  Dogleg step  =  Cauchy step  +  t * femur.
!
  t1 = (t - 1.0D+00 ) * cfact
  v(grdfac) = t1
  t2 = -t * relax
  v(nwtfac) = t2
  v(stppar) = two - t
  v(gtstep) = t1*gnorm**2 + t2*ghinvg
  v(preduc) = -t1*gnorm * ((t2 + 1.0D+00 )*gnorm) &
                  - t2 * ( 1.0D+00 + half*t2)*ghinvg &
                   - half * (v(gthg)*t1)**2

  step(1:n) = t1 * dig(1:n) + t2 * nwtstp(1:n)

  return
end
subroutine deflt ( alg, iv, liv, lv, v )

!*******************************************************************************
!
!! DEFLT: supply default values to IV and V.
!
!  Discussion:
!
!   ALG = 1 means regression constants.
!   ALG = 2 means general unconstrained optimization constants.
!
  integer liv
  integer lv

  integer alg
  integer iv(liv)
  real ( kind = 8 ) v(lv)
  external vdflt
  integer miv, mv
  integer miniv(2), minv(2)
  integer algsav, covprt, covreq, dtype, hc, ierr, inith, inits
  integer ipivot, ivneed, lastiv, lastv, lmat, mxfcal, mxiter
  integer nfcov, ngcov, nvdflt, outlev, parprt, parsav, perm
  integer prunit, qrtyp, rdreq, rmat, solprt, statpr, vneed
  integer vsave, x0prt

  parameter (algsav=51, covprt=14, covreq=15, dtype=16, hc=71 )
  parameter (ierr=75, inith=25, inits=25, ipivot=76, ivneed=3 )
  parameter (lastiv=44, lastv=45, lmat=42, mxfcal=17, mxiter=18 )
  parameter (nfcov=52, ngcov=53, nvdflt=50, outlev=19, parprt=20 )
  parameter (parsav=49, perm=58, prunit=21, qrtyp=80, rdreq=57 )
  parameter (rmat=78, solprt=22, statpr=23, vneed=4, vsave=60 )
  parameter (x0prt=24)

  data miniv(1)/80/, miniv(2)/59/, minv(1)/98/, minv(2)/71/

  if ( alg < 1 .or. 2 < alg ) then
    iv(1) = 67
    return
  end if

  miv = miniv(alg)

  if ( liv < miv ) then
    iv(1) = 15
    return
  end if

  mv = minv(alg)

  if ( lv < mv ) then
    iv(1) = 16
    return
  end if

  call vdflt(alg, lv, v)
  iv(1) = 12
  iv(algsav) = alg
  iv(ivneed) = 0
  iv(lastiv) = miv
  iv(lastv) = mv
  iv(lmat) = mv + 1
  iv(mxfcal) = 200
  iv(mxiter) = 150
  iv(outlev) = 1
  iv(parprt) = 1
  iv(perm) = miv + 1
  iv(prunit) = 6
  iv(solprt) = 1
  iv(statpr) = 1
  iv(vneed) = 0
  iv(x0prt) = 1
!
!  General optimization values.
!
  if ( 2 <= alg ) then

    iv(dtype) = 0
    iv(inith) = 1
    iv(nfcov) = 0
    iv(ngcov) = 0
    iv(nvdflt) = 25
    iv(parsav) = 47
!
!  Regression values.
!
  else

    iv(covprt) = 3
    iv(covreq) = 1
    iv(dtype) = 1
    iv(hc) = 0
    iv(ierr) = 0
    iv(inits) = 0
    iv(ipivot) = 0
    iv(nvdflt) = 32
    iv(parsav) = 67
    iv(qrtyp) = 1
    iv(rdreq) = 3
    iv(rmat) = 0
    iv(vsave) = 58

  end if

  return
end
function dotprd ( p, x, y )

!*******************************************************************************
!
!! DOTPRD returns the inner product of vectors X and Y. 
!
  integer p

  real ( kind = 8 ) dotprd
  integer i
  real ( kind = 8 ) rmdcon
  real ( kind = 8 ), save :: sqteta = 0.0D+00
  real ( kind = 8 ) t
  real ( kind = 8 ) x(p)
  real ( kind = 8 ) y(p)

  dotprd = 0.0D+00

  if ( sqteta == 0.0D+00 ) then
    sqteta = rmdcon(2)
  end if

  do i = 1, p

    t = max ( abs ( x(i) ), abs ( y(i) ) )
    if ( t > 1.0D+00 ) go to 10
    if (t < sqteta) go to 20
    t = (x(i)/sqteta)*y(i)
    if (abs(t) < sqteta) go to 20
 10   dotprd = dotprd + x(i)*y(i)

 20 continue

  end do

  return
end
subroutine dupdu ( d, hdiag, iv, liv, lv, n, v )

!*******************************************************************************
!
!! DUPDU: update scale vector D for HUMSL.
!
!  Modified:
!
!    20 February 2006
!
  integer liv
  integer lv
  integer n

  real ( kind = 8 ) d(n)
  integer d0i
  integer, parameter :: dfac = 41
  integer, parameter :: dtol = 59
  integer, parameter :: dtype = 16
  integer dtoli
  real ( kind = 8 ) hdiag(n)
  integer i
  integer iv(liv)
  integer, parameter :: niter = 31
  real ( kind = 8 ) t
  real ( kind = 8 ) v(lv)
  real ( kind = 8 ) vdfac

  i = iv(dtype)

  if ( i /= 1 ) then
    if ( 0 < iv(niter) ) then
      return
    end if
  end if
   
  dtoli = iv(dtol)
  d0i = dtoli + n
  vdfac = v(dfac)

  do i = 1, n

    t = max ( sqrt ( abs ( hdiag(i) ) ), vdfac * d(i) )

    if ( t < v(dtoli) ) then
      t = max ( v(dtoli), v(d0i) )
    end if

    d(i) = t
    dtoli = dtoli + 1
    d0i = d0i + 1

  end do

  return
end
subroutine gqtst ( d, dig, dihdi, ka, l, p, step, v, w )

!*******************************************************************************
!
!! GQTST: compute Goldfeld-Quandt-Trotter step by More-Hebden technique.
!
!  Discussion:
!
!    Given the (compactly stored) lower triangle of a scaled
!    hessian (approximation) and a nonzero scaled gradient vector,
!    this subroutine computes a goldfeld-quandt-trotter step of
!    approximate length v(radius) by the more-hebden technique.  in
!    other words, step is computed to (approximately) minimize
!    psi(step) = (g**t)*step + 0.5*(step**t)*h*step  such that the
!    2-norm of d*step is at most (approximately) v(radius), where
!    g  is the gradient,  h  is the hessian, and  d  is a diagonal
!    scale matrix whose diagonal is stored in the parameter d.
!    (gqtst assumes  dig = d**-1 * g  and  dihdi = d**-1 * h * d**-1.)
!
!    the desired g-q-t step (ref. 2, 3, 4, 6) satisfies
!    (h + alpha*d**2)*step = -g  for some nonnegative alpha such that
!    h + alpha*d**2 is positive semidefinite.  alpha and step are
!    computed by a scheme analogous to the one described in ref. 5.
!    estimates of the smallest and largest eigenvalues of the hessian
!    are obtained from the gerschgorin circle theorem enhanced by a
!    simple form of the scaling described in ref. 7.  cases in which
!    h + alpha*d**2 is nearly (or exactly) singular are handled by
!    the technique discussed in ref. 2.  in these cases, a step of
!    (exact) length v(radius) is returned for which psi(step) exceeds
!    its optimal value by less than -v(epslon)*psi(step).  the test
!    suggested in ref. 6 for detecting the special case is performed
!    once two matrix factorizations have been done -- doing so sooner
!    seems to degrade the performance of optimization routines that
!    call this routine.
!
!  Reference:
!
!    John Dennis, David Gay, Roy Welsch,
!    An Adaptive Nonlinear Least-squares Algorithm,
!    ACM Transactions on Mathematical Software, 
!    Volume 7, Number 3, 1981.
!
!    David Gay,
!    Computing Optimal Locally Constrained Steps,
!    SIAM Journal on Scientific and Statistical Computing,
!    Volume 2, pages 186-197, 1981.
!
!    S M Goldfeld, R E Quandt, H F Trotter,
!    Maximization by Quadratic Hill-climbing, 
!    Econometrica,
!    Volume 34, pages 541-551, 1966.
!
!    M D Hebden,
!    An Algorithm for Minimization Using Exact Second Derivatives,
!    Report TP 515, Theoretical Physics Division,
!    AERE Harwell, Oxon., England, 1973.
!
!    Jorge More,
!    The Levenberg-Marquardt Algorithm, Implementation and Theory, 
!    Springer Lecture Notes in Mathematics Number 630, pages 105-116,
!    edited by G A Watson, 
!    Springer-Verlag, Berlin and New York, 1978.
!
!    Jorge More and Danny Sorensen,
!    Computing a Trust Region Step, 
!    Technical Report ANL-81-83, 
!    Argonne National Lab, 1981.
!
!    Richard Varga,
!    Minimal Gerschgorin Sets, 
!    Pacific Journal of Mathematics,
!    Volume 15, pages 719-729, 1965.
!
!  Parameters: 
!
!     d (in)  = the scale vector, i.e. the diagonal of the scale
!              matrix  d  mentioned above under purpose.
!   dig (in)  = the scaled gradient vector, d**-1 * g.  if g = 0, then
!              step = 0  and  v(stppar) = 0  are returned.
! dihdi (in)  = lower triangle of the scaled hessian (approximation),
!              i.e., d**-1 * h * d**-1, stored compactly by rows., i.e.,
!              in the order (1,1), (2,1), (2,2), (3,1), (3,2), etc.
!    ka (i/o) = the number of hebden iterations (so far) taken to deter-
!              mine step.  ka < 0 on input means this is the first
!              attempt to determine step (for the present dig and dihdi)
!              -- ka is initialized to 0 in this case.  output with
!              ka = 0  (or v(stppar) = 0)  means  step = -(h**-1)*g.
!     l (i/o) = workspace of length p*(p+1)/2 for cholesky factors.
!     p (in)  = number of parameters -- the hessian is a  p x p  matrix.
!  step (i/o) = the step computed.
!     v (i/o) contains various constants and variables described below.
!     w (i/o) = workspace of length 4*p + 6.
!
!   entries in v 
!
! v(dgnorm) (i/o) = 2-norm of (d**-1)*g.
! v(dstnrm) (output) = 2-norm of d*step.
! v(dst0)   (i/o) = 2-norm of d*(h**-1)*g (for pos. def. h only), or
!             overestimate of smallest eigenvalue of (d**-1)*h*(d**-1).
! v(epslon) (in)  = max. rel. error allowed for psi(step).  for the
!             step returned, psi(step) will exceed its optimal value
!             by less than -v(epslon)*psi(step).  suggested value = 0.1.
! v(gtstep) (out) = inner product between g and step.
! v(nreduc) (out) = psi(-(h**-1)*g) = psi(newton step)  (for pos. def.
!             h only -- v(nreduc) is set to zero otherwise).
! v(phmnfc) (in)  = tol. (together with v(phmxfc)) for accepting step
!             (more*s sigma).  the error v(dstnrm) - v(radius) must lie
!             between v(phmnfc)*v(radius) and v(phmxfc)*v(radius).
! v(phmxfc) (in)  (see v(phmnfc).)
!             suggested values -- v(phmnfc) = -0.25, v(phmxfc) = 0.5.
! v(preduc) (out) = psi(step) = predicted obj. func. reduction for step.
! v(radius) (in)  = radius of current (scaled) trust region.
! v(rad0)   (i/o) = value of v(radius) from previous call.
! v(stppar) (i/o) is normally the marquardt parameter, i.e. the alpha
!             described below under algorithm notes.  if h + alpha*d**2
!             (see algorithm notes) is (nearly) singular, however,
!             then v(stppar) = -alpha.
!
!   usage notes 
!
!     if it is desired to recompute step using a different value of
!     v(radius), then this routine may be restarted by calling it
!     with all parameters unchanged except v(radius).  (this explains
!     why step and w are listed as i/o).  on an initial call (one with
!     ka < 0), step and w need not be initialized and only compo-
!     nents v(epslon), v(stppar), v(phmnfc), v(phmxfc), v(radius), and
!     v(rad0) of v must be initialized.
!
  integer ka, p
  real ( kind = 8 ) d(p), dig(p), dihdi(*), l(*), v(21), step(p), w(*)
!     dimension dihdi(p*(p+1)/2), l(p*(p+1)/2), w(4*p+7)
!
  logical restrt
  integer dggdmx, diag, diag0, dstsav, emax, emin, i, im1, inc, irc
  integer j, k, kalim, kamin, k1, lk0, phipin, q, q0, uk0, x
  real ( kind = 8 ) alphak, aki, akk, delta, dst, eps, gtsta, lk
  real ( kind = 8 ) oldphi
  real ( kind = 8 ) phi, phimax, phimin, psifac, rad, radsq
  real ( kind = 8 ) root, si, sk, sw, t, twopsi, t1, t2, uk, wi
  real ( kind = 8 ) big, dgxfac, epsfac, four, half, kappa, negone
  real ( kind = 8 ) one, p001, six, three, two, zero
  real ( kind = 8 ) dotprd, lsvmin, rmdcon, v2norm
  integer dgnorm, dstnrm, dst0, epslon, gtstep, stppar, nreduc
  integer phmnfc, phmxfc, preduc, radius, rad0

  parameter (dgnorm=1, dstnrm=2, dst0=3, epslon=19, gtstep=4 )
  parameter ( nreduc=6, phmnfc=20, phmxfc=21, preduc=7, radius=8 )
  parameter ( rad0=9, stppar=5)

  parameter (epsfac=50.0d+0, four=4.0d+0, half=0.5d+0 )
  parameter ( kappa=2.0d+0, negone=-1.0d+0, one=1.0d+0, p001=1.0d-3 )
  parameter ( six=6.0d+0, three=3.0d+0, two=2.0d+0, zero=0.0d+0)

  save dgxfac

  data big/0.d+0/, dgxfac/0.d+0/
!
!  Store largest abs. entry in (d**-1)*h*(d**-1) at w(dggdmx).
!
  dggdmx = p + 1
!
!  Store Gerschgorin over- and underestimates of the largest
!  and smallest eigenvalues of (d**-1)*h*(d**-1) at w(emax)
!  and w(emin) respectively.
!
  emax = dggdmx + 1
  emin = emax + 1
!
!  For use in recomputing step, the final values of lk, uk, dst,
!  and the inverse derivative of more*s phi at 0 (for pos. def.
!  h) are stored in w(lk0), w(uk0), w(dstsav), and w(phipin)
!  respectively.
!
  lk0 = emin + 1
  phipin = lk0 + 1
  uk0 = phipin + 1
  dstsav = uk0 + 1
!
!  Store diag of (d**-1)*h*(d**-1) in w(diag),...,w(diag0+p).
!
  diag0 = dstsav
  diag = diag0 + 1
!
!  Store -d*step in w(q),...,w(q0+p).
!
  q0 = diag0 + p
  q = q0 + 1
!
!  Allocate storage for scratch vector x 
!
  x = q + p
  rad = v(radius)
  radsq = rad**2
!
!  phitol = max. error allowed in dst = v(dstnrm) = 2-norm of d*step.
!
  phimax = v(phmxfc) * rad
  phimin = v(phmnfc) * rad
  psifac = two * v(epslon) / (three * (four * (v(phmnfc) + 1.0D+00 ) * &
                    (kappa + 1.0D+00 )  +  kappa  +  two) * rad**2)
!
!  OLDPHI is used to detect limits of numerical accuracy.  if
!  we recompute step and it does not change, then we accept it.
!
  oldphi = 0.0D+00
  eps = v(epslon)
  irc = 0
  restrt = .false.
  kalim = ka + 50
!
!  Start or restart, depending on ka 
!
  if (ka >= 0) go to 290
!
!  fresh start 
!
  k = 0
  uk = negone
  ka = 0
  kalim = 50
  v(dgnorm) = v2norm(p, dig)
  v(nreduc) = 0.0D+00
  v(dst0) = 0.0D+00
  kamin = 3
  if (v(dgnorm) == 0.0D+00 ) kamin = 0
!
!  store diag(dihdi) in w(diag0+1),...,w(diag0+p) 
!
  j = 0
  do i = 1, p
    j = j + i
    k1 = diag0 + i
    w(k1) = dihdi(j)
  end do
!
!  determine w(dggdmx), the largest element of dihdi 
!
  t1 = 0.0D+00
  j = p * (p + 1) / 2
  do i = 1, j
     t = abs(dihdi(i))
     if (t1 < t) t1 = t
  end do
  w(dggdmx) = t1
!
!  try alpha = 0 
!
 30   call lsqrt(1, p, l, dihdi, irc)
  if (irc == 0) go to 50
!
!  indefinite h -- underestimate smallest eigenvalue, use this
!  estimate to initialize lower bound lk on alpha.
!
     j = irc*(irc+1)/2
     t = l(j)
     l(j) = 1.0D+00
     w(1:irc) = 0.0D+00
     w(irc) = one
     call litvmu(irc, w, l, w)
     t1 = v2norm(irc, w)
     lk = -t / t1 / t1
     v(dst0) = -lk
     if (restrt) go to 210
     go to 70
!
!  positive definite h -- compute unmodified newton step. 
!
 50   lk = 0.0D+00
  t = lsvmin(p, l, w(q), w(q))
  if (t >= one) go to 60
     if (big <= 0.0D+00 ) big = rmdcon(6)
     if (v(dgnorm) >= t*t*big) go to 70
 60   call livmul(p, w(q), l, dig)
  gtsta = dotprd(p, w(q), w(q))
  v(nreduc) = half * gtsta
  call litvmu(p, w(q), l, w(q))
  dst = v2norm(p, w(q))
  v(dst0) = dst
  phi = dst - rad
  if (phi <= phimax) go to 260
  if (restrt) go to 210
!
!  Prepare to compute Gerschgorin estimates of largest (and
!  smallest) eigenvalues. 
!
 70   k = 0
  do i = 1, p
     wi = 0.0D+00
     im1 = i - 1
     do j = 1, im1
       k = k + 1
       t = abs(dihdi(k))
       wi = wi + t
       w(j) = w(j) + t
     end do
     w(i) = wi
     k = k + 1
  end do
!
!  (under-)estimate smallest eigenvalue of (d**-1)*h*(d**-1) 
!
  k = 1
  t1 = w(diag) - w(1)

  do i = 2, p
     j = diag0 + i
     t = w(j) - w(i)
     if ( t < t1 ) then
       t1 = t
       k = i
     end if
  end do

  sk = w(k)
  j = diag0 + k
  akk = w(j)
  k1 = k*(k-1)/2 + 1
  inc = 1
  t = 0.0D+00
  do i = 1, p
     if (i == k) go to 130
     aki = abs(dihdi(k1))
     si = w(i)
     j = diag0 + i
     t1 = half * (akk - w(j) + si - aki)
     t1 = t1 + sqrt(t1*t1 + sk*aki)
     if (t < t1) t = t1
     if (i < k) go to 140
 130     inc = i
 140     k1 = k1 + inc
  end do

  w(emin) = akk - t
  uk = v(dgnorm)/rad - w(emin)
  if (v(dgnorm) == 0.0D+00 ) uk = uk + p001 + p001*uk
  if (uk <= 0.0D+00) uk = p001
!
!   compute Gerschgorin overestimate of largest eigenvalue 
!
  k = 1
  t1 = w(diag) + w(1)
  if (p <= 1) go to 170

  do i = 2, p
     j = diag0 + i
     t = w(j) + w(i)
     if (t <= t1) go to 160
          t1 = t
          k = i
 160     continue
  end do

 170  sk = w(k)
  j = diag0 + k
  akk = w(j)
  k1 = k*(k-1)/2 + 1
  inc = 1
  t = 0.0D+00

  do i = 1, p
     if (i == k) go to 180
     aki = abs(dihdi(k1))
     si = w(i)
     j = diag0 + i
     t1 = half * (w(j) + si - aki - akk)
     t1 = t1 + sqrt(t1*t1 + sk*aki)
     if (t < t1) t = t1
     if (i < k) go to 190
 180     inc = i
 190     k1 = k1 + inc
  end do

  w(emax) = akk + t
  lk = max (lk, v(dgnorm)/rad - w(emax))
!
!  alphak = current value of alpha (see alg. notes above).  we
!  use More's scheme for initializing it.
!
  alphak = abs(v(stppar)) * v(rad0)/rad

  if (irc /= 0) go to 210
!
!  Compute l0 for positive definite H. 
!
  call livmul(p, w, l, w(q))
  t = v2norm(p, w)
  w(phipin) = dst / t / t
  lk = max (lk, phi*w(phipin))
!
!  safeguard alphak and add alphak*i to (d**-1)*h*(d**-1) 
!
 210  ka = ka + 1
  if (-v(dst0) >= alphak .or. alphak < lk .or. alphak >= uk) then
    alphak = uk * max (p001, sqrt(lk/uk))
  end if
  if (alphak <= 0.0D+00) alphak = half * uk
  if (alphak <= 0.0D+00) alphak = uk
  k = 0
  do i = 1, p
     k = k + i
     j = diag0 + i
     dihdi(k) = w(j) + alphak
  end do
!
!  Try computing Cholesky decomposition 
!
  call lsqrt(1, p, l, dihdi, irc)
  if (irc == 0) go to 240
!
!  (d**-1)*h*(d**-1) + alphak*i  is indefinite -- overestimate
!  smallest eigenvalue for use in updating lk 
!
  j = (irc*(irc+1))/2
  t = l(j)
  l(j) = one
  w(1:irc) = 0.0D+00
  w(irc) = one
  call litvmu(irc, w, l, w)
  t1 = v2norm(irc, w)
  lk = alphak - t/t1/t1
  v(dst0) = -lk
  go to 210
!
!  Alphak makes (d**-1)*h*(d**-1) positive definite.
!  compute q = -d*step, check for convergence. 
!
 240  call livmul(p, w(q), l, dig)
  gtsta = dotprd(p, w(q), w(q))
  call litvmu(p, w(q), l, w(q))
  dst = v2norm(p, w(q))
  phi = dst - rad
  if (phi <= phimax .and. phi >= phimin) go to 270
  if (phi == oldphi) go to 270
  oldphi = phi
  if (phi < 0.0D+00) go to 330
!
!  unacceptable alphak -- update lk, uk, alphak 
!
 250  if (ka >= kalim) go to 270
!
!  The following dmin1 is necessary because of restarts 
!
  if (phi < 0.0D+00) uk = min (uk, alphak)
!
!  kamin = 0 only iff the gradient vanishes 
!
  if (kamin == 0) go to 210
  call livmul(p, w, l, w(q))
  t1 = v2norm(p, w)
  alphak = alphak  +  (phi/t1) * (dst/t1) * (dst/rad)
  lk = max (lk, alphak)
  go to 210
!
!  Acceptable step on first try.
!
 260  alphak = 0.0D+00
!
!  Successful step in general.  compute step = -(d**-1)*q 
!
 270  continue

  do i = 1, p
    j = q0 + i
    step(i) = -w(j)/d(i)
  end do

  v(gtstep) = -gtsta
  v(preduc) = half * (abs(alphak)*dst*dst + gtsta)
  go to 410
!
!  Restart with new radius 
!
 290  if (v(dst0) <= 0.0D+00 .or. v(dst0) - rad > phimax) go to 310
!
!  Prepare to return Newton step.
!
     restrt = .true.
     ka = ka + 1
     k = 0
     do i = 1, p
       k = k + i
       j = diag0 + i
       dihdi(k) = w(j)
     end do
     uk = negone
     go to 30

 310  kamin = ka + 3
  if (v(dgnorm) == 0.0D+00) kamin = 0
  if (ka == 0) go to 50

  dst = w(dstsav)
  alphak = abs(v(stppar))
  phi = dst - rad
  t = v(dgnorm)/rad
  uk = t - w(emin)
  if (v(dgnorm) == 0.0D+00) uk = uk + p001 + p001*uk
  if (uk <= 0.0D+00) uk = p001
  if (rad > v(rad0)) go to 320
!
!  Smaller radius 
!
     lk = 0.0D+00
     if (alphak > 0.0D+00) lk = w(lk0)
     lk = max (lk, t - w(emax))
     if (v(dst0) > 0.0D+00) lk = max (lk, (v(dst0)-rad)*w(phipin))
     go to 250
!
!  Bigger radius.
!
 320  if (alphak > 0.0D+00) uk = min (uk, w(uk0))
  lk = max (zero, -v(dst0), t - w(emax))
  if (v(dst0) > 0.0D+00) lk = max (lk, (v(dst0)-rad)*w(phipin))
  go to 250
!
!  Decide whether to check for special case... in practice (from
!  the standpoint of the calling optimization code) it seems best
!  not to check until a few iterations have failed -- hence the
!  test on kamin below.
!
 330  delta = alphak + min (zero, v(dst0))
  twopsi = alphak*dst*dst + gtsta
  if (ka >= kamin) go to 340
!
!  if the test in ref. 2 is satisfied, fall through to handle
!  the special case (as soon as the more-sorensen test detects
!  it).
!
  if (delta >= psifac*twopsi) go to 370
!
!  check for the special case of  h + alpha*d**2  (nearly)
!  singular.  use one step of inverse power method with start
!  from lsvmin to obtain approximate eigenvector corresponding
!  to smallest eigenvalue of (d**-1)*h*(d**-1).  lsvmin returns
!  x and w with  l*w = x.
!
 340  t = lsvmin(p, l, w(x), w)
!
!  Normalize w.
!
  w(1:p) = t * w(1:p)
!
!  Complete current inv. power iter. -- replace w by (l**-t)*w.
!
  call litvmu(p, w, l, w)
  t2 = one/v2norm(p, w)

  w(1:p) = t2 * w(1:p)

  t = t2 * t
!
!   now w is the desired approximate (unit) eigenvector and
!   t*x = ((d**-1)*h*(d**-1) + alphak*i)*w.
!
  sw = dotprd(p, w(q), w)
  t1 = (rad + dst) * (rad - dst)
  root = sqrt(sw*sw + t1)
  if (sw < 0.0D+00) root = -root
  si = t1 / (sw + root)
!
!  The actual test for the special case:
!
  if ((t2*si)**2 <= eps*(dst**2 + alphak*radsq)) go to 380
!
!  Update upper bound on smallest eigenvalue (when not positive)
!  (as recommended by more and sorensen) and continue...
!
  if (v(dst0) <= 0.0D+00) v(dst0) = min (v(dst0), t2**2 - alphak)
  lk = max (lk, -v(dst0))
!
!  Check whether we can hope to detect the special case in
!  the available arithmetic.  accept step as it is if not.
!
!  If not yet available, obtain machine dependent value dgxfac.
!
 370  if (dgxfac == 0.0D+00) dgxfac = epsfac * rmdcon(3)

  if (delta > dgxfac*w(dggdmx)) go to 250
     go to 270
!
!  Special case detected... negate alphak to indicate special case
!
 380  alphak = -alphak
  v(preduc) = half * twopsi
!
!  Accept current step if adding si*w would lead to a
!  further relative reduction in psi of less than v(epslon)/3.
!
  t1 = 0.0D+00
  t = si*(alphak*sw - half*si*(alphak + t*dotprd(p,w(x),w)))
  if (t < eps*twopsi/six) go to 390
     v(preduc) = v(preduc) + t
     dst = rad
     t1 = -si
 390  continue

  do i = 1, p
     j = q0 + i
     w(j) = t1*w(i) - w(j)
     step(i) = w(j) / d(i)
  end do

  v(gtstep) = dotprd(p, dig, w(q))
!
!  Save values for use in a possible restart 
!
 410  v(dstnrm) = dst
  v(stppar) = alphak
  w(lk0) = lk
  w(uk0) = uk
  v(rad0) = rad
  w(dstsav) = dst
!
!  Restore diagonal of DIHDI.
!
  j = 0
  do i = 1, p
    j = j + i
    k = diag0 + i
    dihdi(j) = w(k)
  end do

  return
end
subroutine humit ( d, fx, g, h, iv, lh, liv, lv, n, v, x )

!*******************************************************************************
!
!! HUMIT carries out unconstrained minimization iterations for HUMSL.
!
!  Discussion:
!
!    The Hessian matrix is provided by the caller.
!
!    parameters iv, n, v, and x are the same as the corresponding
!    ones to humsl (which see), except that v can be shorter (since
!    the part of v that humsl uses for storing g and h is not needed).
!    moreover, compared with humsl, iv(1) may have the two additional
!    output values 1 and 2, which are explained below, as is the use
!    of iv(toobig) and iv(nfgcal).  the value iv(g), which is an
!    output value from humsl, is not referenced by humit or the
!    subroutines it calls.
!
!  Parameters:
!
! d.... scale vector.
! fx... function value.
! g.... gradient vector.
! h.... lower triangle of the hessian, stored rowwise.
! iv... integer value array.
! lh... length of h = p*(p+1)/2.
! liv.. length of iv (at least 60).
! lv... length of v (at least 78 + n*(n+21)/2).
! n.... number of variables (components in x and g).
! v.... floating-point value array.
! x.... parameter vector.
!
! iv(1) = 1 means the caller should set fx to f(x), the function value
!             at x, and call humit again, having changed none of the
!             other parameters.  an exception occurs if f(x) cannot be
!             computed (e.g. if overflow would occur), which may happen
!             because of an oversized step.  in this case the caller
!             should set iv(toobig) = iv(2) to 1, which will cause
!             humit to ignore fx and try a smaller step.  the para-
!             meter nf that humsl passes to calcf (for possible use by
!             calcgh) is a copy of iv(nfcall) = iv(6).
! iv(1) = 2 means the caller should set g to g(x), the gradient of f at
!             x, and h to the lower triangle of h(x), the hessian of f
!             at x, and call humit again, having changed none of the
!             other parameters except perhaps the scale vector d.
!                  the parameter nf that humsl passes to calcg is
!             iv(nfgcal) = iv(7).  if g(x) and h(x) cannot be evaluated,
!             then the caller may set iv(nfgcal) to 0, in which case
!             humit will return with iv(1) = 65.
!                  note -- humit overwrites h with the lower triangle
!             of  diag(d)**-1 * h(x) * diag(d)**-1.
!
  integer lh
  integer liv
  integer lv
  integer n

  integer iv(liv)
  real ( kind = 8 ) d(n), fx, g(n), h(lh), v(lv), x(n)
  integer dg1, i, j, k, l, lstgst, nn1o2, step1
  integer temp1, w1, x01
  real ( kind = 8 ) t
  real ( kind = 8 ) one, onep2
  logical stopx
  real ( kind = 8 ) dotprd, reldst, v2norm
  integer cnvcod, dg, dgnorm, dinit, dstnrm, dtinit, dtol
  integer dtype, d0init, f, f0, fdif, gtstep, incfac, irc, kagqt
  integer lmat, lmax0, lmaxs, mode, model, mxfcal, mxiter, nextv
  integer nfcall, nfgcal, ngcall, niter, preduc, radfac, radinc
  integer radius, rad0, reldx, restor, step, stglim, stlstg, stppar
  integer toobig, tuner4, tuner5, vneed, w, xirc, x0

  parameter (cnvcod=55, dg=37, dtol=59, dtype=16, irc=29, kagqt=33 )
  parameter ( lmat=42, mode=35, model=5, mxfcal=17, mxiter=18 )
  parameter ( nextv=47, nfcall=6, nfgcal=7, ngcall=30, niter=31 )
  parameter ( radinc=8, restor=9, step=40, stglim=11, stlstg=41 )
  parameter ( toobig=2, vneed=4, w=34, xirc=13, x0=43)
  parameter (dgnorm=1, dinit=38, dstnrm=2, dtinit=39, d0init=40 )
  parameter ( f=10, f0=13, fdif=11, gtstep=4, incfac=23, lmax0=35 )
  parameter ( lmaxs=36, preduc=7, radfac=16, radius=8, rad0=9 )
  parameter ( reldx=17, stppar=5, tuner4=29, tuner5=30)

  parameter (one=1.d+0, onep2=1.2d+0 )

  i = iv(1)
  if (i == 1) go to 30
  if (i == 2) go to 40
!
!   check validity of iv and v input values 
!
  if (iv(1) == 0) call deflt(2, iv, liv, lv, v)
  if (iv(1) == 12 .or. iv(1) == 13) then
   iv(vneed) = iv(vneed) + n*(n+21)/2 + 7
  end if
  call parck(2, d, iv, liv, lv, n, v)
  i = iv(1) - 2
  if (i > 12) then
    return
  end if
  nn1o2 = n * (n + 1) / 2
  if (lh >= nn1o2) go to (210,210,210,210,210,210,160,120,160,10,10,20), i
     iv(1) = 66
     go to 350
!
!   storage allocation 
!
 10   iv(dtol) = iv(lmat) + nn1o2
  iv(x0) = iv(dtol) + 2*n
  iv(step) = iv(x0) + n
  iv(stlstg) = iv(step) + n
  iv(dg) = iv(stlstg) + n
  iv(w) = iv(dg) + n
  iv(nextv) = iv(w) + 4*n + 7
  if (iv(1) /= 13) go to 20
     iv(1) = 14
     return
!
!   initialization 
!
 20   iv(niter) = 0
  iv(nfcall) = 1
  iv(ngcall) = 1
  iv(nfgcal) = 1
  iv(mode) = -1
  iv(model) = 1
  iv(stglim) = 1
  iv(toobig) = 0
  iv(cnvcod) = 0
  iv(radinc) = 0
  v(rad0) = 0.0D+00
  v(stppar) = 0.0D+00
  if (v(dinit) >= 0.0D+00) call vscopy(n, d, v(dinit))
  k = iv(dtol)
  if (v(dtinit) > 0.0D+00) call vscopy(n, v(k), v(dtinit))
  k = k + n
  if (v(d0init) > 0.0D+00) call vscopy(n, v(k), v(d0init))
  iv(1) = 1
  return

 30   v(f) = fx
  if (iv(mode) >= 0) go to 210
  iv(1) = 2
  if (iv(toobig) == 0) then
    return
  end if
     iv(1) = 63
     go to 350
!
!  Make sure gradient could be computed 
!
 40   if (iv(nfgcal) /= 0) go to 50
     iv(1) = 65
     go to 350
!
!  Update the scale vector d 
!
 50   dg1 = iv(dg)
  if (iv(dtype) <= 0) go to 70
  k = dg1
  j = 0
  do i = 1, n
     j = j + i
     v(k) = h(j)
     k = k + 1
  end do

  call dupdu(d, v(dg1), iv, liv, lv, n, v)
!
!  Compute scaled gradient and its norm 
!
 70   dg1 = iv(dg)
  k = dg1
  do i = 1, n
     v(k) = g(i) / d(i)
     k = k + 1
  end do

  v(dgnorm) = v2norm(n, v(dg1))
!
!  Compute scaled hessian 
!
  k = 1
  do i = 1, n
     t = one / d(i)
     do j = 1, i
          h(k) = t * h(k) / d(j)
          k = k + 1
     end do
  end do

  if (iv(cnvcod) /= 0) go to 340
  if (iv(mode) == 0) go to 300
!
!  Allow first step to have scaled 2-norm at most v(lmax0) 
!
  v(radius) = v(lmax0)

  iv(mode) = 0
!
!  Main loop  
!
!  print iteration summary, check iteration limit 
!
 110  call itsum(d, g, iv, liv, lv, n, v, x)
 120  k = iv(niter)
  if (k < iv(mxiter)) go to 130
     iv(1) = 10
     go to 350

 130  iv(niter) = k + 1
!
!  initialize for start of next iteration 
!
  dg1 = iv(dg)
  x01 = iv(x0)
  v(f0) = v(f)
  iv(irc) = 4
  iv(kagqt) = -1
!
!  Copy x to x0 
!
  call vcopy ( n, v(x01), x )
!
!  Update radius 
!
  if (k == 0) go to 150
  step1 = iv(step)
  k = step1
  do i = 1, n
     v(k) = d(i) * v(k)
     k = k + 1
  end do
  v(radius) = v(radfac) * v2norm(n, v(step1))
!
!  Check STOPX and function evaluation limit.
!
 150  if (.not. stopx ( ) ) go to 170
     iv(1) = 11
     go to 180
!
!  Come here when restarting after func. eval. limit or STOPX.
!
 160  if (v(f) >= v(f0)) go to 170
     v(radfac) = one
     k = iv(niter)
     go to 130

 170  if (iv(nfcall) < iv(mxfcal)) go to 190
     iv(1) = 9
 180     if (v(f) >= v(f0)) go to 350
!
!  In case of STOPX or function evaluation limit with
!  improved v(f), evaluate the gradient at x.
!
          iv(cnvcod) = iv(1)
          go to 290
!
!  Compute candidate step  
!
 190  step1 = iv(step)
  dg1 = iv(dg)
  l = iv(lmat)
  w1 = iv(w)
  call gqtst(d, v(dg1), h, iv(kagqt), v(l), n, v(step1), v, v(w1))
  if (iv(irc) == 6) go to 210
!
!  Check whether evaluating f(x0 + step) looks worthwhile 
!
  if (v(dstnrm) <= 0.0D+00) go to 210
  if (iv(irc) /= 5) go to 200
  if (v(radfac) <= one) go to 200
  if (v(preduc) <= onep2 * v(fdif)) go to 210
!
!  Compute f(x0 + step) 
!
 200  x01 = iv(x0)
  step1 = iv(step)
  call vaxpy(n, x, one, v(step1), v(x01))
  iv(nfcall) = iv(nfcall) + 1
  iv(1) = 1
  iv(toobig) = 0
  return
!
!  Assess candidate step.
!
 210  x01 = iv(x0)
  v(reldx) = reldst(n, d, x, v(x01))
  call assst(iv, liv, lv, v)
  step1 = iv(step)
  lstgst = iv(stlstg)
  if (iv(restor) == 1) call vcopy(n, x, v(x01))
  if (iv(restor) == 2) call vcopy(n, v(lstgst), v(step1))
  if (iv(restor) /= 3) go to 220
     call vcopy(n, v(step1), v(lstgst))
     call vaxpy(n, x, one, v(step1), v(x01))
     v(reldx) = reldst(n, d, x, v(x01))

 220  k = iv(irc)
  go to (230,260,260,260,230,240,250,250,250,250,250,250,330,300), k
!
!  Recompute step with new radius 
!
 230     v(radius) = v(radfac) * v(dstnrm)
     go to 150
!
!  Compute step of length v(lmaxs) for singular convergence test.
!
 240  v(radius) = v(lmaxs)
  go to 190
!
!  Convergence or false convergence 
!
 250  iv(cnvcod) = k - 4
  if (v(f) >= v(f0)) go to 340
     if (iv(xirc) == 14) go to 340
          iv(xirc) = 14
!
!  Process acceptable step.
!
 260  if (iv(irc) /= 3) go to 290
     temp1 = lstgst
!
!  prepare for gradient tests 
!  set  temp1 = hessian * step + g(x0)
!  = diag(d) * (h * step + g(x0))
!
!  Use x0 vector as temporary.
!
     k = x01

     do i = 1, n
       v(k) = d(i) * v(step1)
       k = k + 1
       step1 = step1 + 1
     end do

     call slvmul(n, v(temp1), h, v(x01))

     do i = 1, n
       v(temp1) = d(i) * v(temp1) + g(i)
       temp1 = temp1 + 1
     end do
!
!  Compute gradient and hessian.
!
 290  iv(ngcall) = iv(ngcall) + 1
  iv(1) = 2
  return

 300  iv(1) = 2
  if (iv(irc) /= 3) go to 110
!
!  Set v(radfac) by gradient tests.
!
  temp1 = iv(stlstg)
  step1 = iv(step)
!
!  Set  temp1 = diag(d)**-1 * (hessian*step + (g(x0)-g(x))) 
!
  k = temp1
  do i = 1, n
     v(k) = ( v(k) - g(i) ) / d(i)
     k = k + 1
  end do
!
!  Do gradient tests, 
!
  if (v2norm(n, v(temp1)) <= v(dgnorm) * v(tuner4)) go to 320
       if (dotprd(n, g, v(step1)) >= v(gtstep) * v(tuner5))  go to 110
 320            v(radfac) = v(incfac)
            go to 110
!
!  misc. details   
!
!  bad parameters to assess 
!
 330  iv(1) = 64
  go to 350
!
!  Print summary of final iteration and other requested items 
!
 340  iv(1) = iv(cnvcod)
  iv(cnvcod) = 0
 350  call itsum(d, g, iv, liv, lv, n, v, x)

  return
end
subroutine humsl ( n, d, x, calcf, calcgh, iv, liv, lv, v, uiparm, &
  urparm, ufparm )

!*******************************************************************************
!
!! HUMSL minimizes a general unconstrained objective function.
!
!  Discussion:
!
!    The gradient and Hessian are provided by the caller.
!
!    this routine is like sumsl, except that the subroutine para-
!    meter calcg of sumsl (which computes the gradient of the objec-
!    tive function) is replaced by the subroutine parameter calcgh,
!    which computes both the gradient and (lower triangle of the)
!    hessian of the objective function.  
!
!  Reference:
!
!    David Gay,
!    Computing Optimal Locally Constrained Steps,
!    SIAM Journal on Scientific and Statistical Computing,
!    Volume 2, pages 186-197, 1981.
!
!  Parameters:
!
!    the calling sequence is...
!             call calcgh(n, x, nf, g, h, uiparm, urparm, ufparm)
!     parameters n, x, nf, g, uiparm, urparm, and ufparm are the same
!     as for sumsl, while h is an array of length n*(n+1)/2 in which
!     calcgh must store the lower triangle of the hessian at x.  start-
!     ing at h(1), calcgh must store the hessian entries in the order
!     (1,1), (2,1), (2,2), (3,1), (3,2), (3,3), ...
!        the value printed (by itsum) in the column labelled stppar
!     is the levenberg-marquardt used in computing the current step.
!     zero means a full newton step.  if the special case described in
!     ref. 1 is detected, then stppar is negated.  the value printed
!     in the column labelled npreldf is zero if the current hessian
!     is not positive definite.
!        it sometimes proves worthwhile to let d be determined from the
!     diagonal of the hessian matrix by setting iv(dtype) = 1 and
!     v(dinit) = 0.  the following iv and v components are relevant...
!
! iv(dtol)  iv(59) gives the starting subscript in v of the dtol
!             array used when d is updated.  (iv(dtol) can be
!             initialized by calling humsl with iv(1) = 13.)
! iv(dtype).... iv(16) tells how the scale vector d should be chosen.
!             iv(dtype) <= 0 means that d should not be updated, and
!             iv(dtype) >= 1 means that d should be updated as
!             described below with v(dfac).  default = 0.
! v(dfac)  v(41) and the dtol and d0 arrays (see v(dtinit) and
!             v(d0init)) are used in updating the scale vector d when
!             iv(dtype) > 0.  (d is initialized according to
!             v(dinit), described in sumsl.)  let
!                  d1(i) = max(sqrt(abs(h(i,i))), v(dfac)*d(i)),
!             where h(i,i) is the i-th diagonal element of the current
!             hessian.  if iv(dtype) = 1, then d(i) is set to d1(i)
!             unless d1(i) < dtol(i), in which case d(i) is set to
!                  max(d0(i), dtol(i)).
!             if iv(dtype) >= 2, then d is updated during the first
!             iteration as for iv(dtype) = 1 (after any initialization
!             due to v(dinit)) and is left unchanged thereafter.
!             default = 0.6.
! v(dtinit)... v(39), if positive, is the value to which all components
!             of the dtol array (see v(dfac)) are initialized.  if
!             v(dtinit) = 0, then it is assumed that the caller has
!             stored dtol in v starting at v(iv(dtol)).
!             default = 10**-6.
! v(d0init)... v(40), if positive, is the value to which all components
!             of the d0 vector (see v(dfac)) are initialized.  if
!             v(dfac) = 0, then it is assumed that the caller has
!             stored d0 in v starting at v(iv(dtol)+n).  default = 1.0.
!
  integer liv
  integer lv
  integer n

  integer iv(liv)
  integer uiparm(*)
  real ( kind = 8 ) d(n), x(n), v(lv), urparm(*)
!     dimension v(78 + n*(n+12)), uiparm(*), urparm(*)
  external ufparm
  integer g1, h1, iv1, lh, nf
  real ( kind = 8 ) f
  integer g, h, nextv, nfcall, nfgcal, toobig, vneed

  parameter (nextv=47, nfcall=6, nfgcal=7, g=28, h=56, toobig=2, vneed=4)

  lh = n * (n + 1) / 2

  if ( iv(1) == 0 ) then
    call deflt ( 2, iv, liv, lv, v )
  end if

  if (iv(1) == 12 .or. iv(1) == 13) then
    iv(vneed) = iv(vneed) + n*(n+3)/2
  end if

  iv1 = iv(1)
  if (iv1 == 14) go to 10
  if (iv1 > 2 .and. iv1 < 12) go to 10
  g1 = 1
  h1 = 1
  if (iv1 == 12) iv(1) = 13
  go to 20

 10   g1 = iv(g)
  h1 = iv(h)

 20   call humit(d, f, v(g1), v(h1), iv, lh, liv, lv, n, v, x)

  if (iv(1) - 2) 30, 40, 50

 30   nf = iv(nfcall)
  call calcf(n, x, nf, f, uiparm, urparm, ufparm)
  if (nf <= 0) iv(toobig) = 1
  go to 20

 40   call calcgh(n, x, iv(nfgcal), v(g1), v(h1), uiparm, urparm, ufparm)
  go to 20

 50   if (iv(1) /= 14) then
        return
      end if
!
!  storage allocation
!
  iv(g) = iv(nextv)
  iv(h) = iv(g) + n
  iv(nextv) = iv(h) + n*(n+1)/2

  if ( iv1 /= 13 ) then
    go to 10
  end if

  return
end
subroutine itsum ( d, g, iv, liv, lv, p, v, x )

!*******************************************************************************
!
!! ITSUM prints an iteration summary.
!
  integer liv
  integer lv
  integer p

  integer iv(liv)
  real ( kind = 8 ) d(p), g(p), v(lv), x(p)
  integer alg, i, iv1, m, nf, ng, ol, pu
  character*4 model1(6), model2(6)
  real ( kind = 8 ) nreldf, oldf, preldf, reldf
  integer algsav, dstnrm, f, fdif, f0, needhd, nfcall, nfcov, ngcov
  integer ngcall, niter, nreduc, outlev, preduc, prntit, prunit
  integer reldx, solprt, statpr, stppar, sused, x0prt

  parameter (algsav=51, needhd=36, nfcall=6, nfcov=52, ngcall=30 )
  parameter ( ngcov=53, niter=31, outlev=19, prntit=39, prunit=21 )
  parameter ( solprt=22, statpr=23, sused=64, x0prt=24)
  parameter (dstnrm=2, f=10, f0=13, fdif=11, nreduc=6, preduc=7, reldx=17 )
  parameter ( stppar=5)

  data model1/'    ','    ','    ','    ','  g ','  s '/
  data model2/' g  ',' s  ','g-s ','s-g ','-s-g','-g-s'/

  pu = iv(prunit)

  if ( pu == 0 ) then
    return
  end if

  iv1 = iv(1)
  if (iv1 > 62) iv1 = iv1 - 51
  ol = iv(outlev)
  alg = iv(algsav)
  if (iv1 < 2 .or. iv1 > 15) go to 370
  if (iv1 >= 12) go to 120
  if (iv1 == 2 .and. iv(niter) == 0) go to 390
  if (ol == 0) go to 120
  if (iv1 >= 10 .and. iv(prntit) == 0) go to 120
  if (iv1 > 2) go to 10
     iv(prntit) = iv(prntit) + 1
     if (iv(prntit) < iabs(ol)) then
       return
     end if
 10   nf = iv(nfcall) - iabs(iv(nfcov))
  iv(prntit) = 0
  reldf = 0.0D+00
  preldf = 0.0D+00
  oldf = max (abs(v(f0)), abs(v(f)))
  if (oldf <= 0.0D+00) go to 20
     reldf = v(fdif) / oldf
     preldf = v(preduc) / oldf
 20   if (ol > 0) go to 60
!
!  print short summary line 
!
     if (iv(needhd) == 1 .and. alg == 1) write(pu,30)
 30   format(/10h   it   nf,6x,1hf,7x,5hreldf,3x,6hpreldf,3x,5hreldx,2x,13hmodel  stppar)
     if (iv(needhd) == 1 .and. alg == 2) write(pu,40)
 40   format(/11h    it   nf,7x,1hf,8x,5hreldf,4x,6hpreldf,4x,5hreldx,3x,6hstppar)
     iv(needhd) = 0
     if (alg == 2) go to 50
     m = iv(sused)
     write(pu,100) iv(niter), nf, v(f), reldf, preldf, v(reldx), &
       model1(m), model2(m), v(stppar)
     go to 120

 50  write(pu,110) iv(niter), nf, v(f), reldf, preldf, v(reldx), v(stppar)
     go to 120
!
!  print long summary line 
!
 60   if (iv(needhd) == 1 .and. alg == 1) write(pu,70)
 70   format(/11h    it   nf,6x,1hf,7x,5hreldf,3x,6hpreldf,3x,5hreldx, &
        2x,13hmodel  stppar,2x,6hd*step,2x,7hnpreldf)
  if (iv(needhd) == 1 .and. alg == 2) write(pu,80)
 80   format(/11h    it   nf,7x,1hf,8x,5hreldf,4x,6hpreldf,4x,5hreldx, &
        3x,6hstppar,3x,6hd*step,3x,7hnpreldf)
  iv(needhd) = 0
  nreldf = 0.0D+00
  if (oldf > 0.0D+00) nreldf = v(nreduc) / oldf
  if (alg == 2) go to 90
  m = iv(sused)
  write(pu,100) iv(niter), nf, v(f), reldf, preldf, v(reldx), &
              model1(m), model2(m), v(stppar), v(dstnrm), nreldf
  go to 120

 90   write(pu,110) iv(niter), nf, v(f), reldf, preldf, &
             v(reldx), v(stppar), v(dstnrm), nreldf
 100  format(i6,i5,d10.3,2d9.2,d8.1,a3,a4,2d8.1,d9.2)
 110  format(i6,i5,d11.3,2d10.2,3d9.1,d10.2)

 120  if (iv(statpr) < 0) go to 430
  go to (999, 999, 130, 150, 170, 190, 210, 230, 250, 270, 290, 310, 330, 350, 520), iv1

 130  write(pu,140)
 140  format(/' x-convergence')
  go to 430

 150  write ( pu, 160 )
 160  format(/'relative function convergence')
  go to 430

 170  write(pu,180)
 180  format(/'x- and relative function convergence')
  go to 430

 190  write(pu,200)
 200  format(/'Absolute function convergence.')
  go to 430

 210  write(pu,220)
 220  format(/'Singular convergence.')
  go to 430

 230  write(pu,240)
 240  format(/'False convergence.')
  go to 430

 250  write(pu,260)
 260  format(/'Function evaluation limit.')
  go to 430

 270  write(pu,280)
 280  format(/'Iteration limit.')
  go to 430

 290  write(pu,300)
 300  format(/'STOPX')
  go to 430

 310  write(pu,320)
 320  format(/'Initial f(x) cannot be computed.')

  go to 390

 330  write(pu,340)
 340  format(/'Bad parameters to assess.')
  go to 999

 350  write(pu,360)
 360  format(/'Gradient could not be computed.')
  if (iv(niter) > 0) go to 480
  go to 390

 370  write(pu,380) iv(1)
 380  format(/'iv(1) =',i5)
  go to 999
!
!   initial call on itsum 
!
 390  if (iv(x0prt) /= 0) write(pu,400) (i, x(i), d(i), i = 1, p)
 400  format(/23h     i     initial x(i),8x,4hd(i)//(1x,i5,d17.6,d14.3))
!     the following are to avoid undefined variables when the
!     function evaluation limit is 1...
!
  v(dstnrm) = 0.0D+00
  v(fdif) = 0.0D+00
  v(nreduc) = 0.0D+00
  v(preduc) = 0.0D+00
  v(reldx)  = 0.0D+00
  if (iv1 >= 12) go to 999
  iv(needhd) = 0
  iv(prntit) = 0
  if (ol == 0) go to 999
  if (ol < 0 .and. alg == 1) write(pu,30)
  if (ol < 0 .and. alg == 2) write(pu,40)
  if (ol > 0 .and. alg == 1) write(pu,70)
  if (ol > 0 .and. alg == 2) write(pu,80)
  if (alg == 1) write(pu,410) v(f)
  if (alg == 2) write(pu,420) v(f)
 410  format(/11h     0    1,d10.3)
!365  format(/11h     0    1,e11.3)
 420  format(/11h     0    1,d11.3)
  go to 999
!
!  Print various information requested on solution 
!
 430  iv(needhd) = 1
  if (iv(statpr) == 0) go to 480
     oldf = max (abs(v(f0)), abs(v(f)))
     preldf = 0.0D+00
     nreldf = 0.0D+00
     if (oldf <= 0.0D+00) go to 440
          preldf = v(preduc) / oldf
          nreldf = v(nreduc) / oldf
 440     nf = iv(nfcall) - iv(nfcov)
     ng = iv(ngcall) - iv(ngcov)
     write(pu,450) v(f), v(reldx), nf, ng, preldf, nreldf
 450  format(/9h function,d17.6,8h   reldx,d17.3/12h func. evals, &
    i8,9x,11hgrad. evals,i8/7h preldf,d16.3,6x,7hnpreldf,d15.3)

     if (iv(nfcov) > 0) write(pu,460) iv(nfcov)
 460     format(/1x,i4,50h extra func. evals for covariance and diagnostics.)
     if (iv(ngcov) > 0) write(pu,470) iv(ngcov)
 470     format(1x,i4,50h extra grad. evals for covariance and diagnostics.)

 480  if (iv(solprt) == 0) go to 999
     iv(needhd) = 1
     write(pu,490)
 490  format(/22h     i      final x(i),8x,4hd(i),10x,4hg(i)/)
     do i = 1, p
          write(pu,510) i, x(i), d(i), g(i)
     end do
 510     format(1x,i5,d16.6,2d14.3)
  go to 999

 520  write(pu,530)
 530  format(/'Inconsistent dimensions.')
 999  continue

  return
end
subroutine litvmu ( n, x, l, y )

!*******************************************************************************
!
!! LITVMU solves L' * x = y.
!
!  Discussion:
!
!    L is an  n x n  lower triangular
!    matrix stored compactly by rows.  x and y may occupy the same
!    storage. 
!
  integer n

  real ( kind = 8 ) l(*)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)
  integer i, ii, ij, i0, j
  real ( kind = 8 ) xi

  x(1:n) = y(1:n)

  i0 = n*(n+1)/2

  do ii = 1, n
    i = n+1 - ii
    xi = x(i)/l(i0)
    x(i) = xi
    if ( i <= 1 ) then
      exit
    end if
    i0 = i0 - i
    if ( xi /= 0.0D+00 ) then
      do j = 1, i-1
        ij = i0 + j
        x(j) = x(j) - xi*l(ij)
      end do
    end if
  end do

  return
end
subroutine livmul ( n, x, l, y )

!*******************************************************************************
!
!! LIVMUL solves L * x = y.
!
!  Discussion:
!
!    L is an  n x n  lower triangular
!    matrix stored compactly by rows.  x and y may occupy the same
!    storage. 
!
  integer n

  real ( kind = 8 ) x(n), l(*), y(n)
  external dotprd
  real ( kind = 8 ) dotprd
  integer i, j, k
  real ( kind = 8 ) t

  do k = 1, n
    if (y(k) /= 0.0D+00 ) go to 20
    x(k) = 0.0D+00
  end do

  return

20 continue

  j = k*(k+1)/2
  x(k) = y(k) / l(j)

  if (k >= n) then
    return
  end if

  k = k + 1

  do i = k, n
     t = dotprd(i-1, l(j+1), x)
     j = j + i
     x(i) = (y(i) - t)/l(j)
  end do

  return
end
subroutine lsqrt ( n1, n, l, a, irc )

!*******************************************************************************
!
!! LSQRT computes rows N1 through N of the Cholesky factor L.
!
!  Discussion:
!
!   The Cholesky factor L satisfies a = l*(l**t),  where  l  and the 
!   lower triangle of  a  are both stored compactly by rows (and may occupy 
!   the same storage).
!
!   irc = 0 means all went well.  irc = j means the leading
!   principal  j x j  submatrix of  a  is not positive definite --
!   and  l(j*(j+1)/2)  contains the (nonpos.) reduced j-th diagonal.
!
!  Parameters:
!
  integer n1
  integer n, irc
  real ( kind = 8 ) l(*), a(*)
!     dimension l(n*(n+1)/2), a(n*(n+1)/2)
!
  integer i, ij, ik, im1, i0, j, jk, j0, k
  real ( kind = 8 ) t, td

  i0 = n1 * (n1 - 1) / 2

  do i = n1, n
     td = 0.0D+00
     if (i == 1) go to 40
     j0 = 0
     im1 = i - 1
     do j = 1, im1
          t = 0.0D+00
          do k = 1, j-1
            ik = i0 + k
            jk = j0 + k
            t = t + l(ik)*l(jk)
          end do
          ij = i0 + j
          j0 = j0 + j
          t = (a(ij) - t) / l(j0)
          l(ij) = t
          td = td + t*t
       end do
 40    i0 = i0 + i
     t = a(i0) - td
     if (t <= 0.0D+00) go to 60
     l(i0) = sqrt(t)
  end do

  irc = 0
  return

 60   l(i0) = t
  irc = i

  return
end
function lsvmin ( p, l, x, y )

!*******************************************************************************
!
!! LSVMIN estimates the smallest singular value of matrix L.
!
!  Discussion:
!
!    L is a packed lower triangular matrix.
!
!    this function returns a good overestimate of the smallest
!    singular value of the packed lower triangular matrix l.
!
!  Reference: 
!
!    Alan Cline, Cleve Moler, G Stewart, and James Wilkinson,
!    An estimate for the Condition Number of a Matrix, 
!    Report TM-310, 1977,
!    Applied Mathematics Division, 
!    Argonne National Laboratory.
!
!    D C Hoaglin,
!    Theoretical properties of congruential random-number generators,
!    an empirical view,
!    memorandum ns-340, dept. of statistics, harvard univ., 1976.
!
!    D E Knuth,
!    The Art of Computer Programming, 
!    Volume 2: Seminumerical Algorithms, 
!    Addison-wesley, reading, mass., 1969.
!
!    C S Smith,
!    Multiplicative pseudo-random number generators with prime modulus, 
!    Journal of the Association for Computing Machinery,
!    Volume 18, pages 586-593, 1971.
!
!  Parameters:
!
!  p (in)  = the order of l.  l is a  p x p  lower triangular matrix.
!  l (in)  = array holding the elements of  l  in row order, i.e.
!             l(1,1), l(2,1), l(2,2), l(3,1), l(3,2), l(3,3), etc.
!  x (out) if lsvmin returns a positive value, then x is a normalized
!             approximate left singular vector corresponding to the
!             smallest singular value.  this approximation may be very
!             crude.  if lsvmin returns zero, then some components of x
!             are zero and the rest retain their input values.
!  y (out) if lsvmin returns a positive value, then y = (l**-1)*x is an
!             unnormalized approximate right singular vector correspond-
!             ing to the smallest singular value.  this approximation
!             may be crude.  if lsvmin returns zero, then y retains its
!             input value.  the caller may pass the same vector for x
!             and y (nonstandard fortran usage), in which case y over-
!             writes x (for nonzero lsvmin returns).
!
!   algorithm notes 
!
!     the algorithm is based on (1), with the additional provision that
!     lsvmin = 0 is returned if the smallest diagonal element of l
!     (in magnitude) is not more than the unit roundoff times the
!     largest.  the algorithm uses a random number generator proposed
!     in (4), which passes the spectral test with flying colors -- see
!     (2) and (3).
!
  integer p
  real ( kind = 8 ) lsvmin
  real ( kind = 8 ) l(*), x(p), y(p)
!     dimension l(p*(p+1)/2)
!
  integer i, ii, ix, j, ji, jj, jjj, jm1, j0, pm1
  real ( kind = 8 ) b, sminus, splus, t, xminus, xplus
  real ( kind = 8 ) half, one, r9973
  real ( kind = 8 ) dotprd, v2norm

  parameter (half=0.5d+0, one=1.d+0, r9973=9973.d+0 )

  ix = 2
  pm1 = p - 1
!
!  First check whether to return lsvmin = 0 and initialize x 
!
  ii = 0
  j0 = p*pm1/2
  jj = j0 + p
  if (l(jj) == 0.0D+00) go to 110
  ix = mod(3432*ix, 9973)
  b = half*(one + float(ix)/r9973)
  xplus = b / l(jj)
  x(p) = xplus
  if (p <= 1) go to 60
  do i = 1, pm1
     ii = ii + i
     if (l(ii) == 0.0D+00) go to 110
     ji = j0 + i
     x(i) = xplus * l(ji)
  end do
!
!  Solve (l**t)*x = b, where the components of b have randomly
!  chosen magnitudes in (.5,1) with signs chosen to make x large.
!
!     do j = p-1 to 1 by -1...
  do 50 jjj = 1, pm1
     j = p - jjj
!
!  determine x(j) in this iteration. note for i = 1,2,...,j
!  that x(i) holds the current partial sum for row i.
!
     ix = mod(3432*ix, 9973)
     b = half*(one + float(ix)/r9973)
     xplus = (b - x(j))
     xminus = (-b - x(j))
     splus = abs(xplus)
     sminus = abs(xminus)
     jm1 = j - 1
     j0 = j*jm1/2
     jj = j0 + j
     xplus = xplus/l(jj)
     xminus = xminus/l(jj)

     do i = 1, jm1
          ji = j0 + i
          splus = splus + abs(x(i) + l(ji)*xplus)
          sminus = sminus + abs(x(i) + l(ji)*xminus)
     end do

 30      if (sminus > splus) xplus = xminus
     x(j) = xplus
!
!  update partial sums 
!
     if (jm1 > 0) call vaxpy(jm1, x, xplus, l(j0+1), x)
 50      continue
!
!  normalize x 
!
 60   t = one/v2norm(p, x)

  x(1:p) = t * x(1:p)
!
!  solve l*y = x and return lsvmin = 1/twonorm(y) 
!
  do j = 1, p
     jm1 = j - 1
     j0 = j*jm1/2
     jj = j0 + j
     t = 0.0D+00
     if (jm1 > 0) t = dotprd(jm1, l(j0+1), y)
     y(j) = (x(j) - t) / l(jj)
  end do

  lsvmin = one/v2norm(p, y)
  return

 110  lsvmin = 0.0D+00
  return
end
subroutine ltvmul ( n, x, l, y )

!*******************************************************************************
!
!! LTVMUL computes  x = (l**t)*y.
!
!  Discussion:
!
!    L is an  n x n  lower triangular matrix stored compactly by rows.  
!    x and y may occupy the same storage. 
!
  integer n
  real ( kind = 8 ) x(n), l(*), y(n)
!     dimension l(n*(n+1)/2)
  integer i, ij, i0, j
  real ( kind = 8 ) yi

  i0 = 0
  do i = 1, n
    yi = y(i)
    x(i) = 0.0D+00
    do j = 1, i
      ij = i0 + j
      x(j) = x(j) + yi * l(ij)
    end do
    i0 = i0 + i
  end do

  return
end
subroutine lupdat ( beta, gamma, l, lambda, lplus, n, w, z )

!*******************************************************************************
!
!! LUPDAT computes lplus = secant update of L.
!
!  Discussion:
!
!    this routine updates the cholesky factor  l  of a symmetric
!    positive definite matrix to which a secant update is being
!    applied -- it computes a cholesky factor  lplus  of
!    l * (i + z*w**t) * (i + w*z**t) * l**t.  it is assumed that  w
!    and  z  have been chosen so that the updated matrix is strictly
!    positive definite.
!
!    this code uses recurrence 3 of ref. 1 (with d(j) = 1 for all j)
!    to compute  lplus  of the form  l * (i + z*w**t) * q,  where  q
!    is an orthogonal matrix that makes the result lower triangular.
!    lplus may have some negative diagonal elements.
!
!  Reference: 
!
!    D Goldfarb, 
!    Factorized Variable Metric Methods for Unconstrained Optimization, 
!    Mathematics of Computation,
!    Volume 30, pages 796-811, 1976.
!
!  Parameters: 
!
!   beta = scratch vector.
!  gamma = scratch vector.
!      l (input) lower triangular matrix, stored rowwise.
! lambda = scratch vector.
!  lplus (output) lower triangular matrix, stored rowwise, which may
!             occupy the same storage as  l.
!      n (input) length of vector parameters and order of matrices.
!      w (input, destroyed on output) right singular vector of rank 1
!             correction to  l.
!      z (input, destroyed on output) left singular vector of rank 1
!             correction to  l.
!
  integer n
  real ( kind = 8 ) beta(n), gamma(n), l(*), lambda(n), lplus(*), w(n), z(n)
!     dimension l(n*(n+1)/2), lplus(n*(n+1)/2)
!
  integer i, ij, j, jj, jp1, k, nm1
  integer np1
  real ( kind = 8 ) a, b, bj, eta, gj, lj, lij, ljj, nu, s, theta, wj, zj
  real ( kind = 8 ) one

  parameter (one=1.d+0 )

  nu = one
  eta = 0.0D+00
  if (n <= 1) go to 30
  nm1 = n - 1
!
!  temporarily store s(j) = sum over k = j+1 to n of w(k)**2 in
!  lambda(j).
!
  s = 0.0D+00
  do i = 1, nm1
     j = n - i
     s = s + w(j+1)**2
     lambda(j) = s
  end do
!
!  compute lambda, gamma, and beta by goldfarb*s recurrence 3.
!
  do 20 j = 1, nm1
     wj = w(j)
     a = nu*z(j) - eta*wj
     theta = one + a*wj
     s = a*lambda(j)
     lj = sqrt(theta**2 + a*s)
     if (theta > 0.0D+00) lj = -lj
     lambda(j) = lj
     b = theta*wj + s
     gamma(j) = b * nu / lj
     beta(j) = (a - b*eta) / lj
     nu = -nu / lj
     eta = -(eta + (a**2)/(theta - lj)) / lj
 20      continue
 30   lambda(n) = one + (nu*z(n) - eta*w(n))*w(n)
!
!  update l, gradually overwriting  w  and  z  with  l*w  and  l*z.
!
  np1 = n + 1
  jj = n * (n + 1) / 2

  do k = 1, n

     j = np1 - k
     lj = lambda(j)
     ljj = l(jj)
     lplus(jj) = lj * ljj
     wj = w(j)
     w(j) = ljj * wj
     zj = z(j)
     z(j) = ljj * zj
     if (k == 1) go to 50
     bj = beta(j)
     gj = gamma(j)
     ij = jj + j
     jp1 = j + 1

     do i = jp1, n
          lij = l(ij)
          lplus(ij) = lj*lij + bj*w(i) + gj*z(i)
          w(i) = w(i) + lij*wj
          z(i) = z(i) + lij*zj
          ij = ij + i
     end do

 50      jj = jj - j

  end do

  return
end
subroutine lvmul ( n, x, l, y )

!*******************************************************************************
!
!! LVMUL computes x = L * y.
!
!  Discussion:
!
!    L  is an  n x n  lower triangular matrix stored compactly by rows.  
!    x and y may occupy the same storage. 
!
  integer n

  real ( kind = 8 ) x(n), l(*), y(n)
!     dimension l(n*(n+1)/2)
  integer i, ii, ij, i0, j, np1
  real ( kind = 8 ) t

  np1 = n + 1
  i0 = n*(n+1)/2

  do ii = 1, n
    i = np1 - ii
    i0 = i0 - i
    t = 0.0D+00
    do j = 1, i
      ij = i0 + j
      t = t + l(ij)*y(j)
    end do
    x(i) = t
  end do

  return
end
subroutine parck ( alg, d, iv, liv, lv, n, v )

!*******************************************************************************
!
!! PARCK checks parameters, prints changed values.
!
!  Discussion:
!
!    alg = 1 for regression, alg = 2 for general unconstrained opt.
!
  integer alg, liv, lv, n
  integer iv(liv)
  real ( kind = 8 ) d(n), v(lv)
  real ( kind = 8 ) rmdcon
  integer max0
  integer i, ii, iv1, j, k, l, m, miv1, miv2, ndfalt, parsv1, pu
  integer ijmp, jlim(2), miniv(2), ndflt(2)
  character*1 varnm(2), sh(2)
  character*4 cngd(3), dflt(3), vn(2,34), which(3)
  real ( kind = 8 ) big, machep, tiny, vk, vm(34), vx(34)
  integer algsav, dinit, dtype, dtype0, epslon, inits, ivneed
  integer lastiv, lastv, lmat, nextiv, nextv, nvdflt, oldn
  integer parprt, parsav, perm, prunit, vneed

  parameter (algsav=51, dinit=38, dtype=16, dtype0=54, epslon=19 )
  parameter ( inits=25, ivneed=3, lastiv=44, lastv=45, lmat=42 )
  parameter ( nextiv=46, nextv=47, nvdflt=50, oldn=38, parprt=20 )
  parameter ( parsav=49, perm=58, prunit=21, vneed=4)
  save big, machep, tiny

  data big/0.d+0/, machep/-1.d+0/, tiny/1.d+0/

     data vn(1,1),vn(2,1)/'epsl','on..'/
     data vn(1,2),vn(2,2)/'phmn','fc..'/
     data vn(1,3),vn(2,3)/'phmx','fc..'/
     data vn(1,4),vn(2,4)/'decf','ac..'/
     data vn(1,5),vn(2,5)/'incf','ac..'/
     data vn(1,6),vn(2,6)/'rdfc','mn..'/
     data vn(1,7),vn(2,7)/'rdfc','mx..'/
     data vn(1,8),vn(2,8)/'tune','r1..'/
     data vn(1,9),vn(2,9)/'tune','r2..'/
     data vn(1,10),vn(2,10)/'tune','r3..'/
     data vn(1,11),vn(2,11)/'tune','r4..'/
     data vn(1,12),vn(2,12)/'tune','r5..'/
     data vn(1,13),vn(2,13)/'afct','ol..'/
     data vn(1,14),vn(2,14)/'rfct','ol..'/
     data vn(1,15),vn(2,15)/'xcto','l...'/
     data vn(1,16),vn(2,16)/'xfto','l...'/
     data vn(1,17),vn(2,17)/'lmax','0...'/
     data vn(1,18),vn(2,18)/'lmax','s...'/
     data vn(1,19),vn(2,19)/'scto','l...'/
     data vn(1,20),vn(2,20)/'dini','t...'/
     data vn(1,21),vn(2,21)/'dtin','it..'/
     data vn(1,22),vn(2,22)/'d0in','it..'/
     data vn(1,23),vn(2,23)/'dfac','....'/
     data vn(1,24),vn(2,24)/'dltf','dc..'/
     data vn(1,25),vn(2,25)/'dltf','dj..'/
     data vn(1,26),vn(2,26)/'delt','a0..'/
     data vn(1,27),vn(2,27)/'fuzz','....'/
     data vn(1,28),vn(2,28)/'rlim','it..'/
     data vn(1,29),vn(2,29)/'cosm','in..'/
     data vn(1,30),vn(2,30)/'hube','rc..'/
     data vn(1,31),vn(2,31)/'rspt','ol..'/
     data vn(1,32),vn(2,32)/'sigm','in..'/
     data vn(1,33),vn(2,33)/'eta0','....'/
     data vn(1,34),vn(2,34)/'bias','....'/

  data vm(1)/1.0d-3/, vm(2)/-0.99d+0/, vm(3)/1.0d-3/, vm(4)/1.0d-2/
  data vm(5)/1.2d+0/, vm(6)/1.d-2/, vm(7)/1.2d+0/, vm(8)/0.d+0/
  data vm(9)/0.d+0/, vm(10)/1.d-3/, vm(11)/-1.d+0/, vm(13)/0.d+0/
  data vm(15)/0.d+0/, vm(16)/0.d+0/, vm(19)/0.d+0/, vm(20)/-10.d+0/
  data vm(21)/0.d+0/, vm(22)/0.d+0/, vm(23)/0.d+0/, vm(27)/1.01d+0/
  data vm(28)/1.d+10/, vm(30)/0.d+0/, vm(31)/0.d+0/, vm(32)/0.d+0/
  data vm(34)/0.d+0/

  data vx(1)/0.9d+0/, vx(2)/-1.d-3/, vx(3)/1.d+1/, vx(4)/0.8d+0/
  data vx(5)/1.d+2/, vx(6)/0.8d+0/, vx(7)/1.d+2/, vx(8)/0.5d+0/
  data vx(9)/0.5d+0/, vx(10)/1.d+0/, vx(11)/1.d+0/, vx(14)/0.1d+0/
  data vx(15)/1.d+0/, vx(16)/1.d+0/, vx(19)/1.d+0/, vx(23)/1.d+0/
  data vx(24)/1.d+0/, vx(25)/1.d+0/, vx(26)/1.d+0/, vx(27)/1.d+10/
  data vx(29)/1.d+0/, vx(31)/1.d+0/, vx(32)/1.d+0/, vx(33)/1.d+0/
  data vx(34)/1.d+0/

  data varnm(1)/'p'/, varnm(2)/'n'/, sh(1)/'s'/, sh(2)/'h'/
  data cngd(1),cngd(2),cngd(3)/'---c','hang','ed v'/
  data dflt(1),dflt(2),dflt(3)/'nond','efau','lt v'/
  data ijmp/33/, jlim(1)/0/, jlim(2)/24/, ndflt(1)/32/, ndflt(2)/25/
  data miniv(1)/80/, miniv(2)/59/

  pu = 0
  if (prunit <= liv) pu = iv(prunit)
  if (alg < 1 .or. alg > 2) go to 340
  if (iv(1) == 0) call deflt(alg, iv, liv, lv, v)
  iv1 = iv(1)
  if (iv1 /= 13 .and. iv1 /= 12) go to 10
  miv1 = miniv(alg)
  if (perm <= liv) miv1 = max0(miv1, iv(perm) - 1)
  if (ivneed <= liv) miv2 = miv1 + max0(iv(ivneed), 0)
  if (lastiv <= liv) iv(lastiv) = miv2
  if (liv < miv1) go to 300
  iv(ivneed) = 0
  iv(lastv) = max0(iv(vneed), 0) + iv(lmat) - 1
  iv(vneed) = 0
  if (liv < miv2) go to 300
  if (lv < iv(lastv)) go to 320
 10   if (alg == iv(algsav)) go to 30
  if (pu /= 0) write(pu,20) alg, iv(algsav)
 20 format(/39h the first parameter to deflt should be,i3,12h rather than,i3)
     iv(1) = 82
     return
 30   if (iv1 < 12 .or. iv1 > 14) go to 60
     if (n >= 1) go to 50
          iv(1) = 81
          if (pu == 0) then
            return
          end if
          write(pu,40) varnm(alg), n
 40           format(/8h /// bad,a1,2h =,i5)
          return
 50      if (iv1 /= 14) iv(nextiv) = iv(perm)
     if (iv1 /= 14) iv(nextv) = iv(lmat)
     if (iv1 == 13) then
       return
     end if
     k = iv(parsav) - epslon
     call vdflt(alg, lv-k, v(k+1))
     iv(dtype0) = 2 - alg
     iv(oldn) = n
     which(1) = dflt(1)
     which(2) = dflt(2)
     which(3) = dflt(3)
     go to 110
 60   if (n == iv(oldn)) go to 80
     iv(1) = 17
     if (pu == 0) then
       return
     end if
     write(pu,70) varnm(alg), iv(oldn), n
 70      format(/5h /// ,1a1,14h changed from ,i5,4h to ,i5)
     return

 80   if (iv1 <= 11 .and. iv1 >= 1) go to 100
     iv(1) = 80
     if (pu /= 0) write(pu,90) iv1
 90      format(/13h ///  iv(1) =,i5,28h should be between 0 and 14.)
     return

 100  which(1) = cngd(1)
  which(2) = cngd(2)
  which(3) = cngd(3)

 110  if (iv1 == 14) iv1 = 12
  if (big > tiny) go to 120
     tiny = rmdcon(1)
     machep = rmdcon(3)
     big = rmdcon(6)
     vm(12) = machep
     vx(12) = big
     vx(13) = big
     vm(14) = machep
     vm(17) = tiny
     vx(17) = big
     vm(18) = tiny
     vx(18) = big
     vx(20) = big
     vx(21) = big
     vx(22) = big
     vm(24) = machep
     vm(25) = machep
     vm(26) = machep
     vx(28) = rmdcon(5)
     vm(29) = machep
     vx(30) = big
     vm(33) = machep
 120  m = 0
  i = 1
  j = jlim(alg)
  k = epslon
  ndfalt = ndflt(alg)

  do l = 1, ndfalt
    vk = v(k)
    if (vk >= vm(i) .and. vk <= vx(i)) go to 140
      m = k
      if (pu /= 0) write(pu,130) vn(1,i), vn(2,i), k, vk,vm(i), vx(i)
 130  format(/6h ///  ,2a4,5h.. v(,i2,3h) =,d11.3,7h should, &
      11h be between,d11.3,4h and,d11.3)
 140  k = k + 1
     i = i + 1
     if (i == j) i = ijmp
  end do

  if (iv(nvdflt) == ndfalt) go to 170
     iv(1) = 51
     if (pu == 0) then
       return
     end if
     write(pu,160) iv(nvdflt), ndfalt
 160     format(/13h iv(nvdflt) =,i5,13h rather than ,i5)
     return
 170  if ((iv(dtype) > 0 .or. v(dinit) > 0.0D+00) .and. iv1 == 12) then
             go to 200
  end if

  do i = 1, n
     if (d(i) > 0.0D+00) go to 190
          m = 18
          if (pu /= 0) write(pu,180) i, d(i)
 180     format(/8h ///  d(,i3,3h) =,d11.3,19h should be positive)
 190     continue
  end do

 200  if (m == 0) go to 210
     iv(1) = m
     return

 210  if (pu == 0 .or. iv(parprt) == 0) then
        return
      end if
  if (iv1 /= 12 .or. iv(inits) == alg-1) go to 230
     m = 1
     write(pu,220) sh(alg), iv(inits)
 220 format(/22h nondefault values..../5h init,a1,14h      iv(25) =,i3)
 230  if (iv(dtype) == iv(dtype0)) go to 250
     if (m == 0) write(pu,260) which
     m = 1
     write(pu,240) iv(dtype)
 240     format(20h dtype      iv(16) =,i3)
 250  i = 1
  j = jlim(alg)
  k = epslon
  l = iv(parsav)
  ndfalt = ndflt(alg)

  do ii = 1, ndfalt
     if (v(k) == v(l)) go to 280
          if (m == 0) write(pu,260) which
 260          format(/1h ,3a4,9halues..../)
          m = 1
          write(pu,270) vn(1,i), vn(2,i), k, v(k)
 270          format(1x,2a4,5h.. v(,i2,3h) =,d15.7)
 280     k = k + 1
     l = l + 1
     i = i + 1
     if (i == j) i = ijmp
  end do

  iv(dtype0) = iv(dtype)
  parsv1 = iv(parsav)
  call vcopy(iv(nvdflt), v(parsv1), v(epslon))
  return

 300  iv(1) = 15
  if (pu == 0) then
    return
  end if
  write(pu,310) liv, miv2
 310  format(/10h /// liv =,i5,17h must be at least,i5)
  if (liv < miv1) then
    return
  end if
  if (lv < iv(lastv)) go to 320
  return

 320  iv(1) = 16
  if (pu == 0) then
    return
  end if
  write(pu,330) lv, iv(lastv)
 330  format(/9h /// lv =,i5,17h must be at least,i5)
  return

 340  iv(1) = 67
  if (pu == 0) then
    return
  end if
  write(pu,350) alg
 350  format(/10h /// alg =,i5,15h must be 1 or 2)

  return
end
function reldst ( p, d, x, x0 )

!*******************************************************************************
!
!! RELDST computes the relative difference between X and X0.
!
  integer p

  real ( kind = 8 ) reldst
  real ( kind = 8 ) d(p), x(p), x0(p)
  integer i
  real ( kind = 8 ) emax, t, xmax

  emax = 0.0D+00
  xmax = 0.0D+00

  do i = 1, p
    t = abs(d(i) * (x(i) - x0(i)))
    if (emax < t) emax = t
    t = d(i) * (abs(x(i)) + abs(x0(i)))
    if (xmax < t) xmax = t
  end do

  reldst = 0.0D+00
  if ( xmax > 0.0D+00 ) reldst = emax / xmax

  return
end
function rmdcon ( k )

!*******************************************************************************
!
!! RMDCON returns machine dependent constants. 
!
!  Discussion:
!
!    Comments below contain data statements for various machines. 
!    To convert to another machine, place a c in column 1 of the  
!    data statement line(s) that correspond to the current machine
!    and remove the c from column 1 of the data statement line(s) 
!    that correspond to the new machine.                          
!
!    the constant returned depends on k...
!
!         k = 1... smallest pos. eta such that -eta exists.
!         k = 2... square root of eta.
!         k = 3... unit roundoff = smallest pos. no. machep such
!                  that 1 + machep > 1 .and. 1 - machep < 1.
!         k = 4... square root of machep.
!         k = 5... square root of big (see k = 6).
!         k = 6... largest machine no. big such that -big exists.
!
  integer k
  real ( kind = 8 ) rmdcon
  real ( kind = 8 ) big, eta, machep
  integer bigi(4), etai(4), machei(4)
  equivalence (big,bigi(1)), (eta,etai(1)), (machep,machei(1))
!
!  ibm 360, ibm 370, or xerox 
!
!     data big/z7fffffffffffffff/, eta/z0010000000000000/,
!    1     machep/z3410000000000000/
!
!  data general 
!
!     data big/0.7237005577d+76/, eta/0.5397605347d-78/,
!    1     machep/2.22044605d-16/
!
!  dec 11 
!
!     data big/1.7d+38/, eta/2.938735878d-39/, machep/2.775557562d-17/
!
!  hp3000 
!
!     data big/1.157920892d+77/, eta/8.636168556d-78/,
!    1     machep/5.551115124d-17/
!
!  honeywell 
!
!     data big/1.69d+38/, eta/5.9d-39/, machep/2.1680435d-19/
!
!  dec10 
!
!     data big/"377777100000000000000000/,
!    1     eta/"002400400000000000000000/,
!    2     machep/"104400000000000000000000/
!
!  burroughs 
!
!     data big/o0777777777777777,o7777777777777777/,
!    1     eta/o1771000000000000,o7770000000000000/,
!    2     machep/o1451000000000000,o0000000000000000/
!
!  control data 
!
!     data big/37767777777777777777b,37167777777777777777b/,
!    1     eta/00014000000000000000b,00000000000000000000b/,
!    2     machep/15614000000000000000b,15010000000000000000b/
!
!  prime 
!
!     data big/1.0d+9786/, eta/1.0d-9860/, machep/1.4210855d-14/
!
!  univac 
!
!     data big/8.988d+307/, eta/1.2d-308/, machep/1.734723476d-18/
!
!  vax 
!
  data big/1.7d+38/, eta/2.939d-39/, machep/1.3877788d-17/
!
!  cray 1 
!
!     data bigi(1)/577767777777777777777b/,
!    1     bigi(2)/000007777777777777776b/,
!    2     etai(1)/200004000000000000000b/,
!    3     etai(2)/000000000000000000000b/,
!    4     machei(1)/377224000000000000000b/,
!    5     machei(2)/000000000000000000000b/
!
!  port library -- requires more than just a data statement...
!
!     external d1mach
!     real ( kind = 8 ) d1mach, zero
!     data big/0.d+0/, eta/0.d+0/, machep/0.d+0/, zero/0.d+0/
!     if (big > 0.0D+00) go to 1
!        big = d1mach(2)
!        eta = d1mach(1)
!        machep = d1mach(4)
!1    continue
!
! end of port
!
!  body -
!
  go to (10, 20, 30, 40, 50, 60), k

 10   rmdcon = eta
  return

 20   rmdcon = sqrt(256.d+0*eta)/16.d+0
  return

 30   rmdcon = machep
  return

 40   rmdcon = sqrt(machep)
  return

 50   rmdcon = sqrt(big/256.d+0)*16.d+0
  return

 60   rmdcon = big

  return
end
subroutine sgrad2 ( alpha, d, eta0, fx, g, irc, n, w, x )

!*******************************************************************************
!
!! SGRAD2 computes finite difference gradient by Stewart's scheme.
!
!  Discussion:
!
!    This subroutine uses an embellished form of the finite difference 
!    scheme proposed by Stewart to approximate the gradient of the 
!    function f(x), whose values are supplied by reverse communication.
!
!  Reference: 
!
!    G W Stewart,
!    A Modification of Davidon's Minimization Method to Accept Difference 
!    Approximations of Derivatives,
!    Journal of the Association for Computing Machinery,
!    Volume 14, pages. 72-83, 1967.
!
!  Parameters:
!
!  alpha in  (approximate) diagonal elements of the hessian of f(x).
!      d in  scale vector such that d(i)*x(i), i = 1,...,n, are in
!             comparable units.
!   eta0 in  estimated bound on relative error in the function value...
!             (true value) = (computed value)*(1+e),   where
!             abs(e) <= eta0.
!     fx i/o on input,  fx  must be the computed value of f(x).  on
!             output with irc = 0, fx has been restored to its original
!             value, the one it had when sgrad2 was last called with
!             irc = 0.
!      g i/o on input with irc = 0, g should contain an approximation
!             to the gradient of f near x, e.g., the gradient at the
!             previous iterate.  when sgrad2 returns with irc = 0, g is
!             the desired finite-difference approximation to the
!             gradient at x.
!    irc i/o input/return code... before the very first call on sgrad2,
!             the caller must set irc to 0.  whenever sgrad2 returns a
!             nonzero value for irc, it has perturbed some component of
!             x... the caller should evaluate f(x) and call sgrad2
!             again with fx = f(x).
!      n in  the number of variables (components of x) on which f
!             depends.
!      x i/o on input with irc = 0, x is the point at which the
!             gradient of f is desired.  on output with irc nonzero, x
!             is the point at which f should be evaluated.  on output
!             with irc = 0, x has been restored to its original value
!             (the one it had when sgrad2 was last called with irc = 0)
!             and g contains the desired gradient approximation.
!      w i/o work vector of length 6 in which sgrad2 saves certain
!             quantities while the caller is evaluating f(x) at a
!             perturbed x.
!
!      application and usage restrictions 
!
!        this routine is intended for use with quasi-newton routines
!     for unconstrained minimization (in which case  alpha  comes from
!     the diagonal of the quasi-newton hessian approximation).
!
!      algorithm notes 
!
!        this code departs from the scheme proposed by stewart (ref. 1)
!     in its guarding against overly large or small step sizes and its
!     handling of special cases (such as zero components of alpha or g).
!
  integer irc, n
  real ( kind = 8 ) alpha(n), d(n), eta0, fx, g(n), w(6), x(n)
  external rmdcon
  real ( kind = 8 ) rmdcon
  integer fh, fx0, hsave, i, xisave
  real ( kind = 8 ) aai, afx, afxeta, agi, alphai, axi, axibar
  real ( kind = 8 ) discon, eta, gi, h, hmin
  real ( kind = 8 ) c2000, four, hmax0, hmin0, h0, machep, one, p002
  real ( kind = 8 ) three, two

  parameter (c2000=2.0d+3, four=4.0d+0, hmax0=0.02d+0, hmin0=5.0d+1 )
  parameter ( one=1.0d+0, p002=0.002d+0, three=3.0d+0 )
  parameter ( two=2.0d+0 )

  parameter (fh=3, fx0=4, hsave=5, xisave=6)
!
  if (irc) 140, 100, 210
!
!      fresh start -- get machine-dependent constants 
!
!     store machep in w(1) and h0 in w(2), where machep is the unit
!     roundoff (the smallest positive number such that
!     1 + machep > 1  and  1 - machep < 1),  and  h0 is the
!     square-root of machep.
!
 100  w(1) = rmdcon(3)
  w(2) = sqrt(w(1))
!
  w(fx0) = fx
!
!      increment  i  and start computing  g(i) 
!
 110  i = iabs(irc) + 1
  if (i > n) go to 300
     irc = i
     afx = abs(w(fx0))
     machep = w(1)
     h0 = w(2)
     hmin = hmin0 * machep
     w(xisave) = x(i)
     axi = abs(x(i))
     axibar = max (axi, one/d(i))
     gi = g(i)
     agi = abs(gi)
     eta = abs(eta0)
     if (afx > 0.0D+00) eta = max (eta, agi*axi*machep/afx)
     alphai = alpha(i)
     if (alphai == 0.0D+00) go to 170
     if (gi == 0.0D+00 .or. fx == 0.0D+00) go to 180
     afxeta = afx*eta
     aai = abs(alphai)
!
!  compute h = stewart's forward-difference step size.
!
     if (gi**2 <= afxeta*aai) go to 120
          h = two*sqrt(afxeta/aai)
          h = h*(one - aai*h/(three*aai*h + four*agi))
          go to 130
 120     h = two*(afxeta*agi/(aai**2))**(one/three)
     h = h*(one - two*agi/(three*aai*h + four*agi))
!
!  ensure that  h  is not insignificantly small 
!
 130     h = max (h, hmin*axibar)
!
!  use forward difference if bound on truncation error is at
!  most 10**-3.
!
     if (aai*h <= p002*agi) go to 160
!
!  compute h = stewart*s step for central difference.
!
     discon = c2000*afxeta
     h = discon/(agi + sqrt(gi**2 + aai*discon))
!
!  ensure that  h  is neither too small nor too big 
!
     h = max (h, hmin*axibar)
     if (h >= hmax0*axibar) h = axibar * h0**(two/three)
!
!  compute central difference 
!
     irc = -i
     go to 200

 140     h = -w(hsave)
     i = iabs(irc)
     if (h > 0.0D+00) go to 150
     w(fh) = fx
     go to 200

 150     g(i) = (w(fh) - fx) / (two * h)
     x(i) = w(xisave)
     go to 110
!
!  Compute forward differences in various cases 
!
 160     if (h >= hmax0*axibar) h = h0 * axibar
     if (alphai*gi < 0.0D+00) h = -h
     go to 200
 170     h = axibar
     go to 200
 180     h = h0 * axibar

 200     x(i) = w(xisave) + h
     w(hsave) = h
     return
!
!  compute actual forward difference 
!
 210     g(irc) = (fx - w(fx0)) / w(hsave)
     x(irc) = w(xisave)
     go to 110
!
!  Restore fx and indicate that g has been computed 
!
 300  fx = w(fx0)
  irc = 0

  return
end
subroutine slvmul ( p, y, s, x )

!*******************************************************************************
!
!! SLVMUL sets y = S * x.
!
!  Discussion:
!
!    s = p x p symmetric matrix. 
!    lower triangle of  s  stored rowwise.        
!
  integer p
  real ( kind = 8 ) s(*), x(p), y(p)
!     dimension s(p*(p+1)/2)
  integer i, im1, j, k
  real ( kind = 8 ) xi
  real ( kind = 8 ) dotprd

  j = 1

  do i = 1, p
    y(i) = dotprd(i, s(j), x)
    j = j + i
  end do

  if (p <= 1) then
    return
  end if

  j = 1

  do i = 2, p

    xi = x(i)
    im1 = i - 1
    j = j + 1

    do k = 1, im1
      y(k) = y(k) + s(j)*xi
      j = j + 1
    end do

  end do

  return
end
subroutine smsno ( n, d, x, calcf, iv, liv, lv, v, uiparm, urparm, ufparm )

!*******************************************************************************
!
!! SMSNO minimizes a general unconstrained objective function.
!
!  Discussion:
!
!    The routine uses finite-difference gradients and secant hessian
!    approximations.
!
!    This routine interacts with SNOIT in an attempt
!    to find an n-vector  x*  that minimizes the (unconstrained)
!    objective function computed by  calcf.  (often the  x*  found is
!    a local minimizer rather than a global one.)
!
!  Reference: 
!
!    G W Stewart,
!    A Modification of Davidon's Minimization Method to Accept Difference 
!      Approximations of Derivatives,
!    Journal of the Association for Computing Machinery,
!    Volume 14, pages 72-83, 1967.
!
!  Parameters:
!
!        the parameters for smsno are the same as those for sumsl
!     (which see), except that calcg is omitted.  instead of calling
!     calcg to obtain the gradient of the objective function at x,
!     smsno calls sgrad2, which computes an approximation to the
!     gradient by finite (forward and central) differences using the
!     method of ref. 1.  the following input component is of interest
!     in this regard (and is not described in sumsl).
!
! v(eta0)  v(42) is an estimated bound on the relative error in the
!             objective function value computed by calcf...
!                  (true value) = (computed value) * (1 + e),
!             where abs(e) <= v(eta0).  default = machep * 10**3,
!             where machep is the unit roundoff.
!
!        the output values iv(nfcall) and iv(ngcall) have different
!     meanings for smsno than for sumsl...
!
! iv(nfcall)... iv(6) is the number of calls so far made on calcf (i.e.,
!             function evaluations) excluding those made only for
!             computing gradients.  the input value iv(mxfcal) is a
!             limit on iv(nfcall).
! iv(ngcall)... iv(30) is the number of function evaluations made only
!             for computing gradients.  the total number of function
!             evaluations is thus  iv(nfcall) + iv(ngcall).
!
  integer n, liv, lv
  integer iv(liv), uiparm(*)
  real ( kind = 8 ) d(n), x(n), v(lv), urparm(*)
!     dimension v(77 + n*(n+17)/2), uiparm(*), urparm(*)
  external calcf, ufparm
  integer nf
  real ( kind = 8 ) fx
  integer nfcall, toobig
  parameter (nfcall=6, toobig=2)

  do

    call snoit ( d, fx, iv, liv, lv, n, v, x )

    if ( iv(1) > 2 ) then
      exit
    end if

    nf = iv(nfcall)

    call calcf ( n, x, nf, fx, uiparm, urparm, ufparm )

    if ( nf <= 0 ) then
      iv(toobig) = 1
    end if

  end do

  return
end
subroutine snoit ( d, fx, iv, liv, lv, n, v, x )

!*******************************************************************************
!
!! SNOIT is the iteration driver for SMSNO.
!
!  Discussion:
!
!    This routine minimizes a general unconstrained objective function using
!    finite-difference gradients and secant hessian approximations.
!
!    This routine interacts with subroutine  sumit  in an attempt
!    to find an n-vector  x*  that minimizes the (unconstrained)
!    objective function  fx = f(x)  computed by the caller.  (often
!    the  x*  found is a local minimizer rather than a global one.)
!
!  Reference: 
!
!    G W Stewart,
!    A Modification of Davidon's Minimization Method to Accept Difference 
!    Approximations of Derivatives,
!    Journal of the Association for Computing Machinery,
!    Volume 14, pages. 72-83, 1967.
!
!  Parameters:
!
!        the parameters for snoit are the same as those for sumsl
!     (which see), except that calcf, calcg, uiparm, urparm, and ufparm
!     are omitted, and a parameter  fx  for the objective function
!     value at x is added.  instead of calling calcg to obtain the
!     gradient of the objective function at x, snoit calls sgrad2,
!     which computes an approximation to the gradient by finite
!     (forward and central) differences using the method of ref. 1.
!     the following input component is of interest in this regard
!     (and is not described in sumsl).
!
! v(eta0)  v(42) is an estimated bound on the relative error in the
!             objective function value computed by calcf...
!                  (true value) = (computed value) * (1 + e),
!             where abs(e) <= v(eta0).  default = machep * 10**3,
!             where machep is the unit roundoff.
!
!        the output values iv(nfcall) and iv(ngcall) have different
!     meanings for smsno than for sumsl...
!
! iv(nfcall)... iv(6) is the number of calls so far made on calcf (i.e.,
!             function evaluations) excluding those made only for
!             computing gradients.  the input value iv(mxfcal) is a
!             limit on iv(nfcall).
! iv(ngcall)... iv(30) is the number of function evaluations made only
!             for computing gradients.  the total number of function
!             evaluations is thus  iv(nfcall) + iv(ngcall).
!
  integer liv, lv, n
  integer iv(liv)
  real ( kind = 8 ) d(n), fx, x(n), v(lv)
!     dimension v(77 + n*(n+17)/2)
!
  external deflt, dotprd, sgrad2, sumit, vscopy
  real ( kind = 8 ) dotprd
  integer alpha, g1, i, iv1, j, k, w
  real ( kind = 8 ) zero

  integer eta0, f, g, lmat, nextv, nfgcal, ngcall
  integer niter, sgirc, toobig, vneed

  parameter ( eta0=42, f=10, g=28, lmat=42, nextv=47 )
  parameter ( nfgcal=7, ngcall=30, niter=31, sgirc=57 )
  parameter ( toobig=2, vneed=4)

  parameter ( zero=0.d+0)

  iv1 = iv(1)
  if (iv1 == 1) go to 10
  if (iv1 == 2) go to 50
  if (iv(1) == 0) call deflt(2, iv, liv, lv, v)
  iv1 = iv(1)
  if (iv1 == 12 .or. iv1 == 13) iv(vneed) = iv(vneed) + 2*n + 6
  if (iv1 == 14) go to 10
  if (iv1 > 2 .and. iv1 < 12) go to 10
  g1 = 1
  if (iv1 == 12) iv(1) = 13
  go to 20

 10   g1 = iv(g)

 20   call sumit(d, fx, v(g1), iv, liv, lv, n, v, x)
  if (iv(1) - 2) 999, 30, 70
!
!  Compute gradient 
!
 30   if (iv(niter) == 0) call vscopy(n, v(g1), zero)
  j = iv(lmat)
  k = g1 - n

  do i = 1, n
    v(k) = dotprd(i, v(j), v(j))
    k = k + 1
    j = j + i
  end do
!
!  Undo increment of iv(ngcall) done by sumit 
!
  iv(ngcall) = iv(ngcall) - 1
!
!  Store return code from sgrad2 in iv(sgirc) 
!
  iv(sgirc) = 0
!
!  x may have been restored, so copy back fx...
!
  fx = v(f)
  go to 60
!
!  gradient loop 
!
 50   if (iv(toobig) == 0) go to 60
  iv(nfgcal) = 0
  go to 10

 60   g1 = iv(g)
  alpha = g1 - n
  w = alpha - 6
  call sgrad2(v(alpha), d, v(eta0), fx, v(g1), iv(sgirc), n, v(w),x)
  if (iv(sgirc) == 0) go to 10
     iv(ngcall) = iv(ngcall) + 1
     return

 70   if (iv(1) /= 14) then
        return
      end if
!
!  Storage allocation 
!
  iv(g) = iv(nextv) + n + 6
  iv(nextv) = iv(g) + n
  if (iv1 /= 13) go to 10

 999  continue

  return
end
function stopx ( )

!*******************************************************************************
!
!! STOPX checks to see if the BREAK key has been pressed.
!
!  Discussion:
!
!     this function may serve as the stopx (asynchronous interruption)
!     function for the nl2sol (nonlinear least-squares) package at
!     those installations which do not wish to implement a
!     dynamic stopx.
!
!     at installations where the nl2sol system is used
!     interactively, this dummy stopx should be replaced by a
!     function that returns .true. if and only if the interrupt
!     (break) key has been pressed since the last call on stopx.
!
  logical stopx

  stopx = .false.

  return
end
subroutine sumit ( d, fx, g, iv, liv, lv, n, v, x)

!*******************************************************************************
!
!! SUMIT carries out unconstrained minimization iterations for SUMSL.
!
!  Discussion:
!
!    The routine uses double-dogleg/BFGS steps.
!
!    parameters iv, n, v, and x are the same as the corresponding
!    ones to sumsl (which see), except that v can be shorter (since
!    the part of v that sumsl uses for storing g is not needed).
!    moreover, compared with sumsl, iv(1) may have the two additional
!    output values 1 and 2, which are explained below, as is the use
!    of iv(toobig) and iv(nfgcal).  the value iv(g), which is an
!    output value from sumsl (and smsno), is not referenced by
!    sumit or the subroutines it calls.
!
!    fx and g need not have been initialized when sumit is called
!    with iv(1) = 12, 13, or 14.
!
! iv(1) = 1 means the caller should set fx to f(x), the function value
!             at x, and call sumit again, having changed none of the
!             other parameters.  an exception occurs if f(x) cannot be
!             (e.g. if overflow would occur), which may happen because
!             of an oversized step.  in this case the caller should set
!             iv(toobig) = iv(2) to 1, which will cause sumit to ig-
!             nore fx and try a smaller step.  the parameter nf that
!             sumsl passes to calcf (for possible use by calcg) is a
!             copy of iv(nfcall) = iv(6).
! iv(1) = 2 means the caller should set g to g(x), the gradient vector
!             of f at x, and call sumit again, having changed none of
!             the other parameters except possibly the scale vector d
!             when iv(dtype) = 0.  the parameter nf that sumsl passes
!             to calcg is iv(nfgcal) = iv(7).  if g(x) cannot be
!             evaluated, then the caller may set iv(nfgcal) to 0, in
!             which case sumit will return with iv(1) = 65.
!
!  Parameters:  
!
! d.... scale vector.
! fx... function value.
! g.... gradient vector.
! iv... integer value array.
! liv.. length of iv (at least 60).
! lv... length of v (at least 71 + n*(n+13)/2).
! n.... number of variables (components in x and g).
! v.... floating-point value array.
! x.... vector of parameters to be optimized.
!
  integer liv
  integer lv
  integer n

  integer iv(liv)
  real ( kind = 8 ) d(n)
  real ( kind = 8 ) fx
  real ( kind = 8 ) g(n)
  real ( kind = 8 ) v(lv)
  real ( kind = 8 ) x(n)
  integer dg1, g01, i, k, l, lstgst, nwtst1, step1
  integer        temp1, w, x01, z
  real ( kind = 8 ) t
  real ( kind = 8 ) half, negone, one, onep2, zero
  logical stopx
  real ( kind = 8 ) dotprd, reldst, v2norm
  integer cnvcod, dg, dgnorm, dinit, dstnrm, dst0, f, f0, fdif
  integer gthg, gtstep, g0, incfac, inith, irc, kagqt, lmat, lmax0
  integer lmaxs, mode, model, mxfcal, mxiter, nextv, nfcall, nfgcal
  integer ngcall, niter, nreduc, nwtstp, preduc, radfac, radinc
  integer radius, rad0, reldx, restor, step, stglim, stlstg, toobig
  integer tuner4, tuner5, vneed, xirc, x0

  parameter (cnvcod=55, dg=37, g0=48, inith=25, irc=29, kagqt=33 )
  parameter ( mode=35, model=5, mxfcal=17, mxiter=18, nfcall=6 )
  parameter ( nfgcal=7, ngcall=30, niter=31, nwtstp=34, radinc=8 )
  parameter ( restor=9, step=40, stglim=11, stlstg=41, toobig=2 )
  parameter ( vneed=4, xirc=13, x0=43)

  parameter (dgnorm=1, dinit=38, dstnrm=2, dst0=3, f=10, f0=13 )
  parameter ( fdif=11, gthg=44, gtstep=4, incfac=23, lmat=42 )
  parameter ( lmax0=35, lmaxs=36, nextv=47, nreduc=6, preduc=7 )
  parameter ( radfac=16, radius=8, rad0=9, reldx=17, tuner4=29 )
  parameter ( tuner5=30)

  parameter (half=0.5d+0, negone=-1.d+0, one=1.d+0, onep2=1.2d+0, zero=0.d+0)
!
  i = iv(1)
  if (i == 1) go to 50
  if (i == 2) go to 60
!
!   check validity of iv and v input values 
!
  if (iv(1) == 0) call deflt(2, iv, liv, lv, v)
  if (iv(1) == 12 .or. iv(1) == 13) then
    iv(vneed) = iv(vneed) + n*(n+13)/2
  end if
  call parck(2, d, iv, liv, lv, n, v)
  i = iv(1) - 2
  if (i > 12) then
    return
  end if
  go to (180, 180, 180, 180, 180, 180, 120, 90, 120, 10, 10, 20), i
!
!   storage allocation 
!
10    l = iv(lmat)
  iv(x0) = l + n*(n+1)/2
  iv(step) = iv(x0) + n
  iv(stlstg) = iv(step) + n
  iv(g0) = iv(stlstg) + n
  iv(nwtstp) = iv(g0) + n
  iv(dg) = iv(nwtstp) + n
  iv(nextv) = iv(dg) + n
  if (iv(1) /= 13) go to 20
     iv(1) = 14
     return
!
!   initialization 
!
 20   iv(niter) = 0
  iv(nfcall) = 1
  iv(ngcall) = 1
  iv(nfgcal) = 1
  iv(mode) = -1
  iv(model) = 1
  iv(stglim) = 1
  iv(toobig) = 0
  iv(cnvcod) = 0
  iv(radinc) = 0
  v(rad0) = 0.0D+00
  if (v(dinit) >= 0.0D+00) call vscopy(n, d, v(dinit))
  if (iv(inith) /= 1) go to 40
!
!  set the initial hessian approximation to diag(d)**-2 
!
     l = iv(lmat)
     call vscopy(n*(n+1)/2, v(l), zero)
     k = l - 1

     do i = 1, n
       k = k + i
       t = d(i)
       if (t <= 0.0D+00) t = one
       v(k) = t
     end do
!
!  compute initial function value 
!
 40   iv(1) = 1
  return

 50   v(f) = fx
  if (iv(mode) >= 0) go to 180
  iv(1) = 2
  if (iv(toobig) == 0) then
    return
  end if
     iv(1) = 63
     go to 300
!
!   make sure gradient could be computed 
!
 60   if (iv(nfgcal) /= 0) go to 70
     iv(1) = 65
     go to 300

 70   dg1 = iv(dg)
  call vvmulp(n, v(dg1), g, d, -1)
  v(dgnorm) = v2norm(n, v(dg1))

  if (iv(cnvcod) /= 0) go to 290
  if (iv(mode) == 0) go to 250
!
!   allow first step to have scaled 2-norm at most v(lmax0) 
!
  v(radius) = v(lmax0)

  iv(mode) = 0
!
!  main loop  
!
!   print iteration summary, check iteration limit 
!
 80   call itsum(d, g, iv, liv, lv, n, v, x)
 90   k = iv(niter)
  if (k < iv(mxiter)) go to 100
     iv(1) = 10
     go to 300
!
!   update radius 
!
 100  iv(niter) = k + 1
  if(k>0)v(radius) = v(radfac) * v(dstnrm)
!
!   initialize for start of next iteration 
!
  g01 = iv(g0)
  x01 = iv(x0)
  v(f0) = v(f)
  iv(irc) = 4
  iv(kagqt) = -1
!
!      copy x to x0, g to g0 
!
  call vcopy(n, v(x01), x)
  call vcopy(n, v(g01), g)
!
!  Check STOPX and function evaluation limit 
!
 110  if ( .not. stopx ( ) ) go to 130
     iv(1) = 11
     go to 140
!
!  Come here when restarting after func. eval. limit or STOPX.
!
 120  if (v(f) >= v(f0)) go to 130
     v(radfac) = one
     k = iv(niter)
     go to 100

 130  if (iv(nfcall) < iv(mxfcal)) go to 150
     iv(1) = 9
 140     if (v(f) >= v(f0)) go to 300
!
!  in case of STOPX or function evaluation limit with
!  improved v(f), evaluate the gradient at x.
!
          iv(cnvcod) = iv(1)
          go to 240
!
!  Compute candidate step  
!
 150  step1 = iv(step)
  dg1 = iv(dg)
  nwtst1 = iv(nwtstp)
  if (iv(kagqt) >= 0) go to 160
     l = iv(lmat)
     call livmul(n, v(nwtst1), v(l), g)
     v(nreduc) = half * dotprd(n, v(nwtst1), v(nwtst1))
     call litvmu(n, v(nwtst1), v(l), v(nwtst1))
     call vvmulp(n, v(step1), v(nwtst1), d, 1)
     v(dst0) = v2norm(n, v(step1))
     call vvmulp(n, v(dg1), v(dg1), d, -1)
     call ltvmul(n, v(step1), v(l), v(dg1))
     v(gthg) = v2norm(n, v(step1))
     iv(kagqt) = 0
 160  call dbdog(v(dg1), lv, n, v(nwtst1), v(step1), v)
  if (iv(irc) == 6) go to 180
!
!   check whether evaluating f(x0 + step) looks worthwhile 
!
  if (v(dstnrm) <= 0.0D+00) go to 180
  if (iv(irc) /= 5) go to 170
  if (v(radfac) <= one) go to 170
  if (v(preduc) <= onep2 * v(fdif)) go to 180
!
!  Compute f(x0 + step) 
!
 170  x01 = iv(x0)
  step1 = iv(step)
  call vaxpy(n, x, one, v(step1), v(x01))
  iv(nfcall) = iv(nfcall) + 1
  iv(1) = 1
  iv(toobig) = 0
  return
!
!  Assess candidate step.
!
 180  x01 = iv(x0)
  v(reldx) = reldst(n, d, x, v(x01))
  call assst(iv, liv, lv, v)
  step1 = iv(step)
  lstgst = iv(stlstg)
  if (iv(restor) == 1) call vcopy(n, x, v(x01))
  if (iv(restor) == 2) call vcopy(n, v(lstgst), v(step1))
  if (iv(restor) /= 3) go to 190
     call vcopy(n, v(step1), v(lstgst))
     call vaxpy(n, x, one, v(step1), v(x01))
     v(reldx) = reldst(n, d, x, v(x01))

 190  k = iv(irc)
  go to (200,230,230,230,200,210,220,220,220,220,220,220,280,250), k
!
!      recompute step with changed radius 
!
 200     v(radius) = v(radfac) * v(dstnrm)
     go to 110
!
!   compute step of length v(lmaxs) for singular convergence test.
!
 210  v(radius) = v(lmaxs)
  go to 150
!
!   convergence or false convergence 
!
 220  iv(cnvcod) = k - 4
  if (v(f) >= v(f0)) go to 290
     if (iv(xirc) == 14) go to 290
          iv(xirc) = 14
!
!  Process acceptable step.
!
 230  if (iv(irc) /= 3) go to 240
     step1 = iv(step)
     temp1 = iv(stlstg)
!
!      set  temp1 = hessian * step  for use in gradient tests 
!
     l = iv(lmat)
     call ltvmul(n, v(temp1), v(l), v(step1))
     call lvmul(n, v(temp1), v(l), v(temp1))
!
!   compute gradient 
!
 240  iv(ngcall) = iv(ngcall) + 1
  iv(1) = 2
  return
!
!   initializations -- g0 = g - g0, etc. 
!
 250  g01 = iv(g0)
  call vaxpy(n, v(g01), negone, v(g01), g)
  step1 = iv(step)
  temp1 = iv(stlstg)
  if (iv(irc) /= 3) go to 270
!
!   set v(radfac) by gradient tests 
!
!  Set  temp1 = diag(d)**-1 * (hessian*step + (g(x0)-g(x))) 
!
     call vaxpy(n, v(temp1), negone, v(g01), v(temp1))
     call vvmulp(n, v(temp1), v(temp1), d, -1)
!
!  Do gradient tests 
!
     if (v2norm(n, v(temp1)) <= v(dgnorm) * v(tuner4)) then
       go to 260
     end if

     if (dotprd(n, g, v(step1)) >= v(gtstep) * v(tuner5))  then
       go to 270
     end if

 260               v(radfac) = v(incfac)
!
!   update h, loop 
!
 270  w = iv(nwtstp)
  z = iv(x0)
  l = iv(lmat)
  call wzbfgs(v(l), n, v(step1), v(w), v(g01), v(z))
!
!  Use the n-vectors starting at v(step1) and v(g01) for scratch.
!
  call lupdat(v(temp1), v(step1), v(l), v(g01), v(l), n, v(w), v(z))
  iv(1) = 2
  go to 80
!
!   misc. details   
!
!   bad parameters to assess 
!
 280  iv(1) = 64
  go to 300
!
!  Print summary of final iteration and other requested items 
!
 290  iv(1) = iv(cnvcod)
  iv(cnvcod) = 0
 300  call itsum(d, g, iv, liv, lv, n, v, x)

  return
end
subroutine sumsl(n, d, x, calcf, calcg, iv, liv, lv, v, uiparm, urparm, ufparm)

!*******************************************************************************
!
!! SUMSL minimizes a general unconstrained objective function.
!
!  Discussion:
!
!    The routine uses analytic gradient and hessian approximation from 
!    the secant update.
!
!    This routine interacts with subroutine  sumit  in an attempt
!    to find an n-vector  x*  that minimizes the (unconstrained)
!    objective function computed by  calcf.  (often the  x*  found is
!    a local minimizer rather than a global one.)
!
!  Reference:
!
!    J E Dennis, David Gay, and R E Welsch,
!    An Adaptive Nonlinear Least-squares Algorithm,
!    ACM Transactions on Mathematical Software, 
!    Volume 7, Number 3, 1981.
!
!    J E Dennis, H H W Mei,
!    Two New Unconstrained Optimization Algorithms Which Use 
!    Function and Gradient Values, 
!    Journal of Optimization Theory and Applications,
!    Volume 28, pages 453-482, 1979.
!
!    J E Dennis, Jorge More,
!    Quasi-Newton Methods, Motivation and Theory, 
!    SIAM Review, 
!    Volume 19, pages 46-89, 1977.
!
!    D Goldfarb, 
!    Factorized Variable Metric Methods for Unconstrained Optimization,
!    Mathematics of Computation,
!    Volume 30, pages 796-811, 1976.
!
!  Parameters:
!
! n  (input) the number of variables on which  f  depends, i.e.,
!                  the number of components in  x.
! d  (input/output) a scale vector such that  d(i)*x(i),
!                  i = 1,2,...,n,  are all in comparable units.
!                  d can strongly affect the behavior of sumsl.
!                  finding the best choice of d is generally a trial-
!                  and-error process.  choosing d so that d(i)*x(i)
!                  has about the same value for all i often works well.
!                  the defaults provided by subroutine deflt (see iv
!                  below) require the caller to supply d.
! x........ (input/output) before (initially) calling sumsl, the call-
!                  er should set  x  to an initial guess at  x*.  when
!                  sumsl returns,  x  contains the best point so far
!                  found, i.e., the one that gives the least value so
!                  far seen for  f(x).
! calcf.... (input) a subroutine that, given x, computes f(x).  calcf
!                  must be declared external in the calling program.
!                  it is invoked by
!                       call calcf(n, x, nf, f, uiparm, urparm, ufparm)
!                  when calcf is called, nf is the invocation
!                  count for calcf.  nf is included for possible use
!                  with calcg.  if x is out of bounds (e.g., if it
!                  would cause overflow in computing f(x)), then calcf
!                  should set nf to 0.  this will cause a shorter step
!                  to be attempted.  (if x is in bounds, then calcf
!                  should not change nf.)  the other parameters are as
!                  described above and below.  calcf should not change
!                  n, p, or x.
! calcg.... (input) a subroutine that, given x, computes g(x), the gra-
!                  dient of f at x.  calcg must be declared external in
!                  the calling program.  it is invoked by
!                       call calcg(n, x, nf, g, uiparm, urparm, ufaprm)
!                  when calcg is called, nf is the invocation
!                  count for calcf at the time f(x) was evaluated.  the
!                  x passed to calcg is usually the one passed to calcf
!                  on either its most recent invocation or the one
!                  prior to it.  if calcf saves intermediate results
!                  for use by calcg, then it is possible to tell from
!                  nf whether they are valid for the current x (or
!                  which copy is valid if two copies are kept).  if g
!                  cannot be computed at x, then calcg should set nf to
!                  0.  in this case, sumsl will return with iv(1) = 65.
!                  (if g can be computed at x, then calcg should not
!                  changed nf.)  the other parameters to calcg are as
!                  described above and below.  calcg should not change
!                  n or x.
! iv....... (input/output) an integer value array of length liv (see
!                  below) that helps control the sumsl algorithm and
!                  that is used to store various intermediate quanti-
!                  ties.  of particular interest are the initialization/
!                  return code iv(1) and the entries in iv that control
!                  printing and limit the number of iterations and func-
!                  tion evaluations.  see the section on iv input
!                  values below.
! liv...... (input) length of iv array.  must be at least 60.  if liv
!                  is too small, then sumsl returns with iv(1) = 15.
!                  when sumsl returns, the smallest allowed value of
!                  liv is stored in iv(lastiv) -- see the section on
!                  iv output values below.  (this is intended for use
!                  with extensions of sumsl that handle constraints.)
! lv....... (input) length of v array.  must be at least 71+n*(n+15)/2.
!                  (at least 77+n*(n+17)/2 for smsno, at least
!                  78+n*(n+12) for humsl).  if lv is too small, then
!                  sumsl returns with iv(1) = 16.  when sumsl returns,
!                  the smallest allowed value of lv is stored in
!                  iv(lastv) -- see the section on iv output values
!                  below.
! v........ (input/output) a floating-point value array of length lv
!                  (see below) that helps control the sumsl algorithm
!                  and that is used to store various intermediate
!                  quantities.  of particular interest are the entries
!                  in v that limit the length of the first step
!                  attempted (lmax0) and specify convergence tolerances
!                  (afctol, lmaxs, rfctol, sctol, xctol, xftol).
! uiparm... (input) user integer parameter array passed without change
!                  to calcf and calcg.
! urparm... (input) user floating-point parameter array passed without
!                  change to calcf and calcg.
! ufparm... (input) user external subroutine or function passed without
!                  change to calcf and calcg.
!
!   iv input values (from subroutine deflt) 
!
! iv(1)...  on input, iv(1) should have a value between 0 and 14......
!             0 and 12 mean this is a fresh start.  0 means that
!                  deflt(2, iv, liv, lv, v)
!             is to be called to provide all default values to iv and
!             v.  12 (the value that deflt assigns to iv(1)) means the
!             caller has already called deflt and has possibly changed
!             some iv and/or v entries to non-default values.
!             13 means deflt has been called and that sumsl (and
!             sumit) should only do their storage allocation.  that is,
!             they should set the output components of iv that tell
!             where various subarrays arrays of v begin, such as iv(g)
!             (and, for humsl and humit only, iv(dtol)), and return.
!             14 means that a storage has been allocated (by a call
!             with iv(1) = 13) and that the algorithm should be
!             started.  when called with iv(1) = 13, sumsl returns
!             iv(1) = 14 unless liv or lv is too small (or n is not
!             positive).  default = 12.
! iv(inith).... iv(25) tells whether the hessian approximation h should
!             be initialized.  1 (the default) means sumit should
!             initialize h to the diagonal matrix whose i-th diagonal
!             element is d(i)**2.  0 means the caller has supplied a
!             cholesky factor  l  of the initial hessian approximation
!             h = l*(l**t)  in v, starting at v(iv(lmat)) = v(iv(42))
!             (and stored compactly by rows).  note that iv(lmat) may
!             be initialized by calling sumsl with iv(1) = 13 (see
!             the iv(1) discussion above).  default = 1.
! iv(mxfcal)... iv(17) gives the maximum number of function evaluations
!             (calls on calcf) allowed.  if this number does not suf-
!             fice, then sumsl returns with iv(1) = 9.  default = 200.
! iv(mxiter)... iv(18) gives the maximum number of iterations allowed.
!             it also indirectly limits the number of gradient evalua-
!             tions (calls on calcg) to iv(mxiter) + 1.  if iv(mxiter)
!             iterations do not suffice, then sumsl returns with
!             iv(1) = 10.  default = 150.
! iv(outlev)... iv(19) controls the number and length of iteration sum-
!             mary lines printed (by itsum).  iv(outlev) = 0 means do
!             not print any summary lines.  otherwise, print a summary
!             line after each abs(iv(outlev)) iterations.  if iv(outlev)
!             is positive, then summary lines of length 78 (plus carri-
!             age control) are printed, including the following...  the
!             iteration and function evaluation counts, f = the current
!             function value, relative difference in function values
!             achieved by the latest step (i.e., reldf = (f0-v(f))/f01,
!             where f01 is the maximum of abs(v(f)) and abs(v(f0)) and
!             v(f0) is the function value from the previous itera-
!             tion), the relative function reduction predicted for the
!             step just taken (i.e., preldf = v(preduc) / f01, where
!             v(preduc) is described below), the scaled relative change
!             in x (see v(reldx) below), the step parameter for the
!             step just taken (stppar = 0 means a full newton step,
!             between 0 and 1 means a relaxed newton step, between 1
!             and 2 means a double dogleg step, greater than 2 means
!             a scaled down Cauchy step -- see subroutine dbldog), the
!             2-norm of the scale vector d times the step just taken
!             (see v(dstnrm) below), and npreldf, i.e.,
!             v(nreduc)/f01, where v(nreduc) is described below -- if
!             npreldf is positive, then it is the relative function
!             reduction predicted for a newton step (one with
!             stppar = 0).  if npreldf is negative, then it is the
!             negative of the relative function reduction predicted
!             for a step computed with step bound v(lmaxs) for use in
!             testing for singular convergence.
!                  if iv(outlev) is negative, then lines of length 50
!             are printed, including only the first 6 items listed
!             above (through reldx).
!             default = 1.
! iv(parprt)... iv(20) = 1 means print any nondefault v values on a
!             fresh start or any changed v values on a restart.
!             iv(parprt) = 0 means skip this printing.  default = 1.
! iv(prunit)... iv(21) is the output unit number on which all printing
!             is done.  iv(prunit) = 0 means suppress all printing.
!             default = standard output unit (unit 6 on most systems).
! iv(solprt)... iv(22) = 1 means print out the value of x returned (as
!             well as the gradient and the scale vector d).
!             iv(solprt) = 0 means skip this printing.  default = 1.
! iv(statpr)... iv(23) = 1 means print summary statistics upon return-
!             ing.  these consist of the function value, the scaled
!             relative change in x caused by the most recent step (see
!             v(reldx) below), the number of function and gradient
!             evaluations (calls on calcf and calcg), and the relative
!             function reductions predicted for the last step taken and
!             for a newton step (or perhaps a step bounded by v(lmaxs)
!             -- see the descriptions of preldf and npreldf under
!             iv(outlev) above).
!             iv(statpr) = 0 means skip this printing.
!             iv(statpr) = -1 means skip this printing as well as that
!             of the one-line termination reason message.  default = 1.
! iv(x0prt).... iv(24) = 1 means print the initial x and scale vector d
!             (on a fresh start only).  iv(x0prt) = 0 means skip this
!             printing.  default = 1.
!
!   (selected) iv output values 
!
! iv(1)........ on output, iv(1) is a return code....
!             3 = x-convergence.  the scaled relative difference (see
!                  v(reldx)) between the current parameter vector x and
!                  a locally optimal parameter vector is very likely at
!                  most v(xctol).
!             4 = relative function convergence.  the relative differ-
!                  ence between the current function value and its lo-
!                  cally optimal value is very likely at most v(rfctol).
!             5 = both x- and relative function convergence (i.e., the
!                  conditions for iv(1) = 3 and iv(1) = 4 both hold).
!             6 = absolute function convergence.  the current function
!                  value is at most v(afctol) in absolute value.
!             7 = singular convergence.  the hessian near the current
!                  iterate appears to be singular or nearly so, and a
!                  step of length at most v(lmaxs) is unlikely to yield
!                  a relative function decrease of more than v(sctol).
!             8 = false convergence.  the iterates appear to be converg-
!                  ing to a noncritical point.  this may mean that the
!                  convergence tolerances (v(afctol), v(rfctol),
!                  v(xctol)) are too small for the accuracy to which
!                  the function and gradient are being computed, that
!                  there is an error in computing the gradient, or that
!                  the function or gradient is discontinuous near x.
!             9 = function evaluation limit reached without other con-
!                  vergence (see iv(mxfcal)).
!            10 = iteration limit reached without other convergence
!                  (see iv(mxiter)).
!            11 = STOPX returned .true. (external interrupt).  see the
!                  usage notes below.
!            14 = storage has been allocated (after a call with
!                  iv(1) = 13).
!            17 = restart attempted with n changed.
!            18 = d has a negative component and iv(dtype) <= 0.
!            19...43 = v(iv(1)) is out of range.
!            63 = f(x) cannot be computed at the initial x.
!            64 = bad parameters passed to assess (which should not
!                  occur).
!            65 = the gradient could not be computed at x (see calcg
!                  above).
!            67 = bad first parameter to deflt.
!            80 = iv(1) was out of range.
!            81 = n is not positive.
! iv(g)........ iv(28) is the starting subscript in v of the current
!             gradient vector (the one corresponding to x).
! iv(lastiv)... iv(44) is the least acceptable value of liv.  (it is
!             only set if liv is at least 44.)
! iv(lastv).... iv(45) is the least acceptable value of lv.  (it is
!             only set if liv is large enough, at least iv(lastiv).)
! iv(nfcall)... iv(6) is the number of calls so far made on calcf (i.e.,
!             function evaluations).
! iv(ngcall)... iv(30) is the number of gradient evaluations (calls on
!             calcg).
! iv(niter).... iv(31) is the number of iterations performed.
!
!   (selected) v input values (from subroutine deflt) 
!
! v(bias)..... v(43) is the bias parameter used in subroutine dbldog --
!             see that subroutine for details.  default = 0.8.
! v(afctol)... v(31) is the absolute function convergence tolerance.
!             if sumsl finds a point where the function value is less
!             than v(afctol) in absolute value, and if sumsl does not
!             return with iv(1) = 3, 4, or 5, then it returns with
!             iv(1) = 6.  this test can be turned off by setting
!             v(afctol) to zero.  default = max(10**-20, machep**2),
!             where machep is the unit roundoff.
! v(dinit).... v(38), if nonnegative, is the value to which the scale
!             vector d is initialized.  default = -1.
! v(lmax0).... v(35) gives the maximum 2-norm allowed for d times the
!             very first step that sumsl attempts.  this parameter can
!             markedly affect the performance of sumsl.
! v(lmaxs).... v(36) is used in testing for singular convergence -- if
!             the function reduction predicted for a step of length
!             bounded by v(lmaxs) is at most v(sctol) * abs(f0), where
!             f0  is the function value at the start of the current
!             iteration, and if sumsl does not return with iv(1) = 3,
!             4, 5, or 6, then it returns with iv(1) = 7.  default = 1.
! v(rfctol)... v(32) is the relative function convergence tolerance.
!             if the current model predicts a maximum possible function
!             reduction (see v(nreduc)) of at most v(rfctol)*abs(f0)
!             at the start of the current iteration, where  f0  is the
!             then current function value, and if the last step attempt-
!             ed achieved no more than twice the predicted function
!             decrease, then sumsl returns with iv(1) = 4 (or 5).
!             default = max(10**-10, machep**(2/3)), where machep is
!             the unit roundoff.
! v(sctol).... v(37) is the singular convergence tolerance -- see the
!             description of v(lmaxs) above.
! v(tuner1)... v(26) helps decide when to check for false convergence.
!             this is done if the actual function decrease from the
!             current step is no more than v(tuner1) times its predict-
!             ed value.  default = 0.1.
! v(xctol).... v(33) is the x-convergence tolerance.  if a newton step
!             (see v(nreduc)) is tried that has v(reldx) <= v(xctol)
!             and if this step yields at most twice the predicted func-
!             tion decrease, then sumsl returns with iv(1) = 3 (or 5).
!             (see the description of v(reldx) below.)
!             default = machep**0.5, where machep is the unit roundoff.
! v(xftol).... v(34) is the false convergence tolerance.  if a step is
!             tried that gives no more than v(tuner1) times the predict-
!             ed function decrease and that has v(reldx) <= v(xftol),
!             and if sumsl does not return with iv(1) = 3, 4, 5, 6, or
!             7, then it returns with iv(1) = 8.  (see the description
!             of v(reldx) below.)  default = 100*machep, where
!             machep is the unit roundoff.
! v(*)........ deflt supplies to v a number of tuning constants, with
!             which it should ordinarily be unnecessary to tinker.  see
!             section 17 of version 2.2 of the nl2sol usage summary
!             (i.e., the appendix to ref. 1) for details on v(i),
!             i = decfac, incfac, phmnfc, phmxfc, rdfcmn, rdfcmx,
!             tuner2, tuner3, tuner4, tuner5.
!
!   (selected) v output values 
!
! v(dgnorm)... v(1) is the 2-norm of (diag(d)**-1)*g, where g is the
!             most recently computed gradient.
! v(dstnrm)... v(2) is the 2-norm of diag(d)*step, where step is the
!             current step.
! v(f)........ v(10) is the current function value.
! v(f0)....... v(13) is the function value at the start of the current
!             iteration.
! v(nreduc)... v(6), if positive, is the maximum function reduction
!             possible according to the current model, i.e., the func-
!             tion reduction predicted for a newton step (i.e.,
!             step = -h**-1 * g,  where  g  is the current gradient and
!             h is the current hessian approximation).
!                  if v(nreduc) is negative, then it is the negative of
!             the function reduction predicted for a step computed with
!             a step bound of v(lmaxs) for use in testing for singular
!             convergence.
! v(preduc)... v(7) is the function reduction predicted (by the current
!             quadratic model) for the current step.  this (divided by
!             v(f0)) is used in testing for relative function
!             convergence.
! v(reldx).... v(17) is the scaled relative change in x caused by the
!             current step, computed as
!                  max(abs(d(i)*(x(i)-x0(i)), 1 <= i <= p) /
!                     max(d(i)*(abs(x(i))+abs(x0(i))), 1 <= i <= p),
!             where x = x0 + step.
!
!  notes 
!
!   algorithm notes 
!
!        this routine uses a hessian approximation computed from the
!     bfgs update (see ref 3).  only a cholesky factor of the hessian
!     approximation is stored, and this is updated using ideas from
!     ref. 4.  steps are computed by the double dogleg scheme described
!     in ref. 2.  the steps are assessed as in ref. 1.
!
!   usage notes 
!
!        after a return with iv(1) <= 11, it is possible to restart,
!     i.e., to change some of the iv and v input values described above
!     and continue the algorithm from the point where it was interrupt-
!     ed.  iv(1) should not be changed, nor should any entries of iv
!     and v other than the input values (those supplied by deflt).
!        those who do not wish to write a calcg which computes the
!     gradient analytically should call smsno rather than sumsl.
!     smsno uses finite differences to compute an approximate gradient.
!        those who would prefer to provide f and g (the function and
!     gradient) by reverse communication rather than by writing subrou-
!     tines calcf and calcg may call on sumit directly.  see the com-
!     ments at the beginning of sumit.
!        those who use sumsl interactively may wish to supply their
!     own STOPX function, which should return .true. if the break key
!     has been pressed since STOPX was last invoked.  this makes it
!     possible to externally interrupt sumsl (which will return with
!     iv(1) = 11 if STOPX returns .true.).
!        storage for g is allocated at the end of v.  thus the caller
!     may make v longer than specified above and may allow calcg to use
!     elements of g beyond the first n as scratch storage.
!
!   portability notes 
!
!        the sumsl distribution tape contains both single- and double-
!     precision versions of the sumsl source code, so it should be un-
!     necessary to change precisions.
!        only the function rmdcon contains machine-dependent
!     constants.  to change from one machine to another, it should
!     suffice to change the (few) relevant lines in these functions.
!        intrinsic functions are explicitly declared.  on certain com-
!     puters (e.g. univac), it may be necessary to comment out these
!     declarations.  so that this may be done automatically by a simple
!     program, such declarations are preceded by a comment having c/+
!     in columns 1-3 and blanks in columns 4-72 and are followed by
!     a comment having c/ in columns 1 and 2 and blanks in columns 3-72.
!        the sumsl source code is expressed in 1966 ansi standard
!     fortran.  it may be converted to fortran 77 by commenting out all
!     lines that fall between a line having c/6 in columns 1-3 and a
!     line having c/7 in columns 1-3 and by removing (i.e., replacing
!     by a blank) the c in column 1 of the lines that follow the c/7
!     line and precede a line having c/ in columns 1-2 and blanks in
!     columns 3-72.  these changes convert some data statements into
!     parameter statements, convert some variables from real to
!     character*4, and make the data statements that initialize these
!     variables use character strings delimited by primes instead
!     of hollerith constants.  (such variables and data statements
!     appear only in modules itsum and parck.  parameter statements
!     appear nearly everywhere.)  these changes also add save state-
!     ments for variables given machine-dependent constants by rmdcon.
!
  integer n, liv, lv
  integer iv(liv), uiparm(*)
  real ( kind = 8 ) d(n), x(n), v(lv), urparm(*)
!     dimension v(71 + n*(n+15)/2), uiparm(*), urparm(*)
 
  integer g1, iv1, nf
  real ( kind = 8 ) f
  integer nextv, nfcall, nfgcal, g, toobig, vneed

  parameter (nextv=47, nfcall=6, nfgcal=7, g=28, toobig=2, vneed=4)

  external ufparm

  if (iv(1) == 0) call deflt(2, iv, liv, lv, v)
  iv1 = iv(1)
  if (iv1 == 12 .or. iv1 == 13) iv(vneed) = iv(vneed) + n
  if (iv1 == 14) go to 10
  if (iv1 > 2 .and. iv1 < 12) go to 10
  g1 = 1
  if (iv1 == 12) iv(1) = 13
  go to 20

 10   g1 = iv(g)

 20   call sumit(d, f, v(g1), iv, liv, lv, n, v, x)
  if (iv(1) - 2) 30, 40, 50

 30   nf = iv(nfcall)
  call calcf(n, x, nf, f, uiparm, urparm, ufparm)
  if (nf <= 0) iv(toobig) = 1
  go to 20

 40   call calcg(n, x, iv(nfgcal), v(g1), uiparm, urparm, ufparm)
  go to 20

 50   if (iv(1) /= 14) then
        return
      end if
!
!  Storage allocation
!
  iv(g) = iv(nextv)
  iv(nextv) = iv(g) + n
  if (iv1 /= 13) go to 10

  return
end
function v2norm ( p, x )

!*******************************************************************************
!
!! V2NORM returns the 2-norm of the p-vector X.
!
!  Discussion:
!
!    The routine tries to avoid underflow. 
!
!  Parameters:
!
  integer p

  real ( kind = 8 ) x(p)
  integer i, j
  real ( kind = 8 ) r, scale
  real ( kind = 8 ), save :: sqteta = 0.0D+00
  real ( kind = 8 ) t, xi
  real ( kind = 8 ) rmdcon
  real ( kind = 8 ) v2norm

  v2norm = 0.0D+00

  if (p <= 0 ) then
    return
  end if

  if ( all ( x(1:p) == 0.0D+00 ) ) then
    return
  end if

  scale = 0.0D+00
  do i = 1, p
    if ( x(i) /= 0.0D+00 ) then
      scale = abs(x(i))
      exit
    end if
  end do

  if ( scale == 0.0D+00 ) then
    return
  end if

  if ( p <= i ) then
    v2norm = scale
    return
  end if

  t = 1.0D+00
  if ( sqteta == 0.0D+00 ) then
    sqteta = rmdcon(2)
  end if
!
!  sqteta is (slightly larger than) the square root of the
!  smallest positive floating point number on the machine.
!  the tests involving sqteta are done to prevent underflows.
!
  j = i + 1
  do i = j, p
    xi = abs(x(i))
    if (xi <= scale) then
      r = xi / scale
      if (r > sqteta) t = t + r*r
    else
      r = scale / xi
      if (r <= sqteta) r = 0.0D+00
      t = 1.0D+00  +  t * r*r
      scale = xi
    end if
  end do

  v2norm = scale * sqrt(t)

  return
end
subroutine vaxpy ( p, w, a, x, y )

!*******************************************************************************
!
!! VAXPY sets w = a*x + y.
!
!  Discussion:
!
!    w, x, y = p-vectors, a = scalar 
!
!  Parameters:
!
  implicit none

  integer p

  real ( kind = 8 ) a
  real ( kind = 8 ) w(p)
  real ( kind = 8 ) x(p)
  real ( kind = 8 ) y(p)

  w(1:p) = a * x(1:p) + y(1:p)

  return
end
subroutine vcopy ( p, y, x )

!*******************************************************************************
!
!! VCOPY sets y = x.
!
!  Discussion:
!
!    x and y are p-vectors 
!
  implicit none

  integer p
  real ( kind = 8 ) x(p)
  real ( kind = 8 ) y(p)

  y(1:p) = x(1:p)

  return
end
subroutine vdflt ( alg, lv, v )

!*******************************************************************************
!
!! VDFLT supplies default values to V.
!
!  Discussion:
!
!    alg = 1 means regression constants.
!    alg = 2 means general unconstrained optimization constants.
!
  implicit none

  integer alg, lv
  real ( kind = 8 ) v(lv)
  real ( kind = 8 ) rmdcon
  real ( kind = 8 ) machep, mepcrt, one, sqteps, three
  integer afctol, bias, cosmin, decfac, delta0, dfac, dinit, dltfdc
  integer dltfdj, dtinit, d0init, epslon, eta0, fuzz, huberc
  integer incfac, lmax0, lmaxs, phmnfc, phmxfc, rdfcmn, rdfcmx
  integer rfctol, rlimit, rsptol, sctol, sigmin, tuner1, tuner2
  integer tuner3, tuner4, tuner5, xctol, xftol

  parameter (one=1.d+0, three=3.d+0)

  parameter (afctol=31, bias=43, cosmin=47, decfac=22, delta0=44 )
  parameter ( dfac=41, dinit=38, dltfdc=42, dltfdj=43, dtinit=39 )
  parameter ( d0init=40, epslon=19, eta0=42, fuzz=45, huberc=48 )
  parameter ( incfac=23, lmax0=35, lmaxs=36, phmnfc=20, phmxfc=21 )
  parameter ( rdfcmn=24, rdfcmx=25, rfctol=32, rlimit=46, rsptol=49 )
  parameter ( sctol=37, sigmin=50, tuner1=26, tuner2=27, tuner3=28 )
  parameter ( tuner4=29, tuner5=30, xctol=33, xftol=34)

  machep = rmdcon(3)
  v(afctol) = 1.d-20

  if ( machep > 1.d-10 ) then
    v(afctol) = machep**2
  end if

  v(decfac) = 0.5d+0
  sqteps = rmdcon(4)
  v(dfac) = 0.6d+0
  v(delta0) = sqteps
  v(dtinit) = 1.d-6
  mepcrt = machep ** (one/three)
  v(d0init) = 1.d+0
  v(epslon) = 0.1d+0
  v(incfac) = 2.d+0
  v(lmax0) = 1.d+0
  v(lmaxs) = 1.d+0
  v(phmnfc) = -0.1d+0
  v(phmxfc) = 0.1d+0
  v(rdfcmn) = 0.1d+0
  v(rdfcmx) = 4.d+0
  v(rfctol) = max (1.d-10, mepcrt**2)
  v(sctol) = v(rfctol)
  v(tuner1) = 0.1d+0
  v(tuner2) = 1.d-4
  v(tuner3) = 0.75d+0
  v(tuner4) = 0.5d+0
  v(tuner5) = 0.75d+0
  v(xctol) = sqteps
  v(xftol) = 1.d+2 * machep

  if ( alg < 2 ) then
    v(cosmin) = max (1.d-6, 1.d+2 * machep)
    v(dinit) = 0.d+0
    v(dltfdc) = mepcrt
    v(dltfdj) = sqteps
    v(fuzz) = 1.5d+0
    v(huberc) = 0.7d+0
    v(rlimit) = rmdcon(5)
    v(rsptol) = 1.d-3
    v(sigmin) = 1.d-4
  else
    v(bias) = 0.8d+0
    v(dinit) = -1.0d+0
    v(eta0) = 1.0d+3 * machep
  end if

  return
end
subroutine vscopy ( p, y, s )

!*******************************************************************************
!
!! VSCOPY sets the vector Y to scalar S.
!
  implicit none

  integer p

  real ( kind = 8 ) s
  real ( kind = 8 ) y(p)

  y(1:p) = s

  return
end
subroutine vvmulp ( n, x, y, z, k )

!*******************************************************************************
!
!! VVMULP sets x(i) = y(i) * z(i)**k, 1 <= i <= n (for k = 1 or -1) 
!
  implicit none

  integer n

  integer k
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) z(n)

  if ( k < 0 ) then
    x(1:n) = y(1:n) / z(1:n)
  else
    x(1:n) = y(1:n) * z(1:n)
  end if

  return
end
subroutine wzbfgs ( l, n, s, w, y, z )

!*******************************************************************************
!
!! WZBFGS compute Y and Z for LUPDAT corresponding to BFGS update.
!
!  Discussion:
!
!    When S is computed in certain ways, for example by GQTSTP or
!    DBLDOG, it is possible to save N**2/2 operations since L'*S
!    or L*L'*S is then known.
!
!    If the BFGS update to L*L' would reduce its determinant to
!    less than EPS times its old value, then this routine in effect
!    replaces Y by THETA*Y + (1-THETA)*L*L'*S, where THETA
!    (between 0 and 1) is chosen to make the reduction factor = EPS.
!
!  Parameters:
!
!    l (i/o) cholesky factor of hessian, a lower triang. matrix stored
!             compactly by rows.
!
!    n (input) order of  l  and length of  s,  w,  y,  z.
!
!    s (input) the step just taken.
!
!    w (output) right singular vector of rank 1 correction to l.
!
!    y (input) change in gradients corresponding to s.
!
!    z (output) left singular vector of rank 1 correction to l.
!
  implicit none

  integer n

  real ( kind = 8 ) dotprd
  real ( kind = 8 ) cs
  real ( kind = 8 ) cy
  real ( kind = 8 ), parameter :: eps = 0.1D+00
  real ( kind = 8 ) epsrt
  real ( kind = 8 ) l(n*(n+1)/2)
  real ( kind = 8 ) s(n)
  real ( kind = 8 ) shs
  real ( kind = 8 ) theta
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) ys
  real ( kind = 8 ) z(n)

  call ltvmul ( n, w, l, s )
  shs = dotprd ( n, w, w )
  ys = dotprd ( n, y, s )

  if ( ys < eps * shs ) then
    theta = ( 1.0D+00 - eps ) * shs / ( shs - ys )
    epsrt = sqrt ( eps )
    cy = theta / ( shs * epsrt )
    cs = ( 1.0D+00 + ( theta - 1.0D+00 ) / epsrt ) / shs
  else
    cy = 1.0D+00 / ( sqrt ( ys ) * sqrt ( shs ) )
    cs = 1.0D+00 / shs
  end if

  call livmul ( n, z, l, y )

  z(1:n) = cy * z(1:n) - cs * w(1:n)

  return
end
