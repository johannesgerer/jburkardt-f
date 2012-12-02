program main

!*****************************************************************************80
!
!! MAIN is the main program for GEOMPACK3_PRB.
!
!  Discussion:
!
!    GEOMPACK3_PRB tests the routines in GEOMPACK3.
!
!  Modified:
!
!    18 July 2009
!
  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GEOMPACK3_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the GEOMPACK3 library.'

  call test_angle ( )
  call test_area ( )
  call test_dhpsrt ( )
  call test_dtris3 ( )
  call test_ihpsrt ( )
  call test_lu ( )
  call test_meas ( )
  call test_prime ( )
  call test_rotiar ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GEOMPACK3_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test_angle ( )

!*****************************************************************************80
!
!! TEST_ANGLE tests ANGLE.
!
!  Modified:
!
!    25 August 2005
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 6

  real ( kind = 8 ) ang(n)
  real ( kind = 8 ) angle
  real ( kind = 8 ) atr(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) im1
  integer ( kind = 4 ) ip1
  real ( kind = 8 ) xc(n)
  real ( kind = 8 ) yc(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_ANGLE'
  write ( *, '(a)' ) '  ANGLE computes the angles of a polygon.'

  call hexagon_vertices_2d ( xc, yc )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    I    X             Y'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(i5,4g14.6)' ) i, xc(i), yc(i)
  end do

  do i = 1, n

    im1 = i4_wrap ( i - 1, 1, n )
    ip1 = i4_wrap ( i + 1, 1, n )

    ang(i) = angle ( xc(im1), yc(im1), xc(i), yc(i), xc(ip1), yc(ip1) )

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      I       Angle'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i5,2x,g14.6)' ) i, ang(i)
  end do

  return
end
subroutine test_area ( )

!*****************************************************************************80
!
!! TEST_AREA tests AREAPG and AREATR.
!
!  Modified:
!
!    25 August 2005
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 6

  real ( kind = 8 ) area
  real ( kind = 8 ) areapg
  real ( kind = 8 ) areatr
  real ( kind = 8 ) atr(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) sum2
  real ( kind = 8 ) xc(n)
  real ( kind = 8 ) yc(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_AREA'
  write ( *, '(a)' ) '  AREAPG computes the area of a polygon.'
  write ( *, '(a)' ) '  AREATR computes the area of a triangle.'

  call hexagon_vertices_2d ( xc, yc )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I, (X,Y)'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(i5,4g14.6)' ) i, xc(i), yc(i)
  end do

  area = areapg ( n, xc, yc )

  sum2 = 0.0D+00

  atr(1) = 0.0D+00

  do i = 2, n-1

    atr(i) = areatr ( xc(1), yc(1), xc(i), yc(i), xc(i+1), yc(i+1) )

    sum2 = sum2 + atr(i)

  end do

  atr(n) = 0.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I, Triangle'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i5,2x,g14.6)' ) i, atr(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Area computation by:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  AREAPG   AREATR'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,g14.6,2x,g14.6)' ) area, sum2

  return
end
subroutine test_bctr ( )

!*****************************************************************************80
!
!! TEST_BCTR tests ...
!
!     Driver program for testing routines TRIPR3,TRIBFC,DSCONV,INTMVG,
!        BCDTRI,SWPREM,UPDATR,LSRCT3.
!     Routines called:
!        CCSPH,CVDEC3,CVDECF,DSCONV,DSPHDC,DSPHFH,DSPHIH,EMNRTH,GTIME,
!        PRIME,RADRTH,RMCLED,RMCPFC,SANGMN,STATS,TETLST,TRIPR3,
!        VOLCPH,VOLTH,WIDTH3
!
  implicit none

  real ( kind = 8 ) d_pi
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ), parameter :: maxbt = 2000
  integer ( kind = 4 ), parameter :: maxfc = 10000
  integer ( kind = 4 ), parameter :: maxfp = 800
  integer ( kind = 4 ), parameter :: maxfv = 3500
  integer ( kind = 4 ), parameter :: maxhf = 200
  integer ( kind = 4 ), parameter :: maxht = 6011
  integer ( kind = 4 ), parameter :: maxiw = 7500
  integer ( kind = 4 ), parameter :: maxnt = 5000
  integer ( kind = 4 ), parameter :: maxpf = 1200
  integer ( kind = 4 ), parameter :: maxte = 2500
  integer ( kind = 4 ), parameter :: maxvc = 2000
  integer ( kind = 4 ), parameter :: maxvm = 4000
  integer ( kind = 4 ), parameter :: maxwk = 7000
  integer ( kind = 4 ) msglvl
  integer ( kind = 4 ), parameter :: nf = 10
  real ( kind = 8 ) tol

  real ( kind = 8 ) am,av,hm,hv
  parameter (am = 0.0D+00, hm = 0.1D+00)
  parameter (av = -0.1D+00, hv = 0.2D+00)

  integer ( kind = 4 ) btl(3,maxbt),btst(maxfp+1),edno(maxfv)
  integer ( kind = 4 ) edst(maxfv),facep(3,maxfp),factyp(maxfp),fc(7,maxfc)
  integer ( kind = 4 ) fcst(maxfp+1),fvl(6,maxfv),hfl(maxhf),ht(maxht)
  integer ( kind = 4 ) iwk(maxiw),ntete(maxhf),ntetra(maxhf),ntrif(maxhf)
  integer ( kind = 4 ) pfl(2,maxpf),tetra(4,maxte),vm(maxvm),xfhv(3,maxhf+1)
  integer ( kind = 4 ) a,b,c,cnt,ceang,cfvl,chvl,cnrml,crit,d,e,f,htsiz,i,ifach
  integer ( kind = 4 ) imeas,in,ipolh,irdr,irem,j,k,kfc,kt,kvm,nbt,ncface,ncvert
  integer ( kind = 4 ) nface,nfach,nfhol,nfph,nihol,nit,nlo,npf,npolh,nrfe,nt
  integer ( kind = 4 ) nteta,ntetd,ntri,nvc,nvch,nvert,nverth,outmod,p,prime,utet
  real    dectim,holtim,initim,tritim,t0,t1
  real ( kind = 8 ) eang(maxfv),h(maxhf),nrml(3,maxfp),sa(4)
  real ( kind = 8 ) eta(maxnt),rho(maxnt),rvol(maxnt),sig(maxnt)
  real ( kind = 8 ) vcl(3,maxvc),vol(maxhf),wid(maxhf),wk(maxwk)
  real ( kind = 8 ) angacc,angmin,aspc2d,atol2d,emnrth,h3,radrth
  real ( kind = 8 ) radsq,rdacc,sangmn,sf,shrf,sumint,sumvol
  real ( kind = 8 ) volcph,volth,x
  real ( kind = 8 ) efreq(0:nf),emax,emean,emin,estdv
  real ( kind = 8 ) rfreq(0:nf),rmax,rmean,rmin,rstdv
  real ( kind = 8 ) sfreq(0:nf),smax,smean,smin,sstdv
  real ( kind = 8 ) vfreq(0:nf),vmax,vmean,vmin,vstdv
  logical hoflag,holint
  character rgname*60
  real ( kind = 8 ) temp
!
!     OUTMOD = 0 : Output nothing (except for measurements).
!     OUTMOD = 1 : Output polyh decomp, triangulation data structures.
!     OUTMOD = 2 : Output VCL and 3 vertex indices for each triangle,
!        boundary triangles (on faces of decomp) before interior ones.
!     OUTMOD = -1, -2: Same as + case, but also output tetrahedron list
!        to unit UTET.
!     NFHOL <= 0: Simpler routine DSPHDC is called instead of DSPHFH.
!
  irdr = 5
  imeas = 7
  utet = 1
  read (irdr,600) rgname
  read (irdr,*) tol, aspc2d,atol2d,angacc,rdacc,x,x,x,i,ntetd,crit, &
     shrf

  read (irdr,*) nvc,nface,nfhol,npolh,nihol,outmod,msglvl
  write (imeas,800) rgname,tol, aspc2d,atol2d,angacc,rdacc,nfhol, &
     nihol,ntetd,crit,shrf
  hoflag = (nfhol > 0)
  nfhol = max(nfhol,0)
  aspc2d = aspc2d* d_pi ( ) /180.0D+00
  atol2d = atol2d* d_pi ( ) /180.0D+00
  angacc = angacc* d_pi ( ) /180.0D+00
  nfph = nface + nfhol

  if (nvc > maxvc .or. nfph >= maxfp .or. npolh >= maxhf) then
    if (nvc > maxvc) write ( *, 610) 'maxvc',nvc
    if (nfph >= maxfp) write ( *, 610) 'maxfp',nfph+1
    if (npolh >= maxhf) write ( *, 610) 'maxhf',npolh+1
    stop
  end if

  read (irdr,*) ((vcl(j,i),j=1,3),i=1,nvc)
!
!     The outer polygons of NFACE faces should appear first, followed by
!     the inner polygons of NFHOL holes (in order of index of faces
!     containing holes). Face types for NFACE faces should be positive.
!     Face types for NFHOL holes are +F or -F where F is index of face
!     containing hole; positive (negative) sign indicates hole polygon
!     is oriented CCW (CW) in polyhedron when viewed from outside.
!     Holes must only be on boundary faces of polyhedral region.
!
  read (irdr,*) (facep(1,i),i=1,nfph+1)
  nvert = facep(1,nfph+1) - 1
  read (irdr,*) (factyp(i),i=1,nfph)

  if (nvert > maxfv) then
    write ( *, 610) 'maxfv',nvert
    stop
  end if
!
!     Head vertex of each (outer) face must be a strictly convex vertex.
!     Head vertex of each hole may be arbitrary.
!     Positive (negative) sign in PFL(1,*) indicates face is oriented
!     CCW (CW) in polyhedron when viewed from outside.
!     Hole polygons should not be included in PFL.
!
  read (irdr,*) (fvl(1,i),i=1,nvert)
  read (irdr,*) (hfl(i),i=1,npolh+1)
  npf = hfl(npolh+1) - 1

  if (npf + nfhol > maxpf) then
    write ( *, 610) 'maxpf',npf+nfhol
    return
  end if

  read (irdr,*) (pfl(1,i),i=1,npf)
  htsiz = min(prime(nvc+2),maxht)
  call gtime(t0)

  if (hoflag) then
   call dsphfh(aspc2d,atol2d,nvc,nface,nfhol,npolh,maxvc,maxfv, &
        maxiw,maxwk,nvert,npf,vcl,facep,factyp,nrml,fvl,eang,hfl, &
        pfl,htsiz,ht,iwk,wk,ierr)
  else
     call dsphdc(nvc,nface,npolh,vcl,facep,nrml,fvl,eang,hfl,pfl, &
        htsiz,maxiw/4,iwk,ht,ierr)
  end if

  call gtime(initim)
  initim = initim - t0

  if (ierr /= 0) then
    write ( *, 620) ierr
    write (imeas,620) ierr
    return
  end if

  nrfe = 0
  angmin = 2.0D+00 * d_pi ( )
  do i = 1,nvert
    if (eang(i) > d_pi ( ) + tol) nrfe = nrfe + 1
    if (eang(i) > -1.0D+00) angmin = min(angmin,eang(i))
  end do

  angmin = angmin*180.0D+00/ d_pi ( )

  if ( hoflag ) then
     write (imeas,810) 'fholds',nvc,nface,nvert,npolh,npf,nrfe, &
        angmin,'initim',initim
  else
     write (imeas,810) 'initds',nvc,nface,nvert,npolh,npf,nrfe, &
        angmin,'initim',initim
  end if
!
!     Read in description of interior hole polyhedra. The input format
!     is similar to above with no hole faces and only 1 polyhedron at
!     a time, treated as though the region for each polyh is its hole.
!     IPOLH is index of polyhedron containing hole. IFACH is index of
!     `extreme' face of hole used for connection to outer boundary.
!     One interior hole polyhedron per polyhedral region is assumed.
!     A negative sign for IPOLH indicates hole polyh is hole interface.
!
  holtim = 0.0

  do k = 1,nihol

    read (irdr,*) nvch,nfach,ipolh,ifach
    holint = (ipolh < 0)
    ipolh = abs(ipolh)

    if (nvc+nvch > maxvc .or. nface+nfach >= maxfp .or. npf+nfach > maxpf) then
      if (nvc+nvch > maxvc) write ( *, 610) 'maxvc',nvc+nvch
      if (nface+nfach >= maxfp) write ( *, 610) 'maxfp', nface+nfach+1
      if (npf+nfach > maxpf) write ( *, 610) 'maxpf',npf+nfach
      stop
    end if

    read (irdr,*) ((vcl(j,i),j=1,3),i=nvc+1,nvc+nvch)
    read (irdr,*) (facep(1,i),i=nface+1,nface+nfach+1)
    nverth = facep(1,nface+nfach+1) - 1
    read (irdr,*) (factyp(i),i=nface+1,nface+nfach)

    if (nvert+nverth > maxfv) then
      write ( *, 610) 'maxfv',nvert+nverth
      stop
    end if

    read (irdr,*) (fvl(1,i),i=nvert+1,nvert+nverth)
    read (irdr,*) (pfl(1,i),i=npf+1,npf+nfach)
    htsiz = min(prime(nvch+2),maxht)
    call gtime(t0)

    call dsphih(aspc2d,atol2d,angacc,rdacc,nvc,nface,nvert,npolh, &
      npf,nvch,nfach,ipolh,ifach,holint,maxvc,maxfp,maxfv,maxhf, &
      maxpf,maxiw,maxwk,vcl,facep,factyp,nrml,fvl,eang,hfl,pfl, &
      htsiz,ht,iwk,wk,ierr)

    call gtime(t1)
    holtim = holtim + (t1 - t0)

    if (ierr /= 0) then
      write ( *, 620) ierr
      write (imeas,620) ierr
      stop
    end if

  end do

  if (nihol > 0) then
     nrfe = 0
     angmin = 2.0D+00* d_pi ( )
     do i = 1,nvert
       if (eang(i) > d_pi ( ) + tol) nrfe = nrfe + 1
       if (eang(i) > -1.0D+00) angmin = min(angmin,eang(i))
     end do
     angmin = angmin*180.0D+00/ d_pi ( )
     write (imeas,810) 'iholds',nvc,nface,nvert,npolh,npf,nrfe, &
        angmin,'holtim',holtim
  end if
!
!  Decompose polyhedral region into convex parts.
!
  cnt = 1
   40 continue

  call gtime(t0)

  call cvdec3(angacc,rdacc,nvc,nface,nvert,npolh,npf,maxvc,maxfp, &
     maxfv,maxhf,maxpf,maxiw,maxwk,vcl,facep,factyp,nrml,fvl,eang, &
     hfl,pfl,iwk,wk,ierr)

  call gtime(dectim)
  dectim = dectim - t0
  if (ierr /= 0) then
   if (ierr /= 327 .or. cnt >= 3) then
      write ( *, 620) ierr
      write (imeas,620) ierr
      if (ierr /= 327) stop
   end if
  end if
  nrfe = 0
  angmin = 2.0D+00* d_pi ( )
  do i = 1,nvert
    if (eang(i) > d_pi ( ) + tol) nrfe = nrfe + 1
    if (eang(i) > -1.0D+00) angmin = min(angmin,eang(i))
  end do

  angmin = angmin*180.0D+00/ d_pi ( )

  write (imeas,820) nvc,nface,nvert,npolh,npf,nrfe,angmin,dectim
  if (ierr == 327 .and. cnt < 3) then
   angacc = angacc - d_pi ( ) /36.0D+00
   if (angacc > 0.0D+00) then
      rdacc = rdacc*0.95D+00
      ierr = 0
        cnt = cnt + 1
      go to 40
   end if
  else if (ierr == 327) then
   stop
  end if
!
!  Decompose faces of polyhedral region into convex subpolygons.
!
  call gtime(t0)
  call cvdecf(aspc2d,atol2d,nvc,nface,nvert,npf,maxvc,maxfp,maxfv, &
     maxpf,maxiw,maxwk,vcl,facep,factyp,nrml,fvl,eang,hfl,pfl,iwk, &
     wk,ierr)
  call gtime(dectim)
  dectim = dectim - t0
  if (ierr /= 0) then
   write ( *, 620) ierr
   write (imeas,620) ierr
   stop
  end if
  write (imeas,830) nvc,nface,nvert,npolh,npf,dectim
!
!     Determine mesh spacing in each convex polyhedron from volume,
!     width, and NTETD (number of desired tetrahedra in entire region).
!     Mesh distribution function in a polyhedron is 1/width**3.
!     NCFACE = max no. of faces in a polyh and NCVERT = 2 * max no. of
!     edges in a polyh. Array FCST is temporarily used.
!
  ncface = 0
  ncvert = 0
  do f = 1,nface
    k = 0
    i = facep(1,f)
   60    continue
      k = k + 1
      i = fvl(3,i)
    if (i /= facep(1,f)) go to 60
    fcst(f) = k
  end do

  do p = 1,npolh

    j = 0
    k = 0
    i = hfl(p)

80  continue
      f = abs(pfl(1,i))
      j = j + 1
      k = k + fcst(f)
      i = pfl(2,i)
    if (i /= hfl(p)) go to 80

    ncface = max(ncface,j)
    ncvert = max(ncvert,k)

  end do

  chvl = 1
  cfvl = chvl + ncface
  irem = cfvl + 5*ncvert
  cnrml = 1
  ceang = cnrml + 3*ncface

  if (irem + ncvert - 1 > maxiw) then
    write ( *, 610) 'maxiw',irem+ncvert-1
    stop
  else if (ceang + ncvert - 1 > maxwk) then
    write ( *, 610) 'maxwk',ceang+ncvert-1
    stop
  end if

  sumvol = 0.0D+00
  sumint = 0.0D+00

  do i = 1,npolh

    call dsconv(i,hfl(i),facep,nrml,fvl,eang,pfl,ncface,ncvert, &
        iwk(chvl),wk(cnrml),iwk(cfvl),wk(ceang))
    irem = cfvl + 5*ncvert
    call rmcpfc(ncface,ncvert,iwk(chvl),wk(cnrml),iwk(cfvl), &
        wk(ceang),iwk(irem))
    call rmcled(ncface,ncvert,iwk(chvl),iwk(cfvl))
    vol(i) = volcph(ncface,vcl,iwk(chvl),iwk(cfvl))
    irem = cfvl + 5*ncvert
    call width3(ncface,vcl,iwk(chvl),wk(cnrml),iwk(cfvl), &
        maxiw-irem+1,j,k,wid(i),iwk(irem), ierr )
    if (ierr /= 0) then
      write ( *, 620) ierr
      write (imeas,620) ierr
      return
    end if

    sumvol = sumvol + vol(i)
    sumint = sumint + vol(i)/wid(i)**3

  end do

  do i = 1,npolh
    temp = dble ( ntetd ) * vol(i)/wid(i)**3/sumint + 0.5D+00
    ntete(i) = max ( int ( temp ), 1 )
    h(i) = (6.0D+00*vol(i)/ntete(i))**(1.0D+00/3.0D+00)
  end do

  call gtime(t0)
  call tripr3(h,shrf,crit,nvc,nface,nvert,npolh,maxvc,maxbt, &
     maxfc,maxht,maxvm,maxiw,maxwk,vcl,facep,nrml,fvl,eang,hfl,pfl, &
     edst,edno,fcst,btst,btl,fc,ht,vm,xfhv,ntrif,ntetra,iwk,wk, ierr )
  call gtime(tritim)
  tritim = tritim - t0

  if (ierr /= 0) then
    write ( *, 620) ierr
    write (imeas,620) ierr
    return
  end if

  nbt = btst(nface+1)-1

  ntri = 0
  nteta = 0
  nt = 0

  do i = 1,npolh
    ntri = ntri + ntrif(i)
    nteta = nteta + ntetra(i)
    nt = max(nt,ntetra(i))
  end do

  if (npolh > maxiw) then
    write ( *, 610) 'maxiw',npolh
    return
  else if (nt > maxte) then
    write ( *, 610) 'maxte',nt
    return
  else if (nteta > maxnt) then
    write ( *, 610) 'maxnt',nteta
    return
  end if

  if (outmod < 0) then
    write (utet,730) nvc,nteta
    write (utet,740) ((vcl(j,i),j=1,3),i=1,nvc)
    write (utet,730)
  end if
!
!     Collect measurements for number of interior faces that fail local
!     sphere test, and min, max, mean, stdv of relative volume, solid
!     angle, norm. ratio of inradius to circumradius for all tetrahedra.
!
  nit = 0
  nlo = 0
  kt = 0
  sf = 1.5D+00*sqrt(6.0D+00)
  do 160 p = 1,npolh
   k = 0
   j = hfl(p)
  130    continue
      f = abs(pfl(1,j))
      k = k + (btst(f+1) - btst(f))
      j = pfl(2,j)
   if (j /= hfl(p)) go to 130
   nit = nit + (ntrif(p) - k)
   kfc = xfhv(1,p) + k
   iwk(p) = kfc
   kvm = xfhv(3,p) - 1

   do i = kfc,xfhv(1,p+1)-1
      if (fc(1,i) > 0) then
         a = vm(fc(1,i)+kvm)
         b = vm(fc(2,i)+kvm)
         c = vm(fc(3,i)+kvm)
         d = vm(fc(4,i)+kvm)
         e = vm(fc(5,i)+kvm)
         call ccsph(.true.,vcl(1,a),vcl(1,b),vcl(1,c),vcl(1,d), &
              vcl(1,e),sa,radsq,in)
         if (in >= 1) nlo = nlo + 1
      end if
  end do

   call tetlst(xfhv(1,p+1)-xfhv(1,p),vm(xfhv(3,p)), &
        fc(1,xfhv(1,p)),nt,tetra)
   if (nt /= ntetra(p)) then
      write ( *, 750) p,nt,ntetra(p)
      write (imeas,750) p,nt,ntetra(p)
      stop
   end if
   if (outmod < 0) write (utet,730) ((tetra(j,i),j=1,4),i=1,nt)
   h3 = 1.0D+00/h(p)**3
   do i = 1,nt
      kt = kt + 1
      a = tetra(1,i)
      b = tetra(2,i)
      c = tetra(3,i)
      d = tetra(4,i)
      rvol(kt) = h3*volth(vcl(1,a),vcl(1,b),vcl(1,c),vcl(1,d))
      rho(kt) = radrth(vcl(1,a),vcl(1,b),vcl(1,c),vcl(1,d))
      sig(kt) = sangmn(vcl(1,a),vcl(1,b),vcl(1,c),vcl(1,d),sa)*sf
      eta(kt) = emnrth(vcl(1,a),vcl(1,b),vcl(1,c),vcl(1,d))
  end do
  160 continue
  write (imeas,840) fcst(1)-1,fcst(nface+1)-1,nvc,nbt, &
     (xfhv(i,npolh+1)-1,i=1,3),ntri,nteta,nlo,tritim

  call stats(nteta,rvol,av,hv,nf,vmin,vmax,vmean,vstdv,vfreq)
  call stats(nteta,rho,am,hm,nf,rmin,rmax,rmean,rstdv,rfreq)
  call stats(nteta,sig,am,hm,nf,smin,smax,smean,sstdv,sfreq)
  call stats(nteta,eta,am,hm,nf,emin,emax,emean,estdv,efreq)
  write (imeas,850) vmin,vmax,vmean,vstdv
  write (imeas,860) rmin,rmax,rmean,rstdv
  write (imeas,870) smin,smax,smean,sstdv
  write (imeas,880) emin,emax,emean,estdv
  write (imeas,890)

  do i = 0,nf
    write (imeas,900) av+i*hv,vfreq(i),am+i*hm,rfreq(i),am+i*hm, &
        sfreq(i),am+i*hm,efreq(i)
  end do

  if (abs(outmod) == 1) then
     write ( *, 630) nvc,nface,nvert,npolh,npf
     write ( *, 640) (i,(vcl(j,i),j=1,3),i=1,nvc)
     write ( *, 650) (i,(facep(j,i),j=1,3),factyp(i), &
        (nrml(j,i),j=1,3),fcst(i),btst(i),i=1,nface), nface+1, &
        0,0,0,0,0.0D+00,0.0D+00,0.0D+00,fcst(nface+1),btst(nface+1)
     write ( *, 660) (i,(fvl(j,i),j=1,6),eang(i),edst(i),edno(i), &
        i=1,nvert)
     write ( *, 670) (i,hfl(i),(xfhv(j,i),j=1,3),ntrif(i), &
        ntetra(i),ntete(i),h(i),vol(i),wid(i),i=1,npolh), &
        npolh+1,0,(xfhv(j,npolh+1),j=1,3)
     write ( *, 680) (i,(pfl(j,i),j=1,2),i=1,npf)
     write ( *, 690) (i,(btl(j,i),j=1,3),i=1,nbt)
     write ( *, 700) (i,(fc(j,i),j=1,7),i=1,xfhv(1,npolh+1)-1)
     write ( *, 710) (ht(i),i=1,xfhv(2,npolh+1)-1)
     write ( *, 720) (vm(i),i=1,xfhv(3,npolh+1)-1)
  else if (abs(outmod) == 2) then
    write ( *, 730) nvc,nbt,nit
    write ( *, 740) ((vcl(j,i),j=1,3),i=1,nvc)
    write ( *, '(a)' )
    write ( *, 730) (3,(btl(j,i),j=1,3),i=1,nbt)
    write ( *, '(a)' )
    do i = 1,npolh
      kvm = xfhv(3,i) - 1
      do j = iwk(i),xfhv(1,i+1)-1
         if (fc(1,j) > 0) &
              write ( *, 730) 3,(vm(fc(k,j)+kvm),k=1,3)
      end do
    end do
  end if

  600 format (a60)
  610 format (1x,'*** ',a,' must be increased to',i8)
  620 format (/1x,'ierr=',i5)
  630 format (1x,'nvc,nface,nvert,npolh,npf'/1x,5i7)
  640 format (/' vcl'/(1x,i7,3f23.15))
  650 format (/' facep,factyp,nrml,fcst,btst'/(1x,5i5,3f14.6,2i5))
  660 format (/' fvl,eang,edst,edno'/(1x,7i5,f15.7,2i7))
  670 format (/' hfl,xfhv,ntrif,ntetra,ntete,h,vol,wid'/(1x,8i5,3d12.5))
  680 format (/' pfl'/(1x,3i7))
  690 format (/' btl'/(1x,4i7))
  700 format (/' fc'/(1x,8i7))
  710 format (/' ht'/(1x,10i7))
  720 format (/' vm'/(1x,10i7))
  730 format (1x,4i7)
  740 format (/(1x,3f15.7))
  750 format (/1x,'error from tetlst polyh',i5,',   nt != ntetra :',2i7)
  800 format (1x,a60/1x,'input : tol=',d15.7,'   aspc2d=',f9.3, &
     '   atol2d=',f9.3/9x,'angacc=',f9.3,'   rdacc=',f9.5, &
     '   nfhol=',i3,'   nihol=',i3/9x,'ntetd=',i7,'   crit=',i2, &
     '   shrf=',f9.5)
  810 format (1x,a6,': nvc=',i7,'   nface=',i7,'   nvert=',i7, &
     '   npolh=',i7/9x,'npf=',i7,'   nrfe=',i7,'   angmin=',f9.3, &
     3x,a6,'=',f9.3)
  820 format (1x,'decomp: nvc=',i7,'   nface=',i7,'   nvert=',i7, &
     '   npolh=',i7/9x,'npf=',i7,'   nrfe=',i7,'   angmin=',f9.3, &
     '   dectim=',f9.3)
  830 format (1x,'decfac: nvc=',i7,'   nface=',i7,'   nvert=',i7, &
     '   npolh=',i7/9x,'npf=',i7,'   dectim=',f9.3)
  840 format (1x,'triang: nvce=',i7,'   nvcf=',i7,'   nvc=',i7, &
     '   nbt=',i7/9x,'nfc=',i7,'   nht=',i7,'   nvm=',i7,'   ntri=', &
     i7/9x,'nteta=',i7,'   nlo=',i7,'   tritim=',f9.3)
  850 format (1x,'vmin=',f11.7,3x,'vmax=',f11.7,3x,'vmean=',f11.7, &
     3x,'vstdv=',f11.7)
  860 format (1x,'rmin=',f11.7,3x,'rmax=',f11.7,3x,'rmean=',f11.7, &
     3x,'rstdv=',f11.7)
  870 format (1x,'smin=',f11.7,3x,'smax=',f11.7,3x,'smean=',f11.7, &
     3x,'sstdv=',f11.7)
  880 format (1x,'emin=',f11.7,3x,'emax=',f11.7,3x,'emean=',f11.7, &
     3x,'estdv=',f11.7)
  890 format (/3x,'rel volume',9x,'radius ratio',8x,'min solid ang',8x, &
     'mean ratio')
  900 format (1x,f5.2,f9.4,6x,f5.2,f9.4,6x,f5.2,f9.4,6x,f5.2,f9.4)

  return
end
subroutine test_dec3 ( )

!*****************************************************************************80
!
!! TEST_DEC3 tests CVDEC3, RESEDG, CUTFAC, CVDECF.
!
  implicit none

  real ( kind = 8 ) d_pi
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) msglvl
  integer ( kind = 4 ) outmod
  integer ( kind = 4 ) prime
  real ( kind = 8 ) tol

  integer ( kind = 4 ) maxfp,maxfv,maxhf,maxht,maxiw,maxpf,maxvc,maxwk
  parameter (maxfp = 800)
  parameter (maxfv = 3500)
  parameter (maxhf = 200)
  parameter (maxht = 307)
  parameter (maxiw = 5000)
  parameter (maxpf = 2000)
  parameter (maxvc = 600)
  parameter (maxwk = 5000)

  integer ( kind = 4 ) facep(3,maxfp),factyp(maxfp),fvl(6,maxfv),hfl(maxhf)
  integer ( kind = 4 ) ht(0:maxht-1),iwk(maxiw),pfl(2,maxpf)
  integer ( kind = 4 ) cnt,f,htsiz,i,ifach,imeas,ipolh,irdr,j,k,nbf,nface,nfach
  integer ( kind = 4 ) nfhol,nfph,nihol,npf,npolh,nrfe,nvc,nvch,nvert,nverth
  real    dectim,holtim,initim,t0,t1
  real ( kind = 8 ) eang(maxfv),nrml(3,maxfp),vcl(3,maxvc),wk(maxwk)
  real ( kind = 8 ) angacc,angmin,aspc2d,atol2d,rdacc
  logical hoflag,holint
  character rgname*60
!
!     OUTMOD = 1 : Output polyhedral decomposition data structure.
!     OUTMOD = 2 : Output VCL + list of vertex indices for each face.
!     MSGLVL >= 2: Routine DSPHFH prints out non-coplanar vertices of
!        each face (also output from routines CUTFAC, INSFAC).
!     NFHOL < 0: Simpler routine DSPHDC is called instead of DSPHFH.
!
  irdr = 5
  imeas = 7
  read (irdr,600) rgname
  read (irdr,*) tol, aspc2d,atol2d,angacc,rdacc

  read (irdr,*) nvc,nface,nfhol,npolh,nihol,outmod,msglvl
  write (imeas,750) rgname,tol, aspc2d,atol2d,angacc,rdacc,nfhol, &
     nihol
  hoflag = (nfhol >= 0)
  nfhol = max(nfhol,0)
  aspc2d = aspc2d* d_pi ( ) /180.0D+00
  atol2d = atol2d* d_pi ( ) /180.0D+00
  angacc = angacc* d_pi ( ) /180.0D+00
  nfph = nface + nfhol
  if (nvc > maxvc .or. nfph >= maxfp .or. npolh >= maxhf) &
  then
   if (nvc > maxvc) write ( *, 610) 'maxvc',nvc
   if (nfph >= maxfp) write ( *, 610) 'maxfp',nfph+1
   if (npolh >= maxhf) write ( *, 610) 'maxhf',npolh+1
   stop
  end if
  read (irdr,*) ((vcl(j,i),j=1,3),i=1,nvc)
!
!  The outer polygons of NFACE faces should appear first, followed by
!  the inner polygons of NFHOL holes (in order of index of faces
!  containing holes). Face types for NFACE faces should be positive.
!  Face types for NFHOL holes are +F or -F where F is index of face
!  containing hole; positive (negative) sign indicates hole polygon
!  is oriented CCW (CW) in polyhedron when viewed from outside.
!  Holes must only be on boundary faces of polyhedral region.
!
  read (irdr,*) (facep(1,i),i=1,nfph+1)
  nvert = facep(1,nfph+1) - 1
  read (irdr,*) (factyp(i),i=1,nfph)

  if ( maxfv < nvert ) then
    write ( *, 610) 'maxfv',nvert
    stop
  end if
!
!     Head vertex of each (outer) face must be a strictly convex vertex.
!     Head vertex of each hole may be arbitrary.
!     Positive (negative) sign in PFL(1,*) indicates face is oriented
!     CCW (CW) in polyhedron when viewed from outside.
!     Hole polygons should not be included in PFL.
!
  read (irdr,*) (fvl(1,i),i=1,nvert)
  read (irdr,*) (hfl(i),i=1,npolh+1)
  npf = hfl(npolh+1) - 1
  if (npf + nfhol > maxpf) then
   write ( *, 610) 'maxpf',npf+nfhol
   stop
  end if
  read (irdr,*) (pfl(1,i),i=1,npf)
  htsiz = min(prime(nvc+2),maxht)
  call gtime(t0)

  if (hoflag) then
   call dsphfh(aspc2d,atol2d,nvc,nface,nfhol,npolh,maxvc,maxfv, &
        maxiw,maxwk,nvert,npf,vcl,facep,factyp,nrml,fvl,eang,hfl, &
        pfl,htsiz,ht,iwk,wk,ierr)
  else
     call dsphdc(nvc,nface,npolh,vcl,facep,nrml,fvl,eang,hfl,pfl, &
        htsiz,maxiw/4,iwk,ht,ierr)
  end if

  call gtime(initim)
  initim = initim - t0
  if (ierr /= 0) then
   write ( *, 620) ierr
   write (imeas,620) ierr
   stop
  end if
  nrfe = 0
  angmin = 2.0D+00* d_pi ( )

  do i = 1,nvert
    if (eang(i) > d_pi ( ) + tol) nrfe = nrfe + 1
    if (eang(i) > -1.0D+00) angmin = min(angmin,eang(i))
  end do

  angmin = angmin*180.0D+00/ d_pi ( )
  if (hoflag) then
     write (imeas,760) 'fholds',nvc,nface,nvert,npolh,npf,nrfe, &
        angmin,'initim',initim
  else
     write (imeas,760) 'initds',nvc,nface,nvert,npolh,npf,nrfe, &
        angmin,'initim',initim
  end if
!
!     Read in description of interior hole polyhedra. The input format
!     is similar to above with no hole faces and only 1 polyhedron at
!     a time, treated as though the region for each polyh is its hole.
!     IPOLH is index of polyhedron containing hole. IFACH is index of
!     `extreme' face of hole used for connection to outer boundary.
!     One interior hole polyhedron per polyhedral region is assumed.
!     A negative sign for IPOLH indicates hole polyh is hole interface.
!
  holtim = 0.0
  do 20 k = 1,nihol
     read (irdr,*) nvch,nfach,ipolh,ifach
   holint = (ipolh < 0)
   ipolh = abs(ipolh)
     if (nvc+nvch > maxvc .or. nface+nfach >= maxfp .or. &
     npf+nfach > maxpf) then
      if (nvc+nvch > maxvc) write ( *, 610) 'maxvc',nvc+nvch
      if (nface+nfach >= maxfp) write ( *, 610) 'maxfp', &
           nface+nfach+1
      if (npf+nfach > maxpf) write ( *, 610) 'maxpf',npf+nfach
      stop
     end if
     read (irdr,*) ((vcl(j,i),j=1,3),i=nvc+1,nvc+nvch)
     read (irdr,*) (facep(1,i),i=nface+1,nface+nfach+1)
     nverth = facep(1,nface+nfach+1) - 1
     read (irdr,*) (factyp(i),i=nface+1,nface+nfach)
     if (nvert+nverth > maxfv) then
      write ( *, 610) 'maxfv',nvert+nverth
      stop
     end if
   read (irdr,*) (fvl(1,i),i=nvert+1,nvert+nverth)
     read (irdr,*) (pfl(1,i),i=npf+1,npf+nfach)
     htsiz = min(prime(nvch+2),maxht)
     call gtime(t0)
     call dsphih(aspc2d,atol2d,angacc,rdacc,nvc,nface,nvert,npolh, &
        npf,nvch,nfach,ipolh,ifach,holint,maxvc,maxfp,maxfv,maxhf, &
        maxpf,maxiw,maxwk,vcl,facep,factyp,nrml,fvl,eang,hfl,pfl, &
        htsiz,ht,iwk,wk,ierr)

     call gtime(t1)
     holtim = holtim + (t1 - t0)
     if (ierr /= 0) then
      write ( *, 620) ierr
      write (imeas,620) ierr
      stop
     end if
   20 continue
  if (nihol > 0) then
     nrfe = 0
     angmin = 2.0D+00* d_pi ( )
     do 30 i = 1,nvert
      if (eang(i) > d_pi ( ) + tol) nrfe = nrfe + 1
      if (eang(i) > -1.0D+00) angmin = min(angmin,eang(i))
   30    continue
     angmin = angmin*180.0D+00/ d_pi ( )
     write (imeas,760) 'iholds',nvc,nface,nvert,npolh,npf,nrfe, &
        angmin,'holtim',holtim
  end if
!
!     Decompose polyhedral region into convex parts.
!
  cnt = 1
   40 continue
  call gtime(t0)
  call cvdec3(angacc,rdacc,nvc,nface,nvert,npolh,npf,maxvc,maxfp, &
     maxfv,maxhf,maxpf,maxiw,maxwk,vcl,facep,factyp,nrml,fvl,eang, &
     hfl,pfl,iwk,wk,ierr)
  call gtime(dectim)
  dectim = dectim - t0
  if (ierr /= 0) then
   if (ierr /= 327 .or. cnt >= 3) then
      write ( *, 620) ierr
      write (imeas,620) ierr
      if (ierr /= 327) stop
   end if
  end if
  nrfe = 0
  angmin = 2.0D+00* d_pi ( )
  do i = 1,nvert
   if (eang(i) > d_pi ( ) + tol) nrfe = nrfe + 1
   if (eang(i) > -1.0D+00) angmin = min(angmin,eang(i))
  end do
  angmin = angmin*180.0D+00/ d_pi ( )
  write (imeas,770) nvc,nface,nvert,npolh,npf,nrfe,angmin,dectim
  if (ierr == 327 .and. cnt < 3) then
   angacc = angacc - d_pi ( ) /36.0D+00
   if (angacc > 0.0D+00) then
      rdacc = rdacc*0.95D+00
      ierr = 0
        cnt = cnt + 1
      go to 40
   end if
  end if
  if (ierr == 327) stop
!
!  Decompose faces of polyhedral region into convex subpolygons.
!
  call gtime(t0)
  call cvdecf(aspc2d,atol2d,nvc,nface,nvert,npf,maxvc,maxfp,maxfv, &
     maxpf,maxiw,maxwk,vcl,facep,factyp,nrml,fvl,eang,hfl,pfl,iwk, &
     wk,ierr)
  call gtime(dectim)
  dectim = dectim - t0
  if (ierr /= 0) then
   write ( *, 620) ierr
   write (imeas,620) ierr
   stop
  end if
  write (imeas,780) nvc,nface,nvert,npolh,npf,dectim

  if (outmod == 1) then
     write ( *, 630) nvc,nface,nvert,npolh,npf
     write ( *, 640) (i,(vcl(j,i),j=1,3),i=1,nvc)
     write ( *, 650) (i,(facep(j,i),j=1,3),factyp(i), &
        (nrml(j,i),j=1,3),i=1,nface)
     write ( *, 660) (i,(fvl(j,i),j=1,6),eang(i),i=1,nvert)
     write ( *, 670) (hfl(i),i=1,npolh)
     write ( *, 680) (i,(pfl(j,i),j=1,2),i=1,npf)
  else
   nbf = 0
   do f = 1,nface
      if (facep(3,f) == 0) nbf = nbf + 1
   end do
   write ( *, 690) nvc,nbf,nface-nbf
   write ( *, 700) ((vcl(j,i),j=1,3),i=1,nvc)
   do 90 i = 1,2
        write ( *, '(a)' )
      do 80 f = 1,nface
         if (i == 1 .and. facep(3,f) /= 0) go to 80
         if (i == 2 .and. facep(3,f) == 0) go to 80
         k = 0
         j = facep(1,f)
   70          continue
      k = k + 1
      iwk(k) = fvl(1,j)
      j = fvl(3,j)
         if (j /= facep(1,f)) go to 70
         if (k > maxiw) then
            write ( *, 610) 'maxiw',k
            stop
         end if
         write ( *, 710) k,(iwk(j),j=1,k)
   80       continue
   90    continue
  end if

  600 format (a60)
  610 format (1x,'*** ',a,' must be increased to',i8)
  620 format (/1x,'ierr=',i5)
  630 format (1x,'nvc,nface,nvert,npolh,npf'/1x,5i7)
  640 format (/' vcl'/(1x,i7,3f15.7))
  650 format (/' facep,factyp,nrml'/(1x,4i7,i5,3f15.7))
  660 format (/' fvl,eang'/(1x,7i7,f15.7))
  670 format (/' hfl'/(1x,10i7))
  680 format (/' pfl'/(1x,3i7))
  690 format (1x,3i7)
  700 format (/(1x,3f15.7))
  710 format (1x,11i7/(8x,10i7))
  750 format (1x,a60/1x,'input : tol=',d15.7,'   aspc2d=',f9.3, &
     '   atol2d=',f9.3/9x,'angacc=',f9.3,'   rdacc=',f9.5, &
     '   nfhol=',i3,'   nihol=',i3)
  760 format (1x,a6,': nvc=',i7,'   nface=',i7,'   nvert=',i7, &
     '   npolh=',i7/9x,'npf=',i7,'   nrfe=',i7,'   angmin=',f9.3, &
     3x,a6,'=',f9.3)
  770 format (1x,'decomp: nvc=',i7,'   nface=',i7,'   nvert=',i7, &
     '   npolh=',i7/9x,'npf=',i7,'   nrfe=',i7,'   angmin=',f9.3, &
     '   dectim=',f9.3)
  780 format (1x,'decfac: nvc=',i7,'   nface=',i7,'   nvert=',i7, &
     '   npolh=',i7/9x,'npf=',i7,'   dectim=',f9.3)

  return
end
subroutine test_ds ( )

!*****************************************************************************80
!
!! TEST_DS tests DSMCPR, DSPGDC, EDGHT, HOLVRT.
!
  implicit none

  integer ( kind = 4 ) ierr

  integer ( kind = 4 ) incr,maxed,maxho,maxnc,maxnv
  parameter (incr = 10000)
  parameter (maxed = 101)
  parameter (maxho = 50)
  parameter (maxnc = 30)
  parameter (maxnv = 500)

  integer ( kind = 4 ) edge(4,maxed)
  integer ( kind = 4 ) holv(maxho)
  integer ( kind = 4 ) ht(0:maxed-1)
  integer ( kind = 4 ) hvl(maxnc*2)
  integer ( kind = 4 ) icur(maxnc)
  integer ( kind = 4 ) ivrt(maxnv)
  integer ( kind = 4 ) map(maxnc)
  integer ( kind = 4 ) nvbc(maxnc),pvl(4,maxnv*2),regnum(maxnc*2)
  integer ( kind = 4 ) case,htsiz,i,irdr,j,ncur,nhola,nhole,npolg,nsc
  integer ( kind = 4 ) nv,nvc,nvert,prime
  real ( kind = 8 ) iang(maxnv*2),vcl(2,maxnv)
  real ( kind = 8 ) tol
  character rgname*20
!
!     Read in vertices of general polygonal region.
!     CASE = 1 : simple polygon or multiply connected polygonal region
!     CASE = 2 : general polygonal region with holes and interfaces
!
  irdr = 5
  read (irdr,600) rgname
  read (irdr,*) tol
  read (irdr,*) case,nvc,ncur

  if ( maxnv < nvc ) then
    write ( *, 610) 'maxnv',nvc
    stop
  end if

  if (ncur > maxnc) then
    write ( *, 610) 'maxnc',ncur
    stop
  end if

  read (irdr,*) (nvbc(i),i=1,ncur)

  if (case == 2) then
    read (irdr,*) (icur(i),i=1,ncur)
  end if

  read (irdr,*) (vcl(1,i),vcl(2,i),i=1,nvc)
!
!  Call routine DSMCPR or DSPGDC to set data structures in arrays
!  REGNUM, HVL, PVL, IANG, HOLV.
!
  if ( case == 1 ) then

    nhole = ncur - 1

    call dsmcpr(nhole,nvbc,vcl,maxnc*2,maxnv*2,maxho,nvc,npolg, &
        nvert,nhola,regnum,hvl,pvl,iang,holv,ierr)

  else if (case == 2) then

   nv = sum ( nvbc(1:ncur) )

   if (nv > maxnv) then
      write ( *, 610) 'maxnv',nv
      stop
   end if

   read (irdr,*) (ivrt(i),i=1,nv)
   nsc = 0
   do i = 1,nv
     if (ivrt(i) < 0) nsc = nsc + 1
   end do
   nsc = nsc/2
   if (nsc > maxed) then
      write ( *, 610) 'maxed',nsc
      stop
   end if
   htsiz = min(prime(nsc/2),maxed)

   call dspgdc(nvc,vcl,incr,ncur,nvbc,icur,ivrt,maxnc*2,maxnv*2, &
        maxho,npolg,nvert,nhole,nhola,regnum,hvl,pvl,iang,holv, &
        htsiz,nsc,ht,edge,map,ierr)

  end if

  if (ierr /= 0) then
    write ( *, 620) ierr
    return
  end if
!
!  Print out data structures.
!
  write ( *, 630) nvc,(i,vcl(1,i),vcl(2,i),i=1,nvc)
  write ( *, 640) npolg,nhole,(hvl(i),i=1,npolg+nhole)
  write ( *, 650) (regnum(i),i=1,npolg)
  write ( *, 660) nvert,(i,(pvl(j,i),j=1,4),iang(i),i=1,nvert)
  write ( *, 640) nhola,nhole*2+nhola,(holv(i),i=1,nhole*2+nhola)

  600 format (a20)
  610 format (1x,'*** ',a,' must be increased to',i8)
  620 format (/1x,'ierr=',i5)
  630 format (/1x,i7/(1x,i7,2f15.7))
  640 format (/1x,2i7/(1x,10i7))
  650 format (/(1x,10i7))
  660 format (/1x,i7/(1x,5i7,f15.7))

  return
end
subroutine test_ds3 ( )

!*****************************************************************************80
!
!! TEST_DS3 tests DSPHDC,INSVR3,INSED3,INSFAC.
!
  implicit none

  real ( kind = 8 ) d_pi
  integer ( kind = 4 ) ierr
  real ( kind = 8 ) tol

  integer ( kind = 4 ) maxed,maxfp,maxfv,maxhf,maxht,maxpf,maxvc,maxvf
  parameter (maxed = 200)
  parameter (maxfp = 200)
  parameter (maxfv = 1000)
  parameter (maxhf = 100)
  parameter (maxht = 307)
  parameter (maxpf = 200)
  parameter (maxvc = 300)
  parameter (maxvf = 25)

  integer ( kind = 4 ) cedge(2,0:maxvf),edge(4,maxed),facep(3,maxfp)
  integer ( kind = 4 ) factyp(maxfp),fvl(6,maxfv),hfl(maxhf),ht(0:maxht-1)
  integer ( kind = 4 ) pfl(2,maxpf),vindx(maxvf)
  integer ( kind = 4 ) a,b,f,htsiz,i,irdr,j,k,nbf,nce,nface,npf,npolh
  integer ( kind = 4 ) nvc,nvert,outmod,p,prime
  real ( kind = 8 ) eang(maxfv),nrml(3,maxfp),vcl(3,maxvc)
  real ( kind = 8 ) cdang(maxvf),nrmlc(3),sf
  character ch*1,rgname*20
!
!  OUTMOD = 1 : Output polyhedral decomposition data structure.
!  OUTMOD = 2 : Output VCL + list of vertex indices for each face.
!
  irdr = 5
  read (irdr,600) rgname
  read (irdr,*) tol
  read (irdr,*) nvc,nface,i,npolh,i,outmod

  if (nvc > maxvc .or. nface >= maxfp .or. npolh >= maxhf) &
  then
   if (nvc > maxvc) write ( *, 610) 'maxvc',nvc
   if (nface >= maxfp) write ( *, 610) 'maxfp',nface+1
   if (npolh >= maxhf) write ( *, 610) 'maxhf',npolh+1
   stop
  end if
  read (irdr,*) ((vcl(j,i),j=1,3),i=1,nvc)
  read (irdr,*) (facep(1,i),i=1,nface+1)
  nvert = facep(1,nface+1) - 1
!
!  Face types should be positive.
!
  read (irdr,*) (factyp(i),i=1,nface)
  if (nvert > maxfv) then
    write ( *, 610) 'maxfv',nvert
    return
  end if
!
!  Head vertex of each face must be a strictly convex vertex.
!
  read (irdr,*) (fvl(1,i),i=1,nvert)
  read (irdr,*) (hfl(i),i=1,npolh+1)
  npf = hfl(npolh+1) - 1
  if (npf > maxpf) then
    write ( *, 610) 'maxpf',npf
    return
  end if
!
!  Positive (negative) sign in PFL(1,*) indicates face is oriented
!  CCW (CW) in polyhedron when viewed from outside.
!
  read (irdr,*) (pfl(1,i),i=1,npf)
  htsiz = min(prime(nvc+2),maxht)
  call dsphdc(nvc,nface,npolh,vcl,facep,nrml,fvl,eang,hfl,pfl, &
     htsiz,maxed,edge,ht,ierr)

  if (ierr /= 0) then
    write ( *, 620) ierr
    return
  end if
!
!  Read commands to update polyhedral decomposition data structure.
!  CH = 'V': Insert vertex V on edge AB where A is index of FVL and
!               B = FVL(SUCC,A).
!     CH = 'E': Add edge from A to B on a face; A, B are indices of FVL.
!     CH = 'F': Create new face with unit normal NRMLC in polyhedron P
!               with info from CEDGE and CDANG arrays (latter in deg).
!
  sf = d_pi ( ) /180.0D+00

10 continue

  read (irdr,*,end=30) ch

  if (ch == 'v') then
    read (irdr,*) a,(vcl(j,nvc+1),j=1,3)
    call insvr3(a,nvc,nvert,maxfv,vcl,fvl,eang,ierr)
  else if (ch == 'e') then
    read (irdr,*) a,b
    call insed3(a,b,nface,nvert,npf,maxfp,maxfv,maxpf,facep, &
        factyp,nrml,fvl,eang,hfl,pfl,ierr)
  else if (ch == 'f') then
    read (irdr,*) p,nce,(nrmlc(j),j=1,3)
    if (nce > maxvf) then
      write ( *, 610) 'maxvf',nce
      stop
    end if
    k = 0

    do i = 1,nce
      read (irdr,*) cedge(1,i),cedge(2,i),cdang(i)
      cdang(i) = cdang(i)*sf
      if (cedge(1,i) > nvc) then
         k = k + 1
         read (irdr,*) (vcl(j,nvc+k),j=1,3)
      end if
    end do

   if (nvc+k > maxvc) then
      write ( *, 610) 'maxvc',nvc+k
      stop
    end if
    cedge(1,0) = cedge(1,nce)
    call insfac(p,nrmlc,nce,cedge,cdang,nvc,nface,nvert,npolh,npf, &
        maxfp,maxfv,maxhf,maxpf,vcl,facep,factyp,nrml,fvl,eang,hfl, &
        pfl,ierr)
  else
    go to 30
  end if

  if (ierr /= 0) then
   write ( *, 620) ierr
   stop
  end if
  go to 10

   30 continue

  if (outmod == 1) then
     write ( *, 630) nvc,nface,nvert,npolh,npf
     write ( *, 640) (i,(vcl(j,i),j=1,3),i=1,nvc)
     write ( *, 650) (i,(facep(j,i),j=1,3),factyp(i), &
        (nrml(j,i),j=1,3),i=1,nface)
     write ( *, 660) (i,(fvl(j,i),j=1,6),eang(i),i=1,nvert)
     write ( *, 670) (hfl(i),i=1,npolh)
     write ( *, 680) (i,(pfl(j,i),j=1,2),i=1,npf)
  else
   nbf = 0
   do f = 1,nface
      if (facep(3,f) == 0) nbf = nbf + 1
   end do
   write ( *, 690) nvc,nbf,nface-nbf
   write ( *, 700) ((vcl(j,i),j=1,3),i=1,nvc)
   do 70 i = 1,2
        write ( *, '(a)' )
      do 60 f = 1,nface
         if (i == 1 .and. facep(3,f) /= 0) go to 60
         if (i == 2 .and. facep(3,f) == 0) go to 60
         k = 0
         j = facep(1,f)
   50          continue
      k = k + 1
      vindx(k) = fvl(1,j)
      j = fvl(3,j)
         if (j /= facep(1,f)) go to 50
         if (k > maxvf) then
            write ( *, 610) 'maxvf',k
            stop
         end if
         write ( *, 710) k,(vindx(j),j=1,k)
   60       continue
   70    continue
  end if

  600 format (a20)
  610 format (1x,'*** ',a,' must be increased to',i8)
  620 format (/1x,'ierr=',i5)
  630 format (1x,'nvc,nface,nvert,npolh,npf'/1x,5i7)
  640 format (/' vcl'/(1x,i7,3f15.7))
  650 format (/' facep,factyp,nrml'/(1x,4i7,i5,3f15.7))
  660 format (/' fvl,eang'/(1x,7i7,f15.7))
  670 format (/' hfl'/(1x,10i7))
  680 format (/' pfl'/(1x,3i7))
  690 format (1x,3i7)
  700 format (/(1x,3f15.7))
  710 format (1x,11i7/(8x,10i7))

  return
end
subroutine test_dsh3 ( )

!*****************************************************************************80
!
!! TEST_DSH3 tests DSPHFH,SPDECH,INSEH3,RESVRH, FNDSPH,DSPHIH,RESHOL.
!
  implicit none

  real ( kind = 8 ) d_pi
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) msglvl
  real ( kind = 8 ) tol

  integer ( kind = 4 ) maxfp,maxfv,maxhf,maxht,maxiw,maxpf,maxvc,maxwk
  parameter (maxfp = 200)
  parameter (maxfv = 1000)
  parameter (maxhf = 100)
  parameter (maxht = 307)
  parameter (maxiw = 600)
  parameter (maxpf = 200)
  parameter (maxvc = 300)
  parameter (maxwk = 600)

  integer ( kind = 4 ) facep(3,maxfp),factyp(maxfp),fvl(6,maxfv),hfl(maxhf)
  integer ( kind = 4 ) ht(0:maxht-1),iwk(maxiw),pfl(2,maxpf)
  integer ( kind = 4 ) ccw,f,htsiz,i,ifach,ipolh,imeas,irdr,j,k,l,lp,ls,nbf,nface
  integer ( kind = 4 ) nfach,nfhol,nfph,nihol,npf,npolh,nrfe,nvc,nvch,nvert
  integer ( kind = 4 ) nverth,outmod,prime
  real    holtim,t0,t1
  real ( kind = 8 ) eang(maxfv),nrml(3,maxfp),vcl(3,maxvc),wk(maxwk)
  real ( kind = 8 ) angacc,angmin,aspc2d,atol2d,cp,ntol,rdacc
  real ( kind = 8 ) u(3),v(3)
  logical holint
  character rgname*60
!
!     OUTMOD = 0 : Output updated decomposition in format similar to
!        input but with no holes (suitable for DRds3).
!     OUTMOD = 1 : Output polyhedral decomposition data structure.
!     OUTMOD = 2 : Output VCL + list of vertex indices for each face.
!     MSGLVL >= 2: Routine DSPHFH prints out non-coplanar vertices of
!        each face (also output from routines CUTFAC, INSFAC).
!
  irdr = 5
  imeas = 7
  read (irdr,600) rgname
  read (irdr,*) tol, aspc2d,atol2d,angacc,rdacc

  read (irdr,*) nvc,nface,nfhol,npolh,nihol,outmod,msglvl
  nfhol = max(nfhol,0)
  write (imeas,800) rgname,tol, aspc2d,atol2d,angacc,rdacc,nfhol, &
     nihol
  aspc2d = aspc2d* d_pi ( ) /180.0D+00
  atol2d = atol2d* d_pi ( ) /180.0D+00
  angacc = angacc* d_pi ( ) /180.0D+00
  nfph = nface + nfhol

  if (nvc > maxvc .or. nfph >= maxfp .or. npolh >= maxhf) then
   if (nvc > maxvc) write ( *, 610) 'maxvc',nvc
   if (nfph >= maxfp) write ( *, 610) 'maxfp',nfph+1
   if (npolh >= maxhf) write ( *, 610) 'maxhf',npolh+1
   stop
  end if

  read (irdr,*) ((vcl(j,i),j=1,3),i=1,nvc)
!
!     The outer polygons of NFACE faces should appear first, followed by
!     the inner polygons of NFHOL holes (in order of index of faces
!     containing holes). Face types for NFACE faces should be positive.
!     Face types for NFHOL holes are +F or -F where F is index of face
!     containing hole; positive (negative) sign indicates hole polygon
!     is oriented CCW (CW) in polyhedron when viewed from outside.
!     Holes must only be on boundary faces of polyhedral region.
!
  read (irdr,*) (facep(1,i),i=1,nfph+1)
  nvert = facep(1,nfph+1) - 1
  read (irdr,*) (factyp(i),i=1,nfph)
  if (nvert > maxfv) then
   write ( *, 610) 'maxfv',nvert
   stop
  end if
!
!     Head vertex of each (outer) face must be a strictly convex vertex.
!     Head vertex of each hole may be arbitrary.
!     Positive (negative) sign in PFL(1,*) indicates face is oriented
!     CCW (CW) in polyhedron when viewed from outside.
!     Hole polygons should not be included in PFL.
!
  read (irdr,*) (fvl(1,i),i=1,nvert)
  read (irdr,*) (hfl(i),i=1,npolh+1)
  npf = hfl(npolh+1) - 1
  if (npf + nfhol > maxpf) then
   write ( *, 610) 'maxpf',npf+nfhol
   stop
  end if
  read (irdr,*) (pfl(1,i),i=1,npf)
  htsiz = min(prime(nvc+2),maxht)
  call gtime(t0)
  call dsphfh(aspc2d,atol2d,nvc,nface,nfhol,npolh,maxvc,maxfv,maxiw, &
     maxwk,nvert,npf,vcl,facep,factyp,nrml,fvl,eang,hfl,pfl,htsiz, &
     ht,iwk,wk,ierr)

  call gtime(holtim)

  holtim = holtim - t0

  if (ierr /= 0) then
    write ( *, 620) ierr
    write (imeas,620) ierr
    stop
  end if

  nrfe = 0
  angmin = 2.0D+00* d_pi ( )

  do i = 1,nvert
    if (eang(i) > d_pi ( ) + tol) nrfe = nrfe + 1
    if (eang(i) > -1.0D+00) angmin = min(angmin,eang(i))
  end do

  angmin = angmin*180.0D+00/ d_pi ( )
  write (imeas,810) 'f',nvc,nface,nvert,npolh,npf,nrfe,angmin,holtim
!
!     Read in description of interior hole polyhedra. The input format
!     is similar to above with no hole faces and only 1 polyhedron at
!     a time, treated as though the region for each polyh is its hole.
!     IPOLH is index of polyhedron containing hole. IFACH is index of
!     `extreme' face of hole used for connection to outer boundary.
!     One interior hole polyhedron per polyhedral region is assumed.
!     If > 1, this program can be rerun after finding polyhedron index.
!     A negative sign for IPOLH indicates hole polyh is hole interface.
!
  holtim = 0.0
  do 20 k = 1,nihol
     read (irdr,*) nvch,nfach,ipolh,ifach
   holint = (ipolh < 0)
   ipolh = abs(ipolh)
     if (nvc+nvch > maxvc .or. nface+nfach >= maxfp .or. &
     npf+nfach > maxpf) then
      if (nvc+nvch > maxvc) write ( *, 610) 'maxvc',nvc+nvch
      if (nface+nfach >= maxfp) write ( *, 610) 'maxfp', &
           nface+nfach+1
      if (npf+nfach > maxpf) write ( *, 610) 'maxpf',npf+nfach
      stop
     end if
     read (irdr,*) ((vcl(j,i),j=1,3),i=nvc+1,nvc+nvch)
     read (irdr,*) (facep(1,i),i=nface+1,nface+nfach+1)
     nverth = facep(1,nface+nfach+1) - 1
     read (irdr,*) (factyp(i),i=nface+1,nface+nfach)
     if (nvert+nverth > maxfv) then
      write ( *, 610) 'maxfv',nvert+nverth
      stop
     end if
   read (irdr,*) (fvl(1,i),i=nvert+1,nvert+nverth)
     read (irdr,*) (pfl(1,i),i=npf+1,npf+nfach)
     htsiz = min(prime(nvch+2),maxht)
     call gtime(t0)
     call dsphih(aspc2d,atol2d,angacc,rdacc,nvc,nface,nvert,npolh, &
        npf,nvch,nfach,ipolh,ifach,holint,maxvc,maxfp,maxfv,maxhf, &
        maxpf,maxiw,maxwk,vcl,facep,factyp,nrml,fvl,eang,hfl,pfl, &
        htsiz,ht,iwk,wk,ierr)

     call gtime(t1)
     holtim = holtim + (t1 - t0)
     if (ierr /= 0) then
      write ( *, 620) ierr
      write (imeas,620) ierr
      stop
     end if
   20 continue
  if (nihol > 0) then
     nrfe = 0
     angmin = 2.0D+00* d_pi ( )
     do i = 1,nvert
      if (eang(i) > d_pi ( ) + tol) nrfe = nrfe + 1
        if (eang(i) > -1.0D+00) angmin = min(angmin,eang(i))
     end do
     angmin = angmin*180.0D+00/ d_pi ( )
     write (imeas,810) 'i',nvc,nface,nvert,npolh,npf,nrfe,angmin, &
        holtim
  end if

  if (outmod == 0) then
     if (max(nface,npolh) > maxiw) then
      write ( *, 610) 'maxiw',max(nface,npolh)
      stop
     end if
!
!  Set FACEP(1,I) to point to a strictly convex vertex.
!
   do 50 i = 1,nface
      k = 1
      if (abs(nrml(2,i)) > abs(nrml(1,i))) k = 2
      if (abs(nrml(3,i)) > abs(nrml(k,i))) k = 3
      if (facep(2,i) > 0) then
         ccw = 3
      else
         ccw = 4
      end if
      j = facep(1,i)
      l = fvl(1,j)
      lp = fvl(1,fvl(7-ccw,j))
      u(1) = vcl(1,l) - vcl(1,lp)
      u(2) = vcl(2,l) - vcl(2,lp)
      u(3) = vcl(3,l) - vcl(3,lp)
   40       continue
         ls = fvl(1,fvl(ccw,j))
         v(1) = vcl(1,ls) - vcl(1,l)
         v(2) = vcl(2,ls) - vcl(2,l)
         v(3) = vcl(3,ls) - vcl(3,l)
         ntol = tol*max(abs(u(1)), abs(u(2)), abs(u(3)), &
              abs(v(1)), abs(v(2)), abs(v(3)))
         if (k == 1) then
      cp = u(2)*v(3) - u(3)*v(2)
         else if (k == 2) then
      cp = u(3)*v(1) - u(1)*v(3)
         else
      cp = u(1)*v(2) - u(2)*v(1)
         end if
      if (abs(cp) <= ntol .or. cp*nrml(k,i) < 0.0D+00) then
         j = fvl(ccw,j)
         l = ls
         u(1) = v(1)
         u(2) = v(2)
         u(3) = v(3)
         go to 40
      end if
      facep(1,i) = j
   50    continue
     write ( *, 600) rgname
     write ( *, 690) nvc,nface,npolh
   write ( *, '(a)' )
   write ( *, 700) ((vcl(j,i),j=1,3),i=1,nvc)
   write ( *, '(a)' )
   do i = 1,nface
      k = 0
      j = facep(1,i)
   60       continue
         k = k + 1
         j = fvl(3,j)
      if (j /= facep(1,i)) go to 60
      iwk(i) = k
   end do
   iwk(1) = iwk(1) + 1
   do i = 2,nface
      iwk(i) = iwk(i) + iwk(i-1)
   end do
   write ( *, 710) 1,(iwk(i),i=1,nface)
   write ( *, 710) (factyp(i),i=1,nface)
   write ( *, '(a)' )
   do 100 i = 1,nface
      k = 0
      j = facep(1,i)
   90       continue
         k = k + 1
         iwk(k) = fvl(1,j)
         j = fvl(3,j)
      if (j /= facep(1,i)) go to 90
      if (k > maxiw) then
         write ( *, 610) 'maxiw',k
         stop
      end if
      write ( *, 710) (iwk(j),j=1,k)
  100    continue
   write ( *, '(a)' )
   do i = 1,npolh
      k = 0
      j = hfl(i)
  110       continue
         k = k + 1
         j = pfl(2,j)
      if (j /= hfl(i)) go to 110
      iwk(i) = k
   end do
   iwk(1) = iwk(1) + 1
   do i = 2,npolh
      iwk(i) = iwk(i) + iwk(i-1)
   end do
   write ( *, 710) 1,(iwk(i),i=1,npolh)
   write ( *, '(a)' )
   do i = 1,npolh
      k = 0
      j = hfl(i)
  140       continue
         k = k + 1
         iwk(k) = pfl(1,j)
         j = pfl(2,j)
      if (j /= hfl(i)) go to 140
      write ( *, 710) (iwk(j),j=1,k)
    end do
  else if (outmod == 1) then
     write ( *, 630) nvc,nface,nvert,npolh,npf
     write ( *, 640) (i,(vcl(j,i),j=1,3),i=1,nvc)
     write ( *, 650) (i,(facep(j,i),j=1,3),factyp(i), &
        (nrml(j,i),j=1,3),i=1,nface)
     write ( *, 660) (i,(fvl(j,i),j=1,6),eang(i),i=1,nvert)
     write ( *, 670) (hfl(i),i=1,npolh)
     write ( *, 680) (i,(pfl(j,i),j=1,2),i=1,npf)
  else
   nbf = 0
   do f = 1,nface
      if (facep(3,f) == 0) nbf = nbf + 1
   end do
   write ( *, 690) nvc,nbf,nface-nbf
   write ( *, 720) ((vcl(j,i),j=1,3),i=1,nvc)

   do i = 1,2
      write ( *, '(a)' )
      do 180 f = 1,nface
         if (i == 1 .and. facep(3,f) /= 0) go to 180
         if (i == 2 .and. facep(3,f) == 0) go to 180
         k = 0
         j = facep(1,f)
  170          continue
      k = k + 1
      iwk(k) = fvl(1,j)
      j = fvl(3,j)
         if (j /= facep(1,f)) go to 170
         if (k > maxiw) then
      write ( *, 610) 'maxiw',k
      stop
         end if
         write ( *, 730) k,(iwk(j),j=1,k)
  180       continue
    end do
  end if

  600 format (a60)
  610 format (1x,'*** ',a,' must be increased to',i8)
  620 format (/1x,'ierr=',i5)
  630 format (1x,'nvc,nface,nvert,npolh,npf'/1x,5i7)
  640 format (/' vcl'/(1x,i7,3f15.7))
  650 format (/' facep,factyp,nrml'/(1x,4i7,i5,3f15.7))
  660 format (/' fvl,eang'/(1x,7i7,f15.7))
  670 format (/' hfl'/(1x,10i7))
  680 format (/' pfl'/(1x,3i7))
  690 format (3i7)
  700 format ((3f25.15))
  710 format ((10i7))
  720 format (/(1x,3f15.7))
  730 format (1x,11i7/(8x,10i7))
  800 format (1x,a60/1x,'input : tol=',d15.7,'   aspc2d=',f9.3, &
     '   atol2d=',f9.3/9x,'angacc=',f9.3,'   rdacc=',f9.5, &
     '   nfhol=',i3,'   nihol=',i3)
  810 format (1x,a1,'holds: nvc=',i7,'   nface=',i7,'   nvert=',i7, &
     '   npolh=',i7/9x,'npf=',i7,'   nrfe=',i7,'   angmin=',f9.3, &
     '   holtim=',f9.3)

  return
end
subroutine test_dtr3 ( )

!*****************************************************************************80
!
!! TEST_DTR3 tests ...
!
!     Driver program for testing routines DTRIS3,DTRIW3,IMPTR3,ITRIS3,
!        BNSRT3,VBFAC,WALKT3,NWTHOU,NWTHIN,NWTHFC,NWTHED,FRSTET,SWAPES,
!        SWAPMU,TETMU,TETLST,IMPTRF,SWAPTF,FNDMSW,IMPTRD,FNMSWD
!
  implicit none

  integer ( kind = 4 ), parameter :: maxalg = 10
  integer ( kind = 4 ), parameter :: maxbf = 3000
  integer ( kind = 4 ), parameter :: maxfc = 80000
  integer ( kind = 4 ), parameter :: maxht = 8011
  integer ( kind = 4 ), parameter :: maxth = 40000
  integer ( kind = 4 ), parameter :: maxvc = 5000

  logical bndcon
  real ( kind = 8 ) d_pi
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) msglvl
  real ( kind = 8 ) tol

  integer ( kind = 4 ) bf(3,maxbf),crit(maxalg),fc(7,maxfc),ht(0:maxht-1)
  integer ( kind = 4 ) tetra(4,maxth),vm(maxvc)
  integer ( kind = 4 ) a,alg,axis,b,c,d,e,i,ialg,imeas,in,inmode,irdr,j,mlvl,nalg
  integer ( kind = 4 ) nav,nbf,nface,nfc,nlo,npt,nt,ntetra,nx,nz,prime,seed,sizht
  real    t0,tritim
  real ( kind = 8 ) sang(4),scal(3),tran(3),vcl(3,maxvc)
  real ( kind = 8 ) binexp,eta,mineta,minrho,minsig,radsq,rho,sigma
  real ( kind = 8 ) sx,sz, emnrth,radrth,sangmn
!
!  INMODE = 1: read vertices;
!  INMODE = 2: generate random vertices;
!  INMODE = 3: worst case set of vertices.
!
!  ALG = 1: DTRIS3;
!  ALG = 2: DTRIW3;
!  ALG = 3: DTRIW3 with call to BNSRT3 first;
!  ALG = 4: ITRIS3.
!
!  NALG = 1: construct Delaunay triang;
!         2 <= ABS(NALG) <= MAXALG:
!           improve based on local max-min empty circumsphere, solid angle,
!           radius ratio, mean ratio criterion as specified by CRIT(I) = 0,
!           1, 2, 3, resp, 2 <= I <= ABS(NALG), i.e. there are ABS(NALG)-1
!           improvements;
!  NALG < -1: improvements are boundary-constrained.
!
!  MSGLVL = -1: print only measurements and allow multiple runs;
!  MSGLVL =  0: print arrays;
!  MSGLVL =  2: print list of tetrahedra in output units 1, 2, 3;
!  MSGLVL =  4: print new tetrahedra,
!               local optimality tests, local transformations.
!
  irdr = 5
  imeas = 7

10 continue

  read (irdr,*) inmode,alg,nalg,mlvl,npt,tol,binexp

  if ( inmode <= 0 .or. alg <= 0 .or. alg > 4 ) then
    stop
  end if

  if ( maxvc < npt ) then
    write ( *, 600) 'maxvc',npt
    stop
  end if

  msglvl = mlvl
  write (imeas,700) inmode,alg,nalg,npt,tol, binexp
  bndcon = (nalg < 0)
  nalg = min(abs(nalg),maxalg)
  if (nalg > 1) read (irdr,*) (crit(i),i=2,nalg)

  if (alg == 4) then
    crit(1) = -1
  else
    crit(1) = 0
  end if

  if (inmode == 1) then
    read (irdr,*) (vcl(1,i),vcl(2,i),vcl(3,i),i=1,npt)
    sizht = prime(npt*3/2)
  else if (inmode == 2) then
    read (irdr,*) seed,axis,nav,(scal(i),i=1,3),(tran(i),i=1,3)
    write (imeas,710) seed,axis,nav,(scal(i),i=1,3),(tran(i),i=1,3)
    call randpt(3,npt,seed,axis,nav,scal,tran,3,vcl)
    sizht = prime(npt*3/2)
  else
    nx = npt/2
    nz = npt - nx
    sx = 2.0D+00* d_pi ( ) / nx
    sz = 1.0D+00/(nz - 1)

    do i = 1, nx
      vcl(1,i) = cos ( i * sx )
      vcl(2,i) = sin ( i * sx )
      vcl(3,i) = 0.0D+00
    end do

    do i = 1, nz
      vcl(1,nx+i) = 0.0D+00
      vcl(2,nx+i) = 0.0D+00
      vcl(3,nx+i) = (i - 1)*sz
    end do

    sizht = prime(npt**2/20)

  end if

  sizht = min(sizht,maxht)

  do i = 1,npt
    vm(i) = i
  end do

  if (msglvl >= 0) write ( *, 650) npt,(i,vcl(1,i),vcl(2,i), &
     vcl(3,i),i=1,npt)

  do 70 ialg = 1,nalg

    call gtime(t0)

    if (ialg == 1) then

      if (alg == 1) then
        call dtris3(npt,sizht,maxbf,maxfc,vcl,vm,nbf,nfc,nface, &
          ntetra,bf,fc,ht,ierr)
      else if (alg <= 3) then
        if (alg == 3) call bnsrt3(binexp,npt,vcl,vm,fc,ht)
        call dtriw3(npt,sizht,maxbf,maxfc,vcl,vm,nbf,nfc,nface, &
          ntetra,bf,fc,ht,ierr)
      else
        call itris3(npt,sizht,maxbf,maxfc,vcl,vm,nbf,nfc,nface, &
          ntetra,bf,fc,ht,ierr)
      end if
    else
      call imptr3(bndcon,.true.,crit(ialg),npt,sizht,maxfc,vcl,vm, &
        nfc,ntetra,bf,fc,ht,nface,ierr)
    end if

    call gtime(tritim)
    tritim = tritim - t0

    if (ierr /= 0) then
      write ( *, 610) ierr
      write (imeas,610) ierr
      stop
    end if

    if ( maxth < ntetra ) then
      write ( *, 600) 'maxth',ntetra
      stop
    end if

    call tetlst(nfc,vm,fc,nt,tetra)

    if (nt /= ntetra) then
      write ( *, 620) nt,ntetra
      write (imeas,620) nt,ntetra
      stop
    end if

   nlo = 0
   do i = 1,nfc
      if (fc(1,i) > 0 .and. fc(5,i) > 0) then
         a = vm(fc(1,i))
         b = vm(fc(2,i))
         c = vm(fc(3,i))
         d = vm(fc(4,i))
         e = vm(fc(5,i))
         call ccsph(.true.,vcl(1,a),vcl(1,b),vcl(1,c),vcl(1,d), &
              vcl(1,e),sang,radsq,in)
         if (in >= 1) then
      nlo = nlo + 1
      if (msglvl == 4) then
         if (nlo == 1) write ( *, 630)
         write ( *, 640) i,a,b,c,d,e
      end if
         end if
      end if
  end do

   minsig = 2.0D+00
   minrho = 2.0D+00
   mineta = 2.0D+00
   do i = 1,nt
      a = tetra(1,i)
      b = tetra(2,i)
      c = tetra(3,i)
      d = tetra(4,i)
      sigma = sangmn(vcl(1,a),vcl(1,b),vcl(1,c),vcl(1,d),sang)
      rho = radrth(vcl(1,a),vcl(1,b),vcl(1,c),vcl(1,d))
      eta = emnrth(vcl(1,a),vcl(1,b),vcl(1,c),vcl(1,d))
      minsig = min(minsig,sigma)
      minrho = min(minrho,rho)
      mineta = min(mineta,eta)
   end do

   minsig = minsig*1.5D+00*sqrt(6.0D+00)

     write (imeas,720) crit(ialg),nbf,nfc,nface,ntetra,sizht,nlo, &
        minsig,minrho,mineta,tritim
   if (msglvl == 2 .and. ialg <= 4) &
        write (ialg,660) ((tetra(j,i),j=1,4),i=1,nt)
     if (msglvl >= 0) then
      if (ialg == 1) write ( *, 670) 'vm',npt,(vm(i),i=1,npt)
      write ( *, 670) 'ht',sizht,(ht(i),i=0,sizht-1)
      write ( *, 680) nfc,(i,(fc(j,i),j=1,7),i=1,nfc)
      if (ialg == 1 .or. .not. bndcon) &
           write ( *, 690) nbf,(i,(bf(j,i),j=1,3),i=1,nbf)
   end if
   70 continue
  if (msglvl < 0) go to 10

  600 format (1x,'*** ',a,' must be increased to',i8)
  610 format (/1x,'ierr=',i5)
  620 format (/1x,'error from tetlst, nt != ntetra :',2i7)
  630 format (1x,'nonlocally optimal faces - sphere criterion')
  640 format (1x,'face',i7,': a,b,c,d,e=',5i7)
  650 format (1x,'vcl   ',i7/(1x,i7,3f23.15))
  660 format (1x,4i7)
  670 format (/1x,a,3x,i7/(1x,10i7))
  680 format (/1x,'fc   ',i7/(1x,8i7))
  690 format (/1x,'bf   ',i7/(1x,4i7))
  700 format (/1x,'inm=',i1,'   alg=',i1,'   nalg=',i2,'   npt=',i7, &
     '   tol=',d14.7,'   binexp=',f9.7)
  710 format (1x,'seed=',i9,'   axis=',i2,'   nav=',i7/1x, &
     'scal(1:3),tran(1:3)=',6f9.4)
  720 format (1x,'crit=',i2,'   nbf=',i5,'   nfc=',i6,'   nface=',i6, &
     '   ntetra=',i6,'   sizht=',i5/1x,'nlo=',i5,'   msig=',f9.7, &
     '   mrho=',f9.7,'   meta=',f9.7,'   tim=',f9.3)

  return
end
subroutine test_dtri ( )

!*****************************************************************************80
!
!! TEST_DTRI tests DTRIS2, DTRIW2, SWAPEC, VBEDG, WALKT2, BNSRT2.
!
  implicit none

  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) msglvl
  real ( kind = 8 ) tol

  integer ( kind = 4 ) maxnp,maxst
  real ( kind = 8 ) large
  parameter (maxnp = 10000)
  parameter (maxst = 100)
  parameter (large = 1000.0D+00)

  integer ( kind = 4 ) ind(maxnp+3),stack(maxst),til(3,maxnp*2+1)
  integer ( kind = 4 ) tnbr(3,maxnp*2+1)
  integer ( kind = 4 ) a,alg,axis,b,c,d,diaedg,i,imeas,inmode,irdr,j,jp1,jp2,k
  integer ( kind = 4 ) mlvl,nav,nlo,npt,ntri,seed
  real ( kind = 8 ) t0,tritim
  real ( kind = 8 ) scal(2),tran(2),vcl(2,maxnp+3)
  real ( kind = 8 ) binexp
!
!     INMODE = 1: read vertices; INMODE = 2: generate random vertices.
!     ALG = 1: DTRIS2; ALG = 2: DTRIW2; ALG = 3: DTRIW2 with bounding
!        triangle; ALG = 4: DTRIW2 with call to BNSRT2 first.
!     MSGLVL = -1: print only measurements and allow multiple runs;
!     MSGLVL = 0: print arrays; MSGLVL = 4: also print edges as they
!        are created and swapped.
!
  irdr = 5
  imeas = 7

10 continue

  read (irdr,*) inmode,alg,mlvl,npt,tol,binexp

  if (inmode <= 0 .or. alg <= 0 .or. alg > 4) then
    stop
  end if

  if (npt > maxnp) then
    write ( *, 600) 'maxnp',npt
    stop
  end if

  msglvl = mlvl
  write (imeas,700) alg,npt,tol, binexp

  if (inmode == 1) then
    read (irdr,*) (vcl(1,i),vcl(2,i),i=1,npt)
  else
    read (irdr,*) seed,axis,nav,scal(1),scal(2),tran(1),tran(2)
    write (imeas,710) seed,axis,nav,scal(1),scal(2),tran(1),tran(2)
    call randpt(2,npt,seed,axis,nav,scal,tran,2,vcl)
  end if

  if (alg /= 3) then
    do i = 1,npt
      ind(i) = i
    end do
  else
    vcl(1,npt+1) = -large
    vcl(2,npt+1) = -large
    vcl(1,npt+2) = large
    vcl(2,npt+2) = -large
    vcl(1,npt+3) = 0.0D+00
    vcl(2,npt+3) = large
    ind(1) = npt + 1
    ind(2) = npt + 2
    ind(3) = npt + 3

    do i = 1,npt
      ind(i+3) = i
    end do

    npt = npt + 3
  end if

  if (msglvl >= 0) then
    write ( *, 620) msglvl,npt
    write ( *, 630) (i,vcl(1,i),vcl(2,i),i=1,npt)
    if (msglvl == 4) write ( *, 620)
  end if

  call gtime(t0)

  if (alg == 1) then
    call dtris2(npt,maxst,vcl,ind,ntri,til,tnbr,stack,ierr)
  else
    if (alg == 4) call bnsrt2(binexp,npt,vcl,ind,til,tnbr)
    call dtriw2(npt,maxst,vcl,ind,ntri,til,tnbr,stack,ierr)
  end if

  call gtime(tritim)
  tritim = tritim - t0

  if (ierr /= 0) then
    write ( *, 610) ierr
    write (imeas,610) ierr
    stop
  end if

  nlo = 0
  do i = 1,ntri
    do j = 1,3
      k = tnbr(j,i)
      if (k > i) then
         jp1 = j + 1
         if (jp1 > 3) jp1 = 1
         jp2 = jp1 + 1
         if (jp2 > 3) jp2 = 1
         a = til(j,i)
         b = til(jp1,i)
         c = til(jp2,i)
         if (til(1,k) == b) then
           d = til(3,k)
         else if (til(2,k) == b) then
           d = til(1,k)
         else
           d = til(2,k)
         end if
         if (diaedg(vcl(1,c),vcl(2,c),vcl(1,a),vcl(2,a),vcl(1,d), &
              vcl(2,d),vcl(1,b),vcl(2,b)) == 1) nlo = nlo + 1
      end if
    end do
  end do

  write (imeas,720) ntri,nlo,tritim

  if (msglvl >= 0) write ( *, 640) ntri,(i,(til(j,i),j=1,3), &
     (tnbr(j,i),j=1,3),i=1,ntri)
  if (msglvl < 0) go to 10

  600 format (1x,'*** ',a,' must be increased to',i8)
  610 format (/1x,'ierr=',i5)
  620 format (1x,2i7)
  630 format (1x,i7,2f15.7)
  640 format (/1x,i7/(1x,4i7,3i8))
  700 format (/1x,'alg=',i2,'   npt=',i7,'   tol=',d15.7,'   binexp=', &
     f13.7)
  710 format (1x,'seed=',i9,'   axis=',i2,'   nav=',i7/1x, &
     'scal(1:2),tran(1:2)=',4f12.6)
  720 format (1x,'ntri=',i7,'   nlo=',i5,'   tritim=',f9.3)

  return
end
subroutine test_dtris3 ( )

!*****************************************************************************80
!
!! TEST_DTRIS3 tests DTRIS3.
!
  implicit none

  integer ( kind = 4 ), parameter :: bf_max = 100
  integer ( kind = 4 ), parameter :: fc_max = 100
  integer ( kind = 4 ), parameter :: tetra_max = 100
  integer ( kind = 4 ), parameter :: node_max = 100
  integer ( kind = 4 ), parameter :: node_num = 6

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) bf(3,bf_max)
  integer ( kind = 4 ) bf_num
  integer ( kind = 4 ) c
  real ( kind = 8 ) center(3)
  integer ( kind = 4 ) d
  real ( kind = 8 ) d_pi
  integer ( kind = 4 ) e
  real ( kind = 8 ) emnrth
  real ( kind = 8 ) eta
  real ( kind = 8 ) eta_min
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) fc(7,fc_max)
  integer ( kind = 4 ) fc_num
  integer ( kind = 4 ), allocatable, dimension ( : ) :: ht
  integer ( kind = 4 ) ht_num
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) indicator
  integer ( kind = 4 ) j
  integer ( kind = 4 ) :: msglvl = 4
  integer ( kind = 4 ) prime
  real ( kind = 8 ) radrth
  real ( kind = 8 ) radius_sq
  real ( kind = 8 ) rho
  real ( kind = 8 ) rho_min
  real ( kind = 8 ) sang(4)
  real ( kind = 8 ) sangmn
  integer ( kind = 4 ) seed
  real ( kind = 8 ) sigma
  real ( kind = 8 ) sigma_min
  integer ( kind = 4 ) tetra(4,tetra_max)
  integer ( kind = 4 ) tetra_bad_num
  integer ( kind = 4 ) tetra_num
  integer ( kind = 4 ) tetra_num2
  real ( kind = 8 ), dimension(3,node_num) :: vcl = reshape ( (/ &
     0.0D+00,  0.0D+00, -1.0D+00, &
     0.0D+00, -1.0D+00,  0.0D+00, &
     1.0D+00,  0.0D+00,  0.0D+00, &
     0.0D+00,  1.0D+00,  0.0D+00, &
    -1.0D+00,  0.0D+00,  0.0D+00, &
     0.0D+00,  0.0D+00,  1.0D+00 /), (/ 3, node_num /) )
  integer ( kind = 4 ), allocatable, dimension ( : ) :: vm

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_DTRIS3 :'
  write ( *, '(a)' ) '  DTRIS3 computes the Delaunay "triangulation"'
  write ( *, '(a)' ) '  of a set of points in 3D.'

  call r8mat_transpose_print ( 3, node_num, vcl, '  Node coordinates:' )

  ht_num = prime ( node_num * 3 / 2 )

  allocate ( ht(0:ht_num-1) )
  allocate ( vm(1:node_num) )
!
!  Assign the initial ordering to the nodes.
!
  do i = 1, node_num
    vm(i) = i
  end do

  call dtris3 ( node_num, ht_num, bf_max, fc_max, vcl, vm, bf_num, fc_num, &
    face_num, tetra_num, bf, fc, ht, ierr )

  if ( ierr /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) '  IERR = ', ierr
    return
  end if

  if ( tetra_max < tetra_num ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) '  TETRA_MAX must be increased to ', tetra_num
    return
  end if
!
!  Construct a list of the tetrahedra from the FC array.
!
  call tetlst ( fc_num, vm, fc, tetra_num2, tetra )

  if ( tetra_num2 /= tetra_num ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Error from TETLST!'
    write ( *, '(a,i6)' ) &
      '  Expected number of tetrahedra TETRA_NUM = ', tetra_num
    write ( *, '(a,i6)' ) '  Number of tetrahedron found = ', tetra_num2
    return
  end if

  call imat_transpose_print ( 4, tetra_num, tetra, '  Tetrahedron indices:' )
!
!  Find the circumsphere through each tetrahedron, and check
!  whether point E is inside tetrahedron A, B, C, D.
!
  tetra_bad_num = 0

  do i = 1, fc_num

    if ( 0 < fc(1,i) .and. 0 < fc(5,i) ) then

      a = vm(fc(1,i))
      b = vm(fc(2,i))
      c = vm(fc(3,i))
      d = vm(fc(4,i))
      e = vm(fc(5,i))

      call ccsph ( .true., vcl(1,a), vcl(1,b), vcl(1,c), vcl(1,d), &
        vcl(1,e), center, radius_sq, indicator )

      if ( 1 <= indicator ) then

        tetra_bad_num = tetra_bad_num + 1

        if ( msglvl == 4 ) then

          if ( tetra_bad_num == 1 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) '  Nonlocally optimal faces - sphere criterion'
          end if

          write ( *, '(2x,i6,2x,i6,2x,i6,2x,i6,2x,i6,2x,i6)' ) &
            i, a, b, c, d, e

        end if

      end if
    end if

  end do
!
!  Compute the minimum quality measure.
!
  eta_min = 1.0D+00
  rho_min = 1.0D+00
  sigma_min = 1.0D+00

  do i = 1, tetra_num

    a = tetra(1,i)
    b = tetra(2,i)
    c = tetra(3,i)
    d = tetra(4,i)

    eta = emnrth ( vcl(1,a), vcl(1,b), vcl(1,c), vcl(1,d) )

    rho = radrth( vcl(1,a), vcl(1,b), vcl(1,c), vcl(1,d) )

    sigma = sangmn ( vcl(1,a), vcl(1,b), vcl(1,c), vcl(1,d), sang )

    eta_min = min ( eta_min, eta )
    rho_min = min ( rho_min, rho )
    sigma_min = min ( sigma_min, sigma )

  end do

  sigma_min = sigma_min * 1.5D+00 * sqrt ( 6.0D+00 )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  BF_NUM = ', bf_num
  write ( *, '(a,i12)' ) '  FC_NUM = ', fc_num
  write ( *, '(a,i12)' ) '  FACE_NUM = ', face_num
  write ( *, '(a,i12)' ) '  TETRA_NUM = ', tetra_num
  write ( *, '(a,i12)' ) '  TETRA_BAD_NUM = ', tetra_bad_num
  write ( *, '(a,i12)' ) '  HT_NUM = ', ht_num
  write ( *, '(a,g14.6)' ) '  min ( ETA ) = ', eta_min
  write ( *, '(a,g14.6)' ) '  min ( RHO ) = ', rho_min
  write ( *, '(a,g14.6)' ) '  min ( SIGMA ) = ', sigma_min

  if ( msglvl == 2 ) then
    do j = 1, tetra_num
      write (*,'(2x,4i7)' ) tetra(1:4,j)
    end do
  end if

  if (msglvl >= 0) then
      write ( *, 670) 'vm',node_num,(vm(i),i=1,node_num)
    write ( *, 670) 'ht',ht_num,(ht(i),i=0,ht_num-1)
    write ( *, 680) fc_num,(i,(fc(j,i),j=1,7),i=1,fc_num)
    write ( *, 690) bf_num,(i,(bf(j,i),j=1,3),i=1,bf_num)
  end if

  670 format (/1x,a,3x,i7/(1x,10i7))
  680 format (/1x,'FC   ',i7/(1x,8i7))
  690 format (/1x,'BF   ',i7/(1x,4i7))

  deallocate ( ht )
  deallocate ( vm )

  return
end
subroutine test_dtrk ( )

!*****************************************************************************80
!
!! TEST_DTRK
!
!     Driver program for testing routines DTRISK, DTRIWK, BNSRTK,
!        VBFACK, WALKTK, NWSXOU, NWSXIN, NWSXFC, NWSXED, FRSMPX, SWAPDG,
!        SWAPHS, SMPXLS.
!
  implicit none

  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) msglvl
  real ( kind = 8 ) tol

  integer ( kind = 4 ) kmax,maxbf,maxfc,maxht,maxsm,maxvc
  parameter (kmax = 6)
  parameter (maxbf = 25000)
  parameter (maxfc = 110000)
  parameter (maxht = 14011)
  parameter (maxsm = 40000)
  parameter (maxvc = 5000)

  integer ( kind = 4 ) bf(kmax*maxbf),fc((kmax+4)*maxfc),ht(0:maxht-1)
  integer ( kind = 4 ) iwk(6*kmax+2),smplx((kmax+1)*maxsm),vm(maxvc)
  integer ( kind = 4 ) alg,axis,e,i,imeas,in,inmode,irdr,j,k,kp1,kp2,kp4,mlvl,nav
  integer ( kind = 4 ) nbf,nbfac,nface,nfc,nlo,npt,ns,nsmplx,prime,seed,sizht
  real ( kind = 8 ) t0,tritim
  real ( kind = 8 ) urand
  real ( kind = 8 ) scal(kmax),tran(kmax),vcl(kmax*maxvc)
  real ( kind = 8 ) wk(kmax*kmax+2*kmax+1)
  real ( kind = 8 ) binexp,r,radsq
!
!     INMODE = 1: read vertices; INMODE = 2: generate random vertices;
!        INMODE = 3: worst case set of vertices (on moment curve).
!     ALG = 1: DTRISK; ALG = 2: DTRIWK with call to BNSRTK first if
!        BINEXP > 0; ALG = 3: DTRIMK with call to BNSRTK first if
!        BINEXP > 0.
!     MSGLVL = -1: print only measurements and allow multiple runs;
!        MSGLVL = 0: print arrays; MSGLVL = 2: print list of simplices
!        in output unit ALG; MSGLVL = 4: print new simplices, local
!        optimality tests, local transformations or deleted simplices.
!
  irdr = 5
  imeas = 7

10 continue

  read (irdr,*) k,inmode,alg,mlvl,npt,tol, binexp

  if (k < 2 .or. inmode <= 0 .or. alg <= 0 .or. alg > 3) then
    stop
  end if

  if (k > kmax) then
   write ( *, 600) 'kmax',k
   stop
  else if (npt > maxvc) then
   write ( *, 600) 'maxvc',npt
   stop
  end if

  kp1 = k + 1
  kp2 = k + 2
  kp4 = k + 4
  msglvl = mlvl
  write (imeas,800) k,inmode,alg,npt,tol, binexp
  if (inmode == 1) then
   read (irdr,*) (vcl(i),i=1,k*npt)
   j = npt*0.375*(4.5**(k-2)) + 0.5
   sizht = prime(j)
  else if (inmode == 2) then
   read (irdr,*) seed,axis,nav,(scal(i),i=1,k),(tran(i),i=1,k)
     write (imeas,810) seed,axis,nav,(scal(i),i=1,k)
   write (imeas,820) (tran(i),i=1,k)
   call randpt(k,npt,seed,axis,nav,scal,tran,k,vcl)
   j = npt*0.375*(4.5**(k-2)) + 0.5
   sizht = prime(j)
  else
   read (irdr,*) seed,scal(1),tran(1)
     write (imeas,830) seed,scal(1),tran(1)
   j = 0

   do i = 1,npt
      r = urand(seed)*scal(1) + tran(1)
      j = j + 1
      vcl(j) = r
      do e = 2,k
         j = j + 1
         vcl(j) = vcl(j-1)*r
      end do
  end do

   if (k == 2) then
      j = (npt + npt - 3)/8
   else if (k == 3) then
      j = (npt - 2)**2/8
   else if (k == 4) then
      j = (npt - 3)*(3*npt - 10)/16
   else if (k == 5) then
      j = (npt - 3)*(npt - 4)**2/16
   else
      j = (npt - 4)*(npt - 5)*(4*npt - 21)/48
   end if
   sizht = prime(j)
  end if
  sizht = min(sizht,maxht)

  do i = 1,npt
    vm(i) = i
  end do

  if (msglvl >= 0) then
    write ( *, 650) npt
    do i = 1,npt
        if (k <= 4) then
         write ( *, 660) i,(vcl((i-1)*k+j),j=1,k)
      else
         write ( *, 670) i,(vcl((i-1)*k+j),j=1,k)
      end if
    end do
  end if

  call gtime(t0)

  if (alg == 1) then

    call dtrisk(k,npt,sizht,maxbf*kmax/k,maxfc*(kmax+4)/kp4,vcl, &
      vm,nbf,nfc,nbfac,nface,nsmplx,bf,fc,ht,iwk,wk,ierr)

  else if (alg == 2) then

    if ( binexp > 0.0D+00 ) then
      call bnsrtk(k,binexp,npt,vcl,vm,fc,ht,wk)
    end if

    call dtriwk(k,npt,sizht,maxbf*kmax/k,maxfc*(kmax+4)/kp4,vcl, &
        vm,nbf,nfc,nbfac,nface,nsmplx,bf,fc,ht,iwk,wk,ierr)

  else

    if ( binexp > 0.0D+00 ) then
      call bnsrtk(k,binexp,npt,vcl,vm,fc,ht,wk)
    end if

    call dtrimk(k,npt,sizht,maxbf*kmax/k,maxfc*(kmax+4)/kp4,vcl, &
      vm,nbf,nfc,nbfac,nface,nsmplx,bf,fc,ht,iwk,wk,ierr)

  end if

  call gtime(tritim)
  tritim = tritim - t0

  if (ierr /= 0) then
     write ( *, 610) ierr
   write (imeas,610) ierr
   stop
  end if

  if (msglvl == 2) then

    if (nsmplx > maxsm*(kmax+1)/(kp1)) then
      write ( *, 600) 'maxsm',nsmplx
      stop
    end if

    call smpxls(k,nfc,vm,fc,ns,smplx)

    if (ns /= nsmplx) then
      write ( *, 620) ns,nsmplx
      write (imeas,620) ns,nsmplx
      stop
    end if

  end if

  nlo = 0
  do 70 i = 1,nfc
   if (fc((i-1)*kp4+1) > 0 .and. fc((i-1)*kp4+kp2) > 0) then
      do j = 1,kp1
         iwk(j) = vm(fc((i-1)*kp4+j))
      end do
      e = vm(fc((i-1)*kp4+kp2))
      call ccsphk(k,.true.,iwk,vcl,vcl((e-1)*k+1),wk,radsq,in, &
           wk(kp1),iwk(kp2))
      if (in >= 1) then
         nlo = nlo + 1
         if (msglvl == 4) then
      if (nlo == 1) write ( *, 630)
      write ( *, 640) i,(iwk(j),j=1,kp1),e
         end if
      end if
   end if
   70 continue

  write (imeas,840) nbf,nbfac,nfc,nface,nsmplx,sizht,nlo,tritim

  if (msglvl == 2) then
    do 80 i = 1,ns
      write (alg,680) (smplx((i-1)*kp1+j),j=1,kp1)
   80    continue
  end if

  if (msglvl >= 0) then
    write ( *, 690) 'vm',npt,(vm(i),i=1,npt)
    write ( *, 690) 'ht',sizht,(ht(i),i=0,sizht-1)
    write ( *, 700) nfc
    do i = 1,nfc
      write ( *, 680) i,(fc((i-1)*kp4+j),j=1,kp4)
    end do
    write ( *, 710) nbf
    do i = 1,nbf
      write ( *, 680) i,(bf((i-1)*k+j),j=1,k)
    end do
  end if
  if (msglvl < 0) go to 10

  600 format (1x,'*** ',a,' must be increased to',i8)
  610 format (/1x,'ierr=',i5)
  620 format (/1x,'error from smpxls, ns != nsmplx :',2i7)
  630 format (1x,'nonlocally optimal faces')
  640 format (1x,'face',i6,': ab...c, d,e=',8i6)
  650 format (1x,'vcl   ',i7)
  660 format (1x,i7,4f15.7)
  670 format (1x,i7,4f15.7/(8x,4f15.7))
  680 format (1x,11i7)
  690 format (/1x,a,3x,i7/(1x,10i7))
  700 format (/1x,'fc   ',i7)
  710 format (/1x,'bf   ',i7)
  800 format (/1x,'k=',i2,'   inm=',i1,'   alg=',i1,'   npt=',i7, &
     '   tol=',d14.7,'   binexp=',f9.7)
  810 format (1x,'seed=',i9,'   axis=',i2,'   nav=',i7/1x, &
     'scal(1:k)=',6f9.4)
  820 format (1x,'tran(1:k)=',6f9.4)
  830 format (1x,'seed=',i9,'   scal=',f9.4,'   tran=',f9.4)
  840 format (1x,'nbf=',i7,'   nbfac=',i7,'   nfc=',i7,'   nface=',i7/ &
     1x,'nsmplx=',i7,'   sizht=',i7,'   nlo=',i5,'   tim=',f9.3)

  return
end
subroutine test_eqd ( )

!*****************************************************************************80
!
!! TEST_EQD
!
!     Driver program for testing routines EQDIS2, MFDEC2, TRISIZ,
!        INTPG, MMASEP, SEPMDF, SEPSHP, SFDWMF, SFUPMF.
!
  implicit none

  real ( kind = 8 ) d_pi
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) msglvl
  real ( kind = 8 ) tol

  integer ( kind = 4 ) incr,maxed,maxhv,maxiw,maxnc,maxnv,maxpv,maxvc,maxwk
  parameter (incr = 10000)
  parameter (maxed = 101)
  parameter (maxhv = 350)
  parameter (maxiw = 900)
  parameter (maxnc = 30)
  parameter (maxnv = 400)
  parameter (maxpv = 2000)
  parameter (maxvc = 800)
  parameter (maxwk = 1500)

  integer ( kind = 4 ) edge(4,maxed),ht(0:maxed-1)
  integer ( kind = 4 ) hvl(maxhv),icur(maxnc),ivrt(maxnv),iwk(maxiw),map(maxnc)
  integer ( kind = 4 ) ntri(maxhv),nvbc(maxnc),pvl(4,maxpv),regnum(maxhv)
  integer ( kind = 4 ) case,htsiz,i,imeas,irdr,j,ncur,nh,nhola,nhole,nmin,npolg
  integer ( kind = 4 ) nrfv,nsc,ntrid,ntrie,nv,nvc,nvcin,nvert,prime
  real    t0,initim,dectim,eqdtim
  real ( kind = 8 ) area(maxhv),h(maxhv),iang(maxpv),psi(maxhv)
  real ( kind = 8 ) vcl(2,maxvc),wk(maxwk)
  real ( kind = 8 ) angmin,angspc,angtol,dmin,kappa,umdf2
  logical hflag
  character rgname*20
  external umdf2
!
!  Read in vertices of general polygonal region.
!  CASE = 1 : simple polygon or multiply connected polygonal region
!  CASE = 2 : general polygonal region with holes and interfaces
!
  irdr = 5
  imeas = 7
  read (irdr,600) rgname
  read (irdr,*) tol, angspc,angtol,kappa,dmin,nmin,ntrid

  write (imeas,700) rgname,tol, angspc,angtol,kappa,dmin,nmin,ntrid
  angspc = angspc* d_pi ( ) /180.0D+00
  angtol = angtol* d_pi ( ) /180.0D+00
  hflag = (kappa >= 0.0D+00 .and. kappa <= 1.0D+00 )
  read (irdr,*) case,nvc,ncur,msglvl

  if (nvc > maxvc) then
    write ( *, 610) 'maxvc',nvc
    stop
  else if (ncur > maxnc) then
    write ( *, 610) 'maxnc',ncur
    stop
  end if

  read (irdr,*) (nvbc(i),i=1,ncur)
  if (case == 2) read (irdr,*) (icur(i),i=1,ncur)
  read (irdr,*) (vcl(1,i),vcl(2,i),i=1,nvc)
!
!     Call routine DSMCPR or DSPGDC to set data structures in arrays
!     REGNUM, HVL, PVL, IANG, HOLV = IWK.
!
  if (case == 1) then

    nhole = ncur - 1
    call gtime(t0)
    call dsmcpr(nhole,nvbc,vcl,maxhv,maxpv,maxiw,nvc,npolg,nvert, &
        nhola,regnum,hvl,pvl,iang,iwk,ierr)
    call gtime(initim)

  else if (case == 2) then

   nv = sum ( nvbc(1:ncur) )

   if (nv > maxnv) then
      write ( *, 610) 'maxnv',nv
      stop
   end if

   read (irdr,*) (ivrt(i),i=1,nv)
   nsc = 0
   do i = 1,nv
      if (ivrt(i) < 0) nsc = nsc + 1
   end do
   nsc = nsc/2
   if (nsc > maxed) then
      write ( *, 610) 'maxed',nsc
      stop
   end if
   htsiz = min(prime(nsc/2),maxed)
   call gtime(t0)
   call dspgdc(nvc,vcl,incr,ncur,nvbc,icur,ivrt,maxhv,maxpv,maxiw, &
        npolg,nvert,nhole,nhola,regnum,hvl,pvl,iang,iwk,htsiz,nsc, &
        ht,edge,map,ierr)

   call gtime(initim)
  end if
  initim = initim - t0
  if (ierr /= 0) then
   write ( *, 620) ierr
   write (imeas,620) ierr
   stop
  end if
!
!     Print out data structures and measurements after initialization.
!
  nvcin = nvc
  nh = nhole*2 + nhola
  write ( *, 670) msglvl
  write ( *, 630) nvc,(i,vcl(1,i),vcl(2,i),i=1,nvc)
  write ( *, 640) npolg,nhole,(hvl(i),i=1,npolg+nhole)
  write ( *, 650) (regnum(i),i=1,npolg)
  write ( *, 660) nvert,(i,(pvl(j,i),j=1,4),iang(i),i=1,nvert)
  write ( *, 640) nhola,nh,(iwk(i),i=1,nh)
  nrfv = 0
  angmin = 2.0D+00 * d_pi ( )
  do i = 1,nvert
    if (iang(i) > d_pi ( ) + tol) nrfv = nrfv + 1
    angmin = min(angmin,iang(i))
  end do

  angmin = angmin*180.0D+00/ d_pi ( )
  write (imeas,710) nvc,npolg,nvert,nhole,nhola,nrfv,angmin,initim
!
!     Obtain simple and convex polygon decompositions, and print
!     measurements.
!
  if (msglvl == 2) write ( *, 670)
  call gtime(t0)
!
!  Using SPDEC2 from old copy of GEOMPACK.
!
  call spdec2(angspc,angtol,nvc,npolg,nvert,nhole,nhola,maxvc,maxhv, &
     maxpv,maxiw-nh,maxwk,iwk,vcl,regnum,hvl,pvl,iang,iwk(nh+1),wk,ierr)
!
!  Using CVDEC2 from old copy of GEOMPACK.
!
  call cvdec2(angspc,angtol,nvc,npolg,nvert,maxvc,maxhv,maxpv,maxiw, &
     maxwk,vcl,regnum,hvl,pvl,iang,iwk,wk,ierr)

  call gtime(dectim)
  dectim = dectim - t0
  if (ierr /= 0) then
   write ( *, 620) ierr
   write (imeas,620) ierr
   stop
  end if
  angmin = 2.0D+00 * d_pi ( )
  do i = 1,nvert
   angmin = min(angmin,iang(i))
  end do
  angmin = angmin*180.0D+00/ d_pi ( )
  write (imeas,720) nvc,npolg,nvert,angmin,dectim
!
!     Obtain further convex polygon decomposition based on mesh
!     distribution function, and triangle sizes for the polygons.
!
  call gtime(t0)

  call eqdis2(hflag,umdf2,kappa,angspc,angtol,dmin,nmin,ntrid,nvc, &
     npolg,nvert,maxvc,maxhv,maxpv,maxiw,maxwk,vcl,regnum,hvl,pvl, &
     iang,area,psi,h,iwk,wk,ierr)

  call gtime(eqdtim)
  eqdtim = eqdtim - t0
  if (msglvl == 2) write ( *, 670) 0,0,0.0D+00,0.0D+00,0.0D+00,0.0D+00
  if (ierr /= 0) then
   write ( *, 620) ierr
   write (imeas,620) ierr
   stop
  end if
!
!     Print out data structures and measurements after decomposition.
!
  write ( *, 630) nvc,(i,vcl(1,i),vcl(2,i),i=nvcin+1,nvc)
  write ( *, 680) npolg,(hvl(i),i=1,npolg)
  write ( *, 650) (regnum(i),i=1,npolg)
  write ( *, 660) nvert,(i,(pvl(j,i),j=1,4),iang(i),i=1,nvert)
  angmin = 2.0D+00 * d_pi ( )
  do 50 i = 1,nvert
   angmin = min(angmin,iang(i))
   50 continue
  angmin = angmin*180.0D+00/ d_pi ( )
  ntrie = 0
  do 60 i = 1,npolg
   ntri(i) = int(ntrid*psi(i)*area(i) + 0.5D+00 )
   ntrie = ntrie + ntri(i)
   60 continue
  write ( *, 690) (i,area(i),psi(i),h(i),ntri(i),i=1,npolg)
  write (imeas,730) nvc,npolg,nvert,ntrie,angmin,eqdtim
!
  600 format (a20)
  610 format (1x,'*** ',a,' must be increased to',i8)
  620 format (/1x,'ierr=',i5)
  630 format (/1x,i7/(1x,i7,2f15.7))
  640 format (/1x,2i7/(1x,10i7))
  650 format (/(1x,10i7))
  660 format (/1x,i7/(1x,5i7,f15.7))
  670 format (1x,2i7,4f15.7)
  680 format (/1x,i7/(1x,10i7))
  690 format (/(1x,i7,3e15.7,i7))
  700 format (1x,a20/1x,'input : tol=',d15.7,'   angspc=',f9.3, &
     '   angtol=',f9.3/9x,'kappa=',f9.3,'   dmin=',f9.3, &
     '   nmin=',i5,'   ntrid=',i7)
  710 format (1x,'initds: nvc=',i7,'   npolg=',i7,'   nvert=',i7, &
     '   nhole=',i7/9x,'nhola=',i7,'   nrfv=',i7,'   angmin=',f9.3, &
     '   initim=',f9.3)
  720 format (1x,'decomp: nvc=',i7,'   npolg=',i7,'   nvert=',i7/ &
     9x,'angmin=',f9.3,'   dectim=',f9.3)
  730 format (1x,'eqdist: nvc=',i7,'   npolg=',i7,'   nvert=',i7, &
     '   ntrie=',i7/9x,'angmin=',f9.3,'   eqdtim=',f9.3)

  return
end
subroutine test_eqd3 ( )

!*****************************************************************************80
!
!! DREQD3
!
!     Driver program for testing routines EQDIS3,MFDEC3,INTPH,SEPFAC,
!        SFC1MF,SFC2MF,SFCSHP,TETSIZ.
!     Routines called:
!        CVDEC3,CVDECF,DSPHDC,DSPHFH,DSPHIH,EQDIS3,GTIM,PRIME
!
  implicit none

  real ( kind = 8 ) d_pi
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) msglvl
  real ( kind = 8 ) tol

  integer ( kind = 4 ) maxfp,maxfv,maxhf,maxht,maxiw,maxpf,maxvc,maxwk
  parameter (maxfp = 1000)
  parameter (maxfv = 5000)
  parameter (maxhf = 200)
  parameter (maxht = 307)
  parameter (maxiw = 12000)
  parameter (maxpf = 2000)
  parameter (maxvc = 600)
  parameter (maxwk = 20000)

  integer ( kind = 4 ) facep(3,maxfp),factyp(maxfp),fvl(6,maxfv),hfl(maxhf)
  integer ( kind = 4 ) ht(0:maxht-1),iwk(maxiw),pfl(2,maxpf)
  integer ( kind = 4 ) cnt,f,htsiz,i,ifach,imeas,ipolh,irdr,j,k,nbf,nface,nfach
  integer ( kind = 4 ) nfhol,nfph,nihol,nmin,npf,npolh,nrfe,ntetd,nvc,nvch,nvert
  integer ( kind = 4 ) nverth,outmod,prime
  real    dectim,eqdtim,holtim,initim,t0,t1
  real ( kind = 8 ) eang(maxfv),h(maxhf),nrml(3,maxfp),psi(maxhf)
  real ( kind = 8 ) vcl(3,maxvc),vol(maxhf),wk(maxwk)
  real ( kind = 8 ) angacc,angedg,angmin,aspc2d,atol2d,dmin,kappa
  real ( kind = 8 ) rdacc, umdf3
  logical hflag,hoflag,holint,nsflag
  character rgname*60
  external umdf3
!
!     OUTMOD = 1 : Output polyhedral decomposition data structure.
!     OUTMOD = 2 : Output VCL + list of vertex indices for each face.
!     MSGLVL >= 2: Routine DSPHFH prints out non-coplanar vertices of
!        each face (also output from routines CUTFAC, INSFAC).
!     NFHOL <= 0: Simpler routine DSPHDC is called instead of DSPHFH.
!     0 <= KAPPA <= 1: Set HFLAG to TRUE; NMIN < 0: Set NSFLAG to TRUE.
!
  irdr = 5
  imeas = 7
  read (irdr,600) rgname
  read (irdr,*) tol, aspc2d,atol2d,angacc,rdacc,angedg,kappa,dmin, &
     nmin,ntetd

  read (irdr,*) nvc,nface,nfhol,npolh,nihol,outmod,msglvl
  write (imeas,800) rgname,tol, aspc2d,atol2d,angacc,rdacc,nfhol, &
     nihol,angedg,kappa,dmin,nmin,ntetd
  hoflag = (nfhol > 0)
  nfhol = max(nfhol,0)
  hflag = (kappa >= 0.0D+00 .and. kappa <= 1.0D+00 )
  nsflag = (nmin < 0)
  nmin = abs(nmin)
  aspc2d = aspc2d* d_pi ( ) /180.0D+00
  atol2d = atol2d* d_pi ( ) /180.0D+00
  angacc = angacc* d_pi ( ) /180.0D+00
  angedg = angedg* d_pi ( ) /180.0D+00
  nfph = nface + nfhol
  if (nvc > maxvc .or. nfph >= maxfp .or. npolh >= maxhf) &
  then
   if (nvc > maxvc) write ( *, 610) 'maxvc',nvc
   if (nfph >= maxfp) write ( *, 610) 'maxfp',nfph+1
   if (npolh >= maxhf) write ( *, 610) 'maxhf',npolh+1
   stop
  end if
  read (irdr,*) ((vcl(j,i),j=1,3),i=1,nvc)
!
!     The outer polygons of NFACE faces should appear first, followed by
!     the inner polygons of NFHOL holes (in order of index of faces
!     containing holes). Face types for NFACE faces should be positive.
!     Face types for NFHOL holes are +F or -F where F is index of face
!     containing hole; positive (negative) sign indicates hole polygon
!     is oriented CCW (CW) in polyhedron when viewed from outside.
!     Holes must only be on boundary faces of polyhedral region.
!
  read (irdr,*) (facep(1,i),i=1,nfph+1)
  nvert = facep(1,nfph+1) - 1
  read (irdr,*) (factyp(i),i=1,nfph)
  if (nvert > maxfv) then
   write ( *, 610) 'maxfv',nvert
   stop
  end if
!
!     Head vertex of each (outer) face must be a strictly convex vertex.
!     Head vertex of each hole may be arbitrary.
!     Positive (negative) sign in PFL(1,*) indicates face is oriented
!     CCW (CW) in polyhedron when viewed from outside.
!     Hole polygons should not be included in PFL.
!
  read (irdr,*) (fvl(1,i),i=1,nvert)
  read (irdr,*) (hfl(i),i=1,npolh+1)
  npf = hfl(npolh+1) - 1
  if (npf + nfhol > maxpf) then
   write ( *, 610) 'maxpf',npf+nfhol
   stop
  end if
  read (irdr,*) (pfl(1,i),i=1,npf)
  htsiz = min(prime(nvc+2),maxht)
  call gtime(t0)

  if (hoflag) then
   call dsphfh(aspc2d,atol2d,nvc,nface,nfhol,npolh,maxvc,maxfv, &
        maxiw,maxwk,nvert,npf,vcl,facep,factyp,nrml,fvl,eang,hfl, &
        pfl,htsiz,ht,iwk,wk,ierr)
  else
     call dsphdc(nvc,nface,npolh,vcl,facep,nrml,fvl,eang,hfl,pfl, &
        htsiz,maxiw/4,iwk,ht,ierr)
  end if

  call gtime(initim)
  initim = initim - t0
  if (ierr /= 0) then
   write ( *, 620) ierr
   write (imeas,620) ierr
   stop
  end if
  nrfe = 0
  angmin = 2.0D+00 * d_pi ( )

  do i = 1,nvert
    if (eang(i) > d_pi ( ) + tol) nrfe = nrfe + 1
    if (eang(i) > -1.0D+00 ) angmin = min(angmin,eang(i))
  end do

  angmin = angmin*180.0D+00/ d_pi ( )
  if (hoflag) then
     write (imeas,810) 'fholds',nvc,nface,nvert,npolh,npf,nrfe, &
        angmin,'initim',initim
  else
     write (imeas,810) 'initds',nvc,nface,nvert,npolh,npf,nrfe, &
        angmin,'initim',initim
  end if
!
!     Read in description of interior hole polyhedra. The input format
!     is similar to above with no hole faces and only 1 polyhedron at
!     a time, treated as though the region for each polyh is its hole.
!     IPOLH is index of polyhedron containing hole. IFACH is index of
!     `extreme' face of hole used for connection to outer boundary.
!     One interior hole polyhedron per polyhedral region is assumed.
!     A negative sign for IPOLH indicates hole polyh is hole interface.
!
  holtim = 0.0
  do 20 k = 1,nihol
     read (irdr,*) nvch,nfach,ipolh,ifach
   holint = (ipolh < 0)
   ipolh = abs(ipolh)
     if (nvc+nvch > maxvc .or. nface+nfach >= maxfp .or. &
     npf+nfach > maxpf) then
      if (nvc+nvch > maxvc) write ( *, 610) 'maxvc',nvc+nvch
      if (nface+nfach >= maxfp) write ( *, 610) 'maxfp', &
           nface+nfach+1
      if (npf+nfach > maxpf) write ( *, 610) 'maxpf',npf+nfach
      stop
     end if
     read (irdr,*) ((vcl(j,i),j=1,3),i=nvc+1,nvc+nvch)
     read (irdr,*) (facep(1,i),i=nface+1,nface+nfach+1)
     nverth = facep(1,nface+nfach+1) - 1
     read (irdr,*) (factyp(i),i=nface+1,nface+nfach)
     if (nvert+nverth > maxfv) then
      write ( *, 610) 'maxfv',nvert+nverth
      stop
     end if
   read (irdr,*) (fvl(1,i),i=nvert+1,nvert+nverth)
     read (irdr,*) (pfl(1,i),i=npf+1,npf+nfach)
     htsiz = min(prime(nvch+2),maxht)
     call gtime(t0)
     call dsphih(aspc2d,atol2d,angacc,rdacc,nvc,nface,nvert,npolh, &
        npf,nvch,nfach,ipolh,ifach,holint,maxvc,maxfp,maxfv,maxhf, &
        maxpf,maxiw,maxwk,vcl,facep,factyp,nrml,fvl,eang,hfl,pfl, &
        htsiz,ht,iwk,wk,ierr)

     call gtime(t1)
     holtim = holtim + (t1 - t0)
     if (ierr /= 0) then
      write ( *, 620) ierr
      write (imeas,620) ierr
      stop
     end if
   20 continue
  if (nihol > 0) then
     nrfe = 0
     angmin = 2.0D+00 * d_pi ( )
     do 30 i = 1,nvert
      if (eang(i) > d_pi ( ) + tol) nrfe = nrfe + 1
      if (eang(i) > -1.0D+00 ) angmin = min(angmin,eang(i))
   30    continue
     angmin = angmin*180.0D+00/ d_pi ( )
     write (imeas,810) 'iholds',nvc,nface,nvert,npolh,npf,nrfe, &
        angmin,'holtim',holtim
  end if
!
!     Decompose polyhedral region into convex parts.
!
  cnt = 1
   40 continue
  call gtime(t0)
  call cvdec3(angacc,rdacc,nvc,nface,nvert,npolh,npf,maxvc,maxfp, &
     maxfv,maxhf,maxpf,maxiw,maxwk,vcl,facep,factyp,nrml,fvl,eang, &
     hfl,pfl,iwk,wk,ierr)
  call gtime(dectim)
  dectim = dectim - t0
  if (ierr /= 0) then
   if (ierr /= 327 .or. cnt >= 3) then
      write ( *, 620) ierr
      write (imeas,620) ierr
      if (ierr /= 327) stop
   end if
  end if
  nrfe = 0
  angmin = 2.0D+00 * d_pi ( )

  do i = 1,nvert
    if (eang(i) > d_pi ( ) + tol) nrfe = nrfe + 1
    if (eang(i) > -1.0D+00 ) angmin = min(angmin,eang(i))
  end do

  angmin = angmin*180.0D+00/ d_pi ( )
  write (imeas,820) nvc,nface,nvert,npolh,npf,nrfe,angmin,dectim
  if (ierr == 327 .and. cnt < 3) then
   angacc = angacc - d_pi ( ) /36.0D+00
   if (angacc > 0.0D+00) then
      rdacc = rdacc*0.95D+00
      ierr = 0
        cnt = cnt + 1
      go to 40
   end if
  else if (ierr == 327) then
   stop
  end if
!
!     Decompose faces of polyhedral region into convex subpolygons.
!
  call gtime(t0)
  call cvdecf(aspc2d,atol2d,nvc,nface,nvert,npf,maxvc,maxfp,maxfv, &
     maxpf,maxiw,maxwk,vcl,facep,factyp,nrml,fvl,eang,hfl,pfl,iwk, &
     wk,ierr)
  call gtime(dectim)
  dectim = dectim - t0
  if (ierr /= 0) then
   write ( *, 620) ierr
   write (imeas,620) ierr
   stop
  end if
  write (imeas,830) nvc,nface,nvert,npolh,npf,dectim
!
!     Further subdivide convex polyhedra based on mesh distribution
!     function, and determine tetrahedron size for each polyhedron.
!
  call gtime(t0)
  call eqdis3(hflag,umdf3,kappa,angacc,angedg,dmin,nmin,ntetd, &
     nsflag,nvc,nface,nvert,npolh,npf,maxvc,maxfp,maxfv,maxhf,maxpf, &
     maxiw,maxwk,vcl,facep,factyp,nrml,fvl,eang,hfl,pfl,vol,psi,h, &
     iwk,wk,ierr)

  call gtime(eqdtim)
  eqdtim = eqdtim - t0
  if (ierr /= 0) then
   write ( *, 620) ierr
   write (imeas,620) ierr
   stop
  end if
  angmin = 2.0D+00 * d_pi ( )
  do i = 1,nvert
    if (eang(i) > -1.0D+00 ) angmin = min(angmin,eang(i))
  end do
  angmin = angmin*180.0D+00/ d_pi ( )
  write (imeas,840) nvc,nface,nvert,npolh,npf,angmin,eqdtim

  if (outmod == 1) then
     write ( *, 630) nvc,nface,nvert,npolh,npf
     write ( *, 640) (i,(vcl(j,i),j=1,3),i=1,nvc)
     write ( *, 650) (i,(facep(j,i),j=1,3),factyp(i), &
        (nrml(j,i),j=1,3),i=1,nface)
     write ( *, 660) (i,(fvl(j,i),j=1,6),eang(i),i=1,nvert)
     write ( *, 670) (i,hfl(i),vol(i),psi(i),h(i),i=1,npolh)
     write ( *, 680) (i,(pfl(j,i),j=1,2),i=1,npf)
  else
   nbf = 0
   do 70 f = 1,nface
      if (facep(3,f) == 0) nbf = nbf + 1
   70    continue
   write ( *, 690) nvc,nbf,nface-nbf
   write ( *, 700) ((vcl(j,i),j=1,3),i=1,nvc)
   do 100 i = 1,2
        write ( *, '(a)' )
      do 90 f = 1,nface
         if (i == 1 .and. facep(3,f) /= 0) go to 90
         if (i == 2 .and. facep(3,f) == 0) go to 90
         k = 0
         j = facep(1,f)
   80          continue
      k = k + 1
      iwk(k) = fvl(1,j)
      j = fvl(3,j)
         if (j /= facep(1,f)) go to 80
         if (k > maxiw) then
            write ( *, 610) 'maxiw',k
            stop
         end if
         write ( *, 710) k,(iwk(j),j=1,k)
   90       continue
  100    continue
  end if

  600 format (a60)
  610 format (1x,'*** ',a,' must be increased to',i8)
  620 format (/1x,'ierr=',i5)
  630 format (1x,'nvc,nface,nvert,npolh,npf'/1x,5i7)
  640 format (/' vcl'/(1x,i7,3f15.7))
  650 format (/' facep,factyp,nrml'/(1x,4i7,i5,3f15.7))
  660 format (/' fvl,eang'/(1x,7i7,f15.7))
  670 format (/' hfl,vol,psi,h'/(1x,2i7,3f15.7))
  680 format (/' pfl'/(1x,3i7))
  690 format (1x,3i7)
  700 format (/(1x,3f15.7))
  710 format (1x,11i7/(8x,10i7))
  800 format (1x,a60/1x,'input : tol=',d15.7,'   aspc2d=',f9.3, &
     '   atol2d=',f9.3/9x,'angacc=',f9.3,'   rdacc=',f9.5, &
     '   nfhol=',i3,'   nihol=',i3/9x,'angedg=',f9.3,'   kappa=', &
     f7.3,'   dmin=',f6.3,'   nmin=',i5/9x,'ntetd=',i7)
  810 format (1x,a6,': nvc=',i7,'   nface=',i7,'   nvert=',i7, &
     '   npolh=',i7/9x,'npf=',i7,'   nrfe=',i7,'   angmin=',f9.3, &
     3x,a6,'=',f9.3)
  820 format (1x,'decomp: nvc=',i7,'   nface=',i7,'   nvert=',i7, &
     '   npolh=',i7/9x,'npf=',i7,'   nrfe=',i7,'   angmin=',f9.3, &
     '   dectim=',f9.3)
  830 format (1x,'decfac: nvc=',i7,'   nface=',i7,'   nvert=',i7, &
     '   npolh=',i7/9x,'npf=',i7,'   dectim=',f9.3)
  840 format (1x,'eqdist: nvc=',i7,'   nface=',i7,'   nvert=',i7, &
     '   npolh=',i7/9x,'npf=',i7,'   angmin=',f9.3, &
     '   eqdtim=',f9.3)

  return
end
subroutine test_fchk ( )

!*****************************************************************************80
!
!! TEST_FCHK
!
!     Driver program for testing routines ORDERK, AVAILK, HTINSK,
!        HTDELK, HTSRCK, HTSDLK, UPDATK.
!     Routines called:
!        AVAILK, HTDELK, HTINSK, HTSDLK, HTSRCK, UPDATK
!
  implicit none

  integer ( kind = 4 ) kmax,maxfc,n,p
  parameter (kmax = 10)
  parameter (maxfc = 10)
  parameter (n = 10)
  parameter (p = 7)

  integer ( kind = 4 ) ierr

  integer ( kind = 4 ) back,case,d,e,front,hdavfc,htsrck,i,irdr,j,k,nfc,pos
  integer ( kind = 4 ) fc((kmax+4)*maxfc),ht(0:p-1),ind(kmax)

  irdr = 5

  write ( *, '(a)' ) 'enter k:'
  read (irdr,*) k
  if (k <= 1 .or. k > kmax) stop

  ht(0:p-1) = 0

  front = 0
  back = 0
  hdavfc = 0
  nfc = 0
   20 continue
   write ( *, '(a)' ) 'enter case:'
   read (irdr,*) case
   if (case < 1 .or. case > 6) stop
   if (case == 1) then
       write ( *, '(a)' ) 'htinsk - enter ind(1:k),d,e:'
       read (irdr,*) (ind(i),i=1,k),d,e
       call availk(k,hdavfc,nfc,maxfc,fc,pos,ierr)
      if (ierr /= 0) then
         write ( *, '(a)' ) '*** array fc is full'
         ierr = 0
      else
         call htinsk(k,pos,ind,d,e,n,p,fc,ht)
      end if
   else if (case == 2) then
       write ( *, '(a)' ) 'htdelk - enter pos:'
       read (irdr,*) pos
      if (pos <= 0 .or. pos > nfc .or. fc((pos-1)*(k+4)+1) &
        <= 0) then
         write ( *, '(a)' ) '*** invalid index'
      else
         call htdelk(k,pos,n,p,fc,ht)
         fc((pos-1)*(k+4)+1) = -hdavfc
         hdavfc = pos
      end if
   else if (case == 3) then
       write ( *, '(a)' ) 'htsrck - enter ind(1:k):'
       read (irdr,*) (ind(i),i=1,k)
      pos = htsrck(k,ind,n,p,fc,ht)
       write ( *, '(a)' ) 'pos=',pos
   else if (case == 4) then
       write ( *, '(a)' ) 'htsdlk - enter ind(1:k):'
       read (irdr,*) (ind(i),i=1,k)
      call htsdlk(k,ind,n,p,fc,ht,pos)
       write ( *, '(a)' ) 'pos=',pos
      if (pos > 0) then
         fc((pos-1)*(k+4)+1) = -hdavfc
         hdavfc = pos
      end if
   else if (case == 5) then
       write ( *, '(a)' ) 'updatk - enter ind(1:k),d,e,i:'
       read (irdr,*) (ind(j),j=1,k),d,e,i
      call updatk(k,ind,d,e,i,n,p,front,back,fc,ht, ierr )
      if (ierr /= 0) then
         write ( *, '(a)' ) '*** unsuccessful search'
         ierr = 0
      end if
   else
      write ( *, 600) hdavfc,nfc,front,back,(ht(i),i=0,p-1)
      write ( *, 610) 'fc'
      do i = 1,nfc
         write ( *, 620) i,(fc((i-1)*(k+4)+j),j=1,k+4)
      end do
   end if
  go to 20
  600 format (' hdavfc=',i3,'   nfc=',i3,'   front=',i3,'   back=',i3/ &
     ' ht:  ',7i5)
  610 format (1x,a)
  620 format (1x,15i5)

end
subroutine test_fcht ( )

!*****************************************************************************80
!
!! TEST_FCHT tests ORDER3, AVAILF, HTINS, HTDEL, HTSRC, UPDATF.
!
  implicit none

  integer ( kind = 4 ) maxfc,n,p
  parameter (maxfc = 10)
  parameter (n = 10)
  parameter (p = 7)

  integer ( kind = 4 ) ierr

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b,back,c,case,d,e,front,hdavfc,htsrc,i,ind,irdr,nfc
  integer ( kind = 4 ) fc(7,maxfc),ht(0:p-1)

  irdr = 5

  ht(0:p-1) = 0
  front = 0
  back = 0
  hdavfc = 0
  nfc = 0
   20 continue
   write ( *, '(a)' ) 'enter case:'
   read (irdr,*) case
   if (case < 1 .or. case > 5) stop
   if (case == 1) then
       write ( *, '(a)' ) 'htins - enter a,b,c,d,e:'
       read (irdr,*) a,b,c,d,e
       call availf(hdavfc,nfc,maxfc,fc,ind,ierr)
      if (ierr /= 0) then
         write ( *, '(a)' ) '*** array fc is full'
         ierr = 0
      else
         call htins(ind,a,b,c,d,e,n,p,fc,ht)
      end if
   else if (case == 2) then
       write ( *, '(a)' ) 'htdel - enter ind:'
       read (irdr,*) ind
      if (ind <= 0 .or. ind > nfc .or. fc(1,ind) <= 0) &
        then
         write ( *, '(a)' ) '*** invalid index'
      else
         call htdel(ind,n,p,fc,ht)
         fc(1,ind) = -hdavfc
         hdavfc = ind
      end if
   else if (case == 3) then
       write ( *, '(a)' ) 'htsrc - enter a,b,c:'
       read (irdr,*) a,b,c
      ind = htsrc(a,b,c,n,p,fc,ht)
       write ( *, '(a)' ) 'ind=',ind
   else if (case == 4) then
       write ( *, '(a)' ) 'updatf - enter a,b,c,d,e,i:'
       read (irdr,*) a,b,c,d,e,i
      call updatf(a,b,c,d,e,i,n,p,front,back,fc,ht, ierr )
      if (ierr /= 0) then
         write ( *, '(a)' ) '*** unsuccessful search'
         ierr = 0
      end if
   else
      write ( *, 600) hdavfc,nfc,front,back,(ht(i),i=0,p-1), &
           (i,(fc(a,i),a=1,7),i=1,nfc)
   end if

  go to 20

  600 format (' hdavfc=',i3,'   nfc=',i3,'   front=',i3,'   back=',i3/ &
     ' ht:  ',7i5,5x,'fc:'/(1x,8i5))
end
subroutine test_line ( )

!*****************************************************************************80
!
!! TEST_LINE tests LRLINE, XLINE, XEDGE, CMCIRC, DIAEDG.
!
  implicit none

  integer ( kind = 4 ) cmcirc,diaedg,lrline
  integer ( kind = 4 ) case,in,irdr,lr,mode
  real ( kind = 8 ) d12,d34,x0,x1,x2,x3,x4,y0,y1,y2,y3,y4
  logical intsct,parall

  irdr = 5

10 continue

   write ( *, '(a)' ) 'enter case:'
   read (irdr,*) case
   if (case < 1 .or. case > 5) stop
   if (case == 1) then
      write ( *, '(a)' ) 'lrline - enter x0,y0,x1,y1,x2,y2,d12:'
      read (irdr,*) x0,y0,x1,y1,x2,y2,d12
      lr = lrline(x0,y0,x1,y1,x2,y2,d12)
      write ( *, '(a)' ) 'lr=',lr
   else if (case == 2) then
      write ( *, '(a)' ) 'xline - enter x1,y1,x2,y2,x3,y3,x4,y4,', &
           'd12,d34:'
      read (irdr,*) x1,y1,x2,y2,x3,y3,x4,y4,d12,d34
      x0 = 0.0D+00
      y0 = 0.0D+00
      call xline(x1,y1,x2,y2,x3,y3,x4,y4,d12,d34,x0,y0,parall)
      write ( *, '(a)' ) 'parall=',parall,'   x0=',x0,'   y0=',y0
   else if (case == 3) then
      write ( *, '(a)' ) 'xedge - enter mode,x1,y1,x2,y2,x3,y3,x4,y4'
      read (irdr,*) mode,x1,y1,x2,y2,x3,y3,x4,y4
      x0 = 0.0D+00
      y0 = 0.0D+00
      call xedge(mode,x1,y1,x2,y2,x3,y3,x4,y4,x0,y0,intsct)
      write ( *, '(a)' ) 'intsct=',intsct,'   x0=',x0,'   y0=',y0
   else if (case == 4) then
      write ( *, '(a)' ) 'cmcirc - enter x0,y0,x1,y1,x2,y2,x3,y3'
      read (irdr,*) x0,y0,x1,y1,x2,y2,x3,y3
      in = cmcirc(x0,y0,x1,y1,x2,y2,x3,y3)
      write ( *, '(a)' ) 'in=',in
   else
      write ( *, '(a)' ) 'diaedg - enter x0,y0,x1,y1,x2,y2,x3,y3'
      read (irdr,*) x0,y0,x1,y1,x2,y2,x3,y3
      in = diaedg(x0,y0,x1,y1,x2,y2,x3,y3)
      write ( *, '(a)' ) 'in=',in
   end if
  go to 10
end
subroutine test_lu ( )

!*****************************************************************************80
!
!! TEST_LU tests LUFAC and LUSOL.
!
  implicit none

  integer ( kind = 4 ), parameter :: nmax = 200

  real ( kind = 8 ) a(nmax,nmax)
  real ( kind = 8 ) b(nmax)
  real ( kind = 8 ) emax
  real ( kind = 8 ) esum
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ipvt(nmax)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  integer ( kind = 4 ) seed
  logical singlr
  real ( kind = 8 ) t
  real    t0
  real    tf
  real ( kind = 8 ) tol
  real    ts
  real ( kind = 8 ) urand

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_LU'
  write ( *, '(a)' ) '  LUFAC computes the LU factorization of a matrix.'
  write ( *, '(a)' ) '  LUSOL solves a linear system.'

  n = 4
  tol = 0.000001D+00
  seed = 1952

  do j = 1, n
    do i = 1, n
      a(i,j) = urand ( seed ) * 2.0D+00 - 1.0D+00
    end do
  end do

  do i = 1,n
    b(i) = sum ( a(i,1:n) )
  end do

  if ( n <= 4 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Matrix A and right hand side b:'
    do i = 1,n
      write ( *, 610) (a(i,j),j=1,n),b(i)
    end do
  end if

  call gtime(t0)

  call lufac ( a, nmax, n, tol, ipvt, singlr )
  call gtime(tf)
  tf = tf - t0

  if ( singlr ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The matrix is singular'
    return
  end if

  call gtime(t0)
  call lusol ( a, nmax, n, ipvt, b )
  call gtime(ts)
  ts = ts - t0

  if ( n <= 4 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  ipvt, lu'
    ipvt(n) = n
    do i = 1,n
      write ( *, 620) ipvt(i),(a(i,j),j=1,n)
    end do
  end if

  emax = 0.0D+00
  esum = 0.0D+00
  do i = 1, n
    t = abs ( b(i) - 1.0D+00 )
    emax = max(emax,t)
    esum = esum + t
  end do

  write ( *, '(a)' ) ' '
  write ( *, 630) (b(i),i=1,min(4,n))
  write ( *, 640) emax, esum, tf, ts

  600 format (1x,a)
  610 format (1x,5f15.7)
  620 format (1x,i5,4f15.7)
  630 format (1x,'x = ',4f15.7)
  640 format (1x,'emax,esum = ',2e15.7/1x,'tf, ts = ',2f9.3)

  return
end
subroutine test_mdf ( )

!*****************************************************************************80
!
!! TEST_MDF tests DSMDF2, MDF2 and PRMDF2.
!
  implicit none

  real ( kind = 8 ) d_pi
  integer ( kind = 4 ) ierr

  integer ( kind = 4 ) incr,maxed,maxhv,maxiw,maxnc,maxpv,maxvc,maxwk
  parameter (incr = 10000)
  parameter (maxed = 101)
  parameter (maxhv = 200)
  parameter (maxiw = 900)
  parameter (maxnc = 30)
  parameter (maxpv = 1000)
  parameter (maxvc = 500)
  parameter (maxwk = 1500)

  integer ( kind = 4 ) edge(4,maxed),ht(0:maxed-1),hvl(maxhv)
  integer ( kind = 4 ) icur(maxnc),ivrt(maxpv),iwk(maxiw),map(maxnc),nev(maxhv)
  integer ( kind = 4 ) nvbc(maxnc),pvl(4,maxpv),regnum(maxhv),xivrt(maxhv+1)
  integer ( kind = 4 ) case,htsiz,i,ifv,irdr,j,ncur,nh,nhola,nhole
  integer ( kind = 4 ) npolg,nsc,nv,nvc,nvert,prime
  real ( kind = 8 ) area(maxhv),edgval(maxpv),iang(maxpv),val(maxhv)
  real ( kind = 8 ) vcl(2,maxvc),vrtval(maxvc),widsq(maxhv),wk(maxwk)
  real ( kind = 8 ) angspc,angtol,mdf2,tol,x,y
  character rgname*20
!
!     Read in vertices of general polygonal region.
!     CASE = 1 : simple polygon or multiply connected polygonal region
!     CASE = 2 : general polygonal region with holes and interfaces
!
  irdr = 5
  read (irdr,600) rgname
  read (irdr,*) tol, angspc, angtol

  angspc = angspc* d_pi ( ) /180.0D+00
  angtol = angtol* d_pi ( ) /180.0D+00
  read (irdr,*) case,nvc,ncur
  if (nvc > maxvc) then
   write ( *, 610) 'maxvc',nvc
   stop
  else if (ncur > maxnc) then
   write ( *, 610) 'maxnc',ncur
   stop
  end if
  read (irdr,*) (nvbc(i),i=1,ncur)
  if (case == 2) read (irdr,*) (icur(i),i=1,ncur)
  read (irdr,*) (vcl(1,i),vcl(2,i),i=1,nvc)
!
!     Call routine DSMCPR or DSPGDC to set data structures in arrays
!     REGNUM, HVL, PVL, IANG, HOLV = IWK.
!
  if (case == 1) then
   nhole = ncur - 1
   call dsmcpr(nhole,nvbc,vcl,maxhv,maxpv,maxiw,nvc,npolg,nvert, &
        nhola,regnum,hvl,pvl,iang,iwk,ierr)
  else if (case == 2) then

   nv = sum ( nvbc(1:ncur) )

   if (nv > maxpv) then
      write ( *, 610) 'maxpv',nv
      stop
   end if
   read (irdr,*) (ivrt(i),i=1,nv)
   nsc = 0

   do i = 1,nv
     if (ivrt(i) < 0) nsc = nsc + 1
   end do

   nsc = nsc/2
   if (nsc > maxed) then
      write ( *, 610) 'maxed',nsc
      stop
   end if
   htsiz = min(prime(nsc/2),maxed)
   call dspgdc(nvc,vcl,incr,ncur,nvbc,icur,ivrt,maxhv,maxpv,maxiw, &
        npolg,nvert,nhole,nhola,regnum,hvl,pvl,iang,iwk,htsiz,nsc, &
        ht,edge,map,ierr)

  end if
  if (ierr /= 0) then
   write ( *, 620) ierr
   stop
  end if
!
!     Obtain simple and convex polygon decompositions.
!
  nh = nhole*2 + nhola
  call spdec2(angspc,angtol,nvc,npolg,nvert,nhole,nhola,maxvc,maxhv, &
     maxpv,maxiw-nh,maxwk,iwk,vcl,regnum,hvl,pvl,iang,iwk(nh+1),wk)
  call cvdec2(angspc,angtol,nvc,npolg,nvert,maxvc,maxhv,maxpv,maxiw, &
     maxwk,vcl,regnum,hvl,pvl,iang,iwk,wk)
  if (ierr /= 0) then
   write ( *, 620) ierr
   stop
  end if
!
!     Initialize data structure for heuristic mesh distribution function
!     and evaluate function at centroid of each convex polygon.
!
  call dsmdf2(.true.,nvc,npolg,maxwk,vcl,hvl,pvl,iang,ivrt,xivrt, &
     widsq,edgval,vrtval,area,wk,ierr)

  if (ierr /= 0) then
   write ( *, 620) ierr
   stop
  end if
  do 40 i = 1,npolg
   call prmdf2(i,widsq(i),ivrt,xivrt,edgval,vrtval,nev(i),ifv,iwk)
   if (nev(i) == 0) then
      val(i) = widsq(i)
   else
      nv = xivrt(i+1) - xivrt(i)
      x = 0.0D+00
      y = 0.0D+00
      do 30 j = xivrt(i),xivrt(i+1)-1
         x = x + vcl(1,ivrt(j))
         y = y + vcl(2,ivrt(j))
   30       continue
      val(i) = 1.0D+00 /mdf2(x/nv,y/nv,widsq(i),nev(i),ifv,iwk,ivrt, &
           edgval,vrtval,vcl)
   end if
   40 continue
!
!     Print arrays from calls to 3 mdf routines.
!
  write ( *, 630) npolg,(i,xivrt(i),area(i),widsq(i),val(i),nev(i), &
     i=1,npolg)
  write ( *, 640) nvert,nvc,(i,ivrt(i),edgval(i),vrtval(i),i=1,nvc)
  write ( *, 650) (i,ivrt(i),edgval(i),i=nvc+1,nvert)

  600 format (a20)
  610 format (1x,'*** ',a,' must be increased to',i8)
  620 format (/1x,'ierr=',i5)
  630 format (/1x,i7/(1x,2i7,3f15.7,i7))
  640 format (/1x,2i7/(1x,2i7,2f15.7))
  650 format (1x,2i7,f15.7)

  return
end
subroutine test_mdf3 ( )

!*****************************************************************************80
!
!! DRMDF3 tests DSMDF3, MDF3 and PRMDF3.
!
  implicit none

  real ( kind = 8 ) d_pi
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) msglvl
  real ( kind = 8 ) tol

  integer ( kind = 4 ) maxfp,maxfv,maxhf,maxht,maxiw,maxpf,maxvc,maxwk
  parameter (maxfp = 800)
  parameter (maxfv = 3500)
  parameter (maxhf = 200)
  parameter (maxht = 307)
  parameter (maxiw = 2500)
  parameter (maxpf = 1200)
  parameter (maxvc = 600)
  parameter (maxwk = 1000)

  integer ( kind = 4 ) facep(3,maxfp),factyp(maxfp),fvl(6,maxfv),hfl(maxhf)
  integer ( kind = 4 ) ht(0:maxht-1),ifac(maxpf),ivrt(maxfv),iwk(maxiw)
  integer ( kind = 4 ) pfl(2,maxpf),xifac(maxhf+1),xivrt(maxfp+1)
  integer ( kind = 4 ) cnt,f,htsiz,i,ifach,imeas,ipolh,irdr,j,k,n,nbf,ncedge
  integer ( kind = 4 ) ncface,nedev,nface,nfach,nfcev,nfhol,nfph,nihol,npf,npolh
  integer ( kind = 4 ) nrfe,nvc,nvch,nvert,nverth,nvrev,outmod,prime,sume,sumf
  integer ( kind = 4 ) sumv
  real    dectim,holtim,initim,mdftim,t0,t1
  real ( kind = 8 ) eang(maxfv),edgval(maxfv),facval(maxfp)
  real ( kind = 8 ) nrml(3,maxfp),val(maxhf),vcl(3,maxvc)
  real ( kind = 8 ) vrtval(maxvc),wid(maxhf),wk(maxwk)
  real ( kind = 8 ) angacc,angmin,aspc2d,atol2d,mdf3,rdacc
  real ( kind = 8 ) x,y,z
  logical hoflag,holint
  character rgname*60
!
!     OUTMOD = 1 : Output polyhedral decomposition data structure.
!     OUTMOD = 2 : Output VCL + list of vertex indices for each face.
!     MSGLVL >= 2: Routine DSPHFH prints out non-coplanar vertices of
!        each face (also output from routines CUTFAC, INSFAC).
!     NFHOL <= 0: Simpler routine DSPHDC is called instead of DSPHFH.
!
  irdr = 5
  imeas = 7
  read (irdr,600) rgname
  read (irdr,*) tol, aspc2d,atol2d,angacc,rdacc

  read (irdr,*) nvc,nface,nfhol,npolh,nihol,outmod,msglvl
  write (imeas,800) rgname,tol, aspc2d,atol2d,angacc,rdacc,nfhol, &
     nihol
  hoflag = (nfhol > 0)
  nfhol = max(nfhol,0)
  aspc2d = aspc2d* d_pi ( ) /180.0D+00
  atol2d = atol2d* d_pi ( ) /180.0D+00
  angacc = angacc* d_pi ( ) /180.0D+00
  nfph = nface + nfhol
  if (nvc > maxvc .or. nfph >= maxfp .or. npolh >= maxhf) &
  then
   if (nvc > maxvc) write ( *, 610) 'maxvc',nvc
   if (nfph >= maxfp) write ( *, 610) 'maxfp',nfph+1
   if (npolh >= maxhf) write ( *, 610) 'maxhf',npolh+1
   stop
  end if
  read (irdr,*) ((vcl(j,i),j=1,3),i=1,nvc)
!
!     The outer polygons of NFACE faces should appear first, followed by
!     the inner polygons of NFHOL holes (in order of index of faces
!     containing holes). Face types for NFACE faces should be positive.
!     Face types for NFHOL holes are +F or -F where F is index of face
!     containing hole; positive (negative) sign indicates hole polygon
!     is oriented CCW (CW) in polyhedron when viewed from outside.
!     Holes must only be on boundary faces of polyhedral region.
!
  read (irdr,*) (facep(1,i),i=1,nfph+1)
  nvert = facep(1,nfph+1) - 1
  read (irdr,*) (factyp(i),i=1,nfph)
  if (nvert > maxfv) then
   write ( *, 610) 'maxfv',nvert
   stop
  end if
!
!     Head vertex of each (outer) face must be a strictly convex vertex.
!     Head vertex of each hole may be arbitrary.
!     Positive (negative) sign in PFL(1,*) indicates face is oriented
!     CCW (CW) in polyhedron when viewed from outside.
!     Hole polygons should not be included in PFL.
!
  read (irdr,*) (fvl(1,i),i=1,nvert)
  read (irdr,*) (hfl(i),i=1,npolh+1)
  npf = hfl(npolh+1) - 1
  if (npf + nfhol > maxpf) then
   write ( *, 610) 'maxpf',npf+nfhol
   stop
  end if
  read (irdr,*) (pfl(1,i),i=1,npf)
  htsiz = min(prime(nvc+2),maxht)
  call gtime(t0)

  if (hoflag) then
   call dsphfh(aspc2d,atol2d,nvc,nface,nfhol,npolh,maxvc,maxfv, &
        maxiw,maxwk,nvert,npf,vcl,facep,factyp,nrml,fvl,eang,hfl, &
        pfl,htsiz,ht,iwk,wk,ierr)
  else
     call dsphdc(nvc,nface,npolh,vcl,facep,nrml,fvl,eang,hfl,pfl, &
        htsiz,maxiw/4,iwk,ht,ierr)
  end if
  call gtime(initim)
  initim = initim - t0
  if (ierr /= 0) then
   write ( *, 620) ierr
   write (imeas,620) ierr
   stop
  end if
  nrfe = 0
  angmin = 2.0D+00 * d_pi ( )

  do i = 1,nvert
    if (eang(i) > d_pi ( ) + tol) nrfe = nrfe + 1
    if (eang(i) > -1.0D+00 ) angmin = min(angmin,eang(i))
  end do

  angmin = angmin*180.0D+00/ d_pi ( )
  if (hoflag) then
     write (imeas,810) 'fholds',nvc,nface,nvert,npolh,npf,nrfe, &
        angmin,'initim',initim
  else
     write (imeas,810) 'initds',nvc,nface,nvert,npolh,npf,nrfe, &
        angmin,'initim',initim
  end if
!
!     Read in description of interior hole polyhedra. The input format
!     is similar to above with no hole faces and only 1 polyhedron at
!     a time, treated as though the region for each polyh is its hole.
!     IPOLH is index of polyhedron containing hole. IFACH is index of
!     `extreme' face of hole used for connection to outer boundary.
!     One interior hole polyhedron per polyhedral region is assumed.
!     A negative sign for IPOLH indicates hole polyh is hole interface.
!
  holtim = 0.0
  do 20 k = 1,nihol
     read (irdr,*) nvch,nfach,ipolh,ifach
   holint = (ipolh < 0)
   ipolh = abs(ipolh)
     if (nvc+nvch > maxvc .or. nface+nfach >= maxfp .or. &
     npf+nfach > maxpf) then
      if (nvc+nvch > maxvc) write ( *, 610) 'maxvc',nvc+nvch
      if (nface+nfach >= maxfp) write ( *, 610) 'maxfp', &
           nface+nfach+1
      if (npf+nfach > maxpf) write ( *, 610) 'maxpf',npf+nfach
      stop
     end if
     read (irdr,*) ((vcl(j,i),j=1,3),i=nvc+1,nvc+nvch)
     read (irdr,*) (facep(1,i),i=nface+1,nface+nfach+1)
     nverth = facep(1,nface+nfach+1) - 1
     read (irdr,*) (factyp(i),i=nface+1,nface+nfach)
     if (nvert+nverth > maxfv) then
      write ( *, 610) 'maxfv',nvert+nverth
      stop
     end if
   read (irdr,*) (fvl(1,i),i=nvert+1,nvert+nverth)
     read (irdr,*) (pfl(1,i),i=npf+1,npf+nfach)
     htsiz = min(prime(nvch+2),maxht)
     call gtime(t0)
     call dsphih(aspc2d,atol2d,angacc,rdacc,nvc,nface,nvert,npolh, &
        npf,nvch,nfach,ipolh,ifach,holint,maxvc,maxfp,maxfv,maxhf, &
        maxpf,maxiw,maxwk,vcl,facep,factyp,nrml,fvl,eang,hfl,pfl, &
        htsiz,ht,iwk,wk,ierr)

     call gtime(t1)
     holtim = holtim + (t1 - t0)
     if (ierr /= 0) then
      write ( *, 620) ierr
      write (imeas,620) ierr
      stop
     end if
   20 continue
  if (nihol > 0) then
     nrfe = 0
     angmin = 2.0D+00 * d_pi ( )
     do i = 1,nvert
       if (eang(i) > d_pi ( ) + tol) nrfe = nrfe + 1
       if (eang(i) > -1.0D+00 ) angmin = min(angmin,eang(i))
     end do
     angmin = angmin*180.0D+00/ d_pi ( )
     write (imeas,810) 'iholds',nvc,nface,nvert,npolh,npf,nrfe, &
        angmin,'holtim',holtim
  end if
!
!     Decompose polyhedral region into convex parts.
!
  cnt = 1
   40 continue
  call gtime(t0)
  call cvdec3(angacc,rdacc,nvc,nface,nvert,npolh,npf,maxvc,maxfp, &
     maxfv,maxhf,maxpf,maxiw,maxwk,vcl,facep,factyp,nrml,fvl,eang, &
     hfl,pfl,iwk,wk,ierr)
  call gtime(dectim)
  dectim = dectim - t0
  if (ierr /= 0) then
   if (ierr /= 327 .or. cnt >= 3) then
      write ( *, 620) ierr
      write (imeas,620) ierr
      if (ierr /= 327) stop
   end if
  end if
  nrfe = 0
  angmin = 2.0D+00 * d_pi ( )
  do 50 i = 1,nvert
   if (eang(i) > d_pi ( ) + tol) nrfe = nrfe + 1
   if (eang(i) > -1.0D+00 ) angmin = min(angmin,eang(i))
   50 continue
  angmin = angmin*180.0D+00/ d_pi ( )
  write (imeas,820) nvc,nface,nvert,npolh,npf,nrfe,angmin,dectim
  if (ierr == 327 .and. cnt < 3) then
   angacc = angacc - d_pi ( ) /36.0D+00
   if (angacc > 0.0D+00) then
      rdacc = rdacc*0.95D+00
      ierr = 0
        cnt = cnt + 1
      go to 40
   end if
  else if (ierr == 327) then
   stop
  end if
!
!     Decompose faces of polyhedral region into convex subpolygons.
!
  call gtime(t0)
  call cvdecf(aspc2d,atol2d,nvc,nface,nvert,npf,maxvc,maxfp,maxfv, &
     maxpf,maxiw,maxwk,vcl,facep,factyp,nrml,fvl,eang,hfl,pfl,iwk, &
     wk,ierr)
  call gtime(dectim)
  dectim = dectim - t0
  if (ierr /= 0) then
   write ( *, 620) ierr
   write (imeas,620) ierr
   stop
  end if
  write (imeas,830) nvc,nface,nvert,npolh,npf,dectim
!
  if (outmod == 1) then
     write ( *, 630) nvc,nface,nvert,npolh,npf
     write ( *, 640) (i,(vcl(j,i),j=1,3),i=1,nvc)
     write ( *, 650) (i,(facep(j,i),j=1,3),factyp(i), &
        (nrml(j,i),j=1,3),i=1,nface)
     write ( *, 660) (i,(fvl(j,i),j=1,6),eang(i),i=1,nvert)
     write ( *, 670) (hfl(i),i=1,npolh)
     write ( *, 680) (i,(pfl(j,i),j=1,2),i=1,npf)
  else
   nbf = 0
   do f = 1,nface
      if (facep(3,f) == 0) nbf = nbf + 1
   end do
   write ( *, 690) nvc,nbf,nface-nbf
   write ( *, 700) ((vcl(j,i),j=1,3),i=1,nvc)
   do 90 i = 1,2
        write ( *, '(a)' )
      do 80 f = 1,nface
         if (i == 1 .and. facep(3,f) /= 0) go to 80
         if (i == 2 .and. facep(3,f) == 0) go to 80
         k = 0
         j = facep(1,f)
   70          continue
      k = k + 1
      iwk(k) = fvl(1,j)
      j = fvl(3,j)
         if (j /= facep(1,f)) go to 70
         if (k > maxiw) then
            write ( *, 610) 'maxiw',k
            stop
         end if
         write ( *, 710) k,(iwk(j),j=1,k)
   80       continue
   90    continue
  end if
!
!  Initialize data structure for heuristic mesh distribution function
!  and evaluate function at weighted centroid of each polyhedron.
!
  call gtime(t0)

  call dsmdf3(nvc,nface,nvert,npolh,maxiw,maxwk,vcl,facep,nrml,fvl, &
     eang,hfl,pfl,ivrt,xivrt,ifac,xifac,wid,facval,edgval,vrtval, &
     ncface,ncedge,iwk,wk,ierr)

  call gtime(mdftim)

  mdftim = mdftim - t0

  if (ierr /= 0) then
    write ( *, 620) ierr
    write (imeas,620) ierr
    return
  end if

  if (ncface+ncedge+ncedge-2 > maxiw) then
    write ( *, 610) 'maxiw',ncface+ncedge+ncedge-2
    stop
  else if ((ncface+ncedge)*4 > maxwk) then
    write ( *, 610) 'maxwk',(ncface+ncedge)*4
    stop
  end if

  htsiz = min(prime(ncedge),maxht)
  sumf = 0
  sume = 0
  sumv = 0
  do 120 i = 1,npolh
   call prmdf3(i,wid(i),nvc,vcl,nrml,ivrt,xivrt,ifac,xifac,facval, &
        edgval,vrtval,nfcev,nedev,nvrev,factyp,wk,htsiz,maxiw/4,ht, &
        iwk, ierr )
   if (ierr /= 0) then
      write ( *, 620) ierr
      write (imeas,620) ierr
      stop
   end if
   sumf = sumf + nfcev
   sume = sume + nedev
   sumv = sumv + nvrev
   x = 0.0D+00
   y = 0.0D+00
   z = 0.0D+00
   n = 0
   do 110 j = xifac(i),xifac(i+1)-1
      f = abs(ifac(j))
      n = n + (xivrt(f+1) - xivrt(f))
      do k = xivrt(f),xivrt(f+1)-1
         x = x + vcl(1,ivrt(k))
         y = y + vcl(2,ivrt(k))
         z = z + vcl(3,ivrt(k))
      end do
  110    continue
   val(i) = 1.0D+00/mdf3(x/n,y/n,z/n,wid(i),nfcev,nedev,nvrev, &
        factyp,wk,ivrt,facval,edgval,vrtval,vcl)**(1.0D+00/3.0D+00)
  120 continue
!
!     Print arrays from calls to 3 mdf routines.
!
  write (imeas,840) sumf,sume,sumv,mdftim
  write ( *, 720) (i,xifac(i),wid(i),val(i),i=1,npolh)
  write ( *, 730) (i,ifac(i),xivrt(i),facval(i),i=1,nface)
  write ( *, 740) (i,ifac(i),i=nface+1,npf)
  write ( *, 750) (i,ivrt(i),edgval(i),vrtval(i),i=1,nvc)
  write ( *, 760) (i,ivrt(i),edgval(i),i=nvc+1,nvert)
!
  600 format (a60)
  610 format (1x,'*** ',a,' must be increased to',i8)
  620 format (/1x,'ierr=',i5)
  630 format (1x,'nvc,nface,nvert,npolh,npf'/1x,5i7)
  640 format (/' vcl'/(1x,i7,3f15.7))
  650 format (/' facep,factyp,nrml'/(1x,4i7,i5,3f15.7))
  660 format (/' fvl,eang'/(1x,7i7,f15.7))
  670 format (/' hfl'/(1x,10i7))
  680 format (/' pfl'/(1x,3i7))
  690 format (1x,3i7)
  700 format (/(1x,3f15.7))
  710 format (1x,11i7/(8x,10i7))
  720 format (/' xifac,wid,val'/(1x,2i7,2f15.7))
  730 format (/' ifac,xivrt,facval'/(1x,3i7,f15.7))
  740 format ((1x,2i7))
  750 format (/' ivrt,edgval,vrtval'/(1x,2i7,2f15.7))
  760 format ((1x,2i7,f15.7))
  800 format (1x,a60/1x,'input : tol=',d15.7,'   aspc2d=',f9.3, &
     '   atol2d=',f9.3/9x,'angacc=',f9.3,'   rdacc=',f9.5, &
     '   nfhol=',i3,'   nihol=',i3)
  810 format (1x,a6,': nvc=',i7,'   nface=',i7,'   nvert=',i7, &
     '   npolh=',i7/9x,'npf=',i7,'   nrfe=',i7,'   angmin=',f9.3, &
     3x,a6,'=',f9.3)
  820 format (1x,'decomp: nvc=',i7,'   nface=',i7,'   nvert=',i7, &
     '   npolh=',i7/9x,'npf=',i7,'   nrfe=',i7,'   angmin=',f9.3, &
     '   dectim=',f9.3)
  830 format (1x,'decfac: nvc=',i7,'   nface=',i7,'   nvert=',i7, &
     '   npolh=',i7/9x,'npf=',i7,'   dectim=',f9.3)
  840 format (1x,'dsmdf : sumf=',i7,'   sume=',i7,'   sumv=',i7, &
     '   mdftim=',f9.3)

  return
end
subroutine test_meas ( )

!*****************************************************************************80
!
!! TEST_MEAS tests EMNRTH, RADRTH and SANGMN.
!
!  Modified:
!
!    16 August 2005
!
  implicit none

  real ( kind = 8 ), dimension ( 3 ) :: a = (/ 0.0D+00, 0.0D+00, 0.0D+00 /)
  real ( kind = 8 ), dimension ( 3 ) :: b = (/ 1.0D+00, 0.0D+00, 0.0D+00 /)
  real ( kind = 8 ), dimension ( 3 ) :: c = (/ 0.0D+00, 0.0D+00, 0.0D+00 /)
  real ( kind = 8 ), dimension ( 3 ) :: d = (/ 0.0D+00, 0.0D+00, 0.0D+00 /)
  real ( kind = 8 ) emnrth
  real ( kind = 8 ) eta
  integer ( kind = 4 ) irdr
  real ( kind = 8 ) :: l1 = 1.0D+00
  real ( kind = 8 ) :: l2 = 1.0D+00
  real ( kind = 8 ) :: l3 = 1.0D+00
  real ( kind = 8 ) q1
  real ( kind = 8 ) q2
  real ( kind = 8 ) q3
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) r3
  real ( kind = 8 ) radrth
  real ( kind = 8 ) rho
  real ( kind = 8 ) sang(4)
  real ( kind = 8 ) sangmn
  real ( kind = 8 ) sigma
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 5
  real ( kind = 8 ) :: u1 = 0.0D+00
  real ( kind = 8 ) :: u2 = 0.0D+00
  real ( kind = 8 ) :: u3 = 0.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_MEAS'
  write ( *, '(a)' ) '  Compute the three tetrahedron quality measures.'
  write ( *, '(a)' ) '  EMNRTH: eigenvalue measurement;'
  write ( *, '(a)' ) '  RADRTH: inradius / circumradius;'
  write ( *, '(a)' ) '  SANGMN: solid angle measurement'

  do test = 1, test_num

    call random_number ( harvest = c(1:3) )
    call random_number ( harvest = d(1:3) )

    write ( *, '(a)' ) ' '
    write ( *, '(a,3f10.4)' ) '  A = ', a(1:3)
    write ( *, '(a,3f10.4)' ) '  B = ', b(1:3)
    write ( *, '(a,3f10.4)' ) '  C = ', c(1:3)
    write ( *, '(a,3f10.4)' ) '  D = ', d(1:3)

    sigma = sangmn ( a, b, c, d, sang )
    rho = radrth ( a, b, c, d )
    eta = emnrth ( a, b, c, d )

    q1 = rho / eta**3
    r1 = rho / eta**0.75D+00
    q2 = sigma / ( eta * sqrt ( eta ) )
    r2 = sigma / eta**0.75D+00
    q3 = sigma / rho**2
    r3 = sigma / sqrt ( rho )

    l1 = min ( l1, q1 )
    l2 = min ( l2, q2 )
    l3 = min ( l3, q3 )

    u1 = max ( u1, r1 )
    u2 = max ( u2, r2 )
    u3 = max ( u3, r3 )

    write ( *, '(a)' ) ' '
    write ( *, '(a,3e15.7)' ) '  Sigma = ', sigma
    write ( *, '(a,3e15.7)' ) '  Rho =   ', rho
    write ( *, '(a,3e15.7)' ) '  Eta =   ', eta
    write ( *, '(a,3e15.7)' ) '  Q:', q1, q2, q3
    write ( *, '(a,3e15.7)' ) '  U:', r1, r2, r3

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Lower bounds on Q:'
  write ( *, '(2x,3e15.7)' ) l1, l2, l3
  write ( *, '(a)' ) '  Upper bounds on U:'
  write ( *, '(2x,3e15.7)' ) u1, u2, u3

  return
end
subroutine test_prime ( )

!*****************************************************************************80
!
!! TEST_PRIME tests PRIME.
!
!  Modified:
!
!    07 October 2005
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: i_test = 10
  integer ( kind = 4 ) i_uniform
  integer ( kind = 4 ) k
  integer ( kind = 4 ), parameter :: k_max = 15000
  integer ( kind = 4 ) p
  integer ( kind = 4 ) prime
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_PRIME'
  write ( *, '(a)' ) '  PRIME returns a prime number greater than a given'
  write ( *, '(a)' ) '  value I.'
  write ( *, '(a)' ) '  However, its largest stored prime is 14011.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       I  PRIME(I)'
  write ( *, '(a)' ) ' '

  do i = 1, i_test
    k = i_uniform ( 1, k_max, seed )
    p = prime ( k )
    write ( *, '(2x,i6,2x,i6)' ) k, p
  end do

  return
end
subroutine test_ptpolg ( )

!*****************************************************************************80
!
!! TEST_PTPOLG tests PTPOLG.
!
  implicit none

  integer ( kind = 4 ), parameter :: maxn = 100

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  integer ( kind = 4 ) dim
  real ( kind = 8 ) dtol
  integer ( kind = 4 ) i
  integer ( kind = 4 ) inout
  integer ( kind = 4 ) irdr
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  real ( kind = 8 ) nrml(3)
  integer ( kind = 4 ) pgind(0:maxn)
  real ( kind = 8 ) pt(3)
  real ( kind = 8 ) tol
  real ( kind = 8 ) vcl(3,maxn)

  irdr = 5
  tol = 0.000001D+00
  dtol = 10.0D+00 * tol

  read (irdr,*) n,dim,a,b
  if (n < 3) stop

  if (n > maxn) then
    write ( *, 600) 'maxn',n
    stop
  end if

  write ( *, 610) dim,a,b
  write ( *, 620)
  read (irdr,*) (vcl(1,i),vcl(2,i),i=1,n)
  pgind(0) = n

  do i = 1,n
    pgind(i) = i
    if (dim == 3) vcl(3,i) = a*vcl(1,i) + b*vcl(2,i)
    write ( *, 620) i,(vcl(j,i),j=1,dim)
  end do

  if (dim == 3) then
    c = sqrt(1.0D+00 + a**2 + b**2)
    nrml(1) = -a/c
    nrml(2) = -b/c
    nrml(3) = 1.0D+00/c
  end if

  write ( *, 620)

   20 continue
   read (irdr,*,end=30) pt(1),pt(2)
   if (dim == 3) pt(3) = a*pt(1) + b*pt(2)
   call ptpolg(dim,3,n,1,pgind,vcl,pt,nrml,dtol,inout)
     write ( *, 630) inout,(pt(i),i=1,dim)
  go to 20

   30 continue

  return

  600 format (1x,'*** ',a,' must be increased to',i8)
  610 format (1x,'dim=',i2,3x,'a,b=',2f15.7)
  620 format (1x,i5,3f15.7)
  630 format (1x,'inout=',i3,3x,'pt=',3f15.7)
end
subroutine test_rotiar ( )

!*****************************************************************************80
!
!! TEST_ROT tests ROTIAR.
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: shift = 3

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_ROTIAR'
  write ( *, '(a)' ) '  ROTIAR shifts an integer array.'

  do i = 1, n
    a(i) = i
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Initial vector:'
  write ( *, '(a)' ) ' '

  write ( *, '(10i5)' ) a(1:n)

  call rotiar ( n, a, shift )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Vector after applying a shift of ', shift
  write ( *, '(a)' ) ' '

  write ( *, '(10i5)' ) a(1:n)

  return
end
subroutine test_shap ( )

!*****************************************************************************80
!
!! TEST_SHAP tests EMNRTH, RADRTH, and SANGMN.
!
  implicit none

  real ( kind = 8 ) d_pi
  real ( kind = 8 ) tol

  integer ( kind = 4 ) in,irdr
  real ( kind = 8 ) a(3),b(3),c(3),centre(3),d(3),dang(6),sang(4)
  real ( kind = 8 ) e,emnrth,eta,rad,radrth,radsq,rho1,rho2,samin1
  real ( kind = 8 ) samin2,sangmn,xc,xc0,xc1,xcinc,xd,xd0,xd1,xdinc
  real ( kind = 8 ) yc,yc0,yc1,ycinc,yd,yd0,yd1,ydinc,zd,zd0,zd1
  real ( kind = 8 ) zdinc

  data a/0.0D+00,0.0D+00,0.0D+00/
  data b/1.0D+00,0.0D+00,0.0D+00/
  data c(3)/0.0D+00/

  irdr = 5
  tol = 0.000001D+00

  read (irdr,*) xc0,xc1,xcinc,yc0,yc1,ycinc,xd0,xd1,xdinc, &
     yd0,yd1,ydinc,zd0,zd1,zdinc
  write ( *, 600)

  do 50 xc = xc0,xc1,xcinc

    c(1) = xc

    do 40 yc = yc0,yc1,ycinc

      c(2) = yc

      do 30 xd = xd0,xd1,xdinc

        d(1) = xd

        do yd = yd0,yd1,ydinc

          d(2) = yd

          do zd = zd0,zd1,zdinc

            d(3) = zd
            call sdang(a,b,c,d,sang,dang)
            samin1 = min(sang(1),sang(2),sang(3),sang(4)) - d_pi ( )
            samin2 = 2.0D+00 * asin(sangmn(a,b,c,d,sang))
            call insph(a,b,c,d,centre,rad)
            call ccsph(.false.,a,b,c,d,e,centre,radsq,in)
            rho1 = 3.0D+00 * rad/sqrt(abs(radsq))
            rho2 = radrth(a,b,c,d)
            eta = emnrth(a,b,c,d)
            write ( *, 610) c(1),c(2),d(1),d(2),d(3), samin2,rho2,eta

            if (abs(samin1 - samin2) > 100.0D+00*tol) then
              write ( *, 620) 'samin:',samin1,samin2
            end if

            if (abs(rho1 - rho2) > 100.0D+00*tol) then
                    write ( *, 620) 'rho:',rho1,rho2
            end if

          end do

        end do

   30       continue
   40    continue
   50 continue
  600 format (3x,'xc',6x,'yc',6x,'xd',6x,'yd',6x,'zd',7x,'samin',8x, &
     'rho',9x,'eta')
  610 format (5f8.3,3f12.6)
  620 format (10x,a,2f20.14)

  return
end
subroutine test_shr ( )

!*****************************************************************************80
!
!! TEST_SHR tests DIAM2, SHRNK2 and WIDTH2.
!
  implicit none

  integer ( kind = 4 ) ierr

  integer ( kind = 4 ), parameter :: maxn = 100

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) irdr
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nshr
  integer ( kind = 4 ) iedge(0:maxn)
  real ( kind = 8 ) diamsq,tol,widsq
  real ( kind = 8 ) sdist(0:maxn-1)
  real ( kind = 8 ) xc(0:maxn),xs(0:maxn),yc(0:maxn),ys(0:maxn)

  irdr = 5
  read (irdr,*) n, tol

  if (n < 3) then
    return
  end if

  if (n > maxn) then
    write ( *, 600) 'maxn',n
    return
  end if

  read (irdr,*) (xc(i),yc(i),sdist(i),i=0,n-1)
  xc(n) = xc(0)
  yc(n) = yc(0)
  write ( *, 620) n,(xc(i),yc(i),i=0,n)
  xs(0) = 0.0D+00
  ys(0) = 0.0D+00

  call shrnk2 ( n, xc, yc, sdist, nshr, xs, ys, iedge, ierr )

  if (ierr /= 0) then
    write ( *, 610) ierr
    return
  end if

  write ( *, 620) nshr,(xs(i),ys(i),i=0,nshr)

  call diam2(n,xc(1),yc(1),i1,i2,diamsq)

  if (ierr /= 0) then
    write ( *, 610) ierr
    return
  end if

  write ( *, 630) i1,i2,diamsq
  call width2(n,xc(1),yc(1),i1,i2,widsq, ierr )

  if (ierr /= 0) then
    write ( *, 610) ierr
    return
  end if

  write ( *, 630) i1,i2,widsq
  if (nshr < 3) stop
  call diam2(nshr,xs(1),ys(1),i1,i2,diamsq)

  if (ierr /= 0) then
    write ( *, 610) ierr
    return
  end if

  write ( *, 630) i1,i2,diamsq
  call width2(nshr,xs(1),ys(1),i1,i2,widsq, ierr )
  if (ierr /= 0) then
   write ( *, 610) ierr
   return
  end if
  write ( *, 630) i1,i2,widsq

  600 format (1x,'*** ',a,' must be increased to',i8)
  610 format (/1x,'ierr=',i5)
  620 format (/1x,i5/(1x,2f15.7))
  630 format (/1x,2i5,f18.7)

  return
end
subroutine test_shr3 ( )

!*****************************************************************************80
!
!! TEST_SHR3
!
!     Driver program for testing routines DSCPH, DIAM3, WIDTH3, VOLCPH,
!        SHRNK3, XPGHPL, RMCPFC, RMCLED.
!
  implicit none

  integer ( kind = 4 ) ierr
  real ( kind = 8 ) tol

  integer ( kind = 4 ) maxed,maxfv,maxht,maxhv,maxiw,maxvc,maxwk
  parameter (maxed = 100)
  parameter (maxfv = 100)
  parameter (maxht = 53)
  parameter (maxhv = 25)
  parameter (maxiw = 100)
  parameter (maxvc = 100)
  parameter (maxwk = 100)

  integer ( kind = 4 ) edge(4,maxed),fvl(5,maxfv),ht(0:maxed-1),hvl(maxhv+1)
  integer ( kind = 4 ) iwk(maxiw),sfvl(5,maxfv),shvl(maxhv)
  integer ( kind = 4 ) htsiz,i,i1,i2,irdr,j,k,nface,nfold,nsface,nsvc
  integer ( kind = 4 ) nsvert,nvc,nvert,nvold,prime
  real ( kind = 8 ) eang(maxfv),nrml(3,maxhv),svcl(3,maxvc)
  real ( kind = 8 ) vcl(3,maxvc),wk(maxwk)
  real ( kind = 8 ) diamsq,sdist,vol,volcph,wid

  irdr = 5
  read (irdr,*) nvc,nface,sdist,tol

  if (nvc > maxvc) then
   write ( *, 610) 'maxvc',nvc
   stop
  else if (nface > maxhv) then
   write ( *, 610) 'maxhv',nface
   stop
  end if
  read (irdr,*) ((vcl(j,i),j=1,3),i=1,nvc)
  read (irdr,*) (hvl(i),i=1,nface+1)
  nvert = hvl(nface+1) - 1

  if (nvert > maxfv) then
   write ( *, 610) 'maxfv',nvert
   stop
  else if (nvert > maxiw) then
   write ( *, 610) 'maxiw',nvert
   stop
  end if

  read (irdr,*) (fvl(1,i),i=1,nvert)
  htsiz = min(prime(nvc+2),maxht)
  call dscph(nvc,nface,vcl,hvl,nrml,fvl,eang,htsiz,maxed,edge,ht,ierr)

  if (ierr /= 0) then
    write ( *, 620) ierr
    return
  end if

  write ( *, 630) nvc,nface,nvert,tol
  write ( *, 600) 'vcl'
  write ( *, 640) (i,(vcl(j,i),j=1,3),i=1,nvc)
  write ( *, 600) 'hvl,nrml'
  write ( *, 650) (i,hvl(i),(nrml(j,i),j=1,3),i=1,nface)
  write ( *, 600) 'fvl,eang'
  write ( *, 660) (i,(fvl(j,i),j=1,5),eang(i),i=1,nvert)
  nfold = nface
  nvold = nvert
  call rmcpfc(nface,nvert,hvl,nrml,fvl,eang,iwk)
  call rmcled(nface,nvert,hvl,fvl)
  if (nfold > nface .or. nvold > nvert) then
     write ( *, 720) nface,nvert
     write ( *, 600) 'updated hvl,nrml'
     write ( *, 650) (i,hvl(i),(nrml(j,i),j=1,3),i=1,nface)
     write ( *, 600) 'updated fvl,eang'
     write ( *, 660) (i,(fvl(j,i),j=1,5),eang(i),i=1,nvert)
  end if
!
  vol = volcph(nface,vcl,hvl,fvl)
  call diam3(nvc,vcl,i1,i2,diamsq)
  write ( *, 670) vol,sqrt(diamsq),i1,i2
  call width3(nface,vcl,hvl,nrml,fvl,maxiw,i1,i2,wid,iwk, ierr )
  if (ierr /= 0) then
   write ( *, 620) ierr
   stop
  end if
  write ( *, 680) wid,i1,i2
  call shrnk3(sdist,nface,vcl,hvl,nrml,fvl,eang,min(maxfv,maxvc), &
     maxiw,maxwk,nsvc,nsface,nsvert,svcl,shvl,sfvl,iwk,wk, ierr )
  if (ierr /= 0) then
   write ( *, 620) ierr
   stop
  end if
  write ( *, 600) 'shrunken polyhedron'
  write ( *, 690) nsvc,nsface,nsvert,sdist
  if (nsface < 4) stop
  write ( *, 600) 'svcl'
  write ( *, 640) (i,(svcl(j,i),j=1,3),i=1,nsvc)
  write ( *, 600) 'shvl'
  write ( *, 700) (shvl(i),i=1,nface)
  write ( *, 600) 'sfvl'
  write ( *, 710) (i,(sfvl(j,i),j=1,5),i=1,nsvert)
  vol = volcph(nsface,svcl,shvl,sfvl)
  call diam3(nsvc,svcl,i1,i2,diamsq)
  write ( *, 670) vol,sqrt(diamsq),i1,i2
  if (nsface < nface) then
   i = shvl(nsface+1)
   j = nsface + 2
   do 10 k = shvl(nsface+1)+1,nface
      if (j > nface .or. k /= shvl(j)) then
         nrml(1,i) = nrml(1,k)
         nrml(2,i) = nrml(2,k)
         nrml(3,i) = nrml(3,k)
         i = i + 1
      else
         j = j + 1
      end if
   10    continue
  end if
  call width3(nsface,svcl,shvl,nrml,sfvl,maxiw,i1,i2,wid,iwk, ierr )
  if (ierr /= 0) then
   write ( *, 620) ierr
   stop
  end if
  write ( *, 680) wid,i1,i2

  600 format (/1x,a)
  610 format (1x,'*** ',a,' must be increased to',i8)
  620 format (/1x,'ierr=',i5)
  630 format (1x,'nvc=',i5,'   nface=',i5,'   nvert=',i5,'   tol=', &
     d15.7)
  640 format (1x,i5,3f15.7)
  650 format (1x,2i5,3f15.7)
  660 format (1x,6i5,f15.7)
  670 format (/1x,'vol=',f15.7,'   diam=',f15.7,'   i1,i2=',2i5)
  680 format (1x,'wid=',f15.7,'   i1,i2=',2i5)
  690 format (1x,'nsvc=',i5,'   nsface=',i5,'   nsvert=',i5,'   sdist=', &
     f15.7)
  700 format (1x,15i5)
  710 format (1x,6i5)
  720 format (/1x,'nface=',i5,'   nvert=',i5)

  return
end
subroutine test_smpx ( )

!*****************************************************************************80
!
!! TEST_SMPX tests BARYCK, CCSPHK and OPSIDK.
!
  implicit none
  integer ( kind = 4 ) kmax
  parameter (kmax = 10)

  integer ( kind = 4 ) case,i,in,ind(kmax+1),ipvt(kmax-1),irdr,k,op,opsidk
  real ( kind = 8 ) alpha(kmax+1),centre(kmax),mat(kmax*kmax)
  real ( kind = 8 ) pta(kmax),ptb(kmax),radsq,vcl(kmax*(kmax+1))
  logical degen

  irdr = 5

  write ( *, '(a)' ) 'enter k:'
  read (irdr,*) k
  if (k <= 1 .or. k > kmax) stop

  do i = 1,k
    ind(i) = k + 1 - i
  end do

  ind(k+1) = k + 1

   20 continue
   write ( *, '(a)' ) 'enter case:'
   read (irdr,*) case
   if (case < 1 .or. case > 4) stop
   if (case == 1) then
       write ( *, '(a)' ) 'ccsphk - enter v(1:k+1),pt:'
       read (irdr,*) (vcl(i),i=1,k*k+k),(pta(i),i=1,k)
       call ccsphk(k,.true.,ind,vcl,pta,centre,radsq,in,mat,ipvt)
       write ( *, 600) 'centre=',(centre(i),i=1,k)
      if (radsq >= 0.0D+00) radsq = sqrt(radsq)
       write ( *, '(a)' ) 'rad=',radsq,'   in=',in
   else if (case == 2) then
       write ( *, '(a)' ) 'baryck - enter v(1:k+1),pt:'
       read (irdr,*) (vcl(i),i=1,k*k+k),(pta(i),i=1,k)
       call baryck(k,ind,vcl,pta,alpha,degen,mat,ipvt)
       write ( *, '(a)' ) 'degen=',degen
       if (.not. degen) write ( *, 600) 'alpha=', &
           (alpha(ind(i)),i=1,k+1)
   else if (case == 3) then
       write ( *, '(a)' ) 'opsidk - enter v(1:k),pta,ptb:'
       read (irdr,*) (vcl(i),i=1,k*k),(pta(i),i=1,k),(ptb(i),i=1,k)
       op = opsidk(k,ind,vcl,.false.,pta,ptb,mat,alpha)
       write ( *, '(a)' ) 'opsidk=',op
   else
       write ( *, '(a)' ) 'opsidk - enter v(1:k),pta:'
       read (irdr,*) (vcl(i),i=1,k*k),(pta(i),i=1,k)
       op = opsidk(k,ind,vcl,.true.,pta,pta,mat,alpha)
       write ( *, '(a)' ) 'opsidk=',op
   end if

  go to 20
  600 format (a,6f12.6/3x,6f12.6)
end
subroutine test_ihpsrt ( )

!*****************************************************************************80
!
!! TEST_IHPSRT tests IHPSRT.
!
  implicit none

  integer ( kind = 4 ), parameter :: maxk = 4
  integer ( kind = 4 ), parameter :: maxn = 100

  integer ( kind = 4 ) axis
  real ( kind = 8 ) da(maxk,maxn)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(maxk,maxn)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) map(maxn)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nptav
  real ( kind = 8 ), dimension(maxk) :: scale = (/ &
    1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00 /)
  integer ( kind = 4 ) seed
  real ( kind = 8 ), dimension(maxk) :: trans = (/ &
    0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_IHPSRT'
  write ( *, '(a)' ) '  IHPSRT sorts multidimensional integer vectors.'

  k = 2
  n = 20
  seed = 1952
  axis = 4
  nptav = 5

  if ( k < 1 .or. maxk < k .or. n < 1 .or. maxn < n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Bad data.'
    return
  end if

  do j = 1, n
    map(j) = j
  end do

  call randpt ( k, n, seed, axis, nptav, scale, trans, maxk, da )

  ia(1:k,1:n) = int ( n * da(1:k,1:n) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Unsorted array produced by RANDPT:'
  write ( *, '(a)' ) ' '

  write ( *, '(2x,i6,2x,4i6)' ) j, ia(1:k,j)

  call ihpsrt(k,n,maxk,ia,map)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  After sorting by IHPSRT:'
  write ( *, '(a)' ) ' '

  do j = 1, n
    write ( *, '(2x,i6,2x,4i6)' ) map(j), (ia(i,map(j)),i=1,k)
  end do

  return
end
subroutine test_dhpsrt ( )

!*****************************************************************************80
!
!! TEST_DHPSRT tests DHPSRT.
!
  implicit none

  integer ( kind = 4 ), parameter :: maxk = 4
  integer ( kind = 4 ), parameter :: maxn = 100

  integer ( kind = 4 ) axis
  real ( kind = 8 ) da(maxk,maxn)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) map(maxn)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nptav
  real ( kind = 8 ), dimension(maxk) :: scale = (/ &
    1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00 /)
  integer ( kind = 4 ) seed
  real ( kind = 8 ), dimension(maxk) :: trans = (/ &
    0.0D+00,0.0D+00,0.0D+00,0.0D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_DHPSRT'
  write ( *, '(a)' ) '  DHPSRT sorts multidimensional double precision vectors.'

  k = 2
  n = 20
  seed = 1952
  axis = 4
  nptav = 5

  if ( k < 1 .or. maxk < k .or. n < 1 .or. maxn < n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DRDSORT - Fatal error!'
    write ( *, '(a)' ) '  Bad data.'
    return
  end if

  do j = 1, n
    map(j) = j
  end do

  call randpt ( k, n, seed, axis, nptav, scale, trans, maxk, da )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Unsorted array produced by RANDPT:'
  write ( *, '(a)' ) ' '

  do j = 1, n
    write ( *, '(2x,i5,2x,4f15.7)' ) j, da(1:k,j)
  end do

  call dhpsrt ( k, n, maxk, da, map )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  After sorting by DHPSRT:'
  write ( *, '(a)' ) ' '

  do j = 1, n
    write ( *, '(2x,i5,2x,4f15.7)' ) map(j), da(1:k,map(j))
  end do

  return
end
subroutine test_tet ( )

!*****************************************************************************80
!
!! TEST_TET
!
!     Driver program for testing routines CCSPH, BARYTH, OPSIDE, VOLTH,
!        ANGLE3, SDANG, INSPH, RADRTH, EMNRTH, SANGMN, CCRADI.
!
  implicit none

  real ( kind = 8 ) d_pi
  real ( kind = 8 ) tol

  integer ( kind = 4 ) case,i,in,irdr,op,opside
  real ( kind = 8 ) a(3),alpha(4),b(3),c(3),centre(3),d(3),dang(6)
  real ( kind = 8 ) e(3),sang(4),rad,radsq,ratio,vol,volth
  real ( kind = 8 ) ccradi,emnrth,radrth,sangmn
  logical degen

  irdr = 5

10 continue

   write ( *, '(a)' ) 'enter case:'
   read (irdr,*) case
   if (case < 1 .or. case > 7) stop

   if (case == 1) then
       write ( *, '(a)' ) 'ccsph - enter a,b,c,d,e:'
       read (irdr,*) a(1),a(2),a(3),b(1),b(2),b(3),c(1),c(2),c(3), &
           d(1),d(2),d(3),e(1),e(2),e(3)
       call ccsph(.true.,a,b,c,d,e,centre,radsq,in)
       write ( *, 600) 'centre=',(centre(i),i=1,3)
      if (radsq >= 0.0D+00) radsq = sqrt(radsq)
       write ( *, '(a)') 'rad=',radsq,'   in=',in
   else if (case == 2) then
       write ( *, '(a)') 'baryth - enter a,b,c,d,e:'
       read (irdr,*) a(1),a(2),a(3),b(1),b(2),b(3),c(1),c(2),c(3), &
           d(1),d(2),d(3),e(1),e(2),e(3)
       call baryth(a,b,c,d,e,alpha,degen)
       write ( *, 600) 'alpha=',(alpha(i),i=1,4)
       write ( *, '(a)') 'degen=',degen
   else if (case == 3) then
       write ( *, '(a)') 'opside - enter a,b,c,d,e:'
       read (irdr,*) a(1),a(2),a(3),b(1),b(2),b(3),c(1),c(2),c(3), &
           d(1),d(2),d(3),e(1),e(2),e(3)
       op = opside(a,b,c,d,e)
       write ( *, '(a)') 'opside=',op
   else if (case == 4) then
       write ( *, '(a)') 'volth - enter a,b,c,d:'
       read (irdr,*) a(1),a(2),a(3),b(1),b(2),b(3),c(1),c(2),c(3), &
           d(1),d(2),d(3)
       vol = volth(a,b,c,d)/6.0D+00
       write ( *, '(a)') 'volth=',vol
   else if (case == 5) then
       write ( *, '(a)') 'sdang,sangmn - enter a,b,c,d:'
       read (irdr,*) a(1),a(2),a(3),b(1),b(2),b(3),c(1),c(2),c(3), &
           d(1),d(2),d(3)
       call sdang(a,b,c,d,sang,dang)
       write ( *, 600) 'sang=',(sang(i) - d_pi ( ),i=1,4)
       write ( *, 600) 'dang=',(dang(i),i=1,6)
       write ( *, '(a)') 'sin(sangmn/2)=',sangmn(a,b,c,d,sang)
       write ( *, 600) 'sang=',(2.0D+00*asin(sang(i)),i=1,4)
   else if (case == 6) then
       write ( *, '(a)') 'insph,ccradi,radrth - enter a,b,c,d:'
       read (irdr,*) a(1),a(2),a(3),b(1),b(2),b(3),c(1),c(2),c(3), &
           d(1),d(2),d(3)
       call insph(a,b,c,d,centre,rad)
       write ( *, 600) 'incent=',(centre(i),i=1,3)
       radsq = ccradi(a,b,c,d)
      if (radsq == 0.0D+00) then
         ratio = 0.0D+00
      else
         ratio = rad*sqrt(radsq)*12.0D+00
      end if
       write ( *, '(a)') 'inrad=',rad,'   ratio=',ratio
       write ( *, '(a)') 'radrth=',radrth(a,b,c,d)
   else
       write ( *, '(a)') 'emnrth - enter a,b,c,d:'
       read (irdr,*) a(1),a(2),a(3),b(1),b(2),b(3),c(1),c(2),c(3), &
           d(1),d(2),d(3)
      write ( *, '(a)') 'emnrth=',emnrth(a,b,c,d)
   end if

  go to 10
  600 format (a,6f12.6)
end
subroutine test_tri3 ( )

!*****************************************************************************80
!
!! TEST_TRI3 tests EQDIS3, TRIPR3.
!
  implicit none

  real ( kind = 8 ) d_pi
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) msglvl
  real ( kind = 8 ) tol

  integer ( kind = 4 ) maxbt,maxfc,maxfp,maxfv,maxhf,maxht,maxiw,maxnt,maxpf
  integer ( kind = 4 ) maxte,maxvc,maxvm,maxwk
  parameter (maxbt = 10000)
  parameter (maxfc = 45000)
  parameter (maxfp = 1000)
  parameter (maxfv = 5000)
  parameter (maxhf = 200)
  parameter (maxht = 13411)
  parameter (maxiw = 12000)
  parameter (maxnt = 18000)
  parameter (maxpf = 2000)
  parameter (maxte = 8000)
  parameter (maxvc = 5000)
  parameter (maxvm = 8000)
  parameter (maxwk = 19000)
!
!     Parameters for relative frequency tables.
!
  integer ( kind = 4 ) nf
  real ( kind = 8 ) am,av,hm,hv
  parameter (nf = 10)
  parameter (am = 0.0D+00, hm = 0.1D+00)
  parameter (av = -0.1D+00, hv = 0.2D+00)
!
  integer ( kind = 4 ) btl(3,maxbt),btst(maxfp+1),edno(maxfv),edst(maxfv)
  integer ( kind = 4 ) facep(3,maxfp),factyp(maxfp),fc(7,maxfc),fcst(maxfp+1)
  integer ( kind = 4 ) fvl(6,maxfv),hfl(maxhf),ht(maxht),iwk(maxiw),ntetra(maxhf)
  integer ( kind = 4 ) ntrif(maxhf),pfl(2,maxpf),tetra(4,maxte),vm(maxvm)
  integer ( kind = 4 ) xfhv(3,maxhf+1)
  integer ( kind = 4 ) a,b,c,cnt,crit,d,e,f,htsiz,i,ifach,imeas,in,ipolh,irdr,j,k
  integer ( kind = 4 ) kfc,kt,kvm,nbt,nface,nfach,nfhol,nfph,nihol,nit,nlo,nmin
  integer ( kind = 4 ) npf,npolh,nrfe,nt,nteta,ntetd,ntri,nvc,nvch,nvert,nverth
  integer ( kind = 4 ) outmod,p,prime,utet
  real    dectim,eqdtim,holtim,initim,tritim,t0,t1
  real ( kind = 8 ) eang(maxfv),h(maxhf),nrml(3,maxfp),psi(maxhf)
  real ( kind = 8 ) eta(maxnt),rho(maxnt),rvol(maxnt),sig(maxnt)
  real ( kind = 8 ) sa(4),vcl(3,maxvc),vol(maxhf),wk(maxwk)
  real ( kind = 8 ) angacc,angedg,angmin,aspc2d,atol2d,dmin,emnrth
  real ( kind = 8 ) efreq(0:nf),emax,emean,emin,estdv,h3,kappa,radsq
  real ( kind = 8 ) rdacc,rfreq(0:nf),rmax,rmean,rmin,rstdv,sf
  real ( kind = 8 ) sfreq(0:nf),shrf,smax,smean,smin,sstdv,umdf3
  real ( kind = 8 ) vfreq(0:nf),vmax,vmean,vmin,volth,vstdv
  real ( kind = 8 ) radrth,sangmn
  logical hflag,hoflag,holint,nsflag
  character rgname*60
  external umdf3
!
!     OUTMOD = 0 : Output nothing (except for measurements).
!     OUTMOD = 1 : Output polyh decomp, triangulation data structures.
!     OUTMOD = 2 : Output VCL and 3 vertex indices for each triangle,
!        boundary triangles (on faces of decomp) before interior ones.
!     OUTMOD = 3 : Output VCL and 3 vertex indices for each triangle
!        on faces with FACTYP(F) > 0, boundary faces before interior.
!     OUTMOD = -1, -2, -3: Same as + case, but also output tetrahedron
!        list to unit UTET.
!     NFHOL <= 0: Simpler routine DSPHDC is called instead of DSPHFH.
!     0 <= KAPPA <= 1: Set HFLAG to TRUE; NMIN < 0: Set NSFLAG to TRUE.
!
  irdr = 5
  imeas = 7
  utet = 1
  read (irdr,600) rgname
  read (irdr,*) tol, aspc2d,atol2d,angacc,rdacc,angedg,kappa,dmin, &
     nmin,ntetd,crit,shrf

  read (irdr,*) nvc,nface,nfhol,npolh,nihol,outmod,msglvl
  write (imeas,800) rgname,tol, aspc2d,atol2d,angacc,rdacc,nfhol, &
     nihol,angedg,kappa,dmin,nmin,ntetd,crit,shrf
  hoflag = (nfhol > 0)
  nfhol = max(nfhol,0)
  hflag = (kappa >= 0.0D+00 .and. kappa <= 1.0D+00)
  nsflag = (nmin < 0)
  nmin = abs(nmin)
  aspc2d = aspc2d* d_pi ( ) /180.0D+00
  atol2d = atol2d* d_pi ( ) /180.0D+00
  angacc = angacc* d_pi ( ) /180.0D+00
  angedg = angedg* d_pi ( ) /180.0D+00
  nfph = nface + nfhol
  if (nvc > maxvc .or. nfph >= maxfp .or. npolh >= maxhf) &
  then
   if (nvc > maxvc) write ( *, 610) 'maxvc',nvc
   if (nfph >= maxfp) write ( *, 610) 'maxfp',nfph+1
   if (npolh >= maxhf) write ( *, 610) 'maxhf',npolh+1
   stop
  end if
  read (irdr,*) ((vcl(j,i),j=1,3),i=1,nvc)
!
!  The outer polygons of NFACE faces should appear first, followed by
!  the inner polygons of NFHOL holes (in order of index of faces
!  containing holes). Face types for NFACE faces should be positive.
!  Face types for NFHOL holes are +F or -F where F is index of face
!  containing hole; positive (negative) sign indicates hole polygon
!  is oriented CCW (CW) in polyhedron when viewed from outside.
!  Holes must only be on boundary faces of polyhedral region.
!
  read (irdr,*) (facep(1,i),i=1,nfph+1)
  nvert = facep(1,nfph+1) - 1
  read (irdr,*) (factyp(i),i=1,nfph)
  if (nvert > maxfv) then
   write ( *, 610) 'maxfv',nvert
   stop
  end if
!
!     Head vertex of each (outer) face must be a strictly convex vertex.
!     Head vertex of each hole may be arbitrary.
!     Positive (negative) sign in PFL(1,*) indicates face is oriented
!     CCW (CW) in polyhedron when viewed from outside.
!     Hole polygons should not be included in PFL.
!
  read (irdr,*) (fvl(1,i),i=1,nvert)
  read (irdr,*) (hfl(i),i=1,npolh+1)
  npf = hfl(npolh+1) - 1
  if (npf + nfhol > maxpf) then
   write ( *, 610) 'maxpf',npf+nfhol
   stop
  end if
  read (irdr,*) (pfl(1,i),i=1,npf)
  htsiz = min(prime(nvc+2),maxht)
  call gtime(t0)

  if (hoflag) then
   call dsphfh(aspc2d,atol2d,nvc,nface,nfhol,npolh,maxvc,maxfv, &
        maxiw,maxwk,nvert,npf,vcl,facep,factyp,nrml,fvl,eang,hfl, &
        pfl,htsiz,ht,iwk,wk,ierr)
  else
     call dsphdc(nvc,nface,npolh,vcl,facep,nrml,fvl,eang,hfl,pfl, &
        htsiz,maxiw/4,iwk,ht,ierr)
  end if

  call gtime(initim)
  initim = initim - t0

  if (ierr /= 0) then
    write ( *, 620) ierr
    write (imeas,620) ierr
    stop
  end if

  nrfe = 0
  angmin = 2.0D+00* d_pi ( )

  do i = 1,nvert
    if (eang(i) > d_pi ( ) + tol) nrfe = nrfe + 1
    if (eang(i) > -1.0D+00) angmin = min(angmin,eang(i))
  end do

  angmin = angmin*180.0D+00/ d_pi ( )
  if (hoflag) then
     write (imeas,810) 'fholds',nvc,nface,nvert,npolh,npf,nrfe, &
        angmin,'initim',initim
  else
     write (imeas,810) 'initds',nvc,nface,nvert,npolh,npf,nrfe, &
        angmin,'initim',initim
  end if
!
!     Read in description of interior hole polyhedra. The input format
!     is similar to above with no hole faces and only 1 polyhedron at
!     a time, treated as though the region for each polyh is its hole.
!     IPOLH is index of polyhedron containing hole. IFACH is index of
!     `extreme' face of hole used for connection to outer boundary.
!     One interior hole polyhedron per polyhedral region is assumed.
!     A negative sign for IPOLH indicates hole polyh is hole interface.
!
  holtim = 0.0
  do 20 k = 1,nihol
     read (irdr,*) nvch,nfach,ipolh,ifach
   holint = (ipolh < 0)
   ipolh = abs(ipolh)
     if (nvc+nvch > maxvc .or. nface+nfach >= maxfp .or. &
     npf+nfach > maxpf) then
      if (nvc+nvch > maxvc) write ( *, 610) 'maxvc',nvc+nvch
      if (nface+nfach >= maxfp) write ( *, 610) 'maxfp', &
           nface+nfach+1
      if (npf+nfach > maxpf) write ( *, 610) 'maxpf',npf+nfach
      stop
     end if
     read (irdr,*) ((vcl(j,i),j=1,3),i=nvc+1,nvc+nvch)
     read (irdr,*) (facep(1,i),i=nface+1,nface+nfach+1)
     nverth = facep(1,nface+nfach+1) - 1
     read (irdr,*) (factyp(i),i=nface+1,nface+nfach)

     if (nvert+nverth > maxfv) then
       write ( *, 610) 'maxfv',nvert+nverth
       stop
     end if

   read (irdr,*) (fvl(1,i),i=nvert+1,nvert+nverth)
     read (irdr,*) (pfl(1,i),i=npf+1,npf+nfach)
     htsiz = min(prime(nvch+2),maxht)
     call gtime(t0)
     call dsphih(aspc2d,atol2d,angacc,rdacc,nvc,nface,nvert,npolh, &
        npf,nvch,nfach,ipolh,ifach,holint,maxvc,maxfp,maxfv,maxhf, &
        maxpf,maxiw,maxwk,vcl,facep,factyp,nrml,fvl,eang,hfl,pfl, &
        htsiz,ht,iwk,wk,ierr)

     call gtime(t1)
     holtim = holtim + (t1 - t0)
     if (ierr /= 0) then
      write ( *, 620) ierr
      write (imeas,620) ierr
      stop
     end if
   20 continue
  if (nihol > 0) then
     nrfe = 0
     angmin = 2.0D+00* d_pi ( )
     do 30 i = 1,nvert
      if (eang(i) > d_pi ( ) + tol) nrfe = nrfe + 1
      if (eang(i) > -1.0D+00) angmin = min(angmin,eang(i))
   30    continue
     angmin = angmin*180.0D+00/ d_pi ( )
     write (imeas,810) 'iholds',nvc,nface,nvert,npolh,npf,nrfe, &
        angmin,'holtim',holtim
  end if
!
!     Decompose polyhedral region into convex parts.
!
  cnt = 1
   40 continue
  call gtime(t0)
  call cvdec3(angacc,rdacc,nvc,nface,nvert,npolh,npf,maxvc,maxfp, &
     maxfv,maxhf,maxpf,maxiw,maxwk,vcl,facep,factyp,nrml,fvl,eang, &
     hfl,pfl,iwk,wk,ierr)
  call gtime(dectim)
  dectim = dectim - t0
  if (ierr /= 0) then
   if (ierr /= 327 .or. cnt >= 3) then
      write ( *, 620) ierr
      write (imeas,620) ierr
      if (ierr /= 327) stop
   end if
  end if
  nrfe = 0
  angmin = 2.0D+00* d_pi ( )
  do i = 1,nvert
   if (eang(i) > d_pi ( ) + tol) nrfe = nrfe + 1
   if (eang(i) > -1.0D+00) angmin = min(angmin,eang(i))
  end do
  angmin = angmin*180.0D+00/ d_pi ( )
  write (imeas,820) nvc,nface,nvert,npolh,npf,nrfe,angmin,dectim
  if (ierr == 327 .and. cnt < 3) then
   angacc = angacc - d_pi ( ) /36.0D+00
   if (angacc > 0.0D+00) then
      rdacc = rdacc*0.95D+00
      ierr = 0
        cnt = cnt + 1
      go to 40
   end if
  else if (ierr == 327) then
   stop
  end if
!
!     Decompose faces of polyhedral region into convex subpolygons.
!
  call gtime(t0)
  call cvdecf(aspc2d,atol2d,nvc,nface,nvert,npf,maxvc,maxfp,maxfv, &
     maxpf,maxiw,maxwk,vcl,facep,factyp,nrml,fvl,eang,hfl,pfl,iwk, &
     wk,ierr)
  call gtime(dectim)
  dectim = dectim - t0
  if (ierr /= 0) then
   write ( *, 620) ierr
   write (imeas,620) ierr
   stop
  end if
  write (imeas,830) nvc,nface,nvert,npolh,npf,dectim
!
!     Further subdivide convex polyhedra based on mesh distribution
!     function, and determine tetrahedron size for each polyhedron.
!
  call gtime(t0)
  call eqdis3(hflag,umdf3,kappa,angacc,angedg,dmin,nmin,ntetd, &
     nsflag,nvc,nface,nvert,npolh,npf,maxvc,maxfp,maxfv,maxhf,maxpf, &
     maxiw,maxwk,vcl,facep,factyp,nrml,fvl,eang,hfl,pfl,vol,psi,h, &
     iwk,wk,ierr)

  call gtime(eqdtim)
  eqdtim = eqdtim - t0
  if (ierr /= 0) then
   write ( *, 620) ierr
   write (imeas,620) ierr
   stop
  end if
  angmin = 2.0D+00* d_pi ( )
  do 60 i = 1,nvert
   if (eang(i) > -1.0D+00) angmin = min(angmin,eang(i))
   60 continue
  angmin = angmin*180.0D+00/ d_pi ( )
  write (imeas,840) nvc,nface,nvert,npolh,npf,angmin,eqdtim
!
!     Generate mesh vertices on a quasi-uniform grid in each convex
!     polyh, and triangulate them using boundary-constrained triang.
!
  call gtime(t0)
  call tripr3(h,shrf,crit,nvc,nface,nvert,npolh,maxvc,maxbt, &
     maxfc,maxht,maxvm,maxiw,maxwk,vcl,facep,nrml,fvl,eang,hfl,pfl, &
     edst,edno,fcst,btst,btl,fc,ht,vm,xfhv,ntrif,ntetra,iwk,wk, ierr)
  call gtime(tritim)
  tritim = tritim - t0
  if (ierr /= 0) then
   write ( *, 620) ierr
   write (imeas,620) ierr
   stop
  end if
  nbt = btst(nface+1)-1
  ntri = 0
  nteta = 0
  nt = 0
  do 70 i = 1,npolh
   ntri = ntri + ntrif(i)
   nteta = nteta + ntetra(i)
   nt = max(nt,ntetra(i))
   70 continue
  if (npolh > maxiw) then
   write ( *, 610) 'maxiw',npolh
   stop
  else if (nt > maxte) then
   write ( *, 610) 'maxte',nt
   stop
  else if (nteta > maxnt) then
   write ( *, 610) 'maxnt',nteta
   stop
  end if
  if (outmod < 0) then
   write (utet,730) nvc,nteta
   write (utet,740) ((vcl(j,i),j=1,3),i=1,nvc)
   write (utet,730)
  end if
!
!     Collect measurements for number of interior faces that fail local
!     sphere test, and min, max, mean, stdv of relative volume, solid
!     angle, norm. ratio of inradius to circumradius for all tetrahedra.
!
  nit = 0
  nlo = 0
  kt = 0
  sf = 1.5D+00*sqrt(6.0D+00)
  do 110 p = 1,npolh
   k = 0
   j = hfl(p)
   80    continue
      f = abs(pfl(1,j))
      k = k + (btst(f+1) - btst(f))
      j = pfl(2,j)
   if (j /= hfl(p)) go to 80
   nit = nit + (ntrif(p) - k)
   kfc = xfhv(1,p) + k
   iwk(p) = kfc
   kvm = xfhv(3,p) - 1
   do 90 i = kfc,xfhv(1,p+1)-1
      if (fc(1,i) > 0) then
         a = vm(fc(1,i)+kvm)
         b = vm(fc(2,i)+kvm)
         c = vm(fc(3,i)+kvm)
         d = vm(fc(4,i)+kvm)
         e = vm(fc(5,i)+kvm)
         call ccsph(.true.,vcl(1,a),vcl(1,b),vcl(1,c),vcl(1,d), &
              vcl(1,e),sa,radsq,in)
         if (in >= 1) nlo = nlo + 1
      end if
   90    continue
   call tetlst(xfhv(1,p+1)-xfhv(1,p),vm(xfhv(3,p)), &
        fc(1,xfhv(1,p)),nt,tetra)
   if (nt /= ntetra(p)) then
      write ( *, 750) p,nt,ntetra(p)
      write (imeas,750) p,nt,ntetra(p)
      stop
   end if
   if (outmod < 0) write (utet,730) ((tetra(j,i),j=1,4),i=1,nt)
   h3 = 1.0D+00/h(p)**3
   do i = 1,nt
      kt = kt + 1
      a = tetra(1,i)
      b = tetra(2,i)
      c = tetra(3,i)
      d = tetra(4,i)
      rvol(kt) = h3*volth(vcl(1,a),vcl(1,b),vcl(1,c),vcl(1,d))
      rho(kt) = radrth(vcl(1,a),vcl(1,b),vcl(1,c),vcl(1,d))
      sig(kt) = sangmn(vcl(1,a),vcl(1,b),vcl(1,c),vcl(1,d),sa)*sf
      eta(kt) = emnrth(vcl(1,a),vcl(1,b),vcl(1,c),vcl(1,d))
  end do
  110 continue
  write (imeas,850) fcst(1)-1,fcst(nface+1)-1,nvc,nbt, &
     (xfhv(i,npolh+1)-1,i=1,3),ntri,nteta,nlo,tritim
!
  call stats(nteta,rvol,av,hv,nf,vmin,vmax,vmean,vstdv,vfreq)
  call stats(nteta,rho,am,hm,nf,rmin,rmax,rmean,rstdv,rfreq)
  call stats(nteta,sig,am,hm,nf,smin,smax,smean,sstdv,sfreq)
  call stats(nteta,eta,am,hm,nf,emin,emax,emean,estdv,efreq)
  write (imeas,860) vmin,vmax,vmean,vstdv
  write (imeas,870) rmin,rmax,rmean,rstdv
  write (imeas,880) smin,smax,smean,sstdv
  write (imeas,890) emin,emax,emean,estdv
  write (imeas,900)

  do i = 0,nf
    write (imeas,910) av+i*hv,vfreq(i),am+i*hm,rfreq(i),am+i*hm, &
        sfreq(i),am+i*hm,efreq(i)
  end do

  if (abs(outmod) == 1) then
     write ( *, 630) nvc,nface,nvert,npolh,npf
     write ( *, 640) (i,(vcl(j,i),j=1,3),i=1,nvc)
     write ( *, 650) (i,(facep(j,i),j=1,3),factyp(i), &
        (nrml(j,i),j=1,3),fcst(i),btst(i),i=1,nface), nface+1, &
        0,0,0,0,0.0D+00,0.0D+00,0.0D+00,fcst(nface+1),btst(nface+1)
     write ( *, 660) (i,(fvl(j,i),j=1,6),eang(i),edst(i),edno(i), &
        i=1,nvert)
     write ( *, 670) (i,hfl(i),(xfhv(j,i),j=1,3),ntrif(i), &
        ntetra(i),vol(i),psi(i),h(i),i=1,npolh), &
        npolh+1,0,(xfhv(j,npolh+1),j=1,3)
     write ( *, 680) (i,(pfl(j,i),j=1,2),i=1,npf)
     write ( *, 690) (i,(btl(j,i),j=1,3),i=1,nbt)
     write ( *, 700) (i,(fc(j,i),j=1,7),i=1,xfhv(1,npolh+1)-1)
     write ( *, 710) (ht(i),i=1,xfhv(2,npolh+1)-1)
     write ( *, 720) (vm(i),i=1,xfhv(3,npolh+1)-1)
  else if (abs(outmod) == 2) then
   write ( *, 730) nvc,nbt,nit
   write ( *, 740) ((vcl(j,i),j=1,3),i=1,nvc)
   write ( *, '(a)')
   write ( *, 730) (3,(btl(j,i),j=1,3),i=1,nbt)
   write ( *, '(a)')
   do i = 1,npolh
      kvm = xfhv(3,i) - 1
      do j = iwk(i),xfhv(1,i+1)-1
         if (fc(1,j) > 0) &
              write ( *, 730) 3,(vm(fc(k,j)+kvm),k=1,3)
      end do
    end do
  else if (abs(outmod) == 3) then
   nbt = 0
   nit = 0
   do 150 f = 1,nface
      if (factyp(f) <= 0) go to 150
      if (facep(3,f) == 0) then
         nbt = nbt + (btst(f+1) - btst(f))
      else
         nit = nit + (btst(f+1) - btst(f))
      end if
  150    continue
   write ( *, 730) nvc,nbt,nit
   write ( *, 740) ((vcl(j,i),j=1,3),i=1,nvc)
   write ( *, '(a)' )
   do 160 f = 1,nface
      if (factyp(f) <= 0 .or. facep(3,f) /= 0) go to 160
      if (facep(2,f) > 0) then
         write ( *, 730) (3,(btl(j,i),j=1,3),i=btst(f), &
              btst(f+1)-1)
      else
         write ( *, 730) (3,btl(1,i),btl(3,i),btl(2,i), &
              i=btst(f),btst(f+1)-1)
      end if
  160    continue
   write ( *, '(a)' )

   do f = 1,nface
      if (factyp(f) <= 0 .or. facep(3,f) == 0) then
        cycle
      end if
      if (facep(2,f) > 0) then
         write ( *, 730) (3,(btl(j,i),j=1,3),i=btst(f), &
              btst(f+1)-1)
      else
         write ( *, 730) (3,btl(1,i),btl(3,i),btl(2,i), &
              i=btst(f),btst(f+1)-1)
      end if
    end do

  end if

  600 format (a60)
  610 format (1x,'*** ',a,' must be increased to',i8)
  620 format (/1x,'ierr=',i5)
  630 format (1x,'nvc,nface,nvert,npolh,npf'/1x,5i7)
  640 format (/' vcl'/(1x,i7,3f23.15))
  650 format (/' facep,factyp,nrml,fcst,btst'/(1x,5i5,3f14.6,2i5))
  660 format (/' fvl,eang,edst,edno'/(1x,7i5,f15.7,2i7))
  670 format (/' hfl,xfhv,ntrif,ntetra,vol,psi,h'/(1x,2i5,i6,4i5, &
     3d13.6))
  680 format (/' pfl'/(1x,3i7))
  690 format (/' btl'/(1x,4i7))
  700 format (/' fc'/(1x,8i7))
  710 format (/' ht'/(1x,10i7))
  720 format (/' vm'/(1x,10i7))
  730 format (1x,4i7)
  740 format (/(1x,3f15.7))
  750 format (/1x,'error from tetlst polyh',i5,',   nt != ntetra :',2i7)
  800 format (1x,a60/1x,'input : tol=',d15.7,'   aspc2d=',f9.3, &
     '   atol2d=',f9.3/9x,'angacc=',f9.3,'   rdacc=',f9.5, &
     '   nfhol=',i3,'   nihol=',i3/9x,'angedg=',f9.3,'   kappa=', &
     f7.3,'   dmin=',f6.3,'   nmin=',i5/9x,'ntetd=',i7, &
     '   crit=',i2,'   shrf=',f9.5)
  810 format (1x,a6,': nvc=',i7,'   nface=',i7,'   nvert=',i7, &
     '   npolh=',i7/9x,'npf=',i7,'   nrfe=',i7,'   angmin=',f9.3, &
     3x,a6,'=',f9.3)
  820 format (1x,'decomp: nvc=',i7,'   nface=',i7,'   nvert=',i7, &
     '   npolh=',i7/9x,'npf=',i7,'   nrfe=',i7,'   angmin=',f9.3, &
     '   dectim=',f9.3)
  830 format (1x,'decfac: nvc=',i7,'   nface=',i7,'   nvert=',i7, &
     '   npolh=',i7/9x,'npf=',i7,'   dectim=',f9.3)
  840 format (1x,'eqdist: nvc=',i7,'   nface=',i7,'   nvert=',i7, &
     '   npolh=',i7/9x,'npf=',i7,'   angmin=',f9.3, &
     '   eqdtim=',f9.3)
  850 format (1x,'triang: nvce=',i7,'   nvcf=',i7,'   nvc=',i7, &
     '   nbt=',i7/9x,'nfc=',i7,'   nht=',i7,'   nvm=',i7,'   ntri=', &
     i7/9x,'nteta=',i7,'   nlo=',i7,'   tritim=',f9.3)
  860 format (1x,'vmin=',f11.7,3x,'vmax=',f11.7,3x,'vmean=',f11.7, &
     3x,'vstdv=',f11.7)
  870 format (1x,'rmin=',f11.7,3x,'rmax=',f11.7,3x,'rmean=',f11.7, &
     3x,'rstdv=',f11.7)
  880 format (1x,'smin=',f11.7,3x,'smax=',f11.7,3x,'smean=',f11.7, &
     3x,'sstdv=',f11.7)
  890 format (1x,'emin=',f11.7,3x,'emax=',f11.7,3x,'emean=',f11.7, &
     3x,'estdv=',f11.7)
  900 format (/3x,'rel volume',9x,'radius ratio',8x,'min solid ang',8x, &
     'mean ratio')
  910 format (1x,f5.2,f9.4,6x,f5.2,f9.4,6x,f5.2,f9.4,6x,f5.2,f9.4)

  return
end
subroutine test_vnbr ( )

!*****************************************************************************80
!
!! TEST_VNBR tests VISVRT, VORNBR.
!
  implicit none

  integer ( kind = 4 ) ierr

  integer ( kind = 4 ), parameter :: maxn = 200

  integer ( kind = 4 ) i
  integer ( kind = 4 ) irdr,ivert,j,maxnv,n,nvis,nvor,nvrt,nvsvrt
  integer ( kind = 4 ) ivis(0:maxn),ivor(0:maxn)
  real ( kind = 8 ) angle,angspc,phi,tol,xeye,yeye
  real ( kind = 8 ) temp
  real ( kind = 8 ) theta(0:maxn),x(maxn),xc(0:maxn),xvor(0:maxn)
  real ( kind = 8 ) y(maxn),yc(0:maxn),yvor(0:maxn)

  irdr = 5
  read (irdr,*) n,tol
  if (n < 3) stop

  if (n > maxn) then
    write ( *, 600) 'maxn',n
    stop
  end if

  read (irdr,*) (x(i),y(i),i=1,n)
  read (irdr,*) ivert
  read (irdr,*) angspc
!
!     1 <= IVERT <= N is index of polygon vertex for eyepoint.
!     ANGSPC is angle spacing parameter in degrees.
!
  angspc = angspc*acos(-1.0D+00)/180.0D+00
  xeye = x(ivert)
  yeye = y(ivert)
  nvrt = n - 2
  j = -1
  do i = ivert+1,n
    j = j + 1
    xc(j) = x(i)
    yc(j) = y(i)
  end do

  do i = 1,ivert-1
    j = j + 1
    xc(j) = x(i)
    yc(j) = y(i)
  end do

  write ( *, 620) xeye,yeye
  write ( *, 630) nvrt,(xc(i),yc(i),i=0,nvrt)
  call vispol(xeye,yeye,nvrt,xc,yc,nvis,ivis, ierr )
  write ( *, 640) nvis,(xc(i),yc(i),ivis(i),i=0,nvis)
  if (ierr /= 0) go to 30

  phi = angle(xc(nvis),yc(nvis),xeye,yeye,xc(0),yc(0))
  temp = phi / angspc
  maxnv = nvis + int ( temp )

  if (maxnv > maxn) then
    write ( *, 600) 'maxn',maxnv
    stop
  end if

  call visvrt(angspc,xeye,yeye,nvis,xc,yc,ivis,maxnv,nvsvrt,theta)
  write ( *, 650) nvsvrt,(xc(i),yc(i),ivis(i),theta(i),i=0,nvsvrt)
  nvor = -1
  call vornbr(xeye,yeye,nvsvrt,xc,yc,nvor,ivor,xvor,yvor, ierr )
  write ( *, 640) nvor,(xvor(i),yvor(i),ivor(i),i=0,nvor)

   30 continue

  if (ierr /= 0) write ( *, 610) ierr

  600 format (1x,'*** ',a,' must be increased to',i8)
  610 format (/1x,'ierr=',i5)
  620 format (1x,2f15.7)
  630 format (/1x,i5/(1x,2f15.7))
  640 format (/1x,i5/(1x,2f15.7,i5))
  650 format (/1x,i5/(1x,2f15.7,i5,f15.7))

  return
end
subroutine test_vpol ( )

!*****************************************************************************80
!
!! TEST_VPOL tests VISPOL and ROTIPG.
!
  implicit none

  integer ( kind = 4 ) ierr

  integer ( kind = 4 ), parameter :: maxn = 200

  integer ( kind = 4 ) i
  integer ( kind = 4 ) irdr,ivert,j,n,nvis,nvrt,vptype
  integer ( kind = 4 ) ivis(0:maxn+1)
  real ( kind = 8 ) tol,xeye,yeye
  real ( kind = 8 ) x(maxn)
  real ( kind = 8 ) xc(0:maxn+1),y(maxn),yc(0:maxn+1)

  irdr = 5
  read (irdr,*) n,tol
  if (n < 3) stop
  if (n > maxn) then
    write ( *, 600) 'maxn',n
    stop
  end if

  read (irdr,*) (x(i),y(i),i=1,n)
  read (irdr,*) ivert,vptype,xeye,yeye
!
!     1 <= IVERT <= N is index of polygon vertex or IVERT <= 0 if
!     ROTIPG is to be called. BVTYPE = 0, 1, or 2 for boundary,
!     interior, or blocked exterior viewpoint. (XEYE,YEYE) is needed
!     for non-boundary viewpoints only, and must be visible from
!     (X(IVERT),Y(IVERT)) if IVERT > 0.
!
  if (vptype == 0) then
   xeye = x(ivert)
   yeye = y(ivert)
   nvrt = n - 2
   j = -1
   do i = ivert+1,n
      j = j + 1
      xc(j) = x(i)
      yc(j) = y(i)
   end do
   do i = 1,ivert-1
      j = j + 1
      xc(j) = x(i)
      yc(j) = y(i)
    end do
  else if (vptype == 1) then
   nvrt = n
   if (ivert <= 0) then
      do i = 1,n
         xc(i) = x(i)
         yc(i) = y(i)
      end do
      xc(0) = xc(n)
      yc(0) = yc(n)
      call rotipg(xeye,yeye,nvrt,xc,yc, ierr )
   else
      j = -1

      do i = ivert,n
         j = j + 1
         xc(j) = x(i)
         yc(j) = y(i)
      end do

      do i = 1,ivert
         j = j + 1
         xc(j) = x(i)
         yc(j) = y(i)
      end do

   end if
  else
   nvrt = n
   if (ivert <= 0) then
      do i = 1,n
         xc(n-i) = x(i)
         yc(n-i) = y(i)
      end do
      xc(n) = xc(0)
      yc(n) = yc(0)
      call rotipg(xeye,yeye,nvrt,xc,yc, ierr )
   else
      j = -1
      do i = ivert,1,-1
         j = j + 1
         xc(j) = x(i)
         yc(j) = y(i)
      end do
      do i = n,ivert,-1
         j = j + 1
         xc(j) = x(i)
         yc(j) = y(i)
      end do
   end if
  end if
  if (ierr /= 0) go to 90
  write ( *, 620) xeye,yeye
  write ( *, 630) nvrt,(xc(i),yc(i),i=0,nvrt)
  call vispol(xeye,yeye,nvrt,xc,yc,nvis,ivis, ierr )
  write ( *, 640) nvis,(xc(i),yc(i),ivis(i),i=0,nvis)
   90 continue
  if (ierr /= 0) write ( *, 610) ierr

  600 format (1x,'*** ',a,' must be increased to',i8)
  610 format (/1x,'ierr=',i5)
  620 format (1x,2f15.7)
  630 format (/1x,i5/(1x,2f15.7))
  640 format (/1x,i5/(1x,2f15.7,i5))

  return
end
