program main

!*****************************************************************************80
!
!! MAIN is the main program for GEOMPACK2_PRB.
!
!  Discussion:
!
!    GEOMPACK2_PRB runs the GEOMPACK2 tests.
!
!  Modified:
!
!    25 February 2007
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GEOMPACK2_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the GEOMPACK2 library.'

  call test01 ( )
  call test02 ( )
  call test03 ( 'cmos.in' )
  call test04 ( )
  call test05 ( 'annulus.in' )
  call test06 ( )
  call test07 ( )
  call test08 ( 'annulus.in' )
  call test09 ( )

  call test10 ( )
  call test11 ( 'annulus.in' )
  call test12 ( )
  call test13 ( 'ptpg.in' )
  call test14 ( )
  call test15 ( 'shr1.in' )
  call test16 ( )
  call test17 ( 'annulus.in' )
  call test18 ( )
  call test19 ( )

  call test20 ( )
  call test21 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GEOMPACK2_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests ANGLE, AREAPG, AREATR.
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 9

  real ( kind = 8 ) ang(n)
  real ( kind = 8 ) angle
  real ( kind = 8 ) area
  real ( kind = 8 ) areapg
  real ( kind = 8 ) areatr
  real ( kind = 8 ) atr(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) sum2
  real ( kind = 8 ), dimension ( n ) :: xc = (/ &
    5.0D+00, 5.0D+00, 7.0D+00, 9.0D+00, 5.0D+00, &
    1.0D+00, 2.0D+00, 0.0D+00, 3.0D+00 /)
  real ( kind = 8 ), dimension ( n ) :: yc = (/ &
    0.0D+00, 2.0D+00, 2.0D+00, 4.0D+00, 8.0D+00, &
    7.0D+00, 5.0D+00, 3.0D+00, 3.0D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  AREAPG computes the area of a polygon;'
  write ( *, '(a)' ) '  AREATR computes the area of a triangle;'
  write ( *, '(a)' ) '  ANGLE computes the polygonal angle at any vertex.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of polygonal vertices is N = ', n
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I, XC(I), YC(I)'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(i5,2d15.7)' ) i, xc(i), yc(i)
  end do

  area = areapg ( n, xc, yc )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Area computed directly by AREAPG = ', area

  sum2 = 0.0D+00

  do i = 1, n

   if ( i == 1 ) then
     atr(i) = areatr ( xc(1), yc(1), xc(i), yc(i), xc(i+1), yc(i+1) )
     ang(i) = angle ( xc(n), yc(n), xc(i), yc(i), xc(i+1), yc(i+1) )
   else if ( i == n ) then
     atr(i) = areatr ( xc(1), yc(1), xc(i), yc(i), xc(1), yc(1) )
     ang(i) = angle ( xc(i-1), yc(i-1), xc(i), yc(i), xc(1), yc(1) )
   else
     atr(i) = areatr ( xc(1), yc(1), xc(i), yc(i), xc(i+1), yc(i+1) )
     ang(i) = angle ( xc(i-1), yc(i-1), xc(i), yc(i), xc(i+1), yc(i+1) )
   end if

   sum2 = sum2 + atr(i)

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'I, AREATR(I), ANGLE(I)'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(i5,2d15.7)' ) i, atr(i), ang(i)
  end do

  write ( *, '(a,g14.6)' ) &
    'Area computed indirectly by summing AREATR =  ', sum2

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests CMCIRC.
!
  implicit none

  integer ( kind = 4 ) cmcirc
  integer ( kind = 4 ) in
  real ( kind = 8 ), parameter :: x0 = 3.0D+00
  real ( kind = 8 ), parameter :: y0 =  3.0D+00
  real ( kind = 8 ), parameter :: x1 =  5.0D+00
  real ( kind = 8 ), parameter :: y1 =  0.0D+00
  real ( kind = 8 ), parameter :: x2 =  0.0D+00
  real ( kind = 8 ), parameter :: y2 =  5.0D+00
  real ( kind = 8 ), parameter :: x3 = -5.0D+00
  real ( kind = 8 ), parameter :: y3 =  0.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  CMCIRC determines if a point lies in, on'
  write ( *, '(a)' ) '  or outside a circle given by 3 points.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The points defining the circle:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,2g14.6)' ) '    X1,Y1 = ', x1, y1
  write ( *, '(a,2g14.6)' ) '    X2,Y2 = ', x2, y2
  write ( *, '(a,2g14.6)' ) '    X3,Y3 = ', x3, y3
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The point to be tested:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,2g14.6)' ) '    X0,Y0 = ', x0, y0

  in = cmcirc ( x0, y0, x1, y1, x2, y2, x3, y3 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test results:'

  if ( in == 2 ) then
    write ( *, '(a)' ) '    The three points are collinear.'
  else if ( in == 1 ) then
    write ( *, '(a)' ) '    The point is inside the circle.'
  else if ( in == 0 ) then
    write ( *, '(a)' ) '    The point is on the circle.'
  else if ( in == -1 ) then
    write ( *, '(a)' ) '    The point is outside the circle.'
  end if

  return
end
subroutine test03 ( filename )

!*****************************************************************************80
!
!! TEST03 tests CVDEC2, FNDSEP, INSED2, INSVR2, JNHOLE, MINANG, RESVRT, SPDEC2.
!
  implicit none

  integer ( kind = 4 ), parameter :: incr = 10000
  integer ( kind = 4 ), parameter :: inunit = 1
  integer ( kind = 4 ), parameter :: maxed = 101
  integer ( kind = 4 ), parameter :: maxhv = 200
  integer ( kind = 4 ), parameter :: maxiw = 900
  integer ( kind = 4 ), parameter :: maxnc = 30
  integer ( kind = 4 ), parameter :: maxnv = 400
  integer ( kind = 4 ), parameter :: maxpv = 1000
  integer ( kind = 4 ), parameter :: maxvc = 500
  integer ( kind = 4 ), parameter :: maxwk = 1500

  real ( kind = 8 ) angmin
  real ( kind = 8 ) angspc
  real ( kind = 8 ) angtol
  integer ( kind = 4 ) case
  real ( kind = 8 ) degrees_to_radians
  integer ( kind = 4 ) edge(4,maxed)
  character ( len = * ) filename
  integer ( kind = 4 ) ht(0:maxed-1)
  integer ( kind = 4 ) htsiz
  integer ( kind = 4 ) hvl(maxhv)
  integer ( kind = 4 ) i
  real ( kind = 8 ) iang(maxpv)
  integer ( kind = 4 ) icur(maxnc)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ivrt(maxnv)
  integer ( kind = 4 ) iwk(maxiw)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) map(maxnc)
  integer ( kind = 4 ) msglvl
  integer ( kind = 4 ) ncur
  integer ( kind = 4 ) nh
  integer ( kind = 4 ) nhola
  integer ( kind = 4 ) nhole
  integer ( kind = 4 ) npolg
  integer ( kind = 4 ) nrfv
  integer ( kind = 4 ) nsc
  integer ( kind = 4 ) nv
  integer ( kind = 4 ) nvbc(maxnc)
  integer ( kind = 4 ) nvc
  integer ( kind = 4 ) nvcin
  integer ( kind = 4 ) nvert
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ) prime
  integer ( kind = 4 ) pvl(4,maxpv)
  real ( kind = 8 ) radians_to_degrees
  integer ( kind = 4 ) regnum(maxhv)
  character ( len = 20 ) rgname
  real ( kind = 8 ) tol
  real ( kind = 8 ) tolin
  real ( kind = 8 ) vcl(2,maxvc)
  real ( kind = 8 ) wk(maxwk)

  tol = 100.0D+00 * epsilon ( tol )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  Test CVDEC.'
  write ( *, '(a)' ) '  Test FNDSEP.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Reading input file: "' // trim ( filename ) // '".'

  open ( unit = inunit, file = filename, form = 'formatted' )
!
!  Read in vertices of general polygonal region.
!
!  CASE = 1 : simple polygon or multiply connected polygonal region
!  CASE = 2 : general polygonal region with holes and interfaces
!  CASE = -1 or -2: similar, but decompose into simple polygons only
!     and don't obtain convex polygon decomposition
!
  write ( *, '(a)' ) ' '

  read ( inunit, '(a)' ) rgname
  write ( *, '(a)' ) trim ( rgname )

  read ( inunit, * ) tolin, angspc, angtol

  write ( *, '(a,d15.7)' ) '  TOLIN =  ', tolin
  write ( *, '(a,d15.7)' ) '  ANGSPC = ', angspc
  write ( *, '(a,d15.7)' ) '  ANGTOL = ', angtol

  angspc = degrees_to_radians ( angspc )
  angtol = degrees_to_radians ( angtol )

  read ( inunit, * ) case, nvc, ncur, msglvl

  if ( maxvc < nvc ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST03 - Error!'
    write ( *, '(a)' ) '  MAXVC < NVC.'
    return
  end if

  if ( maxnc < ncur ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST03 - Error!'
    write ( *, '(a)' ) '  MAXNC < NCUR.'
    return
  end if

  read ( inunit, * ) nvbc(1:ncur)
  if ( abs(case) == 2 ) then
    read ( inunit, * ) icur(1:ncur)
  end if

  read ( inunit, * ) ( vcl(1,i), vcl(2,i), i = 1, nvc )
!
!  Call DSMCPR or DSPGDC to set data structures in arrays
!  REGNUM, HVL, PVL, IANG, HOLV = IWK.
!
  if ( abs ( case ) == 1 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Calling DSMCPR to set data structures.'

    nhole = ncur - 1

    call dsmcpr ( nhole, nvbc, vcl, maxhv, maxpv, maxiw, nvc, npolg, nvert, &
      nhola, regnum, hvl, pvl, iang, iwk, ierror )

  else if ( abs(case) == 2 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Calling DSPGDC to set data structures.'

    nv = sum ( nvbc(1:ncur) )

    if ( maxnv < nv ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03 - Error!'
      write ( *, '(a)' ) '  MAXNV < NV.'
      return
    end if

    read ( inunit, * ) ivrt(1:nv)
    nsc = 0

    do i = 1, nv
      if ( ivrt(i) < 0 ) then
        nsc = nsc + 1
      end if
    end do

    nsc = nsc / 2

    if ( maxed < nsc ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03 - Fatal error!'
      write ( *, '(a)' ) '  MAXED < NSC.'
      return
    end if

    htsiz = min ( prime ( nsc / 2 ), maxed )

    call dspgdc ( nvc, vcl, incr, ncur, nvbc, icur, ivrt, maxhv, maxpv, maxiw, &
      npolg, nvert, nhole, nhola, regnum, hvl, pvl, iang, iwk, htsiz, nsc, &
      ht, edge, map, ierror )

  end if

  close ( unit = inunit )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST03 - Error!'
    write ( *, '(a)' ) '  IERROR = ', ierror
    return
  end if
!
!  Print out data structures and measurements after initialization.
!
  nvcin = nvc
  nh = nhole * 2 + nhola
  write ( *, '(a,i6)' ) 'MSGLVL = ', msglvl
  write ( *, '(a,i6)' ) 'NVC = ', nvc
  write ( *,630) (i,vcl(1,i),vcl(2,i),i=1,nvc)
  write ( *,640) npolg,nhole,(hvl(i),i=1,npolg+nhole)
  write ( *,650) (regnum(i),i=1,npolg)
  write ( *,660) nvert,(i,(pvl(j,i),j=1,4),iang(i),i=1,nvert)
  write ( *,640) nhola,nh,(iwk(i),i=1,nh)

  nrfv = 0
  do i = 1, nvert
    if ( pi + tol < iang(i) ) then
      nrfv = nrfv + 1
    end if
  end do

  angmin = iang(1)
  do i = 2, nvert
    angmin = min ( angmin, iang(i) )
  end do
  angmin = radians_to_degrees ( angmin )

  write (*,710) nvc,npolg,nvert,nhole,nhola,nrfv,angmin
!
!  Obtain simple (and convex) polygon decompositions.
!
  if ( msglvl == 2 ) then
    write ( *,670)
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Calling SPDEC2'

  call spdec2 ( angspc, angtol, nvc, npolg, nvert, nhole, nhola, maxvc, &
    maxhv, maxpv, maxiw-nh, maxwk, iwk, vcl, regnum, hvl, pvl, iang, &
    iwk(nh+1), wk, ierror )

  if ( 0 < case ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Calling CVDEC2'

    call cvdec2 ( angspc, angtol, nvc, npolg, nvert, maxvc, maxhv, maxpv, &
      maxiw, maxwk, vcl, regnum, hvl, pvl, iang, iwk, wk, ierror )

  end if

  if ( msglvl == 2 ) then
    write ( *,670) 0, 0, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00
  end if

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST03 - Error!'
    write ( *, '(a)' ) '  IERROR = ', ierror
    return
  end if
!
!  Print out data structures and measurements after decomposition.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  NVC = ', nvc
  write ( *,630) (i,vcl(1,i),vcl(2,i),i=nvcin+1,nvc)
  write ( *,680) npolg, hvl(1:npolg)
  write ( *,650) regnum(1:npolg)
  write ( *,660) nvert,(i,(pvl(j,i),j=1,4),iang(i),i=1,nvert)

  angmin = minval ( iang(1:nvert) )
  angmin = radians_to_degrees ( angmin )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DECOMP:'
  write ( *, '(a,i6)' ) '  NVC =    ', nvc
  write ( *, '(a,i6)' ) '  NPOLG =  ', npolg
  write ( *, '(a,i6)' ) '  NVERT =  ', nvert
  write ( *, '(a,g14.6)' ) '  ANGMIN = ', angmin

  630 format ((1x,i7,2f15.7))
  640 format (/1x,2i7/(1x,10i7))
  650 format (/(1x,10i7))
  660 format (/1x,i7/(1x,5i7,f15.7))
  670 format (1x,2i7,4f15.7)
  680 format (/1x,i7/(1x,10i7))

  710 format (1x,'initds: nvc=',i7,'   npolg=',i7,'   nvert=',i7, &
    '   nhole=',i7/9x,'nhola=',i7,'   nrfv=',i7,'   angmin=',f9.3 )

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests DIAEDG.
!
  implicit none

  integer ( kind = 4 ) diaedg
  integer ( kind = 4 ) in
  real ( kind = 8 ) x0
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ) y0
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) y3

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  DIAEDG determines which diagonal of a'
  write ( *, '(a)' ) '  quadrilateral is to be preferred, based on'
  write ( *, '(a)' ) '  the circumcircle criterion.'

  x0 =  0.0D+00
  y0 =  0.0D+00
  x1 =  5.0D+00
  y1 =  0.0D+00
  x2 =  6.0D+00
  y2 =  1.0D+00
  x3 =  1.0D+00
  y3 =  1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The points defining the quadrilateral:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,2g14.6)' ) '  P0: X0,Y0 = ', x0, y0
  write ( *, '(a,2g14.6)' ) '  P1: X1,Y1 = ', x1, y1
  write ( *, '(a,2g14.6)' ) '  P2: X2,Y2 = ', x2, y2
  write ( *, '(a,2g14.6)' ) '  P3: X3,Y3 = ', x3, y3

  in = diaedg ( x0, y0, x1, y1, x2, y2, x3, y3 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  DIAEDG results:'
  write ( *, '(a)' ) ' '

  if ( in == 1 ) then
    write ( *, '(a)' ) '  Use diagonal P0--P2.'
  else if ( in == -1 ) then
    write ( *, '(a)' ) '  Use diagonal P1--P3.'
  else if ( in == 0 ) then
    write ( *, '(a)' ) '  All 4 points lie on a circle.'
    write ( *, '(a)' ) '  Either diagonal can be used.'
  end if

  return
end
subroutine test05 ( filename )

!*****************************************************************************80
!
!! TEST05 tests DSMCPR, DSPGDC, EDGHT, HOLVRT.
!
  implicit none

  integer ( kind = 4 ), parameter :: incr = 10000
  integer ( kind = 4 ), parameter :: maxed = 101
  integer ( kind = 4 ), parameter :: maxho = 50
  integer ( kind = 4 ), parameter :: maxnc = 30
  integer ( kind = 4 ), parameter :: maxnv = 500

  integer ( kind = 4 ) case
  integer ( kind = 4 ) edge(4,maxed)
  character ( len = * ) filename
  integer ( kind = 4 ) holv(maxho)
  integer ( kind = 4 ) ht(0:maxed-1)
  integer ( kind = 4 ) htsiz
  integer ( kind = 4 ) hvl(maxnc*2)
  integer ( kind = 4 ) i
  real ( kind = 8 ) iang(maxnv*2)
  integer ( kind = 4 ) icur(maxnc)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ), parameter :: inunit = 1
  integer ( kind = 4 ) ivrt(maxnv)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) map(maxnc)
  integer ( kind = 4 ) ncur
  integer ( kind = 4 ) nhola
  integer ( kind = 4 ) nhole
  integer ( kind = 4 ) npolg
  integer ( kind = 4 ) nsc
  integer ( kind = 4 ) nv
  integer ( kind = 4 ) nvbc(maxnc)
  integer ( kind = 4 ) nvc
  integer ( kind = 4 ) nvert
  integer ( kind = 4 ) prime
  integer ( kind = 4 ) pvl(4,maxnv*2)
  integer ( kind = 4 ) regnum(maxnc*2)
  character ( len = 20 ) rgname
  real ( kind = 8 ) tolin
  real ( kind = 8 ) vcl(2,maxnv)
!
!  Read in the vertices of a general polygonal region.
!  CASE = 1 : simple polygon or multiply connected polygonal region
!  CASE = 2 : general polygonal region with holes and interfaces
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  Reading input file:'
  write ( *, '(a)' ) filename

  open ( unit = inunit, file = filename, form = 'formatted' )

  read ( inunit, '(a)' ) rgname
  read ( inunit, * ) tolin
  read ( inunit, * ) case, nvc, ncur

  if ( maxnv < nvc ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST05 - Error!'
    write ( *, '(a)' ) '  MAXNV < NVC.'
    return
  end if

  if ( maxnc < ncur ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST05 - Error!'
    write ( *, '(a)' ) '  MAXNC < NCUR.'
    return
  end if

  read ( inunit, * ) nvbc(1:ncur)
  if ( case == 2) then
    read ( inunit, * ) icur(1:ncur)
  end if
  read ( inunit, * ) ( vcl(1,i), vcl(2,i), i = 1, nvc )
!
!  Call DSMCPR or DSPGDC to set data structures in arrays
!  REGNUM, HVL, PVL, IANG, HOLV.
!
  if ( case == 1 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Calling DSMCPR:'

    nhole = ncur - 1

    call dsmcpr ( nhole, nvbc, vcl, maxnc*2, maxnv*2, maxho, nvc, &
      npolg, nvert, nhola, regnum, hvl, pvl, iang, holv, ierror )

  else if ( case == 2 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Calling DSPGDC:'

    nv = sum ( nvbc(1:ncur) )

    if ( maxnv < nv ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST05 - Fatal error!'
      write ( *, '(a)' ) '  MAXNV < NV.'
      return
    end if

    read ( inunit, * ) ivrt(1:nv)

    nsc = 0
    do i = 1, nv
      if ( ivrt(i) < 0 ) then
        nsc = nsc + 1
      end if
    end do

    nsc = nsc / 2

    if ( maxed < nsc ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST05 - Fatal error!'
      write ( *, '(a)' ) '  MAXED < NSC.'
      return
    end if

    htsiz = min ( prime(nsc/2), maxed )

    call dspgdc ( nvc, vcl, incr, ncur, nvbc, icur, ivrt, maxnc*2, maxnv*2, &
      maxho, npolg, nvert, nhole, nhola, regnum, hvl, pvl, iang, holv, &
      htsiz, nsc, ht, edge, map, ierror )

  end if

  close ( unit = inunit )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST05 - Error!'
    write ( *, '(a)' ) '  IERROR = ', ierror
    return
  end if
!
!  Print out the data structures.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  NVC = ', nvc
  write ( *, '(a)' ) ' '
  write ( *, '(i7,2f15.7)' ) (i,vcl(1:2,i),i=1,nvc)
  write ( *, '(a)' ) ' '
  write ( *,640) npolg,nhole,(hvl(i),i=1,npolg+nhole)
  write ( *, '(a)' ) ' '
  write ( *,650) regnum(1:npolg)
  write ( *, '(a)' ) ' '
  write ( *,660) nvert,(i,(pvl(j,i),j=1,4),iang(i),i=1,nvert)
  write ( *, '(a)' ) ' '
  write ( *,640) nhola,nhole*2+nhola,(holv(i),i=1,nhole*2+nhola)

  640 format (1x,2i7/(1x,10i7))
  650 format ((1x,10i7))
  660 format (1x,i7/(1x,5i7,f15.7))

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 tests DTRIS2;
!
  implicit none

  integer ( kind = 4 ), parameter :: g_num = 24

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) c
  integer ( kind = 4 ) d
  integer ( kind = 4 ) diaedg
  character ( len = 255 ) :: file_name = 'test06_triangulation_plot.eps'
  real ( kind = 8 ), dimension (2,g_num) :: g_xy = reshape ( (/ &
     1.0000000D+00,  0.0000000D+00, &
     0.9238795D+00,  0.3826834D+00, &
     0.7071068D+00,  0.7071068D+00, &
     0.3826834D+00,  0.9238795D+00, &
     0.0000000D+00,  1.0000000D+00, &
    -0.3826834D+00,  0.9238795D+00, &
    -0.7071068D+00,  0.7071068D+00, &
    -0.9238795D+00,  0.3826834D+00, &
    -1.0000000D+00,  0.0000000D+00, &
    -0.9238795D+00, -0.3826834D+00, &
    -0.7071068D+00, -0.7071068D+00, &
    -0.3826834D+00, -0.9238795D+00, &
     0.0000000D+00, -1.0000000D+00, &
     0.3826834D+00, -0.9238795D+00, &
     0.7071068D+00, -0.7071068D+00, &
     0.9238795D+00, -0.3826834D+00, &
     0.7500000D+00,  0.0000000D+00, &
     0.6767767D+00,  0.1767767D+00, &
     0.5000000D+00,  0.2500000D+00, &
     0.3232233D+00,  0.1767767D+00, &
     0.2500000D+00,  0.0000000D+00, &
     0.3232233D+00, -0.1767767D+00, &
     0.5000000D+00, -0.2500000D+00, &
     0.6767767D+00, -0.1767767D+00 /), (/ 2, g_num /) )
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ind(g_num)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jp1
  integer ( kind = 4 ) jp2
  integer ( kind = 4 ) k
  integer ( kind = 4 ) nlo
  integer ( kind = 4 ) nod_tri(3,2*g_num)
  integer ( kind = 4 ) tnbr(3,2*g_num)
  integer ( kind = 4 ) tri_num

  do i = 1, g_num
    ind(i) = i
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  DTRIS2 constructs the Delaunay triangulation'
  write ( *, '(a)' ) '  of a set of points in the plane.'

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  The number of points to triangulate is ', g_num
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The coordinates of the points are:'
  write ( *, '(a)' ) ' '
  do i = 1, g_num
    write ( *, '(i5,2f15.7)' ) i, g_xy(1,i), g_xy(2,i)
  end do

  call dtris2 ( g_num, g_xy, ind, tri_num, nod_tri, tnbr, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST06 - Error!'
    write ( *, '(a,i6)' ) '  IERROR = ', ierror
    return
  end if

  nlo = 0

  do i = 1, tri_num

    do j = 1, 3

      k = tnbr(j,i)

      if ( i < k ) then

        jp1 = j + 1
        if ( 3 < jp1 ) then
          jp1 = 1
        end if

        jp2 = jp1 + 1

        if ( 3 < jp2 ) then
          jp2 = 1
        end if

        a = nod_tri(j,i)
        b = nod_tri(jp1,i)
        c = nod_tri(jp2,i)

        if ( nod_tri(1,k) == b ) then
          d = nod_tri(3,k)
        else if ( nod_tri(2,k) == b ) then
          d = nod_tri(1,k)
        else
          d = nod_tri(2,k)
        end if

        if ( diaedg ( g_xy(1,c), g_xy(2,c), g_xy(1,a), g_xy(2,a), g_xy(1,d), &
          g_xy(2,d), g_xy(1,b), g_xy(2,b) ) == 1 ) then
          nlo = nlo + 1
        end if

      end if

    end do

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  The number of triangles = ', tri_num
  write ( *, '(a,i6)' ) '  NLO =  ', nlo

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I, NOD_TRI(1:3,I), TNBR(1:3,I)'
  write ( *, '(a)' ) ' '

  do i = 1, tri_num
    write ( *, '(7i8)' ) i, nod_tri(1:3,i), tnbr(1:3,i)
  end do

  call triangulation_plot_eps ( file_name, g_num, g_xy, tri_num, nod_tri )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  TRIANGULATION_PLOT_EPS has created an EPS file'
  write ( *, '(a)' ) '  containing an image of the triangulation, '
  write ( *, '(a)' ) '  in ' // trim ( file_name )

  return
end
subroutine test07 ( )

!*****************************************************************************80
!
!! TEST07 tests DTRIW2;
!
  implicit none

  real ( kind = 8 ), parameter :: large = 1000.0D+00
  integer ( kind = 4 ), parameter :: maxnp = 10000
  integer ( kind = 4 ), parameter :: maxst = 100

  integer ( kind = 4 ) a
  integer ( kind = 4 ) alg
  integer ( kind = 4 ) b
  real ( kind = 8 ) binexp
  integer ( kind = 4 ) c
  integer ( kind = 4 ) d
  integer ( kind = 4 ) diaedg
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ind(maxnp+3)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jp1
  integer ( kind = 4 ) jp2
  integer ( kind = 4 ) k
  integer ( kind = 4 ) msglvl
  integer ( kind = 4 ) nlo
  integer ( kind = 4 ) npt
  integer ( kind = 4 ) ntri
  integer ( kind = 4 ) stack(maxst)
  integer ( kind = 4 ) til(3,maxnp*2+1)
  integer ( kind = 4 ) tnbr(3,maxnp*2+1)
  real ( kind = 8 ) vcl(2,maxnp+3)
!
!  ALG =
!    2: DTRIW2;
!    3: DTRIW2 with bounding triangle;
!    4: DTRIW2 with call to BNSRT2 first.
!  MSGLVL
!     0: print arrays;
!     4: also print edges as they are created and swapped.
!
!  I HAVE NO IDEA HOW TO CHOOSE BINEXP
!
  msglvl = 0
  binexp = 0.5D+00

  vcl(1,1) =    1.0000000D+00
  vcl(2,1) =    0.0000000D+00
  vcl(1,2) =    0.9238795D+00
  vcl(2,2) =    0.3826834D+00
  vcl(1,3) =    0.7071068D+00
  vcl(2,3) =    0.7071068D+00
  vcl(1,4) =    0.3826834D+00
  vcl(2,4) =    0.9238795D+00
  vcl(1,5) =    0.0000000D+00
  vcl(2,5) =    1.0000000D+00
  vcl(1,6) =  - 0.3826834D+00
  vcl(2,6) =    0.9238795D+00
  vcl(1,7) =  - 0.7071068D+00
  vcl(2,7) =    0.7071068D+00
  vcl(1,8) =  - 0.9238795D+00
  vcl(2,8) =    0.3826834D+00
  vcl(1,9) =  - 1.0000000D+00
  vcl(2,9) =    0.0000000D+00
  vcl(1,10) = - 0.9238795D+00
  vcl(2,10) = - 0.3826834D+00
  vcl(1,11) = - 0.7071068D+00
  vcl(2,11) = - 0.7071068D+00
  vcl(1,12) = - 0.3826834D+00
  vcl(2,12) = - 0.9238795D+00
  vcl(1,13) =   0.0000000D+00
  vcl(2,13) = - 1.0000000D+00
  vcl(1,14) =   0.3826834D+00
  vcl(2,14) = - 0.9238795D+00
  vcl(1,15) =   0.7071068D+00
  vcl(2,15) = - 0.7071068D+00
  vcl(1,16) =   0.9238795D+00
  vcl(2,16) = - 0.3826834D+00
  vcl(1,17) =   0.7500000D+00
  vcl(2,17) =   0.0000000D+00
  vcl(1,18) =   0.6767767D+00
  vcl(2,18) =   0.1767767D+00
  vcl(1,19) =   0.5000000D+00
  vcl(2,19) =   0.2500000D+00
  vcl(1,20) =   0.3232233D+00
  vcl(2,20) =   0.1767767D+00
  vcl(1,21) =   0.2500000D+00
  vcl(2,21) =   0.0000000D+00
  vcl(1,22) =   0.3232233D+00
  vcl(2,22) = - 0.1767767D+00
  vcl(1,23) =   0.5000000D+00
  vcl(2,23) = - 0.2500000D+00
  vcl(1,24) =   0.6767767D+00
  vcl(2,24) = - 0.1767767D+00

  npt = 24

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a,i6)' ) '  MSGLVL = ', msglvl
  write ( *, '(a,i6)' ) '  NPT =    ', npt
  write ( *, '(a,g14.6)' ) '  BINEXP = ', binexp

  if ( maxnp < npt ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST07 - Fatal error!'
    write ( *, '(a)' ) '  MAXNP < NPT.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  The number of points to triangulate is ', npt
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The coordinates of the points are:'
  write ( *, '(a)' ) ' '
  do i = 1, npt
    write ( *, '(i5,2f15.7)' ) i, vcl(1,i), vcl(2,i)
  end do

  do alg = 2, 4

    npt = 24

    write ( *, '(a,i6)' ) 'ALG =    ', alg

    if ( alg /= 3 ) then

      do i = 1, npt
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
      do i = 1, npt
        ind(i+3) = i
      end do

      npt = npt + 3

    end if

    if ( alg == 4 ) then
      call bnsrt2 ( binexp, npt, vcl, ind, til, tnbr )
    end if

    call dtriw2 ( npt, maxst, vcl, ind, ntri, til, tnbr, stack, ierror )

    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST07 - Error!'
      write ( *, '(a,i6)' ) '  IERROR = ', ierror
      return
    end if

    nlo = 0

    do i = 1, ntri

      do j = 1, 3

        k = tnbr(j,i)

        if ( i < k ) then

          jp1 = j + 1
          if ( 3 < jp1 ) then
            jp1 = 1
          end if

          jp2 = jp1 + 1

          if ( 3 < jp2 ) then
            jp2 = 1
          end if

          a = til(j,i)
          b = til(jp1,i)
          c = til(jp2,i)

          if ( til(1,k) == b ) then
            d = til(3,k)
          else if ( til(2,k) == b ) then
            d = til(1,k)
          else
            d = til(2,k)
          end if

          if ( diaedg(vcl(1,c),vcl(2,c),vcl(1,a),vcl(2,a),vcl(1,d), &
               vcl(2,d),vcl(1,b),vcl(2,b)) == 1 ) then
            nlo = nlo + 1
          end if

        end if

      end do

    end do

    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) '  NLO =  ', nlo

    call delaunay_print ( npt, vcl, ntri, til, tnbr )

  end do

  return
end
subroutine test08 ( filename )

!*****************************************************************************80
!
!! TEST08 tests EQDIS2, INTPG, MFDEC2, MMASEP, SEPMDF, SEPSHP, SFDWMF, SFUPMF, TRISIZ.
!
  implicit none

  integer ( kind = 4 ), parameter :: incr = 10000
  integer ( kind = 4 ), parameter :: inunit = 1
  integer ( kind = 4 ), parameter :: maxed = 101
  integer ( kind = 4 ), parameter :: maxhv = 350
  integer ( kind = 4 ), parameter :: maxiw = 900
  integer ( kind = 4 ), parameter :: maxnc = 30
  integer ( kind = 4 ), parameter :: maxnv = 400
  integer ( kind = 4 ), parameter :: maxpv = 2000
  integer ( kind = 4 ), parameter :: maxvc = 800
  integer ( kind = 4 ), parameter :: maxwk = 1500

  real ( kind = 8 ) angmin
  real ( kind = 8 ) angspc
  real ( kind = 8 ) angtol
  real ( kind = 8 ) area(maxhv)
  integer ( kind = 4 ) case
  real ( kind = 8 ) degrees_to_radians
  real ( kind = 8 ) dmin
  integer ( kind = 4 ) edge(4,maxed)
  character ( len = * ) filename
  real ( kind = 8 ) h(maxhv)
  logical hflag
  integer ( kind = 4 ) ht(0:maxed-1)
  integer ( kind = 4 ) htsiz
  integer ( kind = 4 ) hvl(maxhv)
  integer ( kind = 4 ) i
  real ( kind = 8 ) iang(maxpv)
  integer ( kind = 4 ) icur(maxnc)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ivrt(maxnv)
  integer ( kind = 4 ) iwk(maxiw)
  integer ( kind = 4 ) j
  real ( kind = 8 ) kappa
  integer ( kind = 4 ) map(maxnc)
  integer ( kind = 4 ) msglvl
  integer ( kind = 4 ) ncur
  integer ( kind = 4 ) nh
  integer ( kind = 4 ) nhola
  integer ( kind = 4 ) nhole
  integer ( kind = 4 ) nmin
  integer ( kind = 4 ) npolg
  integer ( kind = 4 ) nrfv
  integer ( kind = 4 ) nsc
  integer ( kind = 4 ) ntri(maxhv)
  integer ( kind = 4 ) ntrid
  integer ( kind = 4 ) ntrie
  integer ( kind = 4 ) nv
  integer ( kind = 4 ) nvbc(maxnc)
  integer ( kind = 4 ) nvc
  integer ( kind = 4 ) nvcin
  integer ( kind = 4 ) nvert
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ) prime
  real ( kind = 8 ) psi(maxhv)
  integer ( kind = 4 ) pvl(4,maxpv)
  real ( kind = 8 ) radians_to_degrees
  integer ( kind = 4 ) regnum(maxhv)
  character ( len = 20 ) rgname
  real ( kind = 8 ) tol
  real ( kind = 8 ) tolin
  real ( kind = 8 ) umdf2
  real ( kind = 8 ) vcl(2,maxvc)
  real ( kind = 8 ) wk(maxwk)
!
  external umdf2
!
  tol = 100.0D+00 * epsilon ( tol )
!
!  Read in vertices of general polygonal region.
!  CASE = 1 : simple polygon or multiply connected polygonal region
!  CASE = 2 : general polygonal region with holes and interfaces
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08'
  write ( *, '(a)' ) '  Reading input file:'
  write ( *, '(a)' ) filename

  open ( unit = inunit, file = filename, form = 'formatted' )

  read ( inunit, '(a)' ) rgname
  read ( inunit, * ) tolin,angspc,angtol,kappa,dmin,nmin,ntrid

  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) '  Region: ', rgname
  write ( *, 700) tol,angspc,angtol,kappa,dmin,nmin,ntrid

  angspc = degrees_to_radians ( angspc )
  angtol = degrees_to_radians ( angtol )

  hflag = ( 0.0D+00 <= kappa .and. kappa <= 1.0D+00 )
  read ( inunit, * ) case,nvc,ncur,msglvl

  if ( maxvc < nvc ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST08 - Error'
    write ( *, '(a)' ) '  MAXVC < NVC.'
    return
  end if

  if ( maxnc < ncur ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST08 - Error'
    write ( *, '(a)' ) '  MAXNC < NCUR.'
    return
  end if

  read ( inunit, * ) nvbc(1:ncur)

  if ( case == 2 ) then
    read ( inunit, * ) icur(1:ncur)
  end if

  read ( inunit, * ) (vcl(1,i),vcl(2,i),i=1,nvc)
!
!  Call DSMCPR or DSPGDC to set data structures in arrays
!  REGNUM, HVL, PVL, IANG, HOLV = IWK.
!
  if ( case == 1 ) then

    nhole = ncur - 1

    call dsmcpr ( nhole, nvbc, vcl, maxhv, maxpv, maxiw, nvc, npolg, nvert, &
      nhola, regnum, hvl, pvl, iang, iwk, ierror )

  else if ( case == 2 ) then

    nv = sum ( nvbc(1:ncur) )

    if ( maxnv < nv ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST08 - Error'
      write ( *, '(a)' ) '  MAXNV < NV.'
      return
    end if

    read ( inunit, * ) ivrt(1:nv)

    nsc = 0
    do i = 1, nv
      if ( ivrt(i) < 0 ) then
        nsc = nsc + 1
      end if
    end do
    nsc = nsc / 2

    if ( maxed < nsc ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST08 - Fatal error'
      write ( *, '(a)' ) '  MAXED < NSC.'
      return
    end if

    htsiz = min ( prime ( nsc / 2 ), maxed )

    call dspgdc ( nvc, vcl, incr, ncur, nvbc, icur, ivrt, maxhv, maxpv, maxiw, &
      npolg, nvert, nhole, nhola, regnum, hvl, pvl, iang, iwk, htsiz, nsc, &
      ht, edge, map, ierror )

  end if

  close ( unit = inunit )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST08 - Error!'
    write ( *, '(a,i6)' ) '  IERROR = ', ierror
    return
  end if
!
!  Print out data structures and measurements after initialization.
!
  nvcin = nvc
  nh = nhole*2 + nhola
  write ( *,670) msglvl
  write ( *,630) nvc,(i,vcl(1,i),vcl(2,i),i=1,nvc)
  write ( *,640) npolg,nhole,(hvl(i),i=1,npolg+nhole)
  write ( *,650) (regnum(i),i=1,npolg)
  write ( *,660) nvert,(i,(pvl(j,i),j=1,4),iang(i),i=1,nvert)
  write ( *,640) nhola,nh,(iwk(i),i=1,nh)

  nrfv = 0
  do i = 1, nvert
    if ( pi + tol < iang(i) ) then
      nrfv = nrfv + 1
    end if
  end do

  angmin = minval ( iang(1:nvert) )

  angmin = radians_to_degrees ( angmin )

  write (*,710) nvc,npolg,nvert,nhole,nhola,nrfv,angmin
!
!  Obtain simple and convex polygon decompositions, and print measurements.
!
  if ( msglvl == 2) then
    write ( *,670)
  end if

  call spdec2 ( angspc,angtol,nvc,npolg,nvert,nhole,nhola,maxvc,maxhv, &
    maxpv,maxiw-nh,maxwk,iwk,vcl,regnum,hvl,pvl,iang,iwk(nh+1),wk,ierror)

  call cvdec2 ( angspc, angtol, nvc, npolg, nvert, maxvc, maxhv, maxpv, maxiw, &
    maxwk, vcl, regnum, hvl, pvl, iang, iwk, wk, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST08 - Error!'
    write ( *, '(a,i6)' ) '  IERROR = ', ierror
    return
  end if

  angmin = minval ( iang(1:nvert) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08:'
  write ( *, '(a)' ) '  Before call to EQDIS2:'
  write ( *, '(a,i6)' ) '  NVC =    ', nvc
  write ( *, '(a,i6)' ) '  NPOLG =  ', npolg
  write ( *, '(a,i6)' ) '  NVERT =  ', nvert
  write ( *, '(a,g14.6)' ) '  ANGMIN = ', angmin
!
!  Obtain further convex polygon decomposition based on mesh
!  distribution function, and triangle sizes for the polygons.
!
  write ( *, '(a)' ) 'DEBUG: Call EQDIS2'

  call eqdis2 ( hflag, umdf2, kappa, angspc, angtol, dmin, nmin, ntrid, nvc, &
    npolg, nvert, maxvc, maxhv, maxpv, maxiw, maxwk, vcl, regnum, hvl, pvl, &
    iang, area, psi, h, iwk, wk, ierror )

  if ( msglvl == 2 ) then
    write ( *,670) 0, 0, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00
  end if

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST08 - Error!'
    write ( *, '(a,i6)' ) '  IERROR = ', ierror
    return
  end if
!
!  Print out data structures and measurements after decomposition.
!
  write ( *,630) nvc,(i,vcl(1,i),vcl(2,i),i=nvcin+1,nvc)
  write ( *,680) npolg, hvl(1:npolg)
  write ( *,650) regnum(1:npolg)
  write ( *,660) nvert,(i,(pvl(j,i),j=1,4),iang(i),i=1,nvert)

  angmin = minval ( iang(1:nvert) )
  angmin = radians_to_degrees ( angmin )

  ntrie = 0
  do i = 1, npolg
    ntri(i) = int ( ntrid * psi(i) * area(i) + 0.5D+00 )
    ntrie = ntrie + ntri(i)
  end do

  write ( *,690) (i,area(i),psi(i),h(i),ntri(i),i=1,npolg)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08:'
  write ( *, '(a)' ) '  After call to EQDIS2:'
  write ( *, '(a,i6)' ) '  NVC =    ', nvc
  write ( *, '(a,i6)' ) '  NPOLG =  ', npolg
  write ( *, '(a,i6)' ) '  NVERT =  ', nvert
  write ( *, '(a,i6)' ) '  NTRIE =  ', ntrie
  write ( *, '(a,g14.6)' ) '  ANGMIN = ', angmin

  630 format (/1x,i7/(1x,i7,2f15.7))
  640 format (/1x,2i7/(1x,10i7))
  650 format (/(1x,10i7))
  660 format (/1x,i7/(1x,5i7,f15.7))
  670 format (1x,2i7,4f15.7)
  680 format (/1x,i7/(1x,10i7))
  690 format (/(1x,i7,3e15.7,i7))
  700 format ( 'input : tol=',d15.7,'   angspc=',f9.3, &
    '   angtol=',f9.3/9x,'kappa=',f9.3,'   dmin=',f9.3, &
    '   nmin=',i5,'   ntrid=',i7)
  710 format (1x,'initds: nvc=',i7,'   npolg=',i7,'   nvert=',i7, &
    '   nhole=',i7/9x,'nhola=',i7,'   nrfv=',i7,'   angmin=',f9.3 )

  return
end
subroutine test09 ( )

!*****************************************************************************80
!
!! TEST09 tests LRLINE.
!
  implicit none

  real ( kind = 8 ) dv
  integer ( kind = 4 ) lr
  integer ( kind = 4 ) lrline
  real ( kind = 8 ) xu
  real ( kind = 8 ) xv1
  real ( kind = 8 ) xv2
  real ( kind = 8 ) yu
  real ( kind = 8 ) yv1
  real ( kind = 8 ) yv2

  xu = 0.0D+00
  yu = 0.0D+00

  xv1 =  0.0D+00
  yv1 = -1.0D+00
  xv2 =  1.0D+00
  yv2 =  0.0D+00

  dv = 0.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST09'
  write ( *, '(a)' ) '  LRLINE determines if a point is to the right,'
  write ( *, '(a)' ) '  left, or on a directed line that is a directed'
  write ( *, '(a)' ) '  distance away from a directed line from one'
  write ( *, '(a)' ) '  point to another.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The directed base line goes from'
  write ( *, '(2g14.6)' ) xv1, yv1
  write ( *, '(a)' ) '  to'
  write ( *, '(2g14.6)' ) xv2, yv2
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  The directed line distance is ', dv
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The point to be located is'
  write ( *, '(2g14.6)' ) xu, yu

  lr = lrline ( xu, yu, xv1, yv1, xv2, yv2, dv )

  write ( *, '(a)' ) ' '
  if ( lr == -1 ) then
    write ( *, '(a)' ) '  The point is to the right of the line.'
  else if ( lr == 0 ) then
    write ( *, '(a)' ) '  The point is on the line.'
  else if ( lr == 1 ) then
    write ( *, '(a)' ) '  The point is to the left of the line.'
  end if

  return
end
subroutine test10 ( )

!*****************************************************************************80
!
!! TEST10 tests LUFAC,  LUSOL.
!
  implicit none

  integer ( kind = 4 ), parameter :: nmax = 20

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
  real ( kind = 8 ) tol
  real ( kind = 8 ) urand

  tol = 100.0D+00 * epsilon ( tol )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST10'
  write ( *, '(a)' ) '  LUFAC factors a linear system;'
  write ( *, '(a)' ) '  LUSOL solves a factored linear system.'

  n = 4
  seed = 1952

  b(1:n) = 0.0D+00

  do j = 1, n
    do i = 1, n
      a(i,j) = urand(seed) * 2.0D+00 - 1.0D+00
      b(i) = b(i) + a(i,j)
    end do
  end do

  if ( n <= 4 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'a, b'
    do i = 1, n
      write ( *, '(5f15.7)' ) (a(i,j),j=1,n),b(i)
    end do
  end if

  call lufac ( a, nmax, n, tol, ipvt, singlr )

  if ( singlr ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The matrix is singular'
    return
  end if

  call lusol ( a, nmax, n, ipvt, b )

  if ( n <= 4 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ipvt, lu'
    ipvt(n) = n

    do i = 1, n
      write ( *, '(i5,4f15.7)' ) ipvt(i), (a(i,j),j=1,n)
    end do

  end if

  emax = 0.0D+00
  esum = 0.0D+00
  do i = 1, n
    t = abs ( b(i) - 1.0D+00 )
    emax = max ( emax, t )
    esum = esum + t
  end do

  write ( *,630) (b(i),i=1,min(4,n))
  write ( *,640) emax,esum

  630 format ('  x = ',4f15.7)
  640 format ('  emax,esum = ',2e15.7)

  return
end
subroutine test11 ( filename )

!*****************************************************************************80
!
!! TEST11 tests DSMDF2, MDF2, PRMDF2.
!
  implicit none

  integer ( kind = 4 ), parameter :: incr = 10000
  integer ( kind = 4 ), parameter :: inunit = 1
  integer ( kind = 4 ), parameter :: maxed = 101
  integer ( kind = 4 ), parameter :: maxhv = 200
  integer ( kind = 4 ), parameter :: maxiw = 900
  integer ( kind = 4 ), parameter :: maxnc = 30
  integer ( kind = 4 ), parameter :: maxpv = 1000
  integer ( kind = 4 ), parameter :: maxvc = 500
  integer ( kind = 4 ), parameter :: maxwk = 1500

  real ( kind = 8 ) angspc
  real ( kind = 8 ) angtol
  real ( kind = 8 ) area(maxhv)
  integer ( kind = 4 ) case
  real ( kind = 8 ) degrees_to_radians
  integer ( kind = 4 ) edge(4,maxed)
  real ( kind = 8 ) edgval(maxpv)
  character ( len = * ) filename
  integer ( kind = 4 ) ht(0:maxed-1)
  integer ( kind = 4 ) htsiz
  integer ( kind = 4 ) hvl(maxhv)
  integer ( kind = 4 ) i
  real ( kind = 8 ) iang(maxpv)
  integer ( kind = 4 ) icur(maxnc)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ifv
  integer ( kind = 4 ) ivrt(maxpv)
  integer ( kind = 4 ) iwk(maxiw)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) map(maxnc)
  real ( kind = 8 ) mdf2
  integer ( kind = 4 ) ncur
  integer ( kind = 4 ) nev(maxhv)
  integer ( kind = 4 ) nh
  integer ( kind = 4 ) nhola
  integer ( kind = 4 ) nhole
  integer ( kind = 4 ) npolg
  integer ( kind = 4 ) nsc
  integer ( kind = 4 ) nv
  integer ( kind = 4 ) nvbc(maxnc)
  integer ( kind = 4 ) nvc
  integer ( kind = 4 ) nvert
  integer ( kind = 4 ) prime
  integer ( kind = 4 ) pvl(4,maxpv)
  integer ( kind = 4 ) regnum(maxhv)
  character ( len = 20 ) rgname
  real ( kind = 8 ) tol
  real ( kind = 8 ) tolin
  real ( kind = 8 ) val(maxhv)
  real ( kind = 8 ) vcl(2,maxvc)
  real ( kind = 8 ) vrtval(maxvc)
  real ( kind = 8 ) widsq(maxhv)
  real ( kind = 8 ) wk(maxwk)
  real ( kind = 8 ) x
  integer ( kind = 4 ) xivrt(maxhv+1)
  real ( kind = 8 ) y

  tol = 100.0D+00 * epsilon ( tol )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST11'
  write ( *, '(a)' ) '  Reading input file:'
  write ( *, '(a)' ) filename

  open ( unit = inunit, file = filename, form = 'formatted' )
!
!  Read in vertices of general polygonal region.
!  CASE = 1 : simple polygon or multiply connected polygonal region
!  CASE = 2 : general polygonal region with holes and interfaces
!
  read ( inunit, '(a)' ) rgname
  read ( inunit, * ) tolin,angspc,angtol

  angspc = degrees_to_radians ( angspc )
  angtol = degrees_to_radians ( angtol )

  read ( inunit, * ) case,nvc,ncur

  if ( maxvc < nvc ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST11 - Error'
    write ( *, '(a)' ) '  MAXVC < NVC.'
    return
  end if

  if ( maxnc < ncur ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST11 - Error'
    write ( *, '(a)' ) '  MAXNC < NCUR.'
    return
  end if

  read ( inunit, * ) (nvbc(i),i=1,ncur)
  if ( case == 2 ) then
    read ( inunit, * ) (icur(i),i=1,ncur)
  end if
  read ( inunit, * ) (vcl(1,i),vcl(2,i),i=1,nvc)
!
!  Call DSMCPR or DSPGDC to set data structures in arrays
!  REGNUM, HVL, PVL, IANG, HOLV = IWK.
!
  if ( case == 1 ) then

    nhole = ncur - 1

    call dsmcpr(nhole,nvbc,vcl,maxhv,maxpv,maxiw,nvc,npolg,nvert, &
      nhola,regnum,hvl,pvl,iang,iwk, ierror )

  else if ( case == 2 ) then

    nv = 0
    do i = 1, ncur
      nv = nv + nvbc(i)
    end do

    if ( maxpv < nv ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST11 - Error'
      write ( *, '(a)' ) '  MAXPC < NV.'
      return
    end if

    read ( inunit, * ) ivrt(1:nv)

    nsc = 0
    do i = 1, nv
      if ( ivrt(i) < 0 ) then
        nsc = nsc + 1
      end if
    end do
    nsc = nsc / 2

    if ( maxed < nsc ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST11 - Error'
      write ( *, '(a)' ) '  MAXED < NSC.'
      return
    end if

    htsiz = min ( prime(nsc/2), maxed )

    call dspgdc(nvc,vcl,incr,ncur,nvbc,icur,ivrt,maxhv,maxpv,maxiw, &
      npolg,nvert,nhole,nhola,regnum,hvl,pvl,iang,iwk,htsiz,nsc, &
      ht, edge, map, ierror )

  end if

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST11 - Error!'
    write ( *, '(a,i6)' ) '  IERROR = ', ierror
    return
  end if
!
!  Obtain simple and convex polygon decompositions.
!
  nh = nhole*2 + nhola

  call spdec2(angspc,angtol,nvc,npolg,nvert,nhole,nhola,maxvc,maxhv, &
    maxpv,maxiw-nh,maxwk,iwk,vcl,regnum,hvl,pvl,iang,iwk(nh+1),wk, ierror )

  call cvdec2 ( angspc, angtol, nvc, npolg, nvert, maxvc, maxhv, maxpv, &
    maxiw, maxwk, vcl, regnum, hvl, pvl, iang, iwk, wk, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST11 - Error!'
    write ( *, '(a,i6)' ) '  IERROR = ', ierror
    return
  end if
!
!  Initialize data structure for heuristic mesh distribution function
!  and evaluate function at centroid of each convex polygon.
!
  call dsmdf2 ( .true., nvc, npolg, maxwk, vcl, hvl, pvl, iang, &
    ivrt, xivrt, widsq, edgval, vrtval, area, wk, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST11 - Error!'
    write ( *, '(a,i6)' ) '  IERROR = ', ierror
    return
  end if

  do i = 1, npolg

    call prmdf2(i,widsq(i),ivrt,xivrt,edgval,vrtval,nev(i),ifv,iwk)

    if ( nev(i) == 0 ) then

      val(i) = widsq(i)

    else

      nv = xivrt(i+1) - xivrt(i)
      x = 0.0D+00
      y = 0.0D+00
      do j = xivrt(i), xivrt(i+1)-1
        x = x + vcl(1,ivrt(j))
        y = y + vcl(2,ivrt(j))
      end do

      val(i) = 1.0D+00 / mdf2(x/nv,y/nv,widsq(i),nev(i),ifv,iwk,ivrt, &
        edgval,vrtval,vcl)

    end if

  end do

  close ( unit = inunit )
!
!  Print arrays from calls to 3 mdf routines.
!
  write ( *,630) npolg,(i,xivrt(i),area(i),widsq(i),val(i),nev(i),i=1,npolg)
  write ( *,640) nvert,nvc,(i,ivrt(i),edgval(i),vrtval(i),i=1,nvc)
  write ( *,650) (i,ivrt(i),edgval(i),i=nvc+1,nvert)

  630 format (/1x,i7/(1x,2i7,3f15.7,i7))
  640 format (/1x,2i7/(1x,2i7,2f15.7))
  650 format (1x,2i7,f15.7)

  return
end
subroutine test12 ( )

!***********************************************************************
!
!! TEST12 tests PRIME.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) prime

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST12'
  write ( *, '(a)' ) '  PRIME finds the smallest prime bigger than'
  write ( *, '(a)' ) '  a given value.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I, PRIME(I)'
  write ( *, '(a)' ) ' '

  do i = 100, 500, 100
    write ( *, '(2i10)' ) i, prime ( i )
  end do

  return
end
subroutine test13 ( filename )

!*****************************************************************************80
!
!! TEST13 tests PTPOLG.
!
  implicit none

  integer ( kind = 4 ), parameter :: inunit = 1
  integer ( kind = 4 ), parameter :: maxn = 100

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  integer ( kind = 4 ) dim
  real ( kind = 8 ) dtol
  character ( len = * ) filename
  integer ( kind = 4 ) i
  integer ( kind = 4 ) inout
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  real ( kind = 8 ) nrml(3)
  integer ( kind = 4 ) pgind(0:maxn)
  real ( kind = 8 ) pt(3)
  real ( kind = 8 ) tol
  real ( kind = 8 ) vcl(3,maxn)

  tol = 100.0D+00 * epsilon ( tol )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST13'
  write ( *, '(a)' ) '  Reading input file:'
  write ( *, '(a)' ) filename

  open ( unit = inunit, file = filename, form = 'formatted' )

  dtol = 10.0D+00 * tol
  read ( inunit, * ) n,dim,a,b
  if ( n < 3 ) then
    return
  end if

  if ( maxn < n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DRPTPG - Error'
    write ( *, '(a)' ) '  MAXN < N.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) 'DIM = ', dim
  write ( *, '(a,g14.6)' ) '  A = ', a
  write ( *, '(a,g14.6)' ) '  B = ', b
  write ( *, '(a)' ) ' '

  read ( inunit, * ) (vcl(1,i),vcl(2,i),i=1,n)

  pgind(0) = n
  do i = 1, n
    pgind(i) = i
    if ( dim == 3 ) then
      vcl(3,i) = a * vcl(1,i) + b * vcl(2,i)
    end if
    write ( *, '(i5,3f15.7)' ) i,(vcl(j,i),j=1,dim)
  end do

  if ( dim == 3 ) then
    c = sqrt ( 1.0D+00 + a**2 + b**2)
    nrml(1) = - a / c
    nrml(2) = - b / c
    nrml(3) = 1.0D+00 / c
  end if

  write ( *, '(a)' ) ' '

   20 continue

  read ( inunit, *, end=30 ) pt(1),pt(2)

  if ( dim == 3 ) then
    pt(3) = a * pt(1) + b * pt(2)
  end if

  call ptpolg ( dim, 3, n, 1, pgind, vcl, pt, nrml, dtol, inout )

  write ( *,630) inout,(pt(i),i=1,dim)
  go to 20

   30 continue

  close ( unit = inunit )

  630 format (1x,'inout=',i3,3x,'pt=',3f15.7)

  return
end
subroutine test14 ( )

!*****************************************************************************80
!
!! TEST14 tests ROTIAR.
!
  implicit none

  integer ( kind = 4 ), parameter :: maxn = 50

  integer ( kind = 4 ) a(maxn)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) shift

  n = 10
  shift = 3

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST14'
  write ( *, '(a)' ) '  ROTIAR "rotates" an array.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Using N = ', n
  write ( *, '(a,i6)' ) '  SHIFT = ', shift

  do i = 1, n
    a(i) = i
  end do

  write ( *, '(10i6)' ) a(1:n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Shifted array:'
  write ( *, '(a)' ) ' '

  call rotiar ( n, a, shift )

  write ( *, '(10i6)' ) a(1:n)

  return
end
subroutine test15 ( filename )

!*****************************************************************************80
!
!! TEST15 tests DIAM2, SHRNK2, WIDTH2.
!
  implicit none

  integer ( kind = 4 ), parameter :: inunit = 1
  integer ( kind = 4 ), parameter :: maxn = 100

  real ( kind = 8 ) diamsq
  character ( len = * ) filename
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) iedge(0:maxn)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nshr
  real ( kind = 8 ) sdist(0:maxn-1)
  real ( kind = 8 ) widsq
  real ( kind = 8 ) xc(0:maxn)
  real ( kind = 8 ) xs(0:maxn)
  real ( kind = 8 ) yc(0:maxn)
  real ( kind = 8 ) ys(0:maxn)
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST15'
  write ( *, '(a)' ) '  Reading data file: "' // trim ( filename ) // '".'

  open ( unit = inunit, file = filename, form = 'formatted' )

  read ( inunit, * ) n

  if ( maxn < n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST15 - Error!'
    write ( *, '(a,i6)' ) '  N = ', n
    write ( *, '(a,i6)' ) '  Maximum legal N is ', maxn
    return
  end if

  read ( inunit, * ) (xc(i),yc(i),sdist(i),i=0,n-1)

  close ( unit = inunit )

  xc(n) = xc(0)
  yc(n) = yc(0)
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  N = ', n
  write ( *, '(2f15.7)' ) (xc(i),yc(i),i=0,n)
  xs(0) = 0.0D+00
  ys(0) = 0.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Calling SHRNK2'

  call shrnk2 ( n, xc, yc, sdist, nshr, xs, ys, iedge, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST15 - Error!'
    write ( *, '(a,i6)' ) '  SHRNK2 returns IERROR = ', ierror
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  NSHR = ', n
  write ( *, '(2f15.7)' ) (xs(i),ys(i),i=0,nshr)

  call diam2 ( n, xc(1), yc(1), i1, i2, diamsq, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST15 - Error!'
    write ( *, '(a,i6)' ) '  DIAM2 returns IERROR = ', ierror
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(2i5,f18.7)' ) i1,i2,diamsq

  call width2 ( n, xc(1), yc(1), i1, i2, widsq, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST15 - Error!'
    write ( *, '(a,i6)' ) '  WIDTH2 returns IERROR = ', ierror
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(2i5,f18.7)' ) i1,i2,widsq

  if ( nshr < 3 ) then
    return
  end if

  call diam2 ( nshr, xs(1), ys(1), i1, i2, diamsq, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST15 - Error!'
    write ( *, '(a,i6)' ) '  DIAM2 returns IERROR = ', ierror
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(2i5,f18.7)' ) i1,i2,diamsq

  call width2 ( nshr, xs(1), ys(1), i1, i2, widsq, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST15 - Error!'
    write ( *, '(a,i6)' ) '  WIDTH2 returns IERROR = ', ierror
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(2i5,f18.7)' ) i1,i2,widsq

  return
end
subroutine test16 ( )

!*****************************************************************************80
!
!! TEST16 tests DHPSRT, IHPSRT, RANDPT.
!
  implicit none

  integer ( kind = 4 ), parameter :: maxk = 4
  integer ( kind = 4 ), parameter :: maxn = 100

  integer ( kind = 4 ) axis
  real ( kind = 8 ) da(maxk,maxn)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(maxk,maxn)
  logical iflag
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) map(maxn)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nptav
  real ( kind = 8 ), dimension ( maxk ) :: scale = (/ &
    1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00 /)
  integer ( kind = 4 ) seed
  real ( kind = 8 ), dimension ( maxk ) :: trans = (/ &
    0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST16'

  k = 2
  n = 10
  seed = 1952
  axis = 0
  nptav = 0

  iflag = (k < 0)
  k = abs(k)

  if ( k < 1 .or. maxk < k ) then
    return
  end if

  if ( n < 1 .or. maxn < n ) then
    return
  end if

  do j = 1, n
    map(j) = j
  end do

  call randpt ( k, n, seed, axis, nptav, scale, trans, maxk, da )

  if ( iflag ) then

    do j = 1, n
      do i = 1, k
        ia(i,j) = int ( n * da(i,j) )
      end do
    end do

    do j = 1, n
      write ( *, '(5i5)' ) j,(ia(i,j),i=1,k)
    end do

    call ihpsrt ( k, n, maxk, ia, map )

    write ( *, '(a)' ) ' '
    do j = 1, n
      write ( *, '(5i5)' ) map(j),(ia(i,map(j)),i=1,k)
    end do

  else

    do j = 1, n
      write ( *, '(i5,4f15.7)' ) j,(da(i,j),i=1,k)
    end do

    call dhpsrt ( k, n, maxk, da, map )

    write ( *,'(a)' ) ' '
    do j = 1, n
     write ( *, '(i5,4f15.7)' ) map(j), (da(i,map(j)),i=1,k)
    end do

  end if

  return
end
subroutine test17 ( filename )

!*****************************************************************************80
!
!! TEST17 tests BEDGMV, CVDTRI, FNDTRI, INTTRI, LOP, MTREDG, ROTPG, TMERGE, TRINBR, TRIPR2, TRPOLG.
!
  implicit none

  integer ( kind = 4 ), parameter :: incr = 10000
  integer ( kind = 4 ), parameter :: inunit = 1
  integer ( kind = 4 ), parameter :: maxed = 1201
  integer ( kind = 4 ), parameter :: maxhv = 350
  integer ( kind = 4 ), parameter :: maxiw = 900
  integer ( kind = 4 ), parameter :: maxnc = 30
  integer ( kind = 4 ), parameter :: maxnv = 400
  integer ( kind = 4 ), parameter :: maxpv = 2000
  integer ( kind = 4 ), parameter :: maxti = 8000
  integer ( kind = 4 ), parameter :: maxvc = 5000
  integer ( kind = 4 ), parameter :: maxwk = 1500
  integer ( kind = 4 ), parameter :: nfreq = 12

  integer ( kind = 4 ) a
  integer ( kind = 4 ) afreq(0:nfreq-1)
  real ( kind = 8 ) anga
  real ( kind = 8 ) angb
  real ( kind = 8 ) angc
  real ( kind = 8 ) angle
  real ( kind = 8 ) angmax
  real ( kind = 8 ) angmin
  real ( kind = 8 ) angspc
  real ( kind = 8 ) angtol
  real ( kind = 8 ) area(maxhv)
  integer ( kind = 4 ) b
  integer ( kind = 4 ) c
  integer ( kind = 4 ) case
  real ( kind = 8 ) degrees_to_radians
  real ( kind = 8 ) delta
  real ( kind = 8 ) dmin
  integer ( kind = 4 ) edge(4,maxed)
  character ( len = * ) filename
  real ( kind = 8 ) h(maxhv)
  logical hflag
  integer ( kind = 4 ) ht(0:maxed-1)
  integer ( kind = 4 ) htsiz
  integer ( kind = 4 ) hvl(maxhv)
  integer ( kind = 4 ) i
  real ( kind = 8 ) iang(maxpv)
  integer ( kind = 4 ) icur(maxnc)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ivrt(maxnv)
  integer ( kind = 4 ) iwk(maxiw)
  integer ( kind = 4 ) j
  real ( kind = 8 ) kappa
  integer ( kind = 4 ) map(maxnc)
  integer ( kind = 4 ) msglvl
  integer ( kind = 4 ) ncur
  integer ( kind = 4 ) nh
  integer ( kind = 4 ) nhola
  integer ( kind = 4 ) nhole
  integer ( kind = 4 ) nmin
  integer ( kind = 4 ) npolg
  integer ( kind = 4 ) nrfv
  integer ( kind = 4 ) nsc
  integer ( kind = 4 ) ntri
  integer ( kind = 4 ) ntrid
  integer ( kind = 4 ) nv
  integer ( kind = 4 ) nvbc(maxnc)
  integer ( kind = 4 ) nvc
  integer ( kind = 4 ) nvcin
  integer ( kind = 4 ) nvert
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ) prime
  real ( kind = 8 ) psi(maxhv)
  integer ( kind = 4 ) pvl(4,maxpv)
  real ( kind = 8 ) radians_to_degrees
  integer ( kind = 4 ) regnum(maxhv)
  character ( len = 20 ) rgname
  integer ( kind = 4 ) til(3,maxti)
  integer ( kind = 4 ) tnbr(3,maxti)
  real ( kind = 8 ) tol
  real ( kind = 8 ) tolin
  integer ( kind = 4 ) tstart(maxhv)
  real ( kind = 8 ), external :: umdf2
  real ( kind = 8 ) vcl(2,maxvc)
  integer ( kind = 4 ) vnum(maxpv)
  integer ( kind = 4 ) vstart(maxpv)
  real ( kind = 8 ) wk(maxwk)
!
  tol = 100.0D+00 * epsilon ( tol )
!
!  Read in vertices of general polygonal region.
!  CASE = 1 : simple polygon or multiply connected polygonal region
!  CASE = 2 : general polygonal region with holes and interfaces
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST17'
  write ( *, '(a)' ) '  Reading input file: "' // trim ( filename ) // '".'

  open ( unit = inunit, file = filename, form = 'formatted' )

  read ( inunit, '(a)' ) rgname
  read ( inunit, * ) tolin, angspc, angtol, kappa, dmin, nmin, ntrid

  write ( *, '(a)' ) ' '
  write ( *, 710 ) rgname, tol, angspc, angtol, kappa, dmin, nmin, ntrid
  angspc = degrees_to_radians ( angspc )
  angtol = degrees_to_radians ( angtol )
  hflag = ( 0.0D+00 <= kappa .and. kappa <= 1.0D+00 )

  read ( inunit, * ) case, nvc, ncur, msglvl

  if ( maxvc < nvc ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST17 - Error'
    write ( *, '(a)' ) '  MAXVC < NVC.'
    return
  end if

  if ( maxnc < ncur ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST17 - Error'
    write ( *, '(a)' ) '  MAXNC < NCUR.'
    return
  end if

  read ( inunit, * ) nvbc(1:ncur)

  if ( case == 2 ) then
    read ( inunit, * ) icur(1:ncur)
  end if

  read ( inunit, * ) (vcl(1,i),vcl(2,i),i=1,nvc)
!
!  Call DSMCPR or DSPGDC to set the data structures in arrays
!  REGNUM, HVL, PVL, IANG, HOLV = IWK.
!
  if ( case == 1 ) then

    nhole = ncur - 1

    call dsmcpr ( nhole, nvbc, vcl, maxhv, maxpv, maxiw, nvc, npolg, nvert, &
      nhola, regnum, hvl, pvl, iang, iwk, ierror )

  else if ( case == 2 ) then

    nv = sum ( nvbc(1:ncur) )

    if ( maxnv < nv ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST17 - Error'
      write ( *, '(a)' ) '  MAXNV < NV.'
      return
    end if

    read ( inunit, * ) ivrt(1:nv)

    nsc = 0
    do i = 1, nv
      if ( ivrt(i) < 0 ) then
        nsc = nsc + 1
      end if
    end do

    nsc = nsc / 2

    if ( maxed < nsc ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST17 - Error'
      write ( *, '(a)' ) '  MAXED < NSC.'
      return
    end if

    htsiz = min ( prime(nsc/2), maxed )

    call dspgdc ( nvc, vcl, incr, ncur, nvbc, icur, ivrt, maxhv, maxpv, maxiw, &
      npolg, nvert, nhole, nhola, regnum, hvl, pvl, iang, iwk, htsiz, nsc, &
      ht, edge, map, ierror )

  end if

  close ( unit = inunit )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST17 - Error!'
    write ( *, '(a,i6)' ) '  IERROR = ', ierror
    return
  end if
!
!  Print out the data structures and measurements after initialization.
!
  nvcin = nvc
  nh = nhole * 2 + nhola
  write ( *, '(a,i6)' ) 'MSGLVL = ', msglvl
  write ( *,630) nvc,(i,vcl(1,i),vcl(2,i),i=1,nvc)
  write ( *,640) npolg,nhole,(hvl(i),i=1,npolg+nhole)
  write ( *,650) regnum(1:npolg)
  write ( *,660) nvert, (i,(pvl(j,i),j=1,4), iang(i), i = 1,nvert)
  write ( *,640) nhola, nh, iwk(1:nh)

  nrfv = 0
  do i = 1, nvert
    if ( pi + tol < iang(i) ) then
      nrfv = nrfv + 1
    end if
  end do

  angmin = minval ( iang(1:nvert) )
  angmin = radians_to_degrees ( angmin )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST17: INITDS:'
  write ( *, '(a,i6)' ) '  NVC =    ', nvc
  write ( *, '(a,i6)' ) '  NPOLG =  ', npolg
  write ( *, '(a,i6)' ) '  NVERT =  ', nvert
  write ( *, '(a,i6)' ) '  NHOLE =  ', nhole
  write ( *, '(a,i6)' ) '  NHOLA =  ', nhola
  write ( *, '(a,i6)' ) '  NRFV =   ', nrfv
  write ( *, '(a,g14.6)' ) '  ANGMIN = ', angmin
!
!  Obtain simple and convex polygon decompositions, and print measurements.
!
  if ( msglvl == 2 ) then
    write ( *,670)
  end if

  call spdec2 ( angspc, angtol, nvc, npolg, nvert, nhole, nhola, maxvc, maxhv, &
    maxpv, maxiw-nh, maxwk, iwk, vcl, regnum, hvl, pvl, iang, iwk(nh+1), wk, &
    ierror )

  call cvdec2 ( angspc, angtol, nvc, npolg, nvert, maxvc, maxhv, maxpv, maxiw, &
    maxwk, vcl, regnum, hvl, pvl, iang, iwk, wk, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST17 - Error!'
    write ( *, '(a,i6)' ) '  IERROR = ', ierror
    return
  end if

  angmin = minval ( iang(1:nvert) )

  angmin = radians_to_degrees ( angmin )

  write (*,730) nvc, npolg, nvert, angmin
!
!  Obtain further convex polygon decomposition based on mesh
!  distribution function, and triangle sizes for the polygons.
!  Then print measurements.
!
  call eqdis2 ( hflag, umdf2, kappa, angspc, angtol, dmin, nmin, ntrid, nvc, &
    npolg, nvert, maxvc, maxhv, maxpv, maxiw, maxwk, vcl, regnum, hvl, pvl, &
    iang, area, psi, h, iwk, wk, ierror )

  if ( msglvl == 2 ) then
    write ( *,670) 0, 0, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00
  end if

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST17 - Error!'
    write ( *, '(a,i6)' ) '  IERROR = ', ierror
    return
  end if

  write ( *,630) nvc,(i,vcl(1,i),vcl(2,i),i=nvcin+1,nvc)
  write ( *,680) npolg, hvl(1:npolg)
  write ( *,650) regnum(1:npolg)

  angmin = radians_to_degrees ( minval ( iang(1:nvert) ) )

  write (*,740) nvc, npolg, nvert, angmin
  htsiz = min ( prime(nvcin*2), maxed )
  nvcin = nvc
!
!  Triangulate each convex polygon in the decomposition, according to
!  the mesh spacings in the H array.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Calling TRIPR2'

  call tripr2 ( nvc, npolg, nvert, maxvc, maxti, maxiw, maxwk, h, vcl, hvl, &
    pvl, iang, ntri, til, vstart, vnum, tstart, iwk, wk, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST17 - Error!'
    write ( *, '(a,i6)' ) '  IERROR = ', ierror
    return
  end if
!
!  Compute TNBR array. Then print arrays and measurements.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Calling TRINBR'

  call trinbr ( nvc, ntri, til, tnbr, htsiz, maxed, ht, edge, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST17 - Error!'
    write ( *, '(a,i6)' ) '  IERROR = ', ierror
    return
  end if

  write ( *,690) nvert,(i,(pvl(j,i),j=1,4),iang(i),vstart(i),vnum(i),i=1,nvert)
  write ( *,630) nvc,(i,vcl(1,i),vcl(2,i),i=nvcin+1,nvc)
  write ( *,650) tstart(1:npolg)
  write ( *,700) ntri,(i,(til(j,i),j=1,3),(tnbr(j,i),j=1,3),i=1,ntri)

  angmin = pi
  angmax = 0.0D+00
  delta = pi / real ( nfreq, kind = 8 )
  afreq(0:nfreq-1) = 0

  do i = 1, ntri

    a = til(1,i)
    b = til(2,i)
    c = til(3,i)

    anga = angle(vcl(1,c),vcl(2,c),vcl(1,a),vcl(2,a),vcl(1,b),vcl(2,b))

    angb = angle(vcl(1,a),vcl(2,a),vcl(1,b),vcl(2,b),vcl(1,c),vcl(2,c))

    angc = pi - anga - angb
    angmin = min ( angmin, anga, angb, angc )
    angmax = max ( angmax, anga, angb, angc )

    if ( abs ( anga - 0.5D+00 * pi ) <= tol ) then
      anga = 0.5D+00 * pi - 2.0D+00 * tol
    end if

    j = min ( nfreq - 1, int ( ( anga + tol ) / delta ) )
    afreq(j) = afreq(j) + 1

    if ( abs ( angb - 0.5D+00 * pi ) <= tol ) then
      angb = 0.5D+00 * pi - 2.0D+00 * tol
    end if

    j = min ( nfreq - 1, int ( ( angb + tol ) / delta ) )
    afreq(j) = afreq(j) + 1

    if ( abs ( angc - 0.5D+00 * pi ) <= tol ) then
      angc = 0.5D+00 * pi - 2.0D+00 * tol
    end if

    j = min ( nfreq - 1, int ( ( angc + tol ) / delta ) )
!
!  Having some odd cases where J comes out NEGATIVE.
!
    if ( j < 0 ) then
      j = 0
    end if

    afreq(j) = afreq(j) + 1

  end do

  angmin = radians_to_degrees ( angmin )
  angmax = radians_to_degrees ( angmax )
  write (*,750) nvc,ntri,angmin,angmax
  write (*,760) ( i*180.0D+00 / real ( nfreq, kind = 8 ), &
    afreq(i) / real ( 3*ntri, kind = 8 ),i=0,nfreq-1)

  630 format (/1x,i7/(1x,i7,2f15.7))
  640 format (/1x,2i7/(1x,10i7))
  650 format (/(1x,10i7))
  660 format (/1x,i7/(1x,5i7,f15.7))
  670 format (1x,2i7,4f15.7)
  680 format (/1x,i7/(1x,10i7))
  690 format (/1x,i7/(1x,5i7,f15.7,2i7))
  700 format (/1x,i7/(1x,7i7))
  710 format (1x,a20/1x,'input : tol=',d15.7,'   angspc=',f9.3, &
        '   angtol=',f9.3/9x,'kappa=',f9.3,'   dmin=',f9.3, &
        '   nmin=',i5,'   ntrid=',i7)
  730 format (1x,'decomp: nvc=',i7,'   npolg=',i7,'   nvert=',i7/ &
        9x,'angmin=',f9.3 )
  740 format (1x,'eqdist: nvc=',i7,'   npolg=',i7,'   nvert=',i7/ &
        9x,'angmin=',f9.3 )
  750 format (1x,'triang: nvc=',i7,'   ntri=',i7 / &
        9x,'angmin=',f9.3,'   angmax=',f9.3 )
  760 format (1x,'relative frequency of triangle angles'/4(f8.1,f10.5))

  return
end
subroutine test18 ( )

!*****************************************************************************80
!
!! TEST18 tests VISVRT, VORNBR.
!
  implicit none

  integer ( kind = 4 ), parameter :: maxn = 200

  real ( kind = 8 ) angle
  real ( kind = 8 ) angspc
  real ( kind = 8 ) degrees_to_radians
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ivert
  integer ( kind = 4 ) ivis(0:maxn)
  integer ( kind = 4 ) ivor(0:maxn)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) maxnv
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nvis
  integer ( kind = 4 ) nvor
  integer ( kind = 4 ) nvrt
  integer ( kind = 4 ) nvsvrt
  real ( kind = 8 ) phi
  real ( kind = 8 ) theta(0:maxn)
  real ( kind = 8 ) x(maxn)
  real ( kind = 8 ) xc(0:maxn)
  real ( kind = 8 ) xeye
  real ( kind = 8 ) xvor(0:maxn)
  real ( kind = 8 ) y(maxn)
  real ( kind = 8 ) yc(0:maxn)
  real ( kind = 8 ) yeye
  real ( kind = 8 ) yvor(0:maxn)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST18'
  write ( *, '(a)' ) '  Test being SKIPPED FOR NOW.'
  return

  read ( *, * ) n

  if ( n < 3 ) then
    return
  end if

  if ( maxn < n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST18 - Error'
    write ( *, '(a)' ) '  MAXN < N.'
    return
  end if

  read ( *, * ) ( x(i), y(i), i = 1, n )
  read ( *, * ) ivert
  read ( *, * ) angspc
!
!  1 <= IVERT <= N is index of polygon vertex for eyepoint.
!  ANGSPC is angle spacing parameter in degrees.
!
  angspc = degrees_to_radians ( angspc )
  xeye = x(ivert)
  yeye = y(ivert)
  nvrt = n - 2

  j = -1
  do i = ivert+1, n
    j = j + 1
    xc(j) = x(i)
    yc(j) = y(i)
  end do

  do i = 1, ivert-1
    j = j + 1
    xc(j) = x(i)
    yc(j) = y(i)
  end do

  write ( *, '(a,2g14.6)' ) '(XEYE,YEYE)= ', xeye, yeye
  write ( *,630) nvrt,(xc(i),yc(i),i=0,nvrt)

  call vispol ( xeye, yeye, nvrt, xc, yc, nvis, ivis )

  write ( *,640) nvis,(xc(i),yc(i),ivis(i),i=0,nvis)

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST18 - Fatal error!'
    write ( *, '(a,i6)' ) '  IERROR = ', ierror
    return
  end if

  phi = angle ( xc(nvis), yc(nvis), xeye, yeye, xc(0), yc(0) )
  maxnv = nvis + int ( phi / angspc )

  if ( maxn < maxnv ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST18 - Fatal error!'
    write ( *, '(a)' ) '  MAXN < MAXNV.'
    return
  end if

  call visvrt ( angspc, xeye, yeye, nvis, xc, yc, ivis, maxnv, nvsvrt, theta )

  write ( *,650) nvsvrt,(xc(i),yc(i),ivis(i),theta(i),i=0,nvsvrt)
  nvor = -1

  call vornbr ( xeye, yeye, nvsvrt, xc, yc, nvor, ivor, xvor, yvor, ierror )

  write ( *,640) nvor,(xvor(i),yvor(i),ivor(i),i=0,nvor)

  630 format (/1x,i5/(1x,2f15.7))
  640 format (/1x,i5/(1x,2f15.7,i5))
  650 format (/1x,i5/(1x,2f15.7,i5,f15.7))

  return
end
subroutine test19 ( )

!*****************************************************************************80
!
!! TEST19 tests ROTIPG, VISPOL.
!
  implicit none

  integer ( kind = 4 ), parameter :: maxn = 200

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ivert
  integer ( kind = 4 ) ivis(0:maxn+1)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nvis
  integer ( kind = 4 ) nvrt
  integer ( kind = 4 ) vptype
  real ( kind = 8 ) x(maxn)
  real ( kind = 8 ) xc(0:maxn+1)
  real ( kind = 8 ) xeye
  real ( kind = 8 ) y(maxn)
  real ( kind = 8 ) yc(0:maxn+1)
  real ( kind = 8 ) yeye

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST19'
  write ( *, '(a)' ) '  TEST BEING SKIPPED FOR NOW.'
  return

  read ( *, * ) n
  if ( n < 3 ) then
    return
  end if

  if ( maxn < n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST19 - Error'
    write ( *, '(a)' ) '  MAXN < N.'
    return
  end if

  read ( *,*) (x(i),y(i),i=1,n)
  read ( *,*) ivert,vptype,xeye,yeye
!
!  1 <= IVERT <= N is index of polygon vertex or IVERT <= 0 if
!  ROTIPG is to be called. BVTYPE = 0, 1, or 2 for boundary,
!  interior, or blocked exterior viewpoint. (XEYE,YEYE) is needed
!  for non-boundary viewpoints only, and must be visible from
!  (X(IVERT),Y(IVERT)) if 0 < IVERT.
!
  if ( vptype == 0 ) then

    xeye = x(ivert)
    yeye = y(ivert)
    nvrt = n - 2

    j = -1
    do i = ivert+1, n
      j = j + 1
      xc(j) = x(i)
      yc(j) = y(i)
    end do

    do i = 1, ivert-1
      j = j + 1
      xc(j) = x(i)
      yc(j) = y(i)
    end do

  else if ( vptype == 1 ) then

    nvrt = n

    if ( ivert <= 0 ) then

      do i = 1, n
        xc(i) = x(i)
        yc(i) = y(i)
      end do

      xc(0) = xc(n)
      yc(0) = yc(n)

      call rotipg ( xeye, yeye, nvrt, xc, yc )

    else

      j = -1

      do i = ivert, n
        j = j + 1
        xc(j) = x(i)
        yc(j) = y(i)
      end do

      do i = 1, ivert
        j = j + 1
        xc(j) = x(i)
        yc(j) = y(i)
      end do

    end if

  else

    nvrt = n

    if ( ivert <= 0 ) then

      do i = 1, n
        xc(n-i) = x(i)
        yc(n-i) = y(i)
      end do

      xc(n) = xc(0)
      yc(n) = yc(0)

      call rotipg ( xeye, yeye, nvrt, xc, yc )

    else

      j = -1

      do i = ivert, 1, -1
        j = j + 1
        xc(j) = x(i)
        yc(j) = y(i)
      end do

      do i = n, ivert, -1
        j = j + 1
        xc(j) = x(i)
        yc(j) = y(i)
      end do

    end if

  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a,2g14.6)' ) '(XEYE,YEYE)=', xeye,yeye
  write ( *,630) nvrt,(xc(i),yc(i),i=0,nvrt)

  call vispol ( xeye, yeye, nvrt, xc, yc, nvis, ivis )

  write ( *,640) nvis,(xc(i),yc(i),ivis(i),i=0,nvis)

  630 format (/1x,i5/(1x,2f15.7))
  640 format (/1x,i5/(1x,2f15.7,i5))

  return
end
subroutine pcdec ( )

!*****************************************************************************80
!
!! PCDEC generates PICTEX commands from DRDEC or DREQD output.
!
  implicit none

  integer ( kind = 4 ), parameter :: edgv = 4
  integer ( kind = 4 ), parameter :: loc = 1
  integer ( kind = 4 ), parameter :: maxho = 50
  integer ( kind = 4 ), parameter :: maxnc = 30
  integer ( kind = 4 ), parameter :: maxnv = 500
  integer ( kind = 4 ), parameter :: succ = 3

  integer ( kind = 4 ) holv(maxho)
  integer ( kind = 4 ) hvl(maxnc*2)
  integer ( kind = 4 ) i
  real ( kind = 8 ) iang(maxnv*2)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) l
  integer ( kind = 4 ) l1
  logical mid
  integer ( kind = 4 ) msglvl
  integer ( kind = 4 ) nh
  integer ( kind = 4 ) nhola
  integer ( kind = 4 ) nhole
  integer ( kind = 4 ) npolg
  integer ( kind = 4 ) nvc
  integer ( kind = 4 ) nvert
  integer ( kind = 4 ) pvl(4,maxnv*2)
  integer ( kind = 4 ) regnum(maxnc*2)
  integer ( kind = 4 ) t
  real ( kind = 8 ) vcl(2,maxnv)
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin

  read ( *,*) msglvl
  if ( msglvl /= 2) then
    return
  end if
  read ( *,*) nvc
  read ( *,*) (t,vcl(1,i),vcl(2,i),i=1,nvc)
  read ( *,*) npolg,nhole
  read ( *,*) (hvl(i),i=1,npolg+nhole)
  read ( *,*) (regnum(i),i=1,npolg)
  read ( *,*) nvert
  read ( *,*) (t,(pvl(j,i),j=1,4),iang(i),i=1,nvert)
  read ( *,*) nhola,nh
  read ( *,*) (holv(i),i=1,nh)

  xmin = vcl(1,1)
  ymin = vcl(2,1)
  xmax = xmin
  ymax = ymin
  do i = 2, nvc
    xmin = min(xmin,vcl(1,i))
    xmax = max(xmax,vcl(1,i))
    ymin = min(ymin,vcl(2,i))
    ymax = max(ymax,vcl(2,i))
  end do

  write ( *, '(a)' ) '\\begin{figure}[htbp]'
  write ( *, '(a)' ) '\\begin{center}'
  write ( *, '(a)' ) '\\mbox{\\beginpicture'
  write ( *, '(a)' ) '\\setcoordinatesystem units <0.2cm,0.2cm>'
  write ( *,610) xmin,xmax,ymin,ymax
  write ( *, '(a)' ) '\\setlinear'

  do i = 1, npolg+nhole

    mid = .false.
    j = hvl(i)
    l = pvl(loc,j)

   20   continue

    j1 = pvl(succ,j)
    l1 = pvl(loc,j1)

    if ( pvl(edgv,j) < j ) then

      if ( .not. mid ) then
        write ( *, '(a)' ) '\\plot'
        mid = .true.
        write ( *, '(2f15.7)' ) vcl(1,l),vcl(2,l)
      end if

      write ( *, '(2f15.7)' ) vcl(1,l1),vcl(2,l1)

    else if ( mid ) then

      write ( *, '(a)' ) '/'
      mid = .false.

    end if

    j = j1
    l = l1

    if ( j /= hvl(i) ) then
      go to 20
    end if

    if ( mid ) then
      write ( *, '(a)' ) '/'
    end if

  end do

   40 continue

   read ( *,*) i,j,x1,y1,x2,y2
   if ( i == 0 ) then
     go to 50
   end if
   write ( *,630) x1,y1,x2,y2
   go to 40

   50 continue

  write ( *, '(a)' ) '\\end{picture}'
  write ( *, '(a)' ) '\\end{center}'
  write ( *, '(a)' ) '\\end{figure}'

  610 format ('\\setplotarea x from ',f9.4,' to ',f9.4,', y from ',f9.4, &
    ' to ',f9.4)
  630 format ('\\plot ',4f15.7,' /')

  return
end
subroutine pcds ( )

!*****************************************************************************80
!
!! PCDS generates PICTEX commands from DRDS output.
!
  implicit none

  integer ( kind = 4 ), parameter :: edgv = 4
  integer ( kind = 4 ), parameter :: loc = 1
  integer ( kind = 4 ), parameter :: maxho = 50
  integer ( kind = 4 ), parameter :: maxnc = 30
  integer ( kind = 4 ), parameter :: maxnv = 500
  integer ( kind = 4 ), parameter :: succ = 3

  integer ( kind = 4 ) holv(maxho)
  integer ( kind = 4 ) hvl(maxnc*2)
  integer ( kind = 4 ) i
  real ( kind = 8 ) iang(maxnv*2)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) l
  integer ( kind = 4 ) l1
  logical mid
  integer ( kind = 4 ) nh
  integer ( kind = 4 ) nhola
  integer ( kind = 4 ) nhole
  integer ( kind = 4 ) npolg
  integer ( kind = 4 ) nvc
  integer ( kind = 4 ) nvert
  integer ( kind = 4 ) pvl(4,maxnv*2)
  integer ( kind = 4 ) regnum(maxnc*2)
  integer ( kind = 4 ) t
  real ( kind = 8 ) vcl(2,maxnv)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin
!
  read ( *,*) nvc
  read ( *,*) (t,vcl(1,i),vcl(2,i),i=1,nvc)
  read ( *,*) npolg,nhole
  read ( *,*) (hvl(i),i=1,npolg+nhole)
  read ( *,*) (regnum(i),i=1,npolg)
  read ( *,*) nvert
  read ( *,*) (t,(pvl(j,i),j=1,4),iang(i),i=1,nvert)
  read ( *,*) nhola,nh
  read ( *,*) (holv(i),i=1,nh)

  xmin = vcl(1,1)
  ymin = vcl(2,1)
  xmax = xmin
  ymax = ymin
  do i = 2, nvc
    xmin = min(xmin,vcl(1,i))
    xmax = max(xmax,vcl(1,i))
    ymin = min(ymin,vcl(2,i))
    ymax = max(ymax,vcl(2,i))
  end do

  write ( *, '(a)' ) '\\begin{figure}[htbp]'
  write ( *, '(a)' ) '\\begin{center}'
  write ( *, '(a)' ) '\\mbox{\\beginpicture'
  write ( *, '(a)' ) '\\setcoordinatesystem units <0.2cm,0.2cm>'
  write ( *, 610) xmin,xmax,ymin,ymax
  write ( *, '(a)' ) '\\setlinear'

  do i = 1, npolg+nhole

    mid = .false.
    j = hvl(i)
    l = pvl(loc,j)

   20   continue

    j1 = pvl(succ,j)
    l1 = pvl(loc,j1)

    if ( pvl(edgv,j) < j ) then

      if ( .not. mid ) then
        write ( *, '(a)' ) '\\plot'
        mid = .true.
        write ( *, '(2f15.7)' ) vcl(1,l),vcl(2,l)
      end if

      write ( *, '(2f15.7)' ) vcl(1,l1),vcl(2,l1)

    else if ( mid ) then

      write ( *, '(a)' ) '/'
      mid = .false.

    end if

    j = j1
    l = l1

    if ( j /= hvl(i) ) then
      go to 20
    end if

    if ( mid ) then
      write ( *, '(a)' ) '/'
    end if

  end do

  write ( *, '(a)' ) '\\endpicture}'
  write ( *, '(a)' ) '\\end{center}'
  write ( *, '(a)' ) '\\end{figure}'

  610 format ('\\setplotarea x from ',f9.4,' to ',f9.4,', y from ',f9.4, &
    ' to ',f9.4)

  return
end
subroutine pctri ( )

!*****************************************************************************80
!
!! PCTRI generates PICTEX commands from DRTRI output.
!
  implicit none

  integer ( kind = 4 ), parameter :: maxti = 8000
  integer ( kind = 4 ), parameter :: maxvc = 5000

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jp1
  integer ( kind = 4 ) msglvl
  integer ( kind = 4 ) nh
  integer ( kind = 4 ) nhola
  integer ( kind = 4 ) nhole
  integer ( kind = 4 ) npolg
  integer ( kind = 4 ) ntri
  integer ( kind = 4 ) nvc
  integer ( kind = 4 ) nvcin
  integer ( kind = 4 ) nvert
  integer ( kind = 4 ) t
  integer ( kind = 4 ) til(3,maxti)
  integer ( kind = 4 ) tnbr(3,maxti)
  real ( kind = 8 ) vcl(2,maxvc)
  real ( kind = 8 ) x
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin
!
  read ( *,*,end=40) xmin,xmax,ymin,ymax
  read ( *,*) msglvl
  if ( msglvl /= 0) then
    return
  end if
  read ( *,*) nvc
  read ( *,*) (t,vcl(1,i),vcl(2,i),i=1,nvc)
  read ( *,*) npolg,nhole
  read ( *,*) (t,i=1,npolg+nhole)
  read ( *,*) (t,i=1,npolg)
  read ( *,*) nvert
  read ( *,*) ((t,j=1,5),x,i=1,nvert)
  read ( *,*) nhola,nh
  read ( *,*) (t,i=1,nh)
  nvcin = nvc
  read ( *,*) nvc
  read ( *,*) (t,vcl(1,i),vcl(2,i),i=nvcin+1,nvc)
  read ( *,*) npolg
  read ( *,*) (t,i=1,npolg)
  read ( *,*) (t,i=1,npolg)
  read ( *,*) nvert
  read ( *,*) ((t,j=1,5),x,t,t,i=1,nvert)
  nvcin = nvc
  read ( *,*) nvc
  read ( *,*) (t,vcl(1,i),vcl(2,i),i=nvcin+1,nvc)
  read ( *,*) (t,i=1,npolg)
  read ( *,*) ntri
  read ( *,*) (t,(til(j,i),j=1,3),(tnbr(j,i),j=1,3),i=1,ntri)
!
  write ( *, '(a)' ) '\\begin{figure}[htbp]'
  write ( *, '(a)' ) '\\begin{center}'
  write ( *, '(a)' ) '\\mbox{\\beginpicture'
  write ( *, '(a)' ) '\\setcoordinatesystem units <1.0cm,1.0cm>'
  write ( *,610) xmin,xmax,ymin,ymax
  write ( *, '(a)' ) '\\setlinear'

  do i = 1, ntri

   do j = 1, 3

      if ( tnbr(j,i) < i ) then

        if ( j <= 2 ) then
          jp1 = j + 1
        else
          jp1 = 1
        end if

        a = til(j,i)
        b = til(jp1,i)
        if ( xmin <= vcl(1,a) .and. vcl(1,a) <= xmax .and. &
             ymin <= vcl(2,a) .and. vcl(2,a) <= ymax .and. &
             xmin <= vcl(1,b) .and. vcl(1,b) <= xmax .and. &
             ymin <= vcl(2,b) .and. vcl(2,b) <= ymax) then
              write ( *,620) vcl(1,a),vcl(2,a),vcl(1,b),vcl(2,b)
         end if

       end if

     end do

  end do

  write ( *, '(a)' ) '\\endpicture}'
  write ( *, '(a)' ) '\\end{center}'
  write ( *, '(a)' ) '\\end{figure}'

   40 continue

  610 format ('\\setplotarea x from ',f9.4,' to ',f9.4,', y from ',f9.4, &
    ' to ',f9.4)
  620 format ('\\plot ',4f15.7,' /')

  return
end
subroutine pitri ( )

!*****************************************************************************80
!
!! PITRI generates PIC commands from DRTRI output.
!
  implicit none

  integer ( kind = 4 ), parameter :: maxti = 8000
  integer ( kind = 4 ), parameter :: maxvc = 5000

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jp1
  integer ( kind = 4 ) msglvl
  integer ( kind = 4 ) nh
  integer ( kind = 4 ) nhola
  integer ( kind = 4 ) nhole
  integer ( kind = 4 ) npolg
  integer ( kind = 4 ) ntri
  integer ( kind = 4 ) nvc
  integer ( kind = 4 ) nvcin
  integer ( kind = 4 ) nvert
  integer ( kind = 4 ) t
  integer ( kind = 4 ) til(3,maxti)
  integer ( kind = 4 ) tnbr(3,maxti)
  real ( kind = 8 ) vcl(2,maxvc)
  real ( kind = 8 ) x
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin

  read ( *,*,end=40) xmin,xmax,ymin,ymax
  read ( *,*) msglvl

  if ( msglvl /= 0 ) then
    return
  end if

  read ( *, * ) nvc
  read ( *, * ) (t,vcl(1,i),vcl(2,i),i=1,nvc)
  read ( *, * ) npolg,nhole
  read ( *, * ) (t,i=1,npolg+nhole)
  read ( *, * ) (t,i=1,npolg)
  read ( *, * ) nvert
  read ( *, * ) ((t,j=1,5),x,i=1,nvert)
  read ( *, * ) nhola,nh
  read ( *, * ) (t,i=1,nh)
  nvcin = nvc
  read ( *, * ) nvc
  read ( *, * ) (t,vcl(1,i),vcl(2,i),i=nvcin+1,nvc)
  read ( *, * ) npolg
  read ( *, * ) (t,i=1,npolg)
  read ( *, * ) (t,i=1,npolg)
  read ( *, * ) nvert
  read ( *, * ) ((t,j=1,5),x,t,t,i=1,nvert)
  nvcin = nvc
  read ( *, * ) nvc
  read ( *, * ) (t,vcl(1,i),vcl(2,i),i=nvcin+1,nvc)
  read ( *, * ) (t,i=1,npolg)
  read ( *, * ) ntri
  read ( *, * ) (t,(til(j,i),j=1,3),(tnbr(j,i),j=1,3),i=1,ntri)

  write ( *, '(a)' ) '.ps 5.5i'

  do i = 1, ntri

    do j = 1, 3

      if ( tnbr(j,i) < i ) then

        if ( j <= 2 ) then
          jp1 = j + 1
        else
          jp1 = 1
        end if

        a = til(j,i)
        b = til(jp1,i)

        if ( xmin <= vcl(1,a) .and. vcl(1,a) <= xmax .and. &
             ymin <= vcl(2,a) .and. vcl(2,a) <= ymax .and. &
             xmin <= vcl(1,b) .and. vcl(1,b) <= xmax .and. &
             ymin <= vcl(2,b) .and. vcl(2,b) <= ymax) then
              write ( *,620) -vcl(2,a),vcl(1,a),-vcl(2,b),vcl(1,b)
        end if
!
!  Rotate by 90 degrees in above line.
!

      end if

    end do

  end do

  write ( *, '(a)' ) '.pe'

40    continue

  620 format ('line from ',f11.7,', ',f11.7,' to ', f11.7,', ',f11.7)

  return
end
subroutine test20 ( )

!*****************************************************************************80
!
!! TEST20 tests XEDGE.
!
  implicit none

  logical intsct
  integer ( kind = 4 ) mode
  real ( kind = 8 ) xu
  real ( kind = 8 ) xv1
  real ( kind = 8 ) xv2
  real ( kind = 8 ) xw1
  real ( kind = 8 ) xw2
  real ( kind = 8 ) yu
  real ( kind = 8 ) yv1
  real ( kind = 8 ) yv2
  real ( kind = 8 ) yw1
  real ( kind = 8 ) yw2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST20'
  write ( *, '(a)' ) '  XEDGE determines whether two edges,'
  write ( *, '(a)' ) '  or an edge and a ray, intersect.'
  write ( *, '(a)' ) '  (An edge is a finite line segment.)'
  write ( *, '(a)' ) '  (A ray is a semi-infinite line segment.)'

  mode = 0

  xv1 = 3.0D+00
  yv1 = 0.0D+00
  xv2 = 3.0D+00
  yv2 = 2.0D+00

  xw1 = 0.0D+00
  yw1 = 0.0D+00
  xw2 = 6.0D+00
  yw2 = 2.0D+00

  write ( *, '(a)' ) ' '

  if ( mode == 0 ) then
    write ( *, '(a)' ) '  Edge 1 is from'
    write ( *, * ) '  (', xv1, ',', yv1, ') to'
    write ( *, * ) '  (', xv2, ',', yv2, ').'
  else
    write ( *, '(a)' ) '  Ray 1 is from'
    write ( *, * ) '  (', xv1, ',', yv1, ') through'
    write ( *, * ) '  (', xv2, ',', yv2, ').'
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Edge 2 is from'
  write ( *, * ) '  (', xw1, ',', yw1, ') to'
  write ( *, * ) '  (', xw2, ',', yw2, ').'

  call xedge ( mode, xv1, yv1, xv2, yv2, xw1, yw1, xw2, yw2, xu, yu, intsct )

  if ( .not. intsct ) then
    write ( *, '(a)' ) 'The objects either do not intersect'
    write ( *, '(a)' ) 'or are parallel, or coincide.'
  else
    write ( *, '(a,2g14.6)' ) 'The objects intersect at ', xu, yu
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '(Expecting the answer (3,1) ).'

  return
end
subroutine test21 ( )

!*****************************************************************************80
!
!! TEST21 tests XLINE.
!
  implicit none

  real ( kind = 8 ) dv
  real ( kind = 8 ) dw
  logical parall
  real ( kind = 8 ) xu
  real ( kind = 8 ) xv1
  real ( kind = 8 ) xv2
  real ( kind = 8 ) xw1
  real ( kind = 8 ) xw2
  real ( kind = 8 ) yu
  real ( kind = 8 ) yv1
  real ( kind = 8 ) yv2
  real ( kind = 8 ) yw1
  real ( kind = 8 ) yw2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST21'
  write ( *, '(a)' ) '  XLINE finds the intersection of two lines.'
  write ( *, '(a)' ) '  Each line is defined as the line a given'
  write ( *, '(a)' ) '  distance to the left of a line through two'
  write ( *, '(a)' ) '  points.'

  xv1 = 0.0D+00
  yv1 = 0.0D+00
  xv2 = 0.0D+00
  yv2 = 1.0D+00
  dv = - 6.0D+00

  xw1 = 0.0D+00
  yw1 = 0.0D+00
  xw2 = 3.0D+00
  yw2 = 1.0D+00
  dw = 0.0D+00

  write ( *, '(a)' ) ' '
  write ( *, * ) '  Line 1 is ', dv, ' units left of the line'
  write ( *, * ) '  through (', xv1, ',', yv1, ') and'
  write ( *, * ) '  (', xv2, ',', yv2, ').'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' '
  write ( *, * ) '  Line 2 is ', dw, ' units left of the line'
  write ( *, * ) '  through (', xw1, ',', yw1, ') and'
  write ( *, * ) '  (', xw2, ',', yw2, ').'

  call xline ( xv1, yv1, xv2, yv2, xw1, yw1, xw2, yw2, dv, dw, xu, yu, parall )

  write ( *, '(a)' ) ' '

  if ( parall ) then
    write ( *, '(a)' ) '  The lines are parallel or coincide.'
  else
    write ( *, '(a,2g14.6)' ) '  The lines intersect at ', xu, yu
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  (Expecting the answer (6,2) ).'

  return
end
