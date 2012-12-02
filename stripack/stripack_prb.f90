program main

!*****************************************************************************80
!
!! MAIN is the main program for STRIPACK_PRB.
!
!  Discussion:
!
!    STRIPACK_PRB is a test routine for STRIPACK.
!
!  Modified:
!
!    25 February 2007
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'STRIPACK_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the STRIPACK library.'

  call test01 ( )
  call test02 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'STRIPACK_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 is a test for STRIPACK.
!
!  Discussion:
!
!    This driver tests software package STRIPACK for constructing a 
!    Delaunay triangulation and Voronoi diagram of a set of nodes on 
!    the surface of the unit sphere.
!
!    All STRIPACK subprograms are tested.
!
!    By default, a triangulation is created from a set of N nodes consisting 
!    of the north pole and N-1 points uniformly distributed around the 
!    60-degree parallel (with constant longitudinal separation).  
!
!    The data is stored as RLAT(I), RLON(I), which are the nodal coordinates 
!    in degrees latitude (-90 to 90) and degrees longitude (-180 to 180).
!
!  Modified:
!
!    25 February 2007
!
  implicit none

  integer ( kind = 4 ), parameter :: nmax = 200
  integer ( kind = 4 ), parameter :: nrow = 9

  real    ( kind = 8 ) a
  real    ( kind = 8 ) al
  real    ( kind = 8 ) area
  real    ( kind = 8 ) areas
  real    ( kind = 8 ) ds(nmax)
  real    ( kind = 8 ) elat
  real    ( kind = 8 ) elon
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) iflag
  logical inside
  integer ( kind = 4 ) iwk(2*nmax)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) ksum
  integer ( kind = 4 ) kt
  integer ( kind = 4 ) lbtri(6,nmax)
  integer ( kind = 4 ) lend(nmax)
  integer ( kind = 4 ) list(6*nmax)
  integer ( kind = 4 ) listc(6*nmax)
  integer ( kind = 4 ) lnew
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lpl
  integer ( kind = 4 ), parameter :: lplt = 3
  integer ( kind = 4 ), parameter :: lplv = 4
  integer ( kind = 4 ) lptr(6*nmax)
  integer ( kind = 4 ) ltri(nrow,2*nmax-4)
  integer ( kind = 4 ) lwk
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n0
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n3
  integer ( kind = 4 ) na
  integer ( kind = 4 ) nb
  integer ( kind = 4 ) nearnd
  integer ( kind = 4 ) nn
  integer ( kind = 4 ) nt
  logical numbr
  integer ( kind = 4 ) nv
  real    ( kind = 8 ) p(3)
  real    ( kind = 8 ), parameter :: pltsiz = 7.5D+00
  real    ( kind = 8 ) rc(2*nmax-4)
  real    ( kind = 8 ) rlat(nmax)
  real    ( kind = 8 ) rlon(nmax)
  real    ( kind = 8 ) sc
  character ( len = 80 ) trplot_file_name
  character ( len = 80 ) trplot_title
  real    ( kind = 8 ) v1(3)
  real    ( kind = 8 ) v2(3)
  real    ( kind = 8 ) v3(3)
  real    ( kind = 8 ) vlat
  real    ( kind = 8 ) vlon
  real    ( kind = 8 ) vnrm
  character ( len = 80 ) vrplot_file_name
  character ( len = 80 ) vrplot_title
  real    ( kind = 8 ) x(nmax)
  real    ( kind = 8 ) xc(2*nmax-4)
  real    ( kind = 8 ) y(nmax)
  real    ( kind = 8 ) yc(2*nmax-4)
  real    ( kind = 8 ) z(nmax)
  real    ( kind = 8 ) zc(2*nmax-4)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  TRANS converts Cartesian to spherical coordinates.'
  write ( *, '(a)' ) '  TRMESH creates a triangulation.'
  write ( *, '(a)' ) '  TRPRNT prints out a triangulation.'
  write ( *, '(a)' ) '  TRLIST creates a triangle list.'
  write ( *, '(a)' ) '  TRLPRT prints a triangle list.'
  write ( *, '(a)' ) '  TRPLOT plots a triangulation.'
  write ( *, '(a)' ) '  AREAS computes areas.'
  write ( *, '(a)' ) '  BNODES computes boundary nodes.'
  write ( *, '(a)' ) '  GETNP gets the next nearest node to a given node.'
  write ( *, '(a)' ) '  NEARND returns the nearest node to a given point.'
  write ( *, '(a)' ) '  DELARC removes a boundary arc if possible.'
  write ( *, '(a)' ) '  CRLIST constructs the Voronoi diagram.'
  write ( *, '(a)' ) '  VRPLOT plots the Voronoi diagram.'
  write ( *, '(a)' ) '  SCOORD prints the Voronoi region boundary associated'
  write ( *, '(a)' ) '    with a point.'
  write ( *, '(a)' ) '  INSIDE determines if a point is inside a '
  write ( *, '(a)' ) '    Voronoi region.'
!
!  Generate the default set of nodes as latitudinal and longitudinal
!  coordinates. 
!
  n = 100

  call random_number ( harvest = rlat(1:n) )
  call random_number ( harvest = rlon(1:n) )

  rlat(1:n) = ( ( 1.0D+00 - rlat(1:n) ) * (  -90.0D+00 ) &
              +             rlat(1:n)   *     90.0D+00 )

  rlon(1:n) = ( ( 1.0D+00 - rlon(1:n) ) * ( -180.0D+00 ) &
              +             rlon(1:n)   *    180.0D+00 )

  if ( n < 3 .or. nmax < n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST01 - Fatal error!'
    write ( *, '(a)' ) '  The value of N is illegal.'
    write ( *, '(a,i8,a)' ) '  3 <= N <= NMAX = ', nmax, ' is required.'
    write ( *, '(a,i8)' ) '  Input N = ', n
    stop
  end if
!
!  Set X and Y to the values of RLON and RLAT, respectively,
!  in radians.  (RLON and RLAT are saved for printing by TRPRNT).
!
  sc = atan ( 1.0D+00 ) / 45.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I     RLON     RLAT'
  write ( *, '(a)' ) ' '
  do i = 1, 5
    write ( *, '(2x,i8,2x,f10.6,2x,f10.6)' ) i, rlon(i), rlat(i)
  end do

  x(1:n) = sc * rlon(1:n)
  y(1:n) = sc * rlat(1:n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I     X     Y'
  write ( *, '(a)' ) ' '
  do i = 1, 5
    write ( *, '(2x,i8,2x,f10.6,2x,f10.6)' ) i, x(i), y(i)
  end do
!
!  Transform spherical coordinates X and Y to Cartesian
!  coordinates (X,Y,Z) on the unit sphere (X**2 + Y**2 + Z**2 = 1).
!
  call trans ( n, y, x, x, y, z )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I     X     Y     Z'
  write ( *, '(a)' ) ' '
  do i = 1, 5
    write ( *, '(2x,i8,2x,f10.6,2x,f10.6,2x,f10.6)' ) i, x(i), y(i), z(i)
  end do
!
!  Create the triangulation.
!
  call trmesh ( n, x, y, z, list, lptr, lend, lnew, iwk, iwk(n+1), ds, ier )

  if ( ier == -2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST01 - Warning!'
    write ( *, '(a)' ) '  Error in TRMESH.'
    write ( *, '(a)' ) '  The first three nodes are collinear.'
    stop
  end if

  if ( 0 < ier ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST01 - Fatal error!'
    write ( *, '(a)' ) '  Error in TRMESH.'
    write ( *, '(a)' ) '  Duplicate nodes encountered.'
    stop
  end if
!
!  Print the spherical coordinates and adjacency information.
!
!  0 < IFLAG indicates that RLON and RLAT only are to be printed.
!
  iflag = 1

  call trprnt ( n, rlon, rlat, z, iflag, list, lptr, lend )
!
!  Test TRLIST and TRLPRT by creating and printing a triangle list.
!
  call trlist ( n, list, lptr, lend, nrow, nt, ltri, ier )

  call trlprt ( n, rlon, rlat, z, iflag, nrow, nt, ltri )
!
!  Test TRPLOT by plotting the portion of the triangulation contained 
!  in the hemisphere centered at E = (ELAT,ELON), where ELAT and ELON
!  are taken to be the center of the range of
!  the nodal latitudes and longitudes.
!
  elat = minval ( rlat(1:n) )
  vlat = maxval ( rlat(1:n) )
  elon = minval ( rlon(1:n) )
  vlon = maxval ( rlon(1:n) )

  elat = ( elat + vlat ) / 2.0D+00
  elon = ( elon + vlon ) / 2.0D+00
  a = 90.0D+00
  numbr = n <= 200

  trplot_title = '(Triangulation created by STRIPACK_PRB)'

  trplot_file_name = 'stripack_prb_del.eps'

  open ( lplt, file = trplot_file_name )

  call trplot ( lplt, pltsiz, elat, elon, a, n, x, y, z, list, &
    lptr, lend, trplot_title, numbr, ier )

  if ( ier == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  TRPLOT created the triangulation plot file: "' // &
      trim ( trplot_file_name ) // '".'
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  TRPLOT returned error code ', ier
  end if
!
!  Test AREAS by computing and printing the area of the
!  convex hull of the nodes (sum of triangle
!  areas) relative to the total surface area (4*Pi).
!
  area = 0.0D+00

  do kt = 1, nt
    n1 = ltri(1,kt)
    n2 = ltri(2,kt)
    n3 = ltri(3,kt)
    v1(1) = x(n1)
    v1(2) = y(n1)
    v1(3) = z(n1)
    v2(1) = x(n2)
    v2(2) = y(n2)
    v2(3) = z(n2)
    v3(1) = x(n3)
    v3(2) = y(n3)
    v3(3) = z(n3)
    area = area + areas ( v1, v2, v3 )
  end do

  area = area / ( 16.0D+00 * atan ( 1.0D+00 ) )

  write ( *, '(a)' ) ' '
  write ( *, '(a,f8.2)' ) '  Relative area of convex hull = ', area
!
!  Test BNODES.  The ordered sequence of boundary nodes is stored in IWK.
!
  call bnodes ( n, list, lptr, lend, iwk, nb, na, nt )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Output from BNODES:'
  write ( *, '(a,i8)' ) '  Number of boundary nodes = ', nb
  write ( *, '(a,i8)' ) '  Number of arcs =           ', na
  write ( *, '(a,i8)' ) '  Number of triangles =      ', nt
!
!  Test GETNP by ordering the nodes on distance from N0 and verifying 
!  the ordering.  
!
!  The sequence of nodal indexes is stored in IWK, and the values of an
!  increasing function (the negative cosine) of angular distance is
!  stored in DS.
!
  n0 = n / 2
  iwk(1) = n0
  ds(1) = -1.0D+00
  ksum = n0

  do k = 2, n

    call getnp ( x, y, z, list, lptr, lend, k, iwk, ds(k), ier )

    if ( ier /= 0  .or.  ds(k) < ds(k-1) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01 - Fatal error!'
      write ( *, '(a)' ) '  Error in GETNP.'
      stop
    end if

    ksum = ksum + iwk(k)

  end do
!
!  Test for all nodal indexes included in IWK.
!
  if ( ksum /= ( n * ( n + 1 ) ) / 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST01 - Fatal error!'
    write ( *, '(a)' ) '  Error in GETNP.'
    stop
  end if
!
!  Test NEARND by verifying that the nearest node to K is
!  node K for K = 1 to N.
!
  do k = 1, n

    p(1) = x(k)
    p(2) = y(k)
    p(3) = z(k)

    n0 = nearnd ( p, 1, n, x, y, z, list, lptr, lend, al )

    if ( n0 /= k .or. 0.001D+00 < al ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01 - Fatal error!'
      write ( *, '(a)' ) '  Error in NEARND.'
      stop
    end if

  end do
!
!  Test DELARC by removing a boundary arc if possible.
!  The last two nodes define a boundary arc
!  in the default data set.
!
  n1 = n - 1
  n2 = n
  call delarc ( n, n1, n2, list, lptr, lend, lnew, ier )

  if ( ier == 1  .or.  ier == 4 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST01 - Warning!'
    write ( *, '(a,i8)' ) '  DELARC returned error code ', ier
    stop
  end if

  if ( ier /= 0 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Subroutine DELARC was not tested.'
    write ( *, '(a,i8,a,i8,a)' ) '  Nodes ', n1, ' and ', n2, &
      ' do not form a removable boundary arc.'
  else

    call trmesh ( n, x, y, z, list, lptr, lend, lnew, iwk, iwk(n+1), ds, &
      ier )

  end if
!
!  Test CRLIST, VRPLOT, and SCOORD by constructing and
!  plotting the Voronoi diagram, and printing
!  the Voronoi region boundary (ordered
!  sequence of Voronoi vertices) associated with N0.
!
!  Note that the triangulation data structure
!  is altered if 0 < NB.
!
  call crlist ( n, nmax, x, y, z, list, lend, lptr, lnew, &
    lbtri, listc, nb, xc, yc, zc, rc, ier )

  if ( ier /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST01 - Warning!'
    write ( *, '(a,i8)' ) '  CRLIST returned error code ', ier
    stop
  end if
!
!  Use the same parameter values that were used for the
!  triangulation plot (except the output unit and title).
!
  nt = 2 * n - 4

  vrplot_file_name = 'stripack_prb_vor.eps'

  vrplot_title = '(Voronoi diagram created by STRIPACK_PRB)'

  open ( unit = lplv, file = vrplot_file_name )

  call vrplot ( lplv, pltsiz, elat, elon, a, n, x, y, z, nt, listc, &
    lptr, lend, xc, yc, zc, vrplot_title, numbr, ier )

  if ( ier == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  VRPLOT created the Voronoi plot file: "' // &
      trim ( vrplot_file_name ) // '".'
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST01 - Warning!'
    write ( *, '(a,i8)' ) '  VRPLOT returned error code ', ier
  end if

  n0 = 1

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Voronoi region for node ', n0
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     Triangle     Latitude     Longitude' // &
    '     Circumradius'
  write ( *, '(a)' ) ' '
!
!  Initialize for loop on Voronoi vertices (triangle circumcenters).  
!  The number of vertices is accumulated in NV, and the vertex indexes
!  are stored in IWK.  The vertices are converted to latitude and longitude 
!  in degrees for printing.
!
  nv = 0
  lpl = lend(n0)
  lp = lpl

  do

    lp = lptr(lp)
    kt = listc(lp)
    nv = nv + 1
    iwk(nv) = kt
    call scoord ( xc(kt), yc(kt), zc(kt), vlat, vlon, vnrm )
    vlat = vlat / sc
    vlon = vlon / sc
    write ( *, '(i13,f13.6,f14.6,f17.6)' ) kt, vlat, vlon, rc(kt)

    if ( lp == lpl ) then
      exit
    end if

  end do
!
!  Test INSIDE by checking for node N0 inside its Voronoi region.
!
  p(1) = x(n0)
  p(2) = y(n0)
  p(3) = z(n0)

  if ( .not. inside ( p, 2*nmax-4, xc, yc, zc, nv, iwk, ier ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST01 - Warning!'
    write ( *, '(a)' ) '  Error in INSIDE.'
    write ( *, '(a)' ) '  A node is not contained in its Voronoi region.'
    write ( *, '(a,i8)' ) '  Node index = ', n0
  end if

  if ( ier /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST01 - Fatal error!'
    write ( *, '(a)' ) '  Error in INSIDE.'
    write ( *, '(a,i8)' ) '  IER = ', ier
    stop
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  INSIDE correctly reports that node N0 is'
  write ( *, '(a)' ) '  inside its Voronoi region!'

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests EDGE and DELNOD.
!
!  Modified:
!
!    16 June 2007
!
  implicit none

  integer ( kind = 4 ), parameter :: nmax = 200

  real    ( kind = 8 ) ds(nmax)
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) iwk(2*nmax)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) lend(nmax)
  integer ( kind = 4 ) list(6*nmax)
  integer ( kind = 4 ) listc(6*nmax)
  integer ( kind = 4 ) lnew
  integer ( kind = 4 ) lptr(6*nmax)
  integer ( kind = 4 ) lwk
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) nn
  real    ( kind = 8 ) rlat(nmax)
  real    ( kind = 8 ) rlon(nmax)
  real    ( kind = 8 ) sc
  real    ( kind = 8 ) x(nmax)
  real    ( kind = 8 ) y(nmax)
  real    ( kind = 8 ) z(nmax)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  EDGE can be used to force an edge between two nodes.'
  write ( *, '(a)' ) '  DEL can be used to delete a node.'
!
!  Generate the default set of nodes as latitudinal and longitudinal
!  coordinates. 
!
  n = 9

  rlat(1) = 90.0D+00
  rlat(2:n) = 60.0D+00

  rlon(1) = 0.0D+00
  do k = 2, n
    rlon(k) = real ( k - 2, kind = 8 ) * 360.0D+00 / real ( n - 1, kind = 8 )
  end do

  if ( n < 3 .or. nmax < n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST02 - Fatal error!'
    write ( *, '(a)' ) '  The value of N is illegal.'
    write ( *, '(a,i8,a)' ) '  3 <= N <= NMAX = ', nmax, ' is required.'
    write ( *, '(a,i8)' ) '  Input N = ', n
    stop
  end if
!
!  Set X and Y to the values of RLON and RLAT, respectively,
!  in radians.  (RLON and RLAT are saved for printing by TRPRNT).
!
  sc = atan ( 1.0D+00 ) / 45.0D+00

  x(1:n) = sc * rlon(1:n)
  y(1:n) = sc * rlat(1:n)
!
!  Transform spherical coordinates X and Y to Cartesian
!  coordinates (X,Y,Z) on the unit sphere (X**2 + Y**2 + Z**2 = 1).
!
  call trans ( n, y, x, x, y, z )
!
!  Create the triangulation.
!
  call trmesh ( n, x, y, z, list, lptr, lend, lnew, iwk, iwk(n+1), ds, ier )

  if ( ier == -2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST02 - Warning!'
    write ( *, '(a)' ) '  Error in TRMESH.'
    write ( *, '(a)' ) '  The first three nodes are collinear.'
    return
  else if ( 0 < ier ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST02 - Fatal error!'
    write ( *, '(a)' ) '  Error in TRMESH.'
    write ( *, '(a)' ) '  Duplicate nodes encountered.'
    return
  end if
!
!  Test EDGE by forcing an edge between nodes N1=1 and N2=N. 
!
  n1 = 1
  n2 = n

  call edge ( n1, n2, x, y, z, nmax, iwk, list, lptr, lend, ier )

  if ( ier /= 0 .and. ier /= 5 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST02 - Fatal error!'
    write ( *, '(a)' ) '  Error in EDGE.'
    write ( *, '(a,i8)' ) '  IER = ', ier
    stop
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  EDGE has forced an edge between two nodes.'
!
!  Test DELNOD by removing nodes 4 to N (in reverse order). 
!
  write ( *, '(a)' ) ' '

  if ( n <= 3 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST02:'
    write ( *, '(a)' ) '  Subroutine DELNOD was not tested, because'
    write ( *, '(a)' ) '  the number of nodes N is too small.'

  else

    nn = n
    lwk = nmax

    do

      k = nn

      write ( *, '(a,i8)' ) '  Call DELNOD to delete node ', k

      call delnod ( k, nn, x, y, z, list, lptr, lend, lnew, lwk, iwk, ier )

      if ( ier /= 0 .and. ier /= 5 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TEST02 - Fatal error!'
        write ( *, '(a,i8)' ) '  DELNOD returned IER = ', ier
        stop
      end if

      if ( nn <= 3 ) then
        exit
      end if

    end do

  end if

  return
end
