program main

!*****************************************************************************80
!
!! MAIN is the main program for TRIPACK_PRB.
!
!  Discussion:
!
!    TRIPACK_PRB is a test problem for TRIPACK.
!
!    This driver tests software package TRIPACK for constructing a
!    constrained Delaunay triangulation of a set of points in the plane.
!    All modules other than TRMSHR are tested unless an error is
!    encountered, in which case the program terminates immediately.
!
!    By default, tests are performed on a simple data set
!    consisting of 12 nodes whose convex hull covers the unit
!    square.  The data set includes a single constraint region
!    consisting of four nodes forming a smaller square at the
!    center of the unit square.
!
!    A data set consists of the following sequence of records:
!
!    N = Number of nodes.
!
!    NCC = Number of constraint curves.
!
!    (LCC(I), I = 1,NCC) = Indexes of the first node in each
!             constraint curve (format I4).  1 .LE. LCC(1)
!             and, for 1 < I, LCC(I-1) + 3 .LE. LCC(I)
!             .LE. N-2. (Each constraint curve has at least
!             three nodes.)
!
!    (X(I),Y(I), I = 1,N) = Nodal coordinates with non-
!                constraint nodes followed by the NCC
!                sequences of constraint nodes (format
!                2F13.8).
!
!
!  Local parameters:
!
!   PLTSIZ is the plot size for the triangulation plot.
!
!   TOL is the tolerance for TRMTST and NEARND: an upper bound on squared
!   distances.
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 12
  integer ( kind = 4 ), parameter :: ncc_max = 1
  integer ( kind = 4 ), parameter :: nrow = 9

  integer ( kind = 4 ), parameter :: lwk = 2 * n
  integer ( kind = 4 ), parameter :: n6 = 6 * n
  integer ( kind = 4 ), parameter :: ntmx = 2 * n

  real ( kind = 8 ) a
  real ( kind = 8 ) areap
  real ( kind = 8 ) armax
  real ( kind = 8 ) ds(n)
  real ( kind = 8 ) dsq
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) io1
  integer ( kind = 4 ) io2
  integer ( kind = 4 ) k
  integer ( kind = 4 ) ksum
  integer ( kind = 4 ), dimension ( ncc_max ) :: lcc = (/ 9 /)
  integer ( kind = 4 ) lct(ncc_max)
  integer ( kind = 4 ) lend(n)
  integer ( kind = 4 ), parameter :: lin = 1
  integer ( kind = 4 ) list(n6)
  integer ( kind = 4 ) lnew
  integer ( kind = 4 ), parameter :: lplt = 3
  integer ( kind = 4 ) lptr(n6)
  integer ( kind = 4 ) ltri(nrow,ntmx)
  integer ( kind = 4 ) lw
  integer ( kind = 4 ) n0
  integer ( kind = 4 ) na
  integer ( kind = 4 ) nb
  integer ( kind = 4 ) :: ncc = 1
  integer ( kind = 4 ) nearnd
  integer ( kind = 4 ) nn
  integer ( kind = 4 ) nodes(lwk)
  integer ( kind = 4 ) nt
  logical numbr
  real ( kind = 8 ), parameter :: pltsiz = 7.5D+00
  logical prntx
  character ( len = 80 ) title
  real ( kind = 8 ), parameter :: tol = 0.001D+00
  real ( kind = 8 ) wx1
  real ( kind = 8 ) wx2
  real ( kind = 8 ) wy1
  real ( kind = 8 ) wy2
  real ( kind = 8 ), dimension ( n ) :: x = (/ &
    0.00D+00, 1.00D+00, 0.50D+00, 0.15D+00, 0.85D+00, &
    0.50D+00, 0.00D+00, 1.00D+00, 0.35D+00, 0.65D+00, &
    0.65D+00, 0.35D+00 /)
  real ( kind = 8 ), dimension ( n ) :: y = (/ &
    0.00D+00, 0.00D+00, 0.15D+00, 0.50D+00, 0.50D+00, &
    0.85D+00, 1.00D+00, 1.00D+00, 0.35D+00, 0.35D+00, &
    0.65D+00, 0.65D+00 /)
!
!  Print a heading.
!
  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TRIPACK_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the TRIPACK library.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  The number of nodes is ', n

  if ( n < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRIPACK_PRB - Fatal error!'
    write ( *, '(a)' ) '  N must be at least 3!'
    stop
  end if

  if ( ncc < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRIPACK_PRB - Fatal error!'
    write ( *, '(a)' ) '  NCC must be at least 0!'
    stop
  end if
!
!  Create the Delaunay triangulation (TRMESH), and test
!  for errors (refer to TRMTST below).  NODES and DS are
!  used as work space.
!
  call trmesh ( n, x, y, list, lptr, lend, lnew, nodes, nodes(n+1), ds, ier )

  if ( ier == -2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRIPACK_PRB - Fatal error!'
    write ( *, '(a)' ) '  Error in TRMESH:'
    write ( *, '(a)' ) '  The first three nodes are collinear.'
    stop
  else if ( ier == -4 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRIPACK_PRB - Fatal error!'
    write ( *, '(a)' ) '  Error in TRMESH:'
    write ( *, '(a)' ) '  Invalid triangulation.'
    stop
  else if ( 0 < ier ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRIPACK_PRB - Fatal error!'
    write ( *, '(a)' ) '  Error in TRMESH:'
    write ( *, '(a)' ) '  Duplicate nodes encountered.'
    stop
  end if

  call trmtst ( n, x, y, list, lptr, lend, lnew, tol, armax, ier )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TRIPACK_PRB:'
  write ( *, '(a)' ) '  TRMTST reports that the maximum triangle'
  write ( *, '(a,g14.6)' ) '  aspect ratio is ', armax

  if ( 0 < ier ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRIPACK_PRB - Fatal error!'
    write ( *, '(a)' ) '  Error in TRMTST:'
    stop
  end if
!
!  Add the constraint curves (ADDCST).  Note that edges
!  and triangles are not removed from constraint regions.
!  ADDCST forces the inclusion of triangulation edges
!  connecting the sequences of constraint nodes.  If it
!  is necessary to alter the triangulation, the empty
!  circumcircle property is no longer satisfied.
!
  lw = lwk
  call addcst ( ncc, lcc, n, x, y, lw, nodes, list, lptr, lend, ier )

  if ( ier /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRIPACK_PRB - Fatal error!'
    write ( *, '(a,i6)' ) '  Error in ADDCST, IER = ', ier
    stop
  end if

  if ( lw == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRIPACK_PRB: Note'
    write ( *, '(a)' ) '  Subroutine EDGE was not tested, because'
    write ( *, '(a)' ) '  no edges were swapped by ADDCST.'
  end if
!
!  Test TRPRNT, TRLIST, and TRLPRT, and TRPLOT.
!
  prntx = .true.

  call trprnt ( ncc, lcc, n, x, y, list, lptr, lend, prntx )

  call trlist ( ncc, lcc, n, list, lptr, lend, nrow, nt, ltri, lct, ier )

  call trlprt ( ncc, lct, n, x, y, nrow, nt, ltri, prntx )
!
!  Set the plot window [WX1,WX2] X [WY1,WY2].
!
  wx1 = minval ( x(1:n) )
  wx2 = maxval ( x(1:n) )
  wy1 = minval ( y(1:n) )
  wy2 = maxval ( y(1:n) )
!
!  Create an encapsulated PostScript file image of the triangulation.
!
  open ( unit = lplt, file = 'tripack_prb.eps', status = 'replace' )
!
!  NUMBR = TRUE iff nodal indexes are to be displayed.
!
  numbr = ( n <= 200 )
!
!  Store a plot title.  It must be enclosed in parentheses.
!
  title = '(Triangulation created by TRIPACK_PRB)'

  call trplot ( lplt, pltsiz, wx1, wx2, wy1, wy2, ncc, lcc, n, &
    x, y, list, lptr, lend, title, numbr, ier )

  close ( unit = lplt )

  if ( ier == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRIPACK_PRB:'
    write ( *, '(a)' ) '  TRPLOT has created a triangulation plot.'
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRIPACK_PRB:'
    write ( *, '(a,i6)' ) '  TRPLOT failed with IER = ', ier
  end if
!
!  Test BNODES.
!
  call bnodes ( n, list, lptr, lend, nodes, nb, na, nt )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TRIPACK_PRB:'
  write ( *, '(a)' ) '  Output from BNODES'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of boundary nodes = ', nb
  write ( *, '(a,i6)' ) '  Number of edges          = ', na
  write ( *, '(a,i6)' ) '  Number of triangles      = ', nt
!
!  Test AREAP.
!
  a = areap ( x, y, nb, nodes )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TRIPACK_PRB:'
  write ( *, '(a)' ) '  Output from AREAP'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Area of convex hull is ', a
!
!  Test GETNP by ordering the nodes on distance from N0
!  and verifying the ordering.
!
  n0 = n/2
  nodes(1) = n0
  ds(1) = 0.0D+00
  ksum = n0

  do k = 2, n

    call getnp ( ncc, lcc, n, x, y, list, lptr, lend, &
      k, nodes, ds, ier )

    if ( ier /= 0 .or. ds(k) < ds(k-1) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TRIPACK_PRB - Fatal error!'
      write ( *, '(a)' ) '  Error in GETNP.'
      stop
    end if

    ksum = ksum + nodes(k)

  end do
!
!  Test for all nodal indexes included in NODES.
!
  if ( ksum /= (n*(n+1))/2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRIPACK_PRB - Fatal error!'
    write ( *, '(a)' ) '  Error in GETNP.'
    stop
  end if
!
!  Test NEARND by verifying that the nearest node to K is
!  node K for K = 1 to N.
!
  do k = 1, n

    n0 = nearnd ( x(k), y(k), 1, n, x, y, list, lptr, lend, dsq )

    if ( n0 /= k .or. tol < dsq ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TRIPACK_PRB - Fatal error!'
      write ( *, '(a)' ) '  Error in NEARND.'
      stop
    end if

  end do
!
!  Test DELARC by removing a boundary arc if possible.
!  The first two nodes define a boundary arc
!  in the default data set.
!
  io1 = 1
  io2 = 2
  call delarc ( n, io1, io2, list, lptr, lend, lnew, ier )

  if ( ier == 1 .or. ier == 4 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRIPACK_PRB - Fatal error!'
    write ( *, '(a,i6)' ) '  Error in DELARC, IER = ', ier
    stop
  end if

  if ( ier /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRIPACK_PRB: Note'
    write ( *, '(a)' ) '  Subroutine DELARC was not tested, because'
    write ( *, '(a)' ) '  nodes 1 and 2 do not form a removable'
    write ( *, '(a)' ) '  boundary edge.'
  end if
!
!  Recreate the triangulation without constraints.
!
  call trmesh ( n, x, y, list, lptr, lend, lnew, nodes, &
    nodes(n+1), ds, ier )

  ncc = 0
!
!  Test DELNOD by removing nodes 4 to N (in reverse order).
!
  if ( n <= 3 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRIPACK_PRB:'
    write ( *, '(a)' ) '  Subroutine DELNOD was not tested.'
    write ( *, '(a)' ) '  N cannot be reduced below 3.'

  else

    nn = n

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRIPACK_PRB:'
    write ( *, '(a)' ) '  Test DELNOD.'

    write ( *, '(a)' ) ' '

    do

      k = nn
      lw = lwk / 2

      call delnod ( nn, ncc, lcc, nn, x, y, list, lptr, lend, &
        lnew, lw, nodes, ier )

      if ( ier /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TRIPACK_PRB - Fatal error!'
        write ( *, '(a,i6)' ) '  Error in DELNOD, IER = ', ier
        stop
      end if

      write ( *, '(a,i6)' ) '  DELNOD deleted node ', nn+1

      if ( nn <= 3 ) then
        exit
      end if

    end do

  end if
!
!  Successful test.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TRIPACK_PRB:'
  write ( *, '(a)' ) '  No triangulation errors were encountered.'
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TRIPACK_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
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

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
