subroutine addcst ( ncc, lcc, n, x, y, lwk, iwk, list, lptr, lend, ier )

!*****************************************************************************80
!
!! ADDCST adds constraint curves to a Delaunay triangulation.
!
!  Discussion:
!
!    This subroutine provides for creation of a constrained
!    Delaunay triangulation which, in some sense, covers an
!    arbitrary connected region R rather than the convex hull
!    of the nodes.  This is achieved simply by forcing the
!    presence of certain adjacencies (triangulation arcs) 
!    corresponding to constraint curves.  The union of triangles
!    coincides with the convex hull of the nodes, but triangles
!    in R can be distinguished from those outside of R.  The
!    only modification required to generalize the definition of
!    the Delaunay triangulation is replacement of property 5
!    (refer to TRMESH) by the following:
!
!    5')  If a node is contained in the interior of the 
!         circumcircle of a triangle, then every interior point
!         of the triangle is separated from the node by a
!         constraint arc.
!
!    In order to be explicit, we make the following definitions.  
!    A constraint region is the open interior of a
!    simple closed positively oriented polygonal curve defined
!    by an ordered sequence of three or more distinct nodes
!    (constraint nodes) P(1),P(2),...,P(K), such that P(I) is
!    adjacent to P(I+1) for I = 1,...,K with P(K+1) = P(1).
!    Thus, the constraint region is on the left (and may have
!    nonfinite area) as the sequence of constraint nodes is
!    traversed in the specified order.  The constraint regions
!    must not contain nodes and must not overlap.  The region
!    R is the convex hull of the nodes with constraint regions
!    excluded.
!
!    Note that the terms boundary node and boundary arc are
!    reserved for nodes and arcs on the boundary of the convex
!    hull of the nodes.
!
!    The algorithm is as follows:  given a triangulation
!    which includes one or more sets of constraint nodes, the
!    corresponding adjacencies (constraint arcs) are forced to
!    be present (Subroutine EDGE).  Any additional new arcs
!    required are chosen to be locally optimal (satisfy the
!    modified circumcircle property).
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    Robert Renka,
!    Department of Computer Science,
!    University of North Texas,
!    renka@cs.unt.edu
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 751: TRIPACK, 
!    A Constrained Two-Dimensional Delaunay Triangulation Package,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 1, 1996.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NCC, the number of constraint curves (constraint 
!    regions).  0 <= NCC.
!
!    Input, integer ( kind = 4 ) LCC(NCC) (or dummy array of length 1 if NCC = 0) 
!    containing the index (for X, Y, and LEND) of the first node of 
!    constraint I in LCC(I) for I = 1 to NCC.  Thus, constraint I
!    contains K = LCC(I+1) - LCC(I) nodes, K >= 3, stored in (X,Y) 
!    locations LCC(I), ..., LCC(I+1)-1, where LCC(NCC+1) = N+1.
!
!    Input, integer ( kind = 4 ) N, the number of nodes in the triangulation, including
!    constraint nodes.  3 <= N.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the coordinates of the nodes with
!    non-constraint nodes in the first LCC(1)-1 locations, followed by NCC
!    sequences of constraint nodes.  Only one of these sequences may be
!    specified in clockwise order to represent an exterior constraint curve (a 
!    constraint region with nonfinite area).
!
!    Input/output, integer ( kind = 4 ) LWK.  On input, the length of IWK.  This must be 
!    at least 2*NI where NI is the maximum number of arcs which intersect a
!    constraint arc to be added.  NI is bounded by N-3.  On output, the 
!    required length of IWK unless IER = 1 or IER = 3.  In the case of 
!    IER = 1, LWK is not altered from its input value.
!
!    Output, integer ( kind = 4 ) IWK(LWK), the endpoint indexes of the new arcs which 
!    were swapped in by the last call to subroutine EDGE.
!
!    Input/output, integer ( kind = 4 ) LIST(*), LPTR(*), LEND(N), the data structure 
!    defining the triangulation.  Refer to subroutine TRMESH.  On output,
!    The structure has all constraint arcs present unless IER /= 0.  
!    These arrays are not altered if IER = 1.
!
!    Output, integer ( kind = 4 ) IER = Error indicator:
!    0 if no errors were encountered.
!    1 if NCC, N, or an LCC entry is outside its valid range, or 
!      LWK < 0 on input.
!    2 if more space is required in IWK.
!    3 if the triangulation data structure is invalid, or failure (in 
!      EDGE or OPTIM) was caused by collinear nodes on the convex hull 
!      boundary.  An error message is written to logical unit 6 in
!      this case.
!    4 if intersecting constraint arcs were encountered.
!    5 if a constraint region contains a node.
!
  implicit none

  integer ( kind = 4 ) lwk
  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) ifrst
  integer ( kind = 4 ) ilast
  integer ( kind = 4 ) iwk(lwk)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kbak
  integer ( kind = 4 ) kfor
  integer ( kind = 4 ) kn
  integer ( kind = 4 ) lcc(*)
  integer ( kind = 4 ) lccip1
  integer ( kind = 4 ) lend(n)
  integer ( kind = 4 ) list(*)
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lpb
  integer ( kind = 4 ) lpf
  integer ( kind = 4 ) lpl
  integer ( kind = 4 ) lptr(*)
  integer ( kind = 4 ) lw
  integer ( kind = 4 ) lwd2
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) ncc
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)

  lwd2 = lwk / 2
!
!  Test for errors in input parameters.
!
  ier = 0

  if ( ncc < 0 .or. lwk < 0 ) then
    ier = 1
    return
  end if

  if ( ncc == 0 ) then

    if ( n < 3 ) then
      ier = 1
    else
      ier = 0
      lwk = 0
    end if

    return

  else

    lccip1 = n + 1

    do i = ncc, 1, -1

      if ( lccip1 - lcc(i) < 3 ) then
        ier = 1
        return
      end if

      lccip1 = lcc(i)

    end do

    if ( lccip1 < 1 ) then
      ier = 1
      return
    end if

  end if
!
!  Force the presence of constraint arcs.  The outer loop is
!  on constraints in reverse order.  IFRST and ILAST are
!  the first and last nodes of constraint I.
!
  lwk = 0
  ifrst = n + 1

  do i = ncc, 1, -1

    ilast = ifrst - 1
    ifrst = lcc(i)
!
!  Inner loop on constraint arcs N1-N2 in constraint I.
!
    n1 = ilast

    do n2 = ifrst, ilast

      lw = lwd2

      call edge ( n1, n2, x, y, lw, iwk, list, lptr, lend, ier )

      lwk = max ( lwk, 2 * lw )

      if ( ier == 4 ) then
        ier = 3
      end if

      if ( ier /= 0 ) then
        return
      end if

      n1 = n2
    end do

  end do
!
!  Test for errors.  The outer loop is on constraint I with
!  first and last nodes IFRST and ILAST, and the inner loop
!  is on constraint nodes K with (KBAK,K,KFOR) a subsequence 
!  of constraint I.
!
  ifrst = n + 1

  do i = ncc, 1, -1

    ilast = ifrst - 1
    ifrst = lcc(i)
    kbak = ilast

    do k = ifrst,ilast

      kfor = k + 1

      if ( k == ilast ) then
        kfor = ifrst
      end if
!
!  Find the LIST pointers LPF and LPB of KFOR and KBAK as neighbors of K.
!
      lpf = 0
      lpb = 0
      lpl = lend(k)
      lp = lpl

      do

        lp = lptr(lp)
        kn = abs ( list(lp) )

        if ( kn == kfor ) then
          lpf = lp
        end if

        if ( kn == kbak ) then
          lpb = lp
        end if

        if ( lp == lpl ) then
          exit
        end if

      end do
!
!  A pair of intersecting constraint arcs was encountered
!  if and only if a constraint arc is missing (introduction 
!  of the second caused the first to be swapped out).
!
      if ( lpf == 0 .or. lpb == 0 ) then
        ier = 4
        return
      end if
!
!  Loop on neighbors KN of node K which follow KFOR and
!  precede KBAK.  The constraint region contains no nodes
!  if and only if all such nodes KN are in constraint I.
!
      lp = lpf

      do

        lp = lptr(lp)

        if ( lp == lpb ) then
          exit
        end if

        kn = abs ( list(lp) )

        if ( kn < ifrst .or. ilast < kn ) then
          ier = 5
          return
        end if

      end do

      kbak = k

    end do

  end do

  ier = 0

  return
end
subroutine addnod ( k, xk, yk, ist, ncc, lcc, n, x, y, list, lptr, lend, lnew, &
  ier )

!*****************************************************************************80
!
!! ADDNOD adds a node to a triangulation.
!
!  Discussion:
!
!    Given a triangulation of N nodes in the plane created by
!    subroutine TRMESH or TRMSHR, this subroutine updates the
!    data structure with the addition of a new node in position
!    K.  If node K is inserted into X and Y (K <= N) rather
!    than appended (K = N+1), then a corresponding insertion
!    must be performed in any additional arrays associated
!    with the nodes.  For example, an array of data values Z
!    must be shifted down to open up position K for the new
!    value:  set Z(I+1) to Z(I) for I = N,N-1,...,K.  For
!    optimal efficiency, new nodes should be appended whenever
!    possible.  Insertion is necessary, however, to add a non-
!    constraint node when constraints are present (refer to
!    subroutine ADDCST).
!
!    Note that a constraint node cannot be added by this
!    routine.  In order to insert a constraint node, it is
!    necessary to add the node with no constraints present
!    (call this routine with NCC = 0), update LCC by increment-
!    ing the appropriate entries, and then create (or restore)
!    the constraints by a call to ADDCST.
!
!    The algorithm consists of the following steps:  node K
!    is located relative to the triangulation (TRFIND), its
!    index is added to the data structure (INTADD or BDYADD),
!    and a sequence of swaps (SWPTST and SWAP) are applied to
!    the arcs opposite K so that all arcs incident on node K
!    and opposite node K (excluding constraint arcs) are local-
!    ly optimal (satisfy the circumcircle test).  Thus, if a
!    (constrained) Delaunay triangulation is input, a (con-
!    strained) Delaunay triangulation will result.  All indexes
!    are incremented as necessary for an insertion.
!
!  Modified:
!
!    29 March 2002
!
!  Author:
!
!    Robert Renka,
!    Department of Computer Science,
!    University of North Texas,
!    renka@cs.unt.edu
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 751: TRIPACK, 
!    A Constrained Two-Dimensional Delaunay Triangulation Package,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 1, 1996.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) K, the nodal index (index for X, Y, and LEND) of the
!    new node to be added.  1 <= K <= LCC(1).  (K <= N+1 if NCC=0).
!
!    Input, real ( kind = 8 ) XK, YK, the coordinates of the new node (to be
!    stored in X(K) and Y(K)).  The node must not lie in a constraint region.
!
!    Input, integer ( kind = 4 ) IST, the index of a node at which TRFIND begins the
!    search.  Search time depends on the proximity
!    of this node to node K.  1 <= IST <= N.
!
!    Input, integer ( kind = 4 ) NCC, the number of constraint curves.  0 <= NCC.
!
!    Input/output, integer ( kind = 4 ) LCC(*), list of constraint curve starting indexes 
!    (or dummy array of length 1 if NCC = 0).  Refer to subroutine ADDCST.
!    On output, starting indexes incremented by 1 to reflect the insertion 
!    of K unless NCC = 0 or (IER /= 0 and IER /= -4).
!
!    Input/output, integer ( kind = 4 ) N, the number of nodes in the triangulation.  
!    3 <= N.  Note that N will be incremented following the addition of node K.
!
!    Input, real ( kind = 8 ) X(N+1), real Y(N+1), containing the coordinates 
!    of the nodes in the first N positions with non-constraint nodes
!    in the first LCC(1)-1 locations if 0 < NCC.  On output, updated with
!    the insertion of XK and YK in the K-th positions (node I+1 was node 
!    I before the insertion for I = K to N if K <= N)
!    unless IER /= 0 and IER /= -4.
!
!    Input/output, integer ( kind = 4 ) LIST(*), LPTR(*), LEND(N), LNEW, the data 
!    structure associated with the triangulation of nodes 1 to N.  The 
!    arrays must have sufficient length for N+1 nodes.  Refer to TRMESH.
!    On output, updated with the addition of node K unless
!    IER /= 0 and IER /= -4.
!
!    Output, integer ( kind = 4 ) IER = Error indicator:
!     0 if no errors were encountered.
!    -1 if K, IST, NCC, N, or an LCC entry is outside its valid range on input.
!    -2 if all nodes (including K) are collinear.
!     L if nodes L and K coincide for some L.
!    -3 if K lies in a constraint region.
!    -4 if an error flag is returned by SWAP implying that the triangulation
!      (geometry) was bad on input.
!
  implicit none

  logical crtri
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) ibk
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) indxcc
  integer ( kind = 4 ) io1
  integer ( kind = 4 ) io2
  integer ( kind = 4 ) ist
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kk
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lcc(*)
  integer ( kind = 4 ) lccip1
  integer ( kind = 4 ) lend(*)
  integer ( kind = 4 ) list(*)
  integer ( kind = 4 ) lnew
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lpf
  integer ( kind = 4 ) lpo1
  integer ( kind = 4 ) lptr(*)
  integer ( kind = 4 ) lstptr
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ncc
  integer ( kind = 4 ) nm1
  logical swptst
  real ( kind = 8 ) x(*)
  real ( kind = 8 ) xk
  real ( kind = 8 ) y(*)
  real ( kind = 8 ) yk

  kk = k
!
!  Test for an invalid input parameter.
!
  if ( kk < 1  .or.  ist < 1  .or.  n < ist &
      .or.  ncc < 0  .or.  n < 3 ) then
    ier = -1
    return
  end if

  lccip1 = n + 1

  do i = ncc, 1, -1
    if ( lccip1-lcc(i) < 3 ) then
      ier = -1
      return
    end if
    lccip1 = lcc(i)
  end do

  if ( lccip1 < kk ) then
    ier = -1
    return
  end if
!
!  Find a triangle (I1,I2,I3) containing K or the rightmost
!  (I1) and leftmost (I2) visible boundary nodes as viewed from node K.
!
  call trfind ( ist, xk, yk, n, x, y, list, lptr, lend, i1, i2, i3 )
!
!  Test for collinear nodes, duplicate nodes, and K lying in
!  a constraint region.
!
  if ( i1 == 0 ) then
    ier = -2
    return
  end if

  if ( i3 /= 0 ) then

    l = i1
    if ( xk == x(l)  .and.  yk == y(l) ) then
      ier = l
      return
    end if

    l = i2
    if ( xk == x(l)  .and.  yk == y(l) ) then
      ier = l
      return
    end if

    l = i3
    if ( xk == x(l)  .and.  yk == y(l) ) then
      ier = l
      return
    end if

    if ( 0 < ncc .and.  crtri(ncc,lcc,i1,i2,i3) ) then
      ier = -3
      return
    end if

  else
!
!  K is outside the convex hull of the nodes and lies in a
!  constraint region iff an exterior constraint curve is present.
!
    if ( 0 < ncc .and. indxcc(ncc,lcc,n,list,lend) /= 0 ) then
      ier = -3
      return
    end if

  end if
!
!  No errors encountered.
!
  ier = 0
  nm1 = n
  n = n + 1

  if (kk < n) then
!
!  Open a slot for K in X, Y, and LEND, and increment all
!  nodal indexes which are greater than or equal to K.
!
!  Note that LIST, LPTR, and LNEW are not yet updated with
!  either the neighbors of K or the edges terminating on K.
!
    do ibk = nm1, kk, -1
      x(ibk+1) = x(ibk)
      y(ibk+1) = y(ibk)
      lend(ibk+1) = lend(ibk)
    end do

    do i = 1, ncc
      lcc(i) = lcc(i) + 1
    end do

    l = lnew - 1

    do i = 1, l

      if ( kk <= list(i) ) then
        list(i) = list(i) + 1
      end if

      if ( list(i) <= -kk ) then
        list(i) = list(i) - 1
      end if

    end do

    if ( kk <= i1 ) then
      i1 = i1 + 1
    end if

    if ( kk <= i2 ) then
      i2 = i2 + 1
    end if

    if ( kk <= i3 ) then
      i3 = i3 + 1
    end if

  end if
!
!  Insert K into X and Y, and update LIST, LPTR, LEND, and
!  LNEW with the arcs containing node K.
!
  x(kk) = xk
  y(kk) = yk

  if ( i3 == 0 ) then
    call bdyadd ( kk, i1, i2, list, lptr, lend, lnew )
  else
    call intadd ( kk, i1, i2, i3, list, lptr, lend, lnew )
  end if
!
!  Initialize variables for optimization of the triangulation.
!
  lp = lend(kk)
  lpf = lptr(lp)
  io2 = list(lpf)
  lpo1 = lptr(lpf)
  io1 = abs ( list(lpo1) )
!
!  Begin loop:  find the node opposite K.
!
  do

    lp = lstptr ( lend(io1), io2, list, lptr )

    if ( 0 <= list(lp) ) then

      lp = lptr(lp)
      in1 = abs ( list(lp) )
!
!  Swap test:  if a swap occurs, two new arcs are
!  opposite K and must be tested.
!
      if ( .not. crtri ( ncc, lcc, io1, io2, in1 ) ) then

        if ( swptst(in1,kk,io1,io2,x,y) ) then

          call swap ( in1, kk, io1, io2, list, lptr, lend, lpo1 )

          if ( lpo1 == 0 ) then
            ier = -4
            exit
          end if
  
          io1 = in1

          cycle

        end if

      end if

    end if
!
!  No swap occurred.  Test for termination and reset IO2 and IO1.
!
    if ( lpo1 == lpf .or. list(lpo1) < 0 ) then
      exit
    end if

    io2 = io1
    lpo1 = lptr(lpo1)
    io1 = abs ( list(lpo1) )

  end do

  return
end
function areap ( x, y, nb, nodes )

!*****************************************************************************80
!
!! AREAP computes the signed area of a polygonal curve.
!
!  Discussion:
!
!    Given a sequence of NB points in the plane, this function
!    computes the signed area bounded by the closed polygonal
!    curve which passes through the points in the
!    specified order.  Each simple closed curve is positively
!    oriented (bounds positive area) if and only if the points
!    are specified in counterclockwise order.  The last point
!    of the curve is taken to be the first point specified, and
!    this point should therefore not be specified twice.
!
!    The area of a triangulation may be computed by calling
!    AREAP with values of NB and NODES determined by subroutine
!    BNODES.
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    Robert Renka,
!    Department of Computer Science,
!    University of North Texas,
!    renka@cs.unt.edu
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 751: TRIPACK, 
!    A Constrained Two-Dimensional Delaunay Triangulation Package,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 1, 1996.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X(*), Y(*), the Cartesian coordinates of a set 
!    of points.
!
!    Input, integer ( kind = 4 ) NB, the number of points in the curve.
!
!    Input, integer ( kind = 4 ) NODES(NB), the indices of the points that
!    make up the closed curve.
!
!    Output, real ( kind = 8 ) AREAP, the signed area bounded by the curve.
!
  implicit none

  integer ( kind = 4 ) nb

  real ( kind = 8 ) a
  real ( kind = 8 ) areap
  integer ( kind = 4 ) i
  integer ( kind = 4 ) nd1
  integer ( kind = 4 ) nd2
  integer ( kind = 4 ) nodes(nb)
  real ( kind = 8 ) x(*)
  real ( kind = 8 ) y(*)

  a = 0.0D+00

  if ( nb < 3 ) then
    areap = 0.0D+00
    return
  end if

  nd2 = nodes(nb)
!
!  Loop on line segments NODES(I-1) -> NODES(I), where
!  NODES(0) = NODES(NB), adding twice the signed trapezoid
!  areas (integrals of the linear interpolants) to A.
!
  do i = 1, nb
    nd1 = nd2
    nd2 = nodes(i)
    a = a + ( x(nd2) - x(nd1) ) * ( y(nd1) + y(nd2) )
  end do
!
!  A contains twice the negative signed area of the region.
!
  areap = -a / 2.0D+00

  return
end
subroutine bdyadd ( kk, i1, i2, list, lptr, lend, lnew )

!*****************************************************************************80
!
!! BDYADD adds a boundary node to a triangulation.
!
!  Discussion:
!
!    This subroutine adds a boundary node to a triangulation
!    of a set of points in the plane.  The data structure is
!    updated with the insertion of node KK, but no optimization
!    is performed.
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    Robert Renka,
!    Department of Computer Science,
!    University of North Texas,
!    renka@cs.unt.edu
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 751: TRIPACK, 
!    A Constrained Two-Dimensional Delaunay Triangulation Package,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 1, 1996.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) KK, the index of a node to be connected to the sequence
!    of all visible boundary nodes.  1 <= KK and
!    KK must not be equal to I1 or I2.
!
!    Input, integer ( kind = 4 ) I1, the first (rightmost as viewed from KK) boundary
!    node in the triangulation which is visible from
!    node KK (the line segment KK-I1 intersects no arcs.
!
!    Input, integer ( kind = 4 ) I2, the last (leftmost) boundary node which is visible
!    from node KK.  I1 and I2 may be determined by subroutine TRFIND.
!
!    Input/output, integer ( kind = 4 ) LIST(*), LPTR(*), LEND(N), LNEW.  The 
!    triangulation data structure created by TRMESH or TRMSHR.
!    On input, nodes I1 and I2 must be included in the triangulation.
!    On output, the data structure has been updated with the addition 
!    of node KK.  Node KK is connected to I1, I2, and all boundary 
!    nodes in between.
!
  implicit none

  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kk
  integer ( kind = 4 ) lend(*)
  integer ( kind = 4 ) list(*)
  integer ( kind = 4 ) lnew
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lptr(*)
  integer ( kind = 4 ) lsav
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) next
  integer ( kind = 4 ) nsav

  k = kk
  n1 = i1
  n2 = i2
!
!  Add K as the last neighbor of N1.
!
  lp = lend(n1)
  lsav = lptr(lp)
  lptr(lp) = lnew
  list(lnew) = -k
  lptr(lnew) = lsav
  lend(n1) = lnew
  lnew = lnew + 1
  next = -list(lp)
  list(lp) = next
  nsav = next
!
!  Loop on the remaining boundary nodes between N1 and N2,
!  adding K as the first neighbor.
!
  do

    lp = lend(next)

    call insert ( k, lp, list, lptr, lnew )

    if ( next == n2 ) then
      exit
    end if

    next = -list(lp)
    list(lp) = next

  end do
!
!  Add the boundary nodes between N1 and N2 as neighbors
!  of node K.
!
  lsav = lnew
  list(lnew) = n1
  lptr(lnew) = lnew + 1
  lnew = lnew + 1
  next = nsav

  do

    if ( next == n2 ) then
      exit
    end if

    list(lnew) = next
    lptr(lnew) = lnew + 1
    lnew = lnew + 1
    lp = lend(next)
    next = list(lp)

  end do

  list(lnew) = -n2
  lptr(lnew) = lsav
  lend(k) = lnew
  lnew = lnew + 1

  return
end
subroutine bnodes ( n, list, lptr, lend, nodes, nb, na, nt )

!*****************************************************************************80
!
!! BNODES returns a list of the boundary nodes.
!
!  Discussion:
!
!    Given a triangulation of N points in the plane, this
!    subroutine returns an array containing the indexes, in
!    counterclockwise order, of the nodes on the boundary of
!    the convex hull of the set of points.
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    Robert Renka,
!    Department of Computer Science,
!    University of North Texas,
!    renka@cs.unt.edu
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 751: TRIPACK, 
!    A Constrained Two-Dimensional Delaunay Triangulation Package,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 1, 1996.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of nodes in the triangulation.  3 <= N.
!
!    Input, integer ( kind = 4 ) LIST(*), LPTR(*), LEND(N), the data structure defining 
!    the triangulation.  Refer to subroutine TRMESH.
!
!    Output, integer ( kind = 4 ) NODES(NB), ordered sequence of boundary node indexes
!    in the range 1 to N.
!
!    Output, integer ( kind = 4 ) NB, the number of boundary nodes.
!
!    Output, integer ( kind = 4 ) NA, NT, the number of arcs and triangles, respectively,
!    in the triangulation.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) k
  integer ( kind = 4 ) lend(n)
  integer ( kind = 4 ) list(*)
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lptr(*)
  integer ( kind = 4 ) n0
  integer ( kind = 4 ) na
  integer ( kind = 4 ) nb
  integer ( kind = 4 ) nodes(*)
  integer ( kind = 4 ) nst
  integer ( kind = 4 ) nt
!
!  Set NST to the first boundary node encountered.
!
  nst = 1

  do

    lp = lend(nst)

    if ( list(lp) < 0 ) then
      exit
    end if

    nst = nst + 1

  end do
!
!  Initialization.
!
  nodes(1) = nst
  k = 1
  n0 = nst
!
!  Traverse the boundary in counterclockwise order.
!
  do

    lp = lend(n0)
    lp = lptr(lp)
    n0 = list(lp)

    if ( n0 == nst ) then
      exit
    end if

    k = k + 1
    nodes(k) = n0

  end do
!
!  Termination.
!
  nb = k
  nt = 2 * n - nb - 2
  na = nt + n - 1

  return
end
subroutine circum ( x1, y1, x2, y2, x3, y3, ratio, xc, yc, cr, sa, ar )

!*****************************************************************************80
!
!! CIRCUM determines the circumcenter (and more) of a triangle.
!
!  Discussion:
!
!    Given three vertices defining a triangle, this routine
!    returns the circumcenter, circumradius, signed
!    triangle area, and, optionally, the aspect ratio of the
!    triangle.
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    Robert Renka,
!    Department of Computer Science,
!    University of North Texas,
!    renka@cs.unt.edu
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 751: TRIPACK, 
!    A Constrained Two-Dimensional Delaunay Triangulation Package,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 1, 1996.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X1, Y1, X2, Y2, X3, Y3, the coordinates of
!    the vertices.
!
!    Input, logical RATIO, is TRUE if and only if the aspect ratio is 
!    to be computed.
!
!    Output, real ( kind = 8 ) XC, YC, coordinates of the circumcenter (center
!    of the circle defined by the three points) unless SA = 0, in which XC 
!    and YC are not altered.
!
!    Output, real ( kind = 8 ) CR, the circumradius (radius of the circle
!    defined by the three points) unless SA = 0 (infinite radius), in which
!    case CR is not altered.
!
!    Output, real ( kind = 8 ) SA, the signed triangle area with positive value
!    if and only if the vertices are specified in counterclockwise order:  
!    (X3,Y3) is strictly to the left of the directed line from (X1,Y1)
!    toward (X2,Y2).
!
!    Output, real ( kind = 8 ) AR, the aspect ratio r/CR, where r is the 
!    radius of the inscribed circle, unless RATIO = FALSE, in which case AR
!    is not altered.  AR is in the range 0 to 0.5, with value 0 iff SA = 0 and
!    value 0.5 iff the vertices define an equilateral triangle.
!
  implicit none

  real ( kind = 8 ) ar
  real ( kind = 8 ) cr
  real ( kind = 8 ) ds(3)
  real ( kind = 8 ) fx
  real ( kind = 8 ) fy
  logical ratio
  real ( kind = 8 ) sa
  real ( kind = 8 ) u(3)
  real ( kind = 8 ) v(3)
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ) xc
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) y3
  real ( kind = 8 ) yc
!
!  Set U(K) and V(K) to the x and y components, respectively,
!  of the directed edge opposite vertex K.
!
  u(1) = x3 - x2
  u(2) = x1 - x3
  u(3) = x2 - x1
  v(1) = y3 - y2
  v(2) = y1 - y3
  v(3) = y2 - y1
!
!  Set SA to the signed triangle area.
!
  sa = ( u(1) * v(2) - u(2) * v(1) ) / 2.0D+00

  if ( sa == 0.0D+00 ) then
    if ( ratio ) then
      ar = 0.0D+00
    end if
    return
  end if
!
!  Set DS(K) to the squared distance from the origin to vertex K.
!
  ds(1) = x1 * x1 + y1 * y1
  ds(2) = x2 * x2 + y2 * y2
  ds(3) = x3 * x3 + y3 * y3
!
!  Compute factors of XC and YC.
!
  fx = - dot_product ( ds(1:3), v(1:3) )
  fy =   dot_product ( ds(1:3), u(1:3) )

  xc = fx / ( 4.0D+00 * sa )
  yc = fy / ( 4.0D+00 * sa )
  cr = sqrt ( ( xc - x1 )**2 + ( yc - y1 )**2 )

  if ( .not. ratio ) then
    return
  end if
!
!  Compute the squared edge lengths and aspect ratio.
!
  ds(1:3) = u(1:3)**2 + v(1:3)**2

  ar = 2.0D+00 * abs ( sa ) / &
       ( ( sqrt ( ds(1) ) + sqrt ( ds(2) ) + sqrt ( ds(3) ) ) * cr )

  return
end
function crtri ( ncc, lcc, i1, i2, i3 )

!*****************************************************************************80
!
!! CRTRI determines if a triangle lies in a constraint region.
!
!  Discussion:
!
!    This function returns TRUE if and only if triangle (I1,
!    I2,I3) lies in a constraint region.
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    Robert Renka,
!    Department of Computer Science,
!    University of North Texas,
!    renka@cs.unt.edu
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 751: TRIPACK, 
!    A Constrained Two-Dimensional Delaunay Triangulation Package,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 1, 1996.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NCC, LCC(*), onstraint data structure.  Refer to ADDCST.
!
!    Input, integer ( kind = 4 ) I1, I2, I3, Nodal indexes of the counterclockwise-
!    ordered vertices of a triangle.
!
!    Output, logical CRTRI, is TRUE if and only if (I1,I2,I3) is a 
!    constraint region triangle.
!
  implicit none

  logical crtri
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  integer ( kind = 4 ) lcc(*)
  integer ( kind = 4 ) ncc

  imax = max ( i1, i2, i3 )
!
!  Find the index I of the constraint containing IMAX.
!
  i = ncc + 1

  do

    i = i - 1

    if ( i <= 0 ) then
      crtri = .false.
      return
    end if

    if ( lcc(i) <= imax ) then
      exit
    end if

  end do

  imin = min ( i1, i2, i3 )
!
!  P lies in a constraint region iff I1, I2, and I3 are nodes
!  of the same constraint (LCC(I) <= IMIN), and (IMIN,IMAX)
!  is (I1,I3), (I2,I1), or (I3,I2).
!
  crtri = lcc(i) <= imin .and. ( &
    ( imin == i1 .and. imax == i3 ) .or.  &
    ( imin == i2 .and. imax == i1 ) .or.  &
    ( imin == i3 .and. imax == i2 ) )

  return
end
subroutine delarc ( n, io1, io2, list, lptr, lend, lnew, ier )

!*****************************************************************************80
!
!! DELARC deletes a boundary arc from a triangulation.
!
!  Discussion:
!
!    This subroutine deletes a boundary arc from a triangula-
!    tion.  It may be used to remove a null triangle from the
!    convex hull boundary.  Note, however, that if the union of
!    triangles is rendered nonconvex, Subroutines DELNOD, EDGE,
!    and TRFIND may fail.  Thus, Subroutines ADDCST, ADDNOD,
!    DELNOD, EDGE, and NEARND should not be called following
!    an arc deletion.
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    Robert Renka,
!    Department of Computer Science,
!    University of North Texas,
!    renka@cs.unt.edu
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 751: TRIPACK, 
!    A Constrained Two-Dimensional Delaunay Triangulation Package,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 1, 1996.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of nodes in the triangulation.  4 <= N.
!
!    Input, integer ( kind = 4 ) IO1, IO2, the indexes (in the range 1 to N) of a pair of
!    adjacent boundary nodes defining the arc to be removed.
!
!    Input/output, integer ( kind = 4 ) LIST(*), LPTR(*), LEND(N), LNEW, the
!    triangulation data structure created by TRMESH or TRMSHR.
!    On output, updated with the removal of arc IO1-IO2 unless 0 < IER.
!
!    Output, integer ( kind = 4 ) IER = Error indicator:
!    0, if no errors were encountered.
!    1, if N, IO1, or IO2 is outside its valid range, or IO1 = IO2.
!    2, if IO1-IO2 is not a boundary arc.
!    3, if the node opposite IO1-IO2 is already a boundary node, and 
!      thus IO1 or IO2 has only two neighbors or a deletion would result 
!      in two triangulations sharing a single node.
!    4, if one of the nodes is a neighbor of the other, but not vice 
!      versa, implying an invalid triangulation data structure.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) ier
  integer ( kind = 4 ) io1
  integer ( kind = 4 ) io2
  integer ( kind = 4 ) lend(n)
  integer ( kind = 4 ) list(*)
  integer ( kind = 4 ) lnew
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lph
  integer ( kind = 4 ) lpl
  integer ( kind = 4 ) lptr(*)
  integer ( kind = 4 ) lstptr
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n3

  n1 = io1
  n2 = io2
!
!  Test for errors, and set N1->N2 to the directed boundary
!  edge associated with IO1-IO2:  (N1,N2,N3) is a triangle
!  for some N3.
!
  if ( n < 4  .or.  n1 < 1  .or.  n < n1 .or. &
      n2 < 1  .or.  n < n2 .or.  n1 == n2 ) then
    ier = 1
    return
  end if

  lpl = lend(n2)

  if ( -list(lpl) /= n1 ) then

    n1 = n2
    n2 = io1
    lpl = lend(n2)

    if ( -list(lpl) /= n1 ) then
      ier = 2
      return
    end if

  end if
!
!  Set N3 to the node opposite N1->N2 (the second neighbor
!  of N1), and test for error 3 (N3 already a boundary node).
!
  lpl = lend(n1)
  lp = lptr(lpl)
  lp = lptr(lp)
  n3 = abs ( list(lp) )
  lpl = lend(n3)

  if ( list(lpl) <= 0 ) then
    ier = 3
    return
  end if
!
!  Delete N2 as a neighbor of N1, making N3 the first
!  neighbor, and test for error 4 (N2 not a neighbor
!  of N1).  Note that previously computed pointers may
!  no longer be valid following the call to DELNB.
!
  call delnb ( n1, n2, n, list, lptr, lend, lnew, lph )

  if ( lph < 0 ) then
    ier = 4
    return
  end if
!
!  Delete N1 as a neighbor of N2, making N3 the new last neighbor.
!
  call delnb ( n2, n1, n, list, lptr, lend, lnew, lph )
!
!  Make N3 a boundary node with first neighbor N2 and last neighbor N1.
!
  lp = lstptr ( lend(n3), n1, list, lptr )
  lend(n3) = lp
  list(lp) = -n1
!
!  No errors encountered.
!
  ier = 0

  return
end
subroutine delnb ( n0, nb, n, list, lptr, lend, lnew, lph )

!*****************************************************************************80
!
!! DELNB deletes a neighbor from an adjacency list.
!
!  Discussion:
!
!    This subroutine deletes a neighbor NB from the adjacency
!    list of node N0 (but N0 is not deleted from the adjacency
!    list of NB) and, if NB is a boundary node, makes N0 a
!    boundary node.  For pointer (LIST index) LPH to NB as a
!    neighbor of N0, the empty LIST,LPTR location LPH is filled
!    in with the values at LNEW-1, pointer LNEW-1 (in LPTR and
!    possibly in LEND) is changed to LPH, and LNEW is decremented
!    This requires a search of LEND and LPTR entailing an
!    expected operation count of O(N).
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    Robert Renka,
!    Department of Computer Science,
!    University of North Texas,
!    renka@cs.unt.edu
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 751: TRIPACK, 
!    A Constrained Two-Dimensional Delaunay Triangulation Package,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 1, 1996.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N0, NB, indexes, in the range 1 to N, of a pair of
!    nodes such that NB is a neighbor of N0.
!    (N0 need not be a neighbor of NB.)
!
!    Input, integer ( kind = 4 ) N, the number of nodes in the triangulation.  3 <= N.
!
!    Input/output, integer ( kind = 4 ) LIST(*), LPTR(*), LEND(N), LNEW, the data 
!    structure defining the triangulation.  On output, updated with
!    the removal of NB from the adjacency list of N0 unless LPH < 0.
!
!    Output, integer ( kind = 4 ) LPH, list pointer to the hole (NB as a neighbor of
!    N0) filled in by the values at LNEW-1 or error indicator:
!    >0, if no errors were encountered.
!    -1, if N0, NB, or N is outside its valid range.
!    -2, if NB is not a neighbor of N0.
!
!  Local parameters:
!
!    I =   DO-loop index
!    LNW = LNEW-1 (output value of LNEW)
!    LP =  LIST pointer of the last neighbor of NB
!    LPB = Pointer to NB as a neighbor of N0
!    LPL = Pointer to the last neighbor of N0
!    LPP = Pointer to the neighbor of N0 that precedes NB
!    NN =  Local copy of N
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) lend(n)
  integer ( kind = 4 ) list(*)
  integer ( kind = 4 ) lnew
  integer ( kind = 4 ) lnw
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lpb
  integer ( kind = 4 ) lph
  integer ( kind = 4 ) lpl
  integer ( kind = 4 ) lpp
  integer ( kind = 4 ) lptr(*)
  integer ( kind = 4 ) n0
  integer ( kind = 4 ) nb
  integer ( kind = 4 ) nn

  nn = n
!
!  Test for error 1.
!
  if ( n0 < 1  .or. nn < n0 .or.  nb < 1  .or. &
      nn < nb .or. nn < 3 ) then
    lph = -1
    return
  end if
!
!  Find pointers to neighbors of N0:
!
!  LPL points to the last neighbor,
!  LPP points to the neighbor NP preceding NB, and
!  LPB points to NB.
!
  lpl = lend(n0)
  lpp = lpl
  lpb = lptr(lpp)

  do

    if ( list(lpb) == nb ) then
      go to 2
    end if

    lpp = lpb
    lpb = lptr(lpp)

    if ( lpb == lpl ) then
      exit
    end if

  end do
!
!  Test for error 2 (NB not found).
!
  if ( abs ( list(lpb) ) /= nb ) then
    lph = -2
    return
  end if
!
!  NB is the last neighbor of N0.  Make NP the new last
!  neighbor and, if NB is a boundary node, then make N0
!  a boundary node.
!
  lend(n0) = lpp
  lp = lend(nb)
  if ( list(lp) < 0 ) then
    list(lpp) = -list(lpp)
  end if
  go to 3
!
!  NB is not the last neighbor of N0.  If NB is a boundary
!  node and N0 is not, then make N0 a boundary node with
!  last neighbor NP.
!
2 continue

  lp = lend(nb)

  if ( list(lp) < 0  .and.  0 < list(lpl) ) then
    lend(n0) = lpp
    list(lpp) = -list(lpp)
  end if
!
!  Update LPTR so that the neighbor following NB now fol-
!  lows NP, and fill in the hole at location LPB.
!
3 continue

  lptr(lpp) = lptr(lpb)
  lnw = lnew - 1
  list(lpb) = list(lnw)
  lptr(lpb) = lptr(lnw)

  do i = nn, 1, -1
    if (lend(i) == lnw) then
      lend(i) = lpb
      exit
    end if
  end do

  do i = 1, lnw-1
    if (lptr(i) == lnw) then
      lptr(i) = lpb
    end if
  end do
!
!  No errors encountered.
!
  lnew = lnw
  lph = lpb

  return
end
subroutine delnod ( k, ncc, lcc, n, x, y, list, lptr, lend, lnew, lwk, iwk, &
  ier )

!*****************************************************************************80
!
!! DELNOD deletes a node from a triangulation.
!
!  Discussion:
!
!    This subroutine deletes node K (along with all arcs
!    incident on node K) from a triangulation of N nodes in the
!    plane, and inserts arcs as necessary to produce a triangu-
!    lation of the remaining N-1 nodes.  If a Delaunay triangu-
!    lation is input, a Delaunay triangulation will result, and
!    thus, DELNOD reverses the effect of a call to subroutine
!    ADDNOD.
!
!    Note that a constraint node cannot be deleted by this
!    routine.  In order to delete a constraint node, it is
!    necessary to call this routine with NCC = 0, decrement the
!    appropriate LCC entries (LCC(I) such that K < LCC(I) ), and
!    then create (or restore) the constraints by a call to sub-
!    routine ADDCST.
!
!    Note that the deletion may result in all remaining nodes
!    being collinear.  This situation is not flagged.
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    Robert Renka,
!    Department of Computer Science,
!    University of North Texas,
!    renka@cs.unt.edu
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 751: TRIPACK, 
!    A Constrained Two-Dimensional Delaunay Triangulation Package,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 1, 1996.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) K, the index (for X and Y) of the node to be deleted.
!    1 <= K < LCC(1).  (K <= N if NCC=0).
!
!    Input, integer ( kind = 4 ) NCC, the number of constraint curves.  0 <= NCC.
!
!    Input/output, integer ( kind = 4 ) LCC(*), list of constraint curve starting indexes 
!    (or dummy array of length 1 if NCC = 0).  Refer to subroutine ADDCST.
!    On output, decremented by 1 to reflect the deletion of K unless NCC = 0 
!    or 1 <= IER <= 4.
!
!    Input/output, integer ( kind = 4 ) N, the number of nodes in the triangulation.
!    4 <= N.  Note that N will be decremented following the deletion.
!
!    Input/output, real ( kind = 8 ) X(N), Y(N), the coordinates of the nodes
!    with non-constraint nodes in the first LCC(1)-1 locations if 0 < NCC.
!    On output, updated arrays of length N-1 containing nodal coordinates 
!    (with elements K+1,...,N shifted a position and thus overwriting 
!    element K) unless 1 <= IER <= 4.  (N here denotes the input value.)
!
!    Input/output, integer ( kind = 4 ) LIST(*), LPTR(*), LEND(N), LNEW, the data 
!    structure defining the triangulation.  Refer to subroutine TRMESH.
!    On output, updated to reflect the deletion unless IER /= 0.  Note
!    that the data structure may have been altered if 3 <= IER.
!
!    Input/output, integer ( kind = 4 ) LWK.  On input, the number of columns reserved 
!    for IWK.  LWK must be at least NNB-3, where NNB is the number of
!    neighbors of node K, including an extra pseudo-node if K is a 
!    boundary node.  On output, the number of IWK columns required unless 
!    IER = 1 or IER = 3.
!
!    Output, integer ( kind = 4 ) IWK(2,LWK), indexes of the endpoints of the new arcs added
!    unless LWK = 0 or 1 <= IER <= 4.  (Arcs are associated with columns, or
!    pairs of adjacent elements if IWK is declared as a singly-subscripted
!    array.)
!
!    Output, integer ( kind = 4 ) IER = Error indicator:
!    0, if no errors were encountered.
!    1, if K, NCC, N, or an LCC entry is outside its valid range or LWK < 0 
!      on input.
!    2, if more space is required in IWK.  Refer to LWK.
!    3, if the triangulation data structure is invalid on input.
!    4, if K is an interior node with 4 or more neighbors, and the number of
!      neighbors could not be reduced to 3 by swaps.  This could be caused by
!      floating point errors with collinear nodes or by an invalid data
!      structure.
!    5, if an error flag was returned by OPTIM.  An error message is written
!      to the standard output unit in this event.
!
  implicit none

  logical bdry
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) iwk(2,*)
  integer ( kind = 4 ) iwl
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) lcc(*)
  integer ( kind = 4 ) lccip1
  logical left
  integer ( kind = 4 ) lend(*)
  integer ( kind = 4 ) list(*)
  integer ( kind = 4 ) lnew
  integer ( kind = 4 ) lnw
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lp21
  integer ( kind = 4 ) lpf
  integer ( kind = 4 ) lph
  integer ( kind = 4 ) lpl
  integer ( kind = 4 ) lpl2
  integer ( kind = 4 ) lpn
  integer ( kind = 4 ) lptr(*)
  integer ( kind = 4 ) lstptr
  integer ( kind = 4 ) lwk
  integer ( kind = 4 ) lwkl
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) nbcnt
  integer ( kind = 4 ) ncc
  integer ( kind = 4 ) nfrst
  integer ( kind = 4 ) nit
  integer ( kind = 4 ) nl
  integer ( kind = 4 ) nn
  integer ( kind = 4 ) nnb
  integer ( kind = 4 ) nr
  real ( kind = 8 ) x(*)
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) xl
  real ( kind = 8 ) xr
  real ( kind = 8 ) y(*)
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) yl
  real ( kind = 8 ) yr
!
!  Set N1 to K and NNB to the number of neighbors of N1 (plus
!  one if N1 is a boundary node), and test for errors.  LPF
!  and LPL are LIST indexes of the first and last neighbors
!  of N1, IWL is the number of IWK columns containing arcs,
!  and BDRY is TRUE iff N1 is a boundary node.
!
  ier = 0
  n1 = k
  nn = n

  if ( ncc < 0  .or.  n1 < 1  .or.  nn < 4  .or. lwk < 0 ) then
    ier = 1
    return
  end if

  lccip1 = nn + 1

  do i = ncc, 1, -1

    if ( lccip1 - lcc(i) < 3 ) then
      ier = 1
      return
    end if

    lccip1 = lcc(i)

  end do

  if ( lccip1 <= n1 ) then
    ier = 1
    return
  end if

  lpl = lend(n1)
  lpf = lptr(lpl)
  nnb = nbcnt ( lpl, lptr )
  bdry = list(lpl) < 0

  if ( bdry ) then
    nnb = nnb + 1
  end if

  if ( nnb < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DELNOD - Fatal error!'
    write ( *, '(a)' ) '  NNB < 3.'
    write ( *, '(a,i6)' ) '  NNB = ', nnb
    ier = 3
    return
  end if

  lwkl = lwk
  lwk = nnb - 3

  if ( lwkl < lwk ) then
    ier = 2
    return
  end if

  iwl = 0

  if ( nnb == 3 ) then
    go to 5
  end if
!
!  Initialize for loop on arcs N1-N2 for neighbors N2 of N1,
!  beginning with the second neighbor.  NR and NL are the
!  neighbors preceding and following N2, respectively, and
!  LP indexes NL.  The loop is exited when all possible
!  swaps have been applied to arcs incident on N1.  If N1
!  is interior, the number of neighbors will be reduced
!  to 3.
!
  x1 = x(n1)
  y1 = y(n1)
  nfrst = list(lpf)
  nr = nfrst
  xr = x(nr)
  yr = y(nr)
  lp = lptr(lpf)
  n2 = list(lp)
  x2 = x(n2)
  y2 = y(n2)
  lp = lptr(lp)
!
!  Top of loop:  set NL to the neighbor following N2.
!
2 continue

  nl = abs ( list(lp) )

  if ( nl == nfrst .and. bdry ) then
    go to 5
  end if

  xl = x(nl)
  yl = y(nl)
!
!  Test for a convex quadrilateral.  To avoid an incorrect
!  test caused by collinearity, use the fact that if N1
!  is a boundary node, then N1 LEFT NR->NL and if N2 is
!  a boundary node, then N2 LEFT NL->NR.
!
  lpl2 = lend(n2)

  if ( (bdry  .or.  left(xr,yr,xl,yl,x1,y1))  .and. &
       (list(lpl2) < 0  .or. &
        left(xl,yl,xr,yr,x2,y2)) ) then
    go to 3
  end if
!
!  Nonconvex quadrilateral -- no swap is possible.
!
  nr = n2
  xr = x2
  yr = y2
  go to 4
!
!  The quadrilateral defined by adjacent triangles
!  (N1,N2,NL) and (N2,N1,NR) is convex.  Swap in
!  NL-NR and store it in IWK.  Indexes larger than N1
!  must be decremented since N1 will be deleted from
!  X and Y.
!
3 continue

  call swap ( nl, nr, n1, n2, list, lptr, lend, lp21 )
  iwl = iwl + 1

  if ( nl <= n1 ) then
    iwk(1,iwl) = nl
  else
    iwk(1,iwl) = nl - 1
  end if

  if ( nr <= n1 ) then
    iwk(2,iwl) = nr
  else
    iwk(2,iwl) = nr - 1
  end if
!
!  Recompute the LIST indexes LPL,LP and decrement NNB.
!
  lpl = lend(n1)
  nnb = nnb - 1

  if ( nnb == 3 ) then
    go to 5
  end if

  lp = lstptr ( lpl, nl, list, lptr )

  if ( nr == nfrst ) then
    go to 4
  end if
!
!  NR is not the first neighbor of N1.
!  Back up and test N1-NR for a swap again:  Set N2 to
!  NR and NR to the previous neighbor of N1 -- the
!  neighbor of NR which follows N1.  LP21 points to NL
!  as a neighbor of NR.
!
  n2 = nr
  x2 = xr
  y2 = yr
  lp21 = lptr(lp21)
  lp21 = lptr(lp21)
  nr = abs ( list(lp21) )
  xr = x(nr)
  yr = y(nr)
  go to 2
!
!  Bottom of loop -- test for invalid termination.
!
4 continue

  if ( n2 == nfrst ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DELNOD - Fatal error!'
    write ( *, '(a)' ) '  N2 = NFRST.'
    ier = 4
    return
  end if

  n2 = nl
  x2 = xl
  y2 = yl
  lp = lptr(lp)
  go to 2
!
!  Delete N1 from the adjacency list of N2 for all neighbors
!  N2 of N1.  LPL points to the last neighbor of N1.
!  LNEW is stored in local variable LNW.
!
5 continue

  lp = lpl
  lnw = lnew
!
!  Loop on neighbors N2 of N1, beginning with the first.
!
6 continue

    lp = lptr(lp)
    n2 = abs ( list(lp) )

    call delnb ( n2, n1, n, list, lptr, lend, lnw, lph )

    if ( lph < 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DELNOD - Fatal error!'
      write ( *, '(a)' ) '  LPH < 0.'
      ier = 3
      return
    end if
!
!  LP and LPL may require alteration.
!
    if ( lpl == lnw ) then
      lpl = lph
    end if

    if ( lp == lnw ) then
      lp = lph
    end if

    if ( lp /= lpl ) then
      go to 6
    end if
!
!  Delete N1 from X, Y, and LEND, and remove its adjacency
!  list from LIST and LPTR.  LIST entries (nodal indexes)
!  which are larger than N1 must be decremented.
!
  nn = nn - 1

  if ( nn < n1 ) then
    go to 9
  end if

  do i = n1,nn
    x(i) = x(i+1)
    y(i) = y(i+1)
    lend(i) = lend(i+1)
  end do

  do i = 1,lnw-1

    if (list(i) > n1) then
      list(i) = list(i) - 1
    end if

    if (list(i) < -n1) then
      list(i) = list(i) + 1
    end if

  end do
!
!  For LPN = first to last neighbors of N1, delete the
!  preceding neighbor (indexed by LP).
!
!  Each empty LIST,LPTR location LP is filled in with the
!  values at LNW-1, and LNW is decremented.  All pointers
!  (including those in LPTR and LEND) with value LNW-1
!  must be changed to LP.
!
!  LPL points to the last neighbor of N1.
!
9 continue

  if ( bdry ) then
    nnb = nnb - 1
  end if

  lpn = lpl

  do j = 1, nnb

    lnw = lnw - 1
    lp = lpn
    lpn = lptr(lp)
    list(lp) = list(lnw)
    lptr(lp) = lptr(lnw)

    if ( lptr(lpn) == lnw ) then
      lptr(lpn) = lp
    end if

    if ( lpn == lnw ) then
      lpn = lp
    end if

    do i = nn,1,-1
      if (lend(i) == lnw) then
        lend(i) = lp
        exit
      end if
    end do

    do i = lnw-1, 1, -1
      if ( lptr(i) == lnw ) then
        lptr(i) = lp
      end if
    end do

  end do
!
!  Decrement LCC entries.
!
  lcc(1:ncc) = lcc(1:ncc) - 1
!
!  Update N and LNEW, and optimize the patch of triangles
!  containing K (on input) by applying swaps to the arcs in IWK.
!
  n = nn
  lnew = lnw

  if ( 0 < iwl ) then
    nit = 4 * iwl
    call optim ( x, y, iwl, list, lptr, lend, nit, iwk, ierr )
    if ( ierr /= 0 ) then
      ier = 5
      write (*,100) nit, ierr
  100 format (//5x,'*** error in optim:  nit = ',i4, &
          ', ier = ',i1,' ***'/)
      return
    end if
  end if
!
!  Successful termination.
!
  ier = 0

  return
end
subroutine edge ( in1, in2, x, y, lwk, iwk, list, lptr, lend, ier )

!*****************************************************************************80
!
!! EDGE swaps arcs to force two nodes to be adjacent.
!
!  Discussion:
!
!    Given a triangulation of N nodes and a pair of nodal
!    indexes IN1 and IN2, this routine swaps arcs as necessary
!    to force IN1 and IN2 to be adjacent.  Only arcs which
!    intersect IN1-IN2 are swapped out.  If a Delaunay triangu-
!    lation is input, the resulting triangulation is as close
!    as possible to a Delaunay triangulation in the sense that
!    all arcs other than IN1-IN2 are locally optimal.
!
!    A sequence of calls to EDGE may be used to force the
!    presence of a set of edges defining the boundary of a non-
!    convex and/or multiply connected region (refer to Subrou-
!    tine ADDCST), or to introduce barriers into the triangula-
!    tion.  Note that Subroutine GETNP will not necessarily
!    return closest nodes if the triangulation has been con-
!    strained by a call to EDGE.  However, this is appropriate
!    in some applications, such as triangle-based interpolation
!    on a nonconvex domain.
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    Robert Renka,
!    Department of Computer Science,
!    University of North Texas,
!    renka@cs.unt.edu
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 751: TRIPACK, 
!    A Constrained Two-Dimensional Delaunay Triangulation Package,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 1, 1996.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IN1, IN2, indexes (of X and Y) in the range 1 to N
!    defining a pair of nodes to be connected by an arc.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the coordinates of the nodes.
!
!    Input/output, integer ( kind = 4 ) LWK.  On input, the number of columns reserved
!    for IWK.  This must be at least NI, the number of arcs which intersect
!    IN1-IN2.  (NI is bounded by N-3.)  On output, the number of arcs which
!    intersect IN1-IN2 (but not more than the input value of LWK) unless
!    IER = 1 or IER = 3.  LWK = 0 if and only if IN1 and IN2 were adjacent 
!    (or LWK=0) on input.
!
!    Output, integer ( kind = 4 ) IWK(2*LWK), the indexes of the endpoints of the new 
!    arcs other than IN1-IN2 unless IER > 0 or LWK = 0.  New arcs to the 
!    left of IN2-IN1 are stored in the first K-1 columns (left portion of 
!    IWK), column K contains zeros, and new arcs to the right of IN2-IN1
!    occupy columns K+1,...,LWK.  (K can be determined by searching IWK 
!    for the zeros.)
!
!    Input, integer ( kind = 4 ) LIST(*), LPTR(*), LEND(N), the data structure defining 
!    the triangulation.  Refer to subroutine TRMESH.  On output, updated 
!    if necessary to reflect the presence of an arc connecting IN1 and 
!    IN2 unless IER /= 0.  The data structure has been altered if IER = 4.
!
!    Output, integer ( kind = 4 ) IER = Error indicator:
!    0, if no errors were encountered.
!    1, if IN1 < 1, IN2 .LT. 1, IN1 = IN2, or LWK < 0 on input.
!    2, if more space is required in IWK.
!    3, if IN1 and IN2 could not be connected due to either an invalid 
!      data structure or collinear nodes (and floating point error).
!    4, if an error flag was returned by OPTIM.
!
!  Local parameters:
!
!    DX,DY =   Components of arc N1-N2.
!
!    I =       DO-loop index and column index for IWK
!    IERR =    Error flag returned by Subroutine OPTIM
!    IWC =     IWK index between IWF and IWL -- NL->NR is
!              stored in IWK(1,IWC)->IWK(2,IWC)
!    IWCP1 =   IWC + 1
!    IWEND =   Input or output value of LWK
!    IWF =     IWK (column) index of the first (leftmost) arc
!              which intersects IN1->IN2
!    IWL =     IWK (column) index of the last (rightmost) are
!              which intersects IN1->IN2
!    LFT =     Flag used to determine if a swap results in the
!              new arc intersecting IN1-IN2 -- LFT = 0 iff
!              N0 = IN1, LFT = -1 implies N0 LEFT IN1->IN2,
!              and LFT = 1 implies N0 LEFT IN2->IN1
!    LP21 =    Unused parameter returned by SWAP
!    LP =      List pointer (index) for LIST and LPTR
!    LPL =     Pointer to the last neighbor of IN1 or NL
!    N0 =      Neighbor of N1 or node opposite NR->NL
!    N1,N2 =   Local copies of IN1 and IN2
!    N1FRST =  First neighbor of IN1
!    N1LST =   (Signed) last neighbor of IN1
!    NEXT =    Node opposite NL->NR
!    NIT =     Flag or number of iterations employed by OPTIM
!    NL,NR =   Endpoints of an arc which intersects IN1-IN2
!              with NL LEFT IN1->IN2
!    X0,Y0 =   Coordinates of N0
!    X1,Y1 =   Coordinates of IN1
!    X2,Y2 =   Coordinates of IN2
!
  implicit none

  real ( kind = 8 ) dx
  real ( kind = 8 ) dy
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) iwc
  integer ( kind = 4 ) iwcp1
  integer ( kind = 4 ) iwend
  integer ( kind = 4 ) iwf
  integer ( kind = 4 ) iwk(2,*)
  integer ( kind = 4 ) iwl
  logical left
  integer ( kind = 4 ) lend(*)
  integer ( kind = 4 ) lft
  integer ( kind = 4 ) list(*)
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lp21
  integer ( kind = 4 ) lpl
  integer ( kind = 4 ) lptr(*)
  integer ( kind = 4 ) lwk
  integer ( kind = 4 ) n0
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n1frst
  integer ( kind = 4 ) n1lst
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) next
  integer ( kind = 4 ) nit
  integer ( kind = 4 ) nl
  integer ( kind = 4 ) nr
  real ( kind = 8 ) x(*)
  real ( kind = 8 ) x0
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) y(*)
  real ( kind = 8 ) y0
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
!
!  Store IN1, IN2, and LWK in local variables and test for errors.
!
  n1 = in1
  n2 = in2
  iwend = lwk

  if ( n1 < 1  .or.  n2 < 1  .or.  n1 == n2  .or. iwend < 0 ) then
    ier = 1
    return
  end if
!
!  Test for N2 as a neighbor of N1.  LPL points to the last neighbor of N1.
!
  lpl = lend(n1)
  n0 = abs ( list(lpl) )
  lp = lpl

  do
!
!  IN1 and IN2 were adjacent on input.
!
    if ( n0 == n2 ) then
      ier = 0
      return
    end if

    lp = lptr(lp)
    n0 = list(lp)

    if ( lp == lpl ) then
      exit
    end if

  end do
!
!  Initialize parameters.
!
  iwl = 0
  nit = 0
!
!  Store the coordinates of N1 and N2.
!
2 continue

  x1 = x(n1)
  y1 = y(n1)
  x2 = x(n2)
  y2 = y(n2)
!
!  Set NR and NL to adjacent neighbors of N1 such that
!  NR LEFT N2->N1 and NL LEFT N1->N2,
!  (NR Forward N1->N2 or NL Forward N1->N2), and
!  (NR Forward N2->N1 or NL Forward N2->N1).
!
!  Initialization:  Set N1FRST and N1LST to the first and
!  (signed) last neighbors of N1, respectively, and
!  initialize NL to N1FRST.
!
  lpl = lend(n1)
  n1lst = list(lpl)
  lp = lptr(lpl)
  n1frst = list(lp)
  nl = n1frst

  if ( n1lst < 0 ) then
    go to 4
  end if
!
!  N1 is an interior node.  Set NL to the first candidate
!  for NR (NL LEFT N2->N1).
!
3   continue

    if ( left(x2,y2,x1,y1,x(nl),y(nl)) ) then
      go to 4
    end if

    lp = lptr(lp)
    nl = list(lp)

    if ( nl /= n1frst ) then
      go to 3
    end if
!
!  All neighbors of N1 are strictly left of N1->N2.
!
  go to 5
!
!  NL = LIST(LP) LEFT N2->N1.  Set NR to NL and NL to the
!  following neighbor of N1.
!
4   continue

    nr = nl
    lp = lptr(lp)
    nl = abs ( list(lp) ) 
!
!  NL LEFT N1->N2 and NR LEFT N2->N1.  The Forward tests
!  are employed to avoid an error associated with
!  collinear nodes.
!
    if ( left(x1,y1,x2,y2,x(nl),y(nl)) ) then

      dx = x2-x1
      dy = y2-y1
      if ((dx*(x(nl)-x1)+dy*(y(nl)-y1) >= 0.  .or. &
           dx*(x(nr)-x1)+dy*(y(nr)-y1) >= 0.)  .and. &
          (dx*(x(nl)-x2)+dy*(y(nl)-y2) <= 0.  .or. &
           dx*(x(nr)-x2)+dy*(y(nr)-y2) <= 0.)) go to 6
!
!  NL-NR does not intersect N1-N2.  However, there is
!  another candidate for the first arc if NL lies on
!  the line N1-N2.
!
      if ( .not. left(x2,y2,x1,y1,x(nl),y(nl)) ) go to 5
    end if
!
!  Bottom of loop.
!
    if ( nl /= n1frst ) then
      go to 4
    end if
!
!  Either the triangulation is invalid or N1-N2 lies on the
!  convex hull boundary and an edge NR->NL (opposite N1 and
!  intersecting N1-N2) was not found due to floating point
!  error.  Try interchanging N1 and N2 -- NIT > 0 iff this
!  has already been done.
!
5 continue

  if ( 0 < nit ) then
    go to 33
  end if

  nit = 1
  n1 = n2
  n2 = in1
  go to 2
!
!  Store the ordered sequence of intersecting edges NL->NR in
!  IWK(1,IWL)->IWK(2,IWL).
!
6 continue

  iwl = iwl + 1

  if (iwl > iwend) then
    ier = 2
    return
  end if

  iwk(1,iwl) = nl
  iwk(2,iwl) = nr
!
!  Set NEXT to the neighbor of NL which follows NR.
!
  lpl = lend(nl)
  lp = lptr(lpl)
!
!  Find NR as a neighbor of NL.  The search begins with
!  the first neighbor.
!
7   continue

    if (list(lp) == nr) go to 8
    lp = lptr(lp)
    if (lp /= lpl) go to 7
!
!  NR must be the last neighbor, and NL->NR cannot be a
!  boundary edge.
!
  if (list(lp) /= nr) then
    go to 33
  end if
!
!  Set NEXT to the neighbor following NR, and test for
!  termination of the store loop.
!
8 continue

  lp = lptr(lp)
  next = abs ( list(lp) )

  if (next == n2) then
    go to 9
  end if
!
!  Set NL or NR to NEXT.
!
  if ( left(x1,y1,x2,y2,x(next),y(next)) ) then
    nl = next
  else
    nr = next
  end if

  go to 6
!
!  IWL is the number of arcs which intersect N1-N2.
!  Store LWK.
!
9 continue

  lwk = iwl
  iwend = iwl
!
!  Initialize for edge swapping loop -- all possible swaps
!  are applied (even if the new arc again intersects
!  N1-N2), arcs to the left of N1->N2 are stored in the
!  left portion of IWK, and arcs to the right are stored in
!  the right portion.  IWF and IWL index the first and last
!  intersecting arcs.
!
  iwf = 1
!
!  Top of loop -- set N0 to N1 and NL->NR to the first edge.
!  IWC points to the arc currently being processed.  LFT
!  <= 0 iff N0 LEFT N1->N2.
!
10 continue

  lft = 0
  n0 = n1
  x0 = x1
  y0 = y1
  nl = iwk(1,iwf)
  nr = iwk(2,iwf)
  iwc = iwf
!
!  Set NEXT to the node opposite NL->NR unless IWC is the last arc.
!
11 continue

  if (iwc == iwl) go to 21
  iwcp1 = iwc + 1
  next = iwk(1,iwcp1)
  if (next /= nl) go to 16
  next = iwk(2,iwcp1)
!
!  NEXT RIGHT N1->N2 and IWC < IWL.  Test for a possible swap.
!
  if ( .not. left(x0,y0,x(nr),y(nr),x(next),y(next)) ) then
    go to 14
  end if

  if (lft >= 0) then
    go to 12
  end if

  if ( .not. left(x(nl),y(nl),x0,y0,x(next),y(next)) ) then
    go to 14
  end if
!
!  Replace NL->NR with N0->NEXT.
!
  call swap (next,n0,nl,nr, list,lptr,lend, lp21)
  iwk(1,iwc) = n0
  iwk(2,iwc) = next
  go to 15
!
!  Swap NL-NR for N0-NEXT, shift columns IWC+1,...,IWL to
!  the left, and store N0-NEXT in the right portion of IWK.
!
12 continue

  call swap (next,n0,nl,nr, list,lptr,lend, lp21)

  do i = iwcp1,iwl
    iwk(1,i-1) = iwk(1,i)
    iwk(2,i-1) = iwk(2,i)
  end do

  iwk(1,iwl) = n0
  iwk(2,iwl) = next
  iwl = iwl - 1
  nr = next
  go to 11
!
!  A swap is not possible.  Set N0 to NR.
!
14 continue

  n0 = nr
  x0 = x(n0)
  y0 = y(n0)
  lft = 1
!
!  Advance to the next arc.
!
15 continue

  nr = next
  iwc = iwc + 1
  go to 11
!
!  NEXT LEFT N1->N2, NEXT /= N2, and IWC < IWL.
!  Test for a possible swap.
!
16 continue

  if ( .not. left(x(nl),y(nl),x0,y0,x(next),y(next)) ) then
    go to 19
  end if

  if (lft <= 0) then
    go to 17
  end if

  if ( .not. left(x0,y0,x(nr),y(nr),x(next),y(next)) ) then
    go to 19
  end if
!
!  Replace NL->NR with NEXT->N0.
!
  call swap (next,n0,nl,nr, list,lptr,lend, lp21)
  iwk(1,iwc) = next
  iwk(2,iwc) = n0
  go to 20
!
!  Swap NL-NR for N0-NEXT, shift columns IWF,...,IWC-1 to
!  the right, and store N0-NEXT in the left portion of IWK.
!
17 continue

  call swap (next,n0,nl,nr, list,lptr,lend, lp21)

  do i = iwc-1,iwf,-1
    iwk(1,i+1) = iwk(1,i)
    iwk(2,i+1) = iwk(2,i)
  end do

  iwk(1,iwf) = n0
  iwk(2,iwf) = next
  iwf = iwf + 1
  go to 20
!
!  A swap is not possible.  Set N0 to NL.
!
19 continue

  n0 = nl
  x0 = x(n0)
  y0 = y(n0)
  lft = -1
!
!  Advance to the next arc.
!
20 continue

  nl = next
  iwc = iwc + 1
  go to 11
!
!  N2 is opposite NL->NR (IWC = IWL).
!
21 continue

  if (n0 == n1) go to 24
  if (lft < 0) go to 22
!
!  N0 RIGHT N1->N2.  Test for a possible swap.
!
  if ( .not. left(x0,y0,x(nr),y(nr),x2,y2) ) go to 10
!
!  Swap NL-NR for N0-N2 and store N0-N2 in the right
!  portion of IWK.
!
  call swap (n2,n0,nl,nr, list,lptr,lend, lp21)
  iwk(1,iwl) = n0
  iwk(2,iwl) = n2
  iwl = iwl - 1
  go to 10
!
!  N0 LEFT N1->N2.  Test for a possible swap.
!
22 continue

  if ( .not. left(x(nl),y(nl),x0,y0,x2,y2) ) then
    go to 10
  end if
!
!  Swap NL-NR for N0-N2, shift columns IWF,...,IWL-1 to the
!  right, and store N0-N2 in the left portion of IWK.
!
  call swap ( n2, n0, nl, nr, list, lptr, lend, lp21 )
  i = iwl

23 continue

  iwk(1,i) = iwk(1,i-1)
  iwk(2,i) = iwk(2,i-1)
  i = i - 1

  if (i > iwf) then
    go to 23
  end if

  iwk(1,iwf) = n0
  iwk(2,iwf) = n2
  iwf = iwf + 1
  go to 10
!
!  IWF = IWC = IWL.  Swap out the last arc for N1-N2 and
!  store zeros in IWK.
!
24 continue

  call swap (n2,n1,nl,nr, list,lptr,lend, lp21)
  iwk(1,iwc) = 0
  iwk(2,iwc) = 0
!
!  Optimization procedure.
!
  if ( 1 < iwc ) then
!
!  Optimize the set of new arcs to the left of IN1->IN2.
!
    nit = 3*(iwc-1)
    call optim ( x, y, iwc-1, list, lptr, lend, nit, iwk, ierr )
    if ( ierr /= 0 ) then
      go to 34
    end if
  end if

  if ( iwc < iwend ) then
!
!  Optimize the set of new arcs to the right of IN1->IN2.
!
    nit = 3*(iwend-iwc)
    call optim ( x, y, iwend-iwc, list, lptr, lend, nit, &
      iwk(1,iwc+1), ierr )

    if (ierr /= 0) then
      go to 34
    end if

  end if
!
!  Successful termination.
!
  ier = 0
  return
!
!  Invalid triangulation data structure or collinear nodes
!  on convex hull boundary.
!
   33 ier = 3
  write (*,130) in1, in2
  130 format (//5x,'*** error in edge:  invalid triangula', &
          'tion or null triangles on boundary'/ &
          9x,'in1 =',i4,', in2=',i4/)
  return
!
! Error flag returned by OPTIM.
!
   34 ier = 4
  write (*,140) nit, ierr
  140 format (//5x,'*** error in optim:  nit = ',i4, &
          ', ier = ',i1,' ***'/)
  return
end
subroutine getnp ( ncc, lcc, n, x, y, list, lptr, lend, l, npts, ds, ier )

!*****************************************************************************80
!
!! GETNP sets the next nearest node to a given node.
!
!  Discussion:
!
!    Given a triangulation of N nodes and an array NPTS con-
!    taining the indexes of L-1 nodes ordered by distance from
!    NPTS(1), this subroutine sets NPTS(L) to the index of the
!    next node in the sequence -- the node, other than NPTS(1),
!    ...,NPTS(L-1), which is closest to NPTS(1).  Thus, the
!    ordered sequence of K closest nodes to N1 (including N1)
!    may be determined by K-1 calls to GETNP with NPTS(1) = N1
!    and L = 2,3,...,K for K >= 2.  Note that NPTS must 
!    include constraint nodes as well as non-constraint nodes.
!    Thus, a sequence of K1 closest non-constraint nodes to N1
!    must be obtained as a subset of the closest K2 nodes to N1
!    for some K2 >= K1.
!
!    The terms closest and distance have special definitions
!    when constraint nodes are present in the triangulation.
!    Nodes N1 and N2 are said to be visible from each other if
!    and only if the line segment N1-N2 intersects no constraint
!    arc (except possibly itself) and is not an interi-
!    or constraint arc (arc whose interior lies in a constraint
!    region).  A path from N1 to N2 is an ordered sequence of
!    nodes, with N1 first and N2 last, such that adjacent path
!    elements are visible from each other.  The path length is
!    the sum of the Euclidean distances between adjacent path
!    nodes.  Finally, the distance from N1 to N2 is defined to
!    be the length of the shortest path from N1 to N2.
!
!    The algorithm uses the property of a Delaunay triangulation
!    that the K-th closest node to N1 is a neighbor of one
!    of the K-1 closest nodes to N1.  With the definition of
!    distance used here, this property holds when constraints
!    are present as long as non-constraint arcs are locally
!    optimal.
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    Robert Renka,
!    Department of Computer Science,
!    University of North Texas,
!    renka@cs.unt.edu
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 751: TRIPACK, 
!    A Constrained Two-Dimensional Delaunay Triangulation Package,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 1, 1996.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NCC, the number of constraints.  NCC >= 0.
!
!    Input, integer ( kind = 4 ) LCC(*), a list of constraint curve starting indexes (or
!    dummy array of length 1 if NCC = 0).  Refer to subroutine ADDCST.
!
!    Input, integer ( kind = 4 ) N, the number of nodes in the triangulation.  N >= 3.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the coordinates of the nodes with 
!    non-constraint nodes in the first LCC(1)-1 locations if NCC > 0.
!
!    Input, integer ( kind = 4 ) LIST(*), LPTR(*), LEND(N), the triangulation data 
!    structure.  Refer to subroutine TRMESH.
!
!    Input, integer ( kind = 4 ) L, the number of nodes in the sequence on output.  
!    2 <= L <= N.
!
!    Input/output, integer ( kind = 4 ) NPTS(L), on input, the indexes of the L-1 closest
!    nodes to NPTS(1) in the first L-1 locations.  On output, updated with 
!    the index of the L-th closest node to NPTS(1) in position L unless
!    IER /= 0.
!
!    Input/output, real ( kind = 8 ) DS(L), the distance (defined above) 
!    between NPTS(1) and NPTS(I) in the I-th position for I = 1,...,L-1.  
!    Thus, DS(1) = 0.  On output, updated with the distance between NPTS(1) 
!    and NPTS(L) in position L unless IER /= 0.
!
!    Output, integer ( kind = 4 ) IER = Error indicator:
!     0 if no errors were encountered.
!    -1 if NCC, N, L, or an LCC entry is outside its valid range on input.
!     K if NPTS(K) is not a valid index in the range 1 to N.
!
  implicit none

  integer ( kind = 4 ) l
  integer ( kind = 4 ) n

  real ( kind = 8 ) dc
  real ( kind = 8 ) dl
  real ( kind = 8 ) ds(l)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) ifrst
  integer ( kind = 4 ) ilast
  logical intsec
  logical isw
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) km1
  integer ( kind = 4 ) lcc(*)
  integer ( kind = 4 ) lcc1
  integer ( kind = 4 ) lend(n)
  logical lft1
  logical lft2
  logical lft12
  integer ( kind = 4 ) list(*)
  integer ( kind = 4 ) lm1
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lpcl
  integer ( kind = 4 ) lpk
  integer ( kind = 4 ) lpkl
  integer ( kind = 4 ) lptr(*)
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) nc
  integer ( kind = 4 ) ncc
  logical ncf
  integer ( kind = 4 ) nf1
  integer ( kind = 4 ) nf2
  integer ( kind = 4 ) nj
  logical njf
  integer ( kind = 4 ) nk
  integer ( kind = 4 ) nkbak
  integer ( kind = 4 ) nkfor
  integer ( kind = 4 ) nl
  integer ( kind = 4 ) nn
  integer ( kind = 4 ) npts(l)
  logical skip
  logical sksav
  logical vis
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x1
  real ( kind = 8 ) xc
  real ( kind = 8 ) xj
  real ( kind = 8 ) xk
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) y1
  real ( kind = 8 ) yc
  real ( kind = 8 ) yj
  real ( kind = 8 ) yk
!
!  Store parameters in local variables and test for errors.
!  LCC1 indexes the first constraint node.
!
  ier = -1
  nn = n
  lcc1 = nn + 1
  lm1 = l - 1

  if ( ncc < 0 ) then
    return
  end if

  if ( lm1 < 1  .or.  lm1 >= nn) then
    return
  end if

  if (ncc == 0) then
    if (nn < 3) then
      return
    end if
  else
    do i = ncc,1,-1
      if (lcc1 - lcc(i) < 3) return
      lcc1 = lcc(i)
    end do

    if (lcc1 < 1) then
      return
    end if

  end if
!
!  Test for an invalid index in NPTS.
!
  do k = 1,lm1
    nk = npts(k)
    if ( nk < 1  .or.  nk > nn ) then
      ier = k
      return
    end if
  end do
!
!  Store N1 = NPTS(1) and mark the elements of NPTS.
!
  n1 = npts(1)
  x1 = x(n1)
  y1 = y(n1)

  do k = 1,lm1
    nk = npts(k)
    lend(nk) = -lend(nk)
  end do
!
!  Candidates NC for NL = NPTS(L) are the unmarked visible
!  neighbors of nodes NK in NPTS.  ISW is an initialization
!  switch set to .TRUE. when NL and its distance DL from N1
!  have been initialized with the first candidate encountered.
!
  isw = .false.
  dl = 0.0D+00
!
!  Loop on marked nodes NK = NPTS(K).  LPKL indexes the last
!  neighbor of NK in LIST.
!
  do k = 1,lm1

    km1 = k - 1
    nk = npts(k)
    xk = x(nk)
    yk = y(nk)
    lpkl = -lend(nk)
    nkfor = 0
    nkbak = 0
    vis = .true.
!
!  NK is a constraint node.  Set NKFOR and NKBAK to the
!  constraint nodes which follow and precede NK.  IFRST
!  and ILAST are set to the first and last nodes in the
!  constraint containing NK.
!
    if (nk >= lcc1) then

      ifrst = nn + 1

      do i = ncc,1,-1
        ilast = ifrst - 1
        ifrst = lcc(i)
        if ( nk >= ifrst ) then
          exit
        end if
      end do

      if (nk < ilast) then
        nkfor = nk + 1
      else
        nkfor = ifrst
      end if

      if (nk > ifrst) then
        nkbak = nk - 1
      else
        nkbak = ilast
      end if
!
!  Initialize VIS to TRUE iff NKFOR precedes NKBAK in the
!  adjacency list for NK -- the first neighbor is visible and is not NKBAK.
!
      lpk = lpkl

      do

        lpk = lptr(lpk)
        nc = abs ( list(lpk) )

        if ( nc == nkfor .or. nc == nkbak ) then
          exit
        end if

      end do

      vis = nc == nkfor

    end if
!
!  Loop on neighbors NC of NK, bypassing marked and nonvisible neighbors.
!
    lpk = lpkl

7   continue

    lpk = lptr(lpk)
    nc = abs ( list(lpk) )

    if ( nc == nkbak ) then
      vis = .true.
    end if
!
!  VIS = .FALSE. iff NK-NC is an interior constraint arc
!  (NK is a constraint node and NC lies strictly between
!  NKFOR and NKBAK).
!
    if ( .not. vis ) go to 15

    if ( nc == nkfor ) then
      vis = .false.
    end if

    if ( lend(nc) < 0 ) go to 15
!
!  Initialize distance DC between N1 and NC to Euclidean distance.
!
    xc = x(nc)
    yc = y(nc)
    dc = sqrt((xc-x1)*(xc-x1) + (yc-y1)*(yc-y1))
    if (isw  .and.  dc >= dl) go to 15

    if (k == 1) then
      go to 14
    end if
!
!  K >= 2.  Store the pointer LPCL to the last neighbor of NC.
!
    lpcl = lend(nc)
!
!  Set DC to the length of the shortest path from N1 to NC
!  which has not previously been encountered and which is
!  a viable candidate for the shortest path from N1 to NL.
!  This is Euclidean distance iff NC is visible from N1.
!  Since the shortest path from N1 to NL contains only ele-
!  ments of NPTS which are constraint nodes (in addition to
!  N1 and NL), only these need be considered for the path
!  from N1 to NC.  Thus, for distance function D(A,B) and
!  J = 1,...,K, DC = min(D(N1,NJ) + D(NJ,NC)) over con-
!  straint nodes NJ = NPTS(J) which are visible from NC.
!
    do j = 1,km1

      nj = npts(j)

      if ( 1 < j .and.  nj < lcc1 ) then
        go to 13
      end if
!
!  If NC is a visible neighbor of NJ, a path from N1 to NC
!  containing NJ has already been considered.  Thus, NJ may
!  be bypassed if it is adjacent to NC.
!
        lp = lpcl

8       continue

        lp = lptr(lp)

          if ( nj == abs ( list(lp) ) ) then
            go to 12
          end if

          if (lp /= lpcl) then
            go to 8
          end if
!
!  NJ is a constraint node (unless J=1) not adjacent to NC,
!  and is visible from NC iff NJ-NC is not intersected by
!  a constraint arc.  Loop on constraints I in reverse
!  order.
!
        xj = x(nj)
        yj = y(nj)
        ifrst = nn+1

        do 11 i = ncc,1,-1
          ilast = ifrst - 1
          ifrst = lcc(i)
          nf1 = ilast
          ncf = nf1 == nc
          njf = nf1 == nj
          skip = ncf  .or.  njf
!
!  Loop on boundary constraint arcs NF1-NF2 which contain
!  neither NC nor NJ.  NCF and NJF are TRUE iff NC (or NJ)
!  has been encountered in the constraint, and SKIP =
!  .TRUE. iff NF1 = NC or NF1 = NJ.
!
          do nf2 = ifrst,ilast

            if (nf2 == nc) ncf = .true.
            if (nf2 == nj) njf = .true.
            sksav = skip
            skip = nf2 == nc  .or.  nf2 == nj
!
!  The last constraint arc in the constraint need not be
!  tested if none of the arcs have been skipped.
!
            if ( sksav  .or.  skip  .or. &
                 (nf2 == ilast  .and. &
                 .not. ncf  .and.  .not. njf) ) then
              go to 9
            end if

            if ( intsec(x(nf1),y(nf1),x(nf2),y(nf2), &
                        xc,yc,xj,yj) ) then
              go to 12
            end if

9           continue

            nf1 = nf2

          end do
!
!  NC and NJ are constraint nodes in the same constraint.
!  NC-NJ is intersected by an interior constraint arc iff
!  1)  NC LEFT NF2->NF1 and (NJ LEFT NF1->NC and NJ LEFT NC->NF2) or
!  2)  NC .NOT. LEFT NF2->NF1 and (NJ LEFT NF1->NC or NJ LEFT NC->NF2),
!  where NF1, NC, NF2 are consecutive constraint nodes.
!
          if (.not. ncf  .or.  .not. njf) go to 11

          if (nc /= ifrst) then
            nf1 = nc - 1
          else
            nf1 = ilast
          end if

          if (nc /= ilast) then
            nf2 = nc + 1
          else
            nf2 = ifrst
          end if

          lft1 = (xc-x(nf1))*(yj-y(nf1)) >= (xj-x(nf1))*(yc-y(nf1))
          lft2 = (x(nf2)-xc)*(yj-yc) >= (xj-xc)*(y(nf2)-yc)
          lft12 = (x(nf1)-x(nf2))*(yc-y(nf2)) >= (xc-x(nf2))*(y(nf1)-y(nf2))

          if ( (lft1  .and.  lft2)  .or.  (.not. lft12 &
               .and.  (lft1  .or.  lft2)) ) go to 12

11         continue
!
!  NJ is visible from NC.  Exit the loop with DC = Euclidean
!  distance if J = 1.
!
        if (j == 1) then
          go to 14
        end if

        dc = min(dc,ds(j) + sqrt((xc-xj)*(xc-xj) + &
                    (yc-yj)*(yc-yj)))
        go to 13
!
!  NJ is not visible from NC or is adjacent to NC.  Initialize DC 
!  with D(N1,NK) + D(NK,NC) if J = 1.
!
12      continue

        if (j == 1) then
          dc = ds(k) + sqrt((xc-xk)*(xc-xk) &
                           + (yc-yk)*(yc-yk))
        end if

13      continue

    end do
!
!  Compare DC with DL.
!
      if ( isw  .and.  dc >= dl) then
        go to 15
      end if
!
!  The first (or a closer) candidate for NL has been encountered.
!
14    continue

      nl = nc
      dl = dc
      isw = .true.

15  continue

    if (lpk /= lpkl) then
      go to 7
    end if

  end do
!
!  Unmark the elements of NPTS and store NL and DL.
!
  do k = 1,lm1
    nk = npts(k)
    lend(nk) = -lend(nk)
  end do

  npts(l) = nl
  ds(l) = dl
  ier = 0

  return
end
function indxcc ( ncc, lcc, n, list, lend )

!*****************************************************************************80
!
!! INDXCC returns the index of an exterior constraint curve.
!
!  Discussion:
!
!    Given a constrained Delaunay triangulation, this 
!    function returns the index, if any, of an exterior constraint
!    curve (an unbounded constraint region).  An exterior 
!    constraint curve is assumed to be present if and only if the
!    clockwise-ordered sequence of boundary nodes is a 
!    subsequence of a constraint node sequence.  The triangulation
!    adjacencies corresponding to constraint edges may or may
!    not have been forced by a call to ADDCST, and the 
!    constraint region may or may not be valid (contain no nodes).
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    Robert Renka,
!    Department of Computer Science,
!    University of North Texas,
!    renka@cs.unt.edu
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 751: TRIPACK, 
!    A Constrained Two-Dimensional Delaunay Triangulation Package,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 1, 1996.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NCC, the number of constraints.  NCC >= 0.
!
!    Input, integer ( kind = 4 ) LCC(*), list of constraint curve starting indexes (or
!    dummy array of length 1 if NCC = 0).  Refer to subroutine ADDCST.
!
!    Input, integer ( kind = 4 ) N, the number of nodes in the triangulation.  N >= 3.
!
!    Input, integer ( kind = 4 ) LIST(*), LEND(N), the data structure defining the 
!    triangulation.  Refer to subroutine TRMESH.
!
!    Output, integer ( kind = 4 ) INDXCC, index of the exterior constraint curve, if
!    present, or 0 otherwise.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ifrst
  integer ( kind = 4 ) ilast
  integer ( kind = 4 ) indxcc
  integer ( kind = 4 ) lcc(*)
  integer ( kind = 4 ) lend(n)
  integer ( kind = 4 ) list(*)
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) n0
  integer ( kind = 4 ) ncc
  integer ( kind = 4 ) nst
  integer ( kind = 4 ) nxt

  indxcc = 0

  if ( ncc < 1 ) then
    return
  end if
!
!  Set N0 to the boundary node with smallest index.
!
  n0 = 0

  do

    n0 = n0 + 1
    lp = lend(n0)

    if ( list(lp) <= 0 ) then
      exit
    end if

  end do
!
!  Search in reverse order for the constraint I, if any, that
!  contains N0.  IFRST and ILAST index the first and last
!  nodes in constraint I.
!
  i = ncc
  ilast = n

  do

    ifrst = lcc(i)

    if ( ifrst <= n0 ) then
      exit
    end if

    if ( i == 1 ) then
      return
    end if

    i = i - 1
    ilast = ifrst - 1

  end do
!
!  N0 is in constraint I which indexes an exterior constraint
!  curve iff the clockwise-ordered sequence of boundary
!  node indexes beginning with N0 is increasing and bounded
!  above by ILAST.
!
  nst = n0

  do

    nxt = -list(lp)

    if ( nxt == nst ) then
      exit
    end if

    if ( nxt <= n0  .or. ilast < nxt ) then
      return
    end if

    n0 = nxt
    lp = lend(n0)

  end do
!
!  Constraint I contains the boundary node sequence as a subset.
!
  indxcc = i

  return
end
subroutine insert ( k, lp, list, lptr, lnew )

!*****************************************************************************80
!
!! INSERT inserts K as a neighbor of N1.
!
!  Discussion:
!
!    This subroutine inserts K as a neighbor of N1 following
!    N2, where LP is the LIST pointer of N2 as a neighbor of
!    N1.  Note that, if N2 is the last neighbor of N1, K will
!    become the first neighbor (even if N1 is a boundary node).
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    Robert Renka,
!    Department of Computer Science,
!    University of North Texas,
!    renka@cs.unt.edu
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 751: TRIPACK, 
!    A Constrained Two-Dimensional Delaunay Triangulation Package,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 1, 1996.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) K, the index of the node to be inserted.
!
!    Input, integer ( kind = 4 ) LP, the LIST pointer of N2 as a neighbor of N1.
!
!    Input/output, integer ( kind = 4 ) LIST(*), LPTR(*), LNEW, the data structure 
!    defining the triangulation.  Refer to subroutine TRMESH.  On output,
!    the data structure has been updated to include node K.
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) list(*)
  integer ( kind = 4 ) lnew
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lptr(*)
  integer ( kind = 4 ) lsav

  lsav = lptr(lp)
  lptr(lp) = lnew
  list(lnew) = k
  lptr(lnew) = lsav
  lnew = lnew + 1

  return
end
subroutine intadd ( kk, i1, i2, i3, list, lptr, lend, lnew )

!*****************************************************************************80
!
!! INTADD adds an interior point to a triangulation.
!
!  Discussion:
!
!    This subroutine adds an interior node to a triangulation
!    of a set of points in the plane.  The data structure is
!    updated with the insertion of node KK into the triangle
!    whose vertices are I1, I2, and I3.  No optimization of the
!    triangulation is performed.
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    Robert Renka,
!    Department of Computer Science,
!    University of North Texas,
!    renka@cs.unt.edu
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 751: TRIPACK, 
!    A Constrained Two-Dimensional Delaunay Triangulation Package,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 1, 1996.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) KK, the index of the node to be inserted.  1 <= KK
!    and KK must not be equal to I1, I2, or I3.
!
!    Input, integer ( kind = 4 ) I1, I2, I3, indexes of the counterclockwise-ordered
!    sequence of vertices of a triangle which contains node KK.
!
!    Input/output, integer ( kind = 4 ) LIST(*), LPTR(*), LEND(N), LNEW, the data 
!    structure defining the triangulation.  Refer to subroutine TRMESH. 
!    Triangle (I1,I2,I3) must be included in the triangulation.
!    On output, updated with the addition of node KK.  KK
!    will be connected to nodes I1, I2, and I3.
!
  implicit none

  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kk
  integer ( kind = 4 ) lend(*)
  integer ( kind = 4 ) list(*)
  integer ( kind = 4 ) lnew
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lptr(*)
  integer ( kind = 4 ) lstptr
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n3

  k = kk
!
!  Initialization.
!
  n1 = i1
  n2 = i2
  n3 = i3
!
!  Add K as a neighbor of I1, I2, and I3.
!
  lp = lstptr(lend(n1),n2,list,lptr)
  call insert (k,lp,list,lptr,lnew)
  lp = lstptr(lend(n2),n3,list,lptr)
  call insert (k,lp,list,lptr,lnew)
  lp = lstptr(lend(n3),n1,list,lptr)
  call insert (k,lp,list,lptr,lnew)
!
!  Add I1, I2, and I3 as neighbors of K.
!
  list(lnew) = n1
  list(lnew+1) = n2
  list(lnew+2) = n3
  lptr(lnew) = lnew + 1
  lptr(lnew+1) = lnew + 2
  lptr(lnew+2) = lnew
  lend(k) = lnew + 2
  lnew = lnew + 3

  return
end
function intsec ( x1, y1, x2, y2, x3, y3, x4, y4 )

!*****************************************************************************80
!
!! INTSEC determines if two line segments intersect.
!
!  Discussion:
!
!    Given a pair of line segments P1-P2 and P3-P4, this
!    function returns the value .TRUE. if and only if P1-P2
!    shares one or more points with P3-P4.  The line segments
!    include their endpoints, and the four points need not be
!    distinct.  Thus, either line segment may consist of a
!    single point, and the segments may meet in a V (which is
!    treated as an intersection).  Note that an incorrect
!    decision may result from floating point error if the four
!    endpoints are nearly collinear.
!
!  Modified:
!
!    19 May 2005
!
!  Author:
!
!    Robert Renka,
!    Department of Computer Science,
!    University of North Texas,
!    renka@cs.unt.edu
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 751: TRIPACK, 
!    A Constrained Two-Dimensional Delaunay Triangulation Package,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 1, 1996.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X1, Y1 = Coordinates of P1.
!
!    Input, real ( kind = 8 ) X2, Y2 = Coordinates of P2.
!
!    Input, real ( kind = 8 ) X3, Y3 = Coordinates of P3.
!
!    Input, real ( kind = 8 ) X4, Y4 = Coordinates of P4.
!
!    Output, logical INTSEC, the logical value defined above.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) d
  real ( kind = 8 ) dx12
  real ( kind = 8 ) dx31
  real ( kind = 8 ) dx34
  real ( kind = 8 ) dy12
  real ( kind = 8 ) dy31
  real ( kind = 8 ) dy34
  logical intsec
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ) x4
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) y3
  real ( kind = 8 ) y4
!
!  Test for overlap between the smallest rectangles that
!  contain the line segments and have sides parallel to
!  the axes.
!
  if (( x1 < x3 .and.  x1 < x4 .and.  x2 < x3 .and.  x2 < x4)  .or. &
      ( x3 < x1 .and.  x1 > x4 .and.  x2 > x3 .and.  x2 > x4)  .or. &
      ( y1 < y3 .and.  y1 < y4 .and.  y2 < y3 .and.  y2 < y4)  .or. &
      ( y3 < y1 .and.  y4 < y1 .and.  y2 > y3 .and.  y2 > y4)) then
    intsec = .false.
    return
  end if
!
!  Compute A = P4-P3 X P1-P3, B = P2-P1 X P1-P3, and
!  D = P2-P1 X P4-P3 (Z components).
!
  dx12 = x2 - x1
  dy12 = y2 - y1
  dx34 = x4 - x3
  dy34 = y4 - y3
  dx31 = x1 - x3
  dy31 = y1 - y3
  a = dx34 * dy31 - dx31 * dy34
  b = dx12 * dy31 - dx31 * dy12
  d = dx12 * dy34 - dx34 * dy12
!
!  D /= 0 and the point of intersection of the lines defined by the line
!  segments is P = P1 + (A/D)*(P2-P1) = P3 + (B/D)*(P4-P3).
!
  if ( d /= 0.0D+00 ) then

    intsec = &
      0.0D+00 <= a / d .and. a / d <= 1.0D+00 .and. &
      0.0D+00 <= b / d .and. b / d <= 1.0D+00
!
!  D == 0 and thus either the line segments are parallel,
!  or one (or both) of them is a single point.
!
  else

    intsec = ( a == 0.0D+00 .and. b == 0.0D+00 )

  end if

  return
end
function jrand ( n, ix, iy, iz )

!*****************************************************************************80
!
!! JRAND returns a uniformly distributed random integer between 1 and N.
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    Robert Renka,
!    Department of Computer Science,
!    University of North Texas,
!    renka@cs.unt.edu
!
!  Reference:
!
!    Brian Wichmann, David Hill, 
!    An Efficient and Portable Pseudo-random Number Generator,
!    Applied Statistics, 
!    Volume 31, Number 2, 1982, pages 188-190.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the maximum value to be returned.
!
!    Input/output, integer ( kind = 4 ) IX, IY, IZ, seeds initialized to values in
!    the range 1 to 30,000 before the first call to JRAND, and not altered 
!    by the user between subsequent calls (unless a sequence of random 
!    numbers is to be repeated by reinitializing the seeds).
!
!    Output, integer ( kind = 4 ) JRAND, random integer in the range 1 to N.
!
!  Local parameters:
!
!    U = Pseudo-random number uniformly distributed in the interval (0,1).
!    X = Pseudo-random number in the range 0 to 3 whose fractional part is U.
!
  implicit none

  integer ( kind = 4 ) ix
  integer ( kind = 4 ) iy
  integer ( kind = 4 ) iz
  integer ( kind = 4 ) jrand
  integer ( kind = 4 ) n
  real ( kind = 8 ) u
  real ( kind = 8 ) x

  ix = mod ( 171 * ix, 30269 )
  iy = mod ( 172 * iy, 30307 )
  iz = mod ( 170 * iz, 30323 )

  x = ( real ( ix, kind = 8 ) / 30269.0D+00 ) &
    + ( real ( iy, kind = 8 ) / 30307.0D+00 ) &
    + ( real ( iz, kind = 8 ) / 30323.0D+00 )
 
  u = x - int ( x )
  jrand = real ( n, kind = 8 ) * u + 1.0D+00

  return
end
function left ( x1, y1, x2, y2, x0, y0 )

!*****************************************************************************80
!
!! LEFT determines whether a node is to the left of a line.
!
!  Discussion:
!
!    This function determines whether node N0 is to the left
!    or to the right of the line through N1-N2 as viewed by an
!    observer at N1 facing N2.
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    Robert Renka,
!    Department of Computer Science,
!    University of North Texas,
!    renka@cs.unt.edu
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 751: TRIPACK, 
!    A Constrained Two-Dimensional Delaunay Triangulation Package,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 1, 1996.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X1, Y1, coordinates of N1.
!
!    Input, real ( kind = 8 ) X2, Y2, coordinates of N2.
!
!    Input, real ( kind = 8 ) X0, Y0, coordinates of N0.
!
!    Output, logical LEFT, is .TRUE. if and only if (X0,Y0) is on or 
!    to the left of the directed line N1->N2.
!
!  Local parameters:
!
!    DX1,DY1 = X,Y components of the vector N1->N2
!    DX2,DY2 = X,Y components of the vector N1->N0
!
  implicit none

  real ( kind = 8 ) dx1
  real ( kind = 8 ) dx2
  real ( kind = 8 ) dy1
  real ( kind = 8 ) dy2
  logical left
  real ( kind = 8 ) x0
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) y0
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2

  dx1 = x2 - x1
  dy1 = y2 - y1
  dx2 = x0 - x1
  dy2 = y0 - y1
!
!  If the sign of the vector cross product of N1->N2 and
!  N1->N0 is positive, then sin(A) > 0, where A is the
!  angle between the vectors, and thus A is in the range
!  (0,180) degrees.
!
  left = dx1 * dy2 >= dx2 * dy1

  return
end
function lstptr ( lpl, nb, list, lptr )

!*****************************************************************************80
!
!! LSTPTR returns the index of NB in the adjacency list for N0.
!
!  Discussion:
!
!    This function returns the index (LIST pointer) of NB in
!    the adjacency list for N0, where LPL = LEND(N0).
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    Robert Renka,
!    Department of Computer Science,
!    University of North Texas,
!    renka@cs.unt.edu
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 751: TRIPACK, 
!    A Constrained Two-Dimensional Delaunay Triangulation Package,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 1, 1996.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LPL = LEND(N0).
!
!    Input, integer ( kind = 4 ) NB, the index of the node whose pointer is to be 
!    returned.  NB must be connected to N0.
!
!    Input, integer ( kind = 4 ) LIST(*), LPTR(*), the data structure defining the 
!    triangulation.  Refer to subroutine TRMESH.
!
!    Output, integer ( kind = 4 ) LSTPTR, pointer such that LIST(LSTPTR) = NB or
!    LIST(LSTPTR) = -NB, unless NB is not a neighbor of N0, in which 
!    case LSTPTR = LPL.
!
  implicit none

  integer ( kind = 4 ) list(*)
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lpl
  integer ( kind = 4 ) lptr(*)
  integer ( kind = 4 ) lstptr
  integer ( kind = 4 ) nb
  integer ( kind = 4 ) nd

  lp = lptr(lpl)

  do

    nd = list(lp)

    if ( nd == nb ) then
      exit
    end if

    lp = lptr(lp)

    if ( lp == lpl ) then
      exit
    end if

  end do

  lstptr = lp

  return
end
function nbcnt ( lpl, lptr )

!*****************************************************************************80
!
!! NBCNT returns the number of neighbors of a node.
!
!  Discussion:
!
!    This function returns the number of neighbors of a node
!    in a triangulation created by TRMESH or TRMSHR.
!
!  Modified:
!
!    25 November 2002
!
!  Author:
!
!    Robert Renka,
!    Department of Computer Science,
!    University of North Texas,
!    renka@cs.unt.edu
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 751: TRIPACK, 
!    A Constrained Two-Dimensional Delaunay Triangulation Package,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 1, 1996.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LPL, the LIST pointer to the last neighbor of N0.
!    LPL = LEND(N0).
!
!    Input, integer ( kind = 4 ) LPTR(*), pointers associated with LIST.
!
!    Output, integer ( kind = 4 ) NBCNT, the  number of neighbors of N0.
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lpl
  integer ( kind = 4 ) lptr(*)
  integer ( kind = 4 ) nbcnt

  lp = lpl
  k = 1

  do

    lp = lptr(lp)
    if ( lp == lpl ) then
      exit
    end if

    k = k + 1

  end do

  nbcnt = k

  return
end
function nearnd ( xp, yp, ist, n, x, y, list, lptr, lend, dsq )

!*****************************************************************************80
!
!! NEARND finds the nearest triangulation node to a point.
!
!  Discussion:
!
!    Given a point P in the plane and a Delaunay triangulation created by
!    subroutine TRMESH or TRMSHR, this function returns the index of the 
!    nearest triangulation node to P.
!
!    The algorithm consists of implicitly adding P to the
!    triangulation, finding the nearest neighbor to P, and
!    implicitly deleting P from the triangulation.  Thus, it
!    is based on the fact that, if P is a node in a Delaunay
!    triangulation, the nearest node to P is a neighbor of P.
!
!    Note that the number of candidates for NEARND
!    (neighbors of P) is limited to LMAX defined in
!    the PARAMETER statement below.
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    Robert Renka,
!    Department of Computer Science,
!    University of North Texas,
!    renka@cs.unt.edu
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 751: TRIPACK, 
!    A Constrained Two-Dimensional Delaunay Triangulation Package,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 1, 1996.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XP, YP, the coordinates of the point P.
!
!    Input, integer ( kind = 4 ) IST, the index of a node at which TRFIND begins 
!    the search.  Search time depends on the proximity
!    of this node to P.
!
!    Input, integer ( kind = 4 ) N, the number of nodes in the triangulation.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the coordinates of the nodes.
!
!    Input, integer ( kind = 4 ) LIST(*), LPTR(*), LEND(N), a data structure 
!    defining the triangulation.  Refer to TRMESH.
!
!    Output, real ( kind = 8 ) DSQ, the square of the distance between P and
!    node NEARND.
!
!    Output, integer ( kind = 4 ) NEARND, the index of the nearest node to P, 
!    or 0 if N < 3 or the triangulation data structure is invalid.
!
  implicit none

  integer ( kind = 4 ), parameter :: lmax = 25
  integer ( kind = 4 ) n

  real ( kind = 8 ) cos1
  real ( kind = 8 ) cos2
  real ( kind = 8 ) ds1
  real ( kind = 8 ) dsq
  real ( kind = 8 ) dsr
  real ( kind = 8 ) dx11
  real ( kind = 8 ) dx12
  real ( kind = 8 ) dx21
  real ( kind = 8 ) dx22
  real ( kind = 8 ) dy11
  real ( kind = 8 ) dy12
  real ( kind = 8 ) dy21
  real ( kind = 8 ) dy22
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) ist
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lend(n)
  integer ( kind = 4 ) list(*)
  integer ( kind = 4 ) listp(lmax)
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lp1
  integer ( kind = 4 ) lp2
  integer ( kind = 4 ) lpl
  integer ( kind = 4 ) lptr(*)
  integer ( kind = 4 ) lptrp(lmax)
  integer ( kind = 4 ) lstptr
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n3
  integer ( kind = 4 ) nearnd
  integer ( kind = 4 ) nr
  integer ( kind = 4 ) nst
  real ( kind = 8 ) sin1
  real ( kind = 8 ) sin2
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xp
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) yp
!
!  Store local parameters and test for N invalid.
!
  nearnd = 0
  dsq = -1.0D+00

  if ( n < 3 ) then
    return
  end if

  nst = ist

  if ( nst < 1 .or. n < nst ) then
    nst = 1
  end if
!
!  Find a triangle (I1,I2,I3) containing P, or the rightmost
!  (I1) and leftmost (I2) visible boundary nodes as viewed from P.
!
  call trfind ( nst, xp, yp, n, x, y, list, lptr, lend, i1, i2, i3 )
!
!  Test for collinear nodes.
!
  if ( i1 == 0 ) then
    return
  end if
!
!  Store the linked list of 'neighbors' of P in LISTP and
!  LPTRP.  I1 is the first neighbor, and 0 is stored as
!  the last neighbor if P is not contained in a triangle.
!  L is the length of LISTP and LPTRP, and is limited to LMAX.
!
  if ( i3 /= 0 ) then

    listp(1) = i1
    lptrp(1) = 2
    listp(2) = i2
    lptrp(2) = 3
    listp(3) = i3
    lptrp(3) = 1
    l = 3

  else

    n1 = i1
    l = 1
    lp1 = 2
    listp(l) = n1
    lptrp(l) = lp1
!
!  Loop on the ordered sequence of visible boundary nodes
!  N1 from I1 to I2.
!
    do

      lpl = lend(n1)
      n1 = -list(lpl)
      l = lp1
      lp1 = l+1
      listp(l) = n1
      lptrp(l) = lp1

      if ( n1 == i2 .or. lmax <= lp1 ) then
        exit
      end if

    end do

    l = lp1
    listp(l) = 0
    lptrp(l) = 1

  end if
!
!  Initialize variables for a loop on arcs N1-N2 opposite P
!  in which new 'neighbors' are 'swapped' in.  N1 follows
!  N2 as a neighbor of P, and LP1 and LP2 are the LISTP
!  indexes of N1 and N2.
!
  lp2 = 1
  n2 = i1
  lp1 = lptrp(1)
  n1 = listp(lp1)
!
!  Begin loop:  find the node N3 opposite N1->N2.
!
  do

    lp = lstptr ( lend(n1), n2, list, lptr )

    if ( list(lp) < 0 ) then
      go to 4
    end if

    lp = lptr(lp)
    n3 = abs ( list(lp) )
!
!  Swap test:  Exit the loop if L = LMAX.
!
    if ( lmax <= l ) then
      exit
    end if

    dx11 = x(n1) - x(n3)
    dx12 = x(n2) - x(n3)
    dx22 = x(n2) - xp
    dx21 = x(n1) - xp

    dy11 = y(n1) - y(n3)
    dy12 = y(n2) - y(n3)
    dy22 = y(n2) - yp
    dy21 = y(n1) - yp

    cos1 = dx11 * dx12 + dy11 * dy12
    cos2 = dx22 * dx21 + dy22 * dy21

    if ( 0.0D+00 <= cos1 .and. cos2 >= 0.0D+00 ) then
      go to 4
    end if

    if ( cos1 < 0.0D+00 .and. cos2 < 0.0D+00 ) then
      go to 3
    end if

    sin1 = dx11 * dy12 - dx12 * dy11
    sin2 = dx22 * dy21 - dx21 * dy22

    if ( sin1 * cos2 + cos1 * sin2 >= 0.0D+00 ) then
      go to 4
    end if
!
!  Swap:  Insert N3 following N2 in the adjacency list for P.
!  The two new arcs opposite P must be tested.
!
3   continue

    l = l+1
    lptrp(lp2) = l
    listp(l) = n3
    lptrp(l) = lp1
    lp1 = l
    n1 = n3
    cycle
!
!  No swap:  Advance to the next arc and test for termination
!  on N1 = I1 (LP1 = 1) or N1 followed by 0.
!
4   continue

    if ( lp1 == 1 ) then
      exit
    end if

    lp2 = lp1
    n2 = n1
    lp1 = lptrp(lp1)
    n1 = listp(lp1)

    if ( n1 == 0 ) then
      exit
    end if

  end do
!
!  Set NR and DSR to the index of the nearest node to P and
!  its squared distance from P, respectively.
!
  nr = i1
  dsr = ( x(nr) - xp )**2 + ( y(nr) - yp )**2

  do lp = 2, l

    n1 = listp(lp)

    if ( n1 == 0 ) then
      cycle
    end if

    ds1 = ( x(n1) - xp )**2 + ( y(n1) - yp )**2

    if ( ds1 < dsr ) then
      nr = n1
      dsr = ds1
    end if

  end do

  dsq = dsr
  nearnd = nr

  return
end
subroutine optim ( x, y, na, list, lptr, lend, nit, iwk, ier )

!*****************************************************************************80
!
!! OPTIM optimizes the quadrilateral portion of a triangulation.
!
!  Discussion:
!
!    Given a set of NA triangulation arcs, this subroutine
!    optimizes the portion of the triangulation consisting of
!    the quadrilaterals (pairs of adjacent triangles) which
!    have the arcs as diagonals by applying the circumcircle
!    test and appropriate swaps to the arcs.
!
!    An iteration consists of applying the swap test and
!    swaps to all NA arcs in the order in which they are
!    stored.  The iteration is repeated until no swap occurs
!    or NIT iterations have been performed.  The bound on the
!    number of iterations may be necessary to prevent an
!    infinite loop caused by cycling (reversing the effect of a
!    previous swap) due to floating point inaccuracy when four
!    or more nodes are nearly cocircular.
!
!  Modified:
!
!    20 October 2005
!
!  Author:
!
!    Robert Renka,
!    Department of Computer Science,
!    University of North Texas,
!    renka@cs.unt.edu
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 751: TRIPACK, 
!    A Constrained Two-Dimensional Delaunay Triangulation Package,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 1, 1996.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X(*), Y(*), the nodal coordinates.
!
!    Input, integer ( kind = 4 ) NA, the number of arcs in the set.  0 <= NA.
!
!    Input/output, integer ( kind = 4 ) LIST(*), LPTR(*), LEND(N), data structure 
!    defining the triangulation.  Refer to subroutine TRMESH.
!    On output, updated to reflect the swaps.
!
!    Input/output, integer ( kind = 4 ) NIT.  On input, the maximum number of iterations 
!    to be performed.  A reasonable value is 3*NA.  1 <= NIT.  On output, 
!    the number of iterations performed.
!
!    Input/output, integer ( kind = 4 ) IWK(2,NA), containing the nodal indexes of the 
!    arc endpoints (pairs of endpoints are stored in columns).  On output,
!    the information has been updated to reflect the swaps.
!
!    Output, integer ( kind = 4 ) IER = Error indicator:
!    0, if no errors were encountered.
!    1, if a swap occurred on the last of MAXIT iterations, where MAXIT is 
!      the value of NIT on input.  The new set of arcs in not necessarily
!      optimal in this case.
!    2, if NA < 0 or NIT < 1 on input.
!    3, if IWK(2,I) is not a neighbor of IWK(1,I) for some I in the range 1
!      to NA.  A swap may have occurred in this case.
!    4, if a zero pointer was returned by subroutine SWAP.
!
!  Local parameters:
!
!    I =       Column index for IWK
!    IO1,IO2 = Nodal indexes of the endpoints of an arc in IWK
!    ITER =    Iteration count
!    LP =      LIST pointer
!    LP21 =    Parameter returned by SWAP (not used)
!    LPL =     Pointer to the last neighbor of IO1
!    LPP =     Pointer to the node preceding IO2 as a neighbor of IO1
!    MAXIT =   Input value of NIT
!    N1,N2 =   Nodes opposite IO1->IO2 and IO2->IO1, respectively
!    NNA =     Local copy of NA
!    SWP =     Flag set to TRUE iff a swap occurs in the optimization loop
!
  implicit none

  integer ( kind = 4 ) na

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) io1
  integer ( kind = 4 ) io2
  integer ( kind = 4 ) iter
  integer ( kind = 4 ) iwk(2,na)
  integer ( kind = 4 ) lend(*)
  integer ( kind = 4 ) list(*)
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lp21
  integer ( kind = 4 ) lpl
  integer ( kind = 4 ) lpp
  integer ( kind = 4 ) lptr(*)
  integer ( kind = 4 ) maxit
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) nit
  integer ( kind = 4 ) nna
  logical swp
  logical swptst
  real ( kind = 8 ) x(*)
  real ( kind = 8 ) y(*)

  nna = na
  maxit = nit

  if ( nna < 0 .or. maxit < 1 ) then
    nit = 0
    ier = 2
    return
  end if
!
!  Initialize iteration count ITER and test for NA = 0.
!
  iter = 0

  if ( nna == 0 ) then
    nit = 0
    ier = 0
    return
  end if
!
!  Top of loop.
!  SWP = TRUE iff a swap occurred in the current iteration.
!
  do

    if ( iter == maxit ) then
      nit = maxit
      ier = 1
      return
    end if

    iter = iter + 1
    swp = .false.
!
!  Inner loop on arcs IO1-IO2.
!
    do i = 1, nna

      io1 = iwk(1,i)
      io2 = iwk(2,i)
!
!  Set N1 and N2 to the nodes opposite IO1->IO2 and
!  IO2->IO1, respectively.  Determine the following:
!
!  LPL = pointer to the last neighbor of IO1,
!  LP = pointer to IO2 as a neighbor of IO1, and
!  LPP = pointer to the node N2 preceding IO2.
!
      lpl = lend(io1)
      lpp = lpl
      lp = lptr(lpp)

      do

        if ( list(lp) == io2 ) then
          go to 3
        end if

        lpp = lp
        lp = lptr(lpp)
        if ( lp == lpl ) then
          exit
        end if

      end do
!
!  IO2 should be the last neighbor of IO1.  Test for no
!  arc and bypass the swap test if IO1 is a boundary node.
!
      if ( abs ( list(lp) ) /= io2 ) then
        nit = iter
        ier = 3
        return
      end if

      if ( list(lp) < 0 ) then
        go to 4
      end if
!
!  Store N1 and N2, or bypass the swap test if IO1 is a
!  boundary node and IO2 is its first neighbor.
!
3     continue

      n2 = list(lpp)

      if ( n2 < 0 ) then
        go to 4
      end if

      lp = lptr(lp)
      n1 = abs ( list(lp) )
!
!  Test IO1-IO2 for a swap, and update IWK if necessary.
!
      if ( .not. swptst ( n1, n2, io1, io2, x, y ) ) then
        go to 4
      end if

      call swap ( n1, n2, io1, io2, list, lptr, lend, lp21 )

      if ( lp21 == 0 ) then
        nit = iter
        ier = 4
        return
      end if

      swp = .true.
      iwk(1,i) = n1
      iwk(2,i) = n2

4     continue

    end do

    if ( .not. swp ) then
      exit
    end if

  end do
!
!  Successful termination.
!
5 continue

  nit = iter
  ier = 0

  return
end
function store ( x )

!*****************************************************************************80
!
!! STORE forces its argument to be stored.
!
!  Discussion:
!
!    This function forces its argument X to be stored in a
!    memory location, thus providing a means of determining
!    floating point number characteristics (such as the machine
!    precision) when it is necessary to avoid computation in
!    high precision registers.
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    Robert Renka,
!    Department of Computer Science,
!    University of North Texas,
!    renka@cs.unt.edu
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 751: TRIPACK, 
!    A Constrained Two-Dimensional Delaunay Triangulation Package,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 1, 1996.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the value to be stored.
!
!    Output, real ( kind = 8 ) STORE, the value of X after it has been stored
!    and possibly truncated or rounded to the single precision word length.
!
  implicit none

  real ( kind = 8 ) store
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  common /stcom/ y

  y = x
  store = y

  return
end
subroutine swap ( in1, in2, io1, io2, list, lptr, lend, lp21 )

!*****************************************************************************80
!
!! SWAP adjusts a triangulation by swapping a diagonal arc.
!
!  Discussion:
!
!    Given a triangulation of a set of points on the unit
!    sphere, this subroutine replaces a diagonal arc in a
!    strictly convex quadrilateral (defined by a pair of adja-
!    cent triangles) with the other diagonal.  Equivalently, a
!    pair of adjacent triangles is replaced by another pair
!    having the same union.
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    Robert Renka,
!    Department of Computer Science,
!    University of North Texas,
!    renka@cs.unt.edu
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 751: TRIPACK, 
!    A Constrained Two-Dimensional Delaunay Triangulation Package,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 1, 1996.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IN1, IN2, IO1, IO2, the nodal indexes of the vertices of
!    the quadrilateral.  IO1-IO2 is replaced by IN1-IN2.  (IO1,IO2,IN1)
!    and (IO2,IO1,IN2) must be triangles on input.
!
!    Input/output, integer ( kind = 4 ) LIST(*), LPTR(*), LEND(N), the data structure
!    defining the triangulation.  Refer to subroutine TRMESH.  On output,
!    updated with the swap; triangles (IO1,IO2,IN1) and (IO2,IO1,IN2) are
!    replaced by (IN1,IN2,IO2) and (IN2,IN1,IO1) unless LP21 = 0.
!
!    Output, integer ( kind = 4 ) LP21, the index of IN1 as a neighbor of IN2 after the
!    swap is performed unless IN1 and IN2 are adjacent on input, in which 
!    case LP21 = 0.
!
!  Local parameters:
!
!    LP, LPH, LPSAV = LIST pointers
!
  implicit none

  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) io1
  integer ( kind = 4 ) io2
  integer ( kind = 4 ) lend(*)
  integer ( kind = 4 ) list(*)
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lp21
  integer ( kind = 4 ) lph
  integer ( kind = 4 ) lpsav
  integer ( kind = 4 ) lptr(*)
  integer ( kind = 4 ) lstptr
!
!  Test for IN1 and IN2 adjacent.
!
  lp = lstptr(lend(in1),in2,list,lptr)

  if ( abs ( list(lp) ) == in2 ) then
    lp21 = 0
    return
  end if
!
!  Delete IO2 as a neighbor of IO1.
!
  lp = lstptr(lend(io1),in2,list,lptr)
  lph = lptr(lp)
  lptr(lp) = lptr(lph)
!
!  If IO2 is the last neighbor of IO1, make IN2 the last neighbor.
!
  if ( lend(io1) == lph ) then
    lend(io1) = lp
  end if
!
!  Insert IN2 as a neighbor of IN1 following IO1
!  using the hole created above.
!
  lp = lstptr(lend(in1),io1,list,lptr)
  lpsav = lptr(lp)
  lptr(lp) = lph
  list(lph) = in2
  lptr(lph) = lpsav
!
!  Delete IO1 as a neighbor of IO2.
!
  lp = lstptr(lend(io2),in1,list,lptr)
  lph = lptr(lp)
  lptr(lp) = lptr(lph)
!
!  If IO1 is the last neighbor of IO2, make IN1 the last neighbor.
!
  if ( lend(io2) == lph ) then
    lend(io2) = lp
  end if
!
!  Insert IN1 as a neighbor of IN2 following IO2.
!
  lp = lstptr(lend(in2),io2,list,lptr)
  lpsav = lptr(lp)
  lptr(lp) = lph
  list(lph) = in1
  lptr(lph) = lpsav
  lp21 = lph

  return
end
function swptst ( in1, in2, io1, io2, x, y )

!*****************************************************************************80
!
!! SWPTST applies the circumcircle test to a quadrilateral.
!
!  Discussion:
!
!    This function applies the circumcircle test to a quadri-
!    lateral defined by a pair of adjacent triangles.  The
!    diagonal arc (shared triangle side) should be swapped for
!    the other diagonl if and only if the fourth vertex is
!    strictly interior to the circumcircle of one of the
!    triangles (the decision is independent of the choice of
!    triangle).  Equivalently, the diagonal is chosen to maxi-
!    mize the smallest of the six interior angles over the two
!    pairs of possible triangles (the decision is for no swap
!    if the quadrilateral is not strictly convex).
!
!    When the four vertices are nearly cocircular (the
!    neutral case), the preferred decision is no swap -- in
!    order to avoid unnecessary swaps and, more important, to
!    avoid cycling in subroutine OPTIM which is called by
!    DELNOD and EDGE.  Thus, a tolerance SWTOL (stored in
!    SWPCOM by TRMESH or TRMSHR) is used to define 'nearness'
!    to the neutral case.
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    Robert Renka,
!    Department of Computer Science,
!    University of North Texas,
!    renka@cs.unt.edu
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 751: TRIPACK, 
!    A Constrained Two-Dimensional Delaunay Triangulation Package,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 1, 1996.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IN1, IN2, IO1, IO2, the nodal indexes of the vertices of
!    the quadrilateral.  IO1-IO2 is the triangulation arc (shared triangle
!    side) to be replaced by IN1-IN2 if the decision is to swap.  The
!    triples (IO1,IO2,IN1) and (IO2,IO1,IN2) must define triangles (be
!    in counterclockwise order) on input.
!
!    Input, real ( kind = 8 ) X(*), Y(*), the nodal coordinates.
!
!    Output, logical SWPTST, .TRUE. if and only if the arc connecting
!    IO1 and IO2 is to be replaced.
!
!  Local parameters:
!
!    DX11,DY11 = X,Y components of the vector IN1->IO1
!    DX12,DY12 = X,Y components of the vector IN1->IO2
!    DX22,DY22 = X,Y components of the vector IN2->IO2
!    DX21,DY21 = X,Y components of the vector IN2->IO1
!    SIN1 =      Cross product of the vectors IN1->IO1 and
!                IN1->IO2 -- proportional to sin(T1), where
!                T1 is the angle at IN1 formed by the vectors
!    COS1 =      Inner product of the vectors IN1->IO1 and
!                IN1->IO2 -- proportional to cos(T1)
!    SIN2 =      Cross product of the vectors IN2->IO2 and
!                IN2->IO1 -- proportional to sin(T2), where
!                T2 is the angle at IN2 formed by the vectors
!    COS2 =      Inner product of the vectors IN2->IO2 and
!                IN2->IO1 -- proportional to cos(T2)
!    SIN12 =     SIN1*COS2 + COS1*SIN2 -- proportional to sin(T1+T2)
!
  implicit none

  real ( kind = 8 ) cos1
  real ( kind = 8 ) cos2
  real ( kind = 8 ) dx11
  real ( kind = 8 ) dx12
  real ( kind = 8 ) dx21
  real ( kind = 8 ) dx22
  real ( kind = 8 ) dy11
  real ( kind = 8 ) dy12
  real ( kind = 8 ) dy21
  real ( kind = 8 ) dy22
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) io1
  integer ( kind = 4 ) io2
  real ( kind = 8 ) sin1
  real ( kind = 8 ) sin12
  real ( kind = 8 ) sin2
  logical swptst
  real ( kind = 8 ) swtol
  real ( kind = 8 ) x(*)
  real ( kind = 8 ) y(*)
!
!  Tolerance stored by TRMESH or TRMSHR.
!
  common /swpcom/ swtol
!
!  Compute the vectors containing the angles T1 and T2.
!
  dx11 = x(io1) - x(in1)
  dx12 = x(io2) - x(in1)
  dx22 = x(io2) - x(in2)
  dx21 = x(io1) - x(in2)

  dy11 = y(io1) - y(in1)
  dy12 = y(io2) - y(in1)
  dy22 = y(io2) - y(in2)
  dy21 = y(io1) - y(in2)
!
!  Compute inner products.
!
  cos1 = dx11 * dx12 + dy11 * dy12
  cos2 = dx22 * dx21 + dy22 * dy21
!
!  The diagonals should be swapped iff 180 < (T1+T2)
!  degrees.  The following two tests ensure numerical
!  stability:  the decision must be FALSE when both
!  angles are close to 0, and TRUE when both angles
!  are close to 180 degrees.
!
  if ( 0.0D+00 <= cos1 .and. 0.0D+00 <= cos2 ) then
    swptst = .false.
    return
  end if

  if ( cos1 < 0.0D+00 .and. cos2 < 0.0D+00 ) then
    swptst = .true.
    return
  end if
!
!  Compute vector cross products (Z-components).
!
  sin1 = dx11 * dy12 - dx12 * dy11
  sin2 = dx22 * dy21 - dx21 * dy22
  sin12 = sin1 * cos2 + cos1 * sin2

  if ( -swtol <= sin12 ) then
    swptst = .false.
  else
    swptst = .true.
  end if

  return
end
subroutine trfind ( nst, px, py, n, x, y, list, lptr, lend, i1, i2, i3 )

!*****************************************************************************80
!
!! TRFIND locates a point relative to a triangulation.
!
!  Discussion:
!
!    This subroutine locates a point P relative to a triangu-
!    lation created by subroutine TRMESH or TRMSHR.  If P is
!    contained in a triangle, the three vertex indexes are
!    returned.  Otherwise, the indexes of the rightmost and
!    leftmost visible boundary nodes are returned.
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    Robert Renka,
!    Department of Computer Science,
!    University of North Texas,
!    renka@cs.unt.edu
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 751: TRIPACK, 
!    A Constrained Two-Dimensional Delaunay Triangulation Package,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 1, 1996.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NST, the index of a node at which TRFIND begins the
!    search.  Search time depends on the proximity of this node to P.
!
!    Input, real ( kind = 8 ) PX, PY, the coordinates of the point P to be
!    located.
!
!    Input, integer ( kind = 4 ) N, the number of nodes in the triangulation.  3 <= N.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the coordinates of the nodes in
!    the triangulation.
!
!    Input, integer ( kind = 4 ) LIST(*), LPTR(*), LEND(N), the data structure defining 
!    the triangulation.  Refer to subroutine TRMESH.
!
!    Output, integer ( kind = 4 ) I1, I2, I3, nodal indexes, in counterclockwise order,
!    of the vertices of a triangle containing P if P is contained in a 
!    triangle.  If P is not in the convex hull of the nodes, I1 indexes 
!    the rightmost visible boundary node, I2 indexes the leftmost visible
!    boundary node, and I3 = 0.  Rightmost and leftmost are defined from 
!    the perspective of P, and a pair of points are visible from each 
!    other if and only if the line segment joining them intersects no 
!    triangulation arc.  If P and all of the nodes lie on a common line, 
!    then I1 = I2 = I3 = 0 on output.
!
!  Local parameters:
!
!    B1,B2 =    Unnormalized barycentric coordinates of P with respect 
!               to (N1,N2,N3)
!    IX,IY,IZ = Integer seeds for JRAND
!    LP =       LIST pointer
!    N0,N1,N2 = Nodes in counterclockwise order defining a
!               cone (with vertex N0) containing P
!    N1S,N2S =  Saved values of N1 and N2
!    N3,N4 =    Nodes opposite N1->N2 and N2->N1, respectively
!    NB =       Index of a boundary node -- first neighbor of
!               NF or last neighbor of NL in the boundary traversal loops
!    NF,NL =    First and last neighbors of N0, or first
!               (rightmost) and last (leftmost) nodes
!               visible from P when P is exterior to the triangulation
!    NP,NPP =   Indexes of boundary nodes used in the boundary traversal loops
!    XA,XB,XC = Dummy arguments for FRWRD
!    YA,YB,YC = Dummy arguments for FRWRD
!    XP,YP =    Local variables containing the components of P
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) b1
  real ( kind = 8 ) b2
  logical frwrd
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ), save :: ix = 1
  integer ( kind = 4 ), save :: iy = 2
  integer ( kind = 4 ), save :: iz = 3
  integer ( kind = 4 ) jrand
  logical left
  integer ( kind = 4 ) lend(n)
  integer ( kind = 4 ) list(*)
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lptr(*)
  integer ( kind = 4 ) lstptr
  integer ( kind = 4 ) n0
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n1s
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n2s
  integer ( kind = 4 ) n3
  integer ( kind = 4 ) n4
  integer ( kind = 4 ) nb
  integer ( kind = 4 ) nf
  integer ( kind = 4 ) nl
  integer ( kind = 4 ) np
  integer ( kind = 4 ) npp
  integer ( kind = 4 ) nst
  real ( kind = 8 ) px
  real ( kind = 8 ) py
  real ( kind = 8 ) store
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xa
  real ( kind = 8 ) xb
  real ( kind = 8 ) xc
  real ( kind = 8 ) xp
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) ya
  real ( kind = 8 ) yb
  real ( kind = 8 ) yc
  real ( kind = 8 ) yp
!
!  Statement function:
!
!  FRWRD = TRUE iff C is forward of A->B iff <A->B,A->C> >= 0.
!
  frwrd(xa,ya,xb,yb,xc,yc) = (xb-xa)*(xc-xa) + (yb-ya)*(yc-ya) >= 0.0D+00
!
!  Initialize variables.
!
  xp = px
  yp = py
  n0 = nst

  if ( n0 < 1  .or.  n < n0 ) then
    n0 = jrand ( n, ix, iy, iz )
  end if
!
!  Set NF and NL to the first and last neighbors of N0, and
!  initialize N1 = NF.
!
1 continue

  lp = lend(n0)
  nl = list(lp)
  lp = lptr(lp)
  nf = list(lp)
  n1 = nf
!
!  Find a pair of adjacent neighbors N1,N2 of N0 that define
!  a wedge containing P:  P LEFT N0->N1 and P RIGHT N0->N2.
!
  if ( 0 < nl ) then
    go to 2
  end if
!
!   N0 is a boundary node.  Test for P exterior.
!
  nl = -nl

  if ( .not. left ( x(n0), y(n0), x(nf), y(nf), xp, yp ) ) then
    nl = n0
    go to 9
  end if

  if ( .not. left(x(nl),y(nl),x(n0),y(n0),xp,yp) ) then
    nb = nf
    nf = n0
    np = nl
    npp = n0
    go to 11
  end if

  go to 3
!
!  N0 is an interior node.  Find N1.
!
2 continue

    do

      if ( left(x(n0),y(n0),x(n1),y(n1),xp,yp) ) then
        exit
      end if

      lp = lptr(lp)
      n1 = list(lp)

      if ( n1 == nl ) then
        go to 6
      end if

    end do
!
!  P is to the left of edge N0->N1.  Initialize N2 to the
!  next neighbor of N0.
!
3 continue

    lp = lptr(lp)
    n2 = abs ( list(lp) )

    if ( .not. left(x(n0),y(n0),x(n2),y(n2),xp,yp) ) then
      go to 7
    end if

    n1 = n2
    if ( n1 /= nl ) then
      go to 3
    end if

  if ( .not. left(x(n0),y(n0),x(nf),y(nf),xp,yp) ) then
    go to 6
  end if

  if (xp == x(n0) .and. yp == y(n0)) then
    go to 5
  end if
!
!  P is left of or on edges N0->NB for all neighbors NB of N0.
!  All points are collinear iff P is left of NB->N0 for
!  all neighbors NB of N0.  Search the neighbors of N0.
!  NOTE: N1 = NL and LP points to NL.
!
4   continue

    if ( .not. left(x(n1),y(n1),x(n0),y(n0),xp,yp) ) then
      go to 5
    end if

    lp = lptr(lp)
    n1 = abs ( list(lp) )

    if ( n1 == nl ) then
      i1 = 0
      i2 = 0
      i3 = 0
      return
    end if

    go to 4
!
!  P is to the right of N1->N0, or P=N0.  Set N0 to N1 and start over.
!
5 continue

  n0 = n1
  go to 1
!
!  P is between edges N0->N1 and N0->NF.
!
6 continue

  n2 = nf
!
!  P is contained in the wedge defined by line segments
!  N0->N1 and N0->N2, where N1 is adjacent to N2.  Set
!  N3 to the node opposite N1->N2, and save N1 and N2 to
!  test for cycling.
!
7 continue

  n3 = n0
  n1s = n1
  n2s = n2
!
!  Top of edge hopping loop.  Test for termination.
!
8 continue

  if ( left ( x(n1), y(n1), x(n2), y(n2), xp, yp ) ) then
!
!  P LEFT N1->N2 and hence P is in (N1,N2,N3) unless an
!  error resulted from floating point inaccuracy and
!  collinearity.  Compute the unnormalized barycentric
!  coordinates of P with respect to (N1,N2,N3).
!
    b1 = (x(n3)-x(n2))*(yp-y(n2)) - (xp-x(n2))*(y(n3)-y(n2))
    b2 = (x(n1)-x(n3))*(yp-y(n3)) - (xp-x(n3))*(y(n1)-y(n3))

    if ( store ( b1 + 1.0D+00 ) >= 1.0D+00  .and. &
         store ( b2 + 1.0D+00 ) >= 1.0D+00 ) then
      go to 16
    end if
!
!  Restart with N0 randomly selected.
!
    n0 = jrand ( n, ix, iy, iz )
    go to 1

  end if
!
!  Set N4 to the neighbor of N2 which follows N1 (node
!  opposite N2->N1) unless N1->N2 is a boundary edge.
!
  lp = lstptr(lend(n2),n1,list,lptr)

  if ( list(lp) < 0 ) then
    nf = n2
    nl = n1
    go to 9
  end if

  lp = lptr(lp)
  n4 = abs ( list(lp) )
!
!  Select the new edge N1->N2 which intersects the line
!  segment N0-P, and set N3 to the node opposite N1->N2.
!
  if ( left(x(n0),y(n0),x(n4),y(n4),xp,yp) ) then
    n3 = n1
    n1 = n4
    n2s = n2
    if (n1 /= n1s  .and.  n1 /= n0) go to 8
  else
    n3 = n2
    n2 = n4
    n1s = n1
    if ( n2 /= n2s  .and.  n2 /= n0 ) then
      go to 8
    end if
  end if
!
!  The starting node N0 or edge N1-N2 was encountered
!  again, implying a cycle (infinite loop).  Restart
!  with N0 randomly selected.
!
  n0 = jrand ( n, ix, iy, iz )
  go to 1
!
!  Boundary traversal loops.  NL->NF is a boundary edge and
!  P RIGHT NL->NF.  Save NL and NF.

9 continue

  np = nl
  npp = nf
!
!  Find the first (rightmost) visible boundary node NF.  NB
!  is set to the first neighbor of NF, and NP is the last neighbor.
!
10 continue

  lp = lend(nf)
  lp = lptr(lp)
  nb = list(lp)

  if ( .not. left(x(nf),y(nf),x(nb),y(nb),xp,yp) ) then
    go to 12
  end if
!
!  P LEFT NF->NB and thus NB is not visible unless an error
!  resulted from floating point inaccuracy and collinear-
!  ity of the 4 points NP, NF, NB, and P.
!
11 continue

  if ( frwrd(x(nf),y(nf),x(np),y(np),xp,yp)  .or. &
       frwrd(x(nf),y(nf),x(np),y(np),x(nb),y(nb)) ) then
    i1 = nf
    go to 13
  end if
!
!  Bottom of loop.
!
12 continue

  np = nf
  nf = nb
  go to 10
!
!  Find the last (leftmost) visible boundary node NL.  NB
!  is set to the last neighbor of NL, and NPP is the first
!  neighbor.
!
13 continue

  lp = lend(nl)
  nb = -list(lp)

  if ( .not. left(x(nb),y(nb),x(nl),y(nl),xp,yp) ) then
    go to 14
  end if
!
!  P LEFT NB->NL and thus NB is not visible unless an error
!  resulted from floating point inaccuracy and collinear-
!  ity of the 4 points P, NB, NL, and NPP.
!
  if ( frwrd(x(nl),y(nl),x(npp),y(npp),xp,yp)  .or. &
       frwrd(x(nl),y(nl),x(npp),y(npp),x(nb),y(nb)) ) then
    go to 15
  end if
!
!  Bottom of loop.
!
14 continue

  npp = nl
  nl = nb
  go to 13
!
!  NL is the leftmost visible boundary node.
!
15 continue

  i2 = nl
  i3 = 0
  return
!
!  P is in the triangle (N1,N2,N3).
!
16 continue
 
  i1 = n1
  i2 = n2
  i3 = n3

  return
end
subroutine trlist ( ncc, lcc, n, list, lptr, lend, nrow, nt, ltri, lct, ier )

!*****************************************************************************80
!
!! TRLIST converts a triangulation to triangle list form.
!
!  Discussion:
!
!    This subroutine converts a triangulation data structure
!    from the linked list created by subroutine TRMESH or
!    TRMSHR to a triangle list.
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    Robert Renka,
!    Department of Computer Science,
!    University of North Texas,
!    renka@cs.unt.edu
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 751: TRIPACK, 
!    A Constrained Two-Dimensional Delaunay Triangulation Package,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 1, 1996.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NCC, the number of constraints.  NCC >= 0.
!
!    Input, integer ( kind = 4 ) LCC(*), list of constraint curve starting indexes (or
!    dummy array of length 1 if NCC = 0).  Refer to subroutine ADDCST.
!
!    Input, integer ( kind = 4 ) N, the number of nodes in the triangulation.  N >= 3.
!
!    Input, integer ( kind = 4 ) LIST(*), LPTR(*), LEND(N), linked list data structure 
!    defining the triangulation.  Refer to subroutine TRMESH.
!
!    Input, integer ( kind = 4 ) NROW, the number of rows (entries per triangle) 
!    reserved for the triangle list LTRI.  The value must be 6 if only 
!    the vertex indexes and neighboring triangle indexes are to be
!    stored, or 9 if arc indexes are also to be assigned and stored.  
!    Refer to LTRI.
!
!    Input, integer ( kind = 4 ) LTRI(NROW*NT), where NT is at most 2N-5.  (A sufficient
!    length is 12 * N if NROW=6 or 18*N if NROW=9.)
!
!    Output, integer ( kind = 4 ) NT, the number of triangles in the triangulation unless
!    IER /= 0, in which case NT = 0.  NT = 2N - NB- 2, where NB is the number 
!    of boundary nodes.
!
!    Output, integer ( kind = 4 ) LTRI(NROW,NT), whose J-th column contains the vertex nodal
!    indexes (first three rows), neighboring triangle indexes (second three
!    rows), and, if NROW = 9, arc indexes (last three rows) associated with
!    triangle J for J = 1,...,NT.  The vertices are ordered counterclockwise
!    with the first vertex taken to be the one with smallest index.  Thus,
!    LTRI(2,J) and LTRI(3,J) are larger than LTRI(1,J) and index adjacent
!    neighbors of node LTRI(1,J).  For I = 1,2,3, LTRI(I+3,J) and LTRI(I+6,J)
!    index the triangle and arc, respectively, which are opposite (not shared
!    by) node LTRI(I,J), with LTRI(I+3,J) = 0 if LTRI(I+6,J) indexes a boundary
!    arc.  Vertex indexes range from 1 to N, triangle indexes from 0 to NT,
!    and, if included, arc indexes from 1 to NA = NT+N-1.  The triangles are 
!    ordered on first (smallest) vertex indexes, except that the sets of
!    constraint triangle (triangles contained in the closure of a constraint
!    region) follow the non-constraint triangles.
!
!    Output, integer ( kind = 4 ) LCT(NCC), containing the triangle index of the first
!    triangle of constraint J in LCT(J).  Thus, the number of non-constraint
!    triangles is LCT(1)-1, and constraint J contains LCT(J+1)-LCT(J) 
!    triangles, where LCT(NCC+1) = NT+1.
!
!    Output, integer ( kind = 4 ) IER = Error indicator.
!    0, if no errors were encountered.
!    1, if NCC, N, NROW, or an LCC entry is outside its valid range on input.
!    2, if the triangulation data structure (LIST,LPTR,LEND) is invalid.  
!
!  Local Parameters:
!
!    ARCS = TRUE iff arc indexes are to be stored.
!    KA,KT = Numbers of currently stored arcs and triangles.
!    N1ST = Starting index for the loop on nodes (N1ST = 1 on
!           pass 1, and N1ST = LCC1 on pass 2).
!    NM2 = Upper bound on candidates for N1.
!    PASS2 = TRUE iff constraint triangles are to be stored.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nrow

  logical arcs
  logical cstri
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) isv
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jlast
  integer ( kind = 4 ) ka
  integer ( kind = 4 ) kn
  integer ( kind = 4 ) kt
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lcc(*)
  integer ( kind = 4 ) lcc1
  integer ( kind = 4 ) lct(*)
  integer ( kind = 4 ) lend(n)
  integer ( kind = 4 ) list(*)
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lp2
  integer ( kind = 4 ) lpl
  integer ( kind = 4 ) lpln1
  integer ( kind = 4 ) lptr(*)
  integer ( kind = 4 ) ltri(nrow,*)
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n1st
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n3
  integer ( kind = 4 ) ncc
  integer ( kind = 4 ) nm2
  integer ( kind = 4 ) nn
  integer ( kind = 4 ) nt
  logical pass2
!
!  Test for invalid input parameters and store the index
!  LCC1 of the first constraint node (if any).
!
  nn = n

  if ( ncc < 0 .or. ( nrow /= 6  .and. nrow /= 9 ) ) then
    nt = 0
    ier = 1
    return
  end if

  lcc1 = nn+1

  if (ncc == 0) then

    if ( nn < 3 ) then
      nt = 0
      ier = 1
      return
    end if

  else

    do i = ncc, 1, -1
      if ( lcc1 - lcc(i) < 3 ) then
        nt = 0
        ier = 1
        return
      end if
      lcc1 = lcc(i)
    end do

    if ( lcc1 < 1 ) then
      nt = 0
      ier = 1
      return
    end if

  end if
!
!  Initialize parameters for loop on triangles KT = (N1,N2,
!  N3), where N1 < N2 and N1 < N3.  This requires two
!  passes through the nodes with all non-constraint
!  triangles stored on the first pass, and the constraint
!  triangles stored on the second.
!
  arcs = nrow == 9
  ka = 0
  kt = 0
  n1st = 1
  nm2 = nn - 2
  pass2 = .false.
!
!  Loop on nodes N1:  
!  J = constraint containing N1,
!  JLAST = last node in constraint J.
!
2 continue

  j = 0
  jlast = lcc1 - 1

  do n1 = n1st, nm2

    if ( jlast < n1 ) then
!
!  N1 is the first node in constraint J+1.  Update J and
!  JLAST, and store the first constraint triangle index
!  if in pass 2.
!
      j = j + 1

      if ( j < ncc ) then
        jlast = lcc(j+1) - 1
      else
        jlast = nn
      end if

      if ( pass2 ) then
        lct(j) = kt + 1
      end if

    end if
!
!  Loop on pairs of adjacent neighbors (N2,N3).  LPLN1 points
!  to the last neighbor of N1, and LP2 points to N2.
!
    lpln1 = lend(n1)
    lp2 = lpln1

    3 continue

      lp2 = lptr(lp2)
      n2 = list(lp2)
      lp = lptr(lp2)
      n3 = abs ( list(lp) )

      if ( n2 < n1 .or. n3 < n1 ) then
        go to 10
      end if
!
!  (N1,N2,N3) is a constraint triangle iff the three nodes
!  are in the same constraint and N2 < N3.  Bypass con-
!  straint triangles on pass 1 and non-constraint triangles
!  on pass 2.
!
      cstri = n1 >= lcc1  .and.  n2 < n3  .and. n3 <= jlast

      if ( ( cstri  .and.  .not. pass2 )  .or. &
          ( .not. cstri  .and.  pass2 ) ) then
        go to 10
      end if
!
!  Add a new triangle KT = (N1,N2,N3).
!
      kt = kt + 1
      ltri(1,kt) = n1
      ltri(2,kt) = n2
      ltri(3,kt) = n3
!
!  Loop on triangle sides (I1,I2) with neighboring triangles
!  KN = (I1,I2,I3).
!
      do i = 1,3

        if ( i == 1 ) then
          i1 = n3
          i2 = n2
        else if ( i == 2 ) then
          i1 = n1
          i2 = n3
        else
          i1 = n2
          i2 = n1
        end if
!
!  Set I3 to the neighbor of I1 which follows I2 unless
!  I2->I1 is a boundary arc.
!
        lpl = lend(i1)
        lp = lptr(lpl)

4       continue

          if (list(lp) == i2) then
            go to 5
          end if

          lp = lptr(lp)

          if ( lp /= lpl ) then
            go to 4
          end if
!
!  I2 is the last neighbor of I1 unless the data structure
!  is invalid.  Bypass the search for a neighboring
!  triangle if I2->I1 is a boundary arc.
!
        if ( abs ( list(lp) ) /= i2 ) then
          go to 13
        end if

        kn = 0

        if (list(lp) < 0) then
          go to 8
        end if
!
!  I2->I1 is not a boundary arc, and LP points to I2 as
!  a neighbor of I1.
!
5   continue

        lp = lptr(lp)
        i3 = abs ( list(lp) )
!
!  Find L such that LTRI(L,KN) = I3 (not used if KN > KT),
!  and permute the vertex indexes of KN so that I1 is
!  smallest.
!
        if ( i1 < i2  .and.  i1 < i3 ) then
          l = 3
        else if (i2 < i3) then
          l = 2
          isv = i1
          i1 = i2
          i2 = i3
          i3 = isv
        else
          l = 1
          isv = i1
          i1 = i3
          i3 = i2
          i2 = isv
        end if
!
!  Test for KN > KT (triangle index not yet assigned).
!
        if ( i1 > n1  .and.  .not. pass2 ) then
          go to 9
        end if
! 
!  Find KN, if it exists, by searching the triangle list in
!  reverse order.
!
        do kn = kt-1,1,-1
          if ( ltri(1,kn) == i1  .and.  ltri(2,kn) == &
              i2 .and. ltri(3,kn) == i3 ) then
            go to 7
          end if
        end do

        go to 9
!
!  Store KT as a neighbor of KN.
!
7       continue

        ltri(l+3,kn) = kt
!
!  Store KN as a neighbor of KT, and add a new arc KA.
!
8       continue

        ltri(i+3,kt) = kn

        if (arcs) then
          ka = ka + 1
          ltri(i+6,kt) = ka
          if ( kn /= 0 ) then
            ltri(l+6,kn) = ka
          end if
        end if

9       continue

    end do
! 
!  Bottom of loop on triangles.
!
10  continue

    if ( lp2 /= lpln1 ) then
      go to 3
    end if

  end do
!
!  Bottom of loop on nodes.
!
  if ( .not. pass2 .and. 0 < ncc ) then
    pass2 = .true.
    n1st = lcc1
    go to 2
  end if
!
!  No errors encountered.
!
  nt = kt
  ier = 0
  return
!
!  Invalid triangulation data structure:  I1 is a neighbor of
!  I2, but I2 is not a neighbor of I1.
!
   13 continue

  nt = 0
  ier = 2

  return
end
subroutine trlprt ( ncc, lct, n, x, y, nrow, nt, ltri, prntx )

!*****************************************************************************80
!
!! TRLPRT prints the triangles in a triangulation.
!
!  Discussion:
!
!    Given a triangulation of a set of points in the plane,
!    this subroutine prints the triangle list created by
!    subroutine TRLIST and, optionally, the nodal coordinates
!    on logical unit LOUT.  The numbers of boundary nodes,
!    triangles, and arcs, and the constraint region triangle
!    indexes, if any, are also printed.
!
!    All parameters other than PRNTX should be
!    unaltered from their values on output from TRLIST.
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    Robert Renka,
!    Department of Computer Science,
!    University of North Texas,
!    renka@cs.unt.edu
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 751: TRIPACK, 
!    A Constrained Two-Dimensional Delaunay Triangulation Package,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 1, 1996.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NCC, the number of constraints.
!
!    Input, integer ( kind = 4 ) LCT(NCC), the list of constraint triangle starting 
!    indexes (or dummy array of length 1 if NCC = 0).
!
!    Input, integer ( kind = 4 ) N, the number of nodes in the triangulation.
!    3 <= N <= 9999.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the coordinates of the nodes in the 
!    triangulation; not used unless PRNTX = TRUE.
!
!    Input, integer ( kind = 4 ) NROW, the number of rows (entries per triangle) 
!    reserved for the triangle list LTRI.  The value must be 6 if only 
!    the vertex indexes and neighboring triangle indexes are stored, 
!    or 9 if arc indexes are also stored.
!
!    Input, integer ( kind = 4 ) NT, the number of triangles in the triangulation.
!    1 <= NT <= 9999.
!
!    Input, integer ( kind = 4 ) LTRI(NROW,NT), array whose J-th column contains
!    the vertex nodal indexes (first three rows), neighboring triangle 
!    indexes (second three rows), and, if NROW = 9, arc indexes (last
!    three rows) associated with triangle J for J = 1,...,NT.
!
!    Input, logical PRNTX, is TRUE if and only if X and Y are to be printed.
!
!  Local parameters:
!
!    I = DO-loop, nodal index, and row index for LTRI
!    K = DO-loop and triangle index
!    LUN = Logical unit number for output
!    NA = Number of triangulation arcs
!    NB = Number of boundary nodes
!    NL = Number of lines printed on the current page
!    NLMAX = Maximum number of print lines per page
!    NMAX = Maximum value of N and NT (4-digit format)
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nrow
  integer ( kind = 4 ) nt

  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) lct(*)
  integer ( kind = 4 ) ltri(nrow,nt)
  integer ( kind = 4 ) na
  integer ( kind = 4 ) nb
  integer ( kind = 4 ) nl
  integer ( kind = 4 ), parameter :: nlmax = 60
  integer ( kind = 4 ), parameter :: nmax = 9999
  integer ( kind = 4 ) ncc
  logical prntx
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)
!
!  Print a heading and test for invalid input.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TRIPACK (TRLIST) Output:'

  nl = 1
!
!  Print an error message and bypass the loops.
!
  if ( n < 3  .or.  n > nmax  .or. &
      (nrow /= 6  .and.  nrow /= 9)  .or. &
      nt < 1  .or.  nt > nmax) then
    write (*,110) n, nrow, nt
    go to 3
  end if
!
!  Print X and Y.
!
  if (prntx) then

    write (*,101)
    nl = 6

    do i = 1,n
      if (nl >= nlmax) then
        write (*,106)
        nl = 0
      end if
      write (*,102) i, x(i), y(i)
      nl = nl + 1
    end do

  end if
!
!  Print the triangulation LTRI.
!
  if ( nlmax/2 < nl ) then
    write (*,106)
    nl = 0
  end if

  if ( nrow == 6 ) then
    write (*,103)
  else
    write (*,104)
  end if

  nl = nl + 5

  do k = 1, nt
    if ( nlmax <= nl ) then
      write (*,106)
      nl = 0
    end if
    write (*,105) k, (ltri(i,k), i = 1,nrow)
    nl = nl + 1
  end do
!
!  Print NB, NA, and NT (boundary nodes, arcs, and triangles).
!
  nb = 2 * n - nt - 2
  na = nt + n - 1

  if ( nlmax-6 < nl ) then
    write (*,106)
  end if

  write (*,107) nb, na, nt
!
!  Print NCC and LCT.
!
3 continue

  write ( *, '(a,i3)' ) '  Number of constraint curves, NCC = ', ncc

  if ( 0 < ncc ) then
    write (*,109) lct(1:ncc)
  end if

  return
!
!  Print formats:
!
  101 format (//16x,'node',7x,'x(node)',10x,'y(node)'//)
  102 format (16x,i4,2e17.6)
  103 format (//1x,'triangle',8x,'vertices',12x,'neighbors'/ &
          4x,'kt',7x,'n1',5x,'n2',5x,'n3',4x,'kt1',4x, &
          'kt2',4x,'kt3'/)
  104 format (//1x,'triangle',8x,'vertices',12x,'neighbors', &
          14x,'arcs'/ &
          4x,'kt',7x,'n1',5x,'n2',5x,'n3',4x,'kt1',4x, &
          'kt2',4x,'kt3',4x,'ka1',4x,'ka2',4x,'ka3'/)
  105 format (2x,i4,2x,6(3x,i4),3(2x,i5))
  106 format (///)
  107 format (/1x,'nb = ',i4,' boundary nodes',5x, &
          'na = ',i5,' arcs',5x,'nt = ',i5, &
          ' triangles')
  109 format (1x,9x,14i5)
  110 format (//1x,10x,'*** invalid parameter:  n =',i5, &
          ', nrow =',i5,', nt =',i5,' ***')
end
subroutine trmesh ( n, x, y, list, lptr, lend, lnew, near, next, dist, ier )

!*****************************************************************************80
!
!! TRMESH triangulates a set of points in the plane.
!
!  Discussion:
!
!    This subroutine creates a Delaunay triangulation of a
!    set of N arbitrarily distributed points in the plane
!    referred to as nodes.  The Delaunay triangulation is defined
!    as a set of triangles with the following five properties:
!
!    1)  The triangle vertices are nodes.
!    2)  No triangle contains a node other than its vertices.
!    3)  The interiors of the triangles are pairwise disjoint.
!    4)  The union of triangles is the convex hull of the set
!        of nodes (the smallest convex set which contains
!        the nodes).
!    5)  The interior of the circumcircle of each triangle
!        contains no node.
!
!    The first four properties define a triangulation, and the
!    last property results in a triangulation which is as close
!    as possible to equiangular in a certain sense and which is
!    uniquely defined unless four or more nodes lie on a common
!    circle.  This property makes the triangulation well-suited
!    for solving closest point problems and for triangle-based
!    interpolation.
!
!    The triangulation can be generalized to a constrained
!    Delaunay triangulation by a call to Subroutine ADDCST.
!    This allows for user-specified boundaries defining a non-
!    convex and/or multiply connected region.
!
!    The algorithm for constructing the triangulation has
!    expected time complexity O(N*log(N)) for most nodal dis-
!    tributions.  Also, since the algorithm proceeds by adding
!    nodes incrementally, the triangulation may be updated with
!    the addition (or deletion) of a node very efficiently.
!    The adjacency information representing the triangulation
!    is stored as a linked list requiring approximately 13N
!    storage locations.
!
!
!    The following is a list of the software package modules
!    which a user may wish to call directly:
!
!    ADDCST - Generalizes the Delaunay triangulation to allow
!             for user-specified constraints.
!
!    ADDNOD - Updates the triangulation by appending or
!             inserting a new node.
!
!    AREAP  - Computes the area bounded by a closed polygonal
!             curve such as the boundary of the triangula-
!             tion or of a constraint region.
!
!    BNODES - Returns an array containing the indexes of the
!             boundary nodes in counterclockwise order.
!             Counts of boundary nodes, triangles, and arcs
!             are also returned.
!
!    CIRCUM - Computes the area, circumcenter, circumradius,
!             and, optionally, the aspect ratio of a trian-
!             gle defined by user-specified vertices.
!
!    DELARC - Deletes a boundary arc from the triangulation.
!
!    DELNOD - Updates the triangulation with the deletion of a
!             node.
!
!    EDGE   - Forces a pair of nodes to be connected by an arc
!             in the triangulation.
!
!    GETNP  - Determines the ordered sequence of L closest
!             nodes to a given node, along with the associ-
!             ated distances.  The distance between nodes is
!             taken to be the length of the shortest connec-
!             ting path which intersects no constraint
!             region.
!
!    INTSEC - Determines whether or not an arbitrary pair of
!             line segments share a common point.
!
!    JRAND  - Generates a uniformly distributed pseudo-random
!             integer ( kind = 4 ).
!
!    LEFT   - Locates a point relative to a line.
!
!    NEARND - Returns the index of the nearest node to an
!             arbitrary point, along with its squared
!             distance.
!
!    STORE  - Forces a value to be stored in main memory so
!             that the precision of floating point numbers
!             in memory locations rather than registers is
!             computed.
!
!    TRLIST - Converts the triangulation data structure to a
!             triangle list more suitable for use in a fin-
!             ite element code.
!
!    TRLPRT - Prints the triangle list created by TRLIST.
!
!    TRMESH - Creates a Delaunay triangulation of a set of nodes.
!
!    TRMSHR - Creates a Delaunay triangulation (more effici-
!             ently than TRMESH) of a set of nodes lying at
!             the vertices of a (possibly skewed) rectangu-
!             lar grid.
!
!    TRPLOT - Creates a level-2 Encapsulated Postscript (EPS)
!             file containing a triangulation plot.
!
!    TRPRNT - Prints the triangulation data structure and,
!             optionally, the nodal coordinates.
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    Robert Renka,
!    Department of Computer Science,
!    University of North Texas,
!    renka@cs.unt.edu
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 751: TRIPACK, 
!    A Constrained Two-Dimensional Delaunay Triangulation Package,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 1, 1996.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of nodes in the triangulation.  N >= 3.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the coordinates of the nodes. 
!    (X(K),Y(K)) is referred to as node K, and K is referred to as a nodal
!    index.  The first three nodes must not be collinear.
!
!    Output, integer ( kind = 4 ) LIST(6*N-12), nodal indexes which, along with LPTR,
!    LEND, and LNEW, define the triangulation as a set of N adjacency 
!    lists; counterclockwise-ordered sequences of neighboring nodes such
!    that the first and last neighbors of a boundary node are boundary 
!    nodes (the first neighbor of an interior node is arbitrary).  In
!    order to distinguish between interior and boundary nodes, the last 
!    neighbor of each boundary node is represented by the negative
!    of its index.
!
!    Output, integer ( kind = 4 ) LPTR(6*N-12), pointers (LIST indexes) in one-to-one
!    correspondence with the elements of LIST.  LIST(LPTR(I)) indexes the 
!    node which follows LIST(I) in cyclical counterclockwise order
!    (the first neighbor follows the last neighbor).
!
!    Output, integer ( kind = 4 ) LEND(N), pointers to adjacency lists.  LEND(K)
!    points to the last neighbor of node K for K = 1,...,N.  Thus, 
!    LIST(LEND(K)) < 0 if and only if K is a boundary node.
!
!    Output, integer ( kind = 4 ) LNEW, pointer to the first empty location in LIST
!    and LPTR (list length plus one).  LIST, LPTR, LEND, and LNEW are 
!    not altered if IER < 0, and are incomplete if IER > 0.
!
!    Workspace NEAR(N), NEXT(N), DIST(N).  The space is used to efficiently
!    determine the nearest triangulation node to each unprocessed node for 
!    use by ADDNOD.
!
!    Output, integer ( kind = 4 ) IER = Error indicator:
!     0 if no errors were encountered.
!    -1 if N < 3 on input.
!    -2 if the first three nodes are collinear.
!    -4 if an error flag was returned by a call to SWAP in ADDNOD.  This is 
!      an internal error and should be reported to the programmer.
!     L if nodes L and M coincide for some M > L.  The linked list represents
!      a triangulation of nodes 1 to M-1 in this case.
!
!  Local parameters:
!
!    D =        Squared distance from node K to node I
!    D1,D2,D3 = Squared distances from node K to nodes 1, 2,
!               and 3, respectively
!    EPS =      Half the machine precision
!    I,J =      Nodal indexes
!    I0 =       Index of the node preceding I in a sequence of
!               unprocessed nodes:  I = NEXT(I0)
!    K =        Index of node to be added and DO-loop index: K > 3
!    KM1 =      K-1
!    LCC(1) =   Dummy array
!    LP =       LIST index (pointer) of a neighbor of K
!    LPL =      Pointer to the last neighbor of K
!    NCC =      Number of constraint curves
!    NEXTI =    NEXT(I)
!    NN =       Local copy of N
!    SWTOL =    Tolerance for function SWPTST
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) d
  real ( kind = 8 ) d1
  real ( kind = 8 ) d2
  real ( kind = 8 ) d3
  real ( kind = 8 ) dist(n)
  real ( kind = 8 ) eps
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i0
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) km1
  integer ( kind = 4 ) lcc(1)
  logical left
  integer ( kind = 4 ) lend(n)
  integer ( kind = 4 ) list(*)
  integer ( kind = 4 ) lnew
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lpl
  integer ( kind = 4 ) lptr(*)
  integer ( kind = 4 ) ncc
  integer ( kind = 4 ) near(n)
  integer ( kind = 4 ) next(n)
  integer ( kind = 4 ) nexti
  integer ( kind = 4 ) nn
  real ( kind = 8 ) swtol
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)

  common /swpcom/ swtol

  nn = n

  if ( nn < 3 ) then
    ier = -1
    return
  end if
!
!  Compute a tolerance for function SWPTST:  SWTOL = 10*
!  (machine precision)
!
  eps = epsilon ( eps )

  swtol = eps * 20.0D+00
!
!  Store the first triangle in the linked list.
!
  if ( .not. left ( x(1), y(1), x(2), y(2), x(3), y(3) ) ) then
!
!  The initial triangle is (3,2,1) = (2,1,3) = (1,3,2).
!
    list(1) = 3
    lptr(1) = 2
    list(2) = -2
    lptr(2) = 1
    lend(1) = 2

    list(3) = 1
    lptr(3) = 4
    list(4) = -3
    lptr(4) = 3
    lend(2) = 4

    list(5) = 2
    lptr(5) = 6
    list(6) = -1
    lptr(6) = 5
    lend(3) = 6

  else if ( .not. left(x(2),y(2),x(1),y(1),x(3),y(3)) ) then
!
!  The initial triangle is (1,2,3).
!
    list(1) = 2
    lptr(1) = 2
    list(2) = -3
    lptr(2) = 1
    lend(1) = 2

    list(3) = 3
    lptr(3) = 4
    list(4) = -1
    lptr(4) = 3
    lend(2) = 4

    list(5) = 1
    lptr(5) = 6
    list(6) = -2
    lptr(6) = 5
    lend(3) = 6

  else
!
!  The first three nodes are collinear.
!
    ier = -2
    return
  end if
!
!  Initialize LNEW and test for N = 3.
!
  lnew = 7
  if (nn == 3) then
    ier = 0
    return
  end if
!
!  A nearest-node data structure (NEAR, NEXT, and DIST) is
!  used to obtain an expected-time (N*log(N)) incremental
!  algorithm by enabling constant search time for locating
!  each new node in the triangulation.
!
!  For each unprocessed node K, NEAR(K) is the index of the
!  triangulation node closest to K (used as the starting
!  point for the search in Subroutine TRFIND) and DIST(K)
!  is an increasing function of the distance between nodes
!  K and NEAR(K).
!
!  Since it is necessary to efficiently find the subset of
!  unprocessed nodes associated with each triangulation
!  node J (those that have J as their NEAR entries), the
!  subsets are stored in NEAR and NEXT as follows:  for
!  each node J in the triangulation, I = NEAR(J) is the
!  first unprocessed node in J's set (with I = 0 if the
!  set is empty), L = NEXT(I) (if I > 0) is the second,
!  NEXT(L) (if L > 0) is the third, etc.  The nodes in each
!  set are initially ordered by increasing indexes (which
!  maximizes efficiency) but that ordering is not main-
!  tained as the data structure is updated.
!
!  Initialize the data structure for the single triangle.
!
  near(1) = 0
  near(2) = 0
  near(3) = 0

  do k = nn, 4, -1

    d1 = ( x(k) - x(1) )**2 + ( y(k) - y(1) )**2
    d2 = ( x(k) - x(2) )**2 + ( y(k) - y(2) )**2
    d3 = ( x(k) - x(3) )**2 + ( y(k) - y(3) )**2

    if ( d1 <= d2  .and.  d1 <= d3 ) then
      near(k) = 1
      dist(k) = d1
      next(k) = near(1)
      near(1) = k
    else if (d2 <= d1  .and.  d2 <= d3) then
      near(k) = 2
      dist(k) = d2
      next(k) = near(2)
      near(2) = k
    else
      near(k) = 3
      dist(k) = d3
      next(k) = near(3)
      near(3) = k
    end if

  end do
!
!  Add the remaining nodes.  Parameters for ADDNOD are as follows:
!
!  K = Index of the node to be added.
!  NEAR(K) = Index of the starting node for the search in TRFIND.
!  NCC = Number of constraint curves.
!  LCC = Dummy array (since NCC = 0).
!  KM1 = Number of nodes in the triangulation.
!
  ncc = 0

  do k = 4, nn

    km1 = k-1

    call addnod ( k, x(k), y(k), near(k), ncc, lcc, km1, x, y, &
      list, lptr, lend, lnew, ier )

    if ( ier /= 0 ) then
      return
    end if
!
!  Remove K from the set of unprocessed nodes associated with NEAR(K).
!
    i = near(k)

    if (near(i) == k) then

      near(i) = next(k)

    else

      i = near(i)

      do

        i0 = i
        i = next(i0)
        if (i == k) then
          exit
        end if

      end do

      next(i0) = next(k)

    end if

    near(k) = 0
!
!  Loop on neighbors J of node K.
!
    lpl = lend(k)
    lp = lpl

4   continue

    lp = lptr(lp)
    j = abs ( list(lp) )
!
!  Loop on elements I in the sequence of unprocessed nodes
!  associated with J:  K is a candidate for replacing J
!  as the nearest triangulation node to I.  The next value
!  of I in the sequence, NEXT(I), must be saved before I
!  is moved because it is altered by adding I to K's set.
!
    i = near(j)

5   continue

    if ( i == 0 ) go to 6
    nexti = next(i)
!
!  Test for the distance from I to K less than the distance
!  from I to J.
!
    d = (x(k)-x(i))**2 + (y(k)-y(i))**2
!
!  Replace J by K as the nearest triangulation node to I:
!  update NEAR(I) and DIST(I), and remove I from J's set
!  of unprocessed nodes and add it to K's set.
!
    if ( d < dist(i) ) then
      near(i) = k
      dist(i) = d
      if (i == near(j)) then
        near(j) = nexti
      else
        next(i0) = nexti
      end if
      next(i) = near(k)
      near(k) = i
    else
      i0 = i
    end if
!
!  Bottom of loop on I.
!
    i = nexti
    go to 5
!
!  Bottom of loop on neighbors J.
!
6   continue

    if ( lp /= lpl ) then
      go to 4
    end if

  end do

  return
end
subroutine trmshr ( n, nx, x, y, nit, list, lptr, lend, lnew, ier )

!*****************************************************************************80
!
!! TRMSHR triangulates logically rectangular data.
!
!  Discussion:
!
!    This subroutine creates a Delaunay triangulation of a
!    set of N nodes in the plane, where the nodes are the vert-
!    ices of an NX by NY skewed rectangular grid with the
!    natural ordering.  Thus, N = NX*NY, and the nodes are
!    ordered from left to right beginning at the top row so
!    that adjacent nodes have indexes which differ by 1 in the
!    x-direction and by NX in the y-direction.  A skewed rec-
!    tangular grid is defined as one in which each grid cell is
!    a strictly convex quadrilateral (and is thus the convex
!    hull of its four vertices).  Equivalently, any transfor-
!    mation from a rectangle to a grid cell which is bilinear
!    in both components has an invertible Jacobian.
!
!    If the nodes are not distributed and ordered as defined
!    above, Subroutine TRMESH must be called in place of this
!    routine.  Refer to Subroutine ADDCST for the treatment of
!    constraints.
!
!    The first phase of the algorithm consists of construc-
!    ting a triangulation by choosing a diagonal arc in each
!    grid cell.  If NIT = 0, all diagonals connect lower left
!    to upper right corners and no error checking or additional
!    computation is performed.  Otherwise, each diagonal arc is
!    chosen to be locally optimal, and boundary arcs are added
!    where necessary in order to cover the convex hull of the
!    nodes.  (This is the first iteration.)  If NIT > 1 and no
!    error was detected, the triangulation is then optimized by
!    a sequence of up to NIT-1 iterations in which interior
!    arcs of the triangulation are tested and swapped if appro-
!    priate.  The algorithm terminates when an iteration
!    results in no swaps and/or when the allowable number of
!    iterations has been performed.  NIT = 0 is sufficient to
!    produce a Delaunay triangulation if the original grid is
!    actually rectangular, and NIT = 1 is sufficient if it is
!    close to rectangular.  Note, however, that the ordering
!    and distribution of nodes is not checked for validity in
!    the case NIT = 0, and the triangulation will not be valid
!    unless the rectangular grid covers the convex hull of the
!    nodes.
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    Robert Renka,
!    Department of Computer Science,
!    University of North Texas,
!    renka@cs.unt.edu
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 751: TRIPACK, 
!    A Constrained Two-Dimensional Delaunay Triangulation Package,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 1, 1996.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of nodes in the grid.  N = NX*NY for some
!    NY >= 2.
!
!    Input, integer ( kind = 4 ) NX, the number of grid points in the x-direction.
!    NX >= 2.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the coordinates of the nodes with 
!    the ordering and distribution defined in the header comments above.
!    (X(K),Y(K)) is referred to as node K.
!
!    Input/output, integer ( kind = 4 ) NIT.  On input, the maximum number of iterations
!    to be employed.  Refer to the header comments above.  On output,
!    the actual number of iterations.
!
!    Output, integer ( kind = 4 ) LIST(6*N-12), LPTR(6*N-12), LEND(N), LNEW, data structure
!    defining the triangulation.  Refer to subroutine TRMESH.
!
!    Output, integer ( kind = 4 ) IER = Error indicator:
!     0 if no errors were encountered.
!     K if the grid element with upper left corner at node K is not a strictly
!      convex quadrilateral.  The algorithm is terminated when the first such
!      occurrence is detected.  Note that this test is not performed if
!      NIT = 0 on input.
!    -1 if N, NX, or NIT is outside its valid range on input.
!    -2 if NIT > 1 on input, and the optimization loop failed to converge
!      within the allowable number of iterations.  The triangulation is
!      valid but not optimal in this case.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) eps
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) iter
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kp1
  logical left
  integer ( kind = 4 ) lend(n)
  integer ( kind = 4 ) list(*)
  integer ( kind = 4 ) lnew
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lpf
  integer ( kind = 4 ) lpk
  integer ( kind = 4 ) lpl
  integer ( kind = 4 ) lpp
  integer ( kind = 4 ) lptr(*)
  integer ( kind = 4 ) lstptr
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) m3
  integer ( kind = 4 ) m4
  integer ( kind = 4 ) maxit
  integer ( kind = 4 ) n0
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n3
  integer ( kind = 4 ) n4
  integer ( kind = 4 ) nbcnt
  integer ( kind = 4 ) ni
  integer ( kind = 4 ) nit
  integer ( kind = 4 ) nj
  integer ( kind = 4 ) nm1
  integer ( kind = 4 ) nn
  integer ( kind = 4 ) nnb
  integer ( kind = 4 ) nx
  logical swptst
  real ( kind = 8 ) swtol
  logical tst
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)

  common /swpcom/ swtol
!
!  Store local variables and test for errors in input parameters.
!
  ni = nx
  nj = n / ni
  nn = ni * nj
  maxit = nit
  nit = 0

  if ( n /= nn .or. nj < 2 .or. ni < 2 .or. maxit < 0 ) then
    ier = -1
    return
  end if

  ier = 0
!
!  Compute a tolerance for function SWPTST.
!
  eps = epsilon ( eps )
  swtol = eps * 20.0D+00
!
!  Loop on grid points (I,J) corresponding to nodes K =
!  (J-1)*NI + I.  TST = TRUE iff diagonals are to be
!  chosen by the swap test.  M1, M2, M3, and M4 are the
!  slopes (-1, 0, or 1) of the diagonals in quadrants 1
!  to 4 (counterclockwise beginning with the upper right)
!  for a coordinate system with origin at node K.
!
  tst = maxit > 0
  m1 = 0
  m4 = 0
  lp = 0
  kp1 = 1

  do 6 j = 1,nj

    do 5 i = 1,ni

      m2 = m1
      m3 = m4
      k = kp1
      kp1 = k + 1
      lpf = lp + 1

      if ( j == nj .and. i /= ni ) go to 2

      if ( i /= 1 ) then

        if ( j /= 1 ) then
!
!  K is not in the top row, leftmost column, or bottom row
!  (unless K is the lower right corner).  Take the first
!  neighbor to be the node above K.
!
          lp = lp + 1
          list(lp) = k - ni
          lptr(lp) = lp + 1

          if ( m2 <= 0 ) then
            lp = lp + 1
            list(lp) = k - 1 - ni
            lptr(lp) = lp + 1
          end if

        end if
!
!  K is not in the leftmost column.  The next (or first)
!  neighbor is to the left of K.
!
        lp = lp + 1
        list(lp) = k - 1
        lptr(lp) = lp + 1
        if (j == nj) go to 3

        if (m3 >= 0) then
          lp = lp + 1
          list(lp) = k - 1 + ni
          lptr(lp) = lp + 1
        end if

      end if
!
!  K is not in the bottom row.  The next (or first) neighbor is below K.
!
      lp = lp + 1
      list(lp) = k + ni
      lptr(lp) = lp + 1
!
!  Test for a negative diagonal in quadrant 4 unless K is
!  in the rightmost column.  The quadrilateral associated
!  with the quadrant is tested for strict convexity un-
!  less NIT = 0 on input.
!
      if ( i == ni ) then
        go to 3
      end if

      m4 = 1
      if ( .not. tst ) go to 2

      if ( left(x(kp1),y(kp1),x(k+ni),y(k+ni),x(k),y(k)) .or. &
           left(x(k),y(k),x(kp1+ni),y(kp1+ni),x(k+ni),y(k+ni)) .or. &
           left(x(k+ni),y(k+ni),x(kp1),y(kp1), x(kp1+ni),y(kp1+ni)) .or. &
           left(x(kp1+ni),y(kp1+ni),x(k),y(k), x(kp1),y(kp1)) ) then
        ier = k
        return
      end if

      if ( swptst ( kp1, k+ni, k, kp1+ni, x, y ) ) go to 2

      m4 = -1
      lp = lp + 1
      list(lp) = kp1 + ni
      lptr(lp) = lp + 1
!
!  The next (or first) neighbor is to the right of K.
!
2     continue

      lp = lp + 1
      list(lp) = kp1
      lptr(lp) = lp + 1
!
!  Test for a positive diagonal in quadrant 1 (the neighbor
!  of K-NI which follows K is not K+1) unless K is in the
!  top row.
!
      if (j == 1) go to 3

      if (tst) then
        m1 = -1
        lpk = lstptr(lend(k-ni),k,list,lptr)
        lpk = lptr(lpk)

        if ( list(lpk) /= kp1 ) then
          m1 = 1
          lp = lp + 1
          list(lp) = kp1 - ni
          lptr(lp) = lp + 1
        end if

      end if
!
!  If K is in the leftmost column (and not the top row) or
!  in the bottom row (and not the rightmost column), then
!  the next neighbor is the node above K.
!
      if ( i /= 1 .and. j /= nj ) go to 4

      lp = lp + 1
      list(lp) = k - ni
      lptr(lp) = lp + 1
      if ( i == 1 ) go to 3
!
!  K is on the bottom row (and not the leftmost or rightmost column).
!
      if ( m2 <= 0 ) then
        lp = lp + 1
        list(lp) = k - 1 - ni
        lptr(lp) = lp + 1
      end if

      lp = lp + 1
      list(lp) = k - 1
      lptr(lp) = lp + 1
!
!  K is a boundary node.
!
3     continue

      list(lp) = -list(lp)
!
!  Bottom of loop.  Store LEND and correct LPTR(LP).
!  LPF and LP point to the first and last neighbors of K.
!
4     continue

      lend(k) = lp
      lptr(lp) = lpf

5     continue
6   continue
!
!  Store LNEW, and terminate the algorithm if NIT = 0 on input.
!
  lnew = lp + 1

  if ( maxit == 0 ) then
    return
  end if
!
!  Add boundary arcs where necessary in order to cover the
!  convex hull of the nodes.  N1, N2, and N3 are consecu-
!  tive boundary nodes in counterclockwise order, and N0
!  is the starting point for each loop around the boundary.
!
  n0 = 1
  n1 = n0
  n2 = ni + 1
!
!  TST is set to TRUE if an arc is added.  The boundary
!  loop is repeated until a traversal results in no
!  added arcs.
!
7 continue

  tst = .false.
!
!  Top of boundary loop.  Set N3 to the first neighbor of
!  N2, and test for N3 LEFT N1 -> N2.
!
8     continue

      lpl = lend(n2)

      lp = lptr(lpl)
      n3 = list(lp)

      if ( left(x(n1),y(n1),x(n2),y(n2),x(n3),y(n3)) ) then
         n1 = n2
      end if

      if (n1 /= n2) then
!
!  Add the boundary arc N1-N3.  If N0 = N2, the starting
!  point is changed to N3, since N2 will be removed from
!  the boundary.  N3 is inserted as the first neighbor of
!  N1, N2 is changed to an interior node, and N1 is
!  inserted as the last neighbor of N3.
!
        tst = .true.
        if (n2 == n0) n0 = n3
        lp = lend(n1)
        call insert (n3,lp, list,lptr,lnew )
        list(lpl) = -list(lpl)
        lp = lend(n3)
        list(lp) = n2
        call insert (-n1,lp, list,lptr,lnew )
        lend(n3) = lnew - 1
      end if
!
!  Bottom of loops.  Test for termination.
!
      n2 = n3
      if (n1 /= n0) go to 8
    if (tst) go to 7
!
!  Terminate the algorithm if NIT = 1 on input.
!
  nit = 1

  if ( maxit == 1 ) then
    return
  end if
!
!  Optimize the triangulation by applying the swap test and
!  appropriate swaps to the interior arcs.  The loop is
!  repeated until no swaps are performed or MAXIT itera-
!  tions have been applied.  ITER is the current iteration,
!  and TST is set to TRUE if a swap occurs.
!
  iter = 1
  nm1 = nn - 1

9 continue

  iter = iter + 1
    tst = .false.
!
!  Loop on interior arcs N1-N2, where N2 > N1 and
!  (N1,N2,N3) and (N2,N1,N4) are adjacent triangles.
!
!  Top of loop on nodes N1.
!
    do 11 n1 = 1,nm1

      lpl = lend(n1)
      n4 = list(lpl)
      lpf = lptr(lpl)
      n2 = list(lpf)
      lp = lptr(lpf)
      n3 = list(lp)
      nnb = nbcnt(lpl,lptr)
!
!  Top of loop on neighbors N2 of N1.  NNB is the number of
!  neighbors of N1.
!
      do i = 1,nnb
!
!  Bypass the swap test if N1 is a boundary node and N2 is
!  the first neighbor (N4 < 0), N2 < N1, or N1-N2 is a
!  diagonal arc (already locally optimal) when ITER = 2.
!
        if ( n4 > 0  .and.  n2 > n1  .and. &
            ( iter /= 2  .or.  abs ( n1+ni-n2 ) /= 1 ) ) then

          if (swptst(n3,n4,n1,n2,x,y) ) then
!
!  Swap diagonal N1-N2 for N3-N4, set TST to TRUE, and set
!  N2 to N4 (the neighbor preceding N3).
!
            call swap (n3,n4,n1,n2, list,lptr,lend, lpp)
            if (lpp /= 0) then
              tst = .true.
              n2 = n4
            end if
          end if
        end if
!
!  Bottom of neighbor loop.
!
        if (list(lpl) == -n3) then
          go to 11
        end if

        n4 = n2
        n2 = n3
        lp = lstptr(lpl,n2,list,lptr)
        lp = lptr(lp)
        n3 = abs ( list(lp) )

      end do

11     continue
!
!  Test for termination.
!
  if (tst  .and.  iter < maxit) go to 9

  nit = iter

  if ( tst ) then
    ier = -2
  end if

  return
end
subroutine trmtst ( n, x, y, list, lptr, lend, lnew, tol, armax, ier )

!*****************************************************************************80
!
!! TRMTST tests a data structure representing a Delaunay triangulation.
!
!  Discussion:
!
!    This subroutine tests the validity of the data structure
!    representing a Delaunay triangulation created by subrou-
!    tine TRMESH.  The following properties are tested:
!
!    1)  Each interior node has at least three neighbors, and
!        each boundary node has at least two neighbors.
!
!    2)  abs ( LIST(LP) ) is a valid nodal index in the range
!        1 to N and LIST(LP) > 0 unless LP = LEND(K) for some
!        nodal index K.
!
!    3)  Each pointer LEND(K) for K = 1 to N and LPTR(LP) for
!        LP = 1 to LNEW-1 is a valid LIST index in the range
!        1 to LNEW-1.
!
!    4)  N .GE. NB .GE. 3, NT = 2*N-NB-2, and NA = 3*N-NB-3 =
!        (LNEW-1)/2, where NB, NT, and NA are the numbers of
!        boundary nodes, triangles, and arcs, respectively.
!
!    5)  Each circumcircle defined by the vertices of a tri-
!        angle contains no nodes in its interior.  This prop-
!        erty distinguishes a Delaunay triangulation from an
!        arbitrary triangulation of the nodes.
!
!    Note that no test is made for the property that a triangulation 
!    covers the convex hull of the nodes.
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    Robert Renka,
!    Department of Computer Science,
!    University of North Texas,
!    renka@cs.unt.edu
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 751: TRIPACK, 
!    A Constrained Two-Dimensional Delaunay Triangulation Package,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 1, 1996.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of nodes.  N .GE. 3.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the nodal coordinates.
!
!    Input, integer ( kind = 4 ) LIST(*), LPTR(*), LEND(N), the data structure 
!    containing the triangulation.  Refer to subroutine TRMESH.
!
!    Input, real ( kind = 8 ) TOL, nonnegative tolerance to allow for
!    floating-point errors in the circumcircle test.  An error situation
!    is defined as 
!      (R**2 - D**2) / R**2 > TOL, 
!    where R is the radius of a circumcircle and D is the distance from the
!    circumcenter to the nearest node.  A reasonable value for TOL is 
!    10*EPS, where EPS is the machine precision.  The test is effectively
!    bypassed by making TOL large.  If TOL < 0, the tolerance is taken 
!    to be 0.
!
!    Output, real ( kind = 8 ) ARMAX, maximum aspect ratio (radius of inscribed
!    circle divided by circumradius) of a triangle in the triangulation 
!    unless 0 < IER.
!
!    Output, integer ( kind = 4 ) IER = Error indicator:
!    -1 if one or more null triangles (area = 0) are present but no (other)
!      errors were encountered.  A null triangle is an error only if it 
!      occurs in the the interior.
!     0 if no errors or null triangles were encountered.
!     1 if a node has too few neighbors.
!     2 if a LIST entry is outside its valid range.
!     3 if a LPTR or LEND entry is outside its valid range.
!     4 if the triangulation parameters (N, NB, NT, NA, and LNEW) are
!      inconsistent (or N < 3 or LNEW is invalid).
!     5 if a triangle contains a node interior to its circumcircle.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) ar
  real ( kind = 8 ) armax
  real ( kind = 8 ) cr
  real ( kind = 8 ) cx
  real ( kind = 8 ) cy
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) k
  integer ( kind = 4 ) lend(n)
  integer ( kind = 4 ) list(*)
  integer ( kind = 4 ) lnew
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lpl
  integer ( kind = 4 ) lpmax
  integer ( kind = 4 ) lpn
  integer ( kind = 4 ) lptr(*)
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n3
  integer ( kind = 4 ) na
  integer ( kind = 4 ) nb
  integer ( kind = 4 ) nfail
  integer ( kind = 4 ) nn
  integer ( kind = 4 ) nnb
  integer ( kind = 4 ) nt
  integer ( kind = 4 ) null
  logical ratio
  real ( kind = 8 ) rs
  real ( kind = 8 ) rtol
  real ( kind = 8 ) sa
  real ( kind = 8 ) tol
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)
!
!  Store local variables, test for errors in input, and
!  initialize counts.
!
  nn = n
  lpmax = lnew - 1
  rtol = tol
  rtol = max ( rtol, 0.0D+00 )
  ratio = .true.
  armax = 0.0D+00

  if ( nn < 3 ) then
    go to 14
  end if

  nb = 0
  nt = 0
  null = 0
  nfail = 0
!
!  Loop on triangles (N1,N2,N3) such that N2 and N3 index
!  adjacent neighbors of N1 and are both larger than N1
!  (each triangle is associated with its smallest index).
!  NNB is the neighbor count for N1.
!
  do n1 = 1,nn

    nnb = 0
    lpl = lend(n1)

    if ( lpl < 1  .or.  lpmax < lpl ) then
      lp = lpl
      go to 13
    end if

    lp = lpl
!
!  Loop on neighbors of N1.
!
1   continue

      lp = lptr(lp)
      nnb = nnb + 1

      if (lp < 1  .or.  lp > lpmax) then
        go to 13
      end if

      n2 = list(lp)

      if (n2 < 0) then
        if (lp /= lpl) then
          go to 12
        end if
        if (n2 == 0  .or.  -n2 > nn) go to 12
        nb = nb + 1
        go to 4
      end if

      if (n2 < 1  .or.  n2 > nn) then
        go to 12
      end if

      lpn = lptr(lp)
      n3 = abs ( list(lpn) )
      if (n2 < n1  .or.  n3 < n1) go to 4
      nt = nt + 1
!
!  Compute the coordinates of the circumcenter of (N1,N2,N3).
!
      call circum (x(n1),y(n1),x(n2),y(n2),x(n3),y(n3), &
                   ratio, cx,cy,cr,sa,ar)
      if (sa == 0.) then
        null = null + 1
        go to 4
      end if

      armax = max(armax,ar)
!
!  Test for nodes within the circumcircle.
!
      rs = cr*cr*(1.-rtol)

      do k = 1,nn
        if ( k == n1  .or.  k == n2  .or. &
            k == n3 ) go to 2
        if ((cx-x(k))**2 + (cy-y(k))**2 < rs) go to 3
2       continue
      end do

      go to 4
!
!  Node K is interior to the circumcircle of (N1,N2,N3).
!
3     continue

      nfail = nfail + 1
!
!  Bottom of loop on neighbors.
!
4     continue

      if (lp /= lpl) go to 1
    if (nnb < 2  .or.  (nnb == 2  .and. &
        list(lpl) > 0)) go to 11

  end do
!
!  Test parameters for consistency and check for NFAIL = 0.
!
  na = lpmax/2
  if (nb < 3  .or.  nt /= 2 * nn - nb - 2  .or. &
      na /= 3*nn-nb-3) go to 14
  if (nfail /= 0) go to 15
!
!  No errors were encountered.
!
  ier = 0
  if (null == 0) return
  ier = -1
  write (*,100) null
  100 format (//5x,'*** trmtst -- ',i5,' null triangles ', &
          'are present'/19x,'(null triangles ', &
          'on the boundary are unavoidable) ***'//)
  return
!
!  Node N1 has fewer than three neighbors.
!
11 continue

  ier = 1
  write (*,110) n1, nnb
  110 format (//5x,'*** trmtst -- node ',i5, &
          ' has only ',i5,' neighbors ***'/)
  return
!
!  N2 = LIST(LP) is outside its valid range.
!
12 continue

  ier = 2
  write (*,120) n2, lp, n1
  120 format (//5x,'*** trmtst -- list(lp) =',i5, &
          ', for lp =',i5,','/19x, &
  'is not a valid neighbor of ',i5,' ***'/)
  return
!
!  LIST pointer LP is outside its valid range.
!
13 continue

  ier = 3
  write (*,130) lp, lnew, n1
  130 format (//5x,'*** trmtst -- lp =',i5,' is not in the', &
          ' range 1 to lnew-1 for lnew = ',i5/ &
          19x,'lp points to a neighbor of ',i5, &
          ' ***'/)
  return
!
!  Inconsistent triangulation parameters encountered.
!
14 continue

  ier = 4
  write (*,140) n, lnew, nb, nt, na
  140 format (//5x,'*** trmtst -- inconsistent parameters', &
          ' ***'/19x,'n = ',i5,' nodes',12x,'lnew =',i5/ &
          19x,'nb = ',i5,' boundary nodes'/ &
          19x,'nt = ',i5,' triangles'/ &
          19x,'na = ',i5,' arcs'/)
  return
!
!  Circumcircle test failure.
!
15 continue

  ier = 5
  write (*,150) nfail
  150 format (//5x,'*** trmtst -- ',i5,' circumcircles ', &
          'contain nodes in their interiors ***'/)
  return
end
subroutine trplot ( lun, pltsiz, wx1, wx2, wy1, wy2, ncc, lcc, &
  n, x, y, list, lptr, lend, title, numbr, ier )

!*****************************************************************************80
!
!! TRPLOT plots a triangulation in an EPS file.
!
!  Discussion:
!
!    This subroutine creates a level-2 Encapsulated Postscript (EPS) file 
!    containing a triangulation plot.
!
!    Various plotting options can be controlled by altering
!    the data statement below.
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    Robert Renka,
!    Department of Computer Science,
!    University of North Texas,
!    renka@cs.unt.edu
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 751: TRIPACK, 
!    A Constrained Two-Dimensional Delaunay Triangulation Package,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 1, 1996.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LUN, the logical unit number in the range 0 to 99.
!    The unit should be opened with an appropriate
!    file name before the call to this routine.
!
!    Input, real ( kind = 8 ) PLTSIZ, the plot size in inches.  The window is
!    mapped, with aspect ratio preserved, to a rectangular viewport with maximum
!    side-length equal to .88*PLTSIZ (leaving room for labels outside 
!    the viewport).  The viewport is centered on the 8.5 by 11 inch page, 
!    and its boundary is drawn.  1.0 <= PLTSIZ <= 8.5.
!
!    Input, real ( kind = 8 ) WX1, WX2, WY1, WY2, parameters defining a
!    rectangular window against which the triangulation is clipped.  (Only the
!    portion of the triangulation that lies in the window is drawn.)
!    (WX1,WY1) and (WX2,WY2) are the lower left and upper right 
!    corners, respectively.  WX1 < WX2 and WY1 < WY2.
!
!    Input, integer ( kind = 4 ) NCC, the number of constraint curves.  Refer to 
!    subroutine ADDCST.  NCC >= 0.
!
!    Input, integer ( kind = 4 ) LCC(NCC) (or dummy parameter if NCC = 0) containing the
!    index of the first node of constraint I in LCC(I).  For I = 1 to
!    NCC, LCC(I+1)-LCC(I) >= 3, where LCC(NCC+1) = N+1.
!
!    Input, integer ( kind = 4 ) N, the number of nodes in the triangulation.  3 <= N.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the coordinates of the nodes with
!    non-constraint nodes in the first LCC(1)-1 locations.
!
!    Input, integer ( kind = 4 ) LIST(*), LPTR(*), LEND(N), data structure defining the
!    triangulation.  Refer to subroutine TRMESH.
!
!    Input, character ( len = * ) TITLE, a string to be centered above the 
!    plot.  The string must be enclosed in parentheses; i.e., the first and 
!    last characters must be '(' and ')', respectively, but these are not
!    displayed.  TITLE may have at most 80 characters including the parentheses.
!
!    Input, logical NUMBR, option indicator:  If NUMBR = TRUE, the
!    nodal indexes are plotted next to the nodes.
!
!    Output, integer ( kind = 4 ) IER = Error indicator:
!    0, if no errors were encountered.
!    1, if LUN, PLTSIZ, NCC, or N is outside its valid range.  LCC is not 
!      tested for validity.
!    2, if WX1 >= WX2 or WY1 >= WY2.
!    3, if an error was encountered in writing to unit LUN.
!
!  Local parameters:
!
!    ANNOT =     Logical variable with value TRUE iff the plot
!                is to be annotated with the values of WX1,
!                WX2, WY1, and WY2
!    CNSTR       Logical variable used to flag constraint arcs:
!                TRUE iff N0-N1 lies in a constraint region
!    DASHL =     Length (in points, at 72 points per inch) of
!                dashes and spaces in a dashed line pattern
!                used for drawing constraint arcs
!    DX =        Window width WX2-WX1
!    DY =        Window height WY2-WY1
!    FSIZN =     Font size in points for labeling nodes with
!                their indexes if NUMBR = TRUE
!    FSIZT =     Font size in points for the title (and
!                annotation if ANNOT = TRUE)
!    I =         Constraint index (1 to NCC)
!    IFRST =     Index of the first node in constraint I
!    IH =        Height of the viewport in points
!    ILAST =     Index of the last node in constraint I
!    IPX1,IPY1 = X and y coordinates (in points) of the lower
!                left corner of the bounding box or viewport
!    IPX2,IPY2 = X and y coordinates (in points) of the upper
!                right corner of the bounding box or viewport
!    IW =        Width of the viewport in points
!    LP =        LIST index (pointer)
!    LPL =       Pointer to the last neighbor of N0
!    N0 =        Nodal index and DO-loop index
!    N0BAK =     Predecessor of N0 in a constraint curve
!                (sequence of adjacent constraint nodes)
!    N0FOR =     Successor to N0 in a constraint curve
!    N1 =        Index of a neighbor of N0
!    NLS =       Index of the last non-constraint node
!    PASS1 =     Logical variable used to flag the first pass
!                through the constraint nodes
!    R =         Aspect ratio DX/DY
!    SFX,SFY =   Scale factors for mapping world coordinates
!                (window coordinates in [WX1,WX2] X [WY1,WY2])
!                to viewport coordinates in [IPX1,IPX2] X [IPY1,IPY2]
!    T =         Temporary variable
!    TX,TY =     Translation vector for mapping world coordi-
!                nates to viewport coordinates
!    X0,Y0 =     X(N0),Y(N0) or label location
!
  implicit none

  integer ( kind = 4 ) n

  logical, parameter :: annot = .true.
  logical cnstr
  real ( kind = 8 ), parameter :: dashl = 4.0D+00
  real ( kind = 8 ) dx
  real ( kind = 8 ) dy
  real ( kind = 8 ), parameter :: fsizn = 10.0D+00
  real ( kind = 8 ), parameter :: fsizt = 16.0D+00
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) ifrst
  integer ( kind = 4 ) ih
  integer ( kind = 4 ) ilast
  integer ( kind = 4 ) ipx1
  integer ( kind = 4 ) ipx2
  integer ( kind = 4 ) ipy1
  integer ( kind = 4 ) ipy2
  integer ( kind = 4 ) iw
  integer ( kind = 4 ) lcc(*)
  integer ( kind = 4 ) lend(n)
  integer ( kind = 4 ) list(*)
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lpl
  integer ( kind = 4 ) lptr(*)
  integer ( kind = 4 ) lun
  integer ( kind = 4 ) n0
  integer ( kind = 4 ) n0bak
  integer ( kind = 4 ) n0for
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) ncc
  integer ( kind = 4 ) nls
  logical numbr
  logical pass1
  real ( kind = 8 ) pltsiz
  real ( kind = 8 ) r
  real ( kind = 8 ) sfx
  real ( kind = 8 ) sfy
  real ( kind = 8 ) t
  character ( len = * ) title
  real ( kind = 8 ) tx
  real ( kind = 8 ) ty
  real ( kind = 8 ) wx1
  real ( kind = 8 ) wx2
  real ( kind = 8 ) wy1
  real ( kind = 8 ) wy2
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x0
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) y0
!
!  Test for error 1.
!
  if (lun < 0  .or.  lun > 99  .or. &
      pltsiz < 1.0D+00  .or.  pltsiz > 8.5  .or. &
      ncc < 0  .or.  n < 3) then
    ier = 1
    return
  end if
!
!  Set NLS to the last non-constraint node.
!
  if ( ncc > 0 ) then
    nls = lcc(1)-1
  else
    nls = n
  end if
!
!  Compute the aspect ratio of the window.
!
  dx = wx2 - wx1
  dy = wy2 - wy1

  if ( dx <= 0.0D+00 .or. dy <= 0.0D+00 ) then
    ier = 2
    return
  end if

  r = dx / dy
!
!  Compute the lower left (IPX1,IPY1) and upper right
!  (IPX2,IPY2) corner coordinates of the bounding box.
!  The coordinates, specified in default user space units
!  (points, at 72 points/inch with origin at the lower
!  left corner of the page), are chosen to preserve the
!  aspect ratio R, and to center the plot on the 8.5 by 11
!  inch page.  The center of the page is (306,396), and
!  T = PLTSIZ/2 in points.
!
  t = 36.0D+00 * pltsiz

  if ( 1.0D+00 <= r ) then
    ipx1 = 306 - nint(t)
    ipx2 = 306 + nint(t)
    ipy1 = 396 - nint(t/r)
    ipy2 = 396 + nint(t/r)
  else
    ipx1 = 306 - nint(t*r)
    ipx2 = 306 + nint(t*r)
    ipy1 = 396 - nint(t)
    ipy2 = 396 + nint(t)
  end if
!
!  Output header comments.
!
  write (lun,100,err=13) ipx1, ipy1, ipx2, ipy2
  100 format ('%!ps-adobe-3.0 epsf-3.0'/ &
          '%%boundingbox:',4i4/ &
          '%%title:  triangulation'/ &
          '%%creator:  tripack'/ &
          '%%endcomments')
!
!  Set (IPX1,IPY1) and (IPX2,IPY2) to the corner coordinates
!  of a viewport obtained by shrinking the bounding box by
!  12% in each dimension.
!
  iw = nint ( 0.88D+00 * real ( ipx2 - ipx1, kind = 8 ) )
  ih = nint ( 0.88D+00 * real ( ipy2 - ipy1, kind = 8 ) )
  ipx1 = 306 - iw/2
  ipx2 = 306 + iw/2
  ipy1 = 396 - ih/2
  ipy2 = 396 + ih/2
!
!  Set the line thickness to 2 points, and draw the viewport boundary.
!
  t = 2.0D+00
  write (lun,110,err=13) t
  write (lun,120,err=13) ipx1, ipy1
  write (lun,130,err=13) ipx1, ipy2
  write (lun,130,err=13) ipx2, ipy2
  write (lun,130,err=13) ipx2, ipy1
  write (lun,140,err=13)
  write (lun,150,err=13)
  110 format (f12.6,' setlinewidth')
  120 format (2i4,' moveto')
  130 format (2i4,' lineto')
  140 format ('closepath')
  150 format ('stroke')
!
!  Set up a mapping from the window to the viewport.
!
  sfx = real ( iw, kind = 8 ) / dx
  sfy = real ( ih, kind = 8 ) / dy
  tx = ipx1 - sfx*wx1
  ty = ipy1 - sfy*wy1
  write (lun,160,err=13) tx, ty, sfx, sfy
  160 format (2f12.6,' translate'/ &
          2f12.6,' scale')
!
!  The line thickness (believe it or not) must be
!  changed to reflect the new scaling which is applied to
!  all subsequent output.  Set it to 1.0 point.
!
  t = 2.0D+00 / (sfx+sfy)
  write (lun,110,err=13) t
!
!  Save the current graphics state, and set the clip path to
!  the boundary of the window.
!
  write (lun,170,err=13)
  write (lun,180,err=13) wx1, wy1
  write (lun,190,err=13) wx2, wy1
  write (lun,190,err=13) wx2, wy2
  write (lun,190,err=13) wx1, wy2
  write (lun,200,err=13)
  170 format ('gsave')
  180 format (2f12.6,' moveto')
  190 format (2f12.6,' lineto')
  200 format ('closepath clip newpath')
!
!  Draw the edges N0->N1, where N1 > N0, beginning with a
!  loop on non-constraint nodes N0.  LPL points to the
!  last neighbor of N0.
!
  do n0 = 1,nls

    x0 = x(n0)
    y0 = y(n0)
    lpl = lend(n0)
    lp = lpl
!
!  Loop on neighbors N1 of N0.
!
    do

      lp = lptr(lp)
      n1 = abs ( list(lp) )

      if ( n0 < n1 ) then
        write (lun,210,err=13) x0, y0, x(n1), y(n1)
  210       format (2f12.6,' moveto',2f12.6,' lineto')
      end if

      if ( lp == lpl) then
        exit
      end if

    end do

  end do
!
!  Loop through the constraint nodes twice.  The non-constraint arcs 
!  incident on constraint nodes are drawn (with solid lines) on the first 
!  pass, and the constraint arcs (both boundary and interior, if any)
!  are drawn (with dashed lines) on the second pass.
!
  pass1 = .true.
!
!  Loop on constraint nodes N0 with (N0BAK,N0,N0FOR) a subsequence of 
!  constraint I.  The outer loop is on constraints I with first and last 
!  nodes IFRST and ILAST.
!
4 continue

  ifrst = n+1

  do i = ncc, 1, -1

    ilast = ifrst - 1
    ifrst = lcc(i)
    n0bak = ilast

    do n0 = ifrst, ilast

      n0for = n0 + 1
      if (n0 == ilast) n0for = ifrst
      lpl = lend(n0)
      x0 = x(n0)
      y0 = y(n0)
      lp = lpl
!
!  Loop on neighbors N1 of N0.  CNSTR = TRUE iff N0-N1 is a
!  constraint arc.
!
!  Initialize CNSTR to TRUE iff the first neighbor of N0
!  strictly follows N0FOR and precedes or coincides with
!  N0BAK (in counterclockwise order).
!
      do

        lp = lptr(lp)
        n1 = abs ( list(lp) )
        if ( n1 == n0for .or. n1 == n0bak ) then
          exit
        end if

      end do

      cnstr = n1 == n0bak
      lp = lpl
!
!  Loop on neighbors N1 of N0.  Update CNSTR and test for N1 > N0.
!
6 continue

        lp = lptr(lp)
        n1 = abs ( list(lp) )

        if (n1 == n0for) then
          cnstr = .true.
        end if
!
!  Draw the edge iff (PASS1=TRUE and CNSTR=FALSE) or
!  (PASS1=FALSE and CNSTR=TRUE); i.e., CNSTR and PASS1
!  have opposite values.
!
        if ( n0 < n1 ) then

          if ( cnstr .neqv. pass1 ) then
            write (lun,210,err=13) x0, y0, x(n1), y(n1)
          end if

        end if

        if ( n1 == n0bak ) then
          cnstr = .false.
        end if
!
!  Bottom of loops.
!
        if ( lp /= lpl ) then
          go to 6
        end if

      n0bak = n0

    end do

  end do

  if (pass1) then
!
!  End of first pass:  paint the path and change to dashed
!  lines for subsequent drawing.  Since the scale factors
!  are applied to everything, the dash length must be
!  specified in world coordinates.
!
    pass1 = .false.
    write (lun,150,err=13)
    t = dashl * 2.0D+00 / ( sfx + sfy )
    write (lun,220,err=13) t
  220   format ('[',f12.6,'] 0 setdash')
    go to 4

  end if
!
!  Paint the path and restore the saved graphics state (with
!  no clip path).
!
  write (lun,150,err=13)
  write (lun,230,err=13)
  230 format ('grestore')

  if (numbr) then
!
!  Nodes in the window are to be labeled with their indexes.
!  Convert FSIZN from points to world coordinates, and
!  output the commands to select a font and scale it.
!
    t = fsizn * 2.0D+00 / ( sfx + sfy )
    write (lun,240,err=13) t
  240   format ('/Helvetica findfont'/ &
            f12.6,' scalefont setfont')
!
!  Loop on nodes N0 with coordinates (X0,Y0).
!
    do n0 = 1, n

      x0 = x(n0)
      y0 = y(n0)
!
!  Move to (X0,Y0), and draw the label N0.  The first character will 
!  have its lower left corner about one
!  character width to the right of the nodal position.
!
      if ( x0 >= wx1  .and.  x0 <= wx2  .and. &
           y0 >= wy1  .and.  y0 <= wy2 ) then
        write (lun,180,err=13) x0, y0
        write (lun,250,err=13) n0
  250   format ('(',i3,') show')
      end if

    end do

  end if
!
!  Convert FSIZT from points to world coordinates, and output
!  the commands to select a font and scale it.
!
  t = fsizt * 2.0D+00 / ( sfx + sfy )
  write (lun,240,err=13) t
!
!  Display TITLE centered above the plot:
!
  y0 = wy2 + 3.0D+00 * t
  write (lun,260,err=13) title, ( wx1 + wx2 ) / 2.0D+00, y0
  260 format (a80/'  stringwidth pop 2 div neg ',f12.6, &
          ' add ',f12.6,' moveto')
  write (lun,270,err=13) title
  270 format (a80/'  show')
  if (annot) then
!
!  Display the window extrema below the plot.
!
    x0 = wx1
    y0 = wy1 - 100.0D+00 / ( sfx + sfy )
    write (lun,180,err=13) x0, y0
    write (lun,280,err=13) wx1, wx2
    y0 = y0 - 2.0D+00 * t
    write (lun,290,err=13) x0, y0, wy1, wy2
  280   format ('(window:   wx1 = ',e9.3,',   wx2 = ',e9.3, &
            ') show')
  290   format ('(window:  ) stringwidth pop ',f12.6,' add', &
            f12.6,' moveto'/ &
            '( wy1 = ',e9.3,',   wy2 = ',e9.3,') show')
  end if
!
!  Paint the path and output the showpage command and
!  end-of-file indicator.
!
  write (lun,300,err=13)
  300 format ('stroke'/ &
          'showpage'/ &
          '%%eof')
!
!  HP's interpreters require a one-byte End-of-PostScript-Job
!  indicator (to eliminate a timeout error message): ASCII 4.
!
  write (lun,310,err=13) char(4)
  310 format (a1)
!
!  No error encountered.
!
  ier = 0
  return
!
!  Error writing to unit LUN.
!
   13 ier = 3
  return
end
subroutine trprnt ( ncc, lcc, n, x, y, list, lptr, lend, prntx )

!*****************************************************************************80
!
!! TRPRNT prints information about a planar triangulation.
!
!  Discussion:
!
!    Given a triangulation of a set of points in the plane,
!    this subroutine prints the adjacency lists and, optionally, 
!    the nodal coordinates.  The list of neighbors of a boundary 
!    node is followed by index 0.  The numbers of boundary nodes, 
!    triangles, and arcs, and the constraint curve starting indexes, 
!    if any, are also printed.
!
!  Modified:
!
!    16 June 2007
!
!  Author:
!
!    Robert Renka,
!    Department of Computer Science,
!    University of North Texas,
!    renka@cs.unt.edu
!
!  Reference:
!
!    Robert Renka,
!    Algorithm 751: TRIPACK, 
!    A Constrained Two-Dimensional Delaunay Triangulation Package,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 1, 1996.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NCC, the number of constraints.
!
!    Input, integer ( kind = 4 ) LCC(*), list of constraint curve starting indexes (or
!    dummy array of length 1 if NCC = 0).
!
!    Input, integer ( kind = 4 ) N, the number of nodes in the triangulation.
!    3 <= N <= 9999.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the coordinates of the nodes in the 
!    triangulation; not used unless PRNTX = TRUE.
!
!    Input, integer ( kind = 4 ) LIST(*), LPTR(*), LEND(N), data structure defining 
!    the triangulation.  Refer to subroutine TRMESH.
!
!    Input, logical PRNTX, TRUE if and only if X and Y are to be printed.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) k
  integer ( kind = 4 ) lcc(*)
  integer ( kind = 4 ) lend(n)
  integer ( kind = 4 ) list(*)
  integer ( kind = 4 ) lp
  integer ( kind = 4 ) lpl
  integer ( kind = 4 ) lptr(*)
  integer ( kind = 4 ) na
  integer ( kind = 4 ) nabor(100)
  integer ( kind = 4 ) nb
  integer ( kind = 4 ) ncc
  integer ( kind = 4 ) nd
  integer ( kind = 4 ) nl
  integer ( kind = 4 ), parameter :: nlmax = 60
  integer ( kind = 4 ), parameter :: nmax = 9999
  integer ( kind = 4 ) nn
  integer ( kind = 4 ) node
  integer ( kind = 4 ) nt
  logical prntx
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)

  nn = n
!
!  Print a heading and test the range of N.
!
  write ( *,100) nn
  if ( nn < 3 .or. nmax < nn ) then
!
!  N is outside its valid range.
!
    write (*,110)
    go to 5
  end if
!
!  Initialize NL (the number of lines printed on the current
!  page) and NB (the number of boundary nodes encountered).
!
  nl = 6
  nb = 0
  if (.not. prntx) then
!
!  Print LIST only.  K is the number of neighbors of NODE
!  which are stored in NABOR.
!
    write (*,101)

    do node = 1, nn

      lpl = lend(node)
      lp = lpl
      k = 0

      do

        k = k + 1
        lp = lptr(lp)
        nd = list(lp)
        nabor(k) = nd

        if ( lp == lpl ) then
          exit
        end if

      end do

      if ( nd <= 0 ) then
!
!  NODE is a boundary node.  Correct the sign of the last
!  neighbor, add 0 to the end of the list, and increment NB.
!
        nabor(k) = -nd
        k = k + 1
        nabor(k) = 0
        nb = nb + 1
      end if
!
!  Increment NL and print the list of neighbors.
!
      inc = (k-1) / 14 + 2
      nl = nl + inc

      if ( nlmax < nl ) then
        write (*,106)
        nl = inc
      end if

      write (*,103) node, nabor(1:k)
      if ( k /= 14 ) then
        write (*,105)
      end if

    end do

  else
!
!  Print X, Y, and LIST.
!
    write (*,102)

    do node = 1,nn

      lpl = lend(node)
      lp = lpl
      k = 0

      do

        k = k + 1
        lp = lptr(lp)
        nd = list(lp)
        nabor(k) = nd

        if ( lp == lpl ) then
          exit
        end if

      end do

      if ( nd <= 0 ) then
!
!  NODE is a boundary node.
!
        nabor(k) = -nd
        k = k + 1
        nabor(k) = 0
        nb = nb + 1
      end if
!
!  Increment NL and print X, Y, and NABOR.
!
      inc = (k-1) / 8 + 2
      nl = nl + inc

      if ( nlmax < nl ) then
        write (*,106)
        nl = inc
      end if

      write (*,104) node, x(node), y(node), &
                      (nabor(i), i = 1,k)
      if (k /= 8) write (*,105)

    end do

  end if
!
!  Print NB, NA, and NT (boundary nodes, arcs, and triangles).
!
  nt = 2 * nn - nb - 2
  na = nt + nn - 1

  if ( nlmax - 6 < nl ) then
    write (*,106)
  end if

  write (*,107) nb, na, nt
!
!  Print NCC and LCC.
!
5 continue

  write ( *, '(a)' ) ' '
  write ( *, 108 ) ncc
  write ( *, 109 ) lcc(1:ncc)

  return

  100 format (///,26x,'adjacency sets,    n = ',i5//)
  101 format (1x,'node',32x,'neighbors of node'//)
  102 format (1x,'node',5x,'x(node)',8x,'y(node)', &
          20x,'neighbors of node'//)
  103 format (1x,i4,5x,14i5/(1x,9x,14i5))
  104 format (1x,i4,2e15.6,5x,8i5/(1x,39x,8i5))
  105 format (1x)
  106 format (///)
  107 format (/1x,'nb = ',i4,' boundary nodes',5x, &
          'na = ',i5,' arcs',5x,'nt = ',i5, &
          ' triangles')
  108 format ('ncc =',i3,' constraint curves')
  109 format (1x,9x,14i5)
  110 format (1x,10x,'*** N is outside its valid range ***')
end
