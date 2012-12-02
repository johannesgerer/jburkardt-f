subroutine graph_arc_euler_circ ( nnode, nedge, inode, jnode, loop )

!*****************************************************************************80
!
!! GRAPH_ARC_EULER_CIRC finds an Euler circuit in an Eulerian graph.
!
!  Discussion:
!
!    An Euler circuit of a graph is a path that uses each edge exactly once.
!
!    A graph is Eulerian if it has an Euler circuit.
!
!    An Eulerian graph may have many circuits; this routine only finds one.
!
!  Modified:
!
!    12 October 2010
!
!  Author:
!
!    Hang Tong Lau
!
!  Reference:
!
!    Hang Tong Lau,
!    Combinatorial Heuristic Algorithms in FORTRAN,
!    Springer Verlag, 1986,
!    LC: QA402.5 L37.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes in the graph.
!
!    Input, integer ( kind = 4 ) NEDGE, the number of edges in the graph.
!
!    Input, integer ( kind = 4 ) INODE(NEDGE), JNODE(NEDGE), the two
!    end nodes of each edge.
!
!    Output, integer ( kind = 4 ) LOOP(NEDGE), the Euler circuit, as a
!    series of nodes.
!
  implicit none

  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) nnode

  logical copyon
  logical found
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibase
  integer ( kind = 4 ) iforwd
  integer ( kind = 4 ) inode(nedge)
  integer ( kind = 4 ) insert
  integer ( kind = 4 ) ipivot
  integer ( kind = 4 ) iwork1(nedge)
  integer ( kind = 4 ) iwork2(nedge)
  integer ( kind = 4 ) iwork3(nnode)
  integer ( kind = 4 ) iwork4(nnode)
  integer ( kind = 4 ) iwork5(nnode)
  integer ( kind = 4 ) iwork6(nnode)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jnode(nedge)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) locbas
  integer ( kind = 4 ) loop(nedge)
  integer ( kind = 4 ) nbreak
  integer ( kind = 4 ) ncopy
  integer ( kind = 4 ) numarc
  integer ( kind = 4 ) numnode
!
!  The number of times each node has been visited begins at 0.
!
  iwork3(1:nnode) = 0
  loop(1:nedge) = 0
  iwork1(1:nedge) = 0
  iwork2(1:nedge) = 0
!
!  Begin the Euler circuit with the edge INODE(1), JNODE(1).
!
  numarc = 1
  iwork2(1) = 1

  numnode = 1
  i = inode(1)
  iwork1(numnode) = i
  iwork3(i) = 1

  numnode = numnode + 1
  j = jnode(1)
  iwork1(numnode) = j
  iwork3(j) = 1

  ibase = j
  nbreak = 0
!
!  Look for the next arc.
!
30    continue

  do i = 2, nedge

    if ( iwork2(i) == 0 ) then

      if ( inode(i) == ibase ) then

        found = .true.
        ibase = jnode(i)

      else if ( jnode(i) == ibase ) then

        found = .true.
        ibase = inode(i)

      else

        found = .false.

      end if

      if ( found ) then
        iwork2(i) = 1
        numarc = numarc + 1
        numnode = numnode + 1
        if ( numnode <= nedge ) then
          iwork1(numnode) = ibase
        end if
        iwork3(ibase) = 1
        go to 30
      end if

    end if

  end do
!
!  A cycle has been found.
!
  if ( 0 < nbreak ) then
    numnode = numnode - 1
    iwork5(nbreak) = numnode
  end if

  if ( numarc < nedge )  then

    iwork1(numnode) = ibase
!
!  Find a node in the current Euler circuit.
!
    do i = 2, nedge

      if ( iwork2(i) == 0 ) then

        if ( iwork3(inode(i)) /= 0 ) then

          found = .true.
          j = inode(i)
          k = jnode(i)

        else if ( iwork3(jnode(i)) /= 0 ) then

          found = .true.
          j = jnode(i)
          k = inode(i)

        else

          found = .false.

        end if
!
!  Identify a path which will be added to the circuit.
!
        if ( found ) then
          nbreak = nbreak + 1
          iwork6(nbreak) = j
          ibase = k
          iwork3(k) = 1
          numnode = numnode + 1
          iwork4(nbreak) = numnode
          iwork1(numnode) = ibase
          iwork2(i) = 1
          numarc = numarc + 1
          go to 30
        end if

      end if

    end do

  end if
!
!  Form the Euler circuit.
!
  if ( nbreak == 0 ) then
    numnode = numnode - 1
    loop(1:numnode) = iwork1(1:numnode)
    return
  end if

  insert = 1
  ipivot = iwork6(insert)
  iforwd = 0

  do

    ncopy = 1
    ibase = iwork1(1)
    locbas = 1
    loop(ncopy) = ibase
!
!  A path identified before is added to the circuit.
!
80  continue

    if ( ibase == ipivot ) then

      j = iwork4(insert) + iforwd
      k = iwork5(insert) + iforwd

      do l = j, k
        ncopy = ncopy + 1
        loop(ncopy) = iwork1(l)
        iwork1(l) = 0
      end do

      ncopy = ncopy + 1
!
!  Add the intersecting node to the circuit.
!
      loop(ncopy) = ibase
      iforwd = iforwd + 1

      if ( ncopy < numnode ) then

        do

          if ( nedge <= ncopy ) then
            exit
          end if

          locbas = locbas + 1

          if ( nedge <= locbas ) then
            exit
          end if

          ibase = iwork1(locbas)

          if ( ibase /= 0 ) then
            ncopy = ncopy + 1
            loop(ncopy) = ibase
          end if

        end do

      end if

    else

      ncopy = ncopy + 1

       if ( ncopy <= numnode ) then
         locbas = locbas + 1
         ibase = iwork1(locbas)
         loop(ncopy) = ibase
         go to 80
       end if

     end if
!
!  Check if more paths are to be added to the circuit.
!
    copyon = .false.
    insert = insert + 1

    if ( insert <= nbreak ) then
      copyon = .true.
      ipivot = iwork6(insert)
    end if

    if ( .not. copyon ) then
      exit
    end if

    iwork1(1:nedge) = loop(1:nedge)

  end do

  return
end
subroutine graph_arc_min_path ( nnode, nedge, inode, jnode, cost, &
  start, last, num_path, ispath, xlen )

!*****************************************************************************80
!
!! GRAPH_ARC_MIN_PATH finds the shortest path between two nodes.
!
!  Discussion:
!
!    This routine is used by the approximate Steiner tree algorithm.
!
!  Modified:
!
!    11 September 1999
!
!  Author:
!
!    Hang Tong Lau
!
!  Reference:
!
!    Hang Tong Lau,
!    Combinatorial Heuristic Algorithms in FORTRAN,
!    Springer Verlag, 1986,
!    ISBN: 3540171614,
!    LC: QA402.5.L37.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes in the graph.
!
!    Input, integer ( kind = 4 ) NEDGE, the number of edges in the graph.
!
!    Input, integer ( kind = 4 ) INODE(NEDGE), JNODE(NEDGE), the edges
!    of the graph, describe by pairs of nodes.
!
!    Input, real ( kind = 8 ) COST(NEDGE), the distance or cost of each edge.
!
!    Input, integer ( kind = 4 ) START, LAST, are the two nodes between which a
!    shortest path is desired.
!
!    Output, integer ( kind = 4 ) NUM_PATH, the number of nodes in the shortest
!    path.  NUM_PATH is zero if no path could be found.
!
!    Output, integer ( kind = 4 ) ISPATH(NNODE), lists the nodes in the
!    shortest path, from ISPATH(1) to ISPATH(NUM_PATH).
!
!    Output, real ( kind = 8 ) XLEN, the length of the shortest path
!    from START to LAST.
!
  implicit none

  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) nnode

  real ( kind = 8 ) cost(nedge)
  real ( kind = 8 ) d
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) ient
  logical ifin
  integer ( kind = 4 ) ij
  integer ( kind = 4 ) inode(nedge)
  integer ( kind = 4 ) ispath(nnode)
  logical iwork1(nnode)
  integer ( kind = 4 ) iwork2(nnode)
  integer ( kind = 4 ) iwork3(nedge)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jnode(nedge)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) last
  integer ( kind = 4 ) num_path
  integer ( kind = 4 ) start
  real ( kind = 8 ) wk4(nnode)
  real ( kind = 8 ) xlen

  wk4(1:nnode) = huge ( xlen )
  iwork1(1:nnode) = .true.
  iwork2(1:nnode) = 0

  wk4(start) = 0.0D+00
  i = start
  iwork1(start) = .false.
  xlen = 0
!
!  For each forward arc originating at node I calculate
!  the length of the path to node I.
!
10    continue

  ic = 0

  do k = 1, nedge

    if ( inode(k) == i ) then
      ic = ic + 1
      iwork3(ic) = k
      ispath(ic) = jnode(k)
    end if

    if ( jnode(k) == i ) then
      ic = ic + 1
      iwork3(ic) = k
      ispath(ic) = inode(k)
    end if

  end do

  if ( 0 < ic ) then

    do l = 1, ic
      k = iwork3(l)
      j = ispath(l)
      if ( iwork1(j) ) then
        d = wk4(i) + cost(k)
        if ( d < wk4(j) ) then
          wk4(j) = d
          iwork2(j) = k
        end if
      end if
    end do

  end if
!
!  Find the minimum potential.
!
  d = huge ( d )
  ient = 0
  ifin = .false.

  do i = 1, nnode

    if ( iwork1(i) ) then
      ifin = .true.
      if ( wk4(i) < d ) then
        d = wk4(i)
        ient = i
      end if
    end if

  end do
!
!  Include the node in the current path.
!
  if ( d < huge ( d ) ) then
    iwork1(ient) = .false.
    if ( ient /= last ) then
      i = ient
      go to 10
    end if
  else
    if ( ifin ) then
      num_path = 0
      return
    end if
  end if

  ij = last
  num_path = 1
  ispath(1) = last

  do

    k = iwork2(ij)

    if ( inode(k) == ij ) then
      ij = jnode(k)
    else
      ij = inode(k)
    end if

    num_path = num_path + 1
    ispath(num_path) = ij

    if ( ij == start ) then
      exit
    end if

  end do

  l = num_path / 2
  j = num_path

  do i = 1, l
    call i4_swap ( ispath(i), ispath(j) )
    j = j - 1
  end do

  xlen = wk4(last)

  return
end
subroutine graph_arc_min_span_tree ( nnode, nedge, inode, jnode, cost, &
  itree, jtree, tree_cost )

!*****************************************************************************80
!
!! GRAPH_ARC_MIN_SPAN_TREE finds a minimum spanning tree of a graph.
!
!  Discussion:
!
!    This routine is used by the approximate Steiner tree algorithm.
!
!    The input graph is represented by a list of edges.
!
!  Modified:
!
!    21 July 2000
!
!  Author:
!
!    Hang Tong Lau
!
!  Reference:
!
!    Hang Tong Lau,
!    Combinatorial Heuristic Algorithms in FORTRAN,
!    Springer Verlag, 1986,
!    ISBN: 3540171614,
!    LC: QA402.5.L37.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes in the graph.
!
!    Input, integer ( kind = 4 ) NEDGE, the number of edges in the graph.
!
!    Input, integer ( kind = 4 ) INODE(NEDGE), JNODE(NEDGE), the start and end
!    nodes of the edges.
!
!    Input, real ( kind = 8 ) COST(NEDGE), the cost or length of each edge.
!
!    Output, integer ( kind = 4 ) ITREE(NNODE-1), JTREE(NNODE-1), the pairs of
!    nodes that form the edges of the spanning tree.
!
!    Output, real ( kind = 8 ) TREE_COST, the total cost or length of the
!    spanning tree.
!
  implicit none

  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) best
  real ( kind = 8 ) cost(nedge)
  real ( kind = 8 ) d
  logical free(nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) ij
  integer ( kind = 4 ) inode(nedge)
  integer ( kind = 4 ) itr
  integer ( kind = 4 ) itree(nnode-1)
  integer ( kind = 4 ) iwork1(nnode)
  integer ( kind = 4 ) iwork2(nnode)
  integer ( kind = 4 ) iwork4(nedge)
  integer ( kind = 4 ) iwork5(nedge)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jnode(nedge)
  integer ( kind = 4 ) jtree(nnode-1)
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kk
  integer ( kind = 4 ) l
  real ( kind = 8 ) tree_cost
  real ( kind = 8 ) wk6(nnode)

  wk6(1:nnode) = huge ( tree_cost )
  free(1:nnode) = .true.
  iwork1(1:nnode) = 0
  iwork2(1:nnode) = 0
  itree(1:nnode-1) = 0
  jtree(1:nnode-1) = 0
!
!  Find the first non-zero arc.
!
  do ij = 1, nedge
    if ( inode(ij) /= 0 ) then
      i = inode(ij)
      exit
    end if
  end do

  wk6(i) = 0.0D+00
  free(i) = .false.

  tree_cost = 0.0D+00

  do jj = 1, nnode - 1

    wk6(1:nnode) = huge ( tree_cost )

    do i = 1, nnode
!
!  For each forward arc originating at node I
!  calculate the length of the path to node I.
!
      if ( .not. free(i) ) then

        ic = 0

        do k = 1, nedge

          if ( inode(k) == i ) then
            ic = ic + 1
            iwork5(ic) = k
            iwork4(ic) = jnode(k)
          end if

          if ( jnode(k) == i ) then
            ic = ic + 1
            iwork5(ic) = k
            iwork4(ic) = inode(k)
          end if

        end do

        if ( 0 < ic ) then

          do l = 1, ic

            k = iwork5(l)
            j = iwork4(l)

            if ( free(j) ) then

              d = tree_cost + cost(k)

              if ( d < wk6(j) ) then
                wk6(j) = d
                iwork1(j) = i
                iwork2(j) = k
              end if

            end if

          end do

        end if

      end if

    end do
!
!  Identify the free node of minimum potential.
!
    d = huge ( d )
    best = 0

    do i = 1, nnode

      if ( free(i) ) then
        if ( wk6(i) < d ) then
          d = wk6(i)
          best = i
          itr = iwork1(i)
          kk = iwork2(i)
        end if
      end if

    end do
!
!  Add that node to the tree.
!
    if ( d < huge ( d ) ) then
      free(best) = .false.
      tree_cost = tree_cost + cost(kk)
      itree(jj) = itr
      jtree(jj) = best
    end if

  end do

  return
end
subroutine graph_arc_print ( nedge, inode, jnode, title )

!*****************************************************************************80
!
!! GRAPH_ARC_PRINT prints out a graph from an edge list.
!
!  Modified:
!
!    04 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEDGE, the number of edges.
!
!    Input, integer ( kind = 4 ) INODE(NEDGE), JNODE(NEDGE), the beginning
!    and end nodes of the edges.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) nedge

  integer ( kind = 4 ) i
  integer ( kind = 4 ) inode(nedge)
  integer ( kind = 4 ) jnode(nedge)
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = 1, nedge
    write ( *, '(i8,4x,i8,i8)' ) i, inode(i), jnode(i)
  end do

  return
end
subroutine graph_arc_weight_print ( nedge, inode, jnode, wnode, title )

!*****************************************************************************80
!
!! GRAPH_ARC_WEIGHT_PRINT prints out a weighted graph from an edge list.
!
!  Modified:
!
!    04 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEDGE, the number of edges.
!
!    Input, integer ( kind = 4 ) INODE(NEDGE), JNODE(NEDGE), the beginning
!    and end nodes of the edges.
!
!    Input, real ( kind = 8 ) WNODE(NEDGE), the weights of the edges.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) nedge

  integer ( kind = 4 ) i
  integer ( kind = 4 ) inode(nedge)
  integer ( kind = 4 ) jnode(nedge)
  character ( len = * ) title
  real ( kind = 8 ) wnode(nedge)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = 1, nedge
    write ( *, '(i8,4x,i8,i8,g14.6)' ) i, inode(i), jnode(i), wnode(i)
  end do

  return
end
subroutine graph_dist_min_span_tree3 ( nnode, dist, inode, jnode )

!*****************************************************************************80
!
!! GRAPH_DIST_MIN_SPAN_TREE3 finds a minimum spanning tree.
!
!  Discussion:
!
!    The input graph is represented by a distance matrix.
!
!  Modified:
!
!    03 January 2004
!
!  Author:
!
!    Hang Tong Lau
!
!  Reference:
!
!    Hang Tong Lau,
!    Combinatorial Heuristic Algorithms in FORTRAN,
!    Springer Verlag, 1986,
!    ISBN: 3540171614,
!    LC: QA402.5.L37.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, real ( kind = 8 ) DIST(NNODE,NNODE), the distance matrix.
!    DIST(I,J) is the distance from node I to node J.  The matrix
!    should be symmetric.  If there is no arc from node I to node J,
!    set DIST(I,J) = HUGE(1.0).
!
!    Output, integer ( kind = 4 ) INODE(NNODE), JNODE(NNODE); entries 1
!    through NNODE-1 describe the edges of the spanning tree as pairs of nodes.
!
  implicit none

  integer ( kind = 4 ) nnode

  real ( kind = 8 ) d
  real ( kind = 8 ) dist(nnode,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ient
  integer ( kind = 4 ) ij
  integer ( kind = 4 ) inode(nnode)
  integer ( kind = 4 ) itr
  integer ( kind = 4 ) iwork1(nnode)
  integer ( kind = 4 ) iwork2(nnode)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) jnode(nnode)
  integer ( kind = 4 ) kj
  integer ( kind = 4 ) kk4
  real ( kind = 8 ) tree_length
  real ( kind = 8 ) work(nnode)

  work(1:nnode) = huge ( d )
  iwork1(1:nnode) = 0
  iwork2(1:nnode) = 0
!
!  Find the first non-zero arc.
!
  i = 0

  do ij = 1, nnode
    do kj = 1, nnode
      if ( dist(ij,kj) < huge ( d ) ) then
        i = ij
        exit
      end if
    end do

    if ( i /= 0 ) then
      exit
    end if

  end do

  if ( i == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRAPH_DIST_MIN_SPAN_TREE - Fatal error!'
    write ( *, '(a)' ) '  There are no nonzero arcs in this graph.'
    stop
  end if

  work(i) = 0.0D+00
  iwork1(i) = 1
  tree_length = 0.0D+00
  kk4 = nnode - 1

  do jj = 1, kk4

    work(1:nnode) = huge ( d )

    do i = 1, nnode
!
!  For each forward arc originating at node I calculate
!  the length of the path to node I
!
      if ( iwork1(i) == 1 ) then
        do j = 1, nnode
          if ( dist(i,j) < huge ( 1.0D+00 ) .and. iwork1(j) == 0 ) then
            d = tree_length + dist(i,j)
            if ( d < work(j) ) then
              work(j) = d
              iwork2(j) = i
            end if
          end if
        end do
      end if

    end do
!
!  Find the minimum potential.
!
    d = huge ( d )
    ient = 0

    do i = 1, nnode
      if ( iwork1(i) == 0 .and. work(i) < d ) then
        d = work(i)
        ient = i
        itr = iwork2(i)
      end if
    end do
!
!  Include the node in the current path.
!
    if ( d < huge ( d ) ) then
      iwork1(ient) = 1
      tree_length = tree_length + dist(itr,ient)
      inode(jj) = itr
      jnode(jj) = ient
    end if

  end do

  return
end
subroutine graph_dist_print ( nnode, dist, title )

!*****************************************************************************80
!
!! GRAPH_DIST_PRINT prints a distance matrix.
!
!  Modified:
!
!    03 January 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, real ( kind = 8 ) DIST(NNODE,NNODE), the distance matrix.
!    DIST(I,J) is the distance from node I to node J.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) nnode

  real ( kind = 8 ) dist(nnode,nnode)
  character ( len = * ) title

  call r8mat_print ( nnode, nnode, dist, title )

  return
end
subroutine i4_swap ( i, j )

!*****************************************************************************80
!
!! I4_SWAP switches two I4's.
!
!  Modified:
!
!    30 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) I, J.  On output, the values of I and
!    J have been interchanged.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k

  k = i
  i = j
  j = k

  return
end
subroutine int_lp ( m, n, a, b, c, x, infs )

!*****************************************************************************80
!
!! INT_LP is a heuristic algorithm for the integer linear programming problem.
!
!  Discussion:
!
!    The integer ( kind = 4 ) linear programming problem has the form:
!
!      Maximize the objective function:
!
!        Sum ( 1 <= J <= N ) C(J) * X(J)
!
!      Subject to the constraint equations:
!
!        Sum ( 1 <= J <= N ) A(I,J) * X(J) <= B(I), for 1 <= I <= M,
!
!      and positivity and integer ( kind = 4 ) conditions:
!
!        0 <= X(J),    for 1 <= J <= N,
!        X(J) integer ( kind = 4 ), for 1 <= J <= N.
!
!    The heuristic algorithm described here will usually find a vector
!    X which is feasible, and gives a "reasonably good" objective
!    function value.
!
!  Modified:
!
!    11 September 1999
!
!  Author:
!
!    Hang Tong Lau
!
!  Reference:
!
!    Frederick Hillier,
!    Efficient heuristic procedure for integer linear programming
!    with an interior,
!    Operations Research,
!    Volume 17, pages 600-637, 1969.
!
!    Hang Tong Lau,
!    Combinatorial Heuristic Algorithms in FORTRAN,
!    Springer Verlag, 1986,
!    ISBN: 3540171614,
!    LC: QA402.5.L37.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of constraints.
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) A(M+1,N+M+1), contains the coefficients of the
!    constraints in the first M rows and N columns.  There is also an M+1
!    row used for working storage.  The last M+1 columns are also used as
!    working storage.
!
!    Input, real ( kind = 8 ) B(M), the right hand side of the constraints.
!
!    Input, real ( kind = 8 ) C(N), the coefficients of the objective function.
!
!    Output, real ( kind = 8 ) X(max ( M+1, N ) ), contains the solution.
!
!    Output, integer ( kind = 4 ) INFS, error indicator.
!    0, no error, a solution was found.
!    1, the linear program is infeasible.  No solution was found.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m+1,n+m+1)
  real ( kind = 8 ) b(m)
  real ( kind = 8 ) c(n)
  real ( kind = 8 ) eps
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) infs
  integer ( kind = 4 ) init
  integer ( kind = 4 ) isub
  integer ( kind = 4 ) iwk18(n)
  integer ( kind = 4 ) iwk19(m+1)
  integer ( kind = 4 ) iwk20(n+m)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) jval
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kval
  integer ( kind = 4 ) ll1
  integer ( kind = 4 ) ll2
  integer ( kind = 4 ) ll3
  integer ( kind = 4 ) more
  integer ( kind = 4 ) moveon
  integer ( kind = 4 ) n1
  real ( kind = 8 ) row_norm
  real ( kind = 8 ) temp
  real ( kind = 8 ) wk1(n,n)
  real ( kind = 8 ) wk2(m+1,m+1)
  real ( kind = 8 ) wk3(n)
  real ( kind = 8 ) wk4(n)
  real ( kind = 8 ) wk5(n)
  real ( kind = 8 ) wk6(n)
  real ( kind = 8 ) wk7(n)
  real ( kind = 8 ) wk8(n)
  real ( kind = 8 ) wk9(m)
  real ( kind = 8 ) wk10(m)
  real ( kind = 8 ) wk11(m)
  real ( kind = 8 ) wk12(m)
  real ( kind = 8 ) wk13(m+1)
  real ( kind = 8 ) wk14(m+1)
  real ( kind = 8 ) wk15(m+1)
  real ( kind = 8 ) wk16(n)
  logical              wk17(n)
  real ( kind = 8 ) x(n)

  eps = sqrt ( epsilon ( eps ) )
!
!  Normalize the coefficients.
!
  do i = 1, m

    row_norm = sqrt ( sum ( a(i,1:n)**2 ) )

    if ( row_norm /= 0.0D+00 ) then
      a(i,1:n) = a(i,1:n) / row_norm
      b(i) = b(i) / row_norm
    end if

  end do
!
!  Set up the data for the simplex subroutine.
!
  do i = m, 1, -1
    a(i+1,1:n) = a(i,1:n)
  end do

  a(1,1:n) = - c(1:n)
  a(1,n+1:n+m) = 0.0D+00

  do i = 2, m + 1
    a(i,n+1:n+m) = 0.0D+00
  end do

  do i = 2, m + 1
    j = n + i - 1
    a(i,j) = 1.0D+00
  end do
!
!  Set the right hand side.
!
  wk16(1) = 0.0D+00
  wk16(2:m+1) = b(1:m)
!
!  Use the simplex method to find an optimal noninteger solution.
!
  call simplex ( n, m, a, x, infs, wk2, wk13, wk14, wk16, iwk19, iwk20 )

  if ( infs /= 0 ) then
    return
  end if
!
!  Restore the matrix.
!
  do i = 1, m
    a(i,1:n) = a(i+1,1:n)
  end do

  wk15(1:m+1) = x(1:m+1)

  do j = 1, n
    isub = iwk20(j)
    if ( 0 < isub ) then
      wk7(j) = x(isub)
    else
      wk7(j) = 0.0D+00
    end if
  end do

  jj = 1
  ii = 0

  do j = 1, m + 1
    do i = 1, m + 1
      ii = ii + 1
      if ( m + 1 < ii ) then
        ii = 1
        jj = jj + 1
      end if
    a(i,n+j) = wk2(ii,jj)
    end do
  end do

  do i = 1, m
    do j = 1, m
      wk2(i,j) = a(i+1,n+j+1)
    end do
  end do
!
!  Select a second solution corresponding to the noninteger
!  optimal solution by tightening the binding constraints
!  sufficiently so that a rounded solution still satisfies
!  the original constraints.
!
  do i = 1, m

    isub = iwk20(n+i)

    if ( 0 < isub ) then
      if ( 0.0D+00 < x(isub) ) then
        wk12(i) = b(i)
        go to 220
      end if
    end if

    temp = 0.0D+00

    do j = 1, n
      if ( 0 < iwk20(j) ) then
        temp = temp + abs(a(i,j))
      end if
    end do

    wk12(i) = b(i) - 0.5D+00 * temp

220 continue

  end do
!
!  The second solution is stored in WK8.
!
  wk8(1:n) = 0.0D+00

  do i = 1, m
    isub = iwk19(i+1)
    if ( isub <= n ) then
      wk8(isub) = dot_product ( wk2(i,1:m), wk12(1:m) )
    end if
  end do
!
!  Perform the linear search along the line joined by two points
!  WK7 and WK8 found above.
!
  call parse1 ( n, m, a, b, x, wk4, wk6, wk7, wk8, wk9, wk16, wk17 )
!
!  Search for better feasible solutions.
!
  call parse2 ( n, m, a, b, c, x, wk3, wk5, wk10, iwk18, init )
!
!  Prepare for the interchange of two variables.
!
  if ( 1 < init ) then
    call parse4 ( n, c, wk1, iwk18, init )
  end if
!
!  Increase or decrease one variable by one as long as the
!  resulting better solution is feasible.
!
260   continue

  call parse3 ( n, m, a, c, x, wk10, kval )

  if ( init <= 1 ) then
    return
  end if

  wk6(1:n) = x(1:n)
  jval = 1
  kval = n
!
!  Check whether the changed value of a variable is negative.
!
280   continue

  call parse5 ( n, m, a, x, wk5, wk10, wk11, iwk18, moveon, jval )

  if ( moveon /= 1 ) then
     kval = jval + 1
     go to 320
  end if

290  continue

  isub = iwk18(kval)
  ll1 = 0
  ll2 = 0

  do i = 1, m

    if ( wk11(i) < 0.0D+00 ) then

      if ( a(i,isub) == 0.0D+00 ) then
        go to 310
      else if ( a(i,isub) < 0.0D+00 ) then
        ll1 = 1
      else
        ll2 = 2
      end if

    end if

  end do

  ll3 = ll1 + ll2
!
!  Investigate the possibility of changing X(JVAL) in one favorable
!  direction and seek for a change in another variable X(KVAL).
!
!  Try to decrease X(KVAL) to get a better solution.
!
  if ( ll3 == 2 ) then

    call parse6 ( n, m, a, c, x, wk1, wk3, wk5, wk10, wk11, &
      iwk18, more, jval, kval )

    if ( more == 1 ) then
      go to 280
    end if
!
!  Try to increase X(KVAL) to get a better solution.
!
  else if ( ll3 < 2 ) then

    call parse7 ( n, m, a, c, x, wk1, wk3, wk5, wk10, wk11, &
      iwk18, more, jval, kval )

    if ( more == 1 ) then
      go to 280
    end if

  end if

310   continue

  if ( jval /= kval-1 ) then
    kval = kval - 1
    go to 290
  end if
!
!  Consider changing X(JVAL) by one in the other direction and
!  make a small integer change of X(KVAL).
!
320   continue

  call parse8 ( n, m, a, x, wk1, wk5, wk10, wk11, iwk18, init, &
    jval, kval )

  n1 = min ( n - 1, init )
!
!  If one or more improved solution is found then repeat the search.
!
  if ( jval /= n1 ) then
    jval = jval + 1
    kval = n
    go to 280
  end if

  do i = 1, n
    if ( x(i) /= wk6(i) ) then
      go to 260
    end if
  end do

  return
end
subroutine k_center ( nnode, cost, kmax, knum, kset )

!*****************************************************************************80
!
!! K_CENTER is a heuristic algorithm for the K-center problem.
!
!  Discussion:
!
!    Let V be the set of nodes of a complete undirected graph of
!    NNODE nodes, with the edge from node I to J having nonnegative
!    cost COST(I,J), and COST(I,I) = 0.
!
!    Given an integer ( kind = 4 ) K between 1 and N, the K-center location
!    problem is to find a subset S of the nodes V, of size at most K, such
!    that
!
!      Z(S) = max ( i in V ) min ( j in S ) COST(I,J)
!
!    is minimized.
!
!    In other words, for any subset S of the nodes, we measure the
!    distance from every node in V to a node in S; the maximum of
!    these values is the "cost" of the subset S.  We are seeking the
!    set S with minimum cost.  This is a sort of discrete minimization
!    problem in the L-infinity norm.
!
!    In general, this is a hard problem.  However, in the particular case
!    where the COST function satisfies the triangle inequality, that is,
!
!      COST(I,K) <= COST(I,J) + COST(J,K) for all I, J, K,
!
!    then solutions close to the optimum can be found with this algorithm.
!    In fact, the computed solution is guaranteed to have a value of Z
!    no more than twice that of the optimal solution.
!
!  Modified:
!
!    11 September 1999
!
!  Author:
!
!    Hang Tong Lau
!
!  Reference:
!
!    Dorit Hochbaum, David Shmoys,
!    A best possible heuristic for the K-center problem,
!    Mathematics of Operations Research,
!    Volume 10, pages 180-184, 1985.
!
!    Hang Tong Lau,
!    Combinatorial Heuristic Algorithms in FORTRAN,
!    Springer Verlag, 1986,
!    ISBN: 3540171614,
!    LC: QA402.5.L37.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, real ( kind = 8 ) COST(NNODE,NNODE), the cost matrix, which contains
!    the "cost" for each pair of nodes, or the distance between them.
!
!    Input, integer ( kind = 4 ) KMAX, the maximum size of the subset of nodes
!    to be found.
!
!    Output, integer ( kind = 4 ) KNUM, the size of the subset of nodes found.
!    KNUM will be no greater than KMAX.
!
!    Output, integer ( kind = 4 ) KSET(KMAX), contains in entries 1 through KNUM
!    the list of nodes forming the subset.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) nnode

  real ( kind = 8 ) big
  real ( kind = 8 ) cost(nnode,nnode)
  real ( kind = 8 ) great
  integer ( kind = 4 ) high
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ifirst
  integer ( kind = 4 ) inode
  integer ( kind = 4 ) iwk1(nnode,nnode)
  integer ( kind = 4 ) iwk2(nnode)
  integer ( kind = 4 ) iwk3(nnode)
  integer ( kind = 4 ) iwk4((nnode*(nnode-1))/2)
  integer ( kind = 4 ) iwk5((nnode*(nnode-1))/2)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kmax
  integer ( kind = 4 ) knum
  integer ( kind = 4 ) kset(nnode)
  integer ( kind = 4 ) l
  integer ( kind = 4 ) low
  integer ( kind = 4 ) mid
  integer ( kind = 4 ) ncheck
  integer ( kind = 4 ) num1
  integer ( kind = 4 ) num2
  integer ( kind = 4 ) numk
  real ( kind = 8 ) small
  real ( kind = 8 ) work6((nnode*(nnode-1))/2)

  m = ( nnode * ( nnode - 1 ) ) / 2

  great = 1.0D+00

  do i = 1, nnode - 1
    do j = i + 1, nnode
      great = great + cost(i,j)
    end do
  end do

  small = great

  do i = 1, nnode

    big = - great

    do j = 1, nnode

      if ( j /= i ) then
        if ( big < cost(i,j) ) then
          big = cost(i,j)
          numk = i
        end if
      end if

    end do

    if ( big < small ) then
      small = big
      knum = numk
    end if

  end do

  ifirst = knum
!
!  Return the optimal solution if KMAX = 1.
!
  if ( kmax == 1 ) then
    knum = 1
    kset(1) = ifirst
    return
  end if
!
!  Sort the edges in order of increasing cost.
!
  k = 0
  do i = 1, nnode - 1
    do j = i + 1, nnode
      k = k + 1
      work6(k) = cost(i,j)
    end do
  end do
!
!  Sort the edges in order of increasing cost.
!
  call r8vec_sort_heap_index_a ( m, work6, iwk4 )

  do i = 1, m
    j = iwk4(i)
    iwk5(j) = i
  end do

  k = 0
  do i = 1, nnode - 1
    do j = i + 1, nnode
      k = k + 1
      iwk1(i,j) = iwk5(k)
      iwk1(j,i) = iwk5(k)
    end do
  end do
!
!  Binary search.
!
  low = 1
  high = m

100 continue

  if ( high /= low + 1 ) then

    mid = ( high + low ) / 2
!
!  Restrict to the subgraph with the original N nodes
!  but only having the first MID number of edges.
!
    numk = 0
    ncheck = nnode
    iwk2(1:nnode) = 0
    inode = ifirst

120 continue

    numk = numk + 1
!
!  Include node INODE into the solution set.
!
    iwk3(numk) = inode
!
!  Consider all nodes adjacent to node INODE.
!
    do k = 1, nnode

      if ( k /= inode ) then
        if ( iwk1(k,inode) <= mid ) then
!
!  Node K is adjacent to node INODE;
!  delete node K from the subgraph.
!
          if ( iwk2(k) == 0 ) then
            ncheck = ncheck - 1
            iwk2(k) = 2
          end if
!
!  Delete all nodes adjacent to node K from the subgraph.
!
          do l = 1, nnode

            if ( iwk2(l) == 0 ) then
              if ( l /= k ) then
                if ( iwk1(k,l) <= mid ) then
!
!  Node L is adjacent to node K;
!  delete node L from the subgraph.
!
                  ncheck = ncheck - 1
                  iwk2(l) = 2
                end if
              end if
            end if

          end do

        end if
      end if

    end do
!
!  Mark node INODE as being already selected.
!
    if ( iwk2(inode) == 0 ) then
      ncheck = ncheck - 1
    end if

    iwk2(inode) = 1
!
!  Continue the binary search if the subgraph is nonempty.
!
    if ( 0 < ncheck ) then
!
!  Pick the next center by the greedy heuristic.
!
      if ( ncheck <= 2 ) then

        do i = 1, nnode

          if ( iwk2(i) == 0 ) then
            inode = i
            go to 120
          end if

        end do

      end if

      small = great

      do i = 1, nnode

        if ( iwk2(i) == 0 ) then

          big = - great

          do j = 1, nnode

            if ( iwk2(j) == 0 ) then
              if ( j /= i ) then
                if ( big < cost(i,j) ) then
                  big = cost(i,j)
                  num1 = i
                end if
              end if
            end if

          end do

          if ( big < small ) then
            small = big
            num2 = num1
          end if

        end if

      end do

      inode = num2
      go to 120

    end if
!
!  Store up the temporary solution set.
!
    if ( numk <= kmax ) then

      high = mid
      knum = numk
      kset(1:knum) = iwk3(1:knum)

    else
      low = mid
    end if

    go to 100

  end if

  return
end
subroutine k_median ( m, n, c, k, isol )

!*****************************************************************************80
!
!! K_MEDIAN is a heuristic algorithm for the K-median problem.
!
!  Discussion:
!
!    If C is an M by N matrix, and K is an integer between 1 and M,
!    the K median location problem is to select K rows of C which
!    form a subset S that maximizes the value:
!
!      F(S) = sum ( 1 <= J <= N ) max ( I in S ) C(I,J)
!
!    The heuristic algorithm given here finds a near-optimum solution in
!    time proportional to the number of entries of C.
!
!  Modified:
!
!    11 September 1999
!
!  Author:
!
!    Hang Tong Lau
!
!  Reference:
!
!    Brian Kernighan, Shen Lin,
!    Heuristic solution of a signal design optimization problem,
!    Bell System Technical Journal,
!    Volume 52, pages 1145-1159, 1973.
!
!    Hang Tong Lau,
!    Combinatorial Heuristic Algorithms in FORTRAN,
!    Springer Verlag, 1986,
!    ISBN: 3540171614,
!    LC: QA402.5.L37.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns
!    of the cost matrix.
!
!    Input, real ( kind = 8 ) C(M,N), the M by N cost matrix.
!
!    Input, integer ( kind = 4 ) K, the number of rows to be used
!    in the maximizing set.
!
!    Output, integer ( kind = 4 ) ISOL(M), contains in entries 1 through
!    K the elements of the maximizing set.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) big
  real ( kind = 8 ) c(m,n)
  real ( kind = 8 ) eps
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) ij
  integer ( kind = 4 ) index
  integer ( kind = 4 ) ir
  integer ( kind = 4 ) isol(m)
  integer ( kind = 4 ) iw
  integer ( kind = 4 ) iwk5(n)
  integer ( kind = 4 ) iwk6(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  real ( kind = 8 ) optval
  real ( kind = 8 ) row_sum
  real ( kind = 8 ) temp
  real ( kind = 8 ) tmin
  real ( kind = 8 ) tot1
  real ( kind = 8 ) w
  real ( kind = 8 ) work1(n)
  real ( kind = 8 ) work2(n)
  real ( kind = 8 ) work3(m)
  real ( kind = 8 ) work4(m)
  real ( kind = 8 ) z

  eps = sqrt ( epsilon ( eps ) )
  big = huge ( big )

  work3(1:m) = 0.0D+00
!
!  When K = 1, obtaining the optimal solution is easy.
!
  if ( k == 1 ) then

    optval = -big

    do i = 1, m

      row_sum = sum ( c(i,1:n) )

      if ( optval < row_sum ) then
        optval = row_sum
        index = i
      end if

    end do

    isol(1) = index
    return

  end if
!
!  Find the largest and second-largest elements among the first
!  M elements of column J; store them in WORK1(J) and WORK2(J).
!
  do i = 1, n

    work1(i) = -big

    do j = 1, k
      if ( work1(i) < c(j,i) ) then
        work1(i) = c(j,i)
        iwk5(i) = j
      end if
    end do

    work2(i) = -big

    do j = 1, k
      if ( work2(i) < c(j,i) .and. j /= iwk5(i) ) then
        work2(i) = c(j,i)
        iwk6(i) = j
      end if
    end do

  end do

  do i = 1, k
    tot1 = 0.0D+00
    do j = 1, n
      if ( iwk5(j) == i ) then
        tot1 = tot1 + work1(j) - work2(j)
      end if
    end do
    work3(i) = tot1
  end do

  do i = 1, m
    isol(i) = i
  end do

  iflag = 0
  ir = k

  do

    if ( m < ir + 1 ) then
      ir = k
    end if

    ir = ir + 1

    do i = 1, k
      work4(isol(i)) = 0.0D+00
    end do

    work4(isol(ir)) = 0.0D+00
!
!  Find out which of the K+l rows will make the objective
!  value decrease the least if the row is removed.
!
     do j = 1, n
       z = c(isol(ir),j)
       if ( work2(j) < z .and. z <= work1(j) ) then
         work4(iwk5(j)) = work4(iwk5(j)) - z + work2(j)
       else
         if ( work1(j) < z ) then
           work4(isol(ir)) = work4(isol(ir)) + z - work1(j)
           work4(iwk5(j)) = work4(iwk5(j)) - work1(j) + work2(j)
         end if
      end if
    end do
!
!  Determine which row to remove.
!
    tmin = big
    l = 0

    do ij = 1, k

      i = isol(ij)
      temp = work3(i) + work4(i)

      if ( temp < tmin ) then
        tmin = temp
        l = ij
      end if

    end do

    i = isol(ir)
    temp = work3(i) + work4(i)

    if ( temp <= tmin .or. abs ( temp - tmin ) <= eps ) then
      iflag = iflag + 1
      if ( iflag == m - k ) then
        return
      end if
      cycle
    end if

    iflag = 0

    tot1 = 0.0D+00
    do i = 1, n
      if ( work1(i) < c(isol(ir),i) ) then
        tot1 = tot1 + c(isol(ir),i) - work1(i)
      end if
    end do

    work3(isol(ir)) = tot1
!
!  Update all arrays after exchanging two rows.
!
    do j = 1, n

      z = c(isol(ir),j)

      if ( z <= work2(j) ) then
!
!  Replacing the largest element by something no better than
!  the third largest.
!
        if ( isol(l) == iwk5(j) ) then

          work1(j) = work2(j)
          iwk5(j) = iwk6(j)
          w = - big

          do ij = 1, k
            ii = isol(ij)
            if ( w < c(ii,j) .and. ii /= iwk5(j) ) then
              w = c(ii,j)
              iw = ii
            end if
          end do

          ii = isol(ir)

          if ( w < c(ii,j) .and. ii /= iwk5(j) ) then
            w = c(ii,j)
            iw = ii
          end if

          work3(iwk6(j)) = work3(iwk6(j)) - w + work2(j)
          work2(j) = w
          iwk6(j) = iw

        end if
!
!  Replacing the second largest element by something no
!  better than the third largest.
!
        if ( isol(l) == iwk6(j) ) then

          w = -big

          do ij = 1, k
            ii = isol(ij)
            if ( w < c(ii,j) .and. ii /= iwk5(j) .and. ii /= iwk6(j) ) then
              w = c(ii,j)
              iw = ii
            end if
          end do

          ii = isol(ir)

          if ( w < c(ii,j) .and. ii /= iwk5(j) .and. ii /= iwk6(j) ) then
            w = c(ii,j)
            iw = ii
          end if

          work3(iwk5(j)) = work3(iwk5(j)) - w + work2(j)
          work2(j) = w
          iwk6(j) = iw

        end if

      else

        if ( work2(j) < z .and. z <= work1(j) ) then
!
!  Replacing the largest element by a new and smaller largest element.
!
          if ( isol(l) == iwk5(j) ) then

            work1(j) = z
            iwk5(j) = isol(ir)
            work4(isol(ir)) = work4(isol(ir)) + z - work2(j)
!
!  Z becomes the new second largest element.
!
          else

            work3(iwk5(j)) = work3(iwk5(j)) - z + work2(j)
            work2(j) = z
            iwk6(j) = isol(ir)

          end if

        else

          if ( work1(j) < z ) then
!
!  Replacing the largest element by a new largest element.
!
            if ( isol(l) == iwk5(j) ) then

              work4(isol(ir)) = work4(isol(ir)) + work1(j) - work2(j)
              iwk5(j) = isol(ir)
              work1(j) = z
!
!  Z becomes the largest element.
!
            else

              work3(iwk5(j)) = work3(iwk5(j)) - work1(j) + work2(j)
              work2(j) = work1(j)
              work1(j) = z
              iwk6(j) = iwk5(j)
              iwk5(j) = isol(ir)

            end if
          end if
        end if
      end if

    end do
!
!  Iterate until no improvement by exchange can be found.
!
    call i4_swap ( isol(l), isol(ir) )

    work3(isol(l)) = work4(isol(l))

    do i = k + 1, m
      work3(isol(i)) = 0.0D+00
    end do

  end do

  return
end
subroutine knapsack ( n, a, b, c, eps, m, objval, numsol, isol )

!*****************************************************************************80
!
!! KNAPSACK is a heuristic algorithm for the zero-one knapsack problem.
!
!  Discussion:
!
!    The zero-one knapsack problem seeks to maximize
!
!      Sum ( 1 <= J <= N ) C(J) * X(J)
!
!    subject to
!
!      Sum ( 1 <= J <= N ) A(J) * X(J) <= B
!
!    where
!
!      A(J) and C(J) are nonnegative for each J;
!      B is nonnegative;
!      X(J) = 0 or 1 for each J.
!
!    Intuitively, C(J) is the value of X(J), A(J) is the "size" or "cost"
!    of X(J), and B is an overall size or spending limit.  Thus,
!    the ratio C(J)/A(J) is a heuristic indication of the relative
!    value of X(J).
!
!  Modified:
!
!    03 January 2004
!
!  Author:
!
!    Hang Tong Lau
!
!  Reference:
!
!    Hang Tong Lau,
!    Combinatorial Heuristic Algorithms in FORTRAN,
!    Springer Verlag, 1986,
!    ISBN: 3540171614,
!    LC: QA402.5.L37.
!
!    Oscar Ibarra, Chul Kim,
!    Fast approximation algorithms for the knapsack and sum of
!    subset problems,
!    Journal of the Association for Computing Machinery,
!    Volume 22, pages 463-468, 1975.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) A(N), the coefficients of the constraints.
!
!    Input, real ( kind = 8 ) B, the right hand side of the constraints.
!
!    Input, real ( kind = 8 ) C(N), the coefficients of the objective function.
!
!    Input, real ( kind = 8 ) EPS, a positive real number prescribing the
!    required degree of accuracy in the solution.
!
!    Input, integer ( kind = 4 ) M, the smallest integer greater
!    than ( 3 / EPS )**2.
!
!    Output, real ( kind = 8 ) OBJVAL, the value of the objective function for
!    the suggested solution.
!
!    Output, integer ( kind = 4 ) NUMSOL, the number of nonzero variables in the
!    suggested solution.
!
!    Output, integer ( kind = 4 ) ISOL(N), contains the indices of the nonzero
!    variables in entries 1 through NUMSOL.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) b
  real ( kind = 8 ) big
  real ( kind = 8 ) c(n)
  real ( kind = 8 ) eps
  real ( kind = 8 ) est
  real ( kind = 8 ) gvalue
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icol1
  integer ( kind = 4 ) icol2
  integer ( kind = 4 ) igp1
  integer ( kind = 4 ) igp2
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) index
  integer ( kind = 4 ) irow
  integer ( kind = 4 ) isol(n)
  integer ( kind = 4 ) isub
  integer ( kind = 4 ) iwk1(n)
  integer ( kind = 4 ) iwk2(n)
  integer ( kind = 4 ) iwk3(n)
  integer ( kind = 4 ) iwk4(n)
  logical iwk5(m)
  integer ( kind = 4 ) iwk6(m)
  integer ( kind = 4 ) iwk7(n,m)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) numsol
  real ( kind = 8 ) objval
  real ( kind = 8 ) parm1
  real ( kind = 8 ) parm2
  real ( kind = 8 ) rhs
  logical switch
  real ( kind = 8 ) value
  real ( kind = 8 ) work1(2,m)
  real ( kind = 8 ) work2(2,m)
  real ( kind = 8 ) work3(n)
  real ( kind = 8 ) work4(n)
  real ( kind = 8 ) work5(n)
  real ( kind = 8 ) xtemp
  real ( kind = 8 ) ytemp
!
!  Reorder the variables.
!
  work3(1:n) = c(1:n) / a(1:n)
  big = 1.0D+00 + sum ( c )

  call r8vec_sort_heap_index_d ( n, work3, iwk1 )
!
!  Calculate the initial parameters.
!
  numsol = 0
  objval = 0.0D+00
  rhs = b
  j = 1

  do

    index = iwk1(j)

    if ( rhs < a(index) ) then
      exit
    end if

    rhs = rhs - a(index)
    j = j + 1

    if ( n < j ) then
      exit
    end if

  end do

  if ( n < j ) then

    numsol = n

    do i = 1, n
      objval = objval + c(j)
    end do

    do i = 1, n
      isol(i) = i
    end do

    return

  end if

  est = 0.0D+00
  do i = 1, j
    est = est + c(iwk1(i))
  end do

  xtemp = eps / 3.0D+00
  parm2 = xtemp * est
  parm1 = xtemp * parm2
!
!  Split the variables into two groups.
!
  igp1 = 0
  igp2 = 0

  do i = 1, n

    index = iwk1(i)

    if ( parm2 <= c(index) ) then
      igp1= igp1 + 1
      iwk4(igp1) = c(index) / parm1
      work3(igp1) = a(index)
      iwk2(igp1) = index
    else
      igp2 = igp2 + 1
      work5(igp2) = a(index)
      work4(igp2) = c(index)
      iwk3(igp2) = index
    end if

  end do

  if ( igp1 == 0 .or. m <= 1 ) then

    rhs = b
    j = 1

    do while ( j <= n )

      index = iwk1(j)

      if ( a(index) <= rhs ) then
        numsol = numsol + 1
        isol(numsol) = index
        objval = objval + c(index)
        rhs = rhs - a(index)
      end if

      j = j + 1

    end do

    return

  end if
!
!  Solve the sequence of problems in group one.
!
  work1(1,1) = 0.0D+00
  work1(2,1) = 0.0D+00
  work2(1,1) = 0.0D+00
  work2(2,1) = 0.0D+00

  work1(1,2:m) = big
  work2(1,2:m) = 0.0D+00

  i = iwk4(1) + 1
  work1(1,i) = work3(1)
  index = iwk2(1)
  work2(1,i) = c(index)
  switch = .true.
  icol2 = 1
  k = 2

  do while ( k <= igp1 )

    if ( switch ) then
      icol1 = 1
      icol2 = 2
    else
      icol1 = 2
      icol2 = 1
    end if

    switch = .not. switch

    do i = 2, m

      irow = i - 1 - iwk4(k)

      if ( irow < 0 ) then
        xtemp = big
      else
        irow = irow + 1
        xtemp = work1(icol1,irow)

        if ( xtemp < big ) then
          xtemp = xtemp + work3(k)
        end if

      end if

      ytemp = work1(icol1,i)

      if ( ytemp <= xtemp ) then
        work1(icol2,i) = ytemp
        work2(icol2,i) = work2(icol1,i)
      else
        work1(icol2,i) = xtemp
        index = iwk2(k)
        work2(icol2,i) = work2(icol1,irow) + c(index)
      end if

    end do

    k = k + 1

  end do
!
!  Find the maximum objective value.
!
  do i = 1, m

    value = work1(icol2,i)

    if ( value <= b ) then

      gvalue = 0.0D+00

      if ( 0 < igp2 ) then

        rhs = b - value

        if ( rhs < 0 ) then

          work2(icol2,i) = 0.0D+00

        else

          j = 1

          do while ( j <= igp2 )

            if ( work5(j) <= rhs ) then
              gvalue = gvalue + work4(j)
              rhs = rhs - work5(j)
            end if

            j = j + 1

          end do

        end if
      end if

      gvalue = gvalue + work2(icol2,i)

      if ( objval < gvalue ) then
        objval = gvalue
        isub = i
      end if

    end if

  end do
!
!  Obtain the best solution.
!
  if ( isub == 1 ) then

    if ( 0 < igp2 ) then

      rhs = b
      j = 1

      do while ( j <= igp2 )

        index = iwk3(j)

        if ( work5(j) <= rhs ) then
          numsol = numsol + i
          isol(numsol) = index
          rhs = rhs - work5(j)
        end if

        j = j + 1

      end do

    end if

    return

  end if

  iwk5(1:m) = .false.
  iwk6(1:m) = 0
  work1(1,1) = 0.0D+00
  work1(2,1) = 0.0D+00
  work2(1,1) = 0.0D+00
  work1(1,2:isub) = big

  i = iwk4(1) + 1
  work1(1,i) = work3(1)
  iwk6(i) = 1
  iwk7(1,i) = 1
  switch = .true.
  icol2 = 1
  k = 2

  do while ( k <= igp1 )

    if ( switch ) then
      icol1 = 1
      icol2 = 2
    else
      icol1 = 2
      icol2 = 1
    end if

    switch = .not. switch

    do i = 2, isub

      irow = i - 1 - iwk4(k)

      if ( irow < 0 ) then
        xtemp = big
      else
        irow = irow + 1
        xtemp = work1(icol1,irow)
        if ( xtemp < big ) then
          xtemp = xtemp + work3(k)
        end if
      end if

      ytemp = work1(icol1,i)

      if ( ytemp <= xtemp ) then
        work1(icol2,i) = ytemp
      else
        work1(icol2,i) = xtemp
      end if

      if ( work1(icol1,i) /= work1(icol2,i) ) then
        iwk5(i) = .true.
      end if

    end do

    do ii = 2, isub

      i = isub - ii + 2

      if ( iwk5(i) ) then

        iwk5(i) = .false.
        irow = i - iwk4(k)
        j = iwk6(irow)

        do while ( 0 < j )
          iwk7(j,i) = iwk7(j,irow)
          j = j - 1
        end do

        index = iwk6(irow) + 1
        iwk6(i) = index
        iwk7(index,i) = k

      end if

    end do

    k = k + 1

  end do

  j = iwk6(isub)

  do while ( 0 < j )
    numsol = numsol + 1
    index = iwk7(j,isub)
    isol(numsol) = iwk2(index)
    j = j - 1
  end do

  if ( 0 < igp2 ) then

    rhs = b - work1(icol2,isub)
    j = 1

    do while ( j <= igp2 )

      if ( work5(j) <= rhs ) then
        numsol = numsol + 1
        isol(numsol) = iwk3(j)
        rhs = rhs - work5(j)
      end if

      j = j + 1

    end do

  end if

  return
end
subroutine multi_knap ( m, n, a, b, c, objval, numsol, isol )

!*****************************************************************************80
!
!! MULTI_KNAP: heuristic for the multidimensional 0/1 knapsack problem.
!
!  Discussion:
!
!    Maximize:
!
!      Sum ( 1 <= J <= N ) C(J) * X(J)
!
!    Subject to constraints I = 1 to M:
!
!      Sum ( 1 <= J <= N ) A(I,J) * X(J) <= B(I)
!
!    With
!
!      X(J) = 0 or 1.
!
!    Note that this is NOT the multiple knapsack problem, but rather a
!    form of 0-1 linear programming.
!
!  Modified:
!
!    11 September 1999
!
!  Author:
!
!    Hang Tong Lau
!
!  Reference:
!
!    Hang Tong Lau,
!    Combinatorial Heuristic Algorithms in FORTRAN,
!    Springer Verlag, 1986,
!    ISBN: 3540171614,
!    LC: QA402.5.L37.
!
!    Yoshiaki Toyoda,
!    A simplified algorithm for obtaining approximate solutions
!    to zero-one programming problems,
!    Management Science,
!    Volume 21, pages 1417-1427, 1975.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of constraints.
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix of coefficients of the
!    constraints.
!
!    Input, real ( kind = 8 ) B(M), the right hand sides of the constraints.
!
!    Input, real ( kind = 8 ) C(N), the coefficients of the objective function.
!
!    Output, real ( kind = 8 ) OBJVAL, the value of the objective function
!    for the suggested solution.
!
!    Output, integer ( kind = 4 ) NUMSOL, the number of nonzero variables in the
!    suggested solution.
!
!    Output, integer ( kind = 4 ) ISOL(N), contains the indices of the
!    nonzero variables in entries 1 through NUMSOL.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) b(m)
  real ( kind = 8 ) big
  real ( kind = 8 ) c(n)
  logical check
  real ( kind = 8 ) cmax
  real ( kind = 8 ) col_sum
  real ( kind = 8 ) eps
  real ( kind = 8 ) grad
  real ( kind = 8 ) grdmax
  integer ( kind = 4 ) i
  integer ( kind = 4 ) isol(n)
  integer ( kind = 4 ) iwk3(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jgdmax
  integer ( kind = 4 ) jgm
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) ncan
  integer ( kind = 4 ) numcan
  integer ( kind = 4 ) numsol
  real ( kind = 8 ) objval
  real ( kind = 8 ) orimov
  real ( kind = 8 ) rootm
  real ( kind = 8 ) small
  real ( kind = 8 ) sum2
  real ( kind = 8 ) sumsq
  real ( kind = 8 ) wk1(m)
  real ( kind = 8 ) wk2(m)

  eps = sqrt ( epsilon ( eps ) )
!
!  Compute the machine infinity.
!
  big = huge ( big )
  small = tiny ( small )
!
!  Initialize.
!
  do i = 1, m
    a(i,1:n) = a(i,1:n) / b(i)
  end do

  wk1(1:m) = 0.0D+00
  objval = 0.0D+00
  numsol = 0

  do j = 1, n
    iwk3(j) = j
  end do

  rootm  = sqrt ( real ( m, kind = 8 ) )
  numcan = n

60 continue

  ncan = numcan
  numcan = 0
  grdmax = -big
!
!  When the resource requirement vector is zero.
!
  do j = 1, ncan
!
!  Select variables.
!
    jj = iwk3(j)
    do i = 1, m
      if ( eps < wk1(i) + a(i,jj) - 1.0D+00 ) then
        go to 90
      end if
    end do
    numcan = numcan + 1
    iwk3(numcan) = jj
!
!  Compute the effective gradients.
!
    col_sum = sum ( a(1:m,jj) )

    if ( col_sum <= small ) then
      grad = big
    else
      grad = c(jj) * rootm / col_sum
    end if
!
!  Find the variable whose effective gradient is the largest.
!
    if ( grdmax < grad ) then
      grdmax = grad
      jgdmax = jj
      jgm = numcan
    end if

90  continue

  end do
!
!  Accept the variable whose effective gradient is the largest.
!
  if ( numcan <= 0 ) then
    return
  end if

  if ( numcan == 1 ) then

    objval = objval + c(jgdmax)

    wk1(1:m) = wk1(1:m) + a(1:m,jgdmax)

    numsol = numsol + 1
    isol(numsol) = jgdmax
    return

  end if

  objval = objval + c(jgdmax)
  numsol = numsol + 1
  isol(numsol) = jgdmax
  iwk3(jgm) = iwk3(numcan)
  numcan = numcan - 1

  check = .true.

  wk1(1:m) = wk1(1:m) + a(1:m,jgdmax)

  do i = 1, m
    if ( eps < wk1(i) ) then
      check = .false.
    end if
  end do

  if ( check ) then
    go to 60
  end if

  do

    cmax = maxval ( wk1(1:m) )

    orimov = cmax * cmax

    wk2(1:m) = wk1(1:m) - orimov

    do i = 1, m
      if ( wk2(i) <= eps ) then
        wk2(i) = 0.0D+00
      end if
    end do

    sumsq = sqrt ( sum ( wk2(1:m)**2 ) )

    ncan = numcan
    numcan = 0
    grdmax = - big
!
!  When the resource requirement vector is non-zero.
!
    do j = 1, ncan
!
!  Select variables.
!
      jj = iwk3(j)
      do i = 1, m
        if ( eps < wk1(i) + a(i,jj) - 1.0D+00 ) then
          go to 170
        end if
      end do

      numcan = numcan + 1
      iwk3(numcan) = jj
!
!  Compute the effective gradients.
!
      sum2 = dot_product ( wk2(1:m), a(1:m,jj) )

      if ( sum2 <= small ) then
        grad = big
      else
        grad = c(jj) * sumsq / sum2
      end if
!
!  Find the variable whose effective gradient is the largest.
!
      if ( grdmax < grad ) then
        grdmax = grad
        jgdmax = jj
        jgm = numcan
      end if

170   continue

    end do
!
!  Accept the variable whose effective gradient is the largest.
!
    if ( numcan <= 0 ) then
      exit
    end if

    objval = objval + c(jgdmax)
    numsol = numsol + 1
    isol(numsol) = jgdmax
    wk1(1:m) = wk1(1:m) + a(1:m,jgdmax)

    if ( numcan == 1 ) then
      exit
    end if

    iwk3(jgm) = iwk3(numcan)
    numcan = numcan - 1

  end do

  return
end
subroutine parse1 ( n, m, a, b, x, wk4, wk6, wk7, wk8, wk9, wk16, wk17 )

!*****************************************************************************80
!
!! PARSE1 is used by INT_LP.
!
!  Modified:
!
!    11 September 1999
!
!  Author:
!
!    Hang Tong Lau
!
!  Reference:
!
!    Hang Tong Lau,
!    Combinatorial Heuristic Algorithms in FORTRAN,
!    Springer Verlag, 1986,
!    ISBN: 3540171614,
!    LC: QA402.5.L37.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m+1,n+m+1)
  real ( kind = 8 ) b(m)
  real ( kind = 8 ) delps
  real ( kind = 8 ) delta
  real ( kind = 8 ) eps
  integer ( kind = 4 ) i
  integer ( kind = 4 ) incr
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhigh
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kk
  integer ( kind = 4 ) ll
  real ( kind = 8 ) temp
  logical wk17(n)
  real ( kind = 8 ) wk4(n)
  real ( kind = 8 ) wk6(n)
  real ( kind = 8 ) wk7(n)
  real ( kind = 8 ) wk8(n)
  real ( kind = 8 ) wk9(m)
  real ( kind = 8 ) wk16(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xwork
  real ( kind = 8 ) xx1
  real ( kind = 8 ) xx2
  real ( kind = 8 ) yy1
  real ( kind = 8 ) yy2
!
!  Linear search.
!
  eps = - huge ( eps )
  delps = 0.05D+00
  delta = 0.0D+00

  do

    temp = 1.0D+00 - delta
!
!  Form a linear combination of WK7 and WK8.
!
    wk6(1:n) = wk7(1:n) + delta * ( wk8(1:n) - wk7(1:n) )

    do j = 1, n
      if ( wk6(j) <= 0.0D+00 ) then
        x(j) = 0.0D+00
      else
        x(j) = aint ( wk6(j) + 0.5D+00 )
      end if
    end do
!
!  Compute WK9, the degree of infeasibility.
!
    do i = 1, m
      wk9(i) = dot_product ( a(i,1:n), x(1:n) ) - b(i)
    end do

    xwork = 0.0D+00
    do i = 1, m
      if ( 0.0D+00 < wk9(i) ) then
        xwork = xwork + wk9(i)
      end if
    end do

70  continue

    if ( xwork <= 0.0D+00 ) then
      return
    end if

    ll = m
!
!  Compute the change in solution X that improves the objective value.
!  This is indicated by WK4.
!
    incr = 0
    xx1 = eps
    kk = 0

    do j = 1, n

      yy1 = 0.0D+00
      yy2 = 0.0D+00

      do i = 1, ll
        if ( 0.0D+00 < wk9(i) ) then
          yy1 = yy1 + a(i,j)
        end if
      end do

      xx2 = abs ( yy1 )

      if ( yy1 == 0.0D+00 ) then

        wk4(j) = 0.0D+00
        wk17(j) = .false.

      else

        if ( yy1 < 0.0D+00 ) then

          wk4(j) = 1.0D+00
          do i = 1, ll
            yy1 = wk9(i) + a(i,j)
            if ( 0.0D+00 < yy1 ) then
              yy2 = yy2 + yy1
            end if
          end do

        else if ( 0.0D+00 < yy1 ) then

          if ( 0.0D+00 < x(j) ) then
            wk4(j) = -1.0D+00
            do i = 1, ll
              yy1 = wk9(i) - a(i,j)
              if ( 0.0D+00 < yy1 ) then
                yy2 = yy2 + yy1
              end if
            end do
          end if

        end if

        wk16(j) = yy2

        if ( wk16(j) <= 0.0D+00 ) then

          if ( xx1 < xx2 ) then
            xx1 = xx2
            kk = j
          end if

        else

          if ( yy2 < xwork ) then
            wk17(j) = .true.
            k = j
            incr = incr + 1
          else
            wk17(j) = .false.
          end if

        end if

      end if

    end do

    if ( kk /= 0 ) then
      x(kk) = x(kk) + wk4(kk)
      return
    end if
!
!  Compute the improvement of the objective function.
!
    if ( 0 < incr ) then

      if ( 1 < incr ) then

        temp = eps

        do j = 1, n

          if ( wk17(j) ) then

            if ( temp < xwork - wk16(j) ) then
              jhigh = j
              temp = xwork - wk16(j)
            end if

          end if

        end do

        k = jhigh

      end if

      x(k) = x(k) + wk4(k)
      wk9(1:ll) = wk9(1:ll) + a(1:ll,k) * wk4(k)
      xwork = wk16(k)
      go to 70

    end if
!
!  Continue with the linear search.
!
    delta = delta + delps

    if ( 1.0D+00 < delta ) then
      exit
    end if

  end do

  return
end
subroutine parse2 ( n, m, a, b, c, x, wk3, wk5, wk10, iwk18, init )

!*****************************************************************************80
!
!! PARSE2 is used by INT_LP.
!
!  Modified:
!
!    11 September 1999
!
!  Author:
!
!    Hang Tong Lau
!
!  Reference:
!
!    Hang Tong Lau,
!    Combinatorial Heuristic Algorithms in FORTRAN,
!    Springer Verlag, 1986,
!    ISBN: 3540171614,
!    LC: QA402.5.L37.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m+1,n+m+1)
  real ( kind = 8 ) b(m)
  real ( kind = 8 ) c(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) init
  integer ( kind = 4 ) iwk18(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ll
  real ( kind = 8 ) wk3(n)
  real ( kind = 8 ) wk5(n)
  real ( kind = 8 ) wk10(m)
  real ( kind = 8 ) temp
  real ( kind = 8 ) x(n)
!
!  Improve the feasible solution by changing some variable by one.
!  This is the initialization part.
!
  do j = 1, n
    iwk18(j) = j
  end do

  wk3(1:n) = c(1:n)

  do i = 1, n - 1

    temp = abs ( wk3(i) )

    do j = i + 1, n

      if ( temp < abs ( wk3(j) ) ) then

        temp = abs ( wk3(j) )

        call i4_swap ( iwk18(i), iwk18(j) )
        call r8_swap ( wk3(i), wk3(j) )

      end if

    end do

  end do

  do j = 1, n
    if ( 0.0D+00 < abs ( wk3(j) ) ) then
      init = j
    end if
  end do
!
!  Set WK5 to indicate whether the objective coefficient is positive.
!
  do j = 1, init

    ll = iwk18(j)

    if ( wk3(j) < 0.0D+00 ) then
      wk5(ll) = - 1.0D+00
    else if ( 0.0D+00 < wk3(j) ) then
      wk5(ll) = 1.0D+00
    end if

  end do
!
!  Compute WK10, the feasibility test slacks.
!
  do i = 1, m
    wk10(i) = b(i) - dot_product ( a(i,1:n), x(1:n) )
  end do

  return
end
subroutine parse3 ( n, m, a, c, x, wk10, kval )

!*****************************************************************************80
!
!! PARSE3 is used by INT_LP.
!
!  Modified:
!
!    21 July 2000
!
!  Author:
!
!    Hang Tong Lau
!
!  Reference:
!
!    Hang Tong Lau,
!    Combinatorial Heuristic Algorithms in FORTRAN,
!    Springer Verlag, 1986,
!    ISBN: 3540171614,
!    LC: QA402.5.L37.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m+1,n+m+1)
  real ( kind = 8 ) c(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) kval
  real ( kind = 8 ) temp
  real ( kind = 8 ) wk10(m)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xx1
  real ( kind = 8 ) xx2
!
!  Iteratively change the value of some variable by 1 to improve
!  the objective value.
!
  do

    kval = 1
    xx1 = 0.0D+00
!
!  The variable that gives a large increase in the objective value
!  is chosen to be changed.
!
    do j = 1, n

      if ( c(j) /= 0.0D+00 ) then

        if ( 0.0D+00 < c(j) .or. 0.0D+00 < x(j) ) then

          temp = huge ( 1.0D+00 )

          do i = 1, m

            if ( 0.0D+00 < c(j) * a(i,j) ) then

              xx2 = wk10(i) / abs ( a(i,j) )

              if ( xx2 < 0.0D+00 )  then
                if ( xx2 /= aint ( xx2 ) ) then
                  xx2 = aint ( xx2 ) - 1.0D+00
                end if
              else if ( 0.0D+00 < xx2 ) then
                xx2 = aint ( xx2 )
              end if

              temp = min ( temp, xx2 )

            end if

          end do

          if ( xx1 <= abs ( c(j) ) * temp ) then
            xx1 = abs ( c(j) ) * temp
            kval = j
          end if

        end if

      end if

    end do
!
!  Increase or decrease X(KVAL) by one.
!
    if ( xx1 == 0.0D+00 ) then
      exit
    end if

    if ( c(kval) <= 0.0D+00 ) then
      x(kval) = x(kval) - 1.0D+00
      wk10(1:m) = wk10(1:m) + a(1:m,kval)
    else
      x(kval) = x(kval) + 1.0D+00
      wk10(1:m) = wk10(1:m) - a(1:m,kval)
    end if

  end do

  return
end
subroutine parse4 ( n, c, wk1, iwk18, init )

!*****************************************************************************80
!
!! PARSE4 is used by INT_LP.
!
!  Modified:
!
!    11 September 1999
!
!  Author:
!
!    Hang Tong Lau
!
!  Reference:
!
!    Hang Tong Lau,
!    Combinatorial Heuristic Algorithms in FORTRAN,
!    Springer Verlag, 1986,
!    ISBN: 3540171614,
!    LC: QA402.5.L37.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c(n)
  integer ( kind = 4 ) init
  integer ( kind = 4 ) iwk18(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j3
  integer ( kind = 4 ) k
  real ( kind = 8 ) temp
  real ( kind = 8 ) wk1(n,n)
!
!  Prepare for the interchange of two variables.
!
!  Determine the bounds on the size of changes.
!
!  The ratios of the objective coefficients are computed and stored in WK1.
!
  do j = 1, init - 1

    j2 = iwk18(j)

    do k = j + 1, init

      j3 = iwk18(k)
      temp = abs ( c(j2) / c(j3) )

      if ( temp <= 1.0D+00 ) then
        wk1(j,k) = aint ( temp - 1.0D+00 )
      else
        if ( aint ( temp ) == temp ) then
          wk1(j,k) = aint ( temp - 1.0D+00 )
        else
          wk1(j,k) = aint ( temp - 1.0D+00 ) + 1.0D+00
        end if
      end if

      if ( -1.0D+00 <= temp ) then
        wk1(k,j) = aint ( temp + 1.0D+00 )
      else
        if ( aint ( temp ) == temp ) then
          wk1(k,j) = aint ( temp + 1.0D+00 )
        else
          wk1(k,j) = aint ( temp + 1.0D+00 ) - 1.0D+00
        end if
      end if

    end do

  end do

  return
end
subroutine parse5 ( n, m, a, x, wk5, wk10, wk11, iwk18, moveon, jval )

!*****************************************************************************80
!
!! PARSE5 is used by INT_LP.
!
!  Modified:
!
!    11 September 1999
!
!  Author:
!
!    Hang Tong Lau
!
!  Reference:
!
!    Hang Tong Lau,
!    Combinatorial Heuristic Algorithms in FORTRAN,
!    Springer Verlag, 1986,
!    ISBN: 3540171614,
!    LC: QA402.5.L37.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m+1,n+m+1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icon
  integer ( kind = 4 ) isub
  integer ( kind = 4 ) iwk18(n)
  integer ( kind = 4 ) jval
  integer ( kind = 4 ) moveon
  real ( kind = 8 ) wk5(n)
  real ( kind = 8 ) wk10(m)
  real ( kind = 8 ) wk11(m)
  real ( kind = 8 ) x(n)
!
!  Check if the changed value of X(JVAL) is negative.
!
  isub = iwk18(jval)

  do

    if ( x(isub) < - wk5(isub) ) then
      moveon = 0
      exit
    end if
!
!  If the change is negative, then check whether this change is
!  feasible without changing another variable.
!
    icon = 0

    do i = 1, m
      wk11(i) = wk10(i) - wk5(isub) * a(i,isub)
      if ( wk11(i) < 0.0D+00 ) then
        icon = 1
      end if
    end do

    if ( icon == 1 ) then
      moveon = 1
      exit
    end if

    x(isub) = x(isub) + wk5(isub)
    wk10(1:m) = wk11(1:m)

  end do

  return
end
subroutine parse6 ( n, m, a, c, x, wk1, wk3, wk5, wk10, wk11, &
  iwk18, more, jval, kval )

!*****************************************************************************80
!
!! PARSE6 is used by INT_LP.
!
!  Modified:
!
!    11 September 1999
!
!  Author:
!
!    Hang Tong Lau
!
!  Reference:
!
!    Hang Tong Lau,
!    Combinatorial Heuristic Algorithms in FORTRAN,
!    Springer Verlag, 1986,
!    ISBN: 3540171614,
!    LC: QA402.5.L37.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m+1,n+m+1)
  real ( kind = 8 ) c(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) isub
  integer ( kind = 4 ) iwk18(n)
  integer ( kind = 4 ) jval
  integer ( kind = 4 ) jsub
  integer ( kind = 4 ) kval
  integer ( kind = 4 ) more
  real ( kind = 8 ) temp
  real ( kind = 8 ) wk1(n,n)
  real ( kind = 8 ) wk3(n)
  real ( kind = 8 ) wk5(n)
  real ( kind = 8 ) wk10(m)
  real ( kind = 8 ) wk11(m)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xwk1
  real ( kind = 8 ) xwk2
  real ( kind = 8 ) xwk3
!
!  Check for the change of variables X(JVAL) and X(KVAL) and identify a
!  better solution.
!
  isub = iwk18(kval)

  if ( 0.0D+00 < wk3(kval) ) then
    xwk2 = - min ( x(isub), wk1(jval,kval) )
  else
    xwk2 = - x(isub)
  end if

  xwk3 = huge ( 1.0D+00 )
  isub = iwk18(kval)

  do i = 1, m

    if ( wk11(i) < 0.0D+00 ) then

      temp = wk11(i) / a(i,isub)

      if ( temp /= aint ( temp ) ) then
        temp = aint ( temp ) - 1.0D+00
      end if

      xwk3 = min ( xwk3, temp )

    end if

  end do
!
!  Check whether feasibility can be restored.
!
  if ( xwk3 < xwk2 ) then
    more = 0
    return
  end if

  isub = iwk18(kval)
  xwk1 = - huge ( 1.0D+00 )

  do i = 1, m

    if ( a(i,isub) < 0.0D+00 ) then

      temp = wk11(i) / a(i,isub)

      if ( temp < 0.0D+00 ) then
        temp = aint ( temp )
      else
        if ( 0.0D+00 < temp ) then
          if ( temp /= aint ( temp ) ) then
            temp = aint ( temp ) + 1.0D+00
          end if
        end if

      end if

      xwk1 = max ( xwk1, temp )

    end if

  end do
!
!  Check if an improved solution can be obtained.
!
  xwk2 = max ( xwk2, xwk1 )

  if ( xwk3 < xwk2 ) then
    more = 0
    return
  end if

  isub = iwk18(kval)

  if ( c(isub) <= 0.0D+00 ) then
    jsub = iwk18(jval)
    x(jsub) = x(jsub) + wk5(jsub)
    x(isub) = x(isub) + xwk2
    wk10(1:m) = wk11(1:m) - xwk2 * a(1:m,isub)
  else
    x(isub) = x(isub) + xwk3
    wk10(1:m) = wk11(1:m) - xwk3 * a(1:m,isub)
    isub = iwk18(jval)
    x(isub) = x(isub) + wk5(isub)
  end if

  more = 1

  return
end
subroutine parse7 ( n, m, a, c, x, wk1, wk3, wk5, wk10, wk11, &
  iwk18, more, jval, kval )

!*****************************************************************************80
!
!! PARSE7 is used by INT_LP.
!
!  Modified:
!
!    11 September 1999
!
!  Author:
!
!    Hang Tong Lau
!
!  Reference:
!
!    Hang Tong Lau,
!    Combinatorial Heuristic Algorithms in FORTRAN,
!    Springer Verlag, 1986,
!    ISBN: 3540171614,
!    LC: QA402.5.L37.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m+1,n+m+1)
  real ( kind = 8 ) c(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) isub
  integer ( kind = 4 ) iwk18(n)
  integer ( kind = 4 ) jval
  integer ( kind = 4 ) jsub
  integer ( kind = 4 ) kval
  integer ( kind = 4 ) more
  real ( kind = 8 ) temp
  real ( kind = 8 ) wk1(n,n)
  real ( kind = 8 ) wk3(n)
  real ( kind = 8 ) wk5(n)
  real ( kind = 8 ) wk10(m)
  real ( kind = 8 ) wk11(m)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xwk1
  real ( kind = 8 ) xwk2
  real ( kind = 8 ) xwork
!
!  Try to increase X(KVAL) to get a better solution.
!
  if ( wk3(kval) < 0.0D+00 ) then
    xwk1 = wk1(jval,kval)
  else
    xwk1 = huge ( 1.0D+00 )
  end if

  isub = iwk18(kval)
  temp = - huge ( 1.0D+00 )

  do i = 1, m

    if ( wk11(i) < 0.0D+00 ) then
      xwork = wk11(i) / a(i,isub)
      if ( aint ( xwork ) /= xwork ) then
        xwork = aint ( xwork ) + 1.0D+00
      end if
      temp = max ( temp, xwork )
    end if

  end do
!
!  Check whether feasibility can be restored.
!
  xwk2 = temp

  if ( xwk1 < xwk2 ) then
    more = 0
    return
  end if

  isub = iwk18(kval)
  temp = huge ( 1.0D+00 )

  do i = 1, m

    if ( 0.0D+00 < a(i,isub) ) then

      xwork = wk11(i) / a(i,isub)

      if ( 0.0D+00 < xwork ) then
        xwork = aint ( xwork )
      else if ( xwork < 0.0D+00 ) then
        if ( xwork /= aint ( xwork ) ) then
          xwork = aint ( xwork ) - 1.0D+00
        end if
      end if

      temp = min ( temp, xwork )

    end if

  end do
!
!  Check if an improved solution can be obtained.
!
  xwk1 = min ( xwk1, temp )

  if ( xwk1 < xwk2 ) then
    more = 0
    return
  end if

  isub = iwk18(kval)

  if ( c(isub) <= 0.0D+00 ) then

    jsub = iwk18(jval)
    x(jsub) = x(jsub) + wk5(jsub)
    x(isub) = x(isub) + xwk2
    wk10(1:m) = wk11(1:m) - xwk2 * a(1:m,isub)
!
!  A better solution is found.
!
  else
    x(isub) = x(isub) + xwk1
    wk10(1:m) = wk11(1:m) - xwk1 * a(1:m,isub)
    isub = iwk18(jval)
    x(isub) = x(isub) + wk5(isub)
  end if

  more = 1

  return
end
subroutine parse8 ( n, m, a, x, wk1, wk5, wk10, wk11, iwk18, init, &
  jval, kval )

!*****************************************************************************80
!
!! PARSE8 is used by INT_LP.
!
!  Modified:
!
!    11 September 1999
!
!  Author:
!
!    Hang Tong Lau
!
!  Reference:
!
!    Hang Tong Lau,
!    Combinatorial Heuristic Algorithms in FORTRAN,
!    Springer Verlag, 1986,
!    ISBN: 3540171614,
!    LC: QA402.5.L37.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m+1,n+m+1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) init
  integer ( kind = 4 ) isub
  integer ( kind = 4 ) iwk18(n)
  integer ( kind = 4 ) jval
  integer ( kind = 4 ) kval
  real ( kind = 8 ) wk1(n,n)
  real ( kind = 8 ) wk5(n)
  real ( kind = 8 ) wk10(m)
  real ( kind = 8 ) wk11(m)
  real ( kind = 8 ) x(n)
!
!  Consider changing X(JVAL) and a small integer change of X(KVAL).
!
  do

    isub = iwk18(jval)

    if ( x(isub) < wk5(isub) ) then
      exit
    end if

    wk11(1:m) = wk10(1:m) + a(1:m,isub) * wk5(isub)

30  continue

    if ( init < kval ) then
      return
    end if
!
!  Check whether the simultaneous change of X(JVAL) and X(KVAL) is feasible.
!
    isub = iwk18(kval)

    if ( x(isub) + wk5(isub) * wk1(kval,jval) < 0.0D+00 ) then
      kval = kval + 1
      go to 30
    end if

    isub = iwk18(kval)

    do i = 1, m
      if ( wk11(i) < wk5(isub) * wk1(kval,jval) * a(i,isub) ) then
        kval = kval + 1
        go to 30
      end if
    end do
!
!  Make the simultaneous change of X(JVAL) and X(KVAL).
!
    isub = iwk18(jval)
    x(isub) = x(isub) - wk5(isub)
    isub = iwk18(kval)
    x(isub) = x(isub) + wk5(isub) * wk1(kval,jval)

    wk10(1:m) = wk11(1:m) - a(1:m,isub) * wk5(isub) * wk1(kval,jval)

  end do

  return
end
subroutine partition ( n, cost, init, ip, iq, kp, kq, tcost )

!*****************************************************************************80
!
!! PARTITION is a heuristic algorithm for the graph partitioning problem.
!
!  Discussion:
!
!    Let V be the set of 2 * N nodes of a complete graph with an
!    associated 2*N by 2*N symmetric cost matrix C on its edges.
!    The graph partitioning problem is to partion the nodes into
!    two sets P and Q, each containing N nodes, such that we minimize
!    the total cost of the edges in the cut set:
!
!      Sum ( ( I in P ) and ( J in Q ) ) C(I,J)
!
!  Modified:
!
!    11 September 1999
!
!  Author:
!
!    Hang Tong Lau
!
!  Reference:
!
!    Brian Kernighan, Shen Lin,
!    An efficient heuristic procedure for partitioning graphs,
!    Bell System Technical Journal,
!    Volume 49, pages 291-307, 1970.
!
!    Hang Tong Lau,
!    Combinatorial Heuristic Algorithms in FORTRAN,
!    Springer Verlag, 1986,
!    ISBN: 3540171614,
!    LC: QA402.5.L37.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of nodes in each partition.
!
!    Input, real ( kind = 8 ) COST(2*N,2*N), the edge cost matrix.
!
!    Input, integer ( kind = 4 ) INIT:
!    If TRUE, then the initial partition will be generated by this routine;
!    If FALSE, then the initial partition will be supplied on input.
!
!    Input, integer ( kind = 4 ) IP(N), IQ(N).
!    If INIT is TRUE, then IP and IQ are not input quantities, but simply
!    workspace controlled by the routine.
!    If INIT is FALSE, then the user must supply in IP and IQ and initial
!    partition of the nodes.
!
!    Output, integer ( kind = 4 ) KP(N), KQ(N), the final partition of
!    the nodes.
!
!    Output, real ( kind = 8 ) TCOST, the total cost of the edges in the
!    cut set.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) cost(2*n,2*n)
  real ( kind = 8 ) gain
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia
  integer ( kind = 4 ) ib
  integer ( kind = 4 ) ind1
  integer ( kind = 4 ) ind2
  logical init
  integer ( kind = 4 ) ip(n)
  integer ( kind = 4 ) iq(n)
  logical iwk4(n)
  logical iwk5(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kp(n)
  integer ( kind = 4 ) kq(n)
  real ( kind = 8 ) small
  real ( kind = 8 ) tcost
  real ( kind = 8 ) tmax
  real ( kind = 8 ) tot1
  real ( kind = 8 ) tot2
  real ( kind = 8 ) wk1(n)
  real ( kind = 8 ) wk2(n)
  real ( kind = 8 ) wk3(n)
!
!  Initial partitioning.
!
  if ( init ) then
    do i = 1, n
      ip(i) = i
      iq(i) = i + n
    end do
  end if
!
!  Set the flags.
!
  do

    iwk4(1:n) = .true.
    iwk5(1:n) = .true.

    tcost = 0.0D+00
    do i = 1, n
      do j = 1, n
        tcost = tcost + cost(ip(i),iq(j))
      end do
    end do

    small = - 2.0D+00 * tcost
!
!  Calculate the external cost of each element in the first partition.
!
    do i = 1, n

      tot1 = 0.0D+00
      do j = 1, n
        tot1 = tot1 + cost(ip(i),iq(j))
      end do
!
!  Calculate the internal cost of each element in the first partition.
!
      tot2 = 0.0D+00
      do k = 1, n
        tot2 = tot2 + cost(ip(i),ip(k))
      end do
!
!  Calculate the difference between external and internal costs.
!
      wk1(i) = tot1 - tot2

    end do

    do i = 1, n
!
!  Calculate the external cost of each element in the second partition.
!
      tot1 = 0.0D+00
      do j = 1, n
        tot1 = tot1 + cost(iq(i),ip(j))
      end do
!
!  Calculate the internal cost of each element in the second partition.
!
      tot2 = 0.0D+00
      do k = 1, n
        tot2 = tot2 + cost(iq(i),iq(k))
      end do
!
!  Calculate the difference between external and internal costs.
!
      wk2(i) = tot1 - tot2

    end do

    do i = 1, n
!
!  Choose IA from the first partition and IB from the second
!  partition such that the gain is maximum.
!
      tmax = small

      do j = 1, n

        if ( iwk4(j) )   then

          do k =  1, n

            if ( iwk5(k) ) then

              gain = wk1(j) + wk2(k) - 2.0D+00 * cost(ip(j),iq(k))

              if ( tmax < gain ) then
                tmax = gain
                ia = ip(j)
                ib = iq(k)
                ind1 = j
                ind2 = k
              end if

            end if

          end do

        end if

      end do

      wk3(i) = tmax
      kp(i) = ia
      kq(i) = ib
      iwk4(ind1) = .false.
      iwk5(ind2) = .false.
!
!  Recalculate the cost differences.
!
      do j = 1, n

        if ( iwk4(j) ) then
          wk1(j) = wk1(j) + 2.0D+00 * cost(ip(j),ia) - 2.0*cost(ip(j),ib)
        end if

        if ( iwk5(j) ) then
          wk2(i) = wk2(j) + 2.0D+00 * cost(iq(j),ib) - 2.0 * cost(iq(j),ia)
        end if

      end do

    end do
!
!  Choose K such that WK3(K) is maximal.
!
    tmax = small

    do i = 1, n

      tot1 = 0.0D+00
      do j = 1, i
        tot1 = tot1 + wk3(j)
      end do

      if ( tmax < tot1 ) then
        tmax = tot1
        k = i
      end if

    end do
!
!  Exchange the two elements found above.
!  Iterate until no reduction in cost can be obtained.
!
    if ( tmax <= 0.0D+00 ) then
      exit
    end if

    do i = 1, n

      if ( i <= k ) then
        ip(i) = kq(i)
        iq(i) = kp(i)
      else
        ip(i) = kp(i)
        iq(i) = kq(i)
      end if

    end do

  end do

  return
end
subroutine pmatch ( n, cost, pair )

!*****************************************************************************80
!
!! PMATCH finds a minimum weight perfect matching in a graph.
!
!  Modified:
!
!    15 September 1999
!
!  Author:
!
!    Hang Tong Lau
!
!  Reference:
!
!    Rainer Burkard, Ulrich Derigs,
!    Assignment and Matching Problems: Solution methods with
!    FORTRAN programs,
!    Lecture Notes in Economics and Mathematical Systems,
!    Volume 184,
!    Springer, 1980,
!    ISBN: 0387102671,
!    LC: QA402.5.B86.
!
!    Hang Tong Lau,
!    Combinatorial Heuristic Algorithms in FORTRAN,
!    Springer Verlag, 1986,
!    ISBN: 3540171614,
!    LC: QA402.5.L37.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of nodes in the complete graph.
!    N is assumed to be even.
!
!    Input, real ( kind = 8 ) COST(N*(N-1)/2), the strict upper triangle of the
!    cost matrix, stored by rows.  In other words, the first elements of
!    COST are the costs C(1,2), C(1,3), ..., C(1,N), C(2,3), C(2,4),
!    ..., C(2,N).
!
!    Output, integer ( kind = 4 ) PAIR(N), contains the minimum weight perfect
!    matching.  Node I is connected to node PAIR(I), for I = 1 to N.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) cost((n*(n-1))/2)
  real ( kind = 8 ) cst
  real ( kind = 8 ) cstlow
  real ( kind = 8 ) cswk
  real ( kind = 8 ) cwk2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihead
  integer ( kind = 4 ) index
  integer ( kind = 4 ) isub
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jwk1(n)
  integer ( kind = 4 ) jwk2(n)
  integer ( kind = 4 ) jwk3(n)
  integer ( kind = 4 ) jwk4(n)
  integer ( kind = 4 ) jwk5(n)
  integer ( kind = 4 ) jwk6(n)
  integer ( kind = 4 ) jwk7(n)
  integer ( kind = 4 ) jwk8(n)
  integer ( kind = 4 ) jwk9(n)
  integer ( kind = 4 ) kk1
  integer ( kind = 4 ) kk2
  integer ( kind = 4 ) kk3
  integer ( kind = 4 ) kk4
  integer ( kind = 4 ) kk5
  integer ( kind = 4 ) kk6
  integer ( kind = 4 ) ll1
  integer ( kind = 4 ) ll2
  integer ( kind = 4 ) ll3
  integer ( kind = 4 ) ll4
  integer ( kind = 4 ) ll5
  integer ( kind = 4 ) max
  integer ( kind = 4 ) min
  integer ( kind = 4 ) mm1
  integer ( kind = 4 ) mm2
  integer ( kind = 4 ) mm3
  integer ( kind = 4 ) mm4
  integer ( kind = 4 ) mm5
  integer ( kind = 4 ) nn
  integer ( kind = 4 ) nn2
  integer ( kind = 4 ) pair(n)
  real ( kind = 8 ) xcst
  real ( kind = 8 ) value
  real ( kind = 8 ) work1(n)
  real ( kind = 8 ) work2(n)
  real ( kind = 8 ) work3(n)
  real ( kind = 8 ) work4(n)
  real ( kind = 8 ) xwk2
  real ( kind = 8 ) xwk3
  real ( kind = 8 ) xwork

  if ( mod ( n, 2 ) /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PMATCH - Fatal error!'
    write ( *, '(a)' ) '  The value of N must be even.'
    write ( *, '(a,i8)' ) '  The input value was N = ', n
    stop
  end if
!
!  Initialization
!
  nn2 = ( n * ( n - 1 ) ) / 2

  jwk1(2) = 0

  do i = 3, n
    jwk1(i) = jwk1(i-1) + i - 2
  end do

  ihead = n + 2

  do i = 1, n
    jwk2(i) = i
    jwk3(i) = i
    jwk5(i) = i
  end do

  jwk4(1:n) = 0
  jwk6(1:n) = ihead
  jwk7(1:n) = ihead
  jwk8(1:n) = ihead
  pair(1:n) = ihead
  work1(1:n) = huge ( 1.0D+00 )
  work2(1:n) = 0.0D+00
  work3(1:n) = 0.0D+00
  work4(1:n) = huge ( 1.0D+00 )
!
!  Start procedure.
!
  do i = 1, n

    if ( pair(i) == ihead ) then

      nn = 0
      cwk2 = huge ( 1.0D+00 )

      do j = 1, n

        min = i
        max = j

        if ( i /= j ) then

          if ( j < i ) then
            max = i
            min = j
          end if

          isub = jwk1(max) + min
          xcst = cost(isub)
          cswk = cost(isub) - work2(j)

          if ( cswk <= cwk2 ) then

            if ( cswk == cwk2 ) then
              if ( nn == 0 ) then
                go to 30
              else
                cycle
              end if
            end if

            cwk2 = cswk

            nn = 0

30          continue

            if ( pair(j) == ihead ) then
              nn = j
            end if

          end if

        end if

      end do

      if ( nn /= 0 ) then
        work2(i) = cwk2
        pair(i) = nn
        pair(nn) = i
      end if

    end if

  end do
!
!  Initial labeling.
!
  nn = 0

  do i = 1, n

    if ( pair(i) == ihead ) then

      nn = nn + 1
      jwk6(i) = 0
      work4(i) = 0.0D+00
      xwk2 = work2(i)

      do j = 1, n

        min = i
        max = j

        if ( i /= j ) then

          if ( j < i ) then
            max = i
            min = j
          end if

          isub = jwk1(max) + min
          xcst = cost(isub)
          cswk = cost(isub) - xwk2 - work2(j)

          if ( cswk < work1(j) ) then
            work1(j) = cswk
            jwk4(j) = i
          end if

        end if

      end do

    end if

  end do

  if ( nn <= 1 ) then
    go to 340
  end if
!
!  Examine the labeling and prepare for the next step.
!
80 continue

  cstlow = huge ( 1.0D+00 )

  do i = 1, n

    if ( jwk2(i) == i ) then

      cst = work1(i)

      if ( jwk6(i) < ihead ) then

        cst = 0.5D+00 * ( cst + work4(i) )

        if ( cst <= cstlow ) then
          index = i
          cstlow = cst
        end if

      else

        if ( jwk7(i) < ihead ) then
          if ( jwk3(i) /= i ) then
            cst = cst + work2(i)
            if ( cst < cstlow ) then
              index = i
              cstlow = cst
            end if
          end if
        else
          if ( cst < cstlow ) then
            index = i
            cstlow = cst
          end if
        end if

      end if

    end if

  end do

  if  ( jwk7(index) < ihead ) then
    go to 190
  end if

  if  ( jwk6(index) < ihead ) then

    ll4 = jwk4(index)
    ll5 = jwk5(index)
    kk4 = index
    kk1 = kk4
    kk5 = jwk2(ll4)
    kk2 = kk5

    do

      jwk7(kk1) = kk2
      mm5 = jwk6(kk1)

      if ( mm5 == 0 ) then
        exit
      end if

      kk2 = jwk2(mm5)
      kk1 = jwk7(kk2)
      kk1 = jwk2(kk1)

    end do

    ll2 = kk1
    kk1 = kk5
    kk2 = kk4

110 continue

    if ( ihead <= jwk7(kk1) ) then

      jwk7(kk1) = kk2
      mm5 = jwk6(kk1)

      if ( mm5 == 0 ) then
        go to 280
      end if

      kk2 = jwk2(mm5)
      kk1 = jwk7(kk2)
      kk1 = jwk2(kk1)
      go to 110

    end if

120 continue

    if ( kk1 == ll2 ) then
      go to 130
    end if

    mm5 = jwk7(ll2)
    jwk7(ll2) = ihead
    ll1 = pair(mm5)
    ll2 = jwk2(ll1)
    go to 120

  end if
!
!  Growing an alternating tree, add two edges.
!
  jwk7(index) = jwk4(index)
  jwk8(index) = jwk5(index)
  ll1 = pair(index)
  ll3 = jwk2(ll1)
  work4(ll3) = cstlow
  jwk6(ll3) = pair(ll3)

  call pmatch_sub_b ( ll3, n, nn2, cost, jwk1, jwk2, jwk3, jwk4, &
    jwk5, jwk7, jwk9, work1, work2, work3, work4 )

  go to 80
!
!  Shrink a blossom.
!
130 continue

  xwork = work2(ll2) + cstlow - work4(ll2)
  work2(ll2) = 0.0D+00
  mm1 = ll2

 do

    work3(mm1) = work3(mm1) + xwork
    mm1 = jwk3(mm1)

    if ( mm1 == ll2 ) then
      exit
    end if

  end do

  mm5 = jwk3(ll2)

  if ( ll2 /= kk5 ) then
    go to 160
  end if

150 continue

  kk5 = kk4
  kk2 = jwk7(ll2)

160 continue

  jwk3(mm1) = kk2
  ll1 = pair(kk2)
  jwk6(kk2) = ll1
  xwk2 = work2(kk2) + work1(kk2) - cstlow
  mm1 = kk2

  do

    mm2 = mm1
    work3(mm2) = work3(mm2) + xwk2
    jwk2(mm2) = ll2
    mm1 = jwk3(mm2)

    if ( mm1 == kk2 ) then
      exit
    end if

  end do

  jwk5(kk2) = mm2
  work2(kk2) = xwk2
  kk1 = jwk2(ll1)
  jwk3(mm2) = kk1
  xwk2 = work2(kk1) + cstlow - work4(kk1)
  mm2 = kk1

  do

    mm1 = mm2
    work3(mm1) = work3(mm1) + xwk2
    jwk2(mm1) = ll2
    mm2 = jwk3(mm1)

    if ( mm2 == kk1 ) then
      exit
    end if

  end do

  jwk5(kk1) = mm1
  work2(kk1) = xwk2

  if  ( kk5 /= kk1 ) then
    kk2 = jwk7(kk1)
    jwk7(kk1) = jwk8(kk2)
    jwk8(kk1) = jwk7(kk2)
    go to 160
  end if

  if ( kk5 /= index ) then
    jwk7(kk5) = ll5
    jwk8(kk5) = ll4
    if ( ll2 /= index ) then
      go to 150
    end if
  else
    jwk7(index) = ll4
    jwk8(index) = ll5
  end if

  jwk3(mm1) = mm5
  kk4 = jwk3(ll2)
  jwk4(kk4) = mm5
  work4(kk4) = xwork
  jwk7(ll2) = ihead
  work4(ll2) = cstlow

  call pmatch_sub_b ( ll2, n, nn2, cost, jwk1, jwk2, jwk3, jwk4, &
    jwk5, jwk7, jwk9, work1, work2, work3, work4 )

  go to 80
!
!  Expand a T-labeled blossom.
!
190   continue

  kk4 = jwk3(index)
  kk3 = kk4
  ll4 = jwk4(kk4)
  mm2 = kk4

  do

    mm1 = mm2
    ll5 = jwk5(mm1)
    xwk2 = work2(mm1)

    do

      jwk2(mm2) = mm1
      work3(mm2)= work3(mm2) - xwk2

      if ( mm2 == ll5 ) then
        exit
      end if

      mm2 = jwk3(mm2)

    end do

    mm2 = jwk3(ll5)
    jwk3(ll5) = mm1

    if ( mm2 == ll4 ) then
      exit
    end if

  end do

  xwk2 = work4(kk4)
  work2(index) = xwk2
  jwk3(index) = ll4
  mm2 = ll4

  do

    work3(mm2) = work3(mm2) - xwk2

    if ( mm2 == index ) then
      exit
    end if

    mm2 = jwk3(mm2)

  end do

  mm1 = pair(index)
  kk1 = jwk2(mm1)
  mm2 = jwk6(kk1)
  ll2 = jwk2(mm2)

  if ( ll2 /= index ) then

    kk2 = ll2

    do

      mm5 = jwk7(kk2)
      kk1 = jwk2(mm5)

      if ( kk1 == index ) then
        exit
      end if

      kk2 = jwk6(kk1)
      kk2 = jwk2(kk2)

    end do

    jwk7(ll2) = jwk7(index)
    jwk7(index) = jwk8(kk2)
    jwk8(ll2) = jwk8(index)
    jwk8(index) = mm5
    mm3 = jwk6(ll2)
    kk3 = jwk2(mm3)
    mm4 = jwk6(kk3)
    jwk6(ll2) = ihead
    pair(ll2) = mm1
    kk1 = kk3

    do

      mm1 = jwk7(kk1)
      mm2 = jwk8(kk1)
      jwk7(kk1) = mm4
      jwk8(kk1) = mm3
      jwk6(kk1) = mm1
      pair(kk1) = mm1
      kk2 = jwk2(mm1)
      pair(kk2) = mm2
      mm3 = jwk6(kk2)
      jwk6(kk2) = mm2

      if ( kk2 == index ) then
        exit
      end if

      kk1 = jwk2(mm3)
      mm4 = jwk6(kk1)
      jwk7(kk2) = mm3
      jwk8(kk2) = mm4

    end do

  end if

  mm2 = jwk8(ll2)
  kk1 = jwk2(mm2)
  work1(kk1) = cstlow
  kk4 = 0

  if ( kk1 /= ll2 ) then

    mm1 = jwk7(kk1)
    kk3 = jwk2(mm1)
    jwk7(kk1) = jwk7(ll2)
    jwk8(kk1) = mm2

    do

      mm5 = jwk6(kk1)
      jwk6(kk1) = ihead
      kk2 = jwk2(mm5)
      mm5 = jwk7(kk2)
      jwk7(kk2) = ihead
      kk5 = jwk8(kk2)
      jwk8(kk2) = kk4
      kk4 = kk2
      work4(kk2) = cstlow
      kk1 = jwk2(mm5)
      work1(kk1) = cstlow

      if ( kk1 == ll2 ) then
        exit
      end if

    end do

    jwk7(ll2) = kk5
    jwk8(ll2) = mm5
    jwk6(ll2) = ihead

    if ( kk3 == ll2 ) then
      go to 270
    end if

  end if

  kk1 = 0
  kk2 = kk3

  do

    mm5 = jwk6(kk2)
    jwk6(kk2) = ihead
    jwk7(kk2) = ihead
    jwk8(kk2) = kk1
    kk1 = jwk2(mm5)
    mm5 = jwk7(kk1)
    jwk6(kk1) = ihead
    jwk7(kk1) = ihead
    jwk8(kk1) = kk2
    kk2 = jwk2(mm5)

    if ( kk2 == ll2 ) then
      exit
    end if

  end do

  call pmatch_sub_a ( kk1, n, nn2, cost, jwk1, jwk2, jwk3, jwk4, &
    jwk5, jwk6, jwk8, work1, work2, work3, work4 )

270  continue

  if ( kk4 == 0 ) then
    go to 80
  end if

  ll2 = kk4

  call pmatch_sub_b ( ll2, n, nn2, cost, jwk1, jwk2, jwk3, jwk4, &
    jwk5, jwk7, jwk9, work1, work2, work3, work4 )

  kk4 = jwk8(ll2)
  jwk8(ll2) = ihead
  go to 270
!
!  Augmentation of the matching.
!
!  Exchange the matching and non-matching edges along the augmenting path.
!
 280  continue

  ll2 = kk4
  mm5 = ll4

  do

    kk1 = ll2

    do

      pair(kk1) = mm5
      mm5 = jwk6(kk1)
      jwk7(kk1) = ihead

      if ( mm5 == 0 ) then
        exit
      end if

      kk2 = jwk2(mm5)
      mm1 = jwk7(kk2)
      mm5 = jwk8(kk2)
      kk1 = jwk2(mm1)
      pair(kk2) = mm1

    end do

    if ( ll2 /= kk4 ) then
      exit
    end if

    ll2 = kk5
    mm5 = ll5

  end do
!
!  Remove all labels of non-exposed base nodes.
!
  do i = 1, n

    if ( jwk2(i) == i ) then

      if ( jwk6(i) < ihead ) then

        cst = cstlow - work4(i)
        work2(i) = work2(i) + cst
        jwk6(i) = ihead

        if ( pair(i) /= ihead ) then
          work4(i) = huge ( 1.0D+00 )
        else
          jwk6(i) = 0
          work4(i) = 0.0D+00
        end if

      else

        if ( jwk7(i) < ihead ) then
          cst = work1(i) - cstlow
          work2(i) = work2(i) + cst
          jwk7(i) = ihead
          jwk8(i) = ihead
        end if

        work4(i) = huge ( 1.0D+00 )

      end if

      work1(i) = huge ( 1.0D+00 )

    end if

  end do

  nn = nn - 2
!
!  Determine the new WORK1 values.
!
  if ( nn <= 1 ) then
    go to 340
  end if

  do i = 1, n

    kk1 = jwk2(i)

    if ( jwk6(kk1) == 0 ) then

      xwk2 = work2(kk1)
      xwk3 = work3(i)

      do j = 1, n

        kk2 = jwk2(j)

        if ( kk1 /= kk2 ) then

          min = i
          max = j

          if ( i /= j ) then

            if ( j < i ) then
              max = i
              min = j
            end if

            isub = jwk1(max) + min
            xcst = cost(isub)
            cswk = cost(isub) - xwk2 - xwk3
            cswk = cswk - work2(kk2) - work3(j)

            if ( cswk < work1(kk2) ) then
              jwk4(kk2) = i
              jwk5(kk2) = j
              work1(kk2) = cswk
            end if

          end if
        end if

      end do

    end if

  end do

  go to 80
!
!  Generate the original graph by expanding all shrunken blossoms.
!
340   continue

  value = 0.0D+00

  do i = 1, n

    if ( jwk2(i) == i ) then

      if ( 0 <= jwk6(i) ) then

        kk5 = pair(i)
        kk2 = jwk2(kk5)
        kk4 = pair(kk2)
        jwk6(i) = -1
        jwk6(kk2) = -1
        min = kk4
        max = kk5

        if ( kk4 /= kk5 ) then

          if ( kk5 < kk4 ) then
            max = kk4
            min = kk5
          end if

          isub = jwk1(max) + min
          xcst = cost(isub)
          value = value + xcst
        end if

      end if

    end if

  end do

  do i = 1, n

360 continue

    ll2 = jwk2(i)

    if ( ll2 == i ) then
      cycle
    end if

    mm2 = jwk3(ll2)
    ll4 = jwk4(mm2)
    kk3 = mm2
    xwork = work4(mm2)

    do

      mm1 = mm2
      ll5 = jwk5(mm1)
      xwk2 = work2(mm1)

      do

        jwk2(mm2) = mm1
        work3(mm2) = work3(mm2) - xwk2

        if ( mm2 == ll5 ) then
          exit
        end if

        mm2 = jwk3(mm2)

      end do

      mm2 = jwk3(ll5)
      jwk3(ll5) = mm1

      if ( mm2 == ll4 ) then
        exit
      end if

    end do

    work2(ll2) = xwork
    jwk3(ll2) = ll4
    mm2 = ll4

    do

      work3(mm2) = work3(mm2) - xwork

      if ( mm2 == ll2 ) then
        exit
      end if

      mm2 = jwk3(mm2)

    end do

    mm5 = pair(ll2)
    mm1 = jwk2(mm5)
    mm1 = pair(mm1)
    kk1 = jwk2(mm1)

    if ( ll2 /= kk1 ) then

      pair(kk1) = mm5
      kk3 = jwk7(kk1)
      kk3 = jwk2(kk3)

      do

        mm3 = jwk6(kk1)
        kk2 = jwk2(mm3)
        mm1 = jwk7(kk2)
        mm2 = jwk8(kk2)
        kk1 = jwk2(mm1)
        pair(kk1) = mm2
        pair(kk2) = mm1
        min = mm1
        max = mm2

        if ( mm1 == mm2 ) then
          go to 360
        end if

        if ( mm2 < mm1 ) then
          max = mm1
          min = mm2
        end if

        isub = jwk1(max) + min
        xcst = cost(isub)
        value = value + xcst

        if ( kk1 == ll2 ) then
          exit
        end if

      end do

      if ( kk3 == ll2 ) then
        go to 360
      end if

    end if

    do

      kk5 = jwk6(kk3)
      kk2 = jwk2(kk5)
      kk6 = jwk6(kk2)
      min = kk5
      max = kk6

      if ( kk5 == kk6 ) then
        go to 360
      end if

      if ( kk6 < kk5 ) then
        max = kk5
        min = kk6
      end if

      isub = jwk1(max) + min
      xcst = cost(isub)
      value = value + xcst
      kk6 = jwk7(kk2)
      kk3 = jwk2(kk6)

      if ( kk3 == ll2 ) then
        go to 360
      end if

    end do

  end do

  return
end
subroutine pmatch_sub_a ( kk, n, nn2, cost, jwk1, jwk2, jwk3, jwk4, &
  jwk5, jwk6, jwk8, work1, work2, work3, work4 )

!*****************************************************************************80
!
!! PMATCH_SUB_A is used by PMATCH.
!
!  Modified:
!
!    11 September 1999
!
!  Author:
!
!    Hang Tong Lau
!
!  Reference:
!
!    Hang Tong Lau,
!    Combinatorial Heuristic Algorithms in FORTRAN,
!    Springer Verlag, 1986,
!    ISBN: 3540171614,
!    LC: QA402.5.L37.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nn2

  real ( kind = 8 ) cost(nn2)
  real ( kind = 8 ) cstwk
  real ( kind = 8 ) cswk
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihead
  integer ( kind = 4 ) isub
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj1
  integer ( kind = 4 ) jj2
  integer ( kind = 4 ) jj3
  integer ( kind = 4 ) jj4
  integer ( kind = 4 ) jwk1(n)
  integer ( kind = 4 ) jwk2(n)
  integer ( kind = 4 ) jwk3(n)
  integer ( kind = 4 ) jwk4(n)
  integer ( kind = 4 ) jwk5(n)
  integer ( kind = 4 ) jwk6(n)
  integer ( kind = 4 ) jwk8(n)
  integer ( kind = 4 ) kk
  integer ( kind = 4 ) max
  integer ( kind = 4 ) min
  real ( kind = 8 ) work1(n)
  real ( kind = 8 ) work2(n)
  real ( kind = 8 ) work3(n)
  real ( kind = 8 ) work4(n)
  real ( kind = 8 ) xwk2
  real ( kind = 8 ) xwk3

  ihead = n + 2

  do

    jj1 = kk
    kk  = jwk8(jj1)
    jwk8(jj1) = ihead
    cstwk = huge ( 1.0D+00 )
    jj3 = 0
    jj4 = 0
    j = jj1
    xwk2 = work2(jj1)

    do

      xwk3 = work3(j)

      do i = 1, n

        jj2 = jwk2(i)

        if ( jwk6(jj2) < ihead ) then

          min = j
          max = i

          if ( j /= i ) then

            if ( i < j ) then
              max = j
              min = i
            end if

            isub = jwk1(max) + min

            cswk = cost(isub) - xwk2 - xwk3 - work2(jj2) - work3(i) + work4(jj2)

            if ( cswk < cstwk ) then
              jj3 = i
              jj4 = j
              cstwk = cswk
            end if

          end if

        end if

      end do

      j = jwk3(j)

      if ( j == jj1 ) then
        exit
      end if

    end do

    jwk4(jj1) = jj3
    jwk5(jj1) = jj4
    work1(jj1) = cstwk

    if ( kk == 0 ) then
      exit
    end if

  end do

  return
end
subroutine pmatch_sub_b ( kk, n, nn2, cost, jwk1, jwk2, jwk3, jwk4, &
  jwk5, jwk7, jwk9, work1, work2, work3, work4 )

!*****************************************************************************80
!
!! PMATCH_SUB_B is used by PMATCH.
!
!  Modified:
!
!    11 September 1999
!
!  Author:
!
!    Hang Tong Lau
!
!  Reference:
!
!    Hang Tong Lau,
!    Combinatorial Heuristic Algorithms in FORTRAN,
!    Springer Verlag, 1986,
!    ISBN: 3540171614,
!    LC: QA402.5.L37.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nn2

  real ( kind = 8 ) cost(nn2)
  real ( kind = 8 ) cswk
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihead
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) isub
  integer ( kind = 4 ) jj1
  integer ( kind = 4 ) jj2
  integer ( kind = 4 ) jj3
  integer ( kind = 4 ) jwk1(n)
  integer ( kind = 4 ) jwk2(n)
  integer ( kind = 4 ) jwk3(n)
  integer ( kind = 4 ) jwk4(n)
  integer ( kind = 4 ) jwk5(n)
  integer ( kind = 4 ) jwk7(n)
  integer ( kind = 4 ) jwk9(n)
  integer ( kind = 4 ) kk
  integer ( kind = 4 ) max
  integer ( kind = 4 ) min
  real ( kind = 8 ) work1(n)
  real ( kind = 8 ) work2(n)
  real ( kind = 8 ) work3(n)
  real ( kind = 8 ) work4(n)
  real ( kind = 8 ) xwk1
  real ( kind = 8 ) xwk2

  ihead = n + 2
  xwk1 = work4(kk) - work2(kk)
  work1(kk) = huge ( 1.0D+00 )
  xwk2 = xwk1 - work3(kk)
  jwk7(kk) = 0
  ii = 0

  do i = 1, n

    jj3 = jwk2(i)

    if ( ihead <= jwk7(jj3) ) then

      ii = ii + 1
      jwk9(ii) = i
      min = kk
      max = i

      if ( kk /= i ) then

        if ( i < kk ) then
          max = kk
          min = i
        end if

        isub = jwk1(max) + min

        cswk = cost(isub) + xwk2 - work2(jj3) - work3(i)

        if ( cswk < work1(jj3) ) then
          jwk4(jj3) = kk
          jwk5(jj3) = i
          work1(jj3) = cswk
        end if
      end if
    end if

  end do

  jwk7(kk) = ihead
  jj1 = kk
  jj1 = jwk3(jj1)

  if ( jj1 == kk ) then
    return
  end if

  do

    xwk2 = xwk1 - work3(jj1)

    do i = 1, ii

      jj2 = jwk9(i)
      jj3 = jwk2(jj2)
      min = jj1
      max = jj2

      if ( jj1 /= jj2 ) then

        if ( jj2 < jj1 ) then
          max = jj1
          min = jj2
        end if

        isub = jwk1(max) + min

        cswk = cost(isub) + xwk2 - work2(jj3) - work3(jj2)

        if ( cswk < work1(jj3) ) then
          jwk4(jj3) = jj1
          jwk5(jj3) = jj2
          work1(jj3) = cswk
        end if

      end if

    end do

    jj1 = jwk3(jj1)

    if ( jj1 == kk ) then
      exit
    end if

  end do

  return
end
subroutine r8_swap ( x, y )

!*****************************************************************************80
!
!! R8_SWAP switches two R8's.
!
!  Modified:
!
!    30 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X, Y.  On output, the values of X and
!    Y have been interchanged.
!
  implicit none

  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  z = x
  x = y
  y = z

  return
end
subroutine r8mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8MAT_PRINT prints an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an array of R8 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = * ) title

  call r8mat_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine r8mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8MAT_PRINT_SOME prints some of an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an array of R8 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ( kind = 4 ) ILO, JLO, the first row and column to print.
!
!    Input, integer ( kind = 4 ) IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = 14 ) ctemp(incx)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i8,6x)' ) j
    end do

    write ( *, '(''  Col   '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) ' '

    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi

      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( a(i,j) == real ( int ( a(i,j) ), kind = 8 ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) a(i,j)
        else
          write ( ctemp(j2), '(g14.6)' ) a(i,j)
        end if

      end do

      write ( *, '(i5,a,5a14)' ) i, ':', ( ctemp(j), j = 1, inc )

    end do

  end do

  return
end
subroutine r8vec_sort_heap_index_a ( n, a, indx )

!*****************************************************************************80
!
!! R8VEC_SORT_HEAP_INDEX_A does an indexed heap ascending sort of an R8VEC.
!
!  Discussion:
!
!    The sorting is not actually carried out.  Rather an index array is
!    created which defines the sorting.  This array may be used to sort
!    or index the array, or to sort or index related arrays keyed on the
!    original array.
!
!    Once the index array is computed, the sorting can be carried out
!    "implicitly:
!
!      A(INDX(I)), I = 1 to N is sorted,
!
!    or explicitly, by the call
!
!      call R8VEC_PERMUTE ( N, A, INDX )
!
!    after which A(I), I = 1 to N is sorted.
!
!  Modified:
!
!    17 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, real ( kind = 8 ) A(N), an array to be index-sorted.
!
!    Output, integer ( kind = 4 ) INDX(N), contains the sort index.  The
!    I-th element of the sorted array is A(INDX(I)).
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) aval
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) indxt
  integer ( kind = 4 ) ir
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l

  do i = 1, n
    indx(i) = i
  end do

  l = n / 2 + 1
  ir = n

  do

    if ( 1 < l ) then

      l = l - 1
      indxt = indx(l)
      aval = a(indxt)

    else

      indxt = indx(ir)
      aval = a(indxt)
      indx(ir) = indx(1)
      ir = ir - 1

      if ( ir == 1 ) then
        indx(1) = indxt
        exit
      end if

    end if

    i = l
    j = l + l

    do while ( j <= ir )

      if ( j < ir ) then
        if ( a(indx(j)) < a(indx(j+1)) ) then
          j = j + 1
        end if
      end if

      if ( aval < a(indx(j)) ) then
        indx(i) = indx(j)
        i = j
        j = j + j
      else
        j = ir + 1
      end if

    end do

    indx(i) = indxt

  end do

  return
end
subroutine r8vec_sort_heap_index_d ( n, a, indx )

!*****************************************************************************80
!
!! R8VEC_SORT_HEAP_INDEX_D does an indexed heap descending sort of an R8VEC.
!
!  Discussion:
!
!    The sorting is not actually carried out.  Rather an index array is
!    created which defines the sorting.  This array may be used to sort
!    or index the array, or to sort or index related arrays keyed on the
!    original array.
!
!    Once the index array is computed, the sorting can be carried out
!    "implicitly:
!
!      A(INDX(I)), I = 1 to N is sorted,
!
!    or explicitly, by the call
!
!      call R8VEC_PERMUTE ( N, A, INDX )
!
!    after which A(I), I = 1 to N is sorted.
!
!  Modified:
!
!    07 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, real ( kind = 8 ) A(N), an array to be index-sorted.
!
!    Output, integer ( kind = 4 ) INDX(N), contains the sort index.  The
!    I-th element of the sorted array is A(INDX(I)).
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) aval
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) indxt
  integer ( kind = 4 ) ir
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l

  do i = 1, n
    indx(i) = i
  end do

  l = n / 2 + 1
  ir = n

  do

    if ( 1 < l ) then

      l = l - 1
      indxt = indx(l)
      aval = a(indxt)

    else

      indxt = indx(ir)
      aval = a(indxt)
      indx(ir) = indx(1)
      ir = ir - 1

      if ( ir == 1 ) then
        indx(1) = indxt
        exit
      end if

    end if

    i = l
    j = l + l

    do while ( j <= ir )

      if ( j < ir ) then
        if ( a(indx(j+1)) < a(indx(j)) ) then
          j = j + 1
        end if
      end if

      if ( a(indx(j)) < aval ) then
        indx(i) = indx(j)
        i = j
        j = j + j
      else
        j = ir + 1
      end if

    end do

    indx(i) = indxt

  end do

  return
end
subroutine simplex ( n, m, a, x, infs, wk2, wk13, wk14, wk16, iwk19, iwk20 )

!*****************************************************************************80
!
!! SIMPLEX is used by INT_LP.
!
!  Modified:
!
!    13 September 2000
!
!  Author:
!
!    Hang Tong Lau
!
!  Reference:
!
!    Hang Tong Lau,
!    Combinatorial Heuristic Algorithms in FORTRAN,
!    Springer Verlag, 1986,
!    ISBN: 3540171614,
!    LC: QA402.5.L37.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, integer ( kind = 4 ) M, the number of constraints.
!
!    Input, real ( kind = 8 ) A((M+1)*(N+M+1)), contains the coefficients of the
!    constraints in the first M rows and N columns.  There is also an M+1 row
!    used for working storage.  The last M+1 columns are also used as working
!    storage.
!
!    Output, real ( kind = 8 ) X(max ( M+1, N ) ), contains the solution.
!
!    Output, integer ( kind = 4 ) INFS, error indicator.
!    0, no error, a solution was found.
!    1, the linear program is infeasible.  No solution was found.
!
!    Workspace, real ( kind = 8 ) WK2((M+1)*(M+1)).
!
!    Workspace, real ( kind = 8 ) WK13(M+1).
!
!    Workspace, real ( kind = 8 ) WK14(M+1).
!
!    Workspace, real ( kind = 8 ) WK16(max(M+1,N)).
!
!    Workspace, integer ( kind = 4 ) IWK19(M+1).
!
!    Output, integer ( kind = 4 ) IWK20(N+M).
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a((m+1)*(n+m+1))
  real ( kind = 8 ) beta
  real ( kind = 8 ) delta1
  real ( kind = 8 ) delta2
  real ( kind = 8 ) delta3
  real ( kind = 8 ) delta4
  real ( kind = 8 ) delta5
  real ( kind = 8 ) delta6
  real ( kind = 8 ) delta7
  real ( kind = 8 ) delta8
  real ( kind = 8 ) delta9
  real ( kind = 8 ) eps
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic1
  integer ( kind = 4 ) ic2
  integer ( kind = 4 ) ic3
  integer ( kind = 4 ) ic4
  integer ( kind = 4 ) ic5
  integer ( kind = 4 ) ic6
  integer ( kind = 4 ) ic7
  integer ( kind = 4 ) ic8
  integer ( kind = 4 ) ic9
  integer ( kind = 4 ) inc1
  integer ( kind = 4 ) inc2
  integer ( kind = 4 ) inc3
  integer ( kind = 4 ) inc4
  integer ( kind = 4 ) infs
  integer ( kind = 4 ) ispec
  integer ( kind = 4 ) it
  integer ( kind = 4 ) iter
  integer ( kind = 4 ) iwk19(m+1)
  integer ( kind = 4 ) iwk20(n+m)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jc1
  integer ( kind = 4 ) jc2
  integer ( kind = 4 ) jc3
  integer ( kind = 4 ) jc4
  integer ( kind = 4 ) jc5
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) ll
  integer ( kind = 4 ) loop1
  integer ( kind = 4 ) loop2
  integer ( kind = 4 ) loop3
  integer ( kind = 4 ), parameter :: maxit = 900
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) m99
  integer ( kind = 4 ) n99
  real ( kind = 8 ) wk2(*)
  real ( kind = 8 ) wk13(m+1)
  real ( kind = 8 ) wk14(m+1)
  real ( kind = 8 ) wk16(n)
  real ( kind = 8 ) x(n)
!
!  Compute the machine epsilon.
!
  eps = sqrt ( epsilon ( 1.0D+00 ) )
!
!  Simplex method for linear programming.
!
  n99 = n + m
  ic5 = m + 1
  m99 = m + 1
  ic6 = 2
  ic7 = 1
  inc1 = 0
  k = 0
  iter = 0
  inc2 = 0
  ic1 = 0
  jc5 = 0
  delta1 = eps
  delta2 = eps
  delta3 = - eps * 2.0D+00
  jc1 = 0
  delta4 = 0.0D+00
  delta6 = 0.0D+00
  m2 = m99 + m99
  infs = 1
  ispec = 7777
!
!  Start phase one.
!
  iwk19(1:m99) = 0

  ic9 = 0

  do j = 1, n99

    iwk20(j) = 0
    ic8 = ic9 + ic6
    ll = ic9 + m99
    jc2 = 0

    do l = ic8, ll

      if ( a(l) /= 0.0D+00 ) then
        jc2 = jc1 + 1
        jc1 = l
      end if

    end do

    if ( jc2 == 1 ) then
      jc3 = jc1 - ic9
      if ( iwk19(jc3) == 0 ) then
        if ( 0.0D+00 <= a(jc1) * wk16(jc3) ) then
          iwk19(jc3) = j
          iwk20(j) = jc3
          ic9 = ic9 + ic5
        end if
      end if
    end if

  end do
!
!  Form the inverse from IWK20.
!
40 continue

  loop3 = 1
  loop2 = 1

  if ( jc5 <= 0 ) then
    inc2 = 0
  end if
!
!  Old code:
!
! wk2(1:m2)
!
!  New code:
!
  wk2(1:(m+1)*(m+1)) = 0.0D+00

  ic8 = 1
  do i = 1, m99
    wk2(ic8) = 1.0D+00
    ic8 = ic8 + m99 + 1
  end do

  x(1:m99) = wk16(1:m99)

  do i = ic6, m99
    if ( iwk19(i) /= 0 ) then
      iwk19(i) = ispec
    end if
  end do

  infs = 1
  ic1 = 1
!
!  Form the inverse.
!
80    continue

  if ( iwk20(ic1) == 0 ) then
    go to 110
  else
    go to 250
  end if

90  continue
!
!  Reset the artificials.
!
  delta5 = 0.0D+00

  do i = ic6, m99

    if ( iwk19(i) == ispec ) then
      if ( delta5 < abs ( wk14(i) ) ) then
        ic3 = i
        delta5 = abs ( wk14(i) )
      end if
    end if

  end do

  if ( delta5 < delta1 ) then
    iwk20(ic1) =  0
  else
    iwk19(ic3) = ic1
    iwk20(ic1) = ic3
    go to 350
  end if

110  continue

  ic1 = ic1 + 1

  if ( ic1 <= n99 ) then
    go to 80
  end if

120  continue

  do i = 1, m99
    if ( iwk19(i) == ispec ) then
      iwk19(i) = 0
    end if
  end do

  loop1 = 1
  loop2 = 2
  loop3 = 2
!
!  Perform one iteration.
!
140  continue

  inc3 = 0
  inc4 = 0

  do i = ic6, m99

    if ( abs ( x(i) ) < delta2 ) then
      x(i) = 0.0D+00
    else
      if ( x(i) /= 0.0D+00 ) then
        if ( 0.0D+00 < x(i) ) then
          if ( iwk19(i) == 0 ) then
            inc3 = 1
          end if
        else
          inc4 = 1
          inc3 = 1
        end if

      end if
    end if

  end do
!
!  If infeasible, then invert again.
!
  if ( infs < inc3 ) then
    go to 40
  end if

  if ( inc3 < infs ) then
    infs = 0
    delta4 = 0.0D+00
  end if
!
!  Obtain prices.
!
160  continue


  ic8 = ic7

  do j = 1, m99
    wk13(j) = wk2(ic8)
    ic8 = ic8 + m99
  end do
!
!  Pricing.
!
  if ( infs /= 0 ) then

    wk13(1:m99) = wk13(1:m99) * delta4

    do i = ic6, m99

      ic8 = i

      if ( x(i) < 0.0D+00 ) then

        do j = 1, m99
          wk13(j) = wk13(j) + wk2(ic8)
          ic8 = ic8 + m99
        end do

      else if ( iwk19(i) == 0 ) then

        do j = 1, m99
!         write ( *, '(a,i8)' ) 'J=', j
!         write ( *, '(a,g14.6)' ) 'WK13(j)=', wk13(j)
!         write ( *, '(a,i8)' ) 'IC8=', ic8
!         write ( *, '(a,g14.6)' ) 'WK2(IC8)=', wk2(ic8)
          wk13(j) = wk13(j) - wk2(ic8)
          ic8 = ic8 + m99
        end do

      end if

    end do

  end if
!
!  Select entering column.
!
  ic1 = 0
  delta7 = delta3
  ic2 = 1

220   continue

  if ( iwk20(ic2) == 0 ) then
    go to 460
  else
    go to 240
  end if

230  continue

  if ( delta6 < delta7 ) then
    delta7 = delta6
    ic1 = ic2
  end if

240   continue

  ic2 = ic2 + 1
!
!  All costs are non-negative.
!
  if ( ic2 <= n99 ) then
    go to 220
  end if

  if ( ic1 <= 0 ) then
    k = 3 + infs
    go to 340
  end if
!
!  Multiply the column into the basis inverse.
!
250  continue

  wk14(1:m99) = 0.0D+00
  ic4 = ( ic1 - 1 ) * ic5
  ll = 0

  do i = 1, m99

    ic4 = ic4 + 1

    do j = 1, m99
      ll = ll + 1
      wk14(j) = wk14(j) + a(ic4) * wk2(ll)
    end do

  end do

  if ( loop2 == 1 ) then
    go to 90
  else if ( loop2 == 3 ) then
    return
  end if
!
!  Get the maximum value from the pivot column.
!
290   continue

  ic3 = 0
  delta8 = 0.0D+00
  jc3 = 0

  do i = ic6, m99

    if ( x(i) == 0.0D+00 ) then

      beta = abs ( wk14(i) )

      if ( delta1 < beta ) then

        if ( iwk19(i) == 0 ) then
          go to 300
        end if

        if ( jc3 == 0 ) then

          if ( 0.0D+00 < wk14(i) ) then

300        continue

            if ( jc3 /= 0 ) then

              if ( delta8 < beta ) then
                delta8 = beta
                ic3 = i
              end if

            else

              jc3 = 1
              delta8 = beta
              ic3 = i

            end if

          end if

        end if

      end if

    end if

  end do
!
!  Find the maximum pivot from the positive equations.
!
  if ( ic3 == 0 ) then

    delta8 = huge ( 1.0D+00 )

    do it = ic6, m99

      if ( delta1 < wk14(it) ) then

        if ( 0.0D+00 < x(it) ) then
          delta9 = x(it) / wk14(it)
          if ( delta9 < delta8 ) then
            delta8 = delta9
            ic3 = it
          else
            if ( delta9 == delta8 ) then
              if ( iwk19(it) == 0 ) then
                delta8 = delta9
                ic3 = it
              end if
            end if
          end if
        end if

      end if

    end do
!
!  Find pivot among negative equations.
!
     if ( inc4 /= 0 ) then

       delta7 = - delta1

       do i = ic6, m99

         if ( x(i) < 0.0D+00 ) then
           if ( wk14(i) < delta7 ) then
             if ( wk14(i) * delta8 <= x(i) ) then
               delta7 = wk14(i)
               ic3 = i
             end if
           end if
         end if

       end do

     end if
  end if
!
!  Test pivot.
!
  if ( ic3 <= 0 ) then

    k = 5

340 continue

    if ( delta4 /= 0.0D+00 ) then
      delta4 = 0.0D+00
      go to 160
    end if

    go to 400

  end if
!
!  Terminate if too many iterations.
!
  if ( iter < maxit ) then

350 continue

    beta = - wk14(ic3)
    wk14(ic3) = - 1.0D+00
    ll = 0
!
!  Transform inverse.
!
    do l = ic3, m2, m99
      if ( wk2(l) == 0.0D+00 ) then
        ll = ll + m99
      else
        delta9 = wk2(l) / beta
        wk2(l) = 0.0D+00
        do i = 1, m99
          ll = ll + 1
          wk2(ll) = wk2(ll) + delta9 * wk14(i)
        end do
      end if
    end do
!
!  Transform X.
!
    delta9 = x(ic3) / beta
    x(ic3) = 0.0D+00
    x(1:m99) = x(1:m99) + delta9 * wk14(1:m99)
!
!  Restore the pivot.
!
    wk14(ic3) = - beta

    if ( loop3 == 1 ) then
      go to 120
    end if

390 continue

    jc3 = iwk19(ic3)

    if ( 0 < jc3 ) then
      iwk20(jc3) = 0
    end if

    iwk20(ic1) = ic3
    iwk19(ic3) = ic1
    jc5 = 0
    iter = iter + 1
    inc2 = inc2 + 1

    if ( inc1 == inc2 ) then
      go to 40
    else
      go to 140
    end if

  end if

  k = 6
!
!  Error checking.
!
400 continue

  loop1 = 2
  wk14(1:m99) = - wk16(1:m99)

  do i = 1, m99

    jc4 = iwk19(i)

    if ( jc4 /= 0 ) then
      jc3 = ic5 * ( jc4 - 1 )
      wk14(1:m99) = wk14(1:m99) + x(i) * a(jc3+1:jc3+m99)
    end if

  end do
!
!  Set error flags.
!
  i = 1

440 continue

  ic2 = iwk19(i)

  if ( ic2 == 0 ) then

450   continue

    i = i + 1

    if ( i <= m99 ) then
      go to 440
    end if

    if ( jc5 == 0 ) then
      jc5 = 4
    end if

    if ( k /= 5 ) then
      return
    end if

    loop2 = 3
    go to 250

  end if
!
!  Price out a column.
!
460   continue

  ll = ( ic2 - 1 ) * ic5
  delta6 = 0.0D+00
  do ic8 = 1, m99
    ll = ll + 1
    if ( a(ll) /= 0.0D+00 ) then
      delta6 = delta6 + wk13(ic8) * a(ll)
    end if
  end do

  if ( loop1 == 1 ) then
    go to 230
  else
    go to 450
  end if

end
subroutine steiner ( n, m, inode, jnode, cost, ns, spoint, nsp, &
  istree, jstree, xlen )

!*****************************************************************************80
!
!! STEINER is a heuristic algorithm for the minimal Steiner tree problem.
!
!  Discussion:
!
!    Consider an undirected graph G whose edges have specified lengths.
!
!    Let S be a user-specified subset of the nodes of G, called Steiner
!    points.
!
!    The minimal Steiner tree problem is to find a subset of the edges
!    of G which forms a tree connecting all the nodes of S, and
!    whose edges have minimal total length.
!
!    If N is the number of nodes in G, and P the number of Steiner points,
!    and K the number of leaves in the optimal Steiner tree, then the
!    heuristic algorithm implemented here will find an approximate
!    solution whose total length is no more than
!
!      2 * ( 1 - 1/K )
!
!    times that of the optimal tree, and will do this in time
!    proportional to P * N * N.
!
!  Algorithm:
!
!    Step 1: Construct the complete undirected graph H from G and S,
!    in such a way that the set of nodes in H is equal to S, and
!    for every edge (U,V) in H, the length of (U,V) is set equal to
!    the distance of the shortest path between node U and node V in G.
!
!    Step 2: Find a minimum spanning tree T(H) of H.
!
!    Step 3: Replace each edge (U,V) in T(H) by a shortest path
!    between node U and node V in G.  The resulting graph R is a
!    subgraph of G.
!
!    Step 4: Find a minimum spanning tree T(R) of R.
!
!    Step 5: Delete edges in T(R), if necessary, so that all the
!    leaves in T(R) are elements of S.  The resulting tree is returned
!    as the approximate solution.
!
!  Remarks:
!
!    To estimate how close the heuristic solution comes to an optimal one,
!    let T(OPT) be the optimal Steiner tree, and K the number of leaves
!    in T(OPT).  Denote the total length of the edges of T(OPT) by
!    Z(OPT), and the total length of the edges of the approximate Steiner tree
!    computed by the heuristic algorithm as Z(H).
!
!    If an edge is added in parallel to every edge in T(OPT), then there
!    is an Euler circuit C in T(OPT).  The circuit C can be regarded
!    as composed of K simple paths, each connecting one leaf to another.
!    By deleting from C the longest such simple path, a walk W can be
!    constructed such that the total length of W is no more than
!
!      ( 1 - 1 / K ) * total length of C
!
!    which is equal to
!
!      2 * ( 1 - 1 / K ) * Z(OPT).
!
!    On the other hand, the total length of W is greater than the
!    total length of the edges of T(H), so
!
!      Z(H) <= 2 * ( 1 - 1 / K ) * Z(OPT)
!
!    as was asserted.
!
!  Modified:
!
!    11 September 1999
!
!  Author:
!
!    Hang Tong Lau
!
!  Reference:
!
!    Lawrence Kou, George Markowsky, Leonard Berman,
!    A fast algorithm for Steiner trees,
!    Research report RC 7390,
!    IBM Thomas J Watson Research Center,
!    Yorktown Heights, New York, 1978.
!
!    Hang Tong Lau,
!    Combinatorial Heuristic Algorithms in FORTRAN,
!    Springer Verlag, 1986,
!    ISBN: 3540171614,
!    LC: QA402.5.L37.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of nodes in the graph.
!
!    Input, integer ( kind = 4 ) M, the number of edges in the graph.
!
!    Input, integer ( kind = 4 ) INODE(M), JNODE(M), the edges of the graph,
!    described by pairs of nodes.
!
!    Input, real ( kind = 8 ) COST(M), the length of each edge in the graph.
!
!    Input, integer ( kind = 4 ) NS, the number of nodes which are
!    Steiner points.
!
!    Input, logical SPOINT(N), is TRUE for each node that is a Steiner point.
!
!    Output, integer ( kind = 4 ) NSP, the number of edges in the Steiner tree.
!
!    Output, integer ( kind = 4 ) ISTREE(N), JSTREE(N), the edges of the
!    Steiner tree, stored as pairs of nodes, in entries 1 through NSP.
!
!    Output, real ( kind = 8 ) XLEN, the total length of the edges in the
!    Steiner tree.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ns

  real ( kind = 8 ) cost(m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic1
  integer ( kind = 4 ) ic2
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) ij
  integer ( kind = 4 ) il
  integer ( kind = 4 ) in
  integer ( kind = 4 ) inode(m)
  integer ( kind = 4 ) istree(n)
  integer ( kind = 4 ) it
  integer ( kind = 4 ) iv1
  integer ( kind = 4 ) iv2
  integer ( kind = 4 ) iwk1(n)
  integer ( kind = 4 ) iwk2(n)
  integer ( kind = 4 ) iwk4(n)
  integer ( kind = 4 ) iwk6(n)
  integer ( kind = 4 ) iwk8(m)
  integer ( kind = 4 ) iwk9(m)
  integer ( kind = 4 ) iwk10(m)
  integer ( kind = 4 ) iwk11(m)
  integer ( kind = 4 ) iwk12(ns)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) jk
  integer ( kind = 4 ) jnode(m)
  integer ( kind = 4 ) jstree(n)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kk
  integer ( kind = 4 ) kk1
  integer ( kind = 4 ) l
  integer ( kind = 4 ) ll
  integer ( kind = 4 ) mx
  integer ( kind = 4 ) nsp
  integer ( kind = 4 ) nump
  logical spoint(n)
  real ( kind = 8 ) work14(m)
  real ( kind = 8 ) xlen
!
!  Identify the Steiner points.
!
  it = 0
  do i = 1, n
    if ( spoint(i) ) then
      it = it + 1
      iwk12(it) = i
    end if
  end do
!
!  Construct the complete graph for the Steiner points.
!
  il = 0
  ii = 2

  do i = 1, ns - 1

    do j = ii, ns

      if ( i /= j ) then

        il = il + 1
        iwk10(il) = iwk12(i)
        iwk11(il) = iwk12(j)

        call graph_arc_min_path ( n, m, inode, jnode, cost, iwk12(i), &
          iwk12(j), nump, iwk4, xlen )

        work14(il) = xlen

      endif

    end do

    ii = ii + 1

  end do

  ll = ( ns * ( ns - 1 ) ) / 2

  mx = 0
  do i = 1, it
    mx = max ( mx, iwk12(i) )
  end do
!
!  Construct the complete graph for the Steiner points.
!
!  Find a minimum spanning tree of the complete graph.
!
  call graph_arc_min_span_tree ( mx, ll, iwk10, iwk11, work14, &
    iwk1, iwk2, xlen )

  iwk10(1:m) = 0
!
!  Construct the subgraph by replacing each edge of the
!  minimum spanning tree by its corresponding shortest path.
!
  do i = 1, ns - 1

    ii = iwk1(i)
    jj = iwk2(i)

    call graph_arc_min_path ( n, m, inode, jnode, cost, ii, &
      jj, nump, iwk4, xlen )

    do ij = 1, nump - 1

      iv1 = iwk4(ij)
      iv2 = iwk4(ij+1)

      do jk = 1, m
        if ( ( inode(jk) == iv1 .and. jnode(jk) == iv2 ) .or. &
             ( inode(jk) == iv2 .and. jnode(jk) == iv1 ) ) then
           kk1 = jk
        end if
      end do

      iwk10(kk1) = - 1

    end do
  end do

  do i = 1, m
    if ( iwk10(i) == 0 ) then
      inode(i) = 0
      jnode(i) = 0
    end if
  end do
!
!  Find a minimum spanning tree of the subgraph.
!
  call graph_arc_min_span_tree ( n, m, inode, jnode, cost, &
    iwk1, iwk2, xlen )

  in = 0
  do i = 1, n
    if ( iwk1(i) /= 0 ) then
      in = in + 1
    end if
  end do
!
!  Construct a Steiner tree by deleting edges, if
!  such that all leaves are Steiner points.
!
  do j = 1, in

    iflag = 0

    do i = 1, in

      if ( iwk1(i) /= 0 .and. iwk2(i) /= 0 ) then

        l = iwk2(i)
        ic1 = 0

        do k = 1, in

          if ( iwk1(k) == l ) then
            ic1 = ic1 + 1
            iwk9(ic1) = k
            iwk8(ic1) = iwk2(k)
          endif

          if ( iwk2(k) == l ) then
            ic1 = ic1 + 1
            iwk9(ic1) = k
            iwk8(ic1) = iwk1(k)
          endif

        end do

        l = iwk1(i)
        ic2 = 0

        do k = 1, in

          if ( iwk1(k) == l ) then
            ic2 = ic2 + 1
            iwk9(ic2) = k
            iwk8(ic2) = iwk2(k)
          endif

          if ( iwk2(k) == l ) then
            ic2 = ic2 + 1
            iwk9(ic2) = k
            iwk8(ic2) = iwk1(k)
          endif

        end do

        ii = iwk1(i)
        jj = iwk2(i)

        if ( ( ic1 == 1 .and. .not. spoint(jj) ) .or. &
             ( ic2 == 1 .and. .not. spoint(ii) ) ) then

          do k = 1, m
            if ( inode(k) == iwk1(i) .and. jnode(k) == iwk2(i) ) then
              kk = k
              exit
            end if
          end do

          iwk1(i) = 0
          iwk2(i) = 0
          iflag = 1
          xlen = xlen - cost(kk)

        end if

      end if

    end do

    if ( iflag == 0 ) then
      exit
    end if

  end do
!
!  Store the solution.
!
  i = 0

  do j = 1, in
    if ( iwk1(j) /= 0 .and. iwk2(j) /= 0 ) then
      i = i + 1
      istree(i) = iwk1(j)
      jstree(i) = iwk2(j)
    endif
  end do

  nsp = i

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
  character ( len = 10 ) time
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
subroutine tsp ( nnode, dist, isol )

!*****************************************************************************80
!
!! TSP is a heuristic algorithm for the traveling salesman problem.
!
!  Discussion:
!
!    Let G be a complete graph with an associated distance matrix
!    DIST(I,J) on its edges.  The traveling salesman problem is to
!    start from a node in G, visit all the other nodes exactly once,
!    and return to the starting node in such a way that the total
!    traveled distance is a minimum.
!
!    In general, it is very hard to develop efficient algorithms
!    that yield good approximate solutions to this problem.  However,
!    in the special case where the distance matrix is symmetric:
!
!      DIST(I,J) = DIST(J,I)
!
!    and satisfies the triangle inequality:
!
!      DIST(I,K) <= DIST(I,J) + DIST(J,K),
!
!    solutions close to the optimum value can be found in a relatively
!    short time.  The algorithm here will find a circuit of no worse
!    than 3/2 the optimum length, in polynomial time.
!
!  Modified:
!
!    12 October 2010
!
!  Author:
!
!    Hang Tong Lau
!
!  Reference:
!
!    Nicos Christofides,
!    Worst-case analysis of a new heuristic for the traveling salesman
!    problem,
!    Management Science Research Report Number 388,
!    Carnegie-Mellon University, 1976.
!
!    Hang Tong Lau,
!    Combinatorial Heuristic Algorithms in FORTRAN,
!    Springer Verlag, 1986,
!    ISBN: 3540171614,
!    LC: QA402.5.L37.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, real ( kind = 8 ) DIST(NNODE,NNODE), the distance between
!    each pair of nodes.
!
!    Output, integer ( kind = 4 ) ISOL(NNODE), contains the nodes
!    in the order in which they should be visited.
!
  implicit none

  integer ( kind = 4 ) nnode

  real ( kind = 8 ) dist(nnode,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ict
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) inode((3*nnode)/2)
  integer ( kind = 4 ) isol(nnode)
  integer ( kind = 4 ) iwk10((3*nnode)/2)
  integer ( kind = 4 ) iwk11(nnode)
  integer ( kind = 4 ) iwk18(nnode)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) jnode((3*nnode)/2)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) match
  integer ( kind = 4 ) nedge
  real ( kind = 8 ) wk1( ( nnode * ( nnode - 1 ) ) / 2)
!
!  From the graph described by DIST, construct a minimum length spanning
!  tree, described by INODE, JNODE.
!
  call graph_dist_min_span_tree3 ( nnode, dist, inode, jnode )
!
!  Compute the degree (in the spanning tree) of each node.
!
  iwk10(1:nnode) = 0

  do i = 1, nnode - 1
    ii = inode(i)
    jj = jnode(i)
    iwk10(ii) = iwk10(ii) + 1
    iwk10(jj) = iwk10(jj) + 1
  end do
!
!  Construct a list of the nodes of odd degree.
!
  ict = 0
  do i = 1, nnode
    if ( mod ( iwk10(i), 2 ) /= 0 ) then
      ict = ict + 1
      iwk11(ict) = i
    end if
  end do
!
!  Determine a minimum weight perfect matching for the set of odd-degree nodes.
!
  k = 0
  do i = 2, ict
    do j = 1, i - 1
      k = k + 1
      wk1(k) = dist(iwk11(j),iwk11(i))
    end do
  end do

  call pmatch ( ict, wk1, iwk18 )
!
!  Store up the edges in the perfect matching.
!
  do i = 1, ict
    iwk18(i) = iwk11(iwk18(i))
  end do

  iwk10(1:nnode) = 0

  k = nnode - 1

  do i = 1, ict

    if ( iwk10(iwk11(i)) == 0 ) then
      iwk10(iwk11(i)) = 1
      iwk10(iwk18(i)) = 1
      k = k + 1
      inode(k) = iwk11(i)
      jnode(k) = iwk18(i)
    end if

  end do
!
!  Find an Euler circuit.
!
  nedge = nnode - 1 + ( ict / 2 )
  call graph_arc_euler_circ ( nnode, nedge, inode, jnode, iwk10 )
!
!  Form the Hamiltonian circuit.
!
  isol(1) = iwk10(1)
  j = 2
  isol(j) = iwk10(j)
  k = 2

  do

    match = 0

    do while ( match == 1 )

      j = j + 1

      match = 0

      do i = 1, j - 1
        if ( iwk10(j) == iwk10(i) ) then
          match = 1
          exit
        end if
      end do

    end do

    k = k + 1
    isol(k) = iwk10(j)

    if ( nnode <= k ) then
      exit
    end if

  end do

  return
end
