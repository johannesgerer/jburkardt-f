subroutine color2 ( n, m, level, middle, fail, parm1, paint, arcsec, arctop, &
  examin, arctyp )

!*****************************************************************************80
!
!! COLOR2 carries out a two-coloring for a planarity test.
!
!  Modified:
!
!    02 August 2000
!
!  Author:
!
!    Hang Tong Lau
!
!  Reference:
!
!    Hang Tong Lau,
!    Algorithms on Graphs,
!    Tab Books, 1989,
!    QA166 L38.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of nodes.
!
!    Input, integer ( kind = 4 ) M, the number of edges.
!
!    Input/output, integer ( kind = 4 ) LEVEL, ?
!
!    Input/output, logical MIDDLE, ?
!
!    Output, logical FAIL, is TRUE if the algorithm failed.
!
!    Input/output, integer ( kind = 4 ) PARM1(M), ?
!
!    Input/output, integer ( kind = 4 ) PAINT(M+2-N), ?
!
!    Input/output, integer ( kind = 4 ) ARCSEC(7*M+2-5*N), ?
!
!    Input, integer ( kind = 4 ) ARCTOP(7*M+2-5*N), ?
!
!    Input/output, logical EXAMIN(M+2-N), ?
!
!    Input, logical ARCTYP(7*M+2-5*N), ?
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) arcsec(7*m+2-5*n)
  integer ( kind = 4 ) arctop(7*m+2-5*n)
  logical arctyp(7*m+2-5*n)
  logical dum1
  logical dum2
  logical examin(m+2-n)
  logical fail
  integer ( kind = 4 ) level
  integer ( kind = 4 ) link
  logical middle
  integer ( kind = 4 ) paint(m+2-n)
  integer ( kind = 4 ) parm1(m)
  integer ( kind = 4 ) qnode
  integer ( kind = 4 ) tnode

  fail = .false.

  if ( middle ) then
    level = level - 1
    qnode = parm1(level)
  else
    qnode = parm1(level)
  end if

  do while ( arcsec(qnode) /= 0 )

    link = arcsec(qnode)
    tnode = arctop(link)
    arcsec(qnode) = arcsec(link)

    if ( paint(tnode) == 0 ) then

      if ( arctyp(link) ) then
        paint(tnode) = paint(qnode)
      else
        paint(tnode) = 3 - paint(qnode)
      end if

    else

      dum1 = ( paint(tnode) == paint(qnode) )
      dum2 = .not. arctyp(link)

      if ( dum1 .eqv. dum2 ) then
        fail = .true.
        return
      end if

    end if

    if ( examin(tnode) ) then
      examin(tnode) = .false.
      level = level + 1
      parm1(level) = tnode
      middle = .false.
      return
    end if

  end do

  middle = .true.

  return
end
subroutine decomp ( n, m, level, middle, initp, snode, pnum, nexte, pw18, &
  pw19, pw24, pw25, trail, descp, pindex, parm1, start, finish, first, &
  second, auxpf1, auxpf2, auxpf3, auxpf4, arcsec, arctop, arctyp )

!*****************************************************************************80
!
!! DECOMP does path decomposition for a planarity test.
!
!  Modified:
!
!    03 August 2000
!
!  Author:
!
!    Hang Tong Lau
!
!  Reference:
!
!    Hang Tong Lau,
!    Algorithms on Graphs,
!    Tab Books, 1989,
!    QA166 L38.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of nodes.
!
!    Input, integer ( kind = 4 ) M, the number of edges.
!
!    ?, integer LEVEL, ?
!
!    ?, logical MIDDLE, ?
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) arcsec(7*m+2-5*n)
  integer ( kind = 4 ) arctop(7*m+2-5*n)
  logical arctyp(7*m+2-5*n)
  integer ( kind = 4 ) auxpf1(m+m+2)
  integer ( kind = 4 ) auxpf2(m+m+2)
  integer ( kind = 4 ) auxpf3(m+m+m+3)
  integer ( kind = 4 ) auxpf4(m+m+m+3)
  integer ( kind = 4 ) descp(n)
  integer ( kind = 4 ) finish(m+2-n)
  integer ( kind = 4 ) first(n+m+m)
  logical ind
  integer ( kind = 4 ) initp
  integer ( kind = 4 ) level
  logical middle
  integer ( kind = 4 ) nexte
  integer ( kind = 4 ) node1
  integer ( kind = 4 ) node2
  integer ( kind = 4 ) parm1(m)
  integer ( kind = 4 ) pindex(2*n)
  integer ( kind = 4 ) pnum
  integer ( kind = 4 ) pw18
  integer ( kind = 4 ) pw19
  integer ( kind = 4 ) pw24
  integer ( kind = 4 ) pw25
  integer ( kind = 4 ) qnode
  integer ( kind = 4 ) qnode2
  integer ( kind = 4 ) second(n+m+m)
  integer ( kind = 4 ) snode
  integer ( kind = 4 ) start(m+2-n)
  integer ( kind = 4 ) tnode
  integer ( kind = 4 ) trail(n)

  if ( middle ) go to 20

  qnode = parm1(level)

  do while ( second(qnode) /= 0 )

    tnode = first(second(qnode))
    second(qnode) = second(second(qnode))

    if ( initp == 0 ) then
      initp = qnode
    end if

    if ( qnode < tnode ) then

      descp(tnode) = snode
      trail(tnode) = pnum
      level = level + 1
      parm1(level) = tnode
      middle = .false.
      return

20    continue

      tnode = parm1(level)
      level = level - 1
      qnode = parm1(level)
      snode = tnode - 1
      initp = 0

      do while ( qnode <= auxpf2(pw19+2) )
        pw19 = pw19 - 2
      end do

      do while ( qnode <= auxpf1(pw18+2) )
        pw18 = pw18 - 2
      end do

      do while ( qnode <= auxpf3(pw24+3) )
        pw24 = pw24 - 3
      end do

      do while ( qnode <= auxpf4(pw25+3) )
        pw25 = pw25 - 3
      end do

      ind = .false.
      qnode2 = qnode + qnode

      do while ( auxpf3(pw24+2) < pindex(qnode2-1) .and. &
        qnode < auxpf3(pw24+2) .and. &
        pindex(qnode2) < auxpf3(pw24+1) )

        ind = .true.

        node1 = pindex(qnode2)
        node2 = auxpf3(pw24+1)
        nexte = nexte + 1
        arcsec(nexte) = arcsec(node1)
        arcsec(node1) = nexte
        arctop(nexte) = node2

        node1 = auxpf3(pw24+1)
        node2 = pindex(qnode2)
        nexte = nexte + 1
        arcsec(nexte) = arcsec(node1)
        arcsec(node1) = nexte
        arctop(nexte) = node2

        arctyp(nexte-1) = .false.
        arctyp(nexte) = .false.
        pw24 = pw24 - 3

      end do

      if ( ind ) then
        pw24 = pw24 + 3
      end if

      pindex(qnode2-1) = 0
      pindex(qnode2) = 0

    else

      start(pnum+1) = initp
      finish(pnum+1) = tnode
      ind = .false.

      if ( auxpf1(pw18+2) /= 0 ) then
        pw19 = pw19 + 2
        auxpf2(pw19+1) = auxpf1(pw18+1)
        auxpf2(pw19+2) = auxpf1(pw18+2)
      end if

      if ( finish(auxpf1(pw18+1)+1) /= tnode ) then

        do while ( tnode < auxpf2(pw19+2) )

          node1 = pnum
          node2 = auxpf2(pw19+1)
          nexte = nexte + 1
          arcsec(nexte) = arcsec(node1)
          arcsec(node1) = nexte
          arctop(nexte) = node2

          node1 = auxpf2(pw19+1)
          node2 = pnum
          nexte = nexte + 1
          arcsec(nexte) = arcsec(node1)
          arcsec(node1) = nexte
          arctop(nexte) = node2

          arctyp(nexte-1) = .true.
          arctyp(nexte) = .true.
          ind = .true.
          pw19 = pw19 - 2

        end do

        if ( ind ) then
          pw19 = pw19 + 2
        end if

        ind = .false.

        do while ( tnode < auxpf3(pw24+3) .and. initp < auxpf3(pw24+2) )

          node1 = pnum
          node2 = auxpf3(pw24+1)
          nexte = nexte + 1
          arcsec(nexte) = arcsec(node1)
          arcsec(node1) = nexte
          arctop(nexte) = node2

          node1 = auxpf3(pw24+1)
          node2 = pnum
          nexte = nexte + 1
          arcsec(nexte) = arcsec(node1)
          arcsec(node1) = nexte
          arctop(nexte) = node2

          arctyp(nexte-1) = .false.
          arctyp(nexte) = .false.

         pw24 = pw24 - 3

        end do

        do while ( tnode < auxpf4(pw25+3) .and. initp < auxpf4(pw25+2) )
          pw25 = pw25 - 3
        end do

        if ( pindex(2*tnode-1) < initp ) then
          pindex(2*tnode-1) = initp
          pindex(2*tnode) = pnum
        end if

        pw24 = pw24 + 3
        auxpf3(pw24+1) = pnum
        auxpf3(pw24+2) = initp
        auxpf3(pw24+3) = tnode
!
!  Typo in the original!
!
        pw25 = pw25 + 3
        auxpf3(pw25+1) = pnum
        auxpf3(pw25+2) = initp
        auxpf3(pw25+3) = tnode

      else

        do while ( tnode < auxpf4(pw25+3) .and. initp < auxpf4(pw25+2) .and. &
                   auxpf4(pw25+2) <= descp(initp) )

          ind = .true.

          node1 = pnum
          node2 = auxpf4(pw25+1)
          nexte = nexte + 1
          arcsec(nexte) = arcsec(node1)
          arcsec(node1) = nexte
          arctop(nexte) = node2

          node1 = auxpf4(pw25+1)
          node2 = pnum
          nexte = nexte + 1
          arcsec(nexte) = arcsec(node1)
          arcsec(node1) = nexte
          arctop(nexte) = node2

          arctyp(nexte-1) = .false.
          arctyp(nexte) = .false.
          pw25 = pw25 - 3

        end do

        if ( ind ) then
          pw25 = pw25 + 3
        end if

      end if

      if ( qnode /= initp ) then
        pw18 = pw18 + 2
        auxpf1(pw18+1) = pnum
        auxpf1(pw18+2) = initp
      end if

      pnum = pnum + 1
      initp = 0

    end if

  end do

  middle = .true.

  return
end
subroutine dfs1 ( n, m, level, middle, snum, pw12, mark, small1, small2, &
  parm1, parm2, stacke, first, second )

!*****************************************************************************80
!
!! DFS1 carries out the first depth first search for a planarity test.
!
!  Modified:
!
!    04 August 2000
!
!  Author:
!
!    Hang Tong Lau
!
!  Reference:
!
!    Hang Tong Lau,
!    Algorithms on Graphs,
!    Tab Books, 1989,
!    QA166 L38.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of nodes.
!
!    Input, integer ( kind = 4 ) M, the number of edges.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) first(n+m+m)
  integer ( kind = 4 ) level
  integer ( kind = 4 ) mark(n)
  logical middle
  integer ( kind = 4 ) parm1(m)
  integer ( kind = 4 ) parm2(m)
  integer ( kind = 4 ) pnode
  integer ( kind = 4 ) pw12
  integer ( kind = 4 ) qnode
  integer ( kind = 4 ) second(n+m+m)
  integer ( kind = 4 ) small1(n)
  integer ( kind = 4 ) small2(n)
  integer ( kind = 4 ) snum
  integer ( kind = 4 ) stacke(m+m)
  integer ( kind = 4 ) tnode

  if ( middle ) then

    tnode = parm1(level)
    qnode = parm2(level)
    level = level - 1
    pnode = parm2(level)

    if ( small1(tnode) < small1(qnode) ) then
      small2(qnode) = min ( small2(tnode), small1(qnode) )
      small1(qnode) = small1(tnode)
    else if ( small1(tnode) == small1(qnode) ) then
      small2(qnode) = min ( small2(tnode), small2(qnode) )
    else
      small2(qnode) = min ( small1(tnode), small2(qnode) )
    end if

  else

    qnode = parm1(level)
    pnode = parm2(level)

  end if

  do while ( 0 < second(qnode) )

    tnode = first(second(qnode))
    second(qnode) = second(second(qnode))

    if ( mark(tnode) < mark(qnode) .and. tnode /= pnode ) then

      pw12 = pw12 + 2
      stacke(pw12-1) = qnode
      stacke(pw12) = tnode

      if ( mark(tnode) == 0 ) then

        snum = snum + 1
        mark(tnode) = snum
        level = level + 1
        parm1(level) = tnode
        parm2(level) = qnode
        middle = .false.
        return

      else

        if ( mark(tnode) < small1(qnode) ) then
          small2(qnode) = small1(qnode)
          small1(qnode) = mark(tnode)
        else

          if ( small1(qnode) < mark(tnode) ) then
            small2(qnode) = min ( small2(qnode), mark(tnode) )
          end if

        end if

      end if

    end if

  end do

  middle = .true.

  return
end
subroutine dfs2 ( n, m, level, middle, snum, pw12, mark, parm1, &
  stacke, first, second )

!*****************************************************************************80
!
!! DFS2 carries out the second depth-first search for a planarity test.
!
!  Modified:
!
!    03 August 2000
!
!  Author:
!
!    Hang Tong Lau
!
!  Reference:
!
!    Hang Tong Lau,
!    Algorithms on Graphs,
!    Tab Books, 1989,
!    QA166 L38.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of nodes.
!
!    Input, integer ( kind = 4 ) M, the number of edges.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) first(n+m+m)
  integer ( kind = 4 ) level
  integer ( kind = 4 ) mark(n)
  logical middle
  integer ( kind = 4 ) parm1(m)
  integer ( kind = 4 ) pw12
  integer ( kind = 4 ) qnode
  integer ( kind = 4 ) second(n+m+m)
  integer ( kind = 4 ) snum
  integer ( kind = 4 ) stacke(m+m)
  integer ( kind = 4 ) tnode

  if ( middle ) then
    tnode = parm1(level)
    level = level - 1
    qnode = parm1(level)
    pw12 = pw12 + 2
    stacke(pw12-1) = mark(qnode)
    stacke(pw12) = mark(tnode)
  else
    qnode = parm1(level)
  end if

  do while ( 0 < second(qnode) )

    tnode = first(second(qnode))
    second(qnode) = second(second(qnode))

    if ( mark(tnode) == 0 ) then
      snum = snum + 1
      mark(tnode) = snum
      level = level + 1
      parm1(level) = tnode
      middle = .false.
      return
    end if

    pw12 = pw12 + 2
    stacke(pw12-1) = mark(qnode)
    stacke(pw12) = mark(tnode)

  end do

  middle = .true.

  return
end
subroutine digraph_arc_euler ( nnode, nedge, inode, jnode, success, trail )

!*****************************************************************************80
!
!! DIGRAPH_ARC_EULER returns an Euler circuit in a digraph.
!
!  Discussion:
!
!    An Euler circuit of a digraph is a path which starts and ends at
!    the same node and uses each directed edge exactly once.  A digraph is
!    eulerian if it has an Euler circuit.  The problem is to decide whether
!    a given digraph is eulerian and to find an Euler circuit if the
!    answer is affirmative.
!
!    The digraph is assumed to be connected.
!
!  Method:
!
!    A digraph has an Euler circuit if and only if the number of incoming
!    edges is equal to the number of outgoing edges at each node.
!
!    This characterization gives a straightforward procedure to decide whether
!    a digraph is eulerian.  Furthermore, an Euler circuit in an eulerian
!    digraph G of NEDGE edges can be determined by the following method:
!
!      STEP 1: Choose any node U as the starting node, and traverse any edge
!      ( U, V ) incident to node U, and than traverse any unused edge incident
!      to node U.  Repeat this process of traversing unused edges until the
!      starting node U is reached.  Let P be the resulting walk consisting of
!      all used edges.  If all edges of G are in P, than stop.
!
!      STEP 2: Choose any unused edge ( X,  Y) in G such that X is
!      in P and Y is not in P.  Use node X as the starting node and
!      find another walk Q using all unused edges as in step 1.
!
!      STEP 3: Walk P and walk Q share a common node X, they can be merged
!      to form a walk R by starting at any node S of P and to traverse P
!      until node X is reached; than, detour and traverse all edges of Q
!      until node X is reached and continue to traverse the edges of P until
!      the starting node S is reached.  Set P = R.
!
!      STEP 4: Repeat steps 2 and 3 until all edges are used.
!
!    The running time of the algorithm is O ( NEDGE ).
!
!  Modified:
!
!    25 July 2000
!
!  Author:
!
!    Hang Tong Lau
!
!  Reference:
!
!    Hang Tong Lau,
!    Algorithms on Graphs,
!    Tab Books, 1989,
!    QA166 L38.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) NEDGE, the number of edges.
!
!    Input, integer ( kind = 4 ) INODE(NEDGE), JNODE(NEDGE); the I-th edge
!    starts at node INODE(I) and ends at node JNODE(I).
!
!    Output, logical SUCCESS, is TRUE if an Euler circuit was found,
!    and FALSE otherwise.
!
!    Output, integer ( kind = 4 ) TRAIL(NEDGE).  TRAIL(I) is the edge number
!    of the I-th edge in the Euler circuit.
!
  implicit none

  integer ( kind = 4 ) nedge

  logical candid(nedge)
  integer ( kind = 4 ) endnod(nedge)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) inode(nedge)
  integer ( kind = 4 ) istak
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jnode(nedge)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) len
  integer ( kind = 4 ) lensol
  integer ( kind = 4 ) lenstk
  integer ( kind = 4 ) nnode
  integer ( kind = 4 ) stack(2*nedge)
  logical success
  integer ( kind = 4 ) trail(nedge)
!
!  Check if the digraph is eulerian.
!
  trail(1:nedge) = 0
  endnod(1:nedge) = 0

  do i = 1, nedge
    j = inode(i)
    trail(j) = trail(j) + 1
    j = jnode(i)
    endnod(j) = endnod(j) + 1
  end do

  do i = 1, nnode
    if ( trail(i) /= endnod(i) ) then
      success = .false.
      return
    end if
  end do
!
!  The digraph is eulerian; find an Euler circuit.
!
  success = .true.
  lensol = 1
  lenstk = 0
!
!  Find the next edge.
!
  do

    if ( lensol == 1 ) then

      endnod(1) = inode(1)
      stack(1) = 1
      stack(2) = 1
      lenstk = 2

    else

      l = lensol - 1

      if ( lensol /= 2 ) then
        endnod(l) = inode(trail(l)) + jnode(trail(l)) - endnod(l-1)
      end if

      k = endnod(l)

      do i = 1, nedge
        candid(i) = ( k == jnode(i) )
      end do

      do i = 1, l
        candid(trail(i)) = .false.
      end do

      len = lenstk

      do i = 1, nedge

        if ( candid(i) ) then
          len = len + 1
          stack(len) = i
        end if

      end do

      stack(len+1) = len - lenstk
      lenstk = len + 1

    end if

    do

      istak = stack(lenstk)
      lenstk = lenstk - 1

      if ( istak /= 0 ) then
        exit
      end if

      lensol = lensol - 1

      if ( lensol == 0 ) then
        call i4vec_reverse ( nedge, trail )
        return
      end if

    end do

    trail(lensol) = stack(lenstk)
    stack(lenstk) = istak - 1

    if ( lensol == nedge ) then
      exit
    end if

    lensol = lensol + 1

  end do

  call i4vec_reverse ( nedge, trail )

  return
end
subroutine digraph_arc_find_path ( nnode, nedge, i1, j1, fwdarc, arcfir, &
  mark, pexist )

!*****************************************************************************80
!
!! DIGRAPH_ARC_FIND_PATH determines if a path exists from 11 to J1.
!
!  Discussion:
!
!    Yen's algorithm is used.
!
!    This routine is a utility routine used by DIGRAPH_ARC_MINEQV.
!
!  Modified:
!
!    26 July 2000
!
!  Author:
!
!    Hang Tong Lau
!
!  Reference:
!
!    Hang Tong Lau,
!    Algorithms on Graphs,
!    Tab Books, 1989, pages 47-59,
!    QA166 L38.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) NEDGE, the number of edges.
!
!    Input, integer ( kind = 4 ) I1, J1, the nodes defining the ends of the
!    desired path.
!
!    Input, integer ( kind = 4 ) FWDARC(NEDGE); FWDARC(I) is the ending node
!    of the I-th edge in the forward star representation of the digraph.
!
!    Input, integer ( kind = 4 ) ARCFIR(NNODE+1); ARCFIR(I) is the number of
!    the first edge starting at node I in the forward star representation of
!    the digraph.
!
!    Workspace, logical MARK(NEDGE), a working copy of the array ARCLIS.
!
!    Output, logical PEXIST, is TRUE if a path exists from node I1 to node J1.
!
  implicit none

  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) arcfir(nnode+1)
  integer ( kind = 4 ) fwdarc(nedge)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) index1
  integer ( kind = 4 ) index2
  integer ( kind = 4 ) index3
  integer ( kind = 4 ) iup
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  logical join
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kedge
  integer ( kind = 4 ) low
  integer ( kind = 4 ) next(nnode)
  logical mark(nedge)
  logical nopath(nnode)
  logical pexist
!
!  Initialize.
!
  call i4vec_indicator ( nnode, next )

  nopath(1:nnode) = .FALSE.

  if ( i1 < 1 .or. nnode < i1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DIGRAPH_ARC_FIND_PATH - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of I1 = ', i1
    stop
  end if

  nopath(i1) = .TRUE.
  next(i1) = nnode
  index1 = i1
  index2 = nnode - 1
!
!  Compute the shortest distance labels.
!
  i = 1

20  continue

  j = next(i)
  join =.FALSE.
  low = arcfir(index1)
  iup = arcfir(index1+1) - 1

  do k = low, iup
    if ( fwdarc(k) == j ) then
      join =.TRUE.
      kedge = k
      exit
    end if
  end do

  if ( join ) then
    if ( mark(kedge) ) then
      nopath(j) = .TRUE.
    end if
  end if

  if ( .not. nopath(j) ) then

    i = i + 1

    if ( index2 < i ) then
      pexist = .FALSE.
      return
    end if

    go to 20

  end if

  index3 = i + 1

  if ( index3 <= index2 ) then

    do i2 = index3, index2

      j2 = next(i2)
      join = .FALSE.
      low = arcfir(index1)
      iup = arcfir(index1+1) - 1

      do k = low, iup

        if ( fwdarc(k) == j2 ) then
          join = .TRUE.
          kedge = k
          if ( mark(kedge) ) then
            nopath(j2) = .TRUE.
          end if
          exit
        end if

      end do

    end do

  end if
!
!  Check whether an alternative path exists.
!
  if ( nopath(j1) ) then
    pexist = .TRUE.
    return
  end if

  next(i) = next(index2)
  index1 = j
  index2 = index2 - 1

  if ( 1 < index2 ) then
    go to 20
  end if

  join = .FALSE.
  low = arcfir(index1)
  iup = arcfir(index1+1) - 1

  do k = low, iup
    if ( fwdarc(k) == j1 ) then
      join = .TRUE.
      kedge = k
      exit
    end if
  end do

  pexist = .FALSE.

  if ( join ) then
    if ( mark(kedge) ) then
      pexist = .TRUE.
    end if
  end if

  return
end
subroutine digraph_arc_get_path ( nnode, nedge, k, kpaths, maxque, inode, &
  jnode, pathlen, arcdir, qufirp, qunxtp, crosar, nextrd, nump, length, arcnod )

!*****************************************************************************80
!
!! DIGRAPH_ARC_GET_PATH retrieves the edges of the K-th shortest path.
!
!  Discussion:
!
!    This program is used after DIGRAPH_ARC_KSHORT2.
!
!  Modified:
!
!    23 July 2000
!
!  Author:
!
!    Hang Tong Lau
!
!  Reference:
!
!    Hang Tong Lau,
!    Algorithms on Graphs,
!    Tab Books, 1989,
!    QA166 L38.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) NEDGE, the number of edges.
!
!    Input, integer ( kind = 4 ) K, the number of the path to be generated.
!
!    Input, integer ( kind = 4 ) KPATHS, the number of shortest paths that were
!    generated.
!
!    Input, integer ( kind = 4 ) MAXQUE, the size of the auxilliary storage.
!
!    Input, integer ( kind = 4 ) INODE(NEDGE), JNODE(NEDGE); the I-th edge
!    is directed from INODE(I) to JNODE(I).
!
!    Input, integer ( kind = 4 ) PTHLEN(KPATHS+3); the shortest path lengths are stored
!    in PTHLEN(2), PTHLEN(3), ..., PTHLEN(KPATHS+1).
!
!    Input, integer ( kind = 4 ) ARCDIR(NNODE); array of shortest path tree, as computed
!    by KSHOT2.
!
!    Input, integer ( kind = 4 ) QUFIRP(KPATHS+3); description of the first path in
!    a queue, as computed by KSHOT2.
!
!    Input, integer ( kind = 4 ) QUNXTP(KPATHS+3); next entry in a queue, as computed
!    by KSHOT2.
!
!    Input, integer ( kind = 4 ) CROSAR(MAXQUE); cross edge of the path, as computed
!    by KSHOT2.
!
!    Input, integer ( kind = 4 ) NEXTRD(MAXQUE); next entry of the pool storage, as
!    computed by KSHOT2.
!
!    Output, integer ( kind = 4 ) NUMP, the number of edges in the path.
!
!    Output, integer ( kind = 4 ) LENGTH, the length of the path.
!
!    Output, integer ( kind = 4 ) ARCNOD(NNODE); the edges of the K-th shortest path
!    are stored in ARCNOD(1) through ARCNOD(NUMP).
!
  implicit none

  integer ( kind = 4 ) kpaths
  integer ( kind = 4 ) maxque
  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) arcdir(nnode)
  integer ( kind = 4 ) arcnod(nnode)
  integer ( kind = 4 ) crosar(maxque)
  integer ( kind = 4 ) index1
  integer ( kind = 4 ) index2
  integer ( kind = 4 ) index3
  integer ( kind = 4 ) inode(nedge)
  integer ( kind = 4 ) isub1
  integer ( kind = 4 ) isub2
  integer ( kind = 4 ) isub3
  integer ( kind = 4 ) jnode(nedge)
  integer ( kind = 4 ) jsub
  integer ( kind = 4 ) k
  integer ( kind = 4 ) length
  integer ( kind = 4 ) nextrd(maxque)
  integer ( kind = 4 ) number
  integer ( kind = 4 ) nump
  integer ( kind = 4 ) pathlen(kpaths+3)
  integer ( kind = 4 ) qufirp(kpaths+3)
  integer ( kind = 4 ) qunxtp(kpaths+3)

  number = 0

  if ( k <= 0 .or. qunxtp(1) < k ) then
    nump = number
    return
  end if

  index2 = pathlen(1)
  length = pathlen(k+1)
  isub3 = qufirp(k+1)

  do while ( isub3 /= 0 )

    jsub = abs ( isub3 ) + 1
    index3 = index2

    if ( 0 < crosar(jsub) ) then
      index1 = inode(crosar(jsub))
      index2 = jnode(crosar(jsub))
    else
      index1 = jnode(-crosar(jsub))
      index2 = inode(-crosar(jsub))
    end if

    if ( index2 /= index3 ) then
!
!  Store the edges.
!
      isub2 = nnode
      arcnod(isub2) = crosar(jsub)

      do while ( index1 /= index3 )

        isub1 = arcdir(index1)
        isub2 = isub2 - 1

        if ( 0 < isub2 ) then
          arcnod(isub2) = isub1
        else
          nump = isub1
        end if

        if ( 0 < isub1 ) then
          index1 = inode(isub1)
        else
          index1 = jnode(-isub1)
        end if

      end do

      do while ( isub2 <= nnode )
        number = number + 1
        arcnod(number) = arcnod(isub2)
        isub2 = isub2 + 1
      end do

    end if

    isub3 = nextrd(jsub)

  end do

  nump = number

  return
end
subroutine digraph_arc_hamcyc ( nnode, nedge, inode, jnode, success, hcycle )

!*****************************************************************************80
!
!! DIGRAPH_ARC_HAMCYC finds a Hamiltonian circuit in a digraph.
!
!  Discussion:
!
!    A hamiltonian cycle in a digraph G is a cycle containing every
!    node of G, and a digraph is hamiltonian if it has a hamiltonian cycle.
!    The problem is to decide whether a given digraph is hamiltonian and
!    to find a hamiltonian cycle if the answer is affirmative.
!
!  Method:
!
!    A hamiltonian cycle in a given digraph G of NNODE nodes will be found by
!    an exhaustive search.  Start with a single node, say node 1, as the
!    partially constructed cycle.  The cycle is grown by a backtracking
!    procedure until a hamiltonian cycle is formed.  More precisely, let
!
!      H(1), H(2), ..., H(K-1)
!
!    be the partially constructed cycle at the K-th stage.  Then, the
!    set of candidates, for the next element H(K) is the set of all
!    nodes U in G such that
!
!      if K = 1, than U = 1;
!      if 1 < K, then ( H(K-1), U ) is an edge in G, and U is
!         distinct from H(1), H(2), ..., H(K-1); furthermore,
!      Furthermore, if K = NNODE, then ( U, H(1) ) is an edge in G and
!      U < H(2).
!
!    The search is complete when K = NNODE.  If the set of candidates
!    is empty at any stage, then the digraph is not hamiltonian.
!
!  Modified:
!
!    23 July 2000
!
!  Author:
!
!    Hang Tong Lau
!
!  Reference:
!
!    Hang Tong Lau,
!    Algorithms on Graphs,
!    Tab Books, 1989,
!    QA166 L38.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) NEDGE, the number of edges.
!
!    Input, integer ( kind = 4 ) INODE(NEDGE), JNODE(NEDGE); the I-th edge of
!    the digraph extends from node INODE(I) to JNODE(I).
!
!    Output, logical SUCCESS, is TRUE if a hamiltonian cycle was found,
!    and FALSE otherwise.
!
!    Output, integer ( kind = 4 ) HCYCLE(NNODE); the nodes of the hamiltonian cycle are
!    stored in order in HCYCLE.
!
!    Workspace, integer FWDARC(NEDGE); FWDARC(I) is the ending node of the
!    I-th edge in the forward star representation of the digraph.
!
!    Workspace, integer ARCFIR(NNODE+1); ARCFIR(I) is the number of the first
!    edge starting at node I in the forward star representation of the digraph.
!
!    Workspace, integer STACK(NEDGE).
!
!    Workspace, logical CONNECT(NNODE).
!
  implicit none

  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) arcfir(nnode+1)
  logical connect(nnode)
  integer ( kind = 4 ) fwdarc(nedge)
  integer ( kind = 4 ) hcycle(nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) inode(nedge)
  integer ( kind = 4 ) istak
  integer ( kind = 4 ) iup
  integer ( kind = 4 ) jnode(nedge)
  logical join
  integer ( kind = 4 ) k
  integer ( kind = 4 ) len
  integer ( kind = 4 ) len1
  integer ( kind = 4 ) len2
  integer ( kind = 4 ) lensol
  integer ( kind = 4 ) lenstk
  integer ( kind = 4 ) low
  integer ( kind = 4 ) stack(nedge)
  logical success
!
!  Set up the forward star representation of the digraph.
!
  call digraph_arc_to_star ( nnode, nedge, inode, jnode, arcfir, fwdarc )
!
!  Initialize.
!
  lensol = 1
  lenstk = 0
  hcycle(1:nnode) = 0
!
!  Find the next node.
!
30    continue
!
!  Start with
!    STACK(1) = node 1,
!    STACK(2) = level 1.
!
  if ( lensol == 1 ) then

    stack(1) = 1
    stack(2) = 1
    lenstk = 2

  else

    len1 = lensol - 1
    len2 = hcycle(len1)
!
!  Seek a connection from the last node on the tentative list to any node I.
!
    do i = 1, nnode

      connect(i) = .false.
      low = arcfir(len2)
      iup = arcfir(len2+1) - 1

      do k = low, iup
        if ( fwdarc(k) == i ) then
          connect(i) = .true.
        end if
      end do

    end do
!
!  Disregard connections that take you to nodes already on the tentative list.
!
    do i = 1, len1
      connect ( hcycle(i) ) = .false.
    end do

    len = lenstk

    if ( lensol /= nnode ) then

      do i = 1, nnode
        if ( connect(i) ) then
          len = len + 1
          stack(len) = i
        end if
      end do

      stack(len+1) = len - lenstk
      lenstk = len + 1

    else

      do i = 1, nnode

        if ( connect(i) ) then

          join = .false.
          low = arcfir(i)
          iup = arcfir(i+1) - 1

          do k = low, iup

            if ( fwdarc(k) == 1 ) then
              join = .true.
              exit
            end if

          end do

          if ( join ) then
            lenstk = lenstk + 2
            stack(lenstk-1) = i
            stack(lenstk) = 1
          else
            stack(len+1) = len - lenstk
            lenstk = len + 1
          end if

          go to 110

        end if

      end do

      stack(len+1) = len - lenstk
      lenstk = len + 1

    end if

  end if
!
!  Search further.
!
110   continue

  istak = stack(lenstk)
  lenstk = lenstk - 1

  if ( istak == 0 ) then

    lensol = lensol - 1

    if ( lensol == 0 ) then
      success = .false.
      return
    end if

    go to 110

  else

    hcycle(lensol) = stack(lenstk)
    stack(lenstk) = istak - 1

    if ( lensol == nnode ) then
      success = .true.
      return
    end if

    lensol = lensol + 1
    go to 30

  end if

end
subroutine digraph_arc_kshort1 ( nnode, nedge, k, isorce, iter, inode, jnode, &
  arclen, dist )

!*****************************************************************************80
!
!! DIGRAPH_ARC_KSHORT1 finds the K shortest paths in a digraph.
!
!  Discussion:
!
!    Consider a digraph of NNODE nodes with lengths on edges.
!    It is assumed that all cycles in the digraph have strictly positive
!    lengths.  Given an integer 1 < K and a specified node S,
!    the problem is to find K shortest distinct path lengths from
!    node S to every other node in the digraph.  Furthermore, list
!    the K shortest paths from node S to every other node.  As
!    defined here, the K shortest paths are allowed to contain
!    embedded cycles, that is, some, nodes may be repeated in a path.
!
!  Method:
!
!    The K shortest path problem in a digraph of NNODE nodes with
!    given edge lengths D(I,J) will be solved by a label-correcting
!    method in which an initial guess is given to the K shortest
!    path lengths.  Then the tentative K shortest path lengths will
!    be improved successavely.
!
!    A sequence of alternating forward and backward iterations will be
!    employed.  During the forward iteration, the nodes are examined in
!    the order 1, 2, .., N, and only edges (I,J) with I < J are processed.
!    During the backward iteration, the nodes are examined in the order
!    NNODE, NNODE - 1, . . . , 1, and only edges (I,J) with J < I
!    are processed.
!
!    The alternating forward and backward iterations are continued until the
!    node labels at two consecutive iterations coincide, in which
!    case no improvement is possible.
!
!    The procedure of finding the K shortest path lengths
!    from node I to all other nodes cim be outlined as follows:
!
!    STEP 1. Initialize the required K shortest path lengths, a
!    K-vector
!
!      S(I) = ( S(I,1), S(I,2), ..., S(I,K) ),
!
!    from node 1 to node I by setting
!
!      S(1) = ( 0, Infinity, ..., Infinity ),
!
!    and
!
!      S(I) = ( Infinity, Infinity, ..., Infinity), for 1 < I.
!
!    STEP 2. For J = 2 to N, do the following: For every node I
!    adjacent to node J, where I < J, if
!
!      { S(I,P) + D(I,J): P = 1, 2, ..., K }
!
!    gives a smaller path length than any one of the tentative K
!    shortest path lengths in S(J), then the current K vector S(J)
!    is updated by inclusion of this smaller path length.
!
!    STEP 3. If none of the node labels S(J) has changed in Step 2,
!    then stop.
!
!    STEP 4. For J = NNODE - 1 to 1, do the following.  For every node
!    I adjacent to node J, where J < I, process the edge (I,J) the
!    same way as in Step 2.
!
!    STEP 5. If none of the node labels S(J) has changed in Step 4,
!    then stop; otherwise, return to Step 2.
!
!    After the K shortest distinct path lengths are determined by the
!    above procedure, the paths themselves can be reconstructed from the
!    path length information.  In essence, if a J-th shortest path P of
!    length T from node U to node V passes through node W, then the
!    subpath of P extending from node U to node W is an I-th shortest path,
!    for some I,  1 <= I <= J.
!
!    This can be used to determine the penultimate node W on a J-th shortest
!    path of known length T from node U to node V.  This backtracking method
!    can be applied repeatedly to produce all the K shortest paths from
!    node U to node V.
!
!    It should be noted that several paths can have the same
!    path lengths.  It is therefore possible that more than K paths
!    may be generated from the K distinct shortest path lengths.
!
!    The whole algorithm will take at most O ( N**2 * K**2 * log(N) )
!    operations.
!
!  Modified:
!
!    23 July 2000
!
!  Author:
!
!    Hang Tong Lau
!
!  Reference:
!
!    Hang Tong Lau,
!    Algorithms on Graphs,
!    Tab Books, 1989,
!    QA166 L38.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) NEDGE, the number of edges.
!
!    Input, integer ( kind = 4 ) K, number of shortest paths desired.
!
!    Input, integer ( kind = 4 ) ISORCE, the starting point for the shortest paths.
!
!    Input/output, integer ( kind = 4 ) ITER; if ITER is 0 on input, it is
!    reset to 100.  Otherwise, the input value is used to determine the maximum
!    number of iterations to take.  On output, ITER is the actual number of
!    iterations taken.
!
!    Input, integer ( kind = 4 ) INODE(NEDGE), JNODE(NEDGE).  The I-th edge is directed from
!    INODE(I) to JNODE(I).  It is assumed that the array JNODE is already
!    sorted in nondecreasing order.
!
!    Input, real ( kind = 8 ) ARCLEN(NEDGE); ARCLEN(I) is the length of the
!    I-th edge.
!
!    Output, real ( kind = 8 ) DIST(NNODE,K); DIST(I,J) is the length of the
!    J-th shortest path from ISORCE to node I.  If node I is not reachable
!    from ISORCE in the J-th shortest path, then DIST(I,J) will be set equal
!    to a number greater than the sum of all edge lengths in the input digraph.
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) nnode
  integer ( kind = 4 ) nij
  integer ( kind = 4 ) nji

  real ( kind = 8 ) arclen(nedge)
  logical best
  real ( kind = 8 ) dist(nnode,k)
  real ( kind = 8 ) dist1
  real ( kind = 8 ) dist2
  real ( kind = 8 ) dist3
  integer ( kind = 4 ) i
  real ( kind = 8 ) iaux(k)
  integer ( kind = 4 ) ids1
  integer ( kind = 4 ) ids2
  integer ( kind = 4 ) iedge1(nnode)
  integer ( kind = 4 ) iedge2(nnode)
  real ( kind = 8 ) ilen1(nedge)
  real ( kind = 8 ) ilen2(nedge)
  integer ( kind = 4 ) index
  integer ( kind = 4 ) index1
  integer ( kind = 4 ) index2
  integer ( kind = 4 ) init
  integer ( kind = 4 ) inod1(nedge)
  integer ( kind = 4 ) inod2(nedge)
  integer ( kind = 4 ) inode(nedge)
  integer ( kind = 4 ) isorce
  integer ( kind = 4 ) isub1
  integer ( kind = 4 ) isub2
  integer ( kind = 4 ) iter
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) jnode(nedge)
  real ( kind = 8 ) large
  real ( kind = 8 ) len
  integer ( kind = 4 ) loop1
  integer ( kind = 4 ) loop2
  real ( kind = 8 ) max
  integer ( kind = 4 ) maxitr
  integer ( kind = 4 ) nij1
  integer ( kind = 4 ) nji1
  integer ( kind = 4 ) nodeu
  integer ( kind = 4 ) nodev

  maxitr = iter

  if ( iter <= 0 ) then
    maxitr = 100
  end if

  nij = 0
  nji = 0
  do i = 1, nedge
    if ( inode(i) < jnode(i) ) then
      nij = nij + 1
    else
      nji = nji + 1
    end if
  end do

  index = 0
  init = 0
  isub1 = 0
  isub2 = 0

  large = 1.0D+00 + sum ( arclen(1:nedge) )

  nij1 = 0
  nji1 = 0

  do i = 1, nedge

    nodev = inode(i)
    nodeu = jnode(i)
    len = arclen(i)

    if ( nodeu /= index ) then

      j1 = index + 1
      j2 = nodeu - 1

      iedge1(j1:j2) = 0
      iedge2(j1:j2) = 0

      if ( init /= 0 ) then
        iedge1(index) = isub1
        iedge2(index) = isub2
      end if

      isub1 = 0
      isub2 = 0
      index = nodeu

    end if

    init = init + 1

    if ( nodev <= nodeu ) then
      nji1 = nji1 + 1
      inod2(nji1) = nodev
      ilen2(nji1) = len
      isub2 = isub2 + 1
    else
      nij1 = nij1 + 1
      inod1(nij1) = nodev
      ilen1(nij1) = len
      isub1 = isub1 + 1
    end if

  end do

  iedge1(index) = isub1
  iedge2(index) = isub2

  dist(1:nnode,1:k) = large
  dist(isorce,1) = 0.0D+00

  iter = 1

40    continue

  ids2 = nij1
  best = .true.
  i = nnode - 1

50    continue

  if ( 0 < i ) then

    if ( iedge1(i) /= 0 ) then

      ids1 = ids2 - iedge1(i) + 1
!
!  Matrix multiplication with DIST using the
!  lower triangular part of the edge distance matrix.
!
      iaux(1:k) = dist(i,1:k)

      max = iaux(k)

      do loop1 = ids1, ids2

        index1 = inod1(loop1)
        dist3 = ilen1(loop1)

        do loop2 = 1, k

          dist1 = dist(index1,loop2)

          if ( large <= dist1 ) then
            go to 100
          end if

          dist2 = dist1 + dist3

          if ( max <= dist2 ) then
            go to 100
          end if

          j = k

70            continue

          if ( 2 <= j ) then

            if ( dist2 < iaux(j-1) ) then
              j = j - 1
              go to 70
            end if

            if ( dist2 == iaux(j-1) ) then
              go to 90
            end if

          else

            j = 1

          end if

          index2 = k

          do

            if ( index2 <= j ) then
              exit
            end if
            iaux(index2) = iaux(index2-1)
            index2 = index2 - 1

          end do

          iaux(j) = dist2
          best = .false.
          max = iaux(k)

90        continue

        end do

100     continue

      end do

      if ( .not. best ) then
        dist(i,1:k) = iaux(1:k)
      end if

      ids2 = ids1 - 1

    end if

    i = i - 1
    go to 50

  end if
!
!  If no change between previous iteration and this one, return.
!
  if ( iter /= 1 ) then

    if ( best ) then
      return
    end if

  end if
!
!  Begin the next iteration.
!
  iter = iter + 1
  ids1 = 1
  best = .true.

  do i = 2, nnode

    if ( iedge2(i) /= 0 ) then

      ids2 = ids1 + iedge2(i) - 1
!
!  Matrix multiplication with DIST using the
!  upper triangular part of the edge distance matrix.
!
      iaux(1:k) = dist(i,1:k)

      max = iaux(k)

      do loop1 = ids1, ids2

        index1 = inod2(loop1)
        dist3 = ilen2(loop1)

        do loop2 = 1, k

          dist1 = dist(index1,loop2)

          if ( large <= dist1 ) then
            go to 160
          end if

          dist2 = dist1 + dist3

          if ( max <= dist2 ) then
            go to 160
          end if

          j = k

130       continue

          if ( 2 <= j ) then

            if ( dist2 < iaux(j-1) ) then
              j = j - 1
              go to 130
            else if ( dist2 == iaux(j-1) ) then
              go to 150
            end if

          else

            j = 1

          end if

          index2 = k

          do

            if ( index2 <= j ) then
              exit
            end if

            iaux(index2) = iaux(index2-1)
            index2 = index2 - 1

          end do

          iaux(j) = dist2
          best = .false.
          max = iaux(k)

150       continue

        end do

160     continue

      end do

      if ( .not. best ) then
        dist(i,1:k) = iaux(1:k)
      end if

      ids1 = ids2 + 1

    end if

  end do

  if ( .not. best ) then

    iter = iter + 1

    if ( iter < maxitr ) then
      go to 40
    end if

  end if

  return
end
subroutine digraph_arc_kshort2 ( nnode, nedge, inode, jnode, arclen, kpaths, &
  maxque, isorce, isink, iflag, ipaths, pathlen, arcdir, trdist, arcnod, &
  arcfwd, arcbwd, auxstg, auxdis, auxtre, auxlnk, nxtfwd, nxtbwd, qufirp, &
  qunxtp, crosar, nextrd )

!*****************************************************************************80
!
!! DIGRAPH_ARC_KSHORT2 finds the KPATHS shortest distinct path lengths without repeat nodes.
!
!  Discussion:
!
!    Consider a digraph of NNODE nodes with lengths on edges and which
!    has no cycle of negative length.  Given an integer 0 < K, a specified
!    source node S and a specified sink node T, the problem is to find the
!    K shortest paths from S to T such that each of the K shortest paths
!    does not have any embedded cycle.  That is, each path does not contain
!    any repeated nodes.
!
!    Warning:  This code uses nonstandard FORTRAN, in jumping to
!    statements 260 and 270, which are inside blocks.
!
!  Method:
!
!    Let P(I) be the I-th shortest path from the source S to the sink T.
!    Among the shortest paths P(1), P(2), ..., P(I-1) that have the same
!    initial subpaths from node S to the I-th node, let Q(I,J) be the shortest
!    of the paths that coincide with P(J-1) from node S up to the I-th node
!    and then deviate to a node that is different, from any (I+1)-th node
!    in other paths.
!
!    The procedure of finding the K shortest paths from S to T follows:
!
!    STEP 1: Find P(1), a shortest path from S to T.  If there is only
!    one shortest path, than enter it into a list Ll.  If there is more
!    than one shortest path, then enter one into a list L1, and the
!    rest into another list L2.  Set J = 2.
!
!    STEP 2: For each I = 1, 2, ..., J-1, find all Q(I,J) and place
!    them into the list L2.
!
!    STEP 3:  Take a path from the end of the list L2.  Denote this path
!    by P(J) and move it from L2 to L1.  If J = K then stop; the required
!    list of K shortest paths is in L1; otherwise, set J = J + 1 and
!    return to Step 2.
!
!    The algorithm requires O(kn3) operations if all edge lengths of the
!    input digraph are nonnegative, and 0(kn') operations if some edge lengths
!    are negative.
!
!  Modified:
!
!    23 July 2000
!
!  Author:
!
!    Hang Tong Lau
!
!  Reference:
!
!    Hang Tong Lau,
!    Algorithms on Graphs,
!    Tab Books, 1989,
!    QA166 L38.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) NEDGE, the number of edges.
!
!    Input, integer ( kind = 4 ) INODE(NEDGE), JNODE(NEDGE); the I-th edge
!    is directed from INODE(I) to JNODE(I).
!
!    Input, integer ( kind = 4 ) ARCLEN(NEDGE); ARCLEN(I) is the length of the I-th edge.
!
!    Input, integer ( kind = 4 ) KPATHS, the requested number of shortest paths.
!    KPATHS must be greater than 0.
!
!    Input, integer ( kind = 4 ) MAXQUE, the size of the auxilliary storage.
!    In most cases, MAXQUE = 2 * KPATHS is sufficient.
!
!    Input, integer ( kind = 4 ) ISORCE, ISINK; the paths to be found must go from
!    node ISORCE to node ISINK.
!
!    Output, integer ( kind = 4 ) IFLAG, completion flag.
!    0, for normal termination;
!    1, if no paths are found;
!    2, if the size of the auxiliary arrays CROSAR and NEXTRD are not large
!    enough; in this case, increase their sizes by using a larger value of
!    MAXQUE and rerun the program.
!
!    Output, integer ( kind = 4 ) IPATHS, the actual number of shortest paths found
!    by KSHOT2; IPATHS <= KPATHS.
!
!    Output, integer ( kind = 4 ) PTHLEN(KPATHS+3); the shortest path lengths are stored
!    in PTHLEN(2), PTHLEN(3), ..., PTHLEN(KPATHS+1).
!
!    Workspace, integer ARCDIR(NNODE); array of shortest path tree.
!
!    Workspace, integer TRDIST(NNODE); array of node information.
!
!    Workspace, integer ARCNOD(NNODE); array of node information.
!
!    Workspace, integer ARCFWD(NNODE); first edge in forward star.
!
!    Workspace, integer ARCBWD(NNODE); first edge in backward star.
!
!    Workspace, integer AUXSTG(NNODE); auxiliary storage.
!
!    Workspace, integer AUXDIS(NNODE); auxiliary node array.
!
!    Workspace, integer AUXTRE(NNODE); auxiliary node array.
!
!    Workspace, integer AUXLNK(NNODE); auxiliaty node array.
!
!    Workspace, integer NXTFWD(NEDGE); links of node chain.
!
!    Workspace, integer NXTBWD(NEDGE); links of node chain.
!
!    Workspace, integer QUFIRP(KPATHS+3); description of the first path in
!    a queue.
!
!    Workspace, integer QUNXTP(KPATHS+3); next entry in a queue.
!
!    Workspace, integer CROSAR(MAXQUE); cross edge of the path.
!
!    Workspace, integer NEXTRD(MAXQUE); next entry of the pool storage.
!
  implicit none

  integer ( kind = 4 ) kpaths
  integer ( kind = 4 ) maxque
  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) arcbwd(nnode)
  integer ( kind = 4 ) arcdir(nnode)
  integer ( kind = 4 ) arcfwd(nnode)
  integer ( kind = 4 ) arclen(nedge)
  integer ( kind = 4 ) arcnod(nnode)
  integer ( kind = 4 ) auxdis(nnode)
  integer ( kind = 4 ) auxlnk(nnode)
  integer ( kind = 4 ) auxstg(nnode)
  integer ( kind = 4 ) auxtre(nnode)
  integer ( kind = 4 ) crosar(maxque)
  logical finem1
  logical finem2
  logical forwrd
  logical goon
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iauxd1
  integer ( kind = 4 ) iauxd2
  integer ( kind = 4 ) iauxd3
  integer ( kind = 4 ) ibedge
  integer ( kind = 4 ) icall
  integer ( kind = 4 ) icedge
  integer ( kind = 4 ) icnt
  integer ( kind = 4 ) idedge
  integer ( kind = 4 ) idet
  integer ( kind = 4 ) idet1
  integer ( kind = 4 ) iedge
  integer ( kind = 4 ) iedge1
  integer ( kind = 4 ) iedge2
  integer ( kind = 4 ) iexam
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) ihd1
  integer ( kind = 4 ) ihd2
  integer ( kind = 4 ) ihead
  integer ( kind = 4 ) ihigh
  integer ( kind = 4 ) ilen1
  integer ( kind = 4 ) ilen2
  integer ( kind = 4 ) incrs1
  integer ( kind = 4 ) incrs2
  integer ( kind = 4 ) incrs3
  integer ( kind = 4 ) incrs4
  integer ( kind = 4 ) index1
  integer ( kind = 4 ) index2
  integer ( kind = 4 ) index3
  integer ( kind = 4 ) index4
  logical initsp
  integer ( kind = 4 ) inode(nedge)
  integer ( kind = 4 ) inqlas
  integer ( kind = 4 ) inrd1
  integer ( kind = 4 ) inrd2
  integer ( kind = 4 ) iop1
  integer ( kind = 4 ) iop2
  integer ( kind = 4 ) iop3
  integer ( kind = 4 ) iorder
  integer ( kind = 4 ) iparm
  integer ( kind = 4 ) ipaths
  integer ( kind = 4 ) ipool1
  integer ( kind = 4 ) ipool2
  integer ( kind = 4 ) ipool3
  integer ( kind = 4 ) iqufir
  integer ( kind = 4 ) isink
  integer ( kind = 4 ) isorce
  integer ( kind = 4 ) itail
  integer ( kind = 4 ) itrear
  integer ( kind = 4 ) iupbnd
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j3
  integer ( kind = 4 ) jedge
  integer ( kind = 4 ) jem1
  integer ( kind = 4 ) jem2
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) jnd1
  integer ( kind = 4 ) jnd2
  integer ( kind = 4 ) jnode(nedge)
  integer ( kind = 4 ) jterm
  integer ( kind = 4 ) jump
  integer ( kind = 4 ) large
  logical lasta
  logical lastb
  logical lastno
  integer ( kind = 4 ) length
  integer ( kind = 4 ) lenrst
  integer ( kind = 4 ) lentab
  integer ( kind = 4 ) linkst
  logical loopon
  integer ( kind = 4 ) low
  integer ( kind = 4 ) mark1
  integer ( kind = 4 ) mark2
  integer ( kind = 4 ) mshade
  integer ( kind = 4 ) ncrs1
  integer ( kind = 4 ) ncrs2
  integer ( kind = 4 ) ncrs3
  integer ( kind = 4 ) nextrd(maxque)
  integer ( kind = 4 ) njem1
  integer ( kind = 4 ) njem2
  integer ( kind = 4 ) njem3
  integer ( kind = 4 ) node1
  integer ( kind = 4 ) node2
  integer ( kind = 4 ) nodep1
  integer ( kind = 4 ) nodep2
  integer ( kind = 4 ) nodep3
  logical noroom
  logical nostg
  integer ( kind = 4 ) nqfirs
  integer ( kind = 4 ) nqin
  integer ( kind = 4 ) nqlast
  integer ( kind = 4 ) nqop1
  integer ( kind = 4 ) nqop2
  integer ( kind = 4 ) nqout1
  integer ( kind = 4 ) nqout2
  integer ( kind = 4 ) nqp1
  integer ( kind = 4 ) nqp2
  integer ( kind = 4 ) nqsize
  integer ( kind = 4 ) nqtab(kpaths+3)
  integer ( kind = 4 ) nrdsiz
  integer ( kind = 4 ) nxtbwd(nedge)
  integer ( kind = 4 ) nxtfwd(nedge)
  integer ( kind = 4 ) pathlen(kpaths+3)
  integer ( kind = 4 ) qufirp(kpaths+3)
  integer ( kind = 4 ) qunxtp(kpaths+3)
  logical rdfull
  integer ( kind = 4 ) tqufir
  integer ( kind = 4 ) tqunxt
  integer ( kind = 4 ) tquord
  integer ( kind = 4 ) trdist(nnode)
!
!  Initialize.
!
  finem1 = .false.
  incrs4 = 0
  linkst = 0
  ncrs1 = 0
  ncrs3 = 0
!
!  Set up the network representation.
!
  iflag = 0
  arcfwd(1:nnode) = 0
  arcbwd(1:nnode) = 0

  large = 1 + sum ( arclen )

  do i = 1, nedge
    itail = inode(i)
    ihead = jnode(i)
    nxtfwd(i) = arcfwd(itail)
    arcfwd(itail) = i
    nxtbwd(i) = arcbwd(ihead)
    arcbwd(ihead) = i
  end do
!
!  Initialize.
!
  auxdis(1:nnode) = large
  qufirp(1:kpaths+3) = 0

  call i4vec_indicator ( maxque-1, nextrd )
  nextrd(maxque) = 0
!
!  Build the shortest distance tree.
!
!  TRDIST(I) is used to store the shortest distance of node I from ISORCE;
!  ARCDIR(I) will contain the tree edge coming to node I; it is negative
!  if the direction of the edge is towards ISORCE, and it is zero if not
!  reachable.
!
  trdist(1:nnode) = large
  arcdir(1:nnode) = 0
  arcnod(1:nnode) = 0
  trdist(isorce) = 0
  arcdir(isorce) = isorce
  arcnod(isorce) = isorce
  j = isorce
  node1 = isorce
!
!  Examine neighbors of node J.
!
70    continue

  iedge1 = arcfwd(j)
  forwrd = .true.
  lastno = .false.

80    continue

  if ( iedge1 == 0 ) then

    lastno = .true.

  else

    length = trdist(j) + arclen(iedge1)

    if ( forwrd ) then
      node2 = jnode(iedge1)
      iedge2 = iedge1
    else
      node2 = inode(iedge1)
      iedge2 = - iedge1
    end if

    if ( length < trdist(node2) ) then

      trdist(node2) = length
      arcdir(node2) = iedge2

      if ( arcnod(node2) == 0 ) then

        arcnod(node1) = node2
        arcnod(node2) = node2
        node1 = node2

      else

        if ( arcnod(node2) < 0 ) then

          arcnod(node2) = arcnod(j)
          arcnod(j) = node2

          if ( node1 == j ) then
            node1 = node2
            arcnod(node2) = node2
          end if

        end if

      end if

    end if

    if ( forwrd ) then
      iedge1 = nxtfwd(iedge1)
    else
      iedge1 = nxtbwd(iedge1)
    end if

  end if

  if ( .not. lastno ) then
    go to 80
  end if

  jj = j
  j = arcnod(j)
  arcnod(jj) = - 1

  if ( j /= jj ) then
    go to 70
  end if
!
!  Finish building the shortest distance tree.
!
  ipaths = 0
  noroom = .false.

  if ( arcdir(isink) == 0 ) then
    iflag = 1
    return
  end if
!
!  Initialize the storage pool.
!
  call i4vec_indicator ( kpaths+2, qunxtp )
  qunxtp(kpaths+3) = 0
!
!  Initialize the priority queue.
!
  lentab = kpaths + 3
  low = - large
  ihigh = large
  nqop1 = lentab
  nqop2 = 0
  nqsize = 0
  nqin = 0
  nqout1 = 0
  nqout2 = 0
!
!  Obtain an entry from storage pool.
!
  index1 = qunxtp(1)
  qunxtp(1) = qunxtp(index1+1)
  index2 = qunxtp(1)
  qunxtp(1) = qunxtp(index2+1)
  pathlen(index1+1) = low
  qunxtp(index1+1) = index2
  pathlen(index2+1) = ihigh
  qunxtp(index2+1) = 0
  nqp1 = 0
  nqp2 = 1
  nqtab(1) = index1
  nqtab(2) = index2
  nqfirs = ihigh
  nqlast = low
!
!  Set the shortest path to the queue.
!
  ipool1 = qunxtp(1)
  qunxtp(1) = qunxtp(ipool1+1)
  ipool2 = ipool1
  incrs1 = nextrd(1)
  nextrd(1) = nextrd(incrs1+1)
  crosar(incrs1+1) = arcdir(isink)
  nextrd(incrs1+1) = 0
  pathlen(ipool1+1) = trdist(isink)
  qufirp(ipool1+1) = incrs1
  iparm = ipool1
  icall = 0

100   continue

  iorder = pathlen(iparm+1)
  iop1 = nqp1
  iop2 = nqp2
!
!  Insert IPARM into the priority queue.
!
  do while ( 1 < iop2 - iop1 )

    iop3 = ( iop1 + iop2 ) / 2

    if ( pathlen(nqtab(iop3+1)+1) < iorder ) then
      iop1 = iop3
    else
      iop2 = iop3
    end if

  end do
!
!  Linear search starting from NQTAB(IOP1+1).
!
  index1 = nqtab(iop1+1)

120   continue

  index2 = index1
  index1 = qunxtp(index1+1)

  if ( pathlen(index1+1) <= iorder ) then
    go to 120
  end if
!
!  Insert between INDEX1 and INDEX2.
!
  qunxtp(index2+1) = iparm
  qunxtp(iparm+1) = index1
!
!  Update data in the queue.
!
  nqsize = nqsize + 1
  nqin = nqin + 1
  nqop1 = nqop1 - 1

  if ( nqsize == 1 ) then

    nqfirs = iorder
    nqlast = iorder

  else

    if ( nqlast < iorder ) then

      nqlast = iorder

    else

      if ( iorder < nqfirs ) then
       nqfirs = iorder
      end if

    end if

  end if

  if ( nqop1 <= 0 ) then
!
!  Reorganize.
!
    index1 = nqtab(nqp1+1)
    nqtab(1) = index1
    nqp1 = 0
    index2 = nqtab(nqp2+1)
    j3 = nqsize / lentab
    j2 = j3 + 1
    j1 = mod ( nqsize, lentab )

    do iop2 = 1, j1
      do i = 1, j2
        index1 = qunxtp(index1+1)
      end do
      nqtab(iop2+1) = index1
    end do

    if ( 0 < j3 ) then

      iop2 = j1 + 1

      do while ( iop2 <= lentab-1 )

        do i = 1, j3
          index1 = qunxtp(index1+1)
        end do

        nqtab(iop2+1) = index1
        iop2 = iop2 + 1

      end do

    end if

    nqp2 = iop2
    nqtab(nqp2+1) = index2
    nqop2 = nqop2 + 1
    nqop1 = nqsize / 2

    if ( nqop1 < lentab ) then
      nqop1 = lentab
    end if

  end if

  if ( 0 < icall ) then
    go to 430
  end if

  ilen1 = 0
  mark1 = 0
  initsp = .true.
  arcnod(1:nnode) = 0
!
!  Process the next path.
!
180   continue

  mark1 = mark1 + 2
  mark2 = mark1
  mshade = mark1 + 1
!
!  Obtain the first entry from the priority queue.
!
  if ( 0 < nqsize ) then

    index2 = nqtab(nqp1+1)
    index1 = qunxtp(index2+1)
    qunxtp(index2+1) = qunxtp(index1+1)
    nqfirs = pathlen(qunxtp(index1+1)+1)

    if ( index1 == nqtab(nqp1+2) ) then
      nqp1 = nqp1 + 1
      nqtab(nqp1+1) = index2
    end if

    nqop1 = nqop1 - 1
    nqsize = nqsize - 1
    nqout1 = nqout1 + 1
    ipool3 = index1

  else

    ipool3 = 0

  end if
!
!  No more paths in the queue; stop.
!
  if ( ipool3 == 0 ) then
    noroom = noroom .and. ipaths < kpaths
    go to 450
  end if

  qunxtp(ipool2+1) = ipool3
  ipool2 = ipool3
  ipaths = ipaths + 1

  if ( kpaths < ipaths ) then
    noroom = .false.
    ipaths = ipaths - 1
    go to 450
  end if

  ilen2 = pathlen(ipool3+1)
  iqufir = qufirp(ipool3+1)

  if ( ilen2 < ilen1 ) then
    go to 450
  end if

  ilen1 = ilen2
!
!  Examine the tail of the edge.
!
  incrs2 = iqufir
  ncrs2 = isorce
  nodep1 = nnode + 1
!
!  Obtain data of next path.
!
190   continue

  jump = 1
  go to 500

200   continue

  ilen2 = ilen2 - arclen(incrs4)
  nodep1 = nodep1 - 1
  auxstg(nodep1) = incrs1

  do while ( ncrs1 /= ncrs3 )

    j = abs ( arcdir(ncrs1) )
    ilen2 = ilen2 - arclen(j)

    if ( 0 < arcdir(ncrs1) ) then
      ncrs1 = inode(j)
    else
      ncrs1 = jnode(j)
    end if

  end do

  if ( .not. finem1 ) then
    go to 190
  end if
!
!  Store the tail of the edge.
!
  nodep2 = nodep1
  finem2 = finem1
!
!  Obtain data of next path.
!
  jump = 2
  go to 500

220   continue

  if ( linkst == 2 ) then

    nodep2 = nodep2 - 1
    auxstg(nodep2) = incrs4
    arclen(incrs4) = arclen(incrs4) + large
    finem2 = finem1
!
!  Obtain data of next path.
!
    jump = 2
    go to 500

  end if
!
!  Close the edge on the shortest path.
!
  finem2 = finem2 .and. linkst /= 3

  if ( finem2 ) then
    iedge = abs ( arcdir(ncrs3) )
    nodep2 = nodep2 - 1
    auxstg(nodep2) = iedge
    arclen(iedge) = arclen(iedge) + large
  end if
!
!  Mark more nodes.
!
230   continue

  if ( linkst /= 3 ) then

240     continue

    arcnod(ncrs2) = mark2

    do while ( ncrs1 /= ncrs3 )

      arcnod(ncrs1) = mark2

      if ( 0 < arcdir(ncrs1) ) then
        ncrs1 = inode(arcdir(ncrs1))
      else
        ncrs1 = jnode(-arcdir(ncrs1))
      end if

    end do

    jump = 3
    go to 500

260     continue

    if ( linkst == 1 ) then
      go to 240
    end if

270     continue

    if ( linkst== 2 ) then
     jump = 4
      go to 500
    end if

    go to 230

  end if
!
!  Generate descendants of the tail of the edge.
!
  nodep3 = nodep1
  incrs1 = auxstg(nodep3)
  jnd1 = crosar(incrs1+1)
!
!  Obtain the first node of the edge traversing forward.
!
  if ( jnd1 < 0 ) then
    jnd2 = inode(-jnd1)
  else
    jnd2 = jnode(jnd1)
  end if
!
!  Process a section.
!
280   continue

  nodep3 = nodep3 + 1
  jterm = jnd2
  jedge = jnd1

  if ( nnode < nodep3 ) then

    jnd2 = isorce

  else

    incrs2 = auxstg(nodep3)
    jnd1 = crosar(incrs2+1)

    if ( 0 < -jnd1 ) then
      jnd2 = inode(-jnd1)
    else
      jnd2 = jnode(jnd1)
    end if

  end if
!
!  Process a node.
!
290   continue

  mark1 = mark1 + 2
  itrear = mark1
  iexam = mark1 + 1
  iedge = abs ( jedge )
  arclen(iedge) = arclen(iedge) + large

  if ( initsp ) then
    initsp = nqin < kpaths
  end if

  if ( initsp ) then
    iupbnd = large
  else
    iupbnd = nqlast
  end if
!
!  Obtain the restricted shortest path from ISORCE to JTERM.
!
  lenrst = iupbnd
  ibedge = 0
  auxdis(jterm) = 0
  auxtre(jterm) = 0
  auxlnk(jterm) = 0
  jem1 = jterm
  jem2 = jem1
!
!  Examine the next node.
!
300   continue

  njem1 = jem1
  iauxd1 = auxdis(njem1)
  jem1 = auxlnk(njem1)
  arcnod(njem1) = itrear

  if ( lenrst <= iauxd1 + trdist(njem1) + ilen2 ) then
    go to 360
  end if

  goon = .true.
  lasta = .false.
  iedge1 = arcbwd(njem1)
!
!  Loop through edge from NJEMI.
!
310   continue

  if ( iedge1 == 0 ) then

    lasta = .true.

  else
!
!  Process the edge IEDGE1.
!
    iauxd2 = iauxd1 + arclen(iedge1)

    if ( goon ) then
      njem2 = inode(iedge1)
      iedge2 = iedge1
      iedge1 = nxtbwd(iedge1)
    else
      njem2 = jnode(iedge1)
      iedge2 = - iedge1
      iedge1 = nxtfwd(iedge1)
    end if

    if ( arcnod(njem2) /= mark2 ) then

      iauxd3 = iauxd2 + ilen2 + trdist(njem2)

      if ( lenrst <= iauxd3 ) then
        go to 350
      end if

      if ( arcnod(njem2) < mark2 ) then

        if ( arcdir(njem2) + iedge2 == 0 ) then
          arcnod(njem2) = mshade
          go to 340
        end if
!
!  Examine the status of the path.
!
        loopon = .true.
        njem3 = njem2

        do while ( loopon .and. njem3 /= isorce )

          if ( arcnod(njem3) < mark2 ) then

            j = arcdir(njem3)

            if ( 0 < j ) then
              njem3 = inode(j)
            else
              njem3 = jnode(-j)
            end if

          else

            loopon = .false.

          end if

        end do
!
!  Better path found.
!
        if ( loopon ) then

          lenrst = iauxd3
          ibedge = iedge2
          go to 350

        else

          njem3 = njem2
          lastb = .false.

330           continue

          if ( arcnod(njem3) < mark2 ) then

            arcnod(njem3) = mshade
            j = arcdir(njem3)

            if ( 0 < j ) then
              njem3 = inode(j)
            else
              njem3 = jnode(-j)
            end if

          else

            lastb = .true.

          end if

          if ( .not. lastb ) then
            go to 330
          end if

        end if

      end if

340       continue

      if ( arcnod(njem2) < itrear .or. iauxd2 < auxdis(njem2) ) then
!
!  Update node NJEM2.
!
        auxdis(njem2) = iauxd2
        auxtre(njem2) = iedge2

        if ( arcnod(njem2) /= iexam ) then

          arcnod(njem2) = iexam

          if ( jem1 == 0 ) then

            jem1 = njem2
            jem2 = njem2
            auxlnk(njem2) = 0

          else

            if ( arcnod(njem2) == itrear ) then
              auxlnk(njem2) = jem1
              jem1 = njem2
            else
              auxlnk(njem2) = 0
              auxlnk(jem2) = njem2
              jem2 = njem2
            end if

          end if

        end if

      end if

    end if

  end if

350   continue

  if ( .not. lasta ) then
    go to 310
  end if

360   continue

  if ( 0 < jem1 ) then
    go to 300
  end if

  arcnod(jterm) = mark2
!
!  Finish processing the restricted path.
!
  if ( ibedge /= 0 .and. lenrst < iupbnd ) then

    idet = 0
    icedge = ibedge

370     continue

    if ( 0 < icedge ) then
      idedge = jnode(icedge)
    else
      idedge = inode(-icedge)
    end if

    if ( icedge /= arcdir(idedge) .or. idedge == jterm ) then
      idet = idet + 1
      auxstg(idet) = icedge
    end if

    icedge = auxtre(idedge)

    if ( icedge /= 0 ) then
      go to 370
    end if
!
!  Restore the path data.  NOSTG will be TRUE if the arrays CROSAR
!  and NEXTRD need more space.
!
    idet1 = idet
    inrd1 = nextrd(1)
    inqlas = large
    nostg = .false.

    do while ( 0 < idet1 .and. 0 < inrd1 )
      idet1 = idet1 - 1
      inrd1 = nextrd(inrd1+1)
    end do

    rdfull = ( .not. initsp ) .and. ( kpaths <= ipaths + nqsize )

390     continue

    if ( rdfull .or. ( 0 < idet1 ) ) then
!
!  Remove the last path from the queue.
!
      inqlas = nqlast
      noroom = .true.
      rdfull = .false.
!
!  Get the last entry from the priority queue.
!
      if ( 0 < nqsize ) then

        index4 = nqtab(nqp2+1)
        index3 = nqtab(nqp2)

        if ( qunxtp(index3+1) == index4 ) then
          nqp2 = nqp2 - 1
          nqtab(nqp2+1) = index4
          index3 = nqtab(nqp2)
        end if

        index2 = index3

        do while ( index3 /= index4 )
          index1 = index2
          index2 = index3
          index3 = qunxtp(index3+1)
        end do

        qunxtp(index1+1) = index4
        nqlast = pathlen(index1+1)
        nqop1 = nqop1 - 1
        nrdsiz = index2
        nqsize = nqsize - 1
        nqout2 = nqout2 + 1

      else

        nrdsiz = 0

      end if

      if ( nrdsiz == 0 ) then
        nostg = .true.
        go to 430
      end if

      inrd1 = qufirp(nrdsiz+1)

      do while ( 0 < inrd1 )
        j = inrd1 + 1
        idet1 = idet1 - 1
        inrd2 = inrd1
        inrd1 = nextrd(j)
        nextrd(j) = nextrd(1)
        nextrd(1) = inrd2
      end do
!
!  Put the entry NRDSIZ into storage pool.
!
      qunxtp(nrdsiz+1) = qunxtp(1)
      qunxtp(1) = nrdsiz
      go to 390

    end if
!
!  Build the entries of CROSAR and NEXTRD.
!
    if ( inqlas <= lenrst ) then
      go to 430
    end if

    inrd2 = - incrs1
    idet1 = idet

    do while ( 0 < idet1 )
      inrd1 = nextrd(1)
      nextrd(1) = nextrd(inrd1+1)
      crosar(inrd1+1) = auxstg(idet1)
      nextrd(inrd1+1) = inrd2
      inrd2 = inrd1
      idet1 = idet1 - 1
    end do
!
!  Obtain the entry NRDSlZ from storage pool.
!
    nrdsiz = qunxtp(1)
    qunxtp(1) = qunxtp(nrdsiz+1)
    pathlen(nrdsiz+1) = lenrst
    qufirp(nrdsiz+1) = inrd2
    iparm = nrdsiz
    icall = 1
    go to 100

430     continue

    if ( nostg ) then
      iflag = 2
      go to 450
    end if

  end if

  arclen(iedge) = arclen(iedge) - large
  ilen2 = ilen2 + arclen(iedge)

  if ( jterm /= jnd2 ) then

    if ( 0 < jedge ) then
      jterm = inode(jedge)
    else
      jterm = jnode(-jedge)
    end if

    jedge = arcdir(jterm)

  end if

  if ( jterm /= jnd2 ) then
    go to 290
  end if

  incrs1 = incrs2

  if ( nodep3 <= nnode ) then
    go to 280
  end if
!
!  Restore the join edges.
!
  do while ( nodep2 <= nodep1 - 1 )
    j = auxstg(nodep2)
    arclen(j) = arclen(j) - large
    nodep2 = nodep2 + 1
  end do
!
!  Repeat with the next path.
!
  go to 180
!
!  Sort the paths.
!
450   continue

  ihd2 = ipool1
  icnt = 0

  do

    ihd1 = ihd2
    icnt = icnt + 1
    ihd2 = qunxtp(ihd1+1)
    qunxtp(ihd1+1) = icnt

    if ( ihd1 == ipool2 ) then
      exit
    end if

  end do
!
!  Release all queue entries to the storage pool.
!
  j = nqtab(nqp2+1)
  qunxtp(j+1) = qunxtp(1)
  qunxtp(1) = nqtab(nqp1+1)
  nqp1 = 0
  nqp2 = 0
  ihd2 = 0

  do

    j = ihd2 + 1
    ihd2 = qunxtp(j)
    qunxtp(j) = 0

    if ( ihd2 == 0 ) then
      exit
    end if

  end do
!
!  Exchanging records.
!
  jj = kpaths + 2

  do i = 1, jj

    do while ( 0 < qunxtp(i+1) .and. qunxtp(i+1) /= i )
      tquord = pathlen(i+1)
      tqufir = qufirp(i+1)
      tqunxt = qunxtp(i+1)

      j = qunxtp(i+1) + 1

      pathlen(i+1) = pathlen(j)
      qufirp(i+1) = qufirp(j)
      qunxtp(i+1) = qunxtp(j)

      pathlen(tqunxt+1) = tquord
      qufirp(tqunxt+1) = tqufir
      qunxtp(tqunxt+1) = tqunxt
    end do

  end do

  pathlen(1) = isorce
  qufirp(1) = isink
  qunxtp(1) = ipaths
  return
!
!  Obtain data for the next path.
!
500   continue

  if ( incrs2 == 0 ) then

    linkst = 3

  else

    ncrs3 = ncrs2
    incrs1 = incrs2
    j = abs ( incrs1 ) + 1
    incrs3 = crosar(j)
    incrs2 = nextrd(j)

    if ( 0 < incrs3 ) then
      ncrs1 = inode(incrs3)
      ncrs2 = jnode(incrs3)
      incrs4 = incrs3
    else
      ncrs1 = jnode(-incrs3)
      ncrs2 = inode(-incrs3)
      incrs4 = - incrs3
    end if

    finem1 = incrs2 <= 0

    if ( ncrs2 == ncrs3 ) then
      linkst = 2
    else
      linkst = 1
    end if

  end if

  if ( jump == 1 ) then
    go to 200
  else if ( jump == 2 ) then
    go to 220
  else if ( jump == 3 ) then
    go to 260
  else if ( jump == 4 ) then
    go to 270
  end if

  return
end
subroutine digraph_arc_mineqv ( nnode, nedge, inode, jnode, arclis )

!*****************************************************************************80
!
!! DIGRAPH_ARC_MINEQV finds a minimal equivalent of a strongly connected digraph.
!
!  Warning:
!
!    This routine is failing on the sample problem!
!
!  Description:
!
!    The minimal equivalent digraph problem is to find a digraph H
!    from a given strongly connected digraph G by removing the maximum number
!    of edges from G without affecting its reachability properties.  That is,
!    for any two nodes U and V in G, there is a directed path from node U
!    to node V in H.
!
!  Method:
!
!    Let E be the set of edges in a given strongly connected digraph
!    G of NNODE nodes and NEDGE edges.
!
!    STEP 1. Examine all edges sequentially.  An edge (I,J) is removed
!    from G whenever there exists an alternative path from node I to
!    node J which does not include previously eliminated edges.
!
!    Let F be the set of edges which are removed from G, and S = E - F
!    be the current solution.
!
!    STEP 2. Iteratively execute a combination of backtracking and forward
!    moves.  A backtracking move will remove from F the highest labelled
!    edge - say the K-th one - and add it back in E.  This backtracking move
!    is followed by a forward move which will sequentially consider all edges
!    in E with labels greater than K, removing from E every edge (I,J) for
!    which an alternative path exists, and adding (I,J) back to F.
!
!    STEP 3. If the number of edges in E is less than the number of edges
!    in S, then update the current solution by setting S = E, and return
!    to Step 2.
!
!    The algorithm terminates when no further backtracking is possible
!    in Step 2 (in which case F is empty) or when the number of edges in
!    S is NNODE.
!
!    The number of iterations of the algorithm is bounded by O(2^NEDGE).
!
!  Modified:
!
!    23 July 2000
!
!  Author:
!
!    Hang Tong Lau
!
!  Reference:
!
!    Hang Tong Lau,
!    Algorithms on Graphs,
!    Tab Books, 1989, pages 47-59,
!    QA166 L38.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) NEDGE, the number of edges.
!
!    Input, integer ( kind = 4 ) INODE(NEDGE), JNODE(NEDGE); the I-th edge of the digraph is
!    directed from node INODE(I) to node JNODE(I).  It is assumed that the
!    input digraph is strongly connected.
!
!    Output, logical ARCLIS(NEDGE); ARCLIS(I) has the value TRUE if the I-th
!    edge is in the minimal equivalent digraph; otherwise, it is FALSE.
!
!    Workspace, integer FWDARC(NEDGE); FWDARC(I) is the ending node of the I-th
!    edge in the forward star representation of the digraph.
!
!    Workspace, integer ARCFIR(NNODE+1); ARCFIR(I) is the number of the first
!    edge starting at node I in the forward star representation of the digraph.
!
!    Workspace, integer POINT(NEDGE); pointing to the original edge list in the
!    forward star representation.
!
!    Workspace, logical MARK(NEDGE), a working copy of the array ARCLIS.
!
!    Workspace, integer IDESCN(NNODE); IDESCN(I) is the number of descendants of
!    node I.
!
!    Workspace, integer IANCES(NNODE); IANCES(I) is the number of ancestors of
!    node I.
!
  implicit none

  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) arcfir(nnode+1)
  logical arclis(nedge)
  integer ( kind = 4 ) fwdarc(nedge)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) iances(nnode)
  integer ( kind = 4 ) idescn(nnode)
  integer ( kind = 4 ) iedges
  integer ( kind = 4 ) ierase
  integer ( kind = 4 ) ihigh
  integer ( kind = 4 ) inode(nedge)
  integer ( kind = 4 ) iup
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) jerase
  integer ( kind = 4 ) jnode(nedge)
  logical join
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kedge
  integer ( kind = 4 ) low
  integer ( kind = 4 ) point(nedge)
  logical mark(nedge)
  logical pexist
!
!  Set up the forward star representation of the digraph.
!
  k = 0

  do i = 1, nnode
    arcfir(i) = k + 1
    do j = 1, nedge
      if ( inode(j) == i ) then
        k = k + 1
        point(k) = j
        fwdarc(k) = jnode(j)
        mark(k) = .TRUE.
      end if
    end do
  end do

  arcfir(nnode+1) = nedge + 1
!
!  Compute the number of descendants and ancestors of each node.
!
  idescn(1:nnode) = 0
  iances(1:nnode) = 0

  iedges = 0

  do k = 1, nedge
    i = inode(k)
    j = jnode(k)
    idescn(i) = idescn(i) + 1
    iances(j) = iances(j) + 1
    iedges = iedges + 1
  end do
!
!  If NNODE edges left, then we know we have to stop.
!
  if ( iedges == nnode ) then

    do k = 1, nedge
      arclis(point(k)) = mark(k)
    end do

    return

  end if

  ierase = 0

  do k = 1, nedge

    i = inode(point(k))
    j = jnode(point(k))
!
!  Check for the existence of an alternative path.
!
    if ( idescn(i) /= 1 ) then
      if ( iances(j) /= 1 ) then
        mark(k) = .FALSE.

        call digraph_arc_find_path ( nnode, nedge, i, j, fwdarc, arcfir, &
          mark, pexist )

        if ( pexist ) then
          idescn(i) = idescn(i) - 1
          iances(j) = iances(j) - 1
          ierase = ierase + 1
        else
          mark(k) = .TRUE.
        end if

      end if
    end if

  end do

  if ( ierase == 0 ) then
    do k = 1, nedge
      arclis(point(k)) = mark(k)
    end do
    return
  end if

  ihigh = 0
  i1 = nnode
  j1 = nnode
!
!  Store the current best solution.
!
80    continue

  arclis(1:nedge) = mark(1:nedge)

  jerase = ierase

  if ( iedges - jerase == nnode ) then

    mark(1:nedge) = arclis(1:nedge)

    do k = 1, nedge
      arclis(point(k)) = mark(k)
    end do

    return

  end if
!
!  Forward move.
!
120   continue

  join = .FALSE.
  low = arcfir(i1)
  iup = arcfir(i1+1) - 1

  do k = low, iup
    if ( fwdarc(k) == j1 ) then
      join = .TRUE.
      kedge = k
      exit
    end if
  end do

  if ( join ) then

    if ( .not. mark(kedge) ) then

      mark(kedge) = .TRUE.
      idescn(i1) = idescn(i1) + 1
      iances(j1) = iances(j1) + 1
      ierase = ierase - 1

      if ( jerase < ierase + ihigh - (nnode - i1) ) then
        go to 220
      end if

    end if

    ihigh = ihigh + 1

  end if

150   continue

  if ( j1 /= 1 ) then
    j1 = j1 - 1
    go to 120
  end if

  if ( i1 == 1 ) then

    mark(1:nedge) = arclis(1:nedge)

    do k = 1, nedge
      arclis(point(k)) = mark(k)
    end do

    return

  end if

  i1 = i1 - 1
  j1 = nnode
  go to 120
!
!  Backtrack move.
!
180   continue

  join = .FALSE.
  low = arcfir(i1)
  iup = arcfir(i1+1) - 1

  do k = low, iup
    if ( fwdarc(k) == j1 ) then
      join = .TRUE.
      kedge = k
      exit
    end if
  end do

  if ( join ) then

    ihigh = ihigh - 1

    if ( idescn(i1) /= 1 ) then
      if ( iances(j1) /= 1 ) then
        mark(kedge) = .FALSE.

        call digraph_arc_find_path ( nnode, nedge, i1, j1, fwdarc, &
          arcfir, mark, pexist )

        if ( pexist ) then
          idescn(i1) = idescn(i1) - 1
          iances(j1) = iances(j1) - 1
          ierase = ierase + 1
          go to 210
        end if

        mark(kedge) = .TRUE.

      end if
    end if

    if ( ierase + ihigh - (nnode - i1) <= jerase ) then
      ihigh = ihigh + 1
      go to 150
    end if
!
!  Check for the termination of the forward move.
!
210     continue

    if ( ihigh == nnode - i1 ) then
      go to 80
    end if

  end if

220   continue

  if ( j1 /= nnode ) then
    j1 = j1 + 1
    go to 180
  end if

  if ( i1 /= nnode ) then
    i1 = i1 + 1
    j1 = 1
    go to 180
  end if

  go to 80
end
subroutine digraph_arc_nflow ( nnode, nedge, inode, jnode, capac, isorce, &
  isink, mincut, flow, node_flow )

!*****************************************************************************80
!
!! DIGRAPH_ARC_NFLOW implements Karzanov's network flow algorithm.
!
!  Warning:
!
!    This routine violates some common rules of FORTRAN programming
!    by jumping into blocks of code.  Some of these have been reprogrammed,
!    but there is still a jump into a block, via statement 160.
!
!  Modified:
!
!    23 July 2000
!
!  Author:
!
!    Hang Tong Lau
!
!  Reference:
!
!    Hang Tong Lau,
!    Algorithms on Graphs,
!    Tab Books, 1989,
!    QA166 L38.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) NEDGE, the number of edges in the network.
!
!    Input, integer ( kind = 4 ) INODE(NEDGE), JNODE(NEDGE); these arrays
!    represent the edges of the digraph.  If there is an edge between nodes
!    U and V, it must be represented TWICE.  There will be indices I and J
!    such that
!      INODE(I) = U, JNODE(I) = V, and
!      INODE(J) = V, JNODE(I) = U.
!
!    Input, integer ( kind = 4 ) CAPAC(NEDGE); since a given edge is represented
!    twice, there are two indices of CAPAC that reference the same edge.
!    One of these indices should be set to the capacity of the edge,
!    and the other to zero.
!
!    Input, integer ( kind = 4 ) ISORCE, ISINK, the indices of the source and
!    sink nodes, which must be distinct.
!
!    Output, integer ( kind = 4 ) MINCUT(NNODE); MINCUT(I) = 1 if node I is in
!    the minimal cut set; otherwise, it is 0.
!
!    Output, integer ( kind = 4 ) FLOW(NEDGE); FLOW(I) is the flow on edge I.
!
!    Output, integer ( kind = 4 ) NODE_FLOW(NNODE); NODE_FLOW(I) is the flow
!    through node I.
!
!    Workspace, POINT(NNODE); POINT(I) is the first edge from node I.
!
!    Workspace, IMAP(NNODE), JMAP(N), pointer arrays.
!
  implicit none

  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) capac(nedge)
  logical finish
  integer ( kind = 4 ) flow(nedge)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icont
  integer ( kind = 4 ) iend
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) iflow
  integer ( kind = 4 ) imap(nnode)
  integer ( kind = 4 ) in
  integer ( kind = 4 ) inode(nedge)
  integer ( kind = 4 ) iout
  integer ( kind = 4 ) iparm
  integer ( kind = 4 ) isink
  integer ( kind = 4 ) isorce
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jcont
  integer ( kind = 4 ) jmap(nnode)
  integer ( kind = 4 ) jnode(nedge)
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) medge
  integer ( kind = 4 ) mincut(nnode)
  integer ( kind = 4 ) nodei
  integer ( kind = 4 ) nodej
  integer ( kind = 4 ) nodeu
  integer ( kind = 4 ) nodev
  integer ( kind = 4 ) nodew
  integer ( kind = 4 ) nodex
  integer ( kind = 4 ) nodey
  integer ( kind = 4 ) node_flow(nnode)
  integer ( kind = 4 ) point(nnode)
!
!  Initialize.
!
  icont = 0
  jcont = 0
  m1 = 0
  nodei = 0
  nodev = 0
  point(1:nnode) = 0
  iflow = 0
  flow(1:nedge) = 0

  do i = 1, nedge
    j = inode(i)
    if ( j == isorce ) then
      iflow = iflow + capac(i)
    end if
    point(j) = point(j) + 1
  end do

  node_flow(isorce) = iflow
  nodew = 1

  do i = 1, nnode
    j = point(i)
    point(i) = nodew
    imap(i) = nodew
    nodew = nodew + j
  end do

  finish = .false.
!
!  Sort the edges in lexicographical order.
!
40    continue

  iflag = 0

50    continue

  if ( iflag < 0 ) then

    if ( iflag /= -1 ) then

      if ( nodew < 0 ) then
        nodei = nodei + 1
      end if

      nodej = jcont
      jcont = nodei
      iflag = - 1

    else

      if ( nodew <= 0 ) then

        if ( 1 < icont ) then
          icont = icont - 1
          jcont = icont
          go to 60
        end if

        if ( m1 == 1 ) then
          iflag = 0
        else
          nodei = m1
          m1 = m1 - 1
          nodej = 1
          iflag = 1
        end if

      else

        iflag = 2

      end if

    end if

  else

    if ( 0 < iflag ) then

      if ( iflag <= 1) then
        jcont = icont
      end if

      go to 60

    else

      m1 = nedge
      icont = 1 + nedge / 2
      icont = icont - 1
      jcont = icont

60        continue

      nodei = 2 * jcont

      if ( nodei < m1 ) then

        nodej = nodei + 1
        iflag = - 2

      else

        if ( nodei == m1 ) then

          nodej = jcont
          jcont = nodei
          iflag = - 1

        else

          if ( 1 < icont ) then
            icont = icont - 1
            jcont = icont
            go to 60
          end if

          if ( m1 == 1 ) then
            iflag = 0
          else
            nodei = m1
            m1 = m1 - 1
            nodej = 1
            iflag = 1
          end if

        end if

      end if

    end if

  end if

  if ( iflag < 0 ) then

    nodew = inode(nodei) - inode(nodej)

    if ( nodew == 0 ) then
      nodew = jnode(nodei) - jnode(nodej)
    end if

    go to 50

  else if ( 0 < iflag ) then

    go to 70

  else if ( iflag == 0 ) then

    if ( finish ) then
      return
    end if

  end if
!
!  Set the cross references between edges.
!
  do i = 1, nedge
    nodev = jnode(i)
    inode(i) = imap(nodev)
    imap(nodev) = imap(nodev) + 1
  end do

90    continue

  iflag = 0

  do i = 1, nnode

    if ( i /= isorce ) then
      node_flow(i) = 0
    end if

    jmap(i) = nedge + 1

    if ( i < nnode ) then
      jmap(i) = point(i+1)
    end if

    mincut(i) = 0

  end do

  in = 0
  iout = 1
  imap(1) = isorce
  mincut(isorce) = - 1

110   continue

  in = in + 1

  if ( in <= iout ) then

    nodeu = imap(in)
    medge = jmap(nodeu) - 1
    iend = point(nodeu) - 1

120     continue

    iend = iend + 1

    if ( medge < iend ) then
      go to 110
    end if

    nodev = jnode(iend)
    iflow = capac(iend) - flow(iend)

    if ( mincut(nodev) /= 0 .or. iflow == 0 ) then
      go to 120
    end if

    if ( nodev /= isink ) then
      iout = iout + 1
      imap(iout) = nodev
    end if

    mincut(nodev) = - 1
    go to 120

  end if
!
!  Exit.
!
  if ( mincut(isink) == 0 ) then

    do i = 1, nnode
      mincut(i) = - mincut(i)
    end do

    do i = 1, nedge

      nodeu = jnode(inode(i))

      if ( flow(i) < 0 ) then
        node_flow(nodeu) = node_flow(nodeu) - flow(i)
      end if

      inode(i) = nodeu

    end do

    node_flow(isorce) = node_flow(isink)
    finish = .true.

    go to 40

  end if

  mincut(isink) = 1

150   continue

  in = in - 1

  if ( in /= 0 ) then

    nodeu = imap(in)
    nodei = point(nodeu) - 1
    nodej = jmap(nodeu) - 1

160     continue

    if ( nodei /= nodej ) then

      nodev = jnode(nodej)

      if ( mincut(nodev) <= 0 .or. capac(nodej) == flow(nodej) ) then

        nodej = nodej - 1
        go to 160

      end if

      jnode(nodej) = -nodev
      capac(nodej) = capac(nodej) - flow(nodej)
      flow(nodej) = 0
      nodei = nodei + 1

      if ( nodei < nodej ) then
        inode(inode(nodei)) = nodej
        inode(inode(nodej)) = nodei
        go to 70
      end if

    end if

    if ( point(nodeu) <= nodei ) then
      mincut(nodeu) = nodei
    end if

    go to 150

  end if

  nodex = 0

  do i = 1, iout

    if ( 0 < mincut(imap(i)) ) then
      nodex = nodex + 1
      imap(nodex) = imap(i)
    end if

  end do
!
!  Find a feasible flow.
!
  iflag = - 1
  nodey = 1

180   continue

  nodeu = imap(nodey)

  if ( node_flow(nodeu) <= 0 ) then

190     continue

    nodey = nodey + 1

    if ( nodey <= nodex ) then
      go to 180
    end if

    iparm = 0

200     continue

    nodey = nodey - 1

    if ( nodey /= 1 ) then

      nodeu = imap(nodey)

      if ( node_flow(nodeu) < 0 ) then
        go to 200
      end if
!
!  Accumulating flows.
!
      if ( node_flow(nodeu) == 0 ) then

        medge = nedge + 1

        if ( nodeu < nnode ) then
           medge = point(nodeu+1)
        end if

        iend = jmap(nodeu)
        jmap(nodeu) = medge

210         continue

        if ( iend == medge ) then
          go to 200
        end if

        j = inode(iend)
        iflow = flow(j)
        flow(j) = 0
        capac(j) = capac(j) - iflow
        flow(iend) = flow(iend) - iflow
        iend = iend + 1
        go to 210

      end if

      if ( mincut(nodeu) < point(nodeu) ) then
        iend = jmap(nodeu)
        go to 250
      end if

      iend = mincut(nodeu) + 1
      go to 230

    end if

    do i = 1, nedge

      nodev = - jnode(i)

      if ( 0 <= nodev ) then
        jnode(i) = nodev
        j = inode(i)
        capac(i) = capac(i) - flow(j)
        iflow = flow(i) - flow(j)
        flow(i) = iflow
        flow(j) = - iflow
      end if

    end do

    go to 90

  end if
!
!  An outgoing edge from a node is given maximum flow.
!
  iend = mincut(nodeu) + 1

230   continue

  iend = iend - 1

  if ( point(nodeu) <= iend ) then

    nodev = - jnode(iend)

    if ( node_flow(nodev) < 0 ) then
      go to 230
    end if

    iflow = capac(iend) - flow(iend)

    if ( node_flow(nodeu) < iflow ) then
      iflow = node_flow(nodeu)
    end if

    flow(iend) = flow(iend) + iflow
    node_flow(nodeu) = node_flow(nodeu) - iflow
    node_flow(nodev) = node_flow(nodev) + iflow
    iparm = 1
    nodei = inode(iend)
    nodej = jmap(nodev) - 1

    if ( nodei < nodej ) then
      inode(inode(nodei)) = nodej
      inode(inode(nodej)) = nodei
      go to 70
    end if

    if ( nodei == nodej ) then
      jmap(nodev) = nodej
    end if

240     continue

    if ( 0 < node_flow(nodeu) ) then
      go to 230
    end if

    if ( capac(iend) == flow(iend) ) then
      iend = iend - 1
    end if

  end if

  mincut(nodeu) = iend

  if ( iparm /= 0 ) then
    go to 190
  end if
!
!  Remove excess incoming flows from nodes.
!
  iend = jmap(nodeu)

250   continue

  j = inode(iend)
  iflow = flow(j)

  if ( node_flow(nodeu) < iflow ) then
    iflow = node_flow(nodeu)
  end if

  flow(j) = flow(j) - iflow

  node_flow(nodeu) = node_flow(nodeu) - iflow
  nodev = jnode(iend)
  node_flow(nodev) = node_flow(nodev) + iflow
  iend = iend + 1

  if ( 0 < node_flow(nodeu) ) then
    go to 250
  end if

  node_flow(nodeu) = - 1
  go to 200
!
!  Interchange two edges.
!
70    continue

  nodew = inode(nodei)
  inode(nodei) = inode(nodej)
  inode(nodej) = nodew

  iflow = capac(nodei)
  capac(nodei) = capac(nodej)
  capac(nodej) = iflow

  nodew = jnode(nodei)
  jnode(nodei) = jnode(nodej)
  jnode(nodej) = nodew

  iflow = flow(nodei)
  flow(nodei) = flow(nodej)
  flow(nodej) = iflow

  if ( 0 < iflag ) then
    go to 50
  else if ( iflag == 0 ) then
    go to 160
  end if

  jmap(nodev) = nodej
  go to 240

end
subroutine digraph_arc_print ( nedge, inode, jnode, title )

!*****************************************************************************80
!
!! DIGRAPH_ARC_PRINT prints out a digraph from an edge list.
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
!    Input, integer ( kind = 4 ) INODE(NEDGE), JNODE(NEDGE), the beginning and end
!    nodes of the edges.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) nedge

  integer ( kind = 4 ) i
  integer ( kind = 4 ) inode(nedge)
  integer ( kind = 4 ) jnode(nedge)
  character ( len = * ) title

  if ( len_trim ( title ) /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  do i = 1, nedge
    write ( *, '(2x,i8,4x,2i8)' ) i, inode(i), jnode(i)
  end do

  return
end
subroutine digraph_arc_prt_path ( nnode, nedge, npmax, k, maxpath, isorce, &
  isink, inode, jnode, arclen, dist )

!*****************************************************************************80
!
!! DIGRAPH_ARC_PRT_PATH outputs at most MAXPTH number of K shortest paths between two nodes.
!
!  Discussion:
!
!    Repeated nodes are allowed in the paths.  This routine may only be
!    called after DIGRAPH_ARC_KSHORT1 has computed the lengths of paths
!    between the nodes.
!
!  Modified:
!
!    23 July 2000
!
!  Author:
!
!    Hang Tong Lau
!
!  Reference:
!
!    Hang Tong Lau,
!    Algorithms on Graphs,
!    Tab Books, 1989,
!    QA166 L38.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) NEDGE, the number of edges.
!
!    Input, integer ( kind = 4 ) NPMAX, the maximum number of nodes in a possible
!    path which allows repeated nodes.
!
!    Input, integer ( kind = 4 ) K, the number of shortest paths.
!
!    Input, integer ( kind = 4 ) MAXPTH, the maximum number of paths to be generated.
!
!    Input, integer ( kind = 4 ) ISORCE, ISINK; the paths to be generated will go from
!    node ISORCE to ISINK.  ISORCE and ISINK may be equal.
!
!    Input, integer ( kind = 4 ) INODE(NEDGE), JNODE(NEDGE).  The I-th edge is directed from
!    INODE(I) to JNODE(I).  It is assumed that the array JNODE is already
!    sorted in nondecreasing order.
!
!    Input, real ( kind = 8 ) ARCLEN(NEDGE); ARCLEN(I) is the length of
!    the I-th edge.
!
!    Input, real ( kind = 8 ) DIST(NNODE,K), the output from KSHOT1.
!
!    Workspace, integer FIRARC(NNODE+1); contains the nodes on a path from
!    ISORCE to ISINK.
!
!    Workspace, integer NODINC(NEDGE); array of nodes I incident to node J in
!    the order of increasing J.
!
!    Workspace, integer LENINC(NEDGE); array of edge lengths corresponding
!    to NODINC.
!
!    Workspace, integer PATHND(NPMAX); contains the nodes on a path from
!    ISORCE to ISINK.
!
!    Workspace, integer POINT(NPMAX); POINT(I) is the position of node
!    PATHND(I).
!
!    Workspace, integer PLEN(NPMAX); PLEN(I) is the edge length from node
!    PATHND(I) to node PATHND(I-1).
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) nnode
  integer ( kind = 4 ) npmax

  real ( kind = 8 ) arclen(nedge)
  real ( kind = 8 ) dist(nnode,k)
  real ( kind = 8 ) dist1
  real ( kind = 8 ) dist2
  integer ( kind = 4 ) firarc(nnode+1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ids1
  integer ( kind = 4 ) if
  integer ( kind = 4 ) index
  integer ( kind = 4 ) init
  integer ( kind = 4 ) inode(nedge)
  integer ( kind = 4 ) ipath
  integer ( kind = 4 ) ipt1
  integer ( kind = 4 ) ipt2
  integer ( kind = 4 ) iptag
  integer ( kind = 4 ) isink
  integer ( kind = 4 ) isorce
  integer ( kind = 4 ) isub
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  real ( kind = 8 ) jlen
  integer ( kind = 4 ) jnode(nedge)
  real ( kind = 8 ) large
  real ( kind = 8 ) len
  real ( kind = 8 ) leninc(nedge)
  real ( kind = 8 ) lt
  integer ( kind = 4 ) maxpath
  integer ( kind = 4 ) nd
  integer ( kind = 4 ) nodeu
  integer ( kind = 4 ) nodev
  integer ( kind = 4 ) nodinc(nedge)
  integer ( kind = 4 ) nt
  integer ( kind = 4 ) numpath
  integer ( kind = 4 ) pathnd(npmax)
  integer ( kind = 4 ) point(npmax)
  real ( kind = 8 ) plen(npmax)

  init = 0
  index = 0
  large = 1.0D+00 + sum ( arclen(1:nedge) )

  do i = 1, nedge

    nodev = inode(i)
    nodeu = jnode(i)
    len = arclen(i)

    if ( nodeu /= index ) then

      j1 = index + 1
      j2 = nodeu - 1
      firarc(j1:j2) = 0
      firarc(nodeu) = init + 1
      index = nodeu

    end if

    init = init + 1
    nodinc(init) = nodev
    leninc(init) = len

  end do

  firarc(index+1) = init + 1
  pathnd(1:npmax) = 0
  point(1:npmax) = 0
  plen(1:npmax) = 0.0D+00

  if ( isorce == isink ) then
    ipath = 2
  else
    ipath = 1
  end if

  numpath = 0

  if ( large <= dist(isink,ipath) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DIGRAPH_ARC_PRT_PATH - Warning:'
    write ( *, '(a,i8,a,i8)' ) '  There is no path from ', isorce, ' to ', isink
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DIGRAPH_ARC_PRT_PATH:'
  write ( *, '(a,i8,a,i8)' ) '  Paths from ', isorce, ' to ', isink
  write ( *, '(a)' ) '  Index Length   ------Nodes-----'
  write ( *, '(a)' ) ' '

60    continue

  iptag = 1
  dist1 = dist(isink,ipath)

  if ( dist1 == large ) then
    return
  end if

  dist2 = dist1
  pathnd(1) = isink

70    continue

  ipt1 = 0

80    continue

  nt = pathnd(iptag)
  ids1 = firarc(nt)
  nd = nt
!
!  Find a value of ND for which FIRARC(ND+1) is not 0.
!
  do while ( firarc(nd+1) == 0 .and. nd < nnode )
    nd = nd + 1
  end do

  if = firarc(nd+1) - 1
  ipt2 = ids1 + ipt1

100   continue

  if ( ipt2 <= if ) then

    isub = nodinc(ipt2)
    jlen = leninc(ipt2)
    lt = dist1 - jlen
    j = 1

110     continue

    if ( lt < dist(isub,j) .or. k < j ) then
      ipt2 = ipt2 + 1
      go to 100
    end if

    if ( dist(isub,j) < lt ) then
      j = j + 1
      go to 110
    end if

    iptag = iptag + 1

    if ( npmax < iptag ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DIGRAPH_ARC_PRT_PATH - Warning!'
      write ( *, '(a,i8)' ) '  Number of arcs in path exceeds ', npmax
      write ( *, '(a)' ) '  Remedy:  Increase NPMAX.'
      return
    end if

    pathnd(iptag) = isub
    point(iptag) = ipt2 - ids1 + 1
    plen(iptag) = jlen
    dist1 = lt

    if ( dist1 /= 0.0D+00 ) then
      go to 70
    end if

    if ( isub /= isorce ) then
      go to 70
    end if

    numpath = numpath + 1

    write ( *, '(i4,g14.6)' ) numpath, dist2
    write ( *, '(4x,20i4)' ) ( pathnd(iptag-j+1), j = 1, iptag )

    if ( maxpath <= numpath ) then
      return
    end if

  end if

  ipt1 = point(iptag)

  pathnd(iptag) = 0
  dist1 = dist1 + plen(iptag)
  iptag = iptag - 1

  if ( 0 < iptag ) then
    go to 80
  end if

  ipath = ipath + 1

  if ( ipath <= k ) then
    go to 60
  end if

  return
end
subroutine digraph_arc_shtree ( nnode, nedge, iroot, inode, jnode, arclen, &
  dist, itree, jtree )

!*****************************************************************************80
!
!! DIGRAPH_ARC_SHTREE finds the shortest paths from IROOT to all other nodes.
!
!  Discussion:
!
!    Consider a digraph with lengths associated with the edges.
!    The edge lengths may be nonpositive, but no cycles of negative
!    lengths are present in the digraph.  The problem is to find the
!    shortest paths from a given node to every other node.
!
!  Method:
!
!    Let D(U,V) be the length of the edge from node U to node V.
!    The shortest path from a specified node S to every other node
!    can be found by the following label-correcting method:
!
!    STEP 1.  Let T be a tree consisting of the root S alone.  Set
!    L(S) = 0
!    L(I) = "Infinity" for all I not equal to S.
!
!    STEP 2.  Find an edge (U,V) such that L(U) + D(U,V) < L(V).
!    Set L(V) = L(U) + D(U,V).
!    Include the edge (U,V) in T, and remove any edge in T which
!    ends at V.
!
!    STEP 3.  Repeat Step 2 until LV <= L(U) + D(U,V) for all edges (U,V).
!    The tree T will now contain a shortest path from S to every other
!    node.
!
!    The computation of the method is O(N**3).
!
!  Modified:
!
!    22 July 2000
!
!  Author:
!
!    Hang Tong Lau
!
!  Reference:
!
!    Hang Tong Lau,
!    Algorithms on Graphs,
!    Tab Books, 1989,
!    QA166 L38.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) NEDGE, the number of edges.
!
!    Input, integer ( kind = 4 ) IROOT, the base node, from which distances to all
!    other nodes will be calculated.
!
!    Input, integer ( kind = 4 ) INODE(NEDGE), JNODE(NEDGE); INODE(I), JNODE(I)
!    are the beginning and ending nodes of edge I.
!
!    Input, real ( kind = 8 ) ARCLEN(NEDGE), the length of each edge.
!
!    Output, real ( kind = 8 ) DIST(NNODE).  DIST(I) is the length of
!    the shortest path from node IROOT to node I.
!
!    Output, integer ( kind = 4 ) ITREE(NNODE-1), JTREE(NNODE-1), the arcs in the
!    shortest path tree.
!
  implicit none

  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) arcfir(nnode+1)
  real ( kind = 8 ) arclen(nedge)
  real ( kind = 8 ) dist(nnode)
  integer ( kind = 4 ) fwdarc(nedge)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iedge
  integer ( kind = 4 ) index
  integer ( kind = 4 ) inode(nedge)
  integer ( kind = 4 ) iroot
  integer ( kind = 4 ) istart
  integer ( kind = 4 ) itree(nnode-1)
  integer ( kind = 4 ) iww
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jnode(nedge)
  integer ( kind = 4 ) jtree(nnode-1)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) last
  real ( kind = 8 ) lensum
  real ( kind = 8 ) lenu
  logical mark(nnode)
  integer ( kind = 4 ) nodeu
  integer ( kind = 4 ) nodev
  integer ( kind = 4 ) nodey
  integer ( kind = 4 ) origin(nedge)
  integer ( kind = 4 ) queue(nnode)
  integer ( kind = 4 ) tree_dad(nnode)
!
!  Set up the forward star representation of the digraph.
!
  k = 0
  do i = 1, nnode
    arcfir(i) = k + 1
    do j = 1, nedge
      if ( inode(j) == i ) then
        k = k + 1
        origin(k) = j
        fwdarc(k) = jnode(j)
      end if
    end do
  end do

  arcfir(nnode+1) = nedge + 1

  tree_dad(1:nnode) = 0
  mark(1:nnode) = .true.

  dist(1:nnode) = huge ( 1.0 )
  dist(iroot) = 0

  nodev = 1
  nodey = 1
  nodeu = iroot

50    continue

  lenu = dist(nodeu)
  istart = arcfir(nodeu)

  if ( istart /= 0 ) then

    index = nodeu + 1

    do

      last = arcfir(index) - 1

      if ( -1 < last ) then
        exit
      end if

      index = index + 1

    end do

    do i = istart, last

      iww = fwdarc(i)
      lensum = arclen(origin(i)) + lenu

      if ( lensum < dist(iww) ) then

        dist(iww) = lensum
        tree_dad(iww) = nodeu

        if ( mark(iww) ) then

          mark(iww) = .false.
          queue(nodey) = iww
          nodey = nodey + 1

          if ( nnode < nodey ) then
            nodey = 1
          end if

        end if

      end if

    end do

  end if

  if ( nodev /= nodey ) then

    nodeu = queue(nodev)
    mark(nodeu) = .true.
    nodev = nodev + 1

    if ( nnode < nodev ) then
      nodev = 1
    end if

    go to 50

  end if

  iedge = 0
  do i = 1, nnode
    if ( i /= iroot ) then
      iedge = iedge + 1
      itree(iedge) = tree_dad(i)
      jtree(iedge) = i
    end if
  end do

  return
end
subroutine digraph_arc_stcomp ( nnode, nedge, inode, jnode, numcomp, comp )

!*****************************************************************************80
!
!! DIGRAPH_ARC_STCOMP finds the strongly connected components of a digraph.
!
!  Discussion:
!
!    A strongly connected component of a digraph is a maximal set
!    of nodes in which there is a directed path from any one node in the
!    set to any other node in the set.  The problem is to find the strongly
!    connected components of a given digraph.
!
!  Method:
!
!    A depth-first search will be used to find the strongly connected
!    components of a digraph G of NNODE nodes and NEDGE edges.
!
!    STEP 1. A depth-first search of G is performed by selecting
!    one node S of G as a start node; S is marked "visited."  Each
!    unvisited node adjacent to S is searched in turn, using the
!    depth-first search recursively.  The search of S is completed
!    when all nodes that can be reached from S have been visited.
!    If some nodes remain unvisited, then an unvisited node is
!    arbitrarily selected as a new start node.  This process is
!    repeated until all nodes of G have been visited.
!
!    STEP 2. Let Rl, R2, ..., Rk be the roots in the order in which
!    the depth-first search of the nodes terminated.  Then, the
!    strongly connected component G, with root R, consists of all
!    descendants of R.  Furthermore, for each I, I = 2, 3, ..., K,
!    the strongly connected component G(I) with root R(I) consists of
!    these nodes which are descendants of R(I) but are in none of
!    G(1), G(2), ... , G(I-1).
!
!    The running time of the algorithm is O ( max ( NNODE, NEDGE ) ).
!
!  Modified:
!
!    28 July 2000
!
!  Author:
!
!    Hang Tong Lau
!
!  Reference:
!
!    Hang Tong Lau,
!    Algorithms on Graphs,
!    Tab Books, 1989,
!    QA166 L38.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) NEDGE, the number of edges.
!
!    Input, integer ( kind = 4 ) INODE(NEDGE), JNODE(NEDGE), the I-th edge of
!    the input digraph is directed from node INODE(I) to node JNODE(I).  The
!    digraph does not necessarily have to be connected.
!
!    Output, integer ( kind = 4 ) NUMCOMP, the number of strongly connected components
!    of the digraph
!
!    Output, integer ( kind = 4 ) COMP(NNODE), contains the component of each node.
!
  implicit none

  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) arcfir(nnode+1)
  integer ( kind = 4 ) comp(nnode)
  integer ( kind = 4 ) fwdarc(nedge)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) index
  integer ( kind = 4 ) index1
  integer ( kind = 4 ) index2
  integer ( kind = 4 ) index3
  integer ( kind = 4 ) inode(nedge)
  integer ( kind = 4 ) istack
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j3
  integer ( kind = 4 ) j4
  integer ( kind = 4 ) jnode(nedge)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) lowlnk(nnode)
  integer ( kind = 4 ) nodei
  integer ( kind = 4 ) nodesc(nnode)
  integer ( kind = 4 ) nodeu
  integer ( kind = 4 ) nodev
  integer ( kind = 4 ) number(nnode)
  integer ( kind = 4 ) numcomp
  integer ( kind = 4 ) stack1(nnode)
  integer ( kind = 4 ) stack2(nnode)
!
!  Set up the forward star representation of the digraph.
!
  arcfir(1) = 0
  k = 0

  do i = 1, nnode

    do j = 1, nedge
      if ( inode(j) == i ) then
        k = k + 1
        fwdarc(k) = jnode(j)
      end if
    end do

    arcfir(i+1) = k

  end do

  number(1:nnode) = 0

  numcomp = 0
  j1 = 0
  j2 = 0
  j3 = 0

  do nodei = 1, nnode

    if ( number(nodei) == 0 ) then

      istack = 1
      stack2(1) = nodei

40        continue

      nodev = stack2(istack)
      j1 = j1 + 1
      number(nodev) = j1
      lowlnk(nodev) = j1
      j2 = j2 + 1
      stack1(j2) = nodev
      index3 = arcfir(nodev+1)
      index = arcfir(nodev) + 1

50        continue

      if ( index <= index3 ) then

        nodeu = fwdarc(index)
        index2 = number(nodeu)

        if ( index2 == 0 ) then

          istack = istack + 1
          stack2(istack) = nodeu
          go to 40

        else

          if ( index2 < number(nodev) ) then
            do i = 1, j2
              if ( stack1(i) == nodeu ) then
                if ( index2 < lowlnk(nodev) ) then
                  lowlnk(nodev) = index2
                end if
                index = index + 1
                go to 50
              end if
            end do
          end if

        end if

        index = index + 1
        go to 50

      end if

      index1 = number(nodev)

      if ( lowlnk(nodev) == index1 ) then

        numcomp = numcomp + 1
        j4 = j2

        do while ( 0 < j4 )

          index2 = stack1(j4)

          if ( number(index2) < index1 ) then
            exit
          end if

          j3 = j3 + 1
          nodesc(j3) = index2
          j4 = j4 - 1

        end do

        nodesc(j3) = - nodesc(j3)
        j2 = j4

      end if

      if ( 1 < istack ) then

        nodeu = stack2(istack)
        istack = istack - 1
        nodev = stack2(istack)
        index1 = lowlnk(nodeu)

        if ( index1 < lowlnk(nodev) ) then
          lowlnk(nodev) = index1
        end if

        index = index + 1
        go to 50

      end if

    end if

  end do
!
!  Convert NODESC to component information.
!
  k = 1

  do i = 1, nnode

    j = nodesc(i)
    comp(abs(j)) = k

    if ( j < 0 ) then

      if ( k < numcomp ) then
        k = k + 1
      end if

    end if

  end do

  return
end
subroutine digraph_arc_to_star ( nnode, nedge, inode, jnode, arcfir, fwdarc )

!*****************************************************************************80
!
!! DIGRAPH_ARC_TO_STAR sets up the forward star representation of a digraph.
!
!  Modified:
!
!    04 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) NEDGE, the number of edges.
!
!    Input, integer ( kind = 4 ) INODE(NEDGE), JNODE(NEDGE); the I-th edge of the digraph
!    extends from node INODE(I) to JNODE(I).
!
!    Output, integer ( kind = 4 ) ARCFIR(NNODE+1); ARCFIR(I) is the number of the first
!    edge starting at node I in the forward star representation of the digraph.
!
!    Output, integer ( kind = 4 ) FWDARC(NEDGE); FWDARC(I) is the ending node of
!    the I-th edge in the forward star representation of the digraph.
!
  implicit none

  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) arcfir(nnode+1)
  integer ( kind = 4 ) fwdarc(nedge)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) inode(nedge)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jnode(nedge)
  integer ( kind = 4 ) k
!
!  Set up the forward star representation of the digraph.
!
  k = 0

  do i = 1, nnode

    arcfir(i) = k + 1

    do j = 1, nedge

      if ( inode(j) == i ) then
        k = k + 1
        fwdarc(k) = jnode(j)
      end if

    end do

  end do

  arcfir(nnode+1) = k + 1

  return
end
subroutine digraph_arc_weight_print ( nedge, inode, jnode, wnode, title )

!*****************************************************************************80
!
!! DIGRAPH_ARC_WEIGHT_PRINT prints out a weighted digraph from an edge list.
!
!  Modified:
!
!    23 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEDGE, the number of edges.
!
!    Input, integer ( kind = 4 ) INODE(NEDGE), JNODE(NEDGE), the beginning and end
!    nodes of the edges.
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

  if ( len_trim ( title ) /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  do i = 1, nedge
    write ( *, '(i8,4x,2i8,g14.6)' ) i, inode(i), jnode(i), wnode(i)
  end do

  return
end
subroutine digraph_dist_allpath ( nnode, dist, ndim, next )

!*****************************************************************************80
!
!! DIGRAPH_DIST_ALLPATH finds the shortest paths for all pairs of nodes in a digraph.
!
!  Discussion:
!
!    Find the shortest paths between all pairs of nodes in a digraph
!    with given edge lengths, assumming that there is no cycle of negative
!    length in the digraph.
!
!  Method:
!
!    Let L be a large number. Initially set D(I,J) equal
!    to the length of the edge from node I to node J in the given
!    digraph of a nodes.  It is assumed that D(I,I) = 0, for
!    I = 1, 2, ..., NNODE, and D(I,J) = L if the edge (I,J) is not in
!    the digraph.
!
!    STEP 1: Set K = 0.
!
!    STEP 2: Set K = K + 1.
!
!    STEP 3: For all I =/= K, such that D(I,K) =/= L, and all
!    J =/= K such that D(K,J) =/= L, compute:
!
!      D(I,J) = MIN ( D(I,J), D(I,K) + D(K,J) ).
!
!    STEP 4: If K = NNODE, then stop.  D(I,J) are the lengths of all
!    the shortest paths for all I and J.  Otherwise go to STEP 2.
!
!    The running time of the algorithm is O ( NNODE**3 ).
!
!  Modified:
!
!    20 July 1998
!
!  Author:
!
!    Hang Tong Lau
!
!  Reference:
!
!    Hang Tong Lau,
!    Algorithms on Graphs,
!    Tab Books, 1989,
!    QA166 L38.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input/output, real ( kind = 8 ) DIST(NDIM,NNODE);
!
!    On input, DIST(I,J) is the length of the edge directed from
!    node I to node J, with DIST(I,I) = 0 for all I, and
!    DIST(I,J) = HUGE(1.0) if there is no edge from node I to node J.
!
!    On output, DIST(I,J) is the length of the shortest path from node I
!    to node J.
!
!    Input, integer ( kind = 4 ) NDIM, the row dimension of DIST as specified
!    in the calling program.  NDIM must be at least NNODE.
!
!    Output, integer ( kind = 4 ) NEXT(NNODE,NNODE); NEXT(I,J) is the next-to-last node
!    in the shortest path from node I to node J.  This can be used to trace
!    the shortest path for every pair of nodes.
!
  implicit none

  integer ( kind = 4 ) ndim
  integer ( kind = 4 ) nnode

  real ( kind = 8 ) d
  real ( kind = 8 ) dist(ndim,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) next(nnode,nnode)

  do i = 1, nnode
    next(i,1:nnode) = i
  end do

  do i = 1, nnode
    do j = 1, nnode
      if ( dist(j,i) < huge ( d ) ) then
        do k = 1, nnode
          if ( dist(i,k) < huge ( d ) ) then
            d = dist(j,i) + dist(i,k)
            if ( d < dist(j,k) ) then
              dist(j,k) = d
              next(j,k) = next(i,k)
            end if
          end if
        end do
      end if
    end do
  end do

  return
end
subroutine digraph_dist_fmin ( nnode, distab, stack, key, little, j, level )

!*****************************************************************************80
!
!! DIGRAPH_DIST_FMIN stores values of K such that DISTAB(K) = LITTLE and DISTAB(J) = KEY.
!
!  Discussion:
!
!    The data is stored in STACK(1:LEVEL).
!
!    This routine is a utility routine used by routine SHORTP.
!
!  Modified:
!
!    22 July 2000
!
!  Author:
!
!    Hang Tong Lau
!
!  Reference:
!
!    Hang Tong Lau,
!    Algorithms on Graphs,
!    Tab Books, 1989,
!    QA166 L38.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
  implicit none

  integer ( kind = 4 ) nnode

  real ( kind = 8 ) distab(nnode)
  real ( kind = 8 ) distj
  integer ( kind = 4 ) j
  real ( kind = 8 ) key
  integer ( kind = 4 ) level
  real ( kind = 8 ) little
  integer ( kind = 4 ) stack(nnode)

  distj = distab(j)

  if ( key < distj ) then

    if ( distj < little ) then

      level = 1
      little = distj
      stack(level) = j

    else

      if ( distj == little ) then
        level = level + 1
        stack(level) = j
      end if

    end if

  end if

  return
end
subroutine digraph_dist_print ( dist, lda, nnode, title )

!*****************************************************************************80
!
!! DIGRAPH_DIST_PRINT prints a distance matrix.
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
!    Input, real ( kind = 8 ) DIST(LDA,NNODE), the distance matrix.
!    DIST(I,J) is the distance from node I to node J.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of DIST, which
!    must be at least NNODE.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) nnode

  real ( kind = 8 ) dist(lda,nnode)
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nrow
  character ( len = * ) title

  if ( len_trim ( title ) /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  ilo = 1
  ihi = nnode
  jlo = 1
  jhi = nnode
  ncol = nnode
  nrow = nnode

  call r8mat_print ( dist, ihi, ilo, jhi, jlo, lda, ncol, nrow )

  return
end
subroutine digraph_dist_short_ln ( nnode, dist, ndim, iroot, shdist )

!*****************************************************************************80
!
!! DIGRAPH_DIST_SHORT_LN finds the shortest paths from one node to all others.
!
!  Discussion:
!
!    Find the lengths of all shortest paths from a specified source
!    node to every other node in a given digraph with non-negative
!    edge lengths.
!
!  Method:
!
!    Let D(U,V) be the length of the edge from U to V in a complete
!    digraph of NNODE nodes.  S(I), the length of the shortest path from node 1
!    to node i, can be obtained as follows:
!
!    STEP 1. Set
!      P = 1
!      K = NNODE
!      S(1) = 0
!      L(I) = I, I = 1 to NNODE
!      S(I) = "infinity", I = 2 to NNODE.
!
!    STEP 2.  For I = 2, 3, ..., K do the following:
!      J = L(I)
!      S(J) = min ( S(J), S(P) + D(P,J) ).
!
!      If the value of S(J) is less than the current minimum, say S(Q),
!      during this execution of STEP 2, then set Q = J, T = I.
!
!    STEP 3.
!
!      P = Q,
!      L(T) = L(K),
!      K = K - 1.
!
!      If K = 1, than stop; otherwise, return to Step 2.
!
!    The algorithm requires NNODE**3/2 additions and NNODE**3 comparisoms.
!
!  Modified:
!
!    23 July 2000
!
!  Author:
!
!    Hang Tong Lau
!
!  Reference:
!
!    Hang Tong Lau,
!    Algorithms on Graphs,
!    Tab Books, 1989,
!    QA166 L38.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, real ( kind = 8 ) DIST(NDIM,NNODE);  DIST(I,J) is the nonnegative
!    length of the edge from node I to node J.  DIST(I,I) should be 0.0
!    for all I.
!
!    Input, integer ( kind = 4 ) NDIM, the row dunensdon of matrix DIST exactly as
!    specifled in the dimension statement of the calling program.
!
!    Input, integer ( kind = 4 ) IROOT, the node from which all distances will be
!    measured.
!
!    Output, real ( kind = 8 ) SHDIST(NNODE).  SHDIST(I) is the length of
!    the shortest path from IROOT to node I.
!
  implicit none

  integer ( kind = 4 ) nnode
  integer ( kind = 4 ) ndim

  real ( kind = 8 ) dist(ndim,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iroot
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) locate(nnode)
  integer ( kind = 4 ) minj
  real ( kind = 8 ) minlen
  integer ( kind = 4 ) minv
  integer ( kind = 4 ) nodeu
  integer ( kind = 4 ) nodev
  real ( kind = 8 ) shdist(nnode)
  real ( kind = 8 ) temp
!
!  Shift IROOT to position 1.
!
  if ( iroot /= 1 ) then
!
!  Interchange rows 1 and IROOT.
!
    do i = 1, nnode
      call r8_swap ( dist(1,i), dist(iroot,i) )
    end do
!
!  Interchange columns 1 and IROOT.
!
    do i = 1, nnode
      call r8_swap ( dist(i,1), dist(i,iroot) )
    end do

  end if

  nodeu = 1
  call i4vec_indicator ( nnode, locate )

  do i = 1, nnode
    shdist(i) = dist(nodeu,i)
  end do

  do i = 2, nnode

    k = nnode + 2 - i
    minlen = huge ( minlen )

    do j = 2, k

      nodev = locate(j)
      temp = shdist(nodeu) + dist(nodeu,nodev)

      if ( temp < shdist(nodev) ) then
        shdist(nodev) = temp
      end if

      if ( shdist(nodev) < minlen ) then
        minlen = shdist(nodev)
        minv = nodev
        minj = j
      end if

    end do

    nodeu = minv
    locate(minj) = locate(k)

  end do
!
!  Restore IROOT to its original position.
!
  if ( iroot /= 1 ) then

    shdist(1) = shdist(iroot)
    shdist(iroot) = 0.0D+00
!
!  Interchange rows 1 and IROOT
!
    do i = 1, nnode
      call r8_swap ( dist(1,i), dist(iroot,i) )
    end do
!
!  Interchange columns 1 and IROOT.
!
    do i = 1, nnode
      call r8_swap ( dist(i,1), dist(i,iroot) )
    end do

  end if

  return
end
subroutine digraph_dist_shortp ( nnode, dist, ndim, isorce, isink, numnod, &
  path, lnpath )

!*****************************************************************************80
!
!! DIGRAPH_DIST_SHORTP finds the shortest path from ISORCE to ISINK in a digraph.
!
!  Discussion:
!
!    Find the shortest path from a specified source node to another
!    specified destination node in a digraph with non-negative
!    edge lengths.
!
!  Method:
!
!    As a means of reducing the computational requirements to find the
!    shortest path from the source node S to the sink node T, the paths
!    of both directions from S and into T will be considered.  In general,
!    all paths out of S and into T as far as their adjacent connected nodes
!    are examined simultaneously.  The path which has so far covered the
!    least distance is extended.  This process is repeated until a path
!    is found out of S which has a node an it that already existed on a
!    path into T, or vice versa.  The complete path is then checked
!    to see if the shortest path is obtained.
!
!    If NNODE is the number of nodes in the digraph, than the processing time
!    of the entire algorithm is bounded by O ( NNODE**2 log NNODE**2 ).
!
!  Modified:
!
!    22 July 2000
!
!  Author:
!
!    Hang Tong Lau
!
!  Reference:
!
!    Hang Tong Lau,
!    Algorithms on Graphs,
!    Tab Books, 1989,
!    QA166 L38.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, real ( kind = 8 ) DIST(NDIM,NNODE); DIST(I,J) is the non-negative
!    length of the edge directed from node I to node J, with DIST(I,I) = 0.0
!    for all I.
!
!    Input, integer ( kind = 4 ) NDIM, the row dimension of matrix DIST exactly as
!    specified in the dimension statement of the calling program.
!
!    Input, integer ( kind = 4 ) ISORCE, ISINK, the specified source and destination node.
!
!    Output, integer ( kind = 4 ) NUMNOD, the number of nodes in the shortest path
!    from ISORCE to ISINK.
!
!    Output, integer ( kind = 4 ) PATH(NNODE); the nodes of the shortest path are stored
!    in entries 1 through NUMNOD of PATH.
!
!    Output, real ( kind = 8 ) LNPATH, is the length of the shortest path from
!    ISORCE to ISINK.
!
!    Workspace, integer POINT1(NNODE), POINT2(NNODE); the father of node I in
!    the current path from ISORCE to node I, and the son of node I in
!    the current path from node I to ISINK.
!
!    Workspace, real ( kind = 8 ) LEN1(NNODE), LEN2(NNODE), the lengths of
!    the current path from ISORCE to node I, and from node I to ISINK, as the
!    algorithm progresses.
!
!    Workspace, integer STACK1(NNODE), STACK2(NNODE).
!
  implicit none

  integer ( kind = 4 ) nnode
  integer ( kind = 4 ) ndim

  real ( kind = 8 ) dist(ndim,nnode)
  integer ( kind = 4 ) i
  real ( kind = 8 ) idist
  integer ( kind = 4 ) index
  integer ( kind = 4 ) ipt1
  integer ( kind = 4 ) ipt2
  integer ( kind = 4 ) isink
  integer ( kind = 4 ) isorce
  real ( kind = 8 ) key1
  real ( kind = 8 ) key2
  real ( kind = 8 ) len1(nnode)
  real ( kind = 8 ) len2(nnode)
  real ( kind = 8 ) lensum
  real ( kind = 8 ) lnpath
  real ( kind = 8 ) minln1
  real ( kind = 8 ) minln2
  real ( kind = 8 ) minln3
  integer ( kind = 4 ) mnode
  integer ( kind = 4 ) num1
  integer ( kind = 4 ) num2
  integer ( kind = 4 ) numnod
  integer ( kind = 4 ) path(nnode)
  integer ( kind = 4 ) point1(nnode)
  integer ( kind = 4 ) point2(nnode)
  integer ( kind = 4 ) stack1(nnode)
  integer ( kind = 4 ) stack2(nnode)

  key1 = 0
  key2 = 0

  len1(1:nnode) = dist(isorce,1:nnode)
  len2(1:nnode) = dist(1:nnode,isink)
  point1(1:nnode) = isorce
  point2(1:nnode) = isink
!
!  Find the initial values of MINLN1 and MINLN2 with corresponding INDEX
!  values for LEN1 and LEN2.
!
  ipt1 = 0
  ipt2 = 0
  minln1 = huge ( minln1 )
  minln2 = huge ( minln2 )

  do i = 1, nnode

    call digraph_dist_fmin ( nnode, len1, stack1, key1, minln1, i, ipt1 )

    call digraph_dist_fmin ( nnode, len2, stack2, key2, minln2, i, ipt2 )

  end do
!
!  Reset LEN1.
!
30    continue

  if ( minln1 <= minln2 ) then

    key1 = minln1

    do while ( 0 < ipt1 )

      index = stack1(ipt1)

      do i = 1, nnode

        idist = dist(index,i)
        lensum = minln1 + idist

        if ( lensum < len1(i) ) then
          len1(i) = lensum
          point1(i) = index
        end if

      end do

      ipt1 = ipt1 - 1

    end do
!
!  Find new MINLN1 and INDEX values for LEN1.
!
    minln1 = huge ( minln1 )
    ipt1 = 0

    do i = 1, nnode
      call digraph_dist_fmin ( nnode, len1, stack1, key1, minln1, i, ipt1 )
    end do
!
!  Reset LEN2.
!
  else

    key2 = minln2

    do while ( 0 < ipt2 )

      index = stack2(ipt2)

      do i = 1, nnode

        idist = dist(i,index)
        lensum = minln2 + idist

        if ( lensum < len2(i) ) then
          len2(i) = lensum
          point2(i) = index
        end if

       end do

       ipt2 = ipt2 - 1

      end do
!
!  Find new MINLN2 and INDEX values for LEN2.
!
      minln2 = huge ( minln2 )
      ipt2 = 0

      do i = 1, nnode
        call digraph_dist_fmin ( nnode, len2, stack2, key2, minln2, i, ipt2 )
      end do

    end if
!
!  Compute convergence criterion.
!
     minln3 = huge ( minln3 )

     do i = 1, nnode
       lensum = len1(i) + len2(i)
       if ( lensum < minln3 ) then
         minln3 = lensum
         mnode = i
       end if
     end do

    if ( minln1 + minln2 < minln3 ) then
      go to 30
    end if
!
!  Two ends of a shortest path meet in node mnode; unravel the path.
!
   num1 = mnode
   path(nnode) = mnode

   if ( mnode /= isorce ) then

    numnod = nnode - 1

    do

      num2 = point1(num1)

      if ( num2 == isorce ) then
        exit
      end if

      num1 = num2
      path(numnod) = num2
      numnod = numnod - 1

    end do

  else

    numnod = nnode

  end if

  path(1) = isorce
  num1 = numnod + 1
  numnod = 2

  do while ( num1 <= nnode )
    path(numnod) = path(num1)
    numnod = numnod + 1
    num1 = num1 + 1
  end do

  if ( mnode /= isink ) then

    num1 = mnode

    do

      num2 = point2(num1)

      if ( num2 == isink ) then
        exit
      end if

      num1 = num2
      path(numnod) = num2
      numnod = numnod + 1

    end do

    path(numnod) = isink

  end if

  lnpath = len1(mnode) + len2(mnode)

  return
end
subroutine graph_arc_clique ( nnode, nedge, inode, jnode, cliq )

!*****************************************************************************80
!
!! GRAPH_ARC_CLIQUE finds all the cliques of a graph.
!
!  Discussion:
!
!    Consider a graph G.  An independent set is a subset of the nodes of G
!    such that no two nodes of the set are adjacent in G.  An independent
!    set is maximal if there is no strictly larger independent set that
!    contains it.  A clique is a subset of nodes of G in which every two
!    nodes in the set are adjacent in G.
!
!    The problem is to find all the maximal independent sets
!    and cliques of a given undirected graph.
!
!  Method:
!
!    It is clear that a subset of nodes C of a graph G is a maximal
!    independent set if and only if C is a clique in the complement
!    of G.  Thus, any algorithm which finds the maximal independent
!    sets of a graph can also be used to find its cliques, and vice versa.
!
!    The algorithm for finding all maximal independent sets is essentially
!    an enumerative tree search.  Let P(J) be an independent set at stage J,
!    and Q(J) be the largest set of nodes such that any node from Q(J)
!    added to P(J) will produce an independent set P(J+1).  The set Q(J)
!    can be partitioned into two disjoint sets S(J) and T(J), where S(J) is
!    the set of all nodes which have been used in the search to augment P(J),
!    and T(J) is the set of all nodes which have not been used.  Let E(U)
!    be the set of nodes V such that node U and node V are adjacent in G.
!    The algorithm can be described as follows.
!
!    STEP 1. Set P(0) = empty, S(0) = empty, T(0) = set of
!    nodes of G, and J = 0.
!
!    STEP 2. Perform a forward branch in the tree search by
!    choosing any node X(J) from T(J), and create three new sets:
!
!      P(J+1) = P(J) + {X(J)}
!      S(J+1) = S(J) - E(X(J)),
!      T(J+1) = T(J) - E(X(J)) - {X(J)}.
!
!    Set J = J + 1.
!
!    STEP 3. If there exists a node Y in S(J) such that E(Y) is disjoint
!    from T(J), then go to Step 5.
!
!    STEP 4. If both S(J) and T(J) are empty, then output the
!    maximal independent set P(J) and go to Step 5.
!    If S(J) is nonempty and T(J) is empty, then go to Step 5;
!    otherwise, return to Step 2.
!
!    STEP 5. Backtrack by setting J = J - 1.  Remove node X(J) from P(J+1)
!    to produce P(J).  Remove node X(J) from T(J) and add it to S(J).
!    If J = 0 and T(0) is empty, then stop; otherwise, return
!    to Step 3.
!
!  Modified:
!
!    23 July 2000
!
!  Author:
!
!    Hang Tong Lau
!
!  Reference:
!
!    Hang Tong Lau,
!    Algorithms on Graphs,
!    Tab Books, 1989,
!    QA166 L38.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) NEDGE, the number of edges.
!
!    Input, integer ( kind = 4 ) INODE(NEDGE), JNODE(NEDGE); INODE(I) and JNODE(I) are the
!    two end nodes of the I-th edge in the input graph, which is not
!    necessarily connected.
!
!    Output, integer ( kind = 4 ) CLIQ(NNODE); each time the output WRITE statement in
!    CLIQUE is executed, a clique CLIQ(I), I = 1, 2,. .., CLIQ_NUM, is
!    generated, where CLIQ_NUM is a local variable in CLIQUE.
!
!    Workspace, integer FWDARC(2*NEDGE); FWDARC(I) is the ending node of
!    the I-th edge in the forward star representation of the graph.
!
!    Workspace, integer ARCFIR(NNODE+1); ARCFIR(I) is the number of the first
!    edge starting at node I in the forward star representation of the graph.
!
!    Workspace, integer CAND1(NNODE+1); positions of candidates which have not
!    been used.
!
!    Workspace, integer CAND2(NNODE+1); positions of candidates which have
!    been used.
!
!    Workspace, integer STACK(NNODE+1,NNODE+1); stack of nodes at different
!    stages.
!
!    Workspace, integer NDPS1(NNODE); array of nodes at each stage.
!
!    Workspace, integer NDPS2(NNODE); array of nodes at each stage.
!
  implicit none

  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) arcfir(nnode+1)
  integer ( kind = 4 ) cand1(nnode+1)
  integer ( kind = 4 ) cand2(nnode+1)
  integer ( kind = 4 ) cliq(nnode)
  integer ( kind = 4 ) cliq_num
  integer ( kind = 4 ) fwdarc(2*nedge)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indexv
  integer ( kind = 4 ) inode(nedge)
  integer ( kind = 4 ) ismall
  integer ( kind = 4 ) isub1
  integer ( kind = 4 ) isub2
  integer ( kind = 4 ) isum
  integer ( kind = 4 ) itemp
  integer ( kind = 4 ) iup
  integer ( kind = 4 ) jnode(nedge)
  logical join
  integer ( kind = 4 ) k
  integer ( kind = 4 ) level
  integer ( kind = 4 ) level1
  integer ( kind = 4 ) low
  integer ( kind = 4 ) ndps1(nnode)
  integer ( kind = 4 ) ndps2(nnode)
  integer ( kind = 4 ) nodeu
  integer ( kind = 4 ) nodev
  integer ( kind = 4 ) nodew
  integer ( kind = 4 ) stack(nnode+1,nnode+1)
!
!  Set up the forward star representation of the graph.
!
  call graph_arc_to_star ( nnode, nedge, inode, jnode, arcfir, fwdarc )

  level = 1
  level1 = 2
  do i = 1, nnode
    stack(level,i) = i
  end do

  cliq_num = 0
  cand2(level) = 0
  cand1(level) = nnode

40    continue

  ismall = cand1(level)
  nodeu = 0
  ndps1(level) = 0

50    continue

  nodeu = nodeu + 1

  if ( nodeu <= cand1(level) .and. ismall /= 0 ) then

    isub1 = stack(level,nodeu)
    isum = 0
    nodev = cand2(level)

60      continue

    nodev = nodev + 1

    if ( nodev <= cand1(level) .and. isum < ismall ) then

      itemp = stack(level,nodev)

      if ( itemp == isub1 ) then

        join = .TRUE.

      else

        join = .FALSE.
        low = arcfir(itemp)
        iup = arcfir(itemp + 1) - 1

        do k = low, iup
          if ( fwdarc(k) == isub1 ) then
            join = .TRUE.
            exit
          end if
        end do

      end if
!
!  Store the potential candidate.
!
    if ( .not. join ) then
      isum = isum + 1
      indexv = nodev
    end if

    go to 60

  end if

  if ( isum < ismall ) then

    ndps2(level) = isub1
    ismall = isum

    if ( nodeu <= cand2(level) ) then
      nodew = indexv
    else
      nodew = nodeu
      ndps1(level) = 1
    end if

  end if

  go to 50

  end if
!
!  Backtrack.
!
  ndps1(level) = ndps1(level) + ismall

90    continue

    if ( 0 < ndps1(level) ) then

      isub1 = stack(level,nodew)
      stack(level,nodew) = stack(level,cand2(level) + 1)
      stack(level,cand2(level) + 1) = isub1
      isub2 = isub1
      nodeu = 0
      cand2(level1) = 0

100       continue

      nodeu = nodeu + 1

      if ( nodeu <= cand2(level) ) then

        itemp = stack(level,nodeu)

        if ( itemp == isub2 ) then
          join = .TRUE.
        else
          join = .FALSE.

          low = arcfir(itemp)
          iup = arcfir(itemp+1) - 1

          do k = low, iup
            if ( fwdarc(k) == isub2 ) then
              join = .TRUE.
              exit
            end if
          end do

        end if

        if ( join ) then
          cand2(level1) = cand2(level1) + 1
          stack(level1,cand2(level1)) = stack(level,nodeu)
        end if

        go to 100

      end if

      cand1(level1) = cand2(level1)
      nodeu = cand2(level) + 1

 130      continue

      nodeu = nodeu + 1

      if ( nodeu <= cand1(level) ) then

        itemp = stack(level,nodeu)

        if ( itemp == isub2 ) then

          join = .TRUE.

        else

          join = .FALSE.
          low = arcfir(itemp)
          iup = arcfir(itemp+1) - 1

          do k = low, iup
            if ( fwdarc(k) == isub2 ) then
              join = .TRUE.
              exit
            end if
          end do

        end if

        if ( join ) then
          cand1(level1) = cand1(level1) + 1
          stack(level1,cand1(level1)) = stack(level,nodeu)
        end if

        go to 130

      end if

      cliq_num = cliq_num + 1
      cliq(cliq_num) = isub2

      if ( cand1(level1) == 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(4x,20i3)' ) cliq(1:cliq_num)
      else if ( cand2(level1) < cand1(level1) ) then
        level = level + 1
        level1 = level1 + 1
        go to 40
      end if

170   continue

      cliq_num = cliq_num - 1
      cand2(level) = cand2(level) + 1

      if ( 1 < ndps1(level) ) then

        nodew = cand2(level)
!
!  Look for candidate.
!
180       continue

      nodew = nodew + 1
      itemp = stack(level,nodew)

      if ( itemp == ndps2(level) ) then
        go to 180
      end if

      low = arcfir(itemp)
      iup = arcfir(itemp+1) - 1

      do k = low, iup
        if ( fwdarc(k) == ndps2(level) ) then
          go to 180
        end if
      end do

    end if

    ndps1(level) = ndps1(level) - 1
    go to 90

  end if
!
!  This jump to line 170 is not in good taste...
!
  if ( 1 < level ) then
    level = level - 1
    level1 = level1 - 1
    go to 170
  end if

  return
end
subroutine graph_arc_color_number ( nnode, nedge, inode, jnode, ncolor, color )

!*****************************************************************************80
!
!! GRAPH_ARC_COLOR_NUMBER finds the chromatic number of a graph.
!
!  Discussion:
!
!    A coloring of an undirected graph G is an assignment of colors to the
!    nodes of G such that no adjacent nodes have the same color.  A graph
!    is K-colorable if there is a coloring of G using K colors.  The chromatic
!    number of G is the minimum K for which G is K-colorable.  The problem is
!    to find the chromatic number of a given graph.
!
!    The chromatic number is always between 1 and N.
!
!  Method:
!
!    The coloring is done by a simple implicit enumeration tree
!    search method.  Initially, node 1 is assigned color 1, and the
!    remaining nodes are colored sequentially so that node I is
!    colored with the lowest-numbered color which has not been
!    used so far to color any nodes adjacent to it.
!
!    Let P be the number of colors required by this feasible
!    coloring.  Attempt to generate a feasible coloring using Q < P
!    colors.  To accomplish this, all nodes colored with P must be
!    recolored.  Thus, a backtrack step can be taken up to node U,
!    where node U + 1 is the lowest index assigned color P.
!
!    Attempt to color node U with its smallest feasible alternate
!    color greater than its current color.  If there is no such
!    alternative color which is smaller than P, then backtrack to
!    node U - 1.  Otherwise, recolor node U and proceed forward,
!    sequentially recoloring all nodes U + 1, U + 2, .. ., with
!    the smallest feasible color until either node N is colored or
!    some node V is reached which requires color P.  In the former
!    case, an improved coloring using Q colors has been found;
!    in this case, backtrack and attempt to find a better coloring
!    using less than Q colors.  In the latter case, backtrack from
!    node V and proceed forward as before.
!
!    The algorithm terminates when backtracking reaches node 1.
!
!  Modified:
!
!    27 July 2000
!
!  Author:
!
!    Hang Tong Lau
!
!  Reference:
!
!    Hang Tong Lau,
!    Algorithms on Graphs,
!    Tab Books, 1989,
!    QA166 L38.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) NEDGE, the number of edges.
!
!    Input, integer ( kind = 4 ) INODE(NEDGE), JNODE(NEDGE); the I-th edge of the graph
!    extends from node INODE(I) to JNODE(I).
!
!    Output, integer ( kind = 4 ) NCOLOR, the number of colors used.
!
!    Output, integer ( kind = 4 ) COLOR(NNODE), the colors assigned to each node.
!
  implicit none

  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) arcfir(nnode+1)
  integer ( kind = 4 ) availc(nnode,nnode)
  integer ( kind = 4 ) cmax1(nnode)
  integer ( kind = 4 ) cmax2(nnode)
  integer ( kind = 4 ) color(nnode)
  integer ( kind = 4 ) color1(nnode)
  integer ( kind = 4 ) color2(nnode)
  integer ( kind = 4 ) degree(nnode)
  integer ( kind = 4 ) fwdarc(2*nedge)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) index
  integer ( kind = 4 ) inode(nedge)
  integer ( kind = 4 ) ipaint
  integer ( kind = 4 ) iup
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jnode(nedge)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) loop
  integer ( kind = 4 ) low
  logical more
  integer ( kind = 4 ) ncolor
  integer ( kind = 4 ) newc
  integer ( kind = 4 ) nod
  integer ( kind = 4 ) nodek

  cmax2(1:nnode) = 0
!
!  Set up the forward star representation of the graph.
!
  call graph_arc_to_star ( nnode, nedge, inode, jnode, arcfir, fwdarc )
!
!  Compute the degree of each node.
!
  call graph_arc_degree ( nnode, nedge, inode, jnode, degree )

  color2(1:nnode) = degree(1:nnode) + 1

  do i = 1, nnode

    if ( i < color2(i) ) then
      color2(i) = i
    end if

    cmax1(i) = color2(i)
    k = color2(i)
    availc(i,1:k) = nnode
    availc(i,k+1:nnode) = 0

  end do

  nod = 1
!
!  Color node NOD.
!
  newc = 1
  ncolor = nnode
  ipaint = 0
  more = .true.

  do

    if ( more ) then

      index = cmax1(nod)

      if ( ipaint+1 < index ) then
        index = ipaint + 1
      end if

      do while ( availc(nod,newc) < nod .and. newc <= index )
        newc = newc + 1
      end do
!
!  Node NOD has the color NEWC.
!
      if ( newc == index+1 ) then

        more = .false.

      else
!
!  A new coloring is found.
!
        if ( nod == nnode ) then

          color1(nod) = newc
          color(1:nnode) = color1(1:nnode)

          if ( ipaint < newc ) then
            ipaint = ipaint + 1
          end if

          ncolor = ipaint
!
!  Backtrack to the first node of color NCOLOR.
!
          if ( 2 < ncolor ) then

            index = 1

            do while ( color(index) /= ncolor )
              index = index + 1
            end do

            j = nnode

            do while ( index <= j )

              nod = nod - 1
              newc = color1(nod)
              ipaint = cmax2(nod)
              low = arcfir(nod)
              iup = arcfir(nod+1) - 1

              do k = low, iup

                nodek = fwdarc(k)

                if ( nod < nodek ) then
                  if ( availc(nodek,newc) == nod ) then
                    availc(nodek,newc) = nnode
                    color2(nodek) = color2(nodek) + 1
                  end if
                end if

              end do

              newc = newc + 1
              more = .false.
              j = j - 1

            end do

            ipaint = ncolor - 1

            do i = 1, nnode

              loop = cmax1(i)

              if ( ipaint < loop ) then

                k = ipaint + 1

                do j = k, loop
                  if ( availc(i,j) == nnode ) then
                    color2(i) = color2(i) - 1
                  end if
                end do

                cmax1(i) = ipaint

              end if

            end do

          end if

        else
!
!  NOD is less than NNODE.
!
          low = arcfir(nod)
          iup = arcfir(nod+1) - 1

          if ( low <= iup ) then

            k = low

            do while ( k <= iup .and. more )

              nodek = fwdarc(k)

              if (  nod < nodek ) then
                more = .not. ( ( color2(nodek) == 1 ) .and. &
                 ( nod <= availc(nodek,newc) ) )
              end if

              k = k + 1

            end do

          end if

          if ( more ) then

            color1(nod) = newc
            cmax2(nod) = ipaint

            if ( ipaint < newc ) then
              ipaint = ipaint + 1
            end if

            low = arcfir(nod)
            iup = arcfir(nod+1) - 1

            do k = low, iup

              nodek = fwdarc(k)

              if ( nod < nodek ) then
                if ( nod <= availc(nodek,newc) ) then
                  availc(nodek,newc) = nod
                  color2(nodek) = color2(nodek) - 1
                end if
              end if

            end do

            nod = nod + 1
            newc = 1

          else

            newc = newc + 1

          end if

        end if

      end if

    else

      more = .true.

      if ( cmax1(nod) < newc .or. ipaint+1 < newc ) then

        nod = nod - 1
        newc = color1(nod)
        ipaint = cmax2(nod)

        low = arcfir(nod)
        iup = arcfir(nod+1) - 1

        do k = low, iup

          nodek = fwdarc(k)

          if ( nod < nodek ) then
            if ( availc(nodek,newc) == nod ) then
              availc(nodek,newc) = nnode
              color2(nodek) = color2(nodek) + 1
            end if
          end if

        end do

        newc = newc + 1
        more = .false.

      end if

    end if

    if ( nod == 1 .or. ncolor == 2 ) then
      exit
    end if

  end do

  return
end
subroutine graph_arc_color_poly ( nnode, nedge, inode, jnode, cpoly1, cpoly2, &
  cpoly3 )

!*****************************************************************************80
!
!! GRAPH_ARC_COLOR_POLY computes the chromatic polynomial of a graph.
!
!  Discussion:
!
!    A (node-) coloring of a graph is an assignment of a color to each node,
!    in such a way that nodes sharing an edge have different colors.
!    Letting X be the number of colors used in a particular coloring, then
!    for any X, the number of distinct colorings of a graph can be expressed
!    as a polynomial P(X).
!
!    If a graph has at least one node, it can't have a coloring with no
!    colors, so P must have a factor of X.
!
!    If a graph has any edges at all, it can't have a coloring with 1 color,
!    so P must have a factor of (X-1).
!
!    If the graph has NNODE nodes, then P is a polynomial of degree NNODE-1
!    or less.
!
!    The routine returns the polynomial in three formats:
!
!    * Alternating Sign Standard Form:
!
!      P(X) = A(N) * X^N - A(N-1) * X^(N-1) + ... + (-1)^(N-1) * A(1) * X
!
!    * Tutte or Tree Form:
!
!      P(X) = Sum ( 1 <= I <= N ) (-1)^(N-I) * A(I) * X * ( X - 1 )^(I-1)
!
!    * Factorial Form:
!
!      P(X) = Sum ( 1 <= I <= N ) A(I) * [X}_I
!      where {X}_I = X * ( X - 1 ) * ... * ( X + 1 - I )
!
!  Modified:
!
!    27 July 2000
!
!  Author:
!
!    Hang Tong Lau
!
!  Reference:
!
!    Hang Tong Lau,
!    Algorithms on Graphs,
!    Tab Books, 1989,
!    QA166 L38.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) NEDGE, the number of edges.
!
!    Input, integer ( kind = 4 ) INODE(NEDGE), JNODE(NEDGE); the I-th edge
!    starts at node INODE(I) and ends at node JNODE(I).
!
!    Output, integer ( kind = 4 ) CPOLY1(NNODE), CPOLY2(NNODE), CPOLY3(NNODE), the
!    coefficients of the chromatic polynomial, in Alternating Sign Standard
!    Form, Tutte Form, and Factorial Form.
!
  implicit none

  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) cpoly1(nnode)
  integer ( kind = 4 ) cpoly2(nnode)
  integer ( kind = 4 ) cpoly3(nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ilast
  integer ( kind = 4 ) incr
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) inode(nedge)
  integer ( kind = 4 ) istack((nnode*(nedge+nedge-nnode+1))/2)
  integer ( kind = 4 ) isub1
  integer ( kind = 4 ) isub2
  integer ( kind = 4 ) itop
  integer ( kind = 4 ) ivertx
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) iy
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jlast
  integer ( kind = 4 ) jnode(nedge)
  integer ( kind = 4 ) jstack((nnode*(nedge+nedge-nnode+1))/2)
  integer ( kind = 4 ) jsub1
  integer ( kind = 4 ) jsub2
  integer ( kind = 4 ) jvertx
  integer ( kind = 4 ) loop
  integer ( kind = 4 ) maxmn
  integer ( kind = 4 ) mm
  integer ( kind = 4 ) ncomp
  integer ( kind = 4 ) nn
  integer ( kind = 4 ) nodeu
  integer ( kind = 4 ) nodev
  integer ( kind = 4 ) nodew
  integer ( kind = 4 ) nodex
  integer ( kind = 4 ) nodey
  logical visit

  itop = 0
  cpoly2(1:nnode) = 0
  mm = nedge
  nn = nnode
!
!  Find a spanning tree, and work on the corresponding component.
!
  do

    maxmn = mm + 1

    if ( mm < nn ) then
      maxmn = nn + 1
    end if

    do i = 1, nn
      cpoly3(i) = - i
    end do

    do i = 1, mm
      j = inode(i)
      inode(i) = cpoly3(j)
      cpoly3(j) = - maxmn - i
      j = jnode(i)
      jnode(i) = cpoly3(j)
      cpoly3(j) = - maxmn - maxmn - i
    end do

    ncomp = 0
    indx = 0
    nodew = 0

50  continue

    nodew = nodew + 1

    if ( nodew <= nn ) then

      nodev = cpoly3(nodew)

      if ( 0 < nodev ) then
        go to 50
      end if

      ncomp = ncomp + 1
      cpoly3(nodew) = ncomp

      if ( -nodew <= nodev ) then
        go to 50
      end if

      nodeu = nodew
      nodex = - nodev
      visit = .true.
      isub2 = - nodev / maxmn
      jsub2 = - nodev - isub2 * maxmn

60    continue

      if ( isub2 == 1 ) then
        nodev = inode(jsub2)
      else
        nodev = jnode(jsub2)
      end if

      if ( 0 < nodev ) then
        if ( nodev <= maxmn ) then
          if ( 0 <= cpoly3(nodev) ) then
            nodey = jsub2
          end if
        end if
      end if

      if ( 0 <= nodev ) then

        if ( ix == 1 ) then
          nodex = abs ( inode(iy) )
          inode(iy) = nodeu
        else
          nodex = abs ( jnode(iy) )
          jnode(iy) = nodeu
        end if

        ix = nodex / maxmn
        iy = nodex - ix * maxmn

        if ( ix == 0 ) then

          do

            if ( nodey /= 0 ) then
              indx = indx + 1
              cpoly3(nodeu) = inode(nodey) + jnode(nodey) - nodeu
              inode(nodey) = - indx
              jnode(nodey) = nodeu
            end if

            nodeu = iy

            if ( iy <= 0 ) then
              go to 50
            end if

            nodex = - cpoly3(nodeu)
            ix = nodex / maxmn
            iy = nodex - ix * maxmn

            if ( ix /= 0 ) then
              exit
            end if

          end do

        end if

        isub2 = 3 - ix
        jsub2 = iy
        go to 60

      end if

      if ( nodev < - maxmn ) then

        if ( isub2 == 1 ) then
          inode(jsub2) = - nodev
        else
          jnode(jsub2) = - nodev
        end if

        isub2 = - nodev / maxmn
        jsub2 = - nodev - isub2 * maxmn

      else

        nodev = - nodev

        if ( isub2 == 1 ) then
          inode(jsub2) = 0
        else
          jnode(jsub2) = 0
        end if

        if ( visit ) then

          isub1 = isub2
          jsub1 = jsub2
          nodey = 0
          visit = .false.

        else

          if ( isub1 == 1 ) then
           inode(jsub1) = nodev
          else
           jnode(jsub1) = nodev
          end if

          isub1 = isub2
          jsub1 = jsub2

          if ( ix == 1 ) then
            nodex = abs ( inode(iy) )
            inode(iy) = nodeu
          else
            nodex = abs ( jnode(iy) )
            jnode(iy) = nodeu
          end if

        end if

        ix = nodex / maxmn
        iy = nodex - ix * maxmn

        if ( ix == 0 ) then

          do

            if ( nodey /= 0 ) then
              indx = indx + 1
              cpoly3(nodeu) = inode(nodey) + jnode(nodey) - nodeu
              inode(nodey) = - indx
              jnode(nodey) = nodeu
            end if

            nodeu = iy

            if ( iy <= 0 ) then
              go to 50
            end if

            nodex = - cpoly3(nodeu)
            ix = nodex / maxmn
            iy = nodex - ix * maxmn

            if ( ix /= 0 ) then
              exit
            end if

          end do

        end if

        isub2 = 3 - ix
        jsub2 = iy

      end if

      go to 60

    end if

    do i = 1, mm

      do

        nodey = - inode(i)

        if ( nodey < 0 ) then
          exit
        end if

        nodex = jnode(i)
        jnode(i) = jnode(nodey)
        jnode(nodey) = nodex
        inode(i) = inode(nodey)
        inode(nodey) = cpoly3(jnode(nodey))

      end do

    end do

    do i = 1, indx
      cpoly3(jnode(i)) = cpoly3(inode(i))
    end do
!
!  If NCOMP is not equal to 1, then the graph is not connected.
!
    if ( ncomp /= 1 .or. ( mm < nn .and. itop == 0 ) ) then

      if ( mm < nn .and. itop == 0 ) then
        cpoly2(nn) = cpoly2(nn) + 1
      end if

      cpoly1(1:nnode) = cpoly2(1:nnode)

      do i = 1, nnode
        cpoly3(i) = cpoly2(i) * ( 1 - 2 * mod ( nnode - i, 2 ) )
      end do

      do i = 1, nnode
        jvertx = 0
        do j = i, nnode
          jvertx = cpoly1(nnode+i-j) + jvertx
          cpoly1(nnode+i-j) = jvertx
        end do
      end do

      incr = 0

      do i = 1, nnode
        jvertx = 0
        do j = i, nnode
          jvertx = cpoly3(nnode+i-j) + incr * jvertx
          cpoly3(nnode+i-j) = jvertx
        end do
        incr = incr + 1
      end do

      return

    end if

    if ( mm < nn ) then

      cpoly2(nn) = cpoly2(nn) + 1

      nn = istack(itop)
      mm = jstack(itop)
      itop = itop - mm - 1

      do i = 1, mm
        inode(i) = istack(itop+i)
        jnode(i) = jstack(itop+i)
      end do

      if ( mm == nn ) then
        cpoly2(nn) = cpoly2(nn) + 1
      else
        itop = itop + mm
        istack(itop) = nn
        jstack(itop) = mm - 1
      end if

    else

      if ( mm == nn ) then
        cpoly2(nn) = cpoly2(nn) + 1
      else
        do i = 1, mm
          itop = itop + 1
          istack(itop) = inode(i)
          jstack(itop) = jnode(i)
        end do
        istack(itop) = nn
        jstack(itop) = mm - 1
      end if

    end if

    cpoly1(1:nnode) = 0

    if ( inode(mm) < jnode(mm) ) then
      ivertx = inode(mm)
    else
      ivertx = jnode(mm)
    end if

    jvertx = inode(mm) + jnode(mm) - ivertx
    loop = mm - 1
    mm = 0

    do i = 1, loop

      ilast = inode(i)

      if ( ilast == jvertx ) then
        ilast = ivertx
      end if

      if ( ilast == nn ) then
        ilast = jvertx
      end if

      jlast = jnode(i)
      if ( jlast == jvertx ) then
        jlast = ivertx
      end if

      if ( jlast == nn ) then
        jlast = jvertx
      end if

      if ( ilast == ivertx ) then
        if ( cpoly1(jlast) /= 0 ) then
          cycle
        end if
        cpoly1(jlast) = 1
      end if

      if ( jlast == ivertx ) then
        if ( cpoly1(ilast) /= 0 ) then
          cycle
        end if
        cpoly1(ilast) = 1
      end if

      mm = mm + 1
      inode(mm) = ilast
      jnode(mm) = jlast

    end do

    nn = nn - 1

  end do

  return
end
subroutine graph_arc_conect ( nnode, nedge, inode, jnode, iroot, ncut, nbridg, &
  cutnod, bridge, next )

!*****************************************************************************80
!
!! GRAPH_ARC_CONECT finds the bridges, blocks and cut nodes of an undirected graph.
!
!  Discussion:
!
!    The graph need not be connected.
!
!    Assume that input graph G is connected and has the set of nodes V
!    numbered from 1 to NNODE.  A tree T rooted at an arbitrary node R will be
!    grown to span G.
!
!      P(I) is the unique predecessor of node I in the tree,
!      D(I) is the distance from node I to the root of T,
!      B(I) is the label assigned to the edge (I, P(I)), and
!      H(I) is a Boolean variable for marking edge I.
!
!    STEP 1.
!
!      Set B(I) = 0, H(I) = FALSE for all I.
!      Choose an arbitrary node R as the root.
!      Initially, T consists of the single node R,
!
!        D(R) = 0, X is empty, Y = T, Z = V - {R}.
!
!    STEP 2.
!
!    If Y is nonempty then select the most recent member of Y, say U,
!    delete U from Y and continue from Step 3.
!
!    If Y is empty then count the number of blocks of G by noting that edge
!    (I, P(I)) belongs to block B(I); moreover, if B(I) is equal to -I, then
!    the edge itself is a block of G.
!
!    If G has only one block then stop; otherwise, G has at least one
!    cut node.  The cut nodes are identified by the property that
!    node I is a cut node if and only if there are two or more distinct
!    labels on edges of T through I. Stop.
!
!    STEP 3.
!
!    For each edge (U, V), where V is in Z, do the following:
!
!      add (U, V) to T and transfer V from Z to Y,
!      set P(V) = U, D(V) = D(U) + 1, b(V) = - V.
!
!    For each edge (U, V), where V is in Y, do the following:
!
!      If 0 < B(V) then set H(B(V)) = TRUE,
!      set B(V) = U and L = max (L, D(U) - D(P(V)))
!
!    STEP 4.
!
!    If L = 0 then add U to X and return to Step 2.
!
!    Otherwise, for each edge, say (I, P(I)), on the path of length
!    L from U to the root of the tree, set
!
!      H(B(I)) = TRUE and B(I) = U.
!
!    For each edge (J, P(J)) for which H(B(J)) = TRUE, set
!
!      B(J) = U.
!
!    Add U to X and return to Step 2.
!
!
!    The running time of the algorithm is bounded by O(NNODE**2).
!
!
!    If the input graph is known to be connected, then one call to CONECT
!    rooted at an arbitrary node will find all the cut nodes and bridges
!    of the graph.  If it is not known whether the input graph is connected,
!    then the calling of CONECT can be enclosed in a loop iterating from
!    IROOT = 1 to NNODE.  As usual, the array NEXT should be initialized to zero
!    at the start of the loop.  Each execution of CONECT will assign a
!    positive integer to NEXT(I), for every node I in the same component of
!    IROOT.  If there are K components in the input graph, the loop will
!    be executed K times.  This is illustrated in the test example.
!
!  Modified:
!
!    23 July 2000
!
!  Author:
!
!    Hang Tong Lau
!
!  Reference:
!
!    Hang Tong Lau,
!    Algorithms on Graphs,
!    Tab Books, 1989,
!    QA166 L38.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) NEDGE, the number of edges.
!
!    Input, integer ( kind = 4 ) INODE(NEDGE), JNODE(NEDGE); INODE(I), JNODE(I) are the
!    end nodes of the I-th edge in the input graph.
!
!    Input, integer ( kind = 4 ) IROOT, a specified node for the root of the tree.
!
!    Output, integer ( kind = 4 ) NCUT, the number of cut nodes.
!
!    Output, integer ( kind = 4 ) NBRIDG, the number of bridges.
!
!    Output, integer ( kind = 4 ) CUTNOD(NNODE); during execution of the routine, CUTNOD(I)
!    is the unique predecessor of node I in the tree.  On output, the cut
!    nodes are stored in
!
!      CUTNOD(I), I = 1, 2, ..., NCUT.
!
!    Output, integer ( kind = 4 ) BRIDGE(NEDGE); during execution of the routine, BRIDGE(I)
!    is the end node of the I-th edge in the forward star representation of
!    the graph. On output, if BRIDGE(I) = 1, then (INODE(I),JNODE(J)) is
!    a bridge, and if BRIDGE(I) = 0 then the I-th edge is not a bridge.
!
!    Input/output, integer ( kind = 4 ) NEXT(NNODE);
!    On input for the first component, NEXT(1:N) must be initialized to zero.
!    On output, NEXT(I) is the number of the node following node I on the path
!    from node I to the root of the tree.
!
!    Workspace, integer ARCFIR(NNODE); ARCFIR(I) is the number of the first
!    edge starting at node I in the forward star representation of the graph.
!
!    Workspace, integer LABEL(NNODE); labels assigned to edges in the tree.
!
!    Workspace, LENGTH(NNODE); LENGTH(I) is the distance from node I to the
!    root of the tree.
!
!    Workspace, logical NEW(NNODE); indicates whether the labels are replaced.
!
  implicit none

  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) arcfir(nnode)
  integer ( kind = 4 ) bridge(nedge)
  integer ( kind = 4 ) cutnod(nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iedge
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) index
  integer ( kind = 4 ) inode(nedge)
  integer ( kind = 4 ) iroot
  integer ( kind = 4 ) itemp
  integer ( kind = 4 ) iup
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jnode(nedge)
  logical join
  integer ( kind = 4 ) k
  integer ( kind = 4 ) label(nnode)
  integer ( kind = 4 ) len1
  integer ( kind = 4 ) len2
  integer ( kind = 4 ) length(nnode)
  integer ( kind = 4 ) low
  integer ( kind = 4 ) nblock
  integer ( kind = 4 ) nbridg
  integer ( kind = 4 ) ncut
  logical new(nnode)
  integer ( kind = 4 ) next(nnode)
  integer ( kind = 4 ) node1
  integer ( kind = 4 ) node2
  integer ( kind = 4 ) node3
  integer ( kind = 4 ) node4
  integer ( kind = 4 ) nodeu
  integer ( kind = 4 ) nodev
!
!  Initialize the data.
!
  cutnod(1:nnode) = 0
  label(1:nnode) = 0
  length(1:nnode) = 0
  new(1:nnode) = .false.
  bridge(1:nedge) = 0
!
!  Set up the forward star representation of the graph.
!
  k = 0

  do i = 1, nnode - 1

    arcfir(i) = k + 1

    do j = 1, nedge

      if ( inode(j) == i .and. inode(j) < jnode(j) ) then
        k = k + 1
        bridge(k) = jnode(j)
      else if ( jnode(j) == i .and. jnode(j) < inode(j) ) then
        k = k + 1
        bridge(k) = inode(j)
      end if

    end do

  end do

  arcfir(nnode) = nedge + 1
!
!  Store information about the root.
!
  length(iroot) = 0
  next(iroot) = - 1
  label(iroot) = - iroot
  index = 1
  cutnod(1) = iroot
  iedge = 2

40    continue

  node3 = cutnod(index)
  index = index - 1
  next(node3) = - next(node3)
  len1 = 0
!
!  Big loop.
!
  do node2 = 1, nnode

    join = .false.

    if ( node2 /= node3 ) then

      if ( node2 < node3 ) then
        nodeu = node2
        nodev = node3
      else
        nodeu = node3
        nodev = node2
      end if

      low = arcfir(nodeu)
      iup = arcfir(nodeu+1) - 1

      do k = low, iup

        if ( bridge(k) == nodev ) then
          join = .true.
          exit
        end if

      end do

    end if

    if ( join ) then

      node1 = next(node2)

      if ( node1 == 0 ) then

        next(node2) = - node3
        index = index + 1
        cutnod(index) = node2
        length(node2) = length(node3) + 1
        label(node2) = - node2

      else if ( node1 < 0 ) then
!
!  Next block.
!
        node4 = label(node2)

        if ( 0 < node4 ) then
          new(node4) = .true.
        end if

        label(node2) = node3
        len2 = length(node3) - length(-node1)

        if ( len1 < len2 ) then
          len1 = len2
        end if

      end if

    end if

  end do

  if ( 0 < len1 ) then

    j = node3

80      continue

    len1 = len1 - 1

    if ( 0 <= len1 ) then

      itemp = label(j)

      if ( 0 < itemp ) then
        new(itemp) = .true.
      end if

      label(j) = node3
      j = next(j)
      go to 80

    end if

    do i = 1, nnode

      itemp = label(i)

      if ( 0 < itemp ) then

        if ( new(itemp) ) then
          label(i) = node3
        end if

      end if

    end do

  end if
!
!  And now...
!
  iedge = iedge + 1

  if ( iedge <= nnode .and. 0 < index ) then
    go to 40
  end if

  next(iroot) = 0
  node3 = cutnod(1)
  next(node3) = abs ( next(node3) )
!
!  Count the number of bridges and blocks.
!
  nbridg = 0
  nblock = 0

  do i = 1, nnode

    if ( i /= iroot ) then

      node3 = label(i)

      if ( node3 < 0 ) then

        nblock = nblock + 1
        nbridg = nbridg + 1
        label(i) = nnode + nblock

      else

        if ( node3 <= nnode .and. 0 < node3 ) then

          nblock = nblock + 1
          node4 = nnode + nblock

          do j = i, nnode

            if ( label(j) == node3 ) then
              label(j) = node4
            end if

          end do

        end if

      end if

    end if

  end do
!
!  And now...
!
  do i = 1, nnode
    if ( 0 < label(i) ) then
      label(i) = label(i) - nnode
    end if
  end do
!
!  Find the value of I for which NEXT(I) = IROOT.
!
  i = 1

  do while ( next(i) /= iroot )

    i = i + 1

    if ( nnode < i ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'GRAPH_ARC_CONECT - Fatal error!'
      write ( *, '(a,i8)' ) '  Illegal node index, I = ', i
      stop
    end if

  end do
!
!  Reset the entries of LABEL.
!
  label(iroot) = label(i)

  do i = 1, nnode

    node1 = next(i)

    if ( 0 < node1 ) then

      itemp = abs ( label(node1) )

      if ( abs ( label(i) ) /= itemp ) then
        label(node1) = - itemp
      end if

    end if

  end do
!
!  Count the number of cuts.
!
  ncut = 0

  do i = 1, nnode

    if ( label(i) < 0 ) then
      ncut = ncut + 1
    end if

  end do
!
!  Store the cut nodes.
!
  j = 0

  do i = 1, nnode

    if ( label(i) < 0 ) then
      j = j + 1
      cutnod(j) = i
    end if

  end do
!
!  Find the end nodes.
!
  length(1:nnode) = 0

  do i = 1, nedge
    j = inode(i)
    length(j) = length(j) + 1
    j = jnode(i)
    length(j) = length(j) + 1
  end do

  do i = 1, nnode

    if ( length(i) == 1 ) then

      if ( 0 < label(i) ) then
        label(i) = - label(i)
      end if

    end if

  end do
!
!  Store the bridges.
!
  bridge(1:nedge) = 0

  do i = 1, nnode

    if ( label(i) < 0 ) then

      do j = 1, nedge

        nodev = 0

        if ( i == inode(j) ) then
          nodev = jnode(j)
        end if

        if ( i == jnode(j) ) then
          nodev = inode(j)
        end if

        if ( 0 < nodev ) then

          if ( label(nodev) < 0 ) then

            if ( label(i) /= label(nodev) ) then

              bridge(j) = 1

            else

              k = - label(i)

              do ii = 1, nnode
                if ( ii /= i .and. ii /= nodev ) then
                  if ( label(ii) == k ) then
                    go to 230
                  end if
                end if
              end do

              bridge(j) = 1

            end if

          end if

        end if

      end do

    end if

230     continue

  end do

  return
end
subroutine graph_arc_degree ( nnode, nedge, inode, jnode, degree )

!*****************************************************************************80
!
!! GRAPH_ARC_DEGREE finds the degree of the nodes of an undirected graph.
!
!  Modified:
!
!    04 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) NEDGE, the number of edges.
!
!    Input, integer ( kind = 4 ) INODE(NEDGE), JNODE(NEDGE); the I-th edge of the graph
!    extends from node INODE(I) to JNODE(I).
!
!    Output, integer ( kind = 4 ) DEGREE(NNODE), the degree of each node.
!
  implicit none

  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) degree(nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) inode(nedge)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jnode(nedge)

  do i = 1, nnode

    degree(i) = 0

    do j = 1, nedge

      if ( inode(j) == i ) then
        degree(i) = degree(i) + 1
      end if

      if ( jnode(j) == i ) then
        degree(i) = degree(i) + 1
      end if

    end do

  end do

  return
end
subroutine graph_arc_edge_con ( nnode, nedge, inode, jnode, edge_con )

!*****************************************************************************80
!
!! GRAPH_ARC_EDGE_CON finds the edge-connectivity of a connected graph.
!
!  Method:
!
!    A graph G has edge connectivity K if, given any pair of distinct nodes
!    I and J, there are K paths from I to J, no two of which use a common
!    edge.
!
!    Thus, in particular, if a graph G is Hamiltonian, it must have
!    edge connectivity at least 2.  For we can simply take the Hamiltonian
!    circuit, and use the part from I to J as the first path, and the
!    part from J to I as the second, simply reversing the direction
!    of traversal.
!
!    To determine the edge connectivity, for each J from 2 to NNODE do
!    the following:
!
!      Take node 1 as the source, node J as the sink in G, assign a unit
!      capacity to all edges in both directions, and find the value of the
!      maximum flow G(J) in the resulting network.
!
!    The edge connectivity is then equal to the minimum of G(2:NNODE).
!
!    This routine finds the edge connectivity of a given undirected graph with
!    the help of DIGRAPH_ARC_NFLOW.
!
!    The maximum network flow algorithm requires O(NNODE**3) operations.  The
!    edge connectivity of a graph will therefore be found in O(NNODE**4)
!    operations.
!
!  Modified:
!
!    23 July 2000
!
!  Author:
!
!    Hang Tong Lau
!
!  Reference:
!
!    Hang Tong Lau,
!    Algorithms on Graphs,
!    Tab Books, 1989,
!    QA166 L38.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) NEDGE, the number of edges.
!
!    Input, integer ( kind = 4 ) INODE(NEDGE), JNODE(NEDGE); INODE(I), JNODE(I) are the
!    end nodes of the I-th edge in the graph.
!
!    Output, integer ( kind = 4 ) EDGE_CON, the edge-connectivity of the graph.
!
!    Workspace, integer IEDGE(4*NEDGE), JEDGE(4*NEDGE).
!
!    Workspace, integer CAPAC(4*NEDGE).
!
!    Workspace, integer MINCUT(NNODE).
!
!    Workspace, integer FLOW(4*NEDGE).
!
!    Workspace, integer NODE_FLOW(NNODE).
!
  implicit none

  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) capac(4*nedge)
  integer ( kind = 4 ) edge_con
  integer ( kind = 4 ) flow(4*nedge)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iedge(4*nedge)
  integer ( kind = 4 ) inode(nedge)
  integer ( kind = 4 ) isink
  integer ( kind = 4 ) isorce
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jedge(4*nedge)
  integer ( kind = 4 ) jnode(nedge)
  integer ( kind = 4 ) m4
  integer ( kind = 4 ) mincut(nnode)
  integer ( kind = 4 ) node_flow(nnode)
!
!  Duplicate the edges.
!
  j = 0

  do i = 1, nedge

    j = j + 1
    iedge(j) = inode(i)
    jedge(j) = jnode(i)
    capac(j) = 1

    j = j + 1
    iedge(j) = jnode(i)
    jedge(j) = inode(i)
    capac(j) = 0

    j = j + 1
    iedge(j) = jnode(i)
    jedge(j) = inode(i)
    capac(j) = 1

    j = j + 1
    iedge(j) = inode(i)
    jedge(j) = jnode(i)
    capac(j) = 0

  end do
!
!  Call on the network flow algorithm.
!
  edge_con = nnode
  isorce = 1
  m4 = 4 * nedge

  do isink = 2, nnode

    call digraph_arc_nflow ( nnode, m4, iedge, jedge, capac, isorce, isink, &
      mincut, flow, node_flow )

    if ( node_flow(isorce) < edge_con ) then
      edge_con = node_flow(isorce)
    end if

  end do

  return
end
subroutine graph_arc_euler ( nnode, nedge, inode, jnode, success, trail )

!*****************************************************************************80
!
!! GRAPH_ARC_EULER returns an Euler circuit in an undirected graph.
!
!  Discussion:
!
!    An Euler circuit of an undirected graph is a walk which starts and ends at
!    the same node and uses each edge exactly once.  A graph is eulerian if
!    it has an Euler circuit.  The problem is to decide whether a given
!    connected graph is eulerian and to find an Euler circuit if the answer
!    is affirmative.
!
!  Method:
!
!    A connected graph is eulerian if and only if it has no node of odd degree.
!
!    This characterization gives a straightforward procedure to decide whether
!    a graph is eulerian.  Furthermore, an Euler circuit in an eulerian graph
!    G of NEDGE edges can be determined by the following method:
!
!      STEP 1: Choose any node U as the starting node, and traverse any edge
!      ( U, V ) incident to node U, and than traverse any unused edge incident
!      to node U.  Repeat this process of traversing unused edges until the
!      starting node U is reached.  Let P be the resulting walk consisting of
!      all used edges.  If all edges of G are in P, than stop.
!
!      STEP 2: Choose any unused edge ( X,  Y) in G such that X is
!      in P and Y is not in P.  Use node X as the starting node and
!      find another walk Q using all unused edges as in step 1.
!
!      STEP 3: Walk P and walk Q share a common node X, they can be merged
!      to form a walk R by starting at any node S of P and to traverse P
!      until node X is reached; than, detour and traverse all edges of Q
!      until node X is reached and continue to traverse the edges of P until
!      the starting node S is reached.  Set P = R.
!
!      STEP 4: Repeat steps 2 and 3 until all edges are used.
!
!    The running time of the algorithm is O ( NEDGE ).
!
!  Note:
!
!    The graph is assumed to be connected.
!
!  Modified:
!
!    23 July 2000
!
!  Author:
!
!    Hang Tong Lau
!
!  Reference:
!
!    Hang Tong Lau,
!    Algorithms on Graphs,
!    Tab Books, 1989,
!    QA166 L38.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) NEDGE, the number of edges.
!
!    Input, integer ( kind = 4 ) INODE(NEDGE), JNODE(NEDGE); the I-th edge starts at node
!    INODE(I) and ends at node JNODE(I).
!
!    Output, logical SUCCESS, is TRUE if an Euler circuit was found,
!    and FALSE otherwise.
!
!    Output, integer ( kind = 4 ) TRAIL(NEDGE).  TRAIL(I) is the edge number of the I-th
!    edge in the Euler circuit.
!
  implicit none

  integer ( kind = 4 ) nedge

  logical candid(nedge)
  integer ( kind = 4 ) endnod(nedge)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) inode(nedge)
  integer ( kind = 4 ) istak
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jnode(nedge)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) len
  integer ( kind = 4 ) lensol
  integer ( kind = 4 ) lenstk
  integer ( kind = 4 ) nnode
  integer ( kind = 4 ) stack(2*nedge)
  logical success
  integer ( kind = 4 ) trail(nedge)
!
!  Check if the graph is eulerian.
!
  trail(1:nedge) = 0
  endnod(1:nedge) = 0

  do i = 1, nedge
    j = inode(i)
    endnod(j) = endnod(j) + 1
    j = jnode(i)
    endnod(j) = endnod(j) + 1
  end do

  do i = 1, nnode

    if ( mod ( endnod(i), 2 ) /= 0 ) then
      success = .false.
      return
    end if

  end do
!
!  The graph is eulerian; find an Euler circuit.
!
  success = .true.
  lensol = 1
  lenstk = 0
!
!  Find the next edge.
!
  do

    if ( lensol == 1 ) then

      endnod(1) = inode(1)
      stack(1) = 1
      stack(2) = 1
      lenstk = 2

    else

      l = lensol - 1

      if ( lensol /= 2 ) then
        endnod(l) = inode(trail(l)) + jnode(trail(l)) - endnod(l-1)
      end if

      k = endnod(l)

      do i = 1, nedge
        candid(i) = ( k == inode(i) ) .or. ( k == jnode(i) )
      end do

      do i = 1, l
        candid(trail(i)) = .false.
      end do

      len = lenstk

      do i = 1, nedge

        if ( candid(i) ) then
          len = len + 1
          stack(len) = i
        end if

      end do

      stack(len+1) = len - lenstk
      lenstk = len + 1

    end if

    do

      istak = stack(lenstk)
      lenstk = lenstk - 1

      if ( istak /= 0 ) then
        exit
      end if

      lensol = lensol - 1

      if ( lensol == 0 ) then
        return
      end if

    end do

    trail(lensol) = stack(lenstk)
    stack(lenstk) = istak - 1

    if ( lensol == nedge ) then
      exit
    end if

    lensol = lensol + 1

  end do

  return
end
subroutine graph_arc_fcycle ( nnode, nedge, inode, jnode, numcyc, numcmp, ncyc )

!*****************************************************************************80
!
!! GRAPH_ARC_FCYCLE finds a fundamental set of cycles in an undirected graph.
!
!  Discussion:
!
!    Let T be a spanning tree of an undirected graph G.  The fundamental set
!    of cycles of G corresponding to T is the set of cycles of G consisting
!    of one edge (I,J) of G-T together with the unique path between node I
!    and node J in T.  The problem is to find a fundamental set of cycles
!    in a given undirected graph that is not necessarily connected.
!
!    Note that the graph may be unconnected.
!
!  Method:
!
!    Let G be the given undirected graph of NNODE nodes.  First, find all the
!    connected components of G.  Then, the fundamental set of cycles of G
!    can be found for each component H of G as follows.
!
!    STEP 1. Let E be the set of edges and V the set of nodes of H.
!    Take any node Z from V as the root of the tree consisting of
!    the single node. Set
!
!      T = {Z}, S = V.
!
!    STEP 2. Let X be any node in both T and S. If such a node does not
!    exist, then stop.
!
!    STEP 3. Consider each edge (X,Y) in E.  If Y is in T, then generate the
!    fundamental cycle consisting of edge (X,Y) together with the unique path
!    between X and Y in the tree, and delete the edge (X,Y) from E.  If Y is
!    not in T, then add the edge (X,Y) to the tree, add the node Y to T, and
!    delete the edge (X,Y) from E.
!
!   STEP 4. Remove the node X from S and return to Step 2.
!
!   The processing time of the algorithm is O(N**2).
!
!  Modified:
!
!    23 July 2000
!
!  Author:
!
!    Hang Tong Lau
!
!  Reference:
!
!    Hang Tong Lau,
!    Algorithms on Graphs,
!    Tab Books, 1989,
!    QA166 L38.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) NEDGE, the number of edges.
!
!    Input, integer ( kind = 4 ) INODE(NEDGE), JNODE(NEDGE); INODE(I) and
!    JNODE(I) are the end nodes of the I-th edge in the input graph.
!
!    Output, integer ( kind = 4 ) NUMCYC, the number of independent cycles in the graph.
!
!    Output, integer ( kind = 4 ) NUMCMP, the number of components of the graph.
!
!    Output, integer ( kind = 4 ) NCYC(NNODE); each time the output WRITE statement in FCYCLE
!    is executed, a fundamental cycle
!
!      NCYC(i), i = 1, 2,. . ., LEN
!
!    is generated, where LEN is a local variable in FCYCLE.
!
!    Workspace, integer FWDARC(NEDGE); FWDARC(i) is the ending node of the I-th
!    edge in the forward star representation of the graph.
!
!    Workspace, integer ARCFIR(NNODE); ARCFIR(I) is the number of the first
!    edge starting at node I in the forward star representation.
!
!    Workspace, integer NEXT(NNODE); NEXT(I) is the number of the node following
!    node I on the path from node I to the root of the tree.
!
!    Workspace, integer IPOINT(NNODE); pointer array.
!
  implicit none

  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) arcfir(nnode)
  integer ( kind = 4 ) fwdarc(nedge)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iedge
  integer ( kind = 4 ) index
  integer ( kind = 4 ) inode(nedge)
  integer ( kind = 4 ) ipoint(nnode)
  integer ( kind = 4 ) iroot
  integer ( kind = 4 ) iup
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jnode(nedge)
  logical join
  integer ( kind = 4 ) k
  integer ( kind = 4 ) len
  integer ( kind = 4 ) low
  integer ( kind = 4 ) ncyc(nnode)
  integer ( kind = 4 ) next(nnode)
  integer ( kind = 4 ) node1
  integer ( kind = 4 ) node2
  integer ( kind = 4 ) node3
  integer ( kind = 4 ) nodeu
  integer ( kind = 4 ) nodev
  integer ( kind = 4 ) numcmp
  integer ( kind = 4 ) numcyc
!
!  Set up the forward star representation of the graph.
!
  k = 0

  do i = 1, nnode - 1

    arcfir(i) = k + 1

    do j = 1, nedge

      if ( inode(j) == i .and. inode(j) < jnode(j) ) then
        k = k + 1
        fwdarc(k) = jnode(j)
      else if ( nodev == i .and. jnode(j) < inode(j) ) then
        k = k + 1
        fwdarc(k) = inode(j)
      end if

    end do

  end do

  arcfir(nnode) = nedge + 1
!
!  Initialize the data.
!
  next(1:nnode) = 0
  numcmp = 0
  numcyc = 0
!
!  Loop, using each node successively as the root.
!
  do iroot = 1, nnode

    if ( next(iroot) == 0 ) then

      numcmp = numcmp + 1
      next(iroot) = - 1
      index = 1
      ipoint(1) = iroot
      iedge = 2
!
!  Pick NODE3.
!
40    continue

      node3 = ipoint(index)
      index = index - 1
      next(node3) = - next(node3)
!
!  Decide whether to include an arc from NODE3 to NODE2.
!
      do node2 = 1, nnode

        join = .false.

        if ( node2 /= node3 ) then

          if ( node2 < node3 ) then
            nodeu = node2
            nodev = node3
          else
            nodeu = node3
            nodev = node2
          end if

          low = arcfir(nodeu)
          iup = arcfir(nodeu+1) - 1

          do k = low, iup

            if ( fwdarc(k) == nodev ) then
              join = .true.
              exit
            end if

          end do

        end if

        if ( join ) then

          node1 = next(node2)

          if ( node1 == 0 ) then

            next(node2) = - node3
            index = index + 1
            ipoint(index) = node2

          else

            if ( node1 < 0 ) then
!
!  Generate the next cycle.
!
              numcyc = numcyc + 1
              len = 3
              node1 = - node1
              ncyc(1) = node1
              ncyc(2) = node2
              ncyc(3) = node3
              i = node3

              do

                j = next(i)

                if ( j == node1 ) then
                  exit
                end if

                len = len + 1
                ncyc(len) = j
                i = j

              end do
!
!  Output a cycle.
!
              write ( *, 80 ) numcyc, ncyc(1:len)

80    format ( ' Nodes in cycle number', i3, ':', 1x, 20i4 )

            end if

          end if

        end if

      end do

      iedge = iedge + 1

      if ( iedge <= nnode .and. 0 < index ) then
        go to 40
      end if

      next(iroot) = 0
      node3 = ipoint(1)
      next(node3) = abs ( next(node3) )

    end if

  end do

  return
end
subroutine graph_arc_hamcyc ( nnode, nedge, inode, jnode, success, hcycle )

!*****************************************************************************80
!
!! GRAPH_ARC_HAMCYC finds a Hamiltonian circuit in a graph.
!
!  Discussion:
!
!    A hamiltonian cycle in a graph G is a cycle containing every
!    node of G, and a graph is hamiltonian if it has a hamiltonian cycle.
!    The problem is to decide whether a given graph is hamiltonian and
!    to find a hamiltonian cycle if the answer is affirmative.
!
!  Method:
!
!    A hamiltonian cycle in a given graph G of NNODE nodes will be found by an
!    exhaustive search.  Start with a single node, say node 1, as the
!    partially constructed cycle.  The cycle is grown by a backtracking
!    procedure until a hamiltonian cycle is formed.  More precisely, let
!
!      H(1), H(2), ..., H(K-1)
!
!    be the partially constructed cycle at the K-th stage.  Then, the
!    set of candidates, for the next element H(K) is the set of all
!    nodes U in G such that
!
!      if K = 1, than U = 1;
!      if 1 < K, then ( H(K-1), U ) is an edge in G, and U is
!         distinct from H(1), H(2), ..., H(K-1); furthermore,
!      Furthermore, if K = NNODE, then ( U, H(1) ) is an edge in G and
!      U < H(2).
!
!    The search is complete when K = NNODE.  If the set of candidates
!    is empty at any stage, then the graph is not hamiltonian.
!
!  Modified:
!
!    23 July 2000
!
!  Author:
!
!    Hang Tong Lau
!
!  Reference:
!
!    Hang Tong Lau,
!    Algorithms on Graphs,
!    Tab Books, 1989,
!    QA166 L38.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) NEDGE, the number of edges.
!
!    Input, logical DIRECT, is TRUE if the graph is directed, and FALSE
!    otherwise.
!
!    Input, integer ( kind = 4 ) INODE(NEDGE), JNODE(NEDGE); the I-th edge of
!    the graph extends from node INODE(I) to JNODE(I).
!
!    Output, logical SUCCESS, is TRUE if a hamiltonian cycle was found,
!    and FALSE otherwise.
!
!    Output, integer ( kind = 4 ) HCYCLE(NNODE); the nodes of the hamiltonian cycle are
!    stored in order in HCYCLE.
!
!    Workspace, integer FWDARC(M2); FWDARC(I) is the ending node of the
!    I-th edge in the forward star representation of the graph.
!
!    Workspace, integer ARCFIR(NNODE+1); ARCFIR(I) is the number of the first
!    edge starting at node I in the forward star representation of the graph.
!
!    Workspace, integer STACK(2*NEDGE).
!
!    Workspace, logical CONNECT(NNODE).
!
  implicit none

  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) arcfir(nnode+1)
  logical connect(nnode)
  integer ( kind = 4 ) fwdarc(2*nedge)
  integer ( kind = 4 ) hcycle(nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) inode(nedge)
  integer ( kind = 4 ) istak
  integer ( kind = 4 ) iup
  integer ( kind = 4 ) jnode(nedge)
  logical join
  integer ( kind = 4 ) k
  integer ( kind = 4 ) len
  integer ( kind = 4 ) len1
  integer ( kind = 4 ) len2
  integer ( kind = 4 ) lensol
  integer ( kind = 4 ) lenstk
  integer ( kind = 4 ) low
  integer ( kind = 4 ) stack(2*nedge)
  logical success
!
!  Set up the forward star representation of the graph.
!
  call graph_arc_to_star ( nnode, nedge, inode, jnode, arcfir, fwdarc )
!
!  Initialize.
!
  lensol = 1
  lenstk = 0
  hcycle(1:nnode) = 0
!
!  Find the next node.
!
  do
!
!  Start with
!    STACK(1) = node 1,
!    STACK(2) = level 1.
!
    if ( lensol == 1 ) then

      stack(1) = 1
      stack(2) = 1
      lenstk = 2

    else

      len1 = lensol - 1
      len2 = hcycle(len1)
!
!  Seek a connection from the last node on the tentative list to any node I.
!
      do i = 1, nnode

        connect(i) = .false.
        low = arcfir(len2)
        iup = arcfir(len2+1) - 1

        do k = low, iup
          if ( fwdarc(k) == i ) then
            connect(i) = .true.
          end if
        end do

      end do
!
!  Disregard connections that take you to nodes already on the tentative list.
!
      do i = 1, len1
        connect ( hcycle(i) ) = .false.
      end do

      len = lenstk

      if ( lensol /= nnode ) then

        do i = 1, nnode
          if ( connect(i) ) then
            len = len + 1
            stack(len) = i
          end if
        end do

        stack(len+1) = len - lenstk
        lenstk = len + 1

      else

        do i = 1, nnode

          if ( connect(i) ) then

            if ( hcycle(2) < i ) then
              stack(len+1) = len - lenstk
              lenstk = len + 1
              go to 110
            end if

            join = .false.
            low = arcfir(i)
            iup = arcfir(i+1) - 1

            do k = low, iup

              if ( fwdarc(k) == 1 ) then
                join = .true.
                exit
              end if

            end do

            if ( join ) then
              lenstk = lenstk + 2
              stack(lenstk-1) = i
              stack(lenstk) = 1
            else
              stack(len+1) = len - lenstk
              lenstk = len + 1
            end if

            go to 110

          end if

        end do

        stack(len+1) = len - lenstk
        lenstk = len + 1

      end if

    end if
!
!  Search further.
!
  110   continue

    istak = stack(lenstk)
    lenstk = lenstk - 1

    if ( istak == 0 ) then

      lensol = lensol - 1

      if ( lensol == 0 ) then
        success = .false.
        return
      end if

      go to 110

    end if

    hcycle(lensol) = stack(lenstk)
    stack(lenstk) = istak - 1

    if ( lensol == nnode ) then
      success = .true.
      exit
    end if

    lensol = lensol + 1

  end do

  return
end
subroutine graph_arc_makeg ( nnode, k, inode, jnode )

!*****************************************************************************80
!
!! GRAPH_ARC_MAKEG constructs a graph with NNODE nodes and K-connectivity.
!
!  Discussion:
!
!    Let NNODE and K be given positive integers.  The problem is to construct
!    a K-connected graph G(K,NNODE) on NNODE nodes with as few edges as
!    possible.  Observe that for K = 1, the graph G(1,NNODE) is a spanning
!    tree.  Consequently, it is assumed that K is at least 2.  Moreover,
!    it is known that G(K,NNODE) has exactly [(NNODE*K)/2] edges, where
!    [X] is the smallest integer greater than or equal to X.
!
!  Method:
!
!    Label the nodes of the graph by the integers 0, 1, 2, .. ., NNODE - 1
!
!    CASE 1. K is even.
!
!      Let K = 2*T.  The graph G(2*T,NNODE) is constructed as follows.
!      First, draw an NNODE-gon, that is, add the edges
!
!      (0,1), (1,2), (2,3), ..., (NNODE-2, NNODE-1), (NNODE-1, 0).
!
!      Then join nodes I and J if and only if
!
!        abs ( I - J ) = P (mod NNODE), where 2 <= P <= T.
!
!    CASE 2. K is odd, NNODE is even.
!
!      Let K = 2*T+1.  The graph G(2*T+1,NNODE) is constructed by first drawing
!      G(2*T,NNODE) as above, and then joining node I to node
!
!        I + (NNODE/2), for 0 <= I <= NNODE/2.
!
!    CASE 3. K is odd, NNODE is odd.
!
!      Let K = 2*T + 1.  The graph G(2*T+1,NNODE) is constructed by first
!      drawing  G(2*T,NNODE) as above, and then joining
!
!        node 0 to node (NNODE-1)/2,
!        node 0 to node (NNODE+1)/2,
!        node I to node I + (NNODE+1)/2, for 1 <= I < (NNODE - 1)/2.
!
!  Modified:
!
!    23 July 2000
!
!  Author:
!
!    Hang Tong Lau
!
!  Reference:
!
!    Hang Tong Lau,
!    Algorithms on Graphs,
!    Tab Books, 1989,
!    QA166 L38.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) K, the desired connectivity of the graph.  K
!    must be at least 2.
!
!    Output, integer ( kind = 4 ) INODE((NNODE*K)/2), JNODE((NNODE*K)/2); INODE(I),
!    JNODE(I) are the end nodes of the I-th edge in the graph.
!
  implicit none

  integer ( kind = 4 ) nnode
  integer ( kind = 4 ) k

  logical evenn
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iedge
  integer ( kind = 4 ) inode((nnode*k)/2)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jnode((nnode*k)/2)
  logical join
  integer ( kind = 4 ) l
  integer ( kind = 4 ) nmm
  integer ( kind = 4 ) npp

  iedge = 0
!
!  Make an NNODE-gon.
!
  do i = 1, nnode-1
    iedge = iedge + 1
    inode(iedge) = i
    jnode(iedge) = i + 1
  end do

  iedge = iedge + 1
  inode(iedge) = nnode
  jnode(iedge) = 1

  if ( k == 2 ) then
    return
  end if

  do i = 1, nnode - 1

    do j = i+1, nnode

      join = .false.

      do l = 2, k/2

        if ( mod ( l, nnode ) == j - i .or. j - i + l == nnode ) then

          join = .true.

        end if

      end do

      if ( join ) then
        iedge = iedge + 1
        inode(iedge) = i
        jnode(iedge) = j
      end if

    end do

  end do
!
!  Return if K is even.
!
  if ( mod ( k, 2 ) == 0 ) then
    return
  end if

  evenn = mod ( nnode, 2 ) == 0

  if ( evenn ) then
!
!  K is odd, NNODE is even.
!
    do i = 1, nnode / 2
      iedge = iedge + 1
      inode(iedge) = i
      jnode(iedge) = i + ( nnode / 2 )
    end do

  else
!
!  K is odd, NNODE is odd.
!
    npp = ( nnode + 1 ) / 2
    nmm = ( nnode - 1 ) / 2

    do i = 2, nmm
      iedge = iedge + 1
      inode(iedge) = i
      jnode(iedge) = i + npp
    end do

    iedge = iedge + 1
    inode(iedge) = 1
    jnode(iedge) = nmm + 1

    iedge = iedge + 1
    inode(iedge) = 1
    jnode(iedge) = npp + 1

  end if

  return
end
subroutine graph_arc_match ( nnode, nedge, inode, jnode, pair, notmat )

!*****************************************************************************80
!
!! GRAPH_ARC_MATCH finds a maximum cardinality matching in a graph.
!
!  Discussion:
!
!    A matching in an undirected graph G is a subset S of edges
!    of G in which no two edges in S are adjacent in G.  The
!    problem is to find a matching of maximum cardinality.
!
!    A matching can be used to cover the black and white squares
!    of a checkerboard with dominoes.
!
!    Similarly, if some nodes represent females, and some males,
!    with an arc between two nodes representing a possible marriage,
!    then a maximum matching is an arrangement of marriages that
!    marries off the most number of people.  (A full marriage
!    matching is possible if and only if every set of K women
!    has at least K husbands in mind, for all K).
!
!  Method:
!
!    Let S be a matching in an undirected graph G.  A node that
!    is not matched is called an exposed node.  An alternative path
!    with respect to a given matching S is a path in which the
!    edges are alternately in S and not in S.  An augmenting path
!    is an alternating path which begins with an exposed node
!    and ends with another exposed node.  A fundamental theorem
!    states that a matching S in G is maximum if and only if G
!    has no augmenting path with respect to S.
!
!    The basic method starts with an arbitrary matching Q.  An augmenting
!    path P with respect to Q is found.  Then a new matching is constructed
!    by taking those edges of Q or P that are not in both Q and P.  The
!    process is repeated and the matching is maximal when no augmenting path
!    is found.
!
!    The whole algorithm runs in O ( NNODE**3 ) time, where NNODE is the
!    number of nodes in the graph.
!
!  Modified:
!
!    23 July 2000
!
!  Author:
!
!    Hang Tong Lau
!
!  Reference:
!
!    Hang Tong Lau,
!    Algorithms on Graphs,
!    Tab Books, 1989,
!    QA166 L38.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) NEDGE, the number of edges.
!
!    Input, integer ( kind = 4 ) INODE(NEDGE), JNODE(NEDGE); the I-th edge of
!    the graph extends from node INODE(I) to JNODE(I).
!
!    Output, integer ( kind = 4 ) PAIR(NNODE), describes the matching.  Node I is
!    matched with node PAIR(I).  If PAIR(I) is 0, then node I is
!    unmatched.
!
!    Output, integer ( kind = 4 ) NOTMAT, the number of unmatched nodes.
!
!    Workspace, integer FWDARC(2*NEDGE); FWDARC(I) is the ending node of the
!    I-th edge in the forward star representation of the graph.
!
!    Workspace, integer ARCFIR(NNODE+1); ARCFIR(I) is the number of the first
!    edge starting at node I in the forward star representation of the graph.
!
!    Workspace, integer ANCES(NNODE); ANCES(I) is the grandfather of node I.
!
!    Workspace, integer QUEUE(NNODE).
!
!    Workspace, integer OUTREE(NNODE).
!
  implicit none

  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) ances(nnode)
  integer ( kind = 4 ) arcfir(nnode+1)
  integer ( kind = 4 ) fwdarc(2*nedge)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ifirst
  integer ( kind = 4 ) ilast
  integer ( kind = 4 ) inode(nedge)
  integer ( kind = 4 ) istart
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jnode(nedge)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) neigh1
  integer ( kind = 4 ) neigh2
  logical newnod
  integer ( kind = 4 ) nodei
  integer ( kind = 4 ) nodej
  integer ( kind = 4 ) nodeu
  integer ( kind = 4 ) nodev
  integer ( kind = 4 ) nodew
  logical nopath
  integer ( kind = 4 ) notmat
  logical outree(nnode)
  integer ( kind = 4 ) pair(nnode)
  integer ( kind = 4 ) queue(nnode)
!
!  Set up the forward star representation of the graph.
!
  call graph_arc_to_star ( nnode, nedge, inode, jnode, arcfir, fwdarc )
!
!   All nodes are initially unmatched.
!
  notmat = nnode
  ances(1:nnode) = 0
  pair(1:nnode) = 0

  do i = 1, nnode

    if ( pair(i) == 0 ) then

      j = arcfir(i)
      k = arcfir(i+1) - 1

      do while ( pair(fwdarc(j)) /= 0 .and. j < k )
        j = j + 1
      end do
!
!  Match a pair of nodes.
!
      if ( pair(fwdarc(j)) == 0 ) then
        pair(fwdarc(j)) = i
        pair(i) = fwdarc(j)
        notmat = notmat - 2
      end if

    end if

  end do

  do istart = 1, nnode
!
!  ISTART is not yet matched.
!
    if ( 2 <= notmat .and. pair(istart) == 0 ) then
      outree(1:nnode) = .true.
      outree(istart) = .false.
!
!   Put the root in the queue.
!
      queue(1) = istart
      ifirst = 1
      ilast = 1
      nopath = .true.

70    continue

      nodei = queue(ifirst)
      ifirst = ifirst + 1
      nodeu = arcfir(nodei)
      nodew = arcfir(nodei+1) - 1

80    continue
!
!  Examine the neighbor of NODEI.
!
      if ( nopath .and. nodeu <= nodew ) then

        if ( outree(fwdarc(nodeu)) ) then

          neigh2 = fwdarc(nodeu)
          nodej = pair(neigh2)
!
!  An augmentation path is found.
!
          if ( nodej == 0 ) then

            pair(neigh2) = nodei

            do

              neigh1 = pair(nodei)
              pair(nodei) = neigh2

              if ( neigh1 == 0 ) then
                exit
              end if

              nodei = ances(nodei)
              pair(neigh1) = nodei
              neigh2 = neigh1

            end do

            notmat = notmat - 2
            nopath = .false.

          else

            if ( nodej /= nodei ) then

              if ( nodei == istart ) then
                newnod = .true.
              else

                nodev = ances(nodei)

                do while ( nodev /= istart .and. nodev /= neigh2 )
                  nodev = ances(nodev)
                end do

                if ( nodev == istart ) then
                  newnod = .true.
                else
                  newnod = .false.
                end if

              end if
!
!  Add a tree edge.
!
              if ( newnod ) then
                outree(neigh2)  = .false.
                ances(nodej) = nodei
                ilast = ilast + 1
                queue(ilast) = nodej
              end if

            end if

          end if

        end if

        nodeu = nodeu + 1
        go to 80

      end if

      if ( nopath .and. ifirst <= ilast ) then
        go to 70
      end if

    end if

  end do

  return
end
subroutine graph_arc_mintr2 ( nnode, nedge, inode, jnode, arclen, ntree, &
  itree, jtree )

!*****************************************************************************80
!
!! GRAPH_ARC_MINTR2 finds the minimum spanning tree of a graph, using an edge array.
!
!  Discussion:
!
!    Consider an undirected graph G with given edge lengths.  The minimum
!    spanning tree problem is to find a spanning tree in G such that the
!    sum of the edge lengths in the tree is minimum.
!
!  Method:
!
!    Initially, the set T is empty.  Edges are considered for inclusion in
!    T in the nondecreasing order of their lengths.  An edge is included in
!    T if it does not form a cycle with the edges already in T.  A minimum
!    spanning tree is formed when NNODE - 1 edges are included in T.
!
!    In the implementation, the edges are partially sorted, with the smallest
!    edge at the root of a heap structure (a binary tree in which the weight
!    of every node is not greater than the weight of its sons).
!
!    The running time of the algorithm is O ( NEDGE * log NEDGE ), where
!    NEDGE is the number of edges in the graph.
!
!  Modified:
!
!    24 July 2000
!
!  Author:
!
!    Hang Tong Lau
!
!  Reference:
!
!    Hang Tong Lau,
!    Algorithms on Graphs,
!    Tab Books, 1989,
!    QA166 L38.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) NEDGE, the number of edges.
!
!    Input, integer ( kind = 4 ) INODE(NEDGE), JNODE(NEDGE); the I-th edge of
!    the graph starts at node INODE(I) and goes to node JNODE(I).
!
!    Input, real ( kind = 8 ) ARCLEN(NEDGE), contains the length of each edge.
!
!    Output, integer ( kind = 4 ) NTREE.  If the input graph is connected, then
!    NTREE will have the value NNODE-1, the number of edges in the spanning tree.
!    Otherwise, NTREE will be the number of edges in the minimum spanning
!    forest.
!
!    Output, integer ( kind = 4 ) ITREE(NNODE), JTREE(NNODE).  The edges in the
!    minimum spanning tree or forest are ( ITREE(I), JTREE(I) ).
!
  implicit none

  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) nnode

  real ( kind = 8 ) arclen(nedge)
  real ( kind = 8 ) arclen2(nedge)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) index
  integer ( kind = 4 ) index1
  integer ( kind = 4 ) index2
  integer ( kind = 4 ) inode(nedge)
  integer ( kind = 4 ) inode2(nedge)
  integer ( kind = 4 ) ipred
  integer ( kind = 4 ) ipt
  integer ( kind = 4 ) itree(nnode)
  integer ( kind = 4 ) jnode(nedge)
  integer ( kind = 4 ) jnode2(nedge)
  integer ( kind = 4 ) jtree(nnode)
  integer ( kind = 4 ) md2
  integer ( kind = 4 ) medge
  integer ( kind = 4 ) medge2
  integer ( kind = 4 ) nodeu
  integer ( kind = 4 ) nodev
  integer ( kind = 4 ) ntree
  integer ( kind = 4 ) numarc
  integer ( kind = 4 ) pred(nnode)
!
!  Copy the input data.
!
  inode2(1:nedge) = inode(1:nedge)
  jnode2(1:nedge) = jnode(1:nedge)
  arclen2(1:nedge) = arclen(1:nedge)

  pred(1:nnode) = - 1
!
!  Initialize the heap structure.
!
  i = nedge / 2

  do while ( 0 < i )

    index1 = i
    md2 = nedge / 2

    do while ( index1 <= md2 )

      index = 2 * index1

      if ( index < nedge .and. arclen2(index+1) < arclen2(index) ) then
        index2 = index + 1
      else
        index2 = index
      end if

      if ( arclen2(index2) < arclen2(index1)  ) then

        call i4_swap ( inode2(index1), inode2(index2) )
        call i4_swap ( jnode2(index1), jnode2(index2) )
        call r8_swap ( arclen2(index1), arclen2(index2) )

        index1 = index2

      else

        index1 = nedge

      end if

    end do

    i = i - 1

  end do

  medge = nedge
  ntree = 0
  numarc = 0

  do while ( ntree < nnode-1 .and. numarc < nedge )
!
!  Examine the next edge.
!
    numarc = numarc + 1
    nodeu = inode2(1)
    nodev = jnode2(1)
    ipt = nodeu
!
!  Check if NODEU and NODEV are in the same component.
!
    do while ( 0 < pred(ipt) )
      ipt = pred(ipt)
    end do

    index1 = ipt
    ipt = nodev

    do while ( 0 < pred(ipt) )
      ipt = pred(ipt)
    end do

    index2 = ipt

    if ( index1 /= index2 ) then
!
!  Include NODEU and NODEV in the minimum spanning tree.
!
      ipred = pred(index1) + pred(index2)

      if ( pred(index2) < pred(index1) ) then
        pred(index1) = index2
        pred(index2) = ipred
      else
        pred(index2) = index1
        pred(index1) = ipred
      end if

      ntree = ntree + 1
      itree(ntree) = nodeu
      jtree(ntree) = nodev

    end if
!
!  Restore the heap structure.
!
    inode2(1) = inode2(medge)
    jnode2(1) = jnode2(medge)
    arclen2(1) = arclen2(medge)
    medge = medge - 1
    index1 = 1
    medge2 = medge / 2

    do while ( index1 <= medge2 )

      index = 2 * index1

      if ( index < medge .and. arclen2(index+1) < arclen2(index) ) then
        index2 = index + 1
      else
        index2 = index
      end if

      if ( arclen2(index2) < arclen2(index1) ) then

        call i4_swap ( inode2(index1), inode2(index2) )
        call i4_swap ( jnode2(index1), jnode2(index2) )
        call r8_swap ( arclen2(index1), arclen2(index2) )

        index1 = index2

      else

        index1 = medge

      end if

    end do

  end do

  return
end
subroutine graph_arc_ncolor_print ( nedge, inode, jnode, nnode, color, title )

!*****************************************************************************80
!
!! GRAPH_ARC_NCOLOR_PRINT prints out a node-colored graph from an edge list.
!
!  Discussion:
!
!    The printout is arranged to emphasize the colors of the neighboring nodes.
!
!  Modified:
!
!    23 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEDGE, the number of edges.
!
!    Input, integer ( kind = 4 ) INODE(NEDGE), JNODE(NEDGE), the beginning and end
!    nodes of the edges.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) COLOR(NNODE), the color of each node.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) color(nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) in
  integer ( kind = 4 ) inode(nedge)
  integer ( kind = 4 ) jn
  integer ( kind = 4 ) jnode(nedge)
  character ( len = * ) title

  if ( len_trim ( title ) /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Edge  Node 1  Node 2     Color 1 Color 2'
  write ( *, '(a)' ) ' '

  do i = 1, nedge
    in = inode(i)
    jn = jnode(i)
    write ( *, '(i8,2x,i8,2x,i8,2x,i8,2x,i8)' ) i, in, jn, color(in), color(jn)
  end do

  return
end
subroutine graph_arc_planar ( n, m, inode, jnode, planar )

!*****************************************************************************80
!
!! GRAPH_ARC_PLANAR determines if a graph is planar.
!
!  Algorithm:
!
!    STEP 1. If 3*N - 6 < M, then return the message "the input
!    graph is nonplanar" and stop.
!
!    STEP 2. Perform a depth-first search on graph G to obtain a
!    digraph D so that the edges of G are divided into tree edges
!    and backward edges.  During the search, compute the low
!    point functions SMALL1 and SMALL2, where:
!
!      SMALL1(V) is the lowest node reachable from node V by a
!      sequence of tree edges followed by at most one backward edge,
!
!      SMALL2(V) is the second lowest node reachable from node V
!      in this manner.
!
!    STEP 3. Reorder the adjacency lists of D using a radix sort.
!
!    STEP 4. Use the low point functions computed from Step 2
!    and the ordering of edges from Step 3 to choose a particular
!    adjacency structure so that the nodes of D are numbered in
!    the order they are reached during any depth-first search for
!    D without changing the adjacency structure.
!
!    STEP 5. Perform a second depth-first search to select edges
!    in the reverse order to that given by the adjacency structure.
!    The purpose of this search is to prepare for the path addition
!    process in the next step.
!
!    STEP 6. Perform a third depth-first search to generate one
!    cycle and a number of edge-disjoint paths.  Each generated
!    path is added to a planar embedding which contains the
!    cycle and all the previously generated paths.
!
!    Note that any two paths may not constrain each other, or they may constrain
!    each other to have either the same embedding or the opposite embedding.
!    These dependency relations among paths can be viewed as a dependency
!    graph H.
!
!    STEP 7. This last step makes use of the fact that the dependency
!    graph H is two-colorable if and only if the original graph G
!    is planar.  Use a depth-first search to color the graph H.
!
!    If H is two-colorable, then return the message "the input
!    graph is planar"; otherwise, return the message "the input
!    graph is nonplanar" and stop.
!
!    The computation of the method is bounded by O(N), where N is the number
!    of nodes.
!
!  Modified:
!
!    03 August 2000
!
!  Author:
!
!    Hang Tong Lau
!
!  Reference:
!
!    Hang Tong Lau,
!    Algorithms on Graphs,
!    Tab Books, 1989,
!    QA166 L38.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of nodes.
!
!    Input, integer ( kind = 4 ) M, the number of edges.
!
!    Input, integer ( kind = 4 ) INODE(M), JNODE(M), the pairs of nodes that form
!    the edges of the graph.
!
!    Output, logical PLANAR, is TRUE if the graph is planar.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) arcsec(7*m+2-5*n)
  integer ( kind = 4 ) arctop(7*m+2-5*n)
  logical arctyp(7*m+2-5*n)
  integer ( kind = 4 ) auxpf1(m+m+2)
  integer ( kind = 4 ) auxpf2(m+m+2)
  integer ( kind = 4 ) auxpf3(m+m+m+3)
  integer ( kind = 4 ) auxpf4(m+m+m+3)
  integer ( kind = 4 ) descp(n)
  logical examin(m+2-n)
  logical fail
  integer ( kind = 4 ) finish(m+2-n)
  integer ( kind = 4 ) first(n+m+m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) initp
  integer ( kind = 4 ) inode(m)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jnode(m)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) level
  integer ( kind = 4 ) mark(n)
  logical middle
  integer ( kind = 4 ) mtotal
  integer ( kind = 4 ) nexte
  integer ( kind = 4 ) node1
  integer ( kind = 4 ) node2
  integer ( kind = 4 ) paint(m+2-n)
  integer ( kind = 4 ) parm1(m)
  integer ( kind = 4 ) parm2(m)
  integer ( kind = 4 ) pindex(n+n)
  logical planar
  integer ( kind = 4 ) pnum
  integer ( kind = 4 ) pw12
  integer ( kind = 4 ) pw18
  integer ( kind = 4 ) pw19
  integer ( kind = 4 ) pw20
  integer ( kind = 4 ) pw21
  integer ( kind = 4 ) pw22
  integer ( kind = 4 ) pw23
  integer ( kind = 4 ) pw24
  integer ( kind = 4 ) pw25
  integer ( kind = 4 ) qnode
  integer ( kind = 4 ) second(n+m+m)
  integer ( kind = 4 ) small1(n)
  integer ( kind = 4 ) small2(n)
  integer ( kind = 4 ) snode
  integer ( kind = 4 ) snum
  integer ( kind = 4 ) sortn(n+n+m)
  integer ( kind = 4 ) sortp1(n+n+m)
  integer ( kind = 4 ) sortp2(n+n+m)
  integer ( kind = 4 ) stakc1(m+m+2)
  integer ( kind = 4 ) stakc2(m+m+2)
  integer ( kind = 4 ) stakc3(m+m+2)
  integer ( kind = 4 ) stakc4(m+m+2)
  integer ( kind = 4 ) stacke(m+m)
  integer ( kind = 4 ) start(m+2-n)
  integer ( kind = 4 ) tnode
  integer ( kind = 4 ) tnum
  integer ( kind = 4 ) trail(n)
!
!  If there are too many edges, then we know the graph can't be planar.
!
  if ( 3 * n - 6 < m ) then
    planar = .false.
    return
  end if
!
!  Set up the graph representation.
!
  second(1:n) = 0

  mtotal = n

  do i = 1, m
    node1 = inode(i)
    node2 = jnode(i)
    mtotal = mtotal + 1
    second(mtotal) = second(node1)
    second(node1) = mtotal
    first(mtotal) = node2
    mtotal = mtotal + 1
    second(mtotal) = second(node2)
    second(node2) = mtotal
    first(mtotal) = node1
  end do
!
!  Initial depth-first search; compute the low point functions.
!
  mark(1:n) = 0
  small1(1:n) = n + 1
  small2(1:n) = n + 1

  snum = 1
  pw12 = 0
  mark(1) = 1
  parm1(1) = 1
  parm2(1) = 0
  level = 1
  middle = .false.

  do

    call dfs1 ( n, m, level, middle, snum, pw12, mark, small1, &
      small2, parm1, parm2, stacke, first, second )

    if ( level <= 1 ) then
      exit
    end if

  end do

  do i = 1, n
    if ( mark(i) <= small2(i) ) then
      small2(i) = small1(i)
    end if
  end do
!
!  Radix sort.
!
  mtotal = 2 * n
  k = 2 * n
  sortn(1:2*n) = 0

  do i = 2, 2*m, 2

    k = k + 1
    sortp1(k) = stacke(i-1)
    tnode = stacke(i)
    sortp2(k) = tnode

    if ( mark(tnode) < mark(sortp1(k)) ) then

      j = 2 * mark(tnode) - 1
      sortn(k) = sortn(j)
      sortn(j) = k

    else

      if ( mark(sortp1(k)) <= small2(tnode) ) then
        j = 2 * small1(tnode) - 1
        sortn(k) = sortn(j)
        sortn(j) = k
      else
        j = 2 * small1(tnode)
        sortn(k) = sortn(j)
        sortn(j) = k
      end if

    end if

  end do

  do i = 1, 2*n

    j = sortn(i)

    do while ( j /= 0 )
      node1 = sortp1(j)
      node2 = sortp2(j)
      mtotal = mtotal + 1
      second(mtotal) = second(node1)
      second(node1) = mtotal
      first(mtotal) = node2
      j = sortn(j)
    end do

  end do
!
!  Second depth-first search.
!
  mark(2:n) = 0
  pw12 = 0
  snum = 1
  trail(1) = 1
  parm1(1) = 1
  start(1) = 0
  finish(1) = 0
  level = 1
  middle = .false.

  do

    call dfs2 ( n, m, level, middle, snum, pw12, mark, parm1, &
      stacke, first, second )

    if ( level <= 1 ) then
      exit
    end if

  end do

  mtotal = n

  do i = 1, m
    j = i + i
    node1 = stacke(j-1)
    node2 = stacke(j)
    mtotal = mtotal + 1
    second(mtotal) = second(node1)
    second(node1) = mtotal
    first(mtotal) = node2
  end do
!
!  Path decomposition; construction of the dependency graph.
!
  pw18 = 0
  pw19 = 0
  pw24 = 0
  pw25 = 0
  initp = 0
  pnum = 1
  auxpf1(1) = 0
  auxpf1(2) = 0
  auxpf2(1) = 0
  auxpf2(2) = 0
  auxpf3(1) = 0
  auxpf3(2) = n + 1
  auxpf3(3) = 0
  auxpf4(1) = 0
  auxpf4(2) = n + 1
  auxpf4(3) = 0
  pindex(1:2*n) = 0
  nexte = m - n + 1
  arcsec(1:7*m+2-5*n) = 0
  snode = n
  descp(1) = n
  parm1(1) = 1
  level = 1
  middle = .false.

  do

    call decomp ( n, m, level, middle, initp, snode, pnum, nexte, pw18, &
      pw19, pw24, pw25, trail, descp, pindex, parm1, start, finish, first, &
      second, auxpf1, auxpf2, auxpf3, auxpf4, arcsec, arctop, arctyp )

    if ( level <= 1 ) then
      exit
    end if

  end do
!
!  Perform the two-coloring.
!
  pnum = pnum - 1
  paint(1:m+2-n) = 0
  j = pnum + 1
  examin(2:pnum+1) = .true.
  tnum = 1

  do while ( tnum <= pnum )

    parm1(1) = tnum
    paint(tnum) = 1
    examin(tnum) = .false.
    level = 1
    middle = .false.

    do

      call color2 ( n, m, level, middle, fail, parm1, paint, arcsec, &
        arctop, examin, arctyp )

      if ( fail ) then
        planar = .false.
        return
      end if

      if ( level <= 1 ) then
        exit
      end if

    end do

    do while ( .not. examin(tnum) )
      tnum = tnum + 1
    end do

  end do

  pw20 = 0
  pw21 = 0
  pw22 = 0
  pw23 = 0
  stakc1(1) = 0
  stakc1(2) = 0
  stakc2(1) = 0
  stakc2(2) = 0
  stakc3(1) = 0
  stakc3(2) = 0
  stakc4(1) = 0
  stakc4(2) = 0

  do i = 1, pnum

    qnode = start(i+1)
    tnode = finish(i+1)

    do while ( qnode <= stakc1(pw20+2) )
      pw20 = pw20 - 2
    end do

    do while ( qnode <= stakc2(pw21+2) )
      pw21 = pw21 - 2
    end do

    do while ( qnode <= stakc3(pw22+2) )
      pw22 = pw22 - 2
    end do

    do while ( qnode <= stakc4(pw23+2) )
      pw23 = pw23 - 2
    end do

    if ( paint(i) == 1 ) then

      if ( finish(trail(qnode)+1) /= tnode ) then

        if ( tnode < stakc2(pw21+2) ) then
          planar = .false.
          return
        end if

        if ( tnode < stakc3(pw22+2) ) then
          planar = .false.
          return
        end if

        pw22 = pw22 + 2
        stakc3(pw22+1) = i
        stakc3(pw22+2) = tnode

      else

        if ( tnode < stakc3(pw22+2) .and. &
          start(stakc3(pw22+1)+1) <= descp(qnode) ) then
          planar = .false.
          return
        end if

        pw20 = pw20 + 2

        stakc1(pw20+1) = i
        stakc1(pw20+2) = qnode

      end if

    else

      if ( finish(trail(qnode) + 1) /= tnode ) then

         if ( tnode < stakc1(pw20+2) ) then
           planar = .false.
           return
         end if

         if ( tnode < stakc4(pw23+2) ) then
           planar = .false.
           return
         end if

         pw23 = pw23 + 2
         stakc4(pw23+1) = i
         stakc4(pw23+2) = tnode

        else

         if ( tnode < stakc4(pw23+2) .and. &
           start(stakc4(pw23+1)+1) <= descp(qnode) ) then
           planar = .false.
           return
         end if

         pw21 = pw21 + 2
         stakc2(pw21+1) = i
         stakc2(pw21+2) = qnode

       end if

      end if

    end do

    planar = .true.

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
!    Input, integer ( kind = 4 ) INODE(NEDGE), JNODE(NEDGE), the beginning and end
!    nodes of the edges.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) nedge

  integer ( kind = 4 ) i
  integer ( kind = 4 ) inode(nedge)
  integer ( kind = 4 ) jnode(nedge)
  character ( len = * ) title

  if ( len_trim ( title ) /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  do i = 1, nedge
    write ( *, '(i8,4x,2i8)' ) i, inode(i), jnode(i)
  end do

  return
end
subroutine graph_arc_to_star ( nnode, nedge, inode, jnode, arcfir, fwdarc )

!*****************************************************************************80
!
!! GRAPH_ARC_TO_STAR sets up the forward star representation of an undirected graph.
!
!  Modified:
!
!    04 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) NEDGE, the number of edges.
!
!    Input, integer ( kind = 4 ) INODE(NEDGE), JNODE(NEDGE); the I-th edge of
!    the graph extends from node INODE(I) to JNODE(I).
!
!    Output, integer ( kind = 4 ) ARCFIR(NNODE+1); ARCFIR(I) is the number of
!    the first edge starting at node I in the forward star representation of
!    the graph.
!
!    Output, integer ( kind = 4 ) FWDARC(2*NEDGE); FWDARC(I) is the ending node
!    of the I-th edge in the forward star representation of the graph.
!
  implicit none

  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) arcfir(nnode+1)
  integer ( kind = 4 ) fwdarc(2*nedge)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) inode(nedge)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jnode(nedge)
  integer ( kind = 4 ) k
!
!  Set up the forward star representation of the graph.
!
  k = 0

  do i = 1, nnode

    arcfir(i) = k + 1

    do j = 1, nedge

      if ( inode(j) == i ) then
        k = k + 1
        fwdarc(k) = jnode(j)
      end if

      if ( jnode(j) == i ) then
        k = k + 1
        fwdarc(k) = inode(j)
      end if

    end do

  end do

  arcfir(nnode+1) = k + 1

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
!    Input, integer ( kind = 4 ) INODE(NEDGE), JNODE(NEDGE), the beginning and
!    end nodes of the edges.
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

  if ( len_trim ( title ) /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  do i = 1, nedge
    write ( *, '(i8,4x,2i8,g14.6)' ) i, inode(i), jnode(i), wnode(i)
  end do

  return
end
subroutine graph_dist_mintr1 ( nnode, dist, ndim, itree, jtree )

!*****************************************************************************80
!
!! GRAPH_DIST_MINTR1: minimum spanning tree for a graph using a distance matrix.
!
!  Discussion:
!
!    The graph is assumed to be undirected, and connected.
!
!    Consider an undirected graph G with given edge lengths.  The minimum
!    spanning tree problem is to find a spanning tree in G such that the
!    sum of the edge lengths in the tree is minimum.
!
!  Method:
!
!    Choose any node J and let this single node J be the partially
!    constructed tree T.  Join to T an edge whose length is minimum
!    among all edges with one end in T and the other end not in T.
!    Repeat this process until T becomes a spanning tree of G.
!
!    The running time of the algorithm is O ( NNODE**2 ), where NNODE is
!    the number of nodes in the graph.
!
!  Modified:
!
!    22 July 2000
!
!  Author:
!
!    Hang Tong Lau
!
!  Reference:
!
!    Hang Tong Lau,
!    Algorithms on Graphs,
!    Tab Books, 1989,
!    QA166 L38.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, real ( kind = 8 ) DIST(NDIM,NNODE), the distance matrix.
!    DIST(I,J) is the distance from node I to node J.  DIST should be
!    symmetric, that is, DIST(I,J) = DIST(J,I).  DIST(I,I) should be 0.
!    If there is no edge between nodes I and J, then set DIST(I,J) to
!    HUGE(1.0).
!
!    Input, integer ( kind = 4 ) NDIM, the first dimension of DIST, which must
!    be at least NNODE.
!
!    Output, integer ( kind = 4 ) ITREE(NNODE-1), JTREE(NNODE-1), the pairs of nodes
!    that define the edges of the tree.
!
  implicit none

  integer ( kind = 4 ) nnode
  integer ( kind = 4 ) ndim

  real ( kind = 8 ) d
  real ( kind = 8 ) dist(ndim,nnode)
  real ( kind = 8 ) distmn
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imin
  integer ( kind = 4 ) inode
  integer ( kind = 4 ) itree(nnode-1)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jtree(nnode-1)
  integer ( kind = 4 ) mtree(nnode)

  call i4vec_indicator ( nnode-1, itree )

  mtree(1:nnode-1) = - nnode
  mtree(nnode) = 0

  do i = 1, nnode - 1

    distmn = huge ( distmn )

    do j = 1, nnode - 1

      inode = mtree(j)

      if ( inode <= 0 ) then

        d  = dist(-inode,j)

        if ( d < distmn ) then
          distmn = d
          imin = j
        end if

     end if

    end do

    mtree(imin) = - mtree(imin)

    do j = 1, nnode-1

      inode = mtree(j)

      if ( inode <= 0 ) then

        if ( dist(j,imin) < dist(j,-inode) ) then
          mtree(j) = - imin
        end if

      end if

    end do

  end do

  do i = 1, nnode-1
    jtree(i) = mtree(i)
  end do

  return
end
subroutine graph_dist_print ( dist, lda, nnode, title )

!*****************************************************************************80
!
!! GRAPH_DIST_PRINT prints a distance matrix.
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
!    Input, real ( kind = 8 ) DIST(LDA,NNODE), the distance matrix.
!    DIST(I,J) is the distance from node I to node J.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of DIST, which
!    must be at least NNODE.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) nnode

  real ( kind = 8 ) dist(lda,nnode)
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nrow
  character ( len = * ) title

  if ( len_trim ( title ) /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  ilo = 1
  ihi = nnode
  jlo = 1
  jhi = nnode
  ncol = nnode
  nrow = nnode

  call r8mat_print ( dist, ihi, ilo, jhi, jlo, lda, ncol, nrow )

  return
end
subroutine i4_swap ( i, j )

!*****************************************************************************80
!
!! I4_SWAP swaps two I4's.
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
subroutine i4vec_indicator ( n, a )

!*****************************************************************************80
!
!! I4VEC_INDICATOR sets an I4VEC to the indicator vector.
!
!  Modified:
!
!    09 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Output, integer ( kind = 4 ) A(N), the array to be initialized.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i

  do i = 1, n
    a(i) = i
  end do

  return
end
subroutine i4vec_print ( n, a, title )

!*****************************************************************************80
!
!! I4VEC_PRINT prints an I4VEC.
!
!  Modified:
!
!    16 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, integer ( kind = 4 ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  if ( title /= ' ' ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(i8,i10)' ) i, a(i)
  end do

  return
end
subroutine i4vec_reverse ( n, a )

!*****************************************************************************80
!
!! I4VEC_REVERSE reverses the elements of an I4VEC.
!
!  Example:
!
!    Input:
!
!      N = 5,
!      A = ( 11, 12, 13, 14, 15 ).
!
!    Output:
!
!      A = ( 15, 14, 13, 12, 11 ).
!
!  Modified:
!
!    26 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input/output, integer ( kind = 4 ) A(N), the array to be reversed.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i

  do i = 1, n/2
    call i4_swap ( a(i), a(n+1-i) )
  end do

  return
end
subroutine r8_swap ( x, y )

!*****************************************************************************80
!
!! R8_SWAP swaps two R8's.
!
!  Modified:
!
!    01 May 2000
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
subroutine r8mat_print ( a, ihi, ilo, jhi, jlo, lda, ncol, nrow )

!*****************************************************************************80
!
!! R8MAT_PRINT prints out a portion of an R8MAT.
!
!  Modified:
!
!    01 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A(LDA,NCOL), an NROW by NCOL matrix to be printed.
!
!    Input, integer ( kind = 4 ) IHI, ILO, the first and last rows to print.
!
!    Input, integer ( kind = 4 ) JHI, JLO, the first and last columns to print.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A.
!
!    Input, integer ( kind = 4 ) NCOL, NROW, the number of rows and columns
!    in the matrix A.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) lda
  integer ( kind = 4 ) ncol

  real ( kind = 8 ) a(lda,ncol)
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
  integer ( kind = 4 ) nrow

  write ( *, '(a)' ) ' '

  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, ncol )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i7,7x)' ) j
    end do

    write ( *, '(''Columns '',5a14)' ) ( ctemp(j2), j2 = 1, inc )
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) ' '

    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, nrow )

    do i = i2lo, i2hi

      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( a(i,j) == 0.0D+00 ) then
          ctemp(j2) = '    0.0'
        else
          write ( ctemp(j2), '(g14.6)' ) a(i,j)
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j), j = 1, inc )

    end do

  end do

  return
end
subroutine r8vec_print ( n, a, title )

!*****************************************************************************80
!
!! R8VEC_PRINT prints an R8VEC.
!
!  Modified:
!
!    16 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, real ( kind = 8 ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  if ( title /= ' ' ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(i8,g14.6)' ) i, a(i)
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
