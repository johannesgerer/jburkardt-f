subroutine cdg_code_back ( adj, nnode, code, order )

!*****************************************************************************80
!
!! CDG_CODE_BACK computes a color digraph code via backtracking.
!
!  Discussion:
!
!    The code is the "largest" order code over all possible node
!    orderings.  The lexicographic ordering is used in comparing codes.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,I) is the color of node I; otherwise, ADJ(I,J) is positive
!    if there is an edge from node I to node J.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) CODE(NNODE,NNODE), the code.
!
!    Output, integer ( kind = 4 ) ORDER(NNODE), the ordering of the nodes.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) bestcode(nnode,nnode)
  integer ( kind = 4 ) bestorder(nnode)
  integer ( kind = 4 ) code(nnode,nnode)
  logical, parameter :: debug = .false.
  integer ( kind = 4 ) index
  integer ( kind = 4 ) maxstack
  integer ( kind = 4 ) ncan(nnode)
  integer ( kind = 4 ) ncomp
  integer ( kind = 4 ) nopen
  integer ( kind = 4 ) nswap
  integer ( kind = 4 ) nstack
  integer ( kind = 4 ) order(nnode)
  integer ( kind = 4 ) result
  integer ( kind = 4 ) stack(4*nnode)

  if ( nnode <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CDG_CODE_BACK - Fatal error!'
    write ( *, '(a,i8)' ) '  NNODE = ', nnode
    stop
  end if

  maxstack = 4 * nnode
  nstack = 0
  stack(1) = 0

  ncomp = 0
  nswap = 0
!
!  Start with the identity ordering.
!
  call i4vec_indicator ( nnode, order )
!
!  Compute the corresponding code.
!
  call cdg_order_code ( adj, nnode, nnode, code, order )
!
!  Take this ordering and code as the best so far.
!
  bestorder(1:nnode) = order(1:nnode)
  bestcode(1:nnode,1:nnode) = code(1:nnode,1:nnode)
!
!  Now consider all possible orderings, and their codes.
!
  index = 0

  do

    call i4vec_backtrack ( nnode, order, index, nopen, nstack, stack, &
      maxstack, ncan )
!
!  If the backtracking routine has returned a complete candidate
!  ordering, then compute the resulting code, and see it it is better
!  then our current best.  Then go back for the next backtrack search.
!
    if ( index == 1 ) then

      call cdg_order_code ( adj, nnode, nnode, code, order )

      call cdg_code_compare ( bestcode, code, nnode, nnode, result )

      ncomp = ncomp + 1

      if ( result == - 1 ) then

        nswap = nswap + 1

        bestorder(1:nnode) = order(1:nnode)
        bestcode(1:nnode,1:nnode) = code(1:nnode,1:nnode)

      end if
!
!  If the backtracking routine has a partial reordering,
!  supply candidates for the next item in the ordering.
!
    else if ( index == 2 ) then

      call cdg_code_cand ( adj, bestcode, code, nnode, ncomp, nopen, &
        order, stack, maxstack, nstack, ncan )

    else

      exit

    end if

  end do
!
!  Once we have examined all possibilites, we are done.
!
!  Set the output ordering to the best ordering, and the output
!  code to the corresponding best code.
!
  order(1:nnode) = bestorder(1:nnode)
  code(1:nnode,1:nnode) = bestcode(1:nnode,1:nnode)

  if ( debug ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CDG_CODE_BACK:'
    write ( *, '(a,i8)' ) '  Comparisons: ', ncomp
    write ( *, '(a,i8)' ) '  Swaps:       ', nswap
  end if

  return
end
subroutine cdg_code_brute ( adj, nnode, code, order )

!*****************************************************************************80
!
!! CDG_CODE_BRUTE computes the color digraph code via brute force.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,I) is the color of node I; otherwise, ADJ(I,J) is positive
!    if there is an edge from node I to node J.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) CODE(NNODE,NNODE), the code.
!
!    Output, integer ( kind = 4 ) ORDER(NNODE), the ordering of the nodes.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) bestcode(nnode,nnode)
  integer ( kind = 4 ) bestorder(nnode)
  integer ( kind = 4 ) code(nnode,nnode)
  logical even
  logical more
  integer ( kind = 4 ) ncomp
  integer ( kind = 4 ) nswap
  integer ( kind = 4 ) order(nnode)
  integer ( kind = 4 ) result

  ncomp = 0
  nswap = 0
!
!  Start with the identity ordering.
!
  call i4vec_indicator ( nnode, order )
!
!  Compute the corresponding code.
!
  call cdg_order_code ( adj, nnode, nnode, code, order )
!
!  Take this ordering and code as the best so far.
!
  bestorder(1:nnode) = order(1:nnode)
  bestcode(1:nnode,1:nnode) = code(1:nnode,1:nnode)
!
!  Now consider all possible orderings, and their codes.
!
  more = .false.

  do

    call perm_next ( nnode, order, more, even )

    call cdg_order_code ( adj, nnode, nnode, code, order )

    call cdg_code_compare ( bestcode, code, nnode, nnode, result )

    ncomp = ncomp + 1

    if ( result == - 1 ) then

      nswap = nswap + 1

      bestorder(1:nnode) = order(1:nnode)
      bestcode(1:nnode,1:nnode) = code(1:nnode,1:nnode)

    end if

    if ( .not. more ) then
      exit
    end if

  end do
!
!  Once we have examined all possibilites, we are done.
!
!  Set the output ordering to the best ordering, and the output
!  code to the corresponding best code.
!
  order(1:nnode) = bestorder(1:nnode)
  code(1:nnode,1:nnode) = bestcode(1:nnode,1:nnode)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CDG_CODE_BRUTE:'
  write ( *, '(a,i8)' ) '  Comparisons: ', ncomp
  write ( *, '(a,i8)' ) '  Swaps:       ', nswap

  return
end
subroutine cdg_code_cand ( adj, bestcode, code, nnode, ncomp, nopen, &
  order, stack, maxstack, nstack, ncan )

!*****************************************************************************80
!
!! CDG_CODE_CAND finds candidates for a maximal color digraph code ordering.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,I) is the color of node I; otherwise, ADJ(I,J) is positive
!    if there is an edge from node I to node J.
!
!    Input, integer ( kind = 4 ) BESTCODE(NNODE,NNODE), the best code so far.
!
!    Workspace, integer CODE(NNODE,NNODE).
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input/output, integer ( kind = 4 ) NCOMP, the number of code comparisons.
!    This routine updates NCOMP by 1 each time the routine is called.
!
!    Input, integer ( kind = 4 ) NOPEN, identifies the first open position in ORDER.
!    1 <= NOPEN <= NNODE.
!
!    Input, integer ( kind = 4 ) ORDER(NNODE), contains in entries 1 through NOPEN-1
!    the elements of the current partial list.
!
!    Input/output, integer ( kind = 4 ) STACK(MAXSTACK), used to store the candidates.
!
!    Input, integer ( kind = 4 ) MAXSTACK, the maximum size of the STACK array.
!    A value of NNODE should be sufficient.
!
!    Input/output, integer ( kind = 4 ) NSTACK, the current length of the stack.
!    On output, NSTACK has been increased by the number of
!    candidates found.
!
!    Input/output, integer ( kind = 4 ) NCAN(NNODE), counts the number of candidates
!    for each position.
!
  implicit none

  integer ( kind = 4 ) maxstack
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) bestcode(nnode,nnode)
  integer ( kind = 4 ) code(nnode,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ifree(nnode)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) maxcolor
  integer ( kind = 4 ) ncan(nnode)
  integer ( kind = 4 ) ncomp
  integer ( kind = 4 ) nfree
  integer ( kind = 4 ) ni
  integer ( kind = 4 ) nj
  integer ( kind = 4 ) nopen
  integer ( kind = 4 ) nstack
  integer ( kind = 4 ) order(nnode)
  integer ( kind = 4 ) result
  integer ( kind = 4 ) stack(maxstack)

  if ( nopen < 1 .or. nnode < nopen ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CDG_CODE_CAND - Fatal error 1!'
    write ( *, '(a)' ) '  1 <= NOPEN <= NNODE should be true, but'
    write ( *, '(a,i8)' ) '  NOPEN = ', nopen
    write ( *, '(a,i8)' ) '  NNODE = ', nnode
    stop
  end if
!
!  Start with no candidates.
!
  ncan(nopen) = 0
!
!  If we have fixed at least one entry of the list,
!
!    Compute the partial code;
!
!    Compare the partial code with the corresponding
!    part of the the code for the best ordering so far;
!
!    If the current incomplete ordering is actually LESS than the
!    current best, then bail out now, with zero candidates.
!
  if ( 1 < nopen ) then

    call cdg_order_code ( adj, nnode, nopen-1, code, order )

    call cdg_code_compare ( bestcode, code, nnode, nopen-1, result )

    ncomp = ncomp + 1

    if ( result == + 1 ) then
      ncan(nopen) = 0
      return
    end if

  end if
!
!  Get a list of those nodes which have not been used yet.
!
  nfree = nnode + 1 - nopen
  call perm_free ( order, nopen-1, ifree, nfree )
!
!  Our preferred candidates will be:
!
  do i = 1, nopen-1

    ncan(nopen) = 0

    ni = order(i)
!
!    * for the LOWEST ordered node possible, all unordered OUT neighbors,
!
    do j = 1, nfree

      nj = ifree(j)

      if ( adj(ni,nj) /= 0 ) then

        ncan(nopen) = ncan(nopen) + 1
        nstack = nstack + 1

        if ( maxstack < nstack ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'CDG_CODE_CAND - Fatal error 4!'
          write ( *, '(a)' ) '  MAXSTACK < NSTACK!'
          write ( *, '(a,i8)' ) '  NSTACK =   ', nstack
          write ( *, '(a,i8)' ) '  MAXSTACK = ', maxstack
          stop
        end if

        stack(nstack) = nj

      end if

    end do

    if ( 0 < ncan(nopen) ) then
      return
    end if
!
!    * for the LOWEST ordered node possible, all unordered IN neighbors,
!
    do j = 1, nfree

      nj = ifree(j)

      if ( adj(nj,ni) /= 0 ) then

        ncan(nopen) = ncan(nopen) + 1
        nstack = nstack + 1

        if ( maxstack < nstack ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'CDG_CODE_CAND - Fatal error 4!'
          write ( *, '(a)' ) '  MAXSTACK < NSTACK!'
          write ( *, '(a,i8)' ) '  NSTACK =   ', nstack
          write ( *, '(a,i8)' ) '  MAXSTACK = ', maxstack
          stop
        end if

        stack(nstack) = nj

      end if
    end do

    if ( 0 < ncan(nopen) ) then
      return
    end if

  end do
!
!  NO unordered nodes are connected in any way to ordered nodes.
!  This can happen in two ways:
!
!  * NOPEN = 1; (the list of used nodes is empty)
!  * The graph is disconnected;
!
!  In either case, we must now consider ALL free nodes.
!
!  Compute the maximal color.
!
  maxcolor = 0

  do i = 1, nfree
    ni = ifree(i)
    maxcolor = max ( maxcolor, adj(ni,ni) )
  end do
!
!  Take as candidates every node of color MAXCOLOR.
!
!  We could thin the list down, by looking ahead, and only taking
!  candidates of MAXCOLOR who also happen to have at least one free
!  out neighbor, and so on.
!
  ncan(nopen) = 0

  do i = 1, nfree

    ni = ifree(i)

    if ( adj(ni,ni) == maxcolor ) then

      ncan(nopen) = ncan(nopen) + 1
      nstack = nstack + 1

      if ( maxstack < nstack ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CDG_CODE_CAND - Fatal error 6!'
        write ( *, '(a)' ) '  MAXSTACK < NSTACK!'
        write ( *, '(a,i8)' ) '  NSTACK =   ', nstack
        write ( *, '(a,i8)' ) '  MAXSTACK = ', maxstack
        stop
      end if

      stack(nstack) = ni

    end if

  end do
!
!  This should never happen:
!
  if ( ncan(nopen) == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CDG_CODE_CAND - Fatal error 7!'
    write ( *, '(a)' ) '  No candidates, but there gotta be some!'
    stop
  end if

  return
end
subroutine cdg_code_compare ( code1, code2, nnode, npart, result )

!*****************************************************************************80
!
!! CDG_CODE_COMPARE compares two (partial) color graph codes.
!
!  Discussion:
!
!    CODE1 = CODE2 if every digit of both codes is equal.
!
!    Otherwise, traverse the entries in a funny diagonal way, suggested
!    by this diagram for the first 16 entries:
!
!       1  2  5 10
!       3  4  7 12
!       6  8  9 14
!      11 13 15 16
!
!    As we do that, we examine the corresponding digits of the two codes.
!    For the first entry, (I,J), where the two codes differ, we say:
!
!      if ( CODE1(I,J) < CODE2(I,J) ) then we say
!        CODE1 < CODE2
!      else
!        CODE2 < CODE1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) CODE1(NNODE,NNODE), CODE2(NNODE,NNODE),
!    two codes to be compared.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) NPART, specifies the portion of the codes
!    to compare.  NPART should be between 1 and NNODE.
!
!    If NPART = NNODE, then the full codes are compared.
!
!    If NPART < NNODE, then only entries corresponding to I and J
!    both less than or equal to NPART will be compared.
!
!    Output, integer ( kind = 4 ) RESULT:
!    -1, CODE1 < CODE2;
!     0, CODE1 = CODE2;
!    +1, CODE2 < CODE1.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) code1(nnode,nnode)
  integer ( kind = 4 ) code2(nnode,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) npart
  integer ( kind = 4 ) result

  do j = 1, npart

    do i = 1, j - 1

      if ( code1(i,j) < code2(i,j) ) then

        result = - 1
        return

      else if ( code2(i,j) < code1(i,j) ) then

        result = + 1
        return

      else if ( code1(j,i) < code2(j,i) ) then

        result = - 1
        return

      else if ( code2(j,i) < code1(j,i) ) then

        result = + 1
        return

      end if

    end do

    if ( code1(j,j) < code2(j,j) ) then

      result = - 1
      return

    else if ( code2(j,j) < code1(j,j) ) then

      result = + 1
      return

    end if

  end do

  result = 0

  return
end
subroutine cdg_code_print ( nnode, code, title )

!*****************************************************************************80
!
!! CDG_CODE_PRINT prints a color digraph code.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) CODE(NNODE,NNODE), the code.
!
!    Input, character ( len = * ) TITLE, a title to be printed.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) ck
  integer ( kind = 4 ) code(nnode,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  character ( len = 80 ) string
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  do i = 1, nnode

    do j = 1, nnode

      ck = code(i,j)
      if ( 0 <= ck .and. ck <= 9 ) then
        string(j:j) = char ( 48 + ck )
      else
        string(j:j) = '*'
      end if

    end do

    write ( *, '(2x,i4,2x,a)' ) i, string(1:nnode)

  end do

  return
end
subroutine cdg_color_count ( adj, nnode, mcolor, ncolor )

!*****************************************************************************80
!
!! CDG_COLOR_COUNT counts the number of colors in a color digraph.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,I) is the color of node I; otherwise, ADJ(I,J) is positive
!    if there is an edge from node I to node J.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) MCOLOR, the maximum color index.
!
!    Output, integer ( kind = 4 ) NCOLOR, the number of colors.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) colors(nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) mcolor
  integer ( kind = 4 ) ncolor

  mcolor = 0
  do i = 1, nnode
    mcolor = max ( mcolor, adj(i,i) )
  end do

  do i = 1, nnode
    colors(i) = adj(i,i)
  end do

  call i4vec_sort_heap_d ( nnode, colors )

  call i4vec_sorted_unique_count ( nnode, colors, ncolor )

  return
end
subroutine cdg_color_sequence ( adj, nnode, seq )

!*****************************************************************************80
!
!! CDG_COLOR_SEQUENCE computes the color sequence of a color digraph.
!
!  Discussion:
!
!    The color sequence of a color digraph is constructed by computing the
!    color of each node, and then ordering these values in decreasing order.
!
!    If two color digraphs are isomorphic, they must have the same color
!    sequence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,I) is the color of node I; otherwise, ADJ(I,J) is positive
!    if there is an edge from node I to node J.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) SEQ(NNODE), the color sequence.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) seq(nnode)

  do i = 1, nnode
    seq(i) = adj(i,i)
  end do

  call i4vec_sort_heap_d ( nnode, seq )

  return
end
subroutine cdg_compare ( adj1, nnode1, adj2, nnode2, order1, order2, result )

!*****************************************************************************80
!
!! CDG_COMPARE determines if color digraphs G1 and G2 are isomorphic.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ1(NNODE1,NNODE1), the adjacency information for G1.
!    ADJ1(I,I) is the color of node I; otherwise, ADJ1(I,J) is positive
!    if there is an edge from node I to node J.
!
!    Input, integer ( kind = 4 ) NNODE1, the number of nodes in G1.
!
!    Input, integer ( kind = 4 ) ADJ2(NNODE2,NNODE2), the adjacency information for G2.
!    ADJ2(I,I) is the color of node I; otherwise, ADJ2(I,J) is positive
!    if there is an edge from node I to node J.
!
!    Input, integer ( kind = 4 ) NNODE2, the number of nodes in G2.
!
!    Output, integer ( kind = 4 ) RESULT, is 0 if G1 and G2 are isomorphic,
!    -I if G1 < G2 for test #I, and
!    +I if G2 < G1 for test #I.
!
!    Output, integer ( kind = 4 ) ORDER1(NNODE1), ORDER2(NNODE2).  If RESULT = 0, then
!    ORDER1 and ORDER2 are reorderings of the nodes of G1 and
!    G2 which exhibit the isomorphism.
!
  implicit none

  integer ( kind = 4 ) nnode1
  integer ( kind = 4 ) nnode2

  integer ( kind = 4 ) adj1(nnode1,nnode1)
  integer ( kind = 4 ) adj2(nnode2,nnode2)
  integer ( kind = 4 ) code1(nnode1,nnode1)
  integer ( kind = 4 ) code2(nnode2,nnode2)
  integer ( kind = 4 ) in_seq1(nnode1)
  integer ( kind = 4 ) in_seq2(nnode2)
  integer ( kind = 4 ) mcolor1
  integer ( kind = 4 ) mcolor2
  integer ( kind = 4 ) ncolor1
  integer ( kind = 4 ) ncolor2
  integer ( kind = 4 ) nedge1
  integer ( kind = 4 ) nedge2
  integer ( kind = 4 ) order1(nnode1)
  integer ( kind = 4 ) order2(nnode2)
  integer ( kind = 4 ) out_seq1(nnode1)
  integer ( kind = 4 ) out_seq2(nnode2)
  integer ( kind = 4 ) result
!
!  Test 1: Count the nodes.
!
  if ( nnode1 < nnode2 ) then
    result = - 1
    return
  else if ( nnode2 < nnode1 ) then
    result = + 1
    return
  end if
!
!  Test 2: Count the edges.
!
  call cdg_edge_count ( adj1, nnode1, nedge1 )

  call cdg_edge_count ( adj2, nnode2, nedge2 )

  if ( nedge1 < nedge2 ) then
    result = - 2
    return
  else if ( nedge2 < nedge1 ) then
    result = + 2
    return
  end if
!
!  Tests 3 and 4: Count the colors, and note the maximum color.
!
  call cdg_color_count ( adj1, nnode1, mcolor1, ncolor1 )

  call cdg_color_count ( adj2, nnode2, mcolor2, ncolor2 )

  if ( ncolor1 < ncolor2 ) then
    result = - 3
    return
  else if ( ncolor2 < ncolor1 ) then
    result = + 3
    return
  end if

  if ( mcolor1 < mcolor2 ) then
    result = - 4
    return
  else if ( mcolor2 < mcolor1 ) then
    result = + 4
    return
  end if
!
!  Test 5: Compare the outdegree sequences.
!
  call cdg_degree_seq ( adj1, nnode1, in_seq1, out_seq1 )

  call cdg_degree_seq ( adj2, nnode2, in_seq2, out_seq2 )

  call i4vec_compare ( nnode1, out_seq1, out_seq2, result )

  if ( result < 0 ) then
    result = - 5
    return
  else if ( 0 < result ) then
    result = + 5
    return
  end if
!
!  Test 6: Compare the indegree sequences.
!
  call i4vec_compare ( nnode1, in_seq1, in_seq2, result )

  if ( result < 0 ) then
    result = - 6
    return
  else if ( 0 < result ) then
    result = + 6
    return
  end if
!
!  Test 7: Compare the color sequences.
!
  call cdg_color_sequence ( adj1, nnode1, in_seq1 )

  call cdg_color_sequence ( adj2, nnode2, in_seq2 )

  call i4vec_compare ( nnode1, in_seq1, in_seq2, result )

  if ( result < 0 ) then
    result = - 7
    return
  else if ( 0 < result ) then
    result = + 7
    return
  end if
!
!  Test 8: Compare the codes.
!
  call cdg_code_back ( adj1, nnode1, code1, order1 )

  call cdg_code_back ( adj2, nnode2, code2, order2 )

  call cdg_code_compare ( code1, code2, nnode1, nnode1, result )

  if ( result < 0 ) then
    result = - 8
    return
  else if ( 0 < result ) then
    result = + 8
    return
  end if

  result = 0

  return
end
subroutine cdg_degree ( adj, nnode, indegree, outdegree )

!*****************************************************************************80
!
!! CDG_DEGREE computes the indegree and outdegree of each node.
!
!  Discussion:
!
!    The indegree of a node is the number of directed edges that
!    end at the node.
!
!    The outdegree of a node is the number of directed edges that
!    begin at the node.
!
!    The sum of the indegrees and outdegrees of all the nodes is twice
!    the number of edges.
!
!    The generalized case, where ADJ(I,J) can be greater than 1, indicating
!    the existence of 2 or more distinct edges from node I to node J,
!    will be properly handled by this routine.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information for graph 1.
!    ADJ(I,I) is the color of node I; otherwise, ADJ(I,J) is positive
!    if there is an edge from node I to node J.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) INDEGREE(NNODE), OUTDEGREE(NNODE),
!    the indegree and outdegree of the nodes.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indegree(nnode)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) outdegree(nnode)

  indegree(1:nnode) = 0
  outdegree(1:nnode) = 0

  do i = 1, nnode
    do j = 1, nnode
      if ( i /= j ) then
        if ( adj(i,j) /= 0 ) then
          outdegree(i) = outdegree(i) + adj(i,j)
          indegree(j) = indegree(j) + adj(i,j)
        end if
      end if
    end do
  end do

  return
end
subroutine cdg_degree_seq ( adj, nnode, in_seq, out_seq )

!*****************************************************************************80
!
!! CDG_DEGREE_SEQ computes the degree sequence of a color digraph.
!
!  Discussion:
!
!    The directed degree sequence of a graph is the sequence of indegrees
!    and the sequence of outdegrees, arranged to correspond to nodes of
!    successively decreasing total degree.  For nodes of equal degree, those
!    of higher outdegree take precedence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,I) is the color of node I; otherwise, ADJ(I,J) is positive
!    if there is an edge from node I to node J.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) IN_SEQ(NNODE), OUT_SEQ(NNODE), the degree sequence.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) in_seq(nnode)
  integer ( kind = 4 ) out_seq(nnode)

  call cdg_degree ( adj, nnode, in_seq, out_seq )

  call i4vec2_sort_d ( nnode, out_seq, in_seq )

  return
end
subroutine cdg_edge_count ( adj, nnode, nedge )

!*****************************************************************************80
!
!! CDG_EDGE_COUNT counts the number of edges in a color digraph.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,I) is the color of node I; otherwise, ADJ(I,J) is positive
!    if there is an edge from node I to node J.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) NEDGE, the number of edges.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nedge

  nedge = 0

  do i = 1, nnode
    do j = 1, nnode

      if ( i /= j ) then
        nedge = nedge + adj(i,j)
      end if

    end do
  end do

  return
end
subroutine cdg_example_cube ( adj, nnode )

!*****************************************************************************80
!
!! CDG_EXAMPLE_CUBE sets up the cube color digraph.
!
!  Diagram:
!
!
!    8B----<-----3B
!    |\          /|\
!    | A        V | |
!    |  \      /  | |
!    |  4R-->-7R  | |
!    |   |     |  | |
!    A   A     V  V A
!    |   |     |  | |
!    |   5B-<-2G  | |
!    |  /      \  | |
!    | A        A | |
!    |/          \|/
!    1G----->----6B
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,I) is the color of node I; otherwise, ADJ(I,J) is positive
!    if there is an edge from node I to node J.
!
!    Output, integer ( kind = 4 ) NNODE, the number of nodes.
!
  implicit none

  integer ( kind = 4 ) adj(8,8)
  integer ( kind = 4 ), parameter :: BLUE = 1
  integer ( kind = 4 ), parameter :: GREEN = 2
  integer ( kind = 4 ) nnode
  integer ( kind = 4 ), parameter :: RED = 3

  nnode = 8

  adj(1:nnode,1:nnode) = 0

  adj(1,1) = GREEN
  adj(1,5) = 1
  adj(1,6) = 1
  adj(1,8) = 1

  adj(2,2) = GREEN
  adj(2,5) = 1

  adj(3,3) = BLUE
  adj(3,6) = 1
  adj(3,7) = 1
  adj(3,8) = 1

  adj(4,4) = RED
  adj(4,7) = 1
  adj(4,8) = 1

  adj(5,5) = BLUE
  adj(5,4) = 1

  adj(6,6) = BLUE
  adj(6,2) = 1
  adj(6,3) = 1

  adj(7,7) = RED
  adj(7,2) = 1

  adj(8,8) = BLUE

  return
end
subroutine cdg_example_octo ( example, adj, nnode, seed )

!*****************************************************************************80
!
!! CDG_EXAMPLE_OCTO sets up an 8 node example color digraph.
!
!  Diagram:
!
!      1---2
!     /|   |\
!    8-+---+-3
!    | |   | |
!    7-+---+-4
!     \|   |/
!      6---5
!
!     Graph "A"
!
!    There are 8 graphs to choose from.  They are all on 8 nodes.  The first
!    5 have degree 3 at every node.  Graphs 6 and 7 have degree 5 at every
!    node.  Graph 8 is disconnected.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) EXAMPLE, the index of the example to choose.
!    1 <= EXAMPLE <= 65.
!
!    Output, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information for the graph.
!    ADJ(I,I) is the color of node I.
!    ADJ(I,J) is 1 if nodes I and J are adjacent and 0 otherwise.
!
!    Output, integer ( kind = 4 ) NNODE, the number of nodes, which should be 8.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
  implicit none

  integer ( kind = 4 ) adj(8,8)
  integer ( kind = 4 ), parameter :: BLUE = 1
  integer ( kind = 4 ) example
  integer ( kind = 4 ), parameter :: GREEN = 2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) j
  integer ( kind = 4 ) msave
  integer ( kind = 4 ) nnode
  integer ( kind = 4 ) nsave
  integer ( kind = 4 ), parameter :: RED = 3
  integer ( kind = 4 ) seed
  integer ( kind = 4 ), parameter :: YELLOW = 4

  if ( example <= 0 ) then
    nsave = i4_uniform ( 1, 13, seed )
    msave = i4_uniform ( 1, 5, seed )
  else
    nnode = mod ( example - 1, 65 ) + 1
    msave = ( example - 1 ) / 13 + 1
    nsave = mod ( example - 1, 13 ) + 1
  end if

  nnode = 8

  adj(1:nnode,1:nnode) = 0

  do i = 1, nnode
    j = i + 1
    if ( nnode < j ) then
      j = j - nnode
    end if

    adj(i,j) = 1

  end do
!
!  Underlying graph 1.
!
  if ( nsave == 1 ) then

    adj(1,6) = 1
    adj(2,5) = 1
    adj(3,8) = 1
    adj(4,7) = 1

  else if ( nsave == 2 ) then

    adj(1,6) = 1
    adj(5,2) = 1
    adj(3,8) = 1
    adj(7,4) = 1
!
!  Underlying graph 2.
!  Digraphs 3 and 4 have different indegree/outdegree sequences.
!
  else if ( nsave == 3 ) then

    adj(1,6) = 1
    adj(6,1) = 1
    adj(2,8) = 1
    adj(8,2) = 1
    adj(3,5) = 1
    adj(5,3) = 1
    adj(4,7) = 1
    adj(7,4) = 1

  else if ( nsave == 4 ) then

    adj(1,6) = 1
    adj(2,8) = 1
    adj(3,5) = 1
    adj(4,7) = 1
!
!  Underlying graph 3
!  Digraphs 5 and 6 have the same indegree/outdegree sequences.
!
  else if ( nsave == 5 ) then

    adj(1,5) = 1
    adj(2,6) = 1
    adj(3,7) = 1
    adj(4,8) = 1

  else if ( nsave == 6 ) then

    adj(1:nnode,1:nnode) = 0

    adj(1,8) = 1
    adj(1,5) = 1
    adj(2,1) = 1
    adj(2,3) = 1
    adj(3,4) = 1
    adj(3,7) = 1
    adj(4,5) = 1
    adj(4,8) = 1
    adj(5,6) = 1
    adj(6,2) = 1
    adj(7,6) = 1
    adj(8,7) = 1
!
!  Underlying graph 4
!
  else if ( nsave == 7 ) then

    adj(3,1) = 1
    adj(4,2) = 1
    adj(5,7) = 1
    adj(6,8) = 1

  else if ( nsave == 8 ) then

    adj(3,1) = 1
    adj(4,2) = 1
    adj(5,7) = 1
    adj(8,6) = 1
!
!  Underlying graph 5
!
  else if ( nsave == 9 ) then

    adj(1,4) = 1
    adj(2,6) = 1
    adj(8,3) = 1

    adj(5,7) = 1
    adj(7,5) = 1

  else if ( nsave == 10 ) then

    adj(1,4) = 1
    adj(2,6) = 1
    adj(3,8) = 1

    adj(5,7) = 1
    adj(7,5) = 1
!
!  Underlying graph 6
!
  else if ( nsave == 11 ) then

    adj(1,4) = 1
    adj(1,5) = 1
    adj(1,6) = 1

    adj(2,5) = 1
    adj(2,6) = 1
    adj(2,7) = 1

    adj(3,6) = 1
    adj(3,7) = 1
    adj(3,8) = 1

    adj(4,7) = 1
    adj(4,8) = 1

    adj(5,8) = 1
!
!  Underlying graph 7
!
  else if ( nsave == 12 ) then

    adj(1,3) = 1
    adj(1,5) = 1
    adj(1,7) = 1

    adj(2,4) = 1
    adj(2,6) = 1
    adj(2,8) = 1

    adj(3,5) = 1
    adj(3,7) = 1

    adj(4,6) = 1
    adj(4,8) = 1

    adj(5,7) = 1

    adj(6,8) = 1
!
!  Underlying graph 8.
!
  else if ( nsave == 13 ) then

    adj(1,2) = 1
    adj(3,1) = 1
    adj(2,3) = 1
    adj(3,4) = 1
    adj(5,6) = 1
    adj(6,5) = 1
    adj(5,7) = 1
    adj(6,7) = 1

  end if

  if ( msave == 1 ) then

    adj(1,1) = RED
    adj(2,2) = RED
    adj(3,3) = RED
    adj(4,4) = BLUE
    adj(5,5) = BLUE
    adj(6,6) = BLUE
    adj(7,7) = GREEN
    adj(8,8) = GREEN

  else if ( msave == 2 ) then

    adj(1,1) = RED
    adj(2,2) = RED
    adj(3,3) = RED
    adj(4,4) = BLUE
    adj(5,5) = BLUE
    adj(6,6) = BLUE
    adj(7,7) = GREEN
    adj(8,8) = YELLOW

  else if ( msave == 3 ) then

    adj(1,1) = RED
    adj(2,2) = RED
    adj(3,3) = RED
    adj(4,4) = BLUE
    adj(5,5) = BLUE
    adj(6,6) = BLUE
    adj(7,7) = YELLOW
    adj(8,8) = YELLOW

  else if ( msave == 4 ) then

    adj(1,1) = RED
    adj(2,2) = RED
    adj(3,3) = RED
    adj(4,4) = BLUE
    adj(5,5) = BLUE
    adj(6,6) = GREEN
    adj(7,7) = GREEN
    adj(8,8) = GREEN

  else if ( msave == 5 ) then

    adj(1,1) = RED
    adj(2,2) = BLUE
    adj(3,3) = RED
    adj(4,4) = GREEN
    adj(5,5) = BLUE
    adj(6,6) = RED
    adj(7,7) = BLUE
    adj(8,8) = GREEN

  end if
!
!  Now permute the graph.
!
  call i4mat_perm_random ( nnode, adj, seed )

  return
end
subroutine cdg_order_code ( adj, nnode, npart, code, order )

!*****************************************************************************80
!
!! CDG_ORDER_CODE returns the color digraph code for a specific node ordering.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,I) is the color of node I; otherwise, ADJ(I,J) is positive
!    if there is an edge from node I to node J.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) NPART, the number of nodes to consider.
!    NPART should be between 1 and NNODE.
!
!    If NPART is NNODE, then the usual code is returned.
!
!    If NPART is less than NNODE, then the code is computed as
!    though only NPART nodes existed, namely, those specified in the
!    first NPART entries of order.  This option is provided so that
!    the routine can compute the portion of a code determined
!    by an incomplete ordering of the nodes.
!
!    Output, integer ( kind = 4 ) CODE(NNODE,NNODE), the code for this ordering.
!
!    Input, integer ( kind = 4 ) ORDER(NNODE), the ordering of the nodes.  ORDER(1)
!    is the first node, and so on.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) code(nnode,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ni
  integer ( kind = 4 ) nj
  integer ( kind = 4 ) npart
  integer ( kind = 4 ) order(nnode)

  do i = 1, nnode

    if ( i <= npart ) then

      ni = order(i)

      if ( order(i) < 1 .or. nnode < order(i) ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CDG_ORDER_CODE - Fatal error!'
        write ( *, '(a)' ) '  ORDER is not a proper permutation.'
        stop
      end if

    else
      ni = 0
    end if

    do j = 1, nnode

      if ( j <= npart ) then

        nj = order(j)

        if ( order(j) < 1 .or. nnode < order(j) ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'CDG_ORDER_CODE - Fatal error!'
          write ( *, '(a)' ) '  ORDER is not a proper permutation.'
          stop
        end if

      else
        nj = 0
      end if

      if ( ni == 0 .or. nj == 0 ) then

        code(i,j) = 0

      else

        code(i,j) = adj(ni,nj)

      end if

    end do
  end do

  return
end
subroutine cdg_print ( adj, nnode, title )

!*****************************************************************************80
!
!! CDG_PRINT prints out the adjacency matrix of a color digraph.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,I) is the color of node I; otherwise, ADJ(I,J) is positive
!    if there is an edge from node I to node J.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, character ( len = * ) TITLE, a title to be printed.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  character ( len = 80 ) string
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  do i = 1, nnode

    do j = 1, nnode

      k = (j-1) * 3 + 1
      write ( string(k:k+2), '(i3)' ) adj(i,j)

    end do

    write ( *, '(2x,a)' ) string(1:3*nnode)

  end do

  return
end
subroutine cdg_random ( adj, nnode, ncolor, nedge, seed )

!*****************************************************************************80
!
!! CDG_RANDOM generates a random color graph.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 October 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,I) is the color of node I; otherwise, ADJ(I,J) is positive
!    if there is an edge from node I to node J.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) NCOLOR, the number of colors.
!    Each node is assumed to have an associated color, between 1 and NCOLOR,
!    which will be chosen at random.
!
!    Input, integer ( kind = 4 ) NEDGE, the number of edges, which must be between
!    0 and NNODE*(NNODE-1).
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
  implicit none

  integer ( kind = 4 ) ncolor
  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) icolor
  integer ( kind = 4 ) iwork(nedge)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) maxedge
  integer ( kind = 4 ) perm(ncolor)
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) subset(ncolor)

  if ( nnode <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CDG_RANDOM - Fatal error!'
    write ( *, '(a,i8)' ) '  NNODE = ', nedge
    write ( *, '(a)' ) '  but NNODE must be at least 1.'
    stop
  end if

  maxedge = nnode * ( nnode - 1 )

  if ( nedge < 0 .or. maxedge < nedge ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CDG_RANDOM - Fatal error!'
    write ( *, '(a,i8)' ) '  NEDGE = ', nedge
    write ( *, '(a)') '  but NEDGE must be at least 0, and '
    write ( *, '(a,i8)' ) '  no more than ', maxedge
    stop
  end if
!
!  Start with no edges, no colors.
!
  adj(1:nnode,1:nnode) = 0
!
!  Choose the colors.
!
  call ksub_random ( nnode, ncolor, subset, seed )

  call perm_random ( ncolor, perm, seed )

  do icolor = 1, ncolor
    i = subset(perm(icolor))
    adj(i,i) = icolor
  end do

  do i = 1, nnode
    if ( adj(i,i) == 0 ) then
      adj(i,i) = i4_uniform ( 1, ncolor, seed )
    end if
  end do
!
!  Pick a random NEDGE subset.
!
  call ksub_random ( maxedge, nedge, iwork, seed )
!
!  Mark the potential edges that were chosen.
!
  k = 0
  l = 1

  do i = 1, nnode
    do j = 1, nnode

      if ( i /= j ) then

        k = k + 1
        if ( l <= nedge ) then

          if ( k == iwork(l) ) then
            adj(i,j) = 1
            l = l + 1
          end if

        end if

      end if

    end do
  end do

  return
end
subroutine cdmg_adj_max_max ( adj, nnode, adj_max_max )

!*****************************************************************************80
!
!! CDMG_ADJ_MAX_MAX computes the adjacency maximum maximum of a color dimultigraph.
!
!  Discussion:
!
!    The adjacency maximum maximum of a color dimultigraph may be constructed
!    by computing the maximum entry of the off diagonal entries of the
!    adjacency matrix,
!
!  Example:
!
!    ADJ =
!       3 1 2 3
!       1 9 2 0
!       2 2 2 1
!       3 0 1 7
!
!    ADJ_MAX_MAX = 3
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,I) is the color of node I; otherwise, ADJ(I,J) is
!    the number of edges from node I to node J.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) ADJ_MAX_MAX, the adjacency maximum maximum.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) adj_max_max
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  adj_max_max = 0
  do i = 1, nnode
    do j = 1, nnode
      if ( i /= j ) then
        adj_max_max = max ( adj_max_max, adj(i,j) )
      end if
    end do
  end do

  return
end
subroutine cdmg_adj_max_seq ( adj, nnode, adj_max_seq )

!*****************************************************************************80
!
!! CDMG_ADJ_MAX_SEQ computes the adjacency maximum sequence of a color dimultigraph.
!
!  Discussion:
!
!    The adjacency maximum sequence of a color dimultigraph may be
!    constructed by computing the maximum entry of each row of the
!    off diagonal elements of the adjacency matrix, and then sorting
!    these values in descending order.
!
!  Example:
!
!    ADJ =
!       9 1 2 3
!       1 8 2 0
!       2 2 3 1
!       3 0 1 6
!
!    ADJ_MAX_SEQ =
!
!       3
!       3
!       2
!       2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,I) is the color of node I; otherwise, ADJ(I,J) is
!    the number of edges from node I to node J.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) ADJ_MAX_SEQ(NNODE), the adjacency maximum sequence.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) adj_max_seq(nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
!
!  Copy the adjacency matrix.
!
  do i = 1, nnode
    k = 0
    do j = 1, nnode
      if ( i /= j ) then
        k = max ( k, adj(i,j) )
      end if
    end do
    adj_max_seq(i) = k
  end do
!
!  Sort the elements.
!
  call i4vec_sort_heap_d ( nnode, adj_max_seq )

  return
end
subroutine cdmg_adj_seq_u ( adj, nnode, adj_seq )

!*****************************************************************************80
!
!! CDMG_ADJ_SEQ_U computes the unweighted adjacency sequence of a color dimultigraph.
!
!  Discussion:
!
!    The unweighted adjacency sequence of a color dimultigraph may be
!    constructed by zeroing out the diagonal entries, replacing each nonzero
!    off diagonal entry by 1, sorting the entries of each row in descending
!    order, and then sorting the rows themselves in descending order.
!
!  Example:
!
!    ADJ =
!       5 1 2 3
!       1 7 2 0
!       2 2 8 1
!       3 0 1 9
!
!    ADJ_SEQ =
!
!       1 1 1 0
!       1 1 1 0
!       1 1 0 0
!       1 1 0 0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,I) is the color of node I; otherwise, ADJ(I,J) is
!    the number of edges from node I to node J.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) ADJ_SEQ(NNODE,NNODE), the unweighted adjacency sequence.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) adj_seq(nnode,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
!
!  Copy the adjacency matrix.
!
  do i = 1, nnode
    do j = 1, nnode

      if ( i == j ) then
        adj_seq(i,j) = 0
      else if ( adj(i,j) == 0 ) then
        adj_seq(i,j) = 0
      else
        adj_seq(i,j) = 1
      end if

    end do
  end do
!
!  Sort the elements of each row.
!
  call i4row_sort2_d ( nnode, nnode, adj_seq )
!
!  Sort the rows of the matrix.
!
  call i4row_sort_d ( nnode, nnode, adj_seq )

  return
end
subroutine cdmg_adj_seq_w ( adj, nnode, adj_seq )

!*****************************************************************************80
!
!! CDMG_ADJ_SEQ_W computes the weighted adjacency sequence of a color dimultigraph.
!
!  Discussion:
!
!    The adjacency sequence of a color dimultigraph may be constructed by
!    zeroing out the diagonal entries, sorting the entries of each row of the
!    adjacency matrix in descending order, and then sorting the rows
!    themselves in descending order.
!
!  Example:
!
!    ADJ =
!       8 1 2 3
!       1 7 2 0
!       2 2 5 1
!       3 0 1 6
!
!    ADJ_SEQ =
!
!       3 2 1 0
!       3 1 0 0
!       2 2 1 0
!       2 1 0 0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,I) is the color of node I; otherwise, ADJ(I,J) is
!    the number of edges from node I to node J.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) ADJ_SEQ(NNODE,NNODE), the adjacency sequence.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) adj_seq(nnode,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
!
!  Copy the adjacency matrix.
!
  do i = 1, nnode
    do j = 1, nnode
      if ( i == j ) then
        adj_seq(i,j) = 0
      else
        adj_seq(i,j) = adj(i,j)
      end if
    end do
  end do
!
!  Sort the elements of each row.
!
  call i4row_sort2_d ( nnode, nnode, adj_seq )
!
!  Sort the rows of the matrix.
!
  call i4row_sort_d ( nnode, nnode, adj_seq )

  return
end
subroutine cdmg_code_back ( adj, nnode, code, order )

!*****************************************************************************80
!
!! CDMG_CODE_BACK computes a color dimultigraph code via backtracking.
!
!  Discussion:
!
!    The code is the "largest" order code over all possible node
!    orderings.  The lexicographic ordering is used in comparing codes.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,I) is the color of node I; otherwise, ADJ(I,J) is
!    the number of edges from node I to node J.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) CODE(NNODE,NNODE), the code.
!
!    Output, integer ( kind = 4 ) ORDER(NNODE), the ordering of the nodes.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) bestcode(nnode,nnode)
  integer ( kind = 4 ) bestorder(nnode)
  integer ( kind = 4 ) code(nnode,nnode)
  logical, parameter :: debug = .false.
  integer ( kind = 4 ) index
  integer ( kind = 4 ) maxstack
  integer ( kind = 4 ) ncan(nnode)
  integer ( kind = 4 ) ncomp
  integer ( kind = 4 ) nopen
  integer ( kind = 4 ) nswap
  integer ( kind = 4 ) nstack
  integer ( kind = 4 ) order(nnode)
  integer ( kind = 4 ) result
  integer ( kind = 4 ) stack(4*nnode)

  if ( nnode <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CDMG_CODE_BACK - Fatal error!'
    write ( *, '(a,i8)' ) '  NNODE = ', nnode
    stop
  end if

  maxstack = 4 * nnode
  nstack = 0
  stack(1) = 0

  ncomp = 0
  nswap = 0
!
!  Start with the identity ordering.
!
  call i4vec_indicator ( nnode, order )
!
!  Compute the corresponding code.
!
  call cdmg_order_code ( adj, nnode, nnode, code, order )
!
!  Take this ordering and code as the best so far.
!
  bestorder(1:nnode) = order(1:nnode)
  bestcode(1:nnode,1:nnode) = code(1:nnode,1:nnode)
!
!  Now consider all possible orderings, and their codes.
!
  index = 0

  do

    call i4vec_backtrack ( nnode, order, index, nopen, nstack, stack, &
      maxstack, ncan )
!
!  If the backtracking routine has returned a complete candidate ordering, then
!  compute the resulting code, and see it it is better
!  then our current best.  Then go back for the next backtrack search.
!
    if ( index == 1 ) then

      call cdmg_order_code ( adj, nnode, nnode, code, order )

      call cdmg_code_compare ( bestcode, code, nnode, nnode, result )

      ncomp = ncomp + 1

      if ( result == - 1 ) then

        nswap = nswap + 1

        bestorder(1:nnode) = order(1:nnode)
        bestcode(1:nnode,1:nnode) = code(1:nnode,1:nnode)

      end if
!
!  If the backtracking routine has a partial reordering,
!  supply candidates for the next item in the ordering.
!
    else if ( index == 2 ) then

      call cdmg_code_cand ( adj, bestcode, code, nnode, ncomp, nopen, &
        order, stack, maxstack, nstack, ncan )

    else

      exit

    end if

  end do
!
!  Once we have examined all possibilites, we are done.
!
!  Set the output ordering to the best ordering, and the output
!  code to the corresponding best code.
!
  order(1:nnode) = bestorder(1:nnode)
  code(1:nnode,1:nnode) = bestcode(1:nnode,1:nnode)

  if ( debug ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CDMG_CODE_BACK:'
    write ( *, '(a,i8)' ) '  Comparisons: ', ncomp
    write ( *, '(a,i8)' ) '  Swaps:	  ', nswap
  end if

  return
end
subroutine cdmg_code_brute ( adj, nnode, code, order )

!*****************************************************************************80
!
!! CDMG_CODE_BRUTE computes a color dimultigraph code via brute force.
!
!  Discussion:
!
!    The code is the "largest" order code (in the lexicographic sense)
!    over all possible node orderings.  The brute force method considers
!    every node ordering, computes the corresponding order code, and
!    takes the largest one encountered.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,I) is the color of node I; otherwise, ADJ(I,J) is
!    the number of edges from node I to node J.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) CODE(NNODE,NNODE), the code.
!
!    Output, integer ( kind = 4 ) ORDER(NNODE), the ordering of the nodes.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) bestcode(nnode,nnode)
  integer ( kind = 4 ) bestorder(nnode)
  integer ( kind = 4 ) code(nnode,nnode)
  logical even
  logical more
  integer ( kind = 4 ) ncomp
  integer ( kind = 4 ) nswap
  integer ( kind = 4 ) order(nnode)
  integer ( kind = 4 ) result

  ncomp = 0
  nswap = 0
!
!  Start with the identity ordering.
!
  call i4vec_indicator ( nnode, order )
!
!  Compute the corresponding code.
!
  call cdmg_order_code ( adj, nnode, nnode, code, order )
!
!  Take this ordering and code as the best so far.
!
  bestorder(1:nnode) = order(1:nnode)
  bestcode(1:nnode,1:nnode) = code(1:nnode,1:nnode)
!
!  Now consider all possible orderings, and their codes.
!
  more = .false.

  do

    call perm_next ( nnode, order, more, even )

    call cdmg_order_code ( adj, nnode, nnode, code, order )

    call cdmg_code_compare ( bestcode, code, nnode, nnode, result )

    ncomp = ncomp + 1

    if ( result == - 1 ) then

      nswap = nswap + 1

      bestorder(1:nnode) = order(1:nnode)
      bestcode(1:nnode,1:nnode) = code(1:nnode,1:nnode)

    end if

    if ( .not. more ) then
      exit
    end if

  end do
!
!  Once we have examined all possibilites, we are done.
!
!  Set the output ordering to the best ordering, and the output
!  code to the corresponding best code.
!
  order(1:nnode) = bestorder(1:nnode)
  code(1:nnode,1:nnode) = bestcode(1:nnode,1:nnode)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CDMG_CODE_BRUTE:'
  write ( *, '(a,i8)' ) '  Comparisons: ', ncomp
  write ( *, '(a,i8)' ) '  Swaps:       ', nswap

  return
end
subroutine cdmg_code_cand ( adj, bestcode, code, nnode, ncomp, nopen, &
  order, stack, maxstack, nstack, ncan )

!*****************************************************************************80
!
!! CDMG_CODE_CAND finds candidates for a maximal color dimultigraph code ordering.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,I) is the color of node I; otherwise, ADJ(I,J) is
!    the number of edges from node I to node J.
!
!    Input, integer ( kind = 4 ) BESTCODE(NNODE,NNODE), the best code so far.
!
!    Workspace, integer CODE(NNODE,NNODE).
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input/output, integer ( kind = 4 ) NCOMP, the number of code comparisons.
!    This routine updates NCOMP by 1 each time the routine is called.
!
!    Input, integer ( kind = 4 ) NOPEN, identifies the first open position in ORDER.
!    1 <= NOPEN <= NNODE.
!
!    Input, integer ( kind = 4 ) ORDER(NNODE), contains in entries 1 through NOPEN-1
!    the elements of the current partial list.
!
!    Input/output, integer ( kind = 4 ) STACK(MAXSTACK), used to store the new candidates.
!
!    Input, integer ( kind = 4 ) MAXSTACK, the maximum size of the STACK array.
!    A value of NNODE should be sufficient.
!
!    Input/output, integer ( kind = 4 ) NSTACK, the current length of the stack.
!    On output, NSTACK has been increased by the number of
!    candidates found.
!
!    Input/output, integer ( kind = 4 ) NCAN(NNODE), the number of candidates for
!    each position.
!
  implicit none

  integer ( kind = 4 ) maxstack
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) bestcode(nnode,nnode)
  integer ( kind = 4 ) code(nnode,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ifree(nnode)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) max_adj
  integer ( kind = 4 ) maxcolor
  integer ( kind = 4 ) ncan(nnode)
  integer ( kind = 4 ) ncomp
  integer ( kind = 4 ) nfree
  integer ( kind = 4 ) ni
  integer ( kind = 4 ) nj
  integer ( kind = 4 ) nopen
  integer ( kind = 4 ) nstack
  integer ( kind = 4 ) order(nnode)
  integer ( kind = 4 ) result
  integer ( kind = 4 ) stack(maxstack)

  if ( nopen < 1 .or. nnode < nopen ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CDMG_CODE_CAND - Fatal error 1!'
    write ( *, '(a)' ) '  1 <= NOPEN <= NNODE should be true, but'
    write ( *, '(a,i8)' ) '  NOPEN = ', nopen
    write ( *, '(a,i8)' ) '  NNODE = ', nnode
    stop
  end if
!
!  Start with no candidates.
!
  ncan(nopen) = 0
!
!  If we have fixed at least one entry of the list,
!
!    Compute the partial code;
!
!    Compare the partial code with the corresponding
!    part of the the code for the best ordering so far;
!
!    If the current incomplete ordering is actually LESS than the
!    current best, then bail out now, with zero candidates.
!
  if ( 1 < nopen ) then

    call cdmg_order_code ( adj, nnode, nopen-1, code, order )

    call cdmg_code_compare ( bestcode, code, nnode, nopen-1, result )

    ncomp = ncomp + 1

    if ( result == + 1 ) then
      ncan(nopen) = 0
      return
    end if

  end if
!
!  Get a list of those nodes which have not been used yet.
!
  nfree = nnode + 1 - nopen
  call perm_free ( order, nopen-1, ifree, nfree )
!
!  Our preferred candidates will be:
!
  do i = 1, nopen-1

    ncan(nopen) = 0

    ni = order(i)
!
!  ...note the maximum adjacency FROM NI to any unordered node NJ...
!
    max_adj = 0
    do j = 1, nfree
      nj = ifree(j)
      max_adj = max ( max_adj, adj(ni,nj) )
    end do
!
!   ...and take as candidates all unordered nodes NJ with maximal
!   adjacency FROM NI.
!
!   (We could weed candidates further by only taking the maximal color.)
!
    if ( 0 < max_adj ) then

      do j = 1, nfree

        nj = ifree(j)

        if ( adj(ni,nj) == max_adj ) then

          ncan(nopen) = ncan(nopen) + 1
          nstack = nstack + 1

          if ( maxstack < nstack ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'CDMG_CODE_CAND - Fatal error 2!'
            write ( *, '(a)' ) '  MAXSTACK < NSTACK !'
            write ( *, '(a,i8)' ) '  NSTACK =   ', nstack
            write ( *, '(a,i8)' ) '  MAXSTACK = ', maxstack
            stop
          end if

          stack(nstack) = nj

        end if

      end do

      return

    end if
!
!  Else, note the maximum adjacency TO NI from any unordered node NJ...
!
    max_adj = 0
    do j = 1, nfree
      nj = ifree(j)
      max_adj = max ( max_adj, adj(nj,ni) )
    end do
!
!   ...and take as candidates all unordered nodes NJ with maximal
!   adjacency TO NI.
!
!   (We could weed candidates further by only taking the maximal color.)
!
    if ( 0 < max_adj ) then

      do j = 1, nfree

        nj = ifree(j)

        if ( adj(nj,ni) == max_adj ) then

          ncan(nopen) = ncan(nopen) + 1
          nstack = nstack + 1

          if ( maxstack < nstack ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' )'CDMG_CODE_CAND - Fatal error 2!'
            write ( *, '(a)' ) '  MAXSTACK < NSTACK!'
            write ( *, '(a,i8)' ) '  NSTACK =   ', nstack
            write ( *, '(a,i8)' ) '  MAXSTACK = ', maxstack
            stop
          end if

          stack(nstack) = nj

        end if

      end do

      return

    end if

  end do
!
!  NO unordered nodes are connected in any way to ordered nodes.
!  This can happen in two ways:
!
!  * NOPEN = 1; (the list of used nodes is empty)
!  * The graph is disconnected;
!
!  In either case, we must now consider ALL free nodes.
!
!  Compute the maximal color.
!
  maxcolor = 0

  do i = 1, nfree
    ni = ifree(i)
    maxcolor = max ( maxcolor, adj(ni,ni) )
  end do
!
!  Take as candidates every node of color MAXCOLOR.
!
!  We could thin the list down, by looking ahead, and only taking
!  candidates of MAXCOLOR who also happen to have at least one free
!  out neighbor, and so on.
!
  ncan(nopen) = 0

  do i = 1, nfree

    ni = ifree(i)

    if ( adj(ni,ni) == maxcolor ) then

      ncan(nopen) = ncan(nopen) + 1
      nstack = nstack + 1

      if ( maxstack < nstack ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CDMG_CODE_CAND - Fatal error 6!'
        write ( *, '(a)' ) '  MAXSTACK < NSTACK!'
        write ( *, '(a,i8)' ) '  NSTACK =   ', nstack
        write ( *, '(a,i8)' ) '  MAXSTACK = ', maxstack
        stop
      end if

      stack(nstack) = ni

    end if

  end do
!
!  This should never happen:
!
  if ( ncan(nopen) == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CDMG_CODE_CAND - Fatal error 7!'
    write ( *, '(a)' ) '  No candidates, but there gotta be some!'
    stop
  end if

  return
end
subroutine cdmg_code_compare ( code1, code2, nnode, npart, result )

!*****************************************************************************80
!
!! CDMG_CODE_COMPARE compares two (partial) color dimultigraph codes.
!
!  Discussion:
!
!    CODE1 = CODE2 if every digit of both codes is equal.
!
!    Otherwise, traverse the entries in a funny diagonal way, suggested
!    by this diagram for the first 16 entries:
!
!       1  2  5 10
!       3  4  7 12
!       6  8  9 14
!      11 13 15 16
!
!    As we do that, we examine the corresponding digits of the two codes.
!    For the first entry, (I,J), where the two codes differ, we say:
!
!      if ( CODE1(I,J) < CODE2(I,J) ) then we say
!        CODE1 < CODE2
!      else
!        CODE2 < CODE1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) CODE1(NNODE,NNODE), CODE2(NNODE,NNODE),
!    two codes to be compared.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) NPART, specifies the portion of the codes
!    to compare.  NPART should be between 1 and NNODE.
!
!    If NPART = NNODE, then the full codes are compared.
!
!    If NPART < NNODE, then only entries corresponding to I and J
!    both less than or equal to NPART will be compared.
!
!    Output, integer ( kind = 4 ) RESULT:
!    -1, CODE1 < CODE2;
!     0, CODE1 = CODE2;
!    +1, CODE2 < CODE1.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) code1(nnode,nnode)
  integer ( kind = 4 ) code2(nnode,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) npart
  integer ( kind = 4 ) result

  do j = 1, npart

    do i = 1, j - 1

      if ( code1(i,j) < code2(i,j) ) then

        result = - 1
        return

      else if ( code2(i,j) < code1(i,j) ) then

        result = + 1
        return

      else if ( code1(j,i) < code2(j,i) ) then

        result = - 1
        return

      else if ( code2(j,i) < code1(j,i) ) then

        result = + 1
        return

      end if

    end do

    if ( code1(j,j) < code2(j,j) ) then

      result = - 1
      return

    else if ( code2(j,j) < code1(j,j) ) then

      result = + 1
      return

    end if

  end do

  result = 0

  return
end
subroutine cdmg_code_print ( nnode, code, title )

!*****************************************************************************80
!
!! CDMG_CODE_PRINT prints out a color dimultigraph code.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) CODE(NNODE,NNODE), the code.
!
!    Input, character ( len = * ) TITLE, a title to be printed.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) code(nnode,nnode)
  integer ( kind = 4 ) i
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  do i = 1, nnode

    write ( *, '(2x,78i1)' ) code(i,1:nnode)

  end do

  return
end
subroutine cdmg_color_count ( adj, nnode, mcolor, ncolor )

!*****************************************************************************80
!
!! CDMG_COLOR_COUNT counts the number of colors in a color dimultigraph.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,I) is the color of node I; otherwise, ADJ(I,J) is
!    the number of edges from node I to node J.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) MCOLOR, the maximum color index.
!
!    Output, integer ( kind = 4 ) NCOLOR, the number of colors.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) colors(nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) mcolor
  integer ( kind = 4 ) ncolor

  mcolor = 0
  do i = 1, nnode
    mcolor = max ( mcolor, adj(i,i) )
  end do

  do i = 1, nnode
    colors(i) = adj(i,i)
  end do

  call i4vec_sort_heap_d ( nnode, colors )

  call i4vec_sorted_unique_count ( nnode, colors, ncolor )

  return
end
subroutine cdmg_color_sequence ( adj, nnode, seq )

!*****************************************************************************80
!
!! CDMG_COLOR_SEQUENCE computes the color sequence of a color dimultigraph.
!
!  Discussion:
!
!    The color sequence of a color dimultigraph is constructed by computing the
!    color of each node, and then ordering these values in decreasing order.
!
!    If two color dimultigraphs are isomorphic, they must have the same
!    color sequence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,I) is the color of node I; otherwise, ADJ(I,J) is
!    the number of edges from node I to node J.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) SEQ(NNODE), the color sequence.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) seq(nnode)

  do i = 1, nnode
    seq(i) = adj(i,i)
  end do

  call i4vec_sort_heap_d ( nnode, seq )

  return
end
subroutine cdmg_compare ( adj1, nnode1, adj2, nnode2, order1, &
  order2, result )

!*****************************************************************************80
!
!! CDMG_COMPARE determines if color dimultigraphs G1 and G2 are isomorphic.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ1(NNODE1,NNODE1), the adjacency information
!    for G1.
!
!    Input, integer ( kind = 4 ) NNODE1, the number of nodes in G1.
!
!    Input, integer ( kind = 4 ) ADJ2(NNODE2,NNODE2), the adjacency information
!    for G2.
!
!    Input, integer ( kind = 4 ) NNODE2, the number of nodes in G2.
!
!    Output, integer ( kind = 4 ) RESULT:
!     0 if the dimultigraphs are isomorphic,
!    -I if G1 < G2 for test #I,
!    +I if G2 < G1 for test #I.
!
!    Output, integer ( kind = 4 ) ORDER1(NNODE1), ORDER2(NNODE2).
!    If RESULT = 0, then ORDER1 and ORDER2 are reorderings of the nodes of
!    G1 and G2 which exhibit the isomorphism.
!
  implicit none

  integer ( kind = 4 ) nnode1
  integer ( kind = 4 ) nnode2

  integer ( kind = 4 ) adj_max_max_1
  integer ( kind = 4 ) adj_max_max_2
  integer ( kind = 4 ) adj1(nnode1,nnode1)
  integer ( kind = 4 ) adj2(nnode2,nnode2)
  integer ( kind = 4 ) code1(nnode1,nnode1)
  integer ( kind = 4 ) code2(nnode2,nnode2)
  integer ( kind = 4 ) in_seq1(nnode1)
  integer ( kind = 4 ) in_seq2(nnode2)
  integer ( kind = 4 ) mcolor1
  integer ( kind = 4 ) mcolor2
  integer ( kind = 4 ) ncolor1
  integer ( kind = 4 ) ncolor2
  integer ( kind = 4 ) nedge_u_1
  integer ( kind = 4 ) nedge_u_2
  integer ( kind = 4 ) nedge_w_1
  integer ( kind = 4 ) nedge_w_2
  integer ( kind = 4 ) order1(nnode1)
  integer ( kind = 4 ) order2(nnode2)
  integer ( kind = 4 ) out_seq1(nnode1)
  integer ( kind = 4 ) out_seq2(nnode2)
  integer ( kind = 4 ) result
  integer ( kind = 4 ) seq1(nnode1)
  integer ( kind = 4 ) seq2(nnode2)
!
!  Test 1: Count the nodes.
!
  if ( nnode1 < nnode2 ) then
    result = - 1
    return
  else if ( nnode2 < nnode1 ) then
    result = + 1
    return
  end if
!
!  Test 2: Compare the unweighted edges.
!
  call cdmg_edge_count ( adj1, nnode1, nedge_u_1, nedge_w_1 )

  call cdmg_edge_count ( adj2, nnode2, nedge_u_2, nedge_w_2 )

  if ( nedge_u_1 < nedge_u_2 ) then
    result = - 2
    return
  else if ( nedge_u_2 < nedge_u_1 ) then
    result = + 2
    return
  end if
!
!  Test 3: Compare the weighted edges.
!
  if ( nedge_w_1 < nedge_w_2 ) then
    result = - 3
    return
  else if ( nedge_w_2 < nedge_w_1 ) then
    result = + 3
    return
  end if
!
!  Test 4: Compare the number of colors.
!
  call cdmg_color_count ( adj1, nnode1, mcolor1, ncolor1 )

  call cdmg_color_count ( adj2, nnode2, mcolor2, ncolor2 )

  if ( ncolor1 < ncolor2 ) then
    result = - 4
    return
  else if ( ncolor2 < ncolor1 ) then
    result = + 4
    return
  end if
!
!  Test 5: Compare the maximum color.
!
  if ( mcolor1 < mcolor2 ) then
    result = - 5
    return
  else if ( mcolor2 < mcolor1 ) then
    result = + 5
    return
  end if
!
!  Test 6: Compare the color sequences.
!
  call cdmg_color_sequence ( adj1, nnode1, in_seq1 )

  call cdmg_color_sequence ( adj2, nnode2, in_seq2 )

  call i4vec_compare ( nnode1, in_seq1, in_seq2, result )

  if ( result < 0 ) then
    result = - 6
    return
  else if ( 0 < result ) then
    result = + 6
    return
  end if
!
!  Test 7: Compare the unweighted outdegree sequences.
!
  call cdmg_degree_seq_u ( adj1, nnode1, in_seq1, out_seq1 )

  call cdmg_degree_seq_u ( adj2, nnode2, in_seq2, out_seq2 )

  call i4vec_compare ( nnode1, out_seq1, out_seq2, result )

  if ( result < 0 ) then
    result = - 7
    return
  else if ( 0 < result ) then
    result = + 7
    return
  end if
!
!  Test 8: Compare the unweighted indegree sequences.
!
  call i4vec_compare ( nnode1, in_seq1, in_seq2, result )

  if ( result < 0 ) then
    result = - 8
    return
  else if ( 0 < result ) then
    result = + 8
    return
  end if
!
!  Test 9: Compare the adjacency max max.
!
  call cdmg_adj_max_max ( adj1, nnode1, adj_max_max_1 )

  call cdmg_adj_max_max ( adj2, nnode2, adj_max_max_2 )

  if ( adj_max_max_1 < adj_max_max_2 ) then
    result = - 9
    return
  else if ( adj_max_max_1 < adj_max_max_1 ) then
    result = + 9
    return
  end if
!
!  Test 10: Compare the adjacency max sequences.
!
  call cdmg_adj_max_seq ( adj1, nnode1, seq1 )

  call cdmg_adj_max_seq ( adj2, nnode2, seq2 )

  call i4vec_compare ( nnode1, seq1, seq2, result )

  if ( result < 0 ) then
    result = - 10
    return
  else if ( 0 < result ) then
    result = + 10
    return
  end if
!
!  Test 11: Compare the weighted outdegree sequences.
!
  call cdmg_degree_seq_w ( adj1, nnode1, in_seq1, out_seq1 )

  call cdmg_degree_seq_w ( adj2, nnode2, in_seq2, out_seq2 )

  call i4vec_compare ( nnode1, out_seq1, out_seq2, result )

  if ( result < 0 ) then
    result = - 11
    return
  else if ( 0 < result ) then
    result = + 11
    return
  end if
!
!  Test 12: Compare the weighted indegree sequences.
!
  call i4vec_compare ( nnode1, in_seq1, in_seq2, result )

  if ( result < 0 ) then
    result = - 12
    return
  else if ( 0 < result ) then
    result = + 12
    return
  end if
!
!  Test 13: Compare the weighted adjacency sequences.
!
  call cdmg_adj_seq_w ( adj1, nnode1, code1 )

  call cdmg_adj_seq_w ( adj2, nnode2, code2 )

  call i4mat_row_compare ( nnode1, nnode1, code1, code2, result )

  if ( result < 0 ) then
    result = - 13
    return
  else if ( 0 < result ) then
    result = + 13
    return
  end if
!
!  Test 14: Compare the codes.
!
  call cdmg_code_back ( adj1, nnode1, code1, order1 )

  call cdmg_code_back ( adj2, nnode2, code2, order2 )

  call cdmg_code_compare ( code1, code2, nnode1, nnode1, result )

  if ( result < 0 ) then
    result = - 14
    return
  else if ( 0 < result ) then
    result = + 14
    return
  end if

  result = 0

  return
end
subroutine cdmg_degree_seq_u ( adj, nnode, in_seq, out_seq )

!*****************************************************************************80
!
!! CDMG_DEGREE_SEQ_U: unweighted directed degree sequence of color dimultigraph.
!
!  Discussion:
!
!    The unweighted directed degree sequence is the sequence of indegrees
!    and the sequence of outdegrees, ignoring edge multiplicity, arranged
!    to correspond to nodes of successively decreasing total degree.  For
!    nodes of equal degree, those of higher outdegree take precedence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,I) is the color of node I; otherwise, ADJ(I,J) is
!    the number of edges from node I to node J.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) IN_SEQ(NNODE), OUT_SEQ(NNODE),
!    the unweighted directed degree sequences.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) in_seq(nnode)
  integer ( kind = 4 ) out_seq(nnode)

  call cdmg_degree_u ( adj, nnode, in_seq, out_seq )

  call i4vec2_sort_d ( nnode, out_seq, in_seq )

  return
end
subroutine cdmg_degree_seq_w ( adj, nnode, in_seq, out_seq )

!*****************************************************************************80
!
!! CDMG_DEGREE_SEQ_W: weighted directed degree sequence of a color dimultigraph.
!
!  Discussion:
!
!    The weighted directed degree sequence is the sequence of indegrees
!    and the sequence of outdegrees, with edge multiplicity, arranged
!    to correspond to nodes of successively decreasing total degree.  For
!    nodes of equal degree, those of higher outdegree take precedence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,I) is the color of node I; otherwise, ADJ(I,J) is
!    the number of edges from node I to node J.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) IN_SEQ(NNODE), OUT_SEQ(NNODE),
!    the weighted directed degree sequences.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) in_seq(nnode)
  integer ( kind = 4 ) out_seq(nnode)

  call cdmg_degree_w ( adj, nnode, in_seq, out_seq )

  call i4vec2_sort_d ( nnode, out_seq, in_seq )

  return
end
subroutine cdmg_degree_u ( adj, nnode, indegree, outdegree )

!*****************************************************************************80
!
!! CDMG_DEGREE_U computes the unweighted degrees of a color dimultigraph.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,I) is the color of node I; otherwise, ADJ(I,J) is
!    the number of edges from node I to node J.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) INDEGREE(NNODE), OUTDEGREE(NNODE),
!    the unweighted indegree and outdegree of the nodes.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indegree(nnode)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) outdegree(nnode)

  indegree(1:nnode) = 0
  outdegree(1:nnode) = 0

  do i = 1, nnode
    do j = 1, nnode
      if ( i /= j ) then
        if ( adj(i,j) /= 0 ) then
          outdegree(i) = outdegree(i) + 1
          indegree(j) = indegree(j) + 1
        end if
      end if
    end do
  end do

  return
end
subroutine cdmg_degree_w ( adj, nnode, indegree, outdegree )

!*****************************************************************************80
!
!! CDMG_DEGREE_W computes the weighted degrees of a color dimultigraph.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,I) is the color of node I; otherwise, ADJ(I,J) is
!    the number of edges from node I to node J.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) INDEGREE(NNODE), OUTDEGREE(NNODE),
!    the weighted indegree and outdegree of the nodes.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indegree(nnode)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) outdegree(nnode)

  indegree(1:nnode) = 0
  outdegree(1:nnode) = 0

  do i = 1, nnode
    do j = 1, nnode
      if ( i /= j ) then
        if ( adj(i,j) /= 0 ) then
          outdegree(i) = outdegree(i) + adj(i,j)
          indegree(j) = indegree(j) + adj(i,j)
        end if
      end if
    end do
  end do

  return
end
subroutine cdmg_edge_count ( adj, nnode, nedge_u, nedge_w )

!*****************************************************************************80
!
!! CDMG_EDGE_COUNT counts the number of edges in a color dimultigraph.
!
!  Discussion:
!
!    The number of "unweighted" edges is the number of edges in the
!    underlying digraph, or the number of edges that would be counted
!    if each set of multiple edges was replaced by a single edge.
!
!    The number of "weighted" edges counts separately each edge of a
!    multiple edge.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,I) is the color of node I; otherwise, ADJ(I,J) is
!    the number of edges from node I to node J.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) NEDGE_U, the number of unweighted edges.
!
!    Output, integer ( kind = 4 ) NEDGE_W, the number of weighted edges.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nedge_u
  integer ( kind = 4 ) nedge_w

  nedge_u = 0
  nedge_w = 0

  do i = 1, nnode
    do j = 1, nnode

      if ( i /= j ) then
        nedge_w = nedge_w + adj(i,j)
        if ( 0 < adj(i,j) ) then
          nedge_u = nedge_u + 1
        end if
      end if

    end do
  end do

  return
end
subroutine cdmg_example_octo ( example, adj, nnode, seed )

!*****************************************************************************80
!
!! CDMG_EXAMPLE_OCTO sets up an 8 node example color dimultigraph.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) EXAMPLE, chooses the example, and should be between
!    1 and 12.
!
!    Output, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information for the graph.
!    ADJ(I,I) is the color of node I.
!    ADJ(I,J) is 1 if nodes I and J are adjacent and 0 otherwise.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes, which should be 8.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(8,8)
  integer ( kind = 4 ), parameter :: BLUE = 1
  integer ( kind = 4 ) example
  integer ( kind = 4 ), parameter :: GREEN = 2
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) nsave
  integer ( kind = 4 ), parameter :: RED = 3
  integer ( kind = 4 ) seed
  integer ( kind = 4 ), parameter :: YELLOW = 5
  integer ( kind = 4 ), parameter :: ZIRCON = 4

  if ( example <= 0 ) then
    nsave = i4_uniform ( 1, 12, seed )
  else
    nsave = mod ( example - 1, 12 ) + 1
  end if

  nnode = 8

  adj(1:nnode,1:nnode) = 0
!
!  #1.
!
  if ( nsave == 1 ) then

    adj(1,1) = BLUE
    adj(1,2) = 2
    adj(1,6) = 2
    adj(2,2) = BLUE
    adj(2,3) = 3
    adj(2,5) = 1
    adj(3,3) = BLUE
    adj(3,4) = 1
    adj(3,8) = 4
    adj(4,4) = GREEN
    adj(4,5) = 1
    adj(4,7) = 2
    adj(5,5) = GREEN
    adj(5,6) = 1
    adj(6,6) = RED
    adj(6,7) = 1
    adj(7,7) = RED
    adj(7,8) = 2
    adj(8,1) = 1
    adj(8,8) = YELLOW
!
!  #2, same NNODE, different number of unweighted edges.
!
  else if ( nsave == 2 ) then

    adj(1,1) = BLUE
    adj(1,2) = 2
    adj(1,6) = 2
    adj(2,2) = BLUE
    adj(2,3) = 3
    adj(3,3) = BLUE
    adj(3,4) = 1
    adj(3,8) = 4
    adj(4,4) = GREEN
    adj(4,5) = 1
    adj(4,7) = 2
    adj(5,5) = GREEN
    adj(5,6) = 1
    adj(6,6) = RED
    adj(6,7) = 1
    adj(7,7) = RED
    adj(7,8) = 2
    adj(8,1) = 1
    adj(8,8) = YELLOW
!
!  #3, same NNODE, unweighted edges, different weighted edges.
!
  else if ( nsave == 3 ) then

    adj(1,1) = BLUE
    adj(1,2) = 2
    adj(1,6) = 2
    adj(2,2) = BLUE
    adj(2,3) = 3
    adj(2,5) = 1
    adj(3,3) = BLUE
    adj(3,4) = 1
    adj(3,8) = 4
    adj(4,4) = GREEN
    adj(4,5) = 1
    adj(4,7) = 2
    adj(5,5) = GREEN
    adj(5,6) = 1
    adj(6,6) = RED
    adj(6,7) = 1
    adj(7,7) = RED
    adj(7,8) = 1
    adj(8,1) = 1
    adj(8,8) = YELLOW
!
!  #4, different number of colors
!
  else if ( nsave == 4 ) then

    adj(1,1) = BLUE
    adj(1,2) = 2
    adj(1,6) = 2
    adj(2,2) = BLUE
    adj(2,3) = 3
    adj(2,5) = 1
    adj(3,3) = BLUE
    adj(3,4) = 1
    adj(3,8) = 4
    adj(4,4) = GREEN
    adj(4,5) = 1
    adj(4,7) = 2
    adj(5,5) = GREEN
    adj(5,6) = 1
    adj(6,6) = GREEN
    adj(6,7) = 1
    adj(7,7) = BLUE
    adj(7,8) = 2
    adj(8,1) = 1
    adj(8,8) = YELLOW
!
!  #5, different maximum color index.
!
  else if ( nsave == 5 ) then

    adj(1,1) = BLUE
    adj(1,2) = 2
    adj(1,6) = 2
    adj(2,2) = BLUE
    adj(2,3) = 3
    adj(2,5) = 1
    adj(3,3) = BLUE
    adj(3,4) = 1
    adj(3,8) = 4
    adj(4,4) = GREEN
    adj(4,5) = 1
    adj(4,7) = 2
    adj(5,5) = GREEN
    adj(5,6) = 1
    adj(6,6) = RED
    adj(6,7) = 1
    adj(7,7) = RED
    adj(7,8) = 2
    adj(8,1) = 1
    adj(8,8) = ZIRCON
!
!  #6, different color sequence.
!
  else if ( nsave == 6 ) then

    adj(1,1) = BLUE
    adj(1,2) = 2
    adj(1,6) = 2
    adj(2,2) = BLUE
    adj(2,3) = 3
    adj(2,5) = 1
    adj(3,3) = GREEN
    adj(3,4) = 1
    adj(3,8) = 4
    adj(4,4) = GREEN
    adj(4,5) = 1
    adj(4,7) = 2
    adj(5,5) = GREEN
    adj(5,6) = 1
    adj(6,6) = RED
    adj(6,7) = 1
    adj(7,7) = RED
    adj(7,8) = 2
    adj(8,1) = 1
    adj(8,8) = YELLOW
!
!  #7, unweighted outdegree sequence.
!
  else if ( nsave == 7 ) then

    adj(1,1) = BLUE
    adj(1,2) = 2
    adj(2,2) = BLUE
    adj(2,3) = 3
    adj(2,5) = 1
    adj(2,6) = 2
    adj(3,3) = BLUE
    adj(3,4) = 1
    adj(3,8) = 4
    adj(4,4) = GREEN
    adj(4,5) = 1
    adj(4,7) = 2
    adj(5,5) = GREEN
    adj(5,6) = 1
    adj(6,6) = RED
    adj(6,7) = 1
    adj(7,7) = RED
    adj(7,8) = 2
    adj(8,1) = 1
    adj(8,8) = YELLOW
!
!  #8, unweighted indegree sequence.
!
  else if ( nsave == 8 ) then

    adj(1,1) = BLUE
    adj(1,2) = 2
    adj(1,6) = 2
    adj(2,2) = BLUE
    adj(2,3) = 3
    adj(2,5) = 1
    adj(3,3) = BLUE
    adj(3,4) = 1
    adj(3,8) = 4
    adj(4,4) = GREEN
    adj(4,5) = 1
    adj(4,7) = 2
    adj(5,5) = GREEN
    adj(5,7) = 1
    adj(6,6) = RED
    adj(6,7) = 1
    adj(7,7) = RED
    adj(7,8) = 2
    adj(8,1) = 1
    adj(8,8) = YELLOW
!
!  #9, adjacency max max
!
  else if ( nsave == 9 ) then

    adj(1,1) = BLUE
    adj(1,2) = 2
    adj(1,6) = 2
    adj(2,2) = BLUE
    adj(2,3) = 3
    adj(2,5) = 1
    adj(3,3) = BLUE
    adj(3,4) = 1
    adj(3,8) = 3
    adj(4,4) = GREEN
    adj(4,5) = 1
    adj(4,7) = 3
    adj(5,5) = GREEN
    adj(5,6) = 1
    adj(6,6) = RED
    adj(6,7) = 1
    adj(7,7) = RED
    adj(7,8) = 2
    adj(8,1) = 1
    adj(8,8) = YELLOW
!
!  #10, adjacency max sequence.
!
  else if ( nsave == 10 ) then

    adj(1,1) = BLUE
    adj(1,2) = 2
    adj(1,6) = 2
    adj(2,2) = BLUE
    adj(2,3) = 2
    adj(2,5) = 1
    adj(3,3) = BLUE
    adj(3,4) = 1
    adj(3,8) = 4
    adj(4,4) = GREEN
    adj(4,5) = 1
    adj(4,7) = 2
    adj(5,5) = GREEN
    adj(5,6) = 2
    adj(6,6) = RED
    adj(6,7) = 1
    adj(7,7) = RED
    adj(7,8) = 2
    adj(8,1) = 1
    adj(8,8) = YELLOW
!
!  #11, weighted outdegree sequence
!
  else if ( nsave == 11 ) then

    adj(1,1) = BLUE
    adj(1,2) = 1
    adj(1,6) = 2
    adj(2,2) = BLUE
    adj(2,3) = 3
    adj(2,5) = 2
    adj(3,3) = BLUE
    adj(3,4) = 1
    adj(3,8) = 4
    adj(4,4) = GREEN
    adj(4,5) = 1
    adj(4,7) = 2
    adj(5,5) = GREEN
    adj(5,6) = 1
    adj(6,6) = RED
    adj(6,7) = 1
    adj(7,7) = RED
    adj(7,8) = 2
    adj(8,1) = 1
    adj(8,8) = YELLOW
!
!  #12, weighted indegree sequence.
!
  else if ( nsave == 12 ) then

    adj(1,1) = BLUE
    adj(1,2) = 2
    adj(1,6) = 1
    adj(2,2) = BLUE
    adj(2,3) = 3
    adj(2,5) = 2
    adj(3,3) = BLUE
    adj(3,4) = 1
    adj(3,8) = 4
    adj(4,4) = GREEN
    adj(4,5) = 1
    adj(4,7) = 2
    adj(5,5) = GREEN
    adj(5,6) = 1
    adj(6,6) = RED
    adj(6,7) = 1
    adj(7,7) = RED
    adj(7,8) = 2
    adj(8,1) = 1
    adj(8,8) = YELLOW
!
!  #13: weighted adjacency sequence  NOT SET UP YET
!
  else if ( nsave == 13 ) then

    adj(1,1) = BLUE
    adj(1,2) = 2
    adj(1,6) = 2
    adj(2,2) = BLUE
    adj(2,3) = 3
    adj(2,5) = 1
    adj(3,3) = BLUE
    adj(3,4) = 1
    adj(3,8) = 4
    adj(4,4) = GREEN
    adj(4,5) = 1
    adj(4,7) = 2
    adj(5,5) = GREEN
    adj(5,6) = 1
    adj(6,6) = RED
    adj(6,7) = 1
    adj(7,7) = RED
    adj(7,8) = 2
    adj(8,1) = 1
    adj(8,8) = YELLOW
!
!  #14: code  NOT SET UP YET
!
  else if ( nsave == 14 ) then

    adj(1,1) = BLUE
    adj(1,2) = 2
    adj(1,6) = 2
    adj(2,2) = BLUE
    adj(2,3) = 3
    adj(2,5) = 1
    adj(3,3) = BLUE
    adj(3,4) = 1
    adj(3,8) = 4
    adj(4,4) = GREEN
    adj(4,5) = 1
    adj(4,7) = 2
    adj(5,5) = GREEN
    adj(5,6) = 1
    adj(6,6) = RED
    adj(6,7) = 1
    adj(7,7) = RED
    adj(7,8) = 2
    adj(8,1) = 1
    adj(8,8) = YELLOW

  end if
!
!  Now permute the graph.
!
  call i4mat_perm_random ( nnode, adj, seed )

  return
end
subroutine cdmg_order_code ( adj, nnode, npart, code, order )

!*****************************************************************************80
!
!! CDMG_ORDER_CODE returns the color dimultigraph code for a specific node ordering.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,I) is the color of node I; otherwise, ADJ(I,J) is
!    the number of edges from node I to node J.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) NPART, the number of nodes to consider.
!    NPART should be between 1 and NNODE.
!
!    If NPART is NNODE, then the usual code is returned.
!
!    If NPART is less than NNODE, then the code is computed as
!    though only NPART nodes existed, namely, those specified in the
!    first NPART entries of order.  This option is provided so that
!    the routine can compute the portion of a code determined
!    by an incomplete ordering of the nodes.
!
!    Output, integer ( kind = 4 ) CODE(NNODE,NNODE), the code for this ordering.
!
!    Input, integer ( kind = 4 ) ORDER(NNODE), the ordering of the nodes.  ORDER(1)
!    is the first node, and so on.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) code(nnode,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ni
  integer ( kind = 4 ) nj
  integer ( kind = 4 ) npart
  integer ( kind = 4 ) order(nnode)

  do i = 1, nnode

    if ( i <= npart ) then

      ni = order(i)

      if ( order(i) < 1 .or. nnode < order(i) ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CDMG_ORDER_CODE - Fatal error!'
        write ( *, '(a)' ) '  ORDER is not a proper permutation.'
        stop
      end if

    else
      ni = 0
    end if

    do j = 1, nnode

      if ( j <= npart ) then

        nj = order(j)

        if ( order(j) < 1 .or. nnode < order(j) ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'CDMG_ORDER_CODE - Fatal error!'
          write ( *, '(a)' ) '  ORDER is not a proper permutation.'
          stop
        end if

      else
        nj = 0
      end if

      if ( ni == 0 .or. nj == 0 ) then

        code(i,j) = 0

      else

        code(i,j) = adj(ni,nj)

      end if

    end do
  end do

  return
end
subroutine cdmg_print ( adj, nnode, title )

!*****************************************************************************80
!
!! CDMG_PRINT prints out an adjacency matrix for a color dimultigraph.
!
!  Discussion:
!
!    Color values between 1 and 10 will be printed as
!      'R', 'G', 'B', 'C', 'M', 'Y', 'K', 'W', 'P', 'O'
!
!    Adjacency values between 0 and 9 will be printed as is.
!    Other values will be printed as '*'.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,I) is the color of node I; otherwise, ADJ(I,J) is
!    the number of edges from node I to node J.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, character ( len = * ) TITLE, a title to be printed.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  character, dimension ( 10 ) :: color = &
    (/ 'R', 'G', 'B', 'C', 'M', 'Y', 'K', 'W', 'P', '0' /)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  character ( len = 80 ) string
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  do i = 1, nnode

    jhi = min ( nnode, 80 )

    do j = 1, jhi

      if ( i == j ) then

        if ( 1 <= adj(i,j) .and. adj(i,j) <= 10 ) then
          string(j:j) = color ( adj(i,j) )
        else
          string(j:j) = '*'
        end if

      else

        if ( 0 <= adj(i,j) .and. adj(i,j) <= 9 ) then
          string(j:j) = char ( 48 + adj(i,j) )
        else
          string(j:j) = '*'
        end if

      end if

    end do

    write ( *, '(2x,a)' ) string(1:jhi)

  end do

  return
end
subroutine cdmg_random ( adj, nnode, ncolor, nedge, seed )

!*****************************************************************************80
!
!! CDMG_RANDOM generates a random color dimultigraph.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency matrix.
!    ADJ(I,J) is the number of edges from node I to node J.
!    ADJ(I,I) will always be 0.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) NCOLOR, the number of colors.
!    Each node is assumed to have an associated color, between 1 and NCOLOR,
!    which will be chosen at random.
!
!    Input, integer ( kind = 4 ) NEDGE, the number of edges.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
  implicit none

  integer ( kind = 4 ) ncolor
  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) color_i
  integer ( kind = 4 ) edge_i
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) node_i
  integer ( kind = 4 ) node_j
  integer ( kind = 4 ) perm(ncolor)
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) subset(ncolor)
!
!  Initialize the adjacency matrix.
!
  adj(1:nnode,1:nnode) = 0
!
!  Choose the colors.
!
  call ksub_random ( nnode, ncolor, subset, seed )

  call perm_random ( ncolor, perm, seed )

  do color_i = 1, ncolor
    node_i = subset(perm(color_i))
    adj(node_i,node_i) = color_i
  end do

  do node_i = 1, nnode
    if ( adj(node_i,node_i) == 0 ) then
      adj(node_i,node_i) = i4_uniform ( 1, ncolor, seed )
    end if
  end do
!
!  Essentially, flip a coin NEDGE times to decide where each edge goes.
!
  do edge_i = 1, nedge

    node_i = i4_uniform ( 1, nnode, seed )
    node_j = i4_uniform ( 1, nnode-1, seed )

    if ( node_i <= node_j ) then
      node_j = node_j + 1
    end if

    adj(node_i,node_j) = adj(node_i,node_j) + 1

  end do

  return
end
subroutine cg_code_back ( adj, nnode, code, order )

!*****************************************************************************80
!
!! CG_CODE_BACK computes a color graph code via backtracking.
!
!  Discussion:
!
!    The code is the "largest" order code over all possible node orderings.
!    The lexicographic ordering is used in comparing codes.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,I) is the color of node I; otherwise, ADJ(I,J) is positive
!    if there is an edge between node I and node J.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) CODE(NNODE,NNODE), the code.
!
!    Output, integer ( kind = 4 ) ORDER(NNODE), the ordering of the nodes.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) bestcode(nnode,nnode)
  integer ( kind = 4 ) bestorder(nnode)
  integer ( kind = 4 ) code(nnode,nnode)
  logical, parameter :: debug = .false.
  integer ( kind = 4 ) index
  integer ( kind = 4 ) maxstack
  integer ( kind = 4 ) ncan(nnode)
  integer ( kind = 4 ) ncomp
  integer ( kind = 4 ) nopen
  integer ( kind = 4 ) nstack
  integer ( kind = 4 ) nswap
  integer ( kind = 4 ) order(nnode)
  integer ( kind = 4 ) result
  integer ( kind = 4 ) stack((nnode*(nnode+1))/2)

  if ( nnode <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CG_CODE_BACK - Fatal error!'
    write ( *, '(a,i8)' ) '  NNODE = ', nnode
    stop
  end if

  maxstack = ( nnode * ( nnode + 1 ) ) / 2
  nstack = 0
  stack(1) = 0

  ncomp = 0
  nswap = 0
!
!  Start with the identity ordering.
!
  call i4vec_indicator ( nnode, order )
!
!  Compute the corresponding code.
!
  call cg_order_code ( adj, nnode, nnode, code, order )
!
!  Take this ordering and code as the best so far.
!
  bestorder(1:nnode) = order(1:nnode)
  bestcode(1:nnode,1:nnode) = code(1:nnode,1:nnode)
!
!  Now consider all possible orderings, and their codes.
!
  index = 0

  do

    call i4vec_backtrack ( nnode, order, index, nopen, nstack, stack, &
      maxstack, ncan )
!
!  If the backtrack routine has returned a complete candidate ordering, then
!  compute the resulting code, and see it it is better
!  then our current best.  Then go back for the next backtrack search.
!
    if ( index == 1 ) then

      call cg_order_code ( adj, nnode, nnode, code, order )

      call cg_code_compare ( bestcode, code, nnode, nnode, result )

      ncomp = ncomp + 1

      if ( result == - 1 ) then

        nswap = nswap + 1

        bestorder(1:nnode) = order(1:nnode)
        bestcode(1:nnode,1:nnode) = code(1:nnode,1:nnode)

      end if
!
!  If the backtrack routine has a partial reordering,
!  supply candidates for the next item in the ordering.
!
    else if ( index == 2 ) then

      call cg_code_cand ( adj, bestcode, code, nnode, ncomp, nopen, &
        order, stack, maxstack, nstack, ncan )

    else

      exit

    end if

  end do
!
!  Once we have examined all possibilites, we are done.
!
!  Set the output ordering to the best ordering, and the output
!  code to the corresponding best code.
!
  order(1:nnode) = bestorder(1:nnode)
  code(1:nnode,1:nnode) = bestcode(1:nnode,1:nnode)

  if ( debug ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CG_CODE_BACK:'
    write ( *, '(a,i8)' ) '  Comparisons: ', ncomp
    write ( *, '(a,i8)' ) '  Swaps:       ', nswap
  end if

  return
end
subroutine cg_code_brute ( adj, nnode, code, order )

!*****************************************************************************80
!
!! CG_CODE_BRUTE computes the color graph code via brute force.
!
!  Discussion:
!
!    The code is the "largest" order code over all node orderings.
!    The lexicographic ordering is used in comparing codes.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,I) is the color of node I; otherwise, ADJ(I,J) is positive
!    if there is an edge between node I and node J.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) CODE(NNODE,NNODE), the code.
!
!    Output, integer ( kind = 4 ) ORDER(NNODE), the ordering of the nodes.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) bestcode(nnode,nnode)
  integer ( kind = 4 ) bestorder(nnode)
  integer ( kind = 4 ) code(nnode,nnode)
  logical even
  logical more
  integer ( kind = 4 ) ncomp
  integer ( kind = 4 ) nswap
  integer ( kind = 4 ) order(nnode)
  integer ( kind = 4 ) result

  ncomp = 0
  nswap = 0
!
!  Start with the identity ordering.
!
  call i4vec_indicator ( nnode, order )
!
!  Compute the corresponding code.
!
  call cg_order_code ( adj, nnode, nnode, code, order )
!
!  Take this ordering and code as the best so far.
!
  bestorder(1:nnode) = order(1:nnode)
  bestcode(1:nnode,1:nnode) = code(1:nnode,1:nnode)
!
!  Now consider all possible orderings, and their codes.
!
  more = .false.

  do

    call perm_next ( nnode, order, more, even )

    call cg_order_code ( adj, nnode, nnode, code, order )

    call cg_code_compare ( bestcode, code, nnode, nnode, result )

    ncomp = ncomp + 1

    if ( result == - 1 ) then

      nswap = nswap + 1

      bestorder(1:nnode) = order(1:nnode)
      bestcode(1:nnode,1:nnode) = code(1:nnode,1:nnode)

    end if

    if ( .not. more ) then
      exit
    end if

  end do
!
!  Once we have examined all possibilites, we are done.
!
!  Set the output ordering to the best ordering, and the output
!  code to the corresponding best code.
!
  order(1:nnode) = bestorder(1:nnode)
  code(1:nnode,1:nnode) = bestcode(1:nnode,1:nnode)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CG_CODE_BRUTE:'
  write ( *, '(a,i8)' ) '  Comparisons: ', ncomp
  write ( *, '(a,i8)' ) '  Swaps:       ', nswap

  return
end
subroutine cg_code_cand ( adj, bestcode, code, nnode, ncomp, nopen, &
  order, stack, maxstack, nstack, ncan )

!*****************************************************************************80
!
!! CG_CODE_CAND finds candidates for a maximal color graph code ordering.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,I) is the color of node I; otherwise, ADJ(I,J) is positive
!    if there is an edge between node I and node J.
!
!    Input, integer ( kind = 4 ) BESTCODE(NNODE,NNODE), the best code so far.
!
!    Workspace, integer CODE(NNODE,NNODE).
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input/output, integer ( kind = 4 ) NCOMP, the number of code comparisons.
!    This routine updates NCOMP by 1 each time it is called.
!
!    Input, integer ( kind = 4 ) NOPEN, identifies the first open position in ORDER.
!    1 <= NOPEN <= NNODE.
!
!    Input, integer ( kind = 4 ) ORDER(NNODE), contains in entries 1 through NOPEN-1
!    the elements of the current partial list.
!
!    Input/output, integer ( kind = 4 ) STACK(MAXSTACK), used to store the new candidates.
!
!    Input, integer ( kind = 4 ) MAXSTACK, the maximum size of the STACK array.
!    A value of NNODE should be sufficient.
!
!    Input/output, integer ( kind = 4 ) NSTACK, the current length of the stack.
!    On output, NSTACK has been increased by the number of
!    candidates found.
!
!    Input/output, integer ( kind = 4 ) NCAN(NNODE), the number of candidates for
!    each position.
!
  implicit none

  integer ( kind = 4 ) maxstack
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) bestcode(nnode,nnode)
  integer ( kind = 4 ) code(nnode,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ifree(nnode)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) maxcolor
  integer ( kind = 4 ) ncan(nnode)
  integer ( kind = 4 ) ncomp
  integer ( kind = 4 ) nfree
  integer ( kind = 4 ) ni
  integer ( kind = 4 ) nj
  integer ( kind = 4 ) nopen
  integer ( kind = 4 ) nstack
  integer ( kind = 4 ) order(nnode)
  integer ( kind = 4 ) result
  integer ( kind = 4 ) stack(maxstack)

  if ( nopen < 1 .or. nnode < nopen ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CG_CODE_CAND - Fatal error 1!'
    write ( *, '(a)' ) '  1 <= NOPEN <= NNODE should be true, but'
    write ( *, '(a,i8)' ) '  NOPEN = ', nopen
    write ( *, '(a,i8)' ) '  NNODE = ', nnode
    stop
  end if
!
!  Start with no candidates.
!
  ncan(nopen) = 0
!
!  If we have fixed at least one entry of the list,
!
!    Compute the partial code;
!
!    Compare the partial code with the corresponding
!    part of the the code for the best ordering so far;
!
!    If the current incomplete ordering is actually LESS than the
!    current best, then bail out now, with zero candidates.
!
  if ( 1 < nopen ) then

    call cg_order_code ( adj, nnode, nopen-1, code, order )

    call cg_code_compare ( bestcode, code, nnode, nopen-1, result )

    ncomp = ncomp + 1

    if ( result == + 1 ) then
      ncan(nopen) = 0
      return
    end if

  end if
!
!  Get a list of those nodes which have not been used yet.
!
  nfree = nnode + 1 - nopen
  call perm_free ( order, nopen-1, ifree, nfree )
!
!  Our preferred candidates will be
!    * unused neighbors of the LOWEST ordered node possible.
!
  ncan(nopen) = 0

  do i = 1, nopen-1

    ni = order(i)

    do j = 1, nfree

      nj = ifree(j)

      if ( adj(ni,nj) /= 0 .or. adj(nj,ni) /= 0 ) then

        ncan(nopen) = ncan(nopen) + 1
        nstack = nstack + 1

        if ( maxstack < nstack ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'CG_CODE_CAND - Fatal error 4!'
          write ( *, '(a)' ) '  MAXSTACK < NSTACK!'
          write ( *, '(a,i8)' ) '  NSTACK =   ', nstack
          write ( *, '(a,i8)' ) '  MAXSTACK = ', maxstack
          stop
        end if

        stack(nstack) = nj

      end if

    end do
!
!  If in the middle of this loop, we found unused neighbors of the
!  lowest ordered node possible, then these are the only candidates
!  worth considering.
!
    if ( 0 < ncan(nopen) ) then
      return
    end if

  end do
!
!  If we get here, then NO unused nodes are connected in any way to
!  used nodes.  This can happen in two ways:
!
!  * NOPEN = 1; (the list of used nodes is empty)
!  * The graph is disconnected;
!
!  In either case, we must now consider ALL free nodes.
!
!  Compute the maximal color.
!
  maxcolor = 0

  do i = 1, nfree
    ni = ifree(i)
    maxcolor = max ( maxcolor, adj(ni,ni) )
  end do
!
!  Take as candidates every node of color MAXCOLOR.
!
!  We could thin the list down, by looking ahead, and only taking
!  candidates of MAXCOLOR who also happen to have at least one free
!  out neighbor, and so on.
!
  ncan(nopen) = 0

  do i = 1, nfree

    ni = ifree(i)

    if ( adj(ni,ni) == maxcolor ) then

      ncan(nopen) = ncan(nopen) + 1
      nstack = nstack + 1

      if ( maxstack < nstack ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CG_CODE_CAND - Fatal error 6!'
        write ( *, '(a)' ) '  MAXSTACK < NSTACK!'
        write ( *, '(a,i8)' ) '  NSTACK =   ', nstack
        write ( *, '(a,i8)' ) '  MAXSTACK = ', maxstack
        stop
      end if

      stack(nstack) = ni

    end if

  end do
!
!  This should never happen:
!
  if ( ncan(nopen) == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CG_CODE_CAND - Fatal error 7!'
    write ( *, '(a)' ) '  No candidates, but there gotta be some!'
    stop
  end if

  return
end
subroutine cg_code_compare ( code1, code2, nnode, npart, result )

!*****************************************************************************80
!
!! CG_CODE_COMPARE compares two (partial) color graph codes.
!
!  Discussion:
!
!    CODE1 = CODE2 if every digit of both codes is equal.
!
!    Otherwise, we consider the entries in a special order:
!
!      if ( CODE1(I,J) < CODE2(I,J) ) then we say
!        CODE1 < CODE2
!      else
!        CODE2 < CODE1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) CODE1(NNODE,NNODE), CODE2(NNODE,NNODE), codes to compare.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) NPART, specifies the portion of the codes
!    to compare.  NPART should be between 1 and NNODE.
!
!    If NPART = NNODE, then the full codes are compared.
!
!    If NPART < NNODE, then only entries corresponding to I and J
!    both less than or equal to NPART will be compared.
!
!    Output, integer ( kind = 4 ) RESULT:
!    -1, CODE1 < CODE2;
!     0, CODE1 = CODE2;
!    +1, CODE2 < CODE1.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) code1(nnode,nnode)
  integer ( kind = 4 ) code2(nnode,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) npart
  integer ( kind = 4 ) result

  do j = 1, npart

    do i = 1, j

      if ( code1(i,j) < code2(i,j) ) then

        result = - 1
        return

      else if ( code2(i,j) < code1(i,j) ) then

        result = + 1
        return

      end if

    end do

  end do

  result = 0

  return
end
subroutine cg_code_print ( nnode, code, title )

!*****************************************************************************80
!
!! CG_CODE_PRINT prints a color graph code.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 September
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) CODE(NNODE,NNODE), the code.
!
!    Input, character ( len = * ) TITLE, a title to be printed.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) code(nnode,nnode)
  integer ( kind = 4 ) i
  character ( len = 80 ) string
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  do i = 1, nnode
    string(i:i) = '.'
  end do

  do i = 1, nnode

    write ( *, '(2x,a,80i1)' ) string(1:i-1), code(i,i:nnode)

  end do

  return
end
subroutine cg_color_count ( adj, nnode, mcolor, ncolor )

!*****************************************************************************80
!
!! CG_COLOR_COUNT counts the number of colors in a color graph.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,I) is the color of node I; otherwise, ADJ(I,J) is positive
!    if there is an edge between node I and node J.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) MCOLOR, the maximum color index.
!
!    Output, integer ( kind = 4 ) NCOLOR, the number of colors.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) colors(nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) mcolor
  integer ( kind = 4 ) ncolor

  mcolor = 0
  do i = 1, nnode
    mcolor = max ( mcolor, adj(i,i) )
  end do

  do i = 1, nnode
    colors(i) = adj(i,i)
  end do

  call i4vec_sort_heap_d ( nnode, colors )

  call i4vec_sorted_unique_count ( nnode, colors, ncolor )

  return
end
subroutine cg_color_sequence ( adj, nnode, seq )

!*****************************************************************************80
!
!! CG_COLOR_SEQUENCE computes the color sequence of a color graph.
!
!  Discussion:
!
!    The color sequence of a color graph is constructed by computing the
!    color of each node, and then ordering these values in decreasing order.
!
!    If two color graphs are isomorphic, they must have the same color sequence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,I) is the color of node I; otherwise, ADJ(I,J) is positive
!    if there is an edge between node I and node J.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) SEQ(NNODE), the color sequence.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) seq(nnode)

  do i = 1, nnode
    seq(i) = adj(i,i)
  end do

  call i4vec_sort_heap_d ( nnode, seq )

  return
end
subroutine cg_compare ( adj1, nnode1, adj2, nnode2, order1, &
  order2, result )

!*****************************************************************************80
!
!! CG_COMPARE determines if color graphs G1 and G2 are isomorphic.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ1(NNODE1,NNODE1), the adjacency information for G1.
!    ADJ1(I,I) is the color of node I; otherwise, ADJ1(I,J) is positive
!    if there is an edge between node I and node J.
!
!    Input, integer ( kind = 4 ) NNODE1, the number of nodes in G1.
!
!    Input, integer ( kind = 4 ) ADJ2(NNODE2,NNODE2), the adjacency information for G2.
!    ADJ2(I,I) is the color of node I; otherwise, ADJ2(I,J) is positive
!    if there is an edge between node I and node J.
!
!    Input, integer ( kind = 4 ) NNODE2, the number of nodes in G2.
!
!    Output, integer ( kind = 4 ) RESULT, is 0 if the color graphs are isomorphic,
!    -I if G1 < G2 for test #I, and
!    +I if G2 < G1 for test #I.
!
!    Output, integer ( kind = 4 ) ORDER1(NNODE1), ORDER2(NNODE2).  If RESULT = 0, then
!    ORDER1 and ORDER2 are reorderings of the nodes of G1 and
!    G2 which exhibit the isomorphism.
!
  implicit none

  integer ( kind = 4 ) nnode1
  integer ( kind = 4 ) nnode2

  integer ( kind = 4 ) adj1(nnode1,nnode1)
  integer ( kind = 4 ) adj2(nnode2,nnode2)
  integer ( kind = 4 ) code1(nnode1,nnode1)
  integer ( kind = 4 ) code2(nnode2,nnode2)
  integer ( kind = 4 ) mcolor1
  integer ( kind = 4 ) mcolor2
  integer ( kind = 4 ) ncolor1
  integer ( kind = 4 ) ncolor2
  integer ( kind = 4 ) nedge1
  integer ( kind = 4 ) nedge2
  integer ( kind = 4 ) order1(nnode1)
  integer ( kind = 4 ) order2(nnode2)
  integer ( kind = 4 ) result
  integer ( kind = 4 ) seq1(nnode1)
  integer ( kind = 4 ) seq2(nnode2)
!
!  Test 1: Count the nodes.
!
  if ( nnode1 < nnode2 ) then
    result = - 1
    return
  else if ( nnode2 < nnode1 ) then
    result = + 1
    return
  end if
!
!  Test 2: Count the edges.
!
  call cg_edge_count ( adj1, nnode1, nedge1 )

  call cg_edge_count ( adj2, nnode2, nedge2 )

  if ( nedge1 < nedge2 ) then
    result = - 2
    return
  else if ( nedge2 < nedge1 ) then
    result = + 2
    return
  end if
!
!  Tests 3 and 4: Count the colors, and note the maximum color.
!
  call cg_color_count ( adj1, nnode1, mcolor1, ncolor1 )

  call cg_color_count ( adj2, nnode2, mcolor2, ncolor2 )

  if ( ncolor1 < ncolor2 ) then
    result = - 3
    return
  else if ( ncolor2 < ncolor1 ) then
    result = + 3
    return
  end if

  if ( mcolor1 < mcolor2 ) then
    result = - 4
    return
  else if ( mcolor2 < mcolor1 ) then
    result = + 4
    return
  end if
!
!  Test 5: Compare the degree sequences.
!
  call cg_degree_seq ( adj1, nnode1, seq1 )

  call cg_degree_seq ( adj2, nnode2, seq2 )

  call i4vec_compare ( nnode1, seq1, seq2, result )

  if ( result < 0 ) then
    result = - 5
    return
  else if ( 0 < result ) then
    result = + 5
    return
  end if
!
!  Test 6: Compare the color sequences.
!
  call cg_color_sequence ( adj1, nnode1, seq1 )

  call cg_color_sequence ( adj2, nnode2, seq2 )

  call i4vec_compare ( nnode1, seq1, seq2, result )

  if ( result < 0 ) then
    result = - 6
    return
  else if ( 0 < result ) then
    result = + 6
    return
  end if
!
!  Test 7: Compare the codes.
!
  call cg_code_back ( adj1, nnode1, code1, order1 )

  call cg_code_back ( adj2, nnode2, code2, order2 )

  call cg_code_compare ( code1, code2, nnode1, nnode1, result )

  if ( result < 0 ) then
    result = - 7
    return
  else if ( 0 < result ) then
    result = + 7
    return
  end if

  result = 0

  return
end
subroutine cg_connect_random ( adj, nnode, ncolor, nedge, seed )

!*****************************************************************************80
!
!! CG_CONNECT_RANDOM generates a random connected color graph.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,I) is the color of node I; otherwise, ADJ(I,J) is positive
!    if there is an edge between node I and node J.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) NCOLOR, the number of colors.
!    NCOLOR must be at least 1, and no more than NNODE.
!
!    Input, integer ( kind = 4 ) NEDGE, the number of edges, which must be between
!    NNODE-1 and (NNODE*(NNODE-1))/2.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
  implicit none

  integer ( kind = 4 ) ncolor
  integer ( kind = 4 ) nnode
  integer ( kind = 4 ) nedge

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) code(nnode-2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) icolor
  integer ( kind = 4 ) inode(nnode-1)
  integer ( kind = 4 ) iwork(nedge)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jnode(nnode-1)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) maxedge
  integer ( kind = 4 ) nchoice
  integer ( kind = 4 ) nchoose
  integer ( kind = 4 ) perm(ncolor)
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) subset(ncolor)
!
!  Check.
!
  if ( nnode <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CG_CONNECT_RANDOM - Fatal error!'
    write ( *, '(a,i8)' ) '  NNODE = ', nedge
    write ( *, '(a)' ) '  but NNODE must be at least 1.'
    stop
  end if

  maxedge = ( nnode * ( nnode - 1 ) ) / 2

  if ( nedge < nnode-1 .or. maxedge < nedge ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CG_CONNECT_RANDOM - Fatal error!'
    write ( *, '(a,i8)' ) '  NEDGE = ', nedge
    write ( *, '(a)' ) '  but NEDGE must be at least 0, and '
    write ( *, '(a,i8)' ) '  no more than ', maxedge
    stop
  end if

  if ( ncolor < 1 .or. nnode < ncolor ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CG_CONNECT_RANDOM - Fatal error!'
    write ( *, '(a,i8)' ) '  NCOLOR = ', ncolor
    write ( *, '(a)' ) '  but NCOLOR must be at least 1, and '
    write ( *, '(a,i8)') '  no more than ', nnode
    stop
  end if
!
!  Initialize the adjacency matrix.
!
  adj(1:nnode,1:nnode) = 0
!
!  Choose the colors.
!
  call ksub_random ( nnode, ncolor, subset, seed )

  call perm_random ( ncolor, perm, seed )

  do icolor = 1, ncolor
    i = subset(perm(icolor))
    adj(i,i) = icolor
  end do

  do i = 1, nnode
    if ( adj(i,i) == 0 ) then
      adj(i,i) = i4_uniform ( 1, ncolor, seed )
    end if
  end do
!
!  Pick a random tree.
!
  call tree_arc_random ( nnode, code, inode, jnode, seed )
!
!  Convert information to adjacency form.
!
  call g_arc_to_g_adj ( nnode-1, inode, jnode, adj, nnode )
!
!  Now we have NEDGE - ( NNODE - 1 ) more edges to add.
!
  nchoice = ( nnode * ( nnode - 1 ) ) / 2 - ( nnode - 1 )
  nchoose = nedge - ( nnode - 1 )

  call ksub_random ( nchoice, nchoose, iwork, seed )

  k = 0
  l = 1
  do i = 1, nnode
    do j = i + 1, nnode
      if ( adj(i,j) /= 0 ) then
        k = k + 1

        if ( l <= nchoose ) then
          if ( iwork(l) == k ) then
            adj(i,j) = 1
            adj(j,i) = 1
            l = l + 1
          end if
        end if

      end if
    end do
  end do

  return
end
subroutine cg_degree ( adj, nnode, degree )

!*****************************************************************************80
!
!! CG_DEGREE computes the degree of each node.
!
!  Discussion:
!
!    The degree of a node is the number of edges that are incident on it.
!    The sum of the degrees of the nodes is twice the number of edges.
!
!    The generalized case, where ADJ(I,J) can be greater than 1, indicating
!    the existence of 2 or more distinct edges between nodes I and J,
!    will be properly handled by this routine.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,I) is the color of node I; otherwise, ADJ(I,J) is positive
!    if there is an edge between node I and node J.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) DEGREE(NNODE), the degree of the nodes.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) degree(nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  degree(1:nnode) = 0

  do i = 1, nnode
    do j = 1, nnode
      if ( i /= j ) then
        if ( adj(i,j) /= 0 ) then
          degree(i) = degree(i) + adj(i,j)
        end if
      end if
    end do
  end do

  return
end
subroutine cg_degree_seq ( adj, nnode, seq )

!*****************************************************************************80
!
!! CG_DEGREE_SEQ computes the degree sequence of a color graph.
!
!  Discussion:
!
!    The degree sequence of a color graph is constructed by computing the
!    degree of each node, and then ordering these values in decreasing order.
!
!    If two color graphs are isomorphic, they must have the same
!    degree sequence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,I) is the color of node I; otherwise, ADJ(I,J) is positive
!    if there is an edge between node I and node J.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) SEQ(NNODE), the degree sequence.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) seq(nnode)

  call cg_degree ( adj, nnode, seq )

  call i4vec_sort_heap_d ( nnode, seq )

  return
end
subroutine cg_edge_count ( adj, nnode, nedge )

!*****************************************************************************80
!
!! CG_EDGE_COUNT counts the number of edges in a color graph.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 October 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,I) is the color of node I; otherwise, ADJ(I,J) is positive
!    if there is an edge between node I and node J.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) NEDGE, the number of edges.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nedge

  nedge = 0

  do i = 1, nnode
    do j = 1, nnode

      if ( i /= j ) then
        nedge = nedge + adj(i,j)
      end if

    end do
  end do

  nedge = nedge / 2

  return
end
subroutine cg_example_bush ( adj, nnode )

!*****************************************************************************80
!
!! CG_EXAMPLE_BUSH sets up the bush color graph.
!
!  Diagram:
!
!        6G  3R
!        |   |
!    1B--4G--5W--2R
!        |
!        7W
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,I) is the color of node I; otherwise, ADJ(I,J) is positive
!    if there is an edge between node I and node J.
!
!    Output, integer ( kind = 4 ) NNODE, the number of nodes.
!
  implicit none

  integer ( kind = 4 ) adj(7,7)
  integer ( kind = 4 ), parameter :: BLUE = 1
  integer ( kind = 4 ), parameter :: GREEN = 2
  integer ( kind = 4 ) nnode
  integer ( kind = 4 ), parameter :: RED = 3
  integer ( kind = 4 ), parameter :: WHITE = 4

  nnode = 7

  adj(1:nnode,1:nnode) = 0

  adj(1,1) = BLUE
  adj(1,4) = 1

  adj(2,2) = RED
  adj(2,5) = 1

  adj(3,3) = RED
  adj(3,5) = 1

  adj(4,1) = 1
  adj(4,4) = GREEN
  adj(4,5) = 1
  adj(4,6) = 1
  adj(4,7) = 1

  adj(5,2) = 1
  adj(5,3) = 1
  adj(5,4) = 1
  adj(5,5) = WHITE

  adj(6,4) = 1
  adj(6,6) = GREEN

  adj(7,4) = 1
  adj(7,7) = WHITE

  return
end
subroutine cg_example_cube ( adj, nnode )

!*****************************************************************************80
!
!! CG_EXAMPLE_CUBE sets up the cube color graph.
!
!  Diagram:
!
!      4R----7R
!     /|    /|
!    8B----3B|
!    | |   | |
!    | 5B--|-2G
!    |/    |/
!    1G----6B
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,I) is the color of node I; otherwise, ADJ(I,J) is positive
!    if there is an edge between node I and node J.
!
!    Output, integer ( kind = 4 ) NNODE, the number of nodes, which is 8.
!
  implicit none

  integer ( kind = 4 ) adj(8,8)
  integer ( kind = 4 ), parameter :: BLUE = 1
  integer ( kind = 4 ), parameter :: GREEN = 2
  integer ( kind = 4 ) nnode
  integer ( kind = 4 ), parameter :: RED = 3

  nnode = 8

  adj(1:nnode,1:nnode) = 0

  adj(1,1) = GREEN
  adj(1,5) = 1
  adj(1,6) = 1
  adj(1,8) = 1

  adj(2,2) = GREEN
  adj(2,5) = 1
  adj(2,6) = 1
  adj(2,7) = 1

  adj(3,3) = BLUE
  adj(3,6) = 1
  adj(3,7) = 1
  adj(3,8) = 1

  adj(4,4) = RED
  adj(4,5) = 1
  adj(4,7) = 1
  adj(4,8) = 1

  adj(5,5) = BLUE
  adj(5,1) = 1
  adj(5,2) = 1
  adj(5,4) = 1

  adj(6,6) = BLUE
  adj(6,1) = 1
  adj(6,2) = 1
  adj(6,3) = 1

  adj(7,7) = RED
  adj(7,2) = 1
  adj(7,3) = 1
  adj(7,4) = 1

  adj(8,8) = BLUE
  adj(8,1) = 1
  adj(8,3) = 1
  adj(8,4) = 1

  return
end
subroutine cg_example_octo ( example, adj, nnode, seed )

!*****************************************************************************80
!
!! CG_EXAMPLE_OCTO sets up an 8 node example color graph.
!
!  Diagram:
!
!      1---2
!     /|   |\
!    8-+---+-3
!    | |   | |
!    7-+---+-4
!     \|   |/
!      6---5
!
!     Graph "A"
!
!    There are 8 graphs to choose from.  They are all on 8 nodes.  The first
!    5 have degree 3 at every node.  Graphs 6 and 7 have degree 5 at every
!    node.  Graph 8 is disconnected.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) EXAMPLE, the index of the example, between 1 and 40.
!
!    Output, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,I) is the color of node I.
!    ADJ(I,J) is 1 if nodes I and J are adjacent and 0 otherwise.
!
!    Output, integer ( kind = 4 ) NNODE, the number of nodes, which should be 8.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
  implicit none

  integer ( kind = 4 ) adj(8,8)
  integer ( kind = 4 ), parameter :: BLUE = 1
  integer ( kind = 4 ) example
  integer ( kind = 4 ), parameter :: GREEN = 2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) j
  integer ( kind = 4 ) msave
  integer ( kind = 4 ) nnode
  integer ( kind = 4 ) nsave
  integer ( kind = 4 ), parameter :: RED = 3
  integer ( kind = 4 ) seed
  integer ( kind = 4 ), parameter :: YELLOW = 4

  if ( example <= 0 ) then
    nsave = i4_uniform ( 1, 8, seed )
    msave = i4_uniform ( 1, 5, seed )
  else
    example = mod ( example - 1, 40 ) + 1
    msave = ( ( example - 1 ) / 8 ) + 1
    nsave = mod ( example - 1, 8 ) + 1
  end if

  nnode = 8

  adj(1:nnode,1:nnode) = 0

  do i = 1, nnode
    j = i + 1
    if ( nnode < j ) then
      j = j - nnode
    end if

    adj(i,j) = 1
    adj(j,i) = 1

  end do
!
!  Underlying graph 1.
!
  if ( nsave == 1 ) then

    adj(1,6) = 1
    adj(6,1) = 1
    adj(2,5) = 1
    adj(5,2) = 1
    adj(3,8) = 1
    adj(8,3) = 1
    adj(4,7) = 1
    adj(7,4) = 1
!
!  Underlying graph 2.
!
  else if ( nsave == 2 ) then

    adj(1,6) = 1
    adj(6,1) = 1
    adj(2,8) = 1
    adj(8,2) = 1
    adj(3,5) = 1
    adj(5,3) = 1
    adj(4,7) = 1
    adj(7,4) = 1
!
!  Underlying graph 3.
!
  else if ( nsave == 3 ) then

    adj(1,5) = 1
    adj(5,1) = 1
    adj(2,6) = 1
    adj(6,2) = 1
    adj(3,7) = 1
    adj(7,3) = 1
    adj(4,8) = 1
    adj(8,4) = 1
!
!  Underlying graph 4.
!
  else if ( nsave == 4 ) then

    adj(1,3) = 1
    adj(3,1) = 1
    adj(2,4) = 1
    adj(4,2) = 1
    adj(5,7) = 1
    adj(7,5) = 1
    adj(6,8) = 1
    adj(8,6) = 1
!
!  Underlying graph 5.
!
  else if ( nsave == 5 ) then

    adj(1,4) = 1
    adj(4,1) = 1
    adj(2,6) = 1
    adj(6,2) = 1
    adj(3,8) = 1
    adj(8,3) = 1
    adj(5,7) = 1
    adj(7,5) = 1
!
!  Underlying graph 6.
!
  else if ( nsave == 6 ) then

    adj(1,4) = 1
    adj(4,1) = 1
    adj(1,5) = 1
    adj(5,1) = 1
    adj(1,6) = 1
    adj(6,1) = 1
    adj(2,5) = 1
    adj(5,2) = 1
    adj(2,6) = 1
    adj(6,2) = 1
    adj(2,7) = 1
    adj(7,2) = 1
    adj(3,6) = 1
    adj(6,3) = 1
    adj(3,7) = 1
    adj(7,3) = 1
    adj(3,8) = 1
    adj(8,3) = 1
    adj(4,7) = 1
    adj(7,4) = 1
    adj(4,8) = 1
    adj(8,4) = 1
    adj(5,8) = 1
    adj(8,5) = 1
!
!  Underlying graph 7.
!
  else if ( nsave == 7 ) then

    adj(1,3) = 1
    adj(3,1) = 1
    adj(1,5) = 1
    adj(5,1) = 1
    adj(1,7) = 1
    adj(7,1) = 1
    adj(2,4) = 1
    adj(4,2) = 1
    adj(2,6) = 1
    adj(6,2) = 1
    adj(2,8) = 1
    adj(8,2) = 1
    adj(3,5) = 1
    adj(5,3) = 1
    adj(3,7) = 1
    adj(7,3) = 1
    adj(4,6) = 1
    adj(6,4) = 1
    adj(4,8) = 1
    adj(8,4) = 1
    adj(5,7) = 1
    adj(7,5) = 1
    adj(6,8) = 1
    adj(8,6) = 1

  else if ( nsave == 8 ) then

    adj(1,2) = 1
    adj(2,1) = 1
    adj(1,3) = 1
    adj(3,1) = 1
    adj(2,3) = 1
    adj(3,2) = 1
    adj(3,4) = 1
    adj(4,3) = 1
    adj(5,6) = 1
    adj(6,5) = 1
    adj(5,7) = 1
    adj(7,5) = 1
    adj(6,7) = 1
    adj(7,6) = 1

  end if

  if ( msave == 1 ) then

    adj(1,1) = RED
    adj(2,2) = RED
    adj(3,3) = RED
    adj(4,4) = BLUE
    adj(5,5) = BLUE
    adj(6,6) = BLUE
    adj(7,7) = GREEN
    adj(8,8) = GREEN

  else if ( msave == 2 ) then

    adj(1,1) = RED
    adj(2,2) = RED
    adj(3,3) = RED
    adj(4,4) = BLUE
    adj(5,5) = BLUE
    adj(6,6) = BLUE
    adj(7,7) = GREEN
    adj(8,8) = YELLOW

  else if ( msave == 3 ) then

    adj(1,1) = RED
    adj(2,2) = RED
    adj(3,3) = RED
    adj(4,4) = BLUE
    adj(5,5) = BLUE
    adj(6,6) = BLUE
    adj(7,7) = YELLOW
    adj(8,8) = YELLOW

  else if ( msave == 4 ) then

    adj(1,1) = RED
    adj(2,2) = RED
    adj(3,3) = RED
    adj(4,4) = BLUE
    adj(5,5) = BLUE
    adj(6,6) = GREEN
    adj(7,7) = GREEN
    adj(8,8) = GREEN

  else if ( msave == 5 ) then

    adj(1,1) = RED
    adj(2,2) = BLUE
    adj(3,3) = RED
    adj(4,4) = GREEN
    adj(5,5) = BLUE
    adj(6,6) = RED
    adj(7,7) = BLUE
    adj(8,8) = GREEN

  end if
!
!  Now permute the graph.
!
  call i4mat_perm_random ( nnode, adj, seed )

  return
end
subroutine cg_example_twig ( adj, nnode )

!*****************************************************************************80
!
!! CG_EXAMPLE_TWIG sets up the twig color graph.
!
!  Diagram:
!
!    1R---2R---3B
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,I) is the color of node I; otherwise, ADJ(I,J) is positive
!    if there is an edge between node I and node J.
!
!    Output, integer ( kind = 4 ) NNODE, the number of nodes.
!
  implicit none

  integer ( kind = 4 ) adj(3,3)
  integer ( kind = 4 ), parameter :: BLUE = 1
  integer ( kind = 4 ) nnode
  integer ( kind = 4 ), parameter :: RED = 3

  nnode = 3

  adj(1:nnode,1:nnode) = 0

  adj(1,1) = RED
  adj(1,2) = 1

  adj(2,1) = 1
  adj(2,2) = RED
  adj(2,3) = 1

  adj(3,2) = 1
  adj(3,3) = BLUE

  return
end
subroutine cg_order_code ( adj, nnode, npart, code, order )

!*****************************************************************************80
!
!! CG_ORDER_CODE returns the color graph code for a specific node ordering.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,I) is the color of node I; otherwise, ADJ(I,J) is positive
!    if there is an edge between node I and node J.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) NPART, the number of nodes to consider.
!    NPART should be between 1 and NNODE.
!
!    If NPART is NNODE, then the full code is returned.
!
!    If NPART is less than NNODE, then the code is computed as
!    though only NPART nodes existed, namely, those specified in the
!    first NPART entries of order.  This option is provided so that
!    the routine can compute the portion of a code determined
!    by an incomplete ordering of the nodes.
!
!    Output, integer ( kind = 4 ) CODE(NNODE,NNODE), the code for this ordering.
!
!    Input, integer ( kind = 4 ) ORDER(NNODE), the ordering of the nodes.  ORDER(1)
!    is the first node, and so on.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) code(nnode,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ni
  integer ( kind = 4 ) nj
  integer ( kind = 4 ) npart
  integer ( kind = 4 ) order(nnode)

  do i = 1, nnode

    if ( i <= npart ) then

      ni = order(i)

      if ( order(i) < 1 .or. nnode < order(i) ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CG_ORDER_CODE - Fatal error!'
        write ( *, '(a)' ) '  ORDER is not a proper permutation.'
        stop
      end if

    else
      ni = 0
    end if

    do j = i, nnode

      if ( j <= npart ) then

        nj = order(j)

        if ( order(j) < 1 .or. nnode < order(j) ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'CG_ORDER_CODE - Fatal error!'
          write ( *, '(a)' ) '  ORDER is not a proper permutation.'
          stop
        end if

      else
        nj = 0
      end if

      if ( ni == 0 .or. nj == 0 ) then

        code(i,j) = 0

      else if ( ni <= nj ) then

        code(i,j) = adj(ni,nj)

      else

        code(i,j) = adj(nj,ni)

      end if

    end do
  end do

  return
end
subroutine cg_print ( adj, nnode, title )

!*****************************************************************************80
!
!! CG_PRINT prints out the adjacency matrix of a color graph.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,I) is the color of node I; otherwise, ADJ(I,J) is positive
!    if there is an edge between node I and node J.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, character ( len = * ) TITLE, a title to be printed.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  character ( len = 80 ) string
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  do i = 1, nnode

    do j = 1, nnode

      k = (j-1) * 3 + 1
      write ( string(k:k+2), '(i3)' ) adj(i,j)

    end do

    write ( *, '(2x,a)' ) string(1:3*nnode)

  end do

  return
end
subroutine cg_random ( adj, nnode, ncolor, nedge, seed )

!*****************************************************************************80
!
!! CG_RANDOM generates a random color graph.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,I) is the color of node I; otherwise, ADJ(I,J) is positive
!    if there is an edge between node I and node J.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) NCOLOR, the number of colors.
!    NCOLOR must be at least 1, and no more than NNODE.
!
!    Input, integer ( kind = 4 ) NEDGE, the number of edges, which must be between
!    0 and (NNODE*(NNODE-1))/2.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
  implicit none

  integer ( kind = 4 ) nnode
  integer ( kind = 4 ) nedge

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) icolor
  integer ( kind = 4 ) iwork(nedge)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) maxedge
  integer ( kind = 4 ) ncolor
  integer ( kind = 4 ) perm(ncolor)
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) subset(ncolor)

  if ( nnode <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CG_RANDOM - Fatal error!'
    write ( *, '(a,i8)' ) '  NNODE = ', nedge
    write ( *, '(a)' ) '  but NNODE must be at least 1.'
    stop
  end if

  maxedge = ( nnode * ( nnode - 1 ) ) / 2

  if ( nedge < 0 .or. maxedge < nedge ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CG_RANDOM - Fatal error!'
    write ( *, '(a,i8)' ) '  NEDGE = ', nedge
    write ( *, '(a)' ) '  but NEDGE must be at least 0, and '
    write ( *, '(a,i8)' ) '  no more than ', maxedge
    stop
  end if

  if ( ncolor < 1 .or. nnode < ncolor ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CG_RANDOM - Fatal error!'
    write ( *, '(a)' ) '  Illegal value of NCOLOR.'
    stop
  end if
!
!  Start out with no edges and no colors.
!
  adj(1:nnode,1:nnode) = 0
!
!  Choose the colors.
!
  call ksub_random ( nnode, ncolor, subset, seed )

  call perm_random ( ncolor, perm, seed )

  do icolor = 1, ncolor
    i = subset(perm(icolor))
    adj(i,i) = icolor
  end do

  do i = 1, nnode
    if ( adj(i,i) == 0 ) then
      adj(i,i) = i4_uniform ( 1, ncolor, seed )
    end if
  end do
!
!  Pick a random NEDGE subset of (N*(N-1))/2.
!
  call ksub_random ( maxedge, nedge, iwork, seed )
!
!  The (n*(n-1))/2 spots in the superdiagonal are numbered as follows:
!
!  * 1  2   3  ...  n-1   n
!  * * n+1 n+2 ... 2n-2  2n-1
!  ...
!  * *  *   *  ...   *   (n*(n-1))/2
!  * *  *   *  ...   *    *
!
  k = 0
  l = 1
  do i = 1, nnode-1
    do j = i+1, nnode

      k = k + 1
      if ( l <= nedge ) then

        if ( k == iwork(l) ) then
          adj(i,j) = 1
          adj(j,i) = 1
          l = l + 1
        end if

      end if

    end do
  end do

  return
end
subroutine dg_code_back ( adj, nnode, code, order )

!*****************************************************************************80
!
!! DG_CODE_BACK computes a digraph code via backtracking.
!
!  Discussion:
!
!    The code is the "largest" order code (in the lexicographic
!    sense) over all possible node orderings.  The backtracking method
!    organizes the search of all possible node orderings so that if
!    a partial node ordering is sure to generate an order code
!    less than the best so far, then all the orderings that begin with
!    this partial ordering are immediately dropped from consideration.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,J) is 1 if there is a direct link from node I to node J.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) CODE(NNODE,NNODE), the code.
!
!    Output, integer ( kind = 4 ) ORDER(NNODE), the ordering of the nodes.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) bestcode(nnode,nnode)
  integer ( kind = 4 ) bestorder(nnode)
  integer ( kind = 4 ) code(nnode,nnode)
  logical, parameter :: debug = .false.
  integer ( kind = 4 ) index
  integer ( kind = 4 ) maxstack
  integer ( kind = 4 ) ncan(nnode)
  integer ( kind = 4 ) ncomp
  integer ( kind = 4 ) nopen
  integer ( kind = 4 ) nstack
  integer ( kind = 4 ) nswap
  integer ( kind = 4 ) order(nnode)
  integer ( kind = 4 ) result
  integer ( kind = 4 ) stack(4*nnode)

  if ( nnode <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DG_CODE_BACK - Fatal error!'
    write ( *, '(a,i8)' ) '  NNODE = ', nnode
    stop
  end if

  maxstack = 4 * nnode
  nstack = 0
  stack(1) = 0

  ncomp = 0
  nswap = 0
!
!  Start with the identity ordering.
!
  call i4vec_indicator ( nnode, order )
!
!  Compute the corresponding code.
!
  call dg_order_code ( adj, nnode, nnode, code, order )
!
!  Take this ordering and code as the best so far.
!
  bestorder(1:nnode) = order(1:nnode)
  bestcode(1:nnode,1:nnode) = code(1:nnode,1:nnode)
!
!  Now consider all possible orderings, and their codes.
!
  index = 0

  do

    call i4vec_backtrack ( nnode, order, index, nopen, nstack, stack, &
      maxstack, ncan )
!
!  If the backtrack routine has returned a complete candidate ordering, then
!  compute the resulting code, and see it it is better
!  then our current best.  Then go back for the next backtrack search.
!
    if ( index == 1 ) then

      call dg_order_code ( adj, nnode, nnode, code, order )

      call dg_code_compare ( bestcode, code, nnode, nnode, result )

      ncomp = ncomp + 1

      if ( result == - 1 ) then

        nswap = nswap + 1

        bestorder(1:nnode) = order(1:nnode)
        bestcode(1:nnode,1:nnode) = code(1:nnode,1:nnode)

      end if
!
!  If the backtrack routine has a partial reordering,
!  supply candidates for the next item in the ordering.
!
    else if ( index == 2 ) then

      call dg_code_cand ( adj, bestcode, code, nnode, ncomp, nopen, &
        order, stack, maxstack, nstack, ncan )

    else

      exit

    end if

  end do
!
!  Once we have examined all possibilites, we are done.
!
!  Set the output ordering to the best ordering, and the output
!  code to the corresponding best code.
!
  order(1:nnode) = bestorder(1:nnode)
  code(1:nnode,1:nnode) = bestcode(1:nnode,1:nnode)

  if ( debug ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DG_CODE_BACK:'
    write ( *, '(a,i8)' ) '  Comparisons: ', ncomp
    write ( *, '(a,i8)' ) '  Swaps:	  ', nswap
  end if

  return
end
subroutine dg_code_brute ( adj, nnode, code, order )

!*****************************************************************************80
!
!! DG_CODE_BRUTE computes a digraph code via brute force.
!
!  Discussion:
!
!    The code is the "largest" order code in the lexicographic
!    sense over all node orderings.  The brute force method
!    considers every node ordering, computes the corresponding
!    order code, and takes the largest one encountered.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,J) is 1 if there is a direct link from node I to node J.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) CODE(NNODE,NNODE), the code.
!
!    Output, integer ( kind = 4 ) ORDER(NNODE), the ordering of the nodes.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) bestcode(nnode,nnode)
  integer ( kind = 4 ) bestorder(nnode)
  integer ( kind = 4 ) code(nnode,nnode)
  logical even
  logical more
  integer ( kind = 4 ) ncomp
  integer ( kind = 4 ) nswap
  integer ( kind = 4 ) order(nnode)
  integer ( kind = 4 ) result

  ncomp = 0
  nswap = 0
!
!  Start with the identity ordering.
!
  call i4vec_indicator ( nnode, order )
!
!  Compute the corresponding code.
!
  call dg_order_code ( adj, nnode, nnode, code, order )
!
!  Take this ordering and code as the best so far.
!
  bestorder(1:nnode) = order(1:nnode)
  bestcode(1:nnode,1:nnode) = code(1:nnode,1:nnode)
!
!  Now consider all possible orderings, and their codes.
!
  more = .false.

  do

    call perm_next ( nnode, order, more, even )

    call dg_order_code ( adj, nnode, nnode, code, order )

    call dg_code_compare ( bestcode, code, nnode, nnode, result )

    ncomp = ncomp + 1

    if ( result == - 1 ) then

      nswap = nswap + 1

      bestorder(1:nnode) = order(1:nnode)
      bestcode(1:nnode,1:nnode) = code(1:nnode,1:nnode)

    end if

    if ( .not. more ) then
      exit
    end if

  end do
!
!  Once we have examined all possibilites, we are done.
!
!  Set the output ordering to the best ordering, and the output
!  code to the corresponding best code.
!
  order(1:nnode) = bestorder(1:nnode)
  code(1:nnode,1:nnode) = bestcode(1:nnode,1:nnode)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DG_CODE_BRUTE:'
  write ( *, '(a,i8)' ) '  Comparisons: ', ncomp
  write ( *, '(a,i8)' ) '  Swaps:       ', nswap

  return
end
subroutine dg_code_cand ( adj, bestcode, code, nnode, ncomp, nopen, &
  order, stack, maxstack, nstack, ncan )

!*****************************************************************************80
!
!! DG_CODE_CAND finds candidates for a maximal digraph code ordering.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,J) is 1 if there is a direct link from node I to node J.
!
!    Input, integer ( kind = 4 ) BESTCODE(NNODE,NNODE), the best code so far.
!
!    Workspace, integer CODE(NNODE,NNODE).
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input/output, integer ( kind = 4 ) NCOMP, the number of code comparisons.
!    This routine updates NCOMP by 1 each time the routine is called.
!
!    Input, integer ( kind = 4 ) NOPEN, identifies the first open position in ORDER.
!
!    Input, integer ( kind = 4 ) ORDER(NNODE), contains in entries 1 through NOPEN-1
!    the elements of the current partial list.
!
!    Input/output, integer ( kind = 4 ) STACK(MAXSTACK), used to store the new candidates.
!
!    Input, integer ( kind = 4 ) MAXSTACK, the maximum size of the STACK array.
!
!    Input/output, integer ( kind = 4 ) NSTACK, the current length of the stack.
!    On output, NSTACK has been increased by the number of
!    candidates found.
!
!    Input/output, integer ( kind = 4 ) NCAN(NNODE), the number of candidates for
!    each position.
!
  implicit none

  integer ( kind = 4 ) maxstack
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) bestcode(nnode,nnode)
  integer ( kind = 4 ) code(nnode,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ifree(nnode)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) max_adj(nnode)
  integer ( kind = 4 ) max_max_adj
  integer ( kind = 4 ) ncan(nnode)
  integer ( kind = 4 ) ncomp
  integer ( kind = 4 ) nfree
  integer ( kind = 4 ) ni
  integer ( kind = 4 ) nj
  integer ( kind = 4 ) nopen
  integer ( kind = 4 ) nstack
  integer ( kind = 4 ) order(nnode)
  integer ( kind = 4 ) result
  integer ( kind = 4 ) stack(maxstack)

  if ( nopen < 1 .or. nnode < nopen ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DG_CODE_CAND - Fatal error 1!'
    write ( *, '(a)' ) '  1 <= NOPEN <= NNODE should be true, but'
    write ( *, '(a,i8)' ) '  NOPEN = ', nopen
    write ( *, '(a,i8)' ) '  NNODE = ', nnode
    stop
  end if
!
!  Start with no candidates.
!
  ncan(nopen) = 0
!
!  If we have fixed at least one entry of the list,
!
!    compute the partial code;
!
!    compare the partial code with the corresponding
!    part of the code for the best ordering so far;
!
!    If the current incomplete ordering is actually LESS than the
!    current best, then bail out now, with zero candidates.
!
  if ( 1 < nopen ) then

    call dg_order_code ( adj, nnode, nopen-1, code, order )

    call dg_code_compare ( bestcode, code, nnode, nopen-1, result )

    ncomp = ncomp + 1

    if ( result == + 1 ) then
      ncan(nopen) = 0
      return
    end if

  end if
!
!  Get a list of those nodes which have not been used yet.
!
  nfree = nnode + 1 - nopen
  call perm_free ( order, nopen-1, ifree, nfree )
!
!  Our preferred candidates will be the unused neighbors of the
!  lowest ordered node possible.
!
  ncan(nopen) = 0

  do i = 1, nopen-1

    ncan(nopen) = 0

    ni = order(i)
!
!  First: look for neighbors FROM NI.
!
    do j = 1, nfree

      nj = ifree(j)

      if ( adj(ni,nj) /= 0 ) then
        ncan(nopen) = ncan(nopen) + 1
        nstack = nstack + 1

        if ( maxstack < nstack ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'DG_CODE_CAND - Fatal error 2!'
          write ( *, '(a)' ) '  MAXSTACK < NSTACK!'
          write ( *, '(a,i8)' ) '  NSTACK =   ', nstack
          write ( *, '(a,i8)' ) '  MAXSTACK = ', maxstack
          stop
        end if

        stack(nstack) = nj

      end if

    end do

    if ( 0 < ncan(nopen) ) then
      return
    end if
!
!  Second: look for neighbors TO NI.
!
    do j = 1, nfree

      nj = ifree(j)

      if ( adj(nj,ni) /= 0 ) then

        ncan(nopen) = ncan(nopen) + 1
        nstack = nstack + 1

        if ( maxstack < nstack ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'DG_CODE_CAND - Fatal error 2!'
          write ( *, '(a)' ) '  MAXSTACK < NSTACK!'
          write ( *, '(a,i8)' ) '  NSTACK =   ', nstack
          write ( *, '(a,i8)' ) '  MAXSTACK = ', maxstack
          stop
        end if

        stack(nstack) = nj

      end if

    end do

    if ( 0 < ncan(nopen) ) then
      return
    end if

  end do
!
!  If we get here, no free nodes are connected in any way to
!  used nodes.  This can happen in two ways:
!
!  * NOPEN = 1;
!  * The digraph is disconnected;
!
  max_max_adj = 0

  do i = 1, nfree

    ni = ifree(i)

    max_adj(i) = 0
    do j = 1, nfree
      nj = ifree(j)
      if ( ni /= nj ) then
        max_adj(i) = max ( max_adj(i), adj(ni,nj) )
      end if
    end do

    max_max_adj = max ( max_max_adj, max_adj(i) )

  end do

  ncan(nopen) = 0

  do i = 1, nfree

    if ( max_adj(i) == max_max_adj ) then

      ncan(nopen) = ncan(nopen) + 1
      nstack = nstack + 1

      if ( maxstack < nstack ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'DG_CODE_CAND - Fatal error 2!'
        write ( *, '(a)' ) '  MAXSTACK < NSTACK!'
        write ( *, '(a,i8)' ) '  NSTACK =   ', nstack
        write ( *, '(a,i8)' ) '  MAXSTACK = ', maxstack
        stop
      end if

      stack(nstack) = ifree(i)

    end if
  end do

  return
end
subroutine dg_code_compare ( code1, code2, nnode, npart, result )

!*****************************************************************************80
!
!! DG_CODE_COMPARE compares two partial digraph codes.
!
!  Discussion:
!
!    CODE1 = CODE2 if every digit of both codes is equal.
!
!    Otherwise, traverse the entries in a funny diagonal way, suggested
!    by this diagram for the first 16 entries:
!
!       1  2  5 10
!       3  4  7 12
!       6  8  9 14
!      11 13 15 16
!
!    As we do that, we examine the corresponding digits of the two codes.
!    For the first entry, (I,J), where the two codes differ, we say:
!
!      if ( CODE1(I,J) < CODE2(I,J) ) then we say
!        CODE1 < CODE2
!      else
!        CODE2 < CODE1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) CODE1(NNODE,NNODE), CODE2(NNODE,NNODE), codes to compare.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) NPART, specifies the portion of the codes
!    to compare.  NPART should be between 1 and NNODE.
!
!    If NPART = NNODE, then the full codes are compared.
!
!    If NPART < NNODE, then only entries corresponding to I and J
!    both less than or equal to NPART will be compared.
!
!    Output, integer ( kind = 4 ) RESULT:
!    -1, CODE1 < CODE2;
!     0, CODE1 = CODE2;
!    +1, CODE2 < CODE1.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) code1(nnode,nnode)
  integer ( kind = 4 ) code2(nnode,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) npart
  integer ( kind = 4 ) result

  do j = 2, npart

    do i = 1, j - 1

      if ( code1(i,j) < code2(i,j) ) then

        result = - 1
        return

      else if ( code2(i,j) < code1(i,j) ) then

        result = + 1
        return

      else if ( code1(j,i) < code2(j,i) ) then

        result = - 1
        return

      else if ( code2(j,i) < code1(j,i) ) then

        result = + 1
        return

      end if

    end do

  end do

  result = 0

  return
end
subroutine dg_code_print ( nnode, code, title )

!*****************************************************************************80
!
!! DG_CODE_PRINT prints out a digraph code.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) CODE(NNODE,NNODE), the code.
!
!    Input, character ( len = * ) TITLE, a title to be printed.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) ck
  integer ( kind = 4 ) code(nnode,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  character ( len = 80 ) string
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  do i = 1, nnode

    do j = 1, nnode

      if ( i == j ) then
        string(j:j) = '.'
      else
        ck = code(i,j)
        if ( 0 <= ck .and. ck <= 9 ) then
          string(j:j) = char ( 48 + ck )
        else
          string(j:j) = '*'
        end if
      end if

    end do

    write ( *, '(2x,i4,2x,a)' ) i, string(1:nnode)

  end do

  return
end
subroutine dg_compare ( adj1, nnode1, adj2, nnode2, order1, &
  order2, result )

!*****************************************************************************80
!
!! DG_COMPARE determines if digraphs G1 and G2 are isomorphic.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 October 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ1(NNODE1,NNODE1), the adjacency information for G1.
!
!    Input, integer ( kind = 4 ) NNODE1, the number of nodes in G1.
!
!    Input, integer ( kind = 4 ) ADJ2(NNODE2,NNODE2), the adjacency information for G2.
!
!    Input, integer ( kind = 4 ) NNODE2, the number of nodes in G2.
!
!    Output, integer ( kind = 4 ) RESULT, is 0 if G1 and G2 are isomorphic,
!    -I if G1 < G2 for test #I, and
!    +I if G2 < G1 for test #I.
!
!    Output, integer ( kind = 4 ) ORDER1(NNODE1), ORDER2(NNODE2).  If RESULT = 0, then
!    ORDER1 and ORDER2 are reorderings of the nodes of G1 and G2
!    which exhibit the isomorphism.
!
  implicit none

  integer ( kind = 4 ) nnode1
  integer ( kind = 4 ) nnode2

  integer ( kind = 4 ) adj1(nnode1,nnode1)
  integer ( kind = 4 ) adj2(nnode2,nnode2)
  integer ( kind = 4 ) code1(nnode1,nnode1)
  integer ( kind = 4 ) code2(nnode2,nnode2)
  integer ( kind = 4 ) in_seq1(nnode1)
  integer ( kind = 4 ) in_seq2(nnode2)
  integer ( kind = 4 ) nedge1
  integer ( kind = 4 ) nedge2
  integer ( kind = 4 ) order1(nnode1)
  integer ( kind = 4 ) order2(nnode2)
  integer ( kind = 4 ) out_seq1(nnode1)
  integer ( kind = 4 ) out_seq2(nnode2)
  integer ( kind = 4 ) result
!
!  Test 1: Count the nodes.
!
  if ( nnode1 < nnode2 ) then
    result = - 1
    return
  else if ( nnode2 < nnode1 ) then
    result = + 1
    return
  end if
!
!  Test 2: Count the edges.
!
  call dg_edge_count ( adj1, nnode1, nedge1 )

  call dg_edge_count ( adj2, nnode2, nedge2 )

  if ( nedge1 < nedge2 ) then
    result = - 2
    return
  else if ( nedge2 < nedge1 ) then
    result = + 2
    return
  end if
!
!  Test 3: Compare the outdegree sequences.
!
  call dg_degree_seq ( adj1, nnode1, in_seq1, out_seq1 )

  call dg_degree_seq ( adj2, nnode2, in_seq2, out_seq2 )

  call i4vec_compare ( nnode1, out_seq1, out_seq2, result )

  if ( result < 0 ) then
    result = - 3
    return
  else if ( 0 < result ) then
    result = + 3
    return
  end if
!
!  Test 4: Compare the indegree sequences.
!
  call i4vec_compare ( nnode1, in_seq1, in_seq2, result )

  if ( result < 0 ) then
    result = - 4
    return
  else if ( 0 < result ) then
    result = + 4
    return
  end if
!
!  Test 5: Compare the codes.
!
  call dg_code_back ( adj1, nnode1, code1, order1 )

  call dg_code_back ( adj2, nnode2, code2, order2 )

  call dg_code_compare ( code1, code2, nnode1, nnode1, result )

  if ( result < 0 ) then
    result = - 5
    return
  else if ( 0 < result ) then
    result = + 5
    return
  end if

  result = 0

  return
end
subroutine dg_degree ( adj, nnode, indegree, outdegree )

!*****************************************************************************80
!
!! DG_DEGREE computes the indegree and outdegree of each node.
!
!  Discussion:
!
!    The indegree of a node is the number of directed edges that
!    end at the node.
!
!    The outdegree of a node is the number of directed edges that
!    begin at the node.
!
!    The sum of the indegrees and outdegrees of all the nodes is twice
!    the number of edges.
!
!    The generalized case, where ADJ(I,J) can be greater than 1, indicating
!    the existence of 2 or more distinct edges from node I to node J,
!    will be properly handled by this routine.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 October 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,J) is 1 if there is a direct link from node I to node J.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) INDEGREE(NNODE), OUTDEGREE(NNODE),
!    the indegree and outdegree of the nodes.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indegree(nnode)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) outdegree(nnode)

  indegree(1:nnode) = 0
  outdegree(1:nnode) = 0

  do i = 1, nnode
    do j = 1, nnode
      if ( adj(i,j) /= 0 ) then
        outdegree(i) = outdegree(i) + adj(i,j)
        indegree(j) = indegree(j) + adj(i,j)
      end if
    end do
  end do

  return
end
subroutine dg_degree_max ( adj, nnode, indegree_max, outdegree_max, &
  degree_max )

!*****************************************************************************80
!
!! DG_DEGREE_MAX computes the maximum degrees of a digraph.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 October 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,J) is 1 if there is a direct link from node I to node J.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) INDEGREE_MAX, OUTDEGREE_MAX, the maximum indegree
!    and outdegree, considered independently, which may occur at different
!    nodes.
!
!    Output, integer ( kind = 4 ) DEGREE_MAX, the maximum value of the sum at each
!    node of the indegree and outdegree.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) degree
  integer ( kind = 4 ) degree_max
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indegree
  integer ( kind = 4 ) indegree_max
  integer ( kind = 4 ) j
  integer ( kind = 4 ) outdegree
  integer ( kind = 4 ) outdegree_max

  degree_max = 0
  indegree_max = 0
  outdegree_max = 0

  do i = 1, nnode

    indegree = 0
    outdegree = 0

    do j = 1, nnode
      outdegree = outdegree + adj(i,j)
      indegree = indegree + adj(j,i)
    end do

    degree = indegree + outdegree

    indegree_max = max ( indegree_max, indegree )
    outdegree_max = max ( outdegree_max, outdegree )
    degree_max = max ( degree_max, degree )

  end do

  return
end
subroutine dg_degree_seq ( adj, nnode, in_seq, out_seq )

!*****************************************************************************80
!
!! DG_DEGREE_SEQ computes the directed degree sequence.
!
!  Discussion:
!
!    The directed degree sequence is the sequence of indegrees
!    and the sequence of outdegrees, arranged to correspond to nodes of
!    successively decreasing total degree.  For nodes of equal degree, those
!    of higher outdegree take precedence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 October 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,J) is 1 if there is a direct link from node I to node J.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) IN_SEQ(NNODE), OUT_SEQ(NNODE),
!    the degree sequence.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) in_seq(nnode)
  integer ( kind = 4 ) out_seq(nnode)

  call dg_degree ( adj, nnode, in_seq, out_seq )

  call i4vec2_sort_d ( nnode, out_seq, in_seq )

  return
end
subroutine dg_edge_count ( adj, nnode, nedge )

!*****************************************************************************80
!
!! DG_EDGE_COUNT counts the number of edges in a digraph.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 October 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,J) is 1 if there is an edge from node I to node J.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) NEDGE, the number of edges.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nedge

  nedge = 0

  do i = 1, nnode
    do j = 1, nnode

      nedge = nedge + adj(i,j)

    end do
  end do

  return
end
subroutine dg_example_cycler ( adj, nnode )

!*****************************************************************************80
!
!! DG_EXAMPLE_CYCLER sets up the adjacency information for the cycler digraph.
!
!  Diagram:
!
!           A
!           V
!    9--><--7---<--3--><---4
!    |            /|      /
!    V           A |     /
!    |          /  |    /
!    5----<----1   V   A
!    |        /    |  /
!    V       A     | /
!    |      /      |/
!    2-->---8---<--6
!     \------>----/
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 October 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,J) is 1 if there is a direct link from node I to node J.
!
!    Output, integer ( kind = 4 ) NNODE, the number of nodes.
!
  implicit none

  integer ( kind = 4 ) adj(9,9)
  integer ( kind = 4 ) nnode

  nnode = 9

  adj(1:nnode,1:nnode) = 0

  adj(1,3) = 1
  adj(1,5) = 1

  adj(2,6) = 1
  adj(2,8) = 1

  adj(3,4) = 1
  adj(3,6) = 1
  adj(3,7) = 1

  adj(4,3) = 1

  adj(5,2) = 1

  adj(6,4) = 1
  adj(6,8) = 1

  adj(7,7) = 1
  adj(7,9) = 1

  adj(8,1) = 1

  adj(9,5) = 1
  adj(9,7) = 1

  return
end
subroutine dg_example_octo ( example, adj, nnode, seed )

!*****************************************************************************80
!
!! DG_EXAMPLE_OCTO sets up an 8 node example digraph.
!
!  Diagram:
!
!      1---2
!     /|   |\
!    8-+---+-3
!    | |   | |
!    7-+---+-4
!     \|   |/
!      6---5
!
!     Graph "A"
!
!    There are 12 digraphs to choose from, all on 8 nodes.
!
!    There are 8 underlying graphs.  The first 5 underlying graphs have
!    degree 3 at every node.  Graphs 6 and 7 have degree 5 at every node.
!    Graph 8 is disconnected.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) EXAMPLE, selects the example, and should be
!    between 1 and 13.
!
!    Output, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,J) is 1 if nodes I and J are adjacent and 0 otherwise.
!
!    IOutput, integer ( kind = 4 ) NNODE, the number of nodes, which should be 8.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
  implicit none


  integer ( kind = 4 ) adj(8,8)
  integer ( kind = 4 ) example
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nnode
  integer ( kind = 4 ) nsave
  integer ( kind = 4 ) seed

  if ( example <= 0 ) then
    nsave = i4_uniform ( 1, 13, seed )
  else
    example = mod ( example - 1, 13 ) + 1
    nsave = example
  end if

  nnode = 8

  adj(1:nnode,1:nnode) = 0

  do i = 1, nnode
    j = i + 1
    if ( nnode < j ) then
      j = j - nnode
    end if

    adj(i,j) = 1

  end do
!
!  Underlying graph 1.
!
  if ( nsave == 1 ) then

      adj(1,6) = 1
      adj(2,5) = 1
      adj(3,8) = 1
      adj(4,7) = 1

  else if ( nsave == 2 ) then

      adj(1,6) = 1
      adj(5,2) = 1
      adj(3,8) = 1
      adj(7,4) = 1
!
!  Underlying graph 2.
!  Digraphs 3 and 4 have different indegree/outdegree sequences.
!
  else if ( nsave == 3 ) then

    adj(1,6) = 1
    adj(6,1) = 1
    adj(2,8) = 1
    adj(8,2) = 1
    adj(3,5) = 1
    adj(5,3) = 1
    adj(4,7) = 1
    adj(7,4) = 1

  else if ( nsave == 4 ) then

    adj(1,6) = 1
    adj(2,8) = 1
    adj(3,5) = 1
    adj(4,7) = 1
!
!  Underlying graph 3
!  Digraphs 5 and 6 have the same indegree/outdegree sequences.
!
  else if ( nsave == 5 ) then

    adj(1,5) = 1
    adj(2,6) = 1
    adj(3,7) = 1
    adj(4,8) = 1

  else if ( nsave == 6 ) then

    adj(1:nnode,1:nnode) = 0

    adj(1,8) = 1
    adj(1,5) = 1
    adj(2,1) = 1
    adj(2,3) = 1
    adj(3,4) = 1
    adj(3,7) = 1
    adj(4,5) = 1
    adj(4,8) = 1
    adj(5,6) = 1
    adj(6,2) = 1
    adj(7,6) = 1
    adj(8,7) = 1
!
!  Underlying graph 4
!
  else if ( nsave == 7 ) then

    adj(3,1) = 1
    adj(4,2) = 1
    adj(5,7) = 1
    adj(6,8) = 1

  else if ( nsave == 8 ) then

    adj(3,1) = 1
    adj(4,2) = 1
    adj(5,7) = 1
    adj(8,6) = 1
!
!  Underlying graph 5
!
  else if ( nsave == 9 ) then

    adj(1,4) = 1
    adj(2,6) = 1
    adj(8,3) = 1

    adj(5,7) = 1
    adj(7,5) = 1

  else if ( nsave == 10 ) then

    adj(1,4) = 1
    adj(2,6) = 1
    adj(3,8) = 1

    adj(5,7) = 1
    adj(7,5) = 1
!
!  Underlying graph 6
!
  else if ( nsave == 11 ) then

    adj(1,4) = 1
    adj(1,5) = 1
    adj(1,6) = 1

    adj(2,5) = 1
    adj(2,6) = 1
    adj(2,7) = 1

    adj(3,6) = 1
    adj(3,7) = 1
    adj(3,8) = 1

    adj(4,7) = 1
    adj(4,8) = 1

    adj(5,8) = 1
!
!  Underlying graph 7
!
  else if ( nsave == 12 ) then

    adj(1,3) = 1
    adj(1,5) = 1
    adj(1,7) = 1

    adj(2,4) = 1
    adj(2,6) = 1
    adj(2,8) = 1

    adj(3,5) = 1
    adj(3,7) = 1

    adj(4,6) = 1
    adj(4,8) = 1

    adj(5,7) = 1

    adj(6,8) = 1
!
!  Underlying graph 8.
!
  else if ( nsave == 13 ) then

    adj(1,2) = 1
    adj(3,1) = 1
    adj(2,3) = 1
    adj(3,4) = 1
    adj(5,6) = 1
    adj(6,5) = 1
    adj(5,7) = 1
    adj(6,7) = 1

  end if
!
!  Now permute the graph.
!
  call i4mat_perm_random ( nnode, adj, seed )

  return
end
subroutine dg_example_sixty ( adj, nnode )

!*****************************************************************************80
!
!! DG_EXAMPLE_SIXTY sets up the adjacency information for the sixty digraph.
!
!  Discussion:
!
!    The nodes of the digraph are divisors of 60.  There is a link from I to
!    J if divisor I can be divided by divisor J.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 October 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,J) is 1 if there is a direct link from node I to node J.
!
!    Output, integer ( kind = 4 ) NNODE, the number of nodes.
!
  implicit none

  integer ( kind = 4 ) adj(12,12)
  integer ( kind = 4 ) d(12)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nnode

  nnode = 12

  d(1) = 60
  d(2) = 30
  d(3) = 20
  d(4) = 15
  d(5) = 12
  d(6) = 10
  d(7) = 6
  d(8) = 5
  d(9) = 4
  d(10) = 3
  d(11) = 2
  d(12) = 1

  do i = 1, nnode
    do j = 1, nnode
      if ( i == j ) then
        adj(i,j) = 0
      else if ( mod ( d(i), d(j) ) == 0 ) then
        adj(i,j) = 1
      else
        adj(i,j) = 0
      end if
    end do
  end do

  return
end
subroutine dg_order_code ( adj, nnode, npart, code, order )

!*****************************************************************************80
!
!! DG_ORDER_CODE returns the digraph code for a specific node ordering.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,J) is 1 if there is a direct link from node I to node J.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) NPART, the number of nodes to consider.
!    NPART should be between 1 and NNODE.
!
!    If NPART is NNODE, then the normal code is returned.
!
!    If NPART is less than NNODE, then the code is computed as
!    though only NPART nodes existed, namely, those specified in the
!    first NPART entries of order.  This option is provided so that
!    the routine can compute the portion of a code determined
!    by an incomplete ordering of the nodes.
!
!    Output, integer ( kind = 4 ) CODE(NNODE,NNODE), the code for this ordering.
!
!    Input, integer ( kind = 4 ) ORDER(NNODE), the ordering of the nodes.  ORDER(1)
!    is the first node, and so on.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) code(nnode,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ni
  integer ( kind = 4 ) nj
  integer ( kind = 4 ) npart
  integer ( kind = 4 ) order(nnode)

  do i = 1, nnode

    if ( i <= npart ) then

      ni = order(i)

      if ( order(i) < 1 .or. nnode < order(i) ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'DG_ORDER_CODE - Fatal error!'
        write ( *, '(a)' ) '  ORDER is not a proper permutation.'
        stop
      end if

    else
      ni = 0
    end if

    do j = 1, nnode

      if ( i == j ) then

        code(i,j) = 0

      else

        if ( j <= npart ) then

          nj = order(j)

          if ( order(j) < 1 .or. nnode < order(j) ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'DG_ORDER_CODE - Fatal error!'
            write ( *, '(a)' ) '  ORDER is not a proper permutation.'
            stop
          end if

        else
          nj = 0
        end if

        if ( ni == 0 .or. nj == 0 ) then
          code(i,j) = 0
        else
          code(i,j) = adj(ni,nj)
        end if

      end if

    end do
  end do

  return
end
subroutine dg_random ( adj, nnode, nedge, seed )

!*****************************************************************************80
!
!! DG_RANDOM generates a random digraph.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,J) is 1 if there is a direct link from node I to node J.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) NEDGE, the number of edges, which must be between
!    0 and NNODE*(NNODE-1).
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
  implicit none

  integer ( kind = 4 ) nnode
  integer ( kind = 4 ) nedge

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iwork(nedge)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) maxedge
  integer ( kind = 4 ) seed

  if ( nnode <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DG_RANDOM - Fatal error!'
    write ( *, '(a,i8)' ) '  NNODE = ', nedge
    write ( *, '(a)' ) '  but NNODE must be at least 1.'
    stop
  end if

  maxedge = nnode * ( nnode - 1 )

  if ( nedge < 0 .or. maxedge < nedge ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DG_RANDOM - Fatal error!'
    write ( *, '(a,i8)' ) '  NEDGE = ', nedge
    write ( *, '(a)' ) '  but NEDGE must be at least 0, and '
    write ( *, '(a,i8)' ) '  no more than ', maxedge
    stop
  end if

  adj(1:nnode,1:nnode) = 0
!
!  Pick a random NEDGE subset of NNODE*(NNODE-1).
!
  call ksub_random ( maxedge, nedge, iwork, seed )
!
!  The usable spots in the matrix are numbered as follows:
!
!   *    1    2   3  ...      n-2        n-1
!   n    *   n+1 n+2 ...     2n-1      2(n-1)
!  2n-1  2n   *  ... ... ........  ..........
!  .... ...  ... ... ...     *     (n-1)(n-1)
!  .... ...  ... ... ...   n(n-1)       *
!
  k = 0
  l = 1
  do i = 1, nnode
    do j = 1, nnode

      if ( i /= j ) then

        k = k + 1
        if ( l <= nedge ) then

          if ( k == iwork(l) ) then
            adj(i,j) = 1
            l = l + 1
          end if

        end if

      end if

    end do
  end do

  return
end
subroutine dmg_adj_max_max ( adj, nnode, adj_max_max )

!*****************************************************************************80
!
!! DMG_ADJ_MAX_MAX computes the adjacency maximum maximum of a dimultigraph.
!
!  Discussion:
!
!    The adjacency maximum maximum of a dimultigraph may be constructed by
!    computing the maximum entry of the adjacency matrix,
!
!    If two dimultigraphs are isomorphic, they must have the same
!    adjacency maximum maximum.
!
!  Example:
!
!    ADJ =
!       0 1 2 3
!       1 0 2 0
!       2 2 0 1
!       3 0 1 0
!
!    ADJ_MAX_MAX = 3
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,J) is the number of edges from node I to node J.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) ADJ_MAX_MAX, the adjacency maximum maximum.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) adj_max_max
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  adj_max_max = 0
  do i = 1, nnode
    do j = 1, nnode
      if ( i /= j ) then
        adj_max_max = max ( adj_max_max, adj(i,j) )
      end if
    end do
  end do

  return
end
subroutine dmg_adj_max_seq ( adj, nnode, adj_max_seq )

!*****************************************************************************80
!
!! DMG_ADJ_MAX_SEQ computes the adjacency maximum sequence of a dimultigraph.
!
!  Discussion:
!
!    The adjacency maximum sequence of a dimultigraph may be
!    constructed by computing the maximum entry of each row of the
!    adjacency matrix, and then sorting these values in descending order.
!
!  Example:
!
!    ADJ =
!       0 1 2 3
!       1 0 2 0
!       2 2 0 1
!       3 0 1 0
!
!    ADJ_MAX_SEQ =
!
!       3
!       3
!       2
!       2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,J) is the number of edges from node I to node J.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) ADJ_MAX_SEQ(NNODE), the adjacency maximum sequence.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) adj_max_seq(nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
!
!  Copy the adjacency matrix.
!
  do i = 1, nnode
    k = 0
    do j = 1, nnode
      if ( i /= j ) then
        k = max ( k, adj(i,j) )
      end if
    end do
    adj_max_seq(i) = k
  end do
!
!  Sort the elements.
!
  call i4vec_sort_heap_d ( nnode, adj_max_seq )

  return
end
subroutine dmg_adj_seq_u ( adj, nnode, adj_seq )

!*****************************************************************************80
!
!! DMG_ADJ_SEQ_U computes the unweighted adjacency sequence of a dimultigraph.
!
!  Discussion:
!
!    The unweighted adjacency sequence of a dimultigraph may be constructed
!    by replacing each nonzero entry by 1, sorting the entries of each row
!    in descending order, and then sorting the rows themselves in descending
!    order.
!
!  Example:
!
!    ADJ =
!       0 1 2 3
!       1 0 2 0
!       2 2 0 1
!       3 0 1 0
!
!    ADJ_SEQ =
!
!       1 1 1 0
!       1 1 1 0
!       1 1 0 0
!       1 1 0 0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,J) is the number of edges from node I to node J.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) ADJ_SEQ(NNODE,NNODE), the unweighted adjacency sequence.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) adj_seq(nnode,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
!
!  Copy the adjacency matrix.
!
  do i = 1, nnode
    do j = 1, nnode
      if ( adj(i,j) == 0 ) then
        adj_seq(i,j) = 0
      else
        adj_seq(i,j) = 1
      end if
    end do
  end do
!
!  Sort the elements of each row.
!
  call i4row_sort2_d ( nnode, nnode, adj_seq )
!
!  Sort the rows of the matrix.
!
  call i4row_sort_d ( nnode, nnode, adj_seq )

  return
end
subroutine dmg_adj_seq_w ( adj, nnode, adj_seq )

!*****************************************************************************80
!
!! DMG_ADJ_SEQ_W computes the weighted adjacency sequence of a dimultigraph.
!
!  Discussion:
!
!    The adjacency sequence of a dimultigraph may be constructed by sorting the
!    entries of each row of the adjacency matrix in descending order, and
!    then sorting the rows themselves in descending order.
!
!    If two dimultigraphs are isomorphic, they must have the same adjacency
!    sequence.
!
!  Example:
!
!    ADJ =
!       0 1 2 3
!       1 0 2 0
!       2 2 0 1
!       3 0 1 0
!
!    ADJ_SEQ =
!
!       3 2 1 0
!       3 1 0 0
!       2 2 1 0
!       2 1 0 0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,J) is the number of edges from node I to node J.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) ADJ_SEQ(NNODE,NNODE), the adjacency sequence.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) adj_seq(nnode,nnode)
!
!  Copy the adjacency matrix.
!
  adj_seq(1:nnode,1:nnode) = adj(1:nnode,1:nnode)
!
!  Sort the elements of each row.
!
  call i4row_sort2_d ( nnode, nnode, adj_seq )
!
!  Sort the rows of the matrix.
!
  call i4row_sort_d ( nnode, nnode, adj_seq )

  return
end
subroutine dmg_code_back ( adj, nnode, code, order )

!*****************************************************************************80
!
!! DMG_CODE_BACK computes a dimultigraph code via backtracking.
!
!  Discussion:
!
!    The code is the "largest" order code (in the lexicographic
!    sense) over all possible node orderings.  The backtracking method
!    organizes the search of all possible node orderings so that if
!    a partial node ordering is sure to generate an order code
!    less than the best so far, then all the orderings that begin with
!    this partial ordering are immediately dropped from consideration.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,J) is the number of edges from node I to node J.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) CODE(NNODE,NNODE), the code.
!
!    Output, integer ( kind = 4 ) ORDER(NNODE), the ordering of the nodes.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) bestcode(nnode,nnode)
  integer ( kind = 4 ) bestorder(nnode)
  integer ( kind = 4 ) code(nnode,nnode)
  logical, parameter :: debug = .false.
  integer ( kind = 4 ) index
  integer ( kind = 4 ) maxstack
  integer ( kind = 4 ) ncan(nnode)
  integer ( kind = 4 ) ncomp
  integer ( kind = 4 ) nopen
  integer ( kind = 4 ) nstack
  integer ( kind = 4 ) nswap
  integer ( kind = 4 ) order(nnode)
  integer ( kind = 4 ) result
  integer ( kind = 4 ) stack(4*nnode)

  if ( nnode <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DMG_CODE_BACK - Fatal error!'
    write ( *, '(a,i8)' ) '  NNODE = ', nnode
    stop
  end if

  maxstack = 4 * nnode
  nstack = 0
  stack(1) = 0

  ncomp = 0
  nswap = 0
!
!  Start with the identity ordering.
!
  call i4vec_indicator ( nnode, order )
!
!  Compute the corresponding code.
!
  call dmg_order_code ( adj, nnode, nnode, code, order )
!
!  Take this ordering and code as the best so far.
!
  bestorder(1:nnode) = order(1:nnode)
  bestcode(1:nnode,1:nnode) = code(1:nnode,1:nnode)
!
!  Now consider all possible orderings, and their codes.
!
  index = 0

  do

    call i4vec_backtrack ( nnode, order, index, nopen, nstack, stack, &
      maxstack, ncan )
!
!  If the backtrack routine has returned a complete candidate ordering, then
!  compute the resulting code, and see it it is better
!  then our current best.  Then go back for the next backtrack search.
!
    if ( index == 1 ) then

      call dmg_order_code ( adj, nnode, nnode, code, order )

      call dmg_code_compare ( bestcode, code, nnode, nnode, result )

      ncomp = ncomp + 1

      if ( result == - 1 ) then

        nswap = nswap + 1

        bestorder(1:nnode) = order(1:nnode)
        bestcode(1:nnode,1:nnode) = code(1:nnode,1:nnode)

      end if
!
!  If the backtrack routine has a partial reordering,
!  supply candidates for the next item in the ordering.
!
    else if ( index == 2 ) then

      call dmg_code_cand ( adj, bestcode, code, nnode, ncomp, nopen, &
        order, stack, maxstack, nstack, ncan )

    else

      exit

    end if

  end do
!
!  Once we have examined all possibilites, we are done.
!
!  Set the output ordering to the best ordering, and the output
!  code to the corresponding best code.
!
  order(1:nnode) = bestorder(1:nnode)
  code(1:nnode,1:nnode) = bestcode(1:nnode,1:nnode)

  if ( debug ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DMG_CODE_BACK:'
    write ( *, '(a,i8)' ) '  Comparisons: ', ncomp
    write ( *, '(a,i8)' ) '  Swaps:       ', nswap
  end if

  return
end
subroutine dmg_code_brute ( adj, nnode, code, order )

!*****************************************************************************80
!
!! DMG_CODE_BRUTE computes a dimultigraph code via brute force.
!
!  Discussion:
!
!    The code is the "largest" order code (in the lexicographic sense)
!    over all possible node orderings.  The brute force method considers
!    every node ordering, computes the corresponding order code, and
!    takes the largest one encountered.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,J) is the number of edges from node I to node J.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) CODE(NNODE,NNODE), the code.
!
!    Output, integer ( kind = 4 ) ORDER(NNODE), the ordering of the nodes.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) bestcode(nnode,nnode)
  integer ( kind = 4 ) bestorder(nnode)
  integer ( kind = 4 ) code(nnode,nnode)
  logical even
  logical more
  integer ( kind = 4 ) ncomp
  integer ( kind = 4 ) nswap
  integer ( kind = 4 ) order(nnode)
  integer ( kind = 4 ) result

  ncomp = 0
  nswap = 0
!
!  Start with the identity ordering.
!
  call i4vec_indicator ( nnode, order )
!
!  Compute the corresponding code.
!
  call dmg_order_code ( adj, nnode, nnode, code, order )
!
!  Take this ordering and code as the best so far.
!
  bestorder(1:nnode) = order(1:nnode)
  bestcode(1:nnode,1:nnode) = code(1:nnode,1:nnode)
!
!  Now consider all possible orderings, and their codes.
!
  more = .false.

  do

    call perm_next ( nnode, order, more, even )

    call dmg_order_code ( adj, nnode, nnode, code, order )

    call dmg_code_compare ( bestcode, code, nnode, nnode, result )

    ncomp = ncomp + 1

    if ( result == - 1 ) then

      nswap = nswap + 1

      bestorder(1:nnode) = order(1:nnode)
      bestcode(1:nnode,1:nnode) = code(1:nnode,1:nnode)

    end if

    if ( .not. more ) then
      exit
    end if

  end do
!
!  Once we have examined all possibilites, we are done.
!
!  Set the output ordering to the best ordering, and the output
!  code to the corresponding best code.
!
  order(1:nnode) = bestorder(1:nnode)
  code(1:nnode,1:nnode) = bestcode(1:nnode,1:nnode)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DMG_CODE_BRUTE:'
  write ( *, '(a,i8)' ) '  Comparisons: ', ncomp
  write ( *, '(a,i8)' ) '  Swaps:       ', nswap

  return
end
subroutine dmg_code_cand ( adj, bestcode, code, nnode, ncomp, nopen, &
  order, stack, maxstack, nstack, ncan )

!*****************************************************************************80
!
!! DMG_CODE_CAND finds candidates for a maximal dimultigraph code ordering.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,J) is the number of edges from node I to node J.
!
!    Input, integer ( kind = 4 ) BESTCODE(NNODE,NNODE), the best code so far.
!
!    Workspace, integer CODE(NNODE,NNODE).
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input/output, integer ( kind = 4 ) NCOMP, the number of code comparisons.
!    This routine updates NCOMP by 1 each time the routine is called.
!
!    Input, integer ( kind = 4 ) NOPEN, identifies the first open position in ORDER.
!
!    Input, integer ( kind = 4 ) ORDER(NNODE), contains in entries 1 through
!    NOPEN-1 the elements of the current partial list.
!
!    Input/output, integer ( kind = 4 ) STACK(MAXSTACK), used to store the new
!    candidates.
!
!    Input, integer ( kind = 4 ) MAXSTACK, the maximum size of the STACK array.
!
!    Input/output, integer ( kind = 4 ) NSTACK, the current length of the stack.
!    On output, NSTACK has been increased by the number of
!    candidates found.
!
!    Input/output, integer ( kind = 4 ) NCAN(NNODE), the number of candidates for
!    each position.
!
  implicit none

  integer ( kind = 4 ) maxstack
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) bestcode(nnode,nnode)
  integer ( kind = 4 ) code(nnode,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ifree(nnode)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) max_adj(nnode)
  integer ( kind = 4 ) max_adj_ni
  integer ( kind = 4 ) max_max_adj
  integer ( kind = 4 ) ncan(nnode)
  integer ( kind = 4 ) ncomp
  integer ( kind = 4 ) nfree
  integer ( kind = 4 ) ni
  integer ( kind = 4 ) nj
  integer ( kind = 4 ) nopen
  integer ( kind = 4 ) nstack
  integer ( kind = 4 ) order(nnode)
  integer ( kind = 4 ) result
  integer ( kind = 4 ) stack(maxstack)

  if ( nopen < 1 .or. nnode < nopen ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DMG_CODE_CAND - Fatal error 1!'
    write ( *, '(a)' ) '  1 <= NOPEN <= NNODE should be true, but'
    write ( *, '(a,i8)' ) '  NOPEN = ', nopen
    write ( *, '(a,i8)' ) '  NNODE = ', nnode
    stop
  end if

  ncan(nopen) = 0
!
!  If we have fixed at least one entry of the list, then compare
!  the code of the current incomplete ordering to the
!  code induced by the corresponding part of the current
!  best ordering.
!
!  If the current incomplete ordering is actually LESS than the
!  current best, then bail out with zero candidates.
!
  if ( 1 < nopen ) then

    call dmg_order_code ( adj, nnode, nopen-1, code, order )

    call dmg_code_compare ( bestcode, code, nnode, nopen-1, result )

    ncomp = ncomp + 1

    if ( result == + 1 ) then
      ncan(nopen) = 0
      return
    end if

  end if
!
!  Get a list of those nodes which have not been used yet.
!
  nfree = nnode + 1 - nopen
  call perm_free ( order, nopen-1, ifree, nfree )
!
!  To find candidates, we consider all the ordered nodes.
!  We find the lowest ordered node that has unordered neighbors.
!  We take the unordered neighbors with maximal adjacency.
!
  ncan(nopen) = 0
!
!  For each ordered node NI...
!
  do i = 1, nopen-1

    ni = order(i)
!
!  ...note the maximum adjacency FROM NI to any unordered node NJ...
!
    max_adj_ni = 0
    do j = 1, nfree
      nj = ifree(j)
      max_adj_ni = max ( max_adj_ni, adj(ni,nj) )
    end do
!
!   ...and take as candidates all unordered nodes NJ with maximal
!   adjacency FROM NI.
!
    if ( 0 < max_adj_ni ) then

      do j = 1, nfree

        nj = ifree(j)

        if ( adj(ni,nj) == max_adj_ni ) then
          ncan(nopen) = ncan(nopen) + 1
          nstack = nstack + 1

          if ( maxstack < nstack ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'DMG_CODE_CAND - Fatal error 2!'
            write ( *, '(a)' ) '  MAXSTACK < NSTACK!'
            write ( *, '(a,i8)' ) '  NSTACK =   ', nstack
            write ( *, '(a,i8)' ) '  MAXSTACK = ', maxstack
            stop
          end if

          stack(nstack) = nj
        end if

      end do

      return

    end if
!
!  Else, note the maximum adjacency TO NI from any unordered node NJ...
!
    max_adj_ni = 0
    do j = 1, nfree
      nj = ifree(j)
      max_adj_ni = max ( max_adj_ni, adj(nj,ni) )
    end do
!
!   ...and take as candidates all unordered nodes NJ with maximal
!   adjacency TO NI.
!
    if ( 0 < max_adj_ni ) then

      do j = 1, nfree

        nj = ifree(j)

        if ( adj(nj,ni) == max_adj_ni ) then
          ncan(nopen) = ncan(nopen) + 1
          nstack = nstack + 1

          if ( maxstack < nstack ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' )'DMG_CODE_CAND - Fatal error 2!'
            write ( *, '(a)' ) '  MAXSTACK < NSTACK!'
            write ( *, '(a,i8)' ) '  NSTACK =   ', nstack
            write ( *, '(a,i8)' ) '  MAXSTACK = ', maxstack
            stop
          end if

          stack(nstack) = nj

        end if

      end do

      return

    end if

  end do
!
!  If we get here, no unordered nodes are connected in any way to
!  ordered nodes.  This can happen in two ways:
!
!  * NOPEN = 1;
!  * The dimultigraph is disconnected;
!
!  For each free node, compute the maximum adjacency TO any other free node.
!  Take the maximum of this value over all free nodes.
!  Candidates are free nodes with this maximum value.
!
  max_max_adj = 0

  do i = 1, nfree

    ni = ifree(i)

    max_adj(i) = 0
    do j = 1, nfree
      nj = ifree(j)
      if ( ni /= nj ) then
        max_adj(i) = max ( max_adj(i), adj(ni,nj) )
      end if
    end do

    max_max_adj = max ( max_max_adj, max_adj(i) )

  end do

  ncan(nopen) = 0

  do i = 1, nfree

    if ( max_adj(i) == max_max_adj ) then

      ncan(nopen) = ncan(nopen) + 1
      nstack = nstack + 1

      if ( maxstack < nstack ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'DMG_CODE_CAND - Fatal error 2!'
        write ( *, '(a)' ) '  MAXSTACK < NSTACK!'
        write ( *, '(a,i8)' ) '  NSTACK =   ', nstack
        write ( *, '(a,i8)' ) '  MAXSTACK = ', maxstack
        stop
      end if

      stack(nstack) = ifree(i)

    end if
  end do

  return
end
subroutine dmg_code_compare ( code1, code2, nnode, npart, result )

!*****************************************************************************80
!
!! DMG_CODE_COMPARE compares two partial dimultigraph codes.
!
!  Discussion:
!
!    CODE1 = CODE2 if every digit of both codes is equal.
!
!    Otherwise, traverse the entries in a funny diagonal way, suggested
!    by this diagram for the first 20 entries:
!
!       *  1  3  7 13
!       2  *  5  9 15
!       4  6  * 11 17
!       8 10 12  * 19
!      14 16 18 20  *
!
!    As we do that, we examine the corresponding digits of the two codes.
!    For the first entry, (I,J), where the two codes differ, we say:
!
!      if ( CODE1(I,J) < CODE2(I,J) ) then we say
!        CODE1 < CODE2
!      else
!        CODE2 < CODE1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) CODE1(NNODE,NNODE), CODE2(NNODE,NNODE),
!    codes to compare.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) NPART, specifies the portion of the codes
!    to compare.  NPART should be between 1 and NNODE.
!
!    If NPART = NNODE, then the full codes are compared.
!
!    If NPART < NNODE, then only entries corresponding to I and J
!    both less than or equal to NPART will be compared.
!
!    Output, integer ( kind = 4 ) RESULT:
!    -1, CODE1 < CODE2;
!     0, CODE1 = CODE2;
!    +1, CODE2 < CODE1.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) code1(nnode,nnode)
  integer ( kind = 4 ) code2(nnode,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) npart
  integer ( kind = 4 ) result

  do j = 2, npart

    do i = 1, j-1

      if ( code1(i,j) < code2(i,j) ) then

        result = - 1
        return

      else if ( code2(i,j) < code1(i,j) ) then

        result = + 1
        return

      end if

      if ( code1(j,i) < code2(j,i) ) then

        result = - 1
        return

      else if ( code2(j,i) < code1(j,i) ) then

        result = + 1
        return

      end if

    end do

  end do

  result = 0

  return
end
subroutine dmg_code_print ( nnode, code, title )

!*****************************************************************************80
!
!! DMG_CODE_PRINT prints out a dimultigraph code.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) CODE(NNODE,NNODE), the code.
!
!    Input, character ( len = * ) TITLE, a title to be printed.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) code(nnode,nnode)
  integer ( kind = 4 ) i
  character ( len = 80 ) string
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  do i = 1, nnode

    write ( string, '(80i1)' ) code(i,1:nnode)
    string(i:i) = '.'
    write ( *, '(2x,a)' ) string(1:nnode)

  end do

  return
end
subroutine dmg_compare ( adj1, nnode1, adj2, nnode2, order1, &
  order2, result )

!*****************************************************************************80
!
!! DMG_COMPARE determines if dimultigraphs G1 and G2 are isomorphic.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ1(NNODE1,NNODE1), the adjacency information
!    for G1.
!
!    Input, integer ( kind = 4 ) NNODE1, the number of nodes in G1.
!
!    Input, integer ( kind = 4 ) ADJ2(NNODE2,NNODE2), the adjacency information
!    for G2.
!
!    Input, integer ( kind = 4 ) NNODE2, the number of nodes in G2.
!
!    Output, integer ( kind = 4 ) RESULT, is 0 if the dimultigraphs are isomorphic,
!    -I if G1 < G2 for test #I, and
!    +I if G2 < G1 for test #I.
!
!    Output, integer ( kind = 4 ) ORDER1(NNODE1), ORDER2(NNODE2).
!    If RESULT = 0, then ORDER1 and ORDER2 are reorderings of the nodes of
!    G1 and G2 which exhibit the isomorphism.
!
  implicit none

  integer ( kind = 4 ) nnode1
  integer ( kind = 4 ) nnode2

  integer ( kind = 4 ) adj_max_max_1
  integer ( kind = 4 ) adj_max_max_2
  integer ( kind = 4 ) adj1(nnode1,nnode1)
  integer ( kind = 4 ) adj2(nnode2,nnode2)
  integer ( kind = 4 ) code1(nnode1,nnode1)
  integer ( kind = 4 ) code2(nnode2,nnode2)
  integer ( kind = 4 ) in_seq1(nnode1)
  integer ( kind = 4 ) in_seq2(nnode2)
  integer ( kind = 4 ) nedge_u_1
  integer ( kind = 4 ) nedge_u_2
  integer ( kind = 4 ) nedge_w_1
  integer ( kind = 4 ) nedge_w_2
  integer ( kind = 4 ) order1(nnode1)
  integer ( kind = 4 ) order2(nnode2)
  integer ( kind = 4 ) out_seq1(nnode1)
  integer ( kind = 4 ) out_seq2(nnode2)
  integer ( kind = 4 ) result
  integer ( kind = 4 ) seq1(nnode1)
  integer ( kind = 4 ) seq2(nnode2)
!
!  Test 1: Count the nodes.
!
  if ( nnode1 < nnode2 ) then
    result = - 1
    return
  else if ( nnode2 < nnode1 ) then
    result = + 1
    return
  end if
!
!  Test 2: Compare the unweighted edges.
!
  call dmg_edge_count ( adj1, nnode1, nedge_u_1, nedge_w_1 )

  call dmg_edge_count ( adj2, nnode2, nedge_u_2, nedge_w_2 )

  if ( nedge_u_1 < nedge_u_2 ) then
    result = - 2
    return
  else if ( nedge_u_2 < nedge_u_1 ) then
    result = + 2
    return
  end if
!
!  Test 3: Compare the weighted edges.
!
  if ( nedge_w_1 < nedge_w_2 ) then
    result = - 3
    return
  else if ( nedge_w_2 < nedge_w_1 ) then
    result = + 3
    return
  end if
!
!  Test 4: Compare the unweighted outdegree sequences.
!
  call dmg_degree_seq_u ( adj1, nnode1, in_seq1, out_seq1 )

  call dmg_degree_seq_u ( adj2, nnode2, in_seq2, out_seq2 )

  call i4vec_compare ( nnode1, out_seq1, out_seq2, result )

  if ( result < 0 ) then
    result = - 4
    return
  else if ( 0 < result ) then
    result = + 4
    return
  end if
!
!  Test 5: Compare the unweighted indegree sequences.
!
  call i4vec_compare ( nnode1, in_seq1, in_seq2, result )

  if ( result < 0 ) then
    result = - 5
    return
  else if ( 0 < result ) then
    result = + 5
    return
  end if
!
!  Test 6: Compare the weighted outdegree sequences.
!
  call dmg_degree_seq_w ( adj1, nnode1, in_seq1, out_seq1 )

  call dmg_degree_seq_w ( adj2, nnode2, in_seq2, out_seq2 )

  call i4vec_compare ( nnode1, out_seq1, out_seq2, result )

  if ( result < 0 ) then
    result = - 6
    return
  else if ( 0 < result ) then
    result = + 6
    return
  end if
!
!  Test 7: Compare the weighted indegree sequences.
!
  call i4vec_compare ( nnode1, in_seq1, in_seq2, result )

  if ( result < 0 ) then
    result = - 7
    return
  else if ( 0 < result ) then
    result = + 7
    return
  end if
!
!  Test 8: Compare the adjacency max max.
!
  call dmg_adj_max_max ( adj1, nnode1, adj_max_max_1 )

  call dmg_adj_max_max ( adj2, nnode2, adj_max_max_2 )

  if ( adj_max_max_1 < adj_max_max_2 ) then
    result = - 8
    return
  else if ( adj_max_max_2 < adj_max_max_1 ) then
    result = + 8
    return
  end if
!
!  Test 9: Compare the adjacency max sequences.
!
  call dmg_adj_max_seq ( adj1, nnode1, seq1 )

  call dmg_adj_max_seq ( adj2, nnode2, seq2 )

  call i4vec_compare ( nnode1, seq1, seq2, result )

  if ( result < 0 ) then
    result = - 9
    return
  else if ( 0 < result ) then
    result = + 9
    return
  end if
!
!  Test 10: Compare the unweighted adjacency sequences.
!
  call dmg_adj_seq_u ( adj1, nnode1, code1 )

  call dmg_adj_seq_u ( adj2, nnode2, code2 )

  call i4mat_row_compare ( nnode1, nnode1, code1, code2, result )

  if ( result < 0 ) then
    result = - 10
    return
  else if ( 0 < result ) then
    result = + 10
    return
  end if
!
!  Test 11: Compare the weighted adjacency sequences.
!
  call dmg_adj_seq_w ( adj1, nnode1, code1 )

  call dmg_adj_seq_w ( adj2, nnode2, code2 )

  call i4mat_row_compare ( nnode1, nnode1, code1, code2, result )

  if ( result < 0 ) then
    result = - 11
    return
  else if ( 0 < result ) then
    result = + 11
    return
  end if
!
!  Test 12: Compare the codes.
!
  call dmg_code_back ( adj1, nnode1, code1, order1 )

  call dmg_code_back ( adj2, nnode2, code2, order2 )

  call dmg_code_compare ( code1, code2, nnode1, nnode1, result )

  if ( result < 0 ) then
    result = - 12
    return
  else if ( 0 < result ) then
    result = + 12
    return
  end if

  result = 0

  return
end
subroutine dmg_degree_seq_u ( adj, nnode, in_seq, out_seq )

!*****************************************************************************80
!
!! DMG_DEGREE_SEQ_U: the unweighted directed degree sequence of a dimultigraph.
!
!  Discussion:
!
!    The unweighted directed degree sequence is the sequence of indegrees
!    and the sequence of outdegrees, ignoring edge multiplicity, arranged
!    to correspond to nodes of successively decreasing total degree.  For
!    nodes of equal degree, those of higher outdegree take precedence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,J) is 1 if there is a direct link from node I to node J.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) IN_SEQ(NNODE), OUT_SEQ(NNODE),
!    the unweighted directed degree sequences.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) in_seq(nnode)
  integer ( kind = 4 ) out_seq(nnode)

  call dmg_degree_u ( adj, nnode, in_seq, out_seq )

  call i4vec2_sort_d ( nnode, out_seq, in_seq )

  return
end
subroutine dmg_degree_seq_w ( adj, nnode, in_seq, out_seq )

!*****************************************************************************80
!
!! DMG_DEGREE_SEQ_W: weighted directed degree sequence of a dimultigraph.
!
!  Discussion:
!
!    The weighted directed degree sequence is the sequence of indegrees
!    and the sequence of outdegrees, with edge multiplicity, arranged
!    to correspond to nodes of successively decreasing total degree.  For
!    nodes of equal degree, those of higher outdegree take precedence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,J) is 1 if there is a direct link from node I to node J.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) IN_SEQ(NNODE), OUT_SEQ(NNODE),
!    the weighted directed degree sequences.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) in_seq(nnode)
  integer ( kind = 4 ) out_seq(nnode)

  call dmg_degree_w ( adj, nnode, in_seq, out_seq )

  call i4vec2_sort_d ( nnode, out_seq, in_seq )

  return
end
subroutine dmg_degree_u ( adj, nnode, indegree, outdegree )

!*****************************************************************************80
!
!! DMG_DEGREE_U computes the unweighted degrees of a dimultigraph.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,J) is 1 if there is a direct link from node I to node J.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) INDEGREE(NNODE), OUTDEGREE(NNODE),
!    the unweighted indegree and outdegree of the nodes.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indegree(nnode)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) outdegree(nnode)

  indegree(1:nnode) = 0
  outdegree(1:nnode) = 0

  do i = 1, nnode
    do j = 1, nnode
      if ( adj(i,j) /= 0 ) then
        outdegree(i) = outdegree(i) + 1
        indegree(j) = indegree(j) + 1
      end if
    end do
  end do

  return
end
subroutine dmg_degree_w ( adj, nnode, indegree, outdegree )

!*****************************************************************************80
!
!! DMG_DEGREE_W computes the weighted degrees of a dimultigraph.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,J) is 1 if there is a direct link from node I to node J.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) INDEGREE(NNODE), OUTDEGREE(NNODE),
!    the weighted indegree and outdegree of the nodes.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indegree(nnode)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) outdegree(nnode)

  indegree(1:nnode) = 0
  outdegree(1:nnode) = 0

  do i = 1, nnode
    do j = 1, nnode
      if ( adj(i,j) /= 0 ) then
        outdegree(i) = outdegree(i) + adj(i,j)
        indegree(j) = indegree(j) + adj(i,j)
      end if
    end do
  end do

  return
end
subroutine dmg_edge_count ( adj, nnode, nedge_u, nedge_w )

!*****************************************************************************80
!
!! DMG_EDGE_COUNT counts the number of edges in a dimultigraph.
!
!  Discussion:
!
!    The number of "unweighted" edges is the number of edges in the
!    underlying digraph, or the number of edges that would be counted
!    if each set of multiple edges was replaced by a single edge.
!
!    The number of "weighted" edges counts separately each edge of a
!    multiple edge.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,J) is the number of edges from node I to node J.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) NEDGE_U, the number of unweighted edges.
!
!    Output, integer ( kind = 4 ) NEDGE_W, the number of weighted edges.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nedge_u
  integer ( kind = 4 ) nedge_w

  nedge_u = 0
  nedge_w = 0

  do i = 1, nnode
    do j = 1, nnode

      if ( i /= j ) then
        nedge_w = nedge_w + adj(i,j)
        if ( 0 < adj(i,j) ) then
          nedge_u = nedge_u + 1
        end if
      end if

    end do
  end do

  return
end
subroutine dmg_example_octo ( example, adj, nnode, seed )

!*****************************************************************************80
!
!! DMG_EXAMPLE_OCTO sets up an 8 node example dimultigraph.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) EXAMPLE, selects the example, and should be between
!    1 and 8.
!
!    Output, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,J) is the number of edges from node I to node J.
!
!    Output, integer ( kind = 4 ) NNODE, the number of nodes, which should be 8.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(8,8)
  integer ( kind = 4 ) example
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) nsave
  integer ( kind = 4 ) seed

  if ( example <= 0 ) then
    nsave = i4_uniform ( 1, 8, seed )
  else
    nsave = mod ( example - 1, 8 ) + 1
  end if

  nnode = 8

  adj(1:nnode,1:nnode) = 0
!
!  The "basic" graph.
!
  if ( nsave == 1 ) then

    adj(1,2) = 1
    adj(1,6) = 1
    adj(1,8) = 1
    adj(2,3) = 1
    adj(2,5) = 2
    adj(3,4) = 1
    adj(3,8) = 3
    adj(4,5) = 1
    adj(4,7) = 4
    adj(5,6) = 1
    adj(6,7) = 1
    adj(7,8) = 1
!
!  Same NNODE, different number of unweighted edges.
!
  else if ( nsave == 2 ) then

    adj(1,2) = 1
    adj(1,6) = 1
    adj(1,8) = 1
    adj(2,3) = 1
    adj(2,8) = 2
    adj(3,4) = 1
    adj(3,5) = 3
    adj(4,5) = 1
    adj(4,7) = 3
    adj(5,6) = 1
    adj(6,7) = 1
    adj(7,8) = 1
!
!  Same NNODE, NEDGE_U, different number of weighted edges.
!
  else if ( nsave == 3 ) then

    adj(1,2) = 1
    adj(1,6) = 1
    adj(1,8) = 2
    adj(2,3) = 1
    adj(2,5) = 2
    adj(3,4) = 1
    adj(3,8) = 3
    adj(4,5) = 1
    adj(4,7) = 4
    adj(5,6) = 1
    adj(6,7) = 1
    adj(7,8) = 1
!
!  Same NNODE, NEDGE_U, NEDGE_W, different degree sequence.
!
  else if ( nsave == 4 ) then

    adj(1,2) = 1
    adj(1,5) = 2
    adj(1,8) = 1
    adj(2,3) = 1
    adj(2,6) = 2
    adj(3,4) = 1
    adj(3,7) = 3
    adj(4,5) = 1
    adj(4,8) = 3
    adj(5,6) = 1
    adj(6,7) = 1
    adj(7,8) = 1
!
!  Same NNODE, NEDGE_U, NEDGE_W, degree sequence, different ADJ_MAX_MAX.
!
  else if ( nsave == 5 ) then

    adj(1,2) = 1
    adj(1,7) = 1
    adj(1,8) = 1
    adj(2,3) = 1
    adj(2,5) = 1
    adj(2,8) = 1
    adj(3,4) = 1
    adj(3,7) = 1
    adj(3,8) = 2
    adj(4,5) = 2
    adj(4,6) = 1
    adj(4,7) = 2
    adj(5,6) = 1
    adj(6,7) = 1
    adj(7,8) = 1
!
!  Same NNODE, NEDGE_U, NEDGE_W, degree sequence, ADJ_MAX_MAX, different ADJ_MAX_SEQ.
!
  else if ( nsave == 6 ) then

    adj(1,2) = 1
    adj(1,6) = 1
    adj(1,8) = 1
    adj(2,3) = 1
    adj(2,5) = 2
    adj(3,4) = 2
    adj(3,8) = 2
    adj(3,4) = 2
    adj(4,7) = 4
    adj(5,6) = 1
    adj(5,8) = 1
    adj(6,7) = 1
    adj(7,8) = 1
!
!  Same NNODE, NEDGE_U, NEDGE_W, degree sequence, ADJ_MAX_MAX, ADJ_MAX_SEQ,
!  different ADJ_SEQ.
!
  else if ( nsave == 7 ) then

    adj(1,2) = 4
    adj(1,3) = 2
    adj(2,4) = 2
    adj(3,4) = 3
    adj(5,6) = 2
    adj(5,7) = 1
    adj(5,8) = 1
    adj(6,7) = 1
    adj(6,8) = 1
    adj(7,8) = 1
!
!  Same NNODE, NEDGE_U, NEDGE_W, degree sequence, ADJ_MAX_MAX, ADJ_MAX_SEQ,
!  ADJ_SEQ, different code.
!
  else if ( nsave == 8 ) then

    adj(1,2) = 1
    adj(1,4) = 1
    adj(1,6) = 1
    adj(2,3) = 1
    adj(2,5) = 2
    adj(3,4) = 1
    adj(3,8) = 3
    adj(4,7) = 4
    adj(5,6) = 1
    adj(5,8) = 1
    adj(6,7) = 1
    adj(7,8) = 1

  end if
!
!  Now permute the nodes.
!
  call i4mat_perm_random ( nnode, adj, seed )

  return
end
subroutine dmg_order_code ( adj, nnode, npart, code, order )

!*****************************************************************************80
!
!! DMG_ORDER_CODE returns the dimultigraph code for a specific node ordering.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,J) is the number of edges from node I to node J.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) NPART, the number of nodes to consider.
!    NPART should be between 1 and NNODE.
!
!    If NPART is NNODE, then the full code is returned.
!
!    If NPART is less than NNODE, then the code is computed as
!    though only NPART nodes existed, namely, those specified in the
!    first NPART entries of order.  This option is provided so that
!    the routine can compute the portion of a code determined
!    by an incomplete ordering of the nodes.
!
!    Output, integer ( kind = 4 ) CODE(NNODE,NNODE), the code for this ordering.
!
!    Input, integer ( kind = 4 ) ORDER(NNODE), the ordering of the nodes.  ORDER(1)
!    is the first node, and so on.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) code(nnode,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ni
  integer ( kind = 4 ) nj
  integer ( kind = 4 ) npart
  integer ( kind = 4 ) order(nnode)

  do i = 1, nnode

    if ( i <= npart ) then

      ni = order(i)

      if ( order(i) < 1 .or. nnode < order(i) ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'DMG_ORDER_CODE - Fatal error!'
        write ( *, '(a)' ) '  ORDER is not a proper permutation.'
        stop
      end if

    else
      ni = 0
    end if

    do j = 1, nnode

      if ( j <= npart ) then

        nj = order(j)

        if ( order(j) < 1 .or. nnode < order(j) ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'DMG_ORDER_CODE - Fatal error!'
          write ( *, '(a)' ) '  ORDER is not a proper permutation.'
          stop
        end if

      else
        nj = 0
      end if

      if ( ni == nj ) then

        code(i,j) = 0

      else if ( ni == 0 .or. nj == 0 ) then

        code(i,j) = 0

      else

        code(i,j) = adj(ni,nj)

      end if

    end do
  end do

  return
end
subroutine dmg_print ( adj, nnode, title )

!*****************************************************************************80
!
!! DMG_PRINT prints out an adjacency matrix for a dimultigraph.
!
!  Discussion:
!
!    Values between 0 and 9 will be printed as is.  Other values will
!    be printed as '*'.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency matrix.
!    ADJ(I,J) is the number of edges from node I to node J.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, character ( len = * ) TITLE, a title to be printed.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  character ( len = 80 ) string
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  do i = 1, nnode

    jhi = min ( nnode, 80 )

    do j = 1, jhi

      if ( 0 <= adj(i,j) .and. adj(i,j) <= 9 ) then
        string(j:j) = char ( 48 + adj(i,j) )
      else
        string(j:j) = '*'
      end if

    end do

    write ( *, '(2x,a)' ) string(1:jhi)

  end do

  return
end
subroutine dmg_random ( adj, nnode, nedge, seed )

!*****************************************************************************80
!
!! DMG_RANDOM generates a random dimultigraph on NNODE nodes with NEDGE edges.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency matrix.
!    ADJ(I,J) is the number of edges from node I to node J.
!    ADJ(I,I) will always be 0.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) NEDGE, the number of edges.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
  implicit none

  integer ( kind = 4 ) nnode
  integer ( kind = 4 ) nedge

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
!
!  Initialize the adjacency matrix.
!
  adj(1:nnode,1:nnode) = 0
!
!  Essentially, flip a coin NEDGE times to decide where each edge goes.
!
  do k = 1, nedge

    i = i4_uniform ( 1, nnode,   seed )
    j = i4_uniform ( 1, nnode-1, seed )

    if ( i <= j ) then
      j = j + 1
    end if

    adj(i,j) = adj(i,j) + 1

  end do

  return
end
subroutine g_arc_node_count ( nedge, inode, jnode, mnode, nnode )

!*****************************************************************************80
!
!! G_ARC_NODE_COUNT counts the number of nodes in a graph.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 October 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEDGE, the number of edges.
!
!    Input, integer ( kind = 4 ) INODE(NEDGE), JNODE(NEDGE).  INODE(I) and JNODE(I)
!    are the start and end nodes of the I-th edge.
!
!    Output, integer ( kind = 4 ) MNODE, the maximum node index.
!
!    Output, integer ( kind = 4 ) NNODE, the number of distinct nodes.
!
  implicit none

  integer ( kind = 4 ) nedge

  integer ( kind = 4 ) iedge
  integer ( kind = 4 ) inode(nedge)
  integer ( kind = 4 ) jnode(nedge)
  integer ( kind = 4 ) knode(2*nedge)
  integer ( kind = 4 ) mnode
  integer ( kind = 4 ) nnode

  mnode = max ( maxval ( inode(1:nedge) ), maxval ( jnode(1:nedge) ) )
!
!  Copy all the node labels into KNODE,
!  sort KNODE,
!  count the unique entries.
!
!  That's NNODE.
!
  knode(1:nedge) = inode(1:nedge)

  do iedge = 1, nedge
    knode(nedge+iedge) = jnode(iedge)
  end do

  call i4vec_sort_heap_a ( 2*nedge, knode )

  call i4vec_sorted_unique_count ( 2*nedge, knode, nnode )

  return
end
subroutine g_arc_to_g_adj ( nedge, inode, jnode, adj, nnode )

!*****************************************************************************80
!
!! G_ARC_TO_G_ADJ converts an arc list graph to an adjacency graph.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 October 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEDGE, the number of edges.
!
!    Input, integer ( kind = 4 ) INODE(NEDGE), JNODE(NEDGE), the edge array for an
!    undirected graph.  The I-th edge connects nodes INODE(I) and JNODE(I).
!
!    Output, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!
!    Output, integer ( kind = 4 ) NNODE, the number of nodes.
!
  implicit none

  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) inode(nedge)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jnode(nedge)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) mnode
!
!  Determine the number of nodes.
!
  call g_arc_node_count ( nedge, inode, jnode, mnode, nnode )

  adj(1:nnode,1:nnode) = 0

  do k = 1, nedge
    i = inode(k)
    j = jnode(k)
    adj(i,j) = 1
    adj(j,i) = 1
  end do

  return
end
subroutine g_code_back ( adj, nnode, code, order )

!*****************************************************************************80
!
!! G_CODE_BACK computes a graph code via backtracking.
!
!  Discussion:
!
!    The code is the "largest" order code (in the lexicographic
!    sense) over all possible node orderings.  The backtracking method
!    organizes the search of all possible node orderings so that if
!    a partial node ordering is sure to generate an order code
!    less than the best so far, then all the orderings that begin with
!    this partial ordering are immediately dropped from consideration.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,J) is 1 if nodes I and J are adjacent.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) CODE(NNODE,NNODE), the code.
!
!    Output, integer ( kind = 4 ) ORDER(NNODE), the ordering of the nodes.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) bestcode(nnode,nnode)
  integer ( kind = 4 ) bestorder(nnode)
  integer ( kind = 4 ) code(nnode,nnode)
  logical, parameter :: debug = .false.
  integer ( kind = 4 ) index
  integer ( kind = 4 ) maxstack
  integer ( kind = 4 ) ncan(nnode)
  integer ( kind = 4 ) ncomp
  integer ( kind = 4 ) nopen
  integer ( kind = 4 ) npart
  integer ( kind = 4 ) nstack
  integer ( kind = 4 ) nswap
  integer ( kind = 4 ) order(nnode)
  integer ( kind = 4 ) result
  integer ( kind = 4 ) stack(4*nnode)

  if ( nnode <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'G_CODE_BACK - Fatal error!'
    write ( *, '(a,i8)' ) '  NNODE = ', nnode
    stop
  end if

  if ( debug ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'G_CODE_BACK - DEBUG - Entered routine.'
  end if

  maxstack = 4 * nnode
  nstack = 0
  stack(1) = 0

  ncomp = 0
  nswap = 0
!
!  Start with the identity ordering.
!
  call i4vec_indicator ( nnode, order )
!
!  Compute the corresponding code.
!
  npart = nnode
  call g_order_code ( adj, nnode, npart, code, order )
!
!  Take this ordering and code as the best so far.
!
  bestorder(1:nnode) = order(1:nnode)
  bestcode(1:nnode,1:nnode) = code(1:nnode,1:nnode)

  if ( debug ) then
    write ( *, '(a)' ) 'G_CODE_BACK - DEBUG - Begin backtrack search.'
  end if
!
!  Now consider all possible orderings, and their codes.
!
  index = 0

  do

    call i4vec_backtrack ( nnode, order, index, nopen, nstack, stack, &
      maxstack, ncan )
!
!  If the backtrack routine has returned a complete candidate ordering, then
!  compute the resulting code, and see it it is better
!  then our current best.  Then go back for the next backtrack search.
!
    if ( index == 1 ) then

      call g_order_code ( adj, nnode, npart, code, order )

      call g_code_compare ( bestcode, code, nnode, npart, result )

      ncomp = ncomp + 1

      if ( result == - 1 ) then

        nswap = nswap + 1

        bestorder(1:nnode) = order(1:nnode)
        bestcode(1:nnode,1:nnode) = code(1:nnode,1:nnode)

      end if
!
!  If the backtrack routine has a partial reordering,
!  supply candidates for the next item in the ordering.
!
    else if ( index == 2 ) then

      call g_code_cand ( adj, bestcode, code, nnode, ncomp, nopen, order, &
        stack, maxstack, nstack, ncan )
!
!  If we have examined all possibilities, we are done.
!
    else

      exit

    end if

  end do
!
!  Set the output ordering to the best ordering, and the output
!  code to the corresponding best code.
!
  order(1:nnode) = bestorder(1:nnode)
  code(1:nnode,1:nnode) = bestcode(1:nnode,1:nnode)

  if ( debug ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'G_CODE_BACK:'
    write ( *, '(a,i8)' ) '  Comparisons: ', ncomp
    write ( *, '(a,i8)' ) '  Swaps:       ', nswap
  end if

  return
end
subroutine g_code_brute ( adj, nnode, code, order )

!*****************************************************************************80
!
!! G_CODE_BRUTE computes a graph code via brute force.
!
!  Discussion:
!
!    The code is the "largest" order code in the lexicographic
!    sense over all node orderings.  The brute force method
!    considers every node ordering, computes the corresponding
!    order code, and takes the largest one encountered.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 May 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,J) is 1 if nodes I and J are adjacent.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) CODE(NNODE,NNODE), the code.
!
!    Output, integer ( kind = 4 ) ORDER(NNODE), the ordering of the nodes.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) bestcode(nnode,nnode)
  integer ( kind = 4 ) bestorder(nnode)
  integer ( kind = 4 ) code(nnode,nnode)
  logical even
  logical more
  integer ( kind = 4 ) ncomp
  integer ( kind = 4 ) nswap
  integer ( kind = 4 ) order(nnode)
  integer ( kind = 4 ) result

  ncomp = 0
  nswap = 0
!
!  Start with the identity ordering.
!
  call i4vec_indicator ( nnode, order )
!
!  Compute the corresponding code.
!
  call g_order_code ( adj, nnode, nnode, code, order )
!
!  Take this ordering and code as the best so far.
!
  bestorder(1:nnode) = order(1:nnode)
  bestcode(1:nnode,1:nnode) = code(1:nnode,1:nnode)
!
!  Now consider all possible orderings, and their codes.
!
  more = .false.

  do

    call perm_next ( nnode, order, more, even )

    call g_order_code ( adj, nnode, nnode, code, order )

    call g_code_compare ( bestcode, code, nnode, nnode, result )

    ncomp = ncomp + 1

    if ( result == - 1 ) then

      nswap = nswap + 1

      bestorder(1:nnode) = order(1:nnode)
      bestcode(1:nnode,1:nnode) = code(1:nnode,1:nnode)

    end if

    if ( .not. more ) then
      exit
    end if

  end do
!
!  Once we have examined all possibilites, we are done.
!
!  Set the output ordering to the best ordering, and the output
!  code to the corresponding best code.
!
  order(1:nnode) = bestorder(1:nnode)
  code(1:nnode,1:nnode) = bestcode(1:nnode,1:nnode)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'G_CODE_BRUTE:'
  write ( *, '(a,i8)' ) '  Comparisons: ', ncomp
  write ( *, '(a,i8)' ) '  Swaps:       ', nswap

  return
end
subroutine g_code_cand ( adj, bestcode, code, nnode, ncomp, nopen, &
  order, stack, maxstack, nstack, ncan )

!*****************************************************************************80
!
!! G_CODE_CAND finds candidates for a maximal graph code ordering.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,J) is 1 if nodes I and J are adjacent.
!
!    Input, integer ( kind = 4 ) BESTCODE(NNODE,NNODE), the best code so far.
!
!    Workspace, integer CODE(NNODE,NNODE).
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input/output, integer ( kind = 4 ) NCOMP, the number of code comparisons.
!    This routine updates NCOMP by 1 each time the routine is called.
!
!    Input, integer ( kind = 4 ) NOPEN, identifies the first open position in ORDER.
!
!    Input, integer ( kind = 4 ) ORDER(NNODE), contains in entries 1 through NOPEN-1
!    the elements of the current partial list.
!
!    Input/output, integer ( kind = 4 ) STACK(MAXSTACK), used to store the new candidates.
!
!    Input, integer ( kind = 4 ) MAXSTACK, the maximum size of the STACK array.
!
!    Input/output, integer ( kind = 4 ) NSTACK, the current length of the stack.
!    On output, NSTACK has been increased by the number of
!    candidates found.
!
!    Input/output, integer ( kind = 4 ) NCAN(NNODE), candidates for the positions.
!
  implicit none

  integer ( kind = 4 ) maxstack
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) bestcode(nnode,nnode)
  integer ( kind = 4 ) code(nnode,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ifree(nnode)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) max_adj(nnode)
  integer ( kind = 4 ) max_max_adj
  integer ( kind = 4 ) ncan(nnode)
  integer ( kind = 4 ) ncomp
  integer ( kind = 4 ) nfree
  integer ( kind = 4 ) ni
  integer ( kind = 4 ) nj
  integer ( kind = 4 ) nopen
  integer ( kind = 4 ) nstack
  integer ( kind = 4 ) order(nnode)
  integer ( kind = 4 ) result
  integer ( kind = 4 ) stack(maxstack)

  if ( nopen < 1 .or. nnode < nopen ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'G_CODE_CAND - Fatal error 1!'
    write ( *, '(a)' ) '  1 <= NOPEN <= NNODE should be true, but'
    write ( *, '(a,i8)' ) '  NOPEN = ', nopen
    write ( *, '(a,i8)' ) '  NNODE = ', nnode
    stop
  end if

  ncan(nopen) = 0
!
!  If we have fixed at least one entry of the list, then compare
!  the code of the current incomplete ordering to the
!  code induced by the corresponding part of the current
!  best ordering.
!
!  If the current incomplete ordering is actually LESS than the
!  current best, then bail out with zero candidates.
!
  if ( 1 < nopen ) then

    call g_order_code ( adj, nnode, nopen-1, code, order )

    call g_code_compare ( bestcode, code, nnode, nopen-1, result )

    ncomp = ncomp + 1

    if ( result == + 1 ) then
      ncan(nopen) = 0
      return
    end if

  end if
!
!  Get a list of those nodes which have not been used yet.
!
  nfree = nnode + 1 - nopen
  call perm_free ( order, nopen-1, ifree, nfree )
!
!  Our candidates will be the unused neighbors of the lowest ordered node
!  possible.
!
  ncan(nopen) = 0

  do i = 1, nopen-1

    ni = order(i)

    do j = 1, nfree

      nj = ifree(j)

      if (  adj(ni,nj) /= 0 .or. adj(nj,ni) /= 0 ) then

        ncan(nopen) = ncan(nopen) + 1
        nstack = nstack + 1

        if ( maxstack < nstack ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'G_CODE_CAND - Fatal error 2!'
          write ( *, '(a)' ) '  MAXSTACK < NSTACK!'
          write ( *, '(a,i8)' ) '  NSTACK =   ', nstack
          write ( *, '(a,i8)' ) '  MAXSTACK = ', maxstack
          stop
        end if

        stack(nstack) = nj

      end if

    end do

    if ( 0 < ncan(nopen) ) then
      return
    end if

  end do
!
!  If we get here, no free nodes are connected in any way to
!  used nodes.  This can happen in two ways:
!
!  * NOPEN = 1;
!  * The graph is disconnected;
!
!  In either case, take as candidates the free nodes with at least one
!  neighbor (or maybe zero, if that's the best we can do.)
!
  max_max_adj = 0

  do i = 1, nfree

    ni = ifree(i)

    max_adj(i) = 0
    do j = 1, nfree
      nj = ifree(j)
      if ( ni /= nj ) then
        max_adj(i) = max ( max_adj(i), adj(ni,nj) )
      end if
    end do

    max_max_adj = max ( max_max_adj, max_adj(i) )

  end do

  ncan(nopen) = 0

  do i = 1, nfree

    if ( max_adj(i) == max_max_adj ) then

      ncan(nopen) = ncan(nopen) + 1
      nstack = nstack + 1

      if ( maxstack < nstack ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'G_CODE_CAND - Fatal error 2!'
        write ( *, '(a)' ) '  MAXSTACK < NSTACK!'
        write ( *, '(a,i8)' ) '  NSTACK =   ', nstack
        write ( *, '(a,i8)' ) '  MAXSTACK = ', maxstack
        stop
      end if

      stack(nstack) = ifree(i)

    end if
  end do

  return
end
subroutine g_code_compare ( code1, code2, nnode, npart, result )

!*****************************************************************************80
!
!! G_CODE_COMPARE compares two partial graph codes.
!
!  Discussion:
!
!    CODE1 = CODE2 if every digit of both codes is equal.
!
!    Otherwise, traverse the entries in a funny diagonal way, suggested
!    by this diagram for the first 10 entries:
!
!       *  1  2  4  7
!       *  *  3  5  8
!       *  *  *  6  9
!       *  *  *  * 10
!       *  *  *  *  *
!
!    As we do that, we examine the corresponding digits of the two codes.
!    For the first entry, (I,J), where the two codes differ, we say:
!
!      if ( CODE1(I,J) < CODE2(I,J) ) then we say
!        CODE1 < CODE2
!      else
!        CODE2 < CODE1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) CODE1(NNODE,NNODE), CODE2(NNODE,NNODE), codes to compare.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) NPART, specifies the portion of the codes
!    to compare.  NPART should be between 1 and NNODE.
!
!    If NPART = NNODE, then the full codes are compared.
!
!    If NPART < NNODE, then only entries corresponding to I and J
!    both less than or equal to NPART will be compared.
!
!    Output, integer ( kind = 4 ) RESULT:
!    -1, CODE1 < CODE2;
!     0, CODE1 = CODE2;
!    +1, CODE2 < CODE1.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) code1(nnode,nnode)
  integer ( kind = 4 ) code2(nnode,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) npart
  integer ( kind = 4 ) result

  do j = 2, npart

    do i = 1, j-1

      if ( code1(i,j) < code2(i,j) ) then

        result = - 1
        return

      else if ( code2(i,j) < code1(i,j) ) then

        result = + 1
        return

      end if

    end do

  end do

  result = 0

  return
end
subroutine g_code_print ( nnode, code, title )

!*****************************************************************************80
!
!! G_CODE_PRINT prints out a graph code.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) CODE(NNODE,NNODE), the code.
!
!    Input, character ( len = * ) TITLE, a title to be printed.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) code(nnode,nnode)
  integer ( kind = 4 ) i
  character ( len = 80 ) string
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  do i = 1, nnode
    string(i:i) = '.'
  end do

  do i = 1, nnode

    write ( *, '(2x,a,80i1)' ) string(1:i), code(i,i+1:nnode)

  end do

  return
end
subroutine g_compare ( adj1, nnode1, adj2, nnode2, order1, &
  order2, result )

!*****************************************************************************80
!
!! G_COMPARE determines if graphs G1 and G2 are isomorphic.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 October 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ1(NNODE1,NNODE1), the adjacency information
!    for G1.
!
!    Input, integer ( kind = 4 ) NNODE1, the number of nodes in G1.
!
!    Input, integer ( kind = 4 ) ADJ2(NNODE2,NNODE2), the adjacency information
!    for G2.
!
!    Input, integer ( kind = 4 ) NNODE2, the number of nodes in G2.
!
!    Output, integer ( kind = 4 ) RESULT, is 0 if G1 and G2 are isomorphic,
!    -I if G1 < G2 for test #I, and
!    +I if G2 < G1 for test #I.
!
!    Output, integer ( kind = 4 ) ORDER1(NNODE1), ORDER2(NNODE2).  If
!    RESULT = 0, then ORDER1 and ORDER2 are reorderings of the nodes of G1 and
!    G2 which exhibit the isomorphism.
!
  implicit none

  integer ( kind = 4 ) nnode1
  integer ( kind = 4 ) nnode2

  integer ( kind = 4 ) adj1(nnode1,nnode1)
  integer ( kind = 4 ) adj2(nnode2,nnode2)
  integer ( kind = 4 ) code1(nnode1,nnode1)
  integer ( kind = 4 ) code2(nnode2,nnode2)
  integer ( kind = 4 ) seq1(nnode1)
  integer ( kind = 4 ) seq2(nnode2)
  integer ( kind = 4 ) nedge1
  integer ( kind = 4 ) nedge2
  integer ( kind = 4 ) order1(nnode1)
  integer ( kind = 4 ) order2(nnode2)
  integer ( kind = 4 ) result
!
!  Test 1: Count the nodes.
!
  if ( nnode1 < nnode2 ) then
    result = - 1
    return
  else if ( nnode2 < nnode1 ) then
    result = + 1
    return
  end if
!
!  Test 2: Count the edges.
!
  call g_edge_count ( adj1, nnode1, nedge1 )

  call g_edge_count ( adj2, nnode2, nedge2 )

  if ( nedge1 < nedge2 ) then
    result = - 2
    return
  else if ( nedge2 < nedge1 ) then
    result = + 2
    return
  end if
!
!  Test 3: Compare the degree sequences.
!
  call g_degree_seq ( adj1, nnode1, seq1 )

  call g_degree_seq ( adj2, nnode2, seq2 )

  call i4vec_compare ( nnode1, seq1, seq2, result )

  if ( result < 0 ) then
    result = - 3
    return
  else if ( 0 < result ) then
    result = + 3
    return
  end if
!
!  Test 4: Compare the codes.
!
  call g_code_back ( adj1, nnode1, code1, order1 )

  call g_code_back ( adj2, nnode2, code2, order2 )

  call g_code_compare ( code1, code2, nnode1, nnode1, result )

  if ( result < 0 ) then
    result = - 4
    return
  else if ( 0 < result ) then
    result = + 4
    return
  end if

  result = 0

  return
end
subroutine g_connect_random ( adj, nedge, nnode, seed )

!*****************************************************************************80
!
!! G_CONNECT_RANDOM generates a random connected graph.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency matrix.  ADJ(I,J) is
!    nonzero if there is an edge from node I to node J.  ADJ(I,I) will
!    always be 0.
!
!    Input, integer ( kind = 4 ) NEDGE, the number of edges, which must be between
!    NNODE-1 and (NNODE*(NNODE-1))/2.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
  implicit none

  integer ( kind = 4 ) nnode
  integer ( kind = 4 ) nedge

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) code(nnode-2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) inode(nnode-1)
  integer ( kind = 4 ) iwork(nedge)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jnode(nnode-1)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) maxedge
  integer ( kind = 4 ) nchoice
  integer ( kind = 4 ) nchoose
  integer ( kind = 4 ) seed
!
!  Check.
!
  if ( nnode <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'G_CONNECT_RANDOM - Fatal error!'
    write ( *, '(a,i8)' ) '  NNODE = ', nedge
    write ( *, '(a)' ) '  but NNODE must be at least 1.'
    stop
  end if

  maxedge = ( nnode * ( nnode - 1 ) ) / 2

  if ( nedge < nnode-1 .or. maxedge < nedge ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'G_CONNECT_RANDOM - Fatal error!'
    write ( *, '(a,i8)' ) '  NEDGE = ', nedge
    write ( *, '(a)' ) '  but NEDGE must be at least 0, and '
    write ( *, '(a,i8)' ) '  no more than ', maxedge
    stop
  end if
!
!  Initialize the adjacency matrix.
!
  adj(1:nnode,1:nnode) = 0
!
!  Pick a random tree.
!
  call tree_arc_random ( nnode, code, inode, jnode, seed )
!
!  Convert information to adjacency form.
!
  call g_arc_to_g_adj ( nnode-1, inode, jnode, adj, nnode )
!
!  Now we have NEDGE - ( NNODE - 1 ) more edges to add.
!
  nchoice = ( nnode * ( nnode - 1 ) ) / 2 - ( nnode - 1 )
  nchoose = nedge - ( nnode - 1 )

  call ksub_random ( nchoice, nchoose, iwork, seed )

  k = 0
  l = 1
  do i = 1, nnode
    do j = i + 1, nnode
      if ( adj(i,j) /= 0 .or. adj(j,i) /= 0 ) then
        k = k + 1

        if ( l <= nchoose ) then
          if ( iwork(l) == k ) then
            adj(i,j) = 1
            adj(j,i) = 1
            l = l + 1
          end if
        end if

      end if
    end do
  end do

  return
end
subroutine g_degree ( adj, nnode, degree )

!*****************************************************************************80
!
!! G_DEGREE computes the degree of each node in a graph.
!
!  Discussion:
!
!    The degree of a node in a graph is the number of edges that are
!    incident on it.  The sum of the degrees of the nodes is twice the
!    number of edges.
!
!    The generalized case, where ADJ(I,J) can be greater than 1, indicating
!    the existence of 2 or more distinct edges between nodes I and J,
!    will be properly handled by this routine.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 October 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,J) is 1 if there is an edge from node I to node J.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) DEGREE(NNODE), the degree of the nodes.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) degree(nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  degree(1:nnode) = 0

  do i = 1, nnode
    do j = 1, nnode
      if ( adj(i,j) /= 0 ) then
        degree(i) = degree(i) + adj(i,j)
      end if
    end do
  end do

  return
end
subroutine g_degree_max ( adj, nnode, degree_max )

!*****************************************************************************80
!
!! G_DEGREE_MAX computes the maximum node degree of a graph.
!
!  Discussion:
!
!    The maximum node degree of a graph is the maximum value of the
!    degree of the nodes of the graph.
!
!    If two graphs are isomorphic, they must have the same maximum node degree.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 October 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,J) is 1 if there is an edge from node I to node J.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) DEGREE_MAX, the maximum node degree.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) degree
  integer ( kind = 4 ) degree_max
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  degree_max = 0

  do i = 1, nnode
    degree = 0
    do j = 1, nnode
      if ( adj(i,j) /= 0 ) then
        degree = degree + adj(i,j)
      end if
    end do
    degree_max = max ( degree_max, degree )
  end do

  return
end
subroutine g_degree_seq ( adj, nnode, seq )

!*****************************************************************************80
!
!! G_DEGREE_SEQ computes the degree sequence of a graph.
!
!  Discussion:
!
!    The degree sequence of a graph is constructed by computing the
!    degree of each node, and then ordering these values in decreasing order.
!
!    If two graphs are isomorphic, they must have the same degree sequence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 October 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,J) is 1 if there is an edge from node I to node J.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) SEQ(NNODE), the degree sequence.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) seq(nnode)

  seq(1:nnode) = 0

  do i = 1, nnode
    do j = 1, nnode
      seq(i) = seq(i) + adj(i,j)
    end do
  end do

  call i4vec_sort_heap_d ( nnode, seq )

  return
end
subroutine g_edge_count ( adj, nnode, nedge )

!*****************************************************************************80
!
!! G_EDGE_COUNT counts the number of edges in a graph.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 October 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,J) is 1 if there is an edge from node I to node J.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) NEDGE, the number of edges.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nedge

  nedge = 0

  do i = 1, nnode
    do j = 1, nnode

      if ( i == j ) then
        nedge = nedge + 2 * adj(i,j)
      else
        nedge = nedge + adj(i,j)
      end if

    end do
  end do

  nedge = nedge / 2

  return
end
subroutine g_example_bush ( adj, nnode )

!*****************************************************************************80
!
!! G_EXAMPLE_BUSH sets up the adjacency information for the bush graph.
!
!  Diagram:
!
!        6   3
!        |   |
!    1---4---5---2
!        |
!        7
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 October 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information for the graph.
!    ADJ(I,J) is 1 if nodes I and J are adjacent and 0 otherwise.
!
!    Output, integer ( kind = 4 ) NNODE, the number of nodes, which is 7.
!
  implicit none

  integer ( kind = 4 ) adj(7,7)
  integer ( kind = 4 ) nnode

  nnode = 7

  adj(1:nnode,1:nnode) = 0

  adj(1,4) = 1

  adj(2,5) = 1

  adj(3,5) = 1

  adj(4,1) = 1
  adj(4,5) = 1
  adj(4,6) = 1
  adj(4,7) = 1

  adj(5,2) = 1
  adj(5,3) = 1
  adj(5,4) = 1

  adj(6,4) = 1

  adj(7,4) = 1

  return
end
subroutine g_example_cube ( adj, nnode )

!*****************************************************************************80
!
!! G_EXAMPLE_CUBE sets up the adjacency information for the cube graph.
!
!  Diagram:
!
!      4-----7
!     /|    /|
!    8-----3 |
!    | |   | |
!    | 5---|-2
!    |/    |/
!    1-----6
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 October 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information for the graph.
!    ADJ(I,J) is 1 if nodes I and J are adjacent and 0 otherwise.
!
!    Output, integer ( kind = 4 ) NNODE, the number of nodes.
!
  implicit none

  integer ( kind = 4 ) adj(8,8)
  integer ( kind = 4 ) nnode

  nnode = 8

  adj(1:nnode,1:nnode) = 0

  adj(1,5) = 1
  adj(1,6) = 1
  adj(1,8) = 1

  adj(2,5) = 1
  adj(2,6) = 1
  adj(2,7) = 1

  adj(3,6) = 1
  adj(3,7) = 1
  adj(3,8) = 1

  adj(4,5) = 1
  adj(4,7) = 1
  adj(4,8) = 1

  adj(5,1) = 1
  adj(5,2) = 1
  adj(5,4) = 1

  adj(6,1) = 1
  adj(6,2) = 1
  adj(6,3) = 1

  adj(7,2) = 1
  adj(7,3) = 1
  adj(7,4) = 1

  adj(8,1) = 1
  adj(8,3) = 1
  adj(8,4) = 1

  return
end
subroutine g_example_dodecahedron ( adj, nnode )

!*****************************************************************************80
!
!! G_EXAMPLE_DODECAHEDRON sets adjacency for the dodecahedron graph.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 October 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,J) is 1 if nodes I and J are adjacent and 0 otherwise.
!
!    Output, integer ( kind = 4 ) NNODE, the number of nodes, which is 20.
!
  implicit none

  integer ( kind = 4 ) adj(20,20)
  integer ( kind = 4 ) nnode

  nnode = 20

  adj(1:nnode,1:nnode) = 0

  adj(1,2) = 1
  adj(1,5) = 1
  adj(1,6) = 1

  adj(2,1) = 1
  adj(2,3) = 1
  adj(2,8) = 1

  adj(3,2) = 1
  adj(3,4) = 1
  adj(3,10) = 1

  adj(4,3) = 1
  adj(4,5) = 1
  adj(4,12) = 1

  adj(5,1) = 1
  adj(5,4) = 1
  adj(5,14) = 1

  adj(6,1) = 1
  adj(6,7) = 1
  adj(6,15) = 1

  adj(7,6) = 1
  adj(7,8) = 1
  adj(7,17) = 1

  adj(8,7) = 1
  adj(8,9) = 1
  adj(8,2) = 1

  adj(9,8) = 1
  adj(9,10) = 1
  adj(9,16) = 1

  adj(10,3) = 1
  adj(10,9) = 1
  adj(10,11) = 1

  adj(11,10) = 1
  adj(11,12) = 1
  adj(11,20) = 1

  adj(12,4) = 1
  adj(12,11) = 1
  adj(12,13) = 1

  adj(13,12) = 1
  adj(13,14) = 1
  adj(13,19) = 1

  adj(14,13) = 1
  adj(14,15) = 1
  adj(14,5) = 1

  adj(15,6) = 1
  adj(15,14) = 1
  adj(15,18) = 1

  adj(16,9) = 1
  adj(16,17) = 1
  adj(16,20) = 1

  adj(17,16) = 1
  adj(17,18) = 1
  adj(17,7) = 1

  adj(18,15) = 1
  adj(18,17) = 1
  adj(18,19) = 1

  adj(19,13) = 1
  adj(19,18) = 1
  adj(19,20) = 1

  adj(20,11) = 1
  adj(20,16) = 1
  adj(20,19) = 1

  return
end
subroutine g_example_octo ( example, adj, nnode, seed )

!*****************************************************************************80
!
!! G_EXAMPLE_OCTO sets up an 8 node example graph.
!
!  Diagram:
!
!      1---2
!     /|   |\
!    8-+---+-3
!    | |   | |
!    7-+---+-4
!     \|   |/
!      6---5
!
!     Graph "A"
!
!    There are 8 graphs to choose from.  They are all on 8 nodes.  The first
!    5 have degree 3 at every node.  Graphs 6 and 7 have degree 5 at every
!    node.  Graph 8 is disconnected.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) EXAMPLE, chooses the example, between 1 and 8.
!
!    Output, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,J) is 1 if nodes I and J are adjacent and 0 otherwise.
!
!    Output, integer ( kind = 4 ) NNODE, the number of nodes, which should be 8.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
  implicit none

  integer ( kind = 4 ) adj(8,8)
  integer ( kind = 4 ) example
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nnode
  integer ( kind = 4 ) nsave
  integer ( kind = 4 ) seed

  if ( example <= 0 ) then
    nsave = i4_uniform ( 1, 8, seed )
  else
    example = mod ( example - 1, 8 ) + 1
    nsave = example
  end if

  nnode = 8

  adj(1:nnode,1:nnode) = 0

  do i = 1, nnode
    j = i + 1
    if ( nnode < j ) then
      j = j - nnode
    end if

    adj(i,j) = 1
    adj(j,i) = 1

  end do

  if ( nsave == 1 ) then

    adj(1,6) = 1
    adj(6,1) = 1
    adj(2,5) = 1
    adj(5,2) = 1
    adj(3,8) = 1
    adj(8,3) = 1
    adj(4,7) = 1
    adj(7,4) = 1

  else if ( nsave == 2 ) then

    adj(1,6) = 1
    adj(6,1) = 1
    adj(2,8) = 1
    adj(8,2) = 1
    adj(3,5) = 1
    adj(5,3) = 1
    adj(4,7) = 1
    adj(7,4) = 1

  else if ( nsave == 3 ) then

    adj(1,5) = 1
    adj(5,1) = 1
    adj(2,6) = 1
    adj(6,2) = 1
    adj(3,7) = 1
    adj(7,3) = 1
    adj(4,8) = 1
    adj(8,4) = 1

  else if ( nsave == 4 ) then

    adj(1,3) = 1
    adj(3,1) = 1
    adj(2,4) = 1
    adj(4,2) = 1
    adj(5,7) = 1
    adj(7,5) = 1
    adj(6,8) = 1
    adj(8,6) = 1

  else if ( nsave == 5 ) then

    adj(1,4) = 1
    adj(4,1) = 1
    adj(2,6) = 1
    adj(6,2) = 1
    adj(3,8) = 1
    adj(8,3) = 1
    adj(5,7) = 1
    adj(7,5) = 1

  else if ( nsave == 6 ) then

    adj(1,4) = 1
    adj(4,1) = 1
    adj(1,5) = 1
    adj(5,1) = 1
    adj(1,6) = 1
    adj(6,1) = 1
    adj(2,5) = 1
    adj(5,2) = 1
    adj(2,6) = 1
    adj(6,2) = 1
    adj(2,7) = 1
    adj(7,2) = 1
    adj(3,6) = 1
    adj(6,3) = 1
    adj(3,7) = 1
    adj(7,3) = 1
    adj(3,8) = 1
    adj(8,3) = 1
    adj(4,7) = 1
    adj(7,4) = 1
    adj(4,8) = 1
    adj(8,4) = 1
    adj(5,8) = 1
    adj(8,5) = 1

  else if ( nsave == 7 ) then

    adj(1,3) = 1
    adj(3,1) = 1
    adj(1,5) = 1
    adj(5,1) = 1
    adj(1,7) = 1
    adj(7,1) = 1
    adj(2,4) = 1
    adj(4,2) = 1
    adj(2,6) = 1
    adj(6,2) = 1
    adj(2,8) = 1
    adj(8,2) = 1
    adj(3,5) = 1
    adj(5,3) = 1
    adj(3,7) = 1
    adj(7,3) = 1
    adj(4,6) = 1
    adj(6,4) = 1
    adj(4,8) = 1
    adj(8,4) = 1
    adj(5,7) = 1
    adj(7,5) = 1
    adj(6,8) = 1
    adj(8,6) = 1

  else if ( nsave == 8 ) then

    adj(1,2) = 1
    adj(2,1) = 1
    adj(1,3) = 1
    adj(3,1) = 1
    adj(2,3) = 1
    adj(3,2) = 1
    adj(3,4) = 1
    adj(4,3) = 1
    adj(5,6) = 1
    adj(6,5) = 1
    adj(5,7) = 1
    adj(7,5) = 1
    adj(6,7) = 1
    adj(7,6) = 1

  end if
!
!  Now permute the graph.
!
  call i4mat_perm_random ( nnode, adj, seed )

  return
end
subroutine g_example_twig ( adj, nnode )

!*****************************************************************************80
!
!! G_EXAMPLE_TWIG sets up the adjacency information for the twig graph.
!
!  Diagram:
!
!    1---2---3
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 October 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,J) is 1 if nodes I and J are adjacent.
!
!    Output, integer ( kind = 4 ) NNODE, the number of nodes, which is 3.
!
  implicit none

  integer ( kind = 4 ) adj(3,3)
  integer ( kind = 4 ) nnode

  nnode = 3

  adj(1:nnode,1:nnode) = 0

  adj(1,2) = 1

  adj(2,1) = 1
  adj(2,3) = 1

  adj(3,2) = 1

  return
end
subroutine g_order_code ( adj, nnode, npart, code, order )

!*****************************************************************************80
!
!! G_ORDER_CODE returns the graph code for a specific node ordering.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,J) is 1 if nodes I and J are adjacent.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) NPART, the number of nodes to consider.
!    NPART should be between 1 and NNODE.
!
!    If NPART is NNODE, then the usual code is returned.
!
!    If NPART is less than NNODE, then the code is computed as
!    though only NPART nodes existed, namely, those specified in the
!    first NPART entries of order.  This option is provided so that
!    the routine can compute the portion of a code determined
!    by an incomplete ordering of the nodes.
!
!    Output, integer ( kind = 4 ) CODE(NNODE,NNODE), the code for this ordering.
!
!    Input, integer ( kind = 4 ) ORDER(NNODE), the ordering of the nodes.  ORDER(1)
!    is the first node, and so on.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) code(nnode,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ni
  integer ( kind = 4 ) nj
  integer ( kind = 4 ) npart
  integer ( kind = 4 ) order(nnode)

  do i = 1, nnode

    if ( i <= npart ) then

      ni = order(i)

      if ( order(i) < 1 .or. nnode < order(i) ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'G_ORDER_CODE - Fatal error!'
        write ( *, '(a)' ) '  ORDER is not a proper permutation.'
        stop
      end if

    else
      ni = 0
    end if

    do j = i + 1, nnode

      if ( j <= npart ) then

        nj = order(j)

        if ( order(j) < 1 .or. nnode < order(j) ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'G_ORDER_CODE - Fatal error!'
          write ( *, '(a)' ) '  ORDER is not a proper permutation.'
          stop
        end if

      else
        nj = 0
      end if

      if ( ni == 0 .or. nj == 0 ) then

        code(i,j) = 0
        code(j,i) = 0

      else if ( ( ni < nj .and. adj(ni,nj) /= 0 ) .or. &
                ( nj < ni .and. adj(nj,ni) /= 0 ) ) then

        code(i,j) = 1
        code(j,i) = 1

      else

        code(i,j) = 0
        code(j,i) = 0

      end if

    end do
  end do

  return
end
subroutine g_print ( adj, nnode, title )

!*****************************************************************************80
!
!! G_PRINT prints out an adjacency matrix.
!
!  Discussion:
!
!    This routine actually allows the entries of ADJ to have ANY value.
!    Values between 0 and 9 will be printed as is.  Other values will
!    be printed as '*'.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,J) is 1 if there is a direct connection FROM node I TO node J,
!    and is 0 otherwise.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, character ( len = * ) TITLE, a title to be printed.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  character ( len = 80 ) string
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  do i = 1, nnode

    jhi = min ( nnode, 80 )

    do j = 1, jhi

      if ( 0 <= adj(i,j) .and. adj(i,j) <= 9 ) then
        string(j:j) = char ( 48 + adj(i,j) )
      else
        string(j:j) = '*'
      end if

    end do

    write ( *, '(2x,a)' ) string(1:jhi)

  end do

  return
end
subroutine g_random ( adj, nnode, nedge, seed )

!*****************************************************************************80
!
!! G_RANDOM generates a random graph on NNODE nodes with NEDGE edges.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency matrix.  ADJ(I,J) is
!    nonzero if there is an edge from node I to node J.  ADJ(I,I) will
!    always be 0.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) NEDGE, the number of edges, which must be between
!    0 and (NNODE*(NNODE-1))/2.  (Note that each edge will be listed
!    twice in the adjacency matrix).
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
  implicit none

  integer ( kind = 4 ) nnode
  integer ( kind = 4 ) nedge

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iwork(nedge)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) maxedge
  integer ( kind = 4 ) seed
!
!  Check.
!
  if ( nnode <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'G_RANDOM - Fatal error!'
    write ( *, '(a,i8)' ) '  NNODE = ', nedge
    write ( *, '(a)' ) '  but NNODE must be at least 1.'
    stop
  end if

  maxedge = ( nnode * ( nnode - 1 ) ) / 2

  if ( nedge < 0 .or. maxedge < nedge ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'G_RANDOM - Fatal error!'
    write ( *, '(a,i8)' ) '  NEDGE = ', nedge
    write ( *, '(a)' ) '  but NEDGE must be at least 0, and '
    write ( *, '(a,i8)' ) '  no more than ', maxedge
    stop
  end if
!
!  Initialize the adjacency matrix.
!
  adj(1:nnode,1:nnode) = 0
!
!  Pick a random NEDGE subset of MAXEDGE.
!
  call ksub_random ( maxedge, nedge, iwork, seed )
!
!  The usable spots in the superdiagonal are numbered as follows:
!
!  * 1  2   3  ...  n-1
!  * * n+1 n+2 ... 2n-3
!  ...
!  * *  *   *  ... (n*(n-1))/2
!  * *  *   *  ...   *
!
  k = 0
  l = 1
  do i = 1, nnode-1
    do j = i+1, nnode

      k = k + 1

      if ( l <= nedge ) then

        if ( k == iwork(l) ) then
          adj(i,j) = 1
          adj(j,i) = 1
          l = l + 1
        end if

      end if

    end do
  end do

  return
end
subroutine i4_swap ( i, j )

!*****************************************************************************80
!
!! I4_SWAP switches two I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
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
function i4_uniform ( a, b, seed )

!*****************************************************************************80
!
!! I4_UNIFORM returns a scaled pseudorandom I4.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ) value.
!
!    The pseudorandom number will be scaled to be uniformly distributed
!    between A and B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 November 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley Interscience, page 95, 1998.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A, B, the limits of the interval.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, integer ( kind = 4 ) I4_UNIFORM, a number between A and B.
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) k
  real    ( kind = 4 ) r
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) value

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_UNIFORM - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + 2147483647
  end if

  r = real ( seed, kind = 4 ) * 4.656612875E-10
!
!  Scale R to lie between A-0.5 and B+0.5.
!
  r = ( 1.0E+00 - r ) * ( real ( min ( a, b ), kind = 4 ) - 0.5E+00 ) &
    +             r   * ( real ( max ( a, b ), kind = 4 ) + 0.5E+00 )
!
!  Use rounding to convert R to an integer between A and B.
!
  value = nint ( r, kind = 4 )

  value = max ( value, min ( a, b ) )
  value = min ( value, max ( a, b ) )

  i4_uniform = value

  return
end
subroutine i4mat_perm ( n, a, p )

!*****************************************************************************80
!
!! I4MAT_PERM permutes the rows and columns of a square integer matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input/output, integer ( kind = 4 ) A(N,N).
!    On input, the matrix to be permuted.
!    On output, the permuted matrix.
!
!    Input, integer ( kind = 4 ) P(N), the permutation.  P(I) is the new number of
!    row and column I.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) is
  integer ( kind = 4 ) it
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) k
  integer ( kind = 4 ) lc
  integer ( kind = 4 ) nc
  integer ( kind = 4 ) p(n)

  call perm_check ( n, p, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4MAT_PERM - Fatal error!'
    write ( *, '(a)') '  The input array does not represent'
    write ( *, '(a)' ) '  a proper permutation.  In particular, the'
    write ( *, '(a,i8)' ) '  array is missing the value ', ierror
    stop
  end if

  call perm_cycle ( n, p, is, nc, 1 )

  do i = 1, n

    i1 = - p(i)

    if ( 0 < i1 ) then

      lc = 0

      do

        i1 = p(i1)
        lc = lc + 1

        if ( i1 <= 0 ) then
          exit
        end if

      end do

      i1 = i

      do j = 1, n

        if ( p(j) <= 0 ) then

          j2 = j
          k = lc

          do

            j1 = j2
            it = a(i1,j1)

            do

              i1 = abs ( p(i1) )
              j1 = abs ( p(j1) )

              call i4_swap ( a(i1,j1), it )

              if ( j1 /= j2 ) then
                cycle
              end if

              k = k - 1

              if ( i1 == i ) then
                exit
              end if

            end do

            j2 = abs ( p(j2) )

            if ( k == 0 ) then
              exit
            end if

          end do

        end if

      end do

    end if

  end do
!
!  Restore the positive signs of the data.
!
  p(1:n) = abs ( p(1:n) )

  return
end
subroutine i4mat_perm_random ( n, a, seed )

!*****************************************************************************80
!
!! I4MAT_PERM_RANDOM selects a random permutation of an integer matrix.
!
!  Discussion:
!
!    The matrix is assumed to be square.  A single permutation is
!    applied to both rows and columns.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns in the array.
!
!    Input/output, integer ( kind = 4 ) A(N,N), the array to be permuted.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) j
  integer ( kind = 4 ) seed
!
!  Permute the rows and columns together.
!
  do i = 1, n

    i2 = i4_uniform ( i, n, seed )

    do j = 1, n
      call i4_swap ( a(i2,j), a(i,j) )
    end do

    do j = 1, n
      call i4_swap ( a(j,i2), a(j,i) )
    end do

  end do

  return
end
subroutine i4mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! I4MAT_PRINT prints an integer matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 June 2003
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
!    Input, integer ( kind = 4 ) A(M,N), the matrix to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * ) title

  ilo = 1
  ihi = m
  jlo = 1
  jhi = n

  call i4mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

  return
end
subroutine i4mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! I4MAT_PRINT_SOME prints some of an integer matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 November 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, integer ( kind = 4 ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ( kind = 4 ) ILO, JLO, the first row and column to print.
!
!    Input, integer ( kind = 4 ) IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 10
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  character ( len = 7 ) ctemp(incx)
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

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i7)') j
    end do

    write ( *, '(''  Col '',10a7)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) ' '

    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi

      do j2 = 1, inc

        j = j2lo - 1 + j2

        write ( ctemp(j2), '(i7)' ) a(i,j)

      end do

      write ( *, '(i5,1x,10a7)' ) i, ( ctemp(j), j = 1, inc )

    end do

  end do

  write ( *, '(a)' ) ' '

  return
end
subroutine i4mat_row_compare ( m, n, a1, a2, result )

!*****************************************************************************80
!
!! I4MAT_ROW_COMPARE compares two arrays of row vectors.
!
!  Discussion:
!
!    The arrays are compared by comparing the rows.
!    The rows are compared in the lexicographic sense.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, number of rows in the arrays.
!
!    Input, integer ( kind = 4 ) N, number of columns in the arrays.
!
!    Input, integer ( kind = 4 ) A1(M,N), A2(M,N), the arrays.
!
!    Output, integer ( kind = 4 ) RESULT:
!    -1, A1 < A2,
!     0, A1 = A2,
!    +1, A2 < A1.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a1(m,n)
  integer ( kind = 4 ) a2(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) result

  result = 0

  do i = 1, m

    do j = 1, n

      if ( a1(i,j) < a2(i,j) ) then
        result = - 1
        return
      else if ( a2(i,j) < a1(i,j) ) then
        result = + 1
        return
      end if

    end do

  end do

  return
end
subroutine i4row_compare ( m, n, a, i, j, isgn )

!*****************************************************************************80
!
!! I4ROW_COMPARE compares two rows of a integer array.
!
!  Example:
!
!    Input:
!
!  M = 3, N = 4, I = 2, J = 3
!
!  A = (
!    1  2  3  4
!    5  6  7  8
!    9 10 11 12 )
!
!    Output:
!
!  ISGN = -1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, integer ( kind = 4 ) A(M,N), an array of M rows of vectors of length N.
!
!    Input, integer ( kind = 4 ) I, J, the rows to be compared.
!    I and J must be between 1 and M.
!
!    Output, integer ( kind = 4 ) ISGN, the results of the comparison:
!    -1, row I < row J,
!     0, row I = row J,
!    +1, row J < row I.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
!
!  Check.
!
  if ( i < 1 .or. m < i ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4ROW_COMPARE - Fatal error!'
    write ( *, '(a)' ) '  Row index I is out of bounds.'
    write ( *, '(a,i8)' ) '  I = ', i
    write ( *, '(a,i8)' ) '  M = ', m
    stop
  end if

  if ( j < 1 .or. m < j ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4ROW_COMPARE - Fatal error!'
    write ( *, '(a)' ) '  Row index J is out of bounds.'
    write ( *, '(a,i8)' ) '  J = ', j
    write ( *, '(a,i8)' ) '  M = ', m
    stop
  end if

  isgn = 0

  if ( i == j ) then
    return
  end if

  k = 1

  do while ( k <= n )

    if ( a(i,k) < a(j,k) ) then
      isgn = - 1
      return
    else if ( a(j,k) < a(i,k) ) then
      isgn = + 1
      return
    end if

    k = k + 1

  end do

  return
end
subroutine i4row_sort_d ( m, n, a )

!*****************************************************************************80
!
!! I4ROW_SORT_D descending sorts the rows of an integer array.
!
!  Discussion:
!
!    In lexicographic order, the statement "X < Y", applied to two
!    vectors X and Y of length M, means that there is some index I, with
!    1 <= I <= M, with the property that
!
!      X(J) = Y(J) for J < I, and
!      X(I) < Y(I).
!
!    In other words, the first time they differ, X is smaller.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns of A.
!
!    Input/output, integer ( kind = 4 ) A(M,N).
!    On input, the array of M rows of N-vectors.
!    On output, the rows of A have been sorted in descending
!    lexicographic order.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
!
!  Initialize.
!
  i = 0
  indx = 0
  isgn = 0
  j = 0
!
!  Call the external heap sorter.
!
  do

    call sort_heap_external ( m, indx, i, j, isgn )
!
!  Interchange the I and J objects.
!
    if ( 0 < indx ) then

      call i4row_swap ( m, n, a, i, j )
!
!  Compare the I and J objects.
!
    else if ( indx < 0 ) then

      call i4row_compare ( m, n, a, i, j, isgn )
      isgn = - isgn

    else

      exit

    end if

  end do

  return
end
subroutine i4row_sort2_d ( m, n, a )

!*****************************************************************************80
!
!! I4ROW_SORT2_D descending sorts the elements of each row of an integer array.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of A.
!
!    Input, integer ( kind = 4 ) N, the number of columns of A, and the length
!    of a vector of data.
!
!    Input/output, integer ( kind = 4 ) A(M,N).
!    On input, the array of M rows of N-vectors.
!    On output, the elements of each row of A have been sorted in descending
!    order.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) irow
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
!
!  Initialize.
!
  do irow = 1, m

    i = 0
    indx = 0
    isgn = 0
    j = 0
!
!  Call the external heap sorter.
!
    do

      call sort_heap_external ( n, indx, i, j, isgn )
!
!  Interchange the I and J objects.
!
      if ( 0 < indx ) then

        call i4_swap ( a(irow,i), a(irow,j) )
!
!  Compare the I and J objects.
!
      else if ( indx < 0 ) then

        if ( a(irow,i) < a(irow,j) ) then
          isgn = + 1
        else
          isgn = - 1
        end if

      else

        exit

      end if

    end do

  end do

  return
end
subroutine i4row_swap ( m, n, a, irow1, irow2 )

!*****************************************************************************80
!
!! I4ROW_SWAP swaps two rows of an integer array.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input/output, integer ( kind = 4 ) A(M,N), an array of data.
!
!    Input, integer ( kind = 4 ) IROW1, IROW2, the two rows to swap.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) irow1
  integer ( kind = 4 ) irow2
  integer ( kind = 4 ) j
!
!  Check.
!
  if ( irow1 < 1 .or. m < irow1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4ROW_SWAP - Fatal error!'
    write ( *, '(a)' ) '  IROW1 is out of range.'
    stop
  end if

  if ( irow2 < 1 .or. m < irow2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4ROW_SWAP - Fatal error!'
    write ( *, '(a)' ) '  IROW2 is out of range.'
    stop
  end if

  if ( irow1 == irow2 ) then
    return
  end if

  do j = 1, n
    call i4_swap ( a(irow1,j), a(irow2,j) )
  end do

  return
end
subroutine i4vec_backtrack ( n, x, indx, k, nstack, stack, maxstack, ncan )

!*****************************************************************************80
!
!! I4VEC_BACKTRACK supervises a backtrack search for an integer vector.
!
!  Discussion:
!
!    The routine tries to construct an integer vector one index at a time,
!    using possible candidates as supplied by the user.
!
!    At any time, the partially constructed vector may be discovered to be
!    unsatisfactory, but the routine records information about where the
!    last arbitrary choice was made, so that the search can be
!    carried out efficiently, rather than starting out all over again.
!
!    First, call the routine with INDX = 0 so it can initialize itself.
!
!    Now, on each return from the routine, if INDX is:
!      1, you've just been handed a complete candidate vector;
!         Admire it, analyze it, do what you like.
!      2, please determine suitable candidates for position X(K).
!         Return the number of candidates in NCAN(K), adding each
!         candidate to the end of STACK, and increasing NSTACK.
!      3, you're done.  Stop calling the routine;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 July 2000
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of positions to be filled in the vector.
!
!    Input/output, integer ( kind = 4 ) X(N), the partial or complete candidate vector.
!
!    Input/output, integer ( kind = 4 ) INDX, a communication flag.
!    On input,
!      0 to start a search.
!    On output:
!      1, a complete output vector has been determined and returned in X(1:N);
!      2, candidates are needed for position X(K);
!      3, no more possible vectors exist.
!
!    Output, integer ( kind = 4 ) K, if INDX=2, the current vector index being considered.
!
!    Input/output, integer ( kind = 4 ) NSTACK, the current length of the stack.
!
!    Input, integer ( kind = 4 ) STACK(MAXSTACK), a list of all current candidates for
!    all positions 1 through K.
!
!    Input, integer ( kind = 4 ) MAXSTACK, the maximum length of the stack.
!
!    Input/output, integer ( kind = 4 ) NCAN(N), lists the current number of candidates for
!    positions 1 through K.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) maxstack

  integer ( kind = 4 ) indx
  integer ( kind = 4 ) k
  integer ( kind = 4 ) ncan(n)
  integer ( kind = 4 ) nstack
  integer ( kind = 4 ) stack(maxstack)
  integer ( kind = 4 ) x(n)
!
!  If this is the first call, request a candidate for position 1.
!
  if ( indx == 0 ) then
    k = 1
    nstack = 0
    indx = 2
    return
  end if
!
!  Examine the stack.
!
  do
!
!  If there are candidates for position K, take the first available
!  one off the stack, and increment K.
!
!  This may cause K to reach the desired value of N, in which case
!  we need to signal the user that a complete set of candidates
!  is being returned.
!
    if ( 0 < ncan(k) ) then

      x(k) = stack(nstack)
      nstack = nstack - 1

      ncan(k) = ncan(k) - 1

      if ( k /= n ) then
        k = k + 1
        indx = 2
      else
        indx = 1
      end if

      exit
!
!  If there are no candidates for position K, then decrement K.
!  If K is still positive, repeat the examination of the stack.
!
    else

      k = k - 1

      if ( k <= 0 ) then
        indx = 3
        exit
      end if

    end if

  end do

  return
end
subroutine i4vec_compare ( n, a1, a2, isgn )

!*****************************************************************************80
!
!! I4VEC_COMPARE compares two integer vectors.
!
!  Discussion:
!
!    The lexicographic ordering is used.
!
!  Example:
!
!    Input:
!
!      A1 = ( 2, 6, 2 )
!      A2 = ( 2, 8, 12 )
!
!    Output:
!
!      ISGN = -1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vectors.
!
!    Input, integer ( kind = 4 ) A1(N), A2(N), the vectors to be compared.
!
!    Output, integer ( kind = 4 ) ISGN, the results of the comparison:
!    -1, A1 < A2,
!     0, A1 = A2,
!    +1, A2 < A1.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a1(n)
  integer ( kind = 4 ) a2(n)
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) k

  isgn = 0

  k = 1

  do while ( k <= n )

    if ( a1(k) < a2(k) ) then
      isgn = - 1
      return
    else if ( a2(k) < a1(k) ) then
      isgn = + 1
      return
    end if

    k = k + 1

  end do

  return
end
subroutine i4vec_heap_a ( n, a )

!*****************************************************************************80
!
!! I4VEC_HEAP_A reorders an array of integers into an ascending heap.
!
!  Discussion:
!
!    An ascending heap is an array A with the property that, for every index J,
!    A(J) <= A(2*J) and A(J) <= A(2*J+1), (as long as the indices
!    2*J and 2*J+1 are legal).
!
!  Diagram:
!
!                  A(1)
!                /      \
!            A(2)         A(3)
!          /     \        /  \
!      A(4)       A(5)  A(6) A(7)
!      /  \       /   \
!    A(8) A(9) A(10) A(11)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the input array.
!
!    Input/output, integer ( kind = 4 ) A(N).
!    On input, an unsorted array.
!    On output, the array has been reordered into a heap.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ifree
  integer ( kind = 4 ) key
  integer ( kind = 4 ) m
!
!  Only nodes N/2 down to 1 can be "parent" nodes.
!
  do i = n/2, 1, -1
!
!  Copy the value out of the parent node.
!  Position IFREE is now "open".
!
    key = a(i)
    ifree = i

    do
!
!  Positions 2*IFREE and 2*IFREE + 1 are the descendants of position
!  IFREE.  (One or both may not exist because they exceed N.)
!
      m = 2 * ifree
!
!  Does the first position exist?
!
      if ( n < m ) then
        exit
      end if
!
!  Does the second position exist?
!
      if ( m + 1 <= n ) then
!
!  If both positions exist, take the smaller of the two values,
!  and update M if necessary.
!
        if ( a(m+1) < a(m) ) then
          m = m + 1
        end if

      end if
!
!  If the small descendant is smaller than KEY, move it up,
!  and update IFREE, the location of the free position, and
!  consider the descendants of THIS position.
!
      if ( key <= a(m) ) then
        exit
      end if

      a(ifree) = a(m)
      ifree = m

    end do
!
!  Once there is no more shifting to do, KEY moves into the free spot.
!
    a(ifree) = key

  end do

  return
end
subroutine i4vec_heap_d ( n, a )

!*****************************************************************************80
!
!! I4VEC_HEAP_D reorders an array of integers into an descending heap.
!
!  Discussion:
!
!    A descending heap is an array A with the property that, for every index J,
!    A(J) >= A(2*J) and A(J) >= A(2*J+1), (as long as the indices
!    2*J and 2*J+1 are legal).
!
!  Diagram:
!
!                  A(1)
!                /      \
!            A(2)         A(3)
!          /     \        /  \
!      A(4)       A(5)  A(6) A(7)
!      /  \       /   \
!    A(8) A(9) A(10) A(11)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the input array.
!
!    Input/output, integer ( kind = 4 ) A(N).
!    On input, an unsorted array.
!    On output, the array has been reordered into a heap.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ifree
  integer ( kind = 4 ) key
  integer ( kind = 4 ) m
!
!  Only nodes N/2 down to 1 can be "parent" nodes.
!
  do i = n/2, 1, -1
!
!  Copy the value out of the parent node.
!  Position IFREE is now "open".
!
    key = a(i)
    ifree = i

    do
!
!  Positions 2*IFREE and 2*IFREE + 1 are the descendants of position
!  IFREE.  (One or both may not exist because they exceed N.)
!
      m = 2 * ifree
!
!  Does the first position exist?
!
      if ( n < m ) then
        exit
      end if
!
!  Does the second position exist?
!
      if ( m + 1 <= n ) then
!
!  If both positions exist, take the larger of the two values,
!  and update M if necessary.
!
        if ( a(m) < a(m+1) ) then
          m = m + 1
        end if

      end if
!
!  If the large descendant is larger than KEY, move it up,
!  and update IFREE, the location of the free position, and
!  consider the descendants of THIS position.
!
      if ( a(m) <= key ) then
        exit
      end if

      a(ifree) = a(m)
      ifree = m

    end do
!
!  Once there is no more shifting to do, KEY moves into the free spot IFREE.
!
    a(ifree) = key

  end do

  return
end
subroutine i4vec_indicator ( n, a )

!*****************************************************************************80
!
!! I4VEC_INDICATOR sets an integer vector to the indicator vector.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
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
function i4vec_nonzero_count ( n, a )

!*****************************************************************************80
!
!! I4VEC_NONZERO_COUNT counts the nonzero entries in an integer vector
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the input array.
!
!    Input, integer ( kind = 4 ) A(N), an array.
!
!    Output, integer ( kind = 4 ) I4VEC_NONZERO_COUNT, the number of nonzero entries.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4vec_nonzero_count

  i4vec_nonzero_count = 0

  do i = 1, n
    if ( a(i) /= 0 ) then
      i4vec_nonzero_count = i4vec_nonzero_count + 1
    end if
  end do

  return
end
subroutine i4vec_order_type ( n, a, order )

!*****************************************************************************80
!
!! I4VEC_ORDER_TYPE determines if an integer array is (non)strictly ascending/descending.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries of the array.
!
!    Input, integer ( kind = 4 ) A(N), the array to be checked.
!
!    Output, integer ( kind = 4 ) ORDER, order indicator:
!    -1, no discernable order;
!    0, all entries are equal;
!    1, ascending order;
!    2, strictly ascending order;
!    3, descending order;
!    4, strictly descending order.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) order
!
!  Search for the first value not equal to A(1).
!
  i = 1

  do

    i = i + 1

    if ( n < i ) then
      order = 0
      return
    end if

    if ( a(1) < a(i) ) then

      if ( i == 2 ) then
        order = 2
      else
        order = 1
      end if

      exit

    else if ( a(i) < a(1) ) then

      if ( i == 2 ) then
        order = 4
      else
        order = 3
      end if

      exit

    end if

  end do
!
!  Now we have a "direction".  Examine subsequent entries.
!
  do while ( i < n )

    i = i + 1

    if ( order == 1 ) then

      if ( a(i) < a(i-1) ) then
        order = -1
        exit
      end if

    else if ( order == 2 ) then

      if ( a(i) < a(i-1) ) then
        order = -1
        exit
      else if ( a(i) == a(i-1) ) then
        order = 1
      end if

    else if ( order == 3 ) then

      if ( a(i-1) < a(i) ) then
        order = -1
        exit
      end if

    else if ( order == 4 ) then

      if ( a(i-1) < a(i) ) then
        order = -1
        exit
      else if ( a(i) == a(i-1) ) then
        order = 3
      end if

    end if

  end do

  return
end
subroutine i4vec_perm_random ( n, a, seed )

!*****************************************************************************80
!
!! I4VEC_PERM_RANDOM selects a random permutation of an integer vector.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 May 2005
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of objects to be permuted.
!
!    Input/output, integer ( kind = 4 ) A(N), the vector to be permuted.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) j
  integer ( kind = 4 ) seed

  do i = 1, n

    j = i4_uniform ( i, n, seed )

    call i4_swap ( a(i), a(j) )

  end do

  return
end
subroutine i4vec_print ( n, a, title )

!*****************************************************************************80
!
!! I4VEC_PRINT prints an integer vector.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 November 2000
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
  integer ( kind = 4 ) big
  integer ( kind = 4 ) i
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  big = maxval ( abs ( a(1:n) ) )

  write ( *, '(a)' ) ' '
  if ( big < 1000 ) then
    do i = 1, n
      write ( *, '(2x,i6,2x,i4)' ) i, a(i)
    end do
  else if ( big < 1000000 ) then
    do i = 1, n
      write ( *, '(2x,i6,2x,i7)' ) i, a(i)
    end do
  else
    do i = 1, n
      write ( *, '(2x,i6,2x,i12)' ) i, a(i)
    end do
  end if

  return
end
subroutine i4vec_sort_heap_a ( n, a )

!*****************************************************************************80
!
!! I4VEC_SORT_HEAP_A ascending sorts an integer array using heap sort.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input/output, integer ( kind = 4 ) A(N).
!    On input, the array to be sorted;
!    On output, the array has been sorted.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) n1

  if ( n <= 1 ) then
    return
  end if
!
!  1: Put A into descending heap form.
!
  call i4vec_heap_d ( n, a )
!
!  2: Sort A.
!
!  The largest object in the heap is in A(1).
!  Move it to position A(N).
!
  call i4_swap ( a(1), a(n) )
!
!  Consider the diminished heap of size N1.
!
  do n1 = n-1, 2, -1
!
!  Restore the heap structure of A(1) through A(N1).
!
    call i4vec_heap_d ( n1, a )
!
!  Take the largest object from A(1) and move it to A(N1).
!
    call i4_swap ( a(1), a(n1) )

  end do

  return
end
subroutine i4vec_sort_heap_d ( n, a )

!*****************************************************************************80
!
!! I4VEC_SORT_HEAP_D descending sorts an integer array using heap sort.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input/output, integer ( kind = 4 ) A(N).
!    On input, the array to be sorted;
!    On output, the array has been sorted.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) n1

  if ( n <= 1 ) then
    return
  end if
!
!  1: Put A into ascending heap form.
!
  call i4vec_heap_a ( n, a )
!
!  2: Sort A.
!
!  The smallest object in the heap is in A(1).
!  Move it to position A(N).
!
  call i4_swap ( a(1), a(n) )
!
!  Consider the diminished heap of size N1.
!
  do n1 = n-1, 2, -1
!
!  Restore the heap structure of A(1) through A(N1).
!
    call i4vec_heap_a ( n1, a )
!
!  Take the smallest object from A(1) and move it to A(N1).
!
    call i4_swap ( a(1), a(n1) )

  end do

  return
end
subroutine i4vec_sorted_unique_count ( n, a, nuniq )

!*****************************************************************************80
!
!! I4VEC_SORTED_UNIQUE_COUNT counts unique elements in a sorted integer array.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements in IA.
!
!    Input/output, integer ( kind = 4 ) A(N).  On input, the sorted
!    integer ( kind = 4 ) array.  On output, the unique elements in IA.
!
!    Output, integer ( kind = 4 ) NUNIQ, the number of unique elements in IA.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) itest
  integer ( kind = 4 ) nuniq

  nuniq = 0

  if ( n <= 0 ) then
    return
  end if

  nuniq = 1

  if ( n == 1 ) then
    return
  end if

  itest = 1

  do

    itest = itest + 1

    if ( n < itest ) then
      return
    end if

    if ( a(itest) /= a(nuniq) ) then
      nuniq = nuniq + 1
      a(nuniq) = a(itest)
    end if

  end do

  return
end
subroutine i4vec_uniform ( n, a, b, seed, x )

!*****************************************************************************80
!
!! I4VEC_UNIFORM returns a scaled pseudorandom I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of integer ( kind = 4 ) values.
!
!    The pseudorandom numbers should be scaled to be uniformly distributed
!    between A and B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 November 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the vector.
!
!    Input, integer ( kind = 4 ) A, B, the limits of the interval.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, integer ( kind = 4 ) X(N), a vector of numbers between A and B.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  real    ( kind = 4 ) r
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) value
  integer ( kind = 4 ) x(n)

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4VEC_UNIFORM - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + 2147483647
    end if

    r = real ( seed, kind = 4 ) * 4.656612875E-10
!
!  Scale R to lie between A-0.5 and B+0.5.
!
    r = ( 1.0E+00 - r ) * ( real ( min ( a, b ), kind = 4 ) - 0.5E+00 ) &
      +             r   * ( real ( max ( a, b ), kind = 4 ) + 0.5E+00 )
!
!  Use rounding to convert R to an integer between A and B.
!
    value = nint ( r, kind = 4 )

    value = max ( value, min ( a, b ) )
    value = min ( value, max ( a, b ) )

    x(i) = value

  end do

  return
end
subroutine i4vec2_compare ( n, ivec, jvec, i, j, isgn )

!*****************************************************************************80
!
!! I4VEC2_COMPARE compares pairs of integers stored in two vectors.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 October 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data items.
!
!    Input, integer ( kind = 4 ) IVEC(N), JVEC(N), contain the two components of each item.
!
!    Input, integer ( kind = 4 ) I, J, the items to be compared.
!
!    Output, integer ( kind = 4 ) ISGN, the results of the comparison:
!    -1, item I is less than item J,
!     0, item I is equal to item J,
!    +1, item I is greater than item J.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) ivec(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jvec(n)

  isgn = 0

  if ( ivec(i) < ivec(j) ) then
    isgn = -1
  else if ( ivec(i) == ivec(j) ) then
    if ( jvec(i) < jvec(j) ) then
      isgn = -1
    else if ( jvec(i) < jvec(j) ) then
      isgn = 0
    else if ( jvec(j) < jvec(i) ) then
      isgn = +1
    end if
  else if ( ivec(j) < ivec(i) ) then
    isgn = +1
  end if

  return
end
subroutine i4vec2_sort_d ( n, ivec, jvec )

!*****************************************************************************80
!
!! I4VEC2_SORT_D descending sorts a vector of pairs of integers.
!
!  Discussion:
!
!    Each item to be sorted is a pair of integers (I,J), with the I
!    and J values stored in separate vectors IVEC and JVEC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 October 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items of data.
!
!    Input/output, integer ( kind = 4 ) IVEC(N), JVEC(N), the data to be sorted..
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ivec(n)
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jvec(n)
!
!  Initialize.
!
  i = 0
  indx = 0
  isgn = 0
  j = 0
!
!  Call the external heap sorter.
!
  do

    call sort_heap_external ( n, indx, i, j, isgn )
!
!  Interchange the I and J objects.
!
    if ( 0 < indx ) then

      call i4_swap ( ivec(i), ivec(j) )
      call i4_swap ( jvec(i), jvec(j) )
!
!  Compare the I and J objects.
!
    else if ( indx < 0 ) then

      call i4vec2_compare ( n, ivec, jvec, i, j, isgn )
      isgn = - isgn

    else

      exit

    end if

  end do

  return
end
subroutine i4vec2_uniq ( n, ivec, jvec, nuniq )

!*****************************************************************************80
!
!! I4VEC2_UNIQ keeps the unique elements in a array of pairs of integers.
!
!  Discussion:
!
!    Item I is stored as the pair IVEC(I), JVEC(I).
!
!    The items must have been sorted, or at least it must be the
!    case that equal items are stored in adjacent vector locations.
!
!    If the items were not sorted, then this routine will only
!    replace a string of equal values by a single representative.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 October 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items.
!
!    Input/output, integer ( kind = 4 ) IVEC(N), JVEC(N).
!    On input, the array of N items.
!    On output, an array of NUNIQ unique items.
!
!    Output, integer ( kind = 4 ) NUNIQ, the number of unique items.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) itest
  integer ( kind = 4 ) ivec(n)
  integer ( kind = 4 ) jvec(n)
  integer ( kind = 4 ) nuniq

  nuniq = 0

  if ( n <= 0 ) then
    return
  end if

  nuniq = 1

  if ( n == 1 ) then
    return
  end if

  itest = 1

  do

    itest = itest + 1

    if ( n < itest ) then
      return
    end if

    if ( ivec(itest) /= ivec(nuniq) .or. jvec(itest) /= jvec(nuniq) ) then

      nuniq = nuniq + 1

      ivec(nuniq) = ivec(itest)
      jvec(nuniq) = jvec(itest)

    end if

  end do

  return
end
subroutine ksub_random ( n, k, a, seed )

!*****************************************************************************80
!
!! KSUB_RANDOM selects a random subset of size K from a set of size N.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the set from which subsets are drawn.
!
!    Input, integer ( kind = 4 ) K, number of elements in desired subsets.  K must
!    be between 0 and N.
!
!    Output, integer ( kind = 4 ) A(K).  A(I) is the I-th element of the
!    output set.  The elements of A are in order.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
  implicit none

  integer ( kind = 4 ) k

  integer ( kind = 4 ) a(k)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) ids
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) ir
  integer ( kind = 4 ) is
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) l
  integer ( kind = 4 ) ll
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m0
  integer ( kind = 4 ) n
  integer ( kind = 4 ) seed

  if ( k < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KSUB_RANDOM - Fatal error!'
    write ( *, '(a,i8)' ) '  K = ', k
    write ( *, '(a)' ) '  but 0 <= K is required!'
    stop
  end if

  if ( n < k ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KSUB_RANDOM - Fatal error!'
    write ( *, '(a,i8)' ) '  N = ', n
    write ( *, '(a,i8)' ) '  and K = ', k
    write ( *, '(a)' ) '  K <= N is required!'
    stop
  end if

  if ( k == 0 ) then
    return
  end if

  do i = 1, k
    a(i) = ( ( i - 1 ) * n ) / k
  end do

  do i = 1, k

    do

      ix = i4_uniform ( 1, n, seed )

      l = 1 + ( ix * k - 1 ) / n

      if ( a(l) < ix ) then
        exit
      end if

    end do

    a(l) = a(l) + 1

  end do

  ip = 0
  is = k

  do i = 1, k

    m = a(i)
    a(i) = 0

    if ( m /= ( (i-1) * n ) / k ) then
      ip = ip + 1
      a(ip) = m
    end if

  end do

  ihi = ip

  do i = 1, ihi
    ip = ihi + 1 - i
    l = 1 + ( a(ip) * k - 1 ) / n
    ids = a(ip) - ( ( l - 1 ) * n ) / k
    a(ip) = 0
    a(is) = l
    is = is - ids
  end do

  do ll = 1, k

    l = k + 1 - ll

    if ( a(l) /= 0 ) then
      ir = l
      m0 = 1 + ( ( a(l) - 1 ) * n ) / k
      m = ( a(l) * n ) / k - m0 + 1
    end if

    ix = i4_uniform ( m0, m0 + m - 1, seed )
    i = l + 1

    do while ( i <= ir )

      if ( ix < a(i) ) then
        exit
      end if

      ix = ix + 1
      a(i-1) = a(i)
      i = i + 1

    end do

    a(i-1) = ix
    m = m - 1

  end do

  return
end
subroutine mg_adj_max_max ( adj, nnode, adj_max_max )

!*****************************************************************************80
!
!! MG_ADJ_MAX_MAX computes the adjacency maximum maximum of a multigraph.
!
!  Discussion:
!
!    The adjacency maximum maximum of a multigraph may be constructed by
!    computing the maximum entry of the adjacency matrix,
!
!    If two multigraphs are isomorphic, they must have the same
!    adjacency maximum maximum.
!
!  Example:
!
!    ADJ =
!       0 1 2 3
!       1 0 2 0
!       2 2 0 1
!       3 0 1 0
!
!    ADJ_MAX_MAX = 3
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,J) is the number of edges from node I to node J.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) ADJ_MAX_MAX, the adjacency maximum maximum.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) adj_max_max
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  adj_max_max = 0
  do i = 1, nnode
    do j = 1, nnode
      if ( i /= j ) then
        adj_max_max = max ( adj_max_max, adj(i,j) )
      end if
    end do
  end do

  return
end
subroutine mg_adj_max_seq ( adj, nnode, adj_max_seq )

!*****************************************************************************80
!
!! MG_ADJ_MAX_SEQ computes the adjacency maximum sequence of a multigraph.
!
!  Discussion:
!
!    The adjacency maximum sequence of a multigraph may be constructed by
!    computing the maximum entry of each row of the adjacency matrix,
!    and then sorting these values in descending order.
!
!    If two multigraphs are isomorphic, they must have the same
!    adjacency maximum sequence.
!
!  Example:
!
!    ADJ =
!       0 1 2 3
!       1 0 2 0
!       2 2 0 1
!       3 0 1 0
!
!    ADJ_MAX_SEQ =
!
!       3
!       3
!       2
!       2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,J) is the number of edges from node I to node J.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) ADJ_MAX_SEQ(NNODE), the adjacency maximum sequence.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) adj_max_seq(nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
!
!  Copy the adjacency matrix.
!
  do i = 1, nnode
    k = 0
    do j = 1, nnode
      if ( i /= j ) then
        k = max ( k, adj(i,j) )
      end if
    end do
    adj_max_seq(i) = k
  end do
!
!  Sort the elements.
!
  call i4vec_sort_heap_d ( nnode, adj_max_seq )

  return
end
subroutine mg_adj_seq ( adj, nnode, adj_seq )

!*****************************************************************************80
!
!! MG_ADJ_SEQ computes the adjacency sequence of a multigraph.
!
!  Discussion:
!
!    The adjacency sequence of a multigraph may be constructed by sorting the
!    entries of each row of the adjacency matrix in descending order, and
!    then sorting the rows themselves in descending order.
!
!    If two multigraphs are isomorphic, they must have the same adjacency
!    sequence.
!
!  Example:
!
!    ADJ =
!       0 1 2 3
!       1 0 2 0
!       2 2 0 1
!       3 0 1 0
!
!    ADJ_SEQ =
!
!       3 2 1 0
!       3 1 0 0
!       2 2 1 0
!       2 1 0 0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,J) is the number of edges from node I to node J.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) ADJ_SEQ(NNODE,NNODE), the adjacency sequence.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) adj_seq(nnode,nnode)
!
!  Copy the adjacency matrix.
!
  adj_seq(1:nnode,1:nnode) = adj(1:nnode,1:nnode)
!
!  Sort the elements of each row.
!
  call i4row_sort2_d ( nnode, nnode, adj_seq )
!
!  Sort the rows of the matrix.
!
  call i4row_sort_d ( nnode, nnode, adj_seq )

  return
end
subroutine mg_code_back ( adj, nnode, code, order )

!*****************************************************************************80
!
!! MG_CODE_BACK computes a multigraph code via backtracking.
!
!  Discussion:
!
!    The code is the "largest" order code (in the lexicographic
!    sense) over all possible node orderings.  The backtracking method
!    organizes the search of all possible node orderings so that if
!    a partial node ordering is sure to generate an order code
!    less than the best so far, then all the orderings that begin with
!    this partial ordering are immediately dropped from consideration.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,J) is the number of edges from node I to node J.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) CODE(NNODE,NNODE), the code.
!
!    Output, integer ( kind = 4 ) ORDER(NNODE), the ordering of the nodes.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) bestcode(nnode,nnode)
  integer ( kind = 4 ) bestorder(nnode)
  integer ( kind = 4 ) code(nnode,nnode)
  logical, parameter :: debug = .false.
  integer ( kind = 4 ) index
  integer ( kind = 4 ) maxstack
  integer ( kind = 4 ) ncan(nnode)
  integer ( kind = 4 ) ncomp
  integer ( kind = 4 ) nopen
  integer ( kind = 4 ) nstack
  integer ( kind = 4 ) nswap
  integer ( kind = 4 ) order(nnode)
  integer ( kind = 4 ) result
  integer ( kind = 4 ), allocatable, dimension ( : ) :: stack

  if ( nnode <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)') 'MG_CODE_BACK - Fatal error!'
    write ( *, '(a,i8)' ) '  NNODE = ', nnode
    stop
  end if

  maxstack = 10 * nnode
  allocate ( stack(1:maxstack) )

  nstack = 0
  stack(1) = 0

  ncomp = 0
  nswap = 0
!
!  Start with the identity ordering.
!
  call i4vec_indicator ( nnode, order )
!
!  Compute the corresponding code.
!
  call mg_order_code ( adj, nnode, nnode, code, order )
!
!  Take this ordering and code as the best so far.
!
  bestorder(1:nnode) = order(1:nnode)
  bestcode(1:nnode,1:nnode) = code(1:nnode,1:nnode)
!
!  Now consider all possible orderings, and their codes.
!
  index = 0

  do

    call i4vec_backtrack ( nnode, order, index, nopen, nstack, stack, &
      maxstack, ncan )
!
!  If the backtrack routine has returned a complete candidate ordering, then
!  compute the resulting code, and see it it is better
!  then our current best.  Then go back for the next backtrack search.
!
    if ( index == 1 ) then

      call mg_order_code ( adj, nnode, nnode, code, order )

      call mg_code_compare ( bestcode, code, nnode, nnode, result )

      ncomp = ncomp + 1

      if ( result == -1 ) then

        nswap = nswap + 1

        bestorder(1:nnode) = order(1:nnode)
        bestcode(1:nnode,1:nnode) = code(1:nnode,1:nnode)

      end if
!
!  If the backtrack routine has returned a partial reordering,
!  supply candidates for the next item in the ordering.
!
    else if ( index == 2 ) then

      call mg_code_cand ( adj, bestcode, code, nnode, ncomp, nopen, &
        order, stack, maxstack, nstack, ncan )

    else

      exit

    end if

  end do
!
!  Once we have examined all possibilites, we are done.
!
!  Set the output ordering to the best ordering, and the output
!  code to the corresponding best code.
!
  order(1:nnode) = bestorder(1:nnode)
  code(1:nnode,1:nnode) = bestcode(1:nnode,1:nnode)

  if ( debug ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MG_CODE_BACK:'
    write ( *, '(a,i8)' ) '  Comparisons: ', ncomp
    write ( *, '(a,i8)' ) '  Swaps:       ', nswap
  end if

  deallocate ( stack )

  return
end
subroutine mg_code_brute ( adj, nnode, code, order )

!*****************************************************************************80
!
!! MG_CODE_BRUTE computes a multigraph code via brute force.
!
!  Discussion:
!
!    The code is the "largest" order code (in the lexicographic sense)
!    over all possible node orderings.  The brute force method considers
!    every node ordering, computes the corresponding order code, and
!    takes the largest one encountered.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,J) is the number of edges from node I to node J.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) CODE(NNODE,NNODE), the code.
!
!    Output, integer ( kind = 4 ) ORDER(NNODE), the ordering of the nodes.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) bestcode(nnode,nnode)
  integer ( kind = 4 ) bestorder(nnode)
  integer ( kind = 4 ) code(nnode,nnode)
  logical even
  logical more
  integer ( kind = 4 ) ncomp
  integer ( kind = 4 ) nswap
  integer ( kind = 4 ) order(nnode)
  integer ( kind = 4 ) result

  ncomp = 0
  nswap = 0
!
!  Start with the identity ordering.
!
  call i4vec_indicator ( nnode, order )
!
!  Compute the corresponding code.
!
  call mg_order_code ( adj, nnode, nnode, code, order )
!
!  Take this ordering and code as the best so far.
!
  bestorder(1:nnode) = order(1:nnode)
  bestcode(1:nnode,1:nnode) = code(1:nnode,1:nnode)
!
!  Now consider all possible orderings, and their codes.
!
  more = .false.

  do

    call perm_next ( nnode, order, more, even )

    call mg_order_code ( adj, nnode, nnode, code, order )

    call mg_code_compare ( bestcode, code, nnode, nnode, result )

    ncomp = ncomp + 1

    if ( result == - 1 ) then

      nswap = nswap + 1

      bestorder(1:nnode) = order(1:nnode)
      bestcode(1:nnode,1:nnode) = code(1:nnode,1:nnode)

    end if

    if ( .not. more ) then
      exit
    end if

  end do
!
!  Once we have examined all possibilites, we are done.
!
!  Set the output ordering to the best ordering, and the output
!  code to the corresponding best code.
!
  order(1:nnode) = bestorder(1:nnode)
  code(1:nnode,1:nnode) = bestcode(1:nnode,1:nnode)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MG_CODE_BRUTE:'
  write ( *, '(a,i8)' ) '  Comparisons: ', ncomp
  write ( *, '(a,i8)' ) '  Swaps:       ', nswap

  return
end
subroutine mg_code_cand ( adj, bestcode, code, nnode, ncomp, nopen, &
  order, stack, maxstack, nstack, ncan )

!*****************************************************************************80
!
!! MG_CODE_CAND finds candidates for a maximal multigraph code ordering.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,J) is the number of edges from node I to node J.
!
!    Input, integer ( kind = 4 ) BESTCODE(NNODE,NNODE), the best code so far.
!
!    Workspace, integer CODE(NNODE,NNODE).
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input/output, integer ( kind = 4 ) NCOMP, the number of code comparisons.
!    This routine updates NCOMP by 1 each time the routine is called.
!
!    Input, integer ( kind = 4 ) NOPEN, identifies the first open position in ORDER.
!
!    Input, integer ( kind = 4 ) ORDER(NNODE), contains in entries 1 through NOPEN-1
!    the elements of the current partial list.
!
!    Input/output, integer ( kind = 4 ) STACK(MAXSTACK), used to store the new candidates.
!
!    Input, integer ( kind = 4 ) MAXSTACK, the maximum size of the STACK array.
!
!    Input/output, integer ( kind = 4 ) NSTACK, the current length of the stack.
!    On output, NSTACK has been increased by the number of
!    candidates found.
!
!    Input/output, integer ( kind = 4 ) NCAN(NNODE), the number of candidates
!    for each position.
!
  implicit none

  integer ( kind = 4 ) maxstack
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) bestcode(nnode,nnode)
  integer ( kind = 4 ) code(nnode,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ifree(nnode)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) max_adj(nnode)
  integer ( kind = 4 ) max_adj_ni
  integer ( kind = 4 ) max_max_adj
  integer ( kind = 4 ) ncan(nnode)
  integer ( kind = 4 ) ncomp
  integer ( kind = 4 ) nfree
  integer ( kind = 4 ) ni
  integer ( kind = 4 ) nj
  integer ( kind = 4 ) nopen
  integer ( kind = 4 ) nstack
  integer ( kind = 4 ) order(nnode)
  integer ( kind = 4 ) result
  integer ( kind = 4 ) stack(maxstack)

  if ( nopen < 1 .or. nnode < nopen ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MG_CODE_CAND - Fatal error 1!'
    write ( *, '(a)' ) '  1 <= NOPEN <= NNODE should be true, but'
    write ( *, '(a,i8)' ) '  NOPEN = ', nopen
    write ( *, '(a,i8)' ) '  NNODE = ', nnode
    stop
  end if

  ncan(nopen) = 0
!
!  If we have fixed at least one entry of the list, then compare
!  the code of the current incomplete ordering to the
!  code induced by the corresponding part of the current
!  best ordering.
!
!  If the current incomplete ordering is actually LESS than the
!  current best, then bail out with zero candidates.
!
  if ( 1 < nopen ) then

    call mg_order_code ( adj, nnode, nopen-1, code, order )

    call mg_code_compare ( bestcode, code, nnode, nopen-1, result )

    ncomp = ncomp + 1

    if ( result == +1 ) then
      ncan(nopen) = 0
      return
    end if

  end if
!
!  Get a list of those nodes which have not been used yet.
!
  nfree = nnode + 1 - nopen
  call perm_free ( order, nopen-1, ifree, nfree )
!
!  To find candidates, we consider all the ordered nodes.
!  We find the lowest ordered node that has unordered neighbors.
!  We take the unordered neighbors with maximal adjacency.
!
  ncan(nopen) = 0
!
!  For each ordered node NI...
!
  do i = 1, nopen-1

    ni = order(i)
!
!  ...note the maximum adjacency from NI to any unordered node NJ...
!
    max_adj_ni = 0
    do j = 1, nfree
      nj = ifree(j)
      max_adj_ni = max ( max_adj_ni, adj(ni,nj) )
    end do
!
!   ...and take as candidates all unordered nodes NJ with maximal
!   adjacency from NI.
!
    if ( 0 < max_adj_ni ) then

      do j = 1, nfree

        nj = ifree(j)

        if ( adj(ni,nj) == max_adj_ni ) then

          ncan(nopen) = ncan(nopen) + 1
          nstack = nstack + 1

          if ( maxstack < nstack ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'MG_CODE_CAND - Fatal error 2!'
            write ( *, '(a)' ) '  MAXSTACK < NSTACK!'
            write ( *, '(a,i8)' ) '  NSTACK =   ', nstack
            write ( *, '(a,i8)' ) '  MAXSTACK = ', maxstack
            stop
          end if

          stack(nstack) = nj

        end if

      end do

      return

    end if

  end do
!
!  If we get here, no unordered nodes are connected in any way to
!  ordered nodes.  This can happen in two ways:
!
!  * NOPEN = 1;
!  * The graph is disconnected;
!
!  For each free node, compute the maximum adjacency to any other free node.
!  Take the maximum of this value over all free nodes.
!  Candidates are free nodes with this maximum value.
!
  max_max_adj = 0

  do i = 1, nfree

    ni = ifree(i)

    max_adj(i) = 0
    do j = 1, nfree
      nj = ifree(j)
      if ( ni /= nj ) then
        max_adj(i) = max ( max_adj(i), adj(ni,nj) )
      end if
    end do

    max_max_adj = max ( max_max_adj, max_adj(i) )

  end do
!
!  Now go back and compute the maximum adjacency for each free node.
!  Nodes which achieve the maximum are added to the candidate list.
!
  ncan(nopen) = 0

  do i = 1, nfree

    if ( max_adj(i) == max_max_adj ) then

      ncan(nopen) = ncan(nopen) + 1
      nstack = nstack + 1

      if ( maxstack < nstack ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'MG_CODE_CAND - Fatal error 2!'
        write ( *, '(a)' ) '  MAXSTACK < NSTACK!'
        write ( *, '(a,i8)' ) '  NSTACK =   ', nstack
        write ( *, '(a,i8)' ) '  MAXSTACK = ', maxstack
        stop
      end if

      stack(nstack) = ifree(i)

    end if

  end do

  return
end
subroutine mg_code_compare ( code1, code2, nnode, npart, result )

!*****************************************************************************80
!
!! MG_CODE_COMPARE compares two partial multigraph codes.
!
!  Discussion:
!
!    CODE1 = CODE2 if every digit of both codes is equal.
!
!    Otherwise, traverse the entries in a funny diagonal way, suggested
!    by this diagram for the first 10 entries:
!
!       *  1  2  4  7
!       *  *  3  5  8
!       *  *  *  6  9
!       *  *  *  * 10
!       *  *  *  *  *
!
!    As we do that, we examine the corresponding digits of the two codes.
!    For the first entry, (I,J), where the two codes differ, we say:
!
!      if ( CODE1(I,J) < CODE2(I,J) ) then we say
!        CODE1 < CODE2
!      else
!        CODE2 < CODE1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) CODE1(NNODE,NNODE), CODE2(NNODE,NNODE), codes to compare.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) NPART, specifies the portion of the codes
!    to compare.  NPART should be between 1 and NNODE.
!
!    If NPART = NNODE, then the full codes are compared.
!
!    If NPART < NNODE, then only entries corresponding to I and J
!    both less than or equal to NPART will be compared.
!
!    Output, integer ( kind = 4 ) RESULT:
!    -1, CODE1 < CODE2;
!     0, CODE1 = CODE2;
!    +1, CODE2 < CODE1.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) code1(nnode,nnode)
  integer ( kind = 4 ) code2(nnode,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) npart
  integer ( kind = 4 ) result

  do j = 2, npart

    do i = 1, j-1

      if ( code1(i,j) < code2(i,j) ) then

        result = - 1
        return

      else if ( code2(i,j) < code1(i,j) ) then

        result = + 1
        return

      end if

    end do

  end do

  result = 0

  return
end
subroutine mg_code_print ( nnode, code, title )

!*****************************************************************************80
!
!! MG_CODE_PRINT prints out a multigraph code.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) CODE(NNODE,NNODE), the code.
!
!    Input, character ( len = * ) TITLE, a title to be printed.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) code(nnode,nnode)
  integer ( kind = 4 ) i
  character ( len = 80 ) string
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  do i = 1, nnode
    string(i:i) = '.'
  end do

  do i = 1, nnode

    write ( *, '(2x,a,80i1)' ) string(1:i), code(i,i+1:nnode)

  end do

  return
end
subroutine mg_compare ( adj1, nnode1, adj2, nnode2, order1, &
  order2, result )

!*****************************************************************************80
!
!! MG_COMPARE determines if multigraphs G1 and G2 are isomorphic.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ1(NNODE1,NNODE1), the adjacency information for G1.
!
!    Input, integer ( kind = 4 ) NNODE1, the number of nodes in G1.
!
!    Input, integer ( kind = 4 ) ADJ2(NNODE2,NNODE2), the adjacency information for G2.
!
!    Input, integer ( kind = 4 ) NNODE2, the number of nodes in G2.
!
!    Output, integer ( kind = 4 ) RESULT, is 0 if the multigraphs are isomorphic,
!    -I if G1 < G2 for test #I, and
!    +I if G2 < G1 for test #I.
!
!    Output, integer ( kind = 4 ) ORDER1(NNODE1), ORDER2(NNODE2).  If RESULT = 0, then ORDER1
!    and ORDER2 are reorderings of the nodes of G1 and G2 which
!    exhibit the isomorphism.
!
  implicit none

  integer ( kind = 4 ) nnode1
  integer ( kind = 4 ) nnode2

  integer ( kind = 4 ) adj_max_max_1
  integer ( kind = 4 ) adj_max_max_2
  integer ( kind = 4 ) adj1(nnode1,nnode1)
  integer ( kind = 4 ) adj2(nnode2,nnode2)
  integer ( kind = 4 ) code1(nnode1,nnode1)
  integer ( kind = 4 ) code2(nnode2,nnode2)
  integer ( kind = 4 ) seq1(nnode1)
  integer ( kind = 4 ) seq2(nnode2)
  integer ( kind = 4 ) nedge1
  integer ( kind = 4 ) nedge2
  integer ( kind = 4 ) order1(nnode1)
  integer ( kind = 4 ) order2(nnode2)
  integer ( kind = 4 ) result
!
!  Test 1: Count the nodes.
!
  if ( nnode1 < nnode2 ) then
    result = - 1
    return
  else if ( nnode2 < nnode1 ) then
    result = + 1
    return
  end if
!
!  Test 2: Count the edges.
!
  call mg_edge_count ( adj1, nnode1, nedge1 )

  call mg_edge_count ( adj2, nnode2, nedge2 )

  if ( nedge1 < nedge2 ) then
    result = - 2
    return
  else if ( nedge2 < nedge1 ) then
    result = + 2
    return
  end if
!
!  Test 3: Compare the degree sequences.
!
  call mg_degree_seq ( adj1, nnode1, seq1 )

  call mg_degree_seq ( adj2, nnode2, seq2 )

  call i4vec_compare ( nnode1, seq1, seq2, result )

  if ( result < 0 ) then
    result = - 3
    return
  else if ( 0 < result ) then
    result = + 3
    return
  end if
!
!  Test 4: Compare the adjacency max max.
!
  call mg_adj_max_max ( adj1, nnode1, adj_max_max_1 )

  call mg_adj_max_max ( adj2, nnode2, adj_max_max_2 )

  if ( adj_max_max_1 < adj_max_max_2 ) then
    result = - 4
    return
  else if ( adj_max_max_2 < adj_max_max_1 ) then
    result = + 4
    return
  end if
!
!  Test 5: Compare the adjacency max sequences.
!
  call mg_adj_max_seq ( adj1, nnode1, seq1 )

  call mg_adj_max_seq ( adj2, nnode2, seq2 )

  call i4vec_compare ( nnode1, seq1, seq2, result )

  if ( result < 0 ) then
    result = - 5
    return
  else if ( 0 < result ) then
    result = + 5
    return
  end if
!
!  Test 6: Compare the adjacency sequences.
!
  call mg_adj_seq ( adj1, nnode1, code1 )

  call mg_adj_seq ( adj2, nnode2, code2 )

  call i4mat_row_compare ( nnode1, nnode1, code1, code2, result )

  if ( result < 0 ) then
    result = - 6
    return
  else if ( 0 < result ) then
    result = + 6
    return
  end if
!
!  Test 7: Compare the codes.
!
  call mg_code_back ( adj1, nnode1, code1, order1 )

  call mg_code_back ( adj2, nnode2, code2, order2 )

  call mg_code_compare ( code1, code2, nnode1, nnode1, result )

  if ( result < 0 ) then
    result = - 7
    return
  else if ( 0 < result ) then
    result = + 7
    return
  end if

  result = 0

  return
end
subroutine mg_degree ( adj, nnode, degree )

!*****************************************************************************80
!
!! MG_DEGREE computes the degree of each node of a multigraph.
!
!  Discussion:
!
!    The degree of a node is the number of edges that are incident on it.
!    The sum of the degrees of the nodes is twice the number of edges.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,J) is the numbe of edges from node I to node J.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) DEGREE(NNODE), the degree of the nodes.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) degree(nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  degree(1:nnode) = 0

  do i = 1, nnode
    do j = 1, nnode
      if ( adj(i,j) /= 0 ) then
        degree(i) = degree(i) + adj(i,j)
      end if
    end do
  end do

  return
end
subroutine mg_degree_max ( adj, nnode, degree_max )

!*****************************************************************************80
!
!! MG_DEGREE_MAX computes the maximum node degree of a multigraph.
!
!  Discussion:
!
!    The maximum node degree of a multigraph is the maximum value of the
!    degree of the nodes.
!
!    If two multigraphs are isomorphic, they must have the same
!    maximum node degree.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,J) is the number of edges from node I to node J.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) DEGREE_MAX, the maximum node degree.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) degree
  integer ( kind = 4 ) degree_max
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  degree_max = 0

  do i = 1, nnode
    degree = 0
    do j = 1, nnode
      if ( adj(i,j) /= 0 ) then
        degree = degree + adj(i,j)
      end if
    end do
    degree_max = max ( degree_max, degree )
  end do

  return
end
subroutine mg_degree_seq ( adj, nnode, seq )

!*****************************************************************************80
!
!! MG_DEGREE_SEQ computes the degree sequence of a multigraph.
!
!  Discussion:
!
!    The degree sequence of a multigraph is constructed by computing the
!    degree of each node, and then ordering these values in decreasing order.
!
!    If two multigraphs are isomorphic, they must have the same degree sequence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,J) is the number of edges from node I to node J.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) SEQ(NNODE), the degree sequence.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) seq(nnode)

  do i = 1, nnode
    seq(i) = sum ( adj(i,1:nnode) )
  end do

  call i4vec_sort_heap_d ( nnode, seq )

  return
end
subroutine mg_edge_count ( adj, nnode, nedge )

!*****************************************************************************80
!
!! MG_EDGE_COUNT counts the number of edges in a multigraph.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,J) is the number of edges from node I to node J.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) NEDGE, the number of edges.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nedge

  nedge = 0

  do i = 1, nnode
    do j = 1, nnode

      if ( i == j ) then
        nedge = nedge + 2 * adj(i,j)
      else
        nedge = nedge + adj(i,j)
      end if

    end do
  end do

  nedge = nedge / 2

  return
end
subroutine mg_example_octo ( example, adj, nnode, seed )

!*****************************************************************************80
!
!! MG_EXAMPLE_OCTO sets up an 8 node example multigraph.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) EXAMPLE, chooses the example, and should be
!    between 1 and 7.
!
!    Output, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,J) is the number of edges from node I to node J.
!
!    Output, integer ( kind = 4 ) NNODE, the number of nodes, which should be 8.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
  implicit none

  integer ( kind = 4 ) adj(8,8)
  integer ( kind = 4 ) example
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nnode
  integer ( kind = 4 ) nsave
  integer ( kind = 4 ) seed

  if ( example <= 0 ) then
    nsave = i4_uniform ( 1, 7, seed )
  else
    nsave = mod ( example - 1, 7 ) + 1
  end if

  nnode = 8

  adj(1:nnode,1:nnode) = 0
!
!  The "basic" graph.
!
  if ( nsave == 1 ) then

    adj(1,2) = 1
    adj(1,6) = 1
    adj(1,8) = 1
    adj(2,3) = 1
    adj(2,5) = 2
    adj(3,4) = 1
    adj(3,8) = 3
    adj(4,5) = 1
    adj(4,7) = 4
    adj(5,6) = 1
    adj(6,7) = 1
    adj(7,8) = 1
!
!  Same NNODE, different number of edges.
!
  else if ( nsave == 2 ) then

    adj(1,2) = 1
    adj(1,6) = 1
    adj(1,8) = 1
    adj(2,3) = 1
    adj(2,8) = 2
    adj(3,4) = 1
    adj(3,5) = 3
    adj(4,5) = 1
    adj(4,7) = 3
    adj(5,6) = 1
    adj(6,7) = 1
    adj(7,8) = 1
!
!  Same NNODE, NEDGE, different degree sequence.
!
  else if ( nsave == 3 ) then

    adj(1,2) = 1
    adj(1,5) = 2
    adj(1,8) = 1
    adj(2,3) = 1
    adj(2,6) = 2
    adj(3,4) = 1
    adj(3,7) = 3
    adj(4,5) = 1
    adj(4,8) = 3
    adj(5,6) = 1
    adj(6,7) = 1
    adj(7,8) = 1
!
!  Same NNODE, NEDGE, degree sequence, different ADJ_MAX_MAX.
!
  else if ( nsave == 4 ) then

    adj(1,2) = 1
    adj(1,7) = 1
    adj(1,8) = 1
    adj(2,3) = 1
    adj(2,5) = 1
    adj(2,8) = 1
    adj(3,4) = 1
    adj(3,7) = 1
    adj(3,8) = 2
    adj(4,5) = 2
    adj(4,6) = 1
    adj(4,7) = 2
    adj(5,6) = 1
    adj(6,7) = 1
    adj(7,8) = 1
!
!  Same NNODE, NEDGE, degree sequence, ADJ_MAX_MAX, different ADJ_MAX_SEQ.
!
  else if ( nsave == 5 ) then

    adj(1,2) = 1
    adj(1,6) = 1
    adj(1,8) = 1
    adj(2,3) = 1
    adj(2,5) = 2
    adj(3,4) = 2
    adj(3,8) = 2
    adj(3,4) = 2
    adj(4,7) = 4
    adj(5,6) = 1
    adj(5,8) = 1
    adj(6,7) = 1
    adj(7,8) = 1
!
!  Same NNODE, NEDGE, degree sequence, ADJ_MAX_MAX, ADJ_MAX_SEQ,
!  different ADJ_SEQ.
!
  else if ( nsave == 6 ) then

    adj(1,2) = 4
    adj(1,3) = 2
    adj(2,4) = 2
    adj(3,4) = 3
    adj(5,6) = 2
    adj(5,7) = 1
    adj(5,8) = 1
    adj(6,7) = 1
    adj(6,8) = 1
    adj(7,8) = 1
!
!  Same NNODE, NEDGE, degree sequence, ADJ_MAX_MAX, ADJ_MAX_SEQ,
!  ADJ_SEQ, different code.
!
  else if ( nsave == 7 ) then

    adj(1,2) = 1
    adj(1,4) = 1
    adj(1,6) = 1
    adj(2,3) = 1
    adj(2,5) = 2
    adj(3,4) = 1
    adj(3,8) = 3
    adj(4,7) = 4
    adj(5,6) = 1
    adj(5,8) = 1
    adj(6,7) = 1
    adj(7,8) = 1

  end if
!
!  Copy the upper triangle.
!
  do i = 2, nnode
    do j = 1, i-1
      adj(i,j) = adj(j,i)
    end do
  end do
!
!  Now permute the nodes.
!
  call i4mat_perm_random ( nnode, adj, seed )

  return
end
subroutine mg_order_code ( adj, nnode, npart, code, order )

!*****************************************************************************80
!
!! MG_ORDER_CODE returns the multigraph code for a specific node ordering.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.
!    ADJ(I,J) is the number of edges from node I to node J.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) NPART, the number of nodes to consider.
!    NPART should be between 1 and NNODE.
!
!    If NPART is NNODE, then the full code is returned.
!
!    If NPART is less than NNODE, then the code is computed as
!    though only NPART nodes existed, namely, those specified in the
!    first NPART entries of order.  This option is provided so that
!    the routine can compute the portion of a code determined
!    by an incomplete ordering of the nodes.
!
!    Output, integer ( kind = 4 ) CODE(NNODE,NNODE), the code for this ordering.
!
!    Input, integer ( kind = 4 ) ORDER(NNODE), the ordering of the nodes.  ORDER(1)
!    is the first node, and so on.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) code(nnode,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ni
  integer ( kind = 4 ) nj
  integer ( kind = 4 ) npart
  integer ( kind = 4 ) order(nnode)

  do i = 1, nnode

    if ( i <= npart ) then

      ni = order(i)

      if ( order(i) < 1 .or. nnode < order(i) ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'MG_ORDER_CODE - Fatal error!'
        write ( *, '(a)' ) '  ORDER is not a proper permutation.'
        stop
      end if

    else
      ni = 0
    end if

    do j = i + 1, nnode

      if ( j <= npart ) then

        nj = order(j)

        if ( order(j) < 1 .or. nnode < order(j) ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'MG_ORDER_CODE - Fatal error!'
          write ( *, '(a)' ) '  ORDER is not a proper permutation.'
          stop
        end if

      else
        nj = 0
      end if

      if ( ni == 0 .or. nj == 0 ) then

        code(i,j) = 0
        code(j,i) = 0

      else if ( ( ni < nj .and. adj(ni,nj) /= 0 ) .or. &
                ( nj < ni .and. adj(nj,ni) /= 0 ) ) then

        code(i,j) = adj(ni,nj)
        code(j,i) = adj(nj,ni)

      else

        code(i,j) = 0
        code(j,i) = 0

      end if

    end do
  end do

  return
end
subroutine mg_print ( adj, nnode, title )

!*****************************************************************************80
!
!! MG_PRINT prints out an adjacency matrix for a multigraph.
!
!  Discussion:
!
!    Values between 0 and 9 will be printed as is.  Other values will
!    be printed as '*'.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency matrix.
!    ADJ(I,J) is the number of edges from node I to node J.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, character ( len = * ) TITLE, a title to be printed.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  character ( len = 80 ) string
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  do i = 1, nnode

    jhi = min ( nnode, 80 )

    do j = 1, jhi

      if ( 0 <= adj(i,j) .and. adj(i,j) <= 9 ) then
        string(j:j) = char ( 48 + adj(i,j) )
      else
        string(j:j) = '*'
      end if

    end do

    write ( *, '(2x,a)' ) string(1:jhi)

  end do

  return
end
subroutine mg_random ( adj, nnode, nedge, seed )

!*****************************************************************************80
!
!! MG_RANDOM generates a random multigraph on NNODE nodes with NEDGE edges.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency matrix.
!    ADJ(I,J) is the number of edges from node I to node J.
!    ADJ(I,I) will always be 0.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) NEDGE, the number of edges.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
  implicit none

  integer ( kind = 4 ) nnode
  integer ( kind = 4 ) nedge

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
!
!  Initialize the adjacency matrix.
!
  adj(1:nnode,1:nnode) = 0
!
!  Essentially, flip a coin NEDGE times to decide where each edge goes.
!
  do k = 1, nedge

    i = i4_uniform ( 1, nnode, seed )
    j = i4_uniform ( 1, nnode-1, seed )
    if ( i <= j ) then
      j = j + 1
    end if

    if ( i < 1 .or. nnode < i ) then
      write ( *, '(a,i8)' ) 'I = ', i
      stop
    end if

    if ( j < 1 .or. nnode < j ) then
      write ( *, '(a,i8)' ) 'J = ', j
      stop
    end if

    adj(i,j) = adj(i,j) + 1
    adj(j,i) = adj(j,i) + 1

  end do

  return
end
subroutine node_order_print ( nnode, order, title )

!*****************************************************************************80
!
!! NODE_ORDER_PRINT prints out a node ordering.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) ORDER(NNODE), the node ordering.  ORDER(1) is the label
!    of the node which is to be taken as the first node, and so on.
!
!    Input, character ( len = * ) TITLE, a title to be printed.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) order(nnode)
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  inc = 15

  do ilo = 1, nnode, inc

    ihi = min ( ilo + inc - 1, nnode )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,a,2x,15i4)' ) 'Order:', ( i, i = ilo, ihi )
    write ( *, '(2x,a,2x,15i4)' ) 'Label:', order(ilo:ihi)

  end do

  return
end
subroutine perm_check ( n, p, ierror )

!*****************************************************************************80
!
!! PERM_CHECK checks that a vector represents a permutation.
!
!  Discussion:
!
!    The routine verifies that each of the integers from 1
!    to N occurs among the N entries of the permutation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries.
!
!    Input, integer ( kind = 4 ) P(N), the array to check.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, the array represents a permutation.
!    nonzero, the array does not represent a permutation.  The smallest
!    missing value is equal to IERROR.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ifind
  integer ( kind = 4 ) iseek
  integer ( kind = 4 ) p(n)

  ierror = 0

  do iseek = 1, n

    ierror = iseek

    do ifind = 1, n
      if ( p(ifind) == iseek ) then
        ierror = 0
        exit
      end if
    end do

    if ( ierror /= 0 ) then
      return
    end if

  end do

  return
end
subroutine perm_cycle ( n, p, isgn, ncycle, iopt )

!*****************************************************************************80
!
!! PERM_CYCLE analyzes a permutation.
!
!  Discussion:
!
!    The routine will count cycles, find the sign of a permutation,
!    and tag a permutation.
!
!  Example:
!
!    Input:
!
!      N = 9
!      IOPT = 1
!      P = 2, 3, 9, 6, 7, 8, 5, 4, 1
!
!    Output:
!
!      NCYCLE = 3
!      ISGN = +1
!      P = -2, 3, 9, -6, -7, 8, 5, 4, 1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of objects being permuted.
!
!    Input/output, integer ( kind = 4 ) P(N).  On input, P describes a
!    permutation, in the sense that entry I is to be moved to P(I).
!    If IOPT = 0, then P will not be changed by this routine.
!    If IOPT = 1, then on output, P will be "tagged".  That is,
!    one element of every cycle in P will be negated.  In this way,
!    a user can traverse a cycle by starting at any entry I1 of P
!    which is negative, moving to I2 = ABS(P(I1)), then to
!    P(I2), and so on, until returning to I1.
!
!    Output, integer ( kind = 4 ) ISGN, the "sign" of the permutation, which is
!    +1 if the permutation is even, -1 if odd.  Every permutation
!    may be produced by a certain number of pairwise switches.
!    If the number of switches is even, the permutation itself is
!    called even.
!
!    Output, integer ( kind = 4 ) NCYCLE, the number of cycles in the permutation.
!
!    Input, integer ( kind = 4 ) IOPT, requests tagging.
!    0, the permutation will not be tagged.
!    1, the permutation will be tagged.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) iopt
  integer ( kind = 4 ) is
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) ncycle
  integer ( kind = 4 ) p(n)

  is = 1
  ncycle = n

  do i = 1, n

    i1 = p(i)

    do while ( i < i1 )
      ncycle = ncycle - 1
      i2 = p(i1)
      p(i1) = -i2
      i1 = i2
    end do

    if ( iopt /= 0 ) then
      is = - isign ( 1, p(i) )
    end if

    p(i) = isign ( p(i), is )

  end do

  isgn = 1 - 2 * mod ( n-ncycle, 2 )

  return
end
subroutine perm_free ( ipart, npart, ifree, nfree )

!*****************************************************************************80
!
!! PERM_FREE reports the number of unused items in a partial permutation.
!
!  Discussion:
!
!    It is assumed that the N objects being permuted are the integers
!    from 1 to N, and that IPART contains a "partial" permutation, that
!    is, the NPART entries of IPART represent the beginning of a
!    permutation of all N items.
!
!    The routine returns in IFREE the items that have not been used yet.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IPART(NPART), the partial permutation, which should
!    contain, at most once, some of the integers between 1 and
!    NPART+NFREE.
!
!    Input, integer ( kind = 4 ) NPART, the number of entries in IPART.  NPART may be 0.
!
!    Output, integer ( kind = 4 ) IFREE(NFREE), the integers between 1 and NPART+NFREE
!    that were not used in IPART.
!
!    Input, integer ( kind = 4 ) NFREE, the number of integers that have not been
!    used in IPART.  This is simply N - NPART.  NFREE may be zero.
!
  implicit none

  integer ( kind = 4 ) nfree
  integer ( kind = 4 ) npart

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ifree(nfree)
  integer ( kind = 4 ) ipart(npart)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) match
  integer ( kind = 4 ) n

  n = npart + nfree

  if ( npart < 0 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PERM_FREE - Fatal error!'
    write ( *, '(a)' ) '  NPART < 0.'
    write ( *, '(a,i8)' ) '  NPART = ', npart
    stop

  else if ( npart == 0 ) then

    call i4vec_indicator ( n, ifree )

  else if ( nfree < 0 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PERM_FREE - Fatal error!'
    write ( *, '(a)' ) '  NFREE < 0.'
    write ( *, '(a,i8)' ) '  NFREE = ', nfree
    stop

  else if ( nfree == 0 ) then

    return

  else

    k = 0

    do i = 1, n

      match = 0

      do j = 1, npart
        if ( ipart(j) == i ) then
          match = j
          exit
        end if
      end do

      if ( match == 0 ) then

        k = k + 1

        if ( nfree < k ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)'    ) 'PERM_FREE - Fatal error!'
          write ( *, '(a)'    ) '  The partial permutation is illegal.'
          write ( *, '(a)'    ) '  It should contain, at most once, some of'
          write ( *, '(a,i8)' ) '  the integers between 1 and N = ', n
          write ( *, '(a)'    ) '  The number of integers that have not'
          write ( *, '(a,i8)' ) '  been used is at least K = ', k
          write ( *, '(a,i8)' ) '  This value should be exactly NFREE = ', &
            nfree
          call i4vec_print ( npart, ipart, '  The partial permutation:' )
          stop
        end if

        ifree(k) = i

      end if

    end do

  end if

  return
end
subroutine perm_next ( n, p, more, even )

!*****************************************************************************80
!
!! PERM_NEXT computes all of the permutations on N objects, one at a time.
!
!  Discussion:
!
!    Note that if this routine is called with MORE = .TRUE., any
!    permutation in P, and EVEN = .TRUE., then the successor of the input
!    permutation will be produced, unless P is the last permutation
!    on N letters, in which case P(1) will be set to 0 on return.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 March 2001
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of objects being permuted.
!
!    Input/output, integer ( kind = 4 ) P(N).
!    On input, P contains the previous permutation.
!    On output, P contains the next permutation.
!
!    Input/output, logical MORE.
!    On input, MORE = FALSE means this is the first call.
!    On output, MORE = FALSE means there are no more permutations.
!
!    Output, logical EVEN, is TRUE if the output permutation is even.
!
  implicit none

  integer ( kind = 4 ) n

  logical even
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) ia
  integer ( kind = 4 ) id
  integer ( kind = 4 ) is
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  logical more
  integer ( kind = 4 ) p(n)

  if ( .not. more ) then

    call i4vec_indicator ( n, p )

    more = .true.
    even = .true.

    if ( n == 1 ) then
      more = .false.
      return
    end if

    if ( p(n) /= 1 .or. p(1) /= 2 + mod ( n, 2 ) ) then
      return
    end if

    do i = 1, n-3
      if ( p(i+1) /= p(i) + 1 ) then
        return
      end if
    end do

    more = .false.

  else

    if ( n == 1 ) then
      p(1) = 0
      more = .false.
      return
    end if

    if ( even ) then

      ia = p(1)
      p(1) = p(2)
      p(2) = ia
      even = .false.

      if ( p(n) /= 1 .or. p(1) /= 2 + mod ( n, 2 ) ) then
        return
      end if

      do i = 1, n-3
        if ( p(i+1) /= p(i) + 1 ) then
          return
        end if
      end do

      more = .false.
      return

    else

      is = 0
      more = .false.

      do i1 = 2, n

        ia = p(i1)
        i = i1 - 1
        id = 0

        do j = 1, i
          if ( ia < p(j) ) then
            id = id + 1
          end if
        end do

        is = id + is

        if ( id /= i * mod ( is, 2 ) ) then
          more = .true.
          exit
        end if

      end do

      if ( .not. more ) then
        p(1) = 0
        return
      end if

    end if

    m = mod ( is+1, 2 ) * ( n + 1 )

    do j = 1, i

      if ( sign ( 1, p(j)-ia ) /= sign ( 1, p(j)-m ) ) then
        m = p(j)
        l = j
      end if

    end do

    p(l) = ia
    p(i1) = m
    even = .true.

  end if

  return
end
subroutine perm_random ( n, a, seed )

!*****************************************************************************80
!
!! PERM_RANDOM selects a random permutation of N objects.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of objects to be permuted.
!
!    Output, integer ( kind = 4 ) A(N), the random permutation.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) j
  integer ( kind = 4 ) seed

  call i4vec_indicator ( n, a )

  do i = 1, n
    j = i4_uniform ( i, n, seed )
    call i4_swap ( a(j), a(i) )
  end do

  return
end
subroutine pruefer_to_tree_arc ( nnode, a, inode, jnode )

!*****************************************************************************80
!
!! PRUEFER_TO_TREE_ARC is given a Pruefer code, and computes the tree.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 October 1999
!
!  Author:
!
!    Original FORTRAN77 version by Dennis Stanton, Dennis White.
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Dennis Stanton, Dennis White,
!    Constructive Combinatorics,
!    Springer Verlag, New York, 1986.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes in the tree.
!
!    Input, integer ( kind = 4 ) A(NNODE-2), the Pruefer code of the tree.
!
!    Output, integer ( kind = 4 ) INODE(NNODE-1), JNODE(NNODE-1), the edge array of the
!    tree.  The I-th edge joins nodes INODE(I) and JNODE(I).
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) a(nnode-2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) inode(nnode-1)
  integer ( kind = 4 ) iwork(nnode)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jnode(nnode-1)
!
!  Initialize IWORK(I) to count the number of neighbors of node I.
!  The Pruefer code uses each node one less time than its total
!  number of neighbors.
!
  iwork(1:nnode) = 1

  do i = 1, nnode-2
    iwork(a(i)) = iwork(a(i)) + 1
  end do
!
!  Now process each entry in the Pruefer code.
!
  do i = 1, nnode-2

    ii = 0
    do j = 1, nnode
      if ( iwork(j) == 1 ) then
        ii = j
      end if
    end do

    inode(i) = ii
    jnode(i) = a(i)
    iwork(ii) = 0
    iwork(a(i)) = iwork(a(i)) - 1

  end do

  inode(nnode-1) = a(nnode-2)

  if ( a(nnode-2) /= 1 ) then
    jnode(nnode-1) = 1
  else
    jnode(nnode-1) = 2
  end if

  return
end
function r4_uniform_01 ( seed )

!*****************************************************************************80
!
!! R4_UNIFORM_01 returns a unit pseudorandom R4.
!
!  Discussion:
!
!    An R4 is a real ( kind = 4 ) value.
!
!    This routine implements the recursion
!
!      seed = 16807 * seed mod ( 2**31 - 1 )
!      r4_uniform_01 = seed / ( 2**31 - 1 )
!
!    The integer arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!    If the initial seed is 12345, then the first three computations are
!
!      Input     Output      R4_UNIFORM_01
!      SEED      SEED
!
!         12345   207482415  0.096616
!     207482415  1790989824  0.833995
!    1790989824  2035175616  0.947702
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley Interscience, page 95, 1998.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 4 ) R4_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real    ( kind = 4 ) r4_uniform_01

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + 2147483647
  end if

  r4_uniform_01 = real ( seed, kind = 4 ) * 4.656612875E-10

  return
end
subroutine r8mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8MAT_PRINT prints an R8MAT.
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
!    Input, character ( len = * ) TITLE, a title to be printed.
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
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 March 2005
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
!    Input, character ( len = * ) TITLE, an optional title.
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

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i7,7x)') j
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

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j), j = 1, inc )

    end do

  end do

  write ( *, '(a)' ) ' '

  return
end
function r8_normal_01 ( seed )

!*****************************************************************************80
!
!! R8_NORMAL_01 returns a unit pseudonormal R8.
!
!  Discussion:
!
!    The standard normal probability distribution function (PDF) has
!    mean 0 and standard deviation 1.
!
!    Because this routine uses the Box Muller method, it requires pairs
!    of uniform random values to generate a pair of normal random values.
!    This means that on every other call, essentially, the input value of
!    SEED is ignored, since the code saves the second normal random value.
!
!    If you didn't know this, you might be confused since, usually, the
!    output of a random number generator can be completely controlled by
!    the input value of the SEED.  If I were more careful, I could rewrite
!    this routine so that it would distinguish between cases where the input
!    value of SEED is the output value from the previous call (all is well)
!    and those cases where it is not (the user has decided to do something
!    new.  Restart the uniform random number sequence.)  But I'll leave
!    that for later.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
!    Output, real ( kind = 8 ) R8_NORMAL_01, a sample of the standard normal PDF.
!
  implicit none

  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) r8_normal_01
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  integer ( kind = 4 ), save :: seed2 = 0
  integer ( kind = 4 ), save :: used = 0
  real ( kind = 8 ) x
  real ( kind = 8 ), save :: y = 0.0D+00
!
!  On odd numbered calls, generate two uniforms, create two normals,
!  return the first normal and its corresponding seed.
!
  if ( mod ( used, 2 ) == 0 ) then

    r1 = r8_uniform_01 ( seed )

    if ( r1 == 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8_NORMAL_01 - Fatal error!'
      write ( *, '(a)' ) '  R8_UNIFORM_01 returned a value of 0.'
      stop
    end if

    seed2 = seed
    r2 = r8_uniform_01 ( seed2 )

    x = sqrt ( -2.0D+00 * log ( r1 ) ) * cos ( 2.0D+00 * pi * r2 )
    y = sqrt ( -2.0D+00 * log ( r1 ) ) * sin ( 2.0D+00 * pi * r2 )
!
!  On odd calls, return the second normal and its corresponding seed.
!
  else

    seed = seed2
    x = y

  end if

  used = used + 1

  r8_normal_01 = x

  return
end
function r8_uniform_01 ( seed )

!*****************************************************************************80
!
!! R8_UNIFORM_01 returns a unit pseudorandom R8.
!
!  Discussion:
!
!    An R8 is a real ( kind = 8 ) value.
!
!    For now, the input quantity SEED is an integer ( kind = 4 ) variable.
!
!    This routine implements the recursion
!
!      seed = 16807 * seed mod ( 2**31 - 1 )
!      r8_uniform_01 = seed / ( 2**31 - 1 )
!
!    The integer arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!    If the initial seed is 12345, then the first three computations are
!
!      Input     Output      R8_UNIFORM_01
!      SEED      SEED
!
!         12345   207482415  0.096616
!     207482415  1790989824  0.833995
!    1790989824  2035175616  0.947702
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley Interscience, page 95, 1998.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which should
!    NOT be 0. On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R8_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer ( kind = 4 ) k
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + 2147483647
  end if
!
!  Although SEED can be represented exactly as a 32 bit integer,
!  it generally cannot be represented exactly as a 32 bit real number!
!
  r8_uniform_01 = real ( seed, kind = 8 ) * 4.656612875D-10

  return
end
subroutine sort_heap_external ( n, indx, i, j, isgn )

!*****************************************************************************80
!
!! SORT_HEAP_EXTERNAL externally sorts a list of items into ascending order.
!
!  Discussion:
!
!    The actual list of data is not passed to the routine.  Hence this
!    routine may be used to sort integers, reals, numbers, names,
!    dates, shoe sizes, and so on.  After each call, the routine asks
!    the user to compare or interchange two items, until a special
!    return value signals that the sorting is completed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 February 2004
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items to be sorted.
!
!    Input/output, integer ( kind = 4 ) INDX, the main communication signal.
!
!    The user must set INDX to 0 before the first call.
!    Thereafter, the user should not change the value of INDX until
!    the sorting is done.
!
!    On return, if INDX is
!
!      greater than 0,
!      * interchange items I and J;
!      * call again.
!
!      less than 0,
!      * compare items I and J;
!      * set ISGN = -1 if I < J, ISGN = +1 if J < I;
!      * call again.
!
!      equal to 0, the sorting is done.
!
!    Output, integer ( kind = 4 ) I, J, the indices of two items.
!    On return with INDX positive, elements I and J should be interchanged.
!    On return with INDX negative, elements I and J should be compared, and
!    the result reported in ISGN on the next call.
!
!    Input, integer ( kind = 4 ) ISGN, results of comparison of elements I and J.
!    (Used only when the previous call returned INDX less than 0).
!    ISGN <= 0 means I is less than or equal to J;
!    0 <= ISGN means I is greater than or equal to J.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ), save :: i_save = 0
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  integer ( kind = 4 ), save :: j_save = 0
  integer ( kind = 4 ), save :: k = 0
  integer ( kind = 4 ), save :: k1 = 0
  integer ( kind = 4 ) n
  integer ( kind = 4 ), save :: n1 = 0
!
!  INDX = 0: This is the first call.
!
  if ( indx == 0 ) then

    i_save = 0
    j_save = 0
    k = n / 2
    k1 = k
    n1 = n
!
!  INDX < 0: The user is returning the results of a comparison.
!
  else if ( indx < 0 ) then

    if ( indx == -2 ) then

      if ( isgn < 0 ) then
        i_save = i_save + 1
      end if

      j_save = k1
      k1 = i_save
      indx = -1
      i = i_save
      j = j_save
      return

    end if

    if ( 0 < isgn ) then
      indx = 2
      i = i_save
      j = j_save
      return
    end if

    if ( k <= 1 ) then

      if ( n1 == 1 ) then
        i_save = 0
        j_save = 0
        indx = 0
      else
        i_save = n1
        n1 = n1 - 1
        j_save = 1
        indx = 1
      end if

      i = i_save
      j = j_save
      return

    end if

    k = k - 1
    k1 = k
!
!  0 < INDX, the user was asked to make an interchange.
!
  else if ( indx == 1 ) then

    k1 = k

  end if

  do

    i_save = 2 * k1

    if ( i_save == n1 ) then
      j_save = k1
      k1 = i_save
      indx = -1
      i = i_save
      j = j_save
      return
    else if ( i_save <= n1 ) then
      j_save = i_save + 1
      indx = -2
      i = i_save
      j = j_save
      return
    end if

    if ( k <= 1 ) then
      exit
    end if

    k = k - 1
    k1 = k

  end do

  if ( n1 == 1 ) then
    i_save = 0
    j_save = 0
    indx = 0
    i = i_save
    j = j_save
  else
    i_save = n1
    n1 = n1 - 1
    j_save = 1
    indx = 1
    i = i_save
    j = j_save
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
subroutine tree_arc_random ( nnode, code, inode, jnode, seed )

!*****************************************************************************80
!
!! TREE_ARC_RANDOM selects a random labeled tree and its Pruefer code.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes in the trees.
!
!    Output, integer ( kind = 4 ) CODE(NNODE-2), the Pruefer code for the labeled tree.
!
!    Output, integer ( kind = 4 ) INODE(NNODE-1), JNODE(NNODE-1), the edge array for
!    the tree.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) code(nnode-2)
  integer ( kind = 4 ) inode(nnode-1)
  integer ( kind = 4 ) jnode(nnode-1)
  integer ( kind = 4 ) seed

  call i4vec_uniform ( nnode-2, 0, nnode-1, seed, code )

  code(1:nnode-2) = code(1:nnode-2) + 1

  call pruefer_to_tree_arc ( nnode, code, inode, jnode )

  return
end
subroutine wg_random ( nnode, seed, weight )

!*****************************************************************************80
!
!! WG_RANDOM generates a random weighted graph.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
!    Output, real ( kind = 8 ) WEIGHT(NNODE,NNODE), the symmetric
!    weighted graph.  The diagonal entries are zero.
!
  implicit none

  integer ( kind = 4 ) nnode

  real ( kind = 8 ) r8_normal_01
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) seed
  real ( kind = 8 ) weight(nnode,nnode)

  do i = 1, nnode
    weight(i,i) = 0.0D+00
    do j = i+1, nnode
      weight(i,j) = r8_normal_01 ( seed )
      weight(j,i) = weight(i,j)
    end do
  end do

  return
end
