subroutine balanc ( nm, n, a, low, igh, scale )

!*****************************************************************************80
!
!! BALANC balances a real matrix before eigenvalue calculations.
!
!  Discussion:
!
!    This subroutine balances a real matrix and isolates eigenvalues
!    whenever possible.
!
!    Suppose that the principal submatrix in rows LOW through IGH
!    has been balanced, that P(J) denotes the index interchanged
!    with J during the permutation step, and that the elements
!    of the diagonal matrix used are denoted by D(I,J).  Then
!
!      SCALE(J) = P(J),    J = 1,...,LOW-1,
!               = D(J,J),  J = LOW,...,IGH,
!               = P(J)     J = IGH+1,...,N.
!
!    The order in which the interchanges are made is N to IGH+1,
!    then 1 to LOW-1.
!
!    Note that 1 is returned for LOW if IGH is zero formally.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 December 2008
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe, Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow, 
!    Y Ikebe, V Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NM, the leading dimension of A, which must
!    be at least N.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input/output, real ( kind = 8 ) A(NM,N), the N by N matrix.  On output,
!    the matrix has been balanced.
!
!    Output, integer ( kind = 4 ) LOW, IGH, indicate that A(I,J) is equal to zero if
!    (1) I is greater than J and
!    (2) J=1,...,LOW-1 or I=IGH+1,...,N.
!
!    Output, real ( kind = 8 ) SCALE(N), contains information determining the
!    permutations and scaling factors used.
!
  implicit none

  integer ( kind = 4 ) nm
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(nm,n)
  real ( kind = 8 ) b2
  real ( kind = 8 ) c
  real ( kind = 8 ) f
  real ( kind = 8 ) g
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iexc
  integer ( kind = 4 ) igh
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) low
  integer ( kind = 4 ) m
  logical noconv
  real ( kind = 8 ) r
  real ( kind = 8 ), parameter :: radix = 16.0D+00
  real ( kind = 8 ) s
  real ( kind = 8 ) scale(n)

  iexc = 0
  j = 0
  m = 0

  b2 = radix * radix
  k = 1
  l = n
  go to 100

20 continue

  scale(m) = j

  if ( j /= m ) then

    do i = 1, l
      call r8_swap ( a(i,j), a(i,m) )
    end do

    do i = k, n
      call r8_swap ( a(j,i), a(m,i) )
    end do

  end if

50 continue

  if ( iexc == 2 ) then
    go to 130
  end if
!
!  Search for rows isolating an eigenvalue and push them down.
!
80 continue

  if ( l == 1 ) then
    low = k
    igh = l
    return
  end if

  l = l - 1

100 continue

  do j = l, 1, -1

     do i = 1, l
       if ( i /= j ) then
         if ( a(j,i) /= 0.0D+00 ) then
           go to 120
         end if
       end if
     end do

     m = l
     iexc = 1
     go to 20

120  continue

  end do

  go to 140
!
!  Search for columns isolating an eigenvalue and push them left.
!
130 continue

  k = k + 1

140 continue

  do j = k, l

    do i = k, l
      if ( i /= j ) then
        if ( a(i,j) /= 0.0D+00 ) then
          go to 170
        end if
      end if
    end do

    m = k
    iexc = 2
    go to 20

170 continue

  end do
!
!  Balance the submatrix in rows K to L.
!
  scale(k:l) = 1.0D+00
!
!  Iterative loop for norm reduction.
!
  noconv = .true.

  do while ( noconv )

    noconv = .false.

    do i = k, l

      c = 0.0D+00
      r = 0.0D+00

      do j = k, l
        if ( j /= i ) then
          c = c + abs ( a(j,i) )
          r = r + abs ( a(i,j) )
        end if
      end do
!
!  Guard against zero C or R due to underflow.
!
      if ( c /= 0.0D+00 .and. r /= 0.0D+00 ) then

        g = r / radix
        f = 1.0D+00
        s = c + r

        do while ( c < g )
          f = f * radix
          c = c * b2
        end do

        g = r * radix

        do while ( g <= c )
          f = f / radix
          c = c / b2
        end do
!
!  Balance.
!
        if ( ( c + r ) / f < 0.95D+00 * s ) then

          g = 1.0D+00 / f
          scale(i) = scale(i) * f
          noconv = .true.

          a(i,k:n) = a(i,k:n) * g
          a(1:l,i) = a(1:l,i) * f

        end if

      end if

    end do

  end do

  low = k
  igh = l

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

  character              c
  integer   ( kind = 4 ) itemp

  itemp = ichar ( c )

  if ( 97 <= itemp .and. itemp <= 122 ) then
    c = char ( itemp - 32 )
  end if

  return
end
function ch_eqi ( c1, c2 )

!*****************************************************************************80
!
!! CH_EQI is a case insensitive comparison of two characters for equality.
!
!  Example:
!
!    C_EQI ( 'A', 'a' ) is .TRUE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C1, C2, the characters to compare.
!
!    Output, logical CH_EQI, the result of the comparison.
!
  implicit none

  logical ch_eqi
  character c1
  character c2
  character cc1
  character cc2

  cc1 = c1
  cc2 = c2

  call ch_cap ( cc1 )
  call ch_cap ( cc2 )

  if ( cc1 == cc2 ) then
    ch_eqi = .true.
  else
    ch_eqi = .false.
  end if

  return
end
subroutine ch_to_digit ( c, digit )

!*****************************************************************************80
!
!! CH_TO_DIGIT returns the integer value of a base 10 digit.
!
!  Example:
!
!     C   DIGIT
!    ---  -----
!    '0'    0
!    '1'    1
!    ...  ...
!    '9'    9
!    ' '    0
!    'X'   -1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C, the decimal digit, '0' through '9' or blank
!    are legal.
!
!    Output, integer ( kind = 4 ) DIGIT, the corresponding integer value.  
!    If C was 'illegal', then DIGIT is -1.
!
  implicit none

  character c
  integer ( kind = 4 ) digit

  if ( lge ( c, '0' ) .and. lle ( c, '9' ) ) then

    digit = ichar ( c ) - 48

  else if ( c == ' ' ) then

    digit = 0

  else

    digit = - 1

  end if

  return
end
subroutine catalan ( n, c )

!*****************************************************************************80
!
!! CATALAN computes the Catalan numbers, from C(0) to C(N).
!
!  First values:
!
!     C(0)     1
!     C(1)     1
!     C(2)     2
!     C(3)     5
!     C(4)    14
!     C(5)    42
!     C(6)   132
!     C(7)   429
!     C(8)  1430
!     C(9)  4862
!    C(10) 16796
!
!  Formula:
!
!    C(N) = (2*N)! / ( (N+1) * (N!) * (N!) )
!         = 1 / (N+1) * COMB ( 2N, N )
!         = 1 / (2N+1) * COMB ( 2N+1, N+1).
!
!  Recursion:
!
!    C(N) = 2 * (2*N-1) * C(N-1) / (N+1)
!    C(N) = SUM ( I = 1 to N-1 ) C(I) * C(N-I)
!
!  Comments:
!
!    The Catalan number C(N) counts:
!
!    1) the number of binary trees on N vertices;
!    2) the number of ordered trees on N+1 vertices;
!    3) the number of full binary trees on 2N+1 vertices;
!    4) the number of well formed sequences of 2N parentheses;
!    5) number of ways 2N ballots can be counted, in order,
!       with N positive and N negative, so that the running sum
!       is never negative;
!    6) the number of standard tableaus in a 2 by N rectangular Ferrers diagram;
!    7) the number of monotone functions from [1..N} to [1..N} which
!       satisfy f(i) <= i for all i,
!    8) the number of ways to triangulate a polygon with N+2 vertices.
!
!  Example:
!
!    N = 3
!
!    ()()()
!    ()(())
!    (()())
!    (())()
!    ((()))
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 August 1998
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Dennis Stanton, Dennis White,
!    Constructive Combinatorics,
!    Springer Verlag, New York, 1986.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of Catalan numbers desired.
!
!    Output, integer ( kind = 4 ) C(0:N), the Catalan numbers from C(0) to C(N).
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) c(0:n)

  c(0) = 1
!
!  The extra parentheses ensure that the integer division is
!  done AFTER the integer multiplication.
!
  do i = 1, n
    c(i) = ( c(i-1) * 2 * ( 2 * i - 1 ) ) / ( i + 1 )
  end do

  return
end
subroutine color_digraph_adj_degree ( adj, nnode, indegree, outdegree )

!*****************************************************************************80
!
!! COLOR_DIGRAPH_ADJ_DEGREE computes the indegree and outdegree of each node.
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
        outdegree(i) = outdegree(i) + adj(i,j)
        indegree(j) = indegree(j) + adj(i,j)
      end if
    end do
  end do

  return
end
subroutine color_digraph_adj_degree_seq ( adj, lda, nnode, in_seq, out_seq )

!*****************************************************************************80
!
!! COLOR_DIGRAPH_ADJ_DEGREE_SEQ computes the degree sequence of a color digraph.
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
!    Input, integer ( kind = 4 ) ADJ(LDA,NNODE), the adjacency information.  
!    ADJ(I,I) is the color of node I; otherwise, ADJ(I,J) is positive
!    if there is an edge from node I to node J.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the ADJ array,
!    which must be at least NNODE.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) IN_SEQ(NNODE), OUT_SEQ(NNODE),
!    the degree sequence of the digraph.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) in_seq(nnode)
  integer ( kind = 4 ) out_seq(nnode)

  call color_digraph_adj_degree ( adj, nnode, in_seq, out_seq )

  call i4vec2_sort_d ( nnode, out_seq, in_seq )

  return
end
subroutine color_digraph_adj_edge_count ( adj, lda, nnode, nedge )

!*****************************************************************************80
!
!! COLOR_DIGRAPH_ADJ_EDGE_COUNT counts the number of edges in a color digraph.
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
!    Input, integer ( kind = 4 ) ADJ(LDA,NNODE), the adjacency information.  
!    ADJ(I,I) is the color of node I; otherwise, ADJ(I,J) is positive
!    if there is an edge from node I to node J.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the ADJ array,
!    which must be at least NNODE.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) NEDGE, the number of edges.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(lda,nnode)
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
subroutine color_digraph_adj_example_cube ( adj, lda, nnode )

!*****************************************************************************80
!
!! COLOR_DIGRAPH_ADJ_EXAMPLE_CUBE sets up the cube color digraph.
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
!    Output, integer ( kind = 4 ) ADJ(LDA,NNODE), the adjacency information.  
!    ADJ(I,I) is the color of node I; otherwise, ADJ(I,J) is positive
!    if there is an edge from node I to node J.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the ADJ array,
!    which must be at least NNODE.
!
!    Output, integer ( kind = 4 ) NNODE, the number of nodes.
!
  implicit none

  integer ( kind = 4 ), parameter :: BLUE = 1
  integer ( kind = 4 ), parameter :: GREEN = 2
  integer ( kind = 4 ), parameter :: RED = 3

  integer ( kind = 4 ) lda

  integer ( kind = 4 ) adj(lda,lda)
  integer ( kind = 4 ) nnode

  nnode = 8

  if ( lda < nnode ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'COLOR_DIGRAPH_ADJ_EXAMPLE_CUBE - Fatal error!'
    write ( *, '(a)' ) '  LDA is too small.'
    stop
  end if

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
subroutine color_digraph_adj_example_octo ( lda, example, seed, nnode, adj )

!*****************************************************************************80
!
!! COLOR_DIGRAPH_ADJ_EXAMPLE_OCTO sets up an 8 node example color digraph.
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
!    There are 7 graphs to choose from.  They are all on 8 nodes.  The first
!    5 have degree 3 at every node.  Graphs 6 and 7 have degree 5 at every
!    node.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the ADJ array,
!    which must be at least NNODE.
!
!    Input, integer ( kind = 4 ) EXAMPLE, should be between 1 and 60, and indicates
!    which example graph to pick.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
!    Output, integer ( kind = 4 ) NNODE, the number of nodes, which should be 8.
!
!    Output, integer ( kind = 4 ) ADJ(LDA,LDA), the adjacency information.
!    ADJ(I,I) is the color of node I.
!    ADJ(I,J) is 1 if nodes I and J are adjacent and 0 otherwise.
!
  implicit none

  integer ( kind = 4 ), parameter :: BLUE = 1
  integer ( kind = 4 ), parameter :: GREEN = 2
  integer ( kind = 4 ), parameter :: RED = 3
  integer ( kind = 4 ), parameter :: YELLOW = 4

  integer ( kind = 4 ) lda

  integer ( kind = 4 ) adj(lda,lda)
  integer ( kind = 4 ) example
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) j
  integer ( kind = 4 ) msave
  integer ( kind = 4 ) nnode
  integer ( kind = 4 ) nsave
  integer ( kind = 4 ) seed

  if ( nnode <= 0 ) then
    nsave = i4_uniform ( 1, 12, seed )
    msave = i4_uniform ( 1, 5, seed )
  else
    example = mod ( example - 1, 60 ) + 1
    msave = ( example - 1 ) / 12 + 1
    nsave = mod ( example - 1, 12 ) + 1
  end if

  nnode = 8

  if ( lda < nnode ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'COLOR_DIGRAPH_ADJ_EXAMPLE_OCTO - Fatal error!'
    write ( *, '(a)' ) '  LDA is too small.'
    stop
  end if

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
  call i4mat_perm_random ( lda, nnode, seed, adj )

  return
end
subroutine color_digraph_adj_print ( adj, lda, nnode, title )

!*****************************************************************************80
!
!! COLOR_DIGRAPH_ADJ_PRINT prints out the adjacency matrix of a color digraph.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
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
!    Input, integer ( kind = 4 ) ADJ(LDA,NNODE), the adjacency information.  
!    ADJ(I,I) is the color of node I; otherwise, ADJ(I,J) is positive
!    if there is an edge from node I to node J.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of ADJ, which must be
!    least NNODE.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  character ( len = 80 ) string
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = 1, nnode

    do j = 1, nnode

      k = (j-1) * 3 + 1
      write ( string(k:k+2), '(i3)' ) adj(i,j)

    end do

    write ( *, '(i2,2x,a)' ) i, string(1:3*nnode)

  end do

  return
end
subroutine color_digraph_adj_random ( nnode, ncolor, nedge, seed, adj )

!*****************************************************************************80
!
!! COLOR_DIGRAPH_ADJ_RANDOM generates a random color graph.
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
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) NCOLOR, the number of colors available.  
!    Each node is assumed to have an associated color, between 1 and NCOLOR,
!    which will be chosen at random.
!
!    Input, integer ( kind = 4 ) NEDGE, the number of edges, which must be between
!    0 and NNODE*(NNODE-1).
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
!    Output, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency information.  
!    ADJ(I,I) is the color of node I; otherwise, ADJ(I,J) is positive
!    if there is an edge from node I to node J.
!
  implicit none

  integer ( kind = 4 ) ncolor
  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) color
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) iwork(nedge)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) maxedge
  integer ( kind = 4 ) perm(ncolor)
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) subset(ncolor)

  if ( nnode <= 0  ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'COLOR_DIGRAPH_ADJ_RANDOM - Fatal error!'
    write ( *, '(a,i8)' ) '  NNODE = ', nnode
    write ( *, '(a)' ) '  but NNODE must be at least 1.'
    stop
  end if

  maxedge = nnode * ( nnode - 1 )

  if ( nedge < 0 .or. maxedge < nedge ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'COLOR_DIGRAPH_ADJ_RANDOM - Fatal error!'
    write ( *, '(a,i8)' ) '  NEDGE = ', nedge
    write ( *, '(a)' ) '  but NEDGE must be at least 0, and '
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
  call ksub_random ( nnode, ncolor, seed, subset )

  call perm_random ( ncolor, seed, perm )

  do color = 1, ncolor
    i = subset(perm(color))
    adj(i,i) = color
  end do

  do i = 1, nnode
    if ( adj(i,i) == 0 ) then
      color = i4_uniform ( 1, ncolor, seed )
      adj(i,i) = color
    end if
  end do
!
!  Pick a random NEDGE subset.
!
  call ksub_random ( maxedge, nedge, seed, iwork )
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
subroutine color_graph_adj_color_count ( adj, lda, nnode, mcolor, ncolor )

!*****************************************************************************80
!
!! COLOR_GRAPH_ADJ_COLOR_COUNT counts the number of colors in a color graph.
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
!    Input, integer ( kind = 4 ) ADJ(LDA,NNODE), the adjacency information.  
!    ADJ(I,I) is the color of node I; otherwise, ADJ(I,J) is positive
!    if there is an edge between node I and node J.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the ADJ array,
!    which must be at least NNODE.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) MCOLOR, the maximum color index.
!
!    Output, integer ( kind = 4 ) NCOLOR, the number of colors.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(lda,nnode)
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

  call i4vec_sort_heap_a ( nnode, colors )

  call i4vec_uniq ( nnode, colors, ncolor )

  return
end
subroutine color_graph_adj_color_sequence ( adj, lda, nnode, seq )

!*****************************************************************************80
!
!! COLOR_GRAPH_ADJ_COLOR_SEQUENCE computes the color sequence of a color graph.
!
!  Discussion:
!
!    The color sequence of a color graph is constructed by computing the
!    color of each node, and then ordering these values in decreasing order.
!
!    If two color graphs are isomorphic, they must have the same color sequence.
!
!    If two color graphs have different color sequences, they cannot be
!    isomorphic.
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
!    Input, integer ( kind = 4 ) ADJ(LDA,NNODE), the adjacency information.  
!    ADJ(I,I) is the color of node I; otherwise, ADJ(I,J) is positive
!    if there is an edge between node I and node J.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the ADJ array,
!    which must be at least NNODE.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) SEQ(NNODE), the color sequence.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) seq(nnode)

  do i = 1, nnode
    seq(i) = adj(i,i)
  end do

  call i4vec_sort_heap_d ( nnode, seq )

  return
end
subroutine color_graph_adj_connect_random ( lda, nnode, nedge, &
  ncolor, seed, adj )

!*****************************************************************************80
!
!! COLOR_GRAPH_ADJ_CONNECT_RANDOM generates a random connected color graph.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of LDA, which must be
!    at least NNODE.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) NEDGE, the number of edges, which must be between
!    NNODE-1 and (NNODE*(NNODE-1))/2.  
!
!    Input, integer ( kind = 4 ) NCOLOR, the number of colors available to choose for
!    the nodes.  NCOLOR must be at least 1, and no more than NNODE.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
!    Output, integer ( kind = 4 ) ADJ(LDA,NNODE), the adjacency information.  
!    ADJ(I,I) is the color of node I; otherwise, ADJ(I,J) is positive
!    if there is an edge between node I and node J.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) ncolor
  integer ( kind = 4 ) nnode
  integer ( kind = 4 ) nedge

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) code(nnode-2)
  integer ( kind = 4 ) color
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) inode(nnode-1)
  integer ( kind = 4 ) iwork(nedge)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jnode(nnode-1)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) maxedge
  integer ( kind = 4 ) nchoice
  integer ( kind = 4 ) nchoose
  integer ( kind = 4 ) nnode2
  integer ( kind = 4 ) perm(ncolor)
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) subset(ncolor)
!
!  Check.
!
  if ( nnode <= 0  ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'COLOR_GRAPH_ADJ_CONNECT_RANDOM - Fatal error!'
    write ( *, '(a,i8)' ) '  NNODE = ', nnode
    write ( *, '(a)' ) '  but NNODE must be at least 1.'
    stop
  end if

  if ( lda < nnode ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'COLOR_GRAPH_ADJ_CONNECT_RANDOM - Fatal error!'
    write ( *, '(a,i8)' ) '  LDA = ', lda
    write ( *, '(a,i8)' ) '  but LDA must be at least NNODE = ', nnode
    stop
  end if

  maxedge = ( nnode * ( nnode - 1 ) ) / 2

  if ( nedge < nnode-1 .or. maxedge < nedge ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'COLOR_GRAPH_ADJ_CONNECT_RANDOM - Fatal error!'
    write ( *, '(a,i8)' ) '  NEDGE = ', nedge
    write ( *, '(a)' ) '  but NEDGE must be at least 0, and '
    write ( *, '(a,i8)' ) '  no more than ', maxedge
    stop
  end if

  if ( ncolor < 1 .or. nnode < ncolor ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'COLOR_GRAPH_ADJ_CONNECT_RANDOM - Fatal error!'
    write ( *, '(a,i8)' ) '  NCOLOR = ', ncolor
    write ( *, '(a)' ) '  but NCOLOR must be at least 1, and '
    write ( *, '(a,i8)' ) '  no more than ', nnode
    stop
  end if
!
!  Initialize the adjacency matrix.
!
  adj(1:nnode,1:nnode) = 0
!
!  Choose the colors.
!
  call ksub_random ( nnode, ncolor, seed, subset )

  call perm_random ( ncolor, seed, perm )

  do color = 1, ncolor
    i = subset(perm(color))
    adj(i,i) = color
  end do

  do i = 1, nnode
    if ( adj(i,i) == 0 ) then
      color = i4_uniform ( 1, ncolor, seed )
      adj(i,i) = color
    end if
  end do
!
!  Pick a random tree.
!
  call tree_arc_random ( nnode, seed, code, inode, jnode )
!
!  Convert information to adjacency form.
!
  call graph_arc_to_graph_adj ( nnode-1, inode, jnode, adj, lda, nnode2 )
!
!  Now we have NEDGE - ( NNODE - 1 ) more edges to add.
!
  nchoice = ( nnode * ( nnode - 1 ) ) / 2 - ( nnode - 1 )
  nchoose = nedge - ( nnode - 1 )

  call ksub_random ( nchoice, nchoose, seed, iwork )

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
subroutine color_graph_adj_degree ( adj, lda, nnode, degree )

!*****************************************************************************80
!
!! COLOR_GRAPH_ADJ_DEGREE computes the degree of each node.
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
!    Output, integer ( kind = 4 ) ADJ(LDA,NNODE), the adjacency information.  
!    ADJ(I,I) is the color of node I; otherwise, ADJ(I,J) is positive
!    if there is an edge between node I and node J.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the ADJ array,
!    which must be at least NNODE.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) DEGREE(NNODE), the degree of the nodes.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(lda,nnode)
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
subroutine color_graph_adj_degree_seq ( adj, lda, nnode, seq )

!*****************************************************************************80
!
!! COLOR_GRAPH_ADJ_DEGREE_SEQ computes the degree sequence of a color graph.
!
!  Discussion:
!
!    The degree sequence of a graph is constructed by computing the
!    degree of each node, and then ordering these values in decreasing order.
!
!    If two graphs are isomorphic, they must have the same degree sequence.
!
!    If two graphs have different degree sequences, they cannot be isomorphic.
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
!    Input, integer ( kind = 4 ) ADJ(LDA,NNODE), the adjacency information.  
!    ADJ(I,I) is the color of node I; otherwise, ADJ(I,J) is positive
!    if there is an edge between node I and node J.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the ADJ array,
!    which must be at least NNODE.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) SEQ(NNODE), the degree sequence.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) seq(nnode)

  call color_graph_adj_degree ( adj, lda, nnode, seq )

  call i4vec_sort_heap_d ( nnode, seq )

  return
end
subroutine color_graph_adj_edge_count ( adj, lda, nnode, nedge )

!*****************************************************************************80
!
!! COLOR_GRAPH_ADJ_EDGE_COUNT counts the number of edges in a color graph.
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
!    Input, integer ( kind = 4 ) ADJ(LDA,NNODE), the adjacency information.  
!    ADJ(I,I) is the color of node I; otherwise, ADJ(I,J) is positive
!    if there is an edge between node I and node J.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the ADJ array,
!    which must be at least NNODE.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) NEDGE, the number of edges.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(lda,nnode)
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
subroutine color_graph_adj_example_bush ( adj, lda, nnode )

!*****************************************************************************80
!
!! COLOR_GRAPH_ADJ_EXAMPLE_BUSH sets up the bush color graph.
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
!    Output, integer ( kind = 4 ) ADJ(LDA,NNODE), the adjacency information.  
!    ADJ(I,I) is the color of node I; otherwise, ADJ(I,J) is positive
!    if there is an edge between node I and node J.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the ADJ array,
!    which must be at least NNODE.
!
!    Output, integer ( kind = 4 ) NNODE, the number of nodes.
!
  implicit none

  integer ( kind = 4 ), parameter :: BLUE = 1
  integer ( kind = 4 ), parameter :: GREEN = 2
  integer ( kind = 4 ), parameter :: RED = 3
  integer ( kind = 4 ), parameter :: WHITE = 4

  integer ( kind = 4 ) lda

  integer ( kind = 4 ) adj(lda,lda)
  integer ( kind = 4 ) nnode

  nnode = 7

  if ( lda < nnode ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'COLOR_GRAPH_ADJ_EXAMPLE_BUSH - Fatal error!'
    write ( *, '(a)' ) '  LDA is too small!'
    stop
  end if

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
subroutine color_graph_adj_example_cube ( adj, lda, nnode )

!*****************************************************************************80
!
!! COLOR_GRAPH_ADJ_EXAMPLE_CUBE sets up the cube color graph.
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
!    Output, integer ( kind = 4 ) ADJ(LDA,NNODE), the adjacency information.  
!    ADJ(I,I) is the color of node I; otherwise, ADJ(I,J) is positive
!    if there is an edge between node I and node J.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the ADJ array,
!    which must be at least NNODE.
!
!    Output, integer ( kind = 4 ) NNODE, the number of nodes.
!
  implicit none

  integer ( kind = 4 ), parameter :: BLUE = 1
  integer ( kind = 4 ), parameter :: GREEN = 2
  integer ( kind = 4 ), parameter :: RED = 3

  integer ( kind = 4 ) lda

  integer ( kind = 4 ) adj(lda,lda)
  integer ( kind = 4 ) nnode

  nnode = 8

  if ( lda < nnode ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'COLOR_GRAPH_ADJ_EXAMPLE_CUBE - Fatal error!'
    write ( *, '(a)' ) '  LDA is too small.'
    stop
  end if

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
subroutine color_graph_adj_example_octo ( lda, example, seed, nnode, adj )

!*****************************************************************************80
!
!! COLOR_GRAPH_ADJ_EXAMPLE_OCTO sets up an 8 node example color graph.
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
!    There are 7 graphs to choose from.  They are all on 8 nodes.  The first
!    5 have degree 3 at every node.  Graphs 6 and 7 have degree 5 at every
!    node.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the ADJ array,
!    which must be at least NNODE.
!
!    Input, integer ( kind = 4 ) EXAMPLE, should be between 1 and 35, and indicates
!    which example graph to pick.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
!    Output, integer ( kind = 4 ) NNODE, the number of nodes, which should be 8.
!
!    Output, integer ( kind = 4 ) ADJ(LDA,LDA), the adjacency information.
!    ADJ(I,I) is the color of node I.
!    ADJ(I,J) is 1 if nodes I and J are adjacent and 0 otherwise.
!
  implicit none

  integer ( kind = 4 ), parameter :: BLUE = 1
  integer ( kind = 4 ), parameter :: GREEN = 2
  integer ( kind = 4 ), parameter :: RED = 3
  integer ( kind = 4 ), parameter :: YELLOW = 4

  integer ( kind = 4 ) lda

  integer ( kind = 4 ) adj(lda,lda)
  integer ( kind = 4 ) example
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) j
  integer ( kind = 4 ) msave
  integer ( kind = 4 ) nnode
  integer ( kind = 4 ) nsave
  integer ( kind = 4 ) seed

  if ( example <= 0 ) then
    nsave = i4_uniform ( 1, 7, seed )
    msave = i4_uniform ( 1, 5, seed )
  else
    example = mod ( example - 1, 35 ) + 1
    msave = ( ( example - 1 ) / 7 ) + 1
    nsave = mod ( example - 1, 7 ) + 1
  end if

  nnode = 8

  if ( lda < nnode ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'COLOR_GRAPH_ADJ_EXAMPLE_OCTO - Fatal error!'
    write ( *, '(a)' ) '  LDA is too small.'
    stop
  end if

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
    adj(4,1) = 1

    adj(5,8) = 1
    adj(5,1) = 1
    adj(5,2) = 1

    adj(6,1) = 1
    adj(6,2) = 1
    adj(6,3) = 1

    adj(7,2) = 1
    adj(7,3) = 1
    adj(7,4) = 1

    adj(8,3) = 1
    adj(8,4) = 1
    adj(8,5) = 1
!
!  Underlying graph 7.
!
  else if ( nsave == 7 ) then

    adj(1,3) = 1
    adj(1,5) = 1
    adj(1,7) = 1

    adj(2,4) = 1
    adj(2,6) = 1
    adj(2,8) = 1

    adj(3,5) = 1
    adj(3,7) = 1
    adj(3,1) = 1

    adj(4,6) = 1
    adj(4,8) = 1
    adj(4,2) = 1

    adj(5,7) = 1
    adj(5,1) = 1
    adj(5,3) = 1

    adj(6,8) = 1
    adj(6,2) = 1
    adj(6,4) = 1

    adj(7,1) = 1
    adj(7,3) = 1
    adj(7,5) = 1

    adj(8,2) = 1
    adj(8,4) = 1
    adj(8,6) = 1

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
  call i4mat_perm_random ( lda, nnode, seed, adj )

  return
end
subroutine color_graph_adj_example_twig ( adj, lda, nnode )

!*****************************************************************************80
!
!! COLOR_GRAPH_ADJ_EXAMPLE_TWIG sets up the twig color graph.
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
!    Output, integer ( kind = 4 ) ADJ(LDA,NNODE), the adjacency information.  
!    ADJ(I,I) is the color of node I; otherwise, ADJ(I,J) is positive
!    if there is an edge between node I and node J.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the ADJ array,
!    which must be at least NNODE.
!
!    Output, integer ( kind = 4 ) NNODE, the number of nodes.
!
  implicit none

  integer ( kind = 4 ), parameter :: BLUE = 1
  integer ( kind = 4 ), parameter :: RED = 3

  integer ( kind = 4 ) lda

  integer ( kind = 4 ) adj(lda,lda)
  integer ( kind = 4 ) nnode

  nnode = 3

  if ( lda < nnode ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'COLOR_GRAPH_ADJ_EXAMPLE_TWIG - Fatal error!'
    write ( *, '(a)' ) '  LDA is too small!'
    stop
  end if

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
subroutine color_graph_adj_print ( adj, lda, nnode, title )

!*****************************************************************************80
!
!! COLOR_GRAPH_ADJ_PRINT prints out the adjacency matrix of a color graph.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
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
!    Input, integer ( kind = 4 ) ADJ(LDA,NNODE), the adjacency information.  
!    ADJ(I,I) is the color of node I; otherwise, ADJ(I,J) is positive
!    if there is an edge between node I and node J.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of ADJ, which must be
!    least NNODE.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  character ( len = 80 ) string
  character ( len = * ) title

  if ( len_trim ( title ) /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  do i = 1, nnode

    do j = 1, nnode

      k = (j-1) * 3 + 1
      write ( string(k:k+2), '(i3)' ) adj(i,j)

    end do

    write ( *, '(i2,2x,a)' ) i, string(1:3*nnode)

  end do

  return
end
subroutine color_graph_adj_random ( lda, nnode, ncolor, nedge, seed, adj )

!*****************************************************************************80
!
!! COLOR_GRAPH_ADJ_RANDOM generates a random color graph.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of LDA, which must be
!    at least NNODE.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) NCOLOR, the number of colors available to choose for
!    the nodes.  NCOLOR must be at least 1, and no more than NNODE.
!
!    Input, integer ( kind = 4 ) NEDGE, the number of edges, which must be between
!    0 and (NNODE*(NNODE-1))/2.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
!    Output, integer ( kind = 4 ) ADJ(LDA,NNODE), the adjacency information.  
!    ADJ(I,I) is the color of node I; otherwise, ADJ(I,J) is positive
!    if there is an edge between node I and node J.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) nnode
  integer ( kind = 4 ) nedge

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) color
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) iwork(nedge)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) maxedge
  integer ( kind = 4 ) ncolor
  integer ( kind = 4 ) perm(ncolor)
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) subset(ncolor)

  if ( nnode <= 0  ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'COLOR_GRAPH_ADJ_RANDOM - Fatal error!'
    write ( *, '(a,i8)' ) '  NNODE = ', nnode
    write ( *, '(a)' ) '  but NNODE must be at least 1.'
    stop
  end if

  maxedge = ( nnode * ( nnode - 1 ) ) / 2

  if ( nedge < 0 .or. maxedge < nedge ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'COLOR_GRAPH_ADJ_RANDOM - Fatal error!'
    write ( *, '(a,i8)' ) '  NEDGE = ', nedge
    write ( *, '(a)' ) '  but NEDGE must be at least 0, and '
    write ( *, '(a,i8)' ) '  no more than ', maxedge
    stop
  end if

  if ( ncolor < 1 .or. nnode < ncolor ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'COLOR_GRAPH_ADJ_RANDOM - Fatal error!'
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
  call ksub_random ( nnode, ncolor, seed, subset )

  call perm_random ( ncolor, seed, perm )

  do color = 1, ncolor
    i = subset(perm(color))
    adj(i,i) = color
  end do

  do i = 1, nnode
    if ( adj(i,i) == 0 ) then
      color = i4_uniform ( 1, ncolor, seed )
      adj(i,i) = color
    end if
  end do
!
!  Pick a random NEDGE subset of (N*(N-1))/2.
!
  call ksub_random ( maxedge, nedge, seed, iwork )
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
subroutine degree_seq_is_graphic ( nnode, seq, result )

!*****************************************************************************80
!
!! DEGREE_SEQ_IS_GRAPHIC reports whether a degree sequence represents a graph.
!
!  Discussion:
!
!    The degree sequence of a graph is constructed by computing the
!    degree of each node, and then ordering these values in decreasing order.
!
!    A sequence of NNODE nonnegative integers is said to be "graphic" if
!    there exists a graph whose degree sequence is the given sequence.
!
!    The Havel Hakimi theorem states that 
!
!      s t1 t2 ... ts d1 d2 ... dn
!
!    is graphic if and only if
!
!        t1-1 t2-1 ... ts-1 d1 d2 ... dn
!
!    is graphic (after any necessary resorting and dropping of zeroes).
!    Definitely, the one thing we cannot have is that any nonzero entry
!    is equal to or greater than the number of nonzero entries.
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
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) SEQ(NNODE), the degree sequence to be tested.
!
!    Output, integer ( kind = 4 ) RESULT, the result.
!    0, SEQ is not graphic.
!    1, SEQ is graphic.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) dmax
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4vec_nonzero
  integer ( kind = 4 ) nonzero
  integer ( kind = 4 ) order
  integer ( kind = 4 ) result
  integer ( kind = 4 ) seq(nnode)

  result = 0

  do i = 1, nnode
    if ( seq(i) < 0 ) then
      return
    end if
  end do
!
!  Check that the sequence is decreasing.
!
  call i4vec_order_type ( nnode, seq, order )

  if ( order == -1 .or. order == 1 .or. order == 2 ) then
    return
  end if
!
!  Now apply the Havel Hakimi theorem.
!
  do

    nonzero = i4vec_nonzero ( nnode, seq )

    if ( nonzero == 0 ) then
      result = 1
      exit
    end if

    call i4vec_sort_heap_d ( nnode, seq )

    dmax = seq(1)

    if ( nonzero <= dmax ) then
      result = 0
      exit
    end if

    seq(1) = 0
    do i = 2, dmax + 1
      seq(i) = seq(i) - 1
    end do

  end do
        
  return
end
subroutine degree_seq_to_graph_adj ( nnode, seq, lda, adj, ierror )

!*****************************************************************************80
!
!! DEGREE_SEQ_TO_GRAPH_ADJ computes a graph with the given degree sequence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
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
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) SEQ(NNODE), the degree sequence.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of ADJ.
!
!    Output, integer ( kind = 4 ) ADJ(LDA,NNODE), the adjacency information.  
!    ADJ(I,J) is nonzero if there is an edge from node I to node J.
!
!    Output, integer ( kind = 4 ) IERROR, is nonzero if an error occurred.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) indx(nnode)
  integer ( kind = 4 ) nonzero
  integer ( kind = 4 ) s
  integer ( kind = 4 ) seq(nnode)
  integer ( kind = 4 ) seq2(nnode)

  ierror = 0

  adj(1:nnode,1:nnode) = 0

  seq2(1:nnode) = seq(1:nnode)

  do

    call i4vec_sort_heap_index_d ( nnode, seq2, indx )

    nonzero = 0
    do i = 1, nnode
      if ( seq2(i) /= 0 ) then
        nonzero = nonzero + 1
      end if
    end do

    if ( nonzero == 0 ) then
      exit
    end if

    s = seq2(indx(1))

    if ( nonzero <= s ) then
      ierror = 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DEGREE_SEQ_TO_GRAPH_ADJ - Fatal error!'
      write ( *, '(a)' ) '  The degree sequence is not graphic!'
      return
    end if

    seq2(indx(1)) = 0

    do i = 2, s+1
      adj(indx(i),indx(1)) = 1
      adj(indx(1),indx(i)) = 1
      seq2(indx(i)) = seq2(indx(i)) - 1
    end do

  end do

  return
end
subroutine dge_check ( lda, m, n, ierror )

!*****************************************************************************80
!
!! DGE_CHECK checks the dimensions of a general matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the array.
!    LDA must be at least M.
!
!    Input, integer ( kind = 4 ) M, the number of rows of the matrix.
!    M must be positive.
!
!    Input, integer ( kind = 4 ) N, the number of columns of the matrix.
!    N must be positive.
!
!    Output, integer ( kind = 4 ) IERROR, reports whether any errors were detected.
!    IERROR is set to 0 before the checks are made, and then:
!    IERROR = IERROR + 1 if LDA is illegal;
!    IERROR = IERROR + 2 if M is illegal;
!    IERROR = IERROR + 4 if N is illegal.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) lda
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  ierror = 0

  if ( lda < m ) then
    ierror = ierror + 1
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) 'DGE_CHECK - Illegal LDA = ', lda
  end if

  if ( m < 1 ) then
    ierror = ierror + 2
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) 'DGE_CHECK - Illegal M = ', m
  end if

  if ( n < 1 ) then
    ierror = ierror + 4
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) 'DGE_CHECK - Illegal N = ', n
  end if

  return
end
subroutine dge_det ( lda, n, a, ipivot, det )

!*****************************************************************************80
!
!! DGE_DET computes the determinant of a matrix factored by DGE_FA or DGE_TRF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the array.
!    LDA must be at least N.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Input, real ( kind = 8 ) A(LDA,N), the LU factors computed 
!    by DGE_FA or DGE_TRF.
!
!    Input, integer ( kind = 4 ) IPIVOT(N), as computed by DGE_FA or DGE_TRF.
!
!    Output, real ( kind = 8 ) DET, the determinant of the matrix.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(lda,n)
  real ( kind = 8 ) det
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ipivot(n)
!
!  Check the dimensions.
!
  call dge_check ( lda, n, n, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DGE_DET - Fatal error!'
    write ( *, '(a)' ) '  Illegal dimensions.'
    return
  end if

  det = 1.0D+00

  do i = 1, n
    det = det * a(i,i)
  end do

  do i = 1, n
    if ( ipivot(i) /= i ) then
      det = - det
    end if
  end do

  return
end
subroutine dge_fa ( lda, n, a, ipivot, info )

!*****************************************************************************80
!
!! DGE_FA factors a general matrix.
!
!  Discussion:
!
!    DGE_FA is a simplified version of the LINPACK routine DGEFA.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the array.
!    LDA must be at least N.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Input/output, real ( kind = 8 ) A(LDA,N), the matrix to be factored.
!    On output, A contains an upper triangular matrix and the multipliers
!    which were used to obtain it.  The factorization can be written
!    A = L * U, where L is a product of permutation and unit lower
!    triangular matrices and U is upper triangular.
!
!    Output, integer ( kind = 4 ) IPIVOT(N), a vector of pivot indices.
!
!    Output, integer ( kind = 4 ) INFO, singularity flag.
!    0, no singularity detected.
!    nonzero, the factorization failed on the INFO-th step.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(lda,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipivot(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
!
!  Check the dimensions.
!
  call dge_check ( lda, n, n, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DGE_FA - Fatal error!'
    write ( *, '(a)' ) '  Illegal dimensions.'
    return
  end if

  info = 0

  do k = 1, n-1
!
!  Find L, the index of the pivot row.
!
    l = k
    do i = k+1, n
      if ( abs ( a(l,k) ) < abs ( a(i,k) ) ) then
        l = i
      end if
    end do

    ipivot(k) = l
!
!  If the pivot index is zero, the algorithm has failed.
!
    if ( a(l,k) == 0.0D+00 ) then
      info = k
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DGE_FA - Fatal error!'
      write ( *, '(a,i8)' ) '  Zero pivot on step ', info
      return
    end if
!
!  Interchange rows L and K if necessary.
!
    if ( l /= k ) then
      call r8_swap ( a(l,k), a(k,k) )
    end if
!
!  Normalize the values that lie below the pivot entry A(K,K).
!
    a(k+1:n,k) = -a(k+1:n,k) / a(k,k)
!
!  Row elimination with column indexing.
!
    do j = k+1, n

      if ( l /= k ) then
        call r8_swap ( a(l,j), a(k,j) )
      end if

      a(k+1:n,j) = a(k+1:n,j) + a(k+1:n,k) * a(k,j)

    end do

  end do

  ipivot(n) = n

  if ( a(n,n) == 0.0D+00 ) then
    info = n
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DGE_FA - Fatal error!'
    write ( *, '(a,i8)' ) '  Zero pivot on step ', info
  end if

  return
end
subroutine digraph_adj_closure ( adj, lda, nnode )

!*****************************************************************************80
!
!! DIGRAPH_ADJ_CLOSURE generates the transitive closure of a digraph.
!
!  Discussion:
!
!    The method is due to S Warshall.
!
!  Definition:
!
!    The transitive closure of a graph is a function REACH(I,J) so that
!
!      REACH(I,J) = 0 if node J cannot be reached from node I;
!                   1 if node J can be reached from node I.
!
!    This is an extension of the idea of adjacency.  REACH(I,J)=1 if
!    node J is adjacent to node I, or if node J is adjacent to a node
!    that is adjacent to node I, etc.
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
!  Reference:
!
!    Robert Sedgewick,
!    Algorithms,
!    Addison Wesley, 1983, page 425.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) ADJ(LDA,NNODE).
!
!    On input, ADJ is the adjacency matrix.  ADJ(I,J)
!    is nonzero if there is an edge from node I to node J.
!
!    On output, ADJ is the transitive closure matrix.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of LDA, which must be
!    at least NNODE.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
!
!  You can "reach" a node from itself.
!
  do i = 1, nnode
    adj(i,i) = 1
  end do

  do i = 1, nnode
    do j = 1, nnode
      if ( adj(j,i) /= 0 ) then
        do k = 1, nnode
          if ( adj(i,k) /= 0 ) then
            adj(j,k) = 1
          end if
        end do
      end if
    end do
  end do

  return
end
subroutine digraph_adj_components ( adj, lda, nnode, ncomp, comp, dad, order )

!*****************************************************************************80
!
!! DIGRAPH_ADJ_COMPONENTS finds the strongly connected components of a digraph.
!
!  Discussion:
!
!    A digraph is a directed graph.
!
!    A strongly connected component of a directed graph is the largest
!    set of nodes such that there is a directed path from any node to 
!    any other node in the same component.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 April 1999
!
!  Reference:
!
!    K Thulasiraman, M Swamy,
!    Graph Theory and Algorithms,
!    John Wiley, New York, 1992.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(LDA,NNODE), the adjacency information.
!    ADJ(I,J) is 1 if there is a direct link from node I to node J.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of ADJ.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) NCOMP, the number of strongly connected components.
!
!    Output, integer ( kind = 4 ) COMP(NNODE), lists the connected component to which 
!    each node belongs.
!
!    Output, integer ( kind = 4 ) DAD(NNODE), the father array for the depth first
!    search trees.  DAD(I) = 0 means that node I is the root of 
!    one of the trees.  DAD(I) = J means that the search descended
!    from node J to node I.
!
!    Output, integer ( kind = 4 ) ORDER(NNODE), the order in which the nodes were
!    traversed, from 1 to NNODE.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) comp(nnode)
  integer ( kind = 4 ) dad(nnode)
  integer ( kind = 4 ) iorder
  integer ( kind = 4 ) lowlink(nnode)
  integer ( kind = 4 ) mark(nnode)
  integer ( kind = 4 ) ncomp
  integer ( kind = 4 ) nstack
  integer ( kind = 4 ) order(nnode)
  integer ( kind = 4 ) point(nnode)
  integer ( kind = 4 ) stack(nnode)
  integer ( kind = 4 ) v
  integer ( kind = 4 ) w
  integer ( kind = 4 ) x
!
!  Initialization.
!
  comp(1:nnode) = 0
  dad(1:nnode) = 0
  order(1:nnode) = 0
  lowlink(1:nnode) = 0
  mark(1:nnode) = 0
  point(1:nnode) = 0

  iorder = 0
  nstack = 0
  ncomp = 0
!
!  Select any node V not stored in the stack, that is, with MARK(V) = 0.
!
  do

    v = 0

    do

      v = v + 1

      if ( nnode < v ) then
        adj(1:nnode,1:nnode) = abs ( adj(1:nnode,1:nnode) )
        return
      end if

      if ( mark(v) /= 1 ) then
        exit
      end if

    end do

    iorder = iorder + 1

    order(v) = iorder
    lowlink(v) = iorder
    mark(v) = 1
 
    nstack = nstack + 1
    stack(nstack) = v
    point(v) = 1

30  continue
!
!  Consider each node W.
!
    do w = 1, nnode
!
!  Is there an edge (V,W) and has it not been examined yet?
!
      if ( 0 < adj(v,w) ) then

        adj(v,w) = - adj(v,w)
!
!  Is the node on the other end of the edge undiscovered yet?
!
        if ( mark(w) == 0 ) then

          iorder = iorder + 1
          order(w) = iorder
          lowlink(w) = iorder
          dad(w) = v
          mark(w) = 1

          nstack = nstack + 1
          stack(nstack) = w
          point(w) = 1

          v = w

        else if ( mark(w) == 1 ) then

          if ( order(w) < order(v) .and. point(w) == 1 ) then
            lowlink(v) = min ( lowlink(v), order(w) )
          end if

        end if

        go to 30

      end if

    end do

    if ( lowlink(v) == order(v) ) then

      ncomp = ncomp + 1

      do

        if ( nstack <= 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'DIGRAPH_ADJ_COMPONENTS - Fatal error!'
          write ( *, '(a)' ) '  Illegal stack reference.'
          stop
        end if

        x = stack(nstack)
        nstack = nstack - 1

        point(x) = 0
        comp(x) = ncomp

        if ( x == v ) then
          exit
        end if

      end do

    end if

    if ( dad(v) /= 0 ) then
      lowlink(dad(v)) = min ( lowlink(dad(v)), lowlink(v) )
      v = dad(v)
      go to 30
    end if

  end do

  return
end
subroutine digraph_adj_cycle ( adj, lda, nnode, adj2, dad, order )

!*****************************************************************************80
!
!! DIGRAPH_ADJ_CYCLE searches for cycles in a digraph.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 July 2000
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(LDA,NNODE), the adjacency information.
!    ADJ(I,J) is 1 if there is a direct link from node I to node J.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of ADJ and ADJ2.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) ADJ2(LDA,NNODE), will be one of the following values 
!    depending on the role of the edge from node I to node J:
!       0, no edge,
!       1, neither in a search tree, nor needed to disconnect a cycle;
!      -1, completes a cycle,
!      -2, part of a search tree.
!
!    Output, integer ( kind = 4 ) DAD(NNODE), the father array for the depth first
!    search trees.  DAD(I) = 0 means that node I is the root of 
!    one of the trees.  DAD(I) = J means that the search descended
!    from node J to node I.
!
!    Output, integer ( kind = 4 ) ORDER(NNODE), the order in which the nodes were
!    traversed, from 1 to NNODE.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) adj2(lda,nnode)
  integer ( kind = 4 ) dad(nnode)
  integer ( kind = 4 ) daddy
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) maxstack
  integer ( kind = 4 ) nstack
  integer ( kind = 4 ) order(nnode)
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) stack(2*(nnode-1))
!
!  Initialization.
!
  adj2(1:nnode,1:nnode) = adj(1:nnode,1:nnode)
  dad(1:nnode) = 0
  maxstack = 2 * ( nnode - 1 )
  order(1:nnode) = 0

  rank = 0

  do i = 1, nnode

    if ( order(i) == 0 ) then

      daddy = i
      nstack = 0
!
!  Visit the unvisited node DAD.
!
10    continue

      rank = rank + 1
      order(daddy) = rank
      j = 0
!
!  Consider visiting node J from node DAD.
!
20    continue

      j = j + 1
!
!  If 
!    J is a reasonable value, 
!    J is adjacent to DAD, and 
!    J is unvisited,
!  then 
!    put DAD into the stack, 
!    make J the new value of DAD, and
!    examine J's neighbors.
!
      if ( j <= nnode ) then

        if ( 0 < adj2(daddy,j) ) then

          if ( order(j) == 0 ) then

            adj2(daddy,j) = -2

            if ( nstack+2 <= maxstack ) then
              dad(j) = daddy
              stack(nstack+1) = daddy
              stack(nstack+2) = j
              nstack = nstack + 2
              daddy = j
              go to 10
            else
              write ( *, '(a)' ) ' '
              write ( *, '(a)' ) 'DIGRAPH_ADJ_CYCLE - Fatal error!'
              write ( *, '(a)' ) '  Out of stack space.'
              stop
            end if
!
!  Adjacent node J has already been visited.  If J is actually
!  in the current stack, then we have a cycle.
!
          else

            if ( j == daddy ) then

              adj2(daddy,j) = - 1

            else

              do jj = 1, nstack-1, 2
                if ( stack(jj) == j ) then
                  adj2(daddy,j) = - 1
                end if
              end do

            end if

            go to 20

          end if
!
!  If J is not suitable for a visit, get the next value of J.
!
        else

          go to 20

        end if
!
!  If no more neighbors to consider, back up one node.
!
      else if ( 2 <= nstack ) then

        daddy = stack(nstack-1)
        j = stack(nstack)
        nstack = nstack - 2
        go to 20
!
!  If no more nodes to consider in this tree, bail out.
!
      else

        nstack = 0

      end if

    end if

  end do

  return
end
subroutine digraph_adj_degree ( adj, lda, nnode, indegree, outdegree )

!*****************************************************************************80
!
!! DIGRAPH_ADJ_DEGREE computes the indegree and outdegree of each node.
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
!    Input, integer ( kind = 4 ) ADJ(LDA,NNODE), the adjacency information.
!    ADJ(I,J) is 1 if there is a direct link from node I to node J.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the ADJ array,
!    which must be at least NNODE.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) INDEGREE(NNODE), OUTDEGREE(NNODE), 
!    the indegree and outdegree of the nodes.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(lda,nnode)
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
subroutine digraph_adj_degree_max ( adj, lda, nnode, indegree_max, &
  outdegree_max, degree_max )

!*****************************************************************************80
!
!! DIGRAPH_ADJ_DEGREE_MAX computes the maximum degrees of a digraph.
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
!    Input, integer ( kind = 4 ) ADJ(LDA,NNODE), the adjacency information.
!    ADJ(I,J) is 1 if there is a direct link from node I to node J.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the ADJ array,
!    which must be at least NNODE.
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

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) degree
  integer ( kind = 4 ) degree_max
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indegree
  integer ( kind = 4 ) indegree_max
  integer ( kind = 4 ) outdegree
  integer ( kind = 4 ) outdegree_max

  degree_max = 0
  indegree_max = 0
  outdegree_max = 0

  do i = 1, nnode

    indegree = sum ( adj(1:nnode,i) )
    outdegree = sum ( adj(i,1:nnode) )

    degree = indegree + outdegree

    indegree_max = max ( indegree_max, indegree )
    outdegree_max = max ( outdegree_max, outdegree )
    degree_max = max ( degree_max, degree )
     
  end do

  return
end
subroutine digraph_adj_degree_seq ( adj, lda, nnode, in_seq, out_seq )

!*****************************************************************************80
!
!! DIGRAPH_ADJ_DEGREE_SEQ computes the directed degree sequence.
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
!    22 October 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(LDA,NNODE), the adjacency information.
!    ADJ(I,J) is 1 if there is a direct link from node I to node J.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the ADJ array,
!    which must be at least NNODE.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) IN_SEQ(NNODE), OUT_SEQ(NNODE),
!    the degree sequence of the digraph.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) in_seq(nnode)
  integer ( kind = 4 ) out_seq(nnode)

  call digraph_adj_degree ( adj, lda, nnode, in_seq, out_seq )

  call i4vec2_sort_d ( nnode, out_seq, in_seq )

  return
end
subroutine digraph_adj_edge_count ( adj, lda, nnode, nedge )

!*****************************************************************************80
!
!! DIGRAPH_ADJ_EDGE_COUNT counts the number of edges in a digraph.
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
!    Input, integer ( kind = 4 ) ADJ(LDA,NNODE), the adjacency information.
!    ADJ(I,J) is 1 if there is an edge from node I to node J.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the ADJ array,
!    which must be at least NNODE.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) NEDGE, the number of edges in the digraph.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) nedge

  nedge = sum ( adj(1:nnode,1:nnode) )

  return
end
subroutine digraph_adj_eigen ( adj, lda, nnode, neigen, eigenr, eigeni )

!*****************************************************************************80
!
!! DIGRAPH_ADJ_EIGEN computes the eigenvalues of a digraph from its adjacency matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(LDA,NNODE), the adjacency information.
!    ADJ(I,J) is 1 if there is an edge from node I to node J.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the ADJ array,
!    which must be at least NNODE.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) NEIGEN, the number of eigenvalues computed.
!    Normally, this would be equal to NNODE, unless the algorithm failed.
!
!    Output, real ( kind = 8 ) EIGENR(NNODE), EIGENI(NNODE), contains the real
!    and imaginary parts of the computed eigenvalues.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) nnode

  real ( kind = 8 ) a(nnode,nnode)
  integer ( kind = 4 ) adj(lda,nnode)
  real ( kind = 8 ) eigeni(nnode)
  real ( kind = 8 ) eigenr(nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) igh
  integer ( kind = 4 ) ind(nnode)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) low
  integer ( kind = 4 ) neigen
  real ( kind = 8 ) scale(nnode)

  a(1:nnode,1:nnode) = real ( adj(1:nnode,1:nnode), kind = 8 )

  call balanc ( nnode, nnode, a, low, igh, scale )

  call elmhes ( nnode, nnode, low, igh, a, ind )

  call hqr ( nnode, nnode, low, igh, a, eigenr, eigeni, info )

  if ( info == 0 ) then
    neigen = nnode
  else
    neigen = nnode - info
    do i = 1, neigen
      eigenr(i) = eigenr(i+info)
      eigeni(i) = eigeni(i+info)
    end do
  end if

  return
end
subroutine digraph_adj_example_cycler ( adj, lda, nnode )

!*****************************************************************************80
!
!! DIGRAPH_ADJ_EXAMPLE_CYCLER sets up the adjacency information for the cycler digraph.
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
!    Output, integer ( kind = 4 ) ADJ(LDA,NNODE), the adjacency information.
!    ADJ(I,J) is 1 if there is a direct link from node I to node J.
!
!    Input, integer ( kind = 4 ) LDA, the maximum value of NNODE, which must be at least 9.
!
!    Output, integer ( kind = 4 ) NNODE, the number of nodes.
!
  implicit none

  integer ( kind = 4 ) lda

  integer ( kind = 4 ) adj(lda,lda)
  integer ( kind = 4 ) nnode

  nnode = 9

  if ( lda < nnode ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DIGRAPH_ADJ_EXAMPLE_CYCLER - Fatal error!'
    write ( *, '(a)' ) '  LDA is too small.'
    stop
  end if

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
subroutine digraph_adj_example_octo ( lda, example, seed, nnode, adj )

!*****************************************************************************80
!
!! DIGRAPH_ADJ_EXAMPLE_OCTO sets up an 8 node example digraph.
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
!    There are 12 digraphs to choose from, all on 8 nodes.  There are 7
!    underlying graphs.  The first 5 underlying graphs have degree 3 at 
!    every node.  Graphs 6 and 7 have degree 5 at every node.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the ADJ array,
!    which must be at least NNODE.
!
!    Input, integer ( kind = 4 ) EXAMPLE, should be between 1 and 12, and indicates
!    which example graph to pick.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
!    Output, integer ( kind = 4 ) NNODE, the number of nodes, which should be 8.
!
!    Output, integer ( kind = 4 ) ADJ(LDA,LDA), the adjacency information.
!    ADJ(I,J) is 1 if nodes I and J are adjacent and 0 otherwise.
!
  implicit none

  integer ( kind = 4 ) lda

  integer ( kind = 4 ) adj(lda,lda)
  integer ( kind = 4 ) example
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nnode
  integer ( kind = 4 ) nsave
  integer ( kind = 4 ) seed

  if ( example <= 0 ) then
    nsave = i4_uniform ( 1, 12, seed )
  else
    example = mod ( example - 1, 12 ) + 1
    nsave = example
  end if

  nnode = 8

  if ( lda < nnode ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DIGRAPH_ADJ_EXAMPLE_OCTO - Fatal error!'
    write ( *, '(a)' ) '  LDA is too small.'
    stop
  end if

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

  end if
!
!  Now permute the graph.
!
  call i4mat_perm_random ( lda, nnode, seed, adj )

  return
end
subroutine digraph_adj_example_sixty ( adj, lda, nnode )

!*****************************************************************************80
!
!! DIGRAPH_ADJ_EXAMPLE_SIXTY sets the adjacency matrix for the sixty digraph.
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
!    11 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) ADJ(LDA,NNODE), the adjacency information.
!    ADJ(I,J) is 1 if there is a direct link from node I to node J.
!
!    Input, integer ( kind = 4 ) LDA, the maximum value of NNODE, which must be at least 12.
!
!    Output, integer ( kind = 4 ) NNODE, the number of nodes.
!
  implicit none

  integer ( kind = 4 ) lda

  integer ( kind = 4 ) adj(lda,lda)
  integer ( kind = 4 ) d(12)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nnode

  nnode = 12

  if ( lda < nnode ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DIGRAPH_ADJ_EXAMPLE_SIXTY - Fatal error!'
    write ( *, '(a)' ) '  LDA is too small.'
    stop
  end if

  d(1:12) = (/ 60, 30, 20, 15, 12, 10, 6, 5, 4, 3, 2, 1 /)

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
subroutine digraph_adj_ham_cand ( adj, lda, nnode, circuit, k, nstack, &
  stack, maxstack, ncan )

!*****************************************************************************80
!
!! DIGRAPH_ADJ_HAM_CAND finds candidates for the next node in a Hamiltonian circuit.
!
!  Discussion:
!
!    This routine is used in conjunction with I4VEC_BACKTRACK.  
!
!    A Hamiltonian circuit of a digraph is a path that starts at a given node, 
!    visits every node exactly once, and returns to the starting node.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 August 2000
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
!    Input, integer ( kind = 4 ) ADJ(LDA,NNODE).  ADJ(I,J) = 1 if there is
!    an edge from node I to node J, 0 otherwise.
!
!    Input, integer ( kind = 4 ) LDA, the first dimension of ADJ.
!    LDA must be at least NNODE.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes in the digraph.
!
!    Input, integer ( kind = 4 ) CIRCUIT(NNODE), contains in CIRCUIT(1:K-1) the partial 
!    candidate circuit being constructed.
!
!    Input, integer ( kind = 4 ) K, the index of the next node to be determined for 
!    the circuit.
!
!    Input/output, integer ( kind = 4 ) NSTACK, the current length of the stack.
!
!    Input, integer ( kind = 4 ) STACK(MAXSTACK), candidates for positions 1...K-1.
!
!    Input, integer ( kind = 4 ) MAXSTACK, the dimension of STACK.
!
!    Input/output, integer ( kind = 4 ) NCAN(NNODE), the number of candidates for each 
!    position.  On input, contains values for steps 1 to K-1.  On output, 
!    the value for position K has been determined.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) nnode
  integer ( kind = 4 ) maxstack

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) circuit(nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iwork(nnode)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) ncan(nnode)
  integer ( kind = 4 ) nstack
  integer ( kind = 4 ) stack(maxstack)

  ncan(k) = 0

  if ( k == 1 ) then
    stack(1) = 1
    nstack = 1
    ncan(k) = 1
    return
  end if

  iwork(1:nnode) = adj(circuit(k-1),1:nnode)
 
  iwork(circuit(1:k-1)) = 0
  
  if ( k /= nnode ) then
 
    do i = 1, nnode
      if ( iwork(i) == 1 ) then
        if ( maxstack <= nstack ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'DIGRAPH_ADJ_HAM_CAND - Fatal error!'
          write ( *, '(a)' ) '  MAXSTACK is too small.'
          stop
        end if
        nstack = nstack + 1
        stack(nstack) = i
        ncan(k) = ncan(k) + 1
      end if
    end do
 
    return
 
  else if ( k == nnode ) then
 
    do i = 1, nnode
 
      if ( iwork(i) == 1 ) then
 
        if ( adj(i,1) /= 0 ) then
          if ( maxstack <= nstack ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'DIGRAPH_ADJ_HAM_CAND - Fatal error!'
            write ( *, '(a)' ) '  MAXSTACK is too small.'
            stop
          end if
          nstack = nstack + 1
          stack(nstack) = i
          ncan(k) = ncan(k) + 1
        end if

        return
 
      end if
 
    end do

  end if
 
  return
end
subroutine digraph_adj_ham_next ( adj, lda, nnode, circuit, stack, &
  maxstack, ncan, more )

!*****************************************************************************80
!
!! DIGRAPH_ADJ_HAM_NEXT returns the next Hamilton circuit for a digraph.
!
!  Discussion:
!
!    The routine produces all the Hamilton circuits of a digraph, one at a time.
!
!    A Hamiltonian circuit of a digraph is a path that starts at a given
!    node, visits every node exactly once, and returns to the starting node.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 August 2000
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
!    Input, integer ( kind = 4 ) ADJ(LDA,NNODE), the adjacency matrix of the digraph.  
!    ADJ(I,J) = 1 if there is an edge from node I to node J, 0 otherwise.
!
!    Input, integer ( kind = 4 ) LDA, the first dimension of ADJ as
!    declared in the calling program.  LDA must be at least NNODE.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes in the digraph.
!
!    Input/output, integer ( kind = 4 ) CIRCUIT(NNODE).  On the first call to this routine,
!    the contents of CIRCUIT are irrelevant.  On return, CIRCUIT contains a
!    list of the nodes that form a cirucit.  On each subsequent call, 
!    the input value of CIRCUIT is used to construct the next solution,
!    so the user should not alter the contents of CIRCUIT during a computation.
!
!    Workspace, integer STACK(MAXSTACK).  Candidates for the positions in
!    the circuit.
!
!    Input, integer ( kind = 4 ) MAXSTACK, the dimension of STACK.
!
!    Workspace, integer NCAN(NNODE), a count of the number of candidates for 
!    each step.
!
!    Input/output, logical MORE.
!    On first call, set MORE to .FALSE, and do not alter it after.
!    On return, MORE is TRUE if another circuit has been returned in
!    IARRAY, and FALSE if there are no more circuits.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) nnode
  integer ( kind = 4 ) maxstack

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) circuit(nnode)
  integer ( kind = 4 ), save :: indx = 0
  integer ( kind = 4 ), save :: k = 0
  logical more
  integer ( kind = 4 ) ncan(nnode)
  integer ( kind = 4 ), save :: nstack = 0
  integer ( kind = 4 ) stack(maxstack)

  if ( .not. more ) then
    indx = 0
    k = 0
    more = .true.
    nstack = 0
  end if
 
  do
 
    call i4vec_backtrack ( nnode, circuit, indx, k, nstack, stack, maxstack, &
      ncan )
 
    if ( indx == 1 ) then

      exit

    else if ( indx == 2 ) then

      call digraph_adj_ham_cand ( adj, lda, nnode, circuit, k, nstack, &
        stack, maxstack, ncan )

    else

      more = .false.
      exit

    end if

  end do
 
  return
end
subroutine digraph_adj_ham_next_brute ( adj, lda, nnode, circuit, iset )

!*****************************************************************************80
!
!! DIGRAPH_ADJ_HAM_NEXT_BRUTE finds the next Hamiltonian circuit in a digraph.
!
!  Discussion:
!
!    This is a brute force algorithm, and not suitable for large problems.
!    It is really only useful as a demonstration, and as a check for
!    the backtracking algorithm.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(LDA,NNODE), the adjacency information.
!    ADJ(I,J) is 1 if there is a direct link from node I to node J.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of ADJ, which must be
!    at least NNODE.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input/output, integer ( kind = 4 ) CIRCUIT(NNODE).
!
!    On input, if ISET = 0, then CIRCUIT is not presumed to contain any 
!    information.  If ISET is nonzero, then CIRCUIT contains the circuit 
!    computed on the previous call.
!
!    On output, CIRCUIT contains the circuit computed by this call.
!
!    Input/output, integer ( kind = 4 ) ISET.
!    On input, 0 means this is the first call for this graph.  
!    Any other value means this is a repeated call for more circuits.
!
!    On output, a 0 value means that no more circuits could be computed.
!    Otherwise, ISET is incremented by one on each call.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) circuit(nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ipos
  integer ( kind = 4 ) iset
!
!  If ISET is 0, this is a starting call, and we set CIRCUIT
!  to the lexically first circuit to check.
!
!  Otherwise, we set CIRCUIT to the next permutation.
!
  if ( iset == 0 ) then
    ipos = 0
    circuit(1:nnode) = 0
  else
    ipos = nnode - 1
  end if
 
  do
 
    call perm_inc ( circuit, ipos, nnode )

    if ( ipos <= 0 .or. circuit(1) /= 1 ) then
      iset = 0
      circuit(1:nnode) = 0
      return
    end if
!
!  Check whether the entries of CIRCUIT actually form a circuit.
!  If we find a break in the circuit, store that location in IPOS
!  and move on to try the next permutation.
!
    ipos = 0
    do i = 1, nnode-1
      if ( adj(circuit(i),circuit(i+1)) == 0 ) then
        ipos = i
        exit
      end if
    end do

    if ( ipos /= 0 ) then
      cycle
    end if
!
!  If the circuit connects all the nodes, we only have to check whether
!  the last node connects back to the first one.
!
    if ( adj(circuit(nnode),circuit(1)) /= 0 ) then
      exit
    end if

    ipos = nnode - 1

  end do

  iset = iset + 1

  return
end
subroutine digraph_adj_ham_path_next_brute ( adj, lda, nnode, path, iset )

!*****************************************************************************80
!
!! DIGRAPH_ADJ_HAM_PATH_NEXT_BRUTE finds the next path in a digraph that visits all nodes.
!
!  Discussion:
!
!    The path is not required to be a circuit.  That is, there is no requirement
!    that there be an edge from the last node visited back to the first one.
!
!    This is a brute force algorithm, and not suitable for large problems.
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
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(LDA,NNODE), the adjacency information.
!    ADJ(I,J) is 1 if there is a direct link from node I to node J.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of ADJ, which must be
!    at least NNODE.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input/output, integer ( kind = 4 ) PATH(NNODE).
!
!    On input, if ISET = 0, then PATH is not presumed to contain any
!    information.  If ISET is nonzero, then PATH contains the
!    path computed on the previous call.
!
!    On output, PATH contains the path computed by this call.
!
!    Input/output, integer ( kind = 4 ) ISET.
!
!    On input, a 0 value means this is the first call for this
!    graph.  Any other value means this is a repeated call for more paths.
!
!    On output, a 0 value means that no more paths could be computed.
!    Otherwise, ISET is incremented by one on each call.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ipos
  integer ( kind = 4 ) iset
  integer ( kind = 4 ) path(nnode)
!
!  If ISET is 0, this is a starting call, and we set PATH
!  to the lexically first path to check.
!
!  Otherwise, we set PATH to the next permutation.
!
  if ( iset == 0 ) then
    ipos = 0
    path(1:nnode) = 0
  else
    ipos = nnode - 1
  end if
 
  do
 
    call perm_inc ( path, ipos, nnode )
 
    if ( ipos == 0 ) then
      iset = 0
      path(1:nnode) = 0
      return
    end if
!
!  Check whether the entries of PATH actually form a path.
!
    ipos = 0
    do i = 1, nnode-1
      if ( adj(path(i),path(i+1)) == 0 ) then
        ipos = i
        exit
      end if
    end do

    if ( ipos == 0 ) then
      exit
    end if

  end do 

  iset = iset + 1
 
  return
end
subroutine digraph_adj_is_edge_connected ( adj, lda, nnode, result )

!*****************************************************************************80
!
!! DIGRAPH_ADJ_IS_EDGE_CONNECTED determines if a digraph is edgewise connected.
!
!  Discussion:
!
!    A digraph is edgewise connected if from any edge it is possible to reach
!    any other edge.  An edgewise connected digraph may include isolated nodes.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(LDA,NNODE), the adjacency information.
!    ADJ(I,J) is 1 if there is a direct link from node I to node J.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of LDA, which must be
!    at least NNODE.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) RESULT.
!    0, the digraph is not edgewise connected.
!    1, the digraph is edgewise connected.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) found(nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) list(nnode)
  integer ( kind = 4 ) result
!
!  FOUND(I) is 1 if edge I has been reached.
!  LIST(I) contains a list of the nodes as they are reached.
!
  list(1:nnode) = 0
  found(1:nnode) = 0
!
!  Find an edge.
!
  ilo = 1
  ihi = 0

  do i = 1, nnode
    do j = 1, nnode

      if ( 0 < adj(i,j) ) then

        adj(i,j) = - adj(i,j)
        ihi = ihi + 1
        list(ihi) = i
        found(i) = 1

        if ( i /= j ) then
          ihi = ihi + 1
          list(ihi) = j
          found(j) = 1
        end if

        exit

      end if

    end do

    if ( 0 < ihi ) then
      exit
    end if

  end do
!
!  A digraph with NO edges is edgewise connected!
!
  if ( ihi == 0 ) then
    result = 1
    return
  end if
!
!  From the batch of edge nodes found last time, LIST(ILO:IHI),
!  look for unfound neighbors, and store their indices in LIST(JLO:JHI).
!
  do

    jlo = ihi + 1
    jhi = ihi

    do ii = ilo, ihi

      i = list(ii)

      do j = 1, nnode

        if ( 0 < adj(i,j) ) then

          adj(i,j) = - adj(i,j)

          if ( found(j) == 0 ) then
            jhi = jhi + 1
            list(jhi) = j
            found(j) = 1
          end if

        end if

      end do

    end do

    if ( jhi < jlo ) then
      exit
    end if

    ilo = jlo
    ihi = jhi

  end do
!
!  If any edges were unvisited, then the digraph is not edgewise connected.
!
  result = 1

  do i = 1, nnode
    do j = 1, nnode
      if ( 0 < adj(i,j) ) then
        result = 0
      end if
    end do
  end do
!
!  Restore the positive sign of ADJ.
!
  adj(1:nnode,1:nnode) = abs ( adj(1:nnode,1:nnode) )

  return
end
subroutine digraph_adj_is_eulerian ( adj, lda, nnode, result )

!*****************************************************************************80
!
!! DIGRAPH_ADJ_IS_EULERIAN determines if a digraph is Eulerian from its adjacency ma
!
!  Discussion:
!
!    A digraph is path-Eulerian if there exists a path through the digraph
!    which uses every edge once.
!
!    A digraph is circuit-Eulerian if there exists a path through the digraph
!    which uses every edge once, and which starts and ends on the same node.
!
!    Note that it is NOT necessary for the path or circuit to pass through
!    every node; simply that all the edges can be used exactly once to
!    make a connected path.  This means an Eulerian digraph can have isolated
!    nodes, for instance.
!
!    A digraph is path-Eulerian if and only if it is edge-connected, and 
!    for all but two nodes, the indegree and outdegree are equal, and
!    for those two nodes, the indegree and outdegree, if different, differ
!    by 1.
!
!    A digraph is circuit-Eulerian if and only if it is edge connected and
!    for every node the indegree equals the outdegree.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 October 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(LDA,NNODE), the adjacency information.
!    ADJ(I,J) is 1 if there is a direct link from node I to node J.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of LDA, which must be
!    at least NNODE.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) RESULT.
!    0, the digraph is not Eulerian.
!    1, the digraph is path-Eulerian.
!    2, the digraph is circuit-Eulerian.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indegree(nnode)
  integer ( kind = 4 ) ndiff
  integer ( kind = 4 ) outdegree(nnode)
  integer ( kind = 4 ) result
!
!  First check that the digraph is edgewise connected.
!
  call digraph_adj_is_edge_connected ( adj, lda, nnode, result )

  if ( result == 0 ) then
    return
  end if
!
!  Now look at node degree.
!
  call digraph_adj_degree ( adj, lda, nnode, indegree, outdegree )

  ndiff = 0

  do i = 1, nnode

    if ( indegree(i) /= outdegree(i) ) then

      ndiff = ndiff + 1

      if ( 2 < ndiff ) then
        result = 0
        return
      end if

      if ( 1 < abs ( indegree(i) - outdegree(i) ) ) then
        result = 0
        return
      end if

    end if

  end do

  if ( ndiff == 0 ) then
    result = 2
  else
    result = 1
  end if

  return
end
subroutine digraph_adj_is_strong_connected ( adj, lda, nnode, result )

!*****************************************************************************80
!
!! DIGRAPH_ADJ_IS_STRONG_CONNECTED determines if a digraph is strongly connected.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 November 1999
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(LDA,NNODE), the adjacency information.
!    ADJ(I,J) is 1 if there is a direct link from node I to node J.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of ADJ.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) RESULT,
!    0, the digraph is not strongly connected;
!    1, the digraph is strongly connected.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) dad(nnode)
  integer ( kind = 4 ) iorder
  integer ( kind = 4 ) lowlink(nnode)
  integer ( kind = 4 ) mark(nnode)
  integer ( kind = 4 ) ncomp
  integer ( kind = 4 ) nstack
  integer ( kind = 4 ) order(nnode)
  integer ( kind = 4 ) point(nnode)
  integer ( kind = 4 ) result
  integer ( kind = 4 ) stack(nnode)
  integer ( kind = 4 ) v
  integer ( kind = 4 ) w
  integer ( kind = 4 ) x
!
!  Initialization.
!
  dad(1:nnode) = 0
  order(1:nnode) = 0
  lowlink(1:nnode) = 0
  mark(1:nnode) = 0
  point(1:nnode) = 0

  iorder = 0
  nstack = 0
  ncomp = 0
!
!  Select any node V not stored in the stack, that is, with MARK(V) = 0.
!
  do

    v = 0

    do

      v = v + 1

      if ( nnode < v ) then

        adj(1:nnode,1:nnode) = abs ( adj(1:nnode,1:nnode) )

        if ( 1 < ncomp ) then
          result = 0
        else
          result = 1
        end if

        return
      end if

      if ( mark(v) /= 1 ) then
        exit
      end if

    end do

    iorder = iorder + 1

    order(v) = iorder
    lowlink(v) = iorder
    mark(v) = 1

    nstack = nstack + 1
    stack(nstack) = v
    point(v) = 1

30  continue
!
!  Consider each node W.
!
    do w = 1, nnode
!
!  Is there an edge (V,W) and has it not been examined yet?
!
      if ( 0 < adj(v,w) ) then

        adj(v,w) = - adj(v,w)
!
!  Is the node on the other end of the edge undiscovered yet?
!
        if ( mark(w) == 0 ) then

          iorder = iorder + 1
          order(w) = iorder
          lowlink(w) = iorder
          dad(w) = v
          mark(w) = 1

          nstack = nstack + 1
          stack(nstack) = w
          point(w) = 1

          v = w

        else if ( mark(w) == 1 ) then

          if ( order(w) < order(v) .and. point(w) == 1 ) then
            lowlink(v) = min ( lowlink(v), order(w) )
          end if

        end if

        go to 30

      end if

    end do

    if ( lowlink(v) == order(v) ) then

      ncomp = ncomp + 1

      do

        if ( nstack <= 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'DIGRAPH_ADJ_IS_STRONG_CONNECTED - Fatal error!'
          write ( *, '(a)' ) '  Illegal stack reference.'
          stop
        end if

        x = stack(nstack)
        nstack = nstack - 1

        point(x) = 0

        if ( x == v ) then
          exit
        end if

      end do

    end if

    if ( dad(v) /= 0 ) then
      lowlink(dad(v)) = min ( lowlink(dad(v)), lowlink(v) )
      v = dad(v)
      go to 30
    end if

  end do

  return
end
subroutine digraph_adj_is_tournament ( adj, lda, nnode, result )

!*****************************************************************************80
!
!! DIGRAPH_ADJ_IS_TOURNAMENT determines if a digraph is a tournament.
!
!  Discussion:
!
!    A digraph is a tournament if every pair of distinct nodes is connected by
!    exactly one directed edge.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(LDA,NNODE), the adjacency information.
!    ADJ(I,J) is 1 if there is a direct link from node I to node J.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of LDA, which must 
!    be at least NNODE.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) RESULT.
!    0, the digraph is not a tournament.
!    1, the digraph is a tournament.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) result

  result = 0
!
!  No self links.
!
  do i = 1, nnode
    if ( adj(i,i) /= 0 ) then
      return
    end if
  end do
!
!  Distinct I and J must have exactly one connection.
!
  do i = 1, nnode
    do j = i+1, nnode
      if ( .not. ( adj(i,j) == 0 .and. adj(j,i) == 1 ) .and. &
           .not. ( adj(i,j) == 1 .and. adj(j,i) == 0 ) ) then
        return
      end if
    end do
  end do

  result = 1
 
  return
end
subroutine digraph_adj_is_transitive ( adj, lda, nnode, result )

!*****************************************************************************80
!
!! DIGRAPH_ADJ_IS_TRANSITIVE determines if a digraph is transitive.
!
!  Discussion:
!
!    A digraph is transitive if whenever there's a long way between two
!    nodes, there's an immediate way.  Formally:
!
!      ADJ(I,J) and ADJ(J,K) nonzero imply ADJ(I,K) nonzero.
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
!    Input, integer ( kind = 4 ) ADJ(LDA,NNODE), the adjacency information.
!    ADJ(I,J) is 1 if there is a direct link from node I to node J.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of LDA, which must be
!    at least NNODE.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) RESULT.
!    0, the digraph is not transitive.
!    1, the digraph is transitive.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) result

  result = 0

  do i = 1, nnode
    do j = 1, nnode
      if ( adj(i,j) /= 0 ) then
        do k = 1, nnode
          if ( adj(j,k) /= 0 ) then
            if ( adj(i,k) == 0 ) then
              return
            end if
          end if
        end do
      end if
    end do
  end do

  result = 1

  return
end
subroutine digraph_adj_is_weak_connected ( adj, lda, nnode, result )

!*****************************************************************************80
!
!! DIGRAPH_ADJ_IS_WEAK_CONNECTED determines if a digraph is weakly connected.
!
!  Discussion:
!
!    A digraph is weakly connected if the underlying graph is node connected.
!    In other words, if a graph is constructed from the digraph by replacing
!    every directed edge by an undirected edge, and the it is possible to
!    travel from any node to any other node, then the digraph is weakly
!    connected.
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
!    Input, integer ( kind = 4 ) ADJ(LDA,NNODE), the adjacency matrix for the digraph.  
!    ADJ(I,J) is nonzero if there is an edge from node I to node J.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of LDA, which must be
!    at least NNODE.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) RESULT.
!    0, the digraph is not weakly connected.
!    1, the digraph is weakly connected.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) result

  call graph_adj_is_node_connected ( adj, lda, nnode, result )

  return
end
subroutine digraph_adj_print ( adj, lda, nnode, title )

!*****************************************************************************80
!
!! DIGRAPH_ADJ_PRINT prints out an adjacency matrix for a digraph.
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
!    04 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(LDA,NNODE), the adjacency matrix of a digraph.  
!    ADJ(I,J) is 1 if there is a direct connection FROM node I TO node J,
!    and is 0 otherwise.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of ADJ, which must be
!    at least NNODE.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.  
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  character ( len = 80 ) string
  character ( len = * ) title

  if ( len_trim ( title ) /= 0 ) then
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

    write ( *, '(i2,2x,a)' ) i, string(1:jhi)

  end do

  return
end
subroutine digraph_adj_random ( lda, nnode, nedge, seed, adj )

!*****************************************************************************80
!
!! DIGRAPH_ADJ_RANDOM generates a random digraph.
!
!  Discussion:
!
!    A digraph is a directed graph.
!
!    The user specifies the number of nodes and edges in the digraph.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of LDA, which must be
!    at least NNODE.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) NEDGE, the number of edges, which must be between
!    0 and NNODE*(NNODE-1).
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
!    Output, integer ( kind = 4 ) ADJ(LDA,NNODE), the adjacency information.
!    ADJ(I,J) is 1 if there is a direct link from node I to node J.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) nnode
  integer ( kind = 4 ) nedge

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iwork(nedge)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) maxedge
  integer ( kind = 4 ) seed

  if ( nnode <= 0  ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DIGRAPH_ADJ_RANDOM - Fatal error!'
    write ( *, '(a,i8)' ) '  NNODE = ', nnode
    write ( *, '(a)' ) '  but NNODE must be at least 1.'
    stop
  end if

  maxedge = nnode * ( nnode - 1 )

  if ( nedge < 0 .or. maxedge < nedge ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DIGRAPH_RANDOM - Fatal error!'
    write ( *, '(a,i8)' ) '  NEDGE = ', nedge
    write ( *, '(a)' ) '  but NEDGE must be at least 0, and '
    write ( *, '(a,i8)' ) '  no more than ', maxedge
    stop
  end if

  adj(1:nnode,1:nnode) = 0
!
!  Pick a random NEDGE subset of NNODE*(NNODE-1).
!
  call ksub_random ( maxedge, nedge, seed, iwork )
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
subroutine digraph_adj_reduce ( adj, nnode )

!*****************************************************************************80
!
!! DIGRAPH_ADJ_REDUCE generates a transitive reduction of a digraph.
!
!  Discussion:
!
!    This routine is given an adjacency matrix B, which might be a
!    transitive closure of a graph G.
!
!    The transitive closure graph is generated from a graph G by the 
!    following procedure:
!
!      B(I,J) = 0 if node J cannot be reached from node I in graph G;
!               1 if node J can be reached from node I in graph G.
!
!    The purpose of this routine is to try to find the original, sparser
!    graph G which generated the given transitive closure graph.  Such a
!    graph G is known as a transitive reduction..  In general,
!    there is no unique solution.  In particular, any graph is a transitive
!    reduction of itself.  
!
!    Hence, the real task is to drop as many redundant edges as possible
!    from the given graph, arriving at a graph from which no more edges 
!    may be removed.
!
!  Method:
!
!    One way of explaining the algorithm is based on the adjacency matrix:
!
!    * Zero out the diagonals of the adjacency matrix.
!
!    * Consider row 1.  Any other row that can "reach" row 1 doesn't
!      need a 1 if row 1 has it.  So "subtract" all the 1's in row 1
!      from such rows.  We are done with row 1 and column 1.
!
!    * Repeat for the other rows.
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
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) ADJ(NNODE,NNODE).
!    On input, the adjacency matrix of the transitive closure graph H.
!    On output, the adjacency matrix of a transitive reduction graph G 
!    of the graph H.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
!
!  First discard those useless self-edges.
!
  do i = 1, nnode
    adj(i,i) = 0
  end do
!
!  If you can get from J to I and I to K, you don't need a direct
!  edge from J to K.
!
  do i = 1, nnode
    do j = 1, nnode
      if ( adj(j,i) /= 0 ) then
        do k = 1, nnode
          if ( adj(i,k) /= 0 ) then
            adj(j,k) = 0
          end if
        end do
      end if
    end do
  end do

  return
end
subroutine digraph_adj_to_digraph_arc ( adj, lda, nnode, maxedge, nedge, &
  inode, jnode )

!*****************************************************************************80
!
!! DIGRAPH_ADJ_TO_DIGRAPH_ARC converts a digraph from adjacency to arc list form.
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
!    Input, integer ( kind = 4 ) ADJ(LDA,NNODE), the adjacency information.
!    ADJ(I,J) is 1 if there is a direct link from node I to node J.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of LDA, which must be
!    at least NNODE.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) MAXEDGE, the maximum number of edges.
!
!    Output, integer ( kind = 4 ) NEDGE, the number of edges.
!
!    Output, integer ( kind = 4 ) INODE(MAXEDGE), JNODE(MAXEDGE), the arc list of the
!    digraph.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) maxedge
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) inode(maxedge)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jnode(maxedge)
  integer ( kind = 4 ) nedge

  nedge = 0

  inode(1:maxedge) = 0
  jnode(1:maxedge) = 0

  do j = 1, nnode
    do i = 1, nnode
      if ( adj(i,j) /= 0 ) then
        nedge = nedge + 1
        if ( nedge <= maxedge ) then
          inode(nedge) = i
          jnode(nedge) = j
        else
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'DIGRAPH_ADJ_TO_DIGRAPH_ARC - Fatal error!'
          write ( *, '(a)' ) '  MAXEDGE exceeded.'
          stop
        end if
      end if
    end do
  end do

  return
end
subroutine digraph_adj_to_digraph_inc ( adj, lda, nnode, maxarc, narc, inc )

!*****************************************************************************80
!
!! DIGRAPH_ADJ_TO_DIGRAPH_INC converts an adjacency digraph to an incidence digraph.
!
!  Discussion:
!
!    INC(node,arc) = 0 if NODE is not the beginning or end of ARC, or
!                       if ARC is a loop;
!                     1 if NODE is the beginning of ARC;
!                    -1 if NODE is the end of ARC.
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(LDA,NNODE), the adjacency matrix for the graph.  
!    ADJ(I,J) is nonzero if there is an edge from node I to node J.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of LDA, which must be
!    at least NNODE.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) MAXARC, the maximum number of arcs.
!
!    Output, integer ( kind = 4 ) NARC, the number of arcs.
!
!    Output, integer ( kind = 4 ) INC(LDA,MAXARC), the incidence matrix.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) maxarc
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) inc(lda,maxarc)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) narc

  narc = 0

  do j = 1, nnode
    do i = 1, nnode

      if ( i == j ) then

      else if ( adj(i,j) /= 0 ) then
        narc = narc + 1
        if ( narc <= maxarc ) then
          inc(i,narc) = 1
          inc(j,narc) = -1
        else
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'DIGRAPH_ADJ_TO_DIGRAPH_INC - Fatal error!'
          write ( *, '(a)' ) '  MAXARC exceeded.'
          stop
        end if
      end if
    end do
  end do

  return
end
subroutine digraph_adj_top_sort ( adj, lda, nnode, dad, visit, node_list )

!*****************************************************************************80
!
!! DIGRAPH_ADJ_TOP_SORT finds a reverse topological sorting of a directed acyclic graph.
!
!  Discussion:
!
!    The routine performs a depth first search of the DAG and returns:
!
!    * a list of the order in which the nodes were visited;
!    * a list of the parents of each node in the search trees;
!    * a list of the nodes, in a reverse topological order.
!
!    In a reverse topological sorting of the nodes of a directed
!    acyclic graph, nodes are listed "lowest" first.  That is,
!    if node A precedes node B in the list, then there may or may
!    not be an edge or indirect path from B to A, but there
!    is neither an edge or indirect path from A to B.
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
!    Robert Sedgewick,
!    Algorithms,
!    Addison Wesley, 1983, page 426.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(LDA,NNODE), the adjacency information.
!    ADJ(I,J) is 1 if there is a direct link from node I to node J.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of ADJ, which must be
!    at least NNODE.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) DAD(NNODE), the father array for the depth first
!    search trees.  DAD(I) = 0 means that node I is the root of 
!    one of the trees.  DAD(I) = J means that the search descended
!    from node J to node I.
!
!    Output, integer ( kind = 4 ) VISIT(NNODE), the order in which the nodes were
!    visited, from 1 to NNODE.
!
!    Output, integer ( kind = 4 ) NODE_LIST(NNODE), a list of the nodes, in reverse 
!    topological order.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) dad(nnode)
  integer ( kind = 4 ) daddy
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) maxstack
  integer ( kind = 4 ) nsort
  integer ( kind = 4 ) nstack
  integer ( kind = 4 ) node_list(nnode)
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) stack(2*(nnode-1))
  integer ( kind = 4 ) visit(nnode)

  dad(1:nnode) = 0
  maxstack = 2 * ( nnode - 1 )
  visit(1:nnode) = 0
  node_list(1:nnode) = 0

  rank = 0
  nsort = 0

  do i = 1, nnode
!
!  Find the next unused node and begin a new search tree.
!
    if ( visit(i) == 0 ) then

      daddy = i
      dad(daddy) = 0
      nstack = 0
!
!  Visit node DAD.
!
10    continue

      rank = rank + 1
      visit(daddy) = rank
      j = 0
!
!  Consider visiting node J from node DAD.
!
20    continue

      j = j + 1
!
!  If J is a reasonable value, adjacent to DAD, and unvisited,
!  then put DAD into the stack, make J the new value of DAD,
!  and go to 10.
!
      if ( j <= nnode ) then

        if ( adj(daddy,j) /= 0 .and. visit(j) == 0 ) then

          if ( nstack+2 <= maxstack ) then
            dad(j) = daddy
            stack(nstack+1) = daddy
            stack(nstack+2) = j
            nstack = nstack + 2
            daddy = j
            go to 10
          else
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'DIGRAPH_ADJ_TOP_SORT - Fatal error!'
            write ( *, '(a)' ) '  Out of stack space.'
            stop
          end if
!
!  If J is not suitable for a visit, get the next value of J.
!
        else

          go to 20

        end if
!
!  If no more neighbors to consider, back up one node.
!
      else if ( 2 <= nstack ) then

        nsort = nsort + 1
        node_list(nsort) = daddy

        daddy = stack(nstack-1)
        j = stack(nstack)
        nstack = nstack - 2
        go to 20
!
!  If no more nodes to consider in this tree, bail out.
!
      else

        nsort = nsort + 1
        node_list(nsort) = daddy

        nstack = 0

      end if

    end if

  end do

  return
end
subroutine digraph_adj_tournament_random ( lda, nnode, seed, adj )

!*****************************************************************************80
!
!! DIGRAPH_ADJ_TOURNAMENT_RANDOM generates a random tournament digraph.
!
!  Discussion:
!
!    Definition: A tournament is a directed graph in which every pair 
!    of nodes are joined by exactly one directed edge.
!
!    The user specifies the number of nodes in the digraph.  The number of
!    edges will be (NNODE*(NNODE-1))/2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of LDA, which must be
!    at least NNODE.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
!    Output, integer ( kind = 4 ) ADJ(LDA,NNODE), the adjacency information.
!    ADJ(I,J) is 1 if there is a direct link from node I to node J.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed

  if ( nnode <= 0  ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DIGRAPH_ADJ_TOURNAMENT_RANDOM - Fatal error!'
    write ( *, '(a,i8)' ) '  NNODE = ', nnode
    write ( *, '(a)' ) '  but NNODE must be at least 1.'
    stop
  end if

  adj(1:nnode,1:nnode) = 0

  do i = 1, nnode
    do j = i+1, nnode

      k = i4_uniform ( 1, 2, seed )

      if ( k == 1 ) then
        adj(i,j) = 1
      else
        adj(j,i) = 1
      end if

    end do
  end do

  return
end
subroutine digraph_arc_degree ( nnode, nedge, inode, jnode, indegree, &
  outdegree )

!*****************************************************************************80
!
!! DIGRAPH_ARC_DEGREE determines the degree of the nodes of a digraph.
!
!  Discussion:
!
!    Definition: The degree of a node is the number of edges that 
!    include the node.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
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
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) NEDGE, the number of edges.
!
!    Input, integer ( kind = 4 ) INODE(NEDGE), JNODE(NEDGE), the pairs of nodes
!    that form the edges.
!
!    Output, integer ( kind = 4 ) INDEGREE(NNODE), OUTDEGREE(NNODE), the indegree 
!    and outdegree of each node, that is, the number of edges that end 
!    with the node, and the number of edges that begin with it.
!
  implicit none

  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) i
  integer ( kind = 4 ) indegree(nnode)
  integer ( kind = 4 ) inode(nedge)
  integer ( kind = 4 ) jnode(nedge)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) outdegree(nnode)

  indegree(1:nnode) = 0
  outdegree(1:nnode) = 0

  do i = 1, nedge

    n = inode(i)
    if ( 1 <= n .and. n <= nnode ) then
      outdegree(n) = outdegree(n) + 1
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DIGRAPH_ARC_DEGREE - Fatal error!'
      write ( *, '(a,i8)' ) '  Out-of-range node value = ', n
      stop
    end if

    n = jnode(i)
    if ( 1 <= n .and. n <= nnode ) then
      indegree(n) = indegree(n) + 1
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DIGRAPH_ARC_DEGREE - Fatal error!'
      write ( *, '(a,i8)' ) '  Out-of-range node value = ', n
      stop
    end if

  end do

  return
end
subroutine digraph_arc_edge_sort ( nedge, inode, jnode )

!*****************************************************************************80
!
!! DIGRAPH_ARC_EDGE_SORT sorts the edge array of a graph.
!
!  Discussion:
!
!    The edges are sorted in dictionary order.
!
!  Example:
!
!    Input:
!
!      INODE  JNODE
!
!        3      2
!        2      4
!        4      3
!        2      1
!        1      4
!
!    Output:
!
!      INODE  JNODE
!
!        1      4
!        2      1
!        2      4
!        3      2
!        4      3
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
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEDGE, the number of edges.
!
!    Input/output, integer ( kind = 4 ) INODE(NEDGE), JNODE(NEDGE), the edge array.
!    The I-th edge goes from node INODE(I) to node JNODE(I).
!    On output, the INODE and JNODE arrays have been sorted as described.
!
  implicit none

  integer ( kind = 4 ) nedge

  integer ( kind = 4 ) iedge
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) inode(nedge)
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) jedge
  integer ( kind = 4 ) jnode(nedge)

  if ( nedge <= 1 ) then
    return
  end if
!
!  Sort the edges using an external heap sort.
!
  iedge = 0
  jedge = 0
  indx = 0
  isgn = 0

  do

    call sort_heap_external ( nedge, indx, iedge, jedge, isgn )
!
!  Interchange edges IEDGE and JEDGE.
!
    if ( 0 < indx ) then

      call i4_swap ( inode(iedge), inode(jedge) )
      call i4_swap ( jnode(iedge), jnode(jedge) )
!
!  Compare edges IEDGE and JEDGE.
!
    else if ( indx < 0 ) then

      if ( ( inode(iedge) < inode(jedge) ) .or. &
        ( inode(iedge) == inode(jedge) .and. &
          jnode(iedge) < jnode(jedge) ) ) then
        isgn = -1
      else
        isgn = +1
      end if

    else if ( indx == 0 ) then

      exit

    end if

  end do
 
  return
end
subroutine digraph_arc_euler_circ_cand ( nedge, inode, jnode, circuit, k, &
  nstack, stack, maxstack, ncan, iwork, lwork )

!*****************************************************************************80
!
!! DIGRAPH_ARC_EULER_CIRC_CAND finds candidates for the K-th edge of an Euler circuit.
!
!  Discussion:
!
!    This routine is used in conjunction with I4VEC_BACKTRACK, which directs the 
!    search for a complete Euler circuit.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 August 2000
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
!    Input, integer ( kind = 4 ) NEDGE, the number of edges in the digraph.
!
!    Input, integer ( kind = 4 ) INODE(NEDGE), JNODE(NEDGE), the edge array of the digraph.
!    The I-th edge extends from node INODE(I) to JNODE(I).
!
!    Input, integer ( kind = 4 ) CIRCUIT(NEDGE), CIRCUIT(I) is the I-th edge in the circuit.  
!    A full circuit will have NEDGE edges, but on input we only have K-1.
!
!    Input, integer ( kind = 4 ) K, the index of the next edge to be determined in circuit.
!
!    Input/output, integer ( kind = 4 ) NSTACK, the current length of the stack.
!
!    Input, integer ( kind = 4 ) STACK(MAXSTACK), as yet unused candidates for positions
!    1 to K-1.
!
!    Input, integer ( kind = 4 ) MAXSTACK, the dimension of STACK.
!
!    Workspace, integer IWORK(NEDGE).
!
!    Workspace, logical LWORK(NEDGE).
!
  implicit none

  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) maxstack

  integer ( kind = 4 ) circuit(nedge)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) inode(nedge)
  integer ( kind = 4 ) it
  integer ( kind = 4 ) iwork(nedge)
  integer ( kind = 4 ) jnode(nedge)
  integer ( kind = 4 ) k
  logical lwork(nedge)
  integer ( kind = 4 ) ncan(nedge)
  integer ( kind = 4 ) nstack
  integer ( kind = 4 ) stack(maxstack)

  ncan(k) = 0

  if ( k == 1 ) then
    iwork(1) = jnode(1)
    stack(1) = 1
    nstack = 1
    ncan(k) = 1
    return
  end if
 
  if ( 2 < k ) then
    iwork(k-1) = inode(circuit(k-1)) + jnode(circuit(k-1)) - iwork(k-2)
  end if
 
  it = iwork(k-1)
 
  do i = 1, nedge
    lwork(i) = it == inode(i)
  end do
 
  lwork(circuit(1:k-1)) = .false.
  
  do i = 1, nedge
    if ( lwork(i) ) then
      if ( maxstack <= nstack ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'DIGRAPH_ARC_EULER_CIRC_CAND - Fatal error!'
        write ( *, '(a)' ) '  Stack size exceeded.'
        stop
      end if
      nstack = nstack + 1
      stack(nstack) = i
      ncan(k) = ncan(k) + 1
    end if
  end do
 
  return
end
subroutine digraph_arc_euler_circ_next ( nedge, inode, jnode, circuit, stack, &
  maxstack, ncan, more )

!*****************************************************************************80
!
!! DIGRAPH_ARC_EULER_CIRC_NEXT returns the next Euler circuit for a digraph.
!
!  Discussion:
!
!    The routine produces all the Euler circuits of a digraph, one at a time.
!
!    Definition: An Euler circuit of a digraph is a path starting at some node, 
!    using all the edges of the digraph exactly once, and returning
!    to the starting node.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 August 2000
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
!    Input, integer ( kind = 4 ) NEDGE, the number of edges in the digraph.
!
!    Input, integer ( kind = 4 ) INODE(NEDGE), JNODE(NEDGE), the edge array of the digraph.
!    The I-th edge extends from node INODE(I) to JNODE(I).
!
!    Output, integer ( kind = 4 ) CIRCUIT(NEDGE).  If MORE = TRUE on output, then IARRAY
!    contains the edges, in order, that constitute this circuit.
!
!    Workspace, integer STACK(MAXSTACK).  
!
!    Input, integer ( kind = 4 ) MAXSTACK, the dimension of STACK.
!
!    Input/output, logical MORE.
!    On first call, set MORE to .FALSE, and do not alter it after.
!    On return, MORE is TRUE if another circuit has been returned in
!    IARRAY, and FALSE if there are no more circuits.
!
  implicit none

  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) maxstack

  integer ( kind = 4 ) circuit(nedge)
  integer ( kind = 4 ) inode(nedge)
  integer ( kind = 4 ), save :: indx = 0
  integer ( kind = 4 ) iwork(nedge)
  integer ( kind = 4 ) jnode(nedge)
  integer ( kind = 4 ), save :: k = 0
  logical lwork(nedge)
  logical more
  integer ( kind = 4 ) ncan(nedge)
  integer ( kind = 4 ), save :: nstack = 0
  integer ( kind = 4 ) stack(maxstack)

  if ( .not. more ) then
    indx = 0
    k = 0
    more = .true.
    nstack = 0
  end if
 
  do
 
    call i4vec_backtrack ( nedge, circuit, indx, k, nstack, stack, maxstack, &
      ncan )
 
    if ( indx == 1 ) then

      exit

    else if ( indx == 2 ) then

      call digraph_arc_euler_circ_cand ( nedge, inode, jnode, circuit, k, &
        nstack, stack, maxstack, ncan, iwork, lwork )

    else

      more = .false.
      exit

    end if

  end do
 
  return
end
subroutine digraph_arc_example_cycler ( maxedge, nedge, inode, jnode )

!*****************************************************************************80
!
!! DIGRAPH_ARC_EXAMPLE_CYCLER sets up the arc list information for the cycler digraph.
!
!  Diagram:
!  
!           A
!           |
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
!    26 October 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) MAXEDGE, the maximum number of edges.
!
!    Output, integer ( kind = 4 ) NEDGE, the number of edges.
!
!    Output, integer ( kind = 4 ) INODE(MAXEDGE), JNODE(MAXEDGE), the arc list
!    for the digraph.
!
  implicit none

  integer ( kind = 4 ) maxedge

  integer ( kind = 4 ) inode(maxedge)
  integer ( kind = 4 ) jnode(maxedge)
  integer ( kind = 4 ) nedge

  nedge = 16

  if ( maxedge < nedge ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DIGRAPH_ARC_EXAMPLE_CYCLER - Fatal error!'
    write ( *, '(a)' ) '  MAXEDGE is too small.'
    stop
  end if

  inode(1) = 1
  jnode(1) = 3

  inode(2) = 1
  jnode(2) = 5

  inode(3) = 2
  jnode(3) = 6

  inode(4) = 2
  jnode(4) = 8

  inode(5) = 3
  jnode(5) = 4

  inode(6) = 3
  jnode(6) = 6

  inode(7) = 3
  jnode(7) = 7

  inode(8) = 4
  jnode(8) = 3

  inode(9) = 5
  jnode(9) = 2

  inode(10) = 6
  jnode(10) = 4

  inode(11) = 6
  jnode(11) = 8

  inode(12) = 7
  jnode(12) = 7

  inode(13) = 7
  jnode(13) = 9

  inode(14) = 8
  jnode(14) = 1

  inode(15) = 9
  jnode(15) = 5

  inode(16) = 9
  jnode(16) = 7

  return
end
subroutine digraph_arc_is_eulerian ( nnode, nedge, inode, jnode, indegree, &
  outdegree, result )

!*****************************************************************************80
!
!! DIGRAPH_ARC_IS_EULERIAN determines if a digraph is Eulerian from its edge list.
!
!  Discussion:
!
!    A digraph is Eulerian if there exists a circuit through the graph
!    which uses every edge once.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
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
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) NEDGE, the number of edges.
!
!    Input, integer ( kind = 4 ) INODE(NEDGE), JNODE(NEDGE), the pairs of nodes
!    that form the edges.
!
!    Output, integer ( kind = 4 ) INDEGREE(NNODE), OUTDEGREE(NODE), the indegree and
!    outdegree of each node, that is, the number of edges that end with
!    the node, and that begin the node.
!
!    Output, integer ( kind = 4 ) RESULT.
!    0, the digraph is not Eulerian.
!    1, the digraph is Eulerian, but the starting and ending nodes differ.
!    2, the digraph is Eulerian, and there is a closed Euler circuit.
!
  implicit none

  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) i
  integer ( kind = 4 ) indegree(nnode)
  integer ( kind = 4 ) inode(nedge)
  integer ( kind = 4 ) jnode(nedge)
  integer ( kind = 4 ) n_minus
  integer ( kind = 4 ) n_plus
  integer ( kind = 4 ) outdegree(nnode)
  integer ( kind = 4 ) result

  call digraph_arc_degree ( nnode, nedge, inode, jnode, indegree, outdegree )

  n_plus = 0
  n_minus = 0

  do i = 1, nnode

    if ( indegree(i) == outdegree(i) ) then

    else if ( n_plus == 0 .and. indegree(i) == outdegree(i) + 1 ) then
      n_plus = 1
    else if ( n_minus == 0 .and. indegree(i) == outdegree(i) - 1 ) then
      n_minus = 1
    else
      result = 0
      return
    end if

  end do

  if ( n_plus == 0 .and. n_minus == 0 ) then
    result = 2
  else if ( n_plus == 1 .and. n_minus == 1 ) then
    result = 1
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DIGRAPH_ARC_IS_EULERIAN - Fatal error!'
    write ( *, '(a)' ) '  The algorithm failed.'
    stop
  end if

  return
end
subroutine digraph_arc_print ( nedge, inode, jnode, title )

!*****************************************************************************80
!
!! DIGRAPH_ARC_PRINT prints out a digraph from an edge list.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
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
subroutine digraph_arc_to_digraph_adj ( nedge, inode, jnode, adj, lda, nnode )

!*****************************************************************************80
!
!! DIGRAPH_ARC_TO_DIGRAPH_ADJ converts an arc list digraph to an adjacency digraph.
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
!    Input, integer ( kind = 4 ) NEDGE, the number of edges.
!
!    Input, integer ( kind = 4 ) INODE(NEDGE), JNODE(NEDGE), the edge array.
!    The I-th edge connects nodes INODE(I) and JNODE(I).
!
!    Output, integer ( kind = 4 ) ADJ(LDA,NNODE), the adjacency information.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of ADJ.
!
!    Output, integer ( kind = 4 ) NNODE, the number of nodes.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) nedge

  integer ( kind = 4 ) adj(lda,*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) inode(nedge)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jnode(nedge)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) mnode
  integer ( kind = 4 ) nnode
!
!  Determine the number of nodes.
!
  call graph_arc_node_count ( nedge, inode, jnode, mnode, nnode )

  if ( lda < nnode ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DIGRAPH_ARC_TO_DIGRAPH_ADJ - Fatal error!'
    write ( *, '(a)' ) '  Number of nodes exceeds LDA.'
    stop
  end if

  adj(1:nnode,1:nnode) = 0

  do k = 1, nedge
    i = inode(k)
    j = jnode(k)
    adj(i,j) = 1
  end do

  return
end
subroutine digraph_arc_to_digraph_star ( nnode, nedge, inode, jnode, arcfir, &
  fwdarc )

!*****************************************************************************80
!
!! DIGRAPH_ARC_TO_DIGRAPH_STAR sets up the forward star representation of a digraph.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
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
!    Input, integer ( kind = 4 ) INODE(NEDGE), JNODE(NEDGE); the I-th edge
!    extends from node INODE(I) to JNODE(I).
!
!    Output, integer ( kind = 4 ) ARCFIR(NNODE+1); ARCFIR(I) is the number of the first
!    edge starting at node I in the forward star representation.
!
!    Output, integer ( kind = 4 ) FWDARC(NEDGE); FWDARC(I) is the ending node of
!    the I-th edge in the forward star representation.
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
!  Set up the forward star representation.
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
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
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
subroutine digraph_dist_print ( dist, lda, nnode, title )

!*****************************************************************************80
!
!! DIGRAPH_DIST_PRINT prints the distance matrix defining a digraph.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 July 2000
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
!    Input, integer ( kind = 4 ) LDA, the leading dimension of DIST, which must be at
!    least NNODE.
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
subroutine digraph_inc_print ( lda, nnode, narc, inc, title )

!*****************************************************************************80
!
!! DIGRAPH_INC_PRINT prints the incidence matrix of a digraph.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the array.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) NARC, the number of arcs.
!
!    Input, integer ( kind = 4 ) INC(LDA,NARC), the NNODE by NARC incidence matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) narc

  integer ( kind = 4 ) i
  integer ( kind = 4 ) inc(lda,narc)
  integer ( kind = 4 ) nnode
  character ( len = * ) title

  if ( len_trim ( title ) /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  do i = 1, nnode
    write ( *, '(20i3)' ) inc(i,1:narc)
  end do

  return
end
subroutine edge_add_nodes ( edge, max_edge, num_edge, iface, n1, n2, ierror )

!*****************************************************************************80
!
!! EDGE_ADD_NODES adds the edge defined by two nodes to the edge list.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) EDGE(4,MAX_EDGE), edge information.
!    EDGE(1,I) is the starting node of edge I;
!    EDGE(2,I) is the ending node of edge I;
!    EDGE(3,I) is the positive face;
!    EDGE(4,I) is the negative face, if any.
!
!    Input, integer ( kind = 4 ) MAX_EDGE, the maximum number of edges.
!
!    Input/output, integer ( kind = 4 ) NUM_EDGE, the number of edges.
!
!    Input, integer ( kind = 4 ) IFACE, the face to which the nodes belong.
!
!    Input, integer ( kind = 4 ) N1, N2, two nodes which form an edge.
!
!    Output, integer ( kind = 4 ) IERROR, error flag, 0 = no error, nonzero = error.
!
  implicit none

  integer ( kind = 4 ) max_edge

  integer ( kind = 4 ) edge(4,max_edge)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iface
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) num_edge

  if ( num_edge < max_edge ) then
    num_edge = num_edge + 1
    edge(1,num_edge) = n1
    edge(2,num_edge) = n2
    edge(3,num_edge) = iface
    edge(4,num_edge) = 0
    ierror = 0
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EDGE_ADD_NODES - Fatal error!'
    write ( *, '(a,i8)' ) '  Exceeding MAX_EDGE = ', max_edge
    ierror = 1
  end if

  return
end
subroutine edge_bound ( edge, max_edge, num_edge )

!*****************************************************************************80
!
!! EDGE_BOUND reports the edges which are part of the boundary.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) EDGE(4,MAX_EDGE), edge information.
!    EDGE(1,I) is the starting node of edge I;
!    EDGE(2,I) is the ending node of edge I;
!    EDGE(3,I) is the positive face;
!    EDGE(4,I) is the negative face, if any.
!
!    Input, integer ( kind = 4 ) MAX_EDGE, the maximum number of edges.
!
!    Input, integer ( kind = 4 ) NUM_EDGE, the number of edges.
!
  implicit none

  integer ( kind = 4 ) max_edge

  integer ( kind = 4 ) edge(4,max_edge)
  integer ( kind = 4 ) iedge
  integer ( kind = 4 ) num_bound
  integer ( kind = 4 ) num_edge

  num_bound = 0

  do iedge = 1, num_edge
    if ( ( edge(3,iedge) /= 0 .and. edge(4,iedge) == 0 ) .or. &
         ( edge(3,iedge) == 0 .and. edge(4,iedge) /= 0 ) ) then
      num_bound = num_bound + 1
    end if
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'EDGE_BOUND'
  write ( *, '(a,i8)' ) '  Number of boundary edges = ', num_bound

  return
end
subroutine edge_match_face ( edge, max_edge, num_edge, facelist, n, index )

!*****************************************************************************80
!
!! EDGE_MATCH_FACE seeks an edge common to a face and the edge list.
!
!  Discussion:
!
!    If a common edge is found, then the information in the face node
!    list is adjusted so that the first two entries correspond to the
!    matching edge in EDGE, but in reverse order.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) EDGE(4,MAX_EDGE), edge information.
!    EDGE(1,I) is the starting node of edge I;
!    EDGE(2,I) is the ending node of edge I;
!    EDGE(3,I) is the positive face;
!    EDGE(4,I) is the negative face, if any.
!
!    Input, integer ( kind = 4 ) MAX_EDGE, the maximum number of edges.
!
!    Input, integer ( kind = 4 ) NUM_EDGE, the number of edges.
!
!    Input/output, integer ( kind = 4 ) FACELIST(N), the list of nodes making a face.
!
!    Input, integer ( kind = 4 ) N, the number of nodes in the face.
!
!    Output, integer ( kind = 4 ) INDEX, the results of the search.
!    0, there is no edge common to the face and the EDGE array.
!    nonzero, edge INDEX is common to the face and the EDGE array.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) max_edge

  integer ( kind = 4 ) edge(4,max_edge)
  integer ( kind = 4 ) facelist(n)
  integer ( kind = 4 ) iedge
  integer ( kind = 4 ) index
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jp1
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) num_edge

  index = 0

  if ( n <= 0 ) then
    return
  end if

  if ( num_edge <= 0 ) then
    return
  end if

  do j = 1, n

    if ( j == n ) then
      jp1 = 1
    else
      jp1 = j + 1
    end if

    n1 = facelist(j)
    n2 = facelist(jp1)

    do iedge = 1, num_edge

      if ( edge(1,iedge) == n2 .and. edge(2,iedge) == n1 ) then

        call i4vec_rotate ( n, 1 - j, facelist )

        index = iedge
        return

      else if ( edge(1,iedge) == n1 .and. edge(2,iedge) == n2 ) then

        call i4vec_rotate ( n, n - jp1, facelist )

        call i4vec_reverse ( n, facelist )

        index = iedge
        return

      end if

    end do
   
  end do

  return
end
subroutine edge_match_nodes ( edge, max_edge, num_edge, n1, n2, iedge )

!*****************************************************************************80
!
!! EDGE_MATCH_NODES seeks an edge of the form (N1,N2) or (N2,N1) in EDGE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) EDGE(4,MAX_EDGE), edge information.
!    EDGE(1,I) is the starting node of edge I;
!    EDGE(2,I) is the ending node of edge I;
!    EDGE(3,I) is the positive face;
!    EDGE(4,I) is the negative face, if any.
!
!    Input, integer ( kind = 4 ) MAX_EDGE, the maximum number of edges.
!
!    Input, integer ( kind = 4 ) NUM_EDGE, the number of edges.
!
!    Input, integer ( kind = 4 ) N1, N2, two nodes that form an edge.
!
!    Output, integer ( kind = 4 ) IEDGE, the results of the search.
!    0, no matching edge was found.
!    nonzero, edge IEDGE of the EDGE array matches (N1,N2) or (N2,N1).
!
  implicit none

  integer ( kind = 4 ) max_edge

  integer ( kind = 4 ) edge(4,max_edge)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iedge
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) num_edge

  iedge = 0
  do i = 1, num_edge

    if ( ( n1 == edge(1,i) .and. n2 == edge(2,i) ) .or. &
         ( n2 == edge(1,i) .and. n1 == edge(2,i) ) ) then
      iedge = i
      return
    end if

  end do

  return
end
subroutine edges_to_ps ( plotxmin2, plotymin2, alpha, iunit, inode, jnode, &
  nedge, nnode, x, y, xmin, ymin )

!*****************************************************************************80
!
!! EDGES_TO_PS writes subplot edges to a PostScript file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PLOTXMIN2, PLOTYMIN2, the Postscript origin.
!
!    Input, real ( kind = 8 ) ALPHA, the physical-to-Postscript scale factor.
!
!    Input, integer ( kind = 4 ) IUNIT, the output FORTRAN unit.
!
!    Input, integer ( kind = 4 ) INODE(NEDGE), JNODE(NEDGE), the edge array.
!    The I-th edge connects nodes INODE(I) and JNODE(I).
!
!    Input, integer ( kind = 4 ) NEDGE, the number of edges.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, real ( kind = 8 ) X(NNODE), Y(NNODE), the X and Y components
!    of points.
!
!    Input, real ( kind = 8 ) XMIN, YMIN, the physical origin.
!
  implicit none

  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) nnode

  real ( kind = 8 ) alpha
  integer ( kind = 4 ) i
  integer ( kind = 4 ) inode(nedge)
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) jnode(nedge)
  integer ( kind = 4 ) node
  integer ( kind = 4 ) plotxmin2
  integer ( kind = 4 ) plotymin2
  integer ( kind = 4 ) px1
  integer ( kind = 4 ) px2
  integer ( kind = 4 ) py1
  integer ( kind = 4 ) py2
  real ( kind = 8 ) x(nnode)
  real ( kind = 8 ) xmin
  real ( kind = 8 ) y(nnode)
  real ( kind = 8 ) ymin
!
!  Draw lines.
!
  do i = 1, nedge

    node = inode(i)
    px1 = plotxmin2 + nint ( alpha * ( x(node) - xmin ) )
    py1 = plotymin2 + nint ( alpha * ( y(node) - ymin ) )

    node = jnode(i)
    px2 = plotxmin2 + nint ( alpha * ( x(node) - xmin ) )
    py2 = plotymin2 + nint ( alpha * ( y(node) - ymin ) )

    write ( iunit, '(2i4,a,2i4,a)' ) px1, py1, ' moveto ', px2, py2, &
      ' lineto stroke'

  end do

  return
end
subroutine elmhes ( nm, n, low, igh, a, ind )

!*****************************************************************************80
!
!! ELMHES transforms a real general matrix to upper Hessenberg form.
!
!  Discussion:
!
!    Given a real general matrix, this subroutine reduces a submatrix
!    situated in rows and columns LOW through IGH to upper Hessenberg
!    form by stabilized elementary similarity transformations.
!
!  Reference:
!
!    Martin, James Wilkinson,
!    ELMHES,
!    Numerische Mathematik,
!    Volume 12, pages 349-368, 1968.
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow, 
!    Y Ikebe, V Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NM, the leading dimension of the array A.
!    NM must be at least N.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) LOW, IGH, are determined by the balancing routine
!    BALANC.  If BALANC has not been used, set LOW = 1, IGH = N.
!
!    Input/output, real ( kind = 8 ) A(NM,N).  On input, the matrix to be
!    reduced.  On output, the Hessenberg matrix.  The multipliers
!    which were used in the reduction are stored in the
!    remaining triangle under the Hessenberg matrix.
!
!    Output, integer ( kind = 4 ) IND(N), contains information on the rows and columns
!    interchanged in the reduction.  Only elements LOW through IGH are used.
!
  implicit none

  integer ( kind = 4 ) igh
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nm

  real ( kind = 8 ) a(nm,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ind(igh)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) la
  integer ( kind = 4 ) low
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm1
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  la = igh - 1

  do m = low + 1, la

    mm1 = m - 1
    x = 0.0D+00
    i = m

    do j = m, igh
      if ( abs ( x ) < abs ( a(j,mm1) ) ) then
        x = a(j,mm1)
        i = j
      end if
    end do

    ind(m) = i
!
!  Interchange rows and columns of the matrix.
!
    if ( i /= m ) then

      do j = mm1, n
        call r8_swap ( a(i,j), a(m,j) )
      end do

      do j = 1, igh
        call r8_swap ( a(j,i), a(j,m) )
      end do

    end if

    if ( x /= 0.0D+00 ) then

      do i = m+1, igh

        y = a(i,mm1)

        if ( y /= 0.0D+00 ) then

          y = y / x
          a(i,mm1) = y

          do j = m, n
            a(i,j) = a(i,j) - y * a(m,j)
          end do

          do j = 1, igh
            a(j,m) = a(j,m) + y * a(j,i)
          end do

        end if

      end do

    end if

  end do

  return
end
subroutine face_check ( edge, face, face_object, face_order, face_rank, &
  face_tier, max_edge, max_order, num_edge, num_face, num_object )

!*****************************************************************************80
!
!! FACE_CHECK checks and analyzes a set of faces.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) EDGE(4,MAX_EDGE), edge information.
!    EDGE(1,I) is the starting node of edge I;
!    EDGE(2,I) is the ending node of edge I;
!    EDGE(3,I) is the positive face;
!    EDGE(4,I) is the negative face, or 0 if the edge is used once.
!
!    Input, integer ( kind = 4 ) FACE(MAX_ORDER,NUM_FACE), describes the faces.
!    FACE(I,J) is the index of the I-th node in face J.  It is best
!    if the nodes of all faces are listed in counterclockwise order.
!
!    Output, integer ( kind = 4 ) FACE_OBJECT(NUM_FACE), describes the objects.
!    FACE_OBJECT(I) is the index of the edge-connected "object" to 
!    which face I belongs.
!
!    Input, integer ( kind = 4 ) FACE_ORDER(NUM_FACE), is the number of nodes
!    making up each face.
!
!    Output, integer ( kind = 4 ) FACE_RANK(NUM_FACE), is an ordered list of faces.
!    FACE_RANK(1) is the index of the face in the first tier of the 
!    first object, followed by second tier faces, and so on until
!    object one is complete.  Object two follows, and so on.
!
!    Output, integer ( kind = 4 ) FACE_TIER(NUM_FACE).  FACE_TIER(I) is the "tier"
!    of face I in its object.  The seed of the object has tier 1,
!    the neighbors of the seed have tier 2, and so on.
!
!    Input, integer ( kind = 4 ) MAX_EDGE, the maximum number of edges.
!
!    Input, integer ( kind = 4 ) MAX_ORDER, is the maximum number of nodes that can
!    make up a face, required to dimension FACE.
!
!    Output, integer ( kind = 4 ) NUM_EDGE, the number of edges.
!
!    Input, integer ( kind = 4 ) NUM_FACE, the number of faces.
!
!    Output, integer ( kind = 4 ) NUM_OBJECT, the number of objects.
!
  implicit none

  integer ( kind = 4 ) max_edge
  integer ( kind = 4 ) max_order
  integer ( kind = 4 ) num_face

  integer ( kind = 4 ) edge(4,max_edge)
  integer ( kind = 4 ) face(max_order,num_face)
  integer ( kind = 4 ) face_object(num_face)
  integer ( kind = 4 ) face_order(num_face)
  integer ( kind = 4 ) face_rank(num_face)
  integer ( kind = 4 ) face_tier(num_face)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  integer ( kind = 4 ) num_edge
  integer ( kind = 4 ) num_fix
  integer ( kind = 4 ) num_object
!
!  Organize the faces into layered objects.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Determine edge-connected objects.'

  call object_build ( face, face_object, face_order, face_rank, face_tier, &
    max_order, num_face, num_object )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) 'Number of objects = ', num_object
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Face, Object, Tier'
  write ( *, '(a)' ) ' '

  do i = 1, num_face
    write ( *, '(3i8)' ) i, face_object(i), face_tier(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Preferred order:'
  write ( *, '(a)' ) '  Order, Face'
  write ( *, '(a)' ) ' '
  do i = 1, num_face
    write ( *, '(2i8)' ) i, face_rank(i)
  end do
!
!  Reorder the faces by object and tier.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Reorder the faces.'

  call face_sort ( face, face_object, face_order, face_tier, max_order, &
    num_face )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Face, Label, Object, Tier'
  write ( *, '(a)' ) ' '

  do i = 1, num_face
    write ( *, '(4i8)' ) i, face_rank(i), face_object(i), face_tier(i)
  end do
!
!  Construct the edge list.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Construct the edge list.'
  write ( *, '(a)' ) '(While doing so, check for edges used more'
  write ( *, '(a)' ) 'than twice.)'

  call face_to_edge ( edge, face, face_order, ierror, max_edge, max_order, &
    num_edge, num_face )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FACE_CHECK - Fatal error!'
    write ( *, '(a)' ) '  FACE_TO_EDGE failed.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Edge, Node1, Node2, Face1, Face2, Tier, Object'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' I, node1(i), node2(i), face1(i), face2(i)'
  write ( *, '(a)' ) ' '

  do i = 1, num_edge
    write ( *, '(10i3)' ) i, ( edge(j,i), j = 1, 4 )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Face, Order, Nodes'
  write ( *, '(a)' ) ' '
  do i = 1, num_face
    write ( *, '(10i3)' ) i, face_order(i), ( face(j,i), j = 1, face_order(i) )
  end do
!
!  Now force faces to have a consistent orientation.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Force faces to consistent orientation.'
  
  call face_flip ( edge, face, face_order, max_edge, max_order, num_edge, &
    num_face, num_fix )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Face, Order, Nodes'
  write ( *, '(a)' ) ' '
  do i = 1, num_face
    write ( *, '(10i3)' ) i, face_order(i), ( face(j,i), j = 1, face_order(i) )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'List boundary edges.'

  call edge_bound ( edge, max_edge, num_edge )

  return
end
subroutine face_example_box ( face, face_order, max_face, max_order, num_face )

!*****************************************************************************80
!
!! FACE_EXAMPLE_BOX returns the faces of a simple box.
!
!  Diagram:
!
!    1---------------------------4
!    |\                         /|
!    | \                       / |
!    |  \         1           /  |
!    |   \                   /   |
!    |    2-----------------3    |
!    |    |                 |    |
!    |    |                 |    |
!    |  3 |       4         | 5  |
!    |    |                 |    |
!    |    |                 |    |
!    |    6-----------------7    |
!    |   /                   \   |
!    |  /                     \  |
!    | /          2            \ |
!    |/                         \|
!    5---------------------------8
!
!  Discussion:
!
!    This routine is used to supply some very simple data for the 
!    face checking routines.
!
!    This is "almost" a cube, except that one face is missing.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) FACE(MAX_ORDER,NUM_FACE), describes the faces.
!    FACE(I,J) is the index of the I-th node in face J.  It is best
!    if the nodes of all faces are listed in counterclockwise order.
!
!    Output, integer ( kind = 4 ) FACE_ORDER(NUM_FACE), is the number of nodes
!    making up each face.
!
!    Input, integer ( kind = 4 ) MAX_FACE, the maximum number of faces allowed.
!
!    Input, integer ( kind = 4 ) MAX_ORDER, is the maximum number of nodes that can
!    make up a face, required to dimension FACE.
!
!    Output, integer ( kind = 4 ) NUM_FACE, the number of faces.
!
  implicit none

  integer ( kind = 4 ) max_order
  integer ( kind = 4 ) max_face

  integer ( kind = 4 ) face(max_order,max_face)
  integer ( kind = 4 ) face_order(max_face)
  integer ( kind = 4 ) num_face

  num_face = 5

  if ( max_face < num_face ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FACE_EXAMPLE_OPEN_BOX - Fatal error!'
    write ( *, '(a,i8)' ) '  Increase MAX_FACE to ', num_face
    stop
  end if

  face(1,1) = 1
  face(2,1) = 2
  face(3,1) = 3
  face(4,1) = 4

  face(1,2) = 5
  face(2,2) = 6
  face(3,2) = 7
  face(4,2) = 8

  face(1,3) = 1
  face(2,3) = 2
  face(3,3) = 6
  face(4,3) = 5

  face(1,4) = 6
  face(2,4) = 7
  face(3,4) = 3
  face(4,4) = 2

  face(1,5) = 3
  face(2,5) = 4
  face(3,5) = 8
  face(4,5) = 7

  face_order(1:num_face) = 4

  return
end
subroutine face_example_pieces ( face, face_order, max_face, max_order, &
  num_face )

!*****************************************************************************80
!
!! FACE_EXAMPLE_PIECES returns the faces of a set of three objects.
!
!  Diagram:
!
!    1---------------------------4
!    |\                         /|
!    | \                       / |       9--------10
!    |  \        7            /  |       |         |
!    |   \                   /   |       |   1     |
!    |    2-----------------3    |       |         |
!    |    |                 |    |       |         |
!    |    |                 |    |       11-------12
!    |  3 |       4         | 5  |        \       /
!    |    |                 |    |         \  6  /
!    |    |                 |    |          \   /
!    |    6-----------------7    |           \ /
!    |   /                   \   |           13
!    |  /                     \  |           / \
!    | /          8            \ |          /   \
!    |/                         \|         /  2  \
!    5---------------------------8        /       \
!                                        14-------15
!
!  Discussion:
!
!    THREE_PIECE is used to supply some very simple data for the 
!    face checking routines.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) FACE(MAX_ORDER,MAX_FACE), describes the faces.
!    FACE(I,J) is the index of the I-th node in face J.  It is best
!    if the nodes of all faces are listed in counterclockwise order.
!
!    Output, integer ( kind = 4 ) FACE_ORDER(MAX_FACE), is the number of nodes
!    making up each face.
!
!    Input, integer ( kind = 4 ) MAX_FACE, the maximum number of faces allowed.
!
!    Input, integer ( kind = 4 ) MAX_ORDER, is the maximum number of nodes that can
!    make up a face, required to dimension FACE.
!
!    Output, integer ( kind = 4 ) NUM_FACE, the number of faces.
!
  implicit none

  integer ( kind = 4 ) max_order
  integer ( kind = 4 ) max_face

  integer ( kind = 4 ) face(max_order,max_face)
  integer ( kind = 4 ) face_order(max_face)
  integer ( kind = 4 ) num_face

  num_face = 8

  if ( max_face < num_face ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FACE_EXAMPLE_PIECES - Fatal error!'
    write ( *, '(a)' ) '  MAX_FACE < NUM_FACE!'
    write ( *, '(a,i8)' ) '  NUM_FACE = ', num_face
    write ( *, '(a,i8)' ) '  MAX_FACE = ', max_face
    stop
  end if

  face(1,1) = 9
  face(2,1) = 10
  face(3,1) = 12
  face(4,1) = 11

  face(1,2) = 14
  face(2,2) = 13
  face(3,2) = 15

  face(1,3) = 1
  face(2,3) = 2
  face(3,3) = 6
  face(4,3) = 5

  face(1,4) = 6
  face(2,4) = 7
  face(3,4) = 3
  face(4,4) = 2

  face(1,5) = 3
  face(2,5) = 4
  face(3,5) = 8
  face(4,5) = 7

  face(1,6) = 13
  face(2,6) = 12
  face(3,6) = 11

  face(1,7) = 1
  face(2,7) = 2
  face(3,7) = 3
  face(4,7) = 4

  face(1,8) = 5
  face(2,8) = 6
  face(3,8) = 7
  face(4,8) = 8

  face_order(1) = 4
  face_order(2) = 3
  face_order(3) = 4
  face_order(4) = 4
  face_order(5) = 4
  face_order(6) = 3
  face_order(7) = 4
  face_order(8) = 4

  return
end
subroutine face_flip ( edge, face, face_order, max_edge, max_order, &
  num_edge, num_face, num_fix )

!*****************************************************************************80
!
!! FACE_FLIP flips faces to achieve a consistent orientation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) EDGE(4,MAX_EDGE), edge information.
!    EDGE(1,I) is the starting node of edge I;
!    EDGE(2,I) is the ending node of edge I;
!    EDGE(3,I) is the positive face;
!    EDGE(4,I) is the negative face, if any.
!
!    Input, integer ( kind = 4 ) FACE(MAX_ORDER,NUM_FACE), describes the faces.
!    FACE(I,J) is the index of the I-th node in face J.  It is best
!    if the nodes of all faces are listed in counterclockwise order.
!
!    Input, integer ( kind = 4 ) FACE_ORDER(NUM_FACE), is the number of nodes
!    making up each face.
!
!    Input, integer ( kind = 4 ) MAX_EDGE, the maximum number of edges.
!
!    Input, integer ( kind = 4 ) MAX_ORDER, the maximum number of nodes that can
!    make up a face, required to dimension FACE.
!
!    Input, integer ( kind = 4 ) NUM_EDGE, the number of edges.
!
!    Input, integer ( kind = 4 ) NUM_FACE, the number of faces.
!
!    Output, integer ( kind = 4 ) NUM_FIX, the number of bad faces that were found.
!
  implicit none

  integer ( kind = 4 ) max_edge
  integer ( kind = 4 ) max_order
  integer ( kind = 4 ) num_face

  integer ( kind = 4 ) edge(4,max_edge)
  integer ( kind = 4 ) f1
  integer ( kind = 4 ) f2
  integer ( kind = 4 ) face(max_order,num_face)
  integer ( kind = 4 ) face_order(num_face)
  integer ( kind = 4 ) iedge
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jp1
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) num_edge
  integer ( kind = 4 ) num_fix

  num_fix = 0

  do iedge = 1, num_edge

    n1 = edge(1,iedge)
    n2 = edge(2,iedge)
    f1 = edge(3,iedge)
    f2 = edge(4,iedge)
!
!  For now, just whine unless (N1,N2) is positive in F1 and negative in F2.
!
    if ( f1 /= 0 ) then

      do j = 1, face_order(f1)

        if ( j < face_order(f1) ) then
          jp1 = j + 1
        else
          jp1 = j
        end if

        m1 = face(j,f1)
        m2 = face(jp1,f1)

        if ( m1 == n1 .and. m2 == n2 ) then
          exit
        end if

        if ( m1 == n2 .and. m2 == n1 ) then
          num_fix = num_fix + 1
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'FACE_FLIP - Warning!'
          write ( *, '(a)' ) 'Bad orientation'
          write ( *, '(a,i8)' ) '  Face = ', f1
          write ( *, '(a,i8)' ) '  Side = ', j
          exit
        end if

      end do

    end if

    if ( f2 /= 0 ) then

      do j = 1, face_order(f2)

        if ( j < face_order(f2) ) then
          jp1 = j + 1
        else
          jp1 = j
        end if

        m1 = face(j,f2)
        m2 = face(jp1,f2)

        if ( m1 == n2 .and. m2 == n1 ) then
          exit
        end if

        if ( m1 == n1 .and. m2 == n2 ) then
          num_fix = num_fix + 1
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'FACE_FLIP - Warning!'
          write ( *, '(a)' ) 'Bad orientation'
          write ( *, '(a,i8)' ) '  Face = ', f2
          write ( *, '(a,i8)' ) '  Side = ', j
          exit
        end if

      end do

    end if

  end do

  if ( 0 < num_fix ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FACE_FLIP - Warning:'
    write ( *, '(a,i8)' ) '  Number of badly oriented faces = ', num_fix
  end if

  return
end
subroutine face_to_iv ( file_name, face, face_order, inode, jnode, nedge, &
  maxnode, maxface, maxorder, nnode, nface, x, y, z )

!*****************************************************************************80
!
!! FACE_TO_IV writes some simple graphics data to an Inventor file.
!
!  Example:
!
!     #Inventor V2.0 ascii
!
!     Separator {
!       Separator {
!         LightModel {
!           model PHONG
!         }
!         Material {
!           ambientColor  0.2 0.2 0.2
!           diffuseColor  0.8 0.8 0.8
!           emissiveColor 0.0 0.0 0.0
!           specularColor 0.0 0.0 0.0
!           shininess     0.2
!           transparency  0.0
!         }
!         Coordinate3 {
!           point [
!                8.59816       5.55317      -3.05561,
!                8.59816       2.49756      0.000000D+00,
!                ...etc...
!                2.48695       2.49756      -3.05561,
!           ]
!         }
!         IndexedFaceSet {
!           coordIndex [
!              0,    1,    2,   -1,    3,    4,    5,   -1,    7,    8,    9,
!            ...etc...
!            191,   -1,
!           ]
!         }
!       }
!     }
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the file name.
!
!    Input, integer ( kind = 4 ) FACE(MAX_ORDER,MAX_FACE), the nodes making faces.
!
!    Input, integer ( kind = 4 ) FACE_ORDER(MAX_FACE), the number of nodes per face.
!
!    Input, integer ( kind = 4 ) INODE(NEDGE), JNODE(NEDGE), node pairs for edges.
!
!    Input, integer ( kind = 4 ) NEDGE, the number of edges.
!
!    Input, integer ( kind = 4 ) MAXNODE, the maximum number of nodes.
!
!    Input, integer ( kind = 4 ) MAXFACE, the maximum number of faces.
!
!    Input, integer ( kind = 4 ) MAXORDER, the maximum number of nodes per face.
!
!    Input, integer ( kind = 4 ) NNODE, the number of points.
!
!    Input, integer ( kind = 4 ) NFACE, the number of faces.
!
!    Input, real ( kind = 8 ) X(MAXNODE), Y(MAXNODE), Z(MAXNODE), 
!    the coordinates of points.
!
  implicit none

  integer ( kind = 4 ), parameter :: OFFSET = 1

  integer ( kind = 4 ) maxnode
  integer ( kind = 4 ) maxface
  integer ( kind = 4 ) maxorder
  integer ( kind = 4 ) nedge

  integer ( kind = 4 ) face(maxorder,maxface)
  integer ( kind = 4 ) face_order(maxface)
  character ( len = * ) file_name
  integer ( kind = 4 ) icor3
  integer ( kind = 4 ) iface
  integer ( kind = 4 ) inode(nedge)
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) itemp
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) ivert
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jnode(nedge)
  integer ( kind = 4 ) length
  integer ( kind = 4 ) nnode
  integer ( kind = 4 ) nface
  character ( len = 200 ) text
  character ( len = 20 ) word
  real ( kind = 8 ) x(maxnode)
  real ( kind = 8 ) y(maxnode)
  real ( kind = 8 ) z(maxnode)

  call get_unit ( iunit )

  open ( unit = iunit, file = file_name, status = 'replace', iostat = ios )

  if ( ios /= 0 ) then
    return
  end if

  write ( iunit, '(a)' ) '#Inventor V2.0 ascii'
  write ( iunit, '(a)' ) ' '
  write ( iunit, '(a)' ) 'Separator {'
  write ( iunit, '(a)' ) '  Separator {'
!
!  LightModel:
!
!    BASE_COLOR ignores light sources, and uses only diffuse color
!      and transparency.  Even without normal vector information,
!      the object will show up.  However, you won't get shadow
!      and lighting effects.
!
!    PHONG uses the Phong lighting model, accounting for light sources
!      and surface orientation.  This is the default.  I believe
!      you need accurate normal vector information in order for this
!      option to produce nice pictures.
!
!    DEPTH ignores light sources, and calculates lighting based on
!      the location of the object within the near and far planes
!      of the current camera's view volume.
!
  write ( iunit, '(a)' ) '    LightModel {'
  write ( iunit, '(a)' ) '      model PHONG'
  write ( iunit, '(a)' ) '    }'
!
!  Material
!
  write ( iunit, '(a)' ) '    Material {'
  write ( iunit, '(a)' ) '      ambientColor  0.5 0.2 0.2'
  write ( iunit, '(a)' ) '      diffuseColor  0.5 0.2 0.3'
  write ( iunit, '(a)' ) '      emissiveColor 0.5 0.0 0.0'
  write ( iunit, '(a)' ) '      specularColor 0.5 0.0 0.0'
  write ( iunit, '(a)' ) '      shininess     0.5'
  write ( iunit, '(a)' ) '      transparency  0.0'
  write ( iunit, '(a)' ) '    }'
!
!  Point coordinates.
!
  write ( iunit, '(a)' ) '    Coordinate3 {'
  write ( iunit, '(a)' ) '      point ['

  do icor3 = 1, nnode
    write ( text, '(3f12.4,'','')' ) x(icor3), y(icor3), z(icor3)
    call s_blanks_delete ( text )
    write ( iunit, '(8x,a)' ) trim ( text )
  end do

  write ( iunit, '(a)' ) '      ]'
  write ( iunit, '(a)' ) '    }'
  write ( iunit, '(a)' ) '    IndexedLineSet {'
!
!  IndexedLineSet coordIndex
!
    write ( iunit, '(a)' ) '      coordIndex ['

    do j = 1, nedge
      write ( iunit, '(8x,i8,'','',i8,'','',i8,'','')' ) &
        inode(j) - OFFSET, jnode(j)-offset, -1
    end do

    write ( iunit, '(a)' ) '      ]'

    write ( iunit, '(a)' ) '    }'
!
!  IndexedFaceSet.
!
  if ( 0 < nface ) then

    write ( iunit, '(a)' ) '    IndexedFaceSet {'
!
!  IndexedFaceSet coordIndex
!
    write ( iunit, '(a)' ) '      coordIndex ['

    text = ' '
    length = 0

    do iface = 1, nface

      do ivert = 1, face_order(iface) + 1

        if ( ivert <= face_order(iface) ) then
          itemp = face(ivert,iface) - OFFSET
        else
          itemp = 0 - OFFSET
        end if

        write ( word, '(i8,'','')' ) itemp

        call s_cat ( text, word, text )
        length = length + 1

        if ( itemp == -1 .or. 10 <= length .or. &
          ( iface == nface .and. ivert == face_order(iface) + 1 )  ) then

          call s_blanks_delete ( text )
          write ( iunit, '(8x,a)' ) trim ( text )
          text = ' '
          length = 0

        end if

      end do

    end do

    write ( iunit, '(a)' ) '      ]'

    write ( iunit, '(a)' ) '    }'

  end if
!
!  Close up the Separator node.
!
  write ( iunit, '(a)' ) '  }'
!
!  Close up the Separator node.
!
  write ( iunit, '(a)' ) '}'

  close ( unit = iunit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FACE_TO_IV:'
  write ( *, '(a)' ) '  The data was written to the file: ' &
    // trim ( file_name )

  return
end
subroutine face_sort ( face, face_object, face_order, face_tier, max_order, &
  num_face )

!*****************************************************************************80
!
!! FACE_SORT renumbers the faces in order of object and tier.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) FACE(MAX_ORDER,NUM_FACE), describes the faces.
!    FACE(I,J) is the index of the I-th node in face J.  It is best
!    if the nodes of all faces are listed in counterclockwise order.
!
!    Input/output, integer ( kind = 4 ) FACE_OBJECT(NUM_FACE), describes the objects.
!    FACE_OBJECT(I) is the index of the edge-connected "object" to 
!    which face I belongs.
!
!    Input/output, integer ( kind = 4 ) FACE_ORDER(NUM_FACE), is the number of nodes
!    making up each face.
!
!    Input/output, integer ( kind = 4 ) FACE_TIER(NUM_FACE).  FACE_TIER(I) is the "tier"
!    of face I in its object.  The seed of the object has tier 1,
!    the neighbors of the seed have tier 2, and so on.
!
!    Input, integer ( kind = 4 ) MAX_ORDER, is the maximum number of nodes that can
!    make up a face, required to dimension FACE.
!
!    Input, integer ( kind = 4 ) NUM_FACE, the number of faces.
!
  implicit none

  integer ( kind = 4 ) max_order
  integer ( kind = 4 ) num_face

  integer ( kind = 4 ) face(max_order,num_face)
  integer ( kind = 4 ) face_object(num_face)
  integer ( kind = 4 ) face_order(num_face)
  integer ( kind = 4 ) face_tier(num_face)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iface
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) jface

  iface = 0
  jface = 0
  indx = 0
  isgn = 0

  do

    call sort_heap_external ( num_face, indx, iface, jface, isgn )
!
!  Interchange faces IFACE and JFACE.
!
    if ( 0 < indx ) then

      do i = 1, max_order
        call i4_swap ( face(i,iface), face(i,jface) )
      end do

      call i4_swap ( face_object(iface), face_object(jface) )
      call i4_swap ( face_order(iface), face_order(jface) )
      call i4_swap ( face_tier(iface), face_tier(jface) )
!
!  Compare faces IFACE and JFACE.
!
    else if ( indx < 0 ) then

      if ( ( face_object(iface) < face_object(jface) ) .or. &
           ( face_object(iface) == face_object(jface) .and. &
             face_tier(iface) < face_tier(jface) ) ) then
        isgn = -1
      else
        isgn = +1
      end if

    else

      exit

    end if

  end do

  return
end
subroutine face_to_edge ( edge, face, face_order, ierror, max_edge, &
  max_order, num_edge, num_face )

!*****************************************************************************80
!
!! FACE_TO_EDGE converts face data to edge data.
!
!  Discussion:
!
!    The computation will fail if:
!
!    * More than two faces claim to share an edge (Node1,Node2).
!    * Not enough storage is set aside by MAX_EDGE.
!
!    If is expected that the edge (Node1,Node2) in Face1 is traversed in
!    the opposite sense, as (Node2,Node1), in Face2.  If this is not the
!    case, then some faces may need to be reoriented, but that will not
!    affect the computation.
!    
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) EDGE(4,MAX_EDGE), edge information.
!    EDGE(1,I) is the starting node of edge I;
!    EDGE(2,I) is the ending node of edge I;
!    EDGE(3,I) is the positive face;
!    EDGE(4,I) is the negative face, or 0 if the edge is used once.
!
!    Input, integer ( kind = 4 ) FACE(MAX_ORDER,NUM_FACE), describes the faces.
!    FACE(I,J) is the index of the I-th node in face J.  It is best
!    if the nodes of all faces are listed in counterclockwise order.
!
!    Input, integer ( kind = 4 ) FACE_ORDER(NUM_FACE), is the number of nodes
!    making up each face.
!
!    Output, integer ( kind = 4 ) IERROR, error flag: 0 = no error, nonzero = error.
!
!    Input, integer ( kind = 4 ) MAX_EDGE, the maximum number of edges.
!
!    Input, integer ( kind = 4 ) MAX_ORDER, the maximum number of nodes that can
!    make up a face, required to dimension FACE.
!
!    Output, integer ( kind = 4 ) NUM_EDGE, the number of edges.
!
!    Input, integer ( kind = 4 ) NUM_FACE, the number of faces.
!
  implicit none

  integer ( kind = 4 ) max_edge
  integer ( kind = 4 ) max_order
  integer ( kind = 4 ) num_face

  integer ( kind = 4 ) edge(4,max_edge)
  integer ( kind = 4 ) face(max_order,num_face)
  integer ( kind = 4 ) face_order(num_face)
  integer ( kind = 4 ) iedge
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iface
  integer ( kind = 4 ) index
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jp1
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) num_edge
!
!  Initialize.
!
  ierror = 0

  edge(1:4,1:max_edge) = 0

  num_edge = 0
!
!  Consider face #I.
!
  do iface = 1, num_face
!
!  Seek an edge of face IFACE that already occurs in the edge list.
!  If there is one, then slide and reverse the entries in FACE(*,IFACE)
!  so that that edge occurs first, and in the opposite sense to its
!  occurrence in the edge list.
!
    call edge_match_face ( edge, max_edge, num_edge, face(1,iface), &
      face_order(iface), index )
!
!  Now, in any case, we know that the first two nodes in FACE(*,IFACE)
!  are the negative of an existing edge, or no nodes in FACE(*,IFACE)
!  occur in any existing edge.
!
    do j = 1, face_order(iface)

      n1 = face(j,iface)

      if ( j == face_order(iface) ) then
        jp1 = 1
      else
        jp1 = j + 1
      end if

      n2 = face(jp1,iface)

      call edge_match_nodes ( edge, max_edge, num_edge, n1, n2, iedge )

      if ( iedge == 0 ) then

        call edge_add_nodes ( edge, max_edge, num_edge, iface, n1, n2, ierror )

        if ( ierror /= 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'FACE_TO_EDGE - Fatal error!'
          write ( *, '(a)' ) '  EDGE_ADD_NODES failed.'
          ierror = 1
          return
        end if

      else if ( edge(4,iedge) == 0 ) then

        edge(4,iedge) = iface

      else 

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'FACE_TO_EDGE - Fatal error!'
        write ( *, '(a,2i8)' ) '  Edge between nodes ', &
          edge(1,iedge), edge(2,iedge)
        write ( *, '(a)' ) '  is used at least 3 times, by faces:'
        write ( *, '(3i8)' ) edge(3,iedge), edge(4,iedge), iface
        ierror = 1
        return

      end if

    end do
  end do

  return
end
subroutine face_touch ( face, face_order, max_order, num_face, iface, jface, &
  touch )

!*****************************************************************************80
!
!! FACE_TOUCH reports whether two polygonal faces touch.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) FACE(MAX_ORDER,NUM_FACE), describes the faces.
!    FACE(I,J) is the index of the I-th node in face J.  It is best
!    if the nodes of all faces are listed in counterclockwise order.
!
!    Input, integer ( kind = 4 ) FACE_ORDER(NUM_FACE), is the number of nodes
!    making up each face.
!
!    Input, integer ( kind = 4 ) MAX_ORDER, is the maximum number of nodes that can
!    make up a face, required to dimension FACE.
!
!    Input, integer ( kind = 4 ) NUM_FACE, the number of faces.
!
!    Input, integer ( kind = 4 ) IFACE, JFACE, the faces to be checked.
!
!    Output, integer ( kind = 4 ) TOUCH:
!     0, the faces do not touch;
!    +1, the faces touch, both using an arc in the same direction;
!    -1, the faces touch, using an arc in opposite directions.
!
  implicit none

  integer ( kind = 4 ) max_order
  integer ( kind = 4 ) num_face

  integer ( kind = 4 ) face(max_order,num_face)
  integer ( kind = 4 ) face_order(num_face)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iface
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jface
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mp1
  integer ( kind = 4 ) mm1
  integer ( kind = 4 ) n
  integer ( kind = 4 ) np1
  integer ( kind = 4 ) touch

  touch = 0
!
!  Arc N1-N2 on IFACE must be matched by arc N1-N2 or N2-N1 on JFACE.
!
  do i = 1, face_order(iface)

    n = face(i,iface)
    if ( i < face_order(iface) ) then
      np1 = face(i+1,iface)
    else
      np1 = face(1,iface)
    end if

    do j = 1, face_order(jface)

      m = face(j,jface)
      if ( j < face_order(jface) ) then
        mp1 = face(j+1,jface)
      else
        mp1 = face(1,jface)
      end if

      if ( 1 < j ) then
        mm1 = face(j-1,jface)
      else
        mm1 = face(face_order(jface),jface)
      end if

      if ( n == m ) then
        if ( np1 == mp1 ) then
          touch = + 1
          return
        else if ( np1 == mm1 ) then
          touch = - 1
          return
        end if
      end if

    end do
  end do

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

    if ( i /= 5 .and. i /= 6 ) then

      inquire ( unit = i, opened = lopen, iostat = ios )

      if ( ios == 0 ) then
        if ( .not. lopen ) then
          iunit = i
          return
        end if
      end if

    end if

  end do
!
!  No free unit was found.
!
  iunit = 0

  return
end
subroutine graph_adj_bfs ( adj, lda, nnode, dad, deep, order )

!*****************************************************************************80
!
!! GRAPH_ADJ_BFS carries out a breadth-first traversal of a graph.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 April 1999
!
!  Reference:
!
!    Alan Gibbons,
!    Algorithmic Graph Theory,
!    Cambridge University Press, 1985,
!    ISBN 0-521-28881-9.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(LDA,NNODE), is the adjacency information.
!    ADJ(I,J) is nonzero if there is an edge from node
!    I to node J, and 0 otherwise.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of ADJ, which must be
!    at least NNODE.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) DAD(NNODE), DAD(I) is the node from which
!    node I is visited.  Node 1 is the first node in the search,
!    and has no predecessor, so DAD(1) is zero.  If there is
!    more than one connected component, then there
!    will be other nodes with DAD equal to zero.
!
!    Output, integer ( kind = 4 ) DEEP(NNODE), records the "depth" of the node.
!    The first node, node 1, has depth 1.  All the nodes that
!    can be reached in one step from node 1 have depth 2.  All
!    nodes that can be reached in one step from any of those nodes
!    have depth 3.  If there is more than one connected component,
!    then the depth of nodes in the second component will begin
!    one greater than the greatest depth of the first component,
!    and so on.
!
!    Output, integer ( kind = 4 ) ORDER(NNODE).  ORDER(I) is the step at which
!    node I is visited in the search.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) dad(nnode)
  integer ( kind = 4 ) deep(nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) order(nnode)
  integer ( kind = 4 ) iput
  integer ( kind = 4 ) queue(nnode)
  integer ( kind = 4 ) itake
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jdeep
  integer ( kind = 4 ) k
  integer ( kind = 4 ) nudeep

  deep(1:nnode) = 0
  order(1:nnode) = 0
  dad(1:nnode) = 0
 
  k = 0
  i = 1
  iput = 1
  itake = 1
  nudeep = iput
  queue(iput) = i
  jdeep = 1
  deep(i) = jdeep
  k = k + 1
  order(i) = k
  dad(i) = 0
!
!  Find all sons of this father.
!  Store all sons in the son stack.
!
10    continue
 
  do j = 1, nnode
 
    if ( ( adj(i,j) /= 0 .or. adj(j,i) /= 0 ) .and. order(j) == 0 ) then

      iput = iput + 1

      if ( nnode < iput ) then
        iput = 1
      end if

      queue(iput) = j
      k = k + 1
      dad(j) = i
      order(j) = k
      deep(j) = jdeep + 1

    end if
 
  end do
!
!  Are there more fathers whose sons are to be searched for?
!
  if ( iput /= itake ) then
 
    if ( itake == nudeep ) then
      jdeep = jdeep + 1
      nudeep = iput
    end if
 
    i = queue(itake)
    itake = itake + 1

    if ( nnode < itake ) then
      itake = 1
    end if

    go to 10
!
!  No more fathers, no more sons.  Is there an unvisited component?
!
  else
 
    do i = 1, nnode
 
      if ( order(i) == 0 ) then
        itake = 1
        iput = 1
        queue(iput) = i
        jdeep = jdeep + 1
        nudeep = 1
        k = k + 1
        order(i) = k
        deep(i) = jdeep
        dad(i) = 0
        go to 10
      end if
 
    end do
 
  end if
 
  return
end
subroutine graph_adj_bipartite_random ( lda, nnode1, nnode2, seed, nedge, adj )

!*****************************************************************************80
!
!! GRAPH_ADJ_BIPARTITE_RANDOM generates a random bipartite graph.
!
!  Definition:
!
!    A bipartite graph has the property that its nodes may be divided
!    into two groups, NODE1 and NODE2, with the property that the only
!    edges in the graph are between a node in NODE1 and a node in NODE2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of LDA, which must be
!    at least NNODE1+NNODE2.
!
!    Input, integer ( kind = 4 ) NNODE1, NNODE2, the number of nodes in the first and
!    second groups of nodes.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
!    Output, integer ( kind = 4 ) NEDGE, the number of edges.
!
!    Output, integer ( kind = 4 ) ADJ(LDA,NNODE1+NNODE2), the adjacency matrix.  ADJ(I,J) is
!    nonzero if there is an edge from node I to node J.  ADJ(I,I) will
!    always be 0.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) nnode1
  integer ( kind = 4 ) nnode2
  integer ( kind = 4 ) nedge

  integer ( kind = 4 ) adj(lda,nnode1+nnode2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) nnode
  integer ( kind = 4 ) seed

  if ( nnode1 <= 0  ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRAPH_ADJ_BIPARTITE_RANDOM - Fatal error!'
    write ( *, '(a,i8)' ) '  NNODE1 = ', nnode1
    write ( *, '(a)' ) '  but NNODE1 must be at least 1.'
    stop
  end if

  if ( nnode2 <= 0  ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRAPH_ADJ_BIPARTITE_RANDOM - Fatal error!'
    write ( *, '(a,i8)' ) '  NNODE2 = ', nnode2
    write ( *, '(a)' ) '  but NNODE2 must be at least 1.'
    stop
  end if

  nnode = nnode1 + nnode2
  nedge = 0

  adj(1:nnode,1:nnode) = 0
!
!  For each node in the NODE1 group, 
!  consider a edge to each node in the NODE2 group.
!
  do i = 1, nnode1
    do j = nnode1+1, nnode1+nnode2

      k = i4_uniform ( 0, 1, seed )

      adj(i,j) = k
      adj(j,i) = k
      nedge = nedge + k

    end do
  end do
!
!  Now perform a random permutation of the rows and columns.
!
  call i4mat_perm_random ( lda, nnode, seed, adj )

  return
end
subroutine graph_adj_block ( adj, lda, nnode, dad, order, stack, nblock ) 

!*****************************************************************************80
!
!! GRAPH_ADJ_BLOCK finds the blocks of an undirected graph from its adjacency list.
!
!  Definition:
!
!    A component of a graph is a connected subset of the graph.  If a node
!    is in the component, then all nodes to which it is connected are also
!    in the component.
!
!    An articulation point of a component of a graph is a node whose
!    removal causes the component to no longer be connected.
!
!    A component with no articulation points is called a block.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 April 1999
!
!  Reference:
!
!    Alan Gibbons,
!    Algorithmic Graph Theory,
!    Cambridge University Press, 1985,
!    ISBN 0-521-28881-9.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) ADJ(LDA,NNODE).
!    On input, ADJ is the adjacency matrix.  ADJ(I,J) is
!    positive if there is an edge from node I to node J, and 0 otherwise.
!    On output, each positive entry of ADJ has been replaced
!    by the number of the block that the corresponding edge belongs to.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of ADJ, which must be
!    at least NNODE.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) DAD(NNODE), DAD(I) is the node from which
!    node I is visited.  Node 1 is the first node in the search,
!    and has no predecessor, so DAD(1) is zero.  If there is
!    more than one connected component in the graph, then there
!    will be other nodes with DAD equal to zero.
!
!    Output, integer ( kind = 4 ) ORDER(NNODE).  ORDER(I) records the order
!    in which the node was visited during the depth-first search.
!    The first node, node 1, has ORDER(1) = 1.
!    Note, however, that any node which is an articulation point
!    will have the value of ORDER negated.
!
!    Workspace, integer STACK(NNODE).
!
!    Output, integer ( kind = 4 ) NBLOCK, the number of blocks.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) dad(nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) idir
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) inode(nnode)
  integer ( kind = 4 ) order(nnode)
  integer ( kind = 4 ) iroot
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jedge
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) jnode(nnode)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) label(nnode)
  integer ( kind = 4 ) lstack
  integer ( kind = 4 ) nblock
  integer ( kind = 4 ) stack(nnode)

  dad(1:nnode) = 0
  inode(1:nnode) = 0
  order(1:nnode) = 0
  stack(1:nnode) = 0
  jnode(1:nnode) = 0
  label(1:nnode) = 0
 
  nblock = 0
  k = 0
  i = 1
  lstack = 0
  jedge = 0
!
!  Find all descendants of the parent node in this connected component
!  of the graph.
!
10    continue
 
  iroot = i
  k = k + 1
  order(i) = k
  label(i) = k
  lstack = lstack + 1
  stack(lstack) = i
  idir = + 1
 
30    continue
 
  j = 0
!
!  Check the next neighbor.
!
40    continue
 
  j = j + 1

  if ( nnode < j ) then
    go to 50
  end if
 
  if ( adj(i,j) /= 0 .or. adj(j,i) /= 0 ) then
 
    if ( 0 < adj(i,j) .or. 0 < adj(j,i) ) then
      jedge = jedge + 1
      inode(jedge) = i
      jnode(jedge) = j
    end if
 
    if ( order(j) == 0 ) then
 
      dad(j) = i
      lstack = lstack + 1
      stack(lstack) = j
      idir = + 1
      k = k + 1
      i = j
      order(i) = k
      label(i) = k
      go to 30
 
    else
 
      if ( idir == +1 ) then
        label(i) = min ( label(i), abs ( order(j) ) )
      else
        label(i) = min ( label(i), label(j) )
      end if
 
    end if
 
  end if
 
  go to 40
!
!  Searched all directions from current node.  Back up one node,
!  or, if stack is exhausted, look for a node we haven't visited,
!  which therefore belongs to a new connected component.
!
50    continue
 
  lstack = lstack - 1
  idir = -1
 
  if ( 0 < lstack ) then
 
    j = i
    i = stack(lstack)
 
    if ( abs ( order(i) ) <= label(j) ) then
 
      if ( 0 < order(i) ) then
 
        if ( i /= iroot ) then
          order(i) = - order(i)
        else
          iroot = 0
        end if
 
      end if
 
      nblock = nblock + 1

      do

        ii = inode(jedge)
        jj = jnode(jedge)
        jedge = jedge - 1
        adj(ii,jj) = - nblock
        adj(jj,ii) = - nblock

        if ( ii == i .and. jj == j ) then
          exit
        end if

      end do
 
    end if
 
    go to 40
 
  else
 
    lstack = 0
 
    do l = 1, nnode
      if ( order(l) == 0 ) then
        i = l
        go to 10
      end if
    end do
 
  end if
!
!  Restore the positive sign of the adjacency matrix.
!
  adj(1:nnode,1:nnode) = abs ( adj(1:nnode,1:nnode) )
 
  return
end
subroutine graph_adj_closure ( adj, lda, nnode )

!*****************************************************************************80
!
!! GRAPH_ADJ_CLOSURE generates the transitive closure of a graph.
!
!  Discussion:
!
!    The method is due to S Warshall.
!
!    The transitive closure of a graph is a function REACH(I,J) so that
!
!      REACH(I,J) = 0 if node J cannot be reached from node I;
!                   1 if node J can be reached from node I.
!
!    This is an extension of the idea of adjacency.  REACH(I,J)=1 if
!    node J is adjacent to node I, or if node J is adjacent to a node
!    that is adjacent to node I, etc.
!
!    Note that if a graph is (node) connected, then its transitive closure
!    is the matrix that is 1 everywhere.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Robert Sedgewick,
!    Algorithms,
!    Addison Wesley, 1983, page 425.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) ADJ(LDA,NNODE).
!
!    On input, ADJ is the adjacency matrix for the graph.  ADJ(I,J)
!    is nonzero if there is an edge from node I to node J.
!
!    On output, ADJ is the transitive closure matrix of the graph.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of LDA, which must be
!    at least NNODE.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k

  do i = 1, nnode
    adj(i,i) = 1
  end do

  do i = 1, nnode
    do j = 1, nnode
      if ( adj(j,i) /= 0 .or. adj(i,j) /= 0 ) then
        do k = 1, nnode
          if ( adj(i,k) /= 0 .or. adj(k,i) /= 0 ) then
            adj(j,k) = 1
            adj(k,j) = 1
          end if
        end do
      end if
    end do
  end do

  return
end
subroutine graph_adj_color_cand ( adj, lda, nnode, ncolor, color, k, nstack, &
  stack, maxstack, ncan )

!*****************************************************************************80
!
!! GRAPH_ADJ_COLOR_CAND finds possible colors for a node during a graph coloring.
!
!  Discussion:
!
!    This routine is given a partial coloring of the graph.  
!    The total coloring of the graph must be done in such a way that no 
!    two nodes joined by an edge have the same color.
!
!    This routine must be used in conjunction with I4VEC_BACKTRACK.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 August 2000
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
!    Input, integer ( kind = 4 ) ADJ(LDA,NNODE), the adjacency matrix.  ADJ(I,J) is 
!    nonzero if there is an edge from node I to node J.
!
!    Input, integer ( kind = 4 ) LDA, the first dimension of ADJ as declared in
!    the calling program.  LDA must be at least NNODE.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) NCOLOR, the number of colors available.
!
!    Input, integer ( kind = 4 ) COLOR(NNODE).  COLOR(I) is the color of node I.
!
!    Input, integer ( kind = 4 ) K, node whose possible colors are to be found.
!
!    Input/output, integer ( kind = 4 ) NSTACK, current length of stack.
!
!    Workspace, integer STACK(MAXSTACK), candidates for the colors of the nodes.
!
!    Input, integer ( kind = 4 ) MAXSTACK, dimension of STACK.
!
!    input, integer ( kind = 4 ) NCAN(NNODE), the number of candidates for each position.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) nnode
  integer ( kind = 4 ) maxstack

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) color(nnode)
  integer ( kind = 4 ) i
  logical lwork(nnode)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) ncan(nnode)
  integer ( kind = 4 ) nstack
  integer ( kind = 4 ) ncolor
  integer ( kind = 4 ) stack(maxstack)

  ncan(k) = 0

  if ( k <= 1 ) then

    stack(1) = 1
    nstack = 1
    ncan(k) = 1

  else
 
    lwork(1:ncolor) = .true.
 
    do i = 1, k-1
      if ( adj(i,k) /= 0 .or. adj(k,i) /= 0 ) then
        lwork(color(i)) = .false.
      end if
    end do
  
    do i = 1, ncolor
 
      if ( lwork(i) ) then
        nstack = nstack + 1
        stack(nstack) = i
        ncan(k) = ncan(k) + 1
      end if
 
    end do
 
  end if
 
  return
end
subroutine graph_adj_color_next ( adj, lda, nnode, ncolor, color, stack, &
  maxstack, ncan, more )

!*****************************************************************************80
!
!! GRAPH_ADJ_COLOR_NEXT returns possible colorings of a graph, one at a time.
!
!  Definition:
!
!    A coloring of a graph using NCOLOR colors is an assignment to each
!    node of a label between 1 and NCOLOR, in such a way that no two
!    neighboring nodes have the same label.
!
!  Method:
!
!    This routine uses the backtracking method to produce the colorings.
!    Routine GRAPH_ADJ_COLOR_CAND produces candidates for a partial solution, 
!    and routine I4VEC_BACKTRACK assembles the total solution.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 July 1998
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
!    Input, integer ( kind = 4 ) ADJ(LDA,NNODE), the adjacency matrix.  ADJ(I,J) is nonzero 
!    if there is an edge between node I to node J.
!
!    Input, integer ( kind = 4 ) LDA, the first dimension of ADJ as declared in
!    the calling program.  LDA must be at least NNODE.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) NCOLOR, the number of colors available.
!
!    Output, integer ( kind = 4 ) COLOR(NNODE).  On return with MORE = TRUE, COLOR(I)
!    is the color of node I.
!
!    Workspace, integer STACK(MAXSTACK), candidates for the colors of nodes
!    1 through K-1.
!
!    Input, integer ( kind = 4 ) MAXSTACK, dimension of STACK.
!
!    Workspace, integer NCAN(NNODE), the number of candidates for each position.
!
!    Input/output, logical MORE.
!    On first call, set MORE to .FALSE, and do not alter it after.
!    On return, MORE is TRUE if another coloring has been returned in
!    IARRAY, and FALSE if there are no more colorings.
!
  implicit none

  integer ( kind = 4 ) nnode
  integer ( kind = 4 ) lda
  integer ( kind = 4 ) maxstack

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) color(nnode)
  integer ( kind = 4 ), save :: indx = 0
  integer ( kind = 4 ), save :: k = 0
  logical more
  integer ( kind = 4 ) ncan(nnode)
  integer ( kind = 4 ), save :: nstack = 0
  integer ( kind = 4 ) ncolor
  integer ( kind = 4 ) stack(maxstack)
!
!  First call.
!
  if ( .not. more ) then
    indx = 0
    k = 0
    more = .true.
    nstack = 0
  end if
 
  do
 
    call i4vec_backtrack ( nnode, color, indx, k, nstack, stack, maxstack, &
      ncan )
 
    if ( indx == 1 ) then

      exit

    else if ( indx == 2 ) then

      call graph_adj_color_cand ( adj, lda, nnode, ncolor, color, k, nstack, &
        stack, maxstack, ncan )

    else

      more = .false.
      exit

    end if

  end do
 
  return
end
subroutine graph_adj_complement ( adj, lda, nnode )

!*****************************************************************************80
!
!! GRAPH_ADJ_COMPLEMENT returns the adjacency matrix of the complement of a graph.
!
!  Definition:
!
!    This routine can also handle a directed graph.
!
!    The complement of a graph G is a graph H with the property that
!    nodes u and v are connected in H if and only if they are not
!    connected in G.  However, edges from a node to itself are not
!    allowed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 August 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) ADJ(LDA,NNODE).  On input, the
!    adjacency information for the graph G.  On output, ADJ 
!    the adjacency information for the complement of G.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of ADJ, which must be
!    NNODE or greater.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
!
!  Force the adjacency graph to be symmetric.
!
  call graph_adj_symmetrize ( adj, lda, nnode )

  do i = 1, nnode
    do j = 1, nnode

      if ( i == j ) then
        adj(i,j) = 0
      else if ( adj(i,j) /= 0 ) then
        adj(i,j) = 0
      else if ( adj(i,j) == 0 ) then
        adj(i,j) = 1
      end if

    end do
  end do
 
  return
end
subroutine graph_adj_connect_random ( lda, nnode, nedge, seed, adj )

!*****************************************************************************80
!
!! GRAPH_ADJ_CONNECT_RANDOM generates a random connected graph.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) ADJ(LDA,NNODE), the adjacency matrix.  ADJ(I,J) is
!    nonzero if there is an edge from node I to node J.  ADJ(I,I) will
!    always be 0.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of LDA, which must be
!    at least NNODE.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) NEDGE, the number of edges, which must be between
!    NNODE-1 and (NNODE*(NNODE-1))/2.  
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) nnode
  integer ( kind = 4 ) nedge

  integer ( kind = 4 ) adj(lda,nnode)
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
  integer ( kind = 4 ) nnode2
  integer ( kind = 4 ) seed
!
!  Check.
!
  if ( nnode <= 0  ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRAPH_ADJ_CONNECT_RANDOM - Fatal error!'
    write ( *, '(a,i8)' ) '  NNODE = ', nnode
    write ( *, '(a)' ) '  but NNODE must be at least 1.'
    stop
  end if

  if ( lda < nnode ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRAPH_ADJ_CONNECT_RANDOM - Fatal error!'
    write ( *, '(a,i8)' ) '  LDA = ', lda
    write ( *, '(a,i8)' ) '  but LDA must be at least NNODE = ', nnode
    stop
  end if

  maxedge = ( nnode * ( nnode - 1 ) ) / 2

  if ( nedge < nnode-1 .or. maxedge < nedge ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRAPH_ADJ_CONNECT_RANDOM - Fatal error!'
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
  call tree_arc_random ( nnode, seed, code, inode, jnode )
!
!  Convert information to adjacency form.
!
  call graph_arc_to_graph_adj ( nnode-1, inode, jnode, adj, lda, nnode2 )
!
!  Now we have NEDGE - ( NNODE - 1 ) more edges to add.
!
  nchoice = ( nnode * ( nnode - 1 ) ) / 2 - ( nnode - 1 )
  nchoose = nedge - ( nnode - 1 )

  call ksub_random ( nchoice, nchoose, seed, iwork )

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
subroutine graph_adj_cycle ( adj, lda, nnode, dad, order, maxstack, stack )

!*****************************************************************************80
!
!! GRAPH_ADJ_CYCLE searches for cycles in a graph.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 February 1999
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) ADJ(LDA,NNODE), the adjacency matrix for 
!    the graph.  ADJ(I,J) is 0 if there is no edge from node I to node J.
!
!    On input, ADJ(I,J) should be 1 if there is an edge from node I to node J.
!
!    On output, ADJ(I,J) will be one of the following values:
!      -1 if the edge from node I to node J is part of at least one cycle;
!      -2 if the edge from node I to node J is part of the depth first
!      search trees.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of ADJ, which must
!    be at least NNODE.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) DAD(NNODE), the father array for the depth first
!    search trees.  DAD(I) = 0 means that node I is the root of 
!    one of the trees.  DAD(I) = J means that the search descended
!    from node J to node I.
!
!    Output, integer ( kind = 4 ) ORDER(NNODE), the order in which the nodes were
!    traversed, from 1 to NNODE.
!
!    Input, integer ( kind = 4 ) MAXSTACK, the amount of stack space available.
!    The absolute maximum needed would be 2*(NNODE-1) though less
!    is likely.
!
!    Workspace, integer STACK(MAXSTACK).
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) maxstack
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) dad(nnode)
  integer ( kind = 4 ) daddy
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nstack
  integer ( kind = 4 ) order(nnode)
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) stack(maxstack)

  dad(1:nnode) = 0
  order(1:nnode) = 0

  rank = 0

  do i = 1, nnode

    if ( order(i) == 0 ) then

      daddy = i
      nstack = 0
!
!  Visit node DAD.
!
10    continue

      rank = rank + 1
      order(daddy) = rank
      j = 0
!
!  Consider visiting node J from node DAD.
!
20    continue

      j = j + 1
!
!  If J is a reasonable value, adjacent to DAD, and unvisited,
!  then put DAD into the stack, make J the new value of DAD,
!  and go to 10.
!
      if ( j <= nnode ) then

        if ( 0 < adj(daddy,j) .or. 0 < adj(j,daddy) ) then

          if ( order(j) == 0 ) then

            adj(daddy,j) = - 2
            adj(j,daddy) = - 2

            if ( nstack+2 <= maxstack ) then
              dad(j) = daddy
              stack(nstack+1) = daddy
              stack(nstack+2) = j
              nstack = nstack + 2
              daddy = j
              go to 10
            else
              write ( *, '(a)' ) ' '
              write ( *, '(a)' ) 'GRAPH_CYCLE - Fatal error!'
              write ( *, '(a)' ) '  Out of stack space.'
              stop
            end if
!
!  An adjacent node has already been visited.  This constitutes a cycle.
!
          else

            adj(j,daddy) = - 1
            adj(daddy,j) = - 1
            go to 20

          end if
!
!  If J is not suitable for a visit, get the next value of J.
!
        else

          go to 20

        end if
!
!  If no more neighbors to consider, back up one node.
!
      else if ( 2 <= nstack ) then

        daddy = stack(nstack-1)
        j = stack(nstack)
        nstack = nstack - 2
        go to 20
!
!  If no more nodes to consider in this tree, bail out.
!
      else

        nstack = 0

      end if

    end if

  end do

  return
end
subroutine graph_adj_degree ( adj, lda, nnode, degree )

!*****************************************************************************80
!
!! GRAPH_ADJ_DEGREE computes the degree of each node.
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
!    01 October 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(LDA,NNODE), the adjacency information.
!    ADJ(I,J) is 1 if there is an edge from node I to node J.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the ADJ array,
!    which must be at least NNODE.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) DEGREE(NNODE), the degree of the nodes.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(lda,nnode)
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
subroutine graph_adj_degree_max ( adj, lda, nnode, degree_max )

!*****************************************************************************80
!
!! GRAPH_ADJ_DEGREE_MAX computes the maximum node degree.
!
!  Discussion:
!
!    The maximum node degree of a graph is the maximum value of the
!    degree of the nodes of the graph.
!
!    If two graphs are isomorphic, they must have the same maximum node degree.
!
!    If two graphs have different maximum node degrees, they cannot
!    be isomorphic.
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
!    Input, integer ( kind = 4 ) ADJ(LDA,NNODE), the adjacency information.
!    ADJ(I,J) is 1 if there is an edge from node I to node J.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the ADJ array,
!    which must be at least NNODE.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) DEGREE_MAX, the maximum node degree of the graph.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(lda,nnode)
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
subroutine graph_adj_degree_seq ( adj, lda, nnode, seq )

!*****************************************************************************80
!
!! GRAPH_ADJ_DEGREE_SEQ computes the degree sequence of a graph.
!
!  Discussion:
!
!    The degree sequence of a graph is constructed by computing the
!    degree of each node, and then ordering these values in decreasing order.
!
!    If two graphs are isomorphic, they must have the same degree sequence.
!
!    If two graphs have different degree sequences, they cannot be isomorphic.
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
!    Input, integer ( kind = 4 ) ADJ(LDA,NNODE), the adjacency information.
!    ADJ(I,J) is 1 if there is an edge from node I to node J.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the ADJ array,
!    which must be at least NNODE.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) SEQ(NNODE), the degree sequence of the graph.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(lda,nnode)
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
subroutine graph_adj_dfs ( adj, lda, nnode, dad, order )

!*****************************************************************************80
!
!! GRAPH_ADJ_DFS does a depth first search of a graph.
!
!  Discussion:
!
!    The routine returns:
!
!    * a list of the order in which the nodes were visited,
!    * a list of the parents of each node in the search tree,
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 July 2000
!
!  Reference:
!
!    Robert Sedgewick,
!    Algorithms,
!    Addison Wesley, 1983, page 382.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(LDA,NNODE), the adjacency matrix for the graph.
!    ADJ(I,J) is 0 if node J is not adjacent to node I, and nonzero
!    otherwise.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of ADJ, which must be
!    at least NNODE.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) DAD(NNODE), the father array for the depth first
!    search trees.  DAD(I) = 0 means that node I is the root of 
!    one of the trees.  DAD(I) = J means that the search descended
!    from node J to node I.
!
!    Output, integer ( kind = 4 ) ORDER(NNODE), the order in which the nodes were
!    traversed, from 1 to NNODE.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) dad(nnode)
  integer ( kind = 4 ) daddy
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) maxstack
  integer ( kind = 4 ) nstack
  integer ( kind = 4 ) order(nnode)
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) stack(2*(nnode-1))

  dad(1:nnode) = 0
  maxstack = 2 * ( nnode - 1 )
  order(1:nnode) = 0

  rank = 0

  do i = 1, nnode

    if ( order(i) == 0 ) then

      daddy = i
      nstack = 0
!
!  Visit node DAD.
!
10    continue

      rank = rank + 1
      order(daddy) = rank
      j = 0
!
!  Consider visiting node J from node DAD.
!
20    continue

      j = j + 1
!
!  If J is a reasonable value, adjacent to DAD, and unvisited,
!  then put DAD into the stack, make J the new value of DAD,
!  and go to 10.
!
      if ( j <= nnode ) then

        if ( adj(daddy,j) /= 0 .and. order(j) == 0 ) then

          if ( nstack+2 <= maxstack ) then
            dad(j) = daddy
            stack(nstack+1) = daddy
            stack(nstack+2) = j
            nstack = nstack + 2
            daddy = j
            go to 10
          else
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'GRAPH_ADJ_DFS - Fatal error!'
            write ( *, '(a)' ) '  Out of stack space.'
            stop
          end if
!
!  If J is not suitable for a visit, get the next value of J.
!
        else

          go to 20

        end if
!
!  If no more neighbors to consider, back up one node.
!
      else if ( 2 <= nstack ) then

        daddy = stack(nstack-1)
        j = stack(nstack)
        nstack = nstack - 2
        go to 20
!
!  If no more nodes to consider in this tree, bail out.
!
      else

        nstack = 0

      end if

    end if

  end do

  return
end
subroutine graph_adj_dfs_2 ( adj, lda, nnode, dad, order )

!*****************************************************************************80
!
!! GRAPH_ADJ_DFS_2 does a depth-first search of a graph.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 April 1999
!
!  Reference:
!
!    Alan Gibbons,
!    Algorithmic Graph Theory,
!    Cambridge University Press, 1985,
!    ISBN 0-521-28881-9.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(LDA,NNODE), is the adjacency matrix of
!    the graph.  ADJ(I,J) is nonzero if there is an edge from node
!    I to node J.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of ADJ, which must
!    be at least NNODE.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) DAD(NNODE), DAD(I) is the node from which
!    node I is visited.  Node 1 is the first node in the search,
!    and has no predecessor, so DAD(1) is zero.  If there is
!    more than one connected component in the graph, then there
!    will be other nodes with DAD equal to zero.
!
!    Output, integer ( kind = 4 ) ORDER(NNODE).  ORDER(I) is the step at which
!    node I is visited in the search.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) dad(nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) order(nnode)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) kount
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lstack
  integer ( kind = 4 ) stack(nnode)

  order(1:nnode) = 0
  dad(1:nnode) = 0
  stack(1:nnode) = 0
 
  kount = 0
  i = 1
  lstack = 0
!
!  Find all descendents of the parent node I in this connected component
!  of the graph.
!
10    continue
 
  kount = kount + 1
  dad(i) = 0
  order(i) = kount
  lstack = lstack + 1
  stack(lstack) = i
!
!  Check to see if each node, J, is a "descendant" of node I.
!
30    continue
 
  j = 0
!
!  Check next neighbor, J.
!
40    continue
 
  j = j + 1

  if ( j <= nnode ) then
 
    if ( adj(i,j) /= 0 .and. order(j) == 0 ) then

      lstack = lstack + 1
      stack(lstack) = j
      dad(j) = i
      kount = kount + 1
      order(j) = kount
      i = j

      if ( kount == nnode ) then
        return
      end if

      go to 30

    end if

    go to 40

  end if
!
!  Searched all directions from current node.  Back up one node.
!
  lstack = lstack - 1
 
  if ( 0 < lstack ) then
    j = i
    i = stack(lstack)
    go to 40
  end if
!
!  The stack is exhausted.  It's time to look for another connected
!  component.
!
  lstack = 0
 
  do l = 1, nnode
    if ( order(l) == 0 ) then
      i = l
      go to 10
    end if
  end do
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GRAPH_ADJ_DFS2 - Fatal error!'

  stop
end
subroutine graph_adj_edge_count ( adj, lda, nnode, nedge )

!*****************************************************************************80
!
!! GRAPH_ADJ_EDGE_COUNT counts the number of edges in a graph.
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
!    Input, integer ( kind = 4 ) ADJ(LDA,NNODE), the adjacency information.
!    ADJ(I,J) is 1 if there is an edge from node I to node J.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the ADJ array,
!    which must be at least NNODE.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) NEDGE, the number of edges in the graph.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(lda,nnode)
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
subroutine graph_adj_eigen ( adj, lda, nnode, neigen, eigen )

!*****************************************************************************80
!
!! GRAPH_ADJ_EIGEN computes the eigenvalues of a graph from its adjacency matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(LDA,NNODE), the adjacency information.
!    ADJ(I,J) is 1 if there is an edge from node I to node J.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the ADJ array,
!    which must be at least NNODE.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) NEIGEN, the number of eigenvalues computed.
!    Normally, this would be equal to NNODE, unless the algorithm failed.
!
!    Output, real ( kind = 8 ) EIGEN(NNODE), contains the computed eigenvalues.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) nnode

  real ( kind = 8 ) a(nnode,nnode)
  integer ( kind = 4 ) adj(lda,nnode)
  real ( kind = 8 ) e(nnode)
  real ( kind = 8 ) e2(nnode)
  real ( kind = 8 ) eigen(nnode)
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) neigen

  a(1:nnode,1:nnode) = real ( adj(1:nnode,1:nnode), kind = 8 )

  call tred1 ( nnode, nnode, a, eigen, e, e2 )

  call tqlrat ( nnode, eigen, e2, ierr )

  if ( ierr == 0 ) then
    neigen = nnode
  else
    neigen = ierr - 1
  end if

  return
end
subroutine graph_adj_example_bush ( adj, lda, nnode )

!*****************************************************************************80
!
!! GRAPH_ADJ_EXAMPLE_BUSH sets up the adjacency information for the bush graph.
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
!    Output, integer ( kind = 4 ) ADJ(LDA,LDA), the adjacency information for the graph.  
!    ADJ(I,J) is 1 if nodes I and J are adjacent and 0 otherwise.  
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the ADJ array,
!    which must be at least NNODE.
!
!    Output, integer ( kind = 4 ) NNODE, the number of nodes.
!
  implicit none

  integer ( kind = 4 ) lda

  integer ( kind = 4 ) adj(lda,lda)
  integer ( kind = 4 ) nnode

  nnode = 7

  if ( lda < nnode ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRAPH_ADJ_EXAMPLE_BUSH - Fatal error!'
    write ( *, '(a)' ) '  LDA is too small.'
    stop
  end if

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
subroutine graph_adj_example_cube ( adj, lda, nnode )

!*****************************************************************************80
!
!! GRAPH_ADJ_EXAMPLE_CUBE sets up the adjacency information for the cube graph.
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
!    Output, integer ( kind = 4 ) ADJ(LDA,LDA), the adjacency information for the graph.  
!    ADJ(I,J) is 1 if nodes I and J are adjacent and 0 otherwise.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the ADJ array,
!    which must be at least NNODE.
!
!    Output, integer ( kind = 4 ) NNODE, the number of nodes.
!
  implicit none

  integer ( kind = 4 ) lda

  integer ( kind = 4 ) adj(lda,lda)
  integer ( kind = 4 ) nnode

  nnode = 8

  if ( lda < nnode ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRAPH_ADJ_EXAMPLE_CUBE - Fatal error!'
    write ( *, '(a)' ) '  LDA < NNODE.'
    write ( *, '(a,i8)' ) '  NNODE = ', nnode
    write ( *, '(a,i8)' ) '  LDA = ', lda
    stop
  end if

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
subroutine graph_adj_example_dodecahedron ( adj, lda, nnode )

!*****************************************************************************80
!
!! GRAPH_ADJ_EXAMPLE_DODECAHEDRON: adjacency info for the dodecahedron graph.
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
!    Output, integer ( kind = 4 ) ADJ(LDA,LDA), the adjacency information for 
!    the graph.  ADJ(I,J) is 1 if nodes I and J are adjacent and 0 otherwise.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the ADJ array,
!    which must be at least NNODE.
!
!    Output, integer ( kind = 4 ) NNODE, the number of nodes.
!
  implicit none

  integer ( kind = 4 ) lda

  integer ( kind = 4 ) adj(lda,lda)
  integer ( kind = 4 ) nnode

  nnode = 20

  if ( lda < nnode ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRAPH_ADJ_EXAMPLE_DODECAHEDRON - Fatal error!'
    write ( *, '(a)' ) '  LDA is too small.'
    stop
  end if

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
subroutine graph_adj_example_octo ( lda, example, seed, nnode, adj )

!*****************************************************************************80
!
!! GRAPH_ADJ_EXAMPLE_OCTO sets up an 8 node example graph.
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
!    There are 7 graphs to choose from.  They are all on 8 nodes.  The first
!    5 have degree 3 at every node.  Graphs 6 and 7 have degree 5 at every
!    node.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!

!    Input, integer ( kind = 4 ) LDA, the leading dimension of the ADJ array,
!    which must be at least NNODE.
!
!    Input, integer ( kind = 4 ) EXAMPLE, should be between 1 and 7, and indicates
!    which example graph to pick.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
!    Output, integer ( kind = 4 ) NNODE, the number of nodes, which should be 8.
!
!    Output, integer ( kind = 4 ) ADJ(LDA,LDA), the adjacency information for the graph.  
!    ADJ(I,J) is 1 if nodes I and J are adjacent and 0 otherwise.  
!
  implicit none

  integer ( kind = 4 ) lda

  integer ( kind = 4 ) adj(lda,lda)
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
    example = mod ( example - 1, 7 ) + 1
    nsave = example
  end if

  nnode = 8

  if ( lda < nnode ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRAPH_ADJ_EXAMPLE_OCTO - Fatal error!'
    write ( *, '(a)' ) '  LDA is too small.'
    stop
  end if

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
    adj(4,1) = 1

    adj(5,8) = 1
    adj(5,1) = 1
    adj(5,2) = 1

    adj(6,1) = 1
    adj(6,2) = 1
    adj(6,3) = 1

    adj(7,2) = 1
    adj(7,3) = 1
    adj(7,4) = 1

    adj(8,3) = 1
    adj(8,4) = 1
    adj(8,5) = 1

  else if ( nsave == 7 ) then

    adj(1,3) = 1
    adj(1,5) = 1
    adj(1,7) = 1

    adj(2,4) = 1
    adj(2,6) = 1
    adj(2,8) = 1

    adj(3,5) = 1
    adj(3,7) = 1
    adj(3,1) = 1

    adj(4,6) = 1
    adj(4,8) = 1
    adj(4,2) = 1

    adj(5,7) = 1
    adj(5,1) = 1
    adj(5,3) = 1

    adj(6,8) = 1
    adj(6,2) = 1
    adj(6,4) = 1

    adj(7,1) = 1
    adj(7,3) = 1
    adj(7,5) = 1

    adj(8,2) = 1
    adj(8,4) = 1
    adj(8,6) = 1

  end if
!
!  Now permute the graph.
!
  call i4mat_perm_random ( lda, nnode, seed, adj )

  return
end
subroutine graph_adj_example_twig ( adj, lda, nnode )

!*****************************************************************************80
!
!! GRAPH_ADJ_EXAMPLE_TWIG sets up the adjacency information for the twig graph.
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
!    Output, integer ( kind = 4 ) ADJ(LDA,LDA), the adjacency information for the graph.  
!    ADJ(I,J) is 1 if nodes I and J are adjacent.  
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the ADJ array,
!    which must be at least NNODE.
!
!    Output, integer ( kind = 4 ) NNODE, the number of nodes.
!
  implicit none

  integer ( kind = 4 ) lda

  integer ( kind = 4 ) adj(lda,lda)
  integer ( kind = 4 ) nnode

  nnode = 3

  if ( lda < nnode ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRAPH_ADJ_EXAMPLE_TWIG - Fatal error!'
    write ( *, '(a)' ) '  LDA is too small.'
    stop
  end if

  adj(1:nnode,1:nnode) = 0

  adj(1,2) = 1

  adj(2,1) = 1
  adj(2,3) = 1

  adj(3,2) = 1

  return
end
subroutine graph_adj_ham_cand ( adj, lda, nnode, circuit, k, nstack, &
  stack, maxstack, ncan )

!*****************************************************************************80
!
!! GRAPH_ADJ_HAM_CAND finds candidates for the next node in a Hamiltonian circuit.
!
!  Discussion:
!
!    This routine is used in conjunction with I4VEC_BACKTRACK.  
!
!  Definition:
!
!    A Hamiltonian circuit of a graph is a path that starts at a given node, 
!    visits every node exactly once, and returns to the starting node.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 August 2000
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
!    Input, integer ( kind = 4 ) ADJ(LDA,NNODE).  ADJ(I,J) = 1 if there is
!    an edge from node I to node J, 0 otherwise.
!
!    Input, integer ( kind = 4 ) LDA, the first dimension of ADJ as
!    declared in the calling program.  LDA must be at least NNODE.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes in the graph.
!
!    Input, integer ( kind = 4 ) CIRCUIT(NNODE), the nodes forming the circuit.
!
!    Input, integer ( kind = 4 ) K, index of the next node to be determined for the circuit.
!
!    Input/output, integer ( kind = 4 ) NSTACK, the current length of stack.
!
!    Input, integer ( kind = 4 ) STACK(MAXSTACK), candidates for steps 1...K-1.
!
!    Input, integer ( kind = 4 ) MAXSTACK, the dimension of STACK.
!
!    Workspace, integer NCAN(NNODE), the number of candidates for
!    positions in the circuit.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) nnode
  integer ( kind = 4 ) maxstack

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) circuit(nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iwork(nnode)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) ncan(nnode)
  integer ( kind = 4 ) nstack
  integer ( kind = 4 ) stack(maxstack)

  ncan(k) = 0

  if ( k == 1 ) then
    stack(1) = 1
    nstack = 1
    ncan(k) = 1
    return
  end if
 
  iwork(1:nnode) = adj(circuit(k-1),1:nnode)
 
  iwork(circuit(1:k-1)) = 0
  
  if ( k /= nnode ) then
 
    do i = 1, nnode
      if ( iwork(i) == 1 ) then
        if ( maxstack <= nstack ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'GRAPH_ADJ_HAM_CAND - Fatal error!'
          write ( *, '(a)' ) '  Stack size exceeded.'
          stop
        end if
        nstack = nstack + 1
        stack(nstack) = i
        ncan(k) = ncan(k) + 1
      end if
    end do
 
    return
 
 else if ( k == nnode ) then
 
    do i = 1, nnode
 
      if ( iwork(i) == 1 ) then
 
        if ( circuit(2) < i .or. adj(i,1) == 0 ) then

        else
          if ( maxstack <= nstack ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'GRAPH_ADJ_HAM_CAND - Fatal error!'
            write ( *, '(a)' ) '  Stack size exceeded.'
            stop
          end if
          nstack = nstack + 1
          stack(nstack) = i
          ncan(k) = ncan(k) + 1
        end if

        return
 
      end if
 
    end do

  end if
 
  return
end
subroutine graph_adj_ham_next ( adj, lda, nnode, circuit, stack, maxstack, &
  ncan, more )

!*****************************************************************************80
!
!! GRAPH_ADJ_HAM_NEXT returns the next Hamilton circuit for a graph.
!
!  Discussion:
!
!    The routine produces all the Hamilton circuits of a graph, one at a time.
!
!    A Hamiltonian circuit of a graph is a path that starts at a given
!    node, visits every node exactly once, and returns to the starting node.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 August 2000
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
!    Input, integer ( kind = 4 ) ADJ(LDA,NNODE).  ADJ(I,J) = 1 if there is
!    an edge from node I to node J, 0 otherwise.
!
!    Input, integer ( kind = 4 ) LDA, the first dimension of ADJ as
!    declared in the calling program.  LDA must be at least NNODE.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes in graph.
!
!    Input, integer ( kind = 4 ) CIRCUIT(NNODE).  CIRCUIT(I) is the I-th node
!    in the circuit. 
!
!    Input, integer ( kind = 4 ) K, the index of the next node to be determined
!    for circuit.
!
!    Input/output, integer ( kind = 4 ) NSTACK, the current length of stack.
!
!    Input, integer ( kind = 4 ) STACK(MAXSTACK).  Candidates for steps 1...K-1.
!
!    Input, integer ( kind = 4 ) MAXSTACK, dimension of STACK.
!
!    Workspace, integer NCAN(NNODE), the number of candidates for each
!    position in the circuit.
!
!    Input/output, logical MORE.
!    On first call, set MORE to .FALSE, and do not alter it after.
!    On return, MORE is TRUE if another circuit has been returned in
!    IARRAY, and FALSE if there are no more circuits.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) nnode
  integer ( kind = 4 ) maxstack

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) circuit(nnode)
  integer ( kind = 4 ), save :: indx = 0
  integer ( kind = 4 ), save :: k = 0
  logical more
  integer ( kind = 4 ) ncan(nnode)
  integer ( kind = 4 ), save :: nstack = 0
  integer ( kind = 4 ) stack(maxstack)

  if ( .not. more ) then
    indx = 0
    k = 0
    more = .true.
    nstack = 0
  end if
 
  do
 
    call i4vec_backtrack ( nnode, circuit, indx, k, nstack, stack, maxstack, &
      ncan )
 
    if ( indx == 1 ) then

      exit

    else if ( indx == 2 ) then

      call graph_adj_ham_cand ( adj, lda, nnode, circuit, k, nstack, &
        stack, maxstack, ncan )

    else

      more = .false.
      exit

    end if

  end do
 
  return
end
subroutine graph_adj_ham_next_brute ( adj, lda, nnode, circuit, iset )

!*****************************************************************************80
!
!! GRAPH_ADJ_HAM_NEXT_BRUTE finds the next Hamiltonian circuit in a graph.
!
!  Discussion:
!
!    This is a brute force algorithm, and not suitable for large problems.
!    It is really only useful as a demonstration, and as a check for
!    the backtracking algorithm.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(LDA,NNODE), the adjacency information.
!    ADJ(I,J) is 1 if there is a direct link between nodes I and J.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of ADJ, which must be
!    at least NNODE.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input/output, integer ( kind = 4 ) CIRCUIT(NNODE).
!
!    On input, if ISET = 0, then CIRCUIT is not presumed to contain any 
!    information.  If ISET is nonzero, then CIRCUIT contains the circuit 
!    computed on the previous call.
!
!    On output, CIRCUIT contains the circuit computed by this call.
!
!    Input/output, integer ( kind = 4 ) ISET.
!    On input, 0 means this is the first call for this graph.  
!    Any other value means this is a repeated call for more circuits.
!
!    On output, a 0 value means that no more circuits could be computed.
!    Otherwise, ISET is incremented by one on each call.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) circuit(nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ipos
  integer ( kind = 4 ) iset
!
!  If ISET is 0, this is a starting call, and we set CIRCUIT
!  to the lexically first circuit to check.
!
!  Otherwise, we set CIRCUIT to the next permutation.
!
  if ( iset == 0 ) then
    ipos = 0
    circuit(1:nnode) = 0
  else
    ipos = nnode - 1
  end if
 
10 continue
 
  call perm_inc ( circuit, ipos, nnode )

  if ( ipos <= 0 .or. circuit(1) /= 1 ) then
    iset = 0
    circuit(1:nnode) = 0
    return
  end if
!
!  Check whether the entries of CIRCUIT actually form a circuit.
!  If we find a break in the circuit, store that location in IPOS
!  and move on to try the next permutation.
!
  do i = 1, nnode-1
    ipos = i
    if ( adj(circuit(i),circuit(i+1)) == 0 ) then
      go to 10
    end if
  end do
!
!  If the circuit connects all the nodes, we only have to check whether
!  the last node connects back to the first one.
!
!  To cut down the pairs of equivalent circuits created by going one
!  way or the other over the same set of nodes, we also require that,
!  for 2 < NNODE, the last node be numbered higher than the second one.
!
  if ( adj(circuit(nnode),circuit(1)) == 0 ) then
    go to 10
  end if

  if ( 2 < nnode ) then
    if ( circuit(nnode) < circuit(2) ) then
      go to 10
    end if
  end if
 
  iset = iset + 1

  return
end
subroutine graph_adj_is_bipartite ( adj, lda, nnode, result )

!*****************************************************************************80
!
!! GRAPH_ADJ_IS_BIPARTITE determines if a graph is bipartite.
!
!  Definition:
!
!    A graph is bipartite if its nodes can be divided into two subsets
!    in such a way that every edge joins a node from each subset.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 October 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(LDA,NNODE), the adjacency matrix for the graph.
!    ADJ(I,J) is nonzero if there is an edge from node I to node J.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of LDA, which must be
!    at least NNODE.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) RESULT.
!    0, the graph is not bipartite.
!    1, the graph is bipartite.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) khi
  integer ( kind = 4 ) klo
  integer ( kind = 4 ) lhi
  integer ( kind = 4 ) list(nnode)
  integer ( kind = 4 ) oldset
  integer ( kind = 4 ) result
  integer ( kind = 4 ) set
  integer ( kind = 4 ) subset(nnode)

  result = 1
  subset(1:nnode) = -1
!
!  Node 1 is put in subset 1.
!
  set = 1
  list(1) = 1
  subset(1) = set

  klo = 1
  khi = 1
!
!  Working from the set of nodes found on the previous step, look
!  for all in and out neighbors.  
!
10    continue

  oldset = set
  set = 1 - set

  lhi = khi
!
!  Consider each node I in the previously found set.
!
  do k = klo, khi

    i = list(k)
!
!  Look at all in and out neighbors J.
!
    do j = 1, nnode

      if ( adj(i,j) /= 0 .or. adj(j,i) /= 0 ) then
!
!  If the node is not in any subset, put it in the other one.
!
        if ( subset(j) == -1 ) then

          lhi = lhi + 1
          list(lhi) = j
          subset(j) = set
!
!  But if the node is in the same subset, bipartiteness has failed.
!
        else if ( subset(j) == oldset ) then
          result = 0
          return
        end if

      end if

    end do

  end do
!
!  Assuming we found more nodes, on this sweep, then ...
!
  if ( khi < lhi ) then
    klo = khi + 1
    khi = lhi 
    go to 10
  end if
!
!  Assuming we found no new nodes on this sweep, see if there are any
!  nodes we have missed.  These will be completely isolated from all the
!  nodes we have found so far.
!
  do i = 1, nnode

    if ( subset(i) == -1 ) then
      klo = khi + 1
      khi = klo
      subset(i) = set
      list(klo) = i
      go to 10
    end if

  end do

  result = 1

  return
end
subroutine graph_adj_is_edge_connected ( adj, lda, nnode, result )

!*****************************************************************************80
!
!! GRAPH_ADJ_IS_EDGE_CONNECTED determines if a graph is edgewise connected.
!
!  Definition:
!
!    A graph is edgewise connected if from any edge it is possible to reach
!    any other edge.  An edgewise connected graph may include isolated nodes.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 October 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(LDA,NNODE), the adjacency matrix for the graph.  
!    ADJ(I,J) is nonzero if there is an edge from node I to node J.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of LDA, which must be
!    at least NNODE.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) RESULT.
!    0, the graph is not edgewise connected.
!    1, the graph is edgewise connected.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) found(nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) list(nnode)
  integer ( kind = 4 ) result
!
!  FOUND(I) is 1 if edge I has been reached.
!  LIST(I) contains a list of the nodes as they are reached.
!
  list(1:nnode) = 0
  found(1:nnode) = 0
!
!  Find an edge.
!
  ilo = 1
  ihi = 0

  do i = 1, nnode
    do j = 1, nnode

      if ( 0 < adj(i,j) .or. 0 < adj(j,i) ) then

        adj(i,j) = - abs ( adj(i,j) )
        adj(j,i) = - abs ( adj(j,i) )

        ihi = ihi + 1
        list(ihi) = i
        found(i) = 1
        if ( i /= j ) then
          ihi = ihi + 1
          list(ihi) = j
          found(j) = 1
        end if
        go to 10

      end if

    end do
  end do
!
!  A graph with NO edges is edgewise connected!
!
  result = 1
  return

10    continue
!
!  From the batch of edge nodes found last time, LIST(ILO:IHI),
!  look for unfound neighbors, and store their indices in LIST(JLO:JHI).
!
  jlo = ihi + 1
  jhi = ihi

  do ii = ilo, ihi

    i = list(ii)

    do j = 1, nnode

      if ( 0 < adj(i,j) ) then

        adj(i,j) = - adj(i,j)
        if ( 0 < adj(j,i) ) then
          adj(j,i) = - adj(j,i)
        end if

        if ( found(j) == 0 ) then
          jhi = jhi + 1
          list(jhi) = j
          found(j) = 1
        end if

      end if

    end do

  end do
!
!  If any neighbors were found, go back and find THEIR neighbors.
!
  if ( jlo <= jhi ) then
    ilo = jlo
    ihi = jhi
    go to 10
  end if
!
!  If any edges were unvisited, then the graph is not edgewise connected.
!
  result = 1
  do i = 1, nnode
    do j = 1, nnode
      if ( 0 < adj(i,j) ) then
        result = 0
      end if
    end do
  end do
!
!  Restore the positive sign of ADJ.
!
  adj(1:nnode,1:nnode) = abs ( adj(1:nnode,1:nnode) )

  return
end
subroutine graph_adj_is_eulerian ( adj, lda, nnode, result )

!*****************************************************************************80
!
!! GRAPH_ADJ_IS_EULERIAN determines if a graph is Eulerian from its adjacency matrix.
!
!  Discussion:
!
!    A graph is path-Eulerian if there exists a path through the graph
!    which uses every edge once.
!
!    A graph is circuit-Eulerian if there exists a path through the graph
!    which uses every edge once, and which starts and ends on the same node.
!
!    Note that it is NOT necessary for the path or circuit to pass through
!    every node; simply that all the edges can be used exactly once to
!    make a connected path.  This means an Eulerian graph can have isolated
!    nodes, for instance.
!
!    A graph is path-Eulerian if and only if it is edge connected, and all 
!    but two nodes are of even degree.
!
!    A graph is circuit-Eulerian if and only if it is edge connected and
!    all nodes are of even degree.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(LDA,NNODE), the adjacency matrix for the graph.  
!    ADJ(I,J) is nonzero if there is an edge from node I to node J.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of LDA, which must be
!    at least NNODE.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) RESULT.
!    0, the graph is not Eulerian.
!    1, the graph is path-Eulerian.
!    2, the graph is circuit-Eulerian.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) degree
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nodd
  integer ( kind = 4 ) result
!
!  First check that the graph is edgewise connected.
!
  call graph_adj_is_edge_connected ( adj, lda, nnode, result )

  if ( result == 0 ) then
    return
  end if
!
!  Now look at node degree.
!
  nodd = 0

  do i = 1, nnode

    degree = 0

    do j = 1, nnode
      if ( adj(i,j) /= 0 ) then
        if ( i == j ) then
          degree = degree + 2
        else
          degree = degree + 1
        end if
      end if
    end do

    if ( mod ( degree, 2 ) == 1 ) then
      nodd = nodd + 1
    end if

  end do

  if ( nodd == 0 ) then
    result = 2
  else if ( nodd == 2 ) then
    result = 1
  else
    result = 0
  end if

  return
end
subroutine graph_adj_is_node_connected ( adj, lda, nnode, result )

!*****************************************************************************80
!
!! GRAPH_ADJ_IS_NODE_CONNECTED determines if a graph is nodewise connected.
!
!  Definition:
!
!    A graph is nodewise connected if, from every node, there is a path
!    to any other node.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 October 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(LDA,NNODE), the adjacency matrix for the graph.  
!    ADJ(I,J) is nonzero if there is an edge from node I to node J.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of LDA, which must be
!    at least NNODE.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) RESULT.
!    0, the graph is not nodewise connected.
!    1, the graph is nodewise connected.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) found(nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) list(nnode)
  integer ( kind = 4 ) result
!
!  FOUND(I) is 1 if node I has been reached.
!  LIST(I) contains a list of the nodes as they are reached.
!
  list(1:nnode) = 0
  found(1:nnode) = 0
!
!  Start at node 1.
!
  found(1) = 1
  list(1) = 1
  ilo = 1
  ihi = 1
!
!  From the batch of nodes found last time, LIST(ILO:IHI),
!  look for unfound neighbors, and store their indices in LIST(JLO:JHI).
!
10    continue

  jlo = ihi + 1
  jhi = ihi

  do ii = ilo, ihi

    i = list(ii)

    do j = 1, nnode

      if ( adj(i,j) /= 0 .or. adj(j,i) /= 0 ) then

        if ( found(j) == 0 ) then
          jhi = jhi + 1
          list(jhi) = j
          found(j) = 1
        end if

      end if

    end do

  end do
!
!  If any neighbors were found, go back and find THEIR neighbors.
!
  if ( jlo <= jhi ) then
    ilo = jlo
    ihi = jhi
    go to 10
  end if
!
!  No more neighbors were found.  Have we reached all nodes?
!
  if ( ihi == nnode ) then
    result = 1
  else
    result = 0
  end if

  return
end
subroutine graph_adj_is_tree ( adj, lda, nnode, result )

!*****************************************************************************80
!
!! GRAPH_ADJ_IS_TREE determines whether a graph is a tree.
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
!    Input, integer ( kind = 4 ) ADJ(LDA,NNODE), the adjacency matrix for the graph. 
!    ADJ(I,J) is nonzero if there is an edge from node I to node J.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of LDA, which must be
!    at least NNODE.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) RESULT.
!    0, the graph is not a tree.
!    1, the graph is a tree.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) result

  if ( nnode <= 1 ) then
    result = 1
    return
  end if
!
!  Every node must be connected to every other node.
!
  call graph_adj_is_node_connected ( adj, lda, nnode, result )

  if ( result == 0 ) then
    return
  end if
!
!  There must be exactly NNODE-1 edges.
!
  call graph_adj_edge_count ( adj, lda, nnode, nedge )

  if ( nedge == nnode - 1 ) then
    result = 1
  else
    result = 0
  end if

  return
end
subroutine graph_adj_print ( adj, lda, nnode, title )

!*****************************************************************************80
!
!! GRAPH_ADJ_PRINT prints out an adjacency matrix for a graph.
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
!    04 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(LDA,NNODE), the adjacency matrix of a graph.  
!    ADJ(I,J) is 1 if there is a direct connection FROM node I TO node J,
!    and is 0 otherwise.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of ADJ, which must be
!    at least NNODE.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.  
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  character ( len = 80 ) string
  character ( len = * ) title

  if ( len_trim ( title ) /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  do i = 1, nnode

    jhi = min ( nnode, 74 )

    do j = 1, jhi

      if ( 0 <= adj(i,j) .and. adj(i,j) <= 9 ) then
        string(j:j) = char ( 48 + adj(i,j) )
      else
        string(j:j) = '*'
      end if

    end do

    write ( *, '(2x,i3,2x,a)' ) i, string(1:jhi)

  end do

  return
end
subroutine graph_adj_random ( nnode, nedge, seed, adj )

!*****************************************************************************80
!
!! GRAPH_ADJ_RANDOM generates a random graph on NNODE nodes with NEDGE edges.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 September 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) NEDGE, the number of edges, which must be between
!    0 and (NNODE*(NNODE-1))/2.  (Note that each edge will be listed
!    twice in the adjacency matrix).
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
!    Output, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency matrix.  ADJ(I,J) is
!    nonzero if there is an edge from node I to node J.  ADJ(I,I) will
!    always be 0.
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
  if ( nnode <= 0  ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRAPH_ADJ_RANDOM - Fatal error!'
    write ( *, '(a,i8)' ) '  NNODE = ', nnode
    write ( *, '(a)' ) '  but NNODE must be at least 1.'
    stop
  end if

  maxedge = ( nnode * ( nnode - 1 ) ) / 2

  if ( nedge < 0 .or. maxedge < nedge ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRAPH_ADJ_RANDOM - Fatal error!'
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
  call ksub_random ( maxedge, nedge, seed, iwork )
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
subroutine graph_adj_random2 ( nnode, prob, seed, nedge, adj )

!*****************************************************************************80
!
!! GRAPH_ADJ_RANDOM2 generates a random graph on NNODE nodes with NEDGE edges.
!
!  Discussion:
!
!    The user specifies the probability P that an edge will be generated
!    between any pair of nodes.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 September 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, real ( kind = 8 ) PROB, the probability that an edge will
!    be generated between any given pair of nodes.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
!    Output, integer ( kind = 4 ) NEDGE, the number of edges between distinct pairs
!    of nodes, which must be between 0 and (NNODE*(NNODE-1))/2.
!
!    Output, integer ( kind = 4 ) ADJ(NNODE,NNODE), the adjacency matrix.  ADJ(I,J) is
!    nonzero if there is an edge from node I to node J.  ADJ(I,I) will
!    always be 1.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nedge
  real ( kind = 8 ) p(nnode)
  real ( kind = 8 ) prob
  integer ( kind = 4 ) seed

  nedge = 0

  adj(1:nnode,1:nnode) = 0

  do i = 1, nnode
    adj(i,i) = 1
  end do

  do i = 1, nnode-1

    call r8vec_uniform_01 ( nnode-i, seed, p(i+1:nnode) )

    do j = i+1, nnode
      if ( p(j) <= prob ) then
        nedge = nedge + 1
        adj(i,j) = 1
        adj(j,i) = 1
      end if
    end do

  end do

  return
end
subroutine graph_adj_reduce ( adj, nnode )

!*****************************************************************************80
!
!! GRAPH_ADJ_REDUCE generates a transitive reduction of a graph.
!
!  Discussion:
!
!    This routine is given an adjacency matrix B, which might be a
!    transitive closure of a graph G.
!
!    The transitive closure graph is generated from a graph G by the 
!    following procedure:
!
!      B(I,J) = 0 if node J cannot be reached from node I in graph G;
!               1 if node J can be reached from node I in graph G.
!
!    The purpose of this routine is to try to find the original, sparser
!    graph G which generated the given transitive closure graph.  Such a
!    graph G is known as a transitive reduction..  In general,
!    there is no unique solution.  In particular, any graph is a transitive
!    reduction of itself.  
!
!    Hence, the real task is to drop as many redundant edges as possible
!    from the given graph, arriving at a graph from which no more edges 
!    may be removed.
!
!  Method:
!
!    One way of explaining the algorithm is based on the adjacency matrix:
!
!    * Zero out the diagonals of the adjacency matrix.
!
!    * Consider row 1.  Any other row that can "reach" row 1 doesn't
!      need a 1 if row 1 has it.  So "subtract" all the 1's in row 1
!      from such rows.  We are done with row 1 and column 1.
!
!    * Repeat for the other rows.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 September 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) ADJ(NNODE,NNODE).
!    On input, the adjacency matrix of the transitive closure graph H.
!    On output, the adjacency matrix of a transitive reduction graph G 
!    of the graph H.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(nnode,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
!
!  First discard those useless self-edges.
!
  do i = 1, nnode
    adj(i,i) = 0
  end do
!
!  If you can get from J to I and I to K, you don't need an
!  edge from J to K.
!
  do i = 1, nnode
    do j = 1, nnode
      if ( adj(j,i) /= 0 .or. adj(i,j) /= 0 ) then
        do k = 1, nnode
          if ( adj(i,k) /= 0 .or. adj(k,i) /= 0 ) then
            adj(j,k) = 0
            adj(k,j) = 0
          end if
        end do
      end if
    end do
  end do

  return
end
subroutine graph_adj_span_tree ( adj, lda, nnode, inode, jnode )

!*****************************************************************************80
!
!! GRAPH_ADJ_SPAN_TREE finds a spanning tree of a graph.
!
!  Discussion:
!
!    If the graph is connected, NNODE-1 edges comprise the spanning tree.  
!
!    If the graph is not connected, but divided into NCOMP components, then 
!    NNODE-NCOMP edges will comprise the spanning "forest", and the other 
!    edges will be zero.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 October 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(LDA,NNODE), the adjacency matrix for the graph.  
!    ADJ(I,J) is 0 if there is no edge from node I to node J.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of ADJ, which must
!    be at least NNODE.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) INODE(NNODE-1), JNODE(NNODE-1), the edge list for the
!    spanning tree or forest.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) inode(nnode-1)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jnode(nnode-1)
  integer ( kind = 4 ) label(nnode)
  integer ( kind = 4 ) level
  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) nfound
  integer ( kind = 4 ) nlabel

  label(1:nnode) = 0

  inode(1:nnode-1) = 0
  jnode(1:nnode-1) = 0

  level = 0
  nedge = 0
  nlabel = 0
!
!  Find an unvisited node.
!
  do

    i = 0
 
    do

      i = i + 1

      if ( label(i) == 0 ) then
        exit
      end if

    end do
  
    label(i) = level + 1
    nlabel = nlabel + 1
!
!  Search for all nodes reachable from the node.
!
    do

      level = level + 1
      nfound = 0

      do i = 1, nnode

        if ( label(i) == level ) then

          do j = 1, nnode

            if ( label(j) == 0 ) then

              if ( adj(i,j) /= 0 .or. adj(j,i) /= 0 ) then
                label(j) = level + 1
                nlabel = nlabel + 1
                nfound = nfound + 1
                nedge = nedge + 1
                inode(nedge) = i
                jnode(nedge) = j
              end if

            end if

          end do

        end if

      end do

      if ( nfound <= 0 ) then
        exit
      end if

    end do
!
!  If we have labeled all nodes, exit.
!
    if ( nnode <= nlabel ) then
      exit
    end if

  end do

  return
end
subroutine graph_adj_span_tree_enum ( adj, lda, nnode, tree_num )

!*****************************************************************************80
!
!! GRAPH_ADJ_SPAN_TREE_ENUM enumerates the spanning trees of a graph.
!
!  Discussion:
!
!    If ADJ is the adjacency matrix of the graph, let A be the matrix
!      A = DEG - ADJ
!    where DEG is the diagonal matrix with DEG(I,I) = degree of node I.
!    Then the number of spanning trees of the graph is equal to the
!    determinant of any cofactor of A.  A cofactor of A is obtained by
!    deleting a row and column of A.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ADJ(LDA,NNODE), the adjacency matrix for the graph.
!    ADJ(I,J) is 0 if there is no edge from node I to node J.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of ADJ, which must
!    be at least NNODE.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) TREE_NUM, the number of spanning trees.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) nnode

  real ( kind = 8 ) a(nnode,nnode)
  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) degree(nnode)
  real ( kind = 8 ) det
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipivot(nnode)
  integer ( kind = 4 ) tree_num
!
!  Construct the matrix.
!
  call graph_adj_degree ( adj, lda, nnode, degree )

  a(1:nnode,1:nnode) = - real ( adj(1:nnode,1:nnode), kind = 8 )

  do i = 1, nnode
    a(i,i) = a(i,i) + real ( degree(i), kind = 8 )
  end do
!
!  Factor the NNODE-1 order matrix.
!
  call dge_fa ( nnode, nnode-1, a, ipivot, info )

  if ( info /= 0 ) then
    tree_num = 0
    return
  end if
!
!  Get the determinant.
!
  call dge_det ( nnode, nnode-1, a, ipivot, det )

  tree_num = nint ( det )

  return
end
subroutine graph_adj_symmetrize ( adj, lda, nnode )

!*****************************************************************************80
!
!! GRAPH_ADJ_SYMMETRIZE symmetrizes an adjacency matrix.
!
!  Discussion:
!
!    For a graph, if there is an edge from I to J, there is an edge from
!    J to I.  Therefore, the adjacency matrix should be symmetric.  
!    This routine enforces that condition.  If either ADJ(I,J) or ADJ(J,I)
!    is nonzero, the output adjacency matrix will have both entries nonzero.
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
!    Input/output, integer ( kind = 4 ) ADJ(LDA,NNODE).  On output, the 
!    adjacency information has been symmetrized.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of ADJ, which must 
!    be NNODE or greater.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
!
!  While not perfect, this method does not assume that 1 is the only
!  legal nonzero value in ADJ.
!
  do i = 1, nnode
    do j = i+1, nnode
      if ( adj(i,j) /= 0 ) then
        adj(j,i) = adj(i,j)
      else if ( adj(j,i) /= 0 ) then
        adj(i,j) = adj(j,i)
      end if
    end do
  end do
 
  return
end
subroutine graph_adj_to_graph_arc ( adj, lda, nnode, maxedge, nedge, inode, &
  jnode )

!*****************************************************************************80
!
!! GRAPH_ADJ_TO_GRAPH_ARC converts an adjacency graph to an arc list graph.
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
!    Input, integer ( kind = 4 ) ADJ(LDA,NNODE), the adjacency matrix for the graph.  
!    ADJ(I,J) is nonzero if there is an edge from node I to node J.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of LDA, which must be
!    at least NNODE.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) MAXEDGE, the maximum number of edges.
!
!    Output, integer ( kind = 4 ) NEDGE, the number of edges.
!
!    Output, integer ( kind = 4 ) INODE(MAXEDGE), JNODE(MAXEDGE), the arc list of the
!    graph.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) maxedge
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) inode(maxedge)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jnode(maxedge)
  integer ( kind = 4 ) nedge

  nedge = 0

  inode(1:maxedge) = 0
  jnode(1:maxedge) = 0

  do j = 1, nnode
    do i = j, nnode
      if ( adj(i,j) /= 0 .or. adj(j,i) /= 0 ) then
        nedge = nedge + 1
        if ( nedge <= maxedge ) then
          inode(nedge) = i
          jnode(nedge) = j
        else
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'GRAPH_ADJ_TO_GRAPH_ARC - Fatal error!'
          write ( *, '(a)' ) '  MAXEDGE exceeded.'
          stop
        end if
      end if
    end do
  end do

  return
end
subroutine graph_arc_complement ( inode, jnode, inode2, jnode2, maxedge, &
  nedge, nedge2, nnode )

!*****************************************************************************80
!
!! GRAPH_ARC_COMPLEMENT returns the edge list of the complement of a graph.
!
!  Discussion:
!
!    This routine can also handle a directed graph.
!
!  Definition:
!
!    The complement of a graph G is a graph H with the property that
!    nodes U and V are connected in H if and only if they are not
!    connected in G.  However, edges from a node to itself are not
!    allowed.
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
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) INODE(NEDGE), JNODE(NEDGE).  INODE(I) and JNODE(I) 
!    are the start and end nodes of the I-th edge of the graph G.  On
!    output, the data in INODE and JNODE will have been sorted, but not
!    otherwise disrupted.
!
!    Output, integer ( kind = 4 ) INODE2(MAXEDGE), JNODE2(MAXEDGE).  INODE2(I) and JNODE2(I) 
!    are the start and end nodes of the I-th edge of the complement graph H.
!
!    Input, integer ( kind = 4 ) MAXEDGE, the amount of storage available in INODE2
!    and JNODE2.  MAXEDGE only needs to be as large as NEDGE2, and NEDGE2
!    can be precomputed, assuming that the input value of NEDGE does not
!    count any self edges (edges from a node to itself), and does not
!    count an edge twice (that is, counting the edge from I to J, and
!    the edge from J to I, as distinct).  If this is so, then you can
!    set MAXEDGE = NEDGE2 = 0.5 * ( NNODE * ( NNODE - 1 ) ) - NEDGE.
!
!    Input, integer ( kind = 4 ) NEDGE, the number of edges in the graph G.
!
!    Output, integer ( kind = 4 ) NEDGE2, the number of edges in the complement graph H.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
  implicit none

  integer ( kind = 4 ) maxedge
  integer ( kind = 4 ) nedge

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) inedge
  integer ( kind = 4 ) inode(nedge)
  integer ( kind = 4 ) inode2(maxedge)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) jnode(nedge)
  integer ( kind = 4 ) jnode2(maxedge)
  integer ( kind = 4 ) nedge2
  integer ( kind = 4 ) nnode
!
!  Sort the input edge array.
!
  call graph_arc_edge_sort ( nedge, inode, jnode )
!
!  Compute the complementary edges.
!
  nedge2 = 0
 
  inedge = 0
  i2 = 1
  j2 = 1

  do while ( inedge < nedge )
 
    inedge = inedge + 1
    i1 = i2
    j1 = j2

    if ( inedge <= nedge ) then
      i2 = inode(inedge)
      j2 = jnode(inedge)
    else
      i2 = nnode
      j2 = nnode
    end if
 
    if ( i1 == i2 ) then
 
      do j = j1+1, j2-1
        if ( i1 < j ) then
          nedge2 = nedge2 + 1
          inode2(nedge2) = i2
          jnode2(nedge2) = j
        end if
      end do
 
    else
 
      do j = j1+1, nnode
        if ( i1 < j ) then
          nedge2 = nedge2 + 1
          inode2(nedge2) = i1
          jnode2(nedge2) = j
        end if
      end do
 
      do i = i1+1, i2-1
        do j = 1, nnode
          if ( i < j ) then
            nedge2 = nedge2 + 1
            inode2(nedge2) = i
            jnode2(nedge2) = j
          end if
        end do
      end do

      do j = 1, j2-1
        if ( i2 < j ) then
          nedge2 = nedge2 + 1
          inode2(nedge2) = i2
          jnode2(nedge2) = j
        end if
      end do
 
    end if

  end do
 
  return
end
subroutine graph_arc_degree ( nnode, nedge, inode, jnode, degree )

!*****************************************************************************80
!
!! GRAPH_ARC_DEGREE determines the degree of the nodes of a graph.
!
!  Definition:
!
!    The degree of a node is the number of edges that include the node.
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
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) NEDGE, the number of edges.
!
!    Input, integer ( kind = 4 ) INODE(NEDGE), JNODE(NEDGE), the pairs of nodes
!    that form the edges.
!
!    Output, integer ( kind = 4 ) DEGREE(NNODE), the degree of each node, that is,
!    the number of edges that include the node.
!
  implicit none

  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) degree(nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) inode(nedge)
  integer ( kind = 4 ) jnode(nedge)
  integer ( kind = 4 ) n

  degree(1:nnode) = 0

  do i = 1, nedge

    n = inode(i)
    if ( 1 <= n .and. n <= nnode ) then
      degree(n) = degree(n) + 1
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'GRAPH_ARC_DEGREE - Fatal error!'
      write ( *, '(a,i8)' ) '  Out-of-range node value = ', n
      stop
    end if

    n = jnode(i)
    if ( 1 <= n .and. n <= nnode ) then
      degree(n) = degree(n) + 1
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'GRAPH_ARC_DEGREE - Fatal error!'
      write ( *, '(a,i8)' ) '  Out-of-range node value = ', n
      stop
    end if

  end do

  return
end
subroutine graph_arc_edge_con2 ( nnode, nedge, inode, jnode, edge_con )

!*****************************************************************************80
!
!! GRAPH_ARC_EDGE_CON2 finds the edge-connectivity of a connected graph.
!
!  Method:
!
!    A graph G has edge connectivity K if, given any pair of distinct nodes 
!    I and J, there are K paths from I to J, no two of which use a common edge.
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
!    The edge-connectivity is then equal to the minimum of G(2:NNODE).
!
!    This routine finds the edge-connectivity of a given undirected graph with 
!    the help of a maximum flow algorithm.
!
!    The maximum network flow algorithm requires O(NNODE**3) operations.  The 
!    edge-connectivity of a graph will therefore be found in O(NNODE**4) 
!    operations.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 July 2000
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
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) NEDGE, the number of edges.
!
!    Input, integer ( kind = 4 ) INODE(NEDGE), JNODE(NEDGE), the end nodes of the edges.
!
!    Output, integer ( kind = 4 ) EDGE_CON, the edge-connectivity of the graph.
!
  implicit none

  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) capflo(2,2*nedge)
  integer ( kind = 4 ) edge_con
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icut(nnode)
  integer ( kind = 4 ) iendpt(2,2*nedge)
  integer ( kind = 4 ) inode(nedge)
  integer ( kind = 4 ) isink
  integer ( kind = 4 ) isorce
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jnode(nedge)
  integer ( kind = 4 ) node_flow(nnode)
!
!  Create the network from the graph.
!
  j = 0
 
  do i = 1, nedge

    j = j + 1
    iendpt(1,j) = inode(i)
    iendpt(2,j) = jnode(i)
    capflo(1,j) = 1
    capflo(2,j) = 0

    j = j + 1
    iendpt(1,j) = jnode(i)
    iendpt(2,j) = inode(i)
    capflo(1,j) = 1
    capflo(2,j) = 0

  end do
!
!  Call the network flow algorithm.
!
  edge_con = nnode
  isorce = 1

  do isink = 2, nnode

    call network_flow_max ( nnode, 2*nedge, iendpt, capflo, isorce, isink, &
      icut, node_flow )
 
    if ( node_flow(isorce) < edge_con ) then
      edge_con = node_flow(isorce)
    end if

  end do
 
  return
end
subroutine graph_arc_edge_sort ( nedge, inode, jnode )

!*****************************************************************************80
!
!! GRAPH_ARC_EDGE_SORT sorts the edge array of a graph.
!
!  Discussion:
!
!    The pair of nodes (I,J) representing an edge is reordered so
!    that the smaller node is listed first.
!
!    Then the edges are sorted in dictionary order.
!
!  Example:
!
!    Input:
!
!      INODE  JNODE
!
!        3      2
!        4      3
!        2      1
!        1      4
!
!    Output:
!
!      INODE  JNODE
!
!        1      2
!        1      4
!        2      3
!        3      4
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
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEDGE, the number of edges.
!
!    Input/output, integer ( kind = 4 ) INODE(NEDGE), JNODE(NEDGE), the edge array of a
!    graph.  The I-th edge of the graph connects nodes INODE(I) and
!    JNODE(I).
!
!    On output, the INODE and JNODE arrays have been sorted as described.
!
  implicit none

  integer ( kind = 4 ) nedge

  integer ( kind = 4 ) i
  integer ( kind = 4 ) iedge
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) inode(nedge)
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) jedge
  integer ( kind = 4 ) jnode(nedge)

  if ( nedge <= 1 ) then
    return
  end if
!
!  Sort the node pairs.
!
  do i = 1, nedge
    if ( jnode(i) < inode(i) ) then
      call i4_swap ( inode(i), jnode(i) )
    end if
  end do
!
!  Sort the edges using an external heap sort.
!
  iedge = 0
  jedge = 0
  indx = 0
  isgn = 0

  do

    call sort_heap_external ( nedge, indx, iedge, jedge, isgn )
!
!  Interchange edges IEDGE and JEDGE.
!
    if ( 0 < indx ) then

      call i4_swap ( inode(iedge), inode(jedge) )
      call i4_swap ( jnode(iedge), jnode(jedge) )
!
!  Compare edges IEDGE and JEDGE.
!
    else if ( indx < 0 ) then

      if ( ( inode(iedge) < inode(jedge) ) .or. &
         ( inode(iedge) == inode(jedge) .and. &
           jnode(iedge) < jnode(jedge) ) ) then
        isgn = -1
      else
        isgn = +1
      end if

    else

      exit

    end if

  end do
 
  return
end
subroutine graph_arc_euler_circ ( nnode, nedge, inode, jnode, circuit )

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
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 July 2000
!
!  Reference:
!
!    Hang Tong Lau,
!    Combinatorial Heuristic Algorithms in FORTRAN,
!    Springer Verlag, 1986.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes in the graph.
!
!    Input, integer ( kind = 4 ) NEDGE, the number of edges in the graph.
!
!    Input, integer ( kind = 4 ) INODE(NEDGE), JNODE(NEDGE), the two end nodes of each edge.
!
!    Output, integer ( kind = 4 ) CIRCUIT(NEDGE), the Euler circuit, as a series of nodes.
!
  implicit none

  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) circuit(nedge)
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
  integer ( kind = 4 ) nbreak
  integer ( kind = 4 ) ncopy
  integer ( kind = 4 ) numarc
  integer ( kind = 4 ) numnode
!
!  The number of times each node has been visited begins at 0.
!
  iwork3(1:nnode) = 0
  circuit(1:nedge) = 0
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
    circuit(1:numnode) = iwork1(1:numnode)
    return
  end if

  insert = 1
  ipivot = iwork6(insert)
  iforwd = 0

  do

    ncopy = 1
    ibase = iwork1(1)
    locbas = 1
    circuit(ncopy) = ibase
!
!  A path identified before is added to the circuit.
!
80  continue

    if ( ibase == ipivot ) then

      j = iwork4(insert) + iforwd
      k = iwork5(insert) + iforwd

      do l = j, k
        ncopy = ncopy + 1
        circuit(ncopy) = iwork1(l)
        iwork1(l) = 0
      end do

      ncopy = ncopy + 1
!
!  Add the intersecting node to the circuit.
!
      circuit(ncopy) = ibase
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
            circuit(ncopy) = ibase
          end if

        end do

      end if

    else

      ncopy = ncopy + 1

       if ( ncopy <= numnode ) then
         locbas = locbas + 1
         ibase = iwork1(locbas)
         circuit(ncopy) = ibase
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

    iwork1(1:nedge) = circuit(1:nedge)

  end do

  return
end
subroutine graph_arc_euler_circ_cand ( nedge, inode, jnode, circuit, k, &
  nstack, stack, maxstack, ncan, iwork )

!*****************************************************************************80
!
!! GRAPH_ARC_EULER_CIRC_CAND finds candidates for the K-th edge of an Euler circuit.
!
!  Discussion:
!
!    This routine is used in conjunction with I4VEC_BACKTRACK, which directs the 
!    search for a complete Euler circuit.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 July 2000
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
!    Input, integer ( kind = 4 ) NEDGE, the number of edges in the graph.
!
!    Input, integer ( kind = 4 ) INODE(NEDGE), JNODE(NEDGE), the edge array of the graph.
!    The I-th edge extends from node INODE(I) to JNODE(I).
!
!    Input, integer ( kind = 4 ) CIRCUIT(NEDGE), CIRCUIT(I) is the I-th edge in the circuit.  
!    A full circuit will have NEDGE edges, but on input we only have K-1.
!
!    Input, integer ( kind = 4 ) K, the index of the next edge to be determined in circuit.
!
!    Input/output, integer ( kind = 4 ) NSTACK, the current length of the stack.
!
!    Input, integer ( kind = 4 ) STACK(MAXSTACK).  As yet unused candidates for positions
!    1 to K-1.
!
!    Input, integer ( kind = 4 ) MAXSTACK, the dimension of STACK.
!
!    Input/output, integer ( kind = 4 ) NCAN(NEDGE), the number of candidates for each 
!    position.
!
!    Workspace, integer IWORK(NEDGE).
!
  implicit none

  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) maxstack

  integer ( kind = 4 ) circuit(nedge)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) inode(nedge)
  integer ( kind = 4 ) it
  integer ( kind = 4 ) iwork(nedge)
  integer ( kind = 4 ) jnode(nedge)
  integer ( kind = 4 ) k
  logical lwork(nedge)
  integer ( kind = 4 ) ncan(nedge)
  integer ( kind = 4 ) nstack
  integer ( kind = 4 ) stack(maxstack)

  ncan(k) = 0

  if ( k == 1 ) then
    iwork(1) = jnode(1)
    stack(1) = 1
    nstack = 1
    ncan(k) = 1
    return
  end if
 
  if ( 2 < k ) then
    iwork(k-1) = inode(circuit(k-1)) + jnode(circuit(k-1)) - iwork(k-2)
  end if
 
  it = iwork(k-1)
 
  do i = 1, nedge
    lwork(i) = it == inode(i) .or. it == jnode(i)
  end do
 
  lwork(circuit(1:k-1)) = .false.
 
  do i = 1, nedge
    if ( lwork(i) ) then
      if ( maxstack <= nstack ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'GRAPH_ARC_EULER_CIRC_CAND - Fatal error!'
        write ( *, '(a)' ) '  Stacksize exceeded!'
        stop
      end if
      nstack = nstack + 1
      stack(nstack) = i
      ncan(k) = ncan(k) + 1
    end if
  end do
 
  return
end
subroutine graph_arc_euler_circ_next ( nedge, inode, jnode, circuit, stack, &
  maxstack, ncan, more )

!*****************************************************************************80
!
!! GRAPH_ARC_EULER_CIRC_NEXT returns the next Euler circuit for a graph.
!
!  Discussion:
!
!    The routine produces all the Euler circuits of a graph, one at a time.
!
!    An Euler circuit of a graph is a path starting at some node, 
!    using all the edges of the graph exactly once, and returning
!    to the starting node.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 August 2000
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
!    Input, integer ( kind = 4 ) NEDGE, the number of edges in the graph.
!
!    Input, integer ( kind = 4 ) INODE(NEDGE), JNODE(NEDGE), the edge array of the graph.
!    The I-th edge extends from node INODE(I) to JNODE(I).
!
!    Output, integer ( kind = 4 ) CIRCUIT(NEDGE).  If MORE = TRUE on output, then CIRCUIT
!    contains the edges, in order, that constitute this circuit.
!
!    Workspace, integer STACK(MAXSTACK).  
!
!    Input, integer ( kind = 4 ) MAXSTACK, the dimension of STACK.
!
!    Workspace, integer NCAN(NEDGE), the number of candidates for each position.
!
!    Input/output, logical MORE.
!    On first call, set MORE to .FALSE, and do not alter it after.
!    On return, MORE is TRUE if another circuit has been returned in
!    IARRAY, and FALSE if there are no more circuits.
!
  implicit none

  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) maxstack

  integer ( kind = 4 ) circuit(nedge)
  integer ( kind = 4 ) inode(nedge)
  integer ( kind = 4 ), save :: indx = 0
  integer ( kind = 4 ) iwork(nedge)
  integer ( kind = 4 ) jnode(nedge)
  integer ( kind = 4 ), save :: k = 0
  logical more
  integer ( kind = 4 ) ncan(nedge)
  integer ( kind = 4 ), save :: nstack = 0
  integer ( kind = 4 ) stack(maxstack)

  if ( .not. more ) then
    indx = 0
    k = 0
    more = .true.
    nstack = 0
  end if
 
  do
 
    call i4vec_backtrack ( nedge, circuit, indx, k, nstack, stack, maxstack, &
      ncan )
 
    if ( indx == 1 ) then

      exit

    else if ( indx == 2 ) then

      call graph_arc_euler_circ_cand ( nedge, inode, jnode, circuit, k, &
        nstack, stack, maxstack, ncan, iwork )

    else

      more = .false.
      exit

    end if

  end do
 
  return
end
subroutine graph_arc_example_diamond ( inode, jnode, maxedge, nedge, nnode, &
  x, y, z )

!*****************************************************************************80
!
!! GRAPH_ARC_EXAMPLE_DIAMOND returns the graph of a "diamond" 3D shape.
!
!  Example:
!
!        1
!      /| |\
!     / | | \
!    2--3-4--5--(2)
!    |  | |  |
!    6--7-8--9--(6)
!     \ | | /
!      \| |/
!       10
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters
!
!    Output, integer ( kind = 4 ) INODE(MAXEDGE), JNODE(MAXEDGE), the NEDGE 
!    edges of the graph.  The I-th edge connects nodes INODE(I) and
!    JNODE(I).
!
!    Input, integer ( kind = 4 ) MAXEDGE, the maximum number of edges allocated
!    in the EDGE array.  MAXEDGE should be at least 20.
!
!    Output, integer ( kind = 4 ) NEDGE, the number of edges, which is 20.
!
!    Output, integer ( kind = 4 ) NNODE, the number of nodes, which is 10.
!
!    Output, real ( kind = 8 ) X(NNODE), Y(NNODE), Z(NNODE), the 
!    locations for the nodes.
!
  implicit none

  integer ( kind = 4 ) maxedge

  integer ( kind = 4 ) inode(maxedge)
  integer ( kind = 4 ) jnode(maxedge)
  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) nnode
  real ( kind = 8 ) x(10)
  real ( kind = 8 ) y(10)
  real ( kind = 8 ) z(10)

  nedge = 20
  nnode = 10

  if ( maxedge < nedge ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRAPH_ARC_EXAMPLE_DIAMOND - Fatal error!'
    write ( *, '(a,i8)' ) '  Increase MAXEDGE to at least ', nedge
    stop
  end if

  inode(1) = 1
  jnode(1) = 2

  inode(2) = 1
  jnode(2) = 3

  inode(3) = 1
  jnode(3) = 4

  inode(4) = 1
  jnode(4) = 5

  inode(5) = 2
  jnode(5) = 6

  inode(6) = 3
  jnode(6) = 7

  inode(7) = 4
  jnode(7) = 8

  inode(8) = 5
  jnode(8) = 9

  inode(9) = 6
  jnode(9) = 10

  inode(10) = 7
  jnode(10) = 10

  inode(11) = 8
  jnode(11) = 10

  inode(12) = 9
  jnode(12) = 10

  inode(13) = 2
  jnode(13) = 3

  inode(14) = 3
  jnode(14) = 4

  inode(15) = 4
  jnode(15) = 5

  inode(16) = 5
  jnode(16) = 2

  inode(17) = 6
  jnode(17) = 7

  inode(18) = 7
  jnode(18) = 8

  inode(19) = 8
  jnode(19) = 9

  inode(20) = 9
  jnode(20) = 6

  x(1) =  0.0D+00
  y(1) =  0.0D+00
  z(1) =  2.0D+00

  x(2) =  0.5D+00
  y(2) = -0.5D+00
  z(2) =  1.0D+00

  x(3) =  0.5D+00
  y(3) =  0.5D+00
  z(3) =  1.0D+00

  x(4) = -0.5D+00
  y(4) =  0.5D+00
  z(4) =  1.0D+00

  x(5) = -0.5D+00
  y(5) = -0.5D+00
  z(5) =  1.0D+00

  x(6) =  0.5D+00
  y(6) = -0.5D+00
  z(6) = -1.0D+00

  x(7) =  0.5D+00
  y(7) =  0.5D+00
  z(7) = -1.0D+00

  x(8) = -0.5D+00
  y(8) =  0.5D+00
  z(8) = -1.0D+00

  x(9) = -0.5D+00
  y(9) = -0.5D+00
  z(9) = -1.0D+00

  x(10) =  0.0D+00
  y(10) =  0.0D+00
  z(10) = -2.0D+00

  return
end
subroutine graph_arc_face ( face, face_count, face_order, iface, jface, &
  inode, jnode, maxface, maxorder, nedge, nface, nnode )

!*****************************************************************************80
!
!! GRAPH_ARC_FACE constructs a set of faces for a plane graph.
!
!  Discussion:
!
!    Warning: This is an experimental code.
!
!    The reason this routine was written was to handle the problem of
!    converting certain forms of 3D graphics data from a point and line
!    representation to a list of faces.  While at first glance, this
!    seemed an easy task, it turned out to be one of those problems
!    that becomes harder the longer it is considered.  Particularly
!    vexing was the idea that it might be possible to do this reconstruction
!    without using any of the geometric data supplied with the connectivity
!    data.
!
!    The guiding idea was that a face ought to be a "short" cycle of
!    the graph, and that every edge ought to appear in two separate faces.
!    The resulting method should work for a connected graph which is
!    planar, or merely orientable.  A planar graph will result from a
!    reasonable "triangulation" (meaning decomposition into arbitrary
!    polygons) of the surface of a 3D object that has no holes.
!
!    This algorithm will also handle the case where the graph is not planar,
!    but results from the triangulation of a more complicated 3D object,
!    such as one that includes holes.  Even a Klein bottle, which is
!    a manifold, but not orientable, can be handled, although it may not
!    be possible then to assign a consistent orientation to the faces.
!
!    By the way, this problem is MUCH easier if we can assume that all
!    the faces use the same number of edges, such as in a triangular
!    decomposition.  This algorithm makes no such assumption.
!
!    If the graph is planar, then the decomposition into
!    faces allows us to define the dual graph.  The dual graph H of the
!    planar graph G comprises:
!
!    * nodes V(I), each of which corresponds to a face F(I) of G;
!
!    * edges (V(I), V(J)).  V(I) and V(J) share an edge in H if and only
!      if the faces F(I) and F(J) share an edge in G.  (Thus G and H
!      have the same number of edges).
!
!    In the terminology of this routine, the dual graph has NFACE nodes,
!    NEDGE edges, and the corresponding edge arrays are simply IFACE and
!    JFACE.
!
!  Formula:
!
!    If the graph is actually planar, we can regard it as the flattened
!    triangulation of a connected solid shape, and so we can apply Euler's
!    formula:
!
!      Faces + Vertices = Edges + 2
!
!    This means that we can predict beforehand that the number of faces
!    produced by this routine will be
!
!      NFACE = NEDGE + 2 - NNODE.
!
!  Notes:
!
!    The faces produced by this routine may actually overlap.  Without
!    geometric data, this is surely a possibility, since a graph may
!    have more than one embedding.  For instance, consider the following
!    two embeddings of the same graph:
!
!      A-----B       A-----B
!      |     |       |     |
!      |  E  |       D-----C
!      | / \ |        \   /
!      |/   \|         \ /
!      D-----C          E
!
!    This routine will report the two faces (A,B,C,D) and (C,D,E),
!    although in the first embedding one face seems to be part of
!    another.  This is not as bad as it might seem, since the routine
!    prefers the "smaller" face (A,B,C,D) over (A,B,C,E,D).
!
!
!    A second problem is best illustrated with a simple example.
!    Suppose we have a thin triangular rod, and that we have triangulated
!    the surface of this rod, so that the cross section of the rod
!    is a triangular graph, and the sides are made up of, say, squares.
!    Then this routine will report all the "internal" triangles as
!    faces.  It will still find the "true" faces on the sides, but
!    since it is possible to go around the diameter of the object
!    in a very few steps, the algorithm produces faces we would not
!    expect.
!
!  Restrictions:
!
!    The algorithm will fail if the graph cannot be regarded, at least
!    locally, as the triangulation of a smooth surface.  Smoothness
!    problems will not occur if the graph is planar, or results from
!    the triangulation of a 3D object, which might include holes.
!
!    The graph should be connected.
!
!    There should be no nodes of degree 1.
!
!  Method:
!
!    We have no geometric data from which to deduce physical positions
!    of the nodes.  We are only given that the graph is planar, so that
!    there is at least one embedding of the graph in the plane.
!
!    Our data structure for the method will use arrays called IFACE and JFACE.
!    For each edge I, IFACE(I) and JFACE(I) will eventually hold the
!    indices of the two faces that the edge is part of.  We begin
!    the algorithm by setting all entries of IFACE and JFACE to 0.
!
!    The second step is to find one cycle in the graph, of the shortest
!    length possible.  This cycle constitutes our first face.  We update
!    the appropriate entries of IFACE or JFACE, marking each edge as having
!    been used once.
!
!    The third step is to add one more face to our collection of faces.
!    The new face will be adjacent to the current collection of faces,
!    but will include at least one completely unused edge, if possible.
!
!    To guarantee this, we consider every edge that is part of our
!    collection of faces, and that has only been used once.  We look
!    at the endpoints of each of these edges., and
!
!      We search for an adjacent edge that has not been used.
!      If we find such an edge, then the first two edges of our next face
!      are the edge that was already part of the set of faces, and the
!      unused edge.
!
!      If we cannot find such an edge, then we repeat the search, starting
!      with an edge in the face set that has been used once.  But now
!      when we search for adjacent edges, we will consider using one that
!      has already been used once.
!
!    We then search for a path that will return us to the initial
!    node of the first edge.  Using a breadth-first search, we expect
!    to find the shortest path back to that node, and we assume that
!    this represents a face.  Again, we update the IFACE and JFACE arrays.
!
!    We repeat the third step until there are no more edges in the
!    collection of faces that have not been used twice.  Assuming the
!    graph is connected, this means that every face has been found.
!
!  Improvements:
!
!    It shouldn't be hard to modify the code to handle graphs that are
!    not connected.
!
!    If the edge arrays INODE and JNODE were sorted and indexed, some
!    operations could be done more efficiently.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) FACE(MAXORDER,MAXFACE), contains the list of edges
!    which make up each face.  Face I is made up of the edges
!    FACE(1,I) through FACE(FACE_ORDER(I),I).
!
!    Output, integer ( kind = 4 ) FACE_COUNT(NEDGE).  For each edge I, FACE_COUNT(I)
!    is the number of faces to which the edge belongs.  This value should
!    be 0, 1 or 2.
!
!    Output, integer ( kind = 4 ) IFACE(NEDGE), JFACE(NEDGE).  IFACE(I) and JFACE(I)
!    are the two faces to which edge I belongs.  Either or both may be zero
!    if the algorithm fails.
!
!    Input, integer ( kind = 4 ) INODE(NEDGE), JNODE(NEDGE), the edge list for the graph.
!    The I-th edge connects nodes INODE(I) and JNODE(I).
!
!    Input, integer ( kind = 4 ) MAXFACE, the maximum number of faces for which storage
!    has been set aside in FACE and FACE_ORDER.
!
!    Input, integer ( kind = 4 ) MAXORDER, the maximum number of edges for which storage
!    has been set aside in FACE.
!
!    Input, integer ( kind = 4 ) NEDGE, the number of edges.
!
!    Output, integer ( kind = 4 ) NFACE, the number of faces found by the algorithm.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) FACE_ORDER(MAXFACE).  The number of edges used
!    in constructing face I is stored in FACE_ORDER(I).
!
  implicit none

  logical, parameter :: debug = .false.

  integer ( kind = 4 ) maxface
  integer ( kind = 4 ) maxorder
  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) face(maxorder,maxface)
  integer ( kind = 4 ) face_count(nedge)
  integer ( kind = 4 ) face_order(maxface)
  integer ( kind = 4 ) faceval
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iedge
  integer ( kind = 4 ) iface(nedge)
  integer ( kind = 4 ) inode(nedge)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jface(nedge)
  integer ( kind = 4 ) jnode(nedge)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) length
  integer ( kind = 4 ) nface
  integer ( kind = 4 ) nface_old
  integer ( kind = 4 ) nodes(3)
  integer ( kind = 4 ) nstart
!
!  Initialization.  No arc belongs to any face.
!
  if ( debug ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRAPH_ARC_FACE - Debug:'
    write ( *, '(a)' ) '  Initialization'
  end if

  nface = 0
  face_count(1:nedge) = 0
  iface(1:nedge) = 0
  jface(1:nedge) = 0
  face_order(1:maxface) = 0
  face(1:maxorder,1:maxface) = 0
!
!  We start here.  We may also jump back here if we have used up all the
!  connected parts of a graph, and need to jump to a new piece.
!
!  Find one new face of minimal length.
!
5 continue

  nface_old = nface

  do length = 3, nnode

    do iedge = 1, nedge

      nodes(1) = inode(iedge)
      nodes(2) = jnode(iedge)
      nstart = 2

      call graph_arc_face_next ( face, face_count, face_order, iface, jface, &  
        inode, jnode, maxface, maxorder, nedge, nface, nnode, nodes, nstart )

      if ( nface_old < nface ) then
        go to 10
      end if

    end do

  end do

  if ( nface == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRAPH_ARC_FACE - Note.'
    write ( *, '(a)' ) '  Could not find any starting face.'
  end if

  go to 60
!
!  Find an edge that is in one face, but not two.
!
10    continue

  if ( debug ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRAPH_ARC_FACE - Debug:'
    write ( *, '(a,i8)' ) '  Found starting face #:', nface
    write ( *, '(a,i8)' ) '  Order is ', face_order(nface)
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Vertices:'
    write ( *, '(a)' ) ' '
    do i = 1, face_order(nface)
      write ( *, '(6i8)' ) face(i,nface)
    end do
  end if

  iedge = 0
!
!  Look for an edge with FACE_COUNT of 1.
!
20    continue

  iedge = iedge + 1

  if ( face_count(iedge) == 1 ) then
    go to 30
  else if ( iedge < nedge ) then
    go to 20
  else
    if ( debug ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'GRAPH_ARC_FACE - Debug:'
      write ( *, '(a)' ) '  No more nearby edges to try.'
    end if
!
!  Here, I'd like to be able to jump back and scrounge for other
!  islands of edges, but something's not right.
!
!       go to 5

    go to 60

  end if
!
!  The face will start with the two nodes of edge IEDGE.
!
30    continue

  if ( debug ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRAPH_ARC_FACE - Debug:'
    write ( *, '(a)' ) '  Found a starting edge:'
    write ( *, '(i8)' ) inode(iedge)
    write ( *, '(i8)' ) jnode(iedge)
  end if

  nodes(1) = inode(iedge)
  nodes(2) = jnode(iedge)
!
!  Look for an edge incident to JNODE(IEDGE).  This new edge should have
!  been used just FACEVAL times already.  (FACEVAL is preferably 0, but if
!  we can't find any at that rate, we'll try FACEVAL = 1).
!
  faceval = 0

40    continue

  do i = 1, nedge

    if ( face_count(i) == faceval ) then

      if ( inode(i) == nodes(2) .and. jnode(i) /= nodes(1) ) then
        nodes(3) = jnode(i)
        go to 50
      else if ( jnode(i) == nodes(2) .and. inode(i) /= nodes(1) ) then
        nodes(3) = inode(i)
        go to 50
      end if

    end if

  end do
!
!  If we "fell through" with FACEVAL = 0, then try the search again
!  with FACEVAL = 1.
!
  if ( faceval == 0 ) then

    faceval = 1
    go to 40
!
!  If we fell through with FACEVAL = 1, then we couldn't find any
!  way to use this edge.  Mark it as though it were used, and move on.
!
  else

    if ( debug ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'GRAPH_ARC_FACE - Debug:'
      write ( *, '(a)' ) '  Failure.'
      write ( *, '(a,i8)' ) '  Cannot hook up to edge IEDGE:', iedge
      write ( *, '(2i8)' ) nodes(1), nodes(2)
    end if

    face_count(iedge) = 2
    go to 20

  end if
!
!  Now call FACENEXT to search for the shortest cycle that begins
!  NODES(1), NODES(2), NODES(3), and which involves only edges that
!  have been used once or less.
!
50    continue

  nface_old = nface
  nstart = 3

  call graph_arc_face_next ( face, face_count, face_order, iface, jface, &
    inode, jnode, maxface, maxorder, nedge, nface, nnode, nodes, nstart )

  if ( nface_old < nface ) then

    if ( debug ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'GRAPH_ARC_FACE - Debug'
      write ( *, '(a,i8)' ) '  NFACE_OLD = ', nface_old
      write ( *, '(a,i8)' ) '  NFACE = ', nface
      write ( *, '(a,i8)' ) '  Order is ', face_order(nface)
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Vertices:'
      write ( *, '(a)' ) ' '
      do i = 1, face_order(nface)
        write ( *, '(6i8)' ) face(i,nface)
      end do
      write ( *, '(a)' ) '  Trying the big loop again.'
    end if

    go to 10

  end if
!
!  The algorithm has failed.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GRAPH_ARC_FACE - Error!'
  write ( *, '(a)' ) '  The algorithm has failed.'
  write ( *, '(a)' ) '  Only some of the faces were found.'
!
!  Cleanup
!
60    continue

  do i = 1, nface
    face_order(i) = min ( face_order(i), maxorder )
  end do

  do i = 1, nface
    do j = 1, face_order(i)
      k = face(j,i)
      if ( k < 0 ) then
        face(j,i) = jnode(-k)
      else
        face(j,i) = inode(k)
      end if
    end do
  end do

  return
end
subroutine graph_arc_face_next ( face, face_count, face_order, iface, jface, &
  inode, jnode, maxface, maxorder, nedge, nface, nnode, nodes, nstart )

!*****************************************************************************80
!
!! GRAPH_ARC_FACE_NEXT tries to complete the next face, given a few starting nod
!
!  Discussion:
!
!    This is a utility routine, called by GRAPH_ARC_FACE, which
!    constructs all the faces of a graph.  GRAPH_ARC_FACE finds the first
!    two or three nodes of a face, and then calls this routine, which
!    attempts to complete the face by using a breadth-first search
!    from the final given node of the face.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) FACE(MAXORDER,MAXFACE), contains the list of edges
!    which make up each face.  Face I is made up of the edges
!    FACE(1,I) through FACE(FACE_ORDER(I),I).  If a new face is found, this
!    array is updated.
!
!    Input/output, integer ( kind = 4 ) FACE_COUNT(NEDGE).  For each edge I, FACE_COUNT(I)
!    is the number of faces to which the edge belongs.  This value will
!    be 0, 1 or 2.  If a new face is found, this data is updated.
!
!    Input/output, integer ( kind = 4 ) FACE_ORDER(MAXFACE).  The number of edges used
!    in constructing face I is stored in FACE_ORDER(I).
!
!    Input/output, integer ( kind = 4 ) IFACE(NEDGE), JFACE(NEDGE).  IFACE(I) and JFACE(I)
!    are the two faces to which edge I belongs.  Either or both may be zero
!    if the algorithm fails.  If a new face is found, this data is updated.
!
!    Input, integer ( kind = 4 ) INODE(NEDGE), JNODE(NEDGE), the edge list for the graph.
!    The I-th edge connects nodes INODE(I) and JNODE(I).
!
!    Input, integer ( kind = 4 ) MAXFACE, the maximum number of faces for which storage
!    has been set aside in FACE and FACE_ORDER.
!
!    Input, integer ( kind = 4 ) MAXORDER, the maximum number of edges for which storage
!    has been set aside in FACE.
!
!    Input, integer ( kind = 4 ) NEDGE, the number of edges.
!
!    Input/output, integer ( kind = 4 ) NFACE.  NFACE is the number of faces found so far.
!    This value will be updated by this routine if a new face is found.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) NODES(NSTART), the first few nodes in the partial cycle.
!
!    Input, integer ( kind = 4 ) NSTART, the number of nodes in the partial cycle.
!
!  Workspace:
!
!    Workspace, integer DAD(NNODE), used during the breadth first search
!    of the graph, to point backwards from each node to its predecessor
!    in a path.
!
!    Workspace, integer INDEX(NNODE), used during the breadth first search
!    to label nodes that have been visited.
!
  implicit none

  integer ( kind = 4 ) maxface
  integer ( kind = 4 ) maxorder
  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) nnode
  integer ( kind = 4 ) nstart

  integer ( kind = 4 ) dad(nnode)
  integer ( kind = 4 ) face(maxorder,maxface)
  integer ( kind = 4 ) face_count(nedge)
  integer ( kind = 4 ) face_order(maxface)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iedge2
  integer ( kind = 4 ) iface(nedge)
  integer ( kind = 4 ) index(nnode)
  integer ( kind = 4 ) inode(nedge)
  integer ( kind = 4 ) istart1
  integer ( kind = 4 ) istart2
  integer ( kind = 4 ) itemp
  integer ( kind = 4 ) jface(nedge)
  integer ( kind = 4 ) jnode(nedge)
  integer ( kind = 4 ) kedge
  integer ( kind = 4 ) kk
  integer ( kind = 4 ) nadd
  integer ( kind = 4 ) nface
  integer ( kind = 4 ) nodei
  integer ( kind = 4 ) nodej
  integer ( kind = 4 ) npath
  integer ( kind = 4 ) nodes(nstart)
!
!  Initialization.
!
  index(1:nnode) = 0
  dad(1:nnode) = 0

  istart1 = nodes(1)
  istart2 = nodes(2)

  do i = 1, nstart

    npath = i

    if ( i == 1 ) then
      dad(nodes(i)) = -1
    else
      dad(nodes(i)) = nodes(i-1)
    end if

    index(nodes(i)) = i

  end do
!
!  From the nodes with INDEX = NPATH, consider all neighbors.
!
10    continue

  npath = npath + 1
  nadd = 0

  do iedge2 = 1, nedge

    if ( index(inode(iedge2)) == npath-1 .and. index(jnode(iedge2)) == 0 ) then

      nodei = inode(iedge2)
      nodej = jnode(iedge2)

    else if ( index(jnode(iedge2)) == npath-1 .and. &
      index(inode(iedge2)) == 0 ) then

      nodei = jnode(iedge2)
      nodej = inode(iedge2)

    else if ( index(inode(iedge2)) == npath-1 .and. &
      jnode(iedge2) == istart1 ) then

      nodei = inode(iedge2)
      nodej = jnode(iedge2)

    else if ( index(jnode(iedge2)) == npath-1 .and. &
      inode(iedge2) == istart1 ) then

      nodei = jnode(iedge2)
      nodej = inode(iedge2)

    else

      nodei = 0
      nodej = 0

    end if

    if ( nodei /= 0 .and. nodej /= istart1 ) then

      nadd = nadd + 1
      index(nodej) = npath
      dad(nodej) = nodei
!
!  Success if the marked node is the starting point (except when
!  using the edge (ISTART2,ISTART1)).
!
    else if ( nodej == istart1 .and. nodei == istart2 ) then

    else if ( nodej == istart1 .and. nodei /= istart2 ) then

      nface = nface + 1

20    continue
!
!  Find the edge KK which joins NODEI and NODEJ.
!
      do kk = 1, nedge

        if ( ( inode(kk) == nodei .and. jnode(kk) == nodej ) .or. &
             ( jnode(kk) == nodei .and. inode(kk) == nodej ) ) then

          face_count(kk) = face_count(kk) + 1
          itemp = face_count(kk)

          if ( itemp == 1 ) then
            iface(kk) = nface
          else if ( itemp == 2 ) then
            jface(kk) = nface
          end if

          if ( inode(kk) == nodei ) then
            kedge = kk
          else
            kedge = -kk
          end if

          exit

        end if

      end do

      nodej = nodei
!
!  Add the edge to the face-to-edge list.
!
      if ( nface <= maxface ) then

        if ( face_order(nface) < maxorder ) then
          face_order(nface) = face_order(nface) + 1
        end if

        if ( face_order(nface) <= maxorder ) then
          face(face_order(nface),nface) = kedge
        end if

      end if

      if ( nodej /= istart1 ) then
        nodei = dad(nodej)
        go to 20
      end if

      return

    end if

  end do
!
!  If we were able to proceed another step, and we haven't exceeded
!  our limit, then go back and take another step.
!
  if ( 0 < nadd .and. npath <= nnode ) then
    go to 10
  end if

  return
end
subroutine graph_arc_is_eulerian ( nnode, nedge, inode, jnode, degree, result )

!*****************************************************************************80
!
!! GRAPH_ARC_IS_EULERIAN determines if a graph is Eulerian from its edge list.
!
!  Definition:
!
!    A graph is Eulerian if there exists a circuit through the graph
!    which uses every edge once.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 September 1999
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
!    Input, integer ( kind = 4 ) INODE(NEDGE), JNODE(NEDGE), the pairs of nodes
!    that form the edges.
!
!    Output, integer ( kind = 4 ) DEGREE(NNODE), the degree of each node, that is,
!    the number of edges that include the node.
!
!    Output, integer ( kind = 4 ) RESULT.
!    0, the graph is not Eulerian.
!    1, the graph is Eulerian, but the starting and ending nodes are different.
!    2, the graph is Eulerian, and there is a closed Euler circuit.
!
  implicit none

  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) degree(nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) inode(nedge)
  integer ( kind = 4 ) jnode(nedge)
  integer ( kind = 4 ) nodd
  integer ( kind = 4 ) result

  call graph_arc_degree ( nnode, nedge, inode, jnode, degree )

  nodd = 0

  do i = 1, nnode

    if ( mod ( degree(i), 2 ) == 1 ) then
      nodd = nodd + 1
    end if

  end do

  if ( nodd == 0 ) then
    result = 2
  else if ( nodd == 2 ) then
    result = 1
  else
    result = 0
  end if

  return
end
subroutine graph_arc_match ( nnode, nedge, inode, jnode, type, match )

!*****************************************************************************80
!
!! GRAPH_ARC_MATCH finds a maximum matching in a bipartite graph.
!
!  Discussion:
!
!    The nodes of the graph are assumed to be divided into two groups,
!    and it is desired to determine as matching that is as large as possible.
!    The matching is a set of pairs ( NODE(I), NODE(J) ) with the properties:
!
!    * NODE(I) is in group 1 and NODE(J) is in group 2;
!    * there is an edge between NODE(I) and NODE(J);
!    * NODE(I) and NODE(J) are not used in any other pairing in the matching.
!
!    The user inputs the edge list that defines the graph, and a set of
!    labels that classify the nodes as being in one group or the other.
!    It is not necessary that the graph actually be bipartite; edges between
!    nodes in the same group are allowed, but they will not affect the
!    outcome in any way.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 July 2000
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
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) NEDGE, the number of edges.
!
!    Input, integer ( kind = 4 ) INODE(NEDGE), JNODE(NEDGE), the end nodes of the edges.
!
!    Input, integer ( kind = 4 ) TYPE(NNODE), labels the two types of nodes in the graph.
!    Normally, TYPE(I) would be 0 or 1, but any two distinct values will do.
!
!    Output, integer ( kind = 4 ) MATCH(NNODE), the matching node for each node, or 0
!    if no match was assigned.
!
  implicit none

  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) capflo(2,2*nedge+2*nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icut(nnode+2)
  integer ( kind = 4 ) iendpt(2,2*nedge+2*nnode)
  integer ( kind = 4 ) in
  integer ( kind = 4 ) inode(nedge)
  integer ( kind = 4 ) isink
  integer ( kind = 4 ) isorce
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jn
  integer ( kind = 4 ) jnode(nedge)
  integer ( kind = 4 ) match(nnode)
  integer ( kind = 4 ) nedge2
  integer ( kind = 4 ) nnode2
  integer ( kind = 4 ) node_flow(nnode+2)
  integer ( kind = 4 ) type(nnode)

  match(1:nnode) = 0
!
!  Create a network from the graph, with two extra nodes.
!
  isorce = nnode + 1
  isink = nnode + 2
  nnode2 = nnode + 2

  j = 0

  do i = 1, nedge

    in = inode(i)
    jn = jnode(i)

    if ( type(in) /= type(jn) ) then

      j = j + 1
      iendpt(1,j) = inode(i)
      iendpt(2,j) = jnode(i)
      capflo(1,j) = 1
      capflo(2,j) = 0

      j = j + 1
      iendpt(1,j) = jnode(i)
      iendpt(2,j) = inode(i)
      capflo(1,j) = 1
      capflo(2,j) = 0

    end if

  end do
!
!  Nodes of type 1 are connected to the source, 
!  and nodes of type 2 are connected to the sink.
!
  do i = 1, nnode

    if ( type(i) == type(1) ) then
      j = j + 1
      iendpt(1,j) = isorce
      iendpt(2,j) = i
      capflo(1,j) = 1
      capflo(2,j) = 0
      j = j + 1
      iendpt(1,j) = i
      iendpt(2,j) = isorce
      capflo(1,j) = 1
      capflo(2,j) = 0
    else
      j = j + 1
      iendpt(1,j) = i
      iendpt(2,j) = isink
      capflo(1,j) = 1
      capflo(2,j) = 0
      j = j + 1
      iendpt(1,j) = isink
      iendpt(2,j) = i
      capflo(1,j) = 1
      capflo(2,j) = 0
    end if
  end do
!
!  Determine the maximum flow on the network.
!
!  Then a pair of nodes connected by an edge that has a network flow of 1
!  are part of the maximal matching.
!
  nedge2 = j

  call network_flow_max ( nnode2, nedge2, iendpt, capflo, isorce, isink, &
    icut, node_flow )

  do i = 1, nedge2
    if ( iendpt(1,i) <= nnode .and. &
         iendpt(2,i) <= nnode .and. &
         0 < capflo(1,i) .and. &
         capflo(2,i) == 1 ) then
      in = iendpt(1,i)
      jn = iendpt(2,i)
      match(in) = jn
      match(jn) = in
    end if
  end do

  return
end
subroutine graph_arc_min_path ( nnode, nedge, inode, jnode, arcost, &
  istart, last, num_path, ispath, xlen )

!*****************************************************************************80
!
!! GRAPH_ARC_MIN_PATH finds the shortest path between two nodes.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 September 1999
!
!  Reference:
!
!    Hang Tong Lau,
!    Combinatorial Heuristic Algorithms in FORTRAN,
!    Springer Verlag, 1986.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes in the graph.
!
!    Input, integer ( kind = 4 ) NEDGE, the number of edges in the graph.
!
!    Input, integer ( kind = 4 ) INODE(NEDGE), JNODE(NEDGE), the edges of the graph,
!    describe by pairs of nodes.
!
!    Input, real ( kind = 8 ) ARCOST(NEDGE), the distance or cost of each edge.
!
!    Input, integer ( kind = 4 ) ISTART, LAST, are the two nodes between which a
!    shortest path is desired.
!
!    Output, integer ( kind = 4 ) NUM_PATH, the number of nodes in the shortest path.
!    NUM_PATH is zero if no path could be found.
!
!    Output, integer ( kind = 4 ) ISPATH(NNODE), lists the nodes in the shortest path,
!    from ISPATH(1) to ISPATH(NUM_PATH).
!
!    Output, real ( kind = 8 ) XLEN, the length of the shortest path 
!    from ISTART to LAST.
!
  implicit none

  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) nnode

  real ( kind = 8 ) arcost(nedge)
  real ( kind = 8 ) d
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) ient
  logical ifin
  integer ( kind = 4 ) ij
  integer ( kind = 4 ) inode(nedge)
  integer ( kind = 4 ) ispath(nnode)
  integer ( kind = 4 ) istart
  logical iwork1(nnode)
  integer ( kind = 4 ) iwork2(nnode)
  integer ( kind = 4 ) iwork3(nedge)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jnode(nedge)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) last
  integer ( kind = 4 ) num_path
  real ( kind = 8 ) wk4(nnode)
  real ( kind = 8 ) xlen

  wk4(1:nnode) = huge ( wk4(1) )
  iwork1(1:nnode) = .true.
  iwork2(1:nnode) = 0

  wk4(istart) = 0.0D+00
  i = istart
  iwork1(istart) = .false.
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
        d = wk4(i) + arcost(k)
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

    if ( ij == istart ) then
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
!    The input graph is represented by a list of edges.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 July 2000
!
!  Reference:
!
!    Hang Tong Lau,
!    Combinatorial Heuristic Algorithms in FORTRAN,
!    Springer Verlag, 1986.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes in the graph.
!
!    Input, integer ( kind = 4 ) NEDGE, the number of edges in the graph.
!
!    Input, integer ( kind = 4 ) INODE(NEDGE), JNODE(NEDGE), the start and end nodes
!    of the edges.
!
!    Input, real ( kind = 8 ) COST(NEDGE), the cost or length of each edge.
!
!    Output, integer ( kind = 4 ) ITREE(NNODE-1), JTREE(NNODE-1), the pairs of nodes that 
!    form the edges of the spanning tree.
!
!    Output, real ( kind = 8 ) TREE_COST, the total cost or length 
!    of the spanning tree.
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

  wk6(1:nnode) = huge ( wk6(1) )
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

    wk6(1:nnode) = huge ( wk6(1) )

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
subroutine graph_arc_ncolor_print ( nedge, inode, jnode, nnode, color, title )

!*****************************************************************************80
!
!! GRAPH_ARC_NCOLOR_PRINT prints out a node-colored graph from an edge list.
!
!  Discussion:
!
!    The printout is arranged to emphasize the colors of the neighboring nodes.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
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
subroutine graph_arc_node_count ( nedge, inode, jnode, mnode, nnode )

!*****************************************************************************80
!
!! GRAPH_ARC_NODE_COUNT counts the number of nodes in a graph.
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

  call i4vec_uniq ( 2*nedge, knode, nnode )

  return
end
subroutine graph_arc_print ( nedge, inode, jnode, title )

!*****************************************************************************80
!
!! GRAPH_ARC_PRINT prints out a graph from an edge list.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
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
subroutine graph_arc_to_ps ( file_name, inode, jnode, nedge, nnode, x, y )

!*****************************************************************************80
!
!! GRAPH_ARC_TO_PS writes graph information to a PostScript file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the output file.
!
!    Input, integer ( kind = 4 ) INODE(NEDGE), JNODE(NEDGE), the edge array.
!    The I-th edge connects nodes INODE(I) and JNODE(I).
!
!    Input, integer ( kind = 4 ) NEDGE, the number of edges.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, real ( kind = 8 ) X(NNODE), Y(NNODE), the X and Y components
!    of points.
!
  implicit none

  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) nnode

  real ( kind = 8 ) alpha
  real ( kind = 8 ) blue
  character ( len = 8 ) date
  character ( len = * ) file_name
  real ( kind = 8 ) green
  integer ( kind = 4 ) inode(nedge)
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) jnode(nedge)
  integer ( kind = 4 ) margin
  integer ( kind = 4 ) pagexmax
  integer ( kind = 4 ) pagexmin
  integer ( kind = 4 ) pageymax
  integer ( kind = 4 ) pageymin
  integer ( kind = 4 ) plotxmax
  integer ( kind = 4 ) plotxmin
  integer ( kind = 4 ) plotxmin2
  integer ( kind = 4 ) plotymax
  integer ( kind = 4 ) plotymin
  integer ( kind = 4 ) plotymin2
  real ( kind = 8 ) red
  real ( kind = 8 ) x(nnode)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) y(nnode)
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin
!
!  Bounding box.
!
  xmin = minval ( x(1:nnode) )
  xmax = maxval ( x(1:nnode) )
  ymin = minval ( y(1:nnode) )
  ymax = maxval ( y(1:nnode) )

  if ( xmin == xmax ) then
    xmin = x(1) - 0.5D+00
    xmax = x(1) + 0.5D+00
  end if

  if ( ymin == ymax ) then
    ymin = y(1) - 0.5D+00
    ymax = y(1) + 0.5D+00
  end if
!
!  Compute the scale factor.
!
  pagexmax = 612
  pagexmin = 0
  pageymax = 792
  pageymin = 0

  margin = 36

  plotxmax = pagexmax - margin
  plotxmin = pagexmin + margin
  plotymax = pageymax - margin
  plotymin = pageymin + margin

  alpha = min ( ( plotxmax - plotxmin ) / ( xmax - xmin ), &
                ( plotymax - plotymin ) / ( ymax - ymin ) )
!
!  Adjust PLOTXMIN and PLOTYMIN to center the image.
!
  plotxmin2 = 0.5D+00 * ( plotxmin + plotxmax - alpha * ( xmax - xmin ) )

  plotymin2 = 0.5D+00 * ( plotymin + plotymax - alpha * ( ymax - ymin ) )

  call get_unit ( iunit )

  open ( unit = iunit, file = file_name, status = 'replace', iostat = ios )

  if ( ios /= 0 ) then
    return
  end if
!
!  Write the prolog.
!
  write ( iunit, '(a)' ) '%!PS-Adobe-3.0'
  write ( iunit, '(a)' ) '%%Document-Fonts: Times-Roman'
  write ( iunit, '(a,a)' ) '%%Title: ' , trim ( file_name )
  write ( iunit, '(a)' ) '%%Creator: GRAFPACK(graph_arc_to_ps)'
  call date_and_time ( date )
  write ( iunit, '(a)' ) '%%CreationDate: ' // trim ( date )
  write ( iunit, '(a,4i5)' ) '%%BoundingBox', plotxmin, plotymin, plotxmax, &
    plotymax
  write ( iunit, '(a)' ) '%%LanguageLevel: 2'
  write ( iunit, '(a)' ) '%%EndComments'
  write ( iunit, '(a)' ) '%%BeginProlog'
  write ( iunit, '(a)' ) '%%EndProlog'
!
!  Set the line color.
!
  red = 0.0D+00
  green = 0.0D+00
  blue = 0.0D+00

  write ( iunit, '(3f7.4,a)' ) red, green, blue, ' setrgbcolor'
!
!  Draw lines.
!
  call edges_to_ps ( plotxmin2, plotymin2, alpha, iunit, inode, jnode, &
    nedge, nnode, x, y, xmin, ymin )
!
!  Set the fill color.
!
  red = 0.1
  green = 0.1
  blue = 0.7

  write ( iunit, '(3f7.4,a)' ) red, green, blue, ' setrgbcolor'
!
!  Draw points.
!
  call nodes_to_ps ( plotxmin2, plotymin2, alpha, iunit, nnode, x, y, &
    xmin, ymin )

  write ( iunit, '(a)' ) 'showpage'
!
!  Write the epilog.
!
  write ( iunit, '(a)' ) 'grestore'
  write ( iunit, '(a)' ) '%%Trailer'
  write ( iunit, '(a,i2)' ) '%%Pages: 1'

  close ( iunit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GRAPH_ARC_TO_PS'
  write ( *, '(a)' ) '  The data was written to the file: ' &
    // trim ( file_name )

  return
end
subroutine graph_arc_span_forest ( nnode, nedge, inode, jnode, ncomp, &
  component )

!*****************************************************************************80
!
!! GRAPH_ARC_SPAN_FOREST determines a graph's connectivity and spanning forest.
!
!  Discussion:
!
!    A (connected) component of a graph is a maximal subgraph which
!    is connected.
!
!    A tree is a connected graph containing no cycles.
!
!    A spanning tree of a connected graph is a subgraph which is a 
!    maximal tree.
!
!    A forest is a collection of trees, no two of which share a node.
!
!    A spanning forest of a possibly unconnected graph is a collection
!    containing a single spanning tree for each component of the graph.
!
!    The input graph may be connected or unconnected.
!
!    If the input graph is connected, this routine simply returns a
!    spanning tree for the graph.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 October 1999
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
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) NEDGE, the number of edges in the graph.
!
!    Input/output, integer ( kind = 4 ) INODE(NEDGE), JNODE(NEDGE), the edge list
!    of the graph.  On output, this array has been rearranged.  Edges
!    belonging to the spanning tree of component 1 are first, followed
!    by edges belonging to the other spanning trees, followed last by
!    edges that were not used in any spanning tree.
!
!    Output, integer ( kind = 4 ) NCOMP, the number of connected components of the graph.
!
!    Input/output, integer ( kind = 4 ) IENDPT(2,NEDGE), the edge array of the graph.  
!    IENDPT(1,I) and IENDPT(2,I) are the two nodes that make up edge I.
!  
!    On input, IENDPT describes the graph.
!
!    On output, the input entries of IENDPT have been reordered, so that
!    edges belonging to the spanning forest come first, followed by those
!    edges which are not part of the spanning forest.
!
!    Output, integer ( kind = 4 ) NCOMP, the number of connected components of the graph.
!
!    Output, integer ( kind = 4 ) IARRAY(NNODE).  IARRAY(I) is the component to which
!    node I belongs, and will take on values between 1 and NCOMP.
!
  implicit none

  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) component(nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) inode(nedge)
  integer ( kind = 4 ) inode2(nedge)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jnode(nedge)
  integer ( kind = 4 ) jnode2(nedge)
  integer ( kind = 4 ) left
  integer ( kind = 4 ) ncomp
  integer ( kind = 4 ) nstack
  integer ( kind = 4 ) num
  integer ( kind = 4 ) prev
  integer ( kind = 4 ) r
  integer ( kind = 4 ) right
  integer ( kind = 4 ) stack_node(nnode)
  integer ( kind = 4 ) stack_prev(nnode)
  integer ( kind = 4 ) v
  integer ( kind = 4 ) x_num(nnode)

  left = 0
  right = nedge + 1
  inode2(1:nedge) = 0
  jnode2(1:nedge) = 0
!
!  Part A:
!
  component(1:nnode) = 0
  x_num(1:nnode) = 0
  ncomp = 0
  v = 1
  num = 0

  nstack = 0
!
!  Part B:
!  Scan next V.
!
10    continue

  if ( nnode < v ) then
    inode(1:nedge) = inode2(1:nedge)
    jnode(1:nedge) = jnode2(1:nedge)
    return
  end if

  if ( component(v) /= 0 ) then
    v = v + 1
    go to 10
  end if
!
!  Begin the NCOMP-th component at V.
!
  ncomp = ncomp + 1
  num = num + 1

  component(v) = ncomp
  x_num(v) = num

  nstack = nstack + 1
  stack_node(nstack) = v
  stack_prev(nstack) = 0
!
!  Part C:
!  Is component NCOMP finished?
!
  do

    if ( nstack <= 0 ) then
      v = v + 1
      go to 10
    end if

    j = stack_node(nstack)
    prev = stack_prev(nstack)
    nstack = nstack - 1
!
!  Examine each vertex R that is adjacent to node J.
!
    do i = 1, nedge

      if ( inode(i) == j ) then
        r = jnode(i)
      else if ( jnode(i) == j ) then
        r = inode(i)
      else
        r = 0
      end if

      if ( r /= 0 ) then

        if ( component(r) == 0 ) then

          num = num + 1
          component(r) = ncomp
          x_num(r) = num

          nstack = nstack + 1
          stack_node(nstack) = r
          stack_prev(nstack) = j

          left = left + 1
          inode2(left) = j
          jnode2(left) = r

        else

          if ( r == prev .or. x_num(j) < x_num(r) ) then

          else

            right = right - 1
            inode2(right) = j
            jnode2(right) = r

          end if

        end if

      end if

    end do

  end do

  return
end
subroutine graph_arc_span_tree ( nedge, inode, jnode, nnode, dad )

!*****************************************************************************80
!
!! GRAPH_ARC_SPAN_TREE constructs the spanning tree of a graph.
!
!  Discussion:
!
!    If the graph is connected, then exactly one node will have no
!    parent, and a DAD value of -1.
!
!    If the graph is not connected, but divided into NCOMP components, then 
!    NCOMP nodes will have a DAD value of -1.
!
!    If the graph is connected, then once the tree is computed, the 
!    addition to the tree of any edge not included in the tree will 
!    form a cycle.  Since there are NNODE-1 edges in the tree, this will 
!    normally mean that there are NEDGE-(NNODE-1) "fundamental" cycles 
!    that can be generated in this way.
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
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEDGE, the number of edges.
!
!    Input, integer ( kind = 4 ) INODE(NEDGE), JNODE(NEDGE), the edge array of the graph. 
!    The I-th edge joins nodes INODE(I) and JNODE(I).
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) DAD(NNODE), the "father" array.  If node I is the
!    root of the tree spanning a given component of the graph, then 
!    DAD(I) = -1.  Otherwise, DAD(I) is the index of another node J in
!    the same component, such that the edge (I,J) is the first step
!    in the path along the tree from node I to the root of its component.
! 
  implicit none

  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) dad(nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iedge
  integer ( kind = 4 ) inode(nedge)
  integer ( kind = 4 ) jnode(nedge)
  integer ( kind = 4 ) nodei
  integer ( kind = 4 ) nodej
  integer ( kind = 4 ) nstacki
  integer ( kind = 4 ) nstackj
  integer ( kind = 4 ) stack(nnode)
!
!  Initialize.
!
  nstacki = 0
  nstackj = 0

  dad(1:nnode) = 0
  stack(1:nnode) = 0
!
!  Start at an unvisited node.
!
  do

    i = 0

    do

      i = i + 1

      if ( nnode < i ) then
        return
      end if

      if ( dad(i) == 0 ) then
        exit
      end if

    end do

    nodei = i
    dad(nodei) = - 1

    nstacki = 1
    stack(nstacki) = nodei
!
!  Search for unvisited neighbors of the last set of nodes.
!
    do

      do i = 1, nstacki

        nodei = stack(i)

        do iedge = 1, nedge

          if ( inode(iedge) == nodei ) then
            nodej = jnode(iedge)
          else if ( jnode(iedge) == nodei ) then
            nodej = inode(iedge)
          else
            nodej = 0
          end if
!
!  Store unvisited neighbors in STACK.
!
          if ( nodej /= 0 ) then

            if ( dad(nodej) == 0 ) then
              dad(nodej) = nodei
              nstackj = nstackj + 1
              stack(nstacki+nstackj) = nodej
            end if

          end if
 
        end do

      end do
!
!  If you picked up new neighbors on this pass, then we need to
!  search for THEIR neighbors.
!
      if ( nstackj <= 0 ) then
        exit
      end if

      stack(1:nstackj) = stack(nstacki+1:nstacki+nstackj)
      nstacki = nstackj
      nstackj = 0

    end do

  end do

  return
end
subroutine graph_arc_to_digraph_arc ( iarc, jarc, inode, jnode, maxarc, narc, &
  nedge )

!*****************************************************************************80
!
!! GRAPH_ARC_TO_DIGRAPH_ARC converts an undirected to a directed graph.
!
!  Discussion:
!
!    The intent is that every edge (I,J) of the undirected graph will 
!    become two directed edges or "arcs" (I,J) and (J,I) of the directed
!    graph.  An "arc" (I,J) is a path FROM I TO J, and does not allow
!    passage back from J to I.
!
!    An edge (I,I) results in a single arc (I,I).
!
!    If the input data already includes edges (I,J) and (J,I), then 
!    the code will catch this fact, and will produce two arcs, not four.
!
!    As part of the processing, edges (I,J) in the input array are 
!    reordered if necessary so that I <= J.  Then the edge array is 
!    sorted, and duplicates are removed.  Only then are the arcs 
!    generated.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IARC(MAXARC), JARC(MAXARC), the arcs of a 
!    directed graph, with the property that every edge (I,J) in the undirected 
!    graph corresponds to a pair of arcs (I,J) and (J,I) in the directed 
!    graph, with the exception that an edge (I,I) corresponds to a single 
!    arc (I,I).  The I-th arc goes from IARC(I) to JARC(I).
!
!    Input, integer ( kind = 4 ) INODE(NEDGE), JNODE(NEDGE), the edge array for an 
!    undirected graph.  The I-th edge connects nodes INODE(I) and JNODE(I).
!
!    Input, integer ( kind = 4 ) MAXARC, the maximum number of arcs for which storage
!    has been set aside.  MAXARC = 2*NEDGE is always enough, but less
!    space may be required if there are many duplicate edges, or 
!    edges of the form (I,I).
!
!    Output, integer ( kind = 4 ) NARC, the actual number of arcs constructed for
!    the directed graph.
!
!    Input, integer ( kind = 4 ) NEDGE, the number of edges in the undirected graph.
!
  implicit none

  integer ( kind = 4 ) maxarc
  integer ( kind = 4 ) nedge

  integer ( kind = 4 ) i
  integer ( kind = 4 ) iarc(maxarc)
  integer ( kind = 4 ) inode(nedge)
  integer ( kind = 4 ) jarc(maxarc)
  integer ( kind = 4 ) jnode(nedge)
  integer ( kind = 4 ) narc
  integer ( kind = 4 ) nuniq
!
!  Copy the edge array into the initial part of the arc array.
!
  narc = nedge

  iarc(1:narc) = inode(1:narc)
  jarc(1:narc) = jnode(1:narc)
!
!  Sort the edge array as though it were undirected.
!
  call graph_arc_edge_sort ( narc, iarc, jarc )
!
!  Eliminate duplicates.
!
  call i4vec2_uniq ( narc, iarc, jarc, nuniq )
!
!  Generate the extra arcs.
!
  narc = nuniq

  do i = 1, nuniq

    if ( iarc(i) /= jarc(i) ) then

      narc = narc + 1

      if ( narc <= maxarc ) then
        iarc(narc) = jarc(i)
        jarc(narc) = iarc(i)
      end if

    end if
 
  end do
!
!  Now sort the digraph edge array.
!
  call digraph_arc_edge_sort ( narc, iarc, jarc )

  return
end
subroutine graph_arc_to_graph_adj ( nedge, inode, jnode, adj, lda, nnode )

!*****************************************************************************80
!
!! GRAPH_ARC_TO_GRAPH_ADJ converts an arc list graph to an adjacency graph.
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
!    Output, integer ( kind = 4 ) ADJ(LDA,NNODE), the adjacency information.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of ADJ.
!
!    Output, integer ( kind = 4 ) NNODE, the number of nodes.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) nedge

  integer ( kind = 4 ) adj(lda,*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) inode(nedge)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jnode(nedge)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) mnode
  integer ( kind = 4 ) nnode
!
!  Determine the number of nodes.
!
  call graph_arc_node_count ( nedge, inode, jnode, mnode, nnode )

  if ( lda < nnode ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRAPH_ARC_TO_GRAPH_ADJ - Fatal error!'
    write ( *, '(a)' ) '  Number of nodes exceeds LDA.'
    stop
  end if

  adj(1:nnode,1:nnode) = 0

  do k = 1, nedge
    i = inode(k)
    j = jnode(k)
    adj(i,j) = 1
    adj(j,i) = 1
  end do

  return
end
subroutine graph_arc_to_graph_star ( nnode, nedge, inode, jnode, arcfir, &
  fwdarc )

!*****************************************************************************80
!
!! GRAPH_ARC_TO_GRAPH_STAR sets the forward star form of an undirected graph.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
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
!    Output, integer ( kind = 4 ) ARCFIR(NNODE+1); ARCFIR(I) is the number of the first
!    edge starting at node I in the forward star representation of the graph.
!
!    Output, integer ( kind = 4 ) FWDARC(2*NEDGE); FWDARC(I) is the ending node of
!    the I-th edge in the forward star representation of the graph.
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
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
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
subroutine graph_chro ( nnode, nedge, iendpt, iarray, jarray, karray, stack, &
  maxstack )

!*****************************************************************************80
!
!! GRAPH_CHRO calculates the chromatic polynomial of a connected graph.
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
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) NEDGE, the number of edges.
!
!    Input, integer ( kind = 4 ) IENDPT(2,NEDGE).  IENDPT(1,I) and IENDPT(2,I) are
!    the two nodes which define edge I.  On output, IENDPT has
!    been overwritten.
!
!    Output, integer ( kind = 4 ) IARRAY(NNODE).  Coefficients of the chromatic
!    polynomial in power form:
!
!      P(X) = 
!        IARRAY(N)   * X**NNODE
!      - IARRAY(N-1) * X**NNODE-1
!      + IARRAY(N-2) * X**NNODE-2
!      ...
!      +-IARRAY(1)   * X
!
!    Output, integer ( kind = 4 ) JARRAY(NNODE).  Coefficients of the chromatic
!    polynomial using the Tutte or tree form:
!
!      P(X) = SUM ( I = 1 TO NNODE ) 
!        (-1)**(NNODE-I) * IARRAY(I) * X * (X-1)**(I-1)
!
!    Output, integer ( kind = 4 ) KARRAY(NNODE).  The Stirling or factorial form of
!    chromatic polynomial.
!
!      P(X) = SUM ( I = 1 TO NNODE ) KARRAY(I) * (X)(I)
!
!    Here (X)(I) is meant to represent X*(X-1)*(X-2)...*(X-I+1).
!
!    Workspace, integer STACK(2,MAXSTACK).
!
!    Input, integer ( kind = 4 ) MAXSTACK, dimension of working storage.  An upper
!    estimate for the amount of storage required is
!    NNODE * ( IE - 0.5*(NNODE-1)).
!
  implicit none

  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) nnode
  integer ( kind = 4 ) maxstack

  integer ( kind = 4 ) i
  integer ( kind = 4 ) iarray(nnode)
  integer ( kind = 4 ) nedge1
  integer ( kind = 4 ) ien(2)
  integer ( kind = 4 ) iendpt(2,nedge)
  integer ( kind = 4 ) is
  integer ( kind = 4 ) iu
  integer ( kind = 4 ) iv
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jarray(nnode)
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) k
  integer ( kind = 4 ) karray(nnode)
  integer ( kind = 4 ) l
  integer ( kind = 4 ) nnode1
  integer ( kind = 4 ) stack(2,maxstack)

  is = 0
  jarray(1:nnode) = 0
  nedge1 = nedge
  nnode1 = nnode
 
10    continue
 
  if ( 0 < nnode1 .and. 0 < nedge1 ) then
    call span_forest ( nnode1, nedge1, iendpt, k, karray )
  else
    k = 0
  end if
 
  if ( k /= 1 ) then
    go to 50
  end if

  if ( nedge1 < nnode1 ) then
    go to 40
  end if
 
  if ( nedge1 == nnode1 ) then
 
    jarray(nnode1) = jarray(nnode1) + 1
 
  else
 
    do i = 1, nedge1
      is = is + 1
      stack(1,is) = iendpt(1,i)
      stack(2,is) = iendpt(2,i)
    end do
 
    stack(1,is) = nnode1
    stack(2,is) = nedge1 - 1
 
  end if
 
20    continue
 
  iarray(1:nnode) = 0
  iu = min ( iendpt(1,nedge1), iendpt(2,nedge1) )
  iv = iendpt(1,nedge1) + iendpt(2,nedge1) - iu
  jhi = nedge1 - 1
  nedge1 = 0
 
  do j = 1, jhi
 
    do l = 1, 2

      ien(l) = iendpt(l,j)

      if ( ien(l) == iv ) then
        ien(l) = iu
      end if

      if ( ien(l) == nnode1 ) then
        ien(l) = iv
      end if

    end do
 
    do l = 1, 2
 
      if ( ien(l) == iu ) then
        if ( iarray(ien(3-l)) /= 0 ) then
          go to 30
        end if
        iarray(ien(3-l)) = 1
      end if
 
    end do
 
    nedge1 = nedge1 + 1
 
    iendpt(1,nedge1) = ien(1)
    iendpt(2,nedge1) = ien(2)
 
30      continue
 
  end do
 
  nnode1 = nnode1 - 1
  go to 10
 
40    continue
 
  jarray(nnode1) = jarray(nnode1) + 1
 
  if ( is /= 0 ) then
 
    nnode1 = stack(1,is)
    nedge1 = stack(2,is)
    is = is - nedge1 - 1
 
    do i = 1, nedge1
      iendpt(1,i) = stack(1,is+i)
      iendpt(2,i) = stack(2,is+i)
    end do
 
    if ( nedge1 == nnode1 ) then
      jarray(nnode1) = jarray(nnode1) + 1
    else
      is = is + nedge1
      stack(1,is) = nnode1
      stack(2,is) = nedge1 - 1
    end if
 
    go to 20
 
  end if
 
50    continue
 
  do i = 1, nnode
    iarray(i) = jarray(i)
    karray(i) = ( 1 - 2 * mod ( nnode-i, 2 ) ) * jarray(i)
  end do
 
  call poly ( nnode, iarray, 1, nnode, iv )
  call poly ( nnode, karray, 0, -2, iv )
 
  return
end
subroutine graph_dist_all ( dist, dinfin, lda, nnode, path_dist )

!*****************************************************************************80
!
!! GRAPH_DIST_ALL finds the distance from every node to every other node.
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
!    Alan Gibbons,
!    Algorithmic Graph Theory,
!    Cambridge University Press, 1985,
!    ISBN 0-521-28881-9.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) DIST(LDA,NNODE). 
!
!    On input, DIST(I,J) is the length of the edge FROM node I TO node J.
!    DIST(I,J) = DINFIN if there is no direct edge from I to J.
!
!    On output, DIST has been overwritten by other information.
!
!    Input, real ( kind = 8 ) DINFIN, is a "large" number, larger than any 
!    entry in DIST, which is taken to be "infinity".  DIST(I,J) = DINFIN 
!    means there is no direct edge from node I to node J.  On output,
!    DIST(I,J) = DINFIN means there is no path from node I to node J.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of DIST and PATH_DIST,
!    which must be at least NNODE.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, real ( kind = 8 ) PATH_DIST(LDA,NNODE).  This array contains the
!    lengths of the shortest paths from each node to another node.
!    PATH_DIST(I,J) is the length of the shortest path from node I
!    to node J.  If PATH_DIST(I,J) = DINFIN, there is no path from node
!    I to node J.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) nnode

  real ( kind = 8 ) dist(lda,nnode)
  real ( kind = 8 ) dinfin
  real ( kind = 8 ) path_dist(lda,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k

  do k = 1, nnode
 
    do i = 1, nnode
      do j = 1, nnode
 
        path_dist(i,j) = dist(i,j)

        if ( dist(i,k) /= dinfin .and. dist(k,j) /= dinfin ) then

          path_dist(i,j) = min ( path_dist(i,j), dist(i,k) + dist(k,j) )

        end if
 
      end do
    end do
 
    dist(1:nnode,1:nnode) = path_dist(1:nnode,1:nnode)
 
  end do
 
  return
end
subroutine graph_dist_check ( dist, lda, nnode, ierror )

!*****************************************************************************80
!
!! GRAPH_DIST_CHECK checks a distance matrix for consistency.
!
!  Discussion:
!
!    The checks made are:
!
!      1): DIST(I,I) = 0
!      2): DIST(I,J) > 0 for I different from J
!      3): DIST(I,J) = DIST(J,I) for I different from J.
!      4): DIST(I,J) + DIST(J,K) >= DIST(I,K).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) DIST(LDA,NNODE), the distance matrix.  
!    DIST(I,J) is the distance FROM node I TO node J.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of DIST, which must be at
!    least NNODE.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no errors.
!    1, DIST(I,I) is nonzero for some I;
!    2, DIST(I,J) <= 0 for some distinct I, J
!    3, DIST(I,J) not equal to DIST(J,I) for some distinct I, J.
!    4, DIST(I,J) + DIST(J,K) < DIST(I,K) for some I, J, K.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) nnode

  real ( kind = 8 ) dist(lda,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k

  ierror = 1

  do i = 1, nnode
    if ( dist(i,i) /= 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'GRAPH_DIST_CHECK - Failed test #1:'
      write ( *, '(a,i8)' ) '  DIST(I,I) nonzero for I = ', i
      return
    end if
  end do

  ierror = 2
  do i = 1, nnode
    do j = 1, nnode
      if ( i /= j ) then
        if ( dist(i,j) <= 0.0D+00 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'GRAPH_DIST_CHECK - Failed test #2:'
          write ( *, '(a,2i8)' ) '  DIST(I,J) <= 0 for I, J = ', i, j
          return
        end if
      end if
    end do
  end do

  ierror = 3
  do i = 1, nnode
    do j = 1, i - 1
      if ( dist(i,j) /= dist(j,i) ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'GRAPH_DIST_CHECK - Failed test #3:'
        write ( *, '(a)' ) '  DIST(I,J) is not equal to DIST(J,I)'
        write ( *, '(a,2i8)' ) '  for I, J = ', i, j
        return
      end if  
    end do
  end do

  ierror = 4
  do i = 1, nnode
    do j = 1, nnode
      do k = 1, i - 1
        if ( dist(i,j) + dist(j,k) < dist(i,k) ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'GRAPH_DIST_CHECK - Failed test #4:'
          write ( *, '(a)' ) '  DIST(I,J) + DIST(J,K) < DIST(I,K)'
          write ( *, '(a,3i8)' ) '  I, J, K = ', i, j, k
          write ( *, '(a,g14.6)' ) '  DIST(I,J) = ', dist(i,j)
          write ( *, '(a,g14.6)' ) '  DIST(J,K) = ', dist(j,k)
          write ( *, '(a,g14.6)' ) '  DIST(I,K) = ', dist(i,k)
          return
        end if
      end do
    end do
  end do

  ierror = 0

  return
end
subroutine graph_dist_min_span_tree ( lda, nnode, dist, itree, jtree )

!*****************************************************************************80
!
!! GRAPH_DIST_MIN_SPAN_TREE computes a spanning tree of minimal length.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 July 2000
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
!    Input, integer ( kind = 4 ) LDA, first dimension of DIST in calling program.  
!    LDA must be at least NNODE.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, real ( kind = 8 ) DIST(LDA,NNODE).  DIST(I,J) = distance from node I
!    to node J.
!
!    Output, integer ( kind = 4 ) ITREE(NNODE-1), JTREE(NNODE-1), the pairs of nodes
!    that form the edges of the tree.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) nnode

  real ( kind = 8 ) d
  real ( kind = 8 ) dist(lda,nnode)
  real ( kind = 8 ) dmin
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imin
  integer ( kind = 4 ) it
  integer ( kind = 4 ) itree(nnode-1)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jtree(nnode-1)

  call i4vec_indicator ( nnode-1, itree )

  jtree(1:nnode-1) = -nnode
 
  do j = 1, nnode-1
!
!  Choose the node IMIN whose tree edge ( ITREE(IMIN)=IMIN, JTREE(IMIN) ) 
!  will be set.
!
    dmin = huge ( dmin )
 
    do i = 1, nnode-1
 
      it = jtree(i)
 
      if ( it < 0 ) then
 
        d = dist(-it,i)
 
        if ( d < dmin ) then
          dmin = d
          imin = i
        end if
 
      end if
 
    end do
 
    jtree(imin) = - jtree(imin)

    do i = 1, nnode-1
 
      it = jtree(i)
 
      if ( it < 0 ) then
        if ( dist(i,imin) < dist(i,-it) ) then
          jtree(i) = - imin
        end if
      end if
 
    end do
 
  end do
 
  return
end
subroutine graph_dist_min_span_tree2 ( lda, nnode, dist, class, itree, jtree )

!*****************************************************************************80
!
!! GRAPH_DIST_MIN_SPAN_TREE2 computes a spanning tree of minimal length.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LDA, first dimension of DIST in calling program.
!    LDA must be at least NNODE.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, real ( kind = 8 ) DIST(LDA,NNODE).  DIST(I,J) = distance from node I
!    to node J.
!
!    Output, integer ( kind = 4 ) CLASS(NNODE), lists the nodes in the order in
!    which they joined the tree.
!
!    Output, integer ( kind = 4 ) ITREE(NNODE-1), JTREE(NNODE-1), the pairs of nodes
!    that form the edges of the tree.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) class(nnode)
  real ( kind = 8 ) dist(lda,nnode)
  real ( kind = 8 ) dmin
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) imin
  integer ( kind = 4 ) itree(nnode-1)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) jjmin
  integer ( kind = 4 ) jmin
  integer ( kind = 4 ) jtree(nnode-1)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) npos
  logical smaller
  logical unset

  if ( nnode <= 1 ) then
    return
  end if
!
!  All the nodes start out in the negative class.
!
  npos = 0
  call i4vec_indicator ( nnode, class )
!
!  Find the shortest edge (I,J).
!
  unset = .true.
  dmin = 0.0D+00

  do i = 1, nnode
    do j = i+1, nnode

      if ( unset ) then
        smaller = .true.
      else if ( dist(i,j) < dmin ) then
        smaller = .true.
      else
        smaller = .false.
      end if

      if ( smaller ) then
        imin = i
        jmin = j
        dmin = dist(i,j)
        unset = .false.
      end if

    end do
  end do
!
!  Carry nodes IMIN and JMIN into the positive class.
!
  npos = npos + 1
  call i4_swap ( class(npos), class(imin) )

  npos = npos + 1
  call i4_swap ( class(npos), class(jmin) )

  itree(1) = imin
  jtree(1) = jmin
!
!  Now, repeatedly, find the shortest edge connecting a negative
!  and positive node.  Move the negative node to the positive class and
!  repeat.
!
  do k = 2, nnode-1

    unset = .true.
    dmin = 0.0D+00
    imin = - 99
    jmin = - 99

    do ii = 1, npos

      i = class(ii)

      do jj = npos + 1, nnode

        j = class(jj)

        if ( unset ) then
          smaller = .true.
        else if ( dist(i,j) < dmin ) then
          smaller = .true.
        else
          smaller = .false.
        end if

        if ( smaller ) then
          imin = i
          jmin = j
          jjmin = jj
          dmin = dist(i,j)
          unset = .false.
        end if

      end do

    end do

    npos = npos + 1
    call i4_swap ( class(npos), class(jjmin) )

    itree(k) = imin
    jtree(k) = jmin

  end do

  return
end
subroutine graph_dist_min_span_tree3 ( lda, nnode, dist, inode, jnode )

!*****************************************************************************80
!
!! GRAPH_DIST_MIN_SPAN_TREE3 finds a minimum spanning tree.
!
!  Discussion:
!
!    The input graph is represented by a distance matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 July 2000
!
!  Reference:
!
!    Hang Tong Lau,
!    Combinatorial Heuristic Algorithms in FORTRAN,
!    Springer Verlag, 1986.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of DIST, which should be
!    at least NNODE.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, real ( kind = 8 ) DIST(LDA,NNODE), an NNODE by NNODE distance
!    matrix.  DIST(I,J) is the distance from node I to node J.  The matrix
!    should be symmetric.  If there is no arc from node I to node J,
!    set DIST(I,J) = HUGE(1.0).
!
!    Output, integer ( kind = 4 ) INODE(NNODE), JNODE(NNODE); entries 1 through NNODE-1
!    describe the edges of the spanning tree as pairs of nodes.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) nnode

  real ( kind = 8 ) d
  real ( kind = 8 ) dist(lda,nnode)
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

  work(1:nnode) = huge ( work(1) )
  iwork1(1:nnode) = 0
  iwork2(1:nnode) = 0
!
!  Find the first non-zero arc.
!
  do ij = 1, nnode
    do kj = 1, nnode
      if ( dist(ij,kj) < huge ( dist(1,1) ) ) then
        i = ij
        go to 10
      end if
    end do
  end do

10    continue

  work(i) = 0
  iwork1(i) = 1
  tree_length = 0.0D+00
  kk4 = nnode - 1

  do jj = 1, kk4

    work(1:nnode) = huge ( work(1) )

    do i = 1, nnode
!
!  For each forward arc originating at node I calculate
!  the length of the path to node I
!
      if ( iwork1(i) == 1 ) then
        do j = 1, nnode
          if ( dist(i,j) < huge ( dist(1,1) ) .and. iwork1(j) == 0 ) then
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
subroutine graph_dist_one ( dist, dinfin, path_dist, dad, inode, path, &
  lda, nnode ) 

!*****************************************************************************80
!
!! GRAPH_DIST_ONE computes the distance from one node to all others in a graph.
!
!  Discussion:
!
!    This routine can handle both ordinary graphs and directed graphs.  
!
!    In an ordinary graph, a connection between two nodes is always guaranteed
!    to be "symmetric".  That is, if node I is connected to node J by
!    an edge of length D, then node J is connected to node I, and the
!    distance is again D.
!
!    In a directed graph, if node I is connect to node J by an edge of
!    length D, then nothing is known about a possible connection from 
!    node J back to node I.  In particular, it is possible that:
!
!    * there is no direct edge from node J to node I;
!    * the edge from node J to node I exists, but is a different "length"
!      than the edge from node I to node J.
!
!    The program computes:
!
!    * PATH_DIST, an array of distances from node INODE to all other nodes;
!
!    * DAD, an array which can be used to determine the path from
!      node INODE to any particular node.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 April 1999
!
!  Reference:
!
!    Alan Gibbons,
!    Algorithmic Graph Theory,
!    Cambridge University Press, 1985,
!    ISBN 0-521-28881-9.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) DIST(LDA,NNODE).  DIST contains the weighted
!    adjacency information defining the graph, or directed graph.
!    The diagonal entries of DIST, that is, DIST(I,I), should be set to 0.
!    The value of the typical off diagonal element DIST(I,J) should 
!    represent the length, or weight, of the edge from node I to
!    node J.  If the graph is undirected, then DIST(I,J) should always
!    equal DIST(J,I).  For a directed graph, these quantities may differ.
!    If there is no edge from node I to node J, then it would be natural
!    to set DIST(I,J) to "infinity".  Since this is not computationally
!    possible, the user must specify a special value, called DINFIN,
!    that will be used to mark such entries.  The most natural thing
!    to do would simply be to pick DINFIN to be "very large", such
!    as DINFIN = 10,000.
!    All the entries in DIST should be non-negative.  The algorithm will
!    NOT work correctly if negative edge lengths are input.
!    Off-diagonal elements DIST(I,J) may be set to zero.  This simply
!    means that two nodes are "very close", like St Paul and Minneapolis.
!
!    Input, real ( kind = 8 ) DINFIN, is a "large" number, which should be
!    larger than the length of any edge in the graph, and in fact larger
!    than the length of any reasonable path along the edges of the graph.  
!    The user should have set the DIST matrix so that DIST(I,J) = DINFIN
!    whenever there is no edge from node I to node J.  The program has to 
!    know the value of DINFIN so it can understand this information stored
!    in DIST.
!
!    Output, real ( kind = 8 ) PATH_DIST(NNODE).  On output, for every value
!    of I from 1 to NNODE, PATH_DIST(I) contains the distance from node INODE 
!    to node I in the graph.  Of course, PATH_DIST(INODE) is zero.  Moreover,
!    if PATH_DIST(I) = DINFIN, then this is the program's way of reporting that
!    there is NO path from node INODE to node I.
!
!    Output, integer ( kind = 4 ) DAD(NNODE), information defining the shortest 
!    path from node INODE to any node I, which presumably will be of
!    total distance PATH_DIST(I).
!
!    The path from node I to node INODE, is recorded "in reverse"
!    in DAD.  The last node is INODE, of course.  The previous node
!    is DAD(INODE).  The next node is DAD(DAD(INODE)) and
!    so on, until INODE itself is reached.  
!
!    If the distance from node I to node INODE is "infinity", then
!    DAD will still record a path; it's just probably of no interest.
!
!    Input, integer ( kind = 4 ) INODE, the base node, from which distances to the
!    other nodes are to be calculated.
!
!    Output, integer ( kind = 4 ) PATH(NNODE).  The value of PATH(I) records
!    the step on which the distance from INODE to node I was
!    determined.  There will be NNODE steps, and on each step
!    just one such distance is computed.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of DIST, which must be
!    at least NNODE.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) dad(nnode)
  real ( kind = 8 ) dist(lda,nnode)
  real ( kind = 8 ) dinfin
  real ( kind = 8 ) dmin
  real ( kind = 8 ) dtemp
  integer ( kind = 4 ) imin
  integer ( kind = 4 ) inode
  integer ( kind = 4 ) istep
  integer ( kind = 4 ) j
  integer ( kind = 4 ) path(nnode)
  real ( kind = 8 ) path_dist(nnode)
!
!  Initialize the data.
!
  dad(1:nnode) = inode
  path(1:nnode) = 0
  path_dist(1:nnode) = dist(inode,1:nnode)
!
!  On step 1, we connect node INODE itself.
!
  dad(inode) = inode
  path(inode) = 1
!
!  On steps ISTEP = 2 through NNODE, we try to add just one more node.
!
!  Of all the nodes which are not yet connected to INODE (because PATH
!  is 0 for this node), choose the one whose distance is least.
!
  do istep = 2, nnode
 
    dmin = dinfin
    imin = 0
 
    do j = 1, nnode
 
      if ( path(j) == 0 ) then
        if ( path_dist(j) <= dmin ) then
          dmin = path_dist(j)
          imin = j
        end if
      end if
 
    end do
!
!  If we found no new node to add, then any remaining nodes cannot
!  be connected.
!
    if ( dmin == dinfin ) then
      return
    end if
!
!  Now add the closest node, labeled IMIN, to the list.
!
    path(imin) = istep
!
!  Update the distances of the remaining unconnected nodes.
!
    do j = 1, nnode
 
      if ( path(j) == 0 ) then
 
        dtemp = path_dist(imin) + dist(imin,j)
 
        if ( dtemp < path_dist(j) ) then
          path_dist(j) = dtemp
          dad(j) = imin
        end if
 
      end if
 
    end do
 
  end do
 
  return
end
subroutine graph_dist_print ( dist, lda, nnode, title )

!*****************************************************************************80
!
!! GRAPH_DIST_PRINT prints a distance matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
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
!    Input, integer ( kind = 4 ) LDA, the leading dimension of DIST, which must be at
!    least NNODE.
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
subroutine greedy ( maxit, nodeb, noder, nnode, tol, xb, xr, yb, yr )

!*****************************************************************************80
!
!! GREEDY pairs two sets of nodes using the least total distance.
!
!  Discussion:
!
!    The method is iterative, and is not guaranteed to find the best
!    possible arrangement.  This is particulary true because it is a
!    "local" method, which only considers pairwise switches of the
!    red nodes that reduce the total distance.  This means that a
!    "locally minimizing" pairing might be found which is not the
!    global minimizer.
!
!    On the other hand, in the absence of a theoretical plan for how
!    to reach the global minimizer, the brute force search would
!    require that EVERY possible pairing be considered, and its total
!    distance computed.  This means that a total of NNODE!
!    graphs would have to be generated.
!
!    The approach used here, on each iterative step, looks at a
!    maximum of NNODE * (NNODE-1) graphs, which represents a
!    significantly more efficient method.
!
!    It would not be hard to extend this approach to a method which
!    considers switches of THREE red nodes at a time, though the
!    work there involve looking at NNODE * (NNODE-1) * (NNODE-2)
!    graphs, and as we increase the number of graphs we examine,
!    we begin to approach the NNODE! rate for the brute force
!    algorithm.
!
!    It also would not be hard to extend this method to a case where
!    there are three sets of nodes, arranged in triples, and again
!    the total distance is to be minimized.
!
!
!    If it is suspected that the pairing returned by GREEDY is only
!    a local minimizer, then the user is advised to restart the
!    calculation after randomly permuting the entries of NODER, so that
!    the routine starts from a different point in the space of graphs.
!
!    The routine is given:
!
!      an initial ordering of the black and red nodes, so that
!      ( NODEB(I), NODER(I) ) represents the I-th pair,
!
!      the X and Y coordinates of the black and red nodes,
!
!      a maximum number of iterations, and a relative distance
!      decrease requirement,
!
!    and computes:
!
!      a new ordering of the red nodes, contained in NODER, which should
!      reduce the total distance between corresponding red and black
!      nodes.
!
!
!    GREEDY can be applied to a variety of problems including:
!
!    1) We are given two sets of NNODE points, which we will call the
!       "red" and "black " groups, and the (X,Y) coordinates of each
!       point.  We may imagine these points as forming the two sets of
!       nodes of a bipartite graph lying in the (X,Y) plane.  We wish
!       to choose a pairing of red and black nodes which results in
!       the shortest total arc length.
!
!    2) We are given two sets of NNODE complex quantities, which we
!       believe are approximations to the same (unknown) set of
!       quantities.  We wish to arrange this data into NNODE pairs,
!       each containing a unique element from each set of data, which
!       minimizes the sum of squares of the discrepancies between the
!       pairs.  In particular, the two sets of data might be two
!       separate estimates of the complex eigenvalues of a matrix.
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
!  Parameters:
!
!    Input, integer ( kind = 4 ) MAXIT, the maximum number of iterations allowed.
!    Each iteration considers, one at a time, each black node, and
!    seeks to switch its red neighbor for another red neighbor that
!    reduces the total distance.
!
!    Input, integer ( kind = 4 ) NODEB(NNODE), the "labels" of the black nodes.  
!    You probably want to just set NODEB(I) = I, for i = 1 to NNODE.  
!    The entries in NODEB will not be changed.
!
!    Input/output, integer ( kind = 4 ) NODER(NNODE), the "labels" of the red nodes.  
!    You probably want to just set the input value of NODER(I) = I, 
!    for i = 1 to NNODE.  The entries in NODER WILL be changed.
!
!    At all times, the values of ( NODEB(I), NODER(I) ) contain the
!    labels of the I-th pair of black and red nodes.
!
!    On output, if GREEDY has found a better pairing of the nodes,
!    this will be reflected in the newly permuted values of NODER.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes in the black, and in the 
!    red sets.
!
!    Input, real ( kind = 8 ) TOL.
!    TOL is the relative decrease that the user demands in the
!    total distance, after each iterative step.  If we denote
!    the distance before the iterative step as OLDTOT, and the
!    distance after the iterative step as TOTAL, then the
!    routine will try another iterative step as long as "enough"
!    progress was made on this step.  Enough progress was made
!    whenever:
!
!      OLDTOT - TOTAL < TOL * TOTAL
!
!    Input, real ( kind = 8 ) XB(NNODE), the X coordinates of the black nodes.
!
!    Input, real ( kind = 8 ) XR(NNODE), the X coordinates of the red nodes.
!
!    Input, real ( kind = 8 ) YB(NNODE), the Y coordinates of the black nodes.
!
!    Input, real ( kind = 8 ) YR(NNODE), the Y coordinates of the red nodes.
!
  implicit none

  integer ( kind = 4 ) nnode

  real ( kind = 8 ) dist1
  real ( kind = 8 ) dist2
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) indx1
  integer ( kind = 4 ) indx2
  integer ( kind = 4 ) it
  integer ( kind = 4 ) maxit
  integer ( kind = 4 ) nodeb(nnode)
  integer ( kind = 4 ) nodeb1
  integer ( kind = 4 ) nodeb2
  integer ( kind = 4 ) noder(nnode)
  integer ( kind = 4 ) noder1
  integer ( kind = 4 ) noder2
  integer ( kind = 4 ) nswap
  real ( kind = 8 ) oldtot
  real ( kind = 8 ) temp
  real ( kind = 8 ) tol
  real ( kind = 8 ) total
  real ( kind = 8 ) xb(nnode)
  real ( kind = 8 ) xr(nnode)
  real ( kind = 8 ) yb(nnode)
  real ( kind = 8 ) yr(nnode)
!
!  Compute the total distance of the starting pairing.
!
  total = 0.0D+00
  do indx = 1, nnode

    nodeb1 = nodeb(indx)
    noder1 = noder(indx)

    total = total + sqrt ( &
      ( xb(nodeb1) - xr(noder1) )**2 + ( yb(nodeb1) - yr(noder1) )**2 )

  end do
 
  write ( *, '(a)' ) ' '
!
!  Begin the iterations.
!
  do it = 1, maxit

    if ( total == 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'GREEDY - Early termination.'
      write ( *, '(a)' ) '  Total discrepancy is low enough.'
      return
    end if
!
!  Save the current total distance for comparison at the end of the 
!  iteration.
!
    oldtot = total
    nswap = 0
!
!  Consider each black node, by running through indices INDX1 = 1
!  through NNODE of the NODEB array.
!
    do indx1 = 1, nnode
!
!  Get the actual labels of the current INDX1-th pair of black and 
!  red nodes.
!
      nodeb1 = nodeb(indx1)
      noder1 = noder(indx1)
!
!  Now look at the black node with INDX2 = 1 through NNODE, but ignore
!  the case where INDX1 = INDX2.
!
      do indx2 = 1, nnode
!
!  Get the labels of the current INDX2-th pair of black and red nodes.
!
        nodeb2 = nodeb(indx2)
        noder2 = noder(indx2)
 
        if ( indx2 /= indx1 ) then
!
!  Compute the total distance between (NODEB1,NODER1) and 
!  (NODEB2,NODER2), and compare it to the total where we switch the 
!  red nodes.
!
          dist1 = sqrt ( ( xb(nodeb1) - xr(noder1) )**2 &
                       + ( yb(nodeb1) - yr(noder1) )**2 ) &
                + sqrt ( ( xb(nodeb2) - xr(noder2) )**2 &
                       + ( yb(nodeb2) - yr(noder2) )**2 )

          dist2 = sqrt ( ( xb(nodeb1) - xr(noder2) )**2 &
                       + ( yb(nodeb1) - yr(noder2) )**2 ) &
                + sqrt ( ( xb(nodeb2) - xr(noder1) )**2 &
                       + ( yb(nodeb2) - yr(noder1) )**2 )
! 
!  If the new arrangement is any shorter, take it, by shuffling the
!  red nodes only, and update the total distance.
!
          if ( dist2 < dist1 ) then
            call i4_swap ( noder(indx1), noder(indx2) )
            nswap = nswap + 1
          end if
 
        end if
 
      end do
 
    end do
!
!  Now that we've checked all pairs of nodes,
!  print the new total distance, and see if we may
!  continue, or should give up.
!
    total = 0.0D+00
    do indx1 = 1, nnode

      nodeb1 = nodeb(indx1)
      noder1 = noder(indx1)

      total = total + sqrt ( ( xb(nodeb1) - xr(noder1) )**2 &
                           + ( yb(nodeb1) - yr(noder1) )**2 )

    end do

    write ( *, '(a,i8)' ) '  On step ', it
    write ( *, '(a,g14.6)' ) '  discrepancy =', total
    write ( *, '(a,i8)' ) '  Swaps made was ', nswap
 
    if ( oldtot - total <= tol * oldtot ) then
 
      temp = ( oldtot - total ) / oldtot
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'GREEDY - Warning:'
      write ( *, '(a)' ) '  The relative change in the discrepancy '
      write ( *, '(a,g14.6)' ) '  was only ', temp
      write ( *, '(a,g14.6)' ) '  which is less than the tolerance TOL =',tol
      write ( *, '(a)' ) '  Bailing out of the iteration.'
      write ( *, '(a)' ) ' '
      return
 
    end if
 
  end do
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GREEDY - Note:'
  write ( *, '(a)' ) '  The discrepancy has decreased by at least the'
  write ( *, '(a)' ) '  tolerance on every step.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Increasing the number of iterations might '
  write ( *, '(a)' ) '  provide further improvement at this rate.'
 
  return
end
subroutine grf_read ( file_name, inode, jnode, maxedge, maxnode, nedge, nnode, &
  x, y )

!*****************************************************************************80
!
!! GRF_READ reads a GRF file containing a 2D representation of a graph.
!
!  Example:
!
!    #  A graph where every node has 3 neighbors.
!    #
!    1      0.546  0.956  5      6      2    
!    2      0.144  0.650  7      3      1    
!    3      0.326  0.188  8      4      2    
!    4      0.796  0.188  9      5      3    
!    5      0.988  0.646  10     4      1    
!    6      0.552  0.814  11     12     1    
!    7      0.264  0.616  11     15     2    
!    8      0.404  0.296  15     14     3    
!    9      0.752  0.298  14     13     4    
!    10     0.846  0.624  13     12     5    
!    11     0.430  0.692  16     6      7    
!    12     0.682  0.692  17     10     6    
!    13     0.758  0.492  18     9      10   
!    14     0.566  0.358  19     8      9    
!    15     0.364  0.484  20     7      8    
!    16     0.504  0.602  11     20     17   
!    17     0.608  0.602  12     18     16   
!    18     0.634  0.510  13     19     17   
!    19     0.566  0.444  14     20     18   
!    20     0.480  0.510  15     16     19   
!
!  Discussion:
!
!    The original GRF format has been modified so that a line starting
!    with a # is considered a comment line.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the file name.
!
!    Output, integer ( kind = 4 ) INODE(MAXEDGE), JNODE(MAXEDGE), the edges.  
!    The I-th edge joins nodes INODE(I) and JNODE(I).
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN unit number associated with the
!    graph file, which should already have been opened by the user.
!
!    Input, integer ( kind = 4 ) MAXEDGE, the maximum number of edges.
!
!    Input, integer ( kind = 4 ) MAXNODE, the maximum number of nodes.
!
!    Output, integer ( kind = 4 ) NEDGE, the number of edges that were read.
!
!    Output, integer ( kind = 4 ) NNODE, the number of nodes that were read.
!
!    Output, real ( kind = 8 ) X(MAXNODE), Y(MAXNODE), the coordinates of the
!    nodes.
!
  implicit none

  integer ( kind = 4 ), parameter :: maxchr = 200

  integer ( kind = 4 ) maxedge
  integer ( kind = 4 ) maxnode

  character ( len = * ) file_name
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) inode(maxedge)
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) istring
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) jnode(maxedge)
  integer ( kind = 4 ) lchar
  integer ( kind = 4 ) nbad
  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) nnode
  integer ( kind = 4 ) nodei
  integer ( kind = 4 ) nodej
  integer ( kind = 4 ) ntext
  character ( len = maxchr ) string
  real ( kind = 8 ) x(maxnode)
  real ( kind = 8 ) xval
  real ( kind = 8 ) y(maxnode)
  real ( kind = 8 ) yval

  nbad = 0
  nedge = 0
  nnode = 0
  ntext = 0

  call get_unit ( iunit )

  open ( unit = iunit, file = file_name, status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRF_READ - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input file.'
    return
  end if
!
!  Read information about each node.
!
  do

    read ( iunit, '(a)', iostat = ios ) string

    if ( ios /= 0 ) then
      exit
    end if

    ntext = ntext + 1

    if ( len ( string ) <= 0 ) then
      cycle
    end if

    if ( string(1:1) == '#' ) then
      cycle
    end if

    istring = 1
!
!  Extract the node index, NODEI.
!
    call s_to_i4 ( string(istring:), nodei, ierror, lchar )

    if ( ierror /= 0 .or. lchar == 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'GRF_READ - Fatal error!'
      write ( *, '(a)' ) '  Unreadable node index value.'
      nbad = nbad + 1
      cycle
    end if

    istring = istring + lchar

    if ( nodei < 1 .or. maxnode < nodei ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'GRF_READ - Fatal error!'
      write ( *, '(a,i8)' ) '  Illegal node index value, NODEI = ', nodei
      cycle
    end if

    if ( nodei == nnode + 1 ) then
      nnode = nnode + 1
    else if ( nnode < nodei ) then
      nnode = nodei
    end if
!
!  Extract the X, Y coordinates of the node.
!
    call s_to_r8 ( string(istring:), xval, ierror, lchar )

    if ( ierror /= 0 .or. lchar == 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'GRF_READ - Fatal error!'
      write ( *, '(a)' ) '  Unreadable X coordinate for node.'
      nbad = nbad + 1
      cycle
    end if

    istring = istring + lchar

    call s_to_r8 ( string(istring:), yval, ierror, lchar )

    if ( ierror /= 0 .or. lchar == 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'GRF_READ - Fatal error!'
      write ( *, '(a)' ) '  Unreadable Y coordinate for node.'
      nbad = nbad + 1
      cycle
    end if

    istring = istring + lchar

    x(nodei) = xval
    y(nodei) = yval
!
!  Read the indices of the nodes to which NODEI is connected.
!
    do

      call s_to_i4 ( string(istring:), nodej, ierror, lchar )

      if ( ierror /= 0 .and. ierror /= 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'GRF_READ - Fatal error!'
        write ( *, '(a)' ) '  Unreadable node neighbor value.'
        nbad = nbad + 1
        cycle
      end if

      istring = istring + lchar

      if ( lchar <= 0 ) then
        exit
      end if

      if ( 1 <= nodej .and. nodej <= maxnode ) then

        if ( nedge < maxedge ) then
          nedge = nedge + 1
          inode(nedge) = nodei
          jnode(nedge) = nodej
        end if

      end if

      if ( maxchr < istring ) then
        exit
      end if

    end do

  end do

  close ( unit = iunit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GRF_READ - Input file statistics:'
  write ( *, '(a,i8)' ) '  Text lines:     ', ntext
  write ( *, '(a,i8)' ) '  Bad text lines: ', nbad
  write ( *, '(a,i8)' ) '  Nodes:          ', nnode
  write ( *, '(a,i8)' ) '  Edges:          ', nedge

  return
end
subroutine hqr ( nm, n, low, igh, h, wr, wi, ierr )

!*****************************************************************************80
!
!! HQR computes all eigenvalues of a real upper Hessenberg matrix.
!
!  Discussion:
!
!    This subroutine finds the eigenvalues of a real
!    upper Hessenberg matrix by the QR method.
!
!  Reference:
!
!    Martin, Peters, James Wilkinson,
!    HQR,
!    Numerische Mathematik,
!    Volume 14, pages 219-231, 1970.
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Y Ikebe, V Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NM, the leading dimension of H, which must
!    be at least N.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) LOW, IGH, two integers determined by the routine
!    BALANC.  If BALANC is not used, set LOW=1, IGH=N.
!
!    Input/output, real ( kind = 8 ) H(NM,N), the N by N upper Hessenberg
!    matrix.  Information about the transformations used in the reduction to
!    Hessenberg form by ELMHES or ORTHES, if performed, is stored
!    in the remaining triangle under the Hessenberg matrix.
!    On output, the information in H has been destroyed.
!
!    Output, real ( kind = 8 ) WR(N), WI(N), the real and imaginary parts of the
!    eigenvalues.  The eigenvalues are unordered, except that complex
!    conjugate pairs of values appear consecutively, with the eigenvalue
!    having positive imaginary part listed first.  If an error exit
!    occurred, then the eigenvalues should be correct for indices
!    IERR+1 through N.
!
!    Output, integer ( kind = 4 ) IERR, error flag.
!    0, no error.
!    J, the limit of 30*N iterations was reached while searching for
!       the J-th eigenvalue.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nm

  integer ( kind = 4 ) en
  integer ( kind = 4 ) enm2
  real ( kind = 8 ) h(nm,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) igh
  integer ( kind = 4 ) itn
  integer ( kind = 4 ) its
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) ll
  integer ( kind = 4 ) low
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  integer ( kind = 4 ) na
  real ( kind = 8 ) norm
  logical notlas
  real ( kind = 8 ) p
  real ( kind = 8 ) q
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) t
  real ( kind = 8 ) tst1
  real ( kind = 8 ) tst2
  real ( kind = 8 ) w
  real ( kind = 8 ) wi(n)
  real ( kind = 8 ) wr(n)
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) zz

  ierr = 0
  norm = 0.0D+00
  k = 1
!
!  Store roots isolated by BALANC and compute matrix norm.
!
  do i = 1, n

    do j = k, n
      norm = norm + abs ( h(i,j) )
    end do

    k = i
    if (i < low .or. igh < i ) then
      wr(i) = h(i,i)
      wi(i) = 0.0D+00
    end if

  end do

  en = igh
  t = 0.0D+00
  itn = 30 * n
!
!  Search for next eigenvalues.
!
60 continue

  if ( en < low ) then
    return
  end if

  its = 0
  na = en - 1
  enm2 = na - 1
!
!  Look for a single small sub-diagonal element.
!
70 continue

  do ll = low, en
    l = en + low - ll
    if ( l == low ) then
      exit
    end if
    s = abs ( h(l-1,l-1) ) + abs ( h(l,l) )
    if ( s == 0.0D+00 ) then
      s = norm
    end if
    tst1 = s
    tst2 = tst1 + abs ( h(l,l-1))
    if ( tst2 == tst1 ) then
      exit
    end if
  end do
!
!  Form shift.
!
  x = h(en,en)
  if ( l == en ) go to 270
  y = h(na,na)
  w = h(en,na) * h(na,en)
  if ( l == na ) go to 280

  if ( itn == 0 ) then
    ierr = en
    return
  end if
!
!  Form an exceptional shift.
!
  if ( its == 10 .or. its == 20 ) then

    t = t + x

    do i = low, en
      h(i,i) = h(i,i) - x
    end do

    s = abs ( h(en,na) ) + abs ( h(na,enm2) )
    x = 0.75D+00 * s
    y = x
    w = -0.4375D+00 * s * s

  end if

  its = its + 1
  itn = itn - 1
!
!  Look for two consecutive small sub-diagonal elements.
!
  do mm = l, enm2

    m = enm2 + l - mm
    zz = h(m,m)
    r = x - zz
    s = y - zz
    p = ( r * s - w ) / h(m+1,m) + h(m,m+1)
    q = h(m+1,m+1) - zz - r - s
    r = h(m+2,m+1)
    s = abs ( p ) + abs ( q ) + abs ( r )
    p = p / s
    q = q / s
    r = r / s

    if ( m == l ) then
      exit
    end if

    tst1 = abs ( p ) * ( abs ( h(m-1,m-1) ) + abs ( zz ) + abs ( h(m+1,m+1) ) )
    tst2 = tst1 + abs ( h(m,m-1) ) * ( abs ( q ) + abs ( r ) )

    if ( tst2 == tst1 ) then
      exit
    end if

  end do

  do i = m+2, en
    h(i,i-2) = 0.0D+00
    if ( i /= m+2 ) then
      h(i,i-3) = 0.0D+00
    end if
  end do
!
!  Double QR step involving rows l to EN and columns M to EN.
!
  do k = m, na

     notlas = k /= na

     if ( k == m ) go to 170

     p = h(k,k-1)
     q = h(k+1,k-1)
     if ( notlas ) then
       r = h(k+2,k-1)
     else
       r = 0.0D+00
     end if
     x = abs ( p ) + abs ( q ) + abs ( r )
     if ( x == 0.0D+00 ) go to 260
     p = p / x
     q = q / x
     r = r / x

170  continue

     s = sign ( sqrt ( p**2 + q**2 + r**2 ), p )

     if ( k /= m ) then
       h(k,k-1) = - s * x
     else if ( l /= m ) then
       h(k,k-1) = - h(k,k-1)
     end if

     p = p + s
     x = p / s
     y = q / s
     zz = r / s
     q = q / p
     r = r / p
     if ( notlas ) go to 225
!
!  Row modification.
!
     do j = k, n
       p = h(k,j) + q * h(k+1,j)
       h(k,j) = h(k,j) - p * x
       h(k+1,j) = h(k+1,j) - p * y
     end do

     j = min ( en, k+3 )
!
!  Column modification.
!
     do i = 1, j
       p = x * h(i,k) + y * h(i,k+1)
       h(i,k) = h(i,k) - p
       h(i,k+1) = h(i,k+1) - p * q
     end do

     go to 255

225  continue
!
!  Row modification.
!
     do j = k, n
       p = h(k,j) + q * h(k+1,j) + r * h(k+2,j)
       h(k,j) = h(k,j) - p * x
       h(k+1,j) = h(k+1,j) - p * y
       h(k+2,j) = h(k+2,j) - p * zz
     end do

     j = min ( en, k+3 )
!
!  Column modification.
!
     do i = 1, j
       p = x * h(i,k) + y * h(i,k+1) + zz * h(i,k+2)
       h(i,k) = h(i,k) - p
       h(i,k+1) = h(i,k+1) - p * q
       h(i,k+2) = h(i,k+2) - p * r
     end do

255 continue

260 continue

  end do

  go to 70
!
!  One root found.
!
270 continue

  wr(en) = x + t
  wi(en) = 0.0D+00
  en = na
  go to 60
!
!  Two roots found.
!
280 continue

  p = ( y - x ) / 2.0D+00
  q = p * p + w
  zz = sqrt ( abs ( q ) )
  x = x + t
!
!  Real root, or complex pair.
!
  if ( 0.0D+00 <= q ) then

    zz = p + sign ( zz, p )
    wr(na) = x + zz
    if ( zz == 0.0D+00 ) then
      wr(en) = wr(na)
    else
      wr(en) = x - w / zz
    end if
    wi(na) = 0.0D+00
    wi(en) = 0.0D+00

  else

    wr(na) = x + p
    wr(en) = x + p
    wi(na) = zz
    wi(en) = -zz

  end if

  en = enm2
  go to 60
end
function i4_modp ( i, j )

!*****************************************************************************80
!
!! I4_MODP returns the nonnegative remainder of integer division.
!
!  Discussion:
!
!    If
!      NREM = I4_MODP ( I, J )
!      NMULT = ( I - NREM ) / J
!    then
!      I = J * NMULT + NREM
!    where NREM is always nonnegative.
!
!    The MOD function computes a result with the same sign as the
!    quantity being divided.  Thus, suppose you had an angle A,
!    and you wanted to ensure that it was between 0 and 360.
!    Then mod(A,360.0) would do, if A was positive, but if A
!    was negative, your result would be between -360 and 0.
!
!    On the other hand, I4_MODP(A,360.0) is between 0 and 360, always.
!
!  Example:
!
!        I         J     MOD   I4_MODP   I4_MODP Factorization
!
!      107        50       7       7    107 =  2 *  50 + 7
!      107       -50       7       7    107 = -2 * -50 + 7
!     -107        50      -7      43   -107 = -3 *  50 + 43
!     -107       -50      -7      43   -107 =  3 * -50 + 43
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
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the number to be divided.
!
!    Input, integer ( kind = 4 ) J, the number that divides I.
!
!    Output, integer ( kind = 4 ) I4_MODP, the nonnegative remainder when I is divided by J.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) i4_modp

  if ( j == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_MODP - Fatal error!'
    write ( *, '(a,i8)' ) '  I4_MODP ( I, J ) called with J = ', j
    stop
  end if

  i4_modp = mod ( i, j )

  if ( i4_modp < 0 ) then
    i4_modp = i4_modp + abs ( j )
  end if

  return
end
subroutine i4_swap ( i, j )

!*****************************************************************************80
!
!! I4_SWAP switches two integer values.
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
  real       ( kind = 4 ) r
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
  real       ( kind = 4 ) r
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
subroutine i4col_compare ( lda, m, n, a, i, j, isgn )

!*****************************************************************************80
!
!! I4COL_COMPARE compares columns I and J of a integer array.
!
!  Example:
!
!    Input:
!
!      M = 3, N = 4, I = 2, J = 4
!
!      A = (
!        1  2  3  4
!        5  6  7  8
!        9 10 11 12 )
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
!    30 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the array, which must
!    be at least M.
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, integer ( kind = 4 ) A(LDA,N), an array of N columns of vectors of length M.
!
!    Input, integer ( kind = 4 ) I, J, the columns to be compared.
!    I and J must be between 1 and N.
!
!    Output, integer ( kind = 4 ) ISGN, the results of the comparison:
!    -1, column I < column J,
!     0, column I = column J,
!    +1, column I > column J.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(lda,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
!
!  Check.
!
  if ( i < 1 .or. n < i ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4COL_COMPARE - Fatal error!'
    write ( *, '(a)' ) '  Column index I is out of bounds.'
    stop
  end if

  if ( j < 1 .or. n < j ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4COL_COMPARE - Fatal error!'
    write ( *, '(a)' ) '  Column index J is out of bounds.'
    stop
  end if

  isgn = 0

  if ( i == j ) then
    return
  end if

  k = 1

  do while ( k <= m )

    if ( a(k,i) < a(k,j) ) then
      isgn = - 1
      return
    else if ( a(k,j) < a(k,i) ) then
      isgn = + 1
      return
    end if

    k = k + 1

  end do

  return
end
subroutine i4col_sort_a ( nrow, ncol, ia )

!*****************************************************************************80
!
!! I4COL_SORT_A ascending sorts an I4COL.
!
!  Definition:
!
!    In lexicographic order, the statement "X < Y", applied to two real
!    vectors X and Y of length NROW, means that there is some index I, with
!    1 <= I <= NROW, with the property that
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
!    22 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NROW, the number of rows of A, and the length of
!    a vector of data.
!
!    Input, integer ( kind = 4 ) NCOL, the number of columns of A.
!
!    Input/output, integer ( kind = 4 ) IA(NROW,NCOL).
!    On input, the array of NCOL columns of NROW-vectors.
!    On output, the columns of A have been sorted in lexicographic order.
!
  implicit none

  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nrow

  integer ( kind = 4 ) ia(nrow,ncol)
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

    call sort_heap_external ( ncol, indx, i, j, isgn )
!
!  Interchange the I and J objects.
!
    if ( 0 < indx ) then

      call i4col_swap ( nrow, ncol, ia, i, j )
!
!  Compare the I and J objects.
!
    else if ( indx < 0 ) then

      call i4col_compare ( nrow, nrow, ncol, ia, i, j, isgn )

     else

      exit

    end if

  end do

  return
end
subroutine i4col_swap ( nrow, ncol, ia, i, j )

!*****************************************************************************80
!
!! I4COL_SWAP swaps columns I and J of an I4COL.
!
!  Example:
!
!    Input:
!
!      NROW = 3, NCOL = 4, I = 2, J = 4
!
!      IA = (
!        1  2  3  4
!        5  6  7  8
!        9 10 11 12 )
!
!    Output:
!
!      IA = (
!        1  4  3  2
!        5  8  7  6
!        9 12 11 10 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NROW, NCOL, the number of rows and columns in
!    the table.
!
!    Input, integer ( kind = 4 ) IA(NROW,NCOL), a table of numbers, regarded as
!    NCOL columns of vectors of length NROW.
!
!    Input, integer ( kind = 4 ) I, J, the columns to be swapped.
!
  implicit none

  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nrow

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(nrow,ncol)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k

  if ( 1 <= i .and. i <= ncol .and. 1 <= j .and. j <= ncol ) then

    do k = 1, nrow
      call i4_swap ( ia(k,k), ia(k,j) )
    end do

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4COL_SWAP - Fatal error!'
    write ( *, '(a)' ) '  I or J is out of bounds.'
    write ( *, '(a,i8)' ) '  I =    ', i
    write ( *, '(a,i8)' ) '  J =    ', j
    write ( *, '(a,i8)' ) '  NCOL = ', ncol
    stop

  end if

  return
end
subroutine i4col_uniq ( lda, m, n, a, nuniq )

!*****************************************************************************80
!
!! I4COL_UNIQ keeps the unique elements in a sorted I4COL.
!
!  Discussion:
!
!    The array can be sorted into ascending or descending order.
!    The important point is that identical elements must be stored
!    in adjacent positions.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the array, which
!    must be at least M.
!
!    Input, integer ( kind = 4 ) M, the number of rows of A, and the length of
!    a vector of data.
!
!    Input, integer ( kind = 4 ) N, the number of columns of A.
!
!    Input/output, real ( kind = 8 ) A(LDA,N).
!    On input, the sorted array of N columns of M-vectors.
!    On output, a sorted array of NUNIQ columns of M-vectors.
!
!    Output, integer ( kind = 4 ) NUNIQ, the number of unique columns of A.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(lda,n)
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) itest
  integer ( kind = 4 ) m
  integer ( kind = 4 ) nuniq

  nuniq = 0

  if ( n <= 0 ) then
    return
  end if

  nuniq = 1

  do itest = 2, n

    call i4col_compare ( lda, m, n, a, itest, nuniq, isgn )

    if ( isgn /= 0 ) then
      nuniq = nuniq + 1
      a(1:m,nuniq) = a(1:m,itest)
    end if

  end do

  return
end
subroutine i4mat_perm ( matrix, lda, n, p )

!*****************************************************************************80
!
!! I4MAT_PERM permutes the rows and columns of a square I4MAT.
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
!    Input/output, integer ( kind = 4 ) MATRIX(LDA,N).
!    On input, the matrix to be permuted.
!    On output, the permuted matrix.
!
!    Input, integer ( kind = 4 ) LDA, the declared first dimension of MATRIX.
!    LDA must be at least N.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) P(N), the permutation.  P(I) is the new number of row
!    and column I.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) is
  integer ( kind = 4 ) it
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) k
  integer ( kind = 4 ) lc
  integer ( kind = 4 ) matrix(lda,n)
  integer ( kind = 4 ) nc
  integer ( kind = 4 ) p(n)

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
            it = matrix(i1,j1)

            do

              i1 = abs ( p(i1) )
              j1 = abs ( p(j1) )

              call i4_swap ( matrix(i1,j1), it )

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
subroutine i4mat_perm_random ( lda, n, seed, a )

!*****************************************************************************80
!
!! I4MAT_PERM_RANDOM selects a random permutation of an I4MAT.
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
!    28 March 2005
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
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the array,
!    which must be at least N.
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns in the array.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
!    Input/output, integer ( kind = 4 ) A(LDA,N), the N by N array to be permuted.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(lda,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) j
  real ( kind = 8 ) r
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
subroutine i4mat_print ( lda, m, n, a, title )

!*****************************************************************************80
!
!! I4MAT_PRINT prints an I4MAT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A.
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, integer ( kind = 4 ) A(LDA,N), the matrix to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(lda,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) m
  character ( len = * ) title

  if ( title /= ' ' ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) title
  end if

  do jlo = 1, n, 10
    jhi = min ( jlo + 9, n )
    write ( *, '(a)' ) ' '
    write ( *, '(6x,10(i7))' ) ( j, j = jlo, jhi )
    write ( *, '(a)' ) ' '
    do i = 1, m
      write ( *, '(i8,10i7)' ) i, a(i,jlo:jhi)
    end do
  end do

  return
end
subroutine i4mat_row_compare ( lda, m, n, a, row1, row2, result )

!*****************************************************************************80
!
!! I4MAT_ROW_COMPARE compares two rows of an I4MAT.
!
!  Discussion:
!
!    The rows are compared in the lexicographic sense.  They are equal
!    if every entry is equal.  Otherwise, let I be the first index 
!    where they differ.  Row 1 is less or greater than row 2 as
!    the corresponding indexed values are less or greater.
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
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A.
!    LDA must be at least M.
!
!    Input, integer ( kind = 4 ) M, number of rows in the matrix.
!
!    Input, integer ( kind = 4 ) N, number of columns in the matrix.
!
!    Input, integer ( kind = 4 ) A(LDA,N), the matrix.
!
!    Input, integer ( kind = 4 ) ROW1, ROW2, the indices of the two rows to compare.
!
!    Output, integer ( kind = 4 ) RESULT:
!    -1, ROW1 < ROW2,
!     0, ROW1 = ROW2,
!    +1, ROW1 > ROW2.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(lda,n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) m
  integer ( kind = 4 ) result
  integer ( kind = 4 ) row1
  integer ( kind = 4 ) row2

  if ( row1 < 1 .or. m < row1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4MAT_ROW_COMPARE - Fatal error!'
    write ( *, '(a)' ) '  ROW1 index out of bounds.'
    stop
  end if

  if ( row2 < 1 .or. m < row2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4MAT_ROW_COMPARE - Fatal error!'
    write ( *, '(a)' ) '  ROW2 index out of bounds.'
    stop
  end if

  result = 0

  do j = 1, n

    if ( a(row1,j) < a(row2,j) ) then
      result = -1 
      return
    else if ( a(row2,j) < a(row1,j) ) then
      result = + 1
      return
    end if

  end do

  return
end
subroutine i4mat_row_sort_d ( lda, m, n, a )

!*****************************************************************************80
!
!! I4MAT_ROW_SORT_D sorts the rows of an I4MAT into descending order.
!
!  Discussion:
!
!    Rows are compared lexicographically.  
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
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A, which must be
!    at least M.
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input/output, integer ( kind = 4 ) A(LDA,N).  On input, the M by N matrix to 
!    be row sorted.  On output, the row-sorted matrix.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(lda,n)
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) m
  integer ( kind = 4 ) row1
  integer ( kind = 4 ) row2
!
!  Initialize.
!
  indx = 0
  isgn = 0
  row1 = 0
  row2 = 0
!
!  Call the external heap sorter.
!
  do

    call sort_heap_external ( m, indx, row1, row2, isgn )
!
!  Interchange the objects.
!
    if ( 0 < indx ) then

      call i4mat_row_swap ( lda, m, n, a, row1, row2 )
!
!  Compare the objects.
!
    else if ( indx < 0 ) then

      call i4mat_row_compare ( lda, m, n, a, row1, row2, isgn )
      isgn = - isgn

    else

      exit

    end if

  end do

  return
end
subroutine i4mat_row_swap ( lda, m, n, a, row1, row2 )

!*****************************************************************************80
!
!! I4MAT_ROW_SWAP swaps two rows of an I4MAT.
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
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A.
!    LDA must be at least M.
!
!    Input, integer ( kind = 4 ) M, number of rows in the matrix.
!
!    Input, integer ( kind = 4 ) N, number of columns in the matrix.
!
!    Input/output, integer ( kind = 4 ) A(LDA,N), the matrix.
!
!    Input, integer ( kind = 4 ) ROW1, ROW2, the indices of the two rows to swap.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(lda,n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) m
  integer ( kind = 4 ) row1
  integer ( kind = 4 ) row2

  if ( row1 < 1 .or. m < row1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4MAT_ROW_SWAP - Fatal error!'
    write ( *, '(a)' ) '  ROW1 index out of bounds.'
    stop
  end if

  if ( row2 < 1 .or. m < row2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4MAT_ROW_SWAP - Fatal error!'
    write ( *, '(a)' ) '  ROW2 index out of bounds.'
    stop
  end if

  do j = 1, n
    call i4_swap ( a(row1,j), a(row2,j) )
  end do

  return
end
subroutine i4row_compare ( lda, m, n, a, i, j, isgn )

!*****************************************************************************80
!
!! I4ROW_COMPARE compares two rows of an I4ROW.
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
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the array, which must
!    be at least M.
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, integer ( kind = 4 ) A(LDA,N), an array of M rows of vectors of length N.
!
!    Input, integer ( kind = 4 ) I, J, the rows to be compared.
!    I and J must be between 1 and M.
!
!    Output, integer ( kind = 4 ) ISGN, the results of the comparison:
!    -1, row I < row J,
!     0, row I = row J,
!    +1, row I > row J.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(lda,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
!
!  Check that I and J are legal.
!
  if ( i < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4ROW_COMPARE - Fatal error!'
    write ( *, '(a)' ) '  Row index I is less than 1.'
    write ( *, '(a,i8)' ) '  I = ', i
    stop
  else if ( m < i ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4ROW_COMPARE - Fatal error!'
    write ( *, '(a)' ) '  Row index I is out of bounds.'
    write ( *, '(a,i8)' ) '  I = ', i
    write ( *, '(a,i8)' ) '  Maximum legal value is M = ', m
    stop
  end if

  if ( j < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4ROW_COMPARE - Fatal error!'
    write ( *, '(a)' ) '  Row index J is less than 1.'
    write ( *, '(a,i8)' ) '  J = ', j
    stop
  else if ( m < j ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4ROW_COMPARE - Fatal error!'
    write ( *, '(a)' ) '  Row index J is out of bounds.'
    write ( *, '(a,i8)' ) '  J = ', j
    write ( *, '(a,i8)' ) '  Maximum legal value is M = ', m
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
subroutine i4row_sort_d ( lda, m, n, a )

!*****************************************************************************80
!
!! I4ROW_SORT_D descending sorts the rows of an I4ROW.
!
!  Definition:
!
!    In lexicographic order, the statement "X < Y", applied to two real
!    vectors X and Y of length M, means that there is some index I, with
!    1 <= I <= M, with the property that
!
!      X(J) = Y(J) for J < I,
!    and
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
!    25 September 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the array,
!    which must be at least M.
!
!    Input, integer ( kind = 4 ) M, the number of rows and columns of A.
!
!    Input/output, integer ( kind = 4 ) A(LDA,N).
!    On input, the array of M rows of N-vectors.
!    On output, the rows of A have been sorted in descending
!    lexicographic order.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(lda,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  integer ( kind = 4 ) m

  if ( m <= 1 ) then
    return
  end if

  if ( n <= 0 ) then
    return
  end if
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

      call i4row_swap ( lda, m, n, a, i, j )
!
!  Compare the I and J objects.
!
    else if ( indx < 0 ) then

      call i4row_compare ( lda, m, n, a, i, j, isgn )
      isgn = - isgn

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
subroutine i4row_swap ( lda, m, n, a, row1, row2 )

!*****************************************************************************80
!
!! I4ROW_SWAP swaps two rows of an I4ROW.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the array,
!    which must be at least M.
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input/output, integer ( kind = 4 ) A(LDA,N), an array of data.
!
!    Input, integer ( kind = 4 ) ROW1, ROW2, the two rows to swap.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(lda,n)
  integer ( kind = 4 ) m
  integer ( kind = 4 ) row1
  integer ( kind = 4 ) row2
  integer ( kind = 4 ) row(n)
!
!  Check.
!
  if ( row1 < 1 .or. m < row1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4ROW_SWAP - Fatal error!'
    write ( *, '(a)' ) '  ROW1 is out of range.'
    stop
  end if

  if ( row2 < 1 .or. m < row2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4ROW_SWAP - Fatal error!'
    write ( *, '(a)' ) '  ROW2 is out of range.'
    stop
  end if

  if ( row1 == row2 ) then
    return
  end if

  row(1:n) = a(row1,1:n)
  a(row1,1:n) = a(row2,1:n)
  a(row2,1:n) = row(1:n)

  return
end
function iset2_compare ( x1, y1, x2, y2 )

!*****************************************************************************80
!
!! ISET2_COMPARE compares two I2 sets.
!
!  Discussion:
!
!    The I2 set (X1,Y1) < (X2,Y2) if
!
!      min ( X1, Y1 ) < min ( X2, Y2 ) or
!      min ( X1, Y1 ) = min ( X2, Y2 ) and max ( X1, Y1 ) < max ( X2, Y2 )
!
!    The I2 set (X1,Y1) = (X2,Y2) if
!
!      min ( X1, Y1 ) = min ( X2, Y2 ) and max ( X1, Y1 ) = max ( X2, Y2 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) X1, Y1, the first I2 set.
!
!    Input, integer ( kind = 4 ) X2, Y2, the second I2 set.
!
!    Output, character ISET2_COMPARE: '<', '>' or '=' if the first I2 set
!    is less, greater or equal to the second.
!
  implicit none

  integer ( kind = 4 ) a1
  integer ( kind = 4 ) a2
  integer ( kind = 4 ) b1
  integer ( kind = 4 ) b2
  character c
  character iset2_compare
  integer ( kind = 4 ) x1
  integer ( kind = 4 ) x2
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2

  a1 = min ( x1, y1 )
  b1 = max ( x1, y1 )

  a2 = min ( x2, y2 )
  b2 = max ( x2, y2 )

  if ( a1 < a2 ) then
    c = '<'
  else if ( a1 > a2 ) then
    c = '>'
  else if ( b1 < b2 ) then
    c = '<'
  else if ( b1 > b2 ) then
    c = '>'
  else
    c = '='
  end if

  iset2_compare = c

  return
end
subroutine iset2_index_insert_unique ( maxn, n, x, y, indx, &
  xval, yval, ival, ierror )

!*****************************************************************************80
!
!! ISET2_INDEX_INSERT_UNIQUE inserts a unique I2 set value in an indexed sorted list.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) MAXN, the maximum size of the list.
!
!    Input/output, integer ( kind = 4 ) N, the size of the list.
!
!    Input/output, integer ( kind = 4 ) X(N), Y(N), the list of I2 sets.
!
!    Input, integer ( kind = 4 ) INDX(N), the sort index of the list.
!
!    Input, integer ( kind = 4 ) XVAL, YVAL, the value to be inserted if it is
!    not already in the list.
!
!    Output, integer ( kind = 4 ) IVAL, the index in INDX corresponding to the
!    value XVAL, YVAL.
!
!    Output, integer ( kind = 4 ) IERROR, 0 for no error, 1 if an error occurred.
!
  implicit none

  integer ( kind = 4 ) maxn

  integer ( kind = 4 ) equal
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) indx(maxn)
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) less
  integer ( kind = 4 ) more
  integer ( kind = 4 ) n
  integer ( kind = 4 ) x(maxn)
  integer ( kind = 4 ) xval
  integer ( kind = 4 ) y(maxn)
  integer ( kind = 4 ) yval

  ierror = 0

  if ( n <= 0 ) then

    if ( maxn <= 0 ) then
      ierror = 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ISET2_INDEX_INSERT_UNIQUE - Fatal error!'
      write ( *, '(a)' ) '  Not enough space to store new data.'
      return
    end if

    n = 1
    x(1) = min ( xval, yval )
    y(1) = max ( xval, yval )
    indx(1) = 1
    ival = 1
    return

  end if
!
!  Does ( XVAL, YVAL ) already occur in the list?
!
  call iset2_index_search ( maxn, n, x, y, indx, xval, yval, &
    less, equal, more )

  if ( equal == 0 ) then

    if ( maxn <= n ) then
      ierror = 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ISET2_INDEX_INSERT_UNIQUE - Fatal error!'
      write ( *, '(a)' ) '  Not enough space to store new data.'
      return
    end if

    x(n+1) = min ( xval, yval )
    y(n+1) = max ( xval, yval )
    ival = more
    indx(n+1:more+1:-1) = indx(n:more:-1)
    indx(more) = n + 1
    n = n + 1

  else

    ival = equal

  end if

  return
end
subroutine iset2_index_search ( maxn, n, x, y, indx, xval, yval, &
  less, equal, more )

!*****************************************************************************80
!
!! ISET2_INDEX_SEARCH searches for an I2 set value in an indexed sorted list.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) MAXN, the maximum size of the list.
!
!    Input, integer ( kind = 4 ) N, the size of the current list.
!
!    Input, integer ( kind = 4 ) X(N), Y(N), the list.
!
!    Input, integer ( kind = 4 ) INDX(N), the sort index of the list.
!
!    Input, integer ( kind = 4 ) XVAL, YVAL, the value to be sought.
!
!    Output, integer ( kind = 4 ) LESS, EQUAL, MORE, the indexes in INDX of the
!    list entries that are just less than, equal to, and just greater
!    than the test value.  If the test value does not occur in the list,
!    then EQUAL is zero.  If the test value is the minimum entry of the
!    list, then LESS is 0.  If the test value is the greatest entry of
!    the list, then MORE is N+1.
!
  implicit none

  integer ( kind = 4 ) maxn

  character c
  integer ( kind = 4 ) equal
  integer ( kind = 4 ) hi
  integer ( kind = 4 ) indx(maxn)
  integer ( kind = 4 ) less
  integer ( kind = 4 ) lo
  integer ( kind = 4 ) mid
  integer ( kind = 4 ) more
  integer ( kind = 4 ) n
  character iset2_compare
  integer ( kind = 4 ) x(maxn)
  integer ( kind = 4 ) xhi
  integer ( kind = 4 ) xlo
  integer ( kind = 4 ) xmid
  integer ( kind = 4 ) xval
  integer ( kind = 4 ) y(maxn)
  integer ( kind = 4 ) yhi
  integer ( kind = 4 ) ylo
  integer ( kind = 4 ) ymid
  integer ( kind = 4 ) yval

  if ( n <= 0 ) then
    less = 0
    equal = 0
    more = 0
    return
  end if

  lo = 1
  hi = n

  xlo = x(indx(lo))
  ylo = y(indx(lo))

  xhi = x(indx(hi))
  yhi = y(indx(hi))

  c = iset2_compare ( xval, yval, xlo, ylo )

  if ( c == '<' ) then
    less = 0
    equal = 0
    more = 1
    return
  else if ( c == '=' ) then
    less = 0
    equal = 1
    more = 2
    return
  end if

  c = iset2_compare ( xval, yval, xhi, yhi )

  if ( c == '>' ) then
    less = n
    equal = 0
    more = n + 1
    return
  else if ( c == '=' ) then
    less = n - 1
    equal = n
    more = n + 1
    return
  end if

  do

    if ( lo + 1 == hi ) then
      less = lo
      equal = 0
      more = hi
      return
    end if

    mid = ( lo + hi ) / 2
    xmid = x(indx(mid))
    ymid = y(indx(mid))

    c = iset2_compare ( xval, yval, xmid, ymid )

    if ( c == '=' ) then
      equal = mid
      less = equal - 1
      more = equal + 1
      return
    else if ( c == '<' ) then
      hi = mid
    else if ( c == '>' ) then
      lo = mid
    end if

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
subroutine i4vec_compare ( n, a, b, isgn )

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
!      A = ( 2, 6, 2 )
!      B = ( 2, 8, 12 )
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
!    Input, integer ( kind = 4 ) A(N), B(N), the vectors to be compared.
!
!    Output, integer ( kind = 4 ) ISGN, the results of the comparison:
!    -1, A is lexicographically less than B,
!     0, A is equal to B,
!    +1, A is lexicographically greater than B.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) b(n)
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) k

  isgn = 0

  do k = 1, n

    if ( a(k) < b(k) ) then
      isgn = - 1
      return
    else if ( b(k) < a(k) ) then
      isgn = + 1
      return
    end if

  end do

  return
end
subroutine i4vec_heap_a ( n, a )

!*****************************************************************************80
!
!! I4VEC_HEAP_A reorders an array of integers into an ascending heap.
!
!  Definition:
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

10  continue
!
!  Positions 2*IFREE and 2*IFREE + 1 are the descendants of position
!  IFREE.  (One or both may not exist because they exceed N.)
!
    m = 2 * ifree
!
!  Does the first position exist?
!
    if ( m <= n ) then
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
      if ( a(m) < key ) then
        a(ifree) = a(m)
        ifree = m
        go to 10
      end if

    end if
!
!  Once there is no more shifting to do, the value KEY
!  moves into the free spot IFREE.
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
!  Definition:
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

10  continue
!
!  Positions 2*IFREE and 2*IFREE + 1 are the descendants of position
!  IFREE.  (One or both may not exist because they exceed N.)
!
    m = 2 * ifree
!
!  Does the first position exist?
!
    if ( m <= n ) then
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
      if ( key < a(m) ) then
        a(ifree) = a(m)
        ifree = m
        go to 10
      end if

    end if
!
!  Once there is no more shifting to do, the value KEY
!  moves into the free spot IFREE.
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
function i4vec_nonzero ( n, a )

!*****************************************************************************80
!
!! I4VEC_NONZERO counts the nonzero entries in an integer vector
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
!    Output, integer ( kind = 4 ) I4VEC_NONZERO, the number of nonzero entries.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4vec_nonzero

  i4vec_nonzero = 0

  do i = 1, n
    if ( a(i) /= 0 ) then
      i4vec_nonzero = i4vec_nonzero + 1
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
subroutine i4vec_perm_random ( n, seed, a )

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
!    15 April 1999
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
!    Input, integer ( kind = 4 ) N, the number of objects to be permuted.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
!    Input/output, integer ( kind = 4 ) A(N), the vector to be permuted.
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
!! I4VEC_REVERSE reverses the elements of an integer vector.
!
!  Example:
!
!    Input:
!
!      N = 5, A = ( 11, 12, 13, 14, 15 ).
!
!    Output:
!
!      A = ( 15, 14, 13, 12, 11 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 1998
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
subroutine i4vec_rotate ( n, m, a )

!*****************************************************************************80
!
!! I4VEC_ROTATE rotates an object in place.
!
!  Example:
!
!    Input:
!
!      N = 5, M = 2
!      A = ( 1, 2, 3, 4, 5 )
!
!    Output:
!
!      A = ( 4, 5, 1, 2, 3 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of objects.
!
!    Input, integer ( kind = 4 ) M, the number of positions to the right that
!    each element should be moved.  Elements that shift pass position
!    N "wrap around" to the beginning of the array.
!
!    Input/output, integer ( kind = 4 ) A(N), the array to be rotated.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ) iget
  integer ( kind = 4 ) iput
  integer ( kind = 4 ) istart
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mcopy
  integer ( kind = 4 ) nset
  integer ( kind = 4 ) temp
!
!  Force M to be positive, between 0 and N-1.
!
  mcopy = i4_modp ( m, n )

  if ( mcopy == 0 ) then
    return
  end if

  istart = 0
  nset = 0

  do

    istart = istart + 1

    if ( n < istart ) then
      exit
    end if

    temp = a(istart)
    iget = istart
!
!  Copy the new value into the vacated entry.
!
    do

      iput = iget

      iget = iget - mcopy

      if ( iget < 1 ) then
        iget = iget + n
      end if

      if ( iget == istart ) then
        exit
      end if

      a(iput) = a(iget)
      nset = nset + 1

    end do

    a(iput) = temp
    nset = nset + 1

    if ( n <= nset ) then
      exit
    end if

  end do

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
subroutine i4vec_sort_heap_index_d ( n, a, indx )

!*****************************************************************************80
!
!! I4VEC_SORT_HEAP_INDEX_D does an indexed heap descending sort of an integer vector.
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
!      call I4VEC_PERMUTE ( N, A, INDX )
!
!    after which A(I), I = 1 to N is sorted.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
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
!    Input, integer ( kind = 4 ) A(N), an array to be index-sorted.
!
!    Output, integer ( kind = 4 ) INDX(N), contains the sort index.  The
!    I-th element of the sorted array is A(INDX(I)).
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) aval
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) indxt
  integer ( kind = 4 ) ir
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l

  call i4vec_indicator ( n, indx )

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
        return
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
subroutine i4vec_uniq ( n, a, nuniq )

!*****************************************************************************80
!
!! I4VEC_UNIQ finds the number of unique elements in a sorted integer array.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements in A.
!
!    Input/output, integer ( kind = 4 ) A(N).  On input, the sorted
!    integer ( kind = 4 ) array.  On output, the unique elements in A.
!
!    Output, integer ( kind = 4 ) NUNIQ, the number of unique elements in A.
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
  
  do itest = 2, n

    if ( a(itest) /= a(nuniq) ) then
      nuniq = nuniq + 1
      a(nuniq) = a(itest)
    end if
 
  end do

  return
end
subroutine i4vec2_compare ( n, i4vec, jvec, i, j, isgn )

!*****************************************************************************80
!
!! I4VEC2_COMP compares pairs of integers stored in two vectors.
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
!    Input, integer ( kind = 4 ) I4VEC(N), JVEC(N), contain the two components of each item.
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
  integer ( kind = 4 ) i4vec(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jvec(n)

  isgn = 0

  if ( i4vec(i) < i4vec(j) ) then
    isgn = -1
  else if ( i4vec(i) == i4vec(j) ) then
    if ( jvec(i) < jvec(j) ) then
      isgn = -1
    else if ( jvec(i) < jvec(j) ) then
      isgn = 0
    else if ( jvec(j) < jvec(i) ) then
      isgn = +1
    end if
  else if ( i4vec(j) < i4vec(i) ) then
    isgn = +1
  end if

  return
end
subroutine i4vec2_print ( n, a, b, title )

!*****************************************************************************80
!
!! I4VEC2_PRINT prints a pair of integer vectors.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, integer ( kind = 4 ) A(N), B(N), the vectors to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) b(n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  if ( title /= ' ' ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(i8,2i10)' ) i, a(i), b(i)
  end do

  return
end
subroutine i4vec2_sort_a ( n, a1, a2 )

!*****************************************************************************80
!
!! I4VEC2_SORT_A ascending sorts a vector of pairs of integers.
!
!  Discussion:
!
!    Each item to be sorted is a pair of integers (I,J), with the I
!    and J values stored in separate vectors A1 and A2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items of data.
!
!    Input/output, integer ( kind = 4 ) A1(N), A2(N), the data to be sorted..
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a1(n)
  integer ( kind = 4 ) a2(n)
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

    call sort_heap_external ( n, indx, i, j, isgn )
!
!  Interchange the I and J objects.
!
    if ( 0 < indx ) then

      call i4_swap ( a1(i), a1(j) )
      call i4_swap ( a2(i), a2(j) )
!
!  Compare the I and J objects.
!
    else if ( indx < 0 ) then

      call i4vec2_compare ( n, a1, a2, i, j, isgn )

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
subroutine i4vec2_sort_d ( n, a1, a2 )

!*****************************************************************************80
!
!! I4VEC2_SORT_D descending sorts a vector of pairs of integers.
!
!  Discussion:
!
!    Each item to be sorted is a pair of integers (I,J), with the I
!    and J values stored in separate vectors A1 and A2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items of data.
!
!    Input/output, integer ( kind = 4 ) A1(N), A2(N), the data to be sorted..
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a1(n)
  integer ( kind = 4 ) a2(n)
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

    call sort_heap_external ( n, indx, i, j, isgn )
!
!  Interchange the I and J objects.
!
    if ( 0 < indx ) then

      call i4_swap ( a1(i), a1(j) )
      call i4_swap ( a2(i), a2(j) )
!
!  Compare the I and J objects.
!
    else if ( indx < 0 ) then

      call i4vec2_compare ( n, a1, a2, i, j, isgn )
      isgn = - isgn

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
subroutine i4vec2_uniq ( n, a1, a2, nuniq )

!*****************************************************************************80
!
!! I4VEC2_UNIQ keeps the unique elements in a array of pairs of integers.
!
!  Discussion:
!
!    Item I is stored as the pair A1(I), A2(I).
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
!    09 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items.
!
!    Input/output, integer ( kind = 4 ) A1(N), A2(N).
!    On input, the array of N items.
!    On output, an array of NUNIQ unique items.
!
!    Output, integer ( kind = 4 ) NUNIQ, the number of unique items.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a1(n)
  integer ( kind = 4 ) a2(n)
  integer ( kind = 4 ) itest
  integer ( kind = 4 ) nuniq

  nuniq = 0

  if ( n <= 0 ) then
    return
  end if

  nuniq = 1

  do itest = 2, n

    if ( a1(itest) /= a1(nuniq) .or. a2(itest) /= a2(nuniq) ) then

      nuniq = nuniq + 1

      a1(nuniq) = a1(itest)
      a2(nuniq) = a2(itest)

    end if

  end do

  return
end
subroutine ksub_random ( n, k, seed, iarray )

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
!    28 March 2005
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
!    Input, integer ( kind = 4 ) N, the size of the set from which subsets are drawn.
!
!    Input, integer ( kind = 4 ) K, number of elements in desired subsets.  K must
!    be between 0 and N.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
!    Output, integer ( kind = 4 ) IARRAY(K).  IARRAY(I) is the I-th element of the
!    output set.  The elements of IARRAY are in order.
!
  implicit none

  integer ( kind = 4 ) k

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) iarray(k)
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
  else if ( n < k ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KSUB_RANDOM - Fatal error!'
    write ( *, '(a,i8)' ) '  N = ', n
    write ( *, '(a,i8)' ) '  K = ', k
    write ( *, '(a)' ) '  K <= N is required!'
    stop
  end if

  if ( k == 0 ) then
    return
  end if

  do i = 1, k
    iarray(i) = ( ( i - 1 ) * n ) / k
  end do

  do i = 1, k

    do

      ix = i4_uniform ( 1, n, seed )

      l = 1 + ( ix * k - 1 ) / n

      if ( iarray(l) < ix ) then
        exit
      end if

    end do

    iarray(l) = iarray(l) + 1

  end do

  ip = 0
  is = k

  do i = 1, k

    m = iarray(i)
    iarray(i) = 0

    if ( m /= ( (i-1) * n ) / k ) then
      ip = ip + 1
      iarray(ip) = m
    end if

  end do

  ihi = ip

  do i = 1, ihi
    ip = ihi + 1 - i
    l = 1 + ( iarray(ip) * k - 1 ) / n
    ids = iarray(ip) - ( ( l - 1 ) * n ) / k
    iarray(ip) = 0
    iarray(is) = l
    is = is - ids
  end do

  do ll = 1, k

    l = k + 1 - ll

    if ( iarray(l) /= 0 ) then
      ir = l
      m0 = 1 + ( ( iarray(l) - 1 ) * n ) / k
      m = ( iarray(l) * n ) / k - m0 + 1
    end if

    ix = i4_uniform ( m0, m0+m-1, seed )

    i = l + 1

    do while ( i <= ir )

      if ( ix < iarray(i) ) then
        exit
      end if

      ix = ix + 1
      iarray(i-1) = iarray(i)
      i = i + 1

    end do

    iarray(i-1) = ix
    m = m - 1

  end do

  return
end
subroutine m_graph_adj_edge_seq ( adj, lda, nnode, edge_seq )

!*****************************************************************************80
!
!! M_GRAPH_ADJ_EDGE_SEQ computes the edge sequence of a multigraph.
!
!  Discussion:
!
!    The edge sequence of a multigraph may be constructed by sorting the
!    entries of each row of the adjacency matrix in descending order, and 
!    then sorting the rows themselves in descending order.
!
!    If two multigraphs are isomorphic, they must have the same edge sequence.
!
!    If two multigraphs have different edge sequences, they cannot be
!    isomorphic.
!
!  Example:
!
!    ADJ = 
!       0 1 2 3
!       1 0 2 0
!       2 2 0 1
!       3 0 1 0
!
!    EDGE_SEQ =
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
!    Input, integer ( kind = 4 ) ADJ(LDA,NNODE), the adjacency information.
!    ADJ(I,J) is 1 if there is an edge from node I to node J.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the ADJ array,
!    which must be at least NNODE.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) EDGE_SEQ(LDA,NNODE), the degree sequence of the graph.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) adj(lda,nnode)
  integer ( kind = 4 ) edge_seq(lda,nnode)
!
!  Copy the adjacency matrix.
!
  edge_seq(1:nnode,1:nnode) = adj(1:nnode,1:nnode)
!
!  Descending sort the elements of each row.
!
  call i4row_sort_d ( lda, nnode, nnode, edge_seq )
!
!  Sort the rows of the matrix.
!
  call i4mat_row_sort_d ( lda, nnode, nnode, edge_seq )

  return
end
subroutine maze_diam ( bar, degree, diam, flat, m, n, path, istart, jstart, &
  istop, jstop )

!*****************************************************************************80
!
!! MAZE_DIAM computes the "diameter" of a maze that has no circuits.
!
!  Discussion:
!
!    The routine also returns two cells, (ISTART,JSTART), and (ISTOP,JSTOP)
!    which are separated by a path of length DIAM.
!
!  Definition:
!
!    The diameter is the length of the longest path that never passes
!    through the same cell twice.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 August 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) BAR(M,N+1), records the vertical "bars" in the maze.
!    -1, means "indefinite", that there is no cell of the maze on either
!        side of this position;
!     0, means "wall", that there is a cell on at least one side, and
!        a wall here;
!     1, means "open", that there are cells on both sides (or possibly
!        an opening to the exterior) and the way is open.
!
!    Output, integer ( kind = 4 ) DEGREE(M,N), the degree of each node.
!
!    Output, integer ( kind = 4 ) DIAM, the length of the longest path in the tree.
!
!    Input, integer ( kind = 4 ) FLAT(M+1,N), records the horizontal "flats" in the maze.
!    -1, means "indefinite", that there is no cell of the maze on either
!        side of this position;
!     0, means "wall", that there is a cell on at least one side, and
!        a wall here;
!     1, means "open", that there are cells on both sides (or possibly
!        an opening to the exterior) and the way is open.
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns of cells.
!
!    Output, integer ( kind = 4 ) PATH(M,N), marks the path between the cells
!    (ISTART,JSTART) and (ISTOP,JSTOP).  A cell (I,J) is in the path
!    if PATH(I,J) is 1.
!
!    Output, integer ( kind = 4 ) ISTART, JSTART, are the I and J cell coordinates of the
!    starting cell.
!
!    Output, integer ( kind = 4 ) ISTOP, JSTOP, are the I and J cell coordinates of the
!    goal cell.
!
  implicit none

  integer ( kind = 4 ), parameter :: OPEN = 1
  integer ( kind = 4 ), parameter :: SHUT = 2

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) bar(m,n+1)
  integer ( kind = 4 ) degree(m,n)
  integer ( kind = 4 ) diam
  integer ( kind = 4 ) flat(m+1,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) invals
  integer ( kind = 4 ) istart
  integer ( kind = 4 ) istop
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) jstart
  integer ( kind = 4 ) jstop
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kstep
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) path(m,n)

  if ( m * n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MAZE_DIAM - Fatal error!'
    write ( *, '(a)' ) '  M*N <= 0.'
    stop
  else if ( m * n == 1 ) then
    diam = 0
    return
  end if

  k = 0
  do j = 1, n
    do i = 1, m
      k = k + 1
      path(i,j) = k
    end do
  end do
!
!  On step KSTEP:
!
!    Identify the terminal and interior nodes.
!
!    If there are no interior nodes left,
!
!      then there are just two nodes left at all.  The diameter is 2*K-1,
!      and a maximal path extends between the nodes whose labels are
!      contained in the two remaining terminal nodes.
!
!    Else
!
!      The label of each terminal node is passed to its interior neighbor.
!      If more than one label arrives, take any one.
!
!      The terminal nodes are removed.
!
  kstep = 0

10    continue

  kstep = kstep + 1
!
!  Compute the degree of each node.
!
  do j = 1, n
    do i = 1, m

      degree(i,j) = 0

      if ( flat(i,j) == OPEN ) then
        degree(i,j) = degree(i,j) + 1
      end if
      if ( flat(i+1,j) == OPEN ) then
        degree(i,j) = degree(i,j) + 1
      end if
      if ( bar(i,j) == OPEN ) then
        degree(i,j) = degree(i,j) + 1
      end if
      if ( bar(i,j+1) == OPEN ) then
        degree(i,j) = degree(i,j) + 1
      end if

    end do
  end do
!
!  Count the number of interior nodes.
!
  invals = 0
  do j = 1, n
    do i = 1, m
      if ( 1 < degree(i,j) ) then
        invals = invals + 1
      end if
    end do
  end do
!
!  If there are at least two interior nodes, then chop off the
!  terminal nodes and pass their labels inward.
!
  if ( 2 <= invals ) then

    k = 0

    do j = 1, n

      do i = 1, m

        k = k + 1

        if ( degree(i,j) == 1 ) then

          if ( flat(i,j) == OPEN ) then
            i2 = i - 1
            j2 = j
            flat(i,j) = SHUT
          else if ( flat(i+1,j) == OPEN ) then
            i2 = i + 1
            j2 = j
            flat(i+1,j) = SHUT
          else if ( bar(i,j) == OPEN ) then
            i2 = i
            j2 = j - 1
            bar(i,j) = SHUT
          else if ( bar(i,j+1) == OPEN ) then
            i2 = i
            j2 = j + 1
            bar(i,j+1) = SHUT
          end if

          path(i2,j2) = path(i,j)

        end if

      end do

    end do

    go to 10
!
!  But if there are 1 or 0 interior nodes, it's time to stop.
!
  else if ( invals == 1 ) then

    diam = 2 * kstep + 2

  else if ( invals == 0 ) then

    diam = 2 * kstep + 1

  end if
!
!  Now get the labels from two of the remaining terminal nodes.
!  The nodes represented by these labels will be a diameter apart.
!
  n1 = 0
  n2 = 0

  do j = 1, n
    do i = 1, m

      if ( degree(i,j) == 1 ) then
        if ( n1 == 0 ) then
          n1 = path(i,j)
        else if ( n2 == 0 ) then
          n2 = path(i,j)
        end if
      end if

    end do
  end do
!
!  Set the labels of the interior node (if any) and nodes marked
!  N1 and N2 to 1, and all others to 0.  This will label the nodes on the path.
!
  if ( invals == 1 ) then

    do j = 1, n
      do i = 1, m
        if ( 1 < degree(i,j) ) then
          path(i,j) = 1
        end if
      end do
    end do

  end if

  do j = 1, n
    do i = 1, m

      if ( path(i,j) == n1 .or. path(i,j) == n2 ) then
        path(i,j) = 1
      else
        path(i,j) = 0
      end if

    end do
  end do
!
!  Translate N1 and N2 to row, column.
!
  jstart = ( n1 - 1 ) / m + 1
  istart = n1 - ( jstart - 1 ) * m

  jstop = ( n2 - 1 ) / m + 1
  istop = n2 - ( jstop - 1 ) * m
!
!  Clean up the DEGREE and LINKS arrays.
!
  do i = 1, m
    do j = 1, n+1
      if ( bar(i,j) == SHUT ) then
        bar(i,j) = OPEN
      end if
    end do
  end do

  do i = 1, m+1
    do j = 1, n
      if ( flat(i,j) == SHUT ) then
        flat(i,j) = OPEN
      end if
    end do
  end do

  do j = 1, n
    do i = 1, m

      degree(i,j) = 0

      if ( flat(i,j) == OPEN ) then
        degree(i,j) = degree(i,j) + 1
      end if
      if ( flat(i+1,j) == OPEN ) then
        degree(i,j) = degree(i,j) + 1
      end if
      if ( bar(i,j) == OPEN ) then
        degree(i,j) = degree(i,j) + 1
      end if
      if ( bar(i,j+1) == OPEN ) then
        degree(i,j) = degree(i,j) + 1
      end if

    end do
  end do

  return
end
subroutine maze_path ( bar, flat, m, n, istart, jstart, istop, jstop )

!*****************************************************************************80
!
!! MAZE_PATH finds a path through a maze.
!
!  Warning: 
!
!    This routine has some stupid internal limits which could
!    be fixed by reprogramming.  (Use the BAR and FLAT arrays to record
!    the tentative path, for instance.)
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
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) BAR(M,N+1), records the vertical "bars" in the maze,
!    and on output, the path through open bars:
!
!    -1, means "indefinite", that there is no cell of the maze on either
!        side of this position;
!     0, means "wall", that there is a cell on at least one side, and
!        a wall here;
!     1, means "open", that there are cells on both sides (or possibly
!        an opening to the exterior) and the way is open.
!     2, means the path goes through this open bar.
!
!    Input/output, integer ( kind = 4 ) FLAT(M+1,N), records the horizontal "flats" in the 
!    maze, and on output, the path through open flats:
!
!    -1, means "indefinite", that there is no cell of the maze on either
!        side of this position;
!     0, means "wall", that there is a cell on at least one side, and
!        a wall here;
!     1, means "open", that there are cells on both sides (or possibly
!        an opening to the exterior) and the way is open.
!     2, means the path goes through this open flat.
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns of cells.
!
!    Input, integer ( kind = 4 ) ISTART, JSTART, are the I and J cell coordinates of the 
!    starting cell.
!
!    Input, integer ( kind = 4 ) ISTOP, JSTOP, are the I and J cell coordinates of the 
!    goal cell, which will be required to be a terminal node of the tree.
!
  implicit none

  integer ( kind = 4 ), parameter :: maxpath = 200
  integer ( kind = 4 ), parameter :: maxstack = 500
  integer ( kind = 4 ), parameter :: maxused = 500
  integer ( kind = 4 ), parameter :: OPEN = 1
  integer ( kind = 4 ), parameter :: PATH = 2

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) bar(m,n+1)
  integer ( kind = 4 ) flat(m+1,n)
  integer ( kind = 4 ) ipath
  integer ( kind = 4 ) istart
  integer ( kind = 4 ) istop
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) ival2
  integer ( kind = 4 ) jstart
  integer ( kind = 4 ) jstop
  integer ( kind = 4 ) jval
  integer ( kind = 4 ) jval2
  integer ( kind = 4 ) kval
  integer ( kind = 4 ) kval2
  integer ( kind = 4 ) ncan
  integer ( kind = 4 ) npath
  integer ( kind = 4 ) nstack
  integer ( kind = 4 ) pathlist(maxpath)
  integer ( kind = 4 ) stack(maxstack)
  integer ( kind = 4 ) used(maxused)

  if ( istart < 1 .or. m < istart ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MAZE_PATH - Fatal error!'
    write ( *, '(a,i8)' ) '  ISTART out of range, = ', istart
    write ( *, '(a,i8)' ) '  Must be between 1 and ', m
    stop
  else if ( jstart < 1 .or. n < jstart ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MAZE_PATH - Fatal error!'
    write ( *, '(a,i8)' ) '  JSTART out of range, = ', jstart
    write ( *, '(a,i8)' ) '  Must be between 1 and ', n
    stop
  else if ( istop < 1 .or. m < istop ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MAZE_PATH - Fatal error!'
    write ( *, '(a,i8)' ) '  ISTOP out of range, = ', istop
    write ( *, '(a,i8)' ) '  Must be between 1 and ', m
    stop
  else if ( jstop < 1 .or. n < jstop ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MAZE_PATH - Fatal error!'
    write ( *, '(a,i8)' ) '  JSTOP out of range, = ', jstop
    write ( *, '(a,i8)' ) '  Must be between 1 and ', n
    stop
  end if

  if ( maxused < m * n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MAZE_PATH - Fatal error!'
    write ( *, '(a)' ) '  M * N greater than internal limit MAXUSED.'
    stop
  end if

  used(1:m*n) = 0
  pathlist(1:m*n) = 0
!
!  Begin the path at (ISTART,JSTART).
!
  npath = 1
  ival = istart
  jval = jstart
  kval = ( jval - 1 ) * m + ival
  pathlist(npath) = kval
  used(kval) = npath

  ncan = 0
  nstack = 1
  stack(nstack) = ncan
!
!  Try to take a new step.
!
10    continue
!
!  Find all the accessible never-used neighbors of the current endpoint.
!  Add them to the stack, and set NCAN to the number of candidates.
!
  ncan = 0

  if ( ival /= 1 ) then

    if ( flat ( ival, jval ) == OPEN ) then

      ival2 = ival - 1
      jval2 = jval
      kval2 = ( jval2 - 1 ) * m + ival2

      if ( used(kval2) == 0 ) then
        ncan = ncan + 1
        nstack = nstack + 1
        if ( maxstack < nstack ) then
          go to 100
        end if
        stack(nstack) = kval2
      end if

    end if

  end if

  if ( jval /= n ) then

    if ( bar ( ival, jval+1 ) == OPEN ) then

      ival2 = ival
      jval2 = jval + 1
      kval2 = ( jval2 - 1 ) * m + ival2

      if ( used(kval2) == 0 ) then
        ncan = ncan + 1
        nstack = nstack + 1
        if ( maxstack < nstack ) then
          go to 100
        end if
        stack(nstack) = kval2
      end if

    end if

  end if

  if ( jval /= 1 ) then

    if ( bar ( ival, jval ) == OPEN ) then
      ival2 = ival
      jval2 = jval - 1
      kval2 = ( jval2 - 1 ) * m + ival2

      if ( used(kval2) == 0 ) then
        ncan = ncan + 1
        nstack = nstack + 1
        if ( maxstack < nstack ) then
          go to 100
        end if
        stack(nstack) = kval2
      end if

    end if

  end if

  if ( ival /= m ) then

    if ( flat ( ival+1, jval ) == OPEN ) then

      ival2 = ival + 1
      jval2 = jval
      kval2 = ( jval2 - 1 ) * m + ival2

      if ( used(kval2) == 0 ) then
        ncan = ncan + 1
        nstack = nstack + 1
        if ( maxstack < nstack ) then
          go to 100
        end if
        stack(nstack) = kval2
      end if

    end if

  end if
!
!  Add NCAN to the stack.
!
  nstack = nstack + 1
  if ( maxstack < nstack ) then
    go to 100
  end if
  stack(nstack) = ncan

20    continue
!
!  If NCAN=0, then...
!
  if ( ncan == 0 ) then
!
!  ...if the current cell is the starting point, we've failed.
!
    if ( npath == 1 ) then

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MAZE_PATH - Note'
      write ( *, '(a)' ) '  Could not find a path to the goal.'

      return
!
!  ...Else drop the current endpoint, going back to previous cell
!  on the path, pop the stack one level, (getting new value of NCAN), 
!  go to 20.
!
    else

      used(kval) = - used(kval)

      npath = npath - 1
      kval = pathlist(npath)
      ival = mod ( kval, m )
      jval = 1 + ( kval - ival ) / m 

      nstack = nstack - 1
      ncan = stack(nstack)
      go to 20

    end if
!
!  Else, take one candidate off the stack, add it to the path,
!  mark it as used, set NCAN = NCAN-1.
!
  else 

    kval = stack(nstack-1)

    npath = npath + 1

    if ( maxpath < npath ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MAZE_PATH - Fatal error!'
      write ( *, '(a)' ) '  NPATH exceeds internal limit MAXPATH.'
      stop
    end if

    pathlist(npath) = kval

    used(kval) = npath

    jval = ( kval - 1 ) / m + 1
    ival = kval - ( jval - 1 ) * m

    ncan = ncan - 1
    nstack = nstack-1
    stack(nstack) = ncan
!
!  If the candidate is not the goal, go to 10...
!
    if ( ival /= istop .or. jval /= jstop ) then
      go to 10
    end if
!
!  ...else we're done.
!
    do ipath = 1, npath-1

      kval = pathlist(ipath)
      jval = ( kval - 1 ) / m + 1
      ival = kval - ( jval - 1 ) * m

      kval2 = pathlist(ipath+1)

      if ( kval2 == kval - 1 ) then
        flat(ival,jval) = PATH
      else if ( kval2 == kval + m ) then
        bar(ival,jval+1) = PATH
      else if ( kval2 == kval - m ) then
        bar(ival,jval) = PATH
      else if ( kval2 == kval + 1 ) then
        flat(ival+1,jval) = PATH
      end if

    end do

    return

  end if

100   continue

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MAZE_PATH - Fatal error!'
  write ( *, '(a)' ) '  The size of the internal stack was exceeded!'
  stop
end
subroutine maze_print ( bar, flat, m, n, istart, jstart, istop, jstop, title )

!*****************************************************************************80
!
!! MAZE_PRINT prints out a maze and a path.
!
!  Example:
!
!    +--+--+
!    |*****|$$
!    +**+**+**+
!    |00|*****|
!    +  +--+--+
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
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
!    Input, integer ( kind = 4 ) BAR(M,N+1), records the vertical "bars" in the maze.
!
!    -1, means "indefinite", that there is no cell of the maze on either
!        side of this position;
!     0, means "wall", that there is a cell on at least one side, and
!        a wall here;
!     1, means "open", that there are cells on both sides (or possibly
!        an opening to the exterior) and the way is open.
!     2, means "path", that the way is open, and the path goes this way.
!
!    Input, integer ( kind = 4 ) FLAT(M+1,N), records the horizontal "flats" in the maze.
!
!    -1, means "indefinite", that there is no cell of the maze on either
!        side of this position;
!     0, means "wall", that there is a cell on at least one side, and
!        a wall here;
!     1, means "open", that there are cells on both sides (or possibly
!        an opening to the exterior) and the way is open.
!     2, means "path", that the way is open, and the path goes this way.
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns of cells.  Currently,
!    the program cannot handle a maze with more than 26 columns.
!
!    Input, integer ( kind = 4 ) ISTART, JSTART, are the I and J cell coordinates of the
!    starting cell.  The starting cell will be marked "00".  If no
!    starting cell is to be specified, set ISTART = JSTART = 0.
!
!    Input, integer ( kind = 4 ) ISTOP, JSTOP, are the I and J cell coordinates of the
!    goal cell.  The goal cell will be marked "$$".  If no goal cell
!    is to be specified, set ISTOP = JSTOP = 0.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: NMAX = 26
  integer ( kind = 4 ), parameter :: INDEF = -1
  integer ( kind = 4 ), parameter :: WALL = 0
  integer ( kind = 4 ), parameter :: OPEN = 1
  integer ( kind = 4 ), parameter :: PATH = 2

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) bar(m,n+1)
  integer ( kind = 4 ) flat(m+1,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) istart
  integer ( kind = 4 ) istop
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jstart
  integer ( kind = 4 ) jstop
  integer ( kind = 4 ) nsafe
  character ( len = 3*(NMAX+1) ) string
  character ( len = * ) title

  if ( NMAX < n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MAZE_PRINT - Warning!'
    write ( *, '(a,i8)' ) '  N may not be more than ', NMAX
    write ( *, '(a)' ) '  Only a portion of the maze will be shown.'
  end if

  if ( len_trim ( title ) /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  nsafe = min ( n, NMAX )

  do i = 1, m

    string = ' '

    ilo = 1
    do j = 1, nsafe

      if ( flat(i,j) == WALL ) then

        string(ilo:ilo+3) = '+--+'

      else if ( flat(i,j) == OPEN ) then

        string(ilo:ilo+3) = '+  +'

      else if ( flat(i,j) == PATH ) then

        string(ilo:ilo+3) = '+**+'

      else if ( flat(i,j) == INDEF ) then

      end if

      ilo = ilo + 3

    end do

    write ( *, '(a)' ) string(1:ilo)

    string = ' '

    ilo = 1
    do j = 1, nsafe+1

      if ( bar(i,j) == WALL ) then

        string(ilo:ilo) = '|'

      else if ( bar(i,j) == OPEN ) then

        string(ilo:ilo) = ' '

      else if ( bar(i,j) == PATH ) then

        string(ilo:ilo) = '*'

      else if ( bar(i,j) == INDEF ) then

      end if
!
!  Now fill in the interior of the cell.
!
      if ( i == istart .and. j == jstart ) then

        string(ilo+1:ilo+2) = '00'

      else if ( i== istop .and. j == jstop ) then

        string(ilo+1:ilo+2) = '$$'

      else if ( bar(i,j) == PATH ) then

        string(ilo+1:ilo+2) = '**'

      else if ( j <= n ) then

        if ( flat(i,j) == PATH .or. bar(i,j+1) == PATH .or. &
             flat(i+1,j) == PATH ) then

          string(ilo+1:ilo+2) = '**'

        end if

      end if

      ilo = ilo + 3

    end do

    ilo = ilo - 3

    write ( *, '(a)' ) string(1:ilo)

  end do

  string = ' '
  i = m+1
  ilo = 1
  do j = 1, nsafe

    if ( flat(i,j) == WALL ) then
      string(ilo:ilo+3) = '+--+'
    else if ( flat(i,j) == OPEN ) then
      string(ilo:ilo+3) = '+  +'
    else if ( flat(i,j) == PATH ) then
      string(ilo:ilo+3) = '+**+'
    else if ( flat(i,j) == INDEF ) then

    end if

    ilo = ilo + 3

  end do

  write ( *, '(a)' ) string(1:ilo)

  return
end
subroutine maze_random ( m, n, seed, bar, dad, flat )

!*****************************************************************************80
!
!! MAZE_RANDOM generates a random maze in a rectangular region.
!
!  Discussion:
!
!    The rectangular region is assumed to be made of a grid of M rows
!    and N columns of square cells.  The maze is to be begun in 
!    one cell, and ended in another.  The boundary of the region
!    is walled off, except possibly for entrances to the beginning
!    cell, and an exit from the ending cell.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns of cells.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
!    Output, integer ( kind = 4 ) BAR(M,N+1), records the vertical "bars" in the maze.
!    -1, means "indefinite", that there is no cell of the maze on either
!        side of this position;
!     0, means "wall", that there is a cell on at least one side, and
!        a wall here;
!     1, means "open", that there are cells on both sides (or possibly
!        an opening to the exterior) and the way is open.
!
!    Output, integer ( kind = 4 ) DAD(M,N), a rooted tree representation of
!    the maze.  The root of the tree has DAD(I,J) = 0.  All other cells
!    that are connectable to the root should have DAD(I,J) = K, where
!    K is the cell index K = ( J - 1 ) * M + I, with I and J the row
!    and column indices of the cell.  If the cell is not connectable
!    to the root, then DAD(I,J) is -1.
!
!    Output, integer ( kind = 4 ) FLAT(M+1,N), records the horizontal "flats" in the maze.
!    -1, means "indefinite", that there is no cell of the maze on either
!        side of this position;
!     0, means "wall", that there is a cell on at least one side, and
!        a wall here;
!     1, means "open", that there are cells on both sides (or possibly
!        an opening to the exterior) and the way is open.
!
  implicit none

  integer ( kind = 4 ), parameter :: maxstack = 500
  integer ( kind = 4 ), parameter :: NORTH = 1
  integer ( kind = 4 ), parameter :: EAST = 2
  integer ( kind = 4 ), parameter :: WEST = 3
  integer ( kind = 4 ), parameter :: SOUTH = 4
  integer ( kind = 4 ), parameter :: INDEF = -1
  integer ( kind = 4 ), parameter :: WALL = 0
  integer ( kind = 4 ), parameter :: OPEN = 1

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) bar(m,n+1)
  integer ( kind = 4 ) dad(m,n)
  integer ( kind = 4 ) dir
  integer ( kind = 4 ) flat(m+1,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jval
  integer ( kind = 4 ) nabe
  integer ( kind = 4 ) nbase
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) stack(maxstack)

  if ( m < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MAZE_RANDOM - Fatal error!'
    write ( *, '(a)' ) '  M must be at least 1.'
    stop
  end if

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MAZE_RANDOM - Fatal error!'
    write ( *, '(a)' ) '  N must be at least 1.'
    stop
  end if

  if ( m == 1 .and. n == 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MAZE_RANDOM - Fatal error!'
    write ( *, '(a)' ) '  At least one of M and N must be more than 1.'
    stop
  end if
!
!  Initialize arrays to INDEF.
!
  bar(1:m,1:n+1) = INDEF
  flat(1:m+1,1:n) = INDEF
!
!  Set the boundaries to walls.
!
  flat(1,1:n) = WALL
  flat(m+1,1:n) = WALL
  bar(1:m,1) = WALL
  bar(1:m,n+1) = WALL
!
!  Initialize the tree pointers.
!
  dad(1:m,1:n) = -1
!
!  Pick a random starting point.
!
  ival = i4_uniform ( 1, m, seed )
  jval = i4_uniform ( 1, n, seed )

  dad(ival,jval) = 0
!
!  Count the number of neighbors of the starting cell,
!  choose randomly from the neigbors, and add it.
!
  do

    nabe = 0

    do i = 1, m
      do j = 1, n

        if ( dad(i,j) /= -1 ) then

          if ( flat(i,j) /= WALL ) then

            if ( dad(i-1,j) == -1 ) then

              if ( maxstack < nabe + 3 ) then
                write ( *, '(a)' ) ' '
                write ( *, '(a)' ) 'MAZE_RANDOM - Fatal error!'
                write ( *, '(a)' ) '  Ran out of stack space.'
                return
              end if

              stack(nabe+1) = i
              stack(nabe+2) = j
              stack(nabe+3) = NORTH
              nabe = nabe + 3

            end if

          end if

          if ( bar(i,j+1) /= WALL ) then

            if ( dad(i,j+1) == -1 ) then

              if ( maxstack < nabe + 3 ) then
                write ( *, '(a)' ) ' '
                write ( *, '(a)' ) 'MAZE_RANDOM - Fatal error!'
                write ( *, '(a)' ) '  Ran out of stack space.'
                return
              end if

              stack(nabe+1) = i
              stack(nabe+2) = j
              stack(nabe+3) = EAST
              nabe = nabe + 3

            end if

          end if

          if ( bar(i,j) /= WALL ) then

            if ( dad(i,j-1) == -1 ) then

              if ( maxstack < nabe + 3 ) then
                write ( *, '(a)' ) ' '
                write ( *, '(a)' ) 'MAZE_RANDOM - Fatal error!'
                write ( *, '(a)' ) '  Ran out of stack space.'
                return
              end if

              stack(nabe+1) = i
              stack(nabe+2) = j
              stack(nabe+3) = WEST
              nabe = nabe + 3

            end if

          end if

          if ( flat(i+1,j) /= WALL ) then

            if ( dad(i+1,j) == -1 ) then

              if ( maxstack < nabe + 3 ) then
                write ( *, '(a)' ) ' '
                write ( *, '(a)' ) 'MAZE_RANDOM - Fatal error!'
                write ( *, '(a)' ) '  Ran out of stack space!'
                return
              end if

              stack(nabe+1) = i
              stack(nabe+2) = j
              stack(nabe+3) = SOUTH
              nabe = nabe + 3

            end if

          end if

        end if

      end do
    end do
!
!  If there are accessible neighbors, randomly choose one.
!
    if ( nabe <= 0 ) then
      exit
    end if

    ihi = nabe / 3
 
    ival = i4_uniform ( 1, ihi, seed )
  
    nbase = 3*ival - 3
    i = stack(nbase+1)
    j = stack(nbase+2)
    dir = stack(nbase+3)

    if ( dir == NORTH ) then
      flat(i,j) = OPEN
      dad(i-1,j) = ( j - 1 ) * m + i
    else if ( dir == EAST ) then
      bar(i,j+1) = OPEN
      dad(i,j+1) = ( j - 1 ) * m + i
    else if ( dir == WEST ) then
      bar(i,j) = OPEN
      dad(i,j-1) = ( j - 1 ) * m + i
    else if ( dir == SOUTH ) then
      flat(i+1,j) = OPEN
      dad(i+1,j) = ( j - 1 ) * m + i
    end if

  end do
!
!  Set all remaining INDEF's to WALLS.
!
  do i = 1, m
    do j = 1, n+1
      if ( bar(i,j) == INDEF ) then
        bar(i,j) = WALL
      end if
    end do
  end do

  do i = 1, m+1
    do j = 1, n
      if ( flat(i,j) == INDEF ) then
        flat(i,j) = WALL
      end if
    end do
  end do

  return
end
subroutine network_flow_max ( nnode, nedge, iendpt, icpflo, isorce, isink, &
  icut, node_flow )

!*****************************************************************************80
!
!! NETWORK_FLOW_MAX finds the maximal flow and a minimal cut in a network.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 July 2000
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
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) NEDGE, the number of edges.
!
!    Input/output, integer ( kind = 4 ) IENDPT(2,NEDGE), the edges of the network,
!    defined as pairs of nodes.  Each edge should be listed TWICE,
!    the second time in reverse order.  On output, the edges have
!    been reordered, and so the columns of IENDPT have been rearranged.
!
!    Input/output, integer ( kind = 4 ) ICPFLO(2,NEDGE).  Capacities and flows.
!    On input, ICPFLO(1,I) is the capacity of edge I.  On output,
!    ICPFLO(2,I) is the flow on edge I and ICPFLO(1,I) has
!    been rearranged to match the reordering of IENDPT.
!
!    Input, integer ( kind = 4 ) ISORCE, the designated source node.
!
!    Input, integer ( kind = 4 ) ISINK, the designated sink node.
!
!    Output, integer ( kind = 4 ) ICUT(NNODE).  ICUT(I) = 1 if node I is in the
!    minimal cut set, otherwise 0.
!
!    Output, integer ( kind = 4 ) NODE_FLOW(NNODE).  NODE_FLOW(I) is the value of the flow
!    through node I.
!
  implicit none

  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) i
  integer ( kind = 4 ) iarray(nnode)
  integer ( kind = 4 ) icpflo(2,nedge)
  integer ( kind = 4 ) icut(nnode)
  integer ( kind = 4 ) idel
  integer ( kind = 4 ) ien1
  integer ( kind = 4 ) ien2
  integer ( kind = 4 ) iendpt(2,nedge)
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) iparm
  integer ( kind = 4 ) iq
  integer ( kind = 4 ) ir
  integer ( kind = 4 ) iread
  integer ( kind = 4 ) irite
  integer ( kind = 4 ) is
  integer ( kind = 4 ) isink
  integer ( kind = 4 ) isorce
  integer ( kind = 4 ) isort
  integer ( kind = 4 ) it
  integer ( kind = 4 ) iwork(nnode,2)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) kz
  integer ( kind = 4 ) lst
  integer ( kind = 4 ) m
  integer ( kind = 4 ) node_flow(nnode)

  iarray(1:nnode) = 0
  idel = 0
 
  do i = 1, nedge

    icpflo(2,i) = 0
    ip = iendpt(1,i)

    if ( ip == isorce ) then
      idel = idel + icpflo(1,i)
    end if

    iarray(ip) = iarray(ip) + 1

  end do
 
  node_flow(isorce) = idel
  is = 1
 
  do i = 1, nnode
    it = iarray(i)
    iarray(i) = is
    iwork(i,1) = is
    is = is + it
  end do
 
  isort = 0
  ien1 = 0
  ien2 = 0
 
10    continue
 
  indx = 0
 
50    continue
 
  call sort_heap_external ( nedge, indx, ien1, ien2, is )

  if ( indx < 0 ) then
 
    is = iendpt(1,ien1) - iendpt(1,ien2)

    if ( is == 0 ) then
      is = iendpt(2,ien1) - iendpt(2,ien2)
    end if

  else if ( 0 < indx ) then
 
    do ir = 1, 2
      call i4_swap ( iendpt(ir,ien1), iendpt(ir,ien2) )
      call i4_swap ( icpflo(ir,ien1), icpflo(ir,ien2) )
    end do
 
  else
 
    if ( 0 < isort ) then
      return
    end if
 
    do i = 1, nedge
      iq = iendpt(2,i)
      iendpt(1,i) = iwork(iq,1)
      iwork(iq,1) = iwork(iq,1) + 1
    end do
 
    go to 100
 
  end if
 
  go to 50
 
80    continue
 
  iendpt(1,iendpt(1,ien1)) = ien2
  iendpt(1,iendpt(1,ien2)) = ien1
 
  do ir = 1, 2
    call i4_swap ( iendpt(ir,ien1), iendpt(ir,ien2) )
    call i4_swap ( icpflo(ir,ien1), icpflo(ir,ien2) )
  end do
 
  if ( indx < 0 ) then
    go to 270
  end if

  if ( indx == 0 ) then
    go to 170
  end if

  go to 50
 
100   continue
 
  indx = 0
 
  do i = 1, nnode

    if ( i /= isorce ) then
      node_flow(i) = 0
    end if

    iwork(i,2) = nedge + 1

    if ( i < nnode ) then
      iwork(i,2) = iarray(i+1)
    end if

    icut(i) = 0

  end do
 
  iread = 0
  irite = 1
  iwork(1,1) = isorce
  icut(isorce) = - 1
 
120   continue
 
  iread = iread + 1
 
  if ( iread <= irite ) then
 
    ip = iwork(iread,1)
    lst = iwork(ip,2) - 1
    i = iarray(ip) - 1
 
130     continue
 
    i = i + 1

    if ( lst < i ) then
      go to 120
    end if

    iq = iendpt(2,i)
    idel = icpflo(1,i) - icpflo(2,i)

    if ( icut(iq) /= 0 .or. idel == 0 ) then
      go to 130
    end if
 
    if ( iq /= isink ) then
      irite = irite + 1
      iwork(irite,1) = iq
    end if
 
    icut(iq) = - 1
    go to 130
 
  end if
 
  if ( icut(isink) == 0 ) then
 
    icut(1:nnode) = - icut(1:nnode)
 
    do i = 1, nedge
      ip = iendpt(2,iendpt(1,i))
      if ( icpflo(2,i) < 0 ) then
        node_flow(ip) = node_flow(ip) - icpflo(2,i)
      end if
      iendpt(1,i) = ip
    end do
 
    node_flow(isorce) = node_flow(isink)
    isort = 1
    go to 10
 
  end if
 
  icut(isink) = 1
 
160   continue
 
  iread = iread - 1

  if ( iread == 0 ) then
    go to 180
  end if

  ip = iwork(iread,1)
  ien1 = iarray(ip) - 1
  ien2 = iwork(ip,2) - 1
 
170   continue
 
  if ( ien1 /= ien2 ) then
 
    iq = iendpt(2,ien2)
 
    if ( icut(iq) <= 0 .or. icpflo(1,ien2) == icpflo(2,ien2) ) then
      ien2 = ien2 - 1
      go to 170
    end if
 
    iendpt(2,ien2) = - iq
    icpflo(1,ien2) = icpflo(1,ien2) - icpflo(2,ien2)
    icpflo(2,ien2) = 0
    ien1 = ien1 + 1

    if ( ien1 < ien2 ) then
      go to 80
    end if
 
  end if
 
  if ( iarray(ip) <= ien1 ) then
    icut(ip) = ien1
  end if

  go to 160
 
180   continue
 
  kz = 0
 
  do ir = 1, irite
    if ( 0 < icut(iwork(ir,1)) ) then
      kz = kz + 1
      iwork(kz,1) = iwork(ir,1)
    end if
  end do
 
  indx = - 1
  m = 1
 
200   continue
 
  ip = iwork(m,1)

  if ( 0 < node_flow(ip) ) then
    go to 250
  end if
 
210   continue
 
  m = m + 1

  if ( m <= kz ) then
    go to 200
  end if

  iparm = 0
 
220   continue
 
  m = m - 1
 
  if ( m == 1 ) then
 
    do i = 1, nedge
 
      iq = - iendpt(2,i)
 
      if ( 0 <= iq ) then

        iendpt(2,i) = iq
        j = iendpt(1,i)
        icpflo(1,i) = icpflo(1,i) - icpflo(2,j)

        idel = icpflo(2,i) - icpflo(2,j)
        icpflo(2,i) = idel
        icpflo(2,j) = - idel

      end if
 
    end do
 
    go to 100
 
  end if
 
  ip = iwork(m,1)

  if ( node_flow(ip) < 0 ) then
    go to 220
  end if
 
  if ( node_flow(ip) == 0 ) then

    lst = nedge + 1

    if ( ip < nnode ) then
      lst = iarray(ip+1)
    end if

    i = iwork(ip,2)
    iwork(ip,2) = lst
 
240     continue
 
    if ( i == lst ) then
      go to 220
    end if

    j = iendpt(1,i)
    idel = icpflo(2,j)
    icpflo(2,j) = 0
    icpflo(1,j) = icpflo(1,j) - idel
    icpflo(2,i) = icpflo(2,i) - idel
    i = i + 1
    go to 240
 
  end if
 
  if ( icut(ip) < iarray(ip) ) then
    go to 300
  end if
 
250   continue
 
  i = icut(ip) + 1
 
260   continue

  i = i - 1

  if ( i < iarray(ip) ) then
    go to 290
  end if

  iq = - iendpt(2,i)

  if ( node_flow(iq) < 0 ) then
    go to 260
  end if
 
  idel = icpflo(1,i) - icpflo(2,i)

  if ( node_flow(ip) < idel ) then
    idel = node_flow(ip)
  end if

  icpflo(2,i) = icpflo(2,i) + idel
  node_flow(ip) = node_flow(ip) - idel
  node_flow(iq) = node_flow(iq) + idel
  iparm = 1
  ien1 = iendpt(1,i)
  ien2 = iwork(iq,2) - 1

  if ( ien1 < ien2 ) then
    go to 80
  end if

  if ( ien1 /= ien2 ) then
    go to 280
  end if
 
270   continue
 
  iwork(iq,2) = ien2
 
280   continue
 
  if ( 0 < node_flow(ip) ) then
    go to 260
  end if

  if ( icpflo(1,i) == icpflo(2,i) ) then
    i = i - 1
  end if
 
290   continue
 
  icut(ip) = i

  if ( iparm /= 0 ) then
    go to 210
  end if
 
300   continue
 
  i = iwork(ip,2)
 
310   continue
 
  j = iendpt(1,i)
  idel = icpflo(2,j)

  if ( node_flow(ip) < idel ) then
    idel = node_flow(ip)
  end if

  icpflo(2,j) = icpflo(2,j) - idel
  node_flow(ip) = node_flow(ip) - idel
  iq = iendpt(2,i)
  node_flow(iq) = node_flow(iq) + idel
  i = i + 1

  if ( 0 < node_flow(ip) ) then
    go to 310
  end if

  node_flow(ip) = - 1
  go to 220
 
end
subroutine node_order_print ( nnode, order )

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
!    29 May 1998
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
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) order(nnode)

  inc = 15

  do ilo = 1, nnode, inc

    ihi = min ( ilo + inc - 1, nnode )

    write ( *, '(a)' ) ' '
    write ( *, '(a6,4x,15i4)' ) 'Order:', ( i, i = ilo, ihi )
    write ( *, '(a6,4x,15i4)' ) 'Label:', order(ilo:ihi)

  end do

  return
end
subroutine node_relax ( cor3, cor3_new, cor3_nabe, face, face_order, max_cor3, &
  max_face, max_order, num_cor3, num_face )

!*****************************************************************************80
!
!! NODE_RELAX smooths a shape by an averaging operation on the node positions.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) COR3(3,MAXCOR3), the coordinates of the nodes.
!
!    Output, real ( kind = 8 ) COR3_NEW(3,MAXCOR3), the new, averaged
!    coordinates of the nodes.
!
!    Workspace, integer COR3_NABE(MAXCOR3).  On output, COR3_NABE(I)
!    will contain the number of node neighbors of node I.
!
!    Input, integer ( kind = 4 ) FACE(MAX_ORDER,MAX_FACE), describes the faces.
!    FACE(I,J) is the index of the I-th node in face J. 
!
!    Input, integer ( kind = 4 ) FACE_ORDER(MAX_FACE), is the number of nodes
!    making up each face.
!
!    Input, integer ( kind = 4 ) MAX_FACE, the maximum number of faces.
!
!    Input, integer ( kind = 4 ) MAX_ORDER, is the maximum number of nodes that can
!    make up a face, required to dimension FACE.
!
!    Input, integer ( kind = 4 ) NUM_FACE, the number of faces.
!
  implicit none

  integer ( kind = 4 ) max_cor3
  integer ( kind = 4 ) max_face
  integer ( kind = 4 ) max_order

  real ( kind = 8 ) cor3(3,max_cor3)
  real ( kind = 8 ) cor3_new(3,max_cor3)
  integer ( kind = 4 ) cor3_nabe(max_cor3)
  integer ( kind = 4 ) face(max_order,max_face)
  integer ( kind = 4 ) face_order(max_face)
  integer ( kind = 4 ) icor3
  integer ( kind = 4 ) iface
  integer ( kind = 4 ) inode
  integer ( kind = 4 ) ivert
  integer ( kind = 4 ) jnode
  integer ( kind = 4 ) num_cor3
  integer ( kind = 4 ) num_face
!
!  COR3_NEW will contain the new averaged coordinates.
!
  cor3_nabe(1:num_cor3) = 0
  cor3_new(1:3,1:num_cor3) = 0.0D+00
!
!  Consider each edge.  Essentially, the edge (I,J) is a signal to
!  add the old coordinates of I to the new J coordinates, and vice versa.
!
!  Because we are using a face representation, many, perhaps all the
!  edges, will show up repeatedly, probably twice.  To keep the algorithm
!  simple, for now we will simply use an edge every time it shows up
!  in a face, which means that edges that occur in multiple faces
!  will be weighted more.
!
  do iface = 1, num_face

    inode = face(face_order(iface),iface)

    do ivert = 1, face_order(iface)
      jnode = inode
      inode = face(ivert,iface)
      cor3_nabe(inode) = cor3_nabe(inode) + 1
      cor3_nabe(jnode) = cor3_nabe(jnode) + 1
      cor3_new(1:3,jnode) = cor3_new(1:3,jnode) + cor3(1:3,inode)
      cor3_new(1:3,inode) = cor3_new(1:3,inode) + cor3(1:3,jnode)
    end do

  end do
!
!  Copy the new into the old.
!
  do icor3 = 1, num_cor3

    if ( cor3_nabe(icor3) /= 0 ) then
      cor3_new(1:3,icor3) = cor3_new(1:3,icor3) &
        / real ( cor3_nabe(icor3), kind = 8 )
    end if

  end do

  return
end
subroutine nodes_to_ps ( plotxmin2, plotymin2, alpha, iunit, nnode, x, y, &
  xmin, ymin )

!*****************************************************************************80
!
!! NODES_TO_PS writes subplot nodes to a PostScript file.
!
!  Discussion:
!
!    A small filled circle is placed at each node.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PLOTXMIN2, PLOTYMIN2, the Postscript point corresponding
!    to the physical point XMIN, YMIN.
!
!    Input, real ( kind = 8 ) ALPHA, the physical-to-Postscript scale factor.
!
!    Input, integer ( kind = 4 ) IUNIT, the output FORTRAN unit.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, real ( kind = 8 ) X(NNODE), Y(NNODE), the coordinates of points.
!
!    Input, real ( kind = 8 ) XMIN, YMIN, the coordinates of the physical
!    origin.
!
  implicit none

  integer ( kind = 4 ) nnode

  real ( kind = 8 ) alpha
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) plotxmin2
  integer ( kind = 4 ) plotymin2
  integer ( kind = 4 ) px1
  integer ( kind = 4 ) py1
  integer ( kind = 4 ), parameter :: rad = 10
  real ( kind = 8 ) x(nnode)
  real ( kind = 8 ) xmin
  real ( kind = 8 ) y(nnode)
  real ( kind = 8 ) ymin
!
!  Draw points.
!
  do i = 1, nnode

    px1 = plotxmin2 + nint ( alpha * ( x(i) - xmin ) )
    py1 = plotymin2 + nint ( alpha * ( y(i) - ymin ) )

    write ( iunit, '(3i4,a)' ) px1, py1, rad, ' 0 360 arc closepath fill'

  end do

  return
end
subroutine object_build ( face, face_object, face_order, face_rank, face_tier, &
  max_order, num_face, num_object )

!*****************************************************************************80
!
!! OBJECT_BUILD builds edge-connected "objects" out of polygonal faces.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) FACE(MAX_ORDER,NUM_FACE), describes the faces.
!    FACE(I,J) is the index of the I-th node in face J.  It is best
!    if the nodes of all faces are listed in counterclockwise order.
!
!    Output, integer ( kind = 4 ) FACE_OBJECT(NUM_FACE), describes the objects.
!    FACE_OBJECT(I) is the index of the edge-connected "object" to 
!    which face I belongs.
!
!    Input, integer ( kind = 4 ) FACE_ORDER(NUM_FACE), is the number of nodes
!    making up each face.
!
!    Output, integer ( kind = 4 ) FACE_RANK(NUM_FACE), is an ordered list of faces.
!    FACE_RANK(1) is the index of the face in the first tier of the 
!    first object, followed by second tier faces, and so on until
!    object one is complete.  Object two follows, and so on.
!
!    Output, integer ( kind = 4 ) FACE_TIER(NUM_FACE).  FACE_TIER(I) is the "tier"
!    of face I in its object.  The seed of the object has tier 1,
!    the neighbors of the seed have tier 2, and so on.
!
!    Input, integer ( kind = 4 ) MAX_ORDER, is the maximum number of nodes that can
!    make up a face, required to dimension FACE.
!
!    Input, integer ( kind = 4 ) NUM_FACE, the number of faces.
!
!    Output, integer ( kind = 4 ) NUM_OBJECT, the number of objects.
!
  implicit none

  integer ( kind = 4 ) max_order
  integer ( kind = 4 ) num_face

  integer ( kind = 4 ) face(max_order,num_face)
  integer ( kind = 4 ) face_object(num_face)
  integer ( kind = 4 ) face_order(num_face)
  integer ( kind = 4 ) face_rank(num_face)
  integer ( kind = 4 ) face_tier(num_face)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iface
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ihi_next
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) ilo_next
  integer ( kind = 4 ) irank
  integer ( kind = 4 ) jface
  integer ( kind = 4 ) num_object
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) tier
  integer ( kind = 4 ) touch
!
!  Initialization.
!
  num_object = 0

  if ( num_face <= 0 ) then
    return
  end if

  face_object(1:num_face) = 0
  face_rank(1:num_face) = 0
  face_tier(1:num_face) = 0

  irank = 0

  seed = 1
!
!  Begin the next object, seeded with face SEED.
!
10    continue

  tier = 1

  num_object = num_object + 1
  irank = irank + 1

  face_rank(irank) = seed
  face_tier(seed) = tier
  face_object(seed) = num_object

  ilo = irank
  ihi = irank
!
!  Begin the next tier of faces, which are neighbors of faces we
!  found in the previous tier.
!
20    continue

  tier = tier + 1

  ilo_next = ihi + 1
  ihi_next = ihi

  do jface = 1, num_face

    if ( face_tier(jface) == 0 ) then

      do i = ilo, ihi

        iface = face_rank(i)

        call face_touch ( face, face_order, max_order, num_face, iface, &
          jface, touch )

        if ( touch /= 0 ) then
          ihi_next = ihi_next + 1
          irank = irank + 1
          face_rank(irank) = jface
          face_tier(jface) = tier
          face_object(jface) = num_object
          exit
        end if

      end do

    end if

  end do

  if ( ilo_next <= ihi_next ) then
    ilo = ilo_next
    ihi = ihi_next
    go to 20
  end if
!
!  No neighbors were found, so this object is complete.  
!  Search for an unused face, which will be the seed of the next object.
!
  do iface = 1, num_face

    if ( face_tier(iface) == 0 ) then
      seed = iface
      go to 10
    end if

  end do

  return
end
subroutine perm_cycle ( isig, n, isgn, ncycle, iopt )

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
!      ISIG = 2, 3, 9, 6, 7, 8, 5, 4, 1
!
!    Output:
!
!      NCYCLE = 3
!      ISGN = +1
!      ISIG = -2, 3, 9, -6, -7, 8, 5, 4, 1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 July 2000
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
!    Input/output, integer ( kind = 4 ) ISIG(N).  On input, ISIG describes a
!    permutation, in the sense that entry I is to be moved to ISIG(I).
!    If IOPT = 0, then ISIG will not be changed by this routine.
!    If IOPT = 1, then on output, ISIG will be "tagged".  That is,
!    one element of every cycle in ISIG will be negated.  In this way,
!    a user can traverse a cycle by starting at any entry I1 of ISIG
!    which is negative, moving to I2 = ABS(ISIG(I1)), then to
!    ISIG(I2), and so on, until returning to I1.
!
!    Input, integer ( kind = 4 ) N, the number of objects being permuted.
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
  integer ( kind = 4 ) isig(n)
  integer ( kind = 4 ) ncycle

  is = 1
  ncycle = n

  do i = 1, n

    i1 = isig(i)

    do while ( i < i1 )
      ncycle = ncycle - 1
      i2 = isig(i1)
      isig(i1) = - i2
      i1 = i2
    end do

    if ( iopt /= 0 ) then
      is = - isign ( 1, isig(i) )
    end if

    isig(i) = isign ( isig(i), is )

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
!    02 November 1999
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
  integer ( kind = 4 ) n

  n = npart + nfree

  if ( npart < 0 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PERM_FREE - Fatal error!'
    write ( *, '(a)' ) '  NPART < 0.'
    stop

  else if ( npart == 0 ) then

    call i4vec_indicator ( n, ifree )

  else if ( nfree < 0 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PERM_FREE - Fatal error!'
    write ( *, '(a)' ) '  NFREE < 0.'
    stop

  else if ( nfree == 0 ) then

    return

  else

    k = 0

    do i = 1, n

      do j = 1, npart
        if ( ipart(j) == i ) then
          go to 10
        end if
      end do

      k = k + 1

      if ( nfree < k ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PERM_FREE - Fatal error!'
        write ( *, '(a)' ) '  The partial permutation is illegal.'
        write ( *, '(a)' ) '  It should contain, at most once, some of'
        write ( *, '(a,i8)' ) '  the integers between 1 and ', n
        stop
      end if

      ifree(k) = i

10    continue

    end do

  end if

  return
end
subroutine perm_inc ( iperm, ipos, n )

!*****************************************************************************80
!
!! PERM_INC "increments" a permutation to get the "next" one.
!
!  Discussion:
!
!    The routine is given IPERM, a permutation of the numbers from 1 to N,
!    and a position IPOS between 1 and N.
!
!    It returns the next permutation in the dictionary order which
!    comes after all permutations beginning IPERM(1) through IPERM(IPOS).
!
!  Example:
!
!             PERM              IPOS
!
!    Input    123456789         7
!    Output   123456798         7
!
!    Input    123456789         9
!    Output   213456789         0
!
!    Input    134826795         3
!    Output   134925678         3
!
!    Input    134826795         0
!    Output   123456789         0
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
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) IPERM(N).
!    On input, the current permutation.
!    On output, the "incremented" permutation.
!
!    Input/output, integer ( kind = 4 ) IPOS.
!    On input, IPOS is the location of the end of the string of
!    "digits" in IPERM that form the test string.  That is, the
!    new permutation to be computed must be the very next one,
!    in dictionary order, which succeeds all strings whose first
!    IPOS digits agree with the input value of IPERM.
!
!    On output, IPOS is the position of the last digit of the output
!    value of IPERM which agrees with the input value of IPERM.
!
!    Input, integer ( kind = 4 ) N, is the number of entries in IPERM.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) ipcopy
  integer ( kind = 4 ) iperm(n)
  integer ( kind = 4 ) ipos
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) new

  if ( ipos == 0 ) then
    ipos = n
    call i4vec_indicator ( n, iperm )
    return
  end if
 
  ipcopy = ipos

10    continue
!
!  To get the next permutation, we need to increment the IPOS+1 "digit".
!
!  We do this by finding, if possible, a digit in positions IPOS+2
!  through N that is just larger than the current value IPOS+1 digit.
!  If we find such a digit, it becomes the IPOS+1 digit, and the
!  remaining values are sorted into increasing order.
!
  new = 0
  do j = ipcopy+2, n
    if ( new == 0 ) then
      if ( iperm(ipcopy+1) < iperm(j) ) then
        new = j
      end if
    else
      if ( iperm(ipcopy+1) < iperm(j) .and. iperm(j) < iperm(new) ) then
        new = j
      end if
    end if
  end do
!
!  There is a next candidate that agrees with IPERM through entry I.
!  Swap entries IPOS+1 and NEW, and sort the entries (IPOS+2,...,N).
!
!  The output value of IPOS equals the input value.
!
  if ( new /= 0 ) then

    call i4_swap ( iperm(new), iperm(ipcopy+1) )
 
    do j = ipcopy+2, n
 
      do k = j+1, n
        if ( iperm(k) < iperm(j) ) then
          call i4_swap ( iperm(j), iperm(k) )
        end if
      end do
 
    end do
    return
  end if
!
!  There is no next candidate that agrees with IPERM through entry 
!  IPOS.  Can we decrease IPOS and try for a next candidate that way?
!
  if ( 0 < ipcopy ) then
    ipcopy = ipcopy - 1
    go to 10
  end if
!
!  IPOS is now zero.  There is no successor to the current permutation,
!  so we start again at the first permutation.
!
  ipos = 0
  call i4vec_indicator ( n, iperm )
 
  return
end
subroutine perm_inv ( n, isig )

!*****************************************************************************80
!
!! PERM_INV inverts a permutation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 July 2000
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of objects being permuted.
!
!    Input/output, integer ( kind = 4 ) ISIG(N).
!
!    On input, ISIG describes a permutation.
!
!    ISIG is used to represent a permutation by the convention that
!    the permutation maps the letter I to ISIG(I).  Thus, if ISIG
!    contains the values (4, 1, 3, 2), then the permutation
!    represented permutes 1 to 4, 2 to 1, 3 to 3, and 4 to 2.
!
!    On output, ISIG describes the inverse permutation
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i0
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) is
  integer ( kind = 4 ) isig(n)

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PERM_INV - Fatal error!'
    write ( *, '(a,i8)' ) '  Input value of N = ', n
    stop
  end if

  is = 1

  do i = 1, n

    i1 = isig(i)

    do while ( i < i1 )
      i2 = isig(i1)
      isig(i1) = - i2
      i1 = i2
    end do

    is = - isign ( 1, isig(i) )
    isig(i) = isign ( isig(i), is )

  end do

  do i = 1, n

    i1 = - isig(i)

    if ( 0 <= i1 ) then

      i0 = i

      do

        i2 = isig(i1)
        isig(i1) = i0

        if ( i2 < 0 ) then
          exit
        end if

        i0 = i1
        i1 = i2

      end do

    end if

  end do

  return
end
subroutine perm_next ( n, iarray, more, even )

!*****************************************************************************80
!
!! PERM_NEXT computes all of the permutations on N objects, one at a time.
!
!  Discussion:
!
!    If the routine is called with MORE = .TRUE., any permutation in
!    IARRAY, and EVEN = .TRUE., then the successor of the input
!    permutation will be produced, unless IARRAY is the last permutation
!    on N letters, in which case IARRAY(1) will be set to 0 on return.
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
!    Input, integer ( kind = 4 ) N, the number of objects being permuted.
!
!    Input/output, integer ( kind = 4 ) IARRAY(N).
!
!    If MORE is .TRUE., then IARRAY is assumed to contain the
!    "previous" permutation, and on IARRAY(I) is the value
!    of the I-th object under the next permutation.
!
!    Otherwise, IARRAY(I) will be set to the "first" permutation.
!
!    Input/output, logical MORE.
!
!    Set MORE to FALSE before first calling this routine.
!
!    MORE will be reset to .TRUE. and a permutation will be returned.
!
!    Each new call produces a new permutation until
!    MORE is returned .FALSE.
!
!    Output, logical EVEN.
!
!    EVEN is .TRUE. if the output permutation is even, that is,
!    involves an even number of transpositions.
!
!    EVEN is .FALSE. otherwise.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) ia
  integer ( kind = 4 ) iarray(n)
  integer ( kind = 4 ) id
  integer ( kind = 4 ) is
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  logical more
  logical even

  if ( .not. more ) then

    call i4vec_indicator ( n, iarray )

    more = .true.
    even = .true.

    if ( n == 1 ) then
      more = .false.
      return
    end if

    if ( iarray(n) /= 1 .or. iarray(1) /= 2 + mod ( n, 2 ) ) then
      return
    end if

    do i = 1, n-3
      if ( iarray(i+1) /= iarray(i)+1 ) then
        return
      end if
    end do

    more = .false.

  else

    if ( n == 1 ) then
      iarray(1) = 0
      more = .false.
      return
    end if

    if ( even ) then

      ia = iarray(1)
      iarray(1) = iarray(2)
      iarray(2) = ia
      even = .false.

      if ( iarray(n) /= 1 .or. iarray(1) /= 2 + mod ( n, 2 ) ) then
        return
      end if

      do i = 1, n-3
        if ( iarray(i+1) /= iarray(i)+1 ) then
          return
        end if
      end do

      more = .false.
      return

    else

      is = 0

      do i1 = 2, n

        ia = iarray(i1)
        i = i1-1
        id = 0

        do j = 1, i
          if ( ia < iarray(j) ) then
            id = id + 1
          end if
        end do

        is = id + is

        if ( id /= i * mod ( is, 2 ) ) then
          go to 10
        end if

      end do

      iarray(1) = 0
      more = .false.
      return

    end if

10      continue

    m = mod ( is+1, 2 ) * (n+1)

    do j = 1, i

      if ( isign(1,iarray(j)-ia) /= isign(1,iarray(j)-m) ) then
        m = iarray(j)
        l = j
      end if

    end do

    iarray(l) = ia
    iarray(i1) = m
    even = .true.

  end if

  return
end
subroutine perm_random ( n, seed, iarray )

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
!    02 March 1999
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
!    Input, integer ( kind = 4 ) N, the number of objects to be permuted.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
!    Output, integer ( kind = 4 ) IARRAY(N), the random permutation.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) iarray(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) seed

  call i4vec_indicator ( n, iarray )

  do i = 1, n
    j = i4_uniform ( i, n, seed )
    call i4_swap ( iarray(j), iarray(i) )
  end do

  return
end
subroutine poly ( n, iarray, ix0, iopt, ival )

!*****************************************************************************80
!
!! POLY performs operations on polynomials in power or factorial form.
!
!  Definition:
!
!    The power sum form of a polynomial is
!
!      P(X) = A1+A2*X+A3*X**2+...+(AN+1)*X**N
!
!    The Taylor expansion at C has the form
!
!      P(X) = A1+A2*(X-C)+A3*(X-C)**2+...+(AN+1)*(X-C)**N
!
!    The factorial form of a polynomial is
!
!      P(X) = A1+A2*X+A3*(X)*(X-1)+A4*(X)*(X-1)*(X-2)+...
!        +(AN+1)*(X)*(X-1)*...*(X-N+1)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 April 1999
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of coefficients in the polynomial
!    (in other words, the polynomial degree + 1)
!
!    Input, integer ( kind = 4 ) IOPT, a flag describing which algorithm is to
!    be carried out:
!
!    -3: Reverse Stirling.  Input the coefficients of
!    the polynomial in factorial form, output them in
!    power sum form.
!
!    -2: Stirling.  Input the coefficients in power sum
!    form, output them in factorial form.
!
!    -1: Evaluate a polynomial which has been input
!    in factorial form.
!
!    0:  Evaluate a polynomial input in power sum form.
!
!    1 or more:  Given the coefficients of a polynomial in
!    power sum form, compute the first IOPT coefficients of
!    the polynomial in Taylor expansion form.
!
!    Input, integer ( kind = 4 ) IX0, for IOPT = -1, 0, or positive, the value X of the
!    argument at which the polynomial is to be evaluated, or the
!    Taylor expansion is to be carried out.
!
!    Output, integer ( kind = 4 ) IVAL, for IOPT = -1 or 0, the value of the
!    polynomial at the point IX0.
!
!    Input, integer ( kind = 4 ) IARRAY(N).  Contains the coefficients of the
!    polynomial.  Depending on the option chosen, these coefficients may
!    be overwritten by those of a different form of the polynomial.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) iarray(n)
  integer ( kind = 4 ) ieps
  integer ( kind = 4 ) iopt
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) iw
  integer ( kind = 4 ) ix0
  integer ( kind = 4 ) iz
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n1

  n1 = min ( n, iopt )
  n1 = max ( 1, n1 )
 
  if ( iopt < -1 ) then
    n1 = n
  end if
 
  ieps = mod ( max ( -iopt, 0 ), 2 )
 
  iw = -n * ieps

  if ( -2 < iopt ) then
    iw = iw + ix0
  end if
 
  do m = 1, n1
 
    ival = 0
    iz = iw
 
    do i = m, n
      iz = iz + ieps
      ival = iarray(n+m-i) + iz * ival
      if ( iopt /= 0 .and. iopt /= -1 ) then
        iarray(n+m-i) = ival
      end if
    end do
 
    if ( iopt < 0 ) then
      iw = iw + 1
    end if
 
  end do
 
  return
end
subroutine poly_to_tri ( face, ierror, max_face, max_vert, num_face, num_vert )

!*****************************************************************************80
!
!! POLY_TO_TRI converts a collection of polygons into a collection of triangles.
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
!    Input/output, integer ( kind = 4 ) FACE(MAX_VERT,MAX_FACE), describes 
!    the faces.  FACE(I,J) is the I-th node associated with the J-th face.  
!    This array is updated on return.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error.
!    1, the algorithm failed because MAX_FACE was too small.
!    2, the algorithm failed because there were faces of order < 3.
!    3, the algorithm failed because there were faces of order > MAX_VERT.
!
!    Input, integer ( kind = 4 ) MAX_FACE, the maximum number of faces allowed.
!
!    Input, integer ( kind = 4 ) MAX_VERT, the maximum number of nodes allowed 
!    per face.
!
!    Input/output, integer ( kind = 4 ) NUM_FACE, the number of faces.  
!    This value is updated on return.
!
!    Input/output, integer ( kind = 4 ) NUM_VERT(MAX_FACE), the number of nodes
!    associated with each face.  On successful return, every entry of
!    this array will be 3.
!
  implicit none

  integer ( kind = 4 ) max_face
  integer ( kind = 4 ) max_vert

  integer ( kind = 4 ) face(max_vert,max_face)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iface
  integer ( kind = 4 ) iface_old
  integer ( kind = 4 ) ivert
  integer ( kind = 4 ) k
  integer ( kind = 4 ) num_face
  integer ( kind = 4 ) num_face2
  integer ( kind = 4 ) num_vert(max_face)

  ierror = 0
  num_face2 = 0

  do iface = 1, num_face

    if ( num_vert(iface) < 3 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'POLY_TO_TRI - Fatal error!'
      write ( *, '(a,i8)' ) '  Illegal face ', iface
      write ( *, '(a,i8)' ) '  Number of nodes is ', num_vert(iface)
      ierror = 2
      return
    else if ( max_vert < num_vert(iface) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'POLY_TO_TRI - Fatal error!'
      write ( *, '(a,i8)' ) '  Illegal face ', iface
      write ( *, '(a,i8)' ) '  Number of nodes is ', num_vert(iface)
      write ( *, '(a,i8)' ) '  MAX_VERT is ', max_vert
      ierror = 3
      return
    end if

    do ivert = 3, num_vert(iface)
      num_face2 = num_face2 + 1
    end do

  end do

  if ( max_face < num_face2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'POLY_TO_TRI - Fatal error!'
    write ( *, '(a)' ) '  MAX_FACE is too small to replace all faces'
    write ( *, '(a)' ) '  by triangles.'
    write ( *, '(a,i8)' ) '  MAX_FACE = ', max_face
    write ( *, '(a,i8)' ) '  NUM_FACE2 = ', num_face2
    ierror = 1
    return
  end if

  iface_old = num_face
  k = num_vert(iface_old)

  do iface = num_face2, 1, -1

    if ( k < 3 ) then
      iface_old = iface_old - 1
      k = num_vert(iface_old)
    end if

    num_vert(iface) = 3
    face(1,iface) = face(1,iface_old)
    do ivert = 2, 3
      face(ivert,iface) = face(k+ivert-3,iface_old)
    end do

    k = k - 1

  end do

  num_face = num_face2

  return
end
subroutine pruefer_to_tree_arc ( nnode, iarray, inode, jnode )

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
!  Reference:
!
!    Dennis Stanton, Dennis White,
!    Constructive Combinatorics,
!    Springer Verlag, New York, 1986.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) IARRAY(NNODE-2), the Pruefer code of the tree.
!
!    Output, integer ( kind = 4 ) INODE(NNODE-1), JNODE(NNODE-1), the edge array 
!    of the tree.  The I-th edge joins nodes INODE(I) and JNODE(I).
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) i
  integer ( kind = 4 ) iarray(nnode-2)
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
    iwork(iarray(i)) = iwork(iarray(i)) + 1
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
    jnode(i) = iarray(i)
    iwork(ii) = 0
    iwork(iarray(i)) = iwork(iarray(i)) - 1
 
  end do
 
  inode(nnode-1) = iarray(nnode-2)
 
  if ( iarray(nnode-2) /= 1 ) then
    jnode(nnode-1) = 1
  else
    jnode(nnode-1) = 2
  end if
 
  return
end
subroutine pruefer_to_tree_2 ( nnode, iarray, itree )

!*****************************************************************************80
!
!! PRUEFER_TO_TREE_2 produces the edge list of a tree from its Pruefer code.
!
!  Discussion:
!
!    One can thus exhibit all trees on N nodes, produce
!    one at random, find the M-th one on the list, etc, by
!    manipulating the Pruefer codes.
!
!    For every labeled tree on N nodes, there is a unique N-2 tuple
!    of integers A1 through AN-2, with each A between 1 and N.  There
!    are N**(N-2) such sequences, and each one is associated with exactly
!    one tree.
!
!    This routine apparently assumes that the Pruefer code is
!    generated by taking the LOWEST labeled terminal node each time.
!    This is not consistent with PRUEFER_TO_TREE and TREE_TO_PRUEFER.
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
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Albert Nijenhuis. Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NNODE, number of nodes in desired tree.
!
!    Output, integer ( kind = 4 ) IARRAY(NNODE).  IARRAY(I), I = 1, NNODE-2 is the Pruefer
!    code for the tree.
!
!    Output, integer ( kind = 4 ) ITREE(NNODE); the I-th edge of the tree joins nodes
!    I and ITREE(I).
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ir
  integer ( kind = 4 ) itree(nnode)
  integer ( kind = 4 ) iarray(nnode)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kp
  integer ( kind = 4 ) l

  itree(1:nnode) = 0
 
  do i = nnode-2, 1, -1
 
    l = iarray(i)
 
    if ( itree(l) == 0 ) then
      iarray(i) = - l
      itree(l) = - 1
    end if
 
  end do
 
  iarray(nnode-1) = nnode
!
!  Find next index K so that ITREE(K) is 0.
!
  k = 1
  j = 0
 
10 continue

  do while ( itree(k) /= 0 )
    k = k + 1
  end do
 
  kp = k
 
20    continue
 
  j = j + 1
  ir = abs ( iarray(j) )
  itree(kp) = ir
 
  if ( j /= nnode-1 ) then
 
    if ( 0 < iarray(j) ) then
      go to 10
    end if
 
    if ( k < ir ) then
      itree(ir) = 0
      go to 10
    end if
 
    kp = ir
    go to 20
 
  end if
!
!  Restore the signs of IARRAY.
!
  iarray(1:nnode-2) = abs ( iarray(1:nnode-2) )
 
  return
end
function pythag ( a, b )

!*****************************************************************************80
!
!! PYTHAG computes SQRT ( A**2 + B**2 ) carefully.
!
!  Discussion:
!
!    The formula
!
!      PYTHAG = sqrt ( A**2 + B**2 )
!
!    is reasonably accurate, but the formula can actually fail if
!    for example, A**2 is larger than the machine overflow.  The
!    formula can lose most of its accuracy if the sum of the squares
!    is very large or very small.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 March 2000
!
!  Reference:
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow, 
!    Y Ikebe, V Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the two legs of a right triangle.
!
!    Output, real ( kind = 8 ) PYTHAG, the length of the hypotenuse.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) p
  real ( kind = 8 ) pythag
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) t
  real ( kind = 8 ) u

  p = max ( abs ( a ), abs ( b ) )

  if ( p /= 0.0D+00 ) then

    r = ( min ( abs ( a ), abs ( b ) ) / p )**2

   10   continue

    t = 4.0D+00 + r

    if ( t /= 4.0D+00 ) then
      s = r / t
      u = 1.0D+00 + 2.0D+00 * s
      p = u * p
      r = ( s / u )**2 * r
      go to 10
    end if

  end if

  pythag = p

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
  real ( kind = 4 ) r4_uniform_01

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
subroutine r8_swap ( x, y )

!*****************************************************************************80
!
!! R8_SWAP swaps two double precision values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
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
function r8_uniform_01 ( seed )

!*****************************************************************************80
!
!! R8_UNIFORM_01 returns a unit pseudorandom R8.
!
!  Discussion:
!
!    An R8 is a real ( kind = 4 ) value.
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
subroutine r8col_find ( lda, m, n, a, x, i4col )

!*****************************************************************************80
!
!! R8COL_FIND seeks a table column equal to a real vector.
!
!  Example:
!
!    Input:
!
!      M = 3,
!      N = 4,
!
!      A = (
!        1.  2.  3.  4.
!        5.  6.  7.  8.
!        9. 10. 11. 12. )
!
!      x = ( 3.,
!            7.,
!           11. )
!
!    Output:
!
!      I4COL = 3
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
!    Input, integer ( kind = 4 ) LDA, the leading dimension of the array,
!    which must be at least M.
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(LDA,N), a table of numbers, regarded as
!    N columns of vectors of length M.
!
!    Input, real ( kind = 8 ) X(M), a vector to be matched with a column of A.
!
!    Output, integer ( kind = 4 ) I4COL, the index of the first column of A
!    which exactly matches every entry of X, or 0 if no match
!    could be found.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(lda,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4col
  integer ( kind = 4 ) j
  real ( kind = 8 ) x(m)

  i4col = 0

  do j = 1, n

    i4col = j

    do i = 1, m
      if ( x(i) /= a(i,j) ) then
        i4col = 0
        exit
      end if
    end do

    if ( i4col /= 0 ) then
      return
    end if

  end do

  return
end
subroutine r8mat_print ( a, ihi, ilo, jhi, jlo, lda, ncol, nrow )

!*****************************************************************************80
!
!! R8MAT_PRINT prints out a portion of a dense matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
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
  character ctemp(incx)*14
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

      write ( *, '(i5,1x,5a14)' ) i, ctemp(1:inc)

    end do

  end do

  write ( *, '(a)' ) ' '

  return
end
subroutine r8vec_print ( n, a, title )

!*****************************************************************************80
!
!! R8VEC_PRINT prints an R8VEC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
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
subroutine r8vec_uniform ( n, a, b, seed, r )

!*****************************************************************************80
!
!! R8VEC_UNIFORM returns a scaled pseudorandom R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    For now, the input quantity SEED is an integer variable.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Second Edition,
!    Springer, 1987,
!    ISBN: 0387964673,
!    LC: QA76.9.C65.B73.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, December 1986, pages 362-376.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998,
!    ISBN: 0471134031,
!    LC: T57.62.H37.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, Number 2, 1969, pages 136-143.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of entries in the vector.
!
!    Input, real ( kind = 8 ) A, B, the lower and upper limits.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which 
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R(N), the vector of pseudorandom values.
!
  implicit none

  integer ( kind = 4 ) n

  real       ( kind = 8 ) a
  real       ( kind = 8 ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real       ( kind = 8 ) r(n)

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_UNIFORM - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + i4_huge
    end if

    r(i) = a + ( b - a ) * real ( seed, kind = 8 ) * 4.656612875D-10

  end do

  return
end
subroutine r8vec_uniform_01 ( n, seed, r )

!*****************************************************************************80
!
!! R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    For now, the input quantity SEED is an integer variable.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Second Edition,
!    Springer, 1987,
!    ISBN: 0387964673,
!    LC: QA76.9.C65.B73.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, December 1986, pages 362-376.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998,
!    ISBN: 0471134031,
!    LC: T57.62.H37.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, Number 2, 1969, pages 136-143.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which 
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R(N), the vector of pseudorandom values.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real       ( kind = 8 ) r(n)

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + i4_huge
    end if

    r(i) = real ( seed, kind = 8 ) * 4.656612875D-10

  end do

  return
end
subroutine r8vec2_print ( n, a1, a2, title )

!*****************************************************************************80
!
!! R8VEC2_PRINT prints a pair of R8VEC's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, real ( kind = 8 ) A1(N), A2(N), the vectors to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a1(n)
  real ( kind = 8 ) a2(n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  if ( title /= ' ' ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(i8,2g14.6)' ) i, a1(i), a2(i)
  end do

  return
end
function r8vec3_compare ( x1, y1, z1, x2, y2, z2 )

!*****************************************************************************80
!
!! R8VEC3_COMPARE compares two R8VEC's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X1, Y1, Z1, the first vector.
!
!    Input, real ( kind = 8 ) X2, Y2, Z2, the second vector.
!
!    Output, character R8VEC3_COMPARE: '<', '>' or '=' if the first vector
!    is less, greater or equal to the second.
!
  implicit none

  character c
  character r8vec3_compare
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) z1
  real ( kind = 8 ) z2

  if ( x1 < x2 ) then
    c = '<'
  else if ( x1 > x2 ) then
    c = '>'
  else if ( y1 < y2 ) then
    c = '<'
  else if ( y1 > y2 ) then
    c = '>'
  else if ( z1 < z2 ) then
    c = '<'
  else if ( z1 > z2 ) then
    c = '>'
  else
    c = '='
  end if

  r8vec3_compare = c

  return
end
subroutine r8vec3_index_insert_unique ( maxn, n, x, y, z, indx, &
  xval, yval, zval, ival, ierror )

!*****************************************************************************80
!
!! R8VEC3_INDEX_INSERT_UNIQUE inserts a unique D3 value in an indexed sorted list.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) MAXN, the maximum size of the list.
!
!    Input/output, integer ( kind = 4 ) N, the size of the list.
!
!    Input/output, real ( kind = 8 ) X(N), Y(N), Z(N), the list of R3 vectors.
!
!    Input, integer ( kind = 4 ) INDX(N), the sort index of the list.
!
!    Input, real ( kind = 8 ) XVAL, YVAL, ZVAL, the value to be inserted 
!    if it is not already in the list.
!
!    Output, integer ( kind = 4 ) IVAL, the index in INDX corresponding to the
!    value XVAL, YVAL, ZVAL.
!
!    Output, integer ( kind = 4 ) IERROR, 0 for no error, 1 if an error occurred.
!
  implicit none

  integer ( kind = 4 ) maxn

  integer ( kind = 4 ) equal
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) indx(maxn)
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) less
  integer ( kind = 4 ) more
  integer ( kind = 4 ) n
  real ( kind = 8 ) x(maxn)
  real ( kind = 8 ) xval
  real ( kind = 8 ) y(maxn)
  real ( kind = 8 ) yval
  real ( kind = 8 ) z(maxn)
  real ( kind = 8 ) zval

  ierror = 0

  if ( n <= 0 ) then

    if ( maxn <= 0 ) then
      ierror = 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8VEC3_INDEX_INSERT_UNIQUE - Fatal error!'
      write ( *, '(a)' ) '  Not enough space to store new data.'
      return
    end if

    n = 1
    x(1) = xval
    y(1) = yval
    z(1) = zval
    indx(1) = 1
    ival = 1
    return

  end if
!
!  Does ( XVAL, YVAL, ZVAL ) already occur in ( X, Y, Z)?
!
  call r8vec3_index_search ( maxn, n, x, y, z, indx, xval, yval, zval, &
    less, equal, more )

  if ( equal == 0 ) then

    if ( maxn <= n ) then
      ierror = 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8VEC3_INDEX_INSERT_UNIQUE - Fatal error!'
      write ( *, '(a)' ) '  Not enough space to store new data.'
      return
    end if

    x(n+1) = xval
    y(n+1) = yval
    z(n+1) = zval
    ival = more
    indx(n+1:more+1:-1) = indx(n:more:-1)
    indx(more) = n + 1
    n = n + 1

  else

    ival = equal

  end if

  return
end
subroutine r8vec3_index_search ( maxn, n, x, y, z, indx, xval, yval, &
  zval, less, equal, more )

!*****************************************************************************80
!
!! R8VEC3_INDEX_SEARCH searches for an R3 value in an indexed sorted list.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) MAXN, the maximum size of the list.
!
!    Input, integer ( kind = 4 ) N, the size of the current list.
!
!    Input, real ( kind = 8 ) X(N), Y(N), Z(N), the list.
!
!    Input, integer ( kind = 4 ) INDX(N), the sort index of the list.
!
!    Input, real ( kind = 8 ) XVAL, YVAL, ZVAL, the value to be sought.
!
!    Output, integer ( kind = 4 ) LESS, EQUAL, MORE, the indexes in INDX of the
!    entries of X that are just less than, equal to, and just greater
!    than XVAL.  If XVAL does not occur in X, then EQUAL is zero.
!    If XVAL is the minimum entry of X, then LESS is 0.  If XVAL
!    is the greatest entry of X, then MORE is N+1.
!
  implicit none

  integer ( kind = 4 ) maxn

  character            c
  integer ( kind = 4 ) equal
  integer ( kind = 4 ) hi
  integer ( kind = 4 ) indx(maxn)
  integer ( kind = 4 ) less
  integer ( kind = 4 ) lo
  integer ( kind = 4 ) mid
  integer ( kind = 4 ) more
  integer ( kind = 4 ) n
  character            r8vec3_compare
  real       ( kind = 8 ) x(maxn)
  real       ( kind = 8 ) xhi
  real       ( kind = 8 ) xlo
  real       ( kind = 8 ) xmid
  real       ( kind = 8 ) xval
  real       ( kind = 8 ) y(maxn)
  real       ( kind = 8 ) yhi
  real       ( kind = 8 ) ylo
  real       ( kind = 8 ) ymid
  real       ( kind = 8 ) yval
  real       ( kind = 8 ) z(maxn)
  real       ( kind = 8 ) zhi
  real       ( kind = 8 ) zlo
  real       ( kind = 8 ) zmid
  real       ( kind = 8 ) zval

  if ( n <= 0 ) then
    less = 0
    equal = 0
    more = 0
    return
  end if

  lo = 1
  hi = n

  xlo = x(indx(lo))
  ylo = y(indx(lo))
  zlo = z(indx(lo))

  xhi = x(indx(hi))
  yhi = y(indx(hi))
  zhi = z(indx(hi))

  c = r8vec3_compare ( xval, yval, zval, xlo, ylo, zlo )

  if ( c == '<' ) then
    less = 0
    equal = 0
    more = 1
    return
  else if ( c == '=' ) then
    less = 0
    equal = 1
    more = 2
    return
  end if

  c = r8vec3_compare ( xval, yval, zval, xhi, yhi, zhi )

  if ( c == '>' ) then
    less = n
    equal = 0
    more = n + 1
    return
  else if ( c == '=' ) then
    less = n - 1
    equal = n
    more = n + 1
    return
  end if

  do

    if ( lo + 1 == hi ) then
      less = lo
      equal = 0
      more = hi
      return
    end if

    mid = ( lo + hi ) / 2
    xmid = x(indx(mid))
    ymid = y(indx(mid))
    zmid = z(indx(mid))

    c = r8vec3_compare ( xval, yval, zval, xmid, ymid, zmid )

    if ( c == '=' ) then
      equal = mid
      less = equal - 1
      more = equal + 1
      return
    else if ( c == '<' ) then
      hi = mid
    else if ( c == '>' ) then
      lo = mid
    end if

  end do

  return
end
subroutine s_blanks_delete ( s )

!*****************************************************************************80
!
!! S_BLANKS_DELETE replaces consecutive blanks by one blank.
!
!  Discussion:
!
!    The remaining characters are left justified and right padded with blanks.
!    TAB characters are converted to spaces.
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

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  character newchr
  character oldchr
  character ( len = * ) s
  character, parameter :: TAB = char ( 9 )

  j = 0
  newchr = ' '

  do i = 1, len ( s )

    oldchr = newchr
    newchr = s(i:i)

    if ( newchr == TAB ) then
      newchr = ' '
    end if

    s(i:i) = ' '

    if ( oldchr /= ' ' .or. newchr /= ' ' ) then
      j = j + 1
      s(j:j) = newchr
    end if

  end do

  return
end
subroutine s_cat ( s1, s2, s3 )

!*****************************************************************************80
!
!! S_CAT concatenates two strings to make a third string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S1, the "prefix" string.
!
!    Input, character ( len = * ) S2, the "postfix" string.
!
!    Output, character ( len = * ) S3, the string made by
!    concatenating S1 and S2, ignoring any trailing blanks.
!
  implicit none

  character ( len = * ) s1
  character ( len = * ) s2
  character ( len = * ) s3

  s3 = trim ( s1 ) // trim ( s2 )

  return
end
function s_eqi ( strng1, strng2 )

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
!    Input, character ( len = * ) STRNG1, STRNG2, the strings to compare.
!
!    Output, logical S_EQI, the result of the comparison.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) len1
  integer ( kind = 4 ) len2
  integer ( kind = 4 ) lenc
  logical s_eqi
  character s1
  character s2
  character ( len = * ) strng1
  character ( len = * ) strng2

  len1 = len ( strng1 )
  len2 = len ( strng2 )
  lenc = min ( len1, len2 )

  s_eqi = .false.

  do i = 1, lenc

    s1 = strng1(i:i)
    s2 = strng2(i:i)
    call ch_cap ( s1 )
    call ch_cap ( s2 )

    if ( s1 /= s2 ) then
      return
    end if

  end do

  do i = lenc + 1, len1
    if ( strng1(i:i) /= ' ' ) then
      return
    end if
  end do

  do i = lenc + 1, len2
    if ( strng2(i:i) /= ' ' ) then
      return
    end if
  end do

  s_eqi = .true.

  return
end
subroutine s_to_i4 ( string, ival, ierror, last )

!*****************************************************************************80
!
!! S_TO_I4 reads an I4 from a string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) STRING, a string to be examined.
!
!    Output, integer ( kind = 4 ) IVAL, the integer value read from the string.
!    If STRING is blank, then IVAL will be returned 0.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!    0, no error.
!    1, an error occurred.
!
!    Output, integer ( kind = 4 ) LAST, the last character in STRING that was
!    part of the representation of IVAL.
!
  implicit none

  character c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) istate
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) last
  integer ( kind = 4 ) lens
  character ( len = * ) string

  ierror = 0
  istate = 0

  isgn = 1
  ival = 0

  lens = len ( string )

  i = 0

  do

    i = i + 1

    c = string(i:i)

    if ( istate == 0 ) then

      if ( c == ' ' ) then

      else if ( c == '-' ) then
        istate = 1
        isgn = -1
      else if ( c == '+' ) then
        istate = 1
        isgn = + 1
      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        istate = 2
        ival = ichar ( c ) - ichar ( '0' )
      else
        ierror = 1
        return
      end if

    else if ( istate == 1 ) then

      if ( c == ' ' ) then

      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        istate = 2
        ival = ichar ( c ) - ichar ( '0' )
      else
        ierror = 1
        return
      end if

    else if ( istate == 2 ) then

      if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        ival = 10 * ival + ichar ( c ) - ichar ( '0' )
      else
        istate = 3
      end if

    end if
!
!  Continue or exit?
!
    if ( istate == 3 ) then
      ival = isgn * ival
      last = i - 1
      return
    else if ( lens <= i ) then
      if ( istate == 2 ) then
        ival = isgn * ival
        last = lens
      else
        ierror = 1
        last = 0
      end if
      return
    end if

  end do

  return
end
subroutine s_to_r8 ( s, r, ierror, lchar )

!*****************************************************************************80
!
!! S_TO_R8 reads an R8 from a string.
!
!  Discussion:
!
!    This routine will read as many characters as possible until it reaches
!    the end of the string, or encounters a character which cannot be
!    part of the real number.
!
!    Legal input is:
!
!       1 blanks,
!       2 '+' or '-' sign,
!       2.5 spaces
!       3 integer part,
!       4 decimal point,
!       5 fraction part,
!       6 'E' or 'e' or 'D' or 'd', exponent marker,
!       7 exponent sign,
!       8 exponent integer part,
!       9 exponent decimal point,
!      10 exponent fraction part,
!      11 blanks,
!      12 final comma or semi4colon.
!
!    with most quantities optional.
!
!  Example:
!
!    S                 R
!
!    '1'               1.0
!    '     1   '       1.0
!    '1A'              1.0
!    '12,34,56'        12.0
!    '  34 7'          34.0
!    '-1E2ABCD'        -100.0
!    '-1X2ABCD'        -1.0
!    ' 2E-1'           0.2
!    '23.45'           23.45
!    '-4.2E+2'         -420.0
!    '17d2'            1700.0
!    '-14e-2'         -0.14
!    'e2'              100.0
!    '-12.73e-9.23'   -12.73 * 10.0**(-9.23)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string containing the
!    data to be read.  Reading will begin at position 1 and
!    terminate at the end of the string, or when no more
!    characters can be read to form a legal real.  Blanks,
!    commas, or other nonnumeric data will, in particular,
!    cause the conversion to halt.
!
!    Output, real ( kind = 8 ) R, the real value that was read from the string.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!
!    0, no errors occurred.
!
!    1, 2, 6 or 7, the input number was garbled.  The
!    value of IERROR is the last type of input successfully
!    read.  For instance, 1 means initial blanks, 2 means
!    a plus or minus sign, and so on.
!
!    Output, integer ( kind = 4 ) LCHAR, the number of characters read from
!    the string to form the number, including any terminating
!    characters such as a trailing comma or blanks.
!
  implicit none

  logical ch_eqi
  character c
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ihave
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) iterm
  integer ( kind = 4 ) jbot
  integer ( kind = 4 ) jsgn
  integer ( kind = 4 ) jtop
  integer ( kind = 4 ) lchar
  integer ( kind = 4 ) nchar
  integer ( kind = 4 ) ndig
  real ( kind = 8 ) r
  real ( kind = 8 ) rbot
  real ( kind = 8 ) rexp
  real ( kind = 8 ) rtop
  character ( len = * ) s
  character, parameter :: TAB = char ( 9 )

  nchar = len_trim ( s )
  ierror = 0
  r = 0.0D+00
  lchar = - 1
  isgn = 1
  rtop = 0.0D+00
  rbot = 1.0D+00
  jsgn = 1
  jtop = 0
  jbot = 1
  ihave = 1
  iterm = 0

  do

    lchar = lchar + 1
    c = s(lchar+1:lchar+1)
!
!  Blank or TAB character.
!
    if ( c == ' ' .or. c == TAB ) then

      if ( ihave == 2 ) then

      else if ( ihave == 6 .or. ihave == 7 ) then
        iterm = 1
      else if ( 1 < ihave ) then
        ihave = 11
      end if
!
!  Comma.
!
    else if ( c == ',' .or. c == ';' ) then

      if ( ihave /= 1 ) then
        iterm = 1
        ihave = 12
        lchar = lchar + 1
      end if
!
!  Minus sign.
!
    else if ( c == '-' ) then

      if ( ihave == 1 ) then
        ihave = 2
        isgn = - 1
      else if ( ihave == 6 ) then
        ihave = 7
        jsgn = - 1
      else
        iterm = 1
      end if
!
!  Plus sign.
!
    else if ( c == '+' ) then

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
    else if ( c == '.' ) then

      if ( ihave < 4 ) then
        ihave = 4
      else if ( 6 <= ihave .and. ihave <= 8 ) then
        ihave = 9
      else
        iterm = 1
      end if
!
!  Exponent marker.
!
    else if ( ch_eqi ( c, 'E' ) .or. ch_eqi ( c, 'D' ) ) then

      if ( ihave < 6 ) then
        ihave = 6
      else
        iterm = 1
      end if
!
!  Digit.
!
    else if ( ihave < 11 .and. lge ( c, '0' ) .and. lle ( c, '9' ) ) then

      if ( ihave <= 2 ) then
        ihave = 3
      else if ( ihave == 4 ) then
        ihave = 5
      else if ( ihave == 6 .or. ihave == 7 ) then
        ihave = 8
      else if ( ihave == 9 ) then
        ihave = 10
      end if

      call ch_to_digit ( c, ndig )

      if ( ihave == 3 ) then
        rtop = 10.0D+00 * rtop + real ( ndig, kind = 8 )
      else if ( ihave == 5 ) then
        rtop = 10.0D+00 * rtop + real ( ndig, kind = 8 )
        rbot = 10.0D+00 * rbot
      else if ( ihave == 8 ) then
        jtop = 10 * jtop + ndig
      else if ( ihave == 10 ) then
        jtop = 10 * jtop + ndig
        jbot = 10 * jbot
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
    if ( iterm == 1 .or. nchar <= lchar+1 ) then
      exit
    end if

  end do
!
!  If we haven't seen a terminator, and we have examined the
!  entire string, then we're done, and LCHAR is equal to NCHAR.
!
  if ( iterm /= 1 .and. lchar+1 == nchar ) then
    lchar = nchar
  end if
!
!  Number seems to have terminated.  Have we got a legal number?
!  Not if we terminated in states 1, 2, 6 or 7!
!
  if ( ihave == 1 .or. ihave == 2 .or. ihave == 6 .or. ihave == 7 ) then

    ierror = ihave

    return
  end if
!
!  Number seems OK.  Form it.
!
  if ( jtop == 0 ) then
    rexp = 1.0D+00
  else

    if ( jbot == 1 ) then
      rexp = 10.0D+00**( jsgn * jtop )
    else
      rexp = jsgn * jtop
      rexp = rexp / jbot
      rexp = 10.0D+00**rexp
    end if

  end if

  r = isgn * rexp * rtop / rbot

  return
end
subroutine shape_2d_edges_to_ps ( plotxmin2, plotymin2, alpha, iunit, &
  max_order, nface, nnode, face, face_order, x, y, xmin, ymin )

!*****************************************************************************80
!
!! SHAPE_2D_EDGES_TO_PS writes 2D shape edges to a PostScript file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PLOTXMIN2, PLOTYMIN2, the Postscript origin.
!
!    Input, real ( kind = 8 ) ALPHA, the physical-to-Postscript scale factor.
!
!    Input, integer ( kind = 4 ) IUNIT, the output FORTRAN unit.
!
!    Input, integer ( kind = 4 ) MAX_ORDER, the maximum number of nodes per face.
!
!    Input, integer ( kind = 4 ) NFACE, the number of faces.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) FACE(MAX_ORDER,NFACE), the nodes making faces.
!
!    Input, integer ( kind = 4 ) FACE_ORDER(NFACE), the number of nodes per face.
!
!    Input, real ( kind = 8 ) X(NNODE), Y(NNODE), the coordinates of points.
!
!    Input, real ( kind = 8 ) XMIN, YMIN, the physical origin.
!
  implicit none

  integer ( kind = 4 ) max_order
  integer ( kind = 4 ) nface
  integer ( kind = 4 ) nnode

  real ( kind = 8 ) alpha
  integer ( kind = 4 ) face(max_order,nface)
  integer ( kind = 4 ) face_order(nface)
  integer ( kind = 4 ) iface
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) j
  integer ( kind = 4 ) node
  integer ( kind = 4 ) plotxmin2
  integer ( kind = 4 ) plotymin2
  integer ( kind = 4 ) px
  integer ( kind = 4 ) py
  real ( kind = 8 ) x(nnode)
  real ( kind = 8 ) xmin
  real ( kind = 8 ) y(nnode)
  real ( kind = 8 ) ymin
!
!  Draw faces and fill them.
!
  do iface = 1, nface

    write ( iunit, '(a)' ) 'newpath'

    node = face(face_order(iface),iface)
    px = plotxmin2 + nint ( alpha * ( x(node) - xmin ) )
    py = plotymin2 + nint ( alpha * ( y(node) - ymin ) )
    write ( iunit, '(2i4,a,2i4,a)' ) px, py, ' moveto '

    do j = 1, face_order(iface)

      node = face(j,iface)
      px = plotxmin2 + nint ( alpha * ( x(node) - xmin ) )
      py = plotymin2 + nint ( alpha * ( y(node) - ymin ) )
      write ( iunit, '(2i4,a,2i4,a)' ) px, py, ' lineto '

    end do

    write ( iunit, '(a)' ) 'stroke'
!   write ( iunit, '(a)' ) 'fill'

  end do

  return
end
subroutine shape_2d_faces_to_ps ( plotxmin2, plotymin2, alpha, iunit, &
  max_order, nface, nnode, face, face_order, x, y, xmin, ymin )

!*****************************************************************************80
!
!! SHAPE_2D_FACES_TO_PS writes 2D shape faces to a PostScript file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PLOTXMIN2, PLOTYMIN2, the Postscript origin.
!
!    Input, real ( kind = 8 ) ALPHA, the physical-to-Postscript scale factor.
!
!    Input, integer ( kind = 4 ) IUNIT, the output FORTRAN unit.
!
!    Input, integer ( kind = 4 ) MAX_ORDER, the maximum number of nodes per face.
!
!    Input, integer ( kind = 4 ) NFACE, the number of faces.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) FACE(MAX_ORDER,NFACE), the nodes making faces.
!
!    Input, integer ( kind = 4 ) FACE_ORDER(NFACE), the number of nodes per face.
!
!    Input, real ( kind = 8 ) X(NNODE), Y(NNODE), the coordinates of points.
!
!    Input, real ( kind = 8 ) XMIN, YMIN, the physical origin.
!
  implicit none

  integer ( kind = 4 ) max_order
  integer ( kind = 4 ) nface
  integer ( kind = 4 ) nnode

  real ( kind = 8 ) alpha
  real ( kind = 8 ) blue
  integer ( kind = 4 ) face(max_order,nface)
  real ( kind = 8 ) green
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iface
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) j
  integer ( kind = 4 ) node
  integer ( kind = 4 ) face_order(nface)
  integer ( kind = 4 ) plotxmin2
  integer ( kind = 4 ) plotymin2
  integer ( kind = 4 ) px
  integer ( kind = 4 ) py
  real ( kind = 8 ) red
  real ( kind = 8 ) x(nnode)
  real ( kind = 8 ) xmin
  real ( kind = 8 ) y(nnode)
  real ( kind = 8 ) ymin
!
!  Draw the faces.
!
  do iface = 1, nface

    do i = 1, 2

      if ( i == 1 ) then
        red = 0.9D+00
        green = 0.9D+00
        blue = 1.0D+00
      else
        red = 0.0D+00
        green = 0.0D+00
        blue = 0.0D+00
      end if

      write ( iunit, '(3f7.4,a)' ) red, green, blue, ' setrgbcolor'

      write ( iunit, '(a)' ) 'newpath'

      node = face(face_order(iface),iface)
      px = plotxmin2 + nint ( alpha * ( x(node) - xmin ) )
      py = plotymin2 + nint ( alpha * ( y(node) - ymin ) )
      write ( iunit, '(2i4,a,2i4,a)' ) px, py, ' moveto '

      do j = 1, face_order(iface)

        node = face(j,iface)
        px = plotxmin2 + nint ( alpha * ( x(node) - xmin ) )
        py = plotymin2 + nint ( alpha * ( y(node) - ymin ) )
        write ( iunit, '(2i4,a,2i4,a)' ) px, py, ' lineto '

      end do

      if ( i == 1 ) then
        write ( iunit, '(a)' ) 'fill'
      else
        write ( iunit, '(a)' ) 'stroke'
      end if

    end do

  end do

  return
end
subroutine shape_2d_nodes_to_ps ( plotxmin2, plotymin2, alpha, iunit, &
  nnode, x, y, xmin, ymin )

!*****************************************************************************80
!
!! SHAPE_2D_NODES_TO_PS writes 2D shape nodes to a PostScript file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PLOTXMIN2, PLOTYMIN2, the Postscript origin.
!
!    Input, real ( kind = 8 ) ALPHA, the physical-to-Postscript scale factor.
!
!    Input, integer ( kind = 4 ) IUNIT, the output FORTRAN unit.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, real ( kind = 8 ) X(NNODE), Y(NNODE), the X and Y components
!    of points.
!
!    Input, real ( kind = 8 ) XMIN, YMIN, the physical origin.
!
  implicit none

  integer ( kind = 4 ) nnode

  real ( kind = 8 ) alpha
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) plotxmin2
  integer ( kind = 4 ) plotymin2
  integer ( kind = 4 ) px
  integer ( kind = 4 ) py
  real ( kind = 8 ) x(nnode)
  real ( kind = 8 ) xmin
  real ( kind = 8 ) y(nnode)
  real ( kind = 8 ) ymin
!
!  Draw the nodes.
!
  do i = 1, nnode

    px = plotxmin2 + nint ( alpha * ( x(i) - xmin ) )
    py = plotymin2 + nint ( alpha * ( y(i) - ymin ) )

    write ( iunit, '(a,2i4,a)' ) 'newpath ', px, py, &
      ' 5 0 360 arc closepath stroke'

  end do

  return
end
subroutine shape_3d_edges_to_ps ( file_name, max_order, nface, nnode, &
  face, face_order, x, y, z )

!*****************************************************************************80
!
!! SHAPE_3D_EDGES_TO_PS writes 3D shape edges to a PostScript file.
!
!  Discussion:
!
!    Four views are created in one picture: XY, YZ, ZX, and XYZ.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the output file.
!
!    Input, integer ( kind = 4 ) MAX_ORDER, the maximum number of nodes per face.
!
!    Input, integer ( kind = 4 ) NFACE, the number of faces.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) FACE(MAX_ORDER,NFACE), the nodes making faces.
!
!    Input, integer ( kind = 4 ) FACE_ORDER(NFACE), the number of nodes per face.
!
!    Input, real ( kind = 8 ) X(NNODE), Y(NNODE), Z(NNODE), the X, Y and Z
!    components of points.
!
  implicit none

  integer ( kind = 4 ) max_order
  integer ( kind = 4 ) nface
  integer ( kind = 4 ) nnode

  real ( kind = 8 ) alpha
  real ( kind = 8 ) blue
  character ( len = 8 ) date
  integer ( kind = 4 ) face(max_order,nface)
  integer ( kind = 4 ) face_order(nface)
  character ( len = * ) file_name
  real ( kind = 8 ) green
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ), parameter :: margin = 36
  integer ( kind = 4 ) pagexmax
  integer ( kind = 4 ) pagexmin
  integer ( kind = 4 ) pageymax
  integer ( kind = 4 ) pageymin
  integer ( kind = 4 ) plotxmax
  integer ( kind = 4 ) plotxmin
  integer ( kind = 4 ) plotxmin2
  integer ( kind = 4 ) plotymax
  integer ( kind = 4 ) plotymin
  integer ( kind = 4 ) plotymin2
  integer ( kind = 4 ) px1
  integer ( kind = 4 ) px2
  integer ( kind = 4 ) px3
  integer ( kind = 4 ) px4
  integer ( kind = 4 ) px5
  integer ( kind = 4 ) py1
  integer ( kind = 4 ) py2
  integer ( kind = 4 ) py3
  integer ( kind = 4 ) py4
  integer ( kind = 4 ) py5
  real ( kind = 8 ) red
  real ( kind = 8 ) x(nnode)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) xx(nnode)
  real ( kind = 8 ) y(nnode)
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin
  real ( kind = 8 ) yy(nnode)
  real ( kind = 8 ) z(nnode)
!
!  Open the file.
!
  call get_unit ( iunit )

  open ( unit = iunit, file = file_name, status = 'replace', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SHAPE_3D_EDGES_TO_PS - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file.'
    return
  end if
!
!  Write the prolog.
!
  pagexmax = 612
  pagexmin = 0
  pageymax = 792
  pageymin = 0

  px1 = 0
  px2 = margin
  px3 = pagexmax / 2
  px4 = pagexmax - margin
  px5 = pagexmax

  py1 = 0
  py2 = margin
  py3 = pageymax / 2
  py4 = pageymax - margin
  py5 = pageymax

  write ( iunit, '(a)' ) '%!PS-Adobe-3.0'
  write ( iunit, '(a)' ) '%%Document-Fonts: Times-Roman'
  write ( iunit, '(a,a)' ) '%%Title: ' , trim ( file_name )
  write ( iunit, '(a)' ) '%%Creator: GRAFPACK(shape_3d_edges_to_ps)'
  call date_and_time ( date )
  write ( iunit, '(a)' ) '%%CreationDate: ' // trim ( date )
  write ( iunit, '(a)' ) '%%BoundingBox 0 0 612 794'
  write ( iunit, '(a)' ) '%%LanguageLevel: 2'
  write ( iunit, '(a)' ) '%%EndComments'
  write ( iunit, '(a)' ) '%%BeginProlog'
  write ( iunit, '(a)' ) '%%EndProlog'
!
!  Draw gray lines to separate the boxes.
!
  red = 0.5
  green = 0.5
  blue = 0.5

  write ( iunit, '(3f7.4,a)' ) red, green, blue, ' setrgbcolor'

  write ( iunit, '(2i4,a)' ) px2, py2, ' moveto '
  write ( iunit, '(2i4,a)' ) px4, py2, ' lineto '
  write ( iunit, '(2i4,a)' ) px2, py3, ' moveto '
  write ( iunit, '(2i4,a)' ) px4, py3, ' lineto '
  write ( iunit, '(2i4,a)' ) px2, py4, ' moveto '
  write ( iunit, '(2i4,a)' ) px4, py4, ' lineto '

  write ( iunit, '(2i4,a)' ) px2, py2, ' moveto '
  write ( iunit, '(2i4,a)' ) px2, py4, ' lineto '
  write ( iunit, '(2i4,a)' ) px3, py2, ' moveto '
  write ( iunit, '(2i4,a)' ) px3, py4, ' lineto '
  write ( iunit, '(2i4,a)' ) px4, py2, ' moveto '
  write ( iunit, '(2i4,a)' ) px4, py4, ' lineto '
  write ( iunit, '(a)' ) 'stroke'
!
!  Determine ALPHA, the single scale factor to be used for both
!  directions, and all four plots!
!
  xx(1:nnode) = x(1:nnode)
  yy(1:nnode) = y(1:nnode)
  xmin = minval ( xx(1:nnode) )
  xmax = maxval ( xx(1:nnode) )
  ymin = minval ( yy(1:nnode) )
  ymax = maxval ( yy(1:nnode) )

  alpha = min ( ( px3 - px2 ) / ( xmax - xmin ), &
                ( py4 - py3 ) / ( ymax - ymin ) )

  xx(1:nnode) = y(1:nnode)
  yy(1:nnode) = z(1:nnode)
  xmin = minval ( xx(1:nnode) )
  xmax = maxval ( xx(1:nnode) )
  ymin = minval ( yy(1:nnode) )
  ymax = maxval ( yy(1:nnode) )

  alpha = min ( alpha, &
                ( px4 - px3 ) / ( xmax - xmin ), &
                ( py4 - py3 ) / ( ymax - ymin ) )

  xx(1:nnode) = z(1:nnode)
  yy(1:nnode) = x(1:nnode)
  xmin = minval ( xx(1:nnode) )
  xmax = maxval ( xx(1:nnode) )
  ymin = minval ( yy(1:nnode) )
  ymax = maxval ( yy(1:nnode) )

  alpha = min ( alpha, &
                ( px3 - px2 ) / ( xmax - xmin ), &
                ( py3 - py2 ) / ( ymax - ymin ) )

  xx(1:nnode) = 0.80 * x(1:nnode) - 0.31 * y(1:nnode) + 0.50 * z(1:nnode)
  yy(1:nnode) = 0.50 * x(1:nnode) + 0.80 * y(1:nnode) - 0.31 * z(1:nnode)
  xmin = minval ( xx(1:nnode) )
  xmax = maxval ( xx(1:nnode) )
  ymin = minval ( yy(1:nnode) )
  ymax = maxval ( yy(1:nnode) )

  alpha = min ( alpha, &
                ( px4 - px3 ) / ( xmax - xmin ), &
                ( py3 - py2 ) / ( ymax - ymin ) )
!
!  Set the fill color.
!
  red = 0.9D+00
  green = 0.9D+00
  blue = 1.0D+00

  write ( iunit, '(3f7.4,a)' ) red, green, blue, ' setrgbcolor'
!
!  XY edges.
!
  plotxmin = px2
  plotxmax = px3
  plotymin = py3
  plotymax = py4

  xx(1:nnode) = x(1:nnode)
  yy(1:nnode) = y(1:nnode)
  xmin = minval ( xx(1:nnode) )
  xmax = maxval ( xx(1:nnode) )
  ymin = minval ( yy(1:nnode) )
  ymax = maxval ( yy(1:nnode) )

  plotxmin2 = 0.5 * ( plotxmin + plotxmax - alpha * ( xmax - xmin ) )
  plotymin2 = 0.5 * ( plotymin + plotymax - alpha * ( ymax - ymin ) )

  call shape_2d_edges_to_ps ( plotxmin2, plotymin2, alpha, iunit, &
    max_order, nface, nnode, face, face_order, xx, yy, xmin, ymin )
!
!  YZ edges.
!
  plotxmin = px3
  plotxmax = px4
  plotymin = py3
  plotymax = py4

  xx(1:nnode) = y(1:nnode)
  yy(1:nnode) = z(1:nnode)
  xmin = minval ( xx(1:nnode) )
  xmax = maxval ( xx(1:nnode) )
  ymin = minval ( yy(1:nnode) )
  ymax = maxval ( yy(1:nnode) )

  plotxmin2 = 0.5 * ( plotxmin + plotxmax - alpha * ( xmax - xmin ) )
  plotymin2 = 0.5 * ( plotymin + plotymax - alpha * ( ymax - ymin ) )

  call shape_2d_edges_to_ps ( plotxmin2, plotymin2, alpha, iunit, &
    max_order, nface, nnode, face, face_order, xx, yy, xmin, ymin )
!
!  ZX edges.
!
  plotxmin = px2
  plotxmax = px3
  plotymin = py2
  plotymax = py3

  xx(1:nnode) = z(1:nnode)
  yy(1:nnode) = x(1:nnode)
  xmin = minval ( xx(1:nnode) )
  xmax = maxval ( xx(1:nnode) )
  ymin = minval ( yy(1:nnode) )
  ymax = maxval ( yy(1:nnode) )

  plotxmin2 = 0.5 * ( plotxmin + plotxmax - alpha * ( xmax - xmin ) )
  plotymin2 = 0.5 * ( plotymin + plotymax - alpha * ( ymax - ymin ) )

  call shape_2d_edges_to_ps ( plotxmin2, plotymin2, alpha, iunit, &
    max_order, nface, nnode, face, face_order, xx, yy, xmin, ymin )
!
!  XYZ edges.
!
  plotxmin = px3
  plotxmax = px4
  plotymin = py2
  plotymax = py3

  xx(1:nnode) = 0.80 * x(1:nnode) - 0.31 * y(1:nnode) + 0.50 * z(1:nnode)
  yy(1:nnode) = 0.50 * x(1:nnode) + 0.80 * y(1:nnode) - 0.31 * z(1:nnode)
  xmin = minval ( xx(1:nnode) )
  xmax = maxval ( xx(1:nnode) )
  ymin = minval ( yy(1:nnode) )
  ymax = maxval ( yy(1:nnode) )

  plotxmin2 = 0.5 * ( plotxmin + plotxmax - alpha * ( xmax - xmin ) )
  plotymin2 = 0.5 * ( plotymin + plotymax - alpha * ( ymax - ymin ) )

  call shape_2d_edges_to_ps ( plotxmin2, plotymin2, alpha, iunit, &
    max_order, nface, nnode, face, face_order, xx, yy, xmin, ymin )

  write ( iunit, '(a)' ) 'showpage'
!
!  Write the epilog.
!
  write ( iunit, '(a)' ) 'grestore'
  write ( iunit, '(a)' ) '%%Trailer'
  write ( iunit, '(a,i2)' ) '%%Pages: 1'

  close ( unit = iunit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SHAPE_3D_EDGES_TO_PS'
  write ( *, '(a)' ) '  The data was written to the file: ' &
    // trim ( file_name )

  return
end
subroutine shape_3d_faces_to_ps ( file_name, max_order, nface, nnode, face, &
  face_order, x, y, z )

!*****************************************************************************80
!
!! SHAPE_3D_FACES_TO_PS writes 3D shape faces to a PostScript file.
!
!  Discussion:
!
!    Four views are created in one picture: XY, YZ, ZX, and XYZ.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the output file.
!
!    Input, integer ( kind = 4 ) MAX_ORDER, the maximum number of nodes per face.
!
!    Input, integer ( kind = 4 ) NFACE, the number of faces.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) FACE(MAX_ORDER,NFACE), the nodes making faces.
!
!    Input, integer ( kind = 4 ) FACE_ORDER(NFACE), the number of nodes per face.
!
!    Input, real ( kind = 8 ) X(NNODE), Y(NNODE), Z(NNODE), the X, Y and Z
!    components of points.
!
  implicit none

  integer ( kind = 4 ) max_order
  integer ( kind = 4 ) nface
  integer ( kind = 4 ) nnode

  real ( kind = 8 ) alpha
  real ( kind = 8 ) blue
  character ( len = 8 ) date
  integer ( kind = 4 ) face(max_order,nface)
  integer ( kind = 4 ) face_order(nface)
  character ( len = * ) file_name
  real ( kind = 8 ) green
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ), parameter :: margin = 36
  integer ( kind = 4 ) pagexmax
  integer ( kind = 4 ) pagexmin
  integer ( kind = 4 ) pageymax
  integer ( kind = 4 ) pageymin
  integer ( kind = 4 ) plotxmax
  integer ( kind = 4 ) plotxmin
  integer ( kind = 4 ) plotxmin2
  integer ( kind = 4 ) plotymax
  integer ( kind = 4 ) plotymin
  integer ( kind = 4 ) plotymin2
  integer ( kind = 4 ) px1
  integer ( kind = 4 ) px2
  integer ( kind = 4 ) px3
  integer ( kind = 4 ) px4
  integer ( kind = 4 ) px5
  integer ( kind = 4 ) py1
  integer ( kind = 4 ) py2
  integer ( kind = 4 ) py3
  integer ( kind = 4 ) py4
  integer ( kind = 4 ) py5
  real ( kind = 8 ) red
  real ( kind = 8 ) x(nnode)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) xx(nnode)
  real ( kind = 8 ) y(nnode)
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin
  real ( kind = 8 ) yy(nnode)
  real ( kind = 8 ) z(nnode)
!
!  Open the file.
!
  call get_unit ( iunit )

  open ( unit = iunit, file = file_name, status = 'replace', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SHAPE_3D_EDGES_TO_PS - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file.'
    return
  end if
!
!  Write the prolog.
!
  pagexmax = 612
  pagexmin = 0
  pageymax = 792
  pageymin = 0

  px1 = 0
  px2 = margin
  px3 = pagexmax / 2
  px4 = pagexmax - margin
  px5 = pagexmax

  py1 = 0
  py2 = margin
  py3 = pageymax / 2
  py4 = pageymax - margin
  py5 = pageymax

  write ( iunit, '(a)' ) '%!PS-Adobe-3.0'
  write ( iunit, '(a)' ) '%%Document-Fonts: Times-Roman'
  write ( iunit, '(a,a)' ) '%%Title: ' , trim ( file_name )
  write ( iunit, '(a)' ) '%%Creator: GRAFPACK(shape_3d_edges_to_ps)'
  call date_and_time ( date )
  write ( iunit, '(a)' ) '%%CreationDate: ' // trim ( date )
  write ( iunit, '(a)' ) '%%BoundingBox 0 0 612 794'
  write ( iunit, '(a)' ) '%%LanguageLevel: 2'
  write ( iunit, '(a)' ) '%%EndComments'
  write ( iunit, '(a)' ) '%%BeginProlog'
  write ( iunit, '(a)' ) '%%EndProlog'
!
!  Draw gray lines to separate the boxes.
!
  red = 0.5D+00
  green = 0.5D+00
  blue = 0.5D+00

  write ( iunit, '(3f7.4,a)' ) red, green, blue, ' setrgbcolor'

  write ( iunit, '(2i4,a)' ) px2, py2, ' moveto '
  write ( iunit, '(2i4,a)' ) px4, py2, ' lineto '
  write ( iunit, '(2i4,a)' ) px2, py3, ' moveto '
  write ( iunit, '(2i4,a)' ) px4, py3, ' lineto '
  write ( iunit, '(2i4,a)' ) px2, py4, ' moveto '
  write ( iunit, '(2i4,a)' ) px4, py4, ' lineto '

  write ( iunit, '(2i4,a)' ) px2, py2, ' moveto '
  write ( iunit, '(2i4,a)' ) px2, py4, ' lineto '
  write ( iunit, '(2i4,a)' ) px3, py2, ' moveto '
  write ( iunit, '(2i4,a)' ) px3, py4, ' lineto '
  write ( iunit, '(2i4,a)' ) px4, py2, ' moveto '
  write ( iunit, '(2i4,a)' ) px4, py4, ' lineto '
  write ( iunit, '(a)' ) 'stroke'
!
!  Determine ALPHA, the single scale factor to be used for both
!  directions, and all four plots!
!
  xx(1:nnode) = x(1:nnode)
  yy(1:nnode) = y(1:nnode)
  xmin = minval ( xx(1:nnode) )
  xmax = maxval ( xx(1:nnode) )
  ymin = minval ( yy(1:nnode) )
  ymax = maxval ( yy(1:nnode) )

  alpha = min ( ( px3 - px2 ) / ( xmax - xmin ), &
                ( py4 - py3 ) / ( ymax - ymin ) )

  xx(1:nnode) = y(1:nnode)
  yy(1:nnode) = z(1:nnode)
  xmin = minval ( xx(1:nnode) )
  xmax = maxval ( xx(1:nnode) )
  ymin = minval ( yy(1:nnode) )
  ymax = maxval ( yy(1:nnode) )

  alpha = min ( alpha, &
                ( px4 - px3 ) / ( xmax - xmin ), &
                ( py4 - py3 ) / ( ymax - ymin ) )

  xx(1:nnode) = z(1:nnode)
  yy(1:nnode) = x(1:nnode)
  xmin = minval ( xx(1:nnode) )
  xmax = maxval ( xx(1:nnode) )
  ymin = minval ( yy(1:nnode) )
  ymax = maxval ( yy(1:nnode) )

  alpha = min ( alpha, &
                ( px3 - px2 ) / ( xmax - xmin ), &
                ( py3 - py2 ) / ( ymax - ymin ) )

  xx(1:nnode) = 0.80 * x(1:nnode) - 0.31 * y(1:nnode) + 0.50 * z(1:nnode)
  yy(1:nnode) = 0.50 * x(1:nnode) + 0.80 * y(1:nnode) - 0.31 * z(1:nnode)
  xmin = minval ( xx(1:nnode) )
  xmax = maxval ( xx(1:nnode) )
  ymin = minval ( yy(1:nnode) )
  ymax = maxval ( yy(1:nnode) )

  alpha = min ( alpha, &
                ( px4 - px3 ) / ( xmax - xmin ), &
                ( py3 - py2 ) / ( ymax - ymin ) )
!
!  Set the fill color.
!
  red = 0.9D+00
  green = 0.9D+00
  blue = 1.0D+00

  write ( iunit, '(3f7.4,a)' ) red, green, blue, ' setrgbcolor'
!
!  XY edges.
!
  plotxmin = px2
  plotxmax = px3
  plotymin = py3
  plotymax = py4

  xx(1:nnode) = x(1:nnode)
  yy(1:nnode) = y(1:nnode)
  xmin = minval ( xx(1:nnode) )
  xmax = maxval ( xx(1:nnode) )
  ymin = minval ( yy(1:nnode) )
  ymax = maxval ( yy(1:nnode) )

  plotxmin2 = 0.5 * ( plotxmin + plotxmax - alpha * ( xmax - xmin ) )
  plotymin2 = 0.5 * ( plotymin + plotymax - alpha * ( ymax - ymin ) )

  call shape_2d_faces_to_ps ( plotxmin2, plotymin2, alpha, iunit, &
    max_order, nface, nnode, face, face_order, xx, yy, xmin, ymin )
!
!  YZ edges.
!
  plotxmin = px3
  plotxmax = px4
  plotymin = py3
  plotymax = py4

  xx(1:nnode) = y(1:nnode)
  yy(1:nnode) = z(1:nnode)
  xmin = minval ( xx(1:nnode) )
  xmax = maxval ( xx(1:nnode) )
  ymin = minval ( yy(1:nnode) )
  ymax = maxval ( yy(1:nnode) )

  plotxmin2 = 0.5 * ( plotxmin + plotxmax - alpha * ( xmax - xmin ) )
  plotymin2 = 0.5 * ( plotymin + plotymax - alpha * ( ymax - ymin ) )

  call shape_2d_faces_to_ps ( plotxmin2, plotymin2, alpha, iunit, &
    max_order, nface, nnode, face, face_order, xx, yy, xmin, ymin )
!
!  ZX edges.
!
  plotxmin = px2
  plotxmax = px3
  plotymin = py2
  plotymax = py3

  xx(1:nnode) = z(1:nnode)
  yy(1:nnode) = x(1:nnode)
  xmin = minval ( xx(1:nnode) )
  xmax = maxval ( xx(1:nnode) )
  ymin = minval ( yy(1:nnode) )
  ymax = maxval ( yy(1:nnode) )

  plotxmin2 = 0.5 * ( plotxmin + plotxmax - alpha * ( xmax - xmin ) )
  plotymin2 = 0.5 * ( plotymin + plotymax - alpha * ( ymax - ymin ) )

  call shape_2d_faces_to_ps ( plotxmin2, plotymin2, alpha, iunit, &
    max_order, nface, nnode, face, face_order, xx, yy, xmin, ymin )
!
!  XYZ edges.
!
  plotxmin = px3
  plotxmax = px4
  plotymin = py2
  plotymax = py3

  xx(1:nnode) = 0.80 * x(1:nnode) - 0.31 * y(1:nnode) + 0.50 * z(1:nnode)
  yy(1:nnode) = 0.50 * x(1:nnode) + 0.80 * y(1:nnode) - 0.31 * z(1:nnode)
  xmin = minval ( xx(1:nnode) )
  xmax = maxval ( xx(1:nnode) )
  ymin = minval ( yy(1:nnode) )
  ymax = maxval ( yy(1:nnode) )

  plotxmin2 = 0.5 * ( plotxmin + plotxmax - alpha * ( xmax - xmin ) )
  plotymin2 = 0.5 * ( plotymin + plotymax - alpha * ( ymax - ymin ) )

  call shape_2d_faces_to_ps ( plotxmin2, plotymin2, alpha, iunit, &
    max_order, nface, nnode, face, face_order, xx, yy, xmin, ymin )

  write ( iunit, '(a)' ) 'showpage'
!
!  Write the epilog.
!
  write ( iunit, '(a)' ) 'grestore'
  write ( iunit, '(a)' ) '%%Trailer'
  write ( iunit, '(a,i2)' ) '%%Pages: 1'

  close ( unit = iunit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SHAPE_3D_EDGES_TO_PS'
  write ( *, '(a)' ) '  The data was written to the file: ' &
    // trim ( file_name )

  return
end
subroutine shape_3d_nodes_to_ps ( file_name, nnode, x, y, z )

!*****************************************************************************80
!
!! SHAPE_3D_NODES_TO_PS writes 3D shape nodes to a PostScript file.
!
!  Discussion:
!
!    Four views are created in one picture: XY, YZ, ZX, and XYZ.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the output file.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, real ( kind = 8 ) X(NNODE), Y(NNODE), Z(NNODE), the X, Y and Z
!    components of points.
!
  implicit none

  integer ( kind = 4 ) nnode

  real ( kind = 8 ) alpha
  real ( kind = 8 ) blue
  character ( len = 8 ) date
  character ( len = * ) file_name
  real ( kind = 8 ) green
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ), parameter :: margin = 36
  integer ( kind = 4 ) pagexmax
  integer ( kind = 4 ) pagexmin
  integer ( kind = 4 ) pageymax
  integer ( kind = 4 ) pageymin
  integer ( kind = 4 ) plotxmax
  integer ( kind = 4 ) plotxmin
  integer ( kind = 4 ) plotxmin2
  integer ( kind = 4 ) plotymax
  integer ( kind = 4 ) plotymin
  integer ( kind = 4 ) plotymin2
  integer ( kind = 4 ) px1
  integer ( kind = 4 ) px2
  integer ( kind = 4 ) px3
  integer ( kind = 4 ) px4
  integer ( kind = 4 ) px5
  integer ( kind = 4 ) py1
  integer ( kind = 4 ) py2
  integer ( kind = 4 ) py3
  integer ( kind = 4 ) py4
  integer ( kind = 4 ) py5
  real ( kind = 8 ) red
  real ( kind = 8 ) x(nnode)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) xx(nnode)
  real ( kind = 8 ) y(nnode)
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin
  real ( kind = 8 ) yy(nnode)
  real ( kind = 8 ) z(nnode)
!
!  Open the file.
!
  call get_unit ( iunit )

  open ( unit = iunit, file = file_name, status = 'replace', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SHAPE_3D_NODES_TO_PS - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file.'
    return
  end if
!
!  Write the prolog.
!
  pagexmax = 612
  pagexmin = 0
  pageymax = 792
  pageymin = 0

  px1 = 0
  px2 = margin
  px3 = pagexmax / 2
  px4 = pagexmax - margin
  px5 = pagexmax

  py1 = 0
  py2 = margin
  py3 = pageymax / 2
  py4 = pageymax - margin
  py5 = pageymax

  write ( iunit, '(a)' ) '%!PS-Adobe-3.0'
  write ( iunit, '(a)' ) '%%Document-Fonts: Times-Roman'
  write ( iunit, '(a,a)' ) '%%Title: ' , trim ( file_name )
  write ( iunit, '(a)' ) '%%Creator: GRAFPACK(shape_3d_nodes_to_ps)'
  call date_and_time ( date )
  write ( iunit, '(a)' ) '%%CreationDate: ' // trim ( date )
  write ( iunit, '(a)' ) '%%BoundingBox 0 0 612 794'
  write ( iunit, '(a)' ) '%%LanguageLevel: 2'
  write ( iunit, '(a)' ) '%%EndComments'
  write ( iunit, '(a)' ) '%%BeginProlog'
  write ( iunit, '(a)' ) '%%EndProlog'
!
!  Draw gray lines to separate the boxes.
!
  red = 0.5D+00
  green = 0.5D+00
  blue = 0.5D+00

  write ( iunit, '(3f7.4,a)' ) red, green, blue, ' setrgbcolor'

  write ( iunit, '(2i4,a)' ) px2, py2, ' moveto '
  write ( iunit, '(2i4,a)' ) px4, py2, ' lineto '
  write ( iunit, '(2i4,a)' ) px2, py3, ' moveto '
  write ( iunit, '(2i4,a)' ) px4, py3, ' lineto '
  write ( iunit, '(2i4,a)' ) px2, py4, ' moveto '
  write ( iunit, '(2i4,a)' ) px4, py4, ' lineto '

  write ( iunit, '(2i4,a)' ) px2, py2, ' moveto '
  write ( iunit, '(2i4,a)' ) px2, py4, ' lineto '
  write ( iunit, '(2i4,a)' ) px3, py2, ' moveto '
  write ( iunit, '(2i4,a)' ) px3, py4, ' lineto '
  write ( iunit, '(2i4,a)' ) px4, py2, ' moveto '
  write ( iunit, '(2i4,a)' ) px4, py4, ' lineto '
  write ( iunit, '(a)' ) 'stroke'
!
!  Determine ALPHA, the single scale factor to be used for both
!  directions, and all four plots!
!
  xx(1:nnode) = x(1:nnode)
  yy(1:nnode) = y(1:nnode)
  xmin = minval ( xx(1:nnode) )
  xmax = maxval ( xx(1:nnode) )
  ymin = minval ( yy(1:nnode) )
  ymax = maxval ( yy(1:nnode) )

  alpha = min ( ( px3 - px2 ) / ( xmax - xmin ), &
                ( py4 - py3 ) / ( ymax - ymin ) )

  xx(1:nnode) = y(1:nnode)
  yy(1:nnode) = z(1:nnode)
  xmin = minval ( xx(1:nnode) )
  xmax = maxval ( xx(1:nnode) )
  ymin = minval ( yy(1:nnode) )
  ymax = maxval ( yy(1:nnode) )

  alpha = min ( alpha, &
                ( px4 - px3 ) / ( xmax - xmin ), &
                ( py4 - py3 ) / ( ymax - ymin ) )

  xx(1:nnode) = z(1:nnode)
  yy(1:nnode) = x(1:nnode)
  xmin = minval ( xx(1:nnode) )
  xmax = maxval ( xx(1:nnode) )
  ymin = minval ( yy(1:nnode) )
  ymax = maxval ( yy(1:nnode) )

  alpha = min ( alpha, &
                ( px3 - px2 ) / ( xmax - xmin ), &
                ( py3 - py2 ) / ( ymax - ymin ) )

  xx(1:nnode) = 0.80 * x(1:nnode) - 0.31 * y(1:nnode) + 0.50 * z(1:nnode)
  yy(1:nnode) = 0.50 * x(1:nnode) + 0.80 * y(1:nnode) - 0.31 * z(1:nnode)
  xmin = minval ( xx(1:nnode) )
  xmax = maxval ( xx(1:nnode) )
  ymin = minval ( yy(1:nnode) )
  ymax = maxval ( yy(1:nnode) )

  alpha = min ( alpha, &
                ( px4 - px3 ) / ( xmax - xmin ), &
                ( py3 - py2 ) / ( ymax - ymin ) )
!
!  Set the color.
!
  red = 0.3
  green = 0.3
  blue = 0.3

  write ( iunit, '(3f7.4,a)' ) red, green, blue, ' setrgbcolor'
!
!  XY nodes.
!
  plotxmin = px2
  plotxmax = px3
  plotymin = py3
  plotymax = py4

  xx(1:nnode) = x(1:nnode)
  yy(1:nnode) = y(1:nnode)
  xmin = minval ( xx(1:nnode) )
  xmax = maxval ( xx(1:nnode) )
  ymin = minval ( yy(1:nnode) )
  ymax = maxval ( yy(1:nnode) )

  plotxmin2 = 0.5 * ( plotxmin + plotxmax - alpha * ( xmax - xmin ) )
  plotymin2 = 0.5 * ( plotymin + plotymax - alpha * ( ymax - ymin ) )

  call shape_2d_nodes_to_ps ( plotxmin2, plotymin2, alpha, iunit, &
    nnode, xx, yy, xmin, ymin )
!
!  YZ edges.
!
  plotxmin = px3
  plotxmax = px4
  plotymin = py3
  plotymax = py4

  xx(1:nnode) = y(1:nnode)
  yy(1:nnode) = z(1:nnode)
  xmin = minval ( xx(1:nnode) )
  xmax = maxval ( xx(1:nnode) )
  ymin = minval ( yy(1:nnode) )
  ymax = maxval ( yy(1:nnode) )

  plotxmin2 = 0.5 * ( plotxmin + plotxmax - alpha * ( xmax - xmin ) )
  plotymin2 = 0.5 * ( plotymin + plotymax - alpha * ( ymax - ymin ) )

  call shape_2d_nodes_to_ps ( plotxmin2, plotymin2, alpha, iunit, &
    nnode, xx, yy, xmin, ymin )
!
!  ZX edges.
!
  plotxmin = px2
  plotxmax = px3
  plotymin = py2
  plotymax = py3

  xx(1:nnode) = z(1:nnode)
  yy(1:nnode) = x(1:nnode)
  xmin = minval ( xx(1:nnode) )
  xmax = maxval ( xx(1:nnode) )
  ymin = minval ( yy(1:nnode) )
  ymax = maxval ( yy(1:nnode) )

  plotxmin2 = 0.5 * ( plotxmin + plotxmax - alpha * ( xmax - xmin ) )
  plotymin2 = 0.5 * ( plotymin + plotymax - alpha * ( ymax - ymin ) )

  call shape_2d_nodes_to_ps ( plotxmin2, plotymin2, alpha, iunit, &
    nnode, xx, yy, xmin, ymin )
!
!  XYZ edges.
!
  plotxmin = px3
  plotxmax = px4
  plotymin = py2
  plotymax = py3

  xx(1:nnode) = 0.80 * x(1:nnode) - 0.31 * y(1:nnode) + 0.50 * z(1:nnode)
  yy(1:nnode) = 0.50 * x(1:nnode) + 0.80 * y(1:nnode) - 0.31 * z(1:nnode)
  xmin = minval ( xx(1:nnode) )
  xmax = maxval ( xx(1:nnode) )
  ymin = minval ( yy(1:nnode) )
  ymax = maxval ( yy(1:nnode) )

  plotxmin2 = 0.5 * ( plotxmin + plotxmax - alpha * ( xmax - xmin ) )
  plotymin2 = 0.5 * ( plotymin + plotymax - alpha * ( ymax - ymin ) )

  call shape_2d_nodes_to_ps ( plotxmin2, plotymin2, alpha, iunit, &
    nnode, xx, yy, xmin, ymin )

  write ( iunit, '(a)' ) 'showpage'
!
!  Write the epilog.
!
  write ( iunit, '(a)' ) 'grestore'
  write ( iunit, '(a)' ) '%%Trailer'
  write ( iunit, '(a,i2)' ) '%%Pages: 1'

  close ( unit = iunit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SHAPE_3D_NODES_TO_PS'
  write ( *, '(a)' ) '  The data was written to the file: ' &
    // trim ( file_name )

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
!    routine may be used to sort integers, real ( kind = 8 )s, numbers, names,
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
subroutine span_forest ( nnode, nedge, iendpt, k, component )

!*****************************************************************************80
!
!! SPAN_FOREST determines a graph's connectivity and spanning forest.
!
!  Discussion:
!
!    The input graph may be connected or unconnected.
!
!    If the input graph is connected, this routine simply returns a
!    spanning tree for the graph.
!
!    Definition: A (connected) component of a graph is a maximal subgraph
!    which is connected.
!
!    Definition: A tree is a connected graph containing no cycles.
!
!    Definition: A spanning tree of a connected graph is a subgraph which 
!    is a maximal tree.
!
!    Definition: A forest is a collection of trees, no two of which share 
!    a node.
!
!    Definition: A spanning forest of a possibly unconnected graph 
!    is a collection containing a single spanning tree for each component 
!    of the graph.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 May 1999
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
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) NEDGE, the number of edges in graph.
!
!    Input/output, integer ( kind = 4 ) IENDPT(2,NEDGE), the edge array of the graph.  
!    IENDPT(1,I) and IENDPT(2,I) are the two nodes that make up edge I.
!  
!    On input, IENDPT describes the graph.
!
!    On output, the input entries of IENDPT have been reordered, so that
!    edges belonging to the spanning forest come first, followed by those
!    edges which are not part of the spanning forest.
!
!    Output, integer ( kind = 4 ) K, the number of connected components of the graph.
!
!    Output, integer ( kind = 4 ) COMPONENT(NNODE), the component to which each node belongs.
!
  implicit none

  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) component(nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iendpt(2,nedge)
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) iq
  integer ( kind = 4 ) ir
  integer ( kind = 4 ) iret
  integer ( kind = 4 ) is
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) l0
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) loc
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m0
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) mm

  mm = 1 + max ( nnode, nedge )
 
  do i = 1, nnode
    component(i) = -i
  end do
 
  do m = 1, nedge
    do l = 1, 2
      ip = iendpt(l,m)
      iendpt(l,m) = component(ip)
      component(ip) = - l * mm - m
    end do
  end do
 
  k = 0
  loc = 0
 
10 continue
 
  do i = 1, nnode
 
    iq = component(i)
 
    if ( iq <= 0 ) then
 
      k = k + 1
      component(i) = k
 
      if ( iq + i < 0 ) then
        ip = i
        is = - iq
        iret = 31
        l1 = - iq / mm
        m1 = - iq - l1 * mm
        go to 110
      end if
 
    end if
 
  end do
 
  do m = 1, nedge
 
    do
 
      ir = - iendpt(1,m)
 
      if ( ir < 0 ) then
        exit
      end if

      call i4_swap ( iendpt(2,m), iendpt(2,ir) )

      iendpt(1,m) = iendpt(1,ir)
      iendpt(1,ir) = component(iendpt(2,ir))

    end do
 
  end do
 
  component(iendpt(2,1:loc)) = component(iendpt(1,1:loc))
 
  return
 
90    continue
 
  if ( ir /= 0 ) then
    loc = loc + 1
    component(ip) = iendpt(1,ir) + iendpt(2,ir) - ip
    iendpt(1,ir) = - loc
    iendpt(2,ir) = ip
  end if
 
  ip = m

  if ( m <= 0 ) then
    go to 10
  end if

  is = - component(ip)
 
100   continue
 
  l = is / mm
  m = is - l * mm

  if ( l == 0 ) then
    go to 90
  end if

  l1 = 3 - l
  m1 = m
 
110   continue
 
  iq = iendpt(l1,m1)
 
  if ( 0 < iq ) then
    if ( iq <= mm ) then
      if ( 0 <= component(iq) ) then
        ir = m1
      end if
    end if
  end if
 
  if ( 0 <= iq ) then
    is = abs ( iendpt(l,m) )
    iendpt(l,m) = ip
    go to 100
  end if
 
  if ( -mm <= iq ) then
 
    iq = - iq
    iendpt(l1,m1) = 0
 
    if ( iret == 31 ) then

      l0 = l1
      m0 = m1
      ir = 0
      iret = 43

    else
 
      iendpt(l0,m0) = iq
      l0 = l1
      m0 = m1
      is = abs ( iendpt(l,m) )
      iendpt(l,m) = ip

    end if

    go to 100
 
  end if
 
  iendpt(l1,m1) = - iq
  l1 = - iq / mm
  m1 = - iq - l1 * mm
  go to 110
 
end
subroutine span_tree_cand ( nedge, nnode, iarray, k, nstack, stack, &
  maxstack, iendpt, ierror, ncan )

!*****************************************************************************80
!
!! SPAN_TREE_CAND finds candidates for the K-th edge of a spanning tree of a graph.  
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
!    Input, integer ( kind = 4 ) NEDGE, the number of edges.
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) IARRAY(NNODE).  IARRAY(I) is the I-th edge of
!    the spanning tree.
!
!    Input, integer ( kind = 4 ) K, index of position in IARRAY for which
!    candidates are needed.
!
!    Output, integer ( kind = 4 ) NSTACK, the current size of the stack.
!
!    Output, integer ( kind = 4 ) STACK(MAXSTACK).  List of candidates for all positions.
!
!    Input, integer ( kind = 4 ) MAXSTACK, the maximum size of the stack.
!
!    Input, integer ( kind = 4 ) IENDPT(2,NEDGE).  IENDPT(1,I), IENDPT(2,I) are the
!    two nodes of edge I in graph.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.  0 if no errors, or 1
!    if needed stack size reached available stacksize MAXSTACK.
!    You should increase the dimension of STACK and call again.
!
!    Input/output, integer ( kind = 4 ) NCAN(NNODE-1), the number of candidates
!    for each position.
!
  implicit none

  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) nnode
  integer ( kind = 4 ) maxstack

  integer ( kind = 4 ) i
  integer ( kind = 4 ) iarray(nnode)
  integer ( kind = 4 ) iend(2,nnode)
  integer ( kind = 4 ) iendpt(2,nedge)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iwork(nnode)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) ncan(nnode-1)
  integer ( kind = 4 ) ncomp
  integer ( kind = 4 ) nstack
  integer ( kind = 4 ) stack(maxstack)

  if ( k <= 0 ) then
    ierror = 1
    return
  end if

  ncan(k) = 0

  ierror = 0
 
  if ( k == 1 ) then
 
    nstack = nedge - nnode
 
    if ( maxstack < nstack ) then
      ierror = 1
      return
    end if
 
    call i4vec_indicator ( nstack, stack )
 
    ncan(k) = nedge - nnode
 
  else

    iend(1,1:k-1) = iendpt(1,iarray(1:k-1))
    iend(2,1:k-1) = iendpt(2,iarray(1:k-1))
 
    call span_forest ( nnode, k-1, iend, ncomp, iwork )
 
    do i = iarray(k-1)+1, nedge+k+1-nnode
 
      if ( iwork(iendpt(1,i)) /= iwork(iendpt(2,i)) ) then
 
        nstack = nstack + 1
          
        if ( maxstack < nstack ) then
          ierror = 1
          return
        end if
 
        stack(nstack) = i
        ncan(k) = ncan(k) + 1

      end if
 
    end do
 
  end if
 
  return
end
subroutine span_tree_next ( signal, nnode, nedge, iendpt, iarray, ncan )

!*****************************************************************************80
!
!! SPAN_TREE_NEXT uses backtracking to find spanning forests of a graph.
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
!    Input/output, integer ( kind = 4 ) SIGNAL.
!    On input, 0 means this is the first call for a new problem.
!    On output, 0 means no more solutions exist; 1 means another solution was fo
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) NEDGE, the number of edges.
!
!    Input, integer ( kind = 4 ) IENDPT(2,NEDGE), the edge array of the graph.
!
!    Output, integer ( kind = 4 ) IARRAY(NNODE-1).  If SIGNAL = 1, then IARRAY contains
!    the "next" spanning forest found by the routine, stored as a list of
!    edge indices.
!
!    Workspace, integer NCAN(NNODE-1), the number of candidates for each 
!    position.
!
  implicit none

  integer ( kind = 4 ), parameter :: maxstack = 1000

  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) iarray(nnode-1)
  integer ( kind = 4 ) iendpt(2,nedge)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ), save :: indx
  integer ( kind = 4 ), save :: k
  integer ( kind = 4 ), dimension ( nnode-1) :: ncan
  integer ( kind = 4 ), save :: nstack
  integer ( kind = 4 ) signal
  integer ( kind = 4 ), save, dimension ( maxstack ) :: stack
!
!  First call for this problem.
!
  if ( signal == 0 ) then

    iarray(1:nnode-1) = 0
    indx = 0
    k = 0
    ncan(1:nnode-1) = 0
    nstack = 0
    stack(1:maxstack) = 0

  end if
!
!  Try to extend the current partial solution.
!
  do

    call i4vec_backtrack ( nnode-1, iarray, indx, k, nstack, stack, &
      maxstack, ncan )
!
!  A full solution was found.
!
    if ( indx == 1 ) then

      signal = 1
      exit
!
!  A partial solution was found.  Seek candidates for the next entry.
!
    else if ( indx == 2 ) then

      call span_tree_cand ( nedge, nnode, iarray, k, nstack, stack, &
        maxstack, iendpt, ierror, ncan )

      if ( ierror /= 0 ) then
        signal = 0
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SPAN_TREE_NEXT - Fatal error!'
        write ( *, '(a,i8)' ) '  Stack needs at least ', nstack
        write ( *, '(a,i8)' ) '  Available space is ', maxstack
        exit
      end if
!
!  No more found.
!
    else

      signal = 0
      exit

    end if

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
subroutine tqlrat ( n, d, e2, ierr )

!*****************************************************************************80
!
!! TQLRAT compute all eigenvalues of a real symmetric tridiagonal matrix.
!
!  Discussion:
!
!    This subroutine finds the eigenvalues of a symmetric
!    tridiagonal matrix by the rational QL method.
!
!  Reference:
!
!    Christian Reinsch,
!    Algorithm 464, TQLRAT,
!    Communications of the ACM,
!    Volume 16, page 689, 1973.
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Y Ikebe, V Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input/output, real ( kind = 8 ) D(N).  On input, D contains the diagonal
!    elements of the matrix.  On output, D contains the eigenvalues in ascending
!    order.  If an error exit was made, then the eigenvalues are correct
!    in positions 1 through IERR-1, but may not be the smallest eigenvalues.
!
!    Input/output, real ( kind = 8 ) E2(N), contains in positions 2 through N
!    the squares of the subdiagonal elements of the matrix.  E2(1) is
!    arbitrary.  On output, E2 has been overwritten by workspace
!    information.
!
!    Output, integer ( kind = 4 ) IERR, error flag.
!    0, for no error,
!    J, if the J-th eigenvalue could not be determined after 30 iterations.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) d(n)
  real ( kind = 8 ) e2(n)
  real ( kind = 8 ) f
  real ( kind = 8 ) g
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mml
  real ( kind = 8 ) p
  real ( kind = 8 ) pythag
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) t

  ierr = 0

  if ( n == 1 ) then
    return
  end if

  do i = 2, n
    e2(i-1) = e2(i)
  end do

  f = 0.0D+00
  t = 0.0D+00
  e2(n) = 0.0D+00

  do l = 1, n

     j = 0
     h = abs ( d(l) ) + sqrt ( e2(l) )

     if ( t <= h ) then

       t = h
       b = abs ( t ) * epsilon ( b )
       c = b * b

     end if
!
!  Look for small squared sub-diagonal element.
!
     do m = l, n
       if ( e2(m) <= c ) then
         exit
       end if
     end do

     if ( m == l ) then
       go to 210
     end if

130  continue

     if ( 30 <= j ) then
       ierr = l
       return
     end if

     j = j + 1
!
!  Form shift.
!
     l1 = l + 1
     s = sqrt ( e2(l) )
     g = d(l)
     p = ( d(l1) - g ) / ( 2.0D+00 * s )
     r = pythag ( p, 1.0D+00 )
     d(l) = s / ( p + sign ( r, p ) )
     h = g - d(l)

     do i = l1, n
       d(i) = d(i) - h
     end do

     f = f + h
!
!  Rational QL transformation.
!
     g = d(m)
     if ( g == 0.0D+00 ) g = b
     h = g
     s = 0.0D+00
     mml = m - l

     do ii = 1, mml
       i = m - ii
       p = g * h
       r = p + e2(i)
       e2(i+1) = s * r
       s = e2(i) / r
       d(i+1) = h + s * (h + d(i))
       g = d(i) - e2(i) / g
       if ( g == 0.0D+00 ) then
         g = b
       end if
       h = g * p / r
     end do

     e2(l) = s * g
     d(l) = h
!
!  Guard against underflow in convergence test.
!
     if ( h == 0.0D+00 ) go to 210
     if ( abs ( e2(l) ) <= abs ( c / h ) ) go to 210
     e2(l) = h * e2(l)
     if ( e2(l) /= 0.0D+00 ) go to 130

210  continue

     p = d(l) + f
!
!  Order the eigenvalues.
!
     do ii = 2, l
       i = l + 2 - ii
       if ( d(i-1) <= p ) then
         go to 270
       end if
       d(i) = d(i-1)
     end do

250  continue
     i = 1
270  continue
     d(i) = p

290  continue

  end do

  return
end
subroutine tred1 ( nm, n, a, d, e, e2 )

!*****************************************************************************80
!
!! TRED1 transforms a real symmetric matrix to tridiagonal form.
!
!  Discussion:
!
!    The routine reduces a real symmetric matrix to a symmetric
!    tridiagonal matrix using orthogonal similarity transformations.
!
!  Reference:
!
!    Martin, Reinsch, James Wilkinson,
!    TRED1,
!    Numerische Mathematik,
!    Volume 11, pages 181-195, 1968.
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer Verlag, 1971.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow, 
!    Y Ikebe, V Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NM, the leading dimension of the array A.
!    NM must be at least N.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix A.
!
!    Input/output, real ( kind = 8 ) A(NM,N), on input, contains the real
!    symmetric matrix.  Only the lower triangle of the matrix need be supplied.
!    On output, A contains information about the orthogonal transformations
!    used in the reduction in its strict lower triangle.
!    The full upper triangle of A is unaltered.
!
!    Output, real ( kind = 8 ) D(N), contains the diagonal elements of the
!    tridiagonal matrix.
!
!    Output, real ( kind = 8 ) E(N), contains the subdiagonal elements of the
!    tridiagonal matrix in its last N-1 positions.  E(1) is set to zero.
!
!    Output, real ( kind = 8 ) E2(N), contains the squares of the corresponding
!    elements of E.  E2 may coincide with E if the squares are not needed.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nm

  real ( kind = 8 ) a(nm,n)
  real ( kind = 8 ) d(n)
  real ( kind = 8 ) e(n)
  real ( kind = 8 ) e2(n)
  real ( kind = 8 ) f
  real ( kind = 8 ) g
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  real ( kind = 8 ) scale

  d(1:n) = a(n,1:n)

  do i = 1, n
    a(n,i) = a(i,i)
  end do

  do ii = 1, n

    i = n + 1 - ii
    l = i - 1
    h = 0.0D+00
    scale = 0.0D+00

    if ( l < 1 ) go to 130
!
!  Scale row.
!
    do k = 1, l
      scale = scale + abs ( d(k) )
    end do

    if ( scale /= 0.0D+00 ) go to 140

    do j = 1, l
      d(j) = a(l,j)
      a(l,j) = a(i,j)
      a(i,j) = 0.0D+00
    end do

130 continue

    e(i) = 0.0D+00
    e2(i) = 0.0D+00
    go to 300

140 continue

    do k = 1, l
      d(k) = d(k) / scale
      h = h + d(k) * d(k)
    end do

    e2(i) = scale * scale * h
    f = d(l)
    g = - sign ( sqrt ( h ), f )
    e(i) = scale * g
    h = h - f * g
    d(l) = f - g

    if ( l == 1 ) go to 285
!
!  Form A*U.
!
    e(1:l) = 0.0D+00

    do j = 1, l

      f = d(j)
      g = e(j) + a(j,j) * f

      do k = j+1, l
        g = g + a(k,j) * d(k)
        e(k) = e(k) + a(k,j) * f
      end do

      e(j) = g

    end do
!
!  Form P.
!
    f = 0.0D+00

    do j = 1, l
      e(j) = e(j) / h
      f = f + e(j) * d(j)
    end do

    h = f / ( h + h )
!
!  Form Q.
!
    e(1:l) = e(1:l) - h * d(1:l)
!
!  Form reduced A.
!
    do j = 1, l

      f = d(j)
      g = e(j)

      do k = j, l
        a(k,j) = a(k,j) - f * e(k) - g * d(k)
      end do

    end do

  285 continue

    do j = 1, l
      f = d(j)
      d(j) = a(l,j)
      a(l,j) = a(i,j)
      a(i,j) = f * scale
    end do

300 continue

  end do

  return
end
subroutine tree_arc_center ( nnode, inode, jnode, center, eccent, parity )

!*****************************************************************************
!
!! TREE_ARC_CENTER computes the center, eccentricity, and parity of a tree.
!
!  Definition:
!
!    The edge distance between two nodes I and J is the minimum number of
!    edges that must be traversed in a path from I and J.
!
!    The eccentricity of a node I is the maximum edge distance between
!    node I and the other nodes J in the graph.
!
!    The radius of a graph is the minimum eccentricity over all nodes
!    in the graph.
!
!    The diameter of a graph is the maximum eccentricity over all nodes
!    in the graph.
!
!    The center of a graph is the set of nodes whose eccentricity is 
!    equal to the radius, that is, the set of nodes of minimum eccentricity.
!
!    For a tree, the center is either a single node, or a pair of
!    neighbor nodes.
!
!    The parity of the tree is 1 if the center is a single node, or 2 if
!    the center is 2 nodes.
!
!  Discussion:
!
!    The center of a tree can be found by removing all "leaves", that is,
!    nodes of degree 1.  This step is repeated until only 1 or 2 nodes
!    are left.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) INODE(NNODE-1), JNODE(NNODE-1), the edges of the
!    tree.  Edge I connects nodes INODE(I) and JNODE(I).
!
!    Output, integer ( kind = 4 ) CENTER(2).  CENTER(1) is the index of the first
!    node in the center.  CENTER(2) is 0 if there is only one node
!    in the center, or else the index of the second node.
!
!    Output, integer ( kind = 4 ) ECCENT, the eccentricity of the nodes in the center,
!    and the radius of the the tree.
!
!    Output, integer ( kind = 4 ) PARITY, the parity of the tree, which is normally
!    1 or 2.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) center(2)
  integer ( kind = 4 ) degree(nnode)
  integer ( kind = 4 ) eccent
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iedge
  integer ( kind = 4 ) ileaf
  integer ( kind = 4 ) inode(nnode-1)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jnode(nnode-1)
  integer ( kind = 4 ) list(nnode)
  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) nleaf
  integer ( kind = 4 ) nnode2
  integer ( kind = 4 ) parity

  eccent = 0
  center(1) = 0
  center(2) = 0
  parity = 0

  if ( nnode <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TREE_ARC_CENTER - Fatal error!'
    write ( *, '(a)' ) '  NNODE <= 0.'
    stop
  else if ( nnode == 1 ) then
    eccent = 0
    center(1) = 1
    center(2) = 0
    parity = 1
    return
  else if ( nnode == 2 ) then
    eccent = 1
    center(1) = 1
    center(2) = 2
    parity = 2
    return
  end if
!
!  Compute the degrees.
!
  nedge = nnode - 1
  call graph_arc_degree ( nnode, nedge, inode, jnode, degree )
!
!  Defoliate the tree.
!
  nnode2 = nnode

  do

    eccent = eccent + 1
!
!  Find and mark the leaves.
!
    nleaf = 0

    do i = 1, nnode

      if ( degree(i) == 1 ) then
        nleaf = nleaf + 1
        list(nleaf) = i
      end if

    end do
!
!  Delete the leaves.
!
    do ileaf = 1, nleaf

      i = list(ileaf)

      iedge = 0
      j = 0

      do

        iedge = iedge + 1

        if ( nedge < iedge ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'TREE_ARC_CENTER - Fatal error!'
          write ( *, '(a)' ) '  Data or algorithm failure.'
          stop
        end if

        if ( inode(iedge) .eq. i ) then
          j = jnode(iedge)
          inode(iedge) = - inode(iedge)
          jnode(iedge) = - jnode(iedge)
        else if ( jnode(iedge) == i ) then
          j = inode(iedge)
          inode(iedge) = - inode(iedge)
          jnode(iedge) = - jnode(iedge)
        end if

        if ( j /= 0 ) then
          exit
        end if

      end do

      degree(i) = 0
      nnode2 = nnode2 - 1
      degree(j) = degree(j) - 1
      if ( degree(j) == 0 ) then
        nnode2 = nnode2 - 1
      end if

    end do
!
!  Find the remaining nodes.
!
    nnode2 = 0

    do i = 1, nnode

      if ( 0 < degree(i) ) then
        nnode2 = nnode2 + 1
        list(nnode2) = i
      end if

    end do
!
!  If at least 3, more pruning is required.
!
    if ( nnode2 < 3 ) then
      exit
    end if

  end do
!
!  If only one or two nodes left, we are done.
!
  parity = nnode2

  center(1:nnode2) = list(1:nnode2)
  inode(1:nedge) = abs ( inode(1:nedge) )
  jnode(1:nedge) = abs ( jnode(1:nedge) )

  return
end
subroutine tree_arc_diam ( nnode, inode, jnode, diam, label, n1, n2 )

!*****************************************************************************
!
!! TREE_ARC_DIAM computes the "diameter" of a tree.
!
!  Definition:
!
!    The diameter of a graph is the length of the longest possible
!    path that never repeats an edge.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) INODE(NNODE-1), JNODE(NNODE-1), the edges of the
!    tree.  Edge I connects nodes INODE(I) and JNODE(I).
!
!    Output, integer ( kind = 4 ) DIAM, the length of the longest path in the tree.
!
!    Output, integer ( kind = 4 ) LABEL(NNODE), marks the path between nodes N1 and N2.
!    Node I is in this path if LABEL(I) is 1.
!
!    Output, integer ( kind = 4 ) N1, N2, the indices of two nodes in the tree which
!    are separated by DIAM edges.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) degree(nnode)
  integer ( kind = 4 ) diam
  integer ( kind = 4 ) i
  integer ( kind = 4 ) inode(nnode-1)
  integer ( kind = 4 ) invals
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jnode(nnode-1)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kstep
  integer ( kind = 4 ) label(nnode)
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) nabe
  integer ( kind = 4 ) nedge

  if ( nnode <= 0 ) then
    diam = 0
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TREE_ARC_DIAM - Fatal error!'
    write ( *, '(a)' ) '  NNODE <= 0.'
    stop
  end if

  if ( nnode == 1 ) then
    diam = 0
    return
  end if

  nedge = nnode - 1

  call i4vec_indicator ( nnode, label )
!
!  On step K:
!
!    Identify the terminal and interior nodes.
!
!    If there are no interior nodes left, 
!
!      then there are just two nodes left at all.  The diameter is 2*K-1, 
!      and a maximal path extends between the nodes whose labels are 
!      contained in the two remaining terminal nodes.
!
!    Else
!
!      The label of each terminal node is passed to its interior neighbor.
!      If more than one label arrives, take any one.
!
!      The terminal nodes are removed.
!
  kstep = 0

  do

    kstep = kstep + 1
!
!  Compute the degree of each node.
!
    degree(1:nnode) = 0        

    do j = 1, nedge

      k = inode(j)
      if ( 0 < k ) then
        degree(k) = degree(k) + 1
      end if

      k = jnode(j)
      if ( 0 < k ) then
        degree(k) = degree(k) + 1
      end if

    end do
!
!  Count the number of interior nodes.
!
    invals = 0
    do i = 1, nnode
      if ( 1 < degree(i) ) then
        invals = invals + 1
      end if
    end do
!
!  If there are 1 or 0 interior nodes, it's time to stop.
!
    if ( invals == 1 ) then

      diam = 2 * kstep
      exit
    
    else if ( invals == 0 ) then

      diam = 2 * kstep - 1
      exit

    end if
!
!  If there are at least two interior nodes, then chop off the 
!  terminal nodes and pass their labels inward.
!
    do k = 1, nnode
      if ( degree(k) == 1 ) then

        do j = 1, nedge

          if ( inode(j) == k ) then
            nabe = jnode(j)
            label(nabe) = label(k)
            inode(j) = - inode(j)
            jnode(j) = - jnode(j)
          else if ( jnode(j) == k ) then
            nabe = inode(j)
            label(nabe) = label(k)
            inode(j) = - inode(j)
            jnode(j) = - jnode(j)
          end if

        end do

      end if

    end do

  end do
!
!  Now get the labels from two of the remaining terminal nodes.
!  The nodes represented by these labels will be a diameter apart.
!
  n1 = 0
  n2 = 0

  do i = 1, nnode
    if ( degree(i) == 1 ) then
      if ( n1 == 0 ) then
        n1 = label(i)
      else if ( n2 == 0 ) then
        n2 = label(i)
      end if
    end if
  end do
!
!  Set the labels of the interior node (if any) and nodes marked
!  N1 and N2 to 1, and all others to 0.  This will label the nodes on the path.
!
  if ( invals == 1 ) then

    do i = 1, nnode
      if ( 1 < degree(i) ) then
        label(i) = 1
      end if
    end do

  end if

  do i = 1, nnode
    if ( label(i) == n1 .or. label(i) == n2 ) then
      label(i) = 1
    else
      label(i) = 0
    end if
  end do
!
!  Clean up the arrays.
!
  do j = 1, nedge
    inode(j) = abs ( inode(j) )
    k = inode(j)
    jnode(j) = abs ( jnode(j) )
    k = jnode(j)
  end do

  return
end
subroutine tree_arc_random ( nnode, seed, code, inode, jnode )

!*****************************************************************************
!
!! TREE_ARC_RANDOM selects a random labeled tree and its Pruefer code.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 March 2005
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
!    Output, integer ( kind = 4 ) CODE(NNODE-2), the Pruefer code for the labeled tree.
!
!    Output, integer ( kind = 4 ) INODE(NNODE-1), JNODE(NNODE-1), the edge array for
!    the tree.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) code(nnode-2)
  integer ( kind = 4 ) inode(nnode-1)
  integer ( kind = 4 ) jnode(nnode-1)
  integer ( kind = 4 ) seed

  if ( nnode <= 0  ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TREE_ARC_RANDOM - Fatal error!'
    write ( *, '(a,i8)' ) '  NNODE = ', nnode
    write ( *, '(a)' ) '  but NNODE must be at least 1.'
    stop
  end if

  if ( nnode <= 2 ) then
    return
  end if

  call vec_random ( nnode-2, nnode, seed, code )
 
  code(1:nnode-2) = code(1:nnode-2) + 1
 
  call pruefer_to_tree_arc ( nnode, code, inode, jnode )
 
  return
end
subroutine tree_arc_to_pruefer ( nnode, inode, jnode, code )

!*****************************************************************************80
!
!! TREE_ARC_TO_PRUEFER is given a labeled tree, and computes its Pruefer code.
!
!  Discussion:
!
!    The Pruefer code is a correspondence between all labeled trees of
!    N nodes, and all list of N-2 integers between 1 and N (with repetition
!    allowed).  The number of labeled trees on N nodes is therefore N**(N-2).
!
!  Method:
!
!    The Pruefer code is constructed from the tree as follows:
!
!    A terminal node on the tree is defined as a node with only one neighbor.
!
!    Consider the set of all terminal nodes on the tree.  Take the one
!    with the highest label, I.  Record the label of its neighbor, J.
!    Delete node I and the edge between node I and J.
!
!    J is the first entry in the Pruefer code for the tree.   Repeat
!    the operation a total of N-2 times to get the complete code.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 August 2000
!
!  Reference:
!
!    Dennis Stanton, Dennis White,
!    Constructive Combinatorics,
!    Springer Verlage, New York, 1986.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) INODE(NNODE-1), JNODE(NNODE-1), the edge array of the
!    tree.  The I-th edge joins nodes INODE(I) and JNODE(I).
!
!    Output, integer ( kind = 4 ) CODE(NNODE-2), the Pruefer code of the tree.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) code(nnode-2)
  integer ( kind = 4 ) degree(nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) iterm
  integer ( kind = 4 ) inode(nnode-1)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jnode(nnode-1)
  integer ( kind = 4 ) jsave
  integer ( kind = 4 ) nedge
!
!  Compute the degree of each node.
!
  nedge = nnode - 1
  call graph_arc_degree ( nnode, nedge, inode, jnode, degree )
!
!  Compute the next term of the Pruefer code.
!
  do i = 1, nnode-2
!
!  Find the terminal node with the highest label.
!
    iterm = 0
 
    do j = 1, nnode
      if ( degree(j) == 1 ) then
        iterm = j
      end if
    end do
!
!  Find the edge that includes this node, and note the
!  index of the other node.
!
    do j = 1, nnode-1

      jsave = j
 
      if ( inode(j) == iterm ) then
        i2 = 2
        exit
      else if ( jnode(j) == iterm ) then
        i2 = 1
        exit
      end if
 
    end do
!
!  Delete the edge from the tree.
!
    degree(inode(jsave)) = degree(inode(jsave)) - 1
    degree(jnode(jsave)) = degree(jnode(jsave)) - 1
!
!  Add the neighbor of the node to the Pruefer code.
!
    if ( i2 == 1 ) then
      code(i) = inode(jsave)
    else
      code(i) = jnode(jsave)
    end if
!
!  Negate the nodes in the edge list to mark them as used.
!
    inode(jsave) = - inode(jsave)
    jnode(jsave) = - jnode(jsave)
 
  end do
!
!  Before returning, restore the original form of the edge list.
!
  inode(1:nnode-1) = abs ( inode(1:nnode-1) )
  jnode(1:nnode-1) = abs ( jnode(1:nnode-1) )
 
  return
end
subroutine tree_enum ( nnode, ntree )

!*****************************************************************************80
!
!! TREE_ENUM enumerates the labeled trees on NNODE nodes.
!
!  Example:
!
!    NNODE      NTREE
!
!    0              1
!    1              1
!    2              1
!    3              3
!    4             16
!    5            125
!    6           1296
!    7          16807
!    8         262144
!    9        4782969
!   10      100000000
!
!  Discussion:
!
!    The formula is due to Cauchy.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes in each tree.
!    NNODE must normally be at least 3, but for this routine,
!    any value of NNODE is allowed.  Values of NNODE greater than 10
!    will probably overflow.
!
!    Output, integer ( kind = 4 ) NTREE, the number of distinct labeled trees.
!
  implicit none

  integer ( kind = 4 ) nnode
  integer ( kind = 4 ) ntree

  if ( nnode < 0 ) then
    ntree = 0
  else if ( nnode == 0 ) then
    ntree = 1
  else if ( nnode == 1 ) then
    ntree = 1
  else if ( nnode == 2 ) then
    ntree = 1
  else
    ntree = nnode**( nnode - 2 )
  end if

  return
end
subroutine tree_parent_next ( nnode, iarray, code, itree, more )

!*****************************************************************************
!
!! TREE_PARENT_NEXT generates, one at a time, all labeled trees.
!
!  Discussion:
!
!    The routine also returns the corresponding Pruefer codes.
!
!  Formula:
!
!    There are N**(N-2) labeled trees on N nodes (Cayley's formula).
!
!    The number of trees in which node I has degree D(I) is the
!    multinomial coefficient: ( N-2; D(1)-1, D(2)-1, ..., D(N)-1 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 August 2000
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes to be used in the trees.
!
!    Workspace, integer IARRAY(NNODE).
!
!    Output, integer ( kind = 4 ) CODE(NNODE).  The first NNODE-2 entries of CODE
!    contain the Pruefer code for the given labeled tree.
!
!    Output, integer ( kind = 4 ) ITREE(NNODE).  The first NNODE-1 entries of ITREE
!    describe the edges that go between the nodes.  Each pair
!    (I, ITREE(I)) represents an edge.  Thus if ITREE(5) = 3,
!    there is an edge from node 3 to node 5.
!
!    Input/output, logical MORE.  On the first call only, the
!    user is required to set MORE = .FALSE.  Then call TRENEX, and
!    the program will return information about the first tree
!    as well as setting MORE to the value .TRUE.
!    Keep calling to get another tree until MORE is .FALSE.
!    on return, at which point there are no more trees.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) code(nnode)
  integer ( kind = 4 ), dimension ( nnode ) :: iarray
  integer ( kind = 4 ) itree(nnode)
  logical more

  call vec_next ( nnode-2, iarray, more, nnode )
 
  code(1:nnode-2) = iarray(1:nnode-2) + 1
 
  call pruefer_to_tree_2 ( nnode, code, itree )
 
  return
end
subroutine tree_parent_to_arc ( nnode, parent, nedge, inode, jnode )

!*****************************************************************************80
!
!! TREE_PARENT_TO_ARC converts a tree from parent to arc representation.
!
!  Discussion:
!
!    Parent representation lists the parent node of each node.  For a
!    tree of N nodes, one node has a parent of 0, representing a null link.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes in the tree.
!
!    Input, integer ( kind = 4 ) PARENT(N), the parent node representation of the tree.
!
!    Output, integer ( kind = 4 ) NEDGE, the number of edges, normally NNODE-1.
!
!    Output, integer ( kind = 4 ) INODE(NEDGE), JNODE(NEDGE), pairs of nodes that
!    define the links.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) i
  integer ( kind = 4 ) inode(nnode-1)
  integer ( kind = 4 ) jnode(nnode-1)
  integer ( kind = 4 ) nedge
  integer ( kind = 4 ) parent(nnode)

  nedge = 0

  do i = 1, nnode

    if ( parent(i) /= 0 ) then
      nedge = nedge + 1
      inode(nedge) = i
      jnode(nedge) = parent(i)
    end if

  end do

  return
end
subroutine tree_rb_enum ( n, num )

!*****************************************************************************80
!
!! TREE_RB_ENUM returns the number of rooted binary trees with N nodes.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of nodes in the rooted binary tree.
!    N should be odd.
!
!    Output, integer ( kind = 4 ) NUM, the number of rooted binary trees with N nodes.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) c(0:(n-1)/2)
  integer ( kind = 4 ) num

  if ( n < 0 ) then

    num = 0

  else if ( n == 0 ) then

    num = 1

  else if ( mod ( n, 2 ) == 0 ) then

    num = 0

  else

    call catalan ( ( n - 1 ) / 2, c )

    num = c((n-1)/2)

  end if

  return
end
subroutine tree_rb_lex_next ( n, a, more )

!*****************************************************************************80
!
!! TREE_RB_LEX_NEXT generates rooted binary trees in lexicographic order.
!
!  Discussion:
!
!    The information definining the tree of N nodes is stored in a vector 
!    of 0's and 1's, in preorder traversal form.  Essentially, the
!    shape of the tree is traced out with a pencil that starts at the root,
!    and stops at the very last null leaf.  The first time that a (non-null) 
!    node is encountered, a 1 is added to the vector, and the left 
!    descendant of the node is visited next.  When the path returns from
!    the first descendant, the second descendant is visited.  When then path
!    returns again, the path passes back up from the node to its parent.
!    A null leaf is encountered only once, and causes a zero to be added to 
!    the vector, and the path goes back up to the parent node.  
!
!    The lexicographic order is used to order the vectors of 1's and 0's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 August 2000
!
!  Reference:
!
!    Frank Ruskey,
!    Combinatorial Generation,
!    To appear.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of nodes in the rooted binary tree.
!    N should be odd.
!
!    Input/output, integer ( kind = 4 ) A(N), the preorder traversal form for the 
!    previous/next rooted binary tree.
!
!    Output, logical MORE, is TRUE if the next rooted binary tree was
!    returned on this call, or FALSE if there are no more rooted binary
!    trees, and the output of the previous call was the last one.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) k
  logical more
  integer ( kind = 4 ) p
  integer ( kind = 4 ) q

  if ( .not. more ) then
    a(1:n-2:2) = 1
    a(2:n-1:2) = 0
    a(n) = 0
    more = .true.
    return
  end if
!
!  Find the last 1 in A.
!
  k = n
  do while ( a(k) == 0 )
    k = k - 1
  end do
  q = n - k - 1
!
!  Find the last 0 preceding the last 1 in A.
!  If there is none, then we are done, because 11...1100..00 
!  is the final element.
!
  do 

    if ( k == 1 ) then
      more = .false.
      return
    end if

    if ( a(k) == 0 ) then
      exit
    end if

    k = k - 1

  end do
	
  p = n - k - q - 1
  a(k) = 1
  a(k+1:n-2*p+1) = 0
  a(n-2*p+2:n-2:2) = 1
  a(n-2*p+3:n-1:2) = 0
  a(n) = 0

  return
end
subroutine tree_rb_to_parent ( n, a, parent )

!*****************************************************************************80
!
!! TREE_RB_TO_PARENT converts a rooted binary tree to parent node representation.
!
!  Discussion:
!
!    Parent node representation of a tree assigns to each node a "parent" node,
!    which represents the first link of the path between the node and the 
!    root node.  The root node itself is assigned a parent of 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of nodes in the tree.
!
!    Input, integer ( kind = 4 ) A(N), the preorder traversal form for the rooted 
!    binary tree.
!
!    Output, integer ( kind = 4 ) PARENT(N), the parent node representation of the tree.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) dad
  integer ( kind = 4 ) k
  integer ( kind = 4 ) node
  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) parent(n)
  integer ( kind = 4 ) use(n)

  node = 0
  node_num = 0

  do k = 1, n

    dad = node
    node_num = node_num + 1
    node = node_num
    parent(node) = dad

    if ( a(k) == 1 ) then

      use(node) = 0

    else

      use(node) = 2

      do while ( use(node) == 2 )
        node = dad
        if ( node == 0 ) then
          exit
        end if
        use(node) = use(node) + 1
        dad = parent(node)
      end do

    end if

  end do

  return
end
subroutine tree_rb_yule ( n, seed, a )

!*****************************************************************************80
!
!! TREE_RB_YULE adds two nodes to a rooted binary tree using the Yule model.
!
!  Discussion:
!
!    The Yule model is a simulation of how an evolutionary family tree
!    develops.  We start with a root node.  The internal nodes of the tree 
!    are inactive and never change.  Each pendant or leaf node of the
!    tree represents a biological family that can spontaneously "fission",
!    developing two new distinct sub families.  In graphical terms, the node
!    becomes internal, with two new leaf nodes depending from it.
!
!    The tree is stored in inorder traversal form.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N, the number of nodes in the input tree.
!    On output, this number has been increased, usually by 2.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
!    Input/output, integer ( kind = 4 ) A(*), the preorder traversal form for the rooted 
!    binary tree.  The number of entries in A is N.
!
  implicit none

  integer ( kind = 4 ) a(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) ileaf
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jleaf
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nleaf
  integer ( kind = 4 ) seed

  if ( n <= 0 ) then
    n = 1
    a(1) = 0
    return
  end if
!
!  Count the expected number of leaves, which are the 0 values.
!
  nleaf = ( n + 1 ) / 2
!
!  Choose a random number between 1 and NLEAF.
!
  ileaf = i4_uniform ( 1, nleaf, seed )
!
!  Locate leaf number ILEAF.
!
  j = 0
  jleaf = 0
  do i = 1, n
    if ( a(i) == 0 ) then
      jleaf = jleaf + 1
    end if
    if ( jleaf == ileaf ) then
      j = i
      exit
    end if
  end do
!
!  Replace '0' by '100'
!
  a(n+2:j+2:-1) = a(n:j:-1)
  a(j:j+1) = (/ 1, 0 /)

  n = n + 2

  return
end
subroutine tree_rooted_code ( nnode, parent, code )

!*****************************************************************************80
!
!! TREE_ROOTED_CODE returns the code of a rooted tree.
!
!  Discussion:
!
!    This code for a rooted tree depends on the node ordering, so it's actually
!    the code for a labeled rooted tree.  To eliminate the effects of node
!    labeling, one could choose as the code for a tree the maximum of all
!    the codes associated with the different possible labelings of the tree.
!    There are more effective ways of arriving at this code than simply
!    generating all possible codes and comparing them.  
!
!    For a tree with NNODES, the code is a list of 2*NNODE 0's and 1's,
!    describing a traversal of the tree starting at an imaginary node 0,
!    moving "down" to the root (a code entry of 1), and then moving
!    "down" (1) or "up" (0) as the tree is traversed in a depth first
!    manner.  The final move must be from the root up to the imaginary
!    node 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) PARENT(NNODE), is the parent node of each node.
!    The node with parent 0 is the root.
!
!    Output, integer ( kind = 4 ) CODE(2*NNODE), the code for the tree.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) code(2*nnode)
  integer ( kind = 4 ) father
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) parent(nnode)
  integer ( kind = 4 ) son
!
!  Find the root.
!
  father = 0
  do i = 1, nnode
    if ( parent(i) == 0 ) then
      k = 1
      code(1) = 1
      father = i
      exit
    end if
  end do

  if ( father == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TREE_ROOTED_CODE - Fatal error!'
    write ( *, '(a)' ) '  Could not find the root.'
    stop
  end if

  do while ( father /= 0 ) 

    k = k + 1
    code(k) = 0

    do son = 1, nnode
      if ( parent(son) == father ) then
        code(k) = 1
        father = son
        exit
      end if
    end do

    if ( code(k) == 0 ) then
      parent(father) = - parent(father)
      father = - parent(father)
    end if

  end do

  parent(1:nnode) = - parent(1:nnode)

  return
end
subroutine tree_rooted_code_compare ( nnode, npart, code1, code2, result )

!*****************************************************************************80
!
!! TREE_ROOTED_CODE_COMPARE compares (a portion of) the codes for two rooted trees.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 September 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) NPART, the number of nodes for which the code has
!    been determined.  This determines the portion of the codes to be compared.
!    We expect 0 <= NPART <= NNODE.
!
!    Input, integer ( kind = 4 ) CODE1(2*NNODE), CODE2(2*NNODE), the two rooted tree codes
!    to be compared.
!
!    Output, integer ( kind = 4 ) RESULT, the result of the comparison.
!    -1, CODE1 < CODE2,
!     0, CODE1 = CODE2,
!    +1, CODE1 > CODE2.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) code1(2*nnode)
  integer ( kind = 4 ) code2(2*nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) npart
  integer ( kind = 4 ) result

  result = 0

  if ( npart <= 0 ) then
    return
  end if

  ihi = 2 * min ( npart, nnode )

  do i = 1, ihi

    if ( code1(i) < code2(i) ) then
      result = -1
      return
    else if ( code2(i) < code1(i) ) then
      result = +1
      return
    end if

  end do

  return
end
subroutine tree_rooted_depth ( nnode, parent, depth, depth_node )

!*****************************************************************************80
!
!! TREE_ROOTED_DEPTH returns the depth of a rooted tree.
!
!  Definition:
!
!    The depth of any node of a rooted tree is the number of edges in 
!    the shortest path from the root to the node.
!
!    The depth of the rooted tree is the maximum of the depths
!    of all the nodes.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input, integer ( kind = 4 ) PARENT(NNODE), is the parent node of each node.
!    The node with parent 0 is the root.
!
!    Output, integer ( kind = 4 ) DEPTH, the depth of the tree.
!
!    Output, integer ( kind = 4 ) DEPTH_NODE(NNODE), the depth of each node.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) depth
  integer ( kind = 4 ) depth_node(nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) parent(nnode)
  integer ( kind = 4 ) root
!
!  Find the root.
!
  root = 0
  do i = 1, nnode
    if ( parent(i) == 0 ) then
      root = i
      exit
    end if
  end do

  if ( root == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TREE_ROOTED_DEPTH - Fatal error!'
    write ( *, '(a)' ) '  Could not find the root.'
    stop
  end if
!
!  Determine the depth of each node by moving towards the node.
!  If you reach a node whose depth is already known, stop early.
!
  depth_node(1:nnode) = 0

  do i = 1, nnode

    j = i

    do while ( j /= root )

      depth_node(i) = depth_node(i) + 1
      j = parent(j)

      if ( 0 < depth_node(j) ) then
        depth_node(i) = depth_node(i) + depth_node(j)
        exit
      end if

    end do

  end do
!
!  Determine the maximum depth.
!
  depth = maxval ( depth_node(1:nnode) )

  return
end
subroutine tree_rooted_enum ( nnode, ntree )

!*****************************************************************************80
!
!! TREE_ROOTED_ENUM counts the number of unlabeled rooted trees.
!
!  Example:
!
!    Input    Output
!
!      1         1
!      2         1
!      3         2
!      4         4
!      5         9
!      6        20
!      7        48
!      8       115
!      9       286
!     10       719
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 August 2000
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
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) NTREE(NNODE).  NTREE(I) is the number of rooted,
!    unlabeled trees on I nodes, for I = 1, 2, ... NNODE.
!
  implicit none

  integer ( kind = 4 ) nnode

  integer ( kind = 4 ) i
  integer ( kind = 4 ) id
  integer ( kind = 4 ) isum
  integer ( kind = 4 ) itd
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nlast
  integer ( kind = 4 ) ntree(nnode)

  ntree(1) = 1
 
  do nlast = 2, nnode
 
    isum = 0
 
    do id = 1, nlast-1
 
      i = nlast
      itd = ntree(id) * id
 
      do j = 1, nlast-1

        i = i - id

        if ( i <= 0 ) then
          exit
        end if

        isum = isum + ntree(i) * itd

      end do
 
    end do
 
    ntree(nlast) = isum / ( nlast - 1 )
 
  end do

  return
end
subroutine tree_rooted_random ( nnode, seed, ntree, itree )

!*****************************************************************************80
!
!! TREE_ROOTED_RANDOM selects a random unlabeled rooted tree.
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
!    Input, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, integer ( kind = 4 ) NTREE(NNODE).  NTREE(I) is the number of rooted,
!    unlabeled trees on I nodes, for I = 1, 2, ... NNODE.
!
!    Output, integer ( kind = 4 ) ITREE(NNODE).  (I,ITREE(I)) is the I-th edge of the
!    output tree for I = 2,NNODE.  ITREE(1)=0.
!
  implicit none

  integer ( kind = 4 ) nnode

  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) i
  integer ( kind = 4 ) id
  integer ( kind = 4 ) is1
  integer ( kind = 4 ) is2
  integer ( kind = 4 ) itd
  integer ( kind = 4 ) itree(nnode)
  integer ( kind = 4 ) iz
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l
  integer ( kind = 4 ) ll
  integer ( kind = 4 ) ls
  integer ( kind = 4 ) m
  integer ( kind = 4 ) ntree(nnode)
  integer ( kind = 4 ) nval
  real ( kind = 8 ) r
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) stack(2,nnode)

  if ( nnode <= 0  ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TREE_ROOTED_RANDOM - Fatal error!'
    write ( *, '(a,i8)' ) '  NNODE = ', nnode
    write ( *, '(a)' ) '  but NNODE must be at least 1.'
    stop
  end if
!
!  Compute a table of the number of such trees for a given number of nodes.
!
  call tree_rooted_enum ( nnode, ntree )
!
!  Now select one such tree at random.
!
  l = 0

  nval = nnode
  is1 = 0
  is2 = 0
  
10  continue

  do while ( 2 < nval )
 
    r = r8_uniform_01 ( seed )

    iz = int ( ( nval - 1 ) * ntree(nval) * r )

    id = 0
 
20  continue
 
    id = id + 1
    itd = id * ntree(id)
    m = nval
    j = 0
 
30  continue
 
    j = j + 1
    m = m - id

    if ( m < 1 ) then
      go to 20
    end if

    iz = iz - ntree(m) * itd

    if ( 0 <= iz ) then
      go to 30
    end if

    is1 = is1 + 1
    stack(1,is1) = j
    stack(2,is1) = id
    nval = m
 
  end do
 
  itree(is2+1) = l
  l = is2 + 1
  is2 = is2 + nval

  if ( 1 < nval ) then
    itree(is2) = is2 - 1
  end if
 
  do
 
    nval = stack(2,is1)
 
    if ( nval /= 0 ) then
      stack(2,is1) = 0
      go to 10
    end if
 
    j = stack(1,is1)
    is1 = is1 - 1
    m = is2 - l + 1
    ll = itree(l)
    ls = l + ( j - 1 ) * m - 1
 
    if ( j /= 1 ) then
      do i = l, ls
        itree(i+m) = itree(i) + m
        if ( mod(i-l,m) == 0 ) then
          itree(i+m) = ll
        end if
      end do
    end if
 
    is2 = ls + m
 
    if ( is2 == nnode ) then
      exit
    end if

    l = ll
 
  end do
 
  return
end
subroutine vec_next ( n, iarray, more, ibase )

!*****************************************************************************80
!
!! VEC_NEXT generates all N-vectors of integers modulo a given base.
!
!  Discussion:
!
!    The items are produced one at a time.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 April 1999
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the vectors to be used.
!
!    Output, integer ( kind = 4 ) IARRAY(N).  On each return from VECNEX, IARRAY
!    will contain entries in the range 0 to IBASE-1.
!
!    Input/output, logical MORE.  Set this variable .FALSE. before
!    the first call.  Normally, MORE will be returned .TRUE. but
!    once all the vectors have been generated, MORE will be
!    reset .FALSE. and you should stop calling the program.
!
!    Input, integer ( kind = 4 ) IBASE, the base to be used.  IBASE = 2 will
!    give vectors of 0's and 1's, for instance.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) iarray(n)
  integer ( kind = 4 ) ibase
  integer ( kind = 4 ), save :: kount
  integer ( kind = 4 ), save :: last
  logical more
  integer ( kind = 4 ) nn

  if ( .not. more ) then
 
    kount = 1
    last = ibase**n
    more = .true.
    iarray(1:n) = 0
 
  else
 
    kount = kount + 1

    if ( kount == last ) then
      more = .false.
    end if

    iarray(n) = iarray(n) + 1
 
    do i = 1, n

      nn = n - i

      if ( iarray(nn+1) < ibase ) then
        return
      end if

      iarray(nn+1) = 0

      if ( nn /= 0 ) then
        iarray(nn) = iarray(nn) + 1
      end if

    end do
 
  end if
 
  return
end
subroutine vec_random ( n, base, seed, iarray )

!*****************************************************************************
!
!! VEC_RANDOM selects a random N-vector of integers modulo a given base.
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
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the vector to be generated.
!
!    Input, integer ( kind = 4 ) BASE, the base to be used.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, integer ( kind = 4 ) IARRAY(N), a list of N random values between
!    0 and IBASE-1.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) base
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) iarray(n)
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) seed

  do i = 1, n
    ival = i4_uniform ( 0, base-1, seed )
    iarray(i) = ival
  end do
 
  return
end
subroutine vla_to_graph_arc ( file_name, maxedge, maxnode, nedge, nnode, &
  inode, jnode, x, y, z, ierror )

!*****************************************************************************80
!
!! VLA_TO_GRAPH_ARC reads graphics information from a VLA file.
!
!  Discussion:
!
!    Internal comments begin with a semi4colon in column 1.
!
!    The X, Y, Z coordinates of points begin with a "P" to
!    denote the beginning of a line, and "L" to denote the
!    continuation of a line.  The fourth entry is intensity, which
!    should be between 0.0 and 1.0.
!
!    It is intended that the information read from the file can
!    either start a whole new graphics object, or simply be added
!    to a current graphics object via the '<<' command.
!
!    This is controlled by whether the input values have been zeroed
!    out or not.  This routine simply tacks on the information it
!    finds to the current graphics object.
!
!  Example:
!
!     set comment fish.vla created by IVREAD
!     set comment from data in file fish.iv
!     set comment
!     set intensity EXPLICIT
!     set parametric NON_PARAMETRIC
!     set filecontent LINES
!     set filetype NEW
!     set depthcue 0
!     set defaultdraw stellar
!     set coordsys RIGHT
!     set author IVREAD
!     set site Buhl Planetarium
!     set library_id UNKNOWN
!     ; DXF LINE entity
!     P   8.59816       5.55317      -3.05561       1.00000
!     L   8.59816       2.49756      0.000000D+00   1.00000
!     L   8.59816       2.49756      -3.05561       1.00000
!     L   8.59816       5.55317      -3.05561       1.00000
!     ; DXF LINE entity
!     P   8.59816       5.55317      0.000000D+00   1.00000
!     ...etc...
!     L   2.48695       2.49756      -3.05561       1.00000
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the file name.
!
!    Input, integer ( kind = 4 ) MAXEDGE, the maximum number of edges.
!
!    Input, integer ( kind = 4 ) MAXNODE, the maximum number of nodes.
!
!    Output, integer ( kind = 4 ) NEDGE, the number of edges.
!
!    Output, integer ( kind = 4 ) NNODE, the number of nodes.
!
!    Output, integer ( kind = 4 ) INODE(MAXEDGE), JNODE(MAXEDGE), node pairs 
!    of each edge.
!
!    Output, real ( kind = 8 ) X(MAXNODE), Y(MAXNODE), Z(MAXNODE),
!    the coordinates of nodes.
!
!    Output, integer ( kind = 4 ) IERROR, 0 no error, 1 an error.
!
  implicit none

  integer   ( kind = 4 )  maxnode
  integer   ( kind = 4 )  maxedge

  logical                 done
  character ( len = * )   file_name
  integer   ( kind = 4 )  i
  integer   ( kind = 4 )  icor3
  integer   ( kind = 4 )  icor3_old
  integer   ( kind = 4 )  iedge
  integer   ( kind = 4 )  ierror
  integer   ( kind = 4 )  indx(maxnode)
  integer   ( kind = 4 )  indx_edge(maxedge)
  integer   ( kind = 4 )  inode(maxedge)
  integer   ( kind = 4 )  ios
  integer   ( kind = 4 )  iunit
  integer   ( kind = 4 )  iword
  integer   ( kind = 4 )  jcor3
  integer   ( kind = 4 )  jcor3_old
  integer   ( kind = 4 )  jnode(maxedge)
  integer   ( kind = 4 )  lchar
  integer   ( kind = 4 )  nedge
  integer   ( kind = 4 )  num_bad
  integer   ( kind = 4 )  num_bad_old
  integer   ( kind = 4 )  nnode
  real         ( kind = 8 )  rval
  logical                 s_eqi
  character ( len = 256 ) text
  character ( len = 256 ) word
  character ( len = 256 ) word1
  real         ( kind = 8 )  x(maxnode)
  real         ( kind = 8 )  xval
  real         ( kind = 8 )  y(maxnode)
  real         ( kind = 8 )  yval
  real         ( kind = 8 )  z(maxnode)
  real         ( kind = 8 )  zval

  ierror = 0

  call get_unit ( iunit )

  open ( unit = iunit, file = file_name, status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'VLA_TO_GRAPH_ARC - Fatal error!'
    write ( *, '(a)' ) '  Could not open the file.'
    return
  end if

  ierror = 0
  icor3 = 0
  jcor3 = 0
  num_bad = 0
  num_bad_old = 0
  nedge = 0
  nnode = 0
!
!  Read the next line.
!
  do

    read ( iunit, '(a)', iostat = ios ) text

    if ( ios /= 0 ) then
      exit
    end if

    done = .true.
    iword = 0
!
!  Read the next word.
!
    do

      call word_next_read ( text, word, done )
!
!  If no more words in this line, read a new line.
!
      if ( done ) then
        exit
      end if

      iword = iword + 1
!
!  The first word in the line tells us what's happening.
!
      if ( iword == 1 ) then
        word1 = word
      end if
!
!  If WORD1 is "SET", then we regard this line as comments.
!
      if ( s_eqi ( word1, 'set' ) ) then
!
!  If WORD1 is ";", then we regard this line as comments.
!
      else if ( word1 == ';' ) then
!
!  If WORD1 is "P", then this is the initial point on a line.
!  If WORD1 is "L", then this is a followup point on a line.
!
      else if ( s_eqi ( word1, 'P' ) .or. s_eqi ( word1, 'L' ) ) then
!
!  Read in the point coordinates.
!
        num_bad_old = num_bad

        do i = 1, 3

          call word_next_read ( text, word, done )

          if ( done ) then
            num_bad = num_bad + 1
            exit
          end if

          call s_to_r8 ( word, rval, ierror, lchar )

          if ( ierror /= 0 ) then
            num_bad = num_bad + 1
            exit
          end if

          if ( nnode <= maxnode ) then
            if ( i == 1 ) then
              xval = rval
            else if ( i == 2 ) then
              yval = rval
            else if ( i == 3 ) then
              zval = rval
            end if
          end if

        end do

        if ( num_bad_old < num_bad ) then
          exit
        end if
!
!  Assign a node index to the point.
!
        icor3_old = icor3
        jcor3_old = jcor3
!
!  ICOR3 is the index of the new value.
!  (If such a point already exists, a new one won't be added.)
!
        call r8vec3_index_insert_unique ( maxnode, nnode, x, y, z, indx, &
          xval, yval, zval, icor3, ierror )

        jcor3 = indx(icor3)

        if ( ierror /= 0 ) then
          ierror = 1
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'VLA_TO_GRAPH_ARC - Fatal error!'
          write ( *, '(a)' ) '  R8VEC3_INDEX_INSERT_UNIQUE returned an error!'
          return
        end if
!
!  Define the line as joining JCOR3_OLD to JCOR3.
!  (If such a line already exists, a new copy won't be added.)
!
        if ( s_eqi ( word1, 'L' ) ) then

          call iset2_index_insert_unique ( maxedge, nedge, inode, jnode, &
            indx_edge, jcor3_old, jcor3, iedge, ierror )

        end if

        exit
!
!  If the first word is unrecognized, then skip the whole line.
!
      else

        num_bad = num_bad + 1
        exit

      end if

    end do

  end do

  close ( unit = iunit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'VLA_TO_GRAPH_ARC - Note:'
  write ( *, '(a)' ) '  The graph was read properly.'
  write ( *, '(a,i8)' ) '  Number of nodes = ', nnode
  write ( *, '(a,i8)' ) '  Number of edges = ', nedge

  return
end
subroutine word_next_read ( line, word, done )

!*****************************************************************************80
!
!! WORD_NEXT_READ "reads" words from a string, one at a time.
!
!  Discussion:
!
!    The following characters are considered to be a single word,
!    whether surrounded by spaces or not:
!
!      " ( ) { } [ ]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) LINE, a string, presumably containing words
!    separated by spaces.
!
!    Output, character ( len = * ) WORD.
!
!    If DONE is FALSE, then WORD contains the "next" word read from LINE.
!    If DONE is TRUE, then WORD is blank, because there was no more to read.
!
!    Input/output, logical DONE.
!
!    On input with a fresh value of LINE, set DONE to TRUE.
!
!    On output, the routine sets DONE:
!      FALSE if another word was read from LINE,
!      TRUE if no more words could be read (LINE is exhausted).
!
  implicit none

  logical done
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ), save :: lenc = 0
  character ( len = * ) line
  integer ( kind = 4 ), save :: next = 1
  character TAB
  character ( len = * ) word

  TAB = char ( 9 )
!
!  An input value of DONE = TRUE signals a new line of text to examine.
!
  if ( done ) then

    next = 1
    done = .false.
    lenc = len_trim ( line )

    if ( lenc <= 0 ) then
      done = .true.
      word = ' '
      return
    end if

  end if
!
!  Beginning at index NEXT, search LINE for the next nonblank,
!  which signals the beginning of a word.
!
  ilo = next

10    continue
!
!  ...LINE(NEXT:) is blank.  Return with WORD = ' ' and DONE = TRUE.
!
  if ( lenc < ilo ) then
    word = ' '
    done = .true.
    next = lenc + 1
    return
  end if
!
!  If the current character is blank, skip to the next one.
!
  if ( line(ilo:ilo) == ' ' .or. line(ilo:ilo) == TAB ) then
    ilo = ilo + 1
    go to 10
  end if
!
!  ILO is the index of the next nonblank character in the string.
!
!  If this initial nonblank is a special character,
!  then that's the whole word as far as we're concerned,
!  so return immediately.
!
  if ( line(ilo:ilo) == '"' .or. &
       line(ilo:ilo) == '(' .or. &
       line(ilo:ilo) == ')' .or. &
       line(ilo:ilo) == '{' .or. &
       line(ilo:ilo) == '}' .or. &
       line(ilo:ilo) == '[' .or. &
       line(ilo:ilo) == ']' ) then

    word = line(ilo:ilo)
    next = ilo + 1
    return

  end if
!
!  Now search for the last contiguous character that is not a
!  blank, TAB, or special character.
!
  next = ilo + 1

20    continue

  if ( lenc < next ) then
    word = line(ilo:next-1)
    return
  end if

  if ( line(next:next) /= ' ' .and. &
       line(next:next) /= TAB .and. &
       line(next:next) /= '"' .and. &
       line(next:next) /= '(' .and. &
       line(next:next) /= ')' .and. &
       line(next:next) /= '{' .and. &
       line(next:next) /= '}' .and. &
       line(next:next) /= '[' .and. &
       line(next:next) /= ']' ) then

    next = next + 1
    go to 20

  end if

  word = line(ilo:next-1)

  return
end
