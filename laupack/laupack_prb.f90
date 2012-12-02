program main

!*****************************************************************************80
!
!! MAIN is the main program for LAUPACK_PRB.
!
!  Discussion:
!
!    LAUPACK_PRB calls the LAUPACK test routines.
!
!  Modified:
!
!    11 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LAUPACK_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the LAUPACK library.'

  call test01 ( )
  call test22 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
  call test05 ( )
  call test06 ( )
  call test07 ( )
  call test08 ( )
  call test09 ( )
  call test10 ( )

  call test11 ( )
  call test12 ( )
  call test13 ( )
  call test14 ( )
  call test15 ( )
  call test16 ( )
  call test17 ( )
  call test18 ( )
  call test19 ( )
  call test20 ( )

  call test21 ( )
  call test22 ( )
  call test23 ( )
  call test24 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LAUPACK_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 calls DIGRAPH_ARC_EULER.
!
  implicit none

  integer ( kind = 4 ), parameter :: nedge = 7
  integer ( kind = 4 ), parameter :: nnode = 5

  integer ( kind = 4 ) i
  integer ( kind = 4 ) in
  integer ( kind = 4 ), dimension ( nedge ) :: inode = (/ 2, 1, 2, 1, 3, 5, 4 /)
  integer ( kind = 4 ) j
  integer ( kind = 4 ), dimension ( nedge ) :: jnode = (/ 5, 4, 3, 2, 1, 1, 2 /)
  integer ( kind = 4 ) jp1
  logical success
  integer ( kind = 4 ) trail(nedge)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  DIGRAPH_ARC_EULER finds an Euler circuit of a digraph.'

  call digraph_arc_print ( nedge, inode, jnode, &
    '  The arc list of the digraph:' )

  call digraph_arc_euler ( nnode, nedge, inode, jnode, success, trail )

  if ( success ) then

    call i4vec_print ( nedge, trail, '  The edge list of the Euler circuit:' )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The node list of the Euler circuit:'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '    I  Edge  Node'
    write ( *, '(a)' ) ' '

    do i = 1, nedge

      j = trail(i)

      if ( i == nedge ) then
        jp1 = trail(1)
      else
        jp1 = trail(i+1)
      end if

      if ( jnode(j) == inode(jp1) ) then
        in = jnode(j)
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'The circuit has failed!'
        exit
      end if

      write ( *, '(i8,i8,i8)' ) i, j, in

    end do

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The digraph is not eulerian.'
    write ( *, '(a)' ) ' '

  end if

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 calls DIGRAPH_ARC_GET_PATH, DIGRAPH_ARC_KSHORT2.
!
  implicit none

  integer ( kind = 4 ), parameter :: kpaths = 4
  integer ( kind = 4 ), parameter :: nedge = 11
  integer ( kind = 4 ), parameter :: nnode = 6

  integer ( kind = 4 ), parameter :: maxque = 2 * kpaths

  integer ( kind = 4 ) arcbwd(nnode)
  integer ( kind = 4 ) arcdir(nnode)
  integer ( kind = 4 ) arcfwd(nnode)
  integer ( kind = 4 ), dimension ( nedge ) :: arclen = (/ &
    1, 3, 2, 3, 3, 6,-2, 2, 2, 4, 2 /)
  integer ( kind = 4 ) arcnod(nnode)
  integer ( kind = 4 ) auxdis(nnode)
  integer ( kind = 4 ) auxlnk(nnode)
  integer ( kind = 4 ) auxstg(nnode)
  integer ( kind = 4 ) auxtre(nnode)
  integer ( kind = 4 ) crosar(maxque)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ), dimension ( nedge ) :: inode = &
    (/ 4, 3, 6, 4, 2, 5, 6, 5, 1, 2, 4 /)
  integer ( kind = 4 ) ipaths
  integer ( kind = 4 ) isink
  integer ( kind = 4 ) isorce
  integer ( kind = 4 ), dimension ( nedge ) :: jnode = &
    (/ 2, 1, 5, 3, 1, 4, 1, 2, 3, 6, 1 /)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) length
  integer ( kind = 4 ) nextrd(maxque)
  integer ( kind = 4 ) nump
  integer ( kind = 4 ) nxtbwd(nedge)
  integer ( kind = 4 ) nxtfwd(nedge)
  integer ( kind = 4 ) pathlen(nnode+1)
  integer ( kind = 4 ) qufirp(kpaths+3)
  integer ( kind = 4 ) qunxtp(kpaths+3)
  integer ( kind = 4 ) trdist(nnode)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  DIGRAPH_ARC_KSHORT2 finds the K shortest paths without'
  write ( *, '(a)' ) '  repetition;'
  write ( *, '(a)' ) '  DIGRAPH_ARC_GET_PATH retrieves the computed paths.'

  call digraph_arc_print ( nedge, inode, jnode, &
    '  The arc list of the digraph:' )

  isorce = 5
  isink = 3
!
!  Find the shortest path lengths.
!
  call digraph_arc_kshort2 ( nnode, nedge, inode, jnode, arclen, kpaths, &
    maxque, isorce, isink, iflag, ipaths, pathlen, arcdir, trdist, arcnod, &
    arcfwd, arcbwd, auxstg, auxdis, auxtre, auxlnk, nxtfwd, nxtbwd, qufirp, &
    qunxtp, crosar, nextrd )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  IFLAG = ', iflag
  write ( *, '(a,i8,a)' ) '  The first ', ipaths, ' shortest paths '
  write ( *, '(a,i8,a,i8)' ) '  from ', isorce, ' to ', isink
  write ( *, '(a)' ) '  These path lengths are:'
  write ( *, '(a)' ) ' '

  do i = 2, ipaths + 1
    write ( *, '(i8)' ) pathlen(i)
  end do
!
!  Output the paths
!
  do k = 1, ipaths

    call digraph_arc_get_path ( nnode, nedge, k, kpaths, maxque, inode, jnode, &
      pathlen, arcdir, qufirp, qunxtp, crosar, nextrd, nump, length, arcnod )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8,a,i8,a)' ) '  Path ', k, ' has ', nump, ' edges,'
    write ( *, '(a,i8)' ) '  and a length of ', length
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The edges in the path:'
    write ( *, '(a)' ) ' '
    write ( *, '(10i8)' ) arcnod(1:nump)

  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests DIGRAPH_ARC_HAMCYC.
!
  implicit none

  integer ( kind = 4 ), parameter :: nedge = 16
  integer ( kind = 4 ), parameter :: nnode = 10

  integer ( kind = 4 ) hcycle(nnode)
  integer ( kind = 4 ), dimension ( nedge ) :: inode = &
    (/ 6, 5,  8, 5, 4, 1, 10, 7, 9, 2, 7, 3, 9, 3, 10, 6 /)
  integer ( kind = 4 ), dimension ( nedge ) :: jnode = &
    (/ 2, 3, 10, 2, 9, 8,  5, 4, 5, 4, 1, 8, 7, 7,  6, 1 /)
  logical success

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  DIGRAPH_ARC_HAMCYC finds a hamiltonian circuit.'

  call digraph_arc_print ( nedge, inode, jnode, &
    '  The arc list of the digraph:' )

  call digraph_arc_hamcyc ( nnode, nedge, inode, jnode, success, hcycle )

  if ( success ) then

    call i4vec_print ( nnode, hcycle, '  The Hamiltonian circuit:' )

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The digraph is not hamiltonian.'

  end if

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 calls DIGRAPH_ARC_KSHORT1, DIGRAPH_ARC_PRT_PATH.
!
  implicit none

  integer ( kind = 4 ), parameter :: k = 4
  integer ( kind = 4 ), parameter :: nedge = 14
  integer ( kind = 4 ), parameter :: nnode = 6
  integer ( kind = 4 ), parameter :: npmax = 20

  real ( kind = 8 ), dimension ( nedge ) :: arclen = (/ &
     2.0D+00, &
    -2.0D+00, &
     3.0D+00, &
    -3.0D+00, &
     4.0D+00, &
     6.0D+00, &
     4.0D+00, &
     4.0D+00, &
    -3.0D+00, &
     9.0D+00, &
     6.0D+00, &
     4.0D+00, &
     4.0D+00, &
     2.0D+00 /)
  real ( kind = 8 ) dist(nnode,k)
  integer ( kind = 4 ) i
  integer ( kind = 4 ), dimension ( nedge ) :: inode = &
    (/ 6, 3, 4, 5, 6, 2, 2, 2, 1, 4, 1, 3, 2, 4 /)
  integer ( kind = 4 ) isink
  integer ( kind = 4 ) :: isorce = 2
  integer ( kind = 4 ) iter
  integer ( kind = 4 ), dimension ( nedge ) :: jnode = &
    (/ 1, 1, 2, 3, 3, 3, 3, 4, 4, 5, 5, 5, 6, 6 /)
  integer ( kind = 4 ) maxpath

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  DIGRAPH_ARC_KSHORT1 computes the K shortest path'
  write ( *, '(a)' ) '  lengths from a given node to all others,'
  write ( *, '(a)' ) '  DIGRAPH_ARC_PRT_PATH prints them out.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Input K = ', k
  write ( *, '(a,i8)' ) '  Starting node ', isorce

  call digraph_arc_weight_print ( nedge, inode, jnode, arclen, &
    '  The arc list of the weighted digraph:' )

  iter = 0

  call digraph_arc_kshort1 ( nnode, nedge, k, isorce, iter, inode, jnode, &
    arclen, dist )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Total number of iterations = ', iter
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Distance array.'
  write ( *, '(a)') ' '

  do i = 1, nnode
    write ( *, '(4x,i8,4g14.6)' ) i, dist(i,1:k)
  end do
!
!  Generate at most 10 of the shortest paths from ISORCE to ISINK.
!
  isink = 5
  maxpath = 10

  call digraph_arc_prt_path ( nnode, nedge, npmax, k, maxpath, isorce, &
    isink, inode, jnode, arclen, dist )

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 tests DIGRAPH_ARC_MINEQV.
!
  implicit none

  integer ( kind = 4 ), parameter :: nedge = 10
  integer ( kind = 4 ), parameter :: nnode = 5

  logical arclis(nedge)
  integer ( kind = 4 ) i
  integer ( kind = 4 ), dimension ( nedge ) :: inode = (/ 5,2,2,4,1,3,1,5,1,3 /)
  integer ( kind = 4 ), dimension ( nedge ) :: inode2
  integer ( kind = 4 ), dimension ( nedge ) :: jnode = (/ 2,3,4,5,5,4,2,3,4,1 /)
  integer ( kind = 4 ), dimension ( nedge ) :: jnode2
  integer ( kind = 4 ) nedge2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  DIGRAPH_ARC_MINEQV finds the minimal equivalent'
  write ( *, '(a)' ) '  digraph.'

  call digraph_arc_print ( nedge, inode, jnode, &
    '  The arc list of the digraph:' )

  call digraph_arc_mineqv ( nnode, nedge, inode, jnode, arclis )

  nedge2 = 0

  do i = 1, nedge

    if ( arclis(i) ) then
      nedge2 = nedge2 + 1
      inode2(nedge2) = inode(i)
      jnode2(nedge2) = jnode(i)
    end if

  end do

  call digraph_arc_print ( nedge2, inode2, jnode2, &
    '  The arc list of the minimal equivalent digraph:' )

  nedge2 = 6
  inode2(1) = 5
  jnode2(1) = 2
  inode2(2) = 2
  jnode2(2) = 3
  inode2(3) = 4
  jnode2(3) = 5
  inode2(4) = 1
  jnode2(4) = 5
  inode2(5) = 3
  jnode2(5) = 4
  inode2(6) = 3
  jnode2(6) = 1

  call digraph_arc_print ( nedge2, inode2, jnode2, &
    '  The arc list of the correct answer:' )

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 tests DIGRAPH_ARC_NFLOW
!
  implicit none

  integer ( kind = 4 ), parameter :: nnode = 6
  integer ( kind = 4 ), parameter :: nedge = 20

  integer ( kind = 4 ), dimension ( nedge ) :: capac = &
    (/ 3, 7, 2, 5, 4, 1, 4, 2, 8, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
  integer ( kind = 4 ), dimension ( nedge ) :: flow
  integer ( kind = 4 ), dimension ( nedge ) :: inode = &
    (/ 1, 1, 2, 2, 2, 3, 3, 4, 4, 5, 2, 3, 3, 4, 5, 4, 5, 5, 6, 6 /)
  integer ( kind = 4 ) :: isink = 6
  integer ( kind = 4 ) :: isorce = 1
  integer ( kind = 4 ), dimension ( nedge ) :: jnode = &
    (/ 2, 3, 3, 4, 5, 4, 5, 5, 6, 6, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5 /)
  integer ( kind = 4 ) mincut(nnode)
  integer ( kind = 4 ), dimension ( nnode ) :: nodflow

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  DIGRAPH_ARC_NFLOW finds the maximum flow '
  write ( *, '(a)' ) '  on a network'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The source is node ', isorce
  write ( *, '(a,i8)' ) '  The sink is node   ', isink

  call digraph_arc_print ( nedge, inode, jnode, &
    '  The edge list of the digraph:' )

  call i4vec_print ( nedge, capac, '  The edge capacities:' )

  call digraph_arc_nflow ( nnode, nedge, inode, jnode, capac, isorce, &
    isink, mincut, flow, nodflow )

  call i4vec_print ( nnode, mincut, '  0/1 node cutset array:' )

  call i4vec_print ( nedge, flow, '  Edge flows:' )

  call i4vec_print ( nnode, nodflow, '  Node flows:' )

  return
end
subroutine test07 ( )

!*****************************************************************************80
!
!! TEST07 tests DIGRAPH_ARC_SHTREE.
!
  implicit none

  integer ( kind = 4 ), parameter :: nedge = 11
  integer ( kind = 4 ), parameter :: nnode = 6

  real ( kind = 8 ), dimension ( nedge ) :: arclen = (/ &
     3.0D+00, &
     4.0D+00, &
     2.0D+00, &
     4.0D+00, &
    -2.0D+00, &
    -1.0D+00, &
     5.0D+00, &
     2.0D+00, &
     7.0D+00, &
     6.0D+00, &
     9.0D+00 /)
  real ( kind = 8 ) dist(nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ), dimension ( nedge ) :: inode = &
    (/  3,  2,  5,  4,  1,  1,  6,  6,  1,  3,  6 /)
  integer ( kind = 4 ) :: iroot = 3
  integer ( kind = 4 ), dimension ( nnode - 1 ) :: itree
  integer ( kind = 4 ) j
  integer ( kind = 4 ), dimension ( nedge ) :: jnode = &
    (/  6,  5,  4,  2,  4,  5,  4,  1,  3,  5,  2 /)
  integer ( kind = 4 ), dimension ( nnode - 1 ) :: jtree
  real ( kind = 8 ), dimension ( nnode - 1 ) :: wtree

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a)' ) '  DIGRAPH_ARC_SHTREE constructs a tree of the shortest'
  write ( *, '(a)' ) '  paths from one node to all others.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The base node will be ', iroot

  call digraph_arc_weight_print ( nedge, inode, jnode, arclen, &
    '  The weighted arc list of the digraph:' )

  call digraph_arc_shtree ( nnode, nedge, iroot, inode, jnode, arclen, dist, &
    itree, jtree )

  call r8vec_print ( nnode, dist, '  Distance from base node to other nodes:' )

  do i = 1, nnode-1
    wtree(i) = 0.0D+00
    do j = 1, nedge
      if ( inode(j) == itree(i) .and. jnode(j) == jtree(i) ) then
        wtree(i) = arclen(j)
      end if
    end do
  end do

  call digraph_arc_weight_print ( nnode-1, itree, jtree, wtree, &
    '  The weighted arc list of the shortest path tree:' )

  return
end
subroutine test08 ( )

!*****************************************************************************80
!
!! TEST08 tests DIGRAPH_ARC_STCOMP.
!
!  Find the strongly connected components of the following digraph:
!
!
!  1--><--9-->--4
!
!  5--<--7-->--2--<--8--><--10-->--3-->---6
!    \   |    /                     \    /
!     V  A  V                        A  V
!       \|/                           \/
!        12                           11
!
  implicit none

  integer ( kind = 4 ), parameter :: nedge = 15
  integer ( kind = 4 ), parameter :: nnode = 12

  integer ( kind = 4 ), dimension ( nnode ) :: comp
  integer ( kind = 4 ), dimension ( nedge ) :: inode = &
    (/ 8,11, 9, 5, 3, 8, 9, 6,12,10, 7, 1, 2,10, 7 /)
  integer ( kind = 4 ), dimension ( nedge ) :: jnode = &
    (/ 2, 3, 1,12, 6,10, 4,11, 7, 3, 2, 9,12, 8, 5 /)
  integer ( kind = 4 ) numcomp

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08'
  write ( *, '(a)' ) '  DIGRAPH_ARC_STCOMP finds the strongly connected'
  write ( *, '(a)' ) '  components of a directed graph.'

  call digraph_arc_print ( nedge, inode, jnode, &
    '  The arc list of the digraph:' )

  call digraph_arc_stcomp ( nnode, nedge, inode, jnode, numcomp, comp )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of strong components = ', numcomp

  call i4vec_print ( nnode, comp, '  The component of each node:' )

  return
end
subroutine test09 ( )

!*****************************************************************************80
!
!! TEST09 calls DIGRAPH_DIST_ALLPATH.
!
  implicit none

  integer ( kind = 4 ), parameter :: nnode = 6
  integer ( kind = 4 ), parameter :: lda = nnode

  real ( kind = 8 ) dist(lda,nnode)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) next(nnode,nnode)
  integer ( kind = 4 ) node
  integer ( kind = 4 ) node1
  integer ( kind = 4 ) node2
  integer ( kind = 4 ) num
  integer ( kind = 4 ) path(nnode)

  data dist / &
     0.0D+00,  3.0D+00,  2.0D+00, 99.0D+00, 99.0D+00,  1.0D+00, &
     3.0D+00,  0.0D+00,  8.0D+00, 99.0D+00,  2.0D+00,  2.0D+00, &
    99.0D+00,  8.0D+00,  0.0D+00,  2.0D+00, 99.0D+00, 99.0D+00, &
    99.0D+00, 99.0D+00,  2.0D+00,  0.0D+00,  8.0D+00,  1.0D+00, &
    99.0D+00,  2.0D+00, 99.0D+00,  8.0D+00,  0.0D+00,  9.0D+00, &
     1.0D+00, 99.0D+00, 99.0D+00,  1.0D+00,  9.0D+00,  0.0D+00 /

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST09'
  write ( *, '(a)' ) '  DIGRAPH_DIST_ALLPATH computes the shortest distance '
  write ( *, '(a)' ) '  between all pairs of nodes.'

  call digraph_dist_print ( dist, lda, nnode, '  The distance matrix:' )

  call digraph_dist_allpath ( nnode, dist, lda, next )
!
!  Store the shortest path from node 3 to node 5 in PATH.
!
  node1 = 3
  node2 = 5
  j = node1
  num = 1

  path(num) = node1

  do

    node = next(node2,j)
    num = num + 1
    path(num) = node

    if ( node == node2 ) then
      exit
    end if

    j = node

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The shortest path from '
  write ( *, '(i8,a,i8)' ) node1, ' to ', node2
  write ( *, '(a)' ) ' '
  write ( *, '(4x,20i3)' ) path(1:num)

  return
end
subroutine test10 ( )

!*****************************************************************************80
!
!! TEST10 tests DIGRAPH_DIST_SHORTP.
!
  implicit none

  integer ( kind = 4 ), parameter :: nnode = 6
  integer ( kind = 4 ), parameter :: lda = nnode

  real ( kind = 8 ) dist(lda,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) :: isink = 2
  integer ( kind = 4 ) :: isorce = 3
  integer ( kind = 4 ) j
  real ( kind = 8 ) lnpath
  integer ( kind = 4 ) numnod
  integer ( kind = 4 ) path(nnode)

  data ( ( dist(i,j), i = 1, nnode ), j = 1, nnode ) / &
       0.0D+00, 99.0D+00, 99.0D+00, 99.0D+00, 99.0D+00,  2.0D+00, &
      99.0D+00,  0.0D+00, 99.0D+00,  4.0D+00, 99.0D+00,  9.0D+00, &
       7.0D+00, 99.0D+00,  0.0D+00, 99.0D+00, 99.0D+00, 99.0D+00, &
       2.0D+00, 99.0D+00, 99.0D+00,  0.0D+00,  2.0D+00,  5.0D+00, &
       1.0D+00,  4.0D+00,  6.0D+00, 99.0D+00,  0.0D+00, 99.0D+00, &
      99.0D+00, 99.0D+00,  3.0D+00, 99.0D+00, 99.0D+00,  0.0D+00 /

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST10'
  write ( *, '(a)' ) '  DIGRAPH_DIST_SHORTP finds the shortest path between '
  write ( *, '(a)' ) '  two nodes.'
  write ( *, '(a,i8)' ) '  Start node is  ', isorce
  write ( *, '(a,i8)' ) '  Finish node is ', isink

  call digraph_dist_print ( dist, lda, nnode, '  The distance matrix:' )

  call digraph_dist_shortp ( nnode, dist, lda, isorce, isink, numnod, path, &
    lnpath )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  The path length is ', lnpath

  call i4vec_print ( numnod, path, '  The shortest path:' )

  return
end
subroutine test11 ( )

!*****************************************************************************80
!
!! TEST11 calls DIGRAPH_DIST_SHORT_LN.
!
  implicit none

  integer ( kind = 4 ), parameter :: nnode = 6
  integer ( kind = 4 ), parameter :: lda = nnode

  real ( kind = 8 ) dist(lda,nnode)
  integer ( kind = 4 ) :: root = 5
  real ( kind = 8 ) shdist(nnode)

  data dist / &
     0.0D+00,  3.0D+00,  2.0D+00, 99.0D+00, 99.0D+00,  1.0D+00, &
     3.0D+00,  0.0D+00,  8.0D+00, 99.0D+00,  2.0D+00,  2.0D+00, &
    99.0D+00,  8.0D+00,  0.0D+00,  2.0D+00, 99.0D+00, 99.0D+00, &
    99.0D+00, 99.0D+00,  2.0D+00,  0.0D+00,  8.0D+00,  1.0D+00, &
    99.0D+00,  2.0D+00, 99.0D+00,  8.0D+00,  0.0D+00,  9.0D+00, &
     1.0D+00, 99.0D+00, 99.0D+00,  1.0D+00,  9.0D+00,  0.0D+00 /

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST11'
  write ( *, '(a)' ) '  DIGRAPH_DIST_SHORT_LN finds shortest paths from one'
  write ( *, '(a)' ) '  node to all others.'
  write ( *, '(a,i8)' ) '  The root node will be ', root

  call digraph_dist_print ( dist, lda, nnode, '  The distance matrix:' )

  call digraph_dist_short_ln ( nnode, dist, lda, root, shdist )

  call r8vec_print ( nnode, shdist, '  Root-to-Node distances:' )

  return
end
subroutine test12 ( )

!*****************************************************************************80
!
!! TEST12 tests GRAPH_ARC_CLIQUE.
!
  implicit none

  integer ( kind = 4 ), parameter :: nedge = 11
  integer ( kind = 4 ), parameter :: nnode = 9

  integer ( kind = 4 ), dimension ( nedge ) :: inode = (/ &
    5, 1, 3, 6, 2, 5, 5, 1, 2, 7, 8 /)
  integer ( kind = 4 ), dimension ( nedge ) :: jnode = (/ &
    8, 9, 9, 9, 7, 2, 7, 3, 8, 8, 4 /)
  integer ( kind = 4 ) cliq(nnode)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST12'
  write ( *, '(a)' ) '  GRAPH_ARC_CLIQUE finds cliques in a graph.'

  call graph_arc_print ( nedge, inode, jnode, '  The arc list of the graph:' )

  call graph_arc_clique ( nnode, nedge, inode, jnode, cliq )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The correct answer is :'
  write ( *, '(a)' ) '    (8,2,5,7), (8,4), (3,1,9), (6,9).'

  return
end
subroutine test13 ( )

!*****************************************************************************80
!
!! TEST13 tests GRAPH_ARC_COLOR_NUMBER.
!
!  Here is the graph to be colored.
!
!    6----8----3    7---11
!    |\  / \  /|    |    |
!    | \/   \/ |    |    |
!    | 1----2  |   10----5
!    | |    |  |
!    \ |    | /
!     \|    |/
!      4----9
!
!
  implicit none

  integer ( kind = 4 ), parameter :: nedge = 16
  integer ( kind = 4 ), parameter :: nnode = 11

  integer ( kind = 4 ) color(nnode)
  integer ( kind = 4 ), dimension ( nedge ) :: inode = &
    (/ 3, 6, 7,1,4,10,2,6,8, 5,8,10,1,3,2,9 /)
  integer ( kind = 4 ), dimension ( nedge ) :: jnode = &
    (/ 9, 8,11,2,1, 7,8,4,3,11,1, 5,6,2,9,4 /)
  integer ( kind = 4 ) ncolor

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST13'
  write ( *, '(a)' ) '  GRAPH_ARC_COLOR_NUMBER finds the chromatic number of a'
  write ( *, '(a)' ) '  graph and exhibits a corresponding node coloring.'

  call graph_arc_print ( nedge, inode, jnode, '  The arc list of the graph:' )

  call graph_arc_color_number ( nnode, nedge, inode, jnode, ncolor, color )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Computed answers:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The chromatic number of the graph, that is,'
  write ( *, '(a,i8)' ) '  the number of colors required for the nodes = ', &
    ncolor

  call graph_arc_ncolor_print ( nedge, inode, jnode, nnode, color, &
    '  The arc list and node colors of the graph:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Correct answers:'

  ncolor = 4
  color = (/ 1, 2, 1, 2, 1, 3, 1, 4, 3, 2, 2 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The chromatic number of the graph, that is,'
  write ( *, '(a,i8)' ) '  the number of colors required for the nodes = ', &
    ncolor

  call graph_arc_ncolor_print ( nedge, inode, jnode, nnode, color, &
    '  The arc list and node colors of the graph:' )

  return
end
subroutine test14 ( )

!*****************************************************************************80
!
!! TEST14 tests GRAPH_ARC_COLOR_POLY.
!
  implicit none

  integer ( kind = 4 ), parameter :: nnode = 5
  integer ( kind = 4 ), parameter :: nedge = 7

  integer ( kind = 4 ) cpoly1(nnode)
  integer ( kind = 4 ) cpoly2(nnode)
  integer ( kind = 4 ) cpoly3(nnode)
  integer ( kind = 4 ), dimension ( nedge ) :: inode = &
    (/ 1, 1, 5, 1, 3, 1, 5 /)
  integer ( kind = 4 ), dimension ( nedge ) :: jnode = &
    (/ 2, 3, 2, 4, 2, 5, 4 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST14'
  write ( *, '(a)') '  GRAPH_ARC_COLOR_POLY computes the chromatic polynomial.'

  call graph_arc_color_poly ( nnode, nedge, inode, jnode, cpoly1, cpoly2, &
    cpoly3 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Computed values:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,5i4)' ) '    CPOLY1: ', cpoly1(1:nnode)
  write ( *, '(a,5i4)' ) '    CPOLY2: ', cpoly2(1:nnode)
  write ( *, '(a,5i4)' ) '    CPOLY3: ', cpoly3(1:nnode)

  cpoly1(1:nnode) = (/ 8, 20, 18, 7, 1 /)
  cpoly2(1:nnode) = (/ 0,  1,  3, 3, 1 /)
  cpoly3(1:nnode) = (/ 0,  0,  1, 3, 1 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Correct values:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,5i4)' ) '    CPOLY1: ', cpoly1(1:nnode)
  write ( *, '(a,5i4)' ) '    CPOLY2: ', cpoly2(1:nnode)
  write ( *, '(a,5i4)' ) '    CPOLY3: ', cpoly3(1:nnode)

  return
end
subroutine test15 ( )

!*****************************************************************************80
!
!! TEST15 tests GRAPH_ARC_CONECT.
!
  implicit none

  integer ( kind = 4 ), parameter :: nedge = 17
  integer ( kind = 4 ), parameter :: nnode = 14

  integer ( kind = 4 ) bridge(nedge)
  integer ( kind = 4 ) cutnod(nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ), dimension ( nedge ) :: inode = &
    (/ 3, 4,5,13, 7,9,14,1,3,9, 7,9,5,9,1,6,4 /)
  integer ( kind = 4 ) iroot
  integer ( kind = 4 ), dimension ( nedge ) :: jnode = &
    (/ 6,10,9,11,12,8,7,11,8,1,10,3,2,6,13,8,7 /)
  integer ( kind = 4 ) nbridg
  integer ( kind = 4 ) ncut
  integer ( kind = 4 ) next(nnode)
  integer ( kind = 4 ) numcmp

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST15'
  write ( *, '(a)' ) '  GRAPH_ARC_CONECT finds bridges, blocks, and cut nodes.'

  call graph_arc_print ( nedge, inode, jnode, '  The arc list of the graph:' )

  next(1:nnode) = 0

  numcmp = 0

  do iroot = 1, nnode

    if ( next(iroot) == 0 ) then

      numcmp = numcmp + 1

      call graph_arc_conect ( nnode, nedge, inode, jnode, iroot, ncut, &
        nbridg, cutnod, bridge, next )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Using component with root node:', iroot
      write ( *, '(a,i8)' ) '  The number of cut nodes is ', ncut
      write ( *, '(a,i8)' ) '  The number of bridges is ', nbridg

       if ( ncut > 0 ) then
         write ( *, '(a)' ) ' '
         write ( *, '(a)' ) '  The cut nodes are:'
         write ( *, '(a)' ) ' '
         write ( *, '(10i8)' ) cutnod(1:ncut)
       end if

       if ( nbridg > 0 ) then

         write ( *, '(a)' ) ' '
         write ( *, '(a)' ) '  The bridges are:'
         write ( *, '(a)' ) ' '

         do i = 1, nedge

           if ( bridge(i) > 0 ) then
             write ( *, '(4x,3i3)' ) i, inode(i), jnode(i)
           end if

         end do

       end if

     end if

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of components is ', numcmp

  return
end
subroutine test16 ( )

!*****************************************************************************80
!
!! TEST16 tests GRAPH_ARC_EDGE_CON.
!
  implicit none

  integer ( kind = 4 ), parameter :: nedge = 17
  integer ( kind = 4 ), parameter :: nnode = 9

  integer ( kind = 4 ) edge_con
  integer ( kind = 4 ), dimension ( nedge ) :: inode = &
    (/ 6,2,3,6,7,1,4,7,3,4,9,6,5,4,2,9,4 /)
  integer ( kind = 4 ), dimension ( nedge ) :: jnode = &
    (/ 8,5,1,3,2,8,3,5,8,1,2,1,9,8,6,7,2 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST16'
  write ( *, '(a)' ) '  GRAPH_ARC_EDGE_CON finds graph edge connectivity.'

  call graph_arc_print ( nedge, inode, jnode, '  The arc list of the graph:' )

  call graph_arc_edge_con ( nnode, nedge, inode, jnode, edge_con )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The computed result:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The computed edge connectivity is ', edge_con

  edge_con = 2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The correct result:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The computed edge connectivity is ', edge_con

  return
end
subroutine test17 ( )

!*****************************************************************************80
!
!! TEST17 calls GRAPH_ARC_EULER.
!
  implicit none

  integer ( kind = 4 ), parameter :: nedge = 7
  integer ( kind = 4 ), parameter :: nnode = 5

  integer ( kind = 4 ) i
  integer ( kind = 4 ) in
  integer ( kind = 4 ), dimension ( nedge ) :: inode = (/ 2, 1, 2, 1, 3, 5, 4 /)
  integer ( kind = 4 ) j
  integer ( kind = 4 ), dimension ( nedge ) :: jnode = (/ 5, 4, 3, 2, 1, 1, 2 /)
  integer ( kind = 4 ) jp1
  logical success
  integer ( kind = 4 ) trail(nedge)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST17'
  write ( *, '(a)' ) '  GRAPH_ARC_EULER finds an Euler circuit of a graph.'

  call graph_arc_print ( nedge, inode, jnode, '  The arc list of the graph:' )

  call graph_arc_euler ( nnode, nedge, inode, jnode, success, trail )

  if ( success ) then

    call i4vec_print ( nedge, trail, '  The edge list of the Euler circuit:' )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The node list of the Euler circuit:'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '    I  Edge  Node'
    write ( *, '(a)' ) ' '

    do i = 1, nedge

      j = trail(i)

      if ( i == nedge ) then
        jp1 = trail(1)
      else
        jp1 = trail(i+1)
      end if

      if ( inode(j) == inode(jp1) .or. inode(j) == jnode(jp1) ) then
        in = inode(j)
      else if ( jnode(j) == inode(jp1) .or. jnode(j) == jnode(jp1) ) then
        in = jnode(j)
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'The circuit has failed!'
        exit
      end if

      write ( *, '(3i8)' ) i, j, in

    end do

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The graph is not eulerian.'
    write ( *, '(a)' ) ' '

  end if

  return
end
subroutine test18 ( )

!*****************************************************************************80
!
!! TEST18 tests GRAPH_ARC_FCYCLE.
!
  implicit none

  integer ( kind = 4 ), parameter :: nedge = 14
  integer ( kind = 4 ), parameter :: nnode = 13

  integer ( kind = 4 ), dimension ( nedge ) :: inode = &
    (/ 5, 4, 11, 6, 5,  7,  1, 8, 9, 10, 3,  1,  4, 6 /)
  integer ( kind = 4 ), dimension ( nedge ) :: jnode = &
    (/ 9, 7, 13, 8, 2, 12, 11, 3, 1,  7, 9, 13, 10, 9 /)
  integer ( kind = 4 ) ncyc(nnode)
  integer ( kind = 4 ) numcmp
  integer ( kind = 4 ) numcyc

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST18:'
  write ( *, '(a)' ) '  GRAPH_ARC_FCYCLE finds fundamental cycles of a graph.'

  call graph_arc_print ( nedge, inode, jnode, '  The arc list of the graph:' )

  call graph_arc_fcycle ( nnode, nedge, inode, jnode, numcyc, numcmp, ncyc )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Computed results:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of cycles =     ', numcyc
  write ( *, '(a,i8)' ) '  Number of components = ', numcmp

  numcyc = 3
  numcmp = 2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Correct results:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of cycles =     ', numcyc
  write ( *, '(a,i8)' ) '  Number of components = ', numcmp

  return
end
subroutine test19 ( )

!*****************************************************************************80
!
!! TEST19 tests GRAPH_ARC_HAMCYC.
!
  implicit none

  integer ( kind = 4 ), parameter :: nedge = 16
  integer ( kind = 4 ), parameter :: nnode = 10

  integer ( kind = 4 ) hcycle(nnode)
  integer ( kind = 4 ), dimension ( nedge ) :: inode = &
    (/ 2, 3,  8, 5, 4, 8, 10, 7, 5, 2, 1, 3, 9, 3, 6, 6 /)
  integer ( kind = 4 ), dimension ( nedge ) :: jnode = &
    (/ 6, 5, 10, 2, 9, 1,  5, 4, 9, 4, 7, 8, 7, 7, 10, 1 /)
  logical success

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST19'
  write ( *, '(a)' ) '  GRAPH_ARC_HAMCYC finds a hamiltonian circuit.'

  call graph_arc_print ( nedge, inode, jnode, '  The arc list of the graph:' )

  call graph_arc_hamcyc ( nnode, nedge, inode, jnode, success, hcycle )

  if ( success ) then

    call i4vec_print ( nnode, hcycle, '  The Hamiltonian circuit:' )

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'The graph is not hamiltonian.'

  end if

  return
end
subroutine test20 ( )

!*****************************************************************************80
!
!! TEST20 tests GRAPH_ARC_MAKEG.
!
  implicit none

  integer ( kind = 4 ), parameter :: edge_con = 5
  integer ( kind = 4 ), parameter :: nnode = 8
  integer ( kind = 4 ), parameter :: nedge = ( ( nnode * edge_con ) / 2 )

  integer ( kind = 4 ) edge_con2
  integer ( kind = 4 ) inode(nedge)
  integer ( kind = 4 ) jnode(nedge)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST20:'
  write ( *, '(a)' ) '  GRAPH_ARC_MAKEG makes a graph of given connectivity.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Requested:'
  write ( *, '(a,i8)' ) '    Number of nodes =      ', nnode
  write ( *, '(a,i8)' ) '    Number of edges =      ', nedge
  write ( *, '(a,i8)' ) '    Edge connectivity =    ', edge_con

  call graph_arc_makeg ( nnode, edge_con, inode, jnode )

  call graph_arc_print ( nedge, inode, jnode, '  The arc list of the graph:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Call GRAPH_ARC_EDGE_CON to verify the connectivity.'

  call graph_arc_edge_con ( nnode, nedge, inode, jnode, edge_con2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The computed edge connectivity is ', edge_con2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The correct answer:'
  write ( *, '(a)' ) ' '

  inode = (/ 1, 2, 3, 4, 5, 6, 7, 8, 1, 1, 2, 2, 3, 4, 5, 6, 1, 2, 3, 4 /)
  jnode = (/ 2, 3, 4, 5, 6, 7, 8, 1, 3, 7, 4, 8, 5, 6, 7, 8, 5, 6, 7, 8 /)

  call graph_arc_print ( nedge, inode, jnode, '  The arc list of the graph:' )

  return
end
subroutine test21 ( )

!*****************************************************************************80
!
!  TEST21 tests GRAPH_ARC_MATCH.
!
  implicit none

  integer ( kind = 4 ), parameter :: nedge = 14
  integer ( kind = 4 ), parameter :: nnode = 12

  integer ( kind = 4 ) i
  integer ( kind = 4 ), dimension ( nedge ) :: inode = &
    (/ 6, 9, 3,  4, 11, 6, 4, 5,  6, 10, 3, 4, 1, 3 /)
  integer ( kind = 4 ), dimension ( nedge ) :: jnode = &
    (/ 2, 7, 7, 10,  5, 8, 6, 7, 12,  2, 1, 2, 5, 5 /)
  integer ( kind = 4 ) notmat
  integer ( kind = 4 ) pair(nnode)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST21'
  write ( *, '(a)' ) '  GRAPH_ARC_MATCH finds a maximal matching in a graph.'

  call graph_arc_print ( nedge, inode, jnode, '  The edge list of the graph:' )

  call graph_arc_match ( nnode, nedge, inode, jnode, pair, notmat )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Node, matching node'
  write ( *, '(a)' ) ' '

  do i = 1, nnode
    write ( *, '(4x,i8,i8)' ) i, pair(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of unmatched nodes is ', notmat

  return
end
subroutine test22 ( )

!*****************************************************************************80
!
!! TEST22 tests GRAPH_ARC_MINTR2.
!
  implicit none

  integer ( kind = 4 ), parameter :: nedge = 9
  integer ( kind = 4 ), parameter :: nnode = 6

  real ( kind = 8 ), dimension ( nedge ) :: arclen = (/ &
    7.0D+00, &
    3.0D+00, &
    5.0D+00, &
    1.0D+00, &
    2.0D+00, &
    2.0D+00, &
    3.0D+00, &
    4.0D+00, &
    4.0D+00 /)
  integer ( kind = 4 ) i
  integer ( kind = 4 ), dimension ( nedge ) :: inode = &
    (/ 2, 4, 6, 2, 3, 4, 1, 2, 5 /)
  integer ( kind = 4 ) itree(nnode)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jtree(nnode)
  integer ( kind = 4 ), dimension ( nedge ) :: jnode = &
    (/ 5, 1, 3, 4, 5, 6, 3, 6, 1 /)
  integer ( kind = 4 ) ntree
  real ( kind = 8 ), dimension ( nnode-1 ) :: wtree

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST22'
  write ( *, '(a)' ) '  GRAPH_ARC_MINTR2 finds the minimum spanning tree:'

  call graph_arc_weight_print ( nedge, inode, jnode, arclen, &
    '  The weighted arc list of the graph:' )

  call graph_arc_mintr2 ( nnode, nedge, inode, jnode, arclen, ntree, &
    itree, jtree )

  do i = 1, ntree
    wtree(i) = 0.0
    do j = 1, nedge
      if ( ( inode(j) == itree(i) .and. jnode(j) == jtree(i) ) .or. &
           ( inode(j) == jtree(i) .and. jnode(j) == itree(i) ) ) then
        wtree(i) = arclen(j)
        exit
      end if
    end do
  end do

  call graph_arc_weight_print ( ntree, itree, jtree, wtree, &
    '  The weighted arc list of the tree:' )

  return
end
subroutine test23 ( )

!*****************************************************************************80
!
!! TEST23 tests GRAPH_ARC_PLANAR.
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10
  integer ( kind = 4 ), parameter :: m = 22

  integer ( kind = 4 ), dimension ( m ) :: inode = &
    (/ 9, 3, 10, 2, 6, 7,  7, 1, 8, 1, 10, 5, 1, 4, 1, 3, 9, 3, 2, 2, 5, 3 /)
  integer ( kind = 4 ), dimension ( m ) :: jnode = &
    (/ 8, 2,  5, 4, 5, 9, 10, 3, 6, 2,  2, 8, 9, 7, 4, 5, 4, 6, 7, 5, 9, 8/)
  logical planar

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST23'
  write ( *, '(a)' ) '  GRAPH_ARC_PLANAR determines if a graph is planar.'

  call graph_arc_print ( m, inode, jnode, '  The edge list of the graph.' )

  call graph_arc_planar ( n, m, inode, jnode, planar )

  if ( planar ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The input graph is planar.'
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The input graph is nonplanar.'
  end if

  return
end
subroutine test24 ( )

!*****************************************************************************80
!
!! TEST24 tests GRAPH_DIST_MINTR1.
!
  implicit none

  integer ( kind = 4 ), parameter :: nnode = 6
  integer ( kind = 4 ), parameter :: lda = nnode

  real ( kind = 8 ) dist(lda,nnode)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) itree(nnode-1)
  integer ( kind = 4 ) jtree(nnode-1)
  real ( kind = 8 ) wtree(nnode-1)

  data dist / &
     0.0D+00, 99.0D+00,  3.0D+00,  3.0D+00,  4.0D+00, 99.0D+00, &
    99.0D+00,  0.0D+00, 99.0D+00,  1.0D+00,  7.0D+00,  4.0D+00, &
     3.0D+00, 99.0D+00,  0.0D+00, 99.0D+00,  2.0D+00,  5.0D+00, &
     3.0D+00,  1.0D+00, 99.0D+00,  0.0D+00, 99.0D+00,  2.0D+00, &
     4.0D+00,  7.0D+00,  2.0D+00, 99.0D+00,  0.0D+00, 99.0D+00, &
    99.0D+00,  4.0D+00,  5.0D+00,  2.0D+00, 99.0D+00,  0.0D+00 /

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST24'
  write ( *, '(a)' ) '  GRAPH_DIST_MINTR1 finds the minimum spanning tree:'

  call graph_dist_print ( dist, lda, nnode, '  The distance matrix:' )

  call graph_dist_mintr1 ( nnode, dist, lda, itree, jtree )

  do i = 1, nnode-1
    wtree(i) = dist(itree(i),jtree(i))
  end do

  call graph_arc_weight_print ( nnode-1, itree, jtree, wtree, &
    '  The weighted spanning tree:' )

  return
end
