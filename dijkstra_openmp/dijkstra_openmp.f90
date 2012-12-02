program main

!*****************************************************************************80
!
!! MAIN runs an example of Dijkstra's minimum distance algorithm.
!
!  Discussion:
!
!    Given the distance matrix that defines a graph, we seek a list
!    of the minimum distances between node 0 and all other nodes.
!
!    This program sets up a small example problem and solves it.
!
!    The correct minimum distances are:
!
!      0   35   15   45   49   41
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 July 2010
!
!  Author:
!
!    Original C version by Norm Matloff, CS Dept, UC Davis.
!    FORTRAN90 version by John Burkardt.
!
  use omp_lib

  implicit none

  integer ( kind = 4 ), parameter :: nv = 6

  integer ( kind = 4 ) i
  integer ( kind = 4 ) :: i4_huge = 2147483647
  integer ( kind = 4 ) j
  integer ( kind = 4 ) mind(nv)
  integer ( kind = 4 ) ohd(nv,nv)

  call timestamp ( )
  write ( *, '(a)' )  ' '
  write ( *, '(a)' )  'DIJKSTRA_OPENMP:'
  write ( *, '(a)' )  '  FORTRAN90 version'
  write ( *, '(a)' )  '  Use Dijkstra''s algorithm to determine the minimum'
  write ( *, '(a)' )  '  distance from node 1 to each node in a graph,'
  write ( *, '(a)' )  '  given the distances between each pair of nodes.'
  write ( *, '(a)' )  ' '
  write ( *, '(a)' )  '  Although a very small example is considered, we'
  write ( *, '(a)' )  '  demonstrate the use of OpenMP directives for'
  write ( *, '(a)' )  '  parallel execution.'
!
!  Initialize the problem data.
!
  call init ( nv, ohd )
!
!  Print the distance matrix.
!
  write ( *, '(a)' )  ' '
  write ( *, '(a)' )  '  Distance matrix:'
  write ( *, '(a)' )  ' '
  do i = 1, nv
    do j = 1, nv
      if ( ohd(i,j) == i4_huge ) then
        write ( *, '(2x,a)', advance = 'NO' ) 'Inf'
      else
        write ( *, '(2x,i3)', advance = 'NO' ) ohd(i,j)
      end if
    end do
    write ( *, '(a)', advance = 'yes' )
  end do
!
!  Carry out the algorithm.
!
  call dijkstra_distance ( nv, ohd, mind )
!
!  Print the results.
!
  write ( *, '(a)' )  ' '
  write ( *, '(a)' )  '  Minimum distances from node 1:'
  write ( *, '(a)' )  ' '
  do i = 1, nv
    write ( *, '(2x,i2,2x,i2)' ) i, mind(i)
  end do
!
!  Terminate.
!
  write ( *, '(a)' )  ' '
  write ( *, '(a)' )  'DIJKSTRA_OPENMP:'
  write ( *, '(a)' )  '  Normal end of execution.'

  write ( *, '(a)' )  ' '
  call timestamp ( )

  stop
end
subroutine dijkstra_distance ( nv, ohd, mind )

!*****************************************************************************80
!
!! DIJKSTRA_DISTANCE uses Dijkstra's minimum distance algorithm.
!
!  Discussion:
!
!    We essentially build a tree.  We start with only node 0 connected
!    to the tree, and this is indicated by setting CONNECTED(0) = TRUE.
!
!    We initialize MIND(I) to the one step distance from node 0 to node I.
!    
!    Now we search among the unconnected nodes for the node MV whose minimum
!    distance is smallest, and connect it to the tree.  For each remaining
!    unconnected node I, we check to see whether the distance from 0 to MV
!    to I is less than that recorded in MIND(I), and if so, we can reduce
!    the distance.
!
!    After NV-1 steps, we have connected all the nodes to 0, and computed
!    the correct minimum distances.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 July 2010
!
!  Author:
!
!    Original C version by Norm Matloff, CS Dept, UC Davis.
!    This FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NV, the number of nodes.
!
!    Input, integer ( kind = 4 ) OHD(NV,NV), the distance of the direct
!    link between nodes I and J.
!
!    Output, integer ( kind = 4 ) MIND(NV), the minimum 
!    distance from node 1 to each node.
!
  use omp_lib

  implicit none

  integer ( kind = 4 ) nv

  logical connected(nv)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) :: i4_huge = 2147483647
  integer ( kind = 4 ) md
  integer ( kind = 4 ) mind(nv)
  integer ( kind = 4 ) mv
  integer ( kind = 4 ) my_first
  integer ( kind = 4 ) my_id
  integer ( kind = 4 ) my_last
  integer ( kind = 4 ) my_md
  integer ( kind = 4 ) my_mv
  integer ( kind = 4 ) my_step
  integer ( kind = 4 ) nth
  integer ( kind = 4 ) ohd(nv,nv)
!
!  Start out with only node 1 connected to the tree.
!
  connected(1) = .true.
  connected(2:nv) = .false.
!
!  Initialize the minimum distance to the one-step distance.
!
  mind(1:nv) = ohd(1,1:nv)
!
!  Begin the parallel region.
!
!$omp parallel private ( my_first, my_id, my_last, my_md, my_mv, my_step ) &
!$omp          shared ( connected, md, mind, mv, nth, ohd )

  my_id = omp_get_thread_num ( )
  nth = omp_get_num_threads ( ) 
  my_first =   (   my_id       * nv ) / nth + 1
  my_last  =   ( ( my_id + 1 ) * nv ) / nth
!
!  The SINGLE directive means that the block is to be executed by only
!  one thread, and that thread will be whichever one gets here first.
!
!$omp single
  write ( *, '(a)' ) ' '
  write ( *, * ) '  P', my_id, &
    ': Parallel region begins with ', nth, ' threads.'
  write ( *, '(a)' ) ' '
!$omp end single
  write ( *, * ) '  P', my_id, ':  First=', my_first, '  Last=', my_last
!
!  Attach one more node on each iteration.
!
  do my_step = 2, nv
!
!  Before we compare the results of each thread, set the shared variable 
!  MD to a big value.  Only one thread needs to do this.
!
!$omp single 
    md = i4_huge
    mv = -1
!$omp end single
!
!  Each thread finds the nearest unconnected node in its part of the graph.
!  Some threads might have no unconnected nodes left.
!
    call find_nearest ( my_first, my_last, nv, mind, connected, my_md, my_mv )
!
!  In order to determine the minimum of all the MY_MD's, we must insist
!  that only one thread at a time execute this block!
!
!$omp critical
    if ( my_md < md ) then
      md = my_md
      mv = my_mv
    end if
!$omp end critical
!
!  This barrier means that ALL threads have executed the critical
!  block, and therefore MD and MV have the correct value.  Only then
!  can we proceed.
!
!$omp barrier
!
!  If MV is -1, then NO thread found an unconnected node, so we're done early. 
!  OpenMP does not like to BREAK out of a parallel region, so we'll just have 
!  to let the iteration run to the end, while we avoid doing any more updates.
!
!  Otherwise, we connect the nearest node.
!
!$omp single
    if ( mv /= - 1 ) then
      connected(mv) = .true.
      write ( *, * ) '  P', my_id, ': Connecting node ', mv
    else
      write ( *, * ) '  P', my_id, ': No connecting node on step ', my_step
    end if
!$omp end single
!
!  Again, we don't want any thread to proceed until the value of
!  CONNECTED is updated.
!
!$omp barrier
!
!  Now each thread should update its portion of the MIND vector,
!  by checking to see whether the trip from 0 to MV plus the step
!  from MV to a node is closer than the current record.
!
    if ( mv /= -1 ) then
      call update_mind ( my_first, my_last, nv, connected, ohd, mv, mind )
    end if
!
!  Before starting the next step of the iteration, we need all threads 
!  to complete the updating, so we set a BARRIER here.
!
!$omp barrier

  end do
!
!  Once all the nodes have been connected, we can exit.
!
!$omp single
  write ( *, * ) ' '
  write ( *, * ) '  P', my_id, ': Exiting parallel region.'
!$omp end single

!$omp end parallel

  return
end
subroutine find_nearest ( s, e, nv, mind, connected, d, v )

!*****************************************************************************80
!
!! FIND_NEAREST finds the nearest unconnected node.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 July 2010
!
!  Author:
!
!    Original C version by Norm Matloff, CS Dept, UC Davis.
!    This FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) S, E, the first and last nodes that 
!    are to be checked.
!
!    Input, integer ( kind = 4 ) NV, the number of nodes.
!
!    Input, integer ( kind = 4 ) MIND(NV), the currently computed minimum 
!    distance from node 1 to each node.
!
!    Input, logical CONNECTED(NV), is true for each connected node, whose 
!    minimum distance to node 1 has been determined.
!
!    Output, integer ( kind = 4 ) D, the distance from node 1 to the nearest 
!    unconnected node in the range S to E.
!
!    Output, integer ( kind = 4 ) V, the index of the nearest unconnected node
!    in the range S to E.
!
  implicit none

  integer ( kind = 4 ) nv

  logical connected(nv)
  integer ( kind = 4 ) d
  integer ( kind = 4 ) e
  integer ( kind = 4 ) i
  integer ( kind = 4 ) :: i4_huge = 2147483647
  integer ( kind = 4 ) mind(nv)
  integer ( kind = 4 ) s
  integer ( kind = 4 ) v

  d = i4_huge
  v = -1

  do i = s, e
    if ( .not. connected(i) .and. mind(i) < d ) then
      d = mind(i)
      v = i
    end if
  end do

  return
end
subroutine init ( nv, ohd )

!*****************************************************************************80
!
!! INIT initializes the problem data.
!
!  Discussion:
!
!    The graph uses 6 nodes, and has the following diagram and
!    distance matrix:
!
!    N0--15--N2-100--N3           0   40   15  Inf  Inf  Inf
!      \      |     /            40    0   20   10   25    6
!       \     |    /             15   20    0  100  Inf  Inf
!        40  20  10             Inf   10  100    0  Inf  Inf
!          \  |  /              Inf   25  Inf  Inf    0    8
!           \ | /               Inf    6  Inf  Inf    8    0
!            N1
!            / \
!           /   \
!          6    25
!         /       \
!        /         \
!      N5----8-----N4
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 July 2010
!
!  Author:
!
!    Original C version by Norm Matloff, CS Dept, UC Davis.
!    This FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NV, the number of nodes.
!
!    Output, integer ( kind = 4 ) OHD(NV,NV), the distance of the direct
!    link between nodes I and J.
!
  implicit none

  integer ( kind = 4 ) nv

  integer ( kind = 4 ) i
  integer ( kind = 4 ) :: i4_huge = 2147483647
  integer ( kind = 4 ) ohd(nv,nv)

  ohd(1:nv,1:nv) = i4_huge

  do i = 1, nv
    ohd(i,i) = 0
  end do

  ohd(1,2) = 40
  ohd(1,3) = 15
  ohd(2,3) = 20
  ohd(2,4) = 10
  ohd(2,5) = 25
  ohd(3,4) = 100
  ohd(2,6) = 6
  ohd(5,6) = 8

  ohd(2,1) = 40
  ohd(3,1) = 15
  ohd(3,2) = 20
  ohd(4,2) = 10
  ohd(5,2) = 25
  ohd(4,3) = 100
  ohd(6,2) = 6
  ohd(6,5) = 8

  return
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
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
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
  integer   ( kind = 4 ) d
  integer   ( kind = 4 ) h
  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) mm
  character ( len = 9  ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer   ( kind = 4 ) n
  integer   ( kind = 4 ) s
  integer   ( kind = 4 ) values(8)
  integer   ( kind = 4 ) y

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
subroutine update_mind ( s, e, nv, connected, ohd, mv, mind )

!*****************************************************************************80
!
!! UPDATE_MIND updates the minimum distance vector.
!
!  Discussion:
!
!    We've just determined the minimum distance to node MV.
!
!    For each node I which is not connected yet,
!    check whether the route from node 0 to MV to I is shorter
!    than the currently known minimum distance.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 July 2010
!
!  Author:
!
!    Original C version by Norm Matloff, CS Dept, UC Davis.
!    This FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) S, E, the first and last nodes that are 
!    to be checked.
!
!    Input, integer ( kind = 4 ) NV, the number of nodes.
!
!    Input, logical CONNECTED(NV), is true for each connected node, whose 
!    minimum distance to node 0 has been determined.
!
!    Input, integer ( kind = 4 ) OHD(NV,NV), the distance of the direct link 
!    between nodes I and J.
!
!    Input, integer ( kind = 4 ) MV, the node whose minimum distance to node 20
!    has just been determined.
!
!    Input/output, integer ( kind = 4 ) MIND(NV), the currently computed
!    minimum distances from node 1 to each node.  On output, the values for 
!    nodes S through E have been updated.
!
  implicit none

  integer ( kind = 4 ) nv

  logical              connected(nv)
  integer ( kind = 4 ) e
  integer ( kind = 4 ) i
  integer ( kind = 4 ) :: i4_huge = 2147483647
  integer ( kind = 4 ) mind(nv)
  integer ( kind = 4 ) mv
  integer ( kind = 4 ) ohd(nv,nv)
  integer ( kind = 4 ) s

  do i = s, e
    if ( .not. connected(i) ) then
      if ( ohd(mv,i) < i4_huge ) then
        if ( mind(mv) + ohd(mv,i) < mind(i) ) then
          mind(i) = mind(mv) + ohd(mv,i)
        end if
      end if
    end if
  end do

  return
end



