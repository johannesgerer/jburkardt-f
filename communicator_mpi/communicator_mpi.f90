program main

!*****************************************************************************80
!
!! MAIN is the main program for COMMUNICATOR_MPI.
!
!  Discussion:
!
!    This program demonstrates how an MPI program can start with the
!    default communicator MPI_COMM_WORLD, and create new communicators
!    referencing a subset of the total number of processes.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gropp, Ewing Lusk, Anthony Skjellum,
!    Using MPI: Portable Parallel Programming with the
!    Message-Passing Interface,
!    Second Edition,
!    MIT Press, 1999,
!    ISBN: 0262571323,
!    LC: QA76.642.G76.
!
  use mpi

  integer ( kind = 4 ) even_comm_id
  integer ( kind = 4 ) even_group_id
  integer ( kind = 4 ) even_id
  integer ( kind = 4 ) even_id_sum
  integer ( kind = 4 ) even_p
  integer ( kind = 4 ), allocatable :: even_rank(:)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) id
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) j
  integer ( kind = 4 ) odd_comm_id
  integer ( kind = 4 ) odd_group_id
  integer ( kind = 4 ) odd_id
  integer ( kind = 4 ) odd_id_sum
  integer ( kind = 4 ) odd_p
  integer ( kind = 4 ), allocatable :: odd_rank(:)
  integer ( kind = 4 ) p
  integer ( kind = 4 ) world_group_id
!
!  Initialize MPI.
!
  call MPI_Init ( ierr )
!
!  Get the number of processes.
!
  call MPI_Comm_size ( MPI_COMM_WORLD, p, ierr )
!
!  Get the individual process ID.
!
  call MPI_Comm_rank ( MPI_COMM_WORLD, id, ierr )
!
!  Process 0 prints an introductory message.
!
  if ( id == 0 ) then
    call timestamp ( )
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'COMMUNICATOR_MPI - Master process:'
    write ( *, '(a)' ) '  FORTRAN90/MPI version'
    write ( *, '(a)' ) '  An MPI example program.'
    write ( *, '(a)' ) ' '
    write ( *, '(a,i4)' ) '  The number of processes is ', p
    write ( *, '(a)' ) ' '
  end if
!
!  Every process prints a hello.
!
  write ( *, '(a,i4,a)' ) '  Process ', id, ' says "Hello, world!".'
!
!  Get a group identifier for MPI_COMM_WORLD.
!
  call MPI_Comm_group ( MPI_COMM_WORLD, world_group_id, ierr )
!
!  List the even processes, and create their group.
!
  even_p = ( p + 1 ) / 2
  allocate ( even_rank(1:even_p) )
  j = 0
  do i = 0, p - 1, 2
    j = j + 1
    even_rank(j) = i
  end do
  call MPI_Group_incl ( world_group_id, even_p, even_rank, even_group_id, ierr )

  call MPI_Comm_create ( MPI_COMM_WORLD, even_group_id, even_comm_id, ierr )
!
!  List the odd processes, and create their group.
!
  odd_p = p / 2
  allocate ( odd_rank(1:odd_p) )
  j = 0
  do i = 1, p - 1, 2
    j = j + 1
    odd_rank(j) = i
  end do
  call MPI_Group_incl ( world_group_id, odd_p, odd_rank, odd_group_id, ierr )

  call MPI_Comm_create ( MPI_COMM_WORLD, odd_group_id, odd_comm_id, ierr )
!
!  Try to get ID of each process in both groups.  
!  If a process is not in a communicator, set its ID to -1.
!
  if ( mod ( id, 2 ) == 0 ) then
    call MPI_Comm_rank ( even_comm_id, even_id, ierr )
    odd_id = -1
  else
    call MPI_Comm_rank ( odd_comm_id,  odd_id, ierr )
    even_id = -1
  end if
!
!  Use MPI_Reduce to sum the global ID of each process in the even group.
!  Assuming 4 processes: EVEN_SUM = 0 + 2 = 2
!
  if ( even_id /= -1 ) then
    call MPI_Reduce ( id, even_id_sum, 1, MPI_INTEGER, MPI_SUM, 0, &
      even_comm_id, ierr )
  end if

  if ( even_id == 0 ) then
    write ( *, '(a,i4)' ) &
      '  Number of processes in even communicator = ', even_p
    write ( *, '(a,i4)' ) &
      '  Sum of global ID''s in even communicator = ', even_id_sum
  end if
!
!  Use MPI_Reduce to sum the global ID of each process in the odd group.
!  Assuming 4 processes: ODD_SUM = 1 + 3 = 4
!
  if ( odd_id /= -1 ) then
    call MPI_Reduce ( id, odd_id_sum,  1, MPI_INTEGER, MPI_SUM, 0, &
      odd_comm_id, ierr )
  end if

  if ( odd_id == 0 ) then
    write ( *, '(a,i4)' ) &
      '  Number of processes in odd communicator  = ', odd_p
    write ( *, '(a,i4)' ) &
      '  Sum of global ID''s in odd communicator  = ', odd_id_sum
  end if
!
!  Terminate MPI.
!
  call MPI_Finalize ( ierr )
!
!  Free memory.
!
  deallocate ( even_rank )
  deallocate ( odd_rank )
!
!  Terminate
!
  if ( id == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'COMMUNICATOR_MPI:'
    write ( *, '(a)' ) '  Normal end of execution.'
    write ( *, '(a)' ) ' '
    call timestamp ( )
  end if

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

