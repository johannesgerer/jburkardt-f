program main

!*****************************************************************************80
!
!! TYPE_MPI demonstrates the creation of a user-defined type in MPI.
!
!  Discussion:
!
!    TYPE demonstrates user-defined datatypes in MPI.
!
!    The datatype defined will be a structure that contains three integers.
!
!    Process 0 will set up an example of this structure, and send it
!    to Proces 1, which will alter it and send it back.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2005
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
!    ISBN: 0262571323.
!
  use mpi

  type point
    integer :: x
    integer :: y
    integer :: z
  end type

  integer dest
  integer i
  integer id
  integer ierr
  integer num_procs
  type ( point ) :: my_point
  integer point_type
  integer source
  integer status(MPI_STATUS_SIZE)
  integer tag
!
!  Initialize MPI.
!
  call MPI_Init ( ierr )
!
!  Get the number of processes.
!
  call MPI_Comm_size ( MPI_COMM_WORLD, num_procs, ierr )
!
!  Get the individual process ID.
!
  call MPI_Comm_rank ( MPI_COMM_WORLD, id, ierr )
!
!  Print a message.
!
  if ( id == 0 ) then
    call timestamp ( )
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TYPE_MPI:'
    write ( *, '(a)' ) '  FORTRAN90 version'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  An MPI example program to set up and '
    write ( *, '(a)' ) '  use an MPI datatype.'
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  The number of processes is ', num_procs
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a)' ) '  Process ', id, ' is active.'
!
!  Define and commit the new datatype.
!
  call MPI_Type_contiguous ( 3, MPI_INTEGER, point_type, ierr )

  call MPI_Type_commit ( point_type, ierr )

  if ( id == 0 ) then

    my_point%x = 1
    my_point%y = 2
    my_point%z = 4
    dest = 1
    tag = 1

    call MPI_Send ( my_point, 1, point_type, dest, tag, MPI_COMM_WORLD, ierr )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8,a,3i8)' ) '  Process ', id, &
      ' sent an item of type POINT_TYPE, with value ', &
      my_point%x, my_point%y, my_point%z

    source = 1
    tag = 2
    call MPI_Recv ( my_point, 1, point_type, source, tag, MPI_COMM_WORLD, &
      status, ierr )

    write ( *, '(a,i8,a,3i8)' ) '  Process ', id, &
      ' received a modified item of type POINT_TYPE, with value ', & 
      my_point%x, my_point%y, my_point%z

  else if ( id == 1 ) then

    source = 0
    tag = 1

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8,a,3i8)' ) '  Process ', id, &
      ' expecting an item of type POINT_TYPE.'
 
    call MPI_Recv ( my_point, 1, point_type, source, tag, MPI_COMM_WORLD, &
      status, ierr )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8,a,3i8)' ) '  Process ', id, &
      ' received an item of type POINT_TYPE, with value ', &
      my_point%x, my_point%y, my_point%z
 
    i = my_point%x
    my_point%x = my_point%z * 100
    my_point%y = my_point%y * 10
    my_point%z = i
    dest = 0
    tag = 2
    
    call MPI_Send ( my_point, 1, point_type, dest, tag, MPI_COMM_WORLD, ierr )

    write ( *, '(a,i8,a,3i8)' ) '  Process ', id, &
      ' sent a modified item of type POINT_TYPE, with value ', &
      my_point%x, my_point%y, my_point%z

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8,a)' ) '  Process ', id, &
      ': MPI has nothing for me to do!'

  end if
!
!  Terminate MPI.
!
  call MPI_Finalize ( ierr )
!
!  Terminate.
!
  if ( id == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TYPE_MPI:'
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

  character ( len = 8  ) ampm
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9  ), parameter, dimension(12) :: month = (/ &
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
