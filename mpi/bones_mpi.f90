program main

!*****************************************************************************80
!
!! MAIN is the main program for BONES.
!
!  Discussion:
!
!    BONES is a simple demonstration of the use of MPI by a FORTRAN90 program.
!
!    This program should be run on at least two processes.
!    Only the first two processes will get any work to do.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 October 2002
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

  integer count
  real data(0:99)
  integer dest
  integer i
  integer ierr
  integer num_procs
  integer rank
  integer status(MPI_Status_size)
  integer tag
  real value(200)
!
!  Initialize MPI.
!
  call MPI_Init ( ierr )
!
!  Determine this process's rank.
!
  call MPI_Comm_rank ( MPI_COMM_WORLD, rank, ierr )
!
!  Find out the number of processes available.
!
  call MPI_Comm_size ( MPI_COMM_WORLD, num_procs, ierr )
!
!  Have Process 0 say hello.
!
  if ( rank == 0 ) then
    call timestamp ( )
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BONES:'
    write ( *, '(a)' ) '  FORTRAN90 version'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  An MPI test program.'
    write ( *, '(a,i8)' ) '  The number of processes available is ', num_procs

  end if
!
!  Process 0 expects to receive as much as 200 real values, from any source.
!
  if ( rank == 0 ) then

    tag = 55
    call MPI_Recv ( value, 200, MPI_REAL, MPI_ANY_SOURCE, tag, &
      MPI_COMM_WORLD, status, ierr )

    write ( *, '(a,i1,a,i1)' ) 'P:', rank, ' Got data from processor ', &
      status(MPI_SOURCE)

    call MPI_Get_count ( status, MPI_REAL, count, ierr )

    write ( *, '(a,i1,a,i3,a)' ) 'P:', rank, ' Got ', count, ' elements.'

    write ( *, '(a,i1,a,g14.6)' ) 'P:', rank, ' value(5) = ', value(5)
!
!  Process 1 sends 100 real values to process 0.
!
  else if ( rank == 1 ) then
 
    write ( *, '(a)' ) ' '
    write ( *, '(a,i1,a)' ) 'P:', rank, &
      ' - setting up data to send to process 0.'

    do i = 0, 99
      data(i) = real ( i )
    end do

    dest = 0
    tag = 55
    call MPI_Send ( data, 100, MPI_REAL, dest, tag, MPI_COMM_WORLD, ierr )
  
  else

    write ( *, '(a)' ) ' '
    write ( *, '(a,i1,a)' ) 'P:', rank, ' - MPI has no work for me!'

  end if

  call MPI_Finalize ( ierr )

  if ( rank == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BONES:'
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
