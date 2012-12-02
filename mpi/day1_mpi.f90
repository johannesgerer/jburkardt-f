program main

!*****************************************************************************80
!
!! MAIN is the main program for DAY1.
!
!  Discussion:
!
!    DAY1_EX3 is exercise 3 for first day of the MPI workshop.
!
!    The instructions say:
!
!    Process 1 computes the squares of the first 200 integers.
!    It sends this data to process 3.
!
!    Process 3 should divide the integers between 20 and 119 by 53,
!    getting a real result, and passes this data back to process 1.
!
!    * I presume the first 200 integers are the numbers 0 through 199.
!
!    * The instructions literally mean that process 3 should look
!      at integers whose VALUES are between 20 and 119.  I doubt that
!      is what the instructor meant, but it's more interesting than
!      simply picking the entries with index between 20 and 119,
!      so that's what I'll do.
!
!    * It is also not completely clear whether only the selected data
!      should be sent back, or the entire array.  Again, it is more
!      interesting to send back only part of the data.
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

  integer, parameter :: i_dim = 200
  integer, parameter :: r_dim = 200

  integer count
  integer count2
  integer dest
  integer i
  integer i_buffer(i_dim)
  integer id
  integer ierr
  integer num_procs
  real r_buffer(r_dim)
  integer source
  integer status(MPI_Status_size)
  integer tag
!
!  Initialize MPI.
!
  call MPI_Init ( ierr )
!
!  Determine this process's id.
!
  call MPI_Comm_rank ( MPI_COMM_WORLD, id, ierr )
!
!  Have Process 0 say hello.
!
  if ( id == 0 ) then
    call timestamp ( )
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DAY1:'
    write ( *, '(a)' ) '  FORTRAN90 version'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  MPI exercise #3 for day 1.'
    write ( *, '(a,i8)' ) '  The number of processes available is ', num_procs
  end if
!
!  Get the number of processes.
!
  call MPI_Comm_size ( MPI_COMM_WORLD, num_procs, ierr )
!
!  If we don't have at least 4 processes, then bail out now.
!
  if ( num_procs < 4 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) 'DAY1 - Process ', id
    write ( *, '(a)' ) '  Not enough processes for this task!'
    write ( *, '(a)' ) '  Bailing out now!'
    call MPI_Finalize ( ierr )
    stop
  end if
!
!  Process 1 knows that it will generate 200 integers, and may receive no more
!  than 200 reals.
!
  if ( id == 1 ) then

    count = 200

    do i = 1, count
      i_buffer(i) = ( i - 1 )**2
    end do

    dest = 3
    tag = 1

    call MPI_Send ( i_buffer, count, MPI_INTEGER, dest, tag, &
      MPI_COMM_WORLD, ierr )

    write ( *, '(a,i1,a,i3,a,i1)' ) 'P:', id, ' sent ', count, &
      ' integers to process ', dest

    source = 3
    tag = 2

    call MPI_Recv ( r_buffer, r_dim, MPI_REAL, source, tag, &
      MPI_COMM_WORLD, status, ierr )

    write ( *, '(a,i1,a,i1)' ) 'P:', id, &
      ' received real values from process 3.'

    call MPI_Get_count ( status, MPI_REAL, count, ierr )

    write ( *, '(a,i1,a,i3,a)' ) 'P:', id, &
      ' Number of real values received is ', count

    write ( *, '(a,i1,a,3f8.4)' ) 'P:', id, &
      ' First 3 values = ', r_buffer(1:3)
!
!  Process 3 receives the integer data from process 1, selects some 
!  of the data, does a real computation on it, and sends that part 
!  back to process 1.
!
  else if ( id == 3 ) then
 
    source = 1
    tag = 1

    call MPI_Recv ( i_buffer, i_dim, MPI_INTEGER, source, tag, &
      MPI_COMM_WORLD, status, ierr )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i1,a)' ) 'P:', id, &
      ' received integer values from process 1.'

    call MPI_Get_count ( status, MPI_INTEGER, count, ierr )

    write ( *, '(a,i1,a,i8)' ) 'P:', id, &
      ' - Number of integers received is ', count
    write ( *, '(a,i1,a,3i8)' ) 'P:', id, ' First 3 values = ', i_buffer(1:3)

    count2 = 0

    do i = 1, count
     
      if ( 20 <= i_buffer(i) .and. i_buffer(i) <= 119 ) then

        count2 = count2 + 1
        r_buffer(count2) = real ( i_buffer(i) ) / 53.0E+00

        if ( count2 <= 3 ) then
          write ( *, '(a,i1,a,i8,a,f8.4)' ) 'P:', id, ' Input integer ', &
            i_buffer(i), ' becomes ', r_buffer(count2)
        end if

      end if

    end do

    dest = 1
    tag = 2
  
    call MPI_Send ( r_buffer, count2, MPI_REAL, dest, tag, &
      MPI_COMM_WORLD, ierr )

    write ( *, '(a,i1,a,i3,a,i1)' ) 'P:', id, ' sent ', count2, &
      ' reals to process ', dest

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a,i1,a)' ) 'P:', id, ' - MPI has no work for me!'

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
    write ( *, '(a)' ) 'DAY1:'
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
