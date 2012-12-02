program main

!*****************************************************************************80
!
!! MAIN is the main program for SEARCH_MPI.
!
!  Discussion:
!
!    SEARCH demonstrates the use of MPI routines to carry out a search.
!
!    An array of given size is to be searched for occurrences of a
!    specific value.
!
!    The search is done in parallel.  A master process generates the
!    array and the target value, then distributes the information among
!    a set of worker processes, and waits for them to communicate back
!    the (global) index values at which occurrences of the target value
!    were found.
!
!    An interesting feature of this program is the use of allocatable
!    arrays, which allows the master program to set aside just enough
!    memory for the whole array, and for each worker program to set aside
!    just enough memory for its own part of the array.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 October 2002
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

  integer, allocatable, dimension ( : ) :: a
  integer dest
  real factor
  integer global
  integer i
  integer id
  integer ierr
  integer n
  integer npart
  integer num_procs
  real, allocatable, dimension ( : ) :: r
  integer source
  integer start
  integer status(MPI_STATUS_SIZE)
  integer tag
  integer, parameter :: tag_target = 1
  integer, parameter :: tag_size = 2
  integer, parameter :: tag_data = 3
  integer, parameter :: tag_found = 4
  integer, parameter :: tag_done = 5
  integer target
  integer workers_done
  integer x
!
!  Initialize MPI.
!
  call MPI_Init ( ierr )
!
!  Get this processes's rank.
!
  call MPI_Comm_rank ( MPI_COMM_WORLD, id, ierr )
!
!  Find out how many processes are available.
!
  call MPI_Comm_size ( MPI_COMM_WORLD, num_procs, ierr )

  if ( id == 0 ) then
    call timestamp ( )
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SEARCH_MPI:'
    write ( *, '(a)' ) '  FORTRAN90 version'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  An MPI example program to search an array.'
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  The number of processes is ', num_procs
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a)' ) 'Process ', id, ' is active.'
!
!  Have the master process generate the target and data.  In a more 
!  realistic application, the data might be in a file which the master 
!  process would read.  Here, the master process decides
!
  if ( id == 0 ) then
!
!  Pick the number of data items per process, and set the total.
!
    call random_number ( factor )
    npart = 50 + int ( factor * 100.0E+00 )
    n = npart * num_procs

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SEARCH_MPI - Master process:'
    write ( *, '(a,i8)' ) '  The number of data items per process is ', npart
    write ( *, '(a,i8)' ) '  The total number of data items is       ', n
!
!  Now allocate the master copy of A, fill it with values, and pick 
!  a value for the target.
!
    allocate ( a(1:n) )

    allocate ( r(1:n) )

    call random_number ( r(1:n) )

    factor = real ( n ) / 10
    a(1:n) = int ( factor * r(1:n) )

    deallocate ( r )

    target = a(n/2)

    write ( *, '(a,i8)' ) '  The target value is ', target
!
!  The worker processes need to have the target value, the number of data items,
!  and their individual chunk of the data vector.
!
    do i = 1, num_procs-1

      dest = i
      tag = tag_target

      call MPI_Send ( target, 1, MPI_INTEGER, dest, tag, MPI_COMM_WORLD, &
        ierr )

      tag = tag_size

      call MPI_Send ( npart, 1, MPI_INTEGER, dest, tag, MPI_COMM_WORLD, &
        ierr )

      start = ( i - 1 ) * npart + 1
      tag = tag_data

      call MPI_Send ( a(start), npart, MPI_INTEGER, dest, tag, &
        MPI_COMM_WORLD, ierr )

    end do
!
!  Now the master process simply waits for each worker process to report that 
!  it is done.
!
    workers_done = 0

    do while ( workers_done < num_procs-1 )

      call MPI_Recv ( x, 1, MPI_INTEGER, MPI_ANY_SOURCE, MPI_ANY_TAG, &
        MPI_COMM_WORLD, status, ierr )

      source = status(MPI_SOURCE)
      tag = status(MPI_TAG)
    
      if ( tag == tag_done ) then

        workers_done = workers_done + 1

      else if ( tag == tag_found ) then

        write ( *, '(a,i2,2i8)' ) 'P', source, x, a(x)

      else

        write ( *, '(a,i8)' ) &
          '  Master process received message with unknown tag = ', tag

      end if

    end do
!
!  The master process can throw away A now.
!
    deallocate ( a )
!
!  Each worker process expects to receive the target value, the number of data
!  items, and the data vector.
!
  else 

    tag = tag_target

    call MPI_Recv ( target, 1, MPI_INTEGER, 0, tag, MPI_COMM_WORLD, &
      status, ierr )
 
    tag = tag_size

    call MPI_Recv ( npart, 1, MPI_INTEGER, 0, tag, MPI_COMM_WORLD, &
      status, ierr )

    allocate ( a(1:npart) )

    tag = tag_data

    call MPI_Recv ( a, npart, MPI_INTEGER, 0, tag, MPI_COMM_WORLD, &
      status, ierr )
!
!  The worker simply checks each entry to see if it is equal to the target
!  value.
!
    do i = 1, npart

      if ( a(i) == target ) then

        global = ( id - 1 ) * npart + i
        tag = tag_found

        call MPI_Send ( global, 1, MPI_INTEGER, 0, tag, MPI_COMM_WORLD, &
          ierr )

      end if

    end do  
!
!  When the worker is finished with the loop, it sends a dummy data value with
!  the tag "TAG_DONE" indicating that it is done.
!
    tag = tag_done

    call MPI_Send ( target, 1, MPI_INTEGER, 0, tag, MPI_COMM_WORLD, ierr )

    deallocate ( a )
     
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
    write ( *, '(a)' ) 'SEARCH_MPI:'
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
