program main

!*****************************************************************************80
!
!! MAIN is the main program for MONTE_CARLO.
!
!  Discussion:
!
!    MONTE_CARLO computes PI using Monte Carlo techniques.
!
!    Generate N random points in the unit square.  Count M, the number of
!    points that are in the quarter circle.  Then PI is equal to the
!    ratio 4 * M / N.
!
!    It's important that each processor use DIFFERENT random numbers.
!    One way to ensure this is to have a single master processor
!    generate all the random numbers, and then divide them up.
!
!    (A second way, not explored here, is simply to ensure that each
!    processor uses a different seed, either chosen by a master processor,
!    or generated from the processor ID.)
!
!    The work will be divided as follows:
!
!      Process 0 is in charge.
!
!      Processes 0 through NUM_PROCS - 2 take sets of random points and count
!      the number of points in the quarter circle.
!
!      Process NUM_PROCS - 1 computes the sets of random points.
!
!    The communicator MPI_WORLD_COMM comprises all processes.
!
!    The communicator WORKER_COMM is created to comprise the processes
!    whose world ID's are 1 through NUM_PROCS-2, that is, all but the
!    random number generating process.
!
!    Message tag 0 is a signal from the master process to exit the loop.
!    Message tag 1 is a request for a set of random numbers.
!    Message tag 2 is a reply from the server with a set of random numbers.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 September 2002
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

  integer ( kind = 4 ), parameter :: chunk_size = 1000

  integer ( kind = 4 ) done
  integer ( kind = 4 ) :: dummy = 0
  real ( kind = 8 ) error
  integer ( kind = 4 ) i
  integer ( kind = 4 ) id
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) in_local
  integer ( kind = 4 ) in_total
  integer ( kind = 4 ) num_procs
  integer ( kind = 4 ) out_local
  integer ( kind = 4 ) out_total
  real ( kind = 8 ) pi_est
  real ( kind = 8 ), parameter :: pi_true = 3.141592653589793238462643D+00
  integer ( kind = 4 ) point_local
  integer ( kind = 4 ) :: point_max = 1000000
  integer ( kind = 4 ) point_total
  real ( kind = 8 ) rands(chunk_size)
  integer ( kind = 4 ) ranks(1)
  integer ( kind = 4 ) requestor
  integer ( kind = 4 ) server
  integer ( kind = 4 ) status(MPI_STATUS_SIZE)
  integer ( kind = 4 ) tag
  integer ( kind = 4 ), parameter :: tag_exit = 0
  integer ( kind = 4 ), parameter :: tag_random_send = 1
  integer ( kind = 4 ), parameter :: tag_random_sent = 2
  real ( kind = 8 ) tolerance
  integer ( kind = 4 ) worker_comm
  integer ( kind = 4 ) worker_group
  integer ( kind = 4 ) world_group
  real ( kind = 8 ) x
  real ( kind = 8 ) y
!
!  Initialize MPI.
!
  call MPI_Init ( ierr )
!
!  Get the number of processors.
!
  call MPI_Comm_size ( MPI_COMM_WORLD, num_procs, ierr )
!
!  Get the rank of this processor.
!
  call MPI_Comm_rank ( MPI_COMM_WORLD, id, ierr )
!
!  The master process gets the value of the tolerance...
!
  if ( id == 0 ) then
    call timestamp ( )
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MONTE_CARLO_MPI:'
    write ( *, '(a)' ) '  FORTRAN90 version'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  An MPI example program to'
    write ( *, '(a)' ) '  estimate PI by the Monte Carlo method.'
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  The number of processes is ', num_procs
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Points in the unit square will be tested'
    write ( *, '(a)' ) '  to see if they lie in the unit quarter circle.'
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a,i3,a)' ) '  Process ', id, ' is active.'

  if ( id == 0 ) then
    tolerance = 0.001D+00
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The method continues computing until:'
    write ( *, '(a,g14.6)' ) '  PI is computed to a tolerance of ', tolerance
    write ( *, '(a,i9)' ) '  or the number of points used reaches ', point_max
  end if
!
!  ...and broadcasts it to all other processes.
!
  call MPI_Bcast ( tolerance, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, &
    ierr )
!
!  Now create the WORKER communication group, which excludes process the 
!  SERVER process.
!
!
!  Start by getting the group corresponding to the world communicator.
!
  call MPI_Comm_group ( MPI_COMM_WORLD, world_group, ierr )
!
!  Put SERVER on the list of processes to exclude, and create the new 
!  worker group.
!
  server = num_procs - 1

  if ( id == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MONTE_CARLO_MPI - Master process:'
    write ( *, '(a,i8)' ) '  The random number server process ID is ', server
  end if

  ranks(1) = server

  call MPI_Group_excl ( world_group, 1, ranks, worker_group, ierr )
!
!  Use the worker group to create the new worker communicator.
!
  call MPI_Comm_create ( MPI_COMM_WORLD, worker_group, worker_comm, ierr )
!
!  Since we only needed the worker group to create the worker communicator,
!  we can free the worker group now.
!
  call MPI_Group_free ( worker_group, ierr )
!
!  Here is where the computation is carried out.
!
!  The SERVER process waits to receive a request from any other process.
!
  if ( id == server ) then

    do

      tag = mpi_any_tag

      call MPI_Recv ( dummy, 1, MPI_INTEGER, MPI_ANY_SOURCE, tag, &
        MPI_COMM_WORLD, status, ierr )

      tag = status(mpi_tag)

      if ( tag == tag_exit ) then
        exit
      end if

      requestor = status(mpi_source)

      call random_number ( harvest = rands(1:chunk_size) )

      tag = tag_random_sent

      call MPI_Send ( rands, chunk_size, MPI_DOUBLE_PRECISION, requestor, &
        tag, MPI_COMM_WORLD, ierr )

    end do
!
!  Each worker process sends requests for numbers to the random number server.
!
  else

    in_local = 0
    out_local = 0
    point_local = 0

    do

      tag = tag_random_send

      call MPI_Send ( dummy, 1, MPI_INTEGER, server, tag, MPI_COMM_WORLD, ierr )

      call MPI_Recv ( rands, chunk_size, MPI_DOUBLE_PRECISION, server, &
        MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr )

      do i = 1, chunk_size, 2

        x = rands(i)
        y = rands(i+1)

        point_local = point_local + 1

        if ( x**2 + y**2 <= 1.0D+00 ) then
          in_local = in_local + 1
        else
          out_local = out_local + 1
        end if

      end do

      call MPI_Reduce ( in_local, in_total, 1, MPI_INTEGER, mpi_sum, &
        0, worker_comm, ierr )

      call MPI_Reduce ( out_local, out_total, 1, MPI_INTEGER, mpi_sum, &
        0, worker_comm, ierr )

      call MPI_Reduce ( point_local, point_total, 1, MPI_INTEGER, mpi_sum, &
        0, worker_comm, ierr )
!
!  The Master process now checks the value of PI, and the size of POINT_TOTAL.
!
      if ( id == 0 ) then

        pi_est = 4.0D+00 * in_total / point_total

        error = abs ( pi_est - pi_true )

        done = 1
!
!  If it's time to stop, the Master process informs the random number server.
!
        if ( error < tolerance .or. point_max <= point_total ) then

          done = 0
          tag = tag_exit

          call MPI_Send ( dummy, 1, MPI_INTEGER, server, tag, MPI_COMM_WORLD, &
            ierr )

        end if

      end if
!
!  The Master process broadcasts to all workers whether to quit or go on.
!
      call MPI_Bcast ( done, 1, MPI_INTEGER, 0, worker_comm, ierr )

      if ( done == 0 ) then
        exit
      end if

    end do

  end if

  if ( id == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MONTE_CARLO_MPI - Master process:'
    write ( *, '(a,g24.16)' ) '  Estimate for PI = ', pi_est
    write ( *, '(a,g14.6)' ) '  Error =           ', error
    write ( *, '(a,i8)' ) '  Number of points = ', point_total
  end if
!
!  We free the worker communicator.
!
  call MPI_Comm_free ( worker_comm, ierr )
!
!  Terminate.
!
  call MPI_Finalize ( ierr )
!
!  Terminate.
!
  if ( id == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MONTE_CARLO_MPI - Master process:'
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
