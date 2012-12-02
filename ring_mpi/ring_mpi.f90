program main

!*****************************************************************************80
!
!! MAIN is the main program for RING_MPI.
!
!  Discussion:
!
!    RING_MPI sends messages of various size from process 0 to 1 to 2 to
!    ...the last process and then back to 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Peter Pacheco,
!    Parallel Programming with MPI,
!    Morgan Kaufman, 1996,
!    ISBN: 1558603395,
!    LC: QA76.642.P3.
!
  use mpi

  implicit none

  integer ( kind = 4 ) error
  integer ( kind = 4 ) id
  integer ( kind = 4 ) p
!
!  Initialize MPI.
!
  call MPI_Init ( error )
!
!  Get the number of processes.
!
  call MPI_Comm_size ( MPI_COMM_WORLD, p, error )
!
!  Get the individual process ID.
!
  call MPI_Comm_rank ( MPI_COMM_WORLD, id, error )
!
!  Print a message.
!
  if ( id == 0 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RING_MPI:'
    write ( *, '(a)' ) '  FORTRAN90/MPI version'
    write ( *, '(a)' ) '  Measure time required to transmit data around'
    write ( *, '(a)' ) '  a ring of processes'
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  The number of processes is ', p

  end if

  call ring_io ( p, id )
!
!  Shut down MPI.
!
  call MPI_Finalize ( error )
!
!  Terminate.
!
  if ( id == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RING_MPI:'
    write ( *, '(a)' ) '  Normal end of execution.'
  end if

  stop
end
subroutine ring_io ( p, id )

!*****************************************************************************80
!
!! RING_IO carries out the tasks of process ID, of a total of P processes.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Peter Pacheco,
!    Parallel Programming with MPI,
!    Morgan Kaufman, 1996,
!    ISBN: 1558603395,
!    LC: QA76.642.P3.
!
  use mpi

  implicit none

  integer ( kind = 4 ), parameter :: n_test_num = 5

  integer ( kind = 4 ) dest
  integer ( kind = 4 ) error
  integer ( kind = 4 ) i
  integer ( kind = 4 ) id
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  integer ( kind = 4 ), dimension ( n_test_num ) :: n_test = (/ &
    100, 1000, 10000, 100000, 1000000 /)
  integer ( kind = 4 ) p
  integer ( kind = 4 ) source
  integer ( kind = 4 ) status(MPI_STATUS_SIZE)
  real ( kind = 8 ) tave
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 10
  real ( kind = 8 ) tmax
  real ( kind = 8 ) tmin
  real ( kind = 8 ) wtime
  real ( kind = 8 ), allocatable :: x(:)

  if ( id == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i6,a)' ) '  Timings based on ', test_num, ' experiments'
    write ( *, '(a,i6,a)' ) '  N double precision values were sent'
    write ( *, '(a)' ) '  in a ring transmission starting and ending at process 0'
    write ( *, '(a,i6,a)' ) '  and using a total of ', p, ' processes.'  
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '         N           T min           T ave           T max'
    write ( *, '(a)' ) ' '
  end if
!
!  Choose message size.
!
  do i = 1, n_test_num
    
    n = n_test(i)

    allocate ( x(1:n) )
!
!  Process 0 sends very first message, 
!  then waits to receive the "echo" that has gone around the world.
!
    if ( id == 0 ) then

      dest = 1
      source = p - 1

      tave = 0.0D+00
      tmin = huge ( 1.0D+00 )
      tmax = 0.0D+00

      do test = 1, test_num
!
!  Just in case, set the entries of X in a way that identifies
!  which iteration of the test is being carried out.
!
        do j = 1, n
          x(j) = real ( test + j - 1, kind = 8 )
        end do

        wtime = MPI_Wtime ( )
        call MPI_Send ( x, n, MPI_DOUBLE_PRECISION, dest,   0, MPI_COMM_WORLD,         error )
        call MPI_Recv ( x, n, MPI_DOUBLE_PRECISION, source, 0, MPI_COMM_WORLD, status, error )
        wtime = MPI_Wtime ( ) - wtime
!
!  Record the time it took.
!
        tave =       tave + wtime
        tmin = min ( tmin,  wtime )
        tmax = max ( tmax,  wtime )

      end do

      tave = tave / real ( test_num, kind = 8 )

      write ( *, '(2x,i8,2x,g14.6,2x,g14.6,2x,g14.6)' ) n, tmin, tave, tmax
!
!  Worker ID must receive first from ID-1, then send to ID+1.
!
    else

      source = id - 1
      dest = mod ( id + 1, p )
 
      do test = 1, test_num
        call MPI_Recv ( x, n, MPI_DOUBLE_PRECISION, source, 0, MPI_COMM_WORLD, status, error )
        call MPI_Send ( x, n, MPI_DOUBLE_PRECISION, dest,   0, MPI_COMM_WORLD,         error )
      end do

    end if

    deallocate ( x )

  end do

  return
end
