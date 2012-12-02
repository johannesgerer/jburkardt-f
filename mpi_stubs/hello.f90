program main

!*****************************************************************************80
!
!! MAIN is the main program for HELLO.
!
!  Discussion:
!
!    HELLO is a simple MPI test program.
!
!    Each process prints out a "Hello, world!" message.
!
!    The master process also prints out a short message.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 May 2008
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
  implicit none

  include 'mpif.h'

  integer error
  integer, parameter :: master = 0
  integer num_procs
  integer world_id
!
!  Initialize MPI.
!
  call MPI_Init ( error )
!
!  Get the number of processes.
!
  call MPI_Comm_size ( MPI_COMM_WORLD, num_procs, error )
!
!  Get the individual process ID.
!
  call MPI_Comm_rank ( MPI_COMM_WORLD, world_id, error )
!
!  Print a message.
!
  if ( world_id == master ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HELLO_WORLD - Master process:'
    write ( *, '(a)' ) '  FORTRAN90 version'
    write ( *, '(a)' ) '  An MPI test program.'
    write ( *, '(a)' ) ' '
    write ( *, '(a,i3)' ) '  The number of processes is ', num_procs
    write ( *, '(a)' ) ' '
  end if

  write ( *, '(a,i3,a)' ) '  Process ', world_id, ' says "Hello, world!"'
!
!  Terminate MPI.
!
  call MPI_Finalize ( error )
!
!  Terminate.
!
  if ( world_id == master ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HELLO_WORLD:'
    write ( *, '(a)' ) '  Normal end of execution.'
  end if

  stop
end
