program main

!*****************************************************************************80
!
!! MAIN is the main program for MATMAT.
!
!  Discussion:
!
!    MATMAT uses MPI to multiply two matrices.
!
!    The computation C = A * B is carried out by giving every process
!    a copy of the matrix A.  Then each process is given one column of
!    the matrix B, say "B(*,J)", and computes the product of A times
!    this column, returning the result, which is column J of C.
!
!    This code is "self scheduling", because as soon as a process finishes
!    working on one column, it goes back to the master process and requests
!    another one.  This is an appropriate way to schedule work when the
!    processes are likely to run at different speeds, or individual
!    tasks might take greatly varying times.  Neither problem is likely
!    to occur here, but this is meant as an illustration.
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
!    Snir, Otto, Huss-Lederman, Walker, Dongarra,
!    MPI - The Complete Reference,
!    Volume 1, The MPI Core,
!    second edition,
!    MIT Press, 1998.
!
  use mpi
!
!  NA and MB must be equal.
!
  integer, parameter :: ma = 10
  integer, parameter :: na = 20
  integer, parameter :: mb = 20
  integer, parameter :: nb = 25
!
!  All processes need to store a full copy of A.
!
  real a(ma,na)
!
!  Process 0 needs a full copy of B.
!  The other processes do not.
!  We should really do an allocatable array here.
!
  real b(mb,nb)
  real b_column(mb)
!
!  Process 0 needs a full copy of C.
!  The other processes do not.
!  We should really do an allocatable array here.
!
  real c(ma,nb)
  real c_column(1:ma)
  integer dest
  integer i
  integer id
  integer ierr
  integer j
  integer jhi
  integer num_procs
  integer num_received
  integer num_sent
  integer source
  integer status(mpi_status_size)
  integer tag
  real value
!
!  Initialize MPI.
!
  call MPI_Init ( ierr )
!
!  Get this process's ID.
!
  call MPI_Comm_rank ( mpi_comm_world, id, ierr )
!
!  Get the number of processes.
!
  call MPI_Comm_size (  MPI_COMM_WORLD, num_procs, ierr )

  if ( id == 0 ) then
    call timestamp ( )
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MATMAT'
    write ( *, '(a)' ) '  FORTRAN90 version'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  An MPI example program'
    write ( *, '(a)' ) '  to compute a matrix product C = A * B.'
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  The number of processes is ', num_procs
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Number of rows    of matrix A = ', ma
    write ( *, '(a,i8)' ) '  Number of columns of matrix A = ', na
    write ( *, '(a,i8)' ) '  Number of rows    of matrix B = ', mb
    write ( *, '(a,i8)' ) '  Number of columns of matrix B = ', nb

    if ( num_procs < 2 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MATMAT - Fatal error!'
      write ( *, '(a)' ) '  Must have at least 2 processes!'
    end if

  end if

  if ( num_procs < 2 ) then
    call MPI_Abort (  MPI_COMM_WORLD, 1, ierr )
    stop
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a)' ) '  Process ', id, ' is active.'
!
!  The master process initializes A and B.
!
  if ( id == 0 ) then

    do i = 1, ma
      do j = 1, na
        if ( j == i-1 ) then
          a(i,j) = -1.0E+00
        else if ( j == i ) then
          a(i,j) = 2.0E+00
        else if ( j == i+1 ) then
          a(i,j) = -1.0E+00
        else
          a(i,j) = 0.0E+00
        end if
      end do
    end do

    do i = 1, mb
      do j = 1, nb
        if ( i <= j ) then
          b(i,j) = real ( i * ( nb + 1 - j ) ) / real ( nb + 1 )
        else
          b(i,j) = real ( j * ( mb + 1 - i ) ) / real ( mb + 1 )
        end if
      end do
    end do

  end if
!
!  We have rigged the game so that, if the matrices are square, then B is the 
!  inverse of A, and therefore the product C should be the identity matrix.
!
!  Just so we can do rectangular problems too, and still know we're on the
!  right track, compute the exact result here and print out a chunk of it.
!
  if ( id == 0 ) then

    c(1:ma,1:nb) = matmul ( a(1:ma,1:na), b(1:mb,1:nb) )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MATMAT - Master process:'
    write ( *, '(a)' ) '  Initial 5 x 5 block of exact product matrix C:'
    write ( *, '(a)' ) ' '

    do i = 1, 5
      write ( *, '(5g14.6)' ) c(i,1:5)
    end do

  end if
!
!  The master process broadcasts a copy of the entire A matrix 
!  to all other processes.
!
  call MPI_Bcast ( a, ma*nb, mpi_real, 0, MPI_COMM_WORLD, ierr )
!
!  Now the master process distributes columns of B to the other processes,
!  and collects the products C(1:MA,J) = A(1:MA,1:NA) * B(1:NA,J)
!
  if ( id == 0 ) then

    num_sent = 0
!
!  Start by sending column J of B to process J.
!
    jhi = min ( num_procs-1, nb )

    do j = 1, jhi

      b_column(1:mb) = b(1:mb,j)

      dest = j
      tag = j

      call MPI_Send ( b_column, mb, MPI_REAL, dest, tag, MPI_COMM_WORLD, ierr )

      num_sent = num_sent + 1

    end do
!
!  Process 0 waits to receive a result from any process, and sends that
!  process a new column, or a terminate signal if there are no more columns.
!
!  Once all columns have been received, process 0 exits from the loop.
!
    num_received = 0

    do while ( num_received < nb )

      call MPI_Recv ( c_column, ma, MPI_REAL, MPI_ANY_SOURCE, &
        MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr )

      num_received = num_received + 1
      source = status(mpi_source)
      tag = status(mpi_tag)

      c(1:ma,tag) = c_column(1:ma)

      if ( num_sent < nb ) then

        num_sent = num_sent + 1
        b_column(1:mb) = b(1:mb,num_sent)
        dest = source
        tag = num_sent

        call MPI_Send ( b_column, mb, MPI_REAL, dest, tag, MPI_COMM_WORLD, &
          ierr )

      else

        value = 1.0E+00
        dest = source
        tag = 0

        call MPI_Send ( value, 1,  MPI_REAL, dest, tag, MPI_COMM_WORLD, ierr )

      end if

    end do
!
!  Each worker process receives a column of B.
!
  else if ( id <= nb ) then

    do

      call MPI_Recv ( b_column, mb, MPI_REAL, 0, MPI_ANY_TAG, &
        MPI_COMM_WORLD, status, ierr )
!
!  If the tag is 0, we're actually being told to terminate.
!
      tag = status(mpi_tag)

      if ( tag == 0 ) then
        exit
      end if

      c_column(1:ma) = matmul ( a(1:ma,1:na), b_column(1:mb) )

      call MPI_Send ( c_column, ma, MPI_REAL, 0, tag, MPI_COMM_WORLD, ierr )

    end do

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Process ', id
    write ( *, '(a)' ) '  MPI has no work for me!'

  end if
! 
!  Process 0 prints out a bit of the answer.
!
  if ( id == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MATMAT - Master process:'
    write ( *, '(a)' ) '  Initial 5 x 5 block of computed product matrix C:'
    write ( *, '(a)' ) ' '

    do i = 1, 5
      write ( *, '(5g14.6)' ) c(i,1:5)
    end do

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
    write ( *, '(a)' ) 'MATMAT:'
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
