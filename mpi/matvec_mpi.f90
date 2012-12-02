program main

!*****************************************************************************80
!
!! MAIN is the main program for MATVEC_MPI.
!
!  Discussion:
!
!    MATVEC uses MPI to compute a matrix-vector product b = A * x.
!
!    This is the simple self-scheduling version.  Each worker is given a copy 
!    of x, and then is fed one row of A.  As soon as it computes 
!    B(I) = A(I,1:N)*x(1:N), it is given another column of A, unless there are
!    no more, in which case it is sent a "terminate" message.  Thus, a faster
!    process will be given more work to do.
!
!    By using allocatable arrays, the amount of memory used has been controlled.
!    The master process allocates A and x, but the worker processes only
!    allocate enough memory for one row of A, and x.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 September 2009
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

  real ( kind = 8 ), allocatable :: a(:,:)
  real ( kind = 8 ), allocatable :: a_row(:)
  real ( kind = 8 ) ans
  real ( kind = 8 ), allocatable :: b(:)
  integer ( kind = 4 ) dest
  integer ( kind = 4 ) i
  integer ( kind = 4 ) id
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j_one
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) num_procs
  integer ( kind = 4 ) num_rows
  integer ( kind = 4 ) num_workers
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ) status(MPI_STATUS_SIZE)
  integer ( kind = 4 ) tag
  real ( kind = 8 ), allocatable :: x(:)
!
!  Initialize MPI.
!
  call MPI_Init ( ierr )
!
!  Get this processor's ID.
!
  call MPI_Comm_rank ( MPI_COMM_WORLD, id, ierr )
!
!  Get the number of processors.
!
  call MPI_Comm_size ( MPI_COMM_WORLD, num_procs, ierr )

  if ( id == 0 ) then
    call timestamp ( )
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MATVEC_MPI:'
    write ( *, '(a)' ) '  FORTRAN90 version'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  An MPI example program to compute'
    write ( *, '(a)' ) '  a matrix-vector product b = A * x.'
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  The number of processes is ', num_procs
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a)' ) '  Process ', id, ' is active.'

  m = 100
  n = 50

  if ( id == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  The number of rows is    ', m
    write ( *, '(a,i8)' ) '  The number of columns is ', n
  end if
!
!  The master process allocates and initializes A and X.
!
  if ( id == 0 ) then

    allocate ( a(1:m,1:n) )
    allocate ( x(1:n) )
    allocate ( b(1:m) )

    do i = 1, m
      do j = 1, n
        a(i,j) = sqrt ( 2.0D+00 / real ( n + 1, kind = 8 ) ) &
          * sin ( real ( i * j, kind = 8 ) * pi / real ( n + 1, kind = 8 ) )
      end do
    end do
!
!  X is specially chosen so that b = A * x is known in advance.
!  The value of B will be zero, except that entry J_ONE will be 1.
!  Pick any value of J_ONE between 1 and M.
!
    j_one = 17

    do i = 1, n
      x(i) = sqrt ( 2.0D+00 / real ( n + 1, kind = 8 ) ) &
        * sin ( real ( i * j_one, kind = 8 ) * pi / real ( n + 1, kind = 8 ) )
    end do
!
!  Worker processes set aside room for one row of A, and for the vector X.
!
  else

    allocate ( a_row(1:n) )
    allocate ( x(1:n) )

  end if
!
!  Process 0 broadcasts the vector X to the other processes.
!
  if ( id == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Broadcasting vector X to all processes.'
  end if

  call MPI_Bcast ( x, n, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
!
!  Process 0 sends one row of A to all the other processes.
!
!  Note that the call to MPI_Send uses a FORTRAN90 array section.  Even
!  though the elements of a 2D array row are not contiguous as stored in memory,
!  FORTRAN90 interprets the expression "A(I,1:N)" as requiring it to make
!  a temporary, and contiguous, copy of the indicated elements.
!
  if ( id == 0 ) then

    num_rows = 0

    do i = 1, num_procs - 1

      num_rows = num_rows + 1
      dest = i
      tag = num_rows

      call MPI_Send ( a(num_rows,1:n), n, MPI_DOUBLE_PRECISION, dest, tag, &
        MPI_COMM_WORLD, ierr )

    end do
     
    num_workers = num_procs - 1

    do

      call MPI_Recv ( ans, 1, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, &
        MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr )

      tag = status(MPI_TAG)
      b(tag) = ans

      if ( num_rows < m ) then

        num_rows = num_rows + 1
        dest = status(MPI_SOURCE)
        tag = num_rows

        call MPI_Send ( a(num_rows,1:n), n, MPI_DOUBLE_PRECISION, dest, tag, &
          MPI_COMM_WORLD, ierr )
!
!  Even though we are sending a message for which the TAG is important, not
!  the data, we need to include the right type and amount of data for the
!  message to be accepted.
!
      else

        num_workers = num_workers - 1
        dest = status(MPI_SOURCE)
        tag = m + 1

        call MPI_Send ( a(1,1:n), n, MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr )

        if ( num_workers == 0 ) then
          exit
        end if

      end if

    end do

    deallocate ( a )
    deallocate ( x )
!
!  Each worker process repeatedly receives rows of A (with TAG indicating 
!  which row it is), computes dot products A(I,1:N) * X(1:N) and returns 
!  the result (and TAG), until receiving the "DONE" message.
!
  else

    do

      call MPI_Recv ( a_row, n, MPI_DOUBLE_PRECISION, 0, &
        MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr )

      tag = status(MPI_TAG)

      if ( tag == m + 1 ) then
        write ( *, '(a,i8,a)' ) '  Process ', id, ' shutting down.'
        exit
      end if 

      ans = dot_product ( a_row(1:n), x(1:n) )

      call MPI_Send ( ans, 1, MPI_DOUBLE_PRECISION, 0, tag, &
        MPI_COMM_WORLD, ierr )

    end do

    deallocate ( a_row )
    deallocate ( x )

  end if
! 
!  Print out the answer.
!
  if ( id == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MATVEC_MPI - Master process:'
    write ( *, '(a)' ) '  Product vector b = A * x'
    write ( *, '(a,i8)' ) '  (Should be zero, except for a 1 in entry ', j_one
    write ( *, '(a)' ) ' '
    do i = 1, m
      write ( *, '(i8,g14.6)' )  i, b(i)
    end do

    deallocate ( b )

  end if
!
!  Terminate MPI.
!
  call MPI_Finalize ( ierr )
!
!  Terminate
!
  if ( id == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MATVEC_MPI:'
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

