program main

!*****************************************************************************80
!
!! QUAD_MAIN is the main program for the quadrature example.
!
!  Discussion:
!
!    This example estimates an integral using quadrature.
!
!    The integral of F(X) = 4 / ( 1 + X * X ) from 0 to 1 is PI.
!
!    We break up the interval [0,1] into N subintervals, evaluate
!    F(X) at the midpoint of each subinterval, and multiply the
!    sum of these values by N to get an estimate for the integral.
!
!    If we have M processes available because we are using MPI, then
!    we can ask processes 0, 1, 2, ... M-1 to handle the subintervals
!    in the following order:
!
!          0      1       2            M-1  <-- Process numbers begin at 0
!     ------ ------  ------  -----  ------
!          1      2       3    ...       M
!        M+1    M+2     M+3    ...     2*M
!      2*M+1    2*M+2 2*M+3    ...     3*M
!                              
!    and so on up to subinterval N.  The partial sums collected by 
!    each process are then sent to the master process to be added 
!    together to get the estimated integral.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 September 2011
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
  use, intrinsic :: iso_c_binding

  use mpi

  interface
    function quad_compute ( n ) bind ( c )
      use iso_c_binding
      integer ( c_int ), VALUE :: n
      real ( c_double ) :: quad_compute
    end function quad_compute
  end interface

  real ( kind = 8 ) f
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) id
  integer ( kind = 4 ) ierr
  real ( kind = 8 ) my_part
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_part
  integer ( kind = 4 ) p
  real ( kind = 8 ) q
  real ( kind = 8 ) q_diff
  real ( kind = 8 ), parameter :: q_exact = 3.141592653589793238462643D+00
  real ( kind = 8 ) q_part
  real ( kind = 8 ) sum2
  real ( kind = 8 ) wtime_diff
  real ( kind = 8 ) wtime_end
  real ( kind = 8 ) wtime_start
  real ( kind = 8 ) x

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'QUAD_MAIN:'
  write ( *, '(a)' ) '  FORTRAN90/C++/MPI version'
  write ( *, '(a)' ) '  An example program to estimate an integral.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this example, we have a FORTRAN90 main program,'
  write ( *, '(a)' ) '  and an underlying C++ function, which must cooperate.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Moreover, the main program invokes the MPI library, and'
  write ( *, '(a)' ) '  the C++ function must execute within that environment.'
!
!  Initialize MPI.
!
  call MPI_Init ( ierr )
!
!  Get this process's ID.
!
  call MPI_Comm_rank ( MPI_COMM_WORLD, id, ierr )
!
!  Find out how many processes are available.
!
  call MPI_Comm_size ( MPI_COMM_WORLD, p, ierr )

  if ( id == 0 ) then
    write ( *, '(a,i8)' ) '  The number of processes is ', p
    wtime_start = MPI_Wtime ( )
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a)' ) 'Process ', id, ' is active.'
!
!  Here we call the C++ function.
!
  n = 1000
  q_part = quad_compute ( n )
!
!  Each process sends its value of Q_PART to the master process, to
!  be summed in Q.
!
  call MPI_Reduce ( q_part, q, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, &
    MPI_COMM_WORLD, ierr )

  if ( id == 0 ) then
    wtime_end = MPI_Wtime ( )
  end if
!
!  Finish up.
!
  call MPI_Finalize ( ierr )
!
!  Single process execution.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'QUAD_MAIN:'
  write ( *, '(a,g24.16)' ) '  Integral estimate  ', q
  write ( *, '(a,g24.16)' ) '  Exact value is     ', q_exact
  q_diff = abs ( q_exact - q )
  write ( *, '(a,g24.16)' ) '  Error is           ', q_diff

  wtime_diff = wtime_end - wtime_start

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Elapsed wall clock seconds = ', wtime_diff
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'QUAD_MAIN:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

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
  integer   ( kind = 4 ) d
  integer   ( kind = 4 ) h
  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) mm
  character ( len = 9  ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer   ( kind = 4 ) n
  integer   ( kind = 4 ) s
  integer   ( kind = 4 ) values(8)
  integer   ( kind = 4 ) y

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
