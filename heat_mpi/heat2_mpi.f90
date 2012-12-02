program main

!*****************************************************************************80
!
!!  MAIN is the main program for HEAT_MPI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 April 2008
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
!    ISBN: 0262571323,
!    LC: QA76.642.G76.
!
!    Marc Snir, Steve Otto, Steven Huss-Lederman, David Walker, 
!    Jack Dongarra,
!    MPI: The Complete Reference,
!    Volume I: The MPI Core,
!    Second Edition,
!    MIT Press, 1998,
!    ISBN: 0-262-69216-3,
!     LC: QA76.642.M65.
!
  use mpi

  integer id
  integer ierr
  integer p
  double precision wtime


  call MPI_Init ( ierr )

  call MPI_Comm_rank ( MPI_COMM_WORLD, id, ierr )

  call MPI_Comm_size ( MPI_COMM_WORLD, p, ierr )

  if ( id == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HEAT_MPI:'
    write ( *, '(a)' ) '  FORTRAN90/MPI version.'
    write ( *, '(a)' ) '  Solve the 1D time-dependent heat equation.'
  end if
!
!  Record the starting time.
!
  if ( id == 0 ) then
    wtime = MPI_Wtime ( )
  end if

  call update ( id, p )
!
!  Record the final time.
!
  if ( id == 0 ) then
    wtime = MPI_Wtime ( ) - wtime
    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  Wall clock elapsed seconds = ', wtime
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
    write ( *, '(a)' ) 'HEAT_MPI:'
    write ( *, '(a)' ) '  Normal end of execution.'
  end if

  stop
end
subroutine update ( id, p )

!*****************************************************************************80
!
!! UPDATE computes the solution of the heat equation.
!
!  Discussion:
!
!    If there is only one processor (P .eq. 1 ), then the program writes the
!    values of X and H to files.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 April 2008
!
!  Author:
!
!    John Burkardtc
!
!  Parameters:
!
!    Input, integer ID, the id of this processor.
!
!    Input, integer P, the number of processors.
!
  use mpi

  integer, parameter :: n = 11

  double precision boundary_condition
  double precision cfl
  double precision h(0:n+1)
  integer h_file
  double precision h_new(0:n+1)
  integer i
  integer id
  double precision initial_condition

  integer j
  integer j_max
  integer j_min
  double precision k
  integer p
  double precision rhs
  integer status(MPI_STATUS_SIZE)
  integer tag
  double precision time
  double precision time_delta
  double precision time_max
  double precision time_min
  double precision time_new
  double precision x(0:n+1)
  double precision x_delta
  integer x_file
  double precision x_max
  double precision x_min

  h_file = 11
  j_max = 100
  j_min = 0
  k = 0.002
  x_file = 12
  time_max = 10.0
  time_min = 0.0
  x_max = 1.0
  x_min = 0.0
!
!  Have process 0 print out some information.
!
  if ( id == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Compute an approximate solution to the time dependent'
    write ( *, '(a)' ) '  one dimensional heat equation:'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '    dH/dt - K * d2H/dx2 = f(x,t)'
    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6,a,g14.6)' ) &
    '  for ', x_min, ' = x_min < x < x_max = ', x_max
    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6,a,g14.6)' ) &
    '  and ', time_min, ' = time_min < t <= t_max = ', time_max
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Boundary conditions are specified at x_min and x_max.'
    write ( *, '(a)' ) '  Initial conditions are specified at time_min.'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The finite difference method is used to discretize the'
    write ( *, '(a)' ) '  differential equation.'
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8,a)' ) '  This uses ', p * n, ' equally spaced points in X'
    write ( *, '(a,i8,a)' ) '  and ', j_max, ' equally spaced points in time.'
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8,a)' ) '  Parallel execution is done using ', p, ' processors.'
    write ( *, '(a)' ) '  Domain decomposition is used.'
    write ( *, '(a,i8,a)' ) '  Each processor works on ', n, ' nodes,'
    write ( *, '(a)' ) '  and shares some information with its immediate neighbors.'
  end if
!
!  Set the X coordinates of the N nodes.
!  We don't actually need ghost values of X but we'll throw them in
!  as X(0) and X(N+1).
!
  do i = 0, n + 1
    x(i) = ( dble (         id * n + i - 1 ) * x_max   &
           + dble ( p * n - id * n - i     ) * x_min ) &
           / dble ( p * n              - 1 )
  end do
!
!  In single processor mode, write out the X coordinates for display.
!
  if ( p == 1 ) then

    open ( unit = x_file, file = 'x_data.txt', status = 'unknown' )

    write ( x_file, '(11f14.6)' ) x(1:n)

    close ( unit = x_file )

  end if
!
!  Set the values of H at the initial time.
!
  time = time_min

  h(0) = 0.0
  do i = 1, n
    h(i) = initial_condition ( x(i), time )
  end do
  h(n+1) = 0.0
  
  time_delta = ( time_max - time_min ) / dble ( j_max - j_min )
  x_delta = ( x_max - x_min ) / dble ( p * n - 1 )
!
!  Check the CFL condition, have processor 0 print out its value,
!  and quit if it is too large.
!
  cfl = k * time_delta / x_delta / x_delta

  if ( id == 0 ) then 
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'UPDATE'
    write ( *, '(a,g14.6)' ) '  CFL stability criterion value = ', cfl 
  end if

  if ( 0.5 <= cfl ) then
    if ( id == 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'UPDATE - Warning!'
      write ( *, '(a)' ) '  Computation cancelled!'
      write ( *, '(a)' ) '  CFL condition failed.'
      write ( *, '(a,g14.6)' ) '  0.5 <= K * dT / dX / dX = ', cfl
    end if
    return
  end if
!
!  In single processor mode, write out the values of H.
!
  if ( p == 1 ) then

    open ( unit = h_file, file = 'h_data.txt', status = 'unknown' )

    write ( h_file, '(11f14.6)' ) h(1:n)

  end if
!
!  Compute the values of H at the next time, based on current data.
!
  do j = 1, j_max

    time_new = ( dble (         j - j_min ) * time_max   &
               + dble ( j_max - j         ) * time_min ) &
               / dble ( j_max     - j_min )
!
!  NOTE THAT THE FOLLOWING SECTION OF CODE IS THE ONLY CHANGE MADE!
!
!  THE PREVIOUS VERSION OF THE SEND'S AND RECEIVE'S WAS WRITTEN IN SUCH A WAY
!  THAT ALL THE PROCESSES BUT ONE SENT A MESSAGE.  ONLY PROCESS 0 WAS
!  ABLE TO RECEIVE.  ONCE IT ACKNOWLEDGED RECEIPT, PROCESS 1 COULD RECEIVE ITS
!  MESSAGE, FREEING PROCESS 2 AND SO ON.  SO FOR A CERTAIN TIME, PROCESS P-1
!  HAD TO WAIT FOR THE CASCADE OF MESSAGES TO BE SENT AND ACKNOWLEDGED.
!
!  THIS VERSION OF THE CODE SHOULD AVOID THAT PROBLEM.
!
!  The odd processors send H(1) to ID-1 and the even ones receive it as H(N+1),
!  then
!  the even processors send H(1) to ID-1 and the odd ones receive it as H(N+1).
!
!  Processor 0 does NOT send.
!  Processor P-1 does NOT receive.
!
    tag = 1

    if ( 0 < id .and. mod ( id, 2 ) == 1 ) then
      call MPI_Send ( h(1), 1, MPI_DOUBLE_PRECISION, id-1, tag, &
        MPI_COMM_WORLD, ierr )
    else if ( id < p - 1 .and. mod ( id, 2 ) == 0 ) then
      call MPI_Recv ( h(n+1), 1,  MPI_DOUBLE_PRECISION, id+1, tag, &
        MPI_COMM_WORLD, status, ierr )
    end if

    if ( 0 < id .and. mod ( id, 2 ) == 0 ) then
      call MPI_Send ( h(1), 1, MPI_DOUBLE_PRECISION, id-1, tag, &
        MPI_COMM_WORLD, ierr )
    else if ( id < p - 1 .and. mod ( id, 2 ) == 1 ) then
      call MPI_Recv ( h(n+1), 1,  MPI_DOUBLE_PRECISION, id+1, tag, &
        MPI_COMM_WORLD, status, ierr )
    end if
!
!  The odd processors send H(N) to ID+1 and the even ones receive it as H(0),
!  then
!  the even processors send H(N) to ID+1 and the odd ones receive it as H(0),
!
!  Processor P-1 does NOT send.
!  Processor 0 does NOT receive.
!
    tag = 2

    if ( id < p - 1 .and. mod ( id, 2 ) == 1 ) then
      call MPI_Send ( h(n), 1, MPI_DOUBLE_PRECISION, id+1, tag, &
        MPI_COMM_WORLD, ierr )
    else if ( 0 < id .and. mod ( id, 2 ) == 0 ) then
      call MPI_Recv ( h(0), 1, MPI_DOUBLE_PRECISION, id-1, tag, &
        MPI_COMM_WORLD, status, ierr )
    end if
!
!  Update the temperature based on the four point stencil.
!
    do i = 1, n
      h_new(i) = h(i) &
        + ( time_delta * k / x_delta / x_delta ) &
        * ( h(i-1) - 2.0 * h(i) + h(i+1) ) &
        + time_delta * rhs ( x(i), time )
    end do
!
!  H at the extreme left and right boundaries was incorrectly computed
!  using the differential equation.  Replace that calculation by
!  the boundary conditions.
!
    if ( 0 == id ) then
      h_new(1) = boundary_condition ( x(1), time_new )
    end if

    if ( id == p - 1 ) then
      h_new(n) = boundary_condition ( x(n), time_new )
    end if
!
!  Update time and temperature.
!
    time = time_new

    do i = 1, n
      h(i) = h_new(i)
    end do
!
!  In single processor mode, add current solution data to output file.
!
    if ( p == 1 ) then
      write ( h_file, '(11f14.6)' ) h(1:n)
    end if

  end do

  if ( p == 1 ) then
    close ( unit = h_file )
  end if

  return
end
function boundary_condition ( x, time )

!*****************************************************************************80
!
!! BOUNDARY_CONDITION evaluates the boundary condition of the differential equation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 April 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision X, TIME, the position and time.
!
!    Output, double precision BOUNDARY_CONDITION, the value of the boundary condition.
!
  implicit none

  double precision boundary_condition
  double precision time
  double precision x
!
!  Left condition:
!
  if ( x < 0.5 ) then
    boundary_condition = 100.0 + 10.0 * sin ( time )
  else
    boundary_condition = 75.0
  end if

  return
end
function initial_condition ( x, time )

!*****************************************************************************80
!
!! INITIAL_CONDITION evaluates the initial condition of the differential equation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 April 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision X, TIME, the position and time.
!
!    Output, double precision INITIAL_CONDITION, the value of the initial condition.
!
  implicit none

  double precision initial_condition
  double precision time
  double precision x

  initial_condition = 95.0

  return
end
function rhs ( x, time )

!*****************************************************************************80
!
!! RHS evaluates the right hand side of the differential equation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 April 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision X, TIME, the position and time.
!
!    Output, double precision RHS, the value of the right hand side.
!
  double precision rhs
  double precision time
  double precision x

  rhs = 0.0

  return
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
