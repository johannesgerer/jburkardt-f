program main

!*****************************************************************************80
!
!! MAIN is the main program for FD1D_PREDATOR_PREY.
!
!  Discussion:
!
!    This program sets up and solves a finite difference 1D predator
!    prey system.
!
!    The nondimensional problem has the form
!
!      du/dt =         del u + ( 1 - u ) * u        - v * h(u/alpha)
!
!      dv/dt = delta * del v     - gamma * v + beta * v * h(u/alpha)
!
!    with initial conditions:
!
!      u(x,0) = u0(x)
!      v(x,0) = v0(x)
!
!    and boundary conditions at the left and right endpoints of [A,B]:
!
!      du/dx = 0
!      dv/dx = 0
!
!    The Type II functional response employed here is
!
!      h(eta) = eta / ( 1 + eta )
!
!    The parameters ALPHA, BETA, GAMMA and DELTA are strictly positive.
!
!    The user must input a value H specifying the desired space step
!    to be used in discretizing the space dimension.
!
!    A finite difference scheme is employed to integrate the problem
!    from time 0 to a maximum time T.  The user must input the value
!    T, as well as an appropriate time step DELT.
!
!  Example:
!
!    A typical input for this problem is:
!
!      ALPHA =   0.3
!      BETA  =   2.0
!      GAMMA =   0.8
!      DELTA =   1.0
!      A     =   0.0
!      B     = 200.0
!      H     =   0.5
!      T     =  40.0
!      DELT  =   0.0104
!
!    with the following initial values of U and V supplied in
!    auxiliary subroutines:
!
!      u0(1:n) = exp ( - ( x(1:n) - 100.0 )**2 ) / 5.0
!      v0(1:n) = 2.0 / 5.0
!
!  Modified:
!
!    21 November 2004
!
!  Author:
!
!    Original MATLAB version by Marcus Garvie,
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Marcus Garvie,
!    Computational Algorithms for Spatially Extended Predator-Prey
!    Systems with a Holling Type II Functional Response,
!    To appear.
!
!  Parameters:
!
!    User input, real ( kind = 8 ) ALPHA, the value of the parameter ALPHA.
!
!    User input, real ( kind = 8 ) BETA, the value of the parameter BETA.
!
!    User input, real ( kind = 8 ) GAMMA, the value of the parameter GAMMA.
!
!    User input, real ( kind = 8 ) DELTA, the value of the parameter DELTA.
!
!    User input, real ( kind = 8 ) A, B, the left and right endpoints
!    of the interval.
!
!    User input, real ( kind = 8 ) H, the "space step", the desired
!    spacing between
!    nodes in [A,B].
!
!    User input, real ( kind = 8 ) T, the final time.  The problem is to
!    be integrated from 0 to T.
!
!    User input, real ( kind = 8 ) DELT, the time step to use.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) alpha
  real ( kind = 8 ) b
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: b1
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: b2
  real ( kind = 8 ) beta
  integer bigj
  real ( kind = 8 ) delt
  real ( kind = 8 ) delta
  real ( kind = 8 ), allocatable, dimension ( : ) :: f
  real ( kind = 8 ), allocatable, dimension ( : ) :: g
  real ( kind = 8 ) gamma
  integer gauss
  real ( kind = 8 ) h
  real ( kind = 8 ), allocatable, dimension ( : ) :: hhat
  integer i
  integer info
  integer, parameter :: it_max = 25
  real ( kind = 8 ) mu
  integer n
  integer nt
  real ( kind = 8 ) t
  integer time_steps
  real ( kind = 8 ), allocatable, dimension ( : ) :: u
  real ( kind = 8 ), allocatable, dimension ( : ) :: u_rhs
  real ( kind = 8 ), allocatable, dimension ( : ) :: v
  real ( kind = 8 ), allocatable, dimension ( : ) :: v_rhs
  real ( kind = 8 ), allocatable, dimension ( : ) :: x

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FD1D_PREDATOR_PREY'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A one dimensional finite difference algorithm'
  write ( *, '(a)' ) '  for a predator-prey system.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Enter parameter alpha:'
  read ( *, * ) alpha
  write ( *, '(a)' ) '  Enter parameter beta:'
  read ( *, * ) beta
  write ( *, '(a)' ) '  Enter parameter gamma:'
  read ( *, * ) gamma
  write ( *, '(a)' ) '  Enter parameter delta:'
  read ( *, * ) delta

  write ( *, '(a)' ) '  Enter a in [a,b]:'
  read ( *, * ) a
  write ( *, '(a)' ) '  Enter b in [a,b]:'
  read ( *, * ) b

  write ( *, '(a)' ) '  Enter space-step h:'
  read ( *, * ) h
  write ( *, '(a)' ) '  Enter maximum time t:'
  read ( *, * ) t
  write ( *, '(a)' ) '  Enter time-step Delta t:'
  read ( *, * ) delt
  write ( *, '(a)' ) '  Enter 0 for direct Gauss solution, '
  write ( *, '(a)' ) '        1 for iterative Jacobi:'
  read ( *, * ) gauss

  if ( gauss /= 0 ) then
    gauss = 1
  end if
!
!  Calculate some constants.
!
  mu = delt / ( h * h )
  bigj = nint ( ( b - a ) / h )
!
!  N = number of degrees of freedom for each dependent variable.
!
  n = bigj + 1

  time_steps = nint ( t / delt )
!
!  Initialization
!
  allocate ( b1(3,1:n) )
  allocate ( b2(3,1:n) )
  allocate ( f(1:n) )
  allocate ( g(1:n) )
  allocate ( hhat(1:n) )
  allocate ( u(1:n) )
  allocate ( u_rhs(1:n) )
  allocate ( v(1:n) )
  allocate ( v_rhs(1:n) )
  allocate ( x(1:n) )
!
!  Set the coordinates of the nodes in [A,B].
!
  do i = 1, n
    x(i) = ( real ( n - i     , kind = 8 ) * a   &
           + real (     i - 1 , kind = 8 ) * b ) &
           / real ( n     - 1 , kind = 8 )

  end do
!
!  Call the user-supplied routines to get initial values for U and V.
!
  call u_init ( n, x, u )
  call v_init ( n, x, v )
!
!  Construct the matrix L, without the 1/h^2 factor.
!  For convenience, store it in B1.
!
  b1(1:3,1:n) = 0.0D+00
  b1(1,2) = -2.0D+00
  b1(3,n-1) = -2.0D+00
  b1(1,3:n) = -1.0D+00
  b1(2,1:n) = +2.0D+00
  b1(3,1:n-2) = -1.0D+00
!
!  Construct the matrices B1 & B2.
!
  b1(1:3,1:n) = mu * b1(1:3,1:n)
  b2(1:3,1:n) = delta * b1(1:3,1:n)

  b1(2,1:n) = b1(2,1:n) + 1.0D+00
  b2(2,1:n) = b2(2,1:n) + 1.0D+00
!
!  If we are using Gauss elimination, then LU factor the matrices B1 and B2.
!
  if ( gauss == 0 ) then

    call d3_np_fa ( n, b1, info )

    if ( info /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FD1D_PREDATOR_PREY - Fatal error!'
      write ( *, '(a)' ) '  Matrix B1 is singular.'
      stop
    end if

    call d3_np_fa ( n, b2, info )

    if ( info /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FD1D_PREDATOR_PREY - Fatal error!'
      write ( *, '(a)' ) '  Matrix B2 is singular.'
      stop
    end if

  end if
!
!  March forward in time.
!
  do nt = 1, time_steps
!
!  Evaluate the modified functional response.
!
    hhat(1:n) = u(1:n) / ( alpha + abs ( u(1:n) ) )
!
!  Update the right-hand-side of the linear system.
!
    f(1:n) = u(1:n) - u(1:n) * abs ( u(1:n) ) - v(1:n) * hhat(1:n)
    g(1:n) = beta * v(1:n) * hhat(1:n) - gamma * v(1:n)
    u_rhs(1:n) = u(1:n) + delt * f(1:n)
    v_rhs(1:n) = v(1:n) + delt * g(1:n)
!
!  Solve the linear systems for U and V.
!
    if ( gauss == 0 ) then

      call d3_np_sl ( n, b1, u_rhs, 0 )
      call d3_np_sl ( n, b2, v_rhs, 0 )

      u(1:n) = u_rhs(1:n)
      v(1:n) = v_rhs(1:n)

    else

      call d3_jac_sl ( n, b1, u_rhs, u, it_max, 0 )

      call d3_jac_sl ( n, b2, v_rhs, v, it_max, 0 )

    end if

  end do
!
!  Plot the solution at time level T=N*delt.
!
  call uv_plot_eps ( 'uv.eps', n, x, u, v, 'U and V at final time' )
!
!  Write the solution out to two files.
!
  open (  unit = 1, file = 'u.txt', status = 'replace' )
  write ( 1, '(a)' ) '#  u.txt'
  write ( 1, '(a)' ) '#'
  do i = 1, n
    write ( 1, '(2g14.6)' ) x(i), u(i)
  end do
  close ( unit = 1 )

  open (  unit = 1, file = 'v.txt', status = 'replace' )
  write ( 1, '(a)' ) '#  v.txt'
  write ( 1, '(a)' ) '#'
  do i = 1, n
    write ( 1, '(2g14.6)' ) x(i), v(i)
  end do
  close ( unit = 1 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FD1D_PREDATOR_PREY'
  write ( *, '(a)' ) '  Two text files were created, containing'
  write ( *, '(a)' ) '  values of (X,U) and (X,V) respectively.'
  write ( *, '(a)' ) '  The files are named "u.txt" and "v.txt".'
!
!  Free memory.
!
  deallocate ( b1 )
  deallocate ( b2 )
  deallocate ( f )
  deallocate ( g )
  deallocate ( hhat )
  deallocate ( u )
  deallocate ( u_rhs )
  deallocate ( v )
  deallocate ( v_rhs )
  deallocate ( x )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FD1D_PREDATOR_PREY:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine d3_jac_sl ( n, a, b, x, it_max, job )

!*****************************************************************************80
!
!! D3_JAC_SL solves a D3 system using Jacobi iteration.
!
!  Discussion:
!
!    The D3 storage format is used for a tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:N), the diagonal in
!    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
!    original matrix is "collapsed" vertically into the array.
!
!    This routine simply applies a given number of steps of the
!    iteration to an input approximate solution.  On first call, you can
!    simply pass in the zero vector as an approximate solution.  If
!    the returned value is not acceptable, you may call again, using
!    it as the starting point for additional iterations.
!
!  Example:
!
!    Here is how a D3 matrix of order 5 would be stored:
!
!       *  A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54  *
!
!  Modified:
!
!    20 September 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be at least 2.
!
!    Input, real ( kind = 8 ) A(3,N), the D3 matrix.
!
!    Input, real ( kind = 8 ) B(N), the right hand side of the linear system.
!
!    Input/output, real ( kind = 8 ) X(N), an approximate solution
!    to the system.
!
!    Input, integer IT_MAX, the maximum number of iterations.
!
!    Input, integer JOB, specifies the system to solve.
!    0, solve A * x = b.
!    nonzero, solve A' * x = b.
!
  implicit none

  integer n

  real ( kind = 8 ) a(3,n)
  real ( kind = 8 ) b(n)
  integer i
  integer it_max
  integer it_num
  integer job
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xnew(n)
!
!  No diagonal matrix entry can be zero.
!
  do i = 1, n
    if ( a(2,i) == 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'D3_JAC_SL - Fatal error!'
      write ( *, '(a,i6)' ) '  Zero diagonal entry, index = ', i
      return
    end if
  end do

  if ( job == 0 ) then

    do it_num = 1, it_max

      xnew(1) =   b(1)                   - a(3,1) * x(2)
      do i = 2, n - 1
        xnew(i) = b(i) - a(1,i) * x(i-1) - a(3,i) * x(i+1)
      end do
      xnew(n) =   b(n) - a(1,n) * x(n-1)

      xnew(1:n) = xnew(1:n) / a(2,1:n)

      x(1:n) = xnew(1:n)

    end do

  else

    do it_num = 1, it_max

      xnew(1) =   b(1)                     - a(1,2) * x(2)
      do i = 2, n - 1
        xnew(i) = b(i) - a(3,i-1) * x(i-1) - a(1,i+1) * x(i+1)
      end do
      xnew(n) =   b(n) - a(3,n-1) * x(n-1)

      xnew(1:n) = xnew(1:n) / a(2,1:n)

      x(1:n) = xnew(1:n)

    end do

  end if

  return
end
subroutine d3_np_fa ( n, a, info )

!*****************************************************************************80
!
!! D3_NP_FA factors a D3 matrix without pivoting.
!
!  Discussion:
!
!    The D3 storage format is used for a tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:N), the diagonal in
!    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
!    original matrix is "collapsed" vertically into the array.
!
!    Because this routine does not use pivoting, it can fail even when
!    the matrix is not singular, and it is liable to make larger
!    errors.
!
!    D3_NP_FA and D3_NP_SL may be preferable to the corresponding
!    LINPACK routine DGTSL for tridiagonal systems, which factors and solves
!    in one step, and does not save the factorization.
!
!  Example:
!
!    Here is how a D3 matrix of order 5 would be stored:
!
!       *  A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54  *
!
!  Modified:
!
!    02 November 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be at least 2.
!
!    Input/output, real ( kind = 8 ) A(3,N).
!    On input, the tridiagonal matrix.  On output, factorization information.
!
!    Output, integer INFO, singularity flag.
!    0, no singularity detected.
!    nonzero, the factorization failed on the INFO-th step.
!
  implicit none

  integer n

  real ( kind = 8 ) a(3,n)
  integer i
  integer info

  info = 0

  do i = 1, n-1

    if ( a(2,i) == 0.0D+00 ) then
      info = i
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'D3_NP_FA - Fatal error!'
      write ( *, '(a,i6)' ) '  Zero pivot on step ', info
      return
    end if
!
!  Store the multiplier in L.
!
    a(3,i) = a(3,i) / a(2,i)
!
!  Modify the diagonal entry in the next column.
!
    a(2,i+1) = a(2,i+1) - a(3,i) * a(1,i+1)

  end do

  if ( a(2,n) == 0.0D+00 ) then
    info = n
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'D3_NP_FA - Fatal error!'
    write ( *, '(a,i6)' ) '  Zero pivot on step ', info
    return
  end if

  return
end
subroutine d3_np_sl ( n, a_lu, b, job )

!*****************************************************************************80
!
!! D3_NP_SL solves a D3 system factored by D3_NP_FA.
!
!  Discussion:
!
!    The D3 storage format is used for a tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:N), the diagonal in
!    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
!    original matrix is "collapsed" vertically into the array.
!
!  Example:
!
!    Here is how a D3 matrix of order 5 would be stored:
!
!       *  A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54  *
!
!  Modified:
!
!    02 November 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be at least 2.
!
!    Input, real ( kind = 8 ) A_LU(3,N), the LU factors from D3_NP_FA.
!
!    Input/output, real ( kind = 8 ) B(N).
!    On input, B contains the right hand side of the linear system.
!    On output, B contains the solution of the linear system.
!
!    Input, integer JOB, specifies the system to solve.
!    0, solve A * x = b.
!    nonzero, solve A' * x = b.
!
  implicit none

  integer n

  real ( kind = 8 ) a_lu(3,n)
  real ( kind = 8 ) b(n)
  integer i
  integer job

  if ( job == 0 ) then
!
!  Solve L * Y = B.
!
    do i = 2, n
      b(i) = b(i) - a_lu(3,i-1) * b(i-1)
    end do
!
!  Solve U * X = Y.
!
    do i = n, 1, -1
      b(i) = b(i) / a_lu(2,i)
      if ( 1 < i ) then
        b(i-1) = b(i-1) - a_lu(1,i) * b(i)
      end if
    end do

  else
!
!  Solve U' * Y = B
!
    do i = 1, n
      b(i) = b(i) / a_lu(2,i)
      if ( i < n ) then
        b(i+1) = b(i+1) - a_lu(1,i+1) * b(i)
      end if
    end do
!
!  Solve L' * X = Y.
!
    do i = n-1, 1, -1
      b(i) = b(i) - a_lu(3,i) * b(i+1)
    end do

  end if

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
subroutine uv_plot_eps ( file_name, n, x, u, v, title )

!*****************************************************************************80
!
!! UV_PLOT_EPS creates an EPS file image of U(X) and V(X).
!
!  Modified:
!
!    01 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the file to create.
!
!    Input, integer N, the number of nodes.
!
!    Input, real ( kind = 8 ) X(N), the coordinates of the nodes.
!
!    Input, real ( kind = 8 ) U(N), V(N), the values of U and V at the nodes.
!
!    Input, character ( len = * ) TITLE, a title for the plot.
!
  implicit none

  integer n

  integer, parameter :: eps_unit = 1
  integer eps_x
  integer eps_y
  character ( len = * ) file_name
  integer i
  integer ios
  real ( kind = 8 ) node_x_max
  real ( kind = 8 ) node_x_min
  real ( kind = 8 ) node_y_max
  real ( kind = 8 ) node_y_min
  character ( len = 40 ) string
  character ( len = * ) title
  real ( kind = 8 ) u(n)
  real ( kind = 8 ) v(n)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x_scale
  real ( kind = 8 ) y_scale
!
!  Determine the range.
!
  node_x_min = minval ( x(1:n) )
  node_x_max = maxval ( x(1:n) )
  node_y_min = min ( minval ( u(1:n) ), minval ( v(1:n) ) )
  node_y_max = max ( maxval ( u(1:n) ), maxval ( v(1:n) ) )

  x_scale = node_x_max - node_x_min
  y_scale = node_y_max - node_y_min

  open ( unit = eps_unit, file = file_name, status = 'replace', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'UV_PLOT_EPS - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output EPS file.'
    stop
  end if

  write ( eps_unit, '(a)' ) '%!PS-Adobe-3.0 EPSF-3.0'
  write ( eps_unit, '(a)' ) &
    '%%Creator: uv_plot_eps(fd1d.f90)'
  write ( eps_unit, '(a)' ) '%%Title: ' // trim ( file_name )
  write ( eps_unit, '(a)' ) '%%Pages: 1'
  write ( eps_unit, '(a)' ) '%%BoundingBox:    36    36   576   756'
  write ( eps_unit, '(a)' ) '%%Document-Fonts: Times-Roman'
  write ( eps_unit, '(a)' ) '%%LanguageLevel: 1'
  write ( eps_unit, '(a)' ) '%%EndComments'
  write ( eps_unit, '(a)' ) '%%BeginProlog'
  write ( eps_unit, '(a)' ) '/inch {72 mul} def'
  write ( eps_unit, '(a)' ) '%%EndProlog'
  write ( eps_unit, '(a)' ) '%%Page:      1     1'
  write ( eps_unit, '(a)' ) 'save'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '% Set RGB line color.'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) ' 0.9000 0.9000 0.9000 setrgbcolor'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '% Draw a gray border around the page.'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) 'newpath'
  write ( eps_unit, '(a)' ) '    36   126 moveto'
  write ( eps_unit, '(a)' ) '   576   126 lineto'
  write ( eps_unit, '(a)' ) '   576   666 lineto'
  write ( eps_unit, '(a)' ) '    36   666 lineto'
  write ( eps_unit, '(a)' ) '    36   126 lineto'
  write ( eps_unit, '(a)' ) 'stroke'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '% Set RGB line color.'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) ' 0.0000 0.0000 0.0000 setrgbcolor'

  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '%  Label the plot:'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '/Times-Roman findfont 0.50 inch scalefont setfont'
  write ( eps_unit, '(a)' ) '    36   666 moveto'
  write ( eps_unit, '(a)' ) '(' // trim ( title ) // ') show'

  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '% Define a clipping polygon'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '    36   126 moveto'
  write ( eps_unit, '(a)' ) '   576   126 lineto'
  write ( eps_unit, '(a)' ) '   576   666 lineto'
  write ( eps_unit, '(a)' ) '    36   666 lineto'
  write ( eps_unit, '(a)' ) '    36   126 lineto'
  write ( eps_unit, '(a)' ) 'clip newpath'

  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '%  Draw U in RED'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) ' 0.8000 0.0000 0.0000 setrgbcolor'

  do i = 1, n

    eps_x = int &
      ( ( node_x_max - x(i)              ) *  61.0D+00   &
      + (            + x(i) - node_x_min ) * 551.0D+00 ) &
      / x_scale

    eps_y = int &
      ( ( node_y_max - u(i)              ) * 151.0D+00   &
      + (              u(i) - node_y_min ) * 641.0D+00 ) &
      / y_scale

    if ( i == 1 ) then

      write ( eps_unit, '(a,2x,i4,2x,i4,2x,a)' ) &
        'newpath', eps_x, eps_y, 'moveto'

    else if ( i < n ) then

      write ( eps_unit, '(i4,2x,i4,2x,a)' ) eps_x, eps_y, 'lineto'

    else if ( i == n ) then

      write ( eps_unit, '(i4,2x,i4,2x,a)' ) eps_x, eps_y, 'lineto  stroke'

    end if

  end do

  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '%  Draw V in BLUE'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) ' 0.0000 0.0000 0.8000 setrgbcolor'

  do i = 1, n

    eps_x = int &
      ( ( node_x_max - x(i)              ) *  61.0   &
      + (            + x(i) - node_x_min ) * 551.0 ) &
      / x_scale

    eps_y = int &
      ( ( node_y_max - v(i)              ) * 151.0   &
      + (              v(i) - node_y_min ) * 641.0 ) &
      / y_scale

    if ( i == 1 ) then

      write ( eps_unit, '(a,2x,i4,2x,i4,2x,a)' ) &
        'newpath', eps_x, eps_y, 'moveto'

    else if ( i < n ) then

      write ( eps_unit, '(i4,2x,i4,2x,a)' ) eps_x, eps_y, 'lineto'

    else if ( i == n ) then

      write ( eps_unit, '(i4,2x,i4,2x,a)' ) eps_x, eps_y, 'lineto  stroke'

    end if

  end do

  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) 'restore showpage'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '% End of page'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '%%Trailer'
  write ( eps_unit, '(a)' ) '%%EOF'

  close ( unit = eps_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'UV_PLOT_EPS:'
  write ( *, '(a)' ) '  An encapsulated PostScript file was created'
  write ( *, '(a)' ) '  containing an image of U(X) and V(X).'
  write ( *, '(a)' ) '  The file is named "' // trim ( file_name ) // '".'

  return
end
