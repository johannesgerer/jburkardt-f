program main

!*****************************************************************************80
!
!! MAIN is the main program for POISSON_SERIAL.
!
!  Discussion:
!
!    POISSON_SERIAL is a program for solving the Poisson problem.
!
!    This program runs serially.  Its output is used as a benchmark for
!    comparison with similar programs run in a parallel environment.
!
!    The Poisson equation
!
!      - DEL^2 U(x,y) = F(x,y)
!
!    is solved on the unit square [0,1] x [0,1] using a grid of NX by
!    NX evenly spaced points.  The first and last points in each direction
!    are boundary points.
!
!    The boundary conditions and F are set so that the exact solution is
!
!      U(x,y) = sin ( pi * x * y )
!
!    so that
!
!      - DEL^2 U(x,y) = pi^2 * ( x^2 + y^2 ) * sin ( pi * x * y )
!
!    The Jacobi iteration is repeatedly applied until convergence is detected.
!
!    For convenience in writing the discretized equations, we assume that NX = NY.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 October 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nx = 11
  integer ( kind = 4 ), parameter :: ny = 11

  logical converged
  real ( kind = 8 ) diff
  real ( kind = 8 ) dx
  real ( kind = 8 ) dy
  real ( kind = 8 ) error
  real ( kind = 8 ) f(nx,ny)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) it
  integer ( kind = 4 ), parameter :: it_max = 1000
  integer ( kind = 4 ) j
  real ( kind = 8 ) r8mat_rms
  real ( kind = 8 ), parameter :: tolerance = 0.000001D+00
  real ( kind = 8 ) u(nx,ny)
  real ( kind = 8 ) u_exact
  real ( kind = 8 ) u_norm
  real ( kind = 8 ) udiff(nx,ny)
  real ( kind = 8 ) uexact(nx,ny)
  real ( kind = 8 ) unew(nx,ny)
  real ( kind = 8 ) unew_norm
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  dx = 1.0D+00 / real ( nx - 1, kind = 8 )
  dy = dx

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'POISSON_SERIAL:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  A program for solving the Poisson equation.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  -DEL^2 U = F(X,Y)'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  on the rectangle 0 <= X <= 1, 0 <= Y <= 1.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  F(X,Y) = pi^2 * ( x^2 + y^2 ) * sin ( pi * x * y )'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of interior X grid points is ', nx
  write ( *, '(a,i8)' ) '  The number of interior Y grid points is ', ny
  write ( *, '(a,f10.4)' ) '  The X grid spacing is ', dx
  write ( *, '(a,f10.4)' ) '  The Y grid spacing is ', dy
!
!  Initialize the data.
!
  call rhs ( nx, ny, f )
!
!  Set the initial solution estimate.
!  We are "allowed" to pick up the boundary conditions exactly.
!
  unew(1:nx,1:ny) = 0.0D+00
  unew(1,   1:ny) = f(1,   1:ny)
  unew(  nx,1:ny) = f(  nx,1:ny)
  unew(1:nx,1)    = f(1:nx,1)
  unew(1:nx,  ny) = f(1:nx,  ny)

  unew_norm = r8mat_rms ( nx, ny, unew )
!
!  Set up the exact solution.
!
  do j = 1, ny 
    y = real ( j - 1, kind = 8 ) / real ( ny - 1, kind = 8 )
    do i = 1, nx
      x = real ( i - 1, kind = 8 ) / real ( nx - 1, kind = 8 )
      uexact(i,j) = u_exact ( x, y )
    end do
  end do
  u_norm = r8mat_rms ( nx, ny, uexact )
  write ( *, '(a,g14.6)' ) '  RMS of exact solution = ', u_norm
!
!  Do the iteration.
!
  converged = .false.

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  Step    ||Unew||     ||Unew-U||     ||Unew-Exact||'
  write ( *, '(a)' ) ' '

  udiff = unew(1:nx,1:ny) - uexact(1:nx,1:ny)
  error = r8mat_rms ( nx, ny, udiff )
  write ( *, '(2x,i4,2x,g14.6,2x,14x,2x,g14.6)' ) 0, unew_norm, error

  do it = 1, it_max

    u(1:nx,1:ny) = unew(1:nx,1:ny)
!
!  UNEW is derived from U by one Jacobi step.
!
    call sweep ( nx, ny, dx, dy, f, u, unew )
!
!  Check for convergence.
!
    u_norm = unew_norm
    unew_norm = r8mat_rms ( nx, ny, unew )
    udiff = unew(1:nx,1:ny) - u(1:nx,1:ny)
    diff = r8mat_rms ( nx, ny, udiff )
    udiff = unew(1:nx,1:ny) - uexact(1:nx,1:ny)
    error = r8mat_rms ( nx, ny, udiff )

    write ( *, '(2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) it, unew_norm, diff, error

    if ( diff <= tolerance ) then
      converged = .true.
      exit
    end if

  end do

  if ( converged ) then
    write ( *, '(a)' ) '  The iteration has converged,'
  else
    write ( *, '(a)' ) '  The iteration has NOT converged,'
  end if
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'POISSON_SERIAL:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
function r8mat_rms ( m, n, a )

!*****************************************************************************80
!
!! R8MAT_RMS returns the root mean square of data stored as an R8MAT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns in A.
!
!    Input, real ( kind = 8 ) A(M,N), the data whose RMS is desired.
!
!    Output, real ( kind = 8 ) R8MAT_RMS, the root mean square of A.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) r8mat_rms

  r8mat_rms = sqrt ( sum ( a(1:m,1:n)**2 ) / real ( m * n, kind = 8 ) )

  return
end
subroutine rhs ( nx, ny, f )

!*****************************************************************************80
!
!! RHS initializes the right hand side "vector".
!
!  Discussion:
!
!    It is convenient for us to set up RHS as a 2D array.  However, each
!    entry of RHS is really the right hand side of a linear system of the
!    form
!
!      A * U = F
!
!    In cases where U(I,J) is a boundary value, then the equation is simply
!
!      U(I,J) = F(i,j)
!
!    and F(I,J) holds the boundary data.
!
!    Otherwise, the equation has the form
!
!      (1/DX^2) * ( U(I+1,J)+U(I-1,J)+U(I,J-1)+U(I,J+1)-4*U(I,J) ) = F(I,J)
!
!    where DX is the spacing and F(I,J) is the value at X(I), Y(J) of
!
!      pi^2 * ( x^2 + y^2 ) * sin ( pi * x * y )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 October 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NX, NY, the X and Y grid dimensions.
!
!    Output, real F(NX,NY), the right hand side data.
!
  implicit none

  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny

  real ( kind = 8 ) f(nx,ny)
  real ( kind = 8 ) fnorm
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) r8mat_rms
  real ( kind = 8 ) u_exact
  real ( kind = 8 ) uxxyy_exact
  real ( kind = 8 ) x
  real ( kind = 8 ) y
!
!  The "boundary" entries of F store the boundary values of the solution.
!  The "interior" entries of F store the right hand sides of the Poisson equation.
!
  do j = 1, ny
    y = real ( j - 1, kind = 8 ) / real ( ny - 1, kind = 8 )
    do i = 1, nx
      x = real ( i - 1, kind = 8 ) / real ( nx - 1, kind = 8 )
      if ( i == 1 .or. i == nx .or. j == 1 .or. j == ny ) then
        f(i,j) = u_exact ( x, y )
      else
        f(i,j) = - uxxyy_exact ( x, y )
      end if
    end do
  end do

  fnorm = r8mat_rms ( nx, ny, f ) 

  write ( *, '(a,g14.6)' ) '  RMS of F = ', fnorm

  return
end
subroutine sweep ( nx, ny, dx, dy, f, u, unew )

!*****************************************************************************80
!
!! SWEEP carries out one step of the Jacobi iteration.
!
!  Discussion:
!
!    Assuming DX = DY, we can approximate
!
!      - ( d/dx d/dx + d/dy d/dy ) U(X,Y) 
!
!    by
!
!      ( U(i-1,j) + U(i+1,j) + U(i,j-1) + U(i,j+1) - 4*U(i,j) ) / dx / dy
!
!    The discretization employed below will not be correct in the general
!    case where DX and DY are not equal.  It's only a little more complicated
!    to allow DX and DY to be different, but we're not going to worry about 
!    that right now.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 October 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NX, NY, the X and Y grid dimensions.
!
!    Input, real DX, DY, the spacing between grid points.
!
!    Input, real F(NX,NY), the right hand side data.
!
!    Input, real U(NX,NY), the previous solution estimate.
!
!    Output, real UNEW(NX,NY), the updated solution estimate.
!
  implicit none

  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny

  real ( kind = 8 ) dx
  real ( kind = 8 ) dy
  real ( kind = 8 ) f(nx,ny)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) u(nx,ny)
  real ( kind = 8 ) unew(nx,ny)
!
!  The interior values are set by a Jacobi iteration.
!
  do j = 1, ny
    do i = 1, nx

      if ( i == 1 .or. i == nx .or. j == 1 .or. j == ny ) then
        unew(i,j) = f(i,j)
      else
        unew(i,j) = 0.25D+00 * ( &
            u(i-1,j) + u(i,j+1) + u(i,j-1) + u(i+1,j) + f(i,j) * dx * dy )
      end if

    end do
  end do

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
function u_exact ( x, y )

!*****************************************************************************80
!
!! U_EXACT evaluates the exact solution.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 October 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, Y, the coordinates of a point.
!
!    Output, real ( kind = 8 ) U_EXACT, the value of the exact solution 
!    at (X,Y).
!
  implicit none

  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) u_exact
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  u_exact = sin ( pi * x * y )

  return
end
function uxxyy_exact ( x, y )

!*****************************************************************************80
!
!! UXXYY_EXACT evaluates ( d/dx d/dx + d/dy d/dy ) of the exact solution.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 October 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, Y, the coordinates of a point.
!
!    Output, real ( kind = 8 ) UXXYY_EXACT, the value of 
!    ( d/dx d/dx + d/dy d/dy ) of the exact solution at (X,Y).
!
  implicit none

  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) uxxyy_exact
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  uxxyy_exact = - pi * pi * ( x * x + y * y ) * sin ( pi * x * y )

  return
end
