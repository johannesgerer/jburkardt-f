program main

!*****************************************************************************80
!
!! MAIN is the main program for FD1D_HEAT_IMPLICIT.
!
!  Discussion:
!
!    FD1D_HEAT_IMPLICIT solves the 1D heat equation with an implicit method.
!
!    This program solves
!
!      dUdT - k * d2UdX2 = F(X,T)
!
!    over the interval [A,B] with boundary conditions
!
!      U(A,T) = UA(T),
!      U(B,T) = UB(T),
!
!    over the time interval [T0,T1] with initial conditions
!
!      U(X,T0) = U0(X)
!
!    The code uses the finite difference method to approximate the
!    second derivative in space, and an implicit backward Euler approximation
!    to the first derivative in time.
!
!    The finite difference form can be written as
!
!      U(X,T+dt) - U(X,T)                  ( U(X-dx,T+dt) - 2 U(X,T+dt) + U(X+dx,T+dt) )
!      ------------------ = F(X,T+dt) + k *  --------------------------------------
!               dt                                   dx * dx
!
!    so that we have the following linear system for the values of U at time T+dt:
!
!            -     k * dt / dx / dx   * U(X-dt,T+dt)
!      + ( 1 + 2 * k * dt / dx / dx ) * U(X,   T+dt)
!            -     k * dt / dx / dx   * U(X+dt,T+dt)
!      =               dt             * F(X,   T+dt)
!      +                                U(X,   T)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 May 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable :: a(:,:)
  real ( kind = 8 ), allocatable :: b(:)
  real ( kind = 8 ), allocatable :: fvec(:)
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) info
  integer   ( kind = 4 ) j
  integer   ( kind = 4 ) job
  real ( kind = 8 ) k
  real ( kind = 8 ), allocatable :: t(:)
  real ( kind = 8 ) t_delt
  character ( len = 80 ) t_file
  real ( kind = 8 ) t_max
  real ( kind = 8 ) t_min
  integer   ( kind = 4 ) t_num
  integer   ( kind = 4 ) t_unit
  real ( kind = 8 ), allocatable :: u(:,:)
  character ( len = 80 ) u_file
  integer   ( kind = 4 ) u_unit
  real ( kind = 8 ) w
  real ( kind = 8 ), allocatable :: x(:)
  real ( kind = 8 ) x_delt
  character ( len = 80 ) x_file
  real ( kind = 8 ) x_max
  real ( kind = 8 ) x_min
  integer   ( kind = 4 ) x_num
  integer   ( kind = 4 ) x_unit

  call timestamp ( );
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FD1D_HEAT_IMPLICIT'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Finite difference solution of'
  write ( *, '(a)' ) '  the time dependent 1D heat equation'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    Ut - k * Uxx = F(x,t)'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  for space interval A <= X <= B with boundary conditions'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    U(A,t) = UA(t)'
  write ( *, '(a)' ) '    U(B,t) = UB(t)'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  and time interval T0 <= T <= T1 with initial condition'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    U(X,T0) = U0(X).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A second order difference approximation is used for Uxx.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A first order backward Euler difference approximation'
  write ( *, '(a)' ) '  is used for Ut.'

  k = 5.0D-07
!
!  Set X values.
!
  x_min = 0.0D+00
  x_max = 0.3D+00
  x_num = 11
  x_delt = ( x_max - x_min ) / real ( x_num - 1, kind = 8 )

  allocate ( x(1:x_num) )

  do i = 1, x_num
    x(i) = ( real ( x_num - i,     kind = 8 ) * x_min   &
           + real (         i - 1, kind = 8 ) * x_max ) &
           / real ( x_num     - 1, kind = 8 )
  end do
! 
!  Set T values.
!
  t_min = 0.0D+00
  t_max = 22000.0D+00
  t_num = 51
  t_delt = ( t_max - t_min ) / real ( t_num - 1, kind = 8 )

  allocate ( t(1:t_num) )

  do j = 1, t_num

    t(j) = ( real ( t_num - j,     kind = 8 ) * t_min   &
           + real (         j - 1, kind = 8 ) * t_max ) &
           / real ( t_num     - 1, kind = 8 )
  end do
!
!  Set the initial data, for time T_MIN.
!
  allocate ( u(1:x_num,1:t_num) )

  call u0 ( x_min, x_max, t_min, x_num, x, u(1:x_num,1) )
!
!  The matrix A does not change with time.  We can set it once,
!  factor it once, and solve repeatedly.
!
  w = k * t_delt / x_delt / x_delt

  allocate ( a(1:3,1:x_num) )

  a(1,1) = 0.0D+00

  a(2,1) = 1.0D+00
  a(1,2) = 0.0D+00

  do i = 2, x_num - 1
    a(3,i-1) =                   - w
    a(2,i  ) = 1.0D+00 + 2.0D+00 * w
    a(1,i+1) =                   - w
  end do

  a(3,x_num-1) = 0.0D+00
  a(2,x_num) = 1.0D+00

  a(3,x_num) = 0.0D+00
!
!  Factor the matrix.
!
  call r83_np_fa ( x_num, a, info )

  allocate ( b(1:x_num) )
  allocate ( fvec(1:x_num) )

  do j = 2, t_num
!
!  Set the right hand side B.
!
    call ua ( x_min, x_max, t_min, t(j), b(1) )

    call f ( x_min, x_max, t_min, t(j), x_num, x, fvec )

    b(2:x_num-1) = u(2:x_num-1,j-1) + t_delt * fvec(2:x_num-1)

    call ub ( x_min, x_max, t_min, t(j), b(x_num) )

    job = 0
    call r83_np_sl ( x_num, a, b, job )

    u(1:x_num,j) = b(1:x_num)

  end do

  x_file = 'x.txt'
  call get_unit ( x_unit )
  open ( unit = x_unit, file = x_file, status = 'replace' )

  do i = 1, x_num
    write ( x_unit, '(g14.6)' ) x(i)
  end do

  close ( unit = x_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  X data written to "' // trim ( x_file ) // '".'

  t_file = 't.txt'
  call get_unit ( t_unit )
  open ( unit = t_unit, file = t_file, status = 'replace' )

  do j = 1, t_num
    write ( t_unit, '(g14.6)' ) t(j)
  end do

  close ( unit = t_unit )

  write ( *, '(a)' ) '  T data written to "' // trim ( t_file ) // '".'

  u_file = 'u.txt'
  call get_unit ( u_unit )
  open ( unit = u_unit, file = u_file, status = 'replace' )

  do j = 1, t_num
    write ( u_unit, '(11g14.6)' ) u(1:x_num,j)
  end do

  close ( unit = u_unit )

  write ( *, '(a)' ) '  U data written to "' // trim ( u_file ) // '".'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FD1D_HEAT_IMPLICIT'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) 
  call timestamp ( );

  deallocate ( a )
  deallocate ( b )
  deallocate ( fvec )
  deallocate ( t )
  deallocate ( u  )
  deallocate ( x )

  return
end
subroutine f ( a, b, t0, t, n, x, value )

!*****************************************************************************80
!
!! F returns the right hand side of the heat equation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the left and right endpoints.
!
!    Input, real ( kind = 8 ) T0, the initial time.
!
!    Input, real ( kind = 8 ) T, the current time.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) X(N), the current spatial positions.
!
!    Output, real ( kind = 8 ) VALUE(N), the prescribed value of U(X(:),T0).
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) t
  real ( kind = 8 ) t0
  real ( kind = 8 ) value(n)
  real ( kind = 8 ) x(n)

  value(1:n) = 0.0D+00

  return
end
subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is an integer between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is an integer between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IUNIT, the free unit number.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  logical              lopen

  iunit = 0

  do i = 1, 99

    if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then

      inquire ( unit = i, opened = lopen, iostat = ios )

      if ( ios == 0 ) then
        if ( .not. lopen ) then
          iunit = i
          return
        end if
      end if

    end if

  end do

  return
end
subroutine r83_np_fa ( n, a, info )

!*****************************************************************************80
!
!! R83_NP_FA factors an R83 matrix without pivoting.
!
!  Discussion:
!
!    The R83 storage format is used for a tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:N), the diagonal in
!    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
!    original matrix is "collapsed" vertically into the array.
!
!    Because this routine does not use pivoting, it can fail even when
!    the matrix is not singular, and it is liable to make larger
!    errors.
!
!    R83_NP_FA and R83_NP_SL may be preferable to the corresponding
!    LINPACK routine SGTSL for tridiagonal systems, which factors and solves
!    in one step, and does not save the factorization.
!
!  Example:
!
!    Here is how an R83 matrix of order 5 would be stored:
!
!       *  A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54  *
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
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
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be at least 2.
!
!    Input/output, real ( kind = 8 ) A(3,N).
!    On input, the tridiagonal matrix.  On output, factorization information.
!
!    Output, integer ( kind = 4 ) INFO, singularity flag.
!    0, no singularity detected.
!    nonzero, the factorization failed on the INFO-th step.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(3,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info

  info = 0

  do i = 1, n-1

    if ( a(2,i) == 0.0D+00 ) then
      info = i
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R83_NP_FA - Fatal error!'
      write ( *, '(a,i8)' ) '  Zero pivot on step ', info
      stop
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
    write ( *, '(a)' ) 'R83_NP_FA - Fatal error!'
    write ( *, '(a,i8)' ) '  Zero pivot on step ', info
    stop
  end if

  return
end
subroutine r83_np_sl ( n, a_lu, b, job )

!*****************************************************************************80
!
!! R83_NP_SL solves an R83 system factored by R83_NP_FA.
!
!  Discussion:
!
!    The R83 storage format is used for a tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:N), the diagonal in
!    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
!    original matrix is "collapsed" vertically into the array.
!
!  Example:
!
!    Here is how an R83 matrix of order 5 would be stored:
!
!       *  A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54  *
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
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
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be at least 2.
!
!    Input, real ( kind = 8 ) A_LU(3,N), the LU factors from R83_NP_FA.
!
!    Input/output, real ( kind = 8 ) B(N).
!    On input, B contains the right hand side of the linear system.
!    On output, B contains the solution of the linear system.
!
!    Input, integer ( kind = 4 ) JOB, specifies the system to solve.
!    0, solve A * x = b.
!    nonzero, solve A' * x = b.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a_lu(3,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) job

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
subroutine u0 ( a, b, t0, n, x, value )

!*****************************************************************************80
!
!! U0 returns the initial condition at the starting time.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the left and right endpoints
!
!    Input, real ( kind = 8 ) T0, the initial time.
!
!    Input, real ( kind = 8 ) T, the current time.
!
!    Input, integer ( kind = 4 ) N, the number of points where initial data is needed.
!
!    Input, real ( kind = 8 ) X(N), the positions where initial data is needed.
!
!    Output, real ( kind = 8 ) VALUE(N), the prescribed value of U(X,T0).
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) t
  real ( kind = 8 ) t0
  real ( kind = 8 ) value(n)
  real ( kind = 8 ) x(n)

  value(1:n) = 100.0D+00

  return
end
subroutine ua ( a, b, t0, t, value )

!*****************************************************************************80
!
!! UA returns the Dirichlet boundary condition at the left endpoint.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the left and right endpoints
!
!    Input, real ( kind = 8 ) T0, the initial time.
!
!    Input, real ( kind = 8 ) T, the current time.
!
!    Output, real ( kind = 8 ) VALUE, the prescribed value of U(A,T).
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) t
  real ( kind = 8 ) t0
  real ( kind = 8 ) value

  value = 20.0D+00

  return
end
subroutine ub ( a, b, t0, t, value )

!*****************************************************************************80
!
!! UB returns the Dirichlet boundary condition at the right endpoint.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the left and right endpoints
!
!    Input, real ( kind = 8 ) T0, the initial time.
!
!    Input, real ( kind = 8 ) T, the current time.
!
!    Output, real ( kind = 8 ) VALUE, the prescribed value of U(B,T).
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) t
  real ( kind = 8 ) t0
  real ( kind = 8 ) value

  value = 20.0D+00

  return
end
