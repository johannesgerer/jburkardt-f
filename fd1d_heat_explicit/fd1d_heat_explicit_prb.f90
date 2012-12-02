program main

!*****************************************************************************80
!
!! FD1D_HEAT_EXPLICIT_TEST tests the FD1D_HEAT_EXPLICIT library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 January 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FD1D_HEAT_EXPLICIT_TEST:'
  write ( *, '(a)' ) '  FORTRAN90 version.'
  write ( *, '(a)' ) '  Test the FD1D_HEAT_EXPLICIT library.'

  call fd1d_heat_explicit_test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FD1D_HEAT_EXPLICIT_TEST:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  return
end
subroutine fd1d_heat_explicit_test01 ( )

!*****************************************************************************80
!
!! FD1D_HEAT_EXPLICIT_TEST01 does a simple test problem
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 January 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  external bc_test01
  real ( kind = 8 ) cfl
  real ( kind = 8 ) dt
  real ( kind = 8 ), allocatable :: h(:)
  real ( kind = 8 ), allocatable :: h_new(:)
  real ( kind = 8 ), allocatable :: hmat(:,:)
  integer ( kind = 4 ) j
  real ( kind = 8 ) k
  external rhs_test01
  real ( kind = 8 ), allocatable :: t(:)
  real ( kind = 8 ) t_max
  real ( kind = 8 ) t_min
  integer ( kind = 4 ) t_num
  real ( kind = 8 ), allocatable :: x(:)
  real ( kind = 8 ) x_max
  real ( kind = 8 ) x_min
  integer ( kind = 4 ) x_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FD1D_HEAT_EXPLICIT_TEST01:'
  write ( *, '(a)' ) '  Compute an approximate solution to the time-dependent'
  write ( *, '(a)' ) '  one dimensional heat equation:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    dH/dt - K * d2H/dx2 = f(x,t)'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Run a simple test case.'
!
!  Heat coefficient.
!
  k = 0.002D+00
!
!  X_NUM is the number of equally spaced nodes to use between 0 and 1.
!
  x_num = 21
  x_min = 0.0D+00
  x_max = 1.0D+00
  allocate ( x(1:x_num) )
  call r8vec_linspace ( x_num, x_min, x_max, x )
!
!  T_NUM is the number of equally spaced time points between 0 and 10.0.
!
  t_num = 201
  t_min = 0.0D+00
  t_max = 80.0D+00
  dt = ( t_max - t_min ) / real ( t_num - 1, kind = 8 )
  allocate ( t(1:t_num) )
  call r8vec_linspace ( t_num, t_min, t_max, t )
!
!  Get the CFL coefficient.
!
  call fd1d_heat_explicit_cfl ( k, t_num, t_min, t_max, x_num, x_min, x_max, cfl )
!
!  Running the code produces an array H of temperatures H(t,x),
!  and vectors x and t.
!
  allocate ( h(1:x_num) )
  allocate ( h_new(1:x_num) )

  call ic_test01 ( x_num, x, t(1), h )
  call bc_test01 ( x_num, x, t(1), h )

  allocate ( hmat(1:x_num,1:t_num) )
  hmat(1:x_num,1) = h(1:x_num)

  do j = 2, t_num
    call fd1d_heat_explicit ( x_num, x, t(j-1), dt, cfl, rhs_test01, &
      bc_test01, h, h_new )
    hmat(1:x_num,j) = h_new(1:x_num)
    h(1:x_num) = h_new(1:x_num)
  end do
!
!  Write the data to files.
!
  call r8mat_write ( 'h_test01.txt', x_num, t_num, hmat )
  call r8vec_write ( 't_test01.txt', t_num, t )
  call r8vec_write ( 'x_test01.txt', x_num, x )

  deallocate ( h )
  deallocate ( h_new )
  deallocate ( hmat )
  deallocate ( t )
  deallocate ( x )

  return
end
subroutine bc_test01 ( x_num, x, t, h )

!*****************************************************************************80
!
!! BC_TEST01 evaluates the boundary conditions for problem 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) X_NUM, the number of nodes.
!
!    Input, real ( kind = 8 ) X(X_NUM), the node coordinates.
!
!    Input, real ( kind = 8 ) T, the current time.
!
!    Input, real ( kind = 8 ) H(X_NUM), the current heat values.
!
!    Output, real ( kind = 8 ) H(X_NUM), the current heat values, after boundary
!    conditions have been imposed.
!
  implicit none

  integer ( kind = 4 ) x_num

  real ( kind = 8 ) h(x_num)
  real ( kind = 8 ) t
  real ( kind = 8 ) x(x_num)

  h(1)  = 90.0D+00
  h(x_num) = 70.0D+00

  return
end
subroutine ic_test01 ( x_num, x, t, h  )

!*****************************************************************************80
!
!! IC_TEST01 evaluates the initial condition for problem 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) X_NUM, the number of nodes.
!
!    Input, real ( kind = 8 ) X(X_NUM), the node coordinates.
!
!    Input, real ( kind = 8 ) T, the initial time.
!
!    Output, real ( kind = 8 ) H(X_NUM), the heat values at the initial time.
!
  implicit none

  integer ( kind = 4 ) x_num

  real ( kind = 8 ) h(x_num)
  real ( kind = 8 ) t
  real ( kind = 8 ) x(x_num)

  h(1:x_num) = 50.0D+00

  return
end
subroutine rhs_test01 ( x_num, x, t, value )

!*****************************************************************************80
!
!! RHS_TEST01 evaluates the right hand side for problem 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) X_NUM, the number of nodes.
!
!    Input, real ( kind = 8 ) X(X_NUM), the node coordinates.
!
!    Input, real ( kind = 8 ) T, the current time.
!
!    Output, real ( kind = 8 ) VALUE(X_NUM), the source term.
!
  implicit none

  integer ( kind = 4 ) x_num

  real ( kind = 8 ) t
  real ( kind = 8 ) value(x_num)
  real ( kind = 8 ) x(x_num)

  value(1:x_num) = 0.0D+00

  return
end
