program main

!*****************************************************************************80
!
!! MAIN demonstrates the TEST_ODE test problems.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 March 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_ODE_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the TEST_ODE library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_ODE_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 simply lists the problems with titles and sizes.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 March 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) neqn
  integer ( kind = 4 ) test
  integer ( kind = 4 ) test_num
  character ( len = 80 ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  List the problem titles and sizes.'
!
!  Find out how many test problems are available.
!
  call p00_test_num ( test_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6,a)' ) '  There are ', test_num, ' test problems.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test  Size  Title'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    call p00_title ( test, title )
    call p00_neqn ( test, neqn )
    write ( *, '(2x,i4,2x,i4,2x,a)' ) test, neqn, trim ( title )
  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 solves most of the problems using an Euler method.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 March 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: step_num = 500
  integer ( kind = 4 ) test
  integer ( kind = 4 ) test_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  Solve problems using an Euler method.'
  write ( *, '(a,i8)' ) '  The number of steps taken is ', step_num

  call p00_test_num ( test_num )

  write ( *, '(a,i8)' ) '  The number of tests available is ', test_num
!
!  Solve each problem.
!
  do test = 1, test_num

    if ( test == 32 .or. test == 36 .or. test == 37 ) then

    else
      call euler_test ( test, step_num )
    end if

  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 solves most of the problems using a Runge-Kutta method.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 March 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) order
  integer ( kind = 4 ), parameter :: step_num = 500
  integer ( kind = 4 ) test
  integer ( kind = 4 ) test_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  Solve problems using a Runge-Kutta method.'
  write ( *, '(a,i8)' ) '  The number of steps taken is ', step_num

  call p00_test_num ( test_num )

  write ( *, '(a,i8)' ) '  The number of tests available is ', test_num
!
!  Solve each problem.
!
  order = 3

  do test = 1, test_num

    call rk_test ( test, step_num, order )

  end do

  return
end
subroutine euler_test ( test, step_num )

!*****************************************************************************80
!
!! EULER_TEST uses the Euler method on a test problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer TEST, the number of the problem to be demonstrated.
!
!    Input, integer STEP_NUM, the number of steps to take.
!
  implicit none

  integer ( kind = 4 ) neqn
  logical p00_autonomous
  integer ( kind = 4 ) step
  integer ( kind = 4 ) step_num
  real ( kind = 8 ) t_start
  real ( kind = 8 ) t_stop
  real ( kind = 8 ) t0
  real ( kind = 8 ) t1
  integer ( kind = 4 ) test
  character ( len = 80 ) title
  real ( kind = 8 ), allocatable, dimension ( : ) :: y0
  real ( kind = 8 ), allocatable, dimension ( : ) :: y1
  real ( kind = 8 ) y_ave
  real ( kind = 8 ) y_max
  real ( kind = 8 ) y_min
  real ( kind = 8 ) y_norm
  real ( kind = 8 ), allocatable, dimension ( : ) :: y_start
  real ( kind = 8 ), allocatable, dimension ( : ) :: y_stop

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'EULER_TEST'
  write ( *, '(a,i3)' ) '  Problem number = ', test

  call p00_title ( test, title )

  write ( *, '(2x,a)' ) trim ( title )
!
!  Autonomous?
!
  if ( p00_autonomous ( test ) ) then
    write ( *, '(a)' ) '  The system is autonomous.'
  else
    write ( *, '(a)' ) '  The system is not autonomous.'
  end if
!
!  Get the number of equations.
!
  call p00_neqn ( test, neqn )

  write ( *, '(a,i3)' ) '  Number of equations is ', neqn

  if ( 8 < neqn ) then
    write ( *, '(a)' ) '  The system is large.'
    write ( *, '(a)' ) '  Print only MIN, AVERAGE, MAX, L2NORM'
  end if

  allocate ( y0(1:neqn) )
  allocate ( y1(1:neqn) )
  allocate ( y_start(1:neqn) )
  allocate ( y_stop(1:neqn) )
!
!  Get the starting point.
!
  call p00_start ( test, neqn, t_start, y_start )
!
!  Get the stopping point.
!
  call p00_stop ( test, neqn, t_stop, y_stop )
!
!  Print the stepsize.
!
  write ( *, '(a,g14.6)' ) '  Stepsize H = ', &
    ( t_stop - t_start) / real ( step_num, kind = 8 )

  write ( *, '(a)' ) ' '

  if ( neqn <= 4 ) then
    write ( *, '(g14.6,(4g14.6))' ) t_start, y_start(1:neqn)
  else if ( neqn <= 8 ) then
    write ( *, '(g14.6,4g14.6)' ) t_start, y_start(1:4)
    write ( *, '( 14x, 4g14.6)' )          y_start(5:neqn)
  else
    y_min = minval ( y_start(1:neqn) )
    y_ave = sum ( y_start(1:neqn) ) / real ( neqn, kind = 8 )
    y_max = maxval ( y_start(1:neqn) )
    y_norm = sqrt ( sum ( y_start(1:neqn)**2 ) )
    write ( *, '(g14.6,4g14.6)' ) t_start, y_min, y_ave, y_max, y_norm
  end if

  y0(1:neqn) = y_start(1:neqn)

  do step = 1, step_num

    t0 = ( real ( step_num - step + 1, kind = 8 ) * t_start &
         + real (            step - 1, kind = 8 ) * t_stop ) &
         / real ( step_num,            kind = 8 )

    t1 = ( real ( step_num - step, kind = 8 ) * t_start &
         + real (            step, kind = 8 ) * t_stop ) &
         / real ( step_num,        kind = 8 )

    call p00_euler_step ( test, neqn, t0, y0, t1, y1 )

    if ( mod ( 10 * step, step_num ) == 0 .or. step == step_num ) then
      if ( neqn <= 4 ) then
        write ( *, '(g14.6,(4g14.6))' ) t1, y1(1:neqn)
      else if ( neqn <= 8 ) then
        write ( *, '(g14.6,4g14.6)' ) t1, y1(1:4)
        write ( *, '( 14x, 4g14.6)' )     y1(5:neqn)
      else
        y_min = minval ( y1(1:neqn) )
        y_ave = sum ( y1(1:neqn) ) / real ( neqn, kind = 8 )
        y_max = maxval ( y1(1:neqn) )
        y_norm = sqrt ( sum ( y1(1:neqn)**2 ) )
        write ( *, '(g14.6,4g14.6)' ) t1, y_min, y_ave, y_max, y_norm
      end if
    end if

    y0(1:neqn) = y1(1:neqn)

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Expected final conditions:'
  write ( *, '(a)' ) ' '
  if ( neqn <= 4 ) then
    write ( *, '(g14.6,(4g14.6))' ) t_stop, y_stop(1:neqn)
  else if ( neqn <= 8 ) then
    write ( *, '(g14.6,4g14.6)' ) t_stop, y_stop(1:4)
    write ( *, '( 14x, 4g14.6)' )          y_stop(5:neqn)
  else
    y_min = minval ( y_stop(1:neqn) )
    y_ave = sum ( y_stop(1:neqn) ) / real ( neqn, kind = 8 )
    y_max = maxval ( y_stop(1:neqn) )
    y_norm = sqrt ( sum ( y_stop(1:neqn)**2 ) )
    write ( *, '(g14.6,4g14.6)' ) t_stop, y_min, y_ave, y_max, y_norm
  end if

  deallocate ( y0 )
  deallocate ( y1 )
  deallocate ( y_start )
  deallocate ( y_stop )

  return
end
subroutine rk_test ( test, step_num, order )

!*****************************************************************************80
!
!! RK_TEST uses a Runge-Kutta method on a test problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer TEST, the number of the problem to be demonstrated.
!
!    Input, integer STEP_NUM, the number of steps to take.
!
!    Input, integer ORDER, the order of the Runge-Kutta method to use.
!
  implicit none

  integer ( kind = 4 ) neqn
  integer ( kind = 4 ) order
  logical p00_autonomous
  integer ( kind = 4 ) step
  integer ( kind = 4 ) step_num
  real ( kind = 8 ) t_start
  real ( kind = 8 ) t_stop
  real ( kind = 8 ) t0
  real ( kind = 8 ) t1
  integer ( kind = 4 ) test
  character ( len = 80 ) title
  real ( kind = 8 ), allocatable, dimension ( : ) :: y0
  real ( kind = 8 ), allocatable, dimension ( : ) :: y1
  real ( kind = 8 ) y_ave
  real ( kind = 8 ) y_max
  real ( kind = 8 ) y_min
  real ( kind = 8 ) y_norm
  real ( kind = 8 ), allocatable, dimension ( : ) :: y_start
  real ( kind = 8 ), allocatable, dimension ( : ) :: y_stop

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'RK_TEST'
  write ( *, '(a,i3)' ) '  Problem number =    ', test
  write ( *, '(a,i3)' ) '  Runge Kutta order = ', order

  call p00_title ( test, title )

  write ( *, '(2x,a)' ) trim ( title )
!
!  Autonomous?
!
  if ( p00_autonomous ( test ) ) then
    write ( *, '(a)' ) '  The system is autonomous.'
  else
    write ( *, '(a)' ) '  The system is not autonomous.'
  end if
!
!  Get the number of equations.
!
  call p00_neqn ( test, neqn )

  write ( *, '(a,i3)' ) '  Number of equations is ', neqn

  if ( 8 < neqn ) then
    write ( *, '(a)' ) '  The system is large.'
    write ( *, '(a)' ) '  Print only MIN, AVERAGE, MAX, L2NORM'
  end if

  allocate ( y0(1:neqn) )
  allocate ( y1(1:neqn) )
  allocate ( y_start(1:neqn) )
  allocate ( y_stop(1:neqn) )
!
!  Get the starting point.
!
  call p00_start ( test, neqn, t_start, y_start )
!
!  Get the stopping point.
!
  call p00_stop ( test, neqn, t_stop, y_stop )
!
!  Print the stepsize.
!
  write ( *, '(a,g14.6)' ) '  Stepsize H = ', &
    ( t_stop - t_start) / real ( step_num, kind = 8 )

  write ( *, '(a)' ) ' '

  if ( neqn <= 4 ) then
    write ( *, '(g14.6,(4g14.6))' ) t_start, y_start(1:neqn)
  else if ( neqn <= 8 ) then
    write ( *, '(g14.6,4g14.6)' ) t_start, y_start(1:4)
    write ( *, '( 14x, 4g14.6)' )          y_start(5:neqn)
  else
    y_min = minval ( y_start(1:neqn) )
    y_ave = sum ( y_start(1:neqn) ) / real ( neqn, kind = 8 )
    y_max = maxval ( y_start(1:neqn) )
    y_norm = sqrt ( sum ( y_start(1:neqn)**2 ) )
    write ( *, '(g14.6,4g14.6)' ) t_start, y_min, y_ave, y_max, y_norm
  end if

  y0(1:neqn) = y_start(1:neqn)

  do step = 1, step_num

    t0 = ( real ( step_num - step + 1, kind = 8 ) * t_start &
         + real (            step - 1, kind = 8 ) * t_stop ) &
         / real ( step_num,            kind = 8 )

    t1 = ( real ( step_num - step, kind = 8 ) * t_start &
         + real (            step, kind = 8 ) * t_stop ) &
         / real ( step_num,        kind = 8 )

    call p00_rk_step ( test, neqn, order, t0, y0, t1, y1 )

    if ( mod ( 10 * step, step_num ) == 0 .or. step == step_num ) then
      if ( neqn <= 4 ) then
        write ( *, '(g14.6,(4g14.6))' ) t1, y1(1:neqn)
      else if ( neqn <= 8 ) then
        write ( *, '(g14.6,4g14.6)' ) t1, y1(1:4)
        write ( *, '( 14x, 4g14.6)' )     y1(5:neqn)
      else
        y_min = minval ( y1(1:neqn) )
        y_ave = sum ( y1(1:neqn) ) / real ( neqn, kind = 8 )
        y_max = maxval ( y1(1:neqn) )
        y_norm = sqrt ( sum ( y1(1:neqn)**2 ) )
        write ( *, '(g14.6,4g14.6)' ) t1, y_min, y_ave, y_max, y_norm
      end if
    end if

    y0(1:neqn) = y1(1:neqn)

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Expected final conditions:'
  write ( *, '(a)' ) ' '
  if ( neqn <= 4 ) then
    write ( *, '(g14.6,(4g14.6))' ) t_stop, y_stop(1:neqn)
  else if ( neqn <= 8 ) then
    write ( *, '(g14.6,4g14.6)' ) t_stop, y_stop(1:4)
    write ( *, '( 14x, 4g14.6)' )          y_stop(5:neqn)
  else
    y_min = minval ( y_stop(1:neqn) )
    y_ave = sum ( y_stop(1:neqn) ) / real ( neqn, kind = 8 )
    y_max = maxval ( y_stop(1:neqn) )
    y_norm = sqrt ( sum ( y_stop(1:neqn)**2 ) )
    write ( *, '(g14.6,4g14.6)' ) t_stop, y_min, y_ave, y_max, y_norm
  end if

  deallocate ( y0 )
  deallocate ( y1 )
  deallocate ( y_start )
  deallocate ( y_stop )

  return
end
