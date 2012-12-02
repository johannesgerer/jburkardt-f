program main

!*****************************************************************************80
!
!! MAIN is the main program for RKF45_PRB.
!
!  Discussion:
!
!    RKF45_PRB tests the RKF45 ODE integrator.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 June 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'RKF45_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the RKF45 library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
  call test05 ( )
  call test06 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'RKF45_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 solves a scalar ODE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 June 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: neqn = 1

  real ( kind = 4 ) abserr
  integer ( kind = 4 ) flag
  integer ( kind = 4 ) i_step
  integer ( kind = 4 ) n_step
  external r4_f1
  real ( kind = 4 ) r4_y1x
  real ( kind = 4 ) relerr
  real ( kind = 4 ) t
  real ( kind = 4 ) t_out
  real ( kind = 4 ) t_start
  real ( kind = 4 ) t_stop
  real ( kind = 4 ) y(neqn)
  real ( kind = 4 ) yp(neqn)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Solve a scalar equation using R4_RKF45:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Y'' = 0.25 * Y * ( 1 - Y / 20 )'

  abserr = sqrt ( epsilon ( abserr ) )
  relerr = sqrt ( epsilon ( relerr ) )

  flag = 1

  t_start = 0.0E+00
  t_stop = 20.0E+00

  n_step = 5

  t_out = 0.0E+00
  t = t_out
  y(1) = 1.0E+00
  call r4_f1 ( t, y, yp )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  FLAG     T             Y            Y''           Y_Exact         Error'
  write ( *, '(a)' ) ' '
  write ( *, '(i4,2x,5g14.6)' ) flag, t, y(1), yp(1), r4_y1x ( t ), &
    y(1) - r4_y1x ( t )

  do i_step = 1, n_step

    t = ( real ( n_step - i_step + 1, kind = 4 ) * t_start  &
        + real (          i_step - 1, kind = 4 ) * t_stop ) & 
        / real ( n_step,              kind = 4 )

    t_out = ( real ( n_step - i_step, kind = 4 ) * t_start  &
            + real (          i_step, kind = 4 ) * t_stop ) & 
            / real ( n_step,          kind = 4 )

    call r4_rkf45 ( r4_f1, neqn, y, yp, t, t_out, relerr, abserr, flag )

    write ( *, '(i4,2x,5g14.6)' ) flag, t, y(1), yp(1), r4_y1x ( t ), &
      y(1) - r4_y1x ( t )

  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 solves a vector ODE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 June 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: neqn = 2

  real ( kind = 4 ) abserr
  integer ( kind = 4 ) flag
  integer ( kind = 4 ) i_step
  integer ( kind = 4 ) n_step
  external r4_f2
  real ( kind = 4 ) relerr
  real ( kind = 4 ) t
  real ( kind = 4 ) t_out
  real ( kind = 4 ) t_start
  real ( kind = 4 ) t_stop
  real ( kind = 4 ) y(neqn)
  real ( kind = 4 ) yp(neqn)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  Solve a vector equation using R4_RKF45:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Y''(1) =  Y(2)'
  write ( *, '(a)' ) '  Y''(2) = -Y(1)'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  This system is equivalent to the following'
  write ( *, '(a)' ) '  second order system:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Z" = - Z.'

  abserr = sqrt ( epsilon ( abserr ) )
  relerr = sqrt ( epsilon ( relerr ) )

  flag = 1

  t_start = 0.0E+00
  t_stop = 2.0E+00 * 3.14159265E+00

  n_step = 12

  t = 0.0E+00
  t_out = 0.0E+00

  y(1) = 1.0E+00
  y(2) = 0.0E+00
  call r4_f2 ( t, y, yp )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  FLAG       T          Y(1)          Y(2)'
  write ( *, '(a)' ) ' '
  write ( *, '(i4,2x,4g14.6)' ) flag, t, y(1), y(2)

  do i_step = 1, n_step

    t = ( real ( n_step - i_step + 1, kind = 4 ) * t_start &
        + real (          i_step - 1, kind = 4 ) * t_stop ) &
        / real ( n_step,              kind = 4 )

    t_out = ( real ( n_step - i_step, kind = 4 ) * t_start &
            + real (          i_step, kind = 4 ) * t_stop ) &
            / real ( n_step,          kind = 4 )

    call r4_rkf45 ( r4_f2, neqn, y, yp, t, t_out, relerr, abserr, flag )

    write ( *, '(i4,2x,4g14.6)' ) flag, t, y(1), y(2)

  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 solves a scalar ODE and uses one-step integration.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 June 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: neqn = 1

  real ( kind = 4 ) abserr
  integer ( kind = 4 ) flag
  integer ( kind = 4 ) i_step
  integer ( kind = 4 ) n_step
  external r4_f1
  real ( kind = 4 ) r4_y1x
  real ( kind = 4 ) relerr
  real ( kind = 4 ) t
  real ( kind = 4 ) t_out
  real ( kind = 4 ) t_start
  real ( kind = 4 ) t_stop
  real ( kind = 4 ) y(neqn)
  real ( kind = 4 ) yp(neqn)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  Solve a scalar equation using R4_RKF45:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Y'' = 0.25 * Y * ( 1 - Y / 20 )'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Use the special SINGLE_STEP mode'
  write ( *, '(a)' ) '  which returns after every step.'

  abserr = sqrt ( epsilon ( abserr ) )
  relerr = sqrt ( epsilon ( relerr ) )

  flag = -1

  t_start = 0.0E+00
  t_stop = 20.0E+00

  n_step = 5

  t = 0.0E+00
  t_out = 0.0E+00
  y(1) = 1.0E+00
  call r4_f1 ( t, y, yp )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  FLAG     T             Y           Y''         Y_Exact        Error'
  write ( *, '(a)' ) ' '
  write ( *, '(i4,2x,5g14.6)' ) flag, t, y(1), yp(1), r4_y1x ( t ), &
    y(1) - r4_y1x ( t )

  do i_step = 1, n_step

    t = ( real ( n_step - i_step + 1, kind = 4 ) * t_start  &
        + real (          i_step - 1, kind = 4 ) * t_stop ) & 
        / real ( n_step,              kind = 4 )

    t_out = ( real ( n_step - i_step, kind = 4 ) * t_start  &
            + real (          i_step, kind = 4 ) * t_stop ) & 
            / real ( n_step,          kind = 4 )
!
!  As long as FLAG is negative, we are heading towards T_OUT, but
!  have not reached it!
!
    do while ( flag < 0 )

      call r4_rkf45 ( r4_f1, neqn, y, yp, t, t_out, relerr, abserr, flag )

      write ( *, '(i4,2x,5g14.6)' ) flag, t, y(1), yp(1), r4_y1x ( t ), &
        y(1) - r4_y1x ( t )

    end do
!
!  FLAG is returned as +2 when we reach T_OUT.  Reset it to -2
!  to continue to the next T_OUT in one step mode.
!
    flag = -2

  end do

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 solves a scalar ODE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 June 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: neqn = 1

  real ( kind = 8 ) abserr
  integer ( kind = 4 ) flag
  integer ( kind = 4 ) i_step
  integer ( kind = 4 ) n_step
  external r8_f1
  real ( kind = 8 ) r8_y1x
  real ( kind = 8 ) relerr
  real ( kind = 8 ) t
  real ( kind = 8 ) t_out
  real ( kind = 8 ) t_start
  real ( kind = 8 ) t_stop
  real ( kind = 8 ) y(neqn)
  real ( kind = 8 ) yp(neqn)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  Solve a scalar equation using R8_RKF45:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Y'' = 0.25 * Y * ( 1 - Y / 20 )'

  abserr = sqrt ( epsilon ( abserr ) )
  relerr = sqrt ( epsilon ( relerr ) )

  flag = 1

  t_start = 0.0D+00
  t_stop = 20.0D+00

  n_step = 5

  t_out = 0.0D+00
  t = t_out
  y(1) = 1.0D+00
  call r8_f1 ( t, y, yp )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  FLAG     T             Y            Y''           Y_Exact         Error'
  write ( *, '(a)' ) ' '
  write ( *, '(i4,2x,5g14.6)' ) flag, t, y(1), yp(1), r8_y1x ( t ), &
    y(1) - r8_y1x ( t )

  do i_step = 1, n_step

    t = ( real ( n_step - i_step + 1, kind = 8 ) * t_start  &
        + real (          i_step - 1, kind = 8 ) * t_stop ) & 
        / real ( n_step,              kind = 8 )

    t_out = ( real ( n_step - i_step, kind = 8 ) * t_start  &
            + real (          i_step, kind = 8 ) * t_stop ) & 
            / real ( n_step,          kind = 8 )

    call r8_rkf45 ( r8_f1, neqn, y, yp, t, t_out, relerr, abserr, flag )

    write ( *, '(i4,2x,5g14.6)' ) flag, t, y(1), yp(1), r8_y1x ( t ), &
      y(1) - r8_y1x ( t )

  end do

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 solves a vector ODE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 June 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: neqn = 2

  real ( kind = 8 ) abserr
  integer ( kind = 4 ) flag
  integer ( kind = 4 ) i_step
  integer ( kind = 4 ) n_step
  external r8_f2
  real ( kind = 8 ) relerr
  real ( kind = 8 ) t
  real ( kind = 8 ) t_out
  real ( kind = 8 ) t_start
  real ( kind = 8 ) t_stop
  real ( kind = 8 ) y(neqn)
  real ( kind = 8 ) yp(neqn)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  Solve a vector equation using R8_RKF45:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Y''(1) =  Y(2)'
  write ( *, '(a)' ) '  Y''(2) = -Y(1)'

  abserr = sqrt ( epsilon ( abserr ) )
  relerr = sqrt ( epsilon ( relerr ) )

  flag = 1

  t_start = 0.0D+00
  t_stop = 2.0D+00 * 3.14159265D+00

  n_step = 12

  t = 0.0D+00
  t_out = 0.0D+00

  y(1) = 1.0D+00
  y(2) = 0.0D+00
  call r8_f2 ( t, y, yp )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  FLAG       T          Y(1)          Y(2)'
  write ( *, '(a)' ) ' '
  write ( *, '(i4,2x,4g14.6)' ) flag, t, y(1), y(2)

  do i_step = 1, n_step

    t = ( real ( n_step - i_step + 1, kind = 8 ) * t_start &
        + real (          i_step - 1, kind = 8 ) * t_stop ) &
        / real ( n_step,              kind = 8 )

    t_out = ( real ( n_step - i_step, kind = 8 ) * t_start &
            + real (          i_step, kind = 8 ) * t_stop ) &
            / real ( n_step,          kind = 8 )

    call r8_rkf45 ( r8_f2, neqn, y, yp, t, t_out, relerr, abserr, flag )

    write ( *, '(i4,2x,4g14.6)' ) flag, t, y(1), y(2)

  end do

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 solves a scalar ODE and uses one-step integration.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 June 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: neqn = 1

  real ( kind = 8 ) abserr
  integer ( kind = 4 ) flag
  integer ( kind = 4 ) i_step
  integer ( kind = 4 ) n_step
  external r8_f1
  real ( kind = 8 ) r8_y1x
  real ( kind = 8 ) relerr
  real ( kind = 8 ) t
  real ( kind = 8 ) t_out
  real ( kind = 8 ) t_start
  real ( kind = 8 ) t_stop
  real ( kind = 8 ) y(neqn)
  real ( kind = 8 ) yp(neqn)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  Solve a scalar equation using R8_RKF45:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Y'' = 0.25 * Y * ( 1 - Y / 20 )'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Use the special SINGLE_STEP mode'
  write ( *, '(a)' ) '  which returns after every step.'

  abserr = sqrt ( epsilon ( abserr ) )
  relerr = sqrt ( epsilon ( relerr ) )

  flag = -1

  t_start = 0.0D+00
  t_stop = 20.0D+00

  n_step = 5

  t = 0.0D+00
  t_out = 0.0D+00
  y(1) = 1.0D+00
  call r8_f1 ( t, y, yp )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  FLAG     T             Y           Y''       Y_Exact        Error'
  write ( *, '(a)' ) ' '
  write ( *, '(i4,2x,5g14.6)' ) flag, t, y(1), yp(1), r8_y1x ( t ), &
    y(1) - r8_y1x ( t )

  do i_step = 1, n_step

    t = ( real ( n_step - i_step + 1, kind = 8 ) * t_start  &
        + real (          i_step - 1, kind = 8 ) * t_stop ) & 
        / real ( n_step,              kind = 8 )

    t_out = ( real ( n_step - i_step, kind = 8 ) * t_start  &
            + real (          i_step, kind = 8 ) * t_stop ) & 
            / real ( n_step,          kind = 8 )
!
!  As long as FLAG is negative, we are heading towards T_OUT, but
!  have not reached it!
!
    do while ( flag < 0 )

      call r8_rkf45 ( r8_f1, neqn, y, yp, t, t_out, relerr, abserr, flag )

      write ( *, '(i4,2x,5g14.6)' ) flag, t, y(1), yp(1), r8_y1x ( t ), &
        y(1) - r8_y1x ( t )

    end do
!
!  FLAG is returned as +2 when we reach T_OUT.  Reset it to -2
!  to continue to the next T_OUT in one step mode.
!
    flag = -2

  end do

  return
end
subroutine r4_f1 ( t, y, yp )

!*****************************************************************************80
!
!! R4_F1 evaluates the derivative for the ODE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 4 ) T, the value of the independent variable.
!
!    Input, real ( kind = 4 ) Y(NEQN), the value of the dependent variable.
!
!    Output, real ( kind = 4 ) YP(NEQN), the value of the derivative
!    dY(1:NEQN)/dT.
!
  implicit none

  real ( kind = 4 ) t
  real ( kind = 4 ) y(1)
  real ( kind = 4 ) yp(1)

  yp(1) = 0.25E+00 * y(1) * ( 1.0E+00 - y(1) / 20.0E+00 )

  return
end
function r4_y1x ( t )

!*****************************************************************************80
!
!! R4_Y1X evaluates the exact solution of the ODE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 4 ) T, the value of the independent variable.
!
!    Output, real ( kind = 4 ) Y1X_S, the exact solution.
!
  implicit none

  real ( kind = 4 ) t
  real ( kind = 4 ) r4_y1x

  r4_y1x = 20.0E+00 / ( 1.0E+00 + 19.0E+00 * exp ( - 0.25E+00 * t ) )

  return
end
subroutine r4_f2 ( t, y, yp )

!*****************************************************************************80
!
!! R4_F2 evaluates the derivative for the ODE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 4 ) T, the value of the independent variable.
!
!    Input, real ( kind = 4 ) Y(NEQN), the value of the dependent variable.
!
!    Output, real ( kind = 4 ) YP(NEQN), the value of the derivative
!    dY(1:NEQN)/dT.
!
  implicit none

  real ( kind = 4 ) t
  real ( kind = 4 ) y(2)
  real ( kind = 4 ) yp(2)

  yp(1) =  y(2)
  yp(2) = -y(1)

  return
end
subroutine r8_f1 ( t, y, yp )

!*****************************************************************************80
!
!! R8_F1 evaluates the derivative for the ODE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T, the value of the independent variable.
!
!    Input, real ( kind = 8 ) Y(NEQN), the value of the dependent variable.
!
!    Output, real ( kind = 8 ) YP(NEQN), the value of the derivative
!    dY(1:NEQN)/dT.
!
  implicit none

  real ( kind = 8 ) t
  real ( kind = 8 ) y(1)
  real ( kind = 8 ) yp(1)

  yp(1) = 0.25D+00 * y(1) * ( 1.0D+00 - y(1) / 20.0D+00 )

  return
end
function r8_y1x ( t )

!*****************************************************************************80
!
!! R8_Y1X evaluates the exact solution of the ODE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T, the value of the independent variable.
!
!    Output, real ( kind = 8 ) Y1X_D, the exact solution.
!
  implicit none

  real ( kind = 8 ) t
  real ( kind = 8 ) r8_y1x

  r8_y1x = 20.0D+00 / ( 1.0D+00 + 19.0D+00 * exp ( - 0.25D+00 * t ) )

  return
end
subroutine r8_f2 ( t, y, yp )

!*****************************************************************************80
!
!! R8_F2 evaluates the derivative for the ODE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T, the value of the independent variable.
!
!    Input, real ( kind = 8 ) Y(NEQN), the value of the dependent variable.
!
!    Output, real ( kind = 8 ) YP(NEQN), the value of the derivative
!    dY(1:NEQN)/dT.
!
  implicit none

  real ( kind = 8 ) t
  real ( kind = 8 ) y(2)
  real ( kind = 8 ) yp(2)

  yp(1) =  y(2)
  yp(2) = -y(1)

  return
end
