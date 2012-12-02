subroutine fd1d_wave_alpha ( x_num, x1, x2, t_num, t1, t2, c, alpha )

!*****************************************************************************80
!
!! FD1D_WAVE_ALPHA computes ALPHA for the 1D wave equation.
!
!  Discussion:
!
!    The explicit timestepping procedure uses the quantity ALPHA, which
!    is determined by this function.
!
!    If the spatial region bounds are X1 <= X <= X2, containing X_NUM equally
!    spaced nodes, including the endpoints, and the time domain similarly
!    extends from T1 <= T <= T2 containing T_NUM equally spaced time values,
!    then
!
!      ALPHA = C * DT / DX
!            = C * ( ( T2 - T1 ) / ( T_NUM - 1 ) )
!                / ( ( X2 - X1 ) / ( X_NUM - 1 ) ).
!
!    For a stable computation, it must be the case that ALPHA < 1.
!
!    If ALPHA is greater than 1, then the middle coefficient 1-C^2 DT^2 / DX^2 
!    is negative, and the sum of the magnitudes of the three coefficients 
!    becomes unbounded.  In such a case, the user must reduce the time step 
!    size appropriately.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    George Lindfield, John Penny,
!    Numerical Methods Using MATLAB,
!    Second Edition,
!    Prentice Hall, 1999,
!    ISBN: 0-13-012641-1,
!    LC: QA297.P45.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) X_NUM, the number of nodes in the X direction.
!
!    Input, real ( kind = 8 ) X1, X2, the first and last X coordinates.
!
!    Input, integer ( kind = 4 ) T_NUM, the number of time steps, including the 
!    initial condition.
!
!    Input, real ( kind = 8 ) T1, T2, the first and last T coordinates.
!
!    Input, real ( kind = 8 ) C, a parameter which gives the speed of waves.
!
!    Output, real ( kind = 8 ) ALPHA, the stability coefficient.
!
  implicit none

  real ( kind = 8 ) alpha
  real ( kind = 8 ) c
  real ( kind = 8 ) t_delta
  integer ( kind = 4 ) t_num
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) x_delta
  integer ( kind = 4 ) x_num
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2

  t_delta = ( t2 - t1 ) / real ( t_num - 1, kind = 8 )
  x_delta = ( x2 - x1 ) / real ( x_num - 1, kind = 8 )
  alpha = c * t_delta / x_delta

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Stability condition ALPHA = C * DT / DX = ', alpha

  if ( 1.0D+00 < abs ( alpha ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FD1D_WAVE_ALPHA - Warning!'
    write ( *, '(a)' ) '  The stability condition |ALPHA| <= 1 fails.'
    write ( *, '(a)' ) '  Computed results are liable to be inaccurate.'
  end if

  return
end
subroutine fd1d_wave_start ( x_num, x_vec, t, t_delta, alpha, u_x1, u_x2, &
  ut_t1, u1, u2 )

!*****************************************************************************80
!
!! FD1D_WAVE_START takes the first step for the wave equation.
!
!  Discussion:
!
!    This program solves the 1D wave equation of the form:
!
!      Utt = c^2 Uxx
!
!    over the spatial interval [X1,X2] and time interval [T1,T2],
!    with initial conditions:
!
!      U(T1,X)  = U_T1(X),
!      Ut(T1,X) = UT_T1(X),
!
!    and boundary conditions of Dirichlet type:
!
!      U(T,X1) = U_X1(T),
!      U(T,X2) = U_X2(T).
!
!    The value C represents the propagation speed of waves.
!
!    The program uses the finite difference method, and marches
!    forward in time, solving for all the values of U at the next
!    time step by using the values known at the previous two time steps.
!
!    Central differences may be used to approximate both the time
!    and space derivatives in the original differential equation.
!
!    Thus, assuming we have available the approximated values of U
!    at the current and previous times, we may write a discretized
!    version of the wave equation as follows:
!
!      Uxx(T,X) = ( U(T,   X+dX) - 2 U(T,X) + U(T,   X-dX) ) / dX^2
!      Utt(T,X) = ( U(T+dt,X   ) - 2 U(T,X) + U(T-dt,X   ) ) / dT^2
!
!    If we multiply the first term by C^2 and solve for the single
!    unknown value U(T+dt,X), we have:
!
!      U(T+dT,X) =        (     C^2 * dT^2 / dX^2 ) * U(T,   X+dX)
!                  +  2 * ( 1 - C^2 * dT^2 / dX^2 ) * U(T,   X   )
!                  +      (     C^2 * dT^2 / dX^2 ) * U(T,   X-dX)
!                  -                                  U(T-dT,X   )
!
!    (Equation to advance from time T to time T+dT, except for FIRST step!)
!
!    However, on the very first step, we only have the values of U
!    for the initial time, but not for the previous time step.
!    In that case, we use the initial condition information for dUdT
!    which can be approximated by a central difference that involves
!    U(T+dT,X) and U(T-dT,X):
!
!      dU/dT(T,X) = ( U(T+dT,X) - U(T-dT,X) ) / ( 2 * dT )
!
!    and so we can estimate U(T-dT,X) as
!
!      U(T-dT,X) = U(T+dT,X) - 2 * dT * dU/dT(T,X)
!
!    If we replace the "missing" value of U(T-dT,X) by the known values
!    on the right hand side, we now have U(T+dT,X) on both sides of the
!    equation, so we have to rearrange to get the formula we use
!    for just the first time step:
!
!      U(T+dT,X) =   1/2 * (     C^2 * dT^2 / dX^2 ) * U(T,   X+dX)
!                  +       ( 1 - C^2 * dT^2 / dX^2 ) * U(T,   X   )
!                  + 1/2 * (     C^2 * dT^2 / dX^2 ) * U(T,   X-dX)
!                  +  dT *                         dU/dT(T,   X   )
!
!    (Equation to advance from time T to time T+dT for FIRST step.)
!
!    It should be clear now that the quantity ALPHA = C * DT / DX will affect
!    the stability of the calculation.  If it is greater than 1, then
!    the middle coefficient 1-C^2 DT^2 / DX^2 is negative, and the
!    sum of the magnitudes of the three coefficients becomes unbounded.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    George Lindfield, John Penny,
!    Numerical Methods Using MATLAB,
!    Second Edition,
!    Prentice Hall, 1999,
!    ISBN: 0-13-012641-1,
!    LC: QA297.P45.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) X_NUM, the number of nodes in the X direction.
!
!    Input, real ( kind = 8 ) X_VEC(X_NUM), the coordinates of the nodes.
!
!    Input, real ( kind = 8 ) T, the time after the first step has been taken.
!    In other words, T = T1 + T_DELTA.
!
!    Input, real ( kind = 8 ) T_DELTA, the time step.
!
!    Input, real ( kind = 8 ) ALPHA, the stability coefficient, computed 
!    by FD1D_WAVE_ALPHA.
!
!    Input, real ( kind = 8 ) U_X1(T), U_X2(T), functions for the left and 
!    right boundary conditions.
!
!    Input, real ( kind = 8 ) UT_T1(X), the function that evaluates dUdT at the 
!    initial time.
!
!    Input, real U1(X_NUM), the initial condition.
!
!    Output, real U2(X_NUM), the solution at the first time step.
!
  implicit none

  integer ( kind = 4 ) x_num

  real ( kind = 8 ) alpha
  real ( kind = 8 ) t
  real ( kind = 8 ) t_delta
  external u_x1
  external u_x2
  external ut_t1
  real ( kind = 8 ) u1(x_num)
  real ( kind = 8 ) u2(x_num)
  real ( kind = 8 ) ut(x_num)
  real ( kind = 8 ) x_vec(x_num)

  call ut_t1 ( x_num, x_vec, ut )

  call u_x1 ( t, u2(1) )

  u2(2:x_num-1) =               alpha**2   * u1(3:x_num) / 2.0D+00 &
                  + ( 1.0D+00 - alpha**2 ) * u1(2:x_num-1) &
                  +             alpha**2   * u1(1:x_num-2) / 2.0D+00 &
                  +             t_delta    * ut(2:x_num-1)

  call u_x2 ( t, u2(x_num) )

  return
end
subroutine fd1d_wave_step ( x_num, t, alpha, u_x1, u_x2, u1, u2, u3 )

!*****************************************************************************80
!
!! FD1D_WAVE_STEP computes a step of the 1D wave equation.
!
!  Discussion:
!
!    This program solves the 1D wave equation of the form:
!
!      Utt = c^2 Uxx
!
!    over the spatial interval [X1,X2] and time interval [T1,T2],
!    with initial conditions:
!
!      U(T1,X)  = U_T1(X),
!      Ut(T1,X) = UT_T1(X),
!
!    and boundary conditions of Dirichlet type:
!
!      U(T,X1) = U_X1(T),
!      U(T,X2) = U_X2(T).
!
!    The value C represents the propagation speed of waves.
!
!    The program uses the finite difference method, and marches
!    forward in time, solving for all the values of U at the next
!    time step by using the values known at the previous two time steps.
!
!    Central differences may be used to approximate both the time
!    and space derivatives in the original differential equation.
!
!    Thus, assuming we have available the approximated values of U
!    at the current and previous times, we may write a discretized
!    version of the wave equation as follows:
!
!      Uxx(T,X) = ( U(T,   X+dX) - 2 U(T,X) + U(T,   X-dX) ) / dX^2
!      Utt(T,X) = ( U(T+dt,X   ) - 2 U(T,X) + U(T-dt,X   ) ) / dT^2
!
!    If we multiply the first term by C^2 and solve for the single
!    unknown value U(T+dt,X), we have:
!
!      U(T+dT,X) =        (     C^2 * dT^2 / dX^2 ) * U(T,   X+dX)
!                  +  2 * ( 1 - C^2 * dT^2 / dX^2 ) * U(T,   X   )
!                  +      (     C^2 * dT^2 / dX^2 ) * U(T,   X-dX)
!                  -                                  U(T-dT,X   )
!
!    (Equation to advance from time T to time T+dT, except for FIRST step!)
!
!    However, on the very first step, we only have the values of U
!    for the initial time, but not for the previous time step.
!    In that case, we use the initial condition information for dUdT
!    which can be approximated by a central difference that involves
!    U(T+dT,X) and U(T-dT,X):
!
!      dU/dT(T,X) = ( U(T+dT,X) - U(T-dT,X) ) / ( 2 * dT )
!
!    and so we can estimate U(T-dT,X) as
!
!      U(T-dT,X) = U(T+dT,X) - 2 * dT * dU/dT(T,X)
!
!    If we replace the "missing" value of U(T-dT,X) by the known values
!    on the right hand side, we now have U(T+dT,X) on both sides of the
!    equation, so we have to rearrange to get the formula we use
!    for just the first time step:
!
!      U(T+dT,X) =   1/2 * (     C^2 * dT^2 / dX^2 ) * U(T,   X+dX)
!                  +       ( 1 - C^2 * dT^2 / dX^2 ) * U(T,   X   )
!                  + 1/2 * (     C^2 * dT^2 / dX^2 ) * U(T,   X-dX)
!                  +  dT *                         dU/dT(T,   X   )
!
!    (Equation to advance from time T to time T+dT for FIRST step.)
!
!    It should be clear now that the quantity ALPHA = C * DT / DX will affect
!    the stability of the calculation.  If it is greater than 1, then
!    the middle coefficient 1-C^2 DT^2 / DX^2 is negative, and the
!    sum of the magnitudes of the three coefficients becomes unbounded.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    George Lindfield, John Penny,
!    Numerical Methods Using MATLAB,
!    Second Edition,
!    Prentice Hall, 1999,
!    ISBN: 0-13-012641-1,
!    LC: QA297.P45.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) X_NUM, the number of nodes in the X direction.
!
!    Input, real ( kind = 8 ) T, the new time, that is, the current 
!    time + T_DELTA.
!
!    Input, real ( kind = 8 ) ALPHA, the stability coefficient, computed 
!    by FD1D_WAVE_ALPHA.
!
!    Input, real ( kind = 8 ) U_X1(T), U_X2(T), functions for the left and 
!    right boundary conditions.
!
!    Input, real ( kind = 8 ) U1(X_NUM), the solution at the old time.
!
!    Input, real ( kind = 8 ) U2(X_NUM), the solution at the current time.
!
!    Output, real ( kind = 8 ) U3(X_NUM), the solution at the new time.
!
  implicit none

  integer ( kind = 4 ) x_num

  real ( kind = 8 ) alpha
  real ( kind = 8 ) t
  external u_x1
  external u_x2
  real ( kind = 8 ) u1(x_num)
  real ( kind = 8 ) u2(x_num)
  real ( kind = 8 ) u3(x_num)

  call u_x1 ( t, u3(1) )

  u3(2:x_num-1) =                         alpha**2   * u2(3:x_num) &
                  + 2.0D+00 * ( 1.0D+00 - alpha**2 ) * u2(2:x_num-1) &
                  +                       alpha**2   * u2(1:x_num-2) &
                  -                                    u1(2:x_num-1)

  call u_x2 ( t, u3(x_num) )


  return
end
subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is a value between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is a value between 1 and 99, representing a
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
  logical lopen

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
subroutine piecewise_linear ( nd, xd, yd, nv, xv, yv )

!*****************************************************************************80
!
!! PIECEWISE_LINEAR evaluates a piecewise linear spline.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ND, the number of data points.
!
!    Input, real ( kind = 8 ) XD(ND), YD(ND), the data values.
!
!    Input, integer ( kind = 4 ) NV, the number of evaluation points.
!
!    Input, real ( kind = 8 ) XV(NV), the evaluation arguments.
!
!    Output, real ( kind = 8 ) YV(NV), the values.
!
  implicit none

  integer ( kind = 4 ) nd
  integer ( kind = 4 ) nv

  integer ( kind = 4 ) id
  integer ( kind = 4 ) iv
  real ( kind = 8 ) xd(nd)
  real ( kind = 8 ) xv(nv)
  real ( kind = 8 ) yd(nd)
  real ( kind = 8 ) yv(nv)

  do iv = 1, nv

    if ( xv(iv) < xd(1) ) then
      yv(iv) = yd(1)
    else if ( xd(nd) < xv(iv) ) then
      yv(iv) = yd(nd)
    else 

      do id = 2, nd
        if ( xv(iv) < xd(id) ) then
          yv(iv) = ( ( xd(id) - xv(iv)            ) * yd(id-1) &
                   + (          xv(iv) - xd(id-1) ) * yd(id) ) &
                   / ( xd(id)          - xd(id-1) )
          exit
        end if
      end do

    end if

  end do

  return
end
subroutine r8mat_write ( output_filename, m, n, table )

!*****************************************************************************80
!
!! R8MAT_WRITE writes an R8MAT file.
!
!  Discussion:
!
!    An R8MAT is an array of R8 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) OUTPUT_FILENAME, the output file name.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) TABLE(M,N), the data.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) j
  character ( len = * ) output_filename
  integer ( kind = 4 ) output_status
  integer ( kind = 4 ) output_unit
  character ( len = 30 ) string
  real ( kind = 8 ) table(m,n)
!
!  Open the file.
!
  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_filename, &
    status = 'replace', iostat = output_status )

  if ( output_status /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_WRITE - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the output file "' // &
      trim ( output_filename ) // '" on unit ', output_unit
    output_unit = -1
    stop
  end if
!
!  Create a format string.
!
!  For less precision in the output file, try:
!
!                                            '(', m, 'g', 14, '.', 6, ')'
!
  if ( 0 < m .and. 0 < n ) then

    write ( string, '(a1,i8,a1,i8,a1,i8,a1)' ) '(', m, 'g', 24, '.', 16, ')'
!
!  Write the data.
!
    do j = 1, n
      write ( output_unit, string ) table(1:m,j)
    end do

  end if
!
!  Close the file.
!
  close ( unit = output_unit )

  return
end
subroutine r8vec_linspace ( n, a_first, a_last, a )

!*****************************************************************************80
!
!! R8VEC_LINSPACE creates a vector of linearly spaced values.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 March 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, real ( kind = 8 ) A_FIRST, A_LAST, the first and last entries.
!
!    Output, real ( kind = 8 ) A(N), a vector of linearly spaced data.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) a_first
  real ( kind = 8 ) a_last
  integer ( kind = 4 ) i

  if ( n == 1 ) then

    a(1) = ( a_first + a_last ) / 2.0D+00

  else

    do i = 1, n
      a(i) = ( real ( n - i,     kind = 8 ) * a_first &
             + real (     i - 1, kind = 8 ) * a_last ) &
             / real ( n     - 1, kind = 8 )
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
