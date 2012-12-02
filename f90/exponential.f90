program main

!*****************************************************************************80
!
!! MAIN is the main program for EXPONENTIAL.
!
!  Discussion:
!
!    EXPONENTIAL tries out some ideas for approximating exp(X).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 March 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'EXPONENTIAL:'
  write ( *, '(a)' ) '  FORTRAN90 version.'

  call test01 ( )
  call test02 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'EXPONENTIAL:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '

  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 considers an estimate for the exponential function.
!
!  Discussion:
!
!    EXPONENTIAL approximates exp(X) using a fixed number of terms of
!    a Taylor series.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 March 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real    ( kind = 8 ) factor
  real    ( kind = 8 ) fx
  real    ( kind = 8 ) fx2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ) :: term_num = 5
  real    ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01:'
  write ( *, '(a)' ) '  Try using a Taylor series with fixed number of terms.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '        X             EXP(X)         Approx(X)     Error'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call exp_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = 1.0D+00
    factor = 1.0D+00

    do i = 1, term_num

      factor = 1.0D+00
      do j = 1, i
        factor = factor * j
      end do

      fx2 = fx2 + x**i / factor

    end do

    write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6,2x,g10.4)' ) &
      x, fx, fx2, abs ( fx - fx2 )

  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 considers an estimate for the exponential function.
!
!  Discussion:
!
!    EXPONENTIAL approximates exp(X) using a fixed number of terms of
!    a Taylor series.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 March 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real    ( kind = 8 ) e
  real    ( kind = 8 ) factor
  real    ( kind = 8 ) fx
  real    ( kind = 8 ) fx2
  real    ( kind = 8 ) g
  real    ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ) :: term_num = 5
  real    ( kind = 8 ) x
  real    ( kind = 8 ) x2
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02:'
  write ( *, '(a)' ) '  Try using a Taylor series with fixed number of terms'
  write ( *, '(a)' ) '  and scaling X to be between 0 and 1.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '        X             EXP(X)         Approx(X)     Error'
  write ( *, '(a)' ) ' '

  e = exp ( 1.0D+00 )

  n_data = 0

  do

    call exp_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    x2 = x

    g = 1.0D+00
    do while ( 1.0D+00 < x2 )
      g = g * e
      x2 = x2 - 1.0D+00
    end do

    do while ( x2 < 0.0D+00 )
      g = g / e
      x2 = x2 + 1.0D+00
    end do

    factor = 1.0D+00

    h = 1.0D+00
    do i = 1, term_num

      factor = 1.0D+00
      do j = 1, i
        factor = factor * j
      end do

      h = h + x2**i / factor

    end do

    fx2 = g * h

    write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6,2x,g10.4)' ) &
      x, fx, fx2, abs ( fx - fx2 )

  end do

  return
end
subroutine exp_values ( n_data, x, fx )

!*****************************************************************************80
!
!! EXP_VALUES returns some values of the exponential function.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      Exp[x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 24

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.000045399929762484851536D+00, &
    0.0067379469990854670966D+00, &
    0.36787944117144232160D+00, &
    1.0000000000000000000D+00, &
    1.0000000100000000500D+00, &
    1.0001000050001666708D+00, &
    1.0010005001667083417D+00, &
    1.0100501670841680575D+00, &
    1.1051709180756476248D+00, &
    1.2214027581601698339D+00, &
    1.3498588075760031040D+00, &
    1.4918246976412703178D+00, &
    1.6487212707001281468D+00, &
    1.8221188003905089749D+00, &
    2.0137527074704765216D+00, &
    2.2255409284924676046D+00, &
    2.4596031111569496638D+00, &
    2.7182818284590452354D+00, &
    7.3890560989306502272D+00, &
    23.140692632779269006D+00, &
    148.41315910257660342D+00, &
    22026.465794806716517D+00, &
    4.8516519540979027797D+08, &
    2.3538526683701998541D+17 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
     -10.0D+00, &
      -5.0D+00, &
      -1.0D+00, &
       0.0D+00, &
       0.00000001D+00, &
       0.0001D+00, &
       0.001D+00, &
       0.01D+00, &
       0.1D+00, &
       0.2D+00, &
       0.3D+00, &
       0.4D+00, &
       0.5D+00, &
       0.6D+00, &
       0.7D+00, &
       0.8D+00, &
       0.9D+00, &
       1.0D+00, &
       2.0D+00, &
       3.1415926535897932385D+00, &
       5.0D+00, &
      10.0D+00, &
      20.0D+00, &
      40.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
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
