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
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IUNIT.
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5 and 6).
!    Otherwise, IUNIT is an integer between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  logical ( kind = 4 ) lopen

  iunit = 0

  do i = 1, 99

    if ( i /= 5 .and. i /= 6 ) then

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
subroutine lambert1 ( n, eta )

!*****************************************************************************80
!
!! LAMBERT1 computes the Lambert sequence in 1D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    J P Lambert,
!    Quasi-Random Sequences for Optimization and Numerical Integration,
!    in Numerical Integration,
!    edited by P Keast and G Fairweather,
!    D Reidel, 1987, pages 193-203.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of the sequence
!    to compute.
!
!    Output, real ( kind = 8 ) ETA(1,N), the elements of the sequence.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) eta(1,n)
  integer ( kind = 4 ) j
  real ( kind = 8 ) t
  real ( kind = 8 ) x

  eta(1,1) = 0.0D+00

  do j = 2, n

    t = 1.0D+00

    do

      t = t / 2.0D+00

      if ( x + t < 1.0D+00 ) then
        exit
      end if

    end do

    x = x + t - 1.0D+00

    if ( x < 0.0D+00 ) then
      x = x + 2.0D+00 * t
    end if

    eta(1,j) = x

  end do

  return
end
subroutine lambert2 ( n, eta )

!*****************************************************************************80
!
!! LAMBERT2 computes the Lambert sequence in 2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    J P Lambert,
!    Quasi-Random Sequences for Optimization and Numerical Integration,
!    in Numerical Integration,
!    edited by P Keast and G Fairweather,
!    D Reidel, 1987, pages 193-203.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of the sequence
!    to compute.
!
!    Output, real ( kind = 8 ) ETA(2,N), the elements of the sequence.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) eta(2,n)
  integer ( kind = 4 ) j
  real ( kind = 8 ) t
  real ( kind = 8 ) u
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  eta(1:2,1) = 0.0D+00

  do j = 2, n

    t = 1.0D+00

    do

      t = t / 2.0D+00
      u = 1.0D+00 - t

      if ( x < u .or. t <= y ) then
        exit
      end if

    end do

    x = x - u

    if ( x < 0.0D+00 ) then
      x = x + 2.0D+00 * t
      if ( y < t ) then
        y = y + t
      else
        y = y - t
      end if
    end if

    eta(1:2,j) = (/ x, y /)

  end do

  return
end
subroutine lambert3 ( n, eta )

!*****************************************************************************80
!
!! LAMBERT3 computes the Lambert sequence in 3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    J P Lambert,
!    Quasi-Random Sequences for Optimization and Numerical Integration,
!    in Numerical Integration,
!    edited by P Keast and G Fairweather,
!    D Reidel, 1987, pages 193-203.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of the sequence
!    to compute.
!
!    Output, real ( kind = 8 ) ETA(3,N), the elements of the sequence.
!
  implicit none

  integer n

  real ( kind = 8 ) eta(3,n)
  integer ( kind = 4 ) j
  real ( kind = 8 ) t
  real ( kind = 8 ) u
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  eta(1:3,1) = 0.0D+00

  do j = 2, n

    t = 1.0D+00

    do

      t = t / 2.0D+00
      u = 1.0D+00 - t

      if ( x < u .or. y < u .or. t <= z ) then
        exit
      end if

    end do

    x = x - u
    y = y - u
    z = z - t

    if ( x < 0.0D+00 ) then

      x = x + 2.0D+00 * t
      if ( y < 0.0D+00 ) then
        y = y + 2.0D+00 * t
      end if
      if ( z < 0.0D+00 ) then
        z = z + 2.0D+00 * t
      end if

    else

      if ( y < 0.0D+00 ) then
        if ( z < 0.0D+00 ) then
          y = y + t
          z = z + 2.0D+00 * t
        else if ( y < 0.0D+00 ) then
          y = y + 2.0D+00 * t
          z = z + t
        end if
      else if ( 0.0D+00 <= y ) then
        y = y + t
      end if

    end if

    eta(1:3,j) = (/ x, y, z /)

  end do

  return
end
subroutine lambert4 ( n, eta )

!*****************************************************************************80
!
!! LAMBERT4 computes the Lambert sequence in 4D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    J P Lambert,
!    Quasi-Random Sequences for Optimization and Numerical Integration,
!    in Numerical Integration,
!    edited by P Keast and G Fairweather,
!    D Reidel, 1987, pages 193-203.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of the sequence
!    to compute.
!
!    Output, real ( kind = 8 ) ETA(4,N), the elements of the sequence.
!
  implicit none

  integer n

  real ( kind = 8 ) eta(4,n)
  integer ( kind = 4 ) j
  real ( kind = 8 ) t
  real ( kind = 8 ) u
  real ( kind = 8 ) w
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  eta(1:4,1) = 0.0D+00

  do j = 2, n

    t = 1.0D+00

    do

      t = t / 2.0D+00
      u = 1.0D+00 - t

      if ( x < u .or. y < u .or. z < u .or. t <= w ) then
        exit
      end if

    end do

    x = x - u
    y = y - u
    z = z - u
    w = w - t

    if ( x < 0.0D+00 ) then

      x = x + 2.0D+00 * t

      if ( y < 0.0D+00 ) then
        y = y + 2.0D+00 * t
      end if
      if ( z < 0.0D+00 ) then
        z = z + 2.0D+00 * t
      end if
      if ( w < 0.0D+00 ) then
        w = w + 2.0D+00 * t
      end if

    else if ( y < 0.0D+00 ) then

      if ( z < 0.0D+00 ) then

        if ( w < 0.0D+00 ) then
          y = y + 2.0D+00 * t
          z = z + 2.0D+00 * t
          w = w + t
        else if ( 0.0D+00 <= w ) then
          y = y + t
          z = z + 2.0D+00 * t
        end if

      else if ( 0.0D+00 <= z ) then

        if ( w < 0.0D+00 ) then
          y = y + 2.0D+00 * t
          z = z + t
          w = w + 2.0D+00 * t
        else if ( 0.0D+00 <= w ) then
          y = y + 2.0D+00 * t
          w = w + t
        end if

      end if

    else if ( z < 0.0D+00 ) then

      if ( w < 0.0D+00 ) then
        z = z + t
        w = w + 2.0D+00 * t
      else if ( 0.0D+00 <= w ) then
        z = z + 2.0D+00 * t
        w = w + t
      end if

    else

      y = y + t

    end if

    eta(1:4,j) = (/ x, y, z, w /)

  end do

  return
end
subroutine lambert_write ( m, n, r, file_out_name )

!*****************************************************************************80
!
!! LAMBERT_WRITE writes a Lambert dataset to a file.
!
!  Discussion:
!
!    The initial lines of the file are comments, which begin with a
!    '#' character.
!
!    Thereafter, each line of the file contains the M-dimensional
!    components of the next entry of the Lambert sequence.
!
!  Modified:
!
!    16 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of (successive) points.
!
!    Input, real ( kind = 8 ) R(M,N), the points.
!
!    Input, character ( len = * ) FILE_OUT_NAME, the name of
!    the output file.
!
  implicit none

  integer m
  integer n

  character ( len = * ) file_out_name
  integer ( kind = 4 ) file_out_unit
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) j
  real ( kind = 8 ) r(m,n)

  call get_unit ( file_out_unit )

  open ( unit = file_out_unit, file = file_out_name, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LAMBERT_WRITE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file.'
    stop
  end if

  write ( file_out_unit, '(a)'      ) '#  ' // trim ( file_out_name )
  write ( file_out_unit, '(a)'      ) '#  created by LAMBERT_DATASET.'
  write ( file_out_unit, '(a)'      ) '#'
  write ( file_out_unit, '(a,i6)'   ) '#  Spatial dimension M = ', m
  write ( file_out_unit, '(a,i6)'   ) '#  Number of points N = ', n
  write ( file_out_unit, '(a)'      ) '#'

  do j = 1, n
    write ( file_out_unit, '(20f10.6)' ) r(1:m,j)
  end do

  close ( unit = file_out_unit )

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

  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
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
