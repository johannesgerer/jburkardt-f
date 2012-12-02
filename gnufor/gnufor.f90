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
!    18 September 2005
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
subroutine run_gnuplot ( command_filename )

!*****************************************************************************80
!
!! RUN_GNUPLOT runs GNUPLOT with a given command file.
!
!  Discussion:
!
!    The GNUPLOT program must be available.  To check whether
!    this is so, try typing
!
!      which gnuplot
!
!    If the response is
!
!      gnuplot: command not found
!
!    then you're going to have to make GNUPLOT available.
!
!    You may need to set the environment variable GNUTERM:
!
!      setenv GNUTERM x11
!
!    so that GNUPLOT automatically displays to your X window terminal.
!
!
!    This routine expects that there is a text file containing the appropriate
!    commands to GNUPLOT to display your picture.  There are a number of
!    routines in this package that will do this for simple plotting tasks.
!    Most of them require that you also set up a file of data to be plotted.
!
!    Once this routine invokes GNUPLOT, a graphics window should open
!    up, and the FORTRAN program will pause.  Hitting RETURN should advance
!    to the next picture, or terminate the window at the end, allowing the
!    FORTRAN routine to proceed.
!
!
!    You can look at the data and command files created by the routines.
!    Moreover, you can easily modify the command file to change the options
!    used in GNUPLOT, and then run GNUPLOT interactively, as in:
!
!      gnuplot commands
!
!    In particular, if you want a PostScript version of your graphics files,
!    insert the command "set term postscript" at the beginning of the command
!    file and run gnuplot as follows:
!
!      gnuplot commands > mypicture.ps
!
!    You will also have to hit RETURN once for each plot that is made.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) COMMAND_FILENAME, the name of the
!    command file.
!
  implicit none

  character ( len = 255 ) command
  character ( len = * ) command_filename
  integer ( kind = 4 ) status
  integer ( kind = 4 ) system

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GNUFOR:'
  write ( *, '(a)' ) '  GNUPLOT / FORTRAN90 command interface.'
!
!  Issue a command to the system that will start GNUPLOT, using
!  the file we just wrote as input.
!
!  The "&" will run GNUPLOT in the background, so the FORTRAN program
!  can continue execution independently, and the PERSIST switch tells
!  GNUPLOT that if there are multiple plots, they should each go in
!  a separate window.
!
!  Thanks to Morag Am-Shallem for suggesting these improvements.
!  17 October 2007
!
  write ( command, * ) 'gnuplot -persist ' // trim ( command_filename ) // ' &'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GNUFOR:'
  write ( *, '(a)' ) '  Issuing the command:"' // trim ( command ) // '".'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Press RETURN to proceed.'

  status = system ( trim ( command ) )

  if ( status /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GNUFOR - Fatal error!'
    write ( *, '(a)' ) '  An error code was returned when the GNUPLOT command'
    write ( *, '(a)' ) '  was issued.  Perhaps GNUPLOT is not in your path.'
    write ( *, '(a)' ) '  Type "which gnuplot" to check this.'
    stop
  end if
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GNUFOR:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

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
subroutine write_vector_data ( data_filename, n, x, y, dx, dy, ierror )

!*****************************************************************************80
!
!! WRITE_VECTOR_DATA writes vector data to a file, for plotting by GNUPLOT.
!
!  Discussion:
!
!    Each vector is described by 4 values, X, Y, dX, dY, indicating that
!    a vector is to be drawn from (X,Y) to (X+dX,Y+dY).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) DATA_FILENAME, the name of the data file.
!
!    Input, integer ( kind = 4 ) N, the number of vectors.
!
!    Input, real ( kind = 8 ) X(N), Y(N), DX(N), DY(N), the vector data.
!
!    Output, integer ( kind = 4 ) IERROR, nonzero if an error occurred.
!
  implicit none

  integer n

  character ( len = * ) data_filename
  real ( kind = 8 ) dx(n)
  real ( kind = 8 ) dy(n)
  integer ( kind = 4 ) file_unit
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)

  ierror = 0

  call get_unit ( file_unit )

  if ( file_unit == 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WRITE_VECTOR_DATA - Fatal error!'
    write ( *, '(a)' ) '  Could not get a free FORTRAN unit.'
    return
  end if

  open ( unit = file_unit, file = data_filename, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    ierror = 2
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WRITE_VECTOR_DATA - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file.'
    return
  end if

  do i = 1, n
    write ( file_unit, * ) x(i), y(i), dx(i), dy(i)
  end do

  close ( unit = file_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'WRITE_VECTOR_DATA:'
  write ( *, '(a)' ) '  Wrote the GNUPLOT vector data file "' // &
    trim ( data_filename ) // '"'

  return
end
subroutine write_vector_plot ( command_filename, data_filename, &
  ierror )

!*****************************************************************************80
!
!! WRITE_VECTOR_PLOT writes GNUPLOT commands to plot vectors.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 June 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) COMMAND_FILENAME, the name of the
!    command file.
!
!    Input, character ( len = * ) DATA_FILENAME, the name of the data file.
!
!    Output, integer ( kind = 4 ) IERROR, nonzero if an error occurred.
!
  implicit none

  character ( len = * ) command_filename
  character ( len = * ) data_filename
  integer ( kind = 4 ) file_unit
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
!
!  Write the data file.
!
  ierror = 0

  call get_unit ( file_unit )

  if ( file_unit == 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WRITE_VECTOR_PLOT - Fatal error!'
    write ( *, '(a)' ) '  Could not get a free FORTRAN unit.'
    return
  end if

  open ( unit = file_unit, file = command_filename, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    ierror = 2
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WRITE_VECTOR_PLOT - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file.'
    return
  end if

  write ( file_unit, '(a)' ) 'set title "GNUFOR plot"'
  write ( file_unit, '(a)' ) 'set xlabel "x"'
  write ( file_unit, '(a)' ) 'set ylabel "y"'
  write ( file_unit, '(a)' ) 'set style data vector'
  write ( file_unit, '(a,i2,a)' ) 'plot "' // trim ( data_filename )
  write ( file_unit, '(a)' ) 'pause -1'
  write ( file_unit, '(a)' ) 'q'

  close ( unit = file_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'WRITE_VECTOR_PLOT:'
  write ( *, '(a)' ) '  Wrote the GNUPLOT table plots command file "' // &
    trim ( command_filename ) // '"'

  return
end
subroutine write_xy_data ( data_filename, n, x, y, ierror )

!*****************************************************************************80
!
!! WRITE_XY_DATA writes X(1:N), Y(1:N) data to a file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) DATA_FILENAME, the name of the data file.
!
!    Input, integer ( kind = 4 ) N, the number of data items.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the X and Y data
!
!    Output, integer ( kind = 4 ) IERROR, nonzero if an error occurred.
!
  implicit none

  integer ( kind = 4 ) n

  character ( len = * ) data_filename
  integer ( kind = 4 ) file_unit
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)

  ierror = 0

  call get_unit ( file_unit )

  if ( file_unit == 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WRITE_XY_DATA - Fatal error!'
    write ( *, '(a)' ) '  Could not get a free FORTRAN unit.'
    return
  end if

  open ( unit = file_unit, file = data_filename, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    ierror = 2
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WRITE_XY_DATA - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file.'
    return
  end if

  do i = 1, n
    write ( file_unit, * ) x(i), y(i)
  end do

  close ( unit = file_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'WRITE_XY_DATA:'
  write ( *, '(a)' ) '  Wrote the GNUPLOT XY data file "' // &
    trim ( data_filename ) // '"'

  return
end
subroutine write_xy_plot ( command_filename, data_filename, ierror )

!*****************************************************************************80
!
!! WRITE_XY_PLOT writes GNUPLOT commands to plot X(1:N), Y(1:N) data.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 June 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) COMMAND_FILENAME, the name of the
!    command file.
!
!    Input, character ( len = * ) DATA_FILENAME, the name of the data file.
!
!    Output, integer ( kind = 4 ) IERROR, nonzero if an error occurred.
!
  implicit none

  character ( len = * ) command_filename
  character ( len = * ) data_filename
  integer ( kind = 4 ) file_unit
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
!
!  Write the data file.
!
  ierror = 0

  call get_unit ( file_unit )

  if ( file_unit == 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WRITE_XY_PLOT - Fatal error!'
    write ( *, '(a)' ) '  Could not get a free FORTRAN unit.'
    return
  end if

  open ( unit = file_unit, file = command_filename, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    ierror = 2
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WRITE_XY_PLOT - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file.'
    return
  end if

  write ( file_unit, '(a)' ) 'set title "GNUFOR plot"'
  write ( file_unit, '(a)' ) 'set xlabel "x"'
  write ( file_unit, '(a)' ) 'set ylabel "y"'
  write ( file_unit, '(a,i2,a)' ) 'plot "' // trim ( data_filename ) // &
    '" using 1:2 with lines'
  write ( file_unit, '(a)' ) 'pause -1'
  write ( file_unit, '(a)' ) 'q'

  close ( unit = file_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'WRITE_XY_PLOT:'
  write ( *, '(a)' ) '  Wrote the GNUPLOT XY plot command file "' // &
    trim ( command_filename ) // '"'

  return
end
subroutine write_xyy_data ( data_filename, lda, nrow, ncol, x, ierror )

!*****************************************************************************80
!
!! WRITE_XYY_DATA writes a table of data to a file, for plotting by GNUPLOT.
!
!  Discussion:
!
!    The first column of data is assumed to be the independent variable, X.
!    Separate plots are made of X versus all the other columns of data.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) DATA_FILENAME, the name of the data file.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of X.
!
!    Input, integer ( kind = 4 ) NROW, NCOL, the dimensions of X.
!
!    Input, real ( kind = 8 ) X(LDA,NCOL), the NROW by NCOL data to be written.
!
!    Output, integer ( kind = 4 ) IERROR, nonzero if an error occurred.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) ncol

  character ( len = * ) data_filename
  integer ( kind = 4 ) file_unit
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) nrow
  real ( kind = 8 ) x(lda,ncol)

  ierror = 0

  call get_unit ( file_unit )

  if ( file_unit == 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WRITE_XYY_DATA - Fatal error!'
    write ( *, '(a)' ) '  Could not get a free FORTRAN unit.'
    return
  end if

  open ( unit = file_unit, file = data_filename, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    ierror = 2
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WRITE_XYY_DATA - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file.'
    return
  end if

  do i = 1, nrow
    write ( file_unit, * ) x(i,1:ncol)
  end do

  close ( unit = file_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'WRITE_XYY_DATA:'
  write ( *, '(a)' ) '  Wrote the GNUPLOT XYY data file "' // &
    trim ( data_filename ) // '"'

  return
end
subroutine write_xyy_plots ( command_filename, data_filename, &
  ncol, ierror )

!*****************************************************************************80
!
!! WRITE_XYY_PLOTS writes GNUPLOT commands to make multiple (X,Y) plots.
!
!  Discussion:
!
!    The first column of data is assumed to be the independent variable, X.
!    Separate plots are made of X versus all the other columns of data.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 June 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) COMMAND_FILENAME, the name of the
!    command file.
!
!    Input, character ( len = * ) DATA_FILENAME, the name of the data file.
!
!    Input, integer ( kind = 4 ) NCOL, the number of columns of data.
!
!    Output, integer ( kind = 4 ) IERROR, nonzero if an error occurred.
!
  implicit none

  character ( len = * ) command_filename
  character ( len = * ) data_filename
  integer ( kind = 4 ) file_unit
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) ncol
!
!  Write the data file.
!
  ierror = 0

  call get_unit ( file_unit )

  if ( file_unit == 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WRITE_XYY_PLOTS - Fatal error!'
    write ( *, '(a)' ) '  Could not get a free FORTRAN unit.'
    return
  end if

  open ( unit = file_unit, file = command_filename, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    ierror = 2
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WRITE_XYY_PLOTS - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file.'
    return
  end if

  write ( file_unit, '(a)' ) 'set title "GNUFOR plot"'
  write ( file_unit, '(a)' ) 'set xlabel "x"'
  write ( file_unit, '(a)' ) 'set ylabel "y"'
  do i = 2, ncol
    write ( file_unit, '(a,i2,a)' ) 'plot "' // trim ( data_filename ) // &
      '" using ', i, ' with lines'
    write ( file_unit, '(a)' ) 'pause -1'
  end do
  write ( file_unit, '(a)' ) 'q'

  close ( unit = file_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'WRITE_XYY_PLOTS:'
  write ( *, '(a)' ) '  Wrote the GNUPLOT XYY plots command file "' // &
    trim ( command_filename ) // '"'

  return
end
subroutine write_xyz_data ( data_filename, n, x, y, z, ierror )

!*****************************************************************************80
!
!! WRITE_XYZ_DATA writes X(1:N), Y(1:N), Z(1:N) data to a file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) DATA_FILENAME, the name of the data file.
!
!    Input, integer ( kind = 4 ) N, the number of data items.
!
!    Input, real ( kind = 8 ) X(N), Y(N), Z(N), the X, Y, Z data
!
!    Output, integer ( kind = 4 ) IERROR, nonzero if an error occurred.
!
  implicit none

  integer ( kind = 4 ) n

  character ( len = * ) data_filename
  integer ( kind = 4 ) file_unit
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) z(n)

  ierror = 0

  call get_unit ( file_unit )

  if ( file_unit == 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WRITE_XYZ_DATA - Fatal error!'
    write ( *, '(a)' ) '  Could not get a free FORTRAN unit.'
    return
  end if

  open ( unit = file_unit, file = data_filename, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    ierror = 2
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WRITE_XYZ_DATA - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file.'
    return
  end if

  do i = 1, n
    write ( file_unit, * ) x(i), y(i), z(i)
  end do

  close ( unit = file_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'WRITE_XYZ_DATA:'
  write ( *, '(a)' ) '  Wrote the GNUPLOT XYZ data file "' // &
    trim ( data_filename ) // '"'

  return
end
subroutine write_xyz_plot ( command_filename, data_filename, ierror )

!*****************************************************************************80
!
!! WRITE_XYZ_PLOT writes commands to plot parametric (X,Y,Z) data.
!
!  Discussion:
!
!    This routine tries to write a command file suitable for displaying
!    a 3D arc specified by points (X,Y,Z).  A grid data file, containing
!    values of X, Y and Z, will also be needed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) COMMAND_FILENAME, the name of the
!    command file.
!
!    Input, character ( len = * ) DATA_FILENAME, the name of the data file.
!
!    Output, integer ( kind = 4 ) IERROR, nonzero if an error occurred.
!
  implicit none

  character ( len = * ) command_filename
  character ( len = * ) data_filename
  integer ( kind = 4 ) file_unit
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) ncol
!
!  Write the data file.
!
  ierror = 0

  call get_unit ( file_unit )

  if ( file_unit == 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WRITE_XYZ_PLOT - Fatal error!'
    write ( *, '(a)' ) '  Could not get a free FORTRAN unit.'
    return
  end if

  open ( unit = file_unit, file = command_filename, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    ierror = 2
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WRITE_XYZ_PLOT - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file.'
    return
  end if

  write ( file_unit, '(a)' ) 'set title "GNUFOR plot"'
  write ( file_unit, '(a)' ) 'set xlabel "x"'
  write ( file_unit, '(a)' ) 'set ylabel "y"'
  write ( file_unit, '(a)' ) 'set parametric'
  write ( file_unit, '(a)' ) 'splot "' // trim ( data_filename ) // &
      '" using 1:2:3 with lines'
  write ( file_unit, '(a)' ) 'pause -1'
  write ( file_unit, '(a)' ) 'q'

  close ( unit = file_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'WRITE_XYZ_PLOT:'
  write ( *, '(a)' ) '  Wrote the GNUPLOT SPLOT command file "' // &
    trim ( command_filename ) // '"'

  return
end
subroutine write_xyzgrid_contour ( command_filename, data_filename, ierror )

!*****************************************************************************80
!
!! WRITE_XYZGRID_CONTOUR writes commands to plot contours of Z(X,Y).
!
!  Discussion:
!
!    This routine tries to write a command file suitable for displaying
!    contours of Z(X,Y) gridded data.  A grid data file, containing values
!    of X, Y and Z, will also be needed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 June 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) COMMAND_FILENAME, the name of the
!    command file.
!
!    Input, character ( len = * ) DATA_FILENAME, the name of the data file.
!
!    Output, integer ( kind = 4 ) IERROR, nonzero if an error occurred.
!
  implicit none

  character ( len = * ) command_filename
  character ( len = * ) data_filename
  integer ( kind = 4 ) file_unit
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) ncol
!
!  Write the data file.
!
  ierror = 0

  call get_unit ( file_unit )

  if ( file_unit == 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WRITE_XYZGRID_CONTOUR - Fatal error!'
    write ( *, '(a)' ) '  Could not get a free FORTRAN unit.'
    return
  end if

  open ( unit = file_unit, file = command_filename, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    ierror = 2
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WRITE_XYZGRID_CONTOUR - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file.'
    return
  end if

  write ( file_unit, '(a)' ) 'set title "GNUFOR plot"'
  write ( file_unit, '(a)' ) 'set xlabel "x"'
  write ( file_unit, '(a)' ) 'set ylabel "y"'
  write ( file_unit, '(a)' ) 'set parametric'
  write ( file_unit, '(a)' ) 'set nosurface'
  write ( file_unit, '(a)' ) 'set contour'
  write ( file_unit, '(a)' ) 'set cntrparam levels 10'
  write ( file_unit, '(a)' ) 'splot "' // trim ( data_filename ) // &
    '" using 1:2:3 with lines'

  write ( file_unit, '(a)' ) 'set table "table.txt"'
  write ( file_unit, '(a)' ) 'replot'
  write ( file_unit, '(a)' ) 'unset table'

  write ( file_unit, '(a)' ) 'set term x11'
  write ( file_unit, '(a)' ) 'plot "table.txt" with lines'
  write ( file_unit, '(a)' ) 'pause -1'
  write ( file_unit, '(a)' ) 'q'

  close ( unit = file_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'WRITE_XYZGRID_CONTOUR:'
  write ( *, '(a)' ) &
    '  Wrote the GNUPLOT XYZGRID contour plot command file "' // &
    trim ( command_filename ) // '"'

  return
end
subroutine write_xyzgrid_data ( data_filename, nx, ny, xyz, ierror )

!*****************************************************************************80
!
!! WRITE_XYZGRID_DATA writes a file of XYZ grid data.
!
!  Discussion:
!
!    It is assumed that values of Z are available on a regular NX by NY grid
!    of (X,Y) points.
!
!    The form of the data file requires that all the data for a given value
!    of Y be listed, followed by a blank line, followed by the data for
!    another value of Y.
!
!  Example:
!
!    Here is a grid data file for a 3 by 3 grid, with Z = X + Y.
!
!    0.0 0.0 0.0
!    1.0 0.0 1.0
!    2.0 0.0 2.0
!
!    0.0 1.0 1.0
!    1.0 1.0 2.0
!    2.0 1.0 3.0
!
!    0.0 2.0 2.0
!    1.0 2.0 3.0
!    2.0 2.0 4.0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) DATA_FILENAME, the name of the data file.
!
!    Input, integer ( kind = 4 ) NX, NY, the dimensions of the grid.
!
!    Input, real ( kind = 8 ) XYZ(3,NX,NY), the XYZ grid data to be written.
!
!    Output, integer ( kind = 4 ) IERROR, nonzero if an error occurred.
!
  implicit none

  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny

  character ( len = * ) data_filename
  integer ( kind = 4 ) file_unit
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) j
  real ( kind = 8 ) xyz(3,nx,ny)

  ierror = 0

  call get_unit ( file_unit )

  if ( file_unit == 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WRITE_XYZGRID_DATA - Fatal error!'
    write ( *, '(a)' ) '  Could not get a free FORTRAN unit.'
    return
  end if

  open ( unit = file_unit, file = data_filename, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    ierror = 2
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WRITE_XYZGRID_DATA - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file.'
    return
  end if

  do j = 1, ny
    do i = 1, nx
      write ( file_unit, * ) xyz(1:3,i,j)
    end do
    write ( file_unit, '(a)' )
  end do

  close ( unit = file_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'WRITE_XYZGRID_DATA:'
  write ( *, '(a)' ) '  Wrote the GNUPLOT XYZ grid data file "' // &
    trim ( data_filename ) // '"'

  return
end
subroutine write_xyzgrid_surface ( command_filename, data_filename, ierror )

!*****************************************************************************80
!
!! WRITE_XYZGRID_SURFACE writes a file of GNUPLOT commands to plot a 3D surface.
!
!  Discussion:
!
!    This routine tries to write a command file suitable for displaying
!    a surface Z(X,Y).  A grid data file, containing values of X, Y and Z,
!    will also be needed.  The routine WRITE_XYZGRID_DATA can write this file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) COMMAND_FILENAME, the name of the
!    command file.
!
!    Input, character ( len = * ) DATA_FILENAME, the name of the data file.
!
!    Output, integer ( kind = 4 ) IERROR, nonzero if an error occurred.
!
  implicit none

  character ( len = * ) command_filename
  character ( len = * ) data_filename
  integer ( kind = 4 ) file_unit
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) ncol
!
!  Write the data file.
!
  ierror = 0

  call get_unit ( file_unit )

  if ( file_unit == 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WRITE_XYZGRID_SURFACE - Fatal error!'
    write ( *, '(a)' ) '  Could not get a free FORTRAN unit.'
    return
  end if

  open ( unit = file_unit, file = command_filename, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    ierror = 2
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WRITE_XYZGRID_SURFACE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file.'
    return
  end if

  write ( file_unit, '(a)' ) 'set title "GNUFOR plot"'
  write ( file_unit, '(a)' ) 'set xlabel "x"'
  write ( file_unit, '(a)' ) 'set ylabel "y"'
  write ( file_unit, '(a)' ) 'set parametric'
  write ( file_unit, '(a)' ) 'set hidden3d'
  write ( file_unit, '(a)' ) 'set contour'
  write ( file_unit, '(a)' ) 'splot "' // trim ( data_filename ) // &
    '" using 1:2:3 with lines'
  write ( file_unit, '(a)' ) 'pause -1'
  write ( file_unit, '(a)' ) 'q'

  close ( unit = file_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'WRITE_SURFACE_COMMANDS:'
  write ( *, '(a)' ) '  Wrote the GNUPLOT surface plot command file "' // &
    trim ( command_filename ) // '"'

  return
end
