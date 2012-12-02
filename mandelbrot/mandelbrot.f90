program main

!*****************************************************************************80
!
!! MAIN is the main program for MANDELBROT.
!
!  Discussion:
!
!    MANDELBROT computes an image of the Mandelbrot set.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 August 2009
!
!  Author:
!
!    John Burkardt
!
!  Local Parameters:
!
!    Local, integer ( kind = 4 ) COUNT_MAX, the maximum number of iterations
!    taken for a particular pixel.
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 501

  integer ( kind = 4 ) b(n,n)
  integer ( kind = 4 ) c
  integer ( kind = 4 ) c_max
  integer ( kind = 4 ) count(n,n)
  integer ( kind = 4 ), parameter :: count_max = 400
  character ( len = 255 ) :: filename = 'mandelbrot.ppm'
  integer ( kind = 4 ) g(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) r(n,n)
  real ( kind = 8 ) x
  real ( kind = 8 ) :: x_max =   1.25D+00
  real ( kind = 8 ) :: x_min = - 2.25D+00
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) y
  real ( kind = 8 ) :: y_max =   1.75D+00
  real ( kind = 8 ) :: y_min = - 1.75D+00
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2

  write ( *, '(a)' ) ' '
  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MANDELBROT'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Create an ASCII PPM image of the Mandelbrot set.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  For each point C = X + i*Y'
  write ( *, '(a,g14.6,a,g14.6,a)' ) '  with X range [', x_min, ',', x_max, ']'
  write ( *, '(a,g14.6,a,g14.6,a)' ) '  and  Y range [', y_min, ',', y_max, ']'
  write ( *, '(a,i8,a)' ) '  carry out ', count_max, ' iterations of the map'
  write ( *, '(a)' ) '  Z(n+1) = Z(n)^2 + C.'
  write ( *, '(a)' ) '  If the iterates stay bounded (norm less than 2)'
  write ( *, '(a)' ) '  then C is taken to be a member of the set.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  An ASCII PPM image of the set is created using'
  write ( *, '(a,i8,a)' ) '    N = ', n, ' pixels in the X direction and'
  write ( *, '(a,i8,a)' ) '    N = ', n, ' pixels in the Y direction.'
!
!  Carry out the iteration for each pixel, determining COUNT.
!
  do i = 1, n
    do j = 1, n

      x = ( real (     j - 1, kind = 8 ) * x_max   &
          + real ( n - j,     kind = 8 ) * x_min ) &
          / real ( n     - 1, kind = 8 )

      y = ( real (     i - 1, kind = 8 ) * y_max   &
          + real ( n - i,     kind = 8 ) * y_min ) &
          / real ( n     - 1, kind = 8 )

      count(i,j) = 0

      x1 = x
      y1 = y

      do k = 1, count_max

        x2 = x1 * x1 - y1 * y1 + x
        y2 = 2 * x1 * y1 + y

        if ( x2 < -2.0D+00 .or. &
              2.0D+00 < x2 .or. &
             y2 < -2.0D+00 .or. &
             2.0D+00 < y2 ) then
          count(i,j) = k
          exit
        end if

        x1 = x2
        y1 = y2

      end do

    end do
  end do
!
!  Determine the coloring of each pixel.
!
  c_max = maxval ( count(1:n,1:n) )

  do i = 1, n
    do j = 1, n
      if ( mod ( count(i,j), 2 ) == 1 ) then
        r(i,j) = 255
        g(i,j) = 255
        b(i,j) = 255
      else
        c = int ( 255.0D+00 * sqrt ( sqrt ( sqrt ( &
          ( real ( count(i,j), kind = 8 ) / real ( c_max, kind = 8 ) ) ) ) ) )
        r(i,j) = 3 * c / 5
        g(i,j) = 3 * c / 5
        b(i,j) = c
      end if

    end do
  end do
!
!  Write an image file.
!
  filename = 'mandelbrot.ppm'

  call ppma_write ( filename, n, n, r, g, b, ierror )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  ASCII PPM image data stored in "' // trim ( filename ) // '".'
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MANDELBROT'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
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
subroutine ppma_write ( file_out_name, row_num, col_num, r, g, b, ierror )

!*****************************************************************************80
!
!! PPMA_WRITE writes an ASCII portable pixel map file.
!
!  Discussion:
!
!    PPM files can be viewed by XV.
!
!    Programs to convert files to this format include:
!
!      GIFTOPPM  - GIF file
!      PGMTOPPM  - Portable Gray Map file
!      PICTTOPPM - Macintosh PICT file
!      XPMTOPPM  - X11 pixmap file
!
!    Various programs can convert other formats to PPM format, including:
!
!      BMPTOPPM - Microsoft Windows BMP file.
!
!    A PPM file can also be converted to other formats, by programs:
!
!      PPMTOACAD - AutoCAD file
!      PPMTOGIF  - GIF file
!      PPMTOPGM  - Portable Gray Map file
!      PPMTOPICT - Macintosh PICT file
!      PPMTOPUZZ - X11 puzzle file
!      PPMTORGB3 - 3 Portable Gray Map files
!      PPMTOXPM  - X11 pixmap file
!      PPMTOYUV  - Abekas YUV file
!    
!  Example:
!
!    P3
!    # feep.ppma created by PBMLIB(PPMA_WRITE).
!    4 4
!    15
!     0  0  0    0  0  0    0  0  0   15  0 15
!     0  0  0    0 15  7    0  0  0    0  0  0
!     0  0  0    0  0  0    0 15  7    0  0  0
!    15  0 15    0  0  0    0  0  0    0  0  0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_OUT_NAME, the name of the file 
!    to which the data should be written.
!
!    Input, integer ( kind = 4 ) ROW_NUM, COL_NUM, the number of rows 
!    and columns of data.
!
!    Input, integer ( kind = 4 ) R(ROW_NUM,COL_NUM), G(ROW_NUM,COL_NUM),
!    B(ROW_NUM,COL_NUM), the red, green and blue values of each pixel.  These 
!    should be positive.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!    0, no error.
!    1, the data was illegal.
!    2, the file could not be opened.
!
  implicit none

  integer ( kind = 4 ) col_num
  integer ( kind = 4 ) row_num

  integer ( kind = 4 ) b(row_num,col_num)
  logical, parameter :: debug = .false.
  character ( len = * ) file_out_name
  integer ( kind = 4 ) file_out_unit
  integer ( kind = 4 ) g(row_num,col_num)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) r(row_num,col_num)
  integer ( kind = 4 ) rgb_max

  ierror = 0
!
!  Compute the maximum color value.
!
  rgb_max = max ( &
    maxval ( r(1:row_num,1:col_num) ), &
    maxval ( g(1:row_num,1:col_num) ), &
    maxval ( b(1:row_num,1:col_num) ) )
!
!  Open the file.
!
  call get_unit ( file_out_unit )

  open ( unit = file_out_unit, file = file_out_name, status = 'replace', &
    form = 'formatted', access = 'sequential', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PPMA_WRITE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the file.'
    ierror = 2
    return
  end if
!
!  Write the header.
!
  call ppma_write_header ( file_out_name, file_out_unit, row_num, col_num, &
    rgb_max, ierror )
!
!  Write the data.
!
  call ppma_write_data ( file_out_unit, row_num, col_num, r, g, b, ierror )
!
!  Close the file.
!
  close ( unit = file_out_unit )
!
!  Report
!
  if ( debug ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PPMA_WRITE - Note:'
    write ( *, '(a)' ) '  The data was checked and written.'
    write ( *, '(a,i8)' ) '  Number of data rows ROW_NUM =    ', row_num
    write ( *, '(a,i8)' ) '  Number of data columns COL_NUM = ', col_num
    write ( *, '(a,i8)' ) '  Maximum RGB value RGB_MAX =      ', rgb_max
  end if

  return
end
subroutine ppma_write_data ( file_out_unit, row_num, col_num, r, g, b, ierror )

!*****************************************************************************80
!
!! PPMA_WRITE_DATA writes the data of a PPMA file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) FILE_OUT_UNIT, the output file unit number.
!
!    Input, integer ( kind = 4 ) ROW_NUM, COL_NUM, the number of rows 
!    and columns of data.
!
!    Input, integer ( kind = 4 ) R(ROW_NUM,COL_NUM), G(ROW_NUM,COL_NUM), 
!    B(ROW_NUM,COL_NUM), the red, green and blue values of each pixel.  These
!    should be positive.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!    0, no error.
!    1, the data was illegal.
!    2, the file could not be opened.
!
  implicit none

  integer ( kind = 4 ) col_num
  integer ( kind = 4 ) row_num

  integer ( kind = 4 ) b(row_num,col_num)
  integer ( kind = 4 ) file_out_unit
  integer ( kind = 4 ) g(row_num,col_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) r(row_num,col_num)

  ierror = 0
!
!  Write the header.
!
  do i = 1, row_num
    do jlo = 1, col_num, 4
      jhi = min ( jlo + 3, col_num )
      write ( file_out_unit, '(12i5)' ) ( r(i,j), g(i,j), b(i,j), j = jlo, jhi )
    end do
  end do

  return
end
subroutine ppma_write_header ( file_out_name, file_out_unit, row_num, col_num, &
  rgb_max, ierror )

!*****************************************************************************80
!
!! PPMA_WRITE_HEADER writes the header of a PPMA file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_OUT_NAME, the name of the output file.
!
!    Input, integer ( kind = 4 ) FILE_OUT_UNIT, the output file unit number.
!
!    Input, integer ( kind = 4 ) ROW_NUM, COL_NUM, the number of rows and 
!    columns of data.
!
!    Input, integer ( kind = 4 ) RGB_MAX, the maximum value of any data component.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!    0, no error.
!    1, the data was illegal.
!    2, the file could not be opened.
!
  implicit none

  character ( len = * )  file_out_name
  integer ( kind = 4 ) file_out_unit
  integer ( kind = 4 ) ierror
  character ( len = 2 ) :: magic = 'P3'
  integer ( kind = 4 ) col_num
  integer ( kind = 4 ) row_num
  integer ( kind = 4 ) rgb_max

  ierror = 0
!
!  Write the header.
!
  write ( file_out_unit, '(a2)' ) magic
  write ( file_out_unit, '(a)' ) '# ' // trim ( file_out_name ) &
    // ' created by PPMA_WRITE.F90.'
  write ( file_out_unit, '(i5,2x,i5)' ) col_num, row_num
  write ( file_out_unit, '(i5)' ) rgb_max

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
