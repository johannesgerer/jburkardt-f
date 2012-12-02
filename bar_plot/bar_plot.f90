subroutine angle_to_rgb ( angle, r, g, b )

!*****************************************************************************80
!
!! ANGLE_TO_RGB returns a color on the perimeter of the color hexagon.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 December 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 )ANGLE, the angle in the color hexagon.  
!    The sextants are defined by the following points:
!        0 degrees, 1, 0, 0, red;
!       60 degrees, 1, 1, 0, yellow;
!      120 degrees, 0, 1, 0, green;
!      180 degrees, 0, 1, 1, cyan;
!      240 degrees, 0, 0, 1, blue;
!      300 degrees, 1, 0, 1, magenta.
!
!    Output, real ( kind = 8 )R, G, B, RGB specifications for the color 
!    that lies at the given angle, on the perimeter of the color hexagon.  One
!    value will be 1, and one value will be 0.
!
  implicit none

  real ( kind = 8 ) angle
  real ( kind = 8 ) angle2
  real ( kind = 8 ) b
  real ( kind = 8 ) g
  real ( kind = 8 ), parameter :: degrees_to_radians = &
    3.14159265D+00 / 180.0D+00
  real ( kind = 8 ) r

  angle = mod ( angle, 360.0D+00 )

  if ( angle < 0.0D+00 ) then
    angle = angle + 360.0D+00
  end if

  if ( angle <= 60.0D+00 ) then

    angle2 = degrees_to_radians * 3.0D+00 * angle / 4.0D+00
    r = 1.0D+00
    g = tan ( angle2 )
    b = 0.0D+00

  else if ( angle <= 120.0D+00 ) then

    angle2 = degrees_to_radians * 3.0D+00 * angle / 4.0D+00
    r = cos ( angle2 ) / sin ( angle2 )
    g = 1.0D+00
    b = 0.0D+00

  else if ( angle <= 180.0D+00 ) then

    angle2 = degrees_to_radians * 3.0D+00 * ( angle - 120.0D+00 ) / 4.0D+00
    r = 0.0D+00
    g = 1.0D+00
    b = tan ( angle2 )

  else if ( angle <= 240.0D+00 ) then

    angle2 = degrees_to_radians * 3.0D+00 * ( angle - 120.0D+00 ) / 4.0D+00
    r = 0.0D+00
    g = cos ( angle2 ) / sin ( angle2 )
    b = 1.0D+00

  else if ( angle <= 300.0D+00 ) then

    angle2 = degrees_to_radians * 3.0D+00 * ( angle - 240.0D+00 ) / 4.0D+00
    r = tan ( angle2 )
    g = 0.0D+00
    b = 1.0D+00

  else if ( angle <= 360.0D+00 ) then

    angle2 = degrees_to_radians * 3.0D+00 * ( angle - 240.0D+00 ) / 4.0D+00
    r = 1.0D+00
    g = 0.0D+00
    b = cos ( angle2 ) / sin ( angle2 )

  end if

  return
end
subroutine bar_data_example ( nx, ny, ytab )

!*****************************************************************************80
!
!! BAR_DATA_EXAMPLE returns some sample bar data.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 June 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters
!
!    Input, integer NX, the number of X sample points.
!
!    Input, integer NY, the number of Y divisions.
!
!    Output, real ( kind = 8 ) YTAB(NX,NY), a set of values with the property
!    that, for each I (row index or "X"), the values YTAB(I,1:NY) are
!    all nonnegative, and sum to no more than 1.
!
  implicit none

  integer nx
  integer ny

  integer i
  integer j
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ytab(nx,ny)

  ymax = 0.0D+00

  do i = 1, nx

    x = 2.0D+00 * pi * real ( i - 1, kind = 8 ) &
      / real ( nx - 1, kind = 8 )

    do j = 1, ny
      ytab(i,j) = 1.0D+00 + sin ( j * x )
    end do

    ymax = max ( ymax, sum ( ytab(i,1:ny) ) )

  end do
  
  ytab(1:nx,1:ny) = ytab(1:nx,1:ny) / ymax

  return
end
subroutine bar_data_to_rgb ( nx, ny, ytab, nrow, ncol, r, g, b )

!*****************************************************************************80
!
!! BAR_DATA_TO_RGB makes RGB arrays from bar graph data.
!
!  Discussion:
!
!    This routine makes some fairly strict assumptions.  
!
!    First, no X data is supplied, because it is assumed that the
!    data is sampled at equal spacings, so particular X values are
!    irrelevant.
!
!    Secondly, it is assumed that the Y data is all nonnegative,
!    and scaled in such a way that the maximum sum for any particular
!    I index is 1.  Under this assumption, the program needs to do
!    no scaling of the Y data, and can determine the color of each
!    pixel simply by adding up just enough Y's to exceed the
!    imputed Y value of the pixel.
!
!    Once the RGB arrays have been set up, a plot can be made by
!    calling, for example, the routine PPMA_WRITE from PPMA_IO.
!
!  Example:
!
!    Try really hard, and you can see a bar graph in the "plot" below.
!    We have to imagine that the three symbols $, * and - represent
!    colors.  
!
!    *-----------
!    *-------*---
!    **-----***--
!    $***-***$***
!    $$*****$$$*$
!    $$$**$$$$$$$
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 June 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters
!
!    Input, integer NX, the number of X sample points.
!
!    Input, integer NY, the number of Y divisions.
!
!    Input, real ( kind = 8 )YTAB(NX,NY), a set of values with the property
!    that, for each I (row index or "X"), the values YTAB(I,1:NY) are
!    all nonnegative, and sum to no more than 1.
!
!    Input, integer NROW, NCOL, the number of pixels per row and
!    column of the image.
!
!    Output, integer R(NROW,NCOL), G(NROW,NCOL), B(NROW,NCOL), the
!    RGB arrays for the plot.
!
  implicit none

  integer ncol
  integer nrow
  integer nx
  integer ny

  real ( kind = 8 ) angle
  integer b(nrow,ncol)
  integer g(nrow,ncol)
  integer ip
  integer it
  integer jp
  integer jt
  integer r(nrow,ncol)
  real ( kind = 8 ) rb
  real ( kind = 8 ) rg
  real ( kind = 8 ) rr
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) ysum
  real ( kind = 8 ) ytab(nx,ny)
!
!  Proceed through the columns of the graph.  Each column
!  corresponds to a particular value of X, and a set of Y values.
!
  do jp = 1, ncol

    x = real ( 2 * jp - 1, kind = 8 ) &
      / real ( 2 * ( ncol - 1 ), kind = 8 )

    it = nint ( real ( nx * ( 2 * jp - 1 ) + ncol, kind = 8 ) &
      / real ( 2 * ncol, kind = 8 ) )

    write ( *, * ) jp, it

    jt = 1
    ysum = ytab(it,jt)

    do ip = 1, nrow

      y = real ( 2 * ip - 1, kind = 8 ) / real ( 2 * ( nrow - 1 ), kind = 8 )

      do while ( ysum < y )

        jt = jt + 1

        if ( ny < jt ) then
          exit 
        end if

        ysum = ysum + ytab(it,jt)

      end do

      if ( ny < jt ) then
        rr = 0.80D+00
        rg = 0.80D+00
        rb = 0.85D+00
      else
        call i4_to_angle ( jt-1, angle )
        call angle_to_rgb ( angle, rr, rg, rb )
      end if

      r(nrow+1-ip,jp) = int ( 255 * rr )
      g(nrow+1-ip,jp) = int ( 255 * rg )
      b(nrow+1-ip,jp) = int ( 255 * rb )

    end do

  end do

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
!    18 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer IUNIT, the free unit number.
!
  implicit none

  integer i
  integer ios
  integer iunit
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
function i4_log_2 ( i )

!*****************************************************************************80
!
!! I4_LOG_2 returns the integer part of the logarithm base 2 of |I|.
!
!  Example:
!
!     I  I4_LOG_2
!
!     0  -1
!     1,  0
!     2,  1
!     3,  1
!     4,  2
!     5,  2
!     6,  2
!     7,  2
!     8,  3
!     9,  3
!    10,  3
!   127,  6
!   128,  7
!   129,  7
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 January 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer I, the number whose logarithm base 2 is desired.
!
!    Output, integer I4_LOG_2, the integer part of the logarithm base 2 of
!    the absolute value of I.
!    For positive I4_LOG_2(I), it should be true that
!      2**I4_LOG_2(X) <= |I| < 2**(I4_LOG_2(I)+1).
!    The special case of I4_LOG_2(0) returns -HUGE().
!
  implicit none

  integer i4_log_2
  integer i
  integer i_abs

  if ( i == 0 ) then

    i4_log_2 = - huge ( i4_log_2 )

  else

    i4_log_2 = 0

    i_abs = abs ( i )

    do while ( 2 <= i_abs )
      i_abs = i_abs / 2
      i4_log_2 = i4_log_2 + 1
    end do

  end if

  return
end
subroutine i4_to_angle ( i, angle )

!*****************************************************************************80
!
!! I4_TO_ANGLE maps integers to points on a circle.
!
!  Discussion:
!
!    The angles are intended to be used to select colors on a color
!    hexagon whose 6 vertices are red, yellow, green, cyan, blue,
!    magenta.
!
!  Example:
!
!     I   X      ANGLE
!
!     0   0/3      0
!     1   1/3    120
!     2   2/3    240
!
!     3   1/6     60
!     4   3/6    180
!     5   5/6    300
!
!     6   1/12    30
!     7   3/12    90
!     8   5/12   150
!     9   7/12   210
!    10   9/12   270
!    11  11/12   330
!
!    12   1/24    15
!    13   3/24    45
!    14   5/24    75
!    etc
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 January 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer I, the index of the desired color.
!
!    Output, real ( kind = 8 )ANGLE, an angle, measured in degrees,
!    between 0 and 360.
!
  implicit none

  real ( kind = 8 ) angle
  integer i
  integer i4_log_2
  integer i1
  integer i2
  integer i3
  integer i4

  if ( 0 <= abs ( i ) .and. abs ( i ) <= 2 ) then

    angle = 120.0D+00 * real ( abs ( i ), kind = 8 )

  else

    i1 = i4_log_2 ( abs ( i ) / 3 )
    i2 = abs ( i ) + 1 - 3 * 2**i1
    i3 = 2 * ( i2 - 1 ) + 1
    i4 = 3 * 2**( i1 + 1 )

    angle = 360.0D+00 * real ( i3, kind = 8 ) / real ( i4, kind = 8 )

  end if

  return
end
subroutine ppm_check_data ( r, g, b, ierror, maxcol, ncol, nrow )

!*****************************************************************************80
!
!! PPM_CHECK_DATA checks pixel data.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer R(NROW,NCOL), G(NROW,NCOL), B(NROW,NCOL), contains the
!    RGB pixel data.
!
!    Output, integer IERROR, error flag.
!    0, no error detected.
!    1, the data is illegal.
!
!    Input, integer MAXCOL, the maximum value.
!
!    Input, integer NCOL, NROW, the number of rows and columns of data.
!
  implicit none

  integer ncol
  integer nrow

  integer b(nrow,ncol)
  integer g(nrow,ncol)
  integer i
  integer ierror
  integer j
  integer maxcol
  integer r(nrow,ncol)

  ierror = 0
!
!  Make sure no color is negative.
!
  if ( minval ( r(1:nrow,1:ncol) ) < 0 .or. &
       minval ( g(1:nrow,1:ncol) ) < 0 .or. &
       minval ( b(1:nrow,1:ncol) ) < 0 ) then
    ierror = 1
    return
  end if
!
!  Make sure no color is greater than MAXCOL.
!
  if ( maxcol < maxval ( r(1:nrow,1:ncol) ) .or. &
       maxcol < maxval ( g(1:nrow,1:ncol) ) .or. &
       maxcol < maxval ( b(1:nrow,1:ncol) ) ) then
    ierror = 1
    return
  end if

  return
end
subroutine ppma_write ( file_name, ierror, nrow, ncol, r, g, b )

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
!    07 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the file to which
!    the data should be written.
!
!    Output, integer IERROR, an error flag.
!    0, no error.
!    1, the data was illegal.
!    2, the file could not be opened.
!
!    Input, integer NROW, NCOL, the number of rows and columns of data.
!
!    Input, integer R(NROW,NCOL), G(NROW,NCOL), B(NROW,NCOL), contain
!    the red, green and blue values of each pixel.  These should
!    be positive.
!
  implicit none

  integer ncol
  integer nrow

  integer b(nrow,ncol)
  logical, parameter :: debug = .false.
  character ( len = * ) file_name
  integer g(nrow,ncol)
  integer i
  integer ierror
  integer ios
  integer j
  integer jhi
  integer jlo
  character ( len = 2 ) magic
  integer maxcol
  integer output_unit
  integer r(nrow,ncol)

  ierror = 0
!
!  Compute the maximum color value.
!
  maxcol = max ( &
    maxval ( r(1:nrow,1:ncol) ), &
    maxval ( g(1:nrow,1:ncol) ), &
    maxval ( b(1:nrow,1:ncol) ) )
!
!  Check the data.
!
  call ppm_check_data ( r, g, b, ierror, maxcol, ncol, nrow )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PPMA_WRITE - Fatal error!'
    write ( *, '(a)' ) '  Bad data detected by PPM_CHECK_DATA!'
    ierror = 1
    return
  end if
!
!  Open the file.
!
  call get_unit ( output_unit )

  open ( unit = output_unit, file = file_name, status = 'replace', &
    form = 'formatted', access = 'sequential', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PPMA_WRITE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the file.'
    ierror = 2
    return
  end if
!
!  Write the data.
!
  magic = 'P3'

  write ( output_unit, '(a2)' ) magic
  write ( output_unit, '(a)' ) '# ' // trim ( file_name ) &
    // ' created by PPMLIB(PPMA_WRITE).'
  write ( output_unit, '(i5,2x,i5)' ) ncol, nrow
  write ( output_unit, '(i5)' ) maxcol

  do i = 1, nrow
    do jlo = 1, ncol, 4
      jhi = min ( jlo + 3, ncol )
      write ( output_unit, '(12i5)' ) ( r(i,j), g(i,j), b(i,j), j = jlo, jhi )
    end do
  end do
!
!  Close the file.
!
  close ( unit = output_unit )
!
!  Report
!
  if ( debug ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PPMA_WRITE - Note:'
    write ( *, '(a)' ) '  The data was checked and written.'
    write ( *, '(a,i8)' ) '  Number of data rows NROW =    ', nrow
    write ( *, '(a,i8)' ) '  Number of data columns NCOL = ', ncol
    write ( *, '(a,i8)' ) '  Maximum color value MAXCOL =  ', maxcol
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

  character ( len = 8 ) ampm
  integer d
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer n
  integer s
  integer values(8)
  integer y

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
