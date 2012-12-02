subroutine ch_cap ( c )

!*****************************************************************************80
!
!! CH_CAP capitalizes a single character.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character C, the character to capitalize.
!
  implicit none

  character c
  integer ( kind = 4 ) itemp

  itemp = ichar ( c )

  if ( 97 <= itemp .and. itemp <= 122 ) then
    c = char ( itemp - 32 )
  end if

  return
end
function ch_eqi ( c1, c2 )

!*****************************************************************************80
!
!! CH_EQI is a case insensitive comparison of two characters for equality.
!
!  Example:
!
!    CH_EQI ( 'A', 'a' ) is .TRUE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C1, C2, the characters to compare.
!
!    Output, logical CH_EQI, the result of the comparison.
!
  implicit none

  character c1
  character c1_cap
  character c2
  character c2_cap
  logical ch_eqi

  c1_cap = c1
  c2_cap = c2

  call ch_cap ( c1_cap )
  call ch_cap ( c2_cap )

  if ( c1_cap == c2_cap ) then
    ch_eqi = .true.
  else
    ch_eqi = .false.
  end if

  return
end
subroutine ch_to_digit ( c, digit )

!*****************************************************************************80
!
!! CH_TO_DIGIT returns the value of a base 10 digit.
!
!  Example:
!
!     C   DIGIT
!    ---  -----
!    '0'    0
!    '1'    1
!    ...  ...
!    '9'    9
!    ' '    0
!    'X'   -1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C, the decimal digit, '0' through '9' or blank
!    are legal.
!
!    Output, integer ( kind = 4 ) DIGIT, the corresponding value.  If C was
!    'illegal', then DIGIT is -1.
!
  implicit none

  character c
  integer ( kind = 4 ) digit

  if ( lge ( c, '0' ) .and. lle ( c, '9' ) ) then

    digit = ichar ( c ) - 48

  else if ( c == ' ' ) then

    digit = 0

  else

    digit = -1

  end if

  return
end
function degrees_to_radians ( angle_deg )

!*****************************************************************************80
!
!! DEGREES_TO_RADIANS converts an angle from degrees to radians.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ANGLE_DEG, an angle in degrees.
!
!    Output, real ( kind = 8 ) DEGREES_TO_RADIANS, the equivalent angle
!    in radians.
!
  implicit none

  real ( kind = 8 ) angle_deg
  real ( kind = 8 ) degrees_to_radians
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00

  degrees_to_radians = ( angle_deg / 180.0D+00 ) * pi

  return
end
subroutine dist_table_check ( n, dist_table, check )

!*****************************************************************************80
!
!! DIST_TABLE_CHECK checks a distance table.
!
!  Discussion:
!
!    1) All entries must be nonnegative.
!    2) Diagonal entries must be zero.
!    3) Off-diagonal entries must be symmetric.
!    4) The triangle inequality must be observed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of cities.
!
!    Input, real ( kind = 8 ) DIST_TABLE(N,N), the distance table.
!
!    Output, integer ( kind = 4 ) CHECK, the result of the check.
!    0, the matrix passed the checks.
!    1, Not all entries are nonnegative.
!    2, Not all diagonal entries are zero.
!    3, Not all off-diagonal entries are symmetric.
!    4, Not all entries satisfy the triangle inequality.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) check
  real ( kind = 8 ) dist_table(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k

  if ( any ( dist_table(1:n,1:n) < 0.0D+00 ) ) then
    check = 1
    return
  end if

  do i = 1, n
    if ( dist_table(i,i) /= 0.0D+00 ) then
      check = 2
      return
    end if
  end do

  do i = 1, n
    do j = 1, i - 1
      if ( dist_table(i,j) /= dist_table(j,i) ) then
        check = 3
        return
      end if
    end do
  end do

  do i = 1, n
    do j = 1, n
      do k = 1, n
        if ( dist_table(i,j) + dist_table(j,k) < dist_table(i,k) ) then
          check = 4
          return
        end if
      end do
    end do
  end do

  check = 0

  return
end
subroutine dms_print ( n, lat_dms, long_dms, title )

!*****************************************************************************80
!
!! DMS_PRINT prints the latitude and longitude in degrees/minutes/seconds.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 November 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data items.
!
!    Input, integer ( kind = 4 ) LAT_DMS(4,N), LONG_DMS(4,N), the latitudes
!    and longitudes, in degrees, minutes and seconds.  The fourth
!    argument is +1/-1 for North/South latitude or East/West longitude.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  character lat_char
  integer ( kind = 4 ) lat_dms(4,n)
  character long_char
  integer ( kind = 4 ) long_dms(4,n)
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  #    Latitude            Longitude'
  write ( *, '(a)' ) '       (Deg/Min/Sec)       (Deg/Min/Sec)'
  write ( *, '(a)' ) '---  ---------------      ---------------'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(i3,2x,3i4,2x,a1,5x,3i4,2x,a1)' ) i, &
      lat_dms(1:3,i), lat_char ( lat_dms(4,i) ), &
      long_dms(1:3,i), long_char ( long_dms(4,i) )
  end do

  return
end
subroutine dms_read ( file_name, n, lat_dms, long_dms )

!*****************************************************************************80
!
!! DMS_READ reads DMS data from a file.
!
!  Discussion:
!
!    Blank lines and comment lines (beginning with '#') are ignored.
!
!    Individual data values should be separated by spaces or commas.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 November 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the file to read.
!
!    Input, integer ( kind = 4 ) N, the number of data items.
!
!    Output, integer ( kind = 4 ) LAT_DMS(4,N), LONG_DMS(4,N), the latitude and
!    longitudes, in degrees, minutes and seconds.  The fourth
!    argument is +1/-1 for North/South latitude or East/West longitude. 
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: ncol = 8

  character ( len = * ) file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) input
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) j
  integer ( kind = 4 ) last
  integer ( kind = 4 ) lat_dms(4,n)
  integer ( kind = 4 ) length
  character ( len = 80 ) line
  integer ( kind = 4 ) line_num
  integer ( kind = 4 ) long_dms(4,n)
  logical s_eqi
  integer ( kind = 4 ) value
  integer ( kind = 4 ) vector(ncol)
  character ( len = 80 ) word

  call get_unit ( input )

  open ( unit = input, file = file_name, status = 'old' )

  lat_dms(1:4,1:n) = huge ( lat_dms(1,1) )
  long_dms(1:4,1:n) = huge ( long_dms(1,1) )

  i = 1
  j = 0
  line_num = 0

  do
!
!  Have we read enough data?
!
    if ( i == n .and. j == ncol ) then
      exit
    end if
!
!  Have we read too much data?
!
    if ( n < i .or. ncol < j ) then
      exit
    end if
!
!  Read the next line from the file.
!
    read ( input, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      exit
    end if

    line_num = line_num + 1
!
!  Skip blank lines and comment lines.
!
    if ( len_trim ( line ) == 0 ) then

    else if ( line(1:1) == '#' ) then

    else
!
!  LAST points to the last character associated with the previous 
!  data value read from the line.
!
      last = 0

      do
!
!  Try to read another value from the line.
!  Note that, confusingly, J right now indicates the column
!  of the PREVIOUS data item that was read.
!
        if ( j == 3 ) then

          call s_to_w ( line(last+1:), word, ierror, length )

          if ( ierror /= 0 ) then
            exit
          end if

          if ( s_eqi ( word(1:1), 'N' ) ) then
            value = +1
          else if ( s_eqi ( word(1:1), 'S' ) ) then
            value = -1
          else
            value = 0
          end if

        else if ( j == 7 ) then

          call s_to_w ( line(last+1:), word, ierror, length )

          if ( ierror /= 0 ) then
            exit
          end if

          if ( s_eqi ( word(1:1), 'E' ) ) then
            value = +1
          else if ( s_eqi ( word(1:1), 'W' ) ) then
            value = -1
          else
            value = 0
          end if

        else

          call s_to_i4 ( line(last+1:), value, ierror, length )

          if ( ierror /= 0 ) then
            exit
          end if

        end if
!
!  Update the pointer.
!
        last = last + length
!
!  If we read a new value, where do we put it?
!
        j = j + 1

        if ( ncol < j ) then
          j = 1
          i = i + 1
          if ( n < i ) then
            exit
          end if
        end if

        vector(j) = value
!
!  If you reached the end of the row, it's time to read a new line.
!
        if ( j == ncol ) then
          lat_dms(1:4,i) = vector(1:4)
          long_dms(1:4,i) = vector(5:8)
          exit
        end if

      end do

    end if

  end do

  close ( unit = input )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DMS_READ:'
  write ( *, '(a,i6)' ) '  Number of lines read was ', line_num

  return
end
subroutine dms_to_dist ( n, lat_dms, long_dms, dist_table )

!*****************************************************************************80
!
!! DMS_TO_DIST creates a distance table from latitudes and longitudes.
!
!  Discussion:
!
!    A distance function is used which is appropriate for the earth.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 November 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data items.
!
!    Input, integer ( kind = 4 ) LAT_DMS(4,N), LONG_DMS(4,N), the latitude 
!    and longitude, in degrees, minutes, and seconds, for each point.
!    The fourth argument is +1/-1 for North/South latitude or 
!    East/West longitude.
!
!    Output, real ( kind = 8 ) DIST_TABLE(N,N), the distance matrix.  Distances 
!    are measured in miles.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) dist_table(n,n)
  integer ( kind = 4 ) lat_dms(4,n)
  integer ( kind = 4 ) long_dms(4,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  do i = 1, n
    dist_table(i,i) = 0.0D+0
    do j = i+1, n
      call dms_to_distance_earth ( lat_dms(1,i), long_dms(1,i), &
        lat_dms(1,j), long_dms(1,j), dist_table(i,j) )
      dist_table(j,i) = dist_table(i,j)
    end do
  end do

  return
end
subroutine dms_to_distance_earth ( lat1_dms, long1_dms, lat2_dms, &
  long2_dms, dist )

!*****************************************************************************80
!
!! DMS_TO_DISTANCE_EARTH finds the distance between two points on the earth.
!
!  Discussion:
!
!    The positions of the the points are given as longitude and
!    latitude, measured in degrees, minutes, and seconds.
!
!    The distance is measured on the surface of the earth, which
!    is approximated by a sphere.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 November 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LAT1_DMS(4), LONG1_DMS(4), the latitude and 
!    longitude of the first point.  The fourth
!    argument is +1/-1 for North/South latitude or East/West longitude.
!
!    Input, integer ( kind = 4 ) LAT2_DMS(4), LONG2_DMS(4), the latitude and 
!    longitude of the second point.  The fourth
!    argument is +1/-1 for North/South latitude or East/West longitude.
!
!    Output, real ( kind = 8 ) DIST, the distance between the points, in miles.
!
  implicit none

  real ( kind = 8 ) dist
  real ( kind = 8 ) dms_to_radians
  integer ( kind = 4 ) lat1_dms(4)
  real ( kind = 8 ) lat1_rad
  integer ( kind = 4 ) lat2_dms(4)
  real ( kind = 8 ) lat2_rad
  integer ( kind = 4 ) long1_dms(4)
  real ( kind = 8 ) long1_rad
  integer ( kind = 4 ) long2_dms(4)
  real ( kind = 8 ) long2_rad
  real ( kind = 8 ), parameter :: radius = 3958.89D+00
  real ( kind = 8 ) theta

  lat1_rad = dms_to_radians ( lat1_dms )
  long1_rad = dms_to_radians ( long1_dms )

  lat2_rad = dms_to_radians ( lat2_dms )
  long2_rad = dms_to_radians ( long2_dms )

  theta = acos ( sin ( lat1_rad ) * sin ( lat2_rad ) &
               + cos ( lat1_rad ) * cos ( lat2_rad ) &
               * cos ( long1_rad - long2_rad ) )

  dist = radius * theta

  return
end
function dms_to_radians ( dms )

!*****************************************************************************80
!
!! DMS_TO_RADIANS converts degrees, minutes, seconds to radians.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 November 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DMS(4), the measurement of an angle in
!    degrees, minutes and seconds.  The fourth
!    argument is +1/-1 for North/South latitude or East/West longitude.
!
!    Output, real ( kind = 8 ) DMS_TO_RADIANS, the measurement of the same
!    angle in radians.
!
  implicit none

  real ( kind = 8 ) d_real
  integer ( kind = 4 ) dms(4)
  real ( kind = 8 ) dms_to_radians
  integer ( kind = 4 ) i4_sign
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00

  d_real = real ( i4_sign ( dms(4) ), kind = 8 ) * &
    real ( dms(3) + 60 * dms(2) + 3600 * dms(1), kind = 8 ) / 3600.0D+00

  dms_to_radians = pi * d_real / 180.0D+00

  return
end
subroutine dms_to_xy ( n, lat_dms, long_dms, point_xy )

!*****************************************************************************80
!
!! DMS_TO_XY: Latitude/Longitude in DMS to XY coordinates.
!
!  Discussion:
!
!    Essentially, the latitude and longitude information is treated
!    as though the Earth were a cylinder.  As long as the the
!    data is relatively close on the sphere (and far from either
!    pole!) the distortion will not be too severe.  If the data
!    is closely clustered, and also near the equator, the
!    positions will be relatively accurate.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data items.
!
!    Input, integer ( kind = 4 ) LAT_DMS(4,N), LONG_DMS(4,N), the latitude and
!    longitude, in degrees, minutes, and seconds, for each point.
!    The fourth argument is +1/-1 for North/South latitude or 
!    East/West longitude.
!
!    Output, real ( kind = 8 ) POINT_XY(2,N), the point coordinates, in miles.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) dms_to_radians
  integer ( kind = 4 ) lat_dms(4,n)
  integer ( kind = 4 ) long_dms(4,n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) phi
  real ( kind = 8 ) point_xy(2,n)
  real ( kind = 8 ), parameter :: radius = 3958.89D+00
  real ( kind = 8 ) theta

  do i = 1, n
    theta = dms_to_radians ( long_dms(1,i) )
    phi = dms_to_radians ( lat_dms(1,i) )
    point_xy(1,i) = radius * theta
    point_xy(2,i) = radius * phi
  end do

  return
end
subroutine dms_write ( file_name, n, lat_dms, long_dms )

!*****************************************************************************80
!
!! DMS_WRITE writes a DMS latitude, longitude file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 November 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the file to write.
!
!    Input, integer ( kind = 4 ) N, the number of data items.
!
!    Input, integer ( kind = 4 ) LAT_DMS(4,N), LONG_DMS(4,N), the data values.
!
  implicit none

  integer ( kind = 4 ) n

  character ( len = * )  file_name
  integer ( kind = 4 ) i
  character              lat_char
  integer ( kind = 4 ) lat_dms(4,n)
  character              long_char
  integer ( kind = 4 ) long_dms(4,n)
  integer ( kind = 4 ) output

  call get_unit ( output )

  open ( unit = output, file = file_name, status = 'replace' )

  write ( output, '(a)' ) '# ' // trim ( file_name )
  write ( output, '(a)' ) '#'
  write ( output, '(a)' ) '#  Created by DMS_WRITE.'
  write ( output, '(a)' ) '#'
  write ( output, '(a)' ) '#  Latitude, Longitude in degrees, minutes, seconds'
  write ( output, '(a,i6)' ) '#  Number of points N is ', n
  write ( output, '(a)' ) '#'

  do i = 1,  n

    write ( output, '(3i5,2x,3i5)' ) &
      lat_dms(1:3,i), lat_char ( lat_dms(4,i) ), &
      long_dms(1:3,i), long_char ( long_dms(4,i) )

  end do


  close ( unit = output )

  return
end
subroutine file_column_count ( input_filename, column_num )

!*****************************************************************************80
!
!! FILE_COLUMN_COUNT counts the number of columns in the first line of a file.
!
!  Discussion:
!
!    The file is assumed to be a simple text file.
!
!    Most lines of the file is presumed to consist of COLUMN_NUM words,
!    separated by spaces.  There may also be some blank lines, and some
!    comment lines,
!    which have a "#" in column 1.
!
!    The routine tries to find the first non-comment non-blank line and
!    counts the number of words in that line.
!
!    If all lines are blanks or comments, it goes back and tries to analyze
!    a comment line.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 June 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILENAME, the name of the file.
!
!    Output, integer ( kind = 4 ) COLUMN_NUM, the number of columns in the file.
!
  implicit none

  integer ( kind = 4 ) column_num
  logical got_one
  character ( len = * ) input_filename
  integer ( kind = 4 ) input_status
  integer ( kind = 4 ) input_unit
  character ( len = 255 ) line
!
!  Open the file.
!
  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_filename, status = 'old', &
    form = 'formatted', access = 'sequential', iostat = input_status )

  if ( input_status /= 0 ) then
    column_num = -1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_COLUMN_COUNT - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the input file "' &
      // trim ( input_filename ) // '" on unit ', input_unit
    return
  end if
!
!  Read one line, but skip blank lines and comment lines.
!
  got_one = .false.

  do

    read ( input_unit, '(a)', iostat = input_status ) line

    if ( input_status /= 0 ) then
      exit
    end if

    if ( len_trim ( line ) == 0 ) then
      cycle
    end if

    if ( line(1:1) == '#' ) then
      cycle
    end if

    got_one = .true.
    exit

  end do

  if ( .not. got_one ) then

    rewind ( input_unit )

    do

      read ( input_unit, '(a)', iostat = input_status ) line

      if ( input_status /= 0 ) then
        exit
      end if

      if ( len_trim ( line ) == 0 ) then
        cycle
      end if

      got_one = .true.
      exit

    end do

  end if

  close ( unit = input_unit )

  if ( .not. got_one ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_COLUMN_COUNT - Warning!'
    write ( *, '(a)' ) '  The file does not seem to contain any data.'
    column_num = -1
    return
  end if

  call s_word_count ( line, column_num )

  return
end
function file_exist ( file_name )

!*****************************************************************************80
!
!! FILE_EXIST reports whether a file exists.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the file.
!
!    Output, logical FILE_EXIST, is TRUE if the file exists.
!
  implicit none

  logical file_exist
  character ( len = * ) file_name

  inquire ( file = file_name, exist = file_exist )

  return
end
subroutine file_row_count ( input_filename, row_num )

!*****************************************************************************80
!
!! FILE_ROW_COUNT counts the number of row records in a file.
!
!  Discussion:
!
!    It does not count lines that are blank, or that begin with a
!    comment symbol '#'.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILENAME, the name of the input file.
!
!    Output, integer ( kind = 4 ) ROW_NUM, the number of rows found.
!
  implicit none

  integer ( kind = 4 ) bad_num
  integer ( kind = 4 ) comment_num
  integer ( kind = 4 ) ierror
  character ( len = * ) input_filename
  integer ( kind = 4 ) input_status
  integer ( kind = 4 ) input_unit
  character ( len = 255 ) line
  integer ( kind = 4 ) record_num
  integer ( kind = 4 ) row_num

  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_filename, status = 'old', &
    iostat = input_status )

  if ( input_status /= 0 ) then
    row_num = -1;
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_ROW_COUNT - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the input file "' // &
      trim ( input_filename ) // '" on unit ', input_unit
    stop
  end if

  comment_num = 0
  row_num = 0
  record_num = 0
  bad_num = 0

  do

    read ( input_unit, '(a)', iostat = input_status ) line

    if ( input_status /= 0 ) then
      ierror = record_num
      exit
    end if

    record_num = record_num + 1

    if ( line(1:1) == '#' ) then
      comment_num = comment_num + 1
      cycle
    end if

    if ( len_trim ( line ) == 0 ) then
      comment_num = comment_num + 1
      cycle
    end if

    row_num = row_num + 1

  end do

  close ( unit = input_unit )

  return
end
subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is a valule between 1 and 99 which
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
function i4_sign ( x )

!*****************************************************************************80
!
!! I4_SIGN evaluates the sign of an I4.
!
!  Discussion:
!
!    This function differs from the intrinsic SIGN function, because
!    it returns a value of 0 if the input argument is 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) X, the number whose sign is desired.
!
!    Output, integer ( kind = 4 ) I4_SIGN, a result based on the sign of X:
!
!    -1, if X < 0.
!     0, if X = 0.
!    +1, if X > 0.
!
  implicit none

  integer ( kind = 4 ) i4_sign
  integer ( kind = 4 ) x

  if ( x < 0 ) then
    i4_sign = -1
  else if ( x == 0 ) then
    i4_sign = 0
  else if ( 0 < x ) then
    i4_sign = +1
  end if

  return
end
function i4_to_a ( i )

!*****************************************************************************80
!
!! I4_TO_A returns the I-th alphabetic character.
!
!  Example:
!
!    I  I4_TO_A
!
!   -8  ' '
!    0  ' '
!    1  'A'
!    2  'B'
!   ..
!   26  'Z'
!   27  'a'
!   52  'z'
!   53  ' '
!   99  ' '
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 February 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the index of the letter to be returned.
!    0 is a space;
!    1 through 26 requests 'A' through 'Z', (ASCII 65:90);
!    27 through 52 requests 'a' through 'z', (ASCII 97:122);
!
!    Output, character I4_TO_A, the requested alphabetic letter.
!
  implicit none

  integer ( kind = 4 ), parameter :: cap_shift = 64
  integer ( kind = 4 ) i
  character i4_to_a
  integer ( kind = 4 ), parameter :: low_shift = 96

  if ( i <= 0 ) then
    i4_to_a = ' '
  else if ( 1 <= i .and. i <= 26 ) then
    i4_to_a = char ( cap_shift + i )
  else if ( 27 <= i .and. i <= 52 ) then
    i4_to_a = char ( low_shift + i - 26 )
  else if ( 53 <= i ) then
    i4_to_a = ' '
  end if

  return
end
function lat_char ( i )

!*****************************************************************************80
!
!! LAT_CHAR returns a character for negative or positive latitude.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 November 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, is negative for negative latitude, 
!    and positive for positive latitude.
!
!    Output, character LAT_CHAR, is 'S' for negative latitude, and
!    'N' for positive latitude.
!
  implicit none

  integer ( kind = 4 ) i
  character              lat_char

  if ( i < 0 ) then
    lat_char = 'S'
  else if ( 0 < i ) then
    lat_char = 'N'
  else
    lat_char = '?'
  end if

  return
end
subroutine ll_degrees_to_dist ( n, lat, long, dist_table )

!*****************************************************************************80
!
!! LL_DEGREES_TO_DIST creates a distance table from latitudes and longitudes.
!
!  Discussion:
!
!    A distance function is used which is appropriate for the earth.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 July 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data items.
!
!    Input, real ( kind = 8 ) LAT(N), LONG(N), the latitude and longitude, 
!    in degrees.
!
!    Output, real ( kind = 8 ) DIST_TABLE(N,N), the distance matrix.  Distances 
!    are measured in miles.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) dist_table(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) lat(n)
  real ( kind = 8 ) long(n)

  do i = 1, n
    dist_table(i,i) = 0.0D+0
    do j = i + 1, n
      call ll_degrees_to_distance_earth ( lat(i), long(i), lat(j), long(j), &
        dist_table(i,j) )
      dist_table(j,i) = dist_table(i,j)
    end do
  end do

  return
end
subroutine ll_degrees_to_distance_earth ( lat1, long1, lat2, long2, dist )

!*****************************************************************************80
!
!! LL_DEGREES_TO_DISTANCE_EARTH: distance between two points on the earth.
!
!  Discussion:
!
!    The positions of the points are given as longitude and
!    latitude, measured in degrees.
!
!    The distance is measured on the surface of the earth, which
!    is approximated by a sphere.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 July 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) LAT1, LON1, the latitude and longitude 
!    of the first point.  
!
!    Input, real ( kind = 8 ) LAT2, LON2, the latitude and longitude 
!    of the second point. 
!
!    Output, real ( kind = 8 ) DIST, the distance between the points, in miles.
!
  implicit none

  real ( kind = 8 ) degrees_to_radians
  real ( kind = 8 ) dist
  real ( kind = 8 ) lat1
  real ( kind = 8 ) lat1_rad
  real ( kind = 8 ) lat2
  real ( kind = 8 ) lat2_rad
  real ( kind = 8 ) long1
  real ( kind = 8 ) long1_rad
  real ( kind = 8 ) long2
  real ( kind = 8 ) long2_rad
  real ( kind = 8 ), parameter :: radius = 3958.89D+00
  real ( kind = 8 ) theta

  lat1_rad = degrees_to_radians ( lat1 )
  long1_rad = degrees_to_radians ( long1 )

  lat2_rad = degrees_to_radians ( lat2 )
  long2_rad = degrees_to_radians ( long2 )

  theta = acos ( sin ( lat1_rad ) * sin ( lat2_rad ) &
               + cos ( lat1_rad ) * cos ( lat2_rad ) &
               * cos ( long1_rad - long2_rad ) )

  dist = radius * theta

  return
end
subroutine ll_degrees_to_xy ( n, lat, long, x, y )

!*****************************************************************************80
!
!! LL_DEGREES_TO_XY: Latitude/Longitude in degrees to XY coordinates.
!
!  Discussion:
!
!    Essentially, the latitude and longitude information is treated
!    as though the Earth were a cylinder.  As long as the the
!    data is relatively close on the sphere (and far from either
!    pole!) the distortion will not be too severe.  If the data
!    is closely clustered, and also near the equator, the
!    positions will be relatively accurate.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 July 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data items.
!
!    Input, real ( kind = 8 ) LAT(N), LONG(N), the latitude and longitude.
!
!    Output, real ( kind = 8 ) X(N), Y(N), the point coordinates, in miles.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) degrees_to_radians
  real ( kind = 8 ) lat(n)
  real ( kind = 8 ) long(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) phi
  real ( kind = 8 ), parameter :: radius = 3958.89D+00
  real ( kind = 8 ) theta
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)

  do i = 1, n
    theta = degrees_to_radians ( long(i) )
    phi = degrees_to_radians ( lat(i) )
    x(i) = radius * theta
    y(i) = radius * phi
  end do

  return
end
subroutine ll_rad_dist_sphere ( lat1, long1, lat2, long2, radius, dist )

!*****************************************************************************80
!
!! LL_RAD_DIST_SPHERE: spherical distance, latitude and longitude in radians.
!
!  Discussion:
!
!    On a sphere of given radius, the positions of two points are given as
!    longitude and latitude, in radians.
!
!    This function determines the spherical distance or great circle distance,
!    between the two points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 June 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) LAT1, LONG1, LAT2, LONG2, the latitude and
!    longitude of the two points, in radians.
!
!    Input, real ( kind = 8 ) RADIUS, the radius of the sphere.
!
!    Output, real ( kind = 8 ) DIST, the distance between the points.
!
  implicit none

  real ( kind = 8 ) dist
  real ( kind = 8 ) lat1
  real ( kind = 8 ) lat2
  real ( kind = 8 ) long1
  real ( kind = 8 ) long2
  real ( kind = 8 ) radius
  real ( kind = 8 ) theta

  theta = acos ( sin ( lat1 ) * sin ( lat2 ) &
               + cos ( lat1 ) * cos ( lat2 ) * cos ( long1 - long2 ) )

  dist = radius * theta

  return
end
function long_char ( i )

!*****************************************************************************80
!
!! LONG_CHAR returns a character for negative or positive longitude.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 November 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, is negative for negative longitude, and 
!    positive for positive longitude.
!
!    Output, character LONG_CHAR, is 'W' for negative longitude, and
!    'E' for positive longitude.
!
  implicit none

  integer ( kind = 4 ) i
  character long_char

  if ( i < 0 ) then
    long_char = 'W'
  else if ( 0 < i ) then
    long_char = 'E'
  else
    long_char = '?'
  end if

  return
end
subroutine main_read_code ( file_main, file_code )

!*****************************************************************************80
!
!! MAIN_READ_CODE reads the name of the code file from the main file.
!
!  Discussion:
!
!    FILE_CODE is the name of a file containing short codes for the 
!    cities.
!
!    There MAY be a record in the main file of the form
!
!    "code  key_code.txt"
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 November 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_MAIN, the name of the file to read.
!
!    Output, character ( len = * ) FILE_CODE, the name of the code file,
!    or ' ' if no information was found.
!
  implicit none

  logical done
  character ( len = * ) file_code
  character ( len = * ) file_main
  integer ( kind = 4 ) input
  integer ( kind = 4 ) ios
  character ( len = 80 ) line
  logical s_eqi
  character ( len = 80 ) word

  call get_unit ( input )

  open ( unit = input, file = file_main, status = 'old' )

  file_code = ' '

  do
!
!  Read the next line from the file.
!
    read ( input, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      exit
    end if
!
!  Skip blank lines and comment lines.
!
    if ( len_trim ( line ) == 0 ) then

    else if ( line(1:1) == '#' ) then

    else

      done = .true.

      call word_next_read ( line, word, done )

      if ( done ) then
        cycle
      end if

      if ( .not. s_eqi ( word, 'code' ) ) then
        cycle
      end if

      call word_next_read ( line, file_code, done )

      exit

    end if

  end do

  close ( unit = input )

  return
end
subroutine main_read_dist ( file_main, file_dist )

!*****************************************************************************80
!
!! MAIN_READ_DIST reads the name of the distance file from the main file.
!
!  Discussion:
!
!    FILE_DIST is the name of a file containing the city-to-city
!    distance matrix.
!
!    There MAY be a record in the main file of the form
!
!    "dist  key_dist.txt"
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 November 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_MAIN, the name of the file to read.
!
!    Output, character ( len = * ) FILE_DIST, the name of the distance file,
!    or ' ' if no information was found.
!
  implicit none

  logical done
  character ( len = * ) file_dist
  character ( len = * ) file_main
  integer ( kind = 4 ) input
  integer ( kind = 4 ) ios
  character ( len = 80 ) line
  logical s_eqi
  character ( len = 80 ) word

  call get_unit ( input )

  open ( unit = input, file = file_main, status = 'old' )

  file_dist = ' '

  do
!
!  Read the next line from the file.
!
    read ( input, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      exit
    end if
!
!  Skip blank lines and comment lines.
!
    if ( len_trim ( line ) == 0 ) then

    else if ( line(1:1) == '#' ) then

    else

      done = .true.

      call word_next_read ( line, word, done )

      if ( done ) then
        cycle
      end if

      if ( .not. s_eqi ( word, 'dist' ) ) then
        cycle
      end if

      call word_next_read ( line, file_dist, done )

      exit

    end if

  end do

  close ( unit = input )

  return
end
subroutine main_read_dms ( file_main, file_dms )

!*****************************************************************************80
!
!! MAIN_READ_DMS reads the name of the DMS file from the main file.
!
!  Discussion:
!
!    FILE_DMS is the name of a file containing the longitude and latitude
!    of each city in degrees/minutes/seconds.
!
!    There MAY be a record in the main file of the form
!
!    "dms  key_dms.txt"
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 November 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_MAIN, the name of the file to read.
!
!    Output, character ( len = * ) FILE_DMS, the name of the DMS file,
!    or ' ' if no information was found.
!
  implicit none

  logical done
  character ( len = * ) file_dms
  character ( len = * ) file_main
  integer ( kind = 4 ) input
  integer ( kind = 4 ) ios
  character ( len = 80 ) line
  logical s_eqi
  character ( len = 80 ) word

  call get_unit ( input )

  open ( unit = input, file = file_main, status = 'old' )

  file_dms = ' '

  do
!
!  Read the next line from the file.
!
    read ( input, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      exit
    end if
!
!  Skip blank lines and comment lines.
!
    if ( len_trim ( line ) == 0 ) then

    else if ( line(1:1) == '#' ) then

    else

      done = .true.

      call word_next_read ( line, word, done )

      if ( done ) then
        cycle
      end if

      if ( .not. s_eqi ( word, 'dms' ) ) then
        cycle
      end if

      call word_next_read ( line, file_dms, done )

      exit

    end if

  end do

  close ( unit = input )

  return
end
subroutine main_read_geom ( file_main, geom )

!*****************************************************************************80
!
!! MAIN_READ_GEOM reads the name of the geometry from the main file.
!
!  Discussion:
!
!    GEOM is the name of the geometry of the city data.
!    Typical values include:
!    none - no special geometry
!    plane - the points lie in a plane
!    sphere - the points lie on a sphere
!
!    There MAY be a record in the main file of the form
!
!    "geom  geom_value"
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 November 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_MAIN, the name of the file to read.
!
!    Output, character ( len = * ) GEOM, the name of the geometry,
!    or ' ' if no information was found.
!
  implicit none

  logical done
  character ( len = * ) file_main
  character ( len = * ) geom
  integer ( kind = 4 ) input
  integer ( kind = 4 ) ios
  character ( len = 80 ) line
  logical s_eqi
  character ( len = 80 ) word

  call get_unit ( input )

  open ( unit = input, file = file_main, status = 'old' )

  geom = ' '

  do
!
!  Read the next line from the file.
!
    read ( input, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      exit
    end if
!
!  Skip blank lines and comment lines.
!
    if ( len_trim ( line ) == 0 ) then

    else if ( line(1:1) == '#' ) then

    else

      done = .true.

      call word_next_read ( line, word, done )

      if ( done ) then
        cycle
      end if

      if ( .not. s_eqi ( word, 'geom' ) ) then
        cycle
      end if

      call word_next_read ( line, geom, done )

      exit

    end if

  end do

  close ( unit = input )

  return
end
subroutine main_read_name ( file_main, file_name )

!*****************************************************************************80
!
!! MAIN_READ_NAME reads the name of the name file from the main file.
!
!  Discussion:
!
!    FILE_NAME is the name of a file containing the city names.
!
!    There MAY be a record in the main file of the form
!
!    "name  key_name.txt"
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 November 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_MAIN, the name of the file to read.
!
!    Output, character ( len = * ) FILE_NAME, the name of the name file,
!    or ' ' if no information was found.
!
  implicit none

  logical done
  character ( len = * ) file_main
  character ( len = * ) file_name
  integer ( kind = 4 ) input
  integer ( kind = 4 ) ios
  character ( len = 80 ) line
  logical s_eqi
  character ( len = 80 ) word

  call get_unit ( input )

  open ( unit = input, file = file_main, status = 'old' )

  file_name = ' '

  do
!
!  Read the next line from the file.
!
    read ( input, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      exit
    end if
!
!  Skip blank lines and comment lines.
!
    if ( len_trim ( line ) == 0 ) then

    else if ( line(1:1) == '#' ) then

    else

      done = .true.

      call word_next_read ( line, word, done )

      if ( done ) then
        cycle
      end if

      if ( .not. s_eqi ( word, 'name' ) ) then
        cycle
      end if

      call word_next_read ( line, file_name, done )

      exit

    end if

  end do

  close ( unit = input )

  return
end
subroutine main_read_size ( file_main, n )

!*****************************************************************************80
!
!! MAIN_READ_SIZE reads the problem size N from the main file.
!
!  Discussion:
!
!    The problem size is N, the number of cities.
!
!    There should always be a record in the main file of the form
!
!    "size  7"
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 November 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_MAIN, the name of the file to read.
!
!    Output, integer ( kind = 4 ) N, the problem size.
!
  implicit none

  logical done
  character ( len = * ) file_main
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) input
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) length
  character ( len = 80 ) line
  integer ( kind = 4 ) n
  logical s_eqi
  character ( len = 80 ) word

  call get_unit ( input )

  open ( unit = input, file = file_main, status = 'old' )

  n = 0

  do
!
!  Read the next line from the file.
!
    read ( input, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      exit
    end if
!
!  Skip blank lines and comment lines.
!
    if ( len_trim ( line ) == 0 ) then

    else if ( line(1:1) == '#' ) then

    else

      done = .true.

      call word_next_read ( line, word, done )

      if ( done ) then
        cycle
      end if

      if ( .not. s_eqi ( word, 'size' ) ) then
        cycle
      end if

      call word_next_read ( line, word, done )

      call s_to_i4 ( word, n, ierror, length )

      exit

    end if

  end do

  close ( unit = input )

  return
end
subroutine main_read_xy ( file_main, file_xy )

!*****************************************************************************80
!
!! MAIN_READ_XY reads the name of the XY file from the main file.
!
!  Discussion:
!
!    FILE_XY is the name of a file containing (X,Y) coordinate data.
!
!    There MAY be a record in the main file of the form
!
!    "xy  key_xy.txt"
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 November 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_MAIN, the name of the file to read.
!
!    Output, character ( len = * ) FILE_XY, the name of the XY file,
!    or ' ' if no information was found.
!
  implicit none

  logical done
  character ( len = * ) file_main
  character ( len = * ) file_xy
  integer ( kind = 4 ) input
  integer ( kind = 4 ) ios
  character ( len = 80 ) line
  logical s_eqi
  character ( len = 80 ) word

  call get_unit ( input )

  open ( unit = input, file = file_main, status = 'old' )

  file_xy = ' '

  do
!
!  Read the next line from the file.
!
    read ( input, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      exit
    end if
!
!  Skip blank lines and comment lines.
!
    if ( len_trim ( line ) == 0 ) then

    else if ( line(1:1) == '#' ) then

    else

      done = .true.

      call word_next_read ( line, word, done )

      if ( done ) then
        cycle
      end if

      if ( .not. s_eqi ( word, 'xy' ) ) then
        cycle
      end if

      call word_next_read ( line, file_xy, done )

      exit

    end if

  end do

  close ( unit = input )

  return
end
subroutine point_to_dist_table ( dim_num, point_num, point, dist_table )

!*****************************************************************************80
!
!! POINT_TO_DIST_TABLE creates a distance table from Cartesian coordinates.
!
!  Discussion:
!
!    The euclidean distance is used.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, real ( kind = 8 ) POINT(DIM_NUM,POINT_NUM), the  point coordinates.
!
!    Output, real ( kind = 8 ) DIST_TABLE(POINT_NUM,POINT_NUM), the 
!    distance table.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) dim
  real ( kind = 8 ) dist_table(point_num,point_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) point(dim_num,point_num)

  do i = 1, point_num
    dist_table(i,i) = 0.0D+0
    do j = i + 1, point_num
      dist_table(i,j) = 0.0D+00
      do dim = 1, dim_num 
        dist_table(i,j) = dist_table(i,j) + ( point(dim,i) - point(dim,j) )**2
      end do
      dist_table(i,j) = sqrt ( dist_table(i,j) )
      dist_table(j,i) = dist_table(i,j)
    end do
  end do

  return
end
subroutine r8mat_data_read ( input_filename, m, n, table )

!*****************************************************************************80
!
!! R8MAT_DATA_READ reads data from an R8MAT file.
!
!  Discussion:
!
!    The file may contain more than N points, but this routine will
!    return after reading N of them.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILENAME, the name of the input file.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Output, real ( kind = 8 ) TABLE(M,N), the table data.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) ierror
  character ( len = * ) input_filename
  integer ( kind = 4 ) input_status
  integer ( kind = 4 ) input_unit
  integer ( kind = 4 ) j
  character ( len = 255 ) line
  real ( kind = 8 ) table(m,n)
  real ( kind = 8 ) x(m)

  ierror = 0

  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_filename, status = 'old', &
    iostat = input_status )

  if ( input_status /= 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_DATA_READ - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the input file "' // &
      trim ( input_filename ) // '" on unit ', input_unit
    stop
  end if

  j = 0

  do while ( j < n )

    read ( input_unit, '(a)', iostat = input_status ) line

    if ( input_status /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8MAT_DATA_READ - Fatal error!'
      write ( *, '(a)' ) '  Error while reading lines of data.'
      write ( *, '(a,i8)' ) '  Number of values expected per line M = ', m
      write ( *, '(a,i8)' ) '  Number of data lines read, J =         ', j
      write ( *, '(a,i8)' ) '  Number of data lines needed, N =       ', n
      stop
    end if

    if ( line(1:1) == '#' .or. len_trim ( line ) == 0 ) then
      cycle
    end if

    call s_to_r8vec ( line, m, x, ierror )

    if ( ierror /= 0 ) then
      cycle
    end if

    j = j + 1

    table(1:m,j) = x(1:m)

  end do

  close ( unit = input_unit )

  return
end
subroutine r8mat_header_read ( input_filename, m, n )

!*****************************************************************************80
!
!! R8MAT_HEADER_READ reads the header from an R8MAT file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILENAME, the name of the input file.
!
!    Output, integer ( kind = 4 ) M, spatial dimension.
!
!    Output, integer ( kind = 4 ) N, the number of points.
!
  implicit none

  character ( len = * ) input_filename
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  call file_column_count ( input_filename, m )

  if ( m <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_HEADER_READ - Fatal error!'
    write ( *, '(a)' ) '  There was some kind of I/O problem while trying'
    write ( *, '(a)' ) '  to count the number of data columns in'
    write ( *, '(a)' ) '  the file "' // trim ( input_filename ) // '".'
    stop
  end if

  call file_row_count ( input_filename, n )

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_HEADER_READ - Fatal error!'
    write ( *, '(a)' ) '  There was some kind of I/O problem while trying'
    write ( *, '(a)' ) '  to count the number of data rows in'
    write ( *, '(a)' ) '  the file "' // trim ( input_filename ) // '".'
    stop
  end if

  return
end
subroutine r8mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8MAT_PRINT prints an R8MAT.
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
!    12 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = * )  title

  call r8mat_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine r8mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8MAT_PRINT_SOME prints some of an R8MAT.
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
!    26 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ( kind = 4 ) ILO, JLO, the first row and column to print.
!
!    Input, integer ( kind = 4 ) IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = 14 ) ctemp(incx)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * )  title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i8,6x)' ) j
    end do

    write ( *, '(''  Col   '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) ' '

    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi

      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( a(i,j) == real ( int ( a(i,j) ), kind = 8 ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) a(i,j)
        else
          write ( ctemp(j2), '(g14.6)' ) a(i,j)
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j), j = 1, inc )

    end do

  end do

  return
end
subroutine r8mat_transpose_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8MAT_TRANSPOSE_PRINT prints an R8MAT, transposed.
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
!    14 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), an M by N matrix to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = * )  title

  call r8mat_transpose_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine r8mat_transpose_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8MAT_TRANSPOSE_PRINT_SOME prints some of an R8MAT, transposed.
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
!    14 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ( kind = 4 ) ILO, JLO, the first row and column to print.
!
!    Input, integer ( kind = 4 ) IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = 14 ) ctemp(incx)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * )  title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  do i2lo = max ( ilo, 1 ), min ( ihi, m ), incx

    i2hi = i2lo + incx - 1
    i2hi = min ( i2hi, m )
    i2hi = min ( i2hi, ihi )

    inc = i2hi + 1 - i2lo

    write ( *, '(a)' ) ' '

    do i = i2lo, i2hi
      i2 = i + 1 - i2lo
      write ( ctemp(i2), '(i8,6x)' ) i
    end do

    write ( *, '(''  Row   '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Col'
    write ( *, '(a)' ) ' '

    j2lo = max ( jlo, 1 )
    j2hi = min ( jhi, n )

    do j = j2lo, j2hi

      do i2 = 1, inc
        i = i2lo - 1 + i2
        write ( ctemp(i2), '(g14.6)' ) a(i,j)
      end do

      write ( *, '(i5,1x,5a14)' ) j, ( ctemp(i), i = 1, inc )

    end do

  end do

  return
end
subroutine r8mat_write ( output_filename, m, n, table )

!*****************************************************************************80
!
!! R8MAT_WRITE writes an R8MAT file.
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
!    Input, real ( kind = 8 ) TABLE(M,N), the table data.
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
!  For greater precision in the output file, try:
!
!                                            '(', m, 'g', 24, '.', 16, ')'
!
  if ( 0 < m .and. 0 < n ) then

    write ( string, '(a1,i8,a1,i8,a1,i8,a1)' ) '(', m, 'g', 14, '.', 6, ')'
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
subroutine r8vec_print ( n, a, title )

!*****************************************************************************80
!
!! R8VEC_PRINT prints a real vector.
!
!  Discussion:
!
!    If all the entries are integer ( kind = 4 )s, the data if printed
!    in integer ( kind = 4 ) format.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 November 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, real ( kind = 8 ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  if ( all ( a(1:n) == aint ( a(1:n) ) ) ) then
    do i = 1, n
      write ( *, '(i6,i6)' ) i, int ( a(i) )
    end do
  else if ( all ( abs ( a(1:n) ) < 1000000.0D+00 ) ) then
    do i = 1, n
      write ( *, '(i6,f14.6)' ) i, a(i)
    end do
  else
    do i = 1, n
      write ( *, '(i6,g14.6)' ) i, a(i)
    end do
  end if

  return
end
subroutine r8vec2_data_read ( input_filename, n, x, y )

!*****************************************************************************80
!
!! R8VEC2_DATA_READ reads data from an R8VEC2 file.
!
!  Discussion:
!
!    An R8VEC2 is a pair of R8VEC's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 July 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILENAME, the name of the input file.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Output, real ( kind = 8 ) X(N), Y(N), the data.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) ierror
  character ( len = * ) input_filename
  integer ( kind = 4 ) input_status
  integer ( kind = 4 ) input_unit
  integer ( kind = 4 ) j
  character ( len = 255 ) line
  real ( kind = 8 ) t(2)
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)

  ierror = 0

  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_filename, status = 'old', &
    iostat = input_status )

  if ( input_status /= 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC2_DATA_READ - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the input file "' // &
      trim ( input_filename ) // '" on unit ', input_unit
    stop
  end if

  j = 0

  do while ( j < n )

    read ( input_unit, '(a)', iostat = input_status ) line

    if ( input_status /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8VEC2_DATA_READ - Fatal error!'
      write ( *, '(a)' ) '  Error while reading lines of data.'
      write ( *, '(a,i8)' ) '  Number of data lines read, J =         ', j
      write ( *, '(a,i8)' ) '  Number of data lines needed, N =       ', n
      stop
    end if

    if ( line(1:1) == '#' .or. len_trim ( line ) == 0 ) then
      cycle
    end if

    call s_to_r8vec ( line, 2, t, ierror )

    if ( ierror /= 0 ) then
      cycle
    end if

    j = j + 1

    x(j) = t(1)
    y(j) = t(2)

  end do

  close ( unit = input_unit )

  return
end
subroutine r8vec2_header_read ( input_filename, n )

!*****************************************************************************80
!
!! R8VEC2_HEADER_READ reads the header from an R8VEC2 file.
!
!  Discussion:
!
!    An R8VEC2 is an pair of R8VEC's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 July 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILENAME, the name of the input file.
!
!    Output, integer ( kind = 4 ) N, the number of points.
!
  implicit none

  character ( len = * ) input_filename
  integer ( kind = 4 ) n

  call file_row_count ( input_filename, n )

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC2_HEADER_READ - Fatal error!'
    write ( *, '(a)' ) '  There was some kind of I/O problem while trying'
    write ( *, '(a)' ) '  to count the number of data rows in'
    write ( *, '(a)' ) '  the file "' // trim ( input_filename ) // '".'
    stop
  end if

  return
end
subroutine r8vec2_print ( n, a1, a2, title )

!*****************************************************************************80
!
!! R8VEC2_PRINT prints an R8VEC2.
!
!  Discussion:
!
!    An R8VEC2 is a dataset consisting of N pairs of R8's, stored
!    as two separate vectors A1 and A2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, real ( kind = 8 ) A1(N), A2(N), the vectors to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a1(n)
  real ( kind = 8 ) a2(n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i4,2x,g14.6,2x,g14.6)' ) i, a1(i), a2(i)
  end do

  return
end
subroutine r8vec2_write ( output_filename, n, x, y )

!*****************************************************************************80
!
!! R8VEC2_WRITE writes an R8VEC2 file.
!
!  Discussion:
!
!    An R8VEC2 is a pair of vectors of R8 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 July 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) OUTPUT_FILENAME, the output file name.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the data.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) j
  character ( len = * ) output_filename
  integer ( kind = 4 ) output_status
  integer ( kind = 4 ) output_unit
  real ( kind = 8 ) x(n)
  real      ( kind = 8 ) y(n)
!
!  Open the file.
!
  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_filename, &
    status = 'replace', iostat = output_status )

  if ( output_status /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_WRITE - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the output file "' // &
      trim ( output_filename ) // '" on unit ', output_unit
    output_unit = -1
    stop
  end if

  if ( 0 < n ) then
!
!  Write the data.
!
    do j = 1, n
      write ( output_unit, '(2x,g24.16,2x,g24.16)' ) x(j), y(j)
    end do

  end if
!
!  Close the file.
!
  close ( unit = output_unit )

  return
end
subroutine s_cat ( s1, s2, s3 )

!*****************************************************************************80
!
!! S_CAT concatenates two strings to make a third string.  
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S1, the "prefix" string.
!
!    Input, character ( len = * ) S2, the "postfix" string.
!
!    Output, character ( len = * ) S3, the string made by
!    concatenating S1 and S2, ignoring any trailing blanks.
!
  implicit none

  character ( len = * ) s1
  character ( len = * ) s2
  character ( len = * ) s3

  if ( s1 == ' ' .and. s2 == ' ' ) then
    s3 = ' '
  else if ( s1 == ' ' ) then
    s3 = s2
  else if ( s2 == ' ' ) then
    s3 = s1
  else
    s3 = trim ( s1 ) // trim ( s2 )
  end if
 
  return
end
function s_eqi ( s1, s2 )

!*****************************************************************************80
!
!! S_EQI is a case insensitive comparison of two strings for equality.
!
!  Example:
!
!    S_EQI ( 'Anjana', 'ANJANA' ) is .TRUE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S1, S2, the strings to compare.
!
!    Output, logical S_EQI, the result of the comparison.
!
  implicit none

  character c1
  character c2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) len1
  integer ( kind = 4 ) len2
  integer ( kind = 4 ) lenc
  logical s_eqi
  character ( len = * )  s1
  character ( len = * )  s2

  len1 = len ( s1 )
  len2 = len ( s2 )
  lenc = min ( len1, len2 )

  s_eqi = .false.

  do i = 1, lenc

    c1 = s1(i:i)
    c2 = s2(i:i)
    call ch_cap ( c1 )
    call ch_cap ( c2 )

    if ( c1 /= c2 ) then
      return
    end if

  end do

  do i = lenc + 1, len1
    if ( s1(i:i) /= ' ' ) then
      return
    end if
  end do

  do i = lenc + 1, len2
    if ( s2(i:i) /= ' ' ) then
      return
    end if
  end do

  s_eqi = .true.

  return
end
subroutine s_rep_one ( s1, sub1, sub2, s2 )

!*****************************************************************************80
!
!! S_REP_ONE replaces the first occurrence of SUB1 with SUB2.
!
!  Discussion:
!
!    The input and output strings may coincide.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 November 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S1, the initial string.
!
!    Input, character ( len = * ) SUB1, the string to be replaced.
!
!    Input, character ( len = * ) SUB2, the replacement string.
!
!    Output, character ( len = * ) S2, the final string.
!
  implicit none

  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) i4
  character ( len = * ) s1
  character ( len = * ) s2
  character ( len = 255 ) s3
  character ( len = * ) sub1
  character ( len = * ) sub2

  s3 = ' '

  i1 = index ( s1, sub1 )

  if ( i1 == 0 ) then

    s3 = s1

  else

    s3(1:i1-1) = s1(1:i1-1)

    i2 = len_trim ( sub2 )
    s3(i1:i1+i2-1) = sub2(1:i2)

    i3 = i1 + len_trim ( sub1 )
    i4 = len_trim ( s1 )

    s3(i1+i2:i1+i2+1+i4-i3) = s1(i3:i4)

  end if

  s2 = s3

  return
end
subroutine s_to_i4 ( s, ival, ierror, last )

!*****************************************************************************80
!
!! S_TO_I4 reads an I4 from a string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, a string to be examined.
!
!    Output, integer ( kind = 4 ) IVAL, the value read from the string.
!    If the string is blank, then IVAL will be returned 0.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!    0, no error.
!    1, an error occurred.
!
!    Output, integer ( kind = 4 ) LAST, the last character of S used.
!
  implicit none

  character c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) istate
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) last
  character ( len = * ) s

  ierror = 0
  istate = 0
  isgn = 1
  ival = 0

  do i = 1, len_trim ( s )

    c = s(i:i)
!
!  Haven't read anything.
!
    if ( istate == 0 ) then

      if ( c == ' ' ) then

      else if ( c == '-' ) then
        istate = 1
        isgn = -1
      else if ( c == '+' ) then
        istate = 1
        isgn = + 1
      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        istate = 2
        ival = ichar ( c ) - ichar ( '0' )
      else
        ierror = 1
        return
      end if
!
!  Have read the sign, expecting digits.
!
    else if ( istate == 1 ) then

      if ( c == ' ' ) then

      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        istate = 2
        ival = ichar ( c ) - ichar ( '0' )
      else
        ierror = 1
        return
      end if
!
!  Have read at least one digit, expecting more.
!
    else if ( istate == 2 ) then

      if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        ival = 10 * ival + ichar ( c ) - ichar ( '0' )
      else
        ival = isgn * ival
        last = i - 1
        return
      end if

    end if

  end do
!
!  If we read all the characters in the string, see if we're OK.
!
  if ( istate == 2 ) then
    ival = isgn * ival
    last = len_trim ( s )
  else
    ierror = 1
    last = 0
  end if

  return
end
subroutine s_to_r8 ( s, dval, ierror, length )

!*****************************************************************************80
!
!! S_TO_R8 reads an R8 from a string.
!
!  Discussion:
!
!    The routine will read as many characters as possible until it reaches
!    the end of the string, or encounters a character which cannot be
!    part of the number.
!
!    Legal input is:
!
!       1 blanks,
!       2 '+' or '-' sign,
!       2.5 blanks
!       3 integer part,
!       4 decimal point,
!       5 fraction part,
!       6 'E' or 'e' or 'D' or 'd', exponent marker,
!       7 exponent sign,
!       8 exponent integer part,
!       9 exponent decimal point,
!      10 exponent fraction part,
!      11 blanks,
!      12 final comma or semicolon,
!
!    with most quantities optional.
!
!  Example:
!
!    S                 DVAL
!
!    '1'               1.0
!    '     1   '       1.0
!    '1A'              1.0
!    '12,34,56'        12.0
!    '  34 7'          34.0
!    '-1E2ABCD'        -100.0
!    '-1X2ABCD'        -1.0
!    ' 2E-1'           0.2
!    '23.45'           23.45
!    '-4.2E+2'         -420.0
!    '17d2'            1700.0
!    '-14e-2'         -0.14
!    'e2'              100.0
!    '-12.73e-9.23'   -12.73 * 10.0**(-9.23)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string containing the
!    data to be read.  Reading will begin at position 1 and
!    terminate at the end of the string, or when no more
!    characters can be read to form a legal real.  Blanks,
!    commas, or other nonnumeric data will, in particular,
!    cause the conversion to halt.
!
!    Output, real ( kind = 8 ) DVAL, the value read from the string.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no errors occurred.
!    1, 2, 6 or 7, the input number was garbled.  The
!    value of IERROR is the last type of input successfully
!    read.  For instance, 1 means initial blanks, 2 means
!    a plus or minus sign, and so on.
!
!    Output, integer ( kind = 4 ) LENGTH, the number of characters read
!    to form the number, including any terminating
!    characters such as a trailing comma or blanks.
!
  implicit none

  character c
  logical ch_eqi
  real ( kind = 8 ) dval
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ihave
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) iterm
  integer ( kind = 4 ) jbot
  integer ( kind = 4 ) jsgn
  integer ( kind = 4 ) jtop
  integer ( kind = 4 ) length
  integer ( kind = 4 ) nchar
  integer ( kind = 4 ) ndig
  real ( kind = 8 ) rbot
  real ( kind = 8 ) rexp
  real ( kind = 8 ) rtop
  character ( len = * )  s

  nchar = len_trim ( s )

  ierror = 0
  dval = 0.0D+00
  length = -1
  isgn = 1
  rtop = 0
  rbot = 1
  jsgn = 1
  jtop = 0
  jbot = 1
  ihave = 1
  iterm = 0

  do

    length = length + 1

    if ( nchar < length+1 ) then
      exit
    end if

    c = s(length+1:length+1)
!
!  Blank character.
!
    if ( c == ' ' ) then

      if ( ihave == 2 ) then

      else if ( ihave == 6 .or. ihave == 7 ) then
        iterm = 1
      else if ( 1 < ihave ) then
        ihave = 11
      end if
!
!  Comma.
!
    else if ( c == ',' .or. c == ';' ) then

      if ( ihave /= 1 ) then
        iterm = 1
        ihave = 12
        length = length + 1
      end if
!
!  Minus sign.
!
    else if ( c == '-' ) then

      if ( ihave == 1 ) then
        ihave = 2
        isgn = -1
      else if ( ihave == 6 ) then
        ihave = 7
        jsgn = -1
      else
        iterm = 1
      end if
!
!  Plus sign.
!
    else if ( c == '+' ) then

      if ( ihave == 1 ) then
        ihave = 2
      else if ( ihave == 6 ) then
        ihave = 7
      else
        iterm = 1
      end if
!
!  Decimal point.
!
    else if ( c == '.' ) then

      if ( ihave < 4 ) then
        ihave = 4
      else if ( 6 <= ihave .and. ihave <= 8 ) then
        ihave = 9
      else
        iterm = 1
      end if
!
!  Scientific notation exponent marker.
!
    else if ( ch_eqi ( c, 'E' ) .or. ch_eqi ( c, 'D' ) ) then

      if ( ihave < 6 ) then
        ihave = 6
      else
        iterm = 1
      end if
!
!  Digit.
!
    else if (  ihave < 11 .and. lle ( '0', c ) .and. lle ( c, '9' ) ) then

      if ( ihave <= 2 ) then
        ihave = 3
      else if ( ihave == 4 ) then
        ihave = 5
      else if ( ihave == 6 .or. ihave == 7 ) then
        ihave = 8
      else if ( ihave == 9 ) then
        ihave = 10
      end if

      call ch_to_digit ( c, ndig )

      if ( ihave == 3 ) then
        rtop = 10.0D+00 * rtop + real ( ndig, kind = 8 )
      else if ( ihave == 5 ) then
        rtop = 10.0D+00 * rtop + real ( ndig, kind = 8 )
        rbot = 10.0D+00 * rbot
      else if ( ihave == 8 ) then
        jtop = 10 * jtop + ndig
      else if ( ihave == 10 ) then
        jtop = 10 * jtop + ndig
        jbot = 10 * jbot
      end if
!
!  Anything else is regarded as a terminator.
!
    else
      iterm = 1
    end if
!
!  If we haven't seen a terminator, and we haven't examined the
!  entire string, go get the next character.
!
    if ( iterm == 1 ) then
      exit
    end if

  end do
!
!  If we haven't seen a terminator, and we have examined the
!  entire string, then we're done, and LENGTH is equal to NCHAR.
!
  if ( iterm /= 1 .and. length+1 == nchar ) then
    length = nchar
  end if
!
!  Number seems to have terminated.  Have we got a legal number?
!  Not if we terminated in states 1, 2, 6 or 7!
!
  if ( ihave == 1 .or. ihave == 2 .or. ihave == 6 .or. ihave == 7 ) then
    ierror = ihave
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'S_TO_R8 - Serious error!'
    write ( *, '(a)' ) '  Illegal or nonnumeric input:'
    write ( *, '(a)' ) '    ' // trim ( s )
    return
  end if
!
!  Number seems OK.  Form it.
!
  if ( jtop == 0 ) then
    rexp = 1.0D+00
  else
    if ( jbot == 1 ) then
      rexp = 10.0D+00 ** ( jsgn * jtop )
    else
      rexp = 10.0D+00 ** ( real ( jsgn * jtop, kind = 8 ) &
        / real ( jbot, kind = 8 ) )
    end if
  end if

  dval = real ( isgn, kind = 8 ) * rexp * rtop / rbot

  return
end
subroutine s_to_r8vec ( s, n, rvec, ierror )

!*****************************************************************************80
!
!! S_TO_R8VEC reads an R8VEC from a string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string to be read.
!
!    Input, integer ( kind = 4 ) N, the number of values expected.
!
!    Output, real ( kind = 8 ) RVEC(N), the values read from the string.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no errors occurred.
!    -K, could not read data for entries -K through N.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) lchar
  real ( kind = 8 ) rvec(n)
  character ( len = * )  s

  i = 0
  ierror = 0
  ilo = 1

  do while ( i < n )

    i = i + 1

    call s_to_r8 ( s(ilo:), rvec(i), ierror, lchar )

    if ( ierror /= 0 ) then
      ierror = -i
      exit
    end if

    ilo = ilo + lchar

  end do

  return
end
subroutine s_to_w ( s, w, ierror, last )

!*****************************************************************************80
!
!! S_TO_W reads the next blank-delimited word from a string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 November 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, a string to be examined.
!
!    Output, character ( len = * ) W, the word that was read.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!    0, no error.
!    1, an error occurred.
!
!    Output, integer ( kind = 4 ) LAST, the last character of S used to make W.
!
  implicit none

  character c
  integer ( kind = 4 ) first
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) istate
  integer ( kind = 4 ) last
  character ( len = * ) s
  character ( len = * ) w

  w = ' '
  ierror = 0
  istate = 0
  first = 0
  last = 0
  i = 0

  do

    i = i + 1

    if ( len_trim ( s ) < i ) then

      if ( istate == 0 ) then
        ierror = 1
        last = 0
      else
        last = i-1
        w = s(first:last)
      end if

      exit

    end if

    c = s(i:i)

    if ( istate == 0 ) then

      if ( c /= ' ' ) then
        first = i
        istate = 1
      end if

    else if ( istate == 1 ) then

      if ( c == ' ' ) then
        last = i - 1
        w = s(first:last)
        exit
      end if

    end if

  end do

  return
end
subroutine s_word_count ( s, nword )

!*****************************************************************************80
!
!! S_WORD_COUNT counts the number of "words" in a string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string to be examined.
!
!    Output, integer ( kind = 4 ) NWORD, the number of "words" in the string.
!    Words are presumed to be separated by one or more blanks.
!
  implicit none

  logical blank
  integer ( kind = 4 ) i
  integer ( kind = 4 ) lens
  integer ( kind = 4 ) nword
  character ( len = * ) s

  nword = 0
  lens = len ( s )

  if ( lens <= 0 ) then
    return
  end if

  blank = .true.

  do i = 1, lens

    if ( s(i:i) == ' ' ) then
      blank = .true.
    else if ( blank ) then
      nword = nword + 1
      blank = .false.
    end if

  end do

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    May 31 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 May 2001
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
  character ( len = 8 ) date
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  character ( len = 10 ) time
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y
  character ( len = 5 ) zone

  call date_and_time ( date, time, zone, values )

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

  write ( *, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
subroutine word_next_read ( s, word, done )

!*****************************************************************************80
!
!! WORD_NEXT_READ "reads" words from a string, one at a time.
!
!  Discussion:
!
!    The following characters are considered to be a single word,
!    whether surrounded by spaces or not:
!
!      " ( ) { } [ ]
!
!    Also, if there is a trailing comma on the word, it is stripped off.
!    This is to facilitate the reading of lists.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, a string, presumably containing words
!    separated by spaces.
!
!    Output, character ( len = * ) WORD.
!
!    If DONE is FALSE, then WORD contains the "next" word read.
!    If DONE is TRUE, then WORD is blank, because there was no more to read.
!
!    Input/output, logical DONE.
!
!    On input with a fresh string, set DONE to TRUE.
!
!    On output, the routine sets DONE:
!      FALSE if another word was read,
!      TRUE if no more words could be read.
!
  implicit none

  logical done
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ), save :: lenc = 0
  integer ( kind = 4 ), save :: next = 1
  character ( len = * ) s
  character, parameter :: TAB = char ( 9 )
  character ( len = * ) word
!
!  We "remember" LENC and NEXT from the previous call.
!
!  An input value of DONE = TRUE signals a new line of text to examine.
!
  if ( done ) then

    next = 1
    done = .false.
    lenc = len_trim ( s )

    if ( lenc <= 0 ) then
      done = .true.
      word = ' '
      return
    end if

  end if
!
!  Beginning at index NEXT, search the string for the next nonblank,
!  which signals the beginning of a word.
!
  ilo = next
!
!  ...S(NEXT:) is blank.  Return with WORD = ' ' and DONE = TRUE.
!
  do

    if ( lenc < ilo ) then
      word = ' '
      done = .true.
      next = lenc + 1
      return
    end if
!
!  If the current character is blank, skip to the next one.
!
    if ( s(ilo:ilo) /= ' ' .and. s(ilo:ilo) /= TAB ) then
      exit
    end if

    ilo = ilo + 1

  end do
!
!  ILO is the index of the next nonblank character in the string.
!
!  If this initial nonblank is a special character,
!  then that's the whole word as far as we're concerned,
!  so return immediately.
!
  if ( s(ilo:ilo) == '"' .or. &
       s(ilo:ilo) == '(' .or. &
       s(ilo:ilo) == ')' .or. &
       s(ilo:ilo) == '{' .or. &
       s(ilo:ilo) == '}' .or. &
       s(ilo:ilo) == '[' .or. &
       s(ilo:ilo) == ']' ) then

    word = s(ilo:ilo)
    next = ilo + 1
    return

  end if
!
!  Now search for the last contiguous character that is not a
!  blank, TAB, or special character.
!
  next = ilo + 1

  do while ( next <= lenc )

    if ( s(next:next) == ' ' ) then
      exit
    else if ( s(next:next) == TAB ) then
      exit
    else if ( s(next:next) == '"' ) then
      exit
    else if ( s(next:next) == '(' ) then
      exit
    else if ( s(next:next) == ')' ) then
      exit
    else if ( s(next:next) == '{' ) then
      exit
    else if ( s(next:next) == '}' ) then
      exit
    else if ( s(next:next) == '[' ) then
      exit
    else if ( s(next:next) == ']' ) then
      exit
    end if

    next = next + 1

  end do
!
!  Ignore a trailing comma.
!
  if ( s(next-1:next-1) == ',' ) then
    word = s(ilo:next-2)
  else
    word = s(ilo:next-1)
  end if

  return
end
