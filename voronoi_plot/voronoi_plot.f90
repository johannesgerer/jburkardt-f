program main

!*****************************************************************************80
!
!! MAIN is the main program for VORONOI_PLOT.
!
!  Discussion:
!
!    VORONOI_PLOT manages the process of plotting a Voronoi diagram.
!
!    This program reads in a 2D pointset contained in the unit square,
!    and makes a PPM ASCII file plot of the Voronoi regions.
!
!    The user has a choice of norms within which to work.
!
!  Usage:
!
!    voronoi_plot point_file norm
!
!    Here "norm" is "1" for the L1 norm, "2" for the L2 norm, and so on,
!    but also we have "I" for the LInfinity norm.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Franz Aurenhammer,
!    Voronoi diagrams - a study of a fundamental geometric data structure,
!    ACM Computing Surveys,
!    Volume 23, pages 345-405, September 1991.
!
!  Parameters:
!
!    Command line argument POINT_FILE, is the name of a file containing
!    the coordinates of a (presumably 2-dimensional) point set to be
!    plotted.
!
!    Command line argument NORM, is a real value which indicates
!    the norm in which the Voronoi diagram should be computed, with
!    acceptable values being 1 for L1, 2 for L2 and so on, but
!    with a special value of 0 meaning L infinity.
!
  implicit none

  integer ( kind = 4 ) arg_num
  real ( kind = 8 ), parameter, dimension ( 2 ) :: box_max = (/ &
    1.0D+00, 1.0D+00 /)
  real ( kind = 8 ), parameter, dimension ( 2 ) :: box_min = (/ &
    0.0D+00, 0.0D+00 /)
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: generator
  character ( len = 255 ) generator_file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iarg
  integer ( kind = 4 ) iargc
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  integer ( kind = 4 ) length
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  character ( len = 10 ) norm
  real ( kind = 8 ) p
  character ( len = 255 ) plot_file_name
  integer ( kind = 4 ), parameter :: plot_col = 300
  integer ( kind = 4 ), parameter :: plot_row = 300

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'VORONOI_PLOT:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Given a file of 2D points and a choice of norm,'
  write ( *, '(a)' ) '  plot the Voronoi neighborhoods.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  (Note that the points must lie in the unit square.)'
!
!  Get the number of command line arguments.
!
  arg_num = iargc ( )
!
!  If at least one command line argument, it's the input file name.
!
  if ( 1 <= arg_num ) then

    iarg = 1
    call getarg ( iarg, generator_file_name )

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter the generator file name.'
    read ( *, '(a)' ) generator_file_name

  end if
!
!  If at least two command line arguments, the second one is the norm.
!
  if ( 2 <= arg_num ) then

    iarg = 2
    call getarg ( iarg, norm )

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  For the LP norm, enter the value of P:'
    write ( *, '(a)' ) '  1 for L1,'
    write ( *, '(a)' ) '  2 for L2 and so on, but'
    write ( *, '(a)' ) '  0 means "LInfinity"'
    write ( *, '(a)' ) '  (P can be any positive real value.)'

    read ( *, '(a)' ) norm

  end if

  call s_low ( norm )

  if ( norm(1:1) == 'L' .or. norm(1:1) == 'l' ) then
    norm(1:9) = norm(2:10)
    norm(10:10) = ' '
  end if

  if ( norm(1:1) == 'I' .or. norm(1:1) == 'i' ) then
    p = 0.0D+00
    norm = 'i'
  else
    call s_to_r8 ( norm, p, ierror, length )
  end if
!
!  P < 0 is no good.
!
  if ( p < 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'VORONOI_PLOT: - Fatal error!'
    write ( *, '(a)' ) '  The norm P must be positive.'
    write ( *, '(a,g14.6)' ) '  Input value P = ', p
    stop
  end if
!
!  Measure the size of the data, allocate it, and read it.
!
  call r8mat_header_read ( generator_file_name, m, n )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Spatial dimension M = ', m
  write ( *, '(a,i6)' ) '  Number of points  N = ', n

  allocate ( generator(1:m,1:n) )

  call r8mat_data_read ( generator_file_name, m, n, generator )
!
!  Write the PPMA file.
!
  call file_name_ext_get ( generator_file_name, i, j )

  if ( p == 0 ) then
    norm = 'i'
  end if

  plot_file_name = generator_file_name(1:i-1) // '_l' // &
    trim ( norm ) // '.ppma'

  call region_plot_ppma ( plot_file_name, plot_row, plot_col, m, &
    n, box_min, box_max, generator, p )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Created the graphics file "' &
    // trim ( plot_file_name ) // '".'

  deallocate ( generator )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'VORONOI_PLOT:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
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
!    Input, real ( kind = 8 ) ANGLE, the angle in the color hexagon.  
!    The sextants are defined by the following points:
!        0 degrees, 1, 0, 0, red;
!       60 degrees, 1, 1, 0, yellow;
!      120 degrees, 0, 1, 0, green;
!      180 degrees, 0, 1, 1, cyan;
!      240 degrees, 0, 0, 1, blue;
!      300 degrees, 1, 0, 1, magenta.
!
!    Output, real ( kind = 8 ) R, G, B, RGB specifications for the color 
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

  logical ch_eqi
  character c1
  character c1_cap
  character c2
  character c2_cap

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
subroutine ch_low ( c )

!*****************************************************************************80
!
!! CH_LOW lowercases a single character.
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
!    Input/output, character C, the character to be lowercased.
!
  implicit none

  character c
  integer ( kind = 4 ) i

  i = ichar ( c )

  if ( 65 <= i .and. i <= 90 ) then
    c = char ( i + 32 )
  end if

  return
end
subroutine ch_to_digit ( c, digit )

!*****************************************************************************80
!
!! CH_TO_DIGIT returns the integer value of a base 10 digit.
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
!    Output, integer ( kind = 4 ) DIGIT, the corresponding integer value.  
!    If C was 'illegal', then DIGIT is -1.
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
subroutine digit_to_ch ( digit, c )

!*****************************************************************************80
!
!! DIGIT_TO_CH returns the character representation of a decimal digit.
!
!  Example:
!
!    DIGIT   C
!    -----  ---
!      0    '0'
!      1    '1'
!    ...    ...
!      9    '9'
!     17    '*'
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
!    Input, integer ( kind = 4 ) DIGIT, the digit value between 0 and 9.
!
!    Output, character C, the corresponding character, or '*' if DIGIT
!    was illegal.
!
  implicit none

  character c
  integer ( kind = 4 ) digit

  if ( 0 <= digit .and. digit <= 9 ) then

    c = char ( digit + 48 )

  else

    c = '*'

  end if

  return
end
subroutine file_column_count ( file_in_name, ncolumn )

!*****************************************************************************80
!
!! FILE_COLUMN_COUNT counts the number of columns in the first line of a file.
!
!  Discussion:
!
!    The file is assumed to be a simple text file.
!
!    Most lines of the file is presumed to consist of NCOLUMN words, separated
!    by spaces.  There may also be some blank lines, and some comment lines,
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
!    Input, character ( len = * ) FILE_IN_NAME, the name of the file.
!
!    Output, integer ( kind = 4 ) NCOLUMN, the number of columns assumed to be in the file.
!
  implicit none

  character ( len = * ) file_in_name
  integer ( kind = 4 ) file_in_unit
  logical got_one
  integer ( kind = 4 ) ios
  character ( len = 255 ) line
  integer ( kind = 4 ) ncolumn
!
!  Open the file.
!
  call get_unit ( file_in_unit )

  open ( unit = file_in_unit, file = file_in_name, status = 'old', &
    form = 'formatted', access = 'sequential', iostat = ios )

  if ( ios /= 0 ) then
    ncolumn = - 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_COLUMN_COUNT - Fatal error!'
    write ( *, '(a)' ) '  Could not open the file:'
    write ( *, '(a)' ) '    ' // trim ( file_in_name )
    return
  end if
!
!  Read one line, but skip blank lines and comment lines.
!
  got_one = .false.

  do

    read ( file_in_unit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
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

    rewind ( file_in_unit )

    do

      read ( file_in_unit, '(a)', iostat = ios ) line

      if ( ios /= 0 ) then
        exit
      end if

      if ( len_trim ( line ) == 0 ) then
        cycle
      end if

      got_one = .true.
      exit

    end do

  end if

  close ( unit = file_in_unit )

  if ( .not. got_one ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_COLUMN_COUNT - Warning!'
    write ( *, '(a)' ) '  The file does not seem to contain any data.'
    ncolumn = 0
    return
  end if

  call s_word_count ( line, ncolumn )

  return
end
subroutine file_name_ext_get ( file_name, i, j )

!*****************************************************************************80
!
!! FILE_NAME_EXT_GET determines the "extension" of a file name.
!
!  Discussion:
!
!    The "extension" of a filename is the string of characters
!    that appears after the LAST period in the name.  A file
!    with no period, or with a period as the last character
!    in the name, has a "null" extension.
!
!    Blanks are unusual in filenames.  This routine ignores all
!    trailing blanks, but will treat initial or internal blanks
!    as regular characters acceptable in a file name.
!
!  Example:
!
!    FILE_NAME   I  J
!
!    bob.for     4  7
!    N.B.C.D     6  7
!    Naomi.      6  6
!    Arthur      0  0
!    .com        1  1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, a file name to be examined.
!
!    Output, integer ( kind = 4 ) I, J, the indices of the first and last characters
!    in the file extension.
!
!    If no period occurs in FILE_NAME, then
!      I = J = 0;
!    Otherwise,
!      I is the position of the LAST period in FILE_NAME, and J is the
!      position of the last nonblank character following the period.
!
  implicit none

  character ( len = * ) file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) s_index_last

  i = s_index_last ( file_name, '.' )

  if ( i /= 0 ) then

    j = len_trim ( file_name )

  else

    j = 0

  end if

  return
end
subroutine file_row_count ( file_in_name, row_num )

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
!    Input, character ( len = * ) FILE_IN_NAME, the name of the input file.
!
!    Output, integer ( kind = 4 ) ROW_NUM, the number of rows found.
!
  implicit none

  integer ( kind = 4 ) bad_num
  integer ( kind = 4 ) comment_num
  character ( len = * ) file_in_name
  integer ( kind = 4 ) file_in_unit
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  character ( len = 255 ) line
  integer ( kind = 4 ) record_num
  integer ( kind = 4 ) row_num

  call get_unit ( file_in_unit )

  open ( unit = file_in_unit, file = file_in_name, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_ROW_COUNT - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input file: ' // &
      trim ( file_in_name )
    stop
  end if

  comment_num = 0
  row_num = 0
  record_num = 0
  bad_num = 0

  do

    read ( file_in_unit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
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

  close ( unit = file_in_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FILE_ROW_COUNT:'
  write ( *, '(a,i6)' ) '  Number of records:         ', record_num
  write ( *, '(a,i6)' ) '  Number of data records:    ', row_num
  write ( *, '(a,i6)' ) '  Number of comment records: ', comment_num

  return
end
subroutine find_closest_norm ( m, n, x, generator, p, nearest )

!*****************************************************************************80
!
!! FIND_CLOSEST_NORM finds the generator closest to a point X in a given norm.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 June 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of cell generatorrs.
!
!    Input, real ( kind = 8 ) X(M), the point to be checked.
!
!    Input, real ( kind = 8 ) GENERATOR(M,N), the cell generators.
!
!    Input, real ( kind = 8 ) P, is the value of P in the LP norm.
!    The special value P = 0 is interpreted as a request for the 
!    L infinity norm.
!
!    Output, integer ( kind = 4 ) NEAREST, the index of the nearest 
!    cell generators.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) dist
  real ( kind = 8 ) dist_min
  real ( kind = 8 ) generator(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) nearest
  real ( kind = 8 ) p
  real ( kind = 8 ) p_inverse
  real ( kind = 8 ) x(m)

  nearest = -1
  dist_min = huge ( dist_min )

  if ( p /= 0.0D+00 ) then
    p_inverse = 1.0D+00 / p
  end if

  do i = 1, n

    if ( p == 0 ) then
      dist = maxval ( abs ( generator(1:m,i) - x(1:m) ) )
    else if ( p == 1 ) then
      dist = sum ( abs ( generator(1:m,i) - x(1:m) ) )
    else if ( p == 2 ) then
      dist = sqrt ( sum ( ( generator(1:m,i) - x(1:m) )**2 ) )
    else
      dist = ( sum ( ( abs ( generator(1:m,i) - x(1:m) ) )**p ) )**p_inverse
    end if

    if ( dist < dist_min ) then
      dist_min = dist
      nearest = i
    end if

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
!    Input, integer ( kind = 4 ) I, the number whose logarithm base 2 is desired.
!
!    Output, integer ( kind = 4 ) I4_LOG_2, the integer part of the logarithm 
!    base 2 of the absolute value of I.
!    For positive I4_LOG_2(I), it should be true that
!      2**I4_LOG_2(X) <= |I| < 2**(I4_LOG_2(I)+1).
!    The special case of I4_LOG_2(0) returns -HUGE().
!
  implicit none

  integer ( kind = 4 ) i4_log_2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_abs

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
!    Input, integer ( kind = 4 ) I, the index of the desired color.
!
!    Output, real ( kind = 8 ) ANGLE, an angle, measured in degrees, 
!    between 0 and 360.
!
  implicit none

  real ( kind = 8 ) angle
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_log_2
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) i4

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
subroutine i4_to_rgb ( i, r, g, b )

!*****************************************************************************80
!
!! I4_TO_RGB maps integers to RGB colors.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 January 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the index of the desired color.
!
!    Output, integer ( kind = 4 ) R, G, B, the RGB specifications for a color. 
!
  implicit none

  real ( kind = 8 ) angle
  integer ( kind = 4 ) b
  integer ( kind = 4 ) g
  integer ( kind = 4 ) i
  integer ( kind = 4 ) r
  real ( kind = 8 ) rb
  real ( kind = 8 ) rg
  real ( kind = 8 ) rr
!
!  Red
!
  if ( i == 1 ) then
    r = 255
    g = 0
    b = 0
!
!  Green
!
  else if ( i == 2 ) then
    r = 0
    g = 255
    b = 0
!
!  Blue
!
  else if ( i == 3 ) then
    r = 0
    g = 0
    b = 255
!
!  Cyan
!
  else if ( i == 4 ) then
    r = 0
    g = 255
    b = 255
!
!  Magenta
!
  else if ( i == 5 ) then
    r = 255
    g = 0
    b = 255
!
!  Yellow
!
  else if ( i == 6 ) then
    r = 255
    g = 255
    b = 0
!
!  Brown5
!
  else if ( i == 7 ) then
    r = 139
    g =  35
    b =  35
!
!  Orange
!
  else if ( i == 8 ) then
    r = 255
    g = 165
    b = 0
!
!  Goldenrod5
!
  else if ( i == 9 ) then
    r = 139
    g = 105
    b =  20
!
!  Medium Purple
!
  else if ( i == 10 ) then
    r = 147
    g = 112
    b = 219
!
!  Coral
!
   else if ( i == 11 ) then
    
    r = 255
    g = 127
    b =  80
!
!  Pink5
!
  else if ( i == 12 ) then
    r = 139
    g =  99
    b = 108
!
!  GreenYellow
!
  else if ( i == 13 ) then
    r = 173
    g = 255
    b =  47
!
!  Aquamarine
!
  else if ( i == 14 ) then
    r = 127
    g = 255
    b = 212
!
!  Pale Green3
!
  else if ( i == 15 ) then
    r = 124
    g = 205
    b = 124
!
!  Burlywood
!
  else if ( i == 16 ) then
    r = 222
    g = 184
    b = 135
!
!  Cornsilk3
!
  else if ( i == 17 ) then
    r = 205
    g = 200
    b = 177
!
!  Lemon_Chiffon3
!
  else if ( i == 18 ) then
    r = 205
    g = 201
    b = 165
!
!  Maroon
!
  else if ( i == 19 ) then
    r = 176
    g = 48
    b = 96
!
!  Slate_Blue2
!
  else if ( i == 20 ) then
    r = 131
    g = 111
    b = 255

  else

    call i4_to_angle ( i, angle )

    call angle_to_rgb ( angle, rr, rg, rb )

    r = min ( int ( rr * 255 ), 255 )
    g = min ( int ( rg * 255 ), 255 )
    b = min ( int ( rb * 255 ), 255 )

  end if

  return
end
subroutine i4_to_s_left ( intval, s )

!*****************************************************************************80
!
!! I4_TO_S_LEFT converts an integer to a left-justified string.
!
!  Example:
!
!    Assume that S is 6 characters long:
!
!    INTVAL  S
!
!         1  1
!        -1  -1
!         0  0
!      1952  1952
!    123456  123456
!   1234567  ******  <-- Not enough room!
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
!    Input, integer ( kind = 4 ) INTVAL, an integer to be converted.
!
!    Output, character ( len = * ) S, the representation of the integer.
!    The integer will be left-justified.  If there is not enough space,
!    the string will be filled with stars.
!
  implicit none

  character c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) idig
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) intval
  integer ( kind = 4 ) ipos
  integer ( kind = 4 ) ival
  character ( len = * ) s

  s = ' '

  ilo = 1
  ihi = len ( s )

  if ( ihi <= 0 ) then
    return
  end if
!
!  Make a copy of the integer.
!
  ival = intval
!
!  Handle the negative sign.
!
  if ( ival < 0 ) then

    if ( ihi <= 1 ) then
      s(1:1) = '*'
      return
    end if

    ival = - ival
    s(1:1) = '-'
    ilo = 2

  end if
!
!  The absolute value of the integer goes into S(ILO:IHI).
!
  ipos = ihi
!
!  Find the last digit of IVAL, strip it off, and stick it into the string.
!
  do

    idig = mod ( ival, 10 )
    ival = ival / 10

    if ( ipos < ilo ) then
      do i = 1, ihi
        s(i:i) = '*'
      end do
      return
    end if

    call digit_to_ch ( idig, c )

    s(ipos:ipos) = c
    ipos = ipos - 1

    if ( ival == 0 ) then
      exit
    end if

  end do
!
!  Shift the string to the left.
!
  s(ilo:ilo+ihi-ipos-1) = s(ipos+1:ihi)
  s(ilo+ihi-ipos:ihi) = ' '

  return
end
subroutine ppma_write ( file_out_name, nrow, ncol, r, g, b, ierror )

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
!    Input, integer ( kind = 4 ) NROW, NCOL, the number of rows and columns of data.
!
!    Input, integer ( kind = 4 ) R(NROW,NCOL), G(NROW,NCOL), B(NROW,NCOL), contain
!    the red, green and blue values of each pixel.  These should
!    be positive.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!    0, no error.
!    1, the data was illegal.
!    2, the file could not be opened.
!
  implicit none

  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nrow

  integer ( kind = 4 ) b(nrow,ncol)
  logical, parameter :: debug = .false.
  character ( len = * ) file_out_name
  integer ( kind = 4 ) file_out_unit
  integer ( kind = 4 ) g(nrow,ncol)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) r(nrow,ncol)
  integer ( kind = 4 ) rgb_max

  ierror = 0
!
!  Compute the maximum color value.
!
  rgb_max = max ( &
    maxval ( r(1:nrow,1:ncol) ), &
    maxval ( g(1:nrow,1:ncol) ), &
    maxval ( b(1:nrow,1:ncol) ) )
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
  call ppma_write_header ( file_out_name, file_out_unit, nrow, ncol, &
    rgb_max, ierror )
!
!  Write the data.
!
  call ppma_write_data ( file_out_unit, nrow, ncol, r, g, b, ierror )
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
    write ( *, '(a,i6)' ) '  Number of data rows NROW =    ', nrow
    write ( *, '(a,i6)' ) '  Number of data columns NCOL = ', ncol
    write ( *, '(a,i6)' ) '  Maximum RGB value RGB_MAX =   ', rgb_max
  end if

  return
end
subroutine ppma_write_data ( file_out_unit, nrow, ncol, r, g, b, ierror )

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
!    Input, integer ( kind = 4 ) NROW, NCOL, the number of rows and columns of data.
!
!    Input, integer ( kind = 4 ) R(NROW,NCOL), G(NROW,NCOL), B(NROW,NCOL), contain
!    the red, green and blue values of each pixel.  These should
!    be positive.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!    0, no error.
!    1, the data was illegal.
!    2, the file could not be opened.
!
  implicit none

  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nrow

  integer ( kind = 4 ) b(nrow,ncol)
  integer ( kind = 4 ) file_out_unit
  integer ( kind = 4 ) g(nrow,ncol)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) r(nrow,ncol)

  ierror = 0
!
!  Write the header.
!
  do i = 1, nrow
    do jlo = 1, ncol, 4
      jhi = min ( jlo + 3, ncol )
      write ( file_out_unit, '(12i5)' ) ( r(i,j), g(i,j), b(i,j), j = jlo, jhi )
    end do
  end do

  return
end
subroutine ppma_write_header ( file_out_name, file_out_unit, nrow, ncol, &
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
!    Input, integer ( kind = 4 ) NROW, NCOL, the number of rows and columns of data.
!
!    Input, integer ( kind = 4 ) RGB_MAX, the maximum value of any data component.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!    0, no error.
!    1, the data was illegal.
!    2, the file could not be opened.
!
  implicit none

  character ( len = * ) file_out_name
  integer ( kind = 4 ) file_out_unit
  integer ( kind = 4 ) ierror
  character ( len = 2 ) :: magic = 'P3'
  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nrow
  integer ( kind = 4 ) rgb_max

  ierror = 0
!
!  Write the header.
!
  write ( file_out_unit, '(a2)' ) magic
  write ( file_out_unit, '(a)' ) '# ' // trim ( file_out_name ) &
    // ' created by PPMA_IO::PPMA_WRITE.F90.'
  write ( file_out_unit, '(i5,2x,i5)' ) ncol, nrow
  write ( file_out_unit, '(i5)' ) rgb_max

  return
end
function r82vec_dist_l2 ( a1, a2 )

!*****************************************************************************80
!
!! R82VEC_DIST_L2 returns the L2 distance between a pair of R82VEC's.
!
!  Discussion:
!
!    The vector L2 norm is defined as:
!
!      value = sqrt ( sum ( 1 <= I <= N ) A(I)**2 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 December 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A1(2), A2(2), the vectors whose L2 distance
!    is desired.
!
!    Output, real ( kind = 8 ) R82VEC_DIST_L2, the L2 norm of the distance
!    between A1 and A2.
!
  implicit none

  real ( kind = 8 ) a1(2)
  real ( kind = 8 ) a2(2)
  real ( kind = 8 ) r82vec_dist_l2

  r82vec_dist_l2 = sqrt ( sum ( ( a1(1:2) - a2(1:2) )**2 ) )

  return
end
subroutine r8mat_data_read ( input_filename, m, n, table )

!*****************************************************************************80
!
!! R8MAT_DATA_READ reads data from an R8MAT file.
!
!  Discussion:
!
!    An R8MAT is an array of R8 values.
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
subroutine r8mat_transpose_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8MAT_TRANSPOSE_PRINT prints an R8MAT, transposed.
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
  character ( len = * ) title

  call r8mat_transpose_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine r8mat_transpose_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8MAT_TRANSPOSE_PRINT_SOME prints some of an R8MAT, transposed.
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
  character ( len = * ) title

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
      write ( ctemp(i2), '(i7,7x)') i
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
subroutine region_plot_ppma ( plot_file_name, plot_row, plot_col, m, &
  n, box_min, box_max, generator, p )

!*****************************************************************************80
!
!! REGION_PLOT_PPMA makes a ASCII PPM plot of the CVT regions in the LP norm.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 June 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters
!
!    Input, character ( len = * ) PLOT_FILE_NAME, the name of the file
!    to be created.
!
!    Input, integer ( kind = 4 ) PLOT_ROW, PLOT_COL, the number of pixels per 
!    row and column of the image.
!
!    Input, integer ( kind = 4 ) M, the dimension of the space.
!
!    Input, integer ( kind = 4 ) N, the number of Voronoi cells.
!
!    Input, real ( kind = 8 ) BOX_MIN(2), BOX_MAX(2), the coordinates
!    of the two extreme corners of the bounding box.
!
!    Input, real ( kind = 8 ) GENERATOR(M,N), the Voronoi cell generators.
!
!    Input, real ( kind = 8 ) P, indicates the LP norm being used, with
!    1 = L1, 2 = L2 and so on, but a special value 0 = Linfinity.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) plot_col
  integer ( kind = 4 ) plot_row

  integer ( kind = 4 ) b(plot_row,plot_col)
  real ( kind = 8 ) box_max(2)
  real ( kind = 8 ) box_min(2)
  integer ( kind = 4 ) g(plot_row,plot_col)
  real ( kind = 8 ) generator(m,n)
  real ( kind = 8 ) generator_dist
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nearest
  real ( kind = 8 ) p
  character ( len = * ) plot_file_name
  real ( kind = 8 ) psize
  real ( kind = 8 ) pxsize
  real ( kind = 8 ) pysize
  integer ( kind = 4 ) r(plot_row,plot_col)
  real ( kind = 8 ) r82vec_dist_l2
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) xy(2)
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin
!
!  Define the size of the page.
!
  xmin = box_min(1)
  ymin = box_min(2)

  xmax = box_max(1)
  ymax = box_max(2)

  pysize = ( ymax - ymin ) / real ( plot_row - 1, kind = 8 )
  pxsize = ( xmax - xmin ) / real ( plot_col - 1, kind = 8 )
  psize = min ( pysize, pxsize )

  do i = 1, plot_row

    xy(2) = ( real ( plot_row - i,     kind = 8 ) * ymax   &
            + real (            i - 1, kind = 8 ) * ymin ) &
            / real ( plot_row     - 1, kind = 8 )

    do j = 1, plot_col

      xy(1) = ( real ( plot_col - j,     kind = 8 ) * xmin   &
              + real (            j - 1, kind = 8 ) * xmax ) &
              / real ( plot_col     - 1, kind = 8 )

      call find_closest_norm ( m, n, xy, generator, p, nearest )
       
      generator_dist = r82vec_dist_l2 ( xy, generator(1:2,nearest) )

      if ( generator_dist <= 7.0D+00 * psize ) then
        r(i,j) = 192
        g(i,j) = 192
        b(i,j) = 192
      else
        call i4_to_rgb ( nearest, r(i,j), g(i,j), b(i,j) )
      end if

    end do

  end do
!
!  Now write all this data to a file.
!
  call ppma_write ( plot_file_name, plot_row, plot_col, r, g, b, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'REGION_PLOT_PPMA - Fatal error!'
    write ( *, '(a,i6)' ) '  PPMA_WRITE returns IERROR = ', ierror
  end if

  return
end
function s_index_last ( s, sub )

!*****************************************************************************80
!
!! S_INDEX_LAST finds the LAST occurrence of a given substring.
!
!  Discussion:
!
!    It returns the location in the string at which the substring SUB is
!    first found, or 0 if the substring does not occur at all.
!
!    The routine is also trailing blank insensitive.  This is very
!    important for those cases where you have stored information in
!    larger variables.  If S is of length 80, and SUB is of
!    length 80, then if S = 'FRED' and SUB = 'RED', a match would
!    not be reported by the standard FORTRAN INDEX, because it treats
!    both variables as being 80 characters long!  This routine assumes that
!    trailing blanks represent garbage!
!
!    This means that this routine cannot be used to find, say, the last
!    occurrence of a substring 'A ', since it assumes the blank space
!    was not specified by the user, but is, rather, padding by the
!    system.  However, as a special case, this routine can properly handle
!    the case where either S or SUB is all blanks.
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
!    Input, character ( len = * ) S, the string to be searched.
!
!    Input, character ( len = * ) SUB, the substring to search for.
!
!    Output, integer ( kind = 4 ) S_INDEX_LAST.  0 if SUB does not occur in
!    the string.  Otherwise S_INDEX_LAST = I, where S(I:I+LENS-1) = SUB,
!    where LENS is the length of SUB, and is the last place
!    this happens.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) llen1
  integer ( kind = 4 ) llen2
  character ( len = * ) s
  integer ( kind = 4 ) s_index_last
  character ( len = * ) sub

  s_index_last = 0

  llen1 = len_trim ( s )
  llen2 = len_trim ( sub )
!
!  In case S or SUB is blanks, use LEN
!
  if ( llen1 == 0 ) then
    llen1 = len ( s )
  end if

  if ( llen2 == 0 ) then
    llen2 = len ( sub )
  end if

  if ( llen1 < llen2 ) then
    return
  end if

  do j = 1, llen1+1-llen2

    i = llen1 + 2 - llen2 - j

    if ( s(i:i+llen2-1) == sub ) then
      s_index_last = i
      return
    end if

  end do

  return
end
subroutine s_low ( s )

!*****************************************************************************80
!
!! S_LOW replaces all uppercase letters by lowercase ones.
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
!    Input/output, character ( len = * ) S, the string to be
!    transformed.  On output, the string is all lowercase.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) nchar
  character ( len = * ) s

  nchar = len_trim ( s )

  do i = 1, nchar
    call ch_low ( s(i:i) )
  end do

  return
end
subroutine s_to_i4 ( s, ival, ierror, length )

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
!    Output, integer ( kind = 4 ) IVAL, the integer value read from the string.
!    If the string is blank, then IVAL will be returned 0.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!    0, no error.
!    1, an error occurred.
!
!    Output, integer ( kind = 4 ) LENGTH, the number of characters of S used.
!
  implicit none

  character c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) istate
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) length
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
        length = i - 1
        return
      end if

    end if

  end do
!
!  If we read all the characters in the string, see if we're OK.
!
  if ( istate == 2 ) then
    ival = isgn * ival
    length = len_trim ( s )
  else
    ierror = 1
    length = 0
  end if

  return
end
subroutine s_to_r8 ( s, r, ierror, lchar )

!*****************************************************************************80
!
!! S_TO_R8 reads an R8 from a string.
!
!  Discussion:
!
!    This routine will read as many characters as possible until it reaches
!    the end of the string, or encounters a character which cannot be
!    part of the real number.
!
!    Legal input is:
!
!       1 blanks,
!       2 '+' or '-' sign,
!       2.5 spaces
!       3 integer part,
!       4 decimal point,
!       5 fraction part,
!       6 'E' or 'e' or 'D' or 'd', exponent marker,
!       7 exponent sign,
!       8 exponent integer part,
!       9 exponent decimal point,
!      10 exponent fraction part,
!      11 blanks,
!      12 final comma or semicolon.
!
!    with most quantities optional.
!
!  Example:
!
!    S                 R
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
!    12 February 2001
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
!    Output, real ( kind = 8 ) R, the real value that was read from the string.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!
!    0, no errors occurred.
!
!    1, 2, 6 or 7, the input number was garbled.  The
!    value of IERROR is the last type of input successfully
!    read.  For instance, 1 means initial blanks, 2 means
!    a plus or minus sign, and so on.
!
!    Output, integer ( kind = 4 ) LCHAR, the number of characters read from
!    the string to form the number, including any terminating
!    characters such as a trailing comma or blanks.
!
  implicit none

  logical ch_eqi
  character c
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ihave
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) iterm
  integer ( kind = 4 ) jbot
  integer ( kind = 4 ) jsgn
  integer ( kind = 4 ) jtop
  integer ( kind = 4 ) lchar
  integer ( kind = 4 ) nchar
  integer ( kind = 4 ) ndig
  real ( kind = 8 ) r
  real ( kind = 8 ) rbot
  real ( kind = 8 ) rexp
  real ( kind = 8 ) rtop
  character ( len = * ) s
  character, parameter :: TAB = char ( 9 )

  nchar = len_trim ( s )
  ierror = 0
  r = 0.0D+00
  lchar = - 1
  isgn = 1
  rtop = 0.0D+00
  rbot = 1.0D+00
  jsgn = 1
  jtop = 0
  jbot = 1
  ihave = 1
  iterm = 0

  do

    lchar = lchar + 1
    c = s(lchar+1:lchar+1)
!
!  Blank or TAB character.
!
    if ( c == ' ' .or. c == TAB ) then

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
        lchar = lchar + 1
      end if
!
!  Minus sign.
!
    else if ( c == '-' ) then

      if ( ihave == 1 ) then
        ihave = 2
        isgn = - 1
      else if ( ihave == 6 ) then
        ihave = 7
        jsgn = - 1
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
!  Exponent marker.
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
    else if ( ihave < 11 .and. lge ( c, '0' ) .and. lle ( c, '9' ) ) then

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
    if ( iterm == 1 .or. nchar <= lchar+1 ) then
      exit
    end if

  end do
!
!  If we haven't seen a terminator, and we have examined the
!  entire string, then we're done, and LCHAR is equal to NCHAR.
!
  if ( iterm /= 1 .and. lchar+1 == nchar ) then
    lchar = nchar
  end if
!
!  Number seems to have terminated.  Have we got a legal number?
!  Not if we terminated in states 1, 2, 6 or 7!
!
  if ( ihave == 1 .or. ihave == 2 .or. ihave == 6 .or. ihave == 7 ) then

    ierror = ihave

    return
  end if
!
!  Number seems OK.  Form it.
!
  if ( jtop == 0 ) then
    rexp = 1.0D+00
  else

    if ( jbot == 1 ) then
      rexp = 10.0D+00**( jsgn * jtop )
    else
      rexp = jsgn * jtop
      rexp = rexp / jbot
      rexp = 10.0D+00**rexp
    end if

  end if

  r = isgn * rexp * rtop / rbot

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
!    19 February 2001
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
  character ( len = * ) s

  i = 0

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
!    26 February 2005
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

  write ( *, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
