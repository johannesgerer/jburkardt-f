program main

!*****************************************************************************80
!
!! MAIN is the main program for CIRCLE_TEST.
!
!  Discussion:
!
!    CIRCLE_TEST reads a point dataset, and applies the circle test.
!
!  Usage:
!
!    circle_test point_file
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 April 2003
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) circle_volume
  integer ( kind = 4 ) dim_num
  character ( len = 100 ) input_file_name
  integer ( kind = 4 ) point_num
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: quasi
  logical walls

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CIRCLE_TEST'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Perform the circle packing test.'

  call get_input_file_name ( input_file_name )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The data will be read from "' &
    // trim ( input_file_name ) // '".'

  call file_column_count ( input_file_name, dim_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  The spatial dimension is ', dim_num

  if ( dim_num <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CIRCLE_TEST_MAIN - Fatal error!'
    write ( *, '(a)' ) '  The spatial dimension must be at least 1!'
    stop
  end if

  call file_line_count ( input_file_name, point_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  The number of points is ', point_num

  if ( point_num <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CIRCLE_TEST_MAIN - Fatal error!'
    write ( *, '(a)' ) '  The number of points must be at least 1!'
    stop
  end if

  allocate ( quasi(1:dim_num,1:point_num) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' '

  call read_input_file ( input_file_name, dim_num, point_num, quasi )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The input data file has been read.'

  if ( minval ( quasi ) < 0.0E+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CIRCLE_TEST_MAIN - Fatal error!'
    write ( *, '(a)' ) '  At least one coordinate of a point is less than 0!'
    stop
  else if ( 1.0E+00 < maxval ( quasi ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CIRCLE_TEST_MAIN - Fatal error!'
    write ( *, '(a)' ) '  At least one coordinate of a point is greater than 1!'
    stop
  end if

  walls = .false.
  call circle_test ( dim_num, point_num, quasi, walls, circle_volume )

  walls = .true.
  call circle_test ( dim_num, point_num, quasi, walls, circle_volume )
!
!  Free  memory.
!
  deallocate ( quasi )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CIRCLE_TEST_MAIN:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine circle_test ( dim_num, point_num, quasi, walls, circle_volume )

!*****************************************************************************80
!
!! CIRCLE_TEST performs the circle test on a set of points.
!
!  Discussion:
!
!    This program computes a measure of even spacing in a set of
!    N-dimensional points.  We will discuss the program as though
!    the space is 2-dimensional; the program may be used for general
!    N-dimensional data.
!
!    The points are assumed to lie in the unit square.  
!
!    The program makes a circle-packing measurement on the points
!    by assuming that, at each point, a circle is centered; all
!    the circles start out with zero radius, and then expand
!    together at the same rate.  A circle stops expanding as soon
!    as it touches any other circle.
!
!    The amount of area covered by the circles is compared to the
!    area of the unit square.  This measurement has a certain amount
!    of boundary effect: some circles will naturally extend outside
!    the unit hypercube.  If this is a concern, is possible to restrict 
!    the circles to remain inside the unit hypercube.  In any case,
!    this problem generally goes away as the number of points increases.
!
!    Since we are interested in the coverage of the unit hypercube,
!    it is probably best if the circles are restricted.  This way,
!    computing the area of the circles gives a measure of the even
!    coverage of the region, relative to the presumably best possible
!    covering, by the same number of circles, but of equal radius.
!
!    In the limit, the maximum relative packing density of a 2D 
!    region with equal-sized circles is 0.9069.  In 3D, a density
!    of at least 0.74 can be achieved, and it is known that no
!    greater than 0.7796 is possible.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 April 2003
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
!    Input, real ( kind = 8 ) QUASI(DIM_NUM,POINT_NUM), the points.
!
!    Input, logical WALLS, is TRUE if the circles are restricted
!    to lie within the unit hypercube.
!
!    Output, real ( kind = 8 ) CIRCLE_VOLUME, the amount of volume taken up
!    by the nonintersecting circles of maximum radius around each
!    point.  Ignoring boundary effects, the "ideal" value would be
!    1 (achievable only in 1 dimension), and the maximum value
!    possible is the sphere packing density in the given spatial
!    dimension.  If boundary effects can be ignored, the value of
!    CIRCLE_VOLUME reports how closely the given set of points
!    behaves like a set of close-packed spheres.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  real ( kind = 8 ) circle_volume
  integer ( kind = 4 ) i
  real ( kind = 8 ), dimension ( dim_num, point_num ) :: quasi
  real ( kind = 8 ), dimension ( point_num ) :: radius
  real ( kind = 8 ) radius_ave
  real ( kind = 8 ) radius_max
  real ( kind = 8 ) radius_min
  logical, parameter :: verbose = .true.
  real ( kind = 8 ) volume
  logical walls

  call radius_maximus ( dim_num, point_num, quasi, walls, radius )

  circle_volume = 0.0E+00
  do i = 1, point_num
    call sphere_volume_nd ( dim_num, radius(i), volume )
    circle_volume = circle_volume + volume
  end do

  if ( verbose ) then

    radius_ave = sum ( radius(1:point_num) ) / real ( point_num, kind = 8 )
    radius_min = minval ( radius(1:point_num) )
    radius_max = maxval ( radius(1:point_num) )

    write ( *, '(a)'      ) ' '
    write ( *, '(a,i6)'   ) '  Number of dimensions is ', dim_num
    write ( *, '(a,i6)'   ) '  Number of points is ', point_num
    if ( walls ) then
      write ( *, '(a)' ) '  Circles are required to stay in the unit square.'
    else
      write ( *, '(a)' ) '  Circles are NOT required to stay in the unit square.'
    end if
    write ( *, '(a)'      ) ' '
    write ( *, '(a,f7.4)' ) '  Average radius = ', radius_ave
    write ( *, '(a,f7.4)' ) '  Minimum radius = ', radius_min
    write ( *, '(a,f7.4)' ) '  Maximum radius = ', radius_max
    write ( *, '(a,f7.4)' ) '  Circle volume =  ', circle_volume
  end if

  return
end
function ch_is_digit ( c )

!*****************************************************************************80
!
!! CH_IS_DIGIT returns .TRUE. if a character is a decimal digit.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C, the character to be analyzed.
!
!    Output, logical CH_IS_DIGIT, .TRUE. if C is a digit, .FALSE. otherwise.
!
  implicit none

  character c
  logical ch_is_digit

  if ( lge ( c, '0' ) .and. lle ( c, '9' ) ) then
    ch_is_digit = .true.
  else
    ch_is_digit = .false.
  end if

  return
end
subroutine ch_to_digit ( c, digit )

!*****************************************************************************80
!
!! CH_TO_DIGIT returns the integer ( kind = 4 ) value of a base 10 digit.
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
subroutine digit_inc ( c )

!*****************************************************************************80
!
!! DIGIT_INC increments a decimal digit.
!
!  Example:
!
!    Input  Output
!    -----  ------
!    '0'    '1'
!    '1'    '2'
!    ...
!    '8'    '9'
!    '9'    '0'
!    'A'    'A'
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
!    Input/output, character C, a digit to be incremented.
!
  implicit none

  character c
  integer ( kind = 4 ) digit

  call ch_to_digit ( c, digit )

  if ( digit == -1 ) then
    return
  end if

  digit = digit + 1

  if ( digit == 10 ) then
    digit = 0
  end if

  call digit_to_ch ( digit, c )

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
subroutine file_column_count ( file_name, ncolumn )

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
!    Input, character ( len = * ) FILE_NAME, the name of the file.
!
!    Output, integer ( kind = 4 ) NCOLUMN, the number of columns assumed to be in the file.
!
  implicit none

  character ( len = * ) file_name
  logical got_one
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  character ( len = 256 ) line
  integer ( kind = 4 ) ncolumn
!
!  Open the file.
!
  call get_unit ( iunit )

  open ( unit = iunit, file = file_name, status = 'old', form = 'formatted', &
    access = 'sequential', iostat = ios )

  if ( ios /= 0 ) then
    ncolumn = - 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_COLUMN_COUNT - Fatal error!'
    write ( *, '(a)' ) '  Could not open the file:'
    write ( *, '(4x,a)' ) '"' // trim ( file_name ) // '".'
    stop
  end if
!
!  Read one line, but skip blank lines and comment lines.
!
  got_one = .false.

  do

    read ( iunit, '(a)', iostat = ios ) line

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

    rewind ( iunit )

    do

      read ( iunit, '(a)', iostat = ios ) line

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

  close ( unit = iunit )

  if ( .not. got_one ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_COLUMN_COUNT - Warning!'
    write ( *, '(a)' ) '  The file does not seem to contain any data.'
    ncolumn = 0
    stop
  end if

  call word_count ( line, ncolumn )

  return
end
subroutine file_line_count ( file_name, nline )

!*****************************************************************************80
!
!! FILE_LINE_COUNT counts the number of lines in a file.
!
!  Discussion:
!
!    The file is assumed to be a simple text file.
!
!    Blank lines and comment lines, which begin with '#', are not counted.
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
!    Input, character ( len = * ) FILE_NAME, the name of the file.
!
!    Output, integer ( kind = 4 ) NLINE, the number of lines found in the file.
!
  implicit none

  character ( len = * ) file_name
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  character ( len = 256 ) line
  integer ( kind = 4 ) nline
  logical, parameter :: verbose = .false.

  nline = 0
!
!  Open the file.
!
  call get_unit ( iunit )

  open ( unit = iunit, file = file_name, status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    nline = - 1
    if ( verbose ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FILE_LINE_COUNT - Fatal error!'
      write ( *, '(a)' ) '  Could not open the file:'
      write ( *, '(4x,a)' ) '"' // trim ( file_name ) // '".'
    end if
    return
  end if
!
!  Count the lines.
!
  do

    read ( iunit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      exit
    end if

    if ( len_trim ( line ) == 0 ) then
      cycle
    end if

    if ( line(1:1) == '#' ) then
      cycle
    end if

    nline = nline + 1

  end do

  close ( unit = iunit )

  return
end
subroutine file_name_inc ( file_name )

!*****************************************************************************80
!
!! FILE_NAME_INC generates the next filename in a series.
!
!  Discussion:
!
!    It is assumed that the digits in the name, whether scattered or
!    connected, represent a number that is to be increased by 1 on
!    each call.  If this number is all 9's on input, the output number
!    is all 0's.  Non-numeric letters of the name are unaffected, and
!    if the name contains no digits, then nothing is done.
!
!  Example:
!
!      Input          Output
!      -----          ------
!      a7to11.txt     a7to12.txt
!      a7to99.txt     a8to00.txt
!      a9to99.txt     a0to00.txt
!      cat.txt        cat.txt
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) FILE_NAME.
!    On input, a character string to be incremented.
!    On output, the incremented string.
!
  implicit none

  character c
  logical ch_is_digit
  character ( len = * ) file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) lens

  lens = len_trim ( file_name )

  do i = lens, 1, -1

    c = file_name(i:i)

    if ( ch_is_digit ( c ) ) then

      call digit_inc ( c )

      file_name(i:i) = c

      if ( c /= '0' ) then
        return
      end if

    end if

  end do

  return
end
subroutine get_input_file_name ( input_file_name )

!*****************************************************************************80
!
!! GET_INPUT_FILE_NAME gets the input file name.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) INPUT_FILE_NAME, the name of the input file.
!
  implicit none

  integer ( kind = 4 ) arg_num
  integer ( kind = 4 ) iarg
  integer ( kind = 4 ) iargc
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ilen
  character ( len = * ) input_file_name
  integer ( kind = 4 ) ipxfargc
!
!  Get the number of command line arguments.
!
!  Old style:
!
  arg_num = iargc ( )
!
!  New style:
!
! arg_num = ipxfargc ( )
!
!  If at least one command line argument, it's the input file name.
!
  if ( 1 <= arg_num ) then

    iarg = 1
!
!  Old style:
!
    call getarg ( iarg, input_file_name )
!
!  New style:
!
!   call pxfgetarg ( iarg, input_file_name, ilen, ierror )
!
!   if ( ierror /= 0 ) then
!     write ( *, '(a)' ) ' '
!     write ( *, '(a)' ) 'GET_INPUT_FILE_NAME - Fatal error!'
!     write ( *, '(a)' ) '  Could not read command line argument.'
!     stop
!   end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GET_INPUT_FILE_NAME:'
    write ( *, '(a)' ) '  Please enter the name of the file containing'
    write ( *, '(a)' ) '  the points to be used in the tests.'

    read ( *, '(a)' ) input_file_name

  end if

  return
end
subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is a vlue between 1 and 99 which
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
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5 and 6).
!
!    Otherwise, IUNIT is a vlue between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  logical lopen

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
function r8_pi ( )

!*****************************************************************************80
!
!! R8_PI returns the value of pi.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) R8_PI, the value of pi.
!
  implicit none

  real ( kind = 8 ) r8_pi

  r8_pi = 3.14159265358979323846264338327950288419716939937510E+00

  return
end
subroutine radius_maximus ( dim_num, point_num, coord, walls, radius )

!*****************************************************************************80
!
!! RADIUS_MAXIMUS finds the biggest possible nonintersecting sphere.
!
!  Discussion:
!
!    We are given a set of POINT_NUM points in DIM_NUM space.  We imagine that
!    at each point simultaneously, a sphere begins to expand.
!    Each sphere stops expanding as soon as it touches another sphere.
!    The radius of these spheres is to be computed.
!
!    If WALLS is true, then the spheres must not extend outside the
!    "walls" of the unit hypersquare.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the number of spatial dimensions.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, real ( kind = 8 ) COORD(DIM_NUM,POINT_NUM), the point coordinates.
!    If WALLS is TRUE, these values must be between 0 and 1.
!
!    Input, logical WALLS, is TRUE if the spheres must not extend
!    outside the unit hypercube.  If WALLS is FALSE, then this
!    restriction is not imposed.
!
!    Output, real ( kind = 8 ) RADIUS(POINT_NUM), the radius of the sphere around each point.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  real ( kind = 8 ) coord(dim_num,point_num)
  real ( kind = 8 ) distance(point_num)
  real ( kind = 8 ) distance_j
  real ( kind = 8 ) distance_min
  integer ( kind = 4 ), parameter :: FIXED = 0
  integer ( kind = 4 ), parameter :: FREE = 1
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) next
  real ( kind = 8 ) radius(point_num)
  real ( kind = 8 ) radius_i
  real ( kind = 8 ) radius_min
  integer ( kind = 4 ) status(point_num)
  logical walls

  if ( walls ) then
         if ( any (           coord(1:dim_num,1:point_num) < 0.0E+00 ) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RADIUS_MAXIMUS - Fatal error!'
      write ( *, '(a)' ) '  Some coordinate is less than 0.'
      stop
    else if ( any ( 1.0E+00 < coord(1:dim_num,1:point_num)           ) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RADIUS_MAXIMUS - Fatal error!'
      write ( *, '(a)' ) '  Some coordinate is greater than 1.'
      stop
    end if
  end if
!
!  Initially, all points are "free".
!
  radius(1:point_num) = 0.0E+00
  status(1:point_num) = FREE

  do
!
!  If all points are fixed, we're done.
!
    if ( all ( status(1:point_num) == FIXED ) ) then
      exit
    end if
!
!  Look at all the free points.
!  Imagine an expanding sphere at each free point, and determine
!  which such sphere will first have to stop expanding.
!
    next = 0
    radius_min = huge ( radius_min )

    do i = 1, point_num

      if ( status(i) == FREE ) then

        if ( walls ) then
          radius_i = min ( &
            minval (           coord(1:dim_num,i) ), &
            minval ( 1.0E+00 - coord(1:dim_num,i) ) )
        else
          radius_i = huge ( radius_i )
        end if

        do j = 1, point_num

          if ( j /= i ) then

            distance_j = sqrt ( sum ( &
              ( coord(1:dim_num,i) - coord(1:dim_num,j) )**2 &
            ) )

            if ( status(j) == FREE ) then
              radius_i = min ( radius_i, distance_j / 2.0E+00 )
            else
              radius_i = min ( radius_i, distance_j - radius(j) )
            end if

          end if

        end do

        if ( radius_i < radius_min ) then
          next = i
          radius_min = radius_i
        end if

      end if

    end do

    if ( next == 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RADIUS_MAXIMUS - Fatal error!'
      write ( *, '(a)' ) '  There were points left to handle, but could'
      write ( *, '(a)' ) '  not choose the "next" one to work on.'
      stop
    end if

    i = next
    radius(i) = radius_min
    status(i) = FIXED

  end do

  return
end
subroutine read_input_file ( input_file_name, ndim, npoint, quasi )

!*****************************************************************************80
!
!! READ_INPUT_FILE reads the data from the input file.
!
!  Discusion:
!
!    Blank lines, and lines beginning with a '#' character, are ignored.
!
!    Any other line is assumed to contain the coordinates of a point
!    in NDIM dimensional space.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILE_NAME, the name of the input file.
!
!    Input, integer ( kind = 4 ) NDIM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) NPOINT, the number of points.
!
!    Output, real ( kind = 8 ) QUASI(NDIM,NPOINT), the points.
!
  implicit none

  integer ( kind = 4 ) ndim
  integer ( kind = 4 ) npoint

  character ( len = * ) :: input_file_name
  integer ( kind = 4 ) input_file_unit
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) j
  real ( kind = 8 ) quasi(ndim,npoint)
  character ( len = 256 ) string

  call get_unit ( input_file_unit )

  open ( unit = input_file_unit, file = input_file_name, status = 'old', &
    iostat = ios )

  do j = 1, npoint

    do

      read ( input_file_unit, '(a)', iostat = ios ) string

      if ( ios /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'READ_INPUT_FILE - Fatal error!'
        write ( *, '(a)' ) '  End of file while reading data.'
        stop
      end if

      if ( string(1:1) /= '#' .and. len_trim ( string ) /= 0 ) then
        exit
      end if

      write ( *, '(a)' ) trim ( string )

    end do

    read ( string, * ) quasi(1:ndim,j)

  end do

  close ( unit = input_file_unit )

  return
end
subroutine sphere_unit_volume_nd ( n, volume )

!*****************************************************************************80
!
!! SPHERE_UNIT_VOLUME_ND computes the volume of a unit sphere in ND.
!
!  Discussion:
!
!    A unit sphere in ND satisfies the equation:
!
!      sum ( ( X(1:N) - XC(1:N) )**2 ) = 1
!
!    where XC is the center.
!
!    N  Volume
!
!    2             PI
!    3  (4/3)    * PI
!    4  (1/2)    * PI**2
!    5  (8/15)   * PI**2
!    6  (1/6)    * PI**3
!    7  (16/105) * PI**3
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the space.
!
!    Output, real ( kind = 8 ) VOLUME, the volume of the sphere.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 8 ) r8_pi
  real ( kind = 8 ) volume

  if ( mod ( n, 2 ) == 0 ) then
    m = n / 2
    volume = ( r8_pi ( ) )**m
    do i = 1, m
      volume = volume / real ( i, kind = 8 )
    end do
  else
    m = ( n - 1 ) / 2
    volume = ( r8_pi ( ) )**m * 2.0E+00**n
    do i = m+1, 2*m+1
      volume = volume / real ( i, kind = 8 )
    end do
  end if

  return
end
subroutine sphere_volume_nd ( n, r, volume )

!*****************************************************************************80
!
!! SPHERE_VOLUME_ND computes the volume of a sphere in ND.
!
!  Discussion:
!
!    A sphere in ND satisfies the equation:
!
!      sum ( ( X(1:N) - XC(1:N) )**2 ) = R**2
!
!    where R is the radius and XC is the center.
!
!    N  Volume
!
!    2             PI    * R**2
!    3  (4/3)    * PI    * R**3
!    4  (1/2)    * PI**2 * R**4
!    5  (8/15)   * PI**2 * R**5
!    6  (1/6)    * PI**3 * R**6
!    7  (16/105) * PI**3 * R**7
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the space.
!
!    Input, real ( kind = 8 ) R, the radius of the sphere.
!
!    Output, real ( kind = 8 ) VOLUME, the volume of the sphere.
!
  implicit none

  integer ( kind = 4 ) n
  real ( kind = 8 ) r
  real ( kind = 8 ) volume

  call sphere_unit_volume_nd ( n, volume )

  volume = volume * r**n

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
  character ( len = 10 )  time
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
subroutine word_count ( s, nword )

!*****************************************************************************80
!
!! WORD_COUNT counts the number of "words" in a string.
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
