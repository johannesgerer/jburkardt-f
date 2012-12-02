  integer ( kind = 4 )program main

!*****************************************************************************80
!
!! MAIN is the main program for INTEGRAL_TEST.
!
!  Discussion:
!
!    INTEGRAL_TEST computes test integrals with a given pointset.
!
!    The user may specify a weight file to be used with the point dataset.
!    If no weight file is specified, then uniform weights are used.
!
!  Usage:
!
!    integral_test point_file [weight_file]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 August 2003
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: integral_num1 = 20
  integer ( kind = 4 ), parameter :: integral_num2 = 4
  integer ( kind = 4 ), parameter :: integral_num3 = 10

  real ( kind = 8 ) average_rel_error
  character ( len = 255 ) input_file_name
  integer ( kind = 4 ), dimension ( integral_num1 ) :: integral_index1 = (/ &
     1,  2,  4,  6,  8,  9, 10, 11, 14, 15, &
    16, 17, 18, 19, 24, 25, 26, 28, 30, 31 /)
  integer ( kind = 4 ), dimension ( integral_num2 ) :: integral_index2 = (/ &
    15, 16, 17, 19 /)
  real ( kind = 8 ) integral_rel_error1(integral_num1)
  real ( kind = 8 ) integral_rel_error2(integral_num2)
  real ( kind = 8 ) integral_rel_error3(integral_num3)
  integer ( kind = 4 ) ndim
  integer ( kind = 4 ) ndim1
  integer ( kind = 4 ) npoint
  integer ( kind = 4 ) npoint1
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: quasi
  real ( kind = 8 ), allocatable, dimension ( : ) :: weight
  character ( len = 255 ) weight_file_name

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'INTEGRAL_TEST'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test a file of sample points for use in'
  write ( *, '(a)' ) '  multidimensional quadrature.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The points should lie in the unit hypercube.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The points may have associated weights.'
  write ( *, '(a)' ) '  If no weights are supplied, uniform weights'
  write ( *, '(a)' ) '  are used.'

  call get_input_file_name ( input_file_name, weight_file_name )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The point data will be read from "' &
    // trim ( input_file_name ) // '".'
  if ( 0 < len ( weight_file_name ) ) then
    write ( *, '(a)' ) '  The weight data will be read from "' &
      // trim ( weight_file_name ) // '".'
  else
    write ( *, '(a)' ) &
      '  No weight file was specified, so uniform weights are used.'
  end if

  call file_column_count ( input_file_name, ndim )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The spatial dimension is ', ndim

  if ( ndim <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'INTEGRAL_TEST - Fatal error!'
    write ( *, '(a)' ) '  The spatial dimension must be at least 1!'
    stop
  end if

  call file_line_count ( input_file_name, npoint )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of points is ', npoint

  if ( npoint <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'INTEGRAL_TEST - Fatal error!'
    write ( *, '(a)' ) '  The number of points must be at least 1!'
    stop
  end if

  allocate ( quasi(1:ndim,1:npoint) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' '

  call read_input_file ( input_file_name, ndim, npoint, quasi )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The input data file has been read.'

  if ( minval ( quasi ) < 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'INTEGRAL_TEST - Fatal error!'
    write ( *, '(a)' ) '  At least one coordinate of a point is less than 0!'
    stop
  else if ( 1.0D+00 < maxval ( quasi ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'INTEGRAL_TEST - Fatal error!'
    write ( *, '(a)' ) '  At least one coordinate of a point is greater than 1!'
    stop
  end if
!
!  If a weight file was specified, try to read that now.
!
  if ( 0 < len_trim ( weight_file_name ) ) then

    call file_column_count ( weight_file_name, ndim1 )

    if ( ndim1 /= 1 ) then

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'INTEGRAL_TEST - Fatal error!'
      write ( *, '(a)' ) '  The weight data was not given in column format.'
      write ( *, '(a)' ) '  The weights will not be used.'
      stop

    else

      call file_line_count ( weight_file_name, npoint1 )

      if ( npoint1 /= npoint ) then

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'INTEGRAL_TEST - Fatal error!'
        write ( *, '(a)' ) '  The number of weights does not match'
        write ( *, '(a)' ) '  the number of data points.'
        write ( *, '(a)' ) '  The weights will not be used.'
        stop

      else

        allocate ( weight(1:npoint) )

        call read_input_file ( weight_file_name, 1, npoint, weight )

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  The file of weights was successfully read.'

      end if

    end if

  else

    allocate ( weight(1:npoint) )

    weight(1:npoint) = 1.0D+00 / real ( npoint, kind = 8 )

  end if

  if ( .false. ) then
    call list_titles ( integral_num1, integral_index1 )
  end if

  call integral_test_01 ( ndim, npoint, quasi, integral_num1, &
    integral_index1, weight, integral_rel_error1 )

  if ( .false. ) then
    call list_titles ( integral_num2, integral_index2 )
  end if

  call integral_test_02 ( ndim, npoint, quasi, integral_num2, &
    integral_index2, weight, integral_rel_error2 )

  call integral_test_03 ( ndim, npoint, quasi, integral_num3, &
    weight, integral_rel_error3 )

  deallocate ( quasi )
  deallocate ( weight )
!
!  Averate the relative errors.
!
  average_rel_error = ( &
      sum ( integral_rel_error1(1:integral_num1) ) &
    + sum ( integral_rel_error2(1:integral_num2) ) ) &
    / real ( integral_num1 + integral_num2, kind = 8 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'INTEGRAL_TEST:'
  write ( *, '(a,g14.6)' ) '  Averaged relative error: ', average_rel_error
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'INTEGRAL_TEST:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
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
  character ( len = 255 ) line
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

  call s_word_count ( line, ncolumn )

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
!    If the file could not be opened, then NLINE is returned as -1.
!
  implicit none

  character ( len = * ) file_name
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  character ( len = 255 ) line
  integer ( kind = 4 ) nline
  logical, parameter :: verbose = .false.

  nline = 0
!
!  Open the file.
!
  call get_unit ( iunit )

  open ( unit = iunit, file = file_name, status = 'old', form = 'formatted', &
    access = 'sequential', iostat = ios )

  if ( ios /= 0 ) then

    nline = - 1

    if ( verbose ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FILE_LINE_COUNT - Fatal error!'
      write ( *, '(a)' ) '  Could not open the file:'
      write ( *, '(4x,a)' ) '"' // trim ( file_name ) // '".'
    end if

    stop

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
subroutine get_input_file_name ( input_file_name, weight_file_name )

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
!    20 May 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) INPUT_FILE_NAME, the name of the input file.
!
!    Output, character ( len = * ) WEIGHT_FILE_NAME, the name of the 
!    optional input weight file, if specified on the command line.
!
  implicit none

  integer ( kind = 4 ) arg_num
  integer ( kind = 4 ) iarg
  integer ( kind = 4 ) iargc
  character ( len = * ) input_file_name
  character ( len = * ) weight_file_name
!
!  Get the number of command line arguments.
!
  arg_num = iargc ( )
!
!  If at least one command line argument, it's the input file name.
!
  if ( 1 <= arg_num ) then

    iarg = 1
    call getarg ( iarg, input_file_name )

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GET_INPUT_FILE_NAME:'
    write ( *, '(a)' ) '  Please enter the name of the file containing'
    write ( *, '(a)' ) '  the points to be used in the tests.'

    read ( *, '(a)' ) input_file_name

  end if
!
!  If at least two command line arguments, the second one is the weight 
!  file name.
!
  if ( 2 <= arg_num ) then

    iarg = 2
    call getarg ( iarg, weight_file_name )

  else

    weight_file_name = ' '

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
subroutine integral_test_01 ( ndim, npoint, quasi, integral_num, &
  integral_index, weight, integral_rel_error )

!*****************************************************************************80
!
!! INTEGRAL_TEST_01 does the general integral tests.
!
!  Discussion:
!
!    This program tests the accuracy of numerical integration
!    using sets of points described by files of XY coordinates.
!
!    The data is assumed to be in files with the names 'set_01.txt',
!    'set_02.txt' and so on.  As long as the files are named
!    consecutively, the code will read and process them all.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 June 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NDIM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) NPOINT, the number of points.
!
!    Input, real ( kind = 8 ) QUASI(NDIM,NPOINT), the sample points.
!
!    Input, integer ( kind = 4 ) INTEGRAL_NUM, the number of integrals to approximate.
!
!    Input, integer ( kind = 4 ) INTEGRAL_INDEX(INTEGRAL_NUM), the indices of the
!    integrals to approximate.
!
!    Input, real ( kind = 8 ) WEIGHT(NPOINT), the weights.
!
!    Output, real ( kind = 8 ) INTEGRAL_REL_ERROR(INTEGRAL_NUM), the 
!    relative error in each approximated integral.
!
  implicit none

  integer ( kind = 4 ) integral_num
  integer ( kind = 4 ) ndim
  integer ( kind = 4 ) npoint

  integer ( kind = 4 ) i
  real ( kind = 8 ) integral_error
  real ( kind = 8 ) integral_estimate
  real ( kind = 8 ) integral_exact
  integer ( kind = 4 ), dimension ( integral_num ) :: integral_index
  real ( kind = 8 ), dimension ( integral_num ) :: integral_rel_error
  integer ( kind = 4 ) iprob
  integer ( kind = 4 ) num
  real ( kind = 8 ) p00_f
  real ( kind = 8 ), dimension ( ndim, npoint ) :: quasi
  real ( kind = 8 ) total
  real ( kind = 8 ) volume
  real ( kind = 8 ) weight(npoint)
  real ( kind = 8 ), dimension ( ndim, npoint ) :: x
  real ( kind = 8 ) x0(ndim)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'INTEGRAL_TEST_01'
  write ( *, '(a)' ) '  Using a set of points in [0,1]^N,'
  write ( *, '(a)' ) '  approximate integrals from the TESTNINT set.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  Problem    Region        Correct       Estimated      Absolute      Relative'
  write ( *, '(a)' ) &
    '  Index      Volume        Integral      Integral       Error         Error (%)'
  write ( *, '(a)' ) ' '

  do num = 1, integral_num

    iprob = integral_index ( num )
!
!  Set the real vectors to default values.
!
    call p00_r8vec ( iprob, 'DEFAULT', '*', ndim, x0 )
!
!  Get the volume of the integration domain.
!
    call p00_volume ( iprob, ndim, volume )
!
!  Get the exact value of the integral.
!
    call p00_exact ( iprob, ndim, integral_exact )
!
!  Map the points in [0,1] into the integration domain [A,B].
!
    call p00_remap01 ( iprob, ndim, npoint, quasi, x ) 
!
!  Estimate the integral.
!
    total = 0.0D+00
    do i = 1, npoint
      total = total + weight(i) * p00_f ( iprob, ndim, x(1:ndim,i) )
    end do

    integral_estimate = volume * total

    integral_error = abs ( integral_estimate - integral_exact )

    if ( integral_exact /= 0.0D+00 ) then
      integral_rel_error(num) = 100.0D+00 * abs ( integral_error ) &
        / abs ( integral_exact )
    else
      integral_rel_error(num) = 100.0D+00 * abs ( integral_error )
    end if

    write ( *, '(i8,2x,3g14.6,2g14.6)' ) &
      iprob, volume, integral_exact, integral_estimate, integral_error, &
      integral_rel_error(num)

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'INTEGRAL_TEST_01:'
  write ( *, '(a)' ) '  Normal conclusion.'

  return
end
subroutine integral_test_02 ( ndim, npoint, quasi, integral_num, &
  integral_index, weight, integral_rel_error )

!*****************************************************************************80
!
!! INTEGRAL_TEST_02 does some integrals involving basepoints.
!
!  Discussion:
!
!    Approximate certain test integrals, which depend on a base point X0.
!    The base point will be varied randomly, and the integral evaluated
!    100 times.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NDIM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) NPOINT, the number of points.
!
!    Input, real ( kind = 8 ) QUASI(NDIM,NPOINT), the sample points.
!
!    Input, integer ( kind = 4 ) INTEGRAL_NUM, the number of integrals to approximate.
!
!    Input, integer ( kind = 4 ) INTEGRAL_INDEX(INTEGRAL_NUM), the indices of the
!    integrals to approximate.
!
!    Input, real ( kind = 8 ) WEIGHT(NPOINT), the weights.
!
!    Output, real ( kind = 8 ) INTEGRAL_REL_ERROR(INTEGRAL_NUM), the 
!    relative error in each approximated integral.
!
  implicit none

  integer ( kind = 4 ) integral_num
  integer ( kind = 4 ) ndim
  integer ( kind = 4 ) npoint

  real ( kind = 8 ) era_ave
  real ( kind = 8 ) era_max
  real ( kind = 8 ) era_min
  real ( kind = 8 ) era_now
  real ( kind = 8 ) err_ave
  real ( kind = 8 ) err_max
  real ( kind = 8 ) err_min
  real ( kind = 8 ) err_now
  real ( kind = 8 ) ext_ave
  real ( kind = 8 ) ext_max
  real ( kind = 8 ) ext_min
  real ( kind = 8 ) ext_now
  integer ( kind = 4 ) i
  real ( kind = 8 ) int_ave
  real ( kind = 8 ) int_max
  real ( kind = 8 ) int_min
  real ( kind = 8 ) int_now
  integer ( kind = 4 ) integral_index(integral_num)
  real ( kind = 8 ) integral_rel_error(integral_num)
  integer ( kind = 4 ) iprob
  integer ( kind = 4 ) num
  real ( kind = 8 ) p00_f
  real ( kind = 8 ), dimension ( ndim, npoint ) :: quasi
  integer ( kind = 4 ) rep
  integer ( kind = 4 ), parameter :: repeat = 100
  integer ( kind = 4 ) seed
  real ( kind = 8 ) total
  real ( kind = 8 ) volume
  real ( kind = 8 ) weight(npoint)
  real ( kind = 8 ), dimension ( ndim, npoint ) :: x
  real ( kind = 8 ) x0(ndim)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'INTEGRAL_TEST_02'
  write ( *, '(a)' ) '  Using a set of points in [0,1]^N,'
  write ( *, '(a)' ) '  approximate integrals from the TESTNINT set.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The integrands depend on a base point.'
  write ( *, '(a)' ) '  We randomly vary the base point repeatedly,'
  write ( *, '(a)' ) '  and look at the average behavior.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Because the base point moves, the value of the'
  write ( *, '(a)' ) '  exact integral changes on each repetition!'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of repetitions = ', repeat
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '             Average       Maximum       Average       Maximum'
  write ( *, '(a)' ) '  Integral   Absolute      Absolute      Relative      Relative'
  write ( *, '(a)' ) '             Error         Error         Error %       Error %'
  write ( *, '(a)' ) ' '

  do num = 1, integral_num

    seed = 123456789
    call random_initialize ( seed )

    iprob = integral_index ( num )
!
!  Set the real vectors to default values.
!
    call p00_r8vec ( iprob, 'DEFAULT', '*', ndim, x0 )
!
!  Get the volume.
!
    call p00_volume ( iprob, ndim, volume )
!
!  Map the points in [0,1] into [A,B].
!
    call p00_remap01 ( iprob, ndim, npoint, quasi, x ) 

    era_ave = 0.0D+00
    err_ave = 0.0D+00
    ext_ave = 0.0D+00
    int_ave = 0.0D+00

    era_max = - huge ( era_max )
    err_max = - huge ( err_max )
    ext_max = - huge ( ext_max )
    int_max = - huge ( int_max )

    era_min = huge ( era_min )
    err_min = huge ( err_min )
    ext_min = huge ( ext_min )
    int_min = huge ( int_min )

    integral_rel_error(num) = 0.0D+00

    do rep = 1, repeat
!
!  Randomize the base point x0.
!
      call p00_r8vec ( iprob, 'RANDOMIZE', 'X0', ndim, x0 )
!
!  Get the exact value of the integral (it may depend on X0).
!
      call p00_exact ( iprob, ndim, ext_now )

      ext_ave = ext_ave + ext_now
      ext_max = max ( ext_max, ext_now )
      ext_min = min ( ext_min, ext_now )
!
!  Estimate the integral.
!
      total = 0.0D+00

      do i = 1, npoint
        total = total + weight(i) * p00_f ( iprob, ndim, x(1:ndim,i) )
      end do

      int_now = volume * total

      int_ave = int_ave + int_now
      int_max = max ( int_max, int_now )
      int_min = min ( int_min, int_now )

      era_now = abs ( int_now - ext_now )

      era_ave = era_ave + era_now
      era_max = max ( era_max, era_now )
      era_min = min ( era_min, era_now )

      if ( ext_now == 0.0D+00 ) then
        err_now = abs ( 100.0D+00 * ( int_now - ext_now ) )
      else
        err_now = abs ( 100.0D+00 * ( int_now - ext_now ) ) / abs ( ext_now )
      end if

      integral_rel_error(num) = integral_rel_error(num) + err_now

      err_ave = err_ave + err_now
      err_max = max ( err_max, err_now )
      err_min = min ( err_min, err_now )

    end do

    integral_rel_error(num) = integral_rel_error(num) &
      / real ( repeat, kind = 8 )

    era_ave = era_ave / real ( repeat, kind = 8 )
    err_ave = err_ave / real ( repeat, kind = 8 )
    ext_ave = ext_ave / real ( repeat, kind = 8 )
    int_ave = int_ave / real ( repeat, kind = 8 )

    write ( *, '(i8,2x,4f12.6)' ) iprob, era_ave, era_max, err_ave, err_max

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'INTEGRAL_TEST_02:'
  write ( *, '(a)' ) '  Normal conclusion.'

  return
end
subroutine integral_test_03 ( ndim, npoint, quasi, integral_num, &
  weight, integral_rel_error )

!*****************************************************************************80
!
!! INTEGRAL_TEST_03 does integral tests for one simple integrand.
!
!  Discussion:
!
!    This routine concentrates on a single test integral, of the
!    form 
!
!      F(X(1:NDIM),ORDER) = product ( 1 <= I <= ORDER ) X(MOD(I-1,NDIM)+1)
!
!    We allow ORDER to go from 0 to INTEGRAL_NUM-1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 August 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NDIM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) NPOINT, the number of points.
!
!    Input, real ( kind = 8 ) QUASI(NDIM,NPOINT), the sample points.
!
!    Input, integer ( kind = 4 ) INTEGRAL_NUM, the number of integrals to approximate.
!
!    Input, real ( kind = 8 ) WEIGHT(NPOINT), the weights.
!
!    Output, real ( kind = 8 ) INTEGRAL_REL_ERROR(INTEGRAL_NUM), the 
!    relative error in each approximated integral.
!
  implicit none

  integer ( kind = 4 ) integral_num
  integer ( kind = 4 ) ndim
  integer ( kind = 4 ) npoint

  integer ( kind = 4 ) i
  real ( kind = 8 ) integral_error
  real ( kind = 8 ) integral_estimate
  real ( kind = 8 ) integral_exact
  real ( kind = 8 ), dimension ( integral_num ) :: integral_rel_error
  integer ( kind = 4 ) iprob
  integer ( kind = 4 ) num
  integer ( kind = 4 ) order
  real ( kind = 8 ) p00_f
  real ( kind = 8 ), dimension ( ndim, npoint ) :: quasi
  real ( kind = 8 ) total
  real ( kind = 8 ) volume
  real ( kind = 8 ) weight(npoint)
  real ( kind = 8 ), dimension ( ndim, npoint ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'INTEGRAL_TEST_03'
  write ( *, '(a)' ) '  Using a set of points in [0,1]^N,'
  write ( *, '(a)' ) '  approximate integrals from the TESTNINT set.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Integrand J has the form:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    product ( 1 <= I <= J ) X(mod(I,N)+1)'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  Problem    Region        Correct       Estimated      Absolute      Relative'
  write ( *, '(a)' ) &
    '  Index      Volume        Integral      Integral       Error         Error (%)'
  write ( *, '(a)' ) ' '

  iprob = 33

  do num = 1, integral_num

    order = num - 1

    volume = 1.0D+00

    call p00_i4 ( iprob, 'SET', 'ORDER', order )
!
!  Get the volume of the integration domain.
!
    call p00_volume ( iprob, ndim, volume )
!
!  Get the exact value of the integral.
!
    call p00_exact ( iprob, ndim, integral_exact )
!
!  Map the points in [0,1] into the integration domain [A,B].
!
    call p00_remap01 ( iprob, ndim, npoint, quasi, x ) 
!
!  Estimate the integral.
!
    total = 0.0D+00
    do i = 1, npoint
      total = total + weight(i) * p00_f ( iprob, ndim, x(1:ndim,i) )
    end do

    integral_estimate = volume * total

    integral_error = abs ( integral_estimate - integral_exact )

    if ( integral_exact /= 0.0D+00 ) then
      integral_rel_error(num) = 100.0D+00 * abs ( integral_error ) &
        / abs ( integral_exact )
    else
      integral_rel_error(num) = 100.0D+00 * abs ( integral_error )
    end if

    write ( *, '(i8,2x,3g14.6,2g14.6)' ) &
      num, volume, integral_exact, integral_estimate, integral_error, &
      integral_rel_error(num)

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'INTEGRAL_TEST_03:'
  write ( *, '(a)' ) '  Normal conclusion.'

  return
end
subroutine list_titles ( integral_num, integral_index )

!*****************************************************************************80
!
!! LIST_TITLES lists the titles of the integral tests.
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
!    Input, integer ( kind = 4 ) INTEGRAL_NUM, the number of integrals to approximate.
!
!    Input, integer ( kind = 4 ) INTEGRAL_INDEX(INTEGRAL_NUM), the indices of the
!    integrals to approximate.
!
  implicit none

  integer ( kind = 4 ) integral_num

  integer ( kind = 4 ) i
  integer ( kind = 4 ) index
  integer ( kind = 4 ) integral_index(integral_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LIST_TITLES'
  write ( *, '(a)' ) '  Descriptions of selected test integrals:'

  do i = 1, integral_num

    index = integral_index(i)   

    call p00_title ( index )

  end do

  return
end
subroutine random_initialize ( seed )

!*****************************************************************************80
!
!! RANDOM_INITIALIZE initializes the FORTRAN 90 random number seed.
!
!  Discussion:
!
!    If you don't initialize the random number generator, its behavior
!    is not specified.  If you initialize it simply by:
!
!      call random_seed ( )
!
!    its behavior is not specified.  On the DEC ALPHA, if that's all you
!    do, the same random number sequence is returned.  In order to actually
!    try to scramble up the random number generator a bit, this routine
!    goes through the tedious process of getting the size of the random
!    number seed, making up values based on the current time, and setting
!    the random number seed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED.
!    If SEED is zero on input, then you're asking this routine to come up
!    with a seed value, which is returned as output.
!    If SEED is nonzero on input, then you're asking this routine to
!    use the input value of SEED to initialize the random number generator,
!    and SEED is not changed on output.
!
  implicit none

  integer ( kind = 4 ) count
  integer ( kind = 4 ) count_max
  integer ( kind = 4 ) count_rate
  logical, parameter :: debug = .false.
  integer ( kind = 4 ) i
  integer ( kind = 4 ) seed
  integer ( kind = 4 ), allocatable :: seed_vector(:)
  integer ( kind = 4 ) seed_size
  real ( kind = 8 ) t
!
!  Initialize the random number seed.
!
  call random_seed ( )
!
!  Determine the size of the random number seed.
!
  call random_seed ( size = seed_size )
!
!  Allocate a seed of the right size.
!
  allocate ( seed_vector(seed_size) )

  if ( seed /= 0 ) then

    if ( debug ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RANDOM_INITIALIZE'
      write ( *, '(a,i20)' ) '  Initialize RANDOM_NUMBER, user SEED = ', seed
    end if

  else

    call system_clock ( count, count_rate, count_max )

    seed = count

    if ( debug ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RANDOM_INITIALIZE'
      write ( *, '(a,i20)' ) '  Initialize RANDOM_NUMBER, arbitrary SEED = ', &
        seed
    end if

  end if
!
!  Now set the seed.
!
  seed_vector(1:seed_size) = seed

  call random_seed ( put = seed_vector(1:seed_size) )
!
!  Free up the seed space.
!
  deallocate ( seed_vector )
!
!  Call the random number routine a bunch of times.
!
  do i = 1, 100
    call random_number ( harvest = t )
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
  integer ( kind = 4 )   ( kind = 4 ) d
  integer ( kind = 4 )   ( kind = 4 ) h
  integer ( kind = 4 )   ( kind = 4 ) m
  integer ( kind = 4 )   ( kind = 4 ) mm
  character ( len = 9  ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 )   ( kind = 4 ) n
  integer ( kind = 4 )   ( kind = 4 ) s
  integer ( kind = 4 )   ( kind = 4 ) values(8)
  integer ( kind = 4 )   ( kind = 4 ) y

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
subroutine timestring ( string )

!*****************************************************************************80
!
!! TIMESTRING writes the current YMDHMS date into a string.
!
!  Example:
!
!    STRING = 'May 31 2001   9:45:54.872 AM'
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) STRING, contains the date information.
!    A character length of 40 should always be sufficient.
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
  character ( len = * ) string
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

  write ( string, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
