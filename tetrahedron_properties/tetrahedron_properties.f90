program main

!*****************************************************************************80
!
!! MAIN is the main program for TETRAHEDRON_PROPERTIES.
!
!  Discussion:
!
!    TETRAHEDRON_PROPERTIES reports properties of a tetrahedron.
!
!  Usage:
!
!    tetrahedron_properties filename
!
!    where "filename" is a file containing the coordinates of the vertices.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 August 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) arg_num
  real ( kind = 8 ) centroid(3)
  real ( kind = 8 ) circum_center(3)
  real ( kind = 8 ) circum_radius
  real ( kind = 8 ) dihedral_angles(6)
  integer ( kind = 4 ) dim_num
  real ( kind = 8 ) edge_length(6)
  real ( kind = 8 ) face_angles(3,4)
  real ( kind = 8 ) face_areas(4)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iarg
  integer ( kind = 4 ) iargc
  real ( kind = 8 ) in_center(3)
  real ( kind = 8 ) in_radius
  character ( len = 255 ) node_filename
  integer ( kind = 4 ) node_num
  real ( kind = 8 )  node_xyz(3,4)
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) quality1
  real ( kind = 8 ) quality2
  real ( kind = 8 ) quality3
  real ( kind = 8 ) quality4
  real ( kind = 8 ) solid_angles(4)
  real ( kind = 8 ) volume

  write ( *, '(a)' ) ' '
  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TETRAHEDRON_PROPERTIES:'
  write ( *, '(a)' ) '  FORTRAN90 version:'
  write ( *, '(a)' ) '  Determine properties of a tetrahedron.'
!
!  Get the number of command line arguments.
!
  arg_num = iargc ( )
!
!  Commandline argument #1 is the file name
!
  if ( 1 <= arg_num ) then

    iarg = 1
    call getarg ( iarg, node_filename )

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TETRAHEDRON_PROPERTIES:'
    write ( *, '(a)' ) '  Please enter the name of the node coordinate file.'

    read ( *, '(a)' ) node_filename

  end if
!
!  Read the node data.
!
  call r8mat_header_read ( node_filename, dim_num, node_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Read the header of "' &
    // trim ( node_filename ) //'".'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Spatial dimension DIM_NUM = ', dim_num
  write ( *, '(a,i8)' ) '  Number of points NODE_NUM = ', node_num

  if ( dim_num /= 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TETRAHEDRON_PROPERTIES - Fatal error!'
    write ( *, '(a)' ) '  Dataset must have spatial dimension 3.'
    stop
  end if

  if ( node_num /= 4 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TETRAHEDRON_PROPERTIES - Fatal error!'
    write ( *, '(a)' ) '  Dataset must have 4 nodes.'
    stop
  end if

  call r8mat_data_read ( node_filename, dim_num, node_num, node_xyz )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Read the data in "' &
    // trim ( node_filename ) //'".'

  call r8mat_transpose_print ( dim_num, node_num, node_xyz, '  Node coordinates:' )
!
!  CENTROID
!
  call tetrahedron_centroid_3d ( node_xyz, centroid )

  write ( *, '(a)' ) ' '
  write ( *, '(a,3g14.6)' ) '  CENTROID: ', centroid(1:3)
!
!  CIRCUMSPHERE
!
  call tetrahedron_circumsphere_3d ( node_xyz, circum_radius, circum_center )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  CIRCUM_RADIUS = ', circum_radius
  write ( *, '(a,3g14.6)' ) '  CIRCUM_CENTER: ', circum_center(1:3)
!
!  DIHEDRAL ANGLES
!
  call tetrahedron_dihedral_angles_3d ( node_xyz, dihedral_angles )

  call r8vec_print ( 6, dihedral_angles, '  DIHEDRAL_ANGLES (radians)' )

  dihedral_angles(1:6) = dihedral_angles(1:6) * 180.0D+00 / pi

  call r8vec_print ( 6, dihedral_angles, '  DIHEDRAL_ANGLES (degrees)' )
!
!  EDGE LENGTHS
!
  call tetrahedron_edge_length_3d ( node_xyz, edge_length )

  call r8vec_print ( 6, edge_length, '  EDGE_LENGTHS' )
!
!  FACE ANGLES
!
  call tetrahedron_face_angles_3d ( node_xyz, face_angles )

  call r8mat_transpose_print ( 3, 4, face_angles, '  FACE_ANGLES (radians)' )

  face_angles(1:3,1:4) = face_angles(1:3,1:4) * 180.0D+00 / pi

  call r8mat_transpose_print ( 3, 4, face_angles, '  FACE_ANGLES (degrees)' )
!
!  FACE AREAS
!
  call tetrahedron_face_areas_3d ( node_xyz, face_areas )

  call r8vec_print ( 4, face_areas, '  FACE_AREAS' )
!
!  INSPHERE
!
  call tetrahedron_insphere_3d ( node_xyz, in_radius, in_center )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  IN_RADIUS = ', in_radius
  write ( *, '(a,3g14.6)' ) '  IN_CENTER: ', in_center(1:3)
!
!  QUALITY1
!
  call tetrahedron_quality1_3d ( node_xyz, quality1 )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  QUALITY1 = ', quality1
!
!  QUALITY2
!
  call tetrahedron_quality2_3d ( node_xyz, quality2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  QUALITY2 = ', quality2
!
!  QUALITY3
!
  call tetrahedron_quality3_3d ( node_xyz, quality3 )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  QUALITY3 = ', quality3
!
!  QUALITY4
!
  call tetrahedron_quality4_3d ( node_xyz, quality4 )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  QUALITY4 = ', quality4
!
!  SOLID ANGLES
!
  call tetrahedron_solid_angles_3d ( node_xyz, solid_angles )

  call r8vec_print ( 4, solid_angles, '  SOLID_ANGLES (steradians)' )
!
!  VOLUME
!
  call tetrahedron_volume_3d ( node_xyz, volume )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  VOLUME = ', volume
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TETRAHEDRON_PROPERTIES:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
function arc_cosine ( c )

!*****************************************************************************80
!
!! ARC_COSINE computes the arc cosine function, with argument truncation.
!
!  Discussion:
!
!    If you call your system ACOS routine with an input argument that is
!    even slightly outside the range [-1.0, 1.0 ], you may get an unpleasant
!    surprise (I did).
!
!    This routine simply truncates arguments outside the range.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) C, the argument.
!
!    Output, real ( kind = 8 ) ARC_COSINE, an angle whose cosine is C.
!
  implicit none

  real ( kind = 8 ) arc_cosine
  real ( kind = 8 ) c
  real ( kind = 8 ) c2

  c2 = c
  c2 = max ( c2, -1.0D+00 )
  c2 = min ( c2, +1.0D+00 )

  arc_cosine = acos ( c2 )

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

  character              c
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

  logical   ch_eqi
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

  character              c
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
subroutine file_column_count ( input_file_name, column_num )

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
!    Input, character ( len = * ) INPUT_FILE_NAME, the name of the file.
!
!    Output, integer ( kind = 4 ) COLUMN_NUM, the number of columns in the file.
!
  implicit none

  integer ( kind = 4 )   column_num
  logical                  got_one
  character ( len = * )    input_file_name
  integer ( kind = 4 )   input_status
  integer ( kind = 4 )   input_unit
  character ( len = 255 )  line
!
!  Open the file.
!
  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_file_name, status = 'old', &
    form = 'formatted', access = 'sequential', iostat = input_status )

  if ( input_status /= 0 ) then
    column_num = -1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_COLUMN_COUNT - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the input file "' &
      // trim ( input_file_name ) // '" on unit ', input_unit
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
subroutine file_row_count ( input_file_name, row_num )

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
!    Input, character ( len = * ) INPUT_FILE_NAME, the name of the input file.
!
!    Output, integer ( kind = 4 ) ROW_NUM, the number of rows found.
!
  implicit none

  integer ( kind = 4 )   bad_num
  integer ( kind = 4 )   comment_num
  integer ( kind = 4 )   ierror
  character ( len = * )    input_file_name
  integer ( kind = 4 )   input_status
  integer ( kind = 4 )   input_unit
  character ( len = 255 )  line
  integer ( kind = 4 )   record_num
  integer ( kind = 4 )   row_num

  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_file_name, status = 'old', &
    iostat = input_status )

  if ( input_status /= 0 ) then
    row_num = -1;
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_ROW_COUNT - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the input file "' // &
      trim ( input_file_name ) // '" on unit ', input_unit
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
  logical              lopen

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
subroutine r8_swap ( x, y )

!*****************************************************************************80
!
!! R8_SWAP swaps two R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X, Y.  On output, the values of X and
!    Y have been interchanged.
!
  implicit none

  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  z = x
  x = y
  y = z

  return
end
function r8mat_det_4d ( a )

!*****************************************************************************80
!
!! R8MAT_DET_4D computes the determinant of a 4 by 4 matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A(4,4), the matrix whose determinant is desired.
!
!    Output, real ( kind = 8 ) R8MAT_DET_4D, the determinant of the matrix.
!
  implicit none

  real ( kind = 8 ) a(4,4)
  real ( kind = 8 ) r8mat_det_4d

  r8mat_det_4d = &
      a(1,1) * ( &
        a(2,2) * ( a(3,3) * a(4,4) - a(3,4) * a(4,3) ) &
      - a(2,3) * ( a(3,2) * a(4,4) - a(3,4) * a(4,2) ) &
      + a(2,4) * ( a(3,2) * a(4,3) - a(3,3) * a(4,2) ) ) &
    - a(1,2) * ( &
        a(2,1) * ( a(3,3) * a(4,4) - a(3,4) * a(4,3) ) &
      - a(2,3) * ( a(3,1) * a(4,4) - a(3,4) * a(4,1) ) &
      + a(2,4) * ( a(3,1) * a(4,3) - a(3,3) * a(4,1) ) ) &
    + a(1,3) * ( &
        a(2,1) * ( a(3,2) * a(4,4) - a(3,4) * a(4,2) ) &
      - a(2,2) * ( a(3,1) * a(4,4) - a(3,4) * a(4,1) ) &
      + a(2,4) * ( a(3,1) * a(4,2) - a(3,2) * a(4,1) ) ) &
    - a(1,4) * ( &
        a(2,1) * ( a(3,2) * a(4,3) - a(3,3) * a(4,2) ) &
      - a(2,2) * ( a(3,1) * a(4,3) - a(3,3) * a(4,1) ) &
      + a(2,3) * ( a(3,1) * a(4,2) - a(3,2) * a(4,1) ) )

  return
end
subroutine r8mat_data_read ( input_file_name, m, n, table )

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
!    Input, character ( len = * ) INPUT_FILE_NAME, the name of the input file.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Output, real ( kind = 8 ) TABLE(M,N), the table data.
!
  implicit none

  integer ( kind = 4 )   m
  integer ( kind = 4 )   n

  integer ( kind = 4 )   ierror
  character ( len = * )    input_file_name
  integer ( kind = 4 )   input_status
  integer ( kind = 4 )   input_unit
  integer ( kind = 4 )   j
  character ( len = 255 )  line
  real ( kind = 8 )   table(m,n)
  real ( kind = 8 )   x(m)

  ierror = 0

  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_file_name, status = 'old', &
    iostat = input_status )

  if ( input_status /= 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_DATA_READ - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the input file "' // &
      trim ( input_file_name ) // '" on unit ', input_unit
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
subroutine r8mat_header_read ( input_file_name, m, n )

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
!    Input, character ( len = * ) INPUT_FILE_NAME, the name of the input file.
!
!    Output, integer ( kind = 4 ) M, spatial dimension.
!
!    Output, integer ( kind = 4 ) N, the number of points.
!
  implicit none

  character ( len = * )  input_file_name
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  call file_column_count ( input_file_name, m )

  if ( m <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_HEADER_READ - Fatal error!'
    write ( *, '(a)' ) '  There was some kind of I/O problem while trying'
    write ( *, '(a)' ) '  to count the number of data columns in'
    write ( *, '(a)' ) '  the file "' // trim ( input_file_name ) // '".'
    stop
  end if

  call file_row_count ( input_file_name, n )

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_HEADER_READ - Fatal error!'
    write ( *, '(a)' ) '  There was some kind of I/O problem while trying'
    write ( *, '(a)' ) '  to count the number of data rows in'
    write ( *, '(a)' ) '  the file "' // trim ( input_file_name ) // '".'
    stop
  end if

  return
end
subroutine r8mat_solve ( n, rhs_num, a, info )

!*****************************************************************************80
!
!! R8MAT_SOLVE uses Gauss-Jordan elimination to solve an N by N linear system.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 August 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) RHS_NUM, the number of right hand sides.
!    RHS_NUM must be at least 0.
!
!    Input/output, real ( kind = 8 ) A(N,N+rhs_num), contains in rows and
!    columns 1 to N the coefficient matrix, and in columns N+1 through
!    N+rhs_num, the right hand sides.  On output, the coefficient matrix
!    area has been destroyed, while the right hand sides have
!    been overwritten with the corresponding solutions.
!
!    Output, integer ( kind = 4 ) INFO, singularity flag.
!    0, the matrix was not singular, the solutions were computed;
!    J, factorization failed on step J, and the solutions could not
!    be computed.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) rhs_num

  real ( kind = 8 ) a(n,n+rhs_num)
  real ( kind = 8 ) apivot
  real ( kind = 8 ) factor
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipivot
  integer ( kind = 4 ) j

  info = 0

  do j = 1, n
!
!  Choose a pivot row.
!
    ipivot = j
    apivot = a(j,j)

    do i = j+1, n
      if ( abs ( apivot ) < abs ( a(i,j) ) ) then
        apivot = a(i,j)
        ipivot = i
      end if
    end do

    if ( apivot == 0.0D+00 ) then
      info = j
      return
    end if
!
!  Interchange.
!
    do i = 1, n + rhs_num
      call r8_swap ( a(ipivot,i), a(j,i) )
    end do
!
!  A(J,J) becomes 1.
!
    a(j,j) = 1.0D+00
    a(j,j+1:n+rhs_num) = a(j,j+1:n+rhs_num) / apivot
!
!  A(I,J) becomes 0.
!
    do i = 1, n

      if ( i /= j ) then

        factor = a(i,j)
        a(i,j) = 0.0D+00
        a(i,j+1:n+rhs_num) = a(i,j+1:n+rhs_num) - factor * a(j,j+1:n+rhs_num)

      end if

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
!    Input, character ( len = * ) TITLE, an optional title.
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
!    Input, character ( len = * ) TITLE, an optional title.
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

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

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
subroutine r8vec_angle_3d ( u, v, angle )

!*****************************************************************************80
!
!! R8VEC_ANGLE_3D computes the angle between two vectors in 3D.
!
!  Modified:
!
!    07 July 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) U(3), V(3), the vectors.
!
!    Output, real ( kind = 8 ) ANGLE, the angle between the two vectors.
!
  implicit none

  real ( kind = 8 ) angle
  real ( kind = 8 ) angle_cos
  real ( kind = 8 ) arc_cosine
  real ( kind = 8 ) u(3)
  real ( kind = 8 ) u_norm
  real ( kind = 8 ) uv_dot
  real ( kind = 8 ) v(3)
  real ( kind = 8 ) v_norm

  uv_dot = dot_product ( u(1:3), v(1:3) )

  u_norm = sqrt ( dot_product ( u(1:3), u(1:3) ) )

  v_norm = sqrt ( dot_product ( v(1:3), v(1:3) ) )

  angle_cos = uv_dot / u_norm / v_norm

  angle = arc_cosine ( angle_cos )

  return
end
subroutine r8vec_cross_3d ( v1, v2, v3 )

!*****************************************************************************80
!
!! R8VEC_CROSS_3D computes the cross product of two vectors in 3D.
!
!  Discussion:
!
!    The cross product in 3D can be regarded as the determinant of the
!    symbolic matrix:
!
!          |  i  j  k |
!      det | x1 y1 z1 |
!          | x2 y2 z2 |
!
!      = ( y1 * z2 - z1 * y2 ) * i
!      + ( z1 * x2 - x1 * z2 ) * j
!      + ( x1 * y2 - y1 * x2 ) * k
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) V1(3), V2(3), the two vectors.
!
!    Output, real ( kind = 8 ) V3(3), the cross product vector.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  real ( kind = 8 ) v1(dim_num)
  real ( kind = 8 ) v2(dim_num)
  real ( kind = 8 ) v3(dim_num)

  v3(1) = v1(2) * v2(3) - v1(3) * v2(2)
  v3(2) = v1(3) * v2(1) - v1(1) * v2(3)
  v3(3) = v1(1) * v2(2) - v1(2) * v2(1)

  return
end
function r8vec_length ( dim_num, x )

!*****************************************************************************80
!
!! R8VEC_LENGTH returns the Euclidean length of a vector.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, real ( kind = 8 ) X(DIM_NUM), the vector.
!
!    Output, real ( kind = 8 ) R8VEC_LENGTH, the Euclidean length of the vector.
!
  implicit none

  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) r8vec_length
  real ( kind = 8 ) x(dim_num)

  r8vec_length = sqrt ( sum ( ( x(1:dim_num) )**2 ) )

  return
end
subroutine r8vec_print ( n, a, title )

!*****************************************************************************80
!
!! R8VEC_PRINT prints an R8VEC.
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
!    22 August 2000
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
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  character ( len = * )  title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i8,2x,g16.8)' ) i, a(i)
  end do

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

  character              c
  logical                ch_eqi
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

  logical                blank
  integer ( kind = 4 ) i
  integer ( kind = 4 ) lens
  integer ( kind = 4 ) nword
  character ( len = * )  s

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
subroutine tetrahedron_centroid_3d ( tetra, centroid )

!*****************************************************************************80
!
!! TETRAHEDRON_CENTROID_3D computes the centroid of a tetrahedron in 3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) TETRA(3,4) the tetrahedron vertices.
!
!    Output, real ( kind = 8 ) CENTROID(3), the coordinates of the centroid.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  real ( kind = 8 ) centroid(dim_num)
  integer ( kind = 4 ) i
  real ( kind = 8 ) tetra(dim_num,4)

  do i = 1, dim_num
    centroid(i) = sum ( tetra(i,1:4) ) / 4.0D+00
  end do

  return
end
subroutine tetrahedron_circumsphere_3d ( tetra, r, pc )

!*****************************************************************************80
!
!! TETRAHEDRON_CIRCUMSPHERE_3D computes the circumsphere of a tetrahedron in 3D.
!
!  Discussion:
!
!    The circumsphere, or circumscribed sphere, of a tetrahedron is the
!    sphere that passes through the four vertices.  The circumsphere is
!    not necessarily the smallest sphere that contains the tetrahedron.
!
!    Surprisingly, the diameter of the sphere can be found by solving
!    a 3 by 3 linear system.  This is because the vectors P2 - P1,
!    P3 - P1 and P4 - P1 are secants of the sphere, and each forms a
!    right triangle with the diameter through P1.  Hence, the dot product of
!    P2 - P1 with that diameter is equal to the square of the length
!    of P2 - P1, and similarly for P3 - P1 and P4 - P1.  This determines
!    the diameter vector originating at P1, and hence the radius and
!    center.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Adrian Bowyer, John Woodwark,
!    A Programmer's Geometry,
!    Butterworths, 1983.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) TETRA(3,4) the tetrahedron vertices.
!
!    Output, real ( kind = 8 ) R, PC(3), the center of the
!    circumscribed sphere, and its radius.  If the linear system is
!    singular, then R = -1, PC(1:3) = 0.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3
  integer ( kind = 4 ), parameter :: rhs_num = 1

  real ( kind = 8 ) a(dim_num,dim_num+rhs_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  real ( kind = 8 ) pc(dim_num)
  real ( kind = 8 ) r
  real ( kind = 8 ) tetra(dim_num,4)
!
!  Set up the linear system.
!
  a(1:dim_num,1:3) = transpose ( tetra(1:dim_num,2:4) )

  do j = 1, dim_num
    a(1:dim_num,j) = a(1:dim_num,j) - tetra(j,1)
  end do

  do i = 1, 3
    a(i,4) = sum ( a(i,1:3)**2 )
  end do
!
!  Solve the linear system.
!
  call r8mat_solve ( dim_num, rhs_num, a, info )
!
!  If the system was singular, return a consolation prize.
!
  if ( info /= 0 ) then
    r = -1.0D+00
    pc(1:dim_num) = 0.0D+00
    return
  end if
!
!  Compute the radius and center.
!
  r = 0.5D+00 * sqrt ( sum ( a(1:dim_num,4)**2 ) )

  pc(1:dim_num) = tetra(1:dim_num,1) + 0.5D+00 * a(1:dim_num,4)

  return
end
subroutine tetrahedron_dihedral_angles_3d ( tetra, angle )

!*****************************************************************************80
!
!! TETRAHEDRON_DIHEDRAL_ANGLES_3D computes dihedral angles of a tetrahedron.
!
!  Modified:
!
!    07 July 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) TETRA(3,4), the vertices of the tetrahedron.
!
!    Output, real ( kind = 8 ) ANGLE(6), the dihedral angles along the
!    axes AB, AC, AD, BC, BD and CD, respectively.
!
  implicit none

  real ( kind = 8 ) ab(3)
  real ( kind = 8 ) abc_normal(3)
  real ( kind = 8 ) abd_normal(3)
  real ( kind = 8 ) ac(3)
  real ( kind = 8 ) acd_normal(3)
  real ( kind = 8 ) ad(3)
  real ( kind = 8 ) angle(6)
  real ( kind = 8 ) bc(3)
  real ( kind = 8 ) bcd_normal(3)
  real ( kind = 8 ) bd(3)
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) tetra(3,4)

  ab(1:3) = tetra(1:3,2) - tetra(1:3,1)
  ac(1:3) = tetra(1:3,3) - tetra(1:3,1)
  ad(1:3) = tetra(1:3,4) - tetra(1:3,1)
  bc(1:3) = tetra(1:3,3) - tetra(1:3,2)
  bd(1:3) = tetra(1:3,4) - tetra(1:3,2)

  call r8vec_cross_3d ( ac, ab, abc_normal )
  call r8vec_cross_3d ( ab, ad, abd_normal )
  call r8vec_cross_3d ( ad, ac, acd_normal )
  call r8vec_cross_3d ( bc, bd, bcd_normal )

  call r8vec_angle_3d ( abc_normal, abd_normal, angle(1) )
  call r8vec_angle_3d ( abc_normal, acd_normal, angle(2) )
  call r8vec_angle_3d ( abd_normal, acd_normal, angle(3) )
  call r8vec_angle_3d ( abc_normal, bcd_normal, angle(4) )
  call r8vec_angle_3d ( abd_normal, bcd_normal, angle(5) )
  call r8vec_angle_3d ( acd_normal, bcd_normal, angle(6) )

  angle(1:6) = pi - angle(1:6)

  return
end
subroutine tetrahedron_edge_length_3d ( tetra, edge_length )

!*****************************************************************************80
!
!! TETRAHEDRON_EDGE_LENGTH_3D returns edge lengths of a tetrahedron in 3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) TETRA(3,4), the tetrahedron vertices.
!
!    Output, real ( kind = 8 ) EDGE_LENGTH(6), the length of the edges.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  real ( kind = 8 ) r8vec_length
  real ( kind = 8 ) edge_length(6)
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) k
  real ( kind = 8 ) tetra(dim_num,4)

  k = 0
  do j1 = 1, 3
    do j2 = j1+1, 4
      k = k + 1
      edge_length(k) = r8vec_length ( dim_num, &
        tetra(1:dim_num,j2) - tetra(1:dim_num,j1) )
     end do
  end do

  return
end
subroutine tetrahedron_face_angles_3d ( tetra, angles )

!*****************************************************************************80
!
!! TETRAHEDRON_FACE_ANGLES_3D returns the 12 face angles of a tetrahedron 3D.
!
!  Discussion:
!
!    The tetrahedron has 4 triangular faces.  This routine computes the
!    3 planar angles associated with each face.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 July 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) TETRA(3,4) the tetrahedron vertices.
!
!    Output, real ( kind = 8 ) ANGLES(3,4), the face angles.
!
  implicit none

  real ( kind = 8 ) angles(3,4)
  real ( kind = 8 ) tri(3,3)
  real ( kind = 8 ) tetra(3,4)
!
!  Face 123
!
  tri(1:3,1:3) = tetra(1:3,1:3)
  call triangle_angles_3d ( tri, angles(1:3,1) )
!
!  Face 124
!
  tri(1:3,1:2) = tetra(1:3,1:2)
  tri(1:3,3) = tetra(1:3,4)
  call triangle_angles_3d ( tri, angles(1:3,2) )
!
!  Face 134
!
  tri(1:3,1) = tetra(1:3,1)
  tri(1:3,2:3) = tetra(1:3,3:4)
  call triangle_angles_3d ( tri, angles(1:3,3) )
!
!  Face 234
!
  tri(1:3,1:3) = tetra(1:3,2:4)
  call triangle_angles_3d ( tri, angles(1:3,4) )

  return
end
subroutine tetrahedron_face_areas_3d ( tetra, areas )

!*****************************************************************************80
!
!! TETRAHEDRON_FACE_AREAS_3D returns the 4 face areas of a tetrahedron 3D.
!
!  Discussion:
!
!    The tetrahedron has 4 triangular faces.  This routine computes the
!    area of each face.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 July 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) TETRA(3,4) the tetrahedron vertices.
!
!    Output, real ( kind = 8 ) AREAS(4), the face areas.
!
  implicit none

  real ( kind = 8 ) areas(4)
  real ( kind = 8 ) tri(3,3)
  real ( kind = 8 ) tetra(3,4)
!
!  Face 123
!
  tri(1:3,1:3) = tetra(1:3,1:3)
  call triangle_area_3d ( tri, areas(1) )
!
!  Face 124
!
  tri(1:3,1:2) = tetra(1:3,1:2)
  tri(1:3,3) = tetra(1:3,4)
  call triangle_area_3d ( tri, areas(2) )
!
!  Face 134
!
  tri(1:3,1) = tetra(1:3,1)
  tri(1:3,2:3) = tetra(1:3,3:4)
  call triangle_area_3d ( tri, areas(3) )
!
!  Face 234
!
  tri(1:3,1:3) = tetra(1:3,2:4)
  call triangle_area_3d ( tri, areas(4) )

  return
end
subroutine tetrahedron_insphere_3d ( tetra, r, pc )

!*****************************************************************************80
!
!! TETRAHEDRON_INSPHERE_3D finds the insphere of a tetrahedron in 3D.
!
!  Discussion:
!
!    The insphere of a tetrahedron is the inscribed sphere, which touches
!    each face of the tetrahedron at a single point.
!
!    The points of contact are the centroids of the triangular faces
!    of the tetrahedron.  Therefore, the point of contact for a face
!    can be computed as the average of the vertices of that face.
!
!    The sphere can then be determined as the unique sphere through
!    the four given centroids.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Philip Schneider, David Eberly,
!    Geometric Tools for Computer Graphics,
!    Elsevier, 2002,
!    ISBN: 1558605940,
!    LC: T385.G6974.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) TETRA(3,4), the vertices of the tetrahedron.
!
!    Output, real ( kind = 8 ) R, PC(3), the radius and the center
!    of the sphere.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  real ( kind = 8 ) b(4,4)
  real ( kind = 8 ) r8mat_det_4d
  real ( kind = 8 ) r8vec_length
  real ( kind = 8 ) gamma
  real ( kind = 8 ) l123
  real ( kind = 8 ) l124
  real ( kind = 8 ) l134
  real ( kind = 8 ) l234
  real ( kind = 8 ) n123(1:dim_num)
  real ( kind = 8 ) n124(1:dim_num)
  real ( kind = 8 ) n134(1:dim_num)
  real ( kind = 8 ) n234(1:dim_num)
  real ( kind = 8 ) pc(1:dim_num)
  real ( kind = 8 ) r
  real ( kind = 8 ) tetra(1:dim_num,4)
  real ( kind = 8 ) v21(1:dim_num)
  real ( kind = 8 ) v31(1:dim_num)
  real ( kind = 8 ) v41(1:dim_num)
  real ( kind = 8 ) v32(1:dim_num)
  real ( kind = 8 ) v42(1:dim_num)
  real ( kind = 8 ) v43(1:dim_num)

  v21(1:dim_num) = tetra(1:dim_num,2) - tetra(1:dim_num,1)
  v31(1:dim_num) = tetra(1:dim_num,3) - tetra(1:dim_num,1)
  v41(1:dim_num) = tetra(1:dim_num,4) - tetra(1:dim_num,1)
  v32(1:dim_num) = tetra(1:dim_num,3) - tetra(1:dim_num,2)
  v42(1:dim_num) = tetra(1:dim_num,4) - tetra(1:dim_num,2)
  v43(1:dim_num) = tetra(1:dim_num,4) - tetra(1:dim_num,3)

  call r8vec_cross_3d ( v21, v31, n123 )
  call r8vec_cross_3d ( v41, v21, n124 )
  call r8vec_cross_3d ( v31, v41, n134 )
  call r8vec_cross_3d ( v42, v32, n234 )

  l123 = r8vec_length ( dim_num, n123 )
  l124 = r8vec_length ( dim_num, n124 )
  l134 = r8vec_length ( dim_num, n134 )
  l234 = r8vec_length ( dim_num, n234 )

  pc(1:dim_num) = ( l234 * tetra(1:dim_num,1)   &
                  + l134 * tetra(1:dim_num,2)   &
                  + l124 * tetra(1:dim_num,3)   &
                  + l123 * tetra(1:dim_num,4) ) &
                / ( l234 + l134 + l124 + l123 )

  b(1:dim_num,1:4) = tetra(1:dim_num,1:4)
  b(4,1:4) = 1.0D+00

  gamma = abs ( r8mat_det_4d ( b ) )

! gamma = abs ( &
!     ( tetra(1,2) * tetra(2,3) * tetra(3,4) &
!     - tetra(1,3) * tetra(2,4) * tetra(3,2) &
!     + tetra(1,4) * tetra(2,2) * tetra(3,3) ) &
!   - ( tetra(1,1) * tetra(2,3) * tetra(3,4) &
!     - tetra(1,3) * tetra(2,4) * tetra(3,1) &
!     + tetra(1,4) * tetra(2,1) * tetra(3,3) ) &
!   + ( tetra(1,1) * tetra(2,2) * tetra(3,4) &
!     - tetra(1,2) * tetra(2,4) * tetra(3,1) &
!     + tetra(1,4) * tetra(2,1) * tetra(3,2) ) &
!   - ( tetra(1,1) * tetra(2,2) * tetra(3,3) &
!     - tetra(1,2) * tetra(2,3) * tetra(3,1) &
!     + tetra(1,3) * tetra(2,1) * tetra(3,2) ) )

  r = gamma / ( l234 + l134 + l124 + l123 )

  return
end
subroutine tetrahedron_quality1_3d ( tetra, quality )

!*****************************************************************************80
!
!! TETRAHEDRON_QUALITY1_3D: "quality" of a tetrahedron in 3D.
!
!  Discussion:
!
!    The quality of a tetrahedron is 3 times the ratio of the radius of
!    the inscribed sphere divided by that of the circumscribed sphere.
!
!    An equilateral tetrahredron achieves the maximum possible quality of 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) TETRA(3,4), the tetrahedron vertices.
!
!    Output, real ( kind = 8 ) QUALITY, the quality of the tetrahedron.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  real ( kind = 8 ) pc(dim_num)
  real ( kind = 8 ) quality
  real ( kind = 8 ) r_in
  real ( kind = 8 ) r_out
  real ( kind = 8 ) tetra(dim_num,4)

  call tetrahedron_circumsphere_3d ( tetra, r_out, pc )

  call tetrahedron_insphere_3d ( tetra, r_in, pc )

  quality = 3.0D+00 * r_in / r_out

  return
end
subroutine tetrahedron_quality2_3d ( tetra, quality2 )

!*****************************************************************************80
!
!! TETRAHEDRON_QUALITY2_3D: "quality" of a tetrahedron in 3D.
!
!  Discussion:
!
!    The quality measure #2 of a tetrahedron is:
!
!      QUALITY2 = 2 * sqrt ( 6 ) * RIN / LMAX
!
!    where
!
!      RIN = radius of the inscribed sphere;
!      LMAX = length of longest side of the tetrahedron.
!
!    An equilateral tetrahredron achieves the maximum possible quality of 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Qiang Du, Desheng Wang,
!    The Optimal Centroidal Voronoi Tesselations and the Gersho's
!    Conjecture in the Three-Dimensional Space,
!    Computers and Mathematics with Applications,
!    Volume 49, 2005, pages 1355-1373.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) TETRA(3,4), the tetrahedron vertices.
!
!    Output, real ( kind = 8 ) QUALITY2, the quality of the tetrahedron.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  real ( kind = 8 ) edge_length(6)
  real ( kind = 8 ) l_max
  real ( kind = 8 ) pc(dim_num)
  real ( kind = 8 ) quality2
  real ( kind = 8 ) r_in
  real ( kind = 8 ) tetra(dim_num,4)

  call tetrahedron_edge_length_3d ( tetra, edge_length )

  l_max = maxval ( edge_length(1:6) )

  call tetrahedron_insphere_3d ( tetra, r_in, pc )

  quality2 = 2.0D+00 * sqrt ( 6.0D+00 ) * r_in / l_max

  return
end
subroutine tetrahedron_quality3_3d ( tetra, quality3 )

!*****************************************************************************80
!
!! TETRAHEDRON_QUALITY3_3D computes the mean ratio of a tetrahedron.
!
!  Discussion:
!
!    This routine computes QUALITY3, the eigenvalue or mean ratio of
!    a tetrahedron.
!
!      QUALITY3 = 12 * ( 3 * volume )**(2/3) / (sum of squares of edge lengths).
!
!    This value may be used as a shape quality measure for the tetrahedron.
!
!    For an equilateral tetrahedron, the value of this quality measure
!    will be 1.  For any other tetrahedron, the value will be between
!    0 and 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 August 2005
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) TETRA(3,4), the vertices of the tetrahedron.
!
!    Output, real ( kind = 8 ) QUALITY3, the mean ratio of the tetrahedron.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  real ( kind = 8 ) ab(dim_num)
  real ( kind = 8 ) ac(dim_num)
  real ( kind = 8 ) ad(dim_num)
  real ( kind = 8 ) bc(dim_num)
  real ( kind = 8 ) bd(dim_num)
  real ( kind = 8 ) cd(dim_num)
  real ( kind = 8 ) denom
  real ( kind = 8 ) lab
  real ( kind = 8 ) lac
  real ( kind = 8 ) lad
  real ( kind = 8 ) lbc
  real ( kind = 8 ) lbd
  real ( kind = 8 ) lcd
  real ( kind = 8 ) quality3
  real ( kind = 8 ) tetra(dim_num,4)
  real ( kind = 8 ) volume
!
!  Compute the vectors representing the sides of the tetrahedron.
!
  ab(1:3) = tetra(1:dim_num,2) - tetra(1:dim_num,1)
  ac(1:3) = tetra(1:dim_num,3) - tetra(1:dim_num,1)
  ad(1:3) = tetra(1:dim_num,4) - tetra(1:dim_num,1)
  bc(1:3) = tetra(1:dim_num,3) - tetra(1:dim_num,2)
  bd(1:3) = tetra(1:dim_num,4) - tetra(1:dim_num,2)
  cd(1:3) = tetra(1:dim_num,4) - tetra(1:dim_num,3)
!
!  Compute the squares of the lengths of the sides.
!
  lab = sum ( ab(1:dim_num)**2 )
  lac = sum ( ac(1:dim_num)**2 )
  lad = sum ( ad(1:dim_num)**2 )
  lbc = sum ( bc(1:dim_num)**2 )
  lbd = sum ( bd(1:dim_num)**2 )
  lcd = sum ( cd(1:dim_num)**2 )
!
!  Compute the volume.
!
  volume = abs ( &
      ab(1) * ( ac(2) * ad(3) - ac(3) * ad(2) ) &
    + ab(2) * ( ac(3) * ad(1) - ac(1) * ad(3) ) &
    + ab(3) * ( ac(1) * ad(2) - ac(2) * ad(1) ) ) / 6.0D+00

  denom = lab + lac + lad + lbc + lbd + lcd

  if ( denom == 0.0D+00 ) then
    quality3 = 0.0D+00
  else
    quality3 = 12.0D+00 * ( 3.0D+00 * volume )**( 2.0D+00 / 3.0D+00 ) / denom
  end if

  return
end
subroutine tetrahedron_quality4_3d ( tetra, quality4 )

!*****************************************************************************80
!
!! TETRAHEDRON_QUALITY4_3D computes the minimum solid angle of a tetrahedron.
!
!  Discussion:
!
!    This routine computes a quality measure for a tetrahedron, based
!    on the sine of half the minimum of the four solid angles.
!
!    The quality measure for an equilateral tetrahedron should be 1,
!    since the solid angles of such a tetrahedron are each equal to pi.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 August 2005
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) TETRA(3,4), the vertices of the tetrahedron.
!
!    Output, real ( kind = 8 ) QUALITY4, the value of the quality measure.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  real ( kind = 8 ) ab(dim_num)
  real ( kind = 8 ) ac(dim_num)
  real ( kind = 8 ) ad(dim_num)
  real ( kind = 8 ) bc(dim_num)
  real ( kind = 8 ) bd(dim_num)
  real ( kind = 8 ) cd(dim_num)
  real ( kind = 8 ) denom
  real ( kind = 8 ) l1
  real ( kind = 8 ) l2
  real ( kind = 8 ) l3
  real ( kind = 8 ) lab
  real ( kind = 8 ) lac
  real ( kind = 8 ) lad
  real ( kind = 8 ) lbc
  real ( kind = 8 ) lbd
  real ( kind = 8 ) lcd
  real ( kind = 8 ) quality4
  real ( kind = 8 ) tetra(dim_num,4)
  real ( kind = 8 ) volume
!
!  Compute the vectors that represent the sides.
!
  ab(1:dim_num) = tetra(1:dim_num,2) - tetra(1:dim_num,1)
  ac(1:dim_num) = tetra(1:dim_num,3) - tetra(1:dim_num,1)
  ad(1:dim_num) = tetra(1:dim_num,4) - tetra(1:dim_num,1)
  bc(1:dim_num) = tetra(1:dim_num,3) - tetra(1:dim_num,2)
  bd(1:dim_num) = tetra(1:dim_num,4) - tetra(1:dim_num,2)
  cd(1:dim_num) = tetra(1:dim_num,4) - tetra(1:dim_num,3)
!
!  Compute the lengths of the sides.
!
  lab = sqrt ( sum ( ab(1:dim_num)**2 ) )
  lac = sqrt ( sum ( ac(1:dim_num)**2 ) )
  lad = sqrt ( sum ( ad(1:dim_num)**2 ) )
  lbc = sqrt ( sum ( bc(1:dim_num)**2 ) )
  lbd = sqrt ( sum ( bd(1:dim_num)**2 ) )
  lcd = sqrt ( sum ( cd(1:dim_num)**2 ) )
!
!  Compute the volume
!
  volume = abs ( &
      ab(1) * ( ac(2) * ad(3) - ac(3) * ad(2) ) &
    + ab(2) * ( ac(3) * ad(1) - ac(1) * ad(3) ) &
    + ab(3) * ( ac(1) * ad(2) - ac(2) * ad(1) ) ) / 6.0D+00

  quality4 = 1.0D+00

  l1 = lab + lac
  l2 = lab + lad
  l3 = lac + lad

  denom = ( l1 + lbc ) * ( l1 - lbc ) &
        * ( l2 + lbd ) * ( l2 - lbd ) &
        * ( l3 + lcd ) * ( l3 - lcd )

  if ( denom <= 0.0D+00 ) then
    quality4 = 0.0D+00
  else
    quality4 = min ( quality4, 12.0D+00 * volume / sqrt ( denom ) )
  end if

  l1 = lab + lbc
  l2 = lab + lbd
  l3 = lbc + lbd

  denom = ( l1 + lac ) * ( l1 - lac ) &
        * ( l2 + lad ) * ( l2 - lad ) &
        * ( l3 + lcd ) * ( l3 - lcd )

  if ( denom <= 0.0D+00 ) then
    quality4 = 0.0D+00
  else
    quality4 = min ( quality4, 12.0D+00 * volume / sqrt ( denom ) )
  end if

  l1 = lac + lbc
  l2 = lac + lcd
  l3 = lbc + lcd

  denom = ( l1 + lab ) * ( l1 - lab ) &
        * ( l2 + lad ) * ( l2 - lad ) &
        * ( l3 + lbd ) * ( l3 - lbd )

  if ( denom <= 0.0D+00 ) then
    quality4 = 0.0D+00
  else
    quality4 = min ( quality4, 12.0D+00 * volume / sqrt ( denom ) )
  end if

  l1 = lad + lbd
  l2 = lad + lcd
  l3 = lbd + lcd

  denom = ( l1 + lab ) * ( l1 - lab ) &
        * ( l2 + lac ) * ( l2 - lac ) &
        * ( l3 + lbc ) * ( l3 - lbc )

  if ( denom <= 0.0D+00 ) then
    quality4 = 0.0D+00
  else
    quality4 = min ( quality4, 12.0D+00 * volume / sqrt ( denom ) )
  end if

  quality4 = quality4 * 1.5D+00 * sqrt ( 6.0D+00 )

  return
end
subroutine tetrahedron_solid_angles_3d ( tetra, angle )

!*****************************************************************************80
!
!! TETRAHEDRON_SOLID_ANGLES_3D computes solid angles of a tetrahedron.
!
!  Modified:
!
!    07 July 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) TETRA(3,4), the vertices of the tetrahedron.
!
!    Output, real ( kind = 8 ) ANGLE(4), the solid angles.
!
  implicit none

  real ( kind = 8 ) angle(4)
  real ( kind = 8 ) dihedral_angles(6)
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) tetra(3,4)

  call tetrahedron_dihedral_angles_3d ( tetra, dihedral_angles )

  angle(1) = dihedral_angles(1) + dihedral_angles(2) + dihedral_angles(3) - pi
  angle(2) = dihedral_angles(1) + dihedral_angles(4) + dihedral_angles(5) - pi
  angle(3) = dihedral_angles(2) + dihedral_angles(4) + dihedral_angles(6) - pi
  angle(4) = dihedral_angles(3) + dihedral_angles(5) + dihedral_angles(6) - pi

  return
end
subroutine tetrahedron_volume_3d ( tetra, volume )

!*****************************************************************************80
!
!! TETRAHEDRON_VOLUME_3D computes the volume of a tetrahedron in 3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) TETRA(3,4), the vertices of the tetrahedron.
!
!    Output, real ( kind = 8 ) VOLUME, the volume of the tetrahedron.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  real ( kind = 8 ) a(4,4)
  real ( kind = 8 ) r8mat_det_4d
  real ( kind = 8 ) tetra(dim_num,4)
  real ( kind = 8 ) volume

  a(1:dim_num,1:4) = tetra(1:dim_num,1:4)
  a(4,1:4) = 1.0D+00

  volume = abs ( r8mat_det_4d ( a ) ) / 6.0D+00

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
subroutine triangle_angles_3d ( t, angle )

!*****************************************************************************80
!
!! TRIANGLE_ANGLES_3D computes the angles of a triangle in 3D.
!
!  Discussion:
!
!    The law of cosines is used:
!
!      C * C = A * A + B * B - 2 * A * B * COS ( GAMMA )
!
!    where GAMMA is the angle opposite side C.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T(3,3), the triangle vertices.
!
!    Output, real ( kind = 8 ) ANGLE(3), the angles opposite
!    sides P1-P2, P2-P3 and P3-P1, in radians.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  real ( kind = 8 ) a
  real ( kind = 8 ) angle(3)
  real ( kind = 8 ) arc_cosine
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) t(dim_num,3)
!
!  Compute the length of each side.
!
  a = sqrt ( sum ( ( t(1:dim_num,1) - t(1:dim_num,2) )**2 ) )
  b = sqrt ( sum ( ( t(1:dim_num,2) - t(1:dim_num,3) )**2 ) )
  c = sqrt ( sum ( ( t(1:dim_num,3) - t(1:dim_num,1) )**2 ) )
!
!  Take care of a ridiculous special case.
!
  if ( a == 0.0D+00 .and. b == 0.0D+00 .and. c == 0.0D+00 ) then
    angle(1:3) = 2.0D+00 * pi / 3.0D+00
    return
  end if

  if ( c == 0.0D+00 .or. a == 0.0D+00 ) then
    angle(1) = pi
  else
    angle(1) = arc_cosine ( ( c * c + a * a - b * b ) / ( 2.0D+00 * c * a ) )
  end if

  if ( a == 0.0D+00 .or. b == 0.0D+00 ) then
    angle(2) = pi
  else
    angle(2) = arc_cosine ( ( a * a + b * b - c * c ) / ( 2.0D+00 * a * b ) )
  end if

  if ( b == 0.0D+00 .or. c == 0.0D+00 ) then
    angle(3) = pi
  else
    angle(3) = arc_cosine ( ( b * b + c * c - a * a ) / ( 2.0D+00 * b * c ) )
  end if

  return
end
subroutine triangle_area_3d ( t, area )

!*****************************************************************************80
!
!! TRIANGLE_AREA_3D computes the area of a triangle in 3D.
!
!  Discussion:
!
!    This routine uses the fact that the norm of the cross product
!    of two vectors is the area of the parallelogram they form.
!
!    Therefore, the area of the triangle is half of that value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Adrian Bowyer, John Woodwark,
!    A Programmer's Geometry,
!    Butterworths, 1983.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T(3,3), the triangle vertices.
!
!    Output, real ( kind = 8 ) AREA, the area of the triangle.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  real ( kind = 8 ) area
  real ( kind = 8 ) cross(dim_num)
  real ( kind = 8 ) t(dim_num,3)
!
!  Compute the cross product vector.
!
  cross(1) = ( t(2,2) - t(2,1) ) * ( t(3,3) - t(3,1) ) &
           - ( t(3,2) - t(3,1) ) * ( t(2,3) - t(2,1) )

  cross(2) = ( t(3,2) - t(3,1) ) * ( t(1,3) - t(1,1) ) &
           - ( t(1,2) - t(1,1) ) * ( t(3,3) - t(3,1) )

  cross(3) = ( t(1,2) - t(1,1) ) * ( t(2,3) - t(2,1) ) &
           - ( t(2,2) - t(2,1) ) * ( t(1,3) - t(1,1) )

  area = 0.5D+00 * sqrt ( sum ( cross(1:3)**2 ) )

  return
end
