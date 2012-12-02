program main

!*****************************************************************************80
!
!! MAIN is the main program for TRIANGULATION_REFINE.
!
!  Discussion:
!
!    TRIANGULATION_REFINE refines a triangulation by doubling.
!
!  Usage:
!
!    triangulation_refine prefix
!
!    where 'prefix' is the common filename prefix:
!
!    * prefix_nodes.txt contains the node coordinates,
!    * prefix_elements.txt contains the element definitions.
!    * prefix_ref_nodes.txt will contain the refined nodes;
!    * prefix_ref_elements.txt will contain the refined elements.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 October 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) arg_num
  logical, parameter :: debug = .true.
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: edge_data
  integer ( kind = 4 ) iarg
  integer ( kind = 4 ) iargc
  integer ( kind = 4 ) :: ierror = 0
  character ( len = 255 ) :: node_filename = ' '
  character ( len = 255 ) :: element_filename = ' '
  integer ( kind = 4 ) node_num1
  integer ( kind = 4 ) node_num2
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: node_xy1
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: node_xy2
  character ( len = 255 ) :: node_ref_filename = ' '
  character ( len = 255 ) :: element_ref_filename = ' '
  character ( len = 255 ) prefix
  integer ( kind = 4 ) triangle_num1
  integer ( kind = 4 ) triangle_num2
  integer ( kind = 4 ) triangle_order
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: triangle_node1
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: triangle_node2

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TRIANGULATION_REFINE'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Read a "linear" or "quadratic" triangulation'
  write ( *, '(a)' ) '  and write out a refined triangulation.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In particular:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Read a dataset of NODE_NUM1 points in 2 dimensions.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Read a triangulation dataset of TRIANGLE_NUM1'
  write ( *, '(a)' ) '  triangles using 3 or 6 nodes per triangle.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Subdivide each triangle into 4 triangles,'
  write ( *, '(a)' ) '  generate new nodes as midpoints of current nodes.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Write out the new node and triangulation data.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  If the input triangulation was Delaunay, then'
  write ( *, '(a)' ) '  the output triangulation will be Delaunay.'
!
!  Get the number of command line arguments.
!
  arg_num = iargc ( )

!
!  Argument 1 is the common filename prefix.
!
  if ( 1 <= arg_num ) then

    iarg = 1
    call getarg ( iarg, prefix )

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRIANGULATION_REFINE:'
    write ( *, '(a)' ) '  Please enter the filename prefix.'

    read ( *, '(a)' ) prefix

  end if
!
!  Create the filenames.
!
  node_filename = trim ( prefix ) // '_nodes.txt'
  element_filename = trim ( prefix ) // '_elements.txt'
  node_ref_filename = trim ( prefix ) // '_ref_nodes.txt'
  element_ref_filename = trim ( prefix ) // '_ref_elements.txt'
!
!  Read the node data.
!
  call r8mat_header_read ( node_filename, dim_num, node_num1 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Read the header of "' &
    // trim ( node_filename ) //'".'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Spatial dimension DIM_NUM = ', dim_num
  write ( *, '(a,i8)' ) '  Number of points NODE_NUM1  = ', node_num1

  if ( dim_num /= 2 ) then
    write ( *, * ) ' '
    write ( *, '(a)' ) 'TRIANGULATION_REFINE - Fatal error!'
    write ( *, '(a)' ) '  Spatial dimension DIM_NUM must be 2.'
    stop
  end if

  allocate ( node_xy1(1:dim_num,1:node_num1) )

  call r8mat_data_read ( node_filename, dim_num, node_num1, node_xy1 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Read the data in "' &
    // trim ( node_filename ) //'".'

  call r8mat_transpose_print_some ( dim_num, node_num1, node_xy1, 1, 1, &
    dim_num, 5, '  First 5 nodes:' )
!
!  Read the element data.
!
  call i4mat_header_read ( element_filename, triangle_order, &
    triangle_num1 )

  if ( triangle_order /= 3 .and. triangle_order /= 6 ) then
    write ( *, * ) ' '
    write ( *, '(a)' ) 'TRIANGULATION_REFINE - Fatal error!'
    write ( *, '(a)' ) '  Data is not for a 3 node or 6 node triangulation.'
    stop
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Read the header of "' &
    // trim ( element_filename ) //'".'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Triangle order = ', triangle_order
  write ( *, '(a,i8)' ) '  Number of triangles TRIANGLE_NUM1  = ', triangle_num1

  allocate ( triangle_node1(1:triangle_order,1:triangle_num1) )

  call i4mat_data_read ( element_filename, triangle_order, &
    triangle_num1, triangle_node1 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Read the data in "' &
    // trim ( element_filename ) //'".'

  call i4mat_transpose_print_some ( triangle_order, triangle_num1, &
    triangle_node1, 1, 1, triangle_order, 5, '  First 5 triangles:' )
!
!  Detect and correct 0-based indexing.
!
  call mesh_base_one ( node_num1, triangle_order, triangle_num1, triangle_node1 )
!
!  Determine the size of the refined mesh.
!
  allocate ( edge_data(5,3*triangle_num1) )

  if ( triangle_order == 3 ) then

    call triangulation_order3_refine_size ( node_num1, triangle_num1, &
      triangle_node1, node_num2, triangle_num2, edge_data )

  else if ( triangle_order == 6 ) then

    call triangulation_order6_refine_size ( node_num1, triangle_num1, &
      triangle_node1, node_num2, triangle_num2, edge_data )

  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) &
    '  Number of nodes in refined mesh =      ', node_num2
  write ( *, '(a,i8)' ) &
    '  Number of triangles in refined mesh =  ', triangle_num2

  allocate ( node_xy2(1:dim_num,1:node_num2) )
  allocate ( triangle_node2(1:triangle_order,1:triangle_num2) )
!
!  Compute the refined mesh.
!
  if ( triangle_order == 3 ) then

    call triangulation_order3_refine_compute ( node_num1, triangle_num1, &
      node_xy1, triangle_node1, node_num2, triangle_num2, edge_data, &
      node_xy2, triangle_node2 )

  else if ( triangle_order == 6 ) then

    call triangulation_order6_refine_compute ( node_num1, triangle_num1, &
      node_xy1, triangle_node1, node_num2, triangle_num2, edge_data, &
      node_xy2, triangle_node2 )

  end if

  if ( debug ) then

    call r8mat_transpose_print_some ( dim_num, node_num2, node_xy2, &
      1, 1, dim_num, 5, '  First 5 output nodes:' )

    call i4mat_transpose_print_some ( triangle_order, triangle_num2, &
      triangle_node2, 1, 1, triangle_order, 5, '  First 5 output triangles' )

  end if
!
!  Write out the node and triangle data.
!
  call r8mat_write ( node_ref_filename, dim_num, node_num2, node_xy2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Wrote the refined node data to "' &
    // trim ( node_ref_filename ) //'".'

  call i4mat_write ( element_ref_filename, triangle_order, &
    triangle_num2, triangle_node2 )

  write ( *, '(a)' ) '  Wrote the refined element data to "' &
    // trim ( element_ref_filename ) //'".'
!
!  Release memory.
!
  deallocate ( edge_data )
  deallocate ( node_xy1 )
  deallocate ( node_xy2 )
  deallocate ( triangle_node1 )
  deallocate ( triangle_node2 )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TRIANGULATION_REFINE'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
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
!    Output, integer ( kind = 4 ) DIGIT, the corresponding integer value.  If C was
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
  integer ( kind = 4 ) input_unit
  integer ( kind = 4 ) ios
  character ( len = 255 ) line
!
!  Open the file.
!
  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_filename, status = 'old', &
    form = 'formatted', access = 'sequential', iostat = ios )

  if ( ios /= 0 ) then
    column_num = -1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_COLUMN_COUNT - Fatal error!'
    write ( *, '(a)' ) '  Could not open the file:'
    write ( *, '(a)' ) '    ' // trim ( input_filename )
    return
  end if
!
!  Read one line, but skip blank lines and comment lines.
!
  got_one = .false.

  do

    read ( input_unit, '(a)', iostat = ios ) line

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

    rewind ( input_unit )

    do

      read ( input_unit, '(a)', iostat = ios ) line

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
  integer ( kind = 4 ) input_unit
  integer ( kind = 4 ) ios
  character ( len = 255 ) line
  integer ( kind = 4 ) record_num
  integer ( kind = 4 ) row_num

  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_filename, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    row_num = -1;
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_ROW_COUNT - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input file: ' // &
      trim ( input_filename )
    stop
  end if

  comment_num = 0
  row_num = 0
  record_num = 0
  bad_num = 0

  do

    read ( input_unit, '(a)', iostat = ios ) line

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
subroutine i4col_compare ( m, n, a, i, j, isgn )

!*****************************************************************************80
!
!! I4COL_COMPARE compares columns I and J of an I4COL.
!
!  Example:
!
!    Input:
!
!      M = 3, N = 4, I = 2, J = 4
!
!      A = (
!        1  2  3  4
!        5  6  7  8
!        9 10 11 12 )
!
!    Output:
!
!      ISGN = -1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, integer ( kind = 4 ) A(M,N), an array of N columns of vectors of length M.
!
!    Input, integer ( kind = 4 ) I, J, the columns to be compared.
!    I and J must be between 1 and N.
!
!    Output, integer ( kind = 4 ) ISGN, the results of the comparison:
!    -1, column I < column J,
!     0, column I = column J,
!    +1, column J < column I.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
!
!  Check.
!
  if ( i < 1 .or. n < i ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4COL_COMPARE - Fatal error!'
    write ( *, '(a)' ) '  Column index I is out of bounds.'
    stop
  end if

  if ( j < 1 .or. n < j ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4COL_COMPARE - Fatal error!'
    write ( *, '(a)' ) '  Column index J is out of bounds.'
    stop
  end if

  isgn = 0

  if ( i == j ) then
    return
  end if

  k = 1

  do while ( k <= m )

    if ( a(k,i) < a(k,j) ) then
      isgn = -1
      return
    else if ( a(k,j) < a(k,i) ) then
      isgn = +1
      return
    end if

    k = k + 1

  end do

  return
end
subroutine i4col_sort_a ( m, n, a )

!*****************************************************************************80
!
!! I4COL_SORT_A ascending sorts an I4COL.
!
!  Discussion:
!
!    In lexicographic order, the statement "X < Y", applied to two real
!    vectors X and Y of length M, means that there is some index I, with
!    1 <= I <= M, with the property that
!
!      X(J) = Y(J) for J < I,
!    and
!      X(I) < Y(I).
!
!    In other words, the first time they differ, X is smaller.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 September 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of A, and the length of
!    a vector of data.
!
!    Input, integer ( kind = 4 ) N, the number of columns of A.
!
!    Input/output, integer ( kind = 4 ) A(M,N).
!    On input, the array of N columns of M-vectors.
!    On output, the columns of A have been sorted in ascending
!    lexicographic order.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j

  if ( m <= 0 ) then
    return
  end if

  if ( n <= 1 ) then
    return
  end if
!
!  Initialize.
!
  i = 0
  indx = 0
  isgn = 0
  j = 0
!
!  Call the external heap sorter.
!
  do

    call sort_heap_external ( n, indx, i, j, isgn )
!
!  Interchange the I and J objects.
!
    if ( 0 < indx ) then

      call i4col_swap ( m, n, a, i, j )
!
!  Compare the I and J objects.
!
    else if ( indx < 0 ) then

      call i4col_compare ( m, n, a, i, j, isgn )

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
subroutine i4col_sorted_unique_count ( m, n, a, unique_num )

!*****************************************************************************80
!
!! I4COL_SORTED_UNIQUE_COUNT counts unique elements in an I4COL.
!
!  Discussion:
!
!    The columns of the array may be ascending or descending sorted.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, integer ( kind = 4 ) A(M,N), a sorted array, containing
!    N columns of data.
!
!    Output, integer ( kind = 4 ) UNIQUE_NUM, the number of unique columns.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) unique_num

  if ( n <= 0 ) then
    unique_num = 0
    return
  end if

  unique_num = 1
  j1 = 1

  do j2 = 2, n

    if ( any ( a(1:m,j1) /= a(1:m,j2) ) ) then
      unique_num = unique_num + 1
      j1 = j2
    end if

  end do

  return
end
subroutine i4col_swap ( m, n, a, i, j )

!*****************************************************************************80
!
!! I4COL_SWAP swaps columns I and J of an I4COL.
!
!  Example:
!
!    Input:
!
!      M = 3, N = 4, I = 2, J = 4
!
!      A = (
!        1  2  3  4
!        5  6  7  8
!        9 10 11 12 )
!
!    Output:
!
!      A = (
!        1  4  3  2
!        5  8  7  6
!        9 12 11 10 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns in the array.
!
!    Input, integer ( kind = 4 ) A(M,N), an array of N columns of length M.
!
!    Input, integer ( kind = 4 ) I, J, the columns to be swapped.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) col(m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  if ( i < 1 .or. n < i .or. j < 1 .or. n < j ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4COL_SWAP - Fatal error!'
    write ( *, '(a)' ) '  I or J is out of bounds.'
    write ( *, '(a,i8)' ) '  I =    ', i
    write ( *, '(a,i8)' ) '  J =    ', j
    write ( *, '(a,i8)' ) '  N =    ', n
    stop

  end if

  if ( i == j ) then
    return
  end if

  col(1:m) = a(1:m,i)
  a(1:m,i) = a(1:m,j)
  a(1:m,j) = col(1:m)

  return
end
subroutine i4mat_data_read ( input_filename, m, n, table )

!*****************************************************************************80
!
!! I4MAT_DATA_READ reads data from an I4MAT file.
!
!  Discussion:
!
!    The file may contain more than N points, but this routine
!    will return after reading N points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2005
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
!    Output, integer ( kind = 4 ) TABLE(M,N), the table data.
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
  integer ( kind = 4 ) table(m,n)
  integer ( kind = 4 ) x(m)

  ierror = 0

  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_filename, status = 'old', &
    iostat = input_status )

  if ( input_status /= 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4MAT_DATA_READ - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the input file "' // &
      trim ( input_filename ) // '" on unit ', input_unit
    stop
  end if

  j = 0

  do while ( j < n )

    read ( input_unit, '(a)', iostat = input_status ) line

    if ( input_status /= 0 ) then
      ierror = 2
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'I4MAT_DATA_READ - Fatal error!'
      write ( *, '(a)' ) '  Error while reading lines of data.'
      write ( *, '(a,i8)' ) '  Number of values expected per line M = ', m
      write ( *, '(a,i8)' ) '  Number of data lines read, J =         ', j
      write ( *, '(a,i8)' ) '  Number of data lines needed, N =       ', n
      stop
    end if

    if ( line(1:1) == '#' .or. len_trim ( line ) == 0 ) then
      cycle
    end if

    call s_to_i4vec ( line, m, x, ierror )

    if ( ierror /= 0 ) then
      cycle
    end if

    j = j + 1

    table(1:m,j) = x(1:m)

  end do

  close ( unit = input_unit )

  return
end
subroutine i4mat_header_read ( input_filename, m, n )

!*****************************************************************************80
!
!! I4MAT_HEADER_READ reads the header from an I4MAT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 June 2004
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
    write ( *, '(a)' ) 'I4MAT_HEADER_READ - Fatal error!'
    write ( *, '(a)' ) '  There was some kind of I/O problem while trying'
    write ( *, '(a)' ) '  to count the number of data columns in'
    write ( *, '(a)' ) '  the file "' // trim ( input_filename ) // '".'
    stop
  end if

  call file_row_count ( input_filename, n )

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4MAT_HEADER_READ - Fatal error!'
    write ( *, '(a)' ) '  There was some kind of I/O problem while trying'
    write ( *, '(a)' ) '  to count the number of data rows in'
    write ( *, '(a)' ) '  the file "' // trim ( input_filename ) // '".'
    stop
  end if

  return
end
subroutine i4mat_transpose_print ( m, n, a, title )

!*****************************************************************************80
!
!! I4MAT_TRANSPOSE_PRINT prints an I4MAT, transposed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, integer ( kind = 4 ) A(M,N), an M by N matrix to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  character ( len = * ) title

  call i4mat_transpose_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine i4mat_transpose_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! I4MAT_TRANSPOSE_PRINT_SOME prints some of the transpose of an I4MAT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, integer ( kind = 4 ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ( kind = 4 ) ILO, JLO, the first row and column to print.
!
!    Input, integer ( kind = 4 ) IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 10
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  character ( len = 7 ) ctemp(incx)
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
      write ( ctemp(i2), '(i7)') i
    end do

    write ( *, '(''  Row '',10a7)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Col'
    write ( *, '(a)' ) ' '

    j2lo = max ( jlo, 1 )
    j2hi = min ( jhi, n )

    do j = j2lo, j2hi

      do i2 = 1, inc

        i = i2lo - 1 + i2

        write ( ctemp(i2), '(i7)' ) a(i,j)

      end do

      write ( *, '(i5,1x,10a7)' ) j, ( ctemp(i), i = 1, inc )

    end do

  end do

  return
end
subroutine i4mat_write ( output_filename, m, n, table )

!*****************************************************************************80
!
!! I4MAT_WRITE writes an I4MAT file.
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
!    Input, integer ( kind = 4 ) TABLE(M,N), the table data.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) j
  character ( len = * ) output_filename
  integer ( kind = 4 ) output_status
  integer ( kind = 4 ) output_unit
  character ( len = 30 ) string
  integer ( kind = 4 ) table(m,n)
!
!  Open the file.
!
  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_filename, &
    status = 'replace', iostat = output_status )

  if ( output_status /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4MAT_WRITE - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the output file "' // &
      trim ( output_filename ) // '" on unit ', output_unit
    output_unit = -1
    stop
  end if
!
!  Create a format string.
!
  write ( string, '(a1,i8,a4)' ) '(', m, 'i10)'
!
!  Write the data.
!
  do j = 1, n
    write ( output_unit, string ) table(1:m,j)
  end do
!
!  Close the file.
!
  close ( unit = output_unit )

  return
end
subroutine mesh_base_one ( node_num, element_order, element_num, element_node )

!*****************************************************************************80
!
!! MESH_BASE_ONE ensures that the element definition is one-based.
!
!  Discussion:
!
!    The ELEMENT_NODE array contains nodes indices that form elements.
!    The convention for node indexing might start at 0 or at 1.
!    Since a FORTRAN90 program will naturally assume a 1-based indexing, it is
!    necessary to check a given element definition and, if it is actually
!    0-based, to convert it.
!
!    This function attempts to detect 9-based node indexing and correct it.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 October 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, int NODE_NUM, the number of nodes.
!
!    Input, int ELEMENT_ORDER, the order of the elements.
!
!    Input, int ELEMENT_NUM, the number of elements.
!
!    Input/output, int ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM), the element
!    definitions.
!
  implicit none

  integer ( kind = 4 ) element_num
  integer ( kind = 4 ) element_order

  integer ( kind = 4 ) element
  integer ( kind = 4 ) element_node(element_order,element_num)
  integer ( kind = 4 ) node
  integer ( kind = 4 ) node_max
  integer ( kind = 4 ) node_min
  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) order

  node_min = node_num + 1
  node_max = -1

  node_min = minval ( element_node(1:element_order,1:element_num) )
  node_max = maxval ( element_node(1:element_order,1:element_num) )

  if ( node_min == 0 .and. node_max == node_num - 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' )'MESH_BASE_ONE:'
    write ( *, '(a)' )'  The element indexing appears to be 0-based!'
    write ( *, '(a)' )'  This will be converted to 1-based.'
    element_node(1:element_order,1:element_num) = &
      element_node(1:element_order,1:element_num) + 1
  else if ( node_min == 1 .and. node_max == node_num  ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' )'MESH_BASE_ONE:'
    write ( *, '(a)' )'  The element indexing appears to be 1-based!'
    write ( *, '(a)' )'  No conversion is necessary.'
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MESH_BASE_ONE - Warning!'
    write ( *, '(a)' ) '  The element indexing is not of a recognized type.'
    write ( *, '(a,i8)' ) '  NODE_MIN = ', node_min
    write ( *, '(a,i8)' ) '  NODE_MAX = ', node_max
    write ( *, '(a,i8)' ) '  NODE_NUM = ', node_num
  end if

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
      write ( ctemp(i2), '(i7,7x)') i
    end do

    write ( *, '(''  Row   '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Col'

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
  write ( string, '(a1,i8,a1,i8,a1,i8,a1)' ) '(', m, 'g', 14, '.', 6, ')'
!
!  Write the data.
!
  do j = 1, n
    write ( output_unit, string ) table(1:m,j)
  end do
!
!  Close the file.
!
  close ( unit = output_unit )

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
!    Output, integer ( kind = 4 ) LENGTH, the number of characters of S used to make IVAL.
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
subroutine s_to_i4vec ( s, n, ivec, ierror )

!*****************************************************************************80
!
!! S_TO_I4VEC reads an I4VEC from a string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 October 2003
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
!    Output, integer ( kind = 4 ) IVEC(N), the values read from the string.
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
  integer ( kind = 4 ) ivec(n)
  integer ( kind = 4 ) length
  character ( len = * ) s

  i = 0

  ilo = 1

  do while ( i < n )

    i = i + 1

    call s_to_i4 ( s(ilo:), ivec(i), ierror, length )

    if ( ierror /= 0 ) then
      ierror = -i
      exit
    end if

    ilo = ilo + length

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

  logical ch_eqi
  character c
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
  character ( len = * ) s

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
subroutine sort_heap_external ( n, indx, i, j, isgn )

!*****************************************************************************80
!
!! SORT_HEAP_EXTERNAL externally sorts a list of items into ascending order.
!
!  Discussion:
!
!    The actual list of data is not passed to the routine.  Hence this
!    routine may be used to sort integers, real ( kind = 8 )s, numbers, names,
!    dates, shoe sizes, and so on.  After each call, the routine asks
!    the user to compare or interchange two items, until a special
!    return value signals that the sorting is completed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 February 2004
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Albert Nijenhuis and Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items to be sorted.
!
!    Input/output, integer ( kind = 4 ) INDX, the main communication signal.
!
!    The user must set INDX to 0 before the first call.
!    Thereafter, the user should not change the value of INDX until
!    the sorting is done.
!
!    On return, if INDX is
!
!      greater than 0,
!      * interchange items I and J;
!      * call again.
!
!      less than 0,
!      * compare items I and J;
!      * set ISGN = -1 if I < J, ISGN = +1 if J < I;
!      * call again.
!
!      equal to 0, the sorting is done.
!
!    Output, integer ( kind = 4 ) I, J, the indices of two items.
!    On return with INDX positive, elements I and J should be interchanged.
!    On return with INDX negative, elements I and J should be compared, and
!    the result reported in ISGN on the next call.
!
!    Input, integer ( kind = 4 ) ISGN, results of comparison of elements I and J.
!    (Used only when the previous call returned INDX less than 0).
!    ISGN <= 0 means I is less than or equal to J;
!    0 <= ISGN means I is greater than or equal to J.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ), save :: i_save = 0
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  integer ( kind = 4 ), save :: j_save = 0
  integer ( kind = 4 ), save :: k = 0
  integer ( kind = 4 ), save :: k1 = 0
  integer ( kind = 4 ) n
  integer ( kind = 4 ), save :: n1 = 0
!
!  INDX = 0: This is the first call.
!
  if ( indx == 0 ) then

    i_save = 0
    j_save = 0
    k = n / 2
    k1 = k
    n1 = n
!
!  INDX < 0: The user is returning the results of a comparison.
!
  else if ( indx < 0 ) then

    if ( indx == -2 ) then

      if ( isgn < 0 ) then
        i_save = i_save + 1
      end if

      j_save = k1
      k1 = i_save
      indx = -1
      i = i_save
      j = j_save
      return

    end if

    if ( 0 < isgn ) then
      indx = 2
      i = i_save
      j = j_save
      return
    end if

    if ( k <= 1 ) then

      if ( n1 == 1 ) then
        i_save = 0
        j_save = 0
        indx = 0
      else
        i_save = n1
        n1 = n1 - 1
        j_save = 1
        indx = 1
      end if

      i = i_save
      j = j_save
      return

    end if

    k = k - 1
    k1 = k
!
!  0 < INDX, the user was asked to make an interchange.
!
  else if ( indx == 1 ) then

    k1 = k

  end if

  do

    i_save = 2 * k1

    if ( i_save == n1 ) then
      j_save = k1
      k1 = i_save
      indx = -1
      i = i_save
      j = j_save
      return
    else if ( i_save <= n1 ) then
      j_save = i_save + 1
      indx = -2
      i = i_save
      j = j_save
      return
    end if

    if ( k <= 1 ) then
      exit
    end if

    k = k - 1
    k1 = k

  end do

  if ( n1 == 1 ) then
    i_save = 0
    j_save = 0
    indx = 0
    i = i_save
    j = j_save
  else
    i_save = n1
    n1 = n1 - 1
    j_save = 1
    indx = 1
    i = i_save
    j = j_save
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
subroutine triangulation_order3_refine_compute ( node_num1, triangle_num1, &
  node_xy1, triangle_node1, node_num2, triangle_num2, edge_data, node_xy2, &
  triangle_node2 )

!*****************************************************************************80
!
!! TRIANGULATION_ORDER3_REFINE_COMPUTE computes a refined order 3 triangulation.
!
!  Discussion:
!
!    Given a triangle defined by nodes 1, 2, 3, we need to generate
!    nodes 12, 23, and 13, and create 4 new subtriangles, T1, T2, T3
!    and T4.
!
!    The task is more complicated by the fact that we are working with
!    a mesh of triangles, so that we want to create a node only once,
!    even though it may be shared by other triangles.
!
!          3
!         / \
!        /T3 \
!      13----23
!      / \T4 / \
!     /T1 \ /T2 \
!    1----12-----2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM1, the number of nodes.
!
!    Input, integer ( kind = 4 ) TRIANGLE_NUM1, the number of triangles.
!
!    Input, real ( kind = 8 ) NODE_XY1(2,NODE_NUM1), the nodes.
!
!    Input, integer ( kind = 4 ) TRIANGLE_NODE1(3,TRIANGLE_NUM1), the nodes
!    that make up the triangles.  These should be listed in counterclockwise
!    order.
!
!    Input, integer ( kind = 4 ) NODE_NUM2, the number of nodes in the refined
!    mesh.
!
!    Input, integer ( kind = 4 ) TRIANGLE_NUM2, the number of triangles in the
!    refined mesh.
!
!    Input, integer ( kind = 4 ) EDGE_DATA(5,3*TRIANGLE_NUM1), edge information
!    computed by TRIANGULATION_ORDER3_REFINE_SIZE.
!
!    Output, real ( kind = 8 ) NODE_XY2(2,NODE_NUM2), the refined nodes.
!
!    Output, integer ( kind = 4 ) TRIANGLE_NODE2(3,TRIANGLE_NUM2), the nodes
!    that make up the triangles in the refined mesh.
!
  implicit none

  integer ( kind = 4 ) node_num1
  integer ( kind = 4 ) node_num2
  integer ( kind = 4 ) triangle_num1
  integer ( kind = 4 ) triangle_num2

  integer ( kind = 4 ) edge
  integer ( kind = 4 ) edge_data(5,3*triangle_num1)
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n1_old
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n2_old
  integer ( kind = 4 ) node
  real  ( kind = 8 ) node_xy1(2,node_num1)
  real  ( kind = 8 ) node_xy2(2,node_num2)
  integer ( kind = 4 ) triangle_node1(3,triangle_num1)
  integer ( kind = 4 ) triangle_node2(3,triangle_num2)
  integer ( kind = 4 ) triangle1
  integer ( kind = 4 ) v1
  integer ( kind = 4 ) v2
!
!  Copy the old nodes.
!
  node_xy2(1:2,1:node_num1) = node_xy1(1:2,1:node_num1)

  triangle_node2(1:3,1:triangle_num2) = -1
!
!  We can assign the existing nodes to the new triangles.
!
  do triangle1 = 1, triangle_num1
    triangle_node2(1,(triangle1-1)*4+1) = triangle_node1(1,triangle1)
    triangle_node2(2,(triangle1-1)*4+2) = triangle_node1(2,triangle1)
    triangle_node2(3,(triangle1-1)*4+3) = triangle_node1(3,triangle1)
  end do

  node = node_num1

  n1_old = -1
  n2_old = -1

  do edge = 1, 3 * triangle_num1

    n1 = edge_data(1,edge)
    n2 = edge_data(2,edge)
!
!  If this edge is new, create the coordinates and index for this node.
!
    if ( n1 /= n1_old .or. n2 /= n2_old ) then
      node = node + 1

      if ( node_num2 < node ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TRIANGLE_MESH_ORDER3_REFINE - Fatal error!'
        write ( *, '(a)' ) '  Node index exceeds NODE_NUM2.'
        stop
      end if

      node_xy2(1:2,node) = &
        ( node_xy2(1:2,n1) + node_xy2(1:2,n2) ) / 2.0D+00

      n1_old = n1
      n2_old = n2

    end if
!
!  Assign the node to triangles.
!
    v1 = edge_data(3,edge)
    v2 = edge_data(4,edge)
    triangle1 = edge_data(5,edge)

    if ( v1 == 1 .and. v2 == 2 ) then

      triangle_node2(1,(triangle1-1)*4+2) = node
      triangle_node2(2,(triangle1-1)*4+1) = node
      triangle_node2(3,(triangle1-1)*4+4) = node

    else if ( v1 == 1 .and. v2 == 3 ) then

      triangle_node2(1,(triangle1-1)*4+3) = node
      triangle_node2(2,(triangle1-1)*4+4) = node
      triangle_node2(3,(triangle1-1)*4+1) = node

    else if ( v1 == 2 .and. v2 == 3 ) then

      triangle_node2(1,(triangle1-1)*4+4) = node
      triangle_node2(2,(triangle1-1)*4+3) = node
      triangle_node2(3,(triangle1-1)*4+2) = node

    end if

  end do

  return
end
subroutine triangulation_order3_refine_size ( node_num1, triangle_num1, &
  triangle_node1, node_num2, triangle_num2, edge_data )

!*****************************************************************************80
!
!! TRIANGULATION_ORDER3_REFINE_SIZE sizes a refined order 3 triangulation.
!
!  Discussion:
!
!    Given a triangle defined by nodes 1, 2, 3, we need to generate
!    nodes 12, 23, and 13, and create 4 new subtriangles, T1, T2, T3
!    and T4.
!
!    The task is more complicated by the fact that we are working with
!    a mesh of triangles, so that we want to create a node only once,
!    even though it may be shared by other triangles.
!
!          3
!         / \
!        /T3 \
!      13----23
!      / \T4 / \
!     /T1 \ /T2 \
!    1----12-----2
!
!    This routine simply determines the sizes of the resulting node
!    and triangle arrays.
!
!    The primary amount of work occurs in sorting a list of 3 * TRIANGLE_NUM
!    data items, one item for every edge of every triangle.  Each
!    data item records, for a given edge, the global indices
!    of the two endpoints, the local indices of the two endpoints,
!    and the index of the triangle.
!
!    Through careful sorting, it is possible to arrange this data in
!    a way that allows the proper generation of the interpolated nodes.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM1, the number of nodes in the original mesh.
!
!    Input, integer ( kind = 4 ) TRIANGLE_NUM1, the number of triangles in the
!    original mesh.
!
!    Input, integer ( kind = 4 ) TRIANGLE_NODE1(3,TRIANGLE_NUM1), the indices of the nodes
!    that form the triangles in the input mesh.
!
!    Output, integer ( kind = 4 ) NODE_NUM2, the number of nodes in the refined mesh.
!
!    Output, integer ( kind = 4 ) TRIANGLE_NUM2, the number of triangles in the
!    refined mesh.
!
!    Output, integer ( kind = 4 ) EDGE_DATA(5,3*TRIANGLE_NUM1), edge data that will
!    be needed by TRIANGULATION_ORDER3_REFINE_COMPUTE.
!
  implicit none

  integer ( kind = 4 ) node_num1
  integer ( kind = 4 ) triangle_num1

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) edge
  integer ( kind = 4 ) edge_data(5,3*triangle_num1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n1_old
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n2_old
  integer ( kind = 4 ) node_num2
  integer ( kind = 4 ) triangle
  integer ( kind = 4 ) triangle_node1(3,triangle_num1)
  integer ( kind = 4 ) triangle_num2
!
!  Step 1.
!  From the list of nodes for triangle T, of the form: (I,J,K)
!  construct the edge relations:
!
!    (I,J,1,2,T)
!    (I,K,1,3,T)
!    (J,K,2,3,T)
!
!  In order to make matching easier, we reorder each pair of nodes
!  into ascending order.
!
  do triangle = 1, triangle_num1

    i = triangle_node1(1,triangle)
    j = triangle_node1(2,triangle)
    k = triangle_node1(3,triangle)

    a = min ( i, j )
    b = max ( i, j )

    edge_data(1:5,3*(triangle-1)+1) = (/ a, b, 1, 2, triangle /)

    a = min ( i, k )
    b = max ( i, k )

    edge_data(1:5,3*(triangle-1)+2) = (/ a, b, 1, 3, triangle /)

    a = min ( j, k )
    b = max ( j, k )

    edge_data(1:5,3*(triangle-1)+3) = (/ a, b, 2, 3, triangle /)

  end do
!
!  Step 2. Perform an ascending dictionary sort on the neighbor relations.
!  We only intend to sort on rows 1:2; the routine we call here
!  sorts on the full column but that won't hurt us.
!
!  What we need is to find all cases where triangles share an edge.
!  By sorting the columns of the EDGE_DATA array, we will put shared edges
!  next to each other.
!
  call i4col_sort_a ( 5, 3*triangle_num1, edge_data )
!
!  Step 3. All the triangles which share an edge show up as consecutive
!  columns with identical first two entries.  Figure out how many new
!  nodes there are, and allocate space for their coordinates.
!
  node_num2 = node_num1

  n1_old = -1
  n2_old = -1

  do edge = 1, 3 * triangle_num1
    n1 = edge_data(1,edge)
    n2 = edge_data(2,edge)
    if ( n1 /= n1_old .or. n2 /= n2_old ) then
      node_num2 = node_num2 + 1
      n1_old = n1
      n2_old = n2
    end if
  end do

  triangle_num2 = 4 * triangle_num1

  return
end
subroutine triangulation_order6_refine_compute ( node_num1, triangle_num1, &
  node_xy1, triangle_node1, node_num2, triangle_num2, edge_data, &
  node_xy2, triangle_node2 )

!*****************************************************************************80
!
!! TRIANGULATION_ORDER6_REFINE_COMPUTE computes a refined order 6 triangulation.
!
!  Discussion:
!
!    Given a quadratic triangle defined by nodes 1, 2, 3, 4, 5, 6, we
!    need to generate nodes 14, 16, 24, 25, 35, 36, 45, 46, 56, and 4 new
!    quadratic subtriangles T1, T2, T3 and T4.
!
!    The task is more complicated by the fact that we are working with
!    a mesh of triangles, so that we want to create a node only once,
!    even though it may be shared by other triangles.  (In fact, only
!    the new nodes on the edges can be shared, and then only by at most
!    one other triangle.)
!
!            3
!           / \
!          36 35
!         / T3  \
!        6--56---5
!       / \ T4  / \
!      16 46  45  25
!     / T1  \ / T2  \
!    1--14---4--24---2
!
!    This routine is given sorted information defining the edges, and uses
!    it to build the new node and triangle arrays.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM1, the number of nodes.
!
!    Input, integer ( kind = 4 ) TRIANGLE_NUM1, the number of triangles.
!
!    Input, real ( kind = 8 ) NODE_XY1(2,NODE_NUM1), the nodes.
!
!    Input, integer ( kind = 4 ) TRIANGLE_NODE1(6,TRIANGLE_NUM1), the nodes that make up the
!    triangles.  These should be listed in counterclockwise order.
!
!    Input, integer ( kind = 4 ) NODE_NUM2, the number of nodes in the refined mesh.
!
!    Input, integer ( kind = 4 ) TRIANGLE_NUM2, the number of triangles in the refined mesh.
!
!    Input, integer ( kind = 4 ) EDGE_DATA(5,3*TRIANGLE_NUM1), edge data.
!
!    Output, real ( kind = 8 ) NODE_XY2(2,NODE_NUM2), the refined nodes.
!
!    Output, integer ( kind = 4 ) TRIANGLE_NODE2(6,TRIANGLE_NUM2), the nodes 
!    that make up the triangles in the refined mesh.
!
  implicit none

  integer ( kind = 4 ) node_num1
  integer ( kind = 4 ) node_num2
  integer ( kind = 4 ) triangle_num1
  integer ( kind = 4 ) triangle_num2

  integer ( kind = 4 ) edge
  integer ( kind = 4 ) edge_data(5,3*triangle_num1)
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) l2
  integer ( kind = 4 ) l3
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n1_old
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n2_old
  integer ( kind = 4 ) node
  integer ( kind = 4 ) node_num
  real ( kind = 8 ) node_xy1(2,node_num1)
  real ( kind = 8 ) node_xy2(2,node_num2)
  integer ( kind = 4 ) t1
  integer ( kind = 4 ) t2
  integer ( kind = 4 ) t3
  integer ( kind = 4 ) t4
  integer ( kind = 4 ) triangle_node1(6,triangle_num1)
  integer ( kind = 4 ) triangle_node2(6,triangle_num2)
  integer ( kind = 4 ) triangle1
  integer ( kind = 4 ) v1
  integer ( kind = 4 ) v2
  integer ( kind = 4 ) v3
  integer ( kind = 4 ) v4
  integer ( kind = 4 ) v5
  integer ( kind = 4 ) v6
!
!  Step 1:
!  Copy old nodes.
!
  node_xy2(1:2,1:node_num1) = node_xy1(1:2,1:node_num1)
!
!  Copy indices of existing nodes into new triangle array.
!
  triangle_node2(1:6,1:triangle_num2) = -1

  do triangle1 = 1, triangle_num1

    t1 = ( triangle1 - 1 ) * 4 + 1
    t2 = ( triangle1 - 1 ) * 4 + 2
    t3 = ( triangle1 - 1 ) * 4 + 3
    t4 = ( triangle1 - 1 ) * 4 + 4

    triangle_node2(1,t1) = triangle_node1(1,triangle1)
    triangle_node2(2,t1) = triangle_node1(4,triangle1)
    triangle_node2(3,t1) = triangle_node1(6,triangle1)

    triangle_node2(1,t2) = triangle_node1(4,triangle1)
    triangle_node2(2,t2) = triangle_node1(2,triangle1)
    triangle_node2(3,t2) = triangle_node1(5,triangle1)

    triangle_node2(1,t3) = triangle_node1(6,triangle1)
    triangle_node2(2,t3) = triangle_node1(5,triangle1)
    triangle_node2(3,t3) = triangle_node1(3,triangle1)

    triangle_node2(1,t4) = triangle_node1(5,triangle1)
    triangle_node2(2,t4) = triangle_node1(6,triangle1)
    triangle_node2(3,t4) = triangle_node1(4,triangle1)

  end do
!
!  Step 2.
!  Examine sorted edge information.  The first time an edge is encountered,
!  generate two new nodes, then assign them (usually) to the four subtriangles
!  of the two triangles that share that edge.
!
  node = node_num1

  n1_old = -1
  n2_old = -1

  do edge = 1, 3 * triangle_num1

    n1 = edge_data(1,edge)
    n2 = edge_data(2,edge)

    l1 = edge_data(3,edge)
    l3 = edge_data(4,edge)

    if ( l1 == 1 .and. l3 == 2 ) then
      l2 = 4
    else if ( l1 == 1 .and. l3 == 3 ) then
      l2 = 6
    else if ( l1 == 2 .and. l3 == 3 ) then
      l2 = 5
    end if

    triangle1 = edge_data(5,edge)
!
!  If this is the first time we've encountered this edge,
!  create the new new nodes.
!
    if ( n1 /= n1_old .or. n2 /= n2_old ) then

      n1_old = n1
      n2_old = n2

      v1 = triangle_node1(l1,triangle1)
      v2 = triangle_node1(l2,triangle1)
      v3 = triangle_node1(l3,triangle1)

      node = node + 1
      v4 = node
      node_xy2(1:2,node) = 0.5D+00 * ( node_xy1(1:2,v1) + node_xy1(1:2,v2) )

      node = node + 1
      v5 = node
      node_xy2(1:2,node) = 0.5D+00 * ( node_xy1(1:2,v2) + node_xy1(1:2,v3) )

    end if

    t1 = ( triangle1 - 1 ) * 4 + 1
    t2 = ( triangle1 - 1 ) * 4 + 2
    t3 = ( triangle1 - 1 ) * 4 + 3

    if ( l1 == 1 .and. l3 == 2 ) then

      if ( triangle_node1(1,triangle1) == v1 ) then
        triangle_node2(4,t1) = v4
        triangle_node2(4,t2) = v5
      else
        triangle_node2(4,t1) = v5
        triangle_node2(4,t2) = v4
      end if

    else if ( l1 == 1 .and. l3 == 3 ) then

      if ( triangle_node1(l1,triangle1) == v1 ) then
        triangle_node2(6,t1) = v4
        triangle_node2(6,t3) = v5
      else
        triangle_node2(6,t1) = v5
        triangle_node2(6,t3) = v4
      end if

    else if ( l1 == 2 .and. l3 == 3 ) then

      if ( triangle_node1(l1,triangle1) == v1 ) then
        triangle_node2(5,t3) = v4
        triangle_node2(5,t2) = v5
      else
        triangle_node2(5,t3) = v5
        triangle_node2(5,t2) = v4
      end if

    end if

  end do
!
!  Step 3.
!  Each old triangle has a single central subtriangle, for which we now
!  need to generate three new "interior" nodes.
!
  do triangle1 = 1, triangle_num1

    v4 = triangle_node1(4,triangle1)
    v5 = triangle_node1(5,triangle1)
    v6 = triangle_node1(6,triangle1)

    t1 = ( triangle1 - 1 ) * 4 + 1
    t2 = ( triangle1 - 1 ) * 4 + 2
    t3 = ( triangle1 - 1 ) * 4 + 3
    t4 = ( triangle1 - 1 ) * 4 + 4

    node = node + 1
    node_xy2(1:2,node) = 0.5D+00 * ( node_xy1(1:2,v5) + node_xy1(1:2,v6) )
    triangle_node2(4,t4) = node
    triangle_node2(4,t3) = node

    node = node + 1
    node_xy2(1:2,node) = 0.5D+00 * ( node_xy1(1:2,v6) + node_xy1(1:2,v4) )
    triangle_node2(5,t4) = node
    triangle_node2(5,t1) = node

    node = node + 1
    node_xy2(1:2,node) = 0.5D+00 * ( node_xy1(1:2,v4) + node_xy1(1:2,v5) )
    triangle_node2(6,t4) = node
    triangle_node2(6,t2) = node

  end do

  return
end
subroutine triangulation_order6_refine_size ( node_num1, triangle_num1, &
  triangle_node1, node_num2, triangle_num2, edge_data )

!*****************************************************************************80
!
!! TRIANGULATION_ORDER6_REFINE_SIZE sizes a refined order 6 triangulation.
!
!  Discussion:
!
!    Given a quadratic triangle defined by nodes 1, 2, 3, 4, 5, 6, we
!    need to generate nodes 14, 16, 24, 25, 35, 36, 45, 46, 56, and 4 new
!    quadratic subtriangles T1, T2, T3 and T4.
!
!    The task is more complicated by the fact that we are working with
!    a mesh of triangles, so that we want to create a node only once,
!    even though it may be shared by other triangles.  (In fact, only
!    the new nodes on the edges can be shared, and then only by at most
!    one other triangle.)
!
!            3
!           / \
!          36 35
!         / T3  \
!        6--56---5
!       / \ T4  / \
!      16 46  45  25
!     / T1  \ / T2  \
!    1--14---4--24---2
!
!    This routine determines the sizes of the resulting node and
!    triangles, and constructs an edge array that can be used to
!    properly number the new nodes.
!
!    The primary work occurs in sorting a list related to the edges.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM1, the number of nodes.
!
!    Input, integer ( kind = 4 ) TRIANGLE_NUM1, the number of triangles.
!
!    Input, integer ( kind = 4 ) TRIANGLE_NODE1(6,TRIANGLE_NUM1), the nodes that make up the
!    triangles.  These should be listed in counterclockwise order.
!
!    Input, integer ( kind = 4 ) NODE_NUM2, the number of nodes in the refined mesh.
!
!    Input, integer ( kind = 4 ) TRIANGLE_NUM2, the number of triangles 
!    in the refined mesh.
!
!    Output, integer ( kind = 4 ) EDGE_DATA(5,3*TRIANGLE_NUM1), edge data needed by
!    TRIANGULATION_ORDER6_REFINE_COMPUTE.
!
  implicit none

  integer ( kind = 4 ) node_num1
  integer ( kind = 4 ) triangle_num1

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) edge
  integer ( kind = 4 ) edge_data(5,3*triangle_num1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n1_old
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n2_old
  integer ( kind = 4 ) node_num2
  integer ( kind = 4 ) triangle_num2
  integer ( kind = 4 ) triangle_node1(6,triangle_num1)
  integer ( kind = 4 ) triangle1
!
!  Step 1:
!  From the list of vertices for triangle T, of the form: (I,J,K),
!  construct the edge relations:
!
!    (I,J,1,2,T)
!    (I,K,1,3,T)
!    (J,K,2,3,T)
!
!  To make matching easier, we reorder each pair of nodes into
!  ascending order.
!
  do triangle1 = 1, triangle_num1

    i = triangle_node1(1,triangle1)
    j = triangle_node1(2,triangle1)
    k = triangle_node1(3,triangle1)

    a = min ( i, j )
    b = max ( i, j )

    edge_data(1:5,3*(triangle1-1)+1) = (/ a, b, 1, 2, triangle1 /)

    a = min ( i, k )
    b = max ( i, k )

    edge_data(1:5,3*(triangle1-1)+2) = (/ a, b, 1, 3, triangle1 /)

    a = min ( j, k )
    b = max ( j, k )

    edge_data(1:5,3*(triangle1-1)+3) = (/ a, b, 2, 3, triangle1 /)

  end do
!
!  Step 2: Perform an ascending dictionary sort on the relations.
!
  call i4col_sort_a ( 5, 3*triangle_num1, edge_data )
!
!  Step 3: Each shared edge will show up twice, consecutively,
!  in the EDGE_DATA array.  Each unique edge will generate
!  two new nodes, and each triangle will generate three new nodes.
!
  node_num2 = node_num1

  n1_old = -1
  n2_old = -1

  do edge = 1, 3 * triangle_num1
    n1 = edge_data(1,edge)
    n2 = edge_data(2,edge)
    if ( n1 /= n1_old .or. n2 /= n2_old ) then
      node_num2 = node_num2 + 2
      n1_old = n1
      n2_old = n2
    end if
  end do

  node_num2 = node_num2 + 3 * triangle_num1

  triangle_num2 = 4 * triangle_num1

  return
end
