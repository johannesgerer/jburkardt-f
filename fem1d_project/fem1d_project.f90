program main

!*****************************************************************************80
!
!! MAIN is the main program for FEM1D_PROJECT.
!
!  Discussion:
!
!    FEM1D_PROJECT reads files defining a sampling of a (scalar or vector)
!    function of 1 argument, and a list of nodes and elements to use for
!    a finite element representation of the data.
!
!    It computes a set of finite element coefficients to be associated with
!    the given finite element mesh, and writes that information to a file
!    so that an FEM representation is formed by the node, element and value
!    files.
!
!  Usage:
!
!    fem1d_project sample_prefix fem_prefix
!
!    where 'sample_prefix' is the common prefix for the SAMPLE files:
!
!    * sample_prefix_nodes.txt,  the node coordinates where samples were taken,
!    * sample_prefix_values.txt, the sample values.
!
!    and 'fem_prefix' is the common prefix for the FEM files:
!
!    * fem_prefix_nodes.txt,    the node coordinates.
!    * fem_prefix_elements.txt, the nodes that make up each element;
!    * fem_prefix_values.txt,   the values defined at each node.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 June 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 255 ) fem_element_filename
  integer   ( kind = 4 ), allocatable :: fem_element_node(:,:)
  integer   ( kind = 4 )  fem_element_num
  integer   ( kind = 4 )  fem_element_order
  character ( len = 255 ) fem_prefix
  integer   ( kind = 4 )  fem_node_dim
  character ( len = 255 ) fem_node_filename
  integer   ( kind = 4 )  fem_node_num
  real ( kind = 8 ), allocatable :: fem_node_x(:,:)
  real ( kind = 8 ), allocatable :: fem_value(:,:)
  integer   ( kind = 4 )  fem_value_dim
  integer   ( kind = 4 )  fem_value_num
  character ( len = 255 ) fem_value_filename
  integer   ( kind = 4 )  iarg
  integer   ( kind = 4 )  iargc
  integer   ( kind = 4 )  ios
  integer   ( kind = 4 )  num_arg
  character ( len = 255 ) sample_prefix
  integer   ( kind = 4 )  sample_node_dim
  character ( len = 255 ) sample_node_filename
  integer   ( kind = 4 )  sample_node_num
  real ( kind = 8 ), allocatable :: sample_node_x(:,:)
  integer   ( kind = 4 )  sample_value_dim
  integer   ( kind = 4 )  sample_value_num
  real ( kind = 8 ), allocatable :: sample_value(:,:)
  character ( len = 255 ) sample_value_filename

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FEM1D_PROJECT'
  write ( *, '(a)' ) '  FORTRAN90 version.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Read files defining a sampling of a function of 1 argument.'
  write ( *, '(a)' ) '  Read files defining a finite element mesh.'
  write ( *, '(a)' ) '  Project the sample data onto the mesh, and'
  write ( *, '(a)' ) '  write a file of FEM coefficient values.'
!
!  Get the number of command line arguments.
!
  num_arg = iargc ( )

  if ( num_arg < 1 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Enter the sample file prefix:'
    read ( *, '(a)', iostat = ios ) sample_prefix

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FEM1D_PROJECT - Fatal error!'
      write ( *, '(a)' ) '  Unexpected read error!'
      stop
    end if

  else

    iarg = 1

    call getarg ( iarg, sample_prefix )

  end if

  if ( num_arg < 2 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Enter the FEM file prefix:'
    read ( *, '(a)', iostat = ios ) fem_prefix

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FEM1D_PROJECT - Fatal error!'
      write ( *, '(a)' ) '  Unexpected read error!'
      stop
    end if

  else

    iarg = 2

    call getarg ( iarg, fem_prefix )

  end if
!
!  Create the filenames.
!
  sample_node_filename = trim ( sample_prefix ) // '_nodes.txt'
  sample_value_filename = trim ( sample_prefix ) // '_values.txt'

  fem_node_filename = trim ( fem_prefix ) // '_nodes.txt'
  fem_element_filename = trim ( fem_prefix ) // '_elements.txt'
  fem_value_filename = trim ( fem_prefix ) // '_values.txt'
!
!  Read the SAMPLE data.
!
  call dtable_header_read ( sample_node_filename, sample_node_dim, &
    sample_node_num )

  allocate ( sample_node_x(1:sample_node_dim,1:sample_node_num) )

  call dtable_data_read ( sample_node_filename, sample_node_dim, &
    sample_node_num, sample_node_x )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Sample node spatial dimension is ', sample_node_dim
  write ( *, '(a,i8)' ) '  Sample node number is            ', sample_node_num

  if ( sample_node_dim /= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FEM1D_PROJECT - Fatal error!'
    write ( *, '(a)' ) '  Spatial dimension of the sample nodes is not 1.'
    stop
  end if

  call dtable_header_read ( sample_value_filename, sample_value_dim, &
    sample_value_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Sample value dimension is        ', sample_value_dim
  write ( *, '(a,i8)' ) '  Sample value number is           ', sample_value_num

  if ( sample_value_num /= sample_node_num ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FEM1D_PROJECT - Fatal error!'
    write ( *, '(a)' ) '  Number of sample nodes and values are not equal.'
    stop
  end if

  allocate ( sample_value(1:sample_value_dim,1:sample_value_num) )

  call dtable_data_read ( sample_value_filename, sample_value_dim, &
    sample_value_num, sample_value )
!
!  Read the FEM data.
!
  call dtable_header_read ( fem_node_filename, fem_node_dim, fem_node_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The FEM node dimension is        ', fem_node_dim
  write ( *, '(a,i8)' ) '  The FEM node number is           ', fem_node_num

  if ( fem_node_dim /= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FEM1D_PROJECT - Fatal error!'
    write ( *, '(a)' ) '  Spatial dimension of the nodes is not 1.'
    stop
  end if

  allocate ( fem_node_x(1:fem_node_dim,1:fem_node_num) )

  call dtable_data_read ( fem_node_filename, fem_node_dim, fem_node_num, fem_node_x )

  call itable_header_read ( fem_element_filename, fem_element_order, fem_element_num )

  write ( *, '(a,i8)' ) '  The FEM element order is         ', fem_element_order
  write ( *, '(a,i8)' ) '  The FEM element number is        ', fem_element_num

  allocate ( fem_element_node(1:fem_element_order,1:fem_element_num) )

  call itable_data_read ( fem_element_filename, fem_element_order, fem_element_num, &
    fem_element_node )
!
!  Compute the FEM values.
!
  fem_value_dim = sample_value_dim
  fem_value_num = fem_node_num
  allocate ( fem_value(1:fem_value_dim,1:fem_value_num) )

  call fem1d_approximate ( sample_node_num, sample_value_dim, sample_node_x, &
    sample_value, fem_node_num, fem_node_x, fem_element_order, &
    fem_element_num, fem_value_dim, fem_value_num, fem_value )
!
!  Write the FEM values.
!
  call dtable_write0 ( fem_value_filename, fem_value_dim, &
    fem_value_num, fem_value )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  FEM value data written to "' &
    // trim ( fem_value_filename ) // '".'
!
!  Free memory.
!
  deallocate ( fem_element_node )
  deallocate ( fem_node_x )
  deallocate ( fem_value )
  deallocate ( sample_node_x )
  deallocate ( sample_value )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FEM1D_PROJECT'
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

  character              c
  integer   ( kind = 4 ) itemp

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
  integer   ( kind = 4 ) digit

  if ( lge ( c, '0' ) .and. lle ( c, '9' ) ) then

    digit = ichar ( c ) - 48

  else if ( c == ' ' ) then

    digit = 0

  else

    digit = -1

  end if

  return
end
subroutine dtable_data_read ( input_file_name, m, n, table )

!*****************************************************************************80
!
!! DTABLE_DATA_READ reads data from a DTABLE file.
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

  integer   ( kind = 4 )  m
  integer   ( kind = 4 )  n

  integer   ( kind = 4 )  ierror
  character ( len = * )   input_file_name
  integer   ( kind = 4 )  input_status
  integer   ( kind = 4 )  input_unit
  integer   ( kind = 4 )  j
  character ( len = 255 ) line
  real ( kind = 8 )  table(m,n)
  real ( kind = 8 )  x(m)

  ierror = 0

  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_file_name, status = 'old', &
    iostat = input_status )

  if ( input_status /= 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DTABLE_DATA_READ - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the input file "' // &
      trim ( input_file_name ) // '" on unit ', input_unit
    stop
  end if

  j = 0

  do while ( j < n )

    read ( input_unit, '(a)', iostat = input_status ) line

    if ( input_status /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DTABLE_DATA_READ - Fatal error!'
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
subroutine dtable_header_read ( input_file_name, m, n )

!*****************************************************************************80
!
!! DTABLE_HEADER_READ reads the header from a DTABLE file.
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
  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) n

  call file_column_count ( input_file_name, m )

  if ( m <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DTABLE_HEADER_READ - Fatal error!'
    write ( *, '(a)' ) '  There was some kind of I/O problem while trying'
    write ( *, '(a)' ) '  to count the number of data columns in'
    write ( *, '(a)' ) '  the file "' // trim ( input_file_name ) // '".'
    stop
  end if

  call file_row_count ( input_file_name, n )

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DTABLE_HEADER_READ - Fatal error!'
    write ( *, '(a)' ) '  There was some kind of I/O problem while trying'
    write ( *, '(a)' ) '  to count the number of data rows in'
    write ( *, '(a)' ) '  the file "' // trim ( input_file_name ) // '".'
    stop
  end if

  return
end
subroutine dtable_write0 ( output_file_name, m, n, table )

!*****************************************************************************80
!
!! DTABLE_WRITE0 writes a DTABLE file with no headers.
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
!    Input, character ( len = * ) OUTPUT_FILE_NAME, the output file name.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) TABLE(M,N), the table data.
!
  implicit none

  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) n

  integer   ( kind = 4 ) j
  character ( len = * )  output_file_name
  integer   ( kind = 4 ) output_status
  integer   ( kind = 4 ) output_unit
  character ( len = 30 ) string
  real ( kind = 8 ) table(m,n)
!
!  Open the file.
!
  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_file_name, &
    status = 'replace', iostat = output_status )

  if ( output_status /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DTABLE_WRITE0 - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the output file "' // &
      trim ( output_file_name ) // '" on unit ', output_unit
    output_unit = -1
    stop
  end if
!
!  Create a format string.
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
subroutine fem1d_approximate ( sample_node_num, sample_value_dim, &
  sample_node_x, sample_value, fem_node_num, fem_node_x, fem_element_order, &
  fem_element_num, fem_value_dim, fem_value_num, fem_value )

!*****************************************************************************80
!
!! FEM1D_APPROXIMATE approximates data at sample points with an FEM function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) SAMPLE_NODE_NUM, the number of sample points.
!
!    Input, integer ( kind = 4 ) SAMPLE_VALUE_DIM, the value dimension.
!
!    Input, real ( kind = 8 ) SAMPLE_NODE_X(SAMPLE_NODE_NUM), the sample nodes.
!
!    Input, real ( kind = 8 ) SAMPLE_VALUE(VALUE_DIM,SAMPLE_NODE_NUM),
!    the values at sample nodes.
!
!    Input, integer ( kind = 4 ) FEM_NODE_NUM, the number of FEM nodes.
!
!    Input, real ( kind = 8 ) FEM_NODE_X(FEM_NODE_NUM), the FEM nodes.
!
!    Input, integer ( kind = 4 ) FEM_ELEMENT_ORDER, the element order.
!
!    Input, integer ( kind = 4 ) FEM_ELEMENT_NUM, the number of elements.
!
!    Input, integer ( kind = 4 ) FEM_VALUE_DIM, the FEM value dimension.
!
!    Input, integer ( kind = 4 ) FEM_VALUE_NUM, the number of FEM values.
!
!    Output, real ( kind = 8 ) FEM_VALUE(FEM_VALUE_DIM,FEM_VALUE_NUM),
!    the FEM values.
!
  implicit none

  integer ( kind = 4 ) fem_node_num
  integer ( kind = 4 ) fem_value_dim
  integer ( kind = 4 ) fem_value_num
  integer ( kind = 4 ), parameter :: quad_num = 2
  integer ( kind = 4 ) sample_node_num
  integer ( kind = 4 ) sample_value_dim
  integer ( kind = 4 ) sample_value_num

  real ( kind = 8 ) a(3,fem_node_num)
  real ( kind = 8 ) a1
  real ( kind = 8 ) b(fem_node_num,fem_value_dim)
  real ( kind = 8 ) b1
  integer ( kind = 4 ) dim
  integer ( kind = 4 ) fem_element_num
  integer ( kind = 4 ) fem_element_order
  real ( kind = 8 ) fem_node_x(fem_node_num)
  real ( kind = 8 ) fem_value(fem_value_dim,fem_value_num)
  integer ( kind = 4 ) i
  real ( kind = 8 ) integral
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l
  integer ( kind = 4 ) phi_num
  real ( kind = 8 ) phi_v(3)
  real ( kind = 8 ) phi_x(3)
  real ( kind = 8 ) phil
  real ( kind = 8 ) phir
  real ( kind = 8 ) phis
  integer ( kind = 4 ) quad
  real ( kind = 8 ) :: quad_x(quad_num) = (/ &
    -0.577350269189625764509148780502D+00, &
     0.577350269189625764509148780502D+00 /)
  real ( kind = 8 ) :: quad_w(quad_num) = (/ 1.0D+00, 1.0D+00 /)
  integer ( kind = 4 ) r
  integer ( kind = 4 ) sample
  real ( kind = 8 ) sample_node_x(sample_node_num)
  real ( kind = 8 ) sample_value(sample_value_dim,sample_node_num)
  real ( kind = 8 ) wq
  real ( kind = 8 ) x(fem_value_num,fem_value_dim)
  real ( kind = 8 ) xl
  real ( kind = 8 ) xq
  real ( kind = 8 ) xr
!
!  Set up the matrix A.
!
  a(1:3,1:fem_node_num) = 0.0D+00

  do l = 1, fem_node_num - 1

    r = l + 1
    xl = fem_node_x(l)
    xr = fem_node_x(r)

    do quad = 1, quad_num

      xq = ( ( 1.0D+00 - quad_x(quad) ) * xl   &
           + ( 1.0D+00 + quad_x(quad) ) * xr ) &
           /   2.0D+00

      wq = quad_w(quad) * ( xr - xl ) / 2.0D+00

      phil = (      xq - xr ) &
           / ( xl      - xr )

      phir = ( xl - xq      ) &
           / ( xl      - xr )

      a(2,l) = a(2,l) + wq * phil * phil
      a(3,l) = a(3,l) + wq * phil * phir

      a(1,r) = a(1,r) + wq * phir * phil
      a(2,r) = a(2,r) + wq * phir * phir

    end do

  end do
!
!  Set up the right hand side b.
!
  b(1:fem_node_num,1:fem_value_dim) = 0.0D+00

  do i = 1, fem_node_num

    if ( i == 1 ) then
      phi_num = 2
      phi_x(1) = fem_node_x(1)
      phi_x(2) = fem_node_x(2)
      phi_v(1) = 1.0D+00
      phi_v(2) = 0.0D+00
    else if ( i < fem_node_num ) then
      phi_num = 3
      phi_x(1) = fem_node_x(i-1)
      phi_x(2) = fem_node_x(i)
      phi_x(3) = fem_node_x(i+1)
      phi_v(1) = 0.0D+00
      phi_v(2) = 1.0D+00
      phi_v(3) = 0.0D+00
    else if ( i == fem_node_num ) then
      phi_num = 2
      phi_x(1) = fem_node_x(fem_node_num-1)
      phi_x(2) = fem_node_x(fem_node_num)
      phi_v(1) = 0.0D+00
      phi_v(2) = 1.0D+00
    end if

    a1 = phi_x(1)
    b1 = phi_x(phi_num)

    do dim = 1, fem_value_dim

      call piecewise_linear_product_quad ( a1, b1, phi_num, phi_x, phi_v, &
        sample_node_num, sample_node_x, sample_value(dim,1:sample_node_num), &
        integral )

      b(i,dim) = integral

    end do

  end do
!
!  Solve A * X = B.
!
  call r83_np_fss ( fem_node_num, a, fem_value_dim, b, x )

  fem_value(1:fem_value_dim,1:fem_value_num) = transpose ( x )

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

  integer   ( kind = 4 )  column_num
  logical                 got_one
  character ( len = * )   input_file_name
  integer   ( kind = 4 )  input_status
  integer   ( kind = 4 )  input_unit
  character ( len = 255 ) line
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

  integer   ( kind = 4 )  bad_num
  integer   ( kind = 4 )  comment_num
  integer   ( kind = 4 )  ierror
  character ( len = * )   input_file_name
  integer   ( kind = 4 )  input_status
  integer   ( kind = 4 )  input_unit
  character ( len = 255 ) line
  integer   ( kind = 4 )  record_num
  integer   ( kind = 4 )  row_num

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
subroutine itable_data_read ( input_file_name, m, n, table )

!*****************************************************************************80
!
!! ITABLE_DATA_READ reads data from an ITABLE file.
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
!    Input, character ( len = * ) INPUT_FILE_NAME, the name of the input file.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Output, integer ( kind = 4 ) TABLE(M,N), the table data.
!
  implicit none

  integer   ( kind = 4 )   m
  integer   ( kind = 4 )   n

  integer   ( kind = 4 )   ierror
  character ( len = * )    input_file_name
  integer   ( kind = 4 )   input_status
  integer   ( kind = 4 )   input_unit
  integer   ( kind = 4 )   j
  character ( len = 1023 ) line
  integer   ( kind = 4 )   table(m,n)
  integer   ( kind = 4 )   x(m)

  ierror = 0

  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_file_name, status = 'old', &
    iostat = input_status )

  if ( input_status /= 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ITABLE_DATA_READ - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the input file "' // &
      trim ( input_file_name ) // '" on unit ', input_unit
    stop
  end if

  j = 0

  do while ( j < n )

    read ( input_unit, '(a)', iostat = input_status ) line

    if ( input_status /= 0 ) then
      ierror = 2
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ITABLE_DATA_READ - Fatal error!'
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
subroutine itable_header_read ( input_file_name, m, n )

!*****************************************************************************80
!
!! ITABLE_HEADER_READ reads the header from an integer table file.
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
!    Input, character ( len = * ) INPUT_FILE_NAME, the name of the input file.
!
!    Output, integer ( kind = 4 ) M, spatial dimension.
!
!    Output, integer ( kind = 4 ) N, the number of points.
!
  implicit none

  character ( len = * )  input_file_name
  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) n

  call file_column_count ( input_file_name, m )

  if ( m <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ITABLE_HEADER_READ - Fatal error!'
    write ( *, '(a)' ) '  There was some kind of I/O problem while trying'
    write ( *, '(a)' ) '  to count the number of data columns in'
    write ( *, '(a)' ) '  the file "' // trim ( input_file_name ) // '".'
    stop
  end if

  call file_row_count ( input_file_name, n )

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ITABLE_HEADER_READ - Fatal error!'
    write ( *, '(a)' ) '  There was some kind of I/O problem while trying'
    write ( *, '(a)' ) '  to count the number of data rows in'
    write ( *, '(a)' ) '  the file "' // trim ( input_file_name ) // '".'
    stop
  end if

  return
end
subroutine piecewise_linear_product_quad ( a, b, f_num, f_x, f_v, g_num, &
  g_x, g_v, quad )

!*****************************************************************************80
!
!! PIECEWISE_LINEAR_PRODUCT_QUAD: piecewise linear product integral.
!
!  Discussion:
!
!    We are given two piecewise linear functions F(X) and G(X) and we wish
!    to compute the exact value of the integral
!
!      INTEGRAL = Integral ( A <= X <= B ) F(X) * G(X) dx
!
!    The functions F(X) and G(X) are defined as tables of coordinates X and
!    values V.  A piecewise linear function is evaluated at a point X by
!    evaluating the interpolant to the data at the endpoints of the interval
!    containing X.
!
!    It must be the case that A <= B.
!
!    It must be the case that the node coordinates F_X(*) and G_X(*) are
!    given in ascending order.
!
!    It must be the case that:
!
!      F_X(1) <= A and B <= F_X(F_NUM)
!      G_X(1) <= A and B <= G_X(G_NUM)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the limits of integration.
!
!    Input, integer ( kind = 4 ) F_NUM, the number of nodes for F.
!
!    Input, real ( kind = 8 ) F_X(F_NUM), the node coordinates for F.
!
!    Input, real ( kind = 8 ) F_V(F_NUM), the nodal values for F.
!
!    Input, integer ( kind = 4 ) G_NUM, the number of nodes for G.
!
!    Input, real ( kind = 8 ) G_X(G_NUM), the node coordinates for G.
!
!    Input, real ( kind = 8 ) G_V(G_NUM), the nodal values for G.
!
!    Output, real ( kind = 8 ) QUAD, the integral of F(X) * G(X)
!    from A to B.
!
  implicit none

  integer ( kind = 4 ) f_num
  integer ( kind = 4 ) g_num

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) bit
  integer ( kind = 4 ) f_left
  real ( kind = 8 ) f_v(f_num)
  real ( kind = 8 ) f_x(f_num)
  real ( kind = 8 ) f0
  real ( kind = 8 ) f1
  real ( kind = 8 ) fl
  real ( kind = 8 ) fr
  integer ( kind = 4 ) g_left
  real ( kind = 8 ) g_v(g_num)
  real ( kind = 8 ) g_x(g_num)
  real ( kind = 8 ) g0
  real ( kind = 8 ) g1
  real ( kind = 8 ) gl
  real ( kind = 8 ) gr
  real ( kind = 8 ) h0
  real ( kind = 8 ) h1
  real ( kind = 8 ) h2
  integer ( kind = 4 ) i
  real ( kind = 8 ) quad
  real ( kind = 8 ) xl
  real ( kind = 8 ) xr
  real ( kind = 8 ) xr_max

  quad = 0.0D+00

  if ( f_x(f_num) <= a .or. g_x(g_num) <= a ) then
    return
  end if

  if ( f_num < 2 .or. g_num < 2 ) then
    return
  end if

  xr = a

  f_left = 1
  call r8vec_bracket3 ( f_num, f_x, xr, f_left )
  fr = f_v(f_left) + ( xr - f_x(f_left) ) * ( f_v(f_left+1) - f_v(f_left) ) &
    / ( f_x(f_left+1) - f_x(f_left) )

  g_left = 1
  call r8vec_bracket3 ( g_num, g_x, xr, g_left )
  gr = g_v(g_left) + ( xr - g_x(g_left) ) * ( g_v(g_left+1) - g_v(g_left) ) &
    / ( g_x(g_left+1) - g_x(g_left) )

  xr_max = b
  xr_max = min ( xr_max, f_x(f_num) )
  xr_max = min ( xr_max, g_x(g_num) )

  do while ( xr < xr_max )
!
!  Shift right values to left.
!
    xl = xr
    fl = fr
    gl = gr
!
!  Determine the new right values.
!  The hard part is figuring out how to advance XR some, but not too much.
!
    xr = xr_max

    do i = 1, 2
      if ( f_left + i <= f_num ) then
        if ( xl < f_x(f_left+i) .and. f_x(f_left+i) < xr ) then
          xr = f_x(f_left+i)
          exit
        end if
      end if
    end do

    do i = 1, 2
      if ( g_left + i <= g_num ) then
        if ( xl < g_x(g_left+i) .and. g_x(g_left+i) < xr ) then
          xr = g_x(g_left+i)
          exit
        end if
      end if
    end do

    call r8vec_bracket3 ( f_num, f_x, xr, f_left )
    fr = f_v(f_left) + ( xr - f_x(f_left) ) * ( f_v(f_left+1) - f_v(f_left) ) &
      / ( f_x(f_left+1) - f_x(f_left) )

    call r8vec_bracket3 ( g_num, g_x, xr, g_left )
    gr = g_v(g_left) + ( xr - g_x(g_left) ) * ( g_v(g_left+1) - g_v(g_left) ) &
      / ( g_x(g_left+1) - g_x(g_left) )
!
!  Form the linear polynomials for F(X) and G(X) over [XL,XR],
!  then the product H(X), integrate H(X) and add to the running total.
!
    if ( epsilon ( xl - xr ) <= abs ( xr - xl ) ) then

      f1 = fl - fr
      f0 = fr * xl - fl * xr

      g1 = gl - gr
      g0 = gr * xl - gl * xr

      h2 = f1 * g1
      h1 = f1 * g0 + f0 * g1
      h0 = f0 * g0

      h2 = h2 / 3.0D+00
      h1 = h1 / 2.0D+00

      bit = ( ( h2 * xr + h1 ) * xr + h0 ) * xr &
          - ( ( h2 * xl + h1 ) * xl + h0 ) * xl

      quad = quad + bit / ( xr - xl ) / ( xr - xl )

    end if

  end do

  return
end
subroutine r83_np_fss ( n, a, nb, b, x )

!*****************************************************************************80
!
!! R83_NP_FSS factors and solves multiple R83 systems.
!
!  Discussion:
!
!    The R83 storage format is used for a tridiagonal matrix.
!    The superdiagonal is stored in entries (1,2:N), the diagonal in
!    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
!    original matrix is "collapsed" vertically into the array.
!
!    This algorithm requires that each diagonal entry be nonzero.
!    It does not use pivoting, and so can fail on systems that
!    are actually nonsingular.
!
!  Example:
!
!    Here is how a R83 matrix of order 5 would be stored:
!
!       *  A12 A23 A34 A45
!      A11 A22 A33 A44 A55
!      A21 A32 A43 A54  *
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 November 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the linear system.
!
!    Input/output, real ( kind = 8 ) A(3,N).
!    On input, the tridiagonal matrix.
!    On output, the data in these vectors has been overwritten
!    by factorization information.
!
!    Input, integer ( kind = 4 ) NB, the number of right hand sides.
!
!    Input, real ( kind = 8 ) B(N,NB), the right hand side of the linear system.
!
!    Output, real ( kind = 8 ) X(N,NB), the solution of the linear system.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nb

  real ( kind = 8 ) a(3,n)
  real ( kind = 8 ) b(n,nb)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n,nb)
  real ( kind = 8 ) xmult
!
!  The diagonal entries can't be zero.
!
  do i = 1, n
    if ( a(2,i) == 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R83_NP_FSS - Fatal error!'
      write ( *, '(a,i8,a)' ) '  A(2,', i, ') = 0.'
      return
    end if
  end do

  x(1:n,1:nb) = b(1:n,1:nb)

  do i = 2, n
    xmult = a(3,i-1) / a(2,i-1)
    a(2,i) = a(2,i) - xmult * a(1,i)
    x(i,1:nb)   = x(i,1:nb)   - xmult * x(i-1,1:nb)
  end do

  x(n,1:nb) = x(n,1:nb) / a(2,n)
  do i = n-1, 1, -1
    x(i,1:nb) = ( x(i,1:nb) - a(1,i+1) * x(i+1,1:nb) ) / a(2,i)
  end do

  return
end
subroutine r8vec_bracket3 ( n, t, tval, left )

!*****************************************************************************80
!
!! R8VEC_BRACKET3 finds the interval containing or nearest a given value.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The routine always returns the index LEFT of the sorted array
!    T with the property that either
!    *  T is contained in the interval [ T(LEFT), T(LEFT+1) ], or
!    *  T < T(LEFT) = T(1), or
!    *  T > T(LEFT+1) = T(N).
!
!    The routine is useful for interpolation problems, where
!    the abscissa must be located within an interval of data
!    abscissas for interpolation, or the "nearest" interval
!    to the (extreme) abscissa must be found so that extrapolation
!    can be carried out.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, length of the input array.
!
!    Input, real ( kind = 8 ) T(N), an array that has been sorted
!    into ascending order.
!
!    Input, real ( kind = 8 ) TVAL, a value to be bracketed by entries of T.
!
!    Input/output, integer ( kind = 4 ) LEFT.
!    On input, if 1 <= LEFT <= N-1, LEFT is taken as a suggestion for the
!    interval [ T(LEFT), T(LEFT+1) ] in which TVAL lies.  This interval
!    is searched first, followed by the appropriate interval to the left
!    or right.  After that, a binary search is used.
!    On output, LEFT is set so that the interval [ T(LEFT), T(LEFT+1) ]
!    is the closest to TVAL; it either contains TVAL, or else TVAL
!    lies outside the interval [ T(1), T(N) ].
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) high
  integer ( kind = 4 ) left
  integer ( kind = 4 ) low
  integer ( kind = 4 ) mid
  real ( kind = 8 ) t(n)
  real ( kind = 8 ) tval
!
!  Check the input data.
!
  if ( n < 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_BRACKET3 - Fatal error!'
    write ( *, '(a)' ) '  N must be at least 2.'
    stop
  end if
!
!  If LEFT is not between 1 and N-1, set it to the middle value.
!
  if ( left < 1 .or. n - 1 < left ) then
    left = ( n + 1 ) / 2
  end if
!
!  CASE 1: TVAL < T(LEFT):
!  Search for TVAL in [T(I), T(I+1)] for intervals I = 1 to LEFT-1.
!
  if ( tval < t(left) ) then

    if ( left == 1 ) then
      return
    else if ( left == 2 ) then
      left = 1
      return
    else if ( t(left-1) <= tval ) then
      left = left - 1
      return
    else if ( tval <= t(2) ) then
      left = 1
      return
    end if
!
!  ...Binary search for TVAL in [T(I), T(I+1)] for intervals I = 2 to LEFT-2.
!
    low = 2
    high = left - 2

    do

      if ( low == high ) then
        left = low
        return
      end if

      mid = ( low + high + 1 ) / 2

      if ( t(mid) <= tval ) then
        low = mid
      else
        high = mid - 1
      end if

    end do
!
!  CASE2: T(LEFT+1) < TVAL:
!  Search for TVAL in [T(I),T(I+1)] for intervals I = LEFT+1 to N-1.
!
  else if ( t(left+1) < tval ) then

    if ( left == n - 1 ) then
      return
    else if ( left == n - 2 ) then
      left = left + 1
      return
    else if ( tval <= t(left+2) ) then
      left = left + 1
      return
    else if ( t(n-1) <= tval ) then
      left = n - 1
      return
    end if
!
!  ...Binary search for TVAL in [T(I), T(I+1)] for intervals I = LEFT+2 to N-2.
!
    low = left + 2
    high = n - 2

    do

      if ( low == high ) then
        left = low
        return
      end if

      mid = ( low + high + 1 ) / 2

      if ( t(mid) <= tval ) then
        low = mid
      else
        high = mid - 1
      end if

    end do
!
!  CASE3: T(LEFT) <= TVAL <= T(LEFT+1):
!  T is in [T(LEFT), T(LEFT+1)], as the user said it might be.
!
  else

  end if

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
!    Output, integer ( kind = 4 ) LENGTH, the number of characters of S
!    used to make IVAL.
!
  implicit none

  character              c
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) ierror
  integer   ( kind = 4 ) isgn
  integer   ( kind = 4 ) istate
  integer   ( kind = 4 ) ival
  integer   ( kind = 4 ) length
  character ( len = * )  s

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

  integer   ( kind = 4 ) n

  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) ierror
  integer   ( kind = 4 ) ilo
  integer   ( kind = 4 ) ivec(n)
  integer   ( kind = 4 ) length
  character ( len = * )  s

  i = 0
  ierror = 0
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

  character              c
  logical                ch_eqi
  real ( kind = 8 ) dval
  integer   ( kind = 4 ) ierror
  integer   ( kind = 4 ) ihave
  integer   ( kind = 4 ) isgn
  integer   ( kind = 4 ) iterm
  integer   ( kind = 4 ) jbot
  integer   ( kind = 4 ) jsgn
  integer   ( kind = 4 ) jtop
  integer   ( kind = 4 ) length
  integer   ( kind = 4 ) nchar
  integer   ( kind = 4 ) ndig
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

  integer   ( kind = 4 ) n

  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) ierror
  integer   ( kind = 4 ) ilo
  integer   ( kind = 4 ) lchar
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
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) lens
  integer   ( kind = 4 ) nword
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
  integer   ( kind = 4 ) d
  integer   ( kind = 4 ) h
  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) mm
  character ( len = 9  ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer   ( kind = 4 ) n
  integer   ( kind = 4 ) s
  integer   ( kind = 4 ) values(8)
  integer   ( kind = 4 ) y

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
