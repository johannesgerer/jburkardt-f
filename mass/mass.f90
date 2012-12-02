program main

!*****************************************************************************80
!
!! MAIN is the main program for MASS.
!
!  Discussion:
!
!    MASS is a program to check the mass matrix computation.
!
!    This code was extracted from the POD_BASIS program, in order to
!    build and test the code that computes, factors, and uses
!    the mass matrix.
!
!    This is a tricky process, because the original mass matrix
!    is computed in another executable program, not controlled by
!    me but by H C Lee.  Thus, an appalling number of choices must
!    match, including conventions about the ordering of nodes,
!    the quadrature rule, and so on. 
!
!    After running this program many times, I think we finally
!    worked out the set of parameters to be used.  The information
!    gained from this exercise was incorporated in the "real"
!    program, and this one was retired at the end of July, 2003.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 August 2003
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable, dimension ( :, : ) :: a
  real ( kind = 8 ) a_dot_p
  integer ( kind = 4 ) bandwidth
  character ( len = 255 ) basis_file
  integer ( kind = 4 ) basis_num
  logical ch_eqi
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) distance_type
  character ( len = 255 ) element_file_name
  integer ( kind = 4 ) element_num
  real ( kind = 8 ), allocatable, dimension ( : ) :: energy
  integer ( kind = 4 ) energy_it_max
  logical file_exist
  character ( len = 255 ) gen_file
  integer ( kind = 4 ) gen_unit
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) info
  integer ( kind = 4 ) inunit
  integer ( kind = 4 ) j
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: node
  integer ( kind = 4 ) node_num
  real ( kind = 8 ) norm
  integer ( kind = 4 ) normal
  integer ( kind = 4 ) npe
  integer ( kind = 4 ) nrhs
  real ( kind = 8 ), allocatable, dimension ( :, :) :: point
  integer ( kind = 4 ) point_num
  integer ( kind = 4 ) point_num2
  integer ( kind = 4 ) run_type
  logical s_eqi
  character ( len = 11 ) s_of_i4
  character ( len = 255 ) steady_file
  real ( kind = 8 ) steady_max
  real ( kind = 8 ) steady_norm
  real ( kind = 8 ), allocatable, dimension ( : ) :: sval
  integer ( kind = 4 ) swap_num
  integer ( kind = 4 ) thin
  real ( kind = 8 ), allocatable, dimension ( : ) :: u
  real ( kind = 8 ), allocatable, dimension ( : ) :: u_steady
  character ( len = 255 ) uv_file
  integer ( kind = 4 ) uv_file_num
  character ( len = 255 ) uv0_file
  real ( kind = 8 ), allocatable, dimension ( : ) :: v
  real ( kind = 8 ), allocatable, dimension ( : ) :: v_steady
  real ( kind = 8 ), allocatable, dimension ( : ) :: x
  real ( kind = 8 ) x_max
  real ( kind = 8 ) x_min
  character ( len = 255 ) xy_file
  integer ( kind = 4 ) xy_lines
  integer ( kind = 4 ) xy_values
  real ( kind = 8 ), allocatable, dimension ( : ) :: y
  real ( kind = 8 ) y_max
  real ( kind = 8 ) y_min

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MASS'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)') ' '
  write ( *, '(a)' ) '  Given a PDE for which:'
  write ( *, '(a)' ) '    M is the dimension of each solution vector,'
  write ( *, '(a)' ) '    N is the number of solution vectors,'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Set up A, the M by N matrix of solution vectors,'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Get A = U * S * VT, the singular value decomposition.'
  write ( *, '(a)' ) ' '
!
!  Get the run type
!
  call i4_input ( 'What is the run type?', run_type, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MASS - Fatal error!'
    write ( *, '(a)' ) '  Input error reading the run type.'
    stop
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  RUN_TYPE = ', run_type
!
!  What is the basis size?
!
  call i4_input ( 'What is the requested basis size?', basis_num, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MASS - Fatal error!'
    write ( *, '(a)' ) '  Input error reading the basis size.'
    stop
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  BASIS_NUM = ', basis_num

  if ( run_type == 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  For RUN_TYPE = 1,'
    write ( *, '(a)' ) '  read in the steady state solution'
    write ( *, '(a)' ) '  and, letting SS be the steady state solution,'
    write ( *, '(a)' ) '  subtract 5/3 SS from solutions 1 through 250'
    write ( *, '(a)' ) '  subtract 1/3 SS from solutions 251 through 500.'
  end if
!
!  Get the XY data file.
!
  call s_input ( 'What is the XY data file name?', xy_file, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MASS - Fatal error!'
    write ( *, '(a)' ) '  Input error reading the XY data file name.'
    stop
  end if

  call data_size ( xy_file, xy_lines, xy_values, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MASS - Fatal error!'
    write ( *, '(a)' ) '  Error reading the XY data file.'
    stop
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The file "' // trim ( xy_file ) // '" contains ' &
    // trim ( s_of_i4 ( xy_lines ) ) // ' lines.'
!
!  Allocate space for some arrays.
!
  node_num = xy_lines

  allocate ( u(1:node_num) )
  allocate ( u_steady(1:node_num) )
  allocate ( v(1:node_num) )
  allocate ( v_steady(1:node_num) )
  allocate ( x(1:node_num) )
  allocate ( y(1:node_num) )
!
!  Read in X and Y.
!
  call data_r82_read ( xy_file, node_num, x, y, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MASS - Fatal error!'
    write ( *, '(a)' ) '  Error reading from the XY data file.'
    stop
  end if
!
!  Extract some interesting data.
!
  x_min = minval ( x(1:node_num) )
  x_max = maxval ( x(1:node_num) )
  y_min = minval ( y(1:node_num) )
  y_max = maxval ( y(1:node_num) )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  X minimum : ', x_min
  write ( *, '(a,g14.6)' ) '  X maximum : ', x_max
  write ( *, '(a,g14.6)' ) '  Y minimum : ', y_min
  write ( *, '(a,g14.6)' ) '  Y maximum : ', y_max
!
!  Get the steady state file name.
!
  call s_input ( 'What is the steady state file name?', steady_file, &
    ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'POD_BASIS - Fatal error!'
    write ( *, '(a)' ) '  Input error reading the steady state file name.'
    stop
  end if
!
!  Read in the steady state solution.
!
  call data_r82_read ( steady_file, node_num, u_steady, v_steady, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MASS - Fatal error!'
    write ( *, '(a)' ) '  Error reading the steady state data file.'
    stop
  end if

  steady_max = &
    maxval ( sqrt ( u_steady(1:node_num)**2 + v_steady(1:node_num)**2 ) )
  steady_norm = &
    sqrt ( sum ( u_steady(1:node_num)**2 + v_steady(1:node_num)**2 ) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Steady state information was read from'
  write ( *, '(a)' ) '  the file "' // trim ( steady_file ) // '".'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Steady max norm = ', steady_max
  write ( *, '(a,g14.6)' ) '  Steady l2 norm =  ', steady_norm
!
!  Get the UV0 file name.
!
  call s_input ( 'What is the first solution file name?', uv0_file, &
    ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'POD_BASIS - Fatal error!'
    write ( *, '(a)' ) '  Input error reading the first solution file name.'
    stop
  end if
!
!  Presumably, all the solution files have the same name as the first
!  solution file, but with a numerical increment.  To begin with, simply count
!  the number of files.
!
  uv_file = uv0_file
  uv_file_num = 0

  do

    if ( .not. file_exist ( uv_file ) ) then
      exit
    end if

    uv_file_num = uv_file_num + 1

    call file_name_inc ( uv_file )

  end do

  if ( uv_file_num == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MASS - Fatal error!'
    write ( *, '(a)' ) '  There do not seem to be any solution files;'
    write ( *, '(a)' ) '  that is, files whose names are "incremented"'
    write ( *, '(a)' ) '  versions of the steady state file name.'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The first file we looked for was "' // &
      trim ( uv_file ) // '".'
    stop
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) &
    'We believe the number of solution files is ', uv_file_num
!
!  Now we have enough information to set up a data structure.
!
!  Determine the spatial dimension (columns) and number of points (rows).
!
  dim_num = 2 * node_num

  point_num = uv_file_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The data is stored in an M by N matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  The "spatial" dimension M is   ', dim_num
  write ( *, '(a,i6)' ) '  The number of data points N is ', point_num
!
!  Allocate space for the POINT array.
!
  allocate ( point(1:dim_num,1:point_num) )
!
!  Now read the data from the individual files, process it if necessary,
!  and gather it into a single array called POINT.
!
  uv_file = uv0_file

  do j = 1, point_num

    call data_r82_read ( uv_file, node_num, u, v, ierror )

    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'POD_BASIS - Fatal error!'
      write ( *, '(a)' ) '  Error reading the solution data file.'
      stop
    end if
!
!  For RUN_TYPE = 1,
!  subtract 5/3 of SS from solutions 1-250, 
!  subtract 1/3 of SS from solutions 251-500.
!
    if ( run_type == 1 ) then

      if ( 1 <= j .and. j <= 250 ) then
        u(1:node_num) = u(1:node_num) - 5.0D+00 * u_steady(1:node_num) / 3.0D+00
        v(1:node_num) = v(1:node_num) - 5.0D+00 * v_steady(1:node_num) / 3.0D+00
      else if ( 251 <= j .and. j <= 500 ) then
        u(1:node_num) = u(1:node_num) - u_steady(1:node_num) / 3.0D+00
        v(1:node_num) = v(1:node_num) - v_steady(1:node_num) / 3.0D+00
      end if

    end if

    point(1:2*node_num-1:2,j) = u(1:node_num)
    point(2:2*node_num  :2,j) = v(1:node_num)

    call file_name_inc ( uv_file )

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'All the data has been read into POINT.'
!
!----------------------------------------------------------------------------
!
!  Precondition with the mass matrix M?
!
!  If so, instead of the system X' * X, we wish to study the
!  system X' * M * X. 
!
!  However, it's easy to apply either eigenanalsis to X'*X or
!  SVD to X, however, if we are dealing with the system X' * M * X,
!  it's easy to apply eigenanalysis, but not obvious what system
!  we should apply SVD to.
!
!  In order to make this look like the simpler problem, we need
!  to split M up into appropriate factors and regroup.
!
!  Thus, we will:
!  * compute the (symmetric positive definite) mass matrix M,
!  * determine the Cholesky factorization M = L * L',
!  * premultiply the data matrix X by L'.
!
!  Thus, we regard the system 
!    X' * M * X
!  as the system
!    X' * L * L' * X
!  or
!    (L'*X)' * (L'*X)
!
!  Now it's obvious that we will apply SVD to L'*X.
!
!  THEN (10 July) we must postprocess the vectors we extract from
!  the columns of U by essentially multiplying them by inverse (L').
!
!  So if we overwrite X by L'X, we can proceed as though nothing
!  is really different.
!
!----------------------------------------------------------------------------
!
  call s_input ( &
    'Enter element file for mass matrix preconditioning or "None".', &
    element_file_name, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MASS - Fatal error!'
    write ( *, '(a)' ) '  Input error reading the element file name.'
    stop
  end if

  if ( s_eqi ( element_file_name, 'None' ) ) then

  else

    call data_size ( element_file_name, element_num, npe, ierror )

    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MASS - Fatal error!'
      write ( *, '(a)' ) '  Input error reading the element file.'
      stop
    end if

    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) '  Number of elements = ', element_num
    write ( *, '(a,i6)' ) '  Number of nodes per element = ', npe

    allocate ( node(1:npe,1:element_num) )

    call data_i4vec_read ( element_file_name, npe, element_num, node, ierror )

    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MASS - Fatal error!'
      write ( *, '(a)' ) '  Error reading the element file.'
      stop
    end if

    call bandwidth_determine ( npe, element_num, node, bandwidth )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) '  The bandwidth of the matrix is ', bandwidth
!
!  Allocate storage for A.
!
    allocate ( a(bandwidth+1,node_num) )
!
!  Compute A.
!
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Compute the mass matrix.'

    call mass_matrix ( node_num, npe, element_num, node, bandwidth, x, y, a )
!
!  Print the mass matrix, if small.
!
    call r8pbl_write ( node_num, bandwidth, a, 'mass2.txt' )

    deallocate ( node )

    deallocate ( a )

  end if
!
!  Free memory.
!
  deallocate ( point )
  deallocate ( u )
  deallocate ( u_steady )
  deallocate ( v )
  deallocate ( v_steady )
  deallocate ( x )
  deallocate ( y )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MASS'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine bandwidth_determine ( npe, element_num, node, bandwidth )

!*****************************************************************************80
!
!! BANDWIDTH_DETERMINE computes the lower bandwidth of a finite element matrix.
!
!  Discussion:
!
!    The finite element matrix is assumed to be structured in such a
!    way that row I represents an equation associated with the unknown
!    at node I, and column I stores coefficients associated with
!    the unknown at node I as it appears in various equations.
!
!    Further, it is assumed that the I-th equation, associated with
!    node and unknown I, involves only those nodes and unknowns J
!    with the property that there is an element K that includes both
!    nodes.
!
!    Thus, the (half) bandwidth calculation simply involves finding the
!    greatest difference between two nodes in the same element and adding 1.
!
!    A diagonal matrix will have bandwidth 1.  A tridiagonal matrix
!    will have bandwidth 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NPE, the number of nodes per element.
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of elements.
!
!    Input, integer ( kind = 4 ) NODE(NPE,ELEMENT_NUM), the nodes that 
!    make up each element.
!
!    Output, integer ( kind = 4 ) BANDWIDTH, the (half) bandwidth of the matrix.
!
  implicit none

  integer ( kind = 4 ) element_num
  integer ( kind = 4 ) npe

  integer ( kind = 4 ) bandwidth
  integer ( kind = 4 ) elem
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) node(npe,element_num)

  bandwidth = 1

  do elem = 1, element_num
    do i1 = 1, npe
      n1 = node(i1,elem)
      do i2 = 1, npe
        n2 = node(i2,elem)
        bandwidth = max ( bandwidth, n2 + 1 - n1 )
      end do
    end do
  end do

  return
end
subroutine basis_write ( file_out_name, m, s, x )

!*****************************************************************************80
!
!! BASIS_WRITE writes a basis vector to a file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_OUT_NAME, the name of the file to write.
!
!    Input, integer ( kind = 4 ) M, the number of data items.
!
!    Input, real S, the associated singular value.
!
!    Input, real X(M), the data values.
!
  implicit none

  integer ( kind = 4 ) m

  character ( len = * ) file_out_name
  integer ( kind = 4 ) file_out_unit
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ios
  real ( kind = 8 ) s
  character ( len = 40 ) string
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) x(m)

  call get_unit ( file_out_unit )

  open ( unit = file_out_unit, file = file_out_name, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BASIS_WRITE - Fatal error!'
    write ( *, '(a)' ) '  Error opening output file "' // trim ( file_out_name ) // '".'
    stop
  end if

  write ( file_out_unit, '(a)'       ) '#  ' // trim ( file_out_name )
  write ( file_out_unit, '(a)'       ) '#  created by BASIS_WRITE.'
  write ( file_out_unit, '(a)'       ) '#'
  write ( file_out_unit, '(a,i6)'    ) '#  Velocity vectors.'
  write ( file_out_unit, '(a)'       ) '#'
  write ( file_out_unit, '(a,i6)'    ) '#  Number of records = ', m / 2
  write ( file_out_unit, '(a,g14.6)' ) '#  Singular value S = ', s
  write ( file_out_unit, '(a)'       ) '#'

  do i = 1, m, 2

    if ( abs ( x(i) ) < 1.0D-10 ) then
      u = 0.0D+00
    else
      u = x(i)
    end if

    if ( abs ( x(i+1) ) < 1.0D-10 ) then
      v = 0.0D+00
    else
      v = x(i+1)
    end if

    write ( file_out_unit, '(2e25.15)' ) u, v

  end do

  close ( unit = file_out_unit )

  return
end
subroutine ch_cap ( ch )

!*****************************************************************************80
!
!! CH_CAP capitalizes a single character.
!
!  Discussion:
!
!    Instead of CHAR and ICHAR, we now use the ACHAR and IACHAR functions, 
!    which guarantee the ASCII collating sequence.
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
!    Input/output, character CH, the character to capitalize.
!
  implicit none

  character ch
  integer ( kind = 4 ) itemp

  itemp = iachar ( ch )
 
  if ( 97 <= itemp .and. itemp <= 122 ) then
    ch = achar ( itemp - 32 )
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
!    CH_EQI ( 'A', 'a' ) is TRUE.
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
function ch_is_digit ( ch )

!*****************************************************************************80
!
!! CH_IS_DIGIT is TRUE if a character is a decimal digit.
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
!    Input, character CH, the character to be analyzed.
!
!    Output, logical CH_IS_DIGIT, is TRUE if the character is a digit.
!
  implicit none

  character ch
  logical ch_is_digit

  if ( lle ( '0', ch ) .and. lle ( ch, '9' ) ) then
    ch_is_digit = .true.
  else
    ch_is_digit = .false.
  end if
 
  return
end
subroutine ch_to_digit ( ch, digit )

!*****************************************************************************80
!
!! CH_TO_DIGIT returns the integer value of a base 10 digit.
!
!  Discussion:
!
!    Instead of ICHAR, we now use the IACHAR function, which
!    guarantees the ASCII collating sequence.
!
!  Example:
!
!     CH  DIGIT
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
!    Input, character CH, the decimal digit, '0' through '9' or blank
!    are legal. 
!
!    Output, integer ( kind = 4 ) DIGIT, the corresponding value.  
!    If CH was 'illegal', then DIGIT is -1.
!
  implicit none

  character ch
  integer ( kind = 4 ) digit

  if ( lle ( '0', ch ) .and. lle ( ch, '9' ) ) then
 
    digit = iachar ( ch ) - 48
 
  else if ( ch == ' ' ) then
 
    digit = 0
 
  else

    digit = -1

  end if
 
  return
end
subroutine data_i4vec_read ( file_in_name, m, n, a, ierror )

!*****************************************************************************80
!
!! DATA_I4VEC_READ reads an dataset of integer vectors stored in a file.
!
!  Discussion:
!
!    The data set can be thought of as an integer M by N array.
!
!    Each column of the array corresponds to one data "item".
!
!    The data is stored in a file, one column at a time.
!
!    Each data item begins on a new line.
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
!    22 April 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_IN_NAME, the name of the file to read.
!
!    Input, integer ( kind = 4 ) M, the size of each data item.
!
!    Input, integer ( kind = 4 ) N, the number of data items.
!
!    Output, integer ( kind = 4 ) A(M,N), the data.
!
!    Output, integer ( kind = 4 ) IERROR, is nonzero if an error occurred.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  character ( len = * ) file_in_name
  integer ( kind = 4 ) file_in_unit
  integer ( kind = 4 ) file_in_line
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) last
  integer ( kind = 4 ) length
  character ( len = 80 ) line
  integer ( kind = 4 ) n2

  ierror = 0

  call get_unit ( file_in_unit )

  open ( unit = file_in_unit, file = file_in_name, iostat = ios, &
    status = 'old' )

  if ( ios /= 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DATA_I4VEC_READ - Fatal error!'
    write ( *, '(a)' ) '  Error opening the data file"' // trim ( file_in_name ) // '".'
    return
  end if

  n2 = 0
  file_in_line = 0

  do

    read ( file_in_unit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      ierror = 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DATA_I4VEC_READ - Fatal error!'
      write ( *, '(a)' ) '  Error reading the data file "' // trim ( file_in_name ) // '".'
      return
    end if

    file_in_line = file_in_line + 1

    if ( len_trim ( line ) == 0 ) then

    else if ( line(1:1) == '#' ) then

    else

      n2 = n2 + 1

      last = 0

      do i = 1, m

        call s_to_i4 ( line(last+1:), a(i,n2), ierror, length )

        if ( ierror /= 0 ) then
          ierror = 1
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'DATA_I4VEC_READ - Fatal error!'
          write ( *, '(a)' ) '  Error reading data.'
          write ( *, '(a)' ) '  Reading line number ', file_in_line
          write ( *, '(a)' ) '  Data row number ', i
          write ( *, '(a)' ) '  Data column ', n2
          return
        end if

        last = last + length

      end do

      if ( n <= n2 ) then
        exit
      end if

    end if

  end do

  close ( unit = file_in_unit )

  return
end
subroutine data_r82_read ( file_name, n, x, y, ierror )

!*****************************************************************************80
!
!! DATA_R82_READ reads a data set of pairs of R8 numbers stored in a file.
!
!  Discussion:
!
!    The data set can be thought of as a real M by 2 array.
!
!    Each row of the array corresponds to one data "item".
!
!    The data is stored in a file, one row (pair of values) at a time.
!
!    Each row (pair of values) begins on a new line.
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
!    22 April 2002
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
!    Output, real ( kind = 8 ) X(N), Y(N), the data values.
!
!    Output, integer ( kind = 4 ) IERROR, is nonzero if an error occurred.
!
  implicit none

  integer ( kind = 4 ) n

  character ( len = * ) file_name
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) input
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) last
  integer ( kind = 4 ) length
  character ( len = 80 ) line
  integer ( kind = 4 ) n2
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)

  ierror = 0

  call get_unit ( input )

  open ( unit = input, file = file_name, iostat = ios, status = 'old' )

  if ( ios /= 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DATA_R82_READ - Fatal error!'
    write ( *, '(a)' ) '  Error opening the file "' // trim ( file_name ) // '".'
    return
  end if

  x(1:n) = huge ( x(1) )
  y(1:n) = huge ( y(1) )

  n2 = 0

  do

    read ( input, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      ierror = 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DATA_R82_READ - Fatal error!'
      write ( *, '(a)' ) '  Error reading the file"' // trim ( file_name ) // '".'
      return
    end if

    if ( len_trim ( line ) == 0 ) then

    else if ( line(1:1) == '#' ) then

    else

      n2 = n2 + 1

      last = 0
      call s_to_r8 ( line(last+1:), x(n2), ierror, length )

      if ( ierror /= 0 ) then
        ierror = 1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'DATA_R82_READ - Fatal error!'
        write ( *, '(a)' ) '  Error decoding information in '
        write ( *, '(a)' ) '  the file "' // trim ( file_name ) // '".'
        return
      end if

      last = last + length

      call s_to_r8 ( line(last+1:), y(n2), ierror, length )

      if ( ierror /= 0 ) then
        ierror = 1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'DATA_R82_READ - Fatal error!'
        write ( *, '(a)' ) '  Error decoding information in the file.'
        return
      end if

      if ( n2 == n ) then
        exit
      end if

    end if

  end do

  close ( unit = input )

  return
end
subroutine data_r82_write ( file_name, m, x1, x2 )

!*****************************************************************************80
!
!! DATA_R82_WRITE writes a data set of pairs of real*8 numbers into a file.
!
!  Discussion:
!
!    The data set can be thought of as a real M by 2 array.
!
!    Each row of the array corresponds to one data "item".
!
!    The data is stored in a file, one row (pair of values) at a time.
!
!    Each row (pair of values) begins on a new line.
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
!    20 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the file to write.
!
!    Input, integer ( kind = 4 ) M, the number of data items.
!
!    Input, real ( kind = 8 ) X1(M), X2(M), the data values.
!
  implicit none

  integer ( kind = 4 ) m

  character ( len = * ) file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) output
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) x1(m)
  real ( kind = 8 ) x2(m)

  call get_unit ( output )

  open ( unit = output, file = file_name, status = 'replace' )

  do i = 1, m
!
!  Professor Lee requests D25.15 format.
!
!   write ( output, '(2d25.15)' ) x1(i), x2(i)

    if ( abs ( x1(i) ) < 1.0D-10 ) then
      u = 0.0D+00
    else
      u = x1(i)
    end if

    if ( abs ( x2(i) ) < 1.0D-10 ) then
      v = 0.0D+00
    else
      v = x2(i)
    end if

    write ( output, '(2e25.15)' ) u, v

  end do

  close ( unit = output )

  return
end
subroutine data_size ( file_name, m, n, ierror )

!*****************************************************************************80
!
!! DATA_SIZE counts the size of a data set stored in a file.
!
!  Discussion:
!
!    Blank lines and comment lines (which begin with '#') are ignored).
!
!    All other lines are assumed to represent data items.
!
!    This routine assumes that each data line contains exactly the
!    same number of values, which are separated by spaces.
!
!    (This means this routine cannot handle cases where a data item
!    extends over more than one line, or cases where data is squeezed
!    together with no spaces, or where commas are used as separators,
!    but with space on both sides.)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 April 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the file to read.
!
!    Output, integer ( kind = 4 ) M, the number of nonblank, noncomment lines.
!
!    Output, integer ( kind = 4 ) N, the number of values per line.
!
!    Output, integer ( kind = 4 ) IERROR, is nonzero if an error occurred.
!
  implicit none

  character ( len = * ) file_name
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) input
  integer ( kind = 4 ) ios
  character ( len = 80 ) line
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_max
  integer ( kind = 4 ) n_min
  integer ( kind = 4 ) n_word

  ierror = 0
  m = 0
  n_max = - huge ( n_max )
  n_min = huge ( n_min )

  call get_unit ( input )

  open ( unit = input, file = file_name, status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DATA_SIZE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the file "' // trim ( file_name ) // '".'
    return
  end if

  do

    read ( input, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      exit
    end if

    if ( len_trim ( line ) == 0 ) then

    else if ( line(1:1) == '#' ) then

    else

      m = m + 1

      call s_word_count ( line, n_word )

      n_max = max ( n_max, n_word )
      n_min = min ( n_min, n_word )

    end if

  end do

  if ( n_max /= n_min ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DATA_SIZE - Fatal error!'
    write ( *, '(a)' ) '  Number of words per line varies.'
    write ( *, '(a,i6)' ) '  Minimum is ', n_min
    write ( *, '(a,i6)' ) '  Maximum is ', n_max
    n = 0
  else
    n = n_min
  end if

  close ( unit = input )

  return
end
subroutine digit_inc ( ch )

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
!    Input/output, character CH, a digit to be incremented.
!
  implicit none

  character ch
  integer ( kind = 4 ) digit

  call ch_to_digit ( ch, digit )

  if ( digit == -1 ) then
    return
  end if

  digit = digit + 1

  if ( digit == 10 ) then
    digit = 0
  end if

  call digit_to_ch ( digit, ch )

  return
end
subroutine digit_to_ch ( digit, ch )

!*****************************************************************************80
!
!! DIGIT_TO_CH returns the character representation of a decimal digit.
!
!  Discussion:
!
!    Instead of CHAR, we now use the ACHAR function, which
!    guarantees the ASCII collating sequence.
!
!  Example:
!
!    DIGIT   CH 
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
!    Output, character CH, the corresponding character.
!
  implicit none

  character ch
  integer ( kind = 4 ) digit

  if ( 0 <= digit .and. digit <= 9 ) then

    ch = achar ( digit + 48 )

  else

    ch = '*'

  end if
 
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

  character ( len = * ) file_name
  logical file_exist

  inquire ( file = file_name, exist = file_exist )

  return
end
subroutine file_name_inc ( file_name )

!*****************************************************************************80
!
!! FILE_NAME_INC increments a partially numeric filename.
!
!  Discussion:
!
!    It is assumed that the digits in the name, whether scattered or
!    connected, represent a number that is to be increased by 1 on
!    each call.  If this number is all 9's on input, the output number
!    is all 0's.  Non-numeric letters of the name are unaffected.
!
!    If the name is empty, then the routine stops.
!
!    If the name contains no digits, the empty string is returned.
!
!  Example:
!
!      Input            Output
!      -----            ------
!      'a7to11.txt'     'a7to12.txt'
!      'a7to99.txt'     'a8to00.txt'
!      'a9to99.txt'     'a0to00.txt'
!      'cat.txt'        ' '
!      ' '              STOP!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 September 2005
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
  integer ( kind = 4 ) change
  integer ( kind = 4 ) digit
  character ( len = * ) file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) lens

  lens = len_trim ( file_name )

  if ( lens <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_NAME_INC - Fatal error!'
    write ( *, '(a)' ) '  The input string is empty.'
    stop
  end if

  change = 0

  do i = lens, 1, -1

    c = file_name(i:i)

    if ( lge ( c, '0' ) .and. lle ( c, '9' ) ) then

      change = change + 1

      digit = ichar ( c ) - 48
      digit = digit + 1

      if ( digit == 10 ) then
        digit = 0
      end if

      c = char ( digit + 48 )

      file_name(i:i) = c

      if ( c /= '0' ) then
        return
      end if

    end if

  end do

  if ( change == 0 ) then
    file_name = ' '
    return
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
subroutine i4_input ( string, value, ierror )

!*****************************************************************************80
!
!! I4_INPUT prints a prompt string and reads an I4 from the user.
!
!  Discussion:
!
!    If the input line starts with a comment character ('#') or is
!    blank, the routine ignores that line, and tries to read the next one.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) STRING, the prompt string.
!
!    Output, integer ( kind = 4 ) VALUE, the value input by the user.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag, which is zero 
!    if no error occurred.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) last
  character ( len = 80 ) line
  character ( len = * ) string
  integer ( kind = 4 ) value

  ierror = 0
  value = huge ( value )
!
!  Write the prompt.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( string )

  do

    read ( *, '(a)', iostat = ierror ) line

    if ( ierror /= 0 ) then
      return
    end if
!
!  If the line begins with a comment character, go back and read the next line.
!
    if ( line(1:1) == '#' ) then
      cycle
    end if

    if ( len_trim ( line ) == 0 ) then
      cycle
    end if
!
!  Extract integer information from the string.
!
    call s_to_i4 ( line, value, ierror, last )

    if ( ierror /= 0 ) then
      value = huge ( value )
      return
    end if

    exit

  end do

  return
end
subroutine i4_range_input ( string, value1, value2, ierror )

!*****************************************************************************80
!
!! I4_RANGE_INPUT reads a pair of I4's from the user, representing a range.
!
!  Discussion:
!
!    If the input line starts with a comment character ('#') or is blank,
!    the routine ignores that line, and tries to read the next one.
!
!    The pair of integers may be separated by spaces or a comma or both.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) STRING, the prompt string.
!
!    Output, integer ( kind = 4 ) VALUE1, VALUE2, the values entered by
!    the user.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag, which is zero 
!    if no error occurred.
!
  implicit none

  character, parameter :: comma = ','
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) last
  integer ( kind = 4 ) last2
  character ( len = 80 ) line
  character, parameter :: space = ' '
  character ( len = * ) string
  integer ( kind = 4 ) value1
  integer ( kind = 4 ) value2

  ierror = 0
  value1 = huge ( value1 )
  value2 = huge ( value2 )
!
!  Write the prompt.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( string )

  do

    read ( *, '(a)', iostat = ierror ) line

    if ( ierror /= 0 ) then
      return
    end if
!
!  If the line begins with a comment character, go back and read the next line.
!
    if ( line(1:1) == '#' ) then
      cycle
    end if

    if ( len_trim ( line ) == 0 ) then
      cycle
    end if
!
!  Replace commas by spaces.
!
    call s_replace_ch ( line, comma, space )
!
!  Extract integer information from the string.
!
    call s_to_i4 ( line, value1, ierror, last )

    if ( ierror /= 0 ) then
      value1 = huge ( value1 )
      return
    end if

    call s_to_i4 ( line(last+1:), value2, ierror, last2 )

    if ( ierror /= 0 ) then
      value2 = huge ( value2 )
      return
    end if

    exit

  end do

  return
end
subroutine mass_matrix ( node_num, npe, element_num, node, bandwidth, x, y, a )

!*****************************************************************************80
!
!! MASS_MATRIX computes the mass matrix.
!
!  Discussion:
!
!    I want to compute the mass matrix associated with velocity.
!
!      A(I,J) = integral ( PHI(I)(X,Y) * PHI(J)(X,Y) ) d Region
!
!    where PHI(I) and PHI(J) are the shape functions associated with
!    the I-th and J-th variables.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 May 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, integer ( kind = 4 ) NPE, the number of nodes per element.
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of elements.
!
!    Input, integer ( kind = 4 ) NODE(NPE,ELEMENT_NUM), the nodes that 
!    make up each element.
!
!    Input, integer ( kind = 4 ) BANDWIDTH, the half bandwidth of the matrix.
!
!    Input, real ( kind = 8 ) X(NODE_NUM), Y(NODE_NUM), the coordinates 
!    of the nodes.
!
!    Output, real ( kind = 8 ) A(BANDWIDTH+1,NODE_NUM), the mass matrix.
!
  implicit none

  integer ( kind = 4 ) bandwidth
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) npe
  integer ( kind = 4 ), parameter :: nquad = 13

  real ( kind = 8 ) a(bandwidth+1,node_num)
  real ( kind = 8 ) area
  integer ( kind = 4 ) element
  real ( kind = 8 ) eta
  real ( kind = 8 ) eta_tab(6)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) iq
  integer ( kind = 4 ) iquad
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jp
  integer ( kind = 4 ) jq
  integer ( kind = 4 ) nlband
  integer ( kind = 4 ) node(npe,element_num)
  integer ( kind = 4 ) norder
  integer ( kind = 4 ) p1
  integer ( kind = 4 ) p2
  integer ( kind = 4 ) p3
  real ( kind = 8 ) refqbf
  integer ( kind = 4 ) rule
  real ( kind = 8 ) w(6)
  real ( kind = 8 ) weight(nquad)
  real ( kind = 8 ) wi
  real ( kind = 8 ) wj
  real ( kind = 8 ) x(node_num)
  real ( kind = 8 ) xsi
  real ( kind = 8 ) xsi_tab(6)
  real ( kind = 8 ) xtab(nquad)
  real ( kind = 8 ) y(node_num)
  real ( kind = 8 ) ytab(nquad)
!
!  Zero out the matrix.
!
  a(1:bandwidth+1,1:node_num) = 0.0D+00
!
!  Get the weights and abscissas for a unit triangle.
!
  rule = 12
  call triangle_unit_set ( rule, norder, xtab, ytab, weight )
!
!  For each element.
!
  do element = 1, element_num

    p1 = node(1,element)
    p2 = node(2,element)
    p3 = node(3,element)

    area = 0.5D+00 * abs ( &
        x(p1) * ( y(p2) - y(p3) ) &
      + x(p2) * ( y(p3) - y(p1) ) &
      + x(p3) * ( y(p1) - y(p2) ) )

    if ( area == 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MASS_MATRIX - Fatal error!'
      write ( *, '(a,i6)' ) '  Zero area for element ', element
      stop
    end if
!
!  For each quadrature point in the element...
!
    do iquad = 1, nquad

      xsi = xtab(iquad)
      eta = ytab(iquad)
!
!  For each basis function PHI(I) associated with a node in the element,
!
      do iq = 1, 6

        ip = node(iq,element)
        wi = refqbf ( iq, xsi, eta )
!
!  For each "neighbor" basis function PHI(J) associated with a node in 
!  the element.
!
        do jq = 1, 6

          jp = node(jq,element)

          if ( jp <= ip ) then
            wj = refqbf ( jq, xsi, eta )
            a(ip+1-jp,jp) = a(ip+1-jp,jp) + area * weight(iquad) * wi * wj
          end if

        end do
      end do
    end do
  end do

  return
end
subroutine node_t6 ( r, s )

!*****************************************************************************80
!
!! NODE_T6 returns the basis nodes for a 6 node triangle.
!
!  Diagram:
!
!    |
!    1  3
!    |  |\
!    |  | \
!    S  6  5
!    |  |   \
!    |  |    \
!    0  1--4--2
!    |
!    +--0--R--1-->
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) R(6), S(6), the coordinates of the basis nodes.
!
  implicit none

  real ( kind = 8 ) r(6)
  real ( kind = 8 ) s(6)

  r(1:6) = (/ 0.0D+00, 1.0D+00, 0.0D+00, 0.5D+00, 0.5D+00, 0.0D+00 /)
  s(1:6) = (/ 0.0D+00, 0.0D+00, 1.0D+00, 0.0D+00, 0.5D+00, 0.5D+00 /)

  return
end
subroutine r8blt_check ( n, ml, ierror )

!*****************************************************************************80
!
!! R8BLT_CHECK checks the dimensions of a banded lower triangular matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 May 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) ML, the lower bandwidth of the matrix.
!    ML must be at least 0, and no greater than N-1.
!
!    Output, integer ( kind = 4 ) IERROR, reports whether any errors were detected.
!    IERROR is set to 0 before the checks are made, and then:
!    IERROR = IERROR + 1 if ML is illegal;
!    IERROR = IERROR + 2 if N is illegal.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ml
  integer ( kind = 4 ) n

  ierror = 0

  if ( ml < 0 .or. n - 1 < ml ) then
    ierror = ierror + 1
    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) 'R8BLT_CHECK - Illegal ML = ', ml
  end if

  if ( n <= 0 ) then
    ierror = ierror + 2
    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) 'R8BLT_CHECK - Illegal N = ', n
    return
  end if

  return
end
subroutine r8blt_print ( n, ml, a, title )

!*****************************************************************************80
!
!! R8BLT_PRINT prints a band lower triangular matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 May 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) ML, the lower bandwidth.
!    ML must be nonnegative, and no greater than N-1.
!
!    Input, real ( kind = 8 ) A(ML+1,N), the N by N band matrix, stored 
!    in band lower triangle mode.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) ml
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(ml+1,n)
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  call r8blt_print_some ( n, ml, a, 1, 1, n, n )

  return
end
subroutine r8blt_print_some ( n, ml, a, ilo, jlo, ihi, jhi )

!*****************************************************************************80
!
!! R8BLT_PRINT_SOME prints some of a band lower triangular matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 May 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) ML, the upper (and lower) bandwidth.
!    ML must be nonnegative, and no greater than N-1.
!
!    Input, real ( kind = 8 ) A(ML+1,N), the N by N band lower triangular
!    matrix.
!
!    Input, integer ( kind = 4 ) ILO, JLO, IHI, JHI, designate the first row and
!    column, and the last row and column to be printed.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) ml
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(ml+1,n)
  real ( kind = 8 ) aij
  character ( len = 14 ) ctemp(incx)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
!
!  Check the dimensions.
!
  call r8blt_check ( n, ml, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8BLT_PRINT - Fatal error!'
    write ( *, '(a)' ) '  Illegal dimensions.'
    return
  end if
!
!  Print the columns of the matrix, in strips of 5.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i7,7x)' ) j
    end do

    write ( *, '(a,5a14)' ) 'Columns', ( ctemp(j2), j2 = 1, inc )
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )
    i2lo = max ( i2lo, j2lo )

    i2hi = min ( ihi, n )
    i2hi = min ( i2hi, j2hi + ml )

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( j <= i .and. i <= j + ml ) then
          aij = a(i-j+1,j)
        else
          aij = 0.0D+00
        end if

        if ( ml < i-j .or. 0 < j-i ) then
          ctemp(j2) = '              '
        else
          write ( ctemp(j2), '(g14.6)' ) aij
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  write ( *, '(a)' ) ' '

  return
end
subroutine r8pbl_check ( n, mu, ierror )

!*****************************************************************************80
!
!! R8PBL_CHECK checks the dimensions of a positive definite symmetric band matrix.
!
!  Discussion:
!
!    To save storage, only the diagonal and lower triangle of A is stored,
!    in a compact diagonal format that preserves columns.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 May 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) MU, the upper bandwidth of the matrix.
!    MU must be at least 0, and SHOULD BE no greater than N-1.
!
!    Output, integer ( kind = 4 ) IERROR, reports whether any errors were detected.
!    IERROR is set to 0 before the checks are made, and then:
!    IERROR = IERROR + 1 if MU is illegal;
!    IERROR = IERROR + 2 if N is illegal.
!
  implicit none

  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n

  ierror = 0

  if ( mu < 0 ) then
    ierror = ierror + 1
    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) 'R8PBL_CHECK - Illegal MU < 0 = ', mu
  end if

  if ( n - 1 < mu ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) 'R8PBL_CHECK - Warning!'
    write ( *, '(a,i6)' ) '  Not advisable to have N - 1 < MU!'
    write ( *, '(a,i6)' ) '  MU = ', mu
    write ( *, '(a,i6)' ) '  N =  ', n
  end if

  if ( n <= 0 ) then
    ierror = ierror + 2
    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) 'R8PBL_CHECK - Illegal N = ', n
    return
  end if

  return
end
subroutine r8pbl_write ( n, mu, a, filename )

!*****************************************************************************80
!
!! R8PBL_WRITE writes a symmetric banded matrix to a file.
!
!  Discussion:
!
!    To save storage, only the diagonal and lower triangle of A is stored,
!    in a compact diagonal format that preserves columns.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 May 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Input, integer ( kind = 4 ) MU, the upper (and lower) bandwidth.
!    MU must be nonnegative, and no greater than N-1.
!
!    Input, real ( kind = 8 ) A(MU+1,N), the N by N band matrix, stored 
!    in positive definite symmetric band storage lower triangle mode.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 3
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(mu+1,n)
  real ( kind = 8 ) aij
  character ( len = 14 ) ctemp(incx)
  character ( len = * ) filename
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo

  open ( unit = 1, file = filename, status = 'replace' )
!
!  Determine the range of the rows in this strip.
!
  do i = 1, n

    do jlo = 1, n, 3

      jhi = min ( jlo + 2, n )

      do j = jlo, jhi

        if ( i <= j .and. j <= i + mu ) then
          aij = a(j-i+1,i)
        else if ( j <= i .and. i <= j + mu ) then
          aij = a(i-j+1,j)
        else
          aij = 0.0D+00
        end if

        write ( ctemp(j+1-jlo), '(g14.6)' ) aij

      end do

      write ( 1, '(3a14)' ) ( ctemp(j+1-jlo), j = jlo, jhi )

    end do

  end do

  close ( unit = 1 )

  return
end
function refqbf ( iq, xsi, eta )

!*****************************************************************************80
!
!! REFQBF evaluates a reference element quadratic basis function.
!
!  Discussion:
!
!    There are six possible quadratic basis functions.  This routine
!    evaluates just one of them, and its X and Y derivatives, at a
!    particular point in a particular element, by referring to the
!    reference triangle.
!
!    The point we are interested in is referred to by its coordinates
!    in the reference triangle.  That is, we are given coordinates
!    (XSI, ETA), even though, physically, we are interested
!    in points in (X, Y) space.
!
!    Here is a graph of the (XSI, ETA) reference triangle.
!
!          ^
!          |
!      1.0 +    3
!          |    |\
!      0.5 |    6 5
!          |    |  \
!      0.0 +    1-4-2
!          |
!          +----+---+--->
!               0 0 1
!               . . .
!               0 5 0
!
!                XSI
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 June 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IQ, the local node number, between 1 and
!    6, whose basis function is being evaluated.
!
!    Input, real ( kind = 8 ) XSI, the XSI coordinate of the point.
!
!    Input, real ( kind = 8 ) ETA, the ETA coordinate of the point.
!
!    Output, real ( kind = 8 ) REFQBF, the value of the basis function
!    PSI(IQ)(XSI,ETA).
!
  implicit none

  real ( kind = 8 ) eta
  integer ( kind = 4 ) iq
  real ( kind = 8 ) refqbf
  real ( kind = 8 ) w
  real ( kind = 8 ) xsi
!
!  W(1)(XSI,ETA) = 0 if XSI + ETA = 0.5 or XSI + ETA = 1.
!  W(1)(0.0,0.0) = 1.
!
  if ( iq == 1 ) then

    w = 2.0D+00 * ( 0.5D+00 - xsi - eta ) * ( 1.0D+00 - xsi - eta )
!
!  W(2)(XSI,ETA) = 0 if XSI=0 or XSI=0.5.
!  W(2)(1.0,0.0) = 1.
!
  else if ( iq == 2 ) then

    w = 2.0D+00 * xsi * ( xsi - 0.5D+00 )
!
!  W(3)(XSI,ETA) = 0 if ETA = 0, or ETA = 0.5.
!  W(3)(0.0,1.0) = 1.
!
  else if ( iq == 3 ) then

    w = 2.0D+00 * eta * ( eta - 0.5D+00 )
!
!  W(4)(XSI,ETA) = 0 if XSI = 0 or XSI + ETA = 1.
!  W(4)(0.5,0.0) = 1.
!
  else if ( iq == 4 ) then

    w = 4.0D+00 * xsi * ( 1.0D+00 - xsi - eta )
!
!  W(5)(XSI,ETA) = 0 if ETA = 0 or XSI = 0.
!  W(5)(0.5,0.5) = 1.
!
  else if ( iq == 5 ) then

    w = 4.0D+00 * eta * xsi
!
!  W(6)(XSI,ETA) = 0 if ETA = 0 or XSI + ETA = 1.
!  W(6)(0.0,0.5) = 1.
!
  else if ( iq == 6 ) then

    w = 4.0D+00 * eta * ( 1.0D+00 - xsi - eta )

  end if

  refqbf = w

  return
end
subroutine r8vec_print ( n, a, title )

!*****************************************************************************80
!
!! R8VEC_PRINT prints an R8VEC.
!
!  Discussion:
!
!    If all the entries are integers, the data if printed
!    in integer format.
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

  if ( title /= ' ' ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

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
function s_eqi ( s1, s2 )

!*****************************************************************************80
!
!! S_EQI is a case insensitive comparison of two strings for equality.  
!
!  Example:
!
!    S_EQI ( 'Anjana', 'ANJANA' ) is TRUE.
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
  integer ( kind = 4 ) lenc
  logical s_eqi
  character ( len = * ) s1
  integer ( kind = 4 ) s1_length
  character ( len = * ) s2
  integer ( kind = 4 ) s2_length

  s1_length = len ( s1 )
  s2_length = len ( s2 )
  lenc = min ( s1_length, s2_length )
 
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
 
  do i = lenc + 1, s1_length
    if ( s1(i:i) /= ' ' ) then
      return
    end if
  end do
 
  do i = lenc + 1, s2_length
    if ( s2(i:i) /= ' ' ) then
      return
    end if
  end do
 
  s_eqi = .true.
 
  return
end
subroutine s_input ( string, value, ierror )

!*****************************************************************************80
!
!! S_INPUT prints a prompt string and reads a string from the user.
!
!  Discussion:
!
!    If the input line starts with a comment character ('#'), or is blank,
!    the routine ignores that line, and tries to read the next one.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) STRING, the prompt string.
!
!    Output, character ( len = * ) VALUE, the value input by the user.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag, which is 0
!    if no error occurred.
!
  implicit none

  integer ( kind = 4 ) ierror
  character ( len = * ) string
  character ( len = * ) value

  ierror = 0
  value = ' '
!
!  Write the prompt.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( string )

  do

    read ( *, '(a)', iostat = ierror ) value

    if ( ierror /= 0 ) then
      value = 'S_INPUT: Input error!'
      return
    end if
!
!  If the line begins with a comment character, go back and read the next line.
!
    if ( value(1:1) == '#' ) then
      cycle
    end if

    if ( len_trim ( value ) == 0 ) then
      cycle
    end if

    exit

  end do

  return
end
function s_of_i4 ( i )

!*****************************************************************************80
!
!! S_OF_I4 converts an integer to a left-justified string.
!
!  Example:
!
!         I  S
!
!         1  1
!        -1  -1
!         0  0
!      1952  1952
!    123456  123456
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, an integer to be converted.
!
!    Output, character ( len = 11 ) S_OF_I4, the representation of the
!    integer ( kind = 4 ).  The integer will be left-justified.
!
  implicit none

  character c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) idig
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) ipos
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) j
  character ( len = 11 ) s
  character ( len = 11 ) s_of_i4

  s = ' '

  ilo = 1
  ihi = 11
!
!  Make a copy of the integer.
!
  ival = i
!
!  Handle the negative sign.
!
  if ( ival < 0 ) then

    if ( ihi <= 1 ) then
      s(1:1) = '*'
      return
    end if

    ival = -ival
    s(1:1) = '-'
    ilo = 2

  end if
!
!  The absolute value of the integer goes into S(ILO:IHI).
!
  ipos = ihi
!
!  Find the last digit, strip it off, and stick it into the string.
!
  do

    idig = mod ( ival, 10 )
    ival = ival / 10

    if ( ipos < ilo ) then
      do j = 1, ihi
        s(j:j) = '*'
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

  s_of_i4 = s

  return
end
subroutine s_replace_ch ( s, c1, c2 )

!*****************************************************************************80
!
!! S_REPLACE_CH replaces all occurrences of one character by another.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) S, the string.
!
!    Input, character C1, C2, the character to be replaced, and the
!    replacement character.
!
  implicit none

  character c1
  character c2
  integer ( kind = 4 ) i
  character ( len = * ) s
  integer ( kind = 4 ) s_length

  s_length = len_trim ( s )

  do i = 1, s_length
    if ( s(i:i) == c1 ) then
      s(i:i) = c2
    end if
  end do

  return
end
subroutine s_to_i4 ( s, value, ierror, length )

!*****************************************************************************80
!
!! S_TO_I4 reads an integer value from a string.
!
!  Discussion:
!
!    Instead of ICHAR, we now use the IACHAR function, which
!    guarantees the ASCII collating sequence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 January 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, a string to be examined.
!
!    Output, integer ( kind = 4 ) VALUE, the integer value read from the string.
!    If the string is blank, then VALUE will be returned 0.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!    0, no error.
!    1, an error occurred.
!
!    Output, integer ( kind = 4 ) LENGTH, the number of characters
!    of S used to make the integer.
!
  implicit none

  character c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) length
  character ( len = * ) s
  integer ( kind = 4 ) state
  integer ( kind = 4 ) value

  value = 0
  ierror = 0
  length = 0

  state = 0
  isgn = 1

  do i = 1, len_trim ( s )

    c = s(i:i)
!
!  STATE = 0, haven't read anything.
!
    if ( state == 0 ) then

      if ( c == ' ' ) then

      else if ( c == '-' ) then
        state = 1
        isgn = -1
      else if ( c == '+' ) then
        state = 1
        isgn = +1
      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        state = 2
        value = iachar ( c ) - iachar ( '0' )
      else
        ierror = 1
        return
      end if
!
!  STATE = 1, have read the sign, expecting digits or spaces.
!
    else if ( state == 1 ) then

      if ( c == ' ' ) then

      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        state = 2
        value = iachar ( c ) - iachar ( '0' )
      else
        ierror = 1
        return
      end if
!
!  STATE = 2, have read at least one digit, expecting more.
!
    else if ( state == 2 ) then

      if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then

        value = 10 * value + iachar ( c ) - iachar ( '0' )

      else

        value = isgn * value
        ierror = 0
        length = i - 1
        return

      end if

    end if

  end do
!
!  If we read all the characters in the string, see if we're OK.
!
  if ( state == 2 ) then

    value = isgn * value
    ierror = 0
    length = len_trim ( s )

  else

    value = 0
    ierror = 1
    length = 0

  end if

  return
end
subroutine s_to_r8 ( s, dval, ierror, length )

!*****************************************************************************80
!
!! S_TO_R8 reads an R8 value from a string.
!
!  Discussion:
!
!    An "R8" value is simply a real number to be stored as a
!    variable of type "real ( kind = 8 )".
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
  integer ( kind = 4 ) ndig
  real ( kind = 8 ) rbot
  real ( kind = 8 ) rexp
  real ( kind = 8 ) rtop
  character ( len = * ) s
  integer ( kind = 4 ) s_length

  s_length = len_trim ( s )

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

    if ( s_length < length+1 ) then
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
!  entire string, then we're done, and LENGTH is equal to S_LENGTH.
!
  if ( iterm /= 1 .and. length+1 == s_length ) then
    length = s_length
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
subroutine s_word_count ( s, word_num )

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
!    05 October 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string to be examined.
!
!    Output, integer ( kind = 4 ) WORD_NUM, the number of "words" in the
!    string.  Words are presumed to be separated by one or more blanks.
!
  implicit none

  logical                blank
  integer ( kind = 4 ) i
  character ( len = * )  s
  integer ( kind = 4 ) s_length
  integer ( kind = 4 ) word_num

  word_num = 0
  s_length = len ( s )

  if ( s_length <= 0 ) then
    return
  end if

  blank = .true.

  do i = 1, s_length

    if ( s(i:i) == ' ' ) then
      blank = .true.
    else if ( blank ) then
      word_num = word_num + 1
      blank = .false.
    end if

  end do

  return
end
subroutine singular_vectors ( m, n, basis_num, a, sval )

!*****************************************************************************80
!
!! SINGULAR_VECTORS computes the desired singular values.
!
!  Discussion:
!
!    The LAPACK SVD routine SGESVD is used to compute the singular
!    value decomposition:
!
!      A = U * S * V'
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the number of spatial dimensions.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of data points.
!
!    Input, integer ( kind = 4 ) BASIS_NUM, the number of basis vectors 
!    to be extracted.
!
!    Input, real POINT(DIM_NUM,POINT_NUM), the data points.
!
  implicit none

  integer ( kind = 4 ) basis_num
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) lda
  integer ( kind = 4 ) ldu
  integer ( kind = 4 ) ldvt
  character jobu
  character jobvt
  integer ( kind = 4 ) lwork
  real ( kind = 8 ) s(max(m,n))
  real ( kind = 8 ) sval(basis_num)
  real ( kind = 8 ) u(1,1)
  real ( kind = 8 ) vt(1,1)
  real ( kind = 8 ) work(3*min(m,n)+max(max(m,n),2*min(m,n)))

  lwork = 3 * min ( m, n ) + max ( max ( m, n ), 2 * min ( m, n ) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SINGULAR_VECTORS'
  write ( *, '(a)' ) '  For an MxN matrix A in general storage,'
  write ( *, '(a)' ) '  SGESVD computes the singular value decomposition:'
  write ( *, '(a)' ) '    A = U * S * V'''
  write ( *, '(a)' ) ' '
!
!  Compute the eigenvalues and eigenvectors.
!
  jobu = 'O'
  jobvt = 'N'
  lda = m
  ldu = m
  ldvt = n
  
  call dgesvd ( jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, &
    lwork, info )
 
  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SINGULAR_VECTORS - Warning:'
    write ( *, '(a,i6)' ) '  DGESVD returned nonzero INFO = ', info
    return
  end if

  sval(1:basis_num) = s(1:basis_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Singular values:'
  write ( *, '(a)' ) ' '

  do i = 1, basis_num
    write ( *, '(i4,f10.4)' ) i, sval(i)
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
  character ( len = 8 )  date
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
  character ( len = 5 )  zone

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
subroutine triangle_unit_set ( rule, order, xtab, ytab, weight )

!*****************************************************************************80
!
!! TRIANGLE_UNIT_SET sets a quadrature rule in the unit triangle.
!
!  Discussion:
!
!    The user is responsible for determining the value of ORDER,
!    and appropriately dimensioning the arrays XTAB, YTAB and
!    WEIGHT so that they can accommodate the data.
!
!    The value of ORDER for each rule can be found by invoking
!    the function TRIANGLE_RULE_SIZE.
!
!  Integration region:
!
!      0 <= X,
!    and
!      0 <= Y,
!    and
!      X + Y <= 1.
!
!  Graph:
!
!      ^
!    1 | *
!      | |\
!    Y | | \
!      | |  \
!    0 | *---*
!      +------->
!        0 X 1
!
!   The rules are accessed by an index number, RULE.  The indices,
!   and the descriptions of the corresponding rules, are:
!
!     1, ORDER =  1, precision 1, Zienkiewicz #1.
!     2, ORDER =  2, precision 1, (the "vertex rule").
!     3, ORDER =  3, precision 2, Strang and Fix formula #1.
!     4, ORDER =  3, precision 2, Strang and Fix formula #2,
!                                 Zienkiewicz #2.
!     5, ORDER =  4, precision 3, Strang and Fix formula #3,
!                                 Zienkiewicz #3.
!     6, ORDER =  6, precision 3, Strang and Fix formula #4.
!     7, ORDER =  6, precision 3, Stroud formula T2:3-1.
!     8, ORDER =  6, precision 4, Strang and Fix formula #5.
!     9, ORDER =  7, precision 4, Strang and Fix formula #6.
!    10, ORDER =  7, precision 5, Strang and Fix formula #7,
!                                 Stroud formula T2:5-1,
!                                 Zienkiewicz #4,
!                                 Schwarz Table 2.2.
!    11, ORDER =  9, precision 6, Strang and Fix formula #8.
!    12, ORDER = 12, precision 6, Strang and Fix formula #9.
!    13, ORDER = 13, precision 7, Strang and Fix formula #10.
!        Note that there is a typographical error in Strang and Fix
!        which lists the value of the XSI(3) component of the
!        last generator point as 0.4869... when it should be 0.04869...
!    14, ORDER =  7, precision 3.
!    15, ORDER = 16, precision 7, conical product Gauss, Stroud formula T2:7-1.
!    16, ORDER = 64, precision 15, triangular product Gauss rule.
!    17, ORDER = 19, precision 8, from CUBTRI, ACM TOMS #584.
!    18, ORDER = 19, precision 9, from TRIEX, ACM TOMS #612.
!    19, ORDER = 28, precision 11, from TRIEX, ACM TOMS #612.
!    20, ORDER = 37, precision 13, from ACM TOMS #706.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 April 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jarle Berntsen, Terje Espelid,
!    Algorithm 706,
!    DCUTRI: an algorithm for adaptive cubature over a collection of triangles,
!    ACM Transactions on Mathematical Software,
!    Volume 18, Number 3, September 1992, pages 329-342.
!
!    Elise deDoncker, Ian Robinson,
!    Algorithm 612:
!    Integration over a Triangle Using Nonlinear Extrapolation,
!    ACM Transactions on Mathematical Software,
!    Volume 10, Number 1, March 1984, pages 17-22.
!
!    Dirk Laurie,
!    Algorithm 584,
!    CUBTRI, Automatic Cubature Over a Triangle,
!    ACM Transactions on Mathematical Software,
!    Volume 8, Number 2, 1982, pages 210-218.
!
!    James Lyness, Dennis Jespersen,
!    Moderate Degree Symmetric Quadrature Rules for the Triangle,
!    Journal of the Institute of Mathematics and its Applications,
!    Volume 15, Number 1, February 1975, pages 19-32.
!
!    Hans Rudolf Schwarz,
!    Finite Element Methods,
!    Academic Press, 1988,
!    ISBN: 0126330107,
!    LC: TA347.F5.S3313.
!
!    Gilbert Strang, George Fix,
!    An Analysis of the Finite Element Method,
!    Cambridge, 1973,
!    ISBN: 096140888X,
!    LC: TA335.S77.
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!    Olgierd Zienkiewicz,
!    The Finite Element Method,
!    Sixth Edition,
!    Butterworth-Heinemann, 2005,
!    ISBN: 0750663200,
!    LC: TA640.2.Z54
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) RULE, the index of the rule.
!
!    Input, integer ( kind = 4 ) ORDER, the order of the rule.
!
!    Output, real ( kind = 8 ) XTAB(ORDER), YTAB(ORDER), the abscissas.
!
!    Output, real ( kind = 8 ) WEIGHT(ORDER), the weights of the rule.
!
  implicit none

  integer ( kind = 4 ) order

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) d
  real ( kind = 8 ) e
  real ( kind = 8 ) f
  real ( kind = 8 ) g
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) order2
  real ( kind = 8 ) p
  real ( kind = 8 ) q
  real ( kind = 8 ) r
  integer ( kind = 4 ) rule
  real ( kind = 8 ) s
  real ( kind = 8 ) t
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) w
  real ( kind = 8 ) w1
  real ( kind = 8 ) w2
  real ( kind = 8 ) w3
  real ( kind = 8 ) w4
  real ( kind = 8 ) w5
  real ( kind = 8 ) w6
  real ( kind = 8 ) w7
  real ( kind = 8 ) w8
  real ( kind = 8 ) w9
  real ( kind = 8 ) weight(order)
  real ( kind = 8 ) weight1(8)
  real ( kind = 8 ) weight2(8)
  real ( kind = 8 ) wx
  real ( kind = 8 ) x
  real ( kind = 8 ) xtab(order)
  real ( kind = 8 ) xtab1(8)
  real ( kind = 8 ) xtab2(8)
  real ( kind = 8 ) y
  real ( kind = 8 ) ytab(order)
  real ( kind = 8 ) z
!
!  1 point, precision 1.
!
  if ( rule == 1 ) then

    xtab(1)   = 0.33333333333333333333D+00

    ytab(1)   = 0.33333333333333333333D+00

    weight(1) = 1.00000000000000000000D+00
!
!  3 points, precision 1, the "vertex rule".
!
  else if ( rule == 2 ) then

    xtab(1) =   1.00000000000000000000D+00
    xtab(2) =   0.00000000000000000000D+00
    xtab(3) =   0.00000000000000000000D+00

    ytab(1) =   0.00000000000000000000D+00
    ytab(2) =   1.00000000000000000000D+00
    ytab(3) =   0.00000000000000000000D+00

    weight(1) = 0.33333333333333333333D+00
    weight(2) = 0.33333333333333333333D+00
    weight(3) = 0.33333333333333333333D+00
!
!  3 points, precision 2, Strang and Fix formula #1.
!
  else if ( rule == 3 ) then

    xtab(1)   = 0.66666666666666666667D+00
    xtab(2)   = 0.16666666666666666667D+00
    xtab(3)   = 0.16666666666666666667D+00

    ytab(1)   = 0.16666666666666666667D+00
    ytab(2)   = 0.66666666666666666667D+00
    ytab(3)   = 0.16666666666666666667D+00

    weight(1) = 0.33333333333333333333D+00
    weight(2) = 0.33333333333333333333D+00
    weight(3) = 0.33333333333333333333D+00
!
!  3 points, precision 2, Strang and Fix formula #2.
!
  else if ( rule == 4 ) then

    xtab(1)   = 0.50000000000000000000D+00
    xtab(2)   = 0.50000000000000000000D+00
    xtab(3)   = 0.00000000000000000000D+00

    ytab(1)   = 0.00000000000000000000D+00
    ytab(2)   = 0.50000000000000000000D+00
    ytab(3)   = 0.50000000000000000000D+00

    weight(1) = 0.33333333333333333333D+00
    weight(2) = 0.33333333333333333333D+00
    weight(3) = 0.33333333333333333333D+00
!
!  4 points, precision 3, Strang and Fix formula #3.
!
  else if ( rule == 5 ) then

    a =   6.0D+00
    b =  10.0D+00
    c =  18.0D+00
    d =  25.0D+00
    e = -27.0D+00
    f =  30.0D+00
    g =  48.0D+00

    xtab(1:4) =   (/ b, c, a, a /) / f
    ytab(1:4) =   (/ b, a, c, a /) / f
    weight(1:4) = (/ e, d, d, d /) / g
!
!  6 points, precision 3, Strang and Fix formula #4.
!
  else if ( rule == 6 ) then

    a = 0.659027622374092D+00
    b = 0.231933368553031D+00
    c = 0.109039009072877D+00

    xtab(1:6) =   (/ a, a, b, b, c, c /)
    ytab(1:6) =   (/ b, c, a, c, a, b /)

    weight(1) = 0.16666666666666666667D+00
    weight(2) = 0.16666666666666666667D+00
    weight(3) = 0.16666666666666666667D+00
    weight(4) = 0.16666666666666666667D+00
    weight(5) = 0.16666666666666666667D+00
    weight(6) = 0.16666666666666666667D+00
!
!  6 points, precision 3, Stroud T2:3-1.
!
  else if ( rule == 7 ) then

    a = 0.0D+00
    b = 0.5D+00
    c = 2.0D+00 /  3.0D+00
    d = 1.0D+00 /  6.0D+00
    v = 1.0D+00 / 30.0D+00
    w = 3.0D+00 / 10.0D+00

    xtab(1:6) =   (/ a, b, b, c, d, d /)
    ytab(1:6) =   (/ b, a, b, d, c, d /)
    weight(1:6) = (/ v, v, v, w, w, w /)
!
!  6 points, precision 4, Strang and Fix, formula #5.
!
  else if ( rule == 8 ) then

    a = 0.816847572980459D+00
    b = 0.091576213509771D+00
    c = 0.108103018168070D+00
    d = 0.445948490915965D+00
    v = 0.109951743655322D+00
    w = 0.223381589678011D+00

    xtab(1:6) =   (/ a, b, b, c, d, d /)
    ytab(1:6) =   (/ b, a, b, d, c, d /)
    weight(1:6) = (/ v, v, v, w, w, w /)
!
!  7 points, precision 4, Strang and Fix formula #6.
!
  else if ( rule == 9 ) then

    a = 1.0D+00 / 3.0D+00
    c = 0.736712498968435D+00
    d = 0.237932366472434D+00
    e = 0.025355134551932D+00
    v = 0.375000000000000D+00
    w = 0.104166666666667D+00

    xtab(1:7) =   (/ a, c, c, d, d, e, e /)
    ytab(1:7) =   (/ a, d, e, c, e, c, d /)
    weight(1:7) = (/ v, w, w, w, w, w, w /)
!
!  7 points, precision 5, Strang and Fix formula #7, Stroud T2:5-1
!
  else if ( rule == 10 ) then

    a = 1.0D+00 / 3.0D+00
    b = ( 9.0D+00 + 2.0D+00 * sqrt ( 15.0D+00 ) ) / 21.0D+00
    c = ( 6.0D+00 -           sqrt ( 15.0D+00 ) ) / 21.0D+00
    d = ( 9.0D+00 - 2.0D+00 * sqrt ( 15.0D+00 ) ) / 21.0D+00
    e = ( 6.0D+00 +           sqrt ( 15.0D+00 ) ) / 21.0D+00
    u = 0.225D+00
    v = ( 155.0D+00 - sqrt ( 15.0D+00 ) ) / 1200.0D+00
    w = ( 155.0D+00 + sqrt ( 15.0D+00 ) ) / 1200.0D+00

    xtab(1:7) =   (/ a, b, c, c, d, e, e /)
    ytab(1:7) =   (/ a, c, b, c, e, d, e /)
    weight(1:7) = (/ u, v, v, v, w, w, w /)
!
!  9 points, precision 6, Strang and Fix formula #8.
!
  else if ( rule == 11 ) then

    a = 0.124949503233232D+00
    b = 0.437525248383384D+00
    c = 0.797112651860071D+00
    d = 0.165409927389841D+00
    e = 0.037477420750088D+00

    u = 0.205950504760887D+00
    v = 0.063691414286223D+00

    xtab(1:9) =   (/ a, b, b, c, c, d, d, e, e /)
    ytab(1:9) =   (/ b, a, b, d, e, c, e, c, d /)
    weight(1:9) = (/ u, u, u, v, v, v, v, v, v /)
!
!  12 points, precision 6, Strang and Fix, formula #9.
!
  else if ( rule == 12 ) then

    a = 0.873821971016996D+00
    b = 0.063089014491502D+00
    c = 0.501426509658179D+00
    d = 0.249286745170910D+00
    e = 0.636502499121399D+00
    f = 0.310352451033785D+00
    g = 0.053145049844816D+00

    u = 0.050844906370207D+00
    v = 0.116786275726379D+00
    w = 0.082851075618374D+00

    xtab(1:12) =   (/ a, b, b, c, d, d, e, e, f, f, g, g /)
    ytab(1:12) =   (/ b, a, b, d, c, d, f, g, e, g, e, f /)
    weight(1:12) = (/ u, u, u, v, v, v, w, w, w, w, w, w /)
!
!  13 points, precision 7, Strang and Fix, formula #10.
!
!  Note that there is a typographical error in Strang and Fix
!  which lists the value of the XSI(3) component of the
!  last generator point as 0.4869... when it should be 0.04869...
!
  else if ( rule == 13 ) then

    h = 1.0D+00 / 3.0D+00
    a = 0.479308067841923D+00
    b = 0.260345966079038D+00
    c = 0.869739794195568D+00
    d = 0.065130102902216D+00
    e = 0.638444188569809D+00
    f = 0.312865496004875D+00
    g = 0.048690315425316D+00

    w = -0.149570044467670D+00
    t =  0.175615257433204D+00
    u =  0.053347235608839D+00
    v =  0.077113760890257D+00

    xtab(1:13) =   (/ h, a, b, b, c, d, d, e, e, f, f, g, g /)
    ytab(1:13) =   (/ h, b, a, b, d, c, d, f, g, e, g, e, f /)
    weight(1:13) = (/ w, t, t, t, u, u, u, v, v, v, v, v, v /)
!
!  7 points, precision 3.
!
  else if ( rule == 14 ) then

    a = 1.0D+00 / 3.0D+00
    b = 1.0D+00
    c = 0.5D+00
    z = 0.0D+00

    u = 27.0D+00 / 60.0D+00
    v =  3.0D+00 / 60.0D+00
    w =  8.0D+00 / 60.0D+00

    xtab(1:7) =   (/ a, b, z, z, z, c, c /)
    ytab(1:7) =   (/ a, z, b, z, c, z, c /)
    weight(1:7) = (/ u, v, v, v, w, w, w /)
!
!  16 points, precision 5, Stroud T2:7-1.
!
  else if ( rule == 15 ) then
!
!  Legendre rule of order 4.
!
    order2 = 4

    xtab1(1:4) = (/ &
      -0.861136311594052575223946488893D+00, &
      -0.339981043584856264802665759103D+00, &
       0.339981043584856264802665759103D+00, &
       0.861136311594052575223946488893D+00 /)

    weight1(1:4) = (/ &
      0.347854845137453857373063949222D+00, &
      0.652145154862546142626936050778D+00, &
      0.652145154862546142626936050778D+00, &
      0.347854845137453857373063949222D+00 /)

    xtab1(1:order2) = 0.5D+00 * ( xtab1(1:order2) + 1.0D+00 )

    weight2(1) = 0.1355069134D+00
    weight2(2) = 0.2034645680D+00
    weight2(3) = 0.1298475476D+00
    weight2(4) = 0.0311809709D+00

    xtab2(1) = 0.0571041961D+00
    xtab2(2) = 0.2768430136D+00
    xtab2(3) = 0.5835904324D+00
    xtab2(4) = 0.8602401357D+00

    k = 0
    do i = 1, order2
      do j = 1, order2
        k = k + 1
        xtab(k) = xtab2(j)
        ytab(k) = xtab1(i) * ( 1.0D+00 - xtab2(j) )
        weight(k) = weight1(i) * weight2(j)
      end do
    end do
!
!  64 points, precision 15.
!
  else if ( rule == 16 ) then
!
!  Legendre rule of order 8.
!
    order2 = 8

    xtab1(1) = -0.960289856497536231683560868569D+00
    xtab1(2) = -0.796666477413626739591553936476D+00
    xtab1(3) = -0.525532409916328985817739049189D+00
    xtab1(4) = -0.183434642495649804939476142360D+00
    xtab1(5) =  0.183434642495649804939476142360D+00
    xtab1(6) =  0.525532409916328985817739049189D+00
    xtab1(7) =  0.796666477413626739591553936476D+00
    xtab1(8) =  0.960289856497536231683560868569D+00

    weight1(1) = 0.101228536290376259152531354310D+00
    weight1(2) = 0.222381034453374470544355994426D+00
    weight1(3) = 0.313706645877887287337962201987D+00
    weight1(4) = 0.362683783378361982965150449277D+00
    weight1(5) = 0.362683783378361982965150449277D+00
    weight1(6) = 0.313706645877887287337962201987D+00
    weight1(7) = 0.222381034453374470544355994426D+00
    weight1(8) = 0.101228536290376259152531354310D+00

    weight2(1) = 0.00329519144D+00
    weight2(2) = 0.01784290266D+00
    weight2(3) = 0.04543931950D+00
    weight2(4) = 0.07919959949D+00
    weight2(5) = 0.10604735944D+00
    weight2(6) = 0.11250579947D+00
    weight2(7) = 0.09111902364D+00
    weight2(8) = 0.04455080436D+00

    xtab2(1) = 0.04463395529D+00
    xtab2(2) = 0.14436625704D+00
    xtab2(3) = 0.28682475714D+00
    xtab2(4) = 0.45481331520D+00
    xtab2(5) = 0.62806783542D+00
    xtab2(6) = 0.78569152060D+00
    xtab2(7) = 0.90867639210D+00
    xtab2(8) = 0.98222008485D+00

    k = 0
    do j = 1, order2
      do i = 1, order2
        k = k + 1
        xtab(k) = 1.0D+00 - xtab2(j)
        ytab(k) = 0.5D+00 * ( 1.0D+00 + xtab1(i) ) * xtab2(j)
        weight(k) = weight1(i) * weight2(j)
      end do
    end do
!
!  19 points, precision 8, from CUBTRI.
!
  else if ( rule == 17 ) then

    a = 1.0D+00 / 3.0D+00
    b = ( 9.0D+00 + 2.0D+00 * sqrt ( 15.0D+00 ) ) / 21.0D+00
    c = ( 6.0D+00 -       sqrt ( 15.0D+00 ) ) / 21.0D+00
    d = ( 9.0D+00 - 2.0D+00 * sqrt ( 15.0D+00 ) ) / 21.0D+00
    e = ( 6.0D+00 +       sqrt ( 15.0D+00 ) ) / 21.0D+00
    f = ( 40.0D+00 - 10.0D+00 * sqrt ( 15.0D+00 ) &
      + 10.0D+00 * sqrt ( 7.0D+00 ) + 2.0D+00 * sqrt ( 105.0D+00 ) ) / 90.0D+00
    g = ( 25.0D+00 +  5.0D+00 * sqrt ( 15.0D+00 ) &
      -  5.0D+00 * sqrt ( 7.0D+00 ) - sqrt ( 105.0D+00 ) ) / 90.0D+00
    p = ( 40.0D+00 + 10.0D+00 * sqrt ( 15.0D+00 ) &
      + 10.0D+00 * sqrt ( 7.0D+00 ) - 2.0D+00 * sqrt ( 105.0D+00 ) ) / 90.0D+00
    q = ( 25.0D+00 -  5.0D+00 * sqrt ( 15.0D+00 ) &
      -  5.0D+00 * sqrt ( 7.0D+00 ) + sqrt ( 105.0D+00 ) ) / 90.0D+00
    r = ( 40.0D+00 + 10.0D+00 * sqrt ( 7.0D+00 ) ) / 90.0D+00
    s = ( 25.0D+00 +  5.0D+00 * sqrt ( 15.0D+00 ) - 5.0D+00 * sqrt ( 7.0D+00 ) &
      - sqrt ( 105.0D+00 ) ) / 90.0D+00
    t = ( 25.0D+00 -  5.0D+00 * sqrt ( 15.0D+00 ) - 5.0D+00 * sqrt ( 7.0D+00 ) &
      + sqrt ( 105.0D+00 ) ) / 90.0D+00

    w1 = ( 7137.0D+00 - 1800.0D+00 * sqrt ( 7.0D+00 ) ) / 62720.0D+00
    w2 = -9301697.0D+00 / 4695040.0D+00 - 13517313.0D+00 * sqrt ( 15.0D+00 ) &
      / 23475200.0D+00 + 764885.0D+00 * sqrt ( 7.0D+00 ) / 939008.0D+00 &
      + 198763.0D+00 * sqrt ( 105.0D+00 ) / 939008.0D+00
    w2 = w2 / 3.0D+00
    w3 = -9301697.0D+00 / 4695040.0D+00 + 13517313.0D+00 * sqrt ( 15.0D+00 ) &
      / 23475200.0D+00 &
      + 764885.0D+00 * sqrt ( 7.0D+00 ) / 939008.0D+00 &
      - 198763.0D+00 * sqrt ( 105.0D+00 ) / 939008.0D+00
    w3 = w3 / 3.0D+00
    w4 = ( 102791225.0D+00 - 23876225.0D+00 * sqrt ( 15.0D+00 ) &
      - 34500875.0D+00 * sqrt ( 7.0D+00 ) &
      + 9914825.0D+00 * sqrt ( 105.0D+00 ) ) / 59157504.0D+00
    w4 = w4 / 3.0D+00
    w5 = ( 102791225.0D+00 + 23876225.0D+00 * sqrt ( 15.0D+00 ) &
      - 34500875.0D+00 * sqrt ( 7.0D+00 ) &
      - 9914825D+00 * sqrt ( 105.0D+00 ) ) / 59157504.0D+00
    w5 = w5 / 3.0D+00
    w6 = ( 11075.0D+00 - 3500.0D+00 * sqrt ( 7.0D+00 ) ) / 8064.0D+00
    w6 = w6 / 6.0D+00

    xtab(1:19) =   (/  a,  b,  c,  c,  d,  e,  e,  f,  g,  g,  p,  q,  q, &
                       r,  r,  s,  s,  t,  t /)
    ytab(1:19) =   (/  a,  c,  b,  c,  e,  d,  e,  g,  f,  g,  q,  p,  q, &
                       s,  t,  r,  t,  r,  s /)
    weight(1:19) = (/ w1, w2, w2, w2, w3, w3, w3, w4, w4, w4, w5, w5, w5, &
      w6, w6, w6, w6, w6, w6 /)
!
!  19 points, precision 9.
!  Lyness and Jesperson.
!
  else if ( rule == 18 ) then

    a = 1.0D+00 / 3.0D+00
    b =  0.02063496160252593D+00
    c =  0.4896825191987370D+00
    d =  0.1258208170141290D+00
    e =  0.4370895914929355D+00
    f =  0.6235929287619356D+00
    g =  0.1882035356190322D+00
    r =  0.9105409732110941D+00
    s =  0.04472951339445297D+00
    t =  0.7411985987844980D+00
    u =  0.03683841205473626D+00
    v =  0.22196288916076574D+00

    w1 = 0.09713579628279610D+00
    w2 = 0.03133470022713983D+00
    w3 = 0.07782754100477543D+00
    w4 = 0.07964773892720910D+00
    w5 = 0.02557767565869810D+00
    w6 = 0.04328353937728940D+00

    xtab(1:19) =   (/  a,  b,  c,  c,  d,  e,  e,  f,  g,  g,  r,  s,  s, &
      t, t, u, u, v, v /)
    ytab(1:19) =   (/  a,  c,  b,  c,  e,  d,  e,  g,  f,  g,  s,  r,  s, &
      u, v, t, v, t, u /)
    weight(1:19) = (/ w1, w2, w2, w2, w3, w3, w3, w4, w4, w4, w5, w5, w5, &
      w6, w6, w6, w6, w6, w6 /)
!
!  28 points, precision 11.
!  Lyness and Jesperson.
!
  else if ( rule == 19 ) then

    a = 1.0D+00 / 3.0D+00
    b = 0.9480217181434233D+00
    c = 0.02598914092828833D+00
    d = 0.8114249947041546D+00
    e = 0.09428750264792270D+00
    f = 0.01072644996557060D+00
    g = 0.4946367750172147D+00
    p = 0.5853132347709715D+00
    q = 0.2073433826145142D+00
    r = 0.1221843885990187D+00
    s = 0.4389078057004907D+00
    t = 0.6779376548825902D+00
    u = 0.04484167758913055D+00
    v = 0.27722066752827925D+00
    w = 0.8588702812826364D+00
    x = 0.0D+00
    y = 0.1411297187173636D+00

    w1 = 0.08797730116222190D+00
    w2 = 0.008744311553736190D+00
    w3 = 0.03808157199393533D+00
    w4 = 0.01885544805613125D+00
    w5 = 0.07215969754474100D+00
    w6 = 0.06932913870553720D+00
    w7 = 0.04105631542928860D+00
    w8 = 0.007362383783300573D+00

    xtab(1:28) =   (/  a,  b,  c,  c,  d,  e,  e,  f,  g,  g,  p,  q,  q, &
       r,  s,  s,  t,  t,  u,  u,  v,  v,  w,  w,  x,  x,  y,  y /)
    ytab(1:28) =   (/  a,  c,  b,  c,  e,  d,  e,  g,  f,  g,  q,  p,  q, &
       s,  r,  s,  u,  v,  t,  v,  t,  u,  x,  y,  w,  y,  w,  x /)
    weight(1:28) = (/ w1, w2, w2, w2, w3, w3, w3, w4, w4, w4, w5, w5, w5, &
      w6, w6, w6, w7, w7, w7, w7, w7, w7, w8, w8, w8, w8, w8, w8 /)
!
!  37 points, precision 13.
!
  else if ( rule == 20 ) then

    a = 1.0D+00 / 3.0D+00
    b = 0.950275662924105565450352089520D+00
    c = 0.024862168537947217274823955239D+00
    d = 0.171614914923835347556304795551D+00
    e = 0.414192542538082326221847602214D+00
    f = 0.539412243677190440263092985511D+00
    g = 0.230293878161404779868453507244D+00

    w1 = 0.051739766065744133555179145422D+00
    w2 = 0.008007799555564801597804123460D+00
    w3 = 0.046868898981821644823226732071D+00
    w4 = 0.046590940183976487960361770070D+00
    w5 = 0.031016943313796381407646220131D+00
    w6 = 0.010791612736631273623178240136D+00
    w7 = 0.032195534242431618819414482205D+00
    w8 = 0.015445834210701583817692900053D+00
    w9 = 0.017822989923178661888748319485D+00
    wx = 0.037038683681384627918546472190D+00

    xtab(1:10) =   (/ a, b, c, c, d, e, e, f, g, g /)
    ytab(1:10) =   (/ a, c, b, c, e, d, e, g, f, g /)
    weight(1:37) = (/ w1, w2, w2, w2, w3, w3, w3, w4, w4, w4, w5, w5, w5, &
                      w6, w6, w6, w7, w7, w7, w8, w8, w8, w8, w8, w8, w9, &
                      w9, w9, w9, w9, w9, wx, wx, wx, wx, wx, wx /)

    a = 0.772160036676532561750285570113D+00
    b = 0.113919981661733719124857214943D+00

    xtab(11) = a
    ytab(11) = b

    xtab(12) = b
    ytab(12) = a

    xtab(13) = b
    ytab(13) = b

    a = 0.009085399949835353883572964740D+00
    b = 0.495457300025082323058213517632D+00

    xtab(14) = a
    ytab(14) = b

    xtab(15) = b
    ytab(15) = a

    xtab(16) = b
    ytab(16) = b

    a = 0.062277290305886993497083640527D+00
    b = 0.468861354847056503251458179727D+00

    xtab(17) = a
    ytab(17) = b

    xtab(18) = b
    ytab(18) = a

    xtab(19) = b
    ytab(19) = b

    a = 0.022076289653624405142446876931D+00
    b = 0.851306504174348550389457672223D+00
    c = 0.126617206172027096933163647918263D+00

    xtab(20) = a
    ytab(20) = b

    xtab(21) = a
    ytab(21) = c

    xtab(22) = b
    ytab(22) = a

    xtab(23) = b
    ytab(23) = c

    xtab(24) = c
    ytab(24) = a

    xtab(25) = c
    ytab(25) = b

    a = 0.018620522802520968955913511549D+00
    b = 0.689441970728591295496647976487D+00
    c = 0.291937506468887771754472382212953D+00

    xtab(26) = a
    ytab(26) = b

    xtab(27) = a
    ytab(27) = c

    xtab(28) = b
    ytab(28) = a

    xtab(29) = b
    ytab(29) = c

    xtab(30) = c
    ytab(30) = a

    xtab(31) = c
    ytab(31) = b

    a = 0.096506481292159228736516560903D+00
    b = 0.635867859433872768286976979827D+00
    c = 0.267625659273967961282458816185681D+00

    xtab(32) = a
    ytab(32) = b

    xtab(33) = a
    ytab(33) = c

    xtab(34) = b
    ytab(34) = a

    xtab(35) = b
    ytab(35) = c

    xtab(36) = c
    ytab(36) = a

    xtab(37) = c
    ytab(37) = b

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRIANGLE_UNIT_SET - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of RULE = ', rule
    stop

  end if

  return
end
