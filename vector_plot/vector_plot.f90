program main

!*****************************************************************************80
!
!! MAIN is the main program for VECTOR_PLOT.
!
!  Discussion:
!
!    VECTOR_PLOT displays velocity and direction vector plots.
!
!    Two input files are required:
!
!    * an XY file containing point positions;
!
!    * a UV file containing the velocity component values.
!
!    Both files consist of N rows of pairs of values.
!
!    The program requires the user to specify a "thinning" factor,
!    which describes the amount of data that will be thrown away.
!    This is critical, because most data files have far too many
!    vectors to display properly.  Of course, thinning is a brutal
!    and coarse method, but it works.  Originally, I simply saved
!    every THIN_NUM-th vector.  However, for the T-Cell in particular,
!    this made the rows of data suffer a noticeable "jog" in the
!    transition between the sides and the central well.  Now I try
!    to bin the X and Y coordinate values, and save all those for
!    which the bin index, mod NTHIN, is NTHIN/2.  A little obscure
!    and more prone to error, but produces a more natural looking
!    selection of data for the T-Cell.
!
!    The program writes two files:
!
!    * a VEC file containing velocity vectors;
!
!    * a DIR file containing velocity direction vectors.
!
!  Usage:
!
!    vector_plot xy_file uv_file
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 July 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) arg_num
  logical, parameter :: debug = .false.
  character ( len = 80 ) dir_file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iarg
  integer ( kind = 4 ) iargc
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) node_num
  character ( len = 11 ) s_of_i4
  integer ( kind = 4 ), allocatable, dimension ( : ) :: thin_dex
  integer ( kind = 4 ) :: thin_factor = 1
  integer ( kind = 4 ) thin_num
  real ( kind = 8 ), allocatable, dimension ( : ) :: u
  character ( len = 80 ) uv_file_name
  real ( kind = 8 ) uv_max
  real ( kind = 8 ), allocatable, dimension ( : ) :: v
  character ( len = 80 ) vec_file_name
  real ( kind = 8 ) velocity_scale
  real ( kind = 8 ), allocatable, dimension ( : ) :: x
  real ( kind = 8 ) x_max
  real ( kind = 8 ) x_min
  integer ( kind = 4 ) xy_dim
  character ( len = 80 ) xy_file_name
  integer ( kind = 4 ) xy_lines
  real ( kind = 8 ), allocatable, dimension ( : ) :: y
  real ( kind = 8 ) y_max
  real ( kind = 8 ) y_min

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'VECTOR_PLOT'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Display a vector field ( U(X,Y), V(X,Y) )'
!
!  Get the number of command line arguments.
!
  arg_num = iargc ( )
!
!  If at least one command line argument, it's the XY file name.
!
  if ( 1 <= arg_num ) then

    iarg = 1
    call getarg ( iarg, xy_file_name )

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'VECTOR_PLOT:'
    write ( *, '(a)' ) '  Please enter the name of the XY file.'

    read ( *, '(a)' ) xy_file_name

  end if
!
!  If at least two command line arguments, the second is the UV file name.
!
  if ( 2 <= arg_num ) then

    iarg = 2
    call getarg ( iarg, uv_file_name )

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'VECTOR_PLOT:'
    write ( *, '(a)' ) '  Please enter the name of the UV file.'

    read ( *, '(a)' ) uv_file_name

  end if
!
!  Examine the XY file.
!
  call data_size ( xy_file_name, xy_lines, xy_dim, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'VECTOR_PLOT - Error!'
    write ( *, '(a)' ) '  The XY data file could not be read.'
    write ( *, '(a)' ) '  Possibly the name was incorrect.'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'VECTOR_PLOT'
    write ( *, '(a)' ) '  Abnormal end of execution.'
    write ( *, '(a)' ) ' '
    call timestamp ( )
  end if

  if ( xy_dim /= 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'VECTOR_PLOT - Fatal error!'
    write ( *, '(a)' ) '  The XY data file does not contain 2D data.'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'VECTOR_PLOT'
    write ( *, '(a)' ) '  Abnormal end of execution.'
    write ( *, '(a)' ) ' '
    call timestamp ( )
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The file "' // trim ( xy_file_name ) // '" contains ' &
    // trim ( s_of_i4 ( xy_lines ) ) // ' pairs of data values.'
!
!  Read in X and Y values from the XY data file.
!
  node_num = xy_lines

  allocate ( x(1:node_num) )
  allocate ( y(1:node_num) )
  allocate ( thin_dex(1:node_num) )

  call data_d2_read ( xy_file_name, node_num, x, y, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'VECTOR_PLOT - Fatal error!'
    write ( *, '(a)' ) '  The XY data file could not be read.'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'VECTOR_PLOT'
    write ( *, '(a)' ) '  Abnormal end of execution.'
    write ( *, '(a)' ) ' '
    call timestamp ( )
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
!  Examine the UV file.
!
  allocate ( u(1:node_num) )
  allocate ( v(1:node_num) )
!
!  Read in U, V.
!
  call data_d2_read ( uv_file_name, node_num, u, v, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'VECTOR_PLOT - Error!'
    write ( *, '(a)' ) '  The velocity data file could not be read.'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'VECTOR_PLOT'
    write ( *, '(a)' ) '  Abnormal end of execution.'
    write ( *, '(a)' ) ' '
    call timestamp ( )
  end if

  uv_max = maxval ( sqrt ( u(1:node_num)**2 + v(1:node_num)**2 ) )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  The maximum velocity magnitude = ' , uv_max
!
!  Get the scale factor.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Enter VELOCITY_SCALE, the scale factor.'
  write ( *, '(a)' ) '  1.0 = "normal" length,'
  write ( *, '(a)' ) '  0.5 = half length,'
  write ( *, '(a)' ) '  2.0 = double length.'

  read ( *, * ) velocity_scale
!
!  Get the thinning increment.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Enter THIN_FACTOR, the thinning factor.'
  write ( *, '(a)' ) '  1 = use all the data,'
  write ( *, '(a)' ) '  2 = use 1/2 the data,'
  write ( *, '(a)' ) '  10 = use 1/10 the data, and so on.'

  read ( *, * ) thin_factor
!
!  Set up thinned X, Y, U, V data.
!
  call thin_index ( node_num, x, y, thin_factor, thin_num, thin_dex )

  x(1:thin_num) = x(thin_dex(1:thin_num ) )
  y(1:thin_num) = y(thin_dex(1:thin_num ) )
  u(1:thin_num) = u(thin_dex(1:thin_num ) )
  v(1:thin_num) = v(thin_dex(1:thin_num ) )

  write ( *, * ) ' '
  write ( *, * ) '  Number of nodes before thinning: ', node_num
  write ( *, * ) '  Number of nodes after  thinning: ', thin_num

  node_num = thin_num
!
!  DEBUG:
!
  if ( debug ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '       I      X         Y           U              V'
    write ( *, '(a)' ) ' '
    do i = 1, thin_num
      write ( *, '(2x,i6,2x,f8.4,2x,f8.4,2x,g14.6,2x,g14.6)' ) &
        i, x(i), y(i), u(i), v(i)
    end do
  end if
!
!  Create the vector plot.
!
  vec_file_name = uv_file_name

  call file_name_append ( vec_file_name, '_vec' )
  call file_name_ext_swap ( vec_file_name, 'eps' )

  call vector_field ( node_num, x, y, u, v, vec_file_name, velocity_scale )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'VECTOR_PLOT:'
  write ( *, '(a)' ) '  Created the output vector plot file "' &
    // trim ( vec_file_name ) // '".'
!
!  Create the direction plot.
!
  dir_file_name = uv_file_name

  call file_name_append ( dir_file_name, '_dir' )
  call file_name_ext_swap ( dir_file_name, 'eps' )

  velocity_scale = 1.0D+00

  call direction_field ( node_num, x, y, u, v, dir_file_name, velocity_scale )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'VECTOR_PLOT:'
  write ( *, '(a)' ) '  Created the output vector plot file "' &
    // trim ( dir_file_name ) // '".'
!
!  Free memory.
!
  deallocate ( thin_dex )
  deallocate ( u )
  deallocate ( v )
  deallocate ( x )
  deallocate ( y )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'VECTOR_PLOT'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine box_segment_clip_2d ( p1, p2, pa, pb, ival )

!*****************************************************************************80
!
!! BOX_SEGMENT_CLIP_2D uses a box to clip a line segment in 2D.
!
!  Discussion:
!
!    A box in 2D is a rectangle with sides aligned on coordinate
!    axes.  It can be described by its low and high corners, P1 and P2
!    as the set of points P satisfying:
!
!      P1(1:2) <= P(1:2) <= P2(1:2).
!
!    A line segment is the finite portion of a line that lies between
!    two points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 June 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) P1(2), P2(2), the low and high corners of the box.
!
!    Input/output, real ( kind = 8 ) PA(2), PB(2); on input, the endpoints
!    of a line segment.  On output, the endpoints of the portion of the
!    line segment that lies inside the box.  However, if no part of the
!    initial line segment lies inside the box, the output value is the
!    same as the input value.
!
!    Output, integer ( kind = 4 ) IVAL:
!    -1, no part of the line segment is within the box.
!     0, no clipping was necessary.
!     1, PA was clipped.
!     2, PB was clipped.
!     3, PA and PB were clipped.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  logical clip_a
  logical clip_b
  integer ( kind = 4 ) ival
  real ( kind = 8 ) p1(dim_num)
  real ( kind = 8 ) p2(dim_num)
  real ( kind = 8 ) pa(dim_num)
  real ( kind = 8 ) pb(dim_num)
  real ( kind = 8 ) q(dim_num)

  clip_a = .false.
  clip_b = .false.
!
!  Require that XMIN <= X.
!
  if ( pa(1) < p1(1) .and. pb(1) < p1(1) ) then
    ival = -1
    return
  end if

  if ( pa(1) < p1(1) .and. p1(1) <= pb(1) ) then
    q(1) = p1(1)
    q(2) = pa(2) + ( pb(2) - pa(2) ) * ( q(1) - pa(1) ) / ( pb(1) - pa(1) )
    pa(1:2) = q(1:2)
    clip_a = .true.
  else if ( p1(1) <= pa(1) .and. pb(1) < p1(1) ) then
    q(1) = p1(1)
    q(2) = pa(2) + ( pb(2) - pa(2) ) * ( q(1) - pa(1) ) / ( pb(1) - pa(1) )
    pb(1:2) = q(1:2)
    clip_b = .true.
  end if
!
!  Require that X <= XMAX.
!
  if ( p2(1) < pa(1) .and. p2(1) < pb(1) ) then
    ival = -1
    return
  end if

  if ( p2(1) < pa(1) .and. pb(1) <= p2(1) ) then
    q(1) = p2(1)
    q(2) = pa(2) + ( pb(2) - pa(2) ) * ( q(1) - pa(1) ) / ( pb(1) - pa(1) )
    pa(1:2) = q(1:2)
    clip_a = .true.
  else if ( pa(1) <= p2(1) .and. p2(1) < pb(1) ) then
    q(1) = p2(1)
    q(2) = pa(2) + ( pb(2) - pa(2) ) * ( q(1) - pa(1) ) / ( pb(1) - pa(1) )
    pb(1:2) = q(1:2)
    clip_b = .true.
  end if
!
!  Require that YMIN <= Y.
!
  if ( pa(2) < p1(2) .and. pb(2) < p1(2) ) then
    ival = -1
    return
  end if

  if ( pa(2) < p1(2) .and. p1(2) <= pb(2) ) then
    q(2) = p1(2)
    q(1) = pa(1) + ( pb(1) - pa(1) ) * ( q(2) - pa(2) ) / ( pb(2) - pa(2) )
    pa(1:2) = q(2)
    clip_a = .true.
  else if ( p1(2) <= pa(2) .and. pb(2) < p1(2) ) then
    q(2) = p1(2)
    q(1) = pa(1) + ( pb(1) - pa(1) ) * ( q(2) - pa(2) ) / ( pb(2) - pa(2) )
    pb(1:2) = q(1:2)
    clip_b = .true.
  end if
!
!  Require that Y <= YMAX.
!
  if ( p2(2) < pa(2) .and. p2(2) < pb(2) ) then
    ival = -1
    return
  end if

  if ( p2(2) < pa(2) .and. pb(2) <= p2(2) ) then
    q(2) = p2(2)
    q(1) = pa(1) + ( pb(1) - pa(1) ) * ( q(2) - pa(2) ) / ( pb(2) - pa(2) )
    pa(1:2) = q(1:2)
    clip_a = .true.
  else if ( pa(2) <= p2(2) .and. p2(2) < pb(2) ) then
    q(2) = p2(2)
    q(1) = pa(1) + ( pb(1) - pa(1) ) * ( p2(2) - pa(2) ) / ( pb(2) - pa(2) )
    pb(1:2) = q(1:2)
    clip_b = .true.
  end if

  ival = 0

  if ( clip_a ) then
    ival = ival + 1
  end if

  if ( clip_b ) then
    ival = ival + 2
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
function ch_is_digit ( c )

!*****************************************************************************80
!
!! CH_IS_DIGIT returns TRUE if a character is a decimal digit.
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
!    Input, character C, the character to be analyzed.
!
!    Output, logical CH_IS_DIGIT, TRUE if C is a digit, FALSE otherwise.
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
subroutine data_d2_read ( file_name, n, x, y, ierror )

!*****************************************************************************80
!
!! DATA_D2_READ reads a data set of double precision pairs stored in a file.
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
!    02 December 2002
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
!    Output, integer ( kind = 4 ) IERROR, nonzero on error.
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

  open ( unit = input, file = file_name, status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DATA_D2_READ - Error!'
    write ( *, '(a)' ) '  Could not open the named file.'
    ierror = 1
    return
  end if

  x(1:n) = huge ( x(1) )
  y(1:n) = huge ( y(1) )

  n2 = 0

  do

    read ( input, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      exit
    end if

    if ( len_trim ( line ) == 0 ) then

    else if ( line(1:1) == '#' ) then

    else

      n2 = n2 + 1

      last = 0
      call s_to_r8 ( line(last+1:), x(n2), ierror, length )

      if ( ierror /= 0 ) then
        exit
      end if

      last = last + length

      call s_to_r8 ( line(last+1:), y(n2), ierror, length )

      if ( ierror /= 0 ) then
        exit
      end if

      if ( n2 == n ) then
        exit
      end if

    end if

  end do

  close ( unit = input )

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
!    02 December 2002
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
!    Output, integer ( kind = 4 ) IERROR, nonzero on error.
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
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DATA_SIZE - Error!'
    write ( *, '(a)' ) '  The file could not be opened.'
    write ( *, '(a)' ) '  Possibly the name is wrong.'
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
subroutine direction_field ( n, x, y, u, v, file_name, magnify )

!*****************************************************************************80
!
!! DIRECTION_FIELD plots a vector field.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the location of the points at
!    which vector values are known.
!
!    Input, real ( kind = 8 ) U(N), V(N), the X and Y components of the
!    vector quantity.
!
!    Input, character ( len = * ) FILE_NAME, the name of the file.
!
!    Input, real ( kind = 8 ) MAGNIFY, a magnification factor for the vectors.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) area
  real ( kind = 8 ) blue
  character ( len = * ) file_name
  real ( kind = 8 ) green
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iunit
  real ( kind = 8 ) mag
  real ( kind = 8 ) mag_max
  real ( kind = 8 ) mag_min
  real ( kind = 8 ) magnify
  real ( kind = 8 ) red
  real ( kind = 8 ) scale
  real ( kind = 8 ) u(n)
  real ( kind = 8 ) v(n)
  real ( kind = 8 ) vscale
  real ( kind = 8 ) x(n)
  integer ( kind = 4 ) :: x_ps_max = 612
  integer ( kind = 4 ) :: x_ps_min = 36
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmax2
  real ( kind = 8 ) xmin
  real ( kind = 8 ) xmin2
  real ( kind = 8 ) xscale
  real ( kind = 8 ) y(n)
  integer ( kind = 4 ) :: y_ps_max = 792
  integer ( kind = 4 ) :: y_ps_min = 36
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymax2
  real ( kind = 8 ) ymin
  real ( kind = 8 ) ymin2
!
!  Get the maximum vector length.
!
  mag_max = maxval ( sqrt ( u(1:n)**2 + v(1:n)**2 ) )

  if ( mag_max == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DIRECTION_FIELD - Error!'
    write ( *, '(a)' ) '  The direction field is null.'
    return
  end if
!
!  Open the file.
!
  call get_unit ( iunit )

  call ps_file_open ( file_name, iunit, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DIRECTION_FIELD'
    write ( *, '(a,i6)' ) '  File creation error ', ierror
    return
  end if

  call eps_file_head ( file_name, x_ps_min, y_ps_min, x_ps_max, y_ps_max )
!
!  Define the size of the page.
!
  xmin = minval ( x(1:n) )
  ymin = minval ( y(1:n) )

  xmax = maxval ( x(1:n) )
  ymax = maxval ( y(1:n) )

  area = ( xmax - xmin ) * ( ymax - ymin )

  xscale = sqrt ( area / real ( n, kind = 8 ) )
!
!  Move the boundaries out far enough to catch a typical vector.
!
  xmax = xmax + xscale
  xmin = xmin - xscale
  ymax = ymax + xscale
  ymin = ymin - xscale
!
!  Adjust symmetrically so that X and Y ranges are equal.
!
  if ( ( ymax - ymin ) < ( xmax - xmin ) ) then
    xmax2 = xmax
    xmin2 = xmin
    ymax2 = ymax + 0.5D+00 * ( ( xmax - xmin ) - ( ymax - ymin ) )
    ymin2 = ymin - 0.5D+00 * ( ( xmax - xmin ) - ( ymax - ymin ) )
  else
    xmax2 = xmax + 0.5D+00 * ( ( ymax - ymin ) - ( xmax - xmin ) )
    xmin2 = xmin - 0.5D+00 * ( ( ymax - ymin ) - ( xmax - xmin ) )
    ymax2 = ymax
    ymin2 = ymin
  end if

  call ps_page_head ( xmin2, ymin2, xmax2, ymax2 )
!
!  Mark the points.
!
  call ps_mark_disks ( n, x, y )
!
!  Compute the velocity magnitude range.
!
  mag_min =  huge ( mag_min )
  mag_max = -huge ( mag_max )

  do i = 1, n
    mag = sqrt ( u(i)**2 + v(i)**2 )
    mag_min = min ( mag_min, mag )
    mag_max = max ( mag_max, mag )
  end do

  if ( mag_max == mag_min ) then

    mag = sqrt ( u(1)**2 + v(1)**2 )

    if ( mag /= 0.0D+00 ) then
      mag_min = mag / 2.0D+00
      mag_max = mag + mag_min
    else
      mag_min = 0.0D+00
      mag_max = 1.0D+00
    end if

  end if
!
!  Draw the direction field.
!
  do i = 1, n

    mag = sqrt ( u(i)**2 + v(i)**2 )

    red = ( mag_max - mag ) / ( mag_max - mag_min )
    green = 0.0D+00
    blue = ( mag - mag_min ) / ( mag_max - mag_min )

    call ps_color_line_set ( red, green, blue )

    if ( mag <= 0.0001D+00 * mag_max ) then

    else
      scale = magnify * sqrt ( area / real ( n, kind = 8 ) ) / mag
      call ps_arrow ( x(i), y(i), x(i) + u(i) * scale, y(i) + v(i) * scale )
    end if

  end do
!
!  Close up the page and the file.
!
  call ps_page_tail ( )

  call eps_file_tail ( )

  call ps_file_close ( iunit )

  return
end
subroutine eps_file_head ( file_name, x_ps_min, y_ps_min, x_ps_max, &
  y_ps_max )

!*****************************************************************************80
!
!! EPS_FILE_HEAD writes header information to an encapsulated PostScript file.
!
!  Discussion:
!
!    The file should contain the description of only one page, but this
!    is not currently checked.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the output file.
!
!    Input, integer ( kind = 4 ) X_PS_MIN, Y_PS_MIN, X_PS_MAX, Y_PS_MAX, the minimum
!    and maximum X and Y values of the data, in PostScript units.  Any data
!    that lies outside this range will not show up properly.  A reasonable
!    set of values might be 0, 0, 612, 792, or, for a half inch margin,
!    36, 36, 576, 756.
!
  implicit none

  character ( len = 8 ) date
  character ( len = * ) file_name
  real ( kind = 8 ) line_blue
  real ( kind = 8 ) line_green
  real ( kind = 8 ) line_red
  integer ( kind = 4 ) state
  integer ( kind = 4 ) unit
  integer ( kind = 4 ) x_ps_max
  integer ( kind = 4 ) x_ps_min
  integer ( kind = 4 ) y_ps_max
  integer ( kind = 4 ) y_ps_min
!
!  Determine if the PostScript state is acceptable.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state /= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EPS_FILE_HEAD - Fatal error!'
    write ( *, '(a,i9)' ) '  PostScript state is ', state
    write ( *, '(a)' ) '  PostScript state 1 is required.'
    return
  end if
!
!  Initialization
!
  call ps_default ( )
!
!  Get the unit number.
!
  call ps_setting_int ( 'GET', 'UNIT', unit )

  call date_and_time ( date )
!
!  Write the prolog.
!
  write ( unit, '(a)' )     '%!PS-Adobe-3.0 EPSF-3.0'
  write ( unit, '(a)' )     '%%Creator: ps_write.f90'
  write ( unit, '(a)' )     '%%Title: ' // trim ( file_name )
  write ( unit, '(a)' )     '%%CreationDate: '// trim ( date )
  write ( unit, '(a)' )     '%%Pages: 1'
  write ( unit, '(a,4i6)' ) '%%BoundingBox:', &
    x_ps_min, y_ps_min, x_ps_max, y_ps_max
  write ( unit, '(a)' )     '%%Document-Fonts: Times-Roman'
  write ( unit, '(a)' )     '%%LanguageLevel: 1'
  write ( unit, '(a)' )     '%%EndComments'
  write ( unit, '(a)' )     '%%BeginProlog'
  write ( unit, '(a)' )     '/inch {72 mul} def'
  write ( unit, '(a)' )     '%%EndProlog'
!
!  Set the font.
!
  write ( unit, '(a)' ) '/Times-Roman findfont'
  write ( unit, '(a)' ) '1.00 inch scalefont'
  write ( unit, '(a)' ) 'setfont'
!
!  Set the line color.
!
  line_red = 0.0D+00
  line_green = 0.0D+00
  line_blue = 0.0D+00

  call ps_color_line ( 'SET', line_red, line_green, line_blue )
!
!  Reset the state.
!
  state = 2

  call ps_setting_int ( 'SET', 'STATE', state )

  return
end
subroutine eps_file_tail ( )

!*****************************************************************************80
!
!! EPS_FILE_TAIL writes trailer information to an encapsulated PostScript file.
!
!  Discussion:
!
!    Looks like that penultimate 'end' line is not wanted, so I commented
!    it out.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    None
!
  implicit none

  integer ( kind = 4 ) num_pages
  integer ( kind = 4 ) state
  integer ( kind = 4 ) unit
!
!  Determine if the PostScript state is acceptable.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state == 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EPS_FILE_TAIL - Warning!'
    write ( *, '(a)' ) '  A page was open.  It is being forced closed.'
    state = 2
    call ps_setting_int ( 'SET', 'STATE', state )
  end if

  if ( state /= 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EPS_FILE_TAIL - Fatal error!'
    write ( *, '(a,i9)' ) '  PostScript state is ', state
    write ( *, '(a)' ) '  PostScript state 2 is required.'
    return
  end if
!
!  Get the unit number.
!
  call ps_setting_int ( 'GET', 'UNIT', unit )
!
!  Retrieve the number of pages.
!
  call ps_setting_int ( 'GET', 'NUM_PAGES', num_pages )

  if ( num_pages > 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EPS_FILE_TAIL - Warning!'
    write ( *, '(a)' ) '  An encapsulated PostScript file describes ONE page.'
    write ( *, '(a,i9,a)' ) '  This file describes ', num_pages, ' pages.'
    write ( *, '(a)' ) '  It is not a legal EPS file.'
  end if
!
!  Write the epilog.
!
  write ( unit, '(a)' ) '%%Trailer'
! write ( unit, '(a)' ) 'end'
  write ( unit, '(a)' ) '%%EOF'
!
!  Zero out the number of pages.
!
  num_pages = 0

  call ps_setting_int ( 'SET', 'NUM_PAGES', num_pages )
!
!  Reset the state.
!
  state = 4

  call ps_setting_int ( 'SET', 'STATE', state )

  return
end
subroutine file_name_append ( file_name, append )

!*****************************************************************************80
!
!! FILE_NAME_APPEND appends a string to a filename, before the extension.
!
!  Discussion:
!
!    The "extension" of a filename is the string of characters
!    that appears after the LAST period in the name.
!
!    A file with no period, or with a period as the last character
!    in the name, has a "null" extension.
!
!    The idea is that the string in APPEND is to be appended to
!    the part of the filename that precedes the extension.
!
!    The intended purpose of this routine is to be able to easily
!    generate a filename that indicates its relation to another file.
!
!  Example:
!
!          Input             Output
!    ===================     =========
!    FILE_NAME    APPEND     FILE_NAME
!
!    bob.for      6          bob6.for
!    bob.bob.bob  JOB        bob.bobJOB.bob
!    bob          yak        bobyak
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) FILE_NAME, a file name.
!    On output, the file name has been modified.
!
!    Input, character ( len = * ) APPEND, the string to be appended.
!
  implicit none

  character ( len = * ) append
  character ( len = * ) file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) len_append
  integer ( kind = 4 ) len_name

  call file_name_ext_get ( file_name, i, j )
!
!  If there is no extension, then simply slap APPEND on the end.
!
  if ( i == 0 ) then

    len_name = len_trim ( file_name )
    file_name(len_name+1:) = append
!
!  If there is an extension, then insert APPEND.
!
  else

    len_append = len_trim ( append )
    file_name(i:) = append(1:len_append) // file_name(i:j)

  end if

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
!  Examples:
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
subroutine file_name_ext_swap ( file_name, ext )

!*****************************************************************************80
!
!! FILE_NAME_EXT_SWAP replaces the current "extension" of a file name.
!
!  Discussion:
!
!    The "extension" of a filename is the string of characters
!    that appears after the LAST period in the name.  A file
!    with no period, or with a period as the last character
!    in the name, has a "null" extension.
!
!  Example:
!
!          Input           Output
!    ================     =========
!    FILE_NAME    EXT     FILE_NAME
!
!    bob.for      obj     bob.obj
!    bob.bob.bob  txt     bob.bob.txt
!    bob          yak     bob.yak
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
!    Input/output, character ( len = * ) FILE_NAME, a file name.
!    On output, the extension of the file has been changed.
!
!    Input, character ( len = * ) EXT, the extension to be used on the output
!    copy of FILE_NAME, replacing the current extension if any.
!
  implicit none

  character ( len = * ) ext
  character ( len = * ) file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) len_max
  integer ( kind = 4 ) len_name

  len_max = len ( file_name )
  len_name = len_trim ( file_name )

  call file_name_ext_get ( file_name, i, j )

  if ( i == 0 ) then

    if ( len_max < len_name + 1 ) then
      return
    end if

    len_name = len_name + 1
    file_name(len_name:len_name) = '.'
    i = len_name + 1

  else

    i = i + 1
    file_name(i:j) = ' '

  end if

  file_name(i:) = ext

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
subroutine i4_swap ( i, j )

!*****************************************************************************80
!
!! I4_SWAP swaps two integer values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) I, J.  On output, the values of I and
!    J have been interchanged.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k

  k = i
  i = j
  j = k

  return
end
function point_inside_box_2d ( x1, y1, x2, y2, x, y )

!*****************************************************************************80
!
!! POINT_INSIDE_BOX_2D determines if a point is inside a box in 2D.
!
!  Discussion:
!
!    A "box" is defined by its "left down" corner and its
!    "right up" corner, and all the points between.  It is
!    assumed that the sides of the box align with coordinate directions.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X1, Y1, X2, Y2, the two corners of the box.
!
!    Input, real ( kind = 8 ) X, Y, the point to be checked.
!
!    Output, logical POINT_INSIDE_BOX_2D, is .TRUE. if (X,Y) is inside the
!    box, or on its boundary, and .FALSE. otherwise.
!
  implicit none

  logical point_inside_box_2d
  real ( kind = 8 ) x
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) y
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2

  if ( x1 <= x .and. x <= x2 .and. &
       y1 <= y .and. y <= y2 ) then
    point_inside_box_2d = .true.
  else
    point_inside_box_2d = .false.
  end if

  return
end
subroutine ps_arrow ( x1, y1, x2, y2 )

!*****************************************************************************80
!
!! PS_ARROW draws an arrow from (X1,Y1) to (X2,Y2).
!
!  Discussion:
!
!    The current point is set to (X2,Y2).
!
!    This routine will clip the line, if necessary, so that the line
!    drawn is entirely within the region.
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
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X1, Y1, the starting point of the arrow.
!
!    Input, real ( kind = 8 ) X2, Y2, the ending point of the arrow.
!
  implicit none

  real ( kind = 8 ) alpha
  real ( kind = 8 ) alpha2
  real ( kind = 8 ) alpha3
  real ( kind = 8 ) frac
  integer ( kind = 4 ) ival
  real ( kind = 8 ) p1(2)
  real ( kind = 8 ) p2(2)
  real ( kind = 8 ) pa(2)
  real ( kind = 8 ) pb(2)
  integer ( kind = 4 ) plotxmin2
  integer ( kind = 4 ) plotymin2
  integer ( kind = 4 ) px
  integer ( kind = 4 ) py
  real ( kind = 8 ) r
  real ( kind = 8 ) r2
  integer ( kind = 4 ) state
  integer ( kind = 4 ) unit
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ) x4
  real ( kind = 8 ) x5
  real ( kind = 8 ) x6
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) y3
  real ( kind = 8 ) y4
  real ( kind = 8 ) y5
  real ( kind = 8 ) y6
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin

  if ( x1 == x2 .and. y1 == y2 ) then
    return
  end if
!
!  Determine if the PostScript state is acceptable.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state /= 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_ARROW - Fatal error!'
    write ( *, '(a,i9)' ) '  PostScript state is ', state
    write ( *, '(a)' ) '  PostScript state 3 is required.'
    return
  end if
!
!  Get settings.
!
  call ps_setting_int ( 'GET', 'PXMIN', plotxmin2 )

  call ps_setting_int ( 'GET', 'PYMIN', plotymin2 )

  call ps_setting_int ( 'GET', 'UNIT', unit )

  call ps_setting_real ( 'GET', 'ALPHA', alpha )
  call ps_setting_real ( 'GET', 'XMIN', xmin )
  call ps_setting_real ( 'GET', 'XMAX', xmax )
  call ps_setting_real ( 'GET', 'YMIN', ymin )
  call ps_setting_real ( 'GET', 'YMAX', ymax )
!
!  Clip the line.
!
  p1(1:2) = (/ xmin, ymin /)
  p2(1:2) = (/ xmax, ymax /)
  pa(1:2) = (/ x1, y1 /)
  pb(1:2) = (/ x2, y2 /)

  call box_segment_clip_2d ( p1, p2, pa, pb, ival )

  if ( ival < 0 ) then
    write ( *, * ) ' '
    write ( *, * ) xmin, ymin
    write ( *, * ) xmax, ymax
    write ( *, * ) x1, y1
    write ( *, * ) x2, y2
    return
  end if

  x3 = pa(1)
  y3 = pa(2)
  x4 = pb(1)
  y4 = pb(2)

  r = sqrt ( ( x4 - x3 )**2 + ( y4 - y3 )**2 )

  if ( r == 0.0D+00 ) then
    return
  end if
!
!  FRAC controls the size of the arrow head.  It's specified
!  as a proportion of the length of the line.
!
  frac = 0.1D+00

  r2 = sqrt ( frac**2 + ( 1.0D+00 - frac )**2 ) * r

  alpha2 = atan2 ( y4 - y3, x4 - x3 )
  alpha3 = atan2 ( frac, 1.0D+00 - frac )

  x5 = x3 + r2 * cos ( alpha2 - alpha3 )
  y5 = y3 + r2 * sin ( alpha2 - alpha3 )

  x6 = x3 + r2 * cos ( alpha2 + alpha3 )
  y6 = y3 + r2 * sin ( alpha2 + alpha3 )
!
!  Draw line.
!
  write ( unit, '(a)' ) 'newpath'

  px = plotxmin2 + nint ( alpha * ( x3 - xmin ) )
  py = plotymin2 + nint ( alpha * ( y3 - ymin ) )
  write ( unit, '(2i6,a)' ) px, py, ' moveto'

  px = plotxmin2 + nint ( alpha * ( x4 - xmin ) )
  py = plotymin2 + nint ( alpha * ( y4 - ymin ) )
  write ( unit, '(2i6,a)' ) px, py, ' lineto'
!
!  Draw arrow head.
!
  px = plotxmin2 + nint ( alpha * ( x4 - xmin ) )
  py = plotymin2 + nint ( alpha * ( y4 - ymin ) )
  write ( unit, '(2i6,a)' ) px, py, ' moveto'

  px = plotxmin2 + nint ( alpha * ( x5 - xmin ) )
  py = plotymin2 + nint ( alpha * ( y5 - ymin ) )
  write ( unit, '(2i6,a)' ) px, py, ' lineto'

  px = plotxmin2 + nint ( alpha * ( x4 - xmin ) )
  py = plotymin2 + nint ( alpha * ( y4 - ymin ) )
  write ( unit, '(2i6,a)' ) px, py, ' moveto'

  px = plotxmin2 + nint ( alpha * ( x6 - xmin ) )
  py = plotymin2 + nint ( alpha * ( y6 - ymin ) )
  write ( unit, '(2i6,a)' ) px, py, ' lineto stroke'

  call ps_setting_real ( 'SET', 'XCUR', x2 )
  call ps_setting_real ( 'SET', 'YCUR', y2 )

  return
end
subroutine ps_color_line ( action, r, g, b )

!*****************************************************************************80
!
!! PS_COLOR_LINE handles the line color.
!
!  Discussion:
!
!    By calling this routine, you can temporarily set the line color,
!    draw some lines, and then restore it to whatever it was.
!
!    An earlier version of this routine did not use the SAVE command for
!    the stack arrrays, meaning the stored data was lost.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION, the desired action.
!    'SET', set the line color to RGB.
!    'GET', set RGB to the current line color.
!    'PUSH', push a value onto the RGB stack.
!    'POP', pop the RGB stack.
!
!    Input, real ( kind = 8 ) R, G, B, the RGB values for the new line color.
!
  implicit none

  integer ( kind = 4 ), parameter :: nstack = 10

  character ( len = * ) action
  real ( kind = 8 ) b
  real ( kind = 8 ) b_old
  real ( kind = 8 ), save, dimension ( nstack) :: b_stack
  real ( kind = 8 ) g
  real ( kind = 8 ) g_old
  real ( kind = 8 ), save, dimension ( nstack) :: g_stack
  integer ( kind = 4 ), save :: istack = 0
  real ( kind = 8 ) r
  real ( kind = 8 ) r_old
  real ( kind = 8 ), save, dimension ( nstack) :: r_stack
  logical s_eqi

  if ( s_eqi ( action, 'SET' ) ) then

    call ps_color_line_set ( r, g, b )

  else if ( s_eqi ( action, 'GET' ) ) then

    call ps_setting_real ( 'GET', 'LINE_RED', r )
    call ps_setting_real ( 'GET', 'LINE_GREEN', g )
    call ps_setting_real ( 'GET', 'LINE_BLUE', b )

  else if ( s_eqi ( action, 'POP' ) ) then

    if ( 0 < istack ) then
      r = r_stack(istack)
      g = g_stack(istack)
      b = b_stack(istack)
      istack = istack - 1
    end if

    call ps_color_line_set ( r, g, b )

  else if ( s_eqi ( action, 'PUSH' ) ) then

    call ps_setting_real ( 'GET', 'LINE_RED', r_old )
    call ps_setting_real ( 'GET', 'LINE_GREEN', g_old )
    call ps_setting_real ( 'GET', 'LINE_BLUE', b_old )

    if ( istack <= nstack ) then
      istack = istack + 1
      r_stack(istack) = r_old
      g_stack(istack) = g_old
      b_stack(istack) = b_old
    end if

    call ps_color_line_set ( r, g, b )

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_COLOR_LINE - Fatal error!'
    write ( *, '(a)' ) '  Unexpected ACTION.'
    stop

  end if

  return
end
subroutine ps_color_line_set ( r, g, b )

!*****************************************************************************80
!
!! PS_COLOR_LINE_SET sets the line color.
!
!  Discussion:
!
!    By calling this routine, you guarantee that a check will be made
!    of the current line color.  If the current and new line colors are
!    the same, then we skip the extraneous action of setting the color.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, G, B, the RGB values for the new line color.
!
  implicit none

  real ( kind = 8 ) b
  real ( kind = 8 ) b_old
  real ( kind = 8 ) g
  real ( kind = 8 ) g_old
  real ( kind = 8 ) r
  real ( kind = 8 ) r_old
  integer ( kind = 4 ) state
  integer ( kind = 4 ) unit
!
!  Check the state.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_COLOR_LINE_SET - Fatal error!'
    write ( *, '(a,i9)' ) '  PostScript state is ', state
    write ( *, '(a)' ) '  A PostScript state of at least 1 is required.'
    return
  end if
!
!  Get the current colors.
!
  call ps_setting_real ( 'GET', 'LINE_RED', r_old )
  call ps_setting_real ( 'GET', 'LINE_GREEN', g_old )
  call ps_setting_real ( 'GET', 'LINE_BLUE', b_old )
!
!  If any color has changed, we need to reset them.
!
  if ( r_old /= r .or. g_old /= g .or. b_old /= b ) then

    call ps_setting_int ( 'GET', 'UNIT', unit )

    call ps_comment ( 'Set RGB line color.' )

    write ( unit, '(3f7.4,a)' ) r, g, b, ' setrgbcolor'

    call ps_setting_real ( 'SET', 'LINE_RED', r )
    call ps_setting_real ( 'SET', 'LINE_GREEN', g )
    call ps_setting_real ( 'SET', 'LINE_BLUE', b )

  end if

  return
end
subroutine ps_comment ( string )

!*****************************************************************************80
!
!! PS_COMMENT inserts a comment into the PostScript file.
!
!  Discussion:
!
!    A comment begins with a percent sign in column 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    Input, character ( len = * ) STRING, the comment.
!
  implicit none

  character ( len = * ) string
  integer ( kind = 4 ) unit
!
!  Get the unit number.
!
  call ps_setting_int ( 'GET', 'UNIT', unit )
!
!  Write the comment.
!
  if ( len_trim ( string ) == 0 ) then
    write ( unit, '(a)' ) '%'
  else
    write ( unit, '(a)' ) '%'
    write ( unit, '(a2,a)' ) '% ', trim ( string )
    write ( unit, '(a)' ) '%'
  end if

  return
end
subroutine ps_default ( )

!*****************************************************************************80
!
!! PS_DEFAULT sets the internal settings to their default values
!
!  Discussion:
!
!    Certain variables are not reset, including the number of pages,
!    the unit number, the internal state, and variables relating to
!    the size and shape of the region.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    None
!
  implicit none

  real ( kind = 8 ) fill_blue
  real ( kind = 8 ) fill_green
  real ( kind = 8 ) fill_red
  real ( kind = 8 ) font_size
  real ( kind = 8 ) line_blue
  real ( kind = 8 ) line_green
  real ( kind = 8 ) line_red
  integer ( kind = 4 ) line_width
  integer ( kind = 4 ) marker_size

  line_width = 1
  marker_size = 5

  call ps_setting_int ( 'SET', 'LINE_WIDTH', line_width )
  call ps_setting_int ( 'SET', 'MARKER_SIZE', marker_size )

  fill_blue = 0.7D+00
  fill_green = 0.7D+00
  fill_red = 0.7D+00
  font_size = 0.1D+00
  line_blue = 0.0D+00
  line_green = 0.0D+00
  line_red = 0.0D+00

  call ps_setting_real ( 'SET', 'FILL_BLUE', fill_blue )
  call ps_setting_real ( 'SET', 'FILL_GREEN', fill_green )
  call ps_setting_real ( 'SET', 'FILL_RED', fill_red )
  call ps_setting_real ( 'SET', 'FONT_SIZE', font_size )
  call ps_setting_real ( 'SET', 'LINE_BLUE', line_blue )
  call ps_setting_real ( 'SET', 'LINE_GREEN', line_green )
  call ps_setting_real ( 'SET', 'LINE_RED', line_red )

  return
end
subroutine ps_file_close ( unit )

!*****************************************************************************80
!
!! PS_FILE_CLOSE closes a PostScript file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) UNIT, the FORTRAN unit to which output was written.
!
  implicit none

  integer ( kind = 4 ) state
  integer ( kind = 4 ) unit
!
!  Determine if the PostScript state is acceptable.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state < 1 .or. state > 4 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_FILE_CLOSE - Fatal error!'
    write ( *, '(a,i9)' ) '  PostScript state is ', state
    write ( *, '(a)' ) '  PostScript state 1, 2, 3 or 4 is required.'
    return
  end if

  close ( unit = unit )

  state = 0
  call ps_setting_int ( 'SET', 'STATE', state )

  unit = 0
  call ps_setting_int ( 'SET', 'UNIT', unit )

  return
end
subroutine ps_file_head ( file_name )

!*****************************************************************************80
!
!! PS_FILE_HEAD writes header information to a PostScript file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the output file.
!
  implicit none

  character ( len = 8 ) date
  character ( len = * ) file_name
  real ( kind = 8 ) line_blue
  real ( kind = 8 ) line_green
  real ( kind = 8 ) line_red
  integer ( kind = 4 ) margin
  integer ( kind = 4 ) pagexmax
  integer ( kind = 4 ) pagexmin
  integer ( kind = 4 ) pageymax
  integer ( kind = 4 ) pageymin
  integer ( kind = 4 ) plotxmax
  integer ( kind = 4 ) plotxmin
  integer ( kind = 4 ) plotymax
  integer ( kind = 4 ) plotymin
  integer ( kind = 4 ) state
  integer ( kind = 4 ) unit
!
!  Determine if the PostScript state is acceptable.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state /= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_FILE_HEAD - Fatal error!'
    write ( *, '(a,i9)' ) '  PostScript state is ', state
    write ( *, '(a)' ) '  PostScript state 1 is required.'
    return
  end if
!
!  Initialization
!
  call ps_default ( )
!
!  Compute the scale factor.
!
  pagexmax = 612
  pagexmin = 0
  pageymax = 792
  pageymin = 0

  margin = 36

  plotxmax = pagexmax - margin
  plotxmin = pagexmin + margin
  plotymax = pageymax - margin
  plotymin = pageymin + margin
!
!  Get the unit number.
!
  call ps_setting_int ( 'GET', 'UNIT', unit )

  call date_and_time ( date )
!
!  Write the prolog.
!
  write ( unit, '(a)' )     '%!PS-Adobe-1.0'
  write ( unit, '(a)' )     '%%Creator: ps_write.f90'
  write ( unit, '(a)' )     '%%Title: ' // trim ( file_name )
  write ( unit, '(a)' )     '%%CreationDate: ' // trim ( date )
  write ( unit, '(a)' )     '%%Pages: (atend)'
  write ( unit, '(a,4i6)' ) '%%BoundingBox:', plotxmin, plotymin, plotxmax, &
    plotymax
  write ( unit, '(a)' )     '%%Document-Fonts: Times-Roman'
  write ( unit, '(a)' )     '%%LanguageLevel: 1'
  write ( unit, '(a)' )     '%%EndComments'
  write ( unit, '(a)' )     '%%BeginProlog'
  write ( unit, '(a)' )     '/inch {72 mul} def'
  write ( unit, '(a)' )     '%%EndProlog'
!
!  Set the font.
!
  call ps_comment ( 'Set the font:' )

  write ( unit, '(a)' ) '/Times-Roman findfont'
  write ( unit, '(a)' ) '1.00 inch scalefont'
  write ( unit, '(a)' ) 'setfont'
!
!  Set the line color.
!
  line_red = 0.0D+00
  line_green = 0.0D+00
  line_blue = 0.0D+00

  call ps_color_line ( 'SET', line_red, line_green, line_blue )
!
!  Reset the state.
!
  state = 2

  call ps_setting_int ( 'SET', 'STATE', state )

  return
end
subroutine ps_file_open ( file_name, unit, ierror )

!*****************************************************************************80
!
!! PS_FILE_OPEN opens a new version of a PostScript file with a given name.
!
!  Discussion:
!
!    If a file of the given name already exists, it is deleted.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) UNIT, the FORTRAN unit to which output should
!    be written.
!
!    Input, character ( len = 80 ) FILE_NAME, the name of the output file.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error.
!    nonzero, the file could not be created.
!
  implicit none

  character ( len = * ) file_name
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) state
  integer ( kind = 4 ) unit
!
!  Determine if the PostScript state is acceptable.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_FILE_OPEN - Fatal error!'
    write ( *, '(a,i9)' ) '  PostScript state is ', state
    write ( *, '(a)' ) '  PostScript state 0 is required.'
    write ( *, '(a)' ) '  Call PS_FILE_CLOSE first!'
    return
  end if

  ierror = 0
!
!  Now create a new empty file of the given name.
!
  open ( unit = unit, file = file_name, status = 'replace', iostat = ios )

  if ( ios /= 0 ) then
    ierror = ios
    return
  end if

  state = 1
  call ps_setting_int ( 'SET', 'STATE', state )

  call ps_setting_int ( 'SET', 'UNIT', unit )

  return
end
subroutine ps_file_tail ( )

!*****************************************************************************80
!
!! PS_FILE_TAIL writes trailer information to a PostScript file.
!
!  Discussion:
!
!    Looks like that penultimate 'end' line is not wanted, so
!    I commented it out.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    None
!
  implicit none

  integer ( kind = 4 ) num_pages
  integer ( kind = 4 ) state
  integer ( kind = 4 ) unit
!
!  Determine if the PostScript state is acceptable.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state == 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_FILE_TAIL - Warning!'
    write ( *, '(a)' ) '  A page was open.  It is being forced closed.'
    state = 2
    call ps_setting_int ( 'SET', 'STATE', state )
  end if

  if ( state /= 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_FILE_TAIL - Fatal error!'
    write ( *, '(a,i9)' ) '  PostScript state is ', state
    write ( *, '(a)' ) '  PostScript state 2 is required.'
    return
  end if
!
!  Get the unit number.
!
  call ps_setting_int ( 'GET', 'UNIT', unit )
!
!  Retrieve the number of pages.
!
  call ps_setting_int ( 'GET', 'NUM_PAGES', num_pages )
!
!  Write the epilog.
!
  write ( unit, '(a)' ) '%%Trailer'
  write ( unit, '(a,i6)' ) '%%Pages: ', num_pages
! write ( unit, '(a)' ) 'end'
  write ( unit, '(a)' ) '%%EOF'
!
!  Zero out the number of pages.
!
  num_pages = 0

  call ps_setting_int ( 'SET', 'NUM_PAGES', num_pages )
!
!  Reset the state.
!
  state = 4

  call ps_setting_int ( 'SET', 'STATE', state )

  return
end
subroutine ps_line_closed ( npoint, x, y )

!*****************************************************************************80
!
!! PS_LINE_CLOSED adds the graph of a closed line to a PostScript file.
!
!  Discussion:
!
!    A "closed" line is one in which the last point is connected back
!    to the first one.
!
!    The current point is set to the first (and logically last) point
!    in the list.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NPOINT, the number of points in the line.
!
!    Input, real ( kind = 8 ) X(NPOINT), Y(NPOINT), the X and Y components
!    of the points.
!
  implicit none

  integer ( kind = 4 ) npoint

  real ( kind = 8 ) alpha
  integer ( kind = 4 ) i
  integer ( kind = 4 ) plotxmin2
  integer ( kind = 4 ) plotymin2
  integer ( kind = 4 ) px
  integer ( kind = 4 ) py
  integer ( kind = 4 ) state
  integer ( kind = 4 ) unit
  real ( kind = 8 ) x(npoint)
  real ( kind = 8 ) xmin
  real ( kind = 8 ) y(npoint)
  real ( kind = 8 ) ymin
!
!  Refuse to handle fewer than 2 points.
!
  if ( npoint < 2 ) then
    return
  end if
!
!  Determine if the PostScript state is acceptable.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state /= 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_LINE_CLOSED - Fatal error!'
    write ( *, '(a,i9)' ) '  PostScript state is ', state
    write ( *, '(a)' ) '  PostScript state 3 is required.'
    return
  end if
!
!  Get settings.
!
  call ps_setting_int ( 'GET', 'PXMIN', plotxmin2 )

  call ps_setting_int ( 'GET', 'PYMIN', plotymin2 )

  call ps_setting_int ( 'GET', 'UNIT', unit )

  call ps_setting_real ( 'GET', 'ALPHA', alpha )
  call ps_setting_real ( 'GET', 'XMIN', xmin )
  call ps_setting_real ( 'GET', 'YMIN', ymin )
!
!  Draw lines.
!
  write ( unit, '(a)' ) 'newpath'

  px = plotxmin2 + nint ( alpha * ( x(1) - xmin ) )
  py = plotymin2 + nint ( alpha * ( y(1) - ymin ) )
  write ( unit, '(2i6,a)' ) px, py, ' moveto'

  do i = 2, npoint
    px = plotxmin2 + nint ( alpha * ( x(i) - xmin ) )
    py = plotymin2 + nint ( alpha * ( y(i) - ymin ) )
    write ( unit, '(2i6,a)' ) px, py, ' lineto'
  end do
!
!  Add the final extra segment to the initial point.
!
  px = plotxmin2 + nint ( alpha * ( x(1) - xmin ) )
  py = plotymin2 + nint ( alpha * ( y(1) - ymin ) )
  write ( unit, '(2i6,a)' ) px, py, ' lineto'
!
!  Draw the line.
!
  write ( unit, '(a)' ) 'stroke'

  call ps_setting_real ( 'SET', 'XCUR', x(1) )
  call ps_setting_real ( 'SET', 'YCUR', y(1) )

  return
end
subroutine ps_mark_disks ( n, x, y )

!*****************************************************************************80
!
!! PS_MARK_DISKS marks points with a small filled disk.
!
!  Discussion:
!
!    The current point is set to the center of the last disk.
!
!    The circles are drawn with the current RGB fill colors.
!
!    The circles are drawn the current marker size.
!
!    Points outside the region are not marked.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton, Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the coordinates of the point to mark.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) alpha
  integer ( kind = 4 ) i
  integer ( kind = 4 ) marker_size
  integer ( kind = 4 ) plotxmin2
  integer ( kind = 4 ) plotymin2
  logical point_inside_box_2d
  integer ( kind = 4 ) pxcen
  integer ( kind = 4 ) pycen
  integer ( kind = 4 ) state
  integer ( kind = 4 ) unit
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin
!
!  Determine if the PostScript state is acceptable.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state /= 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_MARK_DISKS - Fatal error!'
    write ( *, '(a,i9)' ) '  PostScript state is ', state
    write ( *, '(a)' ) '  PostScript state 3 is required.'
    return
  end if
!
!  Get settings.
!
!  MARKER_SIZE = 5 seems way too big.
!
  call ps_setting_int ( 'GET', 'MARKER_SIZE', marker_size )
  marker_size = 1

  call ps_setting_int ( 'GET', 'PXMIN', plotxmin2 )
  call ps_setting_int ( 'GET', 'PYMIN', plotymin2 )
  call ps_setting_int ( 'GET', 'UNIT', unit )

  call ps_setting_real ( 'GET', 'ALPHA', alpha )
  call ps_setting_real ( 'GET', 'XMIN', xmin )
  call ps_setting_real ( 'GET', 'YMIN', ymin )
  call ps_setting_real ( 'GET', 'XMAX', xmax )
  call ps_setting_real ( 'GET', 'YMAX', ymax )

  write ( unit, '(a)' ) 'newpath'

  do i = 1, n

    if ( .not. point_inside_box_2d ( xmin, ymin, xmax, ymax, x(i), y(i) ) ) then
      cycle
    end if

    pxcen = plotxmin2 + nint ( alpha * ( x(i) - xmin ) )
    pycen = plotymin2 + nint ( alpha * ( y(i) - ymin ) )
    write ( unit, '(3i6,a)' ) pxcen, pycen, marker_size, &
      ' 0 360 arc closepath fill'

  end do

  call ps_setting_real ( 'SET', 'XCUR', x )
  call ps_setting_real ( 'SET', 'YCUR', y )

  return
end
subroutine ps_page_head ( xmin, ymin, xmax, ymax )

!*****************************************************************************80
!
!! PS_PAGE_HEAD writes header information on a new page.
!
!  Discussion:
!
!    I think an earlier version of this code, which wrote
!    "%% Page:" rather than "%%Page:" may have caused problems
!    for some interpreters.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 February 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XMIN, YMIN, XMAX, YMAX, the minimum and maximum X
!    and Y values of the data to be drawn on this page.
!
  implicit none

  real ( kind = 8 ) alpha
  integer ( kind = 4 ) num_pages
  integer ( kind = 4 ) state
  real ( kind = 8 ) line_blue
  real ( kind = 8 ) line_green
  real ( kind = 8 ) line_red
  integer ( kind = 4 ) margin
  integer ( kind = 4 ) pagexmax
  integer ( kind = 4 ) pagexmin
  integer ( kind = 4 ) pageymax
  integer ( kind = 4 ) pageymin
  integer ( kind = 4 ) plotxmax
  integer ( kind = 4 ) plotxmin
  integer ( kind = 4 ) plotxmin2
  integer ( kind = 4 ) plotymax
  integer ( kind = 4 ) plotymin
  integer ( kind = 4 ) plotymin2
  integer ( kind = 4 ) unit
  real ( kind = 8 ) xcur
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmax2
  real ( kind = 8 ) xmin
  real ( kind = 8 ) xmin2
  real ( kind = 8 ) xvec(4)
  real ( kind = 8 ) ycur
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymax2
  real ( kind = 8 ) ymin
  real ( kind = 8 ) ymin2
  real ( kind = 8 ) yvec(4)
!
!  Determine if the PostScript state is acceptable.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state == 3 ) then
    state = 2
    call ps_setting_int ( 'SET', 'STATE', state )
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_PAGE_HEAD - Warning!'
    write ( *, '(a)' ) '  The current open page is forced closed.'
  end if

  if ( state /= 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_PAGE_HEAD - Fatal error!'
    write ( *, '(a,i9)' ) '  PostScript state is ', state
    write ( *, '(a)' ) '  PostScript state 2 is required.'
    return
  end if
!
!  Get settings.
!
  call ps_setting_int ( 'GET', 'NUM_PAGES', num_pages )

  num_pages = num_pages + 1

  call ps_setting_int ( 'SET', 'NUM_PAGES', num_pages )

  call ps_setting_int ( 'GET', 'UNIT', unit )

  write ( unit, '(a,i6,i6)' ) '%%Page: ', num_pages, num_pages
  write ( unit, '(a)' ) 'save'
!
!  Reset the state.
!
  state = 3

  call ps_setting_int ( 'SET', 'STATE', state )
!
!  Determine and store parameters.
!
  if ( xmax == xmin ) then
    xmax2 = xmax + 1.0D+00
    xmin2 = xmax - 1.0D+00
  else
    xmax2 = xmax
    xmin2 = xmin
  end if

  if ( ymax == ymin ) then
    ymax2 = ymax + 1.0D+00
    ymin2 = ymax - 1.0D+00
  else
    ymax2 = ymax
    ymin2 = ymin
  end if
!
!  Set the value of "current point".
!
  xcur = xmin
  ycur = ymin
!
!  Set the conversion factors.
!
  pagexmax = 612
  pagexmin = 0
  pageymax = 792
  pageymin = 0

  margin = 36

  plotxmax = pagexmax - margin
  plotxmin = pagexmin + margin
  plotymax = pageymax - margin
  plotymin = pageymin + margin

  alpha = min ( real ( plotxmax - plotxmin, kind = 8 ) / ( xmax2 - xmin2 ), &
                real ( plotymax - plotymin, kind = 8 ) / ( ymax2 - ymin2 ) )
!
!  Adjust PLOTXMIN and PLOTYMIN to center the image.
!
  plotxmin2 = nint ( 0.5D+00 * &
    ( real ( plotxmin + plotxmax, kind = 8 ) - alpha * ( xmax2 - xmin2 ) ) )

  plotymin2 = nint ( 0.5D+00 * &
    ( real ( plotymin + plotymax, kind = 8 ) - alpha * ( ymax2 - ymin2 ) ) )
!
!  Store data.
!
  call ps_setting_int ( 'SET', 'PXMIN', plotxmin2 )
  call ps_setting_int ( 'SET', 'PYMIN', plotymin2 )

  call ps_setting_real ( 'SET', 'ALPHA', alpha )
  call ps_setting_real ( 'SET', 'XCUR', xcur )
  call ps_setting_real ( 'SET', 'XMIN', xmin )
  call ps_setting_real ( 'SET', 'XMAX', xmax )
  call ps_setting_real ( 'SET', 'YCUR', ycur )
  call ps_setting_real ( 'SET', 'YMIN', ymin )
  call ps_setting_real ( 'SET', 'YMAX', ymax )
!
!  Draw a gray border around the page.
!
  line_red = 0.9D+00
  line_green = 0.9D+00
  line_blue = 0.9D+00

  call ps_color_line ( 'PUSH', line_red, line_green, line_blue )

  call ps_comment ( 'Draw a gray border around the page.' )

  xvec(1:4) = (/ xmin, xmax, xmax, xmin /)
  yvec(1:4) = (/ ymin, ymin, ymax, ymax /)

  call ps_line_closed ( 4, xvec, yvec )

  call ps_color_line ( 'POP', line_red, line_green, line_blue )

  return
end
subroutine ps_page_tail ( )

!*****************************************************************************80
!
!! PS_PAGE_TAIL writes tail information at the end of a page.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    None
!
  implicit none

  integer ( kind = 4 ) state
  integer ( kind = 4 ) unit
!
!  Determine if the PostScript state is acceptable.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state /= 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_PAGE_TAIL - Fatal error!'
    write ( *, '(a,i9)' ) '  PostScript state is ', state
    write ( *, '(a)' ) '  PostScript state 3 is required.'
    return
  end if
!
!  Get settings.
!
  call ps_setting_int ( 'GET', 'UNIT', unit )

  write ( unit, '(a)' ) 'restore showpage'

  call ps_comment ( 'End of page' )
!
!  Reset the state.
!
  state = 2

  call ps_setting_int ( 'SET', 'STATE', state )

  return
end
subroutine ps_setting_int ( action, variable, value )

!*****************************************************************************80
!
!! PS_SETTING_INT sets, gets, or prints integer internal PS_WRITE parameters.
!
!  Discussion:
!
!    Normally, the user does not call this routine.  It is a utility
!    used by the package.
!
!    I'd like a more sophisticated pop and push.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION, the desired action:
!    'GET' to get the current value of VARIABLE, or
!    'POP' to return the current value and set a new value;
!    'SET' to set a new value of VARIABLE, or
!    'PUSH' to return the current value and set a new value;
!    'PRINT' to print the current value of VARIABLE.
!
!    Input, character ( len = * ) VARIABLE, the variable to get or set:
!    'LINE_WIDTH', the line width.
!      0 is the very thinnest line possible,
!      1 is more usual, 2 is thicker, and so on.
!    'MARKER_SIZE', the size of marker circles and disks, in PostScript points;
!    'NUM_PAGES', the number of pages begun or completed;
!    'PXMIN', the location of the left hand margin of the region
!       in PostScript points;
!    'PYMIN', the location of the lower margin of the region
!       in PostScript points;
!    'STATE', the current internal state,
!      0, file not open,
!      1, file open, no header written, no page open,
!      2, file open, header written, no page open,
!      3, file open, header written, page open.
!      4, file open, header written, trailer written.
!    'UNIT', the FORTRAN output unit associated with the PostScript file.
!
!    Input/output, integer ( kind = 4 ) VALUE.
!    If ACTION = 'GET', then VALUE is an output quantity, and is the
!    current internal value of the variable.
!
!    If ACTION = 'SET', then VALUE is an input quantity, and the
!    current internal value of the variable is set to this value.
!
!    If ACTION = 'PRINT', then VALUE is ignored.
!
  implicit none

  character ( len = * ) action
  integer ( kind = 4 ), save :: line_width = 1
  integer ( kind = 4 ), save :: marker_size = 0
  integer ( kind = 4 ), save :: num_pages = 0
  integer ( kind = 4 ), save :: pxmin = 0
  integer ( kind = 4 ), save :: pymin = 0
  integer ( kind = 4 ), save :: state = 0
  integer ( kind = 4 ), save :: unit = 0
  integer ( kind = 4 ) value
  character ( len = * ) variable

  if ( variable == 'LINE_WIDTH' ) then

    if ( action == 'GET' ) then
      value = line_width
    else if ( action == 'PRINT' ) then
      write ( *, '(a,i9)' ) 'Line width, LINE_WIDTH = ', line_width
    else if ( action == 'SET' ) then
      line_width = value
    else if ( action == 'POP' ) then
      call i4_swap ( line_width, value )
    else if ( action == 'PUSH' ) then
      call i4_swap ( line_width, value )
    end if

  else if ( variable == 'MARKER_SIZE' ) then

    if ( action == 'GET' ) then
      value = marker_size
    else if ( action == 'PRINT' ) then
      write ( *, '(a,i9)' ) 'Marker size, MARKER_SIZE = ', marker_size
    else if ( action == 'SET' ) then
      marker_size = value
    else if ( action == 'POP' ) then
      call i4_swap ( marker_size, value )
    else if ( action == 'PUSH' ) then
      call i4_swap ( marker_size, value )
    end if

  else if ( variable == 'NUM_PAGES' ) then

    if ( action == 'GET' ) then
      value = num_pages
    else if ( action == 'PRINT' ) then
      write ( *, '(a,i9)' ) 'Number of pages, NUM_PAGES = ', num_pages
    else if ( action == 'SET' ) then
      num_pages = value
    end if

  else if ( variable == 'PXMIN' ) then

    if ( action == 'GET' ) then
      value = pxmin
    else if ( action == 'PRINT' ) then
      write ( *, '(a,i9)' ) 'PostScript minimum X point, PXMIN = ', pxmin
    else if ( action == 'SET' ) then
      pxmin = value
    else if ( action == 'POP' ) then
      call i4_swap ( pxmin, value )
    else if ( action == 'PUSH' ) then
      call i4_swap ( pxmin, value )
    end if

  else if ( variable == 'PYMIN' ) then

    if ( action == 'GET' ) then
      value = pymin
    else if ( action == 'PRINT' ) then
      write ( *, '(a,i9)' ) 'PostScript minimum Y point, PYMIN = ', pymin
    else if ( action == 'SET' ) then
      pymin = value
    else if ( action == 'POP' ) then
      call i4_swap ( pymin, value )
    else if ( action == 'PUSH' ) then
      call i4_swap ( pymin, value )
    end if

  else if ( variable == 'STATE' ) then

    if ( action == 'GET' ) then
      value = state
    else if ( action == 'PRINT' ) then
      write ( *, '(a,i9)' ) 'Current internal state, STATE = ', state
    else if ( action == 'SET' ) then
      state = value
    else if ( action == 'POP' ) then
      call i4_swap ( state, value )
    else if ( action == 'PUSH' ) then
      call i4_swap ( state, value )
    end if

  else if ( variable == 'UNIT' ) then

    if ( action == 'GET' ) then
      value = unit
    else if ( action == 'PRINT' ) then
      write ( *, '(a,i9)' ) 'Current FORTRAN unit, UNIT = ', unit
    else if ( action == 'SET' ) then
      unit = value
    else if ( action == 'POP' ) then
      call i4_swap ( unit, value )
    else if ( action == 'PUSH' ) then
      call i4_swap ( unit, value )
    end if

  end if

  return
end
subroutine ps_setting_real ( action, variable, value )

!*****************************************************************************80
!
!! PS_SETTING_REAL sets, gets, or prints real internal PS_WRITE parameters.
!
!  Discussion:
!
!    I'd like a more sophisticated pop and push.
!
!    This routine has been revised to print an error message and stop
!    if the ACTION or VARIABLE is unexpected.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    Input, character ( len = * ) ACTION, is either:
!    'GET' to get the current value, or
!    'POP' to return the current value and set a new one;
!    'PRINT' to print the current value, or
!    'SET' to set the current value or
!    'PUSH' to set a new value and return the current one.
!
!    Input, character ( len = * ) VARIABLE, the variable to get or set:
!    'ALPHA', the scale factor from XY user space to PostScript points;
!    'FILL_BLUE', the intensity of the blue fill color, between 0.0 and 1.0.
!    'FILL_GREEN', the intensity of the green fill color, between 0.0 and 1.0.
!    'FILL_RED', the intensity of the red fill color, between 0.0 and 1.0.
!    'FONT_SIZE', the font size, in inches.
!    'LINE_BLUE', the blue component of the line color, between 0.0 and 1.0.
!    'LINE_GREEN', the green component of the line color, between 0.0 and 1.0.
!    'LINE_RED', the red component of the line color, between 0.0 and 1.0.
!    'XCUR', the current X location.
!    'XMAX', maximum X value of the data.
!    'XMIN', minimum X value of the data.
!    'YCUR', the current Y location.
!    'YMAX', maximum Y value of the data.
!    'YMIN', minimum Y value of the data.
!
!    Input/output, real ( kind = 8 ) VALUE.
!    If ACTION = 'GET', then VALUE is an output quantity, and is the
!    current internal value of the variable.
!
!    If ACTION = 'SET', then VALUE is an input quantity, and the
!    current internal value of the variable is set to this value.
!
!    If ACTION = 'PRINT', then VALUE is ignored.
!
  implicit none

  character ( len = * ) action
  real ( kind = 8 ), save :: alpha = 0.0D+00
  real ( kind = 8 ), save :: fill_blue = 0.7D+00
  real ( kind = 8 ), save :: fill_green = 0.7D+00
  real ( kind = 8 ), save :: fill_red = 0.7D+00
  real ( kind = 8 ), save :: font_size = 0.1D+00
  real ( kind = 8 ), save :: line_blue = 0.0D+00
  real ( kind = 8 ), save :: line_green = 0.0D+00
  real ( kind = 8 ), save :: line_red = 0.0D+00
  real ( kind = 8 ) value
  character ( len = * ) variable
  real ( kind = 8 ), save :: xcur = 0.0D+00
  real ( kind = 8 ), save :: xmax = 1.0D+00
  real ( kind = 8 ), save :: xmin = 0.0D+00
  real ( kind = 8 ), save :: ycur = 0.0D+00
  real ( kind = 8 ), save :: ymax = 0.0D+00
  real ( kind = 8 ), save :: ymin = 0.0D+00

  if ( variable == 'ALPHA' ) then

    if ( action == 'GET' ) then
      value = alpha
    else if ( action == 'PRINT' ) then
      write ( *, '(a,g14.6)' ) 'Scale factor from user to PS, ALPHA = ', alpha
    else if ( action == 'SET' ) then
      alpha = value
    else if ( action == 'POP' ) then
      call r8_swap ( alpha, value )
    else if ( action == 'PUSH' ) then
      call r8_swap ( alpha, value )
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PS_SETTING_REAL - Fatal error'
      write ( *, '(a)' ) '  Unexpected action!'
      stop
    end if

  else if ( variable == 'FILL_BLUE' ) then

    if ( action == 'GET' ) then
      value = fill_blue
    else if ( action == 'PRINT' ) then
      write ( *, '(a,g14.6)' ) 'Blue fill RGB value, FILL_BLUE = ', fill_blue
    else if ( action == 'SET' ) then
      fill_blue = value
    else if ( action == 'POP' ) then
      call r8_swap ( fill_blue, value )
    else if ( action == 'PUSH' ) then
      call r8_swap ( fill_blue, value )
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PS_SETTING_REAL - Fatal error'
      write ( *, '(a)' ) '  Unexpected action!'
      stop
    end if

  else if ( variable == 'FILL_GREEN' ) then

    if ( action == 'GET' ) then
      value = fill_green
    else if ( action == 'PRINT' ) then
      write ( *, '(a,g14.6)' ) 'Green fill RGB value, FILL_GREEN = ', fill_green
    else if ( action == 'SET' ) then
      fill_green = value
    else if ( action == 'POP' ) then
      call r8_swap ( fill_green, value )
    else if ( action == 'PUSH' ) then
      call r8_swap ( fill_green, value )
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PS_SETTING_REAL - Fatal error'
      write ( *, '(a)' ) '  Unexpected action!'
      stop
    end if

  else if ( variable == 'FILL_RED' ) then

    if ( action == 'GET' ) then
      value = fill_red
    else if ( action == 'PRINT' ) then
      write ( *, '(a,g14.6)' ) 'RED fill RGB value, FILL_RED = ', fill_red
    else if ( action == 'SET' ) then
      fill_red = value
    else if ( action == 'POP' ) then
      call r8_swap ( fill_red, value )
    else if ( action == 'PUSH' ) then
      call r8_swap ( fill_red, value )
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PS_SETTING_REAL - Fatal error'
      write ( *, '(a)' ) '  Unexpected action!'
      stop
    end if

  else if ( variable == 'FONT_SIZE' ) then

    if ( action == 'GET' ) then
      value = font_size
    else if ( action == 'PRINT' ) then
      write ( *, '(a,g14.6)' ) 'Font size, FONT_SIZE = ', font_size
    else if ( action == 'SET' ) then
      font_size = value
    else if ( action == 'POP' ) then
      call r8_swap ( font_size, value )
    else if ( action == 'PUSH' ) then
      call r8_swap ( font_size, value )
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PS_SETTING_REAL - Fatal error'
      write ( *, '(a)' ) '  Unexpected action!'
      stop
    end if

  else if ( variable == 'LINE_BLUE' ) then

    if ( action == 'GET' ) then
      value = line_blue
    else if ( action == 'PRINT' ) then
      write ( *, '(a,g14.6)' ) 'Blue line RGB value, LINE_BLUE = ', line_blue
    else if ( action == 'SET' ) then
      line_blue = value
    else if ( action == 'POP' ) then
      call r8_swap ( line_blue, value )
    else if ( action == 'PUSH' ) then
      call r8_swap ( line_blue, value )
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PS_SETTING_REAL - Fatal error'
      write ( *, '(a)' ) '  Unexpected action!'
      stop
    end if

  else if ( variable == 'LINE_GREEN' ) then

    if ( action == 'GET' ) then
      value = line_green
    else if ( action == 'PRINT' ) then
      write ( *, '(a,g14.6)' ) 'Green line RGB value, LINE_GREEN = ', line_green
    else if ( action == 'SET' ) then
      line_green = value
    else if ( action == 'POP' ) then
      call r8_swap ( line_green, value )
    else if ( action == 'PUSH' ) then
      call r8_swap ( line_green, value )
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PS_SETTING_REAL - Fatal error'
      write ( *, '(a)' ) '  Unexpected action!'
      stop
    end if

  else if ( variable == 'LINE_RED' ) then

    if ( action == 'GET' ) then
      value = line_red
    else if ( action == 'PRINT' ) then
      write ( *, '(a,g14.6)' ) 'Red line RGB value, LINE_RED = ', line_red
    else if ( action == 'SET' ) then
      line_red = value
    else if ( action == 'POP' ) then
      call r8_swap ( line_red, value )
    else if ( action == 'PUSH' ) then
      call r8_swap ( line_red, value )
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PS_SETTING_REAL - Fatal error'
      write ( *, '(a)' ) '  Unexpected action!'
      stop
    end if

  else if ( variable == 'XCUR' ) then

    if ( action == 'GET' ) then
      value = xcur
    else if ( action == 'PRINT' ) then
      write ( *, '(a,g14.6)' ) 'Current X location, XCUR = ', xcur
    else if ( action == 'SET' ) then
      xcur = value
    else if ( action == 'POP' ) then
      call r8_swap ( xcur, value )
    else if ( action == 'PUSH' ) then
      call r8_swap ( xcur, value )
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PS_SETTING_REAL - Fatal error'
      write ( *, '(a)' ) '  Unexpected action!'
      stop
    end if

  else if ( variable == 'XMAX' ) then

    if ( action == 'GET' ) then
      value = xmax
    else if ( action == 'PRINT' ) then
      write ( *, '(a,g14.6)' ) 'Maximum X value, XMAX = ', xmax
    else if ( action == 'SET' ) then
      xmax = value
    else if ( action == 'POP' ) then
      call r8_swap ( xmax, value )
    else if ( action == 'PUSH' ) then
      call r8_swap ( xmax, value )
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PS_SETTING_REAL - Fatal error'
      write ( *, '(a)' ) '  Unexpected action!'
      stop
    end if

  else if ( variable == 'XMIN' ) then

    if ( action == 'GET' ) then
      value = xmin
    else if ( action == 'PRINT' ) then
      write ( *, '(a,g14.6)' ) 'Minimum X value, XMIN = ', xmin
    else if ( action == 'SET' ) then
      xmin = value
    else if ( action == 'POP' ) then
      call r8_swap ( xmin, value )
    else if ( action == 'PUSH' ) then
      call r8_swap ( xmin, value )
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PS_SETTING_REAL - Fatal error'
      write ( *, '(a)' ) '  Unexpected action!'
      stop
    end if

  else if ( variable == 'YCUR' ) then

    if ( action == 'GET' ) then
      value = ycur
    else if ( action == 'PRINT' ) then
      write ( *, '(a,g14.6)' ) 'Current Y location, YCUR = ', ycur
    else if ( action == 'SET' ) then
      ycur = value
    else if ( action == 'POP' ) then
      call r8_swap ( ycur, value )
    else if ( action == 'PUSH' ) then
      call r8_swap ( ycur, value )
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PS_SETTING_REAL - Fatal error'
      write ( *, '(a)' ) '  Unexpected action!'
      stop
    end if

  else if ( variable == 'YMAX' ) then

    if ( action == 'GET' ) then
      value = ymax
    else if ( action == 'PRINT' ) then
      write ( *, '(a,g14.6)' ) 'Maximum Y value, YMAX = ', ymax
    else if ( action == 'SET' ) then
      ymax = value
    else if ( action == 'POP' ) then
      call r8_swap ( ymax, value )
    else if ( action == 'PUSH' ) then
      call r8_swap ( ymax, value )
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PS_SETTING_REAL - Fatal error'
      write ( *, '(a)' ) '  Unexpected action!'
      stop
    end if

  else if ( variable == 'YMIN' ) then

    if ( action == 'GET' ) then
      value = ymin
    else if ( action == 'PRINT' ) then
      write ( *, '(a,g14.6)' ) 'Minimum Y value, YMIN = ', ymin
    else if ( action == 'SET' ) then
      ymin = value
    else if ( action == 'POP' ) then
      call r8_swap ( ymin, value )
    else if ( action == 'PUSH' ) then
      call r8_swap ( ymin, value )
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PS_SETTING_REAL - Fatal error'
      write ( *, '(a)' ) '  Unexpected action!'
      stop
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_SETTING_REAL - Fatal error'
    write ( *, '(a)' ) '  Unexpected variable!'
    stop

  end if

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
!    01 May 2000
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
subroutine r8vec_sort_bubble_a ( n, a )

!*****************************************************************************80
!
!! R8VEC_SORT_BUBBLE_A ascending bubble sorts an R8VEC.
!
!  Discussion:
!
!    Bubble sort is simple to program, but inefficient.  It should not
!    be used for large arrays.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input/output, real ( kind = 8 ) A(N).
!    On input, an unsorted array.
!    On output, the array has been sorted.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) temp

  do i = 1, n-1
    do j = i+1, n
      if ( a(j) < a(i) ) then
        temp = a(i)
        a(i) = a(j)
        a(j) = temp
      end if
    end do
  end do

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
  character ( len = * ) s1
  character ( len = * ) s2

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

    ival = - ival
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

  character c
  logical ch_eqi
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
subroutine thin_index ( node_num, x, y, thin_factor, thin_num, thin_dex )

!*****************************************************************************80
!
!! THIN_INDEX determines thinning indices for a X, Y data.
!
!  Discussion:
!
!    A set of X, Y data is given, that is presumably, not too far off
!    from being on a rectangular grid.
!
!    The input value of THIN_FACTOR indicates by how much the data should
!    be thinned.
!
!    The X and Y ranges are computed, and only those data points are
!    retained for which both X and Y lie in an appropriate subrange.
!
!    For instance, a THIN_FACTOR of 2 would essentially save data
!    that lay in the black squares of a checkerboard.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 June
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, real ( kind = 8 ) X(NODE_NUM), Y(NODE_NUM), the X and Y coordinates
!    of the nodes.
!
!    Input, integer ( kind = 4 ) THIN_FACTOR, the thinning factor.
!
!    Output, integer ( kind = 4 ) THIN_NUM, the number of vectors remaining after
!    thinning.
!
!    Output, integer ( kind = 4 ) THIN_DEX(NODE_NUM), contains in (1:THIN_NUM) the
!    indices into X and Y of the vectors to be retained after thinning.
!
  implicit none

  integer ( kind = 4 ) node_num

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) thin_dex(node_num)
  integer ( kind = 4 ) thin_factor
  integer ( kind = 4 ) thin_num
  logical unique
  real ( kind = 8 ) x(node_num)
  integer ( kind = 4 ) x_bin
  real ( kind = 8 ) x_unique(node_num)
  integer ( kind = 4 ) x_unique_num
  real ( kind = 8 ) y(node_num)
  integer ( kind = 4 ) y_bin
  real ( kind = 8 ) y_unique(node_num)
  integer ( kind = 4 ) y_unique_num

  x_unique_num = 0

  do i = 1, node_num

    unique = .true.

    do j = 1, x_unique_num
      if ( x(i) == x_unique(j) ) then
        unique = .false.
        exit
      end if
    end do

    if ( unique ) then
      x_unique_num = x_unique_num + 1
      x_unique(x_unique_num) = x(i)
    end if

  end do

  call r8vec_sort_bubble_a ( x_unique_num, x_unique )

  y_unique_num = 0

  do i = 1, node_num

    unique = .true.

    do j = 1, y_unique_num
      if ( y(i) == y_unique(j) ) then
        unique = .false.
        exit
      end if
    end do

    if ( unique ) then
      y_unique_num = y_unique_num + 1
      y_unique(y_unique_num) = y(i)
    end if

  end do

  call r8vec_sort_bubble_a ( y_unique_num, y_unique )

  thin_num = 0

  do i = 1, node_num

    do j = 1, x_unique_num-1
      if ( x_unique(j) <= x(i) .and. x(i) <= x_unique(j+1) ) then
        x_bin = j
        exit
      end if
    end do

    do j = 1, y_unique_num-1
      if ( y_unique(j) <= y(i) .and. y(i) <= y_unique(j+1) ) then
        y_bin = j
        exit
      end if
    end do

    if ( mod ( y_bin, thin_factor ) == thin_factor / 2 .and. &
         mod ( x_bin, thin_factor ) == thin_factor / 2 ) then

      thin_num = thin_num + 1
      thin_dex(thin_num) = i

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
subroutine vector_field ( n, x, y, u, v, file_name, magnify )

!*****************************************************************************80
!
!! VECTOR_FIELD plots a vector field.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the location of the points at
!    which vector values are known.
!
!    Input, real ( kind = 8 ) U(N), V(N), the X and Y components of the
!    vector quantity.
!
!    Input, character ( len = * ) FILE_NAME, the name of the file.
!
!    Input, real ( kind = 8 ) MAGNIFY, a magnification factor for the vectors.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) area
  real ( kind = 8 ) blue
  character ( len = * ) file_name
  real ( kind = 8 ) green
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iunit
  real ( kind = 8 ) mag
  real ( kind = 8 ) mag_max
  real ( kind = 8 ) mag_min
  real ( kind = 8 ) magnify
  real ( kind = 8 ) red
  real ( kind = 8 ) scale
  real ( kind = 8 ) u(n)
  real ( kind = 8 ) v(n)
  real ( kind = 8 ) vmax
  real ( kind = 8 ) x(n)
  integer ( kind = 4 ) :: x_ps_max = 612
  integer ( kind = 4 ) :: x_ps_min = 36
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmax2
  real ( kind = 8 ) xmin
  real ( kind = 8 ) xmin2
  real ( kind = 8 ) xscale
  real ( kind = 8 ) y(n)
  integer ( kind = 4 ) :: y_ps_max = 792
  integer ( kind = 4 ) :: y_ps_min = 36
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymax2
  real ( kind = 8 ) ymin
  real ( kind = 8 ) ymin2
!
!  Get the maximum vector length.
!
  vmax = maxval ( sqrt ( u(1:n)**2 + v(1:n)**2 ) )

  if ( vmax == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'VECTOR_FIELD - Error!'
    write ( *, '(a)' ) '  The vector field is null.'
    return
  end if
!
!  Open the file.
!
  call get_unit ( iunit )

  call ps_file_open ( file_name, iunit, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'VECTOR_FIELD'
    write ( *, '(a,i6)' ) '  File creation error ', ierror
    return
  end if

  call eps_file_head ( file_name, x_ps_min, y_ps_min, x_ps_max, y_ps_max )
!
!  Define the size of the page.
!
  xmin = minval ( x(1:n) )
  ymin = minval ( y(1:n) )

  xmax = maxval ( x(1:n) )
  ymax = maxval ( y(1:n) )

  area = ( xmax - xmin ) * ( ymax - ymin )

  xscale = sqrt ( area / real ( n, kind = 8 ) )
  scale = magnify * sqrt ( area / real ( n, kind = 8 ) ) / vmax
!
!  Move the boundaries out far enough to catch a typical vector.
!
  xmax = xmax + xscale
  xmin = xmin - xscale
  ymax = ymax + xscale
  ymin = ymin - xscale
!
!  Adjust symmetrically so that X and Y ranges are equal.
!
  if ( ymax - ymin < xmax - xmin ) then
    xmax2 = xmax
    xmin2 = xmin
    ymax2 = ymax + 0.5D+00 * ( ( xmax - xmin ) - ( ymax - ymin ) )
    ymin2 = ymin - 0.5D+00 * ( ( xmax - xmin ) - ( ymax - ymin ) )
  else
    xmax2 = xmax + 0.5D+00 * ( ( ymax - ymin ) - ( xmax - xmin ) )
    xmin2 = xmin - 0.5D+00 * ( ( ymax - ymin ) - ( xmax - xmin ) )
    ymax2 = ymax
    ymin2 = ymin
  end if

  call ps_page_head ( xmin2, ymin2, xmax2, ymax2 )
!
!  Mark the points.
!
  call ps_mark_disks ( n, x, y )
!
!  Compute the velocity magnitude range.
!
  mag_min =  huge ( mag_min )
  mag_max = -huge ( mag_max )

  do i = 1, n
    mag = sqrt ( u(i)**2 + v(i)**2 )
    mag_min = min ( mag_min, mag )
    mag_max = max ( mag_max, mag )
  end do

  if ( mag_max == mag_min ) then

    mag = sqrt ( u(1)**2 + v(1)**2 )

    if ( mag /= 0.0D+00 ) then
      mag_min = mag / 2.0D+00
      mag_max = mag + mag_min
    else
      mag_min = 0.0D+00
      mag_max = 1.0D+00
    end if

  end if
!
!  Draw the vector field.
!
  do i = 1, n

    mag = sqrt ( u(i)**2 + v(i)**2 )

    red = ( mag_max - mag ) / ( mag_max - mag_min )
    green = 0.0D+00
    blue = ( mag - mag_min ) / ( mag_max - mag_min )

    call ps_color_line_set ( red, green, blue )

    call ps_arrow ( x(i), y(i), x(i) + scale * u(i), y(i) + scale * v(i) )

  end do
!
!  Close up the page and the file.
!
  call ps_page_tail ( )

  call eps_file_tail ( )

  call ps_file_close ( iunit )

  return
end
