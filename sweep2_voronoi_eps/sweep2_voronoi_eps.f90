program main

!*****************************************************************************80
!
!! MAIN is the main program for SWEEP2_VORONOI_EPS.
!
!  Discussion:
!
!    SWEEP2_VORONOI_EPS makes an EPS plot of a SWEEP2 Voronoi diagram.
!
!    Two big issues remain:
!
!    * Edges that include a vertex at infinity are not drawn now.
!      I need to figure out how to handle them.
!
!    * The test dataset includes vertices that are far outside the region.
!      This makes the actual plot a mess; there are two long thin spikes,
!      and then a knot of black in the center.  So I would like to include
!      an option to do clipping; or perhaps to make the plot size
!      dependent only on the points, and not the vertices (though I
!      dislike this idea.)
!
!  Invocation:
!
!    sweep2_voronoi_eps  voronoi.txt  plot.eps
!
!  Files:
!
!    The input file "voronoi.txt" contains the description of a Voronoi
!    diagram as computed by the SWEEP2 program.
!
!    The output file "plot.eps" is an Encapsulated PostScript file
!    containing an image of the Voronoi diagram.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 October 2003
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: edge
  integer ( kind = 4 ) edge_num
  character ( len = 100 ) eps_filename
  integer ( kind = 4 ) j
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: line
  integer ( kind = 4 ) line_num
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: point
  integer ( kind = 4 ) point_num
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: vertex
  integer ( kind = 4 ) vertex_num
  character ( len = 100 ) voronoi_filename

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SWEEP2_VORONOI_EPS:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Read'
  write ( *, '(a)' ) '    "voronoi.txt", a file, created by SWEEP2,'
  write ( *, '(a)' ) '    containing a description of a Voronoi diagram;'
  write ( *, '(a)' ) '  Create:'
  write ( *, '(a)' ) '    "plot.eps", an Encapsulated PostScript file'
  write ( *, '(a)' ) '    containing an image of the Voronoi diagram.'

  call get_filenames ( voronoi_filename, eps_filename )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  User specification:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    Voronoi file:  ' // trim ( voronoi_filename )
  write ( *, '(a)' ) '    EPS file:      ' // trim ( eps_filename )
!
!  Read the Voronoi information.
!
  call voronoi_header_read ( voronoi_filename, edge_num, line_num, &
    point_num, vertex_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '    Edges:    ', edge_num
  write ( *, '(a,i8)' ) '    Lines:    ', line_num
  write ( *, '(a,i8)' ) '    Points:   ', point_num
  write ( *, '(a,i8)' ) '    Vertices: ', vertex_num

  if ( point_num <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SWEEP2_VORONOI_EPS - Fatal error!'
    write ( *, '(a)' ) '  The Voronoi file does not contain a positive number'
    write ( *, '(a)' ) '  of points.'
    write ( *, '(a)' ) '  Abnormal end of execution.'
    stop
  end if
!
!  Read the data.
!
  allocate ( point(1:2,1:point_num) )
  allocate ( vertex(1:2,1:vertex_num) )
  allocate ( edge(1:3,1:edge_num) )
  allocate ( line(1:3,1:line_num) )

  call voronoi_data_read ( voronoi_filename, edge_num, line_num, &
    point_num, vertex_num, edge, line, point, vertex )

  write ( *, * ) ' '
  do j = 1, point_num
    write ( *, '(2x,f10.4,2x,f10.4)' ) point(1:2,j)
  end do
  write ( *, * ) ' '
  do j = 1, vertex_num
    write ( *, '(2x,f10.4,2x,f10.4)' ) vertex(1:2,j)
  end do
!
!  Plot the Voronoi diagram.
!
  call voronoi_eps ( eps_filename, edge_num, line_num, &
    point_num, vertex_num, edge, line, point, vertex )
!
!  Free memory.
!
  deallocate ( edge )
  deallocate ( line )
  deallocate ( point )
  deallocate ( vertex )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SWEEP2_VORONOI_EPS:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine box_clip_line_2d ( xmin, ymin, xmax, ymax, x1, y1, x2, y2, x3, y3, &
  x4, y4, ival )

!*****************************************************************************80
!
!! BOX_CLIP_LINE_2D uses a box to clip a line segment in 2D.
!
!  Discussion:
!
!    The box is assumed to be a rectangle with sides aligned on coordinate
!    axes.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 June 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XMIN, YMIN, XMAX, YMAX, the minimum and maximum
!    X and Y values, which define the box.
!
!    Input, real ( kind = 8 ) X1, Y1, X2, Y2, the coordinates of the
!    endpoints of the line segment.
!
!    Output, real ( kind = 8 ) X3, Y3, X4, Y4, the clipped coordinates.
!
!    Output, integer ( kind = 4 ) IVAL:
!    -1, no part of the line segment is within the box.
!     0, no clipping was necessary.  The line segment is entirely within
!        the box.
!     1, (X1,Y1) was clipped.
!     2, (X2,Y2) was clipped.
!     3, (X1,Y1) and (X2,Y2) were clipped.
!
  implicit none

  integer ( kind = 4 ) ival
  logical l1
  logical l2
  real ( kind = 8 ) x
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ) x4
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) y
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) y3
  real ( kind = 8 ) y4
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin

  l1 = .false.
  l2 = .false.

  x3 = x1
  y3 = y1
  x4 = x2
  y4 = y2
!
!  Require that XMIN <= X.
!
  if ( x3 < xmin .and. x4 < xmin ) then
    ival = -1
    return
  end if

  if ( x3 < xmin .and. xmin <= x4 ) then
    x = xmin
    y = y3 + ( y4 - y3 ) * ( x - x3 ) / ( x4 - x3 )
    x3 = x
    y3 = y
    l1 = .true.
  else if ( xmin <= x3 .and. x4 < xmin ) then
    x = xmin
    y = y3 + ( y4 - y3 ) * ( x - x3 ) / ( x4 - x3 )
    x4 = x
    y4 = y
    l2 = .true.
  end if
!
!  Require that X <= XMAX.
!
  if ( xmax < x3 .and. xmax < x4 ) then
    ival = -1
    return
  end if

  if ( xmax < x3 .and. x4 <= xmax ) then
    x = xmax
    y = y3 + ( y4 - y3 ) * ( x - x3 ) / ( x4 - x3 )
    x3 = x
    y3 = y
    l1 = .true.
  else if ( x3 <= xmax .and. xmax < x4 ) then
    x = xmax
    y = y3 + ( y4 - y3 ) * ( x - x3 ) / ( x4 - x3 )
    x4 = x
    y4 = y
    l2 = .true.
  end if
!
!  Require that YMIN <= Y.
!
  if ( y3 < ymin .and. y4 < ymin ) then
    ival = -1
    return
  end if

  if ( y3 < ymin .and. ymin <= y4 ) then
    y = ymin
    x = x3 + ( x4 - x3 ) * ( y - y3 ) / ( y4 - y3 )
    y3 = y
    x3 = x
    l1 = .true.
  else if ( ymin <= y3 .and. y4 < ymin ) then
    y = ymin
    x = x3 + ( x4 - x3 ) * ( y - y3 ) / ( y4 - y3 )
    y4 = y
    x4 = x
    l2 = .true.
  end if
!
!  Require that Y <= YMAX.
!
  if ( ymax < y3 .and. ymax < y4 ) then
    ival = -1
    return
  end if

  if ( ymax < y3 .and. y4 <= ymax ) then
    y = ymax
    x = x3 + ( x4 - x3 ) * ( y - y3 ) / ( y4 - y3 )
    y3 = y
    x3 = x
    l1 = .true.
  else if ( y3 <= ymax .and. ymax < y4 ) then
    y = ymax
    x = x3 + ( x4 - x3 ) * ( y - y3 ) / ( y4 - y3 )
    y4 = y
    x4 = x
    l2 = .true.
  end if

  ival = 0

  if ( l1 ) then
    ival = ival + 1
  end if

  if ( l2 ) then
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
  write ( unit, '(a)' )     '%%Creator: sweep2_voronoi_eps.f90'
  write ( unit, '(a)' )     '%%Title: ' // trim ( file_name )
  write ( unit, '(a)' )     '%%CreationDate: '// trim ( date )
  write ( unit, '(a)' )     '%%Pages: 1'
  write ( unit, '(a,4i8)' ) '%%BoundingBox:', &
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

  if ( 1 < num_pages ) then
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
subroutine get_filenames ( voronoi_filename, eps_filename )

!*****************************************************************************80
!
!! GET_FILENAMES gets the filenames from the command line or the user.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 October 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) VORONOI_FILENAME, EPS_FILENAME,
!    the names of the SWEEP2 Voronoi diagram file, and
!    Encapsulated PostScript graphics file (to be created), as
!    specified by the user on the command line or in response to
!    prompts from this routine.
!
  implicit none

  character ( len = 100 ) eps_filename
  integer ( kind = 4 ) iarg
  integer ( kind = 4 ) iargc
! integer ierror
! integer ilen
  integer ( kind = 4 ) input_unit
! integer ipxfargc
  integer ( kind = 4 ) numarg
  character ( len = 100 ) voronoi_filename
!
!  Count the number of command line arguments.
!
!  New style:
!
! numarg = ipxfargc ( )
!
!  Old style:
!
  numarg = iargc ( )
!
!  Get the point file name.
!
  if ( 1 <= numarg ) then

    iarg = 1
!
!  New style:
!
!   call pxfgetarg ( iarg, voronoi_filename, ilen, ierror )
!
!   if ( ierror /= 0 ) then
!     write ( *, '(a)' ) ' '
!     write ( *, '(a)' ) 'GET_FILENAMES - Fatal error!'
!     write ( *, '(a)' ) '  Could not read commandline argument.'
!     stop
!   end if
!
!  Old style:
!
    call getarg ( iarg, voronoi_filename )

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'What is the name of the Voronoi diagram file?'
    read ( *, '(a)' ) voronoi_filename
    if ( voronoi_filename == ' ' ) then
      stop
    end if

  end if
!
!  Get the eps file name.
!
  if ( 2 <= numarg ) then

    iarg = 2
!
!  New style:
!
!   call pxfgetarg ( iarg, eps_filename, ilen, ierror )
!
!   if ( ierror /= 0 ) then
!     write ( *, '(a)' ) ' '
!     write ( *, '(a)' ) 'GET_FILENAMES - Fatal error!'
!     write ( *, '(a)' ) '  Could not read commandline argument.'
!     stop
!   end if
!
!  Old style:
!
    call getarg ( iarg, eps_filename )

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'What is the name of the output EPS file?'
    read ( *, '(a)' ) eps_filename
    if ( eps_filename == ' ' ) then
      stop
    end if

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
subroutine i4_swap ( i, j )

!*****************************************************************************80
!
!! I4_SWAP swaps two I4's.
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
subroutine ps_clip ( npoint, x, y )

!*****************************************************************************80
!
!! PS_CLIP defines a clipping polygon.
!
!  Discussion:
!
!    Use this routine if you want to draw more than you display.
!    A clipping polygon allows you to define points and lines
!    that lie (partially) outside of the polygon, but only display
!    the portions within the polygon
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 May 2001
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
!    Input, integer ( kind = 4 ) NPOINT, the number of points in the clipping polygon.
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
  if ( npoint < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_CLIP - Warning!'
    write ( *, '(a)' ) '  Clipping polygon has too few sides.'
    write ( *, '(a,i9)' ) '  NPOINT = ', npoint
    return
  end if
!
!  Determine if the PostScript state is acceptable.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state /= 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_CLIP - Fatal error!'
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
  call ps_comment ( 'Define a clipping polygon' )

  px = plotxmin2 + nint ( alpha * ( x(1) - xmin ) )
  py = plotymin2 + nint ( alpha * ( y(1) - ymin ) )
  write ( unit, '(2i8,a)' ) px, py, ' moveto'

  do i = 2, npoint
    px = plotxmin2 + nint ( alpha * ( x(i) - xmin ) )
    py = plotymin2 + nint ( alpha * ( y(i) - ymin ) )
    write ( unit, '(2i8,a)' ) px, py, ' lineto'
  end do
!
!  Add the final extra segment to the initial point.
!
  px = plotxmin2 + nint ( alpha * ( x(1) - xmin ) )
  py = plotymin2 + nint ( alpha * ( y(1) - ymin ) )
  write ( unit, '(2i8,a)' ) px, py, ' lineto'
!
!  Fill the polygon.
!
  write ( unit, '(a)' ) 'clip newpath'

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
  real ( kind = 8 ) b_stack(nstack)
  real ( kind = 8 ) g
  real ( kind = 8 ) g_old
  real ( kind = 8 ) g_stack(nstack)
  integer ( kind = 4 ), save :: istack = 0
  real ( kind = 8 ) r
  real ( kind = 8 ) r_old
  real ( kind = 8 ) r_stack(nstack)
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
!    24 April 2001
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

  if ( state < 1 .or. 4 < state ) then
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
subroutine ps_font_size ( font_size )

!*****************************************************************************80
!
!! PS_FONT_SIZE sets the font size.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 July 1999
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
!    Input, real ( kind = 8 ) FONT_SIZE, the font size, in inches.
!
  implicit none

  real ( kind = 8 ) font_size
  integer ( kind = 4 ) state
  integer ( kind = 4 ) unit
!
!  Determine if the PostScript state is acceptable.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state /= 2 .and. state /= 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_FONT_SIZE - Fatal error!'
    write ( *, '(a,i9)' ) '  PostScript state is ', state
    write ( *, '(a)' ) '  PostScript state 2 or 3 is required.'
    return
  end if

  call ps_setting_int ( 'GET', 'UNIT', unit )

  write ( unit, '(a)' ) '/Times-Roman findfont'
  write ( unit, '(f8.3, a)' ) font_size, ' inch scalefont'
  write ( unit, '(a)' ) 'setfont'

  call ps_setting_real ( 'SET', 'FONT_SIZE', font_size )

  return
end
subroutine ps_label ( string )

!*****************************************************************************80
!
!! PS_LABEL prints a label at the current position.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 July 1999
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
!    Input, character ( len = * ) STRING, the string to be printed.
!
  implicit none

  character ( len = * ) string
  integer ( kind = 4 ) unit

  if ( len_trim ( string ) <= 0 ) then
    return
  end if

  call ps_setting_int ( 'GET', 'UNIT', unit )

  write ( unit, '(a)' ) '(' // trim ( string ) // ') show'

  return
end
subroutine ps_line ( x1, y1, x2, y2 )

!*****************************************************************************80
!
!! PS_LINE draws a line segment from (X1,Y1) to (X2,Y2).
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
!    19 June 2001
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
!    Input, real ( kind = 8 ) X1, Y1, the starting point of the line segment.
!
!    Input, real ( kind = 8 ) X2, Y2, the ending point of the line segment.
!
  implicit none

  real ( kind = 8 ) alpha
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) plotxmin2
  integer ( kind = 4 ) plotymin2
  integer ( kind = 4 ) px
  integer ( kind = 4 ) py
  integer ( kind = 4 ) state
  integer ( kind = 4 ) unit
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ) x4
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) y3
  real ( kind = 8 ) y4
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin
!
!  Determine if the PostScript state is acceptable.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state /= 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_LINE - Fatal error!'
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
   call box_clip_line_2d ( xmin, ymin, xmax, ymax, x1, y1, x2, y2, x3, y3, &
     x4, y4, ival )

   if ( ival < 0 ) then
     return
   end if
!
!  Draw line.
!
  write ( unit, '(a)' ) 'newpath'

  px = plotxmin2 + nint ( alpha * ( x3 - xmin ) )
  py = plotymin2 + nint ( alpha * ( y3 - ymin ) )
  write ( unit, '(2i8,a)' ) px, py, ' moveto'

  px = plotxmin2 + nint ( alpha * ( x4 - xmin ) )
  py = plotymin2 + nint ( alpha * ( y4 - ymin ) )
  write ( unit, '(2i8,a)' ) px, py, ' lineto stroke'

  call ps_setting_real ( 'SET', 'XCUR', x2 )
  call ps_setting_real ( 'SET', 'YCUR', y2 )

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
  write ( unit, '(2i8,a)' ) px, py, ' moveto'

  do i = 2, npoint
    px = plotxmin2 + nint ( alpha * ( x(i) - xmin ) )
    py = plotymin2 + nint ( alpha * ( y(i) - ymin ) )
    write ( unit, '(2i8,a)' ) px, py, ' lineto'
  end do
!
!  Add the final extra segment to the initial point.
!
  px = plotxmin2 + nint ( alpha * ( x(1) - xmin ) )
  py = plotymin2 + nint ( alpha * ( y(1) - ymin ) )
  write ( unit, '(2i8,a)' ) px, py, ' lineto'
!
!  Draw the line.
!
  write ( unit, '(a)' ) 'stroke'

  call ps_setting_real ( 'SET', 'XCUR', x(1) )
  call ps_setting_real ( 'SET', 'YCUR', y(1) )

  return
end
subroutine ps_line_width ( line_width )

!*****************************************************************************80
!
!! PS_LINE_WIDTH sets the line width.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 July 1999
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
!    Input, integer ( kind = 4 ) LINE_WIDTH, the line width.
!    0 is a valid input, and usually produces the thinnest possible line.
!    1 is a more usual line, 2 is thicker, and so on.
!
  implicit none

  integer ( kind = 4 ) line_width
  integer ( kind = 4 ) state
  integer ( kind = 4 ) unit
!
!  Determine if the PostScript state is acceptable.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state /= 2 .and. state /= 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_LINE_WIDTH - Fatal error!'
    write ( *, '(a,i9)' ) '  PostScript state is ', state
    write ( *, '(a)' ) '  PostScript state 2 or 3 is required.'
    return
  end if

  call ps_setting_int ( 'GET', 'UNIT', unit )

  write ( unit, '(i8,a)' ) line_width, ' setlinewidth'

  call ps_setting_int ( 'SET', 'LINE_WIDTH', line_width )

  return
end
subroutine ps_mark_disk ( x, y )

!*****************************************************************************80
!
!! PS_MARK_DISK marks a point with a small filled disk.
!
!  Discussion:
!
!    The current point is set to the center of the disk.
!
!    The circle is drawn with the current RGB fill colors.
!
!    The circle is drawn the current marker size.
!
!    Points outside the region are not marked.
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
!    Henry McGilton and Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, Y, the coordinates of the point to mark.
!
  implicit none

  real ( kind = 8 ) alpha
  integer ( kind = 4 ) marker_size
  integer ( kind = 4 ) plotxmin2
  integer ( kind = 4 ) plotymin2
  logical point_inside_box_2d
  integer ( kind = 4 ) pxcen
  integer ( kind = 4 ) pycen
  integer ( kind = 4 ) state
  integer ( kind = 4 ) unit
  real ( kind = 8 ) x
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) y
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin
!
!  Determine if the PostScript state is acceptable.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state /= 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_MARK_DISK - Fatal error!'
    write ( *, '(a,i9)' ) '  PostScript state is ', state
    write ( *, '(a)' ) '  PostScript state 3 is required.'
    return
  end if
!
!  Get settings.
!
  call ps_setting_int ( 'GET', 'MARKER_SIZE', marker_size )
  call ps_setting_int ( 'GET', 'PXMIN', plotxmin2 )
  call ps_setting_int ( 'GET', 'PYMIN', plotymin2 )
  call ps_setting_int ( 'GET', 'UNIT', unit )

  call ps_setting_real ( 'GET', 'ALPHA', alpha )
  call ps_setting_real ( 'GET', 'XMIN', xmin )
  call ps_setting_real ( 'GET', 'XMAX', xmax )
  call ps_setting_real ( 'GET', 'YMIN', ymin )
  call ps_setting_real ( 'GET', 'YMAX', ymax )
!
!  If the point is outside the plot box, don't draw it.
!
  if ( .not. point_inside_box_2d ( xmin, ymin, xmax, ymax, x, y ) ) then
    return
  end if

  write ( unit, '(a)' ) 'newpath'

  pxcen = plotxmin2 + nint ( alpha * ( x - xmin ) )
  pycen = plotymin2 + nint ( alpha * ( y - ymin ) )

  write ( unit, '(3i8,a)' ) pxcen, pycen, marker_size, &
    ' 0 360 arc closepath fill'

  call ps_setting_real ( 'SET', 'XCUR', x )
  call ps_setting_real ( 'SET', 'YCUR', y )

  return
end
subroutine ps_mark_point ( x, y )

!*****************************************************************************80
!
!! PS_MARK_POINT marks a point with a tiny point.
!
!  Discussion:
!
!    The current point is set to the point.
!
!    The point is drawn with the current RGB line colors.
!
!    If the point is outside the region, the command is ignored.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 May 2001
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
!    Input, real ( kind = 8 ) X, Y, the coordinates of the point to mark.
!
  implicit none

  real ( kind = 8 ) alpha
  integer ( kind = 4 ) marker_size
  integer ( kind = 4 ) plotxmin2
  integer ( kind = 4 ) plotymin2
  logical point_inside_box_2d
  integer ( kind = 4 ) pxcen
  integer ( kind = 4 ) pycen
  integer ( kind = 4 ) state
  integer ( kind = 4 ) unit
  real ( kind = 8 ) x
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) y
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin
!
!  Determine if the PostScript state is acceptable.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state /= 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_MARK_POINT - Fatal error!'
    write ( *, '(a,i9)' ) '  PostScript state is ', state
    write ( *, '(a)' ) '  PostScript state 3 is required.'
    return
  end if
!
!  Get settings.
!
  call ps_setting_int ( 'GET', 'MARKER_SIZE', marker_size )
  call ps_setting_int ( 'GET', 'PXMIN', plotxmin2 )
  call ps_setting_int ( 'GET', 'PYMIN', plotymin2 )
  call ps_setting_int ( 'GET', 'UNIT', unit )

  call ps_setting_real ( 'GET', 'ALPHA', alpha )
  call ps_setting_real ( 'GET', 'XMIN', xmin )
  call ps_setting_real ( 'GET', 'XMAX', xmax )
  call ps_setting_real ( 'GET', 'YMIN', ymin )
  call ps_setting_real ( 'GET', 'YMAX', ymax )
!
!  If the point is outside the plot box, don't draw it.
!
  if ( .not. point_inside_box_2d ( xmin, ymin, xmax, ymax, x, y ) ) then
    return
  end if

  call ps_comment ( 'Draw a point' )

  write ( unit, '(a)' ) 'newpath'

  pxcen = plotxmin2 + nint ( alpha * ( x - xmin ) )
  pycen = plotymin2 + nint ( alpha * ( y - ymin ) )

  write ( unit, '(2i8,a)' ) pxcen, pycen, ' moveto'
  write ( unit, '(2i8,a)' ) pxcen+1, pycen, ' lineto'
  write ( unit, '(a)' ) 'stroke'

  call ps_setting_real ( 'SET', 'XCUR', x )
  call ps_setting_real ( 'SET', 'YCUR', y )

  return
end
subroutine ps_moveto ( x, y )

!*****************************************************************************80
!
!! PS_MOVETO "moves to" a new point, which becomes the current point.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 July 1999
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
!    Input, real ( kind = 8 ) X, Y, the X and Y components of the current point.
!
  implicit none

  real ( kind = 8 ) alpha
  integer ( kind = 4 ) plotxmin2
  integer ( kind = 4 ) plotymin2
  integer ( kind = 4 ) px
  integer ( kind = 4 ) py
  integer ( kind = 4 ) state
  integer ( kind = 4 ) unit
  real ( kind = 8 ) x
  real ( kind = 8 ) xmin
  real ( kind = 8 ) y
  real ( kind = 8 ) ymin
!
!  Determine if the PostScript state is acceptable.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state /= 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_MOVETO - Fatal error!'
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
!  Move to the new point.
!
  px = plotxmin2 + nint ( alpha * ( x - xmin ) )
  py = plotymin2 + nint ( alpha * ( y - ymin ) )
  write ( unit, '(2i8,a)' ) px, py, ' moveto'

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

  write ( unit, '(a,i8,i8)' ) '%%Page: ', num_pages, num_pages
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
    end if

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
subroutine s_to_r8 ( s, r, ierror, length )

!*****************************************************************************80
!
!! S_TO_R8 reads an R8 from a string.
!
!  Discussion:
!
!    This routine will read as many characters as possible until it reaches
!    the end of the string, or encounters a character which cannot be
!    part of the number.
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
!    characters can be read to form a legal number.  Blanks,
!    commas, or other nonnumeric data will, in particular,
!    cause the conversion to halt.
!
!    Output, real ( kind = 8 ) R, the value that was read from the string.
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
!    Output, integer ( kind = 4 ) LENGTH, the number of characters read from
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
  integer ( kind = 4 ) length
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
  length = - 1
  isgn = 1
  rtop = 0.0D+00
  rbot = 1.0D+00
  jsgn = 1
  jtop = 0
  jbot = 1
  ihave = 1
  iterm = 0

  do

    length = length + 1
    c = s(length+1:length+1)
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
        length = length + 1
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
    if ( iterm == 1 .or. nchar <= length+1 ) then
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
!    15 March 2003
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

  character ( len = 40 ) string

  call timestring ( string )

  write ( *, '(a)' ) trim ( string )

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
subroutine voronoi_eps ( eps_filename, edge_num, line_num, point_num, &
  vertex_num, edge, line, point, vertex )

!*****************************************************************************80
!
!! VORONOI_EPS plots a Voronoi diagram.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 October 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) EPS_FILENAME, the name of the output file.
!
!    Input, integer ( kind = 4 ) EDGE_NUM, LINE_NUM, POINT_NUM, VERTEX_NUM, the number
!    of edge, line, point, and vertex records in the file.
!
!    Input, integer ( kind = 4 ) EDGE(3,EDGE_NUM), for each Voronoi edge, the index
!    of a line, and the indices of the two end vertices.
!
!    Input, real ( kind = 8 ) LINE(3,LINE_NUM), for each line, the coefficients
!    of the equation A*X+B*Y=C.
!
!    Input, real ( kind = 8 ) POINT(2,POINT_NUM), the coordinates of each point.
!
!    Input, real ( kind = 8 ) VERTEX(2,POINT_NUM), the coordinates of
!    each vertex.
!
  implicit none

  integer ( kind = 4 ) edge_num
  integer ( kind = 4 ) line_num
  integer ( kind = 4 ) point_num
  integer ( kind = 4 ) vertex_num

  logical, parameter :: debug = .false.
  integer ( kind = 4 ) edge(3,edge_num)
  character ( len = * ) eps_filename
  integer ( kind = 4 ) eps_unit
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) l
  real ( kind = 8 ) line(3,line_num)
  real ( kind = 8 ) plot_max(2)
  real ( kind = 8 ) plot_min(2)
  real ( kind = 8 ) point(2,point_num)
  character ( len = 80 ) title
  integer ( kind = 4 ) v1
  integer ( kind = 4 ) v2
  real ( kind = 8 ) vertex(2,vertex_num)
  real ( kind = 8 ) width
  real ( kind = 8 ) xvec(4)
  real ( kind = 8 ) xx
  real ( kind = 8 ) yvec(4)
  real ( kind = 8 ) yy

  plot_min(1) = min ( minval ( point(1,1:point_num) ), &
                      minval ( vertex(1,1:vertex_num) ) )
  plot_max(1) = max ( maxval ( point(1,1:point_num) ), &
                      maxval ( vertex(1,1:vertex_num) ) )
  width = plot_max(1) - plot_min(1)
  plot_min(1) = plot_min(1) - 0.05D+00 * width
  plot_max(1) = plot_max(1) + 0.05D+00 * width

  plot_min(2) = min ( minval ( point(2,1:point_num) ), &
                      minval ( vertex(2,1:vertex_num) ) )
  plot_max(2) = max ( maxval ( point(2,1:point_num) ), &
                      maxval ( vertex(2,1:vertex_num) ) )

  width = plot_max(2) - plot_min(2)
  plot_min(2) = plot_min(2) - 0.05D+00 * width
  plot_max(2) = plot_max(2) + 0.05D+00 * width

  call get_unit ( eps_unit )

  call ps_file_open ( eps_filename, eps_unit, ierror )

  call eps_file_head ( eps_filename )

  call ps_page_head ( plot_min(1), plot_min(2), plot_max(1), plot_max(2) )
!
!  Initialize the line width to 1.
!
  call ps_line_width ( 1 )
!
!  Print title, if requested.
!
  title = 'Voronoi diagram by SWEEP2'

  if ( 0 < len_trim ( title ) ) then
    call ps_font_size ( 0.30D+00 )
    call ps_moveto ( plot_min(1), plot_max(2) )
    call ps_label ( title )
  end if
!
!  Define a PostScript clipping box.
!
  xvec(1:4) = (/ plot_min(1), plot_max(1), plot_max(1), plot_min(1) /)
  yvec(1:4) = (/ plot_min(2), plot_min(2), plot_max(2), plot_max(2) /)

  call ps_clip ( 4, xvec, yvec )

  xx = plot_min(1)
  yy = plot_min(2)
  call ps_moveto ( xx, yy )

  call ps_comment ( 'Mark the points.' )

  if ( point_num <= 350 ) then
    do i = 1, point_num
      call ps_mark_disk ( point(1,i), point(2,i) )
    end do
  else
    do i = 1, point_num
      call ps_mark_point ( point(1,i), point(2,i) )
    end do
  end if

  call ps_comment ( 'Mark the Voronoi vertices' )

  if ( vertex_num <= 350 ) then
    do i = 1, vertex_num
      call ps_mark_disk ( vertex(1,i), vertex(2,i) )
    end do
  else
    do i = 1, vertex_num
      call ps_mark_point ( vertex(1,i), vertex(2,i) )
    end do
  end if
!
!  Increase the line width to "3".
!
  call ps_comment ( 'Draw the Voronoi edges.' )

  call ps_line_width ( 3 )

  do i = 1, edge_num

    l = edge(1,i) + 1
    v1 = edge(2,i) + 1
    v2 = edge(3,i) + 1

    if ( 1 <= v1 .and. 1 <= v2 ) then
      call ps_line (  vertex(1,v1), vertex(2,v1), vertex(1,v2), vertex(2,v2) )
    end if

  end do

  call ps_page_tail ( )

  call eps_file_tail ( )

  call ps_file_close ( eps_unit )

  return
end
subroutine voronoi_data_read ( voronoi_filename, edge_num, line_num, &
  point_num, vertex_num, edge, line, point, vertex )

!*****************************************************************************80
!
!! VORONOI_DATA_READ reads the data from a Voronoi file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) VORONOI_FILENAME, the SWEEP2 Voronoi
!    filename.
!
!    Input, integer ( kind = 4 ) EDGE_NUM, LINE_NUM, POINT_NUM, VERTEX_NUM, the number
!    of edge, line, point, and vertex records in the file.
!
!    Output, integer ( kind = 4 ) EDGE(3,EDGE_NUM), for each Voronoi edge, the index
!    of a line, and the indices of the two end vertices.
!
!    Output, real ( kind = 8 ) LINE(3,LINE_NUM), for each line, the coefficients
!    of the equation A*X+B*Y=C.
!
!    Output, real ( kind = 8 ) POINT(2,POINT_NUM), the coordinates
!    of each point.
!
!    Output, real ( kind = 8 ) VERTEX(2,POINT_NUM), the coordinates
!    of each vertex.
!
  implicit none

  integer ( kind = 4 ) edge_num
  integer ( kind = 4 ) line_num
  integer ( kind = 4 ) point_num
  integer ( kind = 4 ) vertex_num

  integer ( kind = 4 ) bad_num
  integer ( kind = 4 ) edge(3,edge_num)
  integer ( kind = 4 ) edge_num2
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) length
  real ( kind = 8 ) line(3,line_num)
  integer ( kind = 4 ) line_num2
  real ( kind = 8 ) point(2,point_num)
  integer ( kind = 4 ) point_num2
  character ( len = 255 ) string
  real ( kind = 8 ) vertex(2,vertex_num)
  integer ( kind = 4 ) vertex_num2
  character ( len = * ) voronoi_filename
  integer ( kind = 4 ) voronoi_unit
  character ( len = 255 ) word

  bad_num = 0
  edge_num2 = 0
  line_num2 = 0
  point_num2 = 0
  vertex_num2 = 0

  call get_unit ( voronoi_unit )

  open ( unit = voronoi_unit, file = voronoi_filename, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'VORONOI_DATA_READ - Fatal error!'
    write ( *, '(a)' ) '  Could not open the file "' // &
      trim ( voronoi_filename ) // '".'
    stop
  end if

  do

    read ( voronoi_unit, '(a)', iostat = ios ) string

    if ( ios /= 0 ) then
      exit
    end if

    if ( string(1:1) == '#' .or. len_trim ( string ) == 0 ) then
      cycle
    end if

    call word_extract ( string, word )

    if ( word == 'e' .or. word == 'E' ) then

      edge_num2 = edge_num2 + 1
      call word_extract ( string, word )
      call s_to_i4 ( word, edge(1,edge_num2), ierror, length )
      call word_extract ( string, word )
      call s_to_i4 ( word, edge(2,edge_num2), ierror, length )
      call word_extract ( string, word )
      call s_to_i4 ( word, edge(3,edge_num2), ierror, length )

    else if ( word == 'l' .or. word == 'L' ) then

      line_num2 = line_num2 + 1
      call word_extract ( string, word )
      call s_to_r8 ( word, line(1,line_num2), ierror, length )
      call word_extract ( string, word )
      call s_to_r8 ( word, line(2,line_num2), ierror, length )
      call word_extract ( string, word )
      call s_to_r8 ( word, line(3,line_num2), ierror, length )

    else if ( word == 's' .or. word == 'S' ) then

      point_num2 = point_num2 + 1
      call word_extract ( string, word )
      call s_to_r8 ( word, point(1,point_num2), ierror, length )
      call word_extract ( string, word )
      call s_to_r8 ( word, point(2,point_num2), ierror, length )

    else if ( word == 'v' .or. word == 'V' ) then

      vertex_num2 = vertex_num2 + 1
      call word_extract ( string, word )
      call s_to_r8 ( word, vertex(1,vertex_num2), ierror, length )
      call word_extract ( string, word )
      call s_to_r8 ( word, vertex(2,vertex_num2), ierror, length )

    else
      bad_num = bad_num + 1

    end if

  end do

  close ( unit = voronoi_unit )

  if ( bad_num /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'VORONOI_DATA_READ - Warning!'
    write ( *, '(a,i8)' ) '  Number of bad records was ', bad_num
  end if

  return
end
subroutine voronoi_header_read ( voronoi_filename, edge_num, line_num, &
  point_num, vertex_num )

!*****************************************************************************80
!
!! VORONOI_HEADER_READ reads the header information from a Voronoi file.
!
!  Discussion:
!
!    There is no header information in a SWEEP2 Voronoi file.  What this
!    routine does is read the "implicit" header information, that is,
!    it counts the number of records of each type.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 October 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) VORONOI_FILENAME, the SWEEP2 Voronoi
!    filename.
!
!    Output, integer ( kind = 4 ) EDGE_NUM, LINE_NUM, POINT_NUM, VERTEX_NUM, the number
!    of edge, line, point, and vertex records in the file.
!
  implicit none

  integer ( kind = 4 ) bad_num
  integer ( kind = 4 ) edge_num
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) line_num
  integer ( kind = 4 ) point_num
  character ( len = 255 ) string
  integer ( kind = 4 ) vertex_num
  character ( len = * ) voronoi_filename
  integer ( kind = 4 ) voronoi_unit
  character ( len = 255 ) word

  bad_num = 0
  edge_num = 0
  line_num = 0
  point_num = 0
  vertex_num = 0

  call get_unit ( voronoi_unit )

  open ( unit = voronoi_unit, file = voronoi_filename, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'VORONOI_HEADER_READ - Fatal error!'
    write ( *, '(a)' ) '  Could not open the file "' // &
      trim ( voronoi_filename ) // '".'
    stop
  end if

  do

    read ( voronoi_unit, '(a)', iostat = ios ) string

    if ( ios /= 0 ) then
      exit
    end if

    if ( string(1:1) == '#' .or. len_trim ( string ) == 0 ) then
      cycle
    end if

    call word_extract ( string, word )

    if ( word == 'e' .or. word == 'E' ) then
      edge_num = edge_num + 1
    else if ( word == 'l' .or. word == 'L' ) then
      line_num = line_num + 1
    else if ( word == 's' .or. word == 'S' ) then
      point_num = point_num + 1
    else if ( word == 'v' .or. word == 'V' ) then
      vertex_num = vertex_num + 1
    else
      bad_num = bad_num + 1
    end if

  end do

  close ( unit = voronoi_unit )

  if ( bad_num /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'VORONOI_HEADER_READ - Warning!'
    write ( *, '(a,i8)' ) '  Number of bad records was ', bad_num
  end if

  return
end
subroutine word_extract ( s, w )

!*****************************************************************************80
!
!! WORD_EXTRACT extracts the next word from a string.
!
!  Discussion:
!
!    A "word" is a string of characters terminated by a blank or
!    the end of the string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) S, the string.  On output, the first
!    word has been removed, and the remaining string has been shifted left.
!
!    Output, character ( len = * ) W, the leading word of the string.
!
  implicit none

  integer ( kind = 4 ) get1
  integer ( kind = 4 ) get2
  integer ( kind = 4 ) length
  character ( len = * ) s
  character ( len = * ) w

  w = ' '

  length = len_trim ( s )
!
!  Find the first nonblank.
!
  get1 = 0

  do

    get1 = get1 + 1

    if ( length < get1 ) then
      return
    end if

    if ( s(get1:get1) /= ' ' ) then
      exit
    end if

  end do
!
!  Look for the last contiguous nonblank.
!
  get2 = get1

  do

    if ( length <= get2 ) then
      exit
    end if

    if ( s(get2+1:get2+1) == ' ' ) then
      exit
    end if

    get2 = get2 + 1

  end do
!
!  Copy the word.
!
  w = s(get1:get2)
!
!  Shift the string.
!
  s(1:get2) = ' '
  s = adjustl ( s(get2+1:) )

  return
end
