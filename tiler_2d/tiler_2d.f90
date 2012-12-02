program main

!*****************************************************************************80
!
!! MAIN is the main program for TILER_2D.
!
!  Discussion:
!
!    TILER_2D coordinates the PostScript details and runs the test.
!
!    This program also requires the PS_WRITE routines to run.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gordon and Charles Hall,
!    Construction of Curvilinear Coordinate Systems and Application to
!    Mesh Generation,
!    International Journal of Numerical Methods in Engineering,
!    Volume 7, pages 461-477, 1973.
!
!    Joe Thompson, Bharat Soni, Nigel Weatherill,
!    Handbook of Grid Generation,
!    CRC Press,
!    1999.
!
  implicit none

  character ( len = 80 ) :: file_name = 'tiler_2d.ps'
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iunit
  real ( kind = 8 ) :: xmax =  5.0D+00
  real ( kind = 8 ) :: xmin = -5.0D+00
  real ( kind = 8 ) :: ymax =  5.0D+00
  real ( kind = 8 ) :: ymin =  0.0D+00

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TILER_2D'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Create a PostScript image of a 2D tiling.'

  iunit = 1

  call ps_file_open ( file_name, iunit, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TILER_2D'
    write ( *, '(a,i6)' ) '  File creation error ', ierror
    stop
  end if

  call ps_file_head ( file_name )

  call ps_page_head ( xmin, ymin, xmax, ymax )

  call box_tiler_2d ( iunit )

  call ps_page_tail ( )

  call ps_file_tail ( )

  call ps_file_close ( iunit )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TILER_2D'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine boundary_2d ( umin, vmin, umax, vmax, u, v, x, y )

!*****************************************************************************80
!
!! BOUNDARY_2D maps a point (U,V) to its corresponding (X,Y) coordinates.
!
!  Discussion:
!
!    The routine is only called to evaluate the (X,Y) coordinates of
!    a point (U,V) which lies on the boundary.  That is, one of
!    U and V will always be a maximum or minimum value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gordon and Charles Hall,
!    Construction of Curvilinear Coordinate Systems and Application to
!      Mesh Generation,
!    International Journal of Numerical Methods in Engineering,
!    Volume 7, pages 461-477, 1973.
!
!    Joe Thompson, Bharat Soni, Nigel Weatherill,
!    Handbook of Grid Generation,
!    CRC Press,
!    1999.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) UMIN, VMIN, UMAX, VMAX, the minimum and
!    maximum values of U and V.
!
!    Input, real ( kind = 8 ) U, V, the (U,V) coordinates of a point
!    that lies on the boundary of the region.
!
!    Output, real ( kind = 8 ) X, Y, the (X,Y) coordinates of the image
!    of the point.
!
  implicit none

  real ( kind = 8 ), parameter :: deg2rad = 3.14159265D+00 / 180.0D+00
  real ( kind = 8 ) theta
  real ( kind = 8 ) u
  real ( kind = 8 ) umax
  real ( kind = 8 ) umin
  real ( kind = 8 ) v
  real ( kind = 8 ) v2
  real ( kind = 8 ) vmax
  real ( kind = 8 ) vmin
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  theta = u * deg2rad

       if ( u == umin ) then
    x = v * cos ( theta )
    y = v * sin ( theta )
  else if ( u == umax ) then
    x = v * cos ( theta )
    y = v * sin ( theta )
  else if ( v == vmin ) then
    x = v * cos ( theta )
    y = v * sin ( theta )
  else if ( v == vmax ) then
    v2 = v + sin ( 3.0D+00 * theta )
    x = v2 * cos ( theta )
    y = v2 * sin ( theta )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BOUNDARY_2D - Fatal error!'
    write ( *, '(a)' ) '  Input (U,V) is not on the boundary.'
    stop
  end if

  return
end
subroutine box_tiler_2d ( iunit )

!*****************************************************************************80
!
!! BOX_TILER_2D selects one sub-box at a time to be tiled.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gordon and Charles Hall,
!    Construction of Curvilinear Coordinate Systems and Application to
!    Mesh Generation,
!    International Journal of Numerical Methods in Engineering,
!    Volume 7, pages 461-477, 1973.
!
!    Joe Thompson, Bharat Soni, Nigel Weatherill,
!    Handbook of Grid Generation,
!    CRC Press,
!    1999.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ( kind = 4 ) IUNIT, the FORTRAN unit number associated with the
!    PostScript output file.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) j
  integer ( kind = 4 ), parameter :: ni = 6
  integer ( kind = 4 ), parameter :: nj = 3
  real ( kind = 8 ) u0
  real ( kind = 8 ) u1
  real ( kind = 8 ), parameter :: umax =   30.0D+00
  real ( kind = 8 ), parameter :: umin =  150.0D+00
  real ( kind = 8 ) v0
  real ( kind = 8 ) v1
  real ( kind = 8 ), parameter :: vmax =    5.0D+00
  real ( kind = 8 ), parameter :: vmin =    1.0D+00
!
!  Draw the region:
!
!  Draw the items with horizontal index I:
!
  do i = 1, ni

    u0 = ( real ( ni - i + 1, kind = 8 ) * umin   &
         + real (      i - 1, kind = 8 ) * umax ) &
         / real ( ni,         kind = 8 )

    u1 = ( real ( ni - i, kind = 8 ) * umin   &
         + real (      i, kind = 8 ) * umax ) &
         / real ( ni,     kind = 8     )
!
!  Draw the items with vertical index J:
!
    do j = 1, nj

      v0 = ( real ( nj - j + 1, kind = 8 ) * vmin   &
           + real (      j - 1, kind = 8 ) * vmax ) &
           / real ( nj,         kind = 8 )

      v1 = ( real ( nj - j, kind = 8 ) * vmin   &
           + real (      j, kind = 8 ) * vmax ) &
           / real ( nj,     kind = 8 )
!
!  Now that we have specified a subregion, "draw" the lines associated
!  with the tile.
!
      call sub_box_tiler_2d ( umin, vmin, umax, vmax, u0, v0, u1, v1 )

    end do

  end do

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
!    Input/output, integer ( kind = 4 ) ( kind = 4 ) I, J.  On output, the values of I and
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
!    Henry McGilton, Mary Campione,
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
!    Henry McGilton, Mary Campione,
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
!    Henry McGilton, Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ( kind = 4 ) UNIT, the FORTRAN unit to which output was written.
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
!    Henry McGilton, Mary Campione,
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
!    Henry McGilton, Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ( kind = 4 ) UNIT, the FORTRAN unit to which output should
!    be written.
!
!    Input, character ( len = 80 ) FILE_NAME, the name of the output file.
!
!    Output, integer ( kind = 4 ) ( kind = 4 ) IERROR, error flag.
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
!    Henry McGilton, Mary Campione,
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
!    Henry McGilton, Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ( kind = 4 ) NPOINT, the number of points in the line.
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
subroutine ps_line_open ( npoint, x, y )

!*****************************************************************************80
!
!! PS_LINE_OPEN adds the graph of a line to a PostScript file.
!
!  Discussion:
!
!    The current point is set to the last point in the list.
!
!    This routine does not perform clipping, although it wouldn't be
!    hard to add.
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
!    Henry McGilton, Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ( kind = 4 ) NPOINT, the number of points in the line.
!
!    Input, real ( kind = 8 ) X(NPOINT), Y(NPOINT), the X and Y
!    components of the points.
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
!  Draw the line.
!
  write ( unit, '(a)' ) 'stroke'

  call ps_setting_real ( 'SET', 'XCUR', x(npoint) )
  call ps_setting_real ( 'SET', 'YCUR', y(npoint) )

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
!    Henry McGilton, Mary Campione,
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
!    Henry McGilton, Mary Campione,
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
!! PS_SETTING_INT sets, gets, or prints integer ( kind = 4 ) internal PS_WRITE parameters.
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
!    Henry McGilton, Mary Campione,
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
!    Input/output, integer ( kind = 4 ) ( kind = 4 ) VALUE.
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
!    Henry McGilton, Mary Campione,
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
!! R8_SWAP swaps two real values.
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
subroutine sub_box_tiler_2d ( umin, vmin, umax, vmax, u0, v0, u1, v1 )

!*****************************************************************************80
!
!! SUB_BOX_TILER_2D "tiles" a sub-box with the given pattern.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gordon and Charles Hall,
!    Construction of Curvilinear Coordinate Systems and Application to
!    Mesh Generation,
!    International Journal of Numerical Methods in Engineering,
!    Volume 7, pages 461-477, 1973.
!
!    Joe Thompson, Bharat Soni, Nigel Weatherill,
!    Handbook of Grid Generation,
!    CRC Press,
!    1999.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) UMIN, VMIN, UMAX, VMAX, the minimum and
!    maximum values of U and V.
!
!    Input, real ( kind = 8 ) U0, V0, U1, V1, the (U,V) coordinates of
!    the opposite corners of the sub-box.
!
  implicit none

  integer ( kind = 4 ), parameter :: npoint = 11

  integer ( kind = 4 ) i
  real ( kind = 8 ) r
  real ( kind = 8 ) r_tab(npoint)
  real ( kind = 8 ), parameter :: r0 = 0.0D+00
  real ( kind = 8 ), parameter :: r1 = 1.0D+00
  real ( kind = 8 ), parameter :: s0 = 0.0D+00
  real ( kind = 8 ), parameter :: s1 = 1.0D+00
  real ( kind = 8 ) u
  real ( kind = 8 ) u0
  real ( kind = 8 ) u1
  real ( kind = 8 ) umax
  real ( kind = 8 ) umin
  real ( kind = 8 ) s
  real ( kind = 8 ) s_tab(npoint)
  real ( kind = 8 ) v
  real ( kind = 8 ) v0
  real ( kind = 8 ) v1
  real ( kind = 8 ) vmax
  real ( kind = 8 ) vmin
  real ( kind = 8 ) x(npoint)
  real ( kind = 8 ) x00
  real ( kind = 8 ) x01
  real ( kind = 8 ) x0v
  real ( kind = 8 ) x10
  real ( kind = 8 ) x11
  real ( kind = 8 ) x1v
  real ( kind = 8 ) xu0
  real ( kind = 8 ) xu1
  real ( kind = 8 ) y(npoint)
  real ( kind = 8 ) y00
  real ( kind = 8 ) y01
  real ( kind = 8 ) y0v
  real ( kind = 8 ) y10
  real ( kind = 8 ) y11
  real ( kind = 8 ) y1v
  real ( kind = 8 ) yu0
  real ( kind = 8 ) yu1
!
!  Here is a sequence of points in (R,S) coordinate (between 0 and 1)
!  which describe the "tile" that is to be drawn within each sub-box.
!
  r_tab(1) =  0.3D+00
  s_tab(1) =  0.2D+00
  r_tab(2) =  0.3D+00
  s_tab(2) =  0.8D+00
  r_tab(3) =  0.4D+00
  s_tab(3) =  0.79D+00
  r_tab(4) =  0.5D+00
  s_tab(4) =  0.78D+00
  r_tab(5) =  0.55D+00
  s_tab(5) =  0.75D+00
  r_tab(6) =  0.6D+00
  s_tab(6) =  0.65D+00
  r_tab(7) =  0.55D+00
  s_tab(7) =  0.55D+00
  r_tab(8) =  0.5D+00
  s_tab(8) =  0.42D+00
  r_tab(9) =  0.4D+00
  s_tab(9) =  0.41D+00
  r_tab(10) = 0.3D+00
  s_tab(10) = 0.4D+00
  r_tab(11) = 0.6D+00
  s_tab(11) = 0.2D+00
!
!  Evaluate the (X,Y) coordinates of the corner (U,V) values.
!
  call boundary_2d ( umin, vmin, umax, vmax, umin, vmin, x00, y00 )
  call boundary_2d ( umin, vmin, umax, vmax, umin, vmax, x01, y01 )
  call boundary_2d ( umin, vmin, umax, vmax, umax, vmin, x10, y10 )
  call boundary_2d ( umin, vmin, umax, vmax, umax, vmax, x11, y11 )

  do i = 1, npoint
!
!  Get the point (R,S) in the tile.
!
    r = r_tab(i)
    s = s_tab(i)
!
!  Get the corresponding point (U,V) in the rectangular space.
!
    u = ( ( r1 - r      ) * u0   &
        + (      r - r0 ) * u1 ) &
        / ( r1     - r0 )

    v = ( ( s1 - s      ) * v0   &
        + (      s - s0 ) * v1 ) &
        / ( s1     - s0 )
!
!  Evaluate (X,Y) at (U,0), (U,1), (0,V), (1,V).
!
    call boundary_2d ( umin, vmin, umax, vmax, u, vmin, xu0, yu0 )
    call boundary_2d ( umin, vmin, umax, vmax, u, vmax, xu1, yu1 )
    call boundary_2d ( umin, vmin, umax, vmax, umin, v, x0v, y0v )
    call boundary_2d ( umin, vmin, umax, vmax, umax, v, x1v, y1v )
!
!  Get the corresponding point (X,Y) in the warped space.
!
    x(i) = ( - x00 * ( umax - u    ) * ( vmax - v    ) &
             + x0v * ( umax - u    ) * ( vmax - vmin ) &
             - x01 * ( umax - u    ) * ( v    - vmin ) &
             + xu0 * ( umax - umin ) * ( vmax - v    ) &
             + xu1 * ( umax - umin ) * ( v    - vmin ) &
             - x10 * ( u    - umin ) * ( vmax - v    ) &
             + x1v * ( u    - umin ) * ( vmax - vmin ) &
             - x11 * ( u    - umin ) * ( v    - vmin ) ) &
             / (     ( umax - umin ) * ( vmax - vmin ) )

    y(i) = ( - y00 * ( umax - u    ) * ( vmax - v    ) &
             + y0v * ( umax - u    ) * ( vmax - vmin ) &
             - y01 * ( umax - u    ) * ( v    - vmin ) &
             + yu0 * ( umax - umin ) * ( vmax - v    ) &
             + yu1 * ( umax - umin ) * ( v    - vmin ) &
             - y10 * ( u    - umin ) * ( vmax - v    ) &
             + y1v * ( u    - umin ) * ( vmax - vmin ) &
             - y11 * ( u    - umin ) * ( v    - vmin ) ) &
             / (     ( umax - umin ) * ( vmax - vmin ) )

  end do

  call ps_line_open ( npoint, x, y )

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
