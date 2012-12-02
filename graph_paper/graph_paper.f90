subroutine ch_cap ( ch )

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
!    Input/output, character CH, the character to capitalize.
!
  implicit none

  character ch
  integer ( kind = 4 ) itemp

  itemp = ichar ( ch )

  if ( 97 <= itemp .and. itemp <= 122 ) then
    ch = char ( itemp - 32 )
  end if

  return
end
function degrees_to_radians ( angle )

!*****************************************************************************80
!
!! DEGREES_TO_RADIANS converts an angle from degrees to radians.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ANGLE, an angle in degrees.
!
!    Output, real ( kind = 8 ) DEGREES_TO_RADIANS, the equivalent angle
!    in radians.
!
  implicit none

  real ( kind = 8 ) angle
  real ( kind = 8 ) degrees_to_radians
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00

  degrees_to_radians = ( angle / 180.0D+00 ) * pi

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
!    Input, integer ( kind = 4 ) X_PS_MIN, Y_PS_MIN, X_PS_MAX, Y_PS_MAX, the 
!    minimum and maximum X and Y values of the data, in PostScript units.  Any 
!    data that lies outside this range will not show up properly.  A reasonable
!    set of values might be 0, 0, 612, 792, or, for a half inch margin,
!    36, 36, 576, 756.
!
  implicit none

  character ( len = 8 )  date
  character ( len = * )  file_name
  real ( kind = 8 ) line_blue
  real ( kind = 8 ) line_green
  real ( kind = 8 ) line_red
  integer   ( kind = 4 ) state
  integer   ( kind = 4 ) unit
  integer   ( kind = 4 ) x_ps_max
  integer   ( kind = 4 ) x_ps_min
  integer   ( kind = 4 ) y_ps_max
  integer   ( kind = 4 ) y_ps_min
!
!  Determine if the PostScript state is acceptable.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state /= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EPS_FILE_HEAD - Fatal error!'
    write ( *, '(a,i8)' ) '  PostScript state is ', state
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
    write ( *, '(a,i8)' ) '  PostScript state is ', state
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
    write ( *, '(a,i8,a)' ) '  This file describes ', num_pages, ' pages.'
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
subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is a value between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is a value between 1 and 99, representing a
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
subroutine hexagon_draw ( center, radius, angle, iunit, side )

!*****************************************************************************80
!
!! HEXAGON_DRAW draws a hexagon.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 April 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) CENTER(2), the center of the hexagon, in 
!    PostScript coordinates.
!
!    Input, integer ( kind = 4 ) RADIUS, the radius of the hexagon, in 
!    PostScript coordinates.
!
!    Iniput, integer ( kind = 4 ) ANGLE, the angle through which the hexagon
!    should be turned, in degrees.  If ANGLE is 0, then the first vertex
!    of the hexagon "points" east, and the hexagon is flat on top and bottom.
!
!    Input, integer ( kind = 4 ) IUNIT, the PostScript output unit.
!
!    Input, logical SIDE(6), is true for each side of the hexagon that
!    is to be drawn, starting with the side whose first vertex is ordinarily
!    the easternmost, and proceeding counterclockwise.
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) angle
  integer ( kind = 4 ) center(2)
  integer ( kind = 4 ) h(2,6)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jm1
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ) radius
  logical              side(6)
!
!  Temporarily enforcing limits.
!
  do j = 1, 6

    a = angle + 60 * ( j - 1 )

    h(1,j) = center(1) + int ( radius * cos ( pi * a / 180.0D+00 ) )
!   h(1,j) = min ( max ( h(1,j), 36 ), 576 )
    h(2,j) = center(2) + int ( radius * sin ( pi * a / 180.0D+00 ) )
!   h(2,j) = min ( max ( h(2,j), 36 ), 756 )
  end do

  write ( iunit, '(a)' ) 'newpath'

  jm1 = 6
  do j = 1, 6

    if ( side(j) ) then
      if ( j == 1 ) then
        write ( iunit, '(i4,2x,i4,2x,a)' ) h(1,6), h(2,6), ' moveto'
      else if ( .not. side(j-1) ) then
        write ( iunit, '(i4,2x,i4,2x,a)' ) h(1,j-1), h(2,j-1), ' moveto'
      end if
      write ( iunit, '(i4,2x,i4,2x,a)' ) h(1,j  ), h(2,j  ), ' lineto'
    end if

    jm1 = j

  end do

  write ( iunit, '(a)' ) 'stroke'

  return
end
subroutine hexagonal_1 ( ngrid, file_name )

!*****************************************************************************80
!
!! HEXAGONAL_1 draws a hexagonal grid that roughly forms a triangle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NGRID, the number of "grid lines" drawn in
!    each of the three hexagonal directions.  
!
!    Input, character ( len = * ) FILE_NAME, the name of the file to which the
!    output should be written.
!
  implicit none

  character ( len = * ) file_name
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) ierror
  integer   ( kind = 4 ) iunit
  integer   ( kind = 4 ) j
  integer   ( kind = 4 ) n
  integer   ( kind = 4 ) ngrid
  integer   ( kind = 4 ) pax
  integer   ( kind = 4 ) pay
  integer   ( kind = 4 ) pbx
  integer   ( kind = 4 ) pby
  integer   ( kind = 4 ) pcx
  integer   ( kind = 4 ) pcy
  integer   ( kind = 4 ) px
  integer   ( kind = 4 ) px1
  integer   ( kind = 4 ) px2
  integer   ( kind = 4 ) py
  integer   ( kind = 4 ) py1
  integer   ( kind = 4 ) py2
  integer   ( kind = 4 ) x_ps_max
  integer   ( kind = 4 ) x_ps_min
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  integer   ( kind = 4 ) y_ps_max
  integer   ( kind = 4 ) y_ps_min
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin

  call get_unit ( iunit )
!
!  Open the output file.
!
  call ps_file_open ( file_name, iunit, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HEXAGONAL_1 - Fatal error!'
    write ( *, '(a,i6)' ) '  File creation error ', ierror
    stop
  end if

  x_ps_min = 36
  y_ps_min = 36
  x_ps_max = 576
  y_ps_max = 756

  call eps_file_head ( file_name, x_ps_min, y_ps_min, x_ps_max, &
    y_ps_max )
!
!  Write the header.
!
  xmin = 0.0D+00
  ymin = 0.0D+00

  xmax = 8.5D+00
  ymax = 11.0D+00

  call ps_page_head ( xmin, ymin, xmax, ymax )
!
!  Draw the grid lines.
!
  pax =  36
  pay = 136

  pbx = 576
  pby = 136

  pcx = 306
  pcy = 568

  n = 3 * ngrid

  do i = 0, n-1

    write ( iunit, '(a)' ) ' newpath'

    px1 = int ( real ( ( n - i ) * pax   &
                       +     i   * pbx, kind = 8 ) &
              / real   ( n,             kind = 8 ) )

    py1 = int ( real ( ( n - i ) * pay   &
                           + i   * pby, kind = 8 ) &
              / real   ( n,             kind = 8 ) )

    px2 = int ( real ( ( n - i ) * pcx   &
                           + i   * pbx, kind = 8 ) &
              / real   ( n,             kind = 8 ) )

    py2 = int ( real ( ( n - i ) * pcy   &
                           + i   * pby, kind = 8 ) &
              / real   ( n,             kind = 8 ) )

    do j = mod ( i+1, 3 ), n - i - 1, 3

      px = int ( real ( ( n - i - j ) * px1   &
                                + j   * px2, kind = 8 ) &
               / real   ( n - i,             kind = 8 ) )

      py = int ( real ( ( n - i - j ) * py1   &
                                + j   * py2, kind = 8 ) &
               / real   ( n - i,             kind = 8 ) )

      write ( iunit, '(2i4,a)' ) px, py, ' moveto'

      px = int ( real ( ( n - i - j - 1 ) * px1   &
                      + (         j + 1 ) * px2, kind = 8 ) &
               / real   ( n - i,                 kind = 8 ) )

      py = int ( real ( ( n - i - j - 1 ) * py1   &
                      + (         j + 1 ) * py2, kind = 8 ) &
               / real   ( n - i,                 kind = 8 ) )

      write ( iunit, '(2i4,a)' ) px, py, ' lineto'

    end do

    write ( iunit, '(a)' ) ' stroke'

  end do

  do i = 0, n-1

    write ( iunit, '(a)' ) ' newpath'

    px1 = int ( real ( ( n - i ) * pbx   &
                           + i   * pcx, kind = 8 ) &
              / real   ( n,             kind = 8 ) )

    py1 = int ( real ( ( n - i ) * pby   &
                           + i   * pcy, kind = 8 ) &
              / real   ( n,             kind = 8 ) )

    px2 = int ( real ( ( n - i ) * pax   &
                           + i   * pcx, kind = 8 ) &
              / real   ( n,             kind = 8 ) )

    py2 = int ( real ( ( n - i ) * pay   &
                           + i   * pcy, kind = 8 ) &
              / real   ( n,             kind = 8 ) )

    do j = mod ( i+1, 3 ), n - i - 1, 3

      px = int ( real ( ( n - i - j ) * px1   &
                                + j   * px2, kind = 8  ) &
               / real   ( n - i,             kind = 8 ) )

      py = int ( real ( ( n - i - j ) * py1   &
                                + j   * py2, kind = 8 ) &
               / real   ( n - i,             kind = 8 ) )

      write ( iunit, '(2i4,a)' ) px, py, ' moveto'

      px = int ( real ( ( n - i - j - 1 ) * px1   &
                      + (         j + 1 ) * px2, kind = 8 ) &
               / real   ( n - i,                 kind = 8 ) )

      py = int ( real ( ( n - i - j - 1 ) * py1   &
                      + (         j + 1 ) * py2, kind = 8 ) &
               / real   ( n - i,                 kind = 8 ) )

      write ( iunit, '(2i4,a)' ) px, py, ' lineto'

    end do

    write ( iunit, '(a)' ) ' stroke'

  end do

  do i = 0, n-1

    write ( iunit, '(a)' ) ' newpath'

    px1 = int ( real ( ( n - i ) * pcx   &
                           + i   * pax, kind = 8 ) &
              / real   ( n,             kind = 8 ) )

    py1 = int ( real ( ( n - i ) * pcy   &
                           + i   * pay, kind = 8 ) &
              / real   ( n,             kind = 8 ) )

    px2 = int ( real ( ( n - i ) * pbx   &
                           + i   * pax, kind = 8 ) &
              / real   ( n,             kind = 8 ) )

    py2 = int ( real ( ( n - i ) * pby   &
                           + i   * pay, kind = 8 ) &
              / real   ( n,             kind = 8 ) )

    do j = mod ( i+1, 3 ), n - i - 1, 3

      px = int ( real ( ( n - i - j ) * px1   &
                                + j   * px2, kind = 8 ) &
               / real   ( n - i,             kind = 8 ) )

      py = int ( real ( ( n - i - j ) * py1   &
                                + j   * py2, kind = 8 ) &
               / real   ( n - i,             kind = 8 ) )

      write ( iunit, '(2i4,a)' ) px, py, ' moveto'

      px = int ( real ( ( n - i - j - 1 ) * px1   &
                      + (         j + 1 ) * px2, kind = 8  ) &
               / real   ( n - i,                 kind = 8 ) )

      py = int ( real ( ( n - i - j - 1 ) * py1   &
                      + (         j + 1 ) * py2, kind = 8 ) &
               / real   ( n - i,                 kind = 8 ) )

      write ( iunit, '(2i4,a)' ) px, py, ' lineto'

    end do

    write ( iunit, '(a)' ) ' stroke'

  end do

  call ps_page_tail ( )

  call eps_file_tail ( )

  call ps_file_close ( iunit )

  return
end
subroutine hexagonal_2 ( n, file_name )

!*****************************************************************************80
!
!! HEXAGONAL_2 draws a hexagonal grid of dots (no lines).
!
!  Discussion:
!
!    The hexagonal form is emphasized by omitting what amounts to
!    the central dot of each hexagon.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 January 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, specifies the number of grid lines drawn.
!
!    Input, character ( len = * ) FILE_NAME, the name of the file to which the
!    output should be written.
!
  implicit none

  character ( len = * )  file_name
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) ierror
  integer   ( kind = 4 ) iunit
  integer   ( kind = 4 ) j
  integer   ( kind = 4 ) marker_size
  integer   ( kind = 4 ) n
  integer   ( kind = 4 ) pax
  integer   ( kind = 4 ) pay
  integer   ( kind = 4 ) pbx
  integer   ( kind = 4 ) pby
  integer   ( kind = 4 ) pcx
  integer   ( kind = 4 ) pcy
  integer   ( kind = 4 ) px
  integer   ( kind = 4 ) py
  integer   ( kind = 4 ) x_ps_max
  integer   ( kind = 4 ) x_ps_min
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  integer   ( kind = 4 ) y_ps_max
  integer   ( kind = 4 ) y_ps_min
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin
!
!  Open the output file.
!
  call get_unit ( iunit )

  call ps_file_open ( file_name, iunit, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HEXAGONAL_2 - Fatal error!'
    write ( *, '(a,i6)' ) '  File creation error ', ierror
    stop
  end if

  x_ps_min = 36
  y_ps_min = 36
  x_ps_max = 576
  y_ps_max = 756

  call eps_file_head ( file_name, x_ps_min, y_ps_min, x_ps_max, &
    y_ps_max )
!
!  Write the header.
!
  xmin = 0.0D+00
  ymin = 0.0D+00

  xmax = 8.5D+00
  ymax = 11.0D+00

  call ps_page_head ( xmin, ymin, xmax, ymax )
!
!  Draw the points.
!
  pax =  36
  pay = 136

  pbx = 576
  pby = 568
  
  marker_size = 4

  do j = 0, n

    py = int ( real ( ( n - j ) * pay   &
                          + j   * pby, kind = 8 ) &
               / real ( n,             kind = 8 ) )

    if ( mod ( j, 2 ) == 0 ) then

      do i = 1, n

        if ( mod ( i, 3 ) /= 2 ) then

          px = int ( real ( ( 2 * n - 2 * i + 1 ) * pax   &
                          + (       + 2 * i - 1 ) * pbx, kind = 8 ) &
                   / real   ( 2 * n,                     kind = 8 ) )

          write ( iunit, '(a,3i6,a)' ) 'newpath ', px, py, marker_size, &
            ' 0 360 arc closepath fill'

        end if

      end do

    else

      do i = 0, n

        if ( mod ( i, 3 ) /= 0 ) then

          px = int ( real ( ( n - i ) * pax   &
                                + i   * pbx, kind = 8 ) &
                   / real   ( n,             kind = 8 ) )

          write ( iunit, '(a,3i6,a)' ) 'newpath ', px, py, marker_size, &
            ' 0 360 arc closepath fill'

        end if

      end do

    end if

  end do

  call ps_page_tail ( )

  call eps_file_tail ( )

  call ps_file_close ( iunit )

  return
end
subroutine hexagonal_3 ( n, file_name )

!*****************************************************************************80
!
!! HEXAGONAL_3 draws a hexagonal grid that fills a rectangle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 April 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, specifies the number of grid lines drawn.
!
!    Input, character ( len = * ) FILE_NAME, the name of the file to which the
!    output should be written.
!
  implicit none

  character ( len = * ) file_name
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) ierror
  integer   ( kind = 4 ) iunit
  integer   ( kind = 4 ) j
  integer   ( kind = 4 ) marker_size
  integer   ( kind = 4 ) n
  integer   ( kind = 4 ) pax
  integer   ( kind = 4 ) pay
  integer   ( kind = 4 ) pbx
  integer   ( kind = 4 ) pby
  integer   ( kind = 4 ) pcx
  integer   ( kind = 4 ) pcy
  integer   ( kind = 4 ) px
  integer   ( kind = 4 ) py
  integer   ( kind = 4 ) x_ps_max
  integer   ( kind = 4 ) x_ps_min
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  integer   ( kind = 4 ) y_ps_max
  integer   ( kind = 4 ) y_ps_min
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin
!
!  Open the output file.
!
  call get_unit ( iunit )

  call ps_file_open ( file_name, iunit, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HEXAGONAL_3 - Fatal error!'
    write ( *, '(a,i6)' ) '  File creation error ', ierror
    stop
  end if

  x_ps_min = 36
  y_ps_min = 36
  x_ps_max = 576
  y_ps_max = 756

  call eps_file_head ( file_name, x_ps_min, y_ps_min, x_ps_max, &
    y_ps_max )
!
!  Write the header.
!
  xmin = 0.0D+00
  ymin = 0.0D+00

  xmax = 8.5D+00
  ymax = 11.0D+00

  call ps_page_head ( xmin, ymin, xmax, ymax )
!
!  Draw the points.
!
  pax =  36
  pay = 136

  pbx = 576
  pby = 568
!
!  Horizontal lines.
!
  do j = 0, n+1

    py = int ( real ( ( n - j ) * pay   &
                          + j   * pby, kind = 8 ) &
             / real   ( n,             kind = 8 ) )

    if ( mod ( j, 2 ) == 0 ) then

      do i = 4, n, 3

        px = int ( real ( ( 2 * n - 2 * i + 1 ) * pax   &
                        + (       + 2 * i - 1 ) * pbx, kind = 8 ) &
                 / real   ( 2 * n,                     kind = 8 ) )

        write ( iunit, '(2i4,a)' ) px, py, ' moveto'

        px = int ( real ( ( 2 * n - 2 * i + 3 ) * pax   &
                        + (       + 2 * i - 3 ) * pbx, kind = 8 ) &
                 / real   ( 2 * n,                     kind = 8 ) )

        write ( iunit, '(2i4,a)' ) px, py, ' lineto stroke'

      end do

    else

      do i = 1, n, 3

        if ( mod ( i, 3 ) /= 0 ) then

          px = int ( real ( ( n - i ) * pax   &
                                + i   * pbx, kind = 8 ) &
                   / real   ( n,             kind = 8 ) )

          write ( iunit, '(2i4,a)' ) px, py, ' moveto'

          px = int ( real ( ( n - i - 1 ) * pax   &
                          + (     i + 1 ) * pbx, kind = 8 ) &
                   / real   ( n,                 kind = 8 ) )

          write ( iunit, '(2i4,a)' ) px, py, ' lineto stroke'


        end if

      end do

    end if

  end do
!
!  Slanted lines.
!
  do j = 0, n

    if ( mod ( j, 2 ) == 0 ) then

      do i = 1, n, 3

        px = int ( real ( ( 2 * n - 2 * i + 1 ) * pax   &
                        + (       + 2 * i - 1 ) * pbx, kind = 8 ) &
                 / real   ( 2 * n,                     kind = 8 ) )

        py = int ( real ( ( n - j ) * pay   &
                              + j   * pby, kind = 8 ) &
                 / real   ( n,             kind = 8 ) )
  
        write ( iunit, '(2i4,a)' ) px, py, ' moveto'

        px = int ( real ( ( n - i ) * pax   &
                              + i   * pbx, kind = 8 ) &
                 / real   ( n,             kind = 8 ) )

        py = int ( real ( ( n - j - 1 ) * pay   &
                          +   ( j + 1 ) * pby, kind = 8 ) &
                 / real   ( n,                 kind = 8 ) )

        write ( iunit, '(2i4,a)' ) px, py, ' lineto stroke'

        px = int ( real ( ( 2 * n - 2 * i - 3 ) * pax   &
                        + (       + 2 * i + 3 ) * pbx, kind = 8 ) &
                 / real   ( 2 * n,                     kind = 8 ) )

        py = int ( real ( ( n - j ) * pay   &
                              + j   * pby, kind = 8 ) &
                 / real   ( n,             kind = 8 ) )

        write ( iunit, '(2i4,a)' ) px, py, ' moveto'

        px = int ( real ( ( n - i - 1 ) * pax   &
                        + (     i + 1 ) * pbx, kind = 8 ) &
                 / real   ( n,                 kind = 8 ) )

        py = int ( real ( ( n - j - 1 ) * pay   &
                          +   ( j + 1 ) * pby, kind = 8 ) &
                 / real   ( n,                 kind = 8 ) )

        write ( iunit, '(2i4,a)' ) px, py, ' lineto stroke'

      end do

    else

      do i = 1, n, 3

        px = int ( real ( ( 2 * n - 2 * i + 1 ) * pax   &
                        + (       + 2 * i - 1 ) * pbx, kind = 8 ) &
                 / real   ( 2 * n,                     kind = 8 ) )

        py = int ( real ( ( n - j - 1 ) * pay   &
                          +   ( j + 1 ) * pby, kind = 8 ) &
                 / real   ( n,                 kind = 8 ) )

        write ( iunit, '(2i4,a)' ) px, py, ' moveto'

        px = int ( real ( ( n - i ) * pax   &
                              + i   * pbx, kind = 8 ) &
                 / real   ( n,             kind = 8 ) )

        py = int ( real ( ( n - j ) * pay   &
                              + j   * pby, kind = 8 ) &
                 / real   ( n,             kind = 8 ) )

        write ( iunit, '(2i4,a)' ) px, py, ' lineto stroke'

        px = int ( real ( ( 2 * n - 2 * i - 3 ) * pax   &
                        + (       + 2 * i + 3 ) * pbx, kind = 8 ) &
                 / real   ( 2 * n,                     kind = 8 ) )

        py = int ( real ( ( n - j - 1 ) * pay   &
                          +   ( j + 1 ) * pby, kind = 8 ) &
                 / real   ( n,                 kind = 8 ) )

        write ( iunit, '(2i4,a)' ) px, py, ' moveto'

        px = int ( real ( ( n - i - 1 ) * pax   &
                          +   ( i + 1 ) * pbx, kind = 8 ) &
                 / real   ( n,                 kind = 8 ) )
 
        py = int ( real ( ( n - j ) * pay   &
                              + j   * pby, kind = 8 ) &
                 / real   ( n,             kind = 8 ) )

        write ( iunit, '(2i4,a)' ) px, py, ' lineto stroke'

      end do
    end if

  end do

  call ps_page_tail ( )

  call eps_file_tail ( )

  call ps_file_close ( iunit )

  return
end
subroutine hexagonal_4 ( n, file_name )

!*****************************************************************************80
!
!! HEXAGONAL_4 draws a hexagonal grid of dots (no lines).
!
!  Discussion:
!
!    This routine is similar to HEXAGONAL_2, but does not omit the
!    central dot of each hexagon.  Thus, it should look like a dense
!    array of dots on a hexagonal grid.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, specifies the number of grid lines drawn.
!
!    Input, character ( len = * ) FILE_NAME, the name of the file to which the
!    output should be written.
!
  implicit none

  character ( len = * )  file_name
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) ierror
  integer   ( kind = 4 ) iunit
  integer   ( kind = 4 ) j
  integer   ( kind = 4 ) marker_size
  integer   ( kind = 4 ) n
  integer   ( kind = 4 ) pax
  integer   ( kind = 4 ) pay
  integer   ( kind = 4 ) pbx
  integer   ( kind = 4 ) pby
  integer   ( kind = 4 ) pcx
  integer   ( kind = 4 ) pcy
  integer   ( kind = 4 ) px
  integer   ( kind = 4 ) py
  integer   ( kind = 4 ) x_ps_max
  integer   ( kind = 4 ) x_ps_min
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  integer   ( kind = 4 ) y_ps_max
  integer   ( kind = 4 ) y_ps_min
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin
!
!  Open the output file.
!
  call get_unit ( iunit )

  call ps_file_open ( file_name, iunit, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HEXAGONAL_4 - Fatal error!'
    write ( *, '(a,i6)' ) '  File creation error ', ierror
    stop
  end if

  x_ps_min = 36
  y_ps_min = 36
  x_ps_max = 576
  y_ps_max = 756

  call eps_file_head ( file_name, x_ps_min, y_ps_min, x_ps_max, &
    y_ps_max )
!
!  Write the header.
!
  xmin = 0.0D+00
  ymin = 0.0D+00

  xmax = 8.5D+00
  ymax = 11.0D+00

  call ps_page_head ( xmin, ymin, xmax, ymax )
!
!  Draw the points.
!
  pax =  36
  pay = 136

  pbx = 576
  pby = 568
  
  marker_size = 4

  do j = 0, n

    py = int ( real ( ( n - j ) * pay   &
                          + j   * pby, kind = 8 ) &
               / real ( n,             kind = 8 ) )

    if ( mod ( j, 2 ) == 0 ) then

      do i = 1, n

!       if ( mod ( i, 3 ) /= 2 ) then

          px = int ( real ( ( 2 * n - 2 * i + 1 ) * pax   &
                          + (       + 2 * i - 1 ) * pbx, kind = 8 ) &
                   / real   ( 2 * n,                     kind = 8 ) )

          write ( iunit, '(a,3i6,a)' ) 'newpath ', px, py, marker_size, &
            ' 0 360 arc closepath fill'

!       end if

      end do

    else

      do i = 0, n

!       if ( mod ( i, 3 ) /= 0 ) then

          px = int ( real ( ( n - i ) * pax   &
                                + i   * pbx, kind = 8 ) &
                   / real   ( n,             kind = 8 ) )

          write ( iunit, '(a,3i6,a)' ) 'newpath ', px, py, marker_size, &
            ' 0 360 arc closepath fill'

!       end if

      end do

    end if

  end do

  call ps_page_tail ( )

  call eps_file_tail ( )

  call ps_file_close ( iunit )

  return
end
subroutine hexagonal_5 ( n, file_name )

!*****************************************************************************80
!
!! HEXAGONAL_5 draws a Hex board.
!
!  Discussion:
!
!    Hex is a game that John Nash (and/or Piet Hein) invented.  It could be 
!    played on the hexagonal tiles in a bathroom.  Nash recommended using a 
!    region that was described as a 14x14 hexagonal grid.  You can think of 
!    constructing the board by starting with a strip of 14 hexagons that form
!    a line, and joining 14 such strips in a way that forms an almost 
!    diamond-shaped region.
!
!    There are four sides of the diamond.  The Red player takes two opposing
!    sides, and the Blue player takes the other two.  The object is to construct
!    a bridge (a connected sequence of hexagons belonging to you) joining
!    your two sides.  A turn consists of placing a dot of your color in an 
!    unclaimed hexagon, which claims the hexagon for your side.
!
!    The current version of this routine produces an unlovely image.  Each
!    hexagon is drawn entirely, which results in most lines being drawn twice,
!    but with slightly different positions.  Moreover, the computation of
!    the sizes of the hexagons is imperfect.  All this to be addressed
!    sometime.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 April 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, specifies the number of grid lines drawn.
!
!    Input, character ( len = * ) FILE_NAME, the name of the file to which the
!    output should be written.
!
  implicit none

  integer ( kind = 4 ) angle
  integer ( kind = 4 ) center(2)
  character ( len = * ) file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) marker_size
  integer ( kind = 4 ) n
  integer ( kind = 4 ) pax
  integer ( kind = 4 ) pay
  integer ( kind = 4 ) pbx
  integer ( kind = 4 ) pby
  integer ( kind = 4 ) pcx
  integer ( kind = 4 ) pcy
  integer ( kind = 4 ) px
  integer ( kind = 4 ) py
  integer ( kind = 4 ) radius
  logical              side(6)
  integer ( kind = 4 ) x_ps_max
  integer ( kind = 4 ) x_ps_min
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  integer ( kind = 4 ) y_ps_max
  integer ( kind = 4 ) y_ps_min
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin
!
!  Open the output file.
!
  call get_unit ( iunit )

  call ps_file_open ( file_name, iunit, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HEXAGONAL_5 - Fatal error!'
    write ( *, '(a,i8)' ) '  File creation error ', ierror
    stop
  end if
!
!  Use the whole 8.5x11 page, except for a 36 point margin.
!
  x_ps_min = 36
  y_ps_min = 36
  x_ps_max = 576
  y_ps_max = 756

  call eps_file_head ( file_name, x_ps_min, y_ps_min, x_ps_max, &
    y_ps_max )
!
!  Write the header.
!
  xmin = 0.0D+00
  ymin = 0.0D+00

  xmax = 8.5D+00
  ymax = 11.0D+00

  call ps_page_head ( xmin, ymin, xmax, ymax )
!
!  Draw the points.
!
!  Y: 0 to 792
!  X: 0 to 612
!
!  I counts the row in which the center of a hexagon is found.
!  There will be (N-1)+1+(N-1) = 2*N-1 such rows.
!
  radius = ( 576 - 36 ) / n / 2 / 1.585
  angle = 30
  side(1:6) = .true.

  jlo = n
  jhi = n

  do i = 1, 2*n-1

    center(2) = 360 + ( i - n ) * 1.5 * radius
!
!  J counts the column in which the center of a hexagon is found.
!  There are N+(N-1) = 2*N-1 such columns, staggered.
!
    do j = jlo, jhi

      if ( mod ( i, 2 ) == 1 ) then
        center(1) = 36 + ( 2 * j - 1 ) * ( 576 - 36 ) / 2 / ( 2 * n - 1 )
      else
        center(1) = 36 + (     j - 1 ) * ( 576 - 36 )     / ( 2 * n - 1 )
      end if
!
!  Draw the hexagon at location (I,J).
!  ALWAYS draw 
!    left top, 
!    right top, 
!    right side.
!  If J = JLO and I <= N, draw
!    left bottom
!  If J = JLO and N <= I draw
!    left side
!  If J = JHI and I <= N, draw
!    right bottom
!
      call hexagon_draw ( center, radius, angle, iunit, side )

    end do

    if ( i < n ) then
      if ( mod ( i, 2 ) == 0 ) then
        jlo = jlo - 1
      else
        jhi = jhi + 1
      end if
    else
      if ( mod ( i, 2 ) == 1 ) then
        jlo = jlo + 1
      else
        jhi = jhi - 1
      end if
    end if

  end do

  call ps_page_tail ( )

  call eps_file_tail ( )

  call ps_file_close ( iunit )

  return
end
subroutine i4_swap ( i, j )

!*****************************************************************************80
!
!! I4_SWAP swaps two I4 values.
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
subroutine polar_1 ( angle_num, circle_num, file_name )

!*****************************************************************************80
!
!! POLAR_1 draws a polar grid.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ANGLE_NUM, the number of angular grid lines.
!
!    Input, integer ( kind = 4 ) CIRCLE_NUM, the number of circular grid lines.
!
!    Input, character ( len = * ) FILE_NAME, the name of the file to which the
!    output should be written.
!
  implicit none

  real ( kind = 8 ) angle
  integer ( kind = 4 ) angle_max
  integer ( kind = 4 ) angle_min
  integer ( kind = 4 ) angle_num
  integer ( kind = 4 ) circle
  integer ( kind = 4 ) circle_num
  real ( kind = 8 ) degrees_to_radians
  character ( len = * ) file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) pr
  integer ( kind = 4 ) prmax
  integer ( kind = 4 ) px
  integer ( kind = 4 ) pxcen
  integer ( kind = 4 ) pxmax
  integer ( kind = 4 ) pxmin
  integer ( kind = 4 ) py
  integer ( kind = 4 ) pycen
  integer ( kind = 4 ) pymax
  integer ( kind = 4 ) pymin
  integer ( kind = 4 ) thick
  integer ( kind = 4 ) thin
  integer ( kind = 4 ) x_ps_max
  integer ( kind = 4 ) x_ps_min
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  integer ( kind = 4 ) y_ps_max
  integer ( kind = 4 ) y_ps_min
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin
!
!  Open the output file.
!
  call get_unit ( iunit )

  call ps_file_open ( file_name, iunit, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'POLAR_1 - Fatal error!'
    write ( *, '(a,i6)' ) '  File creation error ', ierror
    stop
  end if

  x_ps_min = 36
  y_ps_min = 36
  x_ps_max = 576
  y_ps_max = 756

  call eps_file_head ( file_name, x_ps_min, y_ps_min, x_ps_max, &
    y_ps_max )
!
!  Write the header.
!
  xmin = 0.0D+00
  ymin = 0.0D+00

  xmax = 8.5D+00
  ymax = 11.0D+00

  call ps_page_head ( xmin, ymin, xmax, ymax )
!
!  Draw the radial grid lines.
!
  pxmin = 0 + 36
  pxmax = 612 - 36

  pymin = 0 + 36
  pymax = 792 - 36

  pxcen = 0.5D+00 * ( pxmin + pxmax )
  pycen = 0.5D+00 * ( pymin + pymax )

  prmax = pxmax - pxcen

  thin = 1
  thick = 2

  do i = 0, angle_num

    angle = real ( ( i - 1 ) * 360, kind = 8 ) / real ( angle_num, kind = 8 )

    write ( iunit, '(a)' ) ' newpath'
    px = pxcen
    py = pycen
    write ( iunit, '(2i4,a)' ) px, py, ' moveto'

    px = pxcen + real ( prmax, kind = 8 ) &
      * cos ( degrees_to_radians ( angle ) )

    py = pycen + real ( prmax, kind = 8 ) &
      * sin ( degrees_to_radians ( angle ) )

    write ( iunit, '(2i4,a)' ) px, py, ' lineto'

    write ( iunit, '(a)' ) ' stroke'

  end do
!
!  Draw the circular grid lines.
!
  angle_min = 0
  angle_max = 360

  do circle = 1, circle_num

    if ( mod ( circle, 2 ) == 0 ) then
      write ( iunit, '(i4,a)' ) thick, ' setlinewidth'
    else 
      write ( iunit, '(i4,a)' ) thin, ' setlinewidth'
    end if

    pr = int ( real ( circle * prmax, kind = 8 ) &
      / real ( circle_num, kind = 8 ) )

    write ( iunit, '(a)' ) ' newpath'
    write ( iunit, '(5i4,a)' ) pxcen, pycen, pr, angle_min, angle_max, ' arc'
    write ( iunit, '(a)' ) ' closepath stroke'

  end do

  call ps_page_tail ( )

  call eps_file_tail ( )

  call ps_file_close ( iunit )

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
    write ( *, '(a,i8)' ) '  PostScript state is ', state
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
!    Input, integer ( kind = 4 ) UNIT, the FORTRAN unit to which output 
!    was written.
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
    write ( *, '(a,i8)' ) '  PostScript state is ', state
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

  character ( len = * )  file_name
  integer   ( kind = 4 ) ierror
  integer   ( kind = 4 ) ios
  integer   ( kind = 4 ) state
  integer   ( kind = 4 ) unit
!
!  Determine if the PostScript state is acceptable.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_FILE_OPEN - Fatal error!'
    write ( *, '(a,i8)' ) '  PostScript state is ', state
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
    write ( *, '(a,i8)' ) '  PostScript state is ', state
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
    write ( *, '(a,i8)' ) '  PostScript state is ', state
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
  write ( *, * ) 'PS_PAGE_HEAD: Call PS_COLOR_LINE ( PUSH, R, G, B)'

  call ps_color_line ( 'PUSH', line_red, line_green, line_blue )

  call ps_comment ( 'Draw a gray border around the page.' )

  xvec(1:4) = (/ xmin, xmax, xmax, xmin /)
  yvec(1:4) = (/ ymin, ymin, ymax, ymax /)

  call ps_line_closed ( 4, xvec, yvec )

  write ( *, * ) 'PS_PAGE_HEAD: Call PS_COLOR_LINE ( POP, R, G, B)'
  call ps_color_line ( 'POP', line_red, line_green, line_blue )
  write ( *, * ) 'PS_PAGE_HEAD: ALL DONE'
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
    write ( *, '(a,i8)' ) '  PostScript state is ', state
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
      write ( *, '(a,i8)' ) 'Line width, LINE_WIDTH = ', line_width
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
      write ( *, '(a,i8)' ) 'Marker size, MARKER_SIZE = ', marker_size
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
      write ( *, '(a,i8)' ) 'Number of pages, NUM_PAGES = ', num_pages
    else if ( action == 'SET' ) then
      num_pages = value
    end if

  else if ( variable == 'PXMIN' ) then

    if ( action == 'GET' ) then
      value = pxmin
    else if ( action == 'PRINT' ) then
      write ( *, '(a,i8)' ) 'PostScript minimum X point, PXMIN = ', pxmin
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
      write ( *, '(a,i8)' ) 'PostScript minimum Y point, PYMIN = ', pymin
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
      write ( *, '(a,i8)' ) 'Current internal state, STATE = ', state
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
      write ( *, '(a,i8)' ) 'Current FORTRAN unit, UNIT = ', unit
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
subroutine staggered_2 ( n, file_name )

!*****************************************************************************80
!
!! STAGGERED_2 draws a staggered grid of dots (no lines).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 January 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, specifies the number of grid lines drawn.
!
!    Input, character ( len = * ) FILE_NAME, the name of the file to which the
!    output should be written.
!
  implicit none

  character ( len = * )  file_name
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) ierror
  integer   ( kind = 4 ) iunit
  integer   ( kind = 4 ) j
  integer   ( kind = 4 ) marker_size
  integer   ( kind = 4 ) n
  integer   ( kind = 4 ) pax
  integer   ( kind = 4 ) pay
  integer   ( kind = 4 ) pbx
  integer   ( kind = 4 ) pby
  integer   ( kind = 4 ) pcx
  integer   ( kind = 4 ) pcy
  integer   ( kind = 4 ) px
  integer   ( kind = 4 ) py
  integer   ( kind = 4 ) x_ps_max
  integer   ( kind = 4 ) x_ps_min
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  integer   ( kind = 4 ) y_ps_max
  integer   ( kind = 4 ) y_ps_min
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin
!
!  Open the output file.
!
  call get_unit ( iunit )

  call ps_file_open ( file_name, iunit, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'STAGGERED_2 - Fatal error!'
    write ( *, '(a,i6)' ) '  File creation error ', ierror
    stop
  end if

  x_ps_min = 36
  y_ps_min = 36
  x_ps_max = 576
  y_ps_max = 756

  call eps_file_head ( file_name, x_ps_min, y_ps_min, x_ps_max, &
    y_ps_max )
!
!  Write the header.
!
  xmin = 0.0D+00
  ymin = 0.0D+00

  xmax = 8.5D+00
  ymax = 11.0D+00

  call ps_page_head ( xmin, ymin, xmax, ymax )
!
!  Draw the points.
!
  pax =  36
  pay = 136

  pbx = 576
  pby = 568
  
  marker_size = 4

  do j = 0, n

    py = int ( real ( ( n - j ) * pay   &
                          + j   * pby, kind = 8 ) &
             / real   ( n,             kind = 8 ) )

    if ( mod ( j, 2 ) == 0 ) then

      do i = 0, n

        px = int ( real ( ( n - i ) * pax   &
                              + i   * pbx, kind = 8 ) &
                / real ( n,             kind = 8 ) )
 
        write ( iunit, '(a,3i6,a)' ) 'newpath ', px, py, marker_size, &
          ' 0 360 arc closepath fill'

      end do

    else

      do i = 1, n

        px = int ( real ( ( 2 * n - 2 * i + 1 ) * pax   &
                        + (       + 2 * i - 1 ) * pbx, kind = 8 ) &
                 / real   ( 2 * n,                     kind = 8 ) )

        write ( iunit, '(a,3i6,a)' ) 'newpath ', px, py, marker_size, &
          ' 0 360 arc closepath fill'

      end do

    end if

  end do

  call ps_page_tail ( )

  call eps_file_tail ( )

  call ps_file_close ( iunit )

  return
end
subroutine sudoku_sheet_blank ( file_name )

!*****************************************************************************80
!
!! SUDOKU_SHEET_BLANK draws a sheet of 12 blank sudoku grids.
!
!  Discussion:
!
!    12 sudoku grids are plotted in a pattern of 4 rows and 3 columns.
!
!    Each grid is a 9 x 9 box, divided by heavy lines into 
!    3 x 3 subboxes.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 June 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the file to which the
!    output should be written.
!
  implicit none

  character ( len = * ) file_name
  integer ( kind = 4 ) grid_column
  integer ( kind = 4 ) grid_px_max
  integer ( kind = 4 ) grid_px_min
  integer ( kind = 4 ) grid_px_num
  integer ( kind = 4 ) grid_px_space
  integer ( kind = 4 ) grid_px_width
  integer ( kind = 4 ) grid_py_max
  integer ( kind = 4 ) grid_py_min
  integer ( kind = 4 ) grid_py_num
  integer ( kind = 4 ) grid_py_space
  integer ( kind = 4 ) grid_py_width
  integer ( kind = 4 ) grid_row
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  integer ( kind = 4 ) px
  integer ( kind = 4 ) pxinc
  integer ( kind = 4 ) pxmax
  integer ( kind = 4 ) pxmin
  integer ( kind = 4 ) py
  integer ( kind = 4 ) pyinc
  integer ( kind = 4 ) pymax
  integer ( kind = 4 ) pymin
  integer ( kind = 4 ) x_ps_max
  integer ( kind = 4 ) x_ps_min
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  integer ( kind = 4 ) y_ps_max
  integer ( kind = 4 ) y_ps_min
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin
!
!  Open the output file.
!
  call get_unit ( iunit )

  call ps_file_open ( file_name, iunit, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SUDOKU_SHEET_BLANK - Fatal error!'
    write ( *, '(a,i8)' ) '  File creation error ', ierror
    stop
  end if

  x_ps_min = 36
  y_ps_min = 36
  x_ps_max = 576
  y_ps_max = 756

  call eps_file_head ( file_name, x_ps_min, y_ps_min, x_ps_max, &
    y_ps_max )
!
!  Write the header.
!
  xmin = 0.0D+00
  ymin = 0.0D+00

  xmax = 8.5D+00
  ymax = 11.0D+00

  call ps_page_head ( xmin, ymin, xmax, ymax )
!
!  Draw the grid lines.
!
  pxmin = 0 + 36
  pxmax = 612 - 36

  pymin = 0 + 36
  pymax = 792 - 36

  grid_px_space = 18
  grid_py_space = 18

  grid_px_num = 3
  grid_py_num = 4

  grid_px_width = ( pxmax - pxmin - ( grid_px_num - 1 ) * grid_px_space ) &
    / grid_px_num

  grid_py_width = ( pymax - pymin - ( grid_py_num - 1 ) * grid_py_space ) &
    / grid_py_num

  grid_py_min = pymax + grid_py_space

  do grid_row = 1, grid_py_num

    grid_py_max = grid_py_min - grid_py_space
    grid_py_min = grid_py_max - grid_py_width

    grid_px_max = pxmin - grid_px_space
 
    do grid_column = 1, grid_px_num

      grid_px_min = grid_px_max + grid_px_space
      grid_px_max = grid_px_min + grid_px_width

      call sudoku_grid_blank ( iunit, grid_px_min, grid_py_min, grid_px_max, &
        grid_py_max )

    end do

  end do

  call ps_page_tail ( )

  call eps_file_tail ( )

  call ps_file_close ( iunit )

  return
end
subroutine sudoku_sheet_filled ( file_name )

!*****************************************************************************80
!
!! SUDOKU_SHEET_FILLED draws a sheet of 12 filled sudoku grids.
!
!  Discussion:
!
!    This routine is still being experimented upon.
!
!    12 sudoku grids are plotted in a pattern of 4 rows and 3 columns.
!
!    Each grid is a 9 x 9 box, divided by heavy lines into 
!    3 x 3 subboxes.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the file to which the
!    output should be written.
!
  implicit none

  character ( len = * ) file_name
  integer ( kind = 4 ) grid_column
  integer ( kind = 4 ) grid_px_max
  integer ( kind = 4 ) grid_px_min
  integer ( kind = 4 ) grid_px_num
  integer ( kind = 4 ) grid_px_space
  integer ( kind = 4 ) grid_px_width
  integer ( kind = 4 ) grid_py_max
  integer ( kind = 4 ) grid_py_min
  integer ( kind = 4 ) grid_py_num
  integer ( kind = 4 ) grid_py_space
  integer ( kind = 4 ) grid_py_width
  integer ( kind = 4 ) grid_row
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  integer ( kind = 4 ) px
  integer ( kind = 4 ) pxinc
  integer ( kind = 4 ) pxmax
  integer ( kind = 4 ) pxmin
  integer ( kind = 4 ) py
  integer ( kind = 4 ) pyinc
  integer ( kind = 4 ) pymax
  integer ( kind = 4 ) pymin
  integer ( kind = 4 ), dimension ( 9, 9 ) :: value = reshape ( &
    (/ &
    0, 6, 0, 0, 7, 9, 4, 0, 0, &
    0, 0, 0, 0, 0, 4, 0, 3, 6, &
    3, 2, 4, 0, 0, 1, 0, 0, 0, &
    0, 0, 0, 9, 5, 6, 0, 0, 4, &
    6, 0, 8, 0, 0, 0, 2, 0, 7, &
    9, 4, 0, 0, 0, 0, 1, 0, 5, &
    5, 0, 2, 8, 0, 0, 0, 4, 0, &
    0, 8, 7, 0, 3, 0, 0, 5, 0, &
    0, 0, 0, 4, 9, 0, 7, 2, 0 /), (/ 9, 9 /) )
  integer ( kind = 4 ) x_ps_max
  integer ( kind = 4 ) x_ps_min
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  integer ( kind = 4 ) y_ps_max
  integer ( kind = 4 ) y_ps_min
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin
!
  value = transpose ( value )
!
!  Open the output file.
!
  call get_unit ( iunit )

  call ps_file_open ( file_name, iunit, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SUDOKU_SHEET_FILLED - Fatal error!'
    write ( *, '(a,i8)' ) '  File creation error ', ierror
    stop
  end if

  x_ps_min = 36
  y_ps_min = 36
  x_ps_max = 576
  y_ps_max = 756

  call eps_file_head ( file_name, x_ps_min, y_ps_min, x_ps_max, &
    y_ps_max )
!
!  Write the header.
!
  xmin = 0.0D+00
  ymin = 0.0D+00

  xmax = 8.5D+00
  ymax = 11.0D+00

  call ps_page_head ( xmin, ymin, xmax, ymax )
!
!  Draw the grid lines.
!
  pxmin = 0 + 36
  pxmax = 612 - 36

  pymin = 0 + 36
  pymax = 792 - 36

  grid_px_space = 18
  grid_py_space = 18

  grid_px_num = 3
  grid_py_num = 4

  grid_px_width = ( pxmax - pxmin - ( grid_px_num - 1 ) * grid_px_space ) &
    / grid_px_num

  grid_py_width = ( pymax - pymin - ( grid_py_num - 1 ) * grid_py_space ) &
    / grid_py_num

  grid_py_min = pymax + grid_py_space

!
!  Set the font.
!
  write ( iunit, '(a)' ) '/Times-Roman findfont'
  write ( iunit, '(a)' ) '0.25 inch scalefont'
  write ( iunit, '(a)' ) 'setfont'

  do grid_row = 1, 1
! do grid_row = 1, grid_py_num

    grid_py_max = grid_py_min - grid_py_space
    grid_py_min = grid_py_max - grid_py_width

    grid_px_max = pxmin - grid_px_space
 
    do grid_column = 1, 1
!   do grid_column = 1, grid_px_num

      grid_px_min = grid_px_max + grid_px_space
      grid_px_max = grid_px_min + grid_px_width

      call sudoku_grid_filled ( iunit, grid_px_min, grid_py_min, grid_px_max, &
        grid_py_max, value )

    end do

  end do

  call ps_page_tail ( )

  call eps_file_tail ( )

  call ps_file_close ( iunit )

  return
end
subroutine sudoku_grid_blank ( iunit, grid_px_min, grid_py_min, grid_px_max, &
  grid_py_max )

!*****************************************************************************80
!
!! SUDOKU_GRID_BLANK draws one blank sudoku grid.
!
!  Discussion:
!
!    This routine should be called by SUDOKU_SHEET_BLANK.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 June 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN output unit.
!
!    Input, integer ( kind = 4 ) GRID_PX_MIN, GRID_PY_MIN, GRID_PX_MAX, 
!    GRID_PY_MAX, the coordinates of the lower left and upper right corners of 
!    the sudoku grid to be drawn.
!
  implicit none

  integer ( kind = 4 ) grid_px_max
  integer ( kind = 4 ) grid_px_min
  integer ( kind = 4 ) grid_py_max
  integer ( kind = 4 ) grid_py_min
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) j
  integer ( kind = 4 ) px
  integer ( kind = 4 ) py

  do j = 1, 10

    if ( j == 4 .or. j == 7 ) then
      write ( iunit, '(i4,a)' ) 2, ' setlinewidth'
    else if ( j == 5 .or. j == 8 ) then
      write ( iunit, '(i4,a)' ) 1, ' setlinewidth'
    end if

    px = ( ( 10 - j ) * grid_px_min + ( j - 1 ) * grid_px_max ) / 9
    write ( iunit, '(a)' ) ' newpath'
    write ( iunit, '(2i4,a)' ) px, grid_py_min, ' moveto'
    write ( iunit, '(2i4,a)' ) px, grid_py_max, ' lineto'
    write ( iunit, '(a)' ) ' stroke'

  end do

  do i = 1, 10

    if ( i == 4 .or. i == 7 ) then
      write ( iunit, '(i4,a)' ) 2, ' setlinewidth'
    else if ( i == 5 .or. i == 8 ) then
      write ( iunit, '(i4,a)' ) 1, ' setlinewidth'
    end if

    py = ( ( 10 - i ) * grid_py_min + ( i - 1 ) * grid_py_max ) / 9
    write ( iunit, '(a)' ) ' newpath'
    write ( iunit, '(2i4,a)' ) grid_px_min, py, ' moveto'
    write ( iunit, '(2i4,a)' ) grid_px_max, py, ' lineto'
    write ( iunit, '(a)' ) ' stroke'

  end do

  return
end
subroutine sudoku_grid_filled ( iunit, grid_px_min, grid_py_min, grid_px_max, &
  grid_py_max, value )

!*****************************************************************************80
!
!! SUDOKU_GRID_FILLED draws one sudoku grid with some values filled.
!
!  Discussion:
!
!    This routine should be called by SUDOKU_SHEET_FILLED.
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
!    Input, integer ( kind = 4 ) IUNIT, the FORTRAN output unit.
!
!    Input, integer ( kind = 4 ) GRID_PX_MIN, GRID_PY_MIN, GRID_PX_MAX, 
!    GRID_PY_MAX, the coordinates of the lower left and upper right corners of 
!    the sudoku grid to be drawn.
!
!    Input, integer ( kind = 4 ) VALUE(9,9), contains 81 values.  Values between
!    1 and 9 will be plotted.  Other values are ignored.
!
  implicit none

  integer ( kind = 4 ) grid_px_max
  integer ( kind = 4 ) grid_px_min
  integer ( kind = 4 ) grid_py_max
  integer ( kind = 4 ) grid_py_min
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) j
  integer ( kind = 4 ) px
  integer ( kind = 4 ) py
  integer ( kind = 4 ) value(9,9)

  do j = 1, 10

    if ( j == 4 .or. j == 7 ) then
      write ( iunit, '(i4,a)' ) 2, ' setlinewidth'
    else if ( j == 5 .or. j == 8 ) then
      write ( iunit, '(i4,a)' ) 1, ' setlinewidth'
    end if

    px = ( ( 10 - j ) * grid_px_min + ( j - 1 ) * grid_px_max ) / 9
    write ( iunit, '(a)' ) ' newpath'
    write ( iunit, '(2i4,a)' ) px, grid_py_min, ' moveto'
    write ( iunit, '(2i4,a)' ) px, grid_py_max, ' lineto'
    write ( iunit, '(a)' ) ' stroke'

  end do

  do i = 1, 10

    if ( i == 4 .or. i == 7 ) then
      write ( iunit, '(i4,a)' ) 2, ' setlinewidth'
    else if ( i == 5 .or. i == 8 ) then
      write ( iunit, '(i4,a)' ) 1, ' setlinewidth'
    end if

    py = ( ( 10 - i ) * grid_py_min + ( i - 1 ) * grid_py_max ) / 9
    write ( iunit, '(a)' ) ' newpath'
    write ( iunit, '(2i4,a)' ) grid_px_min, py, ' moveto'
    write ( iunit, '(2i4,a)' ) grid_px_max, py, ' lineto'
    write ( iunit, '(a)' ) ' stroke'

  end do

  do i = 1, 9

    do j = 1, 9

      if ( 1 <= value(i,j) .and. value(i,j) <= 9 ) then

        px = ( ( 40 - ( 4 * j + 1 )     ) * grid_px_min   &
             + (      ( 4 * j + 1 ) - 4 ) * grid_px_max ) &
               / 36

        py = ( ( 40 - ( 4 * i + 3 )     ) * grid_py_max   &
             + (      ( 4 * i + 3 ) - 4 ) * grid_py_min ) &
               / 36

        write ( iunit, '(2i4,a)' ) px, py, ' moveto'
        write ( iunit, '(a,i1,a)' ) '(', value(i,j), ') show'

      end if
    end do
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

  character ( len = 8 )  ampm
  integer   ( kind = 4 ) d
  integer   ( kind = 4 ) h
  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
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

  write ( *, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
subroutine triangular_1 ( n, file_name )

!*****************************************************************************80
!
!! TRIANGULAR_1 draws a triangular grid.
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
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, specifies the number of grid lines drawn.
!
!    Input, character ( len = * ) FILE_NAME, the name of the file to which the
!    output should be written.
!
  implicit none

  character ( len = * )  file_name
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) ierror
  integer   ( kind = 4 ) iunit
  integer   ( kind = 4 ) n
  integer   ( kind = 4 ) pax
  integer   ( kind = 4 ) pay
  integer   ( kind = 4 ) pbx
  integer   ( kind = 4 ) pby
  integer   ( kind = 4 ) pcx
  integer   ( kind = 4 ) pcy
  integer   ( kind = 4 ) px
  integer   ( kind = 4 ) py
  integer   ( kind = 4 ) x_ps_max
  integer   ( kind = 4 ) x_ps_min
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  integer   ( kind = 4 ) y_ps_max
  integer   ( kind = 4 ) y_ps_min
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin
!
!  Open the output file.
!
  call get_unit ( iunit )

  call ps_file_open ( file_name, iunit, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRIANGULAR_1 - Fatal error!'
    write ( *, '(a,i8)' ) '  File creation error ', ierror
    stop
  end if

  x_ps_min = 36
  y_ps_min = 36
  x_ps_max = 576
  y_ps_max = 756

  call eps_file_head ( file_name, x_ps_min, y_ps_min, x_ps_max, &
    y_ps_max )
!
!  Write the header.
!
  xmin = 0.0D+00
  ymin = 0.0D+00

  xmax = 8.5D+00
  ymax = 11.0D+00

  call ps_page_head ( xmin, ymin, xmax, ymax )
!
!  Draw the grid lines.
!
  pax =  36
  pay = 136

  pbx = 576
  pby = 136

  pcx = 306
  pcy = 568
  
  do i = 0, n

    write ( iunit, '(a)' ) ' newpath'

    px = int ( real ( ( n - i ) * pax   &
                          + i   * pbx, kind = 8 ) &
             / real   ( n,             kind = 8 ) )

    py = int ( real ( ( n - i ) * pay   &
                          + i   * pby, kind = 8 ) &
             / real   ( n,             kind = 8 ) )

    write ( iunit, '(2i4,a)' ) px, py, ' moveto'

    px = int ( real ( ( n - i ) * pcx   &
                          + i   * pbx, kind = 8 ) &
             / real   ( n,             kind = 8 ) )

    py = int ( real ( ( n - i ) * pcy   &
                          + i   * pby, kind = 8 ) &
             / real   ( n,             kind = 8 ) )

    write ( iunit, '(2i4,a)' ) px, py, ' lineto'

    write ( iunit, '(a)' ) ' stroke'

  end do

  do i = 0, n

    write ( iunit, '(a)' ) ' newpath'

    px = int ( real ( ( n - i ) * pax   &
                          + i   * pcx, kind = 8 ) &
             / real   ( n,             kind = 8 ) )

    py = int ( real ( ( n - i ) * pay   &
                          + i   * pcy, kind = 8 ) &
             / real   ( n,             kind = 8 ) )

    write ( iunit, '(2i4,a)' ) px, py, ' moveto'

    px = int ( real ( ( n - i ) * pbx   &
                          + i   * pcx, kind = 8 ) &
             / real   ( n,             kind = 8 ) )

    py = int ( real ( ( n - i ) * pby   &
                          + i   * pcy, kind = 8 ) &
             / real   ( n,             kind = 8 ) )

    write ( iunit, '(2i4,a)' ) px, py, ' lineto'

    write ( iunit, '(a)' ) ' stroke'

  end do

  do i = 0, n

    write ( iunit, '(a)' ) ' newpath'

    px = int ( real ( ( n - i ) * pbx   &
                          + i   * pax, kind = 8 ) &
             / real   ( n,             kind = 8 ) )

    py = int ( real ( ( n - i ) * pby   &
                          + i   * pay, kind = 8 ) &
             / real   ( n,             kind = 8 ) )

    write ( iunit, '(2i4,a)' ) px, py, ' moveto'

    px = int ( real ( ( n - i ) * pcx   &
                          + i   * pax, kind = 8 ) &
             / real   ( n,             kind = 8 ) )

    py = int ( real ( ( n - i ) * pcy   &
                          + i   * pay, kind = 8 ) &
             / real   ( n,             kind = 8 ) )

    write ( iunit, '(2i4,a)' ) px, py, ' lineto'

    write ( iunit, '(a)' ) ' stroke'

  end do

  call ps_page_tail ( )

  call eps_file_tail ( )

  call ps_file_close ( iunit )

  return
end
subroutine triangular_2 ( n, file_name )

!*****************************************************************************80
!
!! TRIANGULAR_2 draws a triangular grid of dots (no lines).
!
!  Discussion:
!
!    Depending on how you view them, the points on this grid are
!    arranged in triangles or hexagons.
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
!    Input, integer ( kind = 4 ) N, specifies the number of grid lines drawn.
!
!    Input, character ( len = * ) FILE_NAME, the name of the file to which the
!    output should be written.
!
  implicit none

  character ( len = * )  file_name
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) ierror
  integer   ( kind = 4 ) iunit
  integer   ( kind = 4 ) j
  integer   ( kind = 4 ) marker_size
  integer   ( kind = 4 ) n
  integer   ( kind = 4 ) pax
  integer   ( kind = 4 ) pay
  integer   ( kind = 4 ) pbx
  integer   ( kind = 4 ) pby
  integer   ( kind = 4 ) pcx
  integer   ( kind = 4 ) pcy
  integer   ( kind = 4 ) px
  integer   ( kind = 4 ) py
  integer   ( kind = 4 ) x_ps_max
  integer   ( kind = 4 ) x_ps_min
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  integer   ( kind = 4 ) y_ps_max
  integer   ( kind = 4 ) y_ps_min
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin
!
!  Open the output file.
!
  call get_unit ( iunit )

  call ps_file_open ( file_name, iunit, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRIANGULAR_2 - Fatal error!'
    write ( *, '(a,i8)' ) '  File creation error ', ierror
    stop
  end if

  x_ps_min = 36
  y_ps_min = 36
  x_ps_max = 576
  y_ps_max = 756

  call eps_file_head ( file_name, x_ps_min, y_ps_min, x_ps_max, &
    y_ps_max )
!
!  Write the header.
!
  xmin = 0.0D+00
  ymin = 0.0D+00

  xmax = 8.5D+00
  ymax = 11.0D+00

  call ps_page_head ( xmin, ymin, xmax, ymax )
!
!  Draw the points.
!
  pax =  36
  pay = 136

  pbx = 576
  pby = 136

  pcx = 306
  pcy = 568
  
  marker_size = 4

  do i = 0, n

    do j = 0, n-i

      px = int ( real ( ( n - i - j ) * pax   &
                            + i       * pbx   &
                                + j   * pcx, kind = 8 ) &
               / real   ( n,                 kind = 8) )

      py = int ( real ( ( n - i - j ) * pay   &
                            + i       * pby   &
                                + j   * pcy, kind = 8 ) &
               / real   ( n,                 kind = 8 ) )

      write ( iunit, '(a,3i6,a)' ) 'newpath ', px, py, marker_size, &
        ' 0 360 arc closepath fill'

    end do

  end do

  call ps_page_tail ( )

  call eps_file_tail ( )

  call ps_file_close ( iunit )

  return
end
subroutine uniform_1 ( file_name )

!*****************************************************************************80
!
!! UNIFORM_1 draws a uniform grid at 1/4 inch intervals.
!
!  Discussion:
!
!    The graph paper is surrounded by a half inch blank margin.
!    Lines are drawn every 1/4 inch.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the file to which the
!    output should be written.
!
  implicit none

  character ( len = * )  file_name
  integer   ( kind = 4 ) ierror
  integer   ( kind = 4 ) iunit
  integer   ( kind = 4 ) px
  integer   ( kind = 4 ) pxinc
  integer   ( kind = 4 ) pxmax
  integer   ( kind = 4 ) pxmin
  integer   ( kind = 4 ) py
  integer   ( kind = 4 ) pyinc
  integer   ( kind = 4 ) pymax
  integer   ( kind = 4 ) pymin
  integer   ( kind = 4 ) x_ps_max
  integer   ( kind = 4 ) x_ps_min
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  integer   ( kind = 4 ) y_ps_max
  integer   ( kind = 4 ) y_ps_min
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin
!
!  Open the output file.
!
  call get_unit ( iunit )

  call ps_file_open ( file_name, iunit, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'UNIFORM_1 - Fatal error!'
    write ( *, '(a,i8)' ) '  File creation error ', ierror
    stop
  end if

  x_ps_min = 36
  y_ps_min = 36
  x_ps_max = 576
  y_ps_max = 756

  call eps_file_head ( file_name, x_ps_min, y_ps_min, x_ps_max, &
    y_ps_max )
!
!  Write the header.
!
  xmin = 0.0D+00
  ymin = 0.0D+00

  xmax = 8.5D+00
  ymax = 11.0D+00

  call ps_page_head ( xmin, ymin, xmax, ymax )
!
!  Draw the grid lines.
!
  pxmin = 0 + 36
  pxmax = 612 - 36
  pxinc = 18

  pymin = 0 + 36
  pymax = 792 - 36
  pyinc = 18

  do px = pxmin, pxmax, pxinc

    write ( iunit, '(a)' ) ' newpath'
    py = pymin
    write ( iunit, '(2i4,a)' ) px, py, ' moveto'

    py = pymax

    write ( iunit, '(2i4,a)' ) px, py, ' lineto'

    write ( iunit, '(a)' ) ' stroke'

  end do

  do py = pymin, pymax, pyinc

    write ( iunit, '(a)' ) ' newpath'
    px = pxmin
    write ( iunit, '(2i4,a)' ) px, py, ' moveto'

    px = pxmax

    write ( iunit, '(2i4,a)' ) px, py, ' lineto'

    write ( iunit, '(a)' ) ' stroke'

  end do

  call ps_page_tail ( )

  call eps_file_tail ( )

  call ps_file_close ( iunit )

  return
end
subroutine uniform_2 ( n, file_name )

!*****************************************************************************80
!
!! UNIFORM_2 draws a uniform grid of dots (no lines).
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
!    Input, integer ( kind = 4 ) N, specifies the number of grid lines drawn.
!
!    Input, character ( len = * ) FILE_NAME, the name of the file to which the
!    output should be written.
!
  implicit none

  character ( len = * )  file_name
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) ierror
  integer   ( kind = 4 ) iunit
  integer   ( kind = 4 ) j
  integer   ( kind = 4 ) marker_size
  integer   ( kind = 4 ) n
  integer   ( kind = 4 ) px
  integer   ( kind = 4 ) pxmax
  integer   ( kind = 4 ) pxmin
  integer   ( kind = 4 ) py
  integer   ( kind = 4 ) pymax
  integer   ( kind = 4 ) pymin
  integer   ( kind = 4 ) x_ps_max
  integer   ( kind = 4 ) x_ps_min
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  integer   ( kind = 4 ) y_ps_max
  integer   ( kind = 4 ) y_ps_min
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin
!
!  Open the output file.
!
  call get_unit ( iunit )

  call ps_file_open ( file_name, iunit, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'UNIFORM_2 - Fatal error!'
    write ( *, '(a,i8)' ) '  File creation error ', ierror
    stop
  end if

  x_ps_min = 36
  y_ps_min = 36
  x_ps_max = 576
  y_ps_max = 756

  call eps_file_head ( file_name, x_ps_min, y_ps_min, x_ps_max, &
    y_ps_max )
!
!  Write the header.
!
  xmin = 0.0D+00
  ymin = 0.0D+00

  xmax = 8.5D+00
  ymax = 11.0D+00

  call ps_page_head ( xmin, ymin, xmax, ymax )
!
!  Draw the points.
!
  pxmin =   0 + 36
  pxmax = 612 - 36

  pymin =   0 + 90 + 36
  pymax = 792 - 90 - 36
  
  marker_size = 4

  do i = 1, n

    do j = 1, n

      px = int ( real ( ( n - i     ) * pxmin   &
                      + (     i - 1 ) * pxmax, kind = 8 ) &
               / real   ( n     - 1,           kind = 8 ) )

      py = int ( real ( ( n - j     ) * pymin   &
                      + (     j - 1 ) * pymax, kind = 8 ) &
               / real   ( n     - 1,           kind = 8 ) )

      write ( iunit, '(a,3i6,a)' ) 'newpath ', px, py, marker_size, &
        ' 0 360 arc closepath fill'

    end do

  end do

  call ps_page_tail ( )

  call eps_file_tail ( )

  call ps_file_close ( iunit )

  return
end
subroutine uniform_3 ( file_name )

!*****************************************************************************80
!
!! UNIFORM_3 draws a uniform 1/4 inch grid; every fifth line is heavy.
!
!  Discussion:
!
!    The graph paper is surrounded by a half inch blank margin.
!    Lines are drawn every 1/4 inch.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the file to which the
!    output should be written.
!
  implicit none

  character ( len = * ) file_name
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) px
  integer ( kind = 4 ) pxinc
  integer ( kind = 4 ) pxmax
  integer ( kind = 4 ) pxmin
  integer ( kind = 4 ) py
  integer ( kind = 4 ) pyinc
  integer ( kind = 4 ) pymax
  integer ( kind = 4 ) pymin
  integer ( kind = 4 ) thick
  integer ( kind = 4 ) thin
  integer ( kind = 4 ) x_ps_max
  integer ( kind = 4 ) x_ps_min
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  integer ( kind = 4 ) y_ps_max
  integer ( kind = 4 ) y_ps_min
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin
!
!  Open the output file.
!
  call get_unit ( iunit )

  call ps_file_open ( file_name, iunit, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'UNIFORM_3 - Fatal error!'
    write ( *, '(a,i8)' ) '  File creation error ', ierror
    stop
  end if

  x_ps_min = 36
  y_ps_min = 36
  x_ps_max = 576
  y_ps_max = 756

  call eps_file_head ( file_name, x_ps_min, y_ps_min, x_ps_max, &
    y_ps_max )
!
!  Write the header.
!
  xmin = 0.0D+00
  ymin = 0.0D+00

  xmax = 8.5D+00
  ymax = 11.0D+00

  call ps_page_head ( xmin, ymin, xmax, ymax )
!
!  Draw the grid lines.
!
  pxmin = 0 + 36
  pxmax = 612 - 36
  pxinc = 18

  pymin = 0 + 36
  pymax = 792 - 36
  pyinc = 18

  thin = 1
  thick = 2

  do px = pxmin, pxmax, pxinc

    if ( mod ( ( px - pxmin ) / pxinc, 5 ) == 0 ) then
      write ( iunit, '(i4,a)' ) thick, ' setlinewidth'
    else if ( mod ( ( px - pxmin ) / pxinc, 5 ) == 1 ) then
      write ( iunit, '(i4,a)' ) thin, ' setlinewidth'
    end if

    write ( iunit, '(a)' ) ' newpath'
    py = pymin
    write ( iunit, '(2i4,a)' ) px, py, ' moveto'

    py = pymax

    write ( iunit, '(2i4,a)' ) px, py, ' lineto'

    write ( iunit, '(a)' ) ' stroke'

  end do

  do py = pymin, pymax, pyinc

    if ( mod ( ( py - pymin ) / pyinc, 5 ) == 0 ) then
      write ( iunit, '(i4,a)' ) thick, ' setlinewidth'
    else if ( mod ( ( py - pymin ) / pyinc, 5 ) == 1 ) then
      write ( iunit, '(i4,a)' ) thin, ' setlinewidth'
    end if

    write ( iunit, '(a)' ) ' newpath'
    px = pxmin
    write ( iunit, '(2i4,a)' ) px, py, ' moveto'

    px = pxmax

    write ( iunit, '(2i4,a)' ) px, py, ' lineto'

    write ( iunit, '(a)' ) ' stroke'

  end do

  call ps_page_tail ( )

  call eps_file_tail ( )

  call ps_file_close ( iunit )

  return
end
subroutine uniform_4 ( n, file_name )

!*****************************************************************************80
!
!! UNIFORM_4 draws a uniform grid with N boxes per inch.
!
!  Discussion:
!
!    The graph paper is surrounded by a half inch blank margin.
!    Lines are drawn every 1/N inch.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of lines or boxes per inch.
!
!    Input, character ( len = * ) FILE_NAME, the name of the file to which the
!    output should be written.
!
  implicit none

  character ( len = * ) file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  integer ( kind = 4 ) px
  integer ( kind = 4 ) pxinc
  integer ( kind = 4 ) pxmax
  integer ( kind = 4 ) pxmin
  integer ( kind = 4 ) py
  integer ( kind = 4 ) pyinc
  integer ( kind = 4 ) pymax
  integer ( kind = 4 ) pymin
  integer ( kind = 4 ) x_ps_max
  integer ( kind = 4 ) x_ps_min
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  integer ( kind = 4 ) y_ps_max
  integer ( kind = 4 ) y_ps_min
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin
!
!  Open the output file.
!
  call get_unit ( iunit )

  call ps_file_open ( file_name, iunit, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'UNIFORM_1 - Fatal error!'
    write ( *, '(a,i8)' ) '  File creation error ', ierror
    stop
  end if

  x_ps_min = 36
  y_ps_min = 36
  x_ps_max = 576
  y_ps_max = 756

  call eps_file_head ( file_name, x_ps_min, y_ps_min, x_ps_max, &
    y_ps_max )
!
!  Write the header.
!
  xmin = 0.0D+00
  ymin = 0.0D+00

  xmax = 8.5D+00
  ymax = 11.0D+00

  call ps_page_head ( xmin, ymin, xmax, ymax )
!
!  Draw the grid lines.
!
  pxmin = 0 + 36
  pxmax = 612 - 36

  pymin = 0 + 36
  pymax = 792 - 36

  do i = 0, ( 15 * n ) / 2

    px = pxmin + nint ( real ( i * 72, kind = 8 ) / real ( n, kind = 8 ) )

    write ( iunit, '(a)' ) ' newpath'
    py = pymin
    write ( iunit, '(2i4,a)' ) px, py, ' moveto'

    py = pymax

    write ( iunit, '(2i4,a)' ) px, py, ' lineto'

    write ( iunit, '(a)' ) ' stroke'

  end do

  do j = 0, 10 * n

    py = pymin + nint ( real ( j * 72, kind = 8 ) / real ( n, kind = 8 ) )

    write ( iunit, '(a)' ) ' newpath'
    px = pxmin
    write ( iunit, '(2i4,a)' ) px, py, ' moveto'

    px = pxmax

    write ( iunit, '(2i4,a)' ) px, py, ' lineto'

    write ( iunit, '(a)' ) ' stroke'

  end do

  call ps_page_tail ( )

  call eps_file_tail ( )

  call ps_file_close ( iunit )

  return
end
