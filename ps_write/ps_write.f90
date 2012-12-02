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
!    Output, integer IVAL:
!    -1, no part of the line segment is within the box.
!     0, no clipping was necessary.  The line segment is entirely within 
!        the box.
!     1, (X1,Y1) was clipped.
!     2, (X2,Y2) was clipped.
!     3, (X1,Y1) and (X2,Y2) were clipped.
!
  implicit none

  integer ival
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
  integer itemp

  itemp = ichar ( c )

  if ( 97 <= itemp .and. itemp <= 122 ) then
    c = char ( itemp - 32 )
  end if

  return
end
subroutine circle_points ( x0, y0, r, n, x, y )

!*****************************************************************************80
!
!! CIRCLE_POINTS returns N equally spaced points on a circle in 2D.
!
!  Discussion:
!
!    The first point is always ( X0 + R, Y0 ), and subsequent points
!    proceed counterclockwise around the circle.
!
!  Modified:
!
!    28 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X0, Y0, the coordinates of the center of 
!    the circle.
!
!    Input, real ( kind = 8 ) R, the radius of the circle.
!
!    Input, integer N, the number of points desired.  N must be at least 1.
!
!    Output, real ( kind = 8 ) X(N), Y(N), the coordinates of points
!    on the circle.
!
  implicit none

  integer n

  real ( kind = 8 ) angle
  integer i
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r
  real ( kind = 8 ) x0
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y0
  real ( kind = 8 ) y(n)

  do i = 1, n
    angle = ( 2.0D+00 * pi * real ( i - 1, kind = 8 ) ) &
      / real ( n, kind = 8 )
    x(i) = x0 + r * cos ( angle )
    y(i) = y0 + r * sin ( angle )
  end do

  return
end
subroutine circle_points_arc ( x0, y0, r, theta1, theta2, n, x, y )

!*****************************************************************************80
!
!! CIRCLE_POINTS_ARC returns N points on a circular arc in 2D.
!
!  Discussion:
!
!    The first point is ( X0 + R * COS ( THETA1 ), Y0 + R * SIN ( THETA1 ) );
!    The last point is  ( X0 + R * COS ( THETA2 ), Y0 + R * SIN ( THETA2 ) );
!    and the intermediate points are evenly spaced in angle between these,
!    and in counterclockwise order.
!
!  Modified:
!
!    29 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X0, Y0, the coordinates of the center
!    of the circle.
!
!    Input, real ( kind = 8 ) R, the radius of the circle.
!
!    Input, real ( kind = 8 ) THETA1, THETA2, the angular coordinates
!    of the first and last points to be drawn, in radians.
!
!    Input, integer N, the number of points desired.  N must be at least 1.
!
!    Output, real ( kind = 8 ) X(N), Y(N), the coordinates of points 
!    on the circle.
!
  implicit none

  integer n

  integer i
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r
  real ( kind = 8 ) r8_modp
  real ( kind = 8 ) theta
  real ( kind = 8 ) theta1
  real ( kind = 8 ) theta2
  real ( kind = 8 ) theta3
  real ( kind = 8 ) x0
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y0
  real ( kind = 8 ) y(n)
!
!  THETA3 is the smallest angle, no less than THETA1, which
!  coincides with THETA2.
!
  theta3 = theta1 + r8_modp ( theta2 - theta1, 2.0D+00 * pi )

  do i = 1, n

    if ( 1 < n ) then
      theta = ( real ( n - i,     kind = 8 ) * theta1 &
              + real (     i - 1, kind = 8 ) * theta3 ) &
              / real ( n     - 1 )
    else
      theta = 0.5D+00 * ( theta1 + theta3 )
    end if

    x(i) = x0 + r * cos ( theta )
    y(i) = y0 + r * sin ( theta )

  end do

  return
end
function degrees_to_radians ( angle )

!*****************************************************************************80
!
!! DEGREES_TO_RADIANS converts an angle from degrees to radians.
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
subroutine ellipse_points ( x0, y0, r1, r2, theta, n, x, y )

!*****************************************************************************80
!
!! ELLIPSE_POINTS returns N points on an ellipse in 2D.
!
!  Discussion:
!
!    The points are "equally spaced" in the angular sense.  They are
!    not equally spaced along the perimeter of the ellipse.
!
!  Modified:
!
!    28 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X0, Y0, the coordinates of the center 
!    of the ellipse.
!
!    Input, real ( kind = 8 ) R1, R2, the "radius" of the ellipse in the major
!    and minor axis directions.  A circle has these values equal.
!
!    Input, real ( kind = 8 ) THETA, the angle that the major axis of the
!    ellipse makes with the X axis.  A value of 0.0 means that the major and
!    minor axes of the ellipse will be the X and Y coordinate axes.
!
!    Input, integer N, the number of points desired.  N must be at least 1.
!
!    Output, real ( kind = 8 ) X(N), Y(N), the coordinates of points on 
!    the ellipse.
!
  implicit none

  integer n

  real ( kind = 8 ) angle
  integer i
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) theta
  real ( kind = 8 ) x0
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y0
  real ( kind = 8 ) y(n)

  do i = 1, n

    angle = ( 2.0D+00 * pi * real ( i - 1, kind = 8 ) ) &
      / real ( n, kind = 8 )

    x(i) = x0 + r1 * cos ( theta ) * cos ( angle ) &
              - r2 * sin ( theta ) * sin ( angle )

    y(i) = y0 + r1 * sin ( theta ) * cos ( angle ) &
              + r2 * cos ( theta ) * sin ( angle )

  end do

  return
end
subroutine ellipse_points_arc ( x0, y0, r1, r2, psi, theta1, theta2, n, x, y )

!*****************************************************************************80
!
!! ELLIPSE_POINTS_ARC returns N points on an elliptical arc in 2D.
!
!  Discussion:
!
!    The points are "equally spaced" in the angular sense.  They are
!    not equally spaced along the perimeter of the ellipse.
!
!  Modified:
!
!    29 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X0, Y0, the coordinates of the center of 
!    the ellipse.
!
!    Input, real ( kind = 8 ) R1, R2, the "radius" of the ellipse in the major
!    and minor axis directions.  A circle has these values equal.
!
!    Input, real ( kind = 8 ) PSI, the angle that the major axis of the ellipse
!    makes with the X axis.  A value of 0.0 means that the major and
!    minor axes of the ellipse will be the X and Y coordinate axes.
!
!    Input, real ( kind = 8 ) THETA1, THETA2, the angular coordinates of the
!    first and last points to be drawn, in radians.  This angle is measured
!    with respect to the (possibly tilted) major axis.
!
!    Input, integer N, the number of points desired.  N must be at least 1.
!
!    Output, real ( kind = 8 ) X(N), Y(N), the coordinates of points 
!    on the ellipse.
!
  implicit none

  integer n

  integer i
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) psi
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) r8_modp
  real ( kind = 8 ) theta
  real ( kind = 8 ) theta1
  real ( kind = 8 ) theta2
  real ( kind = 8 ) theta3
  real ( kind = 8 ) x0
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y0
  real ( kind = 8 ) y(n)
!
!  THETA3 is the smallest angle, no less than THETA1, which
!  coincides with THETA2.
!
  theta3 = theta1 + r8_modp ( theta2 - theta1, 2.0D+00 * pi )

  do i = 1, n

    if ( 1 < n ) then
      theta = ( real ( n - i,     kind = 8 ) * theta1 &
              + real (     i - 1, kind = 8 ) * theta3 ) &
              / real ( n     - 1, kind = 8 )
    else
      theta = 0.5D+00 * ( theta1 + theta3 )
    end if

    x(i) = x0 + r1 * cos ( psi ) * cos ( theta ) &
              - r2 * sin ( psi ) * sin ( theta )

    y(i) = y0 + r1 * sin ( psi ) * cos ( theta ) &
              + r2 * cos ( psi ) * sin ( theta )

  end do

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
!    Henry McGilton, Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the output file.
!
!    Input, integer X_PS_MIN, Y_PS_MIN, X_PS_MAX, Y_PS_MAX, the minimum 
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
  integer state
  integer unit
  integer x_ps_max
  integer x_ps_min
  integer y_ps_max
  integer y_ps_min
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

  integer num_pages
  integer state
  integer unit
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
!    Output, integer IUNIT, the free unit number.
!
  implicit none

  integer i
  integer ios
  integer iunit
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
!    Input/output, integer I, J.  On output, the values of I and
!    J have been interchanged.
!
  implicit none

  integer i
  integer j
  integer k

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
!  Definition:
!
!    A "box" is defined by its "left down" corner and its
!    "right up" corner, and all the points between.  It is
!    assumed that the sides of the box align with coordinate directions.
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
!    Henry McGilton, Mary Campione,
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
  integer ival
  integer plotxmin2
  integer plotymin2
  integer px
  integer py
  real ( kind = 8 ) r
  real ( kind = 8 ) r2
  integer state
  integer unit
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
  call box_clip_line_2d ( xmin, ymin, xmax, ymax, x1, y1, x2, y2, x3, y3, &
    x4, y4, ival )

  if ( ival < 0 ) then
    return
  end if

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
subroutine ps_circle ( x0, y0, r )

!*****************************************************************************80
!
!! PS_CIRCLE draws a circle.
!
!  Discussion:
!
!    As a side effect, the current point is set to the center of the circle.
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
!    Input, real ( kind = 8 ) X0, Y0, the coordinates of the center
!    of the circle.
!
!    Input, real ( kind = 8 ) R, the radius of the circle.
!
  implicit none

  real ( kind = 8 ) alpha
  integer, parameter :: angle_max = 360
  integer, parameter :: angle_min = 0
  integer plotxmin2
  integer plotymin2
  integer pr
  integer pxcen
  integer pycen
  real ( kind = 8 ) r
  integer state
  integer unit
  real ( kind = 8 ) x0
  real ( kind = 8 ) xmin
  real ( kind = 8 ) y0
  real ( kind = 8 ) ymin
!
!  Determine if the PostScript state is acceptable.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state /= 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_CIRCLE - Fatal error!'
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

  pxcen = plotxmin2 + nint ( alpha * ( x0 - xmin ) )
  pycen = plotymin2 + nint ( alpha * ( y0 - ymin ) )
  pr = nint ( alpha * r )

  write ( unit, '(a)' ) 'newpath'
  write ( unit, '(5i6,a)' ) pxcen, pycen, pr, angle_min, angle_max, ' arc'
!
!  Draw the circle.
!
  write ( unit, '(a)' ) 'closepath stroke'

  call ps_setting_real ( 'SET', 'XCUR', x0 )
  call ps_setting_real ( 'SET', 'YCUR', y0 )

  return
end
subroutine ps_circle_arc ( x0, y0, r, theta1, theta2 )

!*****************************************************************************80
!
!! PS_CIRCLE_ARC draws a circular arc.
!
!  Discussion:
!
!    As a side effect, the current point is set to the center of the circle.
!
!  Modified:
!
!    28 May 1999
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
!    Input, real ( kind = 8 ) X0, Y0, the coordinates of the center
!    of the circle.
!
!    Input, real ( kind = 8 ) R, the radius of the circle.
!
!    Input, real ( kind = 8 ) THETA1, THETA2, the angular coordinates
!    of the first and last points on the circular arc to be drawn.  These
!    should be ordered in counter-clockwise order.
!
  implicit none

  real ( kind = 8 ) alpha
  integer angle_max
  integer angle_min
  integer plotxmin2
  integer plotymin2
  integer pr
  integer pxcen
  integer pycen
  real ( kind = 8 ) r
  integer state
  real ( kind = 8 ) theta1
  real ( kind = 8 ) theta2
  integer unit
  real ( kind = 8 ) x0
  real ( kind = 8 ) xmin
  real ( kind = 8 ) y0
  real ( kind = 8 ) ymin
!
!  Determine if the PostScript state is acceptable.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state /= 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_CIRCLE_ARC - Fatal error!'
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

  pxcen = plotxmin2 + nint ( alpha * ( x0 - xmin ) )
  pycen = plotymin2 + nint ( alpha * ( y0 - ymin ) )
  pr = nint ( alpha * r )
  angle_min = nint ( theta1 )
  angle_max = nint ( theta2 )

  write ( unit, '(a)' ) 'newpath'
  write ( unit, '(5i6,a)' ) pxcen, pycen, pr, angle_min, angle_max, ' arc'
!
!  Draw.
!
  write ( unit, '(a)' ) 'stroke'

  call ps_setting_real ( 'SET', 'XCUR', x0 )
  call ps_setting_real ( 'SET', 'YCUR', y0 )

  return
end
subroutine ps_circle_fill ( x0, y0, r )

!*****************************************************************************80
!
!! PS_CIRCLE_FILL draws a filled circle.
!
!  Modified:
!
!    04 October 2000
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
!    Input, real ( kind = 8 ) X0, Y0, the coordinates of the center of the disk.
!
!    Input, real ( kind = 8 ) R, the radius of the disk.
!
  implicit none

  integer, parameter :: n = 32

  real ( kind = 8 ) r
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x0
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) y0

  call circle_points ( x0, y0, r, n, x, y )

  call ps_polygon_fill ( n, x, y )

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
!    Henry McGilton, Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    Input, integer NPOINT, the number of points in the clipping polygon.
!
!    Input, real ( kind = 8 ) X(NPOINT), Y(NPOINT), the X and Y components
!    of the points.
!
  implicit none

  integer npoint

  real ( kind = 8 ) alpha
  integer i
  integer plotxmin2
  integer plotymin2
  integer px
  integer py
  integer state
  integer unit
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
!  Fill the polygon.
!
  write ( unit, '(a)' ) 'clip newpath'

  return
end
subroutine ps_color_fill_set ( r, g, b )

!*****************************************************************************80
!
!! PS_COLOR_FILL_SET sets the fill color.
!
!  Discussion:
!
!    By calling this routine, you guarantee that a check will be made
!    of the current fill color.  If the current and new fill colors are
!    the same, then we skip the extraneous action of setting the color.
!
!  Modified:
!
!    11 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, G, B, the RGB values for the new fill color.
!
  implicit none

  real ( kind = 8 ) b
  real ( kind = 8 ) b_old
  real ( kind = 8 ) g
  real ( kind = 8 ) g_old
  real ( kind = 8 ) r
  real ( kind = 8 ) r_old
  integer state
  integer unit
!
!  Check the state.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_COLOR_FILL_SET - Fatal error!'
    write ( *, '(a,i9)' ) '  PostScript state is ', state
    write ( *, '(a)' ) '  A PostScript state of 1 or more is required.'
    return
  end if
!
!  Get the current colors.
!
  call ps_setting_real ( 'GET', 'FILL_RED', r_old )
  call ps_setting_real ( 'GET', 'FILL_GREEN', g_old )
  call ps_setting_real ( 'GET', 'FILL_BLUE', b_old )
!
!  If any color has changed, we need to reset them.
!
  if ( r_old /= r .or. g_old /= g .or. b_old /= b ) then

    call ps_setting_int ( 'GET', 'UNIT', unit )

    call ps_comment ( 'Set RGB line color.' )

    write ( unit, '(3f7.4,a)' ) r, g, b, ' setrgbcolor'

    call ps_setting_real ( 'SET', 'FILL_RED', r )
    call ps_setting_real ( 'SET', 'FILL_GREEN', g )
    call ps_setting_real ( 'SET', 'FILL_BLUE', b )

  end if

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

  integer, parameter :: nstack = 10

  character ( len = * ) action
  real ( kind = 8 ) b
  real ( kind = 8 ) b_old
  real ( kind = 8 ), save, dimension ( nstack) :: b_stack
  real ( kind = 8 ) g
  real ( kind = 8 ) g_old
  real ( kind = 8 ), save, dimension ( nstack) :: g_stack
  integer, save :: istack = 0
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
  integer state
  integer unit
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
  integer unit
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
  integer line_width
  integer marker_size

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
!    Input, integer UNIT, the FORTRAN unit to which output was written.
!
  implicit none

  integer state
  integer unit
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
  integer margin
  integer pagexmax
  integer pagexmin
  integer pageymax
  integer pageymin
  integer plotxmax
  integer plotxmin
  integer plotymax
  integer plotymin
  integer state
  integer unit
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
!    Input, integer UNIT, the FORTRAN unit to which output should
!    be written.
!
!    Input, character ( len = 80 ) FILE_NAME, the name of the output file.
!
!    Output, integer IERROR, error flag.
!    0, no error.
!    nonzero, the file could not be created.
!
  implicit none

  character ( len = * ) file_name
  integer ierror
  integer ios
  integer state
  integer unit
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

  integer num_pages
  integer state
  integer unit
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
subroutine ps_fill_gray ( fill_gray )

!*****************************************************************************80
!
!! PS_FILL_GRAY sets the gray fill for polygons.
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
!    Input, real ( kind = 8 ) FILL_GRAY, the gray fill used to fill polygons.
!    0.0 is black, 1.0 is white, and values in between represent
!    shades of gray.
!
  implicit none

  real ( kind = 8 ) fill_gray
  integer unit

  call ps_setting_int ( 'GET', 'UNIT', unit )

  write ( unit, '(f8.4,a)' ) fill_gray, ' setgray'

  call ps_setting_real ( 'SET', 'FILL_BLUE', fill_gray )
  call ps_setting_real ( 'SET', 'FILL_GREEN', fill_gray )
  call ps_setting_real ( 'SET', 'FILL_RED', fill_gray )

  return
end
subroutine ps_font_size ( font_size )

!*****************************************************************************80
!
!! PS_FONT_SIZE sets the font size.
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
!    Henry McGilton, Mary Campione,
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
  integer state
  integer unit
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
subroutine ps_grid_cartesian ( xmin, xmax, nx, ymin, ymax, ny )

!*****************************************************************************80
!
!! PS_GRID_CARTESIAN draws a cartesian grid.
!
!  Discussion:
!
!    The current point is not modified by this call.
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
!    Input, real ( kind = 8 ) XMIN, XMAX, the minimum and maximum values
!    at which X grid lines should be drawn.
!
!    Input, integer NX, the number of X grid lines.
!    If NX is not positive, no X grid lines are drawn.
!    If NX is 1, a single grid line is drawn midway.
!
!    Input, real ( kind = 8 ) YMIN, YMAX, the minimum and maximum values 
!    at which Y grid lines should be drawn.
!
!    Input, integer NY, the number of Y grid lines.
!    If NY is not positive, no Y grid lines are drawn.
!    If NY is 1, a single grid line is drawn midway.
!
  implicit none

  real ( kind = 8 ) alpha
  integer i
  integer nx
  integer ny
  integer plotxmin2
  integer plotymin2
  integer px
  integer py
  integer state
  integer unit
  real ( kind = 8 ) x
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) xmin2
  real ( kind = 8 ) y
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin
  real ( kind = 8 ) ymin2
!
!  At least one of NX and NY must be positive.
!
  if ( nx < 1 .and. ny < 1 ) then
    return
  end if
!
!  Determine if the PostScript state is acceptable.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state /= 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_GRID_CARTESIAN - Fatal error!'
    write ( *, '(a,i9)' ) '  PostScript state is ', state
    write ( *, '(a)' ) '  PostScript state 3 is required.'
    return
  end if
!
!  Get the unit number.
!
  call ps_setting_int ( 'GET', 'UNIT', unit )
!
!  Get settings.
!
  alpha = 0.0D+00
  xmin2 = 0.0D+00
  ymin2 = 0.0D+00

  call ps_setting_int ( 'GET', 'PXMIN', plotxmin2 )

  call ps_setting_int ( 'GET', 'PYMIN', plotymin2 )

  call ps_setting_real ( 'GET', 'ALPHA', alpha )
  call ps_setting_real ( 'GET', 'XMIN', xmin2 )
  call ps_setting_real ( 'GET', 'YMIN', ymin2 )
!
!  Draw the vertical (X) grid lines.
!
  do i = 1, nx

    if ( 1 < nx ) then

      x = ( real ( nx - i,     kind = 8 ) * xmin   &
          + real (      i - 1, kind = 8 ) * xmax ) &
          / real ( nx     - 1, kind = 8 )

    else if ( nx == 1 ) then

      x = 0.5D+00 * ( xmin + xmax )

    end if

    px = plotxmin2 + nint ( alpha * ( x - xmin2 ) )

    write ( unit, '(a)' ) 'newpath'

    py = plotymin2 + nint ( alpha * ( ymin - ymin2 ) )
    write ( unit, '(2i6,a)' ) px, py, ' moveto'

    py = plotymin2 + nint ( alpha * ( ymax - ymin2 ) )
    write ( unit, '(2i6,a)' ) px, py, ' lineto'

    write ( unit, '(a)' ) 'stroke'

  end do
!
!  Draw the horizontal (Y) grid lines.
!
  do i = 1, ny

    if ( 1 < ny ) then

      y = ( real ( ny - i,     kind = 8 ) * ymin   &
          + real (      i - 1, kind = 8 ) * ymax ) &
          / real ( ny     - 1, kind = 8 )

    else if ( ny == 1 ) then

      y = 0.5D+00 * ( ymin + ymax )

    end if

    py = plotymin2 + nint ( alpha * ( y - ymin2 ) )

    write ( unit, '(a)' ) 'newpath'

    px = plotxmin2 + nint ( alpha * ( xmin - xmin2 ) )
    write ( unit, '(2i6,a)' ) px, py, ' moveto'

    px = plotxmin2 + nint ( alpha * ( xmax - xmin2 ) )
    write ( unit, '(2i6,a)' ) px, py, ' lineto'

    write ( unit, '(a)' ) 'stroke'

  end do

  return
end
subroutine ps_grid_polar ( x0, y0, nr, r1, r2, nt, theta1, theta2 )

!*****************************************************************************80
!
!! PS_GRID_POLAR draws a polar grid.
!
!  Discussion:
!
!    The current point is not modified by this call.
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
!    Input, real ( kind = 8 ) X0, Y0, the coordinates of the origin.
!
!    Input, integer NR, the number of (circular) grid lines to draw.
!
!    Input, real ( kind = 8 ) R1, R2, the minimum and maximum radii at which
!    a grid line is to be drawn.
!
!    Input, integer NT, the number of grid lines in the angular directions.
!    These are rays emanating from the origin, although only the portion
!    between RMIN and RMAX will be drawn.
!
!    Input, real ( kind = 8 ) THETA1, THETA2, the minimum and maximum angles
!    at which a grid line is to be drawn.  These angles are measured
!    in DEGREES.
!
  implicit none

  real ( kind = 8 ) alpha
  real ( kind = 8 ) degrees_to_radians
  integer i
  integer nr
  integer nt
  integer plotxmin2
  integer plotymin2
  real ( kind = 8 ) psi
  integer px
  integer py
  real ( kind = 8 ) r
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  integer state
  real ( kind = 8 ) theta
  real ( kind = 8 ) theta1
  real ( kind = 8 ) theta2
  integer unit
  real ( kind = 8 ) x
  real ( kind = 8 ) x0
  real ( kind = 8 ) xmin2
  real ( kind = 8 ) y
  real ( kind = 8 ) y0
  real ( kind = 8 ) ymin2
!
!  At least one of NR and NT must be positive.
!
  if ( nr < 1 .and. nt < 1 ) then
    return
  end if
!
!  Determine if the PostScript state is acceptable.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state /= 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_GRID_POLAR - Fatal error!'
    write ( *, '(a,i9)' ) '  PostScript state is ', state
    write ( *, '(a)' ) '  PostScript state 3 is required.'
    return
  end if
!
!  Get the unit number.
!
  call ps_setting_int ( 'GET', 'UNIT', unit )
!
!  Get settings.
!
  alpha = 0.0D+00
  xmin2 = 0.0D+00
  ymin2 = 0.0D+00

  call ps_setting_int ( 'GET', 'PXMIN', plotxmin2 )

  call ps_setting_int ( 'GET', 'PYMIN', plotymin2 )

  call ps_setting_real ( 'GET', 'ALPHA', alpha )
  call ps_setting_real ( 'GET', 'XMIN', xmin2 )
  call ps_setting_real ( 'GET', 'YMIN', ymin2 )
!
!  Draw the circular grid lines.
!
  do i = 1, nr

    if ( nr == 1 ) then
      r = 0.5D+00 * ( r1 + r2 )
    else
      r = ( real ( nr - i,     kind = 8 ) * r1   &
          + real (      i - 1, kind = 8 ) * r2 ) &
          / real ( nr     - 1, kind = 8 )
    end if

    if ( 0.0D+00 < r ) then
      call ps_circle_arc ( x0, y0, r, theta1, theta2 )
    end if

  end do
!
!  Draw the radial grid lines.
!
  do i = 1, nt

    if ( nt == 1 ) then
      theta = 0.5D+00 * ( theta1 + theta2 )
    else
      theta = ( real ( nt - i,     kind = 8 ) * theta1   &
              + real (      i - 1, kind = 8 ) * theta2 ) &
              / real ( nt     - 1, kind = 8 )
    end if

    psi = degrees_to_radians ( theta )

    x = x0 + r1 * cos ( psi )
    y = y0 + r1 * sin ( psi )

    px = plotxmin2 + nint ( alpha * ( x - xmin2 ) )
    py = plotymin2 + nint ( alpha * ( y - ymin2 ) )

    write ( unit, '(a)' ) 'newpath'
    write ( unit, '(2i6,a)' ) px, py, ' moveto'

    x = x0 + r2 * cos ( psi )
    y = y0 + r2 * sin ( psi )

    px = plotxmin2 + nint ( alpha * ( x - xmin2 ) )
    py = plotymin2 + nint ( alpha * ( y - ymin2 ) )

    write ( unit, '(2i6,a)' ) px, py, ' lineto'
    write ( unit, '(a)' ) 'stroke'

  end do

  return
end
subroutine ps_grid_triangular ( x1, y1, x2, y2, x3, y3, n )

!*****************************************************************************80
!
!! PS_GRID_TRIANGULAR draws a simple triangular grid.
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
!    Henry McGilton, Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    Input, input, real ( kind = 8 ) X1, Y1, X2, Y2, X3, Y3, the coordinates
!    of the three corners of the grid.
!
!    Input, integer N, the number of grid lines to draw between each
!    pair of corners.
!
  implicit none

  real ( kind = 8 ) alpha
  integer i
  integer n
  integer pax
  integer pay
  integer pbx
  integer pby
  integer pcx
  integer pcy
  integer plotxmin2
  integer plotymin2
  integer px
  integer py
  integer state
  integer unit
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ) xmin
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) y3
  real ( kind = 8 ) ymin
!
!  Determine if the PostScript state is acceptable.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state /= 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_GRID_TRIANGULAR - Fatal error!'
    write ( *, '(a,i9)' ) '  PostScript state is ', state
    write ( *, '(a)' ) '  PostScript state 3 is required.'
    return
  end if
!
!  Get the unit number.
!
  call ps_setting_int ( 'GET', 'UNIT', unit )
!
!  Get settings.
!
  alpha = 0.0D+00
  xmin = 0.0D+00
  ymin = 0.0D+00

  call ps_setting_int ( 'GET', 'PXMIN', plotxmin2 )

  call ps_setting_int ( 'GET', 'PYMIN', plotymin2 )

  call ps_setting_real ( 'GET', 'ALPHA', alpha )
  call ps_setting_real ( 'GET', 'XMIN', xmin )
  call ps_setting_real ( 'GET', 'YMIN', ymin )
!
!  Get the PostScript coordinates of the corners.
!
  pax = plotxmin2 + nint ( alpha * ( x1 - xmin ) )
  pay = plotymin2 + nint ( alpha * ( y1 - ymin ) )
  pbx = plotxmin2 + nint ( alpha * ( x2 - xmin ) )
  pby = plotymin2 + nint ( alpha * ( y2 - ymin ) )
  pcx = plotxmin2 + nint ( alpha * ( x3 - xmin ) )
  pcy = plotymin2 + nint ( alpha * ( y3 - ymin ) )

  do i = 0, n + 1

    write ( unit, '(a)' ) 'newpath'

    px = int ( &
      ( real ( n + 1 - i, kind = 8 ) * pax   &
      + real (         i, kind = 8 ) * pbx ) &
      / real ( n + 1,     kind = 8 ) )
    py = int ( &
      ( real ( n + 1 - i, kind = 8 ) * pay   &
      + real (         i, kind = 8 ) * pby ) &
      / real ( n + 1,     kind = 8 ) )
    write ( unit, '(2i6,a)' ) px, py, ' moveto'

    px = int ( &
      ( real ( n + 1 - i, kind = 8 ) * pcx &
      + real (         i, kind = 8 ) * pbx ) &
      / real ( n + 1,     kind = 8 ) )
    py = int ( &
      ( real ( n + 1 - i, kind = 8 ) * pcy & 
      + real (         i, kind = 8 ) * pby ) &
      / real ( n + 1,     kind = 8 ) )
    write ( unit, '(2i6,a)' ) px, py, ' lineto'

    write ( unit, '(a)' ) 'stroke'

  end do

  do i = 0, n + 1

    write ( unit, '(a)' ) 'newpath'

    px = int ( &
      ( real ( n + 1 - i, kind = 8  ) * pax &
      + real (         i, kind = 8  ) * pcx ) &
      / real ( n + 1,     kind = 8  ) )
    py = int ( &
      ( real ( n + 1 - i, kind = 8  ) * pay  &
      + real (         i, kind = 8  ) * pcy ) & 
      / real ( n + 1,     kind = 8  ) )
    write ( unit, '(2i6,a)' ) px, py, ' moveto'

    px = int ( &
      ( real ( n + 1 - i, kind = 8  ) * pbx &
      + real (         i, kind = 8  ) * pcx ) &
      / real ( n + 1,     kind = 8  ) )
    py = int ( &
      ( real ( n + 1 - i, kind = 8  ) * pby &
      + real (         i, kind = 8  ) * pcy ) &
      / real ( n + 1,     kind = 8  ) )
    write ( unit, '(2i6,a)' ) px, py, ' lineto'

    write ( unit, '(a)' ) 'stroke'

  end do

  do i = 0, n + 1

    write ( unit, '(a)' ) 'newpath'

    px = int ( &
      ( real ( n + 1 - i, kind = 8  ) * pbx   &
      + real (         i, kind = 8  ) * pax ) &
      / real ( n + 1,     kind = 8  ) )
    py = int ( &
      ( real ( n + 1 - i, kind = 8  ) * pby &
      + real (         i, kind = 8  ) * pay ) &
      / real ( n + 1,     kind = 8  ) )
    write ( unit, '(2i6,a)' ) px, py, ' moveto'

    px = int ( &
      ( real ( n + 1 - i, kind = 8  ) * pcx   &
      + real (         i, kind = 8  ) * pax ) &
      / real ( n + 1,     kind = 8  ) )
    py = int ( &
      ( real ( n + 1 - i, kind = 8  ) * pcy   &
      + real (         i, kind = 8  ) * pay ) &
      / real ( n + 1,     kind = 8  ) )
    write ( unit, '(2i6,a)' ) px, py, ' lineto'

    write ( unit, '(a)' ) 'stroke'

  end do

  return
end
subroutine ps_label ( string )

!*****************************************************************************80
!
!! PS_LABEL prints a label at the current position.
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
!    Henry McGilton, Mary Campione,
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
  integer unit

  if ( len_trim ( string ) <= 0 ) then
    return
  end if

  call ps_setting_int ( 'GET', 'UNIT', unit )

  write ( unit, '(a)' ) '(' // trim ( string ) // ') show'

  return
end
subroutine ps_label_slant ( string, angle )

!*****************************************************************************80
!
!! PS_LABEL_SLANT prints a slanted label at a given position.
!
!  Modified:
!
!    03 August 1999
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
!    Input, character ( len = * ) STRING, the string to be printed.
!
!    Input, real ( kind = 8 ) ANGLE, the angle of rotation, in degrees.
!
  implicit none

  real ( kind = 8 ) angle
  character ( len = * ) string
  integer unit

  if ( len_trim ( string ) <= 0 ) then
    return
  end if

  call ps_setting_int ( 'GET', 'UNIT', unit )

  write ( unit, '(f8.4, a)' ) angle, ' rotate'

  write ( unit, '(a)' ) '(' // trim ( string ) // ') show'

  write ( unit, '(f8.4, a)' ) -angle, ' rotate'

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
!    Henry McGilton, Mary Campione,
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
  integer ival
  integer plotxmin2
  integer plotymin2
  integer px
  integer py
  integer state
  integer unit
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
  write ( unit, '(2i6,a)' ) px, py, ' moveto'

  px = plotxmin2 + nint ( alpha * ( x4 - xmin ) )
  py = plotymin2 + nint ( alpha * ( y4 - ymin ) )
  write ( unit, '(2i6,a)' ) px, py, ' lineto stroke'

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
!    Input, integer NPOINT, the number of points in the line.
!
!    Input, real ( kind = 8 ) X(NPOINT), Y(NPOINT), the X and Y components
!    of the points.
!
  implicit none

  integer npoint

  real ( kind = 8 ) alpha
  integer i
  integer plotxmin2
  integer plotymin2
  integer px
  integer py
  integer state
  integer unit
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
!    Input, integer NPOINT, the number of points in the line.
!
!    Input, real ( kind = 8 ) X(NPOINT), Y(NPOINT), the X and Y 
!    components of the points.
!
  implicit none

  integer npoint

  real ( kind = 8 ) alpha
  integer i
  integer plotxmin2
  integer plotymin2
  integer px
  integer py
  integer state
  integer unit
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
subroutine ps_line_width ( line_width )

!*****************************************************************************80
!
!! PS_LINE_WIDTH sets the line width.
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
!    Henry McGilton, Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    Input, integer LINE_WIDTH, the line width.
!    0 is a valid input, and usually produces the thinnest possible line.
!    1 is a more usual line, 2 is thicker, and so on.
!
  implicit none

  integer line_width
  integer state
  integer unit
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

  write ( unit, '(i6,a)' ) line_width, ' setlinewidth'

  call ps_setting_int ( 'SET', 'LINE_WIDTH', line_width )

  return
end
subroutine ps_lineto ( x, y )

!*****************************************************************************80
!
!! PS_LINETO draws a line from the current point to the given point.
!
!  Discussion:
!
!    The current point is updated to the given point.
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
!    Input, real ( kind = 8 ) X, Y, the X and Y components of the new point.
!
  implicit none

  real ( kind = 8 ) alpha
  integer plotxmin2
  integer plotymin2
  integer px
  integer py
  integer state
  integer unit
  real ( kind = 8 ) x
  real ( kind = 8 ) xcur
  real ( kind = 8 ) xmin
  real ( kind = 8 ) y
  real ( kind = 8 ) ycur
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
  call ps_setting_real ( 'GET', 'XCUR', xcur )
  call ps_setting_real ( 'GET', 'XMIN', xmin )
  call ps_setting_real ( 'GET', 'YCUR', ycur )
  call ps_setting_real ( 'GET', 'YMIN', ymin )
!
!  Draw the line.
!
  write ( unit, '(a)' ) 'newpath'

  px = plotxmin2 + nint ( alpha * ( xcur - xmin ) )
  py = plotymin2 + nint ( alpha * ( ycur - ymin ) )
  write ( unit, '(2i6,a)' ) px, py, ' moveto'

  px = plotxmin2 + nint ( alpha * ( x - xmin ) )
  py = plotymin2 + nint ( alpha * ( y - ymin ) )
  write ( unit, '(2i6,a)' ) px, py, ' lineto'
!
!  Draw the line.
!
  write ( unit, '(a)' ) 'stroke'

  call ps_setting_real ( 'SET', 'XCUR', x )
  call ps_setting_real ( 'SET', 'YCUR', y )

  return
end
subroutine ps_mark_circle ( x, y )

!*****************************************************************************80
!
!! PS_MARK_CIRCLE marks a point with a small open circle.
!
!  Discussion:
!
!    The current point is set to the center of the circle.
!
!    The circle is drawn with the current RGB line colors.
!
!    The circle is drawn the current marker size.
!
!    If the point is outside the region, the command is ignored.
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
!    Input, real ( kind = 8 ) X, Y, the coordinates of the point to mark.
!
  implicit none

  real ( kind = 8 ) alpha
  integer marker_size
  integer plotxmin2
  integer plotymin2
  logical point_inside_box_2d
  integer pxcen
  integer pycen
  integer state
  integer unit
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
    write ( *, '(a)' ) 'PS_MARK_CIRCLE - Fatal error!'
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

  write ( unit, '(3i6,a)' ) pxcen, pycen, marker_size, &
    ' 0 360 arc closepath stroke'

  call ps_setting_real ( 'SET', 'XCUR', x )
  call ps_setting_real ( 'SET', 'YCUR', y )

  return
end
subroutine ps_mark_circles ( n, x, y )

!*****************************************************************************80
!
!! PS_MARK_CIRCLES marks points with a small open circle.
!
!  Discussion:
!
!    The current point is set to the center of the last circle.
!
!    The circles are drawn with the current RGB line colors.
!
!    The circles are drawn the current marker size.
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
!    Henry McGilton, Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    Input, integer N, the number of points.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the coordinates of the points to mark.
!
  implicit none

  integer n

  real ( kind = 8 ) alpha
  integer i
  integer marker_size
  integer plotxmin2
  integer plotymin2
  logical point_inside_box_2d
  integer pxcen
  integer pycen
  integer state
  integer unit
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
    write ( *, '(a)' ) 'PS_MARK_CIRCLE - Fatal error!'
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

  write ( unit, '(a)' ) 'newpath'

  do i = 1, n

    if ( .not. point_inside_box_2d ( xmin, ymin, xmax, ymax, x(i), y(i) ) ) then
      cycle
    end if

    pxcen = plotxmin2 + nint ( alpha * ( x(i) - xmin ) )
    pycen = plotymin2 + nint ( alpha * ( y(i) - ymin ) )
    write ( unit, '(3i6,a)' ) pxcen, pycen, marker_size, &
      ' 0 360 arc closepath stroke'

  end do

  call ps_setting_real ( 'SET', 'XCUR', x )
  call ps_setting_real ( 'SET', 'YCUR', y )

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
!    Henry McGilton, Mary Campione,
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
  integer marker_size
  integer plotxmin2
  integer plotymin2
  logical point_inside_box_2d
  integer pxcen
  integer pycen
  integer state
  integer unit
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

  write ( unit, '(3i6,a)' ) pxcen, pycen, marker_size, &
    ' 0 360 arc closepath fill'

  call ps_setting_real ( 'SET', 'XCUR', x )
  call ps_setting_real ( 'SET', 'YCUR', y )

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
!    Input, integer N, the number of points.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the coordinates of the point to mark.
!
  implicit none

  integer n

  real ( kind = 8 ) alpha
  integer i
  integer marker_size
  integer plotxmin2
  integer plotymin2
  logical point_inside_box_2d
  integer pxcen
  integer pycen
  integer state
  integer unit
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
  call ps_setting_int ( 'GET', 'MARKER_SIZE', marker_size )
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
!    Henry McGilton, Mary Campione,
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
  integer marker_size
  integer plotxmin2
  integer plotymin2
  logical point_inside_box_2d
  integer pxcen
  integer pycen
  integer state
  integer unit
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

  write ( unit, '(2i6,a)' ) pxcen, pycen, ' moveto'
  write ( unit, '(2i6,a)' ) pxcen+1, pycen, ' lineto'
  write ( unit, '(a)' ) 'stroke'

  call ps_setting_real ( 'SET', 'XCUR', x )
  call ps_setting_real ( 'SET', 'YCUR', y )

  return
end
subroutine ps_marker_size ( marker_size )

!*****************************************************************************80
!
!! PS_MARKER_SIZE sets the marker size.
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
!    Input, integer MARKER_SIZE, the marker size.
!    0 is invisible, 1 is a single point.
!    A typical value is 3, 5 or 8.
!
  implicit none

  integer marker_size
  integer state
!
!  Determine if the PostScript state is acceptable.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state /= 2 .and. state /= 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_MARKER_SIZE - Fatal error!'
    write ( *, '(a,i9)' ) '  PostScript state is ', state
    write ( *, '(a)' ) '  PostScript state 2 or 3 is required.'
    return
  end if

  call ps_setting_int ( 'SET', 'MARKER_SIZE', marker_size )

  return
end
subroutine ps_moveto ( x, y )

!*****************************************************************************80
!
!! PS_MOVETO "moves to" a new point, which becomes the current point.
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
!    Henry McGilton, Mary Campione,
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
  integer plotxmin2
  integer plotymin2
  integer px
  integer py
  integer state
  integer unit
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
  write ( unit, '(2i6,a)' ) px, py, ' moveto'

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
  integer num_pages
  integer state
  real ( kind = 8 ) line_blue
  real ( kind = 8 ) line_green
  real ( kind = 8 ) line_red
  integer margin
  integer pagexmax
  integer pagexmin
  integer pageymax
  integer pageymin
  integer plotxmax
  integer plotxmin
  integer plotxmin2
  integer plotymax
  integer plotymin
  integer plotymin2
  integer unit
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

  integer state
  integer unit
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
subroutine ps_polygon_fill ( npoint, x, y )

!*****************************************************************************80
!
!! PS_POLYGON_FILL adds a filled polygon to a PostScript file.
!
!  Discussion:
!
!    A closed polygonal path is the sequence of line segments defined
!    by joining consecutive elements of a list of points; the path is
!    closed because the last point is joined to the first.  A filled
!    polygon is the area "inside" a closed polygonal path.  The meaning of
!    the word "inside" can be ambiguous in some cases.
!
!    The polygon fill color should be set before calling this routine.
!
!    The current point is not affected by this call.
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
!    Input, integer NPOINT, the number of points in the line.
!
!    Input, real ( kind = 8 ) X(NPOINT), Y(NPOINT), the X and Y components
!    of the points.
!
  implicit none

  integer npoint

  real ( kind = 8 ) alpha
  integer i
  integer plotxmin2
  integer plotymin2
  integer px
  integer py
  integer state
  integer unit
  real ( kind = 8 ) x(npoint)
  real ( kind = 8 ) xmin
  real ( kind = 8 ) y(npoint)
  real ( kind = 8 ) ymin
!
!  Refuse to handle fewer than 2 points.
!
  if ( npoint < 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_POLYGON_FILL - Warning!'
    write ( *, '(a)' ) '  Polygon has too few sides.'
    write ( *, '(a,i9)' ) '  NPOINT = ', npoint
    return
  end if
!
!  Determine if the PostScript state is acceptable.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state /= 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_POLYGON_FILL - Fatal error!'
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
  call ps_comment ( 'Draw a polygon' )

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
!  Fill the polygon.
!
  write ( unit, '(a)' ) 'fill'

  return
end
subroutine ps_landscape ( )

!*****************************************************************************80
!
!! PS_LANDSCAPE rotates the page from portrait to landscape.
!
!  Discussion:
!
!    PS_LANDSCAPE must be called AFTER a page has been set up.
!
!  Modified:
!
!    04 September 2001
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

  real ( kind = 8 ), parameter :: angle = 90.0D+00
  integer state
  integer unit
!
!  Determine if the PostScript state is acceptable.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state /= 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_LANDSCAPE - Fatal error!'
    write ( *, '(a,i9)' ) '  PostScript state is ', state
    write ( *, '(a)' ) '  PostScript state 3 is required.'
    return
  end if

  call ps_setting_int ( 'GET', 'UNIT', unit )

  write ( unit, '(f8.4,a)' ) angle, ' rotate'

  write ( unit, '(i4,i6,a)' ) 0, -792, ' translate'

  return
end
subroutine ps_rotate ( angle )

!*****************************************************************************80
!
!! PS_ROTATE rotates the coordinate system by a given angle.
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
!    Henry McGilton, Mary Campione,
!    PostScript by Example,
!    Addison-Wesley,
!    ISBN: 0-201-63228-4
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ANGLE, the angle of rotation, in degrees.
!
  implicit none

  real ( kind = 8 ) angle
  integer state
  integer unit
!
!  Determine if the PostScript state is acceptable.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state /= 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_ROTATE - Fatal error!'
    write ( *, '(a,i9)' ) '  PostScript state is ', state
    write ( *, '(a)' ) '  PostScript state 3 is required.'
    return
  end if

  call ps_setting_int ( 'GET', 'UNIT', unit )

  write ( unit, '(f8.4,a)' ) angle, ' rotate'

  return
end
subroutine ps_sector ( x0, y0, r, angle_min, angle_max )

!*****************************************************************************80
!
!! PS_SECTOR draws a circular sector.
!
!  Discussion:
!
!    As a side effect, the current point is set to the center of the circle.
!
!  Modified:
!
!    17 February 2001
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
!    Input, real ( kind = 8 ) X0, Y0, the coordinates of the center
!    of the circle.
!
!    Input, real ( kind = 8 ) R, the radius of the circle.
!
!    Input, real ( kind = 8 ) ANGLE_MIN, ANGLE_MAX, the minimum and 
!    maximum angles that define the sector, in degrees.
!
  implicit none

  real ( kind = 8 ) alpha
  real ( kind = 8 ) angle_max
  real ( kind = 8 ) angle_min
  integer plotxmin2
  integer plotymin2
  integer pr
  integer pxcen
  integer pycen
  real ( kind = 8 ) r
  integer state
  integer unit
  real ( kind = 8 ) x0
  real ( kind = 8 ) xmin
  real ( kind = 8 ) y0
  real ( kind = 8 ) ymin
!
!  Determine if the PostScript state is acceptable.
!
  call ps_setting_int ( 'GET', 'STATE', state )

  if ( state /= 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_SECTOR - Fatal error!'
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

  pxcen = plotxmin2 + nint ( alpha * ( x0 - xmin ) )
  pycen = plotymin2 + nint ( alpha * ( y0 - ymin ) )
  pr = nint ( alpha * r )

  write ( unit, '(a)' ) 'newpath'
  write ( unit, '(5i6,a)' ) pxcen, pycen, pr, int ( angle_min ), &
    int ( angle_max ), ' arc'
!
!  Draw the circle.
!
  write ( unit, '(a)' ) 'closepath stroke'

  call ps_setting_real ( 'SET', 'XCUR', x0 )
  call ps_setting_real ( 'SET', 'YCUR', y0 )

  return
end
subroutine ps_sector_fill ( x0, y0, r, angle_min, angle_max )

!*****************************************************************************80
!
!! PS_SECTOR_FILL draws a filled circular sector.
!
!  Modified:
!
!    17 February 2001
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
!    Input, real ( kind = 8 ) X0, Y0, the coordinates of the center of the disk.
!
!    Input, real ( kind = 8 ) R, the radius of the disk.
!
!    Input, real ( kind = 8 ) ANGLE_MIN, ANGLE_MAX, the minimum and
!    maximum angles that define the sector, in degrees.
!
  implicit none

  integer, parameter :: n = 33

  real ( kind = 8 ) angle_max
  real ( kind = 8 ) angle_min
  real ( kind = 8 ) r
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x0
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) y0

  call sector_points ( x0, y0, r, angle_min, angle_max, n, x, y )

  call ps_comment ( 'Draw a filled polygon.' )

  call ps_polygon_fill ( n, x, y )

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
!    Input/output, integer VALUE.
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
  integer, save :: line_width = 1
  integer, save :: marker_size = 0
  integer, save :: num_pages = 0
  integer, save :: pxmin = 0
  integer, save :: pymin = 0
  integer, save :: state = 0
  integer, save :: unit = 0
  integer value
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
subroutine ps_setting_print ( )

!*****************************************************************************80
!
!! PS_SETTING_PRINT prints the internal PS_WRITE parameters.
!
!  Modified:
!
!    04 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    NONE
!
  implicit none

  integer i
  real ( kind = 8 ) r

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PS_SETTING_PRINT:'
  write ( *, '(a)' ) '  The current internal PS_WRITE setting values:'
  write ( *, '(a)' ) ' '

  call ps_setting_real ( 'PRINT', 'ALPHA', r )
  call ps_setting_real ( 'PRINT', 'FILL_BLUE', r )
  call ps_setting_real ( 'PRINT', 'FILL_GREEN', r )
  call ps_setting_real ( 'PRINT', 'FILL_RED', r )
  call ps_setting_real ( 'PRINT', 'FONT_SIZE', r )
  call ps_setting_real ( 'PRINT', 'LINE_BLUE', r )
  call ps_setting_real ( 'PRINT', 'LINE_GREEN', r )
  call ps_setting_real ( 'PRINT', 'LINE_RED', r )
  call ps_setting_int ( 'PRINT', 'LINE_WIDTH', i )
  call ps_setting_int ( 'PRINT', 'NUM_PAGES', i )
  call ps_setting_int ( 'PRINT', 'PXMIN', i )
  call ps_setting_int ( 'PRINT', 'PYMIN', i )
  call ps_setting_int ( 'PRINT', 'STATE', i )
  call ps_setting_int ( 'PRINT', 'UNIT', i )
  call ps_setting_real ( 'PRINT', 'XCUR', r )
  call ps_setting_real ( 'PRINT', 'XMIN', r )
  call ps_setting_real ( 'PRINT', 'XMAX', r )
  call ps_setting_real ( 'PRINT', 'YCUR', r )
  call ps_setting_real ( 'PRINT', 'YMIN', r )
  call ps_setting_real ( 'PRINT', 'YMAX', r )

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
subroutine ps_square ( x0, y0, r )

!*****************************************************************************80
!
!! PS_SQUARE draws a square.
!
!  Modified:
!
!    07 October 2000
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
!    Input, real ( kind = 8 ) X0, Y0, the coordinates of the center 
!    of the square.
!
!    Input, real ( kind = 8 ) R, the radius of the square.
!
  implicit none

  integer, parameter :: n = 4

  real ( kind = 8 ) r
  real ( kind = 8 ) xvec(n)
  real ( kind = 8 ) x0
  real ( kind = 8 ) yvec(n)
  real ( kind = 8 ) y0

  xvec(1:n) = (/ x0-r, x0+r, x0+r, x0-r /)
  yvec(1:n) = (/ y0-r, y0-r, y0+r, y0+r /)

  call ps_line_closed ( n, xvec, yvec )

  return
end
subroutine ps_square_fill ( x0, y0, r )

!*****************************************************************************80
!
!! PS_SQUARE_FILL draws a filled square.
!
!  Modified:
!
!    07 October 2000
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
!    Input, real ( kind = 8 ) X0, Y0, the coordinates of the center of the box.
!
!    Input, real ( kind = 8 ) R, the radius of the box.
!
  implicit none

  real ( kind = 8 ) r
  real ( kind = 8 ) x(4)
  real ( kind = 8 ) x0
  real ( kind = 8 ) y(4)
  real ( kind = 8 ) y0

  x(1:4) = (/ x0-r, x0+r, x0+r, x0-r /)
  y(1:4) = (/ y0-r, y0-r, y0+r, y0+r /)

  call ps_polygon_fill ( 4, x, y )

  return
end
subroutine ps_star ( x0, y0, r )

!*****************************************************************************80
!
!! PS_STAR draws an open star of given radius.
!
!  Discussion:
!
!    The radius refers to the circumscribing circle.
!
!  Modified:
!
!    13 October 2000
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
!    Input, real ( kind = 8 ) X0, Y0, the coordinates of the center of the star.
!
!    Input, real ( kind = 8 ) R, the radius of the star.
!
  implicit none

  integer, parameter :: n = 5

  real ( kind = 8 ) angle
  integer i
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r
  real ( kind = 8 ) xvec(n)
  real ( kind = 8 ) x0
  real ( kind = 8 ) yvec(n)
  real ( kind = 8 ) y0

  do i = 1, n
    angle = pi * real ( mod ( ( - 7 + i * 8 ), 20 ), kind = 8 ) / 10.0D+00
    xvec(i) = x0 + r * cos ( angle )
    yvec(i) = y0 + r * sin ( angle )
  end do

  call ps_line_closed ( n, xvec, yvec )

  return
end
subroutine ps_star_circle ( x0, y0, r )

!*****************************************************************************80
!
!! PS_STAR_CIRCLE draws open circles at the 10 points on a star of given radius.
!
!  Discussion:
!
!    The radius refers to the circumscribing circle.
!
!  Modified:
!
!    14 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X0, Y0, the coordinates of the center of the star.
!
!    Input, real ( kind = 8 ) R, the radius of the star.
!
  implicit none

  integer, parameter :: n = 10

  integer i
  real ( kind = 8 ) r
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x0
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) y0

  call ps_star_points ( x0, y0, r, x, y )

  do i = 1, n
    call ps_mark_circle ( x(i), y(i) )
  end do

  return
end
subroutine ps_star_disk ( x0, y0, r )

!*****************************************************************************80
!
!! PS_STAR_DISK draws filled disks at the 10 points on a star of given radius.
!
!  Discussion:
!
!    The radius refers to the circumscribing circle.
!
!  Modified:
!
!    14 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X0, Y0, the coordinates of the center of the star.
!
!    Input, real ( kind = 8 ) R, the radius of the star.
!
  implicit none

  integer, parameter :: n = 10

  integer i
  real ( kind = 8 ) r
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x0
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) y0

  call ps_star_points ( x0, y0, r, x, y )

  do i = 1, n
    call ps_mark_disk ( x(i), y(i) )
  end do

  return
end
subroutine ps_star_points ( x0, y0, r, x, y )

!*****************************************************************************80
!
!! PS_STAR_POINTS returns 10 points on a star of given radius.
!
!  Discussion:
!
!    The radius refers to the circumscribing circle.
!
!  Modified:
!
!    14 October 2000
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
!    Input, real ( kind = 8 ) X0, Y0, the coordinates of the center of the star.
!
!    Input, real ( kind = 8 ) R, the radius of the star.
!
!    Output, real ( kind = 8 ) X(10), Y(10), the coordinates of points 
!    on the star.
!
  implicit none

  integer, parameter :: n = 10

  real ( kind = 8 ) angle
  integer i
  integer j
  real ( kind = 8 ) phi
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x0
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) y0

  phi = ( 1.0D+00 + sqrt ( 5.0D+00 ) ) / 2.0D+00

  j = 0
  do i = 1, 5
    j = j + 1
    angle = pi * real ( mod ( ( - 7 + i * 8 ), 20 ), kind = 8 ) / 10.0D+00
    x(j) = x0 + r * cos ( angle )
    y(j) = y0 + r * sin ( angle )
    j = j + 1
    angle = pi * real ( mod ( ( - 5 + i * 8 ), 20 ), kind = 8 ) / 10.0D+00
    x(j) = x0 + r * cos ( angle ) / phi**2
    y(j) = y0 + r * sin ( angle ) / phi**2
  end do

  return
end
subroutine ps_triangle ( x1, y1, x2, y2, x3, y3 )

!*****************************************************************************80
!
!! PS_TRIANGLE draws an open triangle.
!
!  Modified:
!
!    17 February 2001
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
!    Input, real ( kind = 8 ) X1, Y1, X2, Y2, X3, Y3, the coordinates
!    of the triangle.
!
  implicit none

  integer, parameter :: n = 3

  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ) xvec(n)
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) y3
  real ( kind = 8 ) yvec(n)

  xvec(1:n) = (/ x1, x2, x3 /)
  yvec(1:n) = (/ y1, y2, y3 /)

  call ps_line_closed ( n, xvec, yvec )

  return
end
subroutine ps_triangle_fill ( x1, y1, x2, y2, x3, y3 )

!*****************************************************************************80
!
!! PS_TRIANGLE_FILL draws a filled triangle.
!
!  Modified:
!
!    17 February 2001
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
!    Input, real ( kind = 8 ) X1, Y1, X2, Y2, X3, Y3, the coordinates
!    of the triangle.
!
  implicit none

  integer, parameter :: n = 3

  real ( kind = 8 ) xvec(n)
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ) yvec(n)
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) y3

  xvec(1:n) = (/ x1, x2, x3 /)
  yvec(1:n) = (/ y1, y2, y3 /)

  call ps_polygon_fill ( n, xvec, yvec )

  return
end
function r8_modp ( x, y )

!*****************************************************************************80
!
!! R8_MODP returns the nonnegative remainder of real division.
!
!  Formula:
!
!    If
!      REM = R8_MODP ( X, Y )
!      RMULT = ( X - REM ) / Y
!    then
!      X = Y * RMULT + REM
!    where REM is always nonnegative.
!
!  Comments:
!
!    The MOD function computes a result with the same sign as the
!    quantity being divided.  Thus, suppose you had an angle A,
!    and you wanted to ensure that it was between 0 and 360.
!    Then mod(A,360.0) would do, if A was positive, but if A
!    was negative, your result would be between -360 and 0.
!
!    On the other hand, R8_MODP(A,360.0) is between 0 and 360, always.
!
!  Examples:
!
!        I         J     MOD R8_MODP  R8_MODP Factorization
!
!      107        50       7       7    107 =  2 *  50 + 7
!      107       -50       7       7    107 = -2 * -50 + 7
!     -107        50      -7      43   -107 = -3 *  50 + 43
!     -107       -50      -7      43   -107 =  3 * -50 + 43
!
!  Modified:
!
!    24 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the number to be divided.
!
!    Input, real ( kind = 8 ) Y, the number that divides X.
!
!    Output, real ( kind = 8 ) R8_MODP, the nonnegative remainder 
!    when X is divided by Y.
!
  implicit none

  real ( kind = 8 ) r8_modp
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  if ( y == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_MODP - Fatal error!'
    write ( *, '(a,g14.6)' ) '  R8_MODP ( X, Y ) called with Y = ', y
    stop
  end if

  r8_modp = mod ( x, y )

  if ( r8_modp < 0.0D+00 ) then
    r8_modp = r8_modp + abs ( y )
  end if

  return
end
subroutine r8_swap ( x, y )

!*****************************************************************************80
!
!! R8_SWAP swaps two real values.
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
!  Examples:
!
!    S_EQI ( 'Anjana', 'ANJANA' ) is .TRUE.
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
  integer i
  integer len1
  integer len2
  integer lenc
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
subroutine sector_points ( x0, y0, r, angle_min, angle_max, n, x, y )

!*****************************************************************************80
!
!! SECTOR_POINTS returns N equally spaced points on a circle in 2D.
!
!  Discussion:
!
!    The first point is always ( X0, Y0 ), and the N-1 subsequent points
!    proceed counterclockwise around the circle.
!
!  Modified:
!
!    17 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X0, Y0, the coordinates of the center
!    of the circle.
!
!    Input, real ( kind = 8 ) R, the radius of the circle.
!
!    Input, real ( kind = 8 ) ANGLE_MIN, ANGLE_MAX, the minimum and 
!    maximum angles that define the sector, in degrees.
!
!    Input, integer N, the number of points desired.  N must be at least 3.
!
!    Output, real ( kind = 8 ) X(N), Y(N), the coordinates of points 
!    on the circle.
!
  implicit none

  integer n

  real ( kind = 8 ) angle
  real ( kind = 8 ) angle_max
  real ( kind = 8 ) angle_min
  real ( kind = 8 ) degrees_to_radians
  integer i
  real ( kind = 8 ) r
  real ( kind = 8 ) x0
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y0
  real ( kind = 8 ) y(n)

  x(1) = x0
  y(1) = y0

  do i = 2, n
    angle = ( real ( n - i,     kind = 8 ) * angle_min   &
            + real (     i - 2, kind = 8 ) * angle_max ) &
            / real ( n     - 2, kind = 8  )
    angle = degrees_to_radians ( angle )
    x(i) = x0 + r * cos ( angle )
    y(i) = y0 + r * sin ( angle )
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
  integer d
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer n
  integer s
  integer values(8)
  integer y

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
