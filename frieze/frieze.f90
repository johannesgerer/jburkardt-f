program main

!*****************************************************************************80
!
!! MAIN is the main program for FRIEZE.
!
!  Discussion:
!
!    FRIEZE demonstrates the use of blending to tile a region with a pattern.
!
!    We used this program to see what blending would look like.  We were
!    really interested in 3D models, but wanted to see what was involved
!    with a simple 2D problem first.
!
!    This program requires the PS_WRITE routines in order to create a
!    PostScript file with an image of the tiled region.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 September 1999
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

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FRIEZE'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Demonstrate the use of blending to "tile"'
  write ( *, '(a)' ) '  a region with a pattern.'
!
!  Draw the pattern.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FRIEZE'
  write ( *, '(a)' ) '  Draw the wallpaper pattern:'

  call draw_pattern ( 'pattern.ps' )
!
!  Draw the region
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FRIEZE'
  write ( *, '(a)' ) '  Draw the region:'

  call draw_region ( 'region.ps' )
!
!  Draw the region subdivided into cells
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FRIEZE'
  write ( *, '(a)' ) '  Draw the region divided into cells:'

  call draw_cells ( 'cells.ps' )
!
!  Draw the frieze.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FRIEZE'
  write ( *, '(a)' ) '  Draw the region, divided into cells,'
  write ( *, '(a)' ) '  with each cell covered by a copy of the'
  write ( *, '(a)' ) '  wallpaper pattern.'

  call draw_frieze ( 'frieze.ps' )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FRIEZE'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine blend1d2d ( r, s, x, y )

!*****************************************************************************80
!
!! BLEND1D2D uses transfinite interpolation on a 2D cell.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 July 1999
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
!    Input, real ( kind = 8 ) R, S, the (R,S) coordinates of a point.
!
!    Output, real ( kind = 8 ) X, Y, the (X,Y) coordinates of the point.
!
  implicit none

  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) x
  real ( kind = 8 ) xe
  real ( kind = 8 ) xn
  real ( kind = 8 ) xne
  real ( kind = 8 ) xnw
  real ( kind = 8 ) xs
  real ( kind = 8 ) xse
  real ( kind = 8 ) xsw
  real ( kind = 8 ) xw
  real ( kind = 8 ) y
  real ( kind = 8 ) ye
  real ( kind = 8 ) yn
  real ( kind = 8 ) yne
  real ( kind = 8 ) ynw
  real ( kind = 8 ) ys
  real ( kind = 8 ) yse
  real ( kind = 8 ) ysw
  real ( kind = 8 ) yw
!
!  Find the (X,Y) coordinates of the corners.
!
  call boundary_2d ( 1.0D+00, 0.0D+00, xse, yse )
  call boundary_2d ( 1.0D+00, 1.0D+00, xne, yne )
  call boundary_2d ( 0.0D+00, 1.0D+00, xnw, ynw )
  call boundary_2d ( 0.0D+00, 0.0D+00, xsw, ysw )
!
!  Find the (X,Y) coordinates of corresponding points on the sides.
!
  call boundary_2d ( r, 1.0D+00, xn, yn )
  call boundary_2d ( r, 0.0D+00, xs, ys )
  call boundary_2d ( 0.0D+00, s, xw, yw )
  call boundary_2d ( 1.0D+00, s, xe, ye )
!
!  Now interpolate the (X,Y) coordinates of the point in the interior.
!
  x =                          ( 1.0D+00 - s ) * xs &
           -             r   * ( 1.0D+00 - s ) * xse &
           +             r                     * xe &
           -             r   *             s   * xne &
           +                               s   * xn &
           - ( 1.0D+00 - r ) *             s   * xnw &
           + ( 1.0D+00 - r )                   * xw &
           - ( 1.0D+00 - r ) * ( 1.0D+00 - s ) * xsw

  y =                          ( 1.0D+00 - s ) * ys &
           -             r   * ( 1.0D+00 - s ) * yse &
           +             r                     * ye &
           -             r   *             s   * yne &
           +                               s   * yn &
           - ( 1.0D+00 - r ) *             s   * ynw &
           + ( 1.0D+00 - r )                   * yw &
           - ( 1.0D+00 - r ) * ( 1.0D+00 - s ) * ysw

  return
end
subroutine boundary_2d ( r, s, x, y )

!*****************************************************************************80
!
!! BOUNDARY_2D returns (X,Y) points on a side of the boundary.
!
!  Discussion:
!
!    The boundary is divided into four segments:
!
!      BOTTOM: ( 3 * cos ( ( 3-2*r)*pi/4 ),   3 * sin ( ( 3-2*r)*pi/4 ) );
!      RIGHT:  ( (3+s) * sqrt(2)/2,          (3+s)*sqrt(2)/2            );
!      TOP:    ( 4 * cos ( ( 2*r+1)*pi/4 ),   4 * sin ( (2*r+1)*pi/4 )  );
!      LEFT:   ( -(4-s)*sqrt(2)/2,           (4-s)*sqrt(2)/2            );
!
!    A
!    |
!    1  *-----------*
!    |  |           |
!    S  |           |
!    |  |           |
!    0  *-----------*
!    |
!    +--0-----R-----1--->
!
!    I'm assuming that R and S both go from 0 to 1 exactly.  Since writing
!    this code, I've come to prefer to allow the more general case where
!    the ranges of R and S are allowed to be other values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 December 1998
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
!    Input, real ( kind = 8 ) R, S, the (R,S) coordinates of a point 
!    on the boundary.
!
!    Output, real ( kind = 8 ) X, Y, the (X,Y) coordinates of the point.
!
  implicit none

  real ( kind = 8 ) angle
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r
  real ( kind = 8 ) radius
  real ( kind = 8 ) s
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  if ( s == 0.0D+00 ) then
    angle = ( 3.0D+00 - 2.0D+00 * r ) * pi / 4.0D+00
    x = 3.0D+00 * cos ( angle )
    y = 3.0D+00 * sin ( angle )
  else if ( r == 1.0D+00 ) then
    x = ( 3.0D+00 + s ) * sqrt ( 2.0D+00 ) / 2.0D+00
    y = ( 3.0D+00 + s ) * sqrt ( 2.0D+00 ) / 2.0D+00
  else if ( s == 1.0D+00 ) then
    angle = ( 3.0D+00 - 2.0D+00 * r ) * pi / 4.0D+00
    radius = 4.0D+00 + 0.2D+00 * sin ( 4.0D+00 * pi * r )
    x = radius * cos ( angle )
    y = radius * sin ( angle )
  else if ( r == 0.0D+00 ) then
    x = - ( 3.0D+00 + s ) * sqrt ( 2.0D+00 ) / 2.0D+00
    y =   ( 3.0D+00 + s ) * sqrt ( 2.0D+00 ) / 2.0D+00
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BOUNDARY:'
    write ( *, '(a)' ) '  Illegal side coordinates!'
    write ( *, '(a,2g14.6)' ) '  (R,S) = ', r, s
    stop
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
subroutine draw_cell ( re, rw, sn, ss )

!*****************************************************************************80
!
!! DRAW_CELL draws the borders of a given cell.
!
!  Discussion:
!
!    The internal parameter NPOINT determines how many points are
!    used to draw each of the four lines that define the border.
!    Increase it to draw a wiggly boundary more accurately.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 July 1999
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
!    Input, real ( kind = 8 ) RE, RW, SN, SS, the extreme (east and west)
!    R coordinates and extreme (north and south) S coordinates of the subcell.
!
  implicit none

  integer ( kind = 4 ), parameter :: npoint = 5

  integer ( kind = 4 ) i
  real ( kind = 8 ) r
  real ( kind = 8 ) re
  real ( kind = 8 ) rw
  real ( kind = 8 ) s
  real ( kind = 8 ) sn
  real ( kind = 8 ) ss
  real ( kind = 8 ) xpoint(npoint)
  real ( kind = 8 ) ypoint(npoint)
!
!  East cell border.
!
  r = re
  do i = 1, npoint
    s = ( real ( npoint - i,     kind = 8 ) * ss   &
        + real (          i - 1, kind = 8 ) * sn ) &
        / real ( npoint     - 1, kind = 8 )

    call blend1d2d ( r, s, xpoint(i), ypoint(i) )
  end do

  call ps_line_open ( npoint, xpoint, ypoint )
!
!  North cell border.
!
  s = sn
  do i = 1, npoint

    r = ( real ( npoint - i,     kind = 8 ) * re   &
        + real (          i - 1, kind = 8 ) * rw ) &
        / real ( npoint     - 1, kind = 8 )

    call blend1d2d ( r, s, xpoint(i), ypoint(i) )
  end do

  call ps_line_open ( npoint, xpoint, ypoint )
!
!  West cell border.
!
  r = rw
  do i = 1, npoint
    s = ( real ( npoint - i,     kind = 8 ) * sn   &
        + real (          i - 1, kind = 8 ) * ss ) &
        / real ( npoint     - 1, kind = 8 )
    call blend1d2d ( r, s, xpoint(i), ypoint(i) )
  end do

  call ps_line_open ( npoint, xpoint, ypoint )
!
!  South cell border.
!
  s = ss
  do i = 1, npoint
    r = ( real ( npoint - i,     kind = 8 ) * rw   &
        + real (          i - 1, kind = 8 ) * re ) &
        / real ( npoint     - 1, kind = 8 )
    call blend1d2d ( r, s, xpoint(i), ypoint(i) )
  end do

  call ps_line_open ( npoint, xpoint, ypoint )

  return
end
subroutine draw_cells ( filename )

!*****************************************************************************80
!
!! DRAW_CELLS makes an image of the cells that make up the region.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 September 1999
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4
  integer ( kind = 4 ), parameter :: nfine = 21

  character ( len = * ) filename
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icol
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) irow
  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nrow
  real ( kind = 8 ) r
  real ( kind = 8 ) re
  real ( kind = 8 ) rw
  real ( kind = 8 ) s
  real ( kind = 8 ) sn
  real ( kind = 8 ) ss
  real ( kind = 8 ) x(nfine)
  real ( kind = 8 ) xcval(4)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) xne
  real ( kind = 8 ) xnw
  real ( kind = 8 ) xse
  real ( kind = 8 ) xsw
  real ( kind = 8 ) y(nfine)
  real ( kind = 8 ) ycval(n)
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin
  real ( kind = 8 ) yne
  real ( kind = 8 ) ynw
  real ( kind = 8 ) yse
  real ( kind = 8 ) ysw
!
!  Open the PostScript file.
!
  iunit = 1

  call ps_file_open ( filename, iunit, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DRAW_CELLS - Fatal error!'
    write ( *, '(a,i8)' ) '  PostScript file creation error ', ierror
    stop
  end if

  xmin = -3.0D+00
  ymin =  2.0D+00
  xmax = +3.0D+00
  ymax =  4.5D+00

  call ps_file_head ( filename )

  call ps_page_head ( xmin, ymin, xmax, ymax )
!
!  Query the BOUNDARY routine for the (X,Y) locations of the corners.
!
  call boundary_2d ( 1.0D+00, 0.0D+00, xse, yse )
  call boundary_2d ( 1.0D+00, 1.0D+00, xne, yne )
  call boundary_2d ( 0.0D+00, 1.0D+00, xnw, ynw )
  call boundary_2d ( 0.0D+00, 0.0D+00, xsw, ysw )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '(X,Y) corners of the total region:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,2g14.6)' ) 'XSE, YSE = ', xse, yse
  write ( *, '(a,2g14.6)' ) 'XNE, YNE = ', xne, yne
  write ( *, '(a,2g14.6)' ) 'XNW, YNW = ', xnw, ynw
  write ( *, '(a,2g14.6)' ) 'XSW, YSW = ', xsw, ysw
!
!  Draw the boundary of the region.
!
  call ps_color_line_set ( 0.7D+00, 0.7D+00, 0.7D+00 )

  r = 1.0D+00
  do i = 1, nfine
    s = real ( i - 1, kind = 8 ) / real ( nfine - 1, kind = 8 )
    call boundary_2d ( r, s, x(i), y(i) )
  end do

  call ps_line_open ( nfine, x, y )

  s = 1.0D+00
  do i = nfine, 1, -1
    r = real ( i - 1, kind = 8 ) / real ( nfine - 1, kind = 8 )
    call boundary_2d ( r, s, x(i), y(i) )
  end do

  call ps_line_open ( nfine, x, y )

  r = 0.0D+00
  do i = nfine, 1, -1
    s = real ( i - 1, kind = 8 ) / real ( nfine - 1, kind = 8 )
    call boundary_2d ( r, s, x(i), y(i) )
  end do

  call ps_line_open ( nfine, x, y )

  s = 0.0D+00
  do i = 1, nfine
    r = real ( i - 1, kind = 8 ) / real ( nfine - 1, kind = 8 )
    call boundary_2d ( r, s, x(i), y(i) )
  end do

  call ps_line_open ( nfine, x, y )
!
!  Map the reference element into subregion IROW, JCOL.
!
  ncol = 6
  nrow = 2

  do irow = 1, nrow

    sn = real ( irow, kind = 8 ) / real ( nrow, kind = 8 )
    ss = real ( irow - 1, kind = 8 ) / real ( nrow, kind = 8 )

    do icol = 1, ncol

      re = real ( icol, kind = 8 ) / real ( ncol, kind = 8 )
      rw = real ( icol - 1, kind = 8 ) / real ( ncol, kind = 8 )

      call blend1d2d ( re, ss, xcval(1), ycval(1) )
      call blend1d2d ( re, sn, xcval(2), ycval(2) )
      call blend1d2d ( rw, sn, xcval(3), ycval(3) )
      call blend1d2d ( rw, ss, xcval(4), ycval(4) )

      do i = 1, 4
        write ( *, * ) xcval(i), ycval(i)
      end do

      call ps_color_line_set ( 0.0D+00, 0.0D+00, 0.4D+00 )

      call draw_cell ( re, rw, sn, ss )

    end do
  end do
!
!  Finish up the PostScript file.
!
  call ps_page_tail

  call ps_file_tail

  call ps_file_close ( iunit )

  return
end
subroutine draw_diagonals ( )

!*****************************************************************************80
!
!! DRAW_DIAGONALS draws the diagonals of the entire region.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 July 1999
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

  integer ( kind = 4 ), parameter :: nfine = 21

  integer ( kind = 4 ) i
  real ( kind = 8 ) r
  real ( kind = 8 ) re
  real ( kind = 8 ) rw
  real ( kind = 8 ) s
  real ( kind = 8 ) sn
  real ( kind = 8 ) ss
  real ( kind = 8 ) x(nfine)
  real ( kind = 8 ) y(nfine)

  re = 1.0D+00
  rw = 0.0D+00
  sn = 1.0D+00
  ss = 0.0D+00

  do i = 1, nfine
    r = ( real ( nfine - i,     kind = 8 ) * rw   &
        + real (         i - 1, kind = 8 ) * re ) &
        / real ( nfine     - 1, kind = 8 )

    s = ( real ( nfine - i,     kind = 8 ) * sn   &
        + real (         i - 1, kind = 8 ) * ss ) &
        / real ( nfine     - 1, kind = 8 )

    call blend1d2d ( r, s, x(i), y(i) )
  end do

  call ps_line_open ( nfine, x, y )

  do i = 1, nfine

    r = ( real ( nfine - i,     kind = 8 ) * rw   &
        + real (         i - 1, kind = 8 ) * re ) &
        / real ( nfine     - 1, kind = 8 )

    s = ( real ( nfine - i,     kind = 8 ) * ss   &
        + real (         i - 1, kind = 8 ) * sn ) &
        / real ( nfine     - 1, kind = 8 )

    call blend1d2d ( r, s, x(i), y(i) )
  end do

  call ps_line_open ( nfine, x, y )

  return
end
subroutine draw_frieze ( filename )

!*****************************************************************************80
!
!! DRAW_FRIEZE draws the frieze.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 September 1999
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: max_end = 100
  integer ( kind = 4 ), parameter :: max_pat = 100
  integer ( kind = 4 ), parameter :: n = 4
  integer ( kind = 4 ), parameter :: nfine = 21

  character ( len = * ) filename
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icol
  integer ( kind = 4 ) iend(max_end)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) irow
  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) num_end
  integer ( kind = 4 ) num_pat
  integer ( kind = 4 ) nrow
  real ( kind = 8 ) r
  real ( kind = 8 ) re
  real ( kind = 8 ) rw
  real ( kind = 8 ) s
  real ( kind = 8 ) sn
  real ( kind = 8 ) ss
  real ( kind = 8 ) x(nfine)
  real ( kind = 8 ) xcval(4)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) xne
  real ( kind = 8 ) xnw
  real ( kind = 8 ) xpat(max_pat)
  real ( kind = 8 ) xse
  real ( kind = 8 ) xsw
  real ( kind = 8 ) y(nfine)
  real ( kind = 8 ) ycval(n)
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin
  real ( kind = 8 ) yne
  real ( kind = 8 ) ynw
  real ( kind = 8 ) ypat(max_pat)
  real ( kind = 8 ) yse
  real ( kind = 8 ) ysw
!
!  Get the pattern to be used as a tile.
!
  call pattern ( iend, max_end, num_end, max_pat, num_pat, xpat, ypat )
!
!  Open the PostScript file.
!
  iunit = 1

  call ps_file_open ( filename, iunit, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DRAW_FRIEZE'
    write ( *, '(a,i8)' ) '  PostScript file creation error ', ierror
    stop
  end if

  xmin = -3.0D+00
  ymin =  2.0D+00
  xmax = +3.0D+00
  ymax =  4.5D+00

  call ps_file_head ( filename )

  call ps_page_head ( xmin, ymin, xmax, ymax )
!
!  Query the BOUNDARY routine for the (X,Y) locations of the corners.
!
  call boundary_2d ( 1.0D+00, 0.0D+00, xse, yse )
  call boundary_2d ( 1.0D+00, 1.0D+00, xne, yne )
  call boundary_2d ( 0.0D+00, 1.0D+00, xnw, ynw )
  call boundary_2d ( 0.0D+00, 0.0D+00, xsw, ysw )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '(X,Y) corners of the total region:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,2g14.6)' ) 'XSE, YSE = ', xse, yse
  write ( *, '(a,2g14.6)' ) 'XNE, YNE = ', xne, yne
  write ( *, '(a,2g14.6)' ) 'XNW, YNW = ', xnw, ynw
  write ( *, '(a,2g14.6)' ) 'XSW, YSW = ', xsw, ysw
!
!  Draw the boundary of the region.
!
  call ps_color_line_set ( 0.7D+00, 0.7D+00, 0.7D+00 )

  r = 1.0D+00
  do i = 1, nfine
    s = real ( i - 1, kind = 8 ) / real ( nfine - 1, kind = 8 )
    call boundary_2d ( r, s, x(i), y(i) )
  end do

  call ps_line_open ( nfine, x, y )

  s = 1.0D+00
  do i = nfine, 1, -1
    r = real ( i - 1, kind = 8 ) / real ( nfine - 1, kind = 8 )
    call boundary_2d ( r, s, x(i), y(i) )
  end do

  call ps_line_open ( nfine, x, y )

  r = 0.0D+00
  do i = nfine, 1, -1
    s = real ( i - 1, kind = 8 ) / real ( nfine - 1, kind = 8 )
    call boundary_2d ( r, s, x(i), y(i) )
  end do

  call ps_line_open ( nfine, x, y )

  s = 0.0D+00
  do i = 1, nfine
    r = real ( i - 1, kind = 8 ) / real ( nfine - 1, kind = 8 )
    call boundary_2d ( r, s, x(i), y(i) )
  end do

  call ps_line_open ( nfine, x, y )
!
!  Map the reference element into subregion IROW, JCOL.
!
  ncol = 6
  nrow = 2

  do irow = 1, nrow

    sn = real ( irow, kind = 8 ) / real ( nrow, kind = 8 )
    ss = real ( irow - 1, kind = 8 ) / real ( nrow, kind = 8 )

    do icol = 1, ncol

      re = real ( icol, kind = 8 ) / real ( ncol, kind = 8 )
      rw = real ( icol - 1, kind = 8 ) / real ( ncol, kind = 8 )

      call blend1d2d ( re, ss, xcval(1), ycval(1) )
      call blend1d2d ( re, sn, xcval(2), ycval(2) )
      call blend1d2d ( rw, sn, xcval(3), ycval(3) )
      call blend1d2d ( rw, ss, xcval(4), ycval(4) )

      do i = 1, 4
        write ( *, '(2g14.6)' ) xcval(i), ycval(i)
      end do

      call ps_color_line_set ( 0.0D+00, 0.0D+00, 0.4D+00 )

      call draw_cell ( re, rw, sn, ss )

      call ps_color_line_set ( 1.0D+00, 0.0D+00, 0.0D+00 )

      call draw_pat ( iend, num_end, num_pat, xpat, ypat, re, rw, sn, ss )

    end do
  end do
!
!  Finish up the PostScript file.
!
  call ps_page_tail

  call ps_file_tail

  call ps_file_close ( iunit )

  return
end
subroutine draw_pattern ( filename )

!*****************************************************************************80
!
!! DRAW_PATTERN draws the pattern that will be repeated in the frieze.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 September 1999
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: max_end = 100
  integer ( kind = 4 ), parameter :: max_pat = 100
  integer ( kind = 4 ), parameter :: MAX_POINT = 20
  integer ( kind = 4 ), parameter :: nfine = 21

  character ( len = * ) filename
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iend(max_end)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ), parameter :: iunit = 1
  integer ( kind = 4 ) ipat
  integer ( kind = 4 ) j
  integer ( kind = 4 ) npoint
  integer ( kind = 4 ) num_end
  integer ( kind = 4 ) num_pat
  real ( kind = 8 ) x(nfine)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) xpat(max_pat)
  real ( kind = 8 ) xpoint(MAX_POINT)
  real ( kind = 8 ) y(nfine)
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin
  real ( kind = 8 ) ypat(max_pat)
  real ( kind = 8 ) ypoint(MAX_POINT)
!
!  Get the pattern to be used as a tile.
!
  call pattern ( iend, max_end, num_end, max_pat, num_pat, xpat, ypat )
!
!  Open the PostScript file.
!
  call ps_file_open ( filename, iunit, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DRAW_PATTERN - Fatal error!'
    write ( *, '(a,i8)' ) '  PostScript file creation error ', ierror
    stop
  end if

  xmin =  0.0D+00
  ymin =  0.0D+00
  xmax =  1.0D+00
  ymax =  1.0D+00

  call ps_file_head ( filename )

  call ps_page_head ( xmin, ymin, xmax, ymax )
!
!  Draw the boundary of the region.
!
  call ps_color_line_set ( 0.7D+00, 0.7D+00, 0.7D+00 )

  x(1) = xmin
  y(1) = ymin
  x(2) = xmax
  y(2) = ymin
  x(3) = xmax
  y(3) = ymax
  x(4) = xmin
  y(4) = ymax
  x(5) = xmin
  y(5) = ymin

  call ps_line_open ( 5, x, y )

  call ps_color_line_set ( 1.0D+00, 0.0D+00, 0.0D+00 )

  ipat = 0
  do i = 1, num_end - 1

    npoint = 0
    do j = iend(i)+1, iend(i+1)
      npoint = npoint + 1
      ipat = ipat + 1
      xpoint(npoint) = xpat(ipat)
      ypoint(npoint) = ypat(ipat)
    end do

    call ps_line_open ( npoint, xpoint, ypoint )

  end do
!
!  Finish up the PostScript file.
!
  call ps_page_tail ( )

  call ps_file_tail ( )

  call ps_file_close ( iunit )

  return
end
subroutine draw_region ( filename )

!*****************************************************************************80
!
!! DRAW_REGION makes an image of the region.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 September 1999
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4
  integer ( kind = 4 ), parameter :: nfine = 21

  character ( len = * ) filename
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icol
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) irow
  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) nrow
  real ( kind = 8 ) r
  real ( kind = 8 ) re
  real ( kind = 8 ) rw
  real ( kind = 8 ) s
  real ( kind = 8 ) sn
  real ( kind = 8 ) ss
  real ( kind = 8 ) x(nfine)
  real ( kind = 8 ) xcval(4)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) xne
  real ( kind = 8 ) xnw
  real ( kind = 8 ) xse
  real ( kind = 8 ) xsw
  real ( kind = 8 ) y(nfine)
  real ( kind = 8 ) ycval(n)
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin
  real ( kind = 8 ) yne
  real ( kind = 8 ) ynw
  real ( kind = 8 ) yse
  real ( kind = 8 ) ysw
!
!  Open the PostScript file.
!
  iunit = 1

  call ps_file_open ( filename, iunit, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DRAW_REGION - Fatal error.'
    write ( *, '(a,i8)' ) '  PostScript file creation error ', ierror
    stop
  end if

  xmin = -3.0D+00
  ymin =  2.0D+00
  xmax = +3.0D+00
  ymax =  4.5D+00

  call ps_file_head ( filename )

  call ps_page_head ( xmin, ymin, xmax, ymax )
!
!  Query the BOUNDARY routine for the (X,Y) locations of the corners.
!
  call boundary_2d ( 1.0D+00, 0.0D+00, xse, yse )
  call boundary_2d ( 1.0D+00, 1.0D+00, xne, yne )
  call boundary_2d ( 0.0D+00, 1.0D+00, xnw, ynw )
  call boundary_2d ( 0.0D+00, 0.0D+00, xsw, ysw )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '(X,Y) corners of the total region:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,2g14.6)' ) 'XSE, YSE = ', xse, yse
  write ( *, '(a,2g14.6)' ) 'XNE, YNE = ', xne, yne
  write ( *, '(a,2g14.6)' ) 'XNW, YNW = ', xnw, ynw
  write ( *, '(a,2g14.6)' ) 'XSW, YSW = ', xsw, ysw
!
!  Draw the boundary of the region.
!
  call ps_color_line_set ( 0.1D+00, 0.5D+00, 1.0D+00 )

  r = 1.0D+00
  do i = 1, nfine
    s = real ( i - 1, kind = 8 ) / real ( nfine - 1, kind = 8 )
    call boundary_2d ( r, s, x(i), y(i) )
  end do

  call ps_line_open ( nfine, x, y )

  s = 1.0D+00
  do i = nfine, 1, -1
    r = real ( i - 1, kind = 8 ) / real ( nfine - 1, kind = 8 )
    call boundary_2d ( r, s, x(i), y(i) )
  end do

  call ps_line_open ( nfine, x, y )

  r = 0.0D+00
  do i = nfine, 1, -1
    s = real ( i - 1, kind = 8 ) / real ( nfine - 1, kind = 8 )
    call boundary_2d ( r, s, x(i), y(i) )
  end do

  call ps_line_open ( nfine, x, y )

  s = 0.0D+00
  do i = 1, nfine
    r = real ( i - 1, kind = 8 ) / real ( nfine - 1, kind = 8 )
    call boundary_2d ( r, s, x(i), y(i) )
  end do

  call ps_line_open ( nfine, x, y )
!
!  Finish up the PostScript file.
!
  call ps_page_tail

  call ps_file_tail

  call ps_file_close ( iunit )

  return
end
subroutine draw_pat ( iend, num_end, num_pat, xpat, ypat, re, rw, sn, ss )

!*****************************************************************************80
!
!! DRAW_PAT draws a copy of the pattern.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 July 1999
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
!    Input, integer ( kind = 4 ) IEND(NUM_END), contains the indices of XPAT and YPAT
!    that represent ends of line segments.  The first line segment involves
!    the entries IEND(1) + 1 through IEND(2) of XPAT and YPAT, and the last
!    one involves IEND(NUM_END-1) + 1 through IEND(NUM_END).
!
!    Input, integer ( kind = 4 ) NUM_END, the number of entries used in IEND, and one
!    more than the number of line segments used to draw the pattern.
!
!    Input, integer ( kind = 4 ) NUM_PAT, the number of entries used in XPAT and YPAT.
!
!    Input, real ( kind = 8 ) XPAT(NUM_PAT), YPAT(MAX_PAT), the X and Y
!    coordinates of points that define the line segments used to draw 
!    the pattern.
!
!    Input, real ( kind = 8 ) RE, RW, SN, SS, the extreme (east and west)
!    R coordinates and extreme (north and south) S coordinates of the subcell.
!
  implicit none

  integer ( kind = 4 ), parameter :: MAX_POINT = 20
  integer ( kind = 4 ) num_end
  integer ( kind = 4 ) num_pat

  integer ( kind = 4 ) i
  integer ( kind = 4 ) iend(num_end)
  integer ( kind = 4 ) ipat
  integer ( kind = 4 ) j
  integer ( kind = 4 ) npoint
  real ( kind = 8 ) r
  real ( kind = 8 ) re
  real ( kind = 8 ) rw
  real ( kind = 8 ) s
  real ( kind = 8 ) sn
  real ( kind = 8 ) ss
  real ( kind = 8 ) xpoint(MAX_POINT)
  real ( kind = 8 ) xpat(num_pat)
  real ( kind = 8 ) ypoint(MAX_POINT)
  real ( kind = 8 ) ypat(num_pat)

  ipat = 0
  do i = 1, num_end - 1

    npoint = 0
    do j = iend(i)+1, iend(i+1)
      npoint = npoint + 1
      ipat = ipat + 1
      r = ( xpat(ipat) * re + ( 1.0D+00 - xpat(ipat) ) * rw )
      s = ( ypat(ipat) * sn + ( 1.0D+00 - ypat(ipat) ) * ss )
      call blend1d2d ( r, s, xpoint(npoint), ypoint(npoint) )
    end do

    call ps_line_open ( npoint, xpoint, ypoint )

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
subroutine pattern ( iend, max_end, num_end, max_pat, num_pat, xpat, ypat )

!*****************************************************************************80
!
!! PATTERN defines the tiling pattern.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 July 1999
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
!    Output, integer ( kind = 4 ) IEND(MAX_END), contains the indices of XPAT and YPAT
!    that represent ends of line segments.  The first line segment involves
!    the entries IEND(1) + 1 through IEND(2) of XPAT and YPAT, and the last
!    one involves IEND(NUM_END-1) + 1 through IEND(NUM_END).
!
!    Input, integer ( kind = 4 ) MAX_END, the maximum number of entries in IEND.
!
!    Output, integer ( kind = 4 ) NUM_END, the number of entries used in IEND, and one
!    more than the number of line segments used to draw the pattern.
!
!    Input, integer ( kind = 4 ) MAX_PAT, the maximum number of entries in XPAT and YPAT.
!
!    Output, integer ( kind = 4 ) NUM_PAT, the number of entries used in XPAT and YPAT.
!
!    Output, real ( kind = 8 ) XPAT(MAX_PAT), YPAT(MAX_PAT), the X and Y
!    coordinates of points that define the line segments used to draw 
!    the pattern.
!
  implicit none

  integer ( kind = 4 ) max_end
  integer ( kind = 4 ) max_pat

  integer ( kind = 4 ) iend(max_end)
  integer ( kind = 4 ) num_end
  integer ( kind = 4 ) num_pat
  real ( kind = 8 ) xpat(max_pat)
  real ( kind = 8 ) ypat(max_pat)

  iend(1) = 0

  xpat(1) = 0.0D+00
  ypat(1) = 0.8D+00
  xpat(2) = 0.1D+00
  ypat(2) = 0.83D+00
  xpat(3) = 0.26D+00
  ypat(3) = 0.9D+00
  xpat(4) = 0.3D+00
  ypat(4) = 1.0D+00

  iend(2) = 4

  xpat(5) = 0.1D+00
  ypat(5) = 0.83D+00
  xpat(6) = 0.2D+00
  ypat(6) = 0.6D+00

  iend(3) = 6

  xpat(7) = 0.0D+00
  ypat(7) = 0.4D+00
  xpat(8) = 0.1D+00
  ypat(8) = 0.3D+00
  xpat(9) = 0.2D+00
  ypat(9) = 0.4D+00

  iend(4) = 9

  xpat(10) = 0.1D+00
  ypat(10) = 0.3D+00
  xpat(11) = 0.2D+00
  ypat(11) = 0.3D+00
  xpat(12) = 0.3D+00
  ypat(12) = 0.0D+00

  iend(5) = 12

  xpat(13) = 0.3D+00
  ypat(13) = 0.5D+00
  xpat(14) = 0.4D+00
  ypat(14) = 0.3D+00
  xpat(15) = 0.5D+00
  ypat(15) = 0.3D+00
  xpat(16) = 0.5D+00
  ypat(16) = 0.1D+00
  xpat(17) = 0.6D+00
  ypat(17) = 0.0D+00

  iend(6) = 17

  xpat(18) = 0.4D+00
  ypat(18) = 0.9D+00
  xpat(19) = 0.5D+00
  ypat(19) = 0.8D+00
  xpat(20) = 0.5D+00
  ypat(20) = 0.6D+00
  xpat(21) = 0.6D+00
  ypat(21) = 0.5D+00
  xpat(22) = 0.5D+00
  ypat(22) = 0.3D+00
  xpat(23) = 0.6D+00
  ypat(23) = 0.2D+00
  xpat(24) = 0.7D+00
  ypat(24) = 0.0D+00

  iend(7) = 24

  xpat(25) = 0.7D+00
  ypat(25) = 1.0D+00
  xpat(26) = 0.8D+00
  ypat(26) = 0.8D+00
  xpat(27) = 0.9D+00
  ypat(27) = 0.9D+00
  xpat(28) = 1.0D+00
  ypat(28) = 0.8D+00
  xpat(29) = 0.9D+00
  ypat(29) = 0.7D+00
  xpat(30) = 0.9D+00
  ypat(30) = 0.5D+00
  xpat(31) = 0.9D+00
  ypat(31) = 0.2D+00
  xpat(32) = 1.0D+00
  ypat(32) = 0.4D+00

  iend(8) = 32

  xpat(33) = 0.6D+00
  ypat(33) = 0.5D+00
  xpat(34) = 0.8D+00
  ypat(34) = 0.6D+00
  xpat(35) = 0.9D+00
  ypat(35) = 0.5D+00

  iend(9) = 35

  xpat(36) = 0.5D+00
  ypat(36) = 0.8D+00
  xpat(37) = 0.6D+00
  ypat(37) = 1.0D+00

  iend(10) = 37

  num_pat = 37
  num_end = 10

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
!    Input, integer ( kind = 4 ) NPOINT, the number of points in the line.
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
