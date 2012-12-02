program main

!*****************************************************************************80
!
!! MAIN is the main program for PS_WRITE_PRB.
!
!  Discussion:
!
!    PS_WRITE_PRB is a sample calling program for the PS_WRITE utilities.
!
!  Modified:
!
!    17 November 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PS_WRITE_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the PS_WRITE library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
  call test05 ( )
  call test06 ( )
  call test07 ( )
  call test08 ( )
  call test09 ( )
  call test10 ( )
  call test11 ( )
  call test12 ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PS_WRITE_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( ) 

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 demonstrates how to plot a simple line graph.
!
  implicit none

  real ( kind = 8 ) blue
  character ( len = 80 ) file_name
  real ( kind = 8 ) green
  integer ierror
  integer iunit
  integer npoint
  real ( kind = 8 ) red
  real ( kind = 8 ) x(10)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) y(10)
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Demonstrate a simple line graph.'
!
!  Open the file.
!
  call get_unit ( iunit )

  file_name = 'ps_write_prb01.ps'

  call ps_file_open ( file_name, iunit, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST01'
    write ( *, '(a,i6)' ) '  File creation error ', ierror
    stop
  end if

  call ps_file_head ( file_name )
!
!  Define the size of the page.
!
  xmin = -0.5D+00
  ymin = -0.5D+00

  xmax = 3.0D+00
  ymax = 3.0D+00

  call ps_page_head ( xmin, ymin, xmax, ymax )
!
!  Specify the line color and draw lines.
!
  red = 1.0D+00
  green = 0.0D+00
  blue = 0.0D+00

  call ps_color_line_set ( red, green, blue )

  npoint = 2
  x(1) = 0.0D+00
  y(1) = 0.0D+00
  x(2) = 1.0D+00
  y(2) = 2.0D+00
  call ps_arrow ( x(1), y(1), x(2), y(2) )

  red = 0.0D+00
  green = 1.0D+00
  blue = 0.0D+00

  call ps_color_line_set ( red, green, blue )

  npoint = 3
  x(1) = 0.0D+00
  y(1) = 0.0D+00
  x(2) = 1.0D+00
  y(2) = 1.0D+00
  x(3) = 2.0D+00
  y(3) = 2.0D+00
  call ps_line_open ( npoint, x, y )

  red = 0.0D+00
  green = 0.0D+00
  blue = 1.0D+00

  call ps_color_line_set ( red, green, blue )

  npoint = 2
  x(1) = 1.0D+00
  y(1) = 0.0D+00
  x(2) = 2.0D+00
  y(2) = 2.0D+00
  call ps_line_open ( npoint, x, y )
!
!  Write a label.
!
  red = 0.2D+00
  green = 0.2D+00
  blue = 0.4D+00

  call ps_color_line_set ( red, green, blue )

  call ps_moveto ( 0.5D+00, 0.5D+00 )
  call ps_label ( 'Plot #1' )
!
!  Close up the page and the file.
!
  call ps_page_tail ( )

  call ps_file_tail ( )

  call ps_file_close ( iunit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Normal end of TEST01'

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 plots a curve and marks every 10th point.
!
  implicit none

  integer, parameter :: n = 51

  real ( kind = 8 ) blue
  character ( len = 80 ) file_name
  real ( kind = 8 ) green
  integer i
  integer ierror
  integer iunit
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) red
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  Mark every 10th point of a graph.'

  call get_unit ( iunit )

  file_name = 'ps_write_prb02.ps'

  call ps_file_open ( file_name, iunit, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST02'
    write ( *, '(a,i6)' ) '  File creation error ', ierror
    stop
  end if

  do i = 1, n
    x(i) = ( real ( n - i, kind = 8 ) * 4.0D+00 * pi ) &
      / real ( n - 1, kind = 8 )
  end do

  y(1:n) = sin ( x(1:n) )

  call ps_file_head ( file_name )

  xmin = -1.0D+00
  ymin = - 1.5D+00

  xmax = 4.0D+00 * pi + 1.0D+00
  ymax = + 1.5D+00

  call ps_page_head ( xmin, ymin, xmax, ymax )

  red = 0.0D+00
  green = 0.0D+00
  blue = 0.0D+00

  call ps_color_line_set ( red, green, blue )

  call ps_line_open ( n, x, y )
!
!  Mark every 10-th data point.
!
  do i = 1, n, 10
    call ps_mark_circle ( x(i), y(i) )
  end do

  call ps_page_tail ( )

  call ps_file_tail ( )

  call ps_file_close ( iunit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Normal end of TEST02'

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 plots a curve and a grid.
!
  implicit none

  integer, parameter :: n = 51

  character ( len = 80 ) file_name
  integer i
  integer ierror
  integer iunit
  integer nx
  integer ny
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmaxg
  real ( kind = 8 ) xmin
  real ( kind = 8 ) xming
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymaxg
  real ( kind = 8 ) ymin
  real ( kind = 8 ) yming

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  Plot a curve and a grid.'

  call get_unit ( iunit )

  file_name = 'ps_write_prb03.ps'

  call ps_file_open ( file_name, iunit, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST03'
    write ( *, '(a,i6)' ) '  File creation error ', ierror
    stop
  end if

  do i = 1, n
    x(i) = ( real ( i - 1, kind = 8 ) * 4.0D+00 * pi ) &
      / real ( n - 1, kind = 8 )
  end do

  y(1:n) = sin ( x(1:n) )

  call ps_file_head ( file_name )
!
!  Plot 1: just the curve.
!
  xmin = -1.0D+00
  ymin = - 1.5D+00

  xmax = 4.0D+00 * pi + 1.0D+00
  ymax = + 1.5D+00

  call ps_page_head ( xmin, ymin, xmax, ymax )

  call ps_line_open ( n, x, y )

  call ps_page_tail ( )
!
!  Plot 2: just the grid.
!  And change the mapping.
!
  xmin = 0.0D+00
  ymin = 0.0D+00

  xmax = 4.0D+00 * pi
  ymax = +1.0D+00

  call ps_page_head ( xmin, ymin, xmax, ymax )

  xming = xmin
  xmaxg = xmax
  nx = 21
  yming = -1.0D+00
  ymaxg = +1.0D+00
  ny = 5

  call ps_grid_cartesian ( xming, xmaxg, nx, yming, ymaxg, ny )

  call ps_page_tail
!
!  Plot 3: the curve and the grid.
!
  xmin = -1.0D+00
  ymin = - 1.5D+00

  xmax = 4.0D+00 * pi + 1.0D+00
  ymax = + 1.5D+00

  call ps_page_head ( xmin, ymin, xmax, ymax )

  call ps_line_open ( n, x, y )

  xming = xmin
  xmaxg = xmax
  nx = 21
  yming = -1.0D+00
  ymaxg = +1.0D+00
  ny = 5

  call ps_grid_cartesian ( xming, xmaxg, nx, yming, ymaxg, ny )

  call ps_page_tail ( )

  call ps_file_tail ( )

  call ps_file_close ( iunit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Normal end of TEST03'

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 plots an ellipse and a circle.
!
  implicit none

  integer, parameter :: n = 50

  character ( len = 80 ) file_name
  integer i
  integer ierror
  integer iunit
  integer nx
  integer ny
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) theta
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x0
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmaxg
  real ( kind = 8 ) xmin
  real ( kind = 8 ) xming
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) y0
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymaxg
  real ( kind = 8 ) ymin
  real ( kind = 8 ) yming

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  Draw an ellipse and a circle.'

  call get_unit ( iunit )

  file_name = 'ps_write_prb04.ps'

  call ps_file_open ( file_name, iunit, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST04'
    write ( *, '(a,i6)' ) '  File creation error ', ierror
    stop
  end if

  call ps_file_head ( file_name )

  xmin = - 1.0D+00
  xmax =   5.0D+00
  ymin = - 4.0D+00
  ymax =   2.0D+00

  call ps_page_head ( xmin, ymin, xmax, ymax )
!
!  Get points on the circle, and plot them.
!
  x0 = 2.0D+00
  y0 = -1.0D+00
  r = 2.0D+00

  call circle_points ( x0, y0, r, n, x, y )

  call ps_line_closed ( n, x, y )
!
!  Get points on the ellipse, and plot them.
!
  x0 = 2.0D+00
  y0 = -1.0D+00
  r1 = 3.0D+00
  r2 = 2.0D+00
  theta = pi / 6.0D+00

  call ellipse_points ( x0, y0, r1, r2, theta, n, x, y )

  call ps_line_closed ( n, x, y )

  do i = 1, n, 5
    call ps_mark_circle ( x(i), y(i) )
  end do

  xmin = minval ( x(1:n) )
  xmax = maxval ( x(1:n) )
  ymin = minval ( y(1:n) )
  ymax = maxval ( y(1:n) )
!
!  Draw a grid.
!
  xming = xmin
  xmaxg = xmax
  nx = 21
  yming = ymin
  ymaxg = ymax
  ny = 11

  call ps_grid_cartesian ( xming, xmaxg, nx, yming, ymaxg, ny )

  call ps_page_tail ( )

  call ps_file_tail ( )

  call ps_file_close ( iunit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Normal end of TEST04'

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 plots a graph on a polar grid.
!
  implicit none

  integer, parameter :: n = 50

  real ( kind = 8 ) angle
  character ( len = 80 ) file_name
  real ( kind = 8 ) font_size
  integer i
  integer ierror
  integer iunit
  integer line_width
  integer nr
  integer nt
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  character ( len = 80 ) string
  real ( kind = 8 ) theta
  real ( kind = 8 ) theta1
  real ( kind = 8 ) theta2
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x0
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) y0
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  Draw a bifolium.'

  call get_unit ( iunit )

  file_name = 'ps_write_prb05.ps'

  call ps_file_open ( file_name, iunit, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST05'
    write ( *, '(a,i6)' ) '  File creation error ', ierror
    stop
  end if

  call ps_file_head ( file_name )

  xmin = - 1.25D+00
  xmax =   1.25D+00
  ymin = - 0.25D+00
  ymax =   1.25D+00

  call ps_page_head ( xmin, ymin, xmax, ymax )
!
!  Set the line width.
!
  line_width = 4
  call ps_line_width ( line_width )
!
!  Get points on the cardiod, and plot them.
!
  do i = 1, n
    theta = ( real ( i - 1, kind = 8 ) * pi ) / real ( n, kind = 8 )
    r = 3.0D+00 * sin ( theta ) * ( cos ( theta ) )**2
    x(i) = r * cos ( theta )
    y(i) = r * sin ( theta )
  end do

  call ps_line_closed ( n, x, y )

  line_width = 1
  call ps_line_width ( line_width )

  do i = 1, n, 5
    call ps_mark_circle ( x(i), y(i) )
  end do
!
!  Draw a grid.
!
  line_width = 1
  call ps_line_width ( line_width )

  x0 = 0.0D+00
  y0 = 0.0D+00
  nr = 5
  r1 = 0.2D+00
  r2 = 1.0D+00
  nt = 13
  theta1 = 0.0D+00
  theta2 = 180.0D+00

  call ps_grid_polar ( x0, y0, nr, r1, r2, nt, theta1, theta2 )
!
!  Try sticking a label on.
!
  font_size = 0.50D+00
  call ps_font_size ( font_size )

  x0 = 0.5D+00
  y0 = 0.5D+00
  call ps_moveto ( x0, y0 )

  string = 'Polaris'
  angle = 30.0D+00
  call ps_label_slant ( string, angle )

  call ps_page_tail ( )

  call ps_file_tail ( )

  call ps_file_close ( iunit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Normal end of TEST05'

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 plots a filled ellipse and circle.
!
  implicit none

  integer, parameter :: n = 50

  character ( len = 80 ) file_name
  real ( kind = 8 ) fill_gray
  integer ierror
  integer iunit
  integer line_width
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) theta
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x0
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) y0
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  Draw an ellipse and a circle.'

  call get_unit ( iunit )

  file_name = 'ps_write_prb06.ps'

  call ps_file_open ( file_name, iunit, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST06'
    write ( *, '(a,i6)' ) '  File creation error ', ierror
    stop
  end if

  call ps_file_head ( file_name )

  line_width = 1
  call ps_line_width ( line_width )

  xmin = - 1.0D+00
  xmax =   5.0D+00
  ymin = - 4.0D+00
  ymax =   2.0D+00

  call ps_page_head ( xmin, ymin, xmax, ymax )
!
!  You HAVE to draw the ellipse first, since PostScript just splats
!  each polygon on top of the previous ones.  
!
!  Set the gray fill color.
!
  fill_gray = 0.35D+00
  call ps_fill_gray ( fill_gray )
!
!  Get points on the ellipse, and plot them.
!
  x0 = 2.0D+00
  y0 = -1.0D+00
  r1 = 3.0D+00
  r2 = 2.0D+00
  theta = pi / 6.0D+00

  call ellipse_points ( x0, y0, r1, r2, theta, n, x, y )

  call ps_polygon_fill ( n, x, y )
!
!  Set the gray fill color.
!
  fill_gray = 0.7D+00
  call ps_fill_gray ( fill_gray )
!
!  Get points on the circle, and plot them.
!
  x0 = 2.0D+00
  y0 = -1.0D+00
  r = 2.0D+00

  call circle_points ( x0, y0, r, n, x, y )

  call ps_polygon_fill ( n, x, y )
!
!  Draw a filled circle.
!
  fill_gray = 0.45D+00
  call ps_fill_gray ( fill_gray )

  x0 = 2.0D+00
  y0 = -1.0D+00
  r = 1.0D+00

  call ps_circle_fill ( x0, y0, r )

  call ps_page_tail ( )

  call ps_file_tail ( )

  call ps_file_close ( iunit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Normal end of TEST06'

  return
end
subroutine test07 ( )

!*****************************************************************************80
!
!! TEST07 plots a triangular grid.
!
  implicit none

  character ( len = 80 ) file_name
  real ( kind = 8 ) font_size
  integer ierror
  integer iunit
  integer line_width
  integer n
  character ( len = 100 ) string
  real ( kind = 8 ) x
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) y
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) y3
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a)' ) '  Draw a triangular grid.'

  call get_unit ( iunit )

  file_name = 'ps_write_prb07.ps'

  call ps_file_open ( file_name, iunit, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST07'
    write ( *, '(a,i6)' ) '  File creation error ', ierror
    stop
  end if

  call ps_file_head ( file_name )

  line_width = 1
  call ps_line_width ( line_width )

  xmin =   0.0D+00
  xmax =   2.0D+00
  ymin =   0.0D+00
  ymax =   2.0D+00

  call ps_page_head ( xmin, ymin, xmax, ymax )

  x1 = 0.5D+00
  y1 = 0.1D+00

  x2 = 1.9D+00
  y2 = 1.5D+00

  x3 = 0.1D+00
  y3 = 0.9D+00

  n = 10

  call ps_grid_triangular ( x1, y1, x2, y2, x3, y3, n )
!
!  Try sticking a label on.
!
  font_size = 0.50D+00
  call ps_font_size ( font_size )

  x = 0.5D+00
  y = 0.5D+00
  call ps_moveto ( x, y )

  string = 'Forlorn Hope'
  call ps_label ( string )
!
!  Try writing a smaller label underneath.
!
  font_size = 0.25D+00
  call ps_font_size ( font_size )

  x = 0.5D+00
  y = 0.5D+00 - 2.0D+00 * font_size
  call ps_moveto ( x, y )

  string = '(Well, maybe not)'
  call ps_label ( string )
!
!  Finish up.
!
  call ps_page_tail ( )

  call ps_file_tail ( )

  call ps_file_close ( iunit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Normal end of TEST07'

  return
end
subroutine test08 ( )

!*****************************************************************************80
!
!! TEST08 plots a grid with some filled squares.
!
  implicit none

  integer, parameter :: n = 50

  real ( kind = 8 ) blue
  character ( len = 80 ) file_name
  real ( kind = 8 ) fill_gray
  real ( kind = 8 ) green
  integer i
  integer ierror
  integer iunit
  integer j
  integer line_width
  integer, parameter :: nx = 10
  integer, parameter :: ny = 10
  real ( kind = 8 ) red
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08'
  write ( *, '(a)' ) '  Draw a grid with some boxes filled.'

  call get_unit ( iunit )

  file_name = 'ps_write_prb08.ps'

  call ps_file_open ( file_name, iunit, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST08'
    write ( *, '(a,i6)' ) '  File creation error ', ierror
    stop
  end if

  call ps_file_head ( file_name )

  line_width = 1
  call ps_line_width ( line_width )

  xmin = 0.0D+00
  xmax = 1.0D+00
  ymin = 0.0D+00
  ymax = 1.0D+00

  call ps_page_head ( xmin, ymin, xmax, ymax )
!
!  Set the gray fill color.
!
  fill_gray = 0.20D+00
  call ps_fill_gray ( fill_gray )

  call box ( 1, 1, nx, ny, xmin, xmax, ymin, ymax )
  call box ( 1, 2, nx, ny, xmin, xmax, ymin, ymax )
  call box ( 1, 3, nx, ny, xmin, xmax, ymin, ymax )
  call box ( 2, 3, nx, ny, xmin, xmax, ymin, ymax )
  call box ( 2, 4, nx, ny, xmin, xmax, ymin, ymax )
  call box ( 3, 3, nx, ny, xmin, xmax, ymin, ymax )
  call box ( 3, 4, nx, ny, xmin, xmax, ymin, ymax )
  call box ( 4, 2, nx, ny, xmin, xmax, ymin, ymax )
  call box ( 4, 3, nx, ny, xmin, xmax, ymin, ymax )

  fill_gray = 0.40D+00
  call ps_fill_gray ( fill_gray )

  call box ( 3, 7, nx, ny, xmin, xmax, ymin, ymax )
  call box ( 3, 8, nx, ny, xmin, xmax, ymin, ymax )
  call box ( 3, 9, nx, ny, xmin, xmax, ymin, ymax )
  call box ( 4, 7, nx, ny, xmin, xmax, ymin, ymax )
  call box ( 4, 8, nx, ny, xmin, xmax, ymin, ymax )
  call box ( 4, 9, nx, ny, xmin, xmax, ymin, ymax )
  call box ( 5, 6, nx, ny, xmin, xmax, ymin, ymax )
  call box ( 5, 7, nx, ny, xmin, xmax, ymin, ymax )
  call box ( 5, 9, nx, ny, xmin, xmax, ymin, ymax )
  call box ( 6, 5, nx, ny, xmin, xmax, ymin, ymax )
  call box ( 6, 6, nx, ny, xmin, xmax, ymin, ymax )
  call box ( 6, 9, nx, ny, xmin, xmax, ymin, ymax )
  call box ( 7, 9, nx, ny, xmin, xmax, ymin, ymax )

  fill_gray = 0.60D+00
  call ps_fill_gray ( fill_gray )

  call box ( 8, 3, nx, ny, xmin, xmax, ymin, ymax )
  call box ( 8, 4, nx, ny, xmin, xmax, ymin, ymax )
  call box ( 9, 2, nx, ny, xmin, xmax, ymin, ymax )
  call box ( 9, 3, nx, ny, xmin, xmax, ymin, ymax )
  call box ( 9, 4, nx, ny, xmin, xmax, ymin, ymax )
  call box (10, 4, nx, ny, xmin, xmax, ymin, ymax )
  call box (10, 5, nx, ny, xmin, xmax, ymin, ymax )

  fill_gray = 0.80
  call ps_fill_gray ( fill_gray )

  call box ( 9,10, nx, ny, xmin, xmax, ymin, ymax )
  call box (10, 9, nx, ny, xmin, xmax, ymin, ymax )
  call box (10,10, nx, ny, xmin, xmax, ymin, ymax )
!
!  Draw the grid on top.
!
  red = 0.0D+00
  green = 0.0D+00
  blue = 1.0D+00

  call ps_color_line_set ( red, green, blue )

  do i = 0, nx
    x = real ( i, kind = 8 ) * xmax / real ( nx, kind = 8 )
    y = ymin
    call ps_moveto ( x, y )
    y = ymax
    call ps_lineto ( x, y )
  end do

  do j = 0, ny
    x = xmin
    y = real ( j, kind = 8 ) * ymax / real ( ny, kind = 8 )
    call ps_moveto ( x, y )
    x = xmax
    call ps_lineto ( x, y )
  end do

  call ps_page_tail ( )

  call ps_file_tail ( )

  call ps_file_close ( iunit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Normal end of TEST08'

  return
end
subroutine test09 ( )

!*****************************************************************************80
!
!! TEST09 plots a graph on a polar grid to an encapsulated PostScript file.
!
  implicit none

  integer, parameter :: n = 50

  real ( kind = 8 ) angle
  character ( len = 80 ) file_name
  real ( kind = 8 ) font_size
  integer i
  integer ierror
  integer iunit
  integer line_width
  integer nr
  integer nt
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  character ( len = 80 ) string
  real ( kind = 8 ) theta
  real ( kind = 8 ) theta1
  real ( kind = 8 ) theta2
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) x0
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) y0
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST09'
  write ( *, '(a)' ) '  Draw a bifolium.'
  write ( *, '(a)' ) '  Write it in an encapsulated PostScript file.'

  call get_unit ( iunit )

  file_name = 'ps_write_prb09.eps'

  call ps_file_open ( file_name, iunit, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST09'
    write ( *, '(a,i6)' ) '  File creation error ', ierror
    stop
  end if

  call eps_file_head ( file_name, 36, 36, 576, 756 )

  xmin = - 1.25D+00
  xmax =   1.25D+00
  ymin = - 0.25D+00
  ymax =   1.25D+00

  call ps_page_head ( xmin, ymin, xmax, ymax )
!
!  Set the line width.
!
  line_width = 4
  call ps_line_width ( line_width )
!
!  Get points on the cardioid, and plot them.
!
  do i = 1, n
    theta = ( real ( i - 1, kind = 8 ) * pi ) / real ( n, kind = 8 )
    r = 3.0D+00 * sin ( theta ) * ( cos ( theta ) )**2
    x(i) = r * cos ( theta )
    y(i) = r * sin ( theta )
  end do

  call ps_line_closed ( n, x, y )

  line_width = 1
  call ps_line_width ( line_width )

  do i = 1, n, 5
    call ps_mark_circle ( x(i), y(i) )
  end do
!
!  Draw a grid.
!
  line_width = 1
  call ps_line_width ( line_width )

  x0 = 0.0D+00
  y0 = 0.0D+00
  nr = 5
  r1 = 0.2D+00
  r2 = 1.0D+00
  nt = 13
  theta1 = 0.0D+00
  theta2 = 180.0D+00

  call ps_grid_polar ( x0, y0, nr, r1, r2, nt, theta1, theta2 )
!
!  Try sticking a label on.
!
  font_size = 0.50D+00
  call ps_font_size ( font_size )

  x0 = 0.5D+00
  y0 = 0.5D+00
  call ps_moveto ( x0, y0 )

  string = 'Polaris'
  angle = 30.0D+00
  call ps_label_slant ( string, angle )

  call ps_page_tail ( )

  call eps_file_tail ( )

  call ps_file_close ( iunit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Normal end of TEST09'

  return
end
subroutine test10 ( )

!*****************************************************************************80
!
!! TEST10 tests PS_SETTING_PRINT.
!
  implicit none

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST10'
  write ( *, '(a)' ) '  PS_SETTING_PRINT prints the current PS_WRITE settings.'

  call ps_setting_print ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Normal end of TEST10'

  return
end
subroutine test11 ( )

!*****************************************************************************80
!
!! TEST11 tests PS_MARK_CIRCLES, PS_MARK_DISKS.
!
  implicit none

  integer, parameter :: n = 20

  real ( kind = 8 ) blue
  character ( len = 80 ) file_name
  real ( kind = 8 ) green
  integer i
  integer ierror
  integer iunit
  integer marker_size
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r
  real ( kind = 8 ) red
  real ( kind = 8 ) theta
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST11'
  write ( *, '(a)' ) '  PS_MARK_CIRLES marks points with an open circle.'
  write ( *, '(a)' ) '  PS_MARK_DISKS marks points with an closed disk.'

  call get_unit ( iunit )

  file_name = 'ps_write_prb11.eps'

  call ps_file_open ( file_name, iunit, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST11'
    write ( *, '(a,i6)' ) '  File creation error ', ierror
    stop
  end if

  call eps_file_head ( file_name, 36, 36, 576, 756 )

  xmin = - 2.25D+00
  xmax =   2.25D+00
  ymin = - 2.25D+00
  ymax =   2.25D+00

  call ps_page_head ( xmin, ymin, xmax, ymax )

  do i = 1, n
    theta = ( real ( i - 1, kind = 8 ) * pi ) / real ( n, kind = 8 )
    r = 3.0D+00 * sin ( theta ) * ( cos ( theta ) )**2
    x(i) = r * cos ( theta )
    y(i) = r * sin ( theta )
  end do

  red = 0.0D+00
  green = 0.0D+00
  blue = 1.0D+00

  call ps_color_line_set ( red, green, blue )

  marker_size = 4
  call ps_marker_size ( marker_size )

  call ps_mark_circles ( n, x, y )

  do i = 1, n
    theta = ( real ( 2 * ( i - 1 ), kind = 8 ) * pi ) &
      / real ( n, kind = 8 )
    r = 2.0D+00
    x(i) = r * cos ( theta )
    y(i) = r * sin ( theta )
  end do

  red = 1.0D+00
  green = 0.0D+00
  blue = 0.0D+00

  call ps_color_fill_set ( red, green, blue )

  marker_size = 8
  call ps_marker_size ( marker_size )

  call ps_mark_disks ( n, x, y )

  call ps_page_tail ( )

  call eps_file_tail ( )

  call ps_file_close ( iunit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Normal end of TEST11'

  return
end
subroutine box ( i, j, nx, ny, xmin, xmax, ymin, ymax )

!*****************************************************************************80
!
!! BOX fills in the (I,J) box on a grid.
!
  implicit none

  integer, parameter :: n = 4

  integer i
  integer j
  integer nx
  integer ny
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin

  x(1) = real ( i - 1, kind = 8 ) * xmax / real ( nx, kind = 8 )
  y(1) = real ( j - 1, kind = 8 ) * ymax / real ( ny, kind = 8 )
  x(2) = real ( i,     kind = 8 ) * xmax / real ( nx, kind = 8 )
  y(2) = real ( j - 1, kind = 8 ) * ymax / real ( ny, kind = 8 )
  x(3) = real ( i,     kind = 8 ) * xmax / real ( nx, kind = 8 )
  y(3) = real ( j,     kind = 8 ) * ymax / real ( ny, kind = 8 )
  x(4) = real ( i - 1, kind = 8 ) * xmax / real ( nx, kind = 8 )
  y(4) = real ( j,     kind = 8 ) * ymax / real ( ny, kind = 8 )

  call ps_polygon_fill ( n, x, y )

  return
end
subroutine test12 ( )

!*****************************************************************************80
!
!! TEST12 tests PS_MARK_POINT.
!
  implicit none

  character ( len = 80 ) file_name
  integer i
  integer ierror
  integer iunit
  real ( kind = 8 ) x(2)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST12'	
  write ( *, '(a)' ) '  PS_MARK_POINT marks points with a tiny point.'

  call get_unit ( iunit )

  file_name = 'ps_write_prb12.eps'

  call ps_file_open ( file_name, iunit, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST12'
    write ( *, '(a,i6)' ) '  File creation error ', ierror
    stop
  end if

  call eps_file_head ( file_name, 36, 36, 576, 756 )

  xmin = 0.0D+00
  xmax = 1.0D+00
  ymin = 0.0D+00
  ymax = 1.0D+00

  call ps_page_head ( xmin, ymin, xmax, ymax )

  do i = 1, 10000
    call random_number ( harvest = x(1:2) )
    x(1:2) = sqrt ( x(1:2) )
    call ps_mark_point ( x(1), x(2) )
  end do

  call ps_page_tail ( )

  call eps_file_tail ( )
 
  call ps_file_close ( iunit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Normal end of TEST12'

  return
end
