module AA_point_data_module

!*****************************************************************************80
!
!! AA_POINT_DATA_MODULE keeps track of some point data.
!
!  Discussion:
!
!    The name of this module has 'AA' preprended to it so that
!    when this file is split up, and the pieces compiled, this
!    piece will be compiled before the others.  That's necessary
!    so that routines that use this module will not be compiled
!    first!  Ah, the joys of learning modules!
!
!    The arrays stored here are allocatable.  The documentation
!    lists their intended dimensions, which depend on NODE_NUM
!    and DIM_NUM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Module data, integer ( kind = 4 ) CONNECT(NODE_NUM), connection indicator.
!    If CONNECT(I) is:
!    0, it is an isolated point;
!    1, it is connected to the previous point only;
!    2, it is connected to the next point only;
!    3, it is connected to the previous and next points.
!
!    Module data, real ( kind = 8 ) COORD(DIM_NUM,NODE_NUM), the point
!    coordinates.
!
!    Module data, real ( kind = 8 ) COORD_MIN(DIM_NUM), COORD_MAX(DIM_NUM), 
!    the data ranges.
!
!    Module data, real ( kind = 8 ) RADIUS(NODE_NUM), the radius of a circle
!    around each point.
!
!    Module data, integer ( kind = 4 ) TAG(NODE_NUM), a tag for each point,
!    presumably the cluster index.
!
  implicit none

  integer ( kind = 4 ), allocatable, dimension ( : ) :: connect
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: coord
  real ( kind = 8 ), allocatable, dimension ( : ) :: coord_max
  real ( kind = 8 ), allocatable, dimension ( : ) :: coord_min
  real ( kind = 8 ), allocatable, dimension ( : ) :: radius
  integer ( kind = 4 ), allocatable, dimension ( : ) :: tag

end
program main

!*****************************************************************************80
!
!! MAIN is the main program for PLOT_POINTS.
!
!  Discussion:
!
!    PLOT_POINTS plots the points in a file.
!
!    Most recent change was that SHADOW is no longer a "toggle"
!    command.  Rather, SHADOW turns shadowing on, and NOSHADOW turns
!    it off.
!
!    Also, I started to set up a LINE command, but never got done with it.
!    Any shoemaker's elves with nothing else to do, please come fix this!
!
!  Usage:
!
!    The program can be invoked by:
!
!      plot_points plot_file_name
!
!    or:
!
!      plot_points
!
!    in which case the user will be asked to supply the input file name.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 November 2003
!
!  Author:
!
!    John Burkardt
!
  use AA_point_data_module

  implicit none

  integer ( kind = 4 ), parameter :: line_max = 10

  integer ( kind = 4 ) arg_num
  logical :: box_requested = .false.
  real ( kind = 8 ), dimension ( 2, 2 ) :: box_xy = reshape ( &
    (/ 0.0D+00, 0.0D+00, 1.0D+00, 1.0D+00 /), (/ 2, 2 /) )
  logical :: circle_requested = .false.
  real ( kind = 8 ) :: circle_r = 0.0D+00
  real ( kind = 8 ) :: circle_x = 0.0D+00
  real ( kind = 8 ) :: circle_y = 0.0D+00
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: coord2
  character ( len = 100 ) command
  logical, parameter :: delaunay = .true.
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) dim_num2
  character ( len = 255 ) :: file_in_name = ' '
  integer ( kind = 4 ) file_in_unit
  character ( len = 255 ) :: file_out_name = ' '
  character ( len = 255 ) :: file_out_name_user = ' '
  character ( len = 255 ) :: file_out_name_default = 'plot000.eps'
  integer ( kind = 4 ) first
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iarg
  integer ( kind = 4 ) iargc
  integer ( kind = 4 ) :: ierror = 0
  integer ( kind = 4 ) ilen
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) ipxfargc
  integer ( kind = 4 ) last
  integer ( kind = 4 ) length
  integer ( kind = 4 ) lens
  integer ( kind = 4 ) :: line_num = 0
  real ( kind = 8 ), dimension ( 2, 2, line_max ) :: line_xy
  integer ( kind = 4 ) :: marker_size = 4
  integer ( kind = 4 ) next
  real ( kind = 8 ) plot_max(2)
  real ( kind = 8 ) plot_min(2)
  character ( len = 255 ) point_file
  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) node_num2
  logical s_eqi
  logical :: shadow = .false.
  character ( len = 255 ) :: title = ' '
  real ( kind = 8 ) tol
  logical, parameter :: voronoi = .true.
  logical, parameter :: walls = .true.
  integer ( kind = 4 ) :: x_index = 1
  integer ( kind = 4 ) :: y_index = 2

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PLOT_POINTS'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Given a file of 2D points, make a dot plot, and'
  write ( *, '(a)' ) '  a Delaunay triangulation plot.'
!
!  Get the number of command line arguments.
!
!  Old style:
!
  arg_num = iargc ( )
!
!  New style:
!
! arg_num = ipxfargc ( )
!
!  If at least one command line argument, it's the input file name.
!
  if ( 1 <= arg_num ) then

    iarg = 1
!
!  Old style:
!
    call getarg ( iarg, file_in_name )
!
!  New style:
!
!   call pxfgetarg ( iarg, file_in_name, ilen, ierror )
!
!   if ( ierror /= 0 ) then
!     write ( *, '(a)' ) ' '
!     write ( *, '(a)' ) 'PLOT_POINTS - Fatal error!'
!     write ( *, '(a)' ) '  Could not read command line argument.'
!     stop
!   end if

    call data_read ( dim_num, file_in_name, node_num )

  end if
!
!  Command loop begins here.
!
  do

    if ( command(1:1) /= '#' ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'Enter command (H for help):'
      write ( *, '(a)' ) ' '
    end if

    read ( *, '(a)' ) command
!
!  # is a comment.
!  Echo it in case we're saving the output.
!
    if ( command(1:1) == '#' ) then

      write ( *, '(a)' ) trim ( command )
!
!  BALLOON plot
!  Plot each point with its maximum circle.
!
    else if ( s_eqi ( command(1:3), 'BAL' ) ) then

      if ( len_trim ( file_out_name_user ) == 0 ) then

        call file_name_inc ( file_out_name_default )
        file_out_name = file_out_name_default

      else

        file_out_name = file_out_name_user
        file_out_name_user = ' '

      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Write balloon plot data to "' // &
        trim ( file_out_name ) // '".'

      dim_num2 = 2

      allocate ( coord2(1:2,1:node_num) )

      coord2(1,1:node_num) = coord(x_index,1:node_num)
      coord2(2,1:node_num) = coord(y_index,1:node_num)

      call radius_maximus ( dim_num2, node_num, coord2, walls, radius )

      deallocate ( coord2 )

      call plot_size ( dim_num, box_requested, box_xy, circle_requested, &
        circle_x, circle_y, circle_r, x_index, y_index, plot_min, plot_max )

      call balloon_plot ( box_requested, box_xy, circle_requested, circle_x, &
        circle_y, circle_r, file_out_name, dim_num, &
        node_num, coord, plot_min, plot_max, x_index, y_index, shadow, &
        title, radius )
!
!  BOX x1, y1, x2, y2
!  specify a box to be drawn.
!
    else if ( s_eqi ( command(1:3), 'BOX' ) ) then

      box_requested = .true.

      command(1:3) = '   '
      command = adjustl ( command )

      if ( command(4:4) == '=' ) then
        command(4:4) = ' '
      end if

      if ( len_trim ( command ) /= 0 ) then

        call s_to_r8vec ( command, 4, box_xy(1:2,1:2), ierror )
 
      else

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Enter lower left (x,y), then upper right (x,y):'
        read ( *, *, iostat = ios ) box_xy(1:2,1:2)

        if ( ios /= 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'PLOT_POINTS - Fatal error!'
          write ( *, '(a)' ) '  I/O error reading box coordinates.'
          stop
        end if

      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a,2f9.4)' ) &
        '  Lower left box corner  = ', box_xy(1,1), box_xy(2,1)
      write ( *, '(a,2f9.4)' ) &
       '  Upper right box corner =  ', box_xy(1,2), box_xy(2,2)
!
!  CIRCLE x y r
!  Plot a single circle of given size.
!
    else if ( s_eqi ( command(1:6), 'CIRCLE' ) ) then

      circle_requested = .true.

      command(1:6) = '   '
      command = adjustl ( command )

      if ( command(7:7) == '=' ) then
        command(7:7) = ' '
        command = adjustl ( command )
      end if

      if ( len_trim ( command ) /= 0 ) then

        first = 1
        call s_to_r8 ( command(first:), circle_x, ierror, length )
        first = first + length
        call s_to_r8 ( command(first:), circle_y, ierror, length )
        first = first + length
        call s_to_r8 ( command(first:), circle_r, ierror, length )

      else

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Enter center (x,y), then radius r:'
        read ( *, *, iostat = ios ) circle_x, circle_y, circle_r

        if ( ios /= 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'PLOT_POINTS - Fatal error!'
          write ( *, '(a)' ) '  I/O error reading circle coordinates.'
          stop
        end if

      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a,2f9.4)' ) '  Circle center = ', circle_x, circle_y
      write ( *, '(a,2f9.4)' ) '  Circle radius = ', circle_r
!
!  DASH
!  Create a dash plot of connected points.
!
    else if ( s_eqi ( command(1:3), 'DAS' ) ) then

      if ( len_trim ( file_out_name_user ) == 0 ) then

        call file_name_inc ( file_out_name_default )
        file_out_name = file_out_name_default

      else

        file_out_name = file_out_name_user
        file_out_name_user = ' '

      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Write dash data to "' // &
        trim ( file_out_name ) // '".'
 
      call plot_size ( dim_num, box_requested, box_xy, circle_requested, &
        circle_x, circle_y, circle_r, x_index, y_index, plot_min, plot_max )

      call dash_plot ( box_requested, box_xy, circle_requested, &
        circle_x, circle_y, circle_r, connect, file_out_name, &
        dim_num, node_num, coord, plot_min, plot_max, x_index, y_index, &
        shadow, title )
!
!  DELAUNAY
!  Create a Delaunay plot.
!
    else if ( s_eqi ( command(1:3), 'DEL' ) ) then

      if ( len_trim ( file_out_name_user ) == 0 ) then

        call file_name_inc ( file_out_name_default )
        file_out_name = file_out_name_default

      else

        file_out_name = file_out_name_user
        file_out_name_user = ' '

      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Create Delaunay plot "' // &
        trim ( file_out_name ) // '".'

      call plot_size ( dim_num, box_requested, box_xy, circle_requested, &
        circle_x, circle_y, circle_r, x_index, y_index, plot_min, plot_max )

      call delaunay_plot ( box_requested, box_xy, circle_requested, &
        circle_x, circle_y, circle_r, file_out_name, &
        dim_num, node_num, coord, plot_min, plot_max, x_index, y_index, &
        shadow, title )
!
!  DOT
!  Create a dot plot of points.
!
    else if ( s_eqi ( command(1:3), 'DOT' ) ) then

      if ( len_trim ( file_out_name_user ) == 0 ) then

        call file_name_inc ( file_out_name_default )
        file_out_name = file_out_name_default

      else

        file_out_name = file_out_name_user
        file_out_name_user = ' '

      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Write dot data to "' // &
        trim ( file_out_name ) // '".'

      call plot_size ( dim_num, box_requested, box_xy, circle_requested, &
        circle_x, circle_y, circle_r, x_index, y_index, plot_min, plot_max )

      call dot_plot ( box_requested, box_xy, circle_requested, circle_x, &
        circle_y, circle_r, file_out_name, dim_num, &
        node_num, coord, plot_min, plot_max, x_index, y_index, &
        shadow, title, marker_size )
!
!  H
!
    else if ( s_eqi ( command(1:1), 'H' ) ) then

      write ( *, '(a)'    ) ' '
      write ( *, '(a)'    ) 'BALLOON          a "balloon" plot of points.'
      write ( *, '(a)'    ) 'DASH             dash plot of connected points.'
      write ( *, '(a)'    ) 'DOT              dot plot of points.'
      write ( *, '(a)'    ) 'DEL              Delaunay triangulation plot'
      write ( *, '(a)'    ) 'KM               K-Means plot.'
      write ( *, '(a)'    ) 'TH               thin the points.'
      write ( *, '(a)'    ) 'TV               triangulated Voronoi diagram.'
      write ( *, '(a)'    ) 'VOR              Voronoi diagram plot.'
      write ( *, '(a)'    ) ' '
      write ( *, '(a)'    ) 'Q                quit.'
      write ( *, '(a)'    ) 'READ filename    read another input file.'
      write ( *, '(a)'    ) 'OUTPUT filename  name the next output file.'
      write ( *, '(a)'    ) 'BOX x1 y1 x2 y2  draw a single box in the plot.'
      write ( *, '(a)'    ) 'NOBOX            do not draw a box in the plot.'
      write ( *, '(a)'    ) 'CIRCLE x y r     draw a single circle in the plot.'
      write ( *, '(a)'    ) 'NOCIRCLE         do not draw a circle in the plot.'
      write ( *, '(a)'    ) 'LINE x1 y1 x2 y2 draw a line in the plot.'
      write ( *, '(a)'    ) 'NOLINE           do not draw lines in the plot.'
      write ( *, '(a)'    ) 'MARKER_SIZE = n  Specify marker size.'
      write ( *, '(a,i6)' ) '                 Current value = ', marker_size
      write ( *, '(a)'    ) 'SHADOW           mark X and Y axes for each point.'
      write ( *, '(a)'    ) 'NOSHADOW         cancel shadow command.'
      write ( *, '(a)'    ) 'TITLE            enter a title for the next plot.'
      write ( *, '(a)'    ) 'X = n            specify data index to use as X.'
      write ( *, '(a,i6)' ) '                 Current value = ', x_index
      write ( *, '(a)'    ) 'Y = n            specify data index to use as Y.'
      write ( *, '(a,i6)' ) '                 Current value = ', y_index
      write ( *, '(a)'    ) '# comment        make a comment.'
!
!  KMEANS
!  Create a K-Means plot.
!  We presume the cluster centers are each set off with blank records.
!  We presume the data was tagged with a third cluster coordinate.
!
    else if ( s_eqi ( command(1:1), 'K' ) ) then

      if ( len_trim ( file_out_name_user ) == 0 ) then

        call file_name_inc ( file_out_name_default )
        file_out_name = file_out_name_default

      else

        file_out_name = file_out_name_user
        file_out_name_user = ' '

      end if

      tag(1:node_num) = nint ( coord(3,1:node_num) )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Create a K-Means plot "' // &
        trim ( file_out_name ) // '".'

      call plot_size ( dim_num, box_requested, box_xy, circle_requested, &
        circle_x, circle_y, circle_r, x_index, y_index, plot_min, plot_max )

      call kmeans_plot ( box_requested, box_xy, circle_requested, &
        circle_x, circle_y, circle_r, file_out_name, dim_num, &
        node_num, coord, plot_min, plot_max, x_index, y_index, connect, &
        tag, shadow, title )
!
!  LINE x1, y1, x2, y2
!  specify a line to be drawn.
!
    else if ( s_eqi ( command(1:4), 'LINE' ) ) then

      line_num = line_num + 1

      command(1:4) = '   '
      command = adjustl ( command )

      if ( command(5:5) == '=' ) then
        command(5:5) = ' '
        command = adjustl ( command )
      end if

      if ( len_trim ( command ) /= 0 ) then

        first = 1
        call s_to_r8 ( command(first:), line_xy(1,1,line_num), ierror, length )
        first = first + length
        call s_to_r8 ( command(first:), line_xy(2,1,line_num), ierror, length )
        first = first + length
        call s_to_r8 ( command(first:), line_xy(1,2,line_num), ierror, length )
        first = first + length
        call s_to_r8 ( command(first:), line_xy(2,2,line_num), ierror, length )

      else

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Enter first (x,y), then second (x,y):'
        read ( *, *, iostat = ios ) line_xy(1:2,1:2,line_num)

        if ( ios /= 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'PLOT_POINTS - Fatal error!'
          write ( *, '(a)' ) '  I/O error reading line coordinates.'
          stop
        end if

      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a,2f9.4)' ) &
        '  Line begins at ', line_xy(1,1,line_num), line_xy(2,1,line_num)
      write ( *, '(a,2f9.4)' ) &
       '  Line ends at    ', line_xy(1,2,line_num), line_xy(2,2,line_num)
!
!  MARKER_SIZE = n
!
    else if ( s_eqi ( command(1:11), 'MARKER_SIZE' ) ) then

      call s_blank_delete ( command )

      if ( command(12:12) == '=' ) then
        call s_to_i4 ( command(13:), marker_size, ierror, length )
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'Enter new value for marker size:'
        read ( *, * ) marker_size
      end if

      write ( *, '(a,i6)' ) 'Marker size set to ', marker_size
!
!  NOBOX
!  specify no box to be drawn.
!
    else if ( s_eqi ( command(1:5), 'NOBOX' ) ) then

      box_requested = .false.
!
!  NOCIRCLE
!  specify no circle to be drawn.
!
    else if ( s_eqi ( command(1:8), 'NOCIRCLE' ) ) then

      circle_requested = .false.
!
!  NOLINE
!  specify no lines to be drawn.
!
    else if ( s_eqi ( command(1:6), 'NOLINE' ) ) then

      line_num = 0
!
!  NOSHADOW = Turn "shadow" OFF
!
    else if ( s_eqi ( command(1:4), 'NOSH' ) ) then

      shadow = .false.

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SHADOWING is turned OFF.'
!
!  OUTPUT filename
!
    else if ( s_eqi ( command(1:6), 'OUTPUT' ) ) then

      call s_blank_delete ( command )

      if ( command(7:7) == '=' ) then
        file_out_name_user = command(8:)
      else
        file_out_name_user = command(7:)
      end if

      if ( len_trim ( file_out_name_user ) == 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Enter the name of the next output file.'
        read ( *, '(a)' ) file_out_name_user
      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) &
        '  Next output file is ' // trim ( file_out_name_user )
!
!  Q = Quit
!
    else if ( s_eqi ( command(1:1), 'Q' ) ) then

      exit
!
!  READ
!
    else if ( s_eqi ( command(1:4), 'READ' ) ) then

      call s_blank_delete ( command )

      if ( command(5:5) == '=' ) then
        file_in_name = command(6:)
      else
        file_in_name = command(5:)
      end if

      if ( len_trim ( file_in_name ) == 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Enter the name of the new input file.'
        read ( *, '(a)' ) file_in_name
      end if

      if ( allocated ( connect ) ) then
        deallocate ( connect )
      end if

      if ( allocated ( coord ) ) then
        deallocate ( coord )
      end if

      if ( allocated ( coord_max ) ) then
        deallocate ( coord_max )
      end if

      if ( allocated ( coord_min ) ) then
        deallocate ( coord_min )
      end if

      if ( allocated ( radius ) ) then
        deallocate ( radius )
      end if

      if ( allocated ( tag ) ) then
        deallocate ( tag )
      end if

      dim_num = 0
      node_num = 0

      call data_read ( dim_num, file_in_name, node_num )
!
!  SH = Shadow
!
    else if ( s_eqi ( command(1:2), 'SH' ) ) then

      shadow = .true.

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SHADOWING is turned ON.'
      write ( *, '(a)' ) '  Each point makes a small mark on X and Y axis.'
!
!  TH = Thin
!
    else if ( s_eqi ( command(1:2), 'TH' ) ) then

      tol = 0.01D+00 * maxval ( coord_max(1:dim_num) - coord_min(1:dim_num) )

      call points_thin ( connect, dim_num, node_num, coord, tol, node_num2 )
      node_num = node_num2
!
!  TITLE
!
    else if ( s_eqi ( command(1:5), 'TITLE' ) ) then

      next = 6
      length = len_trim ( command )

      do

        next = next + 1

        if ( length < next ) then
          exit
        end if

        if ( command(next:next) /= ' ' .and. command(next:next) /= '=' ) then
          exit
        end if

      end do

      if ( next <= length ) then
        title = command(next:)
      else
        title = ' '
      end if

      if ( len_trim ( title ) == 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Enter the title:'
        read ( *, '(a)' ) title
      end if

      title = adjustl ( title )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TITLE:'
      write ( *, '(a)' ) '  "' // trim ( title ) // '"'
!
!  TV
!  Create a triangulated Voronoi plot.
!
    else if ( s_eqi ( command(1:2), 'TV' ) ) then

      if ( len_trim ( file_out_name_user ) == 0 ) then

        call file_name_inc ( file_out_name_default )
        file_out_name = file_out_name_default

      else

        file_out_name = file_out_name_user
        file_out_name_user = ' '

      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Create triangulated Voronoi plot "' &
        // trim ( file_out_name ) // '".'

      call plot_size ( dim_num, box_requested, box_xy, circle_requested, &
        circle_x, circle_y, circle_r, x_index, y_index, plot_min, plot_max )

      call trivor_plot ( box_requested, box_xy, circle_requested, &
        circle_x, circle_y, circle_r, file_out_name, dim_num, &
        node_num, coord, plot_min, plot_max, x_index, y_index, &
        shadow, title )
!
!  VORONOI
!  Create a Voronoi plot.
!
    else if ( s_eqi ( command(1:3), 'VOR' ) ) then

      if ( len_trim ( file_out_name_user ) == 0 ) then

        call file_name_inc ( file_out_name_default )
        file_out_name = file_out_name_default

      else

        file_out_name = file_out_name_user
        file_out_name_user = ' '

      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Create Voronoi plot "' // &
        trim ( file_out_name ) // '".'

      call plot_size ( dim_num, box_requested, box_xy, circle_requested, &
        circle_x, circle_y, circle_r, x_index, y_index, plot_min, plot_max )

      call voronoi_plot ( box_requested, box_xy, circle_requested, &
        circle_x, circle_y, circle_r, file_out_name, dim_num, &
        node_num, coord, plot_min, plot_max, x_index, y_index, &
        shadow, title )
!
!  X_INDEX
!
    else if ( s_eqi ( command(1:1), 'X' ) ) then

      call s_blank_delete ( command )

      if ( command(2:2) == '=' ) then
        call s_to_i4 ( command(3:), x_index, ierror, length )
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'Enter new value for X index:'
        read ( *, * ) x_index
      end if

      write ( *, '(a,i6)' ) 'X index set to ', x_index
!
!  Y_INDEX
!
    else if ( s_eqi ( command(1:1), 'Y' ) ) then

      call s_blank_delete ( command )

      if ( command(2:2) == '=' ) then
        call s_to_i4 ( command(3:), y_index, ierror, length )
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'Enter new value for Y index:'
        read ( *, * ) y_index
      end if

      write ( *, '(a,i6)' ) 'Y index set to ', y_index
!
!  Unrecognized command.
!
    else

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PLOT_POINTS: Unrecognized command:'
      write ( *, '(a)' ) '  "'// trim ( command ) // '"'
 
    end if

  end do
!
!  Deallocate arrays.
!
  if ( allocated ( connect ) ) then
    deallocate ( connect )
  end if

  if ( allocated ( coord ) ) then
    deallocate ( coord )
  end if

  if ( allocated ( coord_max ) ) then
    deallocate ( coord_max )
  end if

  if ( allocated ( coord_min ) ) then
    deallocate ( coord_min )
  end if

  if ( allocated ( radius ) ) then
    deallocate ( radius )
  end if

  if ( allocated ( tag ) ) then
    deallocate ( tag )
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PLOT_POINTS'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine angle_contains_ray_2d ( inside, x1, y1, x2, y2, x3, y3, x, y )

!*****************************************************************************80
!
!! ANGLE_CONTAINS_RAY_2D determines if an angle contains a ray, in 2D.
!
!  Discussion:
!
!    The angle is defined by the sequence of points (X1,Y1), (X2,Y2)
!    and (X3,Y3).
!
!    The ray is defined by the sequence of points (X2,Y2), (X,Y).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X1, Y1, X2, Y2, X3, Y3, the X and Y coordinates of
!    the angle.
!
!    Input, real ( kind = 8 ) X, Y, the end point of the ray to be checked.
!    The ray is assumed to have an origin at (X2,Y2).
!
!    Output, logical INSIDE, is .TRUE. if the ray is inside
!    the angle or on its boundary, and .FALSE. otherwise.
!
  implicit none

  real ( kind = 8 ) a1
  real ( kind = 8 ) a2
  real ( kind = 8 ) angle_deg_2d
  logical inside
  real ( kind = 8 ) x
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ) y
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) y3

  a1 = angle_deg_2d ( x1, y1, x2, y2, x, y )
  a2 = angle_deg_2d ( x1, y1, x2, y2, x3, y3 )

  if ( a1 <= a2 ) then
    inside = .true.
  else
    inside = .false.
  end if

  return
end
function angle_deg_2d ( x1, y1, x2, y2, x3, y3 )

!*****************************************************************************80
!
!! ANGLE_DEG_2D returns the angle swept out between two rays in 2D.
!
!  Discussion:
!
!    Except for the zero angle case, it should be true that
!
!      ANGLE_DEG_2D(X1,Y1,X2,Y2,X3,Y3)
!    + ANGLE_DEG_2D(X3,Y3,X2,Y2,X1,Y1) = 360.0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X1, Y1, X2, Y2, X3, Y3, define the rays
!    ( X1-X2, Y1-Y2 ) and ( X3-X2, Y3-Y2 ) which in turn define the
!    angle, counterclockwise from ( X1-X2, Y1-Y2 ).
!
!    Output, real ( kind = 8 ) ANGLE_DEG_2D, the angle swept out by the 
!    rays, measured in degrees.  0 <= ANGLE_DEG_2D < 360.  If either 
!    ray has zero length, then ANGLE_DEG_2D is set to 0.
!
  implicit none

  real ( kind = 8 ) angle_deg_2d
  real ( kind = 8 ) angle_rad_2d
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) radians_to_degrees
  real ( kind = 8 ) x
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ) y
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) y3

  x = ( x1 - x2 ) * ( x3 - x2 ) + ( y1 - y2 ) * ( y3 - y2 )
  y = ( x1 - x2 ) * ( y3 - y2 ) - ( y1 - y2 ) * ( x3 - x2 )

  if ( x == 0.0D+00 .and. y == 0.0D+00 ) then

    angle_deg_2d = 0.0D+00

  else

    angle_rad_2d = atan2 ( y, x )

    if ( angle_rad_2d < 0.0D+00 ) then
      angle_rad_2d = angle_rad_2d + 2.0D+00 * pi
    end if

    angle_deg_2d = radians_to_degrees ( angle_rad_2d )

  end if

  return
end
subroutine angle_to_rgb ( angle, r, g, b )

!*****************************************************************************80
!
!! ANGLE_TO_RGB returns a color on the perimeter of the color hexagon.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 December 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ANGLE, the angle in the color hexagon.
!    The sextants are defined by the following points:
!        0 degrees, 1, 0, 0, red;
!       60 degrees, 1, 1, 0, yellow;
!      120 degrees, 0, 1, 0, green;
!      180 degrees, 0, 1, 1, cyan;
!      240 degrees, 0, 0, 1, blue;
!      300 degrees, 1, 0, 1, magenta.
!
!    Output, real ( kind = 8 ) R, G, B, RGB specifications for the color
!    that lies at the given angle, on the perimeter of the color hexagon.  One
!    value will be 1, and one value will be 0.
!
  implicit none

  real ( kind = 8 ) angle
  real ( kind = 8 ) angle2
  real ( kind = 8 ) b
  real ( kind = 8 ) g
  real ( kind = 8 ), parameter :: degrees_to_radians = &
    3.14159265D+00 / 180.0D+00
  real ( kind = 8 ) r

  angle = mod ( angle, 360.0D+00 )

  if ( angle < 0.0D+00 ) then
    angle = angle + 360.0D+00
  end if

  if ( angle <= 60.0D+00 ) then

    angle2 = degrees_to_radians * 3.0D+00 * angle / 4.0D+00
    r = 1.0D+00
    g = tan ( angle2 )
    b = 0.0D+00

  else if ( angle <= 120.0D+00 ) then

    angle2 = degrees_to_radians * 3.0D+00 * angle / 4.0D+00
    r = cos ( angle2 ) / sin ( angle2 )
    g = 1.0D+00
    b = 0.0D+00

  else if ( angle <= 180.0D+00 ) then

    angle2 = degrees_to_radians * 3.0D+00 * ( angle - 120.0D+00 ) / 4.0D+00
    r = 0.0D+00
    g = 1.0D+00
    b = tan ( angle2 )

  else if ( angle <= 240.0D+00 ) then

    angle2 = degrees_to_radians * 3.0D+00 * ( angle - 120.0D+00 ) / 4.0D+00
    r = 0.0D+00
    g = cos ( angle2 ) / sin ( angle2 )
    b = 1.0D+00

  else if ( angle <= 300.0D+00 ) then

    angle2 = degrees_to_radians * 3.0D+00 * ( angle - 240.0D+00 ) / 4.0D+00
    r = tan ( angle2 )
    g = 0.0D+00
    b = 1.0D+00

  else if ( angle <= 360.0D+00 ) then

    angle2 = degrees_to_radians * 3.0D+00 * ( angle - 240.0D+00 ) / 4.0D+00
    r = 1.0D+00
    g = 0.0D+00
    b = cos ( angle2 ) / sin ( angle2 )

  end if

  return
end
subroutine balloon_plot ( box_requested, box_xy, circle_requested, &
  circle_x, circle_y, circle_r, filedot_name, dim_num, &
  node_num, coord, plot_min, plot_max, x, y, shadow, title, radius )

!*****************************************************************************80
!
!! BALLOON_PLOT plots points with maximal nonintersecting circles.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, logical BOX_REQUESTED, is true if the user has specified
!    a box to be drawn.
!
!    Input, real ( kind = 8 ) BOX_XY(2,2), the coordinates of the lower left and
!    upper right corners of the requested box.
!
!    Input, character ( len = * ) FILEDOT_NAME, the name of the dot file.
!
!    Input, integer ( kind = 4 ) DIM_NUM, the number of spatial dimensions.
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of points.
!
!    Input, real ( kind = 8 ) COORD(DIM_NUM,NODE_NUM), the point coordinates.
!
!    Input, real ( kind = 8 ) PLOT_MIN(2), PLOT_MAX(2), the minimum and maximum
!    X and Y values to plot.
!
!    Input, integer ( kind = 4 ) X, Y, the indices of the two dimensions to plot.
!
!    Input, logical SHADOW, is true if the user would like the
!    X and Y axes to be marked at the shadow of each point.
!
!    Input, character ( len = * ) TITLE, a title for the plot.
!
!    Input, real ( kind = 8 ) RADIUS(NODE_NUM), the radius of the circle
!    around each point.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) node_num

  real ( kind = 8 ) b
  real ( kind = 8 ) bitx
  real ( kind = 8 ) bity
  logical box_requested
  real ( kind = 8 ) box_xy(2,2)
  logical circle_requested
  real ( kind = 8 ) circle_r
  real ( kind = 8 ) circle_x
  real ( kind = 8 ) circle_y
  real ( kind = 8 ) coord(dim_num,node_num)
  character ( len = * ) filedot_name
  integer ( kind = 4 ) filedot_unit
  real ( kind = 8 ) g
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  real ( kind = 8 ) plot_max(2)
  real ( kind = 8 ) plot_min(2)
  real ( kind = 8 ) r
  real ( kind = 8 ) radius(node_num)
  logical shadow
  character ( len = * ) title
  integer ( kind = 4 ) x
  real ( kind = 8 ) xvec(4)
  integer ( kind = 4 ) y
  real ( kind = 8 ) yvec(4)

  call get_unit ( filedot_unit )

  call ps_file_open ( filedot_name, filedot_unit, ierror )

  call eps_file_head ( filedot_name, 36, 36, 576, 756  )

  call ps_page_head ( plot_min(1), plot_min(2), plot_max(1), plot_max(2) )
!
!  Initialize the line width to 1.
!
  call ps_line_width ( 1 )
!
!  Print title, if requested.
!
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
!
!  Because filled objects can obscure lines, draw them first.
!
  r = 0.75D+00
  g = 0.75D+00
  b = 1.00D+00

  call ps_color_fill_set ( r, g, b )

  do i = 1, node_num
    call ps_circle_fill ( coord(x,i), coord(y,i), radius(i) )
  end do

  r = 0.0D+00
  g = 0.0D+00
  b = 0.0D+00

  call ps_color_fill_set ( r, g, b )
!
!  Shadow the points, if requested.
!
  if ( shadow ) then

    bitx = 0.025D+00 * ( plot_max ( 1 ) - plot_min ( 1 ) )
    bity = 0.025D+00 * ( plot_max ( 2 ) - plot_min ( 2 ) )

    do i = 1, node_num
      call ps_line ( plot_min(1), coord(y,i), plot_min(1) + bitx, coord(y,i) )
      call ps_line ( coord(x,i), plot_min(2), coord(x,i), plot_min(2) + bity )
    end do

  end if
!
!  Draw a user-specified box, if requested.
!
  if ( box_requested ) then
    call ps_comment ( 'User-requested box:' )
    call ps_line ( box_xy(1,1), box_xy(2,1), box_xy(1,2), box_xy(2,1) )
    call ps_line ( box_xy(1,2), box_xy(2,1), box_xy(1,2), box_xy(2,2) )
    call ps_line ( box_xy(1,2), box_xy(2,2), box_xy(1,1), box_xy(2,2) )
    call ps_line ( box_xy(1,1), box_xy(2,2), box_xy(1,1), box_xy(2,1) )
  end if
!
!  Draw a user-specified circle, if requested.
!
  if ( circle_requested ) then
    call ps_comment ( 'User-requested circle:' )
    call ps_circle ( circle_x, circle_y, circle_r )
  end if

  r = 0.0D+00
  g = 0.0D+00
  b = 0.0D+00

  call ps_color_fill_set ( r, g, b )

  if ( node_num <= 350 ) then
    do i = 1, node_num
      call ps_mark_disk ( coord(x,i), coord(y,i) )
    end do
  else
    do i = 1, node_num
      call ps_mark_point ( coord(x,i), coord(y,i) )
    end do
  end if

  call ps_page_tail ( )

  call eps_file_tail ( )

  call ps_file_close ( filedot_unit )

  return
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
subroutine box_ray_int_2d ( xmin, ymin, xmax, ymax, xa, ya, xb, yb, xi, yi )

!*****************************************************************************80
!
!! BOX_RAY_INT_2D: intersection ( box, ray ) in 2D.
!
!  Discussion:
!
!    The box is assumed to be a rectangle with sides aligned on coordinate
!    axes.
!
!    The origin of the ray is assumed to be inside the box.  This
!    guarantees that the ray will intersect the box in exactly one point.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 June 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XMIN, YMIN, the lower left corner of the box.
!
!    Input, real ( kind = 8 ) XMAX, YMAX, the upper right corner of the box.
!
!    Input, real ( kind = 8 ) XA, YA, the origin of the ray, which should be
!    inside the box.
!
!    Input, real ( kind = 8 ) XB, YB, a second point on the ray.
!
!    Output, real ( kind = 8 ) XI, YI, the point on the box intersected
!    by the ray.
!
  implicit none

  logical inside
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) side
  real ( kind = 8 ) xa
  real ( kind = 8 ) xb
  real ( kind = 8 ) xc
  real ( kind = 8 ) xd
  real ( kind = 8 ) xi
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) ya
  real ( kind = 8 ) yb
  real ( kind = 8 ) yc
  real ( kind = 8 ) yd
  real ( kind = 8 ) yi
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin

  do side = 1, 4

    if ( side == 1 ) then
      xc = xmin
      yc = ymin
      xd = xmax
      yd = ymin
    else if ( side == 2 ) then
      xc = xmax
      yc = ymin
      xd = xmax
      yd = ymax
    else if ( side == 3 ) then
      xc = xmax
      yc = ymax
      xd = xmin
      yd = ymax
    else if ( side == 4 ) then
      xc = xmin
      yc = ymax
      xd = xmin
      yd = ymin
    end if

    call angle_contains_ray_2d ( inside, xc, yc, xa, ya, xd, yd, xb, yb )

    if ( inside ) then
      exit
    end if

    if ( side == 4 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'BOX_RAY_INT_2D - Fatal error!'
      write ( *, '(a)' ) '  No intersection could be found.'
      stop
    end if

  end do

  call lines_exp_int_2d ( xa, ya, xb, yb, xc, yc, xd, yd, ival, xi, yi )

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
function ch_is_digit ( c )

!*****************************************************************************80
!
!! CH_IS_DIGIT returns .TRUE. if a character is a decimal digit.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C, the character to be analyzed.
!
!    Output, logical CH_IS_DIGIT, .TRUE. if C is a digit, .FALSE. otherwise.
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
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
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
!    Input, integer ( kind = 4 ) N, the number of points desired.  N must be at least 1.
!
!    Output, real ( kind = 8 ) X(N), Y(N), the coordinates of points
!    on the circle.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) angle
  integer ( kind = 4 ) i
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
subroutine dash_plot ( box_requested, box_xy, circle_requested, &
  circle_x, circle_y, circle_r, connect, file_name, dim_num, &
  node_num, coord, plot_min, plot_max, x, y, shadow, title )

!*****************************************************************************80
!
!! DASH_PLOT plots a set of points, and connects them.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) CONNECT(NODE_NUM), connection indicator.
!    If CONNECT(I) is:
!    0, it is an isolated point;
!    1, it is connected to the previous point only;
!    2, it is connected to the next point only;
!    3, it is connected to the previous and next points.
!
!    Input, character ( len = * ) FILE_NAME, the name of the dot file.
!
!    Input, integer ( kind = 4 ) DIM_NUM, the number of spatial dimensions.
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of points.
!
!    Input, real ( kind = 8 ) COORD(DIM_NUM,NODE_NUM), the point coordinates.
!
!    Input, real ( kind = 8 ) PLOT_MIN(2), PLOT_MAX(2), the minimum and maximum
!    X and Y values to plot.
!
!    Input, integer ( kind = 4 ) X, Y, the indices of the two dimensions to plot.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) node_num
!
  real ( kind = 8 ) bitx
  real ( kind = 8 ) bity
  logical box_requested
  real ( kind = 8 ) box_xy(2,2)
  logical circle_requested
  real ( kind = 8 ) circle_r
  real ( kind = 8 ) circle_x
  real ( kind = 8 ) circle_y
  integer ( kind = 4 ) connect(node_num)
  real ( kind = 8 ) coord(dim_num,node_num)
  character ( len = * ) file_name
  integer ( kind = 4 ) file_unit
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  real ( kind = 8 ) plot_max(2)
  real ( kind = 8 ) plot_min(2)
  logical shadow
  character ( len = * ) title
  integer ( kind = 4 ) x
  real ( kind = 8 ) xvec(4)
  integer ( kind = 4 ) y
  real ( kind = 8 ) yvec(4)

  call get_unit ( file_unit )

  call ps_file_open ( file_name, file_unit, ierror )

  call eps_file_head ( file_name, 36, 36, 576, 756  )

  call ps_page_head ( plot_min(1), plot_min(2), plot_max(1), plot_max(2) )
!
!  Initialize the line width to 1.
!
  call ps_line_width ( 1 )
!
!  Print title, if requested.
!
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
!
!  Shadow the points, if requested.
!
  if ( shadow ) then

    bitx = 0.025D+00 * ( plot_max ( 1 ) - plot_min ( 1 ) )
    bity = 0.025D+00 * ( plot_max ( 2 ) - plot_min ( 2 ) )

    do i = 1, node_num
      call ps_line ( plot_min(1), coord(y,i), plot_min(1) + bitx, coord(y,i) )
      call ps_line ( coord(x,i), plot_min(2), coord(x,i), plot_min(2) + bity )
    end do

  end if
!
!  Draw a user-specified box, if requested.
!
  if ( box_requested ) then
    call ps_comment ( 'User-requested box:' )
    call ps_line ( box_xy(1,1), box_xy(2,1), box_xy(1,2), box_xy(2,1) )
    call ps_line ( box_xy(1,2), box_xy(2,1), box_xy(1,2), box_xy(2,2) )
    call ps_line ( box_xy(1,2), box_xy(2,2), box_xy(1,1), box_xy(2,2) )
    call ps_line ( box_xy(1,1), box_xy(2,2), box_xy(1,1), box_xy(2,1) )
  end if
!
!  Draw a user-specified circle, if requested.
!
  if ( circle_requested ) then
    call ps_comment ( 'User-requested circle:' )
    call ps_circle ( circle_x, circle_y, circle_r )
  end if
!
!  Draw the points, with big or little marks as appropriate.
!
  if ( node_num <= 350 ) then
    do i = 1, node_num
      call ps_mark_disk ( coord(x,i), coord(y,i) )
    end do
  else
    do i = 1, node_num
      call ps_mark_point ( coord(x,i), coord(y,i) )
    end do
  end if

  do i = 1, node_num-1
    if ( connect(i) == 2 .or. connect(i) == 3 ) then
      call ps_line ( coord(x,i), coord(y,i), coord(x,i+1), coord(y,i+1) )
    end if
  end do

  if ( connect(node_num) == 2 .or. connect(node_num) == 3 ) then
    call ps_line ( coord(x,node_num), coord(y,node_num), coord(x,1), &
      coord(y,1) )
  end if

  call ps_page_tail ( )

  call eps_file_tail ( )

  call ps_file_close ( file_unit )

  return
end
subroutine data_read ( dim_num, file_in_name, node_num )

!*****************************************************************************80
!
!! DATA_READ reads the point data from a file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) DIM_NUM, the number of spatial dimensions.
!
!    Input, character ( len = * ) FILE_IN_NAME, the file to be read.
!
!    Output, integer ( kind = 4 ) NODE_NUM, the number of points.
!
  use AA_point_data_module

  implicit none

  integer ( kind = 4 ) dim_num
  character ( len = * ) file_in_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) node_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DATA_READ'
  write ( *, '(a)' ) '  Read point data from "' // trim ( file_in_name ) // '".'
!
!  Get the "width" of the file.
!
  call file_column_count ( file_in_name, dim_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Points are assumed to have dimension ', dim_num
!
!  Get the "length" of the file.
!
  call points_count ( file_in_name, dim_num, node_num )
!
!  Now allocate some data.
!
  allocate ( connect(node_num) )
  allocate ( coord(dim_num,node_num) )
  allocate ( coord_min(dim_num) )
  allocate ( coord_max(dim_num) )
  allocate ( radius(node_num) )
  allocate ( tag(node_num) )
!
!  Read the coordinates into COORD.
!
  call points_read ( connect, file_in_name, dim_num, node_num, coord )
!
!  Determine the range.
!
  do i = 1, dim_num
    coord_min(i) = minval ( coord(i,1:node_num) )
    coord_max(i) = maxval ( coord(i,1:node_num) )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' Coord  Min  Max'
  write ( *, '(a)' ) ' '
  do i = 1, dim_num
    write ( *, '(i3,2f10.4)' ) i, coord_min(i), coord_max(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DATA_READ'
  write ( *, '(a)' ) '  The data has been read.'

  return
end
subroutine delaunay_plot ( box_requested, box_xy, circle_requested, &
  circle_x, circle_y, circle_r, file_name, dim_num, &
  node_num, coord, plot_min, plot_max, x, y, shadow, title )

!*****************************************************************************80
!
!! DELAUNAY_PLOT plots the Delaunay triangulation of a pointset.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 January 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, logical BOX_REQUESTED, is true if the user has specified
!    a box to be drawn.
!
!    Input, real ( kind = 8 ) BOX_XY(2,2), the coordinates of the lower left and
!    upper right corners of the requested box.
!
!    Input, character ( len = * ) FILE_NAME, the name of the output file.
!
!    Input, integer ( kind = 4 ) DIM_NUM, the number of spatial dimensions.
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of points.
!
!    Input, real ( kind = 8 ) COORD(DIM_NUM,NODE_NUM), the point coordinates.
!
!    Input, real ( kind = 8 ) PLOT_MIN(2), PLOT_MAX(2), the minimum and maximum
!    X and Y values to plot.
!
!    Input, integer ( kind = 4 ) X, Y, the indices of the two dimensions to plot.
!
!    Input, logical SHADOW, is true if the user would like the
!    X and Y axes to be marked at the shadow of each point.
!
!    Input, character ( len = * ) TITLE, a title for the plot.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) node_num

  real ( kind = 8 ) bitx
  real ( kind = 8 ) bity
  logical box_requested
  real ( kind = 8 ) box_xy(2,2)
  logical circle_requested
  real ( kind = 8 ) circle_r
  real ( kind = 8 ) circle_x
  real ( kind = 8 ) circle_y
  real ( kind = 8 ) coord(dim_num,node_num)
  real ( kind = 8 ) coord2(2,node_num)
  logical, parameter :: debug = .false.
  character ( len = * ) file_name
  integer ( kind = 4 ) file_unit
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j3
  integer ( kind = 4 ) triangle_node(3,2*node_num)
  real ( kind = 8 ) plot_max(2)
  real ( kind = 8 ) plot_min(2)
  logical shadow
  character ( len = * ) title
  integer ( kind = 4 ) triangle_neighbor(3,2*node_num)
  integer ( kind = 4 ) triangle_num
  integer ( kind = 4 ) x
  real ( kind = 8 ) xvec(4)
  real ( kind = 8 ) xx
  integer ( kind = 4 ) y
  real ( kind = 8 ) yvec(4)
  real ( kind = 8 ) yy
!
!  Compute the Delaunay triangulation.
!
  coord2(1,1:node_num) = coord(x,1:node_num)
  coord2(2,1:node_num) = coord(y,1:node_num)

  call dtris2 ( node_num, coord2, triangle_num, triangle_node, triangle_neighbor )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DELAUNAY_PLOT:'
  write ( *, '(a,i6)' ) '  The number of triangles is ', triangle_num
!
!  Print the triangulation.
!
  if ( debug ) then
    call triangulation_order3_print ( node_num, triangle_num, coord2, &
      triangle_node, triangle_neighbor )
  end if
!
!  Plot the Delaunay triangulation.
!
  call get_unit ( file_unit )

  call ps_file_open ( file_name, file_unit, ierror )

  call eps_file_head ( file_name, 36, 36, 576, 756  )

  call ps_page_head ( plot_min(1), plot_min(2), plot_max(1), plot_max(2) )
!
!  Initialize the line width to 1.
!
  call ps_line_width ( 1 )
!
!  Print title, if requested.
!
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
!
!  Shadow the points, if requested.
!
  if ( shadow ) then

    bitx = 0.025D+00 * ( plot_max ( 1 ) - plot_min ( 1 ) )
    bity = 0.025D+00 * ( plot_max ( 2 ) - plot_min ( 2 ) )

    do i = 1, node_num
      call ps_line ( plot_min(1), coord(y,i), plot_min(1) + bitx, coord(y,i) )
      call ps_line ( coord(x,i), plot_min(2), coord(x,i), plot_min(2) + bity )
    end do

  end if
!
!  Draw a user-specified box, if requested.
!
  if ( box_requested ) then
    call ps_comment ( 'User-requested box:' )
    call ps_line ( box_xy(1,1), box_xy(2,1), box_xy(1,2), box_xy(2,1) )
    call ps_line ( box_xy(1,2), box_xy(2,1), box_xy(1,2), box_xy(2,2) )
    call ps_line ( box_xy(1,2), box_xy(2,2), box_xy(1,1), box_xy(2,2) )
    call ps_line ( box_xy(1,1), box_xy(2,2), box_xy(1,1), box_xy(2,1) )
  end if
!
!  Draw a user-specified circle, if requested.
!
  if ( circle_requested ) then
    call ps_comment ( 'User-requested circle:' )
    call ps_circle ( circle_x, circle_y, circle_r )
  end if

  xx = plot_min(1)
  yy = plot_min(2)
  call ps_moveto ( xx, yy )

! call ps_label ( 'Delaunay triangulation' )

  if ( node_num <= 350 ) then
    do i = 1, node_num
      call ps_mark_disk ( coord(x,i), coord(y,i) )
    end do
  else
    do i = 1, node_num
      call ps_mark_point ( coord(x,i), coord(y,i) )
    end do
  end if
!
!  Increase the line width for small problems.
!
  if ( triangle_num < 50 ) then
    call ps_line_width ( 3 )
  else if ( triangle_num < 250 ) then
    call ps_line_width ( 2 )
  else
    call ps_line_width ( 1 )
  end if

  do i = 1, triangle_num

    j1 = triangle_node(1,i)
    j2 = triangle_node(2,i)
    j3 = triangle_node(3,i)

    call ps_triangle ( coord2(1,j1), coord2(2,j1), coord2(1,j2), &
      coord2(2,j2), coord2(1,j3), coord2(2,j3) )

  end do

  call ps_page_tail ( )

  call eps_file_tail ( )

  call ps_file_close ( file_unit )

  return
end
function diaedg ( x0, y0, x1, y1, x2, y2, x3, y3 )

!*****************************************************************************80
!
!! DIAEDG chooses a diagonal edge.
!
!  Discussion:
!
!    The routine determines whether 0--2 or 1--3 is the diagonal edge
!    that should be chosen, based on the circumcircle criterion, where
!    (X0,Y0), (X1,Y1), (X2,Y2), (X3,Y3) are the vertices of a simple
!    quadrilateral in counterclockwise order.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 February 2001
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!      using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X0, Y0, X1, Y1, X2, Y2, X3, Y3, the
!    coordinates of the vertices of a quadrilateral, given in
!    counter clockwise order.
!
!    Output, integer ( kind = 4 ) DIAEDG, chooses a diagonal:
!    +1, if diagonal edge 02 is chosen;
!    -1, if diagonal edge 13 is chosen;
!     0, if the four vertices are cocircular.
!
  implicit none

  real ( kind = 8 ) ca
  real ( kind = 8 ) cb
  integer ( kind = 4 ) diaedg
  real ( kind = 8 ) dx10
  real ( kind = 8 ) dx12
  real ( kind = 8 ) dx30
  real ( kind = 8 ) dx32
  real ( kind = 8 ) dy10
  real ( kind = 8 ) dy12
  real ( kind = 8 ) dy30
  real ( kind = 8 ) dy32
  real ( kind = 8 ) s
  real ( kind = 8 ) tol
  real ( kind = 8 ) tola
  real ( kind = 8 ) tolb
  real ( kind = 8 ) x0
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ) y0
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) y3

  tol = 100.0D+00 * epsilon ( tol )

  dx10 = x1 - x0
  dy10 = y1 - y0
  dx12 = x1 - x2
  dy12 = y1 - y2
  dx30 = x3 - x0
  dy30 = y3 - y0
  dx32 = x3 - x2
  dy32 = y3 - y2

  tola = tol * max ( abs ( dx10 ), abs ( dy10 ), abs ( dx30 ), abs ( dy30 ) )
  tolb = tol * max ( abs ( dx12 ), abs ( dy12 ), abs ( dx32 ), abs ( dy32 ) )

  ca = dx10 * dx30 + dy10 * dy30
  cb = dx12 * dx32 + dy12 * dy32

  if ( tola < ca .and. tolb < cb ) then

    diaedg = -1

  else if ( ca < -tola .and. cb < -tolb ) then

    diaedg = 1

  else

    tola = max ( tola, tolb )
    s = ( dx10 * dy30 - dx30 * dy10 ) * cb + ( dx32 * dy12 - dx12 * dy32 ) * ca

    if ( tola < s ) then
      diaedg = -1
    else if ( s < -tola ) then
      diaedg = 1
    else
      diaedg = 0
    end if

  end if

  return
end
subroutine digit_inc ( c )

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
!    Input/output, character C, a digit to be incremented.
!
  implicit none

  character c
  integer ( kind = 4 ) digit

  call ch_to_digit ( c, digit )

  if ( digit == -1 ) then
    return
  end if

  digit = digit + 1

  if ( digit == 10 ) then
    digit = 0
  end if

  call digit_to_ch ( digit, c )

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
subroutine dot_plot ( box_requested, box_xy, circle_requested, circle_x, &
  circle_y, circle_r, file_name, dim_num, &
  node_num, node_xy, plot_min, plot_max, x, y, shadow, title, marker_size )

!*****************************************************************************80
!
!! DOT_PLOT makes a plot of a set of points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 July 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, logical BOX_REQUESTED, is true if the user has specified
!    a box to be drawn.
!
!    Input, real ( kind = 8 ) BOX_XY(2,2), the coordinates of the lower left and
!    upper right corners of the requested box.
!
!    Input, character ( len = * ) FILE_NAME, the name of the dot file.
!
!    Input, integer ( kind = 4 ) DIM_NUM, the number of spatial dimensions.
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of points.
!
!    Input, real ( kind = 8 ) NODE_XY(DIM_NUM,NODE_NUM), the point coordinates.
!
!    Input, real ( kind = 8 ) PLOT_MIN(2), PLOT_MAX(2), the minimum and maximum
!    X and Y values to plot.
!
!    Input, integer ( kind = 4 ) X, Y, the indices of the two dimensions to plot.
!
!    Input, logical SHADOW, is true if the user would like the
!    X and Y axes to be marked at the shadow of each point.
!
!    Input, character ( len = * ) TITLE, a title for the plot.
!
!    Input, integer ( kind = 4 ) MARKER_SIZE, controls the size of the dots.
!    A value of 3, 4, or 5 is typical.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) node_num

  real ( kind = 8 ) b
  real ( kind = 8 ) bitx
  real ( kind = 8 ) bity
  logical box_requested
  real ( kind = 8 ) box_xy(2,2)
  logical circle_requested
  real ( kind = 8 ) circle_r
  real ( kind = 8 ) circle_x
  real ( kind = 8 ) circle_y
  character ( len = 40 ) date_time
  integer ( kind = 4 ) delta
  character ( len = * ) file_name
  integer ( kind = 4 ) file_unit
  real ( kind = 8 ) g
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) marker_size
  real ( kind = 8 ) node_xy(dim_num,node_num)
  real ( kind = 8 ) plot_max(2)
  real ( kind = 8 ) plot_min(2)
  real ( kind = 8 ) r
  logical shadow
  character ( len = * ) title
  integer ( kind = 4 ) x
  real ( kind = 8 ) x_max
  real ( kind = 8 ) x_min
  integer ( kind = 4 ) x_ps
  integer ( kind = 4 ) :: x_ps_max = 576
  integer ( kind = 4 ) :: x_ps_max_clip = 594
  integer ( kind = 4 ) :: x_ps_min = 36
  integer ( kind = 4 ) :: x_ps_min_clip = 18
  real ( kind = 8 ) x_scale
  real ( kind = 8 ) xvec(4)
  integer ( kind = 4 ) y
  real ( kind = 8 ) y_max
  real ( kind = 8 ) y_min
  integer ( kind = 4 ) y_ps
  integer ( kind = 4 ) :: y_ps_max = 666
  integer ( kind = 4 ) :: y_ps_max_clip = 684
  integer ( kind = 4 ) :: y_ps_min = 126
  integer ( kind = 4 ) :: y_ps_min_clip = 108
  real ( kind = 8 ) y_scale
  real ( kind = 8 ) yvec(4)

  call timestring ( date_time )
!
!  We need to do some figuring here, so that we can determine
!  the range of the data, and hence the height and width
!  of the piece of paper.
!
  x_max = maxval ( node_xy(x,1:node_num) )
  x_min = minval ( node_xy(x,1:node_num) )
  if ( box_requested ) then
    x_max = max ( x_max, box_xy(1,1) )
    x_max = max ( x_max, box_xy(1,2) )
    x_min = min ( x_min, box_xy(1,1) )
    x_min = min ( x_min, box_xy(1,2) )
  end if
  x_scale = x_max - x_min

  x_max = x_max + 0.05D+00 * x_scale
  x_min = x_min - 0.05D+00 * x_scale
  x_scale = x_max - x_min

  y_max = maxval ( node_xy(y,1:node_num) )
  y_min = minval ( node_xy(y,1:node_num) )
  if ( box_requested ) then
    y_max = max ( y_max, box_xy(2,1) )
    y_max = max ( y_max, box_xy(2,2) )
    y_min = min ( y_min, box_xy(2,1) )
    y_min = min ( y_min, box_xy(2,2) )
  end if
  y_scale = y_max - y_min

  y_max = y_max + 0.05D+00 * y_scale
  y_min = y_min - 0.05D+00 * y_scale
  y_scale = y_max - y_min

  if ( x_scale < y_scale ) then

    delta = nint ( real ( x_ps_max - x_ps_min, kind = 8 ) &
      * ( y_scale - x_scale ) / ( 2.0D+00 * y_scale ) )

    x_ps_max = x_ps_max - delta
    x_ps_min = x_ps_min + delta

    x_ps_max_clip = x_ps_max_clip - delta
    x_ps_min_clip = x_ps_min_clip + delta

    x_scale = y_scale

  else if ( y_scale < x_scale ) then

    delta = nint ( real ( y_ps_max - y_ps_min, kind = 8 ) &
      * ( x_scale - y_scale ) / ( 2.0D+00 * x_scale ) )

    y_ps_max = y_ps_max - delta
    y_ps_min = y_ps_min + delta

    y_ps_max_clip = y_ps_max_clip - delta
    y_ps_min_clip = y_ps_min_clip + delta

    y_scale = x_scale

  end if

  call get_unit ( file_unit )

  call ps_file_open ( file_name, file_unit, ierror )

  call eps_file_head ( file_name, x_ps_min, y_ps_min, x_ps_max, y_ps_max )

  call ps_page_head ( x_min, y_min, x_max, y_max )
!
!  Initialize the line width to 1.
!
  call ps_line_width ( 1 )
!
!  Print title, if requested.
!
  if ( 0 < len_trim ( title ) ) then
    call ps_font_size ( 0.30D+00 )
    call ps_moveto ( plot_min(1), plot_max(2) )
    call ps_label ( title )
  end if
!
!  Define a PostScript clipping box.
!
  xvec(1) = x_min
  xvec(2) = x_max
  xvec(3) = x_max
  xvec(4) = x_min

  yvec(1) = y_min
  yvec(2) = y_min
  yvec(3) = y_max
  yvec(4) = y_max

  call ps_clip ( 4, xvec, yvec )
!
!  Shadow the points, if requested.
!
  if ( shadow ) then

    bitx = 0.025D+00 * ( plot_max ( 1 ) - plot_min ( 1 ) )
    bity = 0.025D+00 * ( plot_max ( 2 ) - plot_min ( 2 ) )

    do i = 1, node_num
      call ps_line ( plot_min(1), node_xy(y,i), plot_min(1) + bitx, &
        node_xy(y,i) )
      call ps_line ( node_xy(x,i), plot_min(2), node_xy(x,i), &
        plot_min(2) + bity )
    end do

  end if
!
!  Draw a user-specified box, if requested.
!
  if ( box_requested ) then
    call ps_comment ( 'User-requested box:' )
    call ps_line ( box_xy(1,1), box_xy(2,1), box_xy(1,2), box_xy(2,1) )
    call ps_line ( box_xy(1,2), box_xy(2,1), box_xy(1,2), box_xy(2,2) )
    call ps_line ( box_xy(1,2), box_xy(2,2), box_xy(1,1), box_xy(2,2) )
    call ps_line ( box_xy(1,1), box_xy(2,2), box_xy(1,1), box_xy(2,1) )
  end if
!
!  Draw a user-specified circle, if requested.
!
  if ( circle_requested ) then
    call ps_comment ( 'User-requested circle:' )
    call ps_circle ( circle_x, circle_y, circle_r )
  end if

  call ps_comment ( 'Set color to bluish green:' )

  call ps_setting_int ( 'SET', 'MARKER_SIZE', marker_size )

  r = 0.000D+00
  g = 0.750D+00
  b = 0.150D+00

  call ps_color_fill_set ( r, g, b )

  if ( node_num <= 1100 ) then

    do i = 1, node_num
      call ps_mark_disk ( node_xy(x,i), node_xy(y,i) )
    end do

  else
    do i = 1, node_num
      call ps_mark_point ( node_xy(x,i), node_xy(y,i) )
    end do
  end if

  call ps_page_tail ( )

  call eps_file_tail ( )

  call ps_file_close ( file_unit )

  return
end
subroutine dtris2 ( node_num, node_xy, triangle_num, triangle_node, &
  triangle_neighbor )

!*****************************************************************************80
!
!! DTRIS2 constructs a Delaunay triangulation of 2D vertices.
!
!  Discussion:
!
!    The routine constructs the Delaunay triangulation of a set of 2D vertices
!    using an incremental approach and diagonal edge swaps.  Vertices are
!    first sorted in lexicographically increasing (X,Y) order, and
!    then are inserted one at a time from outside the convex hull.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 August 2001
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input/output, real ( kind = 8 ) NODE_XY(2,NODE_NUM), the coordinates
!    of the nodes.  On output, the vertices have been sorted into 
!    dictionary order.
!
!    Output, integer ( kind = 4 ) TRIANGLE_NUM, the number of triangles in the triangulation;
!    TRIANGLE_NUM is equal to 2*NODE_NUM - NB - 2, where NB is the number
!    of boundary vertices.
!
!    Output, integer ( kind = 4 ) TRIANGLE_NODE(3,TRIANGLE_NUM), the nodes that make up each
!    triangle.  The elements are indices of P.  The vertices of the triangles
!    are in counter clockwise order.
!
!    Output, integer ( kind = 4 ) TRIANGLE_NEIGHBOR(3,TRIANGLE_NUM), the triangle neighbor
!    list.  Positive elements are indices of TIL; negative elements are used
!    for links of a counter clockwise linked list of boundary edges; 
!    LINK = -(3*I + J-1) where I, J = triangle, edge index;
!    TRIANGLE_NEIGHBOR(J,I) refers to the neighbor along edge from vertex J 
!    to J+1 (mod 3).
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ) node_num

  real ( kind = 8 ) cmax
  integer ( kind = 4 ) e
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) indx(node_num)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) ledg
  integer ( kind = 4 ) lr
  integer ( kind = 4 ) lrline
  integer ( kind = 4 ) ltri
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) n
  real ( kind = 8 ) node_xy(dim_num,node_num)
  integer ( kind = 4 ) redg
  integer ( kind = 4 ) rtri
  integer ( kind = 4 ) stack(node_num)
  integer ( kind = 4 ) t
  real ( kind = 8 ) tol
  integer ( kind = 4 ) top
  integer ( kind = 4 ) triangle_neighbor(3,node_num*2)
  integer ( kind = 4 ) triangle_num
  integer ( kind = 4 ) triangle_node(3,node_num*2)

  tol = 100.0D+00 * epsilon ( tol )

  ierr = 0
!
!  Sort the vertices by increasing (x,y).
!
  call r82vec_sort_heap_index_a ( node_num, node_xy, indx )

  call r82vec_permute ( node_num, node_xy, indx )
!
!  Make sure that the data nodes are "reasonably" distinct.
!
  m1 = 1

  do i = 2, node_num

    m = m1
    m1 = i

    k = 0

    do j = 1, dim_num

      cmax = max ( abs ( node_xy(j,m) ), abs ( node_xy(j,m1) ) )

      if ( tol * ( cmax + 1.0D+00 ) &
           < abs ( node_xy(j,m) - node_xy(j,m1) ) ) then
        k = j
        exit
      end if

    end do

    if ( k == 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DTRIS2 - Fatal error!'
      write ( *, '(a,i6)' ) '  Fails for point number I = ', i
      write ( *, '(a,i6)' ) '  M = ', m
      write ( *, '(a,i6)' ) '  M1 = ', m1
      write ( *, '(a,2g14.6)' ) '  NODE_XY(M)  = ', node_xy(1:dim_num,m)
      write ( *, '(a,2g14.6)' ) '  NODE_XY(M1) = ', node_xy(1:dim_num,m1)
      ierr = 224
      stop
    end if

  end do
!
!  Starting from nodes M1 and M2, search for a third point M that
!  makes a "healthy" triangle (M1,M2,M)
!
  m1 = 1
  m2 = 2
  j = 3

  do

    if ( node_num < j ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DTRIS2 - Fatal error!'
      ierr = 225
      stop
    end if

    m = j

    lr = lrline ( node_xy(1,m), node_xy(2,m), node_xy(1,m1), &
      node_xy(2,m1), node_xy(1,m2), node_xy(2,m2), 0.0D+00 )

    if ( lr /= 0 ) then
      exit
    end if

    j = j + 1

  end do
!
!  Set up the triangle information for (M1,M2,M), and for any other
!  triangles you created because points were collinear with M1, M2.
!
  triangle_num = j - 2

  if ( lr == -1 ) then

    triangle_node(1,1) = m1
    triangle_node(2,1) = m2
    triangle_node(3,1) = m
    triangle_neighbor(3,1) = -3

    do i = 2, triangle_num

      m1 = m2
      m2 = i+1

      triangle_node(1,i) = m1
      triangle_node(2,i) = m2
      triangle_node(3,i) = m

      triangle_neighbor(1,i-1) = -3 * i
      triangle_neighbor(2,i-1) = i
      triangle_neighbor(3,i) = i - 1

    end do

    triangle_neighbor(1,triangle_num) = -3 * triangle_num - 1
    triangle_neighbor(2,triangle_num) = -5
    ledg = 2
    ltri = triangle_num

  else

    triangle_node(1,1) = m2
    triangle_node(2,1) = m1
    triangle_node(3,1) = m

    triangle_neighbor(1,1) = -4

    do i = 2, triangle_num

      m1 = m2
      m2 = i+1

      triangle_node(1,i) = m2
      triangle_node(2,i) = m1
      triangle_node(3,i) = m

      triangle_neighbor(3,i-1) = i
      triangle_neighbor(1,i) = -3 * i - 3
      triangle_neighbor(2,i) = i - 1

    end do

    triangle_neighbor(3,triangle_num) = -3 * triangle_num
    triangle_neighbor(2,1) = -3 * triangle_num - 2
    ledg = 2
    ltri = 1

  end if
!
!  Insert the vertices one at a time from outside the convex hull,
!  determine visible boundary edges, and apply diagonal edge swaps until
!  Delaunay triangulation of vertices (so far) is obtained.
!
  top = 0

  do i = j+1, node_num

    m = i
    m1 = triangle_node(ledg,ltri)

    if ( ledg <= 2 ) then
      m2 = triangle_node(ledg+1,ltri)
    else
      m2 = triangle_node(1,ltri)
    end if

    lr = lrline ( node_xy(1,m), node_xy(2,m), node_xy(1,m1), &
      node_xy(2,m1), node_xy(1,m2), node_xy(2,m2), 0.0D+00 )

    if ( 0 < lr ) then
      rtri = ltri
      redg = ledg
      ltri = 0
    else
      l = -triangle_neighbor(ledg,ltri)
      rtri = l / 3
      redg = mod(l,3) + 1
    end if

    call vbedg ( node_xy(1,m), node_xy(2,m), node_num, node_xy, triangle_num, &
      triangle_node, triangle_neighbor, ltri, ledg, rtri, redg )

    n = triangle_num + 1
    l = -triangle_neighbor(ledg,ltri)

    do

      t = l / 3
      e = mod ( l, 3 ) + 1
      l = -triangle_neighbor(e,t)
      m2 = triangle_node(e,t)

      if ( e <= 2 ) then
        m1 = triangle_node(e+1,t)
      else
        m1 = triangle_node(1,t)
      end if

      triangle_num = triangle_num + 1
      triangle_neighbor(e,t) = triangle_num

      triangle_node(1,triangle_num) = m1
      triangle_node(2,triangle_num) = m2
      triangle_node(3,triangle_num) = m

      triangle_neighbor(1,triangle_num) = t
      triangle_neighbor(2,triangle_num) = triangle_num - 1
      triangle_neighbor(3,triangle_num) = triangle_num + 1

      top = top + 1

      if ( node_num < top ) then
        ierr = 8
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'DTRIS2 - Fatal error!'
        write ( *, '(a)' ) '  Stack overflow.'
        stop
      end if

      stack(top) = triangle_num

      if ( t == rtri .and. e == redg ) then
        exit
      end if

    end do

    triangle_neighbor(ledg,ltri) = -3 * n - 1
    triangle_neighbor(2,n) = -3 * triangle_num - 2
    triangle_neighbor(3,triangle_num) = -l

    ltri = n
    ledg = 2

    call swapec ( m, top, ltri, ledg, node_num, node_xy, triangle_num, &
      triangle_node, triangle_neighbor, stack, ierr )

    if ( ierr /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DTRIS2 - Fatal error!'
      write ( *, '(a)' ) '  Error return from SWAPEC.'
      stop
    end if

  end do
!
!  Now account for the sorting that we did.
!
  do i = 1, 3
    do j = 1, triangle_num
      triangle_node(i,j) = indx ( triangle_node(i,j) )
    end do
  end do

  call perm_inv ( node_num, indx )

  call r82vec_permute ( node_num, node_xy, indx )

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
!    Henry McGilton, Mary Campione,
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
subroutine file_column_count ( file_name, ncolumn )

!*****************************************************************************80
!
!! FILE_COLUMN_COUNT counts the number of columns in the first line of a file.
!
!  Discussion:
!
!    The file is assumed to be a simple text file.
!
!    Most lines of the file is presumed to consist of NCOLUMN words, separated
!    by spaces.  There may also be some blank lines, and some comment lines,
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
!    Input, character ( len = * ) FILE_NAME, the name of the file.
!
!    Output, integer ( kind = 4 ) NCOLUMN, the number of columns assumed to be in the file.
!
  implicit none

  character ( len = * ) file_name
  logical got_one
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  character ( len = 256 ) line
  integer ( kind = 4 ) ncolumn
!
!  Open the file.
!
  call get_unit ( iunit )

  open ( unit = iunit, file = file_name, status = 'old', form = 'formatted', &
    access = 'sequential', iostat = ios )

  if ( ios /= 0 ) then
    ncolumn = - 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_COLUMN_COUNT - Fatal error!'
    write ( *, '(a)' ) '  Could not open the file:'
    write ( *, '(a)' ) '    ' // trim ( file_name )
    return
  end if
!
!  Read one line, but skip blank lines and comment lines.
!
  got_one = .false.

  do

    read ( iunit, '(a)', iostat = ios ) line

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

    rewind ( iunit )

    do 

      read ( iunit, '(a)', iostat = ios ) line

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

  close ( unit = iunit )

  if ( .not. got_one ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_COLUMN_COUNT - Warning!'
    write ( *, '(a)' ) '  The file does not seem to contain any data.'
    ncolumn = 0
    return
  end if

  call s_word_count ( line, ncolumn )

  return
end
subroutine file_column_range ( file_name, ncolumn, col_min, col_max )

!*****************************************************************************80
!
!! FILE_COLUMN_RANGE determines the minimum and maximum ranges of each column.
!
!  Discussion:
!
!    The file is assumed to be a simple text file.
!
!    Each line of the file is presumed to consist of NCOLUMN words, separated
!    by spaces.
!
!    The routine computes the range of each column.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the file.
!
!    Input, integer ( kind = 4 ) NCOLUMN, the number of columns assumed to be in the file.
!
!    Output, real ( kind = 8 ) COL_MIN(NCOLUM), COL_MAX(NCOLUMN), the 
!    minimum and maximum for each column.
!
  implicit none

  integer ( kind = 4 ) ncolumn

  real ( kind = 8 ) col_max(ncolumn)
  real ( kind = 8 ) col_min(ncolumn)
  character ( len = * ) file_name
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) j
  character ( len = 256 ) line
  integer ( kind = 4 ) nrow
  real ( kind = 8 ) x(ncolumn)
!
!  Open the file.
!
  call get_unit ( iunit )

  open ( unit = iunit, file = file_name, status = 'old', form = 'formatted', &
    access = 'sequential', iostat = ios )

  if ( ios /= 0 ) then
    ncolumn = - 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_COLUMN_RANGE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the file:'
    write ( *, '(a)' ) '    ' // trim ( file_name )
    return
  end if

  nrow = 0
  col_min(1:ncolumn) = 0.0D+00
  col_max(1:ncolumn) = 0.0D+00

  do

    read ( iunit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      exit
    end if

    call s_to_r8vec ( line, ncolumn, x, ierror )

    if ( ierror /= 0 ) then
      exit
    end if

    nrow = nrow + 1

    if ( nrow == 1 ) then
      col_min(1:ncolumn) = x(1:ncolumn)
      col_max(1:ncolumn) = x(1:ncolumn)
    else
      do j = 1, ncolumn
        col_min(j) = min ( col_min(j), x(j) )
        col_max(j) = max ( col_max(j), x(j) )
      end do
    end if

  end do

  close ( unit = iunit )

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
!  Example:
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
subroutine file_name_inc ( file_name )

!*****************************************************************************80
!
!! FILE_NAME_INC generates the next filename in a series.
!
!  Discussion:
!
!    It is assumed that the digits in the name, whether scattered or
!    connected, represent a number that is to be increased by 1 on
!    each call.  If this number is all 9's on input, the output number
!    is all 0's.  Non-numeric letters of the name are unaffected, and
!    if the name contains no digits, then nothing is done.
!
!  Example:
!
!      Input          Output
!      -----          ------
!      a7to11.txt     a7to12.txt
!      a7to99.txt     a8to00.txt
!      a9to99.txt     a0to00.txt
!      cat.txt        cat.txt
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
!    Input/output, character ( len = * ) FILE_NAME.
!    On input, a character string to be incremented.
!    On output, the incremented string.
!
  implicit none

  character c
  logical ch_is_digit
  character ( len = * ) file_name
  integer ( kind = 4 ) i
  integer ( kind = 4 ) lens

  lens = len_trim ( file_name )

  do i = lens, 1, -1

    c = file_name(i:i)

    if ( ch_is_digit ( c ) ) then

      call digit_inc ( c )

      file_name(i:i) = c

      if ( c /= '0' ) then
        return
      end if

    end if

  end do

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
function i4_log_2 ( i )

!*****************************************************************************80
!
!! I4_LOG_2 returns the integer part of the logarithm base 2 of |I|.
!
!  Example:
!
!     I  I4_LOG_2
!
!     0  -1
!     1,  0
!     2,  1
!     3,  1
!     4,  2
!     5,  2
!     6,  2
!     7,  2
!     8,  3
!     9,  3
!    10,  3
!   127,  6
!   128,  7
!   129,  7
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 January 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the number whose logarithm base 2 is desired.
!
!    Output, integer ( kind = 4 ) I4_LOG_2, the integer part of the logarithm base 2 of
!    the absolute value of I.
!    For positive I4_LOG_2(I), it should be true that
!      2**I4_LOG_2(X) <= |I| < 2**(I4_LOG_2(I)+1).
!    The special case of I4_LOG_2(0) returns -HUGE().
!
  implicit none

  integer ( kind = 4 ) i4_log_2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_abs

  if ( i == 0 ) then

    i4_log_2 = - huge ( i4_log_2 )

  else

    i4_log_2 = 0

    i_abs = abs ( i )

    do while ( 2 <= i_abs )
      i_abs = i_abs / 2
      i4_log_2 = i4_log_2 + 1
    end do

  end if

  return
end
function i4_modp ( i, j )

!*****************************************************************************80
!
!! I4_MODP returns the nonnegative remainder of integer division.
!
!  Discussion:
!
!    If
!      NREM = I4_MODP ( I, J )
!      NMULT = ( I - NREM ) / J
!    then
!      I = J * NMULT + NREM
!    where NREM is always nonnegative.
!
!    The MOD function computes a result with the same sign as the
!    quantity being divided.  Thus, suppose you had an angle A,
!    and you wanted to ensure that it was between 0 and 360.
!    Then mod(A,360) would do, if A was positive, but if A
!    was negative, your result would be between -360 and 0.
!
!    On the other hand, I4_MODP(A,360) is between 0 and 360, always.
!
!  Example:
!
!        I     J     MOD  I4_MODP    Factorization
!
!      107    50       7       7    107 =  2 *  50 + 7
!      107   -50       7       7    107 = -2 * -50 + 7
!     -107    50      -7      43   -107 = -3 *  50 + 43
!     -107   -50      -7      43   -107 =  3 * -50 + 43
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the number to be divided.
!
!    Input, integer ( kind = 4 ) J, the number that divides I.
!
!    Output, integer ( kind = 4 ) I4_MODP, the nonnegative remainder when I is
!    divided by J.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ) j

  if ( j == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_MODP - Fatal error!'
    write ( *, '(a,i6)' ) '  I4_MODP ( I, J ) called with J = ', j
    stop
  end if

  i4_modp = mod ( i, j )

  if ( i4_modp < 0 ) then
    i4_modp = i4_modp + abs ( j )
  end if

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
subroutine i4_to_angle ( i, angle )

!*****************************************************************************80
!
!! I4_TO_ANGLE maps integers to points on a circle.
!
!  Discussion:
!
!    The angles are intended to be used to select colors on a color
!    hexagon whose 6 vertices are red, yellow, green, cyan, blue,
!    magenta.
!
!  Example:
!
!     I   X      ANGLE
!
!     0   0/3      0
!     1   1/3    120
!     2   2/3    240
!
!     3   1/6     60
!     4   3/6    180
!     5   5/6    300
!
!     6   1/12    30
!     7   3/12    90
!     8   5/12   150
!     9   7/12   210
!    10   9/12   270
!    11  11/12   330
!
!    12   1/24    15
!    13   3/24    45
!    14   5/24    75
!    etc
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
!    Input, integer ( kind = 4 ) I, the index of the desired color.
!
!    Output, real ( kind = 8 ) ANGLE, an angle, measured in degrees,
!    between 0 and 360.
!
  implicit none

  real ( kind = 8 ) angle
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_log_2
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) i4

  if ( 0 <= abs ( i ) .and. abs ( i ) <= 2 ) then

    angle = 120.0D+00 * real ( abs ( i ), kind = 8 )

  else

    i1 = i4_log_2 ( abs ( i ) / 3 )
    i2 = abs ( i ) + 1 - 3 * 2**i1
    i3 = 2 * ( i2 - 1 ) + 1
    i4 = 3 * 2**( i1 + 1 )

    angle = 360.0D+00 * real ( i3, kind = 8 ) / real ( i4, kind = 8 )

  end if

  return
end
subroutine i4_to_rgb ( i, r, g, b )

!*****************************************************************************80
!
!! I4_TO_RGB maps integers to RGB colors.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 January 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the index of the desired color.
!
!    Output, integer ( kind = 4 ) R, G, B, the RGB specifications for a color.
!
  implicit none

  real ( kind = 8 ) angle
  integer ( kind = 4 ) b
  integer ( kind = 4 ) g
  integer ( kind = 4 ) i
  integer ( kind = 4 ) r
  real ( kind = 8 ) rb
  real ( kind = 8 ) rg
  real ( kind = 8 ) rr
!
!  Red
!
  if ( i == 1 ) then
    r = 255
    g = 0
    b = 0
!
!  Green
!
  else if ( i == 2 ) then
    r = 0
    g = 255
    b = 0
!
!  Blue
!
  else if ( i == 3 ) then
    r = 0
    g = 0
    b = 255
!
!  Cyan
!
  else if ( i == 4 ) then
    r = 0
    g = 255
    b = 255
!
!  Magenta
!
  else if ( i == 5 ) then
    r = 255
    g = 0
    b = 255
!
!  Yellow
!
  else if ( i == 6 ) then
    r = 255
    g = 255
    b = 0
!
!  Brown5
!
  else if ( i == 7 ) then
    r = 139
    g =  35
    b =  35
!
!  Orange
!
  else if ( i == 8 ) then
    r = 255
    g = 165
    b = 0
!
!  Goldenrod5
!
  else if ( i == 9 ) then
    r = 139
    g = 105
    b =  20
!
!  Medium Purple
!
  else if ( i == 10 ) then
    r = 147
    g = 112
    b = 219
!
!  Coral
!
   else if ( i == 11 ) then

    r = 255
    g = 127
    b =  80
!
!  Pink5
!
  else if ( i == 12 ) then
    r = 139
    g =  99
    b = 108
!
!  GreenYellow
!
  else if ( i == 13 ) then
    r = 173
    g = 255
    b =  47
!
!  Aquamarine
!
  else if ( i == 14 ) then
    r = 127
    g = 255
    b = 212
!
!  Pale Green3
!
  else if ( i == 15 ) then
    r = 124
    g = 205
    b = 124
!
!  Burlywood
!
  else if ( i == 16 ) then
    r = 222
    g = 184
    b = 135
!
!  Cornsilk3
!
  else if ( i == 17 ) then
    r = 205
    g = 200
    b = 177
!
!  Lemon_Chiffon3
!
  else if ( i == 18 ) then
    r = 205
    g = 201
    b = 165
!
!  Maroon
!
  else if ( i == 19 ) then
    r = 176
    g = 48
    b = 96
!
!  Slate_Blue2
!
  else if ( i == 20 ) then
    r = 131
    g = 111
    b = 255

  else

    call i4_to_angle ( i, angle )

    call angle_to_rgb ( angle, rr, rg, rb )

    r = min ( int ( rr * 255 ), 255 )
    g = min ( int ( rg * 255 ), 255 )
    b = min ( int ( rb * 255 ), 255 )

  end if

  return
end
function i4_uniform ( a, b, seed )

!*****************************************************************************80
!
!! I4_UNIFORM returns a scaled pseudorandom I4.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ) value.
!
!    The pseudorandom number will be scaled to be uniformly distributed
!    between A and B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Second Edition,
!    Springer, 1987,
!    ISBN: 0387964673,
!    LC: QA76.9.C65.B73.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, December 1986, pages 362-376.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998,
!    ISBN: 0471134031,
!    LC: T57.62.H37.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, Number 2, 1969, pages 136-143.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A, B, the limits of the interval.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, integer ( kind = 4 ) I4_UNIFORM, a number between A and B.
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) k
  real    ( kind = 4 ) r
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) value

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_UNIFORM - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + i4_huge
  end if

  r = real ( seed, kind = 4 ) * 4.656612875E-10
!
!  Scale R to lie between A-0.5 and B+0.5.
!
  r = ( 1.0E+00 - r ) * ( real ( min ( a, b ), kind = 4 ) - 0.5E+00 ) & 
    +             r   * ( real ( max ( a, b ), kind = 4 ) + 0.5E+00 )
!
!  Use rounding to convert R to an integer between A and B.
!
  value = nint ( r, kind = 4 )

  value = max ( value, min ( a, b ) )
  value = min ( value, max ( a, b ) )

  i4_uniform = value

  return
end
function i4_wrap ( ival, ilo, ihi )

!*****************************************************************************80
!
!! I4_WRAP forces an integer to lie between given limits by wrapping.
!
!  Example:
!
!    ILO = 4, IHI = 8
!
!    I  I4_WRAP
!
!    -2     8
!    -1     4
!     0     5
!     1     6
!     2     7
!     3     8
!     4     4
!     5     5
!     6     6
!     7     7
!     8     8
!     9     4
!    10     5
!    11     6
!    12     7
!    13     8
!    14     4
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IVAL, an integer value.
!
!    Input, integer ( kind = 4 ) ILO, IHI, the desired bounds for the integer value.
!
!    Output, integer ( kind = 4 ) I4_WRAP, a "wrapped" version of IVAL.
!
  implicit none

  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) wide

  wide = ihi + 1 - ilo

  if ( wide == 0 ) then
    i4_wrap = ilo
  else
    i4_wrap = ilo + i4_modp ( ival-ilo, wide )
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
!    Input, character ( len = * ) TITLE, an optional title.
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
!    Input, character ( len = * ) TITLE, an optional title.
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

  write ( *, '(a)' ) ' '

  return
end
subroutine i4vec_heap_d ( n, a )

!*****************************************************************************80
!
!! I4VEC_HEAP_D reorders an array of integers into a descending heap.
!
!  Discussion:
!
!    A descending heap is an array A with the property that, for every index J,
!    A(2*J) <= A(J) and A(2*J+1) <= A(J), (as long as the indices
!    2*J and 2*J+1 are legal).
!
!  Diagram:
!
!                  A(1)
!                /      \
!            A(2)         A(3)
!          /     \        /  \
!      A(4)       A(5)  A(6) A(7)
!      /  \       /   \
!    A(8) A(9) A(10) A(11)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    A Nijenhuis and H Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the input array.
!
!    Input/output, integer ( kind = 4 ) A(N).
!    On input, an unsorted array.
!    On output, the array has been reordered into a heap.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ifree
  integer ( kind = 4 ) key
  integer ( kind = 4 ) m
!
!  Only nodes N/2 down to 1 can be "parent" nodes.
!
  do i = n/2, 1, -1
!
!  Copy the value out of the parent node.
!  Position IFREE is now "open".
!
    key = a(i)
    ifree = i

    do
!
!  Positions 2*IFREE and 2*IFREE + 1 are the descendants of position
!  IFREE.  (One or both may not exist because they exceed N.)
!
      m = 2 * ifree
!
!  Does the first position exist?
!
      if ( n < m ) then
        exit
      end if
!
!  Does the second position exist?
!
      if ( m + 1 <= n ) then
!
!  If both positions exist, take the larger of the two values,
!  and update M if necessary.
!
        if ( a(m) < a(m+1) ) then
          m = m + 1
        end if

      end if
!
!  If the large descendant is larger than KEY, move it up,
!  and update IFREE, the location of the free position, and
!  consider the descendants of THIS position.
!
      if ( a(m) <= key ) then
        exit
      end if

      a(ifree) = a(m)
      ifree = m

    end do
!
!  Once there is no more shifting to do, KEY moves into the free spot IFREE.
!
    a(ifree) = key

  end do

  return
end
subroutine i4vec_indicator ( n, a )

!*****************************************************************************80
!
!! I4VEC_INDICATOR sets an integer vector to the indicator vector.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Output, integer ( kind = 4 ) A(N), the array to be initialized.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i

  do i = 1, n
    a(i) = i
  end do

  return
end
subroutine i4vec_sort_heap_a ( n, a )

!*****************************************************************************80
!
!! I4VEC_SORT_HEAP_A ascending sorts an integer array using heap sort.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    A Nijenhuis and H Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input/output, integer ( kind = 4 ) A(N).
!    On input, the array to be sorted;
!    On output, the array has been sorted.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) n1

  if ( n <= 1 ) then
    return
  end if
!
!  1: Put A into descending heap form.
!
  call i4vec_heap_d ( n, a )
!
!  2: Sort A.
!
!  The largest object in the heap is in A(1).
!  Move it to position A(N).
!
  call i4_swap ( a(1), a(n) )
!
!  Consider the diminished heap of size N1.
!
  do n1 = n-1, 2, -1
!
!  Restore the heap structure of A(1) through A(N1).
!
    call i4vec_heap_d ( n1, a )
!
!  Take the largest object from A(1) and move it to A(N1).
!
    call i4_swap ( a(1), a(n1) )

  end do

  return
end
subroutine i4vec_sorted_unique ( n, a, unique_num )

!*****************************************************************************80
!
!! I4VEC_SORTED_UNIQUE gets the unique elements in a sorted integer array.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements in A.
!
!    Input/output, integer ( kind = 4 ) A(N).  On input, the sorted
!    integer ( kind = 4 ) array.  On output, the unique elements in A.
!
!    Output, integer ( kind = 4 ) UNIQUE_NUM, the number of unique elements in A.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) itest
  integer ( kind = 4 ) unique_num

  unique_num = 0

  if ( n <= 0 ) then
    return
  end if

  unique_num = 1

  do itest = 2, n

    if ( a(itest) /= a(unique_num) ) then
      unique_num = unique_num + 1
      a(unique_num) = a(itest)
    end if

  end do

  return
end
subroutine kmeans_plot ( box_requested, box_xy, circle_requested, &
  circle_x, circle_y, circle_r, file_name, dim_num, &
  node_num, coord, plot_min, plot_max, x, y, connect, tag, shadow, title )

!*****************************************************************************80
!
!! KMEANS_PLOT plots a K-Means diagram.
!
!  Discussion:
!
!    The "unconnected" points, which have CONNECT(I) == 0, are
!    assumed to represent the centers.
!
!    The "connected" points, which have CONNECT(I) /= 0, are assumed
!    to be the data values to be clustered.  The third coordinate of
!    such points is presumed to be an integer that indicates which
!    cluster they belong to.
!
!    A Voronoi diagram is made of the centers.
!
!    Then the connected points are displayed, with colors corresponding
!    to their assigned center.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 July 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the output file.
!
!    Input, integer ( kind = 4 ) DIM_NUM, the number of spatial dimensions.
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of points.
!
!    Input, real ( kind = 8 ) COORD(DIM_NUM,NODE_NUM), the point coordinates.
!
!    Input, real ( kind = 8 ) PLOT_MIN(2), PLOT_MAX(2), the minimum and maximum
!    X and Y values to plot.
!
!    Input, integer ( kind = 4 ) X, Y, the indices of the two dimensions to plot.
!
!    Input, integer ( kind = 4 ) CONNECT(NODE_NUM), connection indicator.
!    If CONNECT(I) is:
!    0, it is an isolated point;
!    1, it is connected to the previous point only;
!    2, it is connected to the next point only;
!    3, it is connected to the previous and next points.
!
!    Input, integer ( kind = 4 ) TAG(NODE_NUM), a tag for each point,
!    presumably the cluster index.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) node_num

  real ( kind = 8 ) b
  integer ( kind = 4 ) b_int
  real ( kind = 8 ) bitx
  real ( kind = 8 ) bity
  logical box_requested
  real ( kind = 8 ) box_xy(2,2)
  integer ( kind = 4 ) center_num
  logical circle_requested
  real ( kind = 8 ) circle_r
  real ( kind = 8 ) circle_x
  real ( kind = 8 ) circle_y
  integer ( kind = 4 ) connect(node_num)
  real ( kind = 8 ) coord(dim_num,node_num)
  real ( kind = 8 ) coord2(2,node_num)
  logical, parameter :: debug = .false.
  character ( len = * ) file_name
  integer ( kind = 4 ) file_unit
  real ( kind = 8 ) g
  integer ( kind = 4 ) g_int
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ind(node_num)
  logical inside
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jp1
  integer ( kind = 4 ) jp2
  integer ( kind = 4 ) k
  real ( kind = 8 ) n1
  real ( kind = 8 ) n2
  real ( kind = 8 ) plot_max(2)
  real ( kind = 8 ) plot_min(2)
  real ( kind = 8 ) r
  integer ( kind = 4 ) r_int
  logical shadow
  real ( kind = 8 ) t(2,3)
  integer ( kind = 4 ) tag(node_num)
  integer ( kind = 4 ) tag2(node_num)
  integer ( kind = 4 ) tag_lo
  integer ( kind = 4 ) tag_hi
  character ( len = * ) title
  real ( kind = 8 ) triangle_center(dim_num,2*node_num)
  integer ( kind = 4 ) triangle_neighbor(3,2*node_num)
  integer ( kind = 4 ) triangle_node(3,2*node_num)
  integer ( kind = 4 ) triangle_num
  integer ( kind = 4 ) x
  real ( kind = 8 ) xvec(4)
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ) x4
  real ( kind = 8 ) x5
  real ( kind = 8 ) x6
  integer ( kind = 4 ) y
  real ( kind = 8 ) yvec(4)
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) y3
  real ( kind = 8 ) y4
  real ( kind = 8 ) y5
  real ( kind = 8 ) y6
!
!  Determine which points represent the centers.
!
  center_num = 0
  
  if ( debug ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CENTERS:'
    write ( *, '(a)' ) ' '
  end if

  do i = 1, node_num

    if ( connect(i) == 0 ) then
      center_num = center_num + 1
      coord2(1,center_num) = coord(x,i)
      coord2(2,center_num) = coord(y,i)
      tag2(center_num) = tag(i)
      if ( debug ) then
        write ( *, '(3i4,2g14.6)' ) &
          center_num, i, tag(i), coord(x,i), coord(y,i)
      end if
    end if

  end do

  if ( center_num < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KMEANS_PLOT - Fatal error!'
    write ( *, '(a)' ) '  There are no centers.'
    return
  end if

  if ( center_num == node_num ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KMEANS_PLOT - Fatal error!'
    write ( *, '(a)' ) '  There are only centers, and no data points.'
    return
  end if
!
!  Compute the Delaunay triangulation of the centers.
!
  call dtris2 ( center_num, coord2, triangle_num, triangle_node, &
    triangle_neighbor )
!
!  Compute the intersection point of the perpendicular bisectors
!  of each Delaunay triangle.
!
  do j = 1, triangle_num

    i1 = triangle_node(1,j)
    i2 = triangle_node(2,j)
    i3 = triangle_node(3,j)

    t(1:2,1) = coord2(1:2,i1)
    t(1:2,2) = coord2(1:2,i2)
    t(1:2,3) = coord2(1:2,i3)

    call triangle_circumcenter_2d ( t, triangle_center(1:2,j) )

  end do

  if ( debug ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KMEANS_PLOT:'
    write ( *, '(a)' ) '  The data has been calculated.'
  end if
!
!  Plot the Voronoi tessellation.
!
  call get_unit ( file_unit )

  call ps_file_open ( file_name, file_unit, ierror )

  call eps_file_head ( file_name, 36, 36, 576, 756 )

  call ps_page_head ( plot_min(1), plot_min(2), plot_max(1), plot_max(2) )
!
!  Initialize the line width to 1.
!
  call ps_line_width ( 1 )
!
!  Print title, if requested.
!
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
!
!  Shadow the points, if requested.
!
  if ( shadow ) then

    bitx = 0.025D+00 * ( plot_max ( 1 ) - plot_min ( 1 ) )
    bity = 0.025D+00 * ( plot_max ( 2 ) - plot_min ( 2 ) )

    do i = 1, node_num
      call ps_line ( plot_min(1), coord(y,i), plot_min(1) + bitx, coord(y,i) )
      call ps_line ( coord(x,i), plot_min(2), coord(x,i), plot_min(2) + bity )
    end do

  end if
!
!  Draw a user-specified box, if requested.
!
  if ( box_requested ) then
    call ps_comment ( 'User-requested box:' )
    call ps_line ( box_xy(1,1), box_xy(2,1), box_xy(1,2), box_xy(2,1) )
    call ps_line ( box_xy(1,2), box_xy(2,1), box_xy(1,2), box_xy(2,2) )
    call ps_line ( box_xy(1,2), box_xy(2,2), box_xy(1,1), box_xy(2,2) )
    call ps_line ( box_xy(1,1), box_xy(2,2), box_xy(1,1), box_xy(2,1) )
  end if
!
!  Draw a user-specified circle, if requested.
!
  if ( circle_requested ) then
    call ps_comment ( 'User-requested circle:' )
    call ps_circle ( circle_x, circle_y, circle_r )
  end if

  x1 = plot_min(1)
  y1 = plot_min(2)
  call ps_moveto ( x1, y1 )

! call ps_label ( 'K-Means diagram' )

  call ps_setting_int ( 'SET', 'MARKER_SIZE', 8 )

  tag_lo = 1
  tag_hi = center_num

  do i = 1, center_num
    call i4_to_rgb ( tag2(i), r_int, g_int, b_int )
    r = real ( r_int, kind = 8 ) / 255.0D+00
    g = real ( g_int, kind = 8 ) / 255.0D+00
    b = real ( b_int, kind = 8 ) / 255.0D+00
    call ps_color_fill_set ( r, g, b )
    call ps_mark_disk ( coord2(1,i), coord2(2,i) )
  end do
!
!  Mark the cluster generators.
!
  call ps_setting_int ( 'SET', 'MARKER_SIZE', 4 )

  do i = 1, node_num
    if ( connect(i) /= 0 ) then
      call i4_to_rgb ( tag(i), r_int, g_int, b_int )
      r = real ( r_int, kind = 8 ) / 255.0D+00
      g = real ( g_int, kind = 8 ) / 255.0D+00
      b = real ( b_int, kind = 8 ) / 255.0D+00
      call ps_color_line_set ( r, g, b )
      call ps_mark_circle ( coord(x,i), coord(y,i) )
    end if
  end do
!
!  Construct the Voronoi region boundaries.
!
  r = 0.0D+00
  g = 0.0D+00
  b = 0.0D+00
  call ps_color_fill_set ( r, g, b )
!
!  For each Delaunay triangle, I
!  For each side J,
!
  do i = 1, triangle_num
    do j = 1, 3
      k = triangle_neighbor(j,i)
!
!  If there is a neighboring triangle K on that side,
!  connect the circumcenters.
!
      if ( 0 < k ) then

        if ( i < k ) then
          call ps_line ( triangle_center(1,i), triangle_center(2,i), &
            triangle_center(1,k), triangle_center(2,k) )
        end if
!
!  If there is no neighboring triangle on that side,
!  extend a line from the circumcenter of I in the direction of the
!  outward normal to that side.
!
      else

        ix = triangle_node(j,i)
        x1 = coord2(1,ix)
        y1 = coord2(2,ix)

        jp1 = i4_wrap ( j+1, 1, 3 )

        ix = triangle_node(jp1,i)
        x2 = coord2(1,ix)
        y2 = coord2(2,ix)

        jp2 = i4_wrap ( j+2, 1, 3 )

        ix = triangle_node(jp2,i)
        x3 = coord2(1,ix)
        y3 = coord2(2,ix)

        x4 = triangle_center(1,i)
        y4 = triangle_center(2,i)

        call line_exp_normal_2d ( x1, y1, x2, y2, n1, n2 )

        x5 = x4 + n1
        y5 = y4 + n2

        call box_ray_int_2d ( plot_min(1), plot_min(2), plot_max(1), &
          plot_max(2), x4, y4, x5, y5, x6, y6 )

        call ps_line ( x4, y4, x6, y6 )

      end if

    end do
  end do

  call ps_color_line ( 'POP', r, g, b )

  call ps_page_tail ( )

  call eps_file_tail ( )

  call ps_file_close ( file_unit )

  return
end
subroutine line_exp2imp_2d ( x1, y1, x2, y2, a, b, c )

!*****************************************************************************80
!
!! LINE_EXP2IMP_2D converts an explicit line to implicit form in 2D.
!
!  Discussion:
!
!    The explicit form of a line in 2D is:
!
!      (X1,Y1), (X2,Y2).
!
!    The implicit form of a line in 2D is:
!
!      A * X + B * Y + C = 0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X1, Y1, X2, Y2.  (X1,Y1) and (X2,Y2) are
!    two points on the line. (X1,Y1) must be different
!    from (X2,Y2).
!
!    Output, real ( kind = 8 ) A, B, C, three coefficients which describe
!    the line that passes through (X1,Y1) and (X2,Y2).
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
!
!  Take care of degenerate cases.
!
  if ( x1 == x2 .and. y1 == y2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LINE_EXP2IMP_2D - Fatal error!'
    write ( *, '(a)' ) '  (X1,Y1) = (X2,Y2)'
    write ( *, '(a,2g14.6)' ) '  (X1,Y1) = ', x1, y1
    write ( *, '(a,2g14.6)' ) '  (X2,Y2) = ', x2, y2
    stop
  end if

  a = y2 - y1
  b = x1 - x2
  c = x2 * y1 - x1 * y2

  return
end
subroutine line_exp_normal_2d ( x1, y1, x2, y2, n1, n2 )

!*****************************************************************************80
!
!! LINE_EXP_NORMAL_2D computes the normal to a line in 2D.
!
!  Discussion:
!
!    The explicit form of a line in 2D is:
!
!      (X1,Y1), (X2,Y2).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 June 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X1, Y1, X2, Y2, two points on the line.
!
!    Output, real ( kind = 8 ) N1, N2, the components of a unit normal
!    vector to the line.  If the two points are equal, then N1 = N2 = 0.
!
  implicit none

  real ( kind = 8 ) n1
  real ( kind = 8 ) n2
  real ( kind = 8 ) norm
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2

  norm = sqrt ( ( x2 - x1 )**2 + ( y2 - y1 )**2 )

  if ( norm == 0.0D+00 ) then
    n1 = 0.0D+00
    n2 = 0.0D+00
    return
  end if

  n1 =   ( y2 - y1 ) / norm
  n2 = - ( x2 - x1 ) / norm

  return
end
subroutine lines_exp_int_2d ( x1, y1, x2, y2, x3, y3, x4, y4, ival, x, y )

!*****************************************************************************80
!
!! LINES_EXP_INT_2D determines where two explicit lines intersect in 2D.
!
!  Discussion:
!
!    The explicit form of a line in 2D is:
!
!      (X1,Y1), (X2,Y2).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 June 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X1, Y1, X2, Y2, define the first line.
!
!    Input, real ( kind = 8 ) X3, Y3, X4, Y4, define the second line.
!
!    Output, integer ( kind = 4 ) IVAL, reports on the intersection:
!    0, no intersection, the lines may be parallel or degenerate.
!    1, one intersection point, returned in X, Y.
!    2, infinitely many intersections, the lines are identical.
!
!    Output, real ( kind = 8 ) X, Y, if IVAl = 1, then X, Y contains
!    the intersection point.  Otherwise, X = 0, Y = 0.
!
  implicit none

  real ( kind = 8 ) a1
  real ( kind = 8 ) a2
  real ( kind = 8 ) b1
  real ( kind = 8 ) b2
  real ( kind = 8 ) c1
  real ( kind = 8 ) c2
  integer ( kind = 4 ) ival
  logical point_1
  logical point_2
  real ( kind = 8 ) x
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ) x4
  real ( kind = 8 ) y
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) y3
  real ( kind = 8 ) y4

  ival = 0
  x = 0.0D+00
  y = 0.0D+00
!
!  Check whether either line is a point.
!
  if ( x1 == x2 .and. y1 == y2 ) then
    point_1 = .true.
  else
    point_1 = .false.
  end if

  if ( x3 == x4 .and. y3 == y4 ) then
    point_2 = .true.
  else
    point_2 = .false.
  end if
!
!  Convert the lines to ABC format.
!
  if ( .not. point_1 ) then
    call line_exp2imp_2d ( x1, y1, x2, y2, a1, b1, c1 )
  end if

  if ( .not. point_2 ) then
    call line_exp2imp_2d ( x3, y3, x4, y4, a2, b2, c2 )
  end if
!
!  Search for intersection of the lines.
!
  if ( point_1 .and. point_2 ) then
    if ( x1 == x3 .and. y1 == y3 ) then
      ival = 1
      x = x1
      y = y1
    end if
  else if ( point_1 ) then
    if ( a2 * x1 + b2 * y1 == c2 ) then
      ival = 1
      x = x1
      y = y1
    end if
  else if ( point_2 ) then
    if ( a1 * x3 + b1 * y3 == c1 ) then
      ival = 1
      x = x3
      y = y3
    end if
  else
    call lines_imp_int_2d ( a1, b1, c1, a2, b2, c2, ival, x, y )
  end if

  return
end
subroutine lines_imp_int_2d ( a1, b1, c1, a2, b2, c2, ival, x, y )

!*****************************************************************************80
!
!! LINES_IMP_INT_2D determines where two implicit lines intersect in 2D.
!
!  Discussion:
!
!    The implicit form of a line in 2D is:
!
!      A * X + B * Y + C = 0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A1, B1, C1, define the first line.
!    At least one of A1 and B1 must be nonzero.
!
!    Input, real ( kind = 8 ) A2, B2, C2, define the second line.
!    At least one of A2 and B2 must be nonzero.
!
!    Output, integer ( kind = 4 ) IVAL, reports on the intersection.
!
!    -1, both A1 and B1 were zero.
!    -2, both A2 and B2 were zero.
!     0, no intersection, the lines are parallel.
!     1, one intersection point, returned in X, Y.
!     2, infinitely many intersections, the lines are identical.
!
!    Output, real ( kind = 8 ) X, Y, if IVAL = 1, then X, Y contains
!    the intersection point.  Otherwise, X = 0, Y = 0.
!
  implicit none

  real ( kind = 8 ) a(2,2)
  real ( kind = 8 ) a1
  real ( kind = 8 ) a2
  real ( kind = 8 ) b(2,2)
  real ( kind = 8 ) b1
  real ( kind = 8 ) b2
  real ( kind = 8 ) c1
  real ( kind = 8 ) c2
  real ( kind = 8 ) det
  integer ( kind = 4 ) ival
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  x = 0.0D+00
  y = 0.0D+00
!
!  Refuse to handle degenerate lines.
!
  if ( a1 == 0.0D+00 .and. b1 == 0.0D+00 ) then
    ival = -1
    return
  else if ( a2 == 0.0D+00 .and. b2 == 0.0D+00 ) then
    ival = -2
    return
  end if
!
!  Set up a linear system, and compute its inverse.
!
  a(1,1) = a1
  a(1,2) = b1
  a(2,1) = a2
  a(2,2) = b2

  call r8mat2_inverse ( a, b, det )
!
!  If the inverse exists, then the lines intersect.
!  Multiply the inverse times -C to get the intersection point.
!
  if ( det /= 0.0D+00 ) then

    ival = 1
    x = - b(1,1) * c1 - b(1,2) * c2
    y = - b(2,1) * c1 - b(2,2) * c2
!
!  If the inverse does not exist, then the lines are parallel
!  or coincident.  Check for parallelism by seeing if the
!  C entries are in the same ratio as the A or B entries.
!
  else

    ival = 0

    if ( a1 == 0.0D+00 ) then
      if ( b2 * c1 == c2 * b1 ) then
        ival = 2
      end if
    else
      if ( a2 * c1 == c2 * a1 ) then
        ival = 2
      end if
    end if

  end if

  return
end
function lrline ( xu, yu, xv1, yv1, xv2, yv2, dv )

!*****************************************************************************80
!
!! LRLINE determines if a point is left of, right or, or on a directed line.
!
!  Discussion:
!
!    The directed line is parallel to, and at a signed distance DV from
!    a directed base line from (XV1,YV1) to (XV2,YV2).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 July 2001
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XU, YU, the coordinates of the point whose
!    position relative to the directed line is to be determined.
!
!    Input, real ( kind = 8 ) XV1, YV1, XV2, YV2, the coordinates of two points
!    that determine the directed base line.
!
!    Input, real ( kind = 8 ) DV, the signed distance of the directed line
!    from the directed base line through the points (XV1,YV1) and (XV2,YV2).
!    DV is positive for a line to the left of the base line.
!
!    Output, integer ( kind = 4 ) LRLINE, the result:
!    +1, the point is to the right of the directed line;
!     0, the point is on the directed line;
!    -1, the point is to the left of the directed line.
!
  implicit none

  real ( kind = 8 ) dv
  real ( kind = 8 ) dx
  real ( kind = 8 ) dxu
  real ( kind = 8 ) dy
  real ( kind = 8 ) dyu
  integer ( kind = 4 ) lrline
  real ( kind = 8 ) t
  real ( kind = 8 ) tol
  real ( kind = 8 ) tolabs
  real ( kind = 8 ) xu
  real ( kind = 8 ) xv1
  real ( kind = 8 ) xv2
  real ( kind = 8 ) yu
  real ( kind = 8 ) yv1
  real ( kind = 8 ) yv2

  tol = 100.0D+00 * epsilon ( tol )

  dx = xv2 - xv1
  dy = yv2 - yv1
  dxu = xu - xv1
  dyu = yu - yv1

  tolabs = tol * max ( abs ( dx ), abs ( dy ), abs ( dxu ), &
    abs ( dyu ), abs ( dv ) )

  t = dy * dxu - dx * dyu + dv * sqrt ( dx * dx + dy * dy )

  if ( tolabs < t ) then
    lrline = 1
  else if ( -tolabs <= t ) then
    lrline = 0
  else
    lrline = -1
  end if

  return
end
subroutine perm_inv ( n, p )

!*****************************************************************************80
!
!! PERM_INV inverts a permutation "in place".
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of objects being permuted.
!
!    Input/output, integer ( kind = 4 ) P(N), the permutation, in standard index form.
!    On output, P describes the inverse permutation
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i0
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) is
  integer ( kind = 4 ) p(n)

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PERM_INV - Fatal error!'
    write ( *, '(a,i6)' ) '  Input value of N = ', n
    stop
  end if

  is = 1

  do i = 1, n

    i1 = p(i)

    do while ( i < i1 )
      i2 = p(i1)
      p(i1) = -i2
      i1 = i2
    end do

    is = -sign ( 1, p(i) )
    p(i) = sign ( p(i), is )

  end do

  do i = 1, n

    i1 = -p(i)

    if ( 0 <= i1 ) then

      i0 = i

      do

        i2 = p(i1)
        p(i1) = i0

        if ( i2 < 0 ) then
          exit
        end if

        i0 = i1
        i1 = i2

      end do

    end if

  end do

  return
end
subroutine plot_size ( dim_num, box_requested, box_xy, circle_requested, &
  circle_x, circle_y, circle_r, x, y, plot_min, plot_max )

!*****************************************************************************80
!
!! PLOT_SIZE determines a plot range.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the number of spatial dimensions.
!
!    Input, logical BOX_REQUESTED, is true if the user has specified
!    a box to be drawn.
!
!    Input, real ( kind = 8 ) BOX_XY(2,2), the coordinates of the lower left and
!    upper right corners of the requested box.
!
!    Input, integer ( kind = 4 ) X, Y, the dimensions to be used for X and Y.
!
!    Output, real ( kind = 8 ) plot_min(2), plot_max(2), the
!    minimum and maximum coordinates to use in the plot.
!
  use AA_point_data_module

  implicit none

  logical box_requested
  real ( kind = 8 ) box_xy(2,2)
  logical circle_requested
  real ( kind = 8 ) circle_r
  real ( kind = 8 ) circle_x
  real ( kind = 8 ) circle_y
  integer ( kind = 4 ) dim_num
  real ( kind = 8 ), parameter :: grace = 0.05D+00
  integer ( kind = 4 ) i
  real ( kind = 8 ) plot_max(2)
  real ( kind = 8 ) plot_min(2)
  real ( kind = 8 ) width
  integer ( kind = 4 ) x
  integer ( kind = 4 ) y
  integer ( kind = 4 ) z

  do i = 1, 2

    if ( i == 1 ) then
      z = x
    else if ( i == 2 ) then
      z = y
    end if

    plot_max(i) = coord_max(z)
    plot_min(i) = coord_min(z)

    if ( box_requested ) then
      plot_max(i) = max ( plot_max(i), box_xy(i,2) )
      plot_min(i) = min ( plot_min(i), box_xy(i,1) )
    end if

    if ( circle_requested ) then
      if ( i == 1 ) then
        plot_max(i) = max ( plot_max(i), circle_x + circle_r )
        plot_min(i) = min ( plot_min(i), circle_x - circle_r )
      else if ( i == 2 ) then
        plot_max(i) = max ( plot_max(i), circle_y + circle_r )
        plot_min(i) = min ( plot_min(i), circle_y - circle_r )
      end if
    end if

    width = plot_max(i) - plot_min(i)

    plot_max(i) = plot_max(i) + grace * width
    plot_min(i) = plot_min(i) - grace * width

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '        Plot Plot'
  write ( *, '(a)' ) ' Coord  Min  Max'
  write ( *, '(a)' ) ' '
  do i = 1, 2
    write ( *, '(i3,2f10.4)' ) i, plot_min(i), plot_max(i)
  end do

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
function points_avoid_point_nd ( dim_num, n, xy_set, xy_test, tol )

!*****************************************************************************80
!
!! POINTS_AVOID_POINT_ND checks if a point is "far" from a set of points in ND.
!
!  Discussion:
!
!    The routine discards points that are too close to other points.
!    The method used to check this is quadratic in the number of points,
!    and may take an inordinate amount of time if there are a large
!    number of points.
!
!    The test point is "far enough" from an accepted point if
!    the Euclidean distance is at least "TOL".  In graphics, you probably
!    can't distinguish items that are 1/100 of the graph size apart,
!    and surely not 1/1000, so if WIDTH is the size of your image, you
!    might try setting TOL to WIDTH/100.
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
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points in the set.
!
!    Input, real ( kind = 8 ) XY_SET(DIM_NUM,N), the set of points 
!    to be avoided.
!
!    Input, real ( kind = 8 ) XY_TEST(DIM_NUM), a point to be tested.
!
!    Input, real ( kind = 8 ) TOL, the tolerance to be used.
!
!    Output, logical POINTS_AVOID_POINT_ND, is TRUE if XY_TEST is
!    "far enough" from all the set of points.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) dim_num

  integer ( kind = 4 ) j
  logical points_avoid_point_nd
  real ( kind = 8 ) tol
  real ( kind = 8 ) xy_set(dim_num,n)
  real ( kind = 8 ) xy_test(dim_num)

  points_avoid_point_nd = .true.

  do j = 1, n

    if ( sum ( ( xy_set(1:dim_num,j) - xy_test(1:dim_num) )**2 ) < tol**2 ) then
      points_avoid_point_nd = .false.
      return
    end if

  end do

  return
end
subroutine points_count ( file_in_name, dim_num, node_num )

!*****************************************************************************80
!
!! POINTS_COUNT counts the valid point coordinates in a file.
!
!  Discussion:
!
!    The routine reads every line, and expects to find DIM_NUM
!    real numbers on the line.  
!
!    It does not count lines that begin with a comment symbol '#'.
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
!  Parameters:
!
!    Input, character ( len = * ) FILE_IN_NAME, the name of the input file.
!
!    Input, integer ( kind = 4 ) DIM_NUM, the number of spatial dimensions.
!
!    Output, integer ( kind = 4 ) NODE_NUM, the number of point coordinate records.
!
  implicit none

  integer ( kind = 4 ) dim_num

  integer ( kind = 4 ) bad_num
  integer ( kind = 4 ) comment_num
  character ( len = * ) file_in_name
  integer ( kind = 4 ) file_in_unit
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  character ( len = 100 ) line
  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) record_num
  real ( kind = 8 ) x(dim_num)

  call get_unit ( file_in_unit )

  open ( unit = file_in_unit, file = file_in_name, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'POINTS_COUNT - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input file: ' // &
      trim ( file_in_name )
    stop
  end if

  comment_num = 0
  node_num = 0
  record_num = 0
  bad_num = 0

  do

    read ( file_in_unit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      ierror = record_num
      exit
    end if

    record_num = record_num + 1

    if ( line(1:1) == '#' ) then
      comment_num = comment_num + 1
      cycle
    end if

    call s_to_r8vec ( line, dim_num, x, ierror )

    if ( ierror == 0 ) then
      node_num = node_num + 1
    else
      bad_num = bad_num + 1
    end if

  end do

  close ( unit = file_in_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'POINTS_COUNT:'
  write ( *, '(a,i6)' ) '  Number of records:         ', record_num
  write ( *, '(a,i6)' ) '  Number of point records:   ', node_num
  write ( *, '(a,i6)' ) '  Number of comment records: ', comment_num
  write ( *, '(a,i6)' ) '  Number of bad records:     ', bad_num

  return
end
subroutine points_read ( connect, file_in_name, dim_num, node_num, coord )

!*****************************************************************************80
!
!! POINTS_READ reads point coordinates from a file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) CONNECT(NODE_NUM), connection indicator.
!    If CONNECT(I) is:
!    0, it is an isolated point;
!    1, it is connected to the previous point only;
!    2, it is connected to the next point only;
!    3, it is connected to the previous and next points.
!
!    Input, character ( len = * ) FILE_IN_NAME, the name of the input file.
!
!    Input, integer ( kind = 4 ) DIM_NUM, the number of spatial dimensions.
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of points.  The program
!    will stop reading data once NODE_NUM values have been read.
!
!    Output, real ( kind = 8 ) COORD(DIM_NUM,NODE_NUM), the point coordinates.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) node_num

  integer ( kind = 4 ) connect(node_num)
  real ( kind = 8 ) coord(dim_num,node_num)
  character ( len = * ) file_in_name
  integer ( kind = 4 ) file_in_unit
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ios
  character ( len = 100 ) line
  integer ( kind = 4 ) p
  real ( kind = 8 ) x(dim_num)

  call get_unit ( file_in_unit )

  open ( unit = file_in_unit, file = file_in_name, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'POINTS_READ - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input file: ' // &
      trim ( file_in_name )
    stop
  end if

  connect(1:node_num) = 0

  i = 0
  p = 0

  do while ( i < node_num )

    read ( file_in_unit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      ierror = i
      exit
    end if

    if ( line(1:1) == '#' .or. len_trim ( line ) == 0 ) then
      p = 0
      cycle
    end if

    call s_to_r8vec ( line, dim_num, x, ierror )

    if ( ierror /= 0 ) then
      p = 0
      cycle
    end if

    i = i + 1

    coord(1:dim_num,i) = x(1:dim_num)
    connect(i) = p

    if ( p == 1 .and. 1 < i ) then
      connect(i-1) = connect(i-1) + 2
    end if

    p = 1

  end do

  close ( unit = file_in_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'POINTS_READ:'
  write ( *, '(a,i6)' ) '  Read coordinate data from file.'

  return
end
subroutine points_thin ( connect, dim_num, node_num, coord, tol, node_num2 )

!*****************************************************************************80
!
!! POINTS_THIN "thins" points that are too close to each other.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) CONNECT(NODE_NUM), connection indicator.
!    If CONNECT(I) is:
!    0, it is an isolated point;
!    1, it is connected to the previous point only;
!    2, it is connected to the next point only;
!    3, it is connected to the previous and next points.
!
!    Input, integer ( kind = 4 ) DIM_NUM, the number of spatial dimensions.
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of points.
!
!    Input/output, real ( kind = 8 ) COORD(DIM_NUM,NODE_NUM), the 
!    point coordinates.
!
!    Input, real ( kind = 8 ) TOL, a tolerance to use for thinning.
!
!    Output, integer ( kind = 4 ) NODE_NUM2, the number of points kept.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) node_num

  integer ( kind = 4 ) check
  integer ( kind = 4 ) connect(node_num)
  real ( kind = 8 ) coord(dim_num,node_num)
  integer ( kind = 4 ) node_num2
  integer ( kind = 4 ) thin_num
  logical points_avoid_point_nd
  real ( kind = 8 ) tol

  node_num2 = 0
  thin_num = 0

  do check = 1, node_num
!
!  Always keep the first point.
!
    if ( check == 1 ) then

      node_num2 = 1
!
!  Keep the next point if it's far enough away.
!
    else if ( points_avoid_point_nd &
      ( dim_num, node_num2, coord, coord(1,check), tol ) ) then

      node_num2 = node_num2 + 1

      coord(1:dim_num,node_num2) = coord(1:dim_num,check)
      connect(node_num2) = connect(check)

    else

      thin_num = thin_num + 1

      if ( connect(check) == 0 ) then
      else if ( connect(check) == 1 ) then
        connect(node_num2) = 1
      else if ( connect(check) == 2 ) then
        connect(check+1) = 2
      else if ( connect(check) == 3 ) then

      end if

    end if

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'POINTS_THIN:'
  write ( *, '(a,g14.6)' ) '  Thinning tolerance was ', tol
  write ( *, '(a,i6)' ) '  Number of points given was   ', node_num
  write ( *, '(a,i6)' ) '  Number of points retained    ', node_num2
  write ( *, '(a,i6)' ) '  Number of points thinned out ', thin_num

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
!    Input, real ( kind = 8 ) X0, Y0, the coordinates of the center
!    of the circle.
!
!    Input, real ( kind = 8 ) R, the radius of the circle.
!
  implicit none

  real ( kind = 8 ) alpha
  integer ( kind = 4 ), parameter :: angle_max = 360
  integer ( kind = 4 ), parameter :: angle_min = 0
  integer ( kind = 4 ) plotxmin2
  integer ( kind = 4 ) plotymin2
  integer ( kind = 4 ) pr
  integer ( kind = 4 ) pxcen
  integer ( kind = 4 ) pycen
  real ( kind = 8 ) r
  integer ( kind = 4 ) state
  integer ( kind = 4 ) unit
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
subroutine ps_circle_fill ( x0, y0, r )

!*****************************************************************************80
!
!! PS_CIRCLE_FILL draws a filled circle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
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

  integer ( kind = 4 ), parameter :: n = 32

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
!    Henry McGilton, Mary Campione,
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
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
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
  integer ( kind = 4 ) state
  integer ( kind = 4 ) unit
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
!    Henry McGilton, Mary Campione,
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

  write ( unit, '(i6,a)' ) line_width, ' setlinewidth'

  call ps_setting_int ( 'SET', 'LINE_WIDTH', line_width )

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

  write ( unit, '(3i6,a)' ) pxcen, pycen, marker_size, &
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

  write ( unit, '(2i6,a)' ) pxcen, pycen, ' moveto'
  write ( unit, '(2i6,a)' ) pxcen+1, pycen, ' lineto'
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
subroutine ps_triangle ( x1, y1, x2, y2, x3, y3 )

!*****************************************************************************80
!
!! PS_TRIANGLE draws an open triangle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
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

  integer ( kind = 4 ), parameter :: n = 3

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
!    22 December 2000
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
subroutine r82vec_permute ( n, a, p )

!*****************************************************************************80
!
!! R82VEC_PERMUTE permutes an R82VEC in place.
!
!  Discussion:
!
!    This routine permutes an array of real "objects", but the same
!    logic can be used to permute an array of objects of any arithmetic
!    type, or an array of objects of any complexity.  The only temporary
!    storage required is enough to store a single object.  The number
!    of data movements made is N + the number of cycles of order 2 or more,
!    which is never more than N + N/2.
!
!  Example:
!
!    Input:
!
!      N = 5
!      P = (   2,    4,    5,    1,    3 )
!      A = ( 1.0,  2.0,  3.0,  4.0,  5.0 )
!          (11.0, 22.0, 33.0, 44.0, 55.0 )
!
!    Output:
!
!      A    = (  2.0,  4.0,  5.0,  1.0,  3.0 )
!             ( 22.0, 44.0, 55.0, 11.0, 33.0 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 January 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of objects.
!
!    Input/output, real ( kind = 8 ) A(2,N), the array to be permuted.
!
!    Input, integer ( kind = 4 ) P(N), the permutation.  P(I) = J means
!    that the I-th element of the output array should be the J-th
!    element of the input array.  P must be a legal permutation
!    of the integers from 1 to N, otherwise the algorithm will
!    fail catastrophically.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(2,n)
  real ( kind = 8 ) a_temp(2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) iget
  integer ( kind = 4 ) iput
  integer ( kind = 4 ) istart
  integer ( kind = 4 ) p(n)
!
!  Search for the next element of the permutation that has not been used.
!
  do istart = 1, n

    if ( p(istart) < 0 ) then

      cycle

    else if ( p(istart) == istart ) then

      p(istart) = - p(istart)
      cycle

    else

      a_temp(1:2) = a(1:2,istart)
      iget = istart
!
!  Copy the new value into the vacated entry.
!
      do

        iput = iget
        iget = p(iget)

        p(iput) = - p(iput)

        if ( iget < 1 .or. n < iget ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'D2VEC_PERMUTE - Fatal error!'
          stop
        end if

        if ( iget == istart ) then
          a(1:2,iput) = a_temp(1:2)
          exit
        end if

        a(1:2,iput) = a(1:2,iget)

      end do

    end if

  end do
!
!  Restore the signs of the entries.
!
  p(1:n) = -p(1:n)

  return
end
subroutine r82vec_sort_heap_index_a ( n, a, indx )

!*****************************************************************************80
!
!! R82VEC_SORT_HEAP_INDEX_A does an indexed heap ascending sort of an R82VEC.
!
!  Discussion:
!
!    The sorting is not actually carried out.  Rather an index array is
!    created which defines the sorting.  This array may be used to sort
!    or index the array, or to sort or index related arrays keyed on the
!    original array.
!
!    Once the index array is computed, the sorting can be carried out
!    "implicitly:
!
!      A(1:2,INDX(I)), I = 1 to N is sorted,
!
!    or explicitly, by the call
!
!      call D2VEC_PERMUTE ( N, A, INDX )
!
!    after which A(1:2,I), I = 1 to N is sorted.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 January 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, real ( kind = 8 ) A(2,N), an array to be index-sorted.
!
!    Output, integer ( kind = 4 ) INDX(N), the sort index.  The
!    I-th element of the sorted array is A(1:2,INDX(I)).
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(2,n)
  real ( kind = 8 ) aval(2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) indxt
  integer ( kind = 4 ) ir
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l

  if ( n < 1 ) then
    return
  end if

  if ( n == 1 ) then
    indx(1) = 1
    return
  end if

  call i4vec_indicator ( n, indx )

  l = n / 2 + 1
  ir = n

  do

    if ( 1 < l ) then

      l = l - 1
      indxt = indx(l)
      aval(1:2) = a(1:2,indxt)

    else

      indxt = indx(ir)
      aval(1:2) = a(1:2,indxt)
      indx(ir) = indx(1)
      ir = ir - 1

      if ( ir == 1 ) then
        indx(1) = indxt
        exit
      end if

    end if

    i = l
    j = l + l

    do while ( j <= ir )

      if ( j < ir ) then
        if (   a(1,indx(j)) <  a(1,indx(j+1)) .or. &
             ( a(1,indx(j)) == a(1,indx(j+1)) .and. &
               a(2,indx(j)) <  a(2,indx(j+1)) ) ) then
          j = j + 1
        end if
      end if

      if (   aval(1) <  a(1,indx(j)) .or. &
           ( aval(1) == a(1,indx(j)) .and. &
             aval(2) <  a(2,indx(j)) ) ) then
        indx(i) = indx(j)
        i = j
        j = j + j
      else
        j = ir + 1
      end if

    end do

    indx(i) = indxt

  end do

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
!    Input, character ( len = * ) TITLE, an optional title.
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
!    Input, character ( len = * ) TITLE, an optional title.
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
    write ( *, '(a)' ) ' '

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

  write ( *, '(a)' ) ' '

  return
end
subroutine r8mat2_inverse ( a, b, det )

!*****************************************************************************80
!
!! R8MAT2_INVERSE inverts a 2 by 2 real matrix using Cramer's rule.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A(2,2), the matrix to be inverted.
!
!    Output, real ( kind = 8 ) B(2,2), the inverse of the matrix A.
!
!    Output, real ( kind = 8 ) DET, the determinant of the matrix A.
!
!    If DET is zero, then A is singular, and does not have an
!    inverse.  In that case, B is simply set to zero, and a
!    message is printed.
!
!    If DET is nonzero, then its value is roughly an estimate
!    of how nonsingular the matrix A is.
!
  implicit none

  real ( kind = 8 ) a(2,2)
  real ( kind = 8 ) b(2,2)
  real ( kind = 8 ) det
!
!  Compute the determinant.
!
  det = a(1,1) * a(2,2) - a(1,2) * a(2,1)
!
!  If the determinant is zero, bail out.
!
  if ( det == 0.0D+00 ) then

    b(1:2,1:2) = 0.0D+00

    return
  end if
!
!  Compute the entries of the inverse matrix using an explicit formula.
!
  b(1,1) = + a(2,2) / det
  b(1,2) = - a(1,2) / det
  b(2,1) = - a(2,1) / det
  b(2,2) = + a(1,1) / det

  return
end
subroutine r8mat_solve ( a, n, nrhs, info )

!*****************************************************************************80
!
!! R8MAT_SOLVE uses Gauss-Jordan elimination to solve an N by N linear system.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) A(N,N+NRHS), contains in rows and columns 1
!    to N the coefficient matrix, and in columns N+1 through
!    N+NRHS, the right hand sides.  On output, the coefficient matrix
!    area has been destroyed, while the right hand sides have
!    been overwritten with the corresponding solutions.
!
!    Input, integer ( kind = 4 ) NRHS, the number of right hand sides.  NRHS
!    must be at least 0.
!
!    Output, integer ( kind = 4 ) INFO, singularity flag.
!    0, the matrix was not singular, the solutions were computed;
!    J, factorization failed on step J, and the solutions could not
!    be computed.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nrhs

  real ( kind = 8 ) a(n,n+nrhs)
  real ( kind = 8 ) apivot
  real ( kind = 8 ) factor
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipivot
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) temp

  info = 0

  do j = 1, n
!
!  Choose a pivot row.
!
    ipivot = j
    apivot = a(j,j)

    do i = j+1, n
      if ( abs ( apivot ) < abs ( a(i,j) ) ) then
        apivot = a(i,j)
        ipivot = i
      end if
    end do

    if ( apivot == 0.0D+00 ) then
      info = j
      return
    end if
!
!  Interchange.
!
    do i = 1, n + nrhs
      call r8_swap ( a(ipivot,i), a(j,i) )
    end do
!
!  A(J,J) becomes 1.
!
    a(j,j) = 1.0D+00
    a(j,j+1:n+nrhs) = a(j,j+1:n+nrhs) / apivot
!
!  A(I,J) becomes 0.
!
    do i = 1, n

      if ( i /= j ) then

        factor = a(i,j)
        a(i,j) = 0.0D+00
        a(i,j+1:n+nrhs) = a(i,j+1:n+nrhs) - factor * a(j,j+1:n+nrhs)

      end if

    end do

  end do

  return
end
function radians_to_degrees ( angle )

!*****************************************************************************80
!
!! RADIANS_TO_DEGREES converts an angle from radians to degrees.
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
!    Input, real ( kind = 8 ) ANGLE, an angle in radians.
!
!    Output, real ( kind = 8 ) RADIANS_TO_DEGREES, the equivalent 
!    angle in degrees.
!
  implicit none

  real ( kind = 8 ) angle
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) radians_to_degrees

  radians_to_degrees = ( angle / pi ) * 180.0D+00

  return
end
subroutine radius_maximus ( dim_num, node_num, coord, walls, radius )

!*****************************************************************************80
!
!! RADIUS_MAXIMUS finds the biggest possible nonintersecting sphere.
!
!  Discussion:
!
!    We are given a set of NODE_NUM points in DIM_NUM space.  We imagine that
!    at each point simultaneously, a sphere begins to expand.
!    Each sphere stops expanding as soon as it touches another sphere.
!    The radius of these spheres is to be computed.
!
!    If WALLS is true, then the spheres must not extend outside the
!    "walls" of the unit hypersquare.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 April 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the number of spatial dimensions.
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of points.
!
!    Input, real ( kind = 8 ) COORD(DIM_NUM,NODE_NUM), the point coordinates.
!
!    Input, logical WALLS, is TRUE if the spheres must not extend
!    outside the unit hypercube.  If WALLS is FALSE, then this
!    restriction is not imposed.
!
!    Output, real ( kind = 8 ) RADIUS(NODE_NUM), the radius of the
!    sphere around each point.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) node_num

  real ( kind = 8 ) coord(dim_num,node_num)
  real ( kind = 8 ) distance(node_num)
  real ( kind = 8 ) distance_j
  real ( kind = 8 ) distance_min
  integer ( kind = 4 ), parameter :: FIXED = 0
  integer ( kind = 4 ), parameter :: FREE = 1
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) next
  real ( kind = 8 ) radius(node_num)
  real ( kind = 8 ) radius_i
  real ( kind = 8 ) radius_min
  integer ( kind = 4 ) status(node_num)
  logical walls

  if ( walls ) then
         if ( any (           coord(1:dim_num,1:node_num) < 0.0D+00 ) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RADIUS_MAXIMUS - Fatal error!'
      write ( *, '(a)' ) '  Some coordinate is less than 0.'
      stop
    else if ( any ( 1.0D+00 < coord(1:dim_num,1:node_num)           ) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RADIUS_MAXIMUS - Fatal error!'
      write ( *, '(a)' ) '  Some coordinate is greater than 1.'
      stop
    end if
  end if
!
!  Initially, all points are "free".
!
  radius(1:node_num) = 0.0D+00
  status(1:node_num) = FREE

  do
!
!  If all points are fixed, we're done.
!
    if ( all ( status(1:node_num) == FIXED ) ) then
      exit
    end if
!
!  Look at all the free points.
!  Imagine an expanding sphere at each free point, and determine
!  which such sphere will first have to stop expanding.
!
    next = 0
    radius_min = huge ( radius_min )

    do i = 1, node_num

      if ( status(i) == FREE ) then

        if ( walls ) then
          radius_i = min ( &
            minval (           coord(1:dim_num,i) ), &
            minval ( 1.0D+00 - coord(1:dim_num,i) ) )
        else
          radius_i = huge ( radius_i )
        end if

        do j = 1, node_num

          if ( j /= i ) then

            distance_j = sqrt ( sum ( &
              ( coord(1:dim_num,i) - coord(1:dim_num,j) )**2 &
            ) )

            if ( status(j) == FREE ) then
              radius_i = min ( radius_i, distance_j / 2.0D+00 )
            else
              radius_i = min ( radius_i, distance_j - radius(j) )
            end if

          end if

        end do

        if ( radius_i < radius_min ) then
          next = i
          radius_min = radius_i
        end if

      end if

    end do

    if ( next == 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RADIUS_MAXIMUS - Fatal error!'
      write ( *, '(a)' ) '  There were points left to handle, but could'
      write ( *, '(a)' ) '  not choose the "next" one to work on.'
      stop
    end if

    i = next
    radius(i) = radius_min
    status(i) = FIXED

  end do

  return
end
subroutine s_blank_delete ( s )

!*****************************************************************************80
!
!! S_BLANK_DELETE removes blanks from a string, left justifying the remainder.
!
!  Discussion:
!
!    All TAB characters are also removed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) S, the string to be transformed.
!
  implicit none

  character c
  integer ( kind = 4 ) iget
  integer ( kind = 4 ) iput
  integer ( kind = 4 ) nchar
  character ( len = * ) s
  character, parameter :: TAB = char ( 9 )

  iput = 0
  nchar = len_trim ( s )

  do iget = 1, nchar

    c = s(iget:iget)

    if ( c /= ' ' .and. c /= TAB ) then
      iput = iput + 1
      s(iput:iput) = c
    end if

  end do

  s(iput+1:nchar) = ' '

  return
end
subroutine s_cap ( s )

!*****************************************************************************80
!
!! S_CAP replaces any lowercase letters by uppercase ones in a string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) S, the string to be transformed.
!
  implicit none

  character c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) nchar
  character ( len = * ) s

  nchar = len ( s )

  do i = 1, nchar

    c = s(i:i)
    call ch_cap ( c )
    s(i:i) = c

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
subroutine s_to_i4 ( s, ival, ierror, last )

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
!    Output, integer ( kind = 4 ) LAST, the last character of S used to make IVAL.
!
  implicit none

  character c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) istate
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) last
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
        last = i - 1
        return
      end if

    end if

  end do
!
!  If we read all the characters in the string, see if we're OK.
!
  if ( istate == 2 ) then
    ival = isgn * ival
    last = len_trim ( s )
  else
    ierror = 1
    last = 0
  end if

  return
end
subroutine s_to_r8 ( s, r, ierror, length )

!******************************************************************************
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
!    Output, integer ( kind = 4 ) LENGTH, the number of characters read from
!    the string to form the number, including any terminating
!    characters such as a trailing comma or blanks.
!
  implicit none

  logical ch_eqi
  character c
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
!    19 February 2001
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
  integer ( kind = 4 ) length
  real ( kind = 8 ) rvec(n)
  character ( len = * ) s

  i = 0

  ilo = 1

  do while ( i < n )

    i = i + 1

    call s_to_r8 ( s(ilo:), rvec(i), ierror, length )

    if ( ierror /= 0 ) then
      ierror = -i
      exit
    end if

    ilo = ilo + length

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
subroutine swapec ( i, top, btri, bedg, node_num, node_xy, triangle_num, &
  triangle_node, triangle_neighbor, stack, ierr )

!*****************************************************************************80
!
!! SWAPEC swaps diagonal edges until all triangles are Delaunay.
!
!  Discussion:
!
!    The routine swaps diagonal edges in a 2D triangulation, based on
!    the empty circumcircle criterion, until all triangles are Delaunay,
!    given that I is the index of the new vertex added to the triangulation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 July 2001
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the index of the new vertex.
!
!    Input/output, integer ( kind = 4 ) TOP, the index of the top of the stack.
!    On output, TOP is zero.
!
!    Input/output, integer ( kind = 4 ) BTRI, BEDG; on input, if positive, are the
!    triangle and edge indices of a boundary edge whose updated indices
!    must be recorded.  On output, these may be updated because of swaps.
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, real ( kind = 8 ) NODE_XY(2,NODE_NUM), the coordinates of the nodes.
!
!    Input, integer ( kind = 4 ) TRIANGLE_NUM, the number of triangles.
!
!    Input/output, integer ( kind = 4 ) TRIANGLE_NODE(3,TRIANGLE_NUM), the triangle incidence
!    list.  May be updated on output because of swaps.
!
!    Input/output, integer ( kind = 4 ) TRIANGLE_NEIGHBOR(3,TRIANGLE_NUM), the triangle
!    neighbor list; negative values are used for links of the counter-clockwise
!    linked list of boundary edges;  May be updated on output because of swaps.
!
!      LINK = -(3*I + J-1) where I, J = triangle, edge index.
!
!    Workspace, integer STACK(MAXST); on input, entries 1 through TOP
!    contain the indices of initial triangles (involving vertex I)
!    put in stack; the edges opposite I should be in interior;  entries
!    TOP+1 through MAXST are used as a stack.
!
!    Output, integer ( kind = 4 ) IERR is set to 8 for abnormal return.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) triangle_num

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) bedg
  integer ( kind = 4 ) btri
  integer ( kind = 4 ) c
  integer ( kind = 4 ) diaedg
  integer ( kind = 4 ) e
  integer ( kind = 4 ) ee
  integer ( kind = 4 ) em1
  integer ( kind = 4 ) ep1
  integer ( kind = 4 ) f
  integer ( kind = 4 ) fm1
  integer ( kind = 4 ) fp1
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) l
  real ( kind = 8 ) node_xy(dim_num,node_num)
  integer ( kind = 4 ) r
  integer ( kind = 4 ) s
  integer ( kind = 4 ) stack(node_num)
  integer ( kind = 4 ) swap
  integer ( kind = 4 ) t
  integer ( kind = 4 ) top
  integer ( kind = 4 ) triangle_node(3,triangle_num)
  integer ( kind = 4 ) triangle_neighbor(3,triangle_num)
  integer ( kind = 4 ) tt
  integer ( kind = 4 ) u
  real ( kind = 8 ) x
  real ( kind = 8 ) y
!
!  Determine whether triangles in stack are Delaunay, and swap
!  diagonal edge of convex quadrilateral if not.
!
  x = node_xy(1,i)
  y = node_xy(2,i)

  do

    if ( top <= 0 ) then
      exit
    end if

    t = stack(top)
    top = top - 1

    if ( triangle_node(1,t) == i ) then
      e = 2
      b = triangle_node(3,t)
    else if ( triangle_node(2,t) == i ) then
      e = 3
      b = triangle_node(1,t)
    else
      e = 1
      b = triangle_node(2,t)
    end if

    a = triangle_node(e,t)
    u = triangle_neighbor(e,t)

    if ( triangle_neighbor(1,u) == t ) then
      f = 1
      c = triangle_node(3,u)
    else if ( triangle_neighbor(2,u) == t ) then
      f = 2
      c = triangle_node(1,u)
    else
      f = 3
      c = triangle_node(2,u)
    end if

    swap = diaedg ( x, y, node_xy(1,a), node_xy(2,a), node_xy(1,c), &
      node_xy(2,c), node_xy(1,b), node_xy(2,b) )

    if ( swap == 1 ) then

      em1 = i4_wrap ( e - 1, 1, 3 )
      ep1 = i4_wrap ( e + 1, 1, 3 )
      fm1 = i4_wrap ( f - 1, 1, 3 )
      fp1 = i4_wrap ( f + 1, 1, 3 )

      triangle_node(ep1,t) = c
      triangle_node(fp1,u) = i

      r = triangle_neighbor(ep1,t)
      s = triangle_neighbor(fp1,u)

      triangle_neighbor(ep1,t) = u
      triangle_neighbor(fp1,u) = t
      triangle_neighbor(e,t) = s
      triangle_neighbor(f,u) = r

      if ( 0 < triangle_neighbor(fm1,u) ) then
        top = top + 1
        stack(top) = u
      end if

      if ( 0 < s ) then

        if ( triangle_neighbor(1,s) == u ) then
          triangle_neighbor(1,s) = t
        else if ( triangle_neighbor(2,s) == u ) then
          triangle_neighbor(2,s) = t
        else
          triangle_neighbor(3,s) = t
        end if

        top = top + 1

        if ( node_num < top ) then
          ierr = 8
          return
        end if

        stack(top) = t

      else

        if ( u == btri .and. fp1 == bedg ) then
          btri = t
          bedg = e
        end if

        l = - ( 3 * t + e - 1 )
        tt = t
        ee = em1

        do while ( 0 < triangle_neighbor(ee,tt) )

          tt = triangle_neighbor(ee,tt)

          if ( triangle_node(1,tt) == a ) then
            ee = 3
          else if ( triangle_node(2,tt) == a ) then
            ee = 1
          else
            ee = 2
          end if

        end do

        triangle_neighbor(ee,tt) = l

      end if

      if ( 0 < r ) then

        if ( triangle_neighbor(1,r) == t ) then
          triangle_neighbor(1,r) = u
        else if ( triangle_neighbor(2,r) == t ) then
          triangle_neighbor(2,r) = u
        else
          triangle_neighbor(3,r) = u
        end if

      else

        if ( t == btri .and. ep1 == bedg ) then
          btri = u
          bedg = f
        end if

        l = - ( 3 * u + f - 1 )
        tt = u
        ee = fm1

        do while ( 0 < triangle_neighbor(ee,tt) )

          tt = triangle_neighbor(ee,tt)

          if ( triangle_node(1,tt) == b ) then
            ee = 3
          else if ( triangle_node(2,tt) == b ) then
            ee = 1
          else
            ee = 2
          end if

        end do

        triangle_neighbor(ee,tt) = l

      end if

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
  integer ( kind = 4 )   d
  integer ( kind = 4 )   h
  integer ( kind = 4 )   m
  integer ( kind = 4 )   mm
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
subroutine triangle_circumcenter_2d ( t, center )

!*****************************************************************************80
!
!! TRIANGLE_CIRCUMCENTER_2D computes the circumcenter of a triangle in 2D.
!
!  Discussion:
!
!    The circumcenter of a triangle is the center of the circumcircle, the
!    circle that passes through the three vertices of the triangle.
!
!    The circumcircle contains the triangle, but it is not necessarily the
!    smallest triangle to do so.
!
!    If all angles of the triangle are no greater than 90 degrees, then
!    the center of the circumscribed circle will lie inside the triangle.
!    Otherwise, the center will lie outside the triangle.
!
!    The circumcenter is the intersection of the perpendicular bisectors
!    of the sides of the triangle.
!
!    In geometry, the circumcenter of a triangle is often symbolized by "O".
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T(2,3), the triangle vertices.
!
!    Output, real ( kind = 8 ) CENTER(2), the circumcenter of the triangle.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) asq
  real ( kind = 8 ) bot
  real ( kind = 8 ) center(dim_num)
  real ( kind = 8 ) csq
  real ( kind = 8 ) t(dim_num,3)
  real ( kind = 8 ) top(dim_num)

  asq = ( t(1,2) - t(1,1) )**2 + ( t(2,2) - t(2,1) )**2
  csq = ( t(1,3) - t(1,1) )**2 + ( t(2,3) - t(2,1) )**2
  
  top(1) =    ( t(2,2) - t(2,1) ) * csq - ( t(2,3) - t(2,1) ) * asq
  top(2) =  - ( t(1,2) - t(1,1) ) * csq + ( t(1,3) - t(1,1) ) * asq

  bot  =  ( t(2,2) - t(2,1) ) * ( t(1,3) - t(1,1) ) &
        - ( t(2,3) - t(2,1) ) * ( t(1,2) - t(1,1) )

  center(1:2) = t(1:2,1) + 0.5D+00 * top(1:2) / bot

  return
end
subroutine triangulation_order3_print ( node_num, triangle_num, node_xy, &
  triangle_node, triangle_neighbor )

!*****************************************************************************80
!
!! TRIANGULATION_ORDER3_PRINT prints out information defining a Delaunay triangulation.
!
!  Discussion:
!
!    Triangulations created by DTRIS2 include extra information encoded
!    in the negative values of TRIANGLE_NEIGHBOR.
!
!    Because some of the nodes counted in NODE_NUM may not actually be
!    used in the triangulation, I needed to compute the true number
!    of vertices.  I added this calculation on 13 October 2001.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 November 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, integer ( kind = 4 ) TRIANGLE_NUM, the number of triangles.
!
!    Input, real ( kind = 8 ) NODE_XY(2,NODE_NUM), the coordinates of the nodes.
!
!    Input, integer ( kind = 4 ) TRIANGLE_NODE(3,TRIANGLE_NUM), the nodes that make up the
!    triangles.
!
!    Input, integer ( kind = 4 ) TRIANGLE_NEIGHBOR(3,TRIANGLE_NUM), the triangle neighbors on
!    each side.  If there is no triangle neighbor on a particular side, the
!    value of TRIANGLE_NEIGHBOR should be negative.  If the triangulation 
!    data was created by DTRIS2, then there is more information encoded
!    in the negative values.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) triangle_num

  integer ( kind = 4 ) boundary_num
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  real ( kind = 8 ) node_xy(dim_num,node_num)
  integer ( kind = 4 ) s
  logical skip
  integer ( kind = 4 ) t
  integer ( kind = 4 ) triangle_node(3,triangle_num)
  integer ( kind = 4 ) triangle_neighbor(3,triangle_num)
  integer ( kind = 4 ), allocatable, dimension ( : ) :: vertex_list
  integer ( kind = 4 ) vertex_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TRIANGULATION_ORDER3_PRINT'
  write ( *, '(a)' ) '  Information defining an order3 triangulation.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  The number of nodes is ', node_num

  call r8mat_transpose_print ( dim_num, node_num, node_xy, '  Node coordinates' )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  The number of triangles is ', triangle_num
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Sets of three nodes are used as vertices of'
  write ( *, '(a)' ) '  the triangles.  For each triangle, the nodes'
  write ( *, '(a)' ) '  are listed in counterclockwise order.'

  call i4mat_transpose_print ( 3, triangle_num, triangle_node, &
    '  Triangle nodes:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  On each side of a given triangle, there is either'
  write ( *, '(a)' ) '  another triangle, or a piece of the convex hull.'
  write ( *, '(a)' ) '  For each triangle, we list the indices of the three'
  write ( *, '(a)' ) '  neighbors, or (if negative) the codes of the'
  write ( *, '(a)' ) '  segments of the convex hull.'

  call i4mat_transpose_print ( 3, triangle_num, triangle_neighbor, &
    '  Triangle neighbors' )
!
!  Determine the number of vertices.  
!
  allocate ( vertex_list(1:3*triangle_num) )

  vertex_list(1:3*triangle_num) = reshape ( triangle_node(1:3,1:triangle_num), &
    (/ 3*triangle_num /) )

  call i4vec_sort_heap_a ( 3*triangle_num, vertex_list )

  call i4vec_sorted_unique ( 3*triangle_num, vertex_list, vertex_num )

  deallocate ( vertex_list )
!
!  Determine the number of boundary points.
!
  boundary_num = 2 * vertex_num - triangle_num - 2

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  The number of boundary points is ', boundary_num
   
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The segments that make up the convex hull can be'
  write ( *, '(a)' ) '  determined from the negative entries of the triangle'
  write ( *, '(a)' ) '  neighbor list.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     #   Tri  Side    N1    N2'
  write ( *, '(a)' ) ' '

  skip = .false.

  k = 0

  do i = 1, triangle_num

    do j = 1, 3

      if ( triangle_neighbor(j,i) < 0 ) then
        s = - triangle_neighbor(j,i)
        t = s / 3

        if ( t < 1 .or. triangle_num < t ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  Sorry, this data does not use the DTRIS2'
          write ( *, '(a)' ) '  convention for convex hull segments.'
          skip = .true.
          exit
        end if

        s = mod ( s, 3 ) + 1
        k = k + 1
        n1 = triangle_node(s,t)
        n2 = triangle_node(i4_wrap(s+1,1,3),t)
        write ( *, '(2x,i4,2x,i4,2x,i4,2x,i4,2x,i4)' ) k, t, s, n1, n2
      end if

    end do

    if ( skip ) then
      exit
    end if

  end do

  return
end
subroutine trivor_plot ( box_requested, box_xy, circle_requested, &
  circle_x, circle_y, circle_r, file_name, dim_num, &
  node_num, coord, plot_min, plot_max, x, y, shadow, title )

!*****************************************************************************80
!
!! TRIVOR_PLOT plots the triangulated Voronoi diagram of a set of points in 2D.
!
!  Discussion:
!
!    The routine first determines the Delaunay triangulation.
!    The Voronoi diagram is then determined from this information.
!    Finally, each Voronoi cell is triangulated, using the cell centroid
!    and the vertices.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the output file.
!
!    Input, integer ( kind = 4 ) DIM_NUM, the number of spatial dimensions.
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of points.
!
!    Input, real ( kind = 8 ) COORD(DIM_NUM,NODE_NUM), the point coordinates.
!
!    Input, real ( kind = 8 ) PLOT_MIN(2), PLOT_MAX(2), the minimum and maximum
!    X and Y values to plot.
!
!    Input, integer ( kind = 4 ) X, Y, the indices of the two dimensions to plot.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) node_num

  real ( kind = 8 ) b
  real ( kind = 8 ) bitx
  real ( kind = 8 ) bity
  logical box_requested
  real ( kind = 8 ) box_xy(2,2)
  logical circle_requested
  real ( kind = 8 ) circle_r
  real ( kind = 8 ) circle_x
  real ( kind = 8 ) circle_y
  real ( kind = 8 ) coord(dim_num,node_num)
  real ( kind = 8 ) coord2(2,node_num)
  logical, parameter :: debug = .false.
  character ( len = * ) file_name
  integer ( kind = 4 ) file_unit
  real ( kind = 8 ) g
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) ierror
  logical inside
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jp1
  integer ( kind = 4 ) jp2
  integer ( kind = 4 ) k
  real ( kind = 8 ) n1
  real ( kind = 8 ) n2
  real ( kind = 8 ) plot_max(2)
  real ( kind = 8 ) plot_min(2)
  real ( kind = 8 ) r
  logical shadow
  real ( kind = 8 ) t(2,3)
  character ( len = * ) title
  real ( kind = 8 ) triangle_center(dim_num,2*node_num)
  integer ( kind = 4 ) triangle_neighbor(3,2*node_num)
  integer ( kind = 4 ) triangle_node(3,2*node_num)
  integer ( kind = 4 ) triangle_num
  integer ( kind = 4 ) x
  real ( kind = 8 ) xvec(4)
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ) x4
  real ( kind = 8 ) x5
  real ( kind = 8 ) x6
  integer ( kind = 4 ) y
  real ( kind = 8 ) yvec(4)
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) y3
  real ( kind = 8 ) y4
  real ( kind = 8 ) y5
  real ( kind = 8 ) y6
!
!  Compute the Delaunay triangulation.
!
  coord2(1,1:node_num) = coord(x,1:node_num)
  coord2(2,1:node_num) = coord(y,1:node_num)

  call dtris2 ( node_num, coord2, triangle_num, triangle_node, &
    triangle_neighbor )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TRIVOR_PLOT:'
  write ( *, '(a,i6)' ) '  The number of triangles is ', triangle_num

  if ( debug ) then
    call r8mat_transpose_print ( 2, node_num, coord2, '  Node coordinates:' )
  end if

  if ( debug ) then
    call i4mat_transpose_print ( 3, triangle_num, triangle_node, '  Triangles:' )
  end if

  if ( debug ) then
    call i4mat_transpose_print ( 3, triangle_num, triangle_neighbor, &
      '  Neighbors:' )
  end if
!
!  Compute the intersection point of the perpendicular bisectors
!  of each Delaunay triangle.
!
  do j = 1, triangle_num

    i1 = triangle_node(1,j)
    i2 = triangle_node(2,j)
    i3 = triangle_node(3,j)

    t(1:2,1) = coord2(1:2,i1)
    t(1:2,2) = coord2(1:2,i2)
    t(1:2,3) = coord2(1:2,i3)

    call triangle_circumcenter_2d ( t, triangle_center(1:2,j) )

  end do

  if ( debug ) then
    call r8mat_transpose_print ( 2, triangle_num, triangle_center, &
      '  Circumcenters:' )
  end if

  call get_unit ( file_unit )

  call ps_file_open ( file_name, file_unit, ierror )

  call eps_file_head ( file_name, 36, 36, 576, 756  )

  call ps_page_head ( plot_min(1), plot_min(2), plot_max(1), plot_max(2) )
!
!  Initialize the line width to 1.
!
  call ps_line_width ( 1 )
!
!  Print title, if requested.
!
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
!
!  Shadow the points, if requested.
!
  if ( shadow ) then

    bitx = 0.025D+00 * ( plot_max ( 1 ) - plot_min ( 1 ) )
    bity = 0.025D+00 * ( plot_max ( 2 ) - plot_min ( 2 ) )

    do i = 1, node_num
      call ps_line ( plot_min(1), coord(y,i), plot_min(1) + bitx, coord(y,i) )
      call ps_line ( coord(x,i), plot_min(2), coord(x,i), plot_min(2) + bity )
    end do

  end if
!
!  Draw a user-specified box, if requested.
!
  if ( box_requested ) then
    call ps_comment ( 'User-requested box:' )
    call ps_line ( box_xy(1,1), box_xy(2,1), box_xy(1,2), box_xy(2,1) )
    call ps_line ( box_xy(1,2), box_xy(2,1), box_xy(1,2), box_xy(2,2) )
    call ps_line ( box_xy(1,2), box_xy(2,2), box_xy(1,1), box_xy(2,2) )
    call ps_line ( box_xy(1,1), box_xy(2,2), box_xy(1,1), box_xy(2,1) )
  end if
!
!  Draw a user-specified circle, if requested.
!
  if ( circle_requested ) then
    call ps_comment ( 'User-requested circle:' )
    call ps_circle ( circle_x, circle_y, circle_r )
  end if

  x1 = plot_min(1)
  y1 = plot_min(2)
  call ps_moveto ( x1, y1 )

  if ( node_num <= 350 ) then
    do i = 1, node_num
      call ps_mark_disk ( coord2(1,i), coord2(2,i) )
    end do
  else
    do i = 1, node_num
      call ps_mark_point ( coord2(1,i), coord2(2,i) )
    end do
  end if

  if ( node_num <= 350 ) then
    do i = 1, triangle_num
      call ps_mark_circle ( triangle_center(1,i), triangle_center(2,i) )
    end do
  end if
!
!  For each Delaunay triangle, I
!  For each side J,
!
  do i = 1, triangle_num
    do j = 1, 3
      k = triangle_neighbor(j,i)
!
!  If there is a neighboring triangle K on that side,
!  connect the circumcenters.
!
      if ( 0 < k ) then
        if ( i < k ) then
          call ps_line ( triangle_center(1,i), triangle_center(2,i), &
            triangle_center(1,k), triangle_center(2,k) )
        end if
!
!  If there is no neighboring triangle on that side,
!  extend a line from the circumcenter of I in the direction of the
!  outward normal to that side.
!
      else

        ix = triangle_node(j,i)
        x1 = coord(x,ix)
        y1 = coord(y,ix)

        jp1 = i4_wrap ( j+1, 1, 3 )

        ix = triangle_node(jp1,i)
        x2 = coord(x,ix)
        y2 = coord(y,ix)

        jp2 = i4_wrap ( j+2, 1, 3 )

        ix = triangle_node(jp2,i)
        x3 = coord(x,ix)
        y3 = coord(y,ix)

        x4 = triangle_center(1,i)
        y4 = triangle_center(2,i)

        call line_exp_normal_2d ( x1, y1, x2, y2, n1, n2 )

        x5 = x4 + n1
        y5 = y4 + n2

        call box_ray_int_2d ( plot_min(1), plot_min(2), plot_max(1), &
          plot_max(2), x4, y4, x5, y5, x6, y6 )

        call ps_line ( x4, y4, x6, y6 )

      end if

    end do
  end do
!
!  For each triangle, join its circumcenter to each of the three nodes,
!  using a gray line.  This will display a triangulation of the Voronoi
!  polygons.
!
  r = 0.8D+00
  g = 0.8D+00
  b = 0.8D+00

  call ps_color_line_set ( r, g, b )

  do i = 1, triangle_num
    x1 = triangle_center(1,i)
    y1 = triangle_center(2,i)
    do j = 1, 3
      k = triangle_node(j,i)
      x2 = coord(x,k)
      y2 = coord(y,k)
      call ps_line ( x1, y1, x2, y2 )
    end do
  end do

  call ps_page_tail ( )

  call eps_file_tail ( )

  call ps_file_close ( file_unit )

  return
end
subroutine vbedg ( x, y, node_num, node_xy, triangle_num, triangle_node, &
  triangle_neighbor, ltri, ledg, rtri, redg )

!*****************************************************************************80
!
!! VBEDG determines which boundary edges are visible to a point.
!
!  Discussion:
!
!    The point (X,Y) is assumed to be outside the convex hull of the
!    region covered by the 2D triangulation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 August 2001
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, Y, the coordinates of a point outside the
!    convex hull of the current triangulation.
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, real ( kind = 8 ) NODE_XY(2,NODE_NUM), the coordinates of the nodes.
!
!    Input, integer ( kind = 4 ) TRIANGLE_NUM, the number of triangles.
!
!    Input, integer ( kind = 4 ) TRIANGLE_NODE(3,TRIANGLE_NUM), the triangle incidence list.
!
!    Input, integer ( kind = 4 ) TRIANGLE_NEIGHBOR(3,TRIANGLE_NUM), the triangle neighbor
!    list; negative values are used for links of a counter clockwise linked
!    list of boundary edges;
!      LINK = -(3*I + J-1) where I, J = triangle, edge index.
!
!    Input/output, integer ( kind = 4 ) LTRI, LEDG.  If LTRI /= 0 then these values are
!    assumed to be already computed and are not changed, else they are updated.
!    On output, LTRI is the index of boundary triangle to the left of the
!    leftmost boundary triangle visible from (X,Y), and LEDG is the boundary
!    edge of triangle LTRI to the left of the leftmost boundary edge visible
!    from (X,Y).  1 <= LEDG <= 3.
!
!    Input/output, integer ( kind = 4 ) RTRI.  On input, the index of the boundary triangle
!    to begin the search at.  On output, the index of the rightmost boundary
!    triangle visible from (X,Y).
!
!    Input/output, integer ( kind = 4 ) REDG, the edge of triangle RTRI that is visible
!    from (X,Y).  1 <= REDG <= 3.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) triangle_num

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) e
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) l
  logical ldone
  integer ( kind = 4 ) ledg
  integer ( kind = 4 ) lr
  integer ( kind = 4 ) lrline
  integer ( kind = 4 ) ltri
  real ( kind = 8 ) node_xy(2,node_num)
  integer ( kind = 4 ) redg
  integer ( kind = 4 ) rtri
  integer ( kind = 4 ) t
  integer ( kind = 4 ) triangle_node(3,triangle_num)
  integer ( kind = 4 ) triangle_neighbor(3,triangle_num)
  real ( kind = 8 ) x
  real ( kind = 8 ) y
!
!  Find the rightmost visible boundary edge using links, then possibly
!  leftmost visible boundary edge using triangle neighbor information.
!
  if ( ltri == 0 ) then
    ldone = .false.
    ltri = rtri
    ledg = redg
  else
    ldone = .true.
  end if

  do

    l = -triangle_neighbor(redg,rtri)
    t = l / 3
    e = mod ( l, 3 ) + 1
    a = triangle_node(e,t)

    if ( e <= 2 ) then
      b = triangle_node(e+1,t)
    else
      b = triangle_node(1,t)
    end if

    lr = lrline ( x, y, node_xy(1,a), node_xy(2,a), node_xy(1,b), &
      node_xy(2,b), 0.0D+00 )

    if ( lr <= 0 ) then
      exit
    end if

    rtri = t
    redg = e

  end do

  if ( ldone ) then
    return
  end if

  t = ltri
  e = ledg

  do

    b = triangle_node(e,t)
    e = i4_wrap ( e-1, 1, 3 )

    do while ( 0 < triangle_neighbor(e,t) )

      t = triangle_neighbor(e,t)

      if ( triangle_node(1,t) == b ) then
        e = 3
      else if ( triangle_node(2,t) == b ) then
        e = 1
      else
        e = 2
      end if

    end do

    a = triangle_node(e,t)

    lr = lrline ( x, y, node_xy(1,a), node_xy(2,a), node_xy(1,b), &
      node_xy(2,b), 0.0D+00 )

    if ( lr <= 0 ) then
      exit
    end if

  end do

  ltri = t
  ledg = e

  return
end
subroutine voronoi_plot ( box_requested, box_xy, circle_requested, &
  circle_x, circle_y, circle_r, file_name, dim_num, &
  node_num, coord, plot_min, plot_max, x, y, shadow, title )

!*****************************************************************************80
!
!! VORONOI_PLOT plots the Voronoi diagram of a set of points in 2D.
!
!  Discussion:
!
!    The routine first determines the Delaunay triangulation.
!
!    The Voronoi diagram is then determined from this information.
!
!    In particular, the circumcenter of each Delaunay triangle
!    is a vertex of a Voronoi polygon.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 October 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the output file.
!
!    Input, integer ( kind = 4 ) DIM_NUM, the number of spatial dimensions.
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of points.
!
!    Input, real ( kind = 8 ) COORD(DIM_NUM,NODE_NUM), the point coordinates.
!
!    Input, real ( kind = 8 ) PLOT_MIN(2), PLOT_MAX(2), the minimum and maximum
!    X and Y values to plot.
!
!    Input, integer ( kind = 4 ) X, Y, the indices of the two dimensions to plot.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) node_num

  real ( kind = 8 ) bitx
  real ( kind = 8 ) bity
  logical box_requested
  real ( kind = 8 ) box_xy(2,2)
  logical circle_requested
  real ( kind = 8 ) circle_r
  real ( kind = 8 ) circle_x
  real ( kind = 8 ) circle_y
  real ( kind = 8 ) coord(dim_num,node_num)
  real ( kind = 8 ) coord2(2,node_num)
  logical, parameter :: debug = .false.
  character ( len = * ) file_name
  integer ( kind = 4 ) file_unit
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) ierror
  logical inside
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jp1
  integer ( kind = 4 ) jp2
  integer ( kind = 4 ) k
  real ( kind = 8 ) n1
  real ( kind = 8 ) n2
  real ( kind = 8 ) plot_max(2)
  real ( kind = 8 ) plot_min(2)
  logical shadow
  real ( kind = 8 ) t(2,3)
  character ( len = * ) title
  real ( kind = 8 ) triangle_center(dim_num,2*node_num)
  integer ( kind = 4 ) triangle_neighbor(3,2*node_num)
  integer ( kind = 4 ) triangle_node(3,2*node_num)
  integer ( kind = 4 ) triangle_num
  integer ( kind = 4 ) x
  real ( kind = 8 ) xvec(4)
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ) x4
  real ( kind = 8 ) x5
  real ( kind = 8 ) x6
  integer ( kind = 4 ) y
  real ( kind = 8 ) yvec(4)
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) y3
  real ( kind = 8 ) y4
  real ( kind = 8 ) y5
  real ( kind = 8 ) y6
!
!  Compute the Delaunay triangulation.
!
  coord2(1,1:node_num) = coord(x,1:node_num)
  coord2(2,1:node_num) = coord(y,1:node_num)

  call dtris2 ( node_num, coord2, triangle_num, triangle_node, &
    triangle_neighbor )

  if ( debug ) then
    call r8mat_transpose_print ( 2, node_num, coord2, '  Node coordinates:' )
  end if

  if ( debug ) then
    call i4mat_transpose_print ( 3, triangle_num, triangle_node, '  Triangles:' )
  end if

  if ( debug ) then
    call i4mat_transpose_print ( 3, triangle_num, triangle_neighbor, &
      '  Neighbors:' )
  end if
!
!  For each Delaunay triangle, compute the intersection point of the
!  perpendicular bisectors of the sides.  
!
!  This is a vertex of a Voronoi polygon.
!
  do j = 1, triangle_num

    i1 = triangle_node(1,j)
    i2 = triangle_node(2,j)
    i3 = triangle_node(3,j)

    t(1:2,1) = coord2(1:2,i1)
    t(1:2,2) = coord2(1:2,i2)
    t(1:2,3) = coord2(1:2,i3)

    call triangle_circumcenter_2d ( t, triangle_center(1:2,j) )

  end do

  if ( debug ) then
    call r8mat_transpose_print ( 2, triangle_num, triangle_center, &
      '  Circumcenters:' )
  end if
!
!  Plot the Voronoi tessellation.
!
  call get_unit ( file_unit )

  call ps_file_open ( file_name, file_unit, ierror )

  call eps_file_head ( file_name, 36, 36, 576, 756  )

  call ps_page_head ( plot_min(1), plot_min(2), plot_max(1), plot_max(2) )
!
!  Initialize the line width to 1.
!
  call ps_line_width ( 1 )
!
!  Print title, if requested.
!
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
!
!  Shadow the points, if requested.
!
  if ( shadow ) then

    bitx = 0.025D+00 * ( plot_max ( 1 ) - plot_min ( 1 ) )
    bity = 0.025D+00 * ( plot_max ( 2 ) - plot_min ( 2 ) )

    do i = 1, node_num
      call ps_line ( plot_min(1), coord(y,i), plot_min(1) + bitx, coord(y,i) )
      call ps_line ( coord(x,i), plot_min(2), coord(x,i), plot_min(2) + bity )
    end do

  end if
!
!  Draw a user-specified box, if requested.
!
  if ( box_requested ) then
    call ps_comment ( 'User-requested box:' )
    call ps_line ( box_xy(1,1), box_xy(2,1), box_xy(1,2), box_xy(2,1) )
    call ps_line ( box_xy(1,2), box_xy(2,1), box_xy(1,2), box_xy(2,2) )
    call ps_line ( box_xy(1,2), box_xy(2,2), box_xy(1,1), box_xy(2,2) )
    call ps_line ( box_xy(1,1), box_xy(2,2), box_xy(1,1), box_xy(2,1) )
  end if
!
!  Draw a user-specified circle, if requested.
!
  if ( circle_requested ) then
    call ps_comment ( 'User-requested circle:' )
    call ps_circle ( circle_x, circle_y, circle_r )
  end if

  x1 = plot_min(1)
  y1 = plot_min(2)
  call ps_moveto ( x1, y1 )

! call ps_label ( 'Voronoi tessellation' )

  if ( node_num <= 350 ) then
    do i = 1, node_num
      call ps_mark_disk ( coord2(1,i), coord2(2,i) )
    end do
  else
    do i = 1, node_num
      call ps_mark_point ( coord2(1,i), coord2(2,i) )
    end do
  end if

  if ( node_num <= 350 ) then
    do i = 1, triangle_num
      call ps_mark_circle ( triangle_center(1,i), triangle_center(2,i) )
    end do
  end if
!
!  Increase the line width for smaller problems.
!
  if ( triangle_num < 50 ) then
    call ps_line_width ( 3 )
  else if ( triangle_num < 250 ) then
    call ps_line_width ( 2 )
  else
    call ps_line_width ( 1 )
  end if
!
!  For each Delaunay triangle, I
!    For each side J,
!
  do i = 1, triangle_num
    do j = 1, 3
      k = triangle_neighbor(j,i)
!
!  If there is a neighboring triangle K on that side,
!  connect the circumcenters of triangles I and K, creating 
!  one edge of the Voronoi polygon.
!
      if ( 0 < k  ) then
        if ( i < k ) then
          call ps_line ( triangle_center(1,i), triangle_center(2,i), &
            triangle_center(1,k), triangle_center(2,k) )
        end if
!
!  If there is no neighboring triangle on that side, then
!  if the circumcenter is outside the plot box, there is nothing to do.
!
      else if ( triangle_center(1,i) <= plot_min(1) .or. &
                plot_max(1) <= triangle_center(1,i) .or. &
                triangle_center(2,i) <= plot_min(2) .or. &
                plot_max(2) <= triangle_center(2,i) ) then
!
!  But if the circumcenter is inside the plot box, then
!  extend a line from the circumcenter of I in the direction of the
!  outward normal to that side.  This is a semi-infinite edge of
!  an infinite Voronoi polygon.
!
      else 

        ix = triangle_node(j,i)
        x1 = coord2(1,ix)
        y1 = coord2(2,ix)
        write ( *, * ) ' '
        write ( *, * ) 'Triangle I = ', i
        write ( *, * ) '  Local node index = ', j
        write ( *, * ) '  Global node index = ', ix
        write ( *, * ) 'X1,Y1=', x1, y1

        jp1 = i4_wrap ( j+1, 1, 3 )

        ix = triangle_node(jp1,i)
        x2 = coord2(1,ix)
        y2 = coord2(2,ix)
        write ( *, * ) '  Local node index = ', jp1
        write ( *, * ) '  Global node index = ', ix
        write ( *, * ) 'X2,Y2=', x2, y2

        x4 = triangle_center(1,i)
        y4 = triangle_center(2,i)

        call line_exp_normal_2d ( x1, y1, x2, y2, n1, n2 )

        x5 = x4 + n1
        y5 = y4 + n2

        write ( *, * ) 'Circumcenter = ', x4, y4
        write ( *, * ) 'Side perp = ', n1, n2
        write ( *, * ) 'New point on ray = ', x5, y5

        call box_ray_int_2d ( plot_min(1), plot_min(2), plot_max(1), &
          plot_max(2), x4, y4, x5, y5, x6, y6 )

        call ps_line ( x4, y4, x6, y6 )

        write ( *, * ) 'Intersects boundary at = ', x6, y6

      end if

    end do
  end do

  call ps_page_tail ( )

  call eps_file_tail ( )

  call ps_file_close ( file_unit )

  return
end
