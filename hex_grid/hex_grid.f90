subroutine box_print_2d ( box )

!*****************************************************************************80
!
!! BOX_PRINT_2D prints information about a coordinate box in 2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) BOX(2,2), the coordinates of the lower left
!    and upper right corners of the box.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) box(dim_num,2)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Coordinate box:'
  write ( *, '(a,f8.4,a,f8.4)' ) '  ', box(1,1), ' <= X <= ', box(1,2)
  write ( *, '(a,f8.4,a,f8.4)' ) '  ', box(2,1), ' <= Y <= ', box(2,2)

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
subroutine hex_grid_01_approximate_h ( h_goal, nodes_per_layer, h )

!*****************************************************************************80
!
!! HEX_GRID_01_APPROXIMATE_H seeks a unit square hex grid with spacing H.
!
!  Discussion:
!
!    The parameter NODES_PER_LAYER controls the number of nodes and the
!    grid spacing, but in a somewhat obscure way.  This routine experiments
!    with various values until it is convinced it has the value 
!    of NODES_PER_LAYER that produces a grid spacing that is no
!    no greater than H.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) H_GOAL, the desired grid spacing.
!
!    Output, integer ( kind = 4 ) NODES_PER_LAYER, the number of nodes per layer
!    which produces a mesh with grid spacing H_GOAL or less.
!
!    Output, real ( kind = 8 ) H, the actual grid spacing.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) h_goal
  real ( kind = 8 ) h
  integer ( kind = 4 ) nodes_per_layer
  integer ( kind = 4 ) nodes_per_layer2

  if ( h_goal <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HEX_GRID_01_APPROXIMATE_H - Fatal error!'
    write ( *, '(a,g14.6)' ) '  Illegal input value of H_GOAL = ', h_goal
    stop
  end if

  nodes_per_layer = 1 + int ( 1.0D+00 / h_goal )
!
!  Check whether roundoff means we could use one less node per layer.
!
  if ( 2 < nodes_per_layer ) then

    nodes_per_layer2 = nodes_per_layer - 1
    h = 1.0D+00 / real ( nodes_per_layer2 - 1, kind = 8 )

    if ( h <= h_goal ) then
      nodes_per_layer = nodes_per_layer2
      return
    end if

  end if

  h = 1.0D+00 / real ( nodes_per_layer - 1, kind = 8 )

  return
end
subroutine hex_grid_01_approximate_n ( n_goal, nodes_per_layer, n )

!*****************************************************************************80
!
!! HEX_GRID_01_APPROXIMATE_N seeks a unit square hex grid of about N nodes.
!
!  Discussion:
!
!    The parameter NODES_PER_LAYER controls the number of nodes, but
!    in a somewhat obscure way.  This routine experiments with various
!    values until it is convinced it has the value of NODES_PER_LAYER
!    that comes as close as possible to producing N nodes.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N_GOAL, the desired number of nodes.
!
!    Output, integer ( kind = 4 ) NODES_PER_LAYER, the number of nodes per layer
!    which produces a mesh with about N_GOAL nodes.
!
!    Output, integer ( kind = 4 ) N, the number of nodes in the mesh.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_goal
  integer ( kind = 4 ) n_high
  integer ( kind = 4 ) n_low
  integer ( kind = 4 ) nodes_per_layer
  integer ( kind = 4 ) nodes_per_layer_high
  integer ( kind = 4 ) nodes_per_layer_low

  if ( n_goal <= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HEX_GRID_01_APPROXIMATE_N - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal input value of N_GOAL = ', n_goal
    stop
  end if

  nodes_per_layer_low = 0
  n_low = 0

  nodes_per_layer = nint ( sqrt ( real ( n_goal, kind = 8 ) ) )

  nodes_per_layer_high = n_goal 
  n_high = n_goal**2

  do 

    call hex_grid_01_n ( nodes_per_layer, n )

    if ( n == n_goal ) then
      exit
    end if

    if ( n < n_goal ) then
      nodes_per_layer_low = nodes_per_layer
      n_low = n
    else
      nodes_per_layer_high = nodes_per_layer
      n_high = n
    end if

    if ( nodes_per_layer_low + 1 == nodes_per_layer_high ) then
      if ( n - n_low <= n_high - n ) then
        nodes_per_layer = nodes_per_layer_high
        n = n_high
      else
        nodes_per_layer = nodes_per_layer_low
        n = n_low
      end if
      exit
    end if

    nodes_per_layer = ( nodes_per_layer_low + nodes_per_layer_high ) / 2

  end do

  return
end
subroutine hex_grid_01_h ( nodes_per_layer, hx, hy )

!*****************************************************************************80
!
!! HEX_GRID_01_H computes the unit square hex grid spacings.
!
!  Discussion:
!
!    This routine determines the values of HX and HY from
!    the fundamental hexagonal grid parameter NODES_PER_LAYER.
!
!    A hexagonal grid is defined in the unit square [0,1] x [0,1].
!
!    All nodes of the grid lie on one of LAYERS horizontal lines.
!    The first of these lines is the X axis, and each successive
!    line is HY units higher.
!
!    On all the odd numbered lines, there are NODES_PER_LAYER points, 
!    equally spaced from 0 to 1, with a spacing of HX.  
!
!    On the even numbered lines, there are NODES_PER_LAYER-1 points,
!    whose values are the midpoints of successive intervals on
!    an odd numbered line.  (The grid is staggered).
!
!    HY = HX * sqrt ( 3 ) / 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODES_PER_LAYER, the number of grid points
!    on the first horizontal layer of points.
!
!    Output, real ( kind = 8 ) HX, the spacing between grid points
!    on a horizontal line.
!
!    Output, real ( kind = 8 ) HY, the spacing between horizontal lines.
!
  implicit none

  real ( kind = 8 ) hx
  real ( kind = 8 ) hy
  integer ( kind = 4 ) nodes_per_layer

  if ( nodes_per_layer < 1 ) then

    hx = 0.0D+00
    hy = 0.0D+00

  else if ( nodes_per_layer == 1 ) then

    hx = 1.0D+00
    hy = 1.0D+00

  else

    hx = 1.0D+00 / real ( nodes_per_layer - 1, kind = 8 )
    hy = hx * sqrt ( 3.0D+00 ) / 2.0D+00

  end if

  return
end
subroutine hex_grid_01_layers ( nodes_per_layer, layers )

!*****************************************************************************80
!
!! HEX_GRID_01_LAYERS computes the unit square hex grid column width.
!
!  Discussion:
!
!    This routine determines the value of LAYERS, the number of
!    layers, from the fundamental hexagonal grid parameter NODES_PER_LAYER.
!
!    A hexagonal grid is defined in the unit square [0,1] x [0,1].
!
!    All nodes of the grid lie on one of LAYERS horizontal lines.
!    The first of these lines is the X axis, and each successive
!    line is HY units higher.
!
!    On all the odd numbered lines, there are NODES_PER_LAYER points, 
!    equally spaced from 0 to 1, with a spacing of HX.  
!
!    On the even numbered lines, there are NODES_PER_LAYER-1 points,
!    whose values are the midpoints of successive intervals on
!    an odd numbered line.  (The grid is staggered).
!
!    HY = HX * sqrt ( 3 ) / 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODES_PER_LAYER, the number of grid points 
!    on the first horizontal layer of points.
!
!    Output, integer ( kind = 4 ) LAYERS, the number of horizontal layers.
!
  implicit none

  real ( kind = 8 ) hx
  real ( kind = 8 ) hy
  integer ( kind = 4 ) layers
  integer ( kind = 4 ) nodes_per_layer

  if ( nodes_per_layer < 1 ) then

    layers = 0
    
  else if ( nodes_per_layer == 1 ) then

    layers = 1

  else

    hx = 1.0D+00 / real ( nodes_per_layer - 1, kind = 8 )
    hy = sqrt ( 3.0D+00 ) * hx / 2.0D+00
    layers = 1 + int ( 1.0D+00 / hy )

  end if

  return
end
subroutine hex_grid_01_n ( nodes_per_layer, n )

!*****************************************************************************80
!
!! HEX_GRID_01_N computes the number of unit square hex grid points.
!
!  Discussion:
!
!    This routine determines the value of N, the number of
!    hex grid points, from the fundamental hexagonal grid 
!    parameter NODES_PER_LAYER.
!
!    A hexagonal grid is defined in the unit square [0,1] x [0,1].
!
!    All nodes of the grid lie on one of LAYERS horizontal lines.
!    The first of these lines is the X axis, and each successive
!    line is HY units higher.
!
!    On all the odd numbered lines, there are NODES_PER_LAYER points, 
!    equally spaced from 0 to 1, with a spacing of HX.  
!
!    On the even numbered lines, there are NODES_PER_LAYER-1 points,
!    whose values are the midpoints of successive intervals on
!    an odd numbered line.  (The grid is staggered).
!
!    HY = HX * sqrt ( 3 ) / 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODES_PER_LAYER, the number of grid points 
!    on the first horizontal layer of points.
!
!    Output, integer ( kind = 4 ) N, the number of hex grid points.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) layers
  integer ( kind = 4 ) nodes_per_layer

  if ( nodes_per_layer < 1 ) then

    n = 0

  else if ( nodes_per_layer == 1 ) then

    n = 1

  else

    call hex_grid_01_layers ( nodes_per_layer, layers )

    n = nodes_per_layer       * ( ( layers + 1 ) / 2 ) + &
      ( nodes_per_layer - 1 ) * ( ( layers     ) / 2 )

  end if

  return
end
subroutine hex_grid_01_points ( nodes_per_layer, layers, n, p )

!*****************************************************************************80
!
!! HEX_GRID_01_POINTS returns unit square hex grid points.
!
!  Discussion:
!
!    This routine determines the coordinates of the elements of
!    a hexagonal grid in the unit square.
!
!    A hexagonal grid is defined in the unit square [0,1] x [0,1].
!
!    All nodes of the grid lie on one of LAYERS horizontal lines.
!    The first of these lines is the X axis, and each successive
!    line is HY units higher.
!
!    On all the odd numbered lines, there are NODES_PER_LAYER points, 
!    equally spaced from 0 to 1, with a spacing of HX.  
!
!    On the even numbered lines, there are NODES_PER_LAYER-1 points,
!    whose values are the midpoints of successive intervals on
!    an odd numbered line.  (The grid is staggered).
!
!    HY = HX * sqrt ( 3 ) / 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODES_PER_LAYER, the number of grid points on 
!    the first horizontal layer of points.
!
!    Input, integer ( kind = 4 ) LAYERS, the number of horizontal layers.
!
!    Input, integer ( kind = 4 ) N, the total number of hex grid points.
!
!    Output, real (kind = 8 ) P(2,N), the coordinates of the
!    mesh points, listed one horizontal layer at a time.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) hx
  real ( kind = 8 ) hy
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jmod
  integer ( kind = 4 ) k
  integer ( kind = 4 ) layers
  integer ( kind = 4 ) nodes_per_layer
  real ( kind = 8 ) p(dim_num,n)
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  if ( nodes_per_layer < 1 ) then
    return
  end if

  if ( nodes_per_layer == 1 ) then
    p(1:dim_num,1) = 0.5D+00
    return
  end if

  call hex_grid_01_h ( nodes_per_layer, hx, hy )

  k = 0

  do j = 1, layers

    y = hy * real ( j - 1, kind = 8 )

    jmod = mod ( j, 2 )

    if ( jmod == 1 ) then

      do i = 1, nodes_per_layer
        x = real ( i - 1, kind = 8 ) / real ( nodes_per_layer - 1, kind = 8 )
        k = k + 1
        if ( k <= n ) then
          p(1,k) = x
          p(2,k) = y
        end if
      end do

    else

      do i = 1, nodes_per_layer-1
        x = real ( 2 * i - 1, kind = 8 ) &
          / real ( 2 * nodes_per_layer - 2, kind = 8 )
        k = k + 1
        if ( k <= n ) then
          p(1,k) = x
          p(2,k) = y
        end if
      end do

    end if

  end do

  return
end
subroutine hex_grid_approximate_h ( box, h_goal, nodes_per_layer, h )

!*****************************************************************************80
!
!! HEX_GRID_APPROXIMATE_H seeks a hex grid with spacing H.
!
!  Discussion:
!
!    The parameter NODES_PER_LAYER controls the number of nodes and the
!    grid spacing, but in a somewhat obscure way.  This routine experiments
!    with various values until it is convinced it has the value 
!    of NODES_PER_LAYER that produces a grid spacing that is no
!    no greater than H.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) BOX(2,2), the lower and upper corners
!    of the rectangular region.
!
!    Input, real ( kind = 8 ) H_GOAL, the desired grid spacing.
!
!    Output, integer ( kind = 4 ) NODES_PER_LAYER, the number of nodes per layer
!    which produces a mesh with grid spacing H_GOAL or less.
!
!    Output, real ( kind = 8 ) H, the actual grid spacing.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) box(dim_num,2)
  real ( kind = 8 ) h_goal
  real ( kind = 8 ) h
  integer ( kind = 4 ) nodes_per_layer
  integer ( kind = 4 ) nodes_per_layer2

  if ( h_goal <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HEX_GRID_APPROXIMATE_H - Fatal error!'
    write ( *, '(a,g14.6)' ) '  Illegal input value of H_GOAL = ', h_goal
    stop
  end if

  nodes_per_layer = 1 + int ( ( box(1,2) - box(1,1) ) / h_goal )
!
!  Check whether roundoff means we could use one less node per layer.
!
  if ( 2 < nodes_per_layer ) then

    nodes_per_layer2 = nodes_per_layer - 1
    h = ( box(1,2) - box(1,1) ) / real ( nodes_per_layer2 - 1, kind = 8 )

    if ( h <= h_goal ) then
      nodes_per_layer = nodes_per_layer2
      return
    end if

  end if

  h = ( box(1,2) - box(1,1) ) / real ( nodes_per_layer - 1, kind = 8 )

  return
end
subroutine hex_grid_approximate_n ( box, n_goal, nodes_per_layer, n )

!*****************************************************************************80
!
!! HEX_GRID_APPROXIMATE_N seeks a hex grid of about N nodes.
!
!  Discussion:
!
!    The parameter NODES_PER_LAYER controls the number of nodes, but
!    in a somewhat obscure way.  This routine experiments with various
!    values until it is convinced it has the value of NODES_PER_LAYER
!    that comes as close as possible to producing N nodes.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) BOX(2,2), the lower and upper corners
!    of the rectangular region.
!
!    Input, integer ( kind = 4 ) N_GOAL, the desired number of nodes.
!
!    Output, integer ( kind = 4 ) NODES_PER_LAYER, the number of nodes per layer
!    which produces a mesh with about N_GOAL nodes.
!
!    Output, integer ( kind = 4 ) N, the number of nodes in the mesh.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) box(dim_num,2)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_goal
  integer ( kind = 4 ) n_high
  integer ( kind = 4 ) n_low
  integer ( kind = 4 ) nodes_per_layer
  integer ( kind = 4 ) nodes_per_layer_high
  integer ( kind = 4 ) nodes_per_layer_low

  if ( n_goal <= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HEX_GRID_APPROXIMATE_N - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal input value of N_GOAL = ', n_goal
    stop
  end if

  nodes_per_layer_low = 0
  n_low = 0

  nodes_per_layer = nint ( sqrt ( real ( n_goal, kind = 8 ) ) )

  nodes_per_layer_high = n_goal  
  n_high = n_goal**2

  do 

    call hex_grid_n ( nodes_per_layer, box, n )

    if ( n == n_goal ) then
      exit
    end if

    if ( n < n_goal ) then
      nodes_per_layer_low = nodes_per_layer
      n_low = n
    else
      nodes_per_layer_high = nodes_per_layer
      n_high = n
    end if

    if ( nodes_per_layer_low + 1 == nodes_per_layer_high ) then
      if ( n - n_low <= n_high - n ) then
        nodes_per_layer = nodes_per_layer_high
        n = n_high
      else
        nodes_per_layer = nodes_per_layer_low
        n = n_low
      end if
      exit
    end if

    nodes_per_layer = ( nodes_per_layer_low + nodes_per_layer_high ) / 2

  end do

  return
end
subroutine hex_grid_h ( nodes_per_layer, box, hx, hy )

!*****************************************************************************80
!
!! HEX_GRID_H computes the coordinate box hex grid spacings.
!
!  Discussion:
!
!    This routine determines the values of HX and HY from
!    the fundamental hexagonal grid parameter NODES_PER_LAYER.
!
!    A hexagonal grid is defined in the coordinate box [A,B] x [C,D].
!
!    All nodes of the grid lie on one of LAYERS horizontal lines.
!    The first of these lines is from (A,C) to (B,C), and each 
!    successive line is HY units higher.
!
!    On all the odd numbered lines, there are NODES_PER_LAYER points, 
!    equally spaced from A to B, with a spacing of HX.  
!
!    On the even numbered lines, there are NODES_PER_LAYER-1 points,
!    whose values are the midpoints of successive intervals on
!    an odd numbered line.  (The grid is staggered).
!
!    HY = HX * sqrt ( 3 ) / 2.
!
!  Modified:
!
!    02 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODES_PER_LAYER, the number of grid points 
!    on the first horizontal layer of points.
!
!    Input, real ( kind = 8 ) BOX(2,2), the values of A, B, C and D
!    that define the coordinate box.
!
!    Output, real ( kind = 8 ) HX, the spacing between grid points
!    on a horizontal line.
!
!    Output, real ( kind = 8 ) HY, the spacing between horizontal lines.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) box(dim_num,2)
  real ( kind = 8 ) hx
  real ( kind = 8 ) hy
  integer ( kind = 4 ) nodes_per_layer

  if ( nodes_per_layer < 1 ) then

    hx = 0.0D+00
    hy = 0.0D+00

  else if ( nodes_per_layer == 1 ) then

    hx = box(1,2) - box(1,1)
    hy = box(2,2) - box(2,1)

  else

    hx = ( box(1,2) - box(1,1) ) / real ( nodes_per_layer - 1, kind = 8 )
    hy = hx * sqrt ( 3.0D+00 ) / 2.0D+00

  end if

  return
end
subroutine hex_grid_layers ( nodes_per_layer, box, layers )

!*****************************************************************************80
!
!! HEX_GRID_LAYERS computes the coordinate box hex grid column width.
!
!  Discussion:
!
!    This routine determines the value of LAYERS, the number of
!    layers, from the fundamental hexagonal grid parameter NODES_PER_LAYER.
!
!    A hexagonal grid is defined in a coordinate box [A,B] x [C,D].
!
!    All nodes of the grid lie on one of LAYERS horizontal lines.
!    The first of these lines is from (A,C) to (B,C), and each 
!    successive line is HY units higher.
!
!    On all the odd numbered lines, there are NODES_PER_LAYER points, 
!    equally spaced from A to B, with a spacing of HX.  
!
!    On the even numbered lines, there are NODES_PER_LAYER-1 points,
!    whose values are the midpoints of successive intervals on
!    an odd numbered line.  (The grid is staggered).
!
!    HY = HX * sqrt ( 3 ) / 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODES_PER_LAYER, the number of grid points 
!    on the first horizontal layer of points.
!
!    Input, real ( kind = 8 ) BOX(2,2), the values of A, B, C and D
!    that define the coordinate box.
!
!    Output, integer ( kind = 4 ) LAYERS, the number of horizontal layers.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) box(dim_num,2)
  real ( kind = 8 ) hx
  real ( kind = 8 ) hy
  integer ( kind = 4 ) layers
  integer ( kind = 4 ) nodes_per_layer

  if ( nodes_per_layer < 1 ) then

    layers = 0
    
  else if ( nodes_per_layer == 1 ) then

    layers = 1

  else

    hx = ( box(1,2) - box(1,1) ) / real ( nodes_per_layer - 1, kind = 8 )
    hy = sqrt ( 3.0D+00 ) * hx / 2.0D+00
    layers = 1 + int ( ( box(2,2) - box(2,1) ) / hy )

  end if

  return
end
subroutine hex_grid_n ( nodes_per_layer, box, n )

!*****************************************************************************80
!
!! HEX_GRID_N computes the number of coordinate box hex grid points.
!
!  Discussion:
!
!    This routine determines the value of N, the number of
!    hex grid points, from the fundamental hexagonal grid 
!    parameter NODES_PER_LAYER.
!
!    A hexagonal grid is defined in the coordinate box [A,B] x [C,D].
!
!    All nodes of the grid lie on one of LAYERS horizontal lines.
!    The first of these lines is from (A,C) to (B,C), and each 
!    successive line is HY units higher.
!
!    On all the odd numbered lines, there are NODES_PER_LAYER points, 
!    equally spaced from A to B, with a spacing of HX.  
!
!    On the even numbered lines, there are NODES_PER_LAYER-1 points,
!    whose values are the midpoints of successive intervals on
!    an odd numbered line.  (The grid is staggered).
!
!    HY = HX * sqrt ( 3 ) / 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODES_PER_LAYER, the number of grid points 
!    on the first horizontal layer of points.
!
!    Input, real ( kind = 8 ) BOX(2,2), the values of A, B, C and D
!    that define the coordinate box.
!
!    Output, integer ( kind = 4 ) N, the number of hex grid points.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) box(dim_num,2)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) layers
  integer ( kind = 4 ) nodes_per_layer

  if ( nodes_per_layer < 1 ) then

    n = 0

  else if ( nodes_per_layer == 1 ) then

    n = 1

  else

    call hex_grid_layers ( nodes_per_layer, box, layers )

    n = nodes_per_layer       * ( ( layers + 1 ) / 2 ) + &
      ( nodes_per_layer - 1 ) * ( ( layers     ) / 2 )

  end if

  return
end
subroutine hex_grid_points ( nodes_per_layer, layers, n, box, p )

!*****************************************************************************80
!
!! HEX_GRID_POINTS returns coordinate box hex grid points.
!
!  Discussion:
!
!    This routine determines the coordinates of the elements of
!    a hexagonal grid in the unit square.
!
!    A hexagonal grid is defined in the coordinate box [A,B] x [C,D].
!
!    All nodes of the grid lie on one of LAYERS horizontal lines.
!    The first of these lines is from (A,C) to (B,C), and each 
!    successive line is HY units higher.
!
!    On all the odd numbered lines, there are NODES_PER_LAYER points, 
!    equally spaced from A to B, with a spacing of HX.  
!
!    On the even numbered lines, there are NODES_PER_LAYER-1 points,
!    whose values are the midpoints of successive intervals on
!    an odd numbered line.  (The grid is staggered).
!
!    HY = HX * sqrt ( 3 ) / 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODES_PER_LAYER, the number of grid points on 
!    the first horizontal layer of points.
!
!    Input, integer ( kind = 4 ) LAYERS, the number of horizontal layers.
!
!    Input, integer ( kind = 4 ) N, the total number of hex grid points.
!
!    Input, real ( kind = 8 ) BOX(2,2), the values of A, B, C and D
!    that define the coordinate box.
!
!    Output, real (kind = 8 ) P(2,N), the coordinates of the
!    mesh points, listed one horizontal layer at a time.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) box(dim_num,2)
  real ( kind = 8 ) hx
  real ( kind = 8 ) hy
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jmod
  integer ( kind = 4 ) k
  integer ( kind = 4 ) layers
  integer ( kind = 4 ) nodes_per_layer
  real ( kind = 8 ) p(dim_num,n)
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  if ( nodes_per_layer < 1 ) then
    return
  end if

  if ( nodes_per_layer == 1 ) then
    do i = 1, dim_num
      p(i,1) = ( box(i,1) + box(i,2) ) / 2.0D+00
    end do
    return
  end if

  call hex_grid_h ( nodes_per_layer, box, hx, hy )

  k = 0

  do j = 1, layers

    y = box(2,1) + hy * real ( j - 1, kind = 8 )

    jmod = mod ( j, 2 )

    if ( jmod == 1 ) then

      do i = 1, nodes_per_layer
        x = box(1,1) + ( box(1,2) - box(1,1) ) &
          * real ( i - 1, kind = 8 ) &
          / real ( nodes_per_layer - 1, kind = 8 )
        k = k + 1
        if ( k <= n ) then
          p(1,k) = x
          p(2,k) = y
        end if
      end do

    else

      do i = 1, nodes_per_layer-1
        x = box(1,1) + ( box(1,2) - box(1,1) ) &
          * real ( 2 * i - 1, kind = 8 ) &
          / real ( 2 * nodes_per_layer - 2, kind = 8 )
        k = k + 1
        if ( k <= n ) then
          p(1,k) = x
          p(2,k) = y
        end if
      end do

    end if

  end do

  return
end
subroutine hex_grid_write ( n, nodes_per_layer, layers, hx, hy, box, r, &
  file_out_name )

!*****************************************************************************80
!
!! HEX_GRID_WRITE writes a hex grid dataset to a file.
!
!  Discussion:
!
!    The initial lines of the file are comments, which begin with a
!    '#' character.
!
!    Thereafter, each line of the file contains the M-dimensional
!    components of the next entry of the sequence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, integer ( kind = 4 ) NODES_PER_LAYER, the number of nodes in the 
!    first layer.
!
!    Input, integer ( kind = 4 ) LAYERS, the number of layers.
!
!    Input, real ( kind = 8 ) HX, HY, the row and column spacings.
!
!    Input, real ( kind = 8 ) BOX(2,2), the values of A, B, C and D
!    that define the coordinate box.
!
!    Input, real ( kind = 8 ) R(M,N), the points.
!
!    Input, character ( len = * ) FILE_OUT_NAME, the name of
!    the output file.
!
  implicit none

  integer   ( kind = 4 ) n
  integer   ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) box(dim_num,2)
  character ( len = * ) file_out_name
  integer ( kind = 4 ) file_out_unit
  real ( kind = 8 ) hx
  real ( kind = 8 ) hy
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) j
  integer ( kind = 4 ) layers
  integer ( kind = 4 ) nodes_per_layer
  real ( kind = 8 ) r(dim_num,n)
  character ( len = 40 ) string

  call get_unit ( file_out_unit )

  open ( unit = file_out_unit, file = file_out_name, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HEX_GRID_WRITE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output file.'
    stop
  end if

  write ( file_out_unit, '(a)'       ) '#  ' // trim ( file_out_name )
  write ( file_out_unit, '(a)'       ) '#  created by HEX_GRID_WRITE.F90.'
  write ( file_out_unit, '(a)'       ) '#'
  write ( file_out_unit, '(a)'       ) '#'
  write ( file_out_unit, '(a,i8)'    ) &
    '#  Spatial dimension DIM_NUM =     ', dim_num
  write ( file_out_unit, '(a,i8)'    ) '#  Number of points N =      ', n
  write ( file_out_unit, '(a,i12)'   ) &
    '#  NODES_PER_LAYER =         ', nodes_per_layer
  write ( file_out_unit, '(a,i12)'   ) '#  Number of rows LAYERS =   ', layers
  write ( file_out_unit, '(a,g14.6)' ) '#  Coordinate box X(1) =     ', box(1,1)
  write ( file_out_unit, '(a,g14.6)' ) '#  Coordinate box X(2) =     ', box(1,2)
  write ( file_out_unit, '(a,g14.6)' ) '#  Coordinate box Y(1) =     ', box(2,1)
  write ( file_out_unit, '(a,g14.6)' ) '#  Coordinate box Y(2) =     ', box(2,2)
  write ( file_out_unit, '(a,g14.6)' ) '#  Node spacing HX =         ', hx
  write ( file_out_unit, '(a,g14.6)' ) '#  Layer spacing HY =        ', hy
  write ( file_out_unit, '(a,g14.6)' ) '#  EPSILON (unit roundoff) = ', &
    epsilon ( r(1,1) )

  write ( file_out_unit, '(a)'       ) '#'

  write ( string, '(a,i3,a)' ) '(', dim_num, '(2x,f12.6))'
  do j = 1, n
    write ( file_out_unit, string ) r(1:dim_num,j)
  end do

  close ( unit = file_out_unit )

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
  integer ( kind = 4 ) get
  integer ( kind = 4 ) put
  integer ( kind = 4 ) nchar
  character ( len = * ) s
  character, parameter :: TAB = char ( 9 )

  put = 0
  nchar = len_trim ( s )

  do get = 1, nchar

    c = s(get:get)

    if ( c /= ' ' .and. c /= TAB ) then
      put = put + 1
      s(put:put) = c
    end if

  end do

  s(put+1:nchar) = ' '

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

  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
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
