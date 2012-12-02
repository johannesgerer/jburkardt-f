program main

!*****************************************************************************80
!
!! MAIN is the main program for HEX_GRID_TRIANGULATE.
!
!  Discussion:
!
!    HEX_GRID_TRIANGULATE uses HEX_GRID to find nodes in a region.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 March 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) test
  integer ( kind = 4 ) test_num

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HEX_GRID_TRIANGULATE'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Find hex grid nodes inside a test triangulation region.'

  call p00_test_num ( test_num )

  do test = 1, test_num
    call test01 ( test )
  end do

  do test = 1, test_num
    call test02 ( test )
  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HEX_GRID_TRiANGULATE'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( test )

!*****************************************************************************80
!
!! TEST01 computes a hexagonal grid in an unknown triangulation region.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 March 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer   ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) box(dim_num,2)
  character ( len = 80 ) file_name
  real ( kind = 8 ) hi(dim_num)
  real ( kind = 8 ) hx
  real ( kind = 8 ) hy
  integer   ( kind = 4 ) i
  logical, allocatable, dimension ( : ) :: inside
  integer   ( kind = 4 ) layers
  real ( kind = 8 ) lo(dim_num)
  integer   ( kind = 4 ) n
  integer   ( kind = 4 ) n2
  integer   ( kind = 4 ) nodes_per_layer
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: p
  integer   ( kind = 4 ) test
!
!  Set the number of nodes per layer.
!
  nodes_per_layer = 40

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Given an unknown triangulation region,'
  write ( *, '(a)' ) '  get its bounding box, compute parameters'
  write ( *, '(a)' ) '  for a hexagonal grid with horizontal count '
  write ( *, '(a,i6)' ) '  NODES_PER_LAYER = ', nodes_per_layer
  write ( *, '(a)' ) '  generate the full grid,'
  write ( *, '(a)' ) '  count those grid points that are in the region,'
  write ( *, '(a)' ) '  and write the grid points to a file.'
!
!  Get the bounding box coordinates.
!
  call p00_box ( test, dim_num, lo, hi )

  box(1:dim_num,1) = lo(1:dim_num)
  box(1:dim_num,2) = hi(1:dim_num)

  call box_print_2d ( box )
!
!  Determine the number of layers.
!
  call hex_grid_layers ( nodes_per_layer, box, layers )

  call hex_grid_h ( nodes_per_layer, box, hx, hy )

  call hex_grid_n ( nodes_per_layer, box, n )

  allocate ( p(1:dim_num,1:n) )
  allocate ( inside(1:n) )

  call hex_grid_points ( nodes_per_layer, layers, n, box, p )

  call p00_inside ( test, dim_num, n, p, inside )

  n2 = 0
  do i = 1, n

    if ( inside(i) ) then
      n2 = n2 + 1
      p(1:dim_num,n2) = p(1:dim_num,i)
    end if

  end do

  write ( file_name, '(a,i2,a,i3,a,i3,a,i6,a)' ) &
    'p_', test, '_', nodes_per_layer, '_', layers, '_', n2, '.txt'

  call s_blank_delete ( file_name )

  write ( *, '(2x,i6,2x,i3,2x,i6,4x,a)' ) &
    nodes_per_layer, layers, n2, trim ( file_name )

  call hex_grid_write ( n2, nodes_per_layer, layers, hx, hy, box, &
    p, file_name )

  deallocate ( p )
  deallocate ( inside )

  return
end
subroutine test02 ( test )

!*****************************************************************************80
!
!! TEST02 tries to compute a hexagonal grid of about 350 points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 March 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer   ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) box(dim_num,2)
  character ( len = 80 ) file_name
  real ( kind = 8 ) hi(dim_num)
  real ( kind = 8 ) hx
  real ( kind = 8 ) hy
  integer   ( kind = 4 ) i
  logical, allocatable, dimension ( : ) :: inside
  integer   ( kind = 4 ) layers
  real ( kind = 8 ) lo(dim_num)
  integer   ( kind = 4 ) n
  integer   ( kind = 4 ), parameter :: n_goal = 350
  integer   ( kind = 4 ) n_high
  integer   ( kind = 4 ) n_low
  integer   ( kind = 4 ) n2
  integer   ( kind = 4 ) nodes_per_layer
  integer   ( kind = 4 ) nodes_per_layer_high
  integer   ( kind = 4 ) nodes_per_layer_low
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: p
  integer   ( kind = 4 ) test
!
!  Set the number of nodes per layer.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  Given an unknown triangulation region,'
  write ( *, '(a)' ) '  get its bounding box, and try to find'
  write ( *, '(a)' ) '  the appropriate value of NODES_PER_LAYER'
  write ( *, '(a,i6,a)' ) '  so that we get about N = ', n_goal, ' nodes.'
!
!  Get the bounding box coordinates.
!
  call p00_box ( test, dim_num, lo, hi )

  box(1:dim_num,1) = lo(1:dim_num)
  box(1:dim_num,2) = hi(1:dim_num)

  call box_print_2d ( box )
!
!  This value we are sure of!
!
  nodes_per_layer_low = 0
  n_low = 0

  nodes_per_layer_high = n_goal
  n_high = n_goal**2

  do

    nodes_per_layer = ( nodes_per_layer_low + nodes_per_layer_high ) / 2

    call hex_grid_layers ( nodes_per_layer, box, layers )

    call hex_grid_n ( nodes_per_layer, box, n )

    call hex_grid_h ( nodes_per_layer, box, hx, hy )

    allocate ( p(1:dim_num,1:n) )
    allocate ( inside(1:n) )

    call hex_grid_points ( nodes_per_layer, layers, n, box, p )

    call p00_inside ( test, dim_num, n, p, inside )

    n2 = 0
    do i = 1, n

      if ( inside(i) ) then
        n2 = n2 + 1
        p(1:dim_num,n2) = p(1:dim_num,i)
      end if

    end do

    if ( n2 == n_goal ) then
      exit
    end if

    if ( n2 < n_goal ) then
      nodes_per_layer_low = nodes_per_layer
      n_low = n2
    else
      nodes_per_layer_high = nodes_per_layer
      n_high = n2
    end if

    if ( nodes_per_layer_low + 1 == nodes_per_layer_high ) then
      if ( n - n_low <= n_high - n ) then
        nodes_per_layer = nodes_per_layer_high
      else
        nodes_per_layer = nodes_per_layer_low
      end if
      exit
    end if

    deallocate ( p )
    deallocate ( inside )

  end do

  write ( file_name, '(a,i2,a,i3,a,i3,a,i6,a)' ) &
    'p_', test, '_', nodes_per_layer, '_', layers, '_', n2, '.txt'

  call s_blank_delete ( file_name )

  write ( *, '(2x,i6,2x,i3,2x,i6,4x,a)' ) &
    nodes_per_layer, layers, n2, trim ( file_name )

  call hex_grid_write ( n2, nodes_per_layer, layers, hx, hy, box, &
    p, file_name )

  deallocate ( p )
  deallocate ( inside )

  return
end
