subroutine triangle_mask ( dim_num, triangle_order, nodes, coord, mask )

!*****************************************************************************80
!
!! TRIANGLE_MASK is a user routine which masks triangles.
!
!  Discussion:
!
!    The region to be considered is the union of two rectangles.
!    The first is  -8 <= X <= 2, -1 <= Y <= 0,
!    the second is -2 <= X <= 8,  0 <= Y <= 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) TRIANGLE_ORDER, the number of nodes in the triangle.
!
!    Input, integer ( kind = 4 ) NODES(TRIANGLE_ORDER), the indices of the nodes.
!
!    Input, real ( kind = 8 ) COORD(DIM_NUM,TRIANGLE_ORDER), the coordinates
!    of the nodes.
!
!    Output, logical MASK, is TRUE if the triangle should be discarded,
!    and FALSE if the triangle should be retained.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) triangle_order

  real ( kind = 8 ) centroid(dim_num)
  real ( kind = 8 ) coord(dim_num,triangle_order)
  integer ( kind = 4 ) dim
  logical mask
  integer ( kind = 4 ) nodes(triangle_order)
!
!  Compute the centroid.
!
  do dim = 1, dim_num
    centroid(dim) = sum ( coord(dim,1:triangle_order) ) &
                  / real ( triangle_order, kind = 8 )
  end do
!
!  MASK = The centroid is outside the region.
!
  if (      -8.0D+00 <= centroid(1)            .and. &
                        centroid(1) <= 2.0D+00 .and. &
            -1.0D+00 <= centroid(2)            .and. &
                        centroid(2) <= 0.0D+00 ) then

    mask = .false.

  else if ( -2.0D+00 <= centroid(1)            .and. &
                        centroid(1) <= 8.0D+00 .and. &
             0.0D+00 <= centroid(2)            .and. &
                        centroid(2) <= 1.0D+00 ) then

    mask = .false.

  else

    mask = .true.

  end if

  return
end
