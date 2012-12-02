subroutine triangle_mask ( dim_num, triangle_order, nodes, coord, mask )

!*****************************************************************************80
!
!! TRIANGLE_MASK is a user routine which masks triangles.
!
!  Discussion:
!
!    The region to be considered is the [0,4]x[0,4] square.
!
!    We want to remove the lower left triangular corner,
!    and part of the upper right triangular corner.
!
!    The following diagram of the 25 nodes indicates by "O" the
!    nodes that should end up being deleted, although the deletion
!    is actually done by triangles.
!
!    Before masking:
!
!      X - X - X - X - X
!      | \ | \ | \ | \ |
!      X - X - X - X - X
!      | \ | \ | \ | \ |
!      X - X - X - X - X
!      | \ | \ | \ | \ |
!      X - X - X - X - X
!      | \ | \ | \ | \ |
!      X - X - X - X - X
!
!    After masking:
!
!      X - X   O   O   O
!      | \ | \          
!      X - X - X   O   O
!      | \ | \ | \      
!      X - X - X - X - X
!        \ | \ | \ | \ |
!      O   X - X - X - X
!            \ | \ | \ |
!      O   O   X - X - X
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 January 2007
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
!  Remove the lower left corner
!
  if ( centroid(1) + centroid(2) < 2.0D+00 ) then

    mask = .true.
!
!  Remove the upper right section.
!
  else if ( 5.0D+00 < centroid(1) + centroid(2) .and. &
            2.0D+00 < centroid(2) ) then

    mask = .true.
!
!  Keep everything else.
!
  else

    mask = .false.

  end if

  return
end
