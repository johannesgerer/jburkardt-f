subroutine bandwidth_mesh ( element_order, element_num, element_node, &
  ml, mu, m )

!*****************************************************************************80
!
!! BANDWIDTH_MESH: bandwidth of finite element mesh.
!
!  Discussion:
!
!    The quantity computed here is the "geometric" bandwidth determined
!    by the finite element mesh alone.
!
!    If a single finite element variable is associated with each node
!    of the mesh, and if the nodes and variables are numbered in the
!    same way, then the geometric bandwidth is the same as the bandwidth
!    of a typical finite element matrix.
!
!    The bandwidth M is defined in terms of the lower and upper bandwidths:
!
!      M = ML + 1 + MU
!
!    where 
!
!      ML = maximum distance from any diagonal entry to a nonzero
!      entry in the same row, but earlier column,
!
!      MU = maximum distance from any diagonal entry to a nonzero
!      entry in the same row, but later column.
!
!    Because the finite element node adjacency relationship is symmetric,
!    we are guaranteed that ML = MU.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 September 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ELEMENT_ORDER, the order of the elements.
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of elements.
!
!    Input, integer ( kind = 4 ) ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM); 
!    ELEMENT_NODE(I,J) is the global index of local node I in element J.
!
!    Output, integer ( kind = 4 ) ML, MU, the lower and upper bandwidths of 
!    the matrix.
!
!    Output, integer ( kind = 4 ) M, the bandwidth of the matrix.
!
  implicit none

  integer ( kind = 4 ) element_num
  integer ( kind = 4 ) element_order

  integer ( kind = 4 ) element
  integer ( kind = 4 ) element_node(element_order,element_num)
  integer ( kind = 4 ) global_i
  integer ( kind = 4 ) global_j
  integer ( kind = 4 ) local_i
  integer ( kind = 4 ) local_j
  integer ( kind = 4 ) m
  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu

  ml = 0
  mu = 0
  
  do element = 1, element_num

    do local_i = 1, element_order
      global_i = element_node(local_i,element)

      do local_j = 1, element_order
        global_j = element_node(local_j,element)

        mu = max ( mu, global_j - global_i )
        ml = max ( ml, global_i - global_j )

      end do
    end do
  end do

  m = ml + 1 + mu

  return
end
subroutine bandwidth_var ( element_order, element_num, element_node, &
  node_num, var_node, var_num, var, ml, mu, m )

!*****************************************************************************80
!
!! BANDWIDTH_VAR determines the bandwidth for finite element variables.
!
!  Discussion:
!
!    We assume that, attached to each node in the finite element mesh
!    there are a (possibly zero) number of finite element variables.
!    We wish to determine the bandwidth necessary to store the stiffness
!    matrix associated with these variables.
!
!    An entry K(I,J) of the stiffness matrix must be zero unless the
!    variables I and J correspond to nodes N(I) and N(J) which are
!    common to some element.
!
!    In order to determine the bandwidth of the stiffness matrix, we
!    essentially seek a nonzero entry K(I,J) for which abs ( I - J )
!    is maximized.
!
!    The bandwidth M is defined in terms of the lower and upper bandwidths:
!
!      M = ML + 1 + MU
!
!    where
!
!      ML = maximum distance from any diagonal entry to a nonzero
!      entry in the same row, but earlier column,
!
!      MU = maximum distance from any diagonal entry to a nonzero
!      entry in the same row, but later column.
!
!    We assume the finite element variable adjacency relationship is 
!    symmetric, so we are guaranteed that ML = MU.
!
!    Note that the user is free to number the variables in any way
!    whatsoever, and to associate variables to nodes in any way,
!    so that some nodes have no variables, some have one, and some
!    have several.  
!
!    The storage of the indices of the variables is fairly simple.
!    In VAR, simply list all the variables associated with node 1, 
!    then all those associated with node 2, and so on.  Then set up
!    the pointer array VAR_NODE so that we can jump to the section of
!    VAR where the list begins for any particular node.
!
!    The routine does not check that each variable is only associated
!    with a single node.  This would normally be the case in a finite
!    element setting.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 September 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ELEMENT_ORDER, the order of the elements.
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of elements.
!
!    Input, integer ( kind = 4 ) ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM);
!    ELEMENT_NODE(I,J) is the global index of local node I in element J.
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, integer ( kind = 4 ) VAR_NODE(NODE_NUM+1), used to find the 
!    variables associated with a given node, which are in VAR in locations 
!    VAR_NODE(NODE) to VAR_NODE(NODE+1)-1.  Note that the last entry of
!    this array points to the location just after the last location in VAR.
!
!    Input, integer ( kind = 4 ) VAR_NUM, the number of variables.
!
!    Input, integer ( kind = 4 ) VAR(VAR_NUM), the indexes of the variables, 
!    which are presumably (but not necessarily) 1, 2, 3, ..., VAR_NUM.
!
!    Output, integer ( kind = 4 ) ML, MU, the lower and upper bandwidths of the 
!    matrix.
!
!    Output, integer ( kind = 4 ) M, the bandwidth of the matrix.
!
  implicit none

  integer ( kind = 4 ) element_num
  integer ( kind = 4 ) element_order
  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) var_num

  integer ( kind = 4 ) element
  integer ( kind = 4 ) element_node(element_order,element_num)
  integer ( kind = 4 ) m
  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
  integer ( kind = 4 ) node_global_i
  integer ( kind = 4 ) node_global_j
  integer ( kind = 4 ) node_local_i
  integer ( kind = 4 ) node_local_j
  integer ( kind = 4 ) var(var_num)
  integer ( kind = 4 ) var_global_i
  integer ( kind = 4 ) var_global_j
  integer ( kind = 4 ) var_local_i
  integer ( kind = 4 ) var_local_j
  integer ( kind = 4 ) var_node(node_num+1)

  ml = 0
  mu = 0

  do element = 1, element_num

    do node_local_i = 1, element_order
      node_global_i = element_node(node_local_i,element)

      do var_local_i = var_node(node_global_i), var_node(node_global_i+1)-1
        var_global_i = var(var_local_i)

        do node_local_j = 1, element_order
          node_global_j = element_node(node_local_j,element)

          do var_local_j = var_node(node_global_j), var_node(node_global_j+1)-1
            var_global_j = var(var_local_j)

            mu = max ( mu, var_global_j - var_global_i )
            ml = max ( ml, var_global_i - var_global_j )

          end do
        end do
      end do
    end do
  end do

  m = ml + 1 + mu

  return
end
subroutine basis_11_q4 ( q, i, p, phi, dphidx, dphidy )

!*****************************************************************************80
!
!! BASIS_11_Q4: one basis at one point for a Q4 element.
!
!  Discussion:
!
!    The routine is given the coordinates of the vertices of a quadrilateral.
!    It works directly with these coordinates, and does not refer to a 
!    reference element.
!
!    The sides of the element are presumed to lie along coordinate axes.
!
!    The routine evaluates the basis functions, and their X and Y derivatives.
!
!  Physical Element Q4:
!
!    |
!    |  4-----3
!    |  |     |
!    |  |     |
!    Y  |     |
!    |  |     |
!    |  |     |
!    |  1-----2
!    |
!    +-----X------>
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 March 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) Q(2,4), the coordinates of the vertices.
!    It is common to list these points in counter clockwise order.
!
!    Input, integer ( kind = 4 ) I, the index of the basis function.
!
!    Input, real ( kind = 8 ) P(2), the evaluation point.
!
!    Output, real ( kind = 8 ) PHI(4), the basis functions 
!    at the evaluation points.
!
!    Output, real ( kind = 8 ) DPHIDX(4), DPHIDY(4), the basis
!    derivatives at the evaluation points.
!
!  Local Parameter:
!
!    Local, real ( kind = 8 ) AREA, the area of the rectangle.
!
  implicit none

  real    ( kind = 8 ) area
  real    ( kind = 8 ) dphidx(4)
  real    ( kind = 8 ) dphidy(4)
  integer ( kind = 4 ) i
  real    ( kind = 8 ) p(2)
  real    ( kind = 8 ) phi
  real    ( kind = 8 ) q(2,4)

  area = ( q(1,3) - q(1,1) ) * ( q(2,3) - q(2,1) )

  if ( i == 1 ) then
    phi    =   ( q(1,3) - p(1) ) * ( q(2,3) - p(2) ) / area
    dphidx = -                     ( q(2,3) - p(2) ) / area
    dphidy = - ( q(1,3) - p(1) )                     / area
  else if ( i == 2 ) then
    phi    =   ( p(1) - q(1,1) ) * ( q(2,3) - p(2) ) / area
    dphidx =                       ( q(2,3) - p(2) ) / area
    dphidy = - ( p(1) - q(1,1) )                     / area
  else if ( i == 3 ) then
    phi    =   ( p(1) - q(1,1) ) * ( p(2) - q(2,1) ) / area
    dphidx =                       ( p(2) - q(2,1) ) / area
    dphidy =   ( p(1) - q(1,1) )                     / area
  else if ( i == 4 ) then
    phi    =   ( q(1,3) - p(1) ) * ( p(2) - q(2,1) ) / area
    dphidx = -                     ( p(2) - q(2,1) ) / area
    dphidy =   ( q(1,3) - p(1) )                     / area
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BASIS_11_Q4 - Fatal error!'
    write ( *, '(a)' ) '  Illegal basis function index.'
    stop
  end if

  return
end
subroutine basis_11_q4_test ( )

!*****************************************************************************80
!
!! BASIS_11_Q4_TEST verifies BASIS_11_Q4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 March 2009
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

  integer ( kind = 4 ), parameter :: node_num = 4

  real    ( kind = 8 ) dphidx(node_num,node_num)
  real    ( kind = 8 ) dphidy(node_num,node_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real    ( kind = 8 ) phi(node_num,node_num)
  real    ( kind = 8 ), dimension ( 2, node_num ) :: q = reshape ( (/ &
    2.0D+00, 0.0D+00, &
    4.0D+00, 4.0D+00, &
    0.0D+00, 3.0D+00, &
    0.0D+00, 1.0D+00 /), (/ 2, node_num /) )
  real    ( kind = 8 ) sum_x
  real    ( kind = 8 ) sum_y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BASIS_11_Q4_TEST:'
  write ( *, '(a)' ) '  Verify basis functions for element Q4.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of nodes = ', node_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Physical Nodes:'
  write ( *, '(a)' ) ' '
  do j = 1, node_num
    write ( *, '(2x,i8,2x,f7.3,2x,f7.3)' ) j, q(1:2,j)
  end do
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The basis function values at basis nodes'
  write ( *, '(a)' ) '  should form the identity matrix.'
  write ( *, '(a)' ) ' '

  do i = 1, node_num
    do j = 1, node_num
      call basis_11_q4 ( q, i, q(1:2,j), phi(i,j), dphidx(i,j), dphidy(i,j) )
    end do
  end do

  do i = 1, node_num
    write ( *, '(2x,10f7.3)' ) phi(i,1:node_num)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The X and Y derivatives should sum to 0.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      dPhidX sum      dPhidY sum'
  write ( *, '(a)' ) ' '
  do j = 1, node_num
    sum_x = sum ( dphidx(1:node_num,j) )
    sum_y = sum ( dphidy(1:node_num,j) )
    write ( *, '(2x,f14.8,2x,f14.8)' ) sum_x, sum_y
  end do

  return
end
subroutine basis_11_t3 ( t, i, p, qi, dqidx, dqidy )

!*****************************************************************************80
!
!! BASIS_11_T3: one basis at one point for the T3 element.
!
!  Discussion:
!
!    The routine is given the coordinates of the nodes of a triangle. 
!        
!           3
!          / \
!         /   \
!        /     \
!       1-------2
!
!    It evaluates the linear basis function Q(I)(X,Y) associated with
!    node I, which has the property that it is a linear function
!    which is 1 at node I and zero at the other two nodes.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T(2,3), the coordinates of the nodes.
!
!    Input, integer ( kind = 4 ) I, the index of the desired basis function.
!    I should be between 1 and 3.
!
!    Input, real ( kind = 8 ) P(2), the coordinates of a point at which the
!    basis function is to be evaluated.
!
!    Output, real ( kind = 8 ) QI, DQIDX, DQIDY, the values of the basis
!    function and its X and Y derivatives.
!
  implicit none

  real    ( kind = 8 ) area
  real    ( kind = 8 ) dqidx
  real    ( kind = 8 ) dqidy
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) ip1
  integer ( kind = 4 ) ip2
  real    ( kind = 8 ) p(2)
  real    ( kind = 8 ) qi
  real    ( kind = 8 ) t(2,3)

  area = abs ( t(1,1) * ( t(2,2) - t(2,3) ) &
             + t(1,2) * ( t(2,3) - t(2,1) ) &
             + t(1,3) * ( t(2,1) - t(2,2) ) )

  if ( area == 0.0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BASIS_11_T3 - Fatal error!'
    write ( *, '(a)' ) '  Element has zero area.'
    stop
  end if

  if ( i < 1 .or. 3 < i ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BASIS_11_T3 - Fatal error!'
    write ( *, '(a)' ) '  Basis index I is not between 1 and 3.'
    write ( *, '(a,i8)' ) '  I = ', i
    stop
  end if

  ip1 = i4_wrap ( i + 1, 1, 3 )
  ip2 = i4_wrap ( i + 2, 1, 3 )

  qi = ( ( t(1,ip2) - t(1,ip1) ) * ( p(2) - t(2,ip1) ) &
       - ( t(2,ip2) - t(2,ip1) ) * ( p(1) - t(1,ip1) ) ) / area

  dqidx = - ( t(2,ip2) - t(2,ip1) ) / area
  dqidy =   ( t(1,ip2) - t(1,ip1) ) / area

  return
end
subroutine basis_11_t3_test ( )

!*****************************************************************************80
!
!! BASIS_11_T3_TEST verifies BASIS_11_T3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 February 2006
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

  integer ( kind = 4 ), parameter :: node_num = 3

  real    ( kind = 8 ) dphidx(node_num,node_num)
  real    ( kind = 8 ) dphidy(node_num,node_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real    ( kind = 8 ) phi(node_num,node_num)
  real    ( kind = 8 ) sum_x
  real    ( kind = 8 ) sum_y
  real    ( kind = 8 ), dimension ( 2, node_num ) :: t = reshape ( (/ &
    2.0D+00, 0.0D+00, &
    4.0D+00, 3.0D+00, &
    0.0D+00, 4.0D+00 /), (/ 2, node_num /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BASIS_11_T3_TEST:'
  write ( *, '(a)' ) '  Verify basis functions for element T3.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of nodes = ', node_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Physical Nodes:'
  write ( *, '(a)' ) ' '
  do j = 1, node_num
    write ( *, '(2x,i8,2x,f7.3,2x,f7.3)' ) j, t(1:2,j)
  end do
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The basis function values at basis nodes'
  write ( *, '(a)' ) '  should form the identity matrix.'
  write ( *, '(a)' ) ' '

  do i = 1, node_num
    do j = 1, node_num
      call basis_11_t3 ( t, i, t(1:2,j), phi(i,j), dphidx(i,j), dphidy(i,j) )
    end do
  end do

  do i = 1, node_num
    write ( *, '(2x,10f7.3)' ) phi(i,1:node_num)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The X and Y derivatives should sum to 0.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      dPhidX sum    dPhidY sum'
  write ( *, '(a)' ) ' '

  do j = 1, node_num

    sum_x = sum ( dphidx(1:node_num,j) )
    sum_y = sum ( dphidy(1:node_num,j) )
    write ( *, '(2x,f14.8,f14.8)' ) sum_x, sum_y

  end do

  return
end
subroutine basis_11_t4 ( t, i, p, phi, dphidx, dphidy )

!*****************************************************************************80
!
!! BASIS_11_T4: one basis at one point for a T4 element.
!
!  Discussion:
!
!    The T4 element is the cubic bubble triangle.
!
!    The routine is given the coordinates of the vertices of a triangle.
!    It works directly with these coordinates, and does not refer to a 
!    reference element.
!
!    The sides of the triangle DO NOT have to lie along a coordinate
!    axis.
!
!    The routine evaluates the basis functions associated with each vertex,
!    and their derivatives with respect to X and Y.
!
!  Physical Element T4: 
!        
!            3
!           / \
!          /   \
!         /  4  \
!        /       \
!       1---------2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 March 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T(2,4), the coordinates of the vertices
!    of the triangle, and the coordinates of the centroid.  
!    It is common to list the first three points in counter clockwise
!    order.
!
!    Input, integer ( kind = 4 ) I, the index of the basis function.
!
!    Input, real ( kind = 8 ) P(2), the points where the basis function
!    is to be evaluated.
!
!    Output, real ( kind = 8 ) PHI, the value of the basis function
!    at the evaluation point.
!
!    Output, real ( kind = 8 ) DPHIDX, DPHIDY, the value of the 
!    derivatives at the evaluation point.
!
!  Local parameters:
!
!    Local, real ( kind = 8 ) AREA, is (twice) the area of the triangle.
!
  implicit none

  real    ( kind = 8 ) area
  real    ( kind = 8 ) dphidx
  real    ( kind = 8 ) dphidy
  real    ( kind = 8 ) dpsidx(4)
  real    ( kind = 8 ) dpsidy(4)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real    ( kind = 8 ) p(2)
  real    ( kind = 8 ) phi
  real    ( kind = 8 ) psi(4)
  real    ( kind = 8 ) t(2,4)

  area = t(1,1) * ( t(2,2) - t(2,3) ) &
       + t(1,2) * ( t(2,3) - t(2,1) ) &
       + t(1,3) * ( t(2,1) - t(2,2) )

  psi(1) =     (   ( t(1,3) - t(1,2) ) * ( p(2) - t(2,2) )     &
                 - ( t(2,3) - t(2,2) ) * ( p(1) - t(1,2) ) ) 
  dpsidx(1) =    - ( t(2,3) - t(2,2) )
  dpsidy(1) =      ( t(1,3) - t(1,2) )

  psi(2) =     (   ( t(1,1) - t(1,3) ) * ( p(2) - t(2,3) )     &
                 - ( t(2,1) - t(2,3) ) * ( p(1) - t(1,3) ) )
  dpsidx(2) =    - ( t(2,1) - t(2,3) )
  dpsidy(2) =      ( t(1,1) - t(1,3) )

  psi(3) =     (   ( t(1,2) - t(1,1) ) * ( p(2) - t(2,1) )     &
                 - ( t(2,2) - t(2,1) ) * ( p(1) - t(1,1) ) )
  dpsidx(3) =    - ( t(2,2) - t(2,1) )
  dpsidy(3) =      ( t(1,2) - t(1,1) )
!
!  Normalize the first three functions.
!
  psi(1:3)    =    psi(1:3) / area
  dpsidx(1:3) = dpsidx(1:3) / area
  dpsidy(1:3) = dpsidy(1:3) / area
!
!  Compute the cubic bubble function.
!
  psi(4) = 27.0D+00 * psi(1) * psi(2) * psi(3)

  dpsidx(4) = 27.0D+00 * ( &
                dpsidx(1) *    psi(2) *    psi(3) &
                +  psi(1) * dpsidx(2) *    psi(3) &
                +  psi(1) *    psi(2) * dpsidx(3) )

  dpsidy(4) = 27.0D+00 * ( &
                dpsidy(1) *    psi(2) *    psi(3) &
                +  psi(1) * dpsidy(2) *    psi(3) &
                +  psi(1) *    psi(2) * dpsidy(3) )
!
!  Subtract 1/3 of the cubic bubble function from each of the three linears.
!
  do j = 1, 3
    psi(j)    =    psi(j) -    psi(4) / 3.0D+00
    dpsidx(j) = dpsidx(j) - dpsidx(4) / 3.0D+00
    dpsidy(j) = dpsidy(j) - dpsidy(4) / 3.0D+00
  end do

  phi    = psi(i)
  dphidx = dpsidx(i)
  dphidy = dpsidy(i)

  return
end
subroutine basis_11_t4_test ( )

!*****************************************************************************80
!
!! BASIS_11_T4_TEST verifies BASIS_11_T4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 March 2009
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

  integer ( kind = 4 ), parameter :: node_num = 4

  real    ( kind = 8 ) dphidx(node_num,node_num)
  real    ( kind = 8 ) dphidy(node_num,node_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real    ( kind = 8 ) phi(node_num,node_num)
  real    ( kind = 8 ) sum_x
  real    ( kind = 8 ) sum_y
  real    ( kind = 8 ), dimension ( 2, node_num ) :: t = reshape ( (/ &
    2.0D+00, 0.0D+00, &
    4.0D+00, 3.0D+00, &
    0.0D+00, 4.0D+00, &
    0.0D+00, 0.0D+00 /), (/ 2, node_num /) )
!
!  The node associated with the fourth basis function is the centroid.
!
  t(1,4) = sum ( t(1,1:3) ) / 3.0D+00
  t(2,4) = sum ( t(2,1:3) ) / 3.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BASIS_11_T4_TEST:'
  write ( *, '(a)' ) '  Verify basis functions for element T4.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of nodes = ', node_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Physical Nodes:'
  write ( *, '(a)' ) ' '
  do j = 1, node_num
    write ( *, '(2x,i8,2x,f7.3,2x,f7.3)' ) j, t(1:2,j)
  end do
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The basis function values at basis nodes'
  write ( *, '(a)' ) '  should form the identity matrix.'
  write ( *, '(a)' ) ' '

  do i = 1, node_num
    do j = 1, node_num
      call basis_11_t4 ( t, i, t(1:2,j), phi(i,j), dphidx(i,j), dphidy(i,j) )
    end do
  end do

  do i = 1, node_num
    write ( *, '(2x,10f7.3)' ) phi(i,1:node_num)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The X and Y derivatives should sum to 0.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      dPhidX sum    dPhidY sum'
  write ( *, '(a)' ) ' '

  do j = 1, node_num

    sum_x = sum ( dphidx(1:node_num,j) )
    sum_y = sum ( dphidy(1:node_num,j) )
    write ( *, '(2x,f14.8,f14.8)' ) sum_x, sum_y

  end do

  return
end
subroutine basis_11_t6 ( t, i, p, bi, dbidx, dbidy )

!*****************************************************************************80
!
!! BASIS_11_T6: one basis at one point for the T6 element.
!
!  Discussion:
!
!    The routine is given the coordinates of the nodes of a triangle. 
!        
!           3
!          / \
!         6   5
!        /     \
!       1---4---2
!
!    It evaluates the quadratic basis function B(I)(X,Y) associated with
!    node I, which has the property that it is a quadratic function
!    which is 1 at node I and zero at the other five nodes.
!
!    This routine assumes that the sides of the triangle are straight,
!    so that the midside nodes fall on the line between two vertices.
!
!    This routine relies on the fact that each basis function can be
!    written as the product of two linear factors, which are easily
!    computed and normalized.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T(2,6), the coordinates of the nodes.
!
!    Input, integer ( kind = 4 ) I, the index of the desired basis function.
!    I should be between 1 and 6.
!
!    Input, real ( kind = 8 ) P(2), the coordinates of a point at which the
!    basis function is to be evaluated.
!
!    Output, real ( kind = 8 ) BI, DBIDX, DBIDY, the values of the basis
!    function and its X and Y derivatives.
!
  implicit none

  real    ( kind = 8 ) bi
  real    ( kind = 8 ) dbidx
  real    ( kind = 8 ) dbidy
  real    ( kind = 8 ) gf
  real    ( kind = 8 ) gn
  real    ( kind = 8 ) hf
  real    ( kind = 8 ) hn
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  real    ( kind = 8 ) p(2)
  real    ( kind = 8 ) t(2,6)

  if ( i < 1 .or. 6 < i ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BASIS_11_T6 - Fatal error!'
    write ( *, '(a)' ) '  Basis index I is not between 1 and 6.'
    write ( *, '(a,i8)' ) '  I = ', i
    stop
  end if
!
!  Determine the pairs of nodes.
!
  if ( i <= 3 ) then
    j1 = i4_wrap ( i + 1, 1, 3 )
    j2 = i4_wrap ( i + 2, 1, 3 )
    k1 = i + 3
    k2 = i4_wrap ( i + 5, 4, 6 )
  else
    j1 = i - 3
    j2 = i4_wrap ( i - 3 + 2, 1, 3 )
    k1 = i4_wrap ( i - 3 + 1, 1, 3 )
    k2 = i4_wrap ( i - 3 + 2, 1, 3 )
  end if
!
!  Evaluate the two linear factors GF and HF, 
!  and their normalizers GN and HN.
!
  gf = ( p(1)    - t(1,j1) ) * ( t(2,j2) - t(2,j1) ) &
     - ( t(1,j2) - t(1,j1) ) * ( p(2)    - t(2,j1) ) 

  gn = ( t(1,i)  - t(1,j1) ) * ( t(2,j2) - t(2,j1) ) &
     - ( t(1,j2) - t(1,j1) ) * ( t(2,i)  - t(2,j1) )   

  hf = ( p(1)    - t(1,k1) ) * ( t(2,k2) - t(2,k1) ) &
     - ( t(1,k2) - t(1,k1) ) * ( p(2)    - t(2,k1) ) 

  hn = ( t(1,i)  - t(1,k1) ) * ( t(2,k2) - t(2,k1) ) &
     - ( t(1,k2) - t(1,k1) ) * ( t(2,i)  - t(2,k1) ) 
!
!  Construct the basis function and its derivatives.
!
  bi =        ( gf                  / gn ) * (   hf                  / hn )

  dbidx =   ( ( t(2,j2) - t(2,j1) ) / gn ) * (   hf                  / hn ) &
          + (   gf                  / gn ) * ( ( t(2,k2) - t(2,k1) ) / hn )

  dbidy = - ( ( t(1,j2) - t(1,j1) ) / gn ) * (   hf                  / hn ) &
          - (   gf                  / gn ) * ( ( t(1,k2) - t(1,k1) ) / hn )

  return
end
subroutine basis_11_t6_test ( )

!*****************************************************************************80
!
!! BASIS_11_T6_TEST verifies BASIS_11_T6.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 February 2006
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

  integer ( kind = 4 ), parameter :: node_num = 6

  real    ( kind = 8 ) dphidx(node_num,node_num)
  real    ( kind = 8 ) dphidy(node_num,node_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real    ( kind = 8 ) phi(node_num,node_num)
  real    ( kind = 8 ) sum_x
  real    ( kind = 8 ) sum_y
  real    ( kind = 8 ), dimension ( 2, node_num ) :: t = reshape ( (/ &
    2.0D+00, 0.0D+00, &
    4.0D+00, 3.0D+00, &
    0.0D+00, 4.0D+00, &
    3.0D+00, 1.5D+00, &
    2.0D+00, 3.5D+00, &
    1.0D+00, 2.0D+00 /), (/ 2, node_num /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BASIS_11_T6_TEST:'
  write ( *, '(a)' ) '  Verify basis functions for element T6.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of nodes = ', node_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Physical Nodes:'
  write ( *, '(a)' ) ' '
  do j = 1, node_num
    write ( *, '(2x,i8,2x,f7.3,2x,f7.3)' ) j, t(1:2,j)
  end do
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The basis function values at basis nodes'
  write ( *, '(a)' ) '  should form the identity matrix.'
  write ( *, '(a)' ) ' '

  do i = 1, node_num
    do j = 1, node_num
      call basis_11_t6 ( t, i, t(1:2,j), phi(i,j), dphidx(i,j), dphidy(i,j) )
    end do
  end do

  do i = 1, node_num
    write ( *, '(2x,10f7.3)' ) phi(i,1:node_num)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The X and Y derivatives should sum to 0.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      dPhidX sum      dPhidY sum'
  write ( *, '(a)' ) ' '
  do j = 1, node_num
    sum_x = sum ( dphidx(1:node_num,j) )
    sum_y = sum ( dphidy(1:node_num,j) )
    write ( *, '(2x,f14.8,2x,f14.8)' ) sum_x, sum_y
  end do

  return
end
subroutine basis_mn_q4 ( q, n, p, phi, dphidx, dphidy )

!*****************************************************************************80
!
!! BASIS_MN_Q4: all bases at N points for a Q4 element.
!
!  Discussion:
!
!    The routine is given the coordinates of the vertices of a quadrilateral.
!    It works directly with these coordinates, and does not refer to a 
!    reference element.
!
!    The sides of the element are presumed to lie along coordinate axes.
!
!    The routine evaluates the basis functions, and their X and Y derivatives.
!
!  Physical Element Q4:
!
!    |
!    |  4-----3
!    |  |     |
!    |  |     |
!    Y  |     |
!    |  |     |
!    |  |     |
!    |  1-----2
!    |
!    +-----X------>
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) Q(2,4), the coordinates of the vertices.
!    It is common to list these points in counter clockwise order.
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) P(2,N), the evaluation points.
!
!    Output, real ( kind = 8 ) PHI(4,N), the basis functions 
!    at the evaluation points.
!
!    Output, real ( kind = 8 ) DPHIDX(4,N), DPHIDY(4,N), the basis
!    derivatives at the evaluation points.
!
!  Local Parameter:
!
!    Local, real ( kind = 8 ) AREA, the area of the rectangle.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) area
  real    ( kind = 8 ) dphidx(4,n)
  real    ( kind = 8 ) dphidy(4,n)
  real    ( kind = 8 ) p(2,n)
  real    ( kind = 8 ) phi(4,n)
  real    ( kind = 8 ) q(2,4)

  area =            ( q(1,3)             - q(1,1) ) &
                  * ( q(2,3)             - q(2,1) )

  phi(1,1:n) =      ( q(1,3) - p(1,1:n)           ) &
                  * ( q(2,3) - p(2,1:n)           )
  phi(2,1:n) =      (          p(1,1:n)  - q(1,1) ) &
                  * ( q(2,3) - p(2,1:n)           )
  phi(3,1:n) =      (          p(1,1:n)  - q(1,1) ) &
                  * (          p(2,1:n)  - q(2,1) )
  phi(4,1:n) =      ( q(1,3) - p(1,1:n)           ) &
                  * (          p(2,1:n)  - q(2,1) )
    
  dphidx(1,1:n) = - ( q(2,3) - p(2,1:n)           )
  dphidx(2,1:n) =   ( q(2,3) - p(2,1:n)           )
  dphidx(3,1:n) =   (          p(2,1:n)  - q(2,1) )
  dphidx(4,1:n) = - (          p(2,1:n)  - q(2,1) )
 
  dphidy(1,1:n) = - ( q(1,3) - p(1,1:n)           )
  dphidy(2,1:n) = - (          p(1,1:n)  - q(1,1) )
  dphidy(3,1:n) =   (          p(1,1:n)  - q(1,1) )
  dphidy(4,1:n) =   ( q(1,3) - p(1,1:n)           )
!
!  Normalize.
!
  phi(1:4,1:n)    = phi(1:4,1:n)    / area
  dphidx(1:4,1:n) = dphidx(1:4,1:n) / area
  dphidy(1:4,1:n) = dphidy(1:4,1:n) / area

  return
end
subroutine basis_mn_q4_test ( )

!*****************************************************************************80
!
!! BASIS_MN_Q4_TEST verifies BASIS_MN_Q4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 February 2006
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

  integer ( kind = 4 ), parameter :: node_num = 4

  real    ( kind = 8 ) dphidx(node_num,node_num)
  real    ( kind = 8 ) dphidy(node_num,node_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real    ( kind = 8 ) phi(node_num,node_num)
  real    ( kind = 8 ), dimension ( 2, node_num ) :: q = reshape ( (/ &
    3.0D+00, 1.0D+00, &
    5.0D+00, 1.0D+00, &
    5.0D+00, 4.0D+00, &
    3.0D+00, 4.0D+00 /), (/ 2, node_num /) )
  real    ( kind = 8 ) sum_x
  real    ( kind = 8 ) sum_y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  BASIS_MN_Q4_TEST: '
  write ( *, '(a)' ) '    Verify basis functions for element Q4.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of nodes = ', node_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Physical Nodes:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       I        X        Y'
  write ( *, '(a)' ) ' '
  do j = 1, node_num
    write ( *, '(2x,i8,2x,f7.3,2x,f7.3)' ) j, q(1:2,j)
  end do
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The basis function values at basis nodes'
  write ( *, '(a)' ) '  should form the identity matrix.'
  write ( *, '(a)' ) ' '

  call basis_mn_q4 ( q, node_num, q, phi, dphidx, dphidy )

  do i = 1, node_num
    write ( *, '(2x,10f7.3)' ) phi(i,1:node_num)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The X and Y derivatives should sum to 0.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      dPhidX sum      dPhidY sum'
  write ( *, '(a)' ) ' '

  do j = 1, node_num

    sum_x = sum ( dphidx(1:node_num,j) )
    sum_y = sum ( dphidy(1:node_num,j) )
    write ( *, '(2x,f14.8,2x,f14.8)' ) sum_x, sum_y

  end do

  return
end
subroutine basis_mn_t3 ( t, n, p, phi, dphidx, dphidy )

!*****************************************************************************80
!
!! BASIS_MN_T3: all bases at N points for a T3 element.
!
!  Discussion:
!
!    The routine is given the coordinates of the vertices of a triangle.
!    It works directly with these coordinates, and does not refer to a 
!    reference element.
!
!    The sides of the triangle DO NOT have to lie along a coordinate
!    axis.
!
!    The routine evaluates the basis functions associated with each vertex,
!    and their derivatives with respect to X and Y.
!
!  Physical Element T3: 
!        
!            3
!           / \
!          /   \
!         /     \
!        /       \
!       1---------2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T(2,3), the coordinates of the vertices
!    of the triangle.  It is common to list these points in counter clockwise
!    order.
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) P(2,N), the points where the basis functions 
!    are to be evaluated.
!
!    Output, real ( kind = 8 ) PHI(3,N), the value of the basis functions 
!    at the evaluation points.
!
!    Output, real ( kind = 8 ) DPHIDX(3,N), DPHIDY(3,N), the value of the 
!    derivatives at the evaluation points.
!
!  Local parameters:
!
!    Local, real ( kind = 8 ) AREA, is (twice) the area of the triangle.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) area
  real    ( kind = 8 ) dphidx(3,n)
  real    ( kind = 8 ) dphidy(3,n)
  real    ( kind = 8 ) p(2,n)
  real    ( kind = 8 ) phi(3,n)
  real    ( kind = 8 ) t(2,3)

  area = t(1,1) * ( t(2,2) - t(2,3) ) &
       + t(1,2) * ( t(2,3) - t(2,1) ) &
       + t(1,3) * ( t(2,1) - t(2,2) )

  if ( area == 0.0D+00 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BASIS_MN_T3 - Fatal error!'
    write ( *, '(a)' ) '  Element has zero area.'
    stop

  end if

  phi(1,1:n) =     (   ( t(1,3) - t(1,2) ) * ( p(2,1:n) - t(2,2) )     &
                     - ( t(2,3) - t(2,2) ) * ( p(1,1:n) - t(1,2) ) ) 
  dphidx(1,1:n) =    - ( t(2,3) - t(2,2) )
  dphidy(1,1:n) =      ( t(1,3) - t(1,2) )

  phi(2,1:n) =     (   ( t(1,1) - t(1,3) ) * ( p(2,1:n) - t(2,3) )     &
                     - ( t(2,1) - t(2,3) ) * ( p(1,1:n) - t(1,3) ) )
  dphidx(2,1:n) =    - ( t(2,1) - t(2,3) )
  dphidy(2,1:n) =      ( t(1,1) - t(1,3) )

  phi(3,1:n) =     (   ( t(1,2) - t(1,1) ) * ( p(2,1:n) - t(2,1) )     &
                     - ( t(2,2) - t(2,1) ) * ( p(1,1:n) - t(1,1) ) )
  dphidx(3,1:n) =    - ( t(2,2) - t(2,1) )
  dphidy(3,1:n) =      ( t(1,2) - t(1,1) )
!
!  Normalize.
!
  phi(1:3,1:n) = phi(1:3,1:n) / area
  dphidx(1:3,1:n) = dphidx(1:3,1:n) / area
  dphidy(1:3,1:n) = dphidy(1:3,1:n) / area

  return
end
subroutine basis_mn_t3_test ( )

!*****************************************************************************80
!
!! BASIS_MN_T3_TEST verifies BASIS_MN_T3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 February 2006
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

  integer ( kind = 4 ), parameter :: node_num = 3

  real    ( kind = 8 ) dphidx(node_num,node_num)
  real    ( kind = 8 ) dphidy(node_num,node_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real    ( kind = 8 ) phi(node_num,node_num)
  real    ( kind = 8 ) sum_x
  real    ( kind = 8 ) sum_y
  real    ( kind = 8 ), dimension ( 2, node_num ) :: t = reshape ( (/ &
    2.0D+00, 0.0D+00, &
    4.0D+00, 3.0D+00, &
    0.0D+00, 4.0D+00 /), (/ 2, node_num /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BASIS_MN_T3_TEST:'
  write ( *, '(a)' ) '  Verify basis functions for element T3.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of nodes = ', node_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Physical Nodes:'
  write ( *, '(a)' ) ' '
  do j = 1, node_num
    write ( *, '(2x,i8,2x,f7.3,2x,f7.3)' ) j, t(1:2,j)
  end do
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The basis function values at basis nodes'
  write ( *, '(a)' ) '  should form the identity matrix.'
  write ( *, '(a)' ) ' '
  call basis_mn_t3 ( t, node_num, t, phi, dphidx, dphidy )

  do i = 1, node_num
    write ( *, '(2x,10f7.3)' ) phi(i,1:node_num)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The X and Y derivatives should sum to 0.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      dPhidX sum    dPhidY sum'
  write ( *, '(a)' ) ' '

  do j = 1, node_num

    sum_x = sum ( dphidx(1:node_num,j) )
    sum_y = sum ( dphidy(1:node_num,j) )
    write ( *, '(2x,f14.8,2x,f14.8)' ) sum_x, sum_y

  end do

  return
end
subroutine basis_mn_t4 ( t, n, p, phi, dphidx, dphidy )

!*****************************************************************************80
!
!! BASIS_MN_T4: all bases at N points for a T4 element.
!
!  Discussion:
!
!    The T4 element is the cubic bubble triangle.
!
!    The routine is given the coordinates of the vertices of a triangle.
!    It works directly with these coordinates, and does not refer to a 
!    reference element.
!
!    The sides of the triangle DO NOT have to lie along a coordinate
!    axis.
!
!    The routine evaluates the basis functions associated with each vertex,
!    and their derivatives with respect to X and Y.
!
!  Physical Element T4: 
!        
!            3
!           / \
!          /   \
!         /  4  \
!        /       \
!       1---------2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T(2,4), the coordinates of the vertices
!    of the triangle, and the coordinates of the centroid.  
!    It is common to list the first three points in counter clockwise
!    order.
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) P(2,N), the points where the basis functions 
!    are to be evaluated.
!
!    Output, real ( kind = 8 ) PHI(4,N), the value of the basis functions 
!    at the evaluation points.
!
!    Output, real ( kind = 8 ) DPHIDX(4,N), DPHIDY(4,N), the value of the 
!    derivatives at the evaluation points.
!
!  Local parameters:
!
!    Local, real ( kind = 8 ) AREA, is (twice) the area of the triangle.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) area
  real    ( kind = 8 ) dphidx(4,n)
  real    ( kind = 8 ) dphidy(4,n)
  integer ( kind = 4 ) i
  real    ( kind = 8 ) p(2,n)
  real    ( kind = 8 ) phi(4,n)
  real    ( kind = 8 ) t(2,4)

  area = t(1,1) * ( t(2,2) - t(2,3) ) &
       + t(1,2) * ( t(2,3) - t(2,1) ) &
       + t(1,3) * ( t(2,1) - t(2,2) )

  phi(1,1:n) =     (   ( t(1,3) - t(1,2) ) * ( p(2,1:n) - t(2,2) )     &
                     - ( t(2,3) - t(2,2) ) * ( p(1,1:n) - t(1,2) ) ) 
  dphidx(1,1:n) =    - ( t(2,3) - t(2,2) )
  dphidy(1,1:n) =      ( t(1,3) - t(1,2) )

  phi(2,1:n) =     (   ( t(1,1) - t(1,3) ) * ( p(2,1:n) - t(2,3) )     &
                     - ( t(2,1) - t(2,3) ) * ( p(1,1:n) - t(1,3) ) )
  dphidx(2,1:n) =    - ( t(2,1) - t(2,3) )
  dphidy(2,1:n) =      ( t(1,1) - t(1,3) )

  phi(3,1:n) =     (   ( t(1,2) - t(1,1) ) * ( p(2,1:n) - t(2,1) )     &
                     - ( t(2,2) - t(2,1) ) * ( p(1,1:n) - t(1,1) ) )
  dphidx(3,1:n) =    - ( t(2,2) - t(2,1) )
  dphidy(3,1:n) =      ( t(1,2) - t(1,1) )
!
!  Normalize the first three functions.
!
  phi(1:3,1:n)    =    phi(1:3,1:n) / area
  dphidx(1:3,1:n) = dphidx(1:3,1:n) / area
  dphidy(1:3,1:n) = dphidy(1:3,1:n) / area
!
!  Compute the cubic bubble function.
!
  phi(4,1:n) = 27.0D+00 * phi(1,1:n) * phi(2,1:n) * phi(3,1:n)

  dphidx(4,1:n) = 27.0D+00 * ( &
                dphidx(1,1:n) *    phi(2,1:n) *    phi(3,1:n) &
                +  phi(1,1:n) * dphidx(2,1:n) *    phi(3,1:n) &
                +  phi(1,1:n) *    phi(2,1:n) * dphidx(3,1:n) )

  dphidy(4,1:n) = 27.0D+00 * ( &
                dphidy(1,1:n) *    phi(2,1:n) *    phi(3,1:n) &
                +  phi(1,1:n) * dphidy(2,1:n) *    phi(3,1:n) &
                +  phi(1,1:n) *    phi(2,1:n) * dphidy(3,1:n) )
!
!  Subtract 1/3 of the cubic bubble function from each of the three linears.
!
  do i = 1, 3
    phi(i,1:n)    =    phi(i,1:n) -    phi(4,1:n) / 3.0D+00
    dphidx(i,1:n) = dphidx(i,1:n) - dphidx(4,1:n) / 3.0D+00
    dphidy(i,1:n) = dphidy(i,1:n) - dphidy(4,1:n) / 3.0D+00
  end do

  return
end
subroutine basis_mn_t4_test ( )

!*****************************************************************************80
!
!! BASIS_MN_T4_TEST verifies BASIS_MN_T4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 February 2006
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

  integer ( kind = 4 ), parameter :: node_num = 4

  real    ( kind = 8 ) dphidx(node_num,node_num)
  real    ( kind = 8 ) dphidy(node_num,node_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real    ( kind = 8 ) phi(node_num,node_num)
  real    ( kind = 8 ) sum_x
  real    ( kind = 8 ) sum_y
  real    ( kind = 8 ), dimension ( 2, node_num ) :: t = reshape ( (/ &
    2.0D+00, 0.0D+00, &
    4.0D+00, 2.0D+00, &
    0.0D+00, 4.0D+00, &
    2.0D+00, 2.0D+00 /), (/ 2, node_num /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BASIS_MN_T4_TEST:'
  write ( *, '(a)' ) '  Verify basis functions for element T4.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of nodes = ', node_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Physical Nodes:'
  write ( *, '(a)' ) ' '
  do j = 1, node_num
    write ( *, '(2x,i8,2x,f7.3,2x,f7.3)' ) j, t(1:2,j)
  end do
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The basis function values at basis nodes'
  write ( *, '(a)' ) '  should form the identity matrix.'
  write ( *, '(a)' ) ' '
  call basis_mn_t4 ( t, node_num, t, phi, dphidx, dphidy )

  do i = 1, node_num
    write ( *, '(2x,10f7.3)' ) phi(i,1:node_num)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The X and Y derivatives should sum to 0.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      dPhidX sum    dPhidY sum'
  write ( *, '(a)' ) ' '

  do j = 1, node_num

    sum_x = sum ( dphidx(1:node_num,j) )
    sum_y = sum ( dphidy(1:node_num,j) )
    write ( *, '(2x,f14.8,2x,f14.8)' ) sum_x, sum_y

  end do

  return
end
subroutine basis_mn_t6 ( t, n, p, phi, dphidx, dphidy )

!*****************************************************************************80
!
!! BASIS_MN_T6: all bases at N points for a T6 element.
!
!  Discussion:
!
!    The routine is given the coordinates of the vertices and midside
!    nodes of a triangle.  It works directly with these coordinates, and does 
!    not refer to a reference element.
!
!    This routine requires that the midside nodes be "in line"
!    with the vertices, that is, that the sides of the triangle be
!    straight.  However, the midside nodes do not actually have to
!    be halfway along the side of the triangle.  
!
!  Physical element T6:
!
!    This picture indicates the assumed ordering of the six nodes
!    of the triangle.
!
!             3
!            / \
!           /   \
!          6     5
!         /       \
!        /         \
!       1-----4-----2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T(2,6), the nodal oordinates of the element.
!    It is common to list these points in counter clockwise order.
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) P(2,N), the coordinates of the point where
!    the basis functions are to be evaluated.
!
!    Output, real ( kind = 8 ) PHI(6,N), the basis functions at the
!    evaluation points.
!
!    Output, real ( kind = 8 ) DPHIDX(6,N), DPHIDY(6,N), the derivatives 
!    of the basis functions at the evaluation points.
!
!  Local Parameters:
!
!    Local, real ( kind = 8 ) AREA, is (twice) the area of the triangle.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) dphidx(6,n)
  real    ( kind = 8 ) dphidy(6,n)
  real    ( kind = 8 ) gn(n)
  real    ( kind = 8 ) gx(n)
  real    ( kind = 8 ) hn(n)
  real    ( kind = 8 ) hx(n)
  real    ( kind = 8 ) p(2,n)
  real    ( kind = 8 ) phi(6,n)
  real    ( kind = 8 ) t(2,6)
!
!  Basis function 1: PHI(X,Y) = G(3,2) * H(6,4) / normalization.
!
  gx(1:n) = ( p(1,1:n) - t(1,2) ) * ( t(2,3)   - t(2,2) ) &
          - ( t(1,3)   - t(1,2) ) * ( p(2,1:n) - t(2,2) )

  gn(1:n) = ( t(1,1)   - t(1,2) ) * ( t(2,3)   - t(2,2) ) &
          - ( t(1,3)   - t(1,2) ) * ( t(2,1)   - t(2,2) )

  hx(1:n) = ( p(1,1:n) - t(1,4) ) * ( t(2,6)   - t(2,4) ) &
          - ( t(1,6)   - t(1,4) ) * ( p(2,1:n) - t(2,4) )

  hn(1:n) = ( t(1,1)   - t(1,4) ) * ( t(2,6)   - t(2,4) ) &
          - ( t(1,6)   - t(1,4) ) * ( t(2,1)   - t(2,4) )

  phi(1,1:n) =     ( gx(1:n) * hx(1:n) ) / ( gn(1:n) * hn(1:n) )
  dphidx(1,1:n) =  (      ( t(2,3) - t(2,2) ) * hx(1:n) &
                   + gx(1:n) * ( t(2,6) - t(2,4) ) ) / ( gn(1:n) * hn(1:n) )
  dphidy(1,1:n) = -(      ( t(1,3) - t(1,2) ) * hx(1:n) &
                   + gx(1:n) * ( t(1,6) - t(1,4) ) ) / ( gn(1:n) * hn(1:n) )
!
!  Basis function 2: PHI(X,Y) = G(3,1) * H(4,5) / normalization.
!
  gx(1:n) = ( p(1,1:n) - t(1,1) ) * ( t(2,3)   - t(2,1) ) &
          - ( t(1,3)   - t(1,1) ) * ( p(2,1:n) - t(2,1) )

  gn(1:n) = ( t(1,2)   - t(1,1) ) * ( t(2,3)   - t(2,1) ) &
          - ( t(1,3)   - t(1,1) ) * ( t(2,2)   - t(2,1) )

  hx(1:n) = ( p(1,1:n) - t(1,5) ) * ( t(2,4)   - t(2,5) ) &
          - ( t(1,4)   - t(1,5) ) * ( p(2,1:n) - t(2,5) )

  hn(1:n) = ( t(1,2)   - t(1,5) ) * ( t(2,4)   - t(2,5) ) &
          - ( t(1,4)   - t(1,5) ) * ( t(2,2)   - t(2,5) )

  phi(2,1:n) = ( gx(1:n) * hx(1:n) ) / ( gn(1:n) * hn(1:n) )
  dphidx(2,1:n) =  (      ( t(2,3) - t(2,1) ) * hx(1:n) &
               + gx(1:n) * ( t(2,4) - t(2,5) ) ) / ( gn(1:n) * hn(1:n) )
  dphidy(2,1:n) = -(      ( t(1,3) - t(1,1) ) * hx(1:n) &
               + gx(1:n) * ( t(1,4) - t(1,5) ) ) / ( gn(1:n) * hn(1:n) )
!
!  Basis function 3: PHI(X,Y) = G(1,2) * H(5,6) / normalization.
!
  gx(1:n) = ( p(1,1:n) - t(1,2) ) * ( t(2,1)   - t(2,2) ) &
          - ( t(1,1)   - t(1,2) ) * ( p(2,1:n) - t(2,2) )

  gn(1:n) = ( t(1,3)   - t(1,2) ) * ( t(2,1)   - t(2,2) ) &
          - ( t(1,1)   - t(1,2) ) * ( t(2,3)   - t(2,2) )

  hx(1:n) = ( p(1,1:n) - t(1,6) ) * ( t(2,5)   - t(2,6) ) &
          - ( t(1,5)   - t(1,6) ) * ( p(2,1:n) - t(2,6) )

  hn(1:n) = ( t(1,3)   - t(1,6) ) * ( t(2,5)   - t(2,6) ) &
          - ( t(1,5)   - t(1,6) ) * ( t(2,3)   - t(2,6) )

  phi(3,1:n) = ( gx(1:n) * hx(1:n) ) / ( gn(1:n) * hn(1:n) )
  dphidx(3,1:n) =  (      ( t(2,1) - t(2,2) ) * hx(1:n) &
               + gx(1:n) * ( t(2,5) - t(2,6) ) ) / ( gn(1:n) * hn(1:n) )
  dphidy(3,1:n) = -(      ( t(1,1) - t(1,2) ) * hx(1:n) &
               + gx(1:n) * ( t(1,5) - t(1,6) ) ) / ( gn(1:n) * hn(1:n) )
!
!  Basis function 4: PHI(X,Y) = G(1,3) * H(2,3) / normalization.
!
  gx(1:n) = ( p(1,1:n) - t(1,3) ) * ( t(2,1)   - t(2,3) ) &
          - ( t(1,1)   - t(1,3) ) * ( p(2,1:n) - t(2,3) )

  gn(1:n) = ( t(1,4)   - t(1,3) ) * ( t(2,1)   - t(2,3) ) &
          - ( t(1,1)   - t(1,3) ) * ( t(2,4)   - t(2,3) )

  hx(1:n) = ( p(1,1:n) - t(1,3) ) * ( t(2,2)   - t(2,3) ) &
          - ( t(1,2)   - t(1,3) ) * ( p(2,1:n) - t(2,3) )

  hn(1:n) = ( t(1,4)   - t(1,3) ) * ( t(2,2)   - t(2,3) ) &
          - ( t(1,2)   - t(1,3) ) * ( t(2,4)   - t(2,3) )

  phi(4,1:n) = ( gx(1:n) * hx(1:n) ) / ( gn(1:n) * hn(1:n) )
  dphidx(4,1:n) =  (      ( t(2,1) - t(2,3) ) * hx(1:n) &
               + gx(1:n) * ( t(2,2) - t(2,3) ) ) / ( gn(1:n) * hn(1:n) )
  dphidy(4,1:n) = -(      ( t(1,1) - t(1,3) ) * hx(1:n) &
               + gx(1:n) * ( t(1,2) - t(1,3) ) ) / ( gn(1:n) * hn(1:n) )
!
!  Basis function 5: PHI(X,Y) = G(2,1) * H(3,1) / normalization.
!
  gx(1:n) = ( p(1,1:n) - t(1,1) ) * ( t(2,2)   - t(2,1) ) &
          - ( t(1,2)   - t(1,1) ) * ( p(2,1:n) - t(2,1) )

  gn(1:n) = ( t(1,5)   - t(1,1) ) * ( t(2,2)   - t(2,1) ) &
          - ( t(1,2)   - t(1,1) ) * ( t(2,5)   - t(2,1) )

  hx(1:n) = ( p(1,1:n) - t(1,1) ) * ( t(2,3)   - t(2,1) ) &
          - ( t(1,3)   - t(1,1) ) * ( p(2,1:n) - t(2,1) )

  hn(1:n) = ( t(1,5)   - t(1,1) ) * ( t(2,3)   - t(2,1) ) &
          - ( t(1,3)   - t(1,1) ) * ( t(2,5)   - t(2,1) )

  phi(5,1:n) = ( gx(1:n) * hx(1:n) ) / ( gn(1:n) * hn(1:n) )
  dphidx(5,1:n) =  (      ( t(2,2) - t(2,1) ) * hx(1:n) &
               + gx(1:n) * ( t(2,3) - t(2,1) ) ) / ( gn(1:n) * hn(1:n) )
  dphidy(5,1:n) = -(      ( t(1,2) - t(1,1) ) * hx(1:n) &
               + gx(1:n) * ( t(1,3) - t(1,1) ) ) / ( gn(1:n) * hn(1:n) )
!
!  Basis function 6: PHI(X,Y) = G(1,2) * H(3,2) / normalization.
!
  gx(1:n) = ( p(1,1:n) - t(1,2) ) * ( t(2,1)   - t(2,2) ) &
          - ( t(1,1)   - t(1,2) ) * ( p(2,1:n) - t(2,2) )

  gn(1:n) = ( t(1,6)   - t(1,2) ) * ( t(2,1)   - t(2,2) ) &
          - ( t(1,1)   - t(1,2) ) * ( t(2,6)   - t(2,2) )

  hx(1:n) = ( p(1,1:n) - t(1,2) ) * ( t(2,3)   - t(2,2) ) &
          - ( t(1,3)   - t(1,2) ) * ( p(2,1:n) - t(2,2) )

  hn(1:n) = ( t(1,6)   - t(1,2) ) * ( t(2,3)   - t(2,2) ) &
          - ( t(1,3)   - t(1,2) ) * ( t(2,6)   - t(2,2) )

  phi(6,1:n) = ( gx(1:n) * hx(1:n) ) / ( gn(1:n) * hn(1:n) )
  dphidx(6,1:n) =  (      ( t(2,1) - t(2,2) ) * hx(1:n) &
               + gx(1:n) * ( t(2,3) - t(2,2) ) ) / ( gn(1:n) * hn(1:n) )
  dphidy(6,1:n) = -(      ( t(1,1) - t(1,2) ) * hx(1:n) &
               + gx(1:n) * ( t(1,3) - t(1,2) ) ) / ( gn(1:n) * hn(1:n) )

  return
end
subroutine basis_mn_t6_test ( )

!*****************************************************************************80
!
!! BASIS_MN_T6_TEST verifies BASIS_MN_T6.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 February 2006
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

  integer ( kind = 4 ), parameter :: node_num = 6

  real    ( kind = 8 ) dphidx(node_num,node_num)
  real    ( kind = 8 ) dphidy(node_num,node_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real    ( kind = 8 ) phi(node_num,node_num)
  real    ( kind = 8 ) sum_x
  real    ( kind = 8 ) sum_y
  real    ( kind = 8 ), dimension ( 2, node_num ) :: t = reshape ( (/ &
    2.0D+00, 0.0D+00, &
    4.0D+00, 3.0D+00, &
    0.0D+00, 4.0D+00, &
    3.0D+00, 1.5D+00, &
    2.0D+00, 3.5D+00, &
    1.0D+00, 2.0D+00 /), (/ 2, node_num /) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BASIS_MN_T6_TEST:'
  write ( *, '(a)' ) '  Verify basis functions for element T6.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of nodes = ', node_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Physical Nodes:'
  write ( *, '(a)' ) ' '
  do j = 1, node_num
    write ( *, '(2x,i8,2x,f7.3,2x,f7.3)' ) j, t(1:2,j)
  end do
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The basis function values at basis nodes'
  write ( *, '(a)' ) '  should form the identity matrix.'
  write ( *, '(a)' ) ' '

  call basis_mn_t6 ( t, node_num, t, phi, dphidx, dphidy )

  do i = 1, node_num
    write ( *, '(2x,10f7.3)' ) phi(i,1:node_num)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The X and Y derivatives should sum to 0.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      dPhidX sum    dPhidY sum'
  write ( *, '(a)' ) ' '
  do j = 1, node_num
    sum_x = sum ( dphidx(1:node_num,j) )
    sum_y = sum ( dphidy(1:node_num,j) )
    write ( *, '(2x,2f14.8)' ) sum_x, sum_y
  end do

  return
end
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

  real    ( kind = 8 ) angle
  real    ( kind = 8 ) degrees_to_radians
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00

  degrees_to_radians = ( angle / 180.0D+00 ) * pi

  return
end
subroutine derivative_average_t3 ( node_num, node_xy, element_num, &
  element_node, c, dcdx, dcdy )

!*****************************************************************************80
!
!! DERIVATIVE_AVERAGE_T3 averages derivatives at the nodes of a T3 mesh.
!
!  Discussion:
!
!    This routine can be used to compute an averaged nodal value of any
!    quantity associated with the finite element function.  At a node 
!    that is shared by several elements, the fundamental function
!    U will be continuous, but its spatial derivatives, for instance,
!    will generally be discontinuous.  This routine computes the
!    value of the spatial derivatives in each element, and averages
!    them, to make a reasonable assignment of a nodal value.
!
!    In this version of the routine, the average is not weighted.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, real ( kind = 8 ) NODE_XY(2,NODE_NUM), the coordinates of 
!    the nodes.
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of elements.
!
!    Input, integer ( kind = 4 ) ELEMENT_NODE(3,ELEMENT_NUM), 
!    the element->node data.
!
!    Input, real ( kind = 8 ) C(NODE_NUM), the finite element coefficient
!    vector.
!
!    Output, real ( kind = 8 ) DCDX(NODE_NUM), DCDY(NODE_NUM), the averaged
!    values of dCdX and dCdY at the nodes.
!
  implicit none

  integer ( kind = 4 ) element_num
  integer ( kind = 4 ), parameter :: element_order = 3
  integer ( kind = 4 ) node_num

  real    ( kind = 8 ) c(node_num)
  real    ( kind = 8 ) dcdx(node_num)
  real    ( kind = 8 ) dcdy(node_num)
  real    ( kind = 8 ) dphidx(element_order,element_order)
  real    ( kind = 8 ) dphidy(element_order,element_order)
  integer ( kind = 4 ) element
  integer ( kind = 4 ) element_node(element_order,element_num)
  integer ( kind = 4 ) node_count(node_num)
  integer ( kind = 4 ) node_global1
  integer ( kind = 4 ) node_global2
  integer ( kind = 4 ) node_local1
  integer ( kind = 4 ) node_local2
  real    ( kind = 8 ) node_xy(2,node_num)
  real    ( kind = 8 ) phi(element_order,element_order)
  real    ( kind = 8 ) t(2,element_order)

  node_count(1:node_num) = 0
  dcdx(1:node_num) = 0.0D+00
  dcdy(1:node_num) = 0.0D+00
!
!  Consider every element.
!
  do element = 1, element_num
!
!  Get the coordinates of the nodes of the element.
!
    t(1:2,1:element_order) = node_xy(1:2,element_node(1:element_order,element))
!
!  Evaluate the X and Y derivatives of the basis functions at the nodes.
!
    call basis_mn_t3 ( t, element_order, t, phi, dphidx, dphidy )
!
!  Evaluate dCdX and dCdY at each node in the element, and add
!  them to the running totals.
!
    do node_local1 = 1, element_order

      node_global1 = element_node(node_local1,element)

      do node_local2 = 1, element_order

        node_global2 = element_node(node_local2,element)

        dcdx(node_global1) = dcdx(node_global1) &
          + c(node_global2) * dphidx(node_local2,node_local1)

        dcdy(node_global1) = dcdy(node_global1) &
          + c(node_global2) * dphidy(node_local2,node_local1)

      end do

      node_count(node_global1) = node_count(node_global1) + 1

    end do

  end do
!
!  Average the running totals.
!
  dcdx(1:node_num) = dcdx(1:node_num) &
    / real ( node_count(1:node_num), kind = 8 )

  dcdy(1:node_num) = dcdy(1:node_num) &
    / real ( node_count(1:node_num), kind = 8 )

  return
end
subroutine derivative_average_t6 ( node_num, node_xy, element_num, &
  element_node, c, dcdx, dcdy )

!*****************************************************************************80
!
!! DERIVATIVE_AVERAGE_T6 averages derivatives at the nodes of a T6 mesh.
!
!  Discussion:
!
!    This routine can be used to compute an averaged nodal value of any
!    quantity associated with the finite element function.  At a node 
!    that is shared by several elements, the fundamental function
!    U will be continuous, but its spatial derivatives, for instance,
!    will generally be discontinuous.  This routine computes the
!    value of the spatial derivatives in each element, and averages
!    them, to make a reasonable assignment of a nodal value.
!
!    In this version of the routine, the average is not weighted.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 April 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, real ( kind = 8 ) NODE_XY(2,NODE_NUM), the coordinates of 
!    the nodes.
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of elements.
!
!    Input, integer ( kind = 4 ) ELEMENT_NODE(6,ELEMENT_NUM), 
!    the element->node data.
!
!    Input, real ( kind = 8 ) C(NODE_NUM), the finite element coefficient
!    vector.
!
!    Output, real ( kind = 8 ) DCDX(NODE_NUM), DCDY(NODE_NUM), the averaged
!    values of dCdX and dCdY at the nodes.
!
  implicit none

  integer ( kind = 4 ) element_num
  integer ( kind = 4 ), parameter :: element_order = 6
  integer ( kind = 4 ) node_num

  real    ( kind = 8 ) c(node_num)
  real    ( kind = 8 ) dcdx(node_num)
  real    ( kind = 8 ) dcdy(node_num)
  real    ( kind = 8 ) dphidx(element_order,element_order)
  real    ( kind = 8 ) dphidy(element_order,element_order)
  integer ( kind = 4 ) element
  integer ( kind = 4 ) element_node(element_order,element_num)
  integer ( kind = 4 ) node_count(node_num)
  integer ( kind = 4 ) node_global1
  integer ( kind = 4 ) node_global2
  integer ( kind = 4 ) node_local1
  integer ( kind = 4 ) node_local2
  real    ( kind = 8 ) node_xy(2,node_num)
  real    ( kind = 8 ) phi(element_order,element_order)
  real    ( kind = 8 ) t(2,element_order)

  node_count(1:node_num) = 0
  dcdx(1:node_num) = 0.0D+00
  dcdy(1:node_num) = 0.0D+00
!
!  Consider every element.
!
  do element = 1, element_num
!
!  Get the coordinates of the nodes of the element.
!
    t(1:2,1:element_order) = node_xy(1:2,element_node(1:element_order,element))
!
!  Evaluate the X and Y derivatives of the basis functions at the nodes.
!
    call basis_mn_t6 ( t, element_order, t, phi, dphidx, dphidy )
!
!  Evaluate dCdX and dCdY at each node in the element, and add
!  them to the running totals.
!
    do node_local1 = 1, element_order

      node_global1 = element_node(node_local1,element)

      do node_local2 = 1, element_order

        node_global2 = element_node(node_local2,element)

        dcdx(node_global1) = dcdx(node_global1) &
          + c(node_global2) * dphidx(node_local2,node_local1)

        dcdy(node_global1) = dcdy(node_global1) &
          + c(node_global2) * dphidy(node_local2,node_local1)

      end do

      node_count(node_global1) = node_count(node_global1) + 1

    end do

  end do
!
!  Average the running totals.
!
  dcdx(1:node_num) = dcdx(1:node_num) &
    / real ( node_count(1:node_num), kind = 8 )

  dcdy(1:node_num) = dcdy(1:node_num) &
    / real ( node_count(1:node_num), kind = 8 )

  return
end
subroutine div_q4 ( m, n, u, v, q, div, vort )

!*****************************************************************************80
!
!! DIV_Q4 estimates the divergence and vorticity of a discrete field.
!
!  Discussion:
!
!    The routine is given the values of a vector field ( U(X,Y), V(X,Y) ) at
!    an array of points ( X(1:M), Y(1:N) ).
!
!    The routine models the vector field over the interior of this region using
!    a bilinear interpolant.  It then uses the interpolant to estimate the
!    value of the divergence:
!
!      DIV(X,Y) = dU/dX + dV/dY
!
!    and the vorticity:
!
!      VORT(X,Y) = dV/dX - dU/dY
!
!    at the center point of each of the bilinear elements.
!
!        |       |       |
!      (3,1)---(3,2)---(3,3)---
!        |       |       |
!        | [2,1] | [2,2] |
!        |       |       |
!      (2,1)---(2,2)---(2,3)---
!        |       |       |
!        | [1,1] | [1,2] |
!        |       |       |
!      (1,1)---(1,2)---(1,3)---
!
!    Here, the nodes labeled with parentheses represent the points at
!    which the original (U,V) data is given, while the nodes labeled
!    with square brackets represent the centers of the bilinear
!    elements, where the approximations to the divergence and vorticity
!    are made.
!
!    The reason for evaluating the divergence and vorticity in this way
!    is that the bilinear interpolant to the (U,V) data is not
!    differentiable at the boundaries of the elements, nor especially at
!    the nodes, but is an (infinitely differentiable) bilinear function
!    in the interior of each element.  If a value at the original nodes
!    is strongly desired, then the average at the four surrounding
!    central nodes may be taken.
!
!  Reference Element Q4:
!
!    |
!    1  4-----3
!    |  |     |
!    |  |     |
!    S  |     |
!    |  |     |
!    |  |     |
!    0  1-----2
!    |
!    +--0--R--1-->
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 June 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of data rows.  
!    M must be at least 2.
!
!    Input, integer ( kind = 4 ) N, the number of data columns.  
!    N must be at least 2.
!
!    Input, real ( kind = 8 ) U(M,N), V(M,N), the value of the components 
!    of a vector quantity whose divergence and vorticity are desired. 
!    A common example would be that U and V are the horizontal and 
!    vertical velocity component of a flow field.
!
!    Input, real ( kind = 8 ) Q(2,4), the coordinates of the nodes of
!    the quadrilateral, in counterclockwise order.
!
!    Output, real ( kind = 8 ) DIV(M-1,N-1), an estimate for
!    the divergence in the bilinear element that lies between
!    data rows I and I+1, and data columns J and J+1.
!
!    Output, real ( kind = 8 ) VORT(M-1,N-1), an estimate for
!    the vorticity in the bilinear element that lies between
!    data rows I and I+1, and data columns J and J+1.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) div(m-1,n-1)
  real    ( kind = 8 ) dphidx(4)
  real    ( kind = 8 ) dphidy(4)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ), parameter :: node_num = 1
  real    ( kind = 8 ) p(2)
  real    ( kind = 8 ) phi(4)
  real    ( kind = 8 ) q(2,4)
  real    ( kind = 8 ) q2(2,4)
  real    ( kind = 8 ) u(m,n)
  real    ( kind = 8 ) v(m,n)
  real    ( kind = 8 ) vort(m-1,n-1)

  if ( m <= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DIV_Q4 - Fatal error!'
    write ( *, '(a)' ) '  M must be at least 2,'
    write ( *, '(a,i8)' ) '  but the input value of M is ', m
    stop
  end if

  if ( n <= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DIV_Q4 - Fatal error!'
    write ( *, '(a)' ) '  N must be at least 2,'
    write ( *, '(a,i8)' ) '  but the input value of N is ', n
    stop
  end if

  do i = 1, m-1

    q2(2,1) =   ( real ( 2 * m - 2 * i,     kind = 8 ) * q(2,1)   &
                + real (         2 * i - 2, kind = 8 ) * q(2,3) ) &
                / real ( 2 * m         - 2, kind = 8 )
    p(2) =      ( real ( 2 * m - 2 * i - 1, kind = 8 ) * q(2,1)   &
                + real (         2 * i - 1, kind = 8 ) * q(2,3) ) &
                / real ( 2 * m         - 2, kind = 8 )
    q2(2,3) =   ( real ( 2 * m - 2 * i - 2, kind = 8 ) * q(2,1)   &
                + real (         2 * i,     kind = 8 ) * q(2,3) ) &
                / real ( 2 * m         - 2, kind = 8 )

    q2(2,2) = q2(2,1)
    q2(2,4) = q2(2,3)

    do j = 1, n-1

      q2(1,1) =   ( real ( 2 * n - 2 * j,     kind = 8 ) * q(1,1)   &
                  + real (         2 * j - 2, kind = 8 ) * q(1,3) ) &
                  / real ( 2 * n         - 2, kind = 8 )
      p(1) =      ( real ( 2 * n - 2 * j - 1, kind = 8 ) * q(1,1)   &
                  + real (         2 * j - 1, kind = 8 ) * q(1,3) ) &
                  / real ( 2 * n         - 2, kind = 8 )
      q2(1,3) =   ( real ( 2 * n - 2 * j - 2, kind = 8 ) * q(1,1)   &
                  + real (         2 * j,     kind = 8 ) * q(1,3) ) &
                  / real ( 2 * n         - 2, kind = 8 )

      q2(1,2) = q2(1,3)
      q2(1,4) = q2(1,1)

      call basis_mn_q4 ( q2, node_num, p, phi, dphidx, dphidy )
!
!  Note the following formula for the value of U and V at the same
!  point that the divergence and vorticity are being evaluated.
!
!         umid =  u(i  ,j  ) * phi(1) &
!               + u(i  ,j+1) * phi(2) &
!               + u(i+1,j+1) * phi(3) &
!               + u(i+1,j  ) * phi(4) 
!
!         vmid =  v(i  ,j  ) * phi(1) &
!               + v(i  ,j+1) * phi(2) &
!               + v(i+1,j+1) * phi(3) &
!               + v(i+1,j  ) * phi(4) 
!
      div(i,j) =  u(i  ,j  ) * dphidx(1) + v(i  ,j  ) * dphidy(1) &
                + u(i  ,j+1) * dphidx(2) + v(i  ,j+1) * dphidy(2) &
                + u(i+1,j+1) * dphidx(3) + v(i+1,j+1) * dphidy(3) &
                + u(i+1,j  ) * dphidx(4) + v(i+1,j  ) * dphidy(4) 

      vort(i,j) =  v(i  ,j  ) * dphidx(1) - u(i  ,j  ) * dphidy(1) &
                 + v(i  ,j+1) * dphidx(2) - u(i  ,j+1) * dphidy(2) &
                 + v(i+1,j+1) * dphidx(3) - u(i+1,j+1) * dphidy(3) &
                 + v(i+1,j  ) * dphidx(4) - u(i+1,j  ) * dphidy(4) 
                  
    end do
  end do

  return
end
subroutine div_t3 ( m, n, u, v, q, div, vor )

!*****************************************************************************80
!
!! DIV_T3 estimates the divergence and vorticity of a discrete field.
!
!  Discussion:
!
!    The routine is given the values of a vector field ( U(X,Y), V(X,Y) ) at
!    a regularly spaced grid of points ( X(1:M), Y(1:N) ).  This grid is 
!    described implicitly by giving the values M, N, and the coordinates
!    Q(2,4) of the bounding quadrilateral.  (Note that Q need not be a 
!    rectangle.)
!
!    The quadrilateral is suggested by the following diagram:
!
!     ^  Q(1:2,4)-----Q(1:2,3)
!     |      |            |
!     N      |            |
!     |      |            |
!     V  Q(1:2,1)-----Q(1:2,2)
!
!              <--(M)--->
!
!    The routine models the vector field over the interior of this region using
!    a linear interpolant over 2*(M-1)*(N-1) triangles.  It then uses the 
!    interpolant to estimate the value of the divergence:
!
!      DIV(X,Y) = dU/dX + dV/dY
!
!    and the vorticity:
!
!      VOR(X,Y) = dV/dX - dU/dY
!
!    at the centroid of each of the triangular elements.
!
!    The grid is (somewhat arbitrarily) subdivided into triangular elements
!    as suggested here:
!
!      (3,1)---(3,2)---(3,3)
!        | \     |  \    |
!        |  \    |   \   |
!        |   \   |    \  |
!        |    \  |     \ |
!      (2,1)---(2,2)---(2,3)
!        | \     |  \    |
!        |  \    |   \   |
!        |   \   |    \  |
!        |    \  |     \ |
!      (1,1)---(1,2)---(1,3)
!
!    In each triangular element, linear functions are used to interpolate
!    the U and V data.  The divergence and vorticity functions are then
!    evaluated at the centroid of each element.
!
!    This means that, given a grid of M X coordinates and N Y coordinates,
!    we will construct 2 * ( M - 1 ) * ( N - 1 ) triangular elements.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 November 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the number of data rows.  M must be at least 2.
!
!    Input, integer N, the number of data columns.  N must be at least 2.
!
!    Input, real ( kind = 8 ) U(M,N), V(M,N), the value of the components 
!    of a vector quantity whose divergence and vorticity are desired. 
!    A common example would be that U and V are the horizontal and 
!    vertical velocity component of a flow field.
!
!    Input, real ( kind = 8 ) Q(2,4), the coordinates of the nodes of
!    the quadrilateral, in counterclockwise order.
!
!    Output, real ( kind = 8 ) DIV(2,M-1,N-1), an estimate for
!    the divergence in the two linear elements that lie between
!    data rows I and I+1, and data columns J and J+1.
!
!    Output, real ( kind = 8 ) VOR(2,M-1,N-1), an estimate for
!    the vorticity in the two linear elements that lie between
!    data rows I and I+1, and data columns J and J+1.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) div(2,m-1,n-1)
  real    ( kind = 8 ) dphidx(3)
  real    ( kind = 8 ) dphidy(3)
  real    ( kind = 8 ) dudx
  real    ( kind = 8 ) dudy
  real    ( kind = 8 ) dvdx
  real    ( kind = 8 ) dvdy
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real    ( kind = 8 ) p(2)
  real    ( kind = 8 ) phi(3)
  real    ( kind = 8 ) q(2,4)
  real    ( kind = 8 ) t(2,3)
  real    ( kind = 8 ) u(m,n)
  real    ( kind = 8 ) v(m,n)
  real    ( kind = 8 ) vor(2,m-1,n-1)
  real    ( kind = 8 ) xlb
  real    ( kind = 8 ) xlt
  real    ( kind = 8 ) xrb
  real    ( kind = 8 ) xrt
  real    ( kind = 8 ) xxlb
  real    ( kind = 8 ) xxlt
  real    ( kind = 8 ) xxrb
  real    ( kind = 8 ) xxrt
  real    ( kind = 8 ) ylb
  real    ( kind = 8 ) ylt
  real    ( kind = 8 ) yrb
  real    ( kind = 8 ) yrt
  real    ( kind = 8 ) yylb
  real    ( kind = 8 ) yylt
  real    ( kind = 8 ) yyrb
  real    ( kind = 8 ) yyrt

  if ( m <= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DIV_T3 - Fatal error!'
    write ( *, '(a)' ) '  M must be at least 2,'
    write ( *, '(a,i8)' ) '  but the input value of M is ', m
    stop
  end if

  if ( n <= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DIV_T3 - Fatal error!'
    write ( *, '(a)' ) '  N must be at least 2,'
    write ( *, '(a,i8)' ) '  but the input value of N is ', n
    stop
  end if
!
!  Consider the data between logical rows I and I + 1.
!
  do i = 1, m - 1
!
!  Consider the data between logical columns J and J + 1.
!
    do j = 1, n - 1

      xlb = q(1,1)
      ylb = q(2,1)
      xrb = q(1,2)
      yrb = q(2,2)
      xrt = q(1,3)
      yrt = q(2,3)
      xlt = q(1,4)
      ylt = q(2,4)
      
      yylb = &
        ( &
            real ( n - j,     kind = 8 )  * ( real ( m - i,     kind = 8 ) * ylb   &
                                            + real (     i - 1, kind = 8 ) * yrb ) &
                                            / real ( m     - 1, kind = 8 )         &
          + real (     j - 1, kind = 8 )  * ( real ( m - i,     kind = 8 ) * ylt   &
                                            + real (     i - 1, kind = 8 ) * yrt ) &
                                            / real ( m     - 1, kind = 8 )         &
        ) / real ( n     - 1, kind = 8 )

      yyrb = ( &
               real ( n - j,     kind = 8 )  * ( real ( m - i - 1, kind = 8 ) * ylb   &
                                               + real (     i,     kind = 8 ) * yrb ) &
                                               / real ( m     - 1, kind = 8 )         &
             + real (     j - 1, kind = 8 )  * ( real ( m - i - 1, kind = 8 ) * ylt   &
                                               + real (     i,     kind = 8 ) * yrt ) &
                                               / real ( m     - 1, kind = 8 )         &
           ) / real ( n     - 1, kind = 8 )

      yylt = ( &
               real ( n - j - 1, kind = 8 )  * ( real ( m - i,     kind = 8 ) * ylb   &
                                               + real (     i - 1, kind = 8 ) * yrb ) &
                                               / real ( m     - 1, kind = 8 )         &
             + real (     j,     kind = 8 )  * ( real ( m - i,     kind = 8 ) * ylt   &
                                               + real (     i - 1, kind = 8 ) * yrt ) &
                                               / real ( m     - 1, kind = 8 )         &
           ) / real ( n     - 1, kind = 8 )       

      yyrt = ( &
               real ( n - j - 1, kind = 8 )  * ( real ( m - i - 1, kind = 8 ) * ylb   &
                                               + real (     i,     kind = 8 ) * yrb ) &
                                               / real ( m     - 1, kind = 8 )         &
             + real (     j,     kind = 8 )  * ( real ( m - i - 1, kind = 8 ) * ylt   &
                                               + real (     i,     kind = 8 ) * yrt ) &
                                               / real ( m     - 1, kind = 8 )         &
           ) / real ( n     - 1, kind = 8 )       

      xxlb = ( &
               real ( n - j,     kind = 8 )  * ( real ( m - i,     kind = 8 ) * xlb   &
                                               + real (     i - 1, kind = 8 ) * xrb ) &
                                               / real ( m     - 1, kind = 8 )         &
             + real (     j - 1, kind = 8 )  * ( real ( m - i,     kind = 8 ) * xlt   &
                                               + real (     i - 1, kind = 8 ) * xrt ) &
                                               / real ( m     - 1, kind = 8 )         &
           ) / real ( n     - 1, kind = 8 )

      xxlt = ( &
               real ( n - j - 1, kind = 8 )  * ( real ( m - i,     kind = 8 ) * xlb   &
                                               + real (     i - 1, kind = 8 ) * xrb ) &
                                               / real ( m     - 1, kind = 8 )         &
             + real (     j,     kind = 8 )  * ( real ( m - i,     kind = 8 ) * xlt   &
                                               + real (     i - 1, kind = 8 ) * xrt ) &
                                               / real ( m     - 1, kind = 8 )         &
           ) / real ( n     - 1, kind = 8 )

      xxrb = ( &
               real ( n - j,     kind = 8 )  * ( real ( m - i - 1, kind = 8 ) * xlb   &
                                               + real (     i,     kind = 8 ) * xrb ) &
                                               / real ( m     - 1, kind = 8 )         & 
             + real (     j - 1, kind = 8 )  * ( real ( m - i - 1, kind = 8 ) * xlt   &
                                               + real (     i,     kind = 8 ) * xrt ) &
                                               / real ( m     - 1, kind = 8 )         &
           ) / real ( n     - 1, kind = 8 )

      xxrt = ( &
               real ( n - j - 1, kind = 8 )  * ( real ( m - i - 1, kind = 8 ) * xlb   &
                                               + real (     i,     kind = 8 ) * xrb ) &
                                               / real ( m     - 1, kind = 8 )         & 
             + real (     j,     kind = 8 )  * ( real ( m - i - 1, kind = 8 ) * xlt   &
                                               + real (     i,     kind = 8 ) * xrt ) &
                                               / real ( m     - 1, kind = 8 )         &
           ) / real ( n     - 1, kind = 8 )

!     write(*,'(i4,i4,8f8.3)')i,j,xxlb, yylb, xxrb, yyrb, xxrt,yyrt, xxlt,yylt
!
!  (I,J+1) = LT-----RT = (I+1,J+1)
!            |\      |
!            | \  T2 |
!            |  \    |
!            |   \   |
!            | T1 \  |
!            |     \ |
!  (I,J)   = LB-----RB = (I+1,J)
!
      t(1:2,1:3) = reshape ( (/ xxlb, yylb, xxrb, yyrb, xxrt, yyrt /), (/ 2, 3 /) )
      p(1:2) = (/ xxlb + xxrb + xxrt, yylb + yyrb + yyrt /) / 3.0D+00
      call basis_mn_t3 ( t, 1, p, phi, dphidx, dphidy )
      dudx = u(i,j) * dphidx(1) + u(i+1,j) * dphidx(2) + u(i+1,j+1) * dphidx(3)
      dudy = u(i,j) * dphidy(1) + u(i+1,j) * dphidy(2) + u(i+1,j+1) * dphidy(3)
      dvdx = v(i,j) * dphidx(1) + v(i+1,j) * dphidx(2) + v(i+1,j+1) * dphidx(3)
      dvdy = v(i,j) * dphidy(1) + v(i+1,j) * dphidy(2) + v(i+1,j+1) * dphidy(3)

      div(1,i,j) = dudx + dvdy
      vor(1,i,j) = dvdx - dudy

!     write ( *, '(4g14.6)' ) p(1), p(2), div(1,i,j), vor(1,i,j)

      t = reshape ( (/ xxrt, yyrt, xxlt, yylt, xxlb, yylb /), (/ 2, 3 /) )
      p(1:2) = (/ xxrt + xxlt + xxlb, yyrt + yylt + yyrb /) / 3.0D+00
      call basis_mn_t3 ( t, 1, p, phi, dphidx, dphidy )
      dudx = u(i+1,j+1) * dphidx(1) + u(i,j+1) * dphidx(2) + u(i,j) * dphidx(3)
      dudy = u(i+1,j+1) * dphidy(1) + u(i,j+1) * dphidy(2) + u(i,j) * dphidy(3)
      dvdx = v(i+1,j+1) * dphidx(1) + v(i,j+1) * dphidx(2) + v(i,j) * dphidx(3)
      dvdy = v(i+1,j+1) * dphidy(1) + v(i,j+1) * dphidy(2) + v(i,j) * dphidy(3)

      div(2,i,j) = dudx + dvdy
      vor(2,i,j) = dvdx - dudy

!     write ( *, '(4g14.6)' ) p(1), p(2), div(2,i,j), vor(2,i,j)

    end do
  end do

  return
end
function element_code ( i )

!*****************************************************************************80
!
!! ELEMENT_CODE returns the code for each element.
!
!  Discussion:
!
!     I  ELEMENT_CODE   Definition
!     -  ------------   ----------
!     1  Q4             4 node linear Lagrange/serendipity quadrilateral;
!     2  Q8             8 node quadratic serendipity quadrilateral;
!     3  Q9             9 node quadratic Lagrange quadrilateral;
!     4  Q12            12 node cubic serendipity quadrilateral;
!     5  Q16            16 node cubic Lagrange quadrilateral;
!     6  QL             6 node linear/quadratic quadrilateral;
!     7  T3             3 node linear triangle;
!     8  T4             4 node cubic bubble triangle
!     9  T6             6 node quadratic triangle;
!    10  T10            10 node cubic triangle.
! 
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the index of the element type.  
!
!    Output, character ( len = 3 ) ELEMENT_CODE, the code for the element type.
!
  implicit none

  character ( len = 3 )  element_code
  integer   ( kind = 4 ) i

  if ( i == 1 ) then
    element_code = 'Q4'
  else if ( i == 2 ) then
    element_code = 'Q8'
  else if ( i == 3 ) then
    element_code = 'Q9'
  else if ( i == 4 ) then
    element_code = 'Q12'
  else if ( i == 5 ) then
    element_code = 'Q16'
  else if ( i == 6 ) then
    element_code = 'QL'
  else if ( i == 7 ) then
    element_code = 'T3'
  else if ( i == 8 ) then
    element_code = 'T4'
  else if ( i == 9 ) then
    element_code = 'T6'
  else if ( i == 10 ) then
    element_code = 'T10'
  else
    element_code = '???'
  end if

  return
end
subroutine elements_eps ( file_name, node_num, node_xy, code, &
  element_order, element_num, element_mask, element_node, title )

!*****************************************************************************80
!
!! ELEMENTS_EPS creates an EPS file image of the elements of a grid.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 April 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the file to create.
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, real ( kind = 8 ) NODE_XY(2,NODE_NUM), the
!    coordinates of the nodes.
!
!    Input, character ( len = * ) CODE, the code for the element.
!
!    Input, integer ( kind = 4 ) ELEMENT_ORDER, the element order.
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of elements.
!
!    Input, logical ELEMENT_MASK(ELEMENT_NUM), a mask for the elements.
!
!    Input, integer ( kind = 4 ) ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM), the 
!    element->node data.
!
!    Input, character ( len = * ) TITLE, a title for the plot.
!
  implicit none

  integer   ( kind = 4 ) element_num
  integer   ( kind = 4 ) element_order
  integer   ( kind = 4 ) node_num

  real ( kind = 8 ) ave_x
  real ( kind = 8 ) ave_y
  integer   ( kind = 4 ), parameter :: circle_size = 3
  character ( len = * )  code
  real ( kind = 8 ) dif
  integer   ( kind = 4 ) element
  logical                element_mask(element_num)
  integer   ( kind = 4 ) element_node(element_order,element_num)
  integer   ( kind = 4 ) eps_unit
  integer   ( kind = 4 ) eps_x
  integer   ( kind = 4 ) eps_y
  character ( len = * )  file_name
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) ios
  integer   ( kind = 4 ) j
  integer   ( kind = 4 ) local
  integer   ( kind = 4 ) next_boundary_node
  integer   ( kind = 4 ) node
  logical                node_mask(node_num)
  real ( kind = 8 ) node_x_max
  real ( kind = 8 ) node_x_min
  real ( kind = 8 ) node_xy(2,node_num)
  real ( kind = 8 ) node_y_max
  real ( kind = 8 ) node_y_min
  real ( kind = 8 ) r8_huge
  real ( kind = 8 ) scale
  character ( len = 40 ) string
  character ( len = * )  title
!
!  Determine the range of the unmasked elements.
!
  node_x_min =  r8_huge ( node_x_min )
  node_x_max = -r8_huge ( node_x_max )
  node_y_min =  r8_huge ( node_y_min )
  node_y_max = -r8_huge ( node_y_max )

  node_mask(1:node_num) = .false.

  do element = 1, element_num
    if ( element_mask(element) ) then
      do j = 1, element_order
        node = element_node(j,element)
        node_mask(node) = .true.
        node_x_min = min ( node_x_min, node_xy(1,node) )
        node_x_max = max ( node_x_max, node_xy(1,node) )
        node_y_min = min ( node_y_min, node_xy(2,node) )
        node_y_max = max ( node_y_max, node_xy(2,node) )
      end do
    end if
  end do

  if ( node_y_max - node_y_min < node_x_max - node_x_min ) then
    scale = node_x_max - node_x_min
    dif = ( node_x_max - node_x_min ) - ( node_y_max - node_y_min )
    node_y_max = node_y_max + 0.5D+00 * dif
    node_y_min = node_y_min - 0.5D+00 * dif
  else
    scale = node_y_max - node_y_min
    dif = ( node_y_max - node_y_min ) - ( node_x_max - node_x_min )
    node_x_max = node_x_max + 0.5D+00 * dif
    node_x_min = node_x_min - 0.5D+00 * dif
  end if

  call get_unit ( eps_unit )

  open ( unit = eps_unit, file = file_name, status = 'replace', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ELEMENTS_EPS - Fatal error!'
    write ( *, '(a)' ) '  Could not open the output EPS file.'
    stop
  end if

  write ( eps_unit, '(a)' ) '%!PS-Adobe-3.0 EPSF-3.0'
  write ( eps_unit, '(a)' ) '%%Creator: elements_eps(fempack.f90)'
  write ( eps_unit, '(a)' ) '%%Title: ' // trim ( file_name )
  write ( eps_unit, '(a)' ) '%%Pages: 1'
  write ( eps_unit, '(a)' ) '%%BoundingBox:    36    36   576   756'
  write ( eps_unit, '(a)' ) '%%Document-Fonts: Times-Roman'
  write ( eps_unit, '(a)' ) '%%LanguageLevel: 1'
  write ( eps_unit, '(a)' ) '%%EndComments'
  write ( eps_unit, '(a)' ) '%%BeginProlog'
  write ( eps_unit, '(a)' ) '/inch {72 mul} def'
  write ( eps_unit, '(a)' ) '%%EndProlog'
  write ( eps_unit, '(a)' ) '%%Page:      1     1'
  write ( eps_unit, '(a)' ) 'save'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '% Set RGB line color.'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) ' 0.9000 0.9000 0.9000 setrgbcolor'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '% Draw a gray border around the page.'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) 'newpath'
  write ( eps_unit, '(a)' ) '    36   126 moveto'
  write ( eps_unit, '(a)' ) '   576   126 lineto'
  write ( eps_unit, '(a)' ) '   576   666 lineto'
  write ( eps_unit, '(a)' ) '    36   666 lineto'
  write ( eps_unit, '(a)' ) '    36   126 lineto'
  write ( eps_unit, '(a)' ) 'stroke'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '% Set RGB line color.'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) ' 0.0000 0.0000 0.0000 setrgbcolor'

  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '%  Label the plot:'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) ' 0.0000 0.0000 0.0000 setrgbcolor'
  write ( eps_unit, '(a)' ) '/Times-Roman findfont 0.50 inch scalefont setfont'
  write ( eps_unit, '(a)' ) '    36   666 moveto'
  write ( eps_unit, '(a)' ) '(' // trim ( title ) // ') show'

  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '% Define a clipping polygon'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '    36   126 moveto'
  write ( eps_unit, '(a)' ) '   576   126 lineto'
  write ( eps_unit, '(a)' ) '   576   666 lineto'
  write ( eps_unit, '(a)' ) '    36   666 lineto'
  write ( eps_unit, '(a)' ) '    36   126 lineto'
  write ( eps_unit, '(a)' ) 'clip newpath'

  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '%  Draw filled dots at each node:'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) ' 0.0000 0.0000 0.9000 setrgbcolor'

  do node = 1, node_num

    if ( node_mask(node) ) then

      eps_x = int &
        ( ( node_x_max - node_xy(1,node)              ) *  61.0D+00   &
        + (            + node_xy(1,node) - node_x_min ) * 551.0D+00 ) &
        / scale

      eps_y = int &
        ( ( node_y_max - node_xy(2,node)              ) * 151.0D+00   &
        + (              node_xy(2,node) - node_y_min ) * 641.0D+00 ) &
        / scale

      write ( eps_unit, '(a,i4,2x,i4,2x,i4,a)' ) &
        'newpath  ', eps_x, eps_y, circle_size, ' 0 360 arc closepath fill'

    end if

  end do

  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '%  Label the nodes:'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) ' 0.0000 0.0000 1.0000 setrgbcolor'
  write ( eps_unit, '(a)' ) '/Times-Roman findfont 0.20 inch scalefont setfont'

  do node = 1, node_num

    if ( node_mask(node) ) then

      eps_x = int &
        ( ( node_x_max - node_xy(1,node)              ) *  61.0D+00   &
        + (            + node_xy(1,node) - node_x_min ) * 551.0D+00 ) &
        / scale

      eps_y = int &
        ( ( node_y_max - node_xy(2,node)              ) * 151.0D+00   &
        + (              node_xy(2,node) - node_y_min ) * 641.0D+00 ) &
        / scale

      write ( string, '(i4)' ) node
      string = adjustl ( string )

      write ( eps_unit, '(i4,2x,i4,a)' ) eps_x, eps_y+5, &
        ' moveto (' // trim ( string ) // ') show'

    end if

  end do

  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '%  Draw the element sides:'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) ' 0.9000 0.0000 0.0000 setrgbcolor'

  do element = 1, element_num

    if ( .not. element_mask(element) ) then
      cycle
    end if

    local = 1
    node = element_node(local,element)

    eps_x = int &
      ( ( node_x_max - node_xy(1,node)              ) *  61.0D+00   &
      + (            + node_xy(1,node) - node_x_min ) * 551.0D+00 ) &
      / scale

    eps_y = int &
      ( ( node_y_max - node_xy(2,node)              ) * 151.0D+00   &
      + (              node_xy(2,node) - node_y_min ) * 641.0D+00 ) &
      / scale

    write ( eps_unit, '(a,i4,2x,i4,a)' ) 'newpath ', eps_x, eps_y, ' moveto'

    do

      local = next_boundary_node ( local, code )
      node = element_node(local,element)

      eps_x = int &
        ( ( node_x_max - node_xy(1,node)              ) *  61.0D+00   &
        + (            + node_xy(1,node) - node_x_min ) * 551.0D+00 ) &
        / scale

      eps_y = int &
        ( ( node_y_max - node_xy(2,node)              ) * 151.0D+00   &
        + (              node_xy(2,node) - node_y_min ) * 641.0D+00 ) &
        / scale

      write ( eps_unit, '(i4,2x,i4,a)' ) eps_x, eps_y, ' lineto'

      if ( local == 1 ) then
        exit
      end if

    end do

    write ( eps_unit, '(a)' ) 'stroke'

  end do

  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '%  Label the elements:'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) ' 1.0000 0.0000 0.0000 setrgbcolor'
  write ( eps_unit, '(a)' ) '/Times-Roman findfont 0.30 inch scalefont setfont'

  do element = 1, element_num

    if ( .not. element_mask(element) ) then
      cycle
    end if

    ave_x = 0.0D+00
    ave_y = 0.0D+00

    do i = 1, element_order

      node = element_node(i,element)

      ave_x = ave_x + node_xy(1,node)
      ave_y = ave_y + node_xy(2,node)

    end do

    ave_x = ave_x / real ( element_order, kind = 8 )
    ave_y = ave_y / real ( element_order, kind = 8 )

    eps_x = int &
      ( ( node_x_max - ave_x              ) *  61.0D+00   &
      + (            + ave_x - node_x_min ) * 551.0D+00 ) &
      / scale

    eps_y = int &
      ( ( node_y_max - ave_y              ) * 151.0D+00   &
      + (              ave_y - node_y_min ) * 641.0D+00 ) &
      / scale

    write ( string, '(i4)' ) element
    string = adjustl ( string )

    write ( eps_unit, '(i4,2x,i4,a)' ) eps_x, eps_y, ' moveto (' &
      // trim ( string ) // ') show'

  end do

  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) 'restore showpage'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '% End of page'
  write ( eps_unit, '(a)' ) '%'
  write ( eps_unit, '(a)' ) '%%Trailer'
  write ( eps_unit, '(a)' ) '%%EOF'

  close ( unit = eps_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  ELEMENTS_EPS: An encapsulated PostScript file'
  write ( *, '(a)' ) '  was created containing an image of the nodes and'
  write ( *, '(a)' ) '  elements.  The file is named "' &
    // trim ( file_name ) // '".'

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
subroutine grid_element ( code, element_order, nelemx, nelemy, element_node )

!*****************************************************************************80
!
!! GRID_ELEMENT returns the element grid associated with any available element.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) CODE, identifies the element desired.
!    Legal values include 'Q4', 'Q8', 'Q9', 'Q12', 'Q16', 'QL', 'T3', 
!    'T4', 'T6' and 'T10'.
!
!    Input, integer ( kind = 4 ) ELEMENT_ORDER, the number of nodes
!    per element.
!
!    Input, integer ( kind = 4 ) NELEMX, NELEMY, the number of quadrilaterals 
!    along the X and Y directions.  The number of elements generated will be
!    NELEMX * NELEMY for quadrilaterals, or 2 * NELEMX * NELEMY for
!    triangles.
!
!    Output, integer ( kind = 4 ) ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM), 
!    the nodes that form each element.  
!
  implicit none

  integer   ( kind = 4 ) element_order

  character ( len = * )  code
  integer   ( kind = 4 ) element_node(element_order,*)
  integer   ( kind = 4 ) nelemx
  integer   ( kind = 4 ) nelemy
  logical                s_eqi

  if ( s_eqi ( code, 'Q4' ) ) then
    call grid_q4_element ( nelemx, nelemy, element_node )
  else if ( s_eqi ( code, 'Q8' ) ) then
    call grid_q8_element ( nelemx, nelemy, element_node )
  else if ( s_eqi ( code, 'Q9' ) ) then
    call grid_q9_element ( nelemx, nelemy, element_node )
  else if ( s_eqi ( code, 'Q12' ) ) then
    call grid_q12_element ( nelemx, nelemy, element_node )
  else if ( s_eqi ( code, 'Q16' ) ) then
    call grid_q16_element ( nelemx, nelemy, element_node )
  else if ( s_eqi ( code, 'QL' ) ) then
    call grid_ql_element ( nelemx, nelemy, element_node )
  else if ( s_eqi ( code, 'T3' ) ) then
    call grid_t3_element ( nelemx, nelemy, element_node )
  else if ( s_eqi ( code, 'T4' ) ) then
    call grid_t4_element ( nelemx, nelemy, element_node )
  else if ( s_eqi ( code, 'T6' ) ) then
    call grid_t6_element ( nelemx, nelemy, element_node )
  else if ( s_eqi ( code, 'T10' ) ) then
    call grid_t10_element ( nelemx, nelemy, element_node )

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRID_ELEMENT - Fatal error!'
    write ( *, '(a)' ) '  Illegal value of CODE = "' // trim ( code ) // '".'
    stop

  end if

  return
end
subroutine grid_element_num ( code, nelemx, nelemy, element_num )

!*****************************************************************************80
!
!! GRID_ELEMENT_NUM returns the number of elements in a grid.
!
!  Discussion:
!
!    The number of elements generated will be NELEMX * NELEMY for
!    quadrilaterals, or 2 * NELEMX * NELEMY for triangles.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) CODE, identifies the element desired.
!    Legal values include 'Q4', 'Q8', 'Q9', 'Q12', 'Q16', 'QL', 'T3', 
!    'T4', 'T6' and 'T10'.
!
!    Input, integer ( kind = 4 ) NELEMX, NELEMY, the number of quadrilaterals 
!    along the X and Y directions.  
!
!    Output, integer ( kind = 4 ) ELEMENT_NUM, the number of elements in 
!    the grid.
!
  implicit none

  character ( len = * )  code
  integer   ( kind = 4 ) element_num
  integer   ( kind = 4 ) nelemx
  integer   ( kind = 4 ) nelemy
  logical                s_eqi

  if ( s_eqi ( code, 'Q4' ) ) then
    call grid_q4_element_num ( nelemx, nelemy, element_num )
  else if ( s_eqi ( code, 'Q8' ) ) then
    call grid_q8_element_num ( nelemx, nelemy, element_num )
  else if ( s_eqi ( code, 'Q9' ) ) then
    call grid_q9_element_num ( nelemx, nelemy, element_num )
  else if ( s_eqi ( code, 'Q12' ) ) then
    call grid_q12_element_num ( nelemx, nelemy, element_num )
  else if ( s_eqi ( code, 'Q16' ) ) then
    call grid_q16_element_num ( nelemx, nelemy, element_num )
  else if ( s_eqi ( code, 'QL' ) ) then
    call grid_ql_element_num ( nelemx, nelemy, element_num )
  else if ( s_eqi ( code, 'T3' ) ) then
    call grid_t3_element_num ( nelemx, nelemy, element_num )
  else if ( s_eqi ( code, 'T4' ) ) then
    call grid_t4_element_num ( nelemx, nelemy, element_num )
  else if ( s_eqi ( code, 'T6' ) ) then
    call grid_t6_element_num ( nelemx, nelemy, element_num )
  else if ( s_eqi ( code, 'T10' ) ) then
    call grid_t10_element_num ( nelemx, nelemy, element_num )

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRID_ELEMENT_NUM - Fatal error!'
    write ( *, '(a)' ) '  Illegal value of CODE = "' // trim ( code ) // '".'
    element_num = -1
    stop

  end if

  return
end
subroutine grid_node_num ( code, nelemx, nelemy, node_num )

!*****************************************************************************80
!
!! GRID_NODE_NUM returns the number of nodes in a grid.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) CODE, identifies the element desired.
!    Legal values include 'Q4', 'Q8', 'Q9', 'Q12', 'Q16', 'QL', 'T3', 
!    'T4', 'T6' and 'T10'.
!
!    Input, integer ( kind = 4 ) NELEMX, NELEMY, the number of quadrilaterals 
!    along the X and Y directions.  
!
!    Output, integer ( kind = 4 ) NODE_NUM, the number of elements in the grid.
!
  implicit none

  character ( len = * )  code
  integer   ( kind = 4 ) node_num
  integer   ( kind = 4 ) nelemx
  integer   ( kind = 4 ) nelemy
  logical                s_eqi

  if ( s_eqi ( code, 'Q4' ) ) then
    call grid_q4_node_num ( nelemx, nelemy, node_num )
  else if ( s_eqi ( code, 'Q8' ) ) then
    call grid_q8_node_num ( nelemx, nelemy, node_num )
  else if ( s_eqi ( code, 'Q9' ) ) then
    call grid_q9_node_num ( nelemx, nelemy, node_num )
  else if ( s_eqi ( code, 'Q12' ) ) then
    call grid_q12_node_num ( nelemx, nelemy, node_num )
  else if ( s_eqi ( code, 'Q16' ) ) then
    call grid_q16_node_num ( nelemx, nelemy, node_num )
  else if ( s_eqi ( code, 'QL' ) ) then
    call grid_ql_node_num ( nelemx, nelemy, node_num )
  else if ( s_eqi ( code, 'T3' ) ) then
    call grid_t3_node_num ( nelemx, nelemy, node_num )
  else if ( s_eqi ( code, 'T4' ) ) then
    call grid_t4_node_num ( nelemx, nelemy, node_num )
  else if ( s_eqi ( code, 'T6' ) ) then
    call grid_t6_node_num ( nelemx, nelemy, node_num )
  else if ( s_eqi ( code, 'T10' ) ) then
    call grid_t10_node_num ( nelemx, nelemy, node_num )

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRID_NODE_NUM - Fatal error!'
    write ( *, '(a)' ) '  Illegal value of CODE = "' // trim ( code ) // '".'
    node_num = -1
    stop

  end if

  return
end
subroutine grid_nodes_01 ( x_num, y_num, node_xy )

!*****************************************************************************80
!
!! GRID_NODES_01 returns an equally spaced rectangular grid in the unit square.
!
!  Example:
!
!    X_NUM = 5
!    Y_NUM = 3
!
!    NODE_XY = 
!    ( 0, 0.25, 0.5, 0.75, 1, 0,   0.25, 0.5, 0.75, 1,   0, 0.25, 0.5, 0.75, 1;
!      0, 0,    0,   0,    0, 0.5, 0.5,  0.5, 0.5,  0.5, 1, 1.0,  1.0, 1.0,  1 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 May 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) X_NUM, Y_NUM, the number of nodes in the 
!    X and Y directions.
!
!    Output, real ( kind = 8 ) NODE_XY(2,X_NUM*Y_NUM), the coordinates of
!    the nodes.
!
  implicit none

  integer ( kind = 4 ) x_num
  integer ( kind = 4 ) y_num

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) node_num
  real    ( kind = 8 ) node_xy(2,x_num*y_num)

  node_num = x_num * y_num

  node_xy(1:2,1:node_num) = 0.0D+00

  if ( x_num == 1 ) then
    node_xy(1,1:node_num) = 0.5D+00
  else
    do i = 1, x_num
      node_xy(1,i:i+(y_num-1)*x_num:x_num) = real ( i     - 1, kind = 8 ) &
                                           / real ( x_num - 1, kind = 8 )
    end do
  end if

  if ( y_num == 1 ) then
    node_xy(2,1:node_num) = 0.5D+00
  else
    do j = 1, y_num
      node_xy(2,1+(j-1)*x_num:j*x_num) = real ( j     - 1, kind = 8 ) &
                                       / real ( y_num - 1, kind = 8 )
    end do
  end if

  return
end
subroutine grid_print ( element_order, element_num, element_node )

!*****************************************************************************80
!
!! GRID_PRINT prints the elements that form a grid.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ELEMENT_ORDER, the number of nodes 
!    per element.
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of elements.
! 
!    Input, integer ( kind = 4 ) ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM), 
!    the nodes that form each element.
!
  implicit none

  integer ( kind = 4 ) element_order
  integer ( kind = 4 ) element_num

  integer ( kind = 4 ) element
  integer ( kind = 4 ) element_node(element_order,element_num)
  integer ( kind = 4 ) i

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  GRID_PRINT: Element -> Node table.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of elements = ', element_num
  write ( *, '(a,i8)' ) '  Element order      = ', element_order
  write ( *, '(a)' ) ' '
  write ( *, '(2x,a,3x,20i3)' ) '  #', ( i, i = 1, element_order )
  write ( *, '(a)' ) ' '

  do element = 1, element_num
    write ( *, '(2x,i3,3x,20i3)' ) &
      element, element_node(1:element_order,element)
  end do

  return
end
subroutine grid_q4_element ( nelemx, nelemy, element_node )

!*****************************************************************************80
!
!! GRID_Q4_ELEMENT produces an element grid of 4 node quadrilaterals.
!
!  Discussion:
!
!    For each element, the nodes are listed in counter-clockwise order.
!
!  Example:
!
!    Input:
!
!      NELEMX = 3, NELEMY = 2
!
!    Output:
!
!      ELEMENT_NODE = 
!         1, 2,  6,  5;
!         2, 3,  7,  6;
!         3, 4,  8,  7;
!         5, 6, 10,  9;
!         6, 7, 11, 10;
!         7, 8, 12, 11.
!
!  Grid:
!
!    9---10---11---12
!    |    |    |    |
!    |    |    |    |
!    |  4 |  5 |  6 |
!    |    |    |    |
!    5----6----7----8
!    |    |    |    |
!    |    |    |    |
!    |  1 |  2 |  3 |
!    |    |    |    |
!    1----2----3----4
!
!  Reference Element Q4:
!
!    |
!    1  4-----3
!    |  |     |
!    |  |     |
!    S  |     |
!    |  |     |
!    |  |     |
!    0  1-----2
!    |
!    +--0--R--1-->
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 April 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NELEMX, NELEMY, the number of elements along 
!    the X and Y directions.  The number of elements generated will be
!    NELEMX * NELEMY.
!
!    Output, integer ( kind = 4 ) ELEMENT_NODE(4,NELEMX*NELEMY), the nodes 
!    that form each element.
!
  implicit none

  integer ( kind = 4 ) nelemx
  integer ( kind = 4 ) nelemy
  integer ( kind = 4 ), parameter :: element_order = 4

  integer ( kind = 4 ) element
  integer ( kind = 4 ) element_node(element_order,nelemx*nelemy)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ne
  integer ( kind = 4 ) nw
  integer ( kind = 4 ) se
  integer ( kind = 4 ) sw
!
!  Node labeling:
!
!    NW---NE
!     |    |
!    SW---SE
!
  element = 0
 
  do j = 1, nelemy
    do i = 1, nelemx

      sw = i     + ( j - 1 ) * ( nelemx + 1 )
      se = i + 1 + ( j - 1 ) * ( nelemx + 1 )
      nw = i     +   j       * ( nelemx + 1 )
      ne = i + 1 +   j       * ( nelemx + 1 )
  
      element = element + 1

      element_node(1,element) = sw
      element_node(2,element) = se
      element_node(3,element) = ne
      element_node(4,element) = nw

    end do
  end do

  return
end
subroutine grid_q4_element_num ( nelemx, nelemy, element_num )

!*****************************************************************************80
!
!! GRID_Q4_ELEMENT_NUM counts the elements in a grid of 4 node quadrilaterals.
!
!  Example:
!
!    Input:
!
!      NELEMX = 3, NELEMY = 2
!
!    Output:
!
!      ELEMENT_NUM = NELEMX * NELEMY = 6
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NELEMX, NELEMY, the number of elements along 
!    the X and Y directions. 
!
!    Output, integer ( kind = 4 ) ELEMENT_NUM, the number of elements in 
!    the grid.
!
  implicit none

  integer ( kind = 4 ) element_num
  integer ( kind = 4 ) nelemx
  integer ( kind = 4 ) nelemy

  element_num = nelemx * nelemy

  return
end
subroutine grid_q4_node_num ( nelemx, nelemy, node_num )

!*****************************************************************************80
!
!! GRID_Q4_NODE_NUM counts the nodes in a grid of 4 node quadrilaterals.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NELEMX, NELEMY, the number of elements along the
!    X and Y directions.
!
!    Output, integer ( kind = 4 ) NODE_NUM, the number of nodes in the grid.
!
  implicit none

  integer ( kind = 4 ) nelemx
  integer ( kind = 4 ) nelemy
  integer ( kind = 4 ) node_num

  node_num = ( nelemx + 1 ) * ( nelemy + 1 )

  return
end
subroutine grid_q8_element ( nelemx, nelemy, element_node )

!*****************************************************************************80
!
!! GRID_Q8_ELEMENT produces an element grid of 8 node quadrilaterals.
!
!  Example:
!
!    Input:
!
!      NELEMX = 3, NELEMY = 2
!
!    Output:
!
!      ELEMENT_NODE =
!         1,  3, 14, 12,  2,  9, 13,  8;
!         3,  5, 16, 14,  4, 10, 15,  9;
!         5,  7, 18, 16,  6, 11, 17, 10;
!        12, 14, 25, 23, 13, 20, 24, 19;
!        14, 16, 27, 25, 15, 21, 26, 20;
!        16, 18, 29, 27, 17, 22, 28, 21.
!
!  Diagram:
!
!   23---24---25---26---27---28---29
!    |         |         |         |
!    |         |         |         |
!   19        20        21        22
!    |         |         |         |
!    | 4       | 5       | 6       |
!   12---13---14---15---16---17---18
!    |         |         |         |
!    |         |         |         |
!    8         9        10        11
!    |         |         |         |
!    | 1       | 2       | 3       |
!    1----2----3----4----5----6----7
!
!  Reference Element Q8:
!
!    |
!    1  4--7--3
!    |  |     |
!    |  |     |
!    S  8     6
!    |  |     |
!    |  |     |
!    0  1--5--2
!    |
!    +--0--R--1-->
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 April 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NELEMX, NELEMY, the number of elements along the
!    X and Y directions.  The number of elements generated will be
!    NELEMX * NELEMY.
!
!    Output, integer ( kind = 4 ) ELEMENT_NODE(8,NELEMX*NELEMY), the nodes that form
!    each element.
!
  implicit none

  integer ( kind = 4 ) nelemx
  integer ( kind = 4 ) nelemy
  integer ( kind = 4 ), parameter :: element_order = 8

  integer ( kind = 4 ) e
  integer ( kind = 4 ) element
  integer ( kind = 4 ) element_node(element_order,nelemx*nelemy)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ne
  integer ( kind = 4 ) nw
  integer ( kind = 4 ) s
  integer ( kind = 4 ) se
  integer ( kind = 4 ) sw
  integer ( kind = 4 ) w
!
!  Node labeling:
!
!    NW----N----NE
!     |          |
!     W   (C)    E
!     |          |
!    SW----S----SE
!

  element = 0

  do j = 1, nelemy
    do i = 1, nelemx

      sw = ( j - 1 )  * ( 3 * nelemx + 2 ) + 2 * i - 1
      w  = sw + 2 * nelemx + 2 - i
      nw = sw + 3 * nelemx + 2

      s =  sw + 1
      n =  sw + ( 3 * nelemx + 2 ) + 1

      se = sw + 2
      e  = sw + 2 * nelemx + 2 - i + 1
      ne = sw + ( 3 * nelemx + 2 ) + 2

      element = element + 1

      element_node(1,element) = sw
      element_node(2,element) = se
      element_node(3,element) = ne
      element_node(4,element) = nw
      element_node(5,element) = s
      element_node(6,element) = e
      element_node(7,element) = n
      element_node(8,element) = w

    end do
  end do

  return
end
subroutine grid_q8_element_num ( nelemx, nelemy, element_num )

!*****************************************************************************80
!
!! GRID_Q8_ELEMENT_NUM counts the elements in a grid of 8 node quadrilaterals.
!
!  Example:
!
!    Input:
!
!      NELEMX = 3, NELEMY = 2
!
!    Output:
!
!      ELEMENT_NUM = NELEMX * NELEMY = 6
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NELEMX, NELEMY, the number of elements along the
!    X and Y directions.  The number of elements generated will be
!    NELEMX * NELEMY.
!
!    Output, integer ( kind = 4 ) ELEMENT_NUM, the number of elements in the grid.
!
  implicit none

  integer ( kind = 4 ) element_num
  integer ( kind = 4 ) nelemx
  integer ( kind = 4 ) nelemy

  element_num = nelemx * nelemy

  return
end
subroutine grid_q8_node_num ( nelemx, nelemy, node_num )

!*****************************************************************************80
!
!! GRID_Q8_NODE_NUM counts the nodes in a grid of 8 node quadrilaterals.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NELEMX, NELEMY, the number of elements along the
!    X and Y directions.
!
!    Output, integer ( kind = 4 ) NODE_NUM, the number of node in the grid.
!
  implicit none

  integer ( kind = 4 ) nelemx
  integer ( kind = 4 ) nelemy
  integer ( kind = 4 ) node_num

  node_num = 3 * nelemx * nelemy + 2 * nelemx + 2 * nelemy + 1

  return
end
subroutine grid_q9_element ( nelemx, nelemy, element_node )

!*****************************************************************************80
!
!! GRID_Q9_ELEMENT produces an element grid of 9 node quadrilaterals.
!
!  Example:
!
!    Input:
!
!      NELEMX = 3, NELEMY = 2
!
!    Output:
!
!      ELEMENT_NODE = 
!         1,  3, 17, 15,  2, 10, 16,  8,  9;
!         3,  5, 19, 17,  4, 12, 18, 10, 11;
!         5,  7, 21, 19,  6, 14, 20, 12, 13;
!        15, 17, 31, 29, 16, 24, 30, 22, 23;
!        17, 19, 33, 31, 18, 26, 32, 24, 25;
!        19, 21, 35, 33, 20, 28, 34, 26, 27.
!
!  Grid:
!
!   29---30---31---32---33---34---35
!    |    .    |    .    |    .    |
!    |    .    |    .    |    .    |
!   22 . 23 . 24 . 25 . 26 . 27 . 28
!    |    .    |    .    |    .    |
!    | 4  .    | 5  .    | 6  .    |
!   15---16---17---18---19---20---21
!    |    .    |    .    |    .    |
!    |    .    |    .    |    .    |
!    8 .  9 . 10 . 11 . 12 . 13 . 14
!    |    .    |    .    |    .    |
!    | 1  .    | 2  .    | 3  .    |
!    1----2----3----4----5----6----7
!
!  Reference Element Q9:
!
!    |
!    1  4--7--3
!    |  |     |
!    |  |     |
!    S  8  9  6
!    |  |     |
!    |  |     |
!    0  1--5--2
!    |
!    +--0--R--1-->
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 April 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NELEMX, NELEMY, the number of elements along 
!    the X and Y directions.  The number of elements generated will be
!    NELEMX * NELEMY.
!
!    Output, integer ( kind = 4 ) ELEMENT_NODE(9,NELEMX*NELEMY), the nodes that form
!    each element.
!
  implicit none

  integer ( kind = 4 ) nelemx
  integer ( kind = 4 ) nelemy
  integer ( kind = 4 ), parameter :: element_order = 9

  integer ( kind = 4 ) c
  integer ( kind = 4 ) e
  integer ( kind = 4 ) element
  integer ( kind = 4 ) element_node(element_order,nelemx*nelemy)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ne
  integer ( kind = 4 ) nw
  integer ( kind = 4 ) s
  integer ( kind = 4 ) se
  integer ( kind = 4 ) sw
  integer ( kind = 4 ) w
!
!  Node labeling:
!
!    NW----N----NE
!     |          |
!     W    C     E
!     |          |
!    SW----S----SE
!
  element = 0

  do j = 1, nelemy
    do i = 1, nelemx

      sw = 2 * ( j - 1 )  * ( 2 * nelemx + 1 ) + 2 * ( i - 1 ) + 1
      w  = sw +               2 * nelemx + 1
      nw = sw +         2 * ( 2 * nelemx + 1 )

      s  = sw + 1
      c  = sw + 1 +               2 * nelemx + 1
      n  = sw + 1 +         2 * ( 2 * nelemx + 1 )

      se = sw + 2
      e  = sw + 2 +               2 * nelemx + 1
      ne = sw + 2 +         2 * ( 2 * nelemx + 1 )

      element = element + 1

      element_node(1,element) = sw
      element_node(2,element) = se
      element_node(3,element) = ne
      element_node(4,element) = nw
      element_node(5,element) = s
      element_node(6,element) = e
      element_node(7,element) = n
      element_node(8,element) = w
      element_node(9,element) = c

    end do
  end do

  return
end
subroutine grid_q9_element_num ( nelemx, nelemy, element_num )

!*****************************************************************************80
!
!! GRID_Q9_ELEMENT_NUM counts the elements in a grid of 9 node quadrilaterals.
!
!  Example:
!
!    Input:
!
!      NELEMX = 3, NELEMY = 2
!
!    Output:
!
!      ELEMENT_NUM = NELEMX * NELEMY = 6
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NELEMX, NELEMY, the number of elements along the
!    X and Y directions.  The number of elements generated will be
!    NELEMX * NELEMY.
!
!    Output, integer ( kind = 4 ) ELEMENT_NUM, the number of elements in the grid.
!
  implicit none

  integer ( kind = 4 ) element_num
  integer ( kind = 4 ) nelemx
  integer ( kind = 4 ) nelemy

  element_num = nelemx * nelemy

  return
end
subroutine grid_q9_node_num ( nelemx, nelemy, node_num )

!*****************************************************************************80
!
!! GRID_Q9_NODE_NUM counts the nodes in a grid of 9 node quadrilaterals.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NELEMX, NELEMY, the number of elements along 
!    the X and Y directions. 
!
!    Output, integer ( kind = 4 ) NODE_NUM, the number of nodes in the grid.
!
  implicit none

  integer ( kind = 4 ) nelemx
  integer ( kind = 4 ) nelemy
  integer ( kind = 4 ) node_num

  node_num = ( 2 * nelemx + 1 ) * ( 2 * nelemy + 1 )

  return
end
subroutine grid_q12_element ( nelemx, nelemy, element_node )

!*****************************************************************************80
!
!! GRID_Q12_ELEMENT produces an element grid of 12 node quadrilaterals.
!
!  Example:
!
!    Input:
!
!      NELEMX = 3, NELEMY = 2
!
!    Output:
!
!      ELEMENT_NODE =
!         1,  2,  3,  4, 11, 12, 15, 16, 19, 20, 21, 22;
!         4,  5,  6,  7, 12, 13, 16, 17, 22, 23, 24, 25;
!         7,  8,  9, 10, 13, 14, 17, 18, 25, 26, 27, 28;
!        19, 20, 21, 22, 29, 30, 33, 34, 37, 38, 39, 40;
!        22, 23, 24, 25, 30, 31, 34, 35, 40, 41, 42, 43;
!        25, 26, 27, 28, 31, 32, 35, 36, 43, 44, 45, 46.
!
!  Grid:
!
!   37-38-39-40-41-42-43-44-45-46
!    |        |        |        |
!   33       34       35       36
!    |        |        |        |
!   29       30       31       32
!    | 4      | 5      | 6      |
!   19-20-21-22-23-24-25-26-27-28
!    |        |        |        |
!   15       16       17       18
!    |        |        |        |
!   11       12       13       14
!    | 1      | 2      | 3      |
!    1--2--3--4--5--6--7--8--9-10
!
!  Reference Element Q12:
!
!    |
!    1  9-10-11-12
!    |  |        |
!    |  7        8
!    S  |        |
!    |  5        6
!    |  |        |
!    0  1--2--3--4
!    |
!    +--0---R---1-->
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NELEMX, NELEMY, the number of elements along the
!    X and Y directions.  The number of elements generated will be
!    NELEMX * NELEMY.
!
!    Output, integer ( kind = 4 ) ELEMENT_NODE(12,NELEMX*NELEMY), the nodes 
!    that form each element.
!
  implicit none

  integer ( kind = 4 ) nelemx
  integer ( kind = 4 ) nelemy
  integer ( kind = 4 ), parameter :: element_order = 12

  integer ( kind = 4 ) base
  integer ( kind = 4 ) element
  integer ( kind = 4 ) element_node(element_order,nelemx*nelemy)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  element = 0

  do j = 1, nelemy
    do i = 1, nelemx

      base = ( j - 1 )  * ( 5 * nelemx + 3 ) + 1

      element = element + 1

      element_node( 1,element) = base + ( i - 1 ) * 3
      element_node( 2,element) = base + ( i - 1 ) * 3 + 1
      element_node( 3,element) = base + ( i - 1 ) * 3 + 2
      element_node( 4,element) = base + ( i - 1 ) * 3 + 3

      element_node( 5,element) = base + 3 * nelemx + i
      element_node( 6,element) = base + 3 * nelemx + i + 1

      element_node( 7,element) = base + 4 * nelemx + i + 1
      element_node( 8,element) = base + 4 * nelemx + i + 2

      element_node( 9,element) = base + 5 * nelemx + 3 * i
      element_node(10,element) = base + 5 * nelemx + 3 * i + 1
      element_node(11,element) = base + 5 * nelemx + 3 * i + 2
      element_node(12,element) = base + 5 * nelemx + 3 * i + 3

    end do
  end do

  return
end
subroutine grid_q12_element_num ( nelemx, nelemy, element_num )

!*****************************************************************************80
!
!! GRID_Q12_ELEMENT_NUM counts the elements in a grid of 12 node quadrilaterals.
!
!  Example:
!
!    Input:
!
!      NELEMX = 3, NELEMY = 2
!
!    Output:
!
!      ELEMENT_NUM = NELEMX * NELEMY = 6
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NELEMX, NELEMY, the number of elements along the
!    X and Y directions.  The number of elements generated will be
!    NELEMX * NELEMY.
!
!    Output, integer ( kind = 4 ) ELEMENT_NUM, the number of elements in 
!    the grid.
!
  implicit none

  integer ( kind = 4 ) element_num
  integer ( kind = 4 ) nelemx
  integer ( kind = 4 ) nelemy

  element_num = nelemx * nelemy

  return
end
subroutine grid_q12_node_num ( nelemx, nelemy, node_num )

!*****************************************************************************80
!
!! GRID_Q12_NODE_NUM counts the nodes in a grid of 12 node quadrilaterals.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NELEMX, NELEMY, the number of elements along the
!    X and Y directions.
!
!    Output, integer ( kind = 4 ) NODE_NUM, the number of node in the grid.
!
  implicit none

  integer ( kind = 4 ) nelemx
  integer ( kind = 4 ) nelemy
  integer ( kind = 4 ) node_num

  node_num = 5 * nelemx * nelemy + 3 * nelemx + 3 * nelemy + 1

  return
end
subroutine grid_q16_element ( nelemx, nelemy, element_node )

!*****************************************************************************80
!
!! GRID_Q16_ELEMENT produces an element grid of 16 node quadrilaterals.
!
!  Example:
!
!    Input:
!
!      NELEMX = 2, NELEMY = 2
!
!    Output:
!
!      ELEMENT_NODE = 
!         1,  2,  3,  4,  8,  9, 10, 11, 15, 16, 17, 18, 22, 23, 24, 25;
!         4,  5,  6,  7, 11, 12, 13, 14, 18, 19, 20, 21, 25, 26, 27, 28;
!        22, 23, 24, 25, 29, 30, 31, 32, 36, 37, 38, 39, 43, 44, 45, 46;
!        25, 26, 27, 28, 32, 33, 34, 35, 39, 40, 41, 42, 46, 47, 48, 49. 
!        
!  Grid:
!
!   43-44-45-46-47-48-49
!    |        |        |
!    |        |        |
!   36 37 38 39 40 41 42
!    |        |        |
!    |        |        |
!   29 30 31 32 33 34 35
!    |        |        |
!    | 3      | 4      |
!   22-23-24-25-26-27-28
!    |        |        |
!    |        |        |
!   15 16 17 18 19 20 21
!    |        |        |
!    |        |        |
!    8  9 10 11 12 13 14
!    |        |        |
!    | 1      | 2      |
!    1--2--3--4--5--6--7
!
!  Reference Element Q16:
!
!    |
!    1 13--14--15--16
!    |  |   :   :   |
!    |  |   :   :   |
!    |  9..10..11..12
!    S  |   :   :   |
!    |  |   :   :   |
!    |  5...6...7...8
!    |  |   :   :   |
!    |  |   :   :   |  
!    0  1---2---3---4
!    |
!    +--0-----R-----1-->
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NELEMX, NELEMY, the number of elements along the
!    X and Y directions.  The number of elements generated will be
!    NELEMX * NELEMY.
!
!    Output, integer ( kind = 4 ) ELEMENT_NODE(16,NELEMX*NELEMY), the nodes 
!    that form each element.
!
  implicit none

  integer ( kind = 4 ) nelemx
  integer ( kind = 4 ) nelemy
  integer ( kind = 4 ), parameter :: element_order = 16

  integer ( kind = 4 ) base
  integer ( kind = 4 ) element
  integer ( kind = 4 ) element_node(element_order,nelemx*nelemy)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  element = 0
 
  do j = 1, nelemy
    do i = 1, nelemx

      base = ( j - 1 ) * 3 * ( 3 * nelemx + 1 ) + 3 * i - 2

      element = element + 1

      element_node( 1,element) = base
      element_node( 2,element) = base                          + 1
      element_node( 3,element) = base                          + 2
      element_node( 4,element) = base                          + 3
      element_node( 5,element) = base +     ( 3 * nelemx + 1 )
      element_node( 6,element) = base +     ( 3 * nelemx + 1 ) + 1
      element_node( 7,element) = base +     ( 3 * nelemx + 1 ) + 2
      element_node( 8,element) = base +     ( 3 * nelemx + 1 ) + 3
      element_node( 9,element) = base + 2 * ( 3 * nelemx + 1 )
      element_node(10,element) = base + 2 * ( 3 * nelemx + 1 ) + 1
      element_node(11,element) = base + 2 * ( 3 * nelemx + 1 ) + 2
      element_node(12,element) = base + 2 * ( 3 * nelemx + 1 ) + 3
      element_node(13,element) = base + 3 * ( 3 * nelemx + 1 )
      element_node(14,element) = base + 3 * ( 3 * nelemx + 1 ) + 1
      element_node(15,element) = base + 3 * ( 3 * nelemx + 1 ) + 2
      element_node(16,element) = base + 3 * ( 3 * nelemx + 1 ) + 3

    end do
  end do

  return
end
subroutine grid_q16_element_num ( nelemx, nelemy, element_num )

!*****************************************************************************80
!
!! GRID_Q16_ELEMENT_NUM counts the elements in a grid of 16 node quadrilaterals.
!
!  Example:
!
!    Input:
!
!      NELEMX = 3, NELEMY = 2
!
!    Output:
!
!      ELEMENT_NUM = NELEMX * NELEMY = 6
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NELEMX, NELEMY, the number of elements along the
!    X and Y directions.  The number of elements generated will be
!    NELEMX * NELEMY.
!
!    Output, integer ( kind = 4 ) ELEMENT_NUM, the number of elements in the grid.
!
  implicit none

  integer ( kind = 4 ) element_num
  integer ( kind = 4 ) nelemx
  integer ( kind = 4 ) nelemy

  element_num = nelemx * nelemy

  return
end
subroutine grid_q16_node_num ( nelemx, nelemy, node_num )

!*****************************************************************************80
!
!! GRID_Q16_NODE_NUM counts the nodes in a grid of 16 node quadrilaterals.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NELEMX, NELEMY, the number of elements along the
!    X and Y directions.
!
!    Output, integer ( kind = 4 ) NODE_NUM, the number of nodes in the grid.
!
  implicit none

  integer ( kind = 4 ) nelemx
  integer ( kind = 4 ) nelemy
  integer ( kind = 4 ) node_num

  node_num = ( 3 * nelemx + 1 ) * ( 3 * nelemy + 1 )

  return
end
subroutine grid_ql_element ( nelemx, nelemy, element_node )

!*****************************************************************************80
!
!! GRID_QL_ELEMENT produces an element grid of 6 node quadratics/linears.
!
!  Example:
!
!    Input:
!
!      NELEMX = 3, NELEMY = 2
!
!    Output:
!
!      ELEMENT_NODE = 
!         1,  2,  3,  8,  9, 10;
!         3,  4,  5, 10, 11, 12;
!         5,  6,  7, 12, 13, 14;
!         8,  9, 10, 15, 16, 17;
!        10, 11, 12, 17, 18, 19;
!        12, 13, 14, 19, 20, 21.
!
!  Grid:
!
!   15---16---17---18---19---20---21
!    |         |         |         |
!    |         |         |         |
!    |    4    |    5    |    6    |
!    |         |         |         |
!    |         |         |         |
!    8----9---10---11---12---13---14
!    |         |         |         |
!    |         |         |         |
!    |    1    |    2    |    3    |
!    |         |         |         |
!    |         |         |         |
!    1----2----3----4----5----6----7
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NELEMX, NELEMY, the number of elements along the
!    X and Y directions.  The number of elements generated will be
!    NELEMX * NELEMY.  X will the the "quadratic direction", and
!    Y will be the "linear direction".
!
!    Output, integer ( kind = 4 ) ELEMENT_NODE(6,NELEMX*NELEMY), the nodes 
!    that form each element.
!
  implicit none

  integer ( kind = 4 ) nelemx
  integer ( kind = 4 ) nelemy
  integer ( kind = 4 ), parameter :: element_order = 6

  integer ( kind = 4 ) base
  integer ( kind = 4 ) element
  integer ( kind = 4 ) element_node(element_order,nelemx*nelemy)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  element = 0
 
  do j = 1, nelemy
    do i = 1, nelemx

      base = ( j - 1 )  * ( 2 * nelemx + 1 ) + 2 * i - 1

      element = element + 1

      element_node(1,element) = base
      element_node(2,element) = base + 1
      element_node(3,element) = base + 2
      element_node(4,element) = base + ( 2 * nelemx + 1 )
      element_node(5,element) = base + ( 2 * nelemx + 1 ) + 1
      element_node(6,element) = base + ( 2 * nelemx + 1 ) + 2

    end do
  end do

  return
end
subroutine grid_ql_element_num ( nelemx, nelemy, element_num )

!*****************************************************************************80
!
!! GRID_QL_ELEMENT_NUM counts the elements in a grid of QL quadrilaterals.
!
!  Example:
!
!    Input:
!
!      NELEMX = 3, NELEMY = 2
!
!    Output:
!
!      ELEMENT_NUM = NELEMX * NELEMY = 6
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NELEMX, NELEMY, the number of elements along the
!    X and Y directions.  The number of elements generated will be
!    NELEMX * NELEMY.
!
!    Output, integer ( kind = 4 ) ELEMENT_NUM, the number of elements in the grid.
!
  implicit none

  integer ( kind = 4 ) element_num
  integer ( kind = 4 ) nelemx
  integer ( kind = 4 ) nelemy

  element_num = nelemx * nelemy

  return
end
subroutine grid_ql_node_num ( nelemx, nelemy, node_num )

!*****************************************************************************80
!
!! GRID_QL_NODE_NUM counts the nodes in a grid of QL quadrilaterals.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NELEMX, NELEMY, the number of elements along the
!    X and Y directions.
!
!    Output, integer ( kind = 4 ) NODE_NUM, the number of nodes in the grid.
!
  implicit none

  integer ( kind = 4 ) nelemx
  integer ( kind = 4 ) nelemy
  integer ( kind = 4 ) node_num

  node_num = 2 * nelemx * nelemy + 2 * nelemx + nelemy + 1

  return
end
subroutine grid_shape_2d ( n, a, n1, n2 )

!*****************************************************************************80
!
!! GRID_SHAPE_2D guesses the shape N1 by N2 of a vector of data.
!
!  Discussion:
!
!    The data vector A is assumed to contain N1 * N2 values, with
!    where each of N2 values is repeated N1 times.
!
!  Example:
!
!    Input:
!
!      A = ( 2, 2, 2, 7, 7, 7 )
!
!    Output:
!
!      N1 = 3, N2 = 2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data values.
!
!    Input, real ( kind = 8 ) A(N), the data, which should have the properties
!    described above.
!
!    Output, integer ( kind = 4 ) N1, N2, the "shape" of the data in the array.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
!
!  Make a guess for N1.
!
  i = 1
  n1 = 1

  do i = 2, n
    if ( a(i) /= a(1) ) then
      exit
    end if
    n1 = n1 + 1
  end do
!
!  Guess that N2 = N / N1.
!
  n2 = n / n1

  return
end
subroutine grid_t3_element ( nelemx, nelemy, element_node )

!*****************************************************************************80
!
!! GRID_T3_ELEMENT produces an element grid of pairs of 3 node triangles.
!
!  Example:
!
!    Input:
!
!      NELEMX = 3, NELEMY = 2
!
!    Output:
!
!      ELEMENT_NODE = 
!         1,  2,  5;
!         6,  5,  2;
!         2,  3,  6;
!         7,  6,  3;
!         3,  4,  7;
!         8,  7,  4;
!         5,  6,  9;
!        10,  9,  6;
!         6,  7, 10;
!        11, 10,  7;
!         7,  8, 11;
!        12, 11,  8.
!
!  Grid:
!
!    9---10---11---12
!    |\ 8 |\10 |\12 |
!    | \  | \  | \  |
!    |  \ |  \ |  \ |
!    |  7\|  9\| 11\|
!    5----6----7----8
!    |\ 2 |\ 4 |\ 6 |
!    | \  | \  | \  |
!    |  \ |  \ |  \ |
!    |  1\|  3\|  5\|
!    1----2----3----4
!
!  Reference Element T3:
!
!    |
!    1  3
!    |  |\
!    |  | \
!    S  |  \
!    |  |   \
!    |  |    \
!    0  1-----2
!    |
!    +--0--R--1-->
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 April 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NELEMX, NELEMY, the number of elements along the
!    X and Y directions.  The number of elements generated will be
!    2 * NELEMX * NELEMY.
!
!    Output, integer ( kind = 4 ) ELEMENT_NODE(3,2*NELEMX*NELEMY), the nodes 
!    that form each element.
!
  implicit none

  integer ( kind = 4 ) nelemx
  integer ( kind = 4 ) nelemy
  integer ( kind = 4 ), parameter :: element_order = 3

  integer ( kind = 4 ) element
  integer ( kind = 4 ) element_node(element_order,2*nelemx*nelemy)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ne
  integer ( kind = 4 ) nw
  integer ( kind = 4 ) se
  integer ( kind = 4 ) sw
!
!  Node labeling:
!
!    NW--NE
!     |\ |
!     | \|
!    SW--SE
!
  element = 0
 
  do j = 1, nelemy
    do i = 1, nelemx

      sw = i     + ( j - 1 ) * ( nelemx + 1 )
      se = i + 1 + ( j - 1 ) * ( nelemx + 1 )
      nw = i     +   j       * ( nelemx + 1 )
      ne = i + 1 +   j       * ( nelemx + 1 )

      element = element + 1

      element_node(1,element) = sw
      element_node(2,element) = se
      element_node(3,element) = nw

      element = element + 1

      element_node(1,element) = ne
      element_node(2,element) = nw
      element_node(3,element) = se

    end do
  end do

  return
end
subroutine grid_t3_element_num ( nelemx, nelemy, element_num )

!*****************************************************************************80
!
!! GRID_T3_ELEMENT_NUM counts the elements in a grid of 3 node triangles.
!
!  Example:
!
!    Input:
!
!      NELEMX = 3, NELEMY = 2
!
!    Output:
!
!      ELEMENT_NUM = 2 * NELEMX * NELEMY = 12
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NELEMX, NELEMY, the number of elements along 
!    the X and Y directions.  The number of elements generated will be
!    NELEMX * NELEMY.
!
!    Output, integer ( kind = 4 ) ELEMENT_NUM, the number of elements in the grid.
!
  implicit none

  integer ( kind = 4 ) element_num
  integer ( kind = 4 ) nelemx
  integer ( kind = 4 ) nelemy

  element_num = 2 * nelemx * nelemy

  return
end
subroutine grid_t3_node_num ( nelemx, nelemy, node_num )

!*****************************************************************************80
!
!! GRID_T3_NODE_NUM counts the nodes in a grid of 3 node triangles.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NELEMX, NELEMY, the number of elements along 
!    the X and Y directions.  The number of elements generated will be
!    NELEMX * NELEMY.
!
!    Output, integer ( kind = 4 ) NODE_NUM, the number of nodes in the grid.
!
  implicit none

  integer ( kind = 4 ) nelemx
  integer ( kind = 4 ) nelemy
  integer ( kind = 4 ) node_num

  node_num = ( nelemx + 1 ) * ( nelemy + 1 )

  return
end
subroutine grid_t4_element ( nelemx, nelemy, element_node )

!*****************************************************************************80
!
!! GRID_T4_ELEMENT produces an element grid of pairs of 4 node triangles.
!
!  Example:
!
!    Input:
!
!      NELEMX = 3, NELEMY = 2
!
!    Output:
!
!      ELEMENT_NODE = 
!         1,  2,  11,  5;
!        12, 11,   2,  8;
!         2,  3,  12,  6;
!        13, 12,   3,  9;
!         3   4   13,  7;
!        14, 13,   4,  10;
!        11, 12,  21,  15;
!        22, 21,  12,  18;
!        12, 13,  22,  16;
!        23, 22,  13,  19;
!        13  14   23,  17;
!        24, 23,  14,  20;
!
!  Grid:
!
!   21---22---23---24
!    |\18 |\19 |\20 |
!    | \  | \  | \  |
!    |  \ |  \ |  \ |
!    | 15\| 16\| 17\|
!   11---12---13---14
!    |\ 8 |\ 9 |\10 |
!    | \  | \  | \  |
!    |  \ |  \ |  \ |
!    | 5 \|  6\|  7\|
!    1----2----3----4
!
!  Reference Element T4:
!
!    |
!    1  3
!    |  |\
!    |  | \
!    S  |  \
!    |  | 4 \
!    |  |    \
!    0  1-----2
!    |
!    +--0--R--1-->
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NELEMX, NELEMY, the number of elements along the
!    X and Y directions.  The number of elements generated will be
!    2 * NELEMX * NELEMY.
!
!    Output, integer ( kind = 4 ) ELEMENT_NODE(4,2*NELEMX*NELEMY), the nodes 
!    that form each element.
!
  implicit none

  integer ( kind = 4 ) nelemx
  integer ( kind = 4 ) nelemy
  integer ( kind = 4 ), parameter :: element_order = 4

  integer ( kind = 4 ) element
  integer ( kind = 4 ) element_node(element_order,2*nelemx*nelemy)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nc
  integer ( kind = 4 ) ne
  integer ( kind = 4 ) nw
  integer ( kind = 4 ) sc
  integer ( kind = 4 ) se
  integer ( kind = 4 ) sw
!
!  Node labeling:
!
!    NW----NE
!     |\   |
!     | \NC|
!     |SC\ |
!     |   \|
!    SW---SE
!
  element = 0
 
  do j = 1, nelemy
    do i = 1, nelemx

      sw = i     + ( j - 1 ) * ( 3 * nelemx + 1 )
      se = sw + 1
      sc = sw +     nelemx + 1
      nc = sw + 2 * nelemx + 1
      nw = sw + 3 * nelemx + 1
      ne = sw + 3 * nelemx + 2

      element = element + 1
      element_node(1,element) = sw
      element_node(2,element) = se
      element_node(3,element) = nw
      element_node(4,element) = sc

      element = element + 1
      element_node(1,element) = ne
      element_node(2,element) = nw
      element_node(3,element) = se
      element_node(4,element) = nc

    end do
  end do

  return
end
subroutine grid_t4_element_num ( nelemx, nelemy, element_num )

!*****************************************************************************80
!
!! GRID_T4_ELEMENT_NUM counts the elements in a grid of 4 node triangles.
!
!  Example:
!
!    Input:
!
!      NELEMX = 3, NELEMY = 2
!
!    Output:
!
!      ELEMENT_NUM = 2 * NELEMX * NELEMY = 12
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NELEMX, NELEMY, the number of elements along 
!    the X and Y directions.  The number of elements generated will be
!    NELEMX * NELEMY.
!
!    Output, integer ( kind = 4 ) ELEMENT_NUM, the number of elements in 
!    the grid.
!
  implicit none

  integer ( kind = 4 ) element_num
  integer ( kind = 4 ) nelemx
  integer ( kind = 4 ) nelemy

  element_num = 2 * nelemx * nelemy

  return
end
subroutine grid_t4_node_num ( nelemx, nelemy, node_num )

!*****************************************************************************80
!
!! GRID_T4_NODE_NUM counts the nodes in a grid of 4 node triangles.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NELEMX, NELEMY, the number of elements along the
!    X and Y directions.
!
!    Output, integer ( kind = 4 ) NODE_NUM, the number of nodes in the grid.
!
  implicit none

  integer ( kind = 4 ) nelemx
  integer ( kind = 4 ) nelemy
  integer ( kind = 4 ) node_num

  node_num = ( nelemx + 1 ) * ( nelemy + 1 ) + 2 * nelemx * nelemy

  return
end
subroutine grid_t6_element ( nelemx, nelemy, element_node )

!*****************************************************************************80
!
!! GRID_T6_ELEMENT produces an element grid of pairs of 6 node triangles.
!
!  Example:
!
!    Input:
!
!      NELEMX = 3, NELEMY = 2
!
!    Output:
!
!      ELEMENT_NODE = 
!         1,  3, 15,  2,  9,  8;
!        17, 15,  3, 16,  9, 10;
!         3,  5, 17,  4, 11, 10;
!        19, 17,  5, 18, 11, 12;
!         5,  7, 19,  6, 13, 12;
!        21, 19,  7, 20, 13, 14;
!        15, 17, 29, 16, 23, 22;
!        31, 29, 17, 30, 23, 24;
!        17, 19, 31, 18, 25, 24;
!        33, 31, 19, 32, 25, 26;
!        19, 21, 33, 20, 27, 26;
!        35, 33, 21, 34, 27, 28.
!
!  Grid:
!
!   29-30-31-32-33-34-35
!    |\ 8  |\10  |\12  |
!    | \   | \   | \   |
!   22 23 24 25 26 27 28
!    |   \ |   \ |   \ |
!    |  7 \|  9 \| 11 \|
!   15-16-17-18-19-20-21
!    |\ 2  |\ 4  |\ 6  |
!    | \   | \   | \   |
!    8  9 10 11 12 13 14
!    |   \ |   \ |   \ |
!    |  1 \|  3 \|  5 \|
!    1--2--3--4--5--6--7
!
!  Reference Element T6:
!
!    |
!    1  3
!    |  |\
!    |  | \
!    S  6  5
!    |  |   \
!    |  |    \
!    0  1--4--2
!    |
!    +--0--R--1-->
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 April 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NELEMX, NELEMY, the number of elements along the
!    X and Y directions.  The number of elements generated will be
!    2 * NELEMX * NELEMY.
!
!    Output, integer ( kind = 4 ) ELEMENT_NODE(6,2*NELEMX*NELEMY), the nodes 
!    that form each element.
!
  implicit none

  integer ( kind = 4 ) nelemx
  integer ( kind = 4 ) nelemy
  integer ( kind = 4 ), parameter :: element_order = 6

  integer ( kind = 4 ) c
  integer ( kind = 4 ) e
  integer ( kind = 4 ) element
  integer ( kind = 4 ) element_node(element_order,2*nelemx*nelemy)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ne
  integer ( kind = 4 ) nw
  integer ( kind = 4 ) s
  integer ( kind = 4 ) se
  integer ( kind = 4 ) sw
  integer ( kind = 4 ) w
!
!  Node labeling:
!
!    NW---N--NE
!     | \     |
!     W   C   E
!     |    \  |
!    SW---S--SE
!
  element = 0
 
  do j = 1, nelemy
    do i = 1, nelemx

      sw = 2 * ( j - 1 )  * ( 2 * nelemx + 1 ) + 2 * ( i - 1 ) + 1
      w  = sw +               2 * nelemx + 1
      nw = sw +         2 * ( 2 * nelemx + 1 )

      s  = sw + 1
      c  = sw + 1 +               2 * nelemx + 1
      n  = sw + 1 +         2 * ( 2 * nelemx + 1 )

      se = sw + 2
      e  = sw + 2 +               2 * nelemx + 1
      ne = sw + 2 +         2 * ( 2 * nelemx + 1 )

      element = element + 1

      element_node(1,element) = sw
      element_node(2,element) = se
      element_node(3,element) = nw
      element_node(4,element) = s
      element_node(5,element) = c
      element_node(6,element) = w

      element = element + 1

      element_node(1,element) = ne
      element_node(2,element) = nw
      element_node(3,element) = se
      element_node(4,element) = n
      element_node(5,element) = c
      element_node(6,element) = e

    end do
  end do

  return
end
subroutine grid_t6_element_num ( nelemx, nelemy, element_num )

!*****************************************************************************80
!
!! GRID_T6_ELEMENT_NUM counts the elements in a grid of 6 node triangles.
!
!  Example:
!
!    Input:
!
!      NELEMX = 3, NELEMY = 2
!
!    Output:
!
!      ELEMENT_NUM = 2 * NELEMX * NELEMY = 12
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NELEMX, NELEMY, the number of elements along the
!    X and Y directions.  The number of elements generated will be
!    NELEMX * NELEMY.
!
!    Output, integer ( kind = 4 ) ELEMENT_NUM, the number of elements in the grid.
!
  implicit none

  integer ( kind = 4 ) element_num
  integer ( kind = 4 ) nelemx
  integer ( kind = 4 ) nelemy

  element_num = 2 * nelemx * nelemy

  return
end
subroutine grid_t6_node_num ( nelemx, nelemy, node_num )

!*****************************************************************************80
!
!! GRID_T6_NODE_NUM counts the nodes in a grid of 6 node triangles.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NELEMX, NELEMY, the number of elements along the
!    X and Y directions.  The number of elements generated will be
!    NELEMX * NELEMY.
!
!    Output, integer ( kind = 4 ) NODE_NUM, the number of nodes in the grid.
!
  implicit none

  integer ( kind = 4 ) nelemx
  integer ( kind = 4 ) nelemy
  integer ( kind = 4 ) node_num

  node_num = ( 2 * nelemx + 1 ) * ( 2 * nelemy + 1 )

  return
end
subroutine grid_t10_element ( nelemx, nelemy, element_node )

!*****************************************************************************80
!
!! GRID_T10_ELEMENT produces an element grid of pairs of 10 node triangles.
!
!  Example:
!
!    Input:
!
!      NELEMX = 2, NELEMY = 2
!
!    Output:
!
!      ELEMENT_NODE = 
!         1,  2,  3,  4, 10, 16, 22, 15,  8,  9;
!        25, 24, 23, 22, 16, 10,  4, 11, 18, 17;
!         4,  5,  6,  7, 13, 19, 25, 18, 11, 12;
!        28, 27, 26, 25, 19, 13,  7, 14, 21, 20;
!        22, 23, 24, 25, 31, 37, 43, 36, 29, 30;
!        46, 45, 44, 43, 37, 31, 25, 32, 39, 38;
!        25, 26, 27, 28, 34, 40, 46, 39, 31, 33;
!        49, 48, 47, 46, 40, 34, 28, 35, 42, 41.
!        
!  Grid:
!
!   43-44-45-46-47-48-49
!    |\     6 |\     8 |
!    | \      | \      |
!   36 37 38 39 40 41 42
!    |   \    |   \    |
!    |    \   |    \   |
!   29 30 31 32 33 34 35
!    |      \ |      \ |
!    | 5     \| 7     \|
!   22-23-24-25-26-27-28
!    |\     2 |\     4 |
!    | \      | \      |
!   15 16 17 18 19 20 21
!    |   \    |   \    |
!    |    \   |    \   |
!    8  9 10 11 12 13 14
!    |      \ |      \ |
!    | 1     \| 3     \|
!    1--2--3--4--5--6--7
!
!  Reference Element T10:
!
!    |
!    1  10
!    |  |\
!    |  | \
!    |  8  9
!    |  |   \
!    S  |    \
!    |  5  6  7
!    |  |      \
!    |  |       \
!    0  1--2--3--4
!    |
!    +--0----R---1-->
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NELEMX, NELEMY, the number of elements along the
!    X and Y directions.  The number of elements generated will be
!    2 * NELEMX * NELEMY.
!
!    Output, integer ( kind = 4 ) ELEMENT_NODE(10,2*NELEMX*NELEMY), the nodes 
!    that form each element.
!
  implicit none

  integer ( kind = 4 ) nelemx
  integer ( kind = 4 ) nelemy
  integer ( kind = 4 ), parameter :: element_order = 10

  integer ( kind = 4 ) base
  integer ( kind = 4 ) element
  integer ( kind = 4 ) element_node(element_order,2*nelemx*nelemy)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  element = 0
 
  do j = 1, nelemy
    do i = 1, nelemx

      base = ( j - 1 ) * 3 * ( 3 * nelemx + 1 ) + 3 * i - 2

      element = element + 1

      element_node( 1,element) = base
      element_node( 2,element) = base                          + 1
      element_node( 3,element) = base                          + 2
      element_node( 4,element) = base                          + 3
      element_node( 5,element) = base +     ( 3 * nelemx + 1 ) + 2
      element_node( 6,element) = base + 2 * ( 3 * nelemx + 1 ) + 1
      element_node( 7,element) = base + 3 * ( 3 * nelemx + 1 )
      element_node( 8,element) = base + 2 * ( 3 * nelemx + 1 )
      element_node( 9,element) = base +     ( 2 * nelemx + 1 ) + 2
      element_node(10,element) = base +     ( 2 * nelemx + 1 ) + 3

      element = element + 1

      element_node( 1,element) = base + 3 * ( 3 * nelemx + 1 ) + 3
      element_node( 2,element) = base + 3 * ( 3 * nelemx + 1 ) + 2
      element_node( 3,element) = base + 3 * ( 3 * nelemx + 1 ) + 1
      element_node( 4,element) = base + 3 * ( 3 * nelemx + 1 )
      element_node( 5,element) = base + 2 * ( 3 * nelemx + 1 ) + 1
      element_node( 6,element) = base +     ( 3 * nelemx + 1 ) + 2
      element_node( 7,element) = base                          + 3
      element_node( 8,element) = base +     ( 3 * nelemx + 1 ) + 3
      element_node( 9,element) = base + 2 * ( 3 * nelemx + 1 ) + 3
      element_node(10,element) = base + 2 * ( 3 * nelemx + 1 ) + 2

    end do
  end do

  return
end
subroutine grid_t10_element_num ( nelemx, nelemy, element_num )

!*****************************************************************************80
!
!! GRID_T10_ELEMENT_NUM counts the elements in a grid of 10 node triangles.
!
!  Example:
!
!    Input:
!
!      NELEMX = 3, NELEMY = 2
!
!    Output:
!
!      ELEMENT_NUM = 2 * NELEMX * NELEMY = 12
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NELEMX, NELEMY, the number of elements along the
!    X and Y directions.  The number of elements generated will be
!    NELEMX * NELEMY.
!
!    Output, integer ( kind = 4 ) ELEMENT_NUM, the number of elements in the grid.
!
  implicit none

  integer ( kind = 4 ) element_num
  integer ( kind = 4 ) nelemx
  integer ( kind = 4 ) nelemy

  element_num = 2 * nelemx * nelemy

  return
end
subroutine grid_t10_node_num ( nelemx, nelemy, node_num )

!*****************************************************************************80
!
!! GRID_T10_NODE_NUM counts the nodes in a grid of 10 node triangles.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NELEMX, NELEMY, the number of elements along the
!    X and Y directions.  The number of elements generated will be
!    NELEMX * NELEMY.
!
!    Output, integer ( kind = 4 ) NODE_NUM, the number of nodes in the grid.
!
  implicit none

  integer ( kind = 4 ) nelemx
  integer ( kind = 4 ) nelemy
  integer ( kind = 4 ) node_num

  node_num = ( 3 * nelemx + 1 ) * ( 3 * nelemy + 1 )

  return
end
subroutine grid_test ( code )

!*****************************************************************************80
!
!! GRID_TEST tests the grid element routines.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) CODE, the code for the element.
!    Legal values include 'Q4', 'Q8', 'Q9', 'Q12', 'Q16', 'QL', 
!    'T3', 'T4', 'T6' and 'T10'.
!
  implicit none
 
  character ( len = * ) code
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: element_node
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ) element_order
  integer ( kind = 4 ) nelemx
  integer ( kind = 4 ) nelemy
  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) order_code
  integer ( kind = 4 ) width
!
!  NODE is defined as a vector rather than a two dimensional array,
!  so that we can handle the various cases using a single array.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  GRID_TEST for element "' // trim ( code ) // '".'

  nelemx = 3
  nelemy = 2

  write ( *, '(a,i8)' ) '  Number of elements in X direction = ', nelemx
  write ( *, '(a,i8)' ) '  Number of elements in Y direction = ', nelemy

  element_order = order_code ( code )

  call grid_node_num ( code, nelemx, nelemy, node_num )

  write ( *, '(a,i8)' ) '  Nodes per element =       ', element_order
  write ( *, '(a,i8)' ) '  Nodes in grid =           ', node_num

  call grid_element_num ( code, nelemx, nelemy, element_num )

  allocate ( element_node(element_order,element_num) )

  call grid_element ( code, element_order, nelemx, nelemy, element_node )

  call grid_print ( element_order, element_num, element_node )

  call grid_width ( element_order, element_num, element_node, width )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Grid width is ', width

  deallocate ( element_node )

  return
end
subroutine grid_width ( element_order, element_num, element_node, width )

!*****************************************************************************80
!
!! GRID_WIDTH computes the width of a given grid.
!
!  Discussion:
!
!    The grid width is defined to the maximum absolute
!    difference of global indices of nodes in the same element.
!
!  Example:
!
!    For the following grid, the grid width is 13.
!
!   23---24---25---26---27---28---29
!    |         |         |         |
!    |         |         |         |
!   19        20        21        22
!    |         |         |         |
!    | 4       | 5       | 6       |
!   12---13---14---15---16---17---18
!    |         |         |         |
!    |         |         |         |
!    8         9        10        11
!    |         |         |         |
!    | 1       | 2       | 3       |
!    1----2----3----4----5----6----7
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ELEMENT_ORDER, the number of nodes per element.
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of elements.
! 
!    Input, integer ( kind = 4 ) ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM), 
!    the nodes that form each element.
!
!    Output, integer ( kind = 4 ) WIDTH, the grid width.
!
  implicit none

  integer ( kind = 4 ) element_num
  integer ( kind = 4 ) element_order

  integer ( kind = 4 ) element
  integer ( kind = 4 ) element_node(element_order,element_num)
  integer ( kind = 4 ) ip1
  integer ( kind = 4 ) ip2
  integer ( kind = 4 ) node1
  integer ( kind = 4 ) node2
  integer ( kind = 4 ) width

  width = 0
 
  do element = 1, element_num
    do node1 = 1, element_order
      ip1 = element_node(node1,element)
      do node2 = 1, element_order
        ip2 = element_node(node2,element)
        width = max ( width, abs ( ip1 - ip2 ) )
      end do
    end do
  end do
 
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
    write ( *, '(a,i8)' ) '  I4_MODP ( I, J ) called with J = ', j
    stop
  end if

  i4_modp = mod ( i, j )

  if ( i4_modp < 0 ) then
    i4_modp = i4_modp + abs ( j )
  end if

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
!    19 August 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IVAL, an integer value.
!
!    Input, integer ( kind = 4 ) ILO, IHI, the desired bounds for the value.
!
!    Output, integer ( kind = 4 ) I4_WRAP, a "wrapped" version of IVAL.
!
  implicit none

  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) wide

  jlo = min ( ilo, ihi )
  jhi = max ( ilo, ihi )

  wide = jhi - jlo + 1

  if ( wide == 1 ) then
    i4_wrap = jlo
  else
    i4_wrap = jlo + i4_modp ( ival - jlo, wide )
  end if

  return
end
subroutine i4mat_write ( output_file_name, m, n, table )

!*****************************************************************************80
!
!! I4MAT_WRITE writes an I4MAT file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) OUTPUT_FILE_NAME, the output file name.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, integer ( kind = 4 ) TABLE(M,N), the table data.
!
  implicit none

  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) n

  integer   ( kind = 4 ) j
  character ( len = * )  output_file_name
  integer   ( kind = 4 ) output_status
  integer   ( kind = 4 ) output_unit
  character ( len = 30 ) string
  integer   ( kind = 4 ) table(m,n)
!
!  Open the file.
!
  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_file_name, &
    status = 'replace', iostat = output_status )

  if ( output_status /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4MAT_WRITE - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the output file "' // &
      trim ( output_file_name ) // '" on unit ', output_unit
    output_unit = -1
    stop
  end if
!
!  Create a format string.
!
  write ( string, '(a1,i8,a4)' ) '(', m, 'i10)'
!
!  Write the data.
!
  do j = 1, n
    write ( output_unit, string ) table(1:m,j)
  end do
!
!  Close the file.
!
  close ( unit = output_unit )

  return
end
subroutine i4vec_print ( n, a, title )

!*****************************************************************************80
!
!! I4VEC_PRINT prints an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, integer ( kind = 4 ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  implicit none

  integer ( kind = 4 ) n

  integer   ( kind = 4 ) a(n)
  integer   ( kind = 4 ) i
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i8,2x,i12)' ) i, a(i)
  end do

  return
end
subroutine interp ( code, element_order, r, s, node_u, u, dudr, duds )

!*****************************************************************************80
!
!! INTERP interpolates a quantity in an element from basis node values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) CODE, identifies the element.
!    Legal values include 'Q4', 'Q8', 'Q9', 'Q12', 'Q16', 'QL', 
!    'T3', 'T4', 'T6' and 'T10'.
!
!    Input, integer ( kind = 4 ) ELEMENT_ORDER, the number of nodes per element.
!
!    Input, real ( kind = 8 ) R, S, the reference coordinates of a point.
!
!    Input, real ( kind = 8 ) NODE_U(ELEMENT_ORDER), the value of the quantity 
!    at the basis nodes.
!
!    Output, real ( kind = 8 ) U, DUDR, DUDS, the interpolated value of the
!    quantity and its derivatives at the point (R,S).
!
  implicit none

  integer   ( kind = 4 ) element_order

  character ( len = * )  code
  real ( kind = 8 ) dtdr(element_order)
  real ( kind = 8 ) dtds(element_order)
  real ( kind = 8 ) dudr
  real ( kind = 8 ) duds
  real ( kind = 8 ) node_u(element_order)
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) t(element_order)
  real ( kind = 8 ) u

  call shape ( code, r, s, t, dtdr, dtds )
 
  u    = dot_product ( node_u(1:element_order), t(1:element_order) )
  dudr = dot_product ( node_u(1:element_order), dtdr(1:element_order) )
  duds = dot_product ( node_u(1:element_order), dtds(1:element_order) )
 
  return
end
subroutine interp_test ( code )

!*****************************************************************************80
!
!! INTERP_TEST tests the interpolation property of an element.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) CODE, identifies the element.
!    Legal values include 'Q4', 'Q8', 'Q9', 'Q12', 'Q16', 'QL', 
!    'T3', 'T4', 'T6' and 'T10'.
!
  implicit none

  real ( kind = 8 ) area
  character ( len = * )  code
  real ( kind = 8 ) dudr
  real ( kind = 8 ) dudr_exact
  real ( kind = 8 ) duds
  real ( kind = 8 ) duds_exact
  integer   ( kind = 4 ) element_order
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) node
  real ( kind = 8 ), allocatable, dimension ( : ) :: node_r
  real ( kind = 8 ), allocatable, dimension ( : ) :: node_s
  real ( kind = 8 ), allocatable, dimension ( : ) :: node_u
  integer   ( kind = 4 ) order_code
  real ( kind = 8 ) r
  real ( kind = 8 ) r_factor
  real ( kind = 8 ) r8_power
  integer   ( kind = 4 ), allocatable, dimension ( : ) :: rexp
  real ( kind = 8 ) s
  real ( kind = 8 ) s_factor
  integer   ( kind = 4 ), allocatable, dimension ( : ) :: sexp
  integer   ( kind = 4 ) seed
  integer   ( kind = 4 ) test
  integer   ( kind = 4 ), parameter :: test_num = 5
  real ( kind = 8 ) u
  real ( kind = 8 ) u_exact

  if ( code == 't4' .or. code == 'T4' ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'INTERP_TEST - Warning!'
    write ( *, '(a)' ) '  Skipping test for element "T4".'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'INTERP_TEST for element "' // trim ( code ) // '".'

  element_order = order_code ( code )

  write ( *, '(a,i8)' ) '  Element order = ', element_order

  allocate ( node_r(element_order) )
  allocate ( node_s(element_order) )
  allocate ( node_u(element_order) )
  allocate ( rexp(element_order) )
  allocate ( sexp(element_order) )
!
!  Get the coordinates of the reference nodes.
!
  call node_reference ( code, node_r, node_s, area )
!
!  Get the monomial exponents for which the element is exact.
!
  call poly ( code, rexp, sexp )

  seed = 123456789

  do i = 1, element_order

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8,a,i8)' ) '  Interpolate R ^ ', rexp(i), ' * S ^ ', sexp(i)
    write ( *, '(a)' ) ' '
!
!  Evaluate R**REXP(I) * S**SEXP(I) at the nodes.  This is our data.
!
    do node = 1, element_order
      r = node_r(node)
      s = node_s(node)
      if ( rexp(i) == 0 ) then
        r_factor = 1.0D+00
      else
        r_factor = r**rexp(i)
      end if
      if ( sexp(i) == 0 ) then
        s_factor = 1.0D+00
      else
        s_factor = s**sexp(i)
      end if
      node_u(node) = r_factor * s_factor
      write ( *, '(a,3g14.6)' ) '  (R,S,U):        ', r, s, node_u(node)
    end do
!
!  Now pick random points in the element, and compute the interpolated
!  value of R**REXP(*) * S**SEXP(I) there.  Mathematically, these
!  values should be exact.
!
    do test = 1, test_num

      call reference_sample ( code, seed, r, s )

      write ( *, '(a)' ) ' '
      write ( *, '(a,2g14.6)' ) '  (R,S):          ', r, s

      u_exact = r8_power ( r, rexp(i) ) * r8_power ( s, sexp(i) )

      dudr_exact = real ( rexp(i), kind = 8 ) &
        * r8_power ( r, rexp(i) - 1 ) * r8_power ( s, sexp(i) )

      duds_exact = r8_power ( r, rexp(i) ) * real ( sexp(i), kind = 8 ) &
        * r8_power ( s, sexp(i) - 1 )

      call interp ( code, element_order, r, s, node_u, u, dudr, duds )

      write ( *, '(a,3g14.6)' ) '  (U,U*,Error):   ', u_exact, u, &
        abs ( u_exact - u )
      write ( *, '(a,3g14.6)' ) '  (Ur,Ur*,Error): ', dudr_exact, dudr, &
        abs ( dudr_exact - dudr )
      write ( *, '(a,3g14.6)' ) '  (Us,Us*,Error): ', duds_exact, duds, &
        abs ( duds_exact - duds )

    end do

  end do

  deallocate ( node_r )
  deallocate ( node_s )
  deallocate ( node_u )
  deallocate ( rexp )
  deallocate ( sexp ) 

  return
end
subroutine legendre_com ( norder, xtab, weight )

!*****************************************************************************80
!
!! LEGENDRE_COM computes abscissas and weights for Gauss-Legendre quadrature.
!
!  Integration interval:
!
!    [ -1, 1 ]
!
!  Weight function:
!
!    1.
!
!  Integral to approximate:
!
!    Integral ( -1 <= X <= 1 ) F(X) dX.
!
!  Approximate integral:
!
!    sum ( 1 <= I <= NORDER ) WEIGHT(I) * F ( XTAB(I) ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 November 1998
!
!  Author:
!
!    Original FORTRAN77 version by Philip Davis, Philip Rabinowitz
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NORDER, the order of the rule.
!    NORDER must be greater than 0.
!
!    Output, real ( kind = 8 ) XTAB(NORDER), the abscissas of the rule.
!
!    Output, real ( kind = 8 ) WEIGHT(NORDER), the weights of the rule.
!    The weights are positive, symmetric, and should sum to 2.
!
  implicit none

  integer ( kind = 4 ) norder

  real    ( kind = 8 ) d1
  real    ( kind = 8 ) d2pn
  real    ( kind = 8 ) d3pn
  real    ( kind = 8 ) d4pn
  real    ( kind = 8 ) dp
  real    ( kind = 8 ) dpn
  real    ( kind = 8 ) e1
  real    ( kind = 8 ) fx
  real    ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iback
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mp1mi
  integer ( kind = 4 ) ncopy
  integer ( kind = 4 ) nmove
  real    ( kind = 8 ) p
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00 
  real    ( kind = 8 ) pk
  real    ( kind = 8 ) pkm1
  real    ( kind = 8 ) pkp1
  real    ( kind = 8 ) t
  real    ( kind = 8 ) u
  real    ( kind = 8 ) v
  real    ( kind = 8 ) x0
  real    ( kind = 8 ) xtab(norder)
  real    ( kind = 8 ) xtemp
  real    ( kind = 8 ) weight(norder)

  if ( norder < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LEGENDRE_COM - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of NORDER = ', norder
    stop
  end if
 
  e1 = real ( norder * ( norder + 1 ), kind = 8 )
 
  m = ( norder + 1 ) / 2
 
  do i = 1, ( norder + 1 ) / 2
 
    mp1mi = m + 1 - i
    t = pi * real ( 4 * i - 1,      kind = 8 ) &
           / real ( 4 * norder + 2, kind = 8 )
    x0 = cos(t) * ( 1.0D+00 - ( 1.0D+00 - 1.0D+00 / &
      real    ( norder, kind = 8 ) ) / real ( 8 * norder * norder, kind = 8 ) )
 
    pkm1 = 1.0D+00
    pk = x0

    do k = 2, norder
      pkp1 = 2.0D+00 * x0 * pk - pkm1 - ( x0 * pk - pkm1 ) &
        / real ( k, kind = 8 )
      pkm1 = pk
      pk = pkp1
    end do
 
    d1 = real ( norder, kind = 8 ) * ( pkm1 - x0 * pk )

    dpn = d1 / ( 1.0D+00 - x0 * x0 )

    d2pn = ( 2.0D+00 * x0 * dpn - e1 * pk ) / ( 1.0D+00 - x0 * x0 )

    d3pn = ( 4.0D+00 * x0 * d2pn + ( 2.0D+00 - e1 ) * dpn ) / &
      ( 1.0D+00 - x0 * x0 )

    d4pn = ( 6.0D+00 * x0 * d3pn + ( 6.0D+00 - e1 ) * d2pn ) &
      / ( 1.0D+00 - x0 * x0 )

    u = pk / dpn
    v = d2pn / dpn
!
!  Initial approximation H:
!
    h = - u * ( 1.0D+00 + 0.5D+00 * u * ( v + u * ( v * v - d3pn &
      / ( 3.0D+00 * dpn ) ) ) )
!
!  Refine H using one step of Newton's method:
!
    p = pk + h * ( dpn + 0.5D+00 * h * ( d2pn + h / 3.0D+00 &
      * ( d3pn + 0.25D+00 * h * d4pn ) ) )

    dp = dpn + h * ( d2pn + 0.5D+00 * h * ( d3pn + h * d4pn / 3.0D+00 ) )

    h = h - p / dp
 
    xtemp = x0 + h

    xtab(mp1mi) = xtemp
 
    fx = d1 - h * e1 * ( pk + 0.5D+00 * h * ( dpn + h / 3.0D+00 &
      * ( d2pn + 0.25D+00 * h * ( d3pn + 0.2D+00 * h * d4pn ) ) ) )

    weight(mp1mi) = 2.0D+00 * ( 1.0D+00 - xtemp * xtemp ) / ( fx * fx )
 
  end do
 
  if ( mod ( norder, 2 ) == 1 ) then
    xtab(1) = 0.0D+00
  end if
!
!  Shift the data up.
!
  nmove = int ( ( norder + 1 ) / 2 )
  ncopy = norder - nmove

  do i = 1, nmove
    iback = norder + 1 - i
    xtab(iback) = xtab(iback-ncopy)
    weight(iback) = weight(iback-ncopy)
  end do
!
!  Reflect values for the negative abscissas.
!
  do i = 1, norder - nmove
    xtab(i) = - xtab(norder+1-i)
    weight(i) = weight(norder+1-i)
  end do
 
  return
end
subroutine legendre_set ( n, x, w )

!*****************************************************************************80
!
!! LEGENDRE_SET sets abscissas and weights for Gauss-Legendre quadrature.
!
!  Discussion:
!
!    The integration interval is [ -1, 1 ].
!
!    The weight function is w(x) = 1.0.
!
!    The integral to approximate:
!
!      Integral ( -1 <= X <= 1 ) F(X) dX
!
!    The quadrature rule:
!
!      Sum ( 1 <= I <= N ) W(I) * F ( X(I) )
!
!    The quadrature rule will integrate exactly all polynomials up to
!    X**(2*N-1).
!
!    The abscissas of the rule are the zeroes of the Legendre polynomial
!    P(N)(X).
!
!    The integral produced by a Gauss-Legendre rule is equal to the
!    integral of the unique polynomial of degree N-1 which
!    agrees with the function at the ORDER abscissas of the rule.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Vladimir Krylov,
!    Approximate Calculation of Integrals,
!    Dover, 2006,
!    ISBN: 0486445798,
!    LC: QA311.K713.
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!    Daniel Zwillinger, editor,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996,
!    ISBN: 0-8493-2479-3,
!    LC: QA47.M315.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the rule.
!    N must be between 1 and 33 or 63, 64, 65, 127 or 255.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!    Output, real ( kind = 8 ) WN), the weights.
!    The weights are positive, symmetric and should sum to 2.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) w(n)
  real    ( kind = 8 ) x(n)

  if ( n == 1 ) then

    x(1) =   0.0D+00

    w(1) = 2.0D+00

  else if ( n == 2 ) then

    x(1) = -0.577350269189625764509148780502D+00
    x(2) =  0.577350269189625764509148780502D+00

    w(1) = 1.0D+00
    w(2) = 1.0D+00

  else if ( n == 3 ) then

    x(1) = -0.774596669241483377035853079956D+00
    x(2) =  0.000000000000000000000000000000D+00
    x(3) =  0.774596669241483377035853079956D+00

    w(1) = 5.0D+00 / 9.0D+00
    w(2) = 8.0D+00 / 9.0D+00
    w(3) = 5.0D+00 / 9.0D+00

  else if ( n == 4 ) then

    x(1) = -0.861136311594052575223946488893D+00
    x(2) = -0.339981043584856264802665759103D+00
    x(3) =  0.339981043584856264802665759103D+00
    x(4) =  0.861136311594052575223946488893D+00

    w(1) = 0.347854845137453857373063949222D+00
    w(2) = 0.652145154862546142626936050778D+00
    w(3) = 0.652145154862546142626936050778D+00
    w(4) = 0.347854845137453857373063949222D+00

  else if ( n == 5 ) then

    x(1) = -0.906179845938663992797626878299D+00
    x(2) = -0.538469310105683091036314420700D+00
    x(3) =  0.000000000000000000000000000000D+00
    x(4) =  0.538469310105683091036314420700D+00
    x(5) =  0.906179845938663992797626878299D+00

    w(1) = 0.236926885056189087514264040720D+00
    w(2) = 0.478628670499366468041291514836D+00
    w(3) = 0.568888888888888888888888888889D+00
    w(4) = 0.478628670499366468041291514836D+00
    w(5) = 0.236926885056189087514264040720D+00

  else if ( n == 6 ) then

    x(1) = - 0.932469514203152027812301554494D+00
    x(2) = - 0.661209386466264513661399595020D+00
    x(3) = - 0.238619186083196908630501721681D+00
    x(4) =   0.238619186083196908630501721681D+00
    x(5) =   0.661209386466264513661399595020D+00
    x(6) =   0.932469514203152027812301554494D+00

    w(1) = 0.171324492379170345040296142173D+00
    w(2) = 0.360761573048138607569833513838D+00
    w(3) = 0.467913934572691047389870343990D+00
    w(4) = 0.467913934572691047389870343990D+00
    w(5) = 0.360761573048138607569833513838D+00
    w(6) = 0.171324492379170345040296142173D+00

  else if ( n == 7 ) then

    x(1) = - 0.949107912342758524526189684048D+00
    x(2) = - 0.741531185599394439863864773281D+00
    x(3) = - 0.405845151377397166906606412077D+00
    x(4) =   0.0D+00
    x(5) =   0.405845151377397166906606412077D+00
    x(6) =   0.741531185599394439863864773281D+00
    x(7) =   0.949107912342758524526189684048D+00

    w(1) = 0.129484966168869693270611432679D+00
    w(2) = 0.279705391489276667901467771424D+00
    w(3) = 0.381830050505118944950369775489D+00
    w(4) = 0.417959183673469387755102040816D+00
    w(5) = 0.381830050505118944950369775489D+00
    w(6) = 0.279705391489276667901467771424D+00
    w(7) = 0.129484966168869693270611432679D+00

  else if ( n == 8 ) then

    x(1) = - 0.960289856497536231683560868569D+00
    x(2) = - 0.796666477413626739591553936476D+00
    x(3) = - 0.525532409916328985817739049189D+00
    x(4) = - 0.183434642495649804939476142360D+00
    x(5) =   0.183434642495649804939476142360D+00
    x(6) =   0.525532409916328985817739049189D+00
    x(7) =   0.796666477413626739591553936476D+00
    x(8) =   0.960289856497536231683560868569D+00

    w(1) = 0.101228536290376259152531354310D+00
    w(2) = 0.222381034453374470544355994426D+00
    w(3) = 0.313706645877887287337962201987D+00
    w(4) = 0.362683783378361982965150449277D+00
    w(5) = 0.362683783378361982965150449277D+00
    w(6) = 0.313706645877887287337962201987D+00
    w(7) = 0.222381034453374470544355994426D+00
    w(8) = 0.101228536290376259152531354310D+00

  else if ( n == 9 ) then

    x(1) = - 0.968160239507626089835576202904D+00
    x(2) = - 0.836031107326635794299429788070D+00
    x(3) = - 0.613371432700590397308702039341D+00
    x(4) = - 0.324253423403808929038538014643D+00
    x(5) =   0.0D+00
    x(6) =   0.324253423403808929038538014643D+00
    x(7) =   0.613371432700590397308702039341D+00
    x(8) =   0.836031107326635794299429788070D+00
    x(9) =   0.968160239507626089835576202904D+00

    w(1) = 0.812743883615744119718921581105D-01
    w(2) = 0.180648160694857404058472031243D+00
    w(3) = 0.260610696402935462318742869419D+00
    w(4) = 0.312347077040002840068630406584D+00
    w(5) = 0.330239355001259763164525069287D+00
    w(6) = 0.312347077040002840068630406584D+00
    w(7) = 0.260610696402935462318742869419D+00
    w(8) = 0.180648160694857404058472031243D+00
    w(9) = 0.812743883615744119718921581105D-01

  else if ( n == 10 ) then

    x(1) =  - 0.973906528517171720077964012084D+00
    x(2) =  - 0.865063366688984510732096688423D+00
    x(3) =  - 0.679409568299024406234327365115D+00
    x(4) =  - 0.433395394129247190799265943166D+00
    x(5) =  - 0.148874338981631210884826001130D+00
    x(6) =    0.148874338981631210884826001130D+00
    x(7) =    0.433395394129247190799265943166D+00
    x(8) =    0.679409568299024406234327365115D+00
    x(9) =    0.865063366688984510732096688423D+00
    x(10) =   0.973906528517171720077964012084D+00

    w(1) =  0.666713443086881375935688098933D-01
    w(2) =  0.149451349150580593145776339658D+00
    w(3) =  0.219086362515982043995534934228D+00
    w(4) =  0.269266719309996355091226921569D+00
    w(5) =  0.295524224714752870173892994651D+00
    w(6) =  0.295524224714752870173892994651D+00
    w(7) =  0.269266719309996355091226921569D+00
    w(8) =  0.219086362515982043995534934228D+00
    w(9) =  0.149451349150580593145776339658D+00
    w(10) = 0.666713443086881375935688098933D-01

  else if ( n == 11 ) then

    x(1) =  - 0.978228658146056992803938001123D+00
    x(2) =  - 0.887062599768095299075157769304D+00
    x(3) =  - 0.730152005574049324093416252031D+00
    x(4) =  - 0.519096129206811815925725669459D+00
    x(5) =  - 0.269543155952344972331531985401D+00
    x(6) =    0.0D+00
    x(7) =    0.269543155952344972331531985401D+00
    x(8) =    0.519096129206811815925725669459D+00
    x(9) =    0.730152005574049324093416252031D+00
    x(10) =   0.887062599768095299075157769304D+00
    x(11) =   0.978228658146056992803938001123D+00

    w(1) =  0.556685671161736664827537204425D-01
    w(2) =  0.125580369464904624634694299224D+00
    w(3) =  0.186290210927734251426097641432D+00
    w(4) =  0.233193764591990479918523704843D+00
    w(5) =  0.262804544510246662180688869891D+00
    w(6) =  0.272925086777900630714483528336D+00
    w(7) =  0.262804544510246662180688869891D+00
    w(8) =  0.233193764591990479918523704843D+00
    w(9) =  0.186290210927734251426097641432D+00
    w(10) = 0.125580369464904624634694299224D+00
    w(11) = 0.556685671161736664827537204425D-01

  else if ( n == 12 ) then

    x(1) =  - 0.981560634246719250690549090149D+00
    x(2) =  - 0.904117256370474856678465866119D+00
    x(3) =  - 0.769902674194304687036893833213D+00
    x(4) =  - 0.587317954286617447296702418941D+00
    x(5) =  - 0.367831498998180193752691536644D+00
    x(6) =  - 0.125233408511468915472441369464D+00
    x(7) =    0.125233408511468915472441369464D+00
    x(8) =    0.367831498998180193752691536644D+00
    x(9) =    0.587317954286617447296702418941D+00
    x(10) =   0.769902674194304687036893833213D+00
    x(11) =   0.904117256370474856678465866119D+00
    x(12) =   0.981560634246719250690549090149D+00

    w(1) =  0.471753363865118271946159614850D-01
    w(2) =  0.106939325995318430960254718194D+00
    w(3) =  0.160078328543346226334652529543D+00
    w(4) =  0.203167426723065921749064455810D+00
    w(5) =  0.233492536538354808760849898925D+00
    w(6) =  0.249147045813402785000562436043D+00
    w(7) =  0.249147045813402785000562436043D+00
    w(8) =  0.233492536538354808760849898925D+00
    w(9) =  0.203167426723065921749064455810D+00
    w(10) = 0.160078328543346226334652529543D+00
    w(11) = 0.106939325995318430960254718194D+00
    w(12) = 0.471753363865118271946159614850D-01

  else if ( n == 13 ) then

    x(1) =  - 0.984183054718588149472829448807D+00
    x(2) =  - 0.917598399222977965206547836501D+00
    x(3) =  - 0.801578090733309912794206489583D+00
    x(4) =  - 0.642349339440340220643984606996D+00
    x(5) =  - 0.448492751036446852877912852128D+00
    x(6) =  - 0.230458315955134794065528121098D+00
    x(7) =    0.0D+00
    x(8) =    0.230458315955134794065528121098D+00
    x(9) =    0.448492751036446852877912852128D+00
    x(10) =   0.642349339440340220643984606996D+00
    x(11) =   0.801578090733309912794206489583D+00
    x(12) =   0.917598399222977965206547836501D+00
    x(13) =   0.984183054718588149472829448807D+00

    w(1) =  0.404840047653158795200215922010D-01
    w(2) =  0.921214998377284479144217759538D-01
    w(3) =  0.138873510219787238463601776869D+00
    w(4) =  0.178145980761945738280046691996D+00
    w(5) =  0.207816047536888502312523219306D+00
    w(6) =  0.226283180262897238412090186040D+00
    w(7) =  0.232551553230873910194589515269D+00
    w(8) =  0.226283180262897238412090186040D+00
    w(9) =  0.207816047536888502312523219306D+00
    w(10) = 0.178145980761945738280046691996D+00
    w(11) = 0.138873510219787238463601776869D+00
    w(12) = 0.921214998377284479144217759538D-01
    w(13) = 0.404840047653158795200215922010D-01

  else if ( n == 14 ) then

    x(1) =  - 0.986283808696812338841597266704D+00
    x(2) =  - 0.928434883663573517336391139378D+00
    x(3) =  - 0.827201315069764993189794742650D+00
    x(4) =  - 0.687292904811685470148019803019D+00
    x(5) =  - 0.515248636358154091965290718551D+00
    x(6) =  - 0.319112368927889760435671824168D+00
    x(7) =  - 0.108054948707343662066244650220D+00
    x(8) =    0.108054948707343662066244650220D+00
    x(9) =    0.319112368927889760435671824168D+00
    x(10) =   0.515248636358154091965290718551D+00
    x(11) =   0.687292904811685470148019803019D+00
    x(12) =   0.827201315069764993189794742650D+00
    x(13) =   0.928434883663573517336391139378D+00
    x(14) =   0.986283808696812338841597266704D+00

    w(1) =  0.351194603317518630318328761382D-01
    w(2) =  0.801580871597602098056332770629D-01
    w(3) =  0.121518570687903184689414809072D+00
    w(4) =  0.157203167158193534569601938624D+00
    w(5) =  0.185538397477937813741716590125D+00
    w(6) =  0.205198463721295603965924065661D+00
    w(7) =  0.215263853463157790195876443316D+00
    w(8) =  0.215263853463157790195876443316D+00
    w(9) =  0.205198463721295603965924065661D+00
    w(10) = 0.185538397477937813741716590125D+00
    w(11) = 0.157203167158193534569601938624D+00
    w(12) = 0.121518570687903184689414809072D+00
    w(13) = 0.801580871597602098056332770629D-01
    w(14) = 0.351194603317518630318328761382D-01

  else if ( n == 15 ) then

    x(1) =  - 0.987992518020485428489565718587D+00
    x(2) =  - 0.937273392400705904307758947710D+00
    x(3) =  - 0.848206583410427216200648320774D+00
    x(4) =  - 0.724417731360170047416186054614D+00
    x(5) =  - 0.570972172608538847537226737254D+00
    x(6) =  - 0.394151347077563369897207370981D+00
    x(7) =  - 0.201194093997434522300628303395D+00
    x(8) =    0.0D+00
    x(9) =    0.201194093997434522300628303395D+00
    x(10) =   0.394151347077563369897207370981D+00
    x(11) =   0.570972172608538847537226737254D+00
    x(12) =   0.724417731360170047416186054614D+00
    x(13) =   0.848206583410427216200648320774D+00
    x(14) =   0.937273392400705904307758947710D+00
    x(15) =   0.987992518020485428489565718587D+00

    w(1) =  0.307532419961172683546283935772D-01
    w(2) =  0.703660474881081247092674164507D-01
    w(3) =  0.107159220467171935011869546686D+00
    w(4) =  0.139570677926154314447804794511D+00
    w(5) =  0.166269205816993933553200860481D+00
    w(6) =  0.186161000015562211026800561866D+00
    w(7) =  0.198431485327111576456118326444D+00
    w(8) =  0.202578241925561272880620199968D+00
    w(9) =  0.198431485327111576456118326444D+00
    w(10) = 0.186161000015562211026800561866D+00
    w(11) = 0.166269205816993933553200860481D+00
    w(12) = 0.139570677926154314447804794511D+00
    w(13) = 0.107159220467171935011869546686D+00
    w(14) = 0.703660474881081247092674164507D-01
    w(15) = 0.307532419961172683546283935772D-01

  else if ( n == 16 ) then

    x(1) =  - 0.989400934991649932596154173450D+00
    x(2) =  - 0.944575023073232576077988415535D+00
    x(3) =  - 0.865631202387831743880467897712D+00
    x(4) =  - 0.755404408355003033895101194847D+00
    x(5) =  - 0.617876244402643748446671764049D+00
    x(6) =  - 0.458016777657227386342419442984D+00
    x(7) =  - 0.281603550779258913230460501460D+00
    x(8) =  - 0.950125098376374401853193354250D-01
    x(9) =    0.950125098376374401853193354250D-01
    x(10) =   0.281603550779258913230460501460D+00
    x(11) =   0.458016777657227386342419442984D+00
    x(12) =   0.617876244402643748446671764049D+00
    x(13) =   0.755404408355003033895101194847D+00
    x(14) =   0.865631202387831743880467897712D+00
    x(15) =   0.944575023073232576077988415535D+00
    x(16) =   0.989400934991649932596154173450D+00

    w(1) =  0.271524594117540948517805724560D-01
    w(2) =  0.622535239386478928628438369944D-01
    w(3) =  0.951585116824927848099251076022D-01
    w(4) =  0.124628971255533872052476282192D+00
    w(5) =  0.149595988816576732081501730547D+00
    w(6) =  0.169156519395002538189312079030D+00
    w(7) =  0.182603415044923588866763667969D+00
    w(8) =  0.189450610455068496285396723208D+00
    w(9) =  0.189450610455068496285396723208D+00
    w(10) = 0.182603415044923588866763667969D+00
    w(11) = 0.169156519395002538189312079030D+00
    w(12) = 0.149595988816576732081501730547D+00
    w(13) = 0.124628971255533872052476282192D+00
    w(14) = 0.951585116824927848099251076022D-01
    w(15) = 0.622535239386478928628438369944D-01
    w(16) = 0.271524594117540948517805724560D-01

  else if ( n == 17 ) then

    x(1) =  - 0.990575475314417335675434019941D+00
    x(2) =  - 0.950675521768767761222716957896D+00
    x(3) =  - 0.880239153726985902122955694488D+00
    x(4) =  - 0.781514003896801406925230055520D+00
    x(5) =  - 0.657671159216690765850302216643D+00
    x(6) =  - 0.512690537086476967886246568630D+00
    x(7) =  - 0.351231763453876315297185517095D+00
    x(8) =  - 0.178484181495847855850677493654D+00
    x(9) =    0.0D+00
    x(10) =   0.178484181495847855850677493654D+00
    x(11) =   0.351231763453876315297185517095D+00
    x(12) =   0.512690537086476967886246568630D+00
    x(13) =   0.657671159216690765850302216643D+00
    x(14) =   0.781514003896801406925230055520D+00
    x(15) =   0.880239153726985902122955694488D+00
    x(16) =   0.950675521768767761222716957896D+00
    x(17) =   0.990575475314417335675434019941D+00

    w(1) =  0.241483028685479319601100262876D-01
    w(2) =  0.554595293739872011294401653582D-01
    w(3) =  0.850361483171791808835353701911D-01
    w(4) =  0.111883847193403971094788385626D+00
    w(5) =  0.135136368468525473286319981702D+00
    w(6) =  0.154045761076810288081431594802D+00
    w(7) =  0.168004102156450044509970663788D+00
    w(8) =  0.176562705366992646325270990113D+00
    w(9) =  0.179446470356206525458265644262D+00
    w(10) = 0.176562705366992646325270990113D+00
    w(11) = 0.168004102156450044509970663788D+00
    w(12) = 0.154045761076810288081431594802D+00
    w(13) = 0.135136368468525473286319981702D+00
    w(14) = 0.111883847193403971094788385626D+00
    w(15) = 0.850361483171791808835353701911D-01
    w(16) = 0.554595293739872011294401653582D-01
    w(17) = 0.241483028685479319601100262876D-01

  else if ( n == 18 ) then

    x(1) =  - 0.991565168420930946730016004706D+00
    x(2) =  - 0.955823949571397755181195892930D+00
    x(3) =  - 0.892602466497555739206060591127D+00
    x(4) =  - 0.803704958972523115682417455015D+00
    x(5) =  - 0.691687043060353207874891081289D+00
    x(6) =  - 0.559770831073947534607871548525D+00
    x(7) =  - 0.411751161462842646035931793833D+00
    x(8) =  - 0.251886225691505509588972854878D+00
    x(9) =  - 0.847750130417353012422618529358D-01
    x(10) =   0.847750130417353012422618529358D-01
    x(11) =   0.251886225691505509588972854878D+00
    x(12) =   0.411751161462842646035931793833D+00
    x(13) =   0.559770831073947534607871548525D+00
    x(14) =   0.691687043060353207874891081289D+00
    x(15) =   0.803704958972523115682417455015D+00
    x(16) =   0.892602466497555739206060591127D+00
    x(17) =   0.955823949571397755181195892930D+00
    x(18) =   0.991565168420930946730016004706D+00

    w(1) =  0.216160135264833103133427102665D-01
    w(2) =  0.497145488949697964533349462026D-01
    w(3) =  0.764257302548890565291296776166D-01
    w(4) =  0.100942044106287165562813984925D+00
    w(5) =  0.122555206711478460184519126800D+00
    w(6) =  0.140642914670650651204731303752D+00
    w(7) =  0.154684675126265244925418003836D+00
    w(8) =  0.164276483745832722986053776466D+00
    w(9) =  0.169142382963143591840656470135D+00
    w(10) = 0.169142382963143591840656470135D+00
    w(11) = 0.164276483745832722986053776466D+00
    w(12) = 0.154684675126265244925418003836D+00
    w(13) = 0.140642914670650651204731303752D+00
    w(14) = 0.122555206711478460184519126800D+00
    w(15) = 0.100942044106287165562813984925D+00
    w(16) = 0.764257302548890565291296776166D-01
    w(17) = 0.497145488949697964533349462026D-01
    w(18) = 0.216160135264833103133427102665D-01

  else if ( n == 19 ) then

    x(1) =  - 0.992406843843584403189017670253D+00
    x(2) =  - 0.960208152134830030852778840688D+00
    x(3) =  - 0.903155903614817901642660928532D+00
    x(4) =  - 0.822714656537142824978922486713D+00
    x(5) =  - 0.720966177335229378617095860824D+00
    x(6) =  - 0.600545304661681023469638164946D+00
    x(7) =  - 0.464570741375960945717267148104D+00
    x(8) =  - 0.316564099963629831990117328850D+00
    x(9) =  - 0.160358645640225375868096115741D+00
    x(10) =   0.0D+00
    x(11) =   0.160358645640225375868096115741D+00
    x(12) =   0.316564099963629831990117328850D+00
    x(13) =   0.464570741375960945717267148104D+00
    x(14) =   0.600545304661681023469638164946D+00
    x(15) =   0.720966177335229378617095860824D+00
    x(16) =   0.822714656537142824978922486713D+00
    x(17) =   0.903155903614817901642660928532D+00
    x(18) =   0.960208152134830030852778840688D+00
    x(19) =   0.992406843843584403189017670253D+00

    w(1) =  0.194617882297264770363120414644D-01
    w(2) =  0.448142267656996003328381574020D-01
    w(3) =  0.690445427376412265807082580060D-01
    w(4) =  0.914900216224499994644620941238D-01
    w(5) =  0.111566645547333994716023901682D+00
    w(6) =  0.128753962539336227675515784857D+00
    w(7) =  0.142606702173606611775746109442D+00
    w(8) =  0.152766042065859666778855400898D+00
    w(9) =  0.158968843393954347649956439465D+00
    w(10) = 0.161054449848783695979163625321D+00
    w(11) = 0.158968843393954347649956439465D+00
    w(12) = 0.152766042065859666778855400898D+00
    w(13) = 0.142606702173606611775746109442D+00
    w(14) = 0.128753962539336227675515784857D+00
    w(15) = 0.111566645547333994716023901682D+00
    w(16) = 0.914900216224499994644620941238D-01
    w(17) = 0.690445427376412265807082580060D-01
    w(18) = 0.448142267656996003328381574020D-01
    w(19) = 0.194617882297264770363120414644D-01

  else if ( n == 20 ) then

    x(1) =  - 0.993128599185094924786122388471D+00
    x(2) =  - 0.963971927277913791267666131197D+00
    x(3) =  - 0.912234428251325905867752441203D+00
    x(4) =  - 0.839116971822218823394529061702D+00
    x(5) =  - 0.746331906460150792614305070356D+00
    x(6) =  - 0.636053680726515025452836696226D+00
    x(7) =  - 0.510867001950827098004364050955D+00
    x(8) =  - 0.373706088715419560672548177025D+00
    x(9) =  - 0.227785851141645078080496195369D+00
    x(10) = - 0.765265211334973337546404093988D-01
    x(11) =   0.765265211334973337546404093988D-01
    x(12) =   0.227785851141645078080496195369D+00
    x(13) =   0.373706088715419560672548177025D+00
    x(14) =   0.510867001950827098004364050955D+00
    x(15) =   0.636053680726515025452836696226D+00
    x(16) =   0.746331906460150792614305070356D+00
    x(17) =   0.839116971822218823394529061702D+00
    x(18) =   0.912234428251325905867752441203D+00
    x(19) =   0.963971927277913791267666131197D+00
    x(20) =   0.993128599185094924786122388471D+00

    w(1) =  0.176140071391521183118619623519D-01
    w(2) =  0.406014298003869413310399522749D-01
    w(3) =  0.626720483341090635695065351870D-01
    w(4) =  0.832767415767047487247581432220D-01
    w(5) =  0.101930119817240435036750135480D+00
    w(6) =  0.118194531961518417312377377711D+00
    w(7) =  0.131688638449176626898494499748D+00
    w(8) =  0.142096109318382051329298325067D+00
    w(9) =  0.149172986472603746787828737002D+00
    w(10) = 0.152753387130725850698084331955D+00
    w(11) = 0.152753387130725850698084331955D+00
    w(12) = 0.149172986472603746787828737002D+00
    w(13) = 0.142096109318382051329298325067D+00
    w(14) = 0.131688638449176626898494499748D+00
    w(15) = 0.118194531961518417312377377711D+00
    w(16) = 0.101930119817240435036750135480D+00
    w(17) = 0.832767415767047487247581432220D-01
    w(18) = 0.626720483341090635695065351870D-01
    w(19) = 0.406014298003869413310399522749D-01
    w(20) = 0.176140071391521183118619623519D-01

  else if ( n == 21 ) then

    x( 1) =  -0.9937521706203896D+00
    x( 2) =  -0.9672268385663063D+00
    x( 3) =  -0.9200993341504008D+00
    x( 4) =  -0.8533633645833173D+00
    x( 5) =  -0.7684399634756779D+00
    x( 6) =  -0.6671388041974123D+00
    x( 7) =  -0.5516188358872198D+00
    x( 8) =  -0.4243421202074388D+00
    x( 9) =  -0.2880213168024011D+00
    x(10) =  -0.1455618541608951D+00
    x(11) =   0.0000000000000000D+00
    x(12) =   0.1455618541608951D+00
    x(13) =   0.2880213168024011D+00
    x(14) =   0.4243421202074388D+00
    x(15) =   0.5516188358872198D+00
    x(16) =   0.6671388041974123D+00
    x(17) =   0.7684399634756779D+00
    x(18) =   0.8533633645833173D+00
    x(19) =   0.9200993341504008D+00
    x(20) =   0.9672268385663063D+00
    x(21) =   0.9937521706203896D+00 
   
    w( 1) =   0.1601722825777420D-01
    w( 2) =   0.3695378977085242D-01
    w( 3) =   0.5713442542685715D-01
    w( 4) =   0.7610011362837928D-01
    w( 5) =   0.9344442345603393D-01
    w( 6) =   0.1087972991671484D+00
    w( 7) =   0.1218314160537285D+00
    w( 8) =   0.1322689386333373D+00
    w( 9) =   0.1398873947910731D+00
    w(10) =   0.1445244039899700D+00
    w(11) =   0.1460811336496904D+00
    w(12) =   0.1445244039899700D+00
    w(13) =   0.1398873947910731D+00
    w(14) =   0.1322689386333373D+00
    w(15) =   0.1218314160537285D+00
    w(16) =   0.1087972991671484D+00
    w(17) =   0.9344442345603393D-01
    w(18) =   0.7610011362837928D-01
    w(19) =   0.5713442542685715D-01
    w(20) =   0.3695378977085242D-01
    w(21) =   0.1601722825777420D-01

  else if ( n == 22 ) then

    x( 1) =  -0.9942945854823994D+00
    x( 2) =  -0.9700604978354287D+00
    x( 3) =  -0.9269567721871740D+00
    x( 4) =  -0.8658125777203002D+00
    x( 5) =  -0.7878168059792081D+00
    x( 6) =  -0.6944872631866827D+00
    x( 7) =  -0.5876404035069116D+00
    x( 8) =  -0.4693558379867570D+00
    x( 9) =  -0.3419358208920842D+00
    x(10) =  -0.2078604266882213D+00
    x(11) =  -0.6973927331972223D-01
    x(12) =   0.6973927331972223D-01
    x(13) =   0.2078604266882213D+00
    x(14) =   0.3419358208920842D+00
    x(15) =   0.4693558379867570D+00
    x(16) =   0.5876404035069116D+00
    x(17) =   0.6944872631866827D+00
    x(18) =   0.7878168059792081D+00
    x(19) =   0.8658125777203002D+00
    x(20) =   0.9269567721871740D+00
    x(21) =   0.9700604978354287D+00
    x(22) =   0.9942945854823994D+00
 
    w( 1) =   0.1462799529827203D-01
    w( 2) =   0.3377490158481413D-01
    w( 3) =   0.5229333515268327D-01
    w( 4) =   0.6979646842452038D-01
    w( 5) =   0.8594160621706777D-01
    w( 6) =   0.1004141444428809D+00
    w( 7) =   0.1129322960805392D+00
    w( 8) =   0.1232523768105124D+00
    w( 9) =   0.1311735047870623D+00
    w(10) =   0.1365414983460152D+00
    w(11) =   0.1392518728556321D+00
    w(12) =   0.1392518728556321D+00
    w(13) =   0.1365414983460152D+00
    w(14) =   0.1311735047870623D+00
    w(15) =   0.1232523768105124D+00
    w(16) =   0.1129322960805392D+00
    w(17) =   0.1004141444428809D+00
    w(18) =   0.8594160621706777D-01
    w(19) =   0.6979646842452038D-01
    w(20) =   0.5229333515268327D-01
    w(21) =   0.3377490158481413D-01
    w(22) =   0.1462799529827203D-01

  else if ( n == 23 ) then

    x( 1) =  -0.9947693349975522D+00
    x( 2) =  -0.9725424712181152D+00
    x( 3) =  -0.9329710868260161D+00
    x( 4) =  -0.8767523582704416D+00
    x( 5) =  -0.8048884016188399D+00
    x( 6) =  -0.7186613631319502D+00
    x( 7) =  -0.6196098757636461D+00
    x( 8) =  -0.5095014778460075D+00
    x( 9) =  -0.3903010380302908D+00
    x(10) =  -0.2641356809703449D+00
    x(11) =  -0.1332568242984661D+00
    x(12) =   0.0000000000000000D+00
    x(13) =   0.1332568242984661D+00
    x(14) =   0.2641356809703449D+00
    x(15) =   0.3903010380302908D+00
    x(16) =   0.5095014778460075D+00
    x(17) =   0.6196098757636461D+00
    x(18) =   0.7186613631319502D+00
    x(19) =   0.8048884016188399D+00
    x(20) =   0.8767523582704416D+00
    x(21) =   0.9329710868260161D+00
    x(22) =   0.9725424712181152D+00
    x(23) =   0.9947693349975522D+00
 
    w( 1) =   0.1341185948714167D-01
    w( 2) =   0.3098800585697944D-01
    w( 3) =   0.4803767173108464D-01
    w( 4) =   0.6423242140852586D-01
    w( 5) =   0.7928141177671895D-01
    w( 6) =   0.9291576606003514D-01
    w( 7) =   0.1048920914645414D+00
    w( 8) =   0.1149966402224114D+00
    w( 9) =   0.1230490843067295D+00
    w(10) =   0.1289057221880822D+00
    w(11) =   0.1324620394046967D+00
    w(12) =   0.1336545721861062D+00
    w(13) =   0.1324620394046967D+00
    w(14) =   0.1289057221880822D+00
    w(15) =   0.1230490843067295D+00
    w(16) =   0.1149966402224114D+00
    w(17) =   0.1048920914645414D+00
    w(18) =   0.9291576606003514D-01
    w(19) =   0.7928141177671895D-01
    w(20) =   0.6423242140852586D-01
    w(21) =   0.4803767173108464D-01
    w(22) =   0.3098800585697944D-01
    w(23) =   0.1341185948714167D-01

  else if ( n == 24 ) then

    x( 1) =  -0.9951872199970213D+00    
    x( 2) =  -0.9747285559713095D+00    
    x( 3) =  -0.9382745520027327D+00    
    x( 4) =  -0.8864155270044011D+00    
    x( 5) =  -0.8200019859739029D+00    
    x( 6) =  -0.7401241915785544D+00    
    x( 7) =  -0.6480936519369755D+00    
    x( 8) =  -0.5454214713888396D+00    
    x( 9) =  -0.4337935076260451D+00    
    x(10) =  -0.3150426796961634D+00    
    x(11) =  -0.1911188674736163D+00    
    x(12) =  -0.6405689286260562D-01
    x(13) =   0.6405689286260562D-01
    x(14) =   0.1911188674736163D+00    
    x(15) =   0.3150426796961634D+00    
    x(16) =   0.4337935076260451D+00    
    x(17) =   0.5454214713888396D+00    
    x(18) =   0.6480936519369755D+00    
    x(19) =   0.7401241915785544D+00    
    x(20) =   0.8200019859739029D+00    
    x(21) =   0.8864155270044011D+00    
    x(22) =   0.9382745520027327D+00    
    x(23) =   0.9747285559713095D+00    
    x(24) =   0.9951872199970213D+00    
 
    w( 1) =   0.1234122979998730D-01
    w( 2) =   0.2853138862893375D-01
    w( 3) =   0.4427743881741982D-01
    w( 4) =   0.5929858491543672D-01
    w( 5) =   0.7334648141108031D-01
    w( 6) =   0.8619016153195320D-01
    w( 7) =   0.9761865210411380D-01
    w( 8) =   0.1074442701159656D+00    
    w( 9) =   0.1155056680537256D+00    
    w(10) =   0.1216704729278035D+00    
    w(11) =   0.1258374563468283D+00    
    w(12) =   0.1279381953467521D+00    
    w(13) =   0.1279381953467521D+00    
    w(14) =   0.1258374563468283D+00    
    w(15) =   0.1216704729278035D+00    
    w(16) =   0.1155056680537256D+00   
    w(17) =   0.1074442701159656D+00    
    w(18) =   0.9761865210411380D-01
    w(19) =   0.8619016153195320D-01
    w(20) =   0.7334648141108031D-01
    w(21) =   0.5929858491543672D-01
    w(22) =   0.4427743881741982D-01
    w(23) =   0.2853138862893375D-01
    w(24) =   0.1234122979998730D-01

  else if ( n == 25 ) then

    x( 1) =  -0.9955569697904981D+00    
    x( 2) =  -0.9766639214595175D+00    
    x( 3) =  -0.9429745712289743D+00    
    x( 4) =  -0.8949919978782754D+00    
    x( 5) =  -0.8334426287608340D+00    
    x( 6) =  -0.7592592630373577D+00    
    x( 7) =  -0.6735663684734684D+00    
    x( 8) =  -0.5776629302412229D+00    
    x( 9) =  -0.4730027314457150D+00    
    x(10) =  -0.3611723058093879D+00    
    x(11) =  -0.2438668837209884D+00    
    x(12) =  -0.1228646926107104D+00    
    x(13) =   0.0000000000000000D+00    
    x(14) =   0.1228646926107104D+00  
    x(15) =   0.2438668837209884D+00    
    x(16) =   0.3611723058093879D+00    
    x(17) =   0.4730027314457150D+00    
    x(18) =   0.5776629302412229D+00    
    x(19) =   0.6735663684734684D+00    
    x(20) =   0.7592592630373577D+00    
    x(21) =   0.8334426287608340D+00    
    x(22) =   0.8949919978782754D+00    
    x(23) =   0.9429745712289743D+00    
    x(24) =   0.9766639214595175D+00    
    x(25) =   0.9955569697904981D+00    
 
    w( 1) =   0.1139379850102617D-01
    w( 2) =   0.2635498661503214D-01
    w( 3) =   0.4093915670130639D-01
    w( 4) =   0.5490469597583517D-01
    w( 5) =   0.6803833381235694D-01
    w( 6) =   0.8014070033500101D-01
    w( 7) =   0.9102826198296370D-01
    w( 8) =   0.1005359490670506D+00    
    w( 9) =   0.1085196244742637D+00    
    w(10) =   0.1148582591457116D+00    
    w(11) =   0.1194557635357847D+00    
    w(12) =   0.1222424429903101D+00    
    w(13) =   0.1231760537267154D+00    
    w(14) =   0.1222424429903101D+00    
    w(15) =   0.1194557635357847D+00    
    w(16) =   0.1148582591457116D+00    
    w(17) =   0.1085196244742637D+00    
    w(18) =   0.1005359490670506D+00    
    w(19) =   0.9102826198296370D-01
    w(20) =   0.8014070033500101D-01
    w(21) =   0.6803833381235694D-01
    w(22) =   0.5490469597583517D-01
    w(23) =   0.4093915670130639D-01
    w(24) =   0.2635498661503214D-01
    w(25) =   0.1139379850102617D-01

  else if ( n == 26 ) then

    x( 1) =  -0.9958857011456169D+00    
    x( 2) =  -0.9783854459564710D+00    
    x( 3) =  -0.9471590666617142D+00    
    x( 4) =  -0.9026378619843071D+00    
    x( 5) =  -0.8454459427884981D+00    
    x( 6) =  -0.7763859488206789D+00    
    x( 7) =  -0.6964272604199573D+00    
    x( 8) =  -0.6066922930176181D+00    
    x( 9) =  -0.5084407148245057D+00    
    x(10) =  -0.4030517551234863D+00    
    x(11) =  -0.2920048394859569D+00    
    x(12) =  -0.1768588203568902D+00    
    x(13) =  -0.5923009342931320D-01
    x(14) =   0.5923009342931320D-01
    x(15) =   0.1768588203568902D+00    
    x(16) =   0.2920048394859569D+00    
    x(17) =   0.4030517551234863D+00    
    x(18) =   0.5084407148245057D+00    
    x(19) =   0.6066922930176181D+00    
    x(20) =   0.6964272604199573D+00    
    x(21) =   0.7763859488206789D+00    
    x(22) =   0.8454459427884981D+00    
    x(23) =   0.9026378619843071D+00    
    x(24) =   0.9471590666617142D+00    
    x(25) =   0.9783854459564710D+00    
    x(26) =   0.9958857011456169D+00    
 
    w( 1) =   0.1055137261734304D-01
    w( 2) =   0.2441785109263173D-01
    w( 3) =   0.3796238329436282D-01
    w( 4) =   0.5097582529714782D-01
    w( 5) =   0.6327404632957484D-01
    w( 6) =   0.7468414976565967D-01
    w( 7) =   0.8504589431348521D-01
    w( 8) =   0.9421380035591416D-01
    w( 9) =   0.1020591610944255D+00    
    w(10) =   0.1084718405285765D+00    
    w(11) =   0.1133618165463197D+00    
    w(12) =   0.1166604434852967D+00    
    w(13) =   0.1183214152792622D+00    
    w(14) =   0.1183214152792622D+00    
    w(15) =   0.1166604434852967D+00    
    w(16) =   0.1133618165463197D+00    
    w(17) =   0.1084718405285765D+00    
    w(18) =   0.1020591610944255D+00    
    w(19) =   0.9421380035591416D-01
    w(20) =   0.8504589431348521D-01
    w(21) =   0.7468414976565967D-01
    w(22) =   0.6327404632957484D-01
    w(23) =   0.5097582529714782D-01
    w(24) =   0.3796238329436282D-01
    w(25) =   0.2441785109263173D-01
    w(26) =   0.1055137261734304D-01

  else if ( n == 27 ) then

    x( 1) =  -0.9961792628889886D+00    
    x( 2) =  -0.9799234759615012D+00    
    x( 3) =  -0.9509005578147051D+00    
    x( 4) =  -0.9094823206774911D+00    
    x( 5) =  -0.8562079080182945D+00    
    x( 6) =  -0.7917716390705082D+00    
    x( 7) =  -0.7170134737394237D+00    
    x( 8) =  -0.6329079719464952D+00    
    x( 9) =  -0.5405515645794569D+00    
    x(10) =  -0.4411482517500269D+00    
    x(11) =  -0.3359939036385089D+00    
    x(12) =  -0.2264593654395369D+00    
    x(13) =  -0.1139725856095300D+00    
    x(14) =   0.0000000000000000D+00    
    x(15) =   0.1139725856095300D+00    
    x(16) =   0.2264593654395369D+00    
    x(17) =   0.3359939036385089D+00    
    x(18) =   0.4411482517500269D+00    
    x(19) =   0.5405515645794569D+00    
    x(20) =   0.6329079719464952D+00    
    x(21) =   0.7170134737394237D+00    
    x(22) =   0.7917716390705082D+00    
    x(23) =   0.8562079080182945D+00    
    x(24) =   0.9094823206774911D+00    
    x(25) =   0.9509005578147051D+00    
    x(26) =   0.9799234759615012D+00    
    x(27) =   0.9961792628889886D+00    
 
    w( 1) =   0.9798996051294232D-02
    w( 2) =   0.2268623159618062D-01
    w( 3) =   0.3529705375741969D-01
    w( 4) =   0.4744941252061504D-01
    w( 5) =   0.5898353685983366D-01
    w( 6) =   0.6974882376624561D-01
    w( 7) =   0.7960486777305781D-01
    w( 8) =   0.8842315854375689D-01
    w( 9) =   0.9608872737002842D-01
    w(10) =   0.1025016378177459D+00    
    w(11) =   0.1075782857885332D+00    
    w(12) =   0.1112524883568452D+00    
    w(13) =   0.1134763461089651D+00    
    w(14) =   0.1142208673789570D+00    
    w(15) =   0.1134763461089651D+00    
    w(16) =   0.1112524883568452D+00    
    w(17) =   0.1075782857885332D+00    
    w(18) =   0.1025016378177459D+00    
    w(19) =   0.9608872737002842D-01
    w(20) =   0.8842315854375689D-01
    w(21) =   0.7960486777305781D-01
    w(22) =   0.6974882376624561D-01
    w(23) =   0.5898353685983366D-01
    w(24) =   0.4744941252061504D-01
    w(25) =   0.3529705375741969D-01
    w(26) =   0.2268623159618062D-01
    w(27) =   0.9798996051294232D-02

  else if ( n == 28 ) then

    x( 1) =  -0.9964424975739544D+00    
    x( 2) =  -0.9813031653708728D+00    
    x( 3) =  -0.9542592806289382D+00    
    x( 4) =  -0.9156330263921321D+00    
    x( 5) =  -0.8658925225743951D+00    
    x( 6) =  -0.8056413709171791D+00    
    x( 7) =  -0.7356108780136318D+00    
    x( 8) =  -0.6566510940388650D+00    
    x( 9) =  -0.5697204718114017D+00    
    x(10) =  -0.4758742249551183D+00    
    x(11) =  -0.3762515160890787D+00    
    x(12) =  -0.2720616276351780D+00    
    x(13) =  -0.1645692821333808D+00    
    x(14) =  -0.5507928988403427D-01
    x(15) =   0.5507928988403427D-01
    x(16) =   0.1645692821333808D+00    
    x(17) =   0.2720616276351780D+00    
    x(18) =   0.3762515160890787D+00    
    x(19) =   0.4758742249551183D+00    
    x(20) =   0.5697204718114017D+00    
    x(21) =   0.6566510940388650D+00    
    x(22) =   0.7356108780136318D+00    
    x(23) =   0.8056413709171791D+00    
    x(24) =   0.8658925225743951D+00    
    x(25) =   0.9156330263921321D+00    
    x(26) =   0.9542592806289382D+00    
    x(27) =   0.9813031653708728D+00    
    x(28) =   0.9964424975739544D+00    
 
    w( 1) =   0.9124282593094672D-02
    w( 2) =   0.2113211259277118D-01
    w( 3) =   0.3290142778230441D-01
    w( 4) =   0.4427293475900429D-01
    w( 5) =   0.5510734567571667D-01
    w( 6) =   0.6527292396699959D-01
    w( 7) =   0.7464621423456877D-01
    w( 8) =   0.8311341722890127D-01
    w( 9) =   0.9057174439303289D-01
    w(10) =   0.9693065799792999D-01
    w(11) =   0.1021129675780608D+00    
    w(12) =   0.1060557659228464D+00    
    w(13) =   0.1087111922582942D+00    
    w(14) =   0.1100470130164752D+00    
    w(15) =   0.1100470130164752D+00    
    w(16) =   0.1087111922582942D+00    
    w(17) =   0.1060557659228464D+00    
    w(18) =   0.1021129675780608D+00   
    w(19) =   0.9693065799792999D-01
    w(20) =   0.9057174439303289D-01
    w(21) =   0.8311341722890127D-01
    w(22) =   0.7464621423456877D-01
    w(23) =   0.6527292396699959D-01
    w(24) =   0.5510734567571667D-01
    w(25) =   0.4427293475900429D-01
    w(26) =   0.3290142778230441D-01
    w(27) =   0.2113211259277118D-01
    w(28) =   0.9124282593094672D-02

  else if ( n == 29 ) then

    x( 1) =  -0.9966794422605966D+00    
    x( 2) =  -0.9825455052614132D+00    
    x( 3) =  -0.9572855957780877D+00    
    x( 4) =  -0.9211802329530588D+00    
    x( 5) =  -0.8746378049201028D+00    
    x( 6) =  -0.8181854876152524D+00    
    x( 7) =  -0.7524628517344771D+00    
    x( 8) =  -0.6782145376026865D+00    
    x( 9) =  -0.5962817971382278D+00    
    x(10) =  -0.5075929551242276D+00    
    x(11) =  -0.4131528881740087D+00    
    x(12) =  -0.3140316378676399D+00    
    x(13) =  -0.2113522861660011D+00    
    x(14) =  -0.1062782301326792D+00    
    x(15) =   0.0000000000000000D+00    
    x(16) =   0.1062782301326792D+00    
    x(17) =   0.2113522861660011D+00    
    x(18) =   0.3140316378676399D+00    
    x(19) =   0.4131528881740087D+00    
    x(20) =   0.5075929551242276D+00    
    x(21) =   0.5962817971382278D+00    
    x(22) =   0.6782145376026865D+00    
    x(23) =   0.7524628517344771D+00    
    x(24) =   0.8181854876152524D+00    
    x(25) =   0.8746378049201028D+00    
    x(26) =   0.9211802329530588D+00    
    x(27) =   0.9572855957780877D+00    
    x(28) =   0.9825455052614132D+00    
    x(29) =   0.9966794422605966D+00    
 
    w( 1) =   0.8516903878746365D-02
    w( 2) =   0.1973208505612276D-01
    w( 3) =   0.3074049220209360D-01
    w( 4) =   0.4140206251868281D-01
    w( 5) =   0.5159482690249799D-01
    w( 6) =   0.6120309065707916D-01
    w( 7) =   0.7011793325505125D-01
    w( 8) =   0.7823832713576385D-01
    w( 9) =   0.8547225736617248D-01
    w(10) =   0.9173775713925882D-01
    w(11) =   0.9696383409440862D-01
    w(12) =   0.1010912737599150D+00    
    w(13) =   0.1040733100777293D+00    
    w(14) =   0.1058761550973210D+00    
    w(15) =   0.1064793817183143D+00    
    w(16) =   0.1058761550973210D+00    
    w(17) =   0.1040733100777293D+00    
    w(18) =   0.1010912737599150D+00    
    w(19) =   0.9696383409440862D-01
    w(20) =   0.9173775713925882D-01
    w(21) =   0.8547225736617248D-01
    w(22) =   0.7823832713576385D-01
    w(23) =   0.7011793325505125D-01
    w(24) =   0.6120309065707916D-01
    w(25) =   0.5159482690249799D-01
    w(26) =   0.4140206251868281D-01
    w(27) =   0.3074049220209360D-01
    w(28) =   0.1973208505612276D-01
    w(29) =   0.8516903878746365D-02

  else if ( n == 30 ) then

    x( 1) =  -0.9968934840746495D+00    
    x( 2) =  -0.9836681232797472D+00    
    x( 3) =  -0.9600218649683075D+00    
    x( 4) =  -0.9262000474292743D+00    
    x( 5) =  -0.8825605357920526D+00    
    x( 6) =  -0.8295657623827684D+00    
    x( 7) =  -0.7677774321048262D+00    
    x( 8) =  -0.6978504947933158D+00    
    x( 9) =  -0.6205261829892429D+00    
    x(10) =  -0.5366241481420199D+00    
    x(11) =  -0.4470337695380892D+00    
    x(12) =  -0.3527047255308781D+00    
    x(13) =  -0.2546369261678899D+00    
    x(14) =  -0.1538699136085835D+00    
    x(15) =  -0.5147184255531770D-01
    x(16) =   0.5147184255531770D-01
    x(17) =   0.1538699136085835D+00    
    x(18) =   0.2546369261678899D+00    
    x(19) =   0.3527047255308781D+00    
    x(20) =   0.4470337695380892D+00    
    x(21) =   0.5366241481420199D+00    
    x(22) =   0.6205261829892429D+00    
    x(23) =   0.6978504947933158D+00    
    x(24) =   0.7677774321048262D+00    
    x(25) =   0.8295657623827684D+00    
    x(26) =   0.8825605357920526D+00    
    x(27) =   0.9262000474292743D+00    
    x(28) =   0.9600218649683075D+00    
    x(29) =   0.9836681232797472D+00    
    x(30) =   0.9968934840746495D+00    
 
    w( 1) =   0.7968192496166648D-02
    w( 2) =   0.1846646831109099D-01
    w( 3) =   0.2878470788332330D-01
    w( 4) =   0.3879919256962704D-01
    w( 5) =   0.4840267283059405D-01
    w( 6) =   0.5749315621761905D-01
    w( 7) =   0.6597422988218052D-01
    w( 8) =   0.7375597473770516D-01
    w( 9) =   0.8075589522942023D-01
    w(10) =   0.8689978720108314D-01
    w(11) =   0.9212252223778619D-01
    w(12) =   0.9636873717464424D-01
    w(13) =   0.9959342058679524D-01
    w(14) =   0.1017623897484056D+00    
    w(15) =   0.1028526528935587D+00    
    w(16) =   0.1028526528935587D+00    
    w(17) =   0.1017623897484056D+00    
    w(18) =   0.9959342058679524D-01
    w(19) =   0.9636873717464424D-01
    w(20) =   0.9212252223778619D-01
    w(21) =   0.8689978720108314D-01
    w(22) =   0.8075589522942023D-01
    w(23) =   0.7375597473770516D-01
    w(24) =   0.6597422988218052D-01
    w(25) =   0.5749315621761905D-01
    w(26) =   0.4840267283059405D-01
    w(27) =   0.3879919256962704D-01
    w(28) =   0.2878470788332330D-01
    w(29) =   0.1846646831109099D-01
    w(30) =   0.7968192496166648D-02

  else if ( n == 31 ) then

    x( 1) =  -0.99708748181947707454263838179654D+00   
    x( 2) =  -0.98468590966515248400211329970113D+00
    x( 3) =  -0.96250392509294966178905249675943D+00
    x( 4) =  -0.93075699789664816495694576311725D+00
    x( 5) =  -0.88976002994827104337419200908023D+00
    x( 6) =  -0.83992032014626734008690453594388D+00
    x( 7) =  -0.78173314841662494040636002019484D+00
    x( 8) =  -0.71577678458685328390597086536649D+00
    x( 9) =  -0.64270672292426034618441820323250D+00
    x(10) =  -0.56324916140714926272094492359516D+00
    x(11) =  -0.47819378204490248044059403935649D+00
    x(12) =  -0.38838590160823294306135146128752D+00
    x(13) =  -0.29471806998170161661790389767170D+00
    x(14) =  -0.19812119933557062877241299603283D+00
    x(15) =  -0.99555312152341520325174790118941D-01
    x(16) =   0.00000000000000000000000000000000D+00  
    x(17) =   0.99555312152341520325174790118941D-01
    x(18) =   0.19812119933557062877241299603283D+00    
    x(19) =   0.29471806998170161661790389767170D+00    
    x(20) =   0.38838590160823294306135146128752D+00    
    x(21) =   0.47819378204490248044059403935649D+00    
    x(22) =   0.56324916140714926272094492359516D+00    
    x(23) =   0.64270672292426034618441820323250D+00    
    x(24) =   0.71577678458685328390597086536649D+00    
    x(25) =   0.78173314841662494040636002019484D+00    
    x(26) =   0.83992032014626734008690453594388D+00    
    x(27) =   0.88976002994827104337419200908023D+00    
    x(28) =   0.93075699789664816495694576311725D+00    
    x(29) =   0.96250392509294966178905249675943D+00    
    x(30) =   0.98468590966515248400211329970113D+00    
    x(31) =   0.99708748181947707454263838179654D+00    
 
    w( 1) =   0.74708315792487746093913218970494D-02
    w( 2) =   0.17318620790310582463552990782414D-01
    w( 3) =   0.27009019184979421800608642617676D-01
    w( 4) =   0.36432273912385464024392008749009D-01
    w( 5) =   0.45493707527201102902315857856518D-01
    w( 6) =   0.54103082424916853711666259085477D-01
    w( 7) =   0.62174786561028426910343543686657D-01
    w( 8) =   0.69628583235410366167756126255124D-01
    w( 9) =   0.76390386598776616426357674901331D-01
    w(10) =   0.82392991761589263903823367431962D-01
    w(11) =   0.87576740608477876126198069695333D-01
    w(12) =   0.91890113893641478215362871607150D-01
    w(13) =   0.95290242912319512807204197487597D-01
    w(14) =   0.97743335386328725093474010978997D-01
    w(15) =   0.99225011226672307874875514428615D-01
    w(16) =   0.99720544793426451427533833734349D-01
    w(17) =   0.99225011226672307874875514428615D-01
    w(18) =   0.97743335386328725093474010978997D-01
    w(19) =   0.95290242912319512807204197487597D-01
    w(20) =   0.91890113893641478215362871607150D-01
    w(21) =   0.87576740608477876126198069695333D-01
    w(22) =   0.82392991761589263903823367431962D-01
    w(23) =   0.76390386598776616426357674901331D-01
    w(24) =   0.69628583235410366167756126255124D-01
    w(25) =   0.62174786561028426910343543686657D-01
    w(26) =   0.54103082424916853711666259085477D-01
    w(27) =   0.45493707527201102902315857856518D-01
    w(28) =   0.36432273912385464024392008749009D-01
    w(29) =   0.27009019184979421800608642617676D-01
    w(30) =   0.17318620790310582463552990782414D-01
    w(31) =   0.74708315792487746093913218970494D-02

  else if ( n == 32 ) then

    x(1) =  - 0.997263861849481563544981128665D+00
    x(2) =  - 0.985611511545268335400175044631D+00
    x(3) =  - 0.964762255587506430773811928118D+00
    x(4) =  - 0.934906075937739689170919134835D+00
    x(5) =  - 0.896321155766052123965307243719D+00
    x(6) =  - 0.849367613732569970133693004968D+00
    x(7) =  - 0.794483795967942406963097298970D+00
    x(8) =  - 0.732182118740289680387426665091D+00
    x(9) =  - 0.663044266930215200975115168663D+00
    x(10) = - 0.587715757240762329040745476402D+00
    x(11) = - 0.506899908932229390023747474378D+00
    x(12) = - 0.421351276130635345364119436172D+00
    x(13) = - 0.331868602282127649779916805730D+00
    x(14) = - 0.239287362252137074544603209166D+00
    x(15) = - 0.144471961582796493485186373599D+00
    x(16) = - 0.483076656877383162348125704405D-01
    x(17) =   0.483076656877383162348125704405D-01
    x(18) =   0.144471961582796493485186373599D+00
    x(19) =   0.239287362252137074544603209166D+00
    x(20) =   0.331868602282127649779916805730D+00
    x(21) =   0.421351276130635345364119436172D+00
    x(22) =   0.506899908932229390023747474378D+00
    x(23) =   0.587715757240762329040745476402D+00
    x(24) =   0.663044266930215200975115168663D+00
    x(25) =   0.732182118740289680387426665091D+00
    x(26) =   0.794483795967942406963097298970D+00
    x(27) =   0.849367613732569970133693004968D+00
    x(28) =   0.896321155766052123965307243719D+00
    x(29) =   0.934906075937739689170919134835D+00
    x(30) =   0.964762255587506430773811928118D+00
    x(31) =   0.985611511545268335400175044631D+00
    x(32) =   0.997263861849481563544981128665D+00

    w(1) =  0.701861000947009660040706373885D-02
    w(2) =  0.162743947309056706051705622064D-01
    w(3) =  0.253920653092620594557525897892D-01
    w(4) =  0.342738629130214331026877322524D-01
    w(5) =  0.428358980222266806568786466061D-01
    w(6) =  0.509980592623761761961632446895D-01
    w(7) =  0.586840934785355471452836373002D-01
    w(8) =  0.658222227763618468376500637069D-01
    w(9) =  0.723457941088485062253993564785D-01
    w(10) = 0.781938957870703064717409188283D-01
    w(11) = 0.833119242269467552221990746043D-01
    w(12) = 0.876520930044038111427714627518D-01
    w(13) = 0.911738786957638847128685771116D-01
    w(14) = 0.938443990808045656391802376681D-01
    w(15) = 0.956387200792748594190820022041D-01
    w(16) = 0.965400885147278005667648300636D-01
    w(17) = 0.965400885147278005667648300636D-01
    w(18) = 0.956387200792748594190820022041D-01
    w(19) = 0.938443990808045656391802376681D-01
    w(20) = 0.911738786957638847128685771116D-01
    w(21) = 0.876520930044038111427714627518D-01
    w(22) = 0.833119242269467552221990746043D-01
    w(23) = 0.781938957870703064717409188283D-01
    w(24) = 0.723457941088485062253993564785D-01
    w(25) = 0.658222227763618468376500637069D-01
    w(26) = 0.586840934785355471452836373002D-01
    w(27) = 0.509980592623761761961632446895D-01
    w(28) = 0.428358980222266806568786466061D-01
    w(29) = 0.342738629130214331026877322524D-01
    w(30) = 0.253920653092620594557525897892D-01
    w(31) = 0.162743947309056706051705622064D-01
    w(32) = 0.701861000947009660040706373885D-02

  else if ( n == 33 ) then

    x( 1) =  -0.9974246942464552D+00    
    x( 2) =  -0.9864557262306425D+00
    x( 3) =  -0.9668229096899927D+00
    x( 4) =  -0.9386943726111684D+00    
    x( 5) =  -0.9023167677434336D+00    
    x( 6) =  -0.8580096526765041D+00    
    x( 7) =  -0.8061623562741665D+00    
    x( 8) =  -0.7472304964495622D+00    
    x( 9) =  -0.6817319599697428D+00    
    x(10) =  -0.6102423458363790D+00    
    x(11) =  -0.5333899047863476D+00    
    x(12) =  -0.4518500172724507D+00    
    x(13) =  -0.3663392577480734D+00    
    x(14) =  -0.2776090971524970D+00    
    x(15) =  -0.1864392988279916D+00    
    x(16) =  -0.09363106585473338D+00
    x(17) =   0.000000000000000D+00
    x(18) =   0.09363106585473338D+00
    x(19) =   0.1864392988279916D+00    
    x(20) =   0.2776090971524970D+00    
    x(21) =   0.3663392577480734D+00    
    x(22) =   0.4518500172724507D+00    
    x(23) =   0.5333899047863476D+00    
    x(24) =   0.6102423458363790D+00    
    x(25) =   0.6817319599697428D+00    
    x(26) =   0.7472304964495622D+00    
    x(27) =   0.8061623562741665D+00    
    x(28) =   0.8580096526765041D+00    
    x(29) =   0.9023167677434336D+00    
    x(30) =   0.9386943726111684D+00    
    x(31) =   0.9668229096899927D+00    
    x(32) =   0.9864557262306425D+00    
    x(33) =   0.9974246942464552D+00    
 
    w( 1) =   0.6606227847587558D-02
    w( 2) =   0.1532170151293465D-01
    w( 3) =   0.2391554810174960D-01
    w( 4) =   0.3230035863232891D-01
    w( 5) =   0.4040154133166965D-01
    w( 6) =   0.4814774281871162D-01
    w( 7) =   0.5547084663166357D-01
    w( 8) =   0.6230648253031755D-01
    w( 9) =   0.6859457281865676D-01
    w(10) =   0.7427985484395420D-01
    w(11) =   0.7931236479488685D-01
    w(12) =   0.8364787606703869D-01
    w(13) =   0.8724828761884425D-01
    w(14) =   0.9008195866063859D-01
    w(15) =   0.9212398664331678D-01
    w(16) =   0.9335642606559612D-01
    w(17) =   0.9376844616020999D-01
    w(18) =   0.9335642606559612D-01
    w(19) =   0.9212398664331678D-01
    w(20) =   0.9008195866063859D-01
    w(21) =   0.8724828761884425D-01
    w(22) =   0.8364787606703869D-01
    w(23) =   0.7931236479488685D-01
    w(24) =   0.7427985484395420D-01
    w(25) =   0.6859457281865676D-01
    w(26) =   0.6230648253031755D-01
    w(27) =   0.5547084663166357D-01
    w(28) =   0.4814774281871162D-01
    w(29) =   0.4040154133166965D-01
    w(30) =   0.3230035863232891D-01
    w(31) =   0.2391554810174960D-01
    w(32) =   0.1532170151293465D-01
    w(33) =   0.6606227847587558D-02

  else if ( n == 63 ) then

    x( 1) =  -0.99928298402912378050701628988630D+00    
    x( 2) =  -0.99622401277797010860209018267357D+00    
    x( 3) =  -0.99072854689218946681089469460884D+00    
    x( 4) =  -0.98280881059372723486251140727639D+00    
    x( 5) =  -0.97248403469757002280196067864927D+00    
    x( 6) =  -0.95977944975894192707035416626398D+00    
    x( 7) =  -0.94472613404100980296637531962798D+00    
    x( 8) =  -0.92736092062184320544703138132518D+00    
    x( 9) =  -0.90772630277853155803695313291596D+00    
    x(10) =  -0.88587032850785342629029845731337D+00    
    x(11) =  -0.86184648236412371953961183943106D+00    
    x(12) =  -0.83571355431950284347180776961571D+00    
    x(13) =  -0.80753549577345676005146598636324D+00    
    x(14) =  -0.77738126299037233556333018991104D+00    
    x(15) =  -0.74532464831784741782932166103759D+00    
    x(16) =  -0.71144409958484580785143153770401D+00    
    x(17) =  -0.67582252811498609013110331596954D+00    
    x(18) =  -0.63854710582136538500030695387338D+00    
    x(19) =  -0.59970905187762523573900892686880D+00    
    x(20) =  -0.55940340948628501326769780007005D+00    
    x(21) =  -0.51772881329003324812447758452632D+00    
    x(22) =  -0.47478724799480439992221230985149D+00    
    x(23) =  -0.43068379879511160066208893391863D+00    
    x(24) =  -0.38552639421224789247761502227440D+00    
    x(25) =  -0.33942554197458440246883443159432D+00    
    x(26) =  -0.29249405858625144003615715555067D+00    
    x(27) =  -0.24484679324595336274840459392483D+00    
    x(28) =  -0.19660034679150668455762745706572D+00    
    x(29) =  -0.14787278635787196856983909655297D+00    
    x(30) =  -0.98783356446945279529703669453922D-01
    x(31) =  -0.49452187116159627234233818051808D-01
    x(32) =    0.0000000000000000000000000000000D+00    
    x(33) =   0.49452187116159627234233818051808D-01
    x(34) =   0.98783356446945279529703669453922D-01
    x(35) =   0.14787278635787196856983909655297D+00    
    x(36) =   0.19660034679150668455762745706572D+00    
    x(37) =   0.24484679324595336274840459392483D+00    
    x(38) =   0.29249405858625144003615715555067D+00    
    x(39) =   0.33942554197458440246883443159432D+00    
    x(40) =   0.38552639421224789247761502227440D+00    
    x(41) =   0.43068379879511160066208893391863D+00    
    x(42) =   0.47478724799480439992221230985149D+00    
    x(43) =   0.51772881329003324812447758452632D+00    
    x(44) =   0.55940340948628501326769780007005D+00    
    x(45) =   0.59970905187762523573900892686880D+00    
    x(46) =   0.63854710582136538500030695387338D+00    
    x(47) =   0.67582252811498609013110331596954D+00    
    x(48) =   0.71144409958484580785143153770401D+00    
    x(49) =   0.74532464831784741782932166103759D+00    
    x(50) =   0.77738126299037233556333018991104D+00    
    x(51) =   0.80753549577345676005146598636324D+00    
    x(52) =   0.83571355431950284347180776961571D+00    
    x(53) =   0.86184648236412371953961183943106D+00    
    x(54) =   0.88587032850785342629029845731337D+00    
    x(55) =   0.90772630277853155803695313291596D+00    
    x(56) =   0.92736092062184320544703138132518D+00    
    x(57) =   0.94472613404100980296637531962798D+00    
    x(58) =   0.95977944975894192707035416626398D+00    
    x(59) =   0.97248403469757002280196067864927D+00    
    x(60) =   0.98280881059372723486251140727639D+00    
    x(61) =   0.99072854689218946681089469460884D+00    
    x(62) =   0.99622401277797010860209018267357D+00    
    x(63) =   0.99928298402912378050701628988630D+00

    w( 1) =   0.18398745955770837880499331680577D-02
    w( 2) =   0.42785083468637618661951422543371D-02
    w( 3) =   0.67102917659601362519069109850892D-02
    w( 4) =   0.91259686763266563540586445877022D-02
    w( 5) =   0.11519376076880041750750606118707D-01
    w( 6) =   0.13884612616115610824866086365937D-01
    w( 7) =   0.16215878410338338882283672974995D-01
    w( 8) =   0.18507464460161270409260545805144D-01
    w( 9) =   0.20753761258039090775341953421471D-01
    w(10) =   0.22949271004889933148942319561770D-01
    w(11) =   0.25088620553344986618630138068443D-01
    w(12) =   0.27166574359097933225189839439413D-01
    w(13) =   0.29178047208280526945551502154029D-01
    w(14) =   0.31118116622219817508215988557189D-01
    w(15) =   0.32982034883779341765683179672459D-01
    w(16) =   0.34765240645355877697180504642788D-01
    w(17) =   0.36463370085457289630452409787542D-01
    w(18) =   0.38072267584349556763638324927889D-01
    w(19) =   0.39587995891544093984807928149202D-01
    w(20) =   0.41006845759666398635110037009072D-01
    w(21) =   0.42325345020815822982505485403028D-01
    w(22) =   0.43540267083027590798964315704401D-01
    w(23) =   0.44648638825941395370332669516813D-01
    w(24) =   0.45647747876292608685885992608542D-01
    w(25) =   0.46535149245383696510395418746953D-01
    w(26) =   0.47308671312268919080604988338844D-01
    w(27) =   0.47966421137995131411052756195132D-01
    w(28) =   0.48506789097883847864090099145802D-01
    w(29) =   0.48928452820511989944709361549215D-01
    w(30) =   0.49230380423747560785043116988145D-01
    w(31) =   0.49411833039918178967039646116705D-01
    w(32) =   0.49472366623931020888669360420926D-01
    w(33) =   0.49411833039918178967039646116705D-01
    w(34) =   0.49230380423747560785043116988145D-01
    w(35) =   0.48928452820511989944709361549215D-01
    w(36) =   0.48506789097883847864090099145802D-01
    w(37) =   0.47966421137995131411052756195132D-01
    w(38) =   0.47308671312268919080604988338844D-01
    w(39) =   0.46535149245383696510395418746953D-01
    w(40) =   0.45647747876292608685885992608542D-01
    w(41) =   0.44648638825941395370332669516813D-01
    w(42) =   0.43540267083027590798964315704401D-01
    w(43) =   0.42325345020815822982505485403028D-01
    w(44) =   0.41006845759666398635110037009072D-01
    w(45) =   0.39587995891544093984807928149202D-01
    w(46) =   0.38072267584349556763638324927889D-01
    w(47) =   0.36463370085457289630452409787542D-01
    w(48) =   0.34765240645355877697180504642788D-01
    w(49) =   0.32982034883779341765683179672459D-01
    w(50) =   0.31118116622219817508215988557189D-01
    w(51) =   0.29178047208280526945551502154029D-01
    w(52) =   0.27166574359097933225189839439413D-01
    w(53) =   0.25088620553344986618630138068443D-01
    w(54) =   0.22949271004889933148942319561770D-01
    w(55) =   0.20753761258039090775341953421471D-01
    w(56) =   0.18507464460161270409260545805144D-01
    w(57) =   0.16215878410338338882283672974995D-01
    w(58) =   0.13884612616115610824866086365937D-01
    w(59) =   0.11519376076880041750750606118707D-01
    w(60) =   0.91259686763266563540586445877022D-02
    w(61) =   0.67102917659601362519069109850892D-02
    w(62) =   0.42785083468637618661951422543371D-02
    w(63) =   0.18398745955770837880499331680577D-02
 
  else if ( n == 64 ) then

    x(1) =  - 0.999305041735772139456905624346D+00
    x(2) =  - 0.996340116771955279346924500676D+00
    x(3) =  - 0.991013371476744320739382383443D+00
    x(4) =  - 0.983336253884625956931299302157D+00
    x(5) =  - 0.973326827789910963741853507352D+00
    x(6) =  - 0.961008799652053718918614121897D+00
    x(7) =  - 0.946411374858402816062481491347D+00
    x(8) =  - 0.929569172131939575821490154559D+00
    x(9) =  - 0.910522137078502805756380668008D+00
    x(10) = - 0.889315445995114105853404038273D+00
    x(11) = - 0.865999398154092819760783385070D+00
    x(12) = - 0.840629296252580362751691544696D+00
    x(13) = - 0.813265315122797559741923338086D+00
    x(14) = - 0.783972358943341407610220525214D+00
    x(15) = - 0.752819907260531896611863774886D+00
    x(16) = - 0.719881850171610826848940217832D+00
    x(17) = - 0.685236313054233242563558371031D+00
    x(18) = - 0.648965471254657339857761231993D+00
    x(19) = - 0.611155355172393250248852971019D+00
    x(20) = - 0.571895646202634034283878116659D+00
    x(21) = - 0.531279464019894545658013903544D+00
    x(22) = - 0.489403145707052957478526307022D+00
    x(23) = - 0.446366017253464087984947714759D+00
    x(24) = - 0.402270157963991603695766771260D+00
    x(25) = - 0.357220158337668115950442615046D+00
    x(26) = - 0.311322871990210956157512698560D+00
    x(27) = - 0.264687162208767416373964172510D+00
    x(28) = - 0.217423643740007084149648748989D+00
    x(29) = - 0.169644420423992818037313629748D+00
    x(30) = - 0.121462819296120554470376463492D+00
    x(31) = - 0.729931217877990394495429419403D-01
    x(32) = - 0.243502926634244325089558428537D-01
    x(33) =   0.243502926634244325089558428537D-01
    x(34) =   0.729931217877990394495429419403D-01
    x(35) =   0.121462819296120554470376463492D+00
    x(36) =   0.169644420423992818037313629748D+00
    x(37) =   0.217423643740007084149648748989D+00
    x(38) =   0.264687162208767416373964172510D+00
    x(39) =   0.311322871990210956157512698560D+00
    x(40) =   0.357220158337668115950442615046D+00
    x(41) =   0.402270157963991603695766771260D+00
    x(42) =   0.446366017253464087984947714759D+00
    x(43) =   0.489403145707052957478526307022D+00
    x(44) =   0.531279464019894545658013903544D+00
    x(45) =   0.571895646202634034283878116659D+00
    x(46) =   0.611155355172393250248852971019D+00
    x(47) =   0.648965471254657339857761231993D+00
    x(48) =   0.685236313054233242563558371031D+00
    x(49) =   0.719881850171610826848940217832D+00
    x(50) =   0.752819907260531896611863774886D+00
    x(51) =   0.783972358943341407610220525214D+00
    x(52) =   0.813265315122797559741923338086D+00
    x(53) =   0.840629296252580362751691544696D+00
    x(54) =   0.865999398154092819760783385070D+00
    x(55) =   0.889315445995114105853404038273D+00
    x(56) =   0.910522137078502805756380668008D+00
    x(57) =   0.929569172131939575821490154559D+00
    x(58) =   0.946411374858402816062481491347D+00
    x(59) =   0.961008799652053718918614121897D+00
    x(60) =   0.973326827789910963741853507352D+00
    x(61) =   0.983336253884625956931299302157D+00
    x(62) =   0.991013371476744320739382383443D+00
    x(63) =   0.996340116771955279346924500676D+00
    x(64) =   0.999305041735772139456905624346D+00

    w(1) =  0.178328072169643294729607914497D-02
    w(2) =  0.414703326056246763528753572855D-02
    w(3) =  0.650445796897836285611736039998D-02
    w(4) =  0.884675982636394772303091465973D-02
    w(5) =  0.111681394601311288185904930192D-01
    w(6) =  0.134630478967186425980607666860D-01
    w(7) =  0.157260304760247193219659952975D-01
    w(8) =  0.179517157756973430850453020011D-01
    w(9) =  0.201348231535302093723403167285D-01
    w(10) = 0.222701738083832541592983303842D-01
    w(11) = 0.243527025687108733381775504091D-01
    w(12) = 0.263774697150546586716917926252D-01
    w(13) = 0.283396726142594832275113052002D-01
    w(14) = 0.302346570724024788679740598195D-01
    w(15) = 0.320579283548515535854675043479D-01
    w(16) = 0.338051618371416093915654821107D-01
    w(17) = 0.354722132568823838106931467152D-01
    w(18) = 0.370551285402400460404151018096D-01
    w(19) = 0.385501531786156291289624969468D-01
    w(20) = 0.399537411327203413866569261283D-01
    w(21) = 0.412625632426235286101562974736D-01
    w(22) = 0.424735151236535890073397679088D-01
    w(23) = 0.435837245293234533768278609737D-01
    w(24) = 0.445905581637565630601347100309D-01
    w(25) = 0.454916279274181444797709969713D-01
    w(26) = 0.462847965813144172959532492323D-01
    w(27) = 0.469681828162100173253262857546D-01
    w(28) = 0.475401657148303086622822069442D-01
    w(29) = 0.479993885964583077281261798713D-01
    w(30) = 0.483447622348029571697695271580D-01
    w(31) = 0.485754674415034269347990667840D-01
    w(32) = 0.486909570091397203833653907347D-01
    w(33) = 0.486909570091397203833653907347D-01
    w(34) = 0.485754674415034269347990667840D-01
    w(35) = 0.483447622348029571697695271580D-01
    w(36) = 0.479993885964583077281261798713D-01
    w(37) = 0.475401657148303086622822069442D-01
    w(38) = 0.469681828162100173253262857546D-01
    w(39) = 0.462847965813144172959532492323D-01
    w(40) = 0.454916279274181444797709969713D-01
    w(41) = 0.445905581637565630601347100309D-01
    w(42) = 0.435837245293234533768278609737D-01
    w(43) = 0.424735151236535890073397679088D-01
    w(44) = 0.412625632426235286101562974736D-01
    w(45) = 0.399537411327203413866569261283D-01
    w(46) = 0.385501531786156291289624969468D-01
    w(47) = 0.370551285402400460404151018096D-01
    w(48) = 0.354722132568823838106931467152D-01
    w(49) = 0.338051618371416093915654821107D-01
    w(50) = 0.320579283548515535854675043479D-01
    w(51) = 0.302346570724024788679740598195D-01
    w(52) = 0.283396726142594832275113052002D-01
    w(53) = 0.263774697150546586716917926252D-01
    w(54) = 0.243527025687108733381775504091D-01
    w(55) = 0.222701738083832541592983303842D-01
    w(56) = 0.201348231535302093723403167285D-01
    w(57) = 0.179517157756973430850453020011D-01
    w(58) = 0.157260304760247193219659952975D-01
    w(59) = 0.134630478967186425980607666860D-01
    w(60) = 0.111681394601311288185904930192D-01
    w(61) = 0.884675982636394772303091465973D-02
    w(62) = 0.650445796897836285611736039998D-02
    w(63) = 0.414703326056246763528753572855D-02
    w(64) = 0.178328072169643294729607914497D-02

  else if ( n == 65 ) then

    x( 1) =  -0.9993260970754129D+00    
    x( 2) =  -0.9964509480618492D+00    
    x( 3) =  -0.9912852761768016D+00    
    x( 4) =  -0.9838398121870350D+00    
    x( 5) =  -0.9741315398335512D+00    
    x( 6) =  -0.9621827547180553D+00    
    x( 7) =  -0.9480209281684076D+00    
    x( 8) =  -0.9316786282287494D+00    
    x( 9) =  -0.9131934405428462D+00    
    x(10) =  -0.8926078805047389D+00    
    x(11) =  -0.8699692949264071D+00    
    x(12) =  -0.8453297528999303D+00    
    x(13) =  -0.8187459259226514D+00    
    x(14) =  -0.7902789574921218D+00    
    x(15) =  -0.7599943224419998D+00    
    x(16) =  -0.7279616763294247D+00    
    x(17) =  -0.6942546952139916D+00    
    x(18) =  -0.6589509061936252D+00    
    x(19) =  -0.6221315090854003D+00    
    x(20) =  -0.5838811896604873D+00    
    x(21) =  -0.5442879248622271D+00    
    x(22) =  -0.5034427804550069D+00    
    x(23) =  -0.4614397015691450D+00    
    x(24) =  -0.4183752966234090D+00    
    x(25) =  -0.3743486151220660D+00    
    x(26) =  -0.3294609198374864D+00    
    x(27) =  -0.2838154539022487D+00    
    x(28) =  -0.2375172033464168D+00    
    x(29) =  -0.1906726556261428D+00    
    x(30) =  -0.1433895546989752D+00    
    x(31) =  -0.9577665320919751D-01
    x(32) =  -0.4794346235317186D-01
    x(33) =    0.000000000000000D+00    
    x(34) =   0.4794346235317186D-01
    x(35) =   0.9577665320919751D-01
    x(36) =   0.1433895546989752D+00    
    x(37) =   0.1906726556261428D+00    
    x(38) =   0.2375172033464168D+00    
    x(39) =   0.2838154539022487D+00    
    x(40) =   0.3294609198374864D+00    
    x(41) =   0.3743486151220660D+00    
    x(42) =   0.4183752966234090D+00    
    x(43) =   0.4614397015691450D+00    
    x(44) =   0.5034427804550069D+00    
    x(45) =   0.5442879248622271D+00    
    x(46) =   0.5838811896604873D+00    
    x(47) =   0.6221315090854003D+00    
    x(48) =   0.6589509061936252D+00    
    x(49) =   0.6942546952139916D+00    
    x(50) =   0.7279616763294247D+00    
    x(51) =   0.7599943224419998D+00    
    x(52) =   0.7902789574921218D+00    
    x(53) =   0.8187459259226514D+00    
    x(54) =   0.8453297528999303D+00    
    x(55) =   0.8699692949264071D+00    
    x(56) =   0.8926078805047389D+00    
    x(57) =   0.9131934405428462D+00    
    x(58) =   0.9316786282287494D+00    
    x(59) =   0.9480209281684076D+00    
    x(60) =   0.9621827547180553D+00    
    x(61) =   0.9741315398335512D+00    
    x(62) =   0.9838398121870350D+00    
    x(63) =   0.9912852761768016D+00    
    x(64) =   0.9964509480618492D+00    
    x(65) =   0.9993260970754129D+00    
 
    w( 1) =   0.1729258251300218D-02
    w( 2) =   0.4021524172003703D-02
    w( 3) =   0.6307942578971821D-02
    w( 4) =   0.8580148266881443D-02
    w( 5) =   0.1083267878959798D-01
    w( 6) =   0.1306031163999490D-01
    w( 7) =   0.1525791214644825D-01
    w( 8) =   0.1742042199767025D-01
    w( 9) =   0.1954286583675005D-01
    w(10) =   0.2162036128493408D-01
    w(11) =   0.2364812969128723D-01
    w(12) =   0.2562150693803776D-01
    w(13) =   0.2753595408845034D-01
    w(14) =   0.2938706778931066D-01
    w(15) =   0.3117059038018911D-01
    w(16) =   0.3288241967636860D-01
    w(17) =   0.3451861839854901D-01
    w(18) =   0.3607542322556527D-01
    w(19) =   0.3754925344825770D-01
    w(20) =   0.3893671920405121D-01
    w(21) =   0.4023462927300549D-01
    w(22) =   0.4143999841724028D-01
    w(23) =   0.4255005424675579D-01
    w(24) =   0.4356224359580051D-01
    w(25) =   0.4447423839508296D-01
    w(26) =   0.4528394102630023D-01
    w(27) =   0.4598948914665173D-01
    w(28) =   0.4658925997223349D-01
    w(29) =   0.4708187401045461D-01
    w(30) =   0.4746619823288551D-01
    w(31) =   0.4774134868124067D-01
    w(32) =   0.4790669250049590D-01
    w(33) =   0.4796184939446662D-01
    w(34) =   0.4790669250049590D-01
    w(35) =   0.4774134868124067D-01
    w(36) =   0.4746619823288551D-01
    w(37) =   0.4708187401045461D-01
    w(38) =   0.4658925997223349D-01
    w(39) =   0.4598948914665173D-01
    w(40) =   0.4528394102630023D-01
    w(41) =   0.4447423839508296D-01
    w(42) =   0.4356224359580051D-01
    w(43) =   0.4255005424675579D-01
    w(44) =   0.4143999841724028D-01
    w(45) =   0.4023462927300549D-01
    w(46) =   0.3893671920405121D-01
    w(47) =   0.3754925344825770D-01
    w(48) =   0.3607542322556527D-01
    w(49) =   0.3451861839854901D-01
    w(50) =   0.3288241967636860D-01
    w(51) =   0.3117059038018911D-01
    w(52) =   0.2938706778931066D-01
    w(53) =   0.2753595408845034D-01
    w(54) =   0.2562150693803776D-01
    w(55) =   0.2364812969128723D-01
    w(56) =   0.2162036128493408D-01
    w(57) =   0.1954286583675005D-01
    w(58) =   0.1742042199767025D-01
    w(59) =   0.1525791214644825D-01
    w(60) =   0.1306031163999490D-01
    w(61) =   0.1083267878959798D-01
    w(62) =   0.8580148266881443D-02
    w(63) =   0.6307942578971821D-02
    w(64) =   0.4021524172003703D-02
    w(65) =   0.1729258251300218D-02
    
  else if ( n == 127 ) then
  
    x(  1) =  -0.99982213041530614629963254927125D+00    
    x(  2) =  -0.99906293435531189513828920479421D+00    
    x(  3) =  -0.99769756618980462107441703193392D+00    
    x(  4) =  -0.99572655135202722663543337085008D+00    
    x(  5) =  -0.99315104925451714736113079489080D+00    
    x(  6) =  -0.98997261459148415760778669967548D+00    
    x(  7) =  -0.98619317401693166671043833175407D+00    
    x(  8) =  -0.98181502080381411003346312451200D+00    
    x(  9) =  -0.97684081234307032681744391886221D+00    
    x( 10) =  -0.97127356816152919228894689830512D+00    
    x( 11) =  -0.96511666794529212109082507703391D+00    
    x( 12) =  -0.95837384942523877114910286998060D+00    
    x( 13) =  -0.95104920607788031054790764659636D+00    
    x( 14) =  -0.94314718462481482734544963026201D+00    
    x( 15) =  -0.93467258232473796857363487794906D+00    
    x( 16) =  -0.92563054405623384912746466814259D+00    
    x( 17) =  -0.91602655919146580931308861741716D+00    
    x( 18) =  -0.90586645826182138280246131760282D+00    
    x( 19) =  -0.89515640941708370896904382642451D+00    
    x( 20) =  -0.88390291468002656994525794802849D+00    
    x( 21) =  -0.87211280599856071141963753428864D+00    
    x( 22) =  -0.85979324109774080981203134414483D+00    
    x( 23) =  -0.84695169913409759845333931085437D+00    
    x( 24) =  -0.83359597615489951437955716480123D+00    
    x( 25) =  -0.81973418036507867415511910167470D+00    
    x( 26) =  -0.80537472720468021466656079404644D+00    
    x( 27) =  -0.79052633423981379994544995252740D+00    
    x( 28) =  -0.77519801587020238244496276354566D+00    
    x( 29) =  -0.75939907785653667155666366659810D+00    
    x( 30) =  -0.74313911167095451292056688997595D+00   
    x( 31) =  -0.72642798867407268553569290153270D+00    
    x( 32) =  -0.70927585412210456099944463906757D+00    
    x( 33) =  -0.69169312100770067015644143286666D+00    
    x( 34) =  -0.67369046373825048534668253831602D+00    
    x( 35) =  -0.65527881165548263027676505156852D+00    
    x( 36) =  -0.63646934240029724134760815684175D+00    
    x( 37) =  -0.61727347512685828385763916340822D+00    
    x( 38) =  -0.59770286357006522938441201887478D+00    
    x( 39) =  -0.57776938897061258000325165713764D+00    
    x( 40) =  -0.55748515286193223292186190687872D+00    
    x( 41) =  -0.53686246972339756745816636353452D+00    
    x( 42) =  -0.51591385950424935727727729906662D+00    
    x( 43) =  -0.49465204002278211739494017368636D+00    
    x( 44) =  -0.47308991924540524164509989939699D+00    
    x( 45) =  -0.45124058745026622733189858020729D+00    
    x( 46) =  -0.42911730928019337626254405355418D+00    
    x( 47) =  -0.40673351568978256340867288124339D+00    
    x( 48) =  -0.38410279579151693577907781452239D+00    
    x( 49) =  -0.36123888860586970607092484346723D+00    
    x( 50) =  -0.33815567472039850137600027657095D+00    
    x( 51) =  -0.31486716786289498148601475374890D+00    
    x( 52) =  -0.29138750639370562079451875284568D+00    
    x( 53) =  -0.26773094472238862088834352027938D+00    
    x( 54) =  -0.24391184465391785797071324453138D+00    
    x( 55) =  -0.21994466666968754245452337866940D+00    
    x( 56) =  -0.19584396114861085150428162519610D+00    
    x( 57) =  -0.17162435953364216500834492248954D+00    
    x( 58) =  -0.14730056544908566938932929319807D+00    
    x( 59) =  -0.12288734577408297172603365288567D+00    
    x( 60) =  -0.98399521677698970751091751509101D-01
    x( 61) =  -0.73851959621048545273440409360569D-01
    x( 62) =  -0.49259562331926630315379321821927D-01
    x( 63) =  -0.24637259757420944614897071846088D-01
    x( 64) =   0.00000000000000000000000000000000D+00    
    x( 65) =   0.24637259757420944614897071846088D-01
    x( 66) =   0.49259562331926630315379321821927D-01
    x( 67) =   0.73851959621048545273440409360569D-01
    x( 68) =   0.98399521677698970751091751509101D-01
    x( 69) =   0.12288734577408297172603365288567D+00    
    x( 70) =   0.14730056544908566938932929319807D+00    
    x( 71) =   0.17162435953364216500834492248954D+00    
    x( 72) =   0.19584396114861085150428162519610D+00    
    x( 73) =   0.21994466666968754245452337866940D+00    
    x( 74) =   0.24391184465391785797071324453138D+00    
    x( 75) =   0.26773094472238862088834352027938D+00    
    x( 76) =   0.29138750639370562079451875284568D+00    
    x( 77) =   0.31486716786289498148601475374890D+00    
    x( 78) =   0.33815567472039850137600027657095D+00    
    x( 79) =   0.36123888860586970607092484346723D+00    
    x( 80) =   0.38410279579151693577907781452239D+00    
    x( 81) =   0.40673351568978256340867288124339D+00    
    x( 82) =   0.42911730928019337626254405355418D+00    
    x( 83) =   0.45124058745026622733189858020729D+00    
    x( 84) =   0.47308991924540524164509989939699D+00    
    x( 85) =   0.49465204002278211739494017368636D+00    
    x( 86) =   0.51591385950424935727727729906662D+00    
    x( 87) =   0.53686246972339756745816636353452D+00    
    x( 88) =   0.55748515286193223292186190687872D+00    
    x( 89) =   0.57776938897061258000325165713764D+00    
    x( 90) =   0.59770286357006522938441201887478D+00    
    x( 91) =   0.61727347512685828385763916340822D+00    
    x( 92) =   0.63646934240029724134760815684175D+00    
    x( 93) =   0.65527881165548263027676505156852D+00    
    x( 94) =   0.67369046373825048534668253831602D+00    
    x( 95) =   0.69169312100770067015644143286666D+00   
    x( 96) =   0.70927585412210456099944463906757D+00    
    x( 97) =   0.72642798867407268553569290153270D+00    
    x( 98) =   0.74313911167095451292056688997595D+00    
    x( 99) =   0.75939907785653667155666366659810D+00    
    x(100) =   0.77519801587020238244496276354566D+00    
    x(101) =   0.79052633423981379994544995252740D+00    
    x(102) =   0.80537472720468021466656079404644D+00    
    x(103) =   0.81973418036507867415511910167470D+00    
    x(104) =   0.83359597615489951437955716480123D+00    
    x(105) =   0.84695169913409759845333931085437D+00    
    x(106) =   0.85979324109774080981203134414483D+00    
    x(107) =   0.87211280599856071141963753428864D+00    
    x(108) =   0.88390291468002656994525794802849D+00    
    x(109) =   0.89515640941708370896904382642451D+00    
    x(110) =   0.90586645826182138280246131760282D+00    
    x(111) =   0.91602655919146580931308861741716D+00    
    x(112) =   0.92563054405623384912746466814259D+00    
    x(113) =   0.93467258232473796857363487794906D+00    
    x(114) =   0.94314718462481482734544963026201D+00    
    x(115) =   0.95104920607788031054790764659636D+00    
    x(116) =   0.95837384942523877114910286998060D+00    
    x(117) =   0.96511666794529212109082507703391D+00    
    x(118) =   0.97127356816152919228894689830512D+00    
    x(119) =   0.97684081234307032681744391886221D+00    
    x(120) =   0.98181502080381411003346312451200D+00    
    x(121) =   0.98619317401693166671043833175407D+00    
    x(122) =   0.98997261459148415760778669967548D+00    
    x(123) =   0.99315104925451714736113079489080D+00    
    x(124) =   0.99572655135202722663543337085008D+00    
    x(125) =   0.99769756618980462107441703193392D+00    
    x(126) =   0.99906293435531189513828920479421D+00    
    x(127) =   0.99982213041530614629963254927125D+00    

    w(  1) =   0.45645726109586654495731936146574D-03
    w(  2) =   0.10622766869538486959954760554099D-02
    w(  3) =   0.16683488125171936761028811985672D-02
    w(  4) =   0.22734860707492547802810838362671D-02
    w(  5) =   0.28772587656289004082883197417581D-02
    w(  6) =   0.34792893810051465908910894094105D-02
    w(  7) =   0.40792095178254605327114733456293D-02
    w(  8) =   0.46766539777779034772638165662478D-02
    w(  9) =   0.52712596565634400891303815906251D-02
    w( 10) =   0.58626653903523901033648343751367D-02
    w( 11) =   0.64505120486899171845442463868748D-02
    w( 12) =   0.70344427036681608755685893032552D-02
    w( 13) =   0.76141028256526859356393930849227D-02
    w( 14) =   0.81891404887415730817235884718726D-02
    w( 15) =   0.87592065795403145773316804234385D-02
    w( 16) =   0.93239550065309714787536985834029D-02
    w( 17) =   0.98830429087554914716648010899606D-02
    w( 18) =   0.10436130863141005225673171997668D-01
    w( 19) =   0.10982883090068975788799657376065D-01
    w( 20) =   0.11522967656921087154811609734510D-01
    w( 21) =   0.12056056679400848183529562144697D-01
    w( 22) =   0.12581826520465013101514365424172D-01
    w( 23) =   0.13099957986718627426172681912499D-01
    w( 24) =   0.13610136522139249906034237533759D-01
    w( 25) =   0.14112052399003395774044161633613D-01
    w( 26) =   0.14605400905893418351737288078952D-01
    w( 27) =   0.15089882532666922992635733981431D-01
    w( 28) =   0.15565203152273955098532590262975D-01
    w( 29) =   0.16031074199309941802254151842763D-01
    w( 30) =   0.16487212845194879399346060358146D-01
    w( 31) =   0.16933342169871654545878815295200D-01
    w( 32) =   0.17369191329918731922164721250350D-01
    w( 33) =   0.17794495722974774231027912900351D-01
    w( 34) =   0.18208997148375106468721469154479D-01
    w( 35) =   0.18612443963902310429440419898958D-01
    w( 36) =   0.19004591238555646611148901044533D-01
    w( 37) =   0.19385200901246454628112623489471D-01
    w( 38) =   0.19754041885329183081815217323169D-01
    w( 39) =   0.20110890268880247225644623956287D-01
    w( 40) =   0.20455529410639508279497065713301D-01
    w( 41) =   0.20787750081531811812652137291250D-01
    w( 42) =   0.21107350591688713643523847921658D-01
    w( 43) =   0.21414136912893259295449693233545D-01
    w( 44) =   0.21707922796373466052301324695331D-01
    w( 45) =   0.21988529885872983756478409758807D-01
    w( 46) =   0.22255787825930280235631416460158D-01
    w( 47) =   0.22509534365300608085694429903050D-01
    w( 48) =   0.22749615455457959852242553240982D-01
    w( 49) =   0.22975885344117206754377437838947D-01
    w( 50) =   0.23188206663719640249922582981729D-01
    w( 51) =   0.23386450514828194170722043496950D-01
    w( 52) =   0.23570496544381716050033676844306D-01
    w( 53) =   0.23740233018760777777714726703424D-01
    w( 54) =   0.23895556891620665983864481754172D-01
    w( 55) =   0.24036373866450369675132086026456D-01
    w( 56) =   0.24162598453819584716522917710986D-01
    w( 57) =   0.24274154023278979833195063936748D-01
    w( 58) =   0.24370972849882214952813561907241D-01
    w( 59) =   0.24452996155301467956140198471529D-01
    w( 60) =   0.24520174143511508275183033290175D-01
    w( 61) =   0.24572466031020653286354137335186D-01
    w( 62) =   0.24609840071630254092545634003360D-01
    w( 63) =   0.24632273575707679066033370218017D-01
    w( 64) =   0.24639752923961094419579417477503D-01
    w( 65) =   0.24632273575707679066033370218017D-01
    w( 66) =   0.24609840071630254092545634003360D-01
    w( 67) =   0.24572466031020653286354137335186D-01
    w( 68) =   0.24520174143511508275183033290175D-01
    w( 69) =   0.24452996155301467956140198471529D-01
    w( 70) =   0.24370972849882214952813561907241D-01
    w( 71) =   0.24274154023278979833195063936748D-01
    w( 72) =   0.24162598453819584716522917710986D-01
    w( 73) =   0.24036373866450369675132086026456D-01
    w( 74) =   0.23895556891620665983864481754172D-01
    w( 75) =   0.23740233018760777777714726703424D-01
    w( 76) =   0.23570496544381716050033676844306D-01
    w( 77) =   0.23386450514828194170722043496950D-01
    w( 78) =   0.23188206663719640249922582981729D-01
    w( 79) =   0.22975885344117206754377437838947D-01
    w( 80) =   0.22749615455457959852242553240982D-01
    w( 81) =   0.22509534365300608085694429903050D-01
    w( 82) =   0.22255787825930280235631416460158D-01
    w( 83) =   0.21988529885872983756478409758807D-01
    w( 84) =   0.21707922796373466052301324695331D-01
    w( 85) =   0.21414136912893259295449693233545D-01
    w( 86) =   0.21107350591688713643523847921658D-01
    w( 87) =   0.20787750081531811812652137291250D-01
    w( 88) =   0.20455529410639508279497065713301D-01
    w( 89) =   0.20110890268880247225644623956287D-01
    w( 90) =   0.19754041885329183081815217323169D-01
    w( 91) =   0.19385200901246454628112623489471D-01
    w( 92) =   0.19004591238555646611148901044533D-01
    w( 93) =   0.18612443963902310429440419898958D-01
    w( 94) =   0.18208997148375106468721469154479D-01
    w( 95) =   0.17794495722974774231027912900351D-01
    w( 96) =   0.17369191329918731922164721250350D-01
    w( 97) =   0.16933342169871654545878815295200D-01
    w( 98) =   0.16487212845194879399346060358146D-01
    w( 99) =   0.16031074199309941802254151842763D-01
    w(100) =   0.15565203152273955098532590262975D-01
    w(101) =   0.15089882532666922992635733981431D-01
    w(102) =   0.14605400905893418351737288078952D-01
    w(103) =   0.14112052399003395774044161633613D-01
    w(104) =   0.13610136522139249906034237533759D-01
    w(105) =   0.13099957986718627426172681912499D-01
    w(106) =   0.12581826520465013101514365424172D-01
    w(107) =   0.12056056679400848183529562144697D-01
    w(108) =   0.11522967656921087154811609734510D-01
    w(109) =   0.10982883090068975788799657376065D-01
    w(110) =   0.10436130863141005225673171997668D-01
    w(111) =   0.98830429087554914716648010899606D-02
    w(112) =   0.93239550065309714787536985834029D-02
    w(113) =   0.87592065795403145773316804234385D-02
    w(114) =   0.81891404887415730817235884718726D-02
    w(115) =   0.76141028256526859356393930849227D-02
    w(116) =   0.70344427036681608755685893032552D-02
    w(117) =   0.64505120486899171845442463868748D-02
    w(118) =   0.58626653903523901033648343751367D-02
    w(119) =   0.52712596565634400891303815906251D-02
    w(120) =   0.46766539777779034772638165662478D-02
    w(121) =   0.40792095178254605327114733456293D-02
    w(122) =   0.34792893810051465908910894094105D-02
    w(123) =   0.28772587656289004082883197417581D-02
    w(124) =   0.22734860707492547802810838362671D-02
    w(125) =   0.16683488125171936761028811985672D-02
    w(126) =   0.10622766869538486959954760554099D-02
    w(127) =   0.45645726109586654495731936146574D-03
 
  else if ( n == 255 ) then

    x(  1) =  -0.9999557053175637D+00
    x(  2) =  -0.9997666213120006D+00
    x(  3) =  -0.9994264746801700D+00
    x(  4) =  -0.9989352412846546D+00
    x(  5) =  -0.9982929861369679D+00
    x(  6) =  -0.9974998041266158D+00
    x(  7) =  -0.9965558144351986D+00
    x(  8) =  -0.9954611594800263D+00
    x(  9) =  -0.9942160046166302D+00
    x( 10) =  -0.9928205380219891D+00
    x( 11) =  -0.9912749706303856D+00
    x( 12) =  -0.9895795360859201D+00
    x( 13) =  -0.9877344906997324D+00
    x( 14) =  -0.9857401134074193D+00
    x( 15) =  -0.9835967057247763D+00
    x( 16) =  -0.9813045917010171D+00
    x( 17) =  -0.9788641178690681D+00
    x( 18) =  -0.9762756531927360D+00
    x( 19) =  -0.9735395890106436D+00
    x( 20) =  -0.9706563389768804D+00
    x( 21) =  -0.9676263389983388D+00
    x( 22) =  -0.9644500471687263D+00
    x( 23) =  -0.9611279436992478D+00
    x( 24) =  -0.9576605308459620D+00
    x( 25) =  -0.9540483328338163D+00
    x( 26) =  -0.9502918957773683D+00
    x( 27) =  -0.9463917875982043D+00
    x( 28) =  -0.9423485979390644D+00
    x( 29) =  -0.9381629380746873D+00
    x( 30) =  -0.9338354408193861D+00
    x( 31) =  -0.9293667604313699D+00
    x( 32) =  -0.9247575725138244D+00
    x( 33) =  -0.9200085739127664D+00
    x( 34) =  -0.9151204826116870D+00
    x( 35) =  -0.9100940376230008D+00
    x( 36) =  -0.9049299988763150D+00
    x( 37) =  -0.8996291471035368D+00
    x( 38) =  -0.8941922837208367D+00
    x( 39) =  -0.8886202307074841D+00
    x( 40) =  -0.8829138304815741D+00
    x( 41) =  -0.8770739457726654D+00
    x( 42) =  -0.8711014594913465D+00
    x( 43) =  -0.8649972745957512D+00
    x( 44) =  -0.8587623139550430D+00
    x( 45) =  -0.8523975202098902D+00
    x( 46) =  -0.8459038556299511D+00
    x( 47) =  -0.8392823019683910D+00
    x( 48) =  -0.8325338603134556D+00
    x( 49) =  -0.8256595509371186D+00
    x( 50) =  -0.8186604131408319D+00
    x( 51) =  -0.8115375050983958D+00
    x( 52) =  -0.8042919036959787D+00
    x( 53) =  -0.7969247043693057D+00
    x( 54) =  -0.7894370209380444D+00
    x( 55) =  -0.7818299854374094D+00
    x( 56) =  -0.7741047479470157D+00
    x( 57) =  -0.7662624764170006D+00
    x( 58) =  -0.7583043564914468D+00
    x( 59) =  -0.7502315913291283D+00
    x( 60) =  -0.7420454014216102D+00
    x( 61) =  -0.7337470244087263D+00
    x( 62) =  -0.7253377148914649D+00
    x( 63) =  -0.7168187442422908D+00
    x( 64) =  -0.7081914004129306D+00
    x( 65) =  -0.6994569877396524D+00
    x( 66) =  -0.6906168267460676D+00
    x( 67) =  -0.6816722539434864D+00
    x( 68) =  -0.6726246216288551D+00
    x( 69) =  -0.6634752976803070D+00
    x( 70) =  -0.6542256653503588D+00
    x( 71) =  -0.6448771230567811D+00
    x( 72) =  -0.6354310841711771D+00
    x( 73) =  -0.6258889768052999D+00
    x( 74) =  -0.6162522435951415D+00
    x( 75) =  -0.6065223414828266D+00
    x( 76) =  -0.5967007414963417D+00
    x( 77) =  -0.5867889285271373D+00
    x( 78) =  -0.5767884011056313D+00
    x( 79) =  -0.5667006711746527D+00
    x( 80) =  -0.5565272638608558D+00
    x( 81) =  -0.5462697172441424D+00
    x( 82) =  -0.5359295821251249D+00
    x( 83) =  -0.5255084217906666D+00
    x( 84) =  -0.5150078117775342D+00
    x( 85) =  -0.5044293396341982D+00
    x( 86) =  -0.4937746046808170D+00
    x( 87) =  -0.4830452177674420D+00
    x( 88) =  -0.4722428010304787D+00
    x( 89) =  -0.4613689876474424D+00
    x( 90) =  -0.4504254215900437D+00
    x( 91) =  -0.4394137573756426D+00
    x( 92) =  -0.4283356598171081D+00
    x( 93) =  -0.4171928037711214D+00
    x( 94) =  -0.4059868738849605D+00
    x( 95) =  -0.3947195643418044D+00
    x( 96) =  -0.3833925786045958D+00
    x( 97) =  -0.3720076291585012D+00
    x( 98) =  -0.3605664372520062D+00
    x( 99) =  -0.3490707326366864D+00
    x(100) =  -0.3375222533056927D+00
    x(101) =  -0.3259227452309905D+00
    x(102) =  -0.3142739620993925D+00
    x(103) =  -0.3025776650474256D+00
    x(104) =  -0.2908356223950708D+00
    x(105) =  -0.2790496093784178D+00
    x(106) =  -0.2672214078812731D+00
    x(107) =  -0.2553528061657641D+00
    x(108) =  -0.2434455986019780D+00
    x(109) =  -0.2315015853966777D+00
    x(110) =  -0.2195225723211354D+00
    x(111) =  -0.2075103704381242D+00
    x(112) =  -0.1954667958281108D+00
    x(113) =  -0.1833936693146885D+00
    x(114) =  -0.1712928161892939D+00
    x(115) =  -0.1591660659352477D+00
    x(116) =  -0.1470152519511620D+00
    x(117) =  -0.1348422112737553D+00
    x(118) =  -0.1226487843001178D+00
    x(119) =  -0.1104368145094688D+00
    x(120) =  -0.9820814818444755D-01
    x(121) =  -0.8596463413198061D-01
    x(122) =  -0.7370812340376778D-01
    x(123) =  -0.6144046901642827D-01
    x(124) =  -0.4916352567134998D-01
    x(125) =  -0.3687914947428402D-01
    x(126) =  -0.2458919765472701D-01
    x(127) =  -0.1229552828513332D-01
    x(128) =    0.000000000000000D+00
    x(129) =   0.1229552828513332D-01
    x(130) =   0.2458919765472701D-01
    x(131) =   0.3687914947428402D-01
    x(132) =   0.4916352567134998D-01
    x(133) =   0.6144046901642827D-01
    x(134) =   0.7370812340376778D-01
    x(135) =   0.8596463413198061D-01
    x(136) =   0.9820814818444755D-01
    x(137) =   0.1104368145094688D+00
    x(138) =   0.1226487843001178D+00
    x(139) =   0.1348422112737553D+00
    x(140) =   0.1470152519511620D+00
    x(141) =   0.1591660659352477D+00
    x(142) =   0.1712928161892939D+00
    x(143) =   0.1833936693146885D+00
    x(144) =   0.1954667958281108D+00
    x(145) =   0.2075103704381242D+00
    x(146) =   0.2195225723211354D+00
    x(147) =   0.2315015853966777D+00
    x(148) =   0.2434455986019780D+00
    x(149) =   0.2553528061657641D+00
    x(150) =   0.2672214078812731D+00
    x(151) =   0.2790496093784178D+00
    x(152) =   0.2908356223950708D+00
    x(153) =   0.3025776650474256D+00
    x(154) =   0.3142739620993925D+00
    x(155) =   0.3259227452309905D+00
    x(156) =   0.3375222533056927D+00
    x(157) =   0.3490707326366864D+00
    x(158) =   0.3605664372520062D+00
    x(159) =   0.3720076291585012D+00
    x(160) =   0.3833925786045958D+00
    x(161) =   0.3947195643418044D+00
    x(162) =   0.4059868738849605D+00
    x(163) =   0.4171928037711214D+00
    x(164) =   0.4283356598171081D+00
    x(165) =   0.4394137573756426D+00
    x(166) =   0.4504254215900437D+00
    x(167) =   0.4613689876474424D+00
    x(168) =   0.4722428010304787D+00
    x(169) =   0.4830452177674420D+00
    x(170) =   0.4937746046808170D+00
    x(171) =   0.5044293396341982D+00
    x(172) =   0.5150078117775342D+00
    x(173) =   0.5255084217906666D+00
    x(174) =   0.5359295821251249D+00
    x(175) =   0.5462697172441424D+00
    x(176) =   0.5565272638608558D+00
    x(177) =   0.5667006711746527D+00
    x(178) =   0.5767884011056313D+00
    x(179) =   0.5867889285271373D+00
    x(180) =   0.5967007414963417D+00
    x(181) =   0.6065223414828266D+00
    x(182) =   0.6162522435951415D+00
    x(183) =   0.6258889768052999D+00
    x(184) =   0.6354310841711771D+00
    x(185) =   0.6448771230567811D+00
    x(186) =   0.6542256653503588D+00
    x(187) =   0.6634752976803070D+00
    x(188) =   0.6726246216288551D+00
    x(189) =   0.6816722539434864D+00
    x(190) =   0.6906168267460676D+00
    x(191) =   0.6994569877396524D+00
    x(192) =   0.7081914004129306D+00
    x(193) =   0.7168187442422908D+00
    x(194) =   0.7253377148914649D+00
    x(195) =   0.7337470244087263D+00
    x(196) =   0.7420454014216102D+00
    x(197) =   0.7502315913291283D+00
    x(198) =   0.7583043564914468D+00
    x(199) =   0.7662624764170006D+00
    x(200) =   0.7741047479470157D+00
    x(201) =   0.7818299854374094D+00
    x(202) =   0.7894370209380444D+00
    x(203) =   0.7969247043693057D+00
    x(204) =   0.8042919036959787D+00
    x(205) =   0.8115375050983958D+00
    x(206) =   0.8186604131408319D+00
    x(207) =   0.8256595509371186D+00
    x(208) =   0.8325338603134556D+00
    x(209) =   0.8392823019683910D+00
    x(210) =   0.8459038556299511D+00
    x(211) =   0.8523975202098902D+00
    x(212) =   0.8587623139550430D+00
    x(213) =   0.8649972745957512D+00
    x(214) =   0.8711014594913465D+00
    x(215) =   0.8770739457726654D+00
    x(216) =   0.8829138304815741D+00
    x(217) =   0.8886202307074841D+00
    x(218) =   0.8941922837208367D+00
    x(219) =   0.8996291471035368D+00
    x(220) =   0.9049299988763150D+00
    x(221) =   0.9100940376230008D+00
    x(222) =   0.9151204826116870D+00
    x(223) =   0.9200085739127664D+00
    x(224) =   0.9247575725138244D+00
    x(225) =   0.9293667604313699D+00
    x(226) =   0.9338354408193861D+00
    x(227) =   0.9381629380746873D+00
    x(228) =   0.9423485979390644D+00
    x(229) =   0.9463917875982043D+00
    x(230) =   0.9502918957773683D+00
    x(231) =   0.9540483328338163D+00
    x(232) =   0.9576605308459620D+00
    x(233) =   0.9611279436992478D+00
    x(234) =   0.9644500471687263D+00
    x(235) =   0.9676263389983388D+00
    x(236) =   0.9706563389768804D+00
    x(237) =   0.9735395890106436D+00
    x(238) =   0.9762756531927360D+00
    x(239) =   0.9788641178690681D+00
    x(240) =   0.9813045917010171D+00
    x(241) =   0.9835967057247763D+00
    x(242) =   0.9857401134074193D+00
    x(243) =   0.9877344906997324D+00
    x(244) =   0.9895795360859201D+00
    x(245) =   0.9912749706303856D+00
    x(246) =   0.9928205380219891D+00
    x(247) =   0.9942160046166302D+00
    x(248) =   0.9954611594800263D+00
    x(249) =   0.9965558144351986D+00
    x(250) =   0.9974998041266158D+00
    x(251) =   0.9982929861369679D+00
    x(252) =   0.9989352412846546D+00
    x(253) =   0.9994264746801700D+00
    x(254) =   0.9997666213120006D+00
    x(255) =   0.9999557053175637D+00

    w(  1) =   0.1136736199914808D-03
    w(  2) =   0.2645938711908564D-03
    w(  3) =   0.4156976252681932D-03
    w(  4) =   0.5667579456482639D-03
    w(  5) =   0.7177364780061286D-03
    w(  6) =   0.8686076661194581D-03
    w(  7) =   0.1019347976427318D-02
    w(  8) =   0.1169934372938800D-02
    w(  9) =   0.1320343990022177D-02
    w( 10) =   0.1470554042778403D-02
    w( 11) =   0.1620541799041545D-02
    w( 12) =   0.1770284570660304D-02
    w( 13) =   0.1919759711713187D-02
    w( 14) =   0.2068944619501569D-02
    w( 15) =   0.2217816736754017D-02
    w( 16) =   0.2366353554396287D-02
    w( 17) =   0.2514532614599710D-02
    w( 18) =   0.2662331513971696D-02
    w( 19) =   0.2809727906820460D-02
    w( 20) =   0.2956699508457498D-02
    w( 21) =   0.3103224098519095D-02
    w( 22) =   0.3249279524294296D-02
    w( 23) =   0.3394843704053401D-02
    w( 24) =   0.3539894630372244D-02
    w( 25) =   0.3684410373449933D-02
    w( 26) =   0.3828369084417135D-02
    w( 27) =   0.3971748998634907D-02
    w( 28) =   0.4114528438981242D-02
    w( 29) =   0.4256685819126112D-02
    w( 30) =   0.4398199646792759D-02
    w( 31) =   0.4539048527006180D-02
    w( 32) =   0.4679211165326077D-02
    w( 33) =   0.4818666371065699D-02
    w( 34) =   0.4957393060495050D-02
    w( 35) =   0.5095370260027839D-02
    w( 36) =   0.5232577109391968D-02
    w( 37) =   0.5368992864783177D-02
    w( 38) =   0.5504596902000804D-02
    w( 39) =   0.5639368719565862D-02
    w( 40) =   0.5773287941820301D-02
    w( 41) =   0.5906334322007422D-02
    w( 42) =   0.6038487745332765D-02
    w( 43) =   0.6169728232005295D-02
    w( 44) =   0.6300035940257733D-02
    w( 45) =   0.6429391169346602D-02
    w( 46) =   0.6557774362530328D-02
    w( 47) =   0.6685166110026254D-02
    w( 48) =   0.6811547151944815D-02
    w( 49) =   0.6936898381201466D-02
    w( 50) =   0.7061200846405536D-02
    w( 51) =   0.7184435754724984D-02
    w( 52) =   0.7306584474728122D-02
    w( 53) =   0.7427628539199977D-02
    w( 54) =   0.7547549647934514D-02
    w( 55) =   0.7666329670501377D-02
    w( 56) =   0.7783950648986801D-02
    w( 57) =   0.7900394800708624D-02
    w( 58) =   0.8015644520904983D-02
    w( 59) =   0.8129682385395602D-02
    w( 60) =   0.8242491153216323D-02
    w( 61) =   0.8354053769225508D-02
    w( 62) =   0.8464353366682819D-02
    w( 63) =   0.8573373269798925D-02
    w( 64) =   0.8681096996256795D-02
    w( 65) =   0.8787508259703609D-02
    w( 66) =   0.8892590972213036D-02
    w( 67) =   0.8996329246717397D-02
    w( 68) =   0.9098707399409718D-02
    w( 69) =   0.9199709952114802D-02
    w( 70) =   0.9299321634629343D-02
    w( 71) =   0.9397527387030594D-02
    w( 72) =   0.9494312361953241D-02
    w( 73) =   0.9589661926834022D-02
    w( 74) =   0.9683561666124043D-02
    w( 75) =   0.9775997383468165D-02
    w( 76) =   0.9866955103851452D-02
    w( 77) =   0.9956421075711706D-02
    w( 78) =   0.1004438177301882D-01
    w( 79) =   0.1013082389731963D-01
    w( 80) =   0.1021573437974821D-01
    w( 81) =   0.1029910038300220D-01
    w( 82) =   0.1038090930328312D-01
    w( 83) =   0.1046114877220228D-01
    w( 84) =   0.1053980665865038D-01
    w( 85) =   0.1061687107063194D-01
    w( 86) =   0.1069233035706287D-01
    w( 87) =   0.1076617310953212D-01
    w( 88) =   0.1083838816402652D-01
    w( 89) =   0.1090896460261843D-01
    w( 90) =   0.1097789175511656D-01
    w( 91) =   0.1104515920067912D-01
    w( 92) =   0.1111075676938929D-01
    w( 93) =   0.1117467454379268D-01
    w( 94) =   0.1123690286039691D-01
    w( 95) =   0.1129743231113249D-01
    w( 96) =   0.1135625374477508D-01
    w( 97) =   0.1141335826832922D-01
    w( 98) =   0.1146873724837283D-01
    w( 99) =   0.1152238231236217D-01
    w(100) =   0.1157428534989815D-01
    w(101) =   0.1162443851395193D-01
    w(102) =   0.1167283422205182D-01
    w(103) =   0.1171946515742932D-01
    w(104) =   0.1176432427012535D-01
    w(105) =   0.1180740477805627D-01
    w(106) =   0.1184870016803913D-01
    w(107) =   0.1188820419677619D-01
    w(108) =   0.1192591089179929D-01
    w(109) =   0.1196181455237226D-01
    w(110) =   0.1199590975035326D-01
    w(111) =   0.1202819133101508D-01
    w(112) =   0.1205865441382472D-01
    w(113) =   0.1208729439318107D-01
    w(114) =   0.1211410693911137D-01
    w(115) =   0.1213908799792579D-01
    w(116) =   0.1216223379283022D-01
    w(117) =   0.1218354082449738D-01
    w(118) =   0.1220300587159574D-01
    w(119) =   0.1222062599127671D-01
    w(120) =   0.1223639851961942D-01
    w(121) =   0.1225032107203351D-01
    w(122) =   0.1226239154361966D-01
    w(123) =   0.1227260810948789D-01
    w(124) =   0.1228096922503318D-01
    w(125) =   0.1228747362616942D-01
    w(126) =   0.1229212032952021D-01
    w(127) =   0.1229490863256759D-01
    w(128) =   0.1229583811375833D-01
    w(129) =   0.1229490863256759D-01
    w(130) =   0.1229212032952021D-01
    w(131) =   0.1228747362616942D-01
    w(132) =   0.1228096922503318D-01
    w(133) =   0.1227260810948789D-01
    w(134) =   0.1226239154361966D-01
    w(135) =   0.1225032107203351D-01
    w(136) =   0.1223639851961942D-01
    w(137) =   0.1222062599127671D-01
    w(138) =   0.1220300587159574D-01
    w(139) =   0.1218354082449738D-01
    w(140) =   0.1216223379283022D-01
    w(141) =   0.1213908799792579D-01
    w(142) =   0.1211410693911137D-01
    w(143) =   0.1208729439318107D-01
    w(144) =   0.1205865441382472D-01
    w(145) =   0.1202819133101508D-01
    w(146) =   0.1199590975035326D-01
    w(147) =   0.1196181455237226D-01
    w(148) =   0.1192591089179929D-01
    w(149) =   0.1188820419677619D-01
    w(150) =   0.1184870016803913D-01
    w(151) =   0.1180740477805627D-01
    w(152) =   0.1176432427012535D-01
    w(153) =   0.1171946515742932D-01
    w(154) =   0.1167283422205182D-01
    w(155) =   0.1162443851395193D-01
    w(156) =   0.1157428534989815D-01
    w(157) =   0.1152238231236217D-01
    w(158) =   0.1146873724837283D-01
    w(159) =   0.1141335826832922D-01
    w(160) =   0.1135625374477508D-01
    w(161) =   0.1129743231113249D-01
    w(162) =   0.1123690286039691D-01
    w(163) =   0.1117467454379268D-01
    w(164) =   0.1111075676938929D-01
    w(165) =   0.1104515920067912D-01
    w(166) =   0.1097789175511656D-01
    w(167) =   0.1090896460261843D-01
    w(168) =   0.1083838816402652D-01
    w(169) =   0.1076617310953212D-01
    w(170) =   0.1069233035706287D-01
    w(171) =   0.1061687107063194D-01
    w(172) =   0.1053980665865038D-01
    w(173) =   0.1046114877220228D-01
    w(174) =   0.1038090930328312D-01
    w(175) =   0.1029910038300220D-01
    w(176) =   0.1021573437974821D-01
    w(177) =   0.1013082389731963D-01
    w(178) =   0.1004438177301882D-01
    w(179) =   0.9956421075711706D-02
    w(180) =   0.9866955103851452D-02
    w(181) =   0.9775997383468165D-02
    w(182) =   0.9683561666124043D-02
    w(183) =   0.9589661926834022D-02
    w(184) =   0.9494312361953241D-02
    w(185) =   0.9397527387030594D-02
    w(186) =   0.9299321634629343D-02
    w(187) =   0.9199709952114802D-02
    w(188) =   0.9098707399409718D-02
    w(189) =   0.8996329246717397D-02
    w(190) =   0.8892590972213036D-02
    w(191) =   0.8787508259703609D-02
    w(192) =   0.8681096996256795D-02
    w(193) =   0.8573373269798925D-02
    w(194) =   0.8464353366682819D-02
    w(195) =   0.8354053769225508D-02
    w(196) =   0.8242491153216323D-02
    w(197) =   0.8129682385395602D-02
    w(198) =   0.8015644520904983D-02
    w(199) =   0.7900394800708624D-02
    w(200) =   0.7783950648986801D-02
    w(201) =   0.7666329670501377D-02
    w(202) =   0.7547549647934514D-02
    w(203) =   0.7427628539199977D-02
    w(204) =   0.7306584474728122D-02
    w(205) =   0.7184435754724984D-02
    w(206) =   0.7061200846405536D-02
    w(207) =   0.6936898381201466D-02
    w(208) =   0.6811547151944815D-02
    w(209) =   0.6685166110026254D-02
    w(210) =   0.6557774362530328D-02
    w(211) =   0.6429391169346602D-02
    w(212) =   0.6300035940257733D-02
    w(213) =   0.6169728232005295D-02
    w(214) =   0.6038487745332765D-02
    w(215) =   0.5906334322007422D-02
    w(216) =   0.5773287941820301D-02
    w(217) =   0.5639368719565862D-02
    w(218) =   0.5504596902000804D-02
    w(219) =   0.5368992864783177D-02
    w(220) =   0.5232577109391968D-02
    w(221) =   0.5095370260027839D-02
    w(222) =   0.4957393060495050D-02
    w(223) =   0.4818666371065699D-02
    w(224) =   0.4679211165326077D-02
    w(225) =   0.4539048527006180D-02
    w(226) =   0.4398199646792759D-02
    w(227) =   0.4256685819126112D-02
    w(228) =   0.4114528438981242D-02
    w(229) =   0.3971748998634907D-02
    w(230) =   0.3828369084417135D-02
    w(231) =   0.3684410373449933D-02
    w(232) =   0.3539894630372244D-02
    w(233) =   0.3394843704053401D-02
    w(234) =   0.3249279524294296D-02
    w(235) =   0.3103224098519095D-02
    w(236) =   0.2956699508457498D-02
    w(237) =   0.2809727906820460D-02
    w(238) =   0.2662331513971696D-02
    w(239) =   0.2514532614599710D-02
    w(240) =   0.2366353554396287D-02
    w(241) =   0.2217816736754017D-02
    w(242) =   0.2068944619501569D-02
    w(243) =   0.1919759711713187D-02
    w(244) =   0.1770284570660304D-02
    w(245) =   0.1620541799041545D-02
    w(246) =   0.1470554042778403D-02
    w(247) =   0.1320343990022177D-02
    w(248) =   0.1169934372938800D-02
    w(249) =   0.1019347976427318D-02
    w(250) =   0.8686076661194581D-03
    w(251) =   0.7177364780061286D-03
    w(252) =   0.5667579456482639D-03
    w(253) =   0.4156976252681932D-03
    w(254) =   0.2645938711908564D-03
    w(255) =   0.1136736199914808D-03

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LEGENDRE_SET - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of N = ', n
    write ( *, '(a)' ) &
      '  Legal values are 1 through 33, 63, 64, 65, 127 and 255.'
    stop

  end if

  return
end
subroutine map ( code, element_order, w )

!*****************************************************************************80
!
!! MAP returns the interpolation matrix for any available element.
!
!  Formula:
!
!    For an element of order ELEMENT_ORDER, we suppose we are given 
!    ELEMENT_ORDER items of data Q associated with the nodes.
!
!   Let PHI(J)(R,S) be the Lagrange basis polynomial associated with 
!   node J.  PHI(J)(R,S) is 1 at node J, and 0 at each of the other nodes.
!
!   Let P(R,S) be the polynomial of ELEMENT_ORDER terms which interpolates the
!   data Q, that is,
!
!      P(R(J),S(J)) = Q(J)
!
!   where the coordinates of node J are (R(J),S(J)).  Then we know
!   that we can write
!
!     P(R,S) = sum ( 1 <= J <= ELEMENT_ORDER ) Q(J) * PHI(J)(R,S)
!
!   But P(R,S) also has a standard representation as
!
!     P(R,S) = sum ( 1 <= I <= ELEMENT_ORDER ) A(I) * R**REXP(I) * S**SEXP(I)
!
!   where REXP(I) and SEXP(I) are the exponents of R and S and
!   the A(I) are the appropriate coefficients.
!
!   The interpolation matrix W allows us to immediately compute
!   the standard basis coefficients A from the data Q to be interpolated
!   using the formula:
!
!      A(I) = sum ( 1 <= J <= ELEMENT_ORDER ) W(I,J) * Q(J)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) CODE, identifies the element.
!    Legal values include 'Q4', 'Q8', 'Q9', 'Q12', 'Q16', 'QL', 
!    'T3', 'T4', 'T6' and 'T10'.
!
!    Input, integer ( kind = 4 ) ELEMENT_ORDER, the number of nodes per element.
!
!    Output, real ( kind = 8 ) W(ELEMENT_ORDER,ELEMENT_ORDER),
!     the interpolation matrix.
!
  implicit none

  integer   ( kind = 4 ) element_order

  real ( kind = 8 ) area
  character ( len = * )  code
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) info
  integer   ( kind = 4 ) j
  integer   ( kind = 4 ) pivot(element_order)
  real ( kind = 8 ) r(element_order)
  integer   ( kind = 4 ) rexp(element_order)
  real ( kind = 8 ) rfact
  real ( kind = 8 ) s(element_order)
  integer   ( kind = 4 ) sexp(element_order)
  real ( kind = 8 ) sfact
  real ( kind = 8 ) w(element_order,element_order)
!
!  Get the (R,S) location of the nodes.
!
  call node_reference ( code, r, s, area )
!
!  Get the associated monomials.
!
  call poly ( code, rexp, sexp )
!
!  Set up the Vandermonde matrix.
!  Factors of the form 0**0 are to be understood as 1.
!
  do i = 1, element_order
    do j = 1, element_order

      if ( rexp(j) == 0 ) then
        rfact = 1.0D+00
      else
        rfact = r(i)**rexp(j)
      end if

      if ( sexp(j) == 0 ) then
        sfact = 1.0D+00
      else
        sfact = s(i)**sexp(j)
      end if

      w(i,j) = rfact * sfact

    end do
  end do
!
!  Factor the Vandermonde matrix.
!
  call r8ge_fa ( element_order, w, pivot, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MAP - Fatal error!'
    write ( *, '(a)' ) '  The Vandermonde matrix is singular.'
    stop
  end if
!
!  Invert the Vandermonde matrix.
!
  call r8ge_inverse ( element_order, w, pivot )

  return
end
subroutine map_test ( code )

!*****************************************************************************80
!
!! MAP_TEST tests the map routines.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) CODE, the code for the element.
!    Legal values include 'Q4', 'Q8', 'Q9', 'Q12', 'Q16', 'QL', 
!    'T3', 'T4', 'T6' and 'T10'.
!
  implicit none

  character ( len = * )  code
  integer   ( kind = 4 ) element_order
  integer   ( kind = 4 ) order_code
  logical                s_eqi
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: w

  if ( s_eqi ( code, 'T4' ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MAP_TEST - Warning!'
    write ( *, '(a)' ) '  Skipping test for element "' &
      // trim ( code ) // '".'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  MAP_TEST: interpolation matrix for element "' &
    // trim ( code ) // '".'
  write ( *, '(a)' ) ' '

  element_order = order_code ( code )

  allocate ( w(1:element_order,1:element_order) )

  call map ( code, element_order, w )

  call r8mat_print ( element_order, element_order, w, &
    '  The interpolation matrix:' );

  deallocate ( w )

  return
end
subroutine mass_matrix_t6 ( node_num, element_num, element_node, node_xy, a )

!*****************************************************************************80
!
!! MASS_MATRIX_T6 computes the mass matrix, using 6-node triangles.
!
!  Discussion:
!
!    The mass matrix to be estimated has the form:
!
!      A(I,J) = integral ( PHI(I)(X,Y) * PHI(J)(X,Y) ) d Region
!
!    where PHI(I) and PHI(J) are the shape functions associated with
!    the I-th and J-th variables.
!
!  Reference Element T6:
!
!    |
!    1  3
!    |  |\
!    |  | \
!    S  6  5
!    |  |   \
!    |  |    \
!    0  1--4--2
!    |
!    +--0--R--1-->
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
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes per element.
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of elements.
!
!    Input, integer ( kind = 4 ) ELEMENT_NODE(6,ELEMENT_NUM), the nodes that 
!    make up each element.
!
!    Input, real ( kind = 8 ) NODE_XY(2,NODE_NUM), the coordinates 
!    of the nodes.
!
!    Output, real ( kind = 8 ) A(NODE_NUM,NODE_NUM), the mass matrix.
!
  implicit none

  integer ( kind = 4 ) element_num
  integer ( kind = 4 ), parameter :: element_order = 6
  integer ( kind = 4 ) node_num

  real    ( kind = 8 ) a(node_num,node_num)
  real    ( kind = 8 ) area
  real    ( kind = 8 ) dwdr(element_order)
  real    ( kind = 8 ) dwds(element_order)
  integer ( kind = 4 ) element
  integer ( kind = 4 ) element_node(element_order,element_num)
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) iq
  integer ( kind = 4 ) jp
  integer ( kind = 4 ) jq
  real    ( kind = 8 ) node_xy(2,node_num)
  integer ( kind = 4 ) p1
  integer ( kind = 4 ) p2
  integer ( kind = 4 ) p3
  integer ( kind = 4 ) quad
  integer ( kind = 4 ) quad_num
  real    ( kind = 8 ) r
  real    ( kind = 8 ), allocatable, dimension ( : ) :: rtab
  integer ( kind = 4 ) rule
  real    ( kind = 8 ) s
  real    ( kind = 8 ), allocatable, dimension ( : ) :: stab
  integer ( kind = 4 ) triangle_unit_size
  real    ( kind = 8 ) w(element_order)
  real    ( kind = 8 ), allocatable, dimension ( : ) :: weight
!
!  Zero out the matrix.
!
  a(1:node_num,1:node_num) = 0.0D+00
!
!  Get the weights and abscissas for a unit triangle.
!
  rule = 12
  quad_num = triangle_unit_size ( rule )

  allocate ( rtab(1:quad_num) )
  allocate ( stab(1:quad_num) )
  allocate ( weight(1:quad_num) )

  call triangle_unit_set ( rule, rtab, stab, weight )
!
!  For each element.
!
  do element = 1, element_num

    p1 = element_node(1,element)
    p2 = element_node(2,element)
    p3 = element_node(3,element)

    area = 0.5D+00 * abs ( &
        node_xy(1,p1) * ( node_xy(2,p2) - node_xy(2,p3) ) &
      + node_xy(1,p2) * ( node_xy(2,p3) - node_xy(2,p1) ) &
      + node_xy(1,p3) * ( node_xy(2,p1) - node_xy(2,p2) ) )

    if ( area == 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MASS_MATRIX_T6 - Fatal error!'
      write ( *, '(a,i8)' ) '  Zero area for element ', element
      stop
    end if
!
!  For each quadrature point in the element...
!
    do quad = 1, quad_num

      r = rtab(quad)
      s = stab(quad)

      call shape_t6 ( r, s, w, dwdr, dwds )
!
!  For each basis function PHI(I) associated with a node in the element,
!
      do iq = 1, element_order

        ip = element_node(iq,element)
!
!  For each "neighbor" basis function PHI(J) associated with a node in
!  the element.
!
        do jq = 1, element_order

          jp = element_node(jq,element)

          a(ip,jp) = a(ip,jp) + area * weight(quad) * w(iq) * w(jq)

        end do
      end do
    end do
  end do

  deallocate ( rtab )
  deallocate ( stab )
  deallocate ( weight )

  return
end
function next_boundary_node ( node, code )

!*****************************************************************************80
!
!! NEXT_BOUNDARY_NODE returns the next boundary node in any element.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE, the index of the current node.  An input
!    value of 0 (or any "unusual" value") indicates that the
!    first edge node is desired.
!
!    Output, integer ( kind = 4 ) NEXT_BOUNDARY_NODE, the index of the next 
!    edge node.
!
  implicit none

  character ( len = * )  code
  integer   ( kind = 4 ) next_boundary_node
  integer   ( kind = 4 ) next_boundary_node_q4
  integer   ( kind = 4 ) next_boundary_node_q8
  integer   ( kind = 4 ) next_boundary_node_q9
  integer   ( kind = 4 ) next_boundary_node_q12
  integer   ( kind = 4 ) next_boundary_node_q16
  integer   ( kind = 4 ) next_boundary_node_ql
  integer   ( kind = 4 ) next_boundary_node_t3
  integer   ( kind = 4 ) next_boundary_node_t4
  integer   ( kind = 4 ) next_boundary_node_t6
  integer   ( kind = 4 ) next_boundary_node_t10
  integer   ( kind = 4 ) node
  logical                s_eqi

  if ( s_eqi ( code, 'Q4' ) ) then
    next_boundary_node = next_boundary_node_q4 ( node )
  else if ( s_eqi ( code, 'Q8' ) ) then
    next_boundary_node = next_boundary_node_q8 ( node )
  else if ( s_eqi ( code, 'Q9' ) ) then
    next_boundary_node = next_boundary_node_q9 ( node )
  else if ( s_eqi ( code, 'Q12' ) ) then
    next_boundary_node = next_boundary_node_q12 ( node )
  else if ( s_eqi ( code, 'Q16' ) ) then
    next_boundary_node = next_boundary_node_q16 ( node )
  else if ( s_eqi ( code, 'QL' ) ) then
    next_boundary_node = next_boundary_node_ql ( node )
  else if ( s_eqi ( code, 'T3' ) ) then
    next_boundary_node = next_boundary_node_t3 ( node )
  else if ( s_eqi ( code, 'T4' ) ) then
    next_boundary_node = next_boundary_node_t4 ( node )
  else if ( s_eqi ( code, 'T6' ) ) then
    next_boundary_node = next_boundary_node_t6 ( node )
  else if ( s_eqi ( code, 'T10' ) ) then
    next_boundary_node = next_boundary_node_t10 ( node )
  else
    next_boundary_node = -1
  end if

  return
end
function next_boundary_node_q4 ( node )

!*****************************************************************************80
!
!! NEXT_BOUNDARY_NODE_Q4 returns the next boundary node in a Q4 element.
!
!  Reference Element Q4:
!
!    |
!    1  4-----3
!    |  |     |
!    |  |     |
!    S  |     |
!    |  |     |
!    |  |     |
!    0  1-----2
!    |
!    +--0--R--1-->
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 April 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE, the index of the current node.  An input
!    value of 0 (or any "unusual" value") indicates that the
!    first edge node is desired.
!
!    Output, integer ( kind = 4 ) NEXT_BOUNDARY_NODE_Q4, the index of the next 
!    edge node.
!
  implicit none

  integer ( kind = 4 ) next_boundary_node_q4
  integer ( kind = 4 ) node

  if ( node == 1 ) then
    next_boundary_node_q4 = 2
  else if ( node == 2 ) then
    next_boundary_node_q4 = 3
  else if ( node == 3 ) then
    next_boundary_node_q4 = 4
  else if ( node == 4 ) then
    next_boundary_node_q4 = 1
  else
    next_boundary_node_q4 = -1
  end if

  return
end
function next_boundary_node_q8 ( node )

!*****************************************************************************80
!
!! NEXT_BOUNDARY_NODE_Q8 returns the next boundary node in a Q8 element.
!
!  Reference Element Q8:
!
!    |
!    1  4--7--3
!    |  |     |
!    |  |     |
!    S  8     6
!    |  |     |
!    |  |     |
!    0  1--5--2
!    |
!    +--0--R--1-->
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 April 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE, the index of the current node.  An input
!    value of 0 (or any "unusual" value") indicates that the
!    first edge node is desired.
!
!    Output, integer ( kind = 4 ) NEXT_BOUNDARY_NODE_Q8, the index of the next 
!    edge node.
!
  implicit none

  integer ( kind = 4 ) next_boundary_node_q8
  integer ( kind = 4 ) node

  if ( node == 1 ) then
    next_boundary_node_q8 = 5
  else if ( node == 5 ) then
    next_boundary_node_q8 = 2
  else if ( node == 2 ) then
    next_boundary_node_q8 = 6
  else if ( node == 6 ) then
    next_boundary_node_q8 = 3
  else if ( node == 3 ) then
    next_boundary_node_q8 = 7
  else if ( node == 7 ) then
    next_boundary_node_q8 = 4
  else if ( node == 4 ) then
    next_boundary_node_q8 = 8
  else if ( node == 8 ) then
    next_boundary_node_q8 = 1
  else
    next_boundary_node_q8 = -1
  end if

  return
end
function next_boundary_node_q9 ( node )

!*****************************************************************************80
!
!! NEXT_BOUNDARY_NODE_Q9 returns the next boundary node in a Q9 element.
!
!  Reference Element Q9:
!
!    |
!    1  4--7--3
!    |  |     |
!    |  |     |
!    S  8  9  6
!    |  |     |
!    |  |     |
!    0  1--5--2
!    |
!    +--0--R--1-->
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 April 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE, the index of the current node.  An input
!    value of 0 (or any "unusual" value") indicates that the
!    first edge node is desired.
!
!    Output, integer ( kind = 4 ) NEXT_BOUNDARY_NODE_Q9, the index of the next 
!    edge node.
!
  implicit none

  integer ( kind = 4 ) next_boundary_node_q9
  integer ( kind = 4 ) node

  if ( node == 1 ) then
    next_boundary_node_q9 = 5
  else if ( node == 5 ) then
    next_boundary_node_q9 = 2
  else if ( node == 2 ) then
    next_boundary_node_q9 = 6
  else if ( node == 6 ) then
    next_boundary_node_q9 = 3
  else if ( node == 3 ) then
    next_boundary_node_q9 = 7
  else if ( node == 7 ) then
    next_boundary_node_q9 = 4
  else if ( node == 4 ) then
    next_boundary_node_q9 = 8
  else if ( node == 8 ) then
    next_boundary_node_q9 = 1
  else 
    next_boundary_node_q9 = -1
  end if

  return
end
function next_boundary_node_q12 ( node )

!*****************************************************************************80
!
!! NEXT_BOUNDARY_NODE_Q12 returns the next boundary node in a Q12 element.
!
!  Reference Element Q12:
!
!    |
!    1  9-10-11-12
!    |  |        |
!    |  7        8
!    S  |        |
!    |  5        6
!    |  |        |
!    0  1--2--3--4
!    |
!    +--0---R---1-->
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 April 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE, the index of the current node.  An input
!    value of 0 (or any "unusual" value") indicates that the
!    first edge node is desired.
!
!    Output, integer ( kind = 4 ) NEXT_BOUNDARY_NODE_Q12, the index of the 
!    next edge node.
!
  implicit none

  integer ( kind = 4 ) next_boundary_node_q12
  integer ( kind = 4 ) node

  if ( node == 1 ) then
    next_boundary_node_q12 = 2
  else if ( node == 2 ) then
    next_boundary_node_q12 = 3
  else if ( node == 3 ) then
    next_boundary_node_q12 = 4
  else if ( node == 4 ) then
    next_boundary_node_q12 = 6
  else if ( node == 6 ) then
    next_boundary_node_q12 = 8
  else if ( node == 8 ) then
    next_boundary_node_q12 = 12
  else if ( node == 12 ) then
    next_boundary_node_q12 = 11
  else if ( node == 11 ) then
    next_boundary_node_q12 = 10
  else if ( node == 10 ) then
    next_boundary_node_q12 = 9
  else if ( node == 9 ) then
    next_boundary_node_q12 = 7
  else if ( node == 7 ) then
    next_boundary_node_q12 = 5
  else if ( node == 5 ) then
    next_boundary_node_q12 = 1
  else
    next_boundary_node_q12 = -1
  end if

  return
end
function next_boundary_node_q16 ( node )

!*****************************************************************************80
!
!! NEXT_BOUNDARY_NODE_Q16 returns the next boundary node in a Q16 element.
!
!  Reference Element Q16:
!
!    |
!    1 13--14--15--16
!    |  |   :   :   |
!    |  |   :   :   |
!    |  9..10..11..12
!    S  |   :   :   |
!    |  |   :   :   |
!    |  5...6...7...8
!    |  |   :   :   |
!    |  |   :   :   |  
!    0  1---2---3---4
!    |
!    +--0-----R-----1-->
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 April 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE, the index of the current node.  An input
!    value of 0 (or any "unusual" value") indicates that the
!    first edge node is desired.
!
!    Output, integer ( kind = 4 ) NEXT_BOUNDARY_NODE_Q16, the index of the next 
!    edge node.
!
  implicit none

  integer ( kind = 4 ) next_boundary_node_q16
  integer ( kind = 4 ) node

  if ( node == 1 ) then
    next_boundary_node_q16 = 2
  else if ( node == 2 ) then
    next_boundary_node_q16 = 3
  else if ( node == 3 ) then
    next_boundary_node_q16 = 4
  else if ( node == 4 ) then
    next_boundary_node_q16 = 8
  else if ( node == 8 ) then
    next_boundary_node_q16 = 12
  else if ( node == 12 ) then
    next_boundary_node_q16 = 16
  else if ( node == 16 ) then
    next_boundary_node_q16 = 15
  else if ( node == 15 ) then
    next_boundary_node_q16 = 14
  else if ( node == 14 ) then
    next_boundary_node_q16 = 13
  else if ( node == 13 ) then
    next_boundary_node_q16 = 9
  else if ( node == 9 ) then
    next_boundary_node_q16 = 5
  else if ( node == 5 ) then
    next_boundary_node_q16 = 1
  else
    next_boundary_node_q16 = -1
  end if

  return
end
function next_boundary_node_ql ( node )

!*****************************************************************************80
!
!! NEXT_BOUNDARY_NODE_QL returns the next boundary node in a QL element.
!
!  Reference Element QL:
!
!    |
!    1  4---5---6
!    |  |       |
!    |  |       |
!    S  |       |
!    |  |       |
!    |  |       |
!    0  1---2---3
!    |
!    +--0---R---1-->
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 April 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE, the index of the current node.  An input
!    value of 0 (or any "unusual" value") indicates that the
!    first edge node is desired.
!
!    Output, integer ( kind = 4 ) NEXT_BOUNDARY_NODE_QL, the index of the next 
!    edge node.
!
  implicit none

  integer ( kind = 4 ) next_boundary_node_ql
  integer ( kind = 4 ) node

  if ( node == 1 ) then
    next_boundary_node_ql = 2
  else if ( node == 2 ) then
    next_boundary_node_ql = 3
  else if ( node == 3 ) then
    next_boundary_node_ql = 6
  else if ( node == 6 ) then
    next_boundary_node_ql = 5
  else if ( node == 5 ) then
    next_boundary_node_ql = 4
  else if ( node == 4 ) then
    next_boundary_node_ql = 1
  else
    next_boundary_node_ql = -1
  end if

  return
end
function next_boundary_node_t3 ( node )

!*****************************************************************************80
!
!! NEXT_BOUNDARY_NODE_T3 returns the next boundary node in a T3 element.
!
!  Reference Element T3:
!
!    |
!    1  3
!    |  |\
!    |  | \
!    S  |  \
!    |  |   \
!    |  |    \
!    0  1-----2
!    |
!    +--0--R--1-->
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 April 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE, the index of the current node.  An input
!    value of 0 (or any "unusual" value") indicates that the
!    first edge node is desired.
!
!    Output, integer ( kind = 4 ) NEXT_BOUNDARY_NODE_T3, the index of the next 
!    edge node.
!
  implicit none

  integer ( kind = 4 ) next_boundary_node_t3
  integer ( kind = 4 ) node

  if ( node == 1 ) then
    next_boundary_node_t3 = 2
  else if ( node == 2 ) then
    next_boundary_node_t3 = 3
  else if ( node == 3 ) then
    next_boundary_node_t3 = 1
  else
    next_boundary_node_t3 = -1
  end if

  return
end
function next_boundary_node_t4 ( node )

!*****************************************************************************80
!
!! NEXT_BOUNDARY_NODE_T4 returns the next boundary node in a T4 element.
!
!  Reference Element T4:
!
!    |
!    1  3
!    |  |\
!    |  | \
!    S  |  \
!    |  | 4 \
!    |  |    \
!    0  1-----2
!    |
!    +--0--R--1-->
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE, the index of the current node.  An input
!    value of 0 (or any "unusual" value") indicates that the
!    first edge node is desired.
!
!    Output, integer ( kind = 4 ) NEXT_BOUNDARY_NODE_T4, the index of the next 
!    edge node.
!
  implicit none

  integer ( kind = 4 ) next_boundary_node_t4
  integer ( kind = 4 ) node

  if ( node == 1 ) then
    next_boundary_node_t4 = 2
  else if ( node == 2 ) then
    next_boundary_node_t4 = 3
  else if ( node == 3 ) then
    next_boundary_node_t4 = 1
  else
    next_boundary_node_t4 = -1
  end if

  return
end
function next_boundary_node_t6 ( node )

!*****************************************************************************80
!
!! NEXT_BOUNDARY_NODE_T6 returns the next boundary node in a T6 element.
!
!  Reference Element T6:
!
!    |
!    1  3
!    |  |\
!    |  | \
!    S  6  5
!    |  |   \
!    |  |    \
!    0  1--4--2
!    |
!    +--0--R--1-->
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 April 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE, the index of the current node.  An input
!    value of 0 (or any "unusual" value") indicates that the
!    first edge node is desired.
!
!    Output, integer ( kind = 4 ) NEXT_BOUNDARY_NODE_T6, the index of the next 
!    edge node.
!
  implicit none

  integer ( kind = 4 ) next_boundary_node_t6
  integer ( kind = 4 ) node

  if ( node == 1 ) then
    next_boundary_node_t6 = 4
  else if ( node == 4 ) then
    next_boundary_node_t6 = 2
  else if ( node == 2 ) then
    next_boundary_node_t6 = 5
  else if ( node == 5 ) then
    next_boundary_node_t6 = 3
  else if ( node == 3 ) then
    next_boundary_node_t6 = 6
  else if ( node == 6 ) then
    next_boundary_node_t6 = 1
  else
    next_boundary_node_t6 = -1
  end if

  return
end
function next_boundary_node_t10 ( node )

!*****************************************************************************80
!
!! NEXT_BOUNDARY_NODE_T10 returns the next boundary node in a T10 element.
!
!  Reference Element T10:
!
!    |
!    1  10
!    |  |\
!    |  | \
!    |  8  9
!    |  |   \
!    S  |    \
!    |  5  6  7
!    |  |      \
!    |  |       \
!    0  1--2--3--4
!    |
!    +--0----R---1-->
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 April 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE, the index of the current node.  An input
!    value of 0 (or any "unusual" value") indicates that the
!    first edge node is desired.
!
!    Output, integer ( kind = 4 ) NEXT_BOUNDARY_NODE_T10, the index of the next 
!    edge node.
!
  implicit none

  integer ( kind = 4 ) next_boundary_node_t10
  integer ( kind = 4 ) node

  if ( node == 1 ) then
    next_boundary_node_t10 = 2
  else if ( node == 2 ) then
    next_boundary_node_t10 = 3
  else if ( node == 3 ) then
    next_boundary_node_t10 = 4
  else if ( node == 4 ) then
    next_boundary_node_t10 = 7
  else if ( node == 7 ) then
    next_boundary_node_t10 = 9
  else if ( node == 9 ) then
    next_boundary_node_t10 = 10
  else if ( node == 10 ) then
    next_boundary_node_t10 = 8
  else if ( node == 8 ) then
    next_boundary_node_t10 = 5
  else if ( node == 5 ) then
    next_boundary_node_t10 = 1
  else
    next_boundary_node_t10 = -1
  end if

  return
end
subroutine node_reference ( code, r, s, area )

!*****************************************************************************80
!
!! NODE_REFERENCE returns the basis nodes for any available element.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) CODE, identifies the element desired.
!    Legal values include 'Q4', 'Q8', 'Q9', 'Q12', 'Q16', 'QL', 
!    'T3', 'T4', 'T6' and 'T10'.
!
!    Output, real ( kind = 8 ) R(NODE_NUM), S(NODE_NUM), the coordinates 
!    of the basis nodes.
!
!    Output, real ( kind = 8 ) AREA, the area of the element.
!
  implicit none

  real ( kind = 8 ) area
  character ( len = * )  code
  real ( kind = 8 ) r(*)
  real ( kind = 8 ) s(*)
  logical                s_eqi

  if ( s_eqi ( code, 'Q4' ) ) then
    call node_reference_q4 ( r, s, area )
  else if ( s_eqi ( code, 'Q8' ) ) then
    call node_reference_q8 ( r, s, area )
  else if ( s_eqi ( code, 'Q9' ) ) then
    call node_reference_q9 ( r, s, area )
  else if ( s_eqi ( code, 'Q12' ) ) then
    call node_reference_q12 ( r, s, area )
  else if ( s_eqi ( code, 'Q16' ) ) then
    call node_reference_q16 ( r, s, area )
  else if ( s_eqi ( code, 'QL' ) ) then
    call node_reference_ql ( r, s, area )
  else if ( s_eqi ( code, 'T3' ) ) then
    call node_reference_t3 ( r, s, area )
  else if ( s_eqi ( code, 'T4' ) ) then
    call node_reference_t4 ( r, s, area )
  else if ( s_eqi ( code, 'T6' ) ) then
    call node_reference_t6 ( r, s, area )
  else if ( s_eqi ( code, 'T10' ) ) then
    call node_reference_t10 ( r, s, area )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'NODE_REFERENCE - Fatal error!'
    write ( *, '(a)' ) '  Illegal value of CODE = ' // trim ( code )
    stop
  end if

  return
end
subroutine node_reference_q4 ( r, s, area )

!*****************************************************************************80
!
!! NODE_REFERENCE_Q4 returns the basis nodes for a 4 node quadrilateral.
!
!  Reference Element Q4:
!
!    |
!    1  4-------3
!    |  |       |
!    |  |       |
!    S  |       |
!    |  |       |
!    |  |       |
!    0  1-------2
!    |
!    +--0---R---1-->
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) R(4), S(4), the coordinates of the basis nodes.
!
!    Output, real ( kind = 8 ) AREA, the area of the element.
!
  implicit none

  real    ( kind = 8 ) area
  real    ( kind = 8 ) r(4)
  real    ( kind = 8 ) s(4)

  r(1:4) = (/ 0.0D+00, 1.0D+00, 1.0D+00, 0.0D+00 /)
  s(1:4) = (/ 0.0D+00, 0.0D+00, 1.0D+00, 1.0D+00 /)

  area = 1.0D+00

  return
end
subroutine node_reference_q8 ( r, s, area )

!*****************************************************************************80
!
!! NODE_REFERENCE_Q8 returns the basis nodes for an 8 node quadrilateral.
!
!  Discussion:
!
!    This element is known as the quadratic "serendipity" element.
!
!  Reference Element Q8:
!
!    |
!    1  4--7--3
!    |  |     |
!    |  |     |
!    S  8     6
!    |  |     |
!    |  |     |
!    0  1--5--2
!    |
!    +--0--R--1-->
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) R(8), S(8), the coordinates of the basis nodes.
!
!    Output, real ( kind = 8 ) AREA, the area of the element.
!
  implicit none

  real    ( kind = 8 ) area
  real    ( kind = 8 ) r(8)
  real    ( kind = 8 ) s(8)

  r(1:8) = (/ 0.0D+00, 1.0D+00, 1.0D+00, 0.0D+00, &
              0.5D+00, 1.0D+00, 0.5D+00, 0.0D+00 /)
  s(1:8) = (/ 0.0D+00, 0.0D+00, 1.0D+00, 1.0D+00, &
              0.0D+00, 0.5D+00, 1.0D+00, 0.5D+00 /)

  area = 1.0D+00

  return
end
subroutine node_reference_q9 ( r, s, area )

!*****************************************************************************80
!
!! NODE_REFERENCE_Q9 returns the basis nodes for a 9 node quadrilateral.
!
!  Reference Element Q9:
!
!    |
!    1  4--7--3
!    |  |     |
!    |  |     |
!    S  8  9  6
!    |  |     |
!    |  |     |
!    0  1--5--2
!    |
!    +--0--R--1-->
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R(9), S(9), the coordinates of the basis nodes.
!
!    Output, real ( kind = 8 ) AREA, the area of the element.
!
  implicit none

  real    ( kind = 8 ) area
  real    ( kind = 8 ) r(9)
  real    ( kind = 8 ) s(9)

  r(1:9) = (/ 0.0D+00, 1.0D+00, 1.0D+00, 0.0D+00, 0.5D+00, &
              1.0D+00, 0.5D+00, 0.0D+00, 0.5D+00 /)
  s(1:9) = (/ 0.0D+00, 0.0D+00, 1.0D+00, 1.0D+00, 0.0D+00, &
              0.5D+00, 1.0D+00, 0.5D+00, 0.5D+00 /)

  area = 1.0D+00

  return
end
subroutine node_reference_q12 ( r, s, area )

!*****************************************************************************80
!
!! NODE_REFERENCE_Q12 returns the basis nodes for a 12 node quadrilateral.
!
!  Discussion:
!
!    This element is known as the cubic "serendipity" element.
!
!  Reference Element Q12:
!
!    |
!    1  9-10-11-12
!    |  |        |
!    |  7        8
!    S  |        |
!    |  5        6
!    |  |        |
!    0  1--2--3--4
!    |
!    +--0---R---1-->
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) R(12), S(12), the coordinates of the basis nodes.
!
!    Output, real ( kind = 8 ) AREA, the area of the element.
!
  implicit none

  real    ( kind = 8 ), parameter :: a = 0.0D+00
  real    ( kind = 8 ) area
  real    ( kind = 8 ), parameter :: b = 1.0D+00 / 3.0D+00
  real    ( kind = 8 ), parameter :: c = 2.0D+00 / 3.0D+00
  real    ( kind = 8 ), parameter :: d = 1.0D+00
  real    ( kind = 8 ) r(12)
  real    ( kind = 8 ) s(12)

  r(1:12) = (/ a, b, c, d, a, d, a, d, a, b, c, d /)
  s(1:12) = (/ a, a, a, a, b, b, c, c, d, d, d, d /)

  area = 1.0D+00

  return
end
subroutine node_reference_q16 ( r, s, area )

!*****************************************************************************80
!
!! NODE_REFERENCE_Q16 returns the basis nodes for a 16 node quadrilateral.
!
!  Reference Element Q16:
!
!    |
!    1 13--14--15--16
!    |  |   :   :   |
!    |  |   :   :   |
!    |  9..10..11..12
!    S  |   :   :   |
!    |  |   :   :   |
!    |  5...6...7...8
!    |  |   :   :   |
!    |  |   :   :   |  
!    0  1---2---3---4
!    |
!    +--0-----R-----1-->
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) R(16), S(16), the coordinates of the basis nodes.
!
!    Output, real ( kind = 8 ) AREA, the area of the element.
!
  implicit none

  real    ( kind = 8 ) area
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real    ( kind = 8 ) r(16)
  real    ( kind = 8 ) s(16)

  k = 0
  do i = 0, 3
    do j = 0, 3
      k = k + 1
      r(k) = real ( j, kind = 8 ) / 3.0D+00
      s(k) = real ( i, kind = 8 ) / 3.0D+00
    end do
  end do

  area = 1.0D+00

  return
end
subroutine node_reference_ql ( r, s, area )

!*****************************************************************************80
!
!! NODE_REFERENCE_QL returns the basis nodes for a quadratic/linear.
!
!  Reference Element QL:
!
!    |
!    1  4---5---6
!    |  |       |
!    |  |       |
!    S  |       |
!    |  |       |
!    |  |       |
!    0  1---2---3
!    |
!    +--0---R---1-->
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) R(6), S(6), the coordinates of the basis nodes.
!
!    Output, real ( kind = 8 ) AREA, the area of the element.
!
  implicit none

  real    ( kind = 8 ) area
  real    ( kind = 8 ) r(6)
  real    ( kind = 8 ) s(6)

  r(1:6) = (/ 0.0D+00, 0.5D+00, 1.0D+00, 0.0D+00, 0.5D+00, 1.0D+00 /)
  s(1:6) = (/ 0.0D+00, 0.0D+00, 0.0D+00, 1.0D+00, 1.0D+00, 1.0D+00 /)

  area = 1.0D+00

  return
end
subroutine node_reference_t3 ( r, s, area )

!*****************************************************************************80
!
!! NODE_REFERENCE_T3 returns the basis nodes for the 3 node triangle.
!
!  Reference Element T3:
!
!    |
!    1  3
!    |  |\
!    |  | \
!    S  |  \
!    |  |   \
!    |  |    \
!    0  1-----2
!    |
!    +--0--R--1-->
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) R(3), S(3), the coordinates of the basis nodes.
!
!    Output, real ( kind = 8 ) AREA, the area of the element.
!
  implicit none

  real    ( kind = 8 ) area
  real    ( kind = 8 ) r(3)
  real    ( kind = 8 ) s(3)

  r(1:3) = (/ 0.0D+00, 1.0D+00, 0.0D+00 /)
  s(1:3) = (/ 0.0D+00, 0.0D+00, 1.0D+00 /)

  area = 0.5D+00

  return
end
subroutine node_reference_t4 ( r, s, area )

!*****************************************************************************80
!
!! NODE_REFERENCE_T4 returns the basis nodes for the 4 node triangle.
!
!  Reference Element T4:
!
!    |
!    1  3
!    |  |\
!    |  | \
!    S  |  \
!    |  | 4 \
!    |  |    \
!    0  1-----2
!    |
!    +--0--R--1-->
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) R(4), S(4), the coordinates of the basis nodes.
!
!    Output, real ( kind = 8 ) AREA, the area of the element.
!
  implicit none

  real    ( kind = 8 ) area
  real    ( kind = 8 ) r(4)
  real    ( kind = 8 ) s(4)

  r(1:4) = (/ 0.0D+00, 1.0D+00, 0.0D+00, 1.0D+00 / 3.0D+00 /)
  s(1:4) = (/ 0.0D+00, 0.0D+00, 1.0D+00, 1.0D+00 / 3.0D+00 /)

  area = 0.5D+00

  return
end
subroutine node_reference_t6 ( r, s, area )

!*****************************************************************************80
!
!! NODE_REFERENCE_T6 returns the basis nodes for a 6 node triangle.
!
!  Reference Element T6:
!
!    |
!    1  3
!    |  |\
!    |  | \
!    S  6  5
!    |  |   \
!    |  |    \
!    0  1--4--2
!    |
!    +--0--R--1-->
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) R(6), S(6), the coordinates of the basis nodes.
!
!    Output, real ( kind = 8 ) AREA, the area of the element.
!
  implicit none

  real    ( kind = 8 ) area
  real    ( kind = 8 ) r(6)
  real    ( kind = 8 ) s(6)

  r(1:6) = (/ 0.0D+00, 1.0D+00, 0.0D+00, 0.5D+00, 0.5D+00, 0.0D+00 /)
  s(1:6) = (/ 0.0D+00, 0.0D+00, 1.0D+00, 0.0D+00, 0.5D+00, 0.5D+00 /)

  area = 0.5D+00

  return
end
subroutine node_reference_t10 ( r, s, area )

!*****************************************************************************80
!
!! NODE_REFERENCE_T10 returns the basis nodes for a 10 node triangle.
!
!  Reference Element T10:
!
!    |
!    1  10
!    |  |\
!    |  | \
!    |  8  9
!    |  |   \
!    S  |    \
!    |  5  6  7
!    |  |      \
!    |  |       \
!    0  1--2--3--4
!    |
!    +--0----R---1-->
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) R(10), S(10), the coordinates of the basis nodes.
!
!    Output, real ( kind = 8 ) AREA, the area of the element.
!
  implicit none

  real    ( kind = 8 ) area
  real    ( kind = 8 ) r(10)
  real    ( kind = 8 ) s(10)

  r(1) = 0.0D+00
  s(1) = 0.0D+00

  r(2) = 1.0D+00 / 3.0D+00
  s(2) = 0.0D+00

  r(3) = 2.0D+00 / 3.0D+00
  s(3) = 0.0D+00

  r(4) = 1.0D+00
  s(4) = 0.0D+00

  r(5) = 0.0D+00
  s(5) = 1.0D+00 / 3.0D+00

  r(6) = 1.0D+00 / 3.0D+00
  s(6) = 1.0D+00 / 3.0D+00

  r(7) = 2.0D+00 / 3.0D+00
  s(7) = 1.0D+00 / 3.0D+00

  r(8) = 0.0D+00
  s(8) = 2.0D+00 / 3.0D+00

  r(9) = 1.0D+00 / 3.0D+00
  s(9) = 2.0D+00 / 3.0D+00

  r(10) = 0.0D+00
  s(10) = 1.0D+00

  area = 0.5D+00

  return
end
subroutine ns_t6_var_count ( element_num, element_node, node_num, var_node, &
  var_num )

!*****************************************************************************80
!
!! NS_T6_VAR_COUNT counts the Navier Stokes variables on a T6 grid.
!
!  Discussion:
!
!    We are given a mesh of T6 elements, and asked to count, in advance,
!    the number of Navier-Stokes variables associated with the grid.
!    In particular, every node has two velocity variables associated with
!    it, but only a node that is a vertex of the element will also have
!    an associated pressure variable.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 September 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of elements.
!
!    Input, integer ( kind = 4 ) ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM); 
!    ELEMENT_NODE(I,J) is the global index of local node I in element J.
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Output, integer ( kind = 4 ) VAR_NODE(NODE_NUM+1), used to find the variables 
!    associated with a given node, which are in VAR in locations 
!    VAR_NODE(NODE) to VAR_NODE(NODE+1)-1.  Note that the last entry of
!    this array points to the location just after the last location in VAR.
!
!    Output, integer ( kind = 4 ) VAR_NUM, the number of variables.
!
  implicit none

  integer ( kind = 4 ) element_num
  integer ( kind = 4 ), parameter :: element_order = 6
  integer ( kind = 4 ) node_num

  integer ( kind = 4 ) count
  integer ( kind = 4 ) element
  integer ( kind = 4 ) element_node(element_order,element_num)
  integer ( kind = 4 ) node
  integer ( kind = 4 ) num
  integer ( kind = 4 ) order
  integer ( kind = 4 ) var_node(node_num+1)
  integer ( kind = 4 ) var_num
!
!  Our job is easy once we determine which nodes are vertices.
!  So to begin with, let VAR_NODE count the number of variables
!  associated with each node.
!
  var_node(1:node_num) = 2

  do element = 1, element_num
    do order = 1, 3
      node = element_node(order,element)
      var_node(node) = 3
    end do
  end do
!
!  Count them.
!
  var_num = sum ( var_node(1:node_num) )
!
!  Make pointers.
!
  count = 1

  do node = 1, node_num
    num = var_node(node)
    var_node(node) = count
    count = count + num
  end do
  var_node(node_num+1) = count

  return
end
subroutine ns_t6_var_set ( element_num, element_node, node_num, var_node, &
  var_num, var )

!*****************************************************************************80
!
!! NS_T6_VAR_SET sets the Navier Stokes variables on a T6 grid.
!
!  Discussion:
!
!    We are given a mesh of T6 elements, and asked to create the natural
!    list of indices for Navier-Stokes variables associated with each node.
!    In particular, every node has two velocity variables associated with
!    it, but only a node that is a vertex of the element will also have
!    an associated pressure variable.
!
!    The hard work has been done for us alread, because the variables
!    have been counted, and the pointers to the occurrence of the
!    first variable associated with each node have been created.
!
!    The indexing of the nodes can be arbitrary, although a bad
!    indexing will result in a miserably large bandwidth (if band
!    storage is being tried for the stiffness matrix).  Here, we
!    simply try to natural ordering, that is, the variables are
!    numbered in order of the node with which they are associated.
!
!    For the Navier Stokes problem on a T6 grid, we take it as
!    understood that each node has either 2 or 3 variables associated
!    with it, that the first two are always the horizontal and
!    then vertical velocity coefficients, and that the third, if
!    present, is a pressure coefficient.
!
!    In other settings, it might be necessary not merely to assign
!    the variables an index, but also to identify them as to type.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 September 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of elements.
!
!    Input, integer ( kind = 4 ) ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM); 
!    ELEMENT_NODE(I,J) is the global index of local node I in element J.
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, integer ( kind = 4 ) VAR_NODE(NODE_NUM+1), used to find the variables 
!    associated with a given node, which are in VAR in locations 
!    VAR_NODE(NODE) to VAR_NODE(NODE+1)-1.  Note that the last entry of
!    this array points to the location just after the last location in VAR.
!
!    Input, integer ( kind = 4 ) VAR_NUM, the number of variables.
!
!    Output, integer ( kind = 4 ) VAR(VAR_NUM), the indexes of the variables, which
!    are simply 1, 2, 3, ..., VAR_NUM.
!
  implicit none

  integer ( kind = 4 ) element_num
  integer ( kind = 4 ), parameter :: element_order = 6
  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) var_num

  integer ( kind = 4 ) element_node(element_order,element_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) var(var_num)
  integer ( kind = 4 ) var_node(node_num+1)

  do i = 1, var_num
    var(i) = i
  end do 
  
  return
end
function order_code ( code )

!*****************************************************************************80
!
!! ORDER_CODE returns the order for each element.
!
!  Discussion:
!
!    CODE  Order  Definition
!    ----  -----  ----------
!    Q4     4     4 node linear Lagrange/serendipity quadrilateral;
!    Q8     8     8 node quadratic serendipity quadrilateral;
!    Q9     9     9 node quadratic Lagrange quadrilateral;
!    Q12   12     12 node cubic serendipity quadrilateral;
!    Q16   16     16 node cubic Lagrange quadrilateral;
!    QL     6     6 node linear/quadratic quadrilateral;
!    T3     3     3 node linear triangle;
!    T4     4     4 node cubic bubble triangle
!    T6     6     6 node quadratic triangle;
!    T10   10     10 node cubic triangle.
! 
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) CODE, the code for the element.
!
!    Output, integer ( kind = 4 ) ORDER_CODE, the order of the element.
!
  implicit none

  character ( len = * ) code
  integer ( kind = 4 ) order_code
  logical s_eqi

  if ( s_eqi ( code, 'Q4' ) ) then
    order_code = 4
  else if ( s_eqi ( code, 'Q8' ) ) then
    order_code = 8
  else if ( s_eqi ( code, 'Q9' ) ) then
    order_code = 9
  else if ( s_eqi ( code, 'Q12' ) ) then
    order_code = 12
  else if ( s_eqi ( code, 'Q16' ) ) then
    order_code = 16
  else if ( s_eqi ( code, 'QL' ) ) then
    order_code = 6
  else if ( s_eqi ( code, 'T3' ) ) then
    order_code = 3
  else if ( s_eqi ( code, 'T4' ) ) then
    order_code = 4
  else if ( s_eqi ( code, 'T6' ) ) then
    order_code = 6
  else if ( s_eqi ( code, 'T10' ) ) then
    order_code = 10
  else
    order_code = -1
  end if

  return
end
subroutine physical_to_reference_t3 ( t, n, phy, ref )

!*****************************************************************************80
!
!! PHYSICAL_TO_REFERENCE_T3 maps physical points to reference points.
!
!  Discussion:
!
!    Given the vertices of an order 3 physical triangle and a point 
!    (X,Y) in the physical triangle, the routine computes the value 
!    of the corresponding image point (XSI,ETA) in reference space.
!
!    This routine is also appropriate for an order 4 triangle, assuming
!    that the fourth node is always the centroid of the triangle.
!
!    This routine may be appropriate for an order 6
!    triangle, if the mapping between reference and physical space
!    is linear.  This implies, in particular, that the sides of the
!    image triangle are straight and that the "midside" nodes in the
!    physical triangle are halfway along the sides of
!    the physical triangle.
!
!  Reference Element T3:
!
!    |
!    1  3
!    |  |\
!    |  | \
!    S  |  \
!    |  |   \
!    |  |    \
!    0  1-----2
!    |
!    +--0--R--1-->
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 June 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T(2,3), the X and Y coordinates
!    of the vertices.  The vertices are assumed to be the images of
!    (0,0), (1,0) and (0,1) respectively.
!
!    Input, integer ( kind = 4 ) N, the number of points to transform.
!
!    Input, real ( kind = 8 ) PHY(2,N), the coordinates of physical points
!    to be transformed.
!
!    Output, real ( kind = 8 ) REF(2,N), the coordinates of the corresponding
!    points in the reference space.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) phy(2,n)
  real    ( kind = 8 ) ref(2,n)
  real    ( kind = 8 ) t(2,3)

  ref(1,1:n) = ( ( t(2,3) - t(2,1) ) * ( phy(1,1:n) - t(1,1) )   &
               - ( t(1,3) - t(1,1) ) * ( phy(2,1:n) - t(2,1) ) ) &
             / ( ( t(2,3) - t(2,1) ) * ( t(1,2)     - t(1,1) )   &
               - ( t(1,3) - t(1,1) ) * ( t(2,2)     - t(2,1) ) )

  ref(2,1:n) = ( ( t(1,2) - t(1,1) ) * ( phy(2,1:n) - t(2,1) )   &
               - ( t(2,2) - t(2,1) ) * ( phy(1,1:n) - t(1,1) ) ) &
             / ( ( t(2,3) - t(2,1) ) * ( t(1,2)     - t(1,1) )   &
               - ( t(1,3) - t(1,1) ) * ( t(2,2)     - t(2,1) ) )

  return
end
subroutine points_plot ( file_name, node_num, node_xy, node_label )

!*****************************************************************************80
!
!! POINTS_PLOT plots a pointset.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the output file.
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of points.
!
!    Input, real ( kind = 8 ) NODE_XY(2,NODE_NUM), the nodes.
!
!    Input, logical NODE_LABEL, is TRUE if the nodes should be labeled.
!
!  Local parameters:
!
!    Local, integer CIRCLE_SIZE, controls the size of the circles depicting
!    the nodes.  Currently set to 5.  3 is pretty small, and 1 is
!    barely visible.
!
  implicit none

  integer   ( kind = 4 ) node_num

  integer   ( kind = 4 ), parameter :: circle_size = 5
  integer   ( kind = 4 ) delta
  character ( len = * )  file_name
  integer   ( kind = 4 ) file_unit
  integer   ( kind = 4 ) ios
  integer   ( kind = 4 ) node
  logical                node_label
  real ( kind = 8 ) node_xy(2,node_num)
  character ( len = 40 ) string
  real ( kind = 8 ) x_max
  real ( kind = 8 ) x_min
  integer   ( kind = 4 ) x_ps
  integer   ( kind = 4 ) :: x_ps_max = 576
  integer   ( kind = 4 ) :: x_ps_max_clip = 594
  integer   ( kind = 4 ) :: x_ps_min = 36
  integer   ( kind = 4 ) :: x_ps_min_clip = 18
  real ( kind = 8 ) x_scale
  real ( kind = 8 ) y_max
  real ( kind = 8 ) y_min
  integer   ( kind = 4 ) y_ps
  integer   ( kind = 4 ) :: y_ps_max = 666
  integer   ( kind = 4 ) :: y_ps_max_clip = 684
  integer   ( kind = 4 ) :: y_ps_min = 126
  integer   ( kind = 4 ) :: y_ps_min_clip = 108
  real ( kind = 8 ) y_scale
!
!  We need to do some figuring here, so that we can determine
!  the range of the data, and hence the height and width
!  of the piece of paper.
!
  x_max = maxval ( node_xy(1,1:node_num) )
  x_min = minval ( node_xy(1,1:node_num) )
  x_scale = x_max - x_min

  x_max = x_max + 0.05D+00 * x_scale
  x_min = x_min - 0.05D+00 * x_scale
  x_scale = x_max - x_min

  y_max = maxval ( node_xy(2,1:node_num) )
  y_min = minval ( node_xy(2,1:node_num) )
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

  open ( unit = file_unit, file = file_name, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'POINTS_PLOT - Fatal error!'
    write ( *, '(a)' ) '  Can not open output file.'
    return
  end if

  write ( file_unit, '(a)' ) '%!PS-Adobe-3.0 EPSF-3.0'
  write ( file_unit, '(a)' ) '%%Creator: points_plot.f90'
  write ( file_unit, '(a)' ) '%%Title: ' // trim ( file_name )
  write ( file_unit, '(a)' ) '%%Pages: 1'
  write ( file_unit, '(a,i3,2x,i3,2x,i3,2x,i3)' ) '%%BoundingBox: ', &
    x_ps_min, y_ps_min, x_ps_max, y_ps_max
  write ( file_unit, '(a)' ) '%%Document-Fonts: Times-Roman'
  write ( file_unit, '(a)' ) '%%LanguageLevel: 1'
  write ( file_unit, '(a)' ) '%%EndComments'
  write ( file_unit, '(a)' ) '%%BeginProlog'
  write ( file_unit, '(a)' ) '/inch {72 mul} def'
  write ( file_unit, '(a)' ) '%%EndProlog'
  write ( file_unit, '(a)' ) '%%Page: 1 1'
  write ( file_unit, '(a)' ) 'save'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%  Set the RGB line color to very light gray.'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '0.900  0.900  0.900 setrgbcolor'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%  Draw a gray border around the page.'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) 'newpath'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', x_ps_min, y_ps_min, ' moveto'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', x_ps_max, y_ps_min, ' lineto'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', x_ps_max, y_ps_max, ' lineto'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', x_ps_min, y_ps_max, ' lineto'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', x_ps_min, y_ps_min, ' lineto'
  write ( file_unit, '(a)' ) 'stroke'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%  Set the RGB line color to black.'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '0.000  0.000  0.000 setrgbcolor'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%  Set the font and its size.'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '/Times-Roman findfont'
  write ( file_unit, '(a)' ) '0.50 inch scalefont'
  write ( file_unit, '(a)' ) 'setfont'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%  Print a title.'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%  210  702  moveto'
  write ( file_unit, '(a)' ) '%  (Pointset)  show'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%  Define a clipping polygon.'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) 'newpath'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', &
    x_ps_min_clip, y_ps_min_clip, ' moveto'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', &
    x_ps_max_clip, y_ps_min_clip, ' lineto'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', &
    x_ps_max_clip, y_ps_max_clip, ' lineto'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', &
    x_ps_min_clip, y_ps_max_clip, ' lineto'
  write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', &
    x_ps_min_clip, y_ps_min_clip, ' lineto'
  write ( file_unit, '(a)' ) 'clip newpath'
!
!  Draw the nodes.
!
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%  Draw filled dots at each node.'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%  Set the RGB color to blue.'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '0.000  0.150  0.750 setrgbcolor'
  write ( file_unit, '(a)' ) '%'

  do node = 1, node_num

    x_ps = int ( &
      ( ( x_max - node_xy(1,node)         ) * real ( x_ps_min, kind = 8 )   &
      + (         node_xy(1,node) - x_min ) * real ( x_ps_max, kind = 8 ) ) &
      / ( x_max                   - x_min ) )

    y_ps = int ( &
      ( ( y_max - node_xy(2,node)         ) * real ( y_ps_min, kind = 8 )   &
      + (         node_xy(2,node) - y_min ) * real ( y_ps_max, kind = 8 ) ) &
      / ( y_max                   - y_min ) )

    write ( file_unit, '(a,i4,2x,i4,2x,i4,2x,a)' ) 'newpath ', x_ps, y_ps, &
      circle_size, '0 360 arc closepath fill'

  end do
!
!  Label the nodes.
!
  if ( node_label ) then

    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '%  Label the nodes:'
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '%  Set the RGB color to darker blue.'
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '0.000  0.250  0.850 setrgbcolor'
    write ( file_unit, '(a)' ) '/Times-Roman findfont'
    write ( file_unit, '(a)' ) '0.20 inch scalefont'
    write ( file_unit, '(a)' ) 'setfont'

    do node = 1, node_num

      x_ps = int ( &
        ( ( x_max - node_xy(1,node)         ) * real ( x_ps_min, kind = 8 )   &
        + (       + node_xy(1,node) - x_min ) * real ( x_ps_max, kind = 8 ) ) &
        / ( x_max                   - x_min ) )

      y_ps = int ( &
        ( ( y_max - node_xy(2,node)         ) * real ( y_ps_min, kind = 8 )   &
        + (         node_xy(2,node) - y_min ) * real ( y_ps_max, kind = 8 ) ) &
        / ( y_max                   - y_min ) )

      write ( string, '(i4)' ) node
      string = adjustl ( string )

      write ( file_unit, '(i4,2x,i4,a)' ) x_ps, y_ps+5, &
        ' moveto (' // trim ( string ) // ') show'

    end do

  end if

  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) 'restore  showpage'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%  End of page.'
  write ( file_unit, '(a)' ) '%'
  write ( file_unit, '(a)' ) '%%Trailer'
  write ( file_unit, '(a)' ) '%%EOF'
  close ( unit = file_unit )

  return
end
subroutine poly ( code, rexp, sexp )

!*****************************************************************************80
!
!! POLY returns the polynomial terms associated with any available element.
!
!  Discussion:
!
!    Given coefficients A(I), the polynomial interpolant at (R,S) is
!
!      P(R,S) = sum ( 1 <= I <= N ) A(I) * R**REXP(I) * S**SEXP(I) 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) CODE, identifies the element desired.
!    Legal values include 'Q4', 'Q8', 'Q9', 'Q12', 'Q16', 'QL', 'T3', 
!    'T4', 'T6' and 'T10'.
!
!    Output, integer ( kind = 4 ) REXP(N), SEXP(N), the powers of R and S 
!    associated with each polynomial.  The value of N, the dimension of these
!    arrays, can be determined by a call to ORDER_CODE.
!
  implicit none

  character ( len = * )  code
  integer   ( kind = 4 ) rexp(*)
  integer   ( kind = 4 ) sexp(*)
  logical                s_eqi

  if ( s_eqi ( code, 'Q4' ) ) then
    call poly_q4 ( rexp, sexp )
  else if ( s_eqi ( code, 'Q8' ) ) then
    call poly_q8 ( rexp, sexp )
  else if ( s_eqi ( code, 'Q9' ) ) then
    call poly_q9 ( rexp, sexp )
  else if ( s_eqi ( code, 'Q12' ) ) then
    call poly_q12 ( rexp, sexp )
  else if ( s_eqi ( code, 'Q16' ) ) then
    call poly_q16 ( rexp, sexp )
  else if ( s_eqi ( code, 'QL' ) ) then
    call poly_ql ( rexp, sexp )
  else if ( s_eqi ( code, 'T3' ) ) then
    call poly_t3 ( rexp, sexp )
  else if ( s_eqi ( code, 'T4' ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'POLY - Fatal error!'
    write ( *, '(a)' ) '  The T4 element does not follow the pattern.' 
    stop
  else if ( s_eqi ( code, 'T6' ) ) then
    call poly_t6 ( rexp, sexp )
  else if ( s_eqi ( code, 'T10' ) ) then
    call poly_t10 ( rexp, sexp )
  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'POLY - Fatal error!'
    write ( *, '(a)' ) '  Illegal value of CODE = "' // trim ( code ) // '".'
    stop

  end if

  return
end
subroutine poly_q4 ( rexp, sexp )

!*****************************************************************************80
!
!! POLY_Q4 returns the monomials associated with a 4 node quadrilateral.
!
!  Reference Element Q4:
!
!    |
!    1  4-----3
!    |  |     |
!    |  |     |
!    S  |     |
!    |  |     |
!    |  |     |
!    0  1-----2
!    |
!    +--0--R--1-->
!
!  Formula:
!
!    Given coefficients A(I), the polynomial interpolant at (R,S) is
!
!      P(R,S) = sum ( 1 <= I <= N ) A(I) * R**REXP(I) * S**SEXP(I) 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) REXP(4), SEXP(4), the powers of R and S 
!    associated with each polynomial.
!
  implicit none

  integer ( kind = 4 ) rexp(4)
  integer ( kind = 4 ) sexp(4)

  rexp(1:4) = (/ 0, 0, 1, 1 /)
  sexp(1:4) = (/ 0, 1, 0, 1 /)

  return
end
subroutine poly_q8 ( rexp, sexp )

!*****************************************************************************80
!
!! POLY_Q8 returns the monomials associated with an 8 node quadrilateral.
!
!  Reference Element Q8:
!
!    |
!    1  4--7--3
!    |  |     |
!    |  |     |
!    S  8     6
!    |  |     |
!    |  |     |
!    0  1--5--2
!    |
!    +--0--R--1-->
!
!  Formula:
!
!    Given coefficients A(I), the polynomial interpolant at (R,S) is
!
!      P(R,S) = sum ( 1 <= I <= N ) A(I) * R**REXP(I) * S**SEXP(I) 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) REXP(8), SEXP(8), the powers of R and S 
!    associated with each monomial.
!
  implicit none

  integer ( kind = 4 ) rexp(8)
  integer ( kind = 4 ) sexp(8)

  rexp(1:8) = (/ 0, 0, 1, 0, 1, 2, 1, 2 /)
  sexp(1:8) = (/ 0, 1, 0, 2, 1, 0, 2, 1 /)

  return
end
subroutine poly_q9 ( rexp, sexp )

!*****************************************************************************80
!
!! POLY_Q9 returns the monomials associated with a 9 node quadrilateral.
!
!  Reference Element Q9:
!
!    |
!    1  4--7--3
!    |  |     |
!    |  |     |
!    S  8  9  6
!    |  |     |
!    |  |     |
!    0  1--5--2
!    |
!    +--0--R--1-->
!
!  Formula:
!
!    Given coefficients A(I), the polynomial interpolant at (R,S) is
!
!      P(R,S) = sum ( 1 <= I <= N ) A(I) * R**REXP(I) * S**SEXP(I) 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) REXP(9), SEXP(9), the powers of R and S 
!    associated with each monomial.
!
  implicit none

  integer ( kind = 4 ) rexp(9)
  integer ( kind = 4 ) sexp(9)

  rexp(1:9) = (/ 0, 0, 1, 0, 1, 2, 1, 2, 2 /)
  sexp(1:9) = (/ 0, 1, 0, 2, 1, 0, 2, 1, 2 /)

  return
end
subroutine poly_q12 ( rexp, sexp )

!*****************************************************************************80
!
!! POLY_Q12 returns the monomials associated with a 12 node quadrilateral.
!
!  Reference Element Q12:
!
!    |
!    1  9-10-11-12
!    |  |        |
!    |  7        8
!    S  |        |
!    |  5        6
!    |  |        |
!    0  1--2--3--4
!    |
!    +--0---R---1-->
!
!  Formula:
!
!    Given coefficients A(I), the polynomial interpolant at (R,S) is
!
!      P(R,S) = sum ( 1 <= I <= N ) A(I) * R**REXP(I) * S**SEXP(I) 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) REXP(12), SEXP(12), the powers of R and S 
!    associated with each monomial.
!
  implicit none

  integer ( kind = 4 ) rexp(12)
  integer ( kind = 4 ) sexp(12)

  rexp(1:12) = (/ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 1, 3 /)
  sexp(1:12) = (/ 0, 1, 0, 2, 1, 0, 3, 2, 1, 0, 3, 1 /)

  return
end
subroutine poly_q16 ( rexp, sexp )

!*****************************************************************************80
!
!! POLY_Q16 returns the monomials associated with a 16 node quadrilateral.
!
!  Reference Element Q16:
!
!    |
!    1 13--14--15--16
!    |  |   :   :   |
!    |  |   :   :   |
!    |  9..10..11..12
!    S  |   :   :   |
!    |  |   :   :   |
!    |  5...6...7...8
!    |  |   :   :   |
!    |  |   :   :   |  
!    0  1---2---3---4
!    |
!    +--0-----R-----1-->
!
!  Formula:
!
!    Given coefficients A(I), the polynomial interpolant at (R,S) is
!
!      P(R,S) = sum ( 1 <= I <= N ) A(I) * R**REXP(I) * S**SEXP(I) 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) REXP(16), SEXP(16), the powers of R and S 
!    associated with each monomial.
!
  implicit none

  integer ( kind = 4 ) rexp(16)
  integer ( kind = 4 ) sexp(16)

  rexp(1:16) = (/ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 1, 2, 3, 2, 3, 3 /)
  sexp(1:16) = (/ 0, 1, 0, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 3, 2, 3 /)

  return
end
subroutine poly_ql ( rexp, sexp )

!*****************************************************************************80
!
!! POLY_QL returns the monomials for a quadratic/linear quadrilateral.
!
!  Reference Element QL:
!
!    |
!    1  4---5---6
!    |  |       |
!    |  |       |
!    S  |       |
!    |  |       |
!    |  |       |
!    0  1---2---3
!    |
!    +--0---R---1-->
!
!  Formula:
!
!    Given coefficients A(I), the polynomial interpolant at (R,S) is
!
!      P(R,S) = sum ( 1 <= I <= N ) A(I) * R**REXP(I) * S**SEXP(I) 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) REXP(N), SEXP(N), the powers of R and S
!    associated with each monomial.
!
  implicit none

  integer ( kind = 4 ) rexp(6)
  integer ( kind = 4 ) sexp(6)

  rexp(1:6) = (/ 0, 0, 1, 1, 2, 2 /)
  sexp(1:6) = (/ 0, 1, 0, 1, 0, 1 /)

  return
end
subroutine poly_t3 ( rexp, sexp )

!*****************************************************************************80
!
!! POLY_T3 returns the monomials associated with a 3 node triangle.
!
!  Reference Element T3:
!
!    |
!    1  3
!    |  |\
!    |  | \
!    S  |  \
!    |  |   \
!    |  |    \
!    0  1-----2
!    |
!    +--0--R--1-->
!
!  Formula:
!
!    Given coefficients A(I), the polynomial interpolant at (R,S) is
!
!      P(R,S) = sum ( 1 <= I <= N ) A(I) * R**REXP(I) * S**SEXP(I) 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) REXP(N), SEXP(N), the powers of R and S 
!    associated with each monomial.
!
  implicit none

  integer ( kind = 4 ) rexp(3)
  integer ( kind = 4 ) sexp(3)

  rexp(1:3) = (/ 0, 0, 1 /)
  sexp(1:3) = (/ 0, 1, 0 /)

  return
end
subroutine poly_t6 ( rexp, sexp )

!*****************************************************************************80
!
!! POLY_T6 returns the monomials associated with a 6 node triangle.
!
!  Reference Element T6:
!
!    |
!    1  3
!    |  |\
!    |  | \
!    S  6  5
!    |  |   \
!    |  |    \
!    0  1--4--2
!    |
!    +--0--R--1-->
!
!  Formula:
!
!    Given coefficients A(I), the polynomial interpolant at (R,S) is
!
!      P(R,S) = sum ( 1 <= I <= N ) A(I) * R**REXP(I) * S**SEXP(I) 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) REXP(6), SEXP(6), the powers of R and S 
!    associated with each monomial.
!
  implicit none

  integer ( kind = 4 ) rexp(6)
  integer ( kind = 4 ) sexp(6)

  rexp(1:6) = (/ 0, 0, 1, 0, 1, 2 /)
  sexp(1:6) = (/ 0, 1, 0, 2, 1, 0 /)

  return
end
subroutine poly_t10 ( rexp, sexp )

!*****************************************************************************80
!
!! POLY_T10 returns the monomials associated with a 10 node triangle.
!
!  Reference Element T10:
!
!    |
!    1  10
!    |  |\
!    |  | \
!    |  8  9
!    |  |   \
!    S  |    \
!    |  5  6  7
!    |  |      \
!    |  |       \
!    0  1--2--3--4
!    |
!    +--0----R---1-->
!
!  Formula:
!
!    Given coefficients A(I), the polynomial interpolant at (R,S) is
!
!      P(R,S) = sum ( 1 <= I <= N ) A(I) * R**REXP(I) * S**SEXP(I) 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) REXP(10), SEXP(10), the powers of R and S 
!    associated with each monomial.
!
  implicit none

  integer ( kind = 4 ) rexp(10)
  integer ( kind = 4 ) sexp(10)

  rexp(1:10) = (/ 0, 0, 1, 0, 1, 2, 0, 1, 2, 3 /)
  sexp(1:10) = (/ 0, 1, 0, 2, 1, 0, 3, 2, 1, 0 /)

  return
end
function r8_huge ( )

!*****************************************************************************80
!
!! R8_HUGE returns the largest legal R8.
!
!  Discussion:
!
!    FORTRAN90 provides a built-in routine HUGE ( X ) that
!    is more suitable for this purpose.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 March 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) R8_HUGE, a "huge" value.
!
  implicit none

  real    ( kind = 8 ) r8_huge

  r8_huge = 1.0D+30

  return
end
function r8_power ( r, p )

!*****************************************************************************80
!
!! R8_POWER computes the P-th power of R.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, the base.
!
!    Input, integer ( kind = 4 ) P, the power, which may be negative.
!
!    Output, real ( kind = 8 ) R8_POWER, the value of the P-th power of R.
!
  implicit none

  integer ( kind = 4 ) p
  real    ( kind = 8 ) r
  real    ( kind = 8 ) r8_power
  real    ( kind = 8 ) value
!
!  Special case.  R^0 = 1.
!
  if ( p == 0 ) then

    value = 1.0D+00
!
!  Special case.  All positive powers of 0 are 0.
!
  else if ( r == 0.0D+00 ) then

    value = 0.0D+00

  else if ( 1 <= p ) then
    value = r**p
  else
    value = 1.0D+00 / r**(-p)
  end if

  r8_power = value

  return
end
subroutine r8_swap ( x, y )

!*****************************************************************************80
!
!! R8_SWAP switches two R8's.
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
!    Input/output, real ( kind = 8 ) X, Y.  On output, the values of X and
!    Y have been interchanged.
!
  implicit none

  real    ( kind = 8 ) x
  real    ( kind = 8 ) y
  real    ( kind = 8 ) z

  z = x
  x = y
  y = z

  return
end
function r8_uniform_01 ( seed )

!*****************************************************************************80
!
!! R8_UNIFORM_01 returns a unit pseudorandom R8.
!
!  Discussion:
!
!    An R8 is a real ( kind = 8 ) value.
!
!    For now, the input quantity SEED is an integer variable.
!
!    This routine implements the recursion
!
!      seed = 16807 * seed mod ( 2**31 - 1 )
!      r8_uniform_01 = seed / ( 2**31 - 1 )
!
!    The integer arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!    If the initial seed is 12345, then the first three computations are
!
!      Input     Output      R8_UNIFORM_01
!      SEED      SEED
!
!         12345   207482415  0.096616
!     207482415  1790989824  0.833995
!    1790989824  2035175616  0.947702
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley Interscience, page 95, 1998.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which should
!    NOT be 0. On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R8_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer ( kind = 4 ) k
  real    ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + 2147483647
  end if
!
!  Although SEED can be represented exactly as a 32 bit integer,
!  it generally cannot be represented exactly as a 32 bit real number!
!
  r8_uniform_01 = real ( seed, kind = 8 ) * 4.656612875D-10

  return
end
subroutine r8ge_fa ( n, a, pivot, info )

!*****************************************************************************80
!
!! R8GE_FA performs a LINPACK style PLU factorization of an R8GE matrix.
!
!  Discussion:
!
!    The R8GE storage format is used for a general M by N matrix.  A storage
!    space is made for each logical entry.  The two dimensional logical
!    array is mapped to a vector, in which storage is by columns.
!
!    R8GE_FA is a simplified version of the LINPACK routine DGEFA.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Input/output, real ( kind = 8 ) A(N,N), the matrix to be factored.
!    On output, A contains an upper triangular matrix and the multipliers
!    which were used to obtain it.  The factorization can be written
!    A = L * U, where L is a product of permutation and unit lower
!    triangular matrices and U is upper triangular.
!
!    Output, integer ( kind = 4 ) PIVOT(N), a vector of pivot indices.
!
!    Output, integer ( kind = 4 ) INFO, singularity flag.
!    0, no singularity detected.
!    nonzero, the factorization failed on the INFO-th step.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) pivot(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l

  info = 0

  do k = 1, n-1
!
!  Find L, the index of the pivot row.
!
    l = k
    do i = k+1, n
      if ( abs ( a(l,k) ) < abs ( a(i,k) ) ) then
        l = i
      end if
    end do

    pivot(k) = l
!
!  If the pivot index is zero, the algorithm has failed.
!
    if ( a(l,k) == 0.0D+00 ) then
      info = k
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8GE_FA - Fatal error!'
      write ( *, '(a,i8)' ) '  Zero pivot on step ', info
      return
    end if
!
!  Interchange rows L and K if necessary.
!
    if ( l /= k ) then
      call r8_swap ( a(l,k), a(k,k) )
    end if
!
!  Normalize the values that lie below the pivot entry A(K,K).
!
    a(k+1:n,k) = -a(k+1:n,k) / a(k,k)
!
!  Row elimination with column indexing.
!
    do j = k+1, n

      if ( l /= k ) then
        call r8_swap ( a(l,j), a(k,j) )
      end if

      a(k+1:n,j) = a(k+1:n,j) + a(k+1:n,k) * a(k,j)

    end do

  end do

  pivot(n) = n

  if ( a(n,n) == 0.0D+00 ) then
    info = n
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8GE_FA - Fatal error!'
    write ( *, '(a,i8)' ) '  Zero pivot on step ', info
  end if

  return
end
subroutine r8ge_inverse ( n, a, pivot )

!*****************************************************************************80
!
!! R8GE_INVERSE computes the inverse of a matrix factored by R8GE_FA.
!
!  Discussion:
!
!    The R8GE storage format is used for a general M by N matrix.  A storage
!    space is made for each logical entry.  The two dimensional logical
!    array is mapped to a vector, in which storage is by columns.
!
!    R8GE_INVERSE is a simplified standalone version of the LINPACK routine
!    DGEDI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 September 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix A.
!
!    Input/output, real ( kind = 8 ) A(N,N).
!    On input, the factor information computed by R8GE_FA.
!    On output, the inverse matrix.
!
!    Input, integer ( kind = 4 ) PIVOT(N), the pivot vector from R8GE_FA.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) pivot(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real    ( kind = 8 ) temp
  real    ( kind = 8 ) work(n)
!
!  Compute Inverse(U).
!
  do k = 1, n

    a(k,k) = 1.0D+00 / a(k,k)
    a(1:k-1,k) = -a(1:k-1,k) * a(k,k)

    do j = k + 1, n

      temp = a(k,j)
      a(k,j) = 0.0D+00
      a(1:k,j) = a(1:k,j) + a(1:k,k) * temp

    end do

  end do
!
!  Form Inverse(U) * Inverse(L).
!
  do k = n - 1, 1, -1

    work(k+1:n) = a(k+1:n,k)
    a(k+1:n,k) = 0.0D+00

    do j = k + 1, n
      a(1:n,k) = a(1:n,k) + a(1:n,j) * work(j)
    end do

    if ( pivot(k) /= k ) then

      do i = 1, n
        call r8_swap ( a(i,k), a(i,pivot(k)) )
      end do

    end if

  end do

  return
end
subroutine r8mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8MAT_PRINT prints an R8MAT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix.
!
!    Input, character ( len = * ) TITLE, a title to be printed.
!
  implicit none

  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = * )  title

  call r8mat_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine r8mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8MAT_PRINT_SOME prints some of an R8MAT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 March 2005
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

  real    ( kind = 8 ) a(m,n)
  character ( len = 14 ) ctemp(incx)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i8,6x)') j
    end do

    write ( *, '(''       Col'',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '       Row'
    write ( *, '(a)' ) ' '

    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi

      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( a(i,j) == real ( int ( a(i,j) ), kind = 8 ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) a(i,j)
        else
          write ( ctemp(j2), '(g14.6)' ) a(i,j)
        end if

      end do

      write ( *, '(2x,i8,5a14)' ) i, ( ctemp(j), j = 1, inc )

    end do

  end do

  return
end
subroutine r8mat_write ( output_file_name, m, n, table )

!*****************************************************************************80
!
!! R8MAT_WRITE writes an R8MAT file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) OUTPUT_FILE_NAME, the output file name.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) TABLE(M,N), the table data.
!
  implicit none

  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) n

  integer   ( kind = 4 ) j
  character ( len = * )  output_file_name
  integer   ( kind = 4 ) output_status
  integer   ( kind = 4 ) output_unit
  character ( len = 30 ) string
  real ( kind = 8 ) table(m,n)
!
!  Open the file.
!
  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_file_name, &
    status = 'replace', iostat = output_status )

  if ( output_status /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_WRITE - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the output file "' // &
      trim ( output_file_name ) // '" on unit ', output_unit
    output_unit = -1
    stop
  end if
!
!  Create a format string.
!
  write ( string, '(a1,i8,a1,i8,a1,i8,a1)' ) '(', m, 'g', 14, '.', 6, ')'
!
!  Write the data.
!
  do j = 1, n
    write ( output_unit, string ) table(1:m,j)
  end do
!
!  Close the file.
!
  close ( unit = output_unit )

  return
end
subroutine reference_sample ( code, seed, r, s )

!*****************************************************************************80
!
!! REFERENCE_SAMPLE samples a reference element.
!
!  Discussion:
!
!    The routine either samples the unit triangle or the unit square.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) CODE, identifies the element desired.
!    Legal values include 'Q4', 'Q8', 'Q9', 'Q12', 'Q16', 'QL', 'T3', 
!    'T4', 'T6' and 'T10'.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) R, S, a random point in the reference element.
!
  implicit none

  character ( len = * )  code
  real ( kind = 8 ) r
  real ( kind = 8 ) r8_uniform_01 
  real ( kind = 8 ) s
  integer   ( kind = 4 ) seed

  if ( code(1:1) == 'Q' .or. code(1:1) == 'q' ) then

    r = r8_uniform_01 ( seed )
    s = r8_uniform_01 ( seed )

  else if ( code(1:1) == 'T' .or. code(1:1) == 't' ) then

    r = r8_uniform_01 ( seed )
    s = r8_uniform_01 ( seed )

    if ( 1.0D+00 < r + s ) then
      r = 1.0D+00 - r
      s = 1.0D+00 - s
    end if

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'REFERENCE_SAMPLE - Fatal error!'
    write ( *, '(a)' ) '  Illegal code = "' // trim ( code ) // '".'
    stop

  end if

  return
end
subroutine reference_to_physical_q4 ( q4, n, rs, xy )

!*****************************************************************************80
!
!! REFERENCE_TO_PHYSICAL_Q4 maps Q4 reference points to physical points.
!
!  Discussion:
!
!    XY(R,S) = XY(0,0) * (1-R) * (1-S)
!            + XY(1,0) *    R  * (1-S)
!            + XY(1,1) *    R  *    S
!            + XY(0,1) * (1-R) *    S
!
!  Reference Element Q4:
!
!    |
!    1  4-----3
!    |  |     |
!    |  |     |
!    S  |     |
!    |  |     |
!    |  |     |
!    0  1-----2
!    |
!    +--0--R--1-->
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 February 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) Q4(2,4), the coordinates of the vertices.
!    The vertices are assumed to be the images of the reference vertices
!    (0,0), (1,0), (1,1) and (0,1) respectively.
!
!    Input, integer ( kind = 4 ) N, the number of points to transform.
!
!    Input, real ( kind = 8 ) RS(2,N), (R,S) points in the reference element.
!
!    Output, real ( kind = 8 ) XY(2,N), (X,Y) points in the physical element.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) psi(4,n)
  real    ( kind = 8 ) q4(2,4)
  real    ( kind = 8 ) rs(2,n)
  real    ( kind = 8 ) xy(2,n)

  psi(1,1:n) = ( 1.0D+00 - rs(1,1:n) ) * ( 1.0D+00 - rs(2,1:n) )
  psi(2,1:n) =             rs(1,1:n)   * ( 1.0D+00 - rs(2,1:n) )
  psi(3,1:n) =             rs(1,1:n)   *             rs(2,1:n)
  psi(4,1:n) = ( 1.0D+00 - rs(1,1:n) ) *             rs(2,1:n)

  xy(1:2,1:n) = matmul ( q4(1:2,1:4), psi(1:4,1:n) )

  return
end
subroutine reference_to_physical_t3 ( t, n, ref, phy )

!*****************************************************************************80
!
!! REFERENCE_TO_PHYSICAL_T3 maps T3 reference points to physical points.
!
!  Discussion:
!
!    Given the vertices of an order 3 physical triangle and a point 
!    (XSI,ETA) in the reference triangle, the routine computes the value 
!    of the corresponding image point (X,Y) in physical space.
!
!    This routine is also appropriate for an order 4 triangle,
!    as long as the fourth node is the centroid of the triangle.
!
!    This routine may also be appropriate for an order 6
!    triangle, if the mapping between reference and physical space
!    is linear.  This implies, in particular, that the sides of the
!    image triangle are straight and that the "midside" nodes in the
!    physical triangle are halfway along the sides of
!    the physical triangle.
!
!  Reference Element T3:
!
!    |
!    1  3
!    |  |\
!    |  | \
!    S  |  \
!    |  |   \
!    |  |    \
!    0  1-----2
!    |
!    +--0--R--1-->
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 June 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T(2,3), the coordinates of the vertices.  
!    The vertices are assumed to be the images of (0,0), (1,0) and 
!    (0,1) respectively.
!
!    Input, integer ( kind = 4 ) N, the number of points to transform.
!
!    Input, real ( kind = 8 ) REF(2,N), points in the reference element.
!
!    Output, real ( kind = 8 ) PHY(2,N), corresponding points in the
!    physical element.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real    ( kind = 8 ) phy(2,n)
  real    ( kind = 8 ) ref(2,n)
  real    ( kind = 8 ) t(2,3)

  do i = 1, 2
    phy(i,1:n) = t(i,1) * ( 1.0D+00 - ref(1,1:n) - ref(2,1:n) ) &
               + t(i,2) *             ref(1,1:n)                &
               + t(i,3) *                          ref(2,1:n)
  end do

  return
end
subroutine reference_to_physical_t6 ( t, n, ref, phy )

!*****************************************************************************80
!
!! REFERENCE_TO_PHYSICAL_T6 maps T6 reference points to physical points.
!
!  Discussion:
!
!    Given the vertices of an order 6 physical triangle and a point 
!    (XSI,ETA) in the reference triangle, the routine computes the value 
!    of the corresponding image point (X,Y) in physical space.
!
!    The mapping from (XSI,ETA) to (X,Y) has the form:
!
!      X(ETA,XSI) = A1 * XSI**2 + B1 * XSI*ETA + C1 * ETA**2
!                 + D1 * XSI    + E1 * ETA     + F1
!
!      Y(ETA,XSI) = A2 * XSI**2 + B2 * XSI*ETA + C2 * ETA**2
!                 + D2 * XSI    + E2 * ETA     + F2
!
!  Reference Element T6:
!
!    |
!    1  3
!    |  |\
!    |  | \
!    S  6  5
!    |  |   \
!    |  |    \
!    0  1--4--2
!    |
!    +--0--R--1-->
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 June 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T(2,6), the coordinates of the vertices.  
!    The vertices are assumed to be the images of (0,0), (1,0),
!    (0,1),(1/2,0), (1/2,1/2) and (0,1/2) respectively.
!
!    Input, integer ( kind = 4 ) N, the number of points to transform.
!
!    Input, real ( kind = 8 ) REF(2,N), points in the reference element.
!
!    Output, real ( kind = 8 ) PHY(2,N), corresponding points in the
!    physical element.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(2)
  real    ( kind = 8 ) b(2)
  real    ( kind = 8 ) c(2)
  real    ( kind = 8 ) d(2)
  real    ( kind = 8 ) e(2)
  real    ( kind = 8 ) f(2)
  integer ( kind = 4 ) i
  real    ( kind = 8 ) phy(2,n)
  real    ( kind = 8 ) ref(2,n)
  real    ( kind = 8 ) t(2,6)

  do i = 1, 2

    a(i) =   2.0D+00 * t(i,1) + 2.0D+00 * t(i,2)                    &
           - 4.0D+00 * t(i,4) 

    b(i) =   4.0D+00 * t(i,1)                                       & 
           - 4.0D+00 * t(i,4) + 4.0D+00 * t(i,5) - 4.0D+00 * t(i,6)

    c(i) =   2.0D+00 * t(i,1)                    + 2.0D+00 * t(i,3) &
                                                 - 4.0D+00 * t(i,6)

    d(i) = - 3.0D+00 * t(i,1) -           t(i,2)                    &
           + 4.0D+00 * t(i,4)

    e(i) = - 3.0D+00 * t(i,1)                    -           t(i,3) &
                                                 + 4.0D+00 * t(i,6)
    f(i) =             t(i,1)

  end do

  do i = 1, 2
    phy(i,1:n) = a(i) * ref(1,1:n) * ref(1,1:n) &
               + b(i) * ref(1,1:n) * ref(2,1:n) &
               + c(i) * ref(2,1:n) * ref(2,1:n) &
               + d(i) * ref(1,1:n) &
               + e(i) * ref(2,1:n) &
               + f(i)
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
function s_eqi ( s1, s2 )

!*****************************************************************************80
!
!! S_EQI is a case insensitive comparison of two strings for equality.
!
!  Example:
!
!    S_EQI ( 'Anjana', 'ANJANA' ) is TRUE.
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
subroutine s_l2norm ( psi_num, element_num, quad_num, element_area, &
  quad_weight, psi_quad, s_coef, l2norm )

!*****************************************************************************80
!
!! S_L2NORM computes the "big" L2 norm of a scalar function over a region.
!
!  Discussion:
!
!    It is assumed that a set of finite element basis functions PSI
!    have been defined over a collection of elements that compose
!    the region.  Moreover, integrals over the region are assumed to
!    be approximated by applying a fixed quadrature rule to all the
!    elements.
!
!    Finally, we assume that we have a scalar function S(X), which
!    is represented as a linear combination of basis functions, and
!    it is desired to determine the L2 norm of S.
!
!    This routine estimates the integral
!
!      Sqrt ( Integral ( X in Omega ) S(X) * S(X) dOmega )
!
!    using the finite element representation of S, and applying the
!    given quadrature rule.
!
!    It assumes that a (probably very large) array of data is available,
!    recording the value of each basis function PSI in every element
!    at every quadrature point.  If this is true, then the computation
!    becomes very simple.
!
!    If your problem is small or sufficient memory is available, this
!    may be an efficient computation.  It requires that the value of
!    all the basis functions be stored for all the elements and all
!    the quadrature points.  That particular information need only
!    be computed once.
!
!    Actually, no assumptions are made here about the dimension of the
!    space, so this same code can handle problems in 1D, 2D, 3D and
!    so on.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 March 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) PSI_NUM, the number of global element functions.
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of elements.
!
!    Input, integer ( kind = 4 ) QUAD_NUM, the number of points in the 
!    quadrature rule.
!
!    Input, real ( kind = 8 ) ELEMENT_AREA(ELEMENT_NUM), the area of 
!    each element.
!
!    Input, real ( kind = 8 ) QUAD_WEIGHT(QUAD_NUM), the quadrature
!    weights associated with the quadrature points.
!
!    Input, real ( kind = 8 ) PSI_QUAD(PSI_NUM,ELEMENT_NUM,QUAD_NUM), the 
!    value of the I-th PSI function in the J-th element at the K-th 
!    quadrature point.
!
!    Input, real ( kind = 8 ) S_COEF(PSI_NUM), the coefficients of the 
!    PSI functions associated with the scalar function S.
!
!    Output, L2NORM, the L2 norm of the scalar function S over the region.
!
  implicit none

  integer ( kind = 4 ) element_num
  integer ( kind = 4 ) psi_num
  integer ( kind = 4 ) quad_num

  real    ( kind = 8 ) element_area(element_num)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real    ( kind = 8 ) l2norm
  real    ( kind = 8 ) psi_quad(psi_num,element_num,quad_num)
  real    ( kind = 8 ) quad_weight(quad_num)
  real    ( kind = 8 ) s(element_num,quad_num)
  real    ( kind = 8 ) s_coef(psi_num)
  real    ( kind = 8 ) t(quad_num)
  real    ( kind = 8 ) u
!
!  #1: Sum over all basis functions to get the value of S in each element
!  at each quadrature point.
!
!  The MATMUL function requires that one of its arguments be shaped
!  like a vector, and one like a 2 dimensional matrix, so we have
!  to insert a loop on the quadrature points.
!
  do j = 1, quad_num

    s(1:element_num,j) = matmul ( &
      s_coef(1:psi_num), psi_quad(1:psi_num,1:element_num,j) )

  end do
!
!  #2: Sum over all elements to get the value of S * S weighted by its element
!  area.  SUM expects to see vector quantities, so we have a loop on
!  quadrature points.
!
  do k = 1, quad_num
    t(k) = sum ( s(1:element_num,k)**2 * element_area(1:element_num) )
  end do
!
!  #3: Sum over all quadrature points weighted by the quadrature weight.
!
  u = dot_product ( t(1:quad_num), quad_weight(1:quad_num) )

  l2norm = sqrt ( u )

  return
end
subroutine serene ( type, ve, vn, vne, vnw, vs, vse, vsw, vw, vterp )

!*****************************************************************************80
!
!! SERENE interpolates data using a Q8 element.
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
!    Input, character ( len = 2 ) TYPE, tells SERENE the geometry of the
!    finite element that surrounds the point of interest.  The options
!    are displayed in the following table, which suggests the meaning
!    of each option by its position:
!
!        |   |
!     NW * N * NE
!        |   |
!     -*-*-*-*-*-
!        |   |
!      W * C * E
!        |   |
!     -*-*-*-*-*-
!        |   |
!     SW * S * SE
!        |   |
!
!    Input, real ( kind = 8 ) VE, VN, VNE, VNW, VS, VSE, VSW, VW,
!    are the values of the function at the nodes to the east,
!    north, northeast, northwest, south, southeast, southwest and
!    west of the point of interest.  If the finite element is of
!    type 'C', then all 8 values are needed.  However, if the
!    finite element is of type 'SE', for instance, then only three
!    values are needed, namely VE, VN, and VNW, since these are
!    the only node positions defined in such a finite element.
!
!    Output, real ( kind = 8 ) VTERP, the interpolated value of the
!    function at the point of interest.
!
  implicit none

  real    ( kind = 8 ) eta
  real    ( kind = 8 ) pe
  real    ( kind = 8 ) pn
  real    ( kind = 8 ) pne
  real    ( kind = 8 ) pnw
  real    ( kind = 8 ) ps
  real    ( kind = 8 ) pse
  real    ( kind = 8 ) psw
  real    ( kind = 8 ) pw
  real    ( kind = 8 ) r8_huge
  character ( len = 2 ) type
  real    ( kind = 8 ) ve
  real    ( kind = 8 ) vn
  real    ( kind = 8 ) vne
  real    ( kind = 8 ) vnw
  real    ( kind = 8 ) vs
  real    ( kind = 8 ) vse
  real    ( kind = 8 ) vsw
  real    ( kind = 8 ) vw
  real    ( kind = 8 ) vterp
  real    ( kind = 8 ) xsi
!
!  To make this routine more general, simply pass in the values of XSI
!  and ETA at which the interpolated value is desired.
!
!  By setting XSI = ETA = 0, we are asking for the interpolated value
!  at the center of the finite element.
!
  xsi = 0.0D+00
  eta = 0.0D+00
!
!  8 node center
!
!  Polynomial space is spanned by:
!
!         1
!       x    y
!    x^2  xy  y^2
!      x^2y xy^2
!
!
!    ^   1    4--7--3
!    |        !     !
!    E        !     !
!    T   0    8  X  6
!    A        !     !
!    |        !     !
!    V  -1    1--5--2
!
!            -1  0  1
!
!           <---XSI--->
!
  if ( type == 'C' ) then

    psw = - 0.25D+00 * ( 1.0D+00 - xsi ) * ( 1.0D+00 - eta ) &
      * ( 1.0D+00 + xsi + eta )
    pse = - 0.25D+00 * ( 1.0D+00 + xsi ) * ( 1.0D+00 - eta ) &
      * ( 1.0D+00 - xsi + eta )
    pne = - 0.25D+00 * ( 1.0D+00 + xsi ) * ( 1.0D+00 + eta ) &
      * ( 1.0D+00 - xsi - eta )
    pnw = - 0.25D+00 * ( 1.0D+00 - xsi ) * ( 1.0D+00 + eta ) &
      * ( 1.0D+00 + xsi - eta )
    ps =    0.50D+00 * ( 1.0D+00 - xsi ) * ( 1.0D+00 + xsi ) &
      * ( 1.0D+00 - eta )
    pe =    0.50D+00 * ( 1.0D+00 + xsi ) * ( 1.0D+00 + eta ) & 
      * ( 1.0D+00 - eta )
    pn =    0.50D+00 * ( 1.0D+00 - xsi ) * ( 1.0D+00 + xsi ) & 
      * ( 1.0D+00 + eta )
    pw =    0.50D+00 * ( 1.0D+00 - xsi ) * ( 1.0D+00 + eta ) & 
      * ( 1.0D+00 - eta )

    vterp = vsw * psw + vse * pse + vne * pne + vnw * pnw &
      + vs * ps + ve * pe + vn * pn + vw * pw
!
!  5 node side
!
!    ^   1
!    |
!    E
!    T   0    8  X  6
!    A        !     !
!    |        !     !
!    V  -1    1--5--2
!
!            -1  0  1
!
!           <---XSI--->
!
  else if ( type == 'N' ) then

    psw =  0.5D+00 * ( xsi - 1.0D+00 ) * ( 1.0D+00 + xsi + eta )
    pse = -0.5D+00 * ( xsi + 1.0D+00 ) * ( 1.0D+00 - xsi + eta )
    ps =  -          ( xsi + 1.0D+00 ) * ( xsi - 1.0D+00 )
    pe =   0.5D+00 * ( xsi + 1.0D+00 ) * ( eta + 1.0D+00 )
    pw =  -0.5D+00 * ( xsi - 1.0D+00 ) * ( eta + 1.0D+00 )

    vterp = vsw * psw + vse * pse + vs * ps + ve * pe + vw * pw
!
!    ^   1    4--7
!    |        !
!    E        !
!    T   0    8  X
!    A        !
!    |        !
!    V  -1    1--5
!
!            -1  0  1
!
!           <---XSI--->
!
  else if ( type == 'E' ) then

    pse =  0.5D+00 * ( eta - 1.0D+00 ) * ( 1.0D+00 + xsi + eta )
    pne = -0.5D+00 * ( eta + 1.0D+00 ) * ( 1.0D+00 + xsi - eta )
    ps =  -0.5D+00 * ( xsi + 1.0D+00 ) * ( eta - 1.0D+00 )
    pn =   0.5D+00 * ( xsi + 1.0D+00 ) * ( eta + 1.0D+00 )
    pw =  -          ( eta + 1.0D+00 ) * ( eta - 1.0D+00 )

    vterp = vse * pse + vne * pne + vs * ps + vn * pn + vw * pw
!
!  5 node side
!
!    ^   1       7--3
!    |              !
!    E              !
!    T   0       X  6
!    A              !
!    |              !
!    V  -1       5--2
!
!            -1  0  1
!
!           <---XSI--->
!
  else if ( type == 'W' ) then

    pse =   0.5D+00 * ( eta - 1.0D+00 ) * ( 1.0D+00 - xsi + eta )
    pne = - 0.5D+00 * ( eta + 1.0D+00 ) * ( 1.0D+00 - xsi - eta )
    ps =    0.5D+00 * ( xsi - 1.0D+00 ) * ( eta - 1.0D+00 )
    pe =  -           ( eta - 1.0D+00 ) * ( eta + 1.0D+00 )
    pn =  - 0.5D+00 * ( xsi - 1.0D+00 ) * ( eta + 1.0D+00 )

    vterp = vse * pse + vne * pne + vs * ps + ve * pe + vn * pn
!
!  5 node side
!
!    ^   1    4--7--3
!    |        !     !
!    E        !     !
!    T   0    8  X  6
!    A
!    |
!    V  -1
!
!            -1  0  1
!
!           <---XSI--->
!
  else if ( type == 'S' ) then

    pne = - 0.5D+00 * ( xsi + 1.0D+00 ) * ( 1.0D+00 - xsi - eta )
    pnw =   0.5D+00 * ( xsi - 1.0D+00 ) * ( 1.0D+00 + xsi - eta )
    pe =  - 0.5D+00 * ( eta - 1.0D+00 ) * ( xsi + 1.0D+00 )
    pn =  -           ( xsi + 1.0D+00 ) * ( xsi - 1.0D+00 )
    pw =    0.5D+00 * ( eta - 1.0D+00 ) * ( xsi - 1.0D+00 )

    vterp = vne * pne + vnw * pnw + ve * pe + vn * pn + vw * pw
!
!  3 node corner
!
!  Polynomial space is spanned by:
!
!         1
!       x    y
!
!
!    ^   1
!    |
!    E
!    T   0    8  X
!    A        !
!    |        !
!    V  -1    1--5
!
!            -1  0  1
!
!           <---XSI--->
!
  else if ( type == 'NE' ) then

    psw = - 1.0D+00 - xsi - eta
    ps =    1.0D+00 + xsi
    pw =    1.0D+00       + eta

    vterp = vsw * psw + vs * ps + vw * pw
!
!  3 node corner
!
!  Polynomial space is spanned by:
!
!         1
!       x    y
!
!    ^   1
!    |
!    E
!    T   0       X  6
!    A              !
!    |              !
!    V  -1       5--2
!
!            -1  0  1
!
!           <---XSI--->
!
  else if ( type == 'NW' ) then

    pse = 1.0D+00 + xsi - eta
    ps =  1.0D+00 - xsi
    pe =  1.0D+00       + eta

    vterp = + vse * pse + vs * ps + ve * pe
!
!  3 node corner
!
!  Polynomial space is spanned by:
!         1
!       x    y
!
!
!    ^   1    4--7
!    |        !
!    E        !
!    T   0    8  X
!    A
!    |
!    V  -1
!
!            -1  0  1
!
!           <---XSI--->
!
  else if ( type == 'SE' ) then

    pnw = - 1.0D+00 - xsi + eta
    pn =    1.0D+00 + xsi
    pw =    1.0D+00       - eta

    vterp = vnw * pnw + vn * pn + vw * pw
!
!  3 node corner
!
!  Polynomial space is spanned by:
!
!         1
!       x    y
!
!    ^   1       7--3
!    |              !
!    E              !
!    T   0       X  6
!    A
!    |
!    V  -1
!
!            -1  0  1
!
!           <---XSI--->
!
  else if ( type == 'SW' ) then

    pne = - 1.0D+00 + xsi + eta
    pe =    1.0D+00       - eta
    pn =    1.0D+00 - xsi

    vterp = vne * pne + ve * pe + vn * pn

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SERENE - Fatal error!'
    write ( *, '(a)' ) '  Illegal value of TYPE = "' // trim ( type ) // '".'
    vterp = - r8_huge ( vterp )
    stop

  end if

  return
end
subroutine shape ( code, r, s, t, dtdr, dtds )

!*****************************************************************************80
!
!! SHAPE evaluates shape functions for any available reference element.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 April 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) CODE, identifies the element.
!    Legal values include 'Q4', 'Q8', 'Q9', 'Q12', 'Q16', 'QL', 
!    'T3', 'T4', 'T6' and 'T10'.
!
!    Input, real ( kind = 8 ) R, S, the reference coordinates of a point.
!
!    Output, real ( kind = 8 ) T(N), the basis functions at the point.
!
!    Output, real ( kind = 8 ) DTDR(N), the R basis derivatives at the point.
!
!    Output, real ( kind = 8 ) DTDS(N), the S basis derivatives at the point.
!
  implicit none

  character ( len = * )  code
  real ( kind = 8 ) dtdr(*)
  real ( kind = 8 ) dtds(*)
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  logical                s_eqi
  real ( kind = 8 ) t(*)

  if ( s_eqi ( code, 'Q4' ) ) then
    call shape_q4 ( r, s, t, dtdr, dtds )
  else if( s_eqi ( code, 'Q8' ) ) then
    call shape_q8 ( r, s, t, dtdr, dtds )
  else if ( s_eqi ( code, 'Q9' ) ) then
    call shape_q9 ( r, s, t, dtdr, dtds )
  else if ( s_eqi ( code, 'Q12' ) ) then
    call shape_q12 ( r, s, t, dtdr, dtds )
  else if ( s_eqi ( code, 'Q16' ) ) then
    call shape_q16 ( r, s, t, dtdr, dtds )
  else if ( s_eqi ( code, 'QL' ) ) then
    call shape_ql ( r, s, t, dtdr, dtds )
  else if ( s_eqi ( code, 'T3' ) ) then
    call shape_t3 ( r, s, t, dtdr, dtds )
  else if ( s_eqi ( code, 'T4' ) ) then
    call shape_t4 ( r, s, t, dtdr, dtds )
  else if ( s_eqi ( code, 'T6' ) ) then
    call shape_t6 ( r, s, t, dtdr, dtds )
  else if ( s_eqi ( code, 'T10' ) ) then
    call shape_t10 ( r, s, t, dtdr, dtds )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SHAPE - Fatal error!'
    write ( *, '(a)' ) '  Unrecognized code = "' // trim ( code ) // '".'
    stop
  end if

  return
end
subroutine shape_q4 ( r, s, t, dtdr, dtds )

!*****************************************************************************80
!
!! SHAPE_Q4 evaluates shape functions for a 4 node reference quadrilateral.
!
!  Reference Element Q4:
!
!    |
!    1  4-----3
!    |  |     |
!    |  |     |
!    S  |     |
!    |  |     |
!    |  |     |
!    0  1-----2
!    |
!    +--0--R--1-->
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 April 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, S, the reference coordinates of a point.
!
!    Output, real ( kind = 8 ) T(4), the basis functions at the point.
!
!    Output, real ( kind = 8 ) DTDR(4), the R basis derivatives at the point.
!
!    Output, real ( kind = 8 ) DTDS(4), the S basis derivatives at the point.
!
  implicit none

  real    ( kind = 8 ) dtdr(4)
  real    ( kind = 8 ) dtds(4)
  real    ( kind = 8 ) r
  real    ( kind = 8 ) s
  real    ( kind = 8 ) t(4)

  t(1) = ( 1.0D+00 - r ) * ( 1.0D+00 - s )
  t(2) =             r   * ( 1.0D+00 - s )
  t(3) =             r   *             s
  t(4) = ( 1.0D+00 - r ) *             s

  dtdr(1) = - 1.0D+00 + s
  dtdr(2) =   1.0D+00 - s     
  dtdr(3) =             s
  dtdr(4) =           - s

  dtds(1) = - 1.0D+00 + r
  dtds(2) =           - r
  dtds(3) =             r
  dtds(4) =   1.0D+00 - r

  return
end
subroutine shape_q8 ( r, s, t, dtdr, dtds )

!*****************************************************************************80
!
!! SHAPE_Q8 evaluates shape functions for an 8 node reference quadrilateral.
!
!  Discussion:
!
!    This element is known as the "serendipity" element.
!
!  Reference Element Q8:
!
!    |
!    1  4--7--3
!    |  |     |
!    |  |     |
!    S  8     6
!    |  |     |
!    |  |     |
!    0  1--5--2
!    |
!    +--0--R--1-->
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 April 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, S, the reference coordinates of a point.
!
!    Output, real ( kind = 8 ) T(8), the basis functions at the point.
!
!    Output, real ( kind = 8 ) DTDR(8), the R basis derivatives at the point.
!
!    Output, real ( kind = 8 ) DTDS(8), the S basis derivatives at the point.
!
  implicit none

  real    ( kind = 8 ) dtdr(8)
  real    ( kind = 8 ) dtds(8)
  real    ( kind = 8 ) r
  real    ( kind = 8 ) s
  real    ( kind = 8 ) t(8)

  t(1) =                 ( r - 1.0D+00 )     * ( s - 1.0D+00 ) &
    * ( 1.0D+00 - 2.0D+00 * r - 2.0D+00 * s )
  t(2) =             r                       * ( s - 1.0D+00 ) &
    * ( 1.0D+00 - 2.0D+00 * r + 2.0D+00 * s )
  t(3) =             r                   * s                   &
    * ( 2.0D+00 * r + 2.0D+00 * s - 3.0D+00 )
  t(4) =                 ( r - 1.0D+00 ) * s                   &
    * ( 2.0D+00 * r - 2.0D+00 * s + 1.0D+00 )
  t(5) =   4.0D+00 * r * ( r - 1.0D+00 )     * ( s - 1.0D+00 )
  t(6) = - 4.0D+00 * r                   * s * ( s - 1.0D+00 )
  t(7) = - 4.0D+00 * r * ( r - 1.0D+00 ) * s    
  t(8) =   4.0D+00 *     ( r - 1.0D+00 ) * s * ( s - 1.0D+00 )

  dtdr(1) = ( s - 1.0D+00 ) * ( - 4.0D+00 * r - 2.0D+00 * s + 3.0D+00 )
  dtdr(2) = ( s - 1.0D+00 ) * ( - 4.0D+00 * r + 2.0D+00 * s + 1.0D+00 )
  dtdr(3) =   s         * (   4.0D+00 * r + 2.0D+00 * s - 3.0D+00 )
  dtdr(4) =   s         * (   4.0D+00 * r - 2.0D+00 * s - 1.0D+00 )
  dtdr(5) =   4.0D+00 * ( 2.0D+00 * r - 1.0D+00 )     * ( s - 1.0D+00 )
  dtdr(6) = - 4.0D+00 *                     s * ( s - 1.0D+00 )
  dtdr(7) = - 4.0D+00 * ( 2.0D+00 * r - 1.0D+00 ) * s
  dtdr(8) =   4.0D+00 *                     s * ( s - 1.0D+00 )

  dtds(1) = ( r - 1.0D+00 ) * ( - 4.0D+00 * s - 2.0D+00 * r + 3.0D+00 )
  dtds(2) =   r *       (   4.0D+00 * s - 2.0D+00 * r - 1.0D+00 )
  dtds(3) =   r *       (   4.0D+00 * s + 2.0D+00 * r - 3.0D+00 )
  dtds(4) = ( r - 1.0D+00 ) * ( - 4.0D+00 * s + 2.0D+00 * r + 1.0D+00 )
  dtds(5) =   4.0D+00 * r * ( r - 1.0D+00 )
  dtds(6) = - 4.0D+00 * r               * ( 2.0D+00 * s - 1.0D+00 )
  dtds(7) = - 4.0D+00 * r * ( r - 1.0D+00 )
  dtds(8) =   4.0D+00 *     ( r - 1.0D+00 ) * ( 2.0D+00 * s - 1.0D+00 )

  return
end
subroutine shape_q9 ( r, s, t, dtdr, dtds )

!*****************************************************************************80
!
!! SHAPE_Q9 evaluates shape functions for a 9 node reference quadrilateral.
!
!  Reference Element Q9:
!
!    |
!    1  4--7--3
!    |  |     |
!    |  |     |
!    S  8  9  6
!    |  |     |
!    |  |     |
!    0  1--5--2
!    |
!    +--0--R--1-->
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 April 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, S, the reference coordinates of a point.
!
!    Output, real ( kind = 8 ) T(9), the basis functions at the point.
!
!    Output, real ( kind = 8 ) DTDR(9), the R basis derivatives at the point.
!
!    Output, real ( kind = 8 ) DTDS(9), the S basis derivatives at the point.
!
  implicit none

  real    ( kind = 8 ) dtdr(9)
  real    ( kind = 8 ) dtds(9)
  real    ( kind = 8 ) r
  real    ( kind = 8 ) s
  real    ( kind = 8 ) t(9)

  t(1) =    4.0D+00 * ( r - 1.0D+00 ) * ( r - 0.5D+00 ) * ( s - 1.0D+00 ) &
            * ( s - 0.5D+00 )
  t(2) =    4.0D+00 * r * ( r - 0.5D+00 ) * ( s - 1.0D+00 ) * ( s - 0.5D+00 )
  t(3) =    4.0D+00 * r * ( r - 0.5D+00 ) * s * ( s - 0.5D+00 )
  t(4) =    4.0D+00 * ( r - 1.0D+00 ) * ( r - 0.5D+00 ) * s * ( s - 0.5D+00 )
  t(5) = -  8.0D+00 * r * ( r - 1.0D+00 ) * ( s - 1.0D+00 ) * ( s - 0.5D+00 )
  t(6) = -  8.0D+00 * r * ( r - 0.5D+00 ) * s * ( s - 1.0D+00 )
  t(7) = -  8.0D+00 * r * ( r - 1.0D+00 ) * s * ( s - 0.5D+00 )
  t(8) = -  8.0D+00 * ( r - 1.0D+00 ) * ( r - 0.5D+00 ) * s * ( s - 1.0D+00 )
  t(9) =   16.0D+00 * r * ( r - 1.0D+00 ) * s * ( s - 1.0D+00 )

  dtdr(1) =   4.0D+00 * ( 2.0D+00 * r - 1.5D+00 ) * ( s - 1.0D+00 ) &
              * ( s - 0.5D+00 )
  dtdr(2) =   4.0D+00 * ( 2.0D+00 * r - 0.5D+00 ) * ( s - 1.0D+00 ) &
              * ( s - 0.5D+00 )
  dtdr(3) =   4.0D+00 * ( 2.0D+00 * r - 0.5D+00 ) * s * ( s - 0.5D+00 )
  dtdr(4) =   4.0D+00 * ( 2.0D+00 * r - 1.5D+00 ) * s * ( s - 0.5D+00 )

  dtdr(5) = - 8.0D+00 * ( 2.0D+00 * r - 1.0D+00 ) * ( s - 1.0D+00 ) &
              * ( s - 0.5D+00 )
  dtdr(6) = - 8.0D+00 * ( 2.0D+00 * r - 0.5D+00 ) * s * ( s - 1.0D+00 )
  dtdr(7) = - 8.0D+00 * ( 2.0D+00 * r - 1.0D+00 ) * s * ( s - 0.5D+00 )
  dtdr(8) = - 8.0D+00 * ( 2.0D+00 * r - 1.5D+00 ) * s * ( s - 1.0D+00 )
  dtdr(9) =  16.0D+00 * ( 2.0D+00 * r - 1.0D+00 ) * s * ( s - 1.0D+00 )

  dtds(1) =   4.0D+00 * ( r - 1.0D+00 ) * ( r - 0.5D+00 ) &
              * ( 2.0D+00 * s - 1.5D+00 )
  dtds(2) =   4.0D+00 * r * ( r - 0.5D+00 ) * ( 2.0D+00 * s - 1.5D+00 )
  dtds(3) =   4.0D+00 * r * ( r - 0.5D+00 ) * ( 2.0D+00 * s - 0.5D+00 )
  dtds(4) =   4.0D+00 * ( r - 1.0D+00 ) * ( r - 0.5D+00 ) &
            * ( 2.0D+00 * s - 0.5D+00 )
  dtds(5) = - 8.0D+00 * r * ( r - 1.0D+00 ) * ( 2.0D+00 * s - 1.5D+00 )
  dtds(6) = - 8.0D+00 * r * ( r - 0.5D+00 ) * ( 2.0D+00 * s - 1.0D+00 )
  dtds(7) = - 8.0D+00 * r * ( r - 1.0D+00 ) * ( 2.0D+00 * s - 0.5D+00 ) 
  dtds(8) = - 8.0D+00 * ( r - 1.0D+00 ) * ( r - 0.5D+00 ) &
            * ( 2.0D+00 * s - 1.0D+00 )
  dtds(9) =  16.0D+00 * r * ( r - 1.0D+00 ) * ( 2.0D+00 * s - 1.0D+00 )

  return
end
subroutine shape_q12 ( r, s, t, dtdr, dtds )

!*****************************************************************************80
!
!! SHAPE_Q12 evaluates shape functions for a 12 node reference quadrilateral.
!
!  Reference Element Q12:
!
!    |
!    1  9-10-11-12
!    |  |        |
!    |  7        8
!    S  |        |
!    |  5        6
!    |  |        |
!    0  1--2--3--4
!    |
!    +--0---R---1-->
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, S, the reference coordinates of a point.
!
!    Output, real ( kind = 8 ) T(12), the basis functions at the point.
!
!    Output, real ( kind = 8 ) DTDR(12), the R basis derivatives at the point.
!
!    Output, real ( kind = 8 ) DTDS(12), the S basis derivatives at the point.
!
  implicit none

  real    ( kind = 8 ) a
  real    ( kind = 8 ) b
  real    ( kind = 8 ) c
  real    ( kind = 8 ) corner
  real    ( kind = 8 ) d
  real    ( kind = 8 ) dcdr
  real    ( kind = 8 ) dcds
  real    ( kind = 8 ) dtdr(12)
  real    ( kind = 8 ) dtds(12)
  real    ( kind = 8 ) r
  real    ( kind = 8 ) s
  real    ( kind = 8 ) t(12)

  a = 0.0D+00
  b = 1.0D+00 / 3.0D+00
  c = 2.0D+00 / 3.0D+00
  d = 1.0D+00

  corner = 9.0D+00 * ( &
      ( 2.0D+00 * r - 1.0D+00 ) * ( 2.0D+00 * r - 1.0D+00 ) &
    + ( 2.0D+00 * s - 1.0D+00 ) * ( 2.0D+00 * s - 1.0D+00 ) ) - 10.0D+00

  t(1) =     0.125D+00  * ( r - d ) * ( s - d ) * corner
  t(2) =  - 13.5D+00    * ( r - a ) * ( r - c ) * ( r - d ) * ( s - d )
  t(3) =    13.5D+00    * ( r - a ) * ( r - b ) * ( r - d ) * ( s - d )
  t(4) =   - 0.125D+00  * ( r - a ) * ( s - d ) * corner
  t(5) =  - 13.5D+00    * ( r - d ) * ( s - a ) * ( s - c ) * ( s - d ) 
  t(6) =    13.5D+00    * ( r - a ) * ( s - a ) * ( s - c ) * ( s - d )
  t(7) =    13.5D+00    * ( r - d ) * ( s - a ) * ( s - b ) * ( s - d )
  t(8) =  - 13.5D+00    * ( r - a ) * ( s - a ) * ( s - b ) * ( s - d )
  t(9) =   - 0.125D+00  * ( r - d ) * ( s - a ) * corner
  t(10) =   13.5D+00    * ( r - a ) * ( r - c ) * ( r - d ) * ( s - a )
  t(11) = - 13.5D+00    * ( r - a ) * ( r - b ) * ( r - d ) * ( s - a )
  t(12) =    0.125D+00  * ( r - a ) * ( s - a ) * corner
 
  dcdr = 36.0D+00 * ( 2.0D+00 * r - 1.0D+00 )

  dtdr(1) =  0.125 * ( s - d ) * ( ( r - d ) * dcdr + corner )
  dtdr(2) =  - 13.5D+00 * ( s - d ) * ( 3.0D+00 * r * r &
    - 2.0D+00 * ( a + c + d ) * r + a * c + c * d + d * a ) 
  dtdr(3) =    13.5D+00 * ( s - d ) * ( 3.0D+00 * r * r &
    - 2.0D+00 * ( a + b + d ) * r + a * b + b * d + d * a )
  dtdr(4) = - 0.125D+00 * ( s - d ) * ( ( r - a ) * dcdr + corner )
  dtdr(5) = - 13.5D+00 * ( s - a ) * ( s - c ) * ( s - d ) 
  dtdr(6) =   13.5D+00 * ( s - a ) * ( s - c ) * ( s - d )
  dtdr(7) =   13.5D+00 * ( s - a ) * ( s - b ) * ( s - d )
  dtdr(8) = - 13.5D+00 * ( s - a ) * ( s - b ) * ( s - d )
  dtdr(9) = - 0.125D+00 * ( s - a ) * ( ( r - d ) * dcdr + corner )
  dtdr(10) =   13.5D+00 * ( s - a ) * ( 3.0D+00 * r * r &
    - 2.0D+00 * ( a + c + d ) * r + a * c + c * d + d * a ) 
  dtdr(11) = - 13.5D+00 * ( s - a ) * ( 3.0D+00 * r * r &
    - 2.0D+00 * ( a + b + d ) * r + a * b + b * d + d * a )
  dtdr(12) = 0.125D+00 * ( s - a ) * ( ( r - a ) * dcdr + corner )

  dcds = 36.0D+00 * ( 2.0D+00 * s - 1.0D+00 )

  dtds(1) =  0.125D+00 * ( r - d ) * ( corner + ( s - d ) * dcds )
  dtds(2) =  - 13.5D+00 * ( r - a ) * ( r - c ) * ( r - d ) 
  dtds(3) =  13.5D+00 * ( r - a ) * ( r - b ) * ( r - d )
  dtds(4) = - 0.125D+00  * ( r - a ) * ( corner + ( s - d ) * dcds )
  dtds(5) =  - 13.5D+00 * ( r - d ) * ( 3.0D+00 * s * s &
    - 2.0D+00 * ( a + c + d ) * s + a * c + c * d + d * a )
  dtds(6) =  13.5D+00 * ( r - a ) * ( 3.0D+00 * s * s &
    - 2.0D+00 * ( a + c + d ) * s + a * c + c * d + d * a )
  dtds(7) =  13.5D+00 * ( r - d ) * ( 3.0D+00 * s * s &
    - 2.0D+00 * ( a + b + d ) * s + a * b + b * d + d * a )
  dtds(8) =  - 13.5D+00 * ( r - a ) * ( 3.0D+00 * s * s &
    - 2.0D+00 * ( a + b + d ) * s + a * b + b * d + d * a )
  dtds(9) =  - 0.125D+00 * ( r - d ) * ( corner + ( s - a ) * dcds )
  dtds(10) = 13.5D+00 * ( r - a ) * ( r - c ) * ( r - d ) 
  dtds(11) = - 13.5D+00 * ( r - a ) * ( r - b ) * ( r - d ) 
  dtds(12) = 0.125D+00 * ( r - a ) * ( corner + ( s - a ) * dcds )

  return
end
subroutine shape_q16 ( r, s, t, dtdr, dtds )

!*****************************************************************************80
!
!! SHAPE_Q16 evaluates shape functions for a 16 node reference quadrilateral.
!
!  Reference Element Q16:
!
!    |
!    1 13--14--15--16
!    |  |   :   :   |
!    |  |   :   :   |
!    |  9..10..11..12
!    S  |   :   :   |
!    |  |   :   :   |
!    |  5...6...7...8
!    |  |   :   :   |
!    |  |   :   :   |  
!    0  1---2---3---4
!    |
!    +--0-----R-----1-->
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, S, the reference coordinates of a point.
!
!    Output, real ( kind = 8 ) T(16), the basis functions at the point.
!
!    Output, real ( kind = 8 ) DTDR(16), the R basis derivatives at the point.
!
!    Output, real ( kind = 8 ) DTDS(16), the S basis derivatives at the point.
!
  implicit none

  real    ( kind = 8 ) dabc
  real    ( kind = 8 ) dabd
  real    ( kind = 8 ) dacd
  real    ( kind = 8 ) dbcd
  real    ( kind = 8 ) dtdr(16)
  real    ( kind = 8 ) dtds(16)
  real    ( kind = 8 ) r
  real    ( kind = 8 ) ra
  real    ( kind = 8 ) rb
  real    ( kind = 8 ) rc
  real    ( kind = 8 ) rd
  real    ( kind = 8 ) s
  real    ( kind = 8 ) sa
  real    ( kind = 8 ) sb
  real    ( kind = 8 ) sc
  real    ( kind = 8 ) sd
  real    ( kind = 8 ) t(16)

  ra = r - 0.0D+00
  rb = r - 1.0D+00 / 3.0D+00
  rc = r - 2.0D+00 / 3.0D+00
  rd = r - 1.0D+00

  sa = s - 0.0D+00
  sb = s - 1.0D+00 / 3.0D+00
  sc = s - 2.0D+00 / 3.0D+00
  sd = s - 1.0D+00

  t(1)  =   (  81.0D+00 / 4.0D+00 ) * rb * rc * rd * sb * sc * sd
  t(2)  = - ( 243.0D+00 / 4.0D+00 ) * ra * rc * rd * sb * sc * sd
  t(3)  =   ( 243.0D+00 / 4.0D+00 ) * ra * rb * rd * sb * sc * sd
  t(4)  = - (  81.0D+00 / 4.0D+00 ) * ra * rb * rc * sb * sc * sd

  t(5)  = - ( 243.0D+00 / 4.0D+00 ) * rb * rc * rd * sa * sc * sd
  t(6)  =   ( 729.0D+00 / 4.0D+00 ) * ra * rc * rd * sa * sc * sd
  t(7)  = - ( 729.0D+00 / 4.0D+00 ) * ra * rb * rd * sa * sc * sd
  t(8)  =   ( 243.0D+00 / 4.0D+00 ) * ra * rb * rc * sa * sc * sd

  t(9)  =   ( 243.0D+00 / 4.0D+00 ) * rb * rc * rd * sa * sb * sd
  t(10) = - ( 729.0D+00 / 4.0D+00 ) * ra * rc * rd * sa * sb * sd
  t(11) =   ( 729.0D+00 / 4.0D+00 ) * ra * rb * rd * sa * sb * sd
  t(12) = - ( 243.0D+00 / 4.0D+00 ) * ra * rb * rc * sa * sb * sd

  t(13) = - (  81.0D+00 / 4.0D+00 ) * rb * rc * rd * sa * sb * sc
  t(14) =   ( 243.0D+00 / 4.0D+00 ) * ra * rc * rd * sa * sb * sc
  t(15) = - ( 243.0D+00 / 4.0D+00 ) * ra * rb * rd * sa * sb * sc
  t(16) =   (  81.0D+00 / 4.0D+00 ) * ra * rb * rc * sa * sb * sc

  dbcd = 3.0D+00 * r * r -  4.0D+00 * r       + 11.0D+00 / 9.0D+00
  dacd = 3.0D+00 * r * r - 10.0D+00 * r / 3.0D+00 +  2.0D+00 / 3.0D+00
  dabd = 3.0D+00 * r * r -  8.0D+00 * r / 3.0D+00 +  1.0D+00 / 3.0D+00
  dabc = 3.0D+00 * r * r -  2.0D+00 * r       +  2.0D+00 / 9.0D+00

  dtdr(1)  =   (  81.0D+00 / 4.0D+00 ) * dbcd * sb * sc * sd
  dtdr(2)  = - ( 243.0D+00 / 4.0D+00 ) * dacd * sb * sc * sd
  dtdr(3)  =   ( 243.0D+00 / 4.0D+00 ) * dabd * sb * sc * sd
  dtdr(4)  = - (  81.0D+00 / 4.0D+00 ) * dabc * sb * sc * sd
  dtdr(5)  = - ( 243.0D+00 / 4.0D+00 ) * dbcd * sa * sc * sd
  dtdr(6)  =   ( 729.0D+00 / 4.0D+00 ) * dacd * sa * sc * sd
  dtdr(7)  = - ( 729.0D+00 / 4.0D+00 ) * dabd * sa * sc * sd
  dtdr(8)  =   ( 243.0D+00 / 4.0D+00 ) * dabc * sa * sc * sd
  dtdr(9)  =   ( 243.0D+00 / 4.0D+00 ) * dbcd * sa * sb * sd
  dtdr(10) = - ( 729.0D+00 / 4.0D+00 ) * dacd * sa * sb * sd
  dtdr(11) =   ( 729.0D+00 / 4.0D+00 ) * dabd * sa * sb * sd
  dtdr(12) = - ( 243.0D+00 / 4.0D+00 ) * dabc * sa * sb * sd
  dtdr(13) = - (  81.0D+00 / 4.0D+00 ) * dbcd * sa * sb * sc
  dtdr(14) =   ( 243.0D+00 / 4.0D+00 ) * dacd * sa * sb * sc
  dtdr(15) = - ( 243.0D+00 / 4.0D+00 ) * dabd * sa * sb * sc
  dtdr(16) =   (  81.0D+00 / 4.0D+00 ) * dabc * sa * sb * sc

  dbcd = 3.0D+00 * s * s -  4.0D+00 * s       + 11.0D+00 / 9.0D+00
  dacd = 3.0D+00 * s * s - 10.0D+00 * s / 3.0D+00 +  2.0D+00 / 3.0D+00
  dabd = 3.0D+00 * s * s -  8.0D+00 * s / 3.0D+00 +  1.0D+00 / 3.0D+00
  dabc = 3.0D+00 * s * s -  2.0D+00 * s       +  2.0D+00 / 9.0D+00

  dtds(1)  =   (  81.0D+00 / 4.0D+00 ) * rb * rc * rd * dbcd
  dtds(2)  = - ( 243.0D+00 / 4.0D+00 ) * ra * rc * rd * dbcd
  dtds(3)  =   ( 243.0D+00 / 4.0D+00 ) * ra * rb * rd * dbcd
  dtds(4)  = - (  81.0D+00 / 4.0D+00 ) * ra * rb * rc * dbcd
  dtds(5)  = - ( 243.0D+00 / 4.0D+00 ) * rb * rc * rd * dacd
  dtds(6)  =   ( 729.0D+00 / 4.0D+00 ) * ra * rc * rd * dacd
  dtds(7)  = - ( 729.0D+00 / 4.0D+00 ) * ra * rb * rd * dacd
  dtds(8)  =   ( 243.0D+00 / 4.0D+00 ) * ra * rb * rc * dacd
  dtds(9)  =   ( 243.0D+00 / 4.0D+00 ) * rb * rc * rd * dabd
  dtds(10) = - ( 729.0D+00 / 4.0D+00 ) * ra * rc * rd * dabd
  dtds(11) =   ( 729.0D+00 / 4.0D+00 ) * ra * rb * rd * dabd
  dtds(12) = - ( 243.0D+00 / 4.0D+00 ) * ra * rb * rc * dabd
  dtds(13) = - (  81.0D+00 / 4.0D+00 ) * rb * rc * rd * dabc
  dtds(14) =   ( 243.0D+00 / 4.0D+00 ) * ra * rc * rd * dabc
  dtds(15) = - ( 243.0D+00 / 4.0D+00 ) * ra * rb * rd * dabc
  dtds(16) =   (  81.0D+00 / 4.0D+00 ) * ra * rb * rc * dabc
  
  return
end
subroutine shape_ql ( r, s, t, dtdr, dtds )

!*****************************************************************************80
!
!! SHAPE_QL evaluates shape functions for a 6 node quadratic/linear.
!
!  Reference Element QL:
!
!    |
!    1  4--5--6
!    |  |     |
!    |  |     |
!    S  |     |
!    |  |     |
!    |  |     |
!    0  1--2--3
!    |
!    +--0--R--1-->
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, S, the reference coordinates of a point.
!
!    Output, real ( kind = 8 ) T(6), the basis functions at the point.
!
!    Output, real ( kind = 8 ) DTDR(6), the R basis derivatives at the point.
!
!    Output, real ( kind = 8 ) DTDS(6), the S basis derivatives at the point.
!
  implicit none

  real    ( kind = 8 ) dtdr(6)
  real    ( kind = 8 ) dtds(6)
  real    ( kind = 8 ) r
  real    ( kind = 8 ) s
  real    ( kind = 8 ) t(6)

  t(1) = - 2.0D+00 *     ( r - 0.5D+00 ) * ( r - 1.0D+00 )     * ( s - 1.0D+00 )
  t(2) =   4.0D+00 * r                   * ( r - 1.0D+00 )     * ( s - 1.0D+00 )
  t(3) = - 2.0D+00 * r * ( r - 0.5D+00 )                       * ( s - 1.0D+00 )
  t(4) =   2.0D+00 *     ( r - 0.5D+00 ) * ( r - 1.0D+00 ) * s
  t(5) = - 4.0D+00 * r                   * ( r - 1.0D+00 ) * s
  t(6) =   2.0D+00 * r * ( r - 0.5D+00 )                   * s

  dtdr(1) = 2.0D+00 * ( - 2.0D+00 * r + 1.5D+00 )     * ( s - 1.0D+00 )
  dtdr(2) = 4.0D+00 * (   2.0D+00 * r - 1.0D+00 )     * ( s - 1.0D+00 )
  dtdr(3) = 2.0D+00 * ( - 2.0D+00 * r + 0.5D+00 )     * ( s - 1.0D+00 ) 
  dtdr(4) = 2.0D+00 * (   2.0D+00 * r - 1.5D+00 ) * s
  dtdr(5) = 4.0D+00 * ( - 2.0D+00 * r + 1.0D+00 ) * s
  dtdr(6) = 2.0D+00 * (   2.0D+00 * r - 0.5D+00 ) * s

  dtds(1) = - 2.0D+00 *     ( r - 0.5D+00 ) * ( r - 1.0D+00 )
  dtds(2) =   4.0D+00 * r                   * ( r - 1.0D+00 )
  dtds(3) = - 2.0D+00 * r * ( r - 0.5D+00 )
  dtds(4) =   2.0D+00 *     ( r - 0.5D+00 ) * ( r - 1.0D+00 )
  dtds(5) = - 4.0D+00 * r                   * ( r - 1.0D+00 )
  dtds(6) =   2.0D+00 * r * ( r - 0.5D+00 )

  return
end
subroutine shape_t3 ( r, s, t, dtdr, dtds )

!*****************************************************************************80
!
!! SHAPE_T3 evaluates shape functions for a 3 node reference triangle.
!
!  Reference Element T3:
!
!    |
!    1  3
!    |  |\
!    |  | \
!    S  |  \
!    |  |   \
!    |  |    \
!    0  1-----2
!    |
!    +--0--R--1-->
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, S, the reference coordinates of a point.
!
!    Output, real ( kind = 8 ) T(3), the basis functions at the point.
!
!    Output, real ( kind = 8 ) DTDR(3), the R basis derivatives at the point.
!
!    Output, real ( kind = 8 ) DTDS(3), the S basis derivatives at the point.
!
  implicit none

  real    ( kind = 8 ) dtdr(3)
  real    ( kind = 8 ) dtds(3)
  real    ( kind = 8 ) r
  real    ( kind = 8 ) s
  real    ( kind = 8 ) t(3)

  t(1) = 1.0D+00 - r - s
  t(2) =           r
  t(3) =               s

  dtdr(1) = -1.0D+00
  dtdr(2) =  1.0D+00
  dtdr(3) =  0.0D+00

  dtds(1) = -1.0D+00
  dtds(2) =  0.0D+00
  dtds(3) =  1.0D+00

  return
end
subroutine shape_t4 ( r, s, t, dtdr, dtds )

!*****************************************************************************80
!
!! SHAPE_T4 evaluates shape functions for a 4 node reference triangle.
!
!  Reference Element T4:
!
!    |
!    1  3
!    |  |\
!    |  | \
!    S  |  \
!    |  | 4 \
!    |  |    \
!    0  1-----2
!    |
!    +--0--R--1-->
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, S, the reference coordinates of a point.
!
!    Output, real ( kind = 8 ) T(4), the basis functions at the point.
!
!    Output, real ( kind = 8 ) DTDR(4), the R basis derivatives at the point.
!
!    Output, real ( kind = 8 ) DTDS(4), the S basis derivatives at the point.
!
  implicit none

  real    ( kind = 8 ) dtdr(4)
  real    ( kind = 8 ) dtds(4)
  real    ( kind = 8 ) r
  real    ( kind = 8 ) s
  real    ( kind = 8 ) t(4)

  t(1) = ( 1.0D+00 - 9.0D+00 * r * s ) * ( 1.0D+00 - r - s )
  t(2) = r * ( 1.0D+00 - 9.0D+00 * ( 1.0D+00 - r - s ) * s )
  t(3) = s * ( 1.0D+00 - 9.0D+00 * ( 1.0D+00 - r - s ) * r )
  t(4) = 27.0D+00 * ( 1.0D+00 - r - s ) * r * s

  dtdr(1) = -1.0D+00 +  9.0D+00 * ( - s + 2.0D+00 * r * s + s**2 )
  dtdr(2) =  1.0D+00 +  9.0D+00 * ( - s + 2.0D+00 * r * s + s**2 )
  dtdr(3) =             9.0D+00 * ( - s + 2.0D+00 * r * s + s**2 )
  dtdr(4) =          - 27.0D+00 * ( - s + 2.0D+00 * r * s + s**2 )

  dtds(1) = -1.0D+00 +  9.0D+00 * ( - r + r**2 + 2.0D+00 * r * s )
  dtds(2) =             9.0D+00 * ( - r + r**2 + 2.0D+00 * r * s )
  dtds(3) =  1.0D+00 +  9.0D+00 * ( - r + r**2 + 2.0D+00 * r * s )
  dtds(4) =          - 27.0D+00 * ( - r + r**2 + 2.0D+00 * r * s )

  return
end
subroutine shape_t6 ( r, s, t, dtdr, dtds )

!*****************************************************************************80
!
!! SHAPE_T6 evaluates shape functions for a 6 node reference triangle.
!
!  Reference Element T6:
!
!    |
!    1  3
!    |  |\
!    |  | \
!    S  6  5
!    |  |   \
!    |  |    \
!    0  1--4--2
!    |
!    +--0--R--1-->
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 April 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, S, the reference coordinates of a point.
!
!    Output, real ( kind = 8 ) T(6), the basis functions at the point.
!
!    Output, real ( kind = 8 ) DTDR(6), the R basis derivatives at the point.
!
!    Output, real ( kind = 8 ) DTDS(6), the S basis derivatives at the point.
!
  implicit none

  real    ( kind = 8 ) dtdr(6)
  real    ( kind = 8 ) dtds(6)
  real    ( kind = 8 ) r
  real    ( kind = 8 ) s
  real    ( kind = 8 ) t(6)

  t(1) = 2.0D+00 *     ( 1.0D+00 - r - s ) * ( 0.5D+00 - r - s )
  t(2) = 2.0D+00 * r * ( r - 0.5D+00 )
  t(3) = 2.0D+00 * s * ( s - 0.5D+00 )
  t(4) = 4.0D+00 * r * ( 1.0D+00 - r - s )
  t(5) = 4.0D+00 * r * s
  t(6) = 4.0D+00 * s * ( 1.0D+00 - r - s )

  dtdr(1) = - 3.0D+00 + 4.0D+00 * r + 4.0D+00 * s
  dtdr(2) = - 1.0D+00 + 4.0D+00 * r
  dtdr(3) =   0.0D+00
  dtdr(4) =   4.0D+00 - 8.0D+00 * r - 4.0D+00 * s
  dtdr(5) =                           4.0D+00 * s
  dtdr(6) =                         - 4.0D+00 * s

  dtds(1) = - 3.0D+00 + 4.0D+00 * r + 4.0D+00 * s
  dtds(2) =   0.0D+00
  dtds(3) = - 1.0D+00               + 4.0D+00 * s
  dtds(4) =           - 4.0D+00 * r
  dtds(5) =             4.0D+00 * r
  dtds(6) =   4.0D+00 - 4.0D+00 * r - 8.0D+00 * s

  return
end
subroutine shape_t10 ( r, s, t, dtdr, dtds )

!*****************************************************************************80
!
!! SHAPE_T10 evaluates shape functions for a 10 node reference triangle.
!
!  Reference Element T10:
!
!    |
!    1  10
!    |  |\
!    |  | \
!    |  8  9
!    |  |   \
!    S  |    \
!    |  5  6  7
!    |  |      \
!    |  |       \
!    0  1--2--3--4
!    |
!    +--0----R---1-->
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, S, the reference coordinates of a point.
!
!    Output, real ( kind = 8 ) T(10), the basis functions at the point.
!
!    Output, real ( kind = 8 ) DTDR(10), the R basis derivatives at the point.
!
!    Output, real ( kind = 8 ) DTDS(10), the S basis derivatives at the point.
!
  implicit none

  real    ( kind = 8 ) a 
  real    ( kind = 8 ) b
  real    ( kind = 8 ) c
  real    ( kind = 8 ) dtdr(10)
  real    ( kind = 8 ) dtds(10)
  real    ( kind = 8 ) r
  real    ( kind = 8 ) s
  real    ( kind = 8 ) t(10)

  a = 1.0D+00 / 3.0D+00
  b = 2.0D+00 / 3.0D+00
  c = 1.0D+00

  t(1)  = 4.5D+00 * ( a - r - s ) * ( b - r - s ) * ( c - r - s )
  t(2)  = 13.5D+00 * r * ( b - r - s ) * ( c - r - s ) 
  t(3)  = - 13.5D+00 * r * ( a - r ) * ( c - r - s )
  t(4)  = 4.5D+00 * r * ( a - r ) * ( b - r )
  t(5)  = 13.5D+00 * s * ( b - r - s ) * ( c - r - s )
  t(6)  = 27.0D+00 * r * s * ( c - r - s )
  t(7)  = 13.5D+00 * r * s * ( r - a )
  t(8)  = 13.5D+00 * s * ( s - a ) * ( c - r - s )
  t(9)  = 13.5D+00 * r * s * ( s - a )
  t(10) = 4.5D+00 * s * ( a - s ) * ( b - s )

  dtdr(1) = 4.5D+00 * ( ( a - s ) * ( 2.0D+00 * r - c - b + 2.0D+00 * s ) &
    - ( s - b ) * ( s - c ) - 2.0D+00 * ( 2.0D+00 * s - b - c ) * r &
    - 3.0D+00 * r * r )
  dtdr(2) = 13.5D+00 * ( &
    ( s - b ) * ( s - c ) + 2.0D+00 * ( 2.0D+00 * s - b - c ) * r &
    + 3.0D+00 * r * r )
  dtdr(3) = - 13.5D+00 * ( a * ( c - s ) + 2.0D+00 * ( s - a - c ) * r &
    + 3.0D+00 * r * r )
  dtdr(4) = 4.5D+00 * ( a * b - 2.0D+00 * ( a + b ) * r + 3.0D+00 * r * r )
  dtdr(5) = 13.5D+00 * s * ( 2.0D+00 * s - b - c + 2.0D+00 * r )
  dtdr(6) = 27.0D+00 * s * ( c - s - 2.0D+00 * r )
  dtdr(7) = 13.5D+00 * s * ( 2.0D+00 * r - a )
  dtdr(8) = - 13.5D+00 * s * ( s - a )
  dtdr(9) = 13.5D+00 * s * ( s - a)
  dtdr(10) = 0.0D+00

  dtds(1) = 4.5D+00 * ( ( a - r ) * ( 2.0D+00 * s - c - b + 2.0D+00 * r ) &
    - ( r - b ) * ( r - c ) - 2.0D+00 * ( 2.0D+00 * r - b - c ) * s &
    - 3.0D+00 * s * s )
  dtds(2) = 13.5D+00 * r * ( 2.0D+00 * s + 2.0D+00 * r - b - c )
  dtds(3) = 13.5D+00 * r * ( a - r )
  dtds(4) = 0.0D+00
  dtds(5) = 13.5D+00 * ( ( r - b ) * ( r - c ) + &
    2.0D+00 * ( 2.0D+00 * r - b - c ) * s + 3.0D+00 * s * s )
  dtds(6) = 27.0D+00 * r * ( c - r - 2.0D+00 * s )
  dtds(7) = 13.5D+00 * r * ( r - a )
  dtds(8) = - 13.5D+00 * ( a * ( c - r ) + 2.0D+00 * ( r - c - a ) * s &
    + 3.0D+00 * s * s )
  dtds(9) = 13.5D+00 * r * ( 2.0D+00 * s - a)
  dtds(10) = 4.5D+00 * ( a * b - 2.0D+00 * ( a + b ) * s + 3.0D+00 * s * s )

  return
end
subroutine shape_test ( code )

!*****************************************************************************80
!
!! SHAPE_TEST verifies the shape function values at the basis nodes.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) CODE, identifies the element to be used.
!    Legal values include 'Q4', 'Q8', 'Q9', 'Q12', 'Q16', 'QL', 
!    'T3', 'T4', 'T6' and 'T10'.
!
  implicit none

  real ( kind = 8 ) area
  character ( len = * )  code
  real ( kind = 8 ), allocatable, dimension ( : ) :: dtdr
  real ( kind = 8 ), allocatable, dimension ( : ) :: dtds
  integer   ( kind = 4 ) element_order
  integer   ( kind = 4 ) node
  integer   ( kind = 4 ) order_code
  real ( kind = 8 ), allocatable, dimension ( : ) :: r
  real ( kind = 8 ) rsum
  real ( kind = 8 ), allocatable, dimension ( : ) :: s
  real ( kind = 8 ) ssum
  real ( kind = 8 ), allocatable, dimension ( : ) :: t

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SHAPE_TEST for "' // trim ( code ) &
    // '" shape functions.'

  element_order = order_code ( code )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Element order = ', element_order

  allocate ( dtdr(1:element_order) )
  allocate ( dtds(1:element_order) )
  allocate ( r(1:element_order) )
  allocate ( s(1:element_order) )
  allocate ( t(1:element_order) )

  call node_reference ( code, r, s, area )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Basis function values at basis nodes'
  write ( *, '(a)' ) '  should form the identity matrix.'
  write ( *, '(a)' ) ' '

  do node = 1, element_order
    call shape ( code, r(node), s(node), t, dtdr, dtds )
    write ( *, '(2x,10f7.3)' ) t(1:element_order)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The R and S derivatives should sum to 0.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '        dTdR sum      dTdS sum'
  write ( *, '(a)' ) ' '
  do node = 1, element_order
    call shape ( code, r(node), s(node), t, dtdr, dtds )
    rsum = sum ( dtdr(1:element_order) )
    ssum = sum ( dtds(1:element_order) )
    write ( *, '(2x,f14.8,f14.8)' ) rsum, ssum
  end do

  deallocate ( dtdr )
  deallocate ( dtds )
  deallocate ( r )
  deallocate ( s )
  deallocate ( t )

  return
end
subroutine sphere_grid_element_num ( code, nelemx, nelemy, element_num )

!*****************************************************************************80
!
!! SPHERE_GRID_ELEMENT_NUM returns the number of elements in a sphere grid.
!
!  Discussion:
!
!    The number of elements generated will be NELEMX * NELEMY for
!    quadrilaterals, or 2 * NELEMX * ( NELEMY - 1 ) for triangles.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) CODE, identifies the element desired.
!    Legal values include 'Q4', 'Q9', 'Q16', 'T3', 'T6'.
!
!    Input, integer ( kind = 4 ) NELEMX, NELEMY, the number of quadrilaterals 
!    along the X and Y directions.  
!
!    Output, integer ( kind = 4 ) ELEMENT_NUM, the number of elements in 
!    the grid.
!
  implicit none

  character ( len = * )  code
  integer   ( kind = 4 ) element_num
  integer   ( kind = 4 ) nelemx
  integer   ( kind = 4 ) nelemy
  logical                s_eqi

  if ( s_eqi ( code, 'Q4' ) ) then
    call sphere_grid_q4_element_num ( nelemx, nelemy, element_num )
  else if ( s_eqi ( code, 'Q9' ) ) then
    call sphere_grid_q9_element_num ( nelemx, nelemy, element_num )
  else if ( s_eqi ( code, 'Q16' ) ) then
    call sphere_grid_q16_element_num ( nelemx, nelemy, element_num )
  else if ( s_eqi ( code, 'T3' ) ) then
    call sphere_grid_t3_element_num ( nelemx, nelemy, element_num )
  else if ( s_eqi ( code, 'T6' ) ) then
    call sphere_grid_t6_element_num ( nelemx, nelemy, element_num )
  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPHERE_GRID_ELEMENT_NUM - Fatal error!'
    write ( *, '(a)' ) '  Illegal value of CODE = "' // trim ( code ) // '".'
    element_num = -1
    stop

  end if

  return
end
subroutine sphere_grid_node_num ( code, nelemx, nelemy, node_num )

!*****************************************************************************80
!
!! SPHERE_GRID_NODE_NUM returns the number of nodes in a sphere grid.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) CODE, identifies the element desired.
!    Legal values include 'Q4', 'Q9', 'Q16', 'T3', 'T6'.
!
!    Input, integer ( kind = 4 ) NELEMX, NELEMY, the number of quadrilaterals 
!    along the X and Y directions.  
!
!    Output, integer ( kind = 4 ) NODE_NUM, the number of elements in the grid.
!
  implicit none

  character ( len = * )  code
  integer   ( kind = 4 ) node_num
  integer   ( kind = 4 ) nelemx
  integer   ( kind = 4 ) nelemy
  logical                s_eqi

  if ( s_eqi ( code, 'Q4' ) ) then
    call sphere_grid_q4_node_num ( nelemx, nelemy, node_num )
  else if ( s_eqi ( code, 'Q9' ) ) then
    call sphere_grid_q9_node_num ( nelemx, nelemy, node_num )
  else if ( s_eqi ( code, 'Q16' ) ) then
    call sphere_grid_q16_node_num ( nelemx, nelemy, node_num )
  else if ( s_eqi ( code, 'T3' ) ) then
    call sphere_grid_t3_node_num ( nelemx, nelemy, node_num )
  else if ( s_eqi ( code, 'T6' ) ) then
    call sphere_grid_t6_node_num ( nelemx, nelemy, node_num )
  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPHERE_GRID_NODE_NUM - Fatal error!'
    write ( *, '(a)' ) '  Illegal value of CODE = "' // trim ( code ) // '".'
    node_num = -1
    stop

  end if

  return
end
subroutine sphere_grid_q4_element ( nelemx, nelemy, element_node )

!*****************************************************************************80
!
!! SPHERE_GRID_Q4_ELEMENT produces a Q4 sphere grid.
!
!  Discussion:
!
!    This would be the same as the grid in a plane, except that all the
!    nodes along the bottom edge are identified (replaced by a single node
!    that is the south pole) and similarly for the top edge, and the
!    nodes on the extreme right edge are identified pairwise with those 
!    on the extreme left edge.
!
!  Example:
!
!    Input:
!
!      NELEMX = 3, NELEMY = 4
!
!    Output:
!
!      ELEMENT_NODE =
!         1,  1,  3,  2;
!         1,  1,  4,  3;
!         1,  1,  2,  4;
!         2,  3,  6,  5;
!         3,  4,  7,  6;
!         4,  2,  5,  7;
!         5,  6,  9,  8;
!         6,  7, 10,  9;
!         7,  5,  8, 10;
!         8,  9, 11, 11;
!         9, 10, 11, 11;
!        10,  8, 11, 11;
!
!  Grid:
!
!   11----11----11----11
!    |     |     |     |
!    | E10 | E11 | E12 |
!    |     |     |     |
!    8-----9----10-----8
!    |     |     |     |
!    | E7  | E8  | E9  |
!    |     |     |     |
!    5-----6-----7-----5
!    |     |     |     |
!    | E4  | E5  | E6  |
!    |     |     |     |
!    2-----3-----4-----2
!    |     |     |     |
!    | E1  | E2  | E3  |
!    |     |     |     |
!    1-----1-----1-----1
!
!  Reference Element Q4:
!
!    |
!    1  4------3
!    |  |      |
!    S  |      |
!    |  |      |
!    |  |      |
!    0  1------2
!    |
!    +--0--R---1-->
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NELEMX, NELEMY, the number of elements along 
!    the X and Y directions.  The number of elements generated will be
!    NELEMX * NELEMY.
!
!    Output, integer ( kind = 4 ) ELEMENT_NODE(4,NELEMX*NELEMY), the nodes that form
!    each element.
!
  implicit none

  integer ( kind = 4 ) nelemx
  integer ( kind = 4 ) nelemy
  integer ( kind = 4 ), parameter :: element_order = 4

  integer ( kind = 4 ) base1
  integer ( kind = 4 ) base2
  integer ( kind = 4 ) element
  integer ( kind = 4 ) element_node(element_order,nelemx*nelemy)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  element = 0

  do j = 1, nelemy

    base1 = ( j - 1 ) * nelemx + 2 - nelemx

    do i = 1, nelemx

      base2 = base1 + i - 1

      element = element + 1

      element_node(1,element) = base2
      if ( i < nelemx ) then
        element_node(2,element) = base2 + 1
      else
        element_node(2,element) = base1
      end if
      element_node(3,element) = element_node(2,element) + nelemx
      element_node(4,element) = element_node(1,element) + nelemx

      if ( j == 1 ) then

        element_node( 1,element) = 1
        element_node( 2,element) = 1

      else if ( j == nelemy ) then

        element_node(3,element) = base1 + nelemx
        element_node(4,element) = base1 + nelemx

      end if

    end do
  end do

  return
end
subroutine sphere_grid_q4_element_num ( nelemx, nelemy, element_num )

!*****************************************************************************80
!
!! SPHERE_GRID_Q4_ELEMENT_NUM counts the elements in a Q4 sphere grid.
!
!  Example:
!
!    Input:
!
!      NELEMX = 3, NELEMY = 2
!
!    Output:
!
!      ELEMENT_NUM = NELEMX * NELEMY = 6
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NELEMX, NELEMY, the number of elements along 
!    the X and Y directions. 
!
!    Output, integer ( kind = 4 ) ELEMENT_NUM, the number of elements in 
!    the grid.
!
  implicit none

  integer ( kind = 4 ) element_num
  integer ( kind = 4 ) nelemx
  integer ( kind = 4 ) nelemy

  element_num = nelemx * nelemy

  return
end
subroutine sphere_grid_q4_node_num ( nelemx, nelemy, node_num )

!*****************************************************************************80
!
!! SPHERE_GRID_Q4_NODE_NUM counts nodes in a Q4 sphere grid.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NELEMX, NELEMY, the number of elements along 
!    the X and Y directions.
!
!    Output, integer ( kind = 4 ) NODE_NUM, the number of nodes in the grid.
!
  implicit none

  integer ( kind = 4 ) nelemx
  integer ( kind = 4 ) nelemy
  integer ( kind = 4 ) node_num

  node_num = nelemx * ( nelemy - 1 ) + 2

  return
end
subroutine sphere_grid_q4_node_xyz ( nelemx, nelemy, node_xyz )

!*****************************************************************************80
!
!! SPHERE_GRID_Q4_NODE_XYZ produces node coordinates for a Q4 sphere grid.
!
!  Discussion:
!
!    The number of nodes to be generated is
!
!      NODE_NUM = NELEMX * ( NELEMY - 1 ) + 2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NELEMX, NELEMY, the number of elements along 
!    the X and Y directions.
!
!    Output, real ( kind = 8 ) NODE_XYZ(3,NODE_NUM), 
!    the node coordinates.
!
  implicit none

  integer ( kind = 4 ) nelemx
  integer ( kind = 4 ) nelemy

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) node
  real    ( kind = 8 ) node_xyz(3,nelemx*(nelemy-1)+2)
  real    ( kind = 8 ) phi
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real    ( kind = 8 ) theta

  node = 0

  node = node + 1      
  node_xyz(1,node) =  0.0D+00
  node_xyz(2,node) =  0.0D+00
  node_xyz(3,node) = -1.0D+00

  do j = nelemy, 2, -1

    phi = real ( j - 1, kind = 8 ) * pi &
        / real ( nelemy, kind = 8 )

    do i = 1, nelemx

      theta = real ( i - 1, kind = 8 ) * 2.0D+00 * pi &
            / real ( nelemx, kind = 8 )

      node = node + 1      
      node_xyz(1,node) = cos ( theta ) * sin ( phi )
      node_xyz(2,node) = sin ( theta ) * sin ( phi )
      node_xyz(3,node) =                 cos ( phi )

    end do
  end do

  node = node + 1      
  node_xyz(1,node) =  0.0D+00
  node_xyz(2,node) =  0.0D+00
  node_xyz(3,node) =  1.0D+00

  return
end
subroutine sphere_grid_q9_element ( nelemx, nelemy, element_node )

!*****************************************************************************80
!
!! SPHERE_GRID_Q9_ELEMENT produces a Q9 sphere grid.
!
!  Discussion:
!
!    This would be the same as the grid in a plane, except that all the
!    nodes along the bottom edge are identified (replaced by a single node
!    that is the south pole) and similarly for the top edge, and the
!    nodes on the extreme right edge are identified pairwise with those 
!    on the extreme left edge.
!
!  Example:
!
!    Input:
!
!      NELEMX = 3, NELEMY = 4
!
!    Output:
!
!      ELEMENT_NODE =
!         1,  1, 10,  8,  1,  4,  9,  2,  3;
!         1,  1, 12, 10,  1,  6, 11,  4,  5;
!         1,  1,  8, 12,  1,  2, 13,  6,  7;
!         8, 10, 22, 20,  9, 16, 21, 14, 15;
!        10, 12, 24, 22, 11, 18, 23, 16, 17;
!        12,  8, 20, 24, 13, 14, 25, 18, 19;
!        20, 22, 34, 32, 21, 28, 33, 26, 27;
!        22, 24, 36, 34, 23, 30, 35, 28, 29;
!        24, 20, 32, 36, 25, 26, 37, 30, 31;
!        32, 34, 44, 44, 33, 40, 44, 38, 39;
!        34, 36, 44, 44, 35, 42, 44, 40, 41;
!        36, 32, 44, 44, 37, 38, 44, 42, 43;
!
!  Grid:
!
!   44-44-44-44-44-44-44
!    |     |     |     |
!   38 39 40 41 42 43 38
!    |     |     |     |
!   32-33-34-35-36-37-32
!    |     |     |     |
!   26 27 28 29 30 31 26
!    |     |     |     |
!   20-21-22-23-24-25-20
!    |     |     |     |
!   14 15 16 17 18 19 14
!    |     |     |     |
!    8--9-10-11-12-13--8
!    |     |     |     |
!    2  3  4  5  6  7  2
!    |     |     |     |
!    1--1--1--1--1--1--1
!
!  Reference Element Q9:
!
!    |
!    1  4--7--3
!    |  |     |
!    |  |     |
!    S  8  9  6
!    |  |     |
!    |  |     |
!    0  1--5--2
!    |
!    +--0--R--1-->
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NELEMX, NELEMY, the number of elements along 
!    the X and Y directions.
!
!    Output, integer ( kind = 4 ) ELEMENT_NODE(9,NELEMX*NELEMY), 
!    the nodes that form each element.
!
  implicit none

  integer ( kind = 4 ) nelemx
  integer ( kind = 4 ) nelemy
  integer ( kind = 4 ), parameter :: element_order = 9

  integer ( kind = 4 ) base1
  integer ( kind = 4 ) base2
  integer ( kind = 4 ) element
  integer ( kind = 4 ) element_node(element_order,4*nelemx*nelemy-2*nelemx+2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  element = 0

  do j = 1, nelemy

    base1 = ( j - 1 ) * 2 * ( 2 * nelemx ) + 2 - 2 * nelemx

    do i = 1, nelemx

      base2 = base1 + 2 * ( i - 1 )

      element = element + 1

      element_node(1,element) = base2
      element_node(5,element) = base2 + 1

      if ( i < nelemx ) then
        element_node(2,element) = base2 + 2
      else
        element_node(2,element) = base1
      end if

      element_node(8,element) = element_node(1,element) + 2 * nelemx
      element_node(9,element) = element_node(5,element) + 2 * nelemx
      element_node(6,element) = element_node(2,element) + 2 * nelemx

      element_node(4,element) = element_node(8,element) + 2 * nelemx
      element_node(7,element) = element_node(9,element) + 2 * nelemx
      element_node(3,element) = element_node(6,element) + 2 * nelemx

      if ( j == 1 ) then

        element_node(1,element) = 1
        element_node(5,element) = 1
        element_node(2,element) = 1

      else if ( j == nelemy ) then

        element_node(4,element) = base1 + 4 * nelemx
        element_node(7,element) = base1 + 4 * nelemx
        element_node(3,element) = base1 + 4 * nelemx

      end if

    end do

  end do

  return
end
subroutine sphere_grid_q9_element_num ( nelemx, nelemy, element_num )

!*****************************************************************************80
!
!! SPHERE_GRID_Q9_ELEMENT_NUM counts the elements in a Q9 sphere grid.
!
!  Example:
!
!    Input:
!
!      NELEMX = 3, NELEMY = 2
!
!    Output:
!
!      ELEMENT_NUM = NELEMX * NELEMY = 6
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NELEMX, NELEMY, the number of elements along 
!    the X and Y directions. 
!
!    Output, integer ( kind = 4 ) ELEMENT_NUM, the number of elements in 
!    the grid.
!
  implicit none

  integer ( kind = 4 ) element_num
  integer ( kind = 4 ) nelemx
  integer ( kind = 4 ) nelemy

  element_num = nelemx * nelemy

  return
end
subroutine sphere_grid_q9_node_num ( nelemx, nelemy, node_num )

!*****************************************************************************80
!
!! SPHERE_GRID_Q9_NODE_NUM counts nodes in a Q9 sphere grid.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NELEMX, NELEMY, the number of elements along the
!    X and Y directions.
!
!    Output, integer ( kind = 4 ) NODE_NUM, the number of nodes in the grid.
!
  implicit none

  integer ( kind = 4 ) nelemx
  integer ( kind = 4 ) nelemy
  integer ( kind = 4 ) node_num

  node_num = 4 * nelemx * nelemy - 2 * nelemx + 2

  return
end
subroutine sphere_grid_q9_node_xyz ( nelemx, nelemy, node_xyz )

!*****************************************************************************80
!
!! SPHERE_GRID_Q9_NODE_XYZ produces node coordinates for a Q9 sphere grid.
!
!  Discussion:
!
!    The number of nodes to be generated is
!
!      NODE_NUM = 4 * NELEMX * NELEMY - 2 * NELEMX + 2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NELEMX, NELEMY, the number of elements along 
!    the X and Y directions.  
!
!    Output, real ( kind = 8 ) NODE_XYZ(3,NODE_NUM), 
!    the node coordinates.
!
  implicit none

  integer ( kind = 4 ) nelemx
  integer ( kind = 4 ) nelemy

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) node
  real    ( kind = 8 ) node_xyz(3,4*nelemx*nelemy-2*nelemx+2)
  real    ( kind = 8 ) phi
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real    ( kind = 8 ) theta

  node = 0

  node = node + 1      
  node_xyz(1,node) =  0.0D+00
  node_xyz(2,node) =  0.0D+00
  node_xyz(3,node) = -1.0D+00

  do j = 2 * nelemy, 2, -1

    phi = real ( j - 1, kind = 8 ) * pi &
        / real ( 2 * nelemy, kind = 8 )

    do i = 1, 2 * nelemx

      theta = real ( i - 1, kind = 8 ) * 2.0D+00 * pi &
            / real ( 2 * nelemx, kind = 8 )

      node = node + 1      
      node_xyz(1,node) = cos ( theta ) * sin ( phi )
      node_xyz(2,node) = sin ( theta ) * sin ( phi )
      node_xyz(3,node) =                 cos ( phi )

    end do
  end do

  node = node + 1      
  node_xyz(1,node) =  0.0D+00
  node_xyz(2,node) =  0.0D+00
  node_xyz(3,node) =  1.0D+00

  return
end
subroutine sphere_grid_q16_element ( nelemx, nelemy, element_node )

!*****************************************************************************80
!
!! SPHERE_GRID_Q16_ELEMENT produces a Q16 sphere grid.
!
!  Discussion:
!
!    This would be the same as the grid in a plane, except that all the
!    nodes along the bottom edge are identified (replaced by a single node
!    that is the south pole) and similarly for the top edge, and the
!    nodes on the extreme right edge are identified pairwise with those 
!    on the extreme left edge.
!
!  Example:
!
!    Input:
!
!      NELEMX = 2, NELEMY = 2
!
!    Output:
!
!      ELEMENT_NODE =
!         1,  1,  1,  1,  2,  3,  4,  5,  8,  9, 10, 11, 14, 15, 16, 17;
!         1,  1,  1,  1,  5,  6,  7,  2, 11, 12, 13,  8, 17, 18, 19, 14;
!        14, 15, 16, 17, 20, 21, 22, 23, 26, 27, 28, 29, 32, 32, 32, 32;
!        17, 18, 19, 14, 23, 24, 25, 20, 29, 30, 31, 26, 32, 32, 32, 32.
!
!  Grid:
!
!   32-32-32-32-32-32-32
!    |        |        |
!    |        |        |
!   26 27 28 29 30 31 26
!    |        |        |
!    |        |        |
!   20 21 22 23 24 25 20
!    |        |        |
!    | E3     | E4     |
!   14-15-16-17-18-19-14
!    |        |        |
!    |        |        |
!    8  9 10 11 12 13  8
!    |        |        |
!    |        |        |
!    2  3  4  5  6  7  2
!    |        |        |
!    | E1     | E2     |
!    1--1--1--1--1--1--1
!
!  Reference Element Q16:
!
!    |
!    1 13--14--15--16
!    |  |   :   :   |
!    |  |   :   :   |
!    |  9..10..11..12
!    S  |   :   :   |
!    |  |   :   :   |
!    |  5...6...7...8
!    |  |   :   :   |
!    |  |   :   :   |
!    0  1---2---3---4
!    |
!    +--0-----R-----1-->
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NELEMX, NELEMY, the number of elements along 
!    the X and Y directions.  The number of elements generated will be
!    NELEMX * NELEMY.
!
!    Output, integer ( kind = 4 ) ELEMENT_NODE(16,NELEMX*NELEMY), the nodes 
!    that form each element.
!
  implicit none

  integer ( kind = 4 ) nelemx
  integer ( kind = 4 ) nelemy
  integer ( kind = 4 ), parameter :: element_order = 16

  integer ( kind = 4 ) base1
  integer ( kind = 4 ) base2
  integer ( kind = 4 ) element
  integer ( kind = 4 ) element_node(element_order,nelemx*nelemy)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  element = 0

  do j = 1, nelemy

    base1 = ( j - 1 ) * 3 * ( 3 * nelemx ) + 2 - 3 * nelemx

    do i = 1, nelemx

      base2 = base1 + 3 * ( i - 1 )

      element = element + 1

      element_node( 1,element) = base2
      element_node( 2,element) = base2 + 1
      element_node( 3,element) = base2 + 2

      if ( i < nelemx ) then
        element_node( 4,element) = base2 + 3
      else
        element_node( 4,element) = base1
      end if

      element_node( 5,element) = element_node( 1,element) + 3 * nelemx
      element_node( 6,element) = element_node( 2,element) + 3 * nelemx
      element_node( 7,element) = element_node( 3,element) + 3 * nelemx
      element_node( 8,element) = element_node( 4,element) + 3 * nelemx

      element_node( 9,element) = element_node( 5,element) + 3 * nelemx
      element_node(10,element) = element_node( 6,element) + 3 * nelemx
      element_node(11,element) = element_node( 7,element) + 3 * nelemx
      element_node(12,element) = element_node( 8,element) + 3 * nelemx

      element_node(13,element) = element_node( 9,element) + 3 * nelemx
      element_node(14,element) = element_node(10,element) + 3 * nelemx
      element_node(15,element) = element_node(11,element) + 3 * nelemx
      element_node(16,element) = element_node(12,element) + 3 * nelemx

      if ( j == 1 ) then

        element_node( 1,element) = 1
        element_node( 2,element) = 1
        element_node( 3,element) = 1
        element_node( 4,element) = 1

      else if ( j == nelemy ) then

        element_node(13,element) = base1 + 9 * nelemx
        element_node(14,element) = base1 + 9 * nelemx
        element_node(15,element) = base1 + 9 * nelemx
        element_node(16,element) = base1 + 9 * nelemx

      end if

    end do
  end do

  return
end
subroutine sphere_grid_q16_element_num ( nelemx, nelemy, element_num )

!*****************************************************************************80
!
!! SPHERE_GRID_Q16_ELEMENT_NUM counts the elements in a Q16 sphere grid.
!
!  Example:
!
!    Input:
!
!      NELEMX = 3, NELEMY = 2
!
!    Output:
!
!      ELEMENT_NUM = NELEMX * NELEMY = 6
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NELEMX, NELEMY, the number of elements along 
!    the X and Y directions. 
!
!    Output, integer ( kind = 4 ) ELEMENT_NUM, the number of elements in 
!    the grid.
!
  implicit none

  integer ( kind = 4 ) element_num
  integer ( kind = 4 ) nelemx
  integer ( kind = 4 ) nelemy

  element_num = nelemx * nelemy

  return
end
subroutine sphere_grid_q16_node_num ( nelemx, nelemy, node_num )

!*****************************************************************************80
!
!! SPHERE_GRID_Q16_NODE_NUM counts nodes in a Q16 sphere grid.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NELEMX, NELEMY, the number of elements along the
!    X and Y directions.
!
!    Output, integer ( kind = 4 ) NODE_NUM, the number of nodes in the grid.
!
  implicit none

  integer ( kind = 4 ) nelemx
  integer ( kind = 4 ) nelemy
  integer ( kind = 4 ) node_num

  node_num = 9 * nelemx * nelemy - 3 * nelemx + 2

  return
end
subroutine sphere_grid_q16_node_xyz ( nelemx, nelemy, node_xyz )

!*****************************************************************************80
!
!! SPHERE_GRID_Q16_NODE_XYZ produces node coordinates for a Q16 sphere grid.
!
!  Discussion:
!
!    The number of nodes to be generated is
!
!      NODE_NUM = 9 * NELEMX * NELEMY - 3 * NELEMX + 2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NELEMX, NELEMY, the number of elements along 
!    the X and Y directions.
!
!    Output, real ( kind = 8 ) NODE_XYZ(3,NODE_NUM), the node coordinates.
!
  implicit none

  integer ( kind = 4 ) nelemx
  integer ( kind = 4 ) nelemy

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) node
  real    ( kind = 8 ) node_xyz(3,9*nelemx*nelemy-3*nelemx+2)
  real    ( kind = 8 ) phi
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real    ( kind = 8 ) theta

  node = 0

  do j = 3 * nelemy + 1, 1, -1

    phi = real ( j - 1, kind = 8 ) * pi &
        / real ( 3 * nelemy, kind = 8 )

    if ( j == 1 ) then

      node = node + 1
      node_xyz(1,node) =  0.0D+00
      node_xyz(2,node) =  0.0D+00
      node_xyz(3,node) =  1.0D+00

    else if ( j < 3 * nelemy + 1 ) then

      do i = 1, 3 * nelemx

        theta = real ( i - 1, kind = 8 ) * 2.0D+00 * pi &
              / real ( 3 * nelemx, kind = 8 )

        node = node + 1      
        node_xyz(1,node) = cos ( theta ) * sin ( phi )
        node_xyz(2,node) = sin ( theta ) * sin ( phi )
        node_xyz(3,node) =                 cos ( phi )

      end do

    else if ( j == 3 * nelemy + 1 ) then

      node = node + 1
      node_xyz(1,node) =  0.0D+00
      node_xyz(2,node) =  0.0D+00
      node_xyz(3,node) = -1.0D+00

    end if

  end do

  return
end
subroutine sphere_grid_t3_element ( nelemx, nelemy, element_node )

!*****************************************************************************80
!
!! SPHERE_GRID_T3_ELEMENT produces a T3 sphere grid.
!
!  Discussion:
!
!    This would be the same as the grid in a plane, except that all the
!    nodes along the bottom edge are identified (replaced by a single node
!    that is the south pole) and similarly for the top edge, and the
!    nodes on the extreme right edge are identified pairwise with those 
!    on the extreme left edge.
!
!  Example:
!
!    Input:
!
!      NELEMX = 3, NELEMY = 4
!
!    Output:
!
!      ELEMENT_NODE =
!         1,  3,  2;
!         1,  4,  3;
!         1,  2,  4;
!         2,  3,  5
!         6,  5,  3
!         3,  4,  6
!         7,  6,  4;
!         4,  2,  7;
!         5,  7,  2;
!         5,  6,  8;
!         9,  8,  6;
!         6,  7,  9;
!        10,  9,  7;
!         7,  5, 10;
!         8, 10,  5;
!         8,  9, 11;
!         9, 10, 11;
!        10,  8, 11;
!
!  Grid:
!
!   11    11    11    11
!    | \   | \   | \   |
!    |  \  |  \  |  \  |
!    |E16\ |E17 \|E18\ |
!    8-----9----10-----8
!    | \E11| \E13| \E15|
!    |  \  |  \  |  \  |
!    |E10\ |E12\ |E14\ |
!    5-----6-----7-----5
!    | \E5 | \E7 | \E9 |
!    |  \  |  \  |  \  |
!    |E4 \ |E6 \ |E8 \ |
!    2-----3-----4-----2
!      \E1 | \E2 | \E3 |
!       \  |  \  |  \  |
!        \ |   \ |   \ |
!          1     1     1
!
!  Reference Element T3:
!
!    |
!    1  3
!    |  |\
!    |  | \
!    S  |  \
!    |  |   \
!    |  |    \
!    0  1-----2
!    |
!    +--0--R--1-->
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NELEMX, NELEMY, the number of elements along 
!    the X and Y directions.
!
!    Output, integer ( kind = 4 ) ELEMENT_NODE(3,2*NELEMX*(NELEMY-1)), the nodes
!    that form each element.
!
  implicit none

  integer ( kind = 4 ) nelemx
  integer ( kind = 4 ) nelemy
  integer ( kind = 4 ), parameter :: element_order = 3

  integer ( kind = 4 ) base1
  integer ( kind = 4 ) base2
  integer ( kind = 4 ) element
  integer ( kind = 4 ) element_node(element_order,nelemx*nelemy)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ne
  integer ( kind = 4 ) nw
  integer ( kind = 4 ) se
  integer ( kind = 4 ) sw

  element = 0

  do j = 1, nelemy

    base1 = ( j - 1 ) * nelemx + 2 - nelemx

    do i = 1, nelemx

      base2 = base1 + i - 1

      sw = base2
      if ( i < nelemx ) then
        se = base2 + 1
      else
        se = base1
      end if
      nw = sw + nelemx
      ne = se + nelemx

      if ( j == 1 ) then
        sw = 1
        se = 1
      else if ( j == nelemx ) then
        nw = base1 + nelemx
        ne = base1 + nelemx
      end if

      if ( 1 < j ) then
        element = element + 1
        element_node(1,element) = sw
        element_node(2,element) = se
        element_node(3,element) = nw
      end if

      if ( j < nelemy ) then
        element = element + 1
        element_node(1,element) = ne
        element_node(2,element) = nw
        element_node(3,element) = se
      end if

    end do
  end do

  return
end
subroutine sphere_grid_t3_element_num ( nelemx, nelemy, element_num )

!*****************************************************************************80
!
!! SPHERE_GRID_T3_ELEMENT_NUM counts the elements in a T3 sphere grid.
!
!  Example:
!
!    Input:
!
!      NELEMX = 6, NELEMY = 4
!
!    Output:
!
!      ELEMENT_NUM = 2 * NELEMX * ( NELEMY - 1 ) = 36
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NELEMX, NELEMY, the number of elements along
!    the X and Y directions. 
!
!    Output, integer ( kind = 4 ) ELEMENT_NUM, the number of elements in 
!    the grid.
!
  implicit none

  integer ( kind = 4 ) element_num
  integer ( kind = 4 ) nelemx
  integer ( kind = 4 ) nelemy

  element_num = 2 * nelemx * ( nelemy - 1 )

  return
end
subroutine sphere_grid_t3_node_num ( nelemx, nelemy, node_num )

!*****************************************************************************80
!
!! SPHERE_GRID_T3_NODE_NUM counts nodes in a T3 sphere grid.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NELEMX, NELEMY, the number of elements along 
!    the X and Y directions.
!
!    Output, integer ( kind = 4 ) NODE_NUM, the number of nodes in the grid.
!
  implicit none

  integer ( kind = 4 ) nelemx
  integer ( kind = 4 ) nelemy
  integer ( kind = 4 ) node_num

  node_num = nelemx * ( nelemy - 1 ) + 2

  return
end
subroutine sphere_grid_t3_node_xyz ( nelemx, nelemy, node_xyz )

!*****************************************************************************80
!
!! SPHERE_GRID_T3_NODE_XYZ produces node coordinates for a T3 sphere grid.
!
!  Discussion:
!
!    The number of nodes to be generated is
!
!      NODE_NUM = NELEMX * ( NELEMY - 1 ) + 2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NELEMX, NELEMY, the number of elements along 
!    the X and Y directions. 
!
!    Output, real ( kind = 8 ) NODE_XYZ(3,NODE_NUM), 
!    the node coordinates.
!
  implicit none

  integer ( kind = 4 ) nelemx
  integer ( kind = 4 ) nelemy

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) node
  real    ( kind = 8 ) node_xyz(3,nelemx*(nelemy-1)+2)
  real    ( kind = 8 ) phi
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real    ( kind = 8 ) theta

  node = 0

  node = node + 1      
  node_xyz(1,node) =  0.0D+00
  node_xyz(2,node) =  0.0D+00
  node_xyz(3,node) = -1.0D+00

  do j = nelemy, 2, -1

    phi = real ( j - 1, kind = 8 ) * pi &
        / real ( nelemy, kind = 8 )

    do i = 1, nelemx

      theta = real ( i - 1, kind = 8 ) * 2.0D+00 * pi &
            / real ( nelemx, kind = 8 )

      node = node + 1      
      node_xyz(1,node) = cos ( theta ) * sin ( phi )
      node_xyz(2,node) = sin ( theta ) * sin ( phi )
      node_xyz(3,node) =                 cos ( phi )

    end do
  end do

  node = node + 1      
  node_xyz(1,node) =  0.0D+00
  node_xyz(2,node) =  0.0D+00
  node_xyz(3,node) =  1.0D+00

  return
end
subroutine sphere_grid_t6_element ( nelemx, nelemy, element_node )

!*****************************************************************************80
!
!! SPHERE_GRID_T6_ELEMENT produces a T6 sphere grid.
!
!  Discussion:
!
!    This would be the same as the grid in a plane, except that all the
!    nodes along the bottom edge are identified (replaced by a single node
!    that is the south pole) and similarly for the top edge, and the
!    nodes on the extreme right edge are identified pairwise with those 
!    on the extreme left edge.
!
!  Example:
!
!    Input:
!
!      NELEMX = 3, NELEMY = 4
!
!    Output:
!
!      ELEMENT_NODE =
!        10,  8,  1,  9,  3,  4;
!        12, 10,  1, 11,  5,  6;
!         8, 12,  1, 13,  7,  2;
!         8, 10, 20,  9, 15, 14;
!        22, 20, 10, 21, 15, 16;
!        10, 12, 22, 11, 17, 16;
!        24, 22, 12, 23, 17, 18;
!        12,  8, 24, 13, 19, 18;
!        20, 24,  8, 25, 19, 14;
!        20, 22, 32, 21, 27, 26;
!        34, 32, 22, 33, 27, 28;
!        22, 24, 34, 23, 29, 28;
!        36, 34, 24, 35, 29, 30;
!        24, 20, 36, 25, 31, 30;
!        32, 36, 20, 37, 31, 26;
!        32, 34, 44, 33, 39, 38;
!        34, 36, 44, 35, 41, 40;
!        36, 32, 44, 37, 43, 42;
!
!  Grid:
!
!   44    44    44
!    |\    |\    |\
!   38 39 40 41 42 43 
!    |    \|    \|    \
!   32-33-34-35-36-37-32
!    |\    |\    |\    |
!   26 27 28 29 30 31 26
!    |    \|    \|    \|
!   20-21-22-23-24-25-20
!    |\    |\    |\    |
!   14 15 16 17 18 19 14
!    |    \|    \|    \|
!    8--9-10-11-12-13--8
!     \    |\    |\    |
!       3  4  5  6  7  2
!         \|    \|    \|
!          1     1     1
!
!  Reference Element T6:
!
!    |
!    1  3
!    |  |\
!    |  | \
!    S  6  5
!    |  |   \
!    |  |    \
!    0  1--4--2
!    |
!    +--0--R--1-->
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NELEMX, NELEMY, the number of elements along the
!    X and Y directions.
!
!    Output, integer ( kind = 4 ) ELEMENT_NODE(6,2*NELEMX*(NELEMY-1)), 
!    the nodes that form each element.
!
  implicit none

  integer ( kind = 4 ) nelemx
  integer ( kind = 4 ) nelemy
  integer ( kind = 4 ), parameter :: element_order = 6

  integer ( kind = 4 ) base1
  integer ( kind = 4 ) base2
  integer ( kind = 4 ) c
  integer ( kind = 4 ) e
  integer ( kind = 4 ) element
  integer ( kind = 4 ) element_node(element_order,4*nelemx*nelemy-2*nelemx+2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ne
  integer ( kind = 4 ) nw
  integer ( kind = 4 ) s
  integer ( kind = 4 ) se
  integer ( kind = 4 ) sw
  integer ( kind = 4 ) w

  element = 0

  do j = 1, nelemy

    base1 = ( j - 1 ) * 2 * ( 2 * nelemx ) + 2 - 2 * nelemx

    do i = 1, nelemx

      base2 = base1 + 2 * ( i - 1 )

      sw = base2
      s = base2 + 1
      if ( i < nelemx ) then
        se = base2 + 2
      else
        se = base1
      end if

      w = sw + 2 * nelemx
      c = s  + 2 * nelemx
      e = se + 2 * nelemx

      nw = w + 2 * nelemx
      n  = c + 2 * nelemx
      ne = e + 2 * nelemx

      if ( j == 1 ) then
        sw = 1
        s  = 1
        se = 1
      else if ( j == nelemy ) then
        nw = base1 + 4 * nelemx
        n  = base1 + 4 * nelemx
        ne = base1 + 4 * nelemx
      end if

      if ( 1 < j ) then
        element = element + 1
        element_node(1,element) = sw
        element_node(2,element) = se
        element_node(3,element) = nw
        element_node(4,element) = s
        element_node(5,element) = c
        element_node(6,element) = w
      end if

      if ( j < nelemy ) then
        element = element + 1
        element_node(1,element) = ne
        element_node(2,element) = nw
        element_node(3,element) = se
        element_node(4,element) = n
        element_node(5,element) = c
        element_node(6,element) = e
      end if


    end do

  end do

  return
end
subroutine sphere_grid_t6_element_num ( nelemx, nelemy, element_num )

!*****************************************************************************80
!
!! SPHERE_GRID_T6_ELEMENT_NUM counts the elements in a T6 sphere grid.
!
!  Example:
!
!    Input:
!
!      NELEMX = 6, NELEMY = 4
!
!    Output:
!
!      ELEMENT_NUM = 2 * NELEMX * ( NELEMY - 1 ) = 36
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NELEMX, NELEMY, the number of elements along the
!    X and Y directions. 
!
!    Output, integer ( kind = 4 ) ELEMENT_NUM, the number of elements in the grid.
!
  implicit none

  integer ( kind = 4 ) element_num
  integer ( kind = 4 ) nelemx
  integer ( kind = 4 ) nelemy

  element_num = 2 * nelemx * ( nelemy - 1 )

  return
end
subroutine sphere_grid_t6_node_num ( nelemx, nelemy, node_num )

!*****************************************************************************80
!
!! SPHERE_GRID_T6_NODE_NUM counts nodes in a T6 sphere grid.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NELEMX, NELEMY, the number of elements along the
!    X and Y directions.
!
!    Output, integer ( kind = 4 ) NODE_NUM, the number of nodes in the grid.
!
  implicit none

  integer ( kind = 4 ) nelemx
  integer ( kind = 4 ) nelemy
  integer ( kind = 4 ) node_num

  node_num = 4 * nelemx * nelemy - 2 * nelemx + 2

  return
end
subroutine sphere_grid_t6_node_xyz ( nelemx, nelemy, node_xyz )

!*****************************************************************************80
!
!! SPHERE_GRID_T6_NODE_XYZ produces node coordinates for a T6 sphere grid.
!
!  Discussion:
!
!    The number of nodes to be generated is
!
!      NODE_NUM = 4 * NELEMX * NELEMY - 2 * NELEMX + 2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NELEMX, NELEMY, the number of elements along the
!    X and Y directions.  
!
!    Output, real ( kind = 8 ) NODE_XYZ(3,NODE_NUM), 
!    the node coordinates.
!
  implicit none

  integer ( kind = 4 ) nelemx
  integer ( kind = 4 ) nelemy

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) node
  real    ( kind = 8 ) node_xyz(3,4*nelemx*nelemy-2*nelemx+2)
  real    ( kind = 8 ) phi
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real    ( kind = 8 ) theta

  node = 0

  node = node + 1      
  node_xyz(1,node) =  0.0D+00
  node_xyz(2,node) =  0.0D+00
  node_xyz(3,node) = -1.0D+00

  do j = 2 * nelemy, 2, -1

    phi = real ( j - 1, kind = 8 ) * pi &
        / real ( 2 * nelemy, kind = 8 )

    do i = 1, 2 * nelemx

      theta = real ( i - 1, kind = 8 ) * 2.0D+00 * pi &
            / real ( 2 * nelemx, kind = 8 )

      node = node + 1      
      node_xyz(1,node) = cos ( theta ) * sin ( phi )
      node_xyz(2,node) = sin ( theta ) * sin ( phi )
      node_xyz(3,node) =                 cos ( phi )

    end do
  end do

  node = node + 1      
  node_xyz(1,node) =  0.0D+00
  node_xyz(2,node) =  0.0D+00
  node_xyz(3,node) =  1.0D+00

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
subroutine triangle_unit_set ( rule, xtab, ytab, weight )

!*****************************************************************************80
!
!! TRIANGLE_UNIT_SET sets a quadrature rule in the unit triangle.
!
!  Discussion:
!
!    The user is responsible for determining the value of ORDER,
!    and appropriately dimensioning the arrays XTAB, YTAB and
!    WEIGHT so that they can accommodate the data.
!
!    The value of ORDER for each rule can be found by invoking
!    the function TRIANGLE_RULE_SIZE.
!
!    The integration region is:
!
!      0 <= X and 0 <= Y and X + Y <= 1.
!
!  Graph:
!
!      ^
!    1 | *
!      | |\
!    Y | | \
!      | |  \
!    0 | *---*
!      +------->
!        0 X 1
!
!   The rules are accessed by an index number, RULE.  The indices,
!   and the descriptions of the corresponding rules, are:
!
!     1, ORDER =  1, precision 1, Zienkiewicz #1.
!     2, ORDER =  2, precision 1, (the "vertex rule").
!     3, ORDER =  3, precision 2, Strang and Fix formula #1.
!     4, ORDER =  3, precision 2, Strang and Fix formula #2,
!                                 Zienkiewicz #2.
!     5, ORDER =  4, precision 3, Strang and Fix formula #3,
!                                 Zienkiewicz #3.
!     6, ORDER =  6, precision 3, Strang and Fix formula #4.
!     7, ORDER =  6, precision 3, Stroud formula T2:3-1.
!     8, ORDER =  6, precision 4, Strang and Fix formula #5.
!     9, ORDER =  7, precision 4, Strang and Fix formula #6.
!    10, ORDER =  7, precision 5, Strang and Fix formula #7,
!                                 Stroud formula T2:5-1,
!                                 Zienkiewicz #4,
!                                 Schwarz Table 2.2.
!    11, ORDER =  9, precision 6, Strang and Fix formula #8.
!    12, ORDER = 12, precision 6, Strang and Fix formula #9.
!    13, ORDER = 13, precision 7, Strang and Fix formula #10.
!        Note that there is a typographical error in Strang and Fix
!        which lists the value of the XSI(3) component of the
!        last generator point as 0.4869... when it should be 0.04869...
!    14, ORDER =  7, precision ?.
!    15, ORDER = 16, precision 7, conical product Gauss, Stroud formula T2:7-1.
!    16, ORDER = 64, precision 15, triangular product Gauss rule.
!    17, ORDER = 19, precision 8, from CUBTRI, ACM TOMS #584.
!    18, ORDER = 19, precision 9, from TRIEX, ACM TOMS #612.
!    19, ORDER = 28, precision 11, from TRIEX, ACM TOMS #612.
!    20, ORDER = 37, precision 13, from ACM TOMS #706.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jarle Berntsen, Terje Espelid,
!    Algorithm 706,
!    DCUTRI: an algorithm for adaptive cubature over a collection of triangles,
!    ACM Transactions on Mathematical Software,
!    Volume 18, Number 3, September 1992, pages 329-342.
!
!    Elise deDoncker, Ian Robinson,
!    Algorithm 612:
!    Integration over a Triangle Using Nonlinear Extrapolation,
!    ACM Transactions on Mathematical Software,
!    Volume 10, Number 1, March 1984, pages 17-22.
!
!    Dirk Laurie,
!    Algorithm 584,
!    CUBTRI, Automatic Cubature Over a Triangle,
!    ACM Transactions on Mathematical Software,
!    Volume 8, Number 2, 1982, pages 210-218.
!
!    James Lyness, Dennis Jespersen,
!    Moderate Degree Symmetric Quadrature Rules for the Triangle,
!    Journal of the Institute of Mathematics and its Applications,
!    Volume 15, Number 1, February 1975, pages 19-32.
!
!    Hans Rudolf Schwarz,
!    Finite Element Methods,
!    Academic Press, 1988,
!    ISBN: 0126330107,
!    LC: TA347.F5.S3313.
!
!    Gilbert Strang, George Fix,
!    An Analysis of the Finite Element Method,
!    Cambridge, 1973,
!    ISBN: 096140888X,
!    LC: TA335.S77.
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!    Olgierd Zienkiewicz,
!    The Finite Element Method,
!    Sixth Edition,
!    Butterworth-Heinemann, 2005,
!    ISBN: 0750663200,
!    LC: TA640.2.Z54
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) RULE, the index of the rule.
!
!    Output, real ( kind = 8 ) XTAB(ORDER), YTAB(ORDER), the abscissas.
!
!    Output, real ( kind = 8 ) WEIGHT(ORDER), the weights of the rule.
!
  implicit none

  real    ( kind = 8 ) a
  real    ( kind = 8 ) b
  real    ( kind = 8 ) c
  real    ( kind = 8 ) d
  real    ( kind = 8 ) e
  real    ( kind = 8 ) f
  real    ( kind = 8 ) g
  real    ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) order2
  real    ( kind = 8 ) p
  real    ( kind = 8 ) q
  real    ( kind = 8 ) r
  integer ( kind = 4 ) rule
  real    ( kind = 8 ) s
  real    ( kind = 8 ) t
  real    ( kind = 8 ) u
  real    ( kind = 8 ) v
  real    ( kind = 8 ) w
  real    ( kind = 8 ) w1
  real    ( kind = 8 ) w2
  real    ( kind = 8 ) w3
  real    ( kind = 8 ) w4
  real    ( kind = 8 ) w5
  real    ( kind = 8 ) w6
  real    ( kind = 8 ) w7
  real    ( kind = 8 ) w8
  real    ( kind = 8 ) w9
  real    ( kind = 8 ) weight(*)
  real    ( kind = 8 ) weight1(8)
  real    ( kind = 8 ) weight2(8)
  real    ( kind = 8 ) wx
  real    ( kind = 8 ) x
  real    ( kind = 8 ) xtab(*)
  real    ( kind = 8 ) xtab1(8)
  real    ( kind = 8 ) xtab2(8)
  real    ( kind = 8 ) y
  real    ( kind = 8 ) ytab(*)
  real    ( kind = 8 ) z
!
!  1 point, precision 1.
!
  if ( rule == 1 ) then

    xtab(1)   = 0.33333333333333333333D+00

    ytab(1)   = 0.33333333333333333333D+00

    weight(1) = 1.00000000000000000000D+00
!
!  3 points, precision 1, the "vertex rule".
!
  else if ( rule == 2 ) then

    xtab(1) =   1.00000000000000000000D+00
    xtab(2) =   0.00000000000000000000D+00
    xtab(3) =   0.00000000000000000000D+00

    ytab(1) =   0.00000000000000000000D+00
    ytab(2) =   1.00000000000000000000D+00
    ytab(3) =   0.00000000000000000000D+00

    weight(1) = 0.33333333333333333333D+00
    weight(2) = 0.33333333333333333333D+00
    weight(3) = 0.33333333333333333333D+00
!
!  3 points, precision 2, Strang and Fix formula #1.
!
  else if ( rule == 3 ) then

    xtab(1)   = 0.66666666666666666667D+00
    xtab(2)   = 0.16666666666666666667D+00
    xtab(3)   = 0.16666666666666666667D+00

    ytab(1)   = 0.16666666666666666667D+00
    ytab(2)   = 0.66666666666666666667D+00
    ytab(3)   = 0.16666666666666666667D+00

    weight(1) = 0.33333333333333333333D+00
    weight(2) = 0.33333333333333333333D+00
    weight(3) = 0.33333333333333333333D+00
!
!  3 points, precision 2, Strang and Fix formula #2.
!
  else if ( rule == 4 ) then

    xtab(1)   = 0.50000000000000000000D+00
    xtab(2)   = 0.50000000000000000000D+00
    xtab(3)   = 0.00000000000000000000D+00

    ytab(1)   = 0.00000000000000000000D+00
    ytab(2)   = 0.50000000000000000000D+00
    ytab(3)   = 0.50000000000000000000D+00

    weight(1) = 0.33333333333333333333D+00
    weight(2) = 0.33333333333333333333D+00
    weight(3) = 0.33333333333333333333D+00
!
!  4 points, precision 3, Strang and Fix formula #3.
!
  else if ( rule == 5 ) then

    a =   6.0D+00
    b =  10.0D+00
    c =  18.0D+00
    d =  25.0D+00
    e = -27.0D+00
    f =  30.0D+00
    g =  48.0D+00

    xtab(1:4) =   (/ b, c, a, a /) / f
    ytab(1:4) =   (/ b, a, c, a /) / f
    weight(1:4) = (/ e, d, d, d /) / g
!
!  6 points, precision 3, Strang and Fix formula #4.
!
  else if ( rule == 6 ) then

    a = 0.659027622374092D+00
    b = 0.231933368553031D+00
    c = 0.109039009072877D+00

    xtab(1:6) =   (/ a, a, b, b, c, c /)
    ytab(1:6) =   (/ b, c, a, c, a, b /)

    weight(1) = 0.16666666666666666667D+00
    weight(2) = 0.16666666666666666667D+00
    weight(3) = 0.16666666666666666667D+00
    weight(4) = 0.16666666666666666667D+00
    weight(5) = 0.16666666666666666667D+00
    weight(6) = 0.16666666666666666667D+00
!
!  6 points, precision 3, Stroud T2:3-1.
!
  else if ( rule == 7 ) then

    a = 0.0D+00
    b = 0.5D+00
    c = 2.0D+00 /  3.0D+00
    d = 1.0D+00 /  6.0D+00
    v = 1.0D+00 / 30.0D+00
    w = 3.0D+00 / 10.0D+00

    xtab(1:6) =   (/ a, b, b, c, d, d /)
    ytab(1:6) =   (/ b, a, b, d, c, d /)
    weight(1:6) = (/ v, v, v, w, w, w /)
!
!  6 points, precision 4, Strang and Fix, formula #5.
!
  else if ( rule == 8 ) then

    a = 0.816847572980459D+00
    b = 0.091576213509771D+00
    c = 0.108103018168070D+00
    d = 0.445948490915965D+00
    v = 0.109951743655322D+00
    w = 0.223381589678011D+00

    xtab(1:6) =   (/ a, b, b, c, d, d /)
    ytab(1:6) =   (/ b, a, b, d, c, d /)
    weight(1:6) = (/ v, v, v, w, w, w /)
!
!  7 points, precision 4, Strang and Fix formula #6.
!
  else if ( rule == 9 ) then

    a = 1.0D+00 / 3.0D+00
    c = 0.736712498968435D+00
    d = 0.237932366472434D+00
    e = 0.025355134551932D+00
    v = 0.375000000000000D+00
    w = 0.104166666666667D+00

    xtab(1:7) =   (/ a, c, c, d, d, e, e /)
    ytab(1:7) =   (/ a, d, e, c, e, c, d /)
    weight(1:7) = (/ v, w, w, w, w, w, w /)
!
!  7 points, precision 5, Strang and Fix formula #7, Stroud T2:5-1
!
  else if ( rule == 10 ) then

    a = 1.0D+00 / 3.0D+00
    b = ( 9.0D+00 + 2.0D+00 * sqrt ( 15.0D+00 ) ) / 21.0D+00
    c = ( 6.0D+00 -           sqrt ( 15.0D+00 ) ) / 21.0D+00
    d = ( 9.0D+00 - 2.0D+00 * sqrt ( 15.0D+00 ) ) / 21.0D+00
    e = ( 6.0D+00 +           sqrt ( 15.0D+00 ) ) / 21.0D+00
    u = 0.225D+00
    v = ( 155.0D+00 - sqrt ( 15.0D+00 ) ) / 1200.0D+00
    w = ( 155.0D+00 + sqrt ( 15.0D+00 ) ) / 1200.0D+00

    xtab(1:7) =   (/ a, b, c, c, d, e, e /)
    ytab(1:7) =   (/ a, c, b, c, e, d, e /)
    weight(1:7) = (/ u, v, v, v, w, w, w /)
!
!  9 points, precision 6, Strang and Fix formula #8.
!
  else if ( rule == 11 ) then

    a = 0.124949503233232D+00
    b = 0.437525248383384D+00
    c = 0.797112651860071D+00
    d = 0.165409927389841D+00
    e = 0.037477420750088D+00

    u = 0.205950504760887D+00
    v = 0.063691414286223D+00

    xtab(1:9) =   (/ a, b, b, c, c, d, d, e, e /)
    ytab(1:9) =   (/ b, a, b, d, e, c, e, c, d /)
    weight(1:9) = (/ u, u, u, v, v, v, v, v, v /)
!
!  12 points, precision 6, Strang and Fix, formula #9.
!
  else if ( rule == 12 ) then

    a = 0.873821971016996D+00
    b = 0.063089014491502D+00
    c = 0.501426509658179D+00
    d = 0.249286745170910D+00
    e = 0.636502499121399D+00
    f = 0.310352451033785D+00
    g = 0.053145049844816D+00

    u = 0.050844906370207D+00
    v = 0.116786275726379D+00
    w = 0.082851075618374D+00

    xtab(1:12) =   (/ a, b, b, c, d, d, e, e, f, f, g, g /)
    ytab(1:12) =   (/ b, a, b, d, c, d, f, g, e, g, e, f /)
    weight(1:12) = (/ u, u, u, v, v, v, w, w, w, w, w, w /)
!
!  13 points, precision 7, Strang and Fix, formula #10.
!
!  Note that there is a typographical error in Strang and Fix
!  which lists the value of the XSI(3) component of the
!  last generator point as 0.4869... when it should be 0.04869...
!
  else if ( rule == 13 ) then

    h = 1.0D+00 / 3.0D+00
    a = 0.479308067841923D+00
    b = 0.260345966079038D+00
    c = 0.869739794195568D+00
    d = 0.065130102902216D+00
    e = 0.638444188569809D+00
    f = 0.312865496004875D+00
    g = 0.048690315425316D+00

    w = -0.149570044467670D+00
    t =  0.175615257433204D+00
    u =  0.053347235608839D+00
    v =  0.077113760890257D+00

    xtab(1:13) =   (/ h, a, b, b, c, d, d, e, e, f, f, g, g /)
    ytab(1:13) =   (/ h, b, a, b, d, c, d, f, g, e, g, e, f /)
    weight(1:13) = (/ w, t, t, t, u, u, u, v, v, v, v, v, v /)
!
!  7 points, precision ?.
!
  else if ( rule == 14 ) then

    a = 1.0D+00 / 3.0D+00
    b = 1.0D+00
    c = 0.5D+00
    z = 0.0D+00

    u = 27.0D+00 / 60.0D+00
    v =  3.0D+00 / 60.0D+00
    w =  8.0D+00 / 60.0D+00

    xtab(1:7) =   (/ a, b, z, z, z, c, c /)
    ytab(1:7) =   (/ a, z, b, z, c, z, c /)
    weight(1:7) = (/ u, v, v, v, w, w, w /)
!
!  16 points, Stroud T2:7-1.
!
  else if ( rule == 15 ) then
!
!  Legendre rule of order 4.
!
    order2 = 4

    xtab1(1:4) = (/ &
      -0.861136311594052575223946488893D+00, &
      -0.339981043584856264802665759103D+00, &
       0.339981043584856264802665759103D+00, &
       0.861136311594052575223946488893D+00 /)

    weight1(1:4) = (/ &
      0.347854845137453857373063949222D+00, &
      0.652145154862546142626936050778D+00, &
      0.652145154862546142626936050778D+00, &
      0.347854845137453857373063949222D+00 /)

    xtab1(1:order2) = 0.5D+00 * ( xtab1(1:order2) + 1.0D+00 )

    weight2(1) = 0.1355069134D+00
    weight2(2) = 0.2034645680D+00
    weight2(3) = 0.1298475476D+00
    weight2(4) = 0.0311809709D+00

    xtab2(1) = 0.0571041961D+00
    xtab2(2) = 0.2768430136D+00
    xtab2(3) = 0.5835904324D+00
    xtab2(4) = 0.8602401357D+00

    k = 0
    do i = 1, order2
      do j = 1, order2
        k = k + 1
        xtab(k) = xtab2(j)
        ytab(k) = xtab1(i) * ( 1.0D+00 - xtab2(j) )
        weight(k) = weight1(i) * weight2(j)
      end do
    end do
!
!  64 points, precision 15.
!
  else if ( rule == 16 ) then
!
!  Legendre rule of order 8.
!
    order2 = 8

    xtab1(1) = -0.960289856497536231683560868569D+00
    xtab1(2) = -0.796666477413626739591553936476D+00
    xtab1(3) = -0.525532409916328985817739049189D+00
    xtab1(4) = -0.183434642495649804939476142360D+00
    xtab1(5) =  0.183434642495649804939476142360D+00
    xtab1(6) =  0.525532409916328985817739049189D+00
    xtab1(7) =  0.796666477413626739591553936476D+00
    xtab1(8) =  0.960289856497536231683560868569D+00

    weight1(1) = 0.101228536290376259152531354310D+00
    weight1(2) = 0.222381034453374470544355994426D+00
    weight1(3) = 0.313706645877887287337962201987D+00
    weight1(4) = 0.362683783378361982965150449277D+00
    weight1(5) = 0.362683783378361982965150449277D+00
    weight1(6) = 0.313706645877887287337962201987D+00
    weight1(7) = 0.222381034453374470544355994426D+00
    weight1(8) = 0.101228536290376259152531354310D+00

    weight2(1) = 0.00329519144D+00
    weight2(2) = 0.01784290266D+00
    weight2(3) = 0.04543931950D+00
    weight2(4) = 0.07919959949D+00
    weight2(5) = 0.10604735944D+00
    weight2(6) = 0.11250579947D+00
    weight2(7) = 0.09111902364D+00
    weight2(8) = 0.04455080436D+00

    xtab2(1) = 0.04463395529D+00
    xtab2(2) = 0.14436625704D+00
    xtab2(3) = 0.28682475714D+00
    xtab2(4) = 0.45481331520D+00
    xtab2(5) = 0.62806783542D+00
    xtab2(6) = 0.78569152060D+00
    xtab2(7) = 0.90867639210D+00
    xtab2(8) = 0.98222008485D+00

    k = 0
    do j = 1, order2
      do i = 1, order2
        k = k + 1
        xtab(k) = 1.0D+00 - xtab2(j)
        ytab(k) = 0.5D+00 * ( 1.0D+00 + xtab1(i) ) * xtab2(j)
        weight(k) = weight1(i) * weight2(j)
      end do
    end do
!
!  19 points, precision 8, from CUBTRI.
!
  else if ( rule == 17 ) then

    a = 1.0D+00 / 3.0D+00
    b = ( 9.0D+00 + 2.0D+00 * sqrt ( 15.0D+00 ) ) / 21.0D+00
    c = ( 6.0D+00 -       sqrt ( 15.0D+00 ) ) / 21.0D+00
    d = ( 9.0D+00 - 2.0D+00 * sqrt ( 15.0D+00 ) ) / 21.0D+00
    e = ( 6.0D+00 +       sqrt ( 15.0D+00 ) ) / 21.0D+00
    f = ( 40.0D+00 - 10.0D+00 * sqrt ( 15.0D+00 ) &
      + 10.0D+00 * sqrt ( 7.0D+00 ) + 2.0D+00 * sqrt ( 105.0D+00 ) ) / 90.0D+00
    g = ( 25.0D+00 +  5.0D+00 * sqrt ( 15.0D+00 ) &
      -  5.0D+00 * sqrt ( 7.0D+00 ) - sqrt ( 105.0D+00 ) ) / 90.0D+00
    p = ( 40.0D+00 + 10.0D+00 * sqrt ( 15.0D+00 ) &
      + 10.0D+00 * sqrt ( 7.0D+00 ) - 2.0D+00 * sqrt ( 105.0D+00 ) ) / 90.0D+00
    q = ( 25.0D+00 -  5.0D+00 * sqrt ( 15.0D+00 ) &
      -  5.0D+00 * sqrt ( 7.0D+00 ) + sqrt ( 105.0D+00 ) ) / 90.0D+00
    r = ( 40.0D+00 + 10.0D+00 * sqrt ( 7.0D+00 ) ) / 90.0D+00
    s = ( 25.0D+00 +  5.0D+00 * sqrt ( 15.0D+00 ) - 5.0D+00 * sqrt ( 7.0D+00 ) &
      - sqrt ( 105.0D+00 ) ) / 90.0D+00
    t = ( 25.0D+00 -  5.0D+00 * sqrt ( 15.0D+00 ) - 5.0D+00 * sqrt ( 7.0D+00 ) &
      + sqrt ( 105.0D+00 ) ) / 90.0D+00

    w1 = ( 7137.0D+00 - 1800.0D+00 * sqrt ( 7.0D+00 ) ) / 62720.0D+00
    w2 = -9301697.0D+00 / 4695040.0D+00 - 13517313.0D+00 * sqrt ( 15.0D+00 ) &
      / 23475200.0D+00 + 764885.0D+00 * sqrt ( 7.0D+00 ) / 939008.0D+00 &
      + 198763.0D+00 * sqrt ( 105.0D+00 ) / 939008.0D+00
    w2 = w2 / 3.0D+00
    w3 = -9301697.0D+00 / 4695040.0D+00 + 13517313.0D+00 * sqrt ( 15.0D+00 ) &
      / 23475200.0D+00 &
      + 764885.0D+00 * sqrt ( 7.0D+00 ) / 939008.0D+00 &
      - 198763.0D+00 * sqrt ( 105.0D+00 ) / 939008.0D+00
    w3 = w3 / 3.0D+00
    w4 = ( 102791225.0D+00 - 23876225.0D+00 * sqrt ( 15.0D+00 ) &
      - 34500875.0D+00 * sqrt ( 7.0D+00 ) &
      + 9914825.0D+00 * sqrt ( 105.0D+00 ) ) / 59157504.0D+00
    w4 = w4 / 3.0D+00
    w5 = ( 102791225.0D+00 + 23876225.0D+00 * sqrt ( 15.0D+00 ) &
      - 34500875.0D+00 * sqrt ( 7.0D+00 ) &
      - 9914825D+00 * sqrt ( 105.0D+00 ) ) / 59157504.0D+00
    w5 = w5 / 3.0D+00
    w6 = ( 11075.0D+00 - 3500.0D+00 * sqrt ( 7.0D+00 ) ) / 8064.0D+00
    w6 = w6 / 6.0D+00

    xtab(1:19) =   (/  a,  b,  c,  c,  d,  e,  e,  f,  g,  g,  p,  q,  q, &
                       r,  r,  s,  s,  t,  t /)
    ytab(1:19) =   (/  a,  c,  b,  c,  e,  d,  e,  g,  f,  g,  q,  p,  q, &
                       s,  t,  r,  t,  r,  s /)
    weight(1:19) = (/ w1, w2, w2, w2, w3, w3, w3, w4, w4, w4, w5, w5, w5, &
      w6, w6, w6, w6, w6, w6 /)
!
!  19 points, precision 9.
!  Lyness and Jesperson.
!
  else if ( rule == 18 ) then

    a = 1.0D+00 / 3.0D+00
    b = 0.02063496160252593D+00
    c = 0.4896825191987370D+00
    d = 0.1258208170141290D+00
    e = 0.4370895914929355D+00
    f = 0.6235929287619356D+00
    g = 0.1882035356190322D+00
    r = 0.9105409732110941D+00
    s = 0.04472951339445297D+00
    t = 0.7411985987844980D+00
    u = 0.03683841205473626D+00
    v = 0.22196288916076574D+00

    w1 = 0.09713579628279610D+00
    w2 = 0.03133470022713983D+00
    w3 = 0.07782754100477543D+00
    w4 = 0.07964773892720910D+00
    w5 = 0.02557767565869810D+00
    w6 = 0.04328353937728940D+00

    xtab(1:19) =   (/  a,  b,  c,  c,  d,  e,  e,  f,  g,  g,  r,  s,  s, &
      t, t, u, u, v, v /)
    ytab(1:19) =   (/  a,  c,  b,  c,  e,  d,  e,  g,  f,  g,  s,  r,  s, &
      u, v, t, v, t, u /)
    weight(1:19) = (/ w1, w2, w2, w2, w3, w3, w3, w4, w4, w4, w5, w5, w5, &
      w6, w6, w6, w6, w6, w6 /)
!
!  28 points, precision 11.
!  Lyness and Jesperson.
!
  else if ( rule == 19 ) then

    a = 1.0D+00 / 3.0D+00
    b = 0.9480217181434233D+00
    c = 0.02598914092828833D+00
    d = 0.8114249947041546D+00
    e = 0.09428750264792270D+00
    f = 0.01072644996557060D+00
    g = 0.4946367750172147D+00
    p = 0.5853132347709715D+00
    q = 0.2073433826145142D+00
    r = 0.1221843885990187D+00
    s = 0.4389078057004907D+00
    t = 0.6779376548825902D+00
    u = 0.04484167758913055D+00
    v = 0.27722066752827925D+00
    w = 0.8588702812826364D+00
    x = 0.0D+00
    y = 0.1411297187173636D+00

    w1 = 0.08797730116222190D+00
    w2 = 0.008744311553736190D+00
    w3 = 0.03808157199393533D+00
    w4 = 0.01885544805613125D+00
    w5 = 0.07215969754474100D+00
    w6 = 0.06932913870553720D+00
    w7 = 0.04105631542928860D+00
    w8 = 0.007362383783300573D+00

    xtab(1:28) =   (/  a,  b,  c,  c,  d,  e,  e,  f,  g,  g,  p,  q,  q, &
       r,  s,  s,  t,  t,  u,  u,  v,  v,  w,  w,  x,  x,  y,  y /)
    ytab(1:28) =   (/  a,  c,  b,  c,  e,  d,  e,  g,  f,  g,  q,  p,  q, &
       s,  r,  s,  u,  v,  t,  v,  t,  u,  x,  y,  w,  y,  w,  x /)
    weight(1:28) = (/ w1, w2, w2, w2, w3, w3, w3, w4, w4, w4, w5, w5, w5, &
      w6, w6, w6, w7, w7, w7, w7, w7, w7, w8, w8, w8, w8, w8, w8 /)
!
!  37 points, precision 13.
!
  else if ( rule == 20 ) then

    a = 1.0D+00 / 3.0D+00
    b = 0.950275662924105565450352089520D+00
    c = 0.024862168537947217274823955239D+00
    d = 0.171614914923835347556304795551D+00
    e = 0.414192542538082326221847602214D+00
    f = 0.539412243677190440263092985511D+00
    g = 0.230293878161404779868453507244D+00

    w1 = 0.051739766065744133555179145422D+00
    w2 = 0.008007799555564801597804123460D+00
    w3 = 0.046868898981821644823226732071D+00
    w4 = 0.046590940183976487960361770070D+00
    w5 = 0.031016943313796381407646220131D+00
    w6 = 0.010791612736631273623178240136D+00
    w7 = 0.032195534242431618819414482205D+00
    w8 = 0.015445834210701583817692900053D+00
    w9 = 0.017822989923178661888748319485D+00
    wx = 0.037038683681384627918546472190D+00

    xtab(1:10) =   (/ a, b, c, c, d, e, e, f, g, g /)
    ytab(1:10) =   (/ a, c, b, c, e, d, e, g, f, g /)
    weight(1:37) = (/ w1, w2, w2, w2, w3, w3, w3, w4, w4, w4, w5, w5, w5, &
                      w6, w6, w6, w7, w7, w7, w8, w8, w8, w8, w8, w8, w9, &
                      w9, w9, w9, w9, w9, wx, wx, wx, wx, wx, wx /)

    a = 0.772160036676532561750285570113D+00
    b = 0.113919981661733719124857214943D+00

    xtab(11) = a
    ytab(11) = b

    xtab(12) = b
    ytab(12) = a

    xtab(13) = b
    ytab(13) = b

    a = 0.009085399949835353883572964740D+00
    b = 0.495457300025082323058213517632D+00

    xtab(14) = a
    ytab(14) = b

    xtab(15) = b
    ytab(15) = a

    xtab(16) = b
    ytab(16) = b

    a = 0.062277290305886993497083640527D+00
    b = 0.468861354847056503251458179727D+00

    xtab(17) = a
    ytab(17) = b

    xtab(18) = b
    ytab(18) = a

    xtab(19) = b
    ytab(19) = b

    a = 0.022076289653624405142446876931D+00
    b = 0.851306504174348550389457672223D+00
    c = 1.0D+00 - a - b

    xtab(20) = a
    ytab(20) = b

    xtab(21) = a
    ytab(21) = c

    xtab(22) = b
    ytab(22) = a

    xtab(23) = b
    ytab(23) = c

    xtab(24) = c
    ytab(24) = a

    xtab(25) = c
    ytab(25) = b

    a = 0.018620522802520968955913511549D+00
    b = 0.689441970728591295496647976487D+00
    c = 1.0D+00 - a - b

    xtab(26) = a
    ytab(26) = b

    xtab(27) = a
    ytab(27) = c

    xtab(28) = b
    ytab(28) = a

    xtab(29) = b
    ytab(29) = c

    xtab(30) = c
    ytab(30) = a

    xtab(31) = c
    ytab(31) = b

    a = 0.096506481292159228736516560903D+00
    b = 0.635867859433872768286976979827D+00
    c = 1.0D+00 - a - b

    xtab(32) = a
    ytab(32) = b

    xtab(33) = a
    ytab(33) = c

    xtab(34) = b
    ytab(34) = a

    xtab(35) = b
    ytab(35) = c

    xtab(36) = c
    ytab(36) = a

    xtab(37) = c
    ytab(37) = b

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRIANGLE_UNIT_SET - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of RULE = ', rule
    stop

  end if

  return
end
function triangle_unit_size ( rule )

!******************************************************************************
!
!! TRIANGLE_UNIT_SIZE returns the "size" of a unit triangle quadrature rule.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Gilbert Strang, George Fix,
!    An Analysis of the Finite Element Method,
!    Prentice Hall, 1973,
!    TA335.S77.
!
!    Olgierd Zienkiewicz,
!    The Finite Element Method,
!    McGraw Hill, Third Edition, 1977, page 202.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) RULE, the index of the rule.
!     1, ORDER =  1, precision 1, Zienkiewicz #1.
!     2, ORDER =  2, precision 1, (the "vertex rule").
!     3, ORDER =  3, precision 2, Strang and Fix formula #1.
!     4, ORDER =  3, precision 2, Strang and Fix formula #2, Zienkiewicz #2.
!     5, ORDER =  4, precision 3, Strang and Fix formula #3, Zienkiewicz #3.
!     6, ORDER =  6, precision 3, Strang and Fix formula #4.
!     7, ORDER =  6, precision 3, Stroud formula T2:3-1.
!     8, ORDER =  6, precision 4, Strang and Fix formula #5.
!     9, ORDER =  7, precision 4, Strang and Fix formula #6.
!    10, ORDER =  7, precision 5, Strang and Fix formula #7,
!        Stroud formula T2:5-1, Zienkiewicz #4, Schwarz Table 2.2.
!    11, ORDER =  9, precision 6, Strang and Fix formula #8.
!    12, ORDER = 12, precision 6, Strang and Fix formula #9.
!    13, ORDER = 13, precision 7, Strang and Fix formula #10.
!    14, ORDER =  7, precision ?.
!    15, ORDER = 16, precision 7, conical product Gauss, Stroud formula T2:7-1.
!    16, ORDER = 64, precision 15, triangular product Gauss rule.
!    17, ORDER = 19, precision 8, from CUBTRI, ACM TOMS #584.
!    18, ORDER = 19, precision 9, from TRIEX, Lyness and Jespersen.
!    19, ORDER = 28, precision 11, from TRIEX, Lyness and Jespersen.
!    20, ORDER = 37, precision 13, from ACM TOMS #706.
!
!    Output, integer ( kind = 4 ) TRIANGLE_UNIT_SIZE, the order of the rule.
!
  implicit none

  integer ( kind = 4 ) rule
  integer ( kind = 4 ) triangle_unit_size

  if ( rule == 1 ) then
    triangle_unit_size = 1
  else if ( rule == 2 ) then
    triangle_unit_size = 3
  else if ( rule == 3 ) then
    triangle_unit_size = 3
  else if ( rule == 4 ) then
    triangle_unit_size = 3
  else if ( rule == 5 ) then
    triangle_unit_size = 4
  else if ( rule == 6 ) then
    triangle_unit_size = 6
  else if ( rule == 7 ) then
    triangle_unit_size = 6
  else if ( rule == 8 ) then
    triangle_unit_size = 6
  else if ( rule == 9 ) then
    triangle_unit_size = 7
  else if ( rule == 10 ) then
    triangle_unit_size = 7
  else if ( rule == 11 ) then
    triangle_unit_size = 9
  else if ( rule == 12 ) then
    triangle_unit_size = 12
  else if ( rule == 13 ) then
    triangle_unit_size = 13
  else if ( rule == 14 ) then
    triangle_unit_size = 7
  else if ( rule == 15 ) then
    triangle_unit_size = 16
  else if ( rule == 16 ) then
    triangle_unit_size = 64
  else if ( rule == 17 ) then
    triangle_unit_size = 19
  else if ( rule == 18 ) then
    triangle_unit_size = 19
  else if ( rule == 19 ) then
    triangle_unit_size = 28
  else if ( rule == 20 ) then
    triangle_unit_size = 37
  else
    triangle_unit_size = -1
  end if

  return
end
