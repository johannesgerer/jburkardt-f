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
subroutine basis_brick8 ( n, p, phi )

!*****************************************************************************80
!
!! BASIS_BRICK8: BRICK8 basis functions at natural coordinates.
!
!  Discussion:
!
!      8------7        t  s
!     /|     /|        | /
!    5------6 |        |/
!    | |    | |        0-------r
!    | 4----|-3        
!    |/     |/        
!    1------2        
!                   
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 March 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) P(3,N), natural coordinates of evaluation
!    points.
!
!    Output, real ( kind = 8 ) PHI(8,N), the basis function values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) p(3,n)
  real ( kind = 8 ) phi(8,n)

  phi(1,1:n) = &
    ( 1.0D+00 - p(1,1:n) ) * ( 1.0D+00 - p(2,1:n) ) * ( 1.0D+00 - p(3,1:n) )
  phi(2,1:n) = &
    ( 1.0D+00 + p(1,1:n) ) * ( 1.0D+00 - p(2,1:n) ) * ( 1.0D+00 - p(3,1:n) )
  phi(3,1:n) = &
    ( 1.0D+00 + p(1,1:n) ) * ( 1.0D+00 + p(2,1:n) ) * ( 1.0D+00 - p(3,1:n) )
  phi(4,1:n) = &
    ( 1.0D+00 - p(1,1:n) ) * ( 1.0D+00 + p(2,1:n) ) * ( 1.0D+00 - p(3,1:n) )
  phi(5,1:n) = &
    ( 1.0D+00 - p(1,1:n) ) * ( 1.0D+00 - p(2,1:n) ) * ( 1.0D+00 + p(3,1:n) )
  phi(6,1:n) = &
    ( 1.0D+00 + p(1,1:n) ) * ( 1.0D+00 - p(2,1:n) ) * ( 1.0D+00 + p(3,1:n) )
  phi(7,1:n) = &
    ( 1.0D+00 + p(1,1:n) ) * ( 1.0D+00 + p(2,1:n) ) * ( 1.0D+00 + p(3,1:n) )
  phi(8,1:n) = &
    ( 1.0D+00 - p(1,1:n) ) * ( 1.0D+00 + p(2,1:n) ) * ( 1.0D+00 + p(3,1:n) )

  phi(1:8,1:n) = phi(1:8,1:n) / 8.0D+00

  return
end
subroutine basis_brick8_test ( )

!*****************************************************************************80
!
!! BASIS_BRICK8_TEST verifies BASIS_BRICK8.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 March 2010
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

  integer ( kind = 4 ), parameter :: node_num = 8

  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  real ( kind = 8 ), allocatable :: p(:,:)
  real ( kind = 8 ), allocatable :: phi(:,:)
  real ( kind = 8 ) phi_sum
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 5

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BASIS_BRICK8_TEST:'
  write ( *, '(a)' ) '  Verify basis functions for element BRICK8.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of nodes = ', node_num
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The basis function values at basis nodes'
  write ( *, '(a)' ) '  should form the identity matrix.'
  write ( *, '(a)' ) ' '

  n = node_num

  allocate ( p(1:3,1:n) )
  allocate ( phi(1:node_num,1:n) )

  call nodes_brick8 ( p )

  call basis_brick8 ( n, p, phi )

  do j = 1, n
    write ( *, '(2x,8f7.3)' ) phi(1:node_num,j)
  end do

  deallocate ( p )
  deallocate ( phi )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The basis function values at ANY point P'
  write ( *, '(a)' ) '  should sum to 1:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    ------------P-------------     PHI_SUM'
  write ( *, '(a)' ) ' '

  n = test_num

  allocate ( p(1:3,1:n) )
  allocate ( phi(1:node_num,1:n) )

  seed = 123456789

  call r8mat_uniform_01 ( 3, n, seed, p )

  call basis_brick8 ( n, p, phi )

  do j = 1, n

    phi_sum = sum ( phi(1:node_num,j) )

    write ( *, '(8(2x,f8.4))' ) p(1:3,j), phi_sum

  end do

  deallocate ( p )
  deallocate ( phi )

  return
end
subroutine basis_brick20 ( n, p, phi )

!*****************************************************************************80
!
!! BASIS_BRICK20: BRICK20 basis functions at natural coordinates.
!
!  Discussion:
!
!        8----19---7
!       /|        /|
!     20 |      18 |        t   s
!     /  16     /  15       |  /
!    5----17---6   |        | /
!    |   |     |   |        |/
!    |   4--11-|---3        0---------r
!   13  /     14  /        
!    | 12      | 10       
!    |/        |/        
!    1----9----2
!                   
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 March 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) P(3,N), natural coordinates of evaluation
!    points.
!
!    Output, real ( kind = 8 ) PHI(20,N), the basis function values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) p(3,n)
  real ( kind = 8 ) phi(20,n)

  phi(1,1:n) = &
    ( 1.0D+00 - p(1,1:n) ) * ( 1.0D+00 - p(2,1:n) ) * ( 1.0D+00 - p(3,1:n) ) &
    * ( - p(1,1:n) - p(2,1:n) - p(3,1:n) - 2.0D+00 ) / 8.0D+00
  phi(2,1:n) = &
    ( 1.0D+00 + p(1,1:n) ) * ( 1.0D+00 - p(2,1:n) ) * ( 1.0D+00 - p(3,1:n) ) &
    * ( + p(1,1:n) - p(2,1:n) - p(3,1:n) - 2.0D+00 ) / 8.0D+00
  phi(3,1:n) = &
    ( 1.0D+00 + p(1,1:n) ) * ( 1.0D+00 + p(2,1:n) ) * ( 1.0D+00 - p(3,1:n) ) &
    * ( + p(1,1:n) + p(2,1:n) - p(3,1:n) - 2.0D+00 ) / 8.0D+00
  phi(4,1:n) = &
    ( 1.0D+00 - p(1,1:n) ) * ( 1.0D+00 + p(2,1:n) ) * ( 1.0D+00 - p(3,1:n) ) &
    * ( - p(1,1:n) + p(2,1:n) - p(3,1:n) - 2.0D+00 ) / 8.0D+00
  phi(5,1:n) = &
    ( 1.0D+00 - p(1,1:n) ) * ( 1.0D+00 - p(2,1:n) ) * ( 1.0D+00 + p(3,1:n) ) &
    * ( - p(1,1:n) - p(2,1:n) + p(3,1:n) - 2.0D+00 ) / 8.0D+00
  phi(6,1:n) = &
    ( 1.0D+00 + p(1,1:n) ) * ( 1.0D+00 - p(2,1:n) ) * ( 1.0D+00 + p(3,1:n) ) &
    * ( + p(1,1:n) - p(2,1:n) + p(3,1:n) - 2.0D+00 ) / 8.0D+00
  phi(7,1:n) = &
    ( 1.0D+00 + p(1,1:n) ) * ( 1.0D+00 + p(2,1:n) ) * ( 1.0D+00 + p(3,1:n) ) &
    * ( + p(1,1:n) + p(2,1:n) + p(3,1:n) - 2.0D+00 ) / 8.0D+00
  phi(8,1:n) = &
    ( 1.0D+00 - p(1,1:n) ) * ( 1.0D+00 + p(2,1:n) ) * ( 1.0D+00 + p(3,1:n) ) &
    * ( - p(1,1:n) + p(2,1:n) + p(3,1:n) - 2.0D+00 ) / 8.0D+00

  phi(9,1:n) =  ( 1.0D+00 + p(1,1:n) ) * ( 1.0D+00 - p(1,1:n) ) &
              * ( 1.0D+00 - p(2,1:n) ) * ( 1.0D+00 - p(3,1:n) ) / 4.0D+00
  phi(10,1:n) = ( 1.0D+00 + p(1,1:n) ) * ( 1.0D+00 + p(2,1:n) ) &
              * ( 1.0D+00 - p(2,1:n) ) * ( 1.0D+00 - p(3,1:n) ) / 4.0D+00
  phi(11,1:n) = ( 1.0D+00 + p(1,1:n) ) * ( 1.0D+00 - p(1,1:n) ) &
              * ( 1.0D+00 + p(2,1:n) ) * ( 1.0D+00 - p(3,1:n) ) / 4.0D+00
  phi(12,1:n) = ( 1.0D+00 - p(1,1:n) ) * ( 1.0D+00 + p(2,1:n) ) &
              * ( 1.0D+00 - p(2,1:n) ) * ( 1.0D+00 - p(3,1:n) ) / 4.0D+00
  phi(13,1:n) = ( 1.0D+00 - p(1,1:n) ) * ( 1.0D+00 - p(2,1:n) ) &
              * ( 1.0D+00 + p(3,1:n) ) * ( 1.0D+00 - p(3,1:n) ) / 4.0D+00
  phi(14,1:n) = ( 1.0D+00 + p(1,1:n) ) * ( 1.0D+00 - p(2,1:n) ) &
              * ( 1.0D+00 + p(3,1:n) ) * ( 1.0D+00 - p(3,1:n) ) / 4.0D+00
  phi(15,1:n) = ( 1.0D+00 + p(1,1:n) ) * ( 1.0D+00 + p(2,1:n) ) &
              * ( 1.0D+00 + p(3,1:n) ) * ( 1.0D+00 - p(3,1:n) ) / 4.0D+00
  phi(16,1:n) = ( 1.0D+00 - p(1,1:n) ) * ( 1.0D+00 + p(2,1:n) ) &
              * ( 1.0D+00 + p(3,1:n) ) * ( 1.0D+00 - p(3,1:n) ) / 4.0D+00
  phi(17,1:n) = ( 1.0D+00 + p(1,1:n) ) * ( 1.0D+00 - p(1,1:n) ) &
              * ( 1.0D+00 - p(2,1:n) ) * ( 1.0D+00 + p(3,1:n) ) / 4.0D+00
  phi(18,1:n) = ( 1.0D+00 + p(1,1:n) ) * ( 1.0D+00 + p(2,1:n) ) &
              * ( 1.0D+00 - p(2,1:n) ) * ( 1.0D+00 + p(3,1:n) ) / 4.0D+00
  phi(19,1:n) = ( 1.0D+00 + p(1,1:n) ) * ( 1.0D+00 - p(1,1:n) ) &
              * ( 1.0D+00 + p(2,1:n) ) * ( 1.0D+00 + p(3,1:n) ) / 4.0D+00
  phi(20,1:n) = ( 1.0D+00 - p(1,1:n) ) * ( 1.0D+00 + p(2,1:n) ) &
              * ( 1.0D+00 - p(2,1:n) ) * ( 1.0D+00 + p(3,1:n) ) / 4.0D+00

  return
end
subroutine basis_brick20_test ( )

!*****************************************************************************80
!
!! BASIS_BRICK20_TEST verifies BASIS_BRICK20.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 March 2010
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

  integer ( kind = 4 ), parameter :: node_num = 20

  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  real ( kind = 8 ), allocatable :: p(:,:)
  real ( kind = 8 ), allocatable :: phi(:,:)
  real ( kind = 8 ) phi_sum
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 5

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BASIS_BRICK20_TEST:'
  write ( *, '(a)' ) '  Verify basis functions for element BRICK20.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of nodes = ', node_num
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The basis function values at basis nodes'
  write ( *, '(a)' ) '  should form the identity matrix.'
  write ( *, '(a)' ) ' '

  n = node_num

  allocate ( p(1:3,1:n) )
  allocate ( phi(1:node_num,1:n) )

  call nodes_brick20 ( p )

  call basis_brick20 ( n, p, phi )

  do j = 1, n
    write ( *, '(2x,20f7.3)' ) phi(1:node_num,j)
  end do

  deallocate ( p )
  deallocate ( phi )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The basis function values at ANY point P'
  write ( *, '(a)' ) '  should sum to 1:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    ------------P-------------     PHI_SUM'
  write ( *, '(a)' ) ' '

  n = test_num

  allocate ( p(1:3,1:n) )
  allocate ( phi(1:node_num,1:n) )

  seed = 123456789

  call r8mat_uniform_01 ( 3, n, seed, p )

  call basis_brick20 ( n, p, phi )

  do j = 1, n

    phi_sum = sum ( phi(1:node_num,j) )

    write ( *, '(8(2x,f8.4))' ) p(1:3,j), phi_sum

  end do

  deallocate ( p )
  deallocate ( phi )

  return
end
subroutine basis_brick27 ( n, p, phi )

!*****************************************************************************80
!
!! BASIS_BRICK20: BRICK20 basis functions at natural coordinates.
!
!  Discussion:
!
!        8----19---7
!       /|         /
!     20 |   26   /|
!     /          / |
!    5----17----6  |
!    |   |      |  |
!    |  16---24-|-15
!    |  /|      | /|
!    |25 |  27  |23|        t
!    |/         |/ |        |   s
!   13----22---14  |        |  /
!    |   |      |  |        | /
!    |   |      |  |        |/
!    |   4--11--|--3        0---------r
!    |  /       | /        
!    | 12   21  |10       
!    |/         |/        
!    1----9-----2
!                   
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 March 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) P(3,N), natural coordinates of evaluation
!    points.
!
!    Output, real ( kind = 8 ) PHI(27,N), the basis function values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) p(3,n)
  real ( kind = 8 ) phi(27,n)
  real ( kind = 8 ) rm(1:n)
  real ( kind = 8 ) rp(1:n)
  real ( kind = 8 ) rz(1:n)
  real ( kind = 8 ) sm(1:n)
  real ( kind = 8 ) sp(1:n)
  real ( kind = 8 ) sz(1:n)
  real ( kind = 8 ) tm(1:n)
  real ( kind = 8 ) tp(1:n)
  real ( kind = 8 ) tz(1:n)


  rm(1:n) = p(1,1:n) + 1.0D+00
  rz(1:n) = p(1,1:n)
  rp(1:n) = p(1,1:n) - 1.0D+00

  sm(1:n) = p(2,1:n) + 1.0D+00
  sz(1:n) = p(2,1:n)
  sp(1:n) = p(2,1:n) - 1.0D+00

  tm(1:n) = p(3,1:n) + 1.0D+00
  tz(1:n) = p(3,1:n)
  tp(1:n) = p(3,1:n) - 1.0D+00

  phi(1,1:n)  =        rz * rp      * sz * sp      * tz * tp / 8.0D+00
  phi(2,1:n)  =   rm * rz           * sz * sp      * tz * tp / 8.0D+00
  phi(3,1:n)  =   rm * rz      * sm * sz           * tz * tp / 8.0D+00
  phi(4,1:n)  =        rz * rp * sm * sz           * tz * tp / 8.0D+00
  phi(5,1:n)  =        rz * rp      * sz * sp * tm * tz      / 8.0D+00
  phi(6,1:n)  =   rm * rz           * sz * sp * tm * tz      / 8.0D+00
  phi(7,1:n)  =   rm * rz      * sm * sz      * tm * tz      / 8.0D+00
  phi(8,1:n)  =        rz * rp * sm * sz      * tm * tz      / 8.0D+00

  phi(9,1:n)  = - rm      * rp      * sz * sp      * tz * tp / 4.0D+00
  phi(10,1:n) = - rm * rz      * sm      * sp      * tz * tp / 4.0D+00
  phi(11,1:n) = - rm      * rp * sm * sz           * tz * tp / 4.0D+00
  phi(12,1:n) = -      rz * rp * sm      * sp      * tz * tp / 4.0D+00
  phi(13,1:n) = -      rz * rp      * sz * sp * tm      * tp / 4.0D+00
  phi(14,1:n) = - rm * rz           * sz * sp * tm      * tp / 4.0D+00
  phi(15,1:n) = - rm * rz      * sm * sz      * tm      * tp / 4.0D+00
  phi(16,1:n) = -      rz * rp * sm * sz      * tm      * tp / 4.0D+00
  phi(17,1:n) = - rm      * rp      * sz * sp * tm * tz      / 4.0D+00
  phi(18,1:n) = - rm * rz      * sm      * sp * tm * tz      / 4.0D+00
  phi(19,1:n) = - rm      * rp * sm * sz      * tm * tz      / 4.0D+00
  phi(20,1:n) = -      rz * rp * sm      * sp * tm * tz      / 4.0D+00

  phi(21,1:n) =   rm      * rp * sm      * sp      * tz * tp / 2.0D+00
  phi(22,1:n) =   rm      * rp      * sz * sp * tm      * tp / 2.0D+00
  phi(23,1:n) =   rm * rz      * sm      * sp * tm      * tp / 2.0D+00
  phi(24,1:n) =   rm      * rp * sm * sz      * tm      * tp / 2.0D+00
  phi(25,1:n) =        rz * rp * sm      * sp * tm      * tp / 2.0D+00
  phi(26,1:n) =   rm      * rp * sm      * sp * tm * tz      / 2.0D+00

  phi(27,1:n) = - rm      * rp * sm      * sp * tm      * tp

  return
end
subroutine basis_brick27_test ( )

!*****************************************************************************80
!
!! BASIS_BRICK27_TEST verifies BASIS_BRICK27.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 March 2010
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

  integer ( kind = 4 ), parameter :: node_num = 27

  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  real ( kind = 8 ), allocatable :: p(:,:)
  real ( kind = 8 ), allocatable :: phi(:,:)
  real ( kind = 8 ) phi_sum
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 5

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BASIS_BRICK27_TEST:'
  write ( *, '(a)' ) '  Verify basis functions for element BRICK27.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of nodes = ', node_num
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The basis function values at basis nodes'
  write ( *, '(a)' ) '  should form the identity matrix.'
  write ( *, '(a)' ) ' '

  n = node_num

  allocate ( p(1:3,1:n) )
  allocate ( phi(1:node_num,1:n) )

  call nodes_brick27 ( p )

  call basis_brick27 ( n, p, phi )

  do j = 1, n
    write ( *, '(2x,27f4.0)' ) phi(1:node_num,j)
  end do

  deallocate ( p )
  deallocate ( phi )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The basis function values at ANY point P'
  write ( *, '(a)' ) '  should sum to 1:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    ------------P-------------     PHI_SUM'
  write ( *, '(a)' ) ' '

  n = test_num

  allocate ( p(1:3,1:n) )
  allocate ( phi(1:node_num,1:n) )

  seed = 123456789

  call r8mat_uniform_01 ( 3, n, seed, p )

  call basis_brick27 ( n, p, phi )

  do j = 1, n

    phi_sum = sum ( phi(1:node_num,j) )

    write ( *, '(8(2x,f8.4))' ) p(1:3,j), phi_sum

  end do

  deallocate ( p )
  deallocate ( phi )

  return
end
subroutine basis_mn_tet4 ( t, n, p, phi )

!*****************************************************************************80
!
!! BASIS_MN_TET4: all bases at N points for a TET4 element.
!
!  Discussion:
!
!    The routine is given the coordinates of the vertices of a tetrahedron.
!
!    It works directly with these coordinates, and does not refer to a
!    reference element.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 August 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Olgierd Zienkiewicz,
!    The Finite Element Method,
!    Sixth Edition,
!    Butterworth-Heinemann, 2005,
!    ISBN: 0750663200,
!    LC: TA640.2.Z54.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T(3,4), the coordinates of the vertices.
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) P(3,N), the points where the basis functions
!    are to be evaluated.
!
!    Output, real ( kind = 8 ) PHI(4,N), the value of the basis functions
!    at the evaluation points.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) p(3,n)
  real ( kind = 8 ) phi(4,n)
  real ( kind = 8 ) t(3,4)
  real ( kind = 8 ) volume
!
!           | x1 x2 x3 x4 |
!  Volume = | y1 y2 y3 y4 |
!           | z1 z2 z3 z4 |
!           |  1  1  1  1 |
!
  volume =                             &
      t(1,1) * (                       &
        t(2,2) * ( t(3,3) - t(3,4) )   &
      - t(2,3) * ( t(3,2) - t(3,4) )   &
      + t(2,4) * ( t(3,2) - t(3,3) ) ) &
    - t(1,2) * (                       &
        t(2,1) * ( t(3,3) - t(3,4) )   &
      - t(2,3) * ( t(3,1) - t(3,4) )   &
      + t(2,4) * ( t(3,1) - t(3,3) ) ) &
    + t(1,3) * (                       &
        t(2,1) * ( t(3,2) - t(3,4) )   &
      - t(2,2) * ( t(3,1) - t(3,4) )   &
      + t(2,4) * ( t(3,1) - t(3,2) ) ) &
    - t(1,4) * (                       &
        t(2,1) * ( t(3,2) - t(3,3) )   &
      - t(2,2) * ( t(3,1) - t(3,3) )   &
      + t(2,3) * ( t(3,1) - t(3,2) ) )

  if ( volume == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BASIS_MN_TET4 - Fatal error!'
    write ( *, '(a)' ) '  Element has zero volume.'
    stop
  end if
!
!             | xp x2 x3 x4 |
!  Phi(1,P) = | yp y2 y3 y4 | / volume
!             | zp z2 z3 z4 |
!             |  1  1  1  1 |
!
  phi(1,1:n) = (                           &
      p(1,1:n) * (                         &
        t(2,2)   * ( t(3,3)   - t(3,4) )   &
      - t(2,3)   * ( t(3,2)   - t(3,4) )   &
      + t(2,4)   * ( t(3,2)   - t(3,3) ) ) &
    - t(1,2) * (                           &
        p(2,1:n) * ( t(3,3)   - t(3,4) )   &
      - t(2,3)   * ( p(3,1:n) - t(3,4) )   &
      + t(2,4)   * ( p(3,1:n) - t(3,3) ) ) &
    + t(1,3) * (                           &
        p(2,1:n) * ( t(3,2)   - t(3,4) )   &
      - t(2,2)   * ( p(3,1:n) - t(3,4) )   &
      + t(2,4)   * ( p(3,1:n) - t(3,2) ) ) &
    - t(1,4) * (                           &
        p(2,1:n) * ( t(3,2)   - t(3,3) )   &
      - t(2,2)   * ( p(3,1:n) - t(3,3) )   &
      + t(2,3)   * ( p(3,1:n) - t(3,2) ) ) ) / volume
!
!             | x1 xp x3 x4 |
!  Phi(2,P) = | y1 yp y3 y4 | / volume
!             | z1 zp z3 z4 |
!             |  1  1  1  1 |
!
  phi(2,1:n) = (                             &
      t(1,1) * (                             &
        p(2,1:n) * ( t(3,3)   - t(3,4) )     &
      - t(2,3)   * ( p(3,1:n) - t(3,4) )     &
      + t(2,4)   * ( p(3,1:n) - t(3,3) ) )   &
    - p(1,1:n)   * (                         &
        t(2,1)   * ( t(3,3)   - t(3,4) )     &
      - t(2,3)   * ( t(3,1)   - t(3,4) )     &
      + t(2,4)   * ( t(3,1)   - t(3,3) ) )   &
    + t(1,3) * (                             &
        t(2,1)   * ( p(3,1:n) - t(3,4) )     &
      - p(2,1:n) * ( t(3,1)   - t(3,4) )     &
      + t(2,4)   * ( t(3,1)   - p(3,1:n) ) ) &
    - t(1,4) * (                             &
        t(2,1)   * ( p(3,1:n) - t(3,3) )     &
      - p(2,1:n) * ( t(3,1)   - t(3,3) )     &
      + t(2,3)   * ( t(3,1)   - p(3,1:n) ) ) ) / volume
!
!             | x1 x2 xp x4 |
!  Phi(3,P) = | y1 y2 yp y4 | / volume
!             | z1 z2 zp z4 |
!             |  1  1  1  1 |
!
  phi(3,1:n) = (                             & 
      t(1,1) * (                             &
        t(2,2)   * ( p(3,1:n) - t(3,4) )     &
      - p(2,1:n) * ( t(3,2)   - t(3,4) )     &
      + t(2,4)   * ( t(3,2)   - p(3,1:n) ) ) &
    - t(1,2) * (                             &
        t(2,1)   * ( p(3,1:n) - t(3,4) )     &
      - p(2,1:n) * ( t(3,1)   - t(3,4) )     &
      + t(2,4)   * ( t(3,1)   - p(3,1:n) ) ) &
    + p(1,1:n) * (                           &
        t(2,1)   * ( t(3,2)   - t(3,4) )     &
      - t(2,2)   * ( t(3,1)   - t(3,4) )     &
      + t(2,4)   * ( t(3,1)   - t(3,2) ) )   &
    - t(1,4) * (                             &
        t(2,1)   * ( t(3,2)   - p(3,1:n) )   &
      - t(2,2)   * ( t(3,1)   - p(3,1:n) )   &
      + p(2,1:n) * ( t(3,1)   - t(3,2) ) ) ) / volume
!
!             | x1 x2 x3 xp |
!  Phi(4,P) = | y1 y2 y3 yp | / volume
!             | z1 z2 z3 zp |
!             |  1  1  1  1 |
!
  phi(4,1:n) = (                             &
      t(1,1) * (                             &
        t(2,2)   * ( t(3,3)   - p(3,1:n) )   &
      - t(2,3)   * ( t(3,2)   - p(3,1:n) )   &
      + p(2,1:n) * ( t(3,2)   - t(3,3) ) )   &
    - t(1,2) * (                             &
        t(2,1)   * ( t(3,3)   - p(3,1:n) )   &
      - t(2,3)   * ( t(3,1)   - p(3,1:n) )   &
      + p(2,1:n) * ( t(3,1)   - t(3,3) ) )   &
    + t(1,3) * (                             &
        t(2,1)   * ( t(3,2)   - p(3,1:n) )   &
      - t(2,2)   * ( t(3,1)   - p(3,1:n) )   &
      + p(2,1:n) * ( t(3,1)   - t(3,2) ) )   &
    - p(1,1:n) * (                           &
        t(2,1)   * ( t(3,2)   - t(3,3) )     &
      - t(2,2)   * ( t(3,1)   - t(3,3) )     &
      + t(2,3)   * ( t(3,1)   - t(3,2) ) ) ) / volume

  return
end
subroutine basis_mn_tet4_test ( )

!*****************************************************************************80
!
!! BASIS_MN_TET4_TEST verifies BASIS_MN_TET4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 August 2009
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

  real ( kind = 8 ) c(4)
  real ( kind = 8 ) c_sum
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) p(3)
  real ( kind = 8 ) phi1(node_num,1)
  real ( kind = 8 ) phi1_sum
  real ( kind = 8 ) phi4(node_num,4)
  integer ( kind = 4 ) seed
  real ( kind = 8 ), dimension ( 3, node_num ) :: t
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 5

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BASIS_MN_TET4_TEST:'
  write ( *, '(a)' ) '  Verify basis functions for element TET4.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of nodes = ', node_num

  call r8mat_uniform_01 ( 3, 4, seed, t )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Tetrahedron Nodes:'
  write ( *, '(a)' ) ' '
  do j = 1, node_num
    write ( *, '(2x,i8,3(2x,f7.3))' ) j, t(1:3,j)
  end do
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The basis function values at basis nodes'
  write ( *, '(a)' ) '  should form the identity matrix.'
  write ( *, '(a)' ) ' '

  call basis_mn_tet4 ( t, node_num, t, phi4 )

  do i = 1, node_num
    write ( *, '(2x,10f7.3)' ) phi4(i,1:node_num)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The basis function values at ANY point P'
  write ( *, '(a)' ) '  should sum to 1:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    ------------P-------------    ' // &
    '-----------------PHI----------------   PHI_SUM'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    call r8vec_uniform_01 ( 4, seed, c )

    c_sum = sum ( c(1:4) )
    c(1:4) = c(1:4) / c_sum
    p(1:3) = matmul ( t(1:3,1:4), c(1:4) )

    call basis_mn_tet4 ( t, 1, p, phi1 )
    phi1_sum = sum ( phi1(1:node_num,1) )

    write ( *, '(8(2x,f8.4))' ) p(1:3), phi1(1:node_num,1), phi1_sum

  end do

  return
end
subroutine basis_mn_tet10 ( t, n, p, phi )

!*****************************************************************************80
!
!! BASIS_MN_TET10: all bases at N points for a TET10 element.
!
!  Discussion:
!
!    The routine is given the coordinates of the vertices of a tetrahedron.
!
!    It works directly with these coordinates, and does not refer to a
!    reference element.
!
!    P1 through P4 are vertices.
!
!    P1 <= P5  <= P2
!    P2 <= P6  <= P3
!    P1 <= P7  <= P3
!    P1 <= P8  <= P4
!    P2 <= P9  <= P4
!    P3 <= P10 <= P4
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 August 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Olgierd Zienkiewicz,
!    The Finite Element Method,
!    Sixth Edition,
!    Butterworth-Heinemann, 2005,
!    ISBN: 0750663200,
!    LC: TA640.2.Z54.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T(3,4), the coordinates of the vertices.
!
!    Input, integer ( kind = 4 ) N, the number of evaluation points.
!
!    Input, real ( kind = 8 ) P(3,N), the points where the basis functions
!    are to be evaluated.
!
!    Output, real ( kind = 8 ) PHI(10,N), the value of the basis functions
!    at the evaluation points.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) p(3,n)
  real ( kind = 8 ) phi(10,n)
  real ( kind = 8 ) phi_linear(4,n)
  real ( kind = 8 ) t(3,4)

  call basis_mn_tet4 ( t, n, p, phi_linear )

  phi( 1,1:n) = ( 2.0D+00 * phi_linear(1,1:n)  - 1.0D+00 ) * phi_linear(1,1:n)
  phi( 2,1:n) = ( 2.0D+00 * phi_linear(2,1:n)  - 1.0D+00 ) * phi_linear(2,1:n)
  phi( 3,1:n) = ( 2.0D+00 * phi_linear(3,1:n)  - 1.0D+00 ) * phi_linear(3,1:n)
  phi( 4,1:n) = ( 2.0D+00 * phi_linear(4,1:n)  - 1.0D+00 ) * phi_linear(4,1:n)
  phi( 5,1:n) =   4.0D+00 * phi_linear(1,1:n)              * phi_linear(2,1:n)
  phi( 6,1:n) =   4.0D+00 * phi_linear(2,1:n)              * phi_linear(3,1:n)
  phi( 7,1:n) =   4.0D+00 * phi_linear(1,1:n)              * phi_linear(3,1:n)
  phi( 8,1:n) =   4.0D+00 * phi_linear(1,1:n)              * phi_linear(4,1:n)
  phi( 9,1:n) =   4.0D+00 * phi_linear(2,1:n)              * phi_linear(4,1:n)
  phi(10,1:n) =   4.0D+00 * phi_linear(3,1:n)              * phi_linear(4,1:n)

  return
end
subroutine basis_mn_tet10_test ( )

!*****************************************************************************80
!
!! BASIS_MN_TET10_TEST verifies BASIS_MN_TET10.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 August 2009
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

  real ( kind = 8 ) c(4)
  real ( kind = 8 ) c_sum
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) p(3)
  real ( kind = 8 ) p10(3,10)
  real ( kind = 8 ) phi1(10,1)
  real ( kind = 8 ) phi1_sum
  real ( kind = 8 ) phi10(10,10)
  integer ( kind = 4 ) seed
  real ( kind = 8 ), dimension ( 3, 4 ) :: t
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 5

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BASIS_MN_TET10_TEST:'
  write ( *, '(a)' ) '  Verify basis functions for element TET10.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of nodes = 10.'

  call r8mat_uniform_01 ( 3, 4, seed, t )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Tetrahedron Nodes:'
  write ( *, '(a)' ) ' '
  do j = 1, 4
    write ( *, '(2x,i8,3(2x,g14.6))' ) j, t(1:3,j)
  end do
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The basis function values at basis nodes'
  write ( *, '(a)' ) '  should form the identity matrix.'
  write ( *, '(a)' ) ' '

  p10(1:3,1:4) = t(1:3,1:4)
  p10(1:3, 5)  = 0.5D+00 * ( t(1:3,1) + t(1:3,2) )
  p10(1:3, 6)  = 0.5D+00 * ( t(1:3,2) + t(1:3,3) )
  p10(1:3, 7)  = 0.5D+00 * ( t(1:3,1) + t(1:3,3) )
  p10(1:3, 8)  = 0.5D+00 * ( t(1:3,1) + t(1:3,4) )
  p10(1:3, 9)  = 0.5D+00 * ( t(1:3,2) + t(1:3,4) )
  p10(1:3,10)  = 0.5D+00 * ( t(1:3,3) + t(1:3,4) )

  call basis_mn_tet10 ( t, 10, p10, phi10 )

  do i = 1, 10
    write ( *, '(2x,10f7.3)' ) phi10(i,1:10)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The basis function values at ANY point P'
  write ( *, '(a)' ) '  should sum to 1:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    ------------P-------------    ' // &
    '----------------------------------------------------' // &
    'PHI-----------------------------------------   PHI_SUM'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    call r8vec_uniform_01 ( 4, seed, c )

    c_sum = sum ( c(1:4) )
    c(1:4) = c(1:4) / c_sum
    p(1:3) = matmul ( t(1:3,1:4), c(1:4) )

    call basis_mn_tet10 ( t, 1, p, phi1 )
    phi1_sum = sum ( phi1(1:10,1) )

    write ( *, '(14(2x,f8.4))' ) p(1:3), phi1(1:10,1), phi1_sum

  end do

  return
end
subroutine nodes_brick8 ( p )

!*****************************************************************************80
!
!! NODES_BRICK8 returns the natural coordinates of the BRICK8 element.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 March 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) P(3,8), the coordinates.
!
  implicit none

  real ( kind = 8 ) p(3,8)
  real ( kind = 8 ), save :: p_save(3,8) = reshape ( (/ &
    -1.0D+00, -1.0D+00, -1.0D+00, &
    +1.0D+00, -1.0D+00, -1.0D+00, &
    +1.0D+00, +1.0D+00, -1.0D+00, &
    -1.0D+00, +1.0D+00, -1.0D+00, &
    -1.0D+00, -1.0D+00, +1.0D+00, &
    +1.0D+00, -1.0D+00, +1.0D+00, &
    +1.0D+00, +1.0D+00, +1.0D+00, &
    -1.0D+00, +1.0D+00, +1.0D+00 /), (/ 3, 8 /) )

  p(1:3,1:8) = p_save(1:3,1:8)

  return
end
subroutine nodes_brick20 ( p )

!*****************************************************************************80
!
!! NODES_BRICK20 returns the natural coordinates of the BRICK20 element.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 March 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) P(3,20), the coordinates.
!
  implicit none

  real ( kind = 8 ) p(3,20)
  real ( kind = 8 ), save :: p_save(3,20) = reshape ( (/ &
    -1.0D+00, -1.0D+00, -1.0D+00, &
    +1.0D+00, -1.0D+00, -1.0D+00, &
    +1.0D+00, +1.0D+00, -1.0D+00, &
    -1.0D+00, +1.0D+00, -1.0D+00, &
    -1.0D+00, -1.0D+00, +1.0D+00, &
    +1.0D+00, -1.0D+00, +1.0D+00, &
    +1.0D+00, +1.0D+00, +1.0D+00, &
    -1.0D+00, +1.0D+00, +1.0D+00, &
     0.0D+00, -1.0D+00, -1.0D+00, &
    +1.0D+00,  0.0D+00, -1.0D+00, &
     0.0D+00, +1.0D+00, -1.0D+00, &
    -1.0D+00,  0.0D+00, -1.0D+00, &
    -1.0D+00, -1.0D+00,  0.0D+00, &
    +1.0D+00, -1.0D+00,  0.0D+00, &
    +1.0D+00, +1.0D+00,  0.0D+00, &
    -1.0D+00, +1.0D+00,  0.0D+00, &
     0.0D+00, -1.0D+00, +1.0D+00, &
    +1.0D+00,  0.0D+00, +1.0D+00, &
     0.0D+00, +1.0D+00, +1.0D+00, &
    -1.0D+00,  0.0D+00, +1.0D+00 /), (/ 3, 20 /) )

  p(1:3,1:20) = p_save(1:3,1:20)

  return
end
subroutine nodes_brick27 ( p )

!*****************************************************************************80
!
!! NODES_BRICK27 returns the natural coordinates of the BRICK27 element.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 March 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) P(3,27), the coordinates.
!
  implicit none

  real ( kind = 8 ) p(3,27)
  real ( kind = 8 ), save :: p_save(3,27) = reshape ( (/ &
    -1.0D+00, -1.0D+00, -1.0D+00, &
    +1.0D+00, -1.0D+00, -1.0D+00, &
    +1.0D+00, +1.0D+00, -1.0D+00, &
    -1.0D+00, +1.0D+00, -1.0D+00, &
    -1.0D+00, -1.0D+00, +1.0D+00, &
    +1.0D+00, -1.0D+00, +1.0D+00, &
    +1.0D+00, +1.0D+00, +1.0D+00, &
    -1.0D+00, +1.0D+00, +1.0D+00, &
     0.0D+00, -1.0D+00, -1.0D+00, &
    +1.0D+00,  0.0D+00, -1.0D+00, &
     0.0D+00, +1.0D+00, -1.0D+00, &
    -1.0D+00,  0.0D+00, -1.0D+00, &
    -1.0D+00, -1.0D+00,  0.0D+00, &
    +1.0D+00, -1.0D+00,  0.0D+00, &
    +1.0D+00, +1.0D+00,  0.0D+00, &
    -1.0D+00, +1.0D+00,  0.0D+00, &
     0.0D+00, -1.0D+00, +1.0D+00, &
    +1.0D+00,  0.0D+00, +1.0D+00, &
     0.0D+00, +1.0D+00, +1.0D+00, &
    -1.0D+00,  0.0D+00, +1.0D+00, &
     0.0D+00,  0.0D+00, -1.0D+00, &
     0.0D+00, -1.0D+00,  0.0D+00, &
    +1.0D+00,  0.0D+00,  0.0D+00, &
     0.0D+00, +1.0D+00,  0.0D+00, &
    -1.0D+00,  0.0D+00,  0.0D+00, &
     0.0D+00,  0.0D+00, +1.0D+00, &
     0.0D+00,  0.0D+00,  0.0D+00 /), (/ 3, 27 /) )

  p(1:3,1:27) = p_save(1:3,1:27)

  return
end
subroutine physical_to_reference_tet4 ( t, n, phy, ref )

!*****************************************************************************80
!
!! PHYSICAL_TO_REFERENCE_TET4 maps physical points to reference points.
!
!  Discussion:
!
!    Given the vertices of an order 4 physical tetrahedron and a point 
!    (X,Y,Z) in the physical tetrahedron, the routine computes the value 
!    of the corresponding point (R,S,T) in the reference tetrahedron.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 August 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T(3,4), the coordinates of the vertices of the
!    physical tetrahedron.  The vertices are assumed to be the images of
!    (1,0,0), (0,1,0), (0,0,1) and (0,0,0) respectively.
!
!    Input, integer ( kind = 4 ) N, the number of points to transform.
!
!    Input, real ( kind = 8 ) PHY(3,N), the coordinates of physical points
!    to be transformed.
!
!    Output, real ( kind = 8 ) REF(3,N), the coordinates of the corresponding
!    points in the reference tetrahedron.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(3,3)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  real ( kind = 8 ) phy(3,n)
  real ( kind = 8 ) ref(3,n)
  real ( kind = 8 ) t(3,4)

  a(1:3,1:3) = t(1:3,1:3)
  do i = 1, 3
    a(i,1:3) = a(i,1:3) - t(i,4)
  end do

  do i = 1, 3
    ref(i,1:n) = phy(i,1:n) - t(i,4)
  end do

  call r8ge_fss ( 3, a, n, ref, info )

  return
end
subroutine r8ge_fss ( n, a, nb, b, info )

!*****************************************************************************80
!
!! R8GE_FSS factors and solves multiple R8GE systems.
!
!  Discussion:
!
!    The R8GE storage format is used for a general M by N matrix.  A storage 
!    space is made for each entry.  The two dimensional logical
!    array can be thought of as a vector of M*N entries, starting with
!    the M entries in the column 1, then the M entries in column 2
!    and so on.  Considered as a vector, the entry A(I,J) is then stored
!    in vector location I+(J-1)*M.
!
!    R8GE storage is used by LINPACK and LAPACK.
!
!    This routine does not save the LU factors of the matrix, and hence cannot
!    be used to efficiently solve multiple linear systems, or even to
!    factor A at one time, and solve a single linear system at a later time.
!
!    This routine uses partial pivoting, but no pivot vector is required.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 June 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!    N must be positive.
!
!    Input/output, real ( kind = 8 ) A(N,N).
!    On input, A is the coefficient matrix of the linear system.
!    On output, A is in unit upper triangular form, and
!    represents the U factor of an LU factorization of the
!    original coefficient matrix.
!
!    Input, integer ( kind = 4 ) NB, the number of right hand sides.
!
!    Input/output, real ( kind = 8 ) B(N,NB).
!    On input, the right hand sides of the linear system.
!    On output, the solutions of the linear systems.
!
!    Output, integer ( kind = 4 ) INFO, singularity flag.
!    0, no singularity detected.
!    nonzero, the factorization failed on the INFO-th step.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nb

  real ( kind = 8 ) a(n,n)
  real ( kind = 8 ) b(n,nb)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipiv
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jcol
  real ( kind = 8 ) piv
  real ( kind = 8 ) row(n)
  real ( kind = 8 ) t(nb)
  real ( kind = 8 ) temp

  info = 0

  do jcol = 1, n
!
!  Find the maximum element in column I.
!
    piv = abs ( a(jcol,jcol) )
    ipiv = jcol
    do i = jcol + 1, n
      if ( piv < abs ( a(i,jcol) ) ) then
        piv = abs ( a(i,jcol) )
        ipiv = i
      end if
    end do

    if ( piv == 0.0D+00 ) then
      info = jcol
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8GE_FSS - Fatal error!'
      write ( *, '(a,i8)' ) '  Zero pivot on step ', info
      stop
    end if
!
!  Switch rows JCOL and IPIV, and B.
!
    if ( jcol /= ipiv ) then

      row(1:n) = a(jcol,1:n)
      a(jcol,1:n) = a(ipiv,1:n)
      a(ipiv,1:n) = row(1:n)

      t(1:nb)      = b(jcol,1:nb)
      b(jcol,1:nb) = b(ipiv,1:nb)
      b(ipiv,1:nb) = t(1:nb)

    end if
!
!  Scale the pivot row.
!
    a(jcol,jcol+1:n) = a(jcol,jcol+1:n) / a(jcol,jcol)
    b(jcol,1:nb) = b(jcol,1:nb) / a(jcol,jcol)
    a(jcol,jcol) = 1.0D+00
!
!  Use the pivot row to eliminate lower entries in that column.
!
    do i = jcol + 1, n
      if ( a(i,jcol) /= 0.0D+00 ) then
        temp = - a(i,jcol)
        a(i,jcol) = 0.0D+00
        a(i,jcol+1:n) = a(i,jcol+1:n) + temp * a(jcol,jcol+1:n)
        b(i,1:nb) = b(i,1:nb) + temp * b(jcol,1:nb)
      end if
    end do

  end do
!
!  Back solve.
!
  do j = 1, nb
    do jcol = n, 2, -1
      b(1:jcol-1,j) = b(1:jcol-1,j) - a(1:jcol-1,jcol) * b(jcol,j)
    end do
  end do

  return
end
function r8mat_det_4d ( a )

!*****************************************************************************80
!
!! R8MAT_DET_4D computes the determinant of a 4 by 4 R8MAT.
!
!  Discussion:
!
!    An R8MAT is a two dimensional matrix of double precision real values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A(4,4), the matrix whose determinant is desired.
!
!    Output, real ( kind = 8 ) R8MAT_DET_4D, the determinant of the matrix.
!
  implicit none

  real ( kind = 8 ) a(4,4)
  real ( kind = 8 ) r8mat_det_4d

  r8mat_det_4d = &
         a(1,1) * ( &
             a(2,2) * ( a(3,3) * a(4,4) - a(3,4) * a(4,3) ) &
           - a(2,3) * ( a(3,2) * a(4,4) - a(3,4) * a(4,2) ) &
           + a(2,4) * ( a(3,2) * a(4,3) - a(3,3) * a(4,2) ) ) &
       - a(1,2) * ( &
             a(2,1) * ( a(3,3) * a(4,4) - a(3,4) * a(4,3) ) &
           - a(2,3) * ( a(3,1) * a(4,4) - a(3,4) * a(4,1) ) &
           + a(2,4) * ( a(3,1) * a(4,3) - a(3,3) * a(4,1) ) ) &
       + a(1,3) * ( &
             a(2,1) * ( a(3,2) * a(4,4) - a(3,4) * a(4,2) ) &
           - a(2,2) * ( a(3,1) * a(4,4) - a(3,4) * a(4,1) ) &
           + a(2,4) * ( a(3,1) * a(4,2) - a(3,2) * a(4,1) ) ) &
       - a(1,4) * ( &
             a(2,1) * ( a(3,2) * a(4,3) - a(3,3) * a(4,2) ) &
           - a(2,2) * ( a(3,1) * a(4,3) - a(3,3) * a(4,1) ) &
           + a(2,3) * ( a(3,1) * a(4,2) - a(3,2) * a(4,1) ) )

  return
end
subroutine r8mat_solve ( n, rhs_num, a, info )

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
!    06 August 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) RHS_NUM, the number of right hand sides.  
!    RHS_NUM must be at least 0.
!
!    Input/output, real ( kind = 8 ) A(N,N+RHS_NUM), contains in rows and
!    columns 1 to N the coefficient matrix, and in columns N+1 through
!    N+rhs_num, the right hand sides.  On output, the coefficient matrix
!    area has been destroyed, while the right hand sides have
!    been overwritten with the corresponding solutions.
!
!    Output, integer ( kind = 4 ) INFO, singularity flag.
!    0, the matrix was not singular, the solutions were computed;
!    J, factorization failed on step J, and the solutions could not
!    be computed.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) rhs_num

  real ( kind = 8 ) a(n,n+rhs_num)
  real ( kind = 8 ) apivot
  real ( kind = 8 ) factor
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipivot
  integer ( kind = 4 ) j
  real ( kind = 8 ) t(n+rhs_num)

  info = 0

  do j = 1, n
!
!  Choose a pivot row.
!
    ipivot = j
    apivot = a(j,j)

    do i = j + 1, n
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
!  The pivot row moves into the J-th row.
!
    if ( ipivot /= j ) then
      t(       1:n+rhs_num) = a(ipivot,1:n+rhs_num)
      a(ipivot,1:n+rhs_num) = a(j,     1:n+rhs_num)
      a(j,     1:n+rhs_num) = t(       1:n+rhs_num)
    end if
!
!  A(J,J) becomes 1.
!
    a(j,j) = 1.0D+00
    a(j,j+1:n+rhs_num) = a(j,j+1:n+rhs_num) / apivot
!
!  A(I,J) becomes 0.
!
    do i = 1, n

      if ( i /= j ) then
        factor = a(i,j)
        a(i,j) = 0.0D+00
        a(i,j+1:n+rhs_num) = a(i,j+1:n+rhs_num) - factor * a(j,j+1:n+rhs_num)
      end if

    end do

  end do

  return
end
subroutine r8mat_uniform_01 ( m, n, seed, r )

!*****************************************************************************80
!
!! R8MAT_UNIFORM_01 returns a unit pseudorandom R8MAT.
!
!  Discussion:
!
!    An R8MAT is an array of R8's.
!
!    For now, the input quantity SEED is an integer variable.
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
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns 
!    in the array.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which 
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R(M,N), the array of pseudorandom values.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real ( kind = 8 ) r(m,n)

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  do j = 1, n

    do i = 1, m

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed < 0 ) then
        seed = seed + i4_huge
      end if

      r(i,j) = real ( seed, kind = 8 ) * 4.656612875D-10

    end do
  end do

  return
end
subroutine r8vec_uniform_01 ( n, seed, r )

!*****************************************************************************80
!
!! R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    For now, the input quantity SEED is an integer variable.
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
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which 
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R(N), the vector of pseudorandom values.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real ( kind = 8 ) r(n)

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + i4_huge
    end if

    r(i) = real ( seed, kind = 8 ) * 4.656612875D-10

  end do

  return
end
subroutine reference_tet4_sample ( n, seed, ref )

!*****************************************************************************80
!
!! REFERENCE_TET4_SAMPLE: sample points in the reference tetrahedron.
!
!  Discussion:
!
!    This sampling method is not uniform.  The algorithm is simple.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 August 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) REF(3,N), points in the reference tetrahedron.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c(4)
  real ( kind = 8 ) c_sum
  integer ( kind = 4 ) j
  real ( kind = 8 ) ref(3,n)
  integer ( kind = 4 ) seed

  do j = 1, n
    call r8vec_uniform_01 ( 4, seed, c )
    c_sum = sum ( c(1:4) )
    ref(1:3,j) = c(1:3) / c_sum
  end do

  return
end
subroutine reference_tet4_uniform ( n, seed, x )

!*****************************************************************************80
!
!! REFERENCE_TET4_UNIFORM: uniform sample points in the reference tetrahedron.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 August 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Reuven Rubinstein,
!    Monte Carlo Optimization, Simulation, and Sensitivity 
!    of Queueing Networks,
!    Krieger, 1992,
!    ISBN: 0894647644,
!    LC: QA298.R79.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) X(3,N), the points.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) e(4)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x(3,n)

  do j = 1, n

    call r8vec_uniform_01 ( 4, seed, e )

    e(1:4) = - log ( e(1:4) )

    x(1:3,j) = e(1:3) / sum ( e(1:4) )

  end do

  return
end
subroutine reference_tet4_uniform2 ( n, seed, x )

!*****************************************************************************80
!
!! REFERENCE_TET4_UNIFORM2: uniform sample points in the reference tetrahedron.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 August 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Claudio Rocchini, Paolo Cignoni,
!    Generating Random Points in a Tetrahedron,
!    Journal of Graphics Tools,
!    Volume 5, Number 5, 2000, pages 9-12.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) X(3,N), the points.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c(4)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) seed
  real ( kind = 8 ) t
  real ( kind = 8 ) x(3,n)

  do j = 1, n

    call r8vec_uniform_01 ( 3, seed, c )

    if ( 1.0D+00 < c(1) + c(2) ) then
      c(1) = 1.0D+00 - c(1)
      c(2) = 1.0D+00 - c(2)
    end if

    if ( 1.0D+00 < c(2) + c(3) ) then
      t = c(3)
      c(3) = 1.0D+00 - c(1) - c(2)
      c(2) = 1.0D+00 - t
    else if ( 1.0D+00 < c(1) + c(2) + c(3) ) then
       t = c(3)
       c(3) = c(1) + c(2) + c(3) - 1.0D+00
       c(1) = 1.0D+00 - c(2) - t
    end if

    c(4) = 1.0D+00 - c(1) - c(2) - c(3)
!
!  C(1:4) are the barycentric coordinates of the point.
!
    x(1:3,j) = c(1:3)

  end do

  return
end
subroutine reference_to_physical_tet4 ( t, n, ref, phy )

!*****************************************************************************80
!
!! REFERENCE_TO_PHYSICAL_TET4 maps TET4 reference points to physical points.
!
!  Discussion:
!
!    Given the vertices of an order 4 physical tetrahedron and a point 
!    (R,S,T) in the reference tetrahedron, the routine computes the value 
!    of the corresponding point (X,Y,Z) in the physical tetrahedron.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 August 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T(3,4), the coordinates of the vertices.  
!    The vertices are assumed to be the images of (1,0,0), (0,1,0),
!    (0,0,1) and (0,0,0) respectively.
!
!    Input, integer ( kind = 4 ) N, the number of points to transform.
!
!    Input, real ( kind = 8 ) REF(3,N), points in the reference tetrahedron.
!
!    Output, real ( kind = 8 ) PHY(3,N), corresponding points in the
!    physical tetrahedron.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real ( kind = 8 ) phy(3,n)
  real ( kind = 8 ) ref(3,n)
  real ( kind = 8 ) t(3,4)

  do i = 1, 3
    phy(i,1:n) =                                                    &
        t(i,1) *             ref(1,1:n)                             &
      + t(i,2) *                          ref(2,1:n)                &
      + t(i,3) *                                       ref(3,1:n)   &
      + t(i,4) * ( 1.0D+00 - ref(1,1:n) - ref(2,1:n) - ref(3,1:n) ) 
  end do

  return
end
subroutine tetrahedron_barycentric ( tetra, p, c )

!*****************************************************************************80
!
!! TETRAHEDRON_BARYCENTRIC: barycentric coordinates of a point.
!
!  Discussion:
!
!    The barycentric coordinates of a point P with respect to
!    a tetrahedron are a set of four values C(1:4), each associated
!    with a vertex of the tetrahedron.  The values must sum to 1.
!    If all the values are between 0 and 1, the point is contained
!    within the tetrahedron.
!
!    The barycentric coordinate of point P related to vertex A can be
!    interpreted as the ratio of the volume of the tetrahedron with 
!    vertex A replaced by vertex P to the volume of the original 
!    tetrahedron.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) TETRA(3,4) the tetrahedron vertices.
!
!    Input, real ( kind = 8 ) P(3), the point to be checked.
!
!    Output, real ( kind = 8 ) C(4), the barycentric coordinates of P with
!    respect to the tetrahedron.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3
  integer ( kind = 4 ), parameter :: rhs_num = 1

  real ( kind = 8 ) a(dim_num,dim_num+rhs_num)
  real ( kind = 8 ) c(dim_num+1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  real ( kind = 8 ) p(dim_num)
  real ( kind = 8 ) tetra(dim_num,4)
!
!  Set up the linear system
!
!    ( X2-X1  X3-X1  X4-X1 ) C2    X - X1
!    ( Y2-Y1  Y3-Y1  Y4-Y1 ) C3  = Y - Y1
!    ( Z2-Z1  Z3-Z1  Z4-Z1 ) C4    Z - Z1
!
!  which is satisfied by the barycentric coordinates of P.
!
  a(1:dim_num,1:3) = tetra(1:dim_num,2:4)
  a(1:dim_num,4) = p(1:dim_num)

  do i = 1, dim_num
    a(i,1:4) = a(i,1:4) - tetra(i,1)
  end do
!
!  Solve the linear system.
!
  call r8mat_solve ( dim_num, rhs_num, a, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TETRAHEDRON_BARYCENTRIC - Fatal error!'
    write ( *, '(a)' ) '  The linear system is singular.'
    write ( *, '(a)' ) '  The input data does not form a proper tetrahedron.'
    stop
  end if

  c(2:4) = a(1:dim_num,4)

  c(1) = 1.0D+00 - sum ( c(2:4) )

  return
end
subroutine tetrahedron_volume ( tet_xyz, volume )

!*****************************************************************************80
!
!! TETRAHEDRON_VOLUME computes the volume of a tetrahedron in 3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) TET_XYZ(3,4), the coordinates of the vertices.
!
!    Output, real ( kind = 8 ) VOLUME, the volume of the tetrahedron.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  real ( kind = 8 ) a(4,4)
  real ( kind = 8 ) r8mat_det_4d
  real ( kind = 8 ) tet_xyz(dim_num,4)
  real ( kind = 8 ) volume

  a(1:dim_num,1:4) = tet_xyz( 1:dim_num,1:4)
  a(4,1:4) = 1.0D+00

  volume = abs ( r8mat_det_4d ( a ) ) / 6.0D+00

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
  integer   ( kind = 4 ) d
  integer   ( kind = 4 ) h
  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) mm
  character ( len = 9  ), parameter, dimension(12) :: month = (/ &
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

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
