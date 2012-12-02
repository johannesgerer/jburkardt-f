subroutine boundary_type ( node_num, node_xy, node_boundary, node_type, &
  node_u_condition, node_v_condition, node_p_condition )

!*******************************************************************************
!
!! BOUNDARY_TYPE determines the type of boundary conditions imposed.
!
!  Discussion:
!
!    On input, the calling program has already determined the "type"
!    of every node (vertex or midside), and whether or not it lies
!    on the boundary.
!
!    The program has also set up an initial guess for the boundary
!    conditions, by setting every boundary node to have Dirichlet
!    conditions for U and V, and by setting a single vertex boundary
!    node to have Dirichlet boundary conditions for P.
!
!    The user is free to adjust these boundary conditions in any
!    reasonable way.  
!
!    The most obvious adjustment is to change some velocity boundary
!    conditions to Neumann conditions.  Keep in mind that, for the moment, 
!    we are only supporting zero Neumann conditions.
!
!    However, it is also possible to constrain ANY variable, whether it
!    is on the boundary or not, or to UNCONSTRAIN any variable that
!    has been tentatively constrained.  You simply have to "warn" the code,
!    by setting U_TYPE, V_TYPE or P_TYPE appropriately, and by supplying
!    a value for the right hand side if you are doing a Dirichlet condition.
!
!
!    For the channel flow, we intend for the Dirichlet boundary conditions 
!    to be applied on the inflow and horizontal walls, and zero Neumann
!    conditions on the outflow.
!
!    The calling program has already found the boundary, and guessed that
!    all boundary velocities are constrained by Dirichlet conditions.  So
!    we only have to switch the velocity variables on the outflow to
!    Neumann type.
!
!
!    The pressure is specified to be zero at a single node, but we let
!    the main program take care of that specification.
!
!                           DIRICHLET
!                           U_BC = V_BC = 0
!                           no slip wall
!
!                       1 +---------------+
!                          
!    Inflow -->                           --->  Outflow
!    DIRICHLET                                  NEUMANN
!    U = Y*(1-Y)        0 +---------------+     dU/dn = dV/dn = 0
!    V = 0
!                         0               3
!
!                           no slip wall
!                           DIRICHLET
!                           U_BC = V_BC = 0
!
!  Modified:
!
!    07 October 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NODE_NUM, the number of nodes
!
!    Input, real ( kind = 8 ) NODE_XY(2,NODE_NUM), the coordinates of the nodes.
!
!    Input, logical NODE_BOUNDARY(NODE_NUM), is TRUE if the node is
!    found to lie on the boundary of the region.
!
!    Input, integer NODE_TYPE(NODE_NUM), determines if the node is a 
!    vertex or midside node.
!    1, the node is a vertex (P, U, V variables are associated with it).
!    2, the node is a midside node (only U and V variables are associated.)
!
!    Input/output, integer NODE_U_CONDITION(NODE_NUM), 
!    indicates the condition used to determine horizontal velocity at a node.
!    0, there is no condition (and no variable) at this node.
!    1, a finite element equation is used;
!    2, a Dirichlet condition is used. 
!    3, a Neumann condition is used. 
!
!    Input/output, integer NODE_V_CONDITION(NODE_NUM), 
!    indicates the condition used to determine vertical velocity at a node.
!    0, there is no condition (and no variable) at this node.
!    1, a finite element equation is used;
!    2, a Dirichlet condition is used.
!    3, a Neumann condition is used.  
!
!    Input/output, integer NODE_P_CONDITION(NODE_NUM), 
!    indicates the condition used to determine pressure at a node.
!    0, there is no condition (and no variable) at this node.
!    1, a finite element equation is used;
!    2, a Dirichlet condition is used. 
!    3, a Neumann condition is used. 
!
  implicit none

  integer node_num

  integer node
  logical node_boundary(node_num)
  integer node_p_condition(node_num)
  integer node_type(node_num)
  integer node_u_condition(node_num)
  integer node_v_condition(node_num)
  real ( kind = 8 ) node_xy(2,node_num)
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  do node = 1, node_num

    x = node_xy(1,node)
    y = node_xy(2,node)
!
!  Reset the boundary condition to Neumann type for velocities
!  at nodes on the outflow.  However, leave the velocities at the
!  very top ( Y = 1 ) and bottom ( Y = 0 ) as Dirichlet.
!
    if ( x == 3.0D+00 .and. 0.0D+00 < y .and. y < 1.0D+00 ) then

      node_u_condition(node) = 3
      node_v_condition(node) = 3

    end if

  end do

  return
end
subroutine dirichlet_condition ( node_num, node_xy, u_bc, v_bc, p_bc )

!*******************************************************************************
!
!! DIRICHLET_CONDITION sets the value of any Dirichlet boundary conditions.
!
!  Discussion:
!
!    Dirichlet boundary conditions might be applied to none, some, or
!    all the boundary nodes, and might apply to any combination of
!    U, V, and P.
!
!    This routine will be asked to supply a right hand side for the
!    Dirichlet conditions for U, V and P at EVERY node.  Simply set
!    the value to zero for nodes and variables at which a Dirichlet
!    condition is not being applied.  But set an appropriate value
!    to U_BC, V_BC or P_BC in cases where a Dirichlet condition is
!    being applied for that degree of freedom at that node.
!
!
!    For the channel flow, we intend for the Dirichlet boundary conditions 
!    to be applied on the inflow and horizontal walls, and zero Neumann
!    conditions on the outflow.
!
!    The pressure is specified to be zero at a single node, but we let
!    the main program take care of that specification.
!
!                           DIRICHLET
!                           U_BC = V_BC = 0
!                           no slip wall
!
!                       1 +---------------+
!                          
!    Inflow -->                           --->  Outflow
!    DIRICHLET                                  NEUMANN
!    U = Y*(1-Y)        0 +---------------+     dU/dn = dV/dn = 0
!    V = 0
!                         0               3
!
!                           no slip wall
!                           DIRICHLET
!                           U_BC = V_BC = 0
!
!  Modified:
!
!    07 October 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, real ( kind = 8 ) NODE_XY(2,NODE_NUM), the coordinates of the nodes.
!
!    Output, real ( kind = 8 ) U_BC(NODE_NUM), V_BC(NODE_NUM), P_BC(NODE_NUM), 
!    the values of the boundary conditions on horizontal velocity, 
!    vertical velocity, and pressure.
!
  implicit none

  integer node_num

  integer node
  real ( kind = 8 ) node_xy(2,node_num)
  real ( kind = 8 ) p_bc(node_num)
  real ( kind = 8 ) u_bc(node_num)
  real ( kind = 8 ) v_bc(node_num)
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  do node = 1, node_num

    x = node_xy(1,node)
    y = node_xy(2,node)
!
!  Inflow:
!
    if ( x == 0.0D+00 ) then

      u_bc(node) = y * ( 1.0D+00 - y )
      v_bc(node) = 0.0D+00
      p_bc(node) = 0.0D+00
!
!  Top or bottom walls.
!
    else if ( y == 0.0D+00 .or. y == 1.0D+00 ) then

      u_bc(node) = 0.0D+00
      v_bc(node) = 0.0D+00
      p_bc(node) = 0.0D+00
!
!  No other Dirichlet conditions are imposed.
!  We set these array values to zero.
!
    else

      u_bc(node) = 0.0D+00
      v_bc(node) = 0.0D+00
      p_bc(node) = 0.0D+00

    end if

  end do

  return
end
subroutine neumann_condition ( node_num, node_xy, u_bc, v_bc, p_bc )

!*******************************************************************************
!
!! NEUMANN_CONDITION sets the value of any Neumann boundary conditions.
!
!  Discussion:
!
!    Note that, at the moment, we are simply trying to implement
!    ZERO Neumann boundary conditions, and are not ready to try
!    to implement the NONZERO case.  Therefore, setting nonzero values
!    here is unlikely to work for some time!
!
!
!    Neumann boundary conditions might be applied to none, some, or
!    all the boundary nodes, and might apply to any combination of
!    U, V, and P.
!
!    This routine will be asked to supply a right hand side for the
!    Neumann conditions for U, V and P at EVERY node.  Simply set
!    the value to zero for nodes and variables at which a Neumann
!    condition is not being applied.  But set an appropriate value
!    to U_BC, V_BC or P_BC in cases where a Neumann condition is
!    being applied for that degree of freedom at that node.
!
!
!    For the channel flow, we intend for the Dirichlet boundary conditions 
!    to be applied on the inflow and horizontal walls, and zero Neumann
!    conditions on the outflow.
!
!    The pressure is specified to be zero at a single node, but we let
!    the main program take care of that specification.
!
!                           DIRICHLET
!                           U_BC = V_BC = 0
!                           no slip wall
!
!                       1 +---------------+
!                          
!    Inflow -->                           --->  Outflow
!    DIRICHLET                                  NEUMANN
!    U = Y*(1-Y)        0 +---------------+     dU/dn = dV/dn = 0
!    V = 0
!                         0               3
!
!                           no slip wall
!                           DIRICHLET
!                           U_BC = V_BC = 0
!
!  Modified:
!
!    07 October 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, real ( kind = 8 ) NODE_XY(2,NODE_NUM), the coordinates of the nodes.
!
!    Output, real ( kind = 8 ) U_BC(NODE_NUM), V_BC(NODE_NUM), P_BC(NODE_NUM),
!    the values of the boundary conditions on horizontal velocity, 
!    vertical velocity, and pressure.
!
  implicit none

  integer node_num

  integer node
  real ( kind = 8 ) node_xy(2,node_num)
  real ( kind = 8 ) p_bc(node_num)
  real ( kind = 8 ) u_bc(node_num)
  real ( kind = 8 ) v_bc(node_num)
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  do node = 1, node_num

    x = node_xy(1,node)
    y = node_xy(2,node)
!
!  Outflow.
!
    if ( x == 3.0D+00 ) then

      u_bc(node) = 0.0D+00
      v_bc(node) = 0.0D+00
      p_bc(node) = 0.0D+00
!
!  No other Neumann conditions are imposed.
!  We simply set these array values to zero.
!
    else

      u_bc(node) = 0.0D+00
      v_bc(node) = 0.0D+00
      p_bc(node) = 0.0D+00

    end if

  end do

  return
end
subroutine rhs ( node_num, node_xy, u_rhs, v_rhs, p_rhs )

!*******************************************************************************
!
!! RHS gives the right-hand side of the differential equation.
!
!  Modified:
!
!    07 September 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, real ( kind = 8 ) NODE_XY(2,NODE_NUM), the coordinates of the nodes.
!
!    Output, real ( kind = 8 ) U_RHS(NODE_NUM), V_RHS(NODE_NUM),
!    P_RHS(NODE_NUM), the right hand sides of the differential equations 
!    at the nodes.
!
  implicit none

  integer node_num

  real ( kind = 8 ) node_xy(2,node_num)
  real ( kind = 8 ) p_rhs(node_num)
  real ( kind = 8 ) u_rhs(node_num)
  real ( kind = 8 ) v_rhs(node_num)

  p_rhs(1:node_num) = 0.0D+00
  u_rhs(1:node_num) = 0.0D+00
  v_rhs(1:node_num) = 0.0D+00

  return 
end
