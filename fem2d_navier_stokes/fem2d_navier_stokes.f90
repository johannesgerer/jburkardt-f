program main

!*****************************************************************************80
!
!! MAIN is the main routine of FEM2D_NAVIER_STOKES.
!
!  Discussion:
!
!    This program solves the steady incompressible Navier Stokes equations
!    for velocity vector W and scalar pressure P:
!
!      -nu * Laplacian W(X,Y) + W dot Grad W + Grad P(X,Y) = F(X,Y)
!
!                                               Div W(X,Y) = G(X,Y)
!
!    in an arbitrary triangulated region in the plane.
!
!    Let U and V denote the scalar components of the velocity vector W.
!
!    Along the boundary of the region, the user controls the type of
!    boundary condition to be imposed, if any.  Currently, these
!    conditions may be of Dirichlet form:
!
!      U(X,Y) = U_BC(X,Y)
!      V(X,Y) = V_BC(X,Y)
!      P(X,Y) = P_BC(X,Y)
!
!    or Neumann form with ZERO right hand side:
!
!      dU/dn(X,Y) = 0
!      dV/dn(X,Y) = 0
!      dP/dn(X,Y) = 0
!
!    The code uses the finite element method.  The Taylor-Hood element
!    is used, in which a single reference triangle is used to define
!    both a piecewise quadratic representation of velocity, and a piecewise
!    linear representation of pressure.
!
!  Geometry specification:
!
!    The user defines the geometry by supplying two data files
!    which list the node coordinates, and list the nodes that make up
!    each element.
!
!  Equation specification:
!
!    The user specifies
!
!    * the kinematic viscosity NU,
!
!    * the type of boundary conditions imposed:
!
!      subroutine boundary_type ( node_num, node_xy, node_boundary, node_type,
!        node_u_condition, node_v_condition, node_p_condition )
!
!    * the right hand side of any Dirichlet boundary conditions:
!
!      subroutine dirichlet_condition ( node_num, node_xy, u_bc, v_bc, p_bc )
!
!    * the right hand side of any Neumann boundary conditions
!      (currently, nonzero values will be ignored):
!
!      subroutine neumann_condition ( node_num, node_xy, u_bc, v_bc, p_bc )
!
!    * the right hand side of the Navier Stokes equation:
!
!      subroutine rhs ( node_num, node_xy, u_rhs, v_rhs, p_rhs )
!
!  Usage:
!
!    fem2d_navier_stokes node_file element_file
!
!    invokes the program:
!
!    * "node_file", the coordinates of the nodes;
!    * "element_file", the indices of nodes that make up each element.
!
!    Graphics files created include:
!
!    * "nodes6.eps", an image of the nodes;
!    * "triangles6.eps", an image of the quadratic triangles;
!
!    Data files created include:
!
!    * "nodes3.txt", the nodes associated with pressure;
!    * "triangles3.txt", the linear triangles associated with pressure;
!    * "stokes_pressure3.txt", the Stokes pressure at the pressure nodes;
!    * "stokes_velocity6.txt", the Stokes velocity at the velocity nodes.
!    * "pressure3.txt", the pressure at the pressure nodes;
!    * "velocity6.txt", the velocity at the velocity nodes.
!
!  Modified:
!
!    07 October 2006
!
!  Author:
!
!    John Burkardt
!
!  Command Line Parameters:
!
!    Command line argument, character ( len = * ) NODE_FILE_NAME,
!    the name of the node file.  If this argument is not supplied,
!    it will be requested.
!
!    Command line argument, character ( len = * ) ELEMENT_FILE_NAME,
!    the name of the element file.  If this argument is not supplied,
!    it will be requested.
!
!  Local parameters:
!
!    Local, real ( kind = 8 ) A(3*IB+1,VARIABLE_NUM), the coefficient matrix.
!
!    Local, integer ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM);
!    ELEMENT_NODE(I,J) is the global index of local node I in element J.
!
!    Local, integer ELEMENT_NUM, the number of elements.
!
!    Local, integer ELEMENT_ORDER, the element order.
!
!    Local, real ( kind = 8 ) F(VARIABLE_NUM), the right hand side.
!
!    Local, integer IB, the half-bandwidth of the matrix.
!
!    Local, integer IT_MAX, the maximum number of Newton iterations allowed.
!    0 is a legal value.  It simply solves the Stokes problem, computes
!    the Navier-Stokes residual, and stops.
!    1 does the above, and then takes one Newton step, and so on.
!
!    Local, logical NODE_BOUNDARY(NODE_NUM), is TRUE if the node is
!    found to lie on the boundary of the region.
!
!    Local, real ( kind = 8 ) NODE_C(VARIABLE_NUM), the finite element
!    coefficients of the current solution estimate.
!
!    Local, real ( kind = 8 ) NODE_C_DEL(VARIABLE_NUM), the correction to
!    the finite element coefficients of the current solution estimate.
!
!    Local, real ( kind = 8 ) NODE_C_OLD(VARIABLE_NUM), the finite element
!    coefficients computed by the prefious step of the iteration.
!
!    Local, integer NODE_NUM, the number of nodes.
!
!    Local, integer NODE_P_CONDITION(NODE_NUM),
!    indicates the condition used to determine pressure at a node.
!    0, there is no condition at this node.
!    1, a finite element equation is used;
!    2, a Dirichlet condition is used.
!    3, a Neumann condition is used.
!
!    Local, integer NODE_P_VARIABLE(NODE_NUM),
!    is the index of the pressure variable associated with the node,
!    or -1 if there is no associated pressure variable.
!
!    Local, integer NODE_TYPE(NODE_NUM), determines if the node is a
!    vertex or midside node.
!    1, the node is a vertex (P, U, V variables are associated with it).
!    2, the node is a midside node (only U and V variables are associated.)
!
!    Local, integer NODE_U_CONDITION(NODE_NUM),
!    indicates the condition used to determine horizontal velocity at a node.
!    0, there is no condition at this node.
!    1, a finite element equation is used;
!    2, a Dirichlet condition is used.
!    3, a Neumann condition is used.
!
!    Local, integer NODE_U_VARIABLE(NODE_NUM),
!    is the index of the horizontal velocity variable associated with the node.
!
!    Local, integer NODE_V_CONDITION(NODE_NUM),
!    indicates the condition used to determine vertical velocity at a node.
!    0, there is no condition at this node.
!    1, a finite element equation is used;
!    2, a Dirichlet condition is used.
!    3, a Neumann condition is used.
!
!    Local, integer NODE_V_VARIABLE(NODE_NUM),
!    is the index of the vertical velocity variable associated with the node.
!
!    Local, real ( kind = 8 ) NODE_XY(2,NODE_NUM), the coordinates of nodes.
!
!    Local, integer NODE3_NUM, the number of pressure nodes.
!
!    Local, integer NODE3_LABEL(NODE_NUM), contains the renumbered
!    labels of pressure nodes, and -1 for nodes that are not pressure nodes.
!
!    Local, real ( kind = 8 ) NU, the kinematic viscosity.
!
!    Local, integer QUAD_NUM, the number of quadrature points used for assembly.
!    This is currently set to 3, the lowest reasonable value.  Legal values
!    are 1, 3, 4, 6, 7, 9, 13, and for some problems, a value of QUAD_NUM
!    greater than 3 may be appropriate.
!
!    Local, integer VARIABLE_NUM, the number of variables.
!
  implicit none

  real ( kind = 8 ), allocatable, dimension (:,:) :: a
  logical, parameter :: debugging = .false.
  integer dim_num
  integer element
  character ( len = 255 ) element_file_name
  integer, allocatable, dimension(:,:) :: element_node
  integer element_num
  integer element_order
  real ( kind = 8 ), allocatable, dimension (:) :: f
  character ( len = 80 ) file_name
  integer i
  integer ib
  integer ierr
  integer ip
  integer, parameter :: it_max = 5
  integer it_num
  integer iu
  integer iv
  integer j
  integer job
  integer neumann_num
  integer node
  logical, allocatable, dimension(:) :: node_boundary
  real ( kind = 8 ), allocatable, dimension (:) :: node_c
  real ( kind = 8 ), allocatable, dimension (:) :: node_c_del
  real ( kind = 8 ) node_c_del_norm
  real ( kind = 8 ), allocatable, dimension (:) :: node_c_old
  integer, allocatable, dimension(:) :: node_condition
  character ( len = 255 ) node_file_name
  logical node_label
  integer node_num
  integer, allocatable, dimension(:) :: node_p_condition
  integer, allocatable, dimension(:) :: node_p_variable
  real ( kind = 8 ), allocatable, dimension (:) :: node_r
  real ( kind = 8 ) :: node_r_norm
  integer node_show
  integer, allocatable, dimension(:) :: node_type
  integer, allocatable, dimension(:) :: node_u_condition
  integer, allocatable, dimension(:) :: node_u_variable
  integer, allocatable, dimension(:) :: node_v_condition
  integer, allocatable, dimension(:) :: node_v_variable
  real ( kind = 8 ), allocatable, dimension(:,:) :: node_xy
  integer node3_num
  integer, allocatable, dimension(:) :: node3_label
  real ( kind = 8 ) :: nu = 1.0D+00
  integer p_node
  integer, allocatable, dimension (:) :: pivot
  integer, parameter :: quad_num = 7
  real ( kind = 8 ) :: res_tol = 0.00001D+00
  integer triangle_show
  integer variable_num

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FEM2D_NAVIER_STOKES'
  write ( *, '(a)' ) '  FORTRAN90 version:'
  write ( *, '(a)' ) '  Finite element solution of the'
  write ( *, '(a)' ) '  the steady incompressible Navier Stokes equations'
  write ( *, '(a)' ) '  on a triangulated grid in 2 dimensions.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  -nu * ( Uxx + Uyy ) + UUx + VUy + dPdx = F1(x,y)'
  write ( *, '(a)' ) '  -nu * ( Vxx + Vyy ) + UVx + VVy + dPdy = F2(x,y)'
  write ( *, '(a)' ) '                         Ux +  Vy        = F3(x,y).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Boundary conditions may be of Dirichlet type:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    U(x,y) = U_BC(x,y)'
  write ( *, '(a)' ) '    V(x,y) = V_BC(x,y)'
  write ( *, '(a)' ) '    P(x,y) = P_BC(x,y)'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  or of Neumann type with zero right hand side:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    dU/dn(x,y) = 0'
  write ( *, '(a)' ) '    dV/dn(x,y) = 0'
  write ( *, '(a)' ) '    dP/dn(x,y) = 0'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The finite element method uses Taylor-Hood'
  write ( *, '(a)' ) '  triangular elements which are linear for pressure'
  write ( *, '(a)' ) '  and quadratic for velocity.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) &
    '  Maximum number of Newton iterations IT_MAX = ', it_max
  write ( *, '(a,i8)' ) '  Quadrature order =          ', quad_num
  write ( *, '(a,g14.6)' ) '  The kinematic viscosity NU = ', nu
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Current status:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  * testing zero Neumann condition option.'
!
!  Get the file names.
!
  call file_name_specification ( node_file_name, element_file_name )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Node file is "' // trim ( node_file_name ) // '".'
  write ( *, '(a)' ) '  Element file is "' // trim ( element_file_name ) &
    // '".'
!
!  Read the node coordinate file.
!
  call r8mat_header_read ( node_file_name, dim_num, node_num )

  write ( *, '(a,i8)' ) '  Number of nodes =          ', node_num

  allocate ( node_boundary(node_num) )
  allocate ( node_p_condition(node_num) )
  allocate ( node_p_variable(node_num) )
  allocate ( node_type(node_num) )
  allocate ( node_u_condition(node_num) )
  allocate ( node_u_variable(node_num) )
  allocate ( node_v_condition(node_num) )
  allocate ( node_v_variable(node_num) )
  allocate ( node_xy(dim_num,node_num) )

  call r8mat_data_read ( node_file_name, dim_num, node_num, node_xy )

  call r8mat_transpose_print_some ( dim_num, node_num, node_xy, 1, 1, 2, 10, &
    '  First 10 nodes' )
!
!  Read the element file.
!
  call i4mat_header_read ( element_file_name, element_order, element_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Element order =            ', element_order
  write ( *, '(a,i8)' ) '  Number of elements =       ', element_num

  if ( element_order /= 6 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FEM2D_NAVIER_STOKES - Fatal error!'
    write ( *, '(a,i8)' ) '  The input triangulation has order ', element_order
    write ( *, '(a)' ) '  But a triangulation of order 6 is required.'
    stop
  end if

  allocate ( element_node(element_order,element_num) )

  call i4mat_data_read ( element_file_name, element_order, element_num, &
    element_node )

  call i4mat_transpose_print_some ( element_order, element_num, &
    element_node, 1, 1, element_order, 10, '  First 10 elements' )
!
!  Determine the "type" of each node.
!  A vertex node, of type 1, has U, V, and P variables.
!  A midside node, of type 2, has U and V only.
!
  node_type(1:node_num) = 1

  do element = 1, element_num
    do j = 4, 6
      node = element_node(j,element)
      node_type(node) = 2
    end do
  end do
!
!  Determine which nodes are boundary nodes.
!
  call triangulation_order6_boundary_node ( node_num, element_order, &
    element_num, element_node, node_boundary )

  if ( .false. ) then
    call lvec_print ( node_num, node_boundary, '    Node  Boundary?' )
  end if
!
!  Determine the node conditions.
!  For now, we'll just assume all boundary nodes are Dirichlet.
!
!  All conditions begin as finite element conditions.
!
  node_p_condition(1:node_num) = 1
  node_u_condition(1:node_num) = 1
  node_v_condition(1:node_num) = 1
!
!  Conditions on velocities associated with a boundary node are Dirichlet
!  conditions.
!
  do node = 1, node_num
    if ( node_boundary(node) ) then
      node_u_condition(node) = 2
      node_v_condition(node) = 2
    end if
  end do
!
!  Midside nodes have no associated pressure variable.
!
  do node = 1, node_num
    if ( node_type(node) == 2 ) then
      node_p_condition(node) = 0
    end if
  end do
!
!  Replace a single finite element pressure condition by a Dirichlet
!  condition.
!
  p_node = -1
  do node = 1, node_num
    if ( node_p_condition(node) == 1 ) then
      node_p_condition(node) = 2
      p_node = node
      exit
    end if
  end do

  if ( p_node == -1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FEM2D_NAVIER_STOKES - Fatal error!'
    write ( *, '(a)' ) '  Unable to find a finite element pressure condition'
    write ( *, '(a)' ) '  suitable for replacement by a Dirichlet condition.'
    stop
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Dirichlet boundary condition on pressure'
  write ( *, '(a,i8)' ) '  will be applied at node ', p_node
!
!  Allow the user to examine and modify the tentative boundary conditions.
!
  call boundary_type ( node_num, node_xy, node_boundary, node_type, &
    node_u_condition, node_v_condition, node_p_condition )

  neumann_num = 0

  do node = 1, node_num

    if ( node_u_condition(node) == 3 ) then
      neumann_num = neumann_num + 1
    end if

    if ( node_v_condition(node) == 3 ) then
      neumann_num = neumann_num + 1
    end if

    if ( node_p_condition(node) == 3 ) then
      neumann_num = neumann_num + 1
    end if

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of Neumann conditions added = ', neumann_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Boundary conditions per node:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      Node    U_cond    V_cond    P_cond'
  write ( *, '(a)' ) ' '
  do node = 1, node_num
    write ( *, '(2x,i8,2x,i8,2x,i8,2x,i8)' ) &
      node, node_u_condition(node), node_v_condition(node), &
      node_p_condition(node)
  end do
!
!  Number the variables.
!
  variable_num = 0

  do node = 1, node_num

    variable_num = variable_num + 1
    node_u_variable(node) = variable_num

    variable_num = variable_num + 1
    node_v_variable(node) = variable_num

    if ( node_type(node) == 1 ) then
      variable_num = variable_num + 1
      node_p_variable(node) = variable_num
    else
      node_p_variable(node) = -1
    end if

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Total number of variables is ', variable_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Variable indices per node:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      Node         U         V         P'
  write ( *, '(a)' ) ' '
  do node = 1, node_num
    write ( *, '(2x,i8,2x,i8,2x,i8,2x,i8)' ) &
      node, node_u_variable(node), node_v_variable(node), node_p_variable(node)
  end do
!
!  Determine the bandwidth of the Stokes stiffness matrix
!  and the Navier-Stokes jacobian.
!
  call bandwidth ( element_order, element_num, element_node, &
    node_num, node_p_variable, node_u_variable, node_v_variable, ib )

  write ( *, '(a)'    ) ' '
  write ( *, '(a,i8)' ) '  The matrix half bandwidth is ', ib
  write ( *, '(a,i8)' ) '  The matrix bandwidth is      ', 2 * ib + 1
  write ( *, '(a,i8)' ) '  The storage bandwidth is     ', 3 * ib + 1
!
!  Plot the nodes.
!
  if ( node_num <= 100 ) then

    file_name = 'nodes6.eps'
    node_label = .true.

    call points_plot ( file_name, node_num, node_xy, node_label )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Order 6 nodes plotted in "' &
      // trim ( file_name ) // '".'

  end if
!
!  Plot the triangles.
!
  if ( node_num <= 100 ) then

    file_name = 'triangles6.eps'
    node_show = 2
    triangle_show = 2

    call triangulation_order6_plot ( file_name, node_num, &
      node_xy, element_num, element_node, node_show, triangle_show )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Order 6 triangles plotted in "' &
      // trim ( file_name ) // '".'

  end if
!
!  Allocate space for the coefficient matrix A and right hand side F.
!
  allocate ( a(3*ib+1,variable_num) )
  allocate ( f(variable_num) )
  allocate ( node_c(variable_num) )
  allocate ( node_c_del(variable_num) )
  allocate ( node_c_old(variable_num) )
  allocate ( node_r(variable_num) )
  allocate ( pivot(variable_num) )
!
!  Get an initial condition, by assembling the Stokes coefficient matrix A
!  and the right-hand side F, and solving.
!
  call assemble_stokes ( node_num, element_num, quad_num, &
    variable_num, node_xy, node_p_variable, node_u_variable, &
    node_v_variable, element_node, nu, ib, a, f )
!
!  Print a tiny portion of the matrix.
!
  if ( debugging ) then

    call dgb_print_some ( variable_num, variable_num, ib, ib, a, 1, 1, 10, 10, &
      '  Part of Stokes matrix:' )

    call r8vec_print_some ( variable_num, f, 1, 10, &
      '  Part of Stokes right hand side:' )

  end if
!
!  Adjust the linear system to account for Dirichlet boundary conditions.
!
  call dirichlet_apply ( node_num, node_xy, node_p_variable, &
    node_u_variable, node_v_variable, node_p_condition, &
    node_u_condition, node_v_condition, variable_num, ib, a, f )

  if ( debugging ) then

    call dgb_print_some ( variable_num, variable_num, ib, ib, a, 1, 1, 10, 10, &
      '  Part of Stokes matrix adjusted for Dirichlet BC:' )

    call r8vec_print_some ( variable_num, f, 1, 10, &
      '  Part of Stokes right hand side adjusted for Dirichlet BC:' )

  end if
!
!  Adjust the linear system to account for Neumann boundary conditions.
!
  call neumann_apply ( node_num, node_xy, node_p_variable, &
    node_u_variable, node_v_variable, node_p_condition, &
    node_u_condition, node_v_condition, variable_num, f )

  if ( .false. ) then

    call dgb_print_some ( variable_num, variable_num, ib, ib, a, 1, 1, 20, 20, &
      '  Part of Stokes matrix adjusted for Neumann BC:' )

    call r8vec_print_some ( variable_num, f, 1, 10, &
      '  Part of Stokes right hand side, adjusted for Neumann BC:' )

  end if
!
!  Solve the linear system using a banded solver.
!
  call dgb_fa ( variable_num, ib, ib, a, pivot, ierr )

  if ( ierr /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FEM2D_NAVIER_STOKES - Fatal error!'
    write ( *, '(a)' ) '  DGB_FA returned an error condition.'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The linear system was not factored, and the'
    write ( *, '(a)' ) '  algorithm cannot proceed.'
    stop
  end if

  job = 0
  node_c(1:variable_num) = f(1:variable_num)

  call dgb_sl ( variable_num, ib, ib, a, pivot, node_c, job )

  if ( debugging ) then
    call r8vec_print_some ( variable_num, node_c, 1, 10, &
      '  Part of the solution vector:' )
  end if
!
!  Print the Stokes solution vector based at nodes.
!
  if ( debugging ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Solution to the STOKES Equations:'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '      Node       U               V               P'
    write ( *, '(a)' ) ' '

    do node = 1, node_num

      iu = node_u_variable(node)
      iv = node_v_variable(node)
      ip = node_p_variable(node)

      if ( 0 < ip ) then
        write ( *, '(2x,i8,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
          node, node_c(iu), node_c(iv), node_c(ip)
      else
        write ( *, '(2x,i8,2x,g14.6,2x,g14.6)' ) &
          node, node_c(iu), node_c(iv)
      end if

    end do

  end if
!
!  Compute a renumbering of the pressure nodes.
!
  node3_num = 0

  do node = 1, node_num
    if ( node_type(node) == 1 ) then
      node3_num = node3_num + 1
    end if
  end do

  allocate ( node3_label(1:node_num) )

  node3_num = 0

  do node = 1, node_num
    if ( node_type(node) == 1 ) then
      node3_num = node3_num + 1
      node3_label(node) = node3_num
    else
      node3_label(node) = -1
    end if
  end do
!
!  Write the pressure nodes to a file.
!
  file_name = 'nodes3.txt'

  call nodes3_write ( file_name, node_num, node_xy, node_type )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Pressure nodes written to "' &
    // trim ( file_name ) // '".'
!
!  Write the pressure triangles to a file.
!
  file_name = 'triangles3.txt'

  call triangles3_write ( file_name, element_num, element_node, &
    node_num, node3_label )

  deallocate ( node3_label )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Pressure triangles written to "' &
    // trim ( file_name ) // '".'
!
!  Write the pressures to a file.
!
  if ( .false. ) then

    file_name = 'stokes_pressure3.txt'

    call pressure3_write ( file_name, node_num, node_p_variable, &
      variable_num, node_c )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Stokes pressures written to "' &
      // trim ( file_name ) // '".'

  end if
!
!  Write the velocities to a file.
!
  if ( .false. ) then

    file_name = 'stokes_velocity6.txt'

    call velocity6_write ( file_name, node_num, node_u_variable, &
      node_v_variable, variable_num, node_c )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Stokes velocities written to "' &
      // trim ( file_name ) // '".'

  end if
!
!  Now we have a solution of the Stokes equations.
!  This is used as a starting value for the Newton iteration
!  to be applied to the Navier Stokes equations.
!
  it_num = 0

  do

    call residual_fem ( node_num, node_xy, element_num, &
      element_node, quad_num, node_u_variable, &
      node_v_variable, node_p_variable, variable_num, nu, node_c, node_r )

    if ( debugging ) then
      call r8vec_print_some ( variable_num, node_r, 1, 10, &
        '  Part of Navier-Stokes FEM residual:' )
    end if

    node_r_norm = sqrt ( sum ( node_r(1:variable_num)**2 ) )

    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  l2-norm of FEM residual = ', node_r_norm

    call residual_adjust_dirichlet ( node_num, node_xy, node_p_variable, &
      node_u_variable, node_v_variable, node_p_condition, &
      node_u_condition, node_v_condition, variable_num, node_c, node_r )

    if ( debugging ) then
      call r8vec_print_some ( variable_num, node_r, 1, 10, &
        '  Part of Navier-Stokes FEM residual adjusted for BC:' )
    end if

    call residual_adjust_neumann ( node_num, node_xy, node_p_variable, &
      node_u_variable, node_v_variable, node_p_condition, &
      node_u_condition, node_v_condition, variable_num, node_c, node_r )

    node_r_norm = sqrt ( sum ( node_r(1:variable_num)**2 ) )

    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  l2-norm of adjusted FEM residual = ', node_r_norm

    if ( node_r_norm <= res_tol ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Convergence.'
      exit
    end if

    if ( it_max <= it_num ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The iteration limit has been reached.'
      exit
    end if

    it_num = it_num + 1
!
!  Compute the finite element jacobian matrix.
!
    call jacobian_fem ( node_num, node_xy, element_num, &
      element_node, quad_num, node_u_variable, &
      node_v_variable, node_p_variable, variable_num, nu, node_c, &
      ib, a )
!
!  Adjust the jacobian for boundary conditions.
!
    call jacobian_adjust_dirichlet ( node_num, node_xy, &
      node_p_variable, node_u_variable, node_v_variable, &
      node_p_condition, node_u_condition, node_v_condition, &
      variable_num, ib, a )
!
!  Factor the jacobian.
!
    call dgb_fa ( variable_num, ib, ib, a, pivot, ierr )

    if ( ierr /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FEM2D_NAVIER_STOKES - Fatal error!'
      write ( *, '(a)' ) '  The Jacobian matrix is singular.'
      stop
    end if
!
!  Set up and solve the Newton system J * dX = - F(x)
!
    node_c_del(1:variable_num) = -node_r(1:variable_num)

    call dgb_sl ( variable_num, ib, ib, a, pivot, node_c_del, job )

    node_c_old(1:variable_num) = node_c(1:variable_num)
    node_c(1:variable_num) = node_c(1:variable_num) + node_c_del(1:variable_num)

    if ( debugging ) then
      call r8vec_print_some ( variable_num, node_c_del, 1, 10, &
        '  Part of Newton correction vector:' )
    end if

    node_c_del_norm = sqrt ( sum ( node_c_del(1:variable_num)**2 ) )

    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  l2-norm of Newton correction = ', &
      node_c_del_norm

  end do
!
!  Print the Navier Stokes solution vector based at nodes.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Solution to the NAVIER STOKES Equations:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      Node       U               V               P'
  write ( *, '(a)' ) ' '

  do node = 1, node_num

    iu = node_u_variable(node)
    iv = node_v_variable(node)
    ip = node_p_variable(node)

    if ( 0 < ip ) then
      write ( *, '(2x,i8,2x,g14.6,2x,g14.6,2x,g14.6)' ) &
       node, node_c(iu), node_c(iv), node_c(ip)
    else
      write ( *, '(2x,i8,2x,g14.6,2x,g14.6)' ) &
       node, node_c(iu), node_c(iv)
    end if

  end do
!
!  Write the pressures to a file.
!
  file_name = 'pressure3.txt'

  call pressure3_write ( file_name, node_num, node_p_variable, &
    variable_num, node_c )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Navier Stokes pressures written to "' &
    // trim ( file_name ) // '".'
!
!  Write the velocities to a file.
!
  file_name = 'velocity6.txt'

  call velocity6_write ( file_name, node_num, node_u_variable, &
    node_v_variable, variable_num, node_c )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Navier Stokes velocities written to "' &
    // trim ( file_name ) // '".'
!
!  Deallocate memory.
!
  deallocate ( a )
  deallocate ( element_node )
  deallocate ( f )
  deallocate ( node_boundary )
  deallocate ( node_c )
  deallocate ( node_c_del )
  deallocate ( node_c_old )
  deallocate ( node_p_condition )
  deallocate ( node_p_variable )
  deallocate ( node_r )
  deallocate ( node_u_condition )
  deallocate ( node_u_variable )
  deallocate ( node_v_condition )
  deallocate ( node_v_variable )
  deallocate ( node_xy )
  deallocate ( pivot )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FEM2D_NAVIER_STOKES:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine assemble_stokes ( node_num, element_num, quad_num, &
  variable_num, node_xy, node_p_variable, node_u_variable, &
  node_v_variable, element_node, nu, ib, a, f )

!*****************************************************************************80
!
!! ASSEMBLE_STOKES assembles the finite element Stokes equations.
!
!  Discussion:
!
!    The matrix is known to be banded.  A special matrix storage format
!    is used to reduce the space required.  Details of this format are
!    discussed in the routine DGB_FA.
!
!    The Stokes equations in weak form are:
!
!      Integral ( nu * ( dBdx(I) * dUdx + dBdy(I) * dUdy )
!        + B(I) * ( dPdx - U_RHS ) ) = 0
!
!      Integral ( nu * ( dBdx(I) * dVdx + dBdy(I) * dVdy )
!        + B(I) * (  dPdy - V_RHS ) ) = 0
!
!      Integral ( Q(I) * ( dUdx + dVdy - P_RHS ) ) = 0
!
!    Once the basic finite element system is set up by this routine, another
!    routine adjusts the system to account for boundary conditions.
!
!  Modified:
!
!    16 September 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, integer ELEMENT_NUM, the number of elements.
!
!    Input, integer QUAD_NUM, the number of quadrature points in an element.
!
!    Input, integer VARIABLE_NUM, the number of unknowns.
!
!    Input, real ( kind = 8 ) NODE_XY(2,NODE_NUM), the coordinates
!    of the nodes.
!
!    Input, integer NODE_P_VARIABLE(NODE_NUM), the index of the pressure
!    variable associated with a node, or -1 if there is none.
!
!    Input, integer NODE_U_VARIABLE(NODE_NUM), the index of the horizontal
!    velocity variable associated with a node, or -1 if there is none.
!
!    Input, integer NODE_V_VARIABLE(NODE_NUM), the index of the vertical
!    velocity variable associated with a node, or -1 if there is none.
!
!    Input, integer ELEMENT_NODE(6,ELEMENT_NUM), the nodes that form each
!    element.  Nodes 1, 2, and 3 are the vertices.  Node 4 is between 1
!    and 2, and so on.
!
!    Input, real ( kind = 8 ) NU, the kinematic viscosity.
!
!    Input, integer IB, the matrix half-bandwidth.
!
!    Output, real ( kind = 8 ) A(3*IB+1,VARIABLE_NUM), the VARIABLE_NUM
!    by VARIABLE_NUM coefficient matrix, stored in a compressed format.
!
!    Output, real ( kind = 8 ) F(VARIABLE_NUM), the right hand side.
!
!  Local parameters:
!
!    Local, real ( kind = 8 ) B(6,QUAD_NUM), DBDX(6,QUAD_NUM),
!    DBDY(6,QUAD_NUM), the values of the quadratic basis functions
!    and first derivatives at the quadrature points.
!
!    Local, real ( kind = 8 ) Q(3,QUAD_NUM), DQDX(3,QUAD_NUM),
!    DQDY(3,QUAD_NUM), the values of the linear basis functions
!    and first derivatives at the quadrature points.
!
!    Local, real ( kind = 8 ) QUAD_W(QUAD_NUM), quadrature weights.
!
!    Local, real ( kind = 8 ) QUAD_XY(2,QUAD_NUM), the quadrature points.
!
  implicit none

  integer element_num
  integer ib
  integer node_num
  integer quad_num
  integer variable_num

  real ( kind = 8 ), dimension(3*ib+1,variable_num) :: a
  real ( kind = 8 ) area
  real ( kind = 8 ) b(6,quad_num)
  real ( kind = 8 ), dimension(node_num) :: c
  real ( kind = 8 ) dbdx(6,quad_num)
  real ( kind = 8 ) dbdy(6,quad_num)
  real ( kind = 8 ) dqdx(3,quad_num)
  real ( kind = 8 ) dqdy(3,quad_num)
  integer element
  integer, dimension(6,element_num) :: element_node
  real ( kind = 8 ), dimension(variable_num) :: f
  integer i
  integer ip(3)
  integer iu(6)
  integer iv(6)
  integer j
  integer node_p_variable(node_num)
  integer node_u_variable(node_num)
  integer node_v_variable(node_num)
  real ( kind = 8 ), dimension(2,node_num) :: node_xy
  real ( kind = 8 ) nu
  real ( kind = 8 ) p_rhs(quad_num)
  real ( kind = 8 ) q(3,quad_num)
  integer quad
  real ( kind = 8 ), dimension(quad_num) :: quad_w
  real ( kind = 8 ), dimension(2,quad_num) :: quad_xy
  real ( kind = 8 ), dimension(2,3) :: t3
  real ( kind = 8 ), dimension(2,6) :: t6
  real ( kind = 8 ) triangle_area_2d
  real ( kind = 8 ) u_rhs(quad_num)
  real ( kind = 8 ) v_rhs(quad_num)
  real ( kind = 8 ) w(quad_num)
  real ( kind = 8 ), dimension(2,quad_num) :: xy
!
!  Initialize the arrays to zero.
!
  f(1:variable_num) = 0.0D+00
  a(1:3*ib+1,1:variable_num) = 0.0D+00
!
!  Get the quadrature weights and nodes.
!
  call quad_rule ( quad_num, quad_w, quad_xy )
!
!  Add up all quantities associated with the ELEMENT-th element.
!
  do element = 1, element_num
!
!  Extract the nodes of the linear and quadratic triangles.
!
    t3(1:2,1:3) = node_xy(1:2,element_node(1:3,element))
    t6(1:2,1:6) = node_xy(1:2,element_node(1:6,element))
!
!  Map the quadrature points QUAD_XY to points XY in the physical triangle.
!
    call reference_to_physical_t6 ( t6, quad_num, quad_xy, xy )
    area = abs ( triangle_area_2d ( t3 ) )
    w(1:quad_num) = area * quad_w(1:quad_num)
    call rhs ( quad_num, xy, u_rhs, v_rhs, p_rhs )
!
!  Evaluate the basis functions at the quadrature points.
!
    call basis_mn_t6 ( t6, quad_num, xy, b, dbdx, dbdy )

    call basis_mn_t3 ( t3, quad_num, xy, q, dqdx, dqdy )
!
!  Extract the indices of the finite element coefficients for this element.
!
    iu(1:6) = node_u_variable(element_node(1:6,element))
    iv(1:6) = node_v_variable(element_node(1:6,element))
    ip(1:3) = node_p_variable(element_node(1:3,element))
!
!  The horizontal momentum equation.
!
    do i = 1, 6

      f(iu(i)) = f(iu(i)) + sum &
      ( w(1:quad_num) * u_rhs(1:quad_num) * b(i,1:quad_num) )

      do j = 1, 6

        a(iu(i)-iu(j)+2*ib+1,iu(j)) = a(iu(i)-iu(j)+2*ib+1,iu(j)) + sum &
        ( w(1:quad_num) *                                               &
          (                                                             &
            nu * ( dbdx(j,1:quad_num) * dbdx(i,1:quad_num)              &
                 + dbdy(j,1:quad_num) * dbdy(i,1:quad_num) )            &
          )                                                             &
        )

      end do

      do j = 1, 3
        a(iu(i)-ip(j)+2*ib+1,ip(j)) = a(iu(i)-ip(j)+2*ib+1,ip(j)) + sum &
          ( w(1:quad_num) * dqdx(j,1:quad_num) * b(i,1:quad_num) )
      end do

    end do
!
!  The vertical momentum equation.
!
    do i = 1, 6

      f(iv(i)) = f(iv(i)) + sum &
      ( w(1:quad_num) * b(i,1:quad_num) * v_rhs(1:quad_num) )

      do j = 1, 6

        a(iv(i)-iv(j)+2*ib+1,iv(j)) = a(iv(i)-iv(j)+2*ib+1,iv(j)) + sum &
        ( w(1:quad_num) *                                               &
          (                                                             &
            nu * ( dbdx(j,1:quad_num) * dbdx(i,1:quad_num)              &
                 + dbdy(j,1:quad_num) * dbdy(i,1:quad_num) )            &
          )                                                             &
        )

      end do

      do j = 1, 3
        a(iv(i)-ip(j)+2*ib+1,ip(j)) = a(iv(i)-ip(j)+2*ib+1,ip(j)) + sum &
        ( w(1:quad_num) * dqdy(j,1:quad_num) * b(i,1:quad_num) )
      end do

    end do
!
!  The pressure equation.
!
    do i = 1, 3

      f(ip(i)) = f(ip(i)) + sum &
      ( w(1:quad_num) * q(i,1:quad_num) * p_rhs(1:quad_num) )

      do j = 1, 6

        a(ip(i)-iu(j)+2*ib+1,iu(j)) = a(ip(i)-iu(j)+2*ib+1,iu(j)) + sum &
        ( w(1:quad_num) * dbdx(j,1:quad_num) * q(i,1:quad_num) )

        a(ip(i)-iv(j)+2*ib+1,iv(j)) = a(ip(i)-iv(j)+2*ib+1,iv(j)) + sum &
        ( w(1:quad_num) * dbdy(j,1:quad_num) * q(i,1:quad_num) )

      end do

    end do

  end do

  return
end
subroutine bandwidth ( element_order, element_num, element_node, &
  node_num, node_p_variable, node_u_variable, node_v_variable, ib )

!*****************************************************************************80
!
!! BANDWIDTH determines the bandwidth of the coefficient matrix.
!
!  Discussion:
!
!    We take the bandwidth to be the maximum difference between the
!    indices of two variables associated with nodes that share an element.
!
!    Therefore, we can compute the bandwidth by examining each element,
!    and finding the maximum difference in indices of any two variables
!    associated with nodes in that element.
!
!  Modified:
!
!    10 June 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ELEMENT_ORDER, the number of nodes per element.
!
!    Input, integer ELEMENT_NUM, the number of elements.
!
!    Input, integer ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM);
!    ELEMENT_NODE(I,J) is the global index of local node I in element J.
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, integer NODE_P_VARIABLE(NODE_NUM),
!    is the index of the pressure variable associated with the node,
!    or -1 if there is no associated pressure variable.
!
!    Input, integer NODE_U_VARIABLE(NODE_NUM),
!    is the index of the horizontal velocity variable associated with the node.
!
!    Input, integer NODE_V_VARIABLE(NODE_NUM),
!    is the index of the vertical velocity variable associated with the node.
!
!    Output, integer IB, the half bandwidth of the matrix.
!
  implicit none

  integer node_num
  integer element_num
  integer element_order

  integer element
  integer element_node(element_order,element_num)
  integer i4_huge
  integer local
  integer ib
  integer node
  integer node_p_variable(node_num)
  integer node_u_variable(node_num)
  integer node_v_variable(node_num)
  integer v
  integer v_max
  integer v_min

  ib = 0

  do element = 1, element_num

    v_max = -i4_huge ( )
    v_min = i4_huge ( )

    do local = 1, element_order

      node = element_node(local,element)

      v = node_u_variable(node)
      v_max = max ( v_max, v )
      v_min = min ( v_min, v )

      v = node_v_variable(node)
      v_max = max ( v_max, v )
      v_min = min ( v_min, v )

      if ( 0 < node_p_variable(node) ) then
        v = node_p_variable(node)
        v_max = max ( v_max, v )
        v_min = min ( v_min, v )
      end if

    end do

    ib = max ( ib, v_max - v_min )

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
!    Input, integer N, the number of evaluation points.
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

  integer n

  real ( kind = 8 ) area
  real ( kind = 8 ) dphidx(3,n)
  real ( kind = 8 ) dphidy(3,n)
  real ( kind = 8 ) p(2,n)
  real ( kind = 8 ) phi(3,n)
  real ( kind = 8 ) t(2,3)

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
!    Input, integer N, the number of evaluation points.
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

  integer n

  real ( kind = 8 ) dphidx(6,n)
  real ( kind = 8 ) dphidy(6,n)
  real ( kind = 8 ) gn(n)
  real ( kind = 8 ) gx(n)
  real ( kind = 8 ) hn(n)
  real ( kind = 8 ) hx(n)
  real ( kind = 8 ) p(2,n)
  real ( kind = 8 ) phi(6,n)
  real ( kind = 8 ) t(2,6)
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
subroutine ch_cap ( c )

!*****************************************************************************80
!
!! CH_CAP capitalizes a single character.
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
  integer itemp

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
!  Examples:
!
!    CH_EQI ( 'A', 'a' ) is .TRUE.
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
!    Output, integer DIGIT, the corresponding integer value.  If C was
!    'illegal', then DIGIT is -1.
!
  implicit none

  character c
  integer digit

  if ( lge ( c, '0' ) .and. lle ( c, '9' ) ) then

    digit = ichar ( c ) - 48

  else if ( c == ' ' ) then

    digit = 0

  else

    digit = -1

  end if

  return
end
subroutine dgb_fa ( n, ml, mu, a, pivot, info )

!*****************************************************************************80
!
!! DGB_FA performs a LINPACK-style PLU factorization of an DGB matrix.
!
!  Discussion:
!
!    The DGB storage format is for an M by N banded matrix, with lower
!    bandwidth ML and upper bandwidth MU.  Storage includes room for ML
!    extra superdiagonals, which may be required to store nonzero entries
!    generated during Gaussian elimination.
!
!    The original M by N matrix is "collapsed" downward, so that diagonals
!    become rows of the storage array, while columns are preserved.  The
!    collapsed array is logically 2*ML+MU+1 by N.
!
!    The following program segment will set up the input.
!
!      m = ml + mu + 1
!      do j = 1, n
!        i1 = max ( 1, j-mu )
!        i2 = min ( n, j+ml )
!        do i = i1, i2
!          k = i - j + m
!          a(k,j) = afull(i,j)
!        end do
!      end do
!
!    This uses rows ML+1 through 2*ML+MU+1 of the array A.
!    In addition, the first ML rows in the array are used for
!    elements generated during the triangularization.
!
!    The ML+MU by ML+MU upper left triangle and the
!    ML by ML lower right triangle are not referenced.
!
!  Modified:
!
!    04 March 1999
!
!  Reference:
!
!    Dongarra, Bunch, Moler, Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N-1.
!
!    Input/output, real ( kind = 8 ) A(2*ML+MU+1,N), on input, the matrix
!    in band storage, on output, information about the LU factorization.
!
!    Output, integer PIVOT(N), the pivot vector.
!
!    Output, integer INFO, singularity flag.
!    0, no singularity detected.
!    nonzero, the factorization failed on the INFO-th step.
!
  implicit none

  integer ml
  integer mu
  integer n

  real ( kind = 8 ) a(2*ml+mu+1,n)
  integer i
  integer i0
  integer info
  integer pivot(n)
  integer j
  integer j0
  integer j1
  integer ju
  integer jz
  integer k
  integer l
  integer lm
  integer m
  integer mm
  real ( kind = 8 ) t
  real ( kind = 8 ) temp

  m = ml + mu + 1
  info = 0
!
!  Zero out the initial fill-in columns.
!
  j0 = mu + 2
  j1 = min ( n, m ) - 1

  do jz = j0, j1
    i0 = m + 1 - jz
    a(i0:ml,jz) = 0.0D+00
  end do

  jz = j1
  ju = 0

  do k = 1, n-1
!
!  Zero out the next fill-in column.
!
    jz = jz + 1
    if ( jz <= n ) then
      a(1:ml,jz) = 0.0D+00
    end if
!
!  Find L = pivot index.
!
    lm = min ( ml, n-k )

    l = m
    do j = m+1, m+lm
      if ( abs ( a(l,k) ) < abs ( a(j,k) ) ) then
        l = j
      end if
    end do

    pivot(k) = l + k - m
!
!  Zero pivot implies this column already triangularized.
!
    if ( a(l,k) == 0.0D+00 ) then
      info = k
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DGB_FA - Fatal error!'
      write ( *, '(a,i8)' ) '  Zero pivot on step ', info
      return
    end if
!
!  Interchange if necessary.
!
    temp = a(l,k)
    a(l,k) = a(m,k)
    a(m,k) = temp
!
!  Compute multipliers.
!
    a(m+1:m+lm,k) = - a(m+1:m+lm,k) / a(m,k)
!
!  Row elimination with column indexing.
!
    ju = max ( ju, mu+pivot(k) )
    ju = min ( ju, n )
    mm = m

    do j = k+1, ju

      l = l - 1
      mm = mm - 1

      if ( l /= mm ) then
        temp = a(l,j)
        a(l,j) = a(mm,j)
        a(mm,j) = temp
      end if

      a(mm+1:mm+lm,j) = a(mm+1:mm+lm,j) + a(mm,j) * a(m+1:m+lm,k)

    end do

  end do

  pivot(n) = n
  if ( a(m,n) == 0.0D+00 ) then
    info = n
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DGB_FA - Fatal error!'
    write ( *, '(a,i8)' ) '  Zero pivot on step ', info
  end if

  return
end
subroutine dgb_print_some ( m, n, ml, mu, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! DGB_PRINT_SOME prints some of a DGB matrix.
!
!  Discussion:
!
!    The DGB storage format is for an M by N banded matrix, with lower
!    bandwidth ML and upper bandwidth MU.  Storage includes room for ML
!    extra superdiagonals, which may be required to store nonzero entries
!    generated during Gaussian elimination.
!
!    The original M by N matrix is "collapsed" downward, so that diagonals
!    become rows of the storage array, while columns are preserved.  The
!    collapsed array is logically 2*ML+MU+1 by N.
!
!    DGB storage is used by LINPACK and LAPACK.
!
!  Modified:
!
!    04 April 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the number of rows of the matrix.
!    M must be positive.
!
!    Input, integer N, the number of columns of the matrix.
!    N must be positive.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than min(M,N)-1..
!
!    Input, real ( kind = 8 ) A(2*ML+MU+1,N), the DGB matrix.
!
!    Input, integer ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed.
!
  implicit none

  integer, parameter :: incx = 5
  integer ml
  integer mu
  integer n

  real ( kind = 8 ) a(2*ml+mu+1,n)
  character ( len = 14 ) ctemp(incx)
  integer i
  integer i2hi
  integer i2lo
  integer ihi
  integer ilo
  integer inc
  integer j
  integer j2
  integer j2hi
  integer j2lo
  integer jhi
  integer jlo
  integer m
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if
!
!  Print the columns of the matrix, in strips of 5.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i7,7x)' ) j
    end do

    write ( *, '(a,5a14)' ) '  Col: ', ( ctemp(j2), j2 = 1, inc )
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )
    i2lo = max ( i2lo, j2lo - mu - ml )
    i2hi = min ( ihi, m )
    i2hi = min ( i2hi, j2hi + ml )

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( i < j - ml - mu  .or. j + ml < i ) then
          ctemp(j2) = '              '
        else
          write ( ctemp(j2), '(g14.6)' ) a(i-j+ml+mu+1,j)
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  write ( *, '(a)' ) ' '

  return
end
subroutine dgb_sl ( n, ml, mu, a, pivot, b, job )

!*****************************************************************************80
!
!! DGB_SL solves a system factored by DGB_FA.
!
!  Discussion:
!
!    The DGB storage format is for an M by N banded matrix, with lower
!    bandwidth ML and upper bandwidth MU.  Storage includes room for ML
!    extra superdiagonals, which may be required to store nonzero entries
!    generated during Gaussian elimination.
!
!    The original M by N matrix is "collapsed" downward, so that diagonals
!    become rows of the storage array, while columns are preserved.  The
!    collapsed array is logically 2*ML+MU+1 by N.
!
!  Modified:
!
!    04 March 1999
!
!  Reference:
!
!    Dongarra, Bunch, Moler, Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input, integer ML, MU, the lower and upper bandwidths.
!    ML and MU must be nonnegative, and no greater than N-1.
!
!    Input, real ( kind = 8 ) A(2*ML+MU+1,N), the LU factors from DGB_FA.
!
!    Input, integer PIVOT(N), the pivot vector from DGB_FA.
!
!    Input/output, real ( kind = 8 )l B(N).
!    On input, the right hand side vector.
!    On output, the solution.
!
!    Input, integer JOB.
!    0, solve A * x = b.
!    nonzero, solve A' * x = b.
!
  implicit none

  integer ml
  integer mu
  integer n

  real ( kind = 8 ) a(2*ml+mu+1,n)
  real ( kind = 8 ) b(n)
  integer pivot(n)
  integer j
  integer job
  integer k
  integer l
  integer la
  integer lb
  integer lm
  integer m
  real ( kind = 8 ) t

  m = mu + ml + 1
!
!  Solve A * x = b.
!
  if ( job == 0 ) then
!
!  Solve L * Y = B.
!
    if ( 1 <= ml ) then

      do k = 1, n-1

        lm = min ( ml, n-k )
        l = pivot(k)

        if ( l /= k ) then
          t    = b(l)
          b(l) = b(k)
          b(k) = t
        end if

        b(k+1:k+lm) = b(k+1:k+lm) + b(k) * a(m+1:m+lm,k)

      end do
    end if
!
!  Solve U * X = Y.
!
    do k = n, 1, -1

      b(k) = b(k) / a(m,k)
      lm = min ( k, m ) - 1
      la = m - lm
      lb = k - lm

      b(lb:lb+lm-1) = b(lb:lb+lm-1) - b(k) * a(la:la+lm-1,k)

    end do
!
!  Solve A' * X = B.
!
  else
!
!  Solve U' * Y = B.
!
    do k = 1, n
      lm = min ( k, m ) - 1
      la = m - lm
      lb = k - lm
      b(k) = ( b(k) - sum ( a(la:la+lm-1,k) * b(lb:lb+lm-1) ) ) &
        / a(m,k)
    end do
!
!  Solve L' * X = Y.
!
    if ( 1 <= ml ) then

      do k = n-1, 1, -1

        lm = min ( ml, n-k )
        b(k) = b(k) + sum ( a(m+1:m+lm,k) * b(k+1:k+lm) )
        l = pivot(k)

        if ( l /= k ) then
          t    = b(l)
          b(l) = b(k)
          b(k) = t
        end if

      end do

    end if

  end if

  return
end
subroutine dirichlet_apply ( node_num, node_xy, node_p_variable, &
  node_u_variable, node_v_variable, node_p_condition, &
  node_u_condition, node_v_condition, variable_num, ib, a, f )

!*****************************************************************************80
!
!! DIRICHLET_APPLY accounts for Dirichlet boundary conditions.
!
!  Modified:
!
!    21 September 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, real ( kind = 8 ) NODE_XY(2,NODE_NUM), the coordinates of nodes.
!
!    Input, integer NODE_P_VARIABLE(NODE_NUM),
!    is the index of the pressure variable associated with the node,
!    or -1 if there is no associated pressure variable.
!
!    Input, integer NODE_U_VARIABLE(NODE_NUM),
!    is the index of the horizontal velocity variable associated with the node.
!
!    Input, integer NODE_V_VARIABLE(NODE_NUM),
!    is the index of the vertical velocity variable associated with the node.
!
!    Input, integer NODE_P_CONDITION(NODE_NUM),
!    indicates the condition used to determine pressure at a node.
!    0, there is no condition at this node.
!    1, a finite element equation is used;
!    2, a Dirichlet condition is used.
!    3, a Neumann condition is used.
!
!    Input, integer NODE_U_CONDITION(NODE_NUM),
!    indicates the condition used to determine horizontal velocity at a node.
!    0, there is no condition at this node.
!    1, a finite element equation is used;
!    2, a Dirichlet condition is used.
!    3, a Neumann condition is used.
!
!    Input, integer NODE_V_CONDITION(NODE_NUM),
!    indicates the condition used to determine vertical velocity at a node.
!    0, there is no condition at this node.
!    1, a finite element equation is used;
!    2, a Dirichlet condition is used.
!    3, a Neumann condition is used.
!
!    Input, integer VARIABLE_NUM, the number of variables.
!
!    Input, integer IB, the half-bandwidth of the matrix.
!
!    Input/output, real ( kind = 8 ) A(3*IB+1,VARIABLE_NUM), the
!    VARIABLE_NUM by VARIABLE_NUM coefficient matrix, stored in a
!    compressed format; on output, the matrix has been adjusted
!    for Dirichlet boundary conditions.
!
!    Input/output, real ( kind = 8 ) F(VARIABLE_NUM), the right hand side.
!    On output, the right hand side has been adjusted for Dirichlet
!    boundary conditions.
!
  implicit none

  integer ib
  integer node_num
  integer variable_num

  real ( kind = 8 ), dimension(3*ib+1,variable_num) :: a
  integer column
  integer column_high
  integer column_low
  integer, parameter :: DIRICHLET = 2
  real ( kind = 8 ), dimension(variable_num) :: f
  integer ip
  integer iu
  integer iv
  integer node
  integer node_p_condition(node_num)
  integer node_p_variable(node_num)
  integer node_u_condition(node_num)
  integer node_u_variable(node_num)
  integer node_v_condition(node_num)
  integer node_v_variable(node_num)
  real ( kind = 8 ), dimension(2,node_num) :: node_xy
  real ( kind = 8 ) p_bc(node_num)
  real ( kind = 8 ) u_bc(node_num)
  real ( kind = 8 ) v_bc(node_num)
  real ( kind = 8 ) value

  call dirichlet_condition ( node_num, node_xy, u_bc, v_bc, p_bc )

  do node = 1, node_num

    iu = node_u_variable(node)
    iv = node_v_variable(node)
    ip = node_p_variable(node)

    if ( node_u_condition(node) == DIRICHLET ) then

      column_low = max ( iu - ib, 1 )
      column_high = min ( iu + ib, variable_num )

      do column = column_low, column_high
        a(iu-column+2*ib+1,column) = 0.0D+00
      end do
      a(2*ib+1,iu) = 1.0D+00

      f(iu) = u_bc(node)

    end if

    if ( node_v_condition(node) == DIRICHLET ) then

      column_low = max ( iv - ib, 1 )
      column_high = min ( iv + ib, variable_num )

      do column = column_low, column_high
        a(iv-column+2*ib+1,column) = 0.0D+00
      end do
      a(2*ib+1,iv) = 1.0D+00

      f(iv) = v_bc(node)


    end if

    if ( 0 < ip ) then

      if ( node_p_condition(node) == DIRICHLET ) then

        column_low = max ( ip - ib, 1 )
        column_high = min ( ip + ib, variable_num )

        do column = column_low, column_high
          a(ip-column+2*ib+1,column) = 0.0D+00
        end do
        a(2*ib+1,ip) = 1.0D+00

        f(ip) = p_bc(node)

      end if

    end if

  end do

  return
end
subroutine file_column_count ( input_file_name, column_num )

!*****************************************************************************80
!
!! FILE_COLUMN_COUNT counts the number of columns in the first line of a file.
!
!  Discussion:
!
!    The file is assumed to be a simple text file.
!
!    Most lines of the file is presumed to consist of COLUMN_NUM words,
!    separated by spaces.  There may also be some blank lines, and some
!    comment lines,
!    which have a "#" in column 1.
!
!    The routine tries to find the first non-comment non-blank line and
!    counts the number of words in that line.
!
!    If all lines are blanks or comments, it goes back and tries to analyze
!    a comment line.
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
!    Input, character ( len = * ) INPUT_FILE_NAME, the name of the file.
!
!    Output, integer COLUMN_NUM, the number of columns in the file.
!
  implicit none

  integer column_num
  logical got_one
  character ( len = * ) input_file_name
  integer input_unit
  integer ios
  character ( len = 256 ) line
!
!  Open the file.
!
  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_file_name, status = 'old', &
    form = 'formatted', access = 'sequential', iostat = ios )

  if ( ios /= 0 ) then
    column_num = -1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_COLUMN_COUNT - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the input file "' &
      // trim ( input_file_name ) // '" on unit ', input_unit
    return
  end if
!
!  Read one line, but skip blank lines and comment lines.
!
  got_one = .false.

  do

    read ( input_unit, '(a)', iostat = ios ) line

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

    rewind ( input_unit )

    do

      read ( input_unit, '(a)', iostat = ios ) line

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

  close ( unit = input_unit )

  if ( .not. got_one ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_COLUMN_COUNT - Warning!'
    write ( *, '(a)' ) '  The file does not seem to contain any data.'
    column_num = -1
    return
  end if

  call s_word_count ( line, column_num )

  return
end
subroutine file_name_specification ( node_file_name, element_file_name )

!*****************************************************************************80
!
!! FILE_NAME_SPECIFICATION determines the names of the input files.
!
!  Modified:
!
!    04 October 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) NODE_FILE_NAME, the name of the node file.
!
!    Output, character ( len = * ) ELEMENT_FILE_NAME, the name
!    of the element file.
!
  implicit none

  integer arg_num
  character ( len = * ) :: element_file_name
  integer iarg
  integer iargc
  character ( len = * ) :: node_file_name
!
!  Get the number of command line arguments.
!
  arg_num = iargc ( )
!
!  If at least one command line argument, it's the node file name.
!
  if ( 1 <= arg_num ) then

    iarg = 1
    call getarg ( iarg, node_file_name )

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_NAME_SPECIFICATION:'
    write ( *, '(a)' ) '  Please enter the name of the node file.'

    read ( *, '(a)' ) node_file_name

  end if
!
!  If at least two command line arguments, the second is the element file.
!
  if ( 2 <= arg_num ) then

    iarg = 2
    call getarg ( iarg, element_file_name )

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_NAME_SPECIFICATION:'
    write ( *, '(a)' ) '  Please enter the name of the element file.'

    read ( *, '(a)' ) element_file_name

  end if

  return
end
subroutine file_row_count ( input_file_name, row_num )

!*****************************************************************************80
!
!! FILE_ROW_COUNT counts the number of row records in a file.
!
!  Discussion:
!
!    It does not count lines that are blank, or that begin with a
!    comment symbol '#'.
!
!  Modified:
!
!    06 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILE_NAME, the name of the input file.
!
!    Output, integer ROW_NUM, the number of rows found.
!
  implicit none

  integer bad_num
  integer comment_num
  integer ierror
  character ( len = * ) input_file_name
  integer input_unit
  integer ios
  character ( len = 100 ) line
  integer record_num
  integer row_num

  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_file_name, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    row_num = -1;
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_ROW_COUNT - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the input file "' // &
      trim ( input_file_name ) // '" on unit ', input_unit
    stop
  end if

  comment_num = 0
  row_num = 0
  record_num = 0
  bad_num = 0

  do

    read ( input_unit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      ierror = record_num
      exit
    end if

    record_num = record_num + 1

    if ( line(1:1) == '#' ) then
      comment_num = comment_num + 1
      cycle
    end if

    if ( len_trim ( line ) == 0 ) then
      comment_num = comment_num + 1
      cycle
    end if

    row_num = row_num + 1

  end do

  close ( unit = input_unit )

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
!    Output, integer IUNIT, the free unit number.
!
  implicit none

  integer i
  integer ios
  integer iunit
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
function i4_huge ( )

!*****************************************************************************80
!
!! I4_HUGE returns a "huge" I4.
!
!  Modified:
!
!    17 April 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer I4_HUGE, a "huge" integer.
!
  implicit none

  integer i4_huge

  i4_huge = huge ( 1 )

  return
end
subroutine i4col_compare ( m, n, a, i, j, isgn )

!*****************************************************************************80
!
!! I4COL_COMPARE compares columns I and J of an I4COL.
!
!  Example:
!
!    Input:
!
!      M = 3, N = 4, I = 2, J = 4
!
!      A = (
!        1  2  3  4
!        5  6  7  8
!        9 10 11 12 )
!
!    Output:
!
!      ISGN = -1
!
!  Modified:
!
!    30 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, integer A(M,N), an array of N columns of vectors of length M.
!
!    Input, integer I, J, the columns to be compared.
!    I and J must be between 1 and N.
!
!    Output, integer ISGN, the results of the comparison:
!    -1, column I < column J,
!     0, column I = column J,
!    +1, column J < column I.
!
  implicit none

  integer m
  integer n

  integer a(m,n)
  integer i
  integer isgn
  integer j
  integer k
!
!  Check.
!
  if ( i < 1 .or. n < i ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4COL_COMPARE - Fatal error!'
    write ( *, '(a)' ) '  Column index I is out of bounds.'
    write ( *, '(a,i8)' ) '  I = ', i
    stop
  end if

  if ( j < 1 .or. n < j ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4COL_COMPARE - Fatal error!'
    write ( *, '(a)' ) '  Column index J is out of bounds.'
    write ( *, '(a,i8)' ) '  J = ', j
    stop
  end if

  isgn = 0

  if ( i == j ) then
    return
  end if

  k = 1

  do while ( k <= m )

    if ( a(k,i) < a(k,j) ) then
      isgn = -1
      return
    else if ( a(k,j) < a(k,i) ) then
      isgn = +1
      return
    end if

    k = k + 1

  end do

  return
end
subroutine i4col_sort_a ( m, n, a )

!*****************************************************************************80
!
!! I4COL_SORT_A ascending sorts an I4COL.
!
!  Discussion:
!
!    In lexicographic order, the statement "X < Y", applied to two real
!    vectors X and Y of length M, means that there is some index I, with
!    1 <= I <= M, with the property that
!
!      X(J) = Y(J) for J < I,
!    and
!      X(I) < Y(I).
!
!    In other words, the first time they differ, X is smaller.
!
!  Modified:
!
!    25 September 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the number of rows of A, and the length of
!    a vector of data.
!
!    Input, integer N, the number of columns of A.
!
!    Input/output, integer A(M,N).
!    On input, the array of N columns of M-vectors.
!    On output, the columns of A have been sorted in ascending
!    lexicographic order.
!
  implicit none

  integer m
  integer n

  integer a(m,n)
  integer i
  integer indx
  integer isgn
  integer j

  if ( m <= 0 ) then
    return
  end if

  if ( n <= 1 ) then
    return
  end if
!
!  Initialize.
!
  i = 0
  indx = 0
  isgn = 0
  j = 0
!
!  Call the external heap sorter.
!
  do

    call sort_heap_external ( n, indx, i, j, isgn )
!
!  Interchange the I and J objects.
!
    if ( 0 < indx ) then

      call i4col_swap ( m, n, a, i, j )
!
!  Compare the I and J objects.
!
    else if ( indx < 0 ) then

      call i4col_compare ( m, n, a, i, j, isgn )

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
subroutine i4col_swap ( m, n, a, i, j )

!*****************************************************************************80
!
!! I4COL_SWAP swaps columns I and J of an I4COL.
!
!  Example:
!
!    Input:
!
!      M = 3, N = 4, I = 2, J = 4
!
!      A = (
!        1  2  3  4
!        5  6  7  8
!        9 10 11 12 )
!
!    Output:
!
!      A = (
!        1  4  3  2
!        5  8  7  6
!        9 12 11 10 )
!
!  Modified:
!
!    04 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns in the array.
!
!    Input/output, integer A(M,N), an array of N columns of length M.
!
!    Input, integer I, J, the columns to be swapped.
!
  implicit none

  integer m
  integer n

  integer a(m,n)
  integer col(m)
  integer i
  integer j

  if ( i < 1 .or. n < i .or. j < 1 .or. n < j ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4COL_SWAP - Fatal error!'
    write ( *, '(a)' ) '  I or J is out of bounds.'
    write ( *, '(a,i8)' ) '  I =    ', i
    write ( *, '(a,i8)' ) '  J =    ', j
    write ( *, '(a,i8)' ) '  N =    ', n
    stop

  end if

  if ( i == j ) then
    return
  end if

  col(1:m) = a(1:m,i)
  a(1:m,i) = a(1:m,j)
  a(1:m,j) = col(1:m)

  return
end
subroutine i4mat_data_read ( input_filename, m, n, table )

!*****************************************************************************80
!
!! I4MAT_DATA_READ reads data from an I4MAT file.
!
!  Discussion:
!
!    An I4MAT is an array of I4's.
!
!    The file may contain more than N points, but this routine
!    will return after reading N points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILENAME, the name of the input file.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Output, integer ( kind = 4 ) TABLE(M,N), the data.
!
  implicit none

  integer   ( kind = 4 )   m
  integer   ( kind = 4 )   n

  integer   ( kind = 4 )   ierror
  character ( len = * )    input_filename
  integer   ( kind = 4 )   input_status
  integer   ( kind = 4 )   input_unit
  integer   ( kind = 4 )   j
  character ( len = 255 )  line
  integer   ( kind = 4 )   table(m,n)
  integer   ( kind = 4 )   x(m)

  ierror = 0

  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_filename, status = 'old', &
    iostat = input_status )

  if ( input_status /= 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4MAT_DATA_READ - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the input file "' // &
      trim ( input_filename ) // '" on unit ', input_unit
    stop
  end if

  j = 0

  do while ( j < n )

    read ( input_unit, '(a)', iostat = input_status ) line

    if ( input_status /= 0 ) then
      ierror = 2
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'I4MAT_DATA_READ - Fatal error!'
      write ( *, '(a)' ) '  Error while reading lines of data.'
      write ( *, '(a,i8)' ) '  Number of values expected per line M = ', m
      write ( *, '(a,i8)' ) '  Number of data lines read, J =         ', j
      write ( *, '(a,i8)' ) '  Number of data lines needed, N =       ', n
      stop
    end if

    if ( line(1:1) == '#' .or. len_trim ( line ) == 0 ) then
      cycle
    end if

    call s_to_i4vec ( line, m, x, ierror )

    if ( ierror /= 0 ) then
      cycle
    end if

    j = j + 1

    table(1:m,j) = x(1:m)

  end do

  close ( unit = input_unit )

  return
end
subroutine i4mat_header_read ( input_filename, m, n )

!*****************************************************************************80
!
!! I4MAT_HEADER_READ reads the header from an I4MAT.
!
!  Discussion:
!
!    An I4MAT is an array of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILENAME, the name of the input file.
!
!    Output, integer ( kind = 4 ) M, spatial dimension.
!
!    Output, integer ( kind = 4 ) N, the number of points.
!
  implicit none

  character ( len = * )  input_filename
  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) n

  call file_column_count ( input_filename, m )

  if ( m <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4MAT_HEADER_READ - Fatal error!'
    write ( *, '(a)' ) '  There was some kind of I/O problem while trying'
    write ( *, '(a)' ) '  to count the number of data columns in'
    write ( *, '(a)' ) '  the file "' // trim ( input_filename ) // '".'
    stop
  end if

  call file_row_count ( input_filename, n )

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4MAT_HEADER_READ - Fatal error!'
    write ( *, '(a)' ) '  There was some kind of I/O problem while trying'
    write ( *, '(a)' ) '  to count the number of data rows in'
    write ( *, '(a)' ) '  the file "' // trim ( input_filename ) // '".'
    stop
  end if

  return
end
subroutine i4mat_transpose_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! I4MAT_TRANSPOSE_PRINT_SOME prints some of the transpose of an I4MAT.
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
!    Input, integer M, N, the number of rows and columns.
!
!    Input, integer A(M,N), an M by N matrix to be printed.
!
!    Input, integer ILO, JLO, the first row and column to print.
!
!    Input, integer IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none

  integer, parameter :: incx = 10
  integer m
  integer n

  integer a(m,n)
  character ( len = 7 ) ctemp(incx)
  integer i
  integer i2
  integer i2hi
  integer i2lo
  integer ihi
  integer ilo
  integer inc
  integer j
  integer j2hi
  integer j2lo
  integer jhi
  integer jlo
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

  return
end
subroutine jacobian_adjust_dirichlet ( node_num, node_xy, node_p_variable, &
  node_u_variable, node_v_variable, node_p_condition, &
  node_u_condition, node_v_condition, variable_num, ib, a )

!*****************************************************************************80
!
!! JACOBIAN_ADJUST_DIRICHLET adjusts the jacobian for Dirichlet conditions.
!
!  Modified:
!
!    21 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, real ( kind = 8 ) NODE_XY(2,NODE_NUM), the coordinates of nodes.
!
!    Input, integer NODE_P_VARIABLE(NODE_NUM),
!    is the index of the pressure variable associated with the node,
!    or -1 if there is no associated pressure variable.
!
!    Input, integer NODE_U_VARIABLE(NODE_NUM),
!    is the index of the horizontal velocity variable associated with the node.
!
!    Input, integer NODE_V_VARIABLE(NODE_NUM),
!    is the index of the vertical velocity variable associated with the node.
!
!    Input, integer NODE_P_CONDITION(NODE_NUM),
!    indicates the condition used to determine pressure at a node.
!    0, there is no condition at this node.
!    1, a finite element equation is used;
!    2, a Dirichlet condition is used.
!    3, a Neumann condition is used.
!
!    Input, integer NODE_U_CONDITION(NODE_NUM),
!    indicates the condition used to determine horizontal velocity at a node.
!    0, there is no condition at this node.
!    1, a finite element equation is used;
!    2, a Dirichlet condition is used.
!    3, a Neumann condition is used.
!
!    Input, integer NODE_V_CONDITION(NODE_NUM),
!    indicates the condition used to determine vertical velocity at a node.
!    0, there is no condition at this node.
!    1, a finite element equation is used;
!    2, a Dirichlet condition is used.
!    3, a Neumann condition is used.
!
!    Input, integer VARIABLE_NUM, the number of variables.
!
!    Input, integer IB, the half-bandwidth of the matrix.
!
!    Input/output, real ( kind = 8 ) A(3*IB+1,VARIABLE_NUM), the VARIABLE_NUM
!    by VARIABLE_NUM coefficient matrix, stored in a compressed format;
!    on output, the matrix has been adjusted for Dirichlet boundary conditions.
!
  implicit none

  integer ib
  integer node_num
  integer variable_num

  real ( kind = 8 ), dimension(3*ib+1,variable_num) :: a
  integer column
  integer column_high
  integer column_low
  integer, parameter :: DIRICHLET = 2
  integer ip
  integer iu
  integer iv
  integer node
  integer node_p_condition(node_num)
  integer node_p_variable(node_num)
  integer node_u_condition(node_num)
  integer node_u_variable(node_num)
  integer node_v_condition(node_num)
  integer node_v_variable(node_num)
  real ( kind = 8 ) node_xy(2,node_num)

  do node = 1, node_num

    iu = node_u_variable(node)
    iv = node_v_variable(node)
    ip = node_p_variable(node)

    if ( node_u_condition(node) == DIRICHLET ) then

      column_low = max ( iu - ib, 1 )
      column_high = min ( iu + ib, variable_num )

      do column = column_low, column_high
        a(iu-column+2*ib+1,column) = 0.0D+00
      end do
      a(2*ib+1,iu) = 1.0D+00

    end if

    if ( node_v_condition(node) == DIRICHLET ) then

      column_low = max ( iv - ib, 1 )
      column_high = min ( iv + ib, variable_num )

      do column = column_low, column_high
        a(iv-column+2*ib+1,column) = 0.0D+00
      end do
      a(2*ib+1,iv) = 1.0D+00

    end if

    if ( 0 < ip ) then

      if ( node_p_condition(node) == DIRICHLET ) then

        column_low = max ( ip - ib, 1 )
        column_high = min ( ip + ib, variable_num )

        do column = column_low, column_high
          a(ip-column+2*ib+1,column) = 0.0D+00
        end do
        a(2*ib+1,ip) = 1.0D+00

      end if

    end if

  end do

  return
end
subroutine jacobian_fem ( node_num, node_xy, element_num, &
  element_node, quad_num, node_u_variable, node_v_variable, &
  node_p_variable, variable_num, nu, c, ib, a )

!*****************************************************************************80
!
!! JACOBIAN_FEM evaluates the Navier Stokes jacobian matrix.
!
!  Discussion:
!
!    The matrix is known to be banded.  A special matrix storage format
!    is used to reduce the space required.  Details of this format are
!    discussed in the routine DGB_FA.
!
!    The Navier Stokes equations in weak form are:
!
!      Integral ( nu * ( dBdx(I) * dUdx + dBdy(I) * dUdy )
!        + B(I) * ( ( U * dUdx + V * dUdy ) + dPdx - U_RHS ) ) = 0
!
!      Integral ( nu * ( dBdx(I) * dVdx + dBdy(I) * dVdy )
!        + B(I) * ( ( U * dVdx + V * dVdy ) + dPdy - V_RHS ) ) = 0
!
!      Integral ( Q(I) * ( dUdx + dVdy - P_RHS ) ) = 0
!
!    This routine sets up the matrix as though every degree of freedom
!    were unconstrained.  Adjustments for boundary conditions and other
!    constraints should be made after calling this routine.
!
!  Modified:
!
!    08 September 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, real ( kind = 8 ) NODE_XY(2,NODE_NUM), the
!    coordinates of nodes.
!
!    Input, integer ELEMENT_NUM, the number of elements.
!
!    Input, integer ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM);
!    ELEMENT_NODE(I,J) is the global index of local node I in element J.
!
!    Input, integer QUAD_NUM, the number of quadrature points used in assembly.
!
!    Input, integer NODE_U_VARIABLE(NODE_NUM),
!    is the index of the horizontal velocity variable associated with the node.
!
!    Input, integer NODE_V_VARIABLE(NODE_NUM),
!    is the index of the vertical velocity variable associated with the node.
!
!    Input, integer NODE_P_VARIABLE(NODE_NUM),
!    is the index of the pressure variable associated with the node,
!    or -1 if there is no associated pressure variable.
!
!    Input, integer VARIABLE_NUM, the number of variables.
!
!    Input, real ( kind = 8 ) NU, the kinematic viscosity.
!
!    Input, real ( kind = 8 ) C(VARIABLE_NUM), the finite element
!    coefficients of an approximate solution of the Navier Stokes equations.
!
!    Input, integer IB, the bandwidth of the jacobian.
!
!    Output, real ( kind = 8 ) A(3*IB+1,VARIABLE_NUM), the VARIABLE_NUM
!    by VARIABLE_NUM Navier Stokes jacobian, stored in a general band
!    matrix format.
!
  implicit none

  integer ib
  integer node_num
  integer quad_num
  integer element_num
  integer variable_num

  real ( kind = 8 ) a(3*ib+1,variable_num)
  real ( kind = 8 ) area
  real ( kind = 8 ) b(6,quad_num)
  real ( kind = 8 ), dimension(node_num) :: c
  real ( kind = 8 ) cp(3)
  real ( kind = 8 ) cu(6)
  real ( kind = 8 ) cv(6)
  real ( kind = 8 ) dbdx(6,quad_num)
  real ( kind = 8 ) dbdy(6,quad_num)
  real ( kind = 8 ) dpdx(quad_num)
  real ( kind = 8 ) dpdy(quad_num)
  real ( kind = 8 ) dqdx(3,quad_num)
  real ( kind = 8 ) dqdy(3,quad_num)
  real ( kind = 8 ) dudx(quad_num)
  real ( kind = 8 ) dudy(quad_num)
  real ( kind = 8 ) dvdx(quad_num)
  real ( kind = 8 ) dvdy(quad_num)
  integer element
  integer, dimension(6,element_num) :: element_node
  integer i
  integer ip(3)
  integer iu(6)
  integer iv(6)
  integer j
  integer node_p_variable(node_num)
  integer node_u_variable(node_num)
  integer node_v_variable(node_num)
  real ( kind = 8 ), dimension(2,node_num) :: node_xy
  real ( kind = 8 ) nu
  real ( kind = 8 ) p(quad_num)
  real ( kind = 8 ) p_rhs(quad_num)
  real ( kind = 8 ) q(3,quad_num)
  integer quad
  real ( kind = 8 ), dimension(quad_num) :: quad_w
  real ( kind = 8 ), dimension(2,quad_num) :: quad_xy
  real ( kind = 8 ), dimension(2,3) :: t3
  real ( kind = 8 ), dimension(2,6) :: t6
  integer test
  real ( kind = 8 ) triangle_area_2d
  real ( kind = 8 ) u(quad_num)
  real ( kind = 8 ) u_rhs(quad_num)
  real ( kind = 8 ) v(quad_num)
  real ( kind = 8 ) v_rhs(quad_num)
  real ( kind = 8 ) w(quad_num)
  real ( kind = 8 ), dimension(2,quad_num) :: xy
!
!  Initialize the jacobian to zero.
!
  a(1:3*ib+1,1:variable_num) = 0.0D+00
!
!  Get the quadrature weights and nodes.
!
  call quad_rule ( quad_num, quad_w, quad_xy )
!
!  Consider all quantities associated with a given ELEMENT.
!
  do element = 1, element_num
!
!  Extract the nodes of the linear and quadratic triangles.
!
    t3(1:2,1:3) = node_xy(1:2,element_node(1:3,element))
    t6(1:2,1:6) = node_xy(1:2,element_node(1:6,element))
!
!  Map the quadrature points QUAD_XY to points XY in the physical element.
!
    call reference_to_physical_t6 ( t6, quad_num, quad_xy, xy )
    area = abs ( triangle_area_2d ( t3 ) )
    w(1:quad_num) = quad_w(1:quad_num) * area
!
!  Evaluate the basis functions at the quadrature points.
!
    call basis_mn_t6 ( t6, quad_num, xy, b, dbdx, dbdy )

    call basis_mn_t3 ( t3, quad_num, xy, q, dqdx, dqdy )
!
!  Extract the indices of the finite element coefficients for this element.
!
    iu(1:6) = node_u_variable(element_node(1:6,element))
    iv(1:6) = node_v_variable(element_node(1:6,element))
    ip(1:3) = node_p_variable(element_node(1:3,element))
!
!  Extract the finite element coefficients for this element.
!
    cu(1:6) = c(iu(1:6))
    cv(1:6) = c(iv(1:6))
    cp(1:3) = c(ip(1:3))
!
!  Evaluate the flowfield at each quadrature point.
!
     u(1:quad_num)   = matmul ( cu(1:6),  b(1:6,1:quad_num) )
    dudx(1:quad_num) = matmul ( cu(1:6), dbdx(1:6,1:quad_num) )
    dudy(1:quad_num) = matmul ( cu(1:6), dbdy(1:6,1:quad_num) )

     v(1:quad_num)   = matmul ( cv(1:6),  b(1:6,1:quad_num) )
    dvdx(1:quad_num) = matmul ( cv(1:6), dbdx(1:6,1:quad_num) )
    dvdy(1:quad_num) = matmul ( cv(1:6), dbdy(1:6,1:quad_num) )

     p(1:quad_num)   = matmul ( cp(1:3),  q(1:3,1:quad_num) )
    dpdx(1:quad_num) = matmul ( cp(1:3), dqdx(1:3,1:quad_num) )
    dpdy(1:quad_num) = matmul ( cp(1:3), dqdy(1:3,1:quad_num) )
!
!  dUeqn/dUcof,
!  dUeqn/dVcof,
!  dUeqn/dPcof.
!
    do i = 1, 6
      do j = 1, 6

        a(iu(i)-iu(j)+2*ib+1,iu(j)) = a(iu(i)-iu(j)+2*ib+1,iu(j)) + sum &
        ( w(1:quad_num) *                                               &
          (                                                             &
            nu * ( dbdx(j,1:quad_num) * dbdx(i,1:quad_num)              &
                 + dbdy(j,1:quad_num) * dbdy(i,1:quad_num) )            &
            +                                                           &
            ( b(j,1:quad_num) * dudx(1:quad_num)                        &
            + u(1:quad_num) * dbdx(j,1:quad_num)                        &
            + v(1:quad_num) * dbdy(j,1:quad_num) ) * b(i,1:quad_num)    &
          )                                                             &
        )

        a(iu(i)-iv(j)+2*ib+1,iv(j)) = a(iu(i)-iv(j)+2*ib+1,iv(j)) + sum &
        ( w(1:quad_num) * b(j,1:quad_num) * dudy(1:quad_num)            &
        * b(i,1:quad_num) )

      end do

      do j = 1, 3
        a(iu(i)-ip(j)+2*ib+1,ip(j)) = a(iu(i)-ip(j)+2*ib+1,ip(j)) + sum &
          ( w(1:quad_num) * dqdx(j,1:quad_num) * b(i,1:quad_num) )
      end do

    end do
!
!  dVeqn/dUcof,
!  dVeqn/dVcof,
!  dVeqn/dPcof.
!
    do i = 1, 6
      do j = 1, 6

        a(iv(i)-iu(j)+2*ib+1,iu(j)) = a(iv(i)-iu(j)+2*ib+1,iu(j)) + sum &
        ( w(1:quad_num) * b(j,1:quad_num) * dvdx(1:quad_num)            &
          * b(i,1:quad_num) )

        a(iv(i)-iv(j)+2*ib+1,iv(j)) = a(iv(i)-iv(j)+2*ib+1,iv(j)) + sum &
        ( w(1:quad_num) *                                               &
          (                                                             &
            nu * ( dbdx(j,1:quad_num) * dbdx(i,1:quad_num)              &
                 + dbdy(j,1:quad_num) * dbdy(i,1:quad_num) )            &
            +                                                           &
            ( u(1:quad_num) * dbdx(j,1:quad_num)                        &
            + b(j,1:quad_num) * dvdy(1:quad_num)                        &
            + v(1:quad_num) * dbdy(j,1:quad_num) )                      &
            * b(i,1:quad_num)                                           &
          )                                                             &
        )

      end do

      do j = 1, 3
        a(iv(i)-ip(j)+2*ib+1,ip(j)) = a(iv(i)-ip(j)+2*ib+1,ip(j)) + sum &
        ( w(1:quad_num) * dqdy(j,1:quad_num) * b(i,1:quad_num) )
      end do

    end do
!
!  dPeqn/dUcof,
!  dPeqn/dVcof,
!
    do i = 1, 3
      do j = 1, 6

        a(ip(i)-iu(j)+2*ib+1,iu(j)) = a(ip(i)-iu(j)+2*ib+1,iu(j)) + sum &
        ( w(1:quad_num) * dbdx(j,1:quad_num) * q(i,1:quad_num) )

        a(ip(i)-iv(j)+2*ib+1,iv(j)) = a(ip(i)-iv(j)+2*ib+1,iv(j)) + sum &
        ( w(1:quad_num) * dbdy(j,1:quad_num) * q(i,1:quad_num) )

      end do
    end do

  end do

  return
end
subroutine lvec_print ( n, a, title )

!*****************************************************************************80
!
!! LVEC_PRINT prints a logical vector.
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
!    Input, integer N, the number of components of the vector.
!
!    Input, logical A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  implicit none

  integer n

  logical a(n)
  integer i
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i8,2x,l1)' ) i, a(i)
  end do

  return
end
subroutine neumann_apply ( node_num, node_xy, node_p_variable, &
  node_u_variable, node_v_variable, node_p_condition, &
  node_u_condition, node_v_condition, variable_num, f )

!*****************************************************************************80
!
!! NEUMANN_APPLY accounts for Neumann boundary conditions.
!
!  Discussion:
!
!    At the moment, this program only allows Neumann boundary conditions
!    of the form
!
!      dU/dn = 0
!      dV/dn = 0
!      dP/dn = 0
!
!    For such conditions, there is NO change necessary to the linear system.
!    So this routine actually does nothing.  It is here as preparation
!    for later treatment of nonzero Neumann conditions.
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
!    Input, real ( kind = 8 ) NODE_XY(2,NODE_NUM), the coordinates of nodes.
!
!    Input, integer NODE_P_VARIABLE(NODE_NUM),
!    is the index of the pressure variable associated with the node,
!    or -1 if there is no associated pressure variable.
!
!    Input, integer NODE_U_VARIABLE(NODE_NUM),
!    is the index of the horizontal velocity variable associated with the node.
!
!    Input, integer NODE_V_VARIABLE(NODE_NUM),
!    is the index of the vertical velocity variable associated with the node.
!
!    Input, integer NODE_P_CONDITION(NODE_NUM),
!    indicates the condition used to determine pressure at a node.
!    0, there is no condition at this node.
!    1, a finite element equation is used;
!    2, a Dirichlet condition is used.
!    3, a Neumann condition is used.
!
!    Input, integer NODE_U_CONDITION(NODE_NUM),
!    indicates the condition used to determine horizontal velocity at a node.
!    0, there is no condition at this node.
!    1, a finite element equation is used;
!    2, a Dirichlet condition is used.
!    3, a Neumann condition is used.
!
!    Input, integer NODE_V_CONDITION(NODE_NUM),
!    indicates the condition used to determine vertical velocity at a node.
!    0, there is no condition at this node.
!    1, a finite element equation is used;
!    2, a Dirichlet condition is used.
!    3, a Neumann condition is used.
!
!    Input, integer VARIABLE_NUM, the number of variables.
!
!    Input/output, real ( kind = 8 ) F(VARIABLE_NUM), the right hand side.
!    On output, the right hand side has been adjusted for Dirichlet
!    boundary conditions.
!
  implicit none

  integer node_num
  integer variable_num

  real ( kind = 8 ), dimension(variable_num) :: f
  integer ip
  integer iu
  integer iv
  integer, parameter :: NEUMANN = 3
  integer node
  integer node_p_condition(node_num)
  integer node_p_variable(node_num)
  integer node_u_condition(node_num)
  integer node_u_variable(node_num)
  integer node_v_condition(node_num)
  integer node_v_variable(node_num)
  real ( kind = 8 ), dimension(2,node_num) :: node_xy
  real ( kind = 8 ) p_bc(node_num)
  real ( kind = 8 ) u_bc(node_num)
  real ( kind = 8 ) v_bc(node_num)
!
!  The user routine supplies a right hand side value for a possible
!  Neumann condition at EVERY node.
!
  call neumann_condition ( node_num, node_xy, u_bc, v_bc, p_bc )

  do node = 1, node_num

    iu = node_u_variable(node)
    iv = node_v_variable(node)
    ip = node_p_variable(node)

    if ( node_u_condition(node) == NEUMANN ) then

!     f(iu) = f(iu) + line integral

    end if

    if ( node_v_condition(node) == NEUMANN ) then

!     f(iv) = f(iv) + line integral

    end if

    if ( 0 < ip ) then

      if ( node_p_condition(node) == NEUMANN ) then

!       f(ip) = f(ip) + line integral

      end if

    end if

  end do

  return
end
subroutine nodes3_write ( file_name, node_num, node_xy, node_type )

!*****************************************************************************80
!
!! NODES3_WRITE writes the pressure nodes to a file.
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
!    Input, character ( len = * ) FILE_NAME, the file name.
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, real ( kind = 8 ) NODE_XY(2,NODE_NUM), the coordinates of the nodes.
!
!    Input, integer NODE_TYPE(NODE_NUM), determines if the node is a
!    vertex or midside node.
!    1, the node is a vertex (P, U, V variables are associated with it).
!    2, the node is a midside node (only U and V variables are associated.)
!
  implicit none

  integer node_num

  character ( len = * ) :: file_name
  integer file_status
  integer file_unit
  integer node
  integer, dimension(node_num) :: node_type
  real ( kind = 8 ), dimension(2,node_num) :: node_xy

  call get_unit ( file_unit )

  open ( unit = file_unit, file = file_name, status = 'replace', &
    iostat = file_status )

  if ( file_status /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'NODES3_WRITE - Warning!'
    write ( *, '(a)' ) '  Could not write the file "' &
      // trim ( file_name ) // '".'
    return
  end if

  do node = 1, node_num

    if ( node_type(node) == 1 ) then
      write ( file_unit, '(2x,g14.6,2x,g14.6)' ) node_xy(1:2,node)
    end if

  end do

  close ( unit = file_unit )

  return
end
subroutine points_plot ( file_name, node_num, node_xy, node_label )

!*****************************************************************************80
!
!! POINTS_PLOT plots a pointset.
!
!  Modified:
!
!    27 September 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the output file.
!
!    Input, integer NODE_NUM, the number of points.
!
!    Input, real ( kind = 8 ) NODE_XY(2,NODE_NUM), the nodes.
!
!    Input, logical NODE_LABEL, is TRUE if the nodes should be labeled.
!
  implicit none

  integer node_num

  integer :: circle_size
  integer delta
  character ( len = * ) file_name
  integer file_unit
  integer i
  integer ios
  integer node
  logical node_label
  real ( kind = 8 ) node_xy(2,node_num)
  character ( len = 40 ) string
  real ( kind = 8 ) x_max
  real ( kind = 8 ) x_min
  integer x_ps
  integer :: x_ps_max = 576
  integer :: x_ps_max_clip = 594
  integer :: x_ps_min = 36
  integer :: x_ps_min_clip = 18
  real ( kind = 8 ) x_scale
  real ( kind = 8 ) y_max
  real ( kind = 8 ) y_min
  integer y_ps
  integer :: y_ps_max = 666
  integer :: y_ps_max_clip = 684
  integer :: y_ps_min = 126
  integer :: y_ps_min_clip = 108
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
  if ( node_num <= 200 ) then
    circle_size = 5
  else if ( node_num <= 500 ) then
    circle_size = 4
  else if ( node_num <= 1000 ) then
    circle_size = 3
  else if ( node_num <= 5000 ) then
    circle_size = 2
  else
    circle_size = 1
  end if

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
subroutine pressure3_write ( file_name, node_num, node_p_variable, &
  variable_num, node_c )

!*****************************************************************************80
!
!! PRESSURE3_WRITE writes the pressures to a file.
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
!    Input, character ( len = * ) FILE_NAME, the file name.
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, integer NODE_P_VARIABLE(NODE_NUM),
!    is the index of the pressure variable associated with the node,
!    or -1 if there is no associated pressure variable.
!
!    Input, integer VARIABLE_NUM, the number of variables.
!
!    Input, real ( kind = 8 ) NODE_C(VARIABLE_NUM), the finite element
!    coefficients.
!
  implicit none

  integer node_num
  integer variable_num

  character ( len = * ) :: file_name
  integer file_status
  integer file_unit
  integer node
  real ( kind = 8 ), dimension(variable_num) :: node_c
  integer, dimension(node_num) :: node_p_variable
  integer variable

  call get_unit ( file_unit )

  open ( unit = file_unit, file = file_name, status = 'replace', &
    iostat = file_status )

  if ( file_status /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PRESSURE3_WRITE - Warning!'
    write ( *, '(a)' ) '  Could not write the file "' &
      // trim ( file_name ) // '".'
    return
  end if

  do node = 1, node_num

    variable = node_p_variable(node)

    if ( 0 <  variable ) then
      write ( file_unit, '(2x,g14.6)' ) node_c(variable)
    end if

  end do

  close ( unit = file_unit )

  return
end
subroutine quad_rule ( quad_num, quad_w, quad_xy )

!*****************************************************************************80
!
!! QUAD_RULE sets the quadrature rule for assembly.
!
!  Discussion:
!
!    The quadrature rule is given for a reference element.
!
!      0 <= X,
!      0 <= Y, and
!      X + Y <= 1.
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
!    The rules have the following precision:
!
!    QUAD_NUM  Precision
!
!     1        1
!     3        2
!     4        3
!     6        4
!     7        5
!     9        6
!    13        7
!
!  Modified:
!
!    18 July 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer QUAD_NUM, the number of quadrature nodes.
!
!    Output, real ( kind = 8 ) QUAD_W(QUAD_NUM), the quadrature weights.
!
!    Output, real ( kind = 8 ) QUAD_XY(2,QUAD_NUM),
!    the coordinates of the quadrature nodes.
!
  implicit none

  integer quad_num

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) d
  real ( kind = 8 ) e
  real ( kind = 8 ) f
  real ( kind = 8 ) g
  real ( kind = 8 ) h
  real ( kind = 8 ), dimension(quad_num) :: quad_w
  real ( kind = 8 ), dimension(2,quad_num) :: quad_xy
  real ( kind = 8 ) t
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) w

  if ( quad_num == 1 ) then

    quad_xy(1:2,1:quad_num) = reshape ( (/ &
      1.0D+00 / 3.0D+00, 1.0D+00 / 3.0D+00 /), (/ 2, quad_num /) )

    quad_w(1:quad_num) = 1.0D+00

  else if ( quad_num == 3 ) then

    quad_xy(1:2,1:quad_num) = reshape ( (/ &
      0.5D+00, 0.0D+00, &
      0.5D+00, 0.5D+00, &
      0.0D+00, 0.5D+00 /), (/ 2, quad_num /) )

    quad_w(1:quad_num) = 1.0D+00 / 3.0D+00

  else if ( quad_num == 4 ) then

    a =   6.0D+00 / 30.0D+00
    b =  10.0D+00 / 30.0D+00
    c =  18.0D+00 / 30.0D+00

    d =  25.0D+00 / 48.0D+00
    e = -27.0D+00 / 48.0D+00

    quad_xy(1:2,1:quad_num) = reshape ( (/ &
      b, b, &
      c, a, &
      a, c, &
      a, a /), (/ 2, quad_num /) )

    quad_w(1:quad_num) = (/ e, d, d, d /)

  else if ( quad_num == 6 ) then

    a = 0.816847572980459D+00
    b = 0.091576213509771D+00
    c = 0.108103018168070D+00
    d = 0.445948490915965D+00
    v = 0.109951743655322D+00
    w = 0.223381589678011D+00

    quad_xy(1:2,1:quad_num) = reshape ( (/ &
      a, b, &
      b, a, &
      b, b, &
      c, d, &
      d, c, &
      d, d /), (/ 2, quad_num /) )

    quad_w(1:quad_num) = (/ v, v, v, w, w, w /)

  else if ( quad_num == 7 ) then

    a = 1.0D+00 / 3.0D+00
    b = ( 9.0D+00 + 2.0D+00 * sqrt ( 15.0D+00 ) ) / 21.0D+00
    c = ( 6.0D+00 -       sqrt ( 15.0D+00 ) ) / 21.0D+00
    d = ( 9.0D+00 - 2.0D+00 * sqrt ( 15.0D+00 ) ) / 21.0D+00
    e = ( 6.0D+00 +       sqrt ( 15.0D+00 ) ) / 21.0D+00
    u = 0.225D+00
    v = ( 155.0D+00 - sqrt ( 15.0D+00 ) ) / 1200.0D+00
    w = ( 155.0D+00 + sqrt ( 15.0D+00 ) ) / 1200.0D+00

    quad_xy(1:2,1:quad_num) = reshape ( (/ &
      a, a, &
      b, c, &
      c, b, &
      c, c, &
      d, e, &
      e, d, &
      e, e /), (/ 2, quad_num /) )

    quad_w(1:quad_num) = (/ u, v, v, v, w, w, w /)

  else if ( quad_num == 9 ) then

    a = 0.124949503233232D+00
    b = 0.437525248383384D+00
    c = 0.797112651860071D+00
    d = 0.165409927389841D+00
    e = 0.037477420750088D+00

    u = 0.205950504760887D+00
    v = 0.063691414286223D+00

    quad_xy(1:2,1:quad_num) = reshape ( (/ &
      a, b, &
      b, a, &
      b, b, &
      c, d, &
      c, e, &
      d, c, &
      d, e, &
      e, c, &
      e, d /), (/ 2, quad_num /) )

    quad_w(1:quad_num) = (/ u, u, u, v, v, v, v, v, v /)

  else if ( quad_num == 13 ) then

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

    quad_xy(1:2,1:quad_num) = reshape ( (/ &
      h, h, &
      a, b, &
      b, a, &
      b, b, &
      c, d, &
      d, c, &
      d, d, &
      e, f, &
      e, g, &
      f, e, &
      f, g, &
      g, e, &
      g, f /), (/ 2, quad_num /) )

    quad_w(1:quad_num) = (/ w, t, t, t, u, u, u, v, v, v, v, v, v /)

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'QUAD_RULE - Fatal error!'
    write ( *, '(a,i8)' ) '  No rule is available of order QUAD_NUM = ', &
      quad_num
    stop

  end if

  return
end
subroutine r8mat_data_read ( input_filename, m, n, table )

!*****************************************************************************80
!
!! R8MAT_DATA_READ reads data from an R8MAT file.
!
!  Discussion:
!
!    An R8MAT is an array of R8 values.
!
!  Discussion:
!
!    The file may contain more than N points, but this routine will
!    return after reading N of them.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILENAME, the name of the input file.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Output, real ( kind = 8 ) TABLE(M,N), the data.
!
  implicit none

  integer   ( kind = 4 )   m
  integer   ( kind = 4 )   n

  integer   ( kind = 4 )   ierror
  character ( len = * )    input_filename
  integer   ( kind = 4 )   input_status
  integer   ( kind = 4 )   input_unit
  integer   ( kind = 4 )   j
  character ( len = 255 )  line
  real      ( kind = 8 )   table(m,n)
  real      ( kind = 8 )   x(m)

  ierror = 0

  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_filename, status = 'old', &
    iostat = input_status )

  if ( input_status /= 0 ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_DATA_READ - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the input file "' // &
      trim ( input_filename ) // '" on unit ', input_unit
    stop
  end if

  j = 0

  do while ( j < n )

    read ( input_unit, '(a)', iostat = input_status ) line

    if ( input_status /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8MAT_DATA_READ - Fatal error!'
      write ( *, '(a)' ) '  Error while reading lines of data.'
      write ( *, '(a,i8)' ) '  Number of values expected per line M = ', m
      write ( *, '(a,i8)' ) '  Number of data lines read, J =         ', j
      write ( *, '(a,i8)' ) '  Number of data lines needed, N =       ', n
      stop
    end if

    if ( line(1:1) == '#' .or. len_trim ( line ) == 0 ) then
      cycle
    end if

    call s_to_r8vec ( line, m, x, ierror )

    if ( ierror /= 0 ) then
      cycle
    end if

    j = j + 1

    table(1:m,j) = x(1:m)

  end do

  close ( unit = input_unit )

  return
end
subroutine r8mat_header_read ( input_filename, m, n )

!*****************************************************************************80
!
!! R8MAT_HEADER_READ reads the header from an R8MAT file.
!
!  Discussion:
!
!    An R8MAT is an array of R8 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILENAME, the name of the input file.
!
!    Output, integer ( kind = 4 ) M, spatial dimension.
!
!    Output, integer ( kind = 4 ) N, the number of points.
!
  implicit none

  character ( len = * )  input_filename
  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) n

  call file_column_count ( input_filename, m )

  if ( m <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_HEADER_READ - Fatal error!'
    write ( *, '(a)' ) '  There was some kind of I/O problem while trying'
    write ( *, '(a)' ) '  to count the number of data columns in'
    write ( *, '(a)' ) '  the file "' // trim ( input_filename ) // '".'
    stop
  end if

  call file_row_count ( input_filename, n )

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_HEADER_READ - Fatal error!'
    write ( *, '(a)' ) '  There was some kind of I/O problem while trying'
    write ( *, '(a)' ) '  to count the number of data rows in'
    write ( *, '(a)' ) '  the file "' // trim ( input_filename ) // '".'
    stop
  end if

  return
end
subroutine r8mat_transpose_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8MAT_TRANSPOSE_PRINT_SOME prints some of an R8MAT, transposed.
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
!    Input, integer M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ILO, JLO, the first row and column to print.
!
!    Input, integer IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none

  integer, parameter :: incx = 5
  integer m
  integer n

  real ( kind = 8 ) a(m,n)
  character ( len = 14 ) ctemp(incx)
  integer i
  integer i2
  integer i2hi
  integer i2lo
  integer ihi
  integer ilo
  integer inc
  integer j
  integer j2hi
  integer j2lo
  integer jhi
  integer jlo
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

  return
end
subroutine r8vec_print_some ( n, a, i_lo, i_hi, title )

!*****************************************************************************80
!
!! R8VEC_PRINT_SOME prints "some" of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8 values.
!
!  Modified:
!
!    16 October 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries of the vector.
!
!    Input, real ( kind = 8 ) A(N), the vector to be printed.
!
!    Input, integer I_LO, I_HI, the first and last indices to print.
!    The routine expects 1 <= I_LO <= I_HI <= N.
!
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none

  integer n

  real ( kind = 8 ) a(n)
  integer i
  integer i_hi
  integer i_lo
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '
  do i = max ( i_lo, 1 ), min ( i_hi, n )
    write ( *, '(2x,i8,2x,g14.8)' ) i, a(i)
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
!    Input, integer N, the number of objects to transform.
!
!    Input, real ( kind = 8 ) REF(2,N), points in the reference triangle.
!
!    Output, real ( kind = 8 ) PHY(2,N), corresponding points in the
!    physical triangle.
!
  implicit none

  integer n

  real ( kind = 8 ) a(2)
  real ( kind = 8 ) b(2)
  real ( kind = 8 ) c(2)
  real ( kind = 8 ) d(2)
  real ( kind = 8 ) e(2)
  real ( kind = 8 ) f(2)
  integer i
  real ( kind = 8 ) phy(2,n)
  real ( kind = 8 ) ref(2,n)
  real ( kind = 8 ) t(2,6)

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
subroutine residual_adjust_dirichlet ( node_num, node_xy, node_p_variable, &
  node_u_variable, node_v_variable, node_p_condition, &
  node_u_condition, node_v_condition, variable_num, node_c, node_r )

!*****************************************************************************80
!
!! RESIDUAL_ADJUST_DIRICHLET adjusts the residual for Dirichlet conditions.
!
!  Modified:
!
!    10 September 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, real ( kind = 8 ) NODE_XY(2,NODE_NUM), the coordinates of nodes.
!
!    Input, integer NODE_P_VARIABLE(NODE_NUM),
!    is the index of the pressure variable associated with the node,
!    or -1 if there is no associated pressure variable.
!
!    Input, integer NODE_U_VARIABLE(NODE_NUM),
!    is the index of the horizontal velocity variable associated with the node.
!
!    Input, integer NODE_V_VARIABLE(NODE_NUM),
!    is the index of the vertical velocity variable associated with the node.
!
!    Input, integer NODE_P_CONDITION(NODE_NUM),
!    indicates the condition used to determine pressure at a node.
!    0, there is no condition at this node.
!    1, a finite element equation is used;
!    2, a Dirichlet condition is used.
!    3, a Neumann condition is used.
!
!    Input, integer NODE_U_CONDITION(NODE_NUM),
!    indicates the condition used to determine horizontal velocity at a node.
!    0, there is no condition at this node.
!    1, a finite element equation is used;
!    2, a Dirichlet condition is used.
!    3, a Neumann condition is used.
!
!    Input, integer NODE_V_CONDITION(NODE_NUM),
!    indicates the condition used to determine vertical velocity at a node.
!    0, there is no condition at this node.
!    1, a finite element equation is used;
!    2, a Dirichlet condition is used.
!    3, a Neumann condition is used.
!
!    Input, integer VARIABLE_NUM, the number of variables.
!
!    Input, real ( kind = 8 ) NODE_C(VARIABLE_NUM), the finite element
!    coefficient vector.
!
!    Input/output, real ( kind = 8 ) NODE_R(VARIABLE_NUM), the
!    residual vector; on output, the residual has been adjusted for
!    Dirichlet boundary conditions.
!
  implicit none

  integer node_num

  integer, parameter :: DIRICHLET = 2
  integer ip
  integer iu
  integer iv
  integer node
  real ( kind = 8 ), dimension(node_num) :: node_c
  integer node_p_condition(node_num)
  integer node_p_variable(node_num)
  real ( kind = 8 ), dimension(node_num) :: node_r
  integer node_u_condition(node_num)
  integer node_u_variable(node_num)
  integer node_v_condition(node_num)
  integer node_v_variable(node_num)
  real ( kind = 8 ), dimension(2,node_num) :: node_xy
  real ( kind = 8 ) p_bc(node_num)
  real ( kind = 8 ) u_bc(node_num)
  real ( kind = 8 ) v_bc(node_num)
  integer variable_num

  call dirichlet_condition ( node_num, node_xy, u_bc, v_bc, p_bc )

  do node = 1, node_num

    iu = node_u_variable(node)
    iv = node_v_variable(node)
    ip = node_p_variable(node)

    if ( node_u_condition(node) == DIRICHLET ) then
      node_r(iu) = node_c(iu) - u_bc(node)
    end if

    if ( node_v_condition(node) == DIRICHLET ) then
      node_r(iv) = node_c(iv) - v_bc(node)
    end if

    if ( 0 < ip ) then
      if ( node_p_condition(node) == DIRICHLET ) then
        node_r(ip) = node_c(ip) - p_bc(node)
      end if
    end if

  end do

  return
end
subroutine residual_adjust_neumann ( node_num, node_xy, node_p_variable, &
  node_u_variable, node_v_variable, node_p_condition, &
  node_u_condition, node_v_condition, variable_num, node_c, node_r )

!*****************************************************************************80
!
!! RESIDUAL_ADJUST_NEUMANN adjusts the residual for Neumann conditions.
!
!  Discussion:
!
!    At the moment, we are only allowing zero Neumann conditions.
!    In that case, no adjustment to the residual is necessary,
!    so this routine is just a placeholder.
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
!    Input, real ( kind = 8 ) NODE_XY(2,NODE_NUM), the coordinates of nodes.
!
!    Input, integer NODE_P_VARIABLE(NODE_NUM),
!    is the index of the pressure variable associated with the node,
!    or -1 if there is no associated pressure variable.
!
!    Input, integer NODE_U_VARIABLE(NODE_NUM),
!    is the index of the horizontal velocity variable associated with the node.
!
!    Input, integer NODE_V_VARIABLE(NODE_NUM),
!    is the index of the vertical velocity variable associated with the node.
!
!    Input, integer NODE_P_CONDITION(NODE_NUM),
!    indicates the condition used to determine pressure at a node.
!    0, there is no condition at this node.
!    1, a finite element equation is used;
!    2, a Dirichlet condition is used.
!    3, a Neumann condition is used.
!
!    Input, integer NODE_U_CONDITION(NODE_NUM),
!    indicates the condition used to determine horizontal velocity at a node.
!    0, there is no condition at this node.
!    1, a finite element equation is used;
!    2, a Dirichlet condition is used.
!    3, a Neumann condition is used.
!
!    Input, integer NODE_V_CONDITION(NODE_NUM),
!    indicates the condition used to determine vertical velocity at a node.
!    0, there is no condition at this node.
!    1, a finite element equation is used;
!    2, a Dirichlet condition is used.
!    3, a Neumann condition is used.
!
!    Input, integer VARIABLE_NUM, the number of variables.
!
!    Input, real ( kind = 8 ) NODE_C(VARIABLE_NUM), the finite element
!    coefficient vector.
!
!    Input/output, real ( kind = 8 ) NODE_R(VARIABLE_NUM), the
!    residual vector; on output, the residual has been adjusted for
!    Neumann boundary conditions.
!
  implicit none

  integer node_num

  integer ip
  integer iu
  integer iv
  integer, parameter :: NEUMANN = 3
  integer node
  real ( kind = 8 ), dimension(node_num) :: node_c
  integer node_p_condition(node_num)
  integer node_p_variable(node_num)
  real ( kind = 8 ), dimension(node_num) :: node_r
  integer node_u_condition(node_num)
  integer node_u_variable(node_num)
  integer node_v_condition(node_num)
  integer node_v_variable(node_num)
  real ( kind = 8 ), dimension(2,node_num) :: node_xy
  real ( kind = 8 ) p_bc(node_num)
  real ( kind = 8 ) u_bc(node_num)
  real ( kind = 8 ) v_bc(node_num)
  integer variable_num

  call neumann_condition ( node_num, node_xy, u_bc, v_bc, p_bc )

  do node = 1, node_num

    iu = node_u_variable(node)
    iv = node_v_variable(node)
    ip = node_p_variable(node)

    if ( node_u_condition(node) == NEUMANN ) then
!     node_r(iu) = node_r(iu) + line integral
    end if

    if ( node_v_condition(node) == NEUMANN ) then
!     node_r(iv) = node_r(iv) + line integral
    end if

    if ( 0 < ip ) then
      if ( node_p_condition(node) == NEUMANN ) then
!       node_r(ip) = node_r(ip) + line integral
      end if
    end if

  end do

  return
end
subroutine residual_fem ( node_num, node_xy, element_num, &
  element_node, quad_num, node_u_variable, node_v_variable, &
  node_p_variable, variable_num, nu, c, r )

!*****************************************************************************80
!
!! RESIDUAL_FEM evaluates the Navier Stokes residual function.
!
!  Discussion:
!
!    This routine computes the Galerkin residual for every basis function.
!    The result must be adjusted (elsewhere) because for certain
!    components, the Galerkin residual is replaced by a boundary condition
!    or other constraint.
!
!    The Navier Stokes equations in weak form are:
!
!      Integral ( nu * ( dBdx(I) * dUdx + dBdy(I) * dUdy )
!        + B(I) * ( ( U * dUdx + V * dUdy ) + dPdx - U_RHS ) ) = 0
!
!      Integral ( nu * ( dBdx(I) * dVdx + dBdy(I) * dVdy )
!        + B(I) * ( ( U * dVdx + V * dVdy ) + dPdy - V_RHS ) ) = 0
!
!      Integral ( Q(I) * ( dUdx + dVdy - P_RHS ) ) = 0
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
!    Input, real ( kind = 8 ) NODE_XY(2,NODE_NUM), the
!    coordinates of nodes.
!
!    Input, integer ELEMENT_NUM, the number of elements.
!
!    Input, integer ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM);
!    ELEMENT_NODE(I,J) is the global index of local node I in element J.
!
!    Input, integer QUAD_NUM, the number of quadrature points used in assembly.
!
!    Input, integer NODE_U_VARIABLE(NODE_NUM),
!    is the index of the horizontal velocity variable associated with the node.
!
!    Input, integer NODE_V_VARIABLE(NODE_NUM),
!    is the index of the vertical velocity variable associated with the node.
!
!    Input, integer NODE_P_VARIABLE(NODE_NUM),
!    is the index of the pressure variable associated with the node,
!    or -1 if there is no associated pressure variable.
!
!    Input, integer VARIABLE_NUM, the number of variables.
!
!    Input, real ( kind = 8 ) NU, the kinematic viscosity.
!
!    Input, real ( kind = 8 ) C(VARIABLE_NUM), the finite element
!    coefficients of an approximate solution of the Navier Stokes equations.
!
!    Output, real ( kind = 8 ) R(VARIABLE_NUM), the Navier Stokes residual.
!
  implicit none

  integer node_num
  integer quad_num
  integer element_num
  integer variable_num

  real ( kind = 8 ) area
  real ( kind = 8 ) b(6,quad_num)
  real ( kind = 8 ), dimension(node_num) :: c
  real ( kind = 8 ) cp(3)
  real ( kind = 8 ) cu(6)
  real ( kind = 8 ) cv(6)
  real ( kind = 8 ) dbdx(6,quad_num)
  real ( kind = 8 ) dbdy(6,quad_num)
  real ( kind = 8 ) dpdx(quad_num)
  real ( kind = 8 ) dpdy(quad_num)
  real ( kind = 8 ) dqdx(3,quad_num)
  real ( kind = 8 ) dqdy(3,quad_num)
  real ( kind = 8 ) dudx(quad_num)
  real ( kind = 8 ) dudy(quad_num)
  real ( kind = 8 ) dvdx(quad_num)
  real ( kind = 8 ) dvdy(quad_num)
  integer element
  integer, dimension(6,element_num) :: element_node
  integer i
  integer ip(3)
  integer iu(6)
  integer iv(6)
  integer j
  integer node_p_variable(node_num)
  integer node_u_variable(node_num)
  integer node_v_variable(node_num)
  real ( kind = 8 ), dimension(2,node_num) :: node_xy
  real ( kind = 8 ) nu
  real ( kind = 8 ) p(quad_num)
  real ( kind = 8 ) p_rhs(quad_num)
  real ( kind = 8 ) q(3,quad_num)
  integer quad
  real ( kind = 8 ), dimension(quad_num) :: quad_w
  real ( kind = 8 ), dimension(2,quad_num) :: quad_xy
  real ( kind = 8 ), dimension(variable_num) :: r
  real ( kind = 8 ), dimension(2,3) :: t3
  real ( kind = 8 ), dimension(2,6) :: t6
  integer test
  real ( kind = 8 ) triangle_area_2d
  real ( kind = 8 ) u(quad_num)
  real ( kind = 8 ) u_rhs(quad_num)
  real ( kind = 8 ) v(quad_num)
  real ( kind = 8 ) v_rhs(quad_num)
  real ( kind = 8 ) w(quad_num)
  real ( kind = 8 ), dimension(2,quad_num) :: xy
!
!  Initialize the residual to zero.
!
  r(1:variable_num) = 0.0D+00
!
!  Get the quadrature weights and nodes.
!
  call quad_rule ( quad_num, quad_w, quad_xy )
!
!  Consider all quantities associated with a given ELEMENT.
!
  do element = 1, element_num
!
!  Extract the nodes of the linear and quadratic triangles.
!
    t3(1:2,1:3) = node_xy(1:2,element_node(1:3,element))
    t6(1:2,1:6) = node_xy(1:2,element_node(1:6,element))
!
!  Map the quadrature points QUAD_XY to points XY in the physical element.
!
    call reference_to_physical_t6 ( t6, quad_num, quad_xy, xy )
    area = abs ( triangle_area_2d ( t3 ) )
    w(1:quad_num) = quad_w(1:quad_num) * area
    call rhs ( quad_num, xy, u_rhs, v_rhs, p_rhs )
!
!  Evaluate the basis functions at the quadrature points.
!
    call basis_mn_t6 ( t6, quad_num, xy, b, dbdx, dbdy )
    call basis_mn_t3 ( t3, quad_num, xy, q, dqdx, dqdy )
!
!  Extract the indices of the finite element coefficients for this element.
!
    iu(1:6) = node_u_variable(element_node(1:6,element))
    iv(1:6) = node_v_variable(element_node(1:6,element))
    ip(1:3) = node_p_variable(element_node(1:3,element))
!
!  Extract the finite element coefficients for this element.
!
    cu(1:6) = c(iu(1:6))
    cv(1:6) = c(iv(1:6))
    cp(1:3) = c(ip(1:3))
!
!  Evaluate the flowfield at each quadrature point.
!
     u(1:quad_num)   = matmul ( cu(1:6),  b(1:6,1:quad_num) )
    dudx(1:quad_num) = matmul ( cu(1:6), dbdx(1:6,1:quad_num) )
    dudy(1:quad_num) = matmul ( cu(1:6), dbdy(1:6,1:quad_num) )

     v(1:quad_num)   = matmul ( cv(1:6),  b(1:6,1:quad_num) )
    dvdx(1:quad_num) = matmul ( cv(1:6), dbdx(1:6,1:quad_num) )
    dvdy(1:quad_num) = matmul ( cv(1:6), dbdy(1:6,1:quad_num) )

     p(1:quad_num)   = matmul ( cp(1:3),  q(1:3,1:quad_num) )
    dpdx(1:quad_num) = matmul ( cp(1:3), dqdx(1:3,1:quad_num) )
    dpdy(1:quad_num) = matmul ( cp(1:3), dqdy(1:3,1:quad_num) )
!
!  The horizontal momentum equation multiplied by Bi
!
    do i = 1, 6
      r(iu(i)) = r(iu(i)) + sum                          &
      ( w(1:quad_num) *                                  &
        (                                                &
          nu * ( dudx(1:quad_num) * dbdx(i,1:quad_num)   &
               + dudy(1:quad_num) * dbdy(i,1:quad_num) ) &
          +                                              &
          ( u(1:quad_num) * dudx(1:quad_num)             &
          + v(1:quad_num) * dudy(1:quad_num)             &
          + dpdx(1:quad_num) - u_rhs(1:quad_num) )       &
          * b(i,1:quad_num)                              &
        )                                                &
      )
    end do
!
!  The vertical momentum equation multiplied by Bi
!
    do i = 1, 6
      r(iv(i)) = r(iv(i)) + sum                          &
      ( w(1:quad_num) *                                  &
        (                                                &
          nu * ( dvdx(1:quad_num) * dbdx(i,1:quad_num)   &
               + dvdy(1:quad_num) * dbdy(i,1:quad_num) ) &
          +                                              &
          ( u(1:quad_num) * dvdx(1:quad_num)             &
          + v(1:quad_num) * dvdy(1:quad_num)             &
          + dpdy(1:quad_num) - v_rhs(1:quad_num) )       &
          * b(i,1:quad_num)                              &
        )                                                &
      )
    end do
!
!  The continuity equation multiplied by Qi
!
    do i = 1, 3
      r(ip(i)) = r(ip(i)) + sum                          &
      ( w(1:quad_num) *                                  &
        (                                                &
          ( dudx(1:quad_num) + dvdy(1:quad_num)          &
          - p_rhs(1:quad_num) )                          &
          * q(i,1:quad_num)                              &
        )                                                &
      )
    end do

  end do

  return
end
subroutine s_to_i4 ( s, ival, ierror, length )

!*****************************************************************************80
!
!! S_TO_I4 reads an I4 from a string.
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
!    Output, integer IVAL, the integer value read from the string.
!    If the string is blank, then IVAL will be returned 0.
!
!    Output, integer IERROR, an error flag.
!    0, no error.
!    1, an error occurred.
!
!    Output, integer LENGTH, the number of characters of S used to make IVAL.
!
  implicit none

  character c
  integer i
  integer ierror
  integer isgn
  integer istate
  integer ival
  integer length
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
        length = i - 1
        return
      end if

    end if

  end do
!
!  If we read all the characters in the string, see if we're OK.
!
  if ( istate == 2 ) then
    ival = isgn * ival
    length = len_trim ( s )
  else
    ierror = 1
    length = 0
  end if

  return
end
subroutine s_to_i4vec ( s, n, ivec, ierror )

!*****************************************************************************80
!
!! S_TO_I4VEC reads an integer vector from a string.
!
!  Modified:
!
!    08 October 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string to be read.
!
!    Input, integer N, the number of values expected.
!
!    Output, integer IVEC(N), the values read from the string.
!
!    Output, integer IERROR, error flag.
!    0, no errors occurred.
!    -K, could not read data for entries -K through N.
!
  implicit none

  integer n

  integer i
  integer ierror
  integer ilo
  integer ivec(n)
  integer length
  character ( len = * ) s

  i = 0
  ierror = 0
  ilo = 1

  do while ( i < n )

    i = i + 1

    call s_to_i4 ( s(ilo:), ivec(i), ierror, length )

    if ( ierror /= 0 ) then
      ierror = -i
      exit
    end if

    ilo = ilo + length

  end do

  return
end
subroutine s_to_r8 ( s, dval, ierror, length )

!*****************************************************************************80
!
!! S_TO_R8 reads an R8 from a string.
!
!  Discussion:
!
!    The routine will read as many characters as possible until it reaches
!    the end of the string, or encounters a character which cannot be
!    part of the number.
!
!    Legal input is:
!
!       1 blanks,
!       2 '+' or '-' sign,
!       2.5 blanks
!       3 integer part,
!       4 decimal point,
!       5 fraction part,
!       6 'E' or 'e' or 'D' or 'd', exponent marker,
!       7 exponent sign,
!       8 exponent integer part,
!       9 exponent decimal point,
!      10 exponent fraction part,
!      11 blanks,
!      12 final comma or semicolon,
!
!    with most quantities optional.
!
!  Examples:
!
!    S                 DVAL
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
!  Modified:
!
!    07 September 2004
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
!    Output, real ( kind = 8 ) DVAL, the value read from the string.
!
!    Output, integer IERROR, error flag.
!    0, no errors occurred.
!    1, 2, 6 or 7, the input number was garbled.  The
!    value of IERROR is the last type of input successfully
!    read.  For instance, 1 means initial blanks, 2 means
!    a plus or minus sign, and so on.
!
!    Output, integer LENGTH, the number of characters read
!    to form the number, including any terminating
!    characters such as a trailing comma or blanks.
!
  implicit none

  logical ch_eqi
  character c
  real ( kind = 8 ) dval
  integer ierror
  integer ihave
  integer isgn
  integer iterm
  integer jbot
  integer jsgn
  integer jtop
  integer length
  integer nchar
  integer ndig
  real ( kind = 8 ) rbot
  real ( kind = 8 ) rexp
  real ( kind = 8 ) rtop
  character ( len = * ) s

  nchar = len_trim ( s )

  ierror = 0
  dval = 0.0D+00
  length = -1
  isgn = 1
  rtop = 0
  rbot = 1
  jsgn = 1
  jtop = 0
  jbot = 1
  ihave = 1
  iterm = 0

  do

    length = length + 1

    if ( nchar < length+1 ) then
      exit
    end if

    c = s(length+1:length+1)
!
!  Blank character.
!
    if ( c == ' ' ) then

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
        isgn = -1
      else if ( ihave == 6 ) then
        ihave = 7
        jsgn = -1
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
!  Scientific notation exponent marker.
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
    else if (  ihave < 11 .and. lle ( '0', c ) .and. lle ( c, '9' ) ) then

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
    if ( iterm == 1 ) then
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
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'S_TO_R8 - Serious error!'
    write ( *, '(a)' ) '  Illegal or nonnumeric input:'
    write ( *, '(a)' ) '    ' // trim ( s )
    return
  end if
!
!  Number seems OK.  Form it.
!
  if ( jtop == 0 ) then
    rexp = 1.0D+00
  else
    if ( jbot == 1 ) then
      rexp = 10.0D+00 ** ( jsgn * jtop )
    else
      rexp = 10.0D+00 ** ( real ( jsgn * jtop, kind = 8 ) &
        / real ( jbot, kind = 8 ) )
    end if
  end if

  dval = real ( isgn, kind = 8 ) * rexp * rtop / rbot

  return
end
subroutine s_to_r8vec ( s, n, rvec, ierror )

!*****************************************************************************80
!
!! S_TO_R8VEC reads an R8VEC from a string.
!
!  Modified:
!
!    07 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string to be read.
!
!    Input, integer N, the number of values expected.
!
!    Output, real ( kind = 8 ) RVEC(N), the values read from the string.
!
!    Output, integer IERROR, error flag.
!    0, no errors occurred.
!    -K, could not read data for entries -K through N.
!
  implicit none

  integer n

  integer i
  integer ierror
  integer ilo
  integer lchar
  real ( kind = 8 ) rvec(n)
  character ( len = * ) s

  i = 0
  ierror = 0
  ilo = 1

  do while ( i < n )

    i = i + 1

    call s_to_r8 ( s(ilo:), rvec(i), ierror, lchar )

    if ( ierror /= 0 ) then
      ierror = -i
      exit
    end if

    ilo = ilo + lchar

  end do

  return
end
subroutine s_word_count ( s, nword )

!*****************************************************************************80
!
!! S_WORD_COUNT counts the number of "words" in a string.
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
!    Output, integer NWORD, the number of "words" in the string.
!    Words are presumed to be separated by one or more blanks.
!
  implicit none

  logical blank
  integer i
  integer lens
  integer nword
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
subroutine sort_heap_external ( n, indx, i, j, isgn )

!*****************************************************************************80
!
!! SORT_HEAP_EXTERNAL externally sorts a list of items into ascending order.
!
!  Discussion:
!
!    The actual list of data is not passed to the routine.  Hence this
!    routine may be used to sort integers, reals, numbers, names,
!    dates, shoe sizes, and so on.  After each call, the routine asks
!    the user to compare or interchange two items, until a special
!    return value signals that the sorting is completed.
!
!  Modified:
!
!    05 February 2004
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
!    Input, integer N, the number of items to be sorted.
!
!    Input/output, integer INDX, the main communication signal.
!
!    The user must set INDX to 0 before the first call.
!    Thereafter, the user should not change the value of INDX until
!    the sorting is done.
!
!    On return, if INDX is
!
!      greater than 0,
!      * interchange items I and J;
!      * call again.
!
!      less than 0,
!      * compare items I and J;
!      * set ISGN = -1 if I < J, ISGN = +1 if J < I;
!      * call again.
!
!      equal to 0, the sorting is done.
!
!    Output, integer I, J, the indices of two items.
!    On return with INDX positive, elements I and J should be interchanged.
!    On return with INDX negative, elements I and J should be compared, and
!    the result reported in ISGN on the next call.
!
!    Input, integer ISGN, results of comparison of elements I and J.
!    (Used only when the previous call returned INDX less than 0).
!    ISGN <= 0 means I is less than or equal to J;
!    0 <= ISGN means I is greater than or equal to J.
!
  implicit none

  integer i
  integer, save :: i_save = 0
  integer indx
  integer isgn
  integer j
  integer, save :: j_save = 0
  integer, save :: k = 0
  integer, save :: k1 = 0
  integer n
  integer, save :: n1 = 0
!
!  INDX = 0: This is the first call.
!
  if ( indx == 0 ) then

    i_save = 0
    j_save = 0
    k = n / 2
    k1 = k
    n1 = n
!
!  INDX < 0: The user is returning the results of a comparison.
!
  else if ( indx < 0 ) then

    if ( indx == -2 ) then

      if ( isgn < 0 ) then
        i_save = i_save + 1
      end if

      j_save = k1
      k1 = i_save
      indx = -1
      i = i_save
      j = j_save
      return

    end if

    if ( 0 < isgn ) then
      indx = 2
      i = i_save
      j = j_save
      return
    end if

    if ( k <= 1 ) then

      if ( n1 == 1 ) then
        i_save = 0
        j_save = 0
        indx = 0
      else
        i_save = n1
        n1 = n1 - 1
        j_save = 1
        indx = 1
      end if

      i = i_save
      j = j_save
      return

    end if

    k = k - 1
    k1 = k
!
!  0 < INDX, the user was asked to make an interchange.
!
  else if ( indx == 1 ) then

    k1 = k

  end if

  do

    i_save = 2 * k1

    if ( i_save == n1 ) then
      j_save = k1
      k1 = i_save
      indx = -1
      i = i_save
      j = j_save
      return
    else if ( i_save <= n1 ) then
      j_save = i_save + 1
      indx = -2
      i = i_save
      j = j_save
      return
    end if

    if ( k <= 1 ) then
      exit
    end if

    k = k - 1
    k1 = k

  end do

  if ( n1 == 1 ) then
    i_save = 0
    j_save = 0
    indx = 0
    i = i_save
    j = j_save
  else
    i_save = n1
    n1 = n1 - 1
    j_save = 1
    indx = 1
    i = i_save
    j = j_save
  end if

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
function triangle_area_2d ( t )

!*****************************************************************************80
!
!! TRIANGLE_AREA_2D computes the area of a triangle in 2D.
!
!  Modified:
!
!    17 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T(2,3), the triangle vertices.
!
!    Output, real ( kind = 8 ) TRIANGLE_AREA_2D, the absolute area of
!    the triangle.
!
  implicit none

  integer, parameter :: dim_num = 2

  real ( kind = 8 ) t(dim_num,3)
  real ( kind = 8 ) triangle_area_2d

  triangle_area_2d = 0.5D+00 * abs ( &
      t(1,1) * ( t(2,2) - t(2,3) ) &
    + t(1,2) * ( t(2,3) - t(2,1) ) &
    + t(1,3) * ( t(2,1) - t(2,2) ) )

  return
end
subroutine triangles3_write ( file_name, element_num, element_node, &
  node_num, node3_label )

!*****************************************************************************80
!
!! TRIANGLES3_WRITE writes the pressure triangles to a file.
!
!  Discussion:
!
!    The first three rows of the array ELEMENT_NODE(6,NODE) contain
!    exactly the nodes that make up the pressure triangles.
!
!    However, we must relabel the nodes!
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
!    Input, character ( len = * ) FILE_NAME, the file name.
!
!    Input, integer ELEMENT_NUM, the number of triangles.
!
!    Input, integer ELEMENT_NODE(6,ELEMENT_NUM), the nodes that
!    make up each triangle.
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, integer NODE3_LABEL(NODE_NUM), contains the renumbered
!    label of order3 nodes, and -1 for nodes that are not order3 nodes.
!
  implicit none

  integer node_num
  integer element_num

  integer element
  integer, dimension(6,element_num) :: element_node
  character ( len = * ) :: file_name
  integer file_status
  integer file_unit
  integer node3_label(node_num)

  call get_unit ( file_unit )

  open ( unit = file_unit, file = file_name, status = 'replace', &
    iostat = file_status )

  if ( file_status /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRIANGLES3_WRITE - Warning!'
    write ( *, '(a)' ) '  Could not write the file "' &
      // trim ( file_name ) // '".'
    return
  end if

  do element = 1, element_num

    write ( file_unit, '(2x,i8,2x,i8,2x,i8)' ) &
      node3_label( element_node(1:3,element) )

  end do

  close ( unit = file_unit )

  return
end
subroutine triangulation_order6_boundary_node ( node_num, element_order, &
  element_num, element_node, node_boundary )

!*****************************************************************************80
!
!! TRIANGULATION_ORDER6_BOUNDARY_NODE indicates which nodes are on the boundary.
!
!  Discussion:
!
!    This routine is given a triangulation, an abstract list of sets of
!    of nodes.  It is assumed that the nodes in each triangle are listed
!    in a counterclockwise order, although the routine should work
!    if the nodes are consistently listed in a clockwise order as well.
!
!    It is assumed that each edge of the triangulation is either
!    * an INTERIOR edge, which is listed twice, once with positive
!      orientation and once with negative orientation, or;
!    * a BOUNDARY edge, which will occur only once.
!
!    This routine should work even if the region has holes - as long
!    as the boundary of the hole comprises more than 3 edges!
!
!  Modified:
!
!    07 June 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, integer ELEMENT_ORDER, the order of the elements.
!
!    Input, integer ELEMENT_NUM, the number of elements.
!
!    Input, integer ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM),
!    the nodes that make up the elements.  These should be listed
!    in counterclockwise order.
!
!    Output, logical NODE_BOUNDARY(NODE_NUM), is TRUE if the node
!    is on a boundary edge.
!
  implicit none

  integer node_num
  integer element_order
  integer element_num

  integer e1(6*element_num)
  integer e2(6*element_num)
  integer edge(2,6*element_num)
  integer i
  integer m
  integer n
  logical node_boundary(node_num)
  integer element_node(element_order,element_num)

  m = 2
  n = 6 * element_num
!
!  Set up the edge array.
!
  edge(1,                 1:  element_num) = element_node(1,1:element_num)
  edge(2,                 1:  element_num) = element_node(4,1:element_num)
  edge(1,    element_num+1:2*element_num) = element_node(4,1:element_num)
  edge(2,    element_num+1:2*element_num) = element_node(2,1:element_num)
  edge(1,  2*element_num+1:3*element_num) = element_node(2,1:element_num)
  edge(2,  2*element_num+1:3*element_num) = element_node(5,1:element_num)
  edge(1,  3*element_num+1:4*element_num) = element_node(5,1:element_num)
  edge(2,  3*element_num+1:4*element_num) = element_node(3,1:element_num)
  edge(1,  4*element_num+1:5*element_num) = element_node(3,1:element_num)
  edge(2,  4*element_num+1:5*element_num) = element_node(6,1:element_num)
  edge(1,  5*element_num+1:6*element_num) = element_node(6,1:element_num)
  edge(2,  5*element_num+1:6*element_num) = element_node(1,1:element_num)
!
!  In each column, force the smaller entry to appear first.
!
  e1(1:n) = minval ( edge(1:2,1:n), dim = 1 )
  e2(1:n) = maxval ( edge(1:2,1:n), dim = 1 )

  edge(1,1:n) = e1(1:n)
  edge(2,1:n) = e2(1:n)
!
!  Ascending sort the column array.
!
  call i4col_sort_a ( m, n, edge )
!
!  Records which appear twice are internal edges and can be ignored.
!
  node_boundary(1:node_num) = .false.

  i = 0

  do while ( i < n )

    i = i + 1

    if ( i == n ) then
      node_boundary(edge(1:m,i)) = .true.
    else if ( all ( edge(1:m,i) == edge(1:m,i+1) ) ) then
      i = i + 1
    else
      node_boundary(edge(1:m,i)) = .true.
    end if

  end do

  return
end
subroutine triangulation_order6_plot ( file_name, node_num, node_xy, &
  element_num, element_node, node_show, triangle_show )

!*****************************************************************************80
!
!! TRIANGULATION_ORDER6_PLOT plots a 6-node triangulation of a set of nodes.
!
!  Discussion:
!
!    The triangulation is most usually a Delaunay triangulation,
!    but this is not necessary.
!
!    In a six node triangulation, it is assumed that nodes 1, 2, and 3
!    are the vertices of the triangles, and that nodes 4, 5, and 6
!    lie between 1 and 2, 2 and 3, and 3 and 1 respectively.
!
!  Modified:
!
!    27 September 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the output file.
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, real ( kind = 8 ) NODE_XY(2,NODE_NUM), the coordinates of the nodes.
!
!    Input, integer ELEMENT_NUM, the number of elements.
!
!    Input, integer ELEMENT_NODE(6,ELEMENT_NUM), lists, for each element,
!    the indices of the nodes that form the vertices of the element.
!
!    Input, integer NODE_SHOW,
!    0, do not show nodes;
!    1, show nodes;
!    2, show nodes and label them.
!
!    Input, integer TRIANGLE_SHOW,
!    0, do not show triangles;
!    1, show triangles;
!    2, show triangles and label them.
!
  implicit none

  integer node_num
  integer element_num

  real ( kind = 8 ) ave_x
  real ( kind = 8 ) ave_y
  integer :: circle_size
  integer delta
  integer element
  integer element_node(6,element_num)
  character ( len = * ) file_name
  integer file_unit
  integer i
  integer ios
  integer node
  integer node_show
  real ( kind = 8 ) node_xy(2,node_num)
  character ( len = 40 ) string
  integer triangle_show
  real ( kind = 8 ) x_max
  real ( kind = 8 ) x_min
  integer x_ps
  integer :: x_ps_max = 576
  integer :: x_ps_max_clip = 594
  integer :: x_ps_min = 36
  integer :: x_ps_min_clip = 18
  real ( kind = 8 ) x_scale
  real ( kind = 8 ) y_max
  real ( kind = 8 ) y_min
  integer y_ps
  integer :: y_ps_max = 666
  integer :: y_ps_max_clip = 684
  integer :: y_ps_min = 126
  integer :: y_ps_min_clip = 108
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

    y_ps_max      = y_ps_max - delta
    y_ps_min      = y_ps_min + delta

    y_ps_max_clip = y_ps_max_clip - delta
    y_ps_min_clip = y_ps_min_clip + delta

    y_scale = x_scale

  end if
!
!  Open the file.
!
  call get_unit ( file_unit )

  open ( unit = file_unit, file = file_name, status = 'replace', &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRIANGULATION_ORDER6_PLOT - Fatal error!'
    write ( *, '(a)' ) '  Can not open output file.'
    return
  end if
!
!  Write the header.
!
  write ( file_unit, '(a)' ) '%!PS-Adobe-3.0 EPSF-3.0'
  write ( file_unit, '(a)' ) '%%Creator: triangulation_order6_plot.f90'
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
  write ( file_unit, '(a)' ) '%  (Triangulation)  show'
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
  if ( node_num <= 200 ) then
    circle_size = 5
  else if ( node_num <= 500 ) then
    circle_size = 4
  else if ( node_num <= 1000 ) then
    circle_size = 3
  else if ( node_num <= 5000 ) then
    circle_size = 2
  else
    circle_size = 1
  end if

  if ( 1 <= node_show ) then
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '%  Draw filled dots at the nodes.'
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

  end if
!
!  Label the nodes.
!
  if ( 2 <= node_show ) then

    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '%  Label the nodes:'
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '%  Set the RGB color to darker blue.'
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '0.100  0.250  0.850 setrgbcolor'
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
!
!  Draw the triangles.
!
  if ( 1 <= triangle_show ) then
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '%  Set the RGB color to red.'
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '0.900  0.200  0.100 setrgbcolor'
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '%  Draw the triangles.'
    write ( file_unit, '(a)' ) '%'

    do element = 1, element_num

      write ( file_unit, '(a)' ) 'newpath'

      node = element_node(6,element)

      x_ps = int ( &
        ( ( x_max - node_xy(1,node)         ) * real ( x_ps_min, kind = 8 )   &
        + (         node_xy(1,node) - x_min ) * real ( x_ps_max, kind = 8 ) ) &
        / ( x_max                   - x_min ) )

      y_ps = int ( &
        ( ( y_max - node_xy(2,node)         ) * real ( y_ps_min, kind = 8 )   &
        + (         node_xy(2,node) - y_min ) * real ( y_ps_max, kind = 8 ) ) &
        / ( y_max                   - y_min ) )

      write ( file_unit, '(i3,2x,i3,2x,a)' ) x_ps, y_ps, ' moveto'

      do i = 1, 3

        node = element_node(i,element)

        x_ps = int ( &
          ( ( x_max - node_xy(1,node)         ) * real ( x_ps_min, kind = 8 )   &
          + (         node_xy(1,node) - x_min ) * real ( x_ps_max, kind = 8 ) ) &
          / ( x_max                   - x_min ) )

        y_ps = int ( &
          ( ( y_max - node_xy(2,node)         ) * real ( y_ps_min, kind = 8 )   &
          + (         node_xy(2,node) - y_min ) * real ( y_ps_max, kind = 8 ) ) &
          / ( y_max                   - y_min ) )

        write ( file_unit, '(i3,2x,i3,2x,a)' ) x_ps, y_ps, ' lineto'

        node = element_node(i+3,element)

        x_ps = int ( &
          ( ( x_max - node_xy(1,node)         ) * real ( x_ps_min, kind = 8 )   &
          + (         node_xy(1,node) - x_min ) * real ( x_ps_max, kind = 8 ) ) &
          / ( x_max                   - x_min ) )

        y_ps = int ( &
          ( ( y_max - node_xy(2,node)         ) * real ( y_ps_min, kind = 8 )   &
          + (         node_xy(2,node) - y_min ) * real ( y_ps_max, kind = 8 ) ) &
          / ( y_max                   - y_min ) )

        write ( file_unit, '(i3,2x,i3,2x,a)' ) x_ps, y_ps, ' lineto'

      end do

      write ( file_unit, '(a)' ) 'stroke'

    end do

  end if
!
!  Label the triangles.
!
  if ( 2 <= triangle_show ) then

    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '%  Label the triangles:'
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '%  Set the RGB color to darker red.'
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '0.950  0.250  0.150 setrgbcolor'
    write ( file_unit, '(a)' ) '/Times-Roman findfont'
    write ( file_unit, '(a)' ) '0.20 inch scalefont'
    write ( file_unit, '(a)' ) 'setfont'

    do element = 1, element_num

      ave_x = 0.0D+00
      ave_y = 0.0D+00

      do i = 1, 6

        node = element_node(i,element)

        ave_x = ave_x + node_xy(1,node)
        ave_y = ave_y + node_xy(2,node)

      end do

      ave_x = ave_x / 6.0D+00
      ave_y = ave_y / 6.0D+00

      x_ps = int ( &
        ( ( x_max - ave_x         ) * real ( x_ps_min, kind = 8 )   &
        + (       + ave_x - x_min ) * real ( x_ps_max, kind = 8 ) ) &
        / ( x_max         - x_min ) )

      y_ps = int ( &
        ( ( y_max - ave_y         ) * real ( y_ps_min, kind = 8 )   &
        + (         ave_y - y_min ) * real ( y_ps_max, kind = 8 ) ) &
        / ( y_max         - y_min ) )

      write ( string, '(i4)' ) element
      string = adjustl ( string )

      write ( file_unit, '(i4,2x,i4,a)' ) x_ps, y_ps, ' moveto (' &
        // trim ( string ) // ') show'

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
subroutine velocity6_write ( file_name, node_num, node_u_variable, &
  node_v_variable, variable_num, node_c )

!*****************************************************************************80
!
!! VELOCITY6_WRITE writes the velocities to a file.
!
!  Modified:
!
!    30 July 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the file name.
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, integer NODE_U_VARIABLE(NODE_NUM), NODE_V_VARIABLE(NODE_NUM),
!    the indices of the horizontal and vertical velocity variables
!    associated with the node, or -1 if there is none.
!
!    Input, integer VARIABLE_NUM, the number of variables.
!
!    Input, real ( kind = 8 ) NODE_C(VARIABLE_NUM), the finite element
!    coefficients.
!
  implicit none

  integer node_num
  integer variable_num

  character ( len = * ) :: file_name
  integer file_status
  integer file_unit
  integer node
  real ( kind = 8 ), dimension(variable_num) :: node_c
  integer, dimension(node_num) :: node_u_variable
  integer, dimension(node_num) :: node_v_variable
  real ( kind = 8 ) u
  integer u_index
  real ( kind = 8 ) v
  integer v_index
  integer variable

  call get_unit ( file_unit )

  open ( unit = file_unit, file = file_name, status = 'replace', &
    iostat = file_status )

  if ( file_status /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'VELOCITY6_WRITE - Warning!'
    write ( *, '(a)' ) '  Could not write the file "' &
      // trim ( file_name ) // '".'
    return
  end if

  do node = 1, node_num

    u_index = node_u_variable(node)

    if ( 0 < u_index ) then
      u = node_c(u_index)
    else
      u = 0.0D+00
    end if

    v_index = node_v_variable(node)

    if ( 0 < v_index ) then
      v = node_c(v_index)
    else
      v = 0.0D+00
    end if

    write ( file_unit, '(2x,g14.6,2x,g14.6)' ) u, v

  end do

  close ( unit = file_unit )

  return
end
