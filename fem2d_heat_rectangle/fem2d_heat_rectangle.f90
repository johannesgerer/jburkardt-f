program main

!*****************************************************************************80
!
!! MAIN is the main program for FEM2D_HEAT_RECTANGLE.
!
!  Discussion:
!
!    FEM2D_HEAT_RECTANGLE solves
!
!      dUdT - Laplacian U(X,Y,T) = F(X,Y,T)
!
!    in a rectangular region in the plane.
!
!    Along the boundary of the region, Dirichlet conditions
!    are imposed:
!
!      U(X,Y,T) = G(X,Y,T)
!
!    At the initial time T_INIT, the value of U is given at all points
!    in the region:
!
!      U(X,Y,T) = H(X,Y,T)
!
!    The code uses continuous piecewise quadratic basis functions on
!    triangles determined by a uniform grid of NX by NY points.
!
!    The backward Euler approximation is used for the time derivatives.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 December 2010
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!  Local parameters:
!
!    Local, real ( kind = 8 ) A(3*IB+1,NODE_NUM), the coefficient matrix.
!
!    Local, real ( kind = 8 ) EH1, the H1 seminorm error.
!
!    Local, real ( kind = 8 ) EL2, the L2 error.
!
!    Local, real ( kind = 8 ) ELEMENT_AREA(ELEMENT_NUM), the area of elements.
!
!    Local, integer ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM);
!    ELEMENT_NODE(I,J) is the global index of local node I in element J.
!
!    Local, integer ELEMENT_NUM, the number of elements.
!
!    Local, real ( kind = 8 ) F(NODE_NUM), the right hand side.
!
!    Local, integer IB, the half-bandwidth of the matrix.
!
!    Local, integer NODE_BOUNDARY(NODE_NUM), is
!    0, if a node is an interior node;
!    1, if a node is a Dirichlet boundary node.
!
!    Local, integer ELEMENT_ORDER, the number of nodes used to form one element.
!
!    Local, real ( kind = 8 ) NODE_XY(2,NODE_NUM), the
!    coordinates of nodes.
!
!    Local, integer NX, the number of points in the X direction.
!
!    Local, integer NY, the number of points in the Y direction.
!
!    Local, integer QUAD_NUM, the number of quadrature points used for assembly.
!
!    Local, real ( kind = 8 ) TIME, the current time.
!
!    Local, real ( kind = 8 ) TIME_FINAL, the final time.
!
!    Local, real ( kind = 8 ) TIME_INIT, the initial time.
!
!    Local, real ( kind = 8 ) TIME_OLD, the time at the previous time step.
!
!    Local, integer TIME_STEP_NUM, the number of time steps to take.
!
!    Local, real ( kind = 8 ) TIME_STEP_SIZE, the size of the time steps.
!
!    Local, real ( kind = 8 ) U(NODE_NUM), the finite element coefficients
!    defining the solution at the current time.
!
!    Local, real ( kind = 8 ) WQ(QUAD_NUM), quadrature weights.
!
!    Local, real ( kind = 8 ) XL, XR, YB, YT, the X coordinates of
!    the left and right sides of the rectangle, and the Y coordinates
!    of the bottom and top of the rectangle.
!
!    Local, real ( kind = 8 ) XQ(QUAD_NUM,ELEMENT_NUM),
!    YQ(QUAD_NUM,ELEMENT_NUM), the
!    coordinates of the quadrature points in each element.
!
  implicit none

  integer, parameter :: element_order = 6
  integer, parameter :: quad_num = 3
  integer, parameter :: nx = 5
  integer, parameter :: ny = 5

  integer, parameter :: element_num = ( nx - 1 ) * ( ny - 1 ) * 2
  integer, parameter :: node_num = ( 2 * nx - 1 ) * ( 2 * ny - 1 )

  real ( kind = 8 ), allocatable, dimension (:,:) :: a
  logical, parameter :: debug = .false.
  real ( kind = 8 ) dt
  real ( kind = 8 ), allocatable, dimension ( : ) :: dudx_exact
  real ( kind = 8 ), allocatable, dimension ( : ) :: dudy_exact
  real ( kind = 8 ) eh1
  real ( kind = 8 ) el2
  real ( kind = 8 ), dimension(element_num) :: element_area
  integer, dimension(element_order,element_num) :: element_node
  real ( kind = 8 ), allocatable, dimension (:) :: f
  integer i
  integer ib
  integer ierr
  integer job
  integer node
  integer, dimension(node_num) :: node_boundary
  character ( len = 255 ) node_eps_file_name
  character ( len = 255 ) node_txt_file_name
  logical node_label
  integer node_show
  real ( kind = 8 ), dimension(2,node_num) :: node_xy
  integer, allocatable, dimension (:) :: pivot
  real ( kind = 8 ) time
  character ( len = 255 ) time_file_name
  real ( kind = 8 ) time_final
  real ( kind = 8 ) time_init
  real ( kind = 8 ) time_old
  integer time_step
  integer time_step_num
  real ( kind = 8 ) time_step_size
  integer time_unit
  integer triangle_show
  character ( len = 255 ) triangulation_eps_file_name
  character ( len = 255 ) triangulation_txt_file_name
  real ( kind = 8 ), allocatable, dimension ( : ) :: u
  real ( kind = 8 ), allocatable, dimension ( : ) :: u_exact
  character ( len = 255 ) u_file_name
  real ( kind = 8 ), allocatable, dimension ( : ) :: u_old
  real ( kind = 8 ), dimension(quad_num) :: wq
  real ( kind = 8 ), parameter :: xl = 0.0D+00
  real ( kind = 8 ), dimension(quad_num,element_num) :: xq
  real ( kind = 8 ), parameter :: xr = 1.0D+00
  real ( kind = 8 ), parameter :: yb = 0.0D+00
  real ( kind = 8 ), dimension(quad_num,element_num) :: yq
  real ( kind = 8 ), parameter :: yt = 1.0D+00

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FEM2D_HEAT_RECTANGLE'
  write ( *, '(a)' ) '  FORTRAN90 version:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Solution of the time dependent heat equation'
  write ( *, '(a)' ) '  on a unit box in 2 dimensions.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Ut - Uxx - Uyy = F(x,y,t) in the box'
  write ( *, '(a)' ) '        U(x,y,t) = G(x,y,t) for (x,y) on the boundary.'
  write ( *, '(a)' ) '        U(x,y,t) = H(x,y,t) for t = T_INIT.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The finite element method is used, with piecewise'
  write ( *, '(a)' ) '  quadratic basis functions on 6 node triangular'
  write ( *, '(a)' ) '  elements.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The backward Euler formula is used for '
  write ( *, '(a)' ) '  the time derivative.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The corner nodes of the triangles are generated by an'
  write ( *, '(a)' ) '  underlying grid whose dimensions are'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  NX =                         ', nx
  write ( *, '(a,i8)' ) '  NY =                         ', ny
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of nodes =            ', node_num
  write ( *, '(a,i8)' ) '  Number of elements =         ', element_num
!
!  Set the coordinates of the nodes.
!
  call xy_set ( nx, ny, node_num, xl, xr, yb, yt, node_xy )
!
!  Organize the nodes into a grid of 6-node triangles.
!
  call grid_t6 ( nx, ny, element_order, element_num, element_node )
!
!  Set the quadrature rule for assembly.
!
  call quad_a ( node_xy, element_node, element_num, node_num, &
    element_order, wq, xq, yq )
!
!  Determine the areas of the elements.
!
  call area_set ( node_num, node_xy, element_order, element_num, &
    element_node, element_area )
!
!  Determine which nodes have an associated finite element unknown.
!
  call node_boundary_set ( nx, ny, node_num, node_boundary )

  if ( .false. ) then
    call i4vec_print_some ( node_num, node_boundary, 10, '  NODE_BOUNDARY:' )
  end if
!
!  Determine the bandwidth of the coefficient matrix.
!
  call bandwidth ( element_order, element_num, element_node, node_num, ib )

  write ( *, '(a,i8)' ) '  The matrix half bandwidth is ', ib
  write ( *, '(a,i8)' ) '  The matrix row size is       ', 3 * ib + 1
!
!  Make an EPS picture of the nodes.
!
  if ( nx <= 10 .and. ny <= 10 ) then

    node_eps_file_name = 'nodes.eps'
    node_label = .true.

    call nodes_plot ( node_eps_file_name, node_num, node_xy, node_label )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FEM2D_HEAT_RECTANGLE:'
    write ( *, '(a)' ) '  Wrote an EPS file'
    write ( *, '(a)' ) '    "' // trim ( node_eps_file_name ) // '"'
    write ( *, '(a)' ) '  containing a picture of the nodes.'

  end if
!
!  Write the nodes to an ASCII file that can be read into MATLAB.
!
  node_txt_file_name = 'nodes.txt'
  call nodes_write ( node_num, node_xy, node_txt_file_name )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FEM2D_HEAT_RECTANGLE:'
  write ( *, '(a)' ) '  Wrote an ASCII node file'
  write ( *, '(a)' ) '    "' // trim ( node_txt_file_name ) // '"'
  write ( *, '(a)' ) '  of the form'
  write ( *, '(a)' ) '    X(I), Y(I)'
  write ( *, '(a)' ) '  which can be used for plotting.'
!
!  Make an EPS picture of the elements.
!
  if ( nx <= 10 .and. ny <= 10 ) then

    triangulation_eps_file_name = 'elements.eps'

    node_show = 2
    triangle_show = 2

    call triangulation_order6_plot ( triangulation_eps_file_name, node_num, &
      node_xy, element_num, element_node, node_show, triangle_show )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FEM2D_HEAT_RECTANGLE:'
    write ( *, '(a)' ) '  Wrote an EPS file'
    write ( *, '(a)' ) '    "' // trim ( triangulation_eps_file_name ) // '"'
    write ( *, '(a)' ) '  containing a picture of the elements.'

  end if
!
!  Write the elements to a file that can be read into MATLAB.
!
  triangulation_txt_file_name = 'elements.txt'

  call element_write ( element_order, element_num, element_node, &
    triangulation_txt_file_name )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FEM2D_HEAT_RECTANGLE:'
  write ( *, '(a)' ) '  Wrote an ASCII element file'
  write ( *, '(a)' ) '    "' // trim ( triangulation_txt_file_name ) // '"'
  write ( *, '(a)' ) '  of the form'
  write ( *, '(a)' ) '    Node(1) Node(2) Node(3) Node(4) Node(5) Node(6)'
  write ( *, '(a)' ) '  which can be used for plotting.'
!
!  Set time stepping quantities.
!
  time_init = 0.0D+00
  time_final = 0.5D+00
  time_step_num = 10
  time_step_size = ( time_final - time_init ) / real ( time_step_num, kind = 8 )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Initial time = ', time_init
  write ( *, '(a,g14.6)' ) '  Final time =   ', time_final
  write ( *, '(a,g14.6)' ) '  Step size =    ', time_step_size
  write ( *, '(a,i8)' ) '  Number of steps = ', time_step_num
!
!  Allocate space for the coefficient matrix A and right hand side F.
!
  allocate ( a(3*ib+1,node_num) )
  allocate ( dudx_exact(node_num) )
  allocate ( dudy_exact(node_num) )
  allocate ( f(node_num) )
  allocate ( pivot(node_num) )
  allocate ( u(node_num) )
  allocate ( u_exact(node_num) )
  allocate ( u_old(node_num) )
!
!  Initialize the names of the time and solution file.
!
  time_file_name = 'time.txt'
  u_file_name = 'u0000.txt'
!
!  Set the value of U at the initial time.
!
  time = time_init
  call exact_u ( node_num, node_xy, time, u_exact, dudx_exact, dudy_exact )
  u(1:node_num) = u_exact(1:node_num)

  call get_unit ( time_unit )
  open ( unit = time_unit, file = time_file_name, status = 'replace' )
  write ( time_unit, '(g14.6)' ) time
  call solution_write ( node_num, u, u_file_name )
!
!  Time looping.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     Time        L2 Error       H1 Error'
  write ( *, '(a)' ) ' '

  do time_step = 1, time_step_num

    time_old = time
    u_old(1:node_num) = u(1:node_num)

    time = ( real ( time_step_num - time_step, kind = 8 ) * time_init    &
           + real (                 time_step, kind = 8 ) * time_final ) &
           / real ( time_step_num,             kind = 8 )
!
!  Assemble the finite element contributions to the coefficient matrix A
!  and the right-hand side F.
!
    call assemble ( node_num, node_xy, element_order, element_num, &
      element_node, quad_num, wq, xq, yq, element_area, ib, time, a, f )

    if ( debug ) then

      call dgb_print_some ( node_num, node_num, ib, ib, a, 10, 1, 12, 25, &
        '  Initial block of matrix A:' )

      call r8vec_print_some ( node_num, f, node_num, &
        '  Part of right hand side F:' )

    end if
!
!  Modify the coefficient matrix and right hand side to account for the dU/dt
!  term, which we are treating using the backward Euler formula.
!
    call adjust_backward_euler ( node_num, node_xy, element_order, &
      element_num, element_node, quad_num, wq, xq, yq, element_area, ib, time, &
      time_step_size, u_old, a, f )

    if ( debug ) then

      call dgb_print_some ( node_num, node_num, ib, ib, a, 10, 1, 12, 25, &
        '  A after DT adjustment:' )

      call r8vec_print_some ( node_num, f, node_num, &
        '  F after DT adjustment:' )

    end if
!
!  Modify the coefficient matrix and right hand side to account for
!  boundary conditions.
!
    call adjust_boundary ( node_num, node_xy, node_boundary, &
      ib, time, a, f )

    if ( debug ) then

      call dgb_print_some ( node_num, node_num, ib, ib, a, 10, 1, 12, 25, &
        '  A after BC adjustment:' )

      call r8vec_print_some ( node_num, f, node_num, &
        '  F after BC adjustment:' )

    end if
!
!  Factor and solve the linear system A * U = F.
!
    call dgb_fa ( node_num, ib, ib, a, pivot, ierr )

    if ( ierr /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FEM2D_HEAT_RECTANGLE - Fatal error!'
      write ( *, '(a)' ) '  DGB_FA returned an error condition.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The linear system was not factored, and the'
      write ( *, '(a)' ) '  algorithm cannot proceed.'
      stop
    end if

    job = 0
    call dgb_sl ( node_num, ib, ib, a, pivot, f, job )

    u(1:node_num) = f(1:node_num)

    if ( debug ) then
      call r8vec_print_some ( node_num, u, node_num, &
        '  Part of the solution vector U:' )
    end if
!
!  Calculate EL2(TIME) and EH1, the L2 and H1 errors at the current time.
!
    call errors ( element_area, element_node, node_xy, u, &
      element_num, element_order, node_num, time, el2, eh1 )
!
!  Compare the approximate and exact solutions at the nodes.
!
    if ( .false. ) then
      call compare ( node_num, node_xy, time, u )
    end if
!
!  Increment the file name, and write the new solution.
!
    write ( time_unit, '(g14.6)' ) time

    call file_name_inc ( u_file_name )

    call solution_write ( node_num, u, u_file_name )

  end do

  close ( unit = time_unit )
!
!  Deallocate memory.
!
  deallocate ( a )
  deallocate ( dudx_exact )
  deallocate ( dudy_exact )
  deallocate ( f )
  deallocate ( pivot )
  deallocate ( u )
  deallocate ( u_exact )
  deallocate ( u_old )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FEM2D_HEAT_RECTANGLE:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine adjust_backward_euler ( node_num, node_xy, element_order, &
  element_num, element_node, quad_num, wq, xq, yq, element_area, ib, time, &
  time_step_size, u_old, a, f )

!*****************************************************************************80
!
!! ADJUST_BACKWARD_EULER adjusts the system for the backward Euler term.
!
!  Discussion:
!
!    The input linear system
!
!      A * U = F
!
!    is appropriate for the equation
!
!      -Uxx - Uyy = RHS
!
!    We need to modify the matrix A and the right hand side F to
!    account for the approximation of the time derivative in
!
!      Ut - Uxx - Uyy = RHS
!
!    by the backward Euler approximation:
!
!      Ut approximately equal to ( U - Uold ) / dT
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 April 2006
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
!    Input, integer ELEMENT_ORDER, the number of nodes used to form one element.
!
!    Input, integer ELEMENT_NUM, the number of elements.
!
!    Input, integer ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM);
!    ELEMENT_NODE(I,J) is the global index of local node I in element J.
!
!    Input, integer QUAD_NUM, the number of quadrature points used in assembly.
!
!    Input, real ( kind = 8 ) WQ(QUAD_NUM), quadrature weights.
!
!    Input, real ( kind = 8 ) XQ(QUAD_NUM,ELEMENT_NUM),
!    YQ(QUAD_NUM,ELEMENT_NUM), the
!    coordinates of the quadrature points in each element.
!
!    Input, real ( kind = 8 ) ELEMENT_AREA(ELEMENT_NUM), the area of elements.
!
!    Input, integer IB, the half-bandwidth of the matrix.
!
!    Input, real ( kind = 8 ) TIME, the current time.
!
!    Input, real ( kind = 8 ) TIME_STEP_SIZE, the size of the time step.
!
!    Input, real ( kind = 8 ) U_OLD(NODE_NUM), the finite element
!    coefficients for the solution at the previous time.
!
!    Input/output, real ( kind = 8 ) A(3*IB+1,NODE_NUM), the NODE_NUM
!    by NODE_NUM coefficient matrix, stored in a compressed format.
!
!    Input/output, real ( kind = 8 ) F(NODE_NUM), the right hand side.
!
  implicit none

  integer ib
  integer element_num
  integer element_order
  integer node_num
  integer quad_num

  real ( kind = 8 ), dimension(3*ib+1,node_num) :: a
  integer basis
  real ( kind = 8 ) bi
  real ( kind = 8 ) bj
  real ( kind = 8 ) dbidx
  real ( kind = 8 ) dbidy
  real ( kind = 8 ) dbjdx
  real ( kind = 8 ) dbjdy
  integer element
  real ( kind = 8 ), dimension(element_num) :: element_area
  integer, dimension(element_order,element_num) :: element_node
  real ( kind = 8 ), dimension(node_num) :: f
  integer i
  integer ip
  integer ipp
  integer j
  integer node
  real ( kind = 8 ), dimension(2,node_num) :: node_xy
  integer quad
  integer test
  real ( kind = 8 ) time
  real ( kind = 8 ) time_step_size
  real ( kind = 8 ) u_old(node_num)
  real ( kind = 8 ) w
  real ( kind = 8 ), dimension(quad_num) :: wq
  real ( kind = 8 ) x
  real ( kind = 8 ), dimension(quad_num,element_num) :: xq
  real ( kind = 8 ) y
  real ( kind = 8 ), dimension(quad_num,element_num) :: yq

  do element = 1, element_num

    do quad = 1, quad_num

      x = xq(quad,element)
      y = yq(quad,element)
      w = element_area(element) * wq(quad)

      do test = 1, element_order

        node = element_node(test,element)

        call qbf ( x, y, element, test, node_xy, element_node, &
          element_num, element_order, node_num, bi, dbidx, dbidy )
!
!  Carry the U_OLD term to the right hand side.
!
        f(node) = f(node) + w * bi * u_old(node) / time_step_size
!
!  Modify the diagonal entries of A.
!
        do basis = 1, element_order

          j = element_node(basis,element)

          call qbf ( x, y, element, basis, node_xy, element_node, &
            element_num, element_order, node_num, bj, dbjdx, dbjdy )

          a(node-j+2*ib+1,j) = a(node-j+2*ib+1,j) &
            + w * bi * bj / time_step_size

        end do
      end do

    end do

  end do

  return
end
subroutine adjust_boundary ( node_num, node_xy, node_boundary, &
  ib, time, a, f )

!*****************************************************************************80
!
!! ADJUST_BOUNDARY modifies the linear system for boundary conditions.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 April 2006
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
!    Input, integer NODE_BOUNDARY(NODE_NUM), is
!    0, if a node is an interior node;
!    1, if a node is a Dirichlet boundary node.
!
!    Input, integer IB, the half-bandwidth of the matrix.
!
!    Input, real ( kind = 8 ) TIME, the current time.
!
!    Input/output, real ( kind = 8 ) A(3*IB+1,NODE_NUM), the NODE_NUM by
!    NODE_NUM coefficient matrix, stored in a compressed format.  On output,
!    A has been adjusted for boundary conditions.
!
!    Input/output, real ( kind = 8 ) F(NODE_NUM), the right hand side.
!    On output, F has been adjusted for boundary conditions.
!
  implicit none

  integer ib
  integer node_num

  real ( kind = 8 ), dimension(3*ib+1,node_num) :: a
  real ( kind = 8 ) dudx_exact(node_num)
  real ( kind = 8 ) dudy_exact(node_num)
  real ( kind = 8 ), dimension(node_num) :: f
  integer i
  integer j
  integer jhi
  integer jlo
  integer node
  integer, dimension(node_num) :: node_boundary
  real ( kind = 8 ), dimension(2,node_num) :: node_xy
  real ( kind = 8 ) time
  real ( kind = 8 ) u_exact(node_num)
!
!  Get the exact solution at every node.
!
  call exact_u ( node_num, node_xy, time, u_exact, dudx_exact, dudy_exact )

  do node = 1, node_num

    if ( node_boundary(node) /= 0 ) then

      jlo = max ( node - ib, 1 )
      jhi = min ( node + ib, node_num )

      do j = jlo, jhi
        a(node-j+2*ib+1,j) = 0.0D+00
      end do

      a(node-node+2*ib+1,node) = 1.0D+00

      f(node) = u_exact(node)

    end if

  end do

  return
end
subroutine area_set ( node_num, node_xy, element_order, element_num, &
  element_node, element_area )

!*****************************************************************************80
!
!! AREA_SET sets the area of each element.
!
!  Discussion:
!
!    The areas of the elements are needed in order to adjust
!    the integral estimates produced by the quadrature formulas.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 April 2006
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!  Parameters:
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, real ( kind = 8 ) NODE_XY(2,NODE_NUM), the
!    coordinates of the nodes.
!
!    Input, integer ELEMENT_ORDER, the number of local nodes per element.
!
!    Input, integer ELEMENT_NUM, the number of elements.
!
!    Input, integer ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM);
!    ELEMENT_NODE(I,J) is the global index of local node I in element J.
!
!    Output, real ( kind = 8 ) ELEMENT_AREA(ELEMENT_NUM), the area of elements.
!
  implicit none

  integer element_num
  integer element_order
  integer node_num

  integer element
  real ( kind = 8 ) element_area(element_num)
  integer element_node(element_order,element_num)
  integer i1
  integer i2
  integer i3
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ) node_xy(2,node_num)
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) y3

  do element = 1, element_num

    i1 = element_node(1,element)
    x1 = node_xy(1,i1)
    y1 = node_xy(2,i1)

    i2 = element_node(2,element)
    x2 = node_xy(1,i2)
    y2 = node_xy(2,i2)

    i3 = element_node(3,element)
    x3 = node_xy(1,i3)
    y3 = node_xy(2,i3)

    element_area(element) = 0.5D+00 * abs &
      ( y1 * ( x2 - x3 ) &
      + y2 * ( x3 - x1 ) &
      + y3 * ( x1 - x2 ) )

  end do

  return
end
subroutine assemble ( node_num, node_xy, element_order, element_num, &
  element_node, quad_num, wq, xq, yq, element_area, ib, time, a, f )

!*****************************************************************************80
!
!! ASSEMBLE assembles the coefficient matrix A and right hand side F.
!
!  Discussion:
!
!    The matrix is known to be banded.  A special matrix storage format
!    is used to reduce the space required.  Details of this format are
!    discussed in the routine DGB_FA.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 April 2006
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!  Parameters:
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, real ( kind = 8 ) NODE_XY(2,NODE_NUM), the coordinates of nodes.
!
!    Input, integer ELEMENT_ORDER, the number of nodes used to form one element.
!
!    Input, integer ELEMENT_NUM, the number of elements.
!
!    Input, integer ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM);
!    ELEMENT_NODE(I,J) is the global index of local node I in element J.
!
!    Input, integer QUAD_NUM, the number of quadrature points used in assembly.
!
!    Input, real ( kind = 8 ) WQ(QUAD_NUM), quadrature weights.
!
!    Input, real ( kind = 8 ) XQ(QUAD_NUM,ELEMENT_NUM),
!    YQ(QUAD_NUM,ELEMENT_NUM), the
!    coordinates of the quadrature points in each element.
!
!    Input, real ( kind = 8 ) ELEMENT_AREA(ELEMENT_NUM), the area of elements.
!
!    Input, integer IB, the half-bandwidth of the matrix.
!
!    Input, real ( kind = 8 ) TIME, the current time.
!
!    Output, real ( kind = 8 ) A(3*IB+1,NODE_NUM), the NODE_NUM by NODE_NUM
!    coefficient matrix, stored in a compressed format.
!
!    Output, real ( kind = 8 ) F(NODE_NUM), the right hand side.
!
!  Local parameters:
!
!    Local, real ( kind = 8 ) BB, BX, BY, the value of some basis function
!    and its first derivatives at a quadrature point.
!
!    Local, real ( kind = 8 ) BBB, BBX, BBY, the value of another basis
!    function and its first derivatives at a quadrature point.
!
  implicit none

  integer ib
  integer element_num
  integer element_order
  integer node_num
  integer quad_num

  real ( kind = 8 ), dimension(3*ib+1,node_num) :: a
  real ( kind = 8 ) aij
  integer basis
  real ( kind = 8 ) bi
  real ( kind = 8 ) bj
  real ( kind = 8 ) dbidx
  real ( kind = 8 ) dbidy
  real ( kind = 8 ) dbjdx
  real ( kind = 8 ) dbjdy
  integer element
  real ( kind = 8 ), dimension(element_num) :: element_area
  integer, dimension(element_order,element_num) :: element_node
  real ( kind = 8 ), dimension(node_num) :: f
  integer i
  integer ip
  integer ipp
  integer j
  integer node
  real ( kind = 8 ), dimension(2,node_num) :: node_xy
  integer quad
  real ( kind = 8 ) :: rhs
  integer test
  real ( kind = 8 ) time
  real ( kind = 8 ) w
  real ( kind = 8 ), dimension(quad_num) :: wq
  real ( kind = 8 ) x
  real ( kind = 8 ), dimension(quad_num,element_num) :: xq
  real ( kind = 8 ) y
  real ( kind = 8 ), dimension(quad_num,element_num) :: yq
!
!  Initialize the arrays to zero.
!
  f(1:node_num) = 0.0D+00
  a(1:3*ib+1,1:node_num) = 0.0D+00
!
!  The actual values of A and F are determined by summing up
!  contributions from all the elements.
!
  do element = 1, element_num
!
!  Consider a quadrature point QUAD, with coordinates (X,Y).
!
    do quad = 1, quad_num

      x = xq(quad,element)
      y = yq(quad,element)
      w = element_area(element) * wq(quad)
!
!  Consider one of the basis functions, which will play the
!  role of test function in the integral.
!
!  We will generate an integral for the test function; that is, a basis
!  function associated with a degree of freedom.
!
!  If the degree of freedom associated with the node is
!  constrained by a boundary condition, then the finite element
!  integral we set up here will be modified or replaced later
!  (see subroutine BOUNDARY).
!
      do test = 1, element_order

        node = element_node(test,element)

        call qbf ( x, y, element, test, node_xy, element_node, &
          element_num, element_order, node_num, bi, dbidx, dbidy )

        f(node) = f(node) + w * rhs ( x, y, time ) * bi
!
!  Consider a basis function, which is used to form the
!  value of the solution function.
!
!  If this basis function is associated with a boundary condition,
!  then subtract the term from the right hand side.
!
!  Otherwise, this term is included in the system matrix.
!
!  Logically, this term goes in entry A(I,J).  Because of the
!  band matrix storage, entry (I,J) is actually stored in
!  A(I-J+2*NHBA+1,J).
!
        do basis = 1, element_order

          j = element_node(basis,element)

          call qbf ( x, y, element, basis, node_xy, element_node, &
            element_num, element_order, node_num, bj, dbjdx, dbjdy )

          aij = dbidx * dbjdx + dbidy * dbjdy

          a(node-j+2*ib+1,j) = a(node-j+2*ib+1,j) + w * aij

        end do

      end do

    end do

  end do

  return
end
subroutine bandwidth ( element_order, element_num, element_node, node_num, &
  nhba )

!*****************************************************************************80
!
!! BANDWIDTH determines the bandwidth of the coefficient matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 April 2006
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!  Parameters:
!
!    Input, integer ELEMENT_ORDER, the number of local nodes per element.
!
!    Input, integer ELEMENT_NUM, the number of elements.
!
!    Input, integer ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM);
!    ELEMENT_NODE(I,J) is the global index of local node I in element J.
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Output, integer NHBA, the half bandwidth of the matrix.
!
  implicit none

  integer element_num
  integer element_order
  integer node_num

  integer element
  integer element_node(element_order,element_num)
  integer i
  integer iln
  integer in
  integer j
  integer jln
  integer jn
  integer nhba

  nhba = 0
  do element = 1, element_num
    do iln = 1, element_order
      i = element_node(iln,element)
      do jln = 1, element_order
        j = element_node(jln,element)
        nhba = max ( nhba, j - i )
      end do
    end do
  end do

  return
end
subroutine compare ( node_num, node_xy, time, u )

!*****************************************************************************80
!
!! COMPARE compares the exact and computed solution at the nodes.
!
!  Discussion:
!
!    This is a rough comparison, done only at the nodes.  Such a pointwise
!    comparison is easy, because the value of the finite element
!    solution is exactly the value of the finite element coefficient
!    associated with that node.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 April 2006
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!  Parameters:
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, real ( kind = 8 ) NODE_XY(2,NODE_NUM), the
!    coordinates of the nodes.
!
!    Input, real ( kind = 8 ) TIME, the curent time.
!
!    Input, real ( kind = 8 ) U(NODE_NUM), the solution vector of the finite
!    element system.
!
  implicit none

  integer node_num

  real ( kind = 8 ) dudx_exact(node_num)
  real ( kind = 8 ) dudy_exact(node_num)
  integer i
  integer node
  real ( kind = 8 ) node_xy(2,node_num)
  real ( kind = 8 ) time
  real ( kind = 8 ) u(node_num)
  real ( kind = 8 ) u_exact(node_num)
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  call exact_u ( node_num, node_xy, time, u_exact, dudx_exact, dudy_exact )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'COMPARE:'
  write ( *, '(a)' ) '  Compare computed and exact solutions at the nodes.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  TIME = ', time
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           Y           U           U'
  write ( *, '(a)' ) '                                exact      computed'
  write ( *, '(a)' ) ' '

  do node = 1, node_num

    x = node_xy(1,node)
    y = node_xy(2,node)

    write ( *, '(2x,4f12.4)' ) x, y, u_exact(node), u(node)

  end do

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
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 April 2006
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
    if ( l /= m ) then
      temp = a(l,k)
      a(l,k) = a(m,k)
      a(m,k) = temp
    end if
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
        temp    = a(l,j)
        a(l,j)  = a(mm,j)
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
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 April 2006
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
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 April 2006
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
  integer job
  integer k
  integer l
  integer la
  integer lb
  integer lm
  integer m
  real ( kind = 8 ) temp

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
          temp = b(l)
          b(l) = b(k)
          b(k) = temp
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
          temp = b(l)
          b(l) = b(k)
          b(k) = temp
        end if

      end do

    end if

  end if

  return
end
subroutine element_write ( element_order, element_num, element_node, file_name )

!*****************************************************************************80
!
!! ELEMENT_WRITE writes the element information to a text file.
!
!  Discussion:
!
!    The element information and the node/solution information can be read
!    into another program to make plots.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ELEMENT_ORDER, the number of local nodes per element.
!
!    Input, integer ELEMENT_NUM, the number of elements.
!
!    Input, integer ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM);
!    ELEMENT_NODE(I,J) is the global index of local node I in element J.
!
!    Input, character ( len = * ) FILE_NAME, the name of the output file.
!
  implicit none

  integer element_num
  integer element_order

  integer element
  integer element_node(element_order,element_num)
  character ( len = * ) file_name
  integer output_status
  integer output_unit

  call get_unit ( output_unit )

  open ( unit = output_unit, file = file_name, status = 'replace', &
    iostat = output_status )

  if ( output_status /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ELEMENTS_WRITE - Warning!'
    write ( *, '(a)' ) '  Could not write the element file.'
    return
  end if

  do element = 1, element_num

    write ( output_unit, '(2x,i4,2x,i4,2x,i4,2x,i4,2x,i4,2x,i4)' ) &
      element_node(1:element_order,element)

  end do

  close ( unit = output_unit )

  return
end
subroutine errors ( element_area, element_node, node_xy, u, &
  element_num, element_order, node_num, time, el2, eh1 )

!*****************************************************************************80
!
!! ERRORS calculates the error in the L2 norm and H1 seminorm.
!
!  Discussion:
!
!    This routine uses a 13 point quadrature rule in each element,
!    in order to estimate the values of
!
!      EL2(t) = Sqrt ( Integral ( U(x,y,t) - Uh(x,y,t) )**2 dx dy )
!
!      EH1(t) = Sqrt ( Integral ( Ux(x,y,t) - Uhx(x,y,t) )**2 +
!                               ( Uy(x,y,t) - Uhy(x,y,t) )**2 dx dy )
!
!    Here U is the exact solution, and Ux and Uy its spatial derivatives,
!    as evaluated by a user-supplied routine.
!
!    Uh, Uhx and Uhy are the computed solution and its spatial derivatives,
!    as specified by the computed finite element solution.
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
!    FORTRAN90 version by John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ELEMENT_AREA(ELEMENT_NUM), the area of elements.
!
!    Input, integer ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM);
!    ELEMENT_NODE(I,J) is the global index of local node I in element J.
!
!    Input, real ( kind = 8 ) NODE_XY(2,NODE_NUM), the X and Y
!    coordinates of nodes.
!
!    Input, real ( kind = 8 ) U(NODE_NUM), the coefficients of the solution.
!
!    Input, integer ELEMENT_NUM, the number of elements.
!
!    Input, integer ELEMENT_ORDER, the number of nodes used to form one element.
!
!    Input, integer QUAD2_NUM, the number of points in the quadrature rule.
!    This is actually fixed at 13.
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, real ( kind = 8 ) TIME, the current time.
!
!    Output, real ( kind = 8 ) EL2, the L2 error.
!
!    Output, real ( kind = 8 ) EH1, the H1 seminorm error.
!
!  Local Parameters:
!
!    Local, real ( kind = 8 ) AR, the weight for a given quadrature point
!    in a given element.
!
!    Local, real ( kind = 8 ) BB, BX, BY, a basis function and its first
!    derivatives evaluated at a particular quadrature point.
!
!    Local, real ( kind = 8 ) EH1, the H1 seminorm error.
!
!    Local, real ( kind = 8 ) EL2, the L2 error.
!
!    Local, real ( kind = 8 ) UEX, UEXX, UEXY, the exact solution and its first
!    derivatives evaluated at a particular quadrature point.
!
!    Local, real ( kind = 8 ) UH, UHX, UHY, the computed solution and its first
!    derivatives evaluated at a particular quadrature point.
!
!    Local, real ( kind = 8 ) WQE(QUAD2_NUM), stores the quadrature weights.
!
!    Local, real ( kind = 8 ) X, Y, the coordinates of a particular
!    quadrature point.
!
!    Local, real ( kind = 8 ) XQE(QUAD2_NUM), YQE(QUAD2_NUM), the
!    quadrature points in a given element.
!
  implicit none

  integer element_num
  integer element_order
  integer node_num
  integer, parameter :: quad2_num = 13

  real ( kind = 8 ) ar
  real ( kind = 8 ) bi
  real ( kind = 8 ) dbidx
  real ( kind = 8 ) dbidy
  real ( kind = 8 ) dudx_exact(1)
  real ( kind = 8 ) dudxh
  real ( kind = 8 ) dudy_exact(1)
  real ( kind = 8 ) dudyh
  real ( kind = 8 ) eh1
  real ( kind = 8 ) el2
  integer element
  real ( kind = 8 ), dimension(element_num) :: element_area
  integer, dimension(element_order,element_num) :: element_node
  integer i
  integer in1
  integer ip
  real ( kind = 8 ), dimension(2,node_num) :: node_xy
  integer quad
  real ( kind = 8 ) time
  real ( kind = 8 ), dimension(node_num) :: u
  real ( kind = 8 ) u_exact(1)
  real ( kind = 8 ) uh
  real ( kind = 8 ), dimension(quad2_num) :: wqe
  real ( kind = 8 ) x
  real ( kind = 8 ), dimension(quad2_num) :: xqe
  real ( kind = 8 ) xy(2)
  real ( kind = 8 ) y
  real ( kind = 8 ), dimension(quad2_num) :: yqe

  el2 = 0.0D+00
  eh1 = 0.0D+00
!
!  For each element, retrieve the nodes, area, quadrature weights,
!  and quadrature points.
!
  do element = 1, element_num

    call quad_e ( node_xy, element_node, element, &
      element_num, element_order, node_num, quad2_num, wqe, xqe, yqe )
!
!  For each quadrature point, evaluate the computed solution and its X and
!  Y derivatives.
!
    do quad = 1, quad2_num

      ar = element_area(element) * wqe(quad)
      x = xqe(quad)
      y = yqe(quad)

      uh = 0.0D+00
      dudxh = 0.0D+00
      dudyh = 0.0D+00

      do in1 = 1, element_order

        i = element_node(in1,element)

        call qbf ( x, y, element, in1, node_xy, element_node, &
          element_num, element_order, node_num, bi, dbidx, dbidy )

        uh    = uh    + bi    * u(i)
        dudxh = dudxh + dbidx * u(i)
        dudyh = dudyh + dbidy * u(i)

      end do
!
!  Evaluate the exact solution and its X and Y derivatives.
!
      xy(1) = x
      xy(2) = y

      call exact_u ( 1, xy, time, u_exact, dudx_exact, dudy_exact )
!
!  Add the weighted value at this quadrature point to the quadrature sum.
!
      el2 = el2 + ar * ( uh    - u_exact(1) )**2

      eh1 = eh1 + ar * ( ( dudxh - dudx_exact(1) )**2 &
                       + ( dudyh - dudy_exact(1) )**2 )

    end do

  end do

  el2 = sqrt ( el2 )
  eh1 = sqrt ( eh1 )

  write ( *, '(3g14.6)' ) time, el2, eh1

  return
end
subroutine exact_u ( node_num, node_xy, time, u, dudx, dudy )

!*****************************************************************************80
!
!! EXACT_U calculates the exact solution.
!
!  Discussion:
!
!    It is assumed that the user knows the exact solution and its
!    derivatives.  This, of course, is NOT true for a real computation.
!    But for this code, we are interested in studying the convergence
!    behavior of the approximations, and so we really need to assume
!    we know the correct solution.
!
!    As a convenience, this single routine is used for several purposes:
!
!    * it supplies the initial value function H(X,Y,T);
!    * it supplies the boundary value function G(X,Y,T);
!    * it is used by the COMPARE routine to make a node-wise comparison
!      of the exact and approximate solutions.
!    * it is used by the ERRORS routine to estimate the integrals of
!      the L2 and H1 errors of approximation.
!
!    DUDX and DUDY are only needed for the ERRORS calculation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 April 2006
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!  Parameters:
!
!    Input, integer NODE_NUM, the number of nodes at which
!    a value is desired.
!
!    Input, real ( kind = 8 ) NODE_XY(2,NODE_NUM), the coordinates of
!    the points where a value is desired.
!
!    Input, real ( kind = 8 ) TIME, the current time.
!
!    Output, real ( kind = 8 ) U(NODE_NUM), the exact solution.
!
!    Output, real ( kind = 8 ) DUDX(NODE_NUM), DUDY(NODE_NUM),
!    the X and Y derivatives of the exact solution.
!
  implicit none

  integer node_num

! real ( kind = 8 ) dudt(node_num)
  real ( kind = 8 ) dudx(node_num)
  real ( kind = 8 ) dudy(node_num)
  integer node
  real ( kind = 8 ) node_xy(2,node_num)
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) time
  real ( kind = 8 ) u(node_num)
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  do node = 1, node_num

    x = node_xy(1,node)
    y = node_xy(2,node)

    u(node)    =      sin ( pi * x ) * sin ( pi * y ) * exp ( - time )
!   dudt(node) =    - sin ( pi * x ) * sin ( pi * y ) * exp ( - time )
    dudx(node) = pi * cos ( pi * x ) * sin ( pi * y ) * exp ( - time )
    dudy(node) = pi * sin ( pi * x ) * cos ( pi * y ) * exp ( - time )

  end do

  return
end
subroutine file_name_inc ( file_name )

!*****************************************************************************80
!
!! FILE_NAME_INC increments a partially numeric filename.
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
!    12 October 2003
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
  integer digit
  character ( len = * ) file_name
  integer i
  integer lens

  lens = len_trim ( file_name )

  do i = lens, 1, -1

    c = file_name(i:i)

    if ( lge ( c, '0' ) .and. lle ( c, '9' ) ) then

      digit = ichar ( c ) - 48
      digit = digit + 1

      if ( digit == 10 ) then
        digit = 0
      end if

      c = char ( digit + 48 )

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
subroutine grid_t6 ( nx, ny, element_order, element_num, element_node )

!*****************************************************************************80
!
!! GRID_T6 produces a grid of pairs of 6 node triangles.
!
!  Example:
!
!    Input:
!
!      NX = 4, NY = 3
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
!  Diagram:
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
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NX, NY, controls the number of elements along the
!    X and Y directions.  The number of elements will be
!    2 * ( NX - 1 ) * ( NY - 1 ).
!
!    Input, integer ELEMENT_ORDER, the number of local nodes per element.
!
!    Input, integer ELEMENT_NUM, the number of elements.
!
!    Output, integer ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM);
!    ELEMENT_NODE(I,J) is the index of the I-th node of the J-th element.
!
  implicit none

  integer element_num
  integer element_order

  integer c
  integer e
  integer element
  integer element_node(element_order,element_num)
  integer i
  integer j
  integer n
  integer ne
  integer nw
  integer nx
  integer ny
  integer s
  integer se
  integer sw
  integer w

  element = 0

  do j = 1, ny - 1
    do i = 1, nx - 1

      sw = ( j - 1 ) * 2 * ( 2 * nx - 1 ) + 2 * i - 1
      w  = sw + 1
      nw = sw + 2

      s  = sw + 2 * nx - 1
      c  = s  + 1
      n  = s  + 2

      se = s  + 2 * nx - 1
      e  = se + 1
      ne = se + 2

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
subroutine i4vec_print_some ( n, a, max_print, title )

!*****************************************************************************80
!
!! I4VEC_PRINT_SOME prints "some" of an I4VEC.
!
!  Discussion:
!
!    The user specifies MAX_PRINT, the maximum number of lines to print.
!
!    If N, the size of the vector, is no more than MAX_PRINT, then
!    the entire vector is printed, one entry per line.
!
!    Otherwise, if possible, the first MAX_PRINT-2 entries are printed,
!    followed by a line of periods suggesting an omission,
!    and the last entry.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries of the vector.
!
!    Input, integer A(N), the vector to be printed.
!
!    Input, integer MAX_PRINT, the maximum number of lines to print.
!
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none

  integer n

  integer a(n)
  integer i
  integer max_print
  character ( len = * ) title

  if ( max_print <= 0 ) then
    return
  end if

  if ( n <= 0 ) then
    return
  end if

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
    write ( *, '(a)' ) ' '
  end if

  if ( n <= max_print ) then

    do i = 1, n
      write ( *, '(i8,2x,i10)' ) i, a(i)
    end do

  else if ( 3 <= max_print ) then

    do i = 1, max_print-2
      write ( *, '(i8,2x,i10)' ) i, a(i)
    end do
    write ( *, '(a)' ) '......  ..............'
    i = n
    write ( *, '(i8,2x,i10)' ) i, a(i)

  else

    do i = 1, max_print - 1
      write ( *, '(i8,2x,i10)' ) i, a(i)
    end do
    i = max_print
    write ( *, '(i8,2x,i10,2x,a)' ) i, a(i), '...more entries...'

  end if

  return
end
subroutine node_boundary_set ( nx, ny, node_num, node_boundary )

!*****************************************************************************80
!
!! NODE_BOUNDARY_SET assigns an unknown value index at each node.
!
!  Discussion:
!
!    Every node is assigned a value which indicates whether it is
!    an interior node, or a boundary node.
!
!  Example:
!
!    On a simple 5 by 5 grid, where the nodes are numbered starting
!    at the lower left, and increasing in X first, we would have the
!    following values of NODE_BOUNDARY:
!
!       1  1  1  1  1
!       1  0  0  0  1
!       1  0  0  0  1
!       1  0  0  0  1
!       1  1  1  1  1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 May 2005
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!  Parameters:
!
!    Input, integer NX, NY, the number of elements in the X and Y directions.
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Output, integer NODE_BOUNDARY(NODE_NUM), is
!    0, if a node is an interior node;
!    1, if a node is a Dirichlet boundary node.
!
  implicit none

  integer node_num

  integer i
  integer in
  integer j
  integer node
  integer node_boundary(node_num)
  integer nx
  integer ny

  node = 0

  do j = 1, 2 * ny - 1

    do i = 1, 2 * nx - 1

      node = node + 1
      if ( j == 1 .or. &
           j == 2 * ny - 1 .or. &
           i == 1 .or. &
           i == 2 * nx - 1 ) then
        node_boundary(node) = 1
      else
        node_boundary(node) = 0
      end if

    end do
  end do

  return
end
subroutine nodes_plot ( file_name, node_num, node_xy, node_label )

!*****************************************************************************80
!
!! NODES_PLOT plots a pointset.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
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
    write ( *, '(a)' ) 'NODES_PLOT - Fatal error!'
    write ( *, '(a)' ) '  Can not open output file.'
    return
  end if

  write ( file_unit, '(a)' ) '%!PS-Adobe-3.0 EPSF-3.0'
  write ( file_unit, '(a)' ) '%%Creator: nodes_plot.f90'
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
subroutine nodes_write ( node_num, node_xy, output_filename )

!*****************************************************************************80
!
!! NODES_WRITE writes the nodes to a file.
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
!  Parameters:
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, real ( kind = 8 ) NODE_XY(2,NODE_NUM), the nodes.
!
!    Input, character ( len = * ) OUTPUT_FILENAME, the name of the file
!    in which the data should be stored.
!
  implicit none

  integer node_num

  integer node
  real ( kind = 8 ), dimension(2,node_num) :: node_xy
  character ( len = * ) :: output_filename
  integer output_status
  integer output_unit

  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_filename, status = 'replace', &
    iostat = output_status )

  if ( output_status /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'NODES_WRITE - Warning!'
    write ( *, '(a)' ) '  Could not write the node file.'
    return
  end if

  do node = 1, node_num
    write ( output_unit, '(f8.4,2x,f8.4,2x,g14.6)' ) node_xy(1:2,node)
  end do

  close ( unit = output_unit )

  return
end
subroutine qbf ( x, y, element, inode, node_xy, element_node, &
  element_num, element_order, node_num, b, dbdx, dbdy )

!*****************************************************************************80
!
!! QBF evaluates the quadratic basis functions.
!
!  Discussion:
!
!    This routine assumes that the "midpoint" nodes are, in fact,
!    exactly the average of the two extreme nodes.  This is NOT true
!    for a general quadratic triangular element.
!
!    Assuming this property of the midpoint nodes makes it easy to
!    determine the values of (R,S) in the reference element that
!    correspond to (X,Y) in the physical element.
!
!    Once we know the (R,S) coordinates, it's easy to evaluate the
!    basis functions and derivatives.
!
!  The physical element T6:
!
!    In this picture, we don't mean to suggest that the bottom of
!    the physical triangle is horizontal.  However, we do assume that
!    each of the sides is a straight line, and that the intermediate
!    points are exactly halfway on each side.
!
!    |
!    |
!    |        3
!    |       / \
!    |      /   \
!    Y     6     5
!    |    /       \
!    |   /         \
!    |  1-----4-----2
!    |
!    +--------X-------->
!
!  Reference element T6:
!
!    In this picture of the reference element, we really do assume
!    that one side is vertical, one horizontal, of length 1.
!
!    |
!    |
!    1  3
!    |  |\
!    |  | \
!    S  6  5
!    |  |   \
!    |  |    \
!    0  1--4--2
!    |
!    +--0--R--1-------->
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, Y, the (global) coordinates of the point
!    at which the basis function is to be evaluated.
!
!    Input, integer ELEMENT, the index of the element which contains the point.
!
!    Input, integer INODE, the local index (between 1 and 6) that
!    specifies which basis function is to be evaluated.
!
!    Input, real ( kind = 8 ) NODE_XY(2,NODE_NUM), the X and Y
!    coordinates of nodes.
!
!    Input, integer ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM);
!    ELEMENT_NODE(I,J) is the global index of local node I in element J.
!
!    Input, integer ELEMENT_NUM, the number of elements.
!
!    Input, integer ELEMENT_ORDER, the number of nodes used to form one element.
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Output, real ( kind = 8 ) B, DBDX, DBDY, the value of the basis function
!    and its X and Y derivatives at (X,Y).
!
  implicit none

  integer element_num
  integer element_order
  integer node_num

  real ( kind = 8 ) b
  real ( kind = 8 ) dbdr
  real ( kind = 8 ) dbds
  real ( kind = 8 ) dbdx
  real ( kind = 8 ) dbdy
  real ( kind = 8 ) det
  real ( kind = 8 ) drdx
  real ( kind = 8 ) drdy
  real ( kind = 8 ) dsdx
  real ( kind = 8 ) dsdy
  integer element
  integer, dimension(element_order,element_num) :: element_node
  integer inode
  real ( kind = 8 ), dimension(2,node_num) :: node_xy
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) x
  real ( kind = 8 ) xn(6)
  real ( kind = 8 ) y
  real ( kind = 8 ) yn(6)

  xn(1:6) = node_xy(1,element_node(1:6,element))
  yn(1:6) = node_xy(2,element_node(1:6,element))
!
!  Determine the (R,S) coordinates corresponding to (X,Y).
!
!  What is happening here is that we are solving the linear system:
!
!    ( X2-X1  X3-X1 ) * ( R ) = ( X - X1 )
!    ( Y2-Y1  Y3-Y1 )   ( S )   ( Y - Y1 )
!
!  by computing the inverse of the coefficient matrix and multiplying
!  it by the right hand side to get R and S.
!
!  The values of dRdX, dRdY, dSdX and dSdY are easily from the formulas
!  for R and S.
!
  det =   ( xn(2) - xn(1) ) * ( yn(3) - yn(1) ) &
        - ( xn(3) - xn(1) ) * ( yn(2) - yn(1) )

  r = ( ( yn(3) - yn(1) ) * ( x     - xn(1) ) &
      + ( xn(1) - xn(3) ) * ( y     - yn(1) ) ) / det

  drdx = ( yn(3) - yn(1) ) / det
  drdy = ( xn(1) - xn(3) ) / det

  s = ( ( yn(1) - yn(2) ) * ( x     - xn(1) ) &
      + ( xn(2) - xn(1) ) * ( y     - yn(1) ) ) / det

  dsdx = ( yn(1) - yn(2) ) / det
  dsdy = ( xn(2) - xn(1) ) / det
!
!  The basis functions can now be evaluated in terms of the
!  reference coordinates R and S.  It's also easy to determine
!  the values of the derivatives with respect to R and S.
!
  if ( inode == 1 ) then

    b    =   2.0D+00 *     ( 1.0D+00 - r - s ) * ( 0.5D+00 - r - s )
    dbdr = - 3.0D+00 + 4.0D+00 * r + 4.0D+00 * s
    dbds = - 3.0D+00 + 4.0D+00 * r + 4.0D+00 * s

  else if ( inode == 2 ) then

    b    =   2.0D+00 * r * ( r - 0.5D+00 )
    dbdr = - 1.0D+00 + 4.0D+00 * r
    dbds =   0.0D+00

  else if ( inode == 3 ) then

    b    =   2.0D+00 * s * ( s - 0.5D+00 )
    dbdr =   0.0D+00
    dbds = - 1.0D+00               + 4.0D+00 * s

  else if ( inode == 4 ) then

    b    =   4.0D+00 * r * ( 1.0D+00 - r - s )
    dbdr =   4.0D+00 - 8.0D+00 * r - 4.0D+00 * s
    dbds =           - 4.0D+00 * r

  else if ( inode == 5 ) then

    b    =   4.0D+00 * r * s
    dbdr =                           4.0D+00 * s
    dbds =             4.0D+00 * r

  else if ( inode == 6 ) then

    b    =   4.0D+00 * s * ( 1.0D+00 - r - s )
    dbdr =                         - 4.0D+00 * s
    dbds =   4.0D+00 - 4.0D+00 * r - 8.0D+00 * s

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'QBF - Fatal error!'
    write ( *, '(a,i8)' ) '  Request for local basis function INODE = ', inode
    stop

  end if
!
!  We need to convert the derivative information from (R(X,Y),S(X,Y))
!  to (X,Y) using the chain rule.
!
  dbdx = dbdr * drdx + dbds * dsdx
  dbdy = dbdr * drdy + dbds * dsdy

  return
end
subroutine quad_a ( node_xy, element_node, element_num, node_num, &
  element_order, wq, xq, yq )

!*****************************************************************************80
!
!! QUAD_A sets the quadrature rule for assembly.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 2004
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) NODE_XY(2,NODE_NUM), the X and
!    Y coordinates of nodes.
!
!    Input, integer ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM);
!    ELEMENT_NODE(I,J) is the global index of local node I in element J.
!
!    Input, integer ELEMENT_NUM, the number of elements.
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, integer ELEMENT_ORDER, the number of nodes used to form one element.
!
!    Output, real ( kind = 8 ) WQ(3), quadrature weights.
!
!    Output, real ( kind = 8 ) XQ(3,ELEMENT_NUM), YQ(3,ELEMENT_NUM), the
!    coordinates of the quadrature points in each element.
!
  implicit none

  integer element_num
  integer element_order
  integer node_num

  integer element
  integer, dimension(element_order,element_num) :: element_node
  integer ip1
  integer ip2
  integer ip3
  real ( kind = 8 ), dimension(2,node_num) :: node_xy
  real ( kind = 8 ), dimension(3) :: wq
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ), dimension(3,element_num) :: xq
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) y3
  real ( kind = 8 ), dimension(3,element_num) :: yq

  wq(1) = 1.0D+00 / 3.0D+00
  wq(2) = wq(1)
  wq(3) = wq(1)

  do element = 1, element_num

    ip1 = element_node(1,element)
    ip2 = element_node(2,element)
    ip3 = element_node(3,element)

    x1 = node_xy(1,ip1)
    x2 = node_xy(1,ip2)
    x3 = node_xy(1,ip3)

    y1 = node_xy(2,ip1)
    y2 = node_xy(2,ip2)
    y3 = node_xy(2,ip3)

    xq(1,element) = 0.5D+00 * ( x1 + x2 )
    xq(2,element) = 0.5D+00 * ( x2 + x3 )
    xq(3,element) = 0.5D+00 * ( x1 + x3 )

    yq(1,element) = 0.5D+00 * ( y1 + y2 )
    yq(2,element) = 0.5D+00 * ( y2 + y3 )
    yq(3,element) = 0.5D+00 * ( y1 + y3 )

  end do

  return
end
subroutine quad_e ( node_xy, element_node, element, element_num, &
  element_order, node_num, quad2_num, wqe, xqe, yqe )

!*****************************************************************************80
!
!! QUAD_E sets up the quadrature rule for error integration.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 2004
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) NODE_XY(2,NODE_NUM), the X and Y
!    coordinates of nodes.
!
!    Input, integer ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM);
!    ELEMENT_NODE(I,J) is the global index of local node I in element J.
!
!    Input, integer ELEMENT, the index of the element for which the quadrature
!    points are to be computed.
!
!    Input, integer ELEMENT_NUM, the number of elements.
!
!    Input, integer ELEMENT_ORDER, the number of nodes used to form one element.
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, integer QUAD2_NUM, the number of points in the quadrature rule.
!    This is actually fixed at 13.
!
!    Output, real ( kind = 8 ) WQE(QUAD2_NUM), the quadrature weights.
!
!    Output, real ( kind = 8 ) XQE(QUAD2_NUM), YQE(QUAD2_NUM), the
!    quadrature points.
!
  implicit none

  integer element_num
  integer element_order
  integer node_num
  integer quad2_num

  integer element
  integer, dimension(element_order,element_num) :: element_node
  integer i
  integer ii
  integer iii
  integer ip1
  integer ip2
  integer ip3
  real ( kind = 8 ), dimension(2,node_num) :: node_xy
  real ( kind = 8 ), dimension(quad2_num) :: wqe
  real ( kind = 8 ), dimension(quad2_num) :: xqe
  real ( kind = 8 ), dimension(quad2_num) :: yqe
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) y3
  real ( kind = 8 ) z1
  real ( kind = 8 ) z2
  real ( kind = 8 ) z3
  real ( kind = 8 ) z4
  real ( kind = 8 ) z5
  real ( kind = 8 ) z6
  real ( kind = 8 ) z7

  do i = 1, 3
    wqe(i) = 0.175615257433204D+00
    ii = i + 3
    wqe(ii) = 0.053347235608839D+00
    ii = i + 6
    iii = ii + 3
    wqe(ii) = 0.077113760890257D+00
    wqe(iii) = wqe(ii)
  end do

  wqe(13) = -0.14957004446767D+00

  z1 = 0.479308067841923D+00
  z2 = 0.260345966079038D+00
  z3 = 0.869739794195568D+00
  z4 = 0.065130102902216D+00
  z5 = 0.638444188569809D+00
  z6 = 0.312865496004875D+00
  z7 = 0.048690315425316D+00

  ip1 = element_node(1,element)
  ip2 = element_node(2,element)
  ip3 = element_node(3,element)
  x1 = node_xy(1,ip1)
  x2 = node_xy(1,ip2)
  x3 = node_xy(1,ip3)
  y1 = node_xy(2,ip1)
  y2 = node_xy(2,ip2)
  y3 = node_xy(2,ip3)

  xqe( 1) = z1 * x1 + z2 * x2 + z2 * x3
  yqe( 1) = z1 * y1 + z2 * y2 + z2 * y3
  xqe( 2) = z2 * x1 + z1 * x2 + z2 * x3
  yqe( 2) = z2 * y1 + z1 * y2 + z2 * y3
  xqe( 3) = z2 * x1 + z2 * x2 + z1 * x3
  yqe( 3) = z2 * y1 + z2 * y2 + z1 * y3
  xqe( 4) = z3 * x1 + z4 * x2 + z4 * x3
  yqe( 4) = z3 * y1 + z4 * y2 + z4 * y3
  xqe( 5) = z4 * x1 + z3 * x2 + z4 * x3
  yqe( 5) = z4 * y1 + z3 * y2 + z4 * y3
  xqe( 6) = z4 * x1 + z4 * x2 + z3 * x3
  yqe( 6) = z4 * y1 + z4 * y2 + z3 * y3
  xqe( 7) = z5 * x1 + z6 * x2 + z7 * x3
  yqe( 7) = z5 * y1 + z6 * y2 + z7 * y3
  xqe( 8) = z5 * x1 + z7 * x2 + z6 * x3
  yqe( 8) = z5 * y1 + z7 * y2 + z6 * y3
  xqe( 9) = z6 * x1 + z5 * x2 + z7 * x3
  yqe( 9) = z6 * y1 + z5 * y2 + z7 * y3
  xqe(10) = z6 * x1 + z7 * x2 + z5 * x3
  yqe(10) = z6 * y1 + z7 * y2 + z5 * y3
  xqe(11) = z7 * x1 + z5 * x2 + z6 * x3
  yqe(11) = z7 * y1 + z5 * y2 + z6 * y3
  xqe(12) = z7 * x1 + z6 * x2 + z5 * x3
  yqe(12) = z7 * y1 + z6 * y2 + z5 * y3
  xqe(13) = ( x1 + x2 + x3 ) / 3.0D+00
  yqe(13) = ( y1 + y2 + y3 ) / 3.0D+00

  return
end
subroutine r8vec_print_some ( n, a, max_print, title )

!*****************************************************************************80
!
!! R8VEC_PRINT_SOME prints "some" of an R8VEC.
!
!  Discussion:
!
!    The user specifies MAX_PRINT, the maximum number of lines to print.
!
!    If N, the size of the vector, is no more than MAX_PRINT, then
!    the entire vector is printed, one entry per line.
!
!    Otherwise, if possible, the first MAX_PRINT-2 entries are printed,
!    followed by a line of periods suggesting an omission,
!    and the last entry.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 September 2003
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
!    Input, integer MAX_PRINT, the maximum number of lines to print.
!
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none

  integer n

  real ( kind = 8 ) a(n)
  integer i
  integer max_print
  character ( len = * ) title

  if ( max_print <= 0 ) then
    return
  end if

  if ( n <= 0 ) then
    return
  end if

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
    write ( *, '(a)' ) ' '
  end if

  if ( n <= max_print ) then

    if ( all ( a(1:n) == aint ( a(1:n) ) ) ) then
      do i = 1, n
        write ( *, '(i8,2x,i8)' ) i, int ( a(i) )
      end do
    else if ( all ( abs ( a(1:n) ) < 1000000.0D+00 ) ) then
      do i = 1, n
        write ( *, '(i8,2x,f14.6)' ) i, a(i)
      end do
    else
      do i = 1, n
        write ( *, '(i8,2x,g14.6)' ) i, a(i)
      end do
    end if

  else if ( 3 <= max_print ) then

    if ( all ( a(1:max_print-2) == aint ( a(1:max_print-2) ) ) ) then
      do i = 1, max_print-2
        write ( *, '(i8,2x,i8)' ) i, int ( a(i) )
      end do
    else if ( all ( abs ( a(1:max_print-2) ) < 1000000.0D+00 ) ) then
      do i = 1, max_print-2
        write ( *, '(i8,2x,f14.6)' ) i, a(i)
      end do
    else
      do i = 1, max_print-2
        write ( *, '(i8,2x,g14.6)' ) i, a(i)
      end do
    end if

    write ( *, '(a)' ) '......  ..............'
    i = n

    if ( a(i) == aint ( a(i) ) ) then
      write ( *, '(i8,2x,i8)' ) i, int ( a(i) )
    else if (  abs ( a(i) ) < 1000000.0D+00 ) then
      write ( *, '(i8,2x,f14.6)' ) i, a(i)
    else
      write ( *, '(i8,2x,g14.6)' ) i, a(i)
    end if

  else

    if ( all ( a(1:max_print-1) == aint ( a(1:max_print-1) ) ) ) then
      do i = 1, max_print-1
        write ( *, '(i8,2x,i8)' ) i, int ( a(i) )
      end do
    else if ( all ( abs ( a(1:max_print-1) ) < 1000000.0D+00 ) ) then
      do i = 1, max_print-1
        write ( *, '(i8,2x,f14.6)' ) i, a(i)
      end do
    else
      do i = 1, max_print-1
        write ( *, '(i8,2x,g14.6)' ) i, a(i)
      end do
    end if

    i = max_print

    if ( a(i) == aint ( a(i) ) ) then
      write ( *, '(i8,2x,i8,a)' ) i, int ( a(i) ), '...more entries...'
    else if (  abs ( a(i) ) < 1000000.0D+00 ) then
      write ( *, '(i8,2x,f14.6,a)' ) i, a(i), '...more entries...'
    else
      write ( *, '(i8,2x,g14.6,a)' ) i, a(i), '...more entries...'
    end if

  end if

  return
end
function rhs ( x, y, time )

!*****************************************************************************80
!
!! RHS gives the right-hand side of the differential equation.
!
!  Discussion:
!
!    The function specified here depends on the problem being
!    solved.  The user must be sure that RHS and EXACT_U are consistent.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 April 2006
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, Y, a point in the region.
!
!    Input, real ( kind = 8 ) TIME, the current time.
!
!    Output, real ( kind = 8 ) RHS, the value of the right
!    hand side of the differential equation at (X,Y,TIME).
!
  implicit none

  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) rhs
  real ( kind = 8 ) time
  real ( kind = 8 ) ut
  real ( kind = 8 ) uxx
  real ( kind = 8 ) uyy
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  ut =            - sin ( pi * x ) * sin ( pi * y ) * exp ( - time )
  uxx = - pi * pi * sin ( pi * x ) * sin ( pi * y ) * exp ( - time )
  uyy = - pi * pi * sin ( pi * x ) * sin ( pi * y ) * exp ( - time )

  rhs = ut - uxx - uyy

  return
end
subroutine solution_write ( node_num, u, output_file_name )

!*****************************************************************************80
!
!! SOLUTION_WRITE writes the solution to a file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, real ( kind = 8 ) U(NODE_NUM), the coefficients of the solution.
!
!    Input, character ( len = * ) OUTPUT_FILE_NAME, the name of the file
!    in which the data should be stored.
!
  implicit none

  integer node_num

  integer node
  character ( len = * ) :: output_file_name
  integer output_status
  integer output_unit
  real ( kind = 8 ), dimension(node_num) :: u

  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_file_name, status = 'replace', &
    iostat = output_status )

  if ( output_status /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SOLUTION_WRITE - Warning!'
    write ( *, '(a)' ) '  Could not write the solution file "' &
      // trim ( output_file_name ) // '".'
    return
  end if

  do node = 1, node_num

    write ( output_unit, '(g14.6)' ) u(node)

  end do

  close ( unit = output_unit )

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
subroutine triangulation_order6_plot ( file_name, node_num, node_xy, tri_num, &
  triangle_node, node_show, triangle_show )

!*****************************************************************************80
!
!! TRIANGULATION_ORDER6_PLOT plots a 6-node triangulation of a pointset.
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
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
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
!    Input, integer TRI_NUM, the number of triangles.
!
!    Input, integer TRIANGLE_NODE(6,TRI_NUM), lists, for each triangle,
!    the indices of the points that form the vertices of the triangle.
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
  integer tri_num

  real ( kind = 8 ) ave_x
  real ( kind = 8 ) ave_y
  integer :: circle_size
  integer delta
  character ( len = * ) file_name
  integer file_unit
  integer i
  integer ios
  integer node
  integer node_show
  real ( kind = 8 ) node_xy(2,node_num)
  character ( len = 40 ) string
  integer triangle
  integer triangle_node(6,tri_num)
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

    do triangle = 1, tri_num

      write ( file_unit, '(a)' ) 'newpath'

      node = triangle_node(6,triangle)

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

        node = triangle_node(i,triangle)

        x_ps = int ( &
          ( ( x_max - node_xy(1,node)         ) * real ( x_ps_min, kind = 8 )   &
          + (         node_xy(1,node) - x_min ) * real ( x_ps_max, kind = 8 ) ) &
          / ( x_max                   - x_min ) )

        y_ps = int ( &
          ( ( y_max - node_xy(2,node)         ) * real ( y_ps_min, kind = 8 )   &
          + (         node_xy(2,node) - y_min ) * real ( y_ps_max, kind = 8 ) ) &
          / ( y_max                   - y_min ) )

        write ( file_unit, '(i3,2x,i3,2x,a)' ) x_ps, y_ps, ' lineto'

        node = triangle_node(i+3,triangle)

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

    do triangle = 1, tri_num

      ave_x = 0.0D+00
      ave_y = 0.0D+00

      do i = 1, 6

        node = triangle_node(i,triangle)

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

      write ( string, '(i4)' ) triangle
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
subroutine xy_set ( nx, ny, node_num, xl, xr, yb, yt, node_xy )

!*****************************************************************************80
!
!! XY_SET sets the XY coordinates of the nodes.
!
!  Discussion:
!
!    The nodes are laid out in an evenly spaced grid, in the unit square.
!
!    The first node is at the origin.  More nodes are created to the
!    right until the value of X = 1 is reached, at which point
!    the next layer is generated starting back at X = 0, and an
!    increased value of Y.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 April 2004
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!  Parameters:
!
!    Input, integer NX, NY, the number of elements in the X and
!    Y direction.
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, real ( kind = 8 ) XL, XR, YB, YT, the X coordinates of
!    the left and right sides of the rectangle, and the Y coordinates
!    of the bottom and top of the rectangle.
!
!    Output, real ( kind = 8 ) NODE_XY(2,NODE_NUM),
!    the coordinates of the nodes.
!
  implicit none

  integer node_num

  integer i
  integer j
  real ( kind = 8 ) node_xy(2,node_num)
  integer nx
  integer ny
  real ( kind = 8 ) xl
  real ( kind = 8 ) xr
  real ( kind = 8 ) yb
  real ( kind = 8 ) yt

  do j = 1, 2 * ny - 1
    do i = 1, 2 * nx - 1

      node_xy(1,i+(j-1)*(2*nx-1)) =         &
        ( real ( 2 * nx - i - 1, kind = 8 ) * xl   &
        + real (          i - 1, kind = 8 ) * xr ) &
        / real ( 2 * nx     - 2, kind = 8 )

      node_xy(2,i+(j-1)*(2*nx-1)) =         &
        ( real ( 2 * ny - j - 1, kind = 8 ) * yb   &
        + real (          j - 1, kind = 8 ) * yt ) &
        / real ( 2 * ny     - 2, kind = 8 )

    end do
  end do

  return
end
