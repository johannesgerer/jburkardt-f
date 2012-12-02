program main

!*****************************************************************************80
!
!! MAIN is the main program for FEM2D_SAMPLE.
!
!  Discussion:
!
!    FEM2D_SAMPLE reads files defining a 2D FEM representation of data,
!    and a set of sample points, and writes out a file containing the 
!    value of the finite element function at the sample points.
!
!  Usage:
!
!    fem2d_sample fem_prefix sample_prefix
!
!    where 'fem_prefix' is the common prefix for the FEM files:
!
!    * fem_prefix_nodes.txt,    the node coordinates.
!    * fem_prefix_elements.txt, the nodes that make up each element;
!    * fem_prefix_values.txt,   the values defined at each node.
!
!    and 'sample_prefix' is the common prefix for the SAMPLE files.
!    (the node file is input, and the values file is created by the program.)
!
!    * sample_prefix_nodes.txt,  the node coordinates where samples are desired.
!    * sample_prefix_values.txt, the values computed at each sample node.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 June 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer   ( kind = 4 )  fem_element
  character ( len = 255 ) fem_element_filename
  integer   ( kind = 4 ), allocatable, dimension ( :, : ) :: fem_element_neighbor
  integer   ( kind = 4 ), allocatable, dimension ( :, : ) :: fem_element_node
  integer   ( kind = 4 )  fem_element_num
  integer   ( kind = 4 )  fem_element_order
  integer   ( kind = 4 )  fem_node_dim
  integer   ( kind = 4 )  fem_node_num
  character ( len = 255 ) fem_node_filename
  real      ( kind = 8 ), allocatable, dimension ( :, : ) :: fem_node_xy
  character ( len = 255 ) fem_prefix
  real      ( kind = 8 ), allocatable, dimension ( :, : ) :: fem_value
  integer   ( kind = 4 )  fem_value_dim
  character ( len = 255 ) fem_value_filename
  integer   ( kind = 4 )  fem_value_num
  integer   ( kind = 4 )  i
  integer   ( kind = 4 )  iarg
  integer   ( kind = 4 )  iargc
  integer   ( kind = 4 )  ios
  integer   ( kind = 4 )  j
  integer   ( kind = 4 )  num_arg
  integer   ( kind = 4 )  sample_node_dim
  character ( len = 255 ) sample_node_filename
  integer   ( kind = 4 )  sample_node_num
  real      ( kind = 8 ), allocatable :: sample_node_xy(:,:)
  character ( len = 255 ) sample_prefix
  real      ( kind = 8 ), allocatable, dimension ( :, : ) :: sample_value
  integer   ( kind = 4 )  sample_value_dim
  character ( len = 255 ) sample_value_filename
  integer   ( kind = 4 )  sample_value_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FEM2D_SAMPLE'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Read files defining an FEM function of 2 arguments.'
  write ( *, '(a)' ) '  Read a file of sample arguments.'
  write ( *, '(a)' ) '  Write a file of function values at the arguments.'
!
!  Get the number of command line arguments.
!
  num_arg = iargc ( )

  if ( num_arg < 1 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Enter the FEM file prefix:'
    read ( *, '(a)', iostat = ios ) fem_prefix

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FEM2D_SAMPLE - Fatal error!'
      write ( *, '(a)' ) '  Unexpected read error!'
      stop
    end if

  else

    iarg = 1

    call getarg ( iarg, fem_prefix )

  end if

  if ( num_arg < 2 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Enter the sample file prefix:'
    read ( *, '(a)', iostat = ios ) sample_prefix

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FEM2D_SAMPLE - Fatal error!'
      write ( *, '(a)' ) '  Unexpected read error!'
      stop
    end if

  else

    iarg = 2

    call getarg ( iarg, sample_prefix )

  end if
!
!  Create the filenames.
!
  fem_node_filename = trim ( fem_prefix ) // '_nodes.txt'
  fem_element_filename = trim ( fem_prefix ) // '_elements.txt'
  fem_value_filename = trim ( fem_prefix ) // '_values.txt'

  sample_node_filename = trim ( sample_prefix ) // '_nodes.txt'
  sample_value_filename = trim ( sample_prefix ) // '_values.txt'
!
!  Read the FEM data.
!
  call r8mat_header_read ( fem_node_filename, fem_node_dim, fem_node_num )

  allocate ( fem_node_xy(1:fem_node_dim,1:fem_node_num) )

  call r8mat_data_read ( fem_node_filename, fem_node_dim, fem_node_num, &
    fem_node_xy )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The FEM node dimension is        ', fem_node_dim
  write ( *, '(a,i8)' ) '  The FEM node number is           ', fem_node_num

  if ( fem_node_dim /= 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FEM2D_SAMPLE - Fatal error!'
    write ( *, '(a)' ) '  Spatial dimension of the nodes is not 2.'
    stop
  end if

  call i4mat_header_read ( fem_element_filename, fem_element_order, fem_element_num )

  allocate ( fem_element_node(1:fem_element_order,1:fem_element_num) )
  allocate ( fem_element_neighbor(1:3,1:fem_element_num) )

  call i4mat_data_read ( fem_element_filename, fem_element_order, fem_element_num, &
    fem_element_node )

  write ( *, '(a,i8)' ) '  The FEM element order is         ', fem_element_order
  write ( *, '(a,i8)' ) '  The FEM element number is        ', fem_element_num

  call r8mat_header_read ( fem_value_filename, fem_value_dim, fem_value_num )

  write ( *, '(a,i8)' ) '  The FEM value order is           ', fem_value_dim
  write ( *, '(a,i8)' ) '  the FEM value number is          ', fem_value_num

  if ( fem_value_num /= fem_node_num ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FEM2D_SAMPLE - Fatal error!'
    write ( *, '(a)' ) '  Number of FEM values and nodes differ.'
    stop
  end if

  allocate ( fem_value(1:fem_value_dim,1:fem_node_num) )

  call r8mat_data_read ( fem_value_filename, fem_value_dim, fem_value_num, fem_value )
!
!  Create the element neighbor array.
!
  if ( fem_element_order == 3 ) then

    call triangulation_order3_neighbor_triangles ( fem_element_num, &
      fem_element_node, fem_element_neighbor )

  else if ( fem_element_order == 6 ) then

    call triangulation_order6_neighbor_triangles ( fem_element_num, &
      fem_element_node, fem_element_neighbor )

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FEM2D_SAMPLE - Fatal error!'
    write ( *, '(a)' ) '  The element order must be 3 or 6.'
    write ( *, '(a,i8)' ) '  But this data has element order = ', &
      fem_element_order
    return

  end if

  write ( *, '(a)' ) '  The element neighbor array has been computed.'

  if ( .false. ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  ELEMENT NEIGHBOR:'
    write ( *, '(a)' ) ' '
    do fem_element = 1, fem_element_num
      write ( *, '(2x,i8,2x,i8,2x,i8,2x,i8)' ) &
        fem_element, fem_element_neighbor(1:3,fem_element)
    end do
  end if
!
!  Read the SAMPLE node data.
!
  call r8mat_header_read ( sample_node_filename, sample_node_dim, &
    sample_node_num )

  allocate ( sample_node_xy(1:sample_node_dim,1:sample_node_num) )

  call r8mat_data_read ( sample_node_filename, sample_node_dim, &
    sample_node_num, sample_node_xy )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Sample node spatial dimension is ', sample_node_dim
  write ( *, '(a,i8)' ) '  Sample node number is            ', sample_node_num

  if ( sample_node_dim /= 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FEM2D_SAMPLE - Fatal error!'
    write ( *, '(a)' ) '  Spatial dimension of the sample nodes is not 2.'
    stop
  end if
!
!  Compute the sample values.
!
  sample_value_dim = fem_value_dim
  sample_value_num = sample_node_num
  allocate ( sample_value(1:sample_value_dim,1:sample_value_num) )

  call fem2d_evaluate ( fem_node_num, fem_node_xy, fem_element_order, &
    fem_element_num, fem_element_node, fem_element_neighbor, fem_value_dim, &
    fem_value, sample_node_num, sample_node_xy, sample_value )
!
!  Write the sample values.
!
  call r8mat_write ( sample_value_filename, sample_value_dim, &
    sample_value_num, sample_value )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Interpolated FEM data written to "' &
    // trim ( sample_value_filename ) // '".'
!
!  Free memory.
!
  deallocate ( fem_element_neighbor )
  deallocate ( fem_element_node )
  deallocate ( fem_node_xy )
  deallocate ( fem_value )
  deallocate ( sample_node_xy )
  deallocate ( sample_value )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FEM2D_SAMPLE'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
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
subroutine ch_cap ( ch )

!*****************************************************************************80
!
!! CH_CAP capitalizes a single character.
!
!  Discussion:
!
!    Instead of CHAR and ICHAR, we now use the ACHAR and IACHAR functions,
!    which guarantee the ASCII collating sequence.
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

  character              ch
  integer   ( kind = 4 ) itemp

  itemp = iachar ( ch )

  if ( 97 <= itemp .and. itemp <= 122 ) then
    ch = achar ( itemp - 32 )
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
!    CH_EQI ( 'A', 'a' ) is TRUE.
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

  character c1
  character c1_cap
  character c2
  character c2_cap
  logical   ch_eqi

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
subroutine ch_to_digit ( ch, digit )

!*****************************************************************************80
!
!! CH_TO_DIGIT returns the integer ( kind = 4 ) value of a base 10 digit.
!
!  Discussion:
!
!    Instead of ICHAR, we now use the IACHAR function, which
!    guarantees the ASCII collating sequence.
!
!  Example:
!
!     CH  DIGIT
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
!    Input, character CH, the decimal digit, '0' through '9' or blank
!    are legal.
!
!    Output, integer ( kind = 4 ) DIGIT, the corresponding value.  If CH was
!    'illegal', then DIGIT is -1.
!
  implicit none

  character              ch
  integer   ( kind = 4 ) digit

  if ( lle ( '0', ch ) .and. lle ( ch, '9' ) ) then

    digit = iachar ( ch ) - 48

  else if ( ch == ' ' ) then

    digit = 0

  else

    digit = -1

  end if

  return
end
subroutine fem2d_evaluate ( fem_node_num, fem_node_xy, fem_element_order, &
  fem_element_num, fem_element_node, fem_element_neighbor, fem_value_dim, &
  fem_value, sample_node_num, sample_node_xy, sample_value )

!*****************************************************************************80
!
!! FEM2D_EVALUATE samples an FEM function on a T3 or T6 triangulation.
!
!  Discussion:
!
!    Note that the sample values returned are true values of the underlying
!    finite element function.  They are NOT produced by constructing some
!    other function that interpolates the data at the finite element nodes
!    (something which MATLAB's griddata function can easily do.)  Instead, 
!    each sampling node is located within one of the associated finite
!    element triangles, and the finite element function is developed and 
!    evaluated there.  
!
!    MATLAB's scattered data interpolation is wonderful, but it cannot
!    be guaranteed to reproduce the finite element function corresponding
!    to nodal data.  This routine can (or at least tries to!).
! 
!    So if you are using finite elements, then using THIS routine
!    (but not MATLAB's griddata function), what you see is what you have!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 June 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) FEM_NODE_NUM, the number of nodes.
!
!    Input, real ( kind = 8 ) FEM_NODE_XY(2,FEM_NODE_NUM), the coordinates 
!    of the nodes.
!
!    Input, integer ( kind = 4 ) FEM_ELEMENT_ORDER, the order of the elements, 
!    either 3 or 6.
!
!    Input, integer ( kind = 4 ) FEM_ELEMENT_NUM, the number of triangles.
!
!    Input, integer ( kind = 4 ) 
!    FEM_ELEMENT_NODE(FEM_ELEMENT_ORDER,FEM_ELEMENT_NUM), the
!    nodes that make up each triangle.
!
!    Input, integer ( kind = 4 ) FEM_ELEMENT_NEIGHBOR(3,FEM_ELEMENT_NUM), the 
!    index of the neighboring triangle on each side, or -1 if no neighbor there.
!
!    Input, integer ( kind = 4 ) FEM_VALUE_DIM, the "dimension" of the values.
!
!    Input, real ( kind = 8 ) FEM_VALUE(FEM_VALUE_DIM,FEM_NODE_NUM), the 
!    finite element coefficient values at each node.
!
!    Input, integer ( kind = 4 ) SAMPLE_NODE_NUM, the number of sample nodes.
!
!    Input, real ( kind = 8 ) SAMPLE_NODE_XY(2,SAMPLE_NODE_NUM), the sample nodes.
!
!    Output, real ( kind = 8 ) SAMPLE_VALUE(FEM_VALUE_DIM,SAMPLE_NODE_NUM),
!    the sampled values.
!
  implicit none

  integer ( kind = 4 ) fem_element_num
  integer ( kind = 4 ) fem_element_order
  integer ( kind = 4 ) fem_node_num
  integer ( kind = 4 ) fem_value_dim
  integer ( kind = 4 ) sample_node_num

  real    ( kind = 8 ) b(fem_element_order)
  real    ( kind = 8 ) dbdx(fem_element_order)
  real    ( kind = 8 ) dbdy(fem_element_order)
  integer ( kind = 4 ) edge
  integer ( kind = 4 ) fem_element_neighbor(3,fem_element_num)
  integer ( kind = 4 ) fem_element_node(fem_element_order,fem_element_num)
  real    ( kind = 8 ) fem_node_xy(2,fem_node_num)
  real    ( kind = 8 ) fem_value(fem_value_dim,fem_node_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real    ( kind = 8 ), dimension ( 2 ) :: p_xy
  real    ( kind = 8 ) sample_node_xy(2,sample_node_num)
  real    ( kind = 8 ) sample_value(fem_value_dim,sample_node_num)
  integer ( kind = 4 ) t
  integer ( kind = 4 ) t_node(fem_element_order)
  real    ( kind = 8 ) t_xy(2,fem_element_order)
!
!  For each sample point: find the eement T that contains it,
!  and evaluate the finite element function there.
!
  do j = 1, sample_node_num

    p_xy(1:2) = sample_node_xy(1:2,j)
!
!  Find the element T that contains the point.
!
    call triangulation_search_delaunay ( fem_node_num, fem_node_xy, &
      fem_element_order, fem_element_num, fem_element_node, &
      fem_element_neighbor, p_xy, t, edge )
!
!  Evaluate the finite element basis functions at the point in T.
!
    t_node(1:fem_element_order) = fem_element_node(1:fem_element_order,t)

    t_xy(1:2,1:fem_element_order) = fem_node_xy(1:2,t_node)

    if ( fem_element_order == 3 ) then
      call basis_mn_t3 ( t_xy, 1, p_xy, b, dbdx, dbdy )
    else if ( fem_element_order == 6 ) then
      call basis_mn_t6 ( t_xy, 1, p_xy, b, dbdx, dbdy )
    end if
!
!  Multiply by the finite element values to get the sample values.
!
    do i = 1, fem_value_dim
      sample_value(i,j) = dot_product ( &
        fem_value(i,t_node(1:fem_element_order)), b(1:fem_element_order) )
    end do

  end do

  return
end
subroutine file_column_count ( input_filename, column_num )

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
!    Input, character ( len = * ) INPUT_FILENAME, the name of the file.
!
!    Output, integer ( kind = 4 ) COLUMN_NUM, the number of columns in the file.
!
  implicit none

  integer   ( kind = 4 )  column_num
  logical                 got_one
  character ( len = * )   input_filename
  integer   ( kind = 4 )  input_unit
  integer   ( kind = 4 )  ios
  character ( len = 255 ) line
!
!  Open the file.
!
  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_filename, status = 'old', &
    form = 'formatted', access = 'sequential', iostat = ios )

  if ( ios /= 0 ) then
    column_num = -1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_COLUMN_COUNT - Fatal error!'
    write ( *, '(a)' ) '  Could not open the file:'
    write ( *, '(a)' ) '    ' // trim ( input_filename )
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
subroutine file_row_count ( input_filename, row_num )

!*****************************************************************************80
!
!! FILE_ROW_COUNT counts the number of row records in a file.
!
!  Discussion:
!
!    It does not count lines that are blank, or that begin with a
!    comment symbol '#'.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
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
!    Input, character ( len = * ) INPUT_FILENAME, the name of the input file.
!
!    Output, integer ( kind = 4 ) ROW_NUM, the number of rows found.
!
  implicit none

  integer   ( kind = 4 )  bad_num
  integer   ( kind = 4 )  comment_num
  integer   ( kind = 4 )  ierror
  character ( len = * )   input_filename
  integer   ( kind = 4 )  input_unit
  integer   ( kind = 4 )  ios
  character ( len = 255 ) line
  integer   ( kind = 4 )  record_num
  integer   ( kind = 4 )  row_num

  call get_unit ( input_unit )

  open ( unit = input_unit, file = input_filename, status = 'old', &
    iostat = ios )

  if ( ios /= 0 ) then
    row_num = -1;
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FILE_ROW_COUNT - Fatal error!'
    write ( *, '(a)' ) '  Could not open the input file: ' // &
      trim ( input_filename )
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
subroutine get_seed ( seed )

!*****************************************************************************80
!
!! GET_SEED returns a seed for the random number generator.
!
!  Discussion:
!
!    The seed depends on the current time, and ought to be (slightly)
!    different every millisecond.  Once the seed is obtained, a random
!    number generator should be called a few times to further process
!    the seed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) SEED, a pseudorandom seed value.
!
  implicit none

  integer   ( kind = 4 ) seed
  real      ( kind = 8 ) temp
  character ( len = 10 ) time
  character ( len = 8 )  today
  integer   ( kind = 4 ) values(8)
  character ( len = 5 )  zone

  call date_and_time ( today, time, zone, values )

  temp = 0.0D+00

  temp = temp + real ( values(2) - 1, kind = 8 ) /  11.0D+00
  temp = temp + real ( values(3) - 1, kind = 8 ) /  30.0D+00
  temp = temp + real ( values(5),     kind = 8 ) /  23.0D+00
  temp = temp + real ( values(6),     kind = 8 ) /  59.0D+00
  temp = temp + real ( values(7),     kind = 8 ) /  59.0D+00
  temp = temp + real ( values(8),     kind = 8 ) / 999.0D+00
  temp = temp                                    /   6.0D+00

  do while ( temp <= 0.0D+00 )
    temp = temp + 1.0D+00
  end do

  do while ( 1.0D+00 < temp )
    temp = temp - 1.0D+00
  end do

  seed = int ( real ( huge ( 1 ), kind = 8 ) * temp )
!
!  Never use a seed of 0 or maximum integer ( kind = 4 ).
!
  if ( seed == 0 ) then
    seed = 1
  end if

  if ( seed == huge ( 1 ) ) then
    seed = seed - 1
  end if

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
!    Otherwise, IUNIT is a a value between 1 and 99, representing a
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
  logical              lopen

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
!    12 November 2006
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
    seed = seed + 2147483647
  end if

  r = real ( seed, kind = 4 ) * 4.656612875E-10
!
!  Scale R to lie between A-0.5 and B+0.5.
!
  r = ( 1.0E+00 - r ) * ( real ( min ( a, b ), kind = 4 ) - 0.5E+00 ) & 
    +             r   * ( real ( max ( a, b ), kind = 4 ) + 0.5E+00 )
!
!  Use rounding to convert R to an integer ( kind = 4 ) between A and B.
!
  value = nint ( r, kind = 4 )

  value = max ( value, min ( a, b ) )
  value = min ( value, max ( a, b ) )

  i4_uniform = value

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
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
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
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, integer ( kind = 4 ) A(M,N), an array of N columns of vectors
!    of length M.
!
!    Input, integer ( kind = 4 ) I, J, the columns to be compared.
!    I and J must be between 1 and N.
!
!    Output, integer ( kind = 4 ) ISGN, the results of the comparison:
!    -1, column I < column J,
!     0, column I = column J,
!    +1, column J < column I.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
!
!  Check.
!
  if ( i < 1 .or. n < i ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4COL_COMPARE - Fatal error!'
    write ( *, '(a)' ) '  Column index I is out of bounds.'
    stop
  end if

  if ( j < 1 .or. n < j ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4COL_COMPARE - Fatal error!'
    write ( *, '(a)' ) '  Column index J is out of bounds.'
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
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
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
!    Input, integer ( kind = 4 ) M, the number of rows of A, and the length of
!    a vector of data.
!
!    Input, integer ( kind = 4 ) N, the number of columns of A.
!
!    Input/output, integer ( kind = 4 ) A(M,N).
!    On input, the array of N columns of M-vectors.
!    On output, the columns of A have been sorted in ascending
!    lexicographic order.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j

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
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
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
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns 
!    in the array.
!
!    Input/output, integer ( kind = 4 ) A(M,N), an array of N columns 
!    of length M.
!
!    Input, integer ( kind = 4 ) I, J, the columns to be swapped.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) col(m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

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
!    Output, integer ( kind = 4 ) TABLE(M,N), the table data.
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
!    Output, real ( kind = 8 ) TABLE(M,N), the table data.
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
subroutine r8mat_write ( output_filename, m, n, table )

!*****************************************************************************80
!
!! R8MAT_WRITE writes an R8MAT file.
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
!    31 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) OUTPUT_FILENAME, the output file name.
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
  character ( len = * )  output_filename
  integer   ( kind = 4 ) output_status
  integer   ( kind = 4 ) output_unit
  character ( len = 30 ) string
  real      ( kind = 8 ) table(m,n)
!
!  Open the file.
!
  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_filename, &
    status = 'replace', iostat = output_status )

  if ( output_status /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_WRITE - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the output file "' // &
      trim ( output_filename ) // '" on unit ', output_unit
    output_unit = -1
    stop
  end if
!
!  Create a format string.
!
!  For less precision in the output file, try:
!
!                                            '(', m, 'g', 14, '.', 6, ')'
!
  if ( 0 < m .and. 0 < n ) then

    write ( string, '(a1,i8,a1,i8,a1,i8,a1)' ) '(', m, 'g', 24, '.', 16, ')'
!
!  Write the data.
!
    do j = 1, n
      write ( output_unit, string ) table(1:m,j)
    end do

  end if
!
!  Close the file.
!
  close ( unit = output_unit )

  return
end
subroutine s_to_i4 ( s, ival, ierror, length )

!*****************************************************************************80
!
!! S_TO_I4 reads an integer ( kind = 4 ) value from a string.
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
!    Output, integer ( kind = 4 ) IVAL, the value read from the string.
!    If the string is blank, then IVAL will be returned 0.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!    0, no error.
!    1, an error occurred.
!
!    Output, integer ( kind = 4 ) LENGTH, the number of characters 
!    used to make IVAL.
!
  implicit none

  character              c
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) ierror
  integer   ( kind = 4 ) isgn
  integer   ( kind = 4 ) istate
  integer   ( kind = 4 ) ival
  integer   ( kind = 4 ) length
  character ( len = * )  s

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
!! S_TO_I4VEC reads an integer ( kind = 4 ) vector from a string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
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
!    Input, integer ( kind = 4 ) N, the number of values expected.
!
!    Output, integer ( kind = 4 ) IVEC(N), the values read from the string.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no errors occurred.
!    -K, could not read data for entries -K through N.
!
  implicit none

  integer   ( kind = 4 ) n

  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) ierror
  integer   ( kind = 4 ) ilo
  integer   ( kind = 4 ) ivec(n)
  integer   ( kind = 4 ) length
  character ( len = * )  s

  i = 0

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
!       3 integer ( kind = 4 ) part,
!       4 decimal point,
!       5 fraction part,
!       6 'E' or 'e' or 'D' or 'd', exponent marker,
!       7 exponent sign,
!       8 exponent integer ( kind = 4 ) part,
!       9 exponent decimal point,
!      10 exponent fraction part,
!      11 blanks,
!      12 final comma or semicolon,
!
!    with most quantities optional.
!
!  Example:
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
!    Input, character ( len = * ) S, the string containing the
!    data to be read.  Reading will begin at position 1 and
!    terminate at the end of the string, or when no more
!    characters can be read to form a legal real.  Blanks,
!    commas, or other nonnumeric data will, in particular,
!    cause the conversion to halt.
!
!    Output, real ( kind = 8 ) DVAL, the value read from the string.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no errors occurred.
!    1, 2, 6 or 7, the input number was garbled.  The
!    value of IERROR is the last type of input successfully
!    read.  For instance, 1 means initial blanks, 2 means
!    a plus or minus sign, and so on.
!
!    Output, integer ( kind = 4 ) LENGTH, the number of characters read
!    to form the number, including any terminating
!    characters such as a trailing comma or blanks.
!
  implicit none

  logical                ch_eqi
  character              c
  real      ( kind = 8 ) dval
  integer   ( kind = 4 ) ierror
  integer   ( kind = 4 ) ihave
  integer   ( kind = 4 ) isgn
  integer   ( kind = 4 ) iterm
  integer   ( kind = 4 ) jbot
  integer   ( kind = 4 ) jsgn
  integer   ( kind = 4 ) jtop
  integer   ( kind = 4 ) length
  integer   ( kind = 4 ) nchar
  integer   ( kind = 4 ) ndig
  real      ( kind = 8 ) rbot
  real      ( kind = 8 ) rexp
  real      ( kind = 8 ) rtop
  character ( len = * )  s

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

  integer   ( kind = 4 ) n

  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) ierror
  integer   ( kind = 4 ) ilo
  integer   ( kind = 4 ) lchar
  real      ( kind = 8 ) rvec(n)
  character ( len = * )  s

  i = 0

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

  logical                blank
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) lens
  integer   ( kind = 4 ) nword
  character ( len = * )  s

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
!    routine may be used to sort integer ( kind = 4 )s, reals, numbers, names,
!    dates, shoe sizes, and so on.  After each call, the routine asks
!    the user to compare or interchange two items, until a special
!    return value signals that the sorting is completed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 February 2004
!
!  Author:
!
!    Original FORTRAN77 version by Nijenhuis, WIlf.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items to be sorted.
!
!    Input/output, integer ( kind = 4 ) INDX, the main communication signal.
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
!    Output, integer ( kind = 4 ) I, J, the indices of two items.
!    On return with INDX positive, elements I and J should be interchanged.
!    On return with INDX negative, elements I and J should be compared, and
!    the result reported in ISGN on the next call.
!
!    Input, integer ( kind = 4 ) ISGN, results of comparison of elements I and J.
!    (Used only when the previous call returned INDX less than 0).
!    ISGN <= 0 means I is less than or equal to J;
!    0 <= ISGN means I is greater than or equal to J.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ), save :: i_save = 0
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  integer ( kind = 4 ), save :: j_save = 0
  integer ( kind = 4 ), save :: k = 0
  integer ( kind = 4 ), save :: k1 = 0
  integer ( kind = 4 ) n
  integer ( kind = 4 ), save :: n1 = 0
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

  character ( len = 8 ) ampm
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

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
subroutine triangulation_order3_neighbor_triangles ( triangle_num, &
  triangle_node, triangle_neighbor )

!*****************************************************************************80
!
!! TRIANGULATION_ORDER3_NEIGHBOR_TRIANGLES determines triangle neighbors.
!
!  Discussion:
!
!    A triangulation of a set of nodes can be completely described by
!    the coordinates of the nodes, and the list of nodes that make up
!    each triangle.  However, in some cases, it is necessary to know
!    triangle adjacency information, that is, which triangle, if any,
!    is adjacent to a given triangle on a particular side.
!
!    This routine creates a data structure recording this information.
!
!    The primary amount of work occurs in sorting a list of 3 * TRIANGLE_NUM
!    data items.
!
!    Note that ROW is a work array allocated dynamically inside this
!    routine.  It is possible, for very large values of TRIANGLE_NUM,
!    that the necessary amount of memory will not be accessible, and the
!    routine will fail.  This is a limitation of the implementation of
!    dynamic arrays in FORTRAN90.  One way to get around this would be
!    to require the user to declare ROW in the calling routine
!    as an allocatable array, get the necessary memory explicitly with
!    an ALLOCATE statement, and then pass ROW into this routine.
!
!    Of course, the point of dynamic arrays was to make it easy to
!    hide these sorts of temporary work arrays from the poor user!
!
!    This routine was revised to store the edge data in a column
!    array rather than a row array.
!
!  Example:
!
!    The input information from TRIANGLE_NODE:
!
!    Triangle   Nodes
!    --------   ---------------
!     1         3      4      1
!     2         3      1      2
!     3         3      2      8
!     4         2      1      5
!     5         8      2     13
!     6         8     13      9
!     7         3      8      9
!     8        13      2      5
!     9         9     13      7
!    10         7     13      5
!    11         6      7      5
!    12         9      7      6
!    13        10      9      6
!    14         6      5     12
!    15        11      6     12
!    16        10      6     11
!
!    The output information in TRIANGLE_NEIGHBOR:
!
!    Triangle  Neighboring Triangles
!    --------  ---------------------
!
!     1        -1     -1      2
!     2         1      4      3
!     3         2      5      7
!     4         2     -1      8
!     5         3      8      6
!     6         5      9      7
!     7         3      6     -1
!     8         5      4     10
!     9         6     10     12
!    10         9      8     11
!    11        12     10     14
!    12         9     11     13
!    13        -1     12     16
!    14        11     -1     15
!    15        16     14     -1
!    16        13     15     -1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) TRIANGLE_NUM, the number of triangles.
!
!    Input, integer ( kind = 4 ) TRIANGLE_NODE(3,TRIANGLE_NUM), the nodes 
!    that make up each triangle.
!
!    Output, integer ( kind = 4 ) TRIANGLE_NEIGHBOR(3,TRIANGLE_NUM), the three
!    triangles that are direct neighbors of a given triangle.  
!    TRIANGLE_NEIGHBOR(1,I) is the index of the triangle which touches side 1, 
!    defined by nodes 2 and 3, and so on.  TRIANGLE_NEIGHBOR(1,I) is negative 
!    if there is no neighbor on that side.  In this case, that side of the 
!    triangle lies on the boundary of the triangulation.
!
  implicit none

  integer ( kind = 4 ) triangle_num
  integer ( kind = 4 ), parameter :: triangle_order = 3

  integer ( kind = 4 ) col(4,3*triangle_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icol
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) side1
  integer ( kind = 4 ) side2
  integer ( kind = 4 ) triangle_neighbor(3,triangle_num)
  integer ( kind = 4 ) tri
  integer ( kind = 4 ) triangle_node(triangle_order,triangle_num)
  integer ( kind = 4 ) tri1
  integer ( kind = 4 ) tri2
!
!  Step 1.
!  From the list of nodes for triangle T, of the form: (I,J,K)
!  construct the three neighbor relations:
!
!    (I,J,1,T) or (J,I,1,T),
!    (J,K,2,T) or (K,J,2,T),
!    (K,I,3,T) or (I,K,3,T)
!
!  where we choose (I,J,1,T) if I < J, or else (J,I,1,T)
!
  do tri = 1, triangle_num

    i = triangle_node(1,tri)
    j = triangle_node(2,tri)
    k = triangle_node(3,tri)

    if ( i < j ) then
      col(1:4,3*(tri-1)+1) = (/ i, j, 1, tri /)
    else
      col(1:4,3*(tri-1)+1) = (/ j, i, 1, tri /)
    end if

    if ( j < k ) then
      col(1:4,3*(tri-1)+2) = (/ j, k, 2, tri /)
    else
      col(1:4,3*(tri-1)+2) = (/ k, j, 2, tri /)
    end if

    if ( k < i ) then
      col(1:4,3*(tri-1)+3) = (/ k, i, 3, tri /)
    else
      col(1:4,3*(tri-1)+3) = (/ i, k, 3, tri /)
    end if

  end do
!
!  Step 2. Perform an ascending dictionary sort on the neighbor relations.
!  We only intend to sort on rows 1 and 2; the routine we call here
!  sorts on rows 1 through 4 but that won't hurt us.
!
!  What we need is to find cases where two triangles share an edge.
!  Say they share an edge defined by the nodes I and J.  Then there are
!  two columns of COL that start out ( I, J, ?, ? ).  By sorting COL,
!  we make sure that these two columns occur consecutively.  That will
!  make it easy to notice that the triangles are neighbors.
!
  call i4col_sort_a ( 4, 3*triangle_num, col )
!
!  Step 3. Neighboring triangles show up as consecutive columns with
!  identical first two entries.  Whenever you spot this happening,
!  make the appropriate entries in TRIANGLE_NEIGHBOR.
!
  triangle_neighbor(1:3,1:triangle_num) = -1

  icol = 1

  do

    if ( 3 * triangle_num <= icol ) then
      exit
    end if

    if ( col(1,icol) /= col(1,icol+1) .or. col(2,icol) /= col(2,icol+1) ) then
      icol = icol + 1
      cycle
    end if

    side1 = col(3,icol)
    tri1 = col(4,icol)
    side2 = col(3,icol+1)
    tri2 = col(4,icol+1)

    triangle_neighbor(side1,tri1) = tri2
    triangle_neighbor(side2,tri2) = tri1

    icol = icol + 2

  end do

  return
end
subroutine triangulation_order6_neighbor_triangles ( triangle_num, &
  triangle_node, triangle_neighbor )

!*****************************************************************************80
!
!! TRIANGULATION_ORDER6_NEIGHBOR_TRIANGLES determines triangle neighbors.
!
!  Discussion:
!
!    A triangulation of a set of nodes can be completely described by
!    the coordinates of the nodes, and the list of nodes that make up
!    each triangle.  However, in some cases, it is necessary to know
!    triangle adjacency information, that is, which triangle, if any,
!    is adjacent to a given triangle on a particular side.
!
!    This routine creates a data structure recording this information.
!
!    The primary amount of work occurs in sorting a list of 3 * TRIANGLE_NUM
!    data items.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) TRIANGLE_NUM, the number of triangles.
!
!    Input, integer ( kind = 4 ) TRIANGLE_ORDER(6,TRIANGLE_NUM), the nodes 
!    that make up each triangle.
!
!    Output, integer ( kind = 4 ) TRIANGLE_NEIGHBOR(3,TRIANGLE_NUM), the three
!    triangles that are direct neighbors of a given triangle.  
!    TRIANGLE_NEIGHBOR(1,I) is the index of the triangle which touches side 1,
!    defined by nodes 2 and 3, and so on.  TRIANGLE_NEIGHBOR(1,I) is negative 
!    if there is no neighbor on that side.  In this case, that side of the 
!    triangle lies on the boundary of the triangulation.
!
  implicit none

  integer ( kind = 4 ) triangle_num
  integer ( kind = 4 ), parameter :: triangle_order = 6

  integer ( kind = 4 ) col(4,3*triangle_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icol
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) side1
  integer ( kind = 4 ) side2
  integer ( kind = 4 ) tri
  integer ( kind = 4 ) triangle_node(triangle_order,triangle_num)
  integer ( kind = 4 ) triangle_neighbor(3,triangle_num)
  integer ( kind = 4 ) tri1
  integer ( kind = 4 ) tri2
!
!  Step 1.
!  From the list of vertices for triangle T, of the form: (I,J,K)
!  construct the three neighbor relations:
!
!    (I,J,1,T) or (J,I,1,T),
!    (J,K,2,T) or (K,J,2,T),
!    (K,I,3,T) or (I,K,3,T)
!
!  where we choose (I,J,1,T) if I < J, or else (J,I,1,T)
!
  do tri = 1, triangle_num

    i = triangle_node(1,tri)
    j = triangle_node(2,tri)
    k = triangle_node(3,tri)

    if ( i < j ) then
      col(1:4,3*(tri-1)+1) = (/ i, j, 1, tri /)
    else
      col(1:4,3*(tri-1)+1) = (/ j, i, 1, tri /)
    end if

    if ( j < k ) then
      col(1:4,3*(tri-1)+2) = (/ j, k, 2, tri /)
    else
      col(1:4,3*(tri-1)+2) = (/ k, j, 2, tri /)
    end if

    if ( k < i ) then
      col(1:4,3*(tri-1)+3) = (/ k, i, 3, tri /)
    else
      col(1:4,3*(tri-1)+3) = (/ i, k, 3, tri /)
    end if

  end do
!
!  Step 2. Perform an ascending dictionary sort on the neighbor relations.
!  We only intend to sort on columns 1 and 2; the routine we call here
!  sorts on columns 1 through 4 but that won't hurt us.
!
!  What we need is to find cases where two triangles share an edge.
!  Say they share an edge defined by the nodes I and J.  Then there are
!  two colums of COL that start out ( I, J, ?, ? ).  By sorting COL,
!  we make sure that these two columns occur consecutively.  That will
!  make it easy to notice that the triangles are neighbors.
!
  call i4col_sort_a ( 3*triangle_num, 4, col )
!
!  Step 3. Neighboring triangles show up as consecutive columns with
!  identical first two entries.  Whenever you spot this happening,
!  make the appropriate entries in TRIANGLE_NEIGHBOR.
!
  triangle_neighbor(1:3,1:triangle_num) = -1

  icol = 1

  do

    if ( 3 * triangle_num <= icol ) then
      exit
    end if

    if ( col(1,icol) /= col(1,icol+1) .or. col(2,icol) /= col(2,icol+1) ) then
      icol = icol + 1
      cycle
    end if

    side1 = col(3,icol)
    tri1 = col(4,icol)
    side2 = col(3,icol+1)
    tri2 = col(4,icol+1)

    triangle_neighbor(side1,tri1) = tri2
    triangle_neighbor(side2,tri2) = tri1

    icol = icol + 2

  end do

  return
end
subroutine triangulation_search_delaunay ( node_num, node_xy, triangle_order, &
  triangle_num, triangle_node, triangle_neighbor, p, triangle_index, edge )

!*****************************************************************************80
!
!! TRIANGULATION_SEARCH_DELAUNAY searches a Delaunay triangulation for a point.
!
!  Discussion:
!
!    The algorithm "walks" from one triangle to its neighboring triangle,
!    and so on, until a triangle is found containing point P, or P is found
!    to be outside the convex hull.
!
!    The algorithm computes the barycentric coordinates of the point with
!    respect to the current triangle.  If all three quantities are positive,
!    the point is contained in the triangle.  If the I-th coordinate is
!    negative, then P lies on the far side of edge I, which is opposite
!    from vertex I.  This gives a hint as to where to search next.
!
!    For a Delaunay triangulation, the search is guaranteed to terminate.
!    For other triangulations, a cycle may occur.
!
!    Note the surprising fact that, even for a Delaunay triangulation of
!    a set of nodes, the nearest node to P need not be one of the
!    vertices of the triangle containing P.
!
!    The code can be called for triangulations of any order, but only
!    the first three nodes in each triangle are considered.  Thus, if
!    higher order triangles are used, and the extra nodes are intended
!    to give the triangle a polygonal shape, these will have no effect,
!    and the results obtained here might be misleading.
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
!    Input, real ( kind = 8 ) NODE_XY(2,NODE_NUM), the coordinates of the nodes.
!
!    Input, integer ( kind = 4 ) TRIANGLE_ORDER, the order of the triangles.
!
!    Input, integer ( kind = 4 ) TRIANGLE_NUM, the number of triangles.
!
!    Input, integer ( kind = 4 ) TRIANGLE_NODE(TRIANGLE_ORDER,TRIANGLE_NUM),
!    the nodes that make up each triangle.
!
!    Input, integer ( kind = 4 ) TRIANGLE_NEIGHBOR(3,TRIANGLE_NUM), the 
!    triangle neighbor list.
!
!    Input, real ( kind = 8 ) P(2), the coordinates of a point.
!
!    Output, integer ( kind = 4 ) TRIANGLE_INDEX, the index of the triangle 
!    where the search ended.  If a cycle occurred, then TRIANGLE_INDEX = -1.
!
!    Output, integer ( kind = 4 ) EDGE, indicates the position of the point P in
!    triangle TRIANGLE_INDEX:
!    0, the interior or boundary of the triangle;
!    -1, outside the convex hull of the triangulation, past edge 1;
!    -2, outside the convex hull of the triangulation, past edge 2;
!    -3, outside the convex hull of the triangulation, past edge 3.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) triangle_num
  integer ( kind = 4 ) triangle_order

  integer ( kind = 4 ) a
  real    ( kind = 8 ) alpha
  integer ( kind = 4 ) b
  real    ( kind = 8 ) beta
  integer ( kind = 4 ) c
  integer ( kind = 4 ) count
  real    ( kind = 8 ) det
  real    ( kind = 8 ) dxp
  real    ( kind = 8 ) dxa
  real    ( kind = 8 ) dxb
  real    ( kind = 8 ) dyp
  real    ( kind = 8 ) dya
  real    ( kind = 8 ) dyb
  integer ( kind = 4 ) edge
  real    ( kind = 8 ) gamma
  integer ( kind = 4 ) i4_uniform
  real    ( kind = 8 ) node_xy(dim_num,node_num)
  real    ( kind = 8 ) p(dim_num)
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) triangle_node(triangle_order,triangle_num)
  integer ( kind = 4 ) triangle_index
  integer ( kind = 4 ) triangle_neighbor(3,triangle_num)

  count = 0
  edge = 0

  call get_seed ( seed )

  triangle_index = i4_uniform ( 1, triangle_num, seed )

  do

    count = count + 1

    if ( triangle_num < count ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TRIANGULATION_SEARCH_DELAUNAY - Fatal error!'
      write ( *, '(a)' ) '  The algorithm seems to be cycling.'
      triangle_index = -1
      edge = -1
      stop
    end if
!
!  Get the nodes of triangle TRIANGLE_INDEX.
!
    a = triangle_node(1,triangle_index)
    b = triangle_node(2,triangle_index)
    c = triangle_node(3,triangle_index)
!
!  Using vertex C as a base, compute the distances to vertices A and B,
!  and the point P.
!
    dxa = node_xy(1,a) - node_xy(1,c)
    dya = node_xy(2,a) - node_xy(2,c)

    dxb = node_xy(1,b) - node_xy(1,c)
    dyb = node_xy(2,b) - node_xy(2,c)

    dxp = p(1)         - node_xy(1,c)
    dyp = p(2)         - node_xy(2,c)

    det = dxa * dyb - dya * dxb
!
!  Compute the barycentric coordinates of the point P with respect
!  to this triangle.
!
    alpha = ( dxp * dyb - dyp * dxb ) / det
    beta =  ( dxa * dyp - dya * dxp ) / det
    gamma = 1.0D+00 - alpha - beta
!
!  If the barycentric coordinates are all positive, then the point
!  is inside the triangle and we're done.
!
    if ( 0.0D+00 <= alpha .and. &
         0.0D+00 <= beta  .and. &
         0.0D+00 <= gamma ) then
      exit
    end if
!
!  At least one barycentric coordinate is negative.
!
!  If there is a negative barycentric coordinate for which there exists
!  an opposing triangle neighbor closer to the point, move to that triangle.
!
!  (Two coordinates could be negative, in which case we could go for the
!  most negative one, or the most negative one normalized by the actual
!  distance it represents).
!
    if ( alpha < 0.0D+00 .and. 0 < triangle_neighbor(2,triangle_index) ) then
      triangle_index = triangle_neighbor(2,triangle_index)
      cycle
    else if ( beta < 0.0D+00 .and. &
      0 < triangle_neighbor(3,triangle_index) ) then
      triangle_index = triangle_neighbor(3,triangle_index)
      cycle
    else if ( gamma < 0.0D+00 .and. &
      0 < triangle_neighbor(1,triangle_index) ) then
      triangle_index = triangle_neighbor(1,triangle_index)
      cycle
    end if
!
!  All negative barycentric coordinates correspond to vertices opposite
!  sides on the convex hull.
!
!  Note the edge and exit.
!
    if ( alpha < 0.0D+00 ) then
      edge = -2
      exit
    else if ( beta < 0.0D+00 ) then
      edge = -3
      exit
    else if ( gamma < 0.0D+00 ) then
      edge = -1
      exit
    end if

  end do

  return
end
