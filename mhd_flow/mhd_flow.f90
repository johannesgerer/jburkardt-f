program main

!*****************************************************************************80
!
!! MAIN is the main program for MHD_FLOW.
!
!  Discussion:
!
!    MHD_FLOW solves a magnetic hydrodynamic system.
!
!    This program solves an electrically conducting fluid flow problem
!    on a rectangular domain.  The equations are time independent
!    and a backward Euler approximation is used.
!
!    The MHD problem is formulated in terms of the primitive flow
!    variables U, V, P and and the magnetic field B = (B1,B2)
!
!    There is a specified parabolic inflow at bottom portion (lower 1/4)
!    of the left boundary and outflow at portion (upper 1/4) of right
!    boundary where grad U dot n = 0.   Consequently, choose the number
!    of intervals in the Y directions to be a multiple of 4.
!
!    This version uses finite element techniques with piecewise linear
!    functions on triangles to approximate the pressure and quadratics
!    on triangles for the velocity (Taylor-Hood element);  quadratics
!    are used for the magnetic field.
!
!    The term
!
!      Integral ( div B * div test function )
!
!    is added to the equation for B since otherwise B is not guaranteed to
!    be exactly weakly divergence-free.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 October 2006
!
!  Author:
!
!    Janet Peterson
!
!  Reference:
!
!    Max Gunzburger, Catalin Trenchea,
!    Analysis and Discretization of an Optimal Control Problem
!    for the Time-Periodic MHD Equations,
!    Journal of Mathematical Analysis and Applications,
!    Volume 308, Number 2, 2005, pages 440-466.
!
!  Parameters:
!
!    Local, real ( kind = 8 ) A(UNKNOWN_NUM,2*N_LBAND+1), the matrix used
!    in the linear system defined by Newton's method, in solving for the
!    corrections to the finite element coefficients.
!
!    Local, real ( kind = 8 ) AREA(ELEMENT_NUM), the area of each
!    element.
!
!    Local, real ( kind = 8 ) CP, the coupling parameter.
!
!    Local, real ( kind = 8 ) DT, the time step.
!
!    Local, integer ELEMENT_MAX, the maximum number of elements.
!
!    Local, integer ELEMENT_NUM, the number of elements.
!
!    Input, integer INDX(5,NODE_NUM), at each node, lists the
!    index of the finite element coefficients for U, V, B1, B2, and P.
!
!    Input, integer ELEMENT_NODE(6,ELEMENT_NUM), the nodes of each element.
!
!    Local, integer NODE_NUM, the number of nodes.
!
!    Local, integer NX, controls the number of elements in the X direction.
!
!    Local, integer NY, controls the number of elements in the Y direction.
!
!    Local, integer QUAD_NUM, the order of the quadrature rule.
!
!    Local, real ( kind = 8 ) REY_MAG, the magnetic Reynolds number.
!
!    Local, real ( kind = 8 ) REYNOLD, the Reynolds number.
!
!    Output, integer UNKNOWN_NUM, the number of unknowns.
!
!    Local, real ( kind = 8 ) XC(NODE_NUM), the X coordinates
!    of nodes.
!
!    Output, real ( kind = 8 ) XQ(ELEMENT_MAX,3), the X coordinates of
!    the quadrature points.
!
!    Input, integer Y_IN_NODE, nodes 1 through Y_IN_NODE are on the
!    inflow boundary.
!
!    Input, real ( kind = 8 ) Y_INFLOW, the maximum Y coordinate at which
!    there is inflow.
!
!    Input, integer Y_OUT_NODE, nodes Y_OUT_NODE through NODE_NUM
!    are on the outflow boundary.
!
!    Input, real ( kind = 8 ) Y_OUTFLOW, the minimum Y coordinate at which
!    there is outflow.
!
!    Local, real ( kind = 8 ) YC(NODE_NUM), the Y coordinates
!    of nodes.
!
!    Output, real ( kind = 8 ) YQ(ELEMENT_MAX,3), the Y coordinates of
!    the quadrature points.
!
!    Input, real ( kind = 8 ) YLENGTH, the maximum Y coordinate at which
!    there is outflow .
!
  implicit none

  integer, parameter :: nx = 9
  integer, parameter :: ny = 9

  integer, parameter :: element_max = 2 * ( nx - 1 ) * ( ny - 1 )
  integer, parameter :: quad_num = 3

  real ( kind = 8 ), allocatable, dimension (:,:) :: a
  real ( kind = 8 ), allocatable, dimension (:) :: area
  real ( kind = 8 ), parameter :: cp = 1.0D+00
  real ( kind = 8 )  :: dt = 0.025D+00
  character ( len = 80 ) :: element_file_name = 'elements.txt'
  integer, allocatable, dimension(:,:) :: element_node
  integer :: element_num
  real ( kind = 8 ), allocatable, dimension (:) :: f_new
  real ( kind = 8 ), allocatable, dimension (:) :: f_prev
  integer, parameter :: flag_write = 2
  integer :: i
  integer, allocatable, dimension ( :, : ) :: indx
  character ( len = 20 ) :: mag_file_name = 'mag000.txt'
  integer, parameter :: max_newton = 12
  integer :: n_lband
  integer, parameter  :: n_simple = 2
  integer :: node_num
  character ( len = 20 ) :: p_file_name = 'p000.txt'
  integer :: plot_steps = 5
  real ( kind = 8 ), parameter :: rey_mag = 1.0D+00
  real ( kind = 8 ), parameter :: reynold = 1.0D+00
  real ( kind = 8 ) :: time_cur
  integer, parameter :: time_step_max = 40
  integer :: time_step_num
  real ( kind = 8 ), parameter :: tol_newton = 0.5D-08
  integer :: unknown_num
  character ( len = 20 ) :: uv_file_name = 'uv000.txt'
  real ( kind = 8 ) :: visc
  real ( kind = 8 ), allocatable, dimension(:) :: xc
  real ( kind = 8 ), parameter :: xlength = 1.0D+00
  real ( kind = 8 ), dimension(element_max,quad_num) :: xq
  character ( len = 20 ) :: xy_file_name = 'xy.txt'
  real ( kind = 8 ) :: y_inflow
  real ( kind = 8 ) :: y_outflow
  integer :: y_in_node
  integer :: y_out_node
  real ( kind = 8 ), allocatable, dimension(:) :: yc
  real ( kind = 8 ), dimension(element_max,quad_num) :: yq
  real ( kind = 8 ), parameter :: ylength = 1.0D+00

  visc = 1.0D+00 / reynold

  write ( *, '(a)' ) ' '
  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MHD_FLOW:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Solution of the incompressible time-dependent'
  write ( *, '(a)' ) '  2D Navier-Stokes equations, coupled with a '
  write ( *, '(a)' ) '  magnetic field.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6,a)' ) '  CP =            ', cp, &
    '  Coupling'
  write ( *, '(a,g14.6,a)' ) '  DT =            ', dt, &
    '  Initial time step'
  write ( *, '(a,i14,a)'   ) '  NX =            ', nx, &
    '  Elements/width'
  write ( *, '(a,i14,a)'   ) '  NY =            ', ny, &
    '  Elements/height'
  write ( *, '(a,g14.6,a)' ) '  REY_MAG =       ', rey_mag, &
    '  Reynolds number'
  write ( *, '(a,g14.6,a)' ) '  REYNOLD =       ', reynold, &
    '  Reynolds number'
  write ( *, '(a,i14,a)'   ) '  TIME_STEP_MAX = ', time_step_max, &
    '  Number of time steps'
  write ( *, '(a,g14.6,a)' ) '  XLENGTH =       ', xlength, &
    '  Domain width'
  write ( *, '(a,g14.6,a)' ) '  YLENGTH =       ', ylength, &
    '  Domain height'
!
!  Node coordinates.
!
  node_num = ( 2 * nx - 1 ) * ( 2 * ny - 1 )
  write ( *, '(a,i14,a)' ) '  NODE_NUM =      ', node_num

  allocate ( indx(5,node_num) )
  allocate ( xc(1:node_num) )
  allocate ( yc(1:node_num) )

  call xy_set ( node_num, nx, ny, xlength, ylength, xc, yc )

  call xy_write ( xy_file_name, node_num, xc, yc )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Node coordinates written to "' &
    // trim ( xy_file_name ) // '".'
!
!  Elements.
!
  element_num = 2 * ( nx - 1 ) * ( ny - 1 )
  write ( *, '(a,i14,a)' ) '  ELEMENT_NUM =   ', element_num

  allocate ( element_node(6,element_num) )

  call element_set ( nx, ny, element_num, element_node )

  call element_write ( element_file_name, element_num, element_node )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Elements written to "' &
    // trim ( element_file_name ) // '".'

  allocate ( area(element_num) )

  call element_area ( element_num, element_node, node_num, xc, yc, area )

  call quad_rule ( element_max, element_num, element_node, node_num, &
    xc, yc, xq, yq )
!
!  Variables.
!
  i = ( ny + ny - 2  ) / 4
  y_in_node = i + 1
  y_out_node = ( 2 * nx - 1 ) * ( 2 * ny - 1 ) - i

  call indx_set ( node_num, nx, ny, y_in_node, y_out_node, indx, unknown_num )
!
!  Bandwidth.
!
  call bandwidth ( element_num, node_num, element_node, indx, n_lband )
!
!  Set position for inflow and outflow boundary locations in y-direction
!  set to 1/4 * ylength at bottom for inflow and 1/4 from top at outflow
!
  y_inflow = yc(y_in_node)
  y_outflow = yc(y_out_node )

  write(*,*)  ' inflow boundary conditions:  u  quadratic from 0 to ', y_inflow
  write(*,*)  ' outflow boundary conditions: grad u dot n=0  from ', y_outflow, ' to ', ylength
  write ( *, '(a)' )  '                   U = 0 elsewhere  '
  write ( *, '(a)' )  '                   V = 0 everywhere '
  write ( *, '(a)' )  '             B dot N = 0 (essential) '
  write ( *, '(a)' )  '      curl B cross N = 0 (natural) '
!
!  Write out geometry information.
!
  write ( *, '(a)' ) ' '
  write ( *, * ) '  Number of nodes is    ', node_num
  write ( *, * ) '  Number of elements is ', element_num
  write ( *, * ) '  Number of unknowns is ', unknown_num
  write ( *, * ) '  Bandwidth is          ', 2 * n_lband + 1
!
!  Allocate arrays for assembly and solver
!
  allocate ( a(unknown_num,2*n_lband+1) )
  allocate ( f_new(unknown_num) )
  allocate ( f_prev(unknown_num) )
!
!  Set initial condition for time iteration.
!
  f_prev(1:unknown_num) = 0.0D+00
!
!  Time stepping loop
!
  time_cur = 0.0D+00

  do time_step_num = 1, time_step_max

    time_cur = time_cur + dt

    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  T = ', time_cur
!
!  Solve the MHD/Navier-Stokes equations at the new time using Newton's method.
!
    call newton ( area, cp, dt, f_prev, indx, element_max, max_newton, &
      element_num,  n_lband, node_num, element_node, &
      quad_num, n_simple, unknown_num, rey_mag, time_cur, &
      tol_newton, visc, xc, xq, yc, yq, ylength, y_inflow, y_outflow,  &
      a, f_new )
!
!  Write data to files.
!
    if ( mod ( time_step_num, plot_steps ) == 0 ) then

      call file_name_inc ( uv_file_name )

      call uv_write ( f_new, indx, node_num, uv_file_name, &
        xc, yc, ylength, y_inflow, y_outflow  )

      call file_name_inc ( mag_file_name )

      call mag_write ( f_new, indx, mag_file_name, node_num )

      call file_name_inc ( p_file_name )

      call p_write ( f_new, indx, p_file_name, element_num, &
        node_num, element_node )

    end if
!
!  Get ready for next time step.
!
    f_prev(1:unknown_num) = f_new(1:unknown_num)

  end do
!
!  Deallocate arrays.
!
  deallocate ( a )
  deallocate ( area )
  deallocate ( element_node )
  deallocate ( f_new )
  deallocate ( f_prev )
  deallocate ( indx )
  deallocate ( xc )
  deallocate ( yc )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MHD_FLOW:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine banded ( a, b, n, m1, m2, maxu )

!*****************************************************************************80
!
!! BANDED solves a banded linear system.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 August 2005
!
!  Inputs
!    a    = banded matrix, with non-zero bands stored as columns;
!              the (i,j) entry in dense mode is stored in
!              the (i,m1+1+j-i) entry of banded storage, where
!              m1=lower half-bandwidth; thus the main diagonal
!              is stored in the column number (m1+1)
!    b    = right-hand-side vector of linear system
!    n    = size of linear system to be solved
!    m1   = lower half-bandwidth
!    m2   = upper half-bandwidth; number of nonzero diagonals is
!              thus (m1+m2+1)
!    maxu = first dimension of a(.,.) in calling program
!
!  Output
!    b = solution vector for linear system
!
  implicit none

  integer :: maxu     ! leading dimension of a
  integer :: n          ! number of equations
  real ( kind = 8 ), dimension (maxu, *) :: a
!         banded matrix, with non-zero bands stored as columns;
!              the (i,j) entry in dense mode is stored in
!              the (i,m1+1+j-i) entry of banded storage, where
!              m1=lower half-bandwidth; thus the main diagonal
!              is stored in the column number (m1+1)
  real ( kind = 8 ), dimension(*) :: b
  real ( kind = 8 ), dimension(n) :: wk

  integer, dimension (n) :: is

  integer :: m1        ! lower bandwidth
  integer :: m2        ! upper bandwidth
!
!  Local variables
  real ( kind = 8 ) :: p, q,  tols, x
  integer :: i, ii, im1, j, je, js,  k, kp1, l, m, nm1, np1

  tols=0.00001D+00
  m=m1+m2+1

  do   i=1,n

    p=abs(a(i,1))

    do  j=2,m

      q=abs(a(i,j))
      if(q.gt.p) p=q

    end do

    wk(i)=p

  end do

  do  i=1,m1

    je=m2+i
    js=m1+1-i
    a(i,1:je) = a(i, 1+js: js +je)
    je=je+1
    a(i, je:m ) = 0.0D+00
  end do

    nm1=n-1
    l=m1
    is(1:n) = 0

    do   k=1,nm1
      kp1=k+1
      if(l.lt.n) l=l+1
      i=k
      p=abs(a(k,1))/wk(k)

      do   j=kp1,l
        q=abs(a(j,1))/wk(j)

        if(q.gt.p) then
          p=q
          i=j
        end if
      end do

      if( tols < abs(p) ) then

          if( i /= k) then

              do   j=1,m
                x=a(k,j)
                a(k,j)=a(i,j)
                a(i,j)=x
              end do

              x=b(k)
              b(k)=b(i)
              b(i)=x

          end if

          do  i=kp1,l
            x=a(i,1)/a(k,1)

            do   j=2,m
              a(i,j-1)=a(i,j)-x*a(k,j)
            end do

            a(i,m)=0.0D+00
            b(i)=b(i)-x*b(k)
          end do

        else

          is(k)=1

          do   i=kp1,l
            a(i, 1:m-1) = a(i,2:m)
            a(i,m)=0.0D+00
          end do

      end if

    end do

  if(dabs(a(n,1)).le.tols) is(n)=1
  np1=n+1
  l=1

  do   ii=1,n
    i=np1-ii

    if(is(i) == 1) then
        b(i)=1.0D+00
      else
        im1=i-1
        x=b(i)

        if(i < n) then
            do k=2,l
              x=x-a(i,k)*b(im1+k)
            end do
        end if

        b(i)=x/a(i,1)

    end if

    if(l.lt.m) l=l+1

  end do

  return
end
subroutine bandwidth ( element_num, node_num, element_node, indx, n_lband )

!*****************************************************************************80
!
!! BANDWIDTH computes the half band width of the system matrix.
!
!  Discussion:
!
!    The matrix is assumed to be symmetric and banded, but with
!    an unknown bandwidth.  This routine determines the half bandwidth.
!    The half bandwidth reports, for any row I of the matrix,
!    the maximum distance from the diagonal entry A(I,I) to any nonzero
!    entry A(I,J).
!
!    Knowing this value, the matrix can be stored and factored in
!    a compact way.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 August 2005
!
!  Parameters:
!
!    Input, integer ELEMENT_NUM, the number of elements.
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, integer ELEMENT_NODE(6,ELEMENT_NUM), the nodes of each element.
!
!    Input, integer INDX(5,NODE_NUM), at each node, lists the
!    index of the finite element coefficients for U, V, B1, B2, and P.
!
  implicit none

  integer element_num
  integer node_num

  integer element
  integer, dimension ( 6, element_num ) :: element_node
  integer i
  integer, dimension(5,node_num) :: indx
  integer ip
  integer ipp
  integer iq
  integer iqq
  integer iuk
  integer iukk
  integer j
  integer n_lband

  n_lband = 0

  do element = 1, element_num

    do iq = 1, 6
      ip = element_node(iq,element)

      do iuk = 1, 5

        i = indx(iuk,ip)

        if ( 0 < i ) then

          do iqq = 1, 6
            ipp = element_node(iqq,element)

            do iukk = 1, 5

              j = indx(iukk,ipp)
              if ( i <= j ) then
                n_lband = max ( n_lband, j-i )
              end if
            end do
          end do
        end if
      end do
    end do
  end do

  return
end
subroutine basis_t3 ( t, i, p, qi, dqidx, dqidy )

!*****************************************************************************80
!
!! BASIS_T3 evaluates a linear basis function.
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
!    17 June 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T(2,3), the coordinates of the nodes.
!
!    Input, integer I, the index of the desired basis function.
!    I should be between 1 and 3.
!
!    Input, real ( kind = 8 ) P(2), the point at which the basis
!    function is to be evaluated.
!
!    Output, real ( kind = 8 ) QI, DQIDX, DQIDY, the values of the basis
!    function and its X and Y derivatives.
!
  implicit none

  real ( kind = 8 ) area
  real ( kind = 8 ) dqidx
  real ( kind = 8 ) dqidy
  integer i
  integer i4_wrap
  integer ip1
  integer ip2
  real ( kind = 8 ) p(2)
  real ( kind = 8 ) qi
  real ( kind = 8 ) t(2,3)

  area = t(1,1) * ( t(2,2) - t(2,3) ) &
       + t(1,2) * ( t(2,3) - t(2,1) ) &
       + t(1,3) * ( t(2,1) - t(2,2) )

! if ( area <= 0.0D+00 ) then
!   write ( *, '(a)' ) ' '
!   write ( *, '(a)' ) 'BASIS_T3 - Fatal error!'
!   write ( *, '(a)' ) '  Element has non-positive area.'
!   write ( *, '(a,g14.6)' ) '  Area = ', area
!   stop
! end if

  if ( i < 1 .or. 3 < i ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BASIS_T3 - Fatal error!'
    write ( *, '(a)' ) '  Basis index I is not between 1 and 3.'
    write ( *, '(a,i6)' ) '  I = ', i
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
subroutine basis_t6 ( t, i, p, bi, dbidx, dbidy )

!*****************************************************************************80
!
!! BASIS_T6 evaluates a quadratic basis function.
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
!    17 June 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T(2,6), the coordinates of the nodes.
!
!    Input, integer I, the index of the desired basis function.
!    T should be between 1 and 6.
!
!    Input, real ( kind = 8 ) P(2), the coordinates of a point at which the basis
!    function is to be evaluated.
!
!    Output, real ( kind = 8 ) BI, DBIDX, DBIDY, the values of the basis function
!    and its X and Y derivatives.
!
  implicit none

  real ( kind = 8 ) bi
  real ( kind = 8 ) dbidx
  real ( kind = 8 ) dbidy
  real ( kind = 8 ) gf
  real ( kind = 8 ) gn
  real ( kind = 8 ) hf
  real ( kind = 8 ) hn
  integer i
  integer i4_wrap
  integer j1
  integer j2
  integer k1
  integer k2
  real ( kind = 8 ) p(2)
  real ( kind = 8 ) t(2,6)

  if ( i < 1 .or. 6 < i ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BASIS_T6 - Fatal error!'
    write ( *, '(a)' ) '  Basis index I is not between 1 and 6.'
    write ( *, '(a,i6)' ) '  I = ', i
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
subroutine element_area ( element_num, element_node, node_num, xc, yc, area )

!*****************************************************************************80
!
!! ELEMENT_AREA sets the area of each element.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 August 2005
!
!  Parameters:
!
!    Input, integer ELEMENT_NUM, the number of elements.
!
!    Input, integer ELEMENT_NODE(6,ELEMENT_NUM), the nodes of each element.
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, real ( kind = 8 ) XC(NODE_NUM), YC(NODE_NUM), the coordinates
!    of nodes.
!
!    Output, real ( kind = 8 ) AREA(ELEMENT_NUM), the area of each element.
!
  implicit none

  integer element_num
  integer node_num

  real ( kind = 8 ), dimension(element_num) :: area
  integer element
  integer, dimension ( 6, element_num ) :: element_node
  real ( kind = 8 ), dimension(node_num) :: xc
  real ( kind = 8 ), dimension(node_num) :: yc

  do element = 1, element_num

    area(element) = 0.5D+00 * abs ( &
        xc(element_node(1,element)) &
      * ( yc(element_node(2,element)) - yc(element_node(3,element)) ) &
      + xc(element_node(2,element)) &
      * ( yc(element_node(3,element)) - yc(element_node(1,element)) ) &
      + xc(element_node(3,element)) &
      * ( yc(element_node(1,element)) - yc(element_node(2,element)) ) )

  end do

  return
end
subroutine element_set ( nx, ny, element_num, element_node )

!*****************************************************************************80
!
!! ELEMENT_SET sets up the elements.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 October 2006
!
!  Parameters:
!
!    Input, integer NX, NY, the number of rows and columns of elements.
!
!    Input, integer ELEMENT_NUM, the number of elements.
!
!    Output, integer ELEMENT_NODE(6,ELEMENT_NUM), the nodes of each element.
!
  implicit none

  integer :: element_num

  integer :: col
  integer :: col_num
  logical :: col_odd
  integer :: element
  integer, dimension ( 6, element_num ) :: element_node
  integer :: node_num
  integer :: nx
  integer :: ny
  integer :: row
  integer :: row_num
  logical :: row_odd
!
!  We start out with no elements.
!
!  We "walk through" the nodes, column by column.
!
!  At node N, if we're in an odd column and an odd row (except for the
!  LAST row or column), we build the two elements which extend
!  from node N to the right by two nodes, and up by two nodes:
!
!    N+2-N+NR+2-N+2*NR+2  odd row
!    |            /|
!    |           / |
!    |          /  |
!    |         /   |
!    |        /    |
!    |       /     |
!    N+1 N+NR+1 N+2*NR+1  even row
!    |     /       |
!    |    /        |
!    |   /         |
!    |  /          |
!    | /           |
!    |/            |
!    N----N+NR--N+2*NR    odd row
!
!   odd   even    odd
!   col   col     col
!
!
!  The local numbering of the nodes is
!
!   2--5--3       2
!   |    /       /|
!   |   /       / |
!   4  6       4  5
!   | /       /   |
!   |/       /    |
!   1       1--6--3
!
  col_num = 2 * nx - 1
  row_num = 2 * ny - 1

  element = 0
  node_num = 0

  do col = 1, col_num

    col_odd = ( mod ( col, 2 ) == 1 )

    do row = 1, row_num

      row_odd = ( mod ( row, 2 ) == 1 )

      node_num = node_num + 1

      if ( col_odd .and. row_odd ) then

        if ( col /= col_num .and. row /= row_num ) then

          element = element + 1
          element_node(1,element) = node_num
          element_node(2,element) = node_num               + 2
          element_node(3,element) = node_num + 2 * row_num + 2
          element_node(4,element) = node_num               + 1
          element_node(5,element) = node_num +     row_num + 2
          element_node(6,element) = node_num +     row_num + 1

          element = element + 1
          element_node(1,element) = node_num
          element_node(2,element) = node_num + 2 * row_num + 2
          element_node(3,element) = node_num + 2 * row_num
          element_node(4,element) = node_num +     row_num + 1
          element_node(5,element) = node_num + 2 * row_num + 1
          element_node(6,element) = node_num +     row_num

        end if

      end if

    end do

  end do

  return
end
subroutine element_write ( element_file_name, element_num, element_node )

!*****************************************************************************80
!
!! ELEMENT_WRITE writes the element data to a file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 October 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) ELEMENT_FILE_NAME, the name of the
!    element file.
!
!    Input, integer ELEMENT_NUM, the number of elements.
!
!    Input, integer ELEMENT_NODE(6,ELEMENT_NUM), the nodes of each element.
!
  implicit none

  integer :: element_num

  integer :: element
  character ( len = * ) element_file_name
  integer, dimension ( 6, element_num ) :: element_node
  integer :: iunit

  call get_unit ( iunit )

  open ( unit = iunit, file = element_file_name, status = 'replace' )

  do element = 1, element_num

    write ( iunit, '(6(2x,i8))' ) element_node(1:6,element)

  end do

  close ( unit = iunit )

  return
end
subroutine eval ( g, indx, element, element_num, node_num, &
  element_node, x, y, xc, yc, ylength, y_inflow, y_outflow, vel, velx, &
  vely, mag, magx, magy )

!*****************************************************************************80
!
!! EVAL evaluates the velocity and magnetic fields at a point.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 August 2005
!
!  Parameters:
!
!    Input, real ( kind = 8 ) G(*), the finite element coefficients.
!
!    Input, integer INDX(5,NODE_NUM), at each node, lists the
!    index of the finite element coefficients for U, V, B1, B2, and P.
!
!    Input, integer ELEMENT, the index of the element.
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, integer ELEMENT_NODE(6,ELEMENT_NUM), the nodes of each element.
!
!    Input, real ( kind = 8 ) X, Y, the coordinates of the point
!    at which the fields are to be evaluated.
!
!    Input, real ( kind = 8 ) XC(NODE_NUM), YC(NODE_NUM), the coordinates
!    of nodes.
!
!    Input, real ( kind = 8 ) YLENGTH, the maximum Y coordinate at which
!    there is outflow .
!
!    Input, real ( kind = 8 ) Y_INFLOW, the maximum Y coordinate at which
!    there is inflow.
!
!    Input, real ( kind = 8 ) Y_OUTFLOW, the minimum Y coordinate at which
!    there is outflow.
!
!    Output, real ( kind = 8 ) VEL(2), VELX(2), VELY(2), the value of the
!    horizonal and vertical components of velocity, and their X and Y
!    derivatives.
!
!    Output, real ( kind = 8 ) MAG(2), MAGX(2), MAGY(2), the value of the
!    horizonal and vertical components of the magnetic field, and their X
!    and Y derivatives.
!
  implicit none

  integer :: element_num
  integer :: node_num

  real ( kind = 8 ) b
  real ( kind = 8 ) dbdx
  real ( kind = 8 ) dbdy
  integer :: element
  integer, dimension ( 6, element_num ) :: element_node
  real ( kind = 8 ), dimension(*) :: g
  integer, dimension(5,node_num) :: indx
  integer :: ip
  integer :: iq
  integer :: iuk
  integer :: iun
  real ( kind = 8 ), dimension(2) :: mag
  real ( kind = 8 ), dimension(2) :: magx
  real ( kind = 8 ), dimension(2) :: magy
  real ( kind = 8 ) p(2)
  real ( kind = 8 ) t6(2,6)
  real ( kind = 8 ) :: ubc
  real ( kind = 8 ) :: ubdry
  real ( kind = 8 ), dimension(2) :: vel
  real ( kind = 8 ), dimension(2) :: velx
  real ( kind = 8 ), dimension(2) :: vely
  real ( kind = 8 ) :: x
  real ( kind = 8 ), dimension(node_num) :: xc
  real ( kind = 8 ) :: y
  real ( kind = 8 ), dimension(node_num) :: yc
  real ( kind = 8 ) :: ylength
  real ( kind = 8 ) :: y_inflow
  real ( kind = 8 ) :: y_outflow
  real ( kind = 8 ) :: yy

  vel(1:2)  = 0.0D+00
  velx(1:2) = 0.0D+00
  vely(1:2) = 0.0D+00

  mag(1:2)  = 0.0D+00
  magx(1:2) = 0.0D+00
  magy(1:2) = 0.0D+00

  t6(1,1:6) = xc(element_node(1:6,element))
  t6(2,1:6) = yc(element_node(1:6,element))
  p(1:2) = (/ x, y /)

  do iq = 1, 6

    call basis_t6 ( t6, iq, p, b, dbdx, dbdy )

    ip = element_node(iq,element)
    yy = yc(ip)

    do iuk = 1, 2

      iun = indx(iuk,ip)
!
!  If 0 < IUN, then the value at this node is stored in a finite element
!  coefficient.
!
!  If 0 == IUN, then the value at this node is zero.
!
!  If IUN < 0, then the value at this node is a boundary value,
!  obtained by calling the UBDRY function.
!
      if ( 0 < iun ) then

        vel(iuk)  = vel(iuk)  +  b   * g(iun)
        velx(iuk) = velx(iuk) + dbdx * g(iun)
        vely(iuk) = vely(iuk) + dbdy * g(iun)

      else if ( iun == 0 ) then

      else if ( iun < 0 ) then

        ubc = ubdry ( iun, yy, ylength, y_inflow, y_outflow )

        vel(iuk)  = vel(iuk)  +  b   * ubc
        velx(iuk) = velx(iuk) + dbdx * ubc
        vely(iuk) = vely(iuk) + dbdy * ubc

      end if

    end do
!
!  All values of the magnetic field are stored in the finite
!  element coefficient vector, or are zero.  (There are no
!  boundary values specified.)
!
    do iuk = 1, 2

      iun = indx(iuk+2,ip)

      if ( 0 < iun ) then

        mag(iuk)  = mag(iuk)  +  b   * g(iun)
        magx(iuk) = magx(iuk) + dbdx * g(iun)
        magy(iuk) = magy(iuk) + dbdy * g(iun)

      else if ( 0 == iun ) then

      end if

    end do

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
function i4_modp ( i, j )

!*****************************************************************************80
!
!! I4_MODP returns the nonnegative remainder of I4 division.
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
!    Input, integer I, the number to be divided.
!
!    Input, integer J, the number that divides I.
!
!    Output, integer I4_MODP, the nonnegative remainder when I is
!    divided by J.
!
  implicit none

  integer i
  integer i4_modp
  integer j

  if ( j == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_MODP - Fatal error!'
    write ( *, '(a,i6)' ) '  I4_MODP ( I, J ) called with J = ', j
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
!! I4_WRAP forces an I4 to lie between given limits by wrapping.
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
!    Input, integer IVAL, an integer value.
!
!    Input, integer ILO, IHI, the desired bounds for the integer value.
!
!    Output, integer I4_WRAP, a "wrapped" version of IVAL.
!
  implicit none

  integer i4_modp
  integer i4_wrap
  integer ihi
  integer ilo
  integer ival
  integer jhi
  integer jlo
  integer wide

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
subroutine indx_set ( node_num, nx, ny, y_in_node, y_out_node, &
  indx, unknown_num )

!*****************************************************************************80
!
!! INDX_SET sets the index array.
!
!  Discussion:
!
!    The INDX array records, for each node, the indices of
!    the unknowns associated with:
!      1) U, the horizontal velocity,
!      2) V, the vertical velocity,
!      3) B1, the horizontal magnetic field component,
!      4) B2, the vertical magnetic field component,
!      5) P, pressure.
!
!    If INDX(variable,node) is positive, then it is reporting the
!    corresponding index in the finite element coefficient array
!    for the variable associated with the node.
!
!    If INDEX(node,variable) is zero, then it is indicating that the
!    corresponding variable at the node is constrained to be zero.
!
!    If INDEX(node,variable) is negative, then it is indicating that the
!    corresponding variable at the node is not unknown, but is determined
!    by a boundary condition.
!
!    The actual boundary conditions for this problem are:
!
!    U: zero on the top and bottom walls.
!       on the left wall, parabolic inflow condition between nodes 1 and
!       Y_INFLOW, otherwise zero;
!       on the right wall, a finite element variable between nodes Y_OUTFLOW
!       and NODE_NUM, otherwise zero;
!
!    V: zero on all walls.
!
!    B1: zero on left and right walls (no normal component of magnetic field)
!        and at corners.
!
!    B2: zero on top and bottom walls (no normal component of magnetic field)
!        and at corners.
!
!    P: finite element variable at nodes that are in odd row and column.
!       otherwise, set INDX(I,J) = 0.  P is not actually zero at those
!       nodes, but the finite element coefficient is not needed there.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 October 2006
!
!  Parameters:
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, integer NX, NY, the number of rows and columns of elements.
!
!    Input, integer Y_IN_NODE, nodes 1 through Y_IN_NODE are on the
!    inflow boundary.
!
!    Input, integer Y_OUT_NODE, nodes Y_OUT_NODE through NODE_NUM
!    are on the outflow boundary.
!
!    Output, integer INDX(5,NODE_NUM), at each node, lists the
!    index of the finite element coefficients for U, V, B1, B2, and P.
!
!    Output, integer UNKNOWN_NUM, the number of unknowns.
!
  implicit none

  integer :: node_num

  integer :: col
  integer :: col_num
  logical :: col_odd
  integer, dimension(5,node_num) :: indx
  integer :: unknown_num
  integer :: nx
  integer :: ny
  integer :: row
  integer :: row_num
  logical :: row_odd
  integer :: y_in_node
  integer :: y_out_node

  col_num = 2 * nx - 1
  row_num = 2 * ny - 1

  unknown_num = 0
  node_num = 0

  do col = 1, col_num

    col_odd = ( mod ( col, 2 ) == 1 )

    do row = 1, row_num

      row_odd = ( mod ( row, 2 ) == 1 )

      node_num = node_num + 1
!
!  Is this node one of the four corner nodes?
!
      if ( ( row == 1 .or. row == row_num ) .and. &
           ( col == 1 .or. col == col_num ) ) then

        indx(1,node_num) = 0
        indx(2,node_num) = 0
        indx(3,node_num) = 0
        indx(4,node_num) = 0
        unknown_num = unknown_num + 1
        indx(5,node_num) = unknown_num
!
!  Is this node in the first or last row (but not the corners)?
!
      else if ( row == 1 .or. row == row_num ) then

        indx(1,node_num) = 0
        indx(2,node_num) = 0
        unknown_num = unknown_num + 1
        indx(3,node_num) = unknown_num
        indx(4,node_num) = 0

        if ( col_odd ) then
          unknown_num = unknown_num + 1
          indx(5,node_num) = unknown_num
        else
          indx(5,node_num) = 0
        end if
!
!  Is this node in the first column (but not the corners)?
!
      else if ( col == 1 ) then

        if ( node_num <= y_in_node ) then
          indx(1,node_num) = -1
        else
          indx(1,node_num) = 0
        end if

        indx(2,node_num) = 0
        indx(3,node_num) = 0
        unknown_num = unknown_num + 1
        indx(4,node_num) = unknown_num

        if ( row_odd ) then
          unknown_num = unknown_num + 1
          indx(5,node_num) = unknown_num
        else
          indx(5,node_num) = 0
        end if
!
!  Is this node in the last column (but not the corners)?
!
      else if ( col == col_num ) then

        if ( y_out_node <= node_num ) then
          unknown_num = unknown_num + 1
          indx(1,node_num) = unknown_num
        else
          indx(1,node_num) = 0
        end if

        indx(2,node_num) = 0
        indx(3,node_num) = 0
        unknown_num = unknown_num + 1
        indx(4,node_num) = unknown_num

        if ( row_odd ) then
          unknown_num = unknown_num + 1
          indx(5,node_num) = unknown_num
        else
          indx(5,node_num) = 0
        end if
!
!  Is this node in an "interior" row and column?
!
      else

        unknown_num = unknown_num + 1
        indx(1,node_num) = unknown_num
        unknown_num = unknown_num + 1
        indx(2,node_num) = unknown_num
        unknown_num = unknown_num + 1
        indx(3,node_num) = unknown_num
        unknown_num = unknown_num + 1
        indx(4,node_num) = unknown_num

        if ( row_odd .and. col_odd ) then
          unknown_num = unknown_num + 1
          indx(5,node_num) = unknown_num
        else
          indx(5,node_num) = 0
        end if

      end if

    end do

  end do

  return
end
subroutine jacobian ( area, cp, dt, f_old, f_prev, indx, element_max, &
  node_num, element_num, n_lband, element_node, quad_num, unknown_num, &
  rey_mag, visc, xc, xq, yc, yq, ylength, y_inflow, y_outflow,  &
  par_sim, a, f_new )

!*****************************************************************************80
!
!! JACOBIAN computes the Jacobian and right hand side of the Newton system.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 August 2005
!
!  Reference:
!
!    Max Gunzburger and Catalin Trenchea,
!    Analysis and Discretization of an Optimal Control Problem
!    for the Time-Periodic MHD Equations,
!    Journal of Mathematical Analysis and Applications,
!    Volume 308, Number 2, 2005, pages 440-466.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) AREA(ELEMENT_NUM), the area of each
!    element.
!
!    Input, real ( kind = 8 ) CP, the coupling parameter.
!
!    Input/output, real ( kind = 8 ) DT, the time step.
!
!    Input, real ( kind = 8 ) F_OLD(), ?
!
!    Input, real ( kind = 8 ) F_PREV(), ?
!
!    Input, integer INDX(5,NODE_NUM), at each node, lists the
!    index of the finite element coefficients for U, V, B1, B2, and P.
!
!    Input, integer ELEMENT_MAX, the maximum number of elements.
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, integer ELEMENT_NUM, the number of elements.
!
!    Input, integer N_LBAND, the half-bandwidth of the matrix.
!
!    Input, integer ELEMENT_NODE(6,ELEMENT_NUM), the nodes of each element.
!
!    Input, integer QUAD_NUM, the order of the quadrature rule.
!
!    Input, integer UNKNOWN_NUM, the number of unknowns.
!
!    Input, real ( kind = 8 ) REY_MAG, the magnetic Reynolds number.
!
!    Input, real ( kind = 8 ) VISC, the fluid viscosity.
!
!    Input, real ( kind = 8 ) XC(NODE_NUM), the X coordinates of nodes.
!
!    Input, real ( kind = 8 ) XQ(ELEMENT_MAX,3), the X coordinates of
!    the quadrature points.
!
!    Input, real ( kind = 8 ) YC(NODE_NUM), the Y coordinates of nodes.
!
!    Output, real ( kind = 8 ) YQ(ELEMENT_MAX,3), the Y coordinates of
!    the quadrature points.
!
!    Input, real ( kind = 8 ) YLENGTH, the maximum Y coordinate at which
!    there is outflow.
!
!    Input, real ( kind = 8 ) Y_INFLOW, the maximum Y coordinate at which
!    there is inflow.
!
!    Input, real ( kind = 8 ) Y_OUTFLOW, the minimum Y coordinate at which
!    there is outflow.
!
!    Input, real ( kind = 8 ) PAR_SIM, is 0 if a simple iteration step
!    is being taken, and 1 if a Newton iteration step is being taken.
!
!    Output, real ( kind = 8 ) A(UNKNOWN_NUM,2*N_LBAND+1), the banded Jacobian
!    matrix needed by the Newton solver.
!
!    Output, real ( kind = 8 ) F_NEW(UNKNOWN_NUM), the current value of
!    the residual.
!
  implicit none

  integer :: element_max
  integer :: element_num
  integer :: node_num
  integer :: unknown_num
  integer :: n_lband

  real ( kind = 8 ), dimension(unknown_num,2*n_lband+1) :: a
  real ( kind = 8 ) :: aij
  real ( kind = 8 ) :: ar
  real ( kind = 8 ), dimension(element_num) :: area
  real ( kind = 8 ) :: bi
  real ( kind = 8 ) :: bj
  real ( kind = 8 ) :: cp
  real ( kind = 8 ) :: dbidx
  real ( kind = 8 ) :: dbidy
  real ( kind = 8 ) :: dbjdx
  real ( kind = 8 ) :: dbjdy
  real ( kind = 8 ) :: dqidx
  real ( kind = 8 ) :: dqidy
  real ( kind = 8 ) :: dqjdx
  real ( kind = 8 ) :: dqjdy
  real ( kind = 8 ) :: dt
  integer :: element
  integer, dimension ( 6, element_num ) :: element_node
  real ( kind = 8 ), dimension(unknown_num) :: f_new
  real ( kind = 8 ), dimension(unknown_num) :: f_old
  real ( kind = 8 ), dimension(unknown_num) :: f_prev
  integer :: i
  integer, dimension(5,node_num) :: indx
  integer :: ip
  integer :: ipp
  integer :: iq
  integer :: iqq
  integer :: iuk
  integer :: iukk
  integer :: iuse
  integer :: j
  real ( kind = 8 ), dimension(2) :: mag
  real ( kind = 8 ), dimension(2) :: mag_prev
  real ( kind = 8 ), dimension(2) :: magx
  real ( kind = 8 ), dimension(2) :: magy
  real ( kind = 8 ) :: p(2)
  real ( kind = 8 ) :: par_sim
  real ( kind = 8 ) :: qi
  real ( kind = 8 ) :: qj
  integer :: quad
  integer :: quad_num
  real ( kind = 8 ) :: rey_mag
  real ( kind = 8 ) :: rhs
  real ( kind = 8 ) :: t3(2,3)
  real ( kind = 8 ) :: t6(2,6)
  real ( kind = 8 ) :: ubdry
  real ( kind = 8 ) :: value
  real ( kind = 8 ), dimension(2) :: vel
  real ( kind = 8 ), dimension(2) :: vel_prev
  real ( kind = 8 ), dimension(2) :: velx
  real ( kind = 8 ), dimension(2) :: vely
  real ( kind = 8 ) :: visc
  real ( kind = 8 ) :: x
  real ( kind = 8 ), dimension(node_num) :: xc
  real ( kind = 8 ), dimension(element_max,3) :: xq
  real ( kind = 8 ) :: y
  real ( kind = 8 ) :: y_inflow
  real ( kind = 8 ) :: y_outflow
  real ( kind = 8 ), dimension(node_num) :: yc
  real ( kind = 8 ) :: ylength
  real ( kind = 8 ), dimension(element_max,3) :: yq
  real ( kind = 8 ) :: yy

  a(1:unknown_num,1:2*n_lband+1) = 0.0D+00

  do element = 1, element_num

    ar = area(element) / 3.0D+00

    t3(1,1:3) = xc(element_node(1:3,element))
    t3(2,1:3) = yc(element_node(1:3,element))

    t6(1,1:6) = xc(element_node(1:6,element))
    t6(2,1:6) = yc(element_node(1:6,element))
!
!  Loop over quadrature points.
!
    do quad = 1, quad_num

      x = xq(element,quad)
      y = yq(element,quad)
      p(1:2) = (/ x, y /)
!
!  Evaluate the velocities and magnetic field at the previous time step.
!
      call eval ( f_prev, indx, element, element_num, node_num,  &
        element_node, x, y, xc, yc, ylength, y_inflow, y_outflow, &
        vel_prev, velx, vely, mag_prev, magx, magy )
!
!  Evaluate velocities and magnetic field at the previous Newton iterate.
!
      call eval ( f_old, indx, element, element_num, node_num,  &
        element_node, x, y, xc, yc, ylength, y_inflow, y_outflow, &
        vel, velx, vely, mag, magx, magy )

      do iq = 1, 6

        ip = element_node(iq,element)

        call basis_t6 ( t6, iq, p, bi, dbidx, dbidy )

        if ( iq <= 3 ) then
          call basis_t3 ( t3, iq, p, qi, dqidx, dqidy )
        end if

        do iuk = 1, 5

          i = indx(iuk,ip)

          if ( 0 < i .and. i /= unknown_num ) then

            if ( iuk == 1 ) then

              f_new(i) = f_new(i) + ar * ( rhs ( iuk, x, y ) &
                + bi * par_sim * ( vel(1) * velx(1) + vel(2) * vely(1) ) &
                + bi * vel_prev(1) / dt &
                + bi * par_sim * cp * mag(2) * ( magx(2) - magy(1) ) )

            else if ( iuk == 2 ) then

              f_new(i) = f_new(i) + ar * ( rhs ( iuk, x, y ) &
                + bi * par_sim * ( vel(1) * velx(2) + vel(2) * vely(2) ) &
                + bi * vel_prev(2) / dt   &
                - bi * par_sim * cp * mag(1) * ( magx(2) - magy(1) ) )

            else if ( iuk == 3 ) then

              f_new(i) = f_new(i) + ar * ( rhs ( iuk, x, y ) &
                + bi * mag_prev(1) / dt       &
                - bi * par_sim * ( vely(1) * mag(2) - vely(2) * mag(1) ) &
                - bi * par_sim * ( vel(1) * magy(2) - vel(2) * magy(1) ) )

            else if ( iuk == 4 ) then

              f_new(i) = f_new(i) +  ar * ( rhs ( iuk, x, y ) &
                + bi * mag_prev(2) / dt   &
                + bi * par_sim * ( velx(1) * mag(2) - velx(2) * mag(1) ) &
                + bi * par_sim * ( vel(1) * magx(2) - vel(2) * magx(1) ) )

            else if ( iuk == 5 ) then

            end if

            do iqq = 1, 6

              ipp = element_node(iqq,element)
              yy = yc(ipp)

              call basis_t6 ( t6, iqq, p, bj, dbjdx, dbjdy )

              if ( iqq <= 3 ) then
                call basis_t3 ( t3, iqq, p, qj, dqjdx, dqjdy )
              end if

              do iukk = 1, 5

                j = indx(iukk,ipp)
!
!  0 < j implies there is an unknown at this node,
!  J<0 implies inhomogeneous Dirichlet
!
                if ( j /= 0 ) then

                  aij = 0.0D+00

                  if ( iuk == 1 ) then

                    if ( iukk == 1 )  then

                      aij = bi * bj / dt &
                        + visc * ( dbidx * dbjdx + dbidy * dbjdy )  &
                        + bi * bj * velx(1) * par_sim &
                        + bi * ( dbjdx * vel(1) + dbjdy * vel(2) )

                    else if ( iukk == 2 ) then

                      aij = par_sim * bi * bj * vely(1)

                    else if ( iukk == 3 ) then

                      aij = -cp * mag(2) * dbjdy * bi

                    else if ( iukk == 4 ) then

                      aij = cp * ( bi * bj &
                        * ( magx(2) - magy(1) ) * par_sim   &
                        + mag(2) * dbjdx * bi  )

                    else if ( iukk == 5 ) then

                      aij = -dbidx * qj

                    end if

                  else if ( iuk == 2 ) then

                    if ( iukk == 1 ) then

                      aij = par_sim * bi * bj * velx(2)

                    else if ( iukk == 2 ) then

                      aij = bi * bj / dt &
                        + visc * ( dbidx * dbjdx + dbidy * dbjdy ) &
                        + bi * bj * vely(2) * par_sim    &
                        + bi * ( dbjdy * vel(2) +  dbjdx * vel(1) )

                    else if ( iukk == 3 ) then

                      aij = -cp * bi * bj &
                        * ( magx(2) - magy(1) ) * par_sim   &
                        + cp * bi * dbjdy * mag(1)

                    else if ( iukk == 4 ) then

                      aij = -cp * bi * dbjdx * mag(1)

                    else if ( iukk == 5 ) then

                      aij = -dbidy * qj

                    end if

                  else if ( iuk == 3 ) then

                    if ( iukk == 1 ) then

                      aij = - bi * dbjdy * mag(2) &
                        - bi * bj * magy(2)

                    else if ( iukk == 2 ) then

                      aij = bi * dbjdy * mag(1) &
                        + bi * bj * magy(1)

                    else if ( iukk == 3 ) then

                      aij = bi * bj / dt  &
                        + dbidy * dbjdy / rey_mag   &
                        + bi * bj * vely(2) * par_sim  &
                        + bi * dbjdy * vel(2) * par_sim &
                        + dbidx * dbjdx / rey_mag

                    else if ( iukk == 4 ) then

                      aij = - dbidy * dbjdx / rey_mag &
                        - bi * bj * vely(1) * par_sim  &
                        - bi * dbjdy * vel(1) * par_sim &
                        + dbidx * dbjdy / rey_mag

                    else if ( iukk == 5 ) then

                      aij = 0.0D+00

                    end if

                  else if ( iuk == 4 ) then

                    if ( iukk == 1 ) then

                      aij = bi * bj * magx(2) &
                        + bi * dbjdx * mag(2)

                    else if ( iukk == 2 ) then

                      aij = - bi * dbjdx * mag(1) &
                        - bi * bj * magx(1)

                    else if ( iukk == 3 ) then

                      aij = - dbidx * dbjdy / rey_mag    &
                        - bi * bj * velx(2) * par_sim &
                        - bi * dbjdx * vel(2) * par_sim &
                        + dbidy * dbjdx / rey_mag

                    else if ( iukk == 4 ) then

                      aij = bi * bj / dt &
                        + dbidx * dbjdx / rey_mag &
                        + bi * bj * velx(1) * par_sim   &
                        + bi * dbjdx * vel(1) * par_sim &
                        + dbidy * dbjdy / rey_mag

                    else if ( iukk == 5 ) then

                      aij = 0.0D+00

                    end if

                  else if ( iuk == 5 ) then

                    if ( iukk == 1 ) then

                      aij = -dbjdx * qi

                    else if ( iukk == 2 ) then

                      aij = -dbjdy * qi

                    else

                      aij = 0.0D+00

                    end if

                  end if
!
!  If J represents an unknown, the coefficient is added to matrix entry A(I,J).
!
!  If J represents a known boundary value, the coefficient is combined with
!  the value, and subtracted from the right hand side F(I).
!
                  if ( 0 < j ) then
                    iuse = j - i + n_lband + 1
                    a(i,iuse) = a(i,iuse) + aij * ar
                  else if  ( j < 0 ) then
                    value = ubdry ( j, yy, ylength, y_inflow, y_outflow )
                    f_new(i) = f_new(i) - ar * value  * aij
                  end if

                end if

              end do
            end do
          end if
        end do
      end do
    end do
  end do
!
!  To avoid singularity of the pressure system, the last pressure equation
!  is replaced by an assignment of the last pressure to the value 0.
!
!  This manipulation takes advantage of the fact that the last pressure
!  is the last variable.
!
  f_new(unknown_num) = 0.0D+00
!
!  Replace the coefficients of the continuity equation.
!
  a(unknown_num,1:2*n_lband+1) = 0.0D+00
  a(unknown_num,    n_lband+1) = 1.0D+00

  return
end
subroutine mag_write ( f, indx, mag_file_name, node_num )

!*****************************************************************************80
!
!! MAG_WRITE writes magnetic data to files.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 September 2005
!
!  Parameters:
!
!    Input, real ( kind = 8 ) F(*), the current solution.
!
!    Input, integer INDX(5,NODE_NUM), at each node, lists the
!    index of the finite element coefficients for U, V, B1, B2, and P.
!
!    Input, character ( len = * ) MAG_FILE_NAME, the name of the
!    output file containing the magnetic field values.
!
!    Input, integer NODE_NUM, the number of nodes.
!
  implicit none

  integer :: node_num

  real ( kind = 8 ) :: b1_value
  real ( kind = 8 ) :: b2_value
  real ( kind = 8 ), dimension(*) :: f
  integer, dimension(5,node_num) :: indx
  integer :: ip
  integer :: iuk_b1
  integer :: iuk_b2
  integer iunit
  character ( len = * ) mag_file_name

  call get_unit ( iunit )

  open ( unit = iunit, file = mag_file_name, status = 'replace' )

  do ip = 1, node_num

    iuk_b1 = indx(3,ip)
    iuk_b2 = indx(4,ip)

    if ( 0 < iuk_b1 ) then
      b1_value = f(iuk_b1)
    else
      b1_value = 0.0D+00
    end if

    if ( 0 < iuk_b2 ) then
      b2_value = f(iuk_b2)
    else
      b2_value = 0.0D+00
    end if

    write ( iunit, '(2x,g14.6,2x,g14.6)' ) b1_value, b2_value

  end do

  close( unit = iunit )

  return
end
subroutine newton ( area, cp, dt, f_prev, indx, element_max, &
  max_newton, element_num,  n_lband, node_num, element_node, &
  quad_num, n_simple, unknown_num, rey_mag, time_cur, &
  tol_newton, visc, xc, xq, yc, yq, ylength, y_inflow, y_outflow,  &
  a, f_new )

!*****************************************************************************80
!
!! NEWTON updates the magneto-hydrodynamic state using Newton's method.
!
!  Discussion:
!
!    Navier Stokes equations coupled with magnetic field
!    using Taylor-Hood for velocity/pressure and quadratics on triangles
!    for magnetic field B
!
!    The f_prev  array contains the solution at previous time step.
!
!    The f_new array contains the right hand side initially and then the
!    current iterate.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 August 2005
!
!  Reference:
!
!    Max Gunzburger and Catalin Trenchea,
!    Analysis and Discretization of an Optimal Control Problem
!    for the Time-Periodic MHD Equations,
!    Journal of Mathematical Analysis and Applications,
!    Volume 308, Number 2, 2005, pages 440-466.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) AREA(ELEMENT_NUM), the area of each
!    element.
!
!    Input, real ( kind = 8 ) CP, the coupling parameter.
!
!    Input/output, real ( kind = 8 ) DT, the time step.
!
!    Input, integer INDX(5,NODE_NUM), at each node, lists the
!    index of the finite element coefficients for U, V, B1, B2, and P.
!
!    Input, integer ELEMENT_MAX, the maximum number of elements.
!
!    Input, integer ELEMENT_NUM, the number of elements.
!
!    Input, integer ELEMENT_NODE(6,ELEMENT_NUM), the nodes of each element.
!
!    Input, real ( kind = 8 ) XC(NODE_NUM), YC(NODE_NUM), the coordinates
!    of nodes.
!
!    Input, real ( kind = 8 ) YLENGTH, the maximum Y coordinate at which
!    there is outflow.
!
!    Input, real ( kind = 8 ) Y_INFLOW, the maximum Y coordinate at which
!    there is inflow.
!
!    Input, real ( kind = 8 ) Y_OUTFLOW, the minimum Y coordinate at which
!    there is outflow.
!
  implicit none

  integer :: element_max
  integer :: element_num
  integer :: node_num
  integer :: unknown_num

  real ( kind = 8 ), dimension(unknown_num,*) :: a
  real ( kind = 8 ), dimension(element_num) :: area
  real ( kind = 8 ) :: cp
  real ( kind = 8 ) :: diff
  real ( kind = 8 ) :: diff_max
  integer :: diff_max_index
  real ( kind = 8 ) :: dt
  integer, dimension ( 6, element_num ) :: element_node
  real ( kind = 8 ), dimension(unknown_num) :: f_new
  real ( kind = 8 ), dimension(unknown_num) :: f_old
  real ( kind = 8 ), dimension(unknown_num) :: f_prev
  integer i
  integer, dimension(5,node_num) :: indx
  integer :: iter
  integer :: max_newton
  integer :: n_lband
  integer :: n_simple
  real ( kind = 8 ) :: normf
  real ( kind = 8 ) :: par_sim
  integer :: quad_num
  real ( kind = 8 ) :: rey_mag
  real ( kind = 8 ) :: time_cur
  real ( kind = 8 ) :: tol_newton
  real ( kind = 8 ) :: visc
  real ( kind = 8 ), dimension(node_num) :: xc
  real ( kind = 8 ), dimension(element_max,3) :: xq
  real ( kind = 8 ) :: y_inflow
  real ( kind = 8 ) :: y_outflow
  real ( kind = 8 ), dimension(node_num) :: yc
  real ( kind = 8 ) :: ylength
  real ( kind = 8 ), dimension(element_max,3) :: yq
!
!  Initialize F_OLD, which stores the previous Newton iterate.
!
  f_old(1:unknown_num) = f_prev(1:unknown_num)
!
!  Try to take a time step of size DT.
!  (If the step fails, we will try again with a smaller DT.)
!
  do

    f_new(1:unknown_num) = 0.0D+00
!
!  For the given time step, seek solution values that satisfy the
!  nonlinear relation, using a hybrid Newton method.
!  We first take N_SIMPLE simple iterations, followed by standard Newton.
!
    do iter = 1, max_newton

      if ( iter <= n_simple ) then
        par_sim = 0.0D+00
      else if ( n_simple < iter ) then
        par_sim = 1.0D+00
      end if

      call jacobian ( area, cp, dt, f_old, f_prev, indx, element_max, &
        node_num, element_num,  n_lband, element_node, &
        quad_num, unknown_num, rey_mag, &
        visc, xc, xq, yc, yq, ylength, y_inflow, y_outflow,  &
        par_sim, a, f_new )

      call banded ( a, f_new, unknown_num, n_lband, n_lband, unknown_num )
!
!  Normalize the pressure to have zero mean.
!
      call p_normalize ( area, f_new, indx, element_num, node_num, &
        element_node )
!
!  Check for convergence by looking at normalized difference in
!  successive iterates.
!
      normf = sqrt ( sum ( f_new(1:unknown_num)**2 ) )

      diff = sqrt ( &
        sum ( ( f_old(1:unknown_num) - f_new(1:unknown_num) )**2 ) )

      diff = diff / normf

      diff_max = -1.0D+00

      do i = 1, unknown_num
        if ( diff_max < abs ( f_new(i) - f_old(i) ) ) then
          diff_max = abs ( f_new(i) - f_old(i) )
          diff_max_index = i
        end if
      end do

      if ( iter == 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '    Iter  L2(Diff)  Max index  Max val'
        write ( *, '(a)' ) ' '
      end if

      write ( *, '(2x,i6,2x,g14.6,2x,i6,2x,g14.6)' ) &
        iter, diff, diff_max_index, diff_max

      if ( diff < tol_newton ) then

        write (*, '(a)' ) '  The Newton iteration has converged.'
        return

      end if

      f_old(1:unknown_num) = f_new(1:unknown_num)
      f_new(1:unknown_num) = 0.0D+00

    end do

    time_cur = time_cur - dt
    dt = dt / 2.0D+00
    write (*, * ) '  Newton did not converge; decreasing delta t to  ', dt

    if ( dt < 0.00011D+00 ) then
      write ( *, * ) '  Reduced delta T too much so stop '
      stop
    end if
!
!  Set initial Newton guess to solution at previous timestep.
!
    time_cur = time_cur + dt
    f_old(1:unknown_num) = f_prev(1:unknown_num)

  end do

  return
end
subroutine p_normalize ( area, f, indx, element_num, node_num, &
  element_node )

!*****************************************************************************80
!
!! P_NORMALIZE normalizes the pressure.
!
!  Discussion:
!
!    The singular pressure system is regularized by requiring that
!    a particular pressure have the value 0.  Once the pressures
!    have been computed, they are renormalized by computing the mean
!    and subtracting that from all the pressures.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 August 2005
!
!  Parameters:
!
!    Input, real ( kind = 8 ) AREA(ELEMENT_NUM), the area of each
!    element.
!
!    Input/output, real ( kind = 8 ) F(*), the finite element coefficients.
!    On output, the coefficients associated with pressure have been
!    adjusted so that the mean pressure is 0.
!
!    Input, integer INDX(5,NODE_NUM), at each node, lists the
!    index of the finite element coefficients for U, V, B1, B2, and P.
!
!    Input, integer ELEMENT_NUM, the number of elements.
!
!    Input, integer NODE_NODE, the number of nodes.
!
!    Input, integer ELEMENT_NODE(6,ELEMENT_NUM), the nodes of each element.
!
  implicit none

  integer :: element_num
  integer :: node_num

  real ( kind = 8 ), dimension(element_num) :: area
  real ( kind = 8 ) area_region
  integer :: element
  integer, dimension ( 6, element_num ) :: element_node
  real ( kind = 8 ), dimension(*) :: f
  integer :: i
  integer, dimension (5,node_num) :: indx
  integer :: ip
  integer :: iq
  real ( kind = 8 ) :: p_mean
  logical, parameter :: verbose = .false.
!
!  Compute the mean pressure.
!
!  To do this requires computing the integral of pressure over the region.
!
!  To compute that integral, we compute the integral of pressure over
!  each element.
!
  p_mean = 0.0D+00
  area_region = 0.0D+00

  do element = 1, element_num

    do iq = 1, 3

      ip = element_node(iq,element)
      i = indx(5,ip)

      if ( 0 < i ) then
        p_mean = p_mean + f(i) * area(element) / 3.0D+00
      else if ( i == 0 ) then
        p_mean = p_mean + 0.0D+00 * area(element) / 3.0D+00
      end if

    end do

    area_region = area_region + area(element)

  end do

  p_mean = p_mean / area_region

  if ( verbose ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) 'P_NORMALIZE: the mean pressure = ', p_mean
  end if
!
!  Modify the finite element coefficients associated with pressure,
!  by subtracting off the mean pressure.
!
!  The resulting data will now have a mean pressure of 0.
!
  do ip = 1, node_num
    i = indx(5,ip)
    if ( 0 < i ) then
      f(i) = f(i) - p_mean
    end if
  end do

  return
end
subroutine p_write ( f, indx, p_file_name, element_num, node_num, element_node )

!*****************************************************************************80
!
!! P_WRITE writes pressure data to files.
!
!  Discussion:
!
!    We want to output the pressure at every node.  At vertices,
!    this value can be read out from the F array.  However, at
!    mid side nodes, pressure is not an unknown.  We get the value
!    of pressure at a midside node by averaging the values at the
!    two incident vertices.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 September 2005
!
!  Parameters:
!
!    Input, real ( kind = 8 ) F(*), the finite element coefficients.
!
!    Input, integer INDX(5,NODE_NUM), at each node, lists the
!    index of the finite element coefficients for U, V, B1, B2, and P.
!
!    Input, character ( len = * ) P_FILE_NAME, the name of the
!    output file containing the pressure field values.
!
!    Input, integer ELEMENT_NUM, the number of elements.
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, integer ELEMENT_NODE(6,ELEMENT_NUM), the nodes of each element.
!
  implicit none

  integer :: element_num
  integer :: node_num

  integer element
  integer, dimension ( 6, element_num ) :: element_node
  real ( kind = 8 ), dimension(*) :: f
  integer :: i1
  integer :: i2
  integer :: i3
  integer, dimension(5,node_num) :: indx
  integer :: ip
  integer iunit
  integer :: n1
  integer :: n2
  integer :: n3
  integer :: n4
  integer :: n5
  integer :: n6
  real ( kind = 8 ), dimension(node_num) :: p
  character ( len = * ) p_file_name
!
!  Extract or compute the pressure at every node.
!
  p(1:node_num) = 0.0D+00

  do element = 1, element_num

    n1 = element_node(1,element)
    n2 = element_node(2,element)
    n3 = element_node(3,element)
    n4 = element_node(4,element)
    n5 = element_node(5,element)
    n6 = element_node(6,element)

    i1 = indx(5,n1)
    i2 = indx(5,n2)
    i3 = indx(5,n3)

    p(n1) = f(i1)
    p(n2) = f(i2)
    p(n3) = f(i3)

    p(n4) = 0.5D+00 * ( p(n1) + p(n2) )
    p(n5) = 0.5D+00 * ( p(n2) + p(n3) )
    p(n6) = 0.5D+00 * ( p(n3) + p(n1) )

  end do
!
!  Write the pressures to a file.
!
  call get_unit ( iunit )

  open ( unit = iunit, file = p_file_name, status = 'replace' )

  do ip = 1, node_num

    write ( iunit, '(2x,g14.6)' ) p(ip)

  end do

  close ( unit = iunit )

  return
end
subroutine quad_rule ( element_max, element_num, element_node, node_num, &
  xc, yc, xq, yq )

!*****************************************************************************80
!
!! QUAD_RULE sets the weights and abscissas for quadrature.
!
!  Discussion:
!
!    The quadrature rule defined here is a simple 3 point rule
!    with equal weights.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 August 2005
!
!  Parameters:
!
!    Input, integer ELEMENT_MAX, the maximum number of elements.
!
!    Input, integer ELEMENT_NUM, the number of elements.
!
!    Input, integer ELEMENT_NODE(6,ELEMENT_NUM), the nodes of each element.
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, real ( kind = 8 ) XC(NODE_NUM), YC(NODE_NUM), the coordinates
!    of nodes.
!
!    Output, real ( kind = 8 ) XQ(ELEMENT_MAX,3), YQ(ELEMENT_MAX,3),
!    the coordinates of the quadrature points.
!
  implicit none

  integer element_max
  integer element_num
  integer node_num

  integer element
  integer, dimension ( 6, element_num ) :: element_node
  integer ip1
  integer ip2
  integer ip3
  real ( kind = 8 ), dimension(node_num) :: xc
  real ( kind = 8 ), dimension(element_max,3) :: xq
  real ( kind = 8 ), dimension(node_num) :: yc
  real ( kind = 8 ), dimension(element_max,3) :: yq

  do element = 1, element_num

    ip1 = element_node(1,element)
    ip2 = element_node(2,element)
    ip3 = element_node(3,element)

    xq(element,1) = 0.5D+00 * ( xc(ip1) + xc(ip2) )
    xq(element,2) = 0.5D+00 * ( xc(ip2) + xc(ip3) )
    xq(element,3) = 0.5D+00 * ( xc(ip3) + xc(ip1) )
    yq(element,1) = 0.5D+00 * ( yc(ip1) + yc(ip2) )
    yq(element,2) = 0.5D+00 * ( yc(ip2) + yc(ip3) )
    yq(element,3) = 0.5D+00 * ( yc(ip3) + yc(ip1) )

  end do

  return
end
function rhs ( id, x, y  )

!*****************************************************************************80
!
!! RHS sets the right hand side of each equation.
!
!  Discussion:
!
!    The vector right hand side for the B equation must be
!    divergence free.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 August 2005
!
!  Parameters:
!
!    Input, integer ID, identifies the equation whose right hand side
!    function is to be evaluated:
!    1, horizontal velocity (momentum) equation;
!    2, vertical velocity (momentum) equation;
!    3, horizontal magnetic field equation;
!    4, vertical magnetic field equation;
!    5, pressure (continuity) equation;
!
!    Input, real ( kind = 8 ) X, Y, the coordinates of the point at which
!    the right hand side function is to be evaluated.
!
!    Output, real ( kind = 8 ) RHS, the right hand side function.
!
  implicit none

  real ( kind = 8 ) :: dphidx
  real ( kind = 8 ) :: dphidy
  integer :: id
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) :: rhs
  real ( kind = 8 ) :: x
  real ( kind = 8 ) :: y
!
!  phi(x,y) = y^2 * sin ( pi * x )
!
  dphidx = pi * y * y * cos ( pi * x )
  dphidy = 2.0D+00 * y * sin( pi * x )
!
!  U equation:
!
  if ( id == 1 ) then

    rhs = 0.0D+00
!
!  V equation:
!
  else if ( id == 2 ) then

    rhs = 0.0D+00
!
!  B1 equation:
!
  else if ( id == 3 ) then

    rhs = dphidy
!
!  B2 equation:
!
  else if ( id == 4) then

    rhs = -dphidx

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
  integer d
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer n
  integer s
  integer values(8)
  integer y

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
function ubdry ( id, y, ylength, y_inflow, y_outflow  )

!*****************************************************************************80
!
!! UBDRY sets the parabolic inflow for velocity.
!
!  Discussion:
!
!    The inflow is parabolic.  It is zero at Y = 0,
!    1 at Y = Y_INFLOW/2, and 0 again at Y = Y_INFLOW.
!
!    The outflow is parabolic.  It is zero at Y = Y_OUTFLOW,
!    -1 at Y = (Y_OUTFLOW+YLENGTH)/2, and 0 again at Y = YLENGTH.
!
!    +---------+
!    |         ->
!    |         -->
!    +         ->
!    ->        +
!    -->       |
!    ->        |
!    +---------+
!
!    Actually, it is typical to specify the inflow, but not the outflow.
!    A different boundary condition may be used to handle the outflow.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 August 2005
!
!  Parameters:
!
!    Input, integer ID, is -1 if the inflow function is to be evaluated,
!    and any other value to evaluate the outflow.
!
!    Input, real ( kind = 8 ) Y, the Y coordinate of a point for which
!    the inflow or outflow function is to be evaluated.
!
!    Input, real ( kind = 8 ) YLENGTH, the maximum Y coordinate at which
!    there is outflow.
!
!    Input, real ( kind = 8 ) Y_INFLOW, the maximum Y coordinate at which
!    there is inflow.
!
!    Input, real ( kind = 8 ) Y_OUTFLOW, the minimum Y coordinate at which
!    there is outflow.
!
  implicit none

  integer :: id
  real ( kind = 8 ) :: ubdry
  real ( kind = 8 ) :: y
  real ( kind = 8 ) :: y_inflow
  real ( kind = 8 ) :: y_outflow
  real ( kind = 8 ) :: ylength

  if ( id == -1 ) then

    ubdry = 4.0D+00 * y * ( y_inflow - y ) / y_inflow**2

  else

    ubdry = -4.0D+00 * ( y - ylength ) * ( y - y_outflow ) &
      / ( ylength - y_outflow )**2

  end if

  return
end
subroutine uv_write ( f, indx, node_num, uv_file_name, &
  xc, yc, ylength, y_inflow, y_outflow )

!*****************************************************************************80
!
!! UV_WRITE writes velocity data to files.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 August 2005
!
!  Parameters:
!
!    Input, real ( kind = 8 ) F(*), the current solution.
!
!    Input, integer INDX(5,NODE_NUM), at each node, lists the
!    index of the finite element coefficients for U, V, B1, B2, and P.
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, character ( len = * ) UV_FILE_NAME, the name of the
!    output file containing the velocity field values.
!
!    Input, real ( kind = 8 ) XC(NODE_NUM), YC(NODE_NUM), the coordinates
!    of nodes.
!
!    Input, real ( kind = 8 ) YLENGTH, the maximum Y coordinate at which
!    there is outflow .
!
!    Input, real ( kind = 8 ) Y_INFLOW, the maximum Y coordinate at which
!    there is inflow.
!
!    Input, real ( kind = 8 ) Y_OUTFLOW, the minimum Y coordinate at which
!    there is outflow.
!
  implicit none

  integer :: node_num

  real ( kind = 8 ), dimension(*) :: f
  integer, dimension(5,node_num) :: indx
  integer :: ip
  integer :: iuk_u
  integer :: iuk_v
  integer iunit
  real ( kind = 8 ) :: u_value
  real ( kind = 8 ) :: ubdry
  character ( len = * ) uv_file_name
  real ( kind = 8 ) :: v_value
  real ( kind = 8 ) :: x
  real ( kind = 8 ), dimension(node_num) :: xc
  real ( kind = 8 ) :: y
  real ( kind = 8 ) :: y_inflow
  real ( kind = 8 ) :: y_outflow
  real ( kind = 8 ), dimension(node_num) :: yc
  real ( kind = 8 ) :: ylength

  call get_unit ( iunit )

  open ( unit = iunit, file = uv_file_name, status = 'replace' )

  do ip = 1, node_num

    y = yc(ip)
    x = xc(ip)
    iuk_u = indx(1,ip)
    iuk_v = indx(2,ip)

    if ( iuk_u < 0 ) then
      u_value = ubdry ( iuk_u,  y, ylength, y_inflow, y_outflow )
    else if ( iuk_u == 0 ) then
      u_value = 0.0D+00
    else
      u_value = f(iuk_u)
    end if

    if ( 0 < iuk_v ) then
      v_value = f(iuk_v)
    else
      v_value = 0.0D+00
    end if

    write ( iunit, '(2x,g14.6,2x,g14.6)' ) u_value, v_value

  end do

  close ( unit = iunit )

  return
end
subroutine xy_set ( node_num, nx, ny, xlength, ylength, xc, yc )

!*****************************************************************************80
!
!! XY_SET sets the node coordinates.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 August 2005
!
!  Parameters:
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, integer NX, NY, the number of elements in the X and Y directions.
!
!    Input, real ( kind = 8 ) XLENGTH, YLENGTH, the maximum X and Y coordinates.
!
!    Output, real ( kind = 8 ) XC(NODE_NUM), YC(NODE_NUM), the coordinates
!    of nodes.
!
  implicit none

  integer :: node_num

  integer :: col
  integer :: col_num
  integer :: node
  integer :: nx
  integer :: ny
  integer :: row
  integer :: row_num
  real ( kind = 8 ), dimension(node_num) :: xc
  real ( kind = 8 ) :: xlength
  real ( kind = 8 ), dimension(node_num) :: yc
  real ( kind = 8 ) :: ylength
!
!  Set the coordinates of the nodes.
!
  col_num = 2 * nx - 1
  row_num = 2 * ny - 1

  node = 0

  do col = 1, col_num

    do row = 1, row_num

      node = node + 1

      xc(node) = xlength * real ( col - 1,     kind = 8 ) &
                         / real ( col_num - 1, kind = 8 )

      yc(node) = ylength * real ( row - 1,     kind = 8 ) &
                         / real ( row_num - 1, kind = 8 )

    end do

  end do

  return
end
subroutine xy_write ( xy_file_name, node_num, xc, yc )

!*****************************************************************************80
!
!! XY_WRITE writes the coordinate data to a file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 September 2005
!
!  Parameters:
!
!    Input, character ( len = * ) XY_FILE_NAME, the name of the
!    output file containing the node coordinate values.
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, real ( kind = 8 ) XC(NODE_NUM), YC(NODE_NUM), the coordinates
!    of nodes.
!
  implicit none

  integer :: node_num

  integer iunit
  integer node
  real ( kind = 8 ), dimension(node_num) :: xc
  character ( len = * ) xy_file_name
  real ( kind = 8 ), dimension(node_num) :: yc

  call get_unit ( iunit )

  open ( unit = iunit, file = xy_file_name, status = 'replace' )

  do node = 1, node_num

    write ( iunit, '(2x,f12.6,2x,f12.6)' ) xc(node), yc(node)

  end do

  close ( unit = iunit )

  return
end
