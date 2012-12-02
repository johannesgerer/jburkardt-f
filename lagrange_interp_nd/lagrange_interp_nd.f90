subroutine cc_compute_points ( n, points )

!*****************************************************************************80
!
!! CC_COMPUTE_POINTS: abscissas of a Clenshaw Curtis rule.
!
!  Discussion:
!
!    Our convention is that the abscissas are numbered from left to right.
!
!    The rule is defined on [-1,1].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order.
!    1 <= N.
!
!    Output, real ( kind = 8 ) POINTS(N), the abscissas.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) points(n)

  if ( n < 1 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CC_COMPUTE_POINTS - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of N = ', n
    stop

  else if ( n == 1 ) then

    points(1) = 0.0D+00

  else

    do i = 1, n
      points(i) = cos ( real ( n - i, kind = 8 ) * pi &
                      / real ( n - 1, kind = 8 ) )
    end do

    points(1) = -1.0D+00
    if ( mod ( n, 2 ) == 1 ) then
      points((n+1)/2) = 0.0D+00
    end if
    points(n) = +1.0D+00

  end if

  return
end
subroutine lagrange_basis_1d ( nd, xd, ni, xi, lb ) 

!*****************************************************************************80
!
!! LAGRANGE_BASIS_1D evaluates a 1D Lagrange basis.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ND, the number of data points.
!
!    Input, real ( kind = 8 ) XD(ND), the interpolation nodes.
!
!    Input, integer ( kind = 4 ) NI, the number of evaluation points.
!
!    Input, real ( kind = 8 ) XI(NI), the evaluation points.
!
!    Output, real ( kind = 8 ) LB(NI,ND), the value, at the I-th point XI, 
!    of the Jth basis function.
!
  implicit none

  integer ( kind = 4 ) nd
  integer ( kind = 4 ) ni

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) lb(ni,nd)
  real ( kind = 8 ) xd(nd)
  real ( kind = 8 ) xi(ni)
  
  do i = 1, ni
    do j = 1, nd
      lb(i,j) = product ( ( xi(i) - xd(1:j-1)  ) / ( xd(j) - xd(1:j-1)  ) ) &
              * product ( ( xi(i) - xd(j+1:nd) ) / ( xd(j) - xd(j+1:nd) ) )
    end do
  end do

  return
end
subroutine lagrange_interp_nd_grid ( m, n_1d, a, b, nd, xd )

!*****************************************************************************80
!
!! LAGRANGE_INTERP_ND_GRID sets an M-dimensional Lagrange interpolant grid.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N_1D(M), the order of the 1D rule to be used
!    in each dimension.
!
!    Input, real ( kind = 8 ) A(M), B(M), the lower and upper limits.
!
!    Input, integer ( kind = 4 ) ND, the number of points in the product grid.
!
!    Output, real ( kind = 8 ) XD(M,ND), the points at which data was sampled.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) nd

  real ( kind = 8 ) a(m)
  real ( kind = 8 ) b(m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_1d(m)
  real ( kind = 8 ), allocatable :: x_1d(:)
  real ( kind = 8 ) xd(m,nd)
!
!  Compute the data points.
!
  xd(1:m,1:nd) = 0.0D+00
  do i = 1, m
    n = n_1d(i)
    allocate ( x_1d(1:n) )
    call cc_compute_points ( n, x_1d )
    x_1d(1:n) = 0.5D+00 * ( ( 1.0D+00 - x_1d(1:n) ) * a(i) &
                          + ( 1.0D+00 + x_1d(1:n) ) * b(i) )
    call r8vec_direct_product ( i, n, x_1d, m, nd, xd )
    deallocate ( x_1d )
  end do

  return
end
subroutine lagrange_interp_nd_grid2 ( m, ind, a, b, nd, xd )

!*****************************************************************************80
!
!! LAGRANGE_INTERP_ND_GRID2 sets an M-dimensional Lagrange interpolant grid.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) IND(M), the index or level of the 1D rule 
!    to be used in each dimension.
!
!    Input, real ( kind = 8 ) A(M), B(M), the lower and upper limits.
!
!    Input, integer ( kind = 4 ) ND, the number of points in the product grid.
!
!    Output, real ( kind = 8 ) XD(M,ND), the points at which data was sampled.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) nd

  real ( kind = 8 ) a(m)
  real ( kind = 8 ) b(m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ind(m)
  integer ( kind = 4 ) n
  real ( kind = 8 ), allocatable :: x_1d(:)
  real ( kind = 8 ) xd(m,nd)
!
!  Compute the data points.
!
  xd(1:m,1:nd) = 0.0D+00
  do i = 1, m
    call order_from_level_135 ( ind(i), n )
    allocate ( x_1d(1:n) )
    call cc_compute_points ( n, x_1d )
    x_1d(1:n) = 0.5D+00 * ( ( 1.0D+00 - x_1d(1:n) ) * a(i) &
                          + ( 1.0D+00 + x_1d(1:n) ) * b(i) )
    call r8vec_direct_product ( i, n, x_1d, m, nd, xd )
    deallocate ( x_1d )
  end do

  return
end
subroutine lagrange_interp_nd_size ( m, n_1d, nd )

!*****************************************************************************80
!
!! LAGRANGE_INTERP_ND_SIZE sizes an M-dimensional Lagrange interpolant.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N_1D(M), the order of the 1D rule to be used
!    in each dimension.
!
!    Output, integer ( kind = 4 ) ND, the number of points in the product grid.
!
  implicit none

  integer ( kind = 4 ) m

  integer ( kind = 4 ) n_1d(m)
  integer ( kind = 4 ) nd
!
!  Determine the number of data points.
!
  nd = product ( n_1d(1:m) )

  return
end
subroutine lagrange_interp_nd_size2 ( m, ind, nd )

!*****************************************************************************80
!
!! LAGRANGE_INTERP_ND_SIZE2 sizes an M-dimensional Lagrange interpolant.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) IND(M), the index or level of the 1D rule 
!    to be used in each dimension.
!
!    Output, integer ( kind = 4 ) ND, the number of points in the product grid.
!
  implicit none

  integer ( kind = 4 ) m

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ind(m)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nd
!
!  Determine the number of data points.
!
  nd = 1
  do i = 1, m
    call order_from_level_135 ( ind(i), n )
    nd = nd * n
  end do

  return
end
subroutine lagrange_interp_nd_value ( m, n_1d, a, b, nd, zd, ni, xi, zi )

!*****************************************************************************80
!
!! LAGRANGE_INTERP_ND_VALUE evaluates an ND Lagrange interpolant.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N_1D(M), the order of the 1D rule to be used
!    in each dimension.
!
!    Input, real ( kind = 8 ) A(M), B(M), the lower and upper limits.
!
!    Input, integer ( kind = 4 ) ND, the number of points in the product grid.
!
!    Input, real ( kind = 8 ) ZD(ND), the function evaluated at the points XD.
!
!    Input, integer ( kind = 4 ) NI, the number of points at which the 
!    interpolant is to be evaluated.
!
!    Input, real ( kind = 8 ) XI(M,NI), the points at which the interpolant is 
!    to be evaluated.
!
!    Output, real ( kind = 8 ) ZI(NI), the interpolant evaluated at the 
!    points XI.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) nd
  integer ( kind = 4 ) ni

  real ( kind = 8 ) a(m)
  real ( kind = 8 ) b(m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_1d(m)
  real ( kind = 8 ), allocatable :: value(:)
  real ( kind = 8 ) w(nd)
  real ( kind = 8 ), allocatable :: x_1d(:)
  real ( kind = 8 ) xi(m,ni)
  real ( kind = 8 ) zd(nd)
  real ( kind = 8 ) zi(ni)

  do j = 1, ni

    w(1:nd) = 1.0D+00

    do i = 1, m
      n = n_1d(i)
      allocate ( x_1d(1:n) )
      allocate ( value(1:n) )
      call cc_compute_points ( n, x_1d )
      x_1d(1:n) = 0.5D+00 * ( ( 1.0D+00 - x_1d(1:n) ) * a(i) &
                            + ( 1.0D+00 + x_1d(1:n) ) * b(i) )
      call lagrange_basis_1d ( n, x_1d, 1, xi(i,j), value )
      call r8vec_direct_product2 ( i, n, value, m, nd, w )
      deallocate ( value )
      deallocate ( x_1d )
    end do

    zi(j) = dot_product ( w, zd )

  end do

  return
end
subroutine lagrange_interp_nd_value2 ( m, ind, a, b, nd, zd, ni, xi, zi )

!*****************************************************************************80
!
!! LAGRANGE_INTERP_ND_VALUE2 evaluates an ND Lagrange interpolant.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) IND(M), the index or level of the 1D rule 
!    to be used in each dimension.
!
!    Input, real ( kind = 8 ) A(M), B(M), the lower and upper limits.
!
!    Input, integer ( kind = 4 ) ND, the number of points in the product grid.
!
!    Input, real ( kind = 8 ) ZD(ND), the function evaluated at the points XD.
!
!    Input, integer ( kind = 4 ) NI, the number of points at which the 
!    interpolant is to be evaluated.
!
!    Input, real ( kind = 8 ) XI(M,NI), the points at which the interpolant 
!    is to be evaluated.
!
!    Output, real ( kind = 8 ) ZI(NI), the interpolant evaluated at the 
!    points XI.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) nd
  integer ( kind = 4 ) ni

  real ( kind = 8 ) a(m)
  real ( kind = 8 ) b(m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ind(m)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  real ( kind = 8 ), allocatable :: value(:)
  real ( kind = 8 ) w(nd)
  real ( kind = 8 ), allocatable :: x_1d(:)
  real ( kind = 8 ) xi(m,ni)
  real ( kind = 8 ) zd(nd)
  real ( kind = 8 ) zi(ni)

  do j = 1, ni

    w(1:nd) = 1.0D+00

    do i = 1, m
      call order_from_level_135 ( ind(i), n )
      allocate ( x_1d(1:n) )
      allocate ( value(1:n) )
      call cc_compute_points ( n, x_1d )
      x_1d(1:n) = 0.5D+00 * ( ( 1.0D+00 - x_1d(1:n) ) * a(i) &
                            + ( 1.0D+00 + x_1d(1:n) ) * b(i) )
      call lagrange_basis_1d ( n, x_1d, 1, xi(i,j), value )
      call r8vec_direct_product2 ( i, n, value, m, nd, w )
      deallocate ( value )
      deallocate ( x_1d )
    end do

    zi(j) = dot_product ( w, zd )

  end do

  return
end
subroutine order_from_level_135 ( l, n )

!*****************************************************************************80
!
!! ORDER_FROM_LEVEL_135 evaluates the 135 level-to-order relationship.
!
!  Discussion:
!
!    Clenshaw Curtis rules, and some others, often use the following
!    scheme:
!
!    L: 0  1  2  3   4   5
!    N: 1  3  5  9  17  33 ... 2^L+1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) L, the level, which should be 0 or greater.
!
!    Output, integer ( kind = 4 ) N, the order.
!
  implicit none

  integer ( kind = 4 ) l
  integer ( kind = 4 ) n

  if ( l < 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'ORDER_FROM_LEVEL_135 - Fatal error!'
    write ( *, '(a)' ) '  Illegal input value of L!'
    stop
  else if ( l == 0 ) then
    n = 1
  else
    n = ( 2 ** l ) + 1
  end if

  return
end
