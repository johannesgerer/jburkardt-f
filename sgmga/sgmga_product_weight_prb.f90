program main

!*****************************************************************************80
!
!! MAIN is the main program for SGMGA_PRODUCT_WEIGHT_PRB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SGMGA_PRODUCT_WEIGHT_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the SGMGA_PRODUCT_WEIGHT function.'

  call sgmga_product_weight_tests ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SGMGA_PRODUCT_WEIGHT_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine sgmga_product_weight_tests ( )

!****************************************************************************80
!
!! SGMGA_PRODUCT_WEIGHT_TESTS calls SGMGA_PRODUCT_WEIGHT_TEST.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 June 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ), allocatable :: np(:)
  integer ( kind = 4 ) np_sum
  integer ( kind = 4 ), allocatable, dimension ( : ) :: order_1d
  integer ( kind = 4 ) order_nd
  real ( kind = 8 ), allocatable :: p(:)
  integer ( kind = 4 ), allocatable, dimension ( : ) :: rule

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SGMGA_PRODUCT_WEIGHT_TESTS'
  write ( *, '(a)' ) '  Call SGMGA_PRODUCT_WEIGHT_TEST with various arguments.'

  dim_num = 2
  allocate ( order_1d(1:dim_num) )
  order_1d = (/ 3, 5 /)
  order_nd = product ( order_1d(1:dim_num) )
  allocate ( rule(1:dim_num) )
  rule = (/ 1, 1 /)
  allocate ( np(1:dim_num) )
  np = (/ 0, 0 /)
  np_sum = sum ( np(1:dim_num) )
  allocate ( p(1:np_sum) )
  call sgmga_product_weight_test ( dim_num, order_1d, order_nd, rule, np, p )
  deallocate ( np )
  deallocate ( order_1d )
  deallocate ( p )
  deallocate ( rule )

  dim_num = 2
  allocate ( order_1d(1:dim_num) )
  order_1d = (/ 3, 7 /)
  order_nd = product ( order_1d(1:dim_num) )
  allocate ( rule(1:dim_num) )
  rule = (/ 1, 5 /)
  allocate ( np(1:dim_num) )
  np = (/ 0, 0 /)
  np_sum = sum ( np(1:dim_num) )
  allocate ( p(1:np_sum) )
  call sgmga_product_weight_test ( dim_num, order_1d, order_nd, rule, np, p )
  deallocate ( np )
  deallocate ( order_1d )
  deallocate ( p )
  deallocate ( rule )

  dim_num = 2
  allocate ( order_1d(1:dim_num) )
  order_1d = (/ 3, 3 /)
  order_nd = product ( order_1d(1:dim_num) )
  allocate ( rule(1:dim_num) )
  rule = (/ 3, 7 /)
  allocate ( np(1:dim_num) )
  np = (/ 0, 0 /)
  np_sum = sum ( np(1:dim_num) )
  allocate ( p(1:np_sum) )
  call sgmga_product_weight_test ( dim_num, order_1d, order_nd, rule, np, p )
  deallocate ( np )
  deallocate ( order_1d )
  deallocate ( p )
  deallocate ( rule )

  dim_num = 2
  allocate ( order_1d(1:dim_num) )
  order_1d = (/ 5, 5 /)
  order_nd = product ( order_1d(1:dim_num) )
  allocate ( rule(1:dim_num) )
  rule = (/ 1, 8 /)
  allocate ( np(1:dim_num) )
  np = (/ 0, 1 /)
  np_sum = sum ( np(1:dim_num) )
  allocate ( p(1:np_sum) )
  p(1) = 1.5D+00
  call sgmga_product_weight_test ( dim_num, order_1d, order_nd, rule, np, p )
  deallocate ( np )
  deallocate ( order_1d )
  deallocate ( p )
  deallocate ( rule )

  dim_num = 2
  allocate ( order_1d(1:dim_num) )
  order_1d = (/ 5, 5 /)
  order_nd = product ( order_1d(1:dim_num) )
  allocate ( rule(1:dim_num) )
  rule = (/ 2, 9 /)
  allocate ( np(1:dim_num) )
  np = (/ 0, 2 /)
  np_sum = sum ( np(1:dim_num) )
  allocate ( p(1:np_sum) )
  p(1) = 0.5D+00
  p(2) = 1.5D+00
  call sgmga_product_weight_test ( dim_num, order_1d, order_nd, rule, np, p )
  deallocate ( np )
  deallocate ( order_1d )
  deallocate ( p )
  deallocate ( rule )

  dim_num = 2
  allocate ( order_1d(1:dim_num) )
  order_1d = (/ 7, 9 /)
  order_nd = product ( order_1d(1:dim_num) )
  allocate ( rule(1:dim_num) )
  rule = (/ 6, 10 /)
  allocate ( np(1:dim_num) )
  np = (/ 1, 0 /)
  np_sum = sum ( np(1:dim_num) )
  allocate ( p(1:np_sum) )
  p(1) = 2.0D+00
  call sgmga_product_weight_test ( dim_num, order_1d, order_nd, rule, np, p )
  deallocate ( np )
  deallocate ( order_1d )
  deallocate ( p )
  deallocate ( rule )

  dim_num = 3
  allocate ( order_1d(1:dim_num) )
  order_1d = (/ 2, 3, 3 /)
  order_nd = product ( order_1d(1:dim_num) )
  allocate ( rule(1:dim_num) )
  rule = (/ 1, 4, 5 /)
  allocate ( np(1:dim_num) )
  np = (/ 0, 0, 0 /)
  np_sum = sum ( np(1:dim_num) )
  allocate ( p(1:np_sum) )
  call sgmga_product_weight_test ( dim_num, order_1d, order_nd, rule, np, p )
  deallocate ( np )
  deallocate ( order_1d )
  deallocate ( p )
  deallocate ( rule )

  return
end
subroutine sgmga_product_weight_test ( dim_num, order_1d, order_nd, &
  rule, np, p )

!****************************************************************************80
!
!! SGMGA_PRODUCT_WEIGHT_TEST: mixed factor product rule test.
!
!  Discussion:
!
!    This routine gets the sparse grid indices and determines the
!    corresponding sparse grid weights.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) ORDER_1D(DIM_NUM), the order of the 1D rules.
!
!    Input, integer ( kind = 4 ) ORDER_ND, the order of the product rule.
!
!    Input, integer ( kind = 4 ) RULE(DIM_NUM), the rule in each dimension.
!     1, "CC",  Clenshaw Curtis, Closed Fully Nested.
!     2, "F2",  Fejer Type 2, Open Fully Nested.
!     3, "GP",  Gauss Patterson, Open Fully Nested.
!     4, "GL",  Gauss Legendre, Open Weakly Nested.
!     5, "GH",  Gauss Hermite, Open Weakly Nested.
!     6, "GGH", Generalized Gauss Hermite, Open Weakly Nested.
!     7, "LG",  Gauss Laguerre, Open Non Nested.
!     8, "GLG", Generalized Gauss Laguerre, Open Non Nested.
!     9, "GJ",  Gauss Jacobi, Open Non Nested.
!    10, "HGK", Hermite Genz-Keister, Open Fully Nested.
!    11, "UO",  User supplied Open, presumably Non Nested.
!    12, "UC",  User supplied Closed, presumably Non Nested.
!
!    Input, integer ( kind = 4 ) NP(DIM_NUM), the number of parameters
!    used by each rule.
!
!    Input, real ( kind = 8 ) P(*), the parameters needed by each rule.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) order_nd

  real ( kind = 8 ) alpha
  real ( kind = 8 ) arg1
  real ( kind = 8 ) arg2
  real ( kind = 8 ) arg3
  real ( kind = 8 ) arg4
  real ( kind = 8 ) beta
  integer ( kind = 4 ) dim
  integer ( kind = 4 ) np(dim_num)
  integer ( kind = 4 ) order_1d(dim_num)
  real ( kind = 8 ) p(*)
  integer ( kind = 4 ) p_index
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r8_gamma
  integer ( kind = 4 ) rule(dim_num)
  real ( kind = 8 ) value1
  real ( kind = 8 ) value2
  real ( kind = 8 ) weight(order_nd)
  real ( kind = 8 ) weight_sum
  real ( kind = 8 ) weight_sum_error
  real ( kind = 8 ) weight_sum_exact

  weight_sum_exact = 1.0D+00

  p_index = 1

  do dim = 1, dim_num

    if ( rule(dim) == 1 ) then
      weight_sum_exact = weight_sum_exact * 2.0D+00
    else if ( rule(dim) == 2 ) then
      weight_sum_exact = weight_sum_exact * 2.0D+00
    else if ( rule(dim) == 3 ) then
      weight_sum_exact = weight_sum_exact * 2.0D+00
    else if ( rule(dim) == 4 ) then
      weight_sum_exact = weight_sum_exact * 2.0D+00
    else if ( rule(dim) == 5 ) then
      weight_sum_exact = weight_sum_exact * sqrt ( pi )
    else if ( rule(dim) == 6 ) then
      alpha = p(p_index)
      weight_sum_exact = weight_sum_exact * r8_gamma ( 0.5D+00 *( alpha + 1.0D+00 ) )
    else if ( rule(dim) == 7 ) then
      weight_sum_exact = weight_sum_exact * 1.0D+00
    else if ( rule(dim) == 8 ) then
      alpha = p(p_index)
      weight_sum_exact = weight_sum_exact * r8_gamma ( alpha + 1.0D+00 )
    else if ( rule(dim) == 9 ) then
      alpha = p(p_index)
      beta = p(p_index+1)
      arg1 = - alpha
      arg2 = 1.0D+00
      arg3 = beta + 2.0D+00
      arg4 = - 1.0D+00
      call r8_hyper_2f1 ( arg1, arg2, arg3, arg4, value1 )
      arg1 = - beta
      arg2 = 1.0D+00
      arg3 = alpha + 2.0D+00
      arg4 = - 1.0D+00
      call r8_hyper_2f1 ( arg1, arg2, arg3, arg4, value2 )
      weight_sum_exact = weight_sum_exact * ( &
        value1 / ( beta + 1.0D+00 ) + value2 / ( alpha + 1.0D+00 ) )
    else if ( rule(dim) == 10 ) then
      weight_sum_exact = weight_sum_exact * sqrt ( pi )
    else if ( rule(dim) == 11 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SGMGA_PRODUCT_WEIGHT_TEST - Fatal error!'
      write ( *, '(a)' ) '  Do not know how to handle rule 11.'
      stop
    else if ( rule(dim) == 12 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SGMGA_PRODUCT_WEIGHT_TEST - Fatal error!'
      write ( *, '(a)' ) '  Do not know how to handle rule 12.'
      stop
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SGMGA_PRODUCT_WEIGHT_TEST - Fatal error!'
      write ( *, '(a,i8)' ) '  Unexpected value of RULE = ', rule(dim)
      stop
    end if

    p_index = p_index + np(dim)

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SGMGA_PRODUCT_WEIGHT_TEST:'
  write ( *, '(a)' ) '  Compute the weights of a mixed factor product grid.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  As a simple test, sum these weights.'
  write ( *, '(a,g24.16)' ) '  They should sum to exactly ', weight_sum_exact
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Spatial dimension DIM_NUM = ', dim_num
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' Dimension      Rule     Order        Parameters'
  write ( *, '(a)' ) ' '

  p_index = 1

  do dim = 1, dim_num

    if ( rule(dim) == 1 ) then
      write ( *, '(2x,i8,2x,i8,2x,i11,2x,g14.6,2x,g14.6)' ) &
        dim, rule(dim), order_1d(dim)
    else if ( rule(dim) ==  2 ) then
      write ( *, '(2x,i8,2x,i8,2x,i11,2x,g14.6,2x,g14.6)' ) &
        dim, rule(dim), order_1d(dim)
    else if ( rule(dim) ==  3 ) then
      write ( *, '(2x,i8,2x,i8,2x,i11,2x,g14.6,2x,g14.6)' ) &
        dim, rule(dim), order_1d(dim)
    else if ( rule(dim) ==  4 ) then
      write ( *, '(2x,i8,2x,i8,2x,i11,2x,g14.6,2x,g14.6)' ) &
        dim, rule(dim), order_1d(dim)
    else if ( rule(dim) ==  5 ) then
      write ( *, '(2x,i8,2x,i8,2x,i11,2x,g14.6,2x,g14.6)' ) &
        dim, rule(dim), order_1d(dim)
    else if ( rule(dim) ==  6 ) then
      write ( *, '(2x,i8,2x,i8,2x,i11,2x,g14.6,2x,g14.6)' ) &
        dim, rule(dim), order_1d(dim), p(p_index)
    else if ( rule(dim) ==  7 ) then
      write ( *, '(2x,i8,2x,i8,2x,i11,2x,g14.6,2x,g14.6)' ) &
        dim, rule(dim), order_1d(dim)
    else if ( rule(dim) ==  8 ) then
      write ( *, '(2x,i8,2x,i8,2x,i11,2x,g14.6,2x,g14.6)' ) &
        dim, rule(dim), order_1d(dim), p(p_index)
    else if ( rule(dim) ==  9 ) then
      write ( *, '(2x,i8,2x,i8,2x,i11,2x,g14.6,2x,g14.6)' ) &
        dim, rule(dim), order_1d(dim), p(p_index), p(p_index+1)
    else if ( rule(dim) == 10 ) then
      write ( *, '(2x,i8,2x,i8,2x,i11,2x,g14.6,2x,g14.6)' ) &
        dim, rule(dim), order_1d(dim)
    else if ( rule(dim) == 11 ) then
      write ( *, '(2x,i8,2x,i8,2x,i11,2x,g14.6,2x,g14.6)' ) &
        dim, rule(dim), order_1d(dim)
    end if

    p_index = p_index + np(dim)

  end do
!
!  Compute the weights.
!
  call sgmga_product_weight ( dim_num, order_1d, order_nd, rule, &
    np, p, weight )
!
!  Sum the weights.
!
  weight_sum = sum ( weight(1:order_nd) )

  weight_sum_error = abs ( weight_sum - weight_sum_exact )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    Weight sum  Expected sum    Difference'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6)' ) &
    weight_sum, weight_sum_exact, weight_sum_error

  return
end

