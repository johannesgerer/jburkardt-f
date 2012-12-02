program main

!****************************************************************************80
!
!! MAIN is the main program for PRODUCT_MIXED_WEIGHT_PRB.
!
!  Discussion:
!
!    PRODUCT_MIXED_WEIGHT_PRB tests PRODUCT_MIXED_WEIGHT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 December 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Fabio Nobile, Raul Tempone, Clayton Webster,
!    A Sparse Grid Stochastic Collocation Method for Partial Differential
!    Equations with Random Input Data,
!    SIAM Journal on Numerical Analysis,
!    Volume 46, Number 5, 2008, pages 2309-2345.
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PRODUCT_MIXED_WEIGHT_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '

  call product_mixed_weight_tests ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PRODUCT_MIXED_WEIGHT_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  return
end
subroutine product_mixed_weight_tests ( )

!****************************************************************************80
!
!! PRODUCT_MIXED_WEIGHT_TESTS calls PRODUCT_MIXED_WEIGHT_TEST with various arguments.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 November 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable, dimension ( : ) :: alpha
  real ( kind = 8 ), allocatable, dimension ( : ) :: beta
  integer   ( kind = 4 ) dim_num
  integer   ( kind = 4 ), allocatable, dimension ( : ) :: order_1d
  integer   ( kind = 4 ) order_nd
  integer   ( kind = 4 ), allocatable, dimension ( : ) :: rule

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PRODUCT_MIXED_WEIGHT_TESTS'
  write ( *, '(a)' ) '  Call PRODUCT_MIXED_WEIGHT_TEST with various arguments.'

  dim_num = 2
  allocate ( alpha(1:dim_num) )
  allocate ( beta(1:dim_num) )
  allocate ( order_1d(1:dim_num) )
  allocate ( rule(1:dim_num) )
  order_1d = (/ 3, 5 /)
  order_nd = product ( order_1d(1:dim_num) )
  rule = (/ 1, 1 /)
  alpha = (/ 0.0D+00, 0.0D+00 /)
  beta = (/ 0.0D+00, 0.0D+00 /)
  call product_mixed_weight_test ( dim_num, order_1d, order_nd, rule, alpha, beta )
  deallocate ( alpha )
  deallocate ( beta )
  deallocate ( order_1d )
  deallocate ( rule )

  dim_num = 2
  allocate ( alpha(1:dim_num) )
  allocate ( beta(1:dim_num) )
  allocate ( order_1d(1:dim_num) )
  allocate ( rule(1:dim_num) )
  order_1d = (/ 3, 7 /)
  order_nd = product ( order_1d(1:dim_num) )
  rule = (/ 1, 5 /)
  alpha = (/ 0.0D+00, 0.0D+00 /)
  beta = (/ 0.0D+00, 0.0D+00 /)
  call product_mixed_weight_test ( dim_num, order_1d, order_nd, rule, alpha, beta )
  deallocate ( alpha )
  deallocate ( beta )
  deallocate ( order_1d )
  deallocate ( rule )

  dim_num = 2
  allocate ( alpha(1:dim_num) )
  allocate ( beta(1:dim_num) )
  allocate ( order_1d(1:dim_num) )
  allocate ( rule(1:dim_num) )
  order_1d = (/ 3, 3 /)
  order_nd = product ( order_1d(1:dim_num) )
  rule = (/ 3, 7 /)
  alpha = (/ 0.0D+00, 0.0D+00 /)
  beta = (/ 0.0D+00, 0.0D+00 /)
  call product_mixed_weight_test ( dim_num, order_1d, order_nd, rule, alpha, beta )
  deallocate ( alpha )
  deallocate ( beta )
  deallocate ( order_1d )
  deallocate ( rule )

  dim_num = 2
  allocate ( alpha(1:dim_num) )
  allocate ( beta(1:dim_num) )
  allocate ( order_1d(1:dim_num) )
  allocate ( rule(1:dim_num) )
  order_1d = (/ 5, 5 /)
  order_nd = product ( order_1d(1:dim_num) )
  rule = (/ 1, 8 /)
  alpha = (/ 0.0D+00, 1.5D+00 /)
  beta = (/ 0.0D+00, 0.0D+00 /)
  call product_mixed_weight_test ( dim_num, order_1d, order_nd, rule, alpha, beta )
  deallocate ( alpha )
  deallocate ( beta )
  deallocate ( order_1d )
  deallocate ( rule )

  dim_num = 2
  allocate ( alpha(1:dim_num) )
  allocate ( beta(1:dim_num) )
  allocate ( order_1d(1:dim_num) )
  allocate ( rule(1:dim_num) )
  order_1d = (/ 5, 5 /)
  order_nd = product ( order_1d(1:dim_num) )
  rule = (/ 2, 9 /)
  alpha = (/ 0.0D+00, 0.5D+00 /)
  beta = (/ 0.0D+00, 1.5D+00 /)
  call product_mixed_weight_test ( dim_num, order_1d, order_nd, rule, alpha, beta )
  deallocate ( alpha )
  deallocate ( beta )
  deallocate ( order_1d )
  deallocate ( rule )

  dim_num = 2
  allocate ( alpha(1:dim_num) )
  allocate ( beta(1:dim_num) )
  allocate ( order_1d(1:dim_num) )
  allocate ( rule(1:dim_num) )
  order_1d = (/ 7, 7 /)
  order_nd = product ( order_1d(1:dim_num) )
  rule = (/ 6, 4 /)
  alpha = (/ 2.0D+00, 0.0D+00 /)
  beta = (/ 0.0D+00, 0.0D+00 /)
  call product_mixed_weight_test ( dim_num, order_1d, order_nd, rule, alpha, beta )
  deallocate ( alpha )
  deallocate ( beta )
  deallocate ( order_1d )
  deallocate ( rule )

  dim_num = 3
  allocate ( alpha(1:dim_num) )
  allocate ( beta(1:dim_num) )
  allocate ( order_1d(1:dim_num) )
  allocate ( rule(1:dim_num) )
  order_1d = (/ 2, 3, 3 /)
  order_nd = product ( order_1d(1:dim_num) )
  rule = (/ 1, 3, 5 /)
  alpha = (/ 0.0D+00, 0.0D+00, 0.0D+00 /)
  beta = (/ 0.0D+00, 0.0D+00, 0.0D+00 /)
  call product_mixed_weight_test ( dim_num, order_1d, order_nd, rule, alpha, beta )
  deallocate ( alpha )
  deallocate ( beta )
  deallocate ( order_1d )
  deallocate ( rule )
!
!  Dimension 2, Rule 13
!
  dim_num = 2
  allocate ( alpha(1:dim_num) )
  allocate ( beta(1:dim_num) )
  allocate ( order_1d(1:dim_num) )
  allocate ( rule(1:dim_num) )
  order_1d = (/ 15, 15 /)
  order_nd = product ( order_1d(1:dim_num) )
  rule = (/ 13, 13 /)
  alpha = (/ 0.0D+00, 0.0D+00 /)
  beta = (/ 0.0D+00, 0.0D+00 /)
  call product_mixed_weight_test ( dim_num, order_1d, order_nd, rule, alpha, beta )
  deallocate ( alpha )
  deallocate ( beta )
  deallocate ( order_1d )
  deallocate ( rule )
!
!  Dimension 2, Rule 16
!
  dim_num = 2
  allocate ( alpha(1:dim_num) )
  allocate ( beta(1:dim_num) )
  allocate ( order_1d(1:dim_num) )
  allocate ( rule(1:dim_num) )
  order_1d = (/ 15, 15 /)
  order_nd = product ( order_1d(1:dim_num) )
  rule = (/ 16, 16 /)
  alpha = (/ 0.0D+00, 0.0D+00 /)
  beta = (/ 0.0D+00, 0.0D+00 /)
  call product_mixed_weight_test ( dim_num, order_1d, order_nd, rule, alpha, beta )
  deallocate ( alpha )
  deallocate ( beta )
  deallocate ( order_1d )
  deallocate ( rule )
!
!  Dimension 2, Rule 17
!
  dim_num = 2
  allocate ( alpha(1:dim_num) )
  allocate ( beta(1:dim_num) )
  allocate ( order_1d(1:dim_num) )
  allocate ( rule(1:dim_num) )
  order_1d = (/ 15, 15 /)
  order_nd = product ( order_1d(1:dim_num) )
  rule = (/ 17, 17 /)
  alpha = (/ 0.0D+00, 0.0D+00 /)
  beta = (/ 0.0D+00, 0.0D+00 /)
  call product_mixed_weight_test ( dim_num, order_1d, order_nd, rule, alpha, beta )
  deallocate ( alpha )
  deallocate ( beta )
  deallocate ( order_1d )
  deallocate ( rule )

  return
end
subroutine product_mixed_weight_test ( dim_num, order_1d, order_nd, rule, &
  alpha, beta )

!****************************************************************************80
!
!! PRODUCT_MIXED_WEIGHT_TEST computes the weights of a mixed factor product rule.
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
!    07 March 2011
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
!    10, "GW",  Golub Welsch, (presumed) Open Non Nested.
!    11, "CC_SE", Clenshaw Curtis Slow Exponential, Closed Fully Nested.
!    12, "F2_SE", Fejer Type 2 Slow Exponential, Closed Fully Nested.
!    13, "GP_SE", Gauss Patterson Slow Exponential, Closed Fully Nested.
!    14, "CC_ME", Clenshaw Curtis Moderate Exponential, Closed Fully Nested.
!    15, "F2_ME", Fejer Type 2 Moderate Exponential, Closed Fully Nested.
!    16, "GP_ME", Gauss Patterson Moderate Exponential, Closed Fully Nested.
!    17, "CCN", Clenshaw Curtis Nested, Linear, Closed Fully Nested rule.
!
!    Input, real ( kind = 8 ) ALPHA(DIM_NUM), BETA(DIM_NUM), parameters used for
!    Generalized Gauss Hermite, Generalized Gauss Laguerre, and Gauss Jacobi rules.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) order_nd

  real ( kind = 8 ) alpha(dim_num)
  real ( kind = 8 ) arg1
  real ( kind = 8 ) arg2
  real ( kind = 8 ) arg3
  real ( kind = 8 ) arg4
  real ( kind = 8 ) beta(dim_num)
  integer ( kind = 4 ) dim
  integer ( kind = 4 ) order_1d(dim_num)
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
      weight_sum_exact = weight_sum_exact * r8_gamma ( 0.5D+00 *( alpha(dim) + 1.0D+00 ) )
    else if ( rule(dim) == 7 ) then
      weight_sum_exact = weight_sum_exact * 1.0D+00
    else if ( rule(dim) == 8 ) then
      weight_sum_exact = weight_sum_exact * r8_gamma ( alpha(dim) + 1.0D+00 )
    else if ( rule(dim) == 9 ) then
      arg1 = - alpha(dim)
      arg2 = 1.0D+00
      arg3 = beta(dim) + 2.0D+00
      arg4 = - 1.0D+00
      call r8_hyper_2f1 ( arg1, arg2, arg3, arg4, value1 )
      arg1 = - beta(dim)
      arg2 = 1.0D+00
      arg3 = alpha(dim) + 2.0D+00
      arg4 = - 1.0D+00
      call r8_hyper_2f1 ( arg1, arg2, arg3, arg4, value2 )
      weight_sum_exact = weight_sum_exact * ( &
        value1 / ( beta(dim) + 1.0D+00 ) + value2 / ( alpha(dim) + 1.0D+00 ) )
    else if ( rule(dim) == 10 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PRODUCT_MIXED_WEIGHT_TEST - Fatal error!'
      write ( *, '(a)' ) '  Do not know how to handle rule 10.'
      stop
    else if ( rule(dim) == 11 ) then
      weight_sum_exact = weight_sum_exact * 2.0D+00
    else if ( rule(dim) == 12 ) then
      weight_sum_exact = weight_sum_exact * 2.0D+00
    else if ( rule(dim) == 13 ) then
      weight_sum_exact = weight_sum_exact * 2.0D+00
    else if ( rule(dim) == 13 ) then
      weight_sum_exact = weight_sum_exact * 2.0D+00
    else if ( rule(dim) == 14 ) then
      weight_sum_exact = weight_sum_exact * 2.0D+00
    else if ( rule(dim) == 15 ) then
      weight_sum_exact = weight_sum_exact * 2.0D+00
    else if ( rule(dim) == 16 ) then
      weight_sum_exact = weight_sum_exact * 2.0D+00
    else if ( rule(dim) == 17 ) then
      weight_sum_exact = weight_sum_exact * 2.0D+00
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PRODUCT_MIXED_WEIGHT_TEST - Fatal error!'
      write ( *, '(a,i8)' ) '  Unexpected value of RULE = ', rule(dim)
      stop
    end if

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PRODUCT_MIXED_WEIGHT_TEST:'
  write ( *, '(a)' ) '  Compute the weights of a mixed factor product grid.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  As a simple test, sum these weights.'
  write ( *, '(a,g24.16)' ) '  They should sum to exactly ', weight_sum_exact
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Spatial dimension DIM_NUM = ', dim_num
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' Dimension      Rule     Order        Alpha          Beta'
  write ( *, '(a)' ) ' '

  do dim = 1, dim_num
    if ( rule(dim) == 1 ) then
      write ( *, '(2x,i8,2x,i8,2x,i8)' ) &
        dim, rule(dim), order_1d(dim)
    else if ( rule(dim) == 2 ) then
      write ( *, '(2x,i8,2x,i8,2x,i8)' ) &
        dim, rule(dim), order_1d(dim)
    else if ( rule(dim) == 3 ) then
      write ( *, '(2x,i8,2x,i8,2x,i8)' ) &
        dim, rule(dim), order_1d(dim)
    else if ( rule(dim) == 4 ) then
      write ( *, '(2x,i8,2x,i8,2x,i8)' ) &
        dim, rule(dim), order_1d(dim)
    else if ( rule(dim) == 5 ) then
      write ( *, '(2x,i8,2x,i8,2x,i8)' ) &
        dim, rule(dim), order_1d(dim)
    else if ( rule(dim) == 6 ) then
      write ( *, '(2x,i8,2x,i8,2x,i8,2x,g14.6)' ) &
        dim, rule(dim), order_1d(dim), alpha(dim)
    else if ( rule(dim) == 7 ) then
      write ( *, '(2x,i8,2x,i8,2x,i8)' ) &
        dim, rule(dim), order_1d(dim)
    else if ( rule(dim) == 8 ) then
      write ( *, '(2x,i8,2x,i8,2x,i8,2x,g14.6)' ) &
        dim, rule(dim), order_1d(dim), alpha(dim)
    else if ( rule(dim) == 9 ) then
      write ( *, '(2x,i8,2x,i8,2x,i8,2x,g14.6,2x,g14.6)' ) &
        dim, rule(dim), order_1d(dim), alpha(dim), beta(dim)
    else if ( rule(dim) == 10 ) then
      write ( *, '(2x,i8,2x,i8,2x,i8)' ) &
        dim, rule(dim), order_1d(dim)
    else if ( rule(dim) == 11 ) then
      write ( *, '(2x,i8,2x,i8,2x,i8)' ) &
        dim, rule(dim), order_1d(dim)
    else if ( rule(dim) == 12 ) then
      write ( *, '(2x,i8,2x,i8,2x,i8)' ) &
        dim, rule(dim), order_1d(dim)
    else if ( rule(dim) == 13 ) then
      write ( *, '(2x,i8,2x,i8,2x,i8)' ) &
        dim, rule(dim), order_1d(dim)
    else if ( rule(dim) == 14 ) then
      write ( *, '(2x,i8,2x,i8,2x,i8)' ) &
        dim, rule(dim), order_1d(dim)
    else if ( rule(dim) == 15 ) then
      write ( *, '(2x,i8,2x,i8,2x,i8)' ) &
        dim, rule(dim), order_1d(dim)
    else if ( rule(dim) == 16 ) then
      write ( *, '(2x,i8,2x,i8,2x,i8)' ) &
        dim, rule(dim), order_1d(dim)
    else if ( rule(dim) == 17 ) then
      write ( *, '(2x,i8,2x,i8,2x,i8)' ) &
        dim, rule(dim), order_1d(dim)
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PRODUCT_MIXED_WEIGHT_TEST - Fatal error!'
      write ( *, '(a,i8)' ) '  Unexpected value of RULE = ', rule(dim)
      stop
    end if
  end do
!
!  Compute the weights and points.
!
  call product_mixed_weight ( dim_num, order_1d, order_nd, rule, &
    alpha, beta, weight )
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

