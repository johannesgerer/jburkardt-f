program main

!****************************************************************************80
!
!! MAIN is the main program for SPARSE_GRID_MIXED_SIZE_TABLE.
!
!  Discussion:
!
!    SPARSE_GRID_MIXED_SIZE_TABLE makes point count tables.
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
!  Reference:
!
!    Fabio Nobile, Raul Tempone, Clayton Webster,
!    A Sparse Grid Stochastic Collocation Method for Partial Differential
!    Equations with Random Input Data,
!    SIAM Journal on Numerical Analysis,
!    Volume 46, Number 5, 2008, pages 2309-2345.
!
  implicit none

  real ( kind = 8 ) alpha_1d
  real ( kind = 8 ) beta_1d
  integer ( kind = 4 ) dim_max
  integer ( kind = 4 ) dim_min
  integer ( kind = 4 ) level_max_max
  integer ( kind = 4 ) level_max_min
  real ( kind = 8 ) r8_epsilon
  integer ( kind = 4 ) rule_1d
  real ( kind = 8 ) tol

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPARSE_GRID_MIXED_SIZE_TABLE'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) &
    '  Tabulate sparse grid size for various rules, levels, dimensions.'
!
!  "Rule 0"
!  A special input that indicates we want a table of the number of polynomials
!  of degree DEGREE or less in a space of dimension DIM.
!
  rule_1d = 0
  alpha_1d = 0.0D+00
  beta_1d = 0.0D+00
  dim_min = 1
  dim_max = 5
  level_max_min = 0
  level_max_max = 7

  call sparse_grid_mixed_size_tabulate ( rule_1d, alpha_1d, beta_1d, dim_min, &
    dim_max, level_max_min, level_max_max )

  write ( *, '(a)' ) ' '
  call timestamp ( )

  rule_1d = 0
  alpha_1d = 0.0D+00
  beta_1d = 0.0D+00
  dim_min = 6
  dim_max = 10
  level_max_min = 0
  level_max_max = 5

  call sparse_grid_mixed_size_tabulate ( rule_1d, alpha_1d, beta_1d, dim_min, &
    dim_max, level_max_min, level_max_max )

  write ( *, '(a)' ) ' '
  call timestamp ( )
!
!  Rule 1
!  Clenshaw Curtis Exponential Growth
!
  rule_1d = 1
  alpha_1d = 0.0D+00
  beta_1d = 0.0D+00
  dim_min = 1
  dim_max = 5
  level_max_min = 0
  level_max_max = 7

  call sparse_grid_mixed_size_tabulate ( rule_1d, alpha_1d, beta_1d, dim_min, &
    dim_max, level_max_min, level_max_max )

  write ( *, '(a)' ) ' '
  call timestamp ( )

  rule_1d = 1
  alpha_1d = 0.0D+00
  beta_1d = 0.0D+00
  dim_min = 6
  dim_max = 10
  level_max_min = 0
  level_max_max = 5

  call sparse_grid_mixed_size_tabulate ( rule_1d, alpha_1d, beta_1d, dim_min, &
    dim_max, level_max_min, level_max_max )

  write ( *, '(a)' ) ' '
  call timestamp ( )

  rule_1d = 1
  alpha_1d = 0.0D+00
  beta_1d = 0.0D+00
  dim_min = 100
  dim_max = 100
  level_max_min = 0
  level_max_max = 2

  call sparse_grid_mixed_size_tabulate ( rule_1d, alpha_1d, beta_1d, dim_min, &
    dim_max, level_max_min, level_max_max )

  write ( *, '(a)' ) ' '
  call timestamp ( )
!
!  Rule 3
!  Gauss Patterson Exponential Growth
!
  rule_1d = 3
  alpha_1d = 0.0D+00
  beta_1d = 0.0D+00
  dim_min = 1
  dim_max = 5
  level_max_min = 0
  level_max_max = 7

  call sparse_grid_mixed_size_tabulate ( rule_1d, alpha_1d, beta_1d, dim_min, &
    dim_max, level_max_min, level_max_max )

  write ( *, '(a)' ) ' '
  call timestamp ( )

  rule_1d = 16
  alpha_1d = 0.0D+00
  beta_1d = 0.0D+00
  dim_min = 6
  dim_max = 10
  level_max_min = 0
  level_max_max = 5

  call sparse_grid_mixed_size_tabulate ( rule_1d, alpha_1d, beta_1d, dim_min, &
    dim_max, level_max_min, level_max_max )

  write ( *, '(a)' ) ' '
  call timestamp ( )
!
!  Rule 11
!  Clenshaw-Curtis Slow Exponential Growth
!
  rule_1d = 11
  alpha_1d = 0.0D+00
  beta_1d = 0.0D+00
  dim_min = 1
  dim_max = 5
  level_max_min = 0
  level_max_max = 7

  call sparse_grid_mixed_size_tabulate ( rule_1d, alpha_1d, beta_1d, dim_min, &
    dim_max, level_max_min, level_max_max )

  write ( *, '(a)' ) ' '
  call timestamp ( )

  rule_1d = 11
  alpha_1d = 0.0D+00
  beta_1d = 0.0D+00
  dim_min = 6
  dim_max = 10
  level_max_min = 0
  level_max_max = 5
  call sparse_grid_mixed_size_tabulate ( rule_1d, alpha_1d, beta_1d, dim_min, &
    dim_max, level_max_min, level_max_max )

  write ( *, '(a)' ) ' '
  call timestamp ( )
!
!  Rule 13
!  Gauss Patterson slow exponential growth
!
  rule_1d = 13
  alpha_1d = 0.0D+00
  beta_1d = 0.0D+00
  dim_min = 1
  dim_max = 5
  level_max_min = 0
  level_max_max = 7

  call sparse_grid_mixed_size_tabulate ( rule_1d, alpha_1d, beta_1d, dim_min, &
    dim_max, level_max_min, level_max_max )

  write ( *, '(a)' ) ' '
  call timestamp ( )

  rule_1d = 13
  alpha_1d = 0.0D+00
  beta_1d = 0.0D+00
  dim_min = 6
  dim_max = 10
  level_max_min = 0
  level_max_max = 5
  call sparse_grid_mixed_size_tabulate ( rule_1d, alpha_1d, beta_1d, dim_min, &
    dim_max, level_max_min, level_max_max )

  write ( *, '(a)' ) ' '
  call timestamp ( )
!
!  Rule 16
!  Gauss Patterson Moderate Exponential Growth
!
  rule_1d = 16
  alpha_1d = 0.0D+00
  beta_1d = 0.0D+00
  dim_min = 1
  dim_max = 5
  level_max_min = 0
  level_max_max = 7

  call sparse_grid_mixed_size_tabulate ( rule_1d, alpha_1d, beta_1d, dim_min, &
    dim_max, level_max_min, level_max_max )

  write ( *, '(a)' ) ' '
  call timestamp ( )

  rule_1d = 16
  alpha_1d = 0.0D+00
  beta_1d = 0.0D+00
  dim_min = 6
  dim_max = 10
  level_max_min = 0
  level_max_max = 5

  call sparse_grid_mixed_size_tabulate ( rule_1d, alpha_1d, beta_1d, dim_min, &
    dim_max, level_max_min, level_max_max )

  write ( *, '(a)' ) ' '
  call timestamp ( )
!
!  Rule 17
!
  rule_1d = 17
  alpha_1d = 0.0D+00
  beta_1d = 0.0D+00
  dim_min = 1
  dim_max = 5
  level_max_min = 0
  level_max_max = 7

  call sparse_grid_mixed_size_tabulate ( rule_1d, alpha_1d, beta_1d, dim_min, &
    dim_max, level_max_min, level_max_max )

  write ( *, '(a)' ) ' '
  call timestamp ( )

  rule_1d = 17
  alpha_1d = 0.0D+00
  beta_1d = 0.0D+00
  dim_min = 6
  dim_max = 10
  level_max_min = 0
  level_max_max = 5

  call sparse_grid_mixed_size_tabulate ( rule_1d, alpha_1d, beta_1d, dim_min, &
    dim_max, level_max_min, level_max_max )

  write ( *, '(a)' ) ' '
  call timestamp ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPARSE_GRID_MIXED_SIZE_TABLE'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  return
end
subroutine sparse_grid_mixed_size_tabulate ( rule_1d, alpha_1d, beta_1d, &
   dim_min, dim_max, level_max_min, level_max_max )

!****************************************************************************80
!
!! SPARSE_GRID_MIXED_SIZE_TABULATE tests SPARSE_GRID_MIXED_SIZE.
!
!  Discussion:
!
!    We do NOT consider mixed rules.  Instead, we are looking at sparse grid
!    rules for which all dimensions use the same 1D rule family.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 February 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) RULE_1D, the 1D rule.
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
!    Input, real ( kind = 8 ) ALPHA_1D, BETA_1D, the optional parameters.
!
!    Input, integer ( kind = 4 ) DIM_MIN, the minimum spatial dimension.
!
!    Input, integer ( kind = 4 ) DIM_MAX, the maximum spatial dimension.
!
!    Input, integer ( kind = 4 ) LEVEL_MAX_MIN, the minimum value of LEVEL_MAX.
!
!    Input, integer ( kind = 4 ) LEVEL_MAX_MAX, the maximum value of LEVEL_MAX.
!
  implicit none

  integer ( kind = 4 ) dim_max
  integer ( kind = 4 ) dim_min

  real ( kind = 8 ), allocatable :: alpha(:)
  real ( kind = 8 ) alpha_1d
  real ( kind = 8 ), allocatable :: beta(:)
  real ( kind = 8 ) beta_1d
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) i4_choose
  integer ( kind = 4 ) level_max
  integer ( kind = 4 ) level_max_max
  integer ( kind = 4 ) level_max_min
  integer ( kind = 4 ) point_num
  integer ( kind = 4 ), allocatable :: point_vec(:)
  real ( kind = 8 ) r8_epsilon
  integer ( kind = 4 ) rule_1d
  integer ( kind = 4 ), allocatable :: rule(:)
  real ( kind = 8 ) tol

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPARSE_GRID_MIXED_SIZE_TABLE'
  write ( *, '(a)' ) '  SPARSE_GRID_MIXED_SIZE returns the number of distinct'
  write ( *, '(a)' ) '  points in a sparse grid.'

  if ( rule_1d == 0 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Here we report the number of polynomials of'
    write ( *, '(a)' ) '  degree DEGREE or less in the given DIM dimensional space.'
    write ( *, '(a)' ) ' '

    allocate ( point_vec(dim_min:dim_max) )

    do dim_num = dim_min, dim_max
      point_vec(dim_num) = dim_num
    end do

    write ( *, '(a8,6(2x,i8))' ) '   DIM: ', point_vec(dim_min:dim_max)
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '   DEGREE'
    write ( *, '(a)' ) ' '

    do level_max = level_max_min, level_max_max

      do dim_num = dim_min, dim_max

        point_vec(dim_num) = i4_choose ( dim_num + level_max, dim_num )

      end do

      write ( *, '(a4,i4,6(2x,i8))' ) '    ', level_max, point_vec(dim_min:dim_max)

    end do

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  We use the same rule in all dimensions, and count the points'
    write ( *, '(a)' ) '  for a range of dimensions and levels.'
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  1D rule index is ', rule_1d
    write ( *, '(a,g14.6)' ) '  ALPHA parameter is ', alpha_1d
    write ( *, '(a,g14.6)' ) '  BETA parameter is  ', beta_1d
    write ( *, '(a)' ) ' '

    tol = sqrt ( r8_epsilon ( ) )

    allocate ( point_vec(dim_min:dim_max) )

    do dim_num = dim_min, dim_max
      point_vec(dim_num) = dim_num
    end do

    write ( *, '(a8,6(2x,i8))' ) '   DIM: ', point_vec(dim_min:dim_max)
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '   LEVEL_MAX'
    write ( *, '(a)' ) ' '

    do level_max = level_max_min, level_max_max

      do dim_num = dim_min, dim_max

        allocate ( rule(1:dim_num) )
        allocate ( alpha(1:dim_num) )
        allocate ( beta(1:dim_num) )

        rule(1:dim_num) = rule_1d
        alpha(1:dim_num) = alpha_1d
        beta(1:dim_num) = beta_1d

        call sparse_grid_mixed_size ( dim_num, level_max, rule, alpha, beta, &
          tol, point_num )

        point_vec(dim_num) = point_num

        deallocate ( alpha )
        deallocate ( beta )
        deallocate ( rule )

      end do

      write ( *, '(a4,i4,6(2x,i8))' ) '    ', level_max, point_vec(dim_min:dim_max)

    end do

  end if

  deallocate ( point_vec )

  return
end
