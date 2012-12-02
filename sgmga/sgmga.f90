subroutine sgmga_aniso_balance ( alpha_max, dim_num, level_weight, &
  level_weight2 )

!*****************************************************************************80
!
!! SGMGA_ANISO_BALANCE "balances" an anisotropic weight vector.
!
!  Discussion:
!
!    The entries in LEVEL_WEIGHT are essentially arbitrary nonnegative numbers.
!
!    The ratio between two entries indicates their relative importance.
!    For instance,
!
!      LEVEL_WEIGHT(1) / LEVEL_WEIGHT(2) = 10
!
!    means that variable 2 is 10 times more important than variable 1.
!    Here, being 10 times more important means that we will generate 10 levels
!    of sparse grid in direction 2 as we generate 1 level in direction 1.
!
!    Under this interpretation, a ratio of 10 already indicates an extreme 
!    imbalanace in the variables, since 10 sparse grid levels in the second
!    variable corresponds roughly to approximating x^1 only, and 
!    all of y^1 through y^10.  A ratio higher than this seems unreasonable.
!
!    Therefore, this function tries to take a somewhat arbitrary level weight
!    vector, and produce a "balanced" level weight vector with the properties
!    that the mininum entry is 1 (representing the item of most importance)
!    and the maximum entry is ALPHA_MAX.  A reasonable value of ALPHA_MAX
!    might be 10 or even 5.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 February 2010
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
!    Fabio Nobile, Raul Tempone, Clayton Webster,
!    An Anisotropic Sparse Grid Stochastic Collocation Method for Partial 
!    Differential Equations with Random Input Data,
!    SIAM Journal on Numerical Analysis,
!    Volume 46, Number 5, 2008, pages 2411-2442.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ALPHA_MAX, the maximum legal value of 
!    LEVEL_WEIGHT, after all entries have been divided by the minimum 
!    nonzero entry.  1 <= ALPHA_MAX.
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, real ( kind = 8 ) LEVEL_WEIGHT(DIM_NUM), the anisotropic weights.  
!    The values must be positive.  
!
!    Output, real ( kind = 8 ) LEVEL_WEIGHT2(DIM_NUM), the balanced 
!    anisotropic weights.  The smallest nonzero entry is 1.0 and 
!    no entry is greater than ALPHA_MAX.
!
  implicit none

  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) alpha_max
  integer ( kind = 4 ) dim
  real ( kind = 8 ) level_weight(dim_num)
  real ( kind = 8 ) level_weight2(dim_num)
  real ( kind = 8 ) level_weight_min
  integer ( kind = 4 ) nonzero_num
  real ( kind = 8 ) r8_huge

  if ( alpha_max < 1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SGMGA_ANISO_BALANCE - Fatal error!'
    write ( *, '(a)' ) '  ALPHA_MAX < 1.0'
    stop
  end if
!
!  Find the smallest nonzero entry.
!
  level_weight_min = r8_huge ( );
  nonzero_num = 0

  do dim = 1, dim_num

    if ( 0.0D+00 < level_weight(dim) ) then
      if ( level_weight(dim) < level_weight_min ) then
        level_weight_min = level_weight(dim)
        nonzero_num = nonzero_num + 1
      end if
    end if

  end do

  if ( nonzero_num == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SGMGA_ANISO_BALANCE - Fatal error!'
    write ( *, '(a)' ) '  Could not find a positive entry in LEVEL_WEIGHT.'
    stop
  end if
!
!  Rescale so the smallest nonzero entry is 1.
!
  level_weight2(1:dim_num) = level_weight(1:dim_num) / level_weight_min
!
!  Set the maximum entry to no more than ALPHA_MAX.
!
  level_weight2(1:dim_num) = min ( alpha_max, level_weight2(1:dim_num) )

  return
end
subroutine sgmga_aniso_normalize ( option, dim_num, level_weight )

!*****************************************************************************80
!
!! SGMGA_ANISO_NORMALIZE normalizes the SGMGA anisotropic weight vector.
!
!  Discussion:
!
!    It is convenient for the user to initialize the anisotropic weight
!    vector with any set of positive values.  These values are to be used
!    as coefficients of the 1D levels, to evaluate an expression which 
!    determines which 1D levels will be included in a given rule.
!
!    This means that a relatively LARGE coefficient forces the corresponding 
!    level to be relatively SMALL.  This is perhaps the opposite of what
!    a user might expect.  If a user wishes to use an importance vector,
!    so that a relatively large importance should correspond to more levels,
!    and hence more points, in that dimension, then the function
!    SGMGA_IMPORTANCE_TO_ANISO should be called first.
!
!    Since the weights only represent the relative importance of the
!    components, they may be multiplied by any (positive) scale factor.
!    Nonetheless, it may be convenient to choose a particular normalization
!    for the weights.  
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 November 2009
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
!    Fabio Nobile, Raul Tempone, Clayton Webster,
!    An Anisotropic Sparse Grid Stochastic Collocation Method for Partial 
!    Differential Equations with Random Input Data,
!    SIAM Journal on Numerical Analysis,
!    Volume 46, Number 5, 2008, pages 2411-2442.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) OPTION, the normalization option.
!    0, no normalization is done at all.
!    1, the minimum nonzero entry will be 1.
!    2, the entries will sum to DIM_NUM.
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input/output, real ( kind = 8 ) LEVEL_WEIGHT(DIM_NUM), the anisotropic
!    weights.  The input values must be positive.  
!    On output, the weights have been normalized.
!
  implicit none

  integer ( kind = 4 ) dim_num

  integer ( kind = 4 ) dim
  real ( kind = 8 ) level_weight(dim_num)
  real ( kind = 8 ) level_weight_min
  real ( kind = 8 ) level_weight_sum
  integer ( kind = 4 ) nonzero_num
  integer ( kind = 4 ) option
  real ( kind = 8 ) r8_huge
!
!  Option 0, no normalization.
!
  if ( option == 0 ) then
!
!  Option 1, minimum (nonzero) entry is 1.
!
  else if ( option == 1 ) then

    level_weight_min = r8_huge ( );
    nonzero_num = 0

    do dim = 1, dim_num

      if ( 0.0D+00 < level_weight(dim) ) then
        if ( level_weight(dim) < level_weight_min ) then
          level_weight_min = level_weight(dim)
          nonzero_num = nonzero_num + 1
        end if
      end if

    end do

    if ( nonzero_num == 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SGMGA_ANISO_NORMALIZE - Fatal error!'
      write ( *, '(a)' ) '  Could not find a positive entry in LEVEL_WEIGHT.'
      stop
    end if

   level_weight(1:dim_num) = level_weight(1:dim_num) / level_weight_min
!
!  Option 2, rescale so sum of weights is DIM_NUM.
!
  else if ( option == 2 ) then

    level_weight_sum = sum ( level_weight(1:dim_num) )

    if ( level_weight_sum <= 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SGMGA_ANISO_NORMALIZE - Fatal error!'
      write ( *, '(a)' ) '  LEVEL_WEIGHT entries have nonpositive sum.'
      stop
    end if

    level_weight(1:dim_num) = ( real ( dim_num, kind = 8 ) &
      * level_weight(1:dim_num) ) / level_weight_sum

  end if

  return
end
subroutine sgmga_importance_to_aniso ( dim_num, importance, level_weight )

!*****************************************************************************80
!
!! SGMGA_IMPORTANCE_TO_ANISO: importance vector to anisotropic weight vector.
!
!  Discussion:
!
!    To specify the anisotropy of a multidimensional problem, the user is
!    allowed to specify an "importance vector".  This vector can contain
!    any set of nonnegative values.  These values represent the relative
!    importance of each dimension.  These values, with a suitable normalization,
!    will be used to evaluate a constraint of the following form:
!
!      QMIN < Level(1) / Importance(1) + Level(2) / Importance(2) + ...
!             Level(N) / Importance(N) <= QMAX
!
!    and a set of levels that satisfies this constraint will then be included
!    in a given anistotropic sparse grid rule.  Thus, increasing the
!    importance value of a particular dimension allows larger level values
!    in that dimension to satisfy the constraint.
!
!    The program actually works with coefficients LEVEL_WEIGHT that are
!    the inverse of the importance vector entries, with a suitable
!    normalization.  This function is supplied to convert between the
!    more natural "importance vector" and the internally useful 
!    "level_weight" vector.
!
!    This function converts the importance vector to an unnormalized 
!    anisotropy weight vector. 
!
!    Note that some (but not all) of the IMPORTANCE vector entries may be zero.
!    This indicates that the corresponding dimension is of "zero" or
!    rather "minimal" importance.  In such a case, only a one-point quadrature
!    rule will be applied for that dimension, no matter what sparse grid
!    level is requested for the overall problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 November 2009
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
!    Fabio Nobile, Raul Tempone, Clayton Webster,
!    An Anisotropic Sparse Grid Stochastic Collocation Method for Partial 
!    Differential Equations with Random Input Data,
!    SIAM Journal on Numerical Analysis,
!    Volume 46, Number 5, 2008, pages 2411-2442.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, real ( kind = 8 ) IMPORTANCE(DIM_NUM), the importance vector.
!    All entries must be nonnegative, and at least one entry must be positive.
!
!    Output, real ( kind = 8 ) LEVEL_WEIGHT(DIM_NUM), the anisotropic
!    weights.
!
  implicit none

  integer ( kind = 4 ) dim_num

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) found
  real ( kind = 8 ) importance(dim_num)
  real ( kind = 8 ) level_weight(dim_num)

  if ( any ( importance(1:dim_num) < 0.0D+00 ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SGMGA_IMPORTANCE_TO_ANISO - Fatal error!'
    write ( *, '(a)' ) '  Some importance entries are negative.'
    stop
  end if

  found = 0

  do dim = 1, dim_num
    if ( 0.0D+00 < importance(dim) ) then
      level_weight(dim) = 1.0D+00 / importance(dim)
      found = found + 1
    else
      level_weight(dim) = 0.0D+00
    end if
  end do

  if ( found == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SGMGA_IMPORTANCE_TO_ANISO - Fatal error!'
    write ( *, '(a)' ) '  No importance entry is positive.'
    stop
  end if

  return
end
subroutine sgmga_index ( dim_num, level_weight, level_max, rule, growth, &
  point_num, point_total_num, sparse_unique_index, sparse_order, sparse_index )

!*****************************************************************************80
!
!! SGMGA_INDEX indexes the unique points in an SGMGA grid.
!
!  Discussion:
!
!    For each "unique" point in the sparse grid, we return its INDEX and ORDER.
!
!    That is, for the I-th unique point P, we determine the product grid which
!    first generated this point, and we return in SPARSE_ORDER the orders of 
!    the 1D rules in that grid, and in SPARSE_INDEX the component indexes in 
!    those rules that generated this specific point.
!
!    For instance, say P was first generated by a rule which was a 3D product
!    of a 9th order CC rule and a 15th order GL rule, and that to generate P,
!    we used the 7-th point of the CC rule and the 3rd point of the GL rule.
!    Then the SPARSE_ORDER information would be (9,15) and the SPARSE_INDEX
!    information would be (7,3).  This, combined with the information in RULE,
!    is enough to regenerate the value of P.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 April 2011
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
!    Fabio Nobile, Raul Tempone, Clayton Webster,
!    An Anisotropic Sparse Grid Stochastic Collocation Method for Partial 
!    Differential Equations with Random Input Data,
!    SIAM Journal on Numerical Analysis,
!    Volume 46, Number 5, 2008, pages 2411-2442.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, real ( kind = 8 ) LEVEL_WEIGHT(DIM_NUM), the anisotropic weights.
!
!    Input, integer ( kind = 4 ) LEVEL_MAX, the maximum value of LEVEL.
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
!    Input, integer ( kind = 4 ) GROWTH(DIM_NUM), the growth in each dimension.
!    0, "DF", default growth associated with this quadrature rule;
!    1, "SL", slow linear, L+1;
!    2  "SO", slow linear odd, O=1+2((L+1)/2)
!    3, "ML", moderate linear, 2L+1;
!    4, "SE", slow exponential;
!    5, "ME", moderate exponential;
!    6, "FE", full exponential.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of unique points 
!    in the grid. 
!
!    Input, integer ( kind = 4 ) POINT_TOTAL_NUM, the total number of points 
!    in the grid.
!
!    Input, integer ( kind = 4 ) SPARSE_UNIQUE_INDEX(POINT_TOTAL_NUM), 
!    associates each point in the grid with its unique representative.
!
!    Output, integer ( kind = 4 ) SPARSE_ORDER(DIM_NUM,POINT_NUM), lists, 
!    for each point, the order of the 1D rules used in the grid that 
!    generated it.
!
!    Output, integer ( kind = 4 ) SPARSE_INDEX(DIM_NUM,POINT_NUM), lists, for 
!    each point, its index in each of the 1D rules in the grid that generated 
!    it.  The indices are 1-based.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num
  integer ( kind = 4 ) point_total_num

  real ( kind = 8 ) coef
  integer ( kind = 4 ) dim
  integer ( kind = 4 ) growth(dim_num)
  integer ( kind = 4 ) level_1d(dim_num)
  integer ( kind = 4 ) level_1d_max(dim_num)
  integer ( kind = 4 ) level_max
  real ( kind = 8 ) level_weight(dim_num)
  real ( kind = 8 ) level_weight_min_pos
  logical more_grids
  logical more_points
  integer ( kind = 4 ) order_1d(dim_num)
  integer ( kind = 4 ) point_count
  integer ( kind = 4 ) point_index(dim_num)
  integer ( kind = 4 ) point_unique
  real ( kind = 8 ) q_max
  real ( kind = 8 ) q_min
  integer ( kind = 4 ) r8_floor
  real ( kind = 8 ) r8vec_min_pos
  integer ( kind = 4 ) rule(dim_num)
  integer ( kind = 4 ) sparse_index(dim_num,point_num)
  integer ( kind = 4 ) sparse_order(dim_num,point_num)
  integer ( kind = 4 ) sparse_unique_index(point_total_num)
!
!  Special cases.
!
  if ( level_max < 0 ) then
    return
  end if

  if ( level_max == 0 ) then
    sparse_order(1:dim_num,1) = 1
    sparse_index(1:dim_num,1) = 1
    return
  end if
!
!  Initialize the INDEX and ORDER arrays to -1 to help catch errors.
!
  sparse_order(1:dim_num,1:point_num) = -1
  sparse_index(1:dim_num,1:point_num) = -1

  point_count = 0
!
!  Initialization.
!
  level_weight_min_pos = r8vec_min_pos ( dim_num, level_weight )
  q_min = real ( level_max, kind = 8 ) * level_weight_min_pos &
    - sum ( level_weight(1:dim_num) )
  q_max = real ( level_max, kind = 8 ) * level_weight_min_pos
  do dim = 1, dim_num
    if ( 0.0D+00 < level_weight(dim) ) then
      level_1d_max(dim) = r8_floor ( q_max / level_weight(dim) ) + 1
      if ( q_max <= ( level_1d_max(dim) - 1 ) * level_weight(dim) ) then
        level_1d_max(dim) = level_1d_max(dim) - 1
      end if
    else
      level_1d_max(dim) = 0
    end if
  end do
  more_grids = .false.
!
!  Seek all vectors LEVEL_1D which satisfy the constraint:
!
!    LEVEL_MAX * LEVEL_WEIGHT_MIN_POS - sum ( LEVEL_WEIGHT ) 
!      < sum ( 1 <= I <= DIM_NUM ) LEVEL_WEIGHT(I) * LEVEL_1D(I)
!      <= LEVEL_MAX * LEVEL_WEIGHT_MIN_POS.
!
  do

    call sgmga_vcn_ordered ( dim_num, level_weight, level_1d_max, &
      level_1d, q_min, q_max, more_grids )

    if ( .not. more_grids ) then
      exit
    end if
!
!  Compute the combinatorial coefficient.
!
    call sgmga_vcn_coef ( dim_num, level_weight, level_1d, q_max, coef )

    if ( coef == 0.0D+00 ) then
      cycle
    end if
!
!  Transform each 1D level to a corresponding 1D order.
!
    call level_growth_to_order ( dim_num, level_1d, rule, growth, order_1d )
!
!  The inner loop generates a POINT of the GRID of the LEVEL.
!
    more_points = .false.

    do

      call vec_colex_next3 ( dim_num, order_1d, point_index, more_points )

      if ( .not. more_points ) then
        exit
      end if

      point_count = point_count + 1
      point_unique = sparse_unique_index(point_count)
      sparse_order(1:dim_num,point_unique) = order_1d(1:dim_num)
      sparse_index(1:dim_num,point_unique) = point_index(1:dim_num)

    end do

  end do

  return
end
subroutine sgmga_point ( dim_num, level_weight, level_max, rule, growth, np, &
  p, point_num, sparse_order, sparse_index, sparse_point )

!*****************************************************************************80
!
!! SGMGA_POINT computes the points of an SGMGA rule.
!
!  Discussion:
!
!    The sparse grid is the logical sum of low degree product rules.
!
!    Each product rule is the product of 1D factor rules.
!
!    The user specifies:
!    * the spatial dimension of the quadrature region,
!    * the weighting of each level.
!    * the level that defines the Smolyak grid.
!    * the quadrature rules.
!    * the number of points.
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
!  Reference:
!
!    Fabio Nobile, Raul Tempone, Clayton Webster,
!    A Sparse Grid Stochastic Collocation Method for Partial Differential
!    Equations with Random Input Data,
!    SIAM Journal on Numerical Analysis,
!    Volume 46, Number 5, 2008, pages 2309-2345.
!
!    Fabio Nobile, Raul Tempone, Clayton Webster,
!    An Anisotropic Sparse Grid Stochastic Collocation Method for Partial 
!    Differential Equations with Random Input Data,
!    SIAM Journal on Numerical Analysis,
!    Volume 46, Number 5, 2008, pages 2411-2442.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, real ( kind = 8 ) LEVEL_WEIGHT(DIM_NUM), the anisotropic weights.
!
!    Input, integer ( kind = 4 ) LEVEL_MAX, controls the size of the final
!    sparse grid.
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
!    Input, integer ( kind = 4 ) GROWTH(DIM_NUM), the growth in each dimension.
!    0, "DF", default growth associated with this quadrature rule;
!    1, "SL", slow linear, L+1;
!    2  "SO", slow linear odd, O=1+2((L+1)/2)
!    3, "ML", moderate linear, 2L+1;
!    4, "SE", slow exponential;
!    5, "ME", moderate exponential;
!    6, "FE", full exponential.
!
!    Input, integer ( kind = 4 ) NP(DIM_NUM), the number of parameters 
!    used by each rule.
!
!    Input, real ( kind = 8 ) P(*), the parameters needed by each rule.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of unique points 
!    in the grid. 
!
!    Input, integer ( kind = 4 ) SPARSE_ORDER(DIM_NUM,POINT_NUM), lists, for 
!    each point, the order of the 1D rules used in the grid that generated it.
!
!    Input, integer ( kind = 4 ) SPARSE_INDEX(DIM_NUM,POINT_NUM), lists, for 
!    each point, its index in each of the 1D rules in the grid that generated 
!    it.  The indices are 1-based.
!
!    Output, real ( kind = 8 ) SPARSE_POINT(DIM_NUM,POINT_NUM), the points.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) dim2
  integer ( kind = 4 ) growth(dim_num)
  integer ( kind = 4 ) level
  integer ( kind = 4 ) level_1d_max(dim_num)
  integer ( kind = 4 ) level_max
  real ( kind = 8 ) level_weight(dim_num)
  real ( kind = 8 ) level_weight_min_pos
  integer ( kind = 4 ) np(dim_num)
  integer ( kind = 4 ) order
  real ( kind = 8 ) p(*)
  integer ( kind = 4 ) p_index
  integer ( kind = 4 ) point
  real ( kind = 8 ), allocatable, dimension ( : ) :: points
  real ( kind = 8 ) q_max
  integer ( kind = 4 ) r8_floor
  real ( kind = 8 ) r8_huge
  real ( kind = 8 ) r8vec_min_pos
  integer ( kind = 4 ) rule(dim_num)
  integer ( kind = 4 ) sparse_index(dim_num,point_num)
  integer ( kind = 4 ) sparse_order(dim_num,point_num)
  real ( kind = 8 ) sparse_point(dim_num,point_num)
!
!  Compute the point coordinates.
!
  sparse_point(1:dim_num,1:point_num) = - r8_huge ( )
  
  level_weight_min_pos = r8vec_min_pos ( dim_num, level_weight )
  q_max = real ( level_max, kind = 8 ) * level_weight_min_pos

  p_index = 1

  do dim = 1, dim_num

    if ( 0.0D+00 < level_weight(dim) ) then
      level_1d_max(dim) = r8_floor ( q_max / level_weight(dim) ) + 1
      if ( q_max <= ( level_1d_max(dim) - 1 ) * level_weight(dim) ) then
        level_1d_max(dim) = level_1d_max(dim) - 1
      end if
    else
      level_1d_max(dim) = 0
    end if

    do level = 0, level_1d_max(dim)

      call level_growth_to_order ( 1, level, rule(dim), growth(dim), order )

      allocate ( points(1:order) )

      if ( rule(dim) == 1 ) then
        call clenshaw_curtis_compute_points_np ( order, np(dim), p(p_index), &
          points )
      else if ( rule(dim) == 2 ) then
        call fejer2_compute_points_np ( order, np(dim), p(p_index), points )
      else if ( rule(dim) == 3 ) then
        call patterson_lookup_points_np ( order, np(dim), p(p_index), points )
      else if ( rule(dim) == 4 ) then
        call legendre_compute_points_np ( order, np(dim), p(p_index), points )
      else if ( rule(dim) == 5 ) then
        call hermite_compute_points_np ( order, np(dim), p(p_index), points )
      else if ( rule(dim) == 6 ) then
        call gen_hermite_compute_points_np ( order, np(dim), p(p_index), &
          points )
      else if ( rule(dim) == 7 ) then
        call laguerre_compute_points_np ( order, np(dim), p(p_index), points )
      else if ( rule(dim) == 8 ) then
        call gen_laguerre_compute_points_np ( order, np(dim), p(p_index), &
          points )
      else if ( rule(dim) == 9 ) then
        call jacobi_compute_points_np ( order, np(dim), p(p_index), points )
      else if ( rule(dim) == 10 ) then
        call hermite_genz_keister_lookup_points_np ( order, np(dim), &
          p(p_index), points )
      else if ( rule(dim) == 11 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SGMGA_POINT - Fatal error!'
        write ( *, '(a)' ) '  Do not know how to assign points for rule 11.'
        stop
      else if ( rule(dim) == 12 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SGMGA_POINT - Fatal error!'
        write ( *, '(a)' ) '  Do not know how to assign points for rule 12.'
        stop
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SGMGA_POINT - Fatal error!'
        write ( *, '(a,i8)' ) '  Unexpected value of RULE = ', rule(dim)
        stop
      end if

      do point = 1, point_num
        if ( sparse_order(dim,point) == order ) then
          sparse_point(dim,point) = points ( sparse_index(dim,point) )
        end if
      end do

      deallocate ( points )

    end do

    p_index = p_index + np(dim)

  end do
!
!  Check to see if we missed any points.
!
  do point = 1, point_num
    do dim = 1, dim_num
      if ( sparse_point(dim,point) == - r8_huge ( ) ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SGMGA_POINT - Fatal error!'
        write ( *, '(a)' ) '  At least one point component was not assigned.'
        write ( *, '(a,i8)' ) '  POINT = ', point
        write ( *, '(a,i8)' ) '  DIM = ', dim
        write ( *, '(a,i8)' ) '  SPARSE_ORDER(DIM,POINT) = ', &
          sparse_order(dim,point)
        write ( *, '(a,g14.6)' ) '  LEVEL_WEIGHT(DIM) = ', level_weight(dim)
        write ( *, '(a,i8)' ) '  LEVEL_1D_MAX(DIM) = ', level_1d_max(dim)
        write ( *, '(a,g14.6)' ) '  Q_MAX = ', q_max
        write ( *, '(a,g14.6)' ) &
          '  LEVEL_WEIGHT_MIN_POS = ', level_weight_min_pos
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  SPARSE_POINT(*,POINT):'
        write ( *, '(a)' ) ' '
        do dim2 = 1, dim_num
          write ( *, '(2x,i8,2x,g14.6)' ) dim2, sparse_point(dim2,point)
        end do
        stop
      end if
    end do 
  end do

  return
end
subroutine sgmga_product_weight ( dim_num, order_1d, order_nd, rule, &
  np, p, weight_nd )

!*****************************************************************************80
!
!! SGMGA_PRODUCT_WEIGHT computes the weights of a mixed product rule.
!
!  Discussion:
!
!    This routine computes the weights for a quadrature rule which is
!    a product of 1D rules of varying order and kind.
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
!  Reference:
!
!    Fabio Nobile, Raul Tempone, Clayton Webster,
!    A Sparse Grid Stochastic Collocation Method for Partial Differential
!    Equations with Random Input Data,
!    SIAM Journal on Numerical Analysis,
!    Volume 46, Number 5, 2008, pages 2309-2345.
!
!    Fabio Nobile, Raul Tempone, Clayton Webster,
!    An Anisotropic Sparse Grid Stochastic Collocation Method for Partial 
!    Differential Equations with Random Input Data,
!    SIAM Journal on Numerical Analysis,
!    Volume 46, Number 5, 2008, pages 2411-2442.
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
!    Output, real ( kind = 8 ) WEIGHT_ND(ORDER_ND), the product rule weights.
!
  implicit none

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) order_nd

  integer ( kind = 4 ) np(dim_num)
  integer ( kind = 4 ) order_1d(dim_num)
  real ( kind = 8 ) p(*)
  integer ( kind = 4 ) p_index
  integer ( kind = 4 ) rule(dim_num)
  real ( kind = 8 ), allocatable, dimension ( : ) :: weight_1d
  real ( kind = 8 ) weight_nd(order_nd)

  weight_nd(1:order_nd) = 1.0D+00

  p_index = 1

  do dim = 1, dim_num

    allocate ( weight_1d(1:order_1d(dim) ) )

    if ( rule(dim) == 1 ) then
      call clenshaw_curtis_compute_weights_np ( order_1d(dim), &
        np(dim), p(p_index), weight_1d )
    else if ( rule(dim) == 2 ) then
      call fejer2_compute_weights_np ( order_1d(dim), np(dim), &
        p(p_index), weight_1d )
    else if ( rule(dim) == 3 ) then
      call patterson_lookup_weights_np ( order_1d(dim), np(dim), &
        p(p_index), weight_1d )
    else if ( rule(dim) == 4 ) then
      call legendre_compute_weights_np ( order_1d(dim), np(dim), &
       p(p_index), weight_1d )
    else if ( rule(dim) == 5 ) then
      call hermite_compute_weights_np ( order_1d(dim), np(dim), &
        p(p_index), weight_1d )
    else if ( rule(dim) == 6 ) then
      call gen_hermite_compute_weights_np ( order_1d(dim), np(dim), &
        p(p_index), weight_1d )
    else if ( rule(dim) == 7 ) then
      call laguerre_compute_weights_np ( order_1d(dim), np(dim), p(p_index), &
        weight_1d )
    else if ( rule(dim) == 8 ) then
      call gen_laguerre_compute_weights_np ( order_1d(dim), np(dim), &
        p(p_index), weight_1d )
    else if ( rule(dim) == 9 ) then
      call jacobi_compute_weights_np ( order_1d(dim), np(dim), p(p_index), &
        weight_1d )
    else if ( rule(dim) == 10 ) then
      call hermite_genz_keister_lookup_weights_np ( order_1d(dim), np(dim), &
        p(p_index), weight_1d )
    else if ( rule(dim) == 11 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SGMGA_PRODUCT_WEIGHT - Fatal error!'
      write ( *, '(a,i8)' ) '  Do not know how to set weights for rule 11.'
      stop
    else if ( rule(dim) == 12 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SGMGA_PRODUCT_WEIGHT - Fatal error!'
      write ( *, '(a,i8)' ) '  Do not know how to set weights for rule 12.'
      stop
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SGMGA_PRODUCT_WEIGHT - Fatal error!'
      write ( *, '(a,i8)' ) '  Unexpected value of RULE = ', rule(dim)
      stop
    end if

    call r8vec_direct_product2 ( dim, order_1d(dim), weight_1d, &
      dim_num, order_nd, weight_nd )

    deallocate ( weight_1d )
 
    p_index = p_index + np(dim)

  end do

  return
end
subroutine sgmga_size ( dim_num, level_weight, level_max, rule, growth, &
  np, p, tol, point_num )

!*****************************************************************************80
!
!! SGMGA_SIZE sizes an SGMGA grid, discounting duplicate points.
!
!  Discussion:
!
!    The sparse grid is the logical sum of product grids that satisfy
!    a particular constraint.
!
!    Depending on the 1D rules involved, there may be many duplicate points
!    in the sparse grid.
!
!    This routine counts the unique points in the sparse grid.  It does this
!    in a straightforward way, by actually generating all the points, and
!    comparing them, with a tolerance for equality.
!
!    This function has been modified to automatically omit points for which
!    the "combinatorial coefficient" is zero, since such points would have
!    a weight of zero in the grid.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 April 2011
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
!    Fabio Nobile, Raul Tempone, Clayton Webster,
!    An Anisotropic Sparse Grid Stochastic Collocation Method for Partial 
!    Differential Equations with Random Input Data,
!    SIAM Journal on Numerical Analysis,
!    Volume 46, Number 5, 2008, pages 2411-2442.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, real ( kind = 8 ) LEVEL_WEIGHT(DIM_NUM), the anisotropic weights.
!
!    Input, integer ( kind = 4 ) LEVEL_MAX, the maximum value of LEVEL.
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
!    Input, integer ( kind = 4 ) GROWTH(DIM_NUM), the growth in each dimension.
!    0, "DF", default growth associated with this quadrature rule;
!    1, "SL", slow linear, L+1;
!    2  "SO", slow linear odd, O=1+2((L+1)/2)
!    3, "ML", moderate linear, 2L+1;
!    4, "SE", slow exponential;
!    5, "ME", moderate exponential;
!    6, "FE", full exponential.
!
!    Input, integer ( kind = 4 ) NP(DIM_NUM), the number of parameters 
!    used by each rule.
!
!    Input, real ( kind = 8 ) P(*), the parameters needed by each rule.
!
!    Input, real ( kind = 8 ) TOL, the tolerance for point equality.
!
!    Output, integer ( kind = 4 ) POINT_NUM, the number of unique points 
!    in the grid. 
!
  implicit none

  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) coef
  integer ( kind = 4 ) dim
  integer ( kind = 4 ) growth(dim_num)
  integer ( kind = 4 ) level
  integer ( kind = 4 ) level_1d(dim_num)
  integer ( kind = 4 ) level_1d_max(dim_num)
  integer ( kind = 4 ) level_max
  real ( kind = 8 ) level_weight(dim_num)
  real ( kind = 8 ) level_weight_min_pos
  logical more_grids
  logical more_points
  integer ( kind = 4 ) np(dim_num)
  integer ( kind = 4 ) order
  integer ( kind = 4 ) order_1d(dim_num)
  real ( kind = 8 ) p(*)
  integer ( kind = 4 ) p_index
  integer ( kind = 4 ) point
  integer ( kind = 4 ) point_index(dim_num)
  integer ( kind = 4 ) point_num
  integer ( kind = 4 ) point_total_num
  integer ( kind = 4 ) point_total_num2
  real ( kind = 8 ), allocatable, dimension ( : ) :: points
  real ( kind = 8 ) q_max
  real ( kind = 8 ) q_min
  integer ( kind = 4 ) r8_floor
  real ( kind = 8 ) r8_huge
  real ( kind = 8 ) r8vec_min_pos
  integer ( kind = 4 ) rule(dim_num)
  integer ( kind = 4 ) seed
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: sparse_total_index
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: sparse_total_order
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: sparse_total_point
  real ( kind = 8 ) tol
!
!  Special cases.
!
  if ( level_max < 0 ) then
    point_num = -1
    return
  end if

  if ( level_max == 0 ) then
    point_num = 1
    return
  end if
!
!  Get total number of points, including duplicates.
!
  call sgmga_size_total ( dim_num, level_weight, level_max, rule, &
    growth, point_total_num )
!
!  Generate SPARSE_TOTAL_ORDER and SPARSE_TOTAL_INDEX arrays 
!  for the TOTAL set of points.
!
  allocate ( sparse_total_order(1:dim_num,1:point_total_num ) )
  allocate ( sparse_total_index(1:dim_num,1:point_total_num ) )

  point_total_num2 = 0
!
!  Initialization.
!
  level_weight_min_pos = r8vec_min_pos ( dim_num, level_weight )
  q_min = real ( level_max, kind = 8 ) * level_weight_min_pos &
    - sum ( level_weight(1:dim_num) )
  q_max = real ( level_max, kind = 8 ) * level_weight_min_pos
  do dim = 1, dim_num
    if ( 0.0D+00 < level_weight(dim) ) then
      level_1d_max(dim) = r8_floor ( q_max / level_weight(dim) ) + 1
      if ( q_max <= ( level_1d_max(dim) - 1 ) * level_weight(dim) ) then
        level_1d_max(dim) = level_1d_max(dim) - 1
      end if
    else
      level_1d_max(dim) = 0
    end if
  end do
  more_grids = .false.
!
!  Seek all vectors LEVEL_1D which satisfy the constraint:
!
!    LEVEL_MAX * LEVEL_WEIGHT_MIN_POS - sum ( LEVEL_WEIGHT ) 
!      < sum ( 1 <= I <= DIM_NUM ) LEVEL_WEIGHT(I) * LEVEL_1D(I)
!      <= LEVEL_MAX * LEVEL_WEIGHT_MIN_POS.
!
  do

    call sgmga_vcn_ordered ( dim_num, level_weight, level_1d_max, &
      level_1d, q_min, q_max, more_grids )

    if ( .not. more_grids ) then
      exit
    end if
!
!  Compute the combinatorial coefficient.
!
!   call sgmga_vcn_coef_naive ( dim_num, level_weight, level_1d_max, level_1d, &
!     q_min, q_max, coef )

    call sgmga_vcn_coef ( dim_num, level_weight, level_1d, q_max, coef )

    if ( coef == 0.0D+00 ) then
      cycle
    end if
!
!  Transform each 1D level to a corresponding 1D order.
!
    call level_growth_to_order ( dim_num, level_1d, rule, growth, order_1d )
!
!  The inner loop generates a POINT of the GRID of the LEVEL.
!
    more_points = .false.

    do

      call vec_colex_next3 ( dim_num, order_1d, point_index, more_points )

      if ( .not. more_points ) then
        exit
      end if

      point_total_num2 = point_total_num2 + 1
      sparse_total_order(1:dim_num,point_total_num2) = order_1d(1:dim_num)
      sparse_total_index(1:dim_num,point_total_num2) = point_index(1:dim_num)

    end do

  end do
!
!  Now compute the coordinates of the TOTAL set of points.
!
  allocate ( sparse_total_point(1:dim_num,1:point_total_num) )

  sparse_total_point(1:dim_num,1:point_total_num) = r8_huge ( )

  level_weight_min_pos = r8vec_min_pos ( dim_num, level_weight )
  q_max = real ( level_max, kind = 8 ) * level_weight_min_pos

  p_index = 1

  do dim = 1, dim_num

    if ( 0.0D+00 < level_weight(dim) ) then
      level_1d_max(dim) = r8_floor ( q_max / level_weight(dim) ) + 1
      if ( q_max <= ( level_1d_max(dim) - 1 ) * level_weight(dim) ) then
        level_1d_max(dim) = level_1d_max(dim) - 1
      end if
    else
      level_1d_max(dim) = 0
    end if

    do level = 0, level_1d_max(dim)

      call level_growth_to_order ( 1, level, rule(dim), growth(dim), order )

      allocate ( points(1:order) )

      if ( rule(dim) == 1 ) then
        call clenshaw_curtis_compute_points_np ( order, np(dim), p(p_index), &
          points )
      else if ( rule(dim) == 2 ) then
        call fejer2_compute_points_np ( order, np(dim), p(p_index), points )
      else if ( rule(dim) == 3 ) then
        call patterson_lookup_points_np ( order, np(dim), p(p_index), points )
      else if ( rule(dim) == 4 ) then
         call legendre_compute_points_np ( order, np(dim), p(p_index), points )
      else if ( rule(dim) == 5 ) then
        call hermite_compute_points_np ( order, np(dim), p(p_index), points )
      else if ( rule(dim) == 6 ) then
        call gen_hermite_compute_points_np ( order, np(dim), p(p_index), &
          points )
      else if ( rule(dim) == 7 ) then
        call laguerre_compute_points_np ( order, np(dim), p(p_index), points )
      else if ( rule(dim) == 8 ) then
        call gen_laguerre_compute_points_np ( order, np(dim), p(p_index), &
          points )
      else if ( rule(dim) == 9 ) then
        call jacobi_compute_points_np ( order, np(dim), p(p_index), points )
      else if ( rule(dim) == 10 ) then
        call hermite_genz_keister_lookup_points_np ( order, np(dim), &
          p(p_index), points )
      else if ( rule(dim) == 11 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SGMGA_SIZE - Fatal error!'
        write ( *, '(a)' ) '  Do not know how to assign points for rule 11.'
        stop
      else if ( rule(dim) == 12 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SGMGA_SIZE - Fatal error!'
        write ( *, '(a)' ) '  Do not know how to assign points for rule 12.'
        stop
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SGMGA_SIZE - Fatal error!'
        write ( *, '(a,i8)' ) '  Unexpected value of RULE = ', rule(dim)
        stop
      end if

      do point = 1, point_total_num
        if ( sparse_total_order(dim,point) == order ) then
          sparse_total_point(dim,point) = &
            points ( sparse_total_index(dim,point) )
        end if
      end do

      deallocate ( points )

    end do

    p_index = p_index + np(dim)

  end do
!
!  Count the tolerably unique points. 
!
  seed = 123456789

  call point_radial_tol_unique_count ( dim_num, point_total_num, &
    sparse_total_point, tol, seed, point_num )

  deallocate ( sparse_total_index )
  deallocate ( sparse_total_order )
  deallocate ( sparse_total_point )

  return
end
subroutine sgmga_size_total ( dim_num, level_weight, level_max, rule, growth, &
  point_total_num )

!*****************************************************************************80
!
!! SGMGA_SIZE_TOTAL sizes an SGMGA grid, counting duplicate points.
!
!  Discussion:
!
!    This routine returns the total point count for an SGMGA
!    ( Sparse Grid of Mixed type with Growth rule and Anisotropic weights).
!
!    The sparse grid is the logical sum of product grids.
!
!    The sparse grid has an associated integer index LEVEL_MAX, whose lowest 
!    value is 0.  LEVEL_MAX = 0 indicates the sparse grid made up of one 
!    product grid, which in turn is the product of 1D factor grids of the 
!    lowest level.  This usually means the sparse grid with LEVEL_MAX equal 
!    to 0 is a one point grid.
!
!    We can assign a level to each factor grid, and hence a LEVEL vector
!    to the corresponding product grid, and a weighted index
!    LEVEL_GRID (which will in general be a real number):
!
!      LEVEL_GRID = sum ( 1 <= I <= DIM_NUM ) LEVEL_WEIGHT(I) * LEVEL(I)
!
!    The product grid will participate in the formation of the sparse grid
!    if it satisfies the following weighted constraint:
!
!      LEVEL_MAX - DIM_NUM < LEVEL_GRID <= LEVEL_MAX
!
!    This routine determines the total number of abscissas in all the 
!    product rules used to form the SGMGA associated with the index LEVEL_MAX.
!    The count disregards duplication.  If the same multidimensional abcsissa
!    occurs in two different product rules that are part of the SGMGA, then
!    that single abcissa is counted twice. 
!
!    This computation is useful in cases where the entire set of abscissas
!    is going to be generated, preparatory to compression to finding, indexing
!    and merging the duplicate abcissass.
!
!    This function has been modified to automatically omit points for which
!    the "combinatorial coefficient" is zero, since such points would have
!    a weight of zero in the grid.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 April 2011
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
!    Fabio Nobile, Raul Tempone, Clayton Webster,
!    An Anisotropic Sparse Grid Stochastic Collocation Method for Partial 
!    Differential Equations with Random Input Data,
!    SIAM Journal on Numerical Analysis,
!    Volume 46, Number 5, 2008, pages 2411-2442.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, real ( kind = 8 ) LEVEL_WEIGHT(DIM_NUM), the anisotropic weights.
!
!    Input, integer ( kind = 4 ) LEVEL_MAX, the maximum value of LEVEL.
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
!    Input, integer ( kind = 4 ) GROWTH(DIM_NUM), the growth in each dimension.
!    0, "DF", default growth associated with this quadrature rule;
!    1, "SL", slow linear, L+1;
!    2  "SO", slow linear odd, O=1+2((L+1)/2)
!    3, "ML", moderate linear, 2L+1;
!    4, "SE", slow exponential;
!    5, "ME", moderate exponential;
!    6, "FE", full exponential.
!
!    Output, integer ( kind = 4 ) POINT_TOTAL_NUM, the total number of points 
!    in the grid. 
!
  implicit none

  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) coef
  integer ( kind = 4 ) dim
  integer ( kind = 4 ) growth(dim_num)
  integer ( kind = 4 ) level_1d(dim_num)
  integer ( kind = 4 ) level_1d_max(dim_num)
  integer ( kind = 4 ) level_max
  real ( kind = 8 ) level_weight(dim_num)
  real ( kind = 8 ) level_weight_min_pos
  logical more_grids
  integer ( kind = 4 ) order_1d(dim_num)
  integer ( kind = 4 ) point_total_num
  real ( kind = 8 ) q_max
  real ( kind = 8 ) q_min
  integer ( kind = 4 ) r8_floor
  real ( kind = 8 ) r8vec_min_pos
  integer ( kind = 4 ) rule(dim_num)
!
!  Special case.
!
  if ( level_max == 0 ) then
    point_total_num = 1
    return
  end if

  point_total_num = 0
!
!  Initialization.
!
  level_weight_min_pos = r8vec_min_pos ( dim_num, level_weight )
  q_min = real ( level_max, kind = 8 ) * level_weight_min_pos &
    - sum ( level_weight(1:dim_num) )
  q_max = real ( level_max, kind = 8 ) * level_weight_min_pos
  do dim = 1, dim_num
    if ( 0.0D+00 < level_weight(dim) ) then
      level_1d_max(dim) = r8_floor ( q_max / level_weight(dim) ) + 1
      if ( q_max <= ( level_1d_max(dim) - 1 ) * level_weight(dim) ) then
        level_1d_max(dim) = level_1d_max(dim) - 1
      end if
    else
      level_1d_max(dim) = 0
    end if
  end do
  more_grids = .false.
!
!  Seek all vectors LEVEL_1D which satisfy the constraint:
!
!    LEVEL_MAX * LEVEL_WEIGHT_MIN_POS - sum ( LEVEL_WEIGHT ) 
!      < sum ( 1 <= I <= DIM_NUM ) LEVEL_WEIGHT(I) * LEVEL_1D(I)
!      <= LEVEL_MAX * LEVEL_WEIGHT_MIN_POS.
!
  do

    call sgmga_vcn_ordered ( dim_num, level_weight, level_1d_max, &
      level_1d, q_min, q_max, more_grids )

    if ( .not. more_grids ) then
      exit
    end if
!
!  Compute the combinatorial coefficient.
!
!   call sgmga_vcn_coef_naive ( dim_num, level_weight, level_1d_max, level_1d, &
!     q_min, q_max, coef )

    call sgmga_vcn_coef ( dim_num, level_weight, level_1d, q_max, coef )

    if ( coef == 0.0D+00 ) then
      cycle
    end if
!
!  Transform each 1D level to a corresponding 1D order.
!
    call level_growth_to_order ( dim_num, level_1d, rule, growth, order_1d )

    point_total_num = point_total_num + product ( order_1d(1:dim_num) )

  end do

  return
end
subroutine sgmga_unique_index ( dim_num, level_weight, level_max, rule, &
  growth, np, p, tol, point_num, point_total_num, sparse_unique_index )

!*****************************************************************************80
!
!! SGMGA_UNIQUE_INDEX maps nonunique to unique points of an SGMGA grid.
!
!  Discussion:
!
!    The sparse grid usually contains many points that occur in more
!    than one product grid.
!
!    When generating the point locations, it is easy to realize that a point
!    has already been generated.
!
!    But when it's time to compute the weights of the sparse grids, it is
!    necessary to handle situations in which weights corresponding to 
!    the same point generated in multiple grids must be collected together.
!
!    This routine generates ALL the points, including their multiplicities,
!    and figures out a mapping from them to the collapsed set of unique points.
!
!    This mapping can then be used during the weight calculation so that
!    a contribution to the weight gets to the right place.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 April 2011
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
!    Fabio Nobile, Raul Tempone, Clayton Webster,
!    An Anisotropic Sparse Grid Stochastic Collocation Method for Partial 
!    Differential Equations with Random Input Data,
!    SIAM Journal on Numerical Analysis,
!    Volume 46, Number 5, 2008, pages 2411-2442.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, real ( kind = 8 ) LEVEL_WEIGHT(DIM_NUM), the anisotropic weights.
!
!    Input, integer ( kind = 4 ) LEVEL_MAX, the maximum value of LEVEL.
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
!    Input, integer ( kind = 4 ) GROWTH(DIM_NUM), the growth in each dimension.
!    0, "DF", default growth associated with this quadrature rule;
!    1, "SL", slow linear, L+1;
!    2  "SO", slow linear odd, O=1+2((L+1)/2)
!    3, "ML", moderate linear, 2L+1;
!    4, "SE", slow exponential;
!    5, "ME", moderate exponential;
!    6, "FE", full exponential.
!
!    Input, integer ( kind = 4 ) NP(DIM_NUM), the number of parameters 
!    used by each rule.
!
!    Input, real ( kind = 8 ) P(*), the parameters needed by each rule.
!
!    Input, real ( kind = 8 ) TOL, the tolerance for point equality.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of unique points in 
!    the grid.
!
!    Input, integer ( kind = 4 ) POINT_TOTAL_NUM, the total number of points 
!    in the grid. 
!
!    Output, integer ( kind = 4 ) SPARSE_UNIQUE_INDEX(POINT_TOTAL_NUM), lists, 
!    for each (nonunique) point, the corresponding index of the same point in 
!    the unique listing.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_total_num

  real ( kind = 8 ) coef
  integer ( kind = 4 ) dim
  integer ( kind = 4 ) growth(dim_num)
  integer ( kind = 4 ) level
  integer ( kind = 4 ) level_1d(dim_num)
  integer ( kind = 4 ) level_1d_max(dim_num)
  integer ( kind = 4 ) level_max
  real ( kind = 8 ) level_weight(dim_num)
  real ( kind = 8 ) level_weight_min_pos
  logical more_grids
  logical more_points
  integer ( kind = 4 ) np(dim_num)
  integer ( kind = 4 ) order
  integer ( kind = 4 ) order_1d(dim_num)
  real ( kind = 8 ) p(*)
  integer ( kind = 4 ) p_index
  integer ( kind = 4 ) point
  integer ( kind = 4 ) point_index(dim_num)
  integer ( kind = 4 ) point_num
  integer ( kind = 4 ) point_total_num2
  real ( kind = 8 ), allocatable, dimension ( : ) :: points
  real ( kind = 8 ) q_max
  real ( kind = 8 ) q_min
  integer ( kind = 4 ) r8_floor
  real ( kind = 8 ) r8_huge
  real ( kind = 8 ) r8vec_min_pos
  integer ( kind = 4 ) rep
  integer ( kind = 4 ) rule(dim_num)
  integer ( kind = 4 ) seed
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: sparse_total_index
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: sparse_total_order
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: sparse_total_point
  integer ( kind = 4 ) sparse_unique_index(point_total_num)
  real ( kind = 8 ) tol
  integer ( kind = 4 ), allocatable, dimension ( : ) :: undx
!
!  Special cases.
!
  if ( level_max < 0 ) then
    return
  end if

  if ( level_max == 0 ) then
    sparse_unique_index(1) = 1
    return
  end if
!
!  Generate SPARSE_TOTAL_ORDER and SPARSE_TOTAL_INDEX arrays 
!  for the TOTAL set of points.
!
  allocate ( sparse_total_order(1:dim_num,1:point_total_num ) )
  allocate ( sparse_total_index(1:dim_num,1:point_total_num ) )

  point_total_num2 = 0
!
!  Initialization.
!
  level_weight_min_pos = r8vec_min_pos ( dim_num, level_weight )
  q_min = real ( level_max, kind = 8 ) * level_weight_min_pos &
    - sum ( level_weight(1:dim_num) )
  q_max = real ( level_max, kind = 8 ) * level_weight_min_pos
  do dim = 1, dim_num
    if ( 0.0D+00 < level_weight(dim) ) then
      level_1d_max(dim) = r8_floor ( q_max / level_weight(dim) ) + 1
      if ( q_max <= ( level_1d_max(dim) - 1 ) * level_weight(dim) ) then
        level_1d_max(dim) = level_1d_max(dim) - 1
      end if
    else
      level_1d_max(dim) = 0
    end if
  end do

  more_grids = .false.
!
!  Seek all vectors LEVEL_1D which satisfy the constraint:
!
!    LEVEL_MAX * LEVEL_WEIGHT_MIN_POS - sum ( LEVEL_WEIGHT ) 
!      < sum ( 1 <= I <= DIM_NUM ) LEVEL_WEIGHT(I) * LEVEL_1D(I)
!      <= LEVEL_MAX * LEVEL_WEIGHT_MIN_POS.
!
  do

    call sgmga_vcn_ordered ( dim_num, level_weight, level_1d_max, &
      level_1d, q_min, q_max, more_grids )

    if ( .not. more_grids ) then
      exit
    end if
!
!  Compute the combinatorial coefficient.
!
    call sgmga_vcn_coef ( dim_num, level_weight, level_1d, q_max, coef )

    if ( coef == 0.0D+00 ) then
      cycle
    end if
!
!  Transform each 1D level to a corresponding 1D order.
!
    call level_growth_to_order ( dim_num, level_1d, rule, growth, order_1d )
!
!  The inner loop generates a POINT of the GRID of the LEVEL.
!
    more_points = .false.

    do

      call vec_colex_next3 ( dim_num, order_1d, point_index, more_points )

      if ( .not. more_points ) then
        exit
      end if

      point_total_num2 = point_total_num2 + 1
      sparse_total_order(1:dim_num,point_total_num2) = order_1d(1:dim_num)
      sparse_total_index(1:dim_num,point_total_num2) = point_index(1:dim_num)

    end do

  end do
!
!  Now compute the coordinates of the TOTAL set of points.
!
  allocate ( sparse_total_point(1:dim_num,1:point_total_num) )
  sparse_total_point(1:dim_num,1:point_total_num) = r8_huge ( )

  level_weight_min_pos = r8vec_min_pos ( dim_num, level_weight )
  q_max = real ( level_max, kind = 8 ) * level_weight_min_pos

  p_index = 1

  do dim = 1, dim_num

    if ( 0.0D+00 < level_weight(dim) ) then
      level_1d_max(dim) = r8_floor ( q_max / level_weight(dim) ) + 1
      if ( q_max <= ( level_1d_max(dim) - 1 ) * level_weight(dim) ) then
        level_1d_max(dim) = level_1d_max(dim) - 1
      end if
    else
      level_1d_max(dim) = 0
    end if

    do level = 0, level_1d_max(dim)

      call level_growth_to_order ( 1, level, rule(dim), growth(dim), order )

      allocate ( points(1:order) )

      if ( rule(dim) == 1 ) then
        call clenshaw_curtis_compute_points_np ( order, np(dim), p(p_index), &
          points )
      else if ( rule(dim) == 2 ) then
        call fejer2_compute_points_np ( order, np(dim), p(p_index), points )
      else if ( rule(dim) == 3 ) then
        call patterson_lookup_points_np ( order, np(dim), p(p_index), points )
      else if ( rule(dim) == 4 ) then
        call legendre_compute_points_np ( order, np(dim), p(p_index), points )
      else if ( rule(dim) == 5 ) then
        call hermite_compute_points_np ( order, np(dim), p(p_index), points )
      else if ( rule(dim) == 6 ) then
        call gen_hermite_compute_points_np ( order, np(dim), p(p_index), &
          points )
      else if ( rule(dim) == 7 ) then
        call laguerre_compute_points_np ( order, np(dim), p(p_index), points )
      else if ( rule(dim) == 8 ) then
        call gen_laguerre_compute_points_np ( order, np(dim), p(p_index), &
          points )
      else if ( rule(dim) == 9 ) then
        call jacobi_compute_points_np ( order, np(dim), p(p_index), points )
      else if ( rule(dim) == 10 ) then
        call hermite_genz_keister_lookup_points_np ( order, np(dim), &
          p(p_index), points )
      else if ( rule(dim) == 11 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SGMGA_UNIQUE_INDEX - Fatal error!'
        write ( *, '(a)' ) '  Do not know how to assign points for rule 11.'
        stop
      else if ( rule(dim) == 12 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SGMGA_UNIQUE_INDEX - Fatal error!'
        write ( *, '(a)' ) '  Do not know how to assign points for rule 12.'
        stop
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SGMGA_UNIQUE_INDEX - Fatal error!'
        write ( *, '(a,i8)' ) '  Unexpected value of RULE = ', rule(dim)
        stop
      end if

      do point = 1, point_total_num
        if ( sparse_total_order(dim,point) == order ) then
          sparse_total_point(dim,point) = &
            points ( sparse_total_index(dim,point) )
        end if
      end do

      deallocate ( points )

    end do

    p_index = p_index + np(dim)

  end do
!
!  Merge points that are too close.
!
  seed = 123456789
 
  allocate ( undx(1:point_num) )

  call point_radial_tol_unique_index ( dim_num, point_total_num, &
    sparse_total_point, tol, seed, point_num, undx, sparse_unique_index )

  do point = 1, point_total_num
    rep = undx(sparse_unique_index(point))
    if ( point /= rep ) then
      sparse_total_point(1:dim_num,point) = sparse_total_point(1:dim_num,rep)
    end if
  end do
!
!  Construct an index that indicates the "rank" of the unique points.
!
  call point_unique_index ( dim_num, point_total_num, sparse_total_point, &
    point_num, undx, sparse_unique_index )

  deallocate ( undx )

  deallocate ( sparse_total_index )
  deallocate ( sparse_total_order )
  deallocate ( sparse_total_point )

  return
end
subroutine sgmga_vcn ( n, w, x, q_min, q_max, more )

!*****************************************************************************80
!
!! SGMGA_VCN returns the next constrained vector.
!
!  Discussion:
!
!    This function is intended to replace the "naive" version, now called
!    SGMGA_VCN_NAIVE, which is too slow for high dimensional problems.
!
!    For nonnegative vectors X of dimension N, and nonnegative
!    weights W, we define:
!
!      Q = sum ( 1 <= I <= N ) W(I) * X(I)
!
!    and seek X satisfying the constraint:
!
!      Q_MIN < Q <= Q_MAX
!
!    This routine returns, one at a time exactly those X which satisfy
!    the constraint.  No attempt is made to return the X values in 
!    any particular order as far as Q goes.  
!
!  Example:
! 
!        W               4.0 3.0 5.0       
!      MIN     16.0       0   0   0
!      ---     ----      -----------
!        1     20.0       5   0   0
!        2     19.0       4   1   0
!        3     18.0       3   2   0
!        4     17.0       2   3   0
!        5     20.0       2   4   0
!        6     19.0       1   5   0
!        7     18.0       0   6   0
!        8     17.0       3   0   1
!        9     20.0       3   1   1
!       10     19.0       2   2   1
!       11     18.0       1   3   1
!       12     17.0       0   4   1
!       13     20.0       0   5   1
!       14     18.0       2   0   2
!       15     17.0       1   1   2
!       16     20.0       1   2   2
!       17     19.0       0   3   2
!       18     19.0       1   0   3
!       19     18.0       0   1   3
!       20     20.0       0   0   4
!      ---     ----      ----------
!      MAX     20.0       6   7   5         
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 May 2010
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
!    Fabio Nobile, Raul Tempone, Clayton Webster,
!    An Anisotropic Sparse Grid Stochastic Collocation Method for Partial 
!    Differential Equations with Random Input Data,
!    SIAM Journal on Numerical Analysis,
!    Volume 46, Number 5, 2008, pages 2411-2442.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the vector.
!
!    Input, real ( kind = 8 ) W(N), the weights, which should be nonnegative.
!    At least one weight must be positive.
!
!    Input/output, integer ( kind = 4 ) X(N).  On first call, with 
!    MORE = FALSE, the input value of X is not important.  On subsequent calls, 
!    the input value of X should be the output value from the previous call.
!    On output, (with MORE = TRUE), the value of X will be the "next"
!    vector in the reverse lexicographical list of vectors that satisfy
!    the condition.  However, on output with MORE = FALSE, the vector
!    X is meaningless, because there are no more vectors in the list.
!
!    Input, real ( kind = 8 ) Q_MIN, Q_MAX, the lower and upper
!    limits on the sum.
!
!    Input/output, logical MORE.  On input, if the user has set MORE
!    FALSE, the user is requesting the initiation of a new sequence
!    of values.  If MORE is TRUE, then the user is requesting "more"
!    values in the current sequence.  On output, if MORE is TRUE,
!    then another value was found and returned in X, but if MORE is
!    FALSE, then there are no more values in the sequence, and X is
!    NOT the next value.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ), save :: dir
  integer ( kind = 4 ) i
  logical              more
  integer ( kind = 4 ), save :: n2
  integer ( kind = 4 ), save :: nstart
  real ( kind = 8 ) q_max
  real ( kind = 8 ) q_min
  integer ( kind = 4 ) r8_ceiling
  integer ( kind = 4 ) r8_floor
  real ( kind = 8 ) w(n)
  integer ( kind = 4 ) x(n)
  integer ( kind = 4 ), save, allocatable :: xmax(:)
  integer ( kind = 4 ), save :: xmin
!
!  Initialization for first call.
!
!  Allocate XMAX to remember the currently maximum possible value for each X.
!
!  Locate NSTART, the index of the first nonzero weight.
!  The algorithm is easier to program if the last index we look at
!  has a nonzero weight, so that it can always make up the remainder.
!
  if ( .not. more ) then

    if ( allocated ( xmax ) ) then
      deallocate ( xmax )
    end if

    allocate ( xmax(1:n) )

    nstart = - 1

    do i = 1, n

      if ( 0.0D+00 < w(i) ) then
        nstart = i
        exit
      end if

    end do
!
!  Theoretically, we could even handle the case where all weights are zero.
!  That case is ruled out elsewhere in this software, so I will not try
!  to deal with it here for now.
!
    if ( nstart == - 1 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SGMGA_VCN - Fatal error!'
      write ( *, '(a)' ) '  No weight is positive.'
      stop
    end if
!
!  Initialize X to zero, even the indices we ignore.
!
    x(1:n) = 0
    xmax(1:n) = 0
!
!  N2 points to our current index of interest.
!
    n2 = n + 1
    dir = - 1

    more = .true.

  end if
!
!  Look for the next solution vector X.
!
  do
!
!  If no more, the search is terminated.
!
    if ( .not. more ) then

      exit
!
!  DIR = -1, decrement N2, and, if possible, set X(N2) to XMIN.
!  DIR =  0, hold N2 at current value, and see if we can increment X(N2).
!
    else if ( dir == - 1 .or. dir == 0 ) then

      if ( dir == - 1 ) then
        n2 = n2 - 1
      end if

      if ( w(n2) == 0.0D+00 ) then

        xmin = 0
        xmax(n2) = 0

      else if ( nstart < n2 ) then

        xmin = 0

        xmax(n2) = r8_floor ( &
          ( q_max - dot_product ( w(n2+1:n), x(n2+1:n) ) ) / w(n2) )

      else if ( n2 == nstart .and. dir == - 1 ) then
         
        xmin = r8_ceiling ( &
          ( q_min - dot_product ( w(n2+1:n), x(n2+1:n) ) ) / w(n2) )
        xmin = max ( xmin, 0 )
        if ( dot_product ( w(1:n2-1), x(1:n2-1) ) + &
                           w(n2) * xmin + &
             dot_product ( w(n2+1:n), x(n2+1:n) ) <= q_min ) then
          xmin = xmin + 1
        end if

        x(n2) = xmin

        xmax(n2) = r8_floor ( &
          ( q_max - dot_product ( w(n2+1:n), x(n2+1:n) ) ) / w(n2) )

      end if

      if ( xmax(n2) < xmin ) then

        dir = + 1

      else

        if ( n2 == nstart ) then

          if ( dir == - 1 ) then
            dir = 0
            exit
          else if ( dir == 0 ) then
            x(n2) = x(n2) + 1
            if ( x(n2) <= xmax(n2) ) then
              exit
            else
              dir = + 1
            end if

          end if

        else

          x(n2) = xmin

        end if

      end if
!
!  DIR = + 1:
!  Try moving backwards to find an index N2 whose X we can increment.
!
    else if ( dir == + 1 ) then

      do

        if ( n2 == n ) then
          dir = 0
          more = .false.
          deallocate ( xmax )
          exit
        end if

        n2 = n2 + 1

        if ( 0.0D+00 < w(n2) ) then

          if ( x(n2) < xmax(n2) ) then
            x(n2) = x(n2) + 1
            dir = - 1
            exit
          end if

        end if

      end do

    end if

  end do 

  return
end
subroutine sgmga_vcn_coef ( dim_num, level_weight, x, q_max, coef )

!*****************************************************************************80
!
!! SGMGA_VCN_COEF returns the "next" constrained vector's coefficient.
!
!  Discussion:
!
!    The related code "SGMGA_VCN_COEF_NAIVE" represents a "naive" approach to
!    this calculation.  This code carries out the same calculation, but tries
!    to avoid the potential explosion in work that is exponential in the
!    spatial dimension for the naive approach.
!
!    We are considering nonnegative integer vectors X of dimension DIM_NUM 
!    for which the functional
!
!      Q(X) = sum ( 1 <= I <= DIM_NUM ) LEVEL_WEIGHT(I) * X(I)
!
!   satisfies the "Q" constraint:
!
!      Q_MIN < Q(X) <= Q_MAX
!
!    where LEVEL_WEIGHT is a vector of (essentially) positive weights.
!    Some, but not all of the entries of LEVEL_WEIGHT might be zero;
!    in that case, the corresponding values of X never vary, and do not
!    play a part in the following computation.
!
!    Supposing we have a suitable vector X, we now wish to count the 
!    number of distinct vectors Y which also satisfy the Q constraint
!    as well as the following "binary" constraint:
!
!      Y(I) = X(I) + B(I)
!
!    where every entry of B is 0 or 1.
!
!    Clearly, there are 2^DIM_NUM vectors Y which satisfy the binary
!    constraint, and a naive calculation would simply generate each 
!    possible Y, evaluate Q(Y), and if Y satisfies the Q constraint,
!    add it to the count.
!
!    But if we are considering even moderately large values of DIM_NUM, 
!    say 20 <= DIM_NUM, then the mere task of generating all possible 
!    Y vectors is burdensome.  If there are in fact likely to be only 
!    a few satisfactory Y vectors, (which depends on the values of 
!    Q_MAX and LEVEL_WEIGHT, of course) then it may be possible to
!    find and count them more rapidly.
!
!    This function attempts a more rapid computation by carrying out the
!    search in a particular order, and realizing that, in certain cases,
!    if a particular value Y* does not satisfy the Q constraint, then
!    a consecutive sequence of Y's following Y* also cannot satisfy the
!    constraint, and hence the search can jump over them.
!
!  Example:
!
!    DIM_NUM = 3
!    LEVEL_WEIGHT    3.0  2.0  1.0
!    Q_MAX    6.0
!
!    U = unsigned count 
!    S =   signed count returned as COEF
!                 
!    #   U  S   X1 X2 X3
!
!    1   8  0    0  0  0
!    2   7  1    0  0  1
!    3   6  0    0  0  2
!    4   5 -1    0  0  3
!    5   3 -1    0  0  4
!    6   2  0    0  0  5
!    7   1  1    0  0  6
!    8   6  0    0  1  0
!    9   5 -1    0  1  1
!   10   3 -1    0  1  2
!   11   2  0    0  1  3
!   12   1  1    0  1  4
!   13   3 -1    0  2  0
!   14   2  0    0  2  1
!   15   1  1    0  2  2
!   16   1  1    0  3  0
!   17   5 -1    1  0  0
!   18   3 -1    1  0  1
!   19   2  0    1  0  2
!   20   1  1    1  0  3
!   21   2  0    1  1  0
!   22   1  1    1  1  1
!   23   1  1    2  0  0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 May 2010
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
!    Fabio Nobile, Raul Tempone, Clayton Webster,
!    An Anisotropic Sparse Grid Stochastic Collocation Method for Partial 
!    Differential Equations with Random Input Data,
!    SIAM Journal on Numerical Analysis,
!    Volume 46, Number 5, 2008, pages 2411-2442.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, real ( kind = 8 ) LEVEL_WEIGHT(DIM_NUM), the weights.
!
!    Input, integer ( kind = 4 ) X(DIM_NUM), satisfies the Q constraint.
!
!    Input, real ( kind = 8 ) Q_MAX, the Q constraint maximum.
!
!    Output, real ( kind = 8 ) COEF, the combinatorial coefficient.
!
  implicit none

  integer ( kind = 4 ) dim_num

  integer ( kind = 4 ) b(dim_num)
  integer ( kind = 4 ) b_sum
  integer ( kind = 4 ) c
  real ( kind = 8 ) coef
  integer ( kind = 4 ) i
  real ( kind = 8 ) level_weight(dim_num)
  real ( kind = 8 ) q_max
  integer ( kind = 4 ) x(dim_num)

  c = 0
  b(1:dim_num) = 0

  do
!
!  Generate the next binary perturbation.
!
    i = 0

    do while ( i < dim_num )

      i = i + 1
!
!  If LEVEL_WEIGHT(I) == 0, B(I) is fixed at 0.  Next I.
!
      if ( level_weight(i) == 0.0D+00 ) then
!
!  If B(I) is 1, set it to 0.  Next I.
!
      else if ( b(i) == 1 ) then

        b(i) = 0
!
!  B(I) is 0.  Convert it to 1.
!
      else 

        b(i) = 1

        do
!
!  Does X + B satisfy the Q_MAX constraint?
!
          if ( dot_product ( level_weight, x + b ) <= q_max ) then          
            exit
          end if
!
!  If Q(X+B) now exceeds QMAX, B is rejected.  But we can also skip
!  all perturbations which agree with B through the I-th position.
!  To skip to the next "interesting" candidate, we essentially carry
!  out binary addition between B and a vector B' which has a single 1
!  in the I-th position.
!
          b(i) = 0

          do while ( i < dim_num )

            i = i + 1

            if ( level_weight(i) == 0.0D+00 ) then

            else if ( b(i) == 1 ) then
              b(i) = 0
            else
              b(i) = 1
              exit
            end if

          end do

        end do

        exit

      end if

    end do

    b_sum = sum ( b(1:dim_num) )
!
!  X+B is another solution to be counted.
!
    c = c + 1 - 2 * mod ( b_sum, 2 )
!
!  We're done if we've got back to 0.
!
    if ( b_sum == 0 ) then
      exit
    end if

  end do

  coef = real ( c, kind = 8 )

  return
end
subroutine sgmga_vcn_coef_naive ( dim_num, level_weight, x_max, x, q_min, &
  q_max, coef )

!*****************************************************************************80
!
!! SGMGA_VCN_COEF_NAIVE returns the "next" constrained vector's coefficient.
!
!  Discussion:
!
!    This function uses a naive approach to the computation, resulting in
!    a set of 2^DIM_NUM tasks.  Hence it is not suitable for cases where
!    DIM_NUM is moderately large.  The function SGMGA_VCN_COEF carries out
!    a more complicated but more efficient algorithm for the same computation.
!
!    We are given a nonnegative vector X of dimension DIM_NUM which satisfies:
!
!      sum ( 1 <= I <= DIM_NUM ) LEVEL_WEIGHT(I) * X(I) <= Q_MAX
!
!    This routine computes the appropriate coefficient for X in the
!    anisotropic sparse grid scheme.
!
!    The coefficient is calculated as follows:
!
!      Let B be a binary vector of length DIM_NUM, and let ||B|| represent
!      the sum of the entries of B.
!
!      Coef = sum ( all B such that X+B satisfies constraints ) (-1)^||B||
!
!    Since X+0 satisfies the constraint, there is always at least one 
!    summand.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 September 2009
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
!    Fabio Nobile, Raul Tempone, Clayton Webster,
!    An Anisotropic Sparse Grid Stochastic Collocation Method for Partial 
!    Differential Equations with Random Input Data,
!    SIAM Journal on Numerical Analysis,
!    Volume 46, Number 5, 2008, pages 2411-2442.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the vector.
!
!    Input, real ( kind = 8 ) LEVEL_WEIGHT(DIM_NUM), the anisotropic weights.
!
!    Input, integer ( kind = 4 ) X_MAX(DIM_NUM), the maximum
!    values allowed in each component.
!
!    Input, integer ( kind = 4 ) X(DIM_NUM), a point which satisifies 
!    the constraints.
!
!    Input, real ( kind = 8 ) Q_MIN, Q_MAX, the lower and upper
!    limits on the sum.
!
!    Output, real ( kind = 8 ) COEF, the combinatorial coefficient.
!
  implicit none

  integer ( kind = 4 ) dim_num

  integer ( kind = 4 ) b(dim_num)
  integer ( kind = 4 ) b_sum
  real ( kind = 8 ) coef
  real ( kind = 8 ) level_weight(dim_num)
  real ( kind = 8 ) q
  real ( kind = 8 ) q_max
  real ( kind = 8 ) q_min
  real ( kind = 8 ) r8_mop
  integer ( kind = 4 ) x(dim_num)
  integer ( kind = 4 ) x_max(dim_num)
  integer ( kind = 4 ) x2(dim_num)

  b(1:dim_num) = 0
  coef = 1.0D+00

  do
!
!  Generate the next binary perturbation.
!
    call binary_vector_next ( dim_num, b )
    b_sum = sum ( b(1:dim_num) )
!
!  We're done if we've got back to 0.
!
    if ( b_sum == 0 ) then
      exit
    end if
!
!  Perturb the vector.
!
    x2(1:dim_num) = x(1:dim_num) + b(1:dim_num)
!
!  Does it satisfy the XMAX constraint?
!  (THIS CHECK IS SURPRISINGLY NECESSARY, PARTLY BECAUSE OF ZERO WEIGHT).
!
    if ( any ( x_max(1:dim_num) < x2(1:dim_num) ) ) then
      cycle
    end if
!
!  Does it satisfy the Q_MAX constraint?
!  (We don't actually have to check Q_MIN!)
!
    q = dot_product ( level_weight(1:dim_num), &
      real ( x2(1:dim_num), kind = 8 ) )

    if ( q <= q_max ) then          
      coef = coef + r8_mop ( b_sum )
    end if

  end do

  return
end
subroutine sgmga_vcn_naive ( dim_num, level_weight, x_max, x, q_min, q_max, &
  more )

!*****************************************************************************80
!
!! SGMGA_VCN_NAIVE returns the next constrained vector.
!
!  Discussion:
!
!    This function uses a naive algorithm, which quickly becomes unsuitable
!    for higher dimensions.  The function SGMGA_VCN is an attempt at a more
!    efficient calculation of the same quantities.
!
!    We consider vectors X of dimension DIM_NUM satisfying:
!
!      0 <= X(1:DIM_NUM) <= X_MAX(1:DIM_NUM).
!
!    and define
!
!      Q = sum ( 1 <= I <= DIM_NUM ) LEVEL_WEIGHT(I) * X(I)
!
!    and seek X satisfying the constraint:
!
!      Q_MIN < Q <= Q_MAX
!
!    For sparse grid applications, we compute
!
!      LEVEL_WEIGHT_MIN_POS = minimum positive entry in LEVEL_WEIGHT
!
!    and assume there is an underlying LEVEL used to index the sets of 
!    constrained vectors, and that 
!
!      Q_MAX = LEVEL * LEVEL_WEIGHT_MIN_POS
!      Q_MIN = LEVEL - LEVEL_WEIGHT_MIN_POS * sum ( LEVEL_WEIGHT(:) )
!      X_MAX(I) = LEVEL * LEVEL_WEIGHT_MIN_POS / LEVEL_WEIGHT(I)
!
!    This routine returns, one at a time exactly those X which satisfy
!    the constraint.  No attempt is made to return the X values in 
!    any particular order as far as Q goes.  
!
!  Example:
!
!    LEVEL_WEIGHT:          1.000000        1.000000
!
!    Q_MIN:        0.000000
!    Q_MAX:        2.000000
!    X_MAX:                         2         2
!
!         1        1.000000         1         0
!         2        2.000000         2         0
!         3        1.000000         0         1
!         4        2.000000         1         1
!         5        2.000000         0         2
!
!    LEVEL_WEIGHT:          1.000000        2.000000
!
!    Q_MIN:       -1.000000
!    Q_MAX:        2.000000
!    X_MAX:                         2         1
!
!         1        0.000000         0         0
!         2        1.000000         1         0
!         3        2.000000         2         0
!         4        2.000000         0         1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 October 2009
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
!    Fabio Nobile, Raul Tempone, Clayton Webster,
!    An Anisotropic Sparse Grid Stochastic Collocation Method for Partial 
!    Differential Equations with Random Input Data,
!    SIAM Journal on Numerical Analysis,
!    Volume 46, Number 5, 2008, pages 2411-2442.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the vector.
!
!    Input, real ( kind = 8 ) LEVEL_WEIGHT(DIM_NUM), the anisotropic weights.
!
!    Input, integer ( kind = 4 ) X_MAX(DIM_NUM), the maximum
!    values allowed in each component off X.
!
!    Input/output, integer ( kind = 4 ) X(DIM_NUM).  On first call, with 
!    MORE = FALSE, the input value of X is not important.  On subsequent calls, 
!    the input value of X should be the output value from the previous call.
!    On output, (with MORE = TRUE), the value of X will be the "next"
!    vector in the reverse lexicographical list of vectors that satisfy
!    the condition.  However, on output with MORE = FALSE, the vector
!    X is meaningless, because there are no more vectors in the list.
!
!    Input, real ( kind = 8 ) Q_MIN, Q_MAX, the lower and upper
!    limits on the sum.
!
!    Input/output, logical MORE.  On input, if the user has set MORE
!    FALSE, the user is requesting the initiation of a new sequence
!    of values.  If MORE is TRUE, then the user is requesting "more"
!    values in the current sequence.  On output, if MORE is TRUE,
!    then another value was found and returned in X, but if MORE is
!    FALSE, then there are no more values in the sequence, and X is
!    NOT the next value.
!
  implicit none

  integer ( kind = 4 ) dim_num

  integer ( kind = 4 ) i
  real ( kind = 8 ) level_weight(dim_num)
  logical              more
  real ( kind = 8 ) q
  real ( kind = 8 ) q_max
  real ( kind = 8 ) q_min
  integer ( kind = 4 ) x(dim_num)
  integer ( kind = 4 ) x_max(dim_num)

  if ( .not. more ) then

    more = .true.

    x(1:dim_num) = 0

    q = dot_product ( level_weight(1:dim_num), &
      real ( x(1:dim_num), kind = 8 ) )

    if ( q_min < q .and. q <= q_max ) then
      return
    end if

  end if

  do

    i = 1

    do

      if ( x(i) < x_max(i) ) then
        exit
      end if

      if ( dim_num <= i ) then
        more = .false.
        return
      end if

      i = i + 1

    end do

    x(i) = x(i) + 1
    x(1:i-1) = 0

    q = dot_product ( level_weight(1:dim_num), &
      real ( x(1:dim_num), kind = 8 ) )

    if ( q_min < q .and. q <= q_max ) then
      exit
    end if

  end do

  return
end
subroutine sgmga_vcn_ordered ( dim_num, level_weight, x_max, x, q_min, q_max, &
  more )

!*****************************************************************************80
!
!! SGMGA_VCN_ORDERED returns the "next" constrained vector, with ordering.
!
!  Discussion:
!
!    We consider vectors X of dimension DIM_NUM satisfying:
!
!      0 <= X(1:DIM_NUM)
!
!    and define
!
!      Q = sum ( 1 <= I <= DIM_NUM ) LEVEL_WEIGHT(I) * X(I)
!
!    and seek X's satisfying the constraint:
!
!      Q_MIN < Q <= Q_MAX
!
!    For sparse grid applications, we compute
!
!      LEVEL_WEIGHT_MIN_POS = minimum positive entry in LEVEL_WEIGHT
!
!    and assume there is an underlying LEVEL used to index the sets of 
!    constrained vectors, and that 
!
!      Q_MAX = LEVEL * LEVEL_WEIGHT_MIN_POS
!      Q_MIN = LEVEL - LEVEL_WEIGHT_MIN_POS * sum ( LEVEL_WEIGHT(:) )
!
!    This function returns, one at a time exactly those X which satisfy
!    the constraint.
!
!    A weak ordering is imposed on the solution vectors.  This function 
!    subdivides the range Q_MIN through Q_MAX into subintervals of width 1, so 
!    that the X vectors returned are roughly sorted (or at least binned) 
!    by Q value.
!
!  Example:
!
!    If the weights are also integral, then the X vectors are in fact sorted 
!    by Q value:
!
!    LEVEL_WEIGHT:          1.000000        1.000000
!
!     Q_MIN:        0.000000
!     Q_MAX:        2.000000
!
!          1        1.000000         1         0
!          2        1.000000         0         1
!          3        2.000000         2         0
!          4        2.000000         1         1
!          5        2.000000         0         2
!
!    When the weights are not integral, then the X values are only BINNED
!    by Q value, that is, we first get all X's with Q values between Q_MIN
!    and Q_MIN+1, then Q_MIN+1 to Q_MIN+2 and so on, as demonstrated here:
!
!    LEVEL_WEIGHT:             1.5               1
!    Q_MIN:  0.5
!    Q_MAX:  3
!
!           1             1.5         1         0
!           2               1         0         1
!           3             2.5         1         1
!           4               2         0         2
!           5               3         2         0
!           6               3         0         3
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 October 2009
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
!    Fabio Nobile, Raul Tempone, Clayton Webster,
!    An Anisotropic Sparse Grid Stochastic Collocation Method for Partial 
!    Differential Equations with Random Input Data,
!    SIAM Journal on Numerical Analysis,
!    Volume 46, Number 5, 2008, pages 2411-2442.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the number of components in 
!    the vector.
!
!    Input, real ( kind = 8 ) LEVEL_WEIGHT(DIM_NUM), the anisotropic weights.
!
!    Input, integer ( kind = 4 ) X_MAX(DIM_NUM), the maximum values allowed in 
!    each component of X.
!
!    Input/output, integer ( kind = 4 ) X(DIM_NUM).  On first call, with 
!    MORE = FALSE, the input value of X is not important.  On subsequent calls, 
!    the input value of X should be the output value from the previous call.
!    On output, (with MORE = TRUE), the value of X will be the "next"
!    vector in the reverse lexicographical list of vectors that satisfy
!    the condition.  However, on output with MORE = FALSE, the vector
!    X is meaningless, because there are no more vectors in the list.
!
!    Input, real ( kind = 8 ) Q_MIN, Q_MAX, the lower and upper
!    limits on the sum.
!
!    Input/output, logical MORE.  On input, if the user has set MORE
!    FALSE, the user is requesting the initiation of a new sequence
!    of values.  If MORE is TRUE, then the user is requesting "more"
!    values in the current sequence.  On output, if MORE is TRUE,
!    then another value was found and returned in X, but if MORE is
!    FALSE, then there are no more values in the sequence, and X is
!    NOT the next value.
!
  implicit none

  integer ( kind = 4 ) dim_num

  logical              more
  real ( kind = 8 ) q_max
  real ( kind = 8 ), save :: q_max2
  real ( kind = 8 ) q_min
  real ( kind = 8 ), save :: q_min2
  real ( kind = 8 ) level_weight(dim_num)
  integer ( kind = 4 ) x(dim_num)
  integer ( kind = 4 ) x_max(dim_num)
!
!  On first call, initialize the subrange.
!
  if ( .not. more ) then
    q_min2 = q_min
    q_max2 = min ( q_min + 1.0D+00, q_max )
  end if
!
!  Call a lower level function to search the subrange.
!
  do

    call sgmga_vcn ( dim_num, level_weight, x, q_min2, q_max2, more )
!
!  If another solution was found, return it.
!
    if ( more ) then
      return
    end if
!
!  If the current subrange is exhausted, try to move to the next one.
!
    if ( q_max2 < q_max ) then
      q_min2 = q_max2
      q_max2 = min ( q_max2 + 1.0D+00, q_max )
!
!  If there are no more subranges, we're done.
!
    else

      exit

    end if

  end do

  return
end
subroutine sgmga_vcn_ordered_naive ( dim_num, level_weight, x_max, x, q_min, &
  q_max, more )

!*****************************************************************************80
!
!! SGMGA_VCN_ORDERED_NAIVE returns the "next" constrained vector, with ordering.
!
!  Discussion:
!
!    We consider vectors X of dimension DIM_NUM satisfying:
!
!      0 <= X(1:DIM_NUM) <= X_MAX(1:DIM_NUM).
!
!    and define
!
!      Q = sum ( 1 <= I <= DIM_NUM ) LEVEL_WEIGHT(I) * X(I)
!
!    and seek X's satisfying the constraint:
!
!      Q_MIN < Q <= Q_MAX
!
!    For sparse grid applications, we compute
!
!      LEVEL_WEIGHT_MIN_POS = minimum positive entry in LEVEL_WEIGHT
!
!    and assume there is an underlying LEVEL used to index the sets of 
!    constrained vectors, and that 
!
!      Q_MAX = LEVEL * LEVEL_WEIGHT_MIN_POS
!      Q_MIN = LEVEL - LEVEL_WEIGHT_MIN_POS * sum ( LEVEL_WEIGHT(:) )
!      X_MAX(I) = LEVEL * LEVEL_WEIGHT_MIN_POS / LEVEL_WEIGHT(I)
!
!    This function returns, one at a time exactly those X which satisfy
!    the constraint.
!
!    A weak ordering is imposed on the solution vectors.  This function 
!    subdivides the range Q_MIN through Q_MAX into subintervals of width 1, so 
!    that the X vectors returned are roughly sorted (or at least binned) 
!    by Q value.
!
!  Example:
!
!    If the weights are also integral, then the X vectors are in fact sorted 
!    by Q value:
!
!    LEVEL_WEIGHT:          1.000000        1.000000
!
!     Q_MIN:        0.000000
!     Q_MAX:        2.000000
!     X_MAX:                         2         2
!
!          1        1.000000         1         0
!          2        1.000000         0         1
!          3        2.000000         2         0
!          4        2.000000         1         1
!          5        2.000000         0         2
!
!    When the weights are not integral, then the X values are only BINNED
!    by Q value, that is, we first get all X's with Q values between Q_MIN
!    and Q_MIN+1, then Q_MIN+1 to Q_MIN+2 and so on, as demonstrated here:
!
!    LEVEL_WEIGHT:             1.5               1
!    Q_MIN:  0.5
!    Q_MAX:  3
!    X_MAX:                           2         3
!
!           1             1.5         1         0
!           2               1         0         1
!           3             2.5         1         1
!           4               2         0         2
!           5               3         2         0
!           6               3         0         3
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 October 2009
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
!    Fabio Nobile, Raul Tempone, Clayton Webster,
!    An Anisotropic Sparse Grid Stochastic Collocation Method for Partial 
!    Differential Equations with Random Input Data,
!    SIAM Journal on Numerical Analysis,
!    Volume 46, Number 5, 2008, pages 2411-2442.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the number of components in 
!    the vector.
!
!    Input, real ( kind = 8 ) LEVEL_WEIGHT(DIM_NUM), the anisotropic weights.
!
!    Input, integer ( kind = 4 ) X_MAX(DIM_NUM), the maximum values allowed in 
!    each component of X.
!
!    Input/output, integer ( kind = 4 ) X(DIM_NUM).  On first call, with 
!    MORE = FALSE, the input value of X is not important.  On subsequent calls, 
!    the input value of X should be the output value from the previous call.
!    On output, (with MORE = TRUE), the value of X will be the "next"
!    vector in the reverse lexicographical list of vectors that satisfy
!    the condition.  However, on output with MORE = FALSE, the vector
!    X is meaningless, because there are no more vectors in the list.
!
!    Input, real ( kind = 8 ) Q_MIN, Q_MAX, the lower and upper
!    limits on the sum.
!
!    Input/output, logical MORE.  On input, if the user has set MORE
!    FALSE, the user is requesting the initiation of a new sequence
!    of values.  If MORE is TRUE, then the user is requesting "more"
!    values in the current sequence.  On output, if MORE is TRUE,
!    then another value was found and returned in X, but if MORE is
!    FALSE, then there are no more values in the sequence, and X is
!    NOT the next value.
!
  implicit none

  integer ( kind = 4 ) dim_num

  logical more
  real ( kind = 8 ) q_max
  real ( kind = 8 ), save :: q_max2
  real ( kind = 8 ) q_min
  real ( kind = 8 ), save :: q_min2
  real ( kind = 8 ) level_weight(dim_num)
  integer ( kind = 4 ) x(dim_num)
  integer ( kind = 4 ) x_max(dim_num)
!
!  On first call, initialize the subrange.
!
  if ( .not. more ) then
    q_min2 = q_min
    q_max2 = min ( q_min + 1.0D+00, q_max )
  end if
!
!  Call a lower level function to search the subrange.
!
  do

    call sgmga_vcn_naive ( dim_num, level_weight, x_max, x, q_min2, &
      q_max2, more )
!
!  If another solution was found, return it.
!
    if ( more ) then
      return
    end if
!
!  If the current subrange is exhausted, try to move to the next one.
!
    if ( q_max2 < q_max ) then
      q_min2 = q_max2
      q_max2 = min ( q_max2 + 1.0D+00, q_max )
!
!  If there are no more subranges, we're done.
!
    else

      exit

    end if

  end do

  return
end
subroutine sgmga_weight ( dim_num, level_weight, level_max, rule, growth, &
  np, p, point_num, point_total_num, sparse_unique_index, sparse_weight )

!*****************************************************************************80
!
!! SGMGA_WEIGHT computes weights for an SGMGA rule.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 April 2011
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
!    Fabio Nobile, Raul Tempone, Clayton Webster,
!    An Anisotropic Sparse Grid Stochastic Collocation Method for Partial 
!    Differential Equations with Random Input Data,
!    SIAM Journal on Numerical Analysis,
!    Volume 46, Number 5, 2008, pages 2411-2442.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, real ( kind = 8 ) LEVEL_WEIGHT(DIM_NUM), the anisotropic weights.
!
!    Input, integer ( kind = 4 ) LEVEL_MAX, the maximum value of LEVEL.
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
!    Input, integer ( kind = 4 ) GROWTH(DIM_NUM), the growth in each dimension.
!    0, "DF", default growth associated with this quadrature rule;
!    1, "SL", slow linear, L+1;
!    2  "SO", slow linear odd, O=1+2((L+1)/2)
!    3, "ML", moderate linear, 2L+1;
!    4, "SE", slow exponential;
!    5, "ME", moderate exponential;
!    6, "FE", full exponential.
!
!    Input, integer ( kind = 4 ) NP(DIM_NUM), the number of parameters 
!    used by each rule.
!
!    Input, real ( kind = 8 ) P(*), the parameters needed by each rule.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of unique points in 
!    the grid.
!
!    Input, integer ( kind = 4 ) POINT_TOTAL_NUM, the total number of points 
!    in the grid. 
!
!    Input, integer ( kind = 4 ) SPARSE_UNIQUE_INDEX(POINT_TOTAL_NUM), lists, 
!    for each (nonunique) point, the corresponding index of the same point in 
!    the unique listing.
!
!    Output, real ( kind = 8 ) SPARSE_WEIGHT(POINT_NUM), the weights
!    associated with the sparse grid points.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num
  integer ( kind = 4 ) point_total_num

  real ( kind = 8 ) coef
  integer ( kind = 4 ) dim
  real ( kind = 8 ), allocatable, dimension ( : ) :: grid_weight
  integer ( kind = 4 ) growth(dim_num)
  integer ( kind = 4 ) level_1d(dim_num)
  integer ( kind = 4 ) level_1d_max(dim_num)
  integer ( kind = 4 ) level_max
  real ( kind = 8 ) level_weight(dim_num)
  real ( kind = 8 ) level_weight_min_pos
  logical more_grids
  integer ( kind = 4 ) np(dim_num)
  integer ( kind = 4 ) order
  integer ( kind = 4 ) order_1d(dim_num)
  integer ( kind = 4 ) order_nd
  real ( kind = 8 ) p(*)
  integer ( kind = 4 ) point_total
  integer ( kind = 4 ) point_unique
  real ( kind = 8 ) q_max
  real ( kind = 8 ) q_min
  integer ( kind = 4 ) r8_floor
  real ( kind = 8 ) r8vec_min_pos
  integer ( kind = 4 ) rule(dim_num)
  integer ( kind = 4 ) sparse_unique_index(point_total_num)
  real ( kind = 8 ) sparse_weight(point_num)

  sparse_weight(1:point_num) = 0.0D+00

  point_total = 0
!
!  Initialization.
!
  level_weight_min_pos = r8vec_min_pos ( dim_num, level_weight )
  q_min = real ( level_max, kind = 8 ) * level_weight_min_pos &
    - sum ( level_weight(1:dim_num) )
  q_max = real ( level_max, kind = 8 ) * level_weight_min_pos
  do dim = 1, dim_num
    if ( 0.0D+00 < level_weight(dim) ) then
      level_1d_max(dim) = r8_floor ( q_max / level_weight(dim) ) + 1
      if ( q_max <= ( level_1d_max(dim) - 1 ) * level_weight(dim) ) then
        level_1d_max(dim) = level_1d_max(dim) - 1
      end if
    else
      level_1d_max(dim) = 0
    end if
  end do
  more_grids = .false.
!
!  Seek all vectors LEVEL_1D which satisfy the constraint:
!
!    LEVEL_MAX * LEVEL_WEIGHT_MIN_POS - sum ( LEVEL_WEIGHT ) 
!      < sum ( 1 <= I <= DIM_NUM ) LEVEL_WEIGHT(I) * LEVEL_1D(I)
!      <= LEVEL_MAX * LEVEL_WEIGHT_MIN_POS.
!
  do

    call sgmga_vcn_ordered ( dim_num, level_weight, level_1d_max, &
      level_1d, q_min, q_max, more_grids )

    if ( .not. more_grids ) then
      exit
    end if
!
!  Compute the combinatorial coefficient.
!
    call sgmga_vcn_coef ( dim_num, level_weight, level_1d, q_max, coef )

    if ( coef == 0.0D+00 ) then
      cycle
    end if
!
!  Transform each 1D level to a corresponding 1D order.
!
    call level_growth_to_order ( dim_num, level_1d, rule, growth, order_1d )
!
!  The product of the 1D orders gives us the number of points in this grid.
!
    order_nd = product ( order_1d(1:dim_num) )
!
!  Compute the weights for this grid.
!
!  The correct transfer of data from the product grid to the sparse grid
!  depends on the fact that the product rule weights are stored under colex
!  order of the points, and this is the same ordering implicitly used in
!  generating the SPARSE_UNIQUE_INDEX array.
!
    allocate ( grid_weight(1:order_nd) )

    call sgmga_product_weight ( dim_num, order_1d, order_nd, rule, &
      np, p, grid_weight )
!
!  Add these weights to the rule.
!
    do order = 1, order_nd

      point_total = point_total + 1

      point_unique = sparse_unique_index(point_total)

      sparse_weight(point_unique) = sparse_weight(point_unique) &
        + coef * grid_weight(order)

    end do

    deallocate ( grid_weight )

  end do

  return
end
subroutine sgmga_write ( dim_num, level_weight, rule, growth, np, p, &
  point_num, sparse_weight, sparse_point, file_name )

!*****************************************************************************80
!
!! SGMGA_WRITE writes an SGMGA rule to several files.
!
!  Discussion:
!
!    The files are:
!    * the "A" file stores the LEVEL_WEIGHT values, as a 1 x DIM_NUM list.
!    * the "N" file stores the NP values, as a 1 x DIM_NUM list.
!    * the "P" file stores the P values, as a 1 x sum(NP(1:DIM_NUM) list.
!    * the "R" file stores the region, as a DIM_NUM x 2 list.
!    * the "W" file stores the weights as a 1 x POINT_NUM list;
!    * the "X" file stores the abscissas as a DIM_NUM x POINT_NUM list;
!
!    The entries in the "N" file are computed locally, although they would
!    normally be done by the user.
!
!    The entries in the "R" file are the two corners of the DIM_NUM dimensional
!    rectangle that constitutes the integration region.  Coordinates that
!    should be infinite are set to 1.0E+30.
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
!  Reference:
!
!    Fabio Nobile, Raul Tempone, Clayton Webster,
!    A Sparse Grid Stochastic Collocation Method for Partial Differential
!    Equations with Random Input Data,
!    SIAM Journal on Numerical Analysis,
!    Volume 46, Number 5, 2008, pages 2309-2345.
!
!    Fabio Nobile, Raul Tempone, Clayton Webster,
!    An Anisotropic Sparse Grid Stochastic Collocation Method for Partial 
!    Differential Equations with Random Input Data,
!    SIAM Journal on Numerical Analysis,
!    Volume 46, Number 5, 2008, pages 2411-2442.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, real ( kind = 8 ) LEVEL_WEIGHT(DIM_NUM), the anisotropic weights.
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
!    Input, integer ( kind = 4 ) GROWTH(DIM_NUM), the growth in each dimension.
!    0, "DF", default growth associated with this quadrature rule;
!    1, "SL", slow linear, L+1;
!    2  "SO", slow linear odd, O=1+2((L+1)/2)
!    3, "ML", moderate linear, 2L+1;
!    4, "SE", slow exponential;
!    5, "ME", moderate exponential;
!    6, "FE", full exponential.
!
!    Input, integer ( kind = 4 ) NP(DIM_NUM), the number of parameters 
!    used by each rule.
!
!    Input, real ( kind = 8 ) P(*), the parameters needed by each rule.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of unique points 
!    in the grid. 
!
!    Input, real ( kind = 8 ) SPARSE_WEIGHT(POINT_NUM), the weights.
!
!    Input, real ( kind = 8 ) SPARSE_POINT(DIM_NUM,POINT_NUM), the points.
!
!    Input, character ( len = * ) FILE_NAME, the main part of the file name.
!
  implicit none

  integer   ( kind = 4 ) dim_num
  integer   ( kind = 4 ) point_num

  integer   ( kind = 4 ) dim
  character ( len = *  ) file_name
  character ( len = 255) file_name_a
  character ( len = 255) file_name_n
  character ( len = 255) file_name_p
  character ( len = 255) file_name_r
  character ( len = 255) file_name_w
  character ( len = 255) file_name_x
  integer   ( kind = 4 ) growth(dim_num)
  real      ( kind = 8 ) level_weight(dim_num)
  integer   ( kind = 4 ) np(dim_num)
  integer   ( kind = 4 ) np_sum
  real      ( kind = 8 ) p(*)
  real      ( kind = 8 ) r8_huge
  integer   ( kind = 4 ) rule(dim_num)
  real      ( kind = 8 ) sparse_point(dim_num,point_num)
  real      ( kind = 8 ) sparse_region(dim_num,2)
  real      ( kind = 8 ) sparse_weight(point_num)
!
!  Construct values needed for the "R" file.
!
  do dim = 1, dim_num
    if ( rule(dim) == 1 ) then
      sparse_region(dim,1) = -1.0D+00
      sparse_region(dim,2) = +1.0D+00
    else if ( rule(dim) == 2 ) then
      sparse_region(dim,1) = -1.0D+00
      sparse_region(dim,2) = +1.0D+00
    else if ( rule(dim) == 3 ) then
      sparse_region(dim,1) = -1.0D+00
      sparse_region(dim,2) = +1.0D+00
    else if ( rule(dim) == 4 ) then
      sparse_region(dim,1) = -1.0D+00
      sparse_region(dim,2) = +1.0D+00
    else if ( rule(dim) == 5 ) then
      sparse_region(dim,1) = - r8_huge ( )
      sparse_region(dim,2) = + r8_huge ( )
    else if ( rule(dim) == 6 ) then
      sparse_region(dim,1) = - r8_huge ( )
      sparse_region(dim,2) = + r8_huge ( )
    else if ( rule(dim) == 7 ) then
      sparse_region(dim,1) = 0.0D+00
      sparse_region(dim,2) = r8_huge ( )
    else if ( rule(dim) == 8 ) then
      sparse_region(dim,1) = 0.0D+00
      sparse_region(dim,2) = r8_huge ( )
    else if ( rule(dim) == 9 ) then
      sparse_region(dim,1) = -1.0D+00
      sparse_region(dim,2) = +1.0D+00
    else if ( rule(dim) == 10 ) then
      sparse_region(dim,1) = - r8_huge ( )
      sparse_region(dim,2) = + r8_huge ( )
    else if ( rule(dim) == 11 ) then
      sparse_region(dim,1) = minval ( sparse_point(dim,1:point_num) )
      sparse_region(dim,2) = maxval ( sparse_point(dim,1:point_num) )
    else if ( rule(dim) == 12 ) then
      sparse_region(dim,1) = minval ( sparse_point(dim,1:point_num) )
      sparse_region(dim,2) = maxval ( sparse_point(dim,1:point_num) )
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SGMGA_WRITE - Fatal error!'
      write ( *, '(a,i8)' ) '  Unexpected value of RULE = ', rule(dim)
      stop
    end if
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SGMGA_WRITE:'

  file_name_a = trim ( file_name ) // '_a.txt'
  call r8mat_write ( file_name_a, 1, dim_num, level_weight )
  write ( *, '(a)' ) '  Wrote the A file = "' // trim ( file_name_a ) // '".'

  file_name_n = trim ( file_name ) // '_n.txt'
  call i4mat_write ( file_name_n, 1, dim_num, np )
  write ( *, '(a)' ) '  Wrote the N file = "' // trim ( file_name_n ) // '".'

  np_sum = sum ( np(1:dim_num) )
  file_name_p = trim ( file_name ) // '_p.txt'
  call r8mat_write ( file_name_p, 1, np_sum, p )
  write ( *, '(a)' ) '  Wrote the P file = "' // trim ( file_name_p ) // '".'

  file_name_r = trim ( file_name ) // '_r.txt'
  call r8mat_write ( file_name_r, dim_num, 2, sparse_region )
  write ( *, '(a)' ) '  Wrote the R file = "' // trim ( file_name_r ) // '".'

  file_name_w = trim ( file_name ) // '_w.txt'
  call r8mat_write ( file_name_w, 1, point_num, sparse_weight )
  write ( *, '(a)' ) '  Wrote the W file = "' // trim ( file_name_w ) // '".'

  file_name_x = trim ( file_name ) // '_x.txt'
  call r8mat_write ( file_name_x, dim_num, point_num, sparse_point )
  write ( *, '(a)' ) '  Wrote the X file = "' // trim ( file_name_x ) // '".'

  return
end
