subroutine comp_next ( n, k, a, more, h, t )

!*****************************************************************************80
!
!! COMP_NEXT computes the compositions of the integer N into K parts.
!
!  Discussion:
!
!    A composition of the integer N into K parts is an ordered sequence
!    of K nonnegative integers which sum to N.  The compositions (1,2,1)
!    and (1,1,2) are considered to be distinct.
!
!    The routine computes one composition on each call until there are no more.
!    For instance, one composition of 6 into 3 parts is
!    3+2+1, another would be 6+0+0.
!
!    On the first call to this routine, set MORE = FALSE.  The routine
!    will compute the first element in the sequence of compositions, and
!    return it, as well as setting MORE = TRUE.  If more compositions
!    are desired, call again, and again.  Each time, the routine will
!    return with a new composition.
!
!    However, when the LAST composition in the sequence is computed 
!    and returned, the routine will reset MORE to FALSE, signaling that
!    the end of the sequence has been reached.
!
!    This routine originally used a SAVE statement to maintain the
!    variables H and T.  I have decided that it is safer
!    to pass these variables as arguments, even though the user should
!    never alter them.  This allows this routine to safely shuffle
!    between several ongoing calculations.
!
!
!    There are 28 compositions of 6 into three parts.  This routine will
!    produce those compositions in the following order:
!
!     I         A
!     -     ---------
!     1     6   0   0
!     2     5   1   0
!     3     4   2   0
!     4     3   3   0
!     5     2   4   0
!     6     1   5   0
!     7     0   6   0
!     8     5   0   1
!     9     4   1   1
!    10     3   2   1
!    11     2   3   1
!    12     1   4   1
!    13     0   5   1
!    14     4   0   2
!    15     3   1   2
!    16     2   2   2
!    17     1   3   2
!    18     0   4   2
!    19     3   0   3
!    20     2   1   3
!    21     1   2   3
!    22     0   3   3
!    23     2   0   4
!    24     1   1   4
!    25     0   2   4
!    26     1   0   5
!    27     0   1   5
!    28     0   0   6
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 July 2008
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms for Computers and Calculators,
!    Second Edition,
!    Academic Press, 1978,
!    ISBN: 0-12-519260-6,
!    LC: QA164.N54.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the integer whose compositions are desired.
!
!    Input, integer ( kind = 4 ) K, the number of parts in the composition.
!
!    Input/output, integer ( kind = 4 ) A(K), the parts of the composition.
!
!    Input/output, logical MORE, set by the user to start the computation,
!    and by the routine to terminate it.
!
!    Input/output, integer ( kind = 4 )  H, T, two internal parameters needed 
!    for the computation.  The user should allocate space for these in the 
!    calling program, include them in the calling sequence, but never alter 
!    them!
!
  implicit none

  integer ( kind = 4 ) k

  integer ( kind = 4 ) a(k)
  integer ( kind = 4 ) h
  logical more
  integer ( kind = 4 ) n
  integer ( kind = 4 ) t
!
!  The first computation.
!
  if ( .not. more ) then

    t = n
    h = 0
    a(1) = n
    a(2:k) = 0
!
!  The next computation.
!
  else
!
!  If the first entry A(1) is positive, then set H to zero, 
!  so that when we increment H, it points to A(1); we will decrement A(1) by 1
!  and increment A(2).
!
    if ( 1 < t ) then
      h = 0
    end if
!
!  Otherwise, A(1) is 0.  Then by H + 1 is the entry we incremented last time.
!  Set H = H + 1, zero A(H), adding all but one of its value to A(1),
!  and incrementing A(H+1) by 1.
!
    h = h + 1
    t = a(h)
    a(h) = 0
    a(1) = t - 1
    a(h+1) = a(h+1) + 1

  end if
!
!  This is the last element of the sequence if all the
!  items are in the last slot.
!
  more = ( a(k) /= n )

  return
end
function i4_choose ( n, k )

!*****************************************************************************80
!
!! I4_CHOOSE computes the binomial coefficient C(N,K).
!
!  Discussion:
!
!    The value is calculated in such a way as to avoid overflow and
!    roundoff.  The calculation is done in integer arithmetic.
!
!    The formula used is:
!
!      C(N,K) = N! / ( K! * (N-K)! )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 June 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    ML Wolfson, HV Wright,
!    Algorithm 160:
!    Combinatorial of M Things Taken N at a Time,
!    Communications of the ACM,
!    Volume 6, Number 4, April 1963, page 161.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, K, are the values of N and K.
!
!    Output, integer ( kind = 4 ) I4_CHOOSE, the number of combinations of N
!    things taken K at a time.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_choose
  integer ( kind = 4 ) k
  integer ( kind = 4 ) mn
  integer ( kind = 4 ) mx
  integer ( kind = 4 ) n
  integer ( kind = 4 ) value

  mn = min ( k, n - k )

  if ( mn < 0 ) then

    value = 0

  else if ( mn == 0 ) then

    value = 1

  else

    mx = max ( k, n - k )
    value = mx + 1

    do i = 2, mn
      value = ( value * ( mx + i ) ) / i
    end do

  end if

  i4_choose = value

  return
end
subroutine sparse_grid_cc_me_size ( dim_num, level_max, point_num )

!*****************************************************************************80
!
!! SPARSE_GRID_CC_ME_SIZE: Clenshaw Curtis, Moderate Exponential Growth.
!
!  Discussion:
!
!    The Clenshaw-Curtis-Moderate family assumes that, for the underlying 1D
!    rules, a precision of 4*L+1 is needed at level L.  Therefore, the
!    lowest possible order Clenshaw-Curtis rule is chosen that will achieve
!    that precision.
!
!    The grid is defined as the sum of the product rules whose LEVEL
!    satisfies:
!
!      0 <= LEVEL <= LEVEL_MAX.
!
!    This calculation is much faster than a previous method.  It simply
!    computes the number of new points that are added at each level in the
!    1D rule, and then counts the new points at a given DIM_NUM dimensional
!    level vector as the product of the new points added in each dimension.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 February 2010
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
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) LEVEL_MAX, the maximum value of LEVEL.
!
!    Output, integer ( kind = 4 ) POINT_NUM, the total number of unique 
!    points in the grids.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) h
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l
  integer ( kind = 4 ) level
  integer ( kind = 4 ), allocatable :: level_1d(:)
  integer ( kind = 4 ) level_max
  logical more
  integer ( kind = 4 ), allocatable :: new_1d(:)
  integer ( kind = 4 ) o
  integer ( kind = 4 ), allocatable :: order_1d(:)
  integer ( kind = 4 ) p
  integer ( kind = 4 ) point_num
  integer ( kind = 4 ) t
!
!  Special case.
!
  if ( level_max < 0 ) then
    point_num = 0
    return
  end if

  if ( level_max == 0 ) then
    point_num = 1
    return
  end if
!
!  Count the points in the 1D rule.
!
  allocate ( order_1d(0:level_max) )

  order_1d(0) = 1
  do level = 1, level_max
    p = 3
    o = 3
    do while ( p < 4 * level + 1 )
      p = 2 * p + 1
      o = 2 * o + 1         
    end do
    order_1d(level) = o
  end do
!
!  Count the new points in the 1D rule.
!
  allocate ( new_1d(0:level_max) )

  new_1d(0) = 1
  do level = 1, level_max
    new_1d(level) = order_1d(level) - order_1d(level-1)
  end do
!
!  Count the number of points by counting the number of new points 
!  associated with each level vector.
!
  allocate ( level_1d(1:dim_num) )

  point_num = 0

  do level = 0, level_max

    more = .false.
    h = 0
    t = 0

    do

      call comp_next ( level, dim_num, level_1d, more, h, t )

      point_num = point_num + product ( new_1d(level_1d(1:dim_num)) )

      if ( .not. more ) then
        exit
      end if

    end do

  end do

  deallocate ( level_1d )
  deallocate ( new_1d )
  deallocate ( order_1d )

  return
end
subroutine sparse_grid_cc_se_size ( dim_num, level_max, point_num )

!*****************************************************************************80
!
!! SPARSE_GRID_CC_SE_SIZE: Clenshaw Curtis Slow Exponential Growth.
!
!  Discussion:
!
!    The sparse grid is defined as the sum of the product rules whose LEVEL
!    satisfies:
!
!      0 <= LEVEL <= LEVEL_MAX.
!
!    This calculation is much faster than a previous method.  It simply
!    computes the number of new points that are added at each level in the
!    1D rule, and then counts the new points at a given DIM_NUM dimensional
!    level vector as the product of the new points added in each dimension.
!
!    This approach will work for nested families, and may be extensible
!    to other families, and to mixed rules.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 December 2009
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
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) LEVEL_MAX, the maximum value of LEVEL.
!
!    Output, integer ( kind = 4 ) POINT_NUM, the total number of unique 
!    points in the grids.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) h
  integer ( kind = 4 ) l
  integer ( kind = 4 ) level
  integer ( kind = 4 ), allocatable :: level_1d(:)
  integer ( kind = 4 ) level_max
  logical more
  integer ( kind = 4 ), allocatable :: new_1d(:)
  integer ( kind = 4 ) o
  integer ( kind = 4 ) p
  integer ( kind = 4 ) point_num
  integer ( kind = 4 ) t
!
!  Special case.
!
  if ( level_max < 0 ) then
    point_num = 0
    return
  end if

  if ( level_max == 0 ) then
    point_num = 1
    return
  end if
!
!  Construct the vector that counts the new points in the 1D rule.
!
  allocate ( new_1d(0:level_max) )

  new_1d(0) = 1
  new_1d(1) = 2

  p = 3
  o = 3

  do l = 2, level_max
    p = 2 * l + 1
    if ( o < p ) then
      new_1d(l) = o - 1
      o = 2 * o - 1
    else
      new_1d(l) = 0
    end if
  end do
!
!  Count the number of points by counting the number of new points 
!  associated with each level vector.
!
  allocate ( level_1d(1:dim_num) )

  point_num = 0

  do level = 0, level_max

    more = .false.
    h = 0
    t = 0

    do

      call comp_next ( level, dim_num, level_1d, more, h, t )

      point_num = point_num + product ( new_1d(level_1d(1:dim_num)) )

      if ( .not. more ) then
        exit
      end if

    end do

  end do

  deallocate ( level_1d )
  deallocate ( new_1d )

  return
end
subroutine sparse_grid_cfn_e_size ( dim_num, level_max, point_num )

!*****************************************************************************80
!
!! SPARSE_GRID_CFN_E_SIZE; Closed Fully Nested, Exponential Growth.
!
!  Discussion:
!
!    The grid is defined as the sum of the product rules whose LEVEL
!    satisfies:
!
!      0 <= LEVEL <= LEVEL_MAX.
!
!    This calculation is much faster than a previous method.  It simply
!    computes the number of new points that are added at each level in the
!    1D rule, and then counts the new points at a given DIM_NUM dimensional
!    level vector as the product of the new points added in each dimension.
!
!    This approach will work for nested families, and may be extensible
!    to other families, and to mixed rules.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 December 2009
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
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) LEVEL_MAX, the maximum value of LEVEL.
!
!    Output, integer ( kind = 4 ) POINT_NUM, the total number of unique 
!    points in the grids.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) h
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l
  integer ( kind = 4 ) level
  integer ( kind = 4 ), allocatable :: level_1d(:)
  integer ( kind = 4 ) level_max
  logical more
  integer ( kind = 4 ), allocatable :: new_1d(:)
  integer ( kind = 4 ) point_num
  integer ( kind = 4 ) t
!
!  Special case.
!
  if ( level_max < 0 ) then
    point_num = 0
    return
  end if

  if ( level_max == 0 ) then
    point_num = 1
    return
  end if
!
!  Construct the vector that counts the new points in the 1D rule.
!
  allocate ( new_1d(0:level_max) )

  new_1d(0) = 1
  new_1d(1) = 2

  j = 1
  do l = 2, level_max
    j = j * 2
    new_1d(l) = j
  end do
!
!  Count the number of points by counting the number of new points 
!  associated with each level vector.
!
  allocate ( level_1d(1:dim_num) )

  point_num = 0

  do level = 0, level_max

    more = .false.
    h = 0
    t = 0

    do

      call comp_next ( level, dim_num, level_1d, more, h, t )

      point_num = point_num + product ( new_1d(level_1d(1:dim_num)) )

      if ( .not. more ) then
        exit
      end if

    end do

  end do

  deallocate ( level_1d )
  deallocate ( new_1d )

  return
end
subroutine sparse_grid_f2_se_size ( dim_num, level_max, point_num )

!*****************************************************************************80
!
!! SPARSE_GRID_F2_SE_SIZE: Fejer Type 2 Slow Exponential Growth.
!
!  Discussion:
!
!    The grid is defined as the sum of the product rules whose LEVEL
!    satisfies:
!
!      0 <= LEVEL <= LEVEL_MAX.
!
!    This calculation is much faster than a previous method.  It simply
!    computes the number of new points that are added at each level in the
!    1D rule, and then counts the new points at a given DIM_NUM dimensional
!    level vector as the product of the new points added in each dimension.
!
!    This approach will work for nested families, and may be extensible
!    to other families, and to mixed rules.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 December 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) LEVEL_MAX, the maximum value of LEVEL.
!
!    Output, integer ( kind = 4 ) POINT_NUM, the total number of unique 
!    points in the grids.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) h
  integer ( kind = 4 ) l
  integer ( kind = 4 ) level
  integer ( kind = 4 ), allocatable :: level_1d(:)
  integer ( kind = 4 ) level_max
  logical more
  integer ( kind = 4 ), allocatable :: new_1d(:)
  integer ( kind = 4 ) o
  integer ( kind = 4 ) p
  integer ( kind = 4 ) point_num
  integer ( kind = 4 ) t
!
!  Special case.
!
  if ( level_max < 0 ) then
    point_num = 0
    return
  end if

  if ( level_max == 0 ) then
    point_num = 1
    return
  end if
!
!  Construct the vector that counts the new points in the 1D rule.
!
  allocate ( new_1d(0:level_max) )

  new_1d(0) = 1

  p = 1
  o = 1
  do l = 1, level_max
    p = 2 * l + 1
    if ( o < p ) then
      new_1d(l) = o + 1
      o = 2 * o + 1
    else
      new_1d(l) = 0
    end if
  end do
!
!  Count the number of points by counting the number of new points 
!  associated with each level vector.
!
  allocate ( level_1d(1:dim_num) )

  point_num = 0

  do level = 0, level_max

    more = .false.
    h = 0
    t = 0

    do

      call comp_next ( level, dim_num, level_1d, more, h, t )

      point_num = point_num + product ( new_1d(level_1d(1:dim_num)) )

      if ( .not. more ) then
        exit
      end if

    end do

  end do

  deallocate ( level_1d )
  deallocate ( new_1d )

  return
end
subroutine sparse_grid_gp_me_size ( dim_num, level_max, point_num )

!*****************************************************************************80
!
!! SPARSE_GRID_GP_ME_SIZE: Gauss Patterson, Moderate Exponential Growth.
!
!  Discussion:
!
!    The Gauss-Patterson-Moderate family assumes that, for the underlying 1D
!    rules, a precision of 4*L+1 is needed at level L.  Therefore, the
!    lowest possible order Gauss-Patterson rule is chosen that will achieve
!    that precision.  This retains a combination of the advantages of
!    nestedness and high accuracy.
!
!    The grid is defined as the sum of the product rules whose LEVEL
!    satisfies:
!
!      0 <= LEVEL <= LEVEL_MAX.
!
!    This calculation is much faster than a previous method.  It simply
!    computes the number of new points that are added at each level in the
!    1D rule, and then counts the new points at a given DIM_NUM dimensional
!    level vector as the product of the new points added in each dimension.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 February 2010
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
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) LEVEL_MAX, the maximum value of LEVEL.
!
!    Output, integer ( kind = 4 ) POINT_NUM, the total number of unique 
!    points in the grids.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) h
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l
  integer ( kind = 4 ) level
  integer ( kind = 4 ), allocatable :: level_1d(:)
  integer ( kind = 4 ) level_max
  logical more
  integer ( kind = 4 ), allocatable :: new_1d(:)
  integer ( kind = 4 ) o
  integer ( kind = 4 ), allocatable :: order_1d(:)
  integer ( kind = 4 ) p
  integer ( kind = 4 ) point_num
  integer ( kind = 4 ) t
!
!  Special case.
!
  if ( level_max < 0 ) then
    point_num = 0
    return
  end if

  if ( level_max == 0 ) then
    point_num = 1
    return
  end if
!
!  Count the points in the 1D rule.
!
  allocate ( order_1d(0:level_max) )

  order_1d(0) = 1
  do level = 1, level_max
    p = 5
    o = 3
    do while ( p < 4 * level + 1 )
      p = 2 * p + 1
      o = 2 * o + 1         
    end do
    order_1d(level) = o
  end do
!
!  Count the new points in the 1D rule.
!
  allocate ( new_1d(0:level_max) )

  new_1d(0) = 1
  do level = 1, level_max
    new_1d(level) = order_1d(level) - order_1d(level-1)
  end do
!
!  Count the number of points by counting the number of new points 
!  associated with each level vector.
!
  allocate ( level_1d(1:dim_num) )

  point_num = 0

  do level = 0, level_max

    more = .false.
    h = 0
    t = 0

    do

      call comp_next ( level, dim_num, level_1d, more, h, t )

      point_num = point_num + product ( new_1d(level_1d(1:dim_num)) )

      if ( .not. more ) then
        exit
      end if

    end do

  end do

  deallocate ( level_1d )
  deallocate ( new_1d )
  deallocate ( order_1d )

  return
end
subroutine sparse_grid_gp_se_size ( dim_num, level_max, point_num )

!*****************************************************************************80
!
!! SPARSE_GRID_GP_SE_SIZE: Gauss Patterson, Slow Exponential Growth.
!
!  Discussion:
!
!    The Gauss-Patterson-Slow-Exponential family assumes that, for the 
!    underlying 1D rules, a precision of 2*L+1 is needed at level L.  Therefore, 
!    the lowest possible order Gauss-Patterson rule is chosen that will achieve
!    that precision.  This retains a combination of the advantages of
!    nestedness and high accuracy.
!
!    The grid is defined as the sum of the product rules whose LEVEL
!    satisfies:
!
!      0 <= LEVEL <= LEVEL_MAX.
!
!    This calculation is much faster than a previous method.  It simply
!    computes the number of new points that are added at each level in the
!    1D rule, and then counts the new points at a given DIM_NUM dimensional
!    level vector as the product of the new points added in each dimension.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 December 2009
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
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) LEVEL_MAX, the maximum value of LEVEL.
!
!    Output, integer ( kind = 4 ) POINT_NUM, the total number of unique 
!    points in the grids.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) h
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l
  integer ( kind = 4 ) level
  integer ( kind = 4 ), allocatable :: level_1d(:)
  integer ( kind = 4 ) level_max
  logical more
  integer ( kind = 4 ), allocatable :: new_1d(:)
  integer ( kind = 4 ) o
  integer ( kind = 4 ), allocatable :: order_1d(:)
  integer ( kind = 4 ) p
  integer ( kind = 4 ) point_num
  integer ( kind = 4 ) t
!
!  Special case.
!
  if ( level_max < 0 ) then
    point_num = 0
    return
  end if

  if ( level_max == 0 ) then
    point_num = 1
    return
  end if
!
!  Count the points in the 1D rule.
!
  allocate ( order_1d(0:level_max) )

  order_1d(0) = 1
  do level = 1, level_max
    p = 5
    o = 3
    do while ( p < 2 * level + 1 )
      p = 2 * p + 1
      o = 2 * o + 1         
    end do
    order_1d(level) = o
  end do
!
!  Count the new points in the 1D rule.
!
  allocate ( new_1d(0:level_max) )

  new_1d(0) = 1
  do level = 1, level_max
    new_1d(level) = order_1d(level) - order_1d(level-1)
  end do
!
!  Count the number of points by counting the number of new points 
!  associated with each level vector.
!
  allocate ( level_1d(1:dim_num) )

  point_num = 0

  do level = 0, level_max

    more = .false.
    h = 0
    t = 0

    do

      call comp_next ( level, dim_num, level_1d, more, h, t )

      point_num = point_num + product ( new_1d(level_1d(1:dim_num)) )

      if ( .not. more ) then
        exit
      end if

    end do

  end do

  deallocate ( level_1d )
  deallocate ( new_1d )
  deallocate ( order_1d )

  return
end
subroutine sparse_grid_ofn_e_size ( dim_num, level_max, point_num )

!*****************************************************************************80
!
!! SPARSE_GRID_OFN_E_SIZE: Open Fully Nested, Exponential Growth.
!
!  Discussion:
!
!    This calculation assumes that an exponential growth rule is being used,
!    that is, that the 1D rules have orders 1, 3, 7, 15, 31, and so on.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 December 2009
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
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) LEVEL_MAX, the maximum value of LEVEL.
!
!    Output, integer ( kind = 4 ) POINT_NUM, the total number of unique 
!    points in the grids.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) h
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l
  integer ( kind = 4 ) level
  integer ( kind = 4 ), allocatable :: level_1d(:)
  integer ( kind = 4 ) level_max
  logical more
  integer ( kind = 4 ), allocatable :: new_1d(:)
  integer ( kind = 4 ) point_num
  integer ( kind = 4 ) t
!
!  Special case.
!
  if ( level_max < 0 ) then
    point_num = 0
    return
  end if

  if ( level_max == 0 ) then
    point_num = 1
    return
  end if
!
!  Construct the vector that counts the new points in the 1D rule.
!
  allocate ( new_1d(0:level_max) )

  new_1d(0) = 1
  do l = 1, level_max
    new_1d(l) = 2 * new_1d(l-1)
  end do
!
!  Count the number of points by counting the number of new points 
!  associated with each level vector.
!
  allocate ( level_1d(1:dim_num) )

  point_num = 0

  do level = 0, level_max

    more = .false.
    h = 0
    t = 0

    do

      call comp_next ( level, dim_num, level_1d, more, h, t )

      point_num = point_num + product ( new_1d(level_1d(1:dim_num)) )

      if ( .not. more ) then
        exit
      end if

    end do

  end do

  deallocate ( level_1d )
  deallocate ( new_1d )

  return
end
subroutine sparse_grid_onn_e_size ( dim_num, level_max, point_num )

!*****************************************************************************80
!
!! SPARSE_GRID_ONN_E_SIZE: Open Non Nested, Exponential Growth.
!
!  Discussion:
!
!    This calculation assumes that an exponential growth rule is being used,
!    that is, that the 1D rules have orders 1, 3, 7, 15, 31, and so on.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 January 2009
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
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) LEVEL_MAX, the maximum value of LEVEL.
!
!    Output, integer ( kind = 4 ) POINT_NUM, the total number of unique 
!    points in the grids.
!
  implicit none

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) h
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l
  integer ( kind = 4 ) level
  integer ( kind = 4 ), allocatable :: level_1d(:)
  integer ( kind = 4 ) level_max
  integer ( kind = 4 ) level_min
  logical more
  integer ( kind = 4 ), allocatable :: new_1d(:)
  integer ( kind = 4 ), allocatable :: order_1d(:)
  integer ( kind = 4 ) point_num
  integer ( kind = 4 ) point_num2
  integer ( kind = 4 ) t
  integer ( kind = 4 ) temp
!
!  Special case.
!
  if ( level_max < 0 ) then
    point_num = 0
    return
  end if

  if ( level_max == 0 ) then
    point_num = 1
    return
  end if
!
!  Construct the 1D order vector.
!
  allocate ( order_1d(0:level_max) )

  temp = 2
  do l = 0, level_max
    order_1d(l) = temp - 1
    temp = temp * 2
  end do

  allocate ( level_1d(1:dim_num) )

  point_num = 0

  level_min = max ( 0, level_max + 1 - dim_num )

  do level = level_min, level_max

    more = .false.
    h = 0
    t = 0

    do

      call comp_next ( level, dim_num, level_1d, more, h, t )

      point_num = point_num + product ( order_1d(level_1d(1:dim_num)) )

      if ( .not. more ) then
        exit
      end if

    end do

  end do

  deallocate ( level_1d )
  deallocate ( order_1d )

  return
end
subroutine sparse_grid_onn_l_size ( dim_num, level_max, point_num )

!*****************************************************************************80
!
!! SPARSE_GRID_ONN_L_SIZE: Open Non Nested, Linear Growth.
!
!  Discussion:
!
!    This calculation assumes that a linear growth rule is being used,
!    that is, that the 1D rules have orders 1, 3, 5, 7, 9, and so on.
!
!    The grid is defined as the sum of the product rules whose LEVEL
!    satisfies:
!
!      0 <= LEVEL <= LEVEL_MAX.
!
!    This calculation is much faster than a previous method.  It simply
!    computes the number of new points that are added at each level in the
!    1D rule, and then counts the new points at a given DIM_NUM dimensional
!    level vector as the product of the new points added in each dimension.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 January 2009
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
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) LEVEL_MAX, the maximum value of LEVEL.
!
!    Output, integer ( kind = 4 ) POINT_NUM, the total number of unique 
!    points in the grids.
!
  implicit none

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) h
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l
  integer ( kind = 4 ) level
  integer ( kind = 4 ), allocatable :: level_1d(:)
  integer ( kind = 4 ) level_max
  integer ( kind = 4 ) level_min
  logical more
  integer ( kind = 4 ), allocatable :: new_1d(:)
  integer ( kind = 4 ), allocatable :: order_1d(:)
  integer ( kind = 4 ) point_num
  integer ( kind = 4 ) point_num2
  integer ( kind = 4 ) t
!
!  Special case.
!
  if ( level_max < 0 ) then
    point_num = 0
    return
  end if

  if ( level_max == 0 ) then
    point_num = 1
    return
  end if
!
!  Construct the 1D order vector.
!
  allocate ( order_1d(0:level_max) )

  do l = 0, level_max
    order_1d(l) = 2 * l + 1
  end do

  allocate ( level_1d(1:dim_num) )

  point_num = 0

  level_min = max ( 0, level_max + 1 - dim_num )

  do level = level_min, level_max

    more = .false.
    h = 0
    t = 0

    do

      call comp_next ( level, dim_num, level_1d, more, h, t )

      point_num = point_num + product ( order_1d(level_1d(1:dim_num)) )

      if ( .not. more ) then
        exit
      end if

    end do

  end do

  deallocate ( level_1d )
  deallocate ( order_1d )

  return
end
subroutine sparse_grid_own_e_size ( dim_num, level_max, point_num )

!*****************************************************************************80
!
!! SPARSE_GRID_OWN_E_SIZE: Open Weakly Nested, Exponential Growth.
!
!  Discussion:
!
!    This calculation assumes that an exponential growth rule is being used,
!    that is, that the 1D rules have orders 1, 3, 7, 15, 31, and so on.
!
!    This calculation assumes that the 1D family of quadrature rules 
!    contains only one repeated point, presumably the value 0.0.
!    This assumption holds for Gauss-Legendre, Gauss-Hermite and 
!    Generalized Gauss-Hermite rules.
!
!    The routine then counts the number of unique abscissas that will
!    be generated for a sparse grid of given dimension and level.
!
!    The computation is complicated.  It starts by counting just those
!    abscissas which have no 0.0 in them.  This is relatively easy, since
!    it is like counting the points in a sparse grid that uses open 
!    non-nested rules, but for which the order of each rule is reduced by 1.
!
!    Then we have to count the abscissas with one 0.0, two 0.0's and so
!    on to DIM_NUM zeros.  We are assuming this is an isotropic grid,
!    so for a particular number K of zeroes we only need to count the case
!    where the first K entries are zero, and multiply by C(DIM_NUM,K).
!
!    To count the number of entries with K zeroes, (and assuming 0 < K),
!    then, we essentially count the number of abscissas in an open 
!    non-nested rule as before, but modifed so that the minimum level is 0,
!    rather than LEVEL_MAX - DIM_NUM + 1.
!
!    I will mention that this was a rather difficult computation to
!    figure out!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 January 2009
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
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) LEVEL_MAX, the maximum value of LEVEL.
!
!    Output, integer ( kind = 4 ) POINT_NUM, the total number of unique 
!    points in the grids.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) dim_num2
  integer ( kind = 4 ) h
  integer ( kind = 4 ) i4_choose
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l
  integer ( kind = 4 ) level
  integer ( kind = 4 ), allocatable :: level_1d(:)
  integer ( kind = 4 ) level_max
  integer ( kind = 4 ) level_min
  logical more
  integer ( kind = 4 ), allocatable :: new_1d(:)
  integer ( kind = 4 ) point_num
  integer ( kind = 4 ) point_num2
  integer ( kind = 4 ) t
  integer ( kind = 4 ) temp
!
!  Special case.
!
  if ( level_max < 0 ) then
    point_num = 0
    return
  end if

  if ( level_max == 0 ) then
    point_num = 1
    return
  end if
!
!  Construct the vector that counts the new points in the "depleted" 1D rule.
!
  allocate ( new_1d(0:level_max) )

  new_1d(0) = 0

  temp = 4
  do l = 1, level_max
    new_1d(l) = temp - 2
    temp = temp * 2
  end do
!
!  Count the nonzero points in the full dimensional table with the usual
!  LEVEL_MIN restriction.
!
!  Then count the points with 1, 2, 3, ... DIM_NUM zeroes, by counting
!  the nonzero points in a DIM_NUM2 table, with LEVEL_MIN set to 0, and
!  multiplying by the appropriate combinatorial coefficient.
!
  point_num = 0

  do dim_num2 = dim_num, 0, -1

    if ( dim_num2 == dim_num ) then
      level_min = max ( 0, level_max - dim_num + 1 )
    else
      level_min = 0
    end if

    if ( dim_num2 == 0 ) then

      point_num2 = 1

    else

      allocate ( level_1d(1:dim_num2 ) )

      point_num2 = 0

      do level = level_min, level_max

        more = .false.
        h = 0
        t = 0

        do

          call comp_next ( level, dim_num2, level_1d, more, h, t )

          point_num2 = point_num2 + product ( new_1d(level_1d(1:dim_num2)) )

          if ( .not. more ) then
            exit
          end if

        end do

      end do

      deallocate ( level_1d )

    end if

    point_num = point_num + i4_choose ( dim_num, dim_num2 ) * point_num2

  end do

  deallocate ( new_1d )

  return
end
subroutine sparse_grid_own_l_size ( dim_num, level_max, point_num )

!*****************************************************************************80
!
!! SPARSE_GRID_OWN_L_SIZE: Open Weakly Nested, Linear Growth.
!
!  Discussion:
!
!    This calculation assumes that a linear growth rule is being used,
!    that is, that the 1D rules have orders 1, 3, 5, 7, 9, and so on.
!
!    This calculation assumes that the 1D family of quadrature rules 
!    contains only one repeated point, presumably the value 0.0.
!    This assumption holds for Gauss-Legendre, Gauss-Hermite and 
!    Generalized Gauss-Hermite rules.
!
!    The routine then counts the number of unique abscissas that will
!    be generated for a sparse grid of given dimension and level.
!
!    The computation is complicated.  It starts by counting just those
!    abscissas which have no 0.0 in them.  This is relatively easy, since
!    it is like counting the points in a sparse grid that uses open 
!    non-nested rules, but for which the order of each rule is reduced by 1.
!
!    Then we have to count the abscissas with one 0.0, two 0.0's and so
!    on to DIM_NUM zeros.  We are assuming this is an isotropic grid,
!    so for a particular number K of zeroes we only need to count the case
!    where the first K entries are zero, and multiply by C(DIM_NUM,K).
!
!    To count the number of entries with K zeroes, (and assuming 0 < K),
!    then, we essentially count the number of abscissas in an open 
!    non-nested rule as before, but modifed so that the minimum level is 0,
!    rather than LEVEL_MAX - DIM_NUM + 1.
!
!    I will mention that this was a rather difficult computation to
!    figure out!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 January 2009
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
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) LEVEL_MAX, the maximum value of LEVEL.
!
!    Output, integer ( kind = 4 ) POINT_NUM, the total number of unique 
!    points in the grids.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) dim_num2
  integer ( kind = 4 ) h
  integer ( kind = 4 ) i4_choose
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l
  integer ( kind = 4 ) level
  integer ( kind = 4 ), allocatable :: level_1d(:)
  integer ( kind = 4 ) level_max
  integer ( kind = 4 ) level_min
  logical more
  integer ( kind = 4 ), allocatable :: new_1d(:)
  integer ( kind = 4 ) point_num
  integer ( kind = 4 ) point_num2
  integer ( kind = 4 ) t
!
!  Special case.
!
  if ( level_max < 0 ) then
    point_num = 0
    return
  end if

  if ( level_max == 0 ) then
    point_num = 1
    return
  end if
!
!  Construct the vector that counts the new points in the "depleted" 1D rule.
!
  allocate ( new_1d(0:level_max) )

  new_1d(0) = 0

  do l = 1, level_max
    new_1d(l) = 2 * l
  end do
!
!  Count the nonzero points in the full dimensional table with the usual
!  LEVEL_MIN restriction.
!
!  Then count the points with 1, 2, 3, ... DIM_NUM zeroes, by counting
!  the nonzero points in a DIM_NUM2 table, with LEVEL_MIN set to 0, and
!  multiplying by the appropriate combinatorial coefficient.
!
  point_num = 0

  do dim_num2 = dim_num, 0, -1

    if ( dim_num2 == dim_num ) then
      level_min = max ( 0, level_max - dim_num + 1 )
    else
      level_min = 0
    end if

    if ( dim_num2 == 0 ) then

      point_num2 = 1

    else

      allocate ( level_1d(1:dim_num2 ) )

      point_num2 = 0

      do level = level_min, level_max

        more = .false.
        h = 0
        t = 0

        do

          call comp_next ( level, dim_num2, level_1d, more, h, t )

          point_num2 = point_num2 + product ( new_1d(level_1d(1:dim_num2)) )

          if ( .not. more ) then
            exit
          end if

        end do

      end do

      deallocate ( level_1d )

    end if

    point_num = point_num + i4_choose ( dim_num, dim_num2 ) * point_num2

  end do

  deallocate ( new_1d )

  return
end
subroutine sparse_grid_own_ls_size ( dim_num, level_max, point_num )

!*****************************************************************************80
!
!! SPARSE_GRID_OWN_LS_SIZE: Open Weakly Nested, Linear (Slow) Growth.
!
!  Discussion:
!
!    This calculation assumes that a slow linear growth rule is being used,
!    that is, that the 1D rules have orders 1, 3, 3, 5, 5, 7, 7, 9, 9, 
!    and so on.
!
!    This repetition of rules is permissible because the sparse grid rule
!    only requires that the 1D rule of level L have precision at least 2*L+1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 September 2012
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
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) LEVEL_MAX, the maximum value of LEVEL.
!
!    Output, integer ( kind = 4 ) POINT_NUM, the total number of unique 
!    points in the grids.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) dim_num2
  integer ( kind = 4 ) h
  integer ( kind = 4 ) i4_choose
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l
  integer ( kind = 4 ) level
  integer ( kind = 4 ), allocatable :: level_1d(:)
  integer ( kind = 4 ) level_max
  integer ( kind = 4 ) level_min
  logical more
  integer ( kind = 4 ), allocatable :: new_1d(:)
  integer ( kind = 4 ) point_num
  integer ( kind = 4 ) point_num2
  integer ( kind = 4 ) t
!
!  Special case.
!
  if ( level_max < 0 ) then
    point_num = 0
    return
  end if

  if ( level_max == 0 ) then
    point_num = 1
    return
  end if

  if ( dim_num == 1 ) then
    point_num = 1 + 2 * ( ( level_max + 1 ) / 2 )
    return
  end if
!
!  Construct the vector that counts the new points in the "depleted" 1D rule.
!
  allocate ( new_1d(0:level_max) )

  new_1d(0:level_max) = 0

  do l = 1, level_max, 2
    new_1d(l) = l + 1
  end do
!
!  Count the nonzero points in the full dimensional table with the usual
!  LEVEL_MIN restriction.
!
!  Then count the points with 1, 2, 3, ... DIM_NUM zeroes, by counting
!  the nonzero points in a DIM_NUM2 table, with LEVEL_MIN set to 0, and
!  multiplying by the appropriate combinatorial coefficient.
!
  point_num = 0

  do dim_num2 = dim_num, 0, -1

    if ( dim_num2 == dim_num ) then
      level_min = max ( 0, level_max - dim_num + 1 )
    else
      level_min = 0
    end if

    if ( dim_num2 == 0 ) then

      point_num2 = 1

    else

      allocate ( level_1d(1:dim_num2 ) )

      point_num2 = 0

      do level = level_min, level_max

        more = .false.
        h = 0
        t = 0

        do

          call comp_next ( level, dim_num2, level_1d, more, h, t )

          point_num2 = point_num2 + product ( new_1d(level_1d(1:dim_num2)) )

          if ( .not. more ) then
            exit
          end if

        end do

      end do

      deallocate ( level_1d )

    end if

    point_num = point_num + i4_choose ( dim_num, dim_num2 ) * point_num2

  end do

  deallocate ( new_1d )

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
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y

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
