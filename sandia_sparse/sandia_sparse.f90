subroutine abscissa_level_closed_nd ( level_max, dim_num, test_num, test_val, &
  test_level )

!*****************************************************************************80
!
!! ABSCISSA_LEVEL_CLOSED_ND: first level at which given abscissa is generated.
!
!  Discussion:
!
!    We assume an underlying product grid.  In each dimension, this product
!    grid has order 2^LEVEL_MAX + 1.
!
!    We will say a sparse grid has total level LEVEL if each point in the
!    grid has a total level of LEVEL or less.
!
!    The "level" of a point is determined as the sum of the levels of the
!    point in each spatial dimension.
!
!    The level of a point in a single spatial dimension I is determined as
!    the level, between 0 and LEVEL_MAX, at which the point's I'th index
!    would have been generated.
!
!
!    This description is terse and perhaps unenlightening.  Keep in mind
!    that the product grid is the product of 1D grids,
!    that the 1D grids are built up by levels, having
!    orders (total number of points ) 1, 3, 5, 9, 17, 33 and so on,
!    and that these 1D grids are nested, so that each point in a 1D grid
!    has a first level at which it appears.
!
!    Our procedure for generating the points of a sparse grid, then, is
!    to choose a value LEVEL_MAX, to generate the full product grid,
!    but then only to keep those points on the full product grid whose
!    LEVEL is less than or equal to LEVEL_MAX.
!
!
!    Note that this routine is really just testing out the idea of
!    determining the level.  Our true desire is to be able to start
!    with a value LEVEL, and determine, in a straightforward manner,
!    all the points that are generated exactly at that level, or
!    all the points that are generated up to and including that level.
!
!    This allows us to generate the new points to be added to one sparse
!    grid to get the next, or to generate a particular sparse grid at once.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LEVEL_MAX, controls the size of the
!    final sparse grid.
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) TEST_NUM, the number of points to be tested.
!
!    Input, integer ( kind = 4 ) TEST_VAL(DIM_NUM,TEST_NUM), the indices of
!    the points to be tested.  Normally, each index would be between 0 and
!    2^LEVEL_MAX.
!
!    Output, integer ( kind = 4 ) TEST_LEVEL(TEST_NUM), the value of LEVEL
!    at which the point would first be generated, assuming that a standard
!    sequence of nested grids is used.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) test_num

  integer ( kind = 4 ) index_to_level_closed
  integer ( kind = 4 ) j
  integer ( kind = 4 ) level_max
  integer ( kind = 4 ) order
  integer ( kind = 4 ) test_level(test_num)
  integer ( kind = 4 ) test_val(dim_num,test_num)
!
!  Special case: LEVEL_MAX = 0.
!
  if ( level_max <= 0 ) then

    test_level(1:test_num) = 0

  else

    order = 2**level_max + 1

    do j = 1, test_num

      test_level(j) = index_to_level_closed ( dim_num, test_val(1:dim_num,j), &
        order, level_max )

    end do

  end if

  return
end
subroutine abscissa_level_open_nd ( level_max, dim_num, test_num, test_val, &
  test_level )

!*****************************************************************************80
!
!! ABSCISSA_LEVEL_OPEN_ND: first level at which given abscissa is generated.
!
!  Discussion:
!
!    We assume an underlying product grid.  In each dimension, this product
!    grid has order 2^(LEVEL_MAX+1) - 1.
!
!    We will say a sparse grid has total level LEVEL if each point in the
!    grid has a total level of LEVEL or less.
!
!    The "level" of a point is determined as the sum of the levels of the
!    point in each spatial dimension.
!
!    The level of a point in a single spatial dimension I is determined as
!    the level, between 0 and LEVEL_MAX, at which the point's I'th index
!    would have been generated.
!
!
!    This description is terse and perhaps unenlightening.  Keep in mind
!    that the product grid is the product of 1D grids,
!    that the 1D grids are built up by levels, having
!    orders (total number of points ) 1, 3, 7, 15, 31 and so on,
!    and that these 1D grids are nested, so that each point in a 1D grid
!    has a first level at which it appears.
!
!    Our procedure for generating the points of a sparse grid, then, is
!    to choose a value LEVEL_MAX, to generate the full product grid,
!    but then only to keep those points on the full product grid whose
!    LEVEL is less than or equal to LEVEL_MAX.
!
!
!    Note that this routine is really just testing out the idea of
!    determining the level.  Our true desire is to be able to start
!    with a value LEVEL, and determine, in a straightforward manner,
!    all the points that are generated exactly at that level, or
!    all the points that are generated up to and including that level.
!
!    This allows us to generate the new points to be added to one sparse
!    grid to get the next, or to generate a particular sparse grid at once.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LEVEL_MAX, controls the size of the final
!    sparse grid.
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) TEST_NUM, the number of points to be tested.
!
!    Input, integer ( kind = 4 ) TEST_VAL(DIM_NUM,TEST_NUM), the indices of
!    the points to be tested.  Normally, each index would be between
!    0 and 2^(LEVEL_MAX+1).
!
!    Output, integer ( kind = 4 ) TEST_LEVEL(TEST_NUM), the value of LEVEL at
!    which the point would first be generated, assuming that a standard
!    sequence of nested grids is used.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) test_num

  integer ( kind = 4 ) index_to_level_open
  integer ( kind = 4 ) j
  integer ( kind = 4 ) level_max
  integer ( kind = 4 ) order
  integer ( kind = 4 ) test_level(test_num)
  integer ( kind = 4 ) test_val(dim_num,test_num)
!
!  Special case: LEVEL_MAX = 0.
!
  if ( level_max == 0 ) then
    test_level(1:test_num) = 0
    return
  end if

  order = 2**( level_max + 1 ) - 1

  do j = 1, test_num

    test_level(j) = index_to_level_open ( dim_num, test_val(1:dim_num,j), &
      order, level_max )

  end do

  return
end
function cc_abscissa ( order, i )

!*****************************************************************************80
!
!! CC_ABSCISSA returns the I-th abscissa for the Clenshaw Curtis rule.
!
!  Discussion:
!
!    Our convention is that the abscissas are numbered from left to
!    right.
!
!    This rule is defined on [-1,1].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ORDER, the order of the Clenshaw Curtis rule.
!    1 <= ORDER.
!
!    Input, integer ( kind = 4 ) I, the index of the desired abscissa.
!    1 <= I <= ORDER.
!
!    Output, real ( kind = 8 ) CC_ABSCISSA, the value of the I-th
!    abscissa in the Clenshaw Curtis rule of order ORDER.
!
  implicit none

  real ( kind = 8 ) cc_abscissa
  integer ( kind = 4 ) i
  integer ( kind = 4 ) order
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r8_huge
  real ( kind = 8 ) value

  if ( order < 1 ) then

    value = - r8_huge ( )

  else if ( i < 1 .or. order < i ) then

    value = - r8_huge ( )

  else if ( order == 1 ) then

    value = 0.0D+00

  else if ( 2 * ( order - i ) == order - 1 ) then

    value = 0.0D+00

  else

    value = cos ( real ( order - i, kind = 8 ) * pi &
                / real ( order - 1, kind = 8 ) )

  end if

  cc_abscissa = value

  return
end
subroutine cc_weights ( n, w )

!*****************************************************************************80
!
!! CC_WEIGHTS computes Clenshaw Curtis weights.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Charles Clenshaw, Alan Curtis,
!    A Method for Numerical Integration on an Automatic Computer,
!    Numerische Mathematik,
!    Volume 2, Number 1, December 1960, pages 197-205.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the rule.
!
!    Output, real ( kind = 8 ) W(N), the weights of the rule.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) :: pi = 3.141592653589793D+00
  real ( kind = 8 ) theta(n)
  real ( kind = 8 ) w(n)

  if ( n == 1 ) then
    w(1) = 2.0D+00
    return
  end if

  do i = 1, n
    theta(i) = real ( i - 1, kind = 8 ) * pi &
             / real ( n - 1, kind = 8 )
  end do

  do i = 1, n

    w(i) = 1.0D+00

    do j = 1, ( n - 1 ) / 2

      if ( 2 * j == ( n - 1 ) ) then
        b = 1.0D+00
      else
        b = 2.0D+00
      end if

      w(i) = w(i) - b * cos ( 2.0D+00 * real ( j, kind = 8 ) * theta(i) ) &
           / real ( 4 * j * j - 1, kind = 8 )

    end do

  end do

  w(1)     =           w(1)     / real ( n - 1, kind = 8 )
  w(2:n-1) = 2.0D+00 * w(2:n-1) / real ( n - 1, kind = 8 )
  w(n)     =           w(n)     / real ( n - 1, kind = 8 )

  return
end
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
!    FORTRAN90 version by John Burkardt
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
!    calling program, include them in the calling sequence, but never
!    alter them!
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

    if ( 1 < t ) then
      h = 0
    end if

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
function f1_abscissa ( order, i )

!*****************************************************************************80
!
!! F1_ABSCISSA returns the I-th abscissa for the Fejer type 1 rule.
!
!  Discussion:
!
!    Our convention is that the abscissas are numbered from left to
!    right.
!
!    This rule is defined on [-1,+1].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ORDER, the order of the Fejer type 1 rule.
!    1 <= ORDER.
!
!    Input, integer ( kind = 4 ) I, the index of the desired abscissa.
!    1 <= I <= ORDER.
!
!    Output, real ( kind = 8 ) F1_ABSCISSA, the value of the I-th
!    abscissa in the Fejer type 1 rule of order ORDER.
!
  implicit none

  real ( kind = 8 ) f1_abscissa
  integer ( kind = 4 ) i
  integer ( kind = 4 ) order
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) value

  if ( order < 1 ) then
    value = - huge ( value )
  else if ( i < 1 .or. order < i ) then
    value = - huge ( value )
  else if ( order == 1 ) then
    value = 0.0D+00
  else if ( 2 * ( 2 * order + 1 - 2 * i ) == 2 * order ) then
    value = 0.0D+00
  else
    value = cos ( real ( 2 * order + 1 - 2 * i, kind = 8 ) * pi &
                / real ( 2 * order,             kind = 8 ) )
  end if

  f1_abscissa = value

  return
end
subroutine f1_weights ( order, w )

!*****************************************************************************80
!
!! F1_WEIGHTS computes weights for a Fejer type 1 rule.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!    Walter Gautschi,
!    Numerical Quadrature in the Presence of a Singularity,
!    SIAM Journal on Numerical Analysis,
!    Volume 4, Number 3, 1967, pages 357-362.
!
!    Joerg Waldvogel,
!    Fast Construction of the Fejer and Clenshaw-Curtis Quadrature Rules,
!    BIT Numerical Mathematics,
!    Volume 43, Number 1, 2003, pages 1-18.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ORDER, the order.
!
!    Output, real ( kind = 8 ) W(ORDER), the weights.
!
  implicit none

  integer ( kind = 4 ) order

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) :: pi = 3.141592653589793D+00
  real ( kind = 8 ) theta(order)
  real ( kind = 8 ) w(order)

  if ( order < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'F1_WEIGHTS - Fatal error!'
    write ( *, '(a)' ) '  ORDER < 1.'
    stop
  end if

  if ( order == 1 ) then
    w(1) = 2.0D+00
    return
  end if

  do i = 1, order
    theta(i) = real ( 2 * ( order + 1 - i ) - 1, kind = 8 ) * pi &
             / real ( 2 * order,     kind = 8 )
  end do

  do i = 1, order
    w(i) = 1.0D+00
    do j = 1, ( order / 2 )
      w(i) = w(i) - 2.0D+00 &
        * cos ( 2.0D+00 * real ( j, kind = 8 ) * theta(i) ) &
        / real ( 4 * j * j - 1, kind = 8 )
    end do
  end do

  w(1:order) = 2.0D+00 * w(1:order) / real ( order, kind = 8 )

  return
end
function f2_abscissa ( order, i )

!*****************************************************************************80
!
!! F2_ABSCISSA returns the I-th abscissa for the Fejer type 2 rule.
!
!  Discussion:
!
!    Our convention is that the abscissas are numbered from left to
!    right.
!
!    This rule is defined on [-1,+1].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ORDER, the order of the Fejer type 2 rule.
!    1 <= ORDER.
!
!    Input, integer ( kind = 4 ) I, the index of the desired abscissa.
!    1 <= I <= ORDER.
!
!    Output, real ( kind = 8 ) F2_ABSCISSA, the value of the I-th
!    abscissa in the Fejer type 2 rule of order ORDER.
!
  implicit none

  real ( kind = 8 ) f2_abscissa
  integer ( kind = 4 ) i
  integer ( kind = 4 ) order
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) value

  if ( order < 1 ) then
    value = - huge ( value )
  else if ( i < 1 .or. order < i ) then
    value = - huge ( value )
  else if ( order == 1 ) then
    value = 0.0D+00
  else if ( 2 * ( order + 1 - i ) == order + 1 ) then
    value = 0.0D+00
  else
    value = cos ( real ( order + 1 - i, kind = 8 ) * pi &
                / real ( order + 1,     kind = 8 ) )
  end if

  f2_abscissa = value

  return
end
subroutine f2_weights ( order, w )

!*****************************************************************************80
!
!! F2_WEIGHTS computes weights for a Fejer type 2 rule.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!    Walter Gautschi,
!    Numerical Quadrature in the Presence of a Singularity,
!    SIAM Journal on Numerical Analysis,
!    Volume 4, Number 3, 1967, pages 357-362.
!
!    Joerg Waldvogel,
!    Fast Construction of the Fejer and Clenshaw-Curtis Quadrature Rules,
!    BIT Numerical Mathematics,
!    Volume 43, Number 1, 2003, pages 1-18.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ORDER, the order.
!
!    Output, real ( kind = 8 ) W(ORDER), the weights.
!
  implicit none

  integer ( kind = 4 ) order

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) p
  real ( kind = 8 ) :: pi = 3.141592653589793D+00
  real ( kind = 8 ) theta(order)
  real ( kind = 8 ) w(order)

  if ( order < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'F2_WEIGHTS - Fatal error!'
    write ( *, '(a)' ) '  ORDER < 1.'
    stop
  end if

  if ( order == 1 ) then
    w(1) = 2.0D+00
    return
  else if ( order == 2 ) then
    w(1:2) = 1.0D+00
    return
  end if

  do i = 1, order
    theta(i) = real ( order + 1 - i, kind = 8 ) * pi &
             / real ( order + 1, kind = 8 )
  end do

  do i = 1, order

    w(i) = 1.0D+00

    do j = 1, ( ( order - 1 ) / 2 )
      w(i) = w(i) - 2.0D+00 &
        * cos ( 2.0D+00 * real ( j, kind = 8 ) * theta(i) ) &
        / real ( 4 * j * j - 1, kind = 8 )
    end do

    if ( 2 < order ) then
      p = 2.0D+00 * real ( ( ( order + 1 ) / 2 ), kind = 8 ) - 1.0D+00
      w(i) = w(i) - cos ( ( p + 1.0D+00 ) * theta(i) ) / p
    end if

  end do

  w(1:order) = 2.0D+00 * w(1:order) / real ( order + 1, kind = 8 )

  return
end
subroutine gh_abscissa ( dim_num, point_num, grid_index, grid_base, &
  grid_point )

!*****************************************************************************80
!
!! GH_ABSCISSA sets abscissas for multidimensional Gauss Hermite quadrature.
!
!  Discussion:
!
!    The "nesting" as it occurs for Gauss Hermite sparse grids simply
!    involves the use of a specified set of permissible orders for the
!    rule.
!
!    The X array lists the (complete) Gauss Hermite abscissas for rules
!    of order 1, 3, 7, 15, 31, 63, and 127 in order.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, integer ( kind = 4 ) GRID_INDEX(DIM_NUM,POINT_NUM), for each
!    point and dimension, the index of the abscissa.
!
!    Input, integer ( kind = 4 ) GRID_BASE(DIM_NUM), the "base" of the
!    rule being used in each dimension.
!
!    Output, real ( kind = 8 ) GRID_POINT(DIM_NUM), the grid points of
!    abscissas.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) grid_base(dim_num)
  integer ( kind = 4 ) grid_index(dim_num,point_num)
  real ( kind = 8 ) grid_point(dim_num,point_num)
  integer ( kind = 4 ) i4_log_2
  integer ( kind = 4 ) level
  integer ( kind = 4 ) point
  integer ( kind = 4 ) pointer
  integer ( kind = 4 ), dimension ( 0:7 ) :: skip = (/ &
    0, 1, 4, 11, 26, 57, 120, 247 /)
  real ( kind = 8 ), dimension ( 247 ) :: x = (/ &
    0.0D+00, &
   -0.122474487139158904909864203735D+01, &
    0.0D+00, &
    0.122474487139158904909864203735D+01, &
   -0.265196135683523349244708200652D+01, &
   -0.167355162876747144503180139830D+01, &
   -0.816287882858964663038710959027D+00, &
    0.0D+00, &
    0.816287882858964663038710959027D+00, &
    0.167355162876747144503180139830D+01, &
    0.265196135683523349244708200652D+01, &
   -0.449999070730939155366438053053D+01, &
   -0.366995037340445253472922383312D+01, &
   -0.296716692790560324848896036355D+01, &
   -0.232573248617385774545404479449D+01, &
   -0.171999257518648893241583152515D+01, &
   -0.113611558521092066631913490556D+01, &
   -0.565069583255575748526020337198D+00, &
    0.0D+00, &
    0.565069583255575748526020337198D+00, &
    0.113611558521092066631913490556D+01, &
    0.171999257518648893241583152515D+01, &
    0.232573248617385774545404479449D+01, &
    0.296716692790560324848896036355D+01, &
    0.366995037340445253472922383312D+01, &
    0.449999070730939155366438053053D+01, &
   -6.9956801237185402753248521473232D+00, &
   -6.2750787049428601427036567812530D+00, &
   -5.6739614446185883296332558789276D+00, &
   -5.1335955771123807045862968913996D+00, &
   -4.6315595063128599420667997654336D+00, &
   -4.1562717558181451724831352315314D+00, &
   -3.7007434032314694224497164589673D+00, &
   -3.2603207323135408104645401509648D+00, &
   -2.8316804533902054557015640151425D+00, &
   -2.4123177054804201051740184582119D+00, &
   -2.0002585489356389657975562598571D+00, &
   -1.5938858604721398261388419455550D+00, &
   -1.1918269983500464260821358649242D+00, &
   -0.79287697691530893968593032998830D+00, &
   -0.39594273647142311094670041663436D+00, &
    0.0000000000000000000000000000000D+00, &
    0.39594273647142311094670041663436D+00, &
    0.79287697691530893968593032998830D+00, &
    1.1918269983500464260821358649242D+00, &
    1.5938858604721398261388419455550D+00, &
    2.0002585489356389657975562598571D+00, &
    2.4123177054804201051740184582119D+00, &
    2.8316804533902054557015640151425D+00, &
    3.2603207323135408104645401509648D+00, &
    3.7007434032314694224497164589673D+00, &
    4.1562717558181451724831352315314D+00, &
    4.6315595063128599420667997654336D+00, &
    5.1335955771123807045862968913996D+00, &
    5.6739614446185883296332558789276D+00, &
    6.2750787049428601427036567812530D+00, &
    6.9956801237185402753248521473232D+00, &
   -10.435499877854168053468115427285D+00, &
   -9.8028759912974963635223935286507D+00, &
   -9.2792019543050391319404745506496D+00, &
   -8.8118581437284546442526628275570D+00, &
   -8.3807683451863219343010651043788D+00, &
   -7.9755950801420373181541806298501D+00, &
   -7.5901395198641066762479783194468D+00, &
   -7.2203167078889678461161324222529D+00, &
   -6.8632544331795368527353285876066D+00, &
   -6.5168348106821160605273395854042D+00, &
   -6.1794379922705969862418461787263D+00, &
   -5.8497884000810673462526582961482D+00, &
   -5.5268572526403031425047575122840D+00, &
   -5.2097979830408354861575136416263D+00, &
   -4.8979018644975742350745099214868D+00, &
   -4.5905665744435190229271294569091D+00, &
   -4.2872733352824404031727616199454D+00, &
   -3.9875699104197157485227052068068D+00, &
   -3.6910577000963465117322810559754D+00, &
   -3.3973817713303911852755941806287D+00, &
   -3.1062230279282566329138616746036D+00, &
   -2.8172919672837977750747135657355D+00, &
   -2.5303236304712010926855221718499D+00, &
   -2.2450734604812066298995918179330D+00, &
   -1.9613138583081485293922008411321D+00, &
   -1.6788312791720137520802800622638D+00, &
   -1.3974237486049625107570752063702D+00, &
   -1.1168987050996462690510970277840D+00, &
   -0.83707109558947615977737795461293D+00, &
   -0.55776166427908221668763665253822D+00, &
   -0.27879538567115223986687628627202D+00, &
    0.00000000000000000000000000000000D+00, &
    0.27879538567115223986687628627202D+00, &
    0.55776166427908221668763665253822D+00, &
    0.83707109558947615977737795461293D+00, &
    1.1168987050996462690510970277840D+00, &
    1.3974237486049625107570752063702D+00, &
    1.6788312791720137520802800622638D+00, &
    1.9613138583081485293922008411321D+00, &
    2.2450734604812066298995918179330D+00, &
    2.5303236304712010926855221718499D+00, &
    2.8172919672837977750747135657355D+00, &
    3.1062230279282566329138616746036D+00, &
    3.3973817713303911852755941806287D+00, &
    3.6910577000963465117322810559754D+00, &
    3.9875699104197157485227052068068D+00, &
    4.2872733352824404031727616199454D+00, &
    4.5905665744435190229271294569091D+00, &
    4.8979018644975742350745099214868D+00, &
    5.2097979830408354861575136416263D+00, &
    5.5268572526403031425047575122840D+00, &
    5.8497884000810673462526582961482D+00, &
    6.1794379922705969862418461787263D+00, &
    6.5168348106821160605273395854042D+00, &
    6.8632544331795368527353285876066D+00, &
    7.2203167078889678461161324222529D+00, &
    7.5901395198641066762479783194468D+00, &
    7.9755950801420373181541806298501D+00, &
    8.3807683451863219343010651043788D+00, &
    8.8118581437284546442526628275570D+00, &
    9.2792019543050391319404745506496D+00, &
    9.8028759912974963635223935286507D+00, &
    10.435499877854168053468115427285D+00, &
   -15.228338148167350978246954433464D+00, &
   -14.669595158833972632746354112896D+00, &
   -14.209085995284870755168244250887D+00, &
   -13.799722290211676634645246746673D+00, &
   -13.423518590070950062438258321855D+00, &
   -13.071208660474601901583995439649D+00, &
   -12.737235652415686338138003924072D+00, &
   -12.417939378869715805445879624069D+00, &
   -12.110749020947747600132123508132D+00, &
   -11.813772198267727195134584136191D+00, &
   -11.525565112572696599167888588564D+00, &
   -11.244994583785543445194384194300D+00, &
   -10.971150569840247423423040263881D+00, &
   -10.703288201027481347670940744690D+00, &
   -10.440787957772772867742591798027D+00, &
   -10.183127473450343888624126450357D+00, &
   -9.9298610495114250736847004273684D+00, &
   -9.6806044412474728038150712732737D+00, &
   -9.4350233389881650135019598506287D+00, &
   -9.1928244988460305715774195052527D+00, &
   -8.9537488108565404323807890169970D+00, &
   -8.7175658087076307363833999548548D+00, &
   -8.4840692689832473326097180339984D+00, &
   -8.2530736454457156579694124243888D+00, &
   -8.0244111514703375578594739796798D+00, &
   -7.7979293513870105420829120455591D+00, &
   -7.5734891556083454022834960763301D+00, &
   -7.3509631392269052701961258043733D+00, &
   -7.1302341220350710668064025713431D+00, &
   -6.9111939615465713197465633109366D+00, &
   -6.6937425208758294190074417381666D+00, &
   -6.4777867811645365448144903821487D+00, &
   -6.2632400742737354345609723857092D+00, &
   -6.0500214161419845694465474482388D+00, &
   -5.8380549248774187386601690807757D+00, &
   -5.6272693105464816659423455794909D+00, &
   -5.4175974259243240722848425872924D+00, &
   -5.2089758693153983587570258372239D+00, &
   -5.0013446320386360038520809107373D+00, &
   -4.7946467843764925009748509930857D+00, &
   -4.5888281947698372951606485031212D+00, &
   -4.3838372778464736294253744407459D+00, &
   -4.1796247675352031349421189892408D+00, &
   -3.9761435120673355916035814195920D+00, &
   -3.7733482881250526721004678400057D+00, &
   -3.5711956317782180447199756485249D+00, &
   -3.3696436841717397896643629240035D+00, &
   -3.1686520501953630191857798261495D+00, &
   -2.9681816685955910267761649521505D+00, &
   -2.7681946921824058801226545958892D+00, &
   -2.5686543769473501723144013022363D+00, &
   -2.3695249790490401080012474645702D+00, &
   -2.1707716587411506879498498083695D+00, &
   -1.9723603904195020079324743227565D+00, &
   -1.7742578780516791584676442103681D+00, &
   -1.5764314753267801315519597621879D+00, &
   -1.3788491099261778091441557053728D+00, &
   -1.1814792113700685848678583598423D+00, &
   -0.98429064194027277726568984213773D+00, &
   -0.78725263021825034151596831878971D+00, &
   -0.59033470680942102142230439346102D+00, &
   -0.39350664185130136568037826200185D+00, &
   -0.19673838392423251964272239737078D+00, &
    0.0000000000000000000000000000000D+00, &
    0.19673838392423251964272239737078D+00, &
    0.39350664185130136568037826200185D+00, &
    0.59033470680942102142230439346102D+00, &
    0.78725263021825034151596831878971D+00, &
    0.98429064194027277726568984213773D+00, &
    1.1814792113700685848678583598423D+00, &
    1.3788491099261778091441557053728D+00, &
    1.5764314753267801315519597621879D+00, &
    1.7742578780516791584676442103681D+00, &
    1.9723603904195020079324743227565D+00, &
    2.1707716587411506879498498083695D+00, &
    2.3695249790490401080012474645702D+00, &
    2.5686543769473501723144013022363D+00, &
    2.7681946921824058801226545958892D+00, &
    2.9681816685955910267761649521505D+00, &
    3.1686520501953630191857798261495D+00, &
    3.3696436841717397896643629240035D+00, &
    3.5711956317782180447199756485249D+00, &
    3.7733482881250526721004678400057D+00, &
    3.9761435120673355916035814195920D+00, &
    4.1796247675352031349421189892408D+00, &
    4.3838372778464736294253744407459D+00, &
    4.5888281947698372951606485031212D+00, &
    4.7946467843764925009748509930857D+00, &
    5.0013446320386360038520809107373D+00, &
    5.2089758693153983587570258372239D+00, &
    5.4175974259243240722848425872924D+00, &
    5.6272693105464816659423455794909D+00, &
    5.8380549248774187386601690807757D+00, &
    6.0500214161419845694465474482388D+00, &
    6.2632400742737354345609723857092D+00, &
    6.4777867811645365448144903821487D+00, &
    6.6937425208758294190074417381666D+00, &
    6.9111939615465713197465633109366D+00, &
    7.1302341220350710668064025713431D+00, &
    7.3509631392269052701961258043733D+00, &
    7.5734891556083454022834960763301D+00, &
    7.7979293513870105420829120455591D+00, &
    8.0244111514703375578594739796798D+00, &
    8.2530736454457156579694124243888D+00, &
    8.4840692689832473326097180339984D+00, &
    8.7175658087076307363833999548548D+00, &
    8.9537488108565404323807890169970D+00, &
    9.1928244988460305715774195052527D+00, &
    9.4350233389881650135019598506287D+00, &
    9.6806044412474728038150712732737D+00, &
    9.9298610495114250736847004273684D+00, &
    10.183127473450343888624126450357D+00, &
    10.440787957772772867742591798027D+00, &
    10.703288201027481347670940744690D+00, &
    10.971150569840247423423040263881D+00, &
    11.244994583785543445194384194300D+00, &
    11.525565112572696599167888588564D+00, &
    11.813772198267727195134584136191D+00, &
    12.110749020947747600132123508132D+00, &
    12.417939378869715805445879624069D+00, &
    12.737235652415686338138003924072D+00, &
    13.071208660474601901583995439649D+00, &
    13.423518590070950062438258321855D+00, &
    13.799722290211676634645246746673D+00, &
    14.209085995284870755168244250887D+00, &
    14.669595158833972632746354112896D+00, &
    15.228338148167350978246954433464D+00  &
         /)

  if ( any ( grid_base(1:dim_num) < 0 ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GH_ABSCISSA - Fatal error!'
    write ( *, '(a)' ) '  Some base values are less than 0.'
    stop
  end if

  if ( any ( 63 < grid_base(1:dim_num) ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GH_ABSCISSA - Fatal error!'
    write ( *, '(a)' ) '  Some base values are greater than 63.'
    stop
  end if

  do point = 1, point_num
    do dim = 1, dim_num

      level = i4_log_2 ( grid_base(dim) + 1 )

      pointer = skip(level) + ( grid_index(dim,point) + grid_base(dim) + 1 )

      grid_point(dim,point) = x(pointer)

    end do
  end do

  return
end
subroutine gh_weights ( order, weight )

!*****************************************************************************80
!
!! GH_WEIGHTS returns weights for certain Gauss Hermite quadrature rules.
!
!  Discussion:
!
!    The allowed orders are 1, 3, 7, 15, 31, 63, and 127.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ORDER, the order of the rule.
!    ORDER must be 1, 3, 7, 15, 31, 63 or 127.
!
!    Output, real ( kind = 8 ) WEIGHT(ORDER), the weights.
!    The weights are positive, symmetric and should sum to SQRT(PI).
!
  implicit none

  integer ( kind = 4 ) order

  real ( kind = 8 ) weight(order)

  if ( order == 1 ) then

    weight(1) = 1.77245385090551602729816748334D+00

  else if ( order == 3 ) then

    weight(1) = 0.295408975150919337883027913890D+00
    weight(2) = 0.118163590060367735153211165556D+01
    weight(3) = 0.295408975150919337883027913890D+00

  else if ( order == 7 ) then

    weight(1) = 0.971781245099519154149424255939D-03
    weight(2) = 0.545155828191270305921785688417D-01
    weight(3) = 0.425607252610127800520317466666D+00
    weight(4) = 0.810264617556807326764876563813D+00
    weight(5) = 0.425607252610127800520317466666D+00
    weight(6) = 0.545155828191270305921785688417D-01
    weight(7) = 0.971781245099519154149424255939D-03

  else if ( order == 15 ) then

    weight(1) =  0.152247580425351702016062666965D-08
    weight(2) =  0.105911554771106663577520791055D-05
    weight(3) =  0.100004441232499868127296736177D-03
    weight(4) =  0.277806884291277589607887049229D-02
    weight(5) =  0.307800338725460822286814158758D-01
    weight(6) =  0.158488915795935746883839384960D+00
    weight(7) =  0.412028687498898627025891079568D+00
    weight(8) =  0.564100308726417532852625797340D+00
    weight(9) =  0.412028687498898627025891079568D+00
    weight(10) = 0.158488915795935746883839384960D+00
    weight(11) = 0.307800338725460822286814158758D-01
    weight(12) = 0.277806884291277589607887049229D-02
    weight(13) = 0.100004441232499868127296736177D-03
    weight(14) = 0.105911554771106663577520791055D-05
    weight(15) = 0.152247580425351702016062666965D-08

  else if ( order == 31 ) then

    weight(  1) =   0.46189683944498305857470556847735D-21
    weight(  2) =   0.51106090079112519643027197715274D-17
    weight(  3) =   0.58995564987355133075257722133966D-14
    weight(  4) =   0.18603735214463569590294465062239D-11
    weight(  5) =   0.23524920032013205739850619940094D-09
    weight(  6) =   0.14611988344865057576066495091513D-07
    weight(  7) =   0.50437125589241034841778074689627D-06
    weight(  8) =   0.10498602757642934202945441341697D-04
    weight(  9) =   0.13952090395003623854995664958146D-03
    weight( 10) =   0.12336833073030489880608311394968D-02
    weight( 11) =   0.74827999140119116765002499116934D-02
    weight( 12) =   0.31847230731201222775249585776902D-01
    weight( 13) =   0.96717948160569462991143316029341D-01
    weight( 14) =   0.21213278866810461318136114862419D+00
    weight( 15) =   0.33877265789305344906000174083214D+00
    weight( 16) =   0.39577855609737786462923720809676D+00
    weight( 17) =   0.33877265789305344906000174083214D+00
    weight( 18) =   0.21213278866810461318136114862419D+00
    weight( 19) =   0.96717948160569462991143316029341D-01
    weight( 20) =   0.31847230731201222775249585776902D-01
    weight( 21) =   0.74827999140119116765002499116934D-02
    weight( 22) =   0.12336833073030489880608311394968D-02
    weight( 23) =   0.13952090395003623854995664958146D-03
    weight( 24) =   0.10498602757642934202945441341697D-04
    weight( 25) =   0.50437125589241034841778074689627E-06
    weight( 26) =   0.14611988344865057576066495091513D-07
    weight( 27) =   0.23524920032013205739850619940094D-09
    weight( 28) =   0.18603735214463569590294465062239D-11
    weight( 29) =   0.58995564987355133075257722133966D-14
    weight( 30) =   0.51106090079112519643027197715274D-17
    weight( 31) =   0.46189683944498305857470556847735D-21

  else if ( order == 63 ) then

    weight(  1) =   0.37099206434787551197827130470031D-47
    weight(  2) =   0.10400778615192299534481914814892D-41
    weight(  3) =   0.19796804708258311251124226474396D-37
    weight(  4) =   0.84687478191640015120141181138947D-34
    weight(  5) =   0.13071305930779945903630127634063D-30
    weight(  6) =   0.93437837175367456929765381518998D-28
    weight(  7) =   0.36027426635173044862245783257252D-25
    weight(  8) =   0.82963863115951789374753323156164D-23
    weight(  9) =   0.12266629909105281472971700203949D-20
    weight( 10) =   0.12288435628797061539461585325494D-18
    weight( 11) =   0.86925536958188009075932426691516D-17
    weight( 12) =   0.44857058689176221240330804981619D-15
    weight( 13) =   0.17335817955735154599902643794700D-13
    weight( 14) =   0.51265062385038307838565047455223D-12
    weight( 15) =   0.11808921844532942490513037158404D-10
    weight( 16) =   0.21508698297808025739828859845140D-09
    weight( 17) =   0.31371929535285447801497640621672D-08
    weight( 18) =   0.37041625984781705796752840204084D-07
    weight( 19) =   0.35734732949879669663960738150956D-06
    weight( 20) =   0.28393114498380927832990899215541D-05
    weight( 21) =   0.18709113003730498008961134765721D-04
    weight( 22) =   0.10284880800653635546698378640623D-03
    weight( 23) =   0.47411702610173128107201781718693D-03
    weight( 24) =   0.18409222622384813438539657470055D-02
    weight( 25) =   0.60436044551187631655712178246467D-02
    weight( 26) =   0.16829299199599730926458559757600D-01
    weight( 27) =   0.39858264027692992170237391875317D-01
    weight( 28) =   0.80467087993950415219587554532823D-01
    weight( 29) =   0.13871950817615293377792092082674D+00
    weight( 30) =   0.20448695346833761570957197160475D+00
    weight( 31) =   0.25799889943058042204920467417642D+00
    weight( 32) =   0.27876694884838411919175686949858D+00
    weight( 33) =   0.25799889943058042204920467417642D+00
    weight( 34) =   0.20448695346833761570957197160475D+00
    weight( 35) =   0.13871950817615293377792092082674D+00
    weight( 36) =   0.80467087993950415219587554532823D-01
    weight( 37) =   0.39858264027692992170237391875317D-01
    weight( 38) =   0.16829299199599730926458559757600D-01
    weight( 39) =   0.60436044551187631655712178246467D-02
    weight( 40) =   0.18409222622384813438539657470055D-02
    weight( 41) =   0.47411702610173128107201781718693D-03
    weight( 42) =   0.10284880800653635546698378640623D-03
    weight( 43) =   0.18709113003730498008961134765721D-04
    weight( 44) =   0.28393114498380927832990899215541D-05
    weight( 45) =   0.35734732949879669663960738150956D-06
    weight( 46) =   0.37041625984781705796752840204084D-07
    weight( 47) =   0.31371929535285447801497640621672D-08
    weight( 48) =   0.21508698297808025739828859845140D-09
    weight( 49) =   0.11808921844532942490513037158404D-10
    weight( 50) =   0.51265062385038307838565047455223D-12
    weight( 51) =   0.17335817955735154599902643794700D-13
    weight( 52) =   0.44857058689176221240330804981619D-15
    weight( 53) =   0.86925536958188009075932426691516D-17
    weight( 54) =   0.12288435628797061539461585325494D-18
    weight( 55) =   0.12266629909105281472971700203949D-20
    weight( 56) =   0.82963863115951789374753323156164D-23
    weight( 57) =   0.36027426635173044862245783257252D-25
    weight( 58) =   0.93437837175367456929765381518998D-28
    weight( 59) =   0.13071305930779945903630127634063E-30
    weight( 60) =   0.84687478191640015120141181138947D-34
    weight( 61) =   0.19796804708258311251124226474396D-37
    weight( 62) =   0.10400778615192299534481914814892D-41
    weight( 63) =   0.37099206434787551197827130470031D-47

  else if ( order == 127 ) then

    weight(  1) =   0.12504497577050595552677230002883D-100
    weight(  2) =   0.17272798059419131415318615789672D-93
    weight(  3) =   0.89321681571986548608031150791499D-88
    weight(  4) =   0.77306185240893578449625186483810D-83
    weight(  5) =   0.20143957652648255497735460506196D-78
    weight(  6) =   0.21503714733610239701351039429345D-74
    weight(  7) =   0.11341924208594594813715533569504D-70
    weight(  8) =   0.33489139011795051950683388483136D-67
    weight(  9) =   0.60486548964016681064424451668405D-64
    weight( 10) =   0.71375092946352177824971347343892D-61
    weight( 11) =   0.57884563374885556636801095624030D-58
    weight( 12) =   0.33581166223858230300409326551248D-55
    weight( 13) =   0.14394641949253923568603163698953D-52
    weight( 14) =   0.46821808383216117724080263903889D-50
    weight( 15) =   0.11817054440684264071348471955361D-47
    weight( 16) =   0.23581659156008927203181682045005D-45
    weight( 17) =   0.37814427940797540210712758405540D-43
    weight( 18) =   0.49411031115771638145610738414006D-41
    weight( 19) =   0.53255303775425059266087298458297D-39
    weight( 20) =   0.47854390680131484999315199332765D-37
    weight( 21) =   0.36191883445952356128627543209554D-35
    weight( 22) =   0.23232083386343554805352497446119D-33
    weight( 23) =   0.12753331411008716683688974281454D-31
    weight( 24) =   0.60277753850758742112436095241270D-30
    weight( 25) =   0.24679773241777200207460855084439D-28
    weight( 26) =   0.88019567691698482573264198727415D-27
    weight( 27) =   0.27482489212040561315005725890593D-25
    weight( 28) =   0.75468218903085486125222816438456D-24
    weight( 29) =   0.18303134636280466270545996891835D-22
    weight( 30) =   0.39355990860860813085582448449811D-21
    weight( 31) =   0.75293161638581191068419292570042D-20
    weight( 32) =   0.12857997786722855037584105682618D-18
    weight( 33) =   0.19659326888445857792541925311450D-17
    weight( 34) =   0.26986511907214101894995783364250D-16
    weight( 35) =   0.33344414303198856330118301113874D-15
    weight( 36) =   0.37173303125150639885726463109574D-14
    weight( 37) =   0.37473954472839737091885387788983D-13
    weight( 38) =   0.34230094493397259538669512076007D-12
    weight( 39) =   0.28385303724993373166810860630552D-11
    weight( 40) =   0.21406920290454669208938772802828D-10
    weight( 41) =   0.14706331273431716244229273183839D-09
    weight( 42) =   0.92173940967434659264335883218167D-09
    weight( 43) =   0.52781663936972714041837056042506D-08
    weight( 44) =   0.27650497044951117835905283127679D-07
    weight( 45) =   0.13267855842539464770913063113371D-06
    weight( 46) =   0.58380944276113062188573331195042D-06
    weight( 47) =   0.23581561724775629112332165335800D-05
    weight( 48) =   0.87524468034280444703919485644809D-05
    weight( 49) =   0.29876790535909012274846532159647D-04
    weight( 50) =   0.93874435720072545206729594267039D-04
    weight( 51) =   0.27170762627931172053444716883938D-03
    weight( 52) =   0.72493929742498358979684249380921D-03
    weight( 53) =   0.17841208326763432884316727108264D-02
    weight( 54) =   0.40524855186046131499765636276283D-02
    weight( 55) =   0.85000263041544110385806705526917D-02
    weight( 56) =   0.16471142241609687824005585301760D-01
    weight( 57) =   0.29499296248213632269675010319119D-01
    weight( 58) =   0.48847387114300011006959603975676D-01
    weight( 59) =   0.74807989768583731416517226905270D-01
    weight( 60) =   0.10598520508090929403834368934301D+00
    weight( 61) =   0.13893945309051540832066283010510D+00
    weight( 62) =   0.16856236074207929740526975049765D+00
    weight( 63) =   0.18927849580120432177170145550076D+00
    weight( 64) =   0.19673340688823289786163676995151D+00
    weight( 65) =   0.18927849580120432177170145550076D+00
    weight( 66) =   0.16856236074207929740526975049765D+00
    weight( 67) =   0.13893945309051540832066283010510D+00
    weight( 68) =   0.10598520508090929403834368934301D+00
    weight( 69) =   0.74807989768583731416517226905270D-01
    weight( 70) =   0.48847387114300011006959603975676D-01
    weight( 71) =   0.29499296248213632269675010319119D-01
    weight( 72) =   0.16471142241609687824005585301760D-01
    weight( 73) =   0.85000263041544110385806705526917D-02
    weight( 74) =   0.40524855186046131499765636276283D-02
    weight( 75) =   0.17841208326763432884316727108264D-02
    weight( 76) =   0.72493929742498358979684249380921D-03
    weight( 77) =   0.27170762627931172053444716883938D-03
    weight( 78) =   0.93874435720072545206729594267039D-04
    weight( 79) =   0.29876790535909012274846532159647D-04
    weight( 80) =   0.87524468034280444703919485644809D-05
    weight( 81) =   0.23581561724775629112332165335800D-05
    weight( 82) =   0.58380944276113062188573331195042D-06
    weight( 83) =   0.13267855842539464770913063113371D-06
    weight( 84) =   0.27650497044951117835905283127679D-07
    weight( 85) =   0.52781663936972714041837056042506D-08
    weight( 86) =   0.92173940967434659264335883218167D-09
    weight( 87) =   0.14706331273431716244229273183839D-09
    weight( 88) =   0.21406920290454669208938772802828D-10
    weight( 89) =   0.28385303724993373166810860630552D-11
    weight( 90) =   0.34230094493397259538669512076007D-12
    weight( 91) =   0.37473954472839737091885387788983D-13
    weight( 92) =   0.37173303125150639885726463109574D-14
    weight( 93) =   0.33344414303198856330118301113874D-15
    weight( 94) =   0.26986511907214101894995783364250D-16
    weight( 95) =   0.19659326888445857792541925311450D-17
    weight( 96) =   0.12857997786722855037584105682618D-18
    weight( 97) =   0.75293161638581191068419292570042D-20
    weight( 98) =   0.39355990860860813085582448449811D-21
    weight( 99) =   0.18303134636280466270545996891835D-22
    weight(100) =   0.75468218903085486125222816438456D-24
    weight(101) =   0.27482489212040561315005725890593D-25
    weight(102) =   0.88019567691698482573264198727415D-27
    weight(103) =   0.24679773241777200207460855084439D-28
    weight(104) =   0.60277753850758742112436095241270D-30
    weight(105) =   0.12753331411008716683688974281454D-31
    weight(106) =   0.23232083386343554805352497446119D-33
    weight(107) =   0.36191883445952356128627543209554D-35
    weight(108) =   0.47854390680131484999315199332765D-37
    weight(109) =   0.53255303775425059266087298458297D-39
    weight(110) =   0.49411031115771638145610738414006D-41
    weight(111) =   0.37814427940797540210712758405540D-43
    weight(112) =   0.23581659156008927203181682045005D-45
    weight(113) =   0.11817054440684264071348471955361D-47
    weight(114) =   0.46821808383216117724080263903889D-50
    weight(115) =   0.14394641949253923568603163698953D-52
    weight(116) =   0.33581166223858230300409326551248D-55
    weight(117) =   0.57884563374885556636801095624030D-58
    weight(118) =   0.71375092946352177824971347343892D-61
    weight(119) =   0.60486548964016681064424451668405D-64
    weight(120) =   0.33489139011795051950683388483136D-67
    weight(121) =   0.11341924208594594813715533569504D-70
    weight(122) =   0.21503714733610239701351039429345D-74
    weight(123) =   0.20143957652648255497735460506196D-78
    weight(124) =   0.77306185240893578449625186483810D-83
    weight(125) =   0.89321681571986548608031150791499D-88
    weight(126) =   0.17272798059419131415318615789672D-93
    weight(127) =   0.12504497577050595552677230002883D-100

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GH_WEIGHTS - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of ORDER = ', order
    write ( *, '(a)' ) '  Legal values are 1, 3, 7, 15, 31, 63 and 127.'
    stop

  end if

  return
end
subroutine gl_abscissa ( dim_num, point_num, grid_index, grid_base, grid_point )

!*****************************************************************************80
!
!! GL_ABSCISSA sets abscissas for multidimensional Gauss Legendre quadrature.
!
!  Discussion:
!
!    The "nesting" as it occurs for Gauss Legendre sparse grids simply
!    involves the use of a specified set of permissible orders for the
!    rule.
!
!    The X array lists the (complete) Gauss Legendre abscissas for rules
!    of order 1, 3, 7, 15, 31, 63, and 127 in order.
!
!    The FORTRAN90 standard does not require more than 99 continuation lines,
!    and the blinking IBM FORTRAN90 compiler enforces that limitation.
!    Hence, instead of a data statement assignment for X, we are reduced
!    to executable assignments to ranges of X.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, integer ( kind = 4 ) GRID_INDEX(DIM_NUM,POINT_NUM), for each
!    point and dimension, the index of the Gauss Legendre abscissa.
!
!    Input, integer ( kind = 4 ) GRID_BASE(DIM_NUM), the "base" of the
!    Gauss Legendre rule being used in each dimension.
!
!    Output, real ( kind = 8 ) GRID_POINT(DIM_NUM), the grid points of
!    Gauss Legendre abscissas.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) grid_base(dim_num)
  integer ( kind = 4 ) grid_index(dim_num,point_num)
  real ( kind = 8 ) grid_point(dim_num,point_num)
  integer ( kind = 4 ) i4_log_2
  integer ( kind = 4 ) level
  integer ( kind = 4 ) point
  integer ( kind = 4 ) pointer
  integer ( kind = 4 ), dimension ( 0:7 ) :: skip = (/ &
    0, 1, 4, 11, 26, 57, 120, 247 /)
  real ( kind = 8 ), dimension ( 247 ) :: x

  x(1:100) = (/ &
       0.0D+00, &
     - 0.774596669241483377035853079956D+00, &
       0.0D+00, &
       0.774596669241483377035853079956D+00, &
     - 0.949107912342758524526189684048D+00, &
     - 0.741531185599394439863864773281D+00, &
     - 0.405845151377397166906606412077D+00, &
       0.0D+00, &
       0.405845151377397166906606412077D+00, &
       0.741531185599394439863864773281D+00, &
       0.949107912342758524526189684048D+00, &
      - 0.987992518020485428489565718587D+00, &
      - 0.937273392400705904307758947710D+00, &
      - 0.848206583410427216200648320774D+00, &
      - 0.724417731360170047416186054614D+00, &
      - 0.570972172608538847537226737254D+00, &
      - 0.394151347077563369897207370981D+00, &
      - 0.201194093997434522300628303395D+00, &
        0.0D+00, &
        0.201194093997434522300628303395D+00, &
       0.394151347077563369897207370981D+00, &
       0.570972172608538847537226737254D+00, &
       0.724417731360170047416186054614D+00, &
       0.848206583410427216200648320774D+00, &
       0.937273392400705904307758947710D+00, &
       0.987992518020485428489565718587D+00, &
      -0.99708748181947707454263838179654D+00, &
      -0.98468590966515248400211329970113D+00, &
      -0.96250392509294966178905249675943D+00, &
      -0.93075699789664816495694576311725D+00, &
      -0.88976002994827104337419200908023D+00, &
      -0.83992032014626734008690453594388D+00, &
      -0.78173314841662494040636002019484D+00, &
      -0.71577678458685328390597086536649D+00, &
      -0.64270672292426034618441820323250D+00, &
      -0.56324916140714926272094492359516D+00, &
      -0.47819378204490248044059403935649D+00, &
      -0.38838590160823294306135146128752D+00, &
      -0.29471806998170161661790389767170D+00, &
      -0.19812119933557062877241299603283D+00, &
      -0.99555312152341520325174790118941D-01, &
       0.00000000000000000000000000000000D+00, &
       0.99555312152341520325174790118941D-01, &
       0.19812119933557062877241299603283D+00, &
       0.29471806998170161661790389767170D+00, &
       0.38838590160823294306135146128752D+00, &
       0.47819378204490248044059403935649D+00, &
       0.56324916140714926272094492359516D+00, &
       0.64270672292426034618441820323250D+00, &
       0.71577678458685328390597086536649D+00, &
       0.78173314841662494040636002019484D+00, &
       0.83992032014626734008690453594388D+00, &
       0.88976002994827104337419200908023D+00, &
       0.93075699789664816495694576311725D+00, &
       0.96250392509294966178905249675943D+00, &
       0.98468590966515248400211329970113D+00, &
       0.99708748181947707454263838179654D+00, &
      -0.99928298402912378050701628988630D+00, &
      -0.99622401277797010860209018267357D+00, &
      -0.99072854689218946681089469460884D+00, &
      -0.98280881059372723486251140727639D+00, &
      -0.97248403469757002280196067864927D+00, &
      -0.95977944975894192707035416626398D+00, &
      -0.94472613404100980296637531962798D+00, &
      -0.92736092062184320544703138132518D+00, &
      -0.90772630277853155803695313291596D+00, &
      -0.88587032850785342629029845731337D+00, &
      -0.86184648236412371953961183943106D+00, &
      -0.83571355431950284347180776961571D+00, &
      -0.80753549577345676005146598636324D+00, &
      -0.77738126299037233556333018991104D+00, &
      -0.74532464831784741782932166103759D+00, &
      -0.71144409958484580785143153770401D+00, &
      -0.67582252811498609013110331596954D+00, &
      -0.63854710582136538500030695387338D+00, &
      -0.59970905187762523573900892686880D+00, &
      -0.55940340948628501326769780007005D+00, &
      -0.51772881329003324812447758452632D+00, &
      -0.47478724799480439992221230985149D+00, &
      -0.43068379879511160066208893391863D+00, &
      -0.38552639421224789247761502227440D+00, &
      -0.33942554197458440246883443159432D+00, &
      -0.29249405858625144003615715555067D+00, &
      -0.24484679324595336274840459392483D+00, &
      -0.19660034679150668455762745706572D+00, &
      -0.14787278635787196856983909655297D+00, &
      -0.98783356446945279529703669453922D-01, &
      -0.49452187116159627234233818051808D-01, &
       0.00000000000000000000000000000000D+00, &
       0.49452187116159627234233818051808D-01, &
       0.98783356446945279529703669453922D-01, &
       0.14787278635787196856983909655297D+00, &
       0.19660034679150668455762745706572D+00, &
       0.24484679324595336274840459392483D+00, &
       0.29249405858625144003615715555067D+00, &
       0.33942554197458440246883443159432D+00, &
       0.38552639421224789247761502227440D+00, &
       0.43068379879511160066208893391863D+00, &
       0.47478724799480439992221230985149D+00, &
       0.51772881329003324812447758452632D+00 /)

  x(101:200) = (/ &
       0.55940340948628501326769780007005D+00, &
       0.59970905187762523573900892686880D+00, &
       0.63854710582136538500030695387338D+00, &
       0.67582252811498609013110331596954D+00, &
       0.71144409958484580785143153770401D+00, &
       0.74532464831784741782932166103759D+00, &
       0.77738126299037233556333018991104D+00, &
       0.80753549577345676005146598636324D+00, &
       0.83571355431950284347180776961571D+00, &
       0.86184648236412371953961183943106D+00, &
       0.88587032850785342629029845731337D+00, &
       0.90772630277853155803695313291596D+00, &
       0.92736092062184320544703138132518D+00, &
       0.94472613404100980296637531962798D+00, &
       0.95977944975894192707035416626398D+00, &
       0.97248403469757002280196067864927D+00, &
       0.98280881059372723486251140727639D+00, &
       0.99072854689218946681089469460884D+00, &
       0.99622401277797010860209018267357D+00, &
       0.99928298402912378050701628988630D+00, &
      -0.99982213041530614629963254927125D+00, &
      -0.99906293435531189513828920479421D+00, &
      -0.99769756618980462107441703193392D+00, &
      -0.99572655135202722663543337085008D+00, &
      -0.99315104925451714736113079489080D+00, &
      -0.98997261459148415760778669967548D+00, &
      -0.98619317401693166671043833175407D+00, &
      -0.98181502080381411003346312451200D+00, &
      -0.97684081234307032681744391886221D+00, &
      -0.97127356816152919228894689830512D+00, &
      -0.96511666794529212109082507703391D+00, &
      -0.95837384942523877114910286998060D+00, &
      -0.95104920607788031054790764659636D+00, &
      -0.94314718462481482734544963026201D+00, &
      -0.93467258232473796857363487794906D+00, &
      -0.92563054405623384912746466814259D+00, &
      -0.91602655919146580931308861741716D+00, &
      -0.90586645826182138280246131760282D+00, &
      -0.89515640941708370896904382642451D+00, &
      -0.88390291468002656994525794802849D+00, &
      -0.87211280599856071141963753428864D+00, &
      -0.85979324109774080981203134414483D+00, &
      -0.84695169913409759845333931085437D+00, &
      -0.83359597615489951437955716480123D+00, &
      -0.81973418036507867415511910167470D+00, &
      -0.80537472720468021466656079404644D+00, &
      -0.79052633423981379994544995252740D+00, &
      -0.77519801587020238244496276354566D+00, &
      -0.75939907785653667155666366659810D+00, &
      -0.74313911167095451292056688997595D+00, &
      -0.72642798867407268553569290153270D+00, &
      -0.70927585412210456099944463906757D+00, &
      -0.69169312100770067015644143286666D+00, &
      -0.67369046373825048534668253831602D+00, &
      -0.65527881165548263027676505156852D+00, &
      -0.63646934240029724134760815684175D+00, &
      -0.61727347512685828385763916340822D+00, &
      -0.59770286357006522938441201887478D+00, &
      -0.57776938897061258000325165713764D+00, &
      -0.55748515286193223292186190687872D+00, &
      -0.53686246972339756745816636353452D+00, &
      -0.51591385950424935727727729906662D+00, &
      -0.49465204002278211739494017368636D+00, &
      -0.47308991924540524164509989939699D+00, &
      -0.45124058745026622733189858020729D+00, &
      -0.42911730928019337626254405355418D+00, &
      -0.40673351568978256340867288124339D+00, &
      -0.38410279579151693577907781452239D+00, &
      -0.36123888860586970607092484346723D+00, &
      -0.33815567472039850137600027657095D+00, &
      -0.31486716786289498148601475374890D+00, &
      -0.29138750639370562079451875284568D+00, &
      -0.26773094472238862088834352027938D+00, &
      -0.24391184465391785797071324453138D+00, &
      -0.21994466666968754245452337866940D+00, &
      -0.19584396114861085150428162519610D+00, &
      -0.17162435953364216500834492248954D+00, &
      -0.14730056544908566938932929319807D+00, &
      -0.12288734577408297172603365288567D+00, &
      -0.98399521677698970751091751509101D-01, &
      -0.73851959621048545273440409360569D-01, &
      -0.49259562331926630315379321821927D-01, &
      -0.24637259757420944614897071846088D-01, &
       0.00000000000000000000000000000000D+00, &
       0.24637259757420944614897071846088D-01, &
       0.49259562331926630315379321821927D-01, &
       0.73851959621048545273440409360569D-01, &
       0.98399521677698970751091751509101D-01, &
       0.12288734577408297172603365288567D+00, &
       0.14730056544908566938932929319807D+00, &
       0.17162435953364216500834492248954D+00, &
       0.19584396114861085150428162519610D+00, &
       0.21994466666968754245452337866940D+00, &
       0.24391184465391785797071324453138D+00, &
       0.26773094472238862088834352027938D+00, &
       0.29138750639370562079451875284568D+00, &
       0.31486716786289498148601475374890D+00, &
       0.33815567472039850137600027657095D+00, &
       0.36123888860586970607092484346723D+00, &
       0.38410279579151693577907781452239D+00 /)

  x(201:247) = (/ &
       0.40673351568978256340867288124339D+00, &
       0.42911730928019337626254405355418D+00, &
       0.45124058745026622733189858020729D+00, &
       0.47308991924540524164509989939699D+00, &
       0.49465204002278211739494017368636D+00, &
       0.51591385950424935727727729906662D+00, &
       0.53686246972339756745816636353452D+00, &
       0.55748515286193223292186190687872D+00, &
       0.57776938897061258000325165713764D+00, &
       0.59770286357006522938441201887478D+00, &
       0.61727347512685828385763916340822D+00, &
       0.63646934240029724134760815684175D+00, &
       0.65527881165548263027676505156852D+00, &
       0.67369046373825048534668253831602D+00, &
       0.69169312100770067015644143286666D+00, &
       0.70927585412210456099944463906757D+00, &
       0.72642798867407268553569290153270D+00, &
       0.74313911167095451292056688997595D+00, &
       0.75939907785653667155666366659810D+00, &
       0.77519801587020238244496276354566D+00, &
       0.79052633423981379994544995252740D+00, &
       0.80537472720468021466656079404644D+00, &
       0.81973418036507867415511910167470D+00, &
       0.83359597615489951437955716480123D+00, &
       0.84695169913409759845333931085437D+00, &
       0.85979324109774080981203134414483D+00, &
       0.87211280599856071141963753428864D+00, &
       0.88390291468002656994525794802849D+00, &
       0.89515640941708370896904382642451D+00, &
       0.90586645826182138280246131760282D+00, &
       0.91602655919146580931308861741716D+00, &
       0.92563054405623384912746466814259D+00, &
       0.93467258232473796857363487794906D+00, &
       0.94314718462481482734544963026201D+00, &
       0.95104920607788031054790764659636D+00, &
       0.95837384942523877114910286998060D+00, &
       0.96511666794529212109082507703391D+00, &
       0.97127356816152919228894689830512D+00, &
       0.97684081234307032681744391886221D+00, &
       0.98181502080381411003346312451200D+00, &
       0.98619317401693166671043833175407D+00, &
       0.98997261459148415760778669967548D+00, &
       0.99315104925451714736113079489080D+00, &
       0.99572655135202722663543337085008D+00, &
       0.99769756618980462107441703193392D+00, &
       0.99906293435531189513828920479421D+00, &
       0.99982213041530614629963254927125D+00  &
         /)

  if ( any ( grid_base(1:dim_num) < 0 ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GL_ABSCISSA - Fatal error!'
    write ( *, '(a)' ) '  Some base values are less than 0.'
    stop
  end if

  if ( any ( 63 < grid_base(1:dim_num) ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GL_ABSCISSA - Fatal error!'
    write ( *, '(a)' ) '  Some base values are greater than 63.'
    stop
  end if

  do point = 1, point_num
    do dim = 1, dim_num

      level = i4_log_2 ( grid_base(dim) + 1 )

      pointer = skip(level) + ( grid_index(dim,point) + grid_base(dim) + 1 )

      grid_point(dim,point) = x(pointer)

    end do
  end do

  return
end
subroutine gl_weights ( order, weight )

!*****************************************************************************80
!
!! GL_WEIGHTS returns weights for certain Gauss Legendre quadrature rules.
!
!  Discussion:
!
!    The allowed orders are 1, 3, 7, 15, 31, 63 and 127.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ORDER, the order of the rule.
!    ORDER must be 1, 3, 7, 15, 31, 63 or 127.
!
!    Output, real ( kind = 8 ) WEIGHT(ORDER), the weights.
!    The weights are positive, symmetric and should sum to 2.
!
  implicit none

  integer ( kind = 4 ) order

  real ( kind = 8 ) weight(order)

  if ( order == 1 ) then

    weight(1) = 2.0D+00

  else if ( order == 3 ) then

    weight(1) = 5.0D+00 / 9.0D+00
    weight(2) = 8.0D+00 / 9.0D+00
    weight(3) = 5.0D+00 / 9.0D+00

  else if ( order == 7 ) then

    weight(1) = 0.129484966168869693270611432679D+00
    weight(2) = 0.279705391489276667901467771424D+00
    weight(3) = 0.381830050505118944950369775489D+00
    weight(4) = 0.417959183673469387755102040816D+00
    weight(5) = 0.381830050505118944950369775489D+00
    weight(6) = 0.279705391489276667901467771424D+00
    weight(7) = 0.129484966168869693270611432679D+00

  else if ( order == 15 ) then

    weight(1) =  0.307532419961172683546283935772D-01
    weight(2) =  0.703660474881081247092674164507D-01
    weight(3) =  0.107159220467171935011869546686D+00
    weight(4) =  0.139570677926154314447804794511D+00
    weight(5) =  0.166269205816993933553200860481D+00
    weight(6) =  0.186161000015562211026800561866D+00
    weight(7) =  0.198431485327111576456118326444D+00
    weight(8) =  0.202578241925561272880620199968D+00
    weight(9) =  0.198431485327111576456118326444D+00
    weight(10) = 0.186161000015562211026800561866D+00
    weight(11) = 0.166269205816993933553200860481D+00
    weight(12) = 0.139570677926154314447804794511D+00
    weight(13) = 0.107159220467171935011869546686D+00
    weight(14) = 0.703660474881081247092674164507D-01
    weight(15) = 0.307532419961172683546283935772D-01

  else if ( order == 31 ) then

    weight( 1) =   0.74708315792487746093913218970494D-02
    weight( 2) =   0.17318620790310582463552990782414D-01
    weight( 3) =   0.27009019184979421800608642617676D-01
    weight( 4) =   0.36432273912385464024392008749009D-01
    weight( 5) =   0.45493707527201102902315857856518D-01
    weight( 6) =   0.54103082424916853711666259085477D-01
    weight( 7) =   0.62174786561028426910343543686657D-01
    weight( 8) =   0.69628583235410366167756126255124D-01
    weight( 9) =   0.76390386598776616426357674901331D-01
    weight(10) =   0.82392991761589263903823367431962D-01
    weight(11) =   0.87576740608477876126198069695333D-01
    weight(12) =   0.91890113893641478215362871607150D-01
    weight(13) =   0.95290242912319512807204197487597D-01
    weight(14) =   0.97743335386328725093474010978997D-01
    weight(15) =   0.99225011226672307874875514428615D-01
    weight(16) =   0.99720544793426451427533833734349D-01
    weight(17) =   0.99225011226672307874875514428615D-01
    weight(18) =   0.97743335386328725093474010978997D-01
    weight(19) =   0.95290242912319512807204197487597D-01
    weight(20) =   0.91890113893641478215362871607150D-01
    weight(21) =   0.87576740608477876126198069695333D-01
    weight(22) =   0.82392991761589263903823367431962D-01
    weight(23) =   0.76390386598776616426357674901331D-01
    weight(24) =   0.69628583235410366167756126255124D-01
    weight(25) =   0.62174786561028426910343543686657D-01
    weight(26) =   0.54103082424916853711666259085477D-01
    weight(27) =   0.45493707527201102902315857856518D-01
    weight(28) =   0.36432273912385464024392008749009D-01
    weight(29) =   0.27009019184979421800608642617676D-01
    weight(30) =   0.17318620790310582463552990782414D-01
    weight(31) =   0.74708315792487746093913218970494D-02

  else if ( order == 63 ) then

    weight( 1) =   0.18398745955770837880499331680577D-02
    weight( 2) =   0.42785083468637618661951422543371D-02
    weight( 3) =   0.67102917659601362519069109850892D-02
    weight( 4) =   0.91259686763266563540586445877022D-02
    weight( 5) =   0.11519376076880041750750606118707D-01
    weight( 6) =   0.13884612616115610824866086365937D-01
    weight( 7) =   0.16215878410338338882283672974995D-01
    weight( 8) =   0.18507464460161270409260545805144D-01
    weight( 9) =   0.20753761258039090775341953421471D-01
    weight(10) =   0.22949271004889933148942319561770D-01
    weight(11) =   0.25088620553344986618630138068443D-01
    weight(12) =   0.27166574359097933225189839439413D-01
    weight(13) =   0.29178047208280526945551502154029D-01
    weight(14) =   0.31118116622219817508215988557189D-01
    weight(15) =   0.32982034883779341765683179672459D-01
    weight(16) =   0.34765240645355877697180504642788D-01
    weight(17) =   0.36463370085457289630452409787542D-01
    weight(18) =   0.38072267584349556763638324927889D-01
    weight(19) =   0.39587995891544093984807928149202D-01
    weight(20) =   0.41006845759666398635110037009072D-01
    weight(21) =   0.42325345020815822982505485403028D-01
    weight(22) =   0.43540267083027590798964315704401D-01
    weight(23) =   0.44648638825941395370332669516813D-01
    weight(24) =   0.45647747876292608685885992608542D-01
    weight(25) =   0.46535149245383696510395418746953D-01
    weight(26) =   0.47308671312268919080604988338844D-01
    weight(27) =   0.47966421137995131411052756195132D-01
    weight(28) =   0.48506789097883847864090099145802D-01
    weight(29) =   0.48928452820511989944709361549215D-01
    weight(30) =   0.49230380423747560785043116988145D-01
    weight(31) =   0.49411833039918178967039646116705D-01
    weight(32) =   0.49472366623931020888669360420926D-01
    weight(33) =   0.49411833039918178967039646116705D-01
    weight(34) =   0.49230380423747560785043116988145D-01
    weight(35) =   0.48928452820511989944709361549215D-01
    weight(36) =   0.48506789097883847864090099145802D-01
    weight(37) =   0.47966421137995131411052756195132D-01
    weight(38) =   0.47308671312268919080604988338844D-01
    weight(39) =   0.46535149245383696510395418746953D-01
    weight(40) =   0.45647747876292608685885992608542D-01
    weight(41) =   0.44648638825941395370332669516813D-01
    weight(42) =   0.43540267083027590798964315704401D-01
    weight(43) =   0.42325345020815822982505485403028D-01
    weight(44) =   0.41006845759666398635110037009072D-01
    weight(45) =   0.39587995891544093984807928149202D-01
    weight(46) =   0.38072267584349556763638324927889D-01
    weight(47) =   0.36463370085457289630452409787542D-01
    weight(48) =   0.34765240645355877697180504642788D-01
    weight(49) =   0.32982034883779341765683179672459D-01
    weight(50) =   0.31118116622219817508215988557189D-01
    weight(51) =   0.29178047208280526945551502154029D-01
    weight(52) =   0.27166574359097933225189839439413D-01
    weight(53) =   0.25088620553344986618630138068443D-01
    weight(54) =   0.22949271004889933148942319561770D-01
    weight(55) =   0.20753761258039090775341953421471D-01
    weight(56) =   0.18507464460161270409260545805144D-01
    weight(57) =   0.16215878410338338882283672974995D-01
    weight(58) =   0.13884612616115610824866086365937D-01
    weight(59) =   0.11519376076880041750750606118707D-01
    weight(60) =   0.91259686763266563540586445877022D-02
    weight(61) =   0.67102917659601362519069109850892D-02
    weight(62) =   0.42785083468637618661951422543371D-02
    weight(63) =   0.18398745955770837880499331680577D-02

  else if ( order == 127 ) then

    weight(  1) =   0.45645726109586654495731936146574D-03
    weight(  2) =   0.10622766869538486959954760554099D-02
    weight(  3) =   0.16683488125171936761028811985672D-02
    weight(  4) =   0.22734860707492547802810838362671D-02
    weight(  5) =   0.28772587656289004082883197417581D-02
    weight(  6) =   0.34792893810051465908910894094105D-02
    weight(  7) =   0.40792095178254605327114733456293D-02
    weight(  8) =   0.46766539777779034772638165662478D-02
    weight(  9) =   0.52712596565634400891303815906251D-02
    weight( 10) =   0.58626653903523901033648343751367D-02
    weight( 11) =   0.64505120486899171845442463868748D-02
    weight( 12) =   0.70344427036681608755685893032552D-02
    weight( 13) =   0.76141028256526859356393930849227D-02
    weight( 14) =   0.81891404887415730817235884718726D-02
    weight( 15) =   0.87592065795403145773316804234385D-02
    weight( 16) =   0.93239550065309714787536985834029D-02
    weight( 17) =   0.98830429087554914716648010899606D-02
    weight( 18) =   0.10436130863141005225673171997668D-01
    weight( 19) =   0.10982883090068975788799657376065D-01
    weight( 20) =   0.11522967656921087154811609734510D-01
    weight( 21) =   0.12056056679400848183529562144697D-01
    weight( 22) =   0.12581826520465013101514365424172D-01
    weight( 23) =   0.13099957986718627426172681912499D-01
    weight( 24) =   0.13610136522139249906034237533759D-01
    weight( 25) =   0.14112052399003395774044161633613D-01
    weight( 26) =   0.14605400905893418351737288078952D-01
    weight( 27) =   0.15089882532666922992635733981431D-01
    weight( 28) =   0.15565203152273955098532590262975D-01
    weight( 29) =   0.16031074199309941802254151842763D-01
    weight( 30) =   0.16487212845194879399346060358146D-01
    weight( 31) =   0.16933342169871654545878815295200D-01
    weight( 32) =   0.17369191329918731922164721250350D-01
    weight( 33) =   0.17794495722974774231027912900351D-01
    weight( 34) =   0.18208997148375106468721469154479D-01
    weight( 35) =   0.18612443963902310429440419898958D-01
    weight( 36) =   0.19004591238555646611148901044533D-01
    weight( 37) =   0.19385200901246454628112623489471D-01
    weight( 38) =   0.19754041885329183081815217323169D-01
    weight( 39) =   0.20110890268880247225644623956287D-01
    weight( 40) =   0.20455529410639508279497065713301D-01
    weight( 41) =   0.20787750081531811812652137291250D-01
    weight( 42) =   0.21107350591688713643523847921658D-01
    weight( 43) =   0.21414136912893259295449693233545D-01
    weight( 44) =   0.21707922796373466052301324695331D-01
    weight( 45) =   0.21988529885872983756478409758807D-01
    weight( 46) =   0.22255787825930280235631416460158D-01
    weight( 47) =   0.22509534365300608085694429903050D-01
    weight( 48) =   0.22749615455457959852242553240982D-01
    weight( 49) =   0.22975885344117206754377437838947D-01
    weight( 50) =   0.23188206663719640249922582981729D-01
    weight( 51) =   0.23386450514828194170722043496950D-01
    weight( 52) =   0.23570496544381716050033676844306D-01
    weight( 53) =   0.23740233018760777777714726703424D-01
    weight( 54) =   0.23895556891620665983864481754172D-01
    weight( 55) =   0.24036373866450369675132086026456D-01
    weight( 56) =   0.24162598453819584716522917710986D-01
    weight( 57) =   0.24274154023278979833195063936748D-01
    weight( 58) =   0.24370972849882214952813561907241D-01
    weight( 59) =   0.24452996155301467956140198471529D-01
    weight( 60) =   0.24520174143511508275183033290175D-01
    weight( 61) =   0.24572466031020653286354137335186D-01
    weight( 62) =   0.24609840071630254092545634003360D-01
    weight( 63) =   0.24632273575707679066033370218017D-01
    weight( 64) =   0.24639752923961094419579417477503D-01
    weight( 65) =   0.24632273575707679066033370218017D-01
    weight( 66) =   0.24609840071630254092545634003360D-01
    weight( 67) =   0.24572466031020653286354137335186D-01
    weight( 68) =   0.24520174143511508275183033290175D-01
    weight( 69) =   0.24452996155301467956140198471529D-01
    weight( 70) =   0.24370972849882214952813561907241D-01
    weight( 71) =   0.24274154023278979833195063936748D-01
    weight( 72) =   0.24162598453819584716522917710986D-01
    weight( 73) =   0.24036373866450369675132086026456D-01
    weight( 74) =   0.23895556891620665983864481754172D-01
    weight( 75) =   0.23740233018760777777714726703424D-01
    weight( 76) =   0.23570496544381716050033676844306D-01
    weight( 77) =   0.23386450514828194170722043496950D-01
    weight( 78) =   0.23188206663719640249922582981729D-01
    weight( 79) =   0.22975885344117206754377437838947D-01
    weight( 80) =   0.22749615455457959852242553240982D-01
    weight( 81) =   0.22509534365300608085694429903050D-01
    weight( 82) =   0.22255787825930280235631416460158D-01
    weight( 83) =   0.21988529885872983756478409758807D-01
    weight( 84) =   0.21707922796373466052301324695331D-01
    weight( 85) =   0.21414136912893259295449693233545D-01
    weight( 86) =   0.21107350591688713643523847921658D-01
    weight( 87) =   0.20787750081531811812652137291250D-01
    weight( 88) =   0.20455529410639508279497065713301D-01
    weight( 89) =   0.20110890268880247225644623956287D-01
    weight( 90) =   0.19754041885329183081815217323169D-01
    weight( 91) =   0.19385200901246454628112623489471D-01
    weight( 92) =   0.19004591238555646611148901044533D-01
    weight( 93) =   0.18612443963902310429440419898958D-01
    weight( 94) =   0.18208997148375106468721469154479D-01
    weight( 95) =   0.17794495722974774231027912900351D-01
    weight( 96) =   0.17369191329918731922164721250350D-01
    weight( 97) =   0.16933342169871654545878815295200D-01
    weight( 98) =   0.16487212845194879399346060358146D-01
    weight( 99) =   0.16031074199309941802254151842763D-01
    weight(100) =   0.15565203152273955098532590262975D-01
    weight(101) =   0.15089882532666922992635733981431D-01
    weight(102) =   0.14605400905893418351737288078952D-01
    weight(103) =   0.14112052399003395774044161633613D-01
    weight(104) =   0.13610136522139249906034237533759D-01
    weight(105) =   0.13099957986718627426172681912499D-01
    weight(106) =   0.12581826520465013101514365424172D-01
    weight(107) =   0.12056056679400848183529562144697D-01
    weight(108) =   0.11522967656921087154811609734510D-01
    weight(109) =   0.10982883090068975788799657376065D-01
    weight(110) =   0.10436130863141005225673171997668D-01
    weight(111) =   0.98830429087554914716648010899606D-02
    weight(112) =   0.93239550065309714787536985834029D-02
    weight(113) =   0.87592065795403145773316804234385D-02
    weight(114) =   0.81891404887415730817235884718726D-02
    weight(115) =   0.76141028256526859356393930849227D-02
    weight(116) =   0.70344427036681608755685893032552D-02
    weight(117) =   0.64505120486899171845442463868748D-02
    weight(118) =   0.58626653903523901033648343751367D-02
    weight(119) =   0.52712596565634400891303815906251D-02
    weight(120) =   0.46766539777779034772638165662478D-02
    weight(121) =   0.40792095178254605327114733456293D-02
    weight(122) =   0.34792893810051465908910894094105D-02
    weight(123) =   0.28772587656289004082883197417581D-02
    weight(124) =   0.22734860707492547802810838362671D-02
    weight(125) =   0.16683488125171936761028811985672D-02
    weight(126) =   0.10622766869538486959954760554099D-02
    weight(127) =   0.45645726109586654495731936146574D-03

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GL_WEIGHTS - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of ORDER = ', order
    write ( *, '(a)' ) '  Legal values are 1, 3, 7, 15, 31, 63 and 127.'
    stop

  end if

  return
end
function gp_abscissa ( order, i )

!*****************************************************************************80
!
!! GP_ABSCISSA returns the I-th abscissa for a Gauss Patterson rule.
!
!  Discussion:
!
!    The rule is specified by its order.
!
!    The number of points in the rule, known as the order, is
!    related to the level by the formula:
!
!      ORDER = 2^(LEVEL+1)-1.
!
!    Only rules of order 1, 3, 7, 15, 31, 63, 127 and 255 are allowed.
!
!    Since the IBM XLF FORTRAN compiler enforces the unreasonable
!    but legal limit on the number of continuation lines, I have
!    had to modify the declaration of the array holding the abscissas.
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
!    Prem Kythe, Michael Schaeferkotter,
!    Handbook of Computational Methods for Integration,
!    Chapman and Hall, 2004,
!    ISBN: 1-58488-428-2,
!    LC: QA299.3.K98.
!
!    Thomas Patterson,
!    The Optimal Addition of Points to Quadrature Formulae,
!    Mathematics of Computation,
!    Volume 22, Number 104, October 1968, pages 847-856.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ORDER, the order.
!    ORDER must be 1, 3, 7, 15, 31, 63, 127 or 255.
!
!    Input, integer ( kind = 4 ) I, the index of the point in the rule.
!
!    Output, real ( kind = 8 ) GP_ABSCISSA, the value of the I-th
!    abscissa in the rule of order ORDER.
!
  implicit none

  real ( kind = 8 ) gp_abscissa
  integer ( kind = 4 ) i
  integer ( kind = 4 ) order
  real ( kind = 8 ) value
  real ( kind = 8 ), save, dimension ( 1 ) :: x_001 = (/ &
     0.0D+00 /)
  real ( kind = 8 ), save, dimension ( 3 ) :: x_003 = (/ &
    -0.77459666924148337704D+00, &
     0.0D+00, &
     0.77459666924148337704D+00 /)
  real ( kind = 8 ), save, dimension ( 7 ) :: x_007 = (/ &
    -0.96049126870802028342D+00, &
    -0.77459666924148337704D+00, &
    -0.43424374934680255800D+00, &
     0.0D+00, &
     0.43424374934680255800D+00, &
     0.77459666924148337704D+00, &
     0.96049126870802028342D+00 /)
  real ( kind = 8 ), save, dimension ( 15 ) :: x_015 = (/ &
    -0.99383196321275502221D+00, &
    -0.96049126870802028342D+00, &
    -0.88845923287225699889D+00, &
    -0.77459666924148337704D+00, &
    -0.62110294673722640294D+00, &
    -0.43424374934680255800D+00, &
    -0.22338668642896688163D+00, &
     0.0D+00, &
     0.22338668642896688163D+00, &
     0.43424374934680255800D+00, &
     0.62110294673722640294D+00, &
     0.77459666924148337704D+00, &
     0.88845923287225699889D+00, &
     0.96049126870802028342D+00, &
     0.99383196321275502221D+00 /)
  real ( kind = 8 ), save, dimension ( 31 ) :: x_031 = (/ &
    -0.99909812496766759766D+00, &
    -0.99383196321275502221D+00, &
    -0.98153114955374010687D+00, &
    -0.96049126870802028342D+00, &
    -0.92965485742974005667D+00, &
    -0.88845923287225699889D+00, &
    -0.83672593816886873550D+00, &
    -0.77459666924148337704D+00, &
    -0.70249620649152707861D+00, &
    -0.62110294673722640294D+00, &
    -0.53131974364437562397D+00, &
    -0.43424374934680255800D+00, &
    -0.33113539325797683309D+00, &
    -0.22338668642896688163D+00, &
    -0.11248894313318662575D+00, &
     0.0D+00, &
     0.11248894313318662575D+00, &
     0.22338668642896688163D+00, &
     0.33113539325797683309D+00, &
     0.43424374934680255800D+00, &
     0.53131974364437562397D+00, &
     0.62110294673722640294D+00, &
     0.70249620649152707861D+00, &
     0.77459666924148337704D+00, &
     0.83672593816886873550D+00, &
     0.88845923287225699889D+00, &
     0.92965485742974005667D+00, &
     0.96049126870802028342D+00, &
     0.98153114955374010687D+00, &
     0.99383196321275502221D+00, &
     0.99909812496766759766D+00 /)
  real ( kind = 8 ), save, dimension ( 63 ) :: x_063 = (/ &
    -0.99987288812035761194D+00, &
    -0.99909812496766759766D+00, &
    -0.99720625937222195908D+00, &
    -0.99383196321275502221D+00, &
    -0.98868475754742947994D+00, &
    -0.98153114955374010687D+00, &
    -0.97218287474858179658D+00, &
    -0.96049126870802028342D+00, &
    -0.94634285837340290515D+00, &
    -0.92965485742974005667D+00, &
    -0.91037115695700429250D+00, &
    -0.88845923287225699889D+00, &
    -0.86390793819369047715D+00, &
    -0.83672593816886873550D+00, &
    -0.80694053195021761186D+00, &
    -0.77459666924148337704D+00, &
    -0.73975604435269475868D+00, &
    -0.70249620649152707861D+00, &
    -0.66290966002478059546D+00, &
    -0.62110294673722640294D+00, &
    -0.57719571005204581484D+00, &
    -0.53131974364437562397D+00, &
    -0.48361802694584102756D+00, &
    -0.43424374934680255800D+00, &
    -0.38335932419873034692D+00, &
    -0.33113539325797683309D+00, &
    -0.27774982202182431507D+00, &
    -0.22338668642896688163D+00, &
    -0.16823525155220746498D+00, &
    -0.11248894313318662575D+00, &
    -0.056344313046592789972D+00, &
     0.0D+00, &
     0.056344313046592789972D+00, &
     0.11248894313318662575D+00, &
     0.16823525155220746498D+00, &
     0.22338668642896688163D+00, &
     0.27774982202182431507D+00, &
     0.33113539325797683309D+00, &
     0.38335932419873034692D+00, &
     0.43424374934680255800D+00, &
     0.48361802694584102756D+00, &
     0.53131974364437562397D+00, &
     0.57719571005204581484D+00, &
     0.62110294673722640294D+00, &
     0.66290966002478059546D+00, &
     0.70249620649152707861D+00, &
     0.73975604435269475868D+00, &
     0.77459666924148337704D+00, &
     0.80694053195021761186D+00, &
     0.83672593816886873550D+00, &
     0.86390793819369047715D+00, &
     0.88845923287225699889D+00, &
     0.91037115695700429250D+00, &
     0.92965485742974005667D+00, &
     0.94634285837340290515D+00, &
     0.96049126870802028342D+00, &
     0.97218287474858179658D+00, &
     0.98153114955374010687D+00, &
     0.98868475754742947994D+00, &
     0.99383196321275502221D+00, &
     0.99720625937222195908D+00, &
     0.99909812496766759766D+00, &
     0.99987288812035761194D+00 /)
  real ( kind = 8 ), save, dimension ( 127 ) :: x_127 = (/ &
    -0.99998243035489159858D+00, &
    -0.99987288812035761194D+00, &
    -0.99959879967191068325D+00, &
    -0.99909812496766759766D+00, &
    -0.99831663531840739253D+00, &
    -0.99720625937222195908D+00, &
    -0.99572410469840718851D+00, &
    -0.99383196321275502221D+00, &
    -0.99149572117810613240D+00, &
    -0.98868475754742947994D+00, &
    -0.98537149959852037111D+00, &
    -0.98153114955374010687D+00, &
    -0.97714151463970571416D+00, &
    -0.97218287474858179658D+00, &
    -0.96663785155841656709D+00, &
    -0.96049126870802028342D+00, &
    -0.95373000642576113641D+00, &
    -0.94634285837340290515D+00, &
    -0.93832039777959288365D+00, &
    -0.92965485742974005667D+00, &
    -0.92034002547001242073D+00, &
    -0.91037115695700429250D+00, &
    -0.89974489977694003664D+00, &
    -0.88845923287225699889D+00, &
    -0.87651341448470526974D+00, &
    -0.86390793819369047715D+00, &
    -0.85064449476835027976D+00, &
    -0.83672593816886873550D+00, &
    -0.82215625436498040737D+00, &
    -0.80694053195021761186D+00, &
    -0.79108493379984836143D+00, &
    -0.77459666924148337704D+00, &
    -0.75748396638051363793D+00, &
    -0.73975604435269475868D+00, &
    -0.72142308537009891548D+00, &
    -0.70249620649152707861D+00, &
    -0.68298743109107922809D+00, &
    -0.66290966002478059546D+00, &
    -0.64227664250975951377D+00, &
    -0.62110294673722640294D+00, &
    -0.59940393024224289297D+00, &
    -0.57719571005204581484D+00, &
    -0.55449513263193254887D+00, &
    -0.53131974364437562397D+00, &
    -0.50768775753371660215D+00, &
    -0.48361802694584102756D+00, &
    -0.45913001198983233287D+00, &
    -0.43424374934680255800D+00, &
    -0.40897982122988867241D+00, &
    -0.38335932419873034692D+00, &
    -0.35740383783153215238D+00, &
    -0.33113539325797683309D+00, &
    -0.30457644155671404334D+00, &
    -0.27774982202182431507D+00, &
    -0.25067873030348317661D+00, &
    -0.22338668642896688163D+00, &
    -0.19589750271110015392D+00, &
    -0.16823525155220746498D+00, &
    -0.14042423315256017459D+00, &
    -0.11248894313318662575D+00, &
    -0.084454040083710883710D+00, &
    -0.056344313046592789972D+00, &
    -0.028184648949745694339D+00, &
     0.0D+00, &
     0.028184648949745694339D+00, &
     0.056344313046592789972D+00, &
     0.084454040083710883710D+00, &
     0.11248894313318662575D+00, &
     0.14042423315256017459D+00, &
     0.16823525155220746498D+00, &
     0.19589750271110015392D+00, &
     0.22338668642896688163D+00, &
     0.25067873030348317661D+00, &
     0.27774982202182431507D+00, &
     0.30457644155671404334D+00, &
     0.33113539325797683309D+00, &
     0.35740383783153215238D+00, &
     0.38335932419873034692D+00, &
     0.40897982122988867241D+00, &
     0.43424374934680255800D+00, &
     0.45913001198983233287D+00, &
     0.48361802694584102756D+00, &
     0.50768775753371660215D+00, &
     0.53131974364437562397D+00, &
     0.55449513263193254887D+00, &
     0.57719571005204581484D+00, &
     0.59940393024224289297D+00, &
     0.62110294673722640294D+00, &
     0.64227664250975951377D+00, &
     0.66290966002478059546D+00, &
     0.68298743109107922809D+00, &
     0.70249620649152707861D+00, &
     0.72142308537009891548D+00, &
     0.73975604435269475868D+00, &
     0.75748396638051363793D+00, &
     0.77459666924148337704D+00, &
     0.79108493379984836143D+00, &
     0.80694053195021761186D+00, &
     0.82215625436498040737D+00, &
     0.83672593816886873550D+00, &
     0.85064449476835027976D+00, &
     0.86390793819369047715D+00, &
     0.87651341448470526974D+00, &
     0.88845923287225699889D+00, &
     0.89974489977694003664D+00, &
     0.91037115695700429250D+00, &
     0.92034002547001242073D+00, &
     0.92965485742974005667D+00, &
     0.93832039777959288365D+00, &
     0.94634285837340290515D+00, &
     0.95373000642576113641D+00, &
     0.96049126870802028342D+00, &
     0.96663785155841656709D+00, &
     0.97218287474858179658D+00, &
     0.97714151463970571416D+00, &
     0.98153114955374010687D+00, &
     0.98537149959852037111D+00, &
     0.98868475754742947994D+00, &
     0.99149572117810613240D+00, &
     0.99383196321275502221D+00, &
     0.99572410469840718851D+00, &
     0.99720625937222195908D+00, &
     0.99831663531840739253D+00, &
     0.99909812496766759766D+00, &
     0.99959879967191068325D+00, &
     0.99987288812035761194D+00, &
     0.99998243035489159858D+00 /)

  real ( kind = 8 ), save, dimension ( 255 ) :: x_255 = (/ &
    -0.99999759637974846462D+00, &
    -0.99998243035489159858D+00, &
    -0.99994399620705437576D+00, &
    -0.99987288812035761194D+00, &
    -0.99976049092443204733D+00, &
    -0.99959879967191068325D+00, &
    -0.99938033802502358193D+00, &
    -0.99909812496766759766D+00, &
    -0.99874561446809511470D+00, &
    -0.99831663531840739253D+00, &
    -0.99780535449595727456D+00, &
    -0.99720625937222195908D+00, &
    -0.99651414591489027385D+00, &
    -0.99572410469840718851D+00, &
    -0.99483150280062100052D+00, &
    -0.99383196321275502221D+00, &
    -0.99272134428278861533D+00, &
    -0.99149572117810613240D+00, &
    -0.99015137040077015918D+00, &
    -0.98868475754742947994D+00, &
    -0.98709252795403406719D+00, &
    -0.98537149959852037111D+00, &
    -0.98351865757863272876D+00, &
    -0.98153114955374010687D+00, &
    -0.97940628167086268381D+00, &
    -0.97714151463970571416D+00, &
    -0.97473445975240266776D+00, &
    -0.97218287474858179658D+00, &
    -0.96948465950245923177D+00, &
    -0.96663785155841656709D+00, &
    -0.96364062156981213252D+00, &
    -0.96049126870802028342D+00, &
    -0.95718821610986096274D+00, &
    -0.95373000642576113641D+00, &
    -0.95011529752129487656D+00, &
    -0.94634285837340290515D+00, &
    -0.94241156519108305981D+00, &
    -0.93832039777959288365D+00, &
    -0.93406843615772578800D+00, &
    -0.92965485742974005667D+00, &
    -0.92507893290707565236D+00, &
    -0.92034002547001242073D+00, &
    -0.91543758715576504064D+00, &
    -0.91037115695700429250D+00, &
    -0.90514035881326159519D+00, &
    -0.89974489977694003664D+00, &
    -0.89418456833555902286D+00, &
    -0.88845923287225699889D+00, &
    -0.88256884024734190684D+00, &
    -0.87651341448470526974D+00, &
    -0.87029305554811390585D+00, &
    -0.86390793819369047715D+00, &
    -0.85735831088623215653D+00, &
    -0.85064449476835027976D+00, &
    -0.84376688267270860104D+00, &
    -0.83672593816886873550D+00, &
    -0.82952219463740140018D+00, &
    -0.82215625436498040737D+00, &
    -0.81462878765513741344D+00, &
    -0.80694053195021761186D+00, &
    -0.79909229096084140180D+00, &
    -0.79108493379984836143D+00, &
    -0.78291939411828301639D+00, &
    -0.77459666924148337704D+00, &
    -0.76611781930376009072D+00, &
    -0.75748396638051363793D+00, &
    -0.74869629361693660282D+00, &
    -0.73975604435269475868D+00, &
    -0.73066452124218126133D+00, &
    -0.72142308537009891548D+00, &
    -0.71203315536225203459D+00, &
    -0.70249620649152707861D+00, &
    -0.69281376977911470289D+00, &
    -0.68298743109107922809D+00, &
    -0.67301883023041847920D+00, &
    -0.66290966002478059546D+00, &
    -0.65266166541001749610D+00, &
    -0.64227664250975951377D+00, &
    -0.63175643771119423041D+00, &
    -0.62110294673722640294D+00, &
    -0.61031811371518640016D+00, &
    -0.59940393024224289297D+00, &
    -0.58836243444766254143D+00, &
    -0.57719571005204581484D+00, &
    -0.56590588542365442262D+00, &
    -0.55449513263193254887D+00, &
    -0.54296566649831149049D+00, &
    -0.53131974364437562397D+00, &
    -0.51955966153745702199D+00, &
    -0.50768775753371660215D+00, &
    -0.49570640791876146017D+00, &
    -0.48361802694584102756D+00, &
    -0.47142506587165887693D+00, &
    -0.45913001198983233287D+00, &
    -0.44673538766202847374D+00, &
    -0.43424374934680255800D+00, &
    -0.42165768662616330006D+00, &
    -0.40897982122988867241D+00, &
    -0.39621280605761593918D+00, &
    -0.38335932419873034692D+00, &
    -0.37042208795007823014D+00, &
    -0.35740383783153215238D+00, &
    -0.34430734159943802278D+00, &
    -0.33113539325797683309D+00, &
    -0.31789081206847668318D+00, &
    -0.30457644155671404334D+00, &
    -0.29119514851824668196D+00, &
    -0.27774982202182431507D+00, &
    -0.26424337241092676194D+00, &
    -0.25067873030348317661D+00, &
    -0.23705884558982972721D+00, &
    -0.22338668642896688163D+00, &
    -0.20966523824318119477D+00, &
    -0.19589750271110015392D+00, &
    -0.18208649675925219825D+00, &
    -0.16823525155220746498D+00, &
    -0.15434681148137810869D+00, &
    -0.14042423315256017459D+00, &
    -0.12647058437230196685D+00, &
    -0.11248894313318662575D+00, &
    -0.098482396598119202090D+00, &
    -0.084454040083710883710D+00, &
    -0.070406976042855179063D+00, &
    -0.056344313046592789972D+00, &
    -0.042269164765363603212D+00, &
    -0.028184648949745694339D+00, &
    -0.014093886410782462614D+00, &
    0.0D+00, &
    0.014093886410782462614D+00, &
    0.028184648949745694339D+00, &
    0.042269164765363603212D+00, &
    0.056344313046592789972D+00, &
    0.070406976042855179063D+00, &
    0.084454040083710883710D+00, &
    0.098482396598119202090D+00, &
    0.11248894313318662575D+00, &
    0.12647058437230196685D+00, &
    0.14042423315256017459D+00, &
    0.15434681148137810869D+00, &
    0.16823525155220746498D+00, &
    0.18208649675925219825D+00, &
    0.19589750271110015392D+00, &
    0.20966523824318119477D+00, &
    0.22338668642896688163D+00, &
    0.23705884558982972721D+00, &
    0.25067873030348317661D+00, &
    0.26424337241092676194D+00, &
    0.27774982202182431507D+00, &
    0.29119514851824668196D+00, &
    0.30457644155671404334D+00, &
    0.31789081206847668318D+00, &
    0.33113539325797683309D+00, &
    0.34430734159943802278D+00, &
    0.35740383783153215238D+00, &
    0.37042208795007823014D+00, &
    0.38335932419873034692D+00, &
    0.39621280605761593918D+00, &
    0.40897982122988867241D+00, &
    0.42165768662616330006D+00, &
    0.43424374934680255800D+00, &
    0.44673538766202847374D+00, &
    0.45913001198983233287D+00, &
    0.47142506587165887693D+00, &
    0.48361802694584102756D+00, &
    0.49570640791876146017D+00, &
    0.50768775753371660215D+00, &
    0.51955966153745702199D+00, &
    0.53131974364437562397D+00, &
    0.54296566649831149049D+00, &
    0.55449513263193254887D+00, &
    0.56590588542365442262D+00, &
    0.57719571005204581484D+00, &
    0.58836243444766254143D+00, &
    0.59940393024224289297D+00, &
    0.61031811371518640016D+00, &
    0.62110294673722640294D+00, &
    0.63175643771119423041D+00, &
    0.64227664250975951377D+00, &
    0.65266166541001749610D+00, &
    0.66290966002478059546D+00, &
    0.67301883023041847920D+00, &
    0.68298743109107922809D+00, &
    0.69281376977911470289D+00, &
    0.70249620649152707861D+00, &
    0.71203315536225203459D+00, &
    0.72142308537009891548D+00, &
    0.73066452124218126133D+00, &
    0.73975604435269475868D+00, &
    0.74869629361693660282D+00, &
    0.75748396638051363793D+00, &
    0.76611781930376009072D+00, &
    0.77459666924148337704D+00, &
    0.78291939411828301639D+00, &
    0.79108493379984836143D+00, &
    0.79909229096084140180D+00, &
    0.80694053195021761186D+00, &
    0.81462878765513741344D+00, &
    0.82215625436498040737D+00, &
    0.82952219463740140018D+00, &
    0.83672593816886873550D+00, &
    0.84376688267270860104D+00, &
    0.85064449476835027976D+00, &
    0.85735831088623215653D+00, &
    0.86390793819369047715D+00, &
    0.87029305554811390585D+00, &
    0.87651341448470526974D+00, &
    0.88256884024734190684D+00, &
    0.88845923287225699889D+00, &
    0.89418456833555902286D+00, &
    0.89974489977694003664D+00, &
    0.90514035881326159519D+00, &
    0.91037115695700429250D+00, &
    0.91543758715576504064D+00, &
    0.92034002547001242073D+00, &
    0.92507893290707565236D+00, &
    0.92965485742974005667D+00, &
    0.93406843615772578800D+00, &
    0.93832039777959288365D+00, &
    0.94241156519108305981D+00, &
    0.94634285837340290515D+00, &
    0.95011529752129487656D+00, &
    0.95373000642576113641D+00, &
    0.95718821610986096274D+00, &
    0.96049126870802028342D+00, &
    0.96364062156981213252D+00, &
    0.96663785155841656709D+00, &
    0.96948465950245923177D+00, &
    0.97218287474858179658D+00, &
    0.97473445975240266776D+00, &
    0.97714151463970571416D+00, &
    0.97940628167086268381D+00, &
    0.98153114955374010687D+00, &
    0.98351865757863272876D+00, &
    0.98537149959852037111D+00, &
    0.98709252795403406719D+00, &
    0.98868475754742947994D+00, &
    0.99015137040077015918D+00, &
    0.99149572117810613240D+00, &
    0.99272134428278861533D+00, &
    0.99383196321275502221D+00, &
    0.99483150280062100052D+00, &
    0.99572410469840718851D+00, &
    0.99651414591489027385D+00, &
    0.99720625937222195908D+00, &
    0.99780535449595727456D+00, &
    0.99831663531840739253D+00, &
    0.99874561446809511470D+00, &
    0.99909812496766759766D+00, &
    0.99938033802502358193D+00, &
    0.99959879967191068325D+00, &
    0.99976049092443204733D+00, &
    0.99987288812035761194D+00, &
    0.99994399620705437576D+00, &
    0.99998243035489159858D+00, &
    0.99999759637974846462D+00 /)

  if ( i < 1 .or. order < i ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GP_ABSCISSA - Fatal error!'
    write ( *, '(a)' ) '  I < 1 or ORDER < I.'
    write ( *, '(a,i12)' ) '  I = ', i
    write ( *, '(a,i12)' ) '  ORDER = ', order
    stop
  end if

  if ( order == 1 ) then
    value = x_001(i)
  else if ( order == 3 ) then
    value = x_003(i)
  else if ( order == 7 ) then
    value = x_007(i)
  else if ( order == 15 ) then
    value = x_015(i)
  else if ( order == 31 ) then
    value = x_031(i)
  else if ( order == 63 ) then
    value = x_063(i)
  else if ( order == 127 ) then
    value = x_127(i)
  else if ( order == 255 ) then
    value = x_255(i)
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GP_ABSCISSA - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal input value of ORDER = ', order
    stop
  end if

  gp_abscissa = value

  return
end
subroutine gp_weights ( order, w )

!*****************************************************************************80
!
!! GP_WEIGHTS sets weights for a Gauss Patterson rule.
!
!  Discussion:
!
!    The zeroth rule, of order 1, is the standard Gauss Legendre rule.
!
!    The first rule, of order 3, is the standard Gauss Legendre rule.
!
!    The second rule, of order 7, includes the abscissas of the previous
!    rule.
!
!    Each subsequent rule is nested in a similar way.  Rules are available
!    of orders 1, 3, 7, 15, 31, 63. 127 or 255.
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
!    Prem Kythe, Michael Schaeferkotter,
!    Handbook of Computational Methods for Integration,
!    Chapman and Hall, 2004,
!    ISBN: 1-58488-428-2,
!    LC: QA299.3.K98.
!
!    Thomas Patterson,
!    The Optimal Addition of Points to Quadrature Formulae,
!    Mathematics of Computation,
!    Volume 22, Number 104, October 1968, pages 847-856.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ORDER, the order.
!    ORDER must be 1, 3, 7, 15, 31, 63, 127 or 255.
!
!    Output, real ( kind = 8 ) W(ORDER), the weights.
!    The weights are positive, symmetric and should sum to 2.
!
  implicit none

  integer ( kind = 4 ) order

  real ( kind = 8 ) w(order)

  if ( order == 1 ) then

    w(1) = 2.0D+00

  else if ( order == 3 ) then

    w(1) = 0.555555555555555555556D+00
    w(2) = 0.888888888888888888889D+00
    w(3) = 0.555555555555555555556D+00

  else if ( order == 7 ) then

    w(1) = 0.104656226026467265194D+00
    w(2) = 0.268488089868333440729D+00
    w(3) = 0.401397414775962222905D+00
    w(4) = 0.450916538658474142345D+00
    w(5) = 0.401397414775962222905D+00
    w(6) = 0.268488089868333440729D+00
    w(7) = 0.104656226026467265194D+00

  else if ( order == 15 ) then

    w( 1) = 0.0170017196299402603390D+00
    w( 2) = 0.0516032829970797396969D+00
    w( 3) = 0.0929271953151245376859D+00
    w( 4) = 0.134415255243784220360D+00
    w( 5) = 0.171511909136391380787D+00
    w( 6) = 0.200628529376989021034D+00
    w( 7) = 0.219156858401587496404D+00
    w( 8) = 0.225510499798206687386D+00
    w( 9) = 0.219156858401587496404D+00
    w(10) = 0.200628529376989021034D+00
    w(11) = 0.171511909136391380787D+00
    w(12) = 0.134415255243784220360D+00
    w(13) = 0.0929271953151245376859D+00
    w(14) = 0.0516032829970797396969D+00
    w(15) = 0.0170017196299402603390D+00

  else if ( order == 31 ) then

    w( 1) = 0.00254478079156187441540D+00
    w( 2) = 0.00843456573932110624631D+00
    w( 3) = 0.0164460498543878109338D+00
    w( 4) = 0.0258075980961766535646D+00
    w( 5) = 0.0359571033071293220968D+00
    w( 6) = 0.0464628932617579865414D+00
    w( 7) = 0.0569795094941233574122D+00
    w( 8) = 0.0672077542959907035404D+00
    w( 9) = 0.0768796204990035310427D+00
    w(10) = 0.0857559200499903511542D+00
    w(11) = 0.0936271099812644736167D+00
    w(12) = 0.100314278611795578771D+00
    w(13) = 0.105669893580234809744D+00
    w(14) = 0.109578421055924638237D+00
    w(15) = 0.111956873020953456880D+00
    w(16) = 0.112755256720768691607D+00
    w(17) = 0.111956873020953456880D+00
    w(18) = 0.109578421055924638237D+00
    w(19) = 0.105669893580234809744D+00
    w(20) = 0.100314278611795578771D+00
    w(21) = 0.0936271099812644736167D+00
    w(22) = 0.0857559200499903511542D+00
    w(23) = 0.0768796204990035310427D+00
    w(24) = 0.0672077542959907035404D+00
    w(25) = 0.0569795094941233574122D+00
    w(26) = 0.0464628932617579865414D+00
    w(27) = 0.0359571033071293220968D+00
    w(28) = 0.0258075980961766535646D+00
    w(29) = 0.0164460498543878109338D+00
    w(30) = 0.00843456573932110624631D+00
    w(31) = 0.00254478079156187441540D+00

  else if ( order == 63 ) then

    w( 1) = 0.000363221481845530659694D+00
    w( 2) = 0.00126515655623006801137D+00
    w( 3) = 0.00257904979468568827243D+00
    w( 4) = 0.00421763044155885483908D+00
    w( 5) = 0.00611550682211724633968D+00
    w( 6) = 0.00822300795723592966926D+00
    w( 7) = 0.0104982469096213218983D+00
    w( 8) = 0.0129038001003512656260D+00
    w( 9) = 0.0154067504665594978021D+00
    w(10) = 0.0179785515681282703329D+00
    w(11) = 0.0205942339159127111492D+00
    w(12) = 0.0232314466399102694433D+00
    w(13) = 0.0258696793272147469108D+00
    w(14) = 0.0284897547458335486125D+00
    w(15) = 0.0310735511116879648799D+00
    w(16) = 0.0336038771482077305417D+00
    w(17) = 0.0360644327807825726401D+00
    w(18) = 0.0384398102494555320386D+00
    w(19) = 0.0407155101169443189339D+00
    w(20) = 0.0428779600250077344929D+00
    w(21) = 0.0449145316536321974143D+00
    w(22) = 0.0468135549906280124026D+00
    w(23) = 0.0485643304066731987159D+00
    w(24) = 0.0501571393058995374137D+00
    w(25) = 0.0515832539520484587768D+00
    w(26) = 0.0528349467901165198621D+00
    w(27) = 0.0539054993352660639269D+00
    w(28) = 0.0547892105279628650322D+00
    w(29) = 0.0554814043565593639878D+00
    w(30) = 0.0559784365104763194076D+00
    w(31) = 0.0562776998312543012726D+00
    w(32) = 0.0563776283603847173877D+00
    w(33) = 0.0562776998312543012726D+00
    w(34) = 0.0559784365104763194076D+00
    w(35) = 0.0554814043565593639878D+00
    w(36) = 0.0547892105279628650322D+00
    w(37) = 0.0539054993352660639269D+00
    w(38) = 0.0528349467901165198621D+00
    w(39) = 0.0515832539520484587768D+00
    w(40) = 0.0501571393058995374137D+00
    w(41) = 0.0485643304066731987159D+00
    w(42) = 0.0468135549906280124026D+00
    w(43) = 0.0449145316536321974143D+00
    w(44) = 0.0428779600250077344929D+00
    w(45) = 0.0407155101169443189339D+00
    w(46) = 0.0384398102494555320386D+00
    w(47) = 0.0360644327807825726401D+00
    w(48) = 0.0336038771482077305417D+00
    w(49) = 0.0310735511116879648799D+00
    w(50) = 0.0284897547458335486125D+00
    w(51) = 0.0258696793272147469108D+00
    w(52) = 0.0232314466399102694433D+00
    w(53) = 0.0205942339159127111492D+00
    w(54) = 0.0179785515681282703329D+00
    w(55) = 0.0154067504665594978021D+00
    w(56) = 0.0129038001003512656260D+00
    w(57) = 0.0104982469096213218983D+00
    w(58) = 0.00822300795723592966926D+00
    w(59) = 0.00611550682211724633968D+00
    w(60) = 0.00421763044155885483908D+00
    w(61) = 0.00257904979468568827243D+00
    w(62) = 0.00126515655623006801137D+00
    w(63) = 0.000363221481845530659694D+00

  else if ( order == 127 ) then

    w(  1) = 0.0000505360952078625176247D+00
    w(  2) = 0.000180739564445388357820D+00
    w(  3) = 0.000377746646326984660274D+00
    w(  4) = 0.000632607319362633544219D+00
    w(  5) = 0.000938369848542381500794D+00
    w(  6) = 0.00128952408261041739210D+00
    w(  7) = 0.00168114286542146990631D+00
    w(  8) = 0.00210881524572663287933D+00
    w(  9) = 0.00256876494379402037313D+00
    w( 10) = 0.00305775341017553113613D+00
    w( 11) = 0.00357289278351729964938D+00
    w( 12) = 0.00411150397865469304717D+00
    w( 13) = 0.00467105037211432174741D+00
    w( 14) = 0.00524912345480885912513D+00
    w( 15) = 0.00584344987583563950756D+00
    w( 16) = 0.00645190005017573692280D+00
    w( 17) = 0.00707248999543355546805D+00
    w( 18) = 0.00770337523327974184817D+00
    w( 19) = 0.00834283875396815770558D+00
    w( 20) = 0.00898927578406413572328D+00
    w( 21) = 0.00964117772970253669530D+00
    w( 22) = 0.0102971169579563555237D+00
    w( 23) = 0.0109557333878379016480D+00
    w( 24) = 0.0116157233199551347270D+00
    w( 25) = 0.0122758305600827700870D+00
    w( 26) = 0.0129348396636073734547D+00
    w( 27) = 0.0135915710097655467896D+00
    w( 28) = 0.0142448773729167743063D+00
    w( 29) = 0.0148936416648151820348D+00
    w( 30) = 0.0155367755558439824399D+00
    w( 31) = 0.0161732187295777199419D+00
    w( 32) = 0.0168019385741038652709D+00
    w( 33) = 0.0174219301594641737472D+00
    w( 34) = 0.0180322163903912863201D+00
    w( 35) = 0.0186318482561387901863D+00
    w( 36) = 0.0192199051247277660193D+00
    w( 37) = 0.0197954950480974994880D+00
    w( 38) = 0.0203577550584721594669D+00
    w( 39) = 0.0209058514458120238522D+00
    w( 40) = 0.0214389800125038672465D+00
    w( 41) = 0.0219563663053178249393D+00
    w( 42) = 0.0224572658268160987071D+00
    w( 43) = 0.0229409642293877487608D+00
    w( 44) = 0.0234067774953140062013D+00
    w( 45) = 0.0238540521060385400804D+00
    w( 46) = 0.0242821652033365993580D+00
    w( 47) = 0.0246905247444876769091D+00
    w( 48) = 0.0250785696529497687068D+00
    w( 49) = 0.0254457699654647658126D+00
    w( 50) = 0.0257916269760242293884D+00
    w( 51) = 0.0261156733767060976805D+00
    w( 52) = 0.0264174733950582599310D+00
    w( 53) = 0.0266966229274503599062D+00
    w( 54) = 0.0269527496676330319634D+00
    w( 55) = 0.0271855132296247918192D+00
    w( 56) = 0.0273946052639814325161D+00
    w( 57) = 0.0275797495664818730349D+00
    w( 58) = 0.0277407021782796819939D+00
    w( 59) = 0.0278772514766137016085D+00
    w( 60) = 0.0279892182552381597038D+00
    w( 61) = 0.0280764557938172466068D+00
    w( 62) = 0.0281388499156271506363D+00
    w( 63) = 0.0281763190330166021307D+00
    w( 64) = 0.0281888141801923586938D+00
    w( 65) = 0.0281763190330166021307D+00
    w( 66) = 0.0281388499156271506363D+00
    w( 67) = 0.0280764557938172466068D+00
    w( 68) = 0.0279892182552381597038D+00
    w( 69) = 0.0278772514766137016085D+00
    w( 70) = 0.0277407021782796819939D+00
    w( 71) = 0.0275797495664818730349D+00
    w( 72) = 0.0273946052639814325161D+00
    w( 73) = 0.0271855132296247918192D+00
    w( 74) = 0.0269527496676330319634D+00
    w( 75) = 0.0266966229274503599062D+00
    w( 76) = 0.0264174733950582599310D+00
    w( 77) = 0.0261156733767060976805D+00
    w( 78) = 0.0257916269760242293884D+00
    w( 79) = 0.0254457699654647658126D+00
    w( 80) = 0.0250785696529497687068D+00
    w( 81) = 0.0246905247444876769091D+00
    w( 82) = 0.0242821652033365993580D+00
    w( 83) = 0.0238540521060385400804D+00
    w( 84) = 0.0234067774953140062013D+00
    w( 85) = 0.0229409642293877487608D+00
    w( 86) = 0.0224572658268160987071D+00
    w( 87) = 0.0219563663053178249393D+00
    w( 88) = 0.0214389800125038672465D+00
    w( 89) = 0.0209058514458120238522D+00
    w( 90) = 0.0203577550584721594669D+00
    w( 91) = 0.0197954950480974994880D+00
    w( 92) = 0.0192199051247277660193D+00
    w( 93) = 0.0186318482561387901863D+00
    w( 94) = 0.0180322163903912863201D+00
    w( 95) = 0.0174219301594641737472D+00
    w( 96) = 0.0168019385741038652709D+00
    w( 97) = 0.0161732187295777199419D+00
    w( 98) = 0.0155367755558439824399D+00
    w( 99) = 0.0148936416648151820348D+00
    w(100) = 0.0142448773729167743063D+00
    w(101) = 0.0135915710097655467896D+00
    w(102) = 0.0129348396636073734547D+00
    w(103) = 0.0122758305600827700870D+00
    w(104) = 0.0116157233199551347270D+00
    w(105) = 0.0109557333878379016480D+00
    w(106) = 0.0102971169579563555237D+00
    w(107) = 0.00964117772970253669530D+00
    w(108) = 0.00898927578406413572328D+00
    w(109) = 0.00834283875396815770558D+00
    w(110) = 0.00770337523327974184817D+00
    w(111) = 0.00707248999543355546805D+00
    w(112) = 0.00645190005017573692280D+00
    w(113) = 0.00584344987583563950756D+00
    w(114) = 0.00524912345480885912513D+00
    w(115) = 0.00467105037211432174741D+00
    w(116) = 0.00411150397865469304717D+00
    w(117) = 0.00357289278351729964938D+00
    w(118) = 0.00305775341017553113613D+00
    w(119) = 0.00256876494379402037313D+00
    w(120) = 0.00210881524572663287933D+00
    w(121) = 0.00168114286542146990631D+00
    w(122) = 0.00128952408261041739210D+00
    w(123) = 0.000938369848542381500794D+00
    w(124) = 0.000632607319362633544219D+00
    w(125) = 0.000377746646326984660274D+00
    w(126) = 0.000180739564445388357820D+00
    w(127) = 0.0000505360952078625176247D+00

  else if ( order == 255 ) then

    w(  1) = 0.69379364324108267170D-05
    w(  2) = 0.25157870384280661489D-04
    w(  3) = 0.53275293669780613125D-04
    w(  4) = 0.90372734658751149261D-04
    w(  5) = 0.13575491094922871973D-03
    w(  6) = 0.18887326450650491366D-03
    w(  7) = 0.24921240048299729402D-03
    w(  8) = 0.31630366082226447689D-03
    w(  9) = 0.38974528447328229322D-03
    w( 10) = 0.46918492424785040975D-03
    w( 11) = 0.55429531493037471492D-03
    w( 12) = 0.64476204130572477933D-03
    w( 13) = 0.74028280424450333046D-03
    w( 14) = 0.84057143271072246365D-03
    w( 15) = 0.94536151685852538246D-03
    w( 16) = 0.10544076228633167722D-02
    w( 17) = 0.11674841174299594077D-02
    w( 18) = 0.12843824718970101768D-02
    w( 19) = 0.14049079956551446427D-02
    w( 20) = 0.15288767050877655684D-02
    w( 21) = 0.16561127281544526052D-02
    w( 22) = 0.17864463917586498247D-02
    w( 23) = 0.19197129710138724125D-02
    w( 24) = 0.20557519893273465236D-02
    w( 25) = 0.21944069253638388388D-02
    w( 26) = 0.23355251860571608737D-02
    w( 27) = 0.24789582266575679307D-02
    w( 28) = 0.26245617274044295626D-02
    w( 29) = 0.27721957645934509940D-02
    w( 30) = 0.29217249379178197538D-02
    w( 31) = 0.30730184347025783234D-02
    w( 32) = 0.32259500250878684614D-02
    w( 33) = 0.33803979910869203823D-02
    w( 34) = 0.35362449977167777340D-02
    w( 35) = 0.36933779170256508183D-02
    w( 36) = 0.38516876166398709241D-02
    w( 37) = 0.40110687240750233989D-02
    w( 38) = 0.41714193769840788528D-02
    w( 39) = 0.43326409680929828545D-02
    w( 40) = 0.44946378920320678616D-02
    w( 41) = 0.46573172997568547773D-02
    w( 42) = 0.48205888648512683476D-02
    w( 43) = 0.49843645647655386012D-02
    w( 44) = 0.51485584789781777618D-02
    w( 45) = 0.53130866051870565663D-02
    w( 46) = 0.54778666939189508240D-02
    w( 47) = 0.56428181013844441585D-02
    w( 48) = 0.58078616599775673635D-02
    w( 49) = 0.59729195655081658049D-02
    w( 50) = 0.61379152800413850435D-02
    w( 51) = 0.63027734490857587172D-02
    w( 52) = 0.64674198318036867274D-02
    w( 53) = 0.66317812429018878941D-02
    w( 54) = 0.67957855048827733948D-02
    w( 55) = 0.69593614093904229394D-02
    w( 56) = 0.71224386864583871532D-02
    w( 57) = 0.72849479805538070639D-02
    w( 58) = 0.74468208324075910174D-02
    w( 59) = 0.76079896657190565832D-02
    w( 60) = 0.77683877779219912200D-02
    w( 61) = 0.79279493342948491103D-02
    w( 62) = 0.80866093647888599710D-02
    w( 63) = 0.82443037630328680306D-02
    w( 64) = 0.84009692870519326354D-02
    w( 65) = 0.85565435613076896192D-02
    w( 66) = 0.87109650797320868736D-02
    w( 67) = 0.88641732094824942641D-02
    w( 68) = 0.90161081951956431600D-02
    w( 69) = 0.91667111635607884067D-02
    w( 70) = 0.93159241280693950932D-02
    w( 71) = 0.94636899938300652943D-02
    w( 72) = 0.96099525623638830097D-02
    w( 73) = 0.97546565363174114611D-02
    w( 74) = 0.98977475240487497440D-02
    w( 75) = 0.10039172044056840798D-01
    w( 76) = 0.10178877529236079733D-01
    w( 77) = 0.10316812330947621682D-01
    w( 78) = 0.10452925722906011926D-01
    w( 79) = 0.10587167904885197931D-01
    w( 80) = 0.10719490006251933623D-01
    w( 81) = 0.10849844089337314099D-01
    w( 82) = 0.10978183152658912470D-01
    w( 83) = 0.11104461134006926537D-01
    w( 84) = 0.11228632913408049354D-01
    w( 85) = 0.11350654315980596602D-01
    w( 86) = 0.11470482114693874380D-01
    w( 87) = 0.11588074033043952568D-01
    w( 88) = 0.11703388747657003101D-01
    w( 89) = 0.11816385890830235763D-01
    w( 90) = 0.11927026053019270040D-01
    w( 91) = 0.12035270785279562630D-01
    w( 92) = 0.12141082601668299679D-01
    w( 93) = 0.12244424981611985899D-01
    w( 94) = 0.12345262372243838455D-01
    w( 95) = 0.12443560190714035263D-01
    w( 96) = 0.12539284826474884353D-01
    w( 97) = 0.12632403643542078765D-01
    w( 98) = 0.12722884982732382906D-01
    w( 99) = 0.12810698163877361967D-01
    w(100) = 0.12895813488012114694D-01
    w(101) = 0.12978202239537399286D-01
    w(102) = 0.13057836688353048840D-01
    w(103) = 0.13134690091960152836D-01
    w(104) = 0.13208736697529129966D-01
    w(105) = 0.13279951743930530650D-01
    w(106) = 0.13348311463725179953D-01
    w(107) = 0.13413793085110098513D-01
    w(108) = 0.13476374833816515982D-01
    w(109) = 0.13536035934956213614D-01
    w(110) = 0.13592756614812395910D-01
    w(111) = 0.13646518102571291428D-01
    w(112) = 0.13697302631990716258D-01
    w(113) = 0.13745093443001896632D-01
    w(114) = 0.13789874783240936517D-01
    w(115) = 0.13831631909506428676D-01
    w(116) = 0.13870351089139840997D-01
    w(117) = 0.13906019601325461264D-01
    w(118) = 0.13938625738306850804D-01
    w(119) = 0.13968158806516938516D-01
    w(120) = 0.13994609127619079852D-01
    w(121) = 0.14017968039456608810D-01
    w(122) = 0.14038227896908623303D-01
    w(123) = 0.14055382072649964277D-01
    w(124) = 0.14069424957813575318D-01
    w(125) = 0.14080351962553661325D-01
    w(126) = 0.14088159516508301065D-01
    w(127) = 0.14092845069160408355D-01
    w(128) = 0.14094407090096179347D-01
    w(129) = 0.14092845069160408355D-01
    w(130) = 0.14088159516508301065D-01
    w(131) = 0.14080351962553661325D-01
    w(132) = 0.14069424957813575318D-01
    w(133) = 0.14055382072649964277D-01
    w(134) = 0.14038227896908623303D-01
    w(135) = 0.14017968039456608810D-01
    w(136) = 0.13994609127619079852D-01
    w(137) = 0.13968158806516938516D-01
    w(138) = 0.13938625738306850804D-01
    w(139) = 0.13906019601325461264D-01
    w(140) = 0.13870351089139840997D-01
    w(141) = 0.13831631909506428676D-01
    w(142) = 0.13789874783240936517D-01
    w(143) = 0.13745093443001896632D-01
    w(144) = 0.13697302631990716258D-01
    w(145) = 0.13646518102571291428D-01
    w(146) = 0.13592756614812395910D-01
    w(147) = 0.13536035934956213614D-01
    w(148) = 0.13476374833816515982D-01
    w(149) = 0.13413793085110098513D-01
    w(150) = 0.13348311463725179953D-01
    w(151) = 0.13279951743930530650D-01
    w(152) = 0.13208736697529129966D-01
    w(153) = 0.13134690091960152836D-01
    w(154) = 0.13057836688353048840D-01
    w(155) = 0.12978202239537399286D-01
    w(156) = 0.12895813488012114694D-01
    w(157) = 0.12810698163877361967D-01
    w(158) = 0.12722884982732382906D-01
    w(159) = 0.12632403643542078765D-01
    w(160) = 0.12539284826474884353D-01
    w(161) = 0.12443560190714035263D-01
    w(162) = 0.12345262372243838455D-01
    w(163) = 0.12244424981611985899D-01
    w(164) = 0.12141082601668299679D-01
    w(165) = 0.12035270785279562630D-01
    w(166) = 0.11927026053019270040D-01
    w(167) = 0.11816385890830235763D-01
    w(168) = 0.11703388747657003101D-01
    w(169) = 0.11588074033043952568D-01
    w(170) = 0.11470482114693874380D-01
    w(171) = 0.11350654315980596602D-01
    w(172) = 0.11228632913408049354D-01
    w(173) = 0.11104461134006926537D-01
    w(174) = 0.10978183152658912470D-01
    w(175) = 0.10849844089337314099D-01
    w(176) = 0.10719490006251933623D-01
    w(177) = 0.10587167904885197931D-01
    w(178) = 0.10452925722906011926D-01
    w(179) = 0.10316812330947621682D-01
    w(180) = 0.10178877529236079733D-01
    w(181) = 0.10039172044056840798D-01
    w(182) = 0.98977475240487497440D-02
    w(183) = 0.97546565363174114611D-02
    w(184) = 0.96099525623638830097D-02
    w(185) = 0.94636899938300652943D-02
    w(186) = 0.93159241280693950932D-02
    w(187) = 0.91667111635607884067D-02
    w(188) = 0.90161081951956431600D-02
    w(189) = 0.88641732094824942641D-02
    w(190) = 0.87109650797320868736D-02
    w(191) = 0.85565435613076896192D-02
    w(192) = 0.84009692870519326354D-02
    w(193) = 0.82443037630328680306D-02
    w(194) = 0.80866093647888599710D-02
    w(195) = 0.79279493342948491103D-02
    w(196) = 0.77683877779219912200D-02
    w(197) = 0.76079896657190565832D-02
    w(198) = 0.74468208324075910174D-02
    w(199) = 0.72849479805538070639D-02
    w(200) = 0.71224386864583871532D-02
    w(201) = 0.69593614093904229394D-02
    w(202) = 0.67957855048827733948D-02
    w(203) = 0.66317812429018878941D-02
    w(204) = 0.64674198318036867274D-02
    w(205) = 0.63027734490857587172D-02
    w(206) = 0.61379152800413850435D-02
    w(207) = 0.59729195655081658049D-02
    w(208) = 0.58078616599775673635D-02
    w(209) = 0.56428181013844441585D-02
    w(210) = 0.54778666939189508240D-02
    w(211) = 0.53130866051870565663D-02
    w(212) = 0.51485584789781777618D-02
    w(213) = 0.49843645647655386012D-02
    w(214) = 0.48205888648512683476D-02
    w(215) = 0.46573172997568547773D-02
    w(216) = 0.44946378920320678616D-02
    w(217) = 0.43326409680929828545D-02
    w(218) = 0.41714193769840788528D-02
    w(219) = 0.40110687240750233989D-02
    w(220) = 0.38516876166398709241D-02
    w(221) = 0.36933779170256508183D-02
    w(222) = 0.35362449977167777340D-02
    w(223) = 0.33803979910869203823D-02
    w(224) = 0.32259500250878684614D-02
    w(225) = 0.30730184347025783234D-02
    w(226) = 0.29217249379178197538D-02
    w(227) = 0.27721957645934509940D-02
    w(228) = 0.26245617274044295626D-02
    w(229) = 0.24789582266575679307D-02
    w(230) = 0.23355251860571608737D-02
    w(231) = 0.21944069253638388388D-02
    w(232) = 0.20557519893273465236D-02
    w(233) = 0.19197129710138724125D-02
    w(234) = 0.17864463917586498247D-02
    w(235) = 0.16561127281544526052D-02
    w(236) = 0.15288767050877655684D-02
    w(237) = 0.14049079956551446427D-02
    w(238) = 0.12843824718970101768D-02
    w(239) = 0.11674841174299594077D-02
    w(240) = 0.10544076228633167722D-02
    w(241) = 0.94536151685852538246D-03
    w(242) = 0.84057143271072246365D-03
    w(243) = 0.74028280424450333046D-03
    w(244) = 0.64476204130572477933D-03
    w(245) = 0.55429531493037471492D-03
    w(246) = 0.46918492424785040975D-03
    w(247) = 0.38974528447328229322D-03
    w(248) = 0.31630366082226447689D-03
    w(249) = 0.24921240048299729402D-03
    w(250) = 0.18887326450650491366D-03
    w(251) = 0.13575491094922871973D-03
    w(252) = 0.90372734658751149261D-04
    w(253) = 0.53275293669780613125D-04
    w(254) = 0.25157870384280661489D-04
    w(255) = 0.69379364324108267170D-05

  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GP_WEIGHTS - Fatal error!'
    write ( *, '(a)' ) '  Illegal input value of ORDER.'
    write ( *, '(a)' ) '  Order must be 1, 3, 7, 15, 31, 63, 127 or 255.'
    write ( *, '(a,i8)' ) '  ORDER = ', order
    stop
  end if

  return
end
function i4_log_2 ( i )

!*****************************************************************************80
!
!! I4_LOG_2 returns the integer part of the logarithm base 2 of an I4.
!
!  Discussion:
!
!    For positive I4_LOG_2(I), it should be true that
!      2**I4_LOG_2(X) <= |I| < 2**(I4_LOG_2(I)+1).
!    The special case of I4_LOG_2(0) returns -HUGE().
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Example:
!
!     I  I4_LOG_2
!
!     0  -1
!     1,  0
!     2,  1
!     3,  1
!     4,  2
!     5,  2
!     6,  2
!     7,  2
!     8,  3
!     9,  3
!    10,  3
!   127,  6
!   128,  7
!   129,  7
!
!  Modified:
!
!    20 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the number whose logarithm base 2
!    is desired.
!
!    Output, integer ( kind = 4 ) I4_LOG_2, the integer part of the
!    logarithm base 2 of the absolute value of I.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_abs
  integer ( kind = 4 ) i4_log_2
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647

  if ( i == 0 ) then

    i4_log_2 = - i4_huge

  else

    i4_log_2 = 0

    i_abs = abs ( i )

    do while ( 2 <= i_abs )
      i_abs = i_abs / 2
      i4_log_2 = i4_log_2 + 1
    end do

  end if

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
!        I     J     MOD I4_MODP    Factorization
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
!    08 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the number to be divided.
!
!    Input, integer ( kind = 4 ) J, the number that divides I.
!
!    Output, integer ( kind = 4 ) I4_MODP, the nonnegative remainder when I is
!    divided by J.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ) j
  integer ( kind = 4 ) value

  if ( j == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_MODP - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal divisor J = ', j
    stop
  end if

  value = mod ( i, j )

  if ( value < 0 ) then
    value = value + abs ( j )
  end if

  i4_modp = value

  return
end
subroutine index_level_own ( level, level_max, dim_num, point_num, &
  grid_index, grid_base, grid_level )

!*****************************************************************************80
!
!! INDEX_LEVEL_OWN: determine first level at which given index is generated.
!
!  Discussion:
!
!    We are constructing a sparse grid based on a 1D OWN rule (Gauss Legendre
!    or Gauss Hermite).  The grid is built up of product grids,
!    with a characteristic LEVEL.
!
!    We are concerned with identifying points in this product grid which
!    have actually been generated previously, on a lower value of LEVEL.
!
!    This routine determines the lowest value of LEVEL at which each of
!    the input points would be generated.
!
!    In 1D, given LEVEL, the number of points is ORDER = 2**(LEVEL+1) + 1,
!    (except that LEVEL = 0 implies ORDER = 1!), the BASE is (ORDER-1)/2,
!    and the point INDEX values range from -BASE to +BASE.
!
!    The values of INDEX and BASE allow us to determine the abstract
!    properties of the point.  In particular, if INDEX is 0, the corresponding
!    abscissa is 0, the special "nested" value we need to take care of.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LEVEL, the level at which these points were
!    generated.  LEVEL_MIN <= LEVEL <= LEVEL_MAX.
!
!    Input, integer ( kind = 4 ) LEVEL_MAX, the maximum level.
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points to be tested.
!
!    Input, integer ( kind = 4 ) GRID_INDEX(DIM_NUM,POINT_NUM), the indices of
!    the points to be tested.
!
!    Input, integer ( kind = 4 ) GRID_BASE(DIM_NUM), the "base", which is
!    essentially the denominator of the index.
!
!    Output, integer ( kind = 4 ) GRID_LEVEL(POINT_NUM), the value of LEVEL at
!    which the point would first be generated.  This will be the same as
!    the input value of LEVEL, unless the point has an INDEX of 0 and
!    a corresponding BASE that is NOT zero.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) grid_base(dim_num)
  integer ( kind = 4 ) grid_index(dim_num,point_num)
  integer ( kind = 4 ) grid_level(point_num)
  integer ( kind = 4 ) level
  integer ( kind = 4 ) level_max
  integer ( kind = 4 ) level_min
  integer ( kind = 4 ) point

  if ( dim_num == 1 ) then
    level_min = level_max
  else
    level_min = 0
  end if
!
!  If a point has a DIM-th component whose INDEX is 0, then the
!  value of LEVEL at which this point would first be generated is
!  less than LEVEL, unless the DIM-th component of GRID_BASE is 0.
!
  do point = 1, point_num

    grid_level(point) = max ( level, level_min )

    do dim = 1, dim_num
      if ( grid_index(dim,point) == 0 ) then
        grid_level(point) = &
          max ( grid_level(point) - grid_base(dim), level_min )
      end if
    end do

  end do

  return
end
function index_to_level_closed ( dim_num, t, order, level_max )

!*****************************************************************************80
!
!! INDEX_TO_LEVEL_CLOSED determines the level of a point given its index.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) T(DIM_NUM), the grid indices of a point
!    in a 1D closed rule.  0 <= T(I) <= ORDER.
!
!    Input, integer ( kind = 4 ) ORDER, the order of the rule.
!
!    Input, integer ( kind = 4 ) LEVEL_MAX, the level with respect to which the
!    index applies.
!
!    Output, integer ( kind = 4 ) INDEX_TO_LEVEL_CLOSED, the first level on
!    which the point associated with the given index will appear.
!
  implicit none

  integer ( kind = 4 ) dim_num

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ) index_to_level_closed
  integer ( kind = 4 ) level
  integer ( kind = 4 ) level_max
  integer ( kind = 4 ) order
  integer ( kind = 4 ) s
  integer ( kind = 4 ) t(dim_num)

  index_to_level_closed = 0

  do dim = 1, dim_num

    s = i4_modp ( t(dim), order )

    if ( s == 0 ) then

      level = 0

    else

      level = level_max

      do while ( mod ( s, 2 ) == 0 )
        s = s / 2
        level = level - 1
      end do

    end if

    if ( level == 0 ) then
      level = 1
    else if ( level == 1 ) then
      level = 0
    end if

    index_to_level_closed = index_to_level_closed + level

  end do

  return
end
function index_to_level_open ( dim_num, t, order, level_max )

!*****************************************************************************80
!
!! INDEX_TO_LEVEL_OPEN determines the level of a point given its index.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) T(DIM_NUM), the grid index of a point.
!
!    Input, integer ( kind = 4 ) ORDER, the order.
!
!    Input, integer ( kind = 4 ) LEVEL_MAX, the level with respect to which the
!    index applies.
!
!    Output, integer ( kind = 4 ) INDEX_TO_LEVEL_OPEN, the first level on which
!    the point associated with the given index will appear.
!
  implicit none

  integer ( kind = 4 ) dim_num

  integer ( kind = 4 ) dim
  integer ( kind = 4 ), parameter :: i4_2 = 2
  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ) index_to_level_open
  integer ( kind = 4 ) level
  integer ( kind = 4 ) level_max
  integer ( kind = 4 ) order
  integer ( kind = 4 ) s
  integer ( kind = 4 ) t(dim_num)

  index_to_level_open = 0

  do dim = 1, dim_num

    s = t(dim)

    s = i4_modp ( s, order )

    if ( s == 0 ) then

      level = 0

    else

      level = level_max

      do while ( mod ( s, i4_2 ) == 0 )
        s = s / 2
        level = level - 1
      end do

    end if

    if ( level == 0 ) then
      level = 1
    else if ( level == 1 ) then
      level = 0
    end if

    index_to_level_open = index_to_level_open + level

  end do

  return
end
subroutine level_to_order_closed ( dim_num, level, order )

!*****************************************************************************80
!
!! LEVEL_TO_ORDER_CLOSED converts a level to an order for closed rules.
!
!  Discussion:
!
!    Sparse grids can naturally be nested.  A natural scheme is to use
!    a series of one-dimensional rules arranged in a series of "levels"
!    whose order roughly doubles with each step.
!
!    The arrangement described here works for the Clenshaw Curtis rule.
!
!    The idea is that we start with LEVEL = 0, ORDER = 1 indicating the single
!    point at the center, and for all values afterwards, we use the
!    relationship
!
!      ORDER = 2**LEVEL + 1
!
!    The following table shows how the growth will occur:
!
!    Level    Order
!
!    0          1
!    1          3 =  2 + 1
!    2          5 =  4 + 1
!    3          9 =  8 + 1
!    4         17 = 16 + 1
!    5         33 = 32 + 1
!
!    For the Clenshaw Curtis rules, the point growth
!    is nested.  If we have ORDER points on a particular LEVEL, the next
!    level includes all these old points, plus ORDER-1 new points, formed
!    in the gaps between successive pairs of old points.
!
!    Level    Order = New + Old
!
!    0          1   =  1  +  0
!    1          3   =  2  +  1
!    2          5   =  2  +  3
!    3          9   =  4  +  5
!    4         17   =  8  +  9
!    5         33   = 16  + 17
!
!    In this routine, we assume that a vector of levels is given,
!    and the corresponding orders are desired.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) LEVEL(DIM_NUM), the nesting levels of the
!    1D rules.
!
!    Output, integer ( kind = 4 ) ORDER(DIM_NUM), the order (number of points)
!    of the 1D rules.
!
  implicit none

  integer ( kind = 4 ) dim_num

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) level(dim_num)
  integer ( kind = 4 ) order(dim_num)

  do dim = 1, dim_num

    if ( level(dim) < 0 ) then
      order(dim) = -1
    else if ( level(dim) == 0 ) then
      order(dim) = 1
    else
      order(dim) = ( 2**level(dim) ) + 1
    end if

  end do

  return
end
subroutine level_to_order_open ( dim_num, level, order )

!*****************************************************************************80
!
!! LEVEL_TO_ORDER_OPEN converts a level to an order for open rules.
!
!  Discussion:
!
!    Sparse grids can naturally be nested.  A natural scheme is to use
!    a series of one-dimensional rules arranged in a series of "levels"
!    whose order roughly doubles with each step.
!
!    The arrangement described here works naturally for the Fejer Type 1,
!    Fejer Type 2, and Gauss Patterson rules.  It also can be used, partially,
!    to describe the growth of Gauss Legendre rules.
!
!    The idea is that we start with LEVEL = 0, ORDER = 1 indicating the single
!    point at the center, and for all values afterwards, we use the
!    relationship
!
!      ORDER = 2^(LEVEL+1) - 1.
!
!    The following table shows how the growth will occur:
!
!    Level    Order
!
!    0          1
!    1          3 =  4 - 1
!    2          7 =  8 - 1
!    3         15 = 16 - 1
!    4         31 = 32 - 1
!    5         63 = 64 - 1
!
!    For the Fejer Type 1, Fejer Type 2, and Gauss Patterson rules, the point
!    growth is nested.  If we have ORDER points on a particular LEVEL, the next
!    level includes all these old points, plus ORDER+1 new points, formed in the
!    gaps between successive pairs of old points plus an extra point at each
!    end.
!
!    Level    Order = New + Old
!
!    0          1   =  1  +  0
!    1          3   =  2  +  1
!    2          7   =  4  +  3
!    3         15   =  8  +  7
!    4         31   = 16  + 15
!    5         63   = 32  + 31
!
!    If we use a series of Gauss Legendre rules, then there is almost no
!    nesting, except that the central point is shared.  If we insist on
!    producing a comparable series of such points, then the "nesting" behavior
!    is as follows:
!
!    Level    Order = New + Old
!
!    0          1   =  1  +  0
!    1          3   =  2  +  1
!    2          7   =  6  +  1
!    3         15   = 14  +  1
!    4         31   = 30  +  1
!    5         63   = 62  +  1
!
!    Moreover, if we consider ALL the points used in such a set of "nested"
!    Gauss Legendre rules, then we must sum the "NEW" column, and we see that
!    we get roughly twice as many points as for the truly nested rules.
!
!    In this routine, we assume that a vector of levels is given,
!    and the corresponding orders are desired.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) LEVEL(DIM_NUM), the nesting levels of the
!    1D rules.
!
!    Output, integer ( kind = 4 ) ORDER(DIM_NUM), the order (number of points)
!    of the 1D rules.
!
  implicit none

  integer ( kind = 4 ) dim_num

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) level(dim_num)
  integer ( kind = 4 ) order(dim_num)

  do dim = 1, dim_num

    if ( level(dim) < 0 ) then
      order(dim) = -1
    else
      order(dim) = 2**( level(dim) + 1 ) - 1
    end if

  end do

  return
end
subroutine levels_index ( dim_num, level_max, rule, point_num, grid_index, &
  grid_base )

!*****************************************************************************80
!
!! LEVELS_INDEX indexes a sparse grid.
!
!  Discussion:
!
!    The sparse grid is the logical sum of product grids with total LEVEL
!    between LEVEL_MIN and LEVEL_MAX.
!
!    The necessary dimensions of GRID_INDEX can be determined by
!    calling LEVELS_INDEX_SIZE first.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 March 2008
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
!    Input, integer ( kind = 4 ) RULE, the index of the rule.
!    1, "CC", Clenshaw Curtis Closed Fully Nested rule.
!    2, "F1", Fejer 1 Open Fully Nested rule.
!    3, "F2", Fejer 2 Open Fully Nested rule.
!    4, "GP", Gauss Patterson Open Fully Nested rule.
!    5, "GL", Gauss Legendre Open Weakly Nested rule.
!    6, "GH", Gauss Hermite Open Weakly Nested rule.
!    7, "LG", Gauss Laguerre Open Non Nested rule.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the total number of points
!    in the grids.
!
!    Output, integer ( kind = 4 ) GRID_INDEX(DIM_NUM,POINT_NUM), a list of
!    point indices, representing a subset of the product grid of level
!    LEVEL_MAX, representing (exactly once) each point that will show up in a
!    sparse grid of level LEVEL_MAX.
!
!    Output, integer ( kind = 4 ) GRID_BASE(DIM_NUM,POINT_NUM), a list of
!    the orders of the rules associated with each point and dimension.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) grid_base(dim_num,point_num)
  integer ( kind = 4 ) grid_index(dim_num,point_num)
  integer ( kind = 4 ) level_max
  integer ( kind = 4 ) rule

  if ( rule == 1 ) then
    call levels_index_cfn ( dim_num, level_max, point_num, grid_index, &
      grid_base )
  else if ( 2 <= rule .and. rule <= 4 ) then
    call levels_index_ofn ( dim_num, level_max, point_num, grid_index, &
      grid_base )
  else if ( 5 <= rule .and. rule <= 6 ) then
    call levels_index_own ( dim_num, level_max, point_num, grid_index, &
      grid_base )
  else if ( 7 == rule ) then
    call levels_index_onn ( dim_num, level_max, point_num, grid_index, &
      grid_base )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LEVELS_INDEX - Fatal error!'
    write ( *, '(a,i8)' ) '  Unrecognized rule number = ', rule
    stop
  end if

  return
end
subroutine levels_index_cfn ( dim_num, level_max, point_num, grid_index, &
  grid_base )

!*****************************************************************************80
!
!! LEVELS_INDEX_CFN indexes a sparse grid made from CFN 1D rules.
!
!  Discussion:
!
!    The sparse grid is presumed to have been created from products
!    of CLOSED FULLY NESTED 1D quadrature rules.
!
!    CFN rules include Clenshaw Curtis rules.
!
!    The sparse grid is the logical sum of product grids with total LEVEL
!    between LEVEL_MIN and LEVEL_MAX.
!
!    The necessary dimensions of GRID_INDEX can be determined by
!    calling LEVELS_INDEX_SIZE_CFN first.
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
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) LEVEL_MAX, the maximum value of LEVEL.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the total number of points
!    in the grids.
!
!    Output, integer ( kind = 4 ) GRID_INDEX(DIM_NUM,POINT_NUM), a list of
!    point indices, representing a subset of the product grid of level
!    LEVEL_MAX, representing (exactly once) each point that will show up in a
!    sparse grid of level LEVEL_MAX.
!
!    Output, integer ( kind = 4 ) GRID_BASE(DIM_NUM,POINT_NUM), a list of
!    the orders of the rules associated with each point and dimension.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) grid_base(dim_num,point_num)
  integer ( kind = 4 ) grid_index(dim_num,point_num)
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: grid_index2
  integer ( kind = 4 ), allocatable, dimension ( : ) :: grid_level
  integer ( kind = 4 ) h
  integer ( kind = 4 ) level
  integer ( kind = 4 ), dimension ( dim_num ) :: level_1d
  integer ( kind = 4 ) level_max
  logical more
  integer ( kind = 4 ), dimension ( dim_num ) :: order_1d
  integer ( kind = 4 ) order_nd
  integer ( kind = 4 ) point
  integer ( kind = 4 ) point_num2
  integer ( kind = 4 ) t
!
!  The outer loop generates LEVELs from 0 to LEVEL_MAX.
!
  point_num2 = 0

  do level = 0, level_max
!
!  The middle loop generates the next partition LEVEL_1D(1:DIM_NUM)
!  that adds up to LEVEL.
!
    more = .false.
    h = 0
    t = 0

    do

      call comp_next ( level, dim_num, level_1d, more, h, t )
!
!  Transform each 1D level to a corresponding 1D order.
!
      call level_to_order_closed ( dim_num, level_1d, order_1d )
!
!  The product of the 1D orders gives us the number of points in this grid.
!
      order_nd = product ( order_1d(1:dim_num) )

      allocate ( grid_index2(1:dim_num,1:order_nd) )
      allocate ( grid_level(1:order_nd) )
!
!  The inner (hidden) loop generates all points corresponding to given grid.
!
      call multigrid_index_cfn ( dim_num, order_1d, order_nd, grid_index2 )
!
!  Adjust these grid indices to reflect LEVEL_MAX.
!
      call multigrid_scale_closed ( dim_num, order_nd, level_max, level_1d, &
        grid_index2 )
!
!  Determine the first level of appearance of each of the points.
!
      call abscissa_level_closed_nd ( level_max, dim_num, order_nd, &
        grid_index2, grid_level )
!
!  Only keep those points which first appear on this level.
!
      do point = 1, order_nd

        if ( grid_level(point) == level ) then

          point_num2 = point_num2 + 1

          if ( point_num < point_num2 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'LEVELS_INDEX_CFN - Fatal error!'
            write ( *, '(a,i8)' ) &
            '  Exceeding maximum point index POINT_NUM = ', point_num
            stop
          end if

          grid_base(1:dim_num,point_num2) = order_1d(1:dim_num)
          grid_index(1:dim_num,point_num2) = grid_index2(1:dim_num,point)

        end if

      end do

      deallocate ( grid_index2 )
      deallocate ( grid_level )

      if ( .not. more ) then
        exit
      end if

    end do

  end do

  if ( point_num2 < point_num ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LEVELS_INDEX_CFN - Fatal error!'
    write ( *, '(a,i8)' ) '  Set fewer points than POINT_NUM = ', point_num
    stop
  end if

  return
end
subroutine levels_index_ofn ( dim_num, level_max, point_num, grid_index, &
  grid_base )

!*****************************************************************************80
!
!! LEVELS_INDEX_OFN indexes a sparse grid made from OFN 1D rules.
!
!  Discussion:
!
!    The sparse grid is presumed to have been created from products
!    of OPEN FULLY NESTED 1D quadrature rules.
!
!    OFN rules include Fejer 1, Fejer 2, and Gauss Patterson rules.
!
!    The sparse grid is the logical sum of product grids with total LEVEL
!    between LEVEL_MIN and LEVEL_MAX.
!
!    The necessary dimensions of GRID_INDEX can be determined by
!    calling LEVELS_INDEX_SIZE_OFN first.
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
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) LEVEL_MAX, the maximum value of LEVEL.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the total number of points in
!    the grids.
!
!    Output, integer ( kind = 4 ) GRID_INDEX(DIM_NUM,POINT_NUM), a list of
!    point indices, representing a subset of the product grid of level
!    LEVEL_MAX, representing (exactly once) each point that will show up in a
!    sparse grid of level LEVEL_MAX.
!
!    Output, integer ( kind = 4 ) GRID_BASE(DIM_NUM,POINT_NUM), a list of
!    the orders of the rules associated with each point and dimension.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) grid_base(dim_num,point_num)
  integer ( kind = 4 ) grid_index(dim_num,point_num)
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: grid_index2
  integer ( kind = 4 ) h
  integer ( kind = 4 ) level
  integer ( kind = 4 ), dimension ( dim_num ) :: level_1d
  integer ( kind = 4 ) level_max
  logical more
  integer ( kind = 4 ), dimension ( dim_num ) :: order_1d
  integer ( kind = 4 ) order_nd
  integer ( kind = 4 ) point
  integer ( kind = 4 ) point_num2
  integer ( kind = 4 ) t
!
!  The outer loop generates LEVELs from 0 to LEVEL_MAX.
!
  point_num2 = 0

  do level = 0, level_max
!
!  The middle loop generates the next partition, LEVEL_1D(1:DIM_NUM),
!  that adds up to LEVEL.
!
    more = .false.
    h = 0
    t = 0

    do

      call comp_next ( level, dim_num, level_1d, more, h, t )
!
!  Transform each 1D level to a corresponding 1D order,
!  ORDER = 2^(LEVEL+1)-1.
!
      call level_to_order_open ( dim_num, level_1d, order_1d )
!
!  The product of the 1D orders gives us the number of points in this grid.
!
      order_nd = product ( order_1d(1:dim_num) )

      allocate ( grid_index2(1:dim_num,1:order_nd) )
!
!  The inner (hidden) loop generates all points corresponding to this grid.
!
      call multigrid_index_ofn ( dim_num, order_1d, order_nd, grid_index2 )
!
!  Only keep those points which first appear on this level.
!  If you keep a point, it is necessary to rescale each of its components
!  so that we save the coordinates as they apply on the final grid.
!
      do point = 1, order_nd

        if ( all ( mod ( grid_index2(1:dim_num,point), 2 ) == 1 ) ) then

          point_num2 = point_num2 + 1

          if ( point_num < point_num2 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'LEVELS_INDEX_OFN - Fatal error!'
            write ( *, '(a,i8)' ) &
            '  Exceeding maximum point index POINT_NUM = ', point_num
            stop
          end if

          grid_base(1:dim_num,point_num2) = order_1d(1:dim_num)

          do dim = 1, dim_num
            grid_index(dim,point_num2) = &
              2**( level_max - level_1d(dim) ) * grid_index2(dim,point)
          end do

        end if

      end do

      deallocate ( grid_index2 )

      if ( .not. more ) then
        exit
      end if

    end do

  end do

  if ( point_num2 < point_num ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LEVELS_INDEX_OFN - Fatal error!'
    write ( *, '(a,i8)' ) '  Set fewer points than POINT_NUM = ', point_num
    stop
  end if

  return
end
subroutine levels_index_onn ( dim_num, level_max, point_num, grid_index, &
  grid_base )

!*****************************************************************************80
!
!! LEVELS_INDEX_ONN indexes a sparse grid made from ONN 1D rules.
!
!  Discussion:
!
!    The sparse grid is presumed to have been created from products
!    of OPEN NON NESTED 1D quadrature rules.
!
!    ONN rules include Gauss Laguerre.
!
!    The sparse grid is the logical sum of product grids with total LEVEL
!    between LEVEL_MIN and LEVEL_MAX.
!
!    The necessary dimensions of GRID_INDEX can be determined by
!    calling LEVELS_INDEX_SIZE_ONN first.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 March 2008
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
!    Input, integer ( kind = 4 ) POINT_NUM, the total number of points in
!    the grids.
!
!    Output, integer ( kind = 4 ) GRID_INDEX(DIM_NUM,POINT_NUM), a list of
!    point indices, representing a subset of the product grid of level
!    LEVEL_MAX, representing (exactly once) each point that will show up in a
!    sparse grid of level LEVEL_MAX.
!
!    Output, integer ( kind = 4 ) GRID_BASE(DIM_NUM,POINT_NUM), a list of
!    the orders of the rules associated with each point and dimension.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) grid_base(dim_num,point_num)
  integer ( kind = 4 ), dimension ( dim_num ) :: grid_base2
  integer ( kind = 4 ) grid_index(dim_num,point_num)
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: grid_index2
  integer ( kind = 4 ) h
  integer ( kind = 4 ) level
  integer ( kind = 4 ), dimension ( dim_num ) :: level_1d
  integer ( kind = 4 ) level_max
  integer ( kind = 4 ) level_min
  logical more
  integer ( kind = 4 ), dimension ( dim_num ) :: order_1d
  integer ( kind = 4 ) order_nd
  integer ( kind = 4 ) point
  integer ( kind = 4 ) point_num2
  integer ( kind = 4 ) t
!
!  The outer loop generates LEVELs from LEVEL_MIN to LEVEL_MAX.
!
  point_num2 = 0

  level_min = max ( 0, level_max + 1 - dim_num )

  do level = level_min, level_max
!
!  The middle loop generates the next partition LEVEL_1D(1:DIM_NUM)
!  that adds up to LEVEL.
!
    more = .false.
    h = 0
    t = 0

    do

      call comp_next ( level, dim_num, level_1d, more, h, t )
!
!  Transform each 1D level to a corresponding 1D order.
!  The relationship is the same as for other OPEN rules.
!  The GL rule differs from the other OPEN rules only in the nesting behavior.
!
      call level_to_order_open ( dim_num, level_1d, order_1d )

      grid_base2(1:dim_num) = order_1d(1:dim_num)
!
!  The product of the 1D orders gives us the number of points in this subgrid.
!
      order_nd = product ( order_1d(1:dim_num) )

      allocate ( grid_index2(1:dim_num,1:order_nd) )
!
!  The inner (hidden) loop generates all points corresponding to given grid.
!  The grid indices will be between 1 and ORDER_1D(DIM).
!
      call multigrid_index_onn ( dim_num, order_1d, order_nd, grid_index2 )

      do point = 1, order_nd

        point_num2 = point_num2 + 1

        if ( point_num < point_num2 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'LEVELS_INDEX_ONN - Fatal error!'
          write ( *, '(a,i8)' ) &
          '  Exceeding maximum point index POINT_NUM = ', point_num
          stop
        end if

        grid_index(1:dim_num,point_num2) = grid_index2(1:dim_num,point)
        grid_base(1:dim_num,point_num2) = grid_base2(1:dim_num)

      end do

      deallocate ( grid_index2 )

      if ( .not. more ) then
        exit
      end if

    end do

  end do

  if ( point_num2 < point_num ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LEVELS_INDEX_ONN - Fatal error!'
    write ( *, '(a,i8)' ) '  Set fewer points than POINT_NUM = ', point_num
    stop
  end if

  return
end
subroutine levels_index_own ( dim_num, level_max, point_num, grid_index, &
  grid_base )

!*****************************************************************************80
!
!! LEVELS_INDEX_OWN indexes a sparse grid made from OWN 1D rules.
!
!  Discussion:
!
!    The sparse grid is presumed to have been created from products
!    of OPEN WEAKLY NESTED 1D quadrature rules.
!
!    OWN rules include Gauss Hermite and Gauss Legendre.
!
!    The sparse grid is the logical sum of product grids with total LEVEL
!    between LEVEL_MIN and LEVEL_MAX.
!
!    The necessary dimensions of GRID_INDEX can be determined by
!    calling LEVELS_INDEX_SIZE_OWN first.
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
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) LEVEL_MAX, the maximum value of LEVEL.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the total number of points in
!    the grids.
!
!    Output, integer ( kind = 4 ) GRID_INDEX(DIM_NUM,POINT_NUM), a list of
!    point indices, representing a subset of the product grid of level
!    LEVEL_MAX, representing (exactly once) each point that will show up in a
!    sparse grid of level LEVEL_MAX.
!
!    Output, integer ( kind = 4 ) GRID_BASE(DIM_NUM,POINT_NUM), a list of
!    the orders of the rules associated with each point and dimension.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) grid_base(dim_num,point_num)
  integer ( kind = 4 ), dimension ( dim_num ) :: grid_base2
  integer ( kind = 4 ) grid_index(dim_num,point_num)
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: grid_index2
  integer ( kind = 4 ), allocatable, dimension ( : ) :: grid_level
  integer ( kind = 4 ) h
  integer ( kind = 4 ) level
  integer ( kind = 4 ), dimension ( dim_num ) :: level_1d
  integer ( kind = 4 ) level_max
  integer ( kind = 4 ) level_min
  logical more
  integer ( kind = 4 ), dimension ( dim_num ) :: order_1d
  integer ( kind = 4 ) order_nd
  integer ( kind = 4 ) point
  integer ( kind = 4 ) point_num2
  integer ( kind = 4 ) t
!
!  The outer loop generates LEVELs from LEVEL_MIN to LEVEL_MAX.
!
  point_num2 = 0

  if ( dim_num == 1 ) then
    level_min = level_max
  else
    level_min = 0
  end if

  do level = level_min, level_max
!
!  The middle loop generates the next partition LEVEL_1D(1:DIM_NUM)
!  that adds up to LEVEL.
!
    more = .false.
    h = 0
    t = 0

    do

      call comp_next ( level, dim_num, level_1d, more, h, t )
!
!  Transform each 1D level to a corresponding 1D order.
!  The relationship is the same as for other OPEN rules.
!  The GL rule differs from the other OPEN rules only in the nesting behavior.
!
      call level_to_order_open ( dim_num, level_1d, order_1d )

      grid_base2(1:dim_num) = ( order_1d(1:dim_num) - 1 ) / 2
!
!  The product of the 1D orders gives us the number of points in this subgrid.
!
      order_nd = product ( order_1d(1:dim_num) )

      allocate ( grid_index2(1:dim_num,1:order_nd) )
      allocate ( grid_level(1:order_nd) )
!
!  The inner (hidden) loop generates all points corresponding to given grid.
!  The grid indices will be between -M to +M, where 2*M + 1 = ORDER_1D(DIM).
!
      call multigrid_index_own ( dim_num, order_1d, order_nd, grid_index2 )
!
!  Determine the first level of appearance of each of the points.
!  This allows us to flag certain points as being repeats of points
!  generated on a grid of lower level.
!
!  This is SLIGHTLY tricky.
!
      call index_level_own ( level, level_max, dim_num, order_nd, &
        grid_index2, grid_base2, grid_level )
!
!  Only keep those points which first appear on this level.
!
      do point = 1, order_nd

        if ( grid_level(point) == level ) then

          point_num2 = point_num2 + 1

          if ( point_num < point_num2 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'LEVELS_INDEX_OWN - Fatal error!'
            write ( *, '(a,i8)' ) &
            '  Exceeding maximum point index POINT_NUM = ', point_num
            stop
          end if

          grid_index(1:dim_num,point_num2) = grid_index2(1:dim_num,point)
          grid_base(1:dim_num,point_num2) = grid_base2(1:dim_num)

        end if

      end do

      deallocate ( grid_index2 )
      deallocate ( grid_level )

      if ( .not. more ) then
        exit
      end if

    end do

  end do

  if ( point_num2 < point_num ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LEVELS_INDEX_OWN - Fatal error!'
    write ( *, '(a,i8)' ) '  Set fewer points than POINT_NUM = ', point_num
    stop
  end if

  return
end
subroutine levels_index_size ( dim_num, level_max, rule, point_num )

!*****************************************************************************80
!
!! LEVELS_INDEX_SIZE sizes a sparse grid.
!
!  Discussion:
!
!    The sparse grid is the logical sum of product grids with total LEVEL
!    between LEVEL_MIN and LEVEL_MAX.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 March 2008
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
!    Input, integer ( kind = 4 ) RULE, the index of the rule.
!    1, "CC", Clenshaw Curtis Closed Fully Nested rule.
!    2, "F1", Fejer 1 Open Fully Nested rule.
!    3, "F2", Fejer 2 Open Fully Nested rule.
!    4, "GP", Gauss Patterson Open Fully Nested rule.
!    5, "GL", Gauss Legendre Open Weakly Nested rule.
!    6, "GH", Gauss Hermite Open Weakly Nested rule.
!    7, "LG", Gauss Laguerre Open Non Nested rule.
!
!    Output, integer ( kind = 4 ) POINT_NUM, the total number of unique
!    points in the grids.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) level_max
  integer ( kind = 4 ) point_num
  integer ( kind = 4 ) rule

  if ( rule == 1 ) then
    call sparse_grid_cc_size ( dim_num, level_max, point_num )
  else if ( 2 <= rule .and. rule <= 4 ) then
    call sparse_grid_ofn_size ( dim_num, level_max, point_num )
  else if ( 5 <= rule .and. rule <= 6 ) then
    call levels_index_size_own ( dim_num, level_max, point_num )
  else if ( 7 == rule ) then
    call levels_index_size_onn ( dim_num, level_max, point_num )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LEVELS_INDEX_SIZE - Fatal error!'
    write ( *, '(a,i8)' ) '  Unrecognized value of RULE = ', rule
    stop
  end if

  return
end
subroutine levels_index_size_onn ( dim_num, level_max, point_num )

!*****************************************************************************80
!
!! LEVELS_INDEX_SIZE_ONN sizes a sparse grid made from ONN 1D rules.
!
!  Discussion:
!
!    The sparse grid is presumed to have been created from products
!    of OPEN NON-NESTED 1D quadrature rules.
!
!    ONN rules include Gauss Laguerre.
!
!    The sparse grid is the logical sum of product grids with total LEVEL
!    between LEVEL_MIN and LEVEL_MAX.
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
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) LEVEL_MAX, the maximum value of LEVEL.
!
!    Output, integer ( kind = 4 ) POINT_NUM, the number of points in the grid.
!
  implicit none

  integer ( kind = 4 ) dim_num

  integer ( kind = 4 ) h
  integer ( kind = 4 ) level
  integer ( kind = 4 ) level_1d(dim_num)
  integer ( kind = 4 ) level_max
  integer ( kind = 4 ) level_min
  logical              more
  integer ( kind = 4 ) order_1d(dim_num)
  integer ( kind = 4 ) point_num
  integer ( kind = 4 ) t
!
!  Special case.
!
  if ( level_max == 0 ) then
    point_num = 1
    return
  end if
!
!  The outer loop generates LEVELs from LEVEL_MIN to LEVEL_MAX.
!
  level_min = max ( 0, level_max + 1 - dim_num )

  point_num = 0

  do level = level_min, level_max
!
!  The middle loop generates the next partition that adds up to LEVEL.
!
    more = .false.
    h = 0
    t = 0

    do

      call comp_next ( level, dim_num, level_1d, more, h, t )
!
!  Transform each 1D level to a corresponding 1D order.
!
      call level_to_order_open ( dim_num, level_1d, order_1d )

      point_num = point_num + product ( order_1d(1:dim_num) )

      if ( .not. more ) then
        exit
      end if

    end do
  end do

  return
end
subroutine levels_index_size_own ( dim_num, level_max, point_num )

!*****************************************************************************80
!
!! LEVELS_INDEX_SIZE_OWN sizes a sparse grid made from OWN 1D rules.
!
!  Discussion:
!
!    The sparse grid is presumed to have been created from products
!    of OPEN WEAKLY NESTED 1D quadrature rules.
!
!    OWN rules include Gauss Hermite and Gauss Legendre.
!
!    The sparse grid is the logical sum of product grids with total LEVEL
!    between LEVEL_MIN and LEVEL_MAX.
!
!    Oddly enough, in order to count the number of points, we will
!    behave as though LEVEL_MIN was zero.  This is because our computation
!    concentrates on throwing away all points generated at lower levels,
!    but, in fact, if we start at a nonzero level, we need to include
!    on that level all the points that would have been generated on lower
!    levels.
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
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) LEVEL_MAX, the maximum value of LEVEL.
!
!    Output, integer ( kind = 4 ) POINT_NUM, the number of points in the grid.
!
  implicit none

  integer ( kind = 4 ) dim_num

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) h
  integer ( kind = 4 ) level
  integer ( kind = 4 ) level_1d(dim_num)
  integer ( kind = 4 ) level_max
  integer ( kind = 4 ) level_min
  logical              more
  integer ( kind = 4 ) order_1d(dim_num)
  integer ( kind = 4 ) point_num
  integer ( kind = 4 ) t
!
!  Special case.
!
  if ( level_max == 0 ) then
    point_num = 1
    return
  end if
!
!  The outer loop generates LEVELs from LEVEL_MIN to LEVEL_MAX.
!
!  The normal definition of LEVEL_MIN:
!
!   level_min = max ( 0, level_max + 1 - dim_num )
!
!  Our somewhat artificial temporary local definition of LEVEL_MIN:
!
  if ( dim_num == 1 ) then
    level_min = level_max
    point_num = 1
  else
    level_min = 0
    point_num = 0
  end if

  do level = level_min, level_max
!
!  The middle loop generates the next partition that adds up to LEVEL.
!
    more = .false.
    h = 0
    t = 0

    do

      call comp_next ( level, dim_num, level_1d, more, h, t )
!
!  Transform each 1D level to a corresponding 1D order.
!
      call level_to_order_open ( dim_num, level_1d, order_1d )

      do dim = 1, dim_num
!
!  Account for the repetition of the center point.
!
        if ( 1 < order_1d(dim) ) then
          order_1d(dim) = order_1d(dim) - 1
        end if

      end do

      point_num = point_num + product ( order_1d(1:dim_num) )

      if ( .not. more ) then
        exit
      end if

    end do
  end do

  return
end
subroutine lg_abscissa ( dim_num, point_num, grid_index, grid_base, &
  grid_point )

!*****************************************************************************80
!
!! LG_ABSCISSA sets abscissas for multidimensional Gauss Laguerre quadrature.
!
!  Discussion:
!
!    The "nesting" as it occurs for Gauss Laguerre sparse grids simply
!    involves the use of a specified set of permissible orders for the
!    rule.
!
!    The X array lists the (complete) Gauss Laguerre abscissas for rules
!    of order 1, 3, 7, 15, 31, 63, and 127 in order.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, integer ( kind = 4 ) GRID_INDEX(DIM_NUM,POINT_NUM), for each
!    point and dimension, the index of the Gauss Laguerre abscissa.
!
!    Input, integer ( kind = 4 ) GRID_BASE(DIM_NUM), the "base" of the
!    Gauss Laguerre rule being used in each dimension.
!
!    Output, real ( kind = 8 ) GRID_POINT(DIM_NUM), the grid points of
!    Gauss Laguerre abscissas.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) grid_base(dim_num)
  integer ( kind = 4 ) grid_index(dim_num,point_num)
  real ( kind = 8 ) grid_point(dim_num,point_num)
  integer ( kind = 4 ) i4_log_2
  integer ( kind = 4 ) level
  integer ( kind = 4 ) point
  integer ( kind = 4 ) pointer
  integer ( kind = 4 ), dimension ( 0:7 ) :: skip = (/ &
    0, 1, 4, 11, 26, 57, 120, 247 /)
  real ( kind = 8 ), dimension ( 247 ) :: x = (/ &
    1.0D+00, &
    0.415774556783479083311533873128D+00, &
    0.229428036027904171982205036136D+01, &
    0.628994508293747919686641576551D+01, &
    0.193043676560362413838247885004D+00, &
    0.102666489533919195034519944317D+01, &
    0.256787674495074620690778622666D+01, &
    0.490035308452648456810171437810D+01, &
    0.818215344456286079108182755123D+01, &
    0.127341802917978137580126424582D+02, &
    0.193957278622625403117125820576D+02, &
    0.933078120172818047629030383672D-01, &
    0.492691740301883908960101791412D+00, &
    0.121559541207094946372992716488D+01, &
    0.226994952620374320247421741375D+01, &
    0.366762272175143727724905959436D+01, &
    0.542533662741355316534358132596D+01, &
    0.756591622661306786049739555812D+01, &
    0.101202285680191127347927394568D+02, &
    0.131302824821757235640991204176D+02, &
    0.166544077083299578225202408430D+02, &
    0.207764788994487667729157175676D+02, &
    0.256238942267287801445868285977D+02, &
    0.314075191697539385152432196202D+02, &
    0.385306833064860094162515167595D+02, &
    0.480260855726857943465734308508D+02, &
    0.45901947621108290743496080275224D-01, &
    0.24198016382477204890408974151714D+00, &
    0.59525389422235073707330165005414D+00, &
    1.1066894995329987162111308789792D+00, &
    1.7775956928747727211593727482675D+00, &
    2.6097034152566806503893375925315D+00, &
    3.6051968023400442698805817554243D+00, &
    4.7667470844717611313629127271123D+00, &
    6.0975545671817409269925429328463D+00, &
    7.6014009492331374229360106942867D+00, &
    9.2827143134708894182536695297710D+00, &
    11.146649755619291358993815629587D+00, &
    13.199189576244998522464925028637D+00, &
    15.447268315549310075809325891801D+00, &
    17.898929826644757646725793817752D+00, &
    20.563526336715822170743048968779D+00, &
    23.451973482011858591050255575933D+00, &
    26.577081352118260459975876986478D+00, &
    29.953990872346445506951917840024D+00, &
    33.600759532902202735410313885784D+00, &
    37.539164407330440882887902558001D+00, &
    41.795830870182219981347945853330D+00, &
    46.403866806411123136029227604386D+00, &
    51.405314476797755161861461088395D+00, &
    56.854992868715843620511922055660D+00, &
    62.826855908786321453677523304806D+00, &
    69.425277191080345623322251656443D+00, &
    76.807047763862732837609972285484D+00, &
    85.230358607545669169387065607043D+00, &
    95.188939891525629981308606853957D+00, &
    107.95224382757871475002440117666D+00, &
    0.22768893732576153785994330248562D-01, &
    0.11998325242727824715771416426383D+00, &
    0.29494185444770149577427738517405D+00, &
    0.54779087896237725363865073775856D+00, &
    0.87869061179931901673895567052285D+00, &
    1.2878464335919706302309207788611D+00, &
    1.7755123815388553763979463268728D+00, &
    2.3419925567085989256055628337716D+00, &
    2.9876423223246473939976731053629D+00, &
    3.7128695992018000346299637413422D+00, &
    4.5181363349503584391105568561550D+00, &
    5.4039601781825946286902599782736D+00, &
    6.3709163787865330220392250891777D+00, &
    7.4196399339311711154888493199004D+00, &
    8.5508280008403328312589048722235D+00, &
    9.7652425999245366807004592977996D+00, &
    11.063713635140661736220550410604D+00, &
    12.447142262356492749798687569289D+00, &
    13.916504641057818562912967008183D+00, &
    15.472856110036296424777143607779D+00, &
    17.117335833863588753116900303886D+00, &
    18.851171974154856850873483787506D+00, &
    20.675687448056515660377265667433D+00, &
    22.592306346311528381292277759986D+00, &
    24.602561094972638883700642760037D+00, &
    26.708100458737343969779087998829D+00, &
    28.910698500451382640177718103234D+00, &
    31.212264631175912885477773820802D+00, &
    33.614854909101154836598842888345D+00, &
    36.120684774484823056306328740825D+00, &
    38.732143442933582145626041607663D+00, &
    41.451810222318741191114726181363D+00, &
    44.282473071479233839358857134636D+00, &
    47.227149784295686898935095231536D+00, &
    50.289112264240695761749021839419D+00, &
    53.471914456788652808348280619542D+00, &
    56.779424636342062213099781057119D+00, &
    60.215862909019862886417550114424D+00, &
    63.785845004235974631701139601836D+00, &
    67.494433702293885830374325695045D+00, &
    71.347199604295266286654803376075D+00, &
    75.350293425653234254290504744279D+00, &
    79.510532629986309149555391354778D+00, &
    83.835506080872257843339817658508D+00, &
    88.333701570354369086112766326498D+00, &
    93.014662728558547405303399037100D+00, &
    97.889184147578140043386727677112D+00, &
    102.96955690741381650783952746778D+00, &
    108.26988161961595392226350967206D+00, &
    113.80647350287462738934485955901D+00, &
    119.59839538830458666962452963285D+00, &
    125.66817255856119431291196303280D+00, &
    132.04277272091165746585590583045D+00, &
    138.75498418103789078167590567526D+00, &
    145.84541318313540358283994248439D+00, &
    153.36548459497863623710815962660D+00, &
    161.38215194813761243562172669592D+00, &
    169.98570600665839438795175301156D+00, &
    179.30366247401580910251827858515D+00, &
    189.52789596532475473668721332981D+00, &
    200.97521159924656741628671841018D+00, &
    214.25368536638788642698056296400D+00, &
    230.93465747089703971246562985079D+00, &
    0.11339635298518611691893169631306D-01, &
    0.59749753435726620281348237057387D-01, &
    0.14685098690746167612388223687431D+00, &
    0.27267590735859553131378008278900D+00, &
    0.43724600644192665554577035869932D+00, &
    0.64058688222566929533576416399983D+00, &
    0.88272968639058364481487653650042D+00, &
    1.1637114160166537661560584700951D+00, &
    1.4835750152834613891313584861012D+00, &
    1.8423694351613565380686320809853D+00, &
    2.2401496839579024244513315656522D+00, &
    2.6769768780141303692167869961238D+00, &
    3.1529182957082825565771508308846D+00, &
    3.6680474360304752540226339926515D+00, &
    4.2224440823301888455977876667425D+00, &
    4.8161943715870502475665535087286D+00, &
    5.4493908694559416755862178908416D+00, &
    6.1221326512997254193944584763155D+00, &
    6.8345253894122668112237994973336D+00, &
    7.5866814466367472174205986836847D+00, &
    8.3787199765932725254842120659452D+00, &
    9.2107670307426558777922506102445D+00, &
    10.082955672528643809166439353647D+00, &
    10.995426098858125429803147358780D+00, &
    11.948325769197725997610605127857D+00, &
    12.941809542585531053723381098192D+00, &
    13.976039822878506520014405668679D+00, &
    15.051186712579523631574796365435D+00, &
    16.167428175612852922977395051768D+00, &
    17.324950209443673446561163712616D+00, &
    18.523947026965688560811711309349D+00, &
    19.764621248611504104071669386884D+00, &
    21.047184105173183606877044020054D+00, &
    22.371855651855542817648123918101D+00, &
    23.738864994122497183652313788712D+00, &
    25.148450525937368234077278385644D+00, &
    26.600860181041749607253384279755D+00, &
    28.096351697964619201753961292129D+00, &
    29.635192899504178910610227138642D+00, &
    31.217661987479759144214467152615D+00, &
    32.844047853610430460522951341338D+00, &
    34.514650407441149149105635947422D+00, &
    36.229780922306804019615388508885D+00, &
    37.989762400399956435968780140278D+00, &
    39.794929958089961778396437141707D+00, &
    41.645631232730180705153990897484D+00, &
    43.542226812286859549950892993822D+00, &
    45.485090689228791137996151336673D+00, &
    47.474610740231964719468766599146D+00, &
    49.511189233379087716728884584381D+00, &
    51.595243364671244443182771266934D+00, &
    53.727205825819316758288140069145D+00, &
    55.907525405447553305830605991732D+00, &
    58.136667626022439197077526025660D+00, &
    60.415115419018590295707192053805D+00, &
    62.743369841051809700207126742685D+00, &
    65.121950833949996311956025417139D+00, &
    67.551398031997886314411872443149D+00, &
    70.032271619884584511229871192030D+00, &
    72.565153245206849090888669416801D+00, &
    75.150646989739935299354362325096D+00, &
    77.789380404085816000647405462136D+00, &
    80.482005610750729205803962926758D+00, &
    83.229200481195914886796120019048D+00, &
    86.031669892953582966798238732643D+00, &
    88.890147073512051099652518544282D+00, &
    91.805395038358177994971250170499D+00, &
    94.778208131331583205387031034825D+00, &
    97.809413676305116411054110115424D+00, &
    100.89987375017285940371939762172D+00, &
    104.05048708821598934704076845022D+00, &
    107.26219113414600428423116401414D+00, &
    110.53596424851500530602771351277D+00, &
    113.87282809075839485348376187652D+00, &
    117.27385019192517774095477886379D+00, &
    120.74014673718880106173978002719D+00, &
    124.27288557955698354259506446928D+00, &
    127.87328950885942645093841745425D+00, &
    131.54263980314366921809377742137D+00, &
    135.28228009311836970132738106369D+00, &
    139.09362057432970013964422086977D+00, &
    142.97814260643601776808227753574D+00, &
    146.93740374437366549441080969072D+00, &
    150.97304325252187127492511437460D+00, &
    155.08678816034612572229641420609D+00, &
    159.28045992663288235401956989889D+00, &
    163.55598178957571104015967182053D+00, &
    167.91538689194360134245547184721D+00, &
    172.36082728473812536838156191681D+00, &
    176.89458392960192176311674993508D+00, &
    181.51907784036813069227528834025D+00, &
    186.23688252828112373861202530357D+00, &
    191.05073794450929196790836610789D+00, &
    195.96356614879879837839002542988D+00, &
    200.97848897600025153696475526130D+00, &
    206.09884802468871112127283042753D+00, &
    211.32822735671655260572377256981D+00, &
    216.67047937658230323477089465777D+00, &
    222.12975445929687246267304963754D+00, &
    227.71053502072232419089132431317D+00, &
    233.41767488282602453367775322563D+00, &
    239.25644498830308620018749667089D+00, &
    245.23258677871567172531254018984D+00, &
    251.35237488718128030005500991754D+00, &
    257.62269123792061413076191882313D+00, &
    264.05111322908240551754377241831D+00, &
    270.64601945722796749299111718606D+00, &
    277.41671750163651071798388218104D+00, &
    284.37359974220870326674402873120D+00, &
    291.52833521346495719581282021650D+00, &
    298.89410837028248600878895615414D+00, &
    306.48591978262611320418112423947D+00, &
    314.32096986471177487400007507615D+00, &
    322.41915589128679683349440361344D+00, &
    330.80372663802405651933847334878D+00, &
    339.50216127832433747735367595958D+00, &
    348.54737559472697355480761787441D+00, &
    357.97942028029845454049007443090D+00, &
    367.84794520076004578858341422871D+00, &
    378.21590623135532818332979188889D+00, &
    389.16539141251004101579475325153D+00, &
    400.80729331451702589996361286427D+00, &
    413.29853681779384418008260081859D+00, &
    426.87579153663675538288509017051D+00, &
    441.93085485310841412460309271842D+00, &
    459.21804639888429981971267313224D+00, &
    480.69378263388373859704269229304D+00  &
         /)

  if ( any ( grid_base(1:dim_num) < 1 ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LG_ABSCISSA - Fatal error!'
    write ( *, '(a)' ) '  Some base values are less than 1.'
    stop
  end if

  if ( any ( 127 < grid_base(1:dim_num) ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LG_ABSCISSA - Fatal error!'
    write ( *, '(a)' ) '  Some base values are greater than 127.'
    stop
  end if

  do point = 1, point_num
    do dim = 1, dim_num

      level = i4_log_2 ( grid_base(dim) + 1 ) - 1

      pointer = skip(level) + grid_index(dim,point)

      if ( pointer < 1 .or. 247 < pointer ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'LG_ABSCISSA - Fatal error!'
        write ( *, '(a)' ) '  POINTER out of bounds.'
        write ( *, '(a,i8)' ) '  POINTER    = ', pointer
        write ( *, '(a,i8)' ) '  POINT      = ', point
        write ( *, '(a,i8)' ) '  DIM        = ', dim
        write ( *, '(a,i8)' ) '  GRID_BASE  = ', grid_base(dim)
        write ( *, '(a,i8)' ) '  LEVEL      = ', level
        write ( *, '(a,i8)' ) '  GRID_INDEX = ', grid_index(dim,point)
        stop
      end if

      grid_point(dim,point) = x(pointer)

    end do
  end do

  return
end
subroutine lg_weights ( order, weight )

!*****************************************************************************80
!
!! LG_WEIGHTS returns weights for certain Gauss Laguerre quadrature rules.
!
!  Discussion:
!
!    The allowed orders are 1, 3, 7, 15, 31, 63, and 127.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ORDER, the order of the rule.
!    ORDER must be 1, 3, 7, 15, 31, 63 or 127.
!
!    Output, real ( kind = 8 ) WEIGHT(ORDER), the weights.
!    The weights are positive, symmetric and should sum to 1.
!
  implicit none

  integer ( kind = 4 ) order

  real ( kind = 8 ) weight(order)

  if ( order == 1 ) then

    weight(1) = 1.0D+00

  else if ( order == 3 ) then

    weight(1) =  0.711093009929173015449590191143D+00
    weight(2) =  0.278517733569240848801444888457D+00
    weight(3) =  0.103892565015861357489649204007D-01

  else if ( order == 7 ) then

    weight(1) =  0.409318951701273902130432880018D+00
    weight(2) =  0.421831277861719779929281005417D+00
    weight(3) =  0.147126348657505278395374184637D+00
    weight(4) =  0.206335144687169398657056149642D-01
    weight(5) =  0.107401014328074552213195962843D-02
    weight(6) =  0.158654643485642012687326223234D-04
    weight(7) =  0.317031547899558056227132215385D-07

  else if ( order == 15 ) then

    weight(1) =  0.218234885940086889856413236448D+00
    weight(2) =  0.342210177922883329638948956807D+00
    weight(3) =  0.263027577941680097414812275022D+00
    weight(4) =  0.126425818105930535843030549378D+00
    weight(5) =  0.402068649210009148415854789871D-01
    weight(6) =  0.856387780361183836391575987649D-02
    weight(7) =  0.121243614721425207621920522467D-02
    weight(8) =  0.111674392344251941992578595518D-03
    weight(9) =  0.645992676202290092465319025312D-05
    weight(10) = 0.222631690709627263033182809179D-06
    weight(11) = 0.422743038497936500735127949331D-08
    weight(12) = 0.392189726704108929038460981949D-10
    weight(13) = 0.145651526407312640633273963455D-12
    weight(14) = 0.148302705111330133546164737187D-15
    weight(15) = 0.160059490621113323104997812370D-19

  else if ( order == 31 ) then

    weight(  1) =   0.11252789550372583820847728082801D+00
    weight(  2) =   0.21552760818089123795222505285045D+00
    weight(  3) =   0.23830825164569654731905788089234D+00
    weight(  4) =   0.19538830929790229249915303390711D+00
    weight(  5) =   0.12698283289306190143635272904602D+00
    weight(  6) =   0.67186168923899300670929441993508D-01
    weight(  7) =   0.29303224993879487404888669311974D-01
    weight(  8) =   0.10597569915295736089529380314433D-01
    weight(  9) =   0.31851272582386980320974842433019D-02
    weight( 10) =   0.79549548307940382922092149012477D-03
    weight( 11) =   0.16480052126636687317862967116412D-03
    weight( 12) =   0.28229237864310816393860971468993D-04
    weight( 13) =   0.39802902551008580387116174900106D-05
    weight( 14) =   0.45931839841801061673729694510289D-06
    weight( 15) =   0.43075545187731100930131457465897D-07
    weight( 16) =   0.32551249938271570855175749257884D-08
    weight( 17) =   0.19620246675410594996247151593142D-09
    weight( 18) =   0.93190499086617587129534716431331D-11
    weight( 19) =   0.34377541819411620520312597898311D-12
    weight( 20) =   0.96795247130446716997405035776206D-14
    weight( 21) =   0.20368066110115247398010624219291D-15
    weight( 22) =   0.31212687280713526831765358632585D-17
    weight( 23) =   0.33729581704161052453395678308350D-19
    weight( 24) =   0.24672796386616696011038363242541D-21
    weight( 25) =   0.11582201904525643634834564576593D-23
    weight( 26) =   0.32472922591425422434798022809020D-26
    weight( 27) =   0.49143017308057432740820076259666D-29
    weight( 28) =   0.34500071104808394132223135953806D-32
    weight( 29) =   0.87663710117162041472932760732881D-36
    weight( 30) =   0.50363643921161490411297172316582D-40
    weight( 31) =   0.19909984582531456482439549080330D-45

  else if ( order == 63 ) then

    weight(  1) =   0.57118633213868979811587283390476D-01
    weight(  2) =   0.12067476090640395283319932036351D+00
    weight(  3) =   0.15925001096581873723870561096472D+00
    weight(  4) =   0.16875178327560799234596192963585D+00
    weight(  5) =   0.15366641977668956696193711310131D+00
    weight(  6) =   0.12368770614716481641086652261948D+00
    weight(  7) =   0.89275098854848671545279150057422D-01
    weight(  8) =   0.58258485446105944957571825725160D-01
    weight(  9) =   0.34546657545992580874717085812508D-01
    weight( 10) =   0.18675685985714656798286552591203D-01
    weight( 11) =   0.92233449044093536528490075241649D-02
    weight( 12) =   0.41671250684839592762582663470209D-02
    weight( 13) =   0.17238120299900582715386728541955D-02
    weight( 14) =   0.65320845029716311169340559359043D-03
    weight( 15) =   0.22677644670909586952405173207471D-03
    weight( 16) =   0.72127674154810668410750270234861D-04
    weight( 17) =   0.21011261180466484598811536851241D-04
    weight( 18) =   0.56035500893357212749181536071292D-05
    weight( 19) =   0.13673642785604888017836641282292D-05
    weight( 20) =   0.30507263930195817240736097189550D-06
    weight( 21) =   0.62180061839309763559981775409241D-07
    weight( 22) =   0.11566529551931711260022448996296D-07
    weight( 23) =   0.19614588267565478081534781863335D-08
    weight( 24) =   0.30286171195709411244334756404054D-09
    weight( 25) =   0.42521344539400686769012963452599D-10
    weight( 26) =   0.54202220578073819334698791381873D-11
    weight( 27) =   0.62627306838597672554166850420603D-12
    weight( 28) =   0.65474443156573322992307089591924D-13
    weight( 29) =   0.61815575808729181846302500000047D-14
    weight( 30) =   0.52592721363507381404263991342633D-15
    weight( 31) =   0.40230920092646484015391506025408D-16
    weight( 32) =   0.27600740511819536505013824207729D-17
    weight( 33) =   0.16936946756968296053322009855265D-18
    weight( 34) =   0.92689146872177087314963772462726D-20
    weight( 35) =   0.45093739060365632939780140603959D-21
    weight( 36) =   0.19435162876132376573629962695374D-22
    weight( 37) =   0.73926270895169207037999639194513D-24
    weight( 38) =   0.24714364154434632615980126000066D-25
    weight( 39) =   0.72288649446741597655145390616476D-27
    weight( 40) =   0.18407617292614039362985209905608D-28
    weight( 41) =   0.40583498566841960105759537058880D-30
    weight( 42) =   0.77000496416438368114463925286343D-32
    weight( 43) =   0.12488505764999334328843314866038D-33
    weight( 44) =   0.17185000226767010697663950619912D-35
    weight( 45) =   0.19896372636672396938013975755522D-37
    weight( 46) =   0.19199671378804058267713164416870D-39
    weight( 47) =   0.15278588285522166920459714708240D-41
    weight( 48) =   0.99054752688842142955854138884590D-44
    weight( 49) =   0.51597523673029211884228858692990D-46
    weight( 50) =   0.21249846664084111245693912887783D-48
    weight( 51) =   0.67903852766852910591172042494884D-51
    weight( 52) =   0.16466654148296177467908300517887D-53
    weight( 53) =   0.29509065402691055027053659375033D-56
    weight( 54) =   0.37838420647571051984882241014675D-59
    weight( 55) =   0.33358130068542431878174667995217D-62
    weight( 56) =   0.19223461022273880981363303073329D-65
    weight( 57) =   0.67812696961083016872779388922288D-69
    weight( 58) =   0.13404752802440604607620468935693D-72
    weight( 59) =   0.13109745101805029757648048223928D-76
    weight( 60) =   0.52624863881401787388694579143866D-81
    weight( 61) =   0.63780013856587414257760666006511D-86
    weight( 62) =   0.12997078942372924566347473916943D-91
    weight( 63) =   0.10008511496968754063443740168421D-98

  else if ( order == 127 ) then

    weight(  1) =   0.28773246692000124355770010301506D-01
    weight(  2) =   0.63817468175134649363480949265236D-01
    weight(  3) =   0.91919669721570571389864194652717D-01
    weight(  4) =   0.11054167914413766381245463002967D+00
    weight(  5) =   0.11879771633375850188328329422643D+00
    weight(  6) =   0.11737818530052695148804451630074D+00
    weight(  7) =   0.10819305984180551488335145581193D+00
    weight(  8) =   0.93827075290489628080377261401107D-01
    weight(  9) =   0.76966450960588843995822485928431D-01
    weight( 10) =   0.59934903912939714332570730063476D-01
    weight( 11) =   0.44417742073889001371708316272923D-01
    weight( 12) =   0.31385080966252320983009372215062D-01
    weight( 13) =   0.21172316041924506411370709025015D-01
    weight( 14) =   0.13650145364230541652171185564626D-01
    weight( 15) =   0.84172852710599172279366657385445D-02
    weight( 16) =   0.49674990059882760515912858620175D-02
    weight( 17) =   0.28069903895001884631961957446400D-02
    weight( 18) =   0.15192951003941952460445341057817D-02
    weight( 19) =   0.78789028751796084086217287140548D-03
    weight( 20) =   0.39156751064868450584507324648999D-03
    weight( 21) =   0.18652434268825860550093566260060D-03
    weight( 22) =   0.85173160415576621908809828160247D-04
    weight( 23) =   0.37285639197853037712145321577724D-04
    weight( 24) =   0.15648416791712993947447805296768D-04
    weight( 25) =   0.62964340695224829035692735524979D-05
    weight( 26) =   0.24288929711328724574541379938222D-05
    weight( 27) =   0.89824607890051007201922871545035D-06
    weight( 28) =   0.31844174740760353710742966328091D-06
    weight( 29) =   0.10821272905566839211861807542741D-06
    weight( 30) =   0.35245076750635536015902779085340D-07
    weight( 31) =   0.11001224365719347407063839761738D-07
    weight( 32) =   0.32904079616717932125329343003261D-08
    weight( 33) =   0.94289145237889976419772700772988D-09
    weight( 34) =   0.25882578904668318184050195309296D-09
    weight( 35) =   0.68047437103370762630942259017560D-10
    weight( 36) =   0.17131398805120837835399564475632D-10
    weight( 37) =   0.41291744524052865469443922304935D-11
    weight( 38) =   0.95264189718807273220707664873469D-12
    weight( 39) =   0.21032604432442425932962942047474D-12
    weight( 40) =   0.44427151938729352860940434285789D-13
    weight( 41) =   0.89760500362833703323319846405449D-14
    weight( 42) =   0.17341511407769287074627948346848D-14
    weight( 43) =   0.32028099548988356631494379835210D-15
    weight( 44) =   0.56531388950793682022660742095189D-16
    weight( 45) =   0.95329672799026591234588044025896D-17
    weight( 46) =   0.15353453477310142565288509437552D-17
    weight( 47) =   0.23608962179467365686057842132176D-18
    weight( 48) =   0.34648742794456611332193876653230D-19
    weight( 49) =   0.48515241897086461320126957663545D-20
    weight( 50) =   0.64786228633519813428137373790678D-21
    weight( 51) =   0.82476020965403242936448553126316D-22
    weight( 52) =   0.10005361880214719793491658282977D-22
    weight( 53) =   0.11561395116207304954233181263632D-23
    weight( 54) =   0.12719342731167922655612134264961D-24
    weight( 55) =   0.13316584714165372967340004160814D-25
    weight( 56) =   0.13261218454678944033646108509198D-26
    weight( 57) =   0.12554995447643949807286074138324D-27
    weight( 58) =   0.11294412178579462703240913107219D-28
    weight( 59) =   0.96491020279562119228500608131696D-30
    weight( 60) =   0.78241846768302099396733076955632D-31
    weight( 61) =   0.60181503542219626658249939076636D-32
    weight( 62) =   0.43882482704961741551510518054138D-33
    weight( 63) =   0.30314137647517256304035802501863D-34
    weight( 64) =   0.19826016543944539545224676057020D-35
    weight( 65) =   0.12267623373665926559013654872402D-36
    weight( 66) =   0.71763931692508888943812834967620D-38
    weight( 67) =   0.39659378833836963584113716149270D-39
    weight( 68) =   0.20688970553868040099581951696677D-40
    weight( 69) =   0.10179587017979517245268418427523D-41
    weight( 70) =   0.47200827745986374625714293679649D-43
    weight( 71) =   0.20606828985553374825744353490744D-44
    weight( 72) =   0.84627575907305987245899032156188D-46
    weight( 73) =   0.32661123687088798658026998931647D-47
    weight( 74) =   0.11833939207883162380564134612682D-48
    weight( 75) =   0.40211209123895013807243250164050D-50
    weight( 76) =   0.12799824394111125389430292847476D-51
    weight( 77) =   0.38123877747548846504399051365162D-53
    weight( 78) =   0.10612057542701156767898551949650D-54
    weight( 79) =   0.27571446947200403594113572720812D-56
    weight( 80) =   0.66772544240928492881306904862856D-58
    weight( 81) =   0.15052438383868234954068178600268D-59
    weight( 82) =   0.31538986800113758526689068500772D-61
    weight( 83) =   0.61326614299483180785237418887960D-63
    weight( 84) =   0.11048510030324810567549119229368D-64
    weight( 85) =   0.18410563538091348076979665543900D-66
    weight( 86) =   0.28323926570052832195543883237652D-68
    weight( 87) =   0.40154409843763655508670978777418D-70
    weight( 88) =   0.52351530215683708779772201956106D-72
    weight( 89) =   0.62634476665005100555787696642851D-74
    weight( 90) =   0.68612210535666530365348093803922D-76
    weight( 91) =   0.68651298840956019297134099761855D-78
    weight( 92) =   0.62581388433728084867318704240915D-80
    weight( 93) =   0.51833271237514904046803469968027D-82
    weight( 94) =   0.38893621571918443533108973497673D-84
    weight( 95) =   0.26357711379476932781525533730623D-86
    weight( 96) =   0.16078851293917979699005509638883D-88
    weight( 97) =   0.87978042070968939637972577886624D-91
    weight( 98) =   0.43013405077495109903408697802188D-93
    weight( 99) =   0.18713435881342838527144321803729D-95
    weight(100) =   0.72125744708060471675805761366523D-98
    weight(101) =   0.24508746062177874383231742333023D-100
    weight(102) =   0.73042094619470875777647865078327D-103
    weight(103) =   0.18983290818383463537886818579820D-105
    weight(104) =   0.42757400244246684123093264825902D-108
    weight(105) =   0.82894681420515755691423485228897D-111
    weight(106) =   0.13729432219324400013067050156048D-113
    weight(107) =   0.19265464126404973222043166489406D-116
    weight(108) =   0.22693344503301354826140809941334D-119
    weight(109) =   0.22209290603717355061909071271535D-122
    weight(110) =   0.17851087685544512662856555121755D-125
    weight(111) =   0.11630931990387164467431190485525D-128
    weight(112) =   0.60524443584652392290952805077893D-132
    weight(113) =   0.24729569115063528647628375096400D-135
    weight(114) =   0.77789065006489410364997205809045D-139
    weight(115) =   0.18409738662712607039570678274636D-142
    weight(116) =   0.31900921131079114970179071968597D-146
    weight(117) =   0.39179487139174199737617666077555D-150
    weight(118) =   0.32782158394188697053774429820559D-154
    weight(119) =   0.17793590713138888062819640128739D-158
    weight(120) =   0.58882353408932623157467835381214D-163
    weight(121) =   0.10957236509071169877747203273886D-167
    weight(122) =   0.10281621114867000898285076975760D-172
    weight(123) =   0.41704725557697758145816510853967D-178
    weight(124) =   0.58002877720316101774638319601971D-184
    weight(125) =   0.18873507745825517106171619101120D-190
    weight(126) =   0.69106601826730911682786705950895D-198
    weight(127) =   0.43506813201105855628383313334402D-207

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LG_WEIGHTS - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of ORDER = ', order
    write ( *, '(a)' ) '  Legal values are 1, 3, 7, 15, 31, 63 and 127.'
    stop

  end if

  return
end
subroutine monomial_integral_hermite ( dim_num, expon, value )

!*****************************************************************************80
!
!! MONOMIAL_INTEGRAL_HERMITE integrates a Hermite monomial.
!
!  Discussion:
!
!    H(d,n) = Integral ( -Infinity < x < Infinity )
!      x1^n1 * x2^n2...*xd^nd * exp(-x1^2-x2^2...-xd^2 ) dx
!
!    H(d,n) is 0 if any n(i) odd.
!
!    H(d,n) = product ( 1 <= i <= d )
!      ( (n(i)-1)!! * sqrt(pi) / 2^(n(i)/2) for all n(i) even.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the integral.
!
!    Input, integer ( kind = 4 ) EXPON(DIM_NUM), the order of the integral.
!    0 <= EXPON(1:DIM_NUM).
!
!    Output, real ( kind = 8 ) VALUE, the value of the integral.
!
  implicit none

  integer ( kind = 4 ) dim_num

  integer ( kind = 4 ) dim
  integer ( kind = 4 ), parameter :: i4_1 = 1
  integer ( kind = 4 ), parameter :: i4_2 = 2
  integer ( kind = 4 ) expon(dim_num)
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r8_factorial2
  real ( kind = 8 ) r8_huge
  real ( kind = 8 ) value

  if ( any ( expon(1:dim_num) < 0 ) ) then

    value = - r8_huge ( )

  else if ( any ( mod ( expon(1:dim_num), i4_2 ) == 1 ) ) then

    value = 0.0D+00

  else

    value = 1.0D+00
    do dim = 1, dim_num
      value = value * r8_factorial2 ( expon(dim) - i4_1 ) * sqrt ( pi ) &
        / 2.0D+00**( expon(dim) / 2 )
    end do

  end if

  return
end
subroutine monomial_integral_laguerre ( dim_num, expon, value )

!*****************************************************************************80
!
!! MONOMIAL_INTEGRAL_LAGUERRE integrates a Laguerre monomial.
!
!  Discussion:
!
!    L(1,n) = Integral ( 0 <= x < Infinity ) x^n exp ( -x ) dx
!           = n!
!
!    L(d,n) = Integral ( 0 <= x(i) < Infinity )
!             x1^n1 * x2^n2...*xd^nd * exp(-x1-x2...-xd ) dx
!           = Product ( n(i)! ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the integral.
!
!    Input, integer ( kind = 4 ) EXPON(DIM_NUM), the exponents.
!
!    Output, real ( kind = 8 ) VALUE, the value of the integral.
!
  implicit none

  integer ( kind = 4 ) dim_num

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) expon(dim_num)
  real ( kind = 8 ) r8_factorial
  real ( kind = 8 ) value

  value = 1.0D+00
  do dim = 1, dim_num
    value = value * r8_factorial ( expon(dim) )
  end do

  return
end
subroutine monomial_integral_legendre ( dim_num, expon, value )

!*****************************************************************************80
!
!! MONOMIAL_INTEGRAL_LEGENDRE integrates a Legendre monomial.
!
!  Discussion:
!
!    This routine returns the exact value of a multidimensional Legendre
!    type integral:
!
!      integral ( -1 <= x(1:n) <= +1 ) f(x) dx
!
!    where f(x) is a monomial of the form:
!
!      product ( 1 <= dim <= dim_num ) x(dim)^expon(dim)
!
!    and the exponents are nonnegative integers.  Note that the combination
!    0^0 is treated as 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) EXPON(DIM_NUM), the exponents.
!
!    Output, real ( kind = 8 ) VALUE, the integral of the monomial.
!
  implicit none

  integer ( kind = 4 ) dim_num

  integer ( kind = 4 ) expon(dim_num)
  real ( kind = 8 ) value

  if ( any ( mod ( expon(1:dim_num), 2 ) == 1 ) ) then
    value = 0.0D+00
  else
    value = 2.0D+00**dim_num &
      / real ( product ( expon(1:dim_num) + 1 ), kind = 8 )
  end if

  return
end
subroutine monomial_quadrature ( dim_num, expon, point_num, weight, &
  x, rule, quad_error )

!*****************************************************************************80
!
!! MONOMIAL_QUADRATURE applies a quadrature rule to a monomial.
!
!  Discussion:
!
!    This routine assumes that the integral being approximated is that of
!    a multidimensional monomial, integrated over the [-1,+1] hypercube,
!    with a Legendre weight (that is, w(x) = 1).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) EXPON(DIM_NUM), the exponents.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points in the rule.
!
!    Input, real ( kind = 8 ) WEIGHT(POINT_NUM), the quadrature weights.
!
!    Input, real ( kind = 8 ) X(DIM_NUM,POINT_NUM), the quadrature points.
!
!    Input, integer ( kind = 4 ) RULE, the index of the rule.
!    1, "CC", Clenshaw Curtis Closed Fully Nested rule.
!    2, "F1", Fejer 1 Open Fully Nested rule.
!    3, "F2", Fejer 2 Open Fully Nested rule.
!    4, "GP", Gauss Patterson Open Fully Nested rule.
!    5, "GL", Gauss Legendre Open Weakly Nested rule.
!    6, "GH", Gauss Hermite Open Weakly Nested rule.
!    7, "LG", Gauss Laguerre Open Non Nested rule.
!
!    Output, real ( kind = 8 ) QUAD_ERROR, the quadrature error.
!
  implicit none

  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) exact
  integer ( kind = 4 ) expon(dim_num)
  integer ( kind = 4 ) point_num
  real ( kind = 8 ) quad
  real ( kind = 8 ) quad_error
  integer ( kind = 4 ) rule
  real ( kind = 8 ) value(point_num)
  real ( kind = 8 ) weight(point_num)
  real ( kind = 8 ) x(dim_num,point_num)
!
!  Get the exact value of the integral of the monomial.
!
  if ( 1 <= rule .and. rule <= 5 ) then
    call monomial_integral_legendre ( dim_num, expon, exact )
  else if ( rule == 6 ) then
    call monomial_integral_hermite ( dim_num, expon, exact )
  else if ( rule == 7 ) then
    call monomial_integral_laguerre ( dim_num, expon, exact )
  end if
!
!  Evaluate the monomial at the quadrature points.
!
  call monomial_value ( dim_num, point_num, x, expon, value )
!
!  Compute the quadrature sum.
!
  quad = dot_product ( weight, value )
!
!  Absolute error if EXACT = 0, relative error otherwise:
!
  if ( exact == 0.0D+00 ) then
    quad_error = abs ( quad - exact )
  else
    quad_error = abs ( quad - exact ) / abs ( exact )
  end if

  return
end
subroutine monomial_value ( dim_num, point_num, x, expon, value )

!*****************************************************************************80
!
!! MONOMIAL_VALUE evaluates a monomial.
!
!  Discussion:
!
!    This routine evaluates a monomial of the form
!
!      product ( 1 <= dim <= dim_num ) x(dim)^expon(dim)
!
!    where the exponents are nonnegative integers.  Note that
!    the combination 0^0 is treated as 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points at which the
!    monomial is to be evaluated.
!
!    Input, real ( kind = 8 ) X(DIM_NUM,POINT_NUM), the point coordinates.
!
!    Input, integer ( kind = 4 ) EXPON(DIM_NUM), the exponents.
!
!    Output, real ( kind = 8 ) VALUE(POINT_NUM), the value of the monomial.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) expon(dim_num)
  real ( kind = 8 ) value(point_num)
  real ( kind = 8 ) x(dim_num,point_num)

  value(1:point_num) = 1.0D+00

  do dim = 1, dim_num
    if ( 0 /= expon(dim) ) then
      value(1:point_num) = value(1:point_num) * x(dim,1:point_num)**expon(dim)
    end if
  end do

  return
end
subroutine multigrid_index_cfn ( dim_num, order_1d, order_nd, indx )

!*****************************************************************************80
!
!! MULTIGRID_INDEX_CFN indexes a sparse grid based on CFN 1D rules.
!
!  Discussion:
!
!    The sparse grid is presumed to have been created from products
!    of CLOSED FULLY NESTED 1D quadrature rules.
!
!    CFN rules include Clenshaw Curtis rules.
!
!    For dimension DIM, the second index of INDX may vary from
!    0 to ORDER_1D(DIM)-1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension of the points.
!
!    Input, integer ( kind = 4 ) ORDER_1D(DIM_NUM), the order of the
!    rule in each dimension.
!
!    Input, integer ( kind = 4 ) ORDER_ND, the product of the entries
!    of ORDER_1D.
!
!    Output, integer ( kind = 4 ) INDX(DIM_NUM,ORDER_ND), the indices of the
!    points in the grid.  The second dimension of this array is equal to the
!    product of the entries of ORDER_1D.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) order_nd

  integer ( kind = 4 ) a(dim_num)
  logical more
  integer ( kind = 4 ) order_1d(dim_num)
  integer ( kind = 4 ) p
  integer ( kind = 4 ) indx(dim_num,order_nd)

  more = .false.
  p = 0

  do

    call vec_colex_next2 ( dim_num, order_1d, a, more )

    if ( .not. more ) then
      exit
    end if

    p = p + 1

    indx(1:dim_num,p) = a(1:dim_num)

  end do

  return
end
subroutine multigrid_index_ofn ( dim_num, order_1d, order_nd, indx )

!*****************************************************************************80
!
!! MULTIGRID_INDEX_OFN indexes a sparse grid based on OFN 1D rules.
!
!  Discussion:
!
!    The sparse grid is presumed to have been created from products
!    of OPEN FULLY NESTED 1D quadrature rules.
!
!    OFN rules include Fejer 1, Fejer 2, and Gauss Patterson rules.
!
!    For dimension DIM, the second index of INDX may vary from
!    1 to ORDER_1D(DIM).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension of the points.
!
!    Input, integer ( kind = 4 ) ORDER_1D(DIM_NUM), the order of the
!    rule in each dimension.
!
!    Input, integer ( kind = 4 ) ORDER_ND, the product of the entries of
!    ORDER_1D.
!
!    Output, integer ( kind = 4 ) INDX(DIM_NUM,ORDER_ND), the indices of the
!    points in the grid.  The second dimension of this array is equal to the
!    product of the entries of ORDER_1D.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) order_nd

  integer ( kind = 4 ) a(dim_num)
  logical more
  integer ( kind = 4 ) order_1d(dim_num)
  integer ( kind = 4 ) p
  integer ( kind = 4 ) indx(dim_num,order_nd)

  more = .false.
  p = 0

  do

    call vec_colex_next2 ( dim_num, order_1d, a, more )

    if ( .not. more ) then
      exit
    end if

    p = p + 1

    indx(1:dim_num,p) = a(1:dim_num) + 1

  end do

  return
end
subroutine multigrid_index_onn ( dim_num, order_1d, order_nd, indx )

!*****************************************************************************80
!
!! MULTIGRID_INDEX_ONN indexes a sparse grid based on ONN 1D rules.
!
!  Discussion:
!
!    The sparse grid is presumed to have been created from products
!    of OPEN NON-NESTED 1D quadrature rules.
!
!    ONN rules include Gauss Laguerre.
!
!    For dimension DIM, the number of points is ORDER_1D(DIM).
!
!    We index the points as
!      1, 2, 3, ..., ORDER_1D(DIM).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension of the points.
!
!    Input, integer ( kind = 4 ) ORDER_1D(DIM_NUM), the order of the
!    rule in each dimension.
!
!    Input, integer ( kind = 4 ) ORDER_ND, the product of the entries
!    of ORDER_1D.
!
!    Output, integer ( kind = 4 ) INDX(DIM_NUM,ORDER_ND), the indices of
!    the points in the grid.  The second dimension of this array is equal
!    to the product of the entries of ORDER_1D.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) order_nd

  integer ( kind = 4 ) a(dim_num)
  logical more
  integer ( kind = 4 ) order_1d(dim_num)
  integer ( kind = 4 ) p
  integer ( kind = 4 ) indx(dim_num,order_nd)

  more = .false.
  p = 0

  do

    call vec_colex_next2 ( dim_num, order_1d, a, more )

    if ( .not. more ) then
      exit
    end if

    p = p + 1
!
!  The values of A(DIM) are between 0 and ORDER_1D(DIM)-1.
!
    indx(1:dim_num,p) = a(1:dim_num) + 1

  end do

  return
end
subroutine multigrid_index_own ( dim_num, order_1d, order_nd, indx )

!*****************************************************************************80
!
!! MULTIGRID_INDEX_OWN indexes a sparse grid based on OWN 1D rules.
!
!  Discussion:
!
!    The sparse grid is presumed to have been created from products
!    of OPEN WEAKLY NESTED 1D quadrature rules.
!
!    OWN rules include Gauss Legendre or Gauss Hermite.
!
!    For dimension DIM, the number of points is ORDER_1D(DIM).
!
!    We assume that ORDER_1D(DIM) is an odd number,
!      ORDER_1D(DIM) = N = 2 * M + 1
!    so that the points have coordinates
!      -M/M, -(M-1)/M, ..., -1/M, 0/M, 1/M, 2/M, 3/M, ..., (M-1)/M, M/M.
!    and we index them as
!      -M,   -(M-1),        -1,   0,   1,   2,   3,   ...,  M-1,    M.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension of the points.
!
!    Input, integer ( kind = 4 ) ORDER_1D(DIM_NUM), the order of the
!    rule in each dimension.
!
!    Input, integer ( kind = 4 ) ORDER_ND, the product of the entries
!    of ORDER_1D.
!
!    Output, integer ( kind = 4 ) INDX(DIM_NUM,ORDER_ND), the indices of
!    the points in the grid.  The second dimension of this array is equal
!    to the product of the entries of ORDER_1D.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) order_nd

  integer ( kind = 4 ) a(dim_num)
  logical              more
  integer ( kind = 4 ) order_1d(dim_num)
  integer ( kind = 4 ) p
  integer ( kind = 4 ) indx(dim_num,order_nd)

  more = .false.
  p = 0

  do

    call vec_colex_next2 ( dim_num, order_1d, a, more )

    if ( .not. more ) then
      exit
    end if

    p = p + 1
!
!  The values of A(DIM) are between 0 and ORDER_1D(DIM)-1 = N - 1 = 2 * M.
!  Subtracting M sets the range to -M to +M, as we wish.
!
    indx(1:dim_num,p) = a(1:dim_num) - ( order_1d(1:dim_num) - 1 ) / 2

  end do

  return
end
subroutine multigrid_scale_closed ( dim_num, order_nd, level_max, level_1d, &
  grid_index )

!*****************************************************************************80
!
!! MULTIGRID_SCALE_CLOSED renumbers a grid as a subgrid on a higher level.
!
!  Discussion:
!
!    This routine takes a grid associated with a given value of
!    LEVEL, and multiplies all the indices by a power of 2, so that
!    the indices reflect the position of the same points, but in
!    a grid of level LEVEL_MAX.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) ORDER_ND, the number of points in the grid.
!
!    Input, integer ( kind = 4 ) LEVEL_MAX, the maximum value of LEVEL.
!
!    Input, integer ( kind = 4 ) LEVEL_1D(DIM_NUM), the level in each dimension.
!
!    Input/output, integer ( kind = 4 ) GRID_INDEX(DIM_NUM,POINT_NUM), the index
!    values for each grid point.  On input, these indices are based in
!    the level for which the grid was generated; on output, the
!    indices are appropriate for the grid as a subgrid of a grid
!    of level LEVEL_MAX.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) order_nd

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) factor
  integer ( kind = 4 ) grid_index(dim_num,order_nd)
  integer ( kind = 4 ) level_1d(dim_num)
  integer ( kind = 4 ) level_max
  integer ( kind = 4 ) order_max

  do dim = 1, dim_num

    if ( level_1d(dim) == 0 ) then

      if ( 0 == level_max ) then
        order_max = 1
      else
        order_max = 2**level_max + 1
      end if

      grid_index(dim,1:order_nd) = ( order_max - 1 ) / 2

    else

      factor = 2**( level_max - level_1d(dim) )

      grid_index(dim,1:order_nd) = grid_index(dim,1:order_nd) * factor

    end if

  end do

  return
end
subroutine multigrid_scale_open ( dim_num, order_nd, level_max, level_1d, &
  grid_index )

!*****************************************************************************80
!
!! MULTIGRID_SCALE_OPEN renumbers a grid as a subgrid on a higher level.
!
!  Discussion:
!
!    This routine takes a grid associated with a given value of
!    LEVEL, and multiplies all the indices by a power of 2, so that
!    the indices reflect the position of the same points, but in
!    a grid of level LEVEL_MAX.
!
!    For an open grid, going from one level to the next, a set of indices
!    will be rescaled by 2*INDEX-1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) ORDER_ND, the number of points in the grid.
!
!    Input, integer ( kind = 4 ) LEVEL_MAX, the maximum value of LEVEL.
!
!    Input, integer ( kind = 4 ) LEVEL_1D(DIM_NUM), the level in each dimension.
!
!    Input/output, integer ( kind = 4 ) GRID_INDEX(DIM_NUM,POINT_NUM), the index
!    values for each grid point.  On input, these indices are based in
!    the level for which the grid was generated; on output, the
!    indices are appropriate for the grid as a subgrid of a grid
!    of level LEVEL_MAX.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) order_nd

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) factor
  integer ( kind = 4 ) grid_index(dim_num,order_nd)
  integer ( kind = 4 ) level_1d(dim_num)
  integer ( kind = 4 ) level_max

  do dim = 1, dim_num

    factor = 2**( level_max - level_1d(dim) )

    grid_index(dim,1:order_nd) = grid_index(dim,1:order_nd) * factor

  end do

  return
end
subroutine product_weights ( dim_num, order_1d, order_nd, rule, w_nd )

!*****************************************************************************80
!
!! PRODUCT_WEIGHTS computes the weights of a product rule.
!
!  Discussion:
!
!    This routine computes the weights for a quadrature rule which is
!    a product of 1D closed rules of varying order.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 March 2008
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
!    Input, integer ( kind = 4 ) RULE, the index of the rule.
!    1, "CC", Clenshaw Curtis Closed Fully Nested rule.
!    2, "F1", Fejer 1 Open Fully Nested rule.
!    3, "F2", Fejer 2 Open Fully Nested rule.
!    4, "GP", Gauss Patterson Open Fully Nested rule.
!    5, "GL", Gauss Legendre Open Weakly Nested rule.
!    6, "GH", Gauss Hermite Open Weakly Nested rule.
!    7, "LG", Gauss Laguerre Open Non Nested rule.
!
!    Output, real ( kind = 8 ) W_ND(DIM_NUM,ORDER_ND), the product rule weights.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) order_nd

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) order_1d(dim_num)
  integer ( kind = 4 ) rule
  real ( kind = 8 ), allocatable, dimension ( : ) :: w_1d
  real ( kind = 8 ) w_nd(order_nd)

  w_nd(1:order_nd) = 1.0D+00

  do dim = 1, dim_num

    allocate ( w_1d(1:order_1d(dim)) )

    if ( rule == 1 ) then
      call cc_weights ( order_1d(dim), w_1d )
    else if ( rule == 2 ) then
      call f1_weights ( order_1d(dim), w_1d )
    else if ( rule == 3 ) then
      call f2_weights ( order_1d(dim), w_1d )
    else if ( rule == 4 ) then
      call gp_weights ( order_1d(dim), w_1d )
    else if ( rule == 5 ) then
      call gl_weights ( order_1d(dim), w_1d )
    else if ( rule == 6 ) then
      call gh_weights ( order_1d(dim), w_1d )
    else if ( rule == 7 ) then
      call lg_weights ( order_1d(dim), w_1d )
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PRODUCT_WEIGHTS - Fatal error!'
      write ( *, '(a,i8)' ) '  Unrecognized rule number = ', rule
      stop
    end if

    call r8vec_direct_product2 ( dim, order_1d(dim), w_1d, dim_num, &
      order_nd, w_nd )

    deallocate ( w_1d )

  end do

  return
end
function r8_choose ( n, k )

!*****************************************************************************80
!
!! R8_CHOOSE computes the binomial coefficient C(N,K) as an R8.
!
!  Discussion:
!
!    The value is calculated in such a way as to avoid overflow and
!    roundoff.  The calculation is done in R8 arithmetic.
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
!    24 March 2008
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
!    Output, real ( kind = 8 ) R8_CHOOSE, the number of combinations of N
!    things taken K at a time.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) mn
  integer ( kind = 4 ) mx
  integer ( kind = 4 ) n
  real ( kind = 8 ) r8_choose
  real ( kind = 8 ) value

  mn = min ( k, n - k )

  if ( mn < 0 ) then

    value = 0.0D+00

  else if ( mn == 0 ) then

    value = 1.0D+00

  else

    mx = max ( k, n - k )
    value = real ( mx + 1, kind = 8 )

    do i = 2, mn
      value = ( value * real ( mx + i, kind = 8 ) ) / real ( i, kind = 8 )
    end do

  end if

  r8_choose = value

  return
end
function r8_factorial ( n )

!*****************************************************************************80
!
!! R8_FACTORIAL computes the factorial function.
!
!  Formula:
!
!    factorial ( N ) = N! = product ( 1 <= I <= N ) I
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the argument of the factorial function.
!    If N is less than 1, the function value is returned as 1.
!
!    Output, real ( kind = 8 ) R8_FACTORIAL, the factorial of N.
!
  implicit none

  real ( kind = 8 ) r8_factorial
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n

  r8_factorial = 1.0D+00

  do i = 1, n
    r8_factorial = r8_factorial * real ( i, kind = 8 )
  end do

  return
end
function r8_factorial2 ( n )

!*****************************************************************************80
!
!! R8_FACTORIAL2 computes the double factorial function.
!
!  Formula:
!
!    FACTORIAL2( N ) = Product ( N * (N-2) * (N-4) * ... * 2 )  (N even)
!                    = Product ( N * (N-2) * (N-4) * ... * 1 )  (N odd)
!
!  Example:
!
!     N    Factorial2(N)
!
!     0     1
!     1     1
!     2     2
!     3     3
!     4     8
!     5    15
!     6    48
!     7   105
!     8   384
!     9   945
!    10  3840
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the argument of the double factorial
!    function.  If N is less than 1, R8_FACTORIAL2 is returned as 1.0.
!
!    Output, real ( kind = 8 ) R8_FACTORIAL2, the value of the double
!    factorial of N.
!
  implicit none

  integer ( kind = 4 ) n
  real ( kind = 8 ) r8_factorial2
  real ( kind = 8 ) r8_n

  if ( n < 1 ) then
    r8_factorial2 = 1.0D+00
    return
  end if

  r8_n = real ( n, kind = 8 )
  r8_factorial2 = 1.0D+00

  do while ( 1.0D+00 < r8_n )
    r8_factorial2 = r8_factorial2 * r8_n
    r8_n = r8_n - 2.0D+00
  end do

  return
end
function r8_huge ( )

!*****************************************************************************80
!
!! R8_HUGE returns a very large R8.
!
!  Discussion:
!
!    The value returned by this function is NOT required to be the
!    maximum representable R8.  This value varies from machine to machine,
!    from compiler to compiler, and may cause problems when being printed.
!    We simply want a "very large" but non-infinite number.
!
!    FORTRAN90 provides a built-in routine HUGE ( X ) that
!    can return the maximum representable number of the same datatype
!    as X, if that is what is really desired.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) R8_HUGE, a "huge" value.
!
  implicit none

  real ( kind = 8 ) r8_huge

  r8_huge = 1.0D+30

  return
end
function r8_mop ( i )

!*****************************************************************************80
!
!! R8_MOP returns the I-th power of -1 as an R8 value.
!
!  Discussion:
!
!    An R8 is a real ( kind = 8 ) value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 November 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the power of -1.
!
!    Output, real ( kind = 8 ) R8_MOP, the I-th power of -1.
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) r8_mop

  if ( mod ( i, 2 ) == 0 ) then
    r8_mop = + 1.0D+00
  else
    r8_mop = - 1.0D+00
  end if

  return
end
subroutine r8vec_direct_product2 ( factor_index, factor_order, factor_value, &
  factor_num, point_num, w )

!*****************************************************************************80
!
!! R8VEC_DIRECT_PRODUCT2 creates a direct product of R8VEC's.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    To explain what is going on here, suppose we had to construct
!    a multidimensional quadrature rule as the product of K rules
!    for 1D quadrature.
!
!    The product rule will be represented as a list of points and weights.
!
!    The J-th item in the product rule will be associated with
!      item J1 of 1D rule 1,
!      item J2 of 1D rule 2,
!      ...,
!      item JK of 1D rule K.
!
!    In particular,
!      X(J) = ( X(1,J1), X(2,J2), ..., X(K,JK))
!    and
!      W(J) = W(1,J1) * W(2,J2) * ... * W(K,JK)
!
!    So we can construct the quadrature rule if we can properly
!    distribute the information in the 1D quadrature rules.
!
!    This routine carries out the task involving the weights W.
!
!    Another way to do this would be to compute, one by one, the
!    set of all possible indices (J1,J2,...,JK), and then index
!    the appropriate information.  An advantage of the method shown
!    here is that you can process the K-th set of information and
!    then discard it.
!
!  Example:
!
!    Rule 1:
!      Order = 4
!      W(1:4) = ( 2, 3, 5, 7 )
!
!    Rule 2:
!      Order = 3
!      W(1:3) = ( 11, 13, 17 )
!
!    Rule 3:
!      Order = 2
!      W(1:2) = ( 19, 23 )
!
!    Product Rule:
!      Order = 24
!      W(1:24) =
!        ( 2 * 11 * 19 )
!        ( 3 * 11 * 19 )
!        ( 4 * 11 * 19 )
!        ( 7 * 11 * 19 )
!        ( 2 * 13 * 19 )
!        ( 3 * 13 * 19 )
!        ( 5 * 13 * 19 )
!        ( 7 * 13 * 19 )
!        ( 2 * 17 * 19 )
!        ( 3 * 17 * 19 )
!        ( 5 * 17 * 19 )
!        ( 7 * 17 * 19 )
!        ( 2 * 11 * 23 )
!        ( 3 * 11 * 23 )
!        ( 5 * 11 * 23 )
!        ( 7 * 11 * 23 )
!        ( 2 * 13 * 23 )
!        ( 3 * 13 * 23 )
!        ( 5 * 13 * 23 )
!        ( 7 * 13 * 23 )
!        ( 2 * 17 * 23 )
!        ( 3 * 17 * 23 )
!        ( 5 * 17 * 23 )
!        ( 7 * 17 * 23 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) FACTOR_INDEX, the index of the factor being
!    processed.  The first factor processed must be factor 1!
!
!    Input, integer ( kind = 4 ) FACTOR_ORDER, the order of the factor.
!
!    Input, real ( kind = 8 ) FACTOR_VALUE(FACTOR_ORDER), the factor values
!    for factor FACTOR_INDEX.
!
!    Input, integer ( kind = 4 ) FACTOR_NUM, the number of factors.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of elements in the
!    direct product.
!
!    Input/output, real ( kind = 8 ) W(POINT_NUM), the elements of the
!    direct product, which are built up gradually.
!
!  Local Parameters:
!
!    Local, integer ( kind = 4 ) START, the first location of a block of values
!    to set.
!
!    Local, integer ( kind = 4 ) CONTIG, the number of consecutive values
!    to set.
!
!    Local, integer SKIP, the distance from the current value of START
!    to the next location of a block of values to set.
!
!    Local, integer REP, the number of blocks of values to set.
!
  implicit none

  integer ( kind = 4 ) factor_num
  integer ( kind = 4 ) factor_order
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ), save :: contig
  integer ( kind = 4 ) factor_index
  real ( kind = 8 ) factor_value(factor_order)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ), save :: rep
  integer ( kind = 4 ), save :: skip
  integer ( kind = 4 ) start
  real ( kind = 8 ) w(point_num)

  if ( factor_index == 1 ) then
    contig = 1
    skip = 1
    rep = point_num
    w(1:point_num) = 1.0D+00
  end if

  rep = rep / factor_order
  skip = skip * factor_order

  do j = 1, factor_order

    start = 1 + ( j - 1 ) * contig

    do k = 1, rep
      w(start:start+contig-1) = w(start:start+contig-1) * factor_value(j)
      start = start + skip
    end do

  end do

  contig = contig * factor_order

  return
end
subroutine sparse_grid ( dim_num, level_max, rule, point_num, &
  grid_weight, grid_point )

!****************************************************************************80
!
!! SPARSE_GRID computes a sparse grid.
!
!  Discussion:
!
!    A Smolyak construction is used to create a multidimensional sparse grid.
!
!    The user specifies:
!    * the spatial dimension of the quadrature region,
!    * the level that defines the Smolyak grid.
!    * the 1D quadrature rule.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) LEVEL_MAX, controls the size of the final
!    sparse grid.
!
!    Input, integer ( kind = 4 ) RULE, the index of the rule.
!    1, "CC", Clenshaw Curtis Closed Fully Nested rule.
!    2, "F1", Fejer 1 Open Fully Nested rule.
!    3, "F2", Fejer 2 Open Fully Nested rule.
!    4, "GP", Gauss Patterson Open Fully Nested rule.
!    5, "GL", Gauss Legendre Open Weakly Nested rule.
!    6, "GH", Gauss Hermite Open Weakly Nested rule.
!    7, "LG", Gauss Laguerre Open Non Nested rule.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points in the grid,
!    as determined by LEVELS_INDEX_SIZE.
!
!    Output, real ( kind = 8 ) GRID_WEIGHT(POINT_NUM), the weights.
!
!    Output, real ( kind = 8 ) GRID_POINT(DIM_NUM,POINT_NUM), the points.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  real ( kind = 8 ) grid_point(dim_num,point_num)
  real ( kind = 8 ) grid_weight(point_num)
  integer ( kind = 4 ) level_max
  integer ( kind = 4 ) rule

  if ( rule == 1 ) then
    call sparse_grid_cfn ( dim_num, level_max, rule, point_num, &
      grid_weight, grid_point )
  else if ( 2 <= rule .and. rule <= 4 ) then
    call sparse_grid_ofn ( dim_num, level_max, rule, point_num, &
      grid_weight, grid_point )
  else if ( 5 <= rule .and. rule <= 6 ) then
    call sparse_grid_own ( dim_num, level_max, rule, point_num, &
      grid_weight, grid_point )
  else if ( 7 == rule ) then
    call sparse_grid_onn ( dim_num, level_max, rule, point_num, &
      grid_weight, grid_point )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPARSE_GRID - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal input rule index = ', rule
    stop
  end if

  return
end
subroutine sparse_grid_cc_size ( dim_num, level_max, point_num )

!*****************************************************************************80
!
!! SPARSE_GRID_CC_SIZE sizes a sparse grid using Clenshaw Curtis rules.
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
  logical              more
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
subroutine sparse_grid_cfn ( dim_num, level_max, rule, point_num, &
  grid_weight, grid_point )

!****************************************************************************80
!
!! SPARSE_GRID_CFN computes a sparse grid based on a CFN 1D rule.
!
!  Discussion:
!
!    The 1D quadrature rule is assumed to be Closed Fully Nested.
!
!    Closed Fully Nested rules include Clenshaw Curtis rules.
!
!    A Smolyak construction is used to create a multidimensional sparse grid.
!
!    The user specifies:
!    * the spatial dimension of the quadrature region,
!    * the level that defines the Smolyak grid.
!    * the quadrature rule.
!    * the number of points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) LEVEL_MAX, controls the size of the final
!    sparse grid.
!
!    Input, integer ( kind = 4 ) RULE, the index of the rule.
!    1, "CC", Clenshaw Curtis Closed Fully Nested rule.
!    2, "F1", Fejer 1 Open Fully Nested rule.
!    3, "F2", Fejer 2 Open Fully Nested rule.
!    4, "GP", Gauss Patterson Open Fully Nested rule.
!    5, "GL", Gauss Legendre Open Weakly Nested rule.
!    6, "GH", Gauss Hermite Open Weakly Nested rule.
!    7, "LG", Gauss Laguerre Open Non Nested rule.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points in the grid,
!    as determined by LEVELS_INDEX_SIZE.
!
!    Output, real ( kind = 8 ) GRID_WEIGHT(POINT_NUM), the weights.
!
!    Output, real ( kind = 8 ) GRID_POINT(DIM_NUM,POINT_NUM), the points.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  real ( kind = 8 ) cc_abscissa
  integer ( kind = 4 ) dim
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: grid_base
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: grid_index
  real ( kind = 8 ) grid_point(dim_num,point_num)
  real ( kind = 8 ) grid_weight(point_num)
  integer ( kind = 4 ) level_max
  integer ( kind = 4 ) order_max
  integer ( kind = 4 ) point
  integer ( kind = 4 ) rule

  if ( rule /= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPARSE_GRID_CFN - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal input rule index = ', rule
    stop
  end if
!
!  Determine the index vector, relative to the full product grid,
!  that identifies the points in the sparse grid.
!
  allocate ( grid_index(dim_num,point_num) )
  allocate ( grid_base(dim_num,point_num) )

  call levels_index_cfn ( dim_num, level_max, point_num, grid_index, grid_base )
!
!  Compute the physical coordinates of the abscissas.
!
  if ( 0 == level_max ) then
    order_max = 1
  else
    order_max = 2**level_max + 1
  end if

  do point = 1, point_num
    do dim = 1, dim_num

      if ( rule == 1 ) then
        grid_point(dim,point) = &
          cc_abscissa ( order_max, grid_index(dim,point) + 1 )
      end if

    end do
  end do
!
!  Gather the weights.
!
  call sparse_grid_weights_cfn ( dim_num, level_max, rule, point_num, &
     grid_index, grid_weight )

  deallocate ( grid_base )
  deallocate ( grid_index )

  return
end
subroutine sparse_grid_ofn ( dim_num, level_max, rule, point_num, &
  grid_weight, grid_point )

!****************************************************************************80
!
!! SPARSE_GRID_OFN computes a sparse grid based on an OFN 1D rule.
!
!  Discussion:
!
!    The 1D quadrature rule is assumed to be Open Fully Nested.
!
!    Open Fully Nested rules include Fejer 1, Fejer 2, and Gauss Patterson.
!
!    A Smolyak construction is used to create a multidimensional sparse grid.
!
!    The user specifies:
!    * the spatial dimension of the quadrature region,
!    * the level that defines the Smolyak grid.
!    * the 1D quadrature rule.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) LEVEL_MAX, controls the size of the final
!    sparse grid.
!
!    Input, integer ( kind = 4 ) RULE, the index of the rule.
!    1, "CC", Clenshaw Curtis Closed Fully Nested rule.
!    2, "F1", Fejer 1 Open Fully Nested rule.
!    3, "F2", Fejer 2 Open Fully Nested rule.
!    4, "GP", Gauss Patterson Open Fully Nested rule.
!    5, "GL", Gauss Legendre Open Weakly Nested rule.
!    6, "GH", Gauss Hermite Open Weakly Nested rule.
!    7, "LG", Gauss Laguerre Open Non Nested rule.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points in the grid,
!    as determined by LEVELS_INDEX_SIZE.
!
!    Output, real ( kind = 8 ) GRID_WEIGHT(POINT_NUM), the weights.
!
!    Output, real ( kind = 8 ) GRID_POINT(DIM_NUM,POINT_NUM), the points.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) dim
  real ( kind = 8 ) f1_abscissa
  real ( kind = 8 ) f2_abscissa
  real ( kind = 8 ) gp_abscissa
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: grid_base
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: grid_index
  real ( kind = 8 ) grid_point(dim_num,point_num)
  real ( kind = 8 ) grid_weight(point_num)
  integer ( kind = 4 ) level_max
  integer ( kind = 4 ) order_max
  integer ( kind = 4 ) point
  integer ( kind = 4 ) rule

  if ( rule < 2 .or. 4 < rule ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPARSE_GRID_OFN - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal input rule index = ', rule
    stop
  end if
!
!  Determine the index vector, relative to the full product grid,
!  that identifies the points in the sparse grid.
!
  allocate ( grid_base(dim_num,point_num) )
  allocate ( grid_index(dim_num,point_num) )

  call levels_index_ofn ( dim_num, level_max, point_num, grid_index, grid_base )
!
!  Compute the physical coordinates of the abscissas.
!
  order_max = 2**( level_max + 1 ) - 1

  do point = 1, point_num
    do dim = 1, dim_num

      if ( rule == 2 ) then
        grid_point(dim,point) = &
          f1_abscissa ( order_max, grid_index(dim,point) )
      else if ( rule == 3 ) then
        grid_point(dim,point) = &
          f2_abscissa ( order_max, grid_index(dim,point) )
      else if ( rule == 4 ) then
        grid_point(dim,point) = &
          gp_abscissa ( order_max, grid_index(dim,point) )
      end if

    end do
  end do
!
!  Gather the weights.
!
  call sparse_grid_weights_ofn ( dim_num, level_max, rule, point_num, &
    grid_index, grid_weight )

  deallocate ( grid_base )
  deallocate ( grid_index )

  return
end
subroutine sparse_grid_ofn_size ( dim_num, level_max, point_num )

!*****************************************************************************80
!
!! SPARSE_GRID_OFN_SIZE sizes a sparse grid using Open Fully Nested rules.
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
  logical              more
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
subroutine sparse_grid_onn ( dim_num, level_max, rule, point_num, &
  grid_weight, grid_point )

!****************************************************************************80
!
!! SPARSE_GRID_ONN computes a sparse grid based on a ONN 1D rule.
!
!  Discussion:
!
!    The 1D quadrature rule is assumed to be Open Non-Nested.
!    Such rules include Gauss Laguerre rules.
!
!    A Smolyak construction is used to create a multidimensional sparse grid.
!
!    The user specifies:
!    * the spatial dimension of the quadrature region,
!    * the level that defines the Smolyak grid.
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
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) LEVEL_MAX, controls the size of the
!    sparse grid.
!
!    Input, integer ( kind = 4 ) RULE, the index of the rule.
!    1, "CC", Clenshaw Curtis Closed Fully Nested rule.
!    2, "F1", Fejer 1 Open Fully Nested rule.
!    3, "F2", Fejer 2 Open Fully Nested rule.
!    4, "GP", Gauss Patterson Open Fully Nested rule.
!    5, "GL", Gauss Legendre Open Weakly Nested rule.
!    6, "GH", Gauss Hermite Open Weakly Nested rule.
!    7, "LG", Gauss Laguerre Open Non Nested rule.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points in the grid,
!    as determined by LEVELS_ONN_INDEX_SIZE.
!
!    Output, real ( kind = 8 ) GRID_WEIGHT(POINT_NUM), the weights.
!
!    Output, real ( kind = 8 ) GRID_POINT(DIM_NUM,POINT_NUM), the points.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  real ( kind = 8 ) coeff
  integer ( kind = 4 ), dimension ( dim_num ) :: grid_base2
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: grid_index2
  integer ( kind = 4 ), allocatable, dimension ( : ) :: grid_level
  real ( kind = 8 ) grid_point(dim_num,point_num)
  real ( kind = 8 ) grid_weight(point_num)
  real ( kind = 8 ), allocatable, dimension ( : ) :: grid_weight2
  integer ( kind = 4 ) h
  integer ( kind = 4 ), parameter :: i4_1 = 1
  integer ( kind = 4 ) level
  integer ( kind = 4 ), dimension ( dim_num ) :: level_1d
  integer ( kind = 4 ) level_max
  integer ( kind = 4 ) level_min
  logical more
  integer ( kind = 4 ), dimension ( dim_num ) :: order_1d
  integer ( kind = 4 ) order_nd
  integer ( kind = 4 ) point
  integer ( kind = 4 ) point_num2
  real ( kind = 8 ) r8_choose
  real ( kind = 8 ) r8_mop
  integer ( kind = 4 ) rule
  integer ( kind = 4 ) t

  grid_weight(1:point_num) = 0.0D+00
!
!  The outer loop generates LEVELs from LEVEL_MIN to LEVEL_MAX.
!
  point_num2 = 0

  level_min = max ( 0, level_max + 1 - dim_num )

  do level = level_min, level_max
!
!  The middle loop generates the next partition LEVEL_1D(1:DIM_NUM)
!  that adds up to LEVEL.
!
    more = .false.
    h = 0
    t = 0

    do

      call comp_next ( level, dim_num, level_1d, more, h, t )
!
!  Transform each 1D level to a corresponding 1D order.
!  The relationship is the same as for other OPEN rules.
!
      call level_to_order_open ( dim_num, level_1d, order_1d )

      grid_base2(1:dim_num) = order_1d(1:dim_num)
!
!  The product of the 1D orders gives us the number of points in this subgrid.
!
      order_nd = product ( order_1d(1:dim_num) )

      allocate ( grid_index2(1:dim_num,1:order_nd) )
      allocate ( grid_level(1:order_nd) )
      allocate ( grid_weight2(1:order_nd) )
!
!  Compute the weights for this product grid.
!
      call product_weights ( dim_num, order_1d, order_nd, rule, grid_weight2 )
!
!  Now determine the coefficient of the weight.
!
      coeff = r8_mop ( level_max - level ) &
        * r8_choose ( dim_num - 1, level_max - level )
!
!  The inner (hidden) loop generates all points corresponding to given grid.
!  The grid indices will be between 1 and ORDER_1D(DIM).
!
      call multigrid_index_onn ( dim_num, order_1d, order_nd, grid_index2 )

      do point = 1, order_nd

        point_num2 = point_num2 + 1

        if ( point_num < point_num2 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'SPARSE_GRID_ONN - Fatal error!'
          write ( *, '(a,i8)' ) &
          '  Exceeding maximum point index POINT_NUM = ', point_num
          stop
        end if

        call lg_abscissa ( dim_num, i4_1, grid_index2(1:dim_num,point), &
          grid_base2(1:dim_num), grid_point(1:dim_num,point_num2) )

        grid_weight(point_num2) = coeff * grid_weight2(point)

      end do

      deallocate ( grid_index2 )
      deallocate ( grid_level )
      deallocate ( grid_weight2 )

      if ( .not. more ) then
        exit
      end if

    end do

  end do

  if ( point_num2 < point_num ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPARSE_GRID_ONN - Fatal error!'
    write ( *, '(a,i8)' ) '  Set fewer points than POINT_NUM = ', point_num
    stop
  end if

  return
end
subroutine sparse_grid_own ( dim_num, level_max, rule, point_num, &
  grid_weight, grid_point )

!****************************************************************************80
!
!! SPARSE_GRID_OWN computes a sparse grid based on an OWN 1D rule.
!
!  Discussion:
!
!    The 1D quadrature rule is assumed to be Open Weakly Nested.
!    Such rules include Gauss Hermite and Gauss Legendre rules.
!
!    A Smolyak construction is used to create a multidimensional sparse grid.
!
!    The user specifies:
!    * the spatial dimension of the quadrature region,
!    * the level that defines the Smolyak grid,
!    * the rule;
!    * the number of points.
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
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) LEVEL_MAX, controls the size of the
!    sparse grid.
!
!    Input, integer ( kind = 4 ) RULE, the index of the rule.
!    1, "CC", Clenshaw Curtis Closed Fully Nested rule.
!    2, "F1", Fejer 1 Open Fully Nested rule.
!    3, "F2", Fejer 2 Open Fully Nested rule.
!    4, "GP", Gauss Patterson Open Fully Nested rule.
!    5, "GL", Gauss Legendre Open Weakly Nested rule.
!    6, "GH", Gauss Hermite Open Weakly Nested rule.
!    7, "LG", Gauss Laguerre Open Non Nested rule.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points in the grid,
!    as determined by LEVELS_INDEX_SIZE_OWN.
!
!    Output, real ( kind = 8 ) GRID_WEIGHT(POINT_NUM), the weights.
!
!    Output, real ( kind = 8 ) GRID_POINT(DIM_NUM,POINT_NUM), the points.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  real ( kind = 8 ) coeff
  integer ( kind = 4 ), dimension ( dim_num ) :: grid_base2
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: grid_index2
  integer ( kind = 4 ), allocatable, dimension ( : ) :: grid_level
  real ( kind = 8 ) grid_point(dim_num,point_num)
  real ( kind = 8 ) grid_point_temp(dim_num)
  real ( kind = 8 ) grid_weight(point_num)
  real ( kind = 8 ), allocatable, dimension ( : ) :: grid_weight2
  integer ( kind = 4 ) h
  integer ( kind = 4 ), parameter :: i4_1 = 1
  integer ( kind = 4 ) level
  integer ( kind = 4 ), dimension ( dim_num ) :: level_1d
  integer ( kind = 4 ) level_max
  integer ( kind = 4 ) level_min
  integer ( kind = 4 ) level_min2
  logical more
  integer ( kind = 4 ), dimension ( dim_num ) :: order_1d
  integer ( kind = 4 ) order_nd
  integer ( kind = 4 ) point
  integer ( kind = 4 ) point2
  integer ( kind = 4 ) point3
  integer ( kind = 4 ) point_num2
  real ( kind = 8 ) r8_choose
  real ( kind = 8 ) r8_mop
  integer ( kind = 4 ) rule
  integer ( kind = 4 ) t

  grid_weight(1:point_num) = 0.0D+00
!
!  The outer loop generates LEVELs from LEVEL_MIN to LEVEL_MAX.
!
  point_num2 = 0

  level_min = max ( 0, level_max + 1 - dim_num )

  if ( dim_num == 1 ) then
    level_min2 = level_min
  else
    level_min2 = 0
  end if

  do level = level_min2, level_max
!
!  The middle loop generates the next partition LEVEL_1D(1:DIM_NUM)
!  that adds up to LEVEL.
!
    more = .false.
    h = 0
    t = 0

    do

      call comp_next ( level, dim_num, level_1d, more, h, t )
!
!  Transform each 1D level to a corresponding 1D order.
!  The relationship is the same as for other OPEN rules.
!  The GL rule differs from the other OPEN rules only in the nesting behavior.
!
      call level_to_order_open ( dim_num, level_1d, order_1d )

      grid_base2(1:dim_num) = ( order_1d(1:dim_num) - 1 ) / 2
!
!  The product of the 1D orders gives us the number of points in this subgrid.
!
      order_nd = product ( order_1d(1:dim_num) )

      allocate ( grid_index2(1:dim_num,1:order_nd) )
      allocate ( grid_level(1:order_nd) )
      allocate ( grid_weight2(1:order_nd) )
!
!  Compute the weights for this product grid.
!
      call product_weights ( dim_num, order_1d, order_nd, rule, grid_weight2 )
!
!  Now determine the coefficient of the weight.
!
      coeff = r8_mop ( level_max - level ) &
        * r8_choose ( dim_num - 1, level_max - level )
!
!  The inner (hidden) loop generates all points corresponding to given grid.
!  The grid indices will be between -M to +M, where 2*M + 1 = ORDER_1D(DIM).
!
      call multigrid_index_own ( dim_num, order_1d, order_nd, grid_index2 )
!
!  Determine the first level of appearance of each of the points.
!  This allows us to flag certain points as being repeats of points
!  generated on a grid of lower level.
!
!  This is SLIGHTLY tricky.
!
      call index_level_own ( level, level_max, dim_num, order_nd, grid_index2, &
        grid_base2, grid_level )
!
!  Only keep those points which first appear on this level.
!
      do point = 1, order_nd
!
!  Either a "new" point (increase count, create point, create weight)
!
        if ( grid_level(point) == level ) then

          point_num2 = point_num2 + 1

          if ( point_num < point_num2 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'SPARSE_GRID_OWN - Fatal error!'
            write ( *, '(a,i8)' ) &
            '  Exceeding maximum point index POINT_NUM = ', point_num
            stop
          end if

          if ( rule == 5 ) then
            call gl_abscissa ( dim_num, i4_1, grid_index2(1:dim_num,point), &
              grid_base2(1:dim_num), grid_point(1:dim_num,point_num2) )
          else if ( rule == 6 ) then
            call gh_abscissa ( dim_num, i4_1, grid_index2(1:dim_num,point), &
              grid_base2(1:dim_num), grid_point(1:dim_num,point_num2) )
          else
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'SPARSE_GRID_OWN - Fatal error!'
            write ( *, '(a,i8)' ) '  Unrecognized rule number = ', rule
            stop
          end if

          if ( level_min <= level ) then
            grid_weight(point_num2) = coeff * grid_weight2(point)
          end if
!
!  or an already existing point (create point temporarily, find match,
!  add weight to matched point's weight).
!
        else

          if ( level_min <= level ) then

            if ( rule == 5 ) then
              call gl_abscissa ( dim_num, i4_1, grid_index2(1:dim_num,point), &
                grid_base2(1:dim_num), grid_point_temp(1:dim_num) )
            else if ( rule == 6 ) then
              call gh_abscissa ( dim_num, i4_1, grid_index2(1:dim_num,point), &
                grid_base2(1:dim_num), grid_point_temp(1:dim_num) )
            else
              write ( *, '(a)' ) ' '
              write ( *, '(a)' ) 'SPARSE_GRID_OWN - Fatal error!'
              write ( *, '(a,i8)' ) '  Unrecognized rule number = ', rule
              stop
            end if

            point3 = -1

            do point2 = 1, point_num2
              if ( all ( grid_point(1:dim_num,point2) == &
                         grid_point_temp(1:dim_num) ) ) then
                point3 = point2
                exit
              end if
            end do

            if ( point3 == -1 ) then
              write ( *, '(a)' ) ' '
              write ( *, '(a)' ) 'SPARSE_GRID_OWN - Fatal error!'
              write ( *, '(a)' ) '  Could not match point.'
              stop
            end if

            grid_weight(point3) = grid_weight(point3) + &
              coeff * grid_weight2(point)

          end if

        end if

      end do

      deallocate ( grid_index2 )
      deallocate ( grid_level )
      deallocate ( grid_weight2 )

      if ( .not. more ) then
        exit
      end if

    end do

  end do

  if ( point_num2 < point_num ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPARSE_GRID_OWN - Fatal error!'
    write ( *, '(a,i8)' ) '  Set fewer points than POINT_NUM = ', point_num
    stop
  end if

  return
end
subroutine sparse_grid_weights_cfn ( dim_num, level_max, rule, point_num, &
  grid_index, grid_weight )

!*****************************************************************************80
!
!! SPARSE_GRID_WEIGHTS_CFN computes sparse grid weights based on a CFN 1D rule.
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
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) LEVEL_MAX, the maximum value of LEVEL.
!
!    Input, integer ( kind = 4 ) RULE, the index of the rule.
!    1, "CC", Clenshaw Curtis Closed Fully Nested rule.
!    2, "F1", Fejer 1 Open Fully Nested rule.
!    3, "F2", Fejer 2 Open Fully Nested rule.
!    4, "GP", Gauss Patterson Open Fully Nested rule.
!    5, "GL", Gauss Legendre Open Weakly Nested rule.
!    6, "GH", Gauss Hermite Open Weakly Nested rule.
!    7, "LG", Gauss Laguerre Open Non Nested rule.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the total number of points in
!    the grids.
!
!    Input, integer ( kind = 4 ) GRID_INDEX(DIM_NUM,POINT_NUM), a list of
!    point indices, representing a subset of the product grid of level
!    LEVEL_MAX, representing (exactly once) each point that will show up in a
!    sparse grid of level LEVEL_MAX.
!
!    Output, real ( kind = 8 ) GRID_WEIGHT(POINT_NUM), the weights
!    associated with the sparse grid points.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  real ( kind = 8 ) coeff
  integer ( kind = 4 ) grid_index(dim_num,point_num)
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: grid_index2
  real ( kind = 8 ) grid_weight(point_num)
  real ( kind = 8 ), allocatable, dimension ( : ) :: grid_weight2
  integer ( kind = 4 ) h
  integer ( kind = 4 ) level
  integer ( kind = 4 ) level_1d(dim_num)
  integer ( kind = 4 ) level_max
  integer ( kind = 4 ) level_min
  logical more
  integer ( kind = 4 ) order_nd
  integer ( kind = 4 ) order_1d(dim_num)
  integer ( kind = 4 ) point
  integer ( kind = 4 ) point2
  real ( kind = 8 ) r8_choose
  real ( kind = 8 ) r8_mop
  integer ( kind = 4 ) rule
  integer ( kind = 4 ) t

  if ( level_max == 0 ) then
    grid_weight(1:point_num) = 2.0D+00**dim_num
    return
  end if

  grid_weight(1:point_num) = 0.0D+00

  level_min = max ( 0, level_max + 1 - dim_num )

  do level = level_min, level_max
!
!  The middle loop generates the next partition LEVEL_1D(1:DIM_NUM)
!  that adds up to LEVEL.
!
    more = .false.
    h = 0
    t = 0

    do

      call comp_next ( level, dim_num, level_1d, more, h, t )
!
!  Transform each 1D level to a corresponding 1D order.
!
      call level_to_order_closed ( dim_num, level_1d, order_1d )
!
!  The product of the 1D orders gives us the number of points in this grid.
!
      order_nd = product ( order_1d(1:dim_num) )

      allocate ( grid_index2(1:dim_num,1:order_nd) )
      allocate ( grid_weight2(1:order_nd) )
!
!  Generate the indices of the points corresponding to the grid.
!
      call multigrid_index_cfn ( dim_num, order_1d, order_nd, grid_index2 )
!
!  Compute the weights for this grid.
!
      call product_weights ( dim_num, order_1d, order_nd, rule, grid_weight2 )
!
!  Adjust the grid indices to reflect LEVEL_MAX.
!
      call multigrid_scale_closed ( dim_num, order_nd, level_max, level_1d, &
        grid_index2 )
!
!  Now determine the coefficient.
!
      coeff = r8_mop ( level_max - level ) &
        * r8_choose ( dim_num - 1, level_max - level )

      do point2 = 1, order_nd

        do point = 1, point_num

          if ( all ( &
            grid_index2(1:dim_num,point2) == grid_index(1:dim_num,point) &
          ) ) then
            grid_weight(point) = grid_weight(point) &
              + coeff * grid_weight2(point2)
            exit
          end if

        end do

      end do

      deallocate ( grid_index2 )
      deallocate ( grid_weight2 )

      if ( .not. more ) then
        exit
      end if

    end do

  end do

  return
end
subroutine sparse_grid_weights_ofn ( dim_num, level_max, rule, point_num, &
  grid_index, grid_weight )

!*****************************************************************************80
!
!! SPARSE_GRID_WEIGHTS_OFN computes sparse grid weights based on a OFN 1D rule.
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
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) LEVEL_MAX, the maximum value of LEVEL.
!
!    Input, integer ( kind = 4 ) RULE, the index of the rule.
!    1, "CC", Clenshaw Curtis Closed Fully Nested rule.
!    2, "F1", Fejer 1 Open Fully Nested rule.
!    3, "F2", Fejer 2 Open Fully Nested rule.
!    4, "GP", Gauss Patterson Open Fully Nested rule.
!    5, "GL", Gauss Legendre Open Weakly Nested rule.
!    6, "GH", Gauss Hermite Open Weakly Nested rule.
!    7, "LG", Gauss Laguerre Open Non Nested rule.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the total number of points in
!    the grids.
!
!    Input, integer ( kind = 4 ) GRID_INDEX(DIM_NUM,POINT_NUM), a list of
!    point indices, representing a subset of the product grid of level
!    LEVEL_MAX, representing (exactly once) each point that will show up in a
!    sparse grid of level LEVEL_MAX.
!
!    Output, real ( kind = 8 ) GRID_WEIGHT(POINT_NUM), the weights
!    associated with the sparse grid points.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  real ( kind = 8 ) coeff
  integer ( kind = 4 ) grid_index(dim_num,point_num)
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: grid_index2
  real ( kind = 8 ) grid_weight(point_num)
  real ( kind = 8 ), allocatable, dimension ( : ) :: grid_weight2
  integer ( kind = 4 ) h
  integer ( kind = 4 ) level
  integer ( kind = 4 ) level_1d(dim_num)
  integer ( kind = 4 ) level_max
  integer ( kind = 4 ) level_min
  logical more
  integer ( kind = 4 ) order_nd
  integer ( kind = 4 ) order_1d(dim_num)
  integer ( kind = 4 ) point
  integer ( kind = 4 ) point2
  real ( kind = 8 ) r8_choose
  real ( kind = 8 ) r8_mop
  integer ( kind = 4 ) rule
  integer ( kind = 4 ) t

  if ( level_max == 0 ) then
    grid_weight(1:point_num) = 2.0D+00**dim_num
    return
  end if

  grid_weight(1:point_num) = 0.0D+00

  level_min = max ( 0, level_max + 1 - dim_num )

  do level = level_min, level_max
!
!  The middle loop generates the next partition LEVEL_1D(1:DIM_NUM)
!  that adds up to LEVEL.
!
    more = .false.
    h = 0
    t = 0

    do

      call comp_next ( level, dim_num, level_1d, more, h, t )
!
!  Transform each 1D level to a corresponding 1D order.
!
      call level_to_order_open ( dim_num, level_1d, order_1d )
!
!  The product of the 1D orders gives us the number of points in this grid.
!
      order_nd = product ( order_1d(1:dim_num) )

      allocate ( grid_index2(1:dim_num,1:order_nd) )
      allocate ( grid_weight2(1:order_nd) )
!
!  Generate the indices of the points corresponding to the grid.
!
      call multigrid_index_ofn ( dim_num, order_1d, order_nd, grid_index2 )
!
!  Compute the weights for this grid.
!
      call product_weights ( dim_num, order_1d, order_nd, rule, grid_weight2 )
!
!  Adjust the grid indices to reflect LEVEL_MAX.
!
      call multigrid_scale_open ( dim_num, order_nd, level_max, level_1d, &
        grid_index2 )
!
!  Now determine the coefficient.
!
      coeff = r8_mop ( level_max - level ) &
        * r8_choose ( dim_num - 1, level_max - level )

      do point2 = 1, order_nd

        do point = 1, point_num

          if ( all ( &
            grid_index2(1:dim_num,point2) == grid_index(1:dim_num,point) &
          ) ) then
            grid_weight(point) = grid_weight(point) &
              + coeff * grid_weight2(point2)
            exit
          end if

        end do

      end do

      deallocate ( grid_index2 )
      deallocate ( grid_weight2 )

      if ( .not. more ) then
        exit
      end if

    end do

  end do

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
!    08 March 2008
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
subroutine vec_colex_next2 ( dim_num, base, a, more )

!*****************************************************************************80
!
!! VEC_COLEX_NEXT2 generates vectors in colex order.
!
!  Discussion:
!
!    The vectors are produced in colexical order, starting with
!
!    (0,        0,        ...,0),
!    (1,        0,        ...,0),
!     ...
!    (BASE(1)-1,0,        ...,0)
!
!    (0,        1,        ...,0)
!    (1,        1,        ...,0)
!    ...
!    (BASE(1)-1,1,        ...,0)
!
!    (0,        2,        ...,0)
!    (1,        2,        ...,0)
!    ...
!    (BASE(1)-1,BASE(2)-1,...,BASE(DIM_NUM)-1).
!
!  Examples:
!
!    DIM_NUM = 2,
!    BASE = ( 3, 3 )
!
!    0   0
!    1   0
!    2   0
!    0   1
!    1   1
!    2   1
!    0   2
!    1   2
!    2   2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 March 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) BASE(DIM_NUM), the bases to be used in each
!    dimension.  In dimension I, entries will range from 0 to BASE(I)-1.
!
!    Input/output, integer ( kind = 4 ) A(DIM_NUM).  On each return, A
!    will contain entries in the range 0 to N-1.
!
!    Input/output, logical MORE.  Set this variable FALSE before
!    the first call.  On return, MORE is TRUE if another vector has
!    been computed.  If MORE is returned FALSE, ignore the output
!    vector and stop calling the routine.
!
  implicit none

  integer ( kind = 4 ) dim_num

  integer ( kind = 4 ) a(dim_num)
  integer ( kind = 4 ) base(dim_num)
  integer ( kind = 4 ) i
  logical more

  if ( .not. more ) then

    a(1:dim_num) = 0
    more = .true.

  else

    do i = 1, dim_num

      a(i) = a(i) + 1

      if ( a(i) < base(i) ) then
        return
      end if

      a(i) = 0

    end do

    more = .false.

  end if

  return
end
