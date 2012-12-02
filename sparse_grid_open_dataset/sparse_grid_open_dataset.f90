 program main

!*****************************************************************************80
!
!! MAIN is the main program for SPARSE_GRID_OPEN_DATASET.
!
!  Discussion:
!
!    This program computes a quadrature rule and writes it to a file.
!
!    The quadrature rule is associated with a sparse grid derived from
!    a Smolyak construction using a open 1D quadrature rule.
!
!    The user specifies:
!    * DIM_NUM, the spatial dimension of the quadrature region,
!    * LEVEL_MAX, the level that defines the Smolyak grid.
!    * RULE, to identify the open 1D quadrature rule to be used.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 August 2009
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

  integer   ( kind = 4 ) arg_num
  integer   ( kind = 4 ) dim
  integer   ( kind = 4 ) dim_num
  character ( len = 2 )  dim_num_string
  real ( kind = 8 ) f2_abscissa
  real ( kind = 8 ) gp_abscissa
  integer   ( kind = 4 ), allocatable, dimension ( :, : ) :: grid_index
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: grid_point
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: grid_region
  real ( kind = 8 ), allocatable, dimension ( : ) :: grid_weight
  real ( kind = 8 ) h
  integer   ( kind = 4 ), parameter :: i4_1 = 1
  integer   ( kind = 4 ), parameter :: i4_10 = 10
  integer   ( kind = 4 ) iarg
  integer   ( kind = 4 ) iargc
  integer   ( kind = 4 ) ierror
  integer   ( kind = 4 ) last
  integer   ( kind = 4 ) level_max
  character ( len = 2 )  level_max_string
  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) n
  real ( kind = 8 ) nco_abscissa
  integer   ( kind = 4 ) order_max
  integer   ( kind = 4 ) point
  integer   ( kind = 4 ) point_num
  integer   ( kind = 4 ) rule
  character ( len = 4 )  rule_string
  character ( len = 255 ) string
  character ( len = 255 ) r_filename
  real ( kind = 8 ) ts_abscissa
  character ( len = 255 ) w_filename
  real ( kind = 8 ) weight_sum
  character ( len = 255 ) x_filename

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPARSE_GRID_OPEN_DATASET'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Compute the abscissas and weights of a quadrature rule'
  write ( *, '(a)' ) '  associated with a sparse grid derived from a Smolyak'
  write ( *, '(a)' ) '  construction based on an open 1D quadrature rule.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Inputs to the program include:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    DIM_NUM, the spatial dimension.'
  write ( *, '(a)' ) '    (typically in the range of 2 to 10)'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    LEVEL_MAX, the "level" of the sparse grid.'
  write ( *, '(a)' ) '    (typically in the range of 0, 1, 2, 3, ...'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    RULE, the 1D quadrature rule'
  write ( *, '(a)' ) '    2: Fejer Type 2 ("F2"), on (-1,1),'
  write ( *, '(a)' ) '    3: Gauss-Patterson ("GP"), on (-1,1),'
  write ( *, '(a)' ) '    4: Newton-Cotes Open ("NCO"), on (-1,1),'
  write ( *, '(a)' ) '    5: Tanh-Sinh ("TS"), on (-1,1),'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Output from the program includes:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    A printed table of the abscissas and weights.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    A set of files defining the quadrature rule:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    "***_d?_level?_x.txt", a file of the abscissas;'
  write ( *, '(a)' ) '    "***_d?_level?_w.txt", a file of the weights;'
  write ( *, '(a)' ) '    "***_d?_level?_r.txt", a file of the ranges.'
!
!  Get the number of command line arguments.
!
  arg_num = iargc ( )
!
!  Get the spatial dimension.
!
  if ( 1 <= arg_num ) then
    iarg = 1
    call getarg ( iarg, string )
    call s_to_i4 ( string, dim_num, ierror, last )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter the value of DIM_NUM (1 or more)'
    read ( *, * ) dim_num
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Spatial dimension requested is = ', dim_num
!
!  Get the level.
!
  if ( 2 <= arg_num ) then
    iarg = 2
    call getarg ( iarg, string )
    call s_to_i4 ( string, level_max, ierror, last )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter the value of LEVEL_MAX (0 or more).'
    read ( *, * ) level_max
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The sparse grid level is = ', level_max
!
!  Get the rule.
!
  if ( 3 <= arg_num ) then
    iarg = 3
    call getarg ( iarg, string )
    call s_to_i4 ( string, rule, ierror, last )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter the value of RULE.'
    write ( *, '(a)' ) '  2 = F2   = Fejer Type 2 Rule,'
    write ( *, '(a)' ) '  3 = GP   = Gauss-Patterson,'
    write ( *, '(a)' ) '  4 = NCO  = Newton-Cotes Open,'
    write ( *, '(a)' ) '  5 = TS   = Tanh-Sinh;'
    read ( *, * ) rule
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The 1D quadrature rule index = ', rule
  if ( rule == 2 ) then
    write ( *, '(a)' ) '  F2:   Fejer Type 2 Rule'
  else if ( rule == 3 ) then
    write ( *, '(a)' ) '  GP:   Gauss-Patterson Rule'
  else if ( rule == 4 ) then
    write ( *, '(a)' ) '  NCO:  Newton-Cotes Open Rule'
  else if ( rule == 5 ) then
    write ( *, '(a)' ) '  TS:   Tanh-Sinh Rule'
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPARSE_GRID_OPEN_DATASET - Fatal error!'
    write ( *, '(a)' ) '  Illegal value of RULE.'
    stop
  end if
!
!  How many distinct points will there be?
!
  call sparse_grid_ofn_size ( dim_num, level_max, point_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The number of distinct abscissas in the '
  write ( *, '(a)' ) '  quadrature rule is determined from the spatial'
  write ( *, '(a)' ) '  dimension DIM_NUM and the level LEVEL_MAX.'
  write ( *, '(a,i8)' ) &
    '  For the given input, this value will be = ', point_num
!
!  Allocate memory.
!
  allocate ( grid_index(dim_num,point_num) )
  allocate ( grid_point(1:dim_num,1:point_num) )
  allocate ( grid_weight(1:point_num) )
  allocate ( grid_region(1:dim_num,2) )
!
!  Determine the index vector, relative to the full product grid,
!  that identifies the points in the sparse grid.
!
  call spgrid_open_index ( dim_num, level_max, point_num, grid_index )

  call i4mat_transpose_print_some ( dim_num, point_num, grid_index, &
    i4_1, i4_1, dim_num, i4_10, '  First 10 entries of grid index:' )
!
!  Compute the physical coordinates of the abscissas.
!
  order_max = 2**( level_max + 1 ) - 1

  if ( rule == 5 ) then

    m = level_max - 3
    n = ( ( order_max + 1 ) / 2 ) - 1
    h = 4.0D+00 / real ( order_max + 1, kind = 8 )

    write ( *, '(a,i8,a,i8,a,i8,a,g14.6)' ) &
      '  M = ', m, &
      '  ORDER_MAX = ', order_max, &
      '  N = ', n, &
      '  H = ', h

  end if

  if ( rule == 2 ) then

    do point = 1, point_num
      do dim = 1, dim_num
        grid_point(dim,point) = &
          f2_abscissa ( order_max, grid_index(dim,point) )
      end do
    end do

  else if ( rule == 3 ) then

    do point = 1, point_num
      do dim = 1, dim_num
        grid_point(dim,point) = &
          gp_abscissa ( order_max, grid_index(dim,point) )
      end do
    end do

  else if ( rule == 4 ) then

    do point = 1, point_num
      do dim = 1, dim_num
        grid_point(dim,point) = &
          nco_abscissa ( order_max, grid_index(dim,point) )
      end do
    end do

  else if ( rule == 5 ) then

    do point = 1, point_num
      do dim = 1, dim_num
        grid_point(dim,point) = &
          ts_abscissa ( order_max, grid_index(dim,point) )
      end do
    end do

  end if

  call r8mat_transpose_print_some ( dim_num, point_num, grid_point, &
    i4_1, i4_1, dim_num, i4_10, '  First 10 entries of grid points:' )
!
!  Gather the weights.
!
  call spgrid_open_weights ( dim_num, level_max, point_num, grid_index, &
    rule, grid_weight )

  call r8vec_print_some ( point_num, grid_weight, i4_1, i4_10, &
    '  First 10 entries of grid weights:' )

  weight_sum = sum ( grid_weight(1:point_num) )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g24.16)' ) '  Weights sum to   ', weight_sum
  write ( *, '(a,g24.16)' ) '  Correct value is ', 2.0**dim_num
!
!  Construct appropriate file names.
!
  if ( rule == 2 ) then
    rule_string = 'f2'
  else if ( rule == 3 ) then
    rule_string = 'gp'
  else if ( rule == 4 ) then
    rule_string = 'nco'
  else if ( rule == 5 ) then
    rule_string = 'ts'
  end if

  if ( level_max < 10 ) then
    write ( level_max_string, '(i1,a)' ) level_max, ' '
  else
    write ( level_max_string, '(i2)' ) level_max
  end if

  if ( dim_num < 10 ) then
    write ( dim_num_string, '(i1,a)' ) dim_num, ' '
  else
    write ( dim_num_string, '(i2)' ) dim_num
  end if

  x_filename = trim ( rule_string )      // '_d' // &
               trim ( dim_num_string )   // '_level' // &
               trim ( level_max_string ) // '_x.txt'

  w_filename = trim ( rule_string )      // '_d' // &
               trim ( dim_num_string )   // '_level' // &
               trim ( level_max_string ) // '_w.txt'

  r_filename = trim ( rule_string )      // '_d' // &
               trim ( dim_num_string )   // '_level' // &
               trim ( level_max_string ) // '_r.txt'
!
!  Write the data to files.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Creating X file = "' // trim ( x_filename ) // '".'

  call r8mat_write ( x_filename, dim_num, point_num, grid_point )

  write ( *, '(a)' ) '  Creating W file = "' // trim ( w_filename ) // '".'

  call r8mat_write ( w_filename, 1, point_num, grid_weight )

  write ( *, '(a)' ) '  Creating R file = "' // trim ( r_filename ) // '".'

  grid_region(1:dim_num,1) = -1.0D+00
  grid_region(1:dim_num,2) = +1.0D+00

  call r8mat_write ( r_filename, dim_num, 2, grid_region )

  deallocate ( grid_index )
  deallocate ( grid_point )
  deallocate ( grid_region )
  deallocate ( grid_weight )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPARSE_GRID_OPEN_DATASET'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
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
!    Input/output, integer ( kind = 4 )  H, T, two internal parameters needed for the
!    computation.  The user should allocate space for these in the calling
!    program, include them in the calling sequence, but never alter them!
!
  implicit none

  integer ( kind = 4 ) k

  integer ( kind = 4 ) a(k)
  integer ( kind = 4 ) h
  logical              more
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
!    09 March 2008
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
  real ( kind = 8 ) :: x_max = +1.0D+00
  real ( kind = 8 ) :: x_min = -1.0D+00

  if ( order < 1 ) then
    f2_abscissa = - huge ( f2_abscissa )
    return
  end if

  if ( i < 1 .or. order < i ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'F2_ABSCISSA - Fatal error!'
    write ( *, '(a)' ) '  1 <= I <= ORDER is required!'
    stop
  end if

  if ( order == 1 ) then
    f2_abscissa = 0.0D+00
    return
  end if

  f2_abscissa = cos ( real ( order + 1 - i, kind = 8 ) * pi &
                    / real ( order + 1,     kind = 8 ) )

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
!    09 March 2008
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
!! GP_WEIGHTS sets weights for a Gauss-Patterson rule.
!
!  Discussion:
!
!    The zeroth rule, of order 1, is the standard Gauss-Legendre rule.
!
!    The first rule, of order 3, is the standard Gauss-Legendre rule.
!
!    The second rule, of order 7, includes the abscissas of the previous
!    rule.
!
!    Each subsequent rule is nested in a similar way.  Rules are available
!    of orders 1, 3, 7, 15, 31, 63, 127 and 255.
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

    w(1) = 0.69379364324108267170D-05
    w(2) = 0.25157870384280661489D-04
    w(3) = 0.53275293669780613125D-04
    w(4) = 0.90372734658751149261D-04
    w(5) = 0.13575491094922871973D-03
    w(6) = 0.18887326450650491366D-03
    w(7) = 0.24921240048299729402D-03
    w(8) = 0.31630366082226447689D-03
    w(9) = 0.38974528447328229322D-03
    w(10) = 0.46918492424785040975D-03
    w(11) = 0.55429531493037471492D-03
    w(12) = 0.64476204130572477933D-03
    w(13) = 0.74028280424450333046D-03
    w(14) = 0.84057143271072246365D-03
    w(15) = 0.94536151685852538246D-03
    w(16) = 0.10544076228633167722D-02
    w(17) = 0.11674841174299594077D-02
    w(18) = 0.12843824718970101768D-02
    w(19) = 0.14049079956551446427D-02
    w(20) = 0.15288767050877655684D-02
    w(21) = 0.16561127281544526052D-02
    w(22) = 0.17864463917586498247D-02
    w(23) = 0.19197129710138724125D-02
    w(24) = 0.20557519893273465236D-02
    w(25) = 0.21944069253638388388D-02
    w(26) = 0.23355251860571608737D-02
    w(27) = 0.24789582266575679307D-02
    w(28) = 0.26245617274044295626D-02
    w(29) = 0.27721957645934509940D-02
    w(30) = 0.29217249379178197538D-02
    w(31) = 0.30730184347025783234D-02
    w(32) = 0.32259500250878684614D-02
    w(33) = 0.33803979910869203823D-02
    w(34) = 0.35362449977167777340D-02
    w(35) = 0.36933779170256508183D-02
    w(36) = 0.38516876166398709241D-02
    w(37) = 0.40110687240750233989D-02
    w(38) = 0.41714193769840788528D-02
    w(39) = 0.43326409680929828545D-02
    w(40) = 0.44946378920320678616D-02
    w(41) = 0.46573172997568547773D-02
    w(42) = 0.48205888648512683476D-02
    w(43) = 0.49843645647655386012D-02
    w(44) = 0.51485584789781777618D-02
    w(45) = 0.53130866051870565663D-02
    w(46) = 0.54778666939189508240D-02
    w(47) = 0.56428181013844441585D-02
    w(48) = 0.58078616599775673635D-02
    w(49) = 0.59729195655081658049D-02
    w(50) = 0.61379152800413850435D-02
    w(51) = 0.63027734490857587172D-02
    w(52) = 0.64674198318036867274D-02
    w(53) = 0.66317812429018878941D-02
    w(54) = 0.67957855048827733948D-02
    w(55) = 0.69593614093904229394D-02
    w(56) = 0.71224386864583871532D-02
    w(57) = 0.72849479805538070639D-02
    w(58) = 0.74468208324075910174D-02
    w(59) = 0.76079896657190565832D-02
    w(60) = 0.77683877779219912200D-02
    w(61) = 0.79279493342948491103D-02
    w(62) = 0.80866093647888599710D-02
    w(63) = 0.82443037630328680306D-02
    w(64) = 0.84009692870519326354D-02
    w(65) = 0.85565435613076896192D-02
    w(66) = 0.87109650797320868736D-02
    w(67) = 0.88641732094824942641D-02
    w(68) = 0.90161081951956431600D-02
    w(69) = 0.91667111635607884067D-02
    w(70) = 0.93159241280693950932D-02
    w(71) = 0.94636899938300652943D-02
    w(72) = 0.96099525623638830097D-02
    w(73) = 0.97546565363174114611D-02
    w(74) = 0.98977475240487497440D-02
    w(75) = 0.10039172044056840798D-01
    w(76) = 0.10178877529236079733D-01
    w(77) = 0.10316812330947621682D-01
    w(78) = 0.10452925722906011926D-01
    w(79) = 0.10587167904885197931D-01
    w(80) = 0.10719490006251933623D-01
    w(81) = 0.10849844089337314099D-01
    w(82) = 0.10978183152658912470D-01
    w(83) = 0.11104461134006926537D-01
    w(84) = 0.11228632913408049354D-01
    w(85) = 0.11350654315980596602D-01
    w(86) = 0.11470482114693874380D-01
    w(87) = 0.11588074033043952568D-01
    w(88) = 0.11703388747657003101D-01
    w(89) = 0.11816385890830235763D-01
    w(90) = 0.11927026053019270040D-01
    w(91) = 0.12035270785279562630D-01
    w(92) = 0.12141082601668299679D-01
    w(93) = 0.12244424981611985899D-01
    w(94) = 0.12345262372243838455D-01
    w(95) = 0.12443560190714035263D-01
    w(96) = 0.12539284826474884353D-01
    w(97) = 0.12632403643542078765D-01
    w(98) = 0.12722884982732382906D-01
    w(99) = 0.12810698163877361967D-01
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
    stop

  end if

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
!  Examples:
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
!    09 March 2008
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
subroutine i4mat_transpose_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! I4MAT_TRANSPOSE_PRINT_SOME prints some of the transpose of an I4MAT.
!
!  Discussion:
!
!    An I4MAT is a rectangular array of integer values.
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
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, integer ( kind = 4 ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ILO, JLO, the first row and column to print.
!
!    Input, integer IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none

  integer   ( kind = 4 ), parameter :: incx = 10
  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) n

  integer   ( kind = 4 ) a(m,n)
  character ( len = 8 )  ctemp(incx)
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) i2
  integer   ( kind = 4 ) i2hi
  integer   ( kind = 4 ) i2lo
  integer   ( kind = 4 ) ihi
  integer   ( kind = 4 ) ilo
  integer   ( kind = 4 ) inc
  integer   ( kind = 4 ) j
  integer   ( kind = 4 ) j2hi
  integer   ( kind = 4 ) j2lo
  integer   ( kind = 4 ) jhi
  integer   ( kind = 4 ) jlo
  character ( len = * )  title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  do i2lo = max ( ilo, 1 ), min ( ihi, m ), incx

    i2hi = i2lo + incx - 1
    i2hi = min ( i2hi, m )
    i2hi = min ( i2hi, ihi )

    inc = i2hi + 1 - i2lo

    write ( *, '(a)' ) ' '

    do i = i2lo, i2hi
      i2 = i + 1 - i2lo
      write ( ctemp(i2), '(i8)' ) i
    end do

    write ( *, '(''  Row '',10a8)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Col'
    write ( *, '(a)' ) ' '

    j2lo = max ( jlo, 1 )
    j2hi = min ( jhi, n )

    do j = j2lo, j2hi

      do i2 = 1, inc

        i = i2lo - 1 + i2

        write ( ctemp(i2), '(i8)' ) a(i,j)

      end do

      write ( *, '(i5,1x,10a8)' ) j, ( ctemp(i), i = 1, inc )

    end do

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
!    The arrangement described here works naturally for the
!    Fejer Type 2, Gauss-Patterson, and Newton Cotes Open rules.
!
!    It also can be used, partially, to describe
!    the growth of Gauss-Legendre rules.
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
!    For the Fejer Type 2, Gauss-Patterson, and Newton Cotes Open rules,
!    the point growth is nested.
!
!    If we have ORDER points on a particular LEVEL, the next level
!    includes all these old points, plus ORDER+1 new points, formed in the
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
!    If we use a series of Gauss-Legendre rules, then there is almost no
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
!    Gauss-Legendre rules, then we must sum the "NEW" column, and we see that
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
      order(dim) = 2**( level(dim) + 1 ) - 1
    end if

  end do

  return
end
subroutine multigrid_index1 ( dim_num, order_1d, order_nd, indx )

!*****************************************************************************80
!
!! MULTIGRID_INDEX1 returns an indexed multidimensional grid.
!
!  Discussion:
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

    indx(1:dim_num,p) = a(1:dim_num) + 1

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
function nco_abscissa ( order, i )

!*****************************************************************************80
!
!! NCO_ABSCISSA returns the I-th abscissa for the Newton Cotes Open rule.
!
!  Discussion:
!
!    Our convention is that the abscissas are numbered from left to
!    right.
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
!    Input, integer ( kind = 4 ) ORDER, the order.
!    1 <= ORDER.
!
!    Input, integer ( kind = 4 ) I, the index of the desired abscissa.
!    1 <= I <= ORDER.
!
!    Output, real ( kind = 8 ) NCO_ABSCISSA, the value of the I-th
!    abscissa in the Newton Cotes open rule of order ORDER.
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) nco_abscissa
  integer ( kind = 4 ) order
  real ( kind = 8 ) :: x_max = +1.0D+00
  real ( kind = 8 ) :: x_min = -1.0D+00

  if ( order < 1 ) then
    nco_abscissa = - huge ( nco_abscissa )
    return
  end if

  if ( i < 1 .or. order < i ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'NCO_ABSCISSA - Fatal error!'
    write ( *, '(a)' ) '  1 <= I <= ORDER is required!'
    stop
  end if

  nco_abscissa = ( real ( order - i + 1, kind = 8 ) * x_min   &
                 + real (         i,     kind = 8 ) * x_max ) &
                 / real ( order     + 1, kind = 8 )

  return
end
subroutine nco_weights ( order, w )

!*****************************************************************************80
!
!! NCO_WEIGHTS computes weights for a Newton Cotes Open rule.
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
!    Input, integer ( kind = 4 ) ORDER, the order.
!
!    Output, real ( kind = 8 ) W(ORDER), the weights.
!
  implicit none

  integer ( kind = 4 ) order

  real ( kind = 8 ) diftab(order)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) w(order)
  real ( kind = 8 ) :: x_max = +1.0D+00
  real ( kind = 8 ) :: x_min = -1.0D+00
  real ( kind = 8 ) x(order)
  real ( kind = 8 ) yvala
  real ( kind = 8 ) yvalb
!
!  Set the abscissas.
!
  do i = 1, order
    x(i) = ( real ( order + 1 - i, kind = 8 ) * x_min   &
           + real (             i, kind = 8 ) * x_max ) &
           / real ( order + 1,     kind = 8 )
  end do
!
!  For the I-th abscissa, compute the Lagrange basis polynomial L(I)(X)
!  which is 1 at X(I), and zero at the other abscissas.
!
  do i = 1, order

    diftab(1:order) = 0.0D+00
    diftab(i) = 1.0D+00

    do j = 2, order
      do k = j, order
        diftab(order+j-k) = ( diftab(order+j-k-1) - diftab(order+j-k) ) &
          / ( x(order+1-k) - x(order+j-k) )
      end do
    end do

    do j = 1, order-1
      do k = 1, order-j
        diftab(order-k) = diftab(order-k) - x(order-k-j+1) * &
          diftab(order-k+1)
      end do
    end do
!
!  Evaluate the antiderivative of L(I)(X) at the left and
!  right endpoints.
!
    yvala = diftab(order) / real ( order, kind = 8 )
    do j = order-1, 1, -1
      yvala = yvala * x_min + diftab(j) / real ( j, kind = 8 )
    end do
    yvala = yvala * x_min

    yvalb = diftab(order) / real ( order, kind = 8 )
    do j = order-1, 1, -1
      yvalb = yvalb * x_max + diftab(j) / real ( j, kind = 8 )
    end do
    yvalb = yvalb * x_max
!
!  The weight W(I) is the integral of L(I)(X) over the interval,
!  that is, the difference between the values of the antiderivative
!  at the right and left endpoints.
!
    w(i) = yvalb - yvala

  end do

  return
end
subroutine product_weights_open ( dim_num, order_1d, order_nd, rule, w_nd )

!*****************************************************************************80
!
!! PRODUCT_WEIGHTS_OPEN: weights for an open product rule.
!
!  Discussion:
!
!    This routine computes the weights for a quadrature rule which is
!    a product of 1D rules of varying order.
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
!    Input, integer ( kind = 4 ) ORDER_1D(DIM_NUM), the order of the 1D rules.
!
!    Input, integer ( kind = 4 ) ORDER_ND, the order of the product rule.
!
!    Input, integer ( kind = 4 ) RULE, the 1D quadrature rule being used.
!    2, Fejer Type 2 Rule;
!    3, Gauss-Patterson Rule,
!    4, Newton-Cotes Open Rule,
!    5, Tanh-Sinh Rule.
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

    if ( rule == 2 ) then
      call f2_weights ( order_1d(dim), w_1d )
    else if ( rule == 3 ) then
      call gp_weights ( order_1d(dim), w_1d )
    else if ( rule == 4 ) then
      call nco_weights ( order_1d(dim), w_1d )
    else if ( rule == 5 ) then
      call ts_weights ( order_1d(dim), w_1d )
    end if

    call r8vec_direct_product2 ( dim, order_1d(dim), w_1d, dim_num, &
      order_nd, w_nd )

    deallocate ( w_1d )

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
!    12 October 2007
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
subroutine r8mat_transpose_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8MAT_TRANSPOSE_PRINT_SOME prints some of an R8MAT, transposed.
!
!  Discussion:
!
!    An R8MAT is a two dimensional matrix of double precision real values.
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
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ( kind = 4 ) ILO, JLO, the first row and column to print.
!
!    Input, integer ( kind = 4 ) IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none

  integer   ( kind = 4 ), parameter :: incx = 5
  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = 14 ) ctemp(incx)
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) i2
  integer   ( kind = 4 ) i2hi
  integer   ( kind = 4 ) i2lo
  integer   ( kind = 4 ) ihi
  integer   ( kind = 4 ) ilo
  integer   ( kind = 4 ) inc
  integer   ( kind = 4 ) j
  integer   ( kind = 4 ) j2hi
  integer   ( kind = 4 ) j2lo
  integer   ( kind = 4 ) jhi
  integer   ( kind = 4 ) jlo
  character ( len = * )  title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  do i2lo = max ( ilo, 1 ), min ( ihi, m ), incx

    i2hi = i2lo + incx - 1
    i2hi = min ( i2hi, m )
    i2hi = min ( i2hi, ihi )

    inc = i2hi + 1 - i2lo

    write ( *, '(a)' ) ' '

    do i = i2lo, i2hi
      i2 = i + 1 - i2lo
      write ( ctemp(i2), '(i8,6x)' ) i
    end do

    write ( *, '(''  Row   '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Col'
    write ( *, '(a)' ) ' '

    j2lo = max ( jlo, 1 )
    j2hi = min ( jhi, n )

    do j = j2lo, j2hi

      do i2 = 1, inc
        i = i2lo - 1 + i2
        write ( ctemp(i2), '(g14.6)' ) a(i,j)
      end do

      write ( *, '(i5,1x,5a14)' ) j, ( ctemp(i), i = 1, inc )

    end do

  end do

  return
end
subroutine r8mat_write ( output_filename, m, n, table )

!*****************************************************************************80
!
!! R8MAT_WRITE writes an R8MAT file.
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
  real ( kind = 8 ) table(m,n)
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
!  For greater precision in the output file, try:
!
!                                            '(', m, 'g', 24, '.', 16, ')'
!
  if ( 0 < m .and. 0 < n ) then

    write ( string, '(a1,i8,a1,i8,a1,i8,a1)' ) '(', m, 'g', 14, '.', 6, ')'
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
subroutine r8vec_print_some ( n, a, i_lo, i_hi, title )

!*****************************************************************************80
!
!! R8VEC_PRINT_SOME prints "some" of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8 values.
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
!    Input, integer ( kind = 4 ) N, the number of entries of the vector.
!
!    Input, real ( kind = 8 ) A(N), the vector to be printed.
!
!    Input, integer ( kind = 4 ) I_LO, I_HI, the first and last indices
!    to print.  The routine expects 1 <= I_LO <= I_HI <= N.
!
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none

  integer   ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) i_hi
  integer   ( kind = 4 ) i_lo
  character ( len = * )  title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '
  do i = max ( i_lo, 1 ), min ( i_hi, n )
    write ( *, '(2x,i8,2x,g14.8)' ) i, a(i)
  end do

  return
end
subroutine s_blank_delete ( s )

!*****************************************************************************80
!
!! S_BLANK_DELETE removes blanks from a string, left justifying the remainder.
!
!  Discussion:
!
!    All TAB characters are also removed.
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
!    Input/output, character ( len = * ) S, the string to be transformed.
!
  implicit none

  character              c
  integer   ( kind = 4 ) get
  integer   ( kind = 4 ) put
  integer   ( kind = 4 ) nchar
  character ( len = * )  s
  character, parameter :: TAB = char ( 9 )

  put = 0
  nchar = len_trim ( s )

  do get = 1, nchar

    c = s(get:get)

    if ( c /= ' ' .and. c /= TAB ) then
      put = put + 1
      s(put:put) = c
    end if

  end do

  s(put+1:nchar) = ' '

  return
end
subroutine s_to_i4 ( s, ival, ierror, length )

!*****************************************************************************80
!
!! S_TO_I4 reads an I4 from a string.
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
!    Input, character ( len = * ) S, a string to be examined.
!
!    Output, integer ( kind = 4 ) IVAL, the integer value read from the string.
!    If the string is blank, then IVAL will be returned 0.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!    0, no error.
!    1, an error occurred.
!
!    Output, integer ( kind = 4 ) LENGTH, the number of characters of S
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
subroutine spgrid_open_index ( dim_num, level_max, point_num, grid_index )

!*****************************************************************************80
!
!! SPGRID_OPEN_INDEX computes open grids with 0 <= LEVEL <= LEVEL_MAX.
!
!  Discussion:
!
!    The necessary dimensions of GRID_INDEX can be determined by
!    calling SPGRID_OPEN_SIZE first.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 July 2008
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
!    Input, integer ( kind = 4 ) POINT_NUM, the total number of points in
!    the grids.
!
!    Output, integer ( kind = 4 ) GRID_INDEX(DIM_NUM,POINT_NUM), a list of
!    point indices, representing a subset of the product grid of level
!    LEVEL_MAX, representing (exactly once) each point that will show up in a
!    sparse grid of level LEVEL_MAX.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) factor
  integer ( kind = 4 ) grid_index(dim_num,point_num)
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: grid_index2
  integer ( kind = 4 ), allocatable, dimension ( : ) :: grid_level
  integer ( kind = 4 ) h
  integer ( kind = 4 ), parameter :: i4_2 = 2
  integer ( kind = 4 ) j
  integer ( kind = 4 ) level
  integer ( kind = 4 ), dimension ( dim_num ) :: level_1d
  integer ( kind = 4 ) level_max
  logical              more
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
      call multigrid_index1 ( dim_num, order_1d, order_nd, grid_index2 )
!
!  Only keep those points which first appear on this level.
!  If you keep a point, it is necessary to rescale each of its components
!  so that we save the coordinates as they apply on the final grid.
!
      do point = 1, order_nd

        if ( all ( mod ( grid_index2(1:dim_num,point), i4_2 ) == 1 ) ) then

          point_num2 = point_num2 + 1

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

  return
end
subroutine spgrid_open_weights ( dim_num, level_max, point_num, grid_index, &
  rule, grid_weight )

!*****************************************************************************80
!
!! SPGRID_OPEN_WEIGHTS gathers the weights.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 July 2008
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
!    Input, integer ( kind = 4 ) POINT_NUM, the total number of points in
!    the grids.
!
!    Input, integer ( kind = 4 ) GRID_INDEX(DIM_NUM,POINT_NUM), a list of
!    point indices, representing a subset of the product grid of level
!    LEVEL_MAX, representing (exactly once) each point that will show up in a
!    sparse grid of level LEVEL_MAX.
!
!    Input, integer ( kind = 4 ) RULE, the 1D quadrature rule being used.
!    2, Fejer Type 2 Rule;
!    3, Gauss-Patterson Rule,
!    4, Newton-Cotes Open Rule,
!
!    Output, real ( kind = 8 ) GRID_WEIGHT(POINT_NUM), the weights
!    associated with the sparse grid points.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) coeff
  integer ( kind = 4 ) dim
  integer ( kind = 4 ) grid_index(dim_num,point_num)
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: grid_index2
  real ( kind = 8 ) grid_weight(point_num)
  real ( kind = 8 ), allocatable, dimension ( : ) :: grid_weight2
  integer ( kind = 4 ) h
  integer ( kind = 4 ) i4_choose
  integer ( kind = 4 ) level
  integer ( kind = 4 ) level_1d(dim_num)
  integer ( kind = 4 ) level_max
  integer ( kind = 4 ) level_min
  integer ( kind = 4 ) match
  logical              more
  integer ( kind = 4 ) order_nd
  integer ( kind = 4 ) order_1d(dim_num)
  integer ( kind = 4 ) point
  integer ( kind = 4 ) point2
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
      call multigrid_index1 ( dim_num, order_1d, order_nd, grid_index2 )
!
!  Compute the weights for this grid.
!
      call product_weights_open ( dim_num, order_1d, order_nd, rule, &
        grid_weight2 )
!
!  Adjust the grid indices to reflect LEVEL_MAX.
!
      call multigrid_scale_open ( dim_num, order_nd, level_max, level_1d, &
        grid_index2 )
!
!  Now determine the coefficient.
!
      coeff = (-1)**( level_max - level ) &
        * i4_choose ( dim_num - 1, level_max - level )

      do point2 = 1, order_nd

        match = -1

        do point = 1, point_num

          if ( all ( &
            grid_index2(1:dim_num,point2) == grid_index(1:dim_num,point) &
          ) ) then
            grid_weight(point) = grid_weight(point) &
              + real ( coeff, kind = 8 ) * grid_weight2(point2)
            match = point
            exit
          end if

        end do

        if ( match == -1 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'SPGRID_OPEN_WEIGHTS - Fatal error!'
          write ( *, '(a)' ) '  Could not match grid index.'
          write ( *, '(a,i8)' ) '  Point index = ', point2
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) ' '
          write ( *, '(a,i8)' ) '  LEVEL = ', level
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  LEVEL_1D:'
          write ( *, '(5i8)' ) level_1d(1:dim_num)
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  ORDER_1D:'
          write ( *, '(5i8)' ) order_1d(1:dim_num)
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  GRID_INDEX2:'
          write ( *, '(5i8)' ) grid_index2(1:dim_num,point2)
          stop
        end if

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
!    09 March 2008
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

  character ( len = 8 )  ampm
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
function ts_abscissa ( order, i )

!*****************************************************************************80
!
!! TS_ABSCISSA returns the I-th abscissa for the tanh-sinh rule.
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
!    01 February 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ORDER, the order of the rule.
!
!    Input, integer ( kind = 4 ) I, the index of the desired abscissa.
!    1 <= I <= ORDER.
!
!    Output, real ( kind = 8 ) TS_ABSCISSA, the value of the I-th abscissa
!    in the rule of order ORDER.
!
  implicit none

  real ( kind = 8 ) ct
  real ( kind = 8 ) ct2
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) order
  real ( kind = 8 ) :: pi = 3.141592653589793D+00
  real ( kind = 8 ) r8_huge
  real ( kind = 8 ) st
  real ( kind = 8 ) t
  real ( kind = 8 ) ts_abscissa
  real ( kind = 8 ) value

  if ( order < 1 ) then
    value = - r8_huge ( )
  else if ( i < 1 .or. order < i ) then
    value = - r8_huge ( )
  else if ( order == 1 ) then
    value = 0.0D+00
  else if ( 2 * i - order - 1 == 0 ) then
    value = 0.0D+00
  else

    h = 4.0D+00 / real ( order + 1, kind = 8 )

    t = real ( 2 * i - order - 1, kind = 8 ) * h / 2.0D+00

    ct = cosh ( t )
    st = sinh ( t )
    ct2 = cosh ( 0.5D+00 * pi * st )

    value = tanh ( 0.5D+00 * pi * st )

  end if

  ts_abscissa = value

  return
end
subroutine ts_weights ( order, w )

!*****************************************************************************80
!
!! TS_WEIGHTS computes weights for a tanh-sinh rule.
!
!  Discussion:
!
!    In the 1D case, a sequence of rules is used of increasing order.
!    For low order, the weights do not sum to 2, but with increasing
!    order, the sum quickly converges to 2.
!
!    However, for sparse grid applications, the lowest order rules are
!    involved in every grid, so it seems it might be useful to force
!    the weights to sum to 2 immediately.  This addresses only one very
!    obvious defect of the lower order rules.  I am not sure what to do
!    about the fact the none of the rules have a definable precision,
!    and the family of rules has not precision but asymptotic accuracy.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 February 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ORDER, the order of the rule.
!
!    Output, real ( kind = 8 ) W(ORDER), the weights of the rule.
!
  integer ( kind = 4 ) order

  real ( kind = 8 ) ct
  real ( kind = 8 ) ct2
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  real ( kind = 8 ) :: pi = 3.141592653589793D+00
  real ( kind = 8 ) st
  real ( kind = 8 ) t
  real ( kind = 8 ) w(order)
  real ( kind = 8 ) w_sum

  if ( order < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TS_WEIGHTS - Fatal error!'
    write ( *, '(a)' ) '  ORDER < 1.'
    stop
  end if

  h = 4.0D+00 / real ( order + 1, kind = 8 )

  do i = 1, order

    t = real ( 2 * i - order - 1, kind = 8 ) * h / 2.0D+00

    ct = cosh ( t )
    st = sinh ( t )
    ct2 = cosh ( 0.5D+00 * pi * st );

    w(i) = 0.5D+00 * pi * h * ct / ct2 / ct2

  end do
!
!  Normalize the weights so that they sum to 2.0.
!
  w_sum = sum ( w(1:order) )

  w(1:order) = 2.0D+00 * w(1:order) / w_sum

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
  logical              more

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
