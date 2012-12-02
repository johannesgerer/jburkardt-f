program main

!*****************************************************************************80
!
!! MAIN is the main program for SPARSE_GRID_CC_DATASET.
!
!  Discussion:
!
!    This program computes a sparse grid quadrature rule based on a 1D
!    Clenshaw-Curtis rule and writes it to a file..
!
!    The user specifies:
!    * the spatial dimension of the quadrature region,
!    * the level that defines the Smolyak grid.
!
!  Licensing:
!
!    This software is released under the GNU LGPL license.
!
!  Modified:
!
!    30 June 2009
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

  integer ( kind = 4 ) arg_num
  integer ( kind = 4 ) dim
  integer ( kind = 4 ) dim_num
  character ( len = 3 ) dim_num_string
  integer ( kind = 4 ) iarg
  integer ( kind = 4 ) iargc
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) last
  integer ( kind = 4 ) level_max
  character ( len = 2 ) level_max_string
  integer ( kind = 4 ) level_min
  integer ( kind = 4 ) point
  integer ( kind = 4 ) point_num
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: r
  character ( len = 255 ) r_filename
  character ( len = 255 ) string
  real ( kind = 8 ), allocatable, dimension ( : ) :: w
  character ( len = 255 ) w_filename
  real ( kind = 8 ) weight_sum
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: x
  character ( len = 255 ) x_filename

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPARSE_GRID_CC_DATASET'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Compute the abscissas and weights of a quadrature rule'
  write ( *, '(a)' ) '  associated with a sparse grid derived from a Smolyak'
  write ( *, '(a)' ) '  construction based on 1D Clenshaw-Curtis rules.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Inputs to the program include:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    DIM_NUM, the spatial dimension.'
  write ( *, '(a)' ) '    (typically in the range of 2 to 10)'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    LEVEL_MAX, the "level" of the sparse grid.'
  write ( *, '(a)' ) '    (typically in the range of 0, 1, 2, 3, ...'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Output from the program includes:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    * A printed table of the abscissas and weights.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    * A set of 3 files that define the quadrature rule.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    (1) "cc_d?_level?_r.txt", the ranges;'
  write ( *, '(a)' ) '    (2) "cc_d?_level?_w.txt", the weights;'
  write ( *, '(a)' ) '    (3) "cc_d?_level?_x.txt", the abscissas.'
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
    write ( *, '(a)' ) '  Enter the value of DIM_NUM (1 or greater)'
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
    write ( *, '(a)' ) '  Enter the value of LEVEL_MAX (0 or greater).'
    read ( *, * ) level_max
  end if

  level_min = max ( 0, level_max + 1 - dim_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  LEVEL_MIN is = ', level_min
  write ( *, '(a,i8)' ) '  LEVEL_MAX is = ', level_max
!
!  How many distinct points will there be?
!
  call sparse_grid_cfn_size ( dim_num, level_max, point_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The number of distinct abscissas in the '
  write ( *, '(a)' ) '  quadrature rule is determined from the spatial'
  write ( *, '(a)' ) '  dimension DIM_NUM and the level LEVEL_MAX.'
  write ( *, '(a,i8)' ) &
    '  For the given input, this value will be = ', point_num
!
!  Allocate memory.
!
  allocate ( r(1:dim_num,2) )
  allocate ( w(1:point_num) )
  allocate ( x(1:dim_num,1:point_num) )
!
!  Compute the weights and points.
!
  r(1:dim_num,1) = -1.0D+00
  r(1:dim_num,2) = +1.0D+00

  call sparse_grid_cc ( dim_num, level_max, point_num, w, x )

  call r8mat_transpose_print_some ( dim_num, point_num, x, &
    1, 1, dim_num, 10, '  First 10 grid points:' )

  call r8vec_print_some ( point_num, w, 1, 10, &
    '  First 10 weights:' )

  weight_sum = sum ( w(1:point_num) )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g24.16)' ) '  Weights sum to   ', weight_sum
  write ( *, '(a,g24.16)' ) '  Correct value is ', 2.0D+00**dim_num
!
!  Construct appropriate file names.
!
  if ( level_max < 10 ) then
    write ( level_max_string, '(i1,a)' ) level_max, ' '
  else
    write ( level_max_string, '(i2)' ) level_max
  end if

  if ( dim_num < 10 ) then
    write ( dim_num_string, '(i1,a)' ) dim_num, ' '
  else if ( dim_num < 100 ) then
    write ( dim_num_string, '(i2)' ) dim_num
  else
    write ( dim_num_string, '(i3)' ) dim_num
  end if

  r_filename = 'cc_d' // &
               trim ( dim_num_string )   // '_level' // &
               trim ( level_max_string ) // '_r.txt'

  w_filename = 'cc_d' // &
               trim ( dim_num_string )   // '_level' // &
               trim ( level_max_string ) // '_w.txt'

  x_filename = 'cc_d' // &
               trim ( dim_num_string )   // '_level' // &
               trim ( level_max_string ) // '_x.txt'
!
!  Write the rule to files.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Creating R file = "' // trim ( r_filename ) // '".'

  call r8mat_write ( r_filename, dim_num, 2, r )

  write ( *, '(a)' ) '  Creating W file = "' // trim ( w_filename ) // '".'

  call r8mat_write ( w_filename, 1, point_num, w )

  write ( *, '(a)' ) '  Creating X file = "' // trim ( x_filename ) // '".'

  call r8mat_write ( x_filename, dim_num, point_num, x )
!
!  Free memory.
!
  deallocate ( r )
  deallocate ( w )
  deallocate ( x )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPARSE_GRID_CC_DATASET:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
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
!    08 November 2007
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
  integer ( kind = 4 ) t
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

  else if ( 2 * i - 1 == order ) then

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
!    08 November 2007
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
function choose ( n, k )

!*****************************************************************************80
!
!! CHOOSE computes the binomial coefficient C(N,K).
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
!    08 November 2007
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
!    Output, integer ( kind = 4 ) CHOOSE, the number of combinations of N
!    things taken K at a time.
!
  implicit none

  integer ( kind = 4 ) choose
  integer ( kind = 4 ) i
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

  choose = value

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
!    variables H and T.  I have decided (based on an wasting an
!    entire morning trying to track down a problem) that it is safer
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
!    Input/output, integer ( kind = 4 ) H, T, two internal parameters needed for the
!    computation.  The user should allocate space for these in the calling
!    program, include them in the calling sequence, but never alter them!
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
!    This software is released under the GNU LGPL license.
!
!  Modified:
!
!    12 November 2007
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
!    This software is released under the GNU LGPL license.
!
!  Modified:
!
!    12 November 2007
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
function index_to_level_closed ( dim_num, t, order, level_max )

!*****************************************************************************80
!
!! INDEX_TO_LEVEL_CLOSED determines the level of a point given its index.
!
!  Licensing:
!
!    This software is released under the GNU LGPL license.
!
!  Modified:
!
!    12 November 2007
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
!    The arrangement described here works naturally for the Clenshaw Curtis
!    and Newton Cotes closed rules.
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
!    For the Clenshaw Curtis and Newton Cotes Closed rules, the point growth
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
!    This software is released under the GNU LGPL license.
!
!  Modified:
!
!    12 November 2007
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
      order(dim) = ( 2**level(dim) ) + 1
    end if

  end do

  return
end
subroutine levels_closed_index ( dim_num, level_max, point_num, grid_index )

!*****************************************************************************80
!
!! LEVELS_CLOSED_INDEX computes closed grids with 0 <= LEVEL <= LEVEL_MAX.
!
!  Discussion:
!
!    The necessary dimensions of GRID_INDEX can be determined by
!    calling LEVELS_CLOSED_INDEX_SIZE first.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 November 2007
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
!    Input, integer ( kind = 4 ) POINT_NUM, the total number of points
!    in the grids.
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
  integer ( kind = 4 ) j
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
      call multigrid_index0 ( dim_num, order_1d, order_nd, grid_index2 )
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

  return
end
subroutine multigrid_index0 ( dim_num, order_1d, order_nd, indx )

!*****************************************************************************80
!
!! MULTIGRID_INDEX0 returns an indexed multidimensional grid.
!
!  Discussion:
!
!    For dimension DIM, the second index of INDX may vary from
!    0 to ORDER_1D(DIM)-1.
!
!  Licensing:
!
!    This software is released under the GNU LGPL license.
!
!  Modified:
!
!    12 November 2007
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
  integer ( kind = 4 ) change
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
!    This software is released under the GNU LGPL license.
!
!  Modified:
!
!    12 November 2007
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
subroutine product_weights_cc ( dim_num, order_1d, order_nd, w_nd )

!*****************************************************************************80
!
!! PRODUCT_WEIGHTS_CC: Clenshaw Curtis product rule weights.
!
!  Discussion:
!
!    This routine computes the weights for a quadrature rule which is
!    a product of 1D closed rules of varying order.
!
!  Licensing:
!
!    This software is released under the GNU LGPL license.
!
!  Modified:
!
!    12 November 2007
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
!    Output, real ( kind = 8 ) W_ND(DIM_NUM,ORDER_ND), the product rule weights.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) order_nd

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) order_1d(dim_num)
  real ( kind = 8 ), allocatable, dimension ( : ) :: w_1d
  real ( kind = 8 ) w_nd(order_nd)

  w_nd(1:order_nd) = 1.0D+00

  do dim = 1, dim_num

    allocate ( w_1d(1:order_1d(dim)) )

    call cc_weights ( order_1d(dim), w_1d )

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
!    This software is released under the GNU LGPL license.
!
!  Modified:
!
!    12 November 2007
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
!    This software is released under the GNU LGPL license.
!
!  Modified:
!
!    12 November 2007
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
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = 14 ) ctemp(incx)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

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

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) j
  character ( len = * ) output_filename
  integer ( kind = 4 ) output_status
  integer ( kind = 4 ) output_unit
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
!    This software is released under the GNU LGPL license.
!
!  Modified:
!
!    12 November 2007
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
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_hi
  integer ( kind = 4 ) i_lo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '
  do i = max ( i_lo, 1 ), min ( i_hi, n )
    write ( *, '(2x,i8,2x,g14.8)' ) i, a(i)
  end do

  return
end
subroutine s_to_i4 ( s, value, ierror, length )

!*****************************************************************************80
!
!! S_TO_I4 reads an integer value from a string.
!
!  Discussion:
!
!    Instead of ICHAR, we now use the IACHAR function, which
!    guarantees the ASCII collating sequence.
!
!  Licensing:
!
!    This software is released under the GNU LGPL license.
!
!  Modified:
!
!    12 November 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, a string to be examined.
!
!    Output, integer ( kind = 4 ) VALUE, the integer value read from the string.
!    If the string is blank, then VALUE will be returned 0.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!    0, no error.
!    1, an error occurred.
!
!    Output, integer ( kind = 4 ) LENGTH, the number of characters
!    of S used to make the integer.
!
  implicit none

  character c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) length
  character ( len = * ) s
  integer ( kind = 4 ) state
  integer ( kind = 4 ) value

  value = 0
  ierror = 0
  length = 0

  state = 0
  isgn = 1

  do i = 1, len_trim ( s )

    c = s(i:i)
!
!  STATE = 0, haven't read anything.
!
    if ( state == 0 ) then

      if ( c == ' ' ) then

      else if ( c == '-' ) then
        state = 1
        isgn = -1
      else if ( c == '+' ) then
        state = 1
        isgn = +1
      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        state = 2
        value = iachar ( c ) - iachar ( '0' )
      else
        ierror = 1
        return
      end if
!
!  STATE = 1, have read the sign, expecting digits or spaces.
!
    else if ( state == 1 ) then

      if ( c == ' ' ) then

      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        state = 2
        value = iachar ( c ) - iachar ( '0' )
      else
        ierror = 1
        return
      end if
!
!  STATE = 2, have read at least one digit, expecting more.
!
    else if ( state == 2 ) then

      if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then

        value = 10 * value + iachar ( c ) - iachar ( '0' )

      else

        value = isgn * value
        ierror = 0
        length = i - 1
        return

      end if

    end if

  end do
!
!  If we read all the characters in the string, see if we're OK.
!
  if ( state == 2 ) then

    value = isgn * value
    ierror = 0
    length = len_trim ( s )

  else

    value = 0
    ierror = 1
    length = 0

  end if

  return
end
subroutine sparse_grid_cc ( dim_num, level_max, point_num, grid_weight, &
  grid_point )

!****************************************************************************80
!
!! SPARSE_GRID_CC computes a sparse grid of Clenshaw Curtis points.
!
!  Discussion:
!
!    This program computes a quadrature rule and writes it to a file.
!
!    The quadrature rule is associated with a sparse grid derived from
!    a Smolyak construction using a closed 1D quadrature rule.
!
!    The user specifies:
!    * the spatial dimension of the quadrature region,
!    * the level that defines the Smolyak grid.
!    * the closed 1D quadrature rule (Clenshaw-Curtis or Newton-Cotes Closed).
!
!  Licensing:
!
!    This software is released under the GNU LGPL license.
!
!  Modified:
!
!    12 November 2007
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
!    Input, integer ( kind = 4 ) LEVEL_MAX, controls the size of the final
!    sparse grid.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points in the grid,
!    as determined by SPARSE_GRID_CC_SIZE.
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
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: grid_index
  real ( kind = 8 ) grid_point(dim_num,point_num)
  real ( kind = 8 ) grid_weight(point_num)
  integer ( kind = 4 ) level_max
  integer ( kind = 4 ) order_max
  integer ( kind = 4 ) point
!
!  Determine the index vector, relative to the full product grid,
!  that identifies the points in the sparse grid.
!
  allocate ( grid_index(dim_num,point_num) )

  call sparse_grid_cc_index ( dim_num, level_max, point_num, grid_index )
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
      grid_point(dim,point) = &
        cc_abscissa ( order_max, grid_index(dim,point) + 1 )
    end do
  end do
!
!  Gather the weights.
!
  call sparse_grid_cc_weights ( dim_num, level_max, point_num, grid_index, &
     grid_weight )

  deallocate ( grid_index )

  return
end
subroutine sparse_grid_cc_index ( dim_num, level_max, point_num, grid_index )

!*****************************************************************************80
!
!! SPARSE_GRID_CC_INDEX indexes the points forming a sparse grid.
!
!  Discussion:
!
!    The points forming the sparse grid are guaranteed to be a subset
!    of a certain product grid.  The product grid is formed by DIM_NUM
!    copies of a 1D rule of fixed order.  The orders of the 1D rule,
!    (called ORDER_1D) and the order of the product grid, (called ORDER)
!    are determined from the value LEVEL_MAX.
!
!    Thus, any point in the product grid can be identified by its grid index,
!    a set of DIM_NUM indices, each between 1 and ORDER_1D.
!
!    This routine creates the GRID_INDEX array, listing (uniquely) the
!    points of the sparse grid.
!
!    An assumption has been made that the 1D rule is closed (includes
!    the interval endpoints) and nested (points that are part of a rule
!    of a given level will be part of every rule of higher level).
!
!    The necessary dimensions of GRID_INDEX can be determined by
!    calling SPARSE_GRID_CC_SIZE first.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 November 2007
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
  integer ( kind = 4 ) j
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
      call multigrid_index0 ( dim_num, order_1d, order_nd, grid_index2 )
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

  return
end
subroutine sparse_grid_cfn_size ( dim_num, level_max, point_num )

!*****************************************************************************80
!
!! SPARSE_GRID_CC_SIZE sizes a sparse grid using Closed Fully Nested rules.
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
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
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
  do i = 2, level_max
    j = j * 2
    new_1d(i) = j
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
subroutine sparse_grid_cc_weights ( dim_num, level_max, point_num, grid_index, &
  grid_weight )

!*****************************************************************************80
!
!! SPARSE_GRID_CC_WEIGHTS gathers the weights.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 November 2007
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
!    Output, real ( kind = 8 ) GRID_WEIGHT(POINT_NUM), the weights
!    associated with the sparse grid points.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) coeff
  integer ( kind = 4 ) choose
  integer ( kind = 4 ) dim
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
      call multigrid_index0 ( dim_num, order_1d, order_nd, grid_index2 )
!
!  Compute the weights for this grid.
!
      call product_weights_cc ( dim_num, order_1d, order_nd, grid_weight2 )
!
!  Adjust the grid indices to reflect LEVEL_MAX.
!
      call multigrid_scale_closed ( dim_num, order_nd, level_max, level_1d, &
        grid_index2 )
!
!  Now determine the coefficient.
!
      coeff = (-1)**( level_max - level ) &
        * choose ( dim_num - 1, level_max - level )

      do point2 = 1, order_nd

        do point = 1, point_num

          if ( all ( &
            grid_index2(1:dim_num,point2) == grid_index(1:dim_num,point) &
          ) ) then
            grid_weight(point) = grid_weight(point) &
              + real ( coeff, kind = 8 ) * grid_weight2(point2)
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
!    This software is released under the GNU LGPL license.
!
!  Modified:
!
!    12 November 2007
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
!    This software is released under the GNU LGPL license.
!
!  Modified:
!
!    12 November 2007
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
