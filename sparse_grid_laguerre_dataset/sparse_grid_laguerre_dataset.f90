program main

!*****************************************************************************80
!
!! MAIN is the main program for SPARSE_GRID_LAGUERRE_DATASET.
!
!  Discussion:
!
!    This program computes a sparse grid quadrature rule based on a 1D
!    Gauss-Laguerre rule and writes it to a file..
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
!    06 October 2007
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
  character ( len = 2 ) dim_num_string
  integer   ( kind = 4 ) iarg
  integer   ( kind = 4 ) iargc
  integer   ( kind = 4 ) ierror
  integer   ( kind = 4 ) last
  integer   ( kind = 4 ) level_max
  character ( len = 2 ) level_max_string
  integer   ( kind = 4 ) level_min
  integer   ( kind = 4 ) point
  integer   ( kind = 4 ) point_num
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: r
  character ( len = 80 ) r_filename
  real ( kind = 8 ) r8_huge
  character ( len = 80 ) string
  real ( kind = 8 ), allocatable, dimension ( : ) :: w
  character ( len = 80 ) w_filename
  real ( kind = 8 ) weight_sum
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: x
  character ( len = 80 ) x_filename

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPARSE_GRID_LAGUERRE_DATASET'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Compute the abscissas and weights of a quadrature rule'
  write ( *, '(a)' ) '  associated with a sparse grid derived from a Smolyak'
  write ( *, '(a)' ) '  construction based on a 1D Gauss-Laguerre rule.'
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
  write ( *, '(a)' ) '    (1) "lag_d?_level?_x.txt", the abscissas;'
  write ( *, '(a)' ) '    (2) "lag_d?_level?_w.txt", the weights;'
  write ( *, '(a)' ) '    (3) "lag_d?_level?_r.txt", the ranges.'
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
  call sparse_grid_laguerre_size ( dim_num, level_max, point_num )

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
  r(1:dim_num,1) =   0.0D+00
  r(1:dim_num,2) = + r8_huge ( )

  call sparse_grid_laguerre ( dim_num, level_max, point_num, w, x )

  call r8mat_transpose_print_some ( dim_num, point_num, x, &
    1, 1, dim_num, 10, '  First 10 grid points:' )

  call r8vec_print_some ( point_num, w, 1, 10, &
    '  First 10 grid weights:' )

  weight_sum = sum ( w(1:point_num) )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g24.16)' ) '  Weights sum to   ', weight_sum
  write ( *, '(a,g24.16)' ) '  Correct value is ', 1.0D+00
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
  else
    write ( dim_num_string, '(i2)' ) dim_num
  end if

  r_filename = 'lag_d' // &
               trim ( dim_num_string )   // '_level' // &
               trim ( level_max_string ) // '_r.txt'

  w_filename = 'lag_d' // &
               trim ( dim_num_string )   // '_level' // &
               trim ( level_max_string ) // '_w.txt'

  x_filename = 'lag_d' // &
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
  write ( *, '(a)' ) 'SPARSE_GRID_LAGUERRE_DATASET'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
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
!    21 May 2007
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
function i4_huge ( )

!*****************************************************************************80
!
!! I4_HUGE returns a "huge" I4.
!
!  Discussion:
!
!    On an IEEE 32 bit machine, I4_HUGE should be 2**31 - 1, and its
!    bit pattern should be
!
!     01111111111111111111111111111111
!
!    In this case, its numerical value is 2147483647.
!
!    Using the Dec/Compaq/HP Alpha FORTRAN compiler FORT, I could
!    use I4_HUGE() and HUGE interchangeably.
!
!    However, when using the G95, the values returned by HUGE were
!    not equal to 2147483647, apparently, and were causing severe
!    and obscure errors in my random number generator, which needs to
!    add I4_HUGE to the seed whenever the seed is negative.  So I
!    am backing away from invoking HUGE, whereas I4_HUGE is under
!    my control.
!
!    Explanation: because under G95 the default integer type is 64 bits!
!    So HUGE ( 1 ) = a very very huge integer indeed, whereas
!    I4_HUGE ( ) = the same old 32 bit big value.
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) I4_HUGE, a "huge" I4.
!
  implicit none

  integer ( kind = 4 ) i4
  integer ( kind = 4 ) i4_huge

  i4_huge = 2147483647

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
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 January 2003
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
  integer ( kind = 4 ) i4_huge

  if ( i == 0 ) then

    i4_log_2 = - i4_huge ( )

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
subroutine laguerre_abscissa ( dim_num, point_num, grid_index, grid_base, &
  grid_point )

!*****************************************************************************80
!
!! LAGUERRE_ABSCISSA sets abscissas for multidimensional Gauss-Laguerre quadrature.
!
!  Discussion:
!
!    The "nesting" as it occurs for Gauss-Laguerre sparse grids simply
!    involves the use of a specified set of permissible orders for the
!    rule.
!
!    The X array lists the (complete) Gauss-Laguerre abscissas for rules
!    of order 1, 3, 7, 15, 31, 63, and 127 in order.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 October 2007
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
!    point and dimension, the index of the Gauss-Laguerre abscissa.
!
!    Input, integer ( kind = 4 ) GRID_BASE(DIM_NUM), the "base" of the
!    Gauss-Laguerre rule being used in each dimension.
!
!    Output, real ( kind = 8 ) GRID_POINT(DIM_NUM), the grid points of
!    Gauss-Laguerre abscissas.
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
  integer ( kind = 4 ), dimension ( 0:7 ) :: skip = (/ 0, 1, 4, 11, 26, 57, 120, 247 /)
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
    write ( *, '(a)' ) 'LAGUERRE_ABSCISSA - Fatal error!'
    write ( *, '(a)' ) '  Some base values are less than 1.'
    stop
  end if

  if ( any ( 127 < grid_base(1:dim_num) ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LAGUERRE_ABSCISSA - Fatal error!'
    write ( *, '(a)' ) '  Some base values are greater than 127.'
    stop
  end if

  do point = 1, point_num
    do dim = 1, dim_num

      level = i4_log_2 ( grid_base(dim) + 1 ) - 1

      pointer = skip(level) + grid_index(dim,point)

      if ( pointer < 1 .or. 247 < pointer ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'LAGUERRE_ABSCISSA - Fatal error!'
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
subroutine laguerre_integral_nd ( dim_num, expon, value )

!*****************************************************************************80
!
!! LAGUERRE_INTEGRAL_ND evaluates a multidimensional Laguerre polynomial integral.
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
!    05 October 2007
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
  integer ( kind = 4 ) expon(dim_num)
  real ( kind = 8 ) r8_factorial
  real ( kind = 8 ) value

  value = 1.0D+00
  do dim = 1, dim_num
    value = value * r8_factorial ( expon(dim) )
  end do

  return
end
subroutine laguerre_weights ( order, weight )

!*****************************************************************************80
!
!! LAGUERRE_WEIGHTS returns weights for certain Gauss-Laguerre quadrature rules.
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
!    02 October 2007
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
    write ( *, '(a)' ) 'LAGUERRE_WEIGHTS - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of ORDER = ', order
    write ( *, '(a)' ) '  Legal values are 1, 3, 7, 15, 31, 63 and 127.'
    stop

  end if

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
!    Fejer Type 2, Gauss-Patterson, Newton Cotes Open, and
!    Newton Cotes Half Open rules.  It also can be used, partially, to describe
!    the growth of Gauss-Legendre and Gauss-Hermite rules.
!
!    The idea is that we start with LEVEL = 0, ORDER = 1 indicating the single
!    point at the center, and for all values afterwards, we use the
!    relationship
!
!      ORDER = 2**(LEVEL+1) - 1.
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
!    For the Fejer Type 1, Fejer Type 2, Gauss-Patterson, Newton Cotes Open,
!    and Newton Cotes Open Half rules, the point growth is
!    nested.  If we have ORDER points on a particular LEVEL, the next level
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
!    If we use a series of Gauss-Legendre or Gauss-Hermite rules, then there is almost no
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
!    Gauss-Hermite or Gauss-Legendre rules, then we must sum the "NEW" column,
!    and we see that we get roughly twice as many points as for the truly nested rules.
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
!    14 April 2007
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
!    Input, integer ( kind = 4 ) LEVEL(DIM_NUM), the nesting levels of the 1D rules.
!
!    Output, integer ( kind = 4 ) ORDER(DIM_NUM), the order (number of points) of the
!    1D rules.
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
subroutine monomial_quadrature ( dim_num, expon, point_num, weight, x, &
  quad_error )

!*****************************************************************************80
!
!! MONOMIAL_QUADRATURE applies a quadrature rule to a monomial.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 August 2007
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
!    Output, real ( kind = 8 ) QUAD_ERROR, the quadrature error.
!
  implicit none

  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) exact
  integer ( kind = 4 ) expon(dim_num)
  integer ( kind = 4 ) point_num
  real ( kind = 8 ) quad
  real ( kind = 8 ) quad_error
  real ( kind = 8 ) value(point_num)
  real ( kind = 8 ) weight(point_num)
  real ( kind = 8 ) x(dim_num,point_num)
!
!  Get the exact value of the integral of the unscaled monomial.
!
  call laguerre_integral_nd ( dim_num, expon, exact )
!
!  Evaluate the monomial at the quadrature points.
!
  call monomial_value ( dim_num, point_num, x, expon, value )
!
!  Compute the weighted sum.
!
  quad = dot_product ( weight, value )
!
!  If exact value is nonzero, use it to scale the data.
!
  if ( exact == 0.0D+00 ) then
    quad_error = abs ( quad )
  else
    quad_error = abs ( ( quad - exact ) / exact )
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
!    if the combination 0^0 is encountered, it should be treated
!    as 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 October 2007
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
  integer ( kind = 4 ) point
  real ( kind = 8 ) value(point_num)
  real ( kind = 8 ) x(dim_num,point_num)

  value(1:point_num) = 1.0D+00

  do dim = 1, dim_num
    do point = 1, point_num
      if ( x(dim,point) /= 0.0D+00 ) then
        value(point) = value(point) * x(dim,point)**expon(dim)
      else if ( expon(dim) == 0 ) then
        value(point) = value(point)
      else
        value(point) = 0.0D+00
      end if
    end do
  end do

  return
end
subroutine multigrid_index_one ( dim_num, order_1d, order_nd, indx )

!*****************************************************************************80
!
!! MULTIGRID_INDEX_ONE returns an indexed multidimensional grid.
!
!  Discussion:
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
!    04 October 2007
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
subroutine product_weight_laguerre ( dim_num, order_1d, order_nd, w_nd )

!*****************************************************************************80
!
!! PRODUCT_WEIGHT_LAGUERRE: weights for a product Gauss-Laguerre rule.
!
!  Discussion:
!
!    This routine computes the weights for a quadrature rule which is
!    a product of 1D Gauss-Laguerre rules of varying order.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 October 2007
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
!    Output, real ( kind = 8 ) W_ND(ORDER_ND), the product rule weights.
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

    call laguerre_weights ( order_1d(dim), w_1d )

    call r8vec_direct_product2 ( dim, order_1d(dim), w_1d, dim_num, &
      order_nd, w_nd )

    deallocate ( w_1d )

  end do

  return
end
function r8_factorial ( n )

!*****************************************************************************80
!
!! R8_FACTORIAL computes the factorial of N, also denoted "N!".
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
!    16 January 1999
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
function r8_huge ( )

!*****************************************************************************80
!
!! R8_HUGE returns the largest legal R8.
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
!  Modified:
!
!    06 October 2007
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
!  Modified:
!
!    14 June 2004
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
!  Modified:
!
!    16 October 2006
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

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_hi
  integer ( kind = 4 ) i_lo
  character ( len = * ) title

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
!  Comment:
!
!    All TAB characters are also removed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 July 1998
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

  character c
  integer ( kind = 4 ) get
  integer ( kind = 4 ) put
  integer ( kind = 4 ) nchar
  character ( len = * ) s
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
!  Modified:
!
!    13 January 2006
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
subroutine sparse_grid_laguerre ( dim_num, level_max, point_num, grid_weight, &
  grid_point )

!****************************************************************************80
!
!! SPARSE_GRID_LAGUERRE computes a sparse grid of Gauss-Laguerre points.
!
!  Discussion:
!
!    The quadrature rule is associated with a sparse grid derived from
!    a Smolyak construction using a 1D Gauss-Laguerre quadrature rule.
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
!    05 July 2008
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
!    Input, integer ( kind = 4 ) LEVEL_MAX, controls the size of the
!    sparse grid.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points in the grid,
!    as determined by SPARSE_GRID_LAGUERRE_SIZE.
!
!    Output, real ( kind = 8 ) GRID_WEIGHT(POINT_NUM), the weights.
!
!    Output, real ( kind = 8 ) GRID_POINT(DIM_NUM,POINT_NUM), the points.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) choose
  integer ( kind = 4 ) coeff
  integer ( kind = 4 ) dim
  integer ( kind = 4 ) factor
  integer ( kind = 4 ), dimension ( dim_num ) :: grid_base2
  integer ( kind = 4 ) grid_index(dim_num,point_num)
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: grid_index2
  integer ( kind = 4 ), allocatable, dimension ( : ) :: grid_level
  real ( kind = 8 ) grid_point(dim_num,point_num)
  real ( kind = 8 ) grid_point_temp(dim_num)
  real ( kind = 8 ) grid_weight(point_num)
  real ( kind = 8 ), allocatable, dimension ( : ) :: grid_weight2
  integer ( kind = 4 ) h
  integer ( kind = 4 ) j
  integer ( kind = 4 ) level
  integer ( kind = 4 ), dimension ( dim_num ) :: level_1d
  integer ( kind = 4 ) level_max
  integer ( kind = 4 ) level_min
  logical more
  integer ( kind = 4 ), dimension ( dim_num ) :: order_1d
  integer ( kind = 4 ) order_nd
  integer ( kind = 4 ) point
  integer ( kind = 4 ) point2
  integer ( kind = 4 ) point3
  integer ( kind = 4 ) point_num2
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
      call product_weight_laguerre ( dim_num, order_1d, order_nd, grid_weight2 )
!
!  Now determine the coefficient of the weight.
!
      coeff = (-1)**( level_max - level ) &
        * choose ( dim_num - 1, level_max - level )
!
!  The inner (hidden) loop generates all points corresponding to given grid.
!  The grid indices will be between -M to +M, where 2*M + 1 = ORDER_1D(DIM).
!
      call multigrid_index_one ( dim_num, order_1d, order_nd, grid_index2 )

      do point = 1, order_nd

        point_num2 = point_num2 + 1

        call laguerre_abscissa ( dim_num, 1, grid_index2(1:dim_num,point), &
          grid_base2(1:dim_num), grid_point(1:dim_num,point_num2) )

        grid_weight(point_num2) = real ( coeff, kind = 8 ) * grid_weight2(point)

      end do

      deallocate ( grid_index2 )
      deallocate ( grid_level )
      deallocate ( grid_weight2 )

      if ( .not. more ) then
        exit
      end if

    end do

  end do

  return
end
subroutine sparse_grid_laguerre_index ( dim_num, level_max, point_num, &
  grid_index, grid_base )

!*****************************************************************************80
!
!! SPARSE_GRID_LAGUERRE_INDEX indexes points in a sparse Gauss-Laguerre grid.
!
!  Discussion:
!
!    The sparse grid is assumed to be formed from 1D Gauss-Laguerre rules.
!
!    The necessary dimensions of GRID_INDEX can be determined by
!    calling SPARSE_GRID_LAGUERRE_SIZE first.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 July 2008
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
!    Output, integer ( kind = 4 ) GRID_BASE(DIM_NUM,POINT_NUM), a list of
!    the orders of the Gauss-Laguerre rules associated with each point
!    and dimension.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) factor
  integer ( kind = 4 ) grid_base(dim_num,point_num)
  integer ( kind = 4 ), dimension ( dim_num ) :: grid_base2
  integer ( kind = 4 ) grid_index(dim_num,point_num)
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: grid_index2
  integer ( kind = 4 ), allocatable, dimension ( : ) :: grid_level
  integer ( kind = 4 ) h
  integer ( kind = 4 ) j
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
!
      call level_to_order_open ( dim_num, level_1d, order_1d )

      grid_base2(1:dim_num) = order_1d(1:dim_num)
!
!  The product of the 1D orders gives us the number of points in this subgrid.
!
      order_nd = product ( order_1d(1:dim_num) )

      allocate ( grid_index2(1:dim_num,1:order_nd) )
      allocate ( grid_level(1:order_nd) )
!
!  The inner (hidden) loop generates all points corresponding to given grid.
!  The grid indices will be between 1 and ORDER_1D(DIM).
!
      call multigrid_index_one ( dim_num, order_1d, order_nd, grid_index2 )
!
!  Only keep those points which first appear on this level.
!
      do point = 1, order_nd

        point_num2 = point_num2 + 1

        grid_index(1:dim_num,point_num2) = grid_index2(1:dim_num,point)
        grid_base(1:dim_num,point_num2) = grid_base2(1:dim_num)

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
subroutine sparse_grid_laguerre_size ( dim_num, level_max, point_num )

!*****************************************************************************80
!
!! SPARSE_GRID_LAGUERRE_SIZE sizes a sparse grid of Gauss-Laguerre points.
!
!  Discussion:
!
!    The grid is defined as the sum of the product rules whose LEVEL
!    satisfies:
!
!      LEVEL_MIN <= LEVEL <= LEVEL_MAX.
!
!    where LEVEL_MAX is user specified, and
!
!      LEVEL_MIN = max ( 0, LEVEL_MAX + 1 - DIM_NUM ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 July 2008
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
  logical more
  integer ( kind = 4 ) num
  integer ( kind = 4 ) order_1d(dim_num)
  integer ( kind = 4 ) order_nd
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
  point_num = 0

  level_min = max ( 0, level_max + 1 - dim_num )

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
!    27 August 2007
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
