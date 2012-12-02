program main

!*****************************************************************************80
!
!! MAIN is the main program for SPARSE_GRID_GL_DATASET.
!
!  Discussion:
!
!    This program computes a sparse grid quadrature rule based on a 1D
!    Gauss-Legendre rule and writes it to a file..
!
!    The user specifies:
!    * the spatial dimension of the quadrature region,
!    * the level that defines the Smolyak grid.
!
!  Modified:
!
!    03 October 2007
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
  integer   ( kind = 4 ), parameter :: i4_1 = 1
  integer   ( kind = 4 ), parameter :: i4_10 = 10
  integer   ( kind = 4 ) iarg
  integer   ( kind = 4 ) iargc
  integer   ( kind = 4 ) ierror
  integer   ( kind = 4 ) last
  integer   ( kind = 4 ) level_max
  character ( len = 2 ) level_max_string
  integer   ( kind = 4 ) level_min
  integer   ( kind = 4 ) point
  integer   ( kind = 4 ) point_num
  real      ( kind = 8 ), allocatable, dimension ( :, : ) :: r
  character ( len = 80 ) r_filename
  character ( len = 80 ) string
  real      ( kind = 8 ), allocatable, dimension ( : ) :: w
  character ( len = 80 ) w_filename
  real      ( kind = 8 ) weight_sum
  real      ( kind = 8 ), allocatable, dimension ( :, : ) :: x
  character ( len = 80 ) x_filename

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPARSE_GRID_GL_DATASET'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Compute the abscissas and weights of a quadrature rule'
  write ( *, '(a)' ) '  associated with a sparse grid derived from a Smolyak'
  write ( *, '(a)' ) '  construction based on 1D Gauss-Legendre rules.'
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
  write ( *, '(a)' ) '    (1) "gl_d?_level?_r.txt", the ranges;'
  write ( *, '(a)' ) '    (2) "gl_d?_level?_w.txt", the weights;'
  write ( *, '(a)' ) '    (3) "gl_d?_level?_x.txt", the abscissas.'
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
  call sparse_grid_gl_size ( dim_num, level_max, point_num )

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

  call sparse_grid_gl ( dim_num, level_max, point_num, w, x )

  call r8mat_transpose_print_some ( dim_num, point_num, x, &
    i4_1, i4_1, dim_num, i4_10, '  First 10 grid points:' )

  call r8vec_print_some ( point_num, w, i4_1, i4_10, &
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
  else
    write ( dim_num_string, '(i2)' ) dim_num
  end if

  r_filename = 'gl_d' // &
               trim ( dim_num_string )   // '_level' // &
               trim ( level_max_string ) // '_r.txt'

  w_filename = 'gl_d' // &
               trim ( dim_num_string )   // '_level' // &
               trim ( level_max_string ) // '_w.txt'

  x_filename = 'gl_d' // &
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
  write ( *, '(a)' ) 'SPARSE_GRID_GL_DATASET:'
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
subroutine gl_abscissa ( dim_num, point_num, grid_index, grid_base, grid_point )

!*****************************************************************************80
!
!! GL_ABSCISSA sets abscissas for multidimensional Gauss-Legendre quadrature.
!
!  Discussion:
!
!    The "nesting" as it occurs for Gauss-Legendre sparse grids simply
!    involves the use of a specified set of permissible orders for the
!    rule.
!
!    The X array lists the (complete) Gauss-Legendre abscissas for rules
!    of order 1, 3, 7, 15, 31, 63, and 127 in order.
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
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, integer ( kind = 4 ) GRID_INDEX(DIM_NUM,POINT_NUM), for each
!    point and dimension, the index of the Gauss-Legendre abscissa.
!
!    Input, integer ( kind = 4 ) GRID_BASE(DIM_NUM), the "base" of the
!    Gauss-Legendre rule being used in each dimension.
!
!    Output, real ( kind = 8 ) GRID_POINT(DIM_NUM), the grid points of
!    Gauss-Legendre abscissas.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) grid_base(dim_num)
  integer ( kind = 4 ) grid_index(dim_num,point_num)
  real    ( kind = 8 ) grid_point(dim_num,point_num)
  integer ( kind = 4 ) i4_log_2
  integer ( kind = 4 ) level
  integer ( kind = 4 ) point
  integer ( kind = 4 ) pointer
  integer ( kind = 4 ), dimension ( 0:7 ) :: skip = (/ &
    0, 1, 4, 11, 26, 57, 120, 247 /)
  real    ( kind = 8 ), dimension ( 247 ) :: x = (/ &
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
       0.51772881329003324812447758452632D+00, &
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
       0.38410279579151693577907781452239D+00, &
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
!! GL_WEIGHTS returns weights for certain Gauss-Legendre quadrature rules.
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
!    The weights are positive, symmetric and should sum to 2.
!
  implicit none

  integer ( kind = 4 ) order

  real    ( kind = 8 ) weight(order)

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
subroutine index_level_gl ( level, level_max, dim_num, point_num, grid_index, &
  grid_base, grid_level )

!*****************************************************************************80
!
!! INDEX_LEVEL_GL: determine first level at which given index is generated.
!
!  Discussion:
!
!    We are constructing a sparse grid of Gauss-Legendre points.  The grid
!    is built up of product grids, with a characteristic LEVEL.
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
!    Gauss-Legendre abscissa is 0, the special "nested" value we need
!    to take care of.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 September 2007
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
!    Input, integer ( kind = 4 ) LEVEL, the level at which these points were
!    generated.  LEVEL_MIN <= LEVEL <= LEVEL_MAX.
!
!    Input, integer ( kind = 4 ) LEVEL_MAX, the maximum level.
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points to be tested.
!
!    Input, integer ( kind = 4 ) GRID_INDEX(DIM_NUM,POINT_NUM), the indices
!    of the points to be tested.
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

  level_min = max ( 0, level_max + 1 - dim_num )
!
!  If a point has a DIM-th component whose INDEX is 0, then the
!  value of LEVEL at which this point would first be generated is
!  less than LEVEL, unless the DIM-th component of GRID_BASE is 0.
!
  do point = 1, point_num

    grid_level(point) = max ( level, level_min )

    do dim = 1, dim_num
      if ( grid_index(dim,point) == 0 ) then
        grid_level(point) = max ( grid_level(point) - grid_base(dim), level_min )
      end if
    end do

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
!    Fejer Type 2, Gauss-Patterson, Newton Cotes Open, and
!    Newton Cotes Half Open rules.  It also can be used, partially, to describe
!    the growth of Gauss-Legendre rules.
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
subroutine multigrid_index_z ( dim_num, order_1d, order_nd, indx )

!*****************************************************************************80
!
!! MULTIGRID_INDEX_Z returns an indexed multidimensional grid.
!
!  Discussion:
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
!    11 September 2007
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
!  The values of A(DIM) are between 0 and ORDER_1D(DIM)-1 = N - 1 = 2 * M.
!  Subtracting M sets the range to -M to +M, as we wish.
!
    indx(1:dim_num,p) = a(1:dim_num) - ( order_1d(1:dim_num) - 1 ) / 2

  end do

  return
end
subroutine product_weight_gl ( dim_num, order_1d, order_nd, w_nd )

!*****************************************************************************80
!
!! PRODUCT_WEIGHT_GL: weights for a product Gauss-Legendre rule.
!
!  Discussion:
!
!    This routine computes the weights for a quadrature rule which is
!    a product of 1D Gauss-Legendre rules of varying order.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 September 2007
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
  real    ( kind = 8 ), allocatable, dimension ( : ) :: w_1d
  real    ( kind = 8 ) w_nd(order_nd)

  w_nd(1:order_nd) = 1.0D+00

  do dim = 1, dim_num

    allocate ( w_1d(1:order_1d(dim)) )

    call gl_weights ( order_1d(dim), w_1d )

    call r8vec_direct_product2 ( dim, order_1d(dim), w_1d, dim_num, &
      order_nd, w_nd )

    deallocate ( w_1d )

  end do

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

  real    ( kind = 8 ) a(m,n)
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
  real      ( kind = 8 ) table(m,n)
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
  real    ( kind = 8 ) factor_value(factor_order)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ), save :: rep
  integer ( kind = 4 ), save :: skip
  integer ( kind = 4 ) start
  real    ( kind = 8 ) w(point_num)

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

  real    ( kind = 8 ) a(n)
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
subroutine sparse_grid_gl ( dim_num, level_max, point_num, grid_weight, &
  grid_point )

!****************************************************************************80
!
!! SPARSE_GRID_GL computes a sparse grid of Gauss-Legendre points.
!
!  Discussion:
!
!    The quadrature rule is associated with a sparse grid derived from
!    a Smolyak construction using a 1D Gauss-Legendre quadrature rule.
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
!    Input, integer ( kind = 4 ) LEVEL_MAX, controls the size of the
!    sparse grid.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points in the grid,
!    as determined by SPARSE_GRID_GL_SIZE.
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
  real    ( kind = 8 ) grid_point(dim_num,point_num)
  real    ( kind = 8 ) grid_point_temp(dim_num)
  real    ( kind = 8 ) grid_weight(point_num)
  real    ( kind = 8 ), allocatable, dimension ( : ) :: grid_weight2
  integer ( kind = 4 ) h
  integer ( kind = 4 ), parameter :: i4_1 = 1
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
      call product_weight_gl ( dim_num, order_1d, order_nd, grid_weight2 )
!
!  Now determine the coefficient of the weight.
!
      coeff = (-1)**( level_max - level ) &
        * choose ( dim_num - 1, level_max - level )
!
!  The inner (hidden) loop generates all points corresponding to given grid.
!  The grid indices will be between -M to +M, where 2*M + 1 = ORDER_1D(DIM).
!
      call multigrid_index_z ( dim_num, order_1d, order_nd, grid_index2 )
!
!  Determine the first level of appearance of each of the points.
!  This allows us to flag certain points as being repeats of points
!  generated on a grid of lower level.
!
!  This is SLIGHTLY tricky.
!
      call index_level_gl ( level, level_max, dim_num, order_nd, grid_index2, &
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

          call gl_abscissa ( dim_num, i4_1, grid_index2(1:dim_num,point), &
            grid_base2(1:dim_num), grid_point(1:dim_num,point_num2) )

          grid_weight(point_num2) = &
            real ( coeff, kind = 8 ) * grid_weight2(point)
!
!  or an already existing point (create point temporarily, find match,
!  add weight to matched point's weight).
!
        else

          call gl_abscissa ( dim_num, i4_1, grid_index2(1:dim_num,point), &
            grid_base2(1:dim_num), grid_point_temp(1:dim_num) )

          point3 = -1

          do point2 = 1, point_num2
            if ( all ( grid_point(1:dim_num,point2) == grid_point_temp(1:dim_num) ) ) then
              point3 = point2
              exit
            end if
          end do

          if ( point3 == -1 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'SPARSE_GRID_GL - Fatal error!'
            write ( *, '(a)' ) '  Could not match point.'
            stop
          end if

          grid_weight(point3) = grid_weight(point3) + &
            real ( coeff, kind = 8 ) * grid_weight2(point)

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

  return
end
subroutine sparse_grid_gl_index ( dim_num, level_max, point_num, grid_index, &
  grid_base )

!*****************************************************************************80
!
!! SPARSE_GRID_GL_INDEX indexes the points forming a sparse grid of GL points.
!
!  Discussion:
!
!    The sparse grid is assumed to be formed from 1D Gauss-Legendre rules
!    of ODD order, which have the property that only the central abscissa,
!    X = 0.0, is "nested".
!
!    The necessary dimensions of GRID_INDEX can be determined by
!    calling SPARSE_GRID_GL_SIZE first.
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
!    Output, integer ( kind = 4 ) GRID_BASE(DIM_NUM,POINT_NUM), a list of
!    the orders of the Gauss-Legendre rules associated with each point
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
      call multigrid_index_z ( dim_num, order_1d, order_nd, grid_index2 )
!
!  Determine the first level of appearance of each of the points.
!  This allows us to flag certain points as being repeats of points
!  generated on a grid of lower level.
!
!  This is SLIGHTLY tricky.
!
      call index_level_gl ( level, level_max, dim_num, order_nd, grid_index2, &
        grid_base2, grid_level )
!
!  Only keep those points which first appear on this level.
!
      do point = 1, order_nd

        if ( grid_level(point) == level ) then

          point_num2 = point_num2 + 1

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

  return
end
subroutine sparse_grid_gl_size ( dim_num, level_max, point_num )

!*****************************************************************************80
!
!! SPARSE_GRID_GL_SIZE sizes a sparse grid of Gauss-Legendre points.
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
!    The grids are only very weakly nested, since Gauss-Legendre rules
!    only have the origin in common.
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

      do dim = 1, dim_num
!
!  If we can reduce the level in this dimension by 1 and
!  still not go below LEVEL_MIN.
!
        if ( level_min < level .and. 1 < order_1d(dim) ) then
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
!    Input, integer ( kind = 4 ) BASE(DIM_NUM), the bases to be used in
!    each dimension.  In dimension I, entries will range from 0 to BASE(I)-1.
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
