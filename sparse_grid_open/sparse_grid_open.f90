subroutine abscissa_level_open_nd ( level_max, dim_num, test_num, test_val, &
  test_level )

!*****************************************************************************80
!
!! ABSCISSA_LEVEL_OPEN_ND: first level at which given abscissa is generated.
!
!  Discussion:
!
!    We assume an underlying product grid.  In each dimension, this product
!    grid has order 2**(LEVEL_MAX+1) - 1.
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
!    18 April 2007
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
!    the points to be tested.  Normally, each index would be between 
!    0 and 2**(LEVEL_MAX+1).
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
!    04 April 2007
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
subroutine gl_abscissa ( dim_num, point_num, grid_index, grid_point )

!*****************************************************************************80
!
!! GL_ABSCISSA sets abscissas for "nested" Gauss-Legendre quadrature.
!
!  Discussion:
!
!    The "nesting" as it occurs for Gauss-Legendre sparse grids simply
!    involves the use of a specified set of permissible orders for the
!    rule.  
!
!    The XTAB array lists the Gauss-Legendre abscissas for rules of order
!    1, 3, 5, 9, 17, 33 and 65, in order. 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 April 2007
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
!    Input, integer ( kind = 4 ) GRID_INDEX(DIM_NUM,POINT_NUM), indices that
!    indicate the Gauss-Legendre abscissa to be used for each component
!    of each point.  Each index should be between 1 and 133, indicating
!    a particular abscissa.
!
!    Output, real ( kind = 8 ) GRID_POINT(DIM_NUM), the grid points of
!    Gauss-Legendre abscissas.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) dim
  integer ( kind = 4 ) grid_index(dim_num,point_num)
  real ( kind = 8 ) grid_point(dim_num,point_num)
  real ( kind = 8 ), dimension ( 133 ) :: xtab = (/ &
       0.0D+00, &
     - 0.774596669241483377035853079956D+00, &
       0.0D+00, &
       0.774596669241483377035853079956D+00, &
     - 0.906179845938663992797626878299D+00, &
     - 0.538469310105683091036314420700D+00, &
       0.0D+00, &
       0.538469310105683091036314420700D+00, &
       0.906179845938663992797626878299D+00, &
     - 0.968160239507626089835576202904D+00, &
     - 0.836031107326635794299429788070D+00, &
     - 0.613371432700590397308702039341D+00, &
     - 0.324253423403808929038538014643D+00, &
       0.0D+00, &
       0.324253423403808929038538014643D+00, &
       0.613371432700590397308702039341D+00, &
       0.836031107326635794299429788070D+00, &
       0.968160239507626089835576202904D+00, &
     - 0.990575475314417335675434019941D+00, &
     - 0.950675521768767761222716957896D+00, &
     - 0.880239153726985902122955694488D+00, &
     - 0.781514003896801406925230055520D+00, &
     - 0.657671159216690765850302216643D+00, &
     - 0.512690537086476967886246568630D+00, &
     - 0.351231763453876315297185517095D+00, &
     - 0.178484181495847855850677493654D+00, &
       0.0D+00, &
       0.178484181495847855850677493654D+00, &
       0.351231763453876315297185517095D+00, &
       0.512690537086476967886246568630D+00, &
       0.657671159216690765850302216643D+00, &
       0.781514003896801406925230055520D+00, &
       0.880239153726985902122955694488D+00, &
       0.950675521768767761222716957896D+00, &
       0.990575475314417335675434019941D+00, &
      -0.9974246942464552D+00, &    
      -0.9864557262306425D+00, &
      -0.9668229096899927D+00, &
      -0.9386943726111684D+00, &    
      -0.9023167677434336D+00, &    
      -0.8580096526765041D+00, &    
      -0.8061623562741665D+00, &    
      -0.7472304964495622D+00, &    
      -0.6817319599697428D+00, &    
      -0.6102423458363790D+00, &    
      -0.5333899047863476D+00, &    
      -0.4518500172724507D+00, &    
      -0.3663392577480734D+00, &    
      -0.2776090971524970D+00, &    
      -0.1864392988279916D+00, &    
      -0.09363106585473338D+00, &
       0.0D+00, &
       0.09363106585473338D+00, &
       0.1864392988279916D+00, &    
       0.2776090971524970D+00, &    
       0.3663392577480734D+00, &    
       0.4518500172724507D+00, &    
       0.5333899047863476D+00, &    
       0.6102423458363790D+00, &    
       0.6817319599697428D+00, &    
       0.7472304964495622D+00, &    
       0.8061623562741665D+00, &    
       0.8580096526765041D+00, &    
       0.9023167677434336D+00, &    
       0.9386943726111684D+00, &    
       0.9668229096899927D+00, &    
       0.9864557262306425D+00, &    
       0.9974246942464552D+00, &    
      -0.9993260970754129D+00, &    
      -0.9964509480618492D+00, &    
      -0.9912852761768016D+00, &    
      -0.9838398121870350D+00, &    
      -0.9741315398335512D+00, &    
      -0.9621827547180553D+00, &    
      -0.9480209281684076D+00, &    
      -0.9316786282287494D+00, &    
      -0.9131934405428462D+00, &    
      -0.8926078805047389D+00, &    
      -0.8699692949264071D+00, &    
      -0.8453297528999303D+00, &    
      -0.8187459259226514D+00, &    
      -0.7902789574921218D+00, &    
      -0.7599943224419998D+00, &    
      -0.7279616763294247D+00, &    
      -0.6942546952139916D+00, &    
      -0.6589509061936252D+00, &    
      -0.6221315090854003D+00, &    
      -0.5838811896604873D+00, &    
      -0.5442879248622271D+00, &    
      -0.5034427804550069D+00, &    
      -0.4614397015691450D+00, &    
      -0.4183752966234090D+00, &    
      -0.3743486151220660D+00, &    
      -0.3294609198374864D+00, &    
      -0.2838154539022487D+00, &    
      -0.2375172033464168D+00, &    
      -0.1906726556261428D+00, &    
      -0.1433895546989752D+00, &    
      -0.9577665320919751D-01, &
      -0.4794346235317186D-01, &
       0.0D+00, &
       0.4794346235317186D-01, &
       0.9577665320919751D-01, &
       0.1433895546989752D+00, &    
       0.1906726556261428D+00, &    
       0.2375172033464168D+00, &    
       0.2838154539022487D+00, &    
       0.3294609198374864D+00, &    
       0.3743486151220660D+00, &    
       0.4183752966234090D+00, &    
       0.4614397015691450D+00, &    
       0.5034427804550069D+00, &    
       0.5442879248622271D+00, &    
       0.5838811896604873D+00, &    
       0.6221315090854003D+00, &    
       0.6589509061936252D+00, &    
       0.6942546952139916D+00, &    
       0.7279616763294247D+00, &    
       0.7599943224419998D+00, &    
       0.7902789574921218D+00, &    
       0.8187459259226514D+00, &    
       0.8453297528999303D+00, &    
       0.8699692949264071D+00, &    
       0.8926078805047389D+00, &    
       0.9131934405428462D+00, &    
       0.9316786282287494D+00, &    
       0.9480209281684076D+00, &    
       0.9621827547180553D+00, &    
       0.9741315398335512D+00, &    
       0.9838398121870350D+00, &    
       0.9912852761768016D+00, &    
       0.9964509480618492D+00, &    
       0.9993260970754129D+00 /)    

  if ( any ( grid_index(1:dim_num,1:point_num) < 1 ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GL_ABSCISSA - Fatal error!'
    write ( *, '(a)' ) '  Some index values are less than 1.'
    stop
  else if ( any ( 127 < grid_index(1:dim_num,1:point_num) ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GL_ABSCISSA - Fatal error!'
    write ( *, '(a)' ) '  Some index values are greater than 127.'
    stop
  end if

  do dim = 1, dim_num
    grid_point(dim,1:point_num) = xtab ( grid_index(dim,1:point_num) ) 
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
!    02 March 1999
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
!    18 April 2007
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
!    Input, integer ( kind = 4 ) ORDER, the order of the rule.
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
!    Fejer Type 2, Gauss-Patterson, and Newton Cotes Open
!    rules.  It also can be used, partially, to describe
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
!    For the Fejer Type 2, Gauss-Patterson, and Newton Cotes Open, 
!    rules, the point growth is
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
subroutine levels_open_index ( dim_num, level_max, point_num, grid_index )

!*****************************************************************************80
!
!! LEVELS_OPEN_INDEX computes open grids with 0 <= LEVEL <= LEVEL_MAX.
!
!  Discussion:
!
!    The necessary dimensions of GRID_INDEX can be determined by 
!    calling LEVELS_OPEN_INDEX_SIZE first.
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
!  Transform each 1D level to a corresponding 1D order, roughly
!  ORDER = 2**(LEVEL+1)+1.
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

        if ( all ( mod ( grid_index2(1:dim_num,point), 2 ) == 1 ) ) then

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
!    23 May 2007
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
!    08 June 2007
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
!    Input, integer ( kind = 4 ) ORDER, the order of the rule.
!    1 <= ORDER.
!
!    Input, integer ( kind = 4 ) I, the index of the desired abscissa.  
!    1 <= I <= ORDER.
!
!    Output, real ( kind = 8 ) NCO_ABSCISSA, the value of the I-th 
!    abscissa in the Newton Cotes Open rule of order ORDER.
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

  character              ch
  integer   ( kind = 4 ) get
  integer   ( kind = 4 ) put
  character ( len = * )  s
  integer   ( kind = 4 ) s_length
  character, parameter :: tab = achar ( 9 )

  put = 0
  s_length = len_trim ( s )

  do get = 1, s_length

    ch = s(get:get)
 
    if ( ch /= ' ' .and. ch /= tab ) then
      put = put + 1
      s(put:put) = ch
    end if
 
  end do
 
  s(put+1:s_length) = ' '
  
  return
end
subroutine sparse_grid_f2s_size ( dim_num, level_max, point_num )

!*****************************************************************************80
!
!! SPARSE_GRID_F2S_SIZE sizes a sparse grid using Fejer Type 2 Slow rules.
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
  logical              more
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
subroutine sparse_grid_gps_size ( dim_num, level_max, point_num )

!*****************************************************************************80
!
!! SPARSE_GRID_GPS_SIZE sizes a sparse grid using Gauss-Patterson-Slow rules.
!
!  Discussion:
!
!    The Gauss-Patterson-Slow family assumes that, for the underlying 1D
!    rules, a precision of 2*L+1 is needed at level L.  Therefore, the
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
  logical              more
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
subroutine sparse_grid_onn_size ( dim_num, level_max, point_num )

!*****************************************************************************80
!
!! SPARSE_GRID_ONN_SIZE sizes a sparse grid using Open Non-Nested rules.
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
!    This routine assumes a linear growth rule is being used.
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
  logical              more
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
subroutine sparse_grid_own_size ( dim_num, level_max, point_num )

!*****************************************************************************80
!
!! SPARSE_GRID_OWN_SIZE sizes a sparse grid using Open Weakly Nested rules.
!
!  Discussion:
!
!    This calculation is much faster than a previous method.  
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
  logical              more
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
!    25 May 2007
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
