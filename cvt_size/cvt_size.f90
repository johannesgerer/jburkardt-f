subroutine angle_to_rgb ( angle, r, g, b )

!*****************************************************************************80
!
!! ANGLE_TO_RGB returns a color on the perimeter of the color hexagon.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 December 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ANGLE, the angle in the color hexagon.
!    The sextants are
!    defined by the following points:
!        0 degrees, 1, 0, 0, red;
!       60 degrees, 1, 1, 0, yellow;
!      120 degrees, 0, 1, 0, green;
!      180 degrees, 0, 1, 1, cyan;
!      240 degrees, 0, 0, 1, blue;
!      300 degrees, 1, 0, 1, magenta.
!
!    Output, real ( kind = 8 ) R, G, B, RGB specifications for the color that
!    lies at the given angle, on the perimeter of the color hexagon.  One
!    value will be 1, and one value will be 0.
!
  implicit none

  real ( kind = 8 ) angle
  real ( kind = 8 ) angle2
  real ( kind = 8 ) b
  real ( kind = 8 ) g
  real ( kind = 8 ), parameter :: degrees_to_radians = &
    3.14159265D+00 / 180.0D+00
  real ( kind = 8 ) r

  angle = mod ( angle, 360.0D+00 )

  if ( angle < 0.0D+00 ) then
    angle = angle + 360.0D+00
  end if

  if ( angle <= 60.0D+00 ) then

    angle2 = degrees_to_radians * 3.0D+00 * angle / 4.0D+00
    r = 1.0D+00
    g = tan ( angle2 )
    b = 0.0D+00

  else if ( angle <= 120.0D+00 ) then

    angle2 = degrees_to_radians * 3.0D+00 * angle / 4.0D+00
    r = cos ( angle2 ) / sin ( angle2 )
    g = 1.0D+00
    b = 0.0D+00

  else if ( angle <= 180.0D+00 ) then

    angle2 = degrees_to_radians * 3.0D+00 * ( angle - 120.0D+00 ) / 4.0D+00
    r = 0.0D+00
    g = 1.0D+00
    b = tan ( angle2 )

  else if ( angle <= 240.0D+00 ) then

    angle2 = degrees_to_radians * 3.0D+00 * ( angle - 120.0D+00 ) / 4.0D+00
    r = 0.0D+00
    g = cos ( angle2 ) / sin ( angle2 )
    b = 1.0D+00

  else if ( angle <= 300.0D+00 ) then

    angle2 = degrees_to_radians * 3.0D+00 * ( angle - 240.0D+00 ) / 4.0D+00
    r = tan ( angle2 )
    g = 0.0D+00
    b = 1.0D+00

  else if ( angle <= 360.0D+00 ) then

    angle2 = degrees_to_radians * 3.0D+00 * ( angle - 240.0D+00 ) / 4.0D+00
    r = 1.0D+00
    g = 0.0D+00
    b = cos ( angle2 ) / sin ( angle2 )

  end if

  return
end
function ch_is_digit ( c )

!*****************************************************************************80
!
!! CH_IS_DIGIT returns .TRUE. if a character is a decimal digit.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C, the character to be analyzed.
!
!    Output, logical CH_IS_DIGIT, .TRUE. if C is a digit, .FALSE. otherwise.
!
  implicit none

  character c
  logical ch_is_digit

  if ( lge ( c, '0' ) .and. lle ( c, '9' ) ) then
    ch_is_digit = .true.
  else
    ch_is_digit = .false.
  end if

  return
end
subroutine ch_to_digit ( c, digit )

!*****************************************************************************80
!
!! CH_TO_DIGIT returns the integer value of a base 10 digit.
!
!  Example:
!
!     C   DIGIT
!    ---  -----
!    '0'    0
!    '1'    1
!    ...  ...
!    '9'    9
!    ' '    0
!    'X'   -1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C, the decimal digit, '0' through '9' or blank
!    are legal.
!
!    Output, integer DIGIT, the corresponding integer value.  If C was
!    'illegal', then DIGIT is -1.
!
  implicit none

  character c
  integer digit

  if ( lge ( c, '0' ) .and. lle ( c, '9' ) ) then

    digit = ichar ( c ) - 48

  else if ( c == ' ' ) then

    digit = 0

  else

    digit = -1

  end if

  return
end
subroutine cvt_size_iteration ( dim_num, box_min, box_max, cell_num, &
  cell_generator, cell_weight, sample_num, cell_centroid, &
  cell_volume, region_volume )

!******************************************************************************
!
!! CVT_SIZE_ITERATION takes one step of the CVT size iteration.
!
!  Discussion:
!
!    The routine is given a set of points, called "generators", which
!    define a tessellation of the region into Voronoi cells.  Each point
!    defines a cell.  Each cell, in turn, has a centroid, but it is
!    unlikely that the centroid and the generator coincide.
!
!    Each time this CVT iteration is carried out, an attempt is made
!    to modify the generators in such a way that they are closer and
!    closer to being the centroids of the Voronoi cells they generate.
!
!    A large number of sample points are generated, and the nearest generator
!    is determined.  A count is kept of how many points were nearest to each
!    generator.  Once the sampling is completed, the centroids are computed.
!
!    If this routine is called in a loop, the centroids from the previous
!    step should be input as the generators for the next step.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 December 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, real ( kind = 8 ) BOX_MIN(DIM_NUM), BOX_MAX(DIM_NUM), the coordinates
!    of the two extreme corners of the bounding box.
!
!    Input, integer CELL_NUM, the number of Voronoi cells.
!
!    Input, real ( kind = 8 ) CELL_GENERATOR(DIM_NUM,CELL_NUM), the Voronoi
!    cell generators.
!
!    Input, real ( kind = 8 ) CELL_WEIGHT(CELL_NUM), the cell weights.
!
!    Input, integer SAMPLE_NUM, the number of sample points.
!
!    Output, real ( kind = 8 ) CELL_CENTROID(DIM_NUM,CELL_NUM), the Voronoi
!    cell centroids.
!
  implicit none

  integer cell_num
  integer dim_num

  real ( kind = 8 ) box_max(1:dim_num)
  real ( kind = 8 ) box_min(1:dim_num)
  integer cell_count(cell_num)
  real ( kind = 8 ) cell_generator(dim_num,cell_num)
  real ( kind = 8 ) cell_centroid(dim_num,cell_num)
  real ( kind = 8 ) cell_volume(cell_num)
  real ( kind = 8 ) cell_weight(cell_num)
  logical, parameter :: debug = .false.
  integer j
  integer nearest
  integer ngen
  real ( kind = 8 ) region_volume
  integer sample_num
  real ( kind = 8 ) x(dim_num)

  cell_centroid(1:dim_num,1:cell_num) = 0.0D+00
  cell_count(1:cell_num) = 0

  do j = 1, sample_num
!
!  Generate a sampling point X.
!
    call region_sampler ( dim_num, box_min, box_max, x, ngen )
!
!  Find the nearest cell generator.
!
    call find_closest_size ( dim_num, x, cell_num, cell_generator, &
      cell_weight, nearest )
!
!  Add X to the averaging data for CELL_GENERATOR(*,NEAREST).
!
    if ( nearest < 1 .or. cell_num < nearest ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CVT_SIZE_ITERATION - Fatal error!'
      write ( *, '(a,i6)' ) 'NEAREST = ',nearest
      stop
    end if

    cell_centroid(1:dim_num,nearest) = &
      cell_centroid(1:dim_num,nearest) + x(1:dim_num)

    cell_count(nearest) = cell_count(nearest) + 1

  end do

  do j = 1, cell_num

    if ( cell_count(j) /= 0 ) then

      cell_centroid(1:dim_num,j) = cell_centroid(1:dim_num,j) &
        / real ( cell_count(j), kind = 8 )

    end if

  end do

  cell_volume(1:cell_num) = cell_count(1:cell_num) * region_volume &
    / real ( sample_num, kind = 8 )

  return
end
subroutine digit_inc ( c )

!*****************************************************************************80
!
!! DIGIT_INC increments a decimal digit.
!
!  Example:
!
!    Input  Output
!    -----  ------
!    '0'    '1'
!    '1'    '2'
!    ...
!    '8'    '9'
!    '9'    '0'
!    'A'    'A'
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character C, a digit to be incremented.
!
  implicit none

  character c
  integer digit

  call ch_to_digit ( c, digit )

  if ( digit == -1 ) then
    return
  end if

  digit = digit + 1

  if ( digit == 10 ) then
    digit = 0
  end if

  call digit_to_ch ( digit, c )

  return
end
subroutine digit_to_ch ( digit, c )

!*****************************************************************************80
!
!! DIGIT_TO_CH returns the character representation of a decimal digit.
!
!  Example:
!
!    DIGIT   C
!    -----  ---
!      0    '0'
!      1    '1'
!    ...    ...
!      9    '9'
!     17    '*'
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIGIT, the digit value between 0 and 9.
!
!    Output, character C, the corresponding character, or '*' if DIGIT
!    was illegal.
!
  implicit none

  character c
  integer digit

  if ( 0 <= digit .and. digit <= 9 ) then

    c = char ( digit + 48 )

  else

    c = '*'

  end if

  return
end
subroutine file_name_inc ( file_name )

!*****************************************************************************80
!
!! FILE_NAME_INC generates the next filename in a series.
!
!  Discussion:
!
!    It is assumed that the digits in the name, whether scattered or
!    connected, represent a number that is to be increased by 1 on
!    each call.  If this number is all 9's on input, the output number
!    is all 0's.  Non-numeric letters of the name are unaffected, and
!    if the name contains no digits, then nothing is done.
!
!  Example:
!
!      Input          Output
!      -----          ------
!      a7to11.txt     a7to12.txt
!      a7to99.txt     a8to00.txt
!      a9to99.txt     a0to00.txt
!      cat.txt        cat.txt
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) FILE_NAME.
!    On input, a character string to be incremented.
!    On output, the incremented string.
!
  implicit none

  character c
  logical ch_is_digit
  character ( len = * ) file_name
  integer i
  integer lens

  lens = len_trim ( file_name )

  do i = lens, 1, -1

    c = file_name(i:i)

    if ( ch_is_digit ( c ) ) then

      call digit_inc ( c )

      file_name(i:i) = c

      if ( c /= '0' ) then
        return
      end if

    end if

  end do

  return
end
subroutine find_closest_size ( dim_num, x, cell_num, cell_generator, &
  cell_weight, nearest )

!******************************************************************************
!
!! FIND_CLOSEST_SIZE finds the Voronoi cell generator closest to a point X.
!
!  Discussion:
!
!    This routine finds the closest Voronoi cell generator by checking every
!    one.  For problems with many cells, this process can take the bulk
!    of the CPU time.  Other approaches, which group the cell generators into
!    bins, can run faster by a large factor.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 September 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, real ( kind = 8 ) X(DIM_NUM), the point to be checked.
!
!    Input, integer CELL_NUM, the number of cell generators.
!
!    Input, real ( kind = 8 ) CELL_GENERATOR(DIM_NUM,CELL_NUM), the cell generators.
!
!    Input, real ( kind = 8 ) CELL_WEIGHT(CELL_NUM), the cell weights.
!
!    Output, integer NEAREST, the index of the nearest cell generators.
!
  implicit none

  integer cell_num
  integer dim_num

  real ( kind = 8 ) cell_generator(dim_num,cell_num)
  real ( kind = 8 ) cell_weight(cell_num)
  real ( kind = 8 ) distance
  real ( kind = 8 ) dist_sq
  integer i
  integer nearest
  real ( kind = 8 ) x(dim_num)

  nearest = 0
  distance = huge ( distance )

  do i = 1, cell_num

    dist_sq = sum ( ( cell_generator(1:dim_num,i) - x(1:dim_num) )**2 )

    dist_sq = sqrt ( dist_sq ) / cell_weight(i)

    if ( dist_sq < distance ) then
      distance = dist_sq
      nearest = i
    end if

  end do

  return
end
subroutine generator_init ( dim_num, box_min, box_max, cell_num, &
  cell_generator )

!******************************************************************************
!
!! GENERATOR_INIT initializes the Voronoi cell generators.
!
!  Discussion:
!
!    The points initialized here will be used to generate a tessellation
!    of the region into Voronoi cells.  Each generator point defines a
!    cell.  The CVT algorithm will try to modify these initial generators
!    in such a way that they are also the centroids of the cells they generate.
!
!    It is probably better to use Halton points for the centroids than
!    uniform random values, in the sense that the algorithm will probably
!    converge more quickly.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 September 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, real ( kind = 8 ) BOX_MIN(DIM_NUM), BOX_MAX(DIM_NUM), the coordinates
!    of the two extreme corners of the bounding box.
!
!    Input, integer CELL_NUM, the number of Voronoi cells.
!
!    Output, real ( kind = 8 ) CELL_GENERATOR(DIM_NUM,CELL_NUM), the
!    Voronoi cell generators.
!
  implicit none

  integer cell_num
  integer dim_num

  real ( kind = 8 ) box_max(dim_num)
  real ( kind = 8 ) box_min(dim_num)
  real ( kind = 8 ) cell_generator(dim_num,cell_num)
  integer i
  integer ngen

  do i = 1, cell_num

    call region_sampler ( dim_num, box_min, box_max, cell_generator(1,i), ngen )

  end do

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
!    Output, integer IUNIT, the free unit number.
!
  implicit none

  integer i
  integer ios
  integer iunit
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
subroutine hexcol ( angle, r, g, b )

!*****************************************************************************80
!
!! HEXCOL returns a color on the perimeter of the color hexagon.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 December 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ANGLE, the angle in the color hexagon.  The
!    sextants are defined by the following points:
!        0 degrees, 1, 0, 0, red;
!       60 degrees, 1, 1, 0, yellow;
!      120 degrees, 0, 1, 0, green;
!      180 degrees, 0, 1, 1, cyan;
!      240 degrees, 0, 0, 1, blue;
!      300 degrees, 1, 0, 1, magenta.
!
!    Output, real ( kind = 8 ) R, G, B, RGB specifications for the color that
!    lies at the given angle, on the perimeter of the color hexagon.  One
!    value will be 1, and one value will be 0.
!
  implicit none

  real ( kind = 8 ) angle
  real ( kind = 8 ) angle2
  real ( kind = 8 ) b
  real ( kind = 8 ) g
  real ( kind = 8 ), parameter :: degrees_to_radians = &
    3.14159265D+00 / 180.0D+00
  real ( kind = 8 ) r

  angle = mod ( angle, 360.0D+00 )

  if ( angle < 0.0D+00 ) then
    angle = angle + 360.0D+00
  end if

  if ( angle <= 60.0D+00 ) then

    angle2 = degrees_to_radians * 3.0D+00 * angle / 4.0D+00
    r = 1.0D+00
    g = tan ( angle2 )
    b = 0.0D+00

  else if ( angle <= 120.0D+00 ) then

    angle2 = degrees_to_radians * 3.0D+00 * angle / 4.0D+00
    r = cos ( angle2 ) / sin ( angle2 )
    g = 1.0D+00
    b = 0.0D+00

  else if ( angle <= 180.0D+00 ) then

    angle2 = degrees_to_radians * 3.0D+00 * ( angle - 120.0D+00 ) / 4.0D+00
    r = 0.0D+00
    g = 1.0D+00
    b = tan ( angle2 )

  else if ( angle <= 240.0D+00 ) then

    angle2 = degrees_to_radians * 3.0D+00 * ( angle - 120.0D+00 ) / 4.0D+00
    r = 0.0D+00
    g = cos ( angle2 ) / sin ( angle2 )
    b = 1.0D+00

  else if ( angle <= 300.0D+00 ) then

    angle2 = degrees_to_radians * 3.0D+00 * ( angle - 240.0D+00 ) / 4.0D+00
    r = tan ( angle2 )
    g = 0.0D+00
    b = 1.0D+00

  else if ( angle <= 360.0D+00 ) then

    angle2 = degrees_to_radians * 3.0D+00 * ( angle - 240.0D+00 ) / 4.0D+00
    r = 1.0D+00
    g = 0.0D+00
    b = cos ( angle2 ) / sin ( angle2 )

  end if

  return
end
function i4_log_2 ( i )

!*****************************************************************************80
!
!! I4_LOG_2 returns the integer part of the logarithm base 2 of |I|.
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
!    Input, integer I, the number whose logarithm base 2 is desired.
!
!    Output, integer I4_LOG_2, the integer part of the logarithm base 2 of
!    the absolute value of I.
!    For positive I4_LOG_2(I), it should be true that
!      2**I4_LOG_2(X) <= |I| < 2**(I4_LOG_2(I)+1).
!    The special case of I4_LOG_2(0) returns -HUGE().
!
  implicit none

  integer i4_log_2
  integer i
  integer i_abs

  if ( i == 0 ) then

    i4_log_2 = - huge ( i4_log_2 )

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
subroutine i4_to_angle ( i, angle )

!*****************************************************************************80
!
!! I4_TO_ANGLE maps integers to points on a circle.
!
!  Discussion:
!
!    The angles are intended to be used to select colors on a color
!    hexagon whose 6 vertices are red, yellow, green, cyan, blue,
!    magenta.
!
!  Example:
!
!     I   X      ANGLE
!
!     0   0/3      0
!     1   1/3    120
!     2   2/3    240
!
!     3   1/6     60
!     4   3/6    180
!     5   5/6    300
!
!     6   1/12    30
!     7   3/12    90
!     8   5/12   150
!     9   7/12   210
!    10   9/12   270
!    11  11/12   330
!
!    12   1/24    15
!    13   3/24    45
!    14   5/24    75
!    etc
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 January 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer I, the index of the desired color.
!
!    Output, real ( kind = 8 ) ANGLE, an angle, measured in degrees,
!    between 0 and 360.
!
  implicit none

  real ( kind = 8 ) angle
  integer i
  integer i4_log_2
  integer i1
  integer i2
  integer i3
  integer i4

  if ( 0 <= abs ( i ) .and. abs ( i ) <= 2 ) then

    angle = 120.0D+00 * real ( abs ( i ), kind = 8 )

  else

    i1 = i4_log_2 ( abs ( i ) / 3 )
    i2 = abs ( i ) + 1 - 3 * 2**i1
    i3 = 2 * ( i2 - 1 ) + 1
    i4 = 3 * 2**( i1 + 1 )

    angle = 360.0D+00 * real ( i3, kind = 8 ) / real ( i4, kind = 8 )

  end if

  return
end
function i4_uniform ( a, b, seed )

!*****************************************************************************80
!
!! I4_UNIFORM returns a scaled pseudorandom I4.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ) value.
!
!    The pseudorandom number will be scaled to be uniformly distributed
!    between A and B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 November 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley Interscience, page 95, 1998.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) A, B, the limits of the interval.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, integer ( kind = 4 ) I4_UNIFORM, a number between A and B.
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) k
  real ( kind = 4 ) r
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) value

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_UNIFORM - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + 2147483647
  end if

  r = real ( seed, kind = 4 ) * 4.656612875E-10
!
!  Scale R to lie between A-0.5 and B+0.5.
!
  r = ( 1.0E+00 - r ) * ( real ( min ( a, b ), kind = 4 ) - 0.5E+00 ) &
    +             r   * ( real ( max ( a, b ), kind = 4 ) + 0.5E+00 )
!
!  Use rounding to convert R to an integer between A and B.
!
  value = nint ( r, kind = 4 )

  value = max ( value, min ( a, b ) )
  value = min ( value, max ( a, b ) )

  i4_uniform = value

  return
end
subroutine ppm_check_data ( r, g, b, ierror, maxcol, ncol, nrow )

!*****************************************************************************80
!
!! PPM_CHECK_DATA checks pixel data.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer R(NROW,NCOL), G(NROW,NCOL), B(NROW,NCOL), contains the
!    RGB pixel data.
!
!    Output, integer IERROR, error flag.
!    0, no error detected.
!    1, the data is illegal.
!
!    Input, integer MAXCOL, the maximum value.
!
!    Input, integer NCOL, NROW, the number of rows and columns of data.
!
  implicit none

  integer ncol
  integer nrow

  integer b(nrow,ncol)
  integer g(nrow,ncol)
  integer i
  integer ierror
  integer j
  integer maxcol
  integer r(nrow,ncol)

  ierror = 0
!
!  Make sure no color is negative.
!
  if ( minval ( r(1:nrow,1:ncol) ) < 0 .or. &
       minval ( g(1:nrow,1:ncol) ) < 0 .or. &
       minval ( b(1:nrow,1:ncol) ) < 0 ) then
    ierror = 1
    return
  end if
!
!  Make sure no color is greater than MAXCOL.
!
  if ( maxcol < maxval ( r(1:nrow,1:ncol) ) .or. &
       maxcol < maxval ( g(1:nrow,1:ncol) ) .or. &
       maxcol < maxval ( b(1:nrow,1:ncol) ) ) then
    ierror = 1
    return
  end if

  return
end
subroutine ppmb_write ( file_name, ierror, nrow, ncol, r, g, b )

!*****************************************************************************80
!
!! PPMB_WRITE writes a binary portable pixel map file.
!
!  Discussion:
!
!    PPM files can be viewed by XV.
!
!    Various programs can convert other formats to PPM format, including:
!
!      BMPTOPPM - Microsoft Windows BMP file.
!
!    A PPM file can also be converted to other formats, by programs:
!
!      PPMTOACAD - AutoCAD file
!      PPMTOGIF  - GIF file
!      PPMTOPGM  - Portable Gray Map file
!      PPMTOPICT - Macintosh PICT file
!      PPMTORGB3 - 3 Portable Gray Map files
!      PPMTOXPM  - X11 pixmap file
!      PPMTOYUV  - Abekas YUV file
!
!    DIRECT ACCESS is used for the output file just so that we can
!    avoid the internal carriage returns and things that FORTRAN
!    seems to want to add.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the file to which
!    the data should be written.
!
!    Output, integer IERROR, an error flag.
!    0, no error.
!    1, the data was illegal.
!    2, the file could not be opened.
!
!    Input, integer NROW, NCOL, the number of rows and columns of data.
!
!    Input, integer R(NROW,NCOL), G(NROW,NCOL), B(NROW,NCOL), contain
!    the red, green and blue values of each pixel.  These should all
!    be values between 0 and 255.
!
  implicit none

  integer ncol
  integer nrow

  integer b(nrow,ncol)
  logical, parameter :: debug = .false.
  character ( len = * ) file_name
  integer g(nrow,ncol)
  integer i
  integer ierror
  integer ios
  integer istring(17)
  integer j
  integer l
  integer maxcol
  integer nchar
  integer nval
  integer output_unit
  integer r(nrow,ncol)
  integer record
  integer record_length
  logical, parameter :: reverse = .false.
  character ( len = 68 ) string

  ierror = 0
!
!  Compute the maximum color value.
!
  maxcol = max ( &
    maxval ( r(1:nrow,1:ncol) ), &
    maxval ( g(1:nrow,1:ncol) ), &
    maxval ( b(1:nrow,1:ncol) ) )
!
!  Check that no color data exceeds 255.
!
  if ( 255 < maxcol ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PPMB_WRITE - Fatal error!'
    write ( *, '(a)' ) '  The color data exceeds 255!'
    write ( *, '(a,i12)' ) '  MAXCOL = ', maxcol
    ierror = 1
    return
  end if
!
!  Check the data.
!
  call ppm_check_data ( r, g, b, ierror, maxcol, ncol, nrow )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PPMB_WRITE - Fatal error!'
    write ( *, '(a)' ) '  Bad data detected by PPM_CHECK_DATA!'
    ierror = 1
    return
  end if
!
!  Open the file.
!
!  The smallest amount of information we can write at a time is
!  1 word = 4 bytes = 32 bits.
!
  call get_unit ( output_unit )
!
!  For the SGI:
!
  record_length = 4
!
!  For the DEC Alpha:
!
! record_length = 1

  open ( unit = output_unit, file = file_name, status = 'replace', &
    form = 'unformatted', access = 'direct', recl = record_length, &
    iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PPMB_WRITE - Fatal error!'
    write ( *, '(a)' ) '  Could not open the file.'
    ierror = 2
    return
  end if

  record = 0
!
!  Write the data.
!
  string = ' '

  string(1:4) = 'P6  '

  write ( string(5:9), '(i5)' ) ncol

  string(10:11) = ' '

  write ( string(12:16), '(i5)' ) nrow

  string(17:18) = ' '

  write ( string(19:23), '(i5)' ) maxcol

  string(24:24) = ' '

  nchar = 24

  do i = 1, nrow
    do j = 1, ncol

      if ( nchar == 68 ) then
        call s_to_i4vec ( string(1:nchar), nval, istring, reverse )
        do l = 1, nval
          record = record + 1
          write ( output_unit, rec = record ) istring(l)
        end do
        string = ' '
        nchar = 0
      end if

      nchar = nchar + 1
      string(nchar:nchar) = char ( r(i,j) )

      if ( nchar == 68 ) then
        call s_to_i4vec ( string(1:nchar), nval, istring, reverse )
        do l = 1, nval
          record = record + 1
          write ( output_unit, rec = record ) istring(l)
        end do
        string = ' '
        nchar = 0
      end if

      nchar = nchar + 1
      string(nchar:nchar) = char ( g(i,j) )

      if ( nchar == 68 ) then
        call s_to_i4vec ( string(1:nchar), nval, istring, reverse )
        do l = 1, nval
          record = record + 1
          write ( output_unit, rec = record ) istring(l)
        end do
        string = ' '
        nchar = 0
      end if

      nchar = nchar + 1
      string(nchar:nchar) = char ( b(i,j) )


    end do
  end do

  if ( 0 < nchar ) then
    call s_to_i4vec ( string(1:nchar), nval, istring, reverse )
    do l = 1, nval
      record = record + 1
      write ( output_unit, rec = record ) istring(l)
    end do
    string = ' '
    nchar = 0
  end if
!
!  Close the file.
!
  close ( unit = output_unit )
!
!  Report
!
  if ( debug ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PPMB_WRITE - Note:'
    write ( *, '(a)' ) '  The data was checked and written.'
    write ( *, '(a,i8)' ) '  Number of words =             ', record
    write ( *, '(a,i8)' ) '  Number of data rows NROW =    ', nrow
    write ( *, '(a,i8)' ) '  Number of data columns NCOL = ', ncol
    write ( *, '(a,i8)' ) '  Maximum color value MAXCOL =  ', maxcol
  end if

  return
end
function r82vec_dist_l2 ( a1, a2 )

!*****************************************************************************80
!
!! R82VEC_DIST_L2 returns the L2 distance between a pair of real 2-vectors.
!
!  Discussion:
!
!    The vector L2 norm is defined as:
!
!      norm = sqrt ( sum ( 1 <= I <= N ) A(I)**2 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 December 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A1(2), A2(2), the vectors whose L2 distance is desired.
!
!    Output, real ( kind = 8 ) R82VEC_DIST_L2, the L2 norm of the distance
!    between A1 and A2.
!
  implicit none

  real ( kind = 8 ) a1(2)
  real ( kind = 8 ) a2(2)
  real ( kind = 8 ) r82vec_dist_l2

  r82vec_dist_l2 = sqrt ( sum ( ( a1(1:2) - a2(1:2) )**2 ) )

  return
end
subroutine random_initialize ( seed )

!*****************************************************************************80
!
!! RANDOM_INITIALIZE initializes the FORTRAN 90 random number seed.
!
!  Discussion:
!
!    If you don't initialize the random number generator, its behavior
!    is not specified.  If you initialize it simply by:
!
!      call random_seed ( )
!
!    its behavior is not specified.  On the DEC ALPHA, if that's all you
!    do, the same random number sequence is returned.  In order to actually
!    try to scramble up the random number generator a bit, this routine
!    goes through the tedious process of getting the size of the random
!    number seed, making up values based on the current time, and setting
!    the random number seed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer SEED.
!    If SEED is zero on input, then you're asking this routine to come up
!    with a seed value, which is returned as output.
!    If SEED is nonzero on input, then you're asking this routine to
!    use the input value of SEED to initialize the random number generator.
!
  implicit none

  integer count
  integer count_max
  integer count_rate
  integer i
  integer seed
  integer, allocatable :: seed_vector(:)
  integer seed_size
  real ( kind = 8 ) t
!
!  Initialize the random number seed.
!
  call random_seed ( )
!
!  Determine the size of the random number seed.
!
  call random_seed ( size = seed_size )
!
!  Allocate a seed of the right size.
!
  allocate ( seed_vector(seed_size) )

  if ( seed /= 0 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RANDOM_INITIALIZE'
    write ( *, '(a,i12)' ) '  Initialize RANDOM_NUMBER with user SEED = ', seed

  else

    call system_clock ( count, count_rate, count_max )

    seed = count

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RANDOM_INITIALIZE'
    write ( *, '(a,i12)' ) &
      '  Initialize RANDOM_NUMBER with arbitrary SEED = ', seed

  end if
!
!  Now set the seed.
!
  seed_vector(1:seed_size) = seed

  call random_seed ( put = seed_vector(1:seed_size) )
!
!  Free up the seed space.
!
  deallocate ( seed_vector )
!
!  Call the random number routine a bunch of times.
!
  do i = 1, 100
    call random_number ( harvest = t )
  end do

  return
end
subroutine region_plot_ppmb ( file_name, nrow, ncol, cell_num, dim_num, &
  box_min, box_max, cell_generator, cell_centroid, cell_weight )

!*****************************************************************************80
!
!! REGION_PLOT_PPMB makes a binary PPM plot of the CVT regions.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 December 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters
!
!    Input, character ( len = * ) FILE_NAME, the name of the file
!    to be created.
!
!    Input, integer NROW, NCOL, the number of pixels per row and
!    column of the image.
!
!    Input, integer CELL_NUM, the number of Voronoi cells.
!
!    Input, integer DIM_NUM, the dimension of the space.
!
!    Input, real ( kind = 8 ) BOX_MIN(DIM_NUM), BOX_MAX(DIM_NUM), the coordinates
!    of the two extreme corners of the bounding box.
!
!    Input, real ( kind = 8 ) CELL_GENERATOR(DIM_NUM,CELL_NUM), the Voronoi
!    cell generators.
!
!    Input, real ( kind = 8 ) CELL_CENTROID(CELL_NUM), the cell centroids.
!
!    Input, real ( kind = 8 ) CELL_WEIGHT(CELL_NUM), the cell weights.
!
  implicit none

  integer cell_num
  integer dim_num
  integer ncol
  integer nrow

  real ( kind = 8 ) angle
  integer b(nrow,ncol)
  real ( kind = 8 ) box_max(2)
  real ( kind = 8 ) box_min(2)
  real ( kind = 8 ) cell_b(cell_num)
  real ( kind = 8 ) cell_centroid(dim_num,cell_num)
  real ( kind = 8 ) cell_g(cell_num)
  real ( kind = 8 ) cell_generator(dim_num,cell_num)
  real ( kind = 8 ) cell_r(cell_num)
  real ( kind = 8 ) cell_weight(cell_num)
  real ( kind = 8 ) centroid_dist
  real ( kind = 8 ) r82vec_dist_l2
  character ( len = * ) file_name
  integer g(nrow,ncol)
  real ( kind = 8 ) generator_dist
  integer i
  integer ierror
  integer j
  logical, save :: loaded = .true.
  integer nearest
  real ( kind = 8 ) psize
  real ( kind = 8 ) pxsize
  real ( kind = 8 ) pysize
  integer r(nrow,ncol)
  real ( kind = 8 ) rb
  real ( kind = 8 ) rg
  real ( kind = 8 ) rr
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) xy(2)
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin

  if ( .not. loaded ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Colors:'
    write ( *, '(a)' ) ' '
    do i = 1, cell_num
      call i4_to_angle ( i-1, angle )
      call angle_to_rgb ( angle, rr, rg, rb )
      write ( *, '(i3,3f8.3)' ) i, rr, rg, rb
      cell_r(i) = rr
      cell_g(i) = rg
      cell_b(i) = rb
    end do

    loaded = .true.

  end if
!
!  Define the size of the page.
!
  xmin = box_min(1)
  ymin = box_min(2)

  xmax = box_max(1)
  ymax = box_max(2)

  pysize = ( ymax - ymin ) / real ( nrow - 1, kind = 8 )
  pxsize = ( xmax - xmin ) / real ( ncol - 1, kind = 8 )
  psize = min ( pysize, pxsize )

  do i = 1, nrow

    xy(2) = ( real ( nrow - i,     kind = 8 ) * ymax   &
            + real (        i - 1, kind = 8 ) * ymin ) &
            / real ( nrow     - 1, kind = 8 )

    do j = 1, ncol

      xy(1) = ( real ( ncol - j,     kind = 8 ) * xmin   &
              + real (        j - 1, kind = 8 ) * xmax ) &
              / real ( ncol     - 1, kind = 8 )

      call find_closest_size ( dim_num, xy, cell_num, cell_generator, &
        cell_weight, nearest )

      centroid_dist = r82vec_dist_l2 ( xy, cell_centroid(1:2,nearest) )

      generator_dist = r82vec_dist_l2 ( xy, cell_generator(1:2,nearest) )

      if ( centroid_dist <= 4.0D+00 * psize ) then
        rr = 0.0D+00
        rg = 0.0D+00
        rb = 0.0D+00
      else if ( generator_dist <= 7.0D+00 * psize ) then
        rr = 0.75D+00
        rg = 0.75D+00
        rb = 0.75D+00
      else

        call i4_to_angle ( nearest-1, angle )
        call angle_to_rgb ( angle, rr, rg, rb )

!       angle = 360.0D+00 * real ( nearest - 1, kind = 8 ) &
!                         / real ( cell_num, kind = 8 )
!       call hexcol ( angle, rr, rg, rb )
!       rr = cell_r(nearest)
!       rg = cell_g(nearest)
!       rb = cell_b(nearest)
      end if

      r(i,j) = int ( 255 * rr )
      g(i,j) = int ( 255 * rg )
      b(i,j) = int ( 255 * rb )

    end do

  end do
!
!  Now write all this data to a file.
!
  call ppmb_write ( file_name, ierror, nrow, ncol, r, g, b )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) 'PPMB_WRITE returns IERROR = ', ierror
  end if

  return
end
subroutine region_sampler ( dim_num, box_min, box_max, x, ngen )

!******************************************************************************
!
!! REGION_SAMPLER returns a sample point in the physical region.
!
!  Discussion:
!
!    The calculations are done in DIM_NUM dimensional space.
!
!    The physical region is enclosed in a bounding box.
!
!    A point is chosen in the bounding box, either by a uniform random
!    number generator, or from a vector Halton sequence.
!
!    If a user-supplied routine determines that this point is
!    within the physical region, this routine returns.  Otherwise,
!    a new random point is chosen.
!
!    The entries of the local vector HALTON_BASE should be distinct primes.
!    Right now, we're assuming DIM_NUM is no greater than 3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 September 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, real ( kind = 8 ) BOX_MIN(DIM_NUM), BOX_MAX(DIM_NUM), the
!    coordinates of the two extreme corners of the bounding box.
!
!    Output, real ( kind = 8 ) X(DIM_NUM), the random point.
!
!    Output, integer NGEN, the number of points that were generated.
!    This is at least 1, but may be larger if some points were rejected.
!
  integer dim_num

  real ( kind = 8 ) box_max(dim_num)
  real ( kind = 8 ) box_min(dim_num)
  integer ival
  integer ngen
  real ( kind = 8 ) r(dim_num)
  real ( kind = 8 ) x(dim_num)

  ngen = 0

  do

    ngen = ngen + 1

    if ( 10000 < ngen ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'REGION_SAMPLER - Fatal error!'
      write ( *, '(a,i8)' ) '  Number of rejected points = ', ngen
      write ( *, '(a)' ) '  There may be a problem with the geometry.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Current random value is:'
      write ( *, '(2x,3g14.6)' ) r(1:dim_num)
      write ( *, '(a)') ' '
      write ( *, '(a)' ) '  Current sample point is:'
      write ( *, '(2x,3g14.6)' ) x(1:dim_num)
      stop
    end if

    call random_number ( r(1:dim_num) )
!
!  Determine a point in the bounding box.
!
    x(1:dim_num) = ( ( 1.0D+00 - r(1:dim_num) ) * box_min(1:dim_num) &
                            + r(1:dim_num)   * box_max(1:dim_num) )

    call test_region ( x, dim_num, ival )

    if ( ival == 1 ) then
      exit
    end if

  end do

  return
end
subroutine s_to_i4vec ( s, n, ivec, reverse )

!*****************************************************************************80
!
!! S_TO_I4VEC converts an string of characters into an I4VEC.
!
!  Discussion:
!
!    This routine can be useful when trying to write character data to an
!    unformatted direct access file.
!
!    Depending on the internal byte ordering used on a particular machine,
!    the parameter REVERSE_ORDER may need to be set TRUE or FALSE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string of characters.
!    Each set of 4 characters is assumed to represent an integer.
!
!    Output, integer N, the number of integers read from the string.
!
!    Output, integer IVEC(*), an array of N integers which contains
!    the information from S.
!
!    Input, logical REVERSE, is TRUE if the bytes in a word need to be
!    reversed.
!
  implicit none

  integer from
  integer frompos
  integer ihi
  integer ilo
  integer ivec(*)
  integer j
  integer, parameter :: length = 8
  integer n
  integer nchar
  logical reverse
  character ( len = * ) s
  integer to
  integer topos

  nchar = len ( s )
  n = 0
  frompos = 0

  do ilo = 1, nchar, 4

    n = n + 1
    ihi = min ( ilo + 3, nchar )
    to = 0

    do j = ilo, ihi

      from = ichar ( s(j:j) )

      if ( reverse ) then
        topos = length * ( j - ilo )
      else
        topos = length * ( ilo + 3 - j )
      end if

      call mvbits ( from, frompos, length, to, topos )

    end do

    ivec(n) = to

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
