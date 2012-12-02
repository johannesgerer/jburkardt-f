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
!    The sextants are defined by the following points:
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
subroutine cvt_iteration ( dim_num, box_min, box_max, cell_num, &
  cell_generator, sample_num, cell_centroid, cell_volume, &
  cutoff_dist, region_volume )

!******************************************************************************
!
!! CVT_ITERATION takes one step of the CVT iteration.
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
!    22 January 2003
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
!    Input, integer CELL_NUM, the number of Voronoi cells.
!
!    Input, real ( kind = 8 ) CELL_GENERATOR(DIM_NUM,CELL_NUM), the Voronoi
!    cell generators. 
!
!    Input, integer SAMPLE_NUM, the number of sample points.
!
!    Output, real ( kind = 8 ) CELL_CENTROID(DIM_NUM,CELL_NUM), the Voronoi
!    cell centroids.
!
!    Output, real ( kind = 8 ) CELL_VOLUME(,CELL_NUM), the Voronoi cell volumes.
!
!    Input, real ( kind = 8 ) CUTOFF_DIST, the maximum influence distance.
!    If the nearest generator is further than this distance,
!    the point is not assigned.
!
!    Input, real ( kind = 8 ) REGION_VOLUME, the volume of the region that
!    is being divided by the CVT algorithm.
!
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
  real ( kind = 8 ) cutoff_dist
  logical, parameter :: debug = .false.
  integer j
  integer nearest(1)
  real ( kind = 8 ) nearest_dist
  integer ngen
  real ( kind = 8 ) region_volume
  real ( kind = 8 ) r8vec_dist_l2
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
    call find_closest ( dim_num, cell_num, 1, x, cell_generator, &
      nearest )
!
!  Add X to the averaging data for CELL_GENERATOR(*,NEAREST).
!
    if ( nearest(1) < 1 .or. cell_num < nearest(1) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CVT_CONSTRAINED_ITERATION - Fatal error!'
      write ( *, '(a,i6)' ) 'NEAREST = ', nearest(1)
      stop
    end if

    nearest_dist = r8vec_dist_l2 ( dim_num, x, &
      cell_generator(1:dim_num,nearest(1)) )

    if ( nearest_dist <= cutoff_dist ) then

      cell_centroid(1:dim_num,nearest(1)) = &
        cell_centroid(1:dim_num,nearest(1)) + x(1:dim_num)

      cell_count(nearest(1)) = cell_count(nearest(1)) + 1

    end if

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
subroutine find_closest ( ndim, n, sample_num, s, r, nearest )

!******************************************************************************
!
!! FIND_CLOSEST finds the nearest R point to each S point.
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
!    02 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NDIM, the spatial dimension.
!
!    Input, integer N, the number of cell generators.
!
!    Input, integer SAMPLE_NUM, the number of sample points.
!
!    Input, real ( kind = 8 ) S(NDIM,SAMPLE_NUM), the points to be checked.
!
!    Input, real ( kind = 8 ) R(NDIM,N), the cell generators.
!
!    Output, integer NEAREST(SAMPLE_NUM), the index of the nearest
!    cell generators.
!
  implicit none

  integer n
  integer ndim
  integer sample_num

  real ( kind = 8 ) dist_sq_min
  real ( kind = 8 ) dist_sq
  integer jr
  integer js
  integer nearest(sample_num)
  real ( kind = 8 ) r(ndim,n)
  real ( kind = 8 ) s(ndim,sample_num)

  do js = 1, sample_num

    dist_sq_min = huge ( dist_sq_min )
    nearest(js) = -1

    do jr = 1, n

      dist_sq = sum ( ( r(1:ndim,jr) - s(1:ndim,js) )**2 )

      if ( dist_sq < dist_sq_min ) then
        dist_sq_min = dist_sq
        nearest(js) = jr
      end if

    end do

  end do

  return
end
subroutine find_closest_size ( dim_num, cell_num, x, cell_generator, &
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
!    Input, integer CELL_NUM, the number of cell generatorrs.
!
!    Input, real ( kind = 8 ) X(DIM_NUM), the point to be checked.
!
!    Input, real ( kind = 8 ) CELL_GENERATOR(DIM_NUM,CELL_NUM), the cell
!    generators.
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
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 January 2003
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
!    Input, integer CELL_NUM, the number of Voronoi cells.
!
!    Output, real ( kind = 8 ) CELL_GENERATOR(DIM_NUM,CELL_NUM), the Voronoi
!    cell generators.
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
!      2**I_LOG_2(X) <= |I| < 2**(I_LOG_2(I)+1).
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
!    Output, real ( kind = 8 ) ANGLE, an angle, measured in degrees, between
!    0 and 360.
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
subroutine i4_to_rgb ( i, r, g, b )

!*****************************************************************************80
!
!! I4_TO_RGB maps integers to RGB colors.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 January 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer I, the index of the desired color.
!
!    Output, integer R, G, B, the RGB specifications for a color. 
!
  implicit none

  real ( kind = 8 ) angle
  integer b
  integer g
  integer i
  integer r
  real ( kind = 8 ) rb
  real ( kind = 8 ) rg
  real ( kind = 8 ) rr
!
!  Red
!
  if ( i == 1 ) then
    r = 255
    g = 0
    b = 0
!
!  Green
!
  else if ( i == 2 ) then
    r = 0
    g = 255
    b = 0
!
!  Blue
!
  else if ( i == 3 ) then
    r = 0
    g = 0
    b = 255
!
!  Cyan
!
  else if ( i == 4 ) then
    r = 0
    g = 255
    b = 255
!
!  Magenta
!
  else if ( i == 5 ) then
    r = 255
    g = 0
    b = 255
!
!  Yellow
!
  else if ( i == 6 ) then
    r = 255
    g = 255
    b = 0
!
!  Brown5
!
  else if ( i == 7 ) then
    r = 139
    g =  35
    b =  35
!
!  Orange
!
  else if ( i == 8 ) then
    r = 255
    g = 165
    b = 0
!
!  Goldenrod5
!
  else if ( i == 9 ) then
    r = 139
    g = 105
    b =  20
!
!  Medium Purple
!
  else if ( i == 10 ) then
    r = 147
    g = 112
    b = 219
!
!  Coral
!
   else if ( i == 11 ) then
    
    r = 255
    g = 127
    b =  80
!
!  Pink5
!
  else if ( i == 12 ) then
    r = 139
    g =  99
    b = 108
!
!  GreenYellow
!
  else if ( i == 13 ) then
    r = 173
    g = 255
    b =  47
!
!  Aquamarine
!
  else if ( i == 14 ) then
    r = 127
    g = 255
    b = 212
!
!  Pale Green3
!
  else if ( i == 15 ) then
    r = 124
    g = 205
    b = 124
!
!  Burlywood
!
  else if ( i == 16 ) then
    r = 222
    g = 184
    b = 135
!
!  Cornsilk3
!
  else if ( i == 17 ) then
    r = 205
    g = 200
    b = 177
!
!  Lemon_Chiffon3
!
  else if ( i == 18 ) then
    r = 205
    g = 201
    b = 165
!
!  Maroon
!
  else if ( i == 19 ) then
    r = 176
    g = 48
    b = 96
!
!  Slate_Blue2
!
  else if ( i == 20 ) then
    r = 131
    g = 111
    b = 255

  else

    call i4_to_angle ( i, angle )

    call angle_to_rgb ( angle, rr, rg, rb )

    r = min ( int ( rr * 255 ), 255 )
    g = min ( int ( rg * 255 ), 255 )
    b = min ( int ( rb * 255 ), 255 )

  end if

  return
end
function r8vec_dist_l2 ( n, a1, a2 )

!*****************************************************************************80
!
!! R8VEC_DIST_L2 returns the L2 distance between a pair of vectors.
!
!  Discussion:
!
!    The vector L2 norm is defined as:
!
!      sqrt ( sum ( 1 <= I <= N ) A(I)**2 ).
!
!    The vector L2 distance is:
!
!      norm ( A1(1:N) - A2(1:N) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 January 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the dimension of the vectors.
!
!    Input, real ( kind = 8 ) A1(N), A2(N), the vectors whose L2 distance is
!    desired.
!
!    Output, real ( kind = 8 ) r8VEC_DIST_L2, the L2 norm of the distance
!    between A1 and A2.
!
  implicit none

  integer  n

  real ( kind = 8 ) a1(n)
  real ( kind = 8 ) a2(n)
  real ( kind = 8 ) r8vec_dist_l2

  r8vec_dist_l2 = sqrt ( sum ( ( a1(1:n) - a2(1:n) )**2 ) )

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
  box_min, box_max, cell_generator, cell_centroid, cell_weight, cutoff_dist )

!*****************************************************************************80
!
!! REGION_PLOT_PPMB makes a binary PPM plot of the CVT regions.
!
!  Discussion:
!
!    This routine is only intended for the case DIM_NUM = 2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 January 2003
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
!    Input, real ( kind = 8 ) BOX_MIN(DIM_NUM), BOX_MAX(DIM_NUM), the
!    coordinates of the two extreme corners of the bounding box.
!
!    Input, real ( kind = 8 ) CELL_GENERATOR(DIM_NUM,CELL_NUM), the Voronoi
!    cell generators.
!
!    Input, real ( kind = 8 ) CELL_CENTROID(CELL_NUM), the cell centroids.
!
!    Input, real ( kind = 8 ) CELL_WEIGHT(CELL_NUM), the cell weights.
!
!    Input, real ( kind = 8 ) CUTOFF_DIST, the maximum influence distance.
!    If the nearest generator is further than this distance,
!    the point is not assigned.
!
  implicit none

  integer cell_num
  integer dim_num
  integer ncol
  integer nrow

  integer b(nrow,ncol)
  real ( kind = 8 ) box_max(2)
  real ( kind = 8 ) box_min(2)
  real ( kind = 8 ) cell_centroid(dim_num,cell_num)
  real ( kind = 8 ) cell_generator(dim_num,cell_num)
  real ( kind = 8 ) cell_weight(cell_num)
  real ( kind = 8 ) centroid_dist
  real ( kind = 8 ) cutoff_dist
  real ( kind = 8 ) r8vec_dist_l2
  character ( len = * ) file_name
  integer g(nrow,ncol)
  real ( kind = 8 ) generator_dist
  integer i
  integer ierror
  integer j
  integer nearest
  real ( kind = 8 ) psize
  real ( kind = 8 ) pxsize
  real ( kind = 8 ) pysize
  integer r(nrow,ncol)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) xy(2)
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin

  if ( dim_num /= 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'REGION_PLOT_PPMB - Fatal error!'
    write ( *, '(a)' ) '  This routine is only intended for cases'
    write ( *, '(a)' ) '  with spatial dimension 2.'
    stop
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

      call find_closest_size ( dim_num, cell_num, xy, cell_generator, &
        cell_weight, nearest )

      centroid_dist = r8vec_dist_l2 ( dim_num, xy, &
        cell_centroid(1:dim_num,nearest) )

      generator_dist = r8vec_dist_l2 ( dim_num, xy, &
        cell_generator(1:dim_num,nearest) )

      if ( cutoff_dist <= generator_dist ) then
        r(i,j) = 255
        g(i,j) = 255
        b(i,j) = 255
      else if ( centroid_dist <= 4.0D+00 * psize ) then
        r(i,j) = 0
        g(i,j) = 0
        b(i,j) = 0
      else if ( generator_dist <= 7.0D+00 * psize ) then
        r(i,j) = 192
        g(i,j) = 192
        b(i,j) = 192
      else
        call i4_to_rgb ( nearest, r(i,j), g(i,j), b(i,j) )
      end if

    end do

  end do
!
!  Now write all this data to a file.
!
  call ppmb_write ( file_name, ierror, nrow, ncol, r, g, b )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'REGION_PLOT_PPMB - Warning!'
    write ( *, '(a,i6)' ) '  PPMB_WRITE returns IERROR = ', ierror
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
      write ( *, '(3g14.6)' ) r(1:dim_num)
      write ( *, '(a)') ' '
      write ( *, '(a)' ) '  Current sample point is:'
      write ( *, '(3g14.6)' ) x(1:dim_num)
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
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    May 31 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 May 2001
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
  character ( len = 8 ) date
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer n
  integer s
  character ( len = 10 )  time
  integer values(8)
  integer y
  character ( len = 5 ) zone

  call date_and_time ( date, time, zone, values )

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

  write ( *, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
