subroutine bin_preprocess ( ndim, box_min, box_max, n, cell_generator, &
  nbin, bin_start, bin_last, bin_next )

!*****************************************************************************80
!
!! BIN_PREPROCESS organizes the preprocessing step for bins.
!
!  Discussion:
!
!    This routine is required to set up the bin data for use in the
!    nearest neighbor algorithm.
!
!    There are separate sets of calls for the 2D and 3D cases.  Although the
!    algorithms are essentially identical, the declarations of data
!    structures, in particular, of NBIN, were a little tricky to 
!    program generically.
!
!    To help with understanding the bin structure, the routine will
!    print out a few statistics about the bins on the first call only.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 April 2001
!
!  Author:
!
!    John Burkardt, Max Gunzburger, Janet Peterson
!
!  Parameters:
!
!    Input, integer NDIM, the spatial dimension.
!
!    Input, real ( kind = 8 ) BOX_MIN(NDIM), BOX_MAX(NDIM), the coordinates
!    of the two extreme corners of the bounding box.
!
!    Input, integer N, the number of Voronoi cells.
!
!    Input, real ( kind = 8 ) CELL_GENERATOR(NDIM,N), the Voronoi
!    cell generators.
!
!    Input, integer NBIN(3) is the number of bins to use in each direction.
!    For 3D problems, set NBIN(3) = 1.
!    For efficiency, these values should be set in such a way that the bins
!    are nearly square or cubical.
!
!    Output, integer BIN_START(NBIN(1),NBIN(2),NBIN(3)), the index of the first
!    cell generator in the bin, or -1 if none.
!
!    Output, integer BIN_LAST(NBIN(1),NBIN(2),NBIN(3)), the index of the last
!    cell generator in the bin, or -1 if none.
!
!    Output, integer BIN_NEXT(N), the index of the next cell generator in
!    the bin containing this cell generator.
!
  implicit none

  integer n
  integer ndim

  integer, dimension ( ndim ) :: nbin

  integer, dimension ( nbin(1),nbin(2),nbin(3) ) :: bin_last
  integer, dimension ( n ) :: bin_next
  integer, dimension ( nbin(1),nbin(2),nbin(3) ) :: bin_start 
  real ( kind = 8 ), dimension ( ndim ) :: box_min
  real ( kind = 8 ), dimension ( ndim ) :: box_max
  real ( kind = 8 ), dimension ( ndim, n ) :: cell_generator
  integer empty
  logical, save :: first_call = .true.
  integer i
  integer j
  integer k
  integer nonempty
  integer total

  if ( ndim == 2 ) then

    call r82vec_bin_even3 ( n, cell_generator, nbin, box_min, box_max, &
      bin_start, bin_last, bin_next )

    call r82vec_binned_reorder2 ( n, cell_generator, nbin, bin_start, &
      bin_last, bin_next )

    call r82vec_binned_sort_a2 ( n, cell_generator, nbin, bin_start, &
      bin_last )

  else if ( ndim == 3 ) then

    call r83vec_bin_even3 ( n, cell_generator, nbin, box_min, box_max, &
      bin_start, bin_last, bin_next )

    call r83vec_binned_reorder2 ( n, cell_generator, nbin, bin_start, &
      bin_last, bin_next )

    call r83vec_binned_sort_a2 ( n, cell_generator, nbin, bin_start, &
      bin_last )

  end if
!
!  On the first call only, prepare a brief statistical report.
!
  if ( first_call ) then

    first_call = .false.

    total = product ( nbin(1:ndim) )

    nonempty = 0
    empty = 0
 
    do i = 1, nbin(1)
      do j = 1, nbin(2)
        do k = 1, nbin(3)
          if ( 0 < bin_start(i,j,k) ) then
            nonempty = nonempty + 1
          else
            empty = empty + 1
          end if
        end do
      end do
    end do

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BIN_PREPROCESS:'
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Number of points =     ', n
    write ( *, '(a,i8)' ) '  Total number of bins = ', total
    write ( *, '(a,i8)' ) '  Number of empty bins = ', empty
    write ( *, '(a,i8)' ) '          nonempy bins = ', nonempty
    write ( *, '(a)' ) ' '
    write ( *, '(a,f5.1)' ) '  Percentage nonempty bins =         ', &
      real ( 100 * nonempty ) / real ( total )
    write ( *, '(a,g14.6)' ) '  Number of points per bin =         ', &
      real ( n ) / real ( total )
    write ( *, '(a,g14.6)' ) '  Number of points per nonempy bin = ', &
      real ( n ) / real ( nonempty )

  end if

  return
end
subroutine bin_to_r8_even2 ( nbin, bin, a, b, cmin, cmax )

!*****************************************************************************80
!
!! BIN_TO_R8_EVEN2 returns the limits for a given "bin" in [A,B].
!
!  Discussion:
!
!    The interval from A to B is divided into NBIN equal subintervals or bins.
!
!  Example:
!
!    NBIN = 5, A = 10, B = 20
!
!    BIN      CMIN  CMAX
!
!    1         10    12
!    2         12    14
!    3         14    16
!    4         16    18
!    5         18    20
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 February 2001
!
!  Author:
!
!    John Burkardt, Max Gunzburger, Janet Peterson
!
!  Parameters:
!
!    Input, integer NBIN, the number of bins.
!
!    Input, integer BIN, the index of the bin to be considered.
!    If BIN is less than 1, or greater than NBIN, the user will get what
!    the user deserves.
!
!    Input, real ( kind = 8 ) A, B, the lower and upper limits of the bin
!    interval.  While A is expected to be less than B, the code should
!    return useful results if A is actually greater than B.
!
!    Output, real ( kind = 8 ) CMIN, CMAX, the minimum and maximum limits
!    on the bin.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer bin
  real ( kind = 8 ) cmax
  real ( kind = 8 ) cmin
  integer nbin
!
!  Compute the bin limits.
!
  if ( bin < 1 ) then
    cmin = - huge ( cmin )
    cmax = a
  else if ( bin <= nbin ) then
    cmin = ( real ( nbin - bin + 1, kind = 8 ) * a &
           + real (        bin - 1, kind = 8 ) * b ) &
           / real ( nbin,           kind = 8 )
    cmax = ( real ( nbin - bin, kind = 8 ) * a &
           + real (        bin, kind = 8 ) * b ) &
           / real ( nbin,       kind = 8 )
  else if ( nbin < bin ) then
    cmin = b
    cmax = huge ( cmax )
  end if

  return
end
subroutine bin_to_r8vec_even3 ( ndim, nbin, bin, a, b, cmin, cmax )

!*****************************************************************************80
!
!! BIN_TO_R8VEC_EVEN3 returns the limits for a given R8VEC "bin" in [A,B].
!
!  Discussion:
!
!    The interval from A(I) to B(I) is divided into NBIN(I) equal
!    subintervals or bins.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 April 2001
!
!  Author:
!
!    John Burkardt, Max Gunzburger, Janet Peterson
!
!  Parameters:
!
!    Input, integer NDIM, the dimension of the space.
!
!    Input, integer NBIN(NDIM), the number of bins in each dimension.
!
!    Input, integer BIN(NDIM), the index of the bin to be considered.
!
!    Input, real ( kind = 8 ) A(NDIM), B(NDIM), the lower and upper limits
!    of the bin interval.  While A(I) is expected to be less than B(I), the
!    code should return useful results if A(I) is actually greater than B(I).
!
!    Output, real ( kind = 8 ) CMIN(NDIM), CMAX(NDIM), the minimum and maximum
!    limits on the bin.
!
  implicit none

  integer ndim

  real ( kind = 8 ) a(ndim)
  real ( kind = 8 ) b(ndim)
  integer bin(ndim)
  real ( kind = 8 ) cmax(ndim)
  real ( kind = 8 ) cmin(ndim)
  integer i
  integer nbin(ndim)

  do i = 1, ndim
    call bin_to_r8_even2 ( nbin(i), bin(i), a(i), b(i), cmin(i), cmax(i) )
  end do

  return
end
subroutine cvt_iteration ( ndim, box_min, box_max, n, cell_generator, ns_cvt, &
  use_diatom, use_bins, dr, updates, nbin, bin_start, bin_last, bin_next )

!*****************************************************************************80
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
!    generator.  Once the sampling is completed, the location of all the
!    generators is adjusted.  This step should decrease the discrepancy
!    between the generators and the centroids.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 April 2001
!
!  Author:
!
!    John Burkardt, Max Gunzburger, Janet Peterson
!
!  Reference:
!
!    John Burkardt, Max Gunzburger, Janet Peterson and Rebecca Brannon,
!    User Manual and Supporting Information for Library of Codes
!    for Centroidal Voronoi Placement and Associated Zeroth,
!    First, and Second Moment Determination,
!    Sandia National Laboratories Technical Report SAND2002-0099,
!    February 2002.
!
!  Parameters:
!
!    Input, integer NDIM, the spatial dimension.
!
!    Input, real ( kind = 8 ) BOX_MIN(NDIM), BOX_MAX(NDIM), the coordinates
!    of the two extreme corners of the bounding box.
!
!    Input, integer N, the number of Voronoi cells.
!
!    Input/output, real ( kind = 8 ) CELL_GENERATOR(NDIM,N), the Voronoi
!    cell generators.  On output, these have been modified
!
!    Input, integer NS_CVT, the number of sample points per generator.
!
!    Input, logical USE_DIATOM, is TRUE if DIATOM is to be called
!    to determine whether a point lies in the physical region; if it is
!    FALSE than a much simplified routine is used.
!
!    Input, logical USE_BINS, is TRUE if the bounding box is to be divided
!    up into bins to speed up the nearest neighbor search;
!    FALSE if the nearest neighbor seach is to be done naively.
!
!    Input, real ( kind = 8 ) DR, a tolerance used by DIATOM when testing
!    whether a point is within, outside of, or on the boundary of the
!    physical region.
!
!    Input/output, integer UPDATES(N), counts the number of times a cell
!    center has been updated.  Before the first call, all the entries of
!    UPDATES should be set to 1.  After each iteration, UPDATES will be
!    incremented by 1 for each cell generator that was updated.  Normally,
!    all of them will be so updated.
!
!    Input, integer NBIN(3) is the number of bins to use in each direction.
!    For 3D problems, set NBIN(3) = 1.
!    For efficiency, these values should be set in such a way that the bins
!    are nearly square or cubical.
!
!    Input, integer BIN_START(NBIN(1),NBIN(2),NBIN(3)), the index of the first
!    cell generator in the bin, or -1 if none.
!
!    Input, integer BIN_LAST(NBIN(1),NBIN(2),NBIN(3)), the index of the last
!    cell generator in the bin, or -1 if none.
!
!    Input, integer BIN_NEXT(N), the index of the next cell generator in
!    the bin containing this cell generator.
!
  implicit none

  integer n
  integer nbin(3)
  integer ndim

  integer, dimension ( nbin(1),nbin(2),nbin(3) ) :: bin_last
  integer, dimension ( n ) :: bin_next
  integer, dimension ( nbin(1),nbin(2),nbin(3) ) :: bin_start
  real ( kind = 8 ) box_max(1:ndim)
  real ( kind = 8 ) box_min(1:ndim)
  real ( kind = 8 ) cell_generator(ndim,n)
  real ( kind = 8 ) cell_generator2(ndim,n)
  integer count(n)
  logical, parameter :: debug = .false.
  real ( kind = 8 ) dr
  integer j
  integer nearest
  integer ngen
  integer ns_cvt
  integer, parameter :: random_sampler = 0
  logical reset
  integer updates(n)
  logical use_bins
  logical use_diatom
  real ( kind = 8 ) x(ndim)

  cell_generator2(1:ndim,1:n) = 0.0D+00
  count(1:n) = 0
  reset = .false.

  do j = 1, ns_cvt * n
!
!  Generate a sampling point X.
!
    call region_sampler ( ndim, box_min, box_max, dr, x, &
     random_sampler, reset, use_diatom, ngen )
!
!  Find the nearest cell generator.
!
    if ( use_bins ) then

      if ( ndim == 2 ) then

        call points_nearest_point_bins3_2d ( n, cell_generator, nbin, &
          box_min, box_max, bin_start, bin_last, bin_next, x, nearest )

      else if ( ndim == 3 ) then

        call points_nearest_point_bins3_3d ( n, cell_generator, nbin, &
          box_min, box_max, bin_start, bin_last, bin_next, x, nearest )

      end if

    else

      call find_closest ( ndim, x, n, cell_generator, nearest )

    end if
!
!  Add X to the averaging data for CELL_GENERATOR(*,NEAREST).
!
    cell_generator2(1:ndim,nearest) = &
      cell_generator2(1:ndim,nearest) + x(1:ndim)

    count(nearest) = count(nearest) + 1

  end do
!
!  Compute the new generators.
!
  do j = 1, n

    if ( count(j) /= 0 ) then

      cell_generator(1:ndim,j) = cell_generator2(1:ndim,j) &
        / real ( count(j), kind = 8 )
      updates(j) = updates(j) + 1

    end if

  end do

  return
end
subroutine find_closest ( ndim, x, n, cell_generator, nearest )

!*****************************************************************************80
!
!! FIND_CLOSEST finds the Voronoi cell generator closest to a point X.
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
!    18 January 2001
!
!  Author:
!
!    John Burkardt, Max Gunzburger, Janet Peterson
!
!  Reference:
!
!    John Burkardt, Max Gunzburger, Janet Peterson and Rebecca Brannon,
!    User Manual and Supporting Information for Library of Codes
!    for Centroidal Voronoi Placement and Associated Zeroth,
!    First, and Second Moment Determination,
!    Sandia National Laboratories Technical Report SAND2002-0099,
!    February 2002.
!
!  Parameters:
!
!    Input, integer NDIM, the spatial dimension.
!
!    Input, real ( kind = 8 ) X(NDIM), the point to be checked.
!
!    Input, integer N, the number of cell generators.
!
!    Input, real ( kind = 8 ) CELL_GENERATOR(NDIM,N), the cell generators.
!
!    Output, integer NEAREST, the index of the nearest cell generators.
!
  implicit none

  integer n
  integer ndim

  real ( kind = 8 ) cell_generator(ndim,n)
  real ( kind = 8 ) distance
  real ( kind = 8 ) dist_sq
  integer i
  integer nearest
  real ( kind = 8 ) x(ndim)

  nearest = 0
  distance = huge ( distance )

  do i = 1, n

    dist_sq = sum ( ( cell_generator(1:ndim,i) - x(1:ndim) )**2 )

    if ( dist_sq < distance ) then
      distance = dist_sq
      nearest = i
    end if

  end do

  distance = sqrt ( distance )

  return
end
subroutine generator_init ( ndim, box_min, box_max, n, cell_generator, &
  use_diatom, dr, random_generator )

!*****************************************************************************80
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
!    02 April 2001
!
!  Author:
!
!    John Burkardt, Max Gunzburger, Janet Peterson
!
!  Reference:
!
!    John Burkardt, Max Gunzburger, Janet Peterson and Rebecca Brannon,
!    User Manual and Supporting Information for Library of Codes
!    for Centroidal Voronoi Placement and Associated Zeroth,
!    First, and Second Moment Determination,
!    Sandia National Laboratories Technical Report SAND2002-0099,
!    February 2002.
!
!  Parameters:
!
!    Input, integer NDIM, the spatial dimension.
!
!    Input, real ( kind = 8 ) BOX_MIN(NDIM), BOX_MAX(NDIM), the coordinates
!    of the two extreme corners of the bounding box.
!
!    Input, integer N, the number of Voronoi cells.
!
!    Output, real ( kind = 8 ) CELL_GENERATOR(NDIM,N), the Voronoi cell
!    generators.
!
!    Input, logical USE_DIATOM, is TRUE if DIATOM is to be called
!    to determine whether a point lies in the physical region; if it is
!    FALSE than a much simplified routine is used.
!
!    Input, real ( kind = 8 ) DR, a tolerance used by DIATOM when testing
!    whether a point is within, outside of, or on the boundary of the
!    physical region.
!
!    Input, integer RANDOM_GENERATOR, specifies how the
!    Voronoi cell generators are to be initialized.
!    0, use the F90 RANDOM_NUMBER routine;
!    1, use the Halton sequence.
!
  implicit none

  integer n
  integer ndim

  real ( kind = 8 ) box_max(ndim)
  real ( kind = 8 ) box_min(ndim)
  real ( kind = 8 ) cell_generator(ndim,n)
  real ( kind = 8 ) dr
  integer i
  integer ngen
  integer random_generator
  logical reset
  logical use_diatom

  reset = .true.

  do i = 1, n
    call region_sampler ( ndim, box_min, box_max, dr, cell_generator(1,i), &
      random_generator, reset, use_diatom, ngen )
  end do

  return
end
subroutine i4_to_halton_vector ( seed, base, ndim, r )

!*****************************************************************************80
!
!! I4_TO_HALTON_VECTOR computes an element of a vector Halton sequence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 February 2001
!
!  Author:
!
!    John Burkardt, Max Gunzburger, Janet Peterson
!
!  Reference:
!
!    John Burkardt, Max Gunzburger, Janet Peterson and Rebecca Brannon,
!    User Manual and Supporting Information for Library of Codes
!    for Centroidal Voronoi Placement and Associated Zeroth,
!    First, and Second Moment Determination,
!    Sandia National Laboratories Technical Report SAND2002-0099,
!    February 2002.
!
!    John Halton,
!    On the efficiency of certain quasi-random sequences of points
!    in evaluating multi-dimensional integrals,
!    Numerische Mathematik,
!    Volume 2, pages 84-90.
!
!  Parameters:
!
!    Input, integer SEED, the index of the desired element.
!    Only the absolute value of SEED is considered.  SEED = 0 is allowed,
!    and returns R = 0.
!
!    Input, integer BASE(NDIM), the Halton bases, which should be
!    distinct prime numbers.  This routine only checks that each base
!    is greater than 1.
!
!    Input, integer NDIM, the dimension of the sequence.
!
!    Output, real ( kind = 8 ) R(NDIM), the SEED-th element of the Halton
!    sequence for the given bases.
!
  implicit none

  integer ndim

  integer base(ndim)
  real ( kind = 8 ) base_inv(ndim)
  integer digit(ndim)
  integer i
  real ( kind = 8 ) r(ndim)
  integer seed
  integer seed2(ndim)

  seed2(1:ndim) = abs ( seed )

  r(1:ndim) = 0.0D+00

  if ( any ( base(1:ndim) <= 1 ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_TO_HALTON_VECTOR - Fatal error!'
    write ( *, '(a)' ) '  An input base BASE is <= 1!'
    do i = 1, ndim
      write ( *, '(2i8)' ) i, base(i)
    end do
    stop
  end if

  base_inv(1:ndim) = 1.0D+00 / real ( base(1:ndim), kind = 8 )

  do while ( any ( seed2(1:ndim) /= 0 ) )
    digit(1:ndim) = mod ( seed2(1:ndim), base(1:ndim) )
    r(1:ndim) = r(1:ndim) + real ( digit(1:ndim), kind = 8 ) * base_inv(1:ndim)
    base_inv(1:ndim) = base_inv(1:ndim) / real ( base(1:ndim), kind = 8 )
    seed2(1:ndim) = seed2(1:ndim) / base(1:ndim)
  end do

  return
end
subroutine index_box2_next_2d ( n1, n2, ic, jc, i, j, more )

!*****************************************************************************80
!
!! INDEX_BOX2_NEXT_2D produces indices on the surface of a box in 2D.
!
!  Discussion:
!
!    The box has center at (IC,JC), and has half-widths N1 and N2.
!    The indices are exactly those which are between (IC-N1,JC-N2) and
!    (IC+N1,JC+N2) with the property that at least one of I and J
!    is an "extreme" value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 May 2001
!
!  Author:
!
!    John Burkardt, Max Gunzburger, Janet Peterson
!
!  Parameters:
!
!    Input, integer N1, N2, the half-widths of the box, that is, the
!    maximum distance allowed between (IC,JC) and (I,J).
!
!    Input, integer IC, JC, the central cell of the box.
!
!    Input/output, integer I, J.  On input, the previous index set.
!    On output, the next index set.  On the first call, MORE should
!    be set to FALSE, and the input values of I and J are ignored.
!
!    Input/output, logical MORE.
!    On the first call for a given box, the user should set MORE to FALSE.
!    On return, the routine sets MORE to TRUE.
!    When there are no more indices, the routine sets MORE to FALSE.
!
  implicit none

  integer i
  integer ic
  integer j
  integer jc
  logical more
  integer n1
  integer n2

  if ( .not. more ) then
    more = .true.
    i = ic - n1
    j = jc - n2
    return
  end if

  if ( i == ic + n1 .and. j == jc + n2 ) then
    more = .false.
    return
  end if
!
!  Increment J.
!
  j = j + 1
!
!  Check J.
!
  if ( jc + n2 < j ) then
    j = jc - n2
    i = i + 1
  else if ( j < jc + n2 .and. ( i == ic - n1 .or. i == ic + n1 ) ) then
    return
  else
    j = jc + n2
    return
  end if

  return
end
subroutine index_box2_next_3d ( n1, n2, n3, ic, jc, kc, i, j, k, more )

!*****************************************************************************80
!
!! INDEX_BOX2_NEXT_3D produces indices on the surface of a box in 3D.
!
!  Discussion:
!
!    The box has a central cell of (IC,JC,KC), with a half widths of
!    (N1,N2,N3).  The indices are exactly those between (IC-N1,JC-N2,KC-N3)
!    and (IC+N1,JC+N2,KC+N3) with the property that at least one of I, J,
!    and K is an "extreme" value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 May 2001
!
!  Author:
!
!    John Burkardt, Max Gunzburger, Janet Peterson
!
!  Parameters:
!
!    Input, integer N1, N2, N3, the "half widths" of the box, that is, the
!    maximum distances from the central cell allowed for I, J and K.
!
!    Input, integer IC, JC, KC, the central cell of the box.
!
!    Input/output, integer I, J, K.  On input, the previous index set.
!    On output, the next index set.  On the first call, MORE should
!    be set to FALSE, and the input values of I, J, and K are ignored.
!
!    Input/output, logical MORE.
!    On the first call for a given box, the user should set MORE to FALSE.
!    On return, the routine sets MORE to TRUE.
!    When there are no more indices, the routine sets MORE to FALSE.
!
  implicit none

  integer i
  integer ic
  integer j
  integer jc
  integer k
  integer kc
  logical more
  integer n1
  integer n2
  integer n3

  if ( .not. more ) then
    more = .true.
    i = ic - n1
    j = jc - n2
    k = kc - n3
    return
  end if

  if ( i == ic + n1 .and. j == jc + n2 .and. k == kc + n3 ) then
    more = .false.
    return
  end if
!
!  Increment K.
!
  k = k + 1
!
!  Check K.
!
  if ( kc + n3 < k ) then
    k = kc - n3
    j = j + 1
  else if ( k < kc + n3 .and. &
    ( i == ic - n1 .or. i == ic + n1 .or. &
      j == jc - n2 .or. j == jc + n2 ) ) then
    return
  else
    k = kc + n3
    return
  end if
!
!  Check J.
!
  if ( jc + n2 < j ) then
    j = jc - n2
    i = i + 1
  else if ( j < jc + n2 .and. &
    ( i == ic - n1 .or. i == ic + n1 .or. &
      k == kc - n3 .or. k == kc + n3 ) ) then
    return
  else
    j = jc + n2
    return
  end if

  return
end
subroutine points_nearest_point_bins3_2d ( nset, pset, nbin, bin_min, bin_max, &
  bin_start, bin_last, bin_next, ptest, i_min )

!*****************************************************************************80
!
!! POINTS_NEAREST_POINT_BINS3_2D finds the nearest point to a given point in 2D.
!
!  Discussion:
!
!    This code differs from POINTS_NEAREST_POINTS_BINS_2D by allowing the
!    user to specify the number of bins in each dimension.
!
!    A set of NSET points with coordinates PSET is assumed to lie within a
!    rectangle.  The limits of this rectangle are given in BIN_MIN and BIN_MAX.
!    The rectangle is divided up into an NBIN(1) by NBIN(2) regular grid of
!    cells.
!
!    The cells are ordered lexicographically, as suggested by the following
!    diagram for NBIN = (/ 5, 4 /)
!
!     *----*----*----*----*----*
!     |    |    |    |    |    |
!     | 16 | 17 | 18 | 19 | 20 |
!     |    |    |    |    |    |
!     *----*----*----*----*----*
!     |    |    |    |    |    |
!     | 11 | 12 | 13 | 14 | 15 |
!     |    |    |    |    |    |
!     *----*----*----*----*----*
!     |    |    |    |    |    |
!     |  6 |  7 |  8 |  9 | 10 |
!     |    |    |    |    |    |
!     *----*----*----*----*----*
!     |    |    |    |    |    |
!     |  1 |  2 |  3 |  4 |  5 |
!     |    |    |    |    |    |
!     *----*----*----*----*----*
!
!    The points in PSET are ordered by cell, and within a cell they
!    are ordered lexicographically.  Thus, for points P1 and P2,
!
!      P1 < P2 implies that
!      * P1 is in a lower ordered cell than P2, or
!      * P1 is in the same cell as P2, but P1.X < P2.X, or
!      * P1 is in the same cell as P2, P1.X = P2.X, but P1.Y < P2.Y.
!
!    The arrays BIN_START(I,J) and BIN_LAST(I,J) are given the coordinates
!    I, J of a cell, and return the lowest and highest index in PSET of a
!    point in the I, J cell.  All indices in between also belong to
!    this cell.  If the cell has no points, then both arrays are -1.
!
!
!    After all this preprocessing, the algorithm for finding the nearest
!    point goes as follows:
!
!    1) for a point PTEST, determine the cell it belongs to.
!       Note that this algorithm will NOT be valid if PTEST lies outside
!       the bin limits.
!
!    2) Search for a cell that has at least one member of PSET in it.
!       We start at the cell containing PTEST, but if there are no members
!       there, we spiral out from the cell, one layer at a time.
!
!    3) Within this cell, find the point nearest to PTEST.  We now know that
!       we don't need to search any cell whose points will all be further
!       from PTEST than this distance.
!
!    4) Now, search in all other cells that could have a point closer to
!       PTEST than what we have found so far.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 July 2008
!
!  Author:
!
!    John Burkardt, Max Gunzburger, Janet Peterson
!
!  Reference:
!
!    Jon Bentley, Bruce Weide, Andrew Yao,
!    Optimal Expected Time Algorithms for Closest Point Problems,
!    ACM Transactions on Mathematical Software,
!    Volume 6, Number 4, December 1980, pages 563-580.
!
!  Parameters:
!
!    Input, integer NSET, the number of points in the set.
!
!    Input, real ( kind = 8 ) PSET(2,NSET), the coordinates of the points
!    in the set.
!
!    Input, integer NBIN(2), the number of cells in the horizontal and
!    vertical directions.
!
!    Input, real ( kind = 8 ) BIN_MIN(2), BIN_MAX(2), the minimum and maximum
!    bin values.
!
!    Input, integer BIN_START(NBIN(1),NBIN(2)), BIN_LAST(NBIN(1),NBIN(2)),
!    indicates the index of the first and last element in the bin, or -1
!    if none.
!
!    Input, integer BIN_NEXT(NSET), the index of the next element of the bin
!    containing this element.
!
!    Input, real ( kind = 8 ) PTEST(2), the coordinates of the test point.
!
!    Output, integer I_MIN, the index of the nearest point in PSET to PTEST.
!
  implicit none

  integer, parameter :: ndim = 2

  integer nbin(ndim)
  integer nset

  integer bin(ndim)
  integer bin_last(nbin(1),nbin(2))
  real ( kind = 8 ) bin_max(ndim)
  real ( kind = 8 ) bin_min(ndim)
  integer bin_start(nbin(1),nbin(2))
  integer bin_next(nset)
  real ( kind = 8 ) c_max(ndim)
  real ( kind = 8 ) c_min(ndim)
  real ( kind = 8 ) d_min
  real ( kind = 8 ) d_min_sq
  real ( kind = 8 ) d_sq
  integer i
  integer i_min
  integer ic
  integer j
  integer jc
  integer layer
  real ( kind = 8 ) layer_width
  logical more_bins
  integer node
  real ( kind = 8 ) pset(ndim,nset)
  real ( kind = 8 ) ptest(ndim)
  real ( kind = 8 ) search_radius
!
!  Special cases.
!
  if ( nset <= 0 ) then
    d_min = huge ( d_min )
    i_min = 0
    return
  end if

  if ( nset == 1 ) then
    d_min = sqrt ( sum ( ( ptest(1:ndim) - pset(1:ndim,1) )**2 ) )
    i_min = 1
    return
  end if
!
!  The efficiency of the code will suffer if the data in the vector
!
!    ( bin_max(1:ndim) - bin_min(1:ndim) ) / real ( nbin(1:ndim), kind = 8 )
!
!  varies significantly.
!
  layer_width = minval ( &
    ( bin_max(1:ndim) - bin_min(1:ndim) ) / real ( nbin(1:ndim), kind = 8 ) &
  )

  d_min_sq = huge ( d_min_sq )
  i_min = 0
  search_radius = 0.0D+00
!
!  Determine the bin coordinates of the point P.
!
  call r8vec_to_bin_even3 ( ndim, nbin, bin_min, bin_max, ptest, bin )
!
!  Determine the radius of the ball of space that will be completely
!  searched after this first bin, "layer 0", is done.
!
  call bin_to_r8vec_even3 ( ndim, nbin, bin, bin_min, bin_max, c_min, c_max )
!
!  Set the central bin of the layers.
!
  ic = bin(1)
  jc = bin(2)

  layer = 0
!
!  Search all legal bins in layer LAYER.
!
  do

    more_bins = .false.
    call index_box2_next_2d ( layer, layer, ic, jc, i, j, more_bins )
!
!  In layer LAYER, search each BIN I, J.
!
    do

      if ( 1 <= i .and. i <= nbin(1) .and. 1 <= j .and. j <= nbin(2) ) then

        node = bin_start(i,j)

        do while ( 0 < node )

          d_sq = sum ( ( ptest(1:ndim) - pset(1:ndim,node) )**2 )

          if ( d_sq < d_min_sq ) then
            d_min_sq = d_sq
            i_min = node
          end if

          node = bin_next(node)

        end do

      end if
!
!  We have now searched all points in bin I, J.
!
!  Determine the next bin in this layer.
!
!  But if it lies outside the region, discard it, and go on to the next one.
!  Once you get past the last bin, exit because you're done the layer.
!
      do

        call index_box2_next_2d ( layer, layer, ic, jc, i, j, more_bins )

        if ( .not. more_bins ) then
          exit
        end if

        if ( 1 <= i .and. i <= nbin(1) .and. &
             1 <= j .and. j <= nbin(2) ) then
          exit
        end if

      end do

      if ( .not. more_bins ) then
        exit
      end if

    end do
!
!  We've completed layer LAYER.
!  Update the radius of the searched area.
!
    if ( layer == 0 ) then
      search_radius = min ( &
        minval ( abs ( ptest(1:ndim) - c_min(1:ndim) ) ), &
        minval ( abs ( ptest(1:ndim) - c_max(1:ndim) ) ) )
    else
      search_radius = search_radius + layer_width
    end if
!
!  We are done if:
!
!    * We've found at least one neighbor;
!    AND
!    * We've searched all legal boxes that contain the circle
!      with PTEST at the center and the nearest neighbor on the circumference.
!
    if ( i_min /= 0 ) then
      d_min = sqrt ( d_min_sq )
      if ( d_min <= search_radius ) then
        exit
      end if
    end if

    layer = layer + 1

  end do

  return
end
subroutine points_nearest_point_bins3_3d ( nset, pset, nbin, bin_min, bin_max, &
  bin_start, bin_last, bin_next, ptest, i_min )

!*****************************************************************************80
!
!! POINTS_NEAREST_POINT_BINS3_3D finds the nearest point to a given point in 3D.
!
!  Discussion:
!
!    This code differs from POINTS_NEAREST_POINTS_BINS_3D by allowing the
!    user to specify the number of bins in each dimension.
!
!    A set of NSET points with coordinates PSET is assumed to lie within a
!    box.  The limits of this rectangle are given in BIN_MIN and BIN_MAX.
!    The rectangle is divided up into an NBIN(1) by NBIN(2) by NBIN(3)
!    regular grid of cells.
!
!    The cells are ordered lexicographically, as suggested by the following
!    diagram for NBIN = (/ 5, 4, 2 /)
!
!             Z LAYER 1                       Z LAYER 2
!     *----*----*----*----*----*     *----*----*----*----*----*
!     |    |    |    |    |    |     |    |    |    |    |    |
!     | 16 | 17 | 18 | 19 | 20 |     | 36 | 37 | 38 | 39 | 40 |
!     |    |    |    |    |    |     |    |    |    |    |    |
!     *----*----*----*----*----*     *----*----*----*----*----*
!     |    |    |    |    |    |     |    |    |    |    |    |
!     | 11 | 12 | 13 | 14 | 15 |     | 31 | 32 | 33 | 34 | 35 |
!     |    |    |    |    |    |     |    |    |    |    |    |
!     *----*----*----*----*----*     *----*----*----*----*----*
!     |    |    |    |    |    |     |    |    |    |    |    |
!     |  6 |  7 |  8 |  9 | 10 |     | 26 | 27 | 28 | 29 | 30 |
!     |    |    |    |    |    |     |    |    |    |    |    |
!     *----*----*----*----*----*     *----*----*----*----*----*
!     |    |    |    |    |    |     |    |    |    |    |    |
!     |  1 |  2 |  3 |  4 |  5 |     | 21 | 22 | 23 | 24 | 25 |
!     |    |    |    |    |    |     |    |    |    |    |    |
!     *----*----*----*----*----*     *----*----*----*----*----*
!
!    The points in PSET are ordered by cell, and within a cell they
!    are ordered lexicographically.  Thus, for points P1 and P2,
!
!      P1 < P2 implies that
!      * P1 is in a lower ordered cell than P2, or
!      * P1 is in the same cell as P2, but P1.X < P2.X, or
!      * P1 is in the same cell as P2, P1.X = P2.X, but P1.Y < P2.Y.
!      * P1 is in the same cell as P2, P1.X = P2.X, P1.Y = P2.Y,
!        but P1.Z < P2.Z
!
!    The arrays BIN_START(I,J,K) and BIN_LAST(I,J,K) are given the coordinates
!    I, J, K of a cell, and return the lowest and highest index in PSET of a
!    point in the I, J, K cell.  All indices in between also belong to
!    this cell.  If the cell has no points, then both arrays are -1.
!
!
!    After all this preprocessing, the algorithm for finding the nearest
!    point goes as follows:
!
!    1) for a point PTEST, determine the cell it belongs to.
!       Note that this algorithm will NOT be valid if PTEST lies outside
!       the bin limits.
!
!    2) Search for a cell that has at least one member of PSET in it.
!       We start at the cell containing PTEST, but if there are no members
!       there, we spiral out from the cell, one layer at a time.
!
!    3) Within this cell, find the point nearest to PTEST.  We now know that
!       we don't need to search any cell whose points will all be further
!       from PTEST than this distance.
!
!    4) Now, search in all other cells that could have a point closer to
!       PTEST than what we have found so far.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 July 2008
!
!  Author:
!
!    John Burkardt, Max Gunzburger, Janet Peterson
!
!  Reference:
!
!    Jon Bentley, Bruce Weide, Andrew Yao,
!    Optimal Expected Time Algorithms for Closest Point Problems,
!    ACM Transactions on Mathematical Software,
!    Volume 6, Number 4, December 1980, pages 563-580.
!
!  Parameters:
!
!    Input, integer NSET, the number of points in the set.
!
!    Input, real ( kind = 8 ) PSET(3,NSET), the coordinates of the points
!    in the set.
!
!    Input, integer NBIN(3), the number of cells in the X, Y and Z directions.
!
!    Input, real ( kind = 8 ) BIN_MIN(3), BIN_MAX(3), the minimum and
!    maximum bin values.
!
!    Input, integer BIN_START(NBIN(1),NBIN(2),NBIN(3)),
!    BIN_LAST(NBIN(1),NBIN(2),NBIN(3)), the index of the first and last
!    element in the bin, or -1 if none.
!
!    Input, integer BIN_NEXT(NSET), the index of the next element of the bin
!    containing this element.
!
!    Input, real ( kind = 8 ) PTEST(3), the coordinates of the test points.
!
!    Output, integer I_MIN, the index of the nearest point in PSET to PTEST.
!
  implicit none

  integer, parameter :: ndim = 3

  integer nbin(ndim)
  integer nset

  integer bin(ndim)
  integer bin_last(nbin(1),nbin(2),nbin(3))
  real ( kind = 8 ) bin_max(ndim)
  real ( kind = 8 ) bin_min(ndim)
  integer bin_start(nbin(1),nbin(2),nbin(3))
  integer bin_next(nset)
  real ( kind = 8 ) c_max(ndim)
  real ( kind = 8 ) c_min(ndim)
  real ( kind = 8 ) d_min
  real ( kind = 8 ) d_min_sq
  real ( kind = 8 ) d_sq
  integer i
  integer i_min
  integer ic
  integer j
  integer jc
  integer k
  integer kc
  integer layer
  real ( kind = 8 ) layer_width
  logical more_bins
  integer node
  real ( kind = 8 ) pset(ndim,nset)
  real ( kind = 8 ) ptest(ndim)
  real ( kind = 8 ) search_radius
!
!  Special cases.
!
  if ( nset <= 0 ) then
    d_min = huge ( d_min )
    i_min = 0
    return
  end if

  if ( nset == 1 ) then
    d_min = sqrt ( sum ( ( ptest(1:ndim) - pset(1:ndim,1) )**2 ) )
    i_min = 1
    return
  end if
!
!  The efficiency of the code will suffer if the data in the vector
!
!    ( bin_max(1:ndim) - bin_min(1:ndim) ) / real ( nbin(1:ndim), kind = 8 )
!
!  varies significantly.
!
  layer_width = minval ( &
    ( bin_max(1:ndim) - bin_min(1:ndim) ) / real ( nbin(1:ndim), kind = 8 ) &
  )

  d_min_sq = huge ( d_min_sq )
  i_min = 0
  search_radius = 0.0D+00
!
!  Determine the bin coordinates of the point P.
!
  call r8vec_to_bin_even3 ( ndim, nbin, bin_min, bin_max, ptest(1), bin )
!
!  Determine the radius of the ball of space that will be completely
!  searched after this first bin, "layer 0", is done.
!
  call bin_to_r8vec_even3 ( ndim, nbin, bin, bin_min, bin_max, c_min, c_max )
!
!  Set the central bin of the layers.
!
  ic = bin(1)
  jc = bin(2)
  kc = bin(3)

  layer = 0
!
!  Search all legal bins in layer LAYER.
!
  do

    more_bins = .false.
    call index_box2_next_3d ( layer, layer, layer, ic, jc, kc, i, j, k, &
      more_bins )
!
!  In layer LAYER, search each BIN I, J, K.
!
    do

      if ( 1 <= i .and. i <= nbin(1) .and. &
           1 <= j .and. j <= nbin(2) .and. &
           1 <= k .and. k <= nbin(3) ) then

        node = bin_start(i,j,k)

        do while ( 0 < node )

          d_sq = sum ( ( ptest(1:ndim) - pset(1:ndim,node) )**2 )

          if ( d_sq < d_min_sq ) then
            d_min_sq = d_sq
            i_min = node
          end if

          node = bin_next(node)

        end do

      end if
!
!  We have now searched all points in bin I, J, K.
!
!  Determine the next bin in this layer.
!
!  But if it lies outside the region, discard it, and go on to the next one.
!  Once you get past the last bin, exit because you're done the layer.
!
      do

        call index_box2_next_3d ( layer, layer, layer, ic, jc, kc, i, j, k, &
          more_bins )

        if ( .not. more_bins ) then
          exit
        end if

        if ( 1 <= i .and. i <= nbin(1) .and. &
             1 <= j .and. j <= nbin(2) .and. &
             1 <= k .and. k <= nbin(3) ) then
          exit
        end if

      end do

      if ( .not. more_bins ) then
        exit
      end if

    end do
!
!  We've completed layer LAYER.
!  Update the radius of the searched area.
!
    if ( layer == 0 ) then
      search_radius = min ( &
        minval ( abs ( ptest(1:ndim) - c_min(1:ndim) ) ), &
        minval ( abs ( ptest(1:ndim) - c_max(1:ndim) ) ) )
    else
      search_radius = search_radius + layer_width
    end if
!
!  We are done with PTEST if:
!
!    * We've found at least one neighbor;
!    AND
!    * We've searched all legal boxes that contain the circle
!      with PTEST at the center and the nearest neighbor on the circumference.
!
    if ( i_min /= 0 ) then
      d_min = sqrt ( d_min_sq )
      if ( d_min <= search_radius ) then
        exit
      end if
    end if

    layer = layer + 1

  end do
!
!  We are now done with all the layers.
!
  return
end
subroutine quality ( ndim, n, cell_moment, cell_volume, region_volume )

!*****************************************************************************80
!
!! QUALITY computes some quality measures for a set of points in a region.
!
!  Discussion:
!
!    The quality measures report on how evenly spread the points are.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 May 2001
!
!  Author:
!
!    John Burkardt, Max Gunzburger, Janet Peterson
!
!  Reference:
!
!    John Burkardt, Max Gunzburger, Janet Peterson and Rebecca Brannon,
!    User Manual and Supporting Information for Library of Codes
!    for Centroidal Voronoi Placement and Associated Zeroth,
!    First, and Second Moment Determination,
!    Sandia National Laboratories Technical Report SAND2002-0099,
!    February 2002.
!
!  Parameters:
!
!    Input, integer NDIM, the spatial dimension.
!
!    Input, integer N, the number of cell generators.
!
!    Input, real ( kind = 8 ) CELL_MOMENT(NDIM,NDIM,N), the second moment
!    matrix for each Voronoi cell.
!
!    Input, real ( kind = 8 ) CELL_VOLUME(N), the Voronoi cell volumes.
!
!    Input, real ( kind = 8 ) REGION_VOLUME, the volume of the region,
!    as input by the user or estimated by VCM.
!
  implicit none

  integer n
  integer ndim

  real ( kind = 8 ) cell_det(n)
  real ( kind = 8 ) cell_moment(ndim,ndim,n)
  real ( kind = 8 ) cell_trace(n)
  real ( kind = 8 ) cell_volume(n)
  real ( kind = 8 ) dd_l1
  real ( kind = 8 ) dd_l2
  real ( kind = 8 ) dd_linf
  real ( kind = 8 ) r8mat_det_2d
  real ( kind = 8 ) r8mat_det_3d
  real ( kind = 8 ) ev
  real ( kind = 8 ) ev_l1
  real ( kind = 8 ) ev_l2
  real ( kind = 8 ) ev_linf
  integer i
  integer j
  real ( kind = 8 ) matrix2(2,2)
  real ( kind = 8 ) matrix3(3,3)
  real ( kind = 8 ) region_volume
  real ( kind = 8 ) tr
  real ( kind = 8 ) tr_l1
  real ( kind = 8 ) tr_l2
  real ( kind = 8 ) tr_linf
!
!  Measure 1: the deviation of the cell volumes from the expected cell volume.
!
  ev = region_volume / real ( n, kind = 8 )
  ev_linf = maxval ( abs ( ev - cell_volume(1:n) ) )
  ev_l1 = sum ( abs ( ev - cell_volume(1:n) ) )
  ev_l2 = sqrt ( sum ( ( ev - cell_volume(1:n) )**2 ) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'QUALITY'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Measure #1:'
  write ( *, '(a)' ) '    ( Cell_Volume - Expected Cell Volume )'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Expected Cell Volume = ', ev
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  L1 norm =              ', ev_l1
  write ( *, '(a,g14.6)' ) '  L2 norm =              ', ev_l2
  write ( *, '(a,g14.6)' ) '  L1 norm / N =          ', ev_l1 &
    / real ( n, kind = 8 )
  write ( *, '(a,g14.6)' ) '  L2 norm / sqrt ( N ) = ', ev_l2 &
    / sqrt ( real ( n, kind = 8 ) )
  write ( *, '(a,g14.6)' ) '  L-Inf norm =           ', ev_linf
!
!  Measure 2: the deviation of the traces of the cell second moment matrices
!  from the average.
!
  do i = 1, n
    cell_trace(i) = 0.0D+00
    do j = 1, ndim
      cell_trace(i) = cell_trace(i) + cell_moment(j,j,i)
    end do
  end do

  tr = sum ( cell_trace(1:n) ) / real ( n, kind = 8 )
  tr_linf = maxval ( abs ( tr - cell_trace(1:n) ) )
  tr_l1 = sum ( abs ( tr - cell_trace(1:n) ) )
  tr_l2 = sqrt ( sum ( ( tr - cell_trace(1:n) )**2 ) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Measure #2:'
  write ( *, '(a)' ) '    ( Cell_Trace - Average Cell Trace )'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Average Cell Trace = ', tr
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  L1 norm =              ', tr_l1
  write ( *, '(a,g14.6)' ) '  L2 norm =              ', tr_l2
  write ( *, '(a,g14.6)' ) '  L1 norm / N =          ', tr_l1 &
    / real ( n, kind = 8 )
  write ( *, '(a,g14.6)' ) '  L2 norm / sqrt ( N ) = ', tr_l2 &
    / sqrt ( real ( n, kind = 8 ) )
  write ( *, '(a,g14.6)' ) '  L-Inf norm =           ', tr_linf
!
!  Measure 3: the determinant of the deviatoric matrix
!
  if ( ndim == 2 ) then

    do i = 1, n
      matrix2(1:ndim,1:ndim) = cell_moment(1:ndim,1:ndim,i)
      do j = 1, ndim
        matrix2(j,j) = matrix2(j,j) - cell_trace(i) / real ( ndim, kind = 8 )
      end do

      cell_det(i) = r8mat_det_2d ( matrix2 )

    end do

  else if ( ndim == 3 ) then

    do i = 1, n
      matrix3(1:ndim,1:ndim) = cell_moment(1:ndim,1:ndim,i)
      do j = 1, ndim
        matrix3(j,j) = matrix3(j,j) - cell_trace(i) / real ( ndim, kind = 8 )
      end do

      cell_det(i) = r8mat_det_3d ( matrix3 )

    end do

  end if

  dd_linf = maxval ( abs ( cell_det(1:n) ) )
  dd_l1 = sum ( abs ( cell_det(1:n) ) )
  dd_l2 = sqrt ( sum ( ( cell_det(1:n) )**2 ) )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Measure #3:'
  write ( *, '(a)' ) '    ( The determinant of the deviatoric matrix )'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  L1 norm = 	     ', dd_l1
  write ( *, '(a,g14.6)' ) '  L2 norm =              ', dd_l2
  write ( *, '(a,g14.6)' ) '  L1 norm / N =          ', dd_l1 &
    / real ( n, kind = 8 )
  write ( *, '(a,g14.6)' ) '  L2 norm / sqrt ( N ) = ', dd_l2 &
    / sqrt ( real ( n, kind = 8 ) )
  write ( *, '(a,g14.6)' ) '  L-Inf norm =           ', dd_linf

  return
end
subroutine r8_to_bin_even2 ( nbin, a, b, c, bin )

!*****************************************************************************80
!
!! R8_TO_BIN_EVEN2 determines the appropriate "bin" for C in [A,B].
!
!  Discussion:
!
!    The interval from A to B is divided into NBIN equal subintervals or bins.
!
!  Example:
!
!    NBIN = 5, A = 5, B = 15
!
!    <-1-+-2-+-3-+-4-+-5->
!    5   7   9  11  13  15
!
!
!    C   BIN
!
!    1    1
!    3    1
!    4.9  1
!    5    1
!    6    1
!    7.1  2
!    8    2
!    9.5  3
!   12    4
!   14    5
!   15    5
!   15.1  5
!   99    5
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 February 2001
!
!  Author:
!
!    John Burkardt, Max Gunzburger, Janet Peterson
!
!  Parameters:
!
!    Input, integer NBIN, the number of bins.
!
!    Input, real ( kind = 8 ) A, B, the lower and upper limits of the bin
!    interval.  While A is expected to be less than B, the code should
!    return useful results if A is actually greater than B.
!
!    Input, real ( kind = 8 ) C, a value to be placed in a bin.
!
!    Output, integer BIN, the index of the bin to which C is assigned.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) a2
  real ( kind = 8 ) b
  real ( kind = 8 ) b2
  integer bin
  real ( kind = 8 ) c
  integer nbin
  logical switch
!
!  Take care of special cases.
!
  if ( nbin < 1 ) then
    bin = 0
    return
  end if

  if ( nbin == 1 ) then
    bin = 1
    return
  end if

  if ( b == a ) then
    bin = 1
    return
  end if
!
!  If the limits are descending, then we switch them now, and
!  unswitch the results at the end.
!
  if ( a < b ) then
    switch = .false.
    a2 = a
    b2 = b
  else
    switch = .true.
    a2 = b
    b2 = a
  end if
!
!  Compute the bin.
!
  if ( c <= a2 ) then
    bin = 1
  else if ( b2 <= c ) then
    bin = nbin
  else
    bin = 1 + int ( real ( nbin, kind = 8 ) * ( c - a2 ) / ( b2 - a2 ) )
    bin = max ( bin, 1 )
    bin = min ( bin, nbin )
  end if
!
!  Reverse the switching.
!
  if ( switch ) then
    bin = nbin + 1 - bin
  end if

  return
end
subroutine r82vec_bin_even3 ( n, a, nbin, bin_min, bin_max, bin_start, &
  bin_last, bin_next )

!*****************************************************************************80
!
!! R82VEC_BIN_EVEN3 bins an R82VEC into evenly spaced bins.
!
!  Discussion:
!
!    A different number of bins may be used in each dimension.
!
!    This is only a partial, indexed, sorting of the data.  To sort
!    the data, it is necessary to build a new array by extracting the
!    data for each bin, sorting that, and appending it to the array.
!
!    There are NBIN(1) 1D bins in the X direction, NBIN(2) for Y, making a
!    total of NBIN(1) * NBIN(2) 2D bins.  Each set of 1D bins begins and
!    ends at user specified mininum and maximum values.
!
!    The 2D bins are indexed by the X and Y bins that construct them,
!    and ordered lexicographically by these indices:
!
!      1,4 | 2,4 | 3,4 | 4,4 | 5,4
!      ----+-----+-----+-----+-----
!      1,3 | 2,3 | 3,3 | 4,3 | 5,3
!      ----+-----+-----+-----+-----
!      1,2 | 2,2 | 3,2 | 4,2 | 5,2
!      ----+-----+-----+-----+-----
!      1,1 | 2,1 | 3,1 | 4,1 | 5,1
!
!    Thus, the 2D bin sequence is (1,1), (1,2), (1,3), (1,4), (2,1),
!    ..., (5,4).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 March 2001
!
!  Author:
!
!    John Burkardt, Max Gunzburger, Janet Peterson
!
!  Reference:
!
!    Thomas Cormen, Charles Leiserson, Ronald Rivest,
!    Introduction to Algorithms,
!    MIT Press, page 180.
!
!  Parameters:
!
!    Input, integer N, the number of points in the data set.
!
!    Input, real ( kind = 8 ) A(2,N), the D2 data to be binned.
!
!    Input, integer NBIN(2), the number of bins in each dimension.
!    NBIN must be at least 1.
!
!    Input, real ( kind = 8 ) BIN_MIN(2), BIN_MAX(2), the bin limits.
!
!    Output, integer BIN_START(NBIN(1),NBIN(2)), BIN_LAST(NBIN(1),NBIN(2)),
!    the index of the first and last elements of A that went into each bin,
!    or -1 if there are no entries in the bin.
!
!    Output, integer BIN_NEXT(N), contains the index of the next element of A
!    that follows this element in the same bin.  A value of 0 means this
!    is the last entry in the particular bin.
!
  implicit none

  integer, parameter :: ndim = 2

  integer n
  integer nbin(ndim)

  real ( kind = 8 ) a(ndim,n)
  integer bin(ndim)
  integer bin_last(nbin(1),nbin(2))
  integer bin_next(n)
  integer bin_start(nbin(1),nbin(2))
  real ( kind = 8 ) bin_max(ndim)
  real ( kind = 8 ) bin_min(ndim)
  integer i1
  integer i2
  integer j
  integer k
!
!  Initialize the bin pointers to -1.
!
  bin_last(1:nbin(1),1:nbin(2)) = -1
  bin_start(1:nbin(1),1:nbin(2)) = -1
  bin_next(1:n) = -1
!
!  Build up linked lists of entries that go into each bin.
!
  do j = 1, n

    call r8vec_to_bin_even3 ( ndim, nbin, bin_min, bin_max, a(1:ndim,j), bin )

    i1 = bin(1)
    i2 = bin(2)

    if ( bin_start(i1,i2) == -1 ) then
      bin_start(i1,i2) = j
    else
      k = bin_last(i1,i2)
      bin_next(k) = j
    end if

    bin_next(j) = 0

    bin_last(i1,i2) = j

  end do

  return
end
subroutine r82vec_binned_reorder2 ( n, a, nbin, bin_start, bin_last, bin_next )

!*****************************************************************************80
!
!! R82VEC_BINNED_REORDER2 reorders a binned R82VEC.
!
!  Discussion:
!
!    This routine allows there to be a different number of bins in
!    each dimension.
!
!    The bin vectors define an implicit ordering of the data
!    vector.  This routine physically rearranges the data vector
!    so that items in the first bin come first, and so on.  The
!    data within a bin is not reordered.
!
!    On output, the BIN_START and BIN_NEXT arrays have also been updated
!    so that they still correspond to the (rearranged) vector A.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 March 2001
!
!  Author:
!
!    John Burkardt, Max Gunzburger, Janet Peterson
!
!  Parameters:
!
!    Input, integer N, the number of data points.
!
!    Input/output, real ( kind = 8 ) A(2,N), the D2 data to be sorted.
!
!    Input, integer NBIN(2), the number of bins in each direction.
!
!    Input/output, integer BIN_START(NBIN(1),NBIN(2)),
!    BIN_LAST(NBIN(1),NBIN(2)), contains the index of the first and last 
!    element of A that went into each bin, or -1 if there are no entries 
!    in the bin.
!
!    Input/output, integer BIN_NEXT(N), contains the index of the next
!    element of A that follows this element in the same bin.  A value of
!    0 means this is the last entry in the particular bin.
!
  implicit none

  integer, parameter :: ndim = 2

  integer n
  integer nbin(ndim)

  real ( kind = 8 ) a(ndim,n)
  real ( kind = 8 ) a2(ndim,n)
  integer bin_last(nbin(1),nbin(2))
  integer bin_next(n)
  integer bin_start(nbin(1),nbin(2))
  integer i
  integer i1
  integer i2
  integer j
  integer k
!
!  Bin by bin, copy the contents of A to A2.
!  The BIN_START array is also updated as we go.
!
  k = 0

  do i1 = 1, nbin(1)

    do i2 = 1, nbin(2)

      j = bin_start(i1,i2)

      if ( 0 < j ) then
        bin_start(i1,i2) = k + 1
      end if

      do while ( 0 < j )
        k = k + 1
        bin_last(i1,i2) = k
        a2(1:ndim,k) = a(1:ndim,j)
        j = bin_next(j)
      end do

    end do

  end do
!
!  Copy A2 back into A.
!
  a(1:ndim,1:n) = a2(1:ndim,1:n)
!
!  Now update the BIN_NEXT array.
!
  do i = 1, n
    bin_next(i) = i+1
  end do

  do i1 = 1, nbin(1)
    do i2 = 1, nbin(2)

      k = bin_last(i1,i2)

      if ( 0 < k ) then
        bin_next(k) = 0
      end if

    end do
  end do

  return
end
subroutine r82vec_binned_sort_a2 ( n, a, nbin, bin_start, bin_last )

!*****************************************************************************80
!
!! R82VEC_BINNED_SORT_A2 sorts each bin of a binned R82VEC.
!
!  Discussion:
!
!    This routine allows a different number of bins in each dimension.
!
!    Presumably, the data vector was first binned by R82VEC_BIN_EVEN3,
!    then reordered by R82VEC_BINNED_REORDER2.  Now, each of the
!    bins of data is sorted one at a time.
!
!    The result is NOT a lexicographically sorted D2 vector.
!
!    What is true is that if I < J, then either the I-th element of A occurs
!    in a lexicographically smaller bin than J, or they share a bin,
!    and the I-th element is lexicographically less than or equal to
!    the J-th element.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 March 2001
!
!  Author:
!
!    John Burkardt, Max Gunzburger, Janet Peterson
!
!  Parameters:
!
!    Input, integer N, the number of data points.
!
!    Input/output, real ( kind = 8 ) A(2,N), the R2 data to be sorted.
!
!    Input, integer NBIN(2), the number of bins in each dimension.
!
!    Input, integer BIN_START(NBIN(1),NBIN(2)), BIN_LAST(NBIN(1),NBIN(2)),
!    the index of the first and last element of A that went into each bin, or -1
!    if there are no entries in this bin.
!
  implicit none

  integer, parameter :: ndim = 2

  integer n
  integer nbin(ndim)

  real ( kind = 8 ) a(ndim,n)
  integer bin_last(nbin(1),nbin(2))
  integer bin_start(nbin(1),nbin(2))
  integer i1
  integer i2
  integer j1
  integer j2
  integer n1

  do i1 = 1, nbin(1)

    do i2 = 1, nbin(2)

      j1 = bin_start(i1,i2)

      if ( 0 < j1 ) then

        j2 = bin_last(i1,i2)

        n1 = j2 + 1 - j1

        if ( 1 < n1 ) then
          call r8col_sort_quick_a ( ndim, n1, a(1:ndim,j1:j2) )
        end if

      end if

    end do

  end do

  return
end
subroutine r83vec_bin_even3 ( n, a, nbin, bin_min, bin_max, bin_start, &
  bin_last, bin_next )

!*****************************************************************************80
!
!! R83VEC_BIN_EVEN3 bins an R83VEC into evenly spaced bins.
!
!  Discussion:
!
!    A different number of bins may be used in each dimension.
!
!    This is only a partial, indexed, sorting of the data.  To sort
!    the data, it is necessary to build a new array by extracting the
!    data for each bin, sorting that, and appending it to the array.
!
!    There are NBIN(1) 1D bins in the X direction, NBIN(2) for Y, and
!    NBIN(3) for Z, making a total of NBIN(1) * NBIN(2) * NBIN(3) 3D bins.
!    Each set of 1D bins begins and ends at user specified mininum and
!    maximum values.
!
!    The 3D bins are indexed by the X, Y and Z bins that construct them,
!    and ordered lexicographically by these indices.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 April 2001
!
!  Author:
!
!    John Burkardt, Max Gunzburger, Janet Peterson
!
!  Reference:
!
!    Thomas Cormen, Charles Leiserson, Ronald Rivest,
!    Introduction to Algorithms,
!    MIT Press, page 180.
!
!  Parameters:
!
!    Input, integer N, the number of points in the data set.
!
!    Input, real ( kind = 8 ) A(3,N), the R3 data to be binned.
!
!    Input, integer NBIN(3), the number of bins in each dimension.
!    NBIN must be at least 1.
!
!    Input, real ( kind = 8 ) BIN_MIN(3), BIN_MAX(3), the bin limits.
!
!    Output, integer BIN_START(NBIN(1),NBIN(2),NBIN(3)),
!    BIN_LAST(NBIN(1),NBIN(2),NBIN(3)), contains the
!    index of the first and last elements of A that went into each bin,
!    or -1 if there are no entries in the bin.
!
!    Output, integer BIN_NEXT(N), contains the index of the next element of A
!    that follows this element in the same bin.  A value of 0 means this
!    is the last entry in the particular bin.
!
  implicit none

  integer, parameter :: ndim = 3

  integer n
  integer nbin(ndim)

  real ( kind = 8 ) a(ndim,n)
  integer bin(ndim)
  integer bin_last(nbin(1),nbin(2),nbin(3))
  integer bin_next(n)
  integer bin_start(nbin(1),nbin(2),nbin(3))
  real ( kind = 8 ) bin_max(ndim)
  real ( kind = 8 ) bin_min(ndim)
  integer i1
  integer i2
  integer i3
  integer j
  integer k
!
!  Initialize the bin pointers to -1.
!
  bin_last(1:nbin(1),1:nbin(2),1:nbin(3)) = -1
  bin_start(1:nbin(1),1:nbin(2),1:nbin(3)) = -1
  bin_next(1:n) = -1
!
!  Build up linked lists of entries that go into each bin.
!
  do j = 1, n

    call r8vec_to_bin_even3 ( ndim, nbin, bin_min, bin_max, a(1:ndim,j), bin )

    i1 = bin(1)
    i2 = bin(2)
    i3 = bin(3)

    if ( bin_start(i1,i2,i3) == -1 ) then
      bin_start(i1,i2,i3) = j
    else
      k = bin_last(i1,i2,i3)
      bin_next(k) = j
    end if

    bin_next(j) = 0

    bin_last(i1,i2,i3) = j

  end do

  return
end
subroutine r83vec_binned_reorder2 ( n, a, nbin, bin_start, bin_last, bin_next )

!*****************************************************************************80
!
!! R83VEC_BINNED_REORDER2 reorders a binned R83VEC.
!
!  Discussion:
!
!    This routine allows there to be a different number of bins in
!    each dimension.
!
!    The bin vectors define an implicit ordering of the data
!    vector.  This routine physically rearranges the data vector
!    so that items in the first bin come first, and so on.  The
!    data within a bin is not reordered.
!
!    On output, the BIN_START and BIN_NEXT arrays have also been updated
!    so that they still correspond to the (rearranged) vector A.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 April 2001
!
!  Author:
!
!    John Burkardt, Max Gunzburger, Janet Peterson
!
!  Parameters:
!
!    Input, integer N, the number of data points.
!
!    Input/output, real ( kind = 8 ) A(3,N), the data to be sorted.
!
!    Input, integer NBIN(3), the number of bins in each dimension.
!
!    Input/output, integer BIN_START(NBIN(1),NBIN(2),NBIN(3)),
!    BIN_LAST(NBIN(1),NBIN(2),NBIN(3)), contains the
!    index of the first and last element of A that went into each bin,
!    or -1 if there are no entries in the bin.
!
!    Input/output, integer BIN_NEXT(N), contains the index of the next
!    element of A that follows this element in the same bin.  A value of 0
!    means this is the last entry in the particular bin.
!
  implicit none

  integer, parameter :: ndim = 3

  integer n
  integer nbin(ndim)

  real ( kind = 8 ) a(ndim,n)
  real ( kind = 8 ) a2(ndim,n)
  integer bin_last(nbin(1),nbin(2),nbin(3))
  integer bin_next(n)
  integer bin_start(nbin(1),nbin(2),nbin(3))
  integer i
  integer i1
  integer i2
  integer i3
  integer j
  integer k
!
!  Bin by bin, copy the contents of A to A2.
!  The BIN_START array is also updated as we go.
!
  k = 0

  do i1 = 1, nbin(1)

    do i2 = 1, nbin(2)

      do i3 = 1, nbin(3)

        j = bin_start(i1,i2,i3)

        if ( 0 < j ) then
          bin_start(i1,i2,i3) = k + 1
        end if

        do while ( 0 < j )
          k = k + 1
          bin_last(i1,i2,i3) = k
          a2(1:ndim,k) = a(1:ndim,j)
          j = bin_next(j)
        end do

      end do

    end do

  end do
!
!  Copy A2 back into A.
!
  a(1:ndim,1:n) = a2(1:ndim,1:n)
!
!  Now update the BIN_NEXT array.
!
  do i = 1, n
    bin_next(i) = i + 1
  end do

  do i1 = 1, nbin(1)
    do i2 = 1, nbin(2)
      do i3 = 1, nbin(3)

        k = bin_last(i1,i2,i3)

        if ( 0 < k ) then
          bin_next(k) = 0
        end if

      end do
    end do
  end do

  return
end
subroutine r83vec_binned_sort_a2 ( n, a, nbin, bin_start, bin_last )

!*****************************************************************************80
!
!! R83VEC_BINNED_SORT_A2 sorts each bin of a binned R83VEC.
!
!  Discussion:
!
!    This routine allows a different number of bins in each dimension.
!
!    Presumably, the data vector was first binned by R83VEC_BIN_EVEN3,
!    then reordered by R83VEC_BINNED_REORDER2.  Now, each of the
!    bins of data is sorted one at a time.
!
!    The result is NOT a lexicographically sorted D3 vector.
!
!    What is true is that if I < J, then either the I-th element of A occurs
!    in a lexicographically smaller bin than J, or they share a bin,
!    and the I-th element is lexicographically less than or equal to
!    the J-th element.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 April 2001
!
!  Author:
!
!    John Burkardt, Max Gunzburger, Janet Peterson
!
!  Parameters:
!
!    Input, integer N, the number of data points.
!
!    Input/output, real ( kind = 8 ) A(3,N), the data to be sorted.
!
!    Input, integer NBIN(3), the number of bins in each dimension.
!
!    Input, integer BIN_START(NBIN(1),NBIN(2),NBIN(3)),
!    BIN_LAST(NBIN(1),NBIN(2),NBIN(3)), contains
!    the index of the first and last element of A that went into each bin,
!    or -1 if there are no entries in this bin.
!
  implicit none

  integer, parameter :: ndim = 3

  integer n
  integer nbin(ndim)

  real ( kind = 8 ) a(ndim,n)
  integer bin_last(nbin(1),nbin(2),nbin(3))
  integer bin_start(nbin(1),nbin(2),nbin(3))
  integer i1
  integer i2
  integer i3
  integer j1
  integer j2
  integer n1

  do i1 = 1, nbin(1)

    do i2 = 1, nbin(2)

      do i3 = 1, nbin(3)

        j1 = bin_start(i1,i2,i3)

        if ( 0 < j1 ) then

          j2 = bin_last(i1,i2,i3)

          n1 = j2 + 1 - j1

          if ( 1 < n1 ) then
            call r8col_sort_quick_a ( ndim, n1, a(1:ndim,j1:j2) )
          end if

        end if

      end do

    end do

  end do

  return
end
subroutine r8col_part_quick_a ( m, n, a, l, r )

!*****************************************************************************80
!
!! R8COL_PART_QUICK_A reorders the columns of an array as part of a quick sort.
!
!  Discussion:
!
!    The routine reorders the columns of A.  Using A(1:M,1) as a
!    key, all entries of A that are less than or equal to the key will
!    precede the key, which precedes all entries that are greater than the key.
!
!  Example:
!
!    Input:
!
!      M = 2, N = 8
!      A = ( 2  8  6  0 10 10  0  5
!            4  8  2  2  6  0  6  8 )
!
!    Output:
!
!      L = 2, R = 4
!
!      A = (  0  0  2  8  6 10 10  4
!             2  6  4  8  2  6  0  8 )
!             ----     -------------
!             LEFT KEY     RIGHT
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 April 2001
!
!  Author:
!
!    John Burkardt, Max Gunzburger, Janet Peterson
!
!  Parameters:
!
!    Input, integer M, the row dimension of A, and the length of a column.
!
!    Input, integer N, the column dimension of A.
!
!    Input/output, real ( kind = 8 ) A(M,N).  On input, the array to be checked.
!    On output, A has been reordered as described above.
!
!    Output, integer L, R, the indices of A that define the three segments.
!    Let KEY = the input value of A(1:M,1).  Then
!    I <= L                 A(1:M,I) < KEY;
!         L < I < R         A(1:M,I) = KEY;
!                 R <= I    A(1:M,I) > KEY.
!
  implicit none

  integer m
  integer n

  real ( kind = 8 ) a(m,n)
  integer i
  integer j
  integer k
  real ( kind = 8 ) key(m)
  integer l
  integer r
  logical r8vec_eq
  logical r8vec_gt
  logical r8vec_lt

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8COL_PART_QUICK_A - Fatal error!'
    write ( *, '(a)' ) '  N < 1.'
    return
  end if

  if ( n == 1 ) then
    l = 0
    r = 2
    return
  end if

  key(1:m) = a(1:m,1)
  k = 1
!
!  The elements of unknown size have indices between L+1 and R-1.
!
  l = 1
  r = n + 1

  do i = 2, n

    if ( r8vec_gt ( m, a(1:m,l+1), key(1:m) ) ) then
      r = r - 1
      call r8vec_swap ( m, a(1:m,r), a(1:m,l+1) )
    else if ( r8vec_eq ( m, a(1:m,l+1), key(1:m) ) ) then
      k = k + 1
      call r8vec_swap ( m, a(1:m,k), a(1:m,l+1) )
      l = l + 1
    else if ( r8vec_lt ( m, a(1:m,l+1), key(1:m) ) ) then
      l = l + 1
    end if

  end do
!
!  Shift small elements to the left.
!
  do j = 1, l - k
    a(1:m,j) = a(1:m,j+k)
  end do
!
!  Shift KEY elements to center.
!
  do j = l-k+1, l
    a(1:m,j) = key(1:m)
  end do
!
!  Update L.
!
  l = l - k

  return
end
subroutine r8col_sort_quick_a ( m, n, a )

!*****************************************************************************80
!
!! R8COL_SORT_QUICK_A ascending sorts the columns of a table using quick sort.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 April 2001
!
!  Author:
!
!    John Burkardt, Max Gunzburger, Janet Peterson
!
!  Parameters:
!
!    Input, integer M, the row order of A, and the length of a column.
!
!    Input, integer N, the number of columns of A.
!
!    Input/output, real ( kind = 8 ) A(M,N).
!    On input, the array to be sorted.
!    On output, the array has been sorted.
!
  implicit none

  integer, parameter :: MAXLEVEL = 25

  integer m
  integer n

  real ( kind = 8 ) a(m,n)
  integer base
  integer l_segment
  integer level
  integer n_segment
  integer rsave(MAXLEVEL)
  integer r_segment

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8COL_SORT_QUICK_A - Fatal error!'
    write ( *, '(a)' ) '  N < 1.'
    stop
  end if

  if ( n == 1 ) then
    return
  end if

  level = 1
  rsave(level) = n + 1
  base = 1
  n_segment = n

  do
!
!  Partition the segment.
!
    call r8col_part_quick_a ( m, n_segment, a(1,base), l_segment, r_segment )
!
!  If the left segment has more than one element, we need to partition it.
!
    if ( 1 < l_segment ) then

      if ( MAXLEVEL < level ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8COL_SORT_QUICK_A - Fatal error!'
        write ( *, '(a,i8)' ) '  Exceeding recursion maximum of ', MAXLEVEL
        stop
      end if

      level = level + 1
      n_segment = l_segment
      rsave(level) = r_segment + base - 1
!
!  The left segment and the middle segment are sorted.
!  Must the right segment be partitioned?
!
    else if ( r_segment < n_segment ) then

      n_segment = n_segment + 1 - r_segment
      base = base + r_segment - 1
!
!  Otherwise, we back up a level if there is an earlier one.
!
    else

      do

        if ( level <= 1 ) then
          return
        end if

        base = rsave(level)
        n_segment = rsave(level-1) - rsave(level)
        level = level - 1

        if ( 0 < n_segment ) then
          exit
        end if

      end do

    end if

  end do

  return
end
function r8mat_det_2d ( a )

!*****************************************************************************80
!
!! R8MAT_DET_2D computes the determinant of a 2 by 2 matrix.
!
!  Discussion:
!
!    The determinant of a 2 by 2 matrix is
!
!      a11 * a22 - a12 * a21.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 March 1999
!
!  Author:
!
!    John Burkardt, Max Gunzburger, Janet Peterson
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A(2,2), the matrix whose determinant is desired.
!
!    Output, real ( kind = 8 ) RMAT_DET_2D, the determinant of the matrix.
!
  implicit none

  real ( kind = 8 ) a(2,2)
  real ( kind = 8 ) r8mat_det_2d

  r8mat_det_2d = a(1,1) * a(2,2) - a(1,2) * a(2,1)

  return
end
function r8mat_det_3d ( a )

!*****************************************************************************80
!
!! R8MAT_DET_3D computes the determinant of a 3 by 3 matrix.
!
!  Discussion:
!
!    The determinant of a 3 by 3 matrix is
!
!        a11 * a22 * a33 - a11 * a23 * a32
!      + a12 * a23 * a31 - a12 * a21 * a33
!      + a13 * a21 * a32 - a13 * a22 * a31
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 March 1999
!
!  Author:
!
!    John Burkardt, Max Gunzburger, Janet Peterson
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A(3,3), the matrix whose determinant is desired.
!
!    Output, real ( kind = 8 ) RMAT_DET_3D, the determinant of the matrix.
!
  implicit none

  real ( kind = 8 ) a(3,3)
  real ( kind = 8 ) r8mat_det_3d

  r8mat_det_3d = &
         a(1,1) * ( a(2,2) * a(3,3) - a(2,3) * a(3,2) ) &
       + a(1,2) * ( a(2,3) * a(3,1) - a(2,1) * a(3,3) ) &
       + a(1,3) * ( a(2,1) * a(3,2) - a(2,2) * a(3,1) )

  return
end
function r8vec_eq ( n, a1, a2 )

!*****************************************************************************80
!
!! R8VEC_EQ is true if every pair of entries in two vectors is equal.
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
!    John Burkardt, Max Gunzburger, Janet Peterson
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vectors.
!
!    Input, real ( kind = 8 ) A1(N), A2(N), two vectors to compare.
!
!    Output, logical R8VEC_EQ.
!    R8VEC_EQ is .TRUE. if every pair of elements A1(I) and A2(I) are equal,
!    and .FALSE. otherwise.
!
  implicit none

  integer n

  real ( kind = 8 ) a1(n)
  real ( kind = 8 ) a2(n)
  logical r8vec_eq

  r8vec_eq = ( all ( a1(1:n) == a2(1:n) ) )

  return
end
function r8vec_gt ( n, a1, a2 )

!*****************************************************************************80
!
!! R8VEC_GT == ( A1 > A2 ) for two R8VEC's.
!
!  Discussion:
!
!    The comparison is lexicographic.
!
!    A1 > A2  <=>                              A1(1) > A2(1) or
!                 ( A1(1)     == A2(1)     and A1(2) > A2(2) ) or
!                 ...
!                 ( A1(1:N-1) == A2(1:N-1) and A1(N) > A2(N)
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
!    John Burkardt, Max Gunzburger, Janet Peterson
!
!  Parameters:
!
!    Input, integer N, the dimension of the vectors.
!
!    Input, real ( kind = 8 ) A1(N), A2(N), the vectors to be compared.
!
!    Output, logical R8VEC_GT, is TRUE if and only if A1 > A2.
!
  implicit none

  integer n

  real ( kind = 8 ) a1(n)
  real ( kind = 8 ) a2(n)
  integer i
  logical r8vec_gt

  r8vec_gt = .false.

  do i = 1, n

    if ( a2(i) < a1(i) ) then
      r8vec_gt = .true.
      exit
    else if ( a1(i) < a2(i) ) then
      r8vec_gt = .false.
      exit
    end if

  end do

  return
end
function r8vec_lt ( n, a1, a2 )

!*****************************************************************************80
!
!! R8VEC_LT == ( A1 < A2 ) for double precision vectors.
!
!  Discussion:
!
!    The comparison is lexicographic.
!
!    A1 < A2  <=>                              A1(1) < A2(1) or
!                 ( A1(1)     == A2(1)     and A1(2) < A2(2) ) or
!                 ...
!                 ( A1(1:N-1) == A2(1:N-1) and A1(N) < A2(N)
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
!    John Burkardt, Max Gunzburger, Janet Peterson
!
!  Parameters:
!
!    Input, integer N, the dimension of the vectors.
!
!    Input, real ( kind = 8 ) A1(N), A2(N), the vectors to be compared.
!
!    Output, logical R8VEC_LT, is TRUE if and only if A1 < A2.
!
  implicit none

  integer n

  real ( kind = 8 ) a1(n)
  real ( kind = 8 ) a2(n)
  integer i
  logical r8vec_lt

  r8vec_lt = .false.

  do i = 1, n

    if ( a1(i) < a2(i) ) then
      r8vec_lt = .true.
      exit
    else if ( a2(i) < a1(i) ) then
      r8vec_lt = .false.
      exit
    end if

  end do

  return
end
subroutine r8vec_swap ( n, a1, a2 )

!*****************************************************************************80
!
!! R8VEC_SWAP swaps the entries of two double precision vectors.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 April 2001
!
!  Author:
!
!    John Burkardt, Max Gunzburger, Janet Peterson
!
!  Parameters:
!
!    Input, integer N, the number of entries in the arrays.
!
!    Input/output, real ( kind = 8 ) A1(N), A2(N), the vectors to swap.
!
  implicit none

  integer n

  real ( kind = 8 ) a1(n)
  real ( kind = 8 ) a2(n)
  real ( kind = 8 ) a3(n)

  a3(1:n) = a1(1:n)
  a1(1:n) = a2(1:n)
  a2(1:n) = a3(1:n)

  return
end
subroutine r8vec_to_bin_even3 ( ndim, nbin, a, b, c, bin )

!*****************************************************************************80
!
!! R8VEC_TO_BIN_EVEN3 determines the appropriate "bin" for a R8VEC value.
!
!  Discussion:
!
!    The interval [A(I),B(I)] is divided into NBIN(I) equal subintervals
!    or bins.
!
!  Example:
!
!    NDIM = 3
!    NBIN = (/ 4, 5, 2 /),
!
!      A(1) = 1,  A(2) = 0,  A(3) = 8
!      B(1) = 17, B(2) = 20, B(3) = 10
!
!
!            8 < Z < 9                    9 < Z < 10
!
!   20 +     +     +     +     +     20 +     +     +     +     +
!        151 | 251 | 351 | 451            152 | 252 | 352 | 452
!   16 +-----+-----+-----+-----+     16 +-----+-----+-----+-----+
!        141 | 241 | 341 | 441            142 | 242 | 342 | 442
!   12 +-----+-----+-----+-----+     12 +-----+-----+-----+-----+
!        131 | 231 | 331 | 431            132 | 232 | 332 | 432
!    8 +-----+-----+-----+-----+      8 +-----+-----+-----+-----+
!        121 | 221 | 321 | 421            122 | 222 | 322 | 422
!    4 +-----+-----+-----+-----+      4 +-----+-----+-----+-----+
!        111 | 211 | 311 | 411            112 | 212 | 312 | 412
!    0 +     +     +     +     +      0 +     +     +     +     +
!      1     5     9    13    17        1     5     9    13    17
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 April 2001
!
!  Author:
!
!    John Burkardt, Max Gunzburger, Janet Peterson
!
!  Parameters:
!
!    Input, integer NDIM, the dimension of the space.
!
!    Input, integer NBIN(NDIM), the number of bins in each dimension.
!
!    Input, real ( kind = 8 ) A(NDIM), B(NDIM), the lower and upper limits of
!    the bin interval.  While A(I) is expected to be less than B(I), the code
!    should return useful results if A(I) is actually greater than B(I).
!
!    Input, real ( kind = 8 ) C(NDIM), a value to be placed in a bin.
!
!    Output, integer BIN(NDIM), the index of the bin to which C is assigned.
!
  implicit none

  integer ndim

  real ( kind = 8 ) a(ndim)
  real ( kind = 8 ) b(ndim)
  integer bin(ndim)
  real ( kind = 8 ) c(ndim)
  integer i
  integer nbin(ndim)

  do i = 1, ndim
    call r8_to_bin_even2 ( nbin(i), a(i), b(i), c(i), bin(i) )
  end do

  return
end
subroutine random_initialize ( seed )

!*****************************************************************************80
!
!! RANDOM_INITIALIZE initializes the FORTRAN90 random number seed.
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
!    John Burkardt, Max Gunzburger, Janet Peterson
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
  real t
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
subroutine region_sampler ( ndim, box_min, box_max, dr, x, &
  random_function, reset, use_diatom, ngen )

!*****************************************************************************80
!
!! REGION_SAMPLER returns a sample point in the physical region.
!
!  Discussion:
!
!    The calculations are done in NDIM dimensional space.
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
!    Right now, we're assuming NDIM is no greater than 3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 April 2001
!
!  Author:
!
!    John Burkardt, Max Gunzburger, Janet Peterson
!
!  Reference:
!
!    John Burkardt, Max Gunzburger, Janet Peterson and Rebecca Brannon,
!    User Manual and Supporting Information for Library of Codes
!    for Centroidal Voronoi Placement and Associated Zeroth,
!    First, and Second Moment Determination,
!    Sandia National Laboratories Technical Report SAND2002-0099,
!    February 2002.
!
!  Parameters:
!
!    Input, integer NDIM, the spatial dimension.
!
!    Input, real ( kind = 8 ) BOX_MIN(NDIM), BOX_MAX(NDIM), the coordinates
!    of the two extreme corners of the bounding box.
!
!    Input, real ( kind = 8 ) DR, a tolerance used by DIATOM when testing
!    whether a point is within, outside of, or on the boundary of the
!    physical region.
!
!    Output, real ( kind = 8 ) X(NDIM), the random point.
!
!    Input, integer RANDOM_FUNCTION, specifies the random function.
!    0, uniform random numbers from F90 RANDOM_NUMBER.
!    1, Halton sequence.
!
!    Input/output, logical RESET.
!    If TRUE, the Halton sequence should be reset.
!
!    Input, logical USE_DIATOM, is TRUE if DIATOM is to be called 
!    to determine whether a point lies in the physical region; if it is
!    FALSE than a much simplified routine is used.
!
!    Output, integer NGEN, the number of points that were generated.
!    This is at least 1, but may be larger if some points were rejected.
!
  implicit none

  integer ndim

  real ( kind = 8 ) box_max(ndim)
  real ( kind = 8 ) box_min(ndim)
  real ( kind = 8 ) dr
  integer, save, dimension ( 3 ) :: halton_base = (/ 2, 3, 5 /)
  integer, save :: halton_seed = 1
  integer ival
  real ( kind = 8 ) mdens
  integer ngen
  real ( kind = 8 ) r(ndim)
  integer random_function
  logical reset
  logical use_diatom
  real ( kind = 8 ) x(ndim)
  real ( kind = 8 ), parameter :: zero = 0.0D+00

  ngen = 0

  if ( reset ) then
    halton_seed = 1
    reset = .false.
  end if

  do

    ngen = ngen + 1

    if ( 10000 < ngen ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'REGION_SAMPLER - Fatal error!'
      write ( *, '(a,i12,a)') '  Generated ', ngen, ' rejected points in a row.'
      write ( *, '(a)' ) &
        '  There may be a problem with the geometry definition.'
      write ( *, '(a)' ) ' '
      if ( random_function == 0 ) then
        write ( *, '(a)' ) '  Using F90 RANDOM_NUMBER.'
      else if ( random_function == 1 ) then
        write ( *, '(a)' ) '  Using Halton sequence.'
        write ( *, '(a,i12)' ) '  Current seed is ', halton_seed
      end if
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Current random value is:'
      write ( *, '(3g14.6)' ) r(1:ndim)
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Current sample point is:'
      write ( *, '(3g14.6)' ) x(1:ndim)
      stop
    end if
!
!  Generate a point at random using:
!  0: a uniformly distributed random value;
!  1: a Halton random value.
!
    if ( random_function == 0 ) then

      call random_number ( r(1:ndim) )

    else if ( random_function == 1 ) then

      call i4_to_halton_vector ( halton_seed, halton_base, ndim, r )
      halton_seed = halton_seed + 1

    else

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'REGION_SAMPLER - Fatal error!'
      write ( *, '(a,i8)' ) '  Illegal value of RANDOM_FUNCTION = ', &
        random_function
      stop

    end if
!
!  Determine a point in the bounding box.
!
    x(1:ndim) = ( ( 1.0D+00 - r(1:ndim) ) * box_min(1:ndim) &
                            + r(1:ndim)   * box_max(1:ndim) )
!
!  Now determine if the point is in the region.
!
    if ( use_diatom ) then

      if ( ndim == 2 ) then
        call diatom_point_test2 ( x(1), x(2), zero, dr, mdens, ival )
      else if ( ndim == 2 ) then
        call diatom_point_test2 ( x(1), x(2), x(3), dr, mdens, ival )
      end if
!
!  Call the routine that has an analytic definition of the region:
!
    else

      call test_region ( x, ndim, ival )

    end if

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
!    John Burkardt, Max Gunzburger, Janet Peterson
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
subroutine vcm ( ndim, box_min, box_max, n, cell_generator, ns_mom, &
  use_diatom, use_bins, region_volume_given, region_volume, dr, nbin, &
  bin_start, bin_last, bin_next, cell_volume, cell_centroid, cell_moment )

!*****************************************************************************80
!
!! VCM calculates Voronoi cell volumes, centroids and second moments.
!
!  Discussion:
!
!    A Monte Carlo sampling is used to estimate the quantities.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 April 2001
!
!  Author:
!
!    John Burkardt, Max Gunzburger, Janet Peterson
!
!  Reference:
!
!    John Burkardt, Max Gunzburger, Janet Peterson and Rebecca Brannon,
!    User Manual and Supporting Information for Library of Codes
!    for Centroidal Voronoi Placement and Associated Zeroth,
!    First, and Second Moment Determination,
!    Sandia National Laboratories Technical Report SAND2002-0099,
!    February 2002.
!
!  Parameters:
!
!    Input, integer NDIM, the spatial dimension.
!
!    Input, real ( kind = 8 ) BOX_MIN(NDIM), BOX_MAX(NDIM), the coordinates
!    of the two extreme corners of the bounding box.
!
!    Input, integer N, the number of cell generators.
!
!    Input, real ( kind = 8 ) CELL_GENERATOR(NDIM,N), the cell generators.
!
!    Input, integer NS_MOM, the number of sample points per cell generator.
!
!    Input, logical USE_DIATOM, is TRUE if DIATOM is to be called 
!    to determine whether a point lies in the physical region; if it is
!    FALSE than a much simplified routine is used.
!
!    Input, logical USE_BINS, is TRUE if the bounding box is to be divided
!    up into bins to speed up the nearest neighbor search;
!    FALSE if the nearest neighbor seach is to be done naively.
!
!    Input, logical REGION_VOLUME_GIVEN,
!    TRUE: the region volume is input in REGION_VOLUME.
!    FALSE: the region volume must be estimated by this routine.
!
!    Input/output, real ( kind = 8 ) REGION_VOLUME, the volume of the region.
!    If REGION_VOLUME_GIVEN is TRUE, then REGION_VOLUME is input by the user.
!    Otherwise, the volume is estimated and output by this routine.
!
!    Input, real ( kind = 8 ) DR, a tolerance used by DIATOM when testing
!    whether a point is within, outside of, or on the boundary of the
!    physical region.
!
!    Input, integer NBIN(3) is the number of bins to use in each direction.
!    For 3D problems, set NBIN(3) = 1.
!    For efficiency, these values should be set in such a way that the bins
!    are nearly square or cubical.
!
!    Input, integer BIN_START(NBIN(1),NBIN(2),NBIN(3)), the index of the first
!    cell generator in the bin, or -1 if none.
!
!    Input, integer BIN_LAST(NBIN(1),NBIN(2),NBIN(3)), the index of the last
!    cell generator in the bin, or -1 if none.
!
!    Input, integer BIN_NEXT(N), the index of the next cell generator in
!    the bin containing this cell generator.
!
!    Output, real ( kind = 8 ) CELL_VOLUME(N), the Voronoi cell volumes.
!
!    Output, real ( kind = 8 ) CELL_CENTROID(NDIM,N), the Voronoi cell
!    centroids.
!
!    Output, real ( kind = 8 ) CELL_MOMENT(NDIM,NDIM,N), the second moment
!    matrix for each Voronoi cell.
!
  implicit none

  integer nbin(3)
  integer ndim
  integer n

  integer bin_last(nbin(1),nbin(2),nbin(3))
  integer bin_next(n)
  integer bin_start(nbin(1),nbin(2),nbin(3))
  real ( kind = 8 ) box_max(ndim)
  real ( kind = 8 ) box_min(ndim)
  real ( kind = 8 ) box_volume
  real ( kind = 8 ) cell_centroid(ndim,n)
  real ( kind = 8 ) cell_generator(ndim,n)
  integer cell_hit(n)
  real ( kind = 8 ) cell_moment(ndim,ndim,n)
  real ( kind = 8 ) cell_volume(n)
  real ( kind = 8 ) dr
  integer i
  integer j
  integer k
  integer nearest
  integer ngen
  integer ns_mom
  integer ntries
  integer, parameter :: random_mom = 0
  real ( kind = 8 ) region_volume
  real ( kind = 8 ) region_volume_estimate
  logical region_volume_given
  logical reset
  logical use_bins
  logical use_diatom
  real ( kind = 8 ) x(ndim)
!
!  Zero out the arrays.
!
  cell_centroid(1:ndim,1:n) = 0.0D+00
  cell_hit(1:n) = 0
  cell_moment(1:ndim,1:ndim,1:n) = 0.0D+00
!
!  Sample the region N * NS_MOM times, and keep track of which cell generator
!  is closest to each sampling point.
!
  ntries = 0
  reset = .true.

  do k = 1, n * ns_mom

    call region_sampler ( ndim, box_min, box_max, dr, x, &
      random_mom, reset, use_diatom, ngen )
    ntries = ntries + ngen

    if ( use_bins ) then

      if ( ndim == 2 ) then

        call points_nearest_point_bins3_2d ( n, cell_generator, nbin, &
          box_min, box_max, bin_start, bin_last, bin_next, x, nearest )

      else if ( ndim == 3 ) then

        call points_nearest_point_bins3_3d ( n, cell_generator, nbin, &
          box_min, box_max, bin_start, bin_last, bin_next, x, nearest )

      end if

    else

      call find_closest ( ndim, x, n, cell_generator, nearest )

    end if

    cell_hit(nearest) = cell_hit(nearest) + 1
    cell_centroid(1:ndim,nearest) = cell_centroid(1:ndim,nearest) + x(1:ndim)

    do i = 1, ndim
      do j = 1, ndim
        cell_moment(i,j,nearest) = cell_moment(i,j,nearest) + x(i) * x(j)
      end do
    end do

  end do
!
!  Estimate the area or volume if it was not given.
!
  box_volume = product ( ( box_max(1:ndim)-box_min(1:ndim) ) )
  region_volume_estimate = real ( n * ns_mom, kind = 8 ) * box_volume &
    / real ( ntries, kind = 8 )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Volume of bounding box is     ', box_volume
  if ( region_volume_given ) then
    write ( *, '(a,g14.6)' ) '  Given volume of region is     ', region_volume
  else
    region_volume = region_volume_estimate
  end if
  write ( *, '(a,g14.6)' ) '  Estimated volume of region is ', &
    region_volume_estimate
!
!  Estimate the geometric integrals for each Voronoi cell.
!
  do k = 1, n

    if ( 0 < cell_hit(k) ) then

      cell_volume(k) = real ( cell_hit(k), kind = 8 ) * region_volume &
        / real ( n * ns_mom, kind = 8 )
      cell_centroid(1:ndim,k) = cell_centroid(1:ndim,k) &
        / real ( cell_hit(k), kind = 8 )
      cell_moment(1:ndim,1:ndim,k) = cell_moment(1:ndim,1:ndim,k) &
        / real ( cell_hit(k), kind = 8 )

      do i = 1, ndim
        do j = 1, ndim
          cell_moment(i,j,k) = cell_moment(i,j,k) &
            - cell_centroid(i,k) * cell_centroid(j,k)
        end do
      end do

    end if

  end do

  return
end
