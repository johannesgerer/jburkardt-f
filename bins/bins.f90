subroutine bin_search_one_2d ( bin, nset, pset, nbin, bin_start, bin_next, &
  ptest, found_a_neighbor, i_min, dist_min_sq, compares )

!*****************************************************************************80
!
!! BIN_SEARCH_ONE_2D searches one cell in a 2D array of bins.
!
!  Discussion:
!
!    The bins are presumed to have been set up by successive calls to:
!
!      R82VEC_BIN_EVEN2,
!      R82VEC_BINNED_REORDER, and
!      R82VEC_BINNED_SORT_A.
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
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) BIN(2), the indices of the cell to be examined.
!
!    Input, integer ( kind = 4 ) NSET, the number of points in the set.
!
!    Input, real ( kind = 8 ) PSET(2,NSET), the coordinates of the points
!    in the set.
!
!    Input, integer ( kind = 4 ) NBIN(2), the number of cells in the horizontal
!    and vertical directions.
!
!    Input, integer ( kind = 4 ) BIN_START(NBIN(1),NBIN(2)), 
!    BIN_LAST(NBIN(1),NBIN(2)), indicates the index of the first and last 
!    element in the bin, or -1 if none.
!
!    Input, integer ( kind = 4 ) BIN_NEXT(NSET), the index of the next element 
!    of the bin containing this element.
!
!    Input, real ( kind = 8 ) PTEST(2), the coordinates of the test point.
!
!    Input/output, logical FOUND_A_NEIGHBOR, is set to TRUE if at least
!    one point of PSET is found in the current bin.  Otherwise, it retains its
!    input value.
!
!    Input/output, integer ( kind = 4 ) I_MIN, the index of the nearest neighbor
!    in PSET to PTEST, if at least one neighbor has been found.
!
!    Input/output, real ( kind = 8 ) DIST_MIN_SQ, the square of the distance
!    from the nearest neighbor in PSET to PTEST, if at least one neighbor
!    has been found.
!
!    Input/output, integer ( kind = 4 ) COMPARES, the number of elements of 
!    PSET whose distance to PTEST has been computed.
!
  implicit none

  integer ( kind = 4 ), parameter :: ndim = 2

  integer ( kind = 4 ) nbin(ndim)
  integer ( kind = 4 ) nset

  integer ( kind = 4 ) bin(ndim)
  integer ( kind = 4 ) bin_next(nset)
  integer ( kind = 4 ) bin_start(nbin(1),nbin(2))
  integer ( kind = 4 ) compares
  real ( kind = 8 ) dist_min_sq
  real ( kind = 8 ) dist_sq
  logical found_a_neighbor
  integer ( kind = 4 ) i_min
  integer ( kind = 4 ) node
  real ( kind = 8 ) pset(ndim,nset)
  real ( kind = 8 ) ptest(ndim)

  node = bin_start(bin(1),bin(2))

  do while ( 0 < node )

    found_a_neighbor = .true.

    dist_sq = sum ( ( ptest(1:ndim) - pset(1:ndim,node) )**2 )
    compares = compares + 1

    if ( dist_sq < dist_min_sq ) then
      dist_min_sq = dist_sq
      i_min = node
    end if

    node = bin_next(node)

  end do

  return
end
subroutine bin_to_r8_even ( nbin, bin, a, b, cmin, cmax )

!*****************************************************************************80
!
!! BIN_TO_R8_EVEN returns the limits for a given "bin" in [A,B].
!
!  Discussion:
!
!    The interval from A to B is divided into NBIN-2 equal subintervals or bins.
!    An initial bin takes everything less than A, and a final bin takes
!    everything greater than B.
!
!  Example:
!
!    NBIN = 7, A = 10, B = 20
!
!    BIN      CMIN  CMAX
!
!    1         -HUGE 10
!    2         10    12
!    3         12    14
!    4         14    16
!    5         16    18
!    6         18    20
!    7         20    HUGE
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NBIN, the number of bins.  NBIN is normally
!    at least 3.  If NBIN is 1 or 2, then everything is assigned to bin 1.
!
!    Input, integer ( kind = 4 ) BIN, the index of the bin to be considered.
!    If BIN is less than 1, or greater than NBIN, the user will get what
!    the user deserves.
!
!    Input, real ( kind = 8 ) A, B, the lower and upper limits of the
!    bin interval.  While A is expected to be less than B, the code
!    should return useful results if A is actually greater than B.
!
!    Output, real ( kind = 8 ) CMIN, CMAX, the minimum and maximum
!    limits on the bin.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) bin
  real ( kind = 8 ) cmax
  real ( kind = 8 ) cmin
  integer ( kind = 4 ) nbin
!
!  Take care of special cases.
!
  if ( nbin <= 2 ) then
    cmin = -huge ( cmin )
    cmax =  huge ( cmax )
    return
  end if

  if ( b == a ) then
    cmin = -huge ( cmin )
    cmax =  huge ( cmax )
    return
  end if
!
!  Compute the bin limits.
!
  if ( bin == 1 ) then
    cmin = -huge ( cmin )
    cmax = a
  else if ( bin < nbin ) then

    cmin = ( real ( nbin - bin,     kind = 8 ) * a   &
           + real (        bin - 2, kind = 8 ) * b ) &
           / real ( nbin       - 2, kind = 8 )

    cmax = ( real ( nbin - bin - 1, kind = 8 ) * a   &
           + real (        bin - 1, kind = 8 ) * b ) &
           / real ( nbin       - 2, kind = 8 )

  else if ( bin == nbin ) then
    cmin = b
    cmax = huge ( cmax )
  else
    cmin = -huge ( cmin )
    cmax =  huge ( cmax )
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
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NBIN, the number of bins.
!
!    Input, integer ( kind = 4 ) BIN, the index of the bin to be considered.
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
  integer ( kind = 4 ) bin
  real ( kind = 8 ) cmax
  real ( kind = 8 ) cmin
  integer ( kind = 4 ) nbin
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
subroutine bin_to_r82_even ( nbin, bin, a, b, cmin, cmax )

!*****************************************************************************80
!
!! BIN_TO_R82_EVEN returns the limits for a given R82 "bin" in [A,B].
!
!  Discussion:
!
!    The interval from A to B is divided into NBIN-2 equal subintervals or bins.
!    An initial bin takes everything less than A, and a final bin takes
!    everything greater than B.
!
!  Example:
!
!    NBIN = 7, A(1) = 5, B(1) = 15
!              A(2) = 0, B(2) = 20
!
!     BIN         CMIN      CMAX
!    ------   -----------  --------
!    1, 1     -HUGE -HUGE   5     0
!    2, 2       5     0     7     4
!    3, 3       7     4     9     8
!    4, 4       9     8    11    12
!    5, 5      11    12    13    16
!    6, 6      13    16    15    20
!    7, 7      15    20    HUGE HUGE
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NBIN, the number of bins.  NBIN is normally
!    at least 3.  If NBIN is 1 or 2, then everything is assigned to bin 1.
!
!    Input, integer ( kind = 4 ) BIN(2), the index of the bin to be considered.
!    If BIN(I) is less than 1, or greater than NBIN, the user will get what
!    the user deserves.
!
!    Input, real ( kind = 8 ) A(2), B(2), the lower and upper limits of the
!    bin interval.  While A(I) is expected to be less than B(I), the code
!    should return useful results if A(I) is actually greater than B(I).
!
!    Output, real ( kind = 8 ) CMIN(2), CMAX(2), the minimum and maximum
!    limits on the bin.
!
  implicit none

  integer ( kind = 4 ), parameter :: ndim = 2

  real ( kind = 8 ) a(ndim)
  real ( kind = 8 ) b(ndim)
  integer ( kind = 4 ) bin(ndim)
  real ( kind = 8 ) cmax(ndim)
  real ( kind = 8 ) cmin(ndim)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) nbin

  do i = 1, ndim
    call bin_to_r8_even ( nbin, bin(i), a(i), b(i), cmin(i), cmax(i) )
  end do

  return
end
subroutine bin_to_r82_even2 ( nbin, bin, a, b, cmin, cmax )

!*****************************************************************************80
!
!! BIN_TO_R82_EVEN2 returns the limits for a given R82 "bin" in [A,B].
!
!  Discussion:
!
!    The interval from A to B is divided into NBIN equal subintervals or bins.
!
!  Example:
!
!    NBIN = 5, A(1) = 5, B(1) = 15
!              A(2) = 0, B(2) = 20
!
!     BIN         CMIN      CMAX
!    ------   -----------  --------
!    1, 1       5     0     7     4
!    2, 2       7     4     9     8
!    3, 3       9     8    11    12
!    4, 4      11    12    13    16
!    5, 5      13    16    15    20
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NBIN, the number of bins.
!
!    Input, integer ( kind = 4 ) BIN(2), the index of the bin to be considered.
!
!    Input, real ( kind = 8 ) A(2), B(2), the lower and upper limits of the
!    bin interval.  While A(I) is expected to be less than B(I), the code
!    should return useful results if A(I) is actually greater than B(I).
!
!    Output, real ( kind = 8 ) CMIN(2), CMAX(2), the minimum and maximum
!    limits on the bin.
!
  implicit none

  integer ( kind = 4 ), parameter :: ndim = 2

  real ( kind = 8 ) a(ndim)
  real ( kind = 8 ) b(ndim)
  integer ( kind = 4 ) bin(ndim)
  real ( kind = 8 ) cmax(ndim)
  real ( kind = 8 ) cmin(ndim)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) nbin

  do i = 1, ndim
    call bin_to_r8_even2 ( nbin, bin(i), a(i), b(i), cmin(i), cmax(i) )
  end do

  return
end
subroutine bin_to_r82_even3 ( nbin, bin, a, b, cmin, cmax )

!*****************************************************************************80
!
!! BIN_TO_R82_EVEN3 returns the limits for a given R82 "bin" in [A,B].
!
!  Discussion:
!
!    The interval from A(I) to B(I) is divided into NBIN(I) equal
!    subintervals or bins.
!
!  Example:
!
!    NBIN = (/ 4, 5, /)
!
!    A(1) = 5, B(1) = 15
!    A(2) = 0, B(2) = 20
!
!     BIN         CMIN      CMAX
!    ------   -----------  --------
!    1, 1       5     0     7     4
!    2, 2       7     4     9     8
!    3, 3       9     8    11    12
!    4, 4      11    12    13    16
!    5, 5      13    16    15    20
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NBIN(2), the number of bins in each dimension.
!
!    Input, integer ( kind = 4 ) BIN(2), the index of the bin to be considered.
!
!    Input, real ( kind = 8 ) A(2), B(2), the lower and upper limits of the
!    bin interval.  While A(I) is expected to be less than B(I), the code
!    should return useful results if A(I) is actually greater than B(I).
!
!    Output, real ( kind = 8 ) CMIN(2), CMAX(2), the minimum and maximum
!    limits on the bin.
!
  implicit none

  integer ( kind = 4 ), parameter :: ndim = 2

  real ( kind = 8 ) a(ndim)
  real ( kind = 8 ) b(ndim)
  integer ( kind = 4 ) bin(ndim)
  real ( kind = 8 ) cmax(ndim)
  real ( kind = 8 ) cmin(ndim)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) nbin(ndim)

  do i = 1, ndim
    call bin_to_r8_even2 ( nbin(i), bin(i), a(i), b(i), cmin(i), cmax(i) )
  end do

  return
end
subroutine bin_to_r83_even2 ( nbin, bin, a, b, cmin, cmax )

!*****************************************************************************80
!
!! BIN_TO_R83_EVEN2 returns the limits for a given R83 "bin" in [A,B].
!
!  Discussion:
!
!    The interval from A to B is divided into NBIN equal subintervals or bins.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NBIN, the number of bins.
!
!    Input, integer ( kind = 4 ) BIN(), the index of the bin to be considered.
!
!    Input, real ( kind = 8 ) A(3), B(3), the lower and upper limits of the
!    bin interval.  While A(I) is expected to be less than B(I), the code
!    should return useful results if A(I) is actually greater than B(I).
!
!    Output, real ( kind = 8 ) CMIN(3), CMAX(3), the minimum and maximum
!    limits on the bin.
!
  implicit none

  integer ( kind = 4 ), parameter :: ndim = 3

  real ( kind = 8 ) a(ndim)
  real ( kind = 8 ) b(ndim)
  integer ( kind = 4 ) bin(ndim)
  real ( kind = 8 ) cmax(ndim)
  real ( kind = 8 ) cmin(ndim)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) nbin

  do i = 1, ndim
    call bin_to_r8_even2 ( nbin, bin(i), a(i), b(i), cmin(i), cmax(i) )
  end do

  return
end
subroutine bin_to_r83_even3 ( nbin, bin, a, b, cmin, cmax )

!*****************************************************************************80
!
!! BIN_TO_R83_EVEN3 returns the limits for a given R83 "bin" in [A,B].
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
!    30 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NBIN(3), the number of bins in each dimension.
!
!    Input, integer ( kind = 4 ) BIN(3), the index of the bin to be considered.
!
!    Input, real ( kind = 8 ) A(3), B(3), the lower and upper limits of the
!    bin interval.  While A(I) is expected to be less than B(I), the code
!    should return useful results if A(I) is actually greater than B(I).
!
!    Output, real ( kind = 8 ) CMIN(3), CMAX(3), the minimum and maximum
!    limits on the bin.
!
  implicit none

  integer ( kind = 4 ), parameter :: ndim = 3

  real ( kind = 8 ) a(ndim)
  real ( kind = 8 ) b(ndim)
  integer ( kind = 4 ) bin(ndim)
  real ( kind = 8 ) cmax(ndim)
  real ( kind = 8 ) cmin(ndim)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) nbin(ndim)

  do i = 1, ndim
    call bin_to_r8_even2 ( nbin(i), bin(i), a(i), b(i), cmin(i), cmax(i) )
  end do

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
!    08 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NDIM, the dimension of the space.
!
!    Input, integer ( kind = 4 ) NBIN(NDIM), the number of bins in 
!    each dimension.
!
!    Input, integer ( kind = 4 ) BIN(NDIM), the index of the bin to be
!    considered.
!
!    Input, real ( kind = 8 ) A(NDIM), B(NDIM), the lower and upper limits
!    of the bin interval.  While A(I) is expected to be less than B(I),
!    the code should return useful results if A(I) is actually greater
!    than B(I).
!
!    Output, real ( kind = 8 ) CMIN(NDIM), CMAX(NDIM), the minimum and
!    maximum limits on the bin.
!
  implicit none

  integer ( kind = 4 ) ndim

  real ( kind = 8 ) a(ndim)
  real ( kind = 8 ) b(ndim)
  integer ( kind = 4 ) bin(ndim)
  real ( kind = 8 ) cmax(ndim)
  real ( kind = 8 ) cmin(ndim)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) nbin(ndim)

  do i = 1, ndim
    call bin_to_r8_even2 ( nbin(i), bin(i), a(i), b(i), cmin(i), cmax(i) )
  end do

  return
end
function diaedg ( x0, y0, x1, y1, x2, y2, x3, y3 )

!*****************************************************************************80
!
!! DIAEDG chooses a diagonal edge.
!
!  Discussion:
!
!    The routine determines whether 0--2 or 1--3 is the diagonal edge
!    that should be chosen, based on the circumcircle criterion, where
!    (X0,Y0), (X1,Y1), (X2,Y2), (X3,Y3) are the vertices of a simple
!    quadrilateral in counterclockwise order.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 February 2001
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X0, Y0, X1, Y1, X2, Y2, X3, Y3, the
!    coordinates of the vertices of a quadrilateral, given in
!    counter clockwise order.
!
!    Output, integer ( kind = 4 ) DIAEDG, chooses a diagonal:
!    +1, if diagonal edge 02 is chosen;
!    -1, if diagonal edge 13 is chosen;
!     0, if the four vertices are cocircular.
!
  implicit none

  real ( kind = 8 ) ca
  real ( kind = 8 ) cb
  integer ( kind = 4 ) diaedg
  real ( kind = 8 ) dx10
  real ( kind = 8 ) dx12
  real ( kind = 8 ) dx30
  real ( kind = 8 ) dx32
  real ( kind = 8 ) dy10
  real ( kind = 8 ) dy12
  real ( kind = 8 ) dy30
  real ( kind = 8 ) dy32
  real ( kind = 8 ) s
  real ( kind = 8 ) tol
  real ( kind = 8 ) tola
  real ( kind = 8 ) tolb
  real ( kind = 8 ) x0
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ) y0
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) y3

  tol = 100.0D+00 * epsilon ( tol )

  dx10 = x1 - x0
  dy10 = y1 - y0
  dx12 = x1 - x2
  dy12 = y1 - y2
  dx30 = x3 - x0
  dy30 = y3 - y0
  dx32 = x3 - x2
  dy32 = y3 - y2

  tola = tol * max ( abs ( dx10 ), abs ( dy10 ), abs ( dx30 ), abs ( dy30 ) )
  tolb = tol * max ( abs ( dx12 ), abs ( dy12 ), abs ( dx32 ), abs ( dy32 ) )

  ca = dx10 * dx30 + dy10 * dy30
  cb = dx12 * dx32 + dy12 * dy32

  if ( tola < ca .and. tolb < cb ) then

    diaedg = -1

  else if ( ca < -tola .and. cb < -tolb ) then

    diaedg = 1

  else

    tola = max ( tola, tolb )
    s = ( dx10 * dy30 - dx30 * dy10 ) * cb + ( dx32 * dy12 - dx12 * dy32 ) * ca

    if ( tola < s ) then
      diaedg = -1
    else if ( s < -tola ) then
      diaedg = 1
    else
      diaedg = 0
    end if

  end if

  return
end
subroutine dtris2 ( point_num, point_xy, tri_num, tri_vert, tri_nabe )

!*****************************************************************************80
!
!! DTRIS2 constructs a Delaunay triangulation of 2D vertices.
!
!  Discussion:
!
!    The routine constructs the Delaunay triangulation of a set of 2D vertices
!    using an incremental approach and diagonal edge swaps.  Vertices are
!    first sorted in lexicographically increasing (X,Y) order, and
!    then are inserted one at a time from outside the convex hull.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 August 2001
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of vertices.
!
!    Input/output, real ( kind = 8 ) POINT_XY(2,POINT_NUM), the coordinates
!    of the vertices.  On output, the vertices have been sorted into
!    dictionary order.
!
!    Output, integer ( kind = 4 ) TRI_NUM, the number of triangles in the
!    triangulation; TRI_NUM is equal to 2*POINT_NUM - NB - 2, where NB is
!    the number of boundary vertices.
!
!    Output, integer ( kind = 4 ) TRI_VERT(3,TRI_NUM), the nodes that make up 
!    each triangle.  The elements are indices of POINT_XY.  The vertices of 
!    the triangles are in counter clockwise order.
!
!    Output, integer ( kind = 4 ) TRI_NABE(3,TRI_NUM), the triangle neighbor 
!    list.  Positive elements are indices of TIL; negative elements are used 
!    for links of a counter clockwise linked list of boundary edges; 
!    LINK = -(3*I + J-1) where I, J = triangle, edge index; TRI_NABE(J,I) refers
!    to the neighbor along edge from vertex J to J+1 (mod 3).
!
  implicit none

  integer ( kind = 4 ) point_num

  real ( kind = 8 ) cmax
  integer ( kind = 4 ) e
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) indx(point_num)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) ledg
  integer ( kind = 4 ) lr
  integer ( kind = 4 ) lrline
  integer ( kind = 4 ) ltri
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) n
  real ( kind = 8 ) point_xy(2,point_num)
  integer ( kind = 4 ) redg
  integer ( kind = 4 ) rtri
  integer ( kind = 4 ) stack(point_num)
  integer ( kind = 4 ) t
  real ( kind = 8 ) tol
  integer ( kind = 4 ) top
  integer ( kind = 4 ) tri_nabe(3,point_num*2)
  integer ( kind = 4 ) tri_num
  integer ( kind = 4 ) tri_vert(3,point_num*2)

  tol = 100.0D+00 * epsilon ( tol )

  ierr = 0
!
!  Sort the vertices by increasing (x,y).
!
  call r82vec_sort_heap_index_a ( point_num, point_xy, indx )

  call r82vec_permute ( point_num, point_xy, indx )
!
!  Make sure that the data points are "reasonably" distinct.
!
  m1 = 1

  do i = 2, point_num

    m = m1
    m1 = i

    k = 0

    do j = 1, 2

      cmax = max ( abs ( point_xy(j,m) ), abs ( point_xy(j,m1) ) )

      if ( tol * ( cmax + 1.0D+00 ) &
           < abs ( point_xy(j,m) - point_xy(j,m1) ) ) then
        k = j
        exit
      end if

    end do

    if ( k == 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DTRIS2 - Fatal error!'
      write ( *, '(a,i6)' ) '  Fails for point number I = ', i
      write ( *, '(a,i6)' ) '  M = ', m
      write ( *, '(a,i6)' ) '  M1 = ', m1
      write ( *, '(a,2g14.6)' ) '  X,Y(M)  = ', point_xy(1,m), point_xy(2,m)
      write ( *, '(a,2g14.6)' ) '  X,Y(M1) = ', point_xy(1,m1), point_xy(2,m1)
      ierr = 224
      return
    end if

  end do
!
!  Starting from points M1 and M2, search for a third point M that
!  makes a "healthy" triangle (M1,M2,M)
!
  m1 = 1
  m2 = 2
  j = 3

  do

    if ( point_num < j ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DTRIS2 - Fatal error!'
      ierr = 225
      return
    end if

    m = j

    lr = lrline ( point_xy(1,m), point_xy(2,m), point_xy(1,m1), &
      point_xy(2,m1), point_xy(1,m2), point_xy(2,m2), 0.0D+00 )

    if ( lr /= 0 ) then
      exit
    end if

    j = j + 1

  end do
!
!  Set up the triangle information for (M1,M2,M), and for any other
!  triangles you created because points were collinear with M1, M2.
!
  tri_num = j - 2

  if ( lr == -1 ) then

    tri_vert(1,1) = m1
    tri_vert(2,1) = m2
    tri_vert(3,1) = m
    tri_nabe(3,1) = -3

    do i = 2, tri_num

      m1 = m2
      m2 = i+1
      tri_vert(1,i) = m1
      tri_vert(2,i) = m2
      tri_vert(3,i) = m
      tri_nabe(1,i-1) = -3 * i
      tri_nabe(2,i-1) = i
      tri_nabe(3,i) = i - 1

    end do

    tri_nabe(1,tri_num) = -3 * tri_num - 1
    tri_nabe(2,tri_num) = -5
    ledg = 2
    ltri = tri_num

  else

    tri_vert(1,1) = m2
    tri_vert(2,1) = m1
    tri_vert(3,1) = m
    tri_nabe(1,1) = -4

    do i = 2, tri_num
      m1 = m2
      m2 = i+1
      tri_vert(1,i) = m2
      tri_vert(2,i) = m1
      tri_vert(3,i) = m
      tri_nabe(3,i-1) = i
      tri_nabe(1,i) = -3 * i - 3
      tri_nabe(2,i) = i - 1
    end do

    tri_nabe(3,tri_num) = -3 * tri_num
    tri_nabe(2,1) = -3 * tri_num - 2
    ledg = 2
    ltri = 1

  end if
!
!  Insert the vertices one at a time from outside the convex hull,
!  determine visible boundary edges, and apply diagonal edge swaps until
!  Delaunay triangulation of vertices (so far) is obtained.
!
  top = 0

  do i = j+1, point_num

    m = i
    m1 = tri_vert(ledg,ltri)

    if ( ledg <= 2 ) then
      m2 = tri_vert(ledg+1,ltri)
    else
      m2 = tri_vert(1,ltri)
    end if

    lr = lrline ( point_xy(1,m), point_xy(2,m), point_xy(1,m1), &
      point_xy(2,m1), point_xy(1,m2), point_xy(2,m2), 0.0D+00 )

    if ( 0 < lr ) then
      rtri = ltri
      redg = ledg
      ltri = 0
    else
      l = -tri_nabe(ledg,ltri)
      rtri = l / 3
      redg = mod(l,3) + 1
    end if

    call vbedg ( point_xy(1,m), point_xy(2,m), point_num, point_xy, tri_num, &
      tri_vert, tri_nabe, ltri, ledg, rtri, redg )

    n = tri_num + 1
    l = -tri_nabe(ledg,ltri)

    do

      t = l / 3
      e = mod ( l, 3 ) + 1
      l = -tri_nabe(e,t)
      m2 = tri_vert(e,t)

      if ( e <= 2 ) then
        m1 = tri_vert(e+1,t)
      else
        m1 = tri_vert(1,t)
      end if

      tri_num = tri_num + 1
      tri_nabe(e,t) = tri_num
      tri_vert(1,tri_num) = m1
      tri_vert(2,tri_num) = m2
      tri_vert(3,tri_num) = m
      tri_nabe(1,tri_num) = t
      tri_nabe(2,tri_num) = tri_num - 1
      tri_nabe(3,tri_num) = tri_num + 1
      top = top + 1

      if ( point_num < top ) then
        ierr = 8
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'DTRIS2 - Fatal error!'
        write ( *, '(a)' ) '  Stack overflow.'
        return
      end if

      stack(top) = tri_num

      if ( t == rtri .and. e == redg ) then
        exit
      end if

    end do

    tri_nabe(ledg,ltri) = -3 * n - 1
    tri_nabe(2,n) = -3 * tri_num - 2
    tri_nabe(3,tri_num) = -l
    ltri = n
    ledg = 2

    call swapec ( m, top, ltri, ledg, point_num, point_xy, tri_num, &
      tri_vert, tri_nabe, stack, ierr )

    if ( ierr /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DTRIS2 - Fatal error!'
      write ( *, '(a)' ) '  Error return from SWAPEC.'
      return
    end if

  end do
!
!  Now account for the sorting that we did.
!
  do i = 1, 3
    do j = 1, tri_num
      tri_vert(i,j) = indx ( tri_vert(i,j) )
    end do
  end do

  call perm_inv ( point_num, indx )

  call r82vec_permute ( point_num, point_xy, indx )

  return
end
subroutine get_seed ( seed )

!*****************************************************************************80
!
!! GET_SEED returns a seed for the random number generator.
!
!  Discussion:
!
!    The seed depends on the current time, and ought to be (slightly)
!    different every millisecond.  Once the seed is obtained, a random
!    number generator should be called a few times to further process
!    the seed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 November 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) SEED, a pseudorandom seed value.
!
  implicit none

  integer ( kind = 4 ) seed
  real ( kind = 8 ) temp
  character ( len = 10 ) time
  character ( len = 8 ) today
  integer ( kind = 4 ) values(8)
  character ( len = 5 ) zone

  call date_and_time ( today, time, zone, values )

  temp = 0.0D+00

  temp = temp + real ( values(2) - 1, kind = 8 ) / 11.0D+00
  temp = temp + real ( values(3) - 1, kind = 8 ) / 30.0D+00
  temp = temp + real ( values(5),     kind = 8 ) / 23.0D+00
  temp = temp + real ( values(6),     kind = 8 ) / 59.0D+00
  temp = temp + real ( values(7),     kind = 8 ) / 59.0D+00
  temp = temp + real ( values(8),     kind = 8 ) / 999.0D+00
  temp = temp / 6.0D+00
!
!  Force 0 < TEMP <= 1.
!
  do while ( temp <= 0.0D+00 )
    temp = temp + 1.0D+00
  end do

  do while ( 1.0D+00 < temp )
    temp = temp - 1.0D+00
  end do

  seed = int ( real ( huge ( 1 ), kind = 8 ) * temp )
!
!  Never use a seed of 0 or maximum integer.
!
  if ( seed == 0 ) then
    seed = 1
  end if

  if ( seed == huge ( 1 ) ) then
    seed = seed - 1
  end if

  return
end
function i4_modp ( i, j )

!*****************************************************************************80
!
!! I4_MODP returns the nonnegative remainder of integer division.
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
!        I     J     MOD  I4_MODP    Factorization
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

  if ( j == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_MODP - Fatal error!'
    write ( *, '(a,i6)' ) '  I4_MODP ( I, J ) called with J = ', j
    stop
  end if

  i4_modp = mod ( i, j )

  if ( i4_modp < 0 ) then
    i4_modp = i4_modp + abs ( j )
  end if

  return
end
subroutine i4_swap ( i, j )

!*****************************************************************************80
!
!! I4_SWAP switches two integer values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) I, J.  On output, the values of I and
!    J have been interchanged.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k

  k = i
  i = j
  j = k

  return
end
function i4_wrap ( ival, ilo, ihi )

!*****************************************************************************80
!
!! I4_WRAP forces an integer to lie between given limits by wrapping.
!
!  Example:
!
!    ILO = 4, IHI = 8
!
!    I  I4_WRAP
!
!    -2     8
!    -1     4
!     0     5
!     1     6
!     2     7
!     3     8
!     4     4
!     5     5
!     6     6
!     7     7
!     8     8
!     9     4
!    10     5
!    11     6
!    12     7
!    13     8
!    14     4
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 August 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IVAL, an integer value.
!
!    Input, integer ( kind = 4 ) ILO, IHI, the desired bounds for the value.
!
!    Output, integer ( kind = 4 ) I4_WRAP, a "wrapped" version of IVAL.
!
  implicit none

  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) wide

  jlo = min ( ilo, ihi )
  jhi = max ( ilo, ihi )

  wide = jhi - jlo + 1

  if ( wide == 1 ) then
    i4_wrap = jlo
  else
    i4_wrap = jlo + i4_modp ( ival - jlo, wide )
  end if

  return
end
subroutine i4mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! I4MAT_PRINT prints an integer matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 June 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, integer ( kind = 4 ) A(M,N), the matrix to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * ) title

  ilo = 1
  ihi = m
  jlo = 1
  jhi = n

  call i4mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

  return
end
subroutine i4mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! I4MAT_PRINT_SOME prints some of an integer matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 November 2003
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
!    Input, integer ( kind = 4 ) ILO, JLO, the first row and column to print.
!
!    Input, integer ( kind = 4 ) IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 10
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  character ( len = 7 ) ctemp(incx)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i7)') j
    end do

    write ( *, '(''  Col '',10a7)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) ' '

    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi

      do j2 = 1, inc

        j = j2lo - 1 + j2

        write ( ctemp(j2), '(i7)' ) a(i,j)

      end do

      write ( *, '(i5,1x,10a7)' ) i, ( ctemp(j), j = 1, inc )

    end do

  end do

  return
end
subroutine i4mat_transpose_print ( m, n, a, title )

!*****************************************************************************80
!
!! I4MAT_TRANSPOSE_PRINT prints an IMAT, transposed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 December 2004
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
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  character ( len = * ) title

  call i4mat_transpose_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine i4mat_transpose_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! I4MAT_TRANSPOSE_PRINT_SOME prints some of the transpose of an integer matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 December 2004
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
!    Input, integer ( kind = 4 ) ILO, JLO, the first row and column to print.
!
!    Input, integer ( kind = 4 ) IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 10
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  character ( len = 7 ) ctemp(incx)
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
      write ( ctemp(i2), '(i7)') i
    end do

    write ( *, '(''  Row '',10a7)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Col'
    write ( *, '(a)' ) ' '

    j2lo = max ( jlo, 1 )
    j2hi = min ( jhi, m )

    do j = j2lo, j2hi

      do i2 = 1, inc

        i = i2lo - 1 + i2

        write ( ctemp(i2), '(i7)' ) a(i,j)

      end do

      write ( *, '(i5,1x,10a7)' ) j, ( ctemp(i), i = 1, inc )

    end do

  end do

  return
end
subroutine i4vec_heap_d ( n, a )

!*****************************************************************************80
!
!! I4VEC_HEAP_D reorders an I4VEC into a descending heap.
!
!  Discussion:
!
!    A descending heap is an array A with the property that, for every index J,
!    A(2*J) <= A(J) and A(2*J+1) <= A(J), (as long as the indices
!    2*J and 2*J+1 are legal).
!
!  Diagram:
!
!                  A(1)
!                /      \
!            A(2)         A(3)
!          /     \        /  \
!      A(4)       A(5)  A(6) A(7)
!      /  \       /   \
!    A(8) A(9) A(10) A(11)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the input array.
!
!    Input/output, integer ( kind = 4 ) A(N).
!    On input, an unsorted array.
!    On output, the array has been reordered into a heap.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ifree
  integer ( kind = 4 ) key
  integer ( kind = 4 ) m
!
!  Only nodes N/2 down to 1 can be "parent" nodes.
!
  do i = n/2, 1, -1
!
!  Copy the value out of the parent node.
!  Position IFREE is now "open".
!
    key = a(i)
    ifree = i

    do
!
!  Positions 2*IFREE and 2*IFREE + 1 are the descendants of position
!  IFREE.  (One or both may not exist because they exceed N.)
!
      m = 2 * ifree
!
!  Does the first position exist?
!
      if ( n < m ) then
        exit
      end if
!
!  Does the second position exist?
!
      if ( m + 1 <= n ) then
!
!  If both positions exist, take the larger of the two values,
!  and update M if necessary.
!
        if ( a(m) < a(m+1) ) then
          m = m + 1
        end if

      end if
!
!  If the large descendant is larger than KEY, move it up,
!  and update IFREE, the location of the free position, and
!  consider the descendants of THIS position.
!
      if ( a(m) <= key ) then
        exit
      end if

      a(ifree) = a(m)
      ifree = m

    end do
!
!  Once there is no more shifting to do, KEY moves into the free spot IFREE.
!
    a(ifree) = key

  end do

  return
end
subroutine i4vec_indicator ( n, a )

!*****************************************************************************80
!
!! I4VEC_INDICATOR sets an I4VEC to the indicator vector.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Output, integer ( kind = 4 ) A(N), the array to be initialized.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i

  do i = 1, n
    a(i) = i
  end do

  return
end
subroutine i4vec_print ( n, a, title )

!*****************************************************************************80
!
!! I4VEC_PRINT prints an I4VEC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, integer ( kind = 4 ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) big
  integer ( kind = 4 ) i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  big = maxval ( abs ( a(1:n) ) )

  write ( *, '(a)' ) ' '
  if ( big < 1000 ) then
    do i = 1, n
      write ( *, '(2x,i6,2x,i4)' ) i, a(i)
    end do
  else if ( big < 1000000 ) then
    do i = 1, n
      write ( *, '(2x,i6,2x,i7)' ) i, a(i)
    end do
  else
    do i = 1, n
      write ( *, '(2x,i6,2x,i12)' ) i, a(i)
    end do
  end if

  return
end
subroutine i4vec_sort_heap_a ( n, a )

!*****************************************************************************80
!
!! I4VEC_SORT_HEAP_A ascending sorts an I4VEC using heap sort.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input/output, integer ( kind = 4 ) A(N).
!    On input, the array to be sorted;
!    On output, the array has been sorted.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) n1

  if ( n <= 1 ) then
    return
  end if
!
!  1: Put A into descending heap form.
!
  call i4vec_heap_d ( n, a )
!
!  2: Sort A.
!
!  The largest object in the heap is in A(1).
!  Move it to position A(N).
!
  call i4_swap ( a(1), a(n) )
!
!  Consider the diminished heap of size N1.
!
  do n1 = n-1, 2, -1
!
!  Restore the heap structure of A(1) through A(N1).
!
    call i4vec_heap_d ( n1, a )
!
!  Take the largest object from A(1) and move it to A(N1).
!
    call i4_swap ( a(1), a(n1) )

  end do

  return
end
subroutine i4vec_sorted_unique ( n, a, nuniq )

!*****************************************************************************80
!
!! I4VEC_SORTED_UNIQUE finds the number of unique elements in a sorted I4VEC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements in A.
!
!    Input/output, integer ( kind = 4 ) A(N).  On input, the sorted
!    integer ( kind = 4 ) array.  On output, the unique elements in A.
!
!    Output, integer ( kind = 4 ) NUNIQ, the number of unique elements in A.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) itest
  integer ( kind = 4 ) nuniq

  nuniq = 0

  if ( n <= 0 ) then
    return
  end if

  nuniq = 1

  do itest = 2, n

    if ( a(itest) /= a(nuniq) ) then
      nuniq = nuniq + 1
      a(nuniq) = a(itest)
    end if

  end do

  return
end
subroutine i4vec2_compare ( n, a1, a2, i, j, isgn )

!*****************************************************************************80
!
!! I4VEC2_COMPARE compares pairs of integers stored in two vectors.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 October 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data items.
!
!    Input, integer ( kind = 4 ) A1(N), A2(N), contain the two components
!    of each item.
!
!    Input, integer ( kind = 4 ) I, J, the items to be compared.
!
!    Output, integer ( kind = 4 ) ISGN, the results of the comparison:
!    -1, item I < item J,
!     0, item I = item J,
!    +1, item J < item I.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a1(n)
  integer ( kind = 4 ) a2(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j

  isgn = 0

       if ( a1(i) < a1(j) ) then

    isgn = -1

  else if ( a1(i) == a1(j) ) then

         if ( a2(i) < a2(j) ) then
      isgn = -1
    else if ( a2(i) < a2(j) ) then
      isgn = 0
    else if ( a2(j) < a2(i) ) then
      isgn = +1
    end if

  else if ( a1(j) < a1(i) ) then

    isgn = +1

  end if

  return
end
subroutine i4vec2_sort_a ( n, a1, a2 )

!*****************************************************************************80
!
!! I4VEC2_SORT_A ascending sorts a vector of pairs of integers.
!
!  Discussion:
!
!    Each item to be sorted is a pair of integers (I,J), with the I
!    and J values stored in separate vectors A1 and A2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items of data.
!
!    Input/output, integer ( kind = 4 ) A1(N), A2(N), the data to be sorted..
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a1(n)
  integer ( kind = 4 ) a2(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
!
!  Initialize.
!
  i = 0
  indx = 0
  isgn = 0
  j = 0
!
!  Call the external heap sorter.
!
  do

    call sort_heap_external ( n, indx, i, j, isgn )
!
!  Interchange the I and J objects.
!
    if ( 0 < indx ) then

      call i4_swap ( a1(i), a1(j) )
      call i4_swap ( a2(i), a2(j) )
!
!  Compare the I and J objects.
!
    else if ( indx < 0 ) then

      call i4vec2_compare ( n, a1, a2, i, j, isgn )

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
subroutine i4vec2_sorted_unique ( n, a1, a2, nuniq )

!*****************************************************************************80
!
!! I4VEC2_SORTED_UNIQUE keeps the unique elements in a array of pairs of integers.
!
!  Discussion:
!
!    Item I is stored as the pair A1(I), A2(I).
!
!    The items must have been sorted, or at least it must be the
!    case that equal items are stored in adjacent vector locations.
!
!    If the items were not sorted, then this routine will only
!    replace a string of equal values by a single representative.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items.
!
!    Input/output, integer ( kind = 4 ) A1(N), A2(N).
!    On input, the array of N items.
!    On output, an array of NUNIQ unique items.
!
!    Output, integer ( kind = 4 ) NUNIQ, the number of unique items.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a1(n)
  integer ( kind = 4 ) a2(n)
  integer ( kind = 4 ) itest
  integer ( kind = 4 ) nuniq

  nuniq = 0

  if ( n <= 0 ) then
    return
  end if

  nuniq = 1

  do itest = 2, n

    if ( a1(itest) /= a1(nuniq) .or. a2(itest) /= a2(nuniq) ) then

      nuniq = nuniq + 1

      a1(nuniq) = a1(itest)
      a2(nuniq) = a2(itest)

    end if

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
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N1, N2, the half-widths of the box, that is, the
!    maximum distance allowed between (IC,JC) and (I,J).
!
!    Input, integer ( kind = 4 ) IC, JC, the central cell of the box.
!
!    Input/output, integer ( kind = 4 ) I, J.  On input, the previous index set.
!    On output, the next index set.  On the first call, MORE should
!    be set to FALSE, and the input values of I and J are ignored.
!
!    Input/output, logical MORE.
!    On the first call for a given box, the user should set MORE to FALSE.
!    On return, the routine sets MORE to TRUE.
!    When there are no more indices, the routine sets MORE to FALSE.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jc
  logical more
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2

  if ( .not. more ) then
    more = .true.
    i = ic - n1
    j = jc - n2
    return
  end if

  if ( i == ic + n1 .and. &
       j == jc + n2 ) then
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
  else if ( j < jc + n2 .and. &
    ( i == ic - n1 .or. i == ic + n1 ) ) then
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
!    (N1,N2,N3).  The indices are exactly those between (IC-N1,JC-N2,KC-N3) and
!    (IC+N1,JC+N2,KC+N3) with the property that at least one of I, J, and K
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
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N1, N2, N3, the "half widths" of the box,
!    that is, the maximum distances from the central cell allowed for I, J and K.
!
!    Input, integer ( kind = 4 ) IC, JC, KC, the central cell of the box.
!
!    Input/output, integer ( kind = 4 ) I, J, K.  On input, the previous
!    index set.  On output, the next index set.  On the first call, MORE should
!    be set to FALSE, and the input values of I, J, and K are ignored.
!
!    Input/output, logical MORE.
!    On the first call for a given box, the user should set MORE to FALSE.
!    On return, the routine sets MORE to TRUE.
!    When there are no more indices, the routine sets MORE to FALSE.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jc
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kc
  logical more
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n3

  if ( .not. more ) then
    more = .true.
    i = ic - n1
    j = jc - n2
    k = kc - n3
    return
  end if

  if ( i == ic + n1 .and. &
       j == jc + n2 .and. &
       k == kc + n3 ) then
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
function lrline ( xu, yu, xv1, yv1, xv2, yv2, dv )

!*****************************************************************************80
!
!! LRLINE determines where a point lies in relation to a directed line.
!
!  Discussion:
!
!    LRLINE determines whether a point is to the left of, right of,
!    or on a directed line parallel to a line through given points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 June 2001
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XU, YU, XV1, YV1, XV2, YV2, are vertex
!    coordinates; the directed line is parallel to and at signed distance
!    DV to the left of the directed line from (XV1,YV1) to (XV2,YV2); (XU,YU)
!    is the vertex for which the position relative to the directed line is
!    to be determined.
!
!    Input, real ( kind = 8 ) DV, the signed distance, positive for left.
!
!    Output, integer ( kind = 4 ) LRLINE, is +1, 0, or -1 depending on
!    whether (XU,YU) is to the right of, on, or left of the directed line.
!    LRLINE is 0 if the line degenerates to a point.
!
  implicit none

  real ( kind = 8 ) dv
  real ( kind = 8 ) dx
  real ( kind = 8 ) dxu
  real ( kind = 8 ) dy
  real ( kind = 8 ) dyu
  integer ( kind = 4 ) lrline
  real ( kind = 8 ) t
  real ( kind = 8 ), parameter :: tol = 0.0000001D+00
  real ( kind = 8 ) tolabs
  real ( kind = 8 ) xu
  real ( kind = 8 ) xv1
  real ( kind = 8 ) xv2
  real ( kind = 8 ) yu
  real ( kind = 8 ) yv1
  real ( kind = 8 ) yv2

  dx = xv2 - xv1
  dy = yv2 - yv1
  dxu = xu - xv1
  dyu = yu - yv1

  tolabs = tol * max ( abs ( dx ), abs ( dy ), abs ( dxu ), abs ( dyu ), &
    abs ( dv ) )

  t = dy * dxu - dx * dyu + dv * sqrt ( dx * dx + dy * dy )

  if ( tolabs < t ) then
    lrline = 1
  else if ( -tolabs <= t ) then
    lrline = 0
  else if ( t < -tolabs ) then
    lrline = -1
  end if

  return
end
subroutine perm_inv ( n, p )

!*****************************************************************************80
!
!! PERM_INV inverts a permutation "in place".
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of objects being permuted.
!
!    Input/output, integer ( kind = 4 ) P(N), the permutation, in standard
!    index form.  On output, P describes the inverse permutation
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i0
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) is
  integer ( kind = 4 ) p(n)

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PERM_INV - Fatal error!'
    write ( *, '(a,i8)' ) '  Input value of N = ', n
    stop
  end if

  is = 1

  do i = 1, n

    i1 = p(i)

    do while ( i < i1 )
      i2 = p(i1)
      p(i1) = -i2
      i1 = i2
    end do

    is = -sign ( 1, p(i) )
    p(i) = sign ( p(i), is )

  end do

  do i = 1, n

    i1 = -p(i)

    if ( 0 <= i1 ) then

      i0 = i

      do

        i2 = p(i1)
        p(i1) = i0

        if ( i2 < 0 ) then
          exit
        end if

        i0 = i1
        i1 = i2

      end do

    end if

  end do

  return
end
subroutine points_nearest_point_bins_2d ( nset, pset, nbin, bin_min, bin_max, &
  bin_start, bin_last, bin_next, p, i_min, dist_min, compares )

!*****************************************************************************80
!
!! POINTS_NEAREST_POINT_BINS_2D finds the nearest point to a given point in 2D.
!
!  Discussion:
!
!    A set of NSET points with coordinates PSET is assumed to lie within a
!    rectangle.  The limits of this rectangle are given in BIN_MIN and BIN_MAX.
!    The rectangle is divided up into an NBIN by NBIN regular grid of cells.
!
!    The cells are ordered lexicographically, as suggested by the following
!    diagram for NBIN = 5:
!
!     *----*----*----*----*----*
!     |    |    |    |    |    |
!     | 21 | 22 | 23 | 24 | 25 |
!     |    |    |    |    |    |
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
!    1) for a point P, determine the cell it belongs to.
!       Note that this algorithm will NOT be valid if P lies outside
!       the bin limits.
!
!    2) Search for a cell that has at least one member of PSET in it.
!       We start at the cell containing P, but if there are no members
!       there, we spiral out from the cell, one layer at a time.
!
!    3) Within this cell, find the point nearest to P.  We now know that
!       we don't need to search any cell whose points will all be further
!       from P than this distance.
!
!    4) Now, search in all other cells that could have a point closer to
!       P than what we have found so far.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 July 2008
!
!  Author:
!
!    John Burkardt
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
!    Input, integer ( kind = 4 ) NSET, the number of points in the set.
!
!    Input, real ( kind = 8 ) PSET(2,NSET), the points in the set.
!
!    Input, integer ( kind = 4 ) NBIN, the number of cells.
!
!    Input, real ( kind = 8 ) BIN_MIN(2), BIN_MAX(2), the minimum and
!    maximum bin values.
!
!    Input, integer ( kind = 4 ) BIN_START(NBIN,NBIN), BIN_LAST(NBIN,NBIN),
!    indicates the index of the first and last element in the bin, or -1 if none.
!
!    Input, integer ( kind = 4 ) BIN_NEXT(NSET), the index of the next element
!    of the bin containing this element.
!
!    Input, real ( kind = 8 ) P(2), the point to be tested.
!
!    Output, integer ( kind = 4 ) I_MIN, the index of the nearest point in
!    PSET to P.
!
!    Output, real ( kind = 8 ) DIST_MIN, the distance between P and
!    PSET(*,I_MIN).
!
!    Output, integer ( kind = 4 ) COMPARES, the number of point-to-point
!    comparisons.
!
  implicit none

  integer ( kind = 4 ) nbin
  integer ( kind = 4 ), parameter :: ndim = 2
  integer ( kind = 4 ) nset

  integer ( kind = 4 ) bin(ndim)
  integer ( kind = 4 ) bin_last(nbin,nbin)
  real ( kind = 8 ) bin_max(ndim)
  real ( kind = 8 ) bin_min(ndim)
  integer ( kind = 4 ) bin_start(nbin,nbin)
  integer ( kind = 4 ) bin_next(nset)
  real ( kind = 8 ) c_max(ndim)
  real ( kind = 8 ) c_min(ndim)
  integer ( kind = 4 ) compares
  real ( kind = 8 ) dist_min
  real ( kind = 8 ) dist_min_sq
  real ( kind = 8 ) dist_sq
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_min
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) il
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jc
  integer ( kind = 4 ) jl
  integer ( kind = 4 ) layer
  real ( kind = 8 ) layer_width
  logical more_bins
  integer ( kind = 4 ) node
  real ( kind = 8 ) p(ndim)
  real ( kind = 8 ) pset(ndim,nset)
  real ( kind = 8 ) search_radius

  compares = 0
!
!  Special cases.
!
  if ( nset <= 0 ) then
    dist_min = huge ( dist_min )
    i_min = 0
    return
  end if

  if ( nset == 1 ) then
    dist_min = sqrt ( sum ( ( p(1:ndim) - pset(1:ndim,1) )**2 ) )
    compares = 1
    i_min = 1
    return
  end if
!
!  Initialize.
!
  dist_min_sq = huge ( dist_min_sq )
  i_min = 0
  search_radius = 0.0D+00

  layer_width = minval ( &
    ( bin_max(1:ndim) - bin_min(1:ndim) ) / real ( nbin, kind = 8 ) &
  )
!
!  Determine the bin coordinates of the point P.
!
  call r82_to_bin_even2 ( nbin, bin_min, bin_max, p, bin )
!
!  Determine the radius of the ball of space that will be completely
!  searched after this first bin, "layer 0", is done.
!
  call bin_to_r82_even2 ( nbin, bin, bin_min, bin_max, c_min, c_max )
!
!  Set the central bin of the layers.
!
  ic = bin(1)
  jc = bin(2)
!
!  Set
!  * the current layer,
!  * the starting bin of the current layer,
!  * the current bin
!
  layer = 0
  il = ic
  jl = jc
  i = il
  j = jl

  do
!
!  Search all legal bins in layer LAYER.
!
    do
!
!  Search BIN I, J.
!
      if ( 1 <= i .and. i <= nbin .and. 1 <= j .and. j <= nbin ) then

        node = bin_start(i,j)

        do while ( 0 < node )

          dist_sq = sum ( ( p(1:ndim) - pset(1:ndim,node) )**2 )
          compares = compares + 1

          if ( dist_sq < dist_min_sq ) then
            dist_min_sq = dist_sq
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
      more_bins = .true.

      do

        if ( i < ic + layer .and. j == jc - layer ) then
          i = i + 1
        else if ( i == ic + layer .and. j < jc + layer ) then
          j = j + 1
        else if ( ic - layer < i .and. j == jc + layer ) then
          i = i - 1
        else if ( i == ic - layer .and. jc - layer + 1 < j ) then
          j = j - 1
        else
          more_bins = .false.
          exit
        end if

        if ( 1 <= i .and. i <= nbin .and. &
             1 <= j .and. j <= nbin ) then
          exit
        end if

      end do

      if ( .not. more_bins ) then
        exit
      end if

    end do
!
!  We've completed this layer.
!  Update the radius of the searched area.
!
    if ( layer == 0 ) then
      search_radius = min ( &
        minval ( abs ( p(1:ndim) - c_min(1:ndim) ) ), &
        minval ( abs ( p(1:ndim) - c_max(1:ndim) ) ) )
    else
      search_radius = search_radius + layer_width
    end if
!
!  We can stop if:
!
!    * We've found at least one neighbor;
!    AND
!    * We've searched all legal boxes that contain the circle
!      with P at the center and the nearest N on the circumference.
!
   if ( i_min /= 0 ) then
     dist_min = sqrt ( dist_min_sq )
     if ( dist_min <= search_radius ) then
       exit
     end if
   end if
!
!  Prepare to search the next layer.
!
    layer = layer + 1

    il = ic - layer
    jl = jc - layer

    i = il
    j = jl

  end do

  return
end
subroutine points_nearest_point_bins_2d_2 ( nset, pset, nbin, bin_min, bin_max, &
  bin_start, bin_last, bin_next, ptest, i_min, dist_min, compares )

!*****************************************************************************80
!
!! POINTS_NEAREST_POINT_BINS_2D_2 finds the nearest point to a given point in 2D.
!
!  Discussion:
!
!    This code differs from POINTS_NEAREST_POINT_BINS_2D by calling
!    a subroutine to compute the next bin index.
!
!    A set of NSET points with coordinates PSET is assumed to lie within a
!    rectangle.  The limits of this rectangle are given in BIN_MIN and BIN_MAX.
!    The rectangle is divided up into an NBIN by NBIN regular grid of cells.
!
!    The cells are ordered lexicographically, as suggested by the following
!    diagram for NBIN = 5:
!
!     *----*----*----*----*----*
!     |    |    |    |    |    |
!     | 21 | 22 | 23 | 24 | 25 |
!     |    |    |    |    |    |
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
!    23 July 2008
!
!  Author:
!
!    John Burkardt
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
!    Input, integer ( kind = 4 ) NSET, the number of points in the set.
!
!    Input, real ( kind = 8 ) PSET(2,NSET), the points in the set.
!
!    Input, integer ( kind = 4 ) NBIN, the number of cells.
!
!    Input, real ( kind = 8 ) BIN_MIN(2), BIN_MAX(2), the minimum and
!    maximum bin values.
!
!    Input, integer ( kind = 4 ) BIN_START(NBIN,NBIN), BIN_LAST(NBIN,NBIN),
!    indicates the index of the first and last element in the bin, or -1 if none.
!
!    Input, integer ( kind = 4 ) BIN_NEXT(NSET), the index of the next element
!    of the bin containing this element.
!
!    Input, real ( kind = 8 ) PTEST(2), the coordinates of the test points.
!
!    Output, integer ( kind = 4 ) I_MIN, the index of the nearest point in
!    PSET to PTEST.
!
!    Output, real ( kind = 8 ) DIST_MIN, the distance between PTEST and
!    PSET(*,I_MIN).
!
!    Output, integer ( kind = 4 ) COMPARES, the number of point-to-point
!    comparisons.
!
  implicit none

  integer ( kind = 4 ) nbin
  integer ( kind = 4 ), parameter :: ndim = 2
  integer ( kind = 4 ) nset

  integer ( kind = 4 ) bin(ndim)
  integer ( kind = 4 ) bin_last(nbin,nbin)
  real ( kind = 8 ) bin_max(ndim)
  real ( kind = 8 ) bin_min(ndim)
  integer ( kind = 4 ) bin_start(nbin,nbin)
  integer ( kind = 4 ) bin_next(nset)
  real ( kind = 8 ) c_max(ndim)
  real ( kind = 8 ) c_min(ndim)
  integer ( kind = 4 ) compares
  real ( kind = 8 ) dist_min
  real ( kind = 8 ) dist_min_sq
  real ( kind = 8 ) dist_sq
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_min
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jc
  integer ( kind = 4 ) layer
  real ( kind = 8 ) layer_width
  logical more_bins
  integer ( kind = 4 ) node
  real ( kind = 8 ) pset(ndim,nset)
  real ( kind = 8 ) ptest(ndim)
  real ( kind = 8 ) search_radius

  compares = 0
!
!  Special cases.
!
  if ( nset <= 0 ) then
    dist_min = huge ( dist_min )
    i_min = 0
    return
  end if

  if ( nset == 1 ) then
    dist_min = sqrt ( sum ( ( ptest(1:ndim) - pset(1:ndim,1) )**2 ) )
    compares = 1
    i_min = 1
    return
  end if

  layer_width = minval ( &
    ( bin_max(1:ndim) - bin_min(1:ndim) ) / real ( nbin, kind = 8 ) &
  )
!
!  Search for each of the NTEST points.
!
  dist_min_sq = huge ( dist_min_sq )
  i_min = 0
  search_radius = 0.0D+00
!
!  Determine the bin coordinates of the point P.
!
  call r82_to_bin_even2 ( nbin, bin_min, bin_max, ptest, bin )
!
!  Determine the radius of the ball of space that will be completely
!  searched after this first bin, "layer 0", is done.
!
  call bin_to_r82_even2 ( nbin, bin, bin_min, bin_max, c_min, c_max )
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

      if ( 1 <= i .and. i <= nbin .and. 1 <= j .and. j <= nbin ) then

        node = bin_start(i,j)

        do while ( 0 < node )

          dist_sq = sum ( ( ptest(1:ndim) - pset(1:ndim,node) )**2 )
          compares = compares + 1

          if ( dist_sq < dist_min_sq ) then
            dist_min_sq = dist_sq
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

        if ( 1 <= i .and. i <= nbin .and. &
             1 <= j .and. j <= nbin ) then
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
      dist_min = sqrt ( dist_min_sq )
      if ( dist_min <= search_radius ) then
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
subroutine points_nearest_point_bins_2d_3 ( nset, pset, nbin, bin_min, &
  bin_max, bin_start, bin_last, bin_next, ptest, i_min, dist_min, compares )

!*****************************************************************************80
!
!! POINTS_NEAREST_POINT_BINS_2D_3 finds the nearest point to a given point in 2D.
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
!    23 July 2008
!
!  Author:
!
!    John Burkardt
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
!    Input, integer ( kind = 4 ) NSET, the number of points in the set.
!
!    Input, real ( kind = 8 ) PSET(2,NSET), the points in the set.
!
!    Input, integer ( kind = 4 ) NBIN(2), the number of cells in the horizontal
!    and vertical directions.
!
!    Input, real ( kind = 8 ) BIN_MIN(2), BIN_MAX(2), the minimum and
!    maximum bin values.
!
!    Input, integer ( kind = 4 ) BIN_START(NBIN(1),NBIN(2)),
!    BIN_LAST(NBIN(1),NBIN(2)), indicate the index of the first and last
!    elements in the bin, or -1 if none.
!
!    Input, integer ( kind = 4 ) BIN_NEXT(NSET), the index of the next element
!    of the bin containing this element.
!
!    Input, real ( kind = 8 ) PTEST(2), the coordinates of the test point.
!
!    Output, integer ( kind = 4 ) I_MIN, the index of the nearest point in
!    PSET to PTEST.
!
!    Output, real ( kind = 8 ) DIST_MIN, the distance between PTEST and
!    PSET(*,I_MIN).
!
!    Output, integer ( kind = 4 ) COMPARES, the number of point-to-point
!    comparisons.
!
  implicit none

  integer ( kind = 4 ), parameter :: ndim = 2

  integer ( kind = 4 ) nbin(ndim)
  integer ( kind = 4 ) nset

  integer ( kind = 4 ) bin(ndim)
  integer ( kind = 4 ) bin_last(nbin(1),nbin(2))
  real ( kind = 8 ) bin_max(ndim)
  real ( kind = 8 ) bin_min(ndim)
  integer ( kind = 4 ) bin_start(nbin(1),nbin(2))
  integer ( kind = 4 ) bin_next(nset)
  real ( kind = 8 ) c_max(ndim)
  real ( kind = 8 ) c_min(ndim)
  integer ( kind = 4 ) compares
  real ( kind = 8 ) dist_min
  real ( kind = 8 ) dist_min_sq
  real ( kind = 8 ) dist_sq
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_min
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jc
  integer ( kind = 4 ) layer
  real ( kind = 8 ) layer_width
  logical more_bins
  integer ( kind = 4 ) node
  real ( kind = 8 ) pset(ndim,nset)
  real ( kind = 8 ) ptest(ndim)
  real ( kind = 8 ) search_radius

  compares = 0
!
!  Special cases.
!
  if ( nset <= 0 ) then
    dist_min = huge ( dist_min )
    i_min = 0
    return
  end if

  if ( nset == 1 ) then
    dist_min = sqrt ( sum ( ( ptest(1:ndim) - pset(1:ndim,1) )**2 ) )
    i_min = 1
    compares = 1
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

  dist_min_sq = huge ( dist_min_sq )
  i_min = 0
  search_radius = 0.0D+00
!
!  Determine the bin coordinates of the point P.
!
  call r82_to_bin_even3 ( nbin, bin_min, bin_max, ptest, bin )
!
!  Determine the radius of the ball of space that will be completely
!  searched after this first bin, "layer 0", is done.
!
  call bin_to_r82_even3 ( nbin, bin, bin_min, bin_max, c_min, c_max )
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

          dist_sq = sum ( ( ptest(1:ndim) - pset(1:ndim,node) )**2 )
          compares = compares + 1

          if ( dist_sq < dist_min_sq ) then
            dist_min_sq = dist_sq
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
      dist_min = sqrt ( dist_min_sq )
      if ( dist_min <= search_radius ) then
        exit
      end if
    end if

    layer = layer + 1

  end do

  return
end
subroutine points_nearest_point_bins_3d_3 ( nset, pset, nbin, bin_min, &
  bin_max, bin_start, bin_last, bin_next, ptest, i_min )

!*****************************************************************************80
!
!! POINTS_NEAREST_POINT_BINS_3D_3 finds the nearest point to a given point in 3D.
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
!       we do not need to search any cell whose points will all be further
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
!    30 August 2008
!
!  Author:
!
!    John Burkardt
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
!    Input, integer ( kind = 4 ) NSET, the number of points in
!    the set.
!
!    Input, real ( kind = 8 ) PSET(3,NSET), the coordinates of the points
!    in the set.
!
!    Input, integer ( kind = 4 ) NBIN(3), the number of cells in the X, Y
!    and Z directions.
!
!    Input, real ( kind = 8 ) BIN_MIN(3), BIN_MAX(3), the minimum and
!    maximum bin values.
!
!    Input, integer ( kind = 4 ) BIN_START(NBIN(1),NBIN(2),NBIN(3)),
!    BIN_LAST(NBIN(1),NBIN(2),NBIN(3)), the index of the first and last
!    element in the bin, or -1 if none.
!
!    Input, integer ( kind = 4 ) BIN_NEXT(NSET), the index of the next element
!    of the bin containing this element.
!
!    Input, real ( kind = 8 ) PTEST(3), the coordinates of the test points.
!
!    Output, integer ( kind = 4 ) I_MIN, the index of the nearest point in
!    PSET to PTEST.
!
  implicit none

  integer ( kind = 4 ), parameter :: ndim = 3

  integer ( kind = 4 ) nbin(ndim)
  integer ( kind = 4 ) nset

  integer ( kind = 4 ) bin(ndim)
  integer ( kind = 4 ) bin_last(nbin(1),nbin(2),nbin(3))
  real    ( kind = 8 ) bin_max(ndim)
  real    ( kind = 8 ) bin_min(ndim)
  integer ( kind = 4 ) bin_start(nbin(1),nbin(2),nbin(3))
  integer ( kind = 4 ) bin_next(nset)
  real    ( kind = 8 ) c_max(ndim)
  real    ( kind = 8 ) c_min(ndim)
  real    ( kind = 8 ) d_min
  real    ( kind = 8 ) d_min_sq
  real    ( kind = 8 ) d_sq
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_min
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jc
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kc
  integer ( kind = 4 ) layer
  real    ( kind = 8 ) layer_width
  logical              more_bins
  integer ( kind = 4 ) node
  real    ( kind = 8 ) pset(ndim,nset)
  real    ( kind = 8 ) ptest(ndim)
  real    ( kind = 8 ) r8_huge
  real    ( kind = 8 ) search_radius
!
!  Special cases.
!
  if ( nset <= 0 ) then
    d_min = r8_huge ( )
    i_min = 0
    return
  end if

  if ( nset == 1 ) then
    d_min = 0.0
    do i = 1, ndim
      d_min = d_min + ( ptest(i) - pset(i,1) )**2
    end do
    d_min = sqrt ( d_min )
    i_min = 1
    return
  end if
!
!  The efficien!y of the code will suffer if the data in the vector
!
!    bin_max(1:ndim) - bin_min(1:ndim) / real ( nbin(1:ndim) )
!
!  varies significantly.
!
  layer_width = r8_huge ( )
  do i = 1, ndim
    layer_width = min ( layer_width, &
      ( bin_max(i) - bin_min(i) ) / real ( nbin(i), kind = 8 ) )
  end do

  d_min_sq = r8_huge ( )
  i_min = 0
  search_radius = 0.0D+00
!
!  Determine the bin coordinates of the point P.
!
  call r8ve!_to_bin_even3 ( ndim, nbin, bin_min, bin_max, ptest(1), bin )
!
!  Determine the radius of the ball of space that will be completely
!  searched after this first bin, "layer 0", is done.
!
  call bin_to_r8ve!_even3 ( ndim, nbin, bin, bin_min, bin_max, c_min, c_max )
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

    call index_box2_next_3d ( layer, layer, layer, ic, jc, kc, i, j, k, more_bins )
!
!  In layer LAYER, search each BIN I, J, K.
!
    do

      if ( 1 <= i .and. i <= nbin(1) .and. &
           1 <= j .and. j <= nbin(2) .and. &
           1 <= k .and. k <= nbin(3) ) then

        node = bin_start(i,j,k)

        do while ( 0 < node )

          d_sq = 0.0
          do i = 1, ndim
            d_sq = d_sq + ( ptest(i) - pset(i,node) )**2
          end do

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
!  Once you get past the last bin, exit because you are done the layer.
!
      do

        call index_box2_next_3d ( layer, layer, layer, ic, jc, &
          kc, i, j, k, more_bins )

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
!  We have completed layer LAYER.
!  Update the radius of the searched area.
!
    if ( layer == 0 ) then
      search_radius = r8_huge ( )
      do i = 1, ndim
        search_radius = min ( search_radius, abs ( ptest(i) - c_min(i) ) )
      end do
      do i = 1, ndim
        search_radius = min ( search_radius, abs ( ptest(i) - c_max(i) ) )
      end do
    else
      search_radius = search_radius + layer_width
    end if
!
!  We are done with PTEST if:
!
!    * We have found at least one neighbor;
!    AND
!    * We have searched all legal boxes that contain the circle
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
subroutine points_nearest_point_del_2d ( point_num, xc, xd, nabes_first, &
  nabes_num, nabes_dim, nabes, nnear, dnear )

!*****************************************************************************80
!
!! POINTS_NEAREST_POINT_DEL_2D: nearest neighbor in a Delaunay triangulation.
!
!  Discussion:
!
!    A set of points XC is given, along with its Delaunay triangulation.
!    Now a new point XD is given, and we need to know the closest point in XC.
!
!    This algorithm starts at a random point in XC, and then repeatedly moves
!    to a neighboring point that is closer to XD.  This is guaranteed to be
!    possible because the triangulation is Delaunay.  Otherwise, it
!    would be possible to reach a vertex which was not the closest,
!    but for which all neighbors were further away.
!
!    This algorithm is able to handle the case where the point XD lies
!    outside the convex hull.
!
!    The algorithm is very simple to code.  In the most likely
!    case, it should have an expected time complexity of O(sqrt(N)).
!
!    Overhead occurs in the development of the vertex adjacency data
!    structure.  The size of this array should be roughly 6*N on average.
!    Given the list of nodes that make up triangles, the vertex adjacency
!    data can be constructed by storing every pair of nodes (I,J) and (J,I),
!    and sorting the data into dictionary order.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 July 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, real ( kind = 8 ) XC(2,POINT_NUM), the set of points.
!
!    Input, real ( kind = 8 ) XD(2), a point whose nearest neighbor is
!    to be found.
!
!    Input, integer ( kind = 4 ) NABES_FIRST(POINT_NUM), the index in NABES of the first
!    neighbor in the list for each node.
!
!    Input, integer ( kind = 4 ) NABES_NUM(POINT_NUM), the number of neighbors
!    of each node.
!
!    Input, integer ( kind = 4 ) NABES_DIM, the dimension of NABES.
!
!    Input, integer ( kind = 4 ) NABES(NABES_DIM), a list of the neighbors of
!    all the nodes.  Neighbors of node 1 are listed first, and so on.
!
!    Output, integer ( kind = 4 ) NNEAR, the nearest node to XD.
!
!    Output, real ( kind = 8 ) DNEAR, the distance of the nearest node to XD.
!
  implicit none

  integer ( kind = 4 ) nabes_dim
  integer ( kind = 4 ) point_num

  real ( kind = 8 ) dist
  real ( kind = 8 ) dnear
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nabes(nabes_dim)
  integer ( kind = 4 ) nabes_first(point_num)
  integer ( kind = 4 ) nabes_num(point_num)
  integer ( kind = 4 ) nnear
  integer ( kind = 4 ) nnear_old
  real ( kind = 8 ) x
  real ( kind = 8 ) x1
  real ( kind = 8 ) xc(2,point_num)
  real ( kind = 8 ) xd(2)
  real ( kind = 8 ) y
  real ( kind = 8 ) y1

  x = xd(1)
  y = xd(2)
!
!  Select a random vertex.
!
  nnear = 1
  x1 = xc(1,nnear)
  y1 = xc(2,nnear)
  dnear = sqrt ( ( x1 - x )**2 + ( y1 - y )**2 )
!
!  From the current vertex, consider all neighboring vertices.
!  For each neighbor, compute the distance to the point.
!  If no neighbor is closer, then the current vertex is the closest.
!  Otherwise, set the current vertex to the neighbor that was closest,
!  and repeat.
!
  do

    nnear_old = nnear

    j = nabes_first(nnear_old)

    do i = 1, nabes_num(nnear_old)

      i1 = nabes(j)
      x1 = xc(1,i1)
      y1 = xc(2,i1)
      dist = sqrt ( ( x1 - x )**2 + ( y1 - y )**2 )

      if ( dist < dnear ) then
        dnear = dist
        nnear = i1
      end if

      j = j + 1

    end do
!
!  If no neighbor was closer, we're done.
!
    if ( nnear == nnear_old ) then
      exit
    end if

  end do

  return
end
subroutine points_nearest_point_naive_2d ( nset, pset, ptest, i_min, dist_min )

!*****************************************************************************80
!
!! POINTS_NEAREST_POINT_NAIVE_2D finds the nearest point to a given point in 2D.
!
!  Discussion:
!
!    A naive algorithm is used.  The distance to every point is calculated,
!    in order to determine the smallest.
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
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NSET, the number of points in the set.
!
!    Input, real ( kind = 8 ) PSET(2,NSET), the points in the set.
!
!    Input, real ( kind = 8 ) PTEST(2), the point whose nearest neighbor
!    is sought.
!
!    Output, integer ( kind = 4 ) I_MIN, the index of the nearest point in
!    PSET to P.
!
!    Output, real ( kind = 8 ) DIST_MIN, the distance between P and PSET(*,I_MIN).
!
  implicit none

  integer ( kind = 4 ) nset
  integer ( kind = 4 ), parameter :: ndim = 2

  real ( kind = 8 ) d
  real ( kind = 8 ) dist_min
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_min
  real ( kind = 8 ) pset(ndim,nset)
  real ( kind = 8 ) ptest(ndim)

  dist_min = huge ( dist_min )
  i_min = -1

  do i = 1, nset
    d = sum ( ( ptest(1:ndim) - pset(1:ndim,i) )**2 )
    if ( d < dist_min ) then
      dist_min = d
      i_min = i
    end if
  end do

  dist_min = sqrt ( dist_min )

  return
end
subroutine points_nearest_point_naive_3d ( nset, pset, ptest, i_min, dist_min )

!*****************************************************************************80
!
!! POINTS_NEAREST_POINT_NAIVE_3D finds the nearest point to a given point in 3D.
!
!  Discussion:
!
!    A naive algorithm is used.  The distance to every point is calculated,
!    in order to determine the smallest.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NSET, the number of points in the set.
!
!    Input, real ( kind = 8 ) PSET(3,NSET), the points in the set.
!
!    Input, real ( kind = 8 ) PTEST(3), the point whose nearest neighbor
!    is sought.
!
!    Output, integer ( kind = 4 ) I_MIN, the index of the nearest point in
!    PSET to P.
!
!    Output, real ( kind = 8 ) DIST_MIN, the distance between P and PSET(*,I_MIN).
!
  implicit none

  integer ( kind = 4 ) nset
  integer ( kind = 4 ), parameter :: ndim = 3

  real ( kind = 8 ) d
  real ( kind = 8 ) dist_min
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_min
  real ( kind = 8 ) pset(ndim,nset)
  real ( kind = 8 ) ptest(ndim)

  dist_min = huge ( dist_min )
  i_min = -1

  do i = 1, nset
    d = sum ( ( ptest(1:ndim) - pset(1:ndim,i) )**2 )
    if ( d < dist_min ) then
      dist_min = d
      i_min = i
    end if
  end do

  dist_min = sqrt ( dist_min )

  return
end
subroutine points_nearest_point_naive_nd ( ndim, n, pset, p, i_min, dist_min )

!*****************************************************************************80
!
!! POINTS_NEAREST_POINT_NAIVE_ND finds the nearest point to a given point in ND.
!
!  Discussion:
!
!    A naive algorithm is used.  No attempt is made to optimize the
!    calculation, so there will be N distance calculations done.
!
!    For a large dataset, it would be better to group the points into
!    clusters, so that far fewer distance calculations need to be done.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NDIM, the dimension of the points.
!
!    Input, integer ( kind = 4 ) N, the number of points in the set.
!
!    Input, real ( kind = 8 ) PSET(NDIM,N), the coordinates of the points
!    in the set.
!
!    Input, real ( kind = 8 ) P(NDIM), the point to be tested.
!
!    Output, integer ( kind = 4 ) I_MIN, the index of the nearest point in
!    PSET to P.
!
!    Output, real ( kind = 8 ) DIST_MIN, the distance between P and PSET(*,I_MIN).
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) ndim

  real ( kind = 8 ) d
  real ( kind = 8 ) dist_min
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_min
  real ( kind = 8 ) p(ndim)
  real ( kind = 8 ) pset(ndim,n)

  dist_min = huge ( dist_min )
  i_min = -1

  do i = 1, n

    d = sum ( ( p(1:ndim) - pset(1:ndim,i) )**2 )

    if ( d < dist_min ) then
      dist_min = d
      i_min = i
    end if

  end do
!
!  We save a little work by waiting til the end to take the square root.
!
  dist_min = sqrt ( dist_min )

  return
end
subroutine points_nearest_points_bins_2d ( nset, pset, nbin, bin_min, bin_max, &
  bin_start, bin_last, bin_next, ntest, ptest, i_min, dist_min, compares )

!*****************************************************************************80
!
!! POINTS_NEAREST_POINTS_BINS_2D finds the nearest point to given points in 2D.
!
!  Discussion:
!
!    A set of NSET points with coordinates PSET is assumed to lie within a
!    rectangle.  The limits of this rectangle are given in BIN_MIN and BIN_MAX.
!    The rectangle is divided up into an NBIN by NBIN regular grid of cells.
!
!    The cells are ordered lexicographically, as suggested by the following
!    diagram for NBIN = 5:
!
!     *----*----*----*----*----*
!     |    |    |    |    |    |
!     | 21 | 22 | 23 | 24 | 25 |
!     |    |    |    |    |    |
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
!    23 July 2008
!
!  Author:
!
!    John Burkardt
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
!    Input, integer ( kind = 4 ) NSET, the number of points in the set.
!
!    Input, real ( kind = 8 ) PSET(2,NSET), the points in the set.
!
!    Input, integer ( kind = 4 ) NBIN, the number of cells.
!
!    Input, real ( kind = 8 ) BIN_MIN(2), BIN_MAX(2), the minimum and
!    maximum bin values.
!
!    Input, integer ( kind = 4 ) BIN_START(NBIN,NBIN), BIN_LAST(NBIN,NBIN),
!    indicates the index of the first and last element in the bin, or -1 
!    if none.
!
!    Input, integer ( kind = 4 ) BIN_NEXT(NSET), the index of the next element
!    of the bin containing this element.
!
!    Input, integer ( kind = 4 ) NTEST, the number of test points.
!
!    Input, real ( kind = 8 ) PTEST(2,NTEST), the test points.
!
!    Output, integer ( kind = 4 ) I_MIN(NTEST), the index of the nearest point 
!    in PSET to PTEST(ITEST).
!
!    Output, real ( kind = 8 ) DIST_MIN(NTEST), the distance between
!    PTEST(*,ITEST) and PSET(*,I_MIN).
!
!    Output, integer ( kind = 4 ) COMPARES(NTEST), the number of 
!    point-to-point comparisons.
!
  implicit none

  integer ( kind = 4 ) nbin
  integer ( kind = 4 ), parameter :: ndim = 2
  integer ( kind = 4 ) nset
  integer ( kind = 4 ) ntest

  integer ( kind = 4 ) bin(ndim)
  integer ( kind = 4 ) bin_last(nbin,nbin)
  real ( kind = 8 ) bin_max(ndim)
  real ( kind = 8 ) bin_min(ndim)
  integer ( kind = 4 ) bin_start(nbin,nbin)
  integer ( kind = 4 ) bin_next(nset)
  real ( kind = 8 ) c_max(ndim)
  real ( kind = 8 ) c_min(ndim)
  integer ( kind = 4 ) compares(ntest)
  real ( kind = 8 ) dist_min(ntest)
  real ( kind = 8 ) dist_min_sq
  real ( kind = 8 ) dist_sq
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_min(ntest)
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) il
  integer ( kind = 4 ) itest
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jc
  integer ( kind = 4 ) jl
  integer ( kind = 4 ) layer
  real ( kind = 8 ) layer_width
  logical more_bins
  integer ( kind = 4 ) node
  real ( kind = 8 ) pset(ndim,nset)
  real ( kind = 8 ) ptest(ndim,ntest)
  real ( kind = 8 ) search_radius

  compares(1:ntest) = 0
!
!  Special cases.
!
  if ( nset <= 0 ) then
    dist_min(1:ntest) = huge ( dist_min )
    i_min(1:ntest) = 0
    return
  end if

  if ( nset == 1 ) then
    do itest = 1, ntest
      dist_min(itest) = sqrt ( &
        sum ( ( ptest(1:ndim,itest) - pset(1:ndim,1) )**2 ) )
    end do
    compares(1:ntest) = 1
    i_min(1:ntest) = 1
    return
  end if

  layer_width = minval ( &
    ( bin_max(1:ndim) - bin_min(1:ndim) ) / real ( nbin, kind = 8 ) &
  )
!
!  Search for each of the NTEST points.
!
  do itest = 1, ntest

    dist_min_sq = huge ( dist_min_sq )
    i_min(itest) = 0
    search_radius = 0.0D+00
!
!  Determine the bin coordinates of the point P.
!
    call r82_to_bin_even2 ( nbin, bin_min, bin_max, ptest(1,itest), bin )
!
!  Determine the radius of the ball of space that will be completely
!  searched after this first bin, "layer 0", is done.
!
    call bin_to_r82_even2 ( nbin, bin, bin_min, bin_max, c_min, c_max )
!
!  Set the central bin of the layers.
!
    ic = bin(1)
    jc = bin(2)
!
!  Set
!  * the current layer,
!  * the starting bin of the current layer,
!  * the current bin
!
    layer = 0
    il = ic
    jl = jc
    i = il
    j = jl

    do
!
!  Search all legal bins in layer LAYER.
!
      do
!
!  Search BIN I, J.
!
        if ( 1 <= i .and. i <= nbin .and. 1 <= j .and. j <= nbin ) then

          node = bin_start(i,j)

          do while ( 0 < node )

            dist_sq = sum ( ( ptest(1:ndim,itest) - pset(1:ndim,node) )**2 )
            compares(itest) = compares(itest) + 1

            if ( dist_sq < dist_min_sq ) then
              dist_min_sq = dist_sq
              i_min(itest) = node
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
        more_bins = .true.

        do

          if ( i < ic + layer .and. j == jc - layer ) then
            i = i + 1
          else if ( i == ic + layer .and. j < jc + layer ) then
            j = j + 1
          else if ( ic - layer < i .and. j == jc + layer ) then
            i = i - 1
          else if ( i == ic - layer .and. jc - layer + 1 < j ) then
            j = j - 1
          else
            more_bins = .false.
            exit
          end if

          if ( 1 <= i .and. i <= nbin .and. &
               1 <= j .and. j <= nbin ) then
            exit
          end if

        end do

        if ( .not. more_bins ) then
          exit
        end if

      end do
!
!  We've completed this layer.
!  Update the radius of the searched area.
!
      if ( layer == 0 ) then
        search_radius = min ( &
          minval ( abs ( ptest(1:ndim,itest) - c_min(1:ndim) ) ), &
          minval ( abs ( ptest(1:ndim,itest) - c_max(1:ndim) ) ) )
      else
        search_radius = search_radius + layer_width
      end if
!
!  We are done with PTEST(*,ITEST) if:
!
!    * We've found at least one neighbor;
!    AND
!    * We've searched all legal boxes that contain the circle
!      with PTEST at the center and the nearest neighbor on the circumference.
!
     if ( i_min(itest) /= 0 ) then
       dist_min(itest) = sqrt ( dist_min_sq )
       if ( dist_min(itest) <= search_radius ) then
         exit
       end if
     end if
!
!  Prepare to search the next layer.
!
      layer = layer + 1

      il = ic - layer
      jl = jc - layer

      i = il
      j = jl

    end do

  end do

  return
end
subroutine points_nearest_points_bins_2d_2 ( nset, pset, nbin, bin_min, &
  bin_max, bin_start, bin_last, bin_next, ntest, ptest, i_min, dist_min, &
  compares )

!*****************************************************************************80
!
!! POINTS_NEAREST_POINTS_BINS_2D_2 finds the nearest point to given points in 2D.
!
!  Discussion:
!
!    This code differs from POINTS_NEAREST_POINTS_BINS_2D by calling
!    a subroutine to compute the next bin index.
!
!    A set of NSET points with coordinates PSET is assumed to lie within a
!    rectangle.  The limits of this rectangle are given in BIN_MIN and BIN_MAX.
!    The rectangle is divided up into an NBIN by NBIN regular grid of cells.
!
!    The cells are ordered lexicographically, as suggested by the following
!    diagram for NBIN = 5:
!
!     *----*----*----*----*----*
!     |    |    |    |    |    |
!     | 21 | 22 | 23 | 24 | 25 |
!     |    |    |    |    |    |
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
!    23 July 2008
!
!  Author:
!
!    John Burkardt
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
!    Input, integer ( kind = 4 ) NSET, the number of points in the set.
!
!    Input, real ( kind = 8 ) PSET(2,NSET), the points in the set.
!
!    Input, integer ( kind = 4 ) NBIN, the number of cells.
!
!    Input, real ( kind = 8 ) BIN_MIN(2), BIN_MAX(2), the minimum and
!    maximum bin values.
!
!    Input, integer ( kind = 4 ) BIN_START(NBIN,NBIN), BIN_LAST(NBIN,NBIN),
!    indicates the index of the first and last element in the bin, or -1 if none.
!
!    Input, integer ( kind = 4 ) BIN_NEXT(NSET), the index of the next element 
!    of the bin containing this element.
!
!    Input, integer ( kind = 4 ) NTEST, the number of test points.
!
!    Input, real ( kind = 8 ) PTEST(2,NTEST), the test points.
!
!    Output, integer ( kind = 4 ) I_MIN(NTEST), the index of the nearest point 
!    in PSET to PTEST(ITEST).
!
!    Output, real ( kind = 8 ) DIST_MIN(NTEST), the distance between
!    PTEST(*,ITEST) and PSET(*,I_MIN).
!
!    Output, integer ( kind = 4 ) COMPARES(NTEST), the number of 
!    point-to-point comparisons.
!
  implicit none

  integer ( kind = 4 ) nbin
  integer ( kind = 4 ), parameter :: ndim = 2
  integer ( kind = 4 ) nset
  integer ( kind = 4 ) ntest

  integer ( kind = 4 ) bin(ndim)
  integer ( kind = 4 ) bin_last(nbin,nbin)
  real ( kind = 8 ) bin_max(ndim)
  real ( kind = 8 ) bin_min(ndim)
  integer ( kind = 4 ) bin_start(nbin,nbin)
  integer ( kind = 4 ) bin_next(nset)
  real ( kind = 8 ) c_max(ndim)
  real ( kind = 8 ) c_min(ndim)
  integer ( kind = 4 ) compares(ntest)
  real ( kind = 8 ) dist_min(ntest)
  real ( kind = 8 ) dist_min_sq
  real ( kind = 8 ) dist_sq
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_min(ntest)
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) itest
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jc
  integer ( kind = 4 ) layer
  real ( kind = 8 ) layer_width
  logical more_bins
  integer ( kind = 4 ) node
  real ( kind = 8 ) pset(ndim,nset)
  real ( kind = 8 ) ptest(ndim,ntest)
  real ( kind = 8 ) search_radius

  compares(1:ntest) = 0
!
!  Special cases.
!
  if ( nset <= 0 ) then
    dist_min(1:ntest) = huge ( dist_min )
    i_min(1:ntest) = 0
    return
  end if

  if ( nset == 1 ) then
    do itest = 1, ntest
      dist_min(itest) = sqrt ( &
        sum ( ( ptest(1:ndim,itest) - pset(1:ndim,1) )**2 ) )
    end do
    compares(1:ntest) = 1
    i_min(1:ntest) = 1
    return
  end if

  layer_width = minval ( &
    ( bin_max(1:ndim) - bin_min(1:ndim) ) / real ( nbin, kind = 8 ) &
  )
!
!  Search for each of the NTEST points.
!
  do itest = 1, ntest

    dist_min_sq = huge ( dist_min_sq )
    i_min(itest) = 0
    search_radius = 0.0D+00
!
!  Determine the bin coordinates of the point P.
!
    call r82_to_bin_even2 ( nbin, bin_min, bin_max, ptest(1,itest), bin )
!
!  Determine the radius of the ball of space that will be completely
!  searched after this first bin, "layer 0", is done.
!
    call bin_to_r82_even2 ( nbin, bin, bin_min, bin_max, c_min, c_max )
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

        if ( 1 <= i .and. i <= nbin .and. 1 <= j .and. j <= nbin ) then

          node = bin_start(i,j)

          do while ( 0 < node )

            dist_sq = sum ( ( ptest(1:ndim,itest) - pset(1:ndim,node) )**2 )
            compares(itest) = compares(itest) + 1

            if ( dist_sq < dist_min_sq ) then
              dist_min_sq = dist_sq
              i_min(itest) = node
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

          if ( 1 <= i .and. i <= nbin .and. &
               1 <= j .and. j <= nbin ) then
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
          minval ( abs ( ptest(1:ndim,itest) - c_min(1:ndim) ) ), &
          minval ( abs ( ptest(1:ndim,itest) - c_max(1:ndim) ) ) )
      else
        search_radius = search_radius + layer_width
      end if
!
!  We are done with PTEST(*,ITEST) if:
!
!    * We've found at least one neighbor;
!    AND
!    * We've searched all legal boxes that contain the circle
!      with PTEST at the center and the nearest neighbor on the circumference.
!
      if ( i_min(itest) /= 0 ) then
        dist_min(itest) = sqrt ( dist_min_sq )
        if ( dist_min(itest) <= search_radius ) then
          exit
        end if
      end if

      layer = layer + 1

    end do
!
!  We are now done with all the layers.
!
  end do
!
!  We are now done with all the test points.
!
  return
end
subroutine points_nearest_points_bins_2d_3 ( nset, pset, nbin, bin_min, &
  bin_max, bin_start, bin_last, bin_next, ntest, ptest, i_min, dist_min, &
  compares )

!*****************************************************************************80
!
!! POINTS_NEAREST_POINTS_BINS_2D_3 finds the nearest point to given points in 2D.
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
!    John Burkardt
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
!    Input, integer ( kind = 4 ) NSET, the number of points in the set.
!
!    Input, real ( kind = 8 ) PSET(2,NSET), the points in the set.
!
!    Input, integer ( kind = 4 ) NBIN(2), the number of cells in the horizontal 
!    and vertical directions.
!
!    Input, real ( kind = 8 ) BIN_MIN(2), BIN_MAX(2), the minimum and
!    maximum bin values.
!
!    Input, integer ( kind = 4 ) BIN_START(NBIN(1),NBIN(2)), 
!    BIN_LAST(NBIN(1),NBIN(2)), indicates the index of the first and last 
!    element in the bin, or -1 if none.
!
!    Input, integer ( kind = 4 ) BIN_NEXT(NSET), the index of the next element 
!    of the bin containing this element.
!
!    Input, integer ( kind = 4 ) NTEST, the number of test points.
!
!    Input, real ( kind = 8 ) PTEST(2,NTEST), the coordinates of the
!    test points.
!
!    Output, integer ( kind = 4 ) I_MIN(NTEST), the index of the nearest point 
!    in PSET to PTEST(ITEST).
!
!    Output, real ( kind = 8 ) DIST_MIN(NTEST), the distance between
!    PTEST(*,ITEST) and PSET(*,I_MIN).
!
!    Output, integer ( kind = 4 ) COMPARES(NTEST), the number of point-to-point
!    comparisons.
!
  implicit none

  integer ( kind = 4 ), parameter :: ndim = 2

  integer ( kind = 4 ) nbin(ndim)
  integer ( kind = 4 ) nset
  integer ( kind = 4 ) ntest

  integer ( kind = 4 ) bin(ndim)
  integer ( kind = 4 ) bin_last(nbin(1),nbin(2))
  real ( kind = 8 ) bin_max(ndim)
  real ( kind = 8 ) bin_min(ndim)
  integer ( kind = 4 ) bin_start(nbin(1),nbin(2))
  integer ( kind = 4 ) bin_next(nset)
  real ( kind = 8 ) c_max(ndim)
  real ( kind = 8 ) c_min(ndim)
  integer ( kind = 4 ) compares(ntest)
  real ( kind = 8 ) dist_min(ntest)
  real ( kind = 8 ) dist_min_sq
  real ( kind = 8 ) dist_sq
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_min(ntest)
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) itest
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jc
  integer ( kind = 4 ) layer
  real ( kind = 8 ) layer_width
  logical more_bins
  integer ( kind = 4 ) node
  real ( kind = 8 ) pset(ndim,nset)
  real ( kind = 8 ) ptest(ndim,ntest)
  real ( kind = 8 ) search_radius

  compares(1:ntest) = 0
!
!  Special cases.
!
  if ( nset <= 0 ) then
    dist_min(1:ntest) = huge ( dist_min(1) )
    i_min(1:ntest) = 0
    return
  end if

  if ( nset == 1 ) then
    do itest = 1, ntest
      dist_min(itest) = &
        sqrt ( sum ( ( ptest(1:ndim,itest) - pset(1:ndim,1) )**2 ) )
    end do
    compares(1:ntest) = 1
    i_min(1:ntest) = 1
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
!
!  Search for each of the NTEST points.
!
  do itest = 1, ntest

    dist_min_sq = huge ( dist_min_sq )
    i_min(itest) = 0
    search_radius = 0.0D+00
!
!  Determine the bin coordinates of the point P.
!
    call r82_to_bin_even3 ( nbin, bin_min, bin_max, ptest(1,itest), bin )
!
!  Determine the radius of the ball of space that will be completely
!  searched after this first bin, "layer 0", is done.
!
    call bin_to_r82_even3 ( nbin, bin, bin_min, bin_max, c_min, c_max )
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

            dist_sq = sum ( ( ptest(1:ndim,itest) - pset(1:ndim,node) )**2 )
            compares(itest) = compares(itest) + 1

            if ( dist_sq < dist_min_sq ) then
              dist_min_sq = dist_sq
              i_min(itest) = node
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
          minval ( abs ( ptest(1:ndim,itest) - c_min(1:ndim) ) ), &
          minval ( abs ( ptest(1:ndim,itest) - c_max(1:ndim) ) ) )
      else
        search_radius = search_radius + layer_width
      end if
!
!  We are done with PTEST(*,ITEST) if:
!
!    * We've found at least one neighbor;
!    AND
!    * We've searched all legal boxes that contain the circle
!      with PTEST at the center and the nearest neighbor on the circumference.
!
      if ( i_min(itest) /= 0 ) then
        dist_min(itest) = sqrt ( dist_min_sq )
        if ( dist_min(itest) <= search_radius ) then
          exit
        end if
      end if

      layer = layer + 1

    end do
!
!  We are now done with all the layers.
!
  end do
!
!  We are now done with all the test points.
!
  return
end
subroutine points_nearest_points_bins_2d_4 ( nset, pset, nbin, bin_min, &
  bin_max, bin_start, bin_last, bin_next, ntest, ptest, i_min, dist_min, &
  compares )

!*****************************************************************************80
!
!! POINTS_NEAREST_POINTS_BINS_2D_4 finds the nearest point to given points in 2D.
!
!  Discussion:
!
!    This code differs from POINTS_NEAREST_POINTS_BINS_2D by allowing the
!    user to specify the number of bins in each dimension.  The main reason
!    for doing this is to efficiently handle problems where the extent
!    of the region varies widely from one dimension to another.
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
!    18 March 2001
!
!  Author:
!
!    John Burkardt
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
!    Input, integer ( kind = 4 ) NSET, the number of points in the set.
!
!    Input, real ( kind = 8 ) PSET(2,NSET), the points in the set.
!
!    Input, integer ( kind = 4 ) NBIN(2), the number of cells in the horizontal
!    and vertical directions.
!
!    Input, real ( kind = 8 ) BIN_MIN(2), BIN_MAX(2), the minimum and
!    maximum bin values.
!
!    Input, integer ( kind = 4 ) BIN_START(NBIN(1),NBIN(2)), 
!    BIN_LAST(NBIN(1),NBIN(2)), indicates the index of the first and last 
!    element in the bin, or -1 if none.
!
!    Input, integer ( kind = 4 ) BIN_NEXT(NSET), the index of the next element 
!    of the bin containing this element.
!
!    Input, integer ( kind = 4 ) NTEST, the number of test points.
!
!    Input, real ( kind = 8 ) PTEST(2,NTEST), the test points.
!
!    Output, integer ( kind = 4 ) I_MIN(NTEST), the index of the nearest point 
!    in PSET to PTEST(ITEST).
!
!    Output, real ( kind = 8 ) DIST_MIN(NTEST), the distance between
!    PTEST(*,ITEST) and PSET(*,I_MIN).
!
!    Output, integer ( kind = 4 ) COMPARES(NTEST), the number of point-to-point 
!    comparisons.
!
  implicit none

  integer ( kind = 4 ), parameter :: ndim = 2

  integer ( kind = 4 ) nbin(ndim)
  integer ( kind = 4 ) nset
  integer ( kind = 4 ) ntest

  integer ( kind = 4 ) bin(ndim)
  integer ( kind = 4 ) bin_last(nbin(1),nbin(2))
  real ( kind = 8 ) bin_max(ndim)
  real ( kind = 8 ) bin_min(ndim)
  integer ( kind = 4 ) bin_start(nbin(1),nbin(2))
  integer ( kind = 4 ) bin_next(nset)
  real ( kind = 8 ) c_max(ndim)
  real ( kind = 8 ) c_min(ndim)
  real ( kind = 8 ) cell_width_i
  integer ( kind = 4 ) compares(ntest)
  real ( kind = 8 ) dist_min(ntest)
  real ( kind = 8 ) dist_min_sq
  real ( kind = 8 ) disc
  real ( kind = 8 ) first_dist
  logical found_a_neighbor
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_min(ntest)
  integer ( kind = 4 ) itest
  integer ( kind = 4 ) k
  real ( kind = 8 ) pset(ndim,nset)
  real ( kind = 8 ) ptest(ndim,ntest)
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) search1(ndim)
  integer ( kind = 4 ) search2(ndim)
  logical searched_everywhere
  integer ( kind = 4 ) searched_hi(ndim)
  integer ( kind = 4 ) searched_hi_new(ndim)
  integer ( kind = 4 ) searched_lo(ndim)
  integer ( kind = 4 ) searched_lo_new(ndim)
  real ( kind = 8 ) w
  real ( kind = 8 ) wall
  real ( kind = 8 ) wall_dist
  real ( kind = 8 ) width
  real ( kind = 8 ) width_i

  compares(1:ntest) = 0
!
!  Special cases.
!
  if ( nset <= 0 ) then
    dist_min(1:ntest) = huge ( dist_min(1) )
    i_min(1:ntest) = 0
    return
  end if

  if ( nset == 1 ) then
    do itest = 1, ntest
      dist_min(itest) = sqrt ( &
        sum ( ( ptest(1:ndim,itest) - pset(1:ndim,1) )**2 ) )
    end do
    compares(1:ntest) = 1
    i_min(1:ntest) = 1
    return
  end if
!
!  Search for each of the NTEST points.
!
  do itest = 1, ntest
!
!  PART ONE: Initialize data.
!
!  Determine the bin coordinates of the test point.
!  Determine the limits of the bin containing P.
!  Search the bin.
!  Set the indices of the searched region to the index of this bin.
!
    dist_min_sq = huge ( dist_min_sq )
    i_min(itest) = 0
    found_a_neighbor = .false.
    searched_everywhere = .false.

    call r82_to_bin_even3 ( nbin, bin_min, bin_max, ptest(1,itest), bin )

    call bin_to_r82_even3 ( nbin, bin, bin_min, bin_max, c_min, c_max )

    call bin_search_one_2d ( bin, nset, pset, nbin, bin_start, bin_next, &
      ptest(1,itest), found_a_neighbor, i_min(itest), dist_min_sq, &
      compares(itest) )

    if ( found_a_neighbor) then
      dist_min(itest) = sqrt ( dist_min_sq )
    end if

    searched_lo(1:ndim) = bin(1:ndim)
    searched_hi(1:ndim) = bin(1:ndim)

    if ( all ( nbin(1:ndim) == 1 ) ) then
      if ( found_a_neighbor ) then
        dist_min(itest) = sqrt ( dist_min_sq )
      else
        dist_min(itest) = huge ( dist_min(1) )
      end if
      searched_everywhere = .true.
      cycle
    end if
!
!  PART TWO. Look for a neighbor.
!
!  Expand the search area in each dimension.
!  Consider push the search limits out in every direction by 1 index.
!  Determine the maximum width of the search region achieved in this way.
!  Where possible, move all smaller indices out to the maximum width.
!
!  Organize the search of the annexed region by dimension.
!  Decrease and increase the first dimension to its new limits.
!  Repeat for each dimension.
!
!  Jump to next section ALMOST as soon as you have found a neighbor;
!  just finish up the search in that coordinate direction, so that the
!  search can be picked up cleanly in the final section.
!
    do while ( .not. found_a_neighbor )

      do i = 1, ndim
        searched_hi_new(i) = min ( searched_hi(i) + 1, nbin(i) )
        searched_lo_new(i) = max ( searched_lo(i) - 1, 1 )
      end do

      width = 0.0D+00
      do i = 1, ndim
        width_i = ( searched_hi_new(i) + 1 - searched_lo_new(i) ) * &
            ( bin_max(i) - bin_min(i) ) / nbin(i)
        width = max ( width, width_i )
      end do

      do i = 1, ndim

        cell_width_i = ( bin_max(i) - bin_min(i) ) / nbin(i)

        width_i = ( searched_hi_new(i) + 1 - searched_lo_new(i) ) * cell_width_i

        disc = width - width_i

        if ( 2.0D+00 * cell_width_i <= disc ) then

          k = int ( disc / ( 2.0D+00 * cell_width_i ) )

          searched_hi_new(i) = searched_hi_new(i) + k
          searched_lo_new(i) = searched_lo_new(i) - k

          searched_hi_new(i) = min ( searched_hi_new(i), nbin(i) )
          searched_lo_new(i) = max ( searched_lo_new(i), 1 )

        end if

      end do

      searched_everywhere = .true.
      do i = 1, ndim
        if ( searched_hi(i) < searched_hi_new(i) ) then
          searched_everywhere = .false.
          exit
        end if
        if ( searched_lo_new(i) < searched_lo(i) ) then
          searched_everywhere = .false.
          exit
        end if
      end do

      if ( searched_everywhere ) then
        exit
      end if

      do i = 1, ndim

        if ( searched_lo_new(i) < searched_lo(i) ) then

          search1(1:ndim) = searched_lo(1:ndim)
          search1(i) = searched_lo_new(i)

          search2(1:i-1) = searched_hi(1:i-1)
          search2(i) = searched_lo(i) - 1
          search2(i+1:ndim) = searched_lo(i+1:ndim)

          rank = 0

          do

            call tuple_next2 ( ndim, search1, search2, bin, rank )

            if ( rank == 0 ) then
              exit
            end if

            call bin_search_one_2d ( bin, nset, pset, nbin, bin_start, &
              bin_next, ptest(1,itest), found_a_neighbor, i_min(itest), &
              dist_min_sq, compares(itest) )

          end do

          searched_lo(i) = searched_lo_new(i)

          if ( found_a_neighbor ) then
            exit
          end if

        end if

        if ( searched_hi(i) < searched_hi_new(i) ) then

          search1(1:i-1) = searched_lo(1:ndim)
          search1(i) = searched_hi(i) + 1
          search1(i+1:ndim) = searched_hi(1:ndim)

          search2(1:ndim) = searched_hi(1:ndim)
          search2(i) = searched_hi_new(i)

          rank = 0

          do

            call tuple_next2 ( ndim, search1, search2, bin, rank )

            if ( rank == 0 ) then
              exit
            end if

            call bin_search_one_2d ( bin, nset, pset, nbin, bin_start, &
              bin_next, ptest(1,itest), found_a_neighbor, i_min(itest), &
              dist_min_sq, compares(itest) )

          end do

          searched_hi(i) = searched_hi_new(i)

          if ( found_a_neighbor ) then
            exit
          end if

        end if

      end do

      dist_min(itest) = sqrt ( dist_min_sq )

    end do
!
!  PART THREE: Final search
!
!  You have found a neighbor in PSET to PTEST.
!  If the neighbor is closer than the nearest wall which might have
!    something on the other side, you're done.
!  Otherwise, expand the region enough in each direction so that, once
!    it is searched, we are sure to be done.
!
    wall_dist = huge ( wall_dist )

    do i = 1, ndim

      if ( 1 < searched_lo(i) ) then

        wall = &
          ( real ( nbin(i) - searched_lo(i) + 1, kind = 8 ) * bin_min(i) &
          + real (           searched_lo(i) - 1, kind = 8 ) * bin_max(i) ) &
          / real ( nbin(i),                      kind = 8 )

        wall_dist = min ( wall_dist, ptest(i,itest) - wall )

      end if

      if ( searched_hi(i) < nbin(i) ) then

        wall = &
          ( real ( nbin(i) - searched_hi(i), kind = 8 ) * bin_min(i) &
          + real (           searched_hi(i), kind = 8 ) * bin_max(i) ) &
          / real ( nbin(i),                  kind = 8 )

        wall_dist = min ( wall_dist, wall - ptest(i,itest) )

      end if

    end do

    first_dist = dist_min(itest)

    if ( first_dist < wall_dist ) then
      cycle
    end if

    do i = 1, ndim

      if ( 1 < searched_lo(i) ) then
!
!  Solve for SEARCH_NEW(I) so that FIRST_DIST < PTEST(I,ITEST) - WALL.
!
        w = ( ptest(i,itest) - first_dist - bin_min(i) ) &
          / real ( bin_max(i) - bin_min(i), kind = 8 )
        k = int ( real ( nbin(i),  kind = 8 ) * w )

        k = max ( k, 1 )

        search1(1:ndim) = searched_lo(1:ndim)
        search1(i) = k
        search2(1:ndim) = searched_hi(1:ndim)
        search2(i) = searched_lo(i) - 1

        if ( search1(i) <= search2(i) ) then

          rank = 0

          do

            call tuple_next2 ( ndim, search1, search2, bin, rank )

            if ( rank == 0 ) then
              exit
            end if

            call bin_search_one_2d ( bin, nset, pset, nbin, bin_start, &
              bin_next, ptest(1,itest), found_a_neighbor, i_min(itest), &
              dist_min_sq, compares(itest) )

          end do

          searched_lo(i) = k

        end if

      end if

      if ( searched_hi(i) < nbin(i) ) then
!
!  Solve for SEARCH_NEW(I) so that FIRST_DIST < WALL - PTEST(I,ITEST).
!
        w = ( ptest(i,itest) + first_dist - bin_min(i) ) &
          / real ( bin_max(i) - bin_min(i), kind = 8 )
        k = 1 + int ( real ( nbin(i), kind = 8 ) * w )

        k = min ( k, nbin(i) )

        search1(1:ndim) = searched_lo(1:ndim)
        search1(i) = searched_hi(i)+1
        search2(1:ndim) = searched_hi(1:ndim)
        search2(i) = k

        if ( search1(i) <= search2(i) ) then

          rank = 0

          do

            call tuple_next2 ( ndim, search1, search2, bin, rank )

            if ( rank == 0 ) then
              exit
            end if

            call bin_search_one_2d ( bin, nset, pset, nbin, bin_start, &
              bin_next, ptest(1,itest), found_a_neighbor, i_min(itest), &
              dist_min_sq, compares(itest) )

          end do

          searched_hi(i) = k

        end if

      end if

    end do

    dist_min(itest) = sqrt ( dist_min_sq )

  end do

  return
end
subroutine points_nearest_points_bins_3d_2 ( nset, pset, nbin, bin_min, &
  bin_max, bin_start, bin_last, bin_next, ntest, ptest, i_min, dist_min, &
  compares )

!*****************************************************************************80
!
!! POINTS_NEAREST_POINTS_BINS_3D_2 finds the nearest point to given points in 3D.
!
!  Discussion:
!
!    A set of NSET points with coordinates PSET is assumed to lie within a
!    rectangle.  The limits of this rectangle are given in BIN_MIN and BIN_MAX.
!    The rectangle is divided up into an NBIN by NBIN regular grid of cells.
!
!    The cells are ordered lexicographically, as suggested by the following
!    diagram for NBIN = 5:
!
!     *----*----*----*----*----*
!     |    |    |    |    |    |
!     | 21 | 22 | 23 | 24 | 25 |
!     |    |    |    |    |    |
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
!    15 February 2001
!
!  Author:
!
!    John Burkardt
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
!    Input, integer ( kind = 4 ) NSET, the number of points in the set.
!
!    Input, real ( kind = 8 ) PSET(3,NSET), the points in the set.
!
!    Input, integer ( kind = 4 ) NBIN, the number of cells.  NBIN must be at
!    least 3.
!
!    Input, real ( kind = 8 ) BIN_MIN(3), BIN_MAX(3), the minimum and
!    maximum bin values.
!
!    Input, integer ( kind = 4 ) BIN_START(NBIN,NBIN,NBIN), 
!    BIN_LAST(NBIN,NBIN,NBIN), the index of the first and last element in 
!    the bin, or -1 if none.
!
!    Input, integer ( kind = 4 ) BIN_NEXT(NSET), the index of the next element 
!    of the bin containing this element.
!
!    Input, integer ( kind = 4 ) NTEST, the number of test points.
!
!    Input, real ( kind = 8 ) PTEST(3,NTEST), the test points.
!
!    Output, integer ( kind = 4 ) I_MIN(NTEST), the index of the nearest point 
!    in PSET to PTEST(ITEST).
!
!    Output, real ( kind = 8 ) DIST_MIN(NTEST), the distance between
!    PTEST(*,ITEST) and PSET(*,I_MIN).
!
!    Output, integer ( kind = 4 ) COMPARES(NTEST), the number of point-to-point 
!    comparisons.
!
  implicit none

  integer ( kind = 4 ) nbin
  integer ( kind = 4 ), parameter :: ndim = 3
  integer ( kind = 4 ) nset
  integer ( kind = 4 ) ntest

  integer ( kind = 4 ) bin(ndim)
  integer ( kind = 4 ) bin_last(nbin,nbin,nbin)
  real ( kind = 8 ) bin_max(ndim)
  real ( kind = 8 ) bin_min(ndim)
  integer ( kind = 4 ) bin_start(nbin,nbin,nbin)
  integer ( kind = 4 ) bin_next(nset)
  real ( kind = 8 ) c_max(ndim)
  real ( kind = 8 ) c_min(ndim)
  integer ( kind = 4 ) compares(ntest)
  real ( kind = 8 ) dist_min(ntest)
  real ( kind = 8 ) dist_min_sq
  real ( kind = 8 ) dist_sq
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_min(ntest)
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) itest
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jc
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kc
  integer ( kind = 4 ) layer
  real ( kind = 8 ) layer_width
  logical more_bins
  integer ( kind = 4 ) node
  real ( kind = 8 ) pset(ndim,nset)
  real ( kind = 8 ) ptest(ndim,ntest)
  real ( kind = 8 ) search_radius

  compares(1:ntest) = 0
!
!  Special cases.
!
  if ( nset <= 0 ) then
    dist_min(1:ntest) = huge ( dist_min )
    i_min(1:ntest) = 0
    return
  end if

  if ( nset == 1 ) then
    do itest = 1, ntest
      dist_min(itest) = sqrt ( &
        sum ( ( ptest(1:ndim,itest) - pset(1:ndim,1) )**2 ) )
    end do
    compares(1:ntest) = 1
    i_min(1:ntest) = 1
    return
  end if

  layer_width = minval ( &
    ( bin_max(1:ndim) - bin_min(1:ndim) ) / real ( nbin, kind = 8 ) &
  )
!
!  Search for each of the NTEST points.
!
  do itest = 1, ntest

    dist_min_sq = huge ( dist_min_sq )
    i_min(itest) = 0
    search_radius = 0.0D+00
!
!  Determine the bin coordinates of the point P.
!
    call r83_to_bin_even2 ( nbin, bin_min, bin_max, ptest(1,itest), bin )
!
!  Determine the radius of the ball of space that will be completely
!  searched after this first bin, "layer 0", is done.
!
    call bin_to_r83_even2 ( nbin, bin, bin_min, bin_max, c_min, c_max )
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

        if ( 1 <= i .and. i <= nbin .and. &
             1 <= j .and. j <= nbin .and. &
             1 <= k .and. k <= nbin ) then

          node = bin_start(i,j,k)

          do while ( 0 < node )

            dist_sq = sum ( ( ptest(1:ndim,itest) - pset(1:ndim,node) )**2 )
            compares(itest) = compares(itest) + 1

            if ( dist_sq < dist_min_sq ) then
              dist_min_sq = dist_sq
              i_min(itest) = node
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

          call index_box2_next_3d ( layer, layer, layer, ic, jc, kc, &
            i, j, k, more_bins )

          if ( .not. more_bins ) then
            exit
          end if

          if ( 1 <= i .and. i <= nbin .and. &
               1 <= j .and. j <= nbin .and. &
               1 <= k .and. k <= nbin ) then
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
          minval ( abs ( ptest(1:ndim,itest) - c_min(1:ndim) ) ), &
          minval ( abs ( ptest(1:ndim,itest) - c_max(1:ndim) ) ) )
      else
        search_radius = search_radius + layer_width
      end if
!
!  We are done with PTEST(*,ITEST) if:
!
!    * We've found at least one neighbor;
!    AND
!    * We've searched all legal boxes that contain the circle
!      with PTEST at the center and the nearest neighbor on the circumference.
!
     if ( i_min(itest) /= 0 ) then
       dist_min(itest) = sqrt ( dist_min_sq )
       if ( dist_min(itest) <= search_radius ) then
         exit
       end if
     end if
!
!  Prepare to search the next layer.
!
      layer = layer + 1

    end do

  end do

  return
end
subroutine points_nearest_points_bins_3d_3 ( nset, pset, nbin, bin_min, &
  bin_max, bin_start, bin_last, bin_next, ntest, ptest, i_min, dist_min, &
  compares )

!*****************************************************************************80
!
!! POINTS_NEAREST_POINTS_BINS_3D_3 finds the nearest point to given points in 3D.
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
!      * P1 is in the same cell as P2, P1.X = P2.X, P1.Y = P2.Y, but P1.Z < P2.Z
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
!    John Burkardt
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
!    Input, integer ( kind = 4 ) NSET, the number of points in the set.
!
!    Input, real ( kind = 8 ) PSET(3,NSET), the points in the set.
!
!    Input, integer ( kind = 4 ) NBIN(3), the number of cells in the X, Y
!    and Z directions.
!
!    Input, real ( kind = 8 ) BIN_MIN(3), BIN_MAX(3), the minimum and
!    maximum bin values.
!
!    Input, integer ( kind = 4 ) BIN_START(NBIN(1),NBIN(2),NBIN(3)),
!    BIN_LAST(NBIN(1),NBIN(2),NBIN(3)),
!    indicates the index of the first and last element in the bin,
!    or -1 if none.
!
!    Input, integer ( kind = 4 ) BIN_NEXT(NSET), the index of the next element 
!    of the bin containing this element.
!
!    Input, integer ( kind = 4 ) NTEST, the number of test points.
!
!    Input, real ( kind = 8 ) PTEST(3,NTEST), the test points.
!
!    Output, integer ( kind = 4 ) I_MIN(NTEST), the index of the nearest point 
!    in PSET to PTEST(ITEST).
!
!    Output, real ( kind = 8 ) DIST_MIN(NTEST), the distance between
!    PTEST(*,ITEST) and PSET(*,I_MIN).
!
!    Output, integer ( kind = 4 ) COMPARES(NTEST), the number of point-to-point 
!    comparisons.
!
  implicit none

  integer ( kind = 4 ), parameter :: ndim = 3

  integer ( kind = 4 ) nbin(ndim)
  integer ( kind = 4 ) nset
  integer ( kind = 4 ) ntest

  integer ( kind = 4 ) bin(ndim)
  integer ( kind = 4 ) bin_last(nbin(1),nbin(2),nbin(3))
  real ( kind = 8 ) bin_max(ndim)
  real ( kind = 8 ) bin_min(ndim)
  integer ( kind = 4 ) bin_start(nbin(1),nbin(2),nbin(3))
  integer ( kind = 4 ) bin_next(nset)
  real ( kind = 8 ) c_max(ndim)
  real ( kind = 8 ) c_min(ndim)
  integer ( kind = 4 ) compares(ntest)
  real ( kind = 8 ) dist_min(ntest)
  real ( kind = 8 ) dist_min_sq
  real ( kind = 8 ) dist_sq
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_min(ntest)
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) itest
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jc
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kc
  integer ( kind = 4 ) layer
  real ( kind = 8 ) layer_width
  logical more_bins
  integer ( kind = 4 ) node
  real ( kind = 8 ) pset(ndim,nset)
  real ( kind = 8 ) ptest(ndim,ntest)
  real ( kind = 8 ) search_radius

  compares(1:ntest) = 0
!
!  Special cases.
!
  if ( nset <= 0 ) then
    dist_min(1:ntest) = huge ( dist_min(1) )
    i_min(1:ntest) = 0
    return
  end if

  if ( nset == 1 ) then
    do itest = 1, ntest
      dist_min(itest) = sqrt ( &
        sum ( ( ptest(1:ndim,itest) - pset(1:ndim,1) )**2 ) )
    end do
    compares(1:ntest) = 1
    i_min(1:ntest) = 1
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
!
!  Search for each of the NTEST points.
!
  do itest = 1, ntest

    dist_min_sq = huge ( dist_min_sq )
    i_min(itest) = 0
    search_radius = 0.0D+00
!
!  Determine the bin coordinates of the point P.
!
    call r83_to_bin_even3 ( nbin, bin_min, bin_max, ptest(1,itest), bin )
!
!  Determine the radius of the ball of space that will be completely
!  searched after this first bin, "layer 0", is done.
!
    call bin_to_r83_even3 ( nbin, bin, bin_min, bin_max, c_min, c_max )
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

            dist_sq = sum ( ( ptest(1:ndim,itest) - pset(1:ndim,node) )**2 )
            compares(itest) = compares(itest) + 1

            if ( dist_sq < dist_min_sq ) then
              dist_min_sq = dist_sq
              i_min(itest) = node
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
          minval ( abs ( ptest(1:ndim,itest) - c_min(1:ndim) ) ), &
          minval ( abs ( ptest(1:ndim,itest) - c_max(1:ndim) ) ) )
      else
        search_radius = search_radius + layer_width
      end if
!
!  We are done with PTEST(*,ITEST) if:
!
!    * We've found at least one neighbor;
!    AND
!    * We've searched all legal boxes that contain the circle
!      with PTEST at the center and the nearest neighbor on the circumference.
!
      if ( i_min(itest) /= 0 ) then
        dist_min(itest) = sqrt ( dist_min_sq )
        if ( dist_min(itest) <= search_radius ) then
          exit
        end if
      end if

      layer = layer + 1

    end do
!
!  We are now done with all the layers.
!
  end do
!
!  We are now done with all the test points.
!
  return
end
subroutine points_nearest_points_naive_2d ( nset, pset, ntest, ptest, i_min, &
  dist_min )

!*****************************************************************************80
!
!! POINTS_NEAREST_POINTS_NAIVE_2D finds the nearest point to given points in 2D.
!
!  Discussion:
!
!    A naive algorithm is used.  The distance to every point is calculated,
!    in order to determine the smallest.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NSET, the number of points in the set.
!
!    Input, real ( kind = 8 ) PSET(2,NSET), the points in the set.
!
!    Input, integer ( kind = 4 ) NTEST, the number of test points.
!
!    Input, real ( kind = 8 ) PTEST(2,NTEST), the test points.
!
!    Output, integer ( kind = 4 ) I_MIN(NTEST), the index of the nearest
!    point in PSET to PTEST(ITEST).
!
!    Output, real ( kind = 8 ) DIST_MIN(NTEST), the distance between 
!    PTEST(*,ITEST) and PSET(*,I_MIN).
!
  implicit none

  integer ( kind = 4 ) nset
  integer ( kind = 4 ), parameter :: ndim = 2
  integer ( kind = 4 ) ntest

  real ( kind = 8 ) d
  real ( kind = 8 ) dist_min(ntest)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_min(ntest)
  integer ( kind = 4 ) itest
  real ( kind = 8 ) pset(ndim,nset)
  real ( kind = 8 ) ptest(ndim,ntest)

  do itest = 1, ntest

    dist_min(itest) = huge ( dist_min )
    i_min(itest) = -1

    do i = 1, nset
      d = sum ( ( ptest(1:ndim,itest) - pset(1:ndim,i) )**2 )
      if ( d < dist_min(itest) ) then
        dist_min(itest) = d
        i_min(itest) = i
      end if
    end do

    dist_min(itest) = sqrt ( dist_min(itest) )

  end do

  return
end
subroutine points_nearest_points_naive_3d ( nset, pset, ntest, ptest, i_min, &
  dist_min )

!*****************************************************************************80
!
!! POINTS_NEAREST_POINTS_NAIVE_3D finds the nearest point to given points in 3D.
!
!  Discussion:
!
!    A naive algorithm is used.  The distance to every point is calculated,
!    in order to determine the smallest.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NSET, the number of points in the set.
!
!    Input, real ( kind = 8 ) PSET(3,NSET), the points in the set.
!
!    Input, integer ( kind = 4 ) NTEST, the number of test points.
!
!    Input, real ( kind = 8 ) PTEST(3,NTEST), the test points.
!
!    Output, integer ( kind = 4 ) I_MIN(NTEST), the index of the nearest point
!    in PSET to PTEST(ITEST).
!
!    Output, real ( kind = 8 ) DIST_MIN(NTEST), the distance between 
!    PTEST(*,ITEST) and PSET(*,I_MIN).
!
  implicit none

  integer ( kind = 4 ) nset
  integer ( kind = 4 ), parameter :: ndim = 3
  integer ( kind = 4 ) ntest

  real ( kind = 8 ) d
  real ( kind = 8 ) dist_min(ntest)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_min(ntest)
  integer ( kind = 4 ) itest
  real ( kind = 8 ) pset(ndim,nset)
  real ( kind = 8 ) ptest(ndim,ntest)

  do itest = 1, ntest

    dist_min(itest) = huge ( dist_min )
    i_min(itest) = -1

    do i = 1, nset
      d = sum ( ( ptest(1:ndim,itest) - pset(1:ndim,i) )**2 )
      if ( d < dist_min(itest) ) then
        dist_min(itest) = d
        i_min(itest) = i
      end if
    end do

    dist_min(itest) = sqrt ( dist_min(itest) )

  end do

  return
end
subroutine r8_to_bin_even ( nbin, a, b, c, bin )

!*****************************************************************************80
!
!! R8_TO_BIN_EVEN determines the appropriate "bin" for C in [A,B].
!
!  Discussion:
!
!    The interval from A to B is divided into NBIN-2 equal subintervals or bins.
!    An initial bin takes everything less than A, and a final bin takes
!    everything greater than B.
!
!  Example:
!
!    NBIN = 7, A = 5, B = 15
!
!    C   BIN
!
!    1    1
!    3    1
!    4.9  1
!    5    2
!    6    2
!    7    3
!    8    3
!    9.5  4
!   13    6
!   14    6
!   15    6
!   15.1  7
!   99    7
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NBIN, the number of bins.  NBIN is normally at 
!    least 3.  If NBIN is 1 or 2, then everything is assigned to bin 1.
!
!    Input, real ( kind = 8 ) A, B, the lower and upper limits of the bin
!    interval.  While A is expected to be less than B, the code should
!    return useful results if A is actually greater than B.
!
!    Input, real ( kind = 8 ) C, a value to be placed in a bin.
!
!    Output, integer ( kind = 4 ) BIN, the index of the bin to which C is 
!    assigned.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) a2
  real ( kind = 8 ) b
  real ( kind = 8 ) b2
  integer ( kind = 4 ) bin
  real ( kind = 8 ) c
  integer ( kind = 4 ) nbin
  logical switch
!
!  Take care of special cases.
!
  if ( nbin < 1 ) then
    bin = 0
    return
  else if ( nbin == 1 .or. nbin == 2 ) then
    bin = 1
    return
  end if

  if ( b == a ) then
    bin = 0
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
  if ( c < a2 ) then
    bin = 1
  else if ( c == a2 ) then
    bin = 2
  else if ( c == b2 ) then
    bin = nbin - 1
  else if ( b2 < c ) then
    bin = nbin
  else
    bin = 2 + int ( real ( nbin - 2, kind = 8 ) * ( c - a2 ) / ( b2 - a2 ) )
    bin = max ( bin, 2 )
    bin = min ( bin, nbin-1 )
  end if
!
!  Reverse the switching.
!
  if ( switch ) then
    bin = nbin + 1 - bin
  end if

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
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NBIN, the number of bins.
!
!    Input, real ( kind = 8 ) A, B, the lower and upper limits of the bin
!    interval.  While A is expected to be less than B, the code should
!    return useful results if A is actually greater than B.
!
!    Input, real ( kind = 8 ) C, a value to be placed in a bin.
!
!    Output, integer ( kind = 4 ) BIN, the index of the bin to which C 
!    is assigned.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) a2
  real ( kind = 8 ) b
  real ( kind = 8 ) b2
  integer ( kind = 4 ) bin
  real ( kind = 8 ) c
  integer ( kind = 4 ) nbin
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
subroutine r8_to_bin_uneven ( nbin, xbin, x, bin )

!*****************************************************************************80
!
!! R8_TO_BIN_UNEVEN places X in one of several unevenly spaced bins.
!
!  Discussion:
!
!    The XBIN array is assumed to be sorted.
!
!  Example:
!
!    NBIN = 5
!    XBIN(1:4) = (/ 0.0, 2.0, 8.0, 9.0 /)
!
!    so bins are
!
!    1  ( -Inf,   0 )
!    2  (    0,   2 )
!    3  (    2,   8 )
!    4  (    8,   9 )
!    5  (    9, Inf )
!
!    X   BIN
!
!   -7    1
!   -3    1
!    0    1
!    0.1  2
!    1    2
!    3    3
!    8    3
!    9.5  5
!   13    5
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NBIN, the number of bins.
!
!    Input, real ( kind = 8 ) XBIN(NBIN-1), the dividing values for the bins.
!
!    Input, real ( kind = 8 ) X, a value to be placed in a bin.
!
!    Output, integer ( kind = 4 ) BIN, the index of the bin to which X 
!    is assigned.
!
  implicit none

  integer ( kind = 4 ) nbin

  integer ( kind = 4 ) bin
  real ( kind = 8 ) x
  real ( kind = 8 ) xbin(nbin-1)

  bin = 1

  do while ( bin < nbin )

    if ( x <= xbin(bin) ) then
      return
    end if

    bin = bin + 1

  end do

  return
end
function r8_uniform ( a, b, seed )

!*****************************************************************************80
!
!! R8_UNIFORM returns a scaled pseudorandom R8.
!
!  Discussion:
!
!    An R8 is a real ( kind = 8 ) value.
!
!    For now, the input quantity SEED is an integer ( kind = 4 ) variable.
!
!    The pseudorandom number should be uniformly distributed
!    between A and B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the limits of the interval.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which should
!    NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R8_UNIFORM, a number strictly between A and B.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) k
  real ( kind = 8 ) r8_uniform
  integer ( kind = 4 )seed

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_UNIFORM - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + 2147483647
  end if

  r8_uniform = a + ( b - a ) * real ( seed, kind = 8 ) * 4.656612875D-10

  return
end
function r8_uniform_01 ( seed )

!*****************************************************************************80
!
!! R8_UNIFORM_01 returns a unit pseudorandom R8.
!
!  Discussion:
!
!    An R8 is a real ( kind = 8 ) value.
!
!    For now, the input quantity SEED is an integer ( kind = 4 ) variable.
!
!    This routine implements the recursion
!
!      seed = 16807 * seed mod ( 2**31 - 1 )
!      r8_uniform_01 = seed / ( 2**31 - 1 )
!
!    The integer arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!    If the initial seed is 12345, then the first three computations are
!
!      Input     Output      R8_UNIFORM_01
!      SEED      SEED
!
!         12345   207482415  0.096616
!     207482415  1790989824  0.833995
!    1790989824  2035175616  0.947702
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 July 2006
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
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which should
!    NOT be 0. On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R8_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer ( kind = 4 ) k
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + 2147483647
  end if
!
!  Although SEED can be represented exactly as a 32 bit integer,
!  it generally cannot be represented exactly as a 32 bit real number!
!
  r8_uniform_01 = real ( seed, kind = 8 ) * 4.656612875D-10

  return
end
subroutine r82_to_bin_even ( nbin, a, b, c, bin )

!*****************************************************************************80
!
!! R82_TO_BIN_EVEN determines the appropriate "bin" for an R82 value.
!
!  Discussion:
!
!    The intervals [A(1),B(1)] and [A(2),B(2)] are each divided into NBIN-2
!    equal subintervals or bins.  Boundary bins take care of extreme values.
!
!  Example:
!
!    NBIN = 7, A(1) = 5,  A(2) = 0,
!              B(1) = 15, B(2) = 20.
!
!      C      BIN
!   ------  ------
!    8 -2    3  1
!    0  1    1  2
!    6  9    2  4
!   10 11    4  4
!   14 23    6  7
!   25 13    7  5
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NBIN, the number of bins in each dimension.
!    NBIN is normally at least 3.  If NBIN is 1 or 2, then everything
!    is assigned to bin 1.
!
!    Input, real ( kind = 8 ) A(2), B(2), the lower and upper limits of
!    the bin interval.  While A(I) is expected to be less than B(I), the
!    code should return useful results if A(I) is actually greater than B(I).
!
!    Input, real ( kind = 8 ) C(2), a value to be placed in a bin.
!
!    Output, integer ( kind = 4 ) BIN(2), the index of the bin to which C
!    is assigned.
!
  implicit none

  integer ( kind = 4 ), parameter :: ndim = 2

  real ( kind = 8 ) a(ndim)
  real ( kind = 8 ) b(ndim)
  integer ( kind = 4 ) bin(ndim)
  real ( kind = 8 ) c(ndim)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) nbin

  do i = 1, ndim
    call r8_to_bin_even ( nbin, a(i), b(i), c(i), bin(i) )
  end do

  return
end
subroutine r82_to_bin_even2 ( nbin, a, b, c, bin )

!*****************************************************************************80
!
!! R82_TO_BIN_EVEN2 determines the appropriate "bin" for an R82 value.
!
!  Discussion:
!
!    The intervals [A(1),B(1)] and [A(2),B(2)] are each divided into NBIN
!    equal subintervals or bins.  Boundary bins take care of extreme values.
!
!  Example:
!
!    NBIN = 5, A(1) = 5,  A(2) = 0,
!              B(1) = 15, B(2) = 20.
!
!   20 +    +    +    +    +    +
!        15 | 25 | 35 | 45 | 55
!   16 +----+----+----+----+----+
!        14 | 24 | 34 | 44 | 54
!   12 +----+----+----+----+----+
!        13 | 23 | 33 | 43 | 53
!    8 +----+----+----+----+----+
!        12 | 22 | 32 | 42 | 52
!    4 +----+----+----+----+----+
!        11 | 21 | 31 | 41 | 51
!    0 +    +    +    +    +    +
!      5    7    9   11   13   15
!
!      C      BIN
!   ------  ------
!    8 -2    2  1
!    0  1    1  1
!    6  9    1  3
!   10 11    3  3
!   14 23    5  5
!   25 13    5  4
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
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NBIN, the number of bins in each dimension.
!    NBIN must be at least 1.
!
!    Input, real ( kind = 8 ) A(2), B(2), the lower and upper limits of the
!    bin interval.  While A(I) is expected to be less than B(I), the code
!    should return useful results if A(I) is actually greater than B(I).
!
!    Input, real ( kind = 8 ) C(2), a value to be placed in a bin.
!
!    Output, integer ( kind = 4 ) BIN(2), the index of the bin to which C 
!    is assigned.
!
  implicit none

  integer ( kind = 4 ), parameter :: ndim = 2

  real ( kind = 8 ) a(ndim)
  real ( kind = 8 ) b(ndim)
  integer ( kind = 4 ) bin(ndim)
  real ( kind = 8 ) c(ndim)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) nbin

  do i = 1, ndim
    call r8_to_bin_even2 ( nbin, a(i), b(i), c(i), bin(i) )
  end do

  return
end
subroutine r82_to_bin_even3 ( nbin, a, b, c, bin )

!*****************************************************************************80
!
!! R82_TO_BIN_EVEN3 determines the appropriate "bin" for an R82 value.
!
!  Discussion:
!
!    The interval [A(I),B(I)] is divided into NBIN(I) equal subintervals
!    or bins.
!
!  Example:
!
!    NBIN = (/ 4, 5 /),
!
!      A(1) = 1,  A(2) = 0,
!      B(1) = 17, B(2) = 20.
!
!   20 +    +    +    +    +
!        15 | 25 | 35 | 45
!   16 +----+----+----+----+
!        14 | 24 | 34 | 44
!   12 +----+----+----+----+
!        13 | 23 | 33 | 43
!    8 +----+----+----+----+
!        12 | 22 | 32 | 42
!    4 +----+----+----+----+
!        11 | 21 | 31 | 41
!    0 +    +    +    +    +
!      1    5    9   13   17
!
!      C      BIN
!   ------  ------
!    8 -2    2  1
!    0  1    1  1
!    6  9    2  3
!   10 11    3  3
!   14 23    4  5
!   25 13    4  4
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NBIN(2), the number of bins in each dimension.
!
!    Input, real ( kind = 8 ) A(2), B(2), the lower and upper limits of the
!    bin interval.  While A(I) is expected to be less than B(I), the code
!    should return useful results if A(I) is actually greater than B(I).
!
!    Input, real ( kind = 8 ) C(2), a value to be placed in a bin.
!
!    Output, integer ( kind = 4 ) BIN(2), the index of the bin to which C is assigned.
!
  implicit none

  integer ( kind = 4 ), parameter :: ndim = 2

  real ( kind = 8 ) a(ndim)
  real ( kind = 8 ) b(ndim)
  integer ( kind = 4 ) bin(ndim)
  real ( kind = 8 ) c(ndim)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) nbin(ndim)

  do i = 1, ndim
    call r8_to_bin_even2 ( nbin(i), a(i), b(i), c(i), bin(i) )
  end do

  return
end
subroutine r82_uniform ( rlo, rhi, seed, r )

!*****************************************************************************80
!
!! R82_UNIFORM returns a random R82 value in a given range.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) RLO(2), RHI(2), the minimum and maximum values.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
!    Output, real ( kind = 8 ) R(2), the randomly chosen value.
!
  implicit none

  integer ( kind = 4 ), parameter :: ndim = 2

  integer ( kind = 4 ) i
  real ( kind = 8 ) r(ndim)
  real ( kind = 8 ) r8_uniform
  real ( kind = 8 ) rhi(ndim)
  real ( kind = 8 ) rlo(ndim)
  integer ( kind = 4 ) seed

  do i = 1, ndim
    r(i) = r8_uniform ( rlo(i), rhi(i), seed )
  end do

  return
end
subroutine r82vec_bin_even ( n, a, nbin, bin_min, bin_max, bin_start, &
  bin_last, bin_next )

!*****************************************************************************80
!
!! R82VEC_BIN_EVEN bins an R82VEC into evenly spaced bins.
!
!  Discussion:
!
!    This is only a partial, indexed, sorting of the data.  To sort
!    the data, it is necessary to build a new array by extracting the
!    data for each bin, sorting that, and appending it to the array.
!
!    There are NBIN 1D bins in both X and Y directions, making a total
!    of NBIN**2 2D bins.  Each set of 1D bins begins and ends with an
!    "open-ended" bin that catches extreme values.
!
!    The 2D bins are indexed by the X and Y bins that construct them,
!    and ordered lexicographically by these indices:
!
!      1,4 | 2,4 | 3,4 | 4,4
!      ----+-----+-----+----
!      1,3 | 2,3 | 3,3 | 4,3
!      ----+-----+-----+----
!      1,2 | 2,2 | 3,2 | 4,2
!      ----+-----+-----+----
!      1,1 | 2,1 | 3,1 | 4,1
!
!    Thus, the 2D bin sequence is (1,1), (1,2), (1,3), (1,4), (2,1), ..., (4,4).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Thomas Cormen, Charles Leiserson, Ronald Rivest,
!    Introduction to Algorithms,
!    MIT Press, page 180.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of points in the data set.
!
!    Input, real ( kind = 8 ) A(2,N), the R82 data to be binned.
!
!    Input, integer ( kind = 4 ) NBIN, the (square root of) the number of bins.
!    NBIN must be at least 3.
!
!    Input, real ( kind = 8 ) BIN_MIN(2), BIN_MAX(2), the bin limits.
!
!    Output, integer ( kind = 4 ) BIN_START(NBIN,NBIN), BIN_LAST(NBIN,NBIN), the index
!    of the first and last elements of A that went into each bin, or
!    -1 if there were no entries in the bin.
!
!    Output, integer ( kind = 4 ) BIN_NEXT(N), the index of the next element of A
!    that follows this element in the same bin.  A value of 0 means this
!    is the last entry in the particular bin.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nbin
  integer ( kind = 4 ), parameter :: ndim = 2

  real ( kind = 8 ) a(ndim,n)
  integer ( kind = 4 ) bin(ndim)
  integer ( kind = 4 ) bin_last(nbin,nbin)
  integer ( kind = 4 ) bin_next(n)
  integer ( kind = 4 ) bin_start(nbin,nbin)
  real ( kind = 8 ) bin_max(ndim)
  real ( kind = 8 ) bin_min(ndim)
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
!
!  Initialize the bin pointers to -1.
!
  bin_last(1:nbin,1:nbin) = -1
  bin_start(1:nbin,1:nbin) = -1
  bin_next(1:n) = -1
!
!  Build up linked lists of entries that go into each bin.
!
  do j = 1, n

    call r82_to_bin_even ( nbin, bin_min, bin_max, a(1:ndim,j), bin )

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
subroutine r82vec_bin_even2 ( n, a, nbin, bin_min, bin_max, bin_start, &
  bin_last, bin_next )

!*****************************************************************************80
!
!! R82VEC_BIN_EVEN2 bins an R82VEC into evenly spaced bins.
!
!  Discussion:
!
!    This is only a partial, indexed, sorting of the data.  To sort
!    the data, it is necessary to build a new array by extracting the
!    data for each bin, sorting that, and appending it to the array.
!
!    There are NBIN 1D bins in both X and Y directions, making a total
!    of NBIN**2 2D bins.  Each set of 1D bins begins and ends at user
!    specified mininum and maximum values.
!
!    The 2D bins are indexed by the X and Y bins that construct them,
!    and ordered lexicographically by these indices:
!
!      1,4 | 2,4 | 3,4 | 4,4
!      ----+-----+-----+----
!      1,3 | 2,3 | 3,3 | 4,3
!      ----+-----+-----+----
!      1,2 | 2,2 | 3,2 | 4,2
!      ----+-----+-----+----
!      1,1 | 2,1 | 3,1 | 4,1
!
!    Thus, the 2D bin sequence is (1,1), (1,2), (1,3), (1,4), (2,1), ..., (4,4).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Thomas Cormen, Charles Leiserson, Ronald Rivest,
!    Introduction to Algorithms,
!    MIT Press, page 180.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of points in the data set.
!
!    Input, real ( kind = 8 ) A(2,N), the R82 data to be binned.
!
!    Input, integer ( kind = 4 ) NBIN, the (square root of) the number of bins.
!    NBIN must be at least 1.
!
!    Input, real ( kind = 8 ) BIN_MIN(2), BIN_MAX(2), the bin limits.
!
!    Output, integer ( kind = 4 ) BIN_START(NBIN,NBIN), BIN_LAST(NBIN,NBIN), the
!    index of the first and last elements of A that went into each bin,
!    or -1 if there are no entries in the bin.
!
!    Output, integer ( kind = 4 ) BIN_NEXT(N), the index of the next element of A
!    that follows this element in the same bin.  A value of 0 means this
!    is the last entry in the particular bin.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nbin
  integer ( kind = 4 ), parameter :: ndim = 2

  real ( kind = 8 ) a(ndim,n)
  integer ( kind = 4 ) bin(ndim)
  integer ( kind = 4 ) bin_last(nbin,nbin)
  integer ( kind = 4 ) bin_next(n)
  integer ( kind = 4 ) bin_start(nbin,nbin)
  real ( kind = 8 ) bin_max(ndim)
  real ( kind = 8 ) bin_min(ndim)
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
!
!  Initialize the bin pointers to -1.
!
  bin_last(1:nbin,1:nbin) = -1
  bin_start(1:nbin,1:nbin) = -1
  bin_next(1:n) = -1
!
!  Build up linked lists of entries that go into each bin.
!
  do j = 1, n

    call r82_to_bin_even2 ( nbin, bin_min, bin_max, a(1:ndim,j), bin )

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
!    Thus, the 2D bin sequence is (1,1), (1,2), (1,3), (1,4), (2,1), ..., (5,4).
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
!    John Burkardt
!
!  Reference:
!
!    Thomas Cormen, Charles Leiserson, Ronald Rivest,
!    Introduction to Algorithms,
!    MIT Press, page 180.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of points in the data set.
!
!    Input, real ( kind = 8 ) A(2,N), the R82 data to be binned.
!
!    Input, integer ( kind = 4 ) NBIN(2), the number of bins in each dimension.
!    NBIN must be at least 1.
!
!    Input, real ( kind = 8 ) BIN_MIN(2), BIN_MAX(2), the bin limits.
!
!    Output, integer ( kind = 4 )  BIN_START(NBIN(1),NBIN(2)), 
!    BIN_LAST(NBIN(1),NBIN(2)), the index of the first and last elements of A 
!    that went into each bin, or -1 if there are no entries in the bin.
!
!    Output, integer ( kind = 4 ) BIN_NEXT(N), the index of the next element 
!    of A that follows this element in the same bin.  A value of 0 means this
!    is the last entry in the particular bin.
!
  implicit none

  integer ( kind = 4 ), parameter :: ndim = 2

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nbin(ndim)

  real ( kind = 8 ) a(ndim,n)
  integer ( kind = 4 ) bin(ndim)
  integer ( kind = 4 ) bin_last(nbin(1),nbin(2))
  integer ( kind = 4 ) bin_next(n)
  integer ( kind = 4 ) bin_start(nbin(1),nbin(2))
  real ( kind = 8 ) bin_max(ndim)
  real ( kind = 8 ) bin_min(ndim)
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
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

    call r82_to_bin_even3 ( nbin, bin_min, bin_max, a(1:ndim,j), bin )

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
subroutine r82vec_binned_reorder ( n, a, nbin, bin_start, bin_last, bin_next )

!*****************************************************************************80
!
!! R82VEC_BINNED_REORDER reorders a binned R82 data vector.
!
!  Discussion:
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
!    09 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data points.
!
!    Input/output, real ( kind = 8 ) A(2,N), the R82 data to be sorted.
!
!    Input, integer ( kind = 4 ) NBIN, the (square root of the) number of bins.
!
!    Input/output, integer ( kind = 4 ) BIN_START(NBIN,NBIN), 
!    BIN_LAST(NBIN,NBIN), the index of the first and last element of A that 
!    went into each bin, or -1 if there are no entries in the bin.
!
!    Input/output, integer ( kind = 4 ) BIN_NEXT(N), the index of the next
!    element of A that follows this element in the same bin.  A value of 0
!    means this is the last entry in the particular bin.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nbin
  integer ( kind = 4 ), parameter :: ndim = 2

  real ( kind = 8 ) a(ndim,n)
  real ( kind = 8 ) a2(ndim,n)
  integer ( kind = 4 ) bin_last(nbin,nbin)
  integer ( kind = 4 ) bin_next(n)
  integer ( kind = 4 ) bin_start(nbin,nbin)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
!
!  Bin by bin, copy the contents of A to A2.
!  The BIN_START array is also updated as we go.
!
  k = 0

  do i1 = 1, nbin

    do i2 = 1, nbin

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

  do i1 = 1, nbin
    do i2 = 1, nbin

      k = bin_last(i1,i2)

      if ( 0 < k ) then
        bin_next(k) = 0
      end if

    end do
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
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data points.
!
!    Input/output, real ( kind = 8 ) A(2,N), the R82 data to be sorted.
!
!    Input, integer ( kind = 4 ) NBIN(2), the number of bins in each direction.
!
!    Input/output, integer ( kind = 4 ) BIN_START(NBIN(1),NBIN(2)),
!    BIN_LAST(NBIN(1),NBIN(2)), the index of the first and last element of A
!    that went into each bin, or -1 if there are no entries in the bin.
!
!    Input/output, integer ( kind = 4 ) BIN_NEXT(N), the index of the next 
!    element of A that follows this element in the same bin.  A value of 0 
!    means this is the last entry in the particular bin.
!
  implicit none

  integer ( kind = 4 ), parameter :: ndim = 2

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nbin(ndim)

  real ( kind = 8 ) a(ndim,n)
  real ( kind = 8 ) a2(ndim,n)
  integer ( kind = 4 ) bin_last(nbin(1),nbin(2))
  integer ( kind = 4 ) bin_next(n)
  integer ( kind = 4 ) bin_start(nbin(1),nbin(2))
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
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
subroutine r82vec_binned_sort_a ( n, a, nbin, bin_start, bin_last )

!*****************************************************************************80
!
!! R82VEC_BINNED_SORT_A sorts each bin of a binned R82VEC.
!
!  Discussion:
!
!    Presumably, the data vector was first binned by R82VEC_BIN_EVEN,
!    then reordered by R82VEC_BINNED_REORDER.  Now, each of the
!    bins of data is sorted one at a time.
!
!    The result is NOT a lexicographically sorted R82VEC.
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
!    09 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data points.
!
!    Input/output, real ( kind = 8 ) A(2,N), the R82 data to be sorted.
!
!    Input, integer ( kind = 4 ) NBIN, the (square root of the) number of bins.
!
!    Input, integer ( kind = 4 ) BIN_START(NBIN,NBIN), BIN_LAST(NBIN,NBIN), 
!    the index of the first and last element of A that went into each bin, or -1
!    if there are no entries in this bin.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nbin
  integer ( kind = 4 ), parameter :: ndim = 2

  real ( kind = 8 ) a(ndim,n)
  integer ( kind = 4 ) bin_last(nbin,nbin)
  integer ( kind = 4 ) bin_start(nbin,nbin)
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) n1

  do i1 = 1, nbin

    do i2 = 1, nbin

      j1 = bin_start(i1,i2)

      if ( 0 < j1 ) then

        j2 = bin_last(i1,i2)

        n1 = j2 + 1 - j1

        if ( 1 < n1 ) then
          call r82vec_sort_quick_a ( n1, a(1:ndim,j1:j2) )
        end if

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
!    The result is NOT a lexicographically sorted R82 vector.
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
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data points.
!
!    Input/output, real ( kind = 8 ) A(2,N), the R82 data to be sorted.
!
!    Input, integer ( kind = 4 ) NBIN(2), the number of bins in each dimension.
!
!    Input, integer ( kind = 4 ) BIN_START(NBIN(1),NBIN(2)),
!    BIN_LAST(NBIN(1),NBIN(2)), the index of the first and last elements of A
!    that went into each bin, or -1 if there are no entries in this bin.
!
  implicit none

  integer ( kind = 4 ), parameter :: ndim = 2

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nbin(ndim)

  real ( kind = 8 ) a(ndim,n)
  integer ( kind = 4 ) bin_last(nbin(1),nbin(2))
  integer ( kind = 4 ) bin_start(nbin(1),nbin(2))
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) n1

  do i1 = 1, nbin(1)

    do i2 = 1, nbin(2)

      j1 = bin_start(i1,i2)

      if ( 0 < j1 ) then

        j2 = bin_last(i1,i2)

        n1 = j2 + 1 - j1

        if ( 1 < n1 ) then
          call r82vec_sort_quick_a ( n1, a(1:ndim,j1:j2) )
        end if

      end if

    end do

  end do

  return
end
subroutine r82vec_part_quick_a ( n, a, l, r )

!*****************************************************************************80
!
!! R82VEC_PART_QUICK_A reorders an R82VEC as part of a quick sort.
!
!  Discussion:
!
!    The routine reorders the entries of A.  Using A(1:2,1) as a
!    key, all entries of A that are less than or equal to the key will
!    precede the key, which precedes all entries that are greater than the key.
!
!  Example:
!
!    Input:
!
!      N = 8
!
!      A = ( (2,4), (8,8), (6,2), (0,2), (10,6), (10,0), (0,6), (4,8) )
!
!    Output:
!
!      L = 2, R = 4
!
!      A = ( (0,2), (0,6), (2,4), (8,8), (6,2), (10,6), (10,0), (4,8) )
!             -----------          ----------------------------------
!             LEFT          KEY    RIGHT
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries of A.
!
!    Input/output, real ( kind = 8 ) A(2,N).  On input, the array to be checked.
!    On output, A has been reordered as described above.
!
!    Output, integer ( kind = 4 ) L, R, the indices of A that define the
!    three segments.  Let KEY = the input value of A(1:2,1).  Then
!    I <= L                 A(1:2,I) < KEY;
!         L < I < R         A(1:2,I) = KEY;
!                 R <= I    A(1:2,I) > KEY.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: ndim = 2

  real ( kind = 8 ) a(ndim,n)
  logical r8vec_eq
  logical r8vec_gt
  logical r8vec_lt
  integer ( kind = 4 ) i
  real ( kind = 8 ) key(ndim)
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) r

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R82VEC_PART_QUICK_A - Fatal error!'
    write ( *, '(a)' ) '  N < 1.'
    stop
  else if ( n == 1 ) then
    l = 0
    r = 2
    return
  end if

  key(1:ndim) = a(1:ndim,1)
  m = 1
!
!  The elements of unknown size have indices between L+1 and R-1.
!
  l = 1
  r = n + 1

  do i = 2, n

    if ( r8vec_gt ( ndim, a(1:ndim,l+1), key(1:ndim) ) ) then
      r = r - 1
      call r8vec_swap ( ndim, a(1:ndim,r), a(1:ndim,l+1) )
    else if ( r8vec_eq ( ndim, a(1:ndim,l+1), key(1:ndim) ) ) then
      m = m + 1
      call r8vec_swap ( ndim, a(1:ndim,m), a(1:ndim,l+1) )
      l = l + 1
    else if ( r8vec_lt ( ndim, a(1:ndim,l+1), key(1:ndim) ) ) then
      l = l + 1
    end if

  end do
!
!  Now shift small elements to the left, and KEY elements to center.
!
  do i = 1, l - m
    a(1:ndim,i) = a(1:ndim,i+m)
  end do

  l = l - m

  do i = 1, ndim
    a(i,l+1:l+m) = key(i)
  end do

  return
end
subroutine r82vec_permute ( n, a, p )

!*****************************************************************************80
!
!! R82VEC_PERMUTE permutes an R82VEC in place.
!
!  Discussion:
!
!    This routine permutes an array of real "objects", but the same
!    logic can be used to permute an array of objects of any arithmetic
!    type, or an array of objects of any complexity.  The only temporary
!    storage required is enough to store a single object.  The number
!    of data movements made is N + the number of cycles of order 2 or more,
!    which is never more than N + N/2.
!
!  Example:
!
!    Input:
!
!      N = 5
!      P = (   2,    4,    5,    1,    3 )
!      A = ( 1.0,  2.0,  3.0,  4.0,  5.0 )
!          (11.0, 22.0, 33.0, 44.0, 55.0 )
!
!    Output:
!
!      A    = (  2.0,  4.0,  5.0,  1.0,  3.0 )
!             ( 22.0, 44.0, 55.0, 11.0, 33.0 ).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 January 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of objects.
!
!    Input/output, real ( kind = 8 ) A(2,N), the array to be permuted.
!
!    Input, integer ( kind = 4 ) P(N), the permutation.  P(I) = J means
!    that the I-th element of the output array should be the J-th
!    element of the input array.  P must be a legal permutation
!    of the integers from 1 to N, otherwise the algorithm will
!    fail catastrophically.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(2,n)
  real ( kind = 8 ) a_temp(2)
  integer ( kind = 4 ) iget
  integer ( kind = 4 ) iput
  integer ( kind = 4 ) istart
  integer ( kind = 4 ) p(n)
!
!  Search for the next element of the permutation that has not been used.
!
  do istart = 1, n

    if ( p(istart) < 0 ) then

      cycle

    else if ( p(istart) == istart ) then

      p(istart) = - p(istart)
      cycle

    else

      a_temp(1:2) = a(1:2,istart)
      iget = istart
!
!  Copy the new value into the vacated entry.
!
      do

        iput = iget
        iget = p(iget)

        p(iput) = - p(iput)

        if ( iget < 1 .or. n < iget ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R82VEC_PERMUTE - Fatal error!'
          stop
        end if

        if ( iget == istart ) then
          a(1:2,iput) = a_temp(1:2)
          exit
        end if

        a(1:2,iput) = a(1:2,iget)

      end do

    end if

  end do
!
!  Restore the signs of the entries.
!
  p(1:n) = -p(1:n)

  return
end
subroutine r82vec_print ( n, a, title )

!*****************************************************************************80
!
!! R82VEC_PRINT prints an R82VEC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, real ( kind = 8 ) A(2,N), the R82 vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: ndim = 2

  real ( kind = 8 ) a(ndim,n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  if ( title /= ' ' ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(i6,(5g14.6))' ) i, a(1:ndim,i)
  end do

  return
end
subroutine r82vec_sort_heap_index_a ( n, a, indx )

!*****************************************************************************80
!
!! R82VEC_SORT_HEAP_INDEX_A does an indexed heap ascending sort of an R82VEC.
!
!  Discussion:
!
!    The sorting is not actually carried out.  Rather an index array is
!    created which defines the sorting.  This array may be used to sort
!    or index the array, or to sort or index related arrays keyed on the
!    original array.
!
!    Once the index array is computed, the sorting can be carried out
!    "implicitly:
!
!      A(1:2,INDX(I)), I = 1 to N is sorted,
!
!    or explicitly, by the call
!
!      call R82VEC_PERMUTE ( N, A, INDX )
!
!    after which A(1:2,I), I = 1 to N is sorted.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 January 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, real ( kind = 8 ) A(2,N), an array to be index-sorted.
!
!    Output, integer ( kind = 4 ) INDX(N), the sort index.  The
!    I-th element of the sorted array is A(1:2,INDX(I)).
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(2,n)
  real ( kind = 8 ) aval(2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) indxt
  integer ( kind = 4 ) ir
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l

  if ( n < 1 ) then
    return
  end if

  if ( n == 1 ) then
    indx(1) = 1
    return
  end if

  call i4vec_indicator ( n, indx )

  l = n / 2 + 1
  ir = n

  do

    if ( 1 < l ) then

      l = l - 1
      indxt = indx(l)
      aval(1:2) = a(1:2,indxt)

    else

      indxt = indx(ir)
      aval(1:2) = a(1:2,indxt)
      indx(ir) = indx(1)
      ir = ir - 1

      if ( ir == 1 ) then
        indx(1) = indxt
        exit
      end if

    end if

    i = l
    j = l + l

    do while ( j <= ir )

      if ( j < ir ) then
        if (   a(1,indx(j)) <  a(1,indx(j+1)) .or. &
             ( a(1,indx(j)) == a(1,indx(j+1)) .and. &
               a(2,indx(j)) <  a(2,indx(j+1)) ) ) then
          j = j + 1
        end if
      end if

      if (   aval(1) <  a(1,indx(j)) .or. &
           ( aval(1) == a(1,indx(j)) .and. &
             aval(2) <  a(2,indx(j)) ) ) then
        indx(i) = indx(j)
        i = j
        j = j + j
      else
        j = ir + 1
      end if

    end do

    indx(i) = indxt

  end do

  return
end
subroutine r82vec_sort_quick_a ( n, a )

!*****************************************************************************80
!
!! R82VEC_SORT_QUICK_A ascending sorts an R82VEC using quick sort.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input/output, real ( kind = 8 ) A(2,N).
!    On input, the array to be sorted.
!    On output, the array has been sorted.
!
  implicit none

  integer ( kind = 4 ), parameter :: maxlevel = 25
  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: ndim = 2

  real ( kind = 8 ) a(ndim,n)
  integer ( kind = 4 ) base
  integer ( kind = 4 ) l_segment
  integer ( kind = 4 ) level
  integer ( kind = 4 ) n_segment
  integer ( kind = 4 ) rsave(maxlevel)
  integer ( kind = 4 ) r_segment

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R82VEC_SORT_QUICK_A - Fatal error!'
    write ( *, '(a)' ) '  N < 1.'
    stop
  else if ( n == 1 ) then
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
    call r82vec_part_quick_a ( n_segment, a(1,base), l_segment, r_segment )
!
!  If the left segment has more than one element, we need to partition it.
!
    if ( 1 < l_segment ) then

      if ( maxlevel < level ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R82VEC_SORT_QUICK_A - Fatal error!'
        write ( *, '(a,i6)' ) '  Exceeding recursion maximum of ', maxlevel
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
subroutine r82vec_uniform ( n, alo, ahi, seed, a )

!*****************************************************************************80
!
!! R82VEC_UNIFORM returns a random R82VEC in a given range.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 April 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, real ( kind = 8 ) ALO(2), AHI(2), the minimum and maximum
!    values allowed for A(1,1:N) and A(2,1:N).
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number generator.
!
!    Output, real ( kind = 8 ) A(2,N), the vector of randomly chosen values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(2,n)
  real ( kind = 8 ) ahi(2)
  real ( kind = 8 ) alo(2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) r8_uniform
  integer ( kind = 4 ) seed

  do i = 1, 2
    do j = 1, n
      a(i,j) = r8_uniform ( alo(i), ahi(i), seed )
    end do
  end do

  return
end
subroutine r83_to_bin_even2 ( nbin, a, b, c, bin )

!*****************************************************************************80
!
!! R83_TO_BIN_EVEN2 determines the appropriate "bin" for an R83 value.
!
!  Discussion:
!
!    The intervals [A(I),B(I)] are each divided into NBIN
!    equal subintervals or bins.  Boundary bins take care of extreme values.
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
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NBIN, the number of bins in each dimension.
!    NBIN must be at least 1.
!
!    Input, real ( kind = 8 ) A(3), B(3), the lower and upper limits of the
!    bin interval.  While A(I) is expected to be less than B(I), the code
!    should return useful results if A(I) is actually greater than B(I).
!
!    Input, real ( kind = 8 ) C(3), a value to be placed in a bin.
!
!    Output, integer ( kind = 4 ) BIN(3), the index of the bin to which C
!    is assigned.
!
  implicit none

  integer ( kind = 4 ), parameter :: ndim = 3

  real ( kind = 8 ) a(ndim)
  real ( kind = 8 ) b(ndim)
  integer ( kind = 4 ) bin(ndim)
  real ( kind = 8 ) c(ndim)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) nbin

  do i = 1, ndim
    call r8_to_bin_even2 ( nbin, a(i), b(i), c(i), bin(i) )
  end do

  return
end
subroutine r83_to_bin_even3 ( nbin, a, b, c, bin )

!*****************************************************************************80
!
!! R83_TO_BIN_EVEN3 determines the appropriate "bin" for an R83 value.
!
!  Discussion:
!
!    The interval [A(I),B(I)] is divided into NBIN(I) equal subintervals
!    or bins.
!
!  Example:
!
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
!    05 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NBIN(3), the number of bins in each dimension.
!
!    Input, real ( kind = 8 ) A(3), B(3), the lower and upper limits of the
!    bin interval.  While A(I) is expected to be less than B(I), the code
!    should return useful results if A(I) is actually greater than B(I).
!
!    Input, real ( kind = 8 ) C(3), a value to be placed in a bin.
!
!    Output, integer ( kind = 4 ) BIN(3), the index of the bin to which C is assigned.
!
  implicit none

  integer ( kind = 4 ), parameter :: ndim = 3

  real ( kind = 8 ) a(ndim)
  real ( kind = 8 ) b(ndim)
  integer ( kind = 4 ) bin(ndim)
  real ( kind = 8 ) c(ndim)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) nbin(ndim)

  do i = 1, ndim
    call r8_to_bin_even2 ( nbin(i), a(i), b(i), c(i), bin(i) )
  end do

  return
end
subroutine r83vec_bin_even2 ( n, a, nbin, bin_min, bin_max, bin_start, &
  bin_last, bin_next )

!*****************************************************************************80
!
!! R83VEC_BIN_EVEN2 bins an R83VEC into evenly spaced bins.
!
!  Discussion:
!
!    This is only a partial, indexed, sorting of the data.  To sort
!    the data, it is necessary to build a new array by extracting the
!    data for each bin, sorting that, and appending it to the array.
!
!    There are NBIN 1D bins in each coordinate, making a total
!    of NBIN**NDIM bins.  Each set of 1D bins begins and ends at user
!    specified mininum and maximum values.
!
!    The bins are indexed by the 1D bins that construct them,
!    and ordered lexicographically by these indices:
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Thomas Cormen, Charles Leiserson, Ronald Rivest,
!    Introduction to Algorithms,
!    MIT Press, page 180.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of points in the data set.
!
!    Input, real ( kind = 8 ) A(3,N), the data to be binned.
!
!    Input, integer ( kind = 4 ) NBIN, the (cube root of) the number of bins.
!    NBIN must be at least 1.
!
!    Input, real ( kind = 8 ) BIN_MIN(3), BIN_MAX(3), the bin limits.
!
!    Output, integer ( kind = 4 ) BIN_START(NBIN,NBIN,NBIN), BIN_LAST(NBIN,NBIN,NBIN),
!    the index of the first and last elements of A that went into each bin,
!    or -1 if there are no entries in the bin.
!
!    Output, integer ( kind = 4 ) BIN_NEXT(N), the index of the next element of A
!    that follows this element in the same bin.  A value of 0 means this
!    is the last entry in the particular bin.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nbin
  integer ( kind = 4 ), parameter :: ndim = 3

  real ( kind = 8 ) a(ndim,n)
  integer ( kind = 4 ) bin(ndim)
  integer ( kind = 4 ) bin_last(nbin,nbin,nbin)
  integer ( kind = 4 ) bin_next(n)
  integer ( kind = 4 ) bin_start(nbin,nbin,nbin)
  real ( kind = 8 ) bin_max(ndim)
  real ( kind = 8 ) bin_min(ndim)
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
!
!  Initialize the bin pointers to -1.
!
  bin_last(1:nbin,1:nbin,1:nbin) = -1
  bin_start(1:nbin,1:nbin,1:nbin) = -1
  bin_next(1:n) = -1
!
!  Build up linked lists of entries that go into each bin.
!
  do j = 1, n

    call r83_to_bin_even2 ( nbin, bin_min, bin_max, a(1:ndim,j), bin )

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
subroutine r83vec_bin_even3 ( n, a, nbin, bin_min, bin_max, bin_start, &
  bin_last, bin_next )

!*****************************************************************************80
!
!! R83VEC_BIN_EVEN3 bins an R83 array into evenly spaced bins.
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
!    07 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Thomas Cormen, Charles Leiserson, Ronald Rivest,
!    Introduction to Algorithms,
!    MIT Press, page 180.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of points in the data set.
!
!    Input, real ( kind = 8 ) A(3,N), the R83 data to be binned.
!
!    Input, integer ( kind = 4 ) NBIN(3), the number of bins in each dimension.
!    NBIN must be at least 1.
!
!    Input, real ( kind = 8 ) BIN_MIN(3), BIN_MAX(3), the bin limits.
!
!    Output, integer ( kind = 4 ) BIN_START(NBIN(1),NBIN(2),NBIN(3)),
!    BIN_LAST(NBIN(1),NBIN(2),NBIN(3)), the
!    index of the first and last elements of A that went into each bin,
!    or -1 if there are no entries in the bin.
!
!    Output, integer ( kind = 4 ) BIN_NEXT(N), the index of the next element of A
!    that follows this element in the same bin.  A value of 0 means this
!    is the last entry in the particular bin.
!
  implicit none

  integer ( kind = 4 ), parameter :: ndim = 3

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nbin(ndim)

  real ( kind = 8 ) a(ndim,n)
  integer ( kind = 4 ) bin(ndim)
  integer ( kind = 4 ) bin_last(nbin(1),nbin(2),nbin(3))
  integer ( kind = 4 ) bin_next(n)
  integer ( kind = 4 ) bin_start(nbin(1),nbin(2),nbin(3))
  real ( kind = 8 ) bin_max(ndim)
  real ( kind = 8 ) bin_min(ndim)
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
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

    call r83_to_bin_even3 ( nbin, bin_min, bin_max, a(1:ndim,j), bin )

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
subroutine r83vec_binned_reorder ( n, a, nbin, bin_start, bin_last, bin_next )

!*****************************************************************************80
!
!! R83VEC_BINNED_REORDER reorders a binned R83 data vector.
!
!  Discussion:
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
!    09 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data points.
!
!    Input/output, real ( kind = 8 ) A(3,N), the data to be sorted.
!
!    Input, integer ( kind = 4 ) NBIN, the (cube root of the) number of bins.
!
!    Input/output, integer ( kind = 4 ) BIN_START(NBIN,NBIN,NBIN),
!    BIN_LAST(NBIN,NBIN,NBIN), the index of the first and last element of A
!    that went into each bin, or -1 if there are no entries in the bin.
!
!    Input/output, integer ( kind = 4 ) BIN_NEXT(N), the index of the next element of A
!    that follows this element in the same bin.  A value of 0 means this
!    is the last entry in the particular bin.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nbin
  integer ( kind = 4 ), parameter :: ndim = 3

  real ( kind = 8 ) a(ndim,n)
  real ( kind = 8 ) a2(ndim,n)
  integer ( kind = 4 ) bin_last(nbin,nbin,nbin)
  integer ( kind = 4 ) bin_next(n)
  integer ( kind = 4 ) bin_start(nbin,nbin,nbin)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
!
!  Bin by bin, copy the contents of A to A2.
!  The BIN_START array is also updated as we go.
!
  k = 0

  do i1 = 1, nbin

    do i2 = 1, nbin

      do i3 = 1, nbin

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
    bin_next(i) = i+1
  end do

  do i1 = 1, nbin
    do i2 = 1, nbin
      do i3 = 1, nbin
        k = bin_last(i1,i2,i3)

        if ( 0 < k ) then
          bin_next(k) = 0
        end if

      end do
    end do
  end do

  return
end
subroutine r83vec_binned_reorder2 ( n, a, nbin, bin_start, bin_last, bin_next )

!*****************************************************************************80
!
!! R83VEC_BINNED_REORDER2 reorders a binned R83 data vector.
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
!    30 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data points.
!
!    Input/output, real ( kind = 8 ) A(3,N), the data to be sorted.
!
!    Input, integer ( kind = 4 ) NBIN(3), the number of bins in each dimension.
!
!    Input/output, integer ( kind = 4 ) BIN_START(NBIN(1),NBIN(2),NBIN(3)),
!    BIN_LAST(NBIN(1),NBIN(2),NBIN(3)), the
!    index of the first and last element of A that went into each bin,
!    or -1 if there are no entries in the bin.
!
!    Input/output, integer ( kind = 4 ) BIN_NEXT(N), the index of the next element of A
!    that follows this element in the same bin.  A value of 0 means this
!    is the last entry in the particular bin.
!
  implicit none

  integer ( kind = 4 ), parameter :: ndim = 3

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nbin(ndim)

  real ( kind = 8 ) a(ndim,n)
  real ( kind = 8 ) a2(ndim,n)
  integer ( kind = 4 ) bin_last(nbin(1),nbin(2),nbin(3))
  integer ( kind = 4 ) bin_next(n)
  integer ( kind = 4 ) bin_start(nbin(1),nbin(2),nbin(3))
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
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
    bin_next(i) = i+1
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
subroutine r83vec_binned_sort_a ( n, a, nbin, bin_start, bin_last )

!*****************************************************************************80
!
!! R83VEC_BINNED_SORT_A sorts each bin of a binned R83VEC.
!
!  Discussion:
!
!    Presumably, the data vector was first binned by R83VEC_BIN_EVEN,
!    then reordered by R83VEC_BINNED_REORDER.  Now, each of the
!    bins of data is sorted one at a time.
!
!    The result is NOT a lexicographically sorted R83 vector.
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
!    09 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data points.
!
!    Input/output, real ( kind = 8 ) A(3,N), the data to be sorted.
!
!    Input, integer ( kind = 4 ) NBIN, the (cube root of the) number of bins.
!
!    Input, integer ( kind = 4 ) BIN_START(NBIN,NBIN,NBIN), BIN_LAST(NBIN,NBIN,NBIN),
!    the index of the first and last element of A that went into each bin,
!    or -1 if there are no entries in this bin.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nbin
  integer ( kind = 4 ), parameter :: ndim = 3

  real ( kind = 8 ) a(ndim,n)
  integer ( kind = 4 ) bin_last(nbin,nbin,nbin)
  integer ( kind = 4 ) bin_start(nbin,nbin,nbin)
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) n1

  do i1 = 1, nbin

    do i2 = 1, nbin

      do i3 = 1, nbin

        j1 = bin_start(i1,i2,i3)

        if ( 0 < j1 ) then

          j2 = bin_last(i1,i2,i3)

          n1 = j2 + 1 - j1

          if ( 1 < n1 ) then
            call r83vec_sort_quick_a ( n1, a(1:ndim,j1:j2) )
          end if

        end if

      end do

    end do

  end do

  return
end
subroutine r83vec_binned_sort_a2 ( n, a, nbin, bin_start, bin_last )

!*****************************************************************************80
!
!! R83VEC_BINNED_SORT_A2 sorts each bin of an R83 binned data vector.
!
!  Discussion:
!
!    This routine allows a different number of bins in each dimension.
!
!    Presumably, the data vector was first binned by R83VEC_BIN_EVEN3,
!    then reordered by R83VEC_BINNED_REORDER2.  Now, each of the
!    bins of data is sorted one at a time.
!
!    The result is NOT a lexicographically sorted R83 vector.
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
!    07 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data points.
!
!    Input/output, real ( kind = 8 ) A(3,N), the data to be sorted.
!
!    Input, integer ( kind = 4 ) NBIN(3), the number of bins in each dimension.
!
!    Input, integer ( kind = 4 ) BIN_START(NBIN(1),NBIN(2),NBIN(3)),
!    BIN_LAST(NBIN(1),NBIN(2),NBIN(3)),
!    the index of the first and last element of A that went into each bin,
!    or -1 if there are no entries in this bin.
!
  implicit none

  integer ( kind = 4 ), parameter :: ndim = 3

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nbin(ndim)

  real ( kind = 8 ) a(ndim,n)
  integer ( kind = 4 ) bin_last(nbin(1),nbin(2),nbin(3))
  integer ( kind = 4 ) bin_start(nbin(1),nbin(2),nbin(3))
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) n1

  if ( n <= 1 ) then
    return
  end if

  do i1 = 1, nbin(1)

    do i2 = 1, nbin(2)

      do i3 = 1, nbin(3)

        j1 = bin_start(i1,i2,i3)

        if ( 0 < j1 ) then

          j2 = bin_last(i1,i2,i3)

          n1 = j2 + 1 - j1

          if ( 1 < n1 ) then
            call r83vec_sort_quick_a ( n1, a(1:ndim,j1:j2) )
          end if

        end if

      end do

    end do

  end do

  return
end
subroutine r83vec_part_quick_a ( n, a, l, r )

!*****************************************************************************80
!
!! R83VEC_PART_QUICK_A reorders an R83VEC as part of a quick sort.
!
!  Discussion:
!
!    The routine reorders the entries of A.  Using A(1:3,1) as a
!    key, all entries of A that are less than or equal to the key will
!    precede the key, which precedes all entries that are greater than the key.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries of A.
!
!    Input/output, real ( kind = 8 ) A(3,N).  On input, the array to be checked.
!    On output, A has been reordered as described above.
!
!    Output, integer ( kind = 4 ) L, R, the indices of A that define the three segments.
!    Let KEY = the input value of A(1:3,1).  Then
!    I <= L                 A(1:3,I) < KEY;
!         L < I < R         A(1:3,I) = KEY;
!                 R <= I    A(1:3,I) > KEY.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: ndim = 3

  real ( kind = 8 ) a(ndim,n)
  logical r8vec_eq
  logical r8vec_gt
  logical r8vec_lt
  integer ( kind = 4 ) i
  real ( kind = 8 ) key(ndim)
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) r

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R83VEC_PART_QUICK_A - Fatal error!'
    write ( *, '(a)' ) '  N < 1.'
    stop
  else if ( n == 1 ) then
    l = 0
    r = 2
    return
  end if

  key(1:ndim) = a(1:ndim,1)
  m = 1
!
!  The elements of unknown size have indices between L+1 and R-1.
!
  l = 1
  r = n + 1

  do i = 2, n

    if ( r8vec_gt ( ndim, a(1:ndim,l+1), key(1:ndim) ) ) then
      r = r - 1
      call r8vec_swap ( ndim, a(1:ndim,r), a(1:ndim,l+1) )
    else if ( r8vec_eq ( ndim, a(1:ndim,l+1), key(1:ndim) ) ) then
      m = m + 1
      call r8vec_swap ( ndim, a(1:ndim,m), a(1:ndim,l+1) )
      l = l + 1
    else if ( r8vec_lt ( ndim, a(1:ndim,l+1), key(1:ndim) ) ) then
      l = l + 1
    end if

  end do
!
!  Now shift small elements to the left, and KEY elements to center.
!
  do i = 1, l - m
    a(1:ndim,i) = a(1:ndim,i+m)
  end do

  l = l - m

  do i = 1, ndim
    a(i,l+1:l+m) = key(i)
  end do

  return
end
subroutine r83vec_print ( n, a, title )

!*****************************************************************************80
!
!! R83VEC_PRINT prints an R83VEC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 September 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, real ( kind = 8 ) A(3,N), the R83 vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: ndim = 3

  real ( kind = 8 ) a(ndim,n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  if ( title /= ' ' ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(i6,(5g14.6))' ) i, a(1:ndim,i)
  end do

  return
end
subroutine r83vec_sort_quick_a ( n, a )

!*****************************************************************************80
!
!! R83VEC_SORT_QUICK_A ascending sorts an R83VEC using quick sort.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input/output, real ( kind = 8 ) A(3,N).
!    On input, the array to be sorted.
!    On output, the array has been sorted.
!
  implicit none

  integer ( kind = 4 ), parameter :: maxlevel = 25
  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: ndim = 3

  real ( kind = 8 ) a(ndim,n)
  integer ( kind = 4 ) base
  integer ( kind = 4 ) l_segment
  integer ( kind = 4 ) level
  integer ( kind = 4 ) n_segment
  integer ( kind = 4 ) rsave(maxlevel)
  integer ( kind = 4 ) r_segment

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R83VEC_SORT_QUICK_A - Fatal error!'
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
    call r83vec_part_quick_a ( n_segment, a(1,base), l_segment, r_segment )
!
!  If the left segment has more than one element, we need to partition it.
!
    if ( 1 < l_segment ) then

      if ( maxlevel < level ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R83VEC_SORT_QUICK_A - Fatal error!'
        write ( *, '(a,i6)' ) '  Exceeding recursion maximum of ', maxlevel
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
subroutine r83vec_uniform ( n, alo, ahi, seed, a )

!*****************************************************************************80
!
!! R83VEC_UNIFORM returns a random R83VEC in a given range.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 April 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, real ( kind = 8 ) ALO(3), AHI(3), the minimum and maximum values.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number
!    generator.
!
!    Output, real ( kind = 8 ) A(3,N), the vector of randomly chosen values.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: ndim = 3

  real ( kind = 8 ) a(ndim,n)
  real ( kind = 8 ) ahi(ndim)
  real ( kind = 8 ) alo(ndim)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) r8_uniform
  integer ( kind = 4 ) seed

  do i = 1, ndim
    do j = 1, n
      a(i,j) = r8_uniform ( alo(i), ahi(i), seed )
    end do
  end do

  return
end
subroutine r8mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8MAT_PRINT prints an R8MAT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 May 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix.
!
!    Input, character ( len = * ) TITLE, a title to be printed.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = * ) title

  call r8mat_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine r8mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8MAT_PRINT_SOME prints some of an R8MAT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 November 2003
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
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i7,7x)') j
    end do

    write ( *, '(''  Col   '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) ' '

    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi

      do j2 = 1, inc

        j = j2lo - 1 + j2

        write ( ctemp(j2), '(g14.6)' ) a(i,j)

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j), j = 1, inc )

    end do

  end do

  return
end
subroutine r8mat_transpose_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8MAT_TRANSPOSE_PRINT prints an R8MAT, transposed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
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
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = * ) title

  call r8mat_transpose_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine r8mat_transpose_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8MAT_TRANSPOSE_PRINT_SOME prints some of an R8MAT, transposed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
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
      write ( ctemp(i2), '(i7,7x)') i
    end do

    write ( *, '(''  Row   '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Col'

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
subroutine r8vec_bin ( n, a, nbin, bin_min, bin_max, bin, bin_limit )

!*****************************************************************************80
!
!! R8VEC_BIN bins an R8VEC, returning the population of each bin.
!
!  Discussion:
!
!    The user specifies minimum and maximum bin values, BIN_MIN and
!    BIN_MAX, and the number of bins, NBIN.  This determines a
!    "bin width":
!
!      H = ( BIN_MAX - BIN_MIN ) / NBIN
!
!    so that bin I will count all entries X(J) such that
!
!      BIN_LIMIT(I-1) <= A(J) < BIN_LIMIT(I).
!
!    The array does NOT have to be sorted.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries of A.
!
!    Input, real ( kind = 8 ) A(N), an (unsorted) array to be binned.
!
!    Input, integer ( kind = 4 ) NBIN, the number of bins.  Two extra bins, #0 and
!    #NBIN+1, count extreme values.
!
!    Input, real ( kind = 8 ) BIN_MIN, BIN_MAX, define the range and size
!    of the bins.  BIN_MIN and BIN_MAX must be distinct.
!    Normally, BIN_MIN < BIN_MAX, and the documentation will assume
!    this, but proper results will be computed if BIN_MAX < BIN_MIN.
!
!    Output, integer ( kind = 4 ) BIN(0:NBIN+1).
!    BIN(0) counts entries of A less than BIN_MIN.
!    BIN(NBIN+1) counts entries greater than or equal to BIN_MAX.
!    For 1 <= I <= NBIN, BIN(I) counts the entries X(J) such that
!      BIN_LIMIT(I-1) <= A(J) < BIN_LIMIT(I).
!    where H is the bin spacing.
!
!    Output, real ( kind = 8 ) BIN_LIMIT(0:NBIN), the "limits" of the bins.
!    BIN(I) counts the number of entries X(J) such that
!      BIN_LIMIT(I-1) <= A(J) < BIN_LIMIT(I).
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nbin

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) bin(0:nbin+1)
  real ( kind = 8 ) bin_limit(0:nbin)
  real ( kind = 8 ) bin_max
  real ( kind = 8 ) bin_min
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) t

  if ( bin_max == bin_min ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_BIN - Fatal error!'
    write ( *, '(a)' ) '  BIN_MIN = BIN_MAX.'
    stop
  end if

  bin(0:nbin+1) = 0

  do i = 1, n

    t = ( a(i) - bin_min ) / ( bin_max - bin_min )

    if ( t < 0.0D+00 ) then
      j = 0
    else if ( 1.0D+00 <= t ) then
      j = nbin + 1
    else
      j = 1 + int ( real ( nbin, kind = 8 ) * t )
    end if

    bin(j) = bin(j) + 1

  end do
!
!  Compute the bin limits.
!
  do i = 0, nbin
    bin_limit(i) = ( real ( nbin - i, kind = 8 ) * bin_min   &
                   + real (        i, kind = 8 ) * bin_max ) &
                   / real ( nbin,     kind = 8 )
  end do

  return
end
subroutine r8vec_bin_even ( n, a, nbin, bin_min, bin_max, bin_start, &
  bin_last, bin_next )

!*****************************************************************************80
!
!! R8VEC_BIN_EVEN bins an R8VEC into evenly spaced bins.
!
!  Discussion:
!
!    This is only a partial, indexed, sorting of the data.  To sort
!    the data, it is necessary to build a new array by extracting the
!    data for each bin, sorting that, and appending it to the array.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Thomas Cormen, Charles Leiserson, Ronald Rivest,
!    Introduction to Algorithms,
!    MIT Press, page 180.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of points in the set.
!
!    Input, real ( kind = 8 ) A(N), the data to be binned.
!
!    Input, integer ( kind = 4 ) NBIN, the number of bins.
!
!    Input, real ( kind = 8 ) BIN_MIN, BIN_MAX, the bin limits.
!
!    Output, integer ( kind = 4 ) BIN_START(NBIN), BIN_LAST(NBIN), the index of the
!    first and last element of A that went into each bin, or -1 if there
!    are no entries in this bin.
!
!    Output, integer ( kind = 4 ) BIN_NEXT(N), the index of the next element of A
!    that follows this element in the same bin.  A value of 0 means this
!    is the last entry in the particular bin.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nbin

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) bin_last(nbin)
  integer ( kind = 4 ) bin_next(n)
  integer ( kind = 4 ) bin_start(nbin)
  real ( kind = 8 ) bin_max
  real ( kind = 8 ) bin_min
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
!
!  Initialize the bin pointers to -1.
!
  bin_last(1:nbin) = -1
  bin_start(1:nbin) = -1
  bin_next(1:n) = -1
!
!  Build up linked lists of entries that go into each bin.
!
  do i = 1, n

    call r8_to_bin_even ( nbin, bin_min, bin_max, a(i), j )

    if ( bin_start(j) == -1 ) then
      bin_start(j) = i
    else
      k = bin_last(j)
      bin_next(k) = i
    end if

    bin_next(i) = 0

    bin_last(j) = i

  end do

  return
end
subroutine r8vec_binned_reorder ( n, a, nbin, bin_start, bin_last, bin_next )

!*****************************************************************************80
!
!! R8VEC_BINNED_REORDER reorders a binned R8VEC.
!
!  Discussion:
!
!    The bin vectors define an implicit ordering of the data
!    vector.  This routine physically rearranges the data vector
!    so that items in the first bin come first, and so on.  The
!    data within a bin is not reordered.
!
!    On output, the BIN_START, BIN_LAST and BIN_NEXT arrays have also been
!    updated so that they still correspond to the (rearranged) vector A.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data points.
!
!    Input/output, real ( kind = 8 ) A(N), the data to be sorted.
!
!    Input, integer ( kind = 4 ) NBIN, the number of bins.
!
!    Input/output, integer ( kind = 4 ) BIN_START(NBIN), BIN_LAST(NBIN), the index of
!    the first and last element of A that went into each bin, or -1 if
!    there are no entries in the bin.
!
!    Input/output, integer ( kind = 4 ) BIN_NEXT(N), the index of the next element of A
!    that follows this element in the same bin.  A value of 0 means this
!    is the last entry in the particular bin.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nbin

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) a2(n)
  integer ( kind = 4 ) bin_last(nbin)
  integer ( kind = 4 ) bin_next(n)
  integer ( kind = 4 ) bin_start(nbin)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k

  k = 0

  do i = 1, nbin

    j = bin_start(i)

    if ( 0 < j ) then
      bin_start(i) = k + 1
    end if

    do while ( 0 < j )
      k = k + 1
      bin_last(i) = k
      a2(k) = a(j)
      j = bin_next(j)
    end do

  end do

  a(1:n) = a2(1:n)
!
!  Now update the BIN_NEXT array.
!
  do i = 1, n
    bin_next(i) = i+1
  end do

  do i = 1, nbin
    k = bin_last(i)
    if ( 0 < k ) then
      bin_next(k) = 0
    end if
  end do

  return
end
subroutine r8vec_binned_sort_a ( n, a, nbin, bin_start, bin_last )

!*****************************************************************************80
!
!! R8VEC_BINNED_SORT_A ascending sorts a binned reordered R8VEC.
!
!  Discussion:
!
!    Presumably, the data vector was first binned by R8VEC_BIN_EVEN,
!    then reordered by R8VEC_BINNED_REORDER.  Now, each of the
!    bins of data is sorted one at a time, which results in sorting
!    the entire vector.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data points.
!
!    Input/output, real ( kind = 8 ) A(N), the data to be sorted.
!
!    Input, integer ( kind = 4 ) NBIN, the number of bins.
!
!    Input, integer ( kind = 4 ) BIN_START(NBIN), BIN_LAST(NBIN), the index of the first
!    and last element of A that went into each bin, or -1 if there were no
!    elements in the bin.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nbin

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) bin_last(nbin)
  integer ( kind = 4 ) bin_start(nbin)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) n1

  if ( n <= 1 ) then
    return
  end if

  do i = 1, nbin

    i1 = bin_start(i)

    if ( 0 < i1 ) then

      i2 = bin_last(i)

      n1 = i2 + 1 - i1

      call r8vec_sort_quick_a ( n1, a(i1:i2) )

    end if

  end do

  return
end
subroutine r8vec_bracket ( n, x, xval, left, right )

!*****************************************************************************80
!
!! R8VEC_BRACKET searches a sorted array for successive brackets of a value.
!
!  Discussion:
!
!    If the values in the vector are thought of as defining intervals
!    on the real line, then this routine searches for the interval
!    nearest to or containing the given value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, length of input array.
!
!    Input, real ( kind = 8 ) X(N), an array that has been sorted into
!    ascending order.
!
!    Input, real ( kind = 8 ) XVAL, a value to be bracketed.
!
!    Output, integer ( kind = 4 ) LEFT, RIGHT, the results of the search.
!    Either:
!      XVAL < X(1), when LEFT = 1, RIGHT = 2;
!      X(N) < XVAL, when LEFT = N-1, RIGHT = N;
!    or
!      X(LEFT) <= XVAL <= X(RIGHT).
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) left
  integer ( kind = 4 ) right
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xval

  do i = 2, n - 1

    if ( xval < x(i) ) then
      left = i - 1
      right = i
      return
    end if

   end do

  left = n - 1
  right = n

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
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vectors.
!
!    Input, real ( kind = 8 ) A1(N), A2(N), two vectors to compare.
!
!    Output, logical R8VEC_EQ.
!    R8VEC_EQ is TRUE if every pair of elements A1(I) and A2(I) are equal,
!    and FALSE otherwise.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a1(n)
  real ( kind = 8 ) a2(n)
  logical r8vec_eq

  r8vec_eq = ( all ( a1(1:n) == a2(1:n) ) )

  return
end
function r8vec_gt ( n, a1, a2 )

!*****************************************************************************80
!
!! R8VEC_GT == ( A1 > A2 ) for real vectors.
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
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the vectors.
!
!    Input, real ( kind = 8 ) A1(N), A2(N), the vectors to be compared.
!
!    Output, logical R8VEC_GT, is TRUE if and only if A1 > A2.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a1(n)
  real ( kind = 8 ) a2(n)
  integer ( kind = 4 ) i
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
!! R8VEC_LT == ( A1 < A2 ) for real vectors.
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
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the vectors.
!
!    Input, real ( kind = 8 ) A1(N), A2(N), the vectors to be compared.
!
!    Output, logical R8VEC_LT, is TRUE if and only if A1 < A2.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a1(n)
  real ( kind = 8 ) a2(n)
  integer ( kind = 4 ) i
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
subroutine r8vec_part_quick_a ( n, a, l, r )

!*****************************************************************************80
!
!! R8VEC_PART_QUICK_A reorders an R8VEC as part of a quick sort.
!
!  Discussion:
!
!    The routine reorders the entries of A.  Using A(1) as the key,
!    all entries of A that are less than or equal to the key will
!    precede the key which precedes all entries that are greater than the key.
!
!  Example:
!
!    Input:
!
!      N = 8
!
!      A = ( 6, 7, 3, 1, 6, 8, 2, 9 )
!
!    Output:
!
!      L = 3, R = 6
!
!      A = ( 3, 1, 2, 6, 6, 8, 9, 7 )
!            -------        -------
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries of A.
!
!    Input/output, real ( kind = 8 ) A(N).  On input, the array to be checked.
!    On output, A has been reordered as described above.
!
!    Output, integer ( kind = 4 ) L, R, the indices of A that define the 
!    three segments.  Let KEY = the input value of A(1).  Then
!    I <= L                 A(I) < KEY;
!         L < I < R         A(I) = KEY;
!                 R <= I    KEY < A(I).
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) key
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) r
  real ( kind = 8 ) temp

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_PART_QUICK_A - Fatal error!'
    write ( *, '(a)' ) '  N < 1.'
    stop
  else if ( n == 1 ) then
    l = 0
    r = 2
    return
  end if

  key = a(1)
  m = 1
!
!  The elements of unknown size have indices between L+1 and R-1.
!
  l = 1
  r = n + 1

  do i = 2, n

    if ( key < a(l+1) ) then
      r = r - 1
      temp = a(r)
      a(r) = a(l+1)
      a(l+1) = temp
    else if ( a(l+1) == key ) then
      m = m + 1
      temp = a(m)
      a(m) = a(l+1)
      a(l+1) = temp
      l = l + 1
    else if ( a(l+1) < key ) then
      l = l + 1
    end if

  end do
!
!  Now shift small elements to the left, and KEY elements to center.
!
  do i = 1, l - m
    a(i) = a(i+m)
  end do
!
!  Out of bounds here, occasionally
!
  l = l - m

  a(l+1:l+m) = key

  return
end
subroutine r8vec_print ( n, a, title )

!*****************************************************************************80
!
!! R8VEC_PRINT prints an R8VEC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, real ( kind = 8 ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i6,g16.8)' ) i, a(i)
  end do

  return
end
subroutine r8vec_sort_quick_a ( n, a )

!*****************************************************************************80
!
!! R8VEC_SORT_QUICK_A ascending sorts an R8VEC using quick sort.
!
!  Example:
!
!    Input:
!
!      N = 7
!      A = ( 6, 7, 3, 2, 9, 1, 8 )
!
!    Output:
!
!      A = ( 1, 2, 3, 6, 7, 8, 9 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input/output, real ( kind = 8 ) A(N).
!    On input, the array to be sorted.
!    On output, the array has been sorted.
!
  implicit none

  integer ( kind = 4 ), parameter :: level_max = 25
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) base
  integer ( kind = 4 ) l_segment
  integer ( kind = 4 ) level
  integer ( kind = 4 ) n_segment
  integer ( kind = 4 ) rsave(level_max)
  integer ( kind = 4 ) r_segment

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_SORT_QUICK_A - Fatal error!'
    write ( *, '(a)' ) '  N < 1.'
    stop
  else if ( n == 1 ) then
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
    call r8vec_part_quick_a ( n_segment, a(base), l_segment, r_segment )
!
!  If the left segment has more than one element, we need to partition it.
!
    if ( 1 < l_segment ) then

      if ( level_max < level ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8VEC_SORT_QUICK_A - Fatal error!'
        write ( *, '(a,i6)' ) '  Exceeding recursion maximum of ', level_max
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
subroutine r8vec_swap ( n, a1, a2 )

!*****************************************************************************80
!
!! R8VEC_SWAP swaps the entries of two real vectors.
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
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the arrays.
!
!    Input/output, real ( kind = 8 ) A1(N), A2(N), the vectors to swap.
!
  implicit none

  integer ( kind = 4 ) n

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
!! R8VEC_TO_BIN_EVEN3 determines the appropriate "bin" for an R8VEC.
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
!    08 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NDIM, the dimension of the space.
!
!    Input, integer ( kind = 4 ) NBIN(NDIM), the number of bins in each dimension.
!
!    Input, real ( kind = 8 ) A(NDIM), B(NDIM), the lower and upper limits
!    of the bin interval.  While A(I) is expected to be less than B(I),
!    the code should return useful results if A(I) is actually greater
!    than B(I).
!
!    Input, real ( kind = 8 ) C(NDIM), a value to be placed in a bin.
!
!    Output, integer ( kind = 4 ) BIN(NDIM), the index of the bin to which 
!    C is assigned.
!
  implicit none

  integer ( kind = 4 ) ndim

  real ( kind = 8 ) a(ndim)
  real ( kind = 8 ) b(ndim)
  integer ( kind = 4 ) bin(ndim)
  real ( kind = 8 ) c(ndim)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) nbin(ndim)

  do i = 1, ndim
    call r8_to_bin_even2 ( nbin(i), a(i), b(i), c(i), bin(i) )
  end do

  return
end
subroutine r8vec_uniform ( n, alo, ahi, seed, a )

!*****************************************************************************80
!
!! R8VEC_UNIFORM returns a random R8VEC in a given range.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 April 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, real ( kind = 8 ) ALO, AHI, the range allowed for the entries.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = 8 ) A(N), the vector of randomly chosen values.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) ahi
  real ( kind = 8 ) alo
  integer ( kind = 4 ) i
  real ( kind = 8 ) r8_uniform
  integer ( kind = 4 ) seed

  do i = 1, n
    a(i) = r8_uniform ( alo, ahi, seed )
  end do

  return
end
subroutine sort_heap_external ( n, indx, i, j, isgn )

!*****************************************************************************80
!
!! SORT_HEAP_EXTERNAL externally sorts a list of items into ascending order.
!
!  Discussion:
!
!    The actual list of data is not passed to the routine.  Hence this
!    routine may be used to sort integers, real ( kind = 8 )s, numbers, names,
!    dates, shoe sizes, and so on.  After each call, the routine asks
!    the user to compare or interchange two items, until a special
!    return value signals that the sorting is completed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 February 2004
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items to be sorted.
!
!    Input/output, integer ( kind = 4 ) INDX, the main communication signal.
!
!    The user must set INDX to 0 before the first call.
!    Thereafter, the user should not change the value of INDX until
!    the sorting is done.
!
!    On return, if INDX is
!
!      greater than 0,
!      * interchange items I and J;
!      * call again.
!
!      less than 0,
!      * compare items I and J;
!      * set ISGN = -1 if I < J, ISGN = +1 if J < I;
!      * call again.
!
!      equal to 0, the sorting is done.
!
!    Output, integer ( kind = 4 ) I, J, the indices of two items.
!    On return with INDX positive, elements I and J should be interchanged.
!    On return with INDX negative, elements I and J should be compared, and
!    the result reported in ISGN on the next call.
!
!    Input, integer ( kind = 4 ) ISGN, results of comparison of elements
!    I and J.  (Used only when the previous call returned INDX less than 0).
!    ISGN <= 0 means I is less than or equal to J;
!    0 <= ISGN means I is greater than or equal to J.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ), save :: i_save = 0
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  integer ( kind = 4 ), save :: j_save = 0
  integer ( kind = 4 ), save :: k = 0
  integer ( kind = 4 ), save :: k1 = 0
  integer ( kind = 4 ) n
  integer ( kind = 4 ), save :: n1 = 0
!
!  INDX = 0: This is the first call.
!
  if ( indx == 0 ) then

    i_save = 0
    j_save = 0
    k = n / 2
    k1 = k
    n1 = n
!
!  INDX < 0: The user is returning the results of a comparison.
!
  else if ( indx < 0 ) then

    if ( indx == -2 ) then

      if ( isgn < 0 ) then
        i_save = i_save + 1
      end if

      j_save = k1
      k1 = i_save
      indx = -1
      i = i_save
      j = j_save
      return

    end if

    if ( 0 < isgn ) then
      indx = 2
      i = i_save
      j = j_save
      return
    end if

    if ( k <= 1 ) then

      if ( n1 == 1 ) then
        i_save = 0
        j_save = 0
        indx = 0
      else
        i_save = n1
        n1 = n1 - 1
        j_save = 1
        indx = 1
      end if

      i = i_save
      j = j_save
      return

    end if

    k = k - 1
    k1 = k
!
!  0 < INDX, the user was asked to make an interchange.
!
  else if ( indx == 1 ) then

    k1 = k

  end if

  do

    i_save = 2 * k1

    if ( i_save == n1 ) then
      j_save = k1
      k1 = i_save
      indx = -1
      i = i_save
      j = j_save
      return
    else if ( i_save <= n1 ) then
      j_save = i_save + 1
      indx = -2
      i = i_save
      j = j_save
      return
    end if

    if ( k <= 1 ) then
      exit
    end if

    k = k - 1
    k1 = k

  end do

  if ( n1 == 1 ) then
    i_save = 0
    j_save = 0
    indx = 0
    i = i_save
    j = j_save
  else
    i_save = n1
    n1 = n1 - 1
    j_save = 1
    indx = 1
    i = i_save
    j = j_save
  end if

  return
end
subroutine swapec ( i, top, btri, bedg, point_num, point_xy, tri_num, &
  tri_vert, tri_nabe, stack, ierr )

!*****************************************************************************80
!
!! SWAPEC swaps diagonal edges until all triangles are Delaunay.
!
!  Discussion:
!
!    The routine swaps diagonal edges in a 2D triangulation, based on
!    the empty circumcircle criterion, until all triangles are Delaunay,
!    given that I is the index of the new vertex added to the triangulation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 July 2001
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the index of the new vertex.
!
!    Input/output, integer ( kind = 4 ) TOP, the index of the top of the stack.
!    On output, TOP is zero.
!
!    Input/output, integer ( kind = 4 ) BTRI, BEDG; on input, if positive,
!    are the triangle and edge indices of a boundary edge whose updated indices
!    must be recorded.  On output, these may be updated because of swaps.
!
!    Input, intger POINT_NUM, the number of points.
!
!    Input, real ( kind = 8 ) POINT_XY(2,POINT_NUM), the coordinates of
!    the points.
!
!    Input, integer ( kind = 4 ) TRI_NUM, the number of triangles.
!
!    Input/output, integer ( kind = 4 ) TRI_VERT(3,TRI_NUM), the triangle
!    incidence list.  May be updated on output because of swaps.
!
!    Input/output, integer ( kind = 4 ) TRI_NABE(3,TRI_NUM), the triangle 
!    neighbor list; negative values are used for links of the counter-clockwise 
!    linked list of boundary edges;  May be updated on output because of swaps.
!
!      LINK = -(3*I + J-1) where I, J = triangle, edge index.
!
!    Workspace, integer STACK(MAXST); on input, entries 1 through TOP
!    contain the indices of initial triangles (involving vertex I)
!    put in stack; the edges opposite I should be in interior;  entries
!    TOP+1 through MAXST are used as a stack.
!
!    Output, integer ( kind = 4 ) IERR is set to 8 for abnormal return.
!
  implicit none

  integer ( kind = 4 ) point_num
  integer ( kind = 4 ) tri_num

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) bedg
  integer ( kind = 4 ) btri
  integer ( kind = 4 ) c
  integer ( kind = 4 ) diaedg
  integer ( kind = 4 ) e
  integer ( kind = 4 ) ee
  integer ( kind = 4 ) em1
  integer ( kind = 4 ) ep1
  integer ( kind = 4 ) f
  integer ( kind = 4 ) fm1
  integer ( kind = 4 ) fp1
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) l
  integer ( kind = 4 ) r
  integer ( kind = 4 ) s
  integer ( kind = 4 ) stack(point_num)
  integer ( kind = 4 ) swap
  integer ( kind = 4 ) t
  integer ( kind = 4 ) top
  integer ( kind = 4 ) tri_nabe(3,tri_num)
  integer ( kind = 4 ) tri_vert(3,tri_num)
  integer ( kind = 4 ) tt
  integer ( kind = 4 ) u
  real ( kind = 8 ) point_xy(2,point_num)
  real ( kind = 8 ) x
  real ( kind = 8 ) y
!
!  Determine whether triangles in stack are Delaunay, and swap
!  diagonal edge of convex quadrilateral if not.
!
  x = point_xy(1,i)
  y = point_xy(2,i)

  do

    if ( top <= 0 ) then
      exit
    end if

    t = stack(top)
    top = top - 1

    if ( tri_vert(1,t) == i ) then
      e = 2
      b = tri_vert(3,t)
    else if ( tri_vert(2,t) == i ) then
      e = 3
      b = tri_vert(1,t)
    else
      e = 1
      b = tri_vert(2,t)
    end if

    a = tri_vert(e,t)
    u = tri_nabe(e,t)

    if ( tri_nabe(1,u) == t ) then
      f = 1
      c = tri_vert(3,u)
    else if ( tri_nabe(2,u) == t ) then
      f = 2
      c = tri_vert(1,u)
    else
      f = 3
      c = tri_vert(2,u)
    end if

    swap = diaedg ( x, y, point_xy(1,a), point_xy(2,a), point_xy(1,c), &
      point_xy(2,c), point_xy(1,b), point_xy(2,b) )

    if ( swap == 1 ) then

      em1 = i4_wrap ( e - 1, 1, 3 )
      ep1 = i4_wrap ( e + 1, 1, 3 )
      fm1 = i4_wrap ( f - 1, 1, 3 )
      fp1 = i4_wrap ( f + 1, 1, 3 )

      tri_vert(ep1,t) = c
      tri_vert(fp1,u) = i
      r = tri_nabe(ep1,t)
      s = tri_nabe(fp1,u)
      tri_nabe(ep1,t) = u
      tri_nabe(fp1,u) = t
      tri_nabe(e,t) = s
      tri_nabe(f,u) = r

      if ( 0 < tri_nabe(fm1,u) ) then
        top = top + 1
        stack(top) = u
      end if

      if ( 0 < s ) then

        if ( tri_nabe(1,s) == u ) then
          tri_nabe(1,s) = t
        else if ( tri_nabe(2,s) == u ) then
          tri_nabe(2,s) = t
        else
          tri_nabe(3,s) = t
        end if

        top = top + 1

        if ( point_num < top ) then
          ierr = 8
          return
        end if

        stack(top) = t

      else

        if ( u == btri .and. fp1 == bedg ) then
          btri = t
          bedg = e
        end if

        l = - ( 3 * t + e - 1 )
        tt = t
        ee = em1

        do while ( 0 < tri_nabe(ee,tt) )

          tt = tri_nabe(ee,tt)

          if ( tri_vert(1,tt) == a ) then
            ee = 3
          else if ( tri_vert(2,tt) == a ) then
            ee = 1
          else
            ee = 2
          end if

        end do

        tri_nabe(ee,tt) = l

      end if

      if ( 0 < r ) then

        if ( tri_nabe(1,r) == t ) then
          tri_nabe(1,r) = u
        else if ( tri_nabe(2,r) == t ) then
          tri_nabe(2,r) = u
        else
          tri_nabe(3,r) = u
        end if

      else

        if ( t == btri .and. ep1 == bedg ) then
          btri = u
          bedg = f
        end if

        l = - ( 3 * u + f - 1 )
        tt = u
        ee = fm1

        do while ( 0 < tri_nabe(ee,tt) )

          tt = tri_nabe(ee,tt)

          if ( tri_vert(1,tt) == b ) then
            ee = 3
          else if ( tri_vert(2,tt) == b ) then
            ee = 1
          else
            ee = 2
          end if

        end do

        tri_nabe(ee,tt) = l

      end if

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
subroutine triangle_area_2d ( t, area )

!*****************************************************************************80
!
!! TRIANGLE_AREA_2D computes the area of a triangle in 2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T(2,3), the triangle vertices.
!
!    Output, real ( kind = 8 ) AREA, the absolute area of the triangle.
!
  implicit none

  real ( kind = 8 ) area
  real ( kind = 8 ) t(2,3)

  area = 0.5D+00 * abs ( &
      t(1,1) * ( t(2,2) - t(2,3) ) &
    + t(1,2) * ( t(2,3) - t(2,1) ) &
    + t(1,3) * ( t(2,1) - t(2,2) ) )

  return
end
subroutine triangle_sample_2d ( t, seed, p )

!*****************************************************************************80
!
!! TRIANGLE_SAMPLE_2D returns a random point in a triangle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T(2,3), the triangle vertices.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = 8 ) P(2), a random point in the triangle.
!
  implicit none

  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  real ( kind = 8 ) p(2)
  real ( kind = 8 ) r
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) t(2,3)
  real ( kind = 8 ) x12
  real ( kind = 8 ) x13
  real ( kind = 8 ) y12
  real ( kind = 8 ) y13

  r = r8_uniform_01 ( seed )
!
!  Interpret R as a percentage of the triangle's area.
!
!  Imagine a line L, parallel to side 1, so that the area between
!  vertex 1 and line L is R percent of the full triangle's area.
!
!  The line L will intersect sides 2 and 3 at a fraction
!  ALPHA = SQRT ( R ) of the distance from vertex 1 to vertices 2 and 3.
!
  alpha = sqrt ( r )
!
!  Determine the coordinates of the points on sides 2 and 3 intersected
!  by line L.
!
  x12 = alpha * t(1,1) + ( 1.0D+00 - alpha ) * t(1,2)
  y12 = alpha * t(2,1) + ( 1.0D+00 - alpha ) * t(2,2)

  x13 = alpha * t(1,1) + ( 1.0D+00 - alpha ) * t(1,3)
  y13 = alpha * t(2,1) + ( 1.0D+00 - alpha ) * t(2,3)
!
!  Now choose, uniformly at random, a point on the line L.
!
  beta = r8_uniform_01 ( seed )

  p(1) = beta * x12 + ( 1.0D+00 - beta ) * x13
  p(2) = beta * y12 + ( 1.0D+00 - beta ) * y13

  return
end
subroutine triangulation_nabe_nodes ( point_num, tri_num, tri_vert, &
  nabes_first, nabes_num, nabes_max, nabes_dim, nabes )

!*****************************************************************************80
!
!! TRIANGULATION_NABE_NODES determines the neighbors of triangulation nodes.
!
!  Example:
!
!    On input, the triangle data structure is:
!
!    Triangle  Nodes
!    --------  ----------
!     1        3,   4,   1
!     2        3,   1,   2
!     3        3,   2,   6
!     4        2,   1,   5
!     5        6,   2,   5
!
!  On output, the auxilliary neighbor arrays are:
!
!    Node  Num  First
!    ----  ---  -----
!     1     4     1
!     2     4     5
!     3     4     9
!     4     2    13
!     5     3    15
!     6     3    18
!
!  and the neighbor array is:
!
!    Position  Node
!    --------  ----
!
!     1        2
!     2        3
!     3        4
!     4        5
!    -----------
!     5        1
!     6        3
!     7        5
!     8        6
!    -----------
!     9        1
!    10        2
!    11        4
!    12        6
!    -----------
!    13        1
!    14        3
!    -----------
!    15        1
!    16        2
!    17        6
!    -----------
!    18        2
!    19        3
!    20        5
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 July 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, integer ( kind = 4 ) TRI_NUM, the number of triangles.
!
!    Input, integer ( kind = 4 ) TRI_VERT(3,TRI_NUM), the nodes that make
!    up each triangle.
!
!    Output, integer ( kind = 4 ) NABES_FIRST(POINT_NUM), the index in NABES
!    of the first neighbor in the list for each node.
!
!    Output, integer ( kind = 4 ) NABES_NUM(POINT_NUM), the number of neighbors
!    of each node.
!
!    Input, integer ( kind = 4 ) NABES_MAX, the maximum dimension of NABES.
!
!    Output, integer ( kind = 4 ) NABES_DIM, the dimension of NABES.
!
!    Output, integer ( kind = 4 ) NABES(NABES_DIM), a list of the neighbors
!    of all the nodes.  Neighbors of node 1 are listed first, and so on.
!
  implicit none

  integer ( kind = 4 ) nabes_max
  integer ( kind = 4 ) point_num
  integer ( kind = 4 ) tri_num

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_current
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) nabe
  integer ( kind = 4 ) nabes(nabes_max)
  integer ( kind = 4 ) nabes1(nabes_max)
  integer ( kind = 4 ) nabes_dim
  integer ( kind = 4 ) nabes_first(point_num)
  integer ( kind = 4 ) nabes_num(point_num)
  integer ( kind = 4 ) nuniq
  integer ( kind = 4 ) tri
  integer ( kind = 4 ) tri_vert(3,tri_num)
!
!  Step 1.  From the triangle list (I,J,K)
!  construct the neighbor relations: (I,J), (J,K), (K,I), (J,I), (K,J), (I,K).
!
  nabes_dim = 0
  do tri = 1, tri_num
    i = tri_vert(1,tri)
    j = tri_vert(2,tri)
    k = tri_vert(3,tri)
    nabes1(nabes_dim+1:nabes_dim+6) = (/ i, i, j, j, k, k /)
    nabes(nabes_dim+1:nabes_dim+6) = (/ j, k, i, k, i, j /)
    nabes_dim = nabes_dim + 6
  end do
!
!  Step 2. Dictionary sort the neighbor relations.
!
  call i4vec2_sort_a ( nabes_dim, nabes1, nabes )
!
!  Step 3. Remove duplicate entries.
!
  call i4vec2_sorted_unique ( nabes_dim, nabes1, nabes, nuniq )

  nabes_dim = nuniq
!
!  Step 4. Construct the NABES_NUM and NABES_FIRST data.
!
  nabes_num(1:point_num) = 0
  nabes_first(1:point_num) = 0
  i_current = 0
  do nabe = 1, nabes_dim
    i = nabes1(nabe)
    if ( i == i_current ) then
      nabes_num(i) = nabes_num(i) + 1
    else
      i_current = i
      nabes_first(i) = nabe
      nabes_num(i) = 1
    end if
  end do

  return
end
subroutine triangulation_print ( point_num, tri_num, xc, tri_vert, tri_nabe )

!*****************************************************************************80
!
!! TRIANGULATION_PRINT prints out information defining a Delaunay triangulation.
!
!  Discussion:
!
!    Triangulations created by RTRIS include extra information encoded
!    in the negative values of TRI_NABE.
!
!    Because some of the nodes counted in POINT_NUM may not actually be
!    used in the triangulation, I needed to compute the true number
!    of vertices.  I added this calculation on 13 October 2001.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 November 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, integer ( kind = 4 ) TRI_NUM, the number of triangles.
!
!    Input, real ( kind = 8 ) XC(2,POINT_NUM), the point coordinates.
!
!    Input, integer ( kind = 4 ) TRI_VERT(3,TRI_NUM), the nodes that make
!    up the triangles.
!
!    Input, integer ( kind = 4 ) TRI_NABE(3,TRI_NUM), the triangle neighbors 
!    on each side.  If there is no triangle neighbor on a particular side, the 
!    value of TRI_NABE should be negative.  If the triangulation data was 
!    created by DTRIS2, then there is more information encoded in the negative
!    values.
!
  implicit none

  integer ( kind = 4 ) point_num
  integer ( kind = 4 ) tri_num

  integer ( kind = 4 ) boundary_num
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) s
  logical skip
  integer ( kind = 4 ) t
  integer ( kind = 4 ) tri_nabe(3,tri_num)
  integer ( kind = 4 ) tri_vert(3,tri_num)
  integer ( kind = 4 ), allocatable, dimension ( : ) :: vertex_list
  integer ( kind = 4 ) vertex_num
  real ( kind = 8 ) xc(2,point_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TRIANGULATION_PRINT'
  write ( *, '(a)' ) '  Information defining a triangulation.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  The number of points is ', point_num

  call r8mat_transpose_print ( 2, point_num, xc, '  Point coordinates' )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  The number of triangles is ', tri_num
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Sets of three points are used as vertices of'
  write ( *, '(a)' ) '  the triangles.  For each triangle, the points'
  write ( *, '(a)' ) '  are listed in counterclockwise order.'

  call i4mat_transpose_print ( 3, tri_num, tri_vert, '  Triangle nodes:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  On each side of a given triangle, there is either'
  write ( *, '(a)' ) '  another triangle, or a piece of the convex hull.'
  write ( *, '(a)' ) '  For each triangle, we list the indices of the three'
  write ( *, '(a)' ) '  neighbors, or (if negative) the codes of the'
  write ( *, '(a)' ) '  segments of the convex hull.'

  call i4mat_transpose_print ( 3, tri_num, tri_nabe, '  Triangle neighbors' )
!
!  Determine the number of vertices.  This is not the same as the
!  number of points!
!
  allocate ( vertex_list(1:3*tri_num) )

  vertex_list(1:3*tri_num) = reshape ( tri_vert(1:3,1:tri_num), &
    (/ 3*tri_num /) )

  call i4vec_sort_heap_a ( 3*tri_num, vertex_list )

  call i4vec_sorted_unique ( 3*tri_num, vertex_list, vertex_num )

  deallocate ( vertex_list )
!
!  Determine the number of boundary points.
!
  boundary_num = 2 * vertex_num - tri_num - 2

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  The number of boundary points is ', boundary_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The segments that make up the convex hull can be'
  write ( *, '(a)' ) '  determined from the negative entries of the triangle'
  write ( *, '(a)' ) '  neighbor list.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  # Tri Side  N1  N2'
  write ( *, '(a)' ) ' '

  skip = .false.

  k = 0

  do i = 1, tri_num

    do j = 1, 3

      if ( tri_nabe(j,i) < 0 ) then
        s = - tri_nabe(j,i)
        t = s / 3

        if ( t < 1 .or. tri_num < t ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  Sorry, this data does not use the DTRIS2'
          write ( *, '(a)' ) '  convention for convex hull segments.'
          skip = .true.
          exit
        end if

        s = mod ( s, 3 ) + 1
        k = k + 1
        n1 = tri_vert(s,t)
        n2 = tri_vert(i4_wrap(s+1,1,3),t)
        write ( *, '(5i4)' ) k, t, s, n1, n2
      end if

    end do

    if ( skip ) then
      exit
    end if

  end do

  return
end
subroutine triangulation_sample ( point_num, xc, tri_num, tri_vert, &
  num_ran, seed, xd, td )

!*****************************************************************************80
!
!! TRIANGULATION_SAMPLE returns random points in a triangulation.
!
!  Discussion:
!
!    It is assumed that the triangulation consists of a set of non-overlapping
!    triangles.
!
!    The point is chosen uniformly in the area covered by the triangulation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 July 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points used in
!    the triangulation.
!
!    Input, real ( kind = 8 ) XC(2,POINT_NUM), the point coordinates.
!
!    Input, integer ( kind = 4 ) TRI_NUM, the number of triangles.
!
!    Input, integer ( kind = 4 ) TRI_VERT(3,TRI_NUM), the nodes that make
!    up the triangles.
!
!    Input, integer ( kind = 4 ) NUM_RAN, the number of points to sample.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = 8 ) XD(2,NUM_RAN), the sample points.
!
!    Output, integer ( kind = 4 ) TD(NUM_RAN), the triangle to which each
!    sample point belongs.
!
  implicit none

  integer ( kind = 4 ) point_num
  integer ( kind = 4 ) num_ran
  integer ( kind = 4 ) tri_num

  real ( kind = 8 ) area
  real ( kind = 8 ) area_cum(0:tri_num)
  real ( kind = 8 ) area_total
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) left
  real ( kind = 8 ) r
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) right
  integer ( kind = 4 ) seed
  real ( kind = 8 ) t(2,3)
  integer ( kind = 4 ) td(num_ran)
  integer ( kind = 4 ) tri_vert(3,tri_num)
  real ( kind = 8 ) xc(2,point_num)
  real ( kind = 8 ) xd(2,num_ran)
!
!  Compute the areas of the triangles.
!  Build a cumulative area vector.
!  Convert it to a relative cumulative area vector.
!
  area_cum(0) = 0.0D+00

  do i = 1, tri_num

    i1 = tri_vert(1,i)
    i2 = tri_vert(2,i)
    i3 = tri_vert(3,i)

    t(1:2,1) = xc(1:2,i1)
    t(1:2,2) = xc(1:2,i2)
    t(1:2,3) = xc(1:2,i3)

    call triangle_area_2d ( t, area )

    area_cum(i) = area_cum(i-1) + area

  end do

  area_total = area_cum(tri_num)

  area_cum(0:tri_num) = area_cum(0:tri_num) / area_total
!
!  Pick random values.  A random value R indicates the corresponding triangle
!  whose cumulative relative area contains R.
!
!  Bracket the random value in the cumulative relative areas,
!  indicating a triangle.
!
!  Pick a random point in the triangle.
!
  do i = 1, num_ran

    r = r8_uniform_01 ( seed )

    call r8vec_bracket ( tri_num+1, area_cum, r, left, right )

    td(i) = right - 1

    i1 = tri_vert(1,td(i))
    i2 = tri_vert(2,td(i))
    i3 = tri_vert(3,td(i))

    t(1:2,1) = xc(1:2,i1)
    t(1:2,2) = xc(1:2,i2)
    t(1:2,3) = xc(1:2,i3)

    call triangle_sample_2d ( t, seed, xd(1:2,i) )

  end do

  return
end
subroutine tuple_next2 ( n, xmin, xmax, x, rank )

!*****************************************************************************80
!
!! TUPLE_NEXT2 computes the next element of an integer tuple space.
!
!  Discussion:
!
!    The elements X are N vectors.
!
!    Each entry X(I) is constrained to lie between XMIN(I) and XMAX(I).
!
!    The elements are produced one at a time.
!
!    The first element is
!      (XMIN(1), XMIN(2), ..., XMIN(N)),
!    the second is (probably)
!      (XMIN(1), XMIN(2), ..., XMIN(N)+1),
!    and the last element is
!      (XMAX(1), XMAX(2), ..., XMAX(N))
!
!    Intermediate elements are produced in a lexicographic order, with
!    the first index more important than the last, and the ordering of
!    values at a fixed index implicitly defined by the sign of
!    XMAX(I) - XMIN(I).
!
!  Example:
!
!    N = 2,
!    XMIN = (/ 1, 10 /)
!    XMAX = (/ 3,  8 /)
!
!    RANK    X
!    ----  -----
!      1   1 10
!      2   1  9
!      3   1  8
!      4   2 10
!      5   2  9
!      6   2  8
!      7   3 10
!      8   3  9
!      9   3  8
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components.
!
!    Input, integer ( kind = 4 ) XMIN(N), XMAX(N), the "minimum" and "maximum"
!    entry values.  These values are minimum and maximum only in the sense
!    of the lexicographic ordering.  In fact, XMIN(I) may be less than,
!    equal to, or greater than XMAX(I).
!
!    Input/output, integer ( kind = 4 ) X(N), on input the previous tuple.
!    On output, the next tuple.
!
!    Input/output, integer ( kind = 4 ) RANK, the rank of the item.  On first call,
!    set RANK to 0 to start up the sequence.  On return, if RANK is zero,
!    there are no more items in the sequence.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) x(n)
  integer ( kind = 4 ) xmin(n)
  integer ( kind = 4 ) xmax(n)

  if ( rank < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TUPLE_NEXT2 - Fatal error!'
    write ( *, '(a,i6)' ) '  Illegal value of RANK = ', rank
    stop
  end if

  if ( product ( 1 + abs ( xmax(1:n) - xmin(1:n) ) ) < rank ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TUPLE_NEXT2 - Fatal error!'
    write ( *, '(a,i6)' ) '  Illegal value of RANK = ', rank
    stop
  end if

  if ( rank == 0 ) then
    x(1:n) = xmin(1:n)
    rank = 1
    return
  end if

  rank = rank + 1
  i = n

  do

    if ( x(i) /= xmax(i) ) then
      x(i) = x(i) + sign ( 1, xmax(i) - xmin(i) )
      exit
    end if

    x(i) = xmin(i)

    if ( i == 1 ) then
      rank = 0
      exit
    end if

    i = i - 1

  end do

  return
end
subroutine vbedg ( x, y, point_num, point_xy, tri_num, tri_vert, tri_nabe, &
  ltri, ledg, rtri, redg )

!*****************************************************************************80
!
!! VBEDG determines which boundary edges are visible to a point.
!
!  Discussion:
!
!    The point (X,Y) is assumed to be outside the convex hull of the
!    region covered by the 2D triangulation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 July 2009
!
!  Author:
!
!    Original FORTRAN77 version by Barry Joe.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, Y, the coordinates of a point outside the
!    convex hull of the current triangulation.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, real ( kind = 8 ) POINT_XY(2,POINT_NUM), the coordinates of the
!    vertices.
!
!    Input, integer ( kind = 4 ) TRI_NUM, the number of triangles.
!
!    Input, integer ( kind = 4 ) TRI_VERT(3,TRI_NUM), the triangle incidence list.
!
!    Input, integer ( kind = 4 ) TRI_NABE(3,TRI_NUM), the triangle neighbor
!    list; negative values are used for links of a counter clockwise linked
!    list of boundary edges;
!      LINK = -(3*I + J-1) where I, J = triangle, edge index.
!
!    Input/output, integer ( kind = 4 ) LTRI, LEDG.  If LTRI /= 0 then these
!    values are assumed to be already computed and are not changed, else they
!    are updated.  On output, LTRI is the index of boundary triangle to the
!    left of the leftmost boundary triangle visible from (X,Y), and LEDG is
!    the boundary edge of triangle LTRI to the left of the leftmost boundary
!    edge visible from (X,Y).  1 <= LEDG <= 3.
!
!    Input/output, integer ( kind = 4 ) RTRI.  On input, the index of the 
!    boundary triangle to begin the search at.  On output, the index of the 
!    rightmost boundary triangle visible from (X,Y).
!
!    Input/output, integer ( kind = 4 ) REDG, the edge of triangle RTRI that
!    is visible from (X,Y).  1 <= REDG <= 3.
!
  implicit none

  integer ( kind = 4 ) point_num
  integer ( kind = 4 ) tri_num

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) e
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) l
  logical ldone
  integer ( kind = 4 ) ledg
  integer ( kind = 4 ) lr
  integer ( kind = 4 ) lrline
  integer ( kind = 4 ) ltri
  real ( kind = 8 ) point_xy(2,point_num)
  integer ( kind = 4 ) redg
  integer ( kind = 4 ) rtri
  integer ( kind = 4 ) t
  integer ( kind = 4 ) tri_nabe(3,tri_num)
  integer ( kind = 4 ) tri_vert(3,tri_num)
  real ( kind = 8 ) x
  real ( kind = 8 ) y
!
!  Find the rightmost visible boundary edge using links, then possibly
!  leftmost visible boundary edge using triangle neighbor information.
!
  if ( ltri == 0 ) then
    ldone = .false.
    ltri = rtri
    ledg = redg
  else
    ldone = .true.
  end if

  do

    l = -tri_nabe(redg,rtri)
    t = l / 3
    e = mod ( l, 3 ) + 1
    a = tri_vert(e,t)

    if ( e <= 2 ) then
      b = tri_vert(e+1,t)
    else
      b = tri_vert(1,t)
    end if

    lr = lrline ( x, y, point_xy(1,a), point_xy(2,a), point_xy(1,b), &
      point_xy(2,b), 0.0D+00 )

    if ( lr <= 0 ) then
      exit
    end if

    rtri = t
    redg = e

  end do

  if ( ldone ) then
    return
  end if

  t = ltri
  e = ledg

  do

    b = tri_vert(e,t)
    e = i4_wrap ( e-1, 1, 3 )

    do while ( 0 < tri_nabe(e,t) )

      t = tri_nabe(e,t)

      if ( tri_vert(1,t) == b ) then
        e = 3
      else if ( tri_vert(2,t) == b ) then
        e = 1
      else
        e = 2
      end if

    end do

    a = tri_vert(e,t)

    lr = lrline ( x, y, point_xy(1,a), point_xy(2,a), point_xy(1,b), &
       point_xy(2,b), 0.0D+00 )

    if ( lr <= 0 ) then
      exit
    end if

  end do

  ltri = t
  ledg = e

  return
end
