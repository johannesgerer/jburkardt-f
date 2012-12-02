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
!    counterclockwise order.
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

    diaedg = - 1

  else if ( ca < - tola .and. cb < - tolb ) then

    diaedg = 1

  else

    tola = max ( tola, tolb )

    s = ( dx10 * dy30 - dx30 * dy10 ) * cb &
      + ( dx32 * dy12 - dx12 * dy32 ) * ca

    if ( tola < s ) then
      diaedg = - 1
    else if ( s < - tola ) then
      diaedg = 1
    else
      diaedg = 0
    end if

  end if

  return
end
function i4_sign ( x )

!*****************************************************************************80
!
!! I4_SIGN evaluates the sign of an I4.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ) value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) X, the number whose sign is desired.
!
!    Output, integer ( kind = 4 ) I4_SIGN, the sign of X:
!
  implicit none

  integer ( kind = 4 ) i4_sign
  integer ( kind = 4 ) x

  if ( x < 0 ) then
    i4_sign = -1
  else
    i4_sign = +1
  end if

  return
end
subroutine i4_swap ( i, j )

!*****************************************************************************80
!
!! I4_SWAP switches two I4's.
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
subroutine i4mat_transpose_print ( m, n, a, title )

!*****************************************************************************80
!
!! I4MAT_TRANSPOSE_PRINT prints an I4MAT, transposed.
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
!! I4MAT_TRANSPOSE_PRINT_SOME prints some of the transpose of an I4MAT.
!
!  Modified:
!
!    09 February 2005
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
    j2hi = min ( jhi, n )

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
subroutine i4vec_min ( n, a, amin )

!*****************************************************************************80
!
!! I4VEC_MIN computes the minimum element of an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, integer ( kind = 4 ) A(N), the array.
!
!    Output, integer ( kind = 4 ) AMIN, the value of the smallest entry.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) amin

  amin = minval ( a(1:n) )

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
subroutine i4vec_sorted_unique ( n, a, unique_num )

!*****************************************************************************80
!
!! I4VEC_SORTED_UNIQUE gets the unique elements in a sorted I4VEC.
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
!    integer array.  On output, the unique elements in A.
!
!    Output, integer ( kind = 4 ) UNIQUE_NUM, the number of unique elements
!    in A.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) itest
  integer ( kind = 4 ) unique_num

  unique_num = 0

  if ( n <= 0 ) then
    return
  end if

  unique_num = 1

  do itest = 2, n

    if ( a(itest) /= a(unique_num) ) then
      unique_num = unique_num + 1
      a(unique_num) = a(itest)
    end if

  end do

  return
end
function lrline ( xu, yu, xv1, yv1, xv2, yv2, dv )

!*****************************************************************************80
!
!! LRLINE determines if a point is left of, right or, or on a directed line.
!
!  Discussion:
!
!    The directed line is parallel to, and at a signed distance DV from
!    a directed base line from (XV1,YV1) to (XV2,YV2).
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
!    Input, real ( kind = 8 ) XU, YU, the coordinates of the point whose
!    position relative to the directed line is to be determined.
!
!    Input, real ( kind = 8 ) XV1, YV1, XV2, YV2, the coordinates of two points
!    that determine the directed base line.
!
!    Input, real ( kind = 8 ) DV, the signed distance of the directed line
!    from the directed base line through the points (XV1,YV1) and (XV2,YV2).
!    DV is positive for a line to the left of the base line.
!
!    Output, integer ( kind = 4 ) LRLINE, the result:
!    +1, the point is to the right of the directed line;
!     0, the point is on the directed line;
!    -1, the point is to the left of the directed line.
!
  implicit none

  real ( kind = 8 ) dv
  real ( kind = 8 ) dx
  real ( kind = 8 ) dxu
  real ( kind = 8 ) dy
  real ( kind = 8 ) dyu
  integer ( kind = 4 ) lrline
  real ( kind = 8 ) t
  real ( kind = 8 ) tol
  real ( kind = 8 ) tolabs
  real ( kind = 8 ) xu
  real ( kind = 8 ) xv1
  real ( kind = 8 ) xv2
  real ( kind = 8 ) yu
  real ( kind = 8 ) yv1
  real ( kind = 8 ) yv2

  tol = 100.0D+00 * epsilon ( tol )

  dx = xv2 - xv1
  dy = yv2 - yv1
  dxu = xu - xv1
  dyu = yu - yv1

  tolabs = tol * max ( abs ( dx ), abs ( dy ), abs ( dxu ), &
    abs ( dyu ), abs ( dv ) )

  t = dy * dxu - dx * dyu + dv * sqrt ( dx * dx + dy * dy )

  if ( tolabs < t ) then
    lrline = 1
  else if ( -tolabs <= t ) then
    lrline = 0
  else
    lrline = -1
  end if

  return
end
subroutine perm_check2 ( n, p, base, ierror )

!*****************************************************************************80
!
!! PERM_CHECK2 checks that a vector represents a permutation.
!
!  Discussion:
!
!    The routine verifies that each of the integers from BASE to
!    to BASE+N-1 occurs among the N entries of the permutation.
!
!    Set the input quantity BASE to 0, if P is a 0-based permutation,
!    or to 1 if P is a 1-based permutation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries.
!
!    Input, integer ( kind = 4 ) P(N), the array to check.
!
!    Input, integer ( kind = 4 ) BASE, the index base.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, the array represents a permutation.
!    nonzero, the array does not represent a permutation.  The smallest
!    missing value is equal to IERROR.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) base
  integer ( kind = 4 ) find
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) p(n)
  integer ( kind = 4 ) seek

  ierror = 0

  do seek = base, base + n - 1

    ierror = 1

    do find = 1, n
      if ( p(find) == seek ) then
        ierror = 0
        exit
      end if
    end do

    if ( ierror /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PERM_CHECK2 - Fatal error!'
      write ( *, '(a)' ) '  The input array does not represent'
      write ( *, '(a)' ) '  a proper permutation.'
      stop
    end if

  end do

  return
end
subroutine perm_inverse ( n, p )

!*****************************************************************************80
!
!! PERM_INVERSE inverts a permutation "in place".
!
!  Discussion:
!
!    This algorithm assumes that the entries in the permutation vector are
!    strictly positive.  In particular, the value 0 must not occur.
!
!    When necessary, this function shifts the data temporarily so that
!    this requirement is satisfied.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 June 2009
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

  integer ( kind = 4 ) base
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i0
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i4_sign
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) is
  integer ( kind = 4 ) p(n)
  integer ( kind = 4 ) p_min

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PERM_INVERSE - Fatal error!'
    write ( *, '(a,i8)' ) '  Input value of N = ', n
    stop
  end if
!
!  Find the least value, and shift data so it begins at 1.
!
  call i4vec_min ( n, p, p_min )
  base = 1
  p(1:n) = p(1:n) - p_min + base
!
!  Check the permutation.
!
  call perm_check2 ( n, p, base, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PERM_INVERSE - Fatal error!'
    write ( *, '(a)' ) '  PERM_CHECK rejects this permutation.'
    stop
  end if
!
!  Invert the permutation.
!
  is = 1

  do i = 1, n

    i1 = p(i)

    do while ( i < i1 )
      i2 = p(i1)
      p(i1) = -i2
      i1 = i2
    end do

    is = - i4_sign ( p(i) )
    p(i) = is * abs ( p(i) )

  end do

  do i = 1, n

    i1 = - p(i)

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
!
!  Reverse the shift.
!
  p(1:n) = p(1:n) + p_min - base

  return
end
subroutine pwl_interp_2d_scattered_value ( nd, xyd, zd, t_num, t, t_neighbor, &
  ni, xyi, zi )

!*****************************************************************************80
!
!! PWL_INTERP_2D_SCATTERED_VALUE evaluates a 2d interpolant of scattered data
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ND, the number of data points.
!
!    Input, real ( kind = 8 ) XYD(2,ND), the data point coordinates.
!
!    Input, real ( kind = 8 ) ZD(ND), the data values.
!
!    Input, integer ( kind = 4 ) T_NUM, the number of triangles.
!
!    Input, integer ( kind = 4 ) T(3,T_NUM), the triangle information.
!
!    Input, integer ( kind = 4 ) T_NEIGHBOR(3,T_NUM), the triangle neighbors.
!
!    Input, integer ( kind = 4 ) NI, the number of interpolation points.
!
!    Input, real ( kind = 8 ) XYI(2,NI), the interpolation point coordinates.
!
!    Output, real ( kind = 8 ) ZI(NI), the interpolated values.
!
  implicit none

  integer ( kind = 4 ) nd
  integer ( kind = 4 ) ni
  integer ( kind = 4 ) t_num

  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  real ( kind = 8 ) gamma
  integer ( kind = 4 ) edge
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) step_num
  integer ( kind = 4 ) t(3,t_num)
  integer ( kind = 4 ) t_neighbor(3,t_num)
  real ( kind = 8 ) xyd(2,nd)
  real ( kind = 8 ) xyi(2,ni)
  real ( kind = 8 ) zd(nd)
  real ( kind = 8 ) zi(ni)

  do i = 1, ni

    call triangulation_search_delaunay ( nd, xyd, 3, t_num, t, t_neighbor, &
      xyi(1:2,i), j, alpha, beta, gamma, edge, step_num )

    if ( j == -1 ) then
      zi(i) = -1.0D+00
    end if

    zi(i) = alpha * zd(t(1,j)) &
          + beta  * zd(t(2,j)) &
          + gamma * zd(t(3,j))

  end do

  return
end
subroutine r8tris2 ( node_num, node_xy, element_num, element_node, &
  element_neighbor )

!*****************************************************************************80
!
!! R8TRIS2 constructs a Delaunay triangulation of 2D vertices.
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
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input/output, real ( kind = 8 ) NODE_XY(2,NODE_NUM), the coordinates
!    of the nodes.  On output, the vertices have been sorted into
!    dictionary order.
!
!    Output, integer ( kind = 4 ) ELEMENT_NUM, the number of triangles in the
!    triangulation;  ELEMENT_NUM is equal to 2*NODE_NUM - NB - 2, where NB is
!    the number of boundary vertices.
!
!    Output, integer ( kind = 4 ) ELEMENT_NODE(3,ELEMENT_NUM), the nodes that
!    make up each triangle.  The elements are indices of P.  The vertices of
!    the triangles are in counter clockwise order.
!
!    Output, integer ( kind = 4 ) ELEMENT_NEIGHBOR(3,ELEMENT_NUM), the
!    triangle neighbor list.  Positive elements are indices of TIL; negative
!    elements are used for links of a counter clockwise linked list of boundary
!    edges;  LINK = -(3*I + J-1) where I, J = triangle, edge index;
!    ELEMENT_NEIGHBOR(J,I) refers to the neighbor along edge from vertex J
!    to J+1 (mod 3).
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ) node_num

  real ( kind = 8 ) cmax
  integer ( kind = 4 ) e
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) indx(node_num)
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
  real ( kind = 8 ) node_xy(dim_num,node_num)
  integer ( kind = 4 ) redg
  integer ( kind = 4 ) rtri
  integer ( kind = 4 ) stack(node_num)
  integer ( kind = 4 ) t
  real ( kind = 8 ) tol
  integer ( kind = 4 ) top
  integer ( kind = 4 ) element_neighbor(3,node_num*2)
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ) element_node(3,node_num*2)

  tol = 100.0D+00 * epsilon ( tol )

  ierr = 0
!
!  Sort the vertices by increasing (x,y).
!
  call r82vec_sort_heap_index_a ( node_num, node_xy, indx )

  call r82vec_permute ( node_num, indx, node_xy )
!
!  Make sure that the data nodes are "reasonably" distinct.
!
  m1 = 1

  do i = 2, node_num

    m = m1
    m1 = i

    k = 0

    do j = 1, dim_num

      cmax = max ( abs ( node_xy(j,m) ), abs ( node_xy(j,m1) ) )

      if ( tol * ( cmax + 1.0D+00 ) &
           < abs ( node_xy(j,m) - node_xy(j,m1) ) ) then
        k = j
        exit
      end if

    end do

    if ( k == 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8TRIS2 - Fatal error!'
      write ( *, '(a,i8)' ) '  Fails for point number I = ', i
      write ( *, '(a,i8)' ) '  M = ', m
      write ( *, '(a,i8)' ) '  M1 = ', m1
      write ( *, '(a,2g14.6)' ) '  NODE_XY(M)  = ', node_xy(1:dim_num,m)
      write ( *, '(a,2g14.6)' ) '  NODE_XY(M1) = ', node_xy(1:dim_num,m1)
      ierr = 224
      stop
    end if

  end do
!
!  Starting from nodes M1 and M2, search for a third point M that
!  makes a "healthy" triangle (M1,M2,M)
!
  m1 = 1
  m2 = 2
  j = 3

  do

    if ( node_num < j ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8TRIS2 - Fatal error!'
      ierr = 225
      stop
    end if

    m = j

    lr = lrline ( node_xy(1,m), node_xy(2,m), node_xy(1,m1), &
      node_xy(2,m1), node_xy(1,m2), node_xy(2,m2), 0.0D+00 )

    if ( lr /= 0 ) then
      exit
    end if

    j = j + 1

  end do
!
!  Set up the triangle information for (M1,M2,M), and for any other
!  triangles you created because points were collinear with M1, M2.
!
  element_num = j - 2

  if ( lr == -1 ) then

    element_node(1,1) = m1
    element_node(2,1) = m2
    element_node(3,1) = m
    element_neighbor(3,1) = -3

    do i = 2, element_num

      m1 = m2
      m2 = i+1

      element_node(1,i) = m1
      element_node(2,i) = m2
      element_node(3,i) = m

      element_neighbor(1,i-1) = -3 * i
      element_neighbor(2,i-1) = i
      element_neighbor(3,i) = i - 1

    end do

    element_neighbor(1,element_num) = -3 * element_num - 1
    element_neighbor(2,element_num) = -5
    ledg = 2
    ltri = element_num

  else

    element_node(1,1) = m2
    element_node(2,1) = m1
    element_node(3,1) = m

    element_neighbor(1,1) = -4

    do i = 2, element_num

      m1 = m2
      m2 = i+1

      element_node(1,i) = m2
      element_node(2,i) = m1
      element_node(3,i) = m

      element_neighbor(3,i-1) = i
      element_neighbor(1,i) = -3 * i - 3
      element_neighbor(2,i) = i - 1

    end do

    element_neighbor(3,element_num) = -3 * element_num
    element_neighbor(2,1) = -3 * element_num - 2
    ledg = 2
    ltri = 1

  end if
!
!  Insert the vertices one at a time from outside the convex hull,
!  determine visible boundary edges, and apply diagonal edge swaps until
!  Delaunay triangulation of vertices (so far) is obtained.
!
  top = 0

  do i = j+1, node_num

    m = i
    m1 = element_node(ledg,ltri)

    if ( ledg <= 2 ) then
      m2 = element_node(ledg+1,ltri)
    else
      m2 = element_node(1,ltri)
    end if

    lr = lrline ( node_xy(1,m), node_xy(2,m), node_xy(1,m1), &
      node_xy(2,m1), node_xy(1,m2), node_xy(2,m2), 0.0D+00 )

    if ( 0 < lr ) then
      rtri = ltri
      redg = ledg
      ltri = 0
    else
      l = -element_neighbor(ledg,ltri)
      rtri = l / 3
      redg = mod ( l, 3 ) + 1
    end if

    call vbedg ( node_xy(1,m), node_xy(2,m), node_num, node_xy, element_num, &
      element_node, element_neighbor, ltri, ledg, rtri, redg )

    n = element_num + 1
    l = -element_neighbor(ledg,ltri)

    do

      t = l / 3
      e = mod ( l, 3 ) + 1
      l = -element_neighbor(e,t)
      m2 = element_node(e,t)

      if ( e <= 2 ) then
        m1 = element_node(e+1,t)
      else
        m1 = element_node(1,t)
      end if

      element_num = element_num + 1
      element_neighbor(e,t) = element_num

      element_node(1,element_num) = m1
      element_node(2,element_num) = m2
      element_node(3,element_num) = m

      element_neighbor(1,element_num) = t
      element_neighbor(2,element_num) = element_num - 1
      element_neighbor(3,element_num) = element_num + 1

      top = top + 1

      if ( node_num < top ) then
        ierr = 8
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8TRIS2 - Fatal error!'
        write ( *, '(a)' ) '  Stack overflow.'
        stop
      end if

      stack(top) = element_num

      if ( t == rtri .and. e == redg ) then
        exit
      end if

    end do

    element_neighbor(ledg,ltri) = -3 * n - 1
    element_neighbor(2,n) = -3 * element_num - 2
    element_neighbor(3,element_num) = -l

    ltri = n
    ledg = 2

    call swapec ( m, top, ltri, ledg, node_num, node_xy, element_num, &
      element_node, element_neighbor, stack, ierr )

    if ( ierr /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8TRIS2 - Fatal error!'
      write ( *, '(a)' ) '  Error return from SWAPEC.'
      stop
    end if

  end do
!
!  Now account for the sorting that we did.
!
  do i = 1, 3
    do j = 1, element_num
      element_node(i,j) = indx ( element_node(i,j) )
    end do
  end do

  call perm_inverse ( node_num, indx )

  call r82vec_permute ( node_num, indx, node_xy )

  return
end
subroutine swapec ( i, top, btri, bedg, node_num, node_xy, element_num, &
  element_node, element_neighbor, stack, ierr )

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
!    Input/output, integer ( kind = 4 ) BTRI, BEDG; on input, if positive, are
!    the triangle and edge indices of a boundary edge whose updated indices
!    must be recorded.  On output, these may be updated because of swaps.
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, real ( kind = 8 ) NODE_XY(2,NODE_NUM), the coordinates of the nodes.
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of triangles.
!
!    Input/output, integer ( kind = 4 ) ELEMENT_NODE(3,ELEMENT_NUM), the
!    triangle incidence list.  May be updated on output because of swaps.
!
!    Input/output, integer ( kind = 4 ) ELEMENT_NEIGHBOR(3,ELEMENT_NUM), the
!    triangle neighbor list; negative values are used for links of the
!    counter-clockwise linked list of boundary edges;  May be updated on output
!    because of swaps.
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

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) element_num

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
  real ( kind = 8 ) node_xy(dim_num,node_num)
  integer ( kind = 4 ) r
  integer ( kind = 4 ) s
  integer ( kind = 4 ) stack(node_num)
  integer ( kind = 4 ) swap
  integer ( kind = 4 ) t
  integer ( kind = 4 ) top
  integer ( kind = 4 ) element_node(3,element_num)
  integer ( kind = 4 ) element_neighbor(3,element_num)
  integer ( kind = 4 ) tt
  integer ( kind = 4 ) u
  real ( kind = 8 ) x
  real ( kind = 8 ) y
!
!  Determine whether triangles in stack are Delaunay, and swap
!  diagonal edge of convex quadrilateral if not.
!
  x = node_xy(1,i)
  y = node_xy(2,i)

  do

    if ( top <= 0 ) then
      exit
    end if

    t = stack(top)
    top = top - 1

    if ( element_node(1,t) == i ) then
      e = 2
      b = element_node(3,t)
    else if ( element_node(2,t) == i ) then
      e = 3
      b = element_node(1,t)
    else
      e = 1
      b = element_node(2,t)
    end if

    a = element_node(e,t)
    u = element_neighbor(e,t)

    if ( element_neighbor(1,u) == t ) then
      f = 1
      c = element_node(3,u)
    else if ( element_neighbor(2,u) == t ) then
      f = 2
      c = element_node(1,u)
    else
      f = 3
      c = element_node(2,u)
    end if

    swap = diaedg ( x, y, node_xy(1,a), node_xy(2,a), node_xy(1,c), &
      node_xy(2,c), node_xy(1,b), node_xy(2,b) )

    if ( swap == 1 ) then

      em1 = e - 1
      em1 = i4_wrap ( em1, 1, 3 )
      ep1 = e + 1
      ep1 = i4_wrap ( ep1, 1, 3 )
      fm1 = f - 1
      fm1 = i4_wrap ( fm1, 1, 3 )
      fp1 = f + 1
      fp1 = i4_wrap ( fp1, 1, 3 )

      element_node(ep1,t) = c
      element_node(fp1,u) = i

      r = element_neighbor(ep1,t)
      s = element_neighbor(fp1,u)

      element_neighbor(ep1,t) = u
      element_neighbor(fp1,u) = t
      element_neighbor(e,t) = s
      element_neighbor(f,u) = r

      if ( 0 < element_neighbor(fm1,u) ) then
        top = top + 1
        stack(top) = u
      end if

      if ( 0 < s ) then

        if ( element_neighbor(1,s) == u ) then
          element_neighbor(1,s) = t
        else if ( element_neighbor(2,s) == u ) then
          element_neighbor(2,s) = t
        else
          element_neighbor(3,s) = t
        end if

        top = top + 1

        if ( node_num < top ) then
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

        do while ( 0 < element_neighbor(ee,tt) )

          tt = element_neighbor(ee,tt)

          if ( element_node(1,tt) == a ) then
            ee = 3
          else if ( element_node(2,tt) == a ) then
            ee = 1
          else
            ee = 2
          end if

        end do

        element_neighbor(ee,tt) = l

      end if

      if ( 0 < r ) then

        if ( element_neighbor(1,r) == t ) then
          element_neighbor(1,r) = u
        else if ( element_neighbor(2,r) == t ) then
          element_neighbor(2,r) = u
        else
          element_neighbor(3,r) = u
        end if

      else

        if ( t == btri .and. ep1 == bedg ) then
          btri = u
          bedg = f
        end if

        l = - ( 3 * u + f - 1 )
        tt = u
        ee = fm1

        do while ( 0 < element_neighbor(ee,tt) )

          tt = element_neighbor(ee,tt)

          if ( element_node(1,tt) == b ) then
            ee = 3
          else if ( element_node(2,tt) == b ) then
            ee = 1
          else
            ee = 2
          end if

        end do

        element_neighbor(ee,tt) = l

      end if

    end if

  end do

  return
end
subroutine triangulation_order3_print ( node_num, element_num, node_xy, &
  element_node, element_neighbor )

!*****************************************************************************80
!
!! TRIANGULATION_ORDER3_PRINT prints information about a triangulation.
!
!  Discussion:
!
!    Triangulations created by R8TRIS2 include extra information encoded
!    in the negative values of ELEMENT_NEIGHBOR.
!
!    Because some of the nodes counted in NODE_NUM may not actually be
!    used in the triangulation, I needed to compute the true number
!    of vertices.  I added this calculation on 13 October 2001.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 June 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of triangles.
!
!    Input, real ( kind = 8 ) NODE_XY(2,NODE_NUM), the coordinates of the nodes.
!
!    Input, integer ( kind = 4 ) ELEMENT_NODE(3,ELEMENT_NUM), the nodes
!    that make up the triangles.
!
!    Input, integer ( kind = 4 ) ELEMENT_NEIGHBOR(3,ELEMENT_NUM), the
!    triangle neighbors on each side.  If there is no triangle neighbor on
!    a particular side, the value of ELEMENT_NEIGHBOR should be negative.
!    If the triangulation data was created by R8TRIS22, then there is more
!    information encoded in the negative values.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ), parameter :: element_order = 3

  integer ( kind = 4 ) boundary_num
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  real ( kind = 8 ) node_xy(dim_num,node_num)
  integer ( kind = 4 ) s
  logical skip
  integer ( kind = 4 ) sp1
  integer ( kind = 4 ) t
  integer ( kind = 4 ) element_node(element_order,element_num)
  integer ( kind = 4 ) element_neighbor(3,element_num)
  integer ( kind = 4 ), allocatable, dimension ( : ) :: vertex_list
  integer ( kind = 4 ) vertex_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TRIANGULATION_ORDER3_PRINT'
  write ( *, '(a)' ) '  Information defining an order3 triangulation.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of nodes is ', node_num

  call r8mat_transpose_print ( dim_num, node_num, node_xy, &
    '  Node coordinates' )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of triangles is ', element_num
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Sets of three nodes are used as vertices of'
  write ( *, '(a)' ) '  the triangles.  For each triangle, the nodes'
  write ( *, '(a)' ) '  are listed in counterclockwise order.'

  call i4mat_transpose_print ( 3, element_num, element_node, &
    '  Triangle nodes:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  On each side of a given triangle, there is either'
  write ( *, '(a)' ) '  another triangle, or a piece of the convex hull.'
  write ( *, '(a)' ) '  For each triangle, we list the indices of the three'
  write ( *, '(a)' ) '  neighbors, or (if negative) the codes of the'
  write ( *, '(a)' ) '  segments of the convex hull.'

  call i4mat_transpose_print ( 3, element_num, element_neighbor, &
    '  Triangle neighbors' )
!
!  Determine the number of vertices.
!
  allocate ( vertex_list(1:3*element_num) )

  vertex_list(1:3*element_num) = reshape ( element_node(1:3,1:element_num), &
    (/ 3 * element_num /) )

  call i4vec_sort_heap_a ( 3*element_num, vertex_list )

  call i4vec_sorted_unique ( 3*element_num, vertex_list, vertex_num )

  deallocate ( vertex_list )
!
!  Determine the number of boundary points.
!
  boundary_num = 2 * vertex_num - element_num - 2

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of boundary points is ', boundary_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The segments that make up the convex hull can be'
  write ( *, '(a)' ) '  determined from the negative entries of the triangle'
  write ( *, '(a)' ) '  neighbor list.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     #   Tri  Side    N1    N2'
  write ( *, '(a)' ) ' '

  skip = .false.

  k = 0

  do i = 1, element_num

    do j = 1, 3

      if ( element_neighbor(j,i) < 0 ) then
        s = - element_neighbor(j,i)
        t = s / 3

        if ( t < 1 .or. element_num < t ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  Sorry, this data does not use the R8TRIS2'
          write ( *, '(a)' ) '  convention for convex hull segments.'
          skip = .true.
          exit
        end if

        s = mod ( s, 3 ) + 1
        k = k + 1
        n1 = element_node(s,t)
        sp1 = s + 1
        sp1 = i4_wrap ( sp1, 1, 3 )
        n2 = element_node(sp1,t)
        write ( *, '(2x,i4,2x,i4,2x,i4,2x,i4,2x,i4)' ) k, t, s, n1, n2
      end if

    end do

    if ( skip ) then
      exit
    end if

  end do

  return
end
subroutine triangulation_search_delaunay ( node_num, node_xy, element_order, &
  element_num, element_node, element_neighbor, p, triangle_index, alpha, &
  beta, gamma, edge, step_num )

!*****************************************************************************80
!
!! TRIANGULATION_SEARCH_DELAUNAY searches a Delaunay triangulation for a point.
!
!  Discussion:
!
!    The algorithm "walks" from one triangle to its neighboring triangle,
!    and so on, until a triangle is found containing point P, or P is found
!    to be outside the convex hull.
!
!    The algorithm computes the barycentric coordinates of the point with
!    respect to the current triangle.  If all three quantities are positive,
!    the point is contained in the triangle.  If the I-th coordinate is
!    negative, then P lies on the far side of edge I, which is opposite
!    from vertex I.  This gives a hint as to where to search next.
!
!    For a Delaunay triangulation, the search is guaranteed to terminate.
!    For other triangulations, a cycle may occur.
!
!    Note the surprising fact that, even for a Delaunay triangulation of
!    a set of nodes, the nearest node to P need not be one of the
!    vertices of the triangle containing P.
!
!    The code can be called for triangulations of any order, but only
!    the first three nodes in each triangle are considered.  Thus, if
!    higher order triangles are used, and the extra nodes are intended
!    to give the triangle a polygonal shape, these will have no effect,
!    and the results obtained here might be misleading.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2012
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
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, real ( kind = 8 ) NODE_XY(2,NODE_NUM), the coordinates of the nodes.
!
!    Input, integer ( kind = 4 ) ELEMENT_ORDER, the order of the triangles.
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of triangles.
!
!    Input, integer ( kind = 4 ) ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM),
!    the nodes that make up each triangle.
!
!    Input, integer ( kind = 4 ) ELEMENT_NEIGHBOR(3,ELEMENT_NUM), the
!    triangle neighbor list.
!
!    Input, real ( kind = 8 ) P(2), the coordinates of a point.
!
!    Output, integer ( kind = 4 ) TRIANGLE_INDEX, the index of the triangle
!    where the search ended.  If a cycle occurred, then TRIANGLE_INDEX = -1.
!
!    Output, real ( kind = 8 ) ALPHA, BETA, GAMMA, the barycentric
!    coordinates of the point relative to triangle TRIANGLE_INDEX.
!
!    Output, integer ( kind = 4 ) EDGE, indicates the position of the point P in
!    triangle TRIANGLE_INDEX:
!    0, the interior or boundary of the triangle;
!    -1, outside the convex hull of the triangulation, past edge 1;
!    -2, outside the convex hull of the triangulation, past edge 2;
!    -3, outside the convex hull of the triangulation, past edge 3.
!
!    Output, integer ( kind = 4 ) STEP_NUM, the number of steps.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ) element_order

  integer ( kind = 4 ) a
  real ( kind = 8 ) alpha
  integer ( kind = 4 ) b
  real ( kind = 8 ) beta
  integer ( kind = 4 ) c
  real ( kind = 8 ) det
  real ( kind = 8 ) dxp
  real ( kind = 8 ) dxa
  real ( kind = 8 ) dxb
  real ( kind = 8 ) dyp
  real ( kind = 8 ) dya
  real ( kind = 8 ) dyb
  integer ( kind = 4 ) edge
  real ( kind = 8 ) gamma
  real ( kind = 8 ) node_xy(dim_num,node_num)
  real ( kind = 8 ) p(dim_num)
  integer ( kind = 4 ) step_num
  integer ( kind = 4 ) element_node(element_order,element_num)
  integer ( kind = 4 ) triangle_index
  integer ( kind = 4 ), save :: triangle_index_save = -1
  integer ( kind = 4 ) element_neighbor(3,element_num)
!
!  If possible, start with the previous successful value of TRIANGLE_INDEX.
!
  if ( triangle_index_save < 1 .or. element_num < triangle_index_save ) then
    triangle_index = ( element_num + 1 ) / 2
  else
    triangle_index = triangle_index_save
  end if

  step_num = - 1
  edge = 0

  do

    step_num = step_num + 1

    if ( element_num < step_num ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TRIANGULATION_SEARCH_DELAUNAY - Fatal error!'
      write ( *, '(a)' ) '  The algorithm seems to be cycling.'
      triangle_index = -1
      edge = -1
      stop
    end if
!
!  Get the nodes of triangle TRIANGLE_INDEX.
!
    a = element_node(1,triangle_index)
    b = element_node(2,triangle_index)
    c = element_node(3,triangle_index)
!
!  Using vertex C as a base, compute the distances to vertices A and B,
!  and the point P.
!
    dxa = node_xy(1,a) - node_xy(1,c)
    dya = node_xy(2,a) - node_xy(2,c)

    dxb = node_xy(1,b) - node_xy(1,c)
    dyb = node_xy(2,b) - node_xy(2,c)

    dxp = p(1)         - node_xy(1,c)
    dyp = p(2)         - node_xy(2,c)

    det = dxa * dyb - dya * dxb
!
!  Compute the barycentric coordinates of the point P with respect
!  to this triangle.
!
    alpha = ( dxp * dyb - dyp * dxb ) / det
    beta =  ( dxa * dyp - dya * dxp ) / det
    gamma = 1.0D+00 - alpha - beta
!
!  If the barycentric coordinates are all positive, then the point
!  is inside the triangle and we're done.
!
    if ( 0.0D+00 <= alpha .and. &
         0.0D+00 <= beta  .and. &
         0.0D+00 <= gamma ) then
      exit
    end if
!
!  At least one barycentric coordinate is negative.
!
!  If there is a negative barycentric coordinate for which there exists
!  an opposing triangle neighbor closer to the point, move to that triangle.
!
!  (Two coordinates could be negative, in which case we could go for the
!  most negative one, or the most negative one normalized by the actual
!  distance it represents).
!
    if ( alpha < 0.0D+00 .and. 0 < element_neighbor(2,triangle_index) ) then
      triangle_index = element_neighbor(2,triangle_index)
      cycle
    else if ( beta < 0.0D+00 .and. &
      0 < element_neighbor(3,triangle_index) ) then
      triangle_index = element_neighbor(3,triangle_index)
      cycle
    else if ( gamma < 0.0D+00 .and. &
      0 < element_neighbor(1,triangle_index) ) then
      triangle_index = element_neighbor(1,triangle_index)
      cycle
    end if
!
!  All negative barycentric coordinates correspond to vertices opposite
!  sides on the convex hull.
!
!  Note the edge and exit.
!
    if ( alpha < 0.0D+00 ) then
      edge = -2
      exit
    else if ( beta < 0.0D+00 ) then
      edge = -3
      exit
    else if ( gamma < 0.0D+00 ) then
      edge = -1
      exit
    end if

  end do

  triangle_index_save = triangle_index

  return
end
subroutine vbedg ( x, y, node_num, node_xy, element_num, element_node, &
  element_neighbor, ltri, ledg, rtri, redg )

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
!    Input, real ( kind = 8 ) X, Y, the coordinates of a point outside the
!    convex hull of the current triangulation.
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, real ( kind = 8 ) NODE_XY(2,NODE_NUM), the coordinates of the nodes.
!
!    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of triangles.
!
!    Input, integer ( kind = 4 ) ELEMENT_NODE(3,ELEMENT_NUM), the triangle
!    incidence list.
!
!    Input, integer ( kind = 4 ) ELEMENT_NEIGHBOR(3,ELEMENT_NUM), the
!    triangle neighbor list; negative values are used for links of a
!    counter clockwise linked list of boundary edges;
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

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) element_num

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) e
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) l
  logical              ldone
  integer ( kind = 4 ) ledg
  integer ( kind = 4 ) lr
  integer ( kind = 4 ) lrline
  integer ( kind = 4 ) ltri
  real ( kind = 8 ) node_xy(2,node_num)
  integer ( kind = 4 ) redg
  integer ( kind = 4 ) rtri
  integer ( kind = 4 ) t
  integer ( kind = 4 ) element_node(3,element_num)
  integer ( kind = 4 ) element_neighbor(3,element_num)
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

    l = -element_neighbor(redg,rtri)
    t = l / 3
    e = mod ( l, 3 ) + 1
    a = element_node(e,t)

    if ( e <= 2 ) then
      b = element_node(e+1,t)
    else
      b = element_node(1,t)
    end if

    lr = lrline ( x, y, node_xy(1,a), node_xy(2,a), node_xy(1,b), &
      node_xy(2,b), 0.0D+00 )

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

    b = element_node(e,t)
    e = e - 1
    e = i4_wrap ( e, 1, 3 )

    do while ( 0 < element_neighbor(e,t) )

      t = element_neighbor(e,t)

      if ( element_node(1,t) == b ) then
        e = 3
      else if ( element_node(2,t) == b ) then
        e = 1
      else
        e = 2
      end if

    end do

    a = element_node(e,t)

    lr = lrline ( x, y, node_xy(1,a), node_xy(2,a), node_xy(1,b), &
      node_xy(2,b), 0.0D+00 )

    if ( lr <= 0 ) then
      exit
    end if

  end do

  ltri = t
  ledg = e

  return
end
