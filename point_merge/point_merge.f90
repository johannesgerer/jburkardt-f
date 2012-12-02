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
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) k
  real    ( kind = 4 ) r
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
    seed = seed + i4_huge
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
subroutine i4vec_print ( n, a, title )

!*****************************************************************************80
!
!! I4VEC_PRINT prints an I4VEC.
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
!    02 May 2010
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

  integer   ( kind = 4 ) n

  integer   ( kind = 4 ) a(n)
  integer   ( kind = 4 ) i
  character ( len = * )  title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i8,a,2x,i12)' ) i, ':', a(i)
  end do

  return
end
subroutine point_radial_tol_unique_count ( m, n, a, tol, seed, unique_num )

!*****************************************************************************80
!
!! POINT_RADIAL_TOL_UNIQUE_COUNT counts the tolerably unique points.
!
!  Discussion:
!
!    The input data is an M x N array A, representing the M-dimensional
!    coordinates of N points.
!
!    The output is the number of tolerably unique points in the list.
!
!    This program performs the same task as POINT_TOL_UNIQUE_COUNT.
!    But that program is guaranteed to use N^2 comparisons.
!
!    It is hoped that this function, on the other hand, will tend
!    to use O(N) comparisons after an O(NLog(N)) sort.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 July 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows.
!
!    Input, integer ( kind = 4 ) N, the number of columns.
!
!    Input, real ( kind = 8 ) A(M,N), the array of N columns of data.
!
!    Input, real ( kind = 8 ) TOL, a tolerance for equality.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, integer ( kind = 4 ) UNIQUE_NUM, the number of tolerably
!    unique points.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(m,n)
  real    ( kind = 8 ) dist
  integer ( kind = 4 ) hi
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) j
  real    ( kind = 8 ) r(n)
  integer ( kind = 4 ) seed
  real    ( kind = 8 ) tol
  logical              unique(n)
  integer ( kind = 4 ) unique_num
  real    ( kind = 8 ) w(n)
  real    ( kind = 8 ) z(m)

  if ( n <= 0 ) then
    unique_num = 0
    return
  end if
!
!  Assign a base point Z randomly in the convex hull.
!
  call r8vec_uniform_01 ( n, seed, w )

  do i = 1, m
    z(i) = dot_product ( w(1:n), a(i,1:n) ) &
         / sum ( w(1:n) )
  end do
!
!  Compute the radial distance R of each point to Z.
!
  do j = 1, n
    r(j) = sqrt ( sum ( ( a(1:m,j) - z(1:m) )**2 ) )
  end do
!
!  Implicitly sort the R array.
!
  call r8vec_sort_heap_index_a ( n, r, indx )
!
!  To determine if a point I is tolerably unique, we only have to check
!  whether it is distinct from all points J such that R(I) <= R(J) <= R(J)+TOL.
!
  unique_num = 0
  unique(1:n) = .true.

  do i = 1, n

    if ( unique(indx(i)) ) then
!
!  Point INDX(I) is unique, in that no earlier point is near it.
!
      unique_num = unique_num + 1
!
!  Look for later points which are close to point INDX(I)
!  in terms of R.
!
      hi = i

      do while ( hi < n )
        if ( r(indx(i)) + tol < r(indx(hi+1)) ) then
          exit
        end if
        hi = hi + 1
      end do
!
!  Points INDX(I+1) through INDX(HI) have an R value close to 
!  point INDX(I).  Are they truly close to point INDEX(I)?
!
      do j = i + 1, hi
        if ( unique(indx(j)) ) then
          dist = sqrt ( sum ( ( a(1:m,indx(i)) - a(1:m,indx(j)) )**2 ) )
          if ( dist <= tol ) then
            unique(indx(j)) = .false.
          end if
        end if
      end do

    end if

  end do

  return
end
subroutine point_radial_tol_unique_index ( m, n, a, tol, seed, unique_num, &
  undx, xdnu )

!*****************************************************************************80
!
!! POINT_RADIAL_TOL_UNIQUE_INDEX indexes the tolerably unique points.
!
!  Discussion:
!
!    The input data is an M x N array A, representing the M-dimensional
!    coordinates of N points.
!
!    The output is:
!    * the number of tolerably unique points in the list;
!    * the index, in the list of unique items, of the representatives 
!      of each point;
!    * the index, in A, of the tolerably unique representatives.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 July 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows.
!
!    Input, integer ( kind = 4 ) N, the number of columns.
!
!    Input, real ( kind = 8 ) A(M,N), the array of N columns of data.
!
!    Input, real ( kind = 8 ) TOL, a tolerance for equality.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, integer ( kind = 4 ) UNIQUE_NUM, the number of tolerably
!    unique points.
!
!    Output, integer ( kind = 4 ) UNDX(UNIQUE_NUM), the index, in A, of the 
!    tolerably unique points.
!
!    Output, integer ( kind = 4 ) XDNU(N), the index, in UNDX, of the 
!    tolerably unique point that "represents" this point.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(m,n)
  real    ( kind = 8 ) dist
  integer ( kind = 4 ) hi
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) j
  real    ( kind = 8 ) r(n)
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) sui(n)
  real    ( kind = 8 ) tol
  logical              unique(n)
  integer ( kind = 4 ) undx(n)
  integer ( kind = 4 ) unique_num
  real    ( kind = 8 ) w(n)
  integer ( kind = 4 ) xdnu(n)
  real    ( kind = 8 ) z(m)

  if ( n <= 0 ) then
    unique_num = 0
    return
  end if
!
!  Assign a base point Z randomly in the convex hull.
!
  call r8vec_uniform_01 ( n, seed, w )

  do i = 1, m
    z(i) = dot_product ( w(1:n), a(i,1:n) ) &
         / sum ( w(1:n) )
  end do
!
!  Compute the radial distance R of each point to Z.
!
  do j = 1, n
    r(j) = sqrt ( sum ( ( a(1:m,j) - z(1:m) )**2 ) )
  end do
!
!  Implicitly sort the R array.
!
  call r8vec_sort_heap_index_a ( n, r, indx )
!
!  To determine if a point I is tolerably unique, we only have to check
!  whether it is distinct from all points J such that R(I) <= R(J) <= R(J)+TOL.
!
  unique_num = 0
  unique(1:n) = .true.

   do i = 1, n

    if ( unique(indx(i)) ) then
!
!  Point INDX(I) is unique, in that no earlier point is near it.
!
      unique_num = unique_num + 1
      xdnu(indx(i)) = unique_num
      undx(unique_num) = indx(i)
!
!  Look for later points which are close to point INDX(I)
!  in terms of R.
!
      hi = i

      do while ( hi < n )
        if ( r(indx(i)) + tol < r(indx(hi+1)) ) then
          exit
        end if
        hi = hi + 1
      end do
!
!  Points INDX(I+1) through INDX(HI) have an R value close to 
!  point INDX(I).  Are they truly close to point INDEX(I)?
!
      do j = i + 1, hi
        if ( unique(indx(j)) ) then
          dist = sqrt ( sum ( ( a(1:m,indx(i)) - a(1:m,indx(j)) )**2 ) )
          if ( dist <= tol ) then
            unique(indx(j)) = .false.
            xdnu(indx(j)) = xdnu(indx(i))
          end if
        end if
      end do

    end if

  end do

  return
end
subroutine point_radial_unique_count ( m, n, a, seed, unique_num )

!*****************************************************************************80
!
!! POINT_RADIAL_UNIQUE_COUNT counts the unique points.
!
!  Discussion:
!
!    The input data is an M x N array A, representing the M-dimensional
!    coordinates of N points.
!
!    The output is the number of unique points in the list.
!
!    This program performs the same task as POINT_UNIQUE_COUNT, and
!    carries out more work.  Hence, it is not a substitute for
!    POINT_UNIQUE_COUNT.  Instead, it is intended to be a starting point
!    for a similar program which includes a tolerance.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 July 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows.
!
!    Input, integer ( kind = 4 ) N, the number of columns.
!
!    Input, real ( kind = 8 ) A(M,N), the array of N columns of data.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, integer ( kind = 4 ) UNIQUE_NUM, the number of unique points.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) hi
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) lo
  real    ( kind = 8 ) r(n)
  integer ( kind = 4 ) seed
  logical              unique
  integer ( kind = 4 ) unique_index
  integer ( kind = 4 ) unique_num
  real    ( kind = 8 ) w(n)
  real    ( kind = 8 ) z(m)

  if ( n <= 0 ) then
    unique_num = 0
    return
  end if
!
!  Assign a base point Z randomly in the convex hull.
!
  call r8vec_uniform_01 ( n, seed, w )

  do i = 1, m
    z(i) = dot_product ( w(1:n), a(i,1:n) ) &
         / sum ( w(1:n) )
  end do
!
!  Compute the radial distance R of each point to Z.
!
  do j = 1, n
    r(j) = sqrt ( sum ( ( a(1:m,j) - z(1:m) )**2 ) )
  end do
!
!  Implicitly sort the R array.
!
  call r8vec_sort_heap_index_a ( n, r, indx )
!
!  To determine if a point is unique, we only have to check
!  whether it is distinct from all points with the same
!  R value and lower ordering.
!
  unique_num = 0
  hi = 0

  do while ( hi < n )
!
!  Advance LO.
!
    lo = hi + 1
!
!  Extend HI.
!
    hi = lo

    do while ( hi < n )
      if ( r(indx(hi+1)) == r(indx(lo)) ) then
        hi = hi + 1
      else
        exit
      end if
    end do
!
!  Points INDX(LO) through INDX(HI) have same R value.
!
!  Find the unique ones.
!
    unique_num = unique_num + 1

    do j1 = lo + 1, hi
      unique = .true.
      do j2 = lo, j1 - 1
        if ( all ( a(1:m,indx(j2)) == a(1:m,indx(j1)) ) ) then
          unique = .false.
          exit
        end if
      end do
      if ( unique ) then
        unique_num = unique_num + 1
      end if
    end do

  end do

  return
end
subroutine point_tol_unique_count ( m, n, a, tol, unique_num )

!*****************************************************************************80
!
!! POINT_TOL_UNIQUE_COUNT counts the tolerably unique points.
!
!  Discussion:
!
!    The input data is an M x N array A, representing the M-dimensional
!    coordinates of N points.
!
!    This function uses a simple but expensive approach.  The first point
!    is accepted as unique.  Each subsequent point is accepted as unique
!    only if it is at least a tolerance away from all accepted unique points.
!    This means the expected amount of work is O(N^2).
!
!    The output is the number of unique points in the list.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 July 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows.
!
!    Input, integer ( kind = 4 ) N, the number of columns.
!
!    Input, real ( kind = 8 ) A(M,N), the array of N columns of data.
!
!    Input, real ( kind = 8 ) TOL, a tolerance.
!
!    Output, integer ( kind = 4 ) UNIQUE_NUM, the number of unique points.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(m,n)
  real    ( kind = 8 ) dist
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real    ( kind = 8 ) tol
  logical              unique(n)
  integer ( kind = 4 ) unique_num

  unique(1:n) = .true.
  unique_num = n

  do i = 2, n

    do j = 1, i - 1
      if ( unique(j) ) then
        dist = sqrt ( sum ( ( a(1:m,i) - a(1:m,j) )**2 ) )
        if ( dist <= tol ) then
          unique(i) = .false.
          unique_num = unique_num - 1
          exit
        end if
      end if
    end do

  end do

  return
end
subroutine point_tol_unique_index ( m, n, a, tol, unique_num, xdnu )

!*****************************************************************************80
!
!! POINT_TOL_UNIQUE_INDEX indexes the tolerably unique points.
!
!  Discussion:
!
!    This routine uses an algorithm that is O(N^2).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 July 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the dimension of the data values.
!
!    Input, integer ( kind = 4 ) N, the number of data values.
!
!    Input, real ( kind = 8 ) A(M,N), the data values.
!
!    Input, real ( kind = 8 ) TOL, a tolerance for equality.
!
!    Output, integer ( kind = 4 ) UNIQUE_NUM, the number of tolerably
!    unique points.
!
!    Output, integer ( kind = 4 ) XDNU(N), the index, in A, of the tolerably 
!    unique point that "represents" this point.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(m,n)
  real    ( kind = 8 ) dist
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real    ( kind = 8 ) tol
  logical              unique(n)
  integer ( kind = 4 ) unique_num
  integer ( kind = 4 ) xdnu(n)

  unique(1:n) = .true.
  do i = 1, n
    xdnu(i) = i
  end do
  unique_num = n

  i = 1
  xdnu(i) = 1

  do i = 2, n

    do j = 1, i - 1
      if ( unique(j) ) then
        dist = sqrt ( sum ( ( a(1:m,i) - a(1:m,j) )**2 ) )
        if ( dist <= tol ) then
          unique(i) = .false.
          unique_num = unique_num - 1
          xdnu(i) = j
          exit
        end if
      end if
    end do

  end do

  return
end
subroutine point_unique_count ( m, n, a, unique_num )

!*****************************************************************************80
!
!! POINT_UNIQUE_COUNT counts the unique points.
!
!  Discussion:
!
!    The input data is an M x N array A, representing the M-dimensional
!    coordinates of N points.
!
!    The algorithm relies on the fact that, in a sorted list, points that
!    are exactly equal must occur consecutively.
!
!    The output is the number of unique points in the list.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 July 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows.
!
!    Input, integer ( kind = 4 ) N, the number of columns.
!
!    Input, real ( kind = 8 ) A(M,N), the array of N columns of data.
!
!    Output, integer ( kind = 4 ) UNIQUE_NUM, the number of unique points.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) unique_index
  integer ( kind = 4 ) unique_num

  if ( n <= 0 ) then
    unique_num = 0
    return
  end if
!
!  Implicitly sort the array.
!
  call r8col_sort_heap_index_a ( m, n, a, indx )
!
!  Two points are considered equal only if they exactly match.
!  In that case, equal points can only occur as consecutive items
!  in the sorted list.   This makes counting easy.
!
  unique_num = 1
  unique_index = indx(1)

  do j = 2, n

    if ( any ( a(1:m,unique_index) /= a(1:m,indx(j)) ) ) then
      unique_num = unique_num + 1
      unique_index = indx(j)
    end if

  end do

  return
end
subroutine point_unique_index ( m, n, a, unique_num, undx, xdnu )

!*****************************************************************************80
!
!! POINT_UNIQUE_INDEX indexes unique points.
!
!  Discussion:
!
!    An R8COL is an M by N array of R8's, regarded as an array of N columns,
!    each of length M.
!
!    The goal of this routine is to determine a vector UNDX,
!    which points, to the unique elements of A, in sorted order,
!    and a vector XDNU, which identifies, for each entry of A, the index of
!    the unique sorted element of A.
!
!    This is all done with index vectors, so that the elements of
!    A are never moved.
!
!    The first step of the algorithm requires the indexed sorting
!    of A, which creates arrays INDX and XDNI.  (If all the entries
!    of A are unique, then these arrays are the same as UNDX and XDNU.)
!
!    We then use INDX to examine the entries of A in sorted order,
!    noting the unique entries, creating the entries of XDNU and
!    UNDX as we go.
!
!    Once this process has been completed, the object X could be
!    replaced by a compressed object XU, containing the unique entries
!    of X in sorted order, using the formula
!
!      XU(*) = A(UNDX(*)).
!
!    We could then, if we wished, reconstruct the entire vector A, or
!    any element of it, by index, as follows:
!
!      A(I) = XU(XDNU(I)).
!
!    We could then replace A by the combination of XU and XDNU.
!
!    Later, when we need the I-th entry of A, we can locate it as
!    the XDNU(I)-th entry of XU.
!
!    Here is an example of a vector A, the sort and inverse sort
!    index vectors, and the unique sort and inverse unique sort vectors
!    and the compressed unique sorted vector.
!
!      I    A   Indx  Xdni      XU   Undx  Xdnu
!    ----+-----+-----+-----+--------+-----+-----+
!      1 | 11.     1     1 |    11.     1     1
!      2 | 22.     3     5 |    22.     2     2
!      3 | 11.     6     2 |    33.     4     1
!      4 | 33.     9     8 |    55.     5     3
!      5 | 55.     2     9 |                  4
!      6 | 11.     7     3 |                  1
!      7 | 22.     8     6 |                  2
!      8 | 22.     4     7 |                  2
!      9 | 11.     5     4 |                  1
!
!    INDX(2) = 3 means that sorted item(2) is A(3).
!    XDNI(2) = 5 means that A(2) is sorted item(5).
!
!    UNDX(3) = 4 means that unique sorted item(3) is at A(4).
!    XDNU(8) = 2 means that A(8) is at unique sorted item(2).
!
!    XU(XDNU(I))) = A(I).
!    XU(I)        = A(UNDX(I)).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 July 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the dimension of the data values.
!
!    Input, integer ( kind = 4 ) N, the number of data values.
!
!    Input, real ( kind = 8 ) A(M,N), the data values.
!
!    Input, integer ( kind = 4 ) UNIQUE_NUM, the number of unique values 
!    in A.  This value is only required for languages in which the size of
!    UNDX must be known in advance.
!
!    Output, integer ( kind = 4 ) UNDX(UNIQUE_NUM), the UNDX vector.
!
!    Output, integer ( kind = 4 ) XDNU(N), the XDNU vector.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) unique_num

  real    ( kind = 8 ) a(m,n)
  real    ( kind = 8 ) diff
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) undx(unique_num)
  integer ( kind = 4 ) xdnu(n)
!
!  Implicitly sort the array.
!
  call r8col_sort_heap_index_a ( m, n, a, indx )
!
!  Walk through the implicitly sorted array.
!
  i = 1
  j = 1
  undx(j) = indx(i)
  xdnu(indx(i)) = j

  do i = 2, n

    diff = maxval ( abs ( a(1:m,indx(i)) - a(1:m,undx(j)) ) )

    if ( 0.0D+00 < diff ) then
      j = j + 1
      undx(j) = indx(i)
    end if

    xdnu(indx(i)) = j

  end do

  return
end
subroutine r8col_duplicates ( m, n, n_unique, seed, a )

!*****************************************************************************80
!
!! R8COL_DUPLICATES generates an R8COL with some duplicate columns.
!
!  Discussion:
!
!    An R8COL is an M by N array of R8's, regarded as an array of N columns,
!    each of length M.
!
!    This routine generates a random R8COL with a specified number of
!    duplicate columns.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 July 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in each column of A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, integer ( kind = 4 ) N_UNIQUE, the number of unique columns in A.
!    1 <= N_UNIQUE <= N.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = 8 ) A(M,N), the array.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) n_unique
  integer ( kind = 4 ) seed
  real    ( kind = 8 ) temp(m)

  if ( n_unique < 1 .or. n < n_unique ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8COL_DUPLICATES - Fatal error!'
    write ( *, '(a)' ) '  1 <= N_UNIQUE <= N is required.'
    stop
  end if

  call r8mat_uniform_01 ( m, n_unique, seed, a )
!
!  Randomly copy unique columns.
!
  do j1 = n_unique + 1, n
    j2 = i4_uniform ( 1, n_unique, seed )
    a(1:m,j1) = a(1:m,j2)
  end do
!
!  Permute the columns.
!
  do j1 = 1, n
    j2 = i4_uniform ( j1, n, seed )
    temp(1:m) = a(1:m,j1)
    a(1:m,j1) = a(1:m,j2)
    a(1:m,j2) = temp(1:m)
  end do

  return
end
subroutine r8col_sort_heap_index_a ( m, n, a, indx )

!*****************************************************************************80
!
!! R8COL_SORT_HEAP_INDEX_A does an indexed heap ascending sort of an R8COL.
!
!  Discussion:
!
!    An R8COL is an M by N array of R8's, regarded as an array of N columns,
!    each of length M.
!
!    The sorting is not actually carried out.  Rather an index array is
!    created which defines the sorting.  This array may be used to sort
!    or index the array, or to sort or index related arrays keyed on the
!    original array.
!
!    A(*,J1) < A(*,J2) if the first nonzero entry of A(*,J1)-A(*,J2) 
!    is negative.
!
!    Once the index array is computed, the sorting can be carried out
!    "implicitly:
!
!      A(*,INDX(*)) is sorted,
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in each column of A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 8 ) A(M,N), the array.
!
!    Output, integer ( kind = 4 ) INDX(N), the sort index.  The I-th element
!    of the sorted array is column INDX(I).
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(m,n)
  real    ( kind = 8 ) column(m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) indxt
  integer ( kind = 4 ) ir
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l

  if ( n < 1 ) then
    return
  end if

  do i = 1, n
    indx(i) = i
  end do

  if ( n == 1 ) then
    return
  end if

  l = ( n / 2 ) + 1
  ir = n

  do

    if ( 1 < l ) then

      l = l - 1
      indxt = indx(l)
      column(1:m) = a(1:m,indxt)

    else

      indxt = indx(ir)
      column(1:m) = a(1:m,indxt)
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

        call r8vec_compare ( m, a(1:m,indx(j)), a(1:m,indx(j+1)), isgn )

        if ( isgn < 0 ) then
          j = j + 1
        end if

      end if

      call r8vec_compare ( m, column, a(1:m,indx(j)), isgn )

      if ( isgn < 0 ) then
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
subroutine r8mat_transpose_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8MAT_TRANSPOSE_PRINT prints an R8MAT, transposed.
!
!  Discussion:
!
!    An R8MAT is an array of R8 values.
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
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) n

  real      ( kind = 8 ) a(m,n)
  character ( len = * )  title

  call r8mat_transpose_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine r8mat_transpose_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8MAT_TRANSPOSE_PRINT_SOME prints some of an R8MAT, transposed.
!
!  Discussion:
!
!    An R8MAT is an array of R8 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 September 2009
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

  integer   ( kind = 4 ), parameter :: incx = 5
  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) n

  real      ( kind = 8 ) a(m,n)
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

      write ( *, '(i5,a,5a14)' ) j, ':', ( ctemp(i), i = 1, inc )

    end do

  end do

  return
end
subroutine r8mat_uniform_01 ( m, n, seed, r )

!*****************************************************************************80
!
!! R8MAT_UNIFORM_01 returns a unit pseudorandom R8MAT.
!
!  Discussion:
!
!    An R8MAT is an array of R8's.
!
!    For now, the input quantity SEED is an integer variable.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Second Edition,
!    Springer, 1987,
!    ISBN: 0387964673,
!    LC: QA76.9.C65.B73.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, December 1986, pages 362-376.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998,
!    ISBN: 0471134031,
!    LC: T57.62.H37.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, Number 2, 1969, pages 136-143.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns 
!    in the array.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which 
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R(M,N), the array of pseudorandom values.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real    ( kind = 8 ) r(m,n)

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  do j = 1, n

    do i = 1, m

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed < 0 ) then
        seed = seed + i4_huge
      end if

      r(i,j) = real ( seed, kind = 8 ) * 4.656612875D-10

    end do
  end do

  return
end
subroutine r8vec_compare ( n, a1, a2, isgn )

!*****************************************************************************80
!
!! R8VEC_COMPARE compares two R8VEC's.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The lexicographic ordering is used.
!
!  Example:
!
!    Input:
!
!      A1 = ( 2.0, 6.0, 2.0 )
!      A2 = ( 2.0, 8.0, 12.0 )
!
!    Output:
!
!      ISGN = -1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vectors.
!
!    Input, real ( kind = 8 ) A1(N), A2(N), the vectors to be compared.
!
!    Output, integer ( kind = 4 ) ISGN, the results of the comparison:
!    -1, A1 < A2,
!     0, A1 = A2,
!    +1, A1 > A2.
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a1(n)
  real    ( kind = 8 ) a2(n)
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) k

  isgn = 0

  k = 1

  do while ( k <= n )

    if ( a1(k) < a2(k) ) then
      isgn = -1
      return
    else if ( a2(k) < a1(k) ) then
      isgn = + 1
      return
    end if

    k = k + 1

  end do

  return
end
subroutine r8vec_print ( n, a, title )

!*****************************************************************************80
!
!! R8VEC_PRINT prints an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
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
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer   ( kind = 4 ) n

  real      ( kind = 8 ) a(n)
  integer   ( kind = 4 ) i
  character ( len = * )  title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i8,a,1x,g16.8)' ) i, ':', a(i)
  end do

  return
end
subroutine r8vec_sort_heap_index_a ( n, a, indx )

!*****************************************************************************80
!
!! R8VEC_SORT_HEAP_INDEX_A does an indexed heap ascending sort of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    The sorting is not actually carried out.  Rather an index array is
!    created which defines the sorting.  This array may be used to sort
!    or index the array, or to sort or index related arrays keyed on the
!    original array.
!
!    Once the index array is computed, the sorting can be carried out
!    "implicitly:
!
!      A(INDX(I:N)) is sorted,
!
!    or explicitly, by the call
!
!      call r8vec_permute ( n, indx, a )
!
!    after which A(1:N) is sorted.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the array.
!
!    Input, real ( kind = 8 ) A(N), an array to be index-sorted.
!
!    Output, integer ( kind = 4 ) INDX(N), the sort index.  The
!    I-th element of the sorted array is A(INDX(I)).
!
  implicit none

  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(n)
  real    ( kind = 8 ) aval
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) indxt
  integer ( kind = 4 ) ir
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l

  if ( n < 1 ) then
    return
  end if

  do i = 1, n
    indx(i) = i
  end do

  if ( n == 1 ) then
    return
  end if

  l = n / 2 + 1
  ir = n

  do

    if ( 1 < l ) then

      l = l - 1
      indxt = indx(l)
      aval = a(indxt)

    else

      indxt = indx(ir)
      aval = a(indxt)
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
        if ( a(indx(j)) < a(indx(j+1)) ) then
          j = j + 1
        end if
      end if

      if ( aval < a(indx(j)) ) then
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
subroutine r8vec_uniform_01 ( n, seed, r )

!*****************************************************************************80
!
!! R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    For now, the input quantity SEED is an integer variable.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Second Edition,
!    Springer, 1987,
!    ISBN: 0387964673,
!    LC: QA76.9.C65.B73.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, December 1986, pages 362-376.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998,
!    ISBN: 0471134031,
!    LC: T57.62.H37.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, Number 2, 1969, pages 136-143.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which 
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R(N), the vector of pseudorandom values.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real    ( kind = 8 ) r(n)

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + i4_huge
    end if

    r(i) = real ( seed, kind = 8 ) * 4.656612875D-10

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

  character ( len = 8  ) ampm
  integer   ( kind = 4 ) d
  integer   ( kind = 4 ) h
  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) mm
  character ( len = 9  ), parameter, dimension(12) :: month = (/ &
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
