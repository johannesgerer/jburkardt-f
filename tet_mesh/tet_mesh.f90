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
subroutine i4col_compare ( m, n, a, i, j, isgn )

!*****************************************************************************80
!
!! I4COL_COMPARE compares columns I and J of an I4COL.
!
!  Example:
!
!    Input:
!
!      M = 3, N = 4, I = 2, J = 4
!
!      A = (
!        1  2  3  4
!        5  6  7  8
!        9 10 11 12 )
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
!    12 June 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, integer ( kind = 4 ) A(M,N), an array of N columns of vectors 
!    of length M.
!
!    Input, integer ( kind = 4 ) I, J, the columns to be compared.
!    I and J must be between 1 and N.
!
!    Output, integer ( kind = 4 ) ISGN, the results of the comparison:
!    -1, column I < column J,
!     0, column I = column J,
!    +1, column J < column I.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
!
!  Check.
!
  if ( i < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4COL_COMPARE - Fatal error!'
    write ( *, '(a,i6,a)' ) '  Column index I = ', i, ' is less than 1.'
    stop
  end if

  if ( n < i ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4COL_COMPARE - Fatal error!'
    write ( *, '(a,i6,a)' ) '  N = ', n, ' is less than column index I = ', i
    stop
  end if

  if ( j < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4COL_COMPARE - Fatal error!'
    write ( *, '(a,i6,a)' ) '  Column index J = ', j, ' is less than 1.'
    stop
  end if

  if ( n < j ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4COL_COMPARE - Fatal error!'
    write ( *, '(a,i6,a)' ) '  N = ', n, ' is less than column index J = ', j
    stop
  end if

  isgn = 0

  if ( i == j ) then
    return
  end if

  k = 1

  do while ( k <= m )

    if ( a(k,i) < a(k,j) ) then
      isgn = -1
      return
    else if ( a(k,j) < a(k,i) ) then
      isgn = +1
      return
    end if

    k = k + 1

  end do

  return
end
subroutine i4col_sort_a ( m, n, a )

!*****************************************************************************80
!
!! I4COL_SORT_A ascending sorts an I4COL.
!
!  Discussion:
!
!    In lexicographic order, the statement "X < Y", applied to two real
!    vectors X and Y of length M, means that there is some index I, with
!    1 <= I <= M, with the property that
!
!      X(J) = Y(J) for J < I,
!    and
!      X(I) < Y(I).
!
!    In other words, the first time they differ, X is smaller.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 September 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of A, and the length of
!    a vector of data.
!
!    Input, integer ( kind = 4 ) N, the number of columns of A.
!
!    Input/output, integer ( kind = 4 ) A(M,N).
!    On input, the array of N columns of M-vectors.
!    On output, the columns of A have been sorted in ascending
!    lexicographic order.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j

  if ( m <= 0 ) then
    return
  end if

  if ( n <= 1 ) then
    return
  end if
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

      call i4col_swap ( m, n, a, i, j )
!
!  Compare the I and J objects.
!
    else if ( indx < 0 ) then

      call i4col_compare ( m, n, a, i, j, isgn )

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
subroutine i4col_sort2_a ( m, n, a )

!*****************************************************************************80
!
!! I4COL_SORT2_A ascending sorts the elements of each column of an I4COL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of A.
!
!    Input, integer ( kind = 4 ) N, the number of columns of A, and the length
!    of a vector of data.
!
!    Input/output, integer ( kind = 4 ) A(M,N).
!    On input, the array of N columns of M vectors.
!    On output, the elements of each column of A have been sorted in ascending
!    order.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) col
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  integer ( kind = 4 ) row
  integer ( kind = 4 ) t

  if ( m <= 1 ) then
    return
  end if

  if ( n <= 0 ) then
    return
  end if
!
!  Initialize.
!
  do col = 1, n

    i = 0
    indx = 0
    isgn = 0
    j = 0
!
!  Call the external heap sorter.
!
    do

      call sort_heap_external ( m, indx, i, j, isgn )
!
!  Interchange the I and J objects.
!
      if ( 0 < indx ) then

        t        = a(i,col)
        a(i,col) = a(j,col)
        a(j,col) = t
!
!  Compare the I and J objects.
!
      else if ( indx < 0 ) then

        if ( a(j,col) < a(i,col) ) then
          isgn = +1
        else
          isgn = -1
        end if

      else if ( indx == 0 ) then

        exit

      end if

    end do

  end do

  return
end
subroutine i4col_sorted_unique_count ( m, n, a, unique_num )

!*****************************************************************************80
!
!! I4COL_SORTED_UNIQUE_COUNT counts unique elements in an I4COL.
!
!  Discussion:
!
!    An I4COL is an M by N array of integer values, regarded
!    as an array of N columns of length M.
!
!    The columns of the array may be ascending or descending sorted.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, integer ( kind = 4 ) A(M,N), a sorted array, containing
!    N columns of data.
!
!    Output, integer ( kind = 4 ) UNIQUE_NUM, the number of unique columns.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) unique_num

  if ( n <= 0 ) then
    unique_num = 0
    return
  end if

  unique_num = 1
  j1 = 1

  do j2 = 2, n

    if ( any ( a(1:m,j1) /= a(1:m,j2) ) ) then
      unique_num = unique_num + 1
      j1 = j2
    end if

  end do

  return
end
subroutine i4col_swap ( m, n, a, j1, j2 )

!*****************************************************************************80
!
!! I4COL_SWAP swaps columns J1 and J2 of an I4COL.
!
!  Example:
!
!    Input:
!
!      M = 3, N = 4, J1 = 2, J2 = 4
!
!      A = (
!        1  2  3  4
!        5  6  7  8
!        9 10 11 12 )
!
!    Output:
!
!      A = (
!        1  4  3  2
!        5  8  7  6
!        9 12 11 10 )
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
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns in 
!    the array.
!
!    Input/output, integer ( kind = 4 ) A(M,N), an array of N columns 
!    of length M.
!
!    Input, integer ( kind = 4 ) J1, J2, the columns to be swapped.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) col(m)
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2

  if ( j1 < 1 .or. n < j1 .or. j2 < 1 .or. n < j2 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4COL_SWAP - Fatal error!'
    write ( *, '(a)' ) '  J1 or J2 is out of bounds.'
    write ( *, '(a,i6)' ) '  J1 =    ', j1
    write ( *, '(a,i6)' ) '  J2 =    ', j2
    write ( *, '(a,i6)' ) '  N =     ', n
    stop

  end if

  if ( j1 == j2 ) then
    return
  end if

  col(1:m)  = a(1:m,j1)
  a(1:m,j1) = a(1:m,j2)
  a(1:m,j2) = col(1:m)

  return
end
subroutine i4i4_sort_a ( i1, i2, j1, j2 )

!*****************************************************************************80
!
!! I4I4_SORT_A ascending sorts a pair of I4's.
!
!  Discussion:
!
!    The program allows the reasonable call:
!
!      call i4i4_sort_a ( i1, i2, i1, i2 )
!
!    and this will return the reasonable result.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I1, I2, the values to sort.
!
!    Output, integer ( kind = 4 ) J1, J2, the sorted values.
!
  implicit none

  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
!
!  Copy arguments, so that the user can make "reasonable" calls like:
!
!    call i4i4_sort_a ( i1, i2, i1, i2 )
!
  k1 = i1
  k2 = i2

  j1 = min ( k1, k2 )
  j2 = max ( k1, k2 )

  return
end
subroutine i4i4i4_sort_a ( i1, i2, i3, j1, j2, j3 )

!*****************************************************************************80
!
!! I4I4I4_SORT_A ascending sorts a triple of I4's.
!
!  Discussion:
!
!    The program allows the reasonable call:
!
!      call i4i4i4_sort_a ( i1, i2, i3, i1, i2, i3 )
! 
!    and this will return the reasonable result.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I1, I2, I3, the values to sort.
!
!    Output, integer ( kind = 4 ) J1, J2, J3, the sorted values.
!
  implicit none

  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j3
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  integer ( kind = 4 ) k3
!
!  Copy arguments, so that the user can make "reasonable" calls like:
!
!    call i4i4i4_sort_a ( i1, i2, i3, i1, i2, i3 )
!
  k1 = i1
  k2 = i2
  k3 = i3

  j1 = min ( min ( k1, k2 ), min ( k2, k3 ) )
  j2 = min ( max ( k1, k2 ), &
       min ( max ( k2, k3 ), max ( k3, k1 ) ) )
  j3 = max ( max ( k1, k2 ), max ( k2, k3 ) )

  return
end
subroutine i4mat_transpose_print ( m, n, a, title )

!*****************************************************************************80
!
!! I4MAT_TRANSPOSE_PRINT prints an I4MAT, transposed.
!
!  Discussion:
!
!    An I4MAT is a rectangular array of I4 values.
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
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  character ( len = * )  title

  call i4mat_transpose_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine i4mat_transpose_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! I4MAT_TRANSPOSE_PRINT_SOME prints some of the transpose of an I4MAT.
!
!  Discussion:
!
!    An I4MAT is a rectangular array of integer values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
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
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 10
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  character ( len = 8 )  ctemp(incx)
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
      write ( ctemp(i2), '(i8)' ) i
    end do

    write ( *, '(''  Row '',10a8)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Col'
    write ( *, '(a)' ) ' '

    j2lo = max ( jlo, 1 )
    j2hi = min ( jhi, n )

    do j = j2lo, j2hi

      do i2 = 1, inc

        i = i2lo - 1 + i2

        write ( ctemp(i2), '(i8)' ) a(i,j)

      end do

      write ( *, '(i5,1x,10a8)' ) j, ( ctemp(i), i = 1, inc )

    end do

  end do

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
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i8,2x,i12)' ) i, a(i)
  end do

  return
end
subroutine mesh_base_one ( node_num, element_order, element_num, element_node )

!*****************************************************************************80
!
!! MESH_BASE_ONE ensures that the element definition is one-based.
!
!  Discussion:
!
!    The ELEMENT_NODE array contains nodes indices that form elements.
!    The convention for node indexing might start at 0 or at 1.
!    Since a FORTRAN90 program will naturally assume a 1-based indexing, it is
!    necessary to check a given element definition and, if it is actually
!    0-based, to convert it.
!
!    This function attempts to detect 9-based node indexing and correct it.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 October 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, int NODE_NUM, the number of nodes.
!
!    Input, int ELEMENT_ORDER, the order of the elements.
!
!    Input, int ELEMENT_NUM, the number of elements.
!
!    Input/output, int ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM), the element
!    definitions.
!
  implicit none

  integer ( kind = 4 ) element_num
  integer ( kind = 4 ) element_order

  integer ( kind = 4 ) element
  integer ( kind = 4 ) element_node(element_order,element_num)
  integer ( kind = 4 ) node
  integer ( kind = 4 ) node_max
  integer ( kind = 4 ) node_min
  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) order

  node_min = node_num + 1
  node_max = -1

  node_min = minval ( element_node(1:element_order,1:element_num) )
  node_max = maxval ( element_node(1:element_order,1:element_num) )

  if ( node_min == 0 .and. node_max == node_num - 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' )'MESH_BASE_ONE:'
    write ( *, '(a)' )'  The element indexing appears to be 0-based!'
    write ( *, '(a)' )'  This will be converted to 1-based.'
    element_node(1:element_order,1:element_num) = &
      element_node(1:element_order,1:element_num) + 1
  else if ( node_min == 1 .and. node_max == node_num  ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' )'MESH_BASE_ONE:'
    write ( *, '(a)' )'  The element indexing appears to be 1-based!'
    write ( *, '(a)' )'  No conversion is necessary.'
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MESH_BASE_ONE - Warning!'
    write ( *, '(a)' ) '  The element indexing is not of a recognized type.'
    write ( *, '(a,i8)' ) '  NODE_MIN = ', node_min
    write ( *, '(a,i8)' ) '  NODE_MAX = ', node_max
    write ( *, '(a,i8)' ) '  NODE_NUM = ', node_num
  end if

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
!    Second Edition,
!    Springer, 1987,
!    ISBN: 0387964673.
!
!    Pierre LEcuyer,
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
function r8mat_det_4d ( a )

!*****************************************************************************80
!
!! R8MAT_DET_4D computes the determinant of a 4 by 4 R8MAT.
!
!  Discussion:
!
!    An R8MAT is a two dimensional matrix of double precision real values.
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
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A(4,4), the matrix whose determinant is desired.
!
!    Output, real ( kind = 8 ) R8MAT_DET_4D, the determinant of the matrix.
!
  implicit none

  real ( kind = 8 ) a(4,4)
  real ( kind = 8 ) r8mat_det_4d

  r8mat_det_4d = &
         a(1,1) * ( &
             a(2,2) * ( a(3,3) * a(4,4) - a(3,4) * a(4,3) ) &
           - a(2,3) * ( a(3,2) * a(4,4) - a(3,4) * a(4,2) ) &
           + a(2,4) * ( a(3,2) * a(4,3) - a(3,3) * a(4,2) ) ) &
       - a(1,2) * ( &
             a(2,1) * ( a(3,3) * a(4,4) - a(3,4) * a(4,3) ) &
           - a(2,3) * ( a(3,1) * a(4,4) - a(3,4) * a(4,1) ) &
           + a(2,4) * ( a(3,1) * a(4,3) - a(3,3) * a(4,1) ) ) &
       + a(1,3) * ( &
             a(2,1) * ( a(3,2) * a(4,4) - a(3,4) * a(4,2) ) &
           - a(2,2) * ( a(3,1) * a(4,4) - a(3,4) * a(4,1) ) &
           + a(2,4) * ( a(3,1) * a(4,2) - a(3,2) * a(4,1) ) ) &
       - a(1,4) * ( &
             a(2,1) * ( a(3,2) * a(4,3) - a(3,3) * a(4,2) ) &
           - a(2,2) * ( a(3,1) * a(4,3) - a(3,3) * a(4,1) ) &
           + a(2,3) * ( a(3,1) * a(4,2) - a(3,2) * a(4,1) ) )

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
!    12 September 2004
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
!    26 March 2005
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

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

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

        if ( a(i,j) == real ( int ( a(i,j) ), kind = 8 ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) a(i,j)
        else
          write ( ctemp(j2), '(g14.6)' ) a(i,j)
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j), j = 1, inc )

    end do

  end do

  return
end
subroutine r8mat_solve ( n, rhs_num, a, info )

!*****************************************************************************80
!
!! R8MAT_SOLVE uses Gauss-Jordan elimination to solve an N by N linear system.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 August 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) RHS_NUM, the number of right hand sides.  
!    RHS_NUM must be at least 0.
!
!    Input/output, real ( kind = 8 ) A(N,N+RHS_NUM), contains in rows and
!    columns 1 to N the coefficient matrix, and in columns N+1 through
!    N+rhs_num, the right hand sides.  On output, the coefficient matrix
!    area has been destroyed, while the right hand sides have
!    been overwritten with the corresponding solutions.
!
!    Output, integer ( kind = 4 ) INFO, singularity flag.
!    0, the matrix was not singular, the solutions were computed;
!    J, factorization failed on step J, and the solutions could not
!    be computed.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) rhs_num

  real ( kind = 8 ) a(n,n+rhs_num)
  real ( kind = 8 ) apivot
  real ( kind = 8 ) factor
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipivot
  integer ( kind = 4 ) j
  real ( kind = 8 ) t(n+rhs_num)

  info = 0

  do j = 1, n
!
!  Choose a pivot row.
!
    ipivot = j
    apivot = a(j,j)

    do i = j + 1, n
      if ( abs ( apivot ) < abs ( a(i,j) ) ) then
        apivot = a(i,j)
        ipivot = i
      end if
    end do

    if ( apivot == 0.0D+00 ) then
      info = j
      return
    end if
!
!  The pivot row moves into the J-th row.
!
    if ( ipivot /= j ) then
      t(       1:n+rhs_num) = a(ipivot,1:n+rhs_num)
      a(ipivot,1:n+rhs_num) = a(j,     1:n+rhs_num)
      a(j,     1:n+rhs_num) = t(       1:n+rhs_num)
    end if
!
!  A(J,J) becomes 1.
!
    a(j,j) = 1.0D+00
    a(j,j+1:n+rhs_num) = a(j,j+1:n+rhs_num) / apivot
!
!  A(I,J) becomes 0.
!
    do i = 1, n

      if ( i /= j ) then
        factor = a(i,j)
        a(i,j) = 0.0D+00
        a(i,j+1:n+rhs_num) = a(i,j+1:n+rhs_num) - factor * a(j,j+1:n+rhs_num)
      end if

    end do

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
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
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
!    An R8MAT is a two dimensional matrix of double precision real values.
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
subroutine r8mat_uniform_01 ( m, n, seed, r )

!*****************************************************************************80
!
!! R8MAT_UNIFORM_01 fills an R8MAT with unit pseudorandom numbers.
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
!    11 August 2004
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
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns in 
!    the array.
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
  real ( kind = 8 ) r(m,n)

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
subroutine r8vec_cross_3d ( v1, v2, v3 )

!*****************************************************************************80
!
!! R8VEC_CROSS_3D computes the cross product of two R8VEC's in 3D.
!
!  Discussion:
!
!    The cross product in 3D can be regarded as the determinant of the
!    symbolic matrix:
!
!          |  i  j  k |
!      det | x1 y1 z1 |
!          | x2 y2 z2 |
!
!      = ( y1 * z2 - z1 * y2 ) * i
!      + ( z1 * x2 - x1 * z2 ) * j
!      + ( x1 * y2 - y1 * x2 ) * k
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) V1(3), V2(3), the two vectors.
!
!    Output, real ( kind = 8 ) V3(3), the cross product vector.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  real ( kind = 8 ) v1(dim_num)
  real ( kind = 8 ) v2(dim_num)
  real ( kind = 8 ) v3(dim_num)

  v3(1) = v1(2) * v2(3) - v1(3) * v2(2)
  v3(2) = v1(3) * v2(1) - v1(1) * v2(3)
  v3(3) = v1(1) * v2(2) - v1(2) * v2(1)

  return
end
function r8vec_length ( dim_num, x )

!*****************************************************************************80
!
!! R8VEC_LENGTH returns the Euclidean length of an R8VEC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, real ( kind = 8 ) X(DIM_NUM), the vector.
!
!    Output, real ( kind = 8 ) R8VEC_LENGTH, the Euclidean length of the vector.
!
  implicit none

  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) r8vec_length
  real ( kind = 8 ) x(dim_num)

  r8vec_length = sqrt ( sum ( ( x(1:dim_num) )**2 ) )

  return
end
subroutine r8vec_mean ( n, a, mean )

!*****************************************************************************80
!
!! R8VEC_MEAN returns the mean of an R8VEC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, real ( kind = 8 ) A(N), the vector whose mean is desired.
!
!    Output, real ( kind = 8 ) MEAN, the mean of the vector entries.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) mean

  mean = sum ( a(1:n) ) / real ( n, kind = 8 )

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
    write ( *, '(2x,i8,2x,g16.8)' ) i, a(i)
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
!    For now, the input quantity SEED is an integer ( kind = 4 ) variable.
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
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real ( kind = 8 ) r(n)

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
      seed = seed + 2147483647
    end if

    r(i) = real ( seed, kind = 8 ) * 4.656612875D-10

  end do

  return
end
subroutine r8vec_variance ( n, a, variance )

!*****************************************************************************80
!
!! R8VEC_VARIANCE returns the variance of an R8VEC.
!
!  Discussion:
!
!    The variance of a vector X of length N is defined as
!
!      mean ( X(1:n) ) = sum ( X(1:n) ) / n
!
!      var ( X(1:n) ) = sum ( ( X(1:n) - mean )**2 ) / ( n - 1 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!    N should be at least 2.
!
!    Input, real ( kind = 8 ) A(N), the vector.
!
!    Output, real ( kind = 8 ) VARIANCE, the variance of the vector.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) mean
  real ( kind = 8 ) variance

  if ( n < 2 ) then

    variance = 0.0D+00

  else

    mean = sum ( a(1:n) ) / real ( n, kind = 8 )

    variance = sum ( ( a(1:n) - mean )**2 )

    variance = variance / real ( n - 1, kind = 8 )

  end if

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
!    routine may be used to sort integers, reals, numbers, names,
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
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert WIlf.
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
!    Input, integer ( kind = 4 ) ISGN, results of comparison of elements I 
!    and J.  (Used only when the previous call returned INDX less than 0).
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
subroutine tet_mesh_neighbor_tets ( tet_order, tet_num, tet_node, &
  tet_neighbor )

!*****************************************************************************80
!
!! TET_MESH_NEIGHBOR_TETS determines tetrahedron neighbors.
!
!  Discussion:
!
!    A tet mesh of a set of nodes can be completely described by
!    the coordinates of the nodes, and the list of nodes that make up
!    each tetrahedron.  In the most common case, four nodes are used.
!    There is also a 10 node case, where nodes are also placed on
!    the midsides of the tetrahedral edges.
!
!    This routine can handle 4 or 10-node tetrahedral meshes.  The
!    10-node case is handled simply by ignoring the six midside nodes,
!    which are presumed to be listed after the vertices.
!
!    The tetrahedron adjacency information records which tetrahedron
!    is adjacent to a given tetrahedron on a particular face.
!
!    This routine creates a data structure recording this information.
!
!    The primary amount of work occurs in sorting a list of 4 * TET_NUM
!    data items.
!
!    The neighbor tetrahedrons are indexed by the face they share with
!    the tetrahedron.
!
!    Each face of the tetrahedron is indexed by the node which is NOT
!    part of the face.  That is:
!
!    * Neighbor 1 shares face 1 defined by nodes 2, 3, 4.
!    * Neighbor 2 shares face 2 defined by nodes 1, 3, 4;
!    * Neighbor 3 shares face 3 defined by nodes 1, 2, 4;
!    * Neighbor 4 shares face 4 defined by nodes 1, 2, 3.
!
!    For instance, if the (transposed) TET_NODE array was:
!
!    Row       1      2      3      4
!    Col
!
!      1       4      3      5      1
!      2       4      2      5      1
!      3       4      7      3      5
!      4       4      7      8      5
!      5       4      6      2      5
!      6       4      6      8      5
!
!    then the (transposed) TET_NEIGHBOR array should be:
!
!    Row       1      2      3      4
!    Col
!
!      1      -1      2     -1      3
!      2      -1      1     -1      5
!      3      -1      1      4     -1
!      4      -1      6      3     -1
!      5      -1      2      6     -1
!      6      -1      4      5     -1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) TET_ORDER, the order of the tetrahedrons.
!
!    Input, integer ( kind = 4 ) TET_NUM, the number of tetrahedrons.
!
!    Input, integer ( kind = 4 ) TET_NODE(TET_ORDER,TET_NUM), the 
!    indices of the nodes.
!
!    Output, integer ( kind = 4 ) TET_NEIGHBOR(4,TET_NUM), the four
!    tetrahedrons that are direct neighbors of a given tetrahedron.  If 
!    there is no neighbor sharing a given face, the index is set to -1.
!
  implicit none

  integer ( kind = 4 ) tet_num
  integer ( kind = 4 ) tet_order

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) c
  integer ( kind = 4 ) face
  integer ( kind = 4 ) face1
  integer ( kind = 4 ) face2
  integer ( kind = 4 ) faces(5,4*tet_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) tet
  integer ( kind = 4 ) tet_neighbor(4,tet_num)
  integer ( kind = 4 ) tet_node(tet_order,tet_num)
  integer ( kind = 4 ) tet1
  integer ( kind = 4 ) tet2
!
!  Step 1.
!  From the list of nodes for tetrahedron T, of the form: (I,J,K,L)
!  construct the four face relations:
!
!    (J,K,L,1,T)
!    (I,K,L,2,T)
!    (I,J,L,3,T)
!    (I,J,K,4,T)
!
!  In order to make matching easier, we reorder each triple of nodes
!  into ascending order.
!
  do tet = 1, tet_num

    i = tet_node(1,tet)
    j = tet_node(2,tet)
    k = tet_node(3,tet)
    l = tet_node(4,tet)

    call i4i4i4_sort_a ( j, k, l, a, b, c )

    faces(1:5,4*(tet-1)+1) = (/ a, b, c, 1, tet /)

    call i4i4i4_sort_a ( i, k, l, a, b, c )

    faces(1:5,4*(tet-1)+2) = (/ a, b, c, 2, tet /)

    call i4i4i4_sort_a ( i, j, l, a, b, c )

    faces(1:5,4*(tet-1)+3) = (/ a, b, c, 3, tet /)

    call i4i4i4_sort_a ( i, j, k, a, b, c )

    faces(1:5,4*(tet-1)+4) = (/ a, b, c, 4, tet /)

  end do
!
!  Step 2. Perform an ascending dictionary sort on the neighbor relations.
!  We only intend to sort on rows 1:3; the routine we call here
!  sorts on rows 1 through 5 but that won't hurt us.
!
!  What we need is to find cases where two tetrahedrons share a face.
!  By sorting the columns of the FACES array, we will put shared faces
!  next to each other.
!
  call i4col_sort_a ( 5, 4*tet_num, faces )
!
!  Step 3. Neighboring tetrahedrons show up as consecutive columns with
!  identical first three entries.  Whenever you spot this happening,
!  make the appropriate entries in TET_NEIGHBOR.
!
  tet_neighbor(1:4,1:tet_num) = -1

  face = 1

  do

    if ( 4 * tet_num <= face ) then
      exit
    end if

    if ( all ( faces(1:3,face) == faces(1:3,face+1) ) ) then
      face1 = faces(4,face)
      tet1 = faces(5,face)
      face2 = faces(4,face+1)
      tet2 = faces(5,face+1)
      tet_neighbor(face1,tet1) = tet2
      tet_neighbor(face2,tet2) = tet1
      face = face + 2
    else
      face = face + 1
    end if

  end do

  return
end
subroutine tet_mesh_node_order ( tet_order, tet_num, tet_node, &
  node_num, node_order )

!*****************************************************************************80
!
!! TET_MESH_NODE_ORDER: determine the order of nodes in a tet mesh.
!
!  Discussion:
!
!    The order of a node is the number of tetrahedrons that use that node
!    as a vertex.
!
!    Tetrahedrons of order 4 or 10 are allowed as input.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) TET_ORDER, the order of the mesh, either 
!    4 or 10.
!
!    Input, integer ( kind = 4 ) TET_NUM, the number of tetrahedrons.
!
!    Input, integer ( kind = 4 ) TET_NODE(TET_ORDER,TET_NUM), the indices
!    of the nodes.
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Output, integer ( kind = 4 ) NODE_ORDER(NODE_NUM), the order of each node.
!
  implicit none

  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) tet_num
  integer ( kind = 4 ) tet_order

  integer ( kind = 4 ) i
  integer ( kind = 4 ) node
  integer ( kind = 4 ) node_order(node_num)
  integer ( kind = 4 ) tet
  integer ( kind = 4 ) tet_node(tet_order,tet_num)

  node_order(1:node_num) = 0

  do tet = 1, tet_num
    do i = 1, tet_order
      node = tet_node(i,tet)
      if ( node < 1 .or. node_num < node ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TET_MESH_NODE_ORDER - Fatal error!'
        write ( *, '(a)' ) '  Illegal entry in TET_NODE.'
        stop
      else
        node_order(node) = node_order(node) + 1
      end if
    end do
  end do

  return
end
subroutine tet_mesh_order4_adj_count ( node_num, tet_num, tet_node, &
  adj_num, adj_row )

!*****************************************************************************80
!
!! TET_MESH_ORDER4_ADJ_COUNT counts the number of nodal adjacencies.
!
!  Discussion:
!
!    Assuming that the tet mesh is to be used in a finite element
!    computation, we declare that two distinct nodes are "adjacent" if and
!    only if they are both included in some tetrahedron.
!
!    It is the purpose of this routine to determine the number of
!    such adjacency relationships.
!
!    The initial count gets only the (I,J) relationships, for which
!    node I is strictly less than node J.  This value is doubled
!    to account for symmetry.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 November 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, integer ( kind = 4 ) TET_NUM, the number of tetrahedrons.
!
!    Input, integer ( kind = 4 ) TET_NODE(4,TET_NUM), the indices of 
!    the nodes.
!
!    Output, integer ( kind = 4 ) ADJ_NUM, the total number of adjacency
!    relationships,
!
!    Output, integer ( kind = 4 ) ADJ_ROW(NODE_NUM+1), the ADJ pointer array.
!
  implicit none

  integer ( kind = 4 ) tet_num
  integer ( kind = 4 ) node_num

  integer ( kind = 4 ) adj_num
  integer ( kind = 4 ) adj_row(node_num+1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) pair(2,6*tet_num)
  integer ( kind = 4 ) pair_num
  integer ( kind = 4 ) pair_unique_num
  integer ( kind = 4 ) tet_node(4,tet_num)
!
!  Each order 4 tetrahedron defines 6 adjacency pairs.
!
  pair(1,            1:  tet_num) = tet_node(1,1:tet_num)
  pair(2,            1:  tet_num) = tet_node(2,1:tet_num)

  pair(1,  tet_num+1:2*tet_num) = tet_node(1,1:tet_num)
  pair(2,  tet_num+1:2*tet_num) = tet_node(3,1:tet_num)

  pair(1,2*tet_num+1:3*tet_num) = tet_node(1,1:tet_num)
  pair(2,2*tet_num+1:3*tet_num) = tet_node(4,1:tet_num)

  pair(1,3*tet_num+1:4*tet_num) = tet_node(2,1:tet_num)
  pair(2,3*tet_num+1:4*tet_num) = tet_node(3,1:tet_num)

  pair(1,4*tet_num+1:5*tet_num) = tet_node(2,1:tet_num)
  pair(2,4*tet_num+1:5*tet_num) = tet_node(4,1:tet_num)

  pair(1,5*tet_num+1:6*tet_num) = tet_node(3,1:tet_num)
  pair(2,5*tet_num+1:6*tet_num) = tet_node(4,1:tet_num)

  pair_num = 6 * tet_num
!
!  Force the nodes of each pair to be listed in ascending order.
!
  call i4col_sort2_a ( 2, pair_num, pair )
!
!  Rearrange the columns in ascending order.
!
  call i4col_sort_a ( 2, pair_num, pair )
!
!  Get the number of unique columns.
!
  call i4col_sorted_unique_count ( 2, pair_num, pair, pair_unique_num )
!
!  The number of adjacencies is TWICE this value, plus the number of nodes.
!
  adj_num = 2 * pair_unique_num
!
!  Now set up the ADJ_ROW counts.
!
  adj_row(1:node_num) = 0

  do k = 1, pair_num

    if ( 1 < k ) then
      if ( pair(1,k-1) == pair(1,k) .and. &
           pair(2,k-1) == pair(2,k) ) then
        cycle
      end if
    end if

    i = pair(1,k)
    j = pair(2,k)

    adj_row(i) = adj_row(i) + 1
    adj_row(j) = adj_row(j) + 1

  end do
!
!  We used ADJ_ROW to count the number of entries in each row.
!  Convert it to pointers into the ADJ array.
!
  adj_row(2:node_num+1) = adj_row(1:node_num)

  adj_row(1) = 1
  do i = 2, node_num+1
    adj_row(i) = adj_row(i-1) + adj_row(i)
  end do

  return
end
subroutine tet_mesh_order4_adj_set ( node_num, tet_num, tet_node, &
  adj_num, adj_row, adj )

!*****************************************************************************80
!
!! TET_MESH_ORDER4_ADJ_SET sets the nodal adjacency matrix.
!
!  Discussion:
!
!    A compressed format is used for the nodal adjacency matrix.
!
!    It is assumed that we know ADJ_NUM, the number of adjacency entries
!    and the ADJ_ROW array, which keeps track of the list of slots
!    in ADJ where we can store adjacency information for each row.
!
!    We essentially repeat the work of TET_MESH_ORDER4_ADJ_COUNT, but
!    now we have a place to store the adjacency information.
!
!    A copy of the ADJ_ROW array is useful, as we can use it to keep track
!    of the next available entry in ADJ for adjacencies associated with
!    a given row.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 November 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, integer ( kind = 4 ) TET_NUM, the number of tetrahedrons.
!
!    Input, integer ( kind = 4 ) TET_NODE(4,TET_NUM), the indices of 
!    the nodes.
!
!    Input, integer ( kind = 4 ) ADJ_NUM, the total number of adjacency 
!    relationships,
!
!    Input, integer ( kind = 4 ) ADJ_ROW(NODE_NUM+1), the ADJ pointer array.
!
!    Output, integer ( kind = 4 ) ADJ(ADJ_NUM), the adjacency information.
!
  implicit none

  integer ( kind = 4 ) adj_num
  integer ( kind = 4 ) tet_num
  integer ( kind = 4 ) node_num

  integer ( kind = 4 ) adj(adj_num)
  integer ( kind = 4 ) adj_row(node_num+1)
  integer ( kind = 4 ) adj_row_copy(node_num+1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) pair(2,6*tet_num)
  integer ( kind = 4 ) pair_num
  integer ( kind = 4 ) tet_node(4,tet_num)
!
!  Each order 4 tetrahedron defines 6 adjacency pairs.
!
  pair(1,            1:  tet_num) = tet_node(1,1:tet_num)
  pair(2,            1:  tet_num) = tet_node(2,1:tet_num)

  pair(1,  tet_num+1:2*tet_num) = tet_node(1,1:tet_num)
  pair(2,  tet_num+1:2*tet_num) = tet_node(3,1:tet_num)

  pair(1,2*tet_num+1:3*tet_num) = tet_node(1,1:tet_num)
  pair(2,2*tet_num+1:3*tet_num) = tet_node(4,1:tet_num)

  pair(1,3*tet_num+1:4*tet_num) = tet_node(2,1:tet_num)
  pair(2,3*tet_num+1:4*tet_num) = tet_node(3,1:tet_num)

  pair(1,4*tet_num+1:5*tet_num) = tet_node(2,1:tet_num)
  pair(2,4*tet_num+1:5*tet_num) = tet_node(4,1:tet_num)

  pair(1,5*tet_num+1:6*tet_num) = tet_node(3,1:tet_num)
  pair(2,5*tet_num+1:6*tet_num) = tet_node(4,1:tet_num)

  pair_num = 6 * tet_num
!
!  Force the nodes of each pair to be listed in ascending order.
!
  call i4col_sort2_a ( 2, pair_num, pair )
!
!  Rearrange the columns in ascending order.
!
  call i4col_sort_a ( 2, pair_num, pair )
!
!  Mark all entries of ADJ so we will know later if we missed one.
!
  adj(1:adj_num) = -1
!
!  Copy the ADJ_ROW array and use it to keep track of the next
!  free entry for each row.
!
  adj_row_copy(1:node_num) = adj_row(1:node_num)
!
!  Now set up the ADJ_ROW counts.
!
  do k = 1, pair_num

    if ( 1 < k ) then
      if ( pair(1,k-1) == pair(1,k) .and. &
           pair(2,k-1) == pair(2,k) ) then
        cycle
      end if
    end if

    i = pair(1,k)
    j = pair(2,k)

    adj(adj_row_copy(i)) = j
    adj_row_copy(i) = adj_row_copy(i) + 1
    adj(adj_row_copy(j)) = i
    adj_row_copy(j) = adj_row_copy(j) + 1

  end do

  return
end
subroutine tet_mesh_order4_boundary_face_count ( tet_num, tet_node, &
  boundary_face_num )

!*****************************************************************************80
!
!! TET_MESH_ORDER4_BOUNDARY_FACE_COUNT counts the number of boundary faces.
!
!  Discussion:
!
!    This routine is given a tet mesh, an abstract list of 
!    quadruples of nodes.  It is assumed that the nodes forming each 
!    face of each tetrahedron are listed in a counterclockwise order, 
!    although the routine should work if the nodes are consistently 
!    listed in a clockwise order as well.
!
!    It is assumed that each face of the tet mesh is either
!    * an INTERIOR face, which is listed twice, once with positive
!      orientation and once with negative orientation, or;
!    * a BOUNDARY face, which will occur only once.
!
!    This routine should work even if the region has holes.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 November 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) TET_NUM, the number of tetrahedrons.
!
!    Input, integer ( kind = 4 ) TET_NODE(4,TET_NUM), the indices of 
!    the nodes.
!
!    Output, integer ( kind = 4 ) BOUNDARY_FACE_NUM, the number of boundary 
!    faces.
!
  implicit none

  integer ( kind = 4 ) tet_num

  integer ( kind = 4 ) boundary_face_num
  integer ( kind = 4 ) face(3,4*tet_num)
  integer ( kind = 4 ) face_num
  integer ( kind = 4 ) interior_face_num
  integer ( kind = 4 ) m
  integer ( kind = 4 ) tet_node(4,tet_num)
  integer ( kind = 4 ) unique_face_num

  m = 3
  face_num = 4 * tet_num
!
!  Set up the face array:
!  (Omit node 1)
!  (Omit node 2)
!  (Omit node 3)
!  (Omit node 4)
!
  face(1:3,            1:  tet_num) = tet_node(2:4,1:tet_num)

  face(1,    tet_num+1:2*tet_num) = tet_node(1,  1:tet_num)
  face(2:3,  tet_num+1:2*tet_num) = tet_node(3:4,1:tet_num)

  face(1:2,2*tet_num+1:3*tet_num) = tet_node(1:2,1:tet_num)
  face(3,  2*tet_num+1:3*tet_num) = tet_node(4,  1:tet_num)

  face(1:3,3*tet_num+1:4*tet_num) = tet_node(1:3,1:tet_num)
!
!  Force the nodes of each face to be listed in ascending order.
!
  call i4col_sort2_a ( m, face_num, face )
!
!  Ascending sort the columns.
!
  call i4col_sort_a ( m, face_num, face )
!
!  Get the number of unique columns.
!
  call i4col_sorted_unique_count ( m, face_num, face, unique_face_num )
!
!  Determine the number of interior and boundary faces.
!
  interior_face_num = 4 * tet_num - unique_face_num

  boundary_face_num = 4 * tet_num - 2 * interior_face_num

  return
end
subroutine tet_mesh_order4_edge_count ( tet_num, tet_node, edge_num )

!*****************************************************************************80
!
!! TET_MESH_ORDER4_EDGE_COUNT counts the number of edges.
!
!  Discussion:
!
!    This routine is given a tet mesh, an abstract list of
!    quadruples of nodes.  Each tetrahedron defines 6 edges; however,
!    assuming that tetrahedrons are touching each other, most edges
!    will be used more than once.  This routine determines the actual
!    number of "geometric" edges associated with the tet mesh.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 November 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) TET_NUM, the number of tetrahedrons.
!
!    Input, integer ( kind = 4 ) TET_NODE(4,TET_NUM), the indices of 
!    the nodes.
!
!    Output, integer ( kind = 4 ) EDGE_NUM, the number of edges.
!
  implicit none

  integer ( kind = 4 ) tet_num

  integer ( kind = 4 ) edge(2,6*tet_num)
  integer ( kind = 4 ) edge_num
  integer ( kind = 4 ) edge_num_raw
  integer ( kind = 4 ) m
  integer ( kind = 4 ) tet_node(4,tet_num)

  m = 3
  edge_num_raw = 6 * tet_num
!
!  Set up the raw edge array:
!
  edge(1,            1:  tet_num) = tet_node(1,1:tet_num)
  edge(2,            1:  tet_num) = tet_node(2,1:tet_num)

  edge(1,  tet_num+1:2*tet_num) = tet_node(1,1:tet_num)
  edge(2,  tet_num+1:2*tet_num) = tet_node(3,1:tet_num)

  edge(1,2*tet_num+1:3*tet_num) = tet_node(1,1:tet_num)
  edge(2,2*tet_num+1:3*tet_num) = tet_node(4,1:tet_num)

  edge(1,3*tet_num+1:4*tet_num) = tet_node(2,1:tet_num)
  edge(2,3*tet_num+1:4*tet_num) = tet_node(3,1:tet_num)

  edge(1,4*tet_num+1:5*tet_num) = tet_node(2,1:tet_num)
  edge(2,4*tet_num+1:5*tet_num) = tet_node(4,1:tet_num)

  edge(1,5*tet_num+1:6*tet_num) = tet_node(3,1:tet_num)
  edge(2,5*tet_num+1:6*tet_num) = tet_node(4,1:tet_num)
!
!  Force the nodes of each face to be listed in ascending order.
!
  call i4col_sort2_a ( m, edge_num_raw, edge )
!
!  Ascending sort the columns.
!
  call i4col_sort_a ( m, edge_num_raw, edge )
!
!  Get the number of unique columns.
!
  call i4col_sorted_unique_count ( m, edge_num_raw, edge, edge_num )

  return
end
subroutine tet_mesh_order4_example_set ( node_num, tet_num, &
  node_xyz, tet_node )

!*****************************************************************************80
!
!! TET_MESH_ORDER4_EXAMPLE_SET sets an example linear tet mesh.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 August 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, integer ( kind = 4 ) TET_NUM, the number of tetrahedrons.
!
!    Output, real ( kind = 8 ) NODE_XYZ(3,NODE_NUM), the node coordinates.
!
!    Output, integer ( kind = 4 ) TET_NODE(4,TET_NUM), the nodes 
!    forming each tet.
!
  implicit none

  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) tet_num

  integer ( kind = 4 ), dimension ( 4, tet_num ) :: tet_node
  real ( kind = 8 ), dimension ( 3, node_num ) :: node_xyz

  tet_node(1:4,1:tet_num) = reshape ( (/  &
     1,   2,   4,  10, &
     2,   4,   5,  10, &
     2,   5,  10,  11, &
     2,   3,   5,  11, &
     4,   5,  10,  13, &
     3,   5,   6,  11, &
     5,  10,  11,  13, &
     4,   5,   7,  13, &
     5,   6,   8,  14, &
     5,   7,   8,  13, &
     6,   8,   9,  14, &
    11,  13,  14,  19, &
    12,  14,  15,  20, &
     3,   6,  11,  12, &
     5,   6,  11,  14, &
     6,   9,  14,  15, &
     6,  11,  12,  14, &
     6,  12,  14,  15, &
     7,   8,  13,  16, &
     5,   8,  13,  14, &
    10,  11,  13,  19, &
     8,   9,  14,  17, &
    11,  12,  14,  20, &
     5,  11,  13,  14, &
     8,  13,  14,  16, &
     9,  14,  15,  17, &
    13,  14,  16,  22, &
     8,  14,  16,  17, &
    14,  15,  17,  23, &
    14,  16,  17,  22, &
     9,  15,  17,  18, &
    15,  17,  18,  23, &
    14,  17,  22,  23, &
    13,  14,  19,  22, &
    11,  14,  19,  20, &
    14,  15,  20,  23, &
    15,  20,  21,  23, &
    21,  23,  24,  29, &
    20,  22,  23,  28, &
    14,  19,  20,  22, &
    15,  18,  23,  24, &
    12,  15,  20,  21, &
    15,  21,  23,  24, &
    16,  17,  22,  25, &
    19,  20,  22,  28, &
    17,  18,  23,  26, &
    20,  21,  23,  29, &
    14,  20,  22,  23, &
    17,  22,  23,  25, &
    18,  23,  24,  26, &
    22,  23,  25,  31, &
    17,  23,  25,  26, &
    23,  24,  26,  32, &
    23,  25,  26,  31, &
    18,  24,  26,  27, &
    24,  26,  27,  32, &
    23,  26,  31,  32, &
    22,  23,  28,  31, &
    20,  23,  28,  29, &
    23,  24,  29,  32, &
    24,  29,  30,  32, &
    30,  32,  33,  38, &
    29,  31,  32,  37, &
    23,  28,  29,  31, &
    24,  27,  32,  33, &
    21,  24,  29,  30, &
    24,  30,  32,  33, &
    25,  26,  31,  34, &
    28,  29,  31,  37, &
    26,  27,  32,  35, &
    29,  30,  32,  38, &
    23,  29,  31,  32, &
    26,  31,  32,  34, &
    27,  32,  33,  35, &
    31,  32,  34,  40, &
    26,  32,  34,  35, &
    32,  33,  35,  41, &
    32,  34,  35,  40, &
    27,  33,  35,  36, &
    33,  35,  36,  41, &
    32,  35,  40,  41, &
    31,  32,  37,  40, &
    29,  32,  37,  38, &
    32,  33,  38,  41, &
    33,  38,  39,  41, &
    39,  41,  42,  47, &
    38,  40,  41,  46, &
    32,  37,  38,  40, &
    33,  36,  41,  42, &
    30,  33,  38,  39, &
    33,  39,  41,  42, &
    34,  35,  40,  43, &
    37,  38,  40,  46, &
    35,  36,  41,  44, &
    38,  39,  41,  47, &
    32,  38,  40,  41, &
    35,  40,  41,  43, &
    36,  41,  42,  44, &
    40,  41,  43,  49, &
    35,  41,  43,  44, &
    41,  42,  44,  50, &
    41,  43,  44,  49, &
    36,  42,  44,  45, &
    42,  44,  45,  50, &
    41,  44,  49,  50, &
    40,  41,  46,  49, &
    38,  41,  46,  47, &
    41,  42,  47,  50, &
    42,  47,  48,  50, &
    48,  50,  51,  56, &
    47,  49,  50,  55, &
    41,  46,  47,  49, &
    42,  45,  50,  51, &
    39,  42,  47,  48, &
    42,  48,  50,  51, &
    43,  44,  49,  52, &
    46,  47,  49,  55, &
    44,  45,  50,  53, &
    47,  48,  50,  56, &
    41,  47,  49,  50, &
    44,  49,  50,  52, &
    45,  50,  51,  53, &
    49,  50,  52,  58, &
    44,  50,  52,  53, &
    50,  51,  53,  59, &
    50,  52,  53,  58, &
    45,  51,  53,  54, &
    51,  53,  54,  59, &
    50,  53,  58,  59, &
    49,  50,  55,  58, &
    47,  50,  55,  56, &
    50,  51,  56,  59, &
    51,  56,  57,  59, &
    50,  55,  56,  58, &
    51,  54,  59,  60, &
    48,  51,  56,  57, &
    51,  57,  59,  60, &
    52,  53,  58,  61, &
    53,  54,  59,  62, &
    50,  56,  58,  59, &
    53,  58,  59,  61, &
    54,  59,  60,  62, &
    53,  59,  61,  62, &
    54,  60,  62,  63  &
  /), (/ 4, tet_num /) )

  node_xyz(1:3,1:node_num) = reshape ( (/ &
  0.0D+00,  0.0D+00,  0.0D+00, &
  0.0D+00,  0.0D+00,  0.5D+00, &
  0.0D+00,  0.0D+00,  1.0D+00, &
  0.0D+00,  0.5D+00,  0.0D+00, &
  0.0D+00,  0.5D+00,  0.5D+00, &
  0.0D+00,  0.5D+00,  1.0D+00, &
  0.0D+00,  1.0D+00,  0.0D+00, &
  0.0D+00,  1.0D+00,  0.5D+00, &
  0.0D+00,  1.0D+00,  1.0D+00, &
  0.5D+00,  0.0D+00,  0.0D+00, &
  0.5D+00,  0.0D+00,  0.5D+00, &
  0.5D+00,  0.0D+00,  1.0D+00, &
  0.5D+00,  0.5D+00,  0.0D+00, &
  0.5D+00,  0.5D+00,  0.5D+00, &
  0.5D+00,  0.5D+00,  1.0D+00, &
  0.5D+00,  1.0D+00,  0.0D+00, &
  0.5D+00,  1.0D+00,  0.5D+00, &
  0.5D+00,  1.0D+00,  1.0D+00, &
  1.0D+00,  0.0D+00,  0.0D+00, &
  1.0D+00,  0.0D+00,  0.5D+00, &
  1.0D+00,  0.0D+00,  1.0D+00, &
  1.0D+00,  0.5D+00,  0.0D+00, &
  1.0D+00,  0.5D+00,  0.5D+00, &
  1.0D+00,  0.5D+00,  1.0D+00, &
  1.0D+00,  1.0D+00,  0.0D+00, &
  1.0D+00,  1.0D+00,  0.5D+00, &
  1.0D+00,  1.0D+00,  1.0D+00, &
  1.5D+00,  0.0D+00,  0.0D+00, &
  1.5D+00,  0.0D+00,  0.5D+00, &
  1.5D+00,  0.0D+00,  1.0D+00, &
  1.5D+00,  0.5D+00,  0.0D+00, &
  1.5D+00,  0.5D+00,  0.5D+00, &
  1.5D+00,  0.5D+00,  1.0D+00, &
  1.5D+00,  1.0D+00,  0.0D+00, &
  1.5D+00,  1.0D+00,  0.5D+00, &
  1.5D+00,  1.0D+00,  1.0D+00, &
  2.0D+00,  0.0D+00,  0.0D+00, &
  2.0D+00,  0.0D+00,  0.5D+00, &
  2.0D+00,  0.0D+00,  1.0D+00, &
  2.0D+00,  0.5D+00,  0.0D+00, &
  2.0D+00,  0.5D+00,  0.5D+00, &
  2.0D+00,  0.5D+00,  1.0D+00, &
  2.0D+00,  1.0D+00,  0.0D+00, &
  2.0D+00,  1.0D+00,  0.5D+00, &
  2.0D+00,  1.0D+00,  1.0D+00, &
  2.5D+00,  0.0D+00,  0.0D+00, &
  2.5D+00,  0.0D+00,  0.5D+00, &
  2.5D+00,  0.0D+00,  1.0D+00, &
  2.5D+00,  0.5D+00,  0.0D+00, &
  2.5D+00,  0.5D+00,  0.5D+00, &
  2.5D+00,  0.5D+00,  1.0D+00, &
  2.5D+00,  1.0D+00,  0.0D+00, &
  2.5D+00,  1.0D+00,  0.5D+00, &
  2.5D+00,  1.0D+00,  1.0D+00, &
  3.0D+00,  0.0D+00,  0.0D+00, &
  3.0D+00,  0.0D+00,  0.5D+00, &
  3.0D+00,  0.0D+00,  1.0D+00, &
  3.0D+00,  0.5D+00,  0.0D+00, &
  3.0D+00,  0.5D+00,  0.5D+00, &
  3.0D+00,  0.5D+00,  1.0D+00, &
  3.0D+00,  1.0D+00,  0.0D+00, &
  3.0D+00,  1.0D+00,  0.5D+00, &
  3.0D+00,  1.0D+00,  1.0D+00  &  
  /), (/ 3, node_num /) )

  return
end
subroutine tet_mesh_order4_example_size ( node_num, tet_num )

!*****************************************************************************80
!
!! TET_MESH_ORDER4_EXAMPLE_SIZE sizes an example linear tet mesh.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 August 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Output, integer ( kind = 4 ) TET_NUM, the number of tetrahedrons.
!
  implicit none

  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) tet_num

  node_num = 63
  tet_num = 144

  return
end
subroutine tet_mesh_order4_refine_compute ( node_num1, tet_num1, node_xyz1, &
  tet_node1, node_num2, tet_num2, edge_data, node_xyz2, tet_node2 )

!*****************************************************************************80
!
!! TET_MESH_ORDER4_REFINE_COMPUTE computes a refined order4 tet mesh.
!
!  Discussion:
!
!    A refined 4-node tet mesh can be derived from a given
!    4-node tet mesh by interpolating nodes at the midpoint of
!    every edge of the mesh.
!
!    The mesh is described indirectly, as the sum of individual
!    tetrahedrons.  A single physical edge may be a logical edge of
!    any number of tetrahedrons.  It is important, however, that a
!    new node be created exactly once for each edge, assigned an index,
!    and associated with every tetrahedron that shares this edge. 
!
!    This routine handles that problem.
!
!    The primary amount of work occurs in sorting a list of 6 * TET_NUM
!    data items, one item for every edge of every tetrahedron.  Each
!    data item records, for a given tetrahedron edge, the global indices
!    of the two endpoints, the local indices of the two endpoints,
!    and the index of the tetrahedron.
!
!    Through careful sorting, it is possible to arrange this data in
!    a way that allows the proper generation of the interpolated nodes.
!
!    Let us add the new nodes and temporarily assign them local indices
!    5 through X, based on the following ordering:
!
!      1, 2, 3, 4, (1+2), (1+3), (1+4), (2+3), (2+4), (3+4).
!
!    Then let us assign these nodes to eight subtetrahedrons as follows:
!
!      1, 5, 6, 7
!      2, 5, 8, 9
!      3, 6, 8, 9
!      4, 7, 9, X
!      5, 6, 7, 9
!      5, 6, 8, 9
!      6, 7, 9, X
!      6, 8, 9, X
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Anwei Liu, Barry Joe,
!    Quality Local Refinement of Tetrahedral Meshes Based
!    on 8-Subtetrahedron Subdivision,
!    Mathematics of Computation,
!    Volume 65, Number 215, July 1996, pages 1183-1200.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM1, the number of nodes in the input
!    mesh.
!
!    Input, integer ( kind = 4 ) TET_NUM1, the number of tetrahedrons in 
!    the input mesh.
!
!    Input, real ( kind = 8 ) NODE_XYZ1(3,NODE_NUM1), the coordinates of
!    the nodes that make up the input mesh.
!
!    Input, integer ( kind = 4 ) TET_NODE1(4,TET_NUM1), the indices of 
!    the nodes in the input mesh.
!
!    Input, integer ( kind = 4 ) NODE_NUM2, the number of nodes for the 
!    refined mesh.
!
!    Input, integer ( kind = 4 ) TET_NUM2, the number of tetrahedrons in the
!    refined mesh.
!
!    Input, integer ( kind = 4 ) EDGE_DATA(5,6*TET_NUM), edge data.
!
!    Output, real ( kind = 8 ) NODE_XYZ2(3,NODE_NUM2), the coordinates of
!    the nodes that make up the output mesh.
!
!    Output, integer ( kind = 4 ) TET_NODE2(4,TET_NUM2), the indices of 
!    the nodes in the output mesh.
!
  implicit none

  integer ( kind = 4 ) node_num1
  integer ( kind = 4 ) node_num2
  integer ( kind = 4 ) tet_num1
  integer ( kind = 4 ) tet_num2

  integer ( kind = 4 ) edge
  integer ( kind = 4 ) edge_data(5,6*tet_num1)
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n1_old
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n2_old
  integer ( kind = 4 ) node
  real ( kind = 8 ) node_xyz1(3,node_num1)
  real ( kind = 8 ) node_xyz2(3,node_num2)
  integer ( kind = 4 ) tet_node1(4,tet_num1)
  integer ( kind = 4 ) tet_node2(4,tet_num2)
  integer ( kind = 4 ) tet1
  integer ( kind = 4 ) tet2
  integer ( kind = 4 ) v
  integer ( kind = 4 ) v1
  integer ( kind = 4 ) v2
!
!  Generate the index and coordinates of the new midside nodes, 
!  and update the tetradehron-node data.
!
  node_xyz2(1:3,1:node_num1) = node_xyz1(1:3,1:node_num1)

  tet_node2(1:4,1:tet_num2) = -1
!
!  The vertices of the input tetrahedron can be assigned now.
!
  do tet1 = 1, tet_num1
    tet_node2(1,(tet1-1)*8+1) = tet_node1(1,tet1)
    tet_node2(1,(tet1-1)*8+2) = tet_node1(2,tet1)
    tet_node2(1,(tet1-1)*8+3) = tet_node1(3,tet1)    
    tet_node2(1,(tet1-1)*8+4) = tet_node1(4,tet1)
  end do

  node = node_num1

  n1_old = -1
  n2_old = -1

  do edge = 1, 6 * tet_num1
!
!  Read the data defining the edge.
!
    n1 = edge_data(1,edge)
    n2 = edge_data(2,edge)
!
!  If this edge is new, create the coordinates and index.
!
    if ( n1 /= n1_old .or. n2 /= n2_old ) then

      node = node + 1

      if ( node_num2 < node ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TET_MESH_ORDER4_REFINE_COMPUTE - Fatal error!'
        write ( *, '(a)' ) '  Node index exceeds NODE_NUM2.'
        stop
      end if

      node_xyz2(1:3,node) = &
        ( node_xyz2(1:3,n1) + node_xyz2(1:3,n2) ) / 2.0D+00

      n1_old = n1
      n2_old = n2

    end if
!
!  Assign the node to the tetrahedron.
!
    v1 = edge_data(3,edge)
    v2 = edge_data(4,edge)
    tet1 = edge_data(5,edge)
!
!  We know the two vertices that bracket this new node.
!  This tells us whether it is new node number 5, 6, 7, 8, 9 or 10.
!  This tells us which of the new subtetrahedrons it belongs to,
!  and what position it occupies.
!
    if ( v1 == 1 .and. v2 == 2 ) then

      tet_node2(2,(tet1-1)*8+1) = node
      tet_node2(2,(tet1-1)*8+2) = node
      tet_node2(1,(tet1-1)*8+5) = node
      tet_node2(1,(tet1-1)*8+6) = node

    else if ( v1 == 1 .and. v2 == 3 ) then

      tet_node2(3,(tet1-1)*8+1) = node
      tet_node2(2,(tet1-1)*8+3) = node
      tet_node2(2,(tet1-1)*8+5) = node
      tet_node2(2,(tet1-1)*8+6) = node
      tet_node2(1,(tet1-1)*8+7) = node
      tet_node2(1,(tet1-1)*8+8) = node

    else if ( v1 == 1 .and. v2 == 4 ) then

      tet_node2(4,(tet1-1)*8+1) = node
      tet_node2(2,(tet1-1)*8+4) = node
      tet_node2(3,(tet1-1)*8+5) = node
      tet_node2(2,(tet1-1)*8+7) = node

    else if ( v1 == 2 .and. v2 == 3 ) then

      tet_node2(3,(tet1-1)*8+2) = node
      tet_node2(3,(tet1-1)*8+3) = node
      tet_node2(3,(tet1-1)*8+6) = node
      tet_node2(2,(tet1-1)*8+8) = node

    else if ( v1 == 2 .and. v2 == 4 ) then

      tet_node2(4,(tet1-1)*8+2) = node
      tet_node2(4,(tet1-1)*8+3) = node
      tet_node2(3,(tet1-1)*8+4) = node
      tet_node2(4,(tet1-1)*8+5) = node
      tet_node2(4,(tet1-1)*8+6) = node
      tet_node2(3,(tet1-1)*8+7) = node
      tet_node2(3,(tet1-1)*8+8) = node

    else if ( v1 == 3 .and. v2 == 4 ) then

      tet_node2(4,(tet1-1)*8+4) = node
      tet_node2(4,(tet1-1)*8+7) = node
      tet_node2(4,(tet1-1)*8+8) = node

    end if

  end do

  return
end
subroutine tet_mesh_order4_refine_size ( node_num1, tet_num1, tet_node1, &
  node_num2, tet_num2, edge_data )

!*****************************************************************************80
!
!! TET_MESH_ORDER4_REFINE_SIZE sizes a refined order 4 tet mesh.
!
!  Discussion:
!
!    A refined tet mesh can be derived from an existing one by interpolating 
!    nodes at the midpoint of every edge of the mesh.
!
!    The mesh is described indirectly, as the sum of individual
!    tetrahedrons.  A single physical edge may be a logical edge of
!    any number of tetrahedrons.  It is important, however, that a
!    new node be created exactly once for each edge, assigned an index,
!    and associated with every tetrahedron that shares this edge. 
!
!    This routine handles that problem.
!
!    The primary amount of work occurs in sorting a list of 6 * TET_NUM
!    data items, one item for every edge of every tetrahedron.  Each
!    data item records, for a given tetrahedron edge, the global indices
!    of the two endpoints, the local indices of the two endpoints,
!    and the index of the tetrahedron.
!
!    Through careful sorting, it is possible to arrange this data in
!    a way that allows the proper generation of the interpolated nodes.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM1, the number of nodes in the 
!    original mesh.
!
!    Input, integer ( kind = 4 ) TET_NUM1, the number of tetrahedrons in the
!    original mesh.
!
!    Input, integer ( kind = 4 ) TET_NODE1(4,TET_NUM1), the indices of 
!    the nodes that form the tetrahedrons in the input mesh.
!
!    Output, integer ( kind = 4 ) NODE_NUM2, the number of nodes in the 
!    refined mesh.
!
!    Output, integer ( kind = 4 ) TET_NUM2, the number of tetrahedrons in the
!    refined mesh.
!
!    Output, integer ( kind = 4 ) EDGE_DATA(5,6*TET_NUM1), edge data.
!
  implicit none

  integer ( kind = 4 ) node_num1
  integer ( kind = 4 ) tet_num1

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) edge
  integer ( kind = 4 ) edge_data(5,6*tet_num1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n1_old
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n2_old
  integer ( kind = 4 ) node_num2
  integer ( kind = 4 ) tet
  integer ( kind = 4 ) tet_node1(4,tet_num1)
  integer ( kind = 4 ) tet_num2
!
!  Step 1.
!  From the list of nodes for tetrahedron T, of the form: (I,J,K,L)
!  construct the six edge relations:
!
!    (I,J,1,2,T)
!    (I,K,1,3,T)
!    (I,L,1,4,T)
!    (J,K,2,3,T)
!    (J,L,2,4,T)
!    (K,L,3,4,T)
!
!  In order to make matching easier, we reorder each pair of nodes
!  into ascending order.
!
  do tet = 1, tet_num1

    i = tet_node1(1,tet)
    j = tet_node1(2,tet)
    k = tet_node1(3,tet)
    l = tet_node1(4,tet)

    call i4i4_sort_a ( i, j, a, b )

    edge_data(1:5,6*(tet-1)+1) = (/ a, b, 1, 2, tet /)

    call i4i4_sort_a ( i, k, a, b )

    edge_data(1:5,6*(tet-1)+2) = (/ a, b, 1, 3, tet /)

    call i4i4_sort_a ( i, l, a, b )

    edge_data(1:5,6*(tet-1)+3) = (/ a, b, 1, 4, tet /)

    call i4i4_sort_a ( j, k, a, b )

    edge_data(1:5,6*(tet-1)+4) = (/ a, b, 2, 3, tet /)

    call i4i4_sort_a ( j, l, a, b )

    edge_data(1:5,6*(tet-1)+5) = (/ a, b, 2, 4, tet /)

    call i4i4_sort_a ( k, l, a, b )

    edge_data(1:5,6*(tet-1)+6) = (/ a, b, 3, 4, tet /)

  end do
!
!  Step 2. Perform an ascending dictionary sort on the neighbor relations.
!  We only intend to sort on rows 1:2; the routine we call here
!  sorts on the full column but that won't hurt us.
!
!  What we need is to find all cases where tetrahedrons share an edge.
!  By sorting the columns of the EDGE_DATA array, we will put shared edges
!  next to each other.
!
  call i4col_sort_a ( 5, 6*tet_num1, edge_data )
!
!  Step 3. All the tetrahedrons which share an edge show up as consecutive
!  columns with identical first two entries.  Figure out how many new
!  nodes there are, and allocate space for their coordinates.
!
  node_num2 = node_num1

  n1_old = -1
  n2_old = -1

  do edge = 1, 6 * tet_num1
    n1 = edge_data(1,edge)
    n2 = edge_data(2,edge)
    if ( n1 /= n1_old .or. n2 /= n2_old ) then
      node_num2 = node_num2 + 1
      n1_old = n1
      n2_old = n2
    end if
  end do

  tet_num2 = 8 * tet_num1

  return
end
subroutine tet_mesh_order4_to_order10_compute ( tet_num, tet_node1, &
  node_num1, node_xyz1, edge_data, tet_node2, node_num2, node_xyz2 )

!*****************************************************************************80
!
!! TET_MESH_ORDER4_TO_ORDER10_COMPUTE: quadratic tet mesh from a linear one.
!
!  Discussion:
!
!    A quadratic (10 node) tet mesh can be derived from a linear
!    (4 node) tet mesh by interpolating nodes at the midpoint of
!    every edge of the mesh.
!
!    The mesh is described indirectly, as the sum of individual
!    tetrahedrons.  A single physical edge may be a logical edge of
!    any number of tetrahedrons.  It is important, however, that a
!    new node be created exactly once for each edge, assigned an index,
!    and associated with every tetrahedron that shares this edge. 
!
!    This routine handles that problem.
!
!    The primary amount of work occurs in sorting a list of 6 * TET_NUM
!    data items, one item for every edge of every tetrahedron.  Each
!    data item records, for a given tetrahedron edge, the global indices
!    of the two endpoints, the local indices of the two endpoints,
!    and the index of the tetrahedron.
!
!    Through careful sorting, it is possible to arrange this data in
!    a way that allows the proper generation of the interpolated nodes.
!
!    The node ordering for the quadratic tetrahedron is somewhat
!    arbitrary.  In the current scheme, the vertices are listed
!    first, followed by the 6 midside nodes.  Each midside node
!    may be identified by the two vertices that bracket it.  Thus,
!    the node ordering may be suggested by:
!
!      1  2  3  4 (1+2) (1+3) (1+4) (2+3) (2+4) (3+4)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) TET_NUM, the number of tetrahedrons in the
!    linear mesh.
!
!    Input, integer ( kind = 4 ) TET_NODE1(4,TET_NUM), the indices of 
!    the nodes in the linear mesh.
!
!    Input, integer ( kind = 4 ) NODE_NUM1, the number of nodes for the 
!    linear mesh.
!
!    Input, real ( kind = 8 ) NODE_XYZ1(3,NODE_NUM1), the coordinates of
!    the nodes that make up the linear mesh.
!
!    Input, integer ( kind = 4 ) EDGE_DATA(5,6*TET_NUM), edge data.
!
!    Output, integer ( kind = 4 ) TET_NODE2(10,TET_NUM), the indices of 
!    the nodes in the quadratic mesh.
!
!    Input, integer ( kind = 4 ) NODE_NUM2, the number of nodes for the 
!    quadratic mesh.
!
!    Output, real ( kind = 8 ) NODE_XYZ2(3,NODE_NUM2), the coordinates of
!    the nodes that make up the quadratic mesh.
!
  implicit none

  integer ( kind = 4 ) node_num1
  integer ( kind = 4 ) node_num2
  integer ( kind = 4 ) tet_num

  integer ( kind = 4 ) edge
  integer ( kind = 4 ) edge_data(5,6*tet_num)
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n1_old
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n2_old
  integer ( kind = 4 ) node
  real ( kind = 8 ) node_xyz1(3,node_num1)
  real ( kind = 8 ) node_xyz2(3,node_num2)
  integer ( kind = 4 ) tet
  integer ( kind = 4 ) tet_node1(4,tet_num)
  integer ( kind = 4 ) tet_node2(10,tet_num)
  integer ( kind = 4 ) v
  integer ( kind = 4 ) v1
  integer ( kind = 4 ) v2
!
!  Generate the index and coordinates of the new midside nodes, 
!  and update the tetradehron node data.
!
  node_xyz2(1:3,1:node_num1) = node_xyz1(1:3,1:node_num1)

  tet_node2(1:4,1:tet_num) = tet_node1(1:4,1:tet_num)

  node = node_num1

  n1_old = -1
  n2_old = -1

  do edge = 1, 6 * tet_num
!
!  Read the data defining the edge.
!
    n1 = edge_data(1,edge)
    n2 = edge_data(2,edge)
!
!  If this edge is new, create the coordinates and index.
!
    if ( n1 /= n1_old .or. n2 /= n2_old ) then

      node = node + 1

      if ( node_num2 < node ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TET_MESH_ORDER4_TO_ORDER10_COMPUTE - Fatal error!'
        write ( *, '(a)' ) '  Node index exceeds NODE_NUM2.'
        stop
      end if

      node_xyz2(1:3,node) = &
        ( node_xyz2(1:3,n1) + node_xyz2(1:3,n2) ) / 2.0D+00

      n1_old = n1
      n2_old = n2

    end if
!
!  Assign the node to the tetrahedron.
!
    v1 = edge_data(3,edge)
    v2 = edge_data(4,edge)
!
!  Here is where the local ordering of the nodes is effected:
!
    if ( v1 == 1 .and. v2 == 2 ) then
      v = 5
    else if ( v1 == 1 .and. v2 == 3 ) then
      v = 6
    else if ( v1 == 1 .and. v2 == 4 ) then
      v = 7
    else if ( v1 == 2 .and. v2 == 3 ) then
      v = 8
    else if ( v1 == 2 .and. v2 == 4 ) then
      v = 9
    else if ( v1 == 3 .and. v2 == 4 ) then
      v = 10
    end if

    tet = edge_data(5,edge)

    tet_node2(v,tet) = node

  end do

  return
end
subroutine tet_mesh_order4_to_order10_size ( tet_num, tet_node1, &
  node_num1, edge_data, node_num2 )

!*****************************************************************************80
!
!! TET_MESH_ORDER4_TO_ORDER10_SIZE sizes a quadratic tet mesh from a linear one.
!
!  Discussion:
!
!    A quadratic (10 node) tet mesh can be derived from a linear
!    (4 node) tet mesh by interpolating nodes at the midpoint of
!    every edge of the mesh.
!
!    The mesh is described indirectly, as the sum of individual
!    tetrahedrons.  A single physical edge may be a logical edge of
!    any number of tetrahedrons.  It is important, however, that a
!    new node be created exactly once for each edge, assigned an index,
!    and associated with every tetrahedron that shares this edge. 
!
!    This routine handles that problem.
!
!    The primary amount of work occurs in sorting a list of 6 * TET_NUM
!    data items, one item for every edge of every tetrahedron.  Each
!    data item records, for a given tetrahedron edge, the global indices
!    of the two endpoints, the local indices of the two endpoints,
!    and the index of the tetrahedron.
!
!    Through careful sorting, it is possible to arrange this data in
!    a way that allows the proper generation of the interpolated nodes.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) TET_NUM, the number of tetrahedrons in the
!    linear mesh.
!
!    Input, integer ( kind = 4 ) TET_NODE1(4,TET_NUM), the indices of 
!    the nodes in the linear mesh.
!
!    Input, integer ( kind = 4 ) NODE_NUM1, the number of nodes for the 
!    linear mesh.
!
!    Output, integer ( kind = 4 ) EDGE_DATA(5,6*TET_NUM), edge data.
!
!    Output, integer ( kind = 4 ) NODE_NUM2, the number of nodes for the
!    quadratic mesh.
!
  implicit none

  integer ( kind = 4 ) node_num1
  integer ( kind = 4 ) tet_num

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) edge
  integer ( kind = 4 ) edge_data(5,6*tet_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n1_old
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n2_old
  integer ( kind = 4 ) node_num2
  integer ( kind = 4 ) tet
  integer ( kind = 4 ) tet_node1(4,tet_num)
!
!  Step 1.
!  From the list of nodes for tetrahedron T, of the form: (I,J,K,L)
!  construct the six edge relations:
!
!    (I,J,1,2,T)
!    (I,K,1,3,T)
!    (I,L,1,4,T)
!    (J,K,2,3,T)
!    (J,L,2,4,T)
!    (K,L,3,4,T)
!
!  In order to make matching easier, we reorder each pair of nodes
!  into ascending order.
!
  do tet = 1, tet_num

    i = tet_node1(1,tet)
    j = tet_node1(2,tet)
    k = tet_node1(3,tet)
    l = tet_node1(4,tet)

    call i4i4_sort_a ( i, j, a, b )

    edge_data(1:5,6*(tet-1)+1) = (/ a, b, 1, 2, tet /)

    call i4i4_sort_a ( i, k, a, b )

    edge_data(1:5,6*(tet-1)+2) = (/ a, b, 1, 3, tet /)

    call i4i4_sort_a ( i, l, a, b )

    edge_data(1:5,6*(tet-1)+3) = (/ a, b, 1, 4, tet /)

    call i4i4_sort_a ( j, k, a, b )

    edge_data(1:5,6*(tet-1)+4) = (/ a, b, 2, 3, tet /)

    call i4i4_sort_a ( j, l, a, b )

    edge_data(1:5,6*(tet-1)+5) = (/ a, b, 2, 4, tet /)

    call i4i4_sort_a ( k, l, a, b )

    edge_data(1:5,6*(tet-1)+6) = (/ a, b, 3, 4, tet /)

  end do
!
!  Step 2. Perform an ascending dictionary sort on the neighbor relations.
!  We only intend to sort on rows 1:2; the routine we call here
!  sorts on the full column but that won't hurt us.
!
!  What we need is to find all cases where tetrahedrons share an edge.
!  By sorting the columns of the EDGE_DATA array, we will put shared edges
!  next to each other.
!
  call i4col_sort_a ( 5, 6*tet_num, edge_data )
!
!  Step 3. All the tetrahedrons which share an edge show up as consecutive
!  columns with identical first two entries.  Figure out how many new
!  nodes there are, and allocate space for their coordinates.
!
  node_num2 = node_num1

  n1_old = -1
  n2_old = -1

  do edge = 1, 6 * tet_num
    n1 = edge_data(1,edge)
    n2 = edge_data(2,edge)
    if ( n1 /= n1_old .or. n2 /= n2_old ) then
      node_num2 = node_num2 + 1
      n1_old = n1
      n2_old = n2
    end if
  end do

  return
end
subroutine tet_mesh_order10_example_set ( node_num, tet_num, &
  node_xyz, tet_node )

!*****************************************************************************80
!
!! TET_MESH_ORDER10_EXAMPLE_SET sets an example quadratic tet mesh.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 August 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, integer ( kind = 4 ) TET_NUM, the number of tetrahedrons.
!
!    Output, real ( kind = 8 ) NODE_XYZ(3,NODE_NUM), the node coordinates.
!
!    Output, integer ( kind = 4 ) TET_NODE(10,TET_NUM), the nodes 
!    forming each tet.
!
  implicit none

  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) tet_num

  integer ( kind = 4 ), dimension ( 10, tet_num ) :: tet_node
  real ( kind = 8 ), dimension ( 3, node_num ) :: node_xyz

  tet_node(1:10,1:tet_num) = reshape ( (/ &
    4,   3,   5,   1,  16,  19,  17,  11,  10,  12, &
    4,   2,   5,   1,  13,  19,  14,  11,   9,  12, &
    4,   7,   3,   5,  21,  16,  18,  19,  24,  17, &
    4,   7,   8,   5,  21,  22,  27,  19,  24,  25, &
    4,   6,   2,   5,  20,  13,  15,  19,  23,  14, &
    4,   6,   8,   5,  20,  22,  26,  19,  23,  25  &
  /), (/ 10, tet_num /) )

  node_xyz(1:3,1:node_num) = reshape ( (/ &
   0.0D+00,  0.0D+00,  0.0D+00, &   
   0.0D+00,  0.0D+00,  1.0D+00, &    
   0.0D+00,  1.0D+00,  0.0D+00, &    
   0.0D+00,  1.0D+00,  1.0D+00, &    
   1.0D+00,  0.0D+00,  0.0D+00, &  
   1.0D+00,  0.0D+00,  1.0D+00, &
   1.0D+00,  1.0D+00,  0.0D+00, &    
   1.0D+00,  1.0D+00,  1.0D+00, &    
   0.0D+00,  0.0D+00,  0.5D+00, &    
   0.0D+00,  0.5D+00,  0.0D+00, &   
   0.0D+00,  0.5D+00,  0.5D+00, &   
   0.5D+00,  0.0D+00,  0.0D+00, &   
   0.0D+00,  0.5D+00,  1.0D+00, &    
   0.5D+00,  0.0D+00,  0.5D+00, &    
   0.5D+00,  0.0D+00,  1.0D+00, &   
   0.0D+00,  1.0D+00,  0.5D+00, &    
   0.5D+00,  0.5D+00,  0.0D+00, &   
   0.5D+00,  1.0D+00,  0.0D+00, &   
   0.5D+00,  0.5D+00,  0.5D+00, &    
   0.5D+00,  0.5D+00,  1.0D+00, &   
   0.5D+00,  1.0D+00,  0.5D+00, &    
   0.5D+00,  1.0D+00,  1.0D+00, &    
   1.0D+00,  0.0D+00,  0.5D+00, &    
   1.0D+00,  0.5D+00,  0.0D+00, &    
   1.0D+00,  0.5D+00,  0.5D+00, &    
   1.0D+00,  0.5D+00,  1.0D+00, &   
   1.0D+00,  1.0D+00,  0.5D+00  &    
  /), (/ 3, node_num /) )

  return
end
subroutine tet_mesh_order10_example_size ( node_num, tet_num )

!*****************************************************************************80
!
!! TET_MESH_ORDER10_EXAMPLE_SIZE sizes an example quadratic tet mesh.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 August 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Output, integer ( kind = 4 ) TET_NUM, the number of tetrahedrons.
!
  implicit none

  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) tet_num

  node_num = 27
  tet_num = 6

  return
end
subroutine tet_mesh_order10_to_order4_compute ( tet_num1, tet_node1, &
  tet_num2, tet_node2 )

!*****************************************************************************80
!
!! TET_MESH_ORDER10_TO_ORDER4_COMPUTE linearizes a quadratic tet mesh.
!
!  Discussion:
!
!    A quadratic tet mesh is assumed to consist of 10-node
!    tetrahedrons.
!
!    This routine rearranges the information so as to define a 4-node
!    tet mesh.
!
!    The same nodes are used, but there are 8 times as many
!    tetrahedrons.
!
!    The node ordering for the quadratic tetrahedron is somewhat
!    arbitrary.  In the current scheme, the vertices are listed
!    first, followed by the 6 midside nodes.  Each midside node
!    may be identified by the two vertices that bracket it.  Thus,
!    the node ordering may be suggested by:
!
!      1  2  3  4 (1+2) (1+3) (1+4) (2+3) (2+4) (3+4)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Anwei Liu, Barry Joe,
!    Quality Local Refinement of Tetrahedral Meshes Based
!    on 8-Subtetrahedron Subdivision,
!    Mathematics of Computation,
!    Volume 65, Number 215, July 1996, pages 1183-1200.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) TET_NUM1, the number of tetrahedrons in 
!    the quadratic tet mesh.
!
!    Input, integer ( kind = 4 ) TET_NODE1(10,TET_NUM1), the indices of 
!    the nodes in the quadratic tet mesh.
!
!    Input, integer ( kind = 4 ) TET_NUM2, the number of tetrahedrons in 
!    the linear tet mesh.  TET_NUM2 = 8 * TET_NUM1.
!
!    Output, integer ( kind = 4 ) TET_NODE2(4,TET_NUM2), the indices of 
!    the nodes in the linear tet mesh.
!
  implicit none

  integer ( kind = 4 ) tet_num1
  integer ( kind = 4 ) tet_num2

  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n3
  integer ( kind = 4 ) n4
  integer ( kind = 4 ) n5
  integer ( kind = 4 ) n6
  integer ( kind = 4 ) n7
  integer ( kind = 4 ) n8
  integer ( kind = 4 ) n9
  integer ( kind = 4 ) nx
  integer ( kind = 4 ) tet1
  integer ( kind = 4 ) tet2
  integer ( kind = 4 ) tet_node1(10,tet_num1)
  integer ( kind = 4 ) tet_node2(4,tet_num2)

  tet2 = 0

  do tet1 = 1, tet_num1

    n1 = tet_node1(1,tet1)
    n2 = tet_node1(2,tet1)
    n3 = tet_node1(3,tet1)
    n4 = tet_node1(4,tet1)
    n5 = tet_node1(5,tet1)
    n6 = tet_node1(6,tet1)
    n7 = tet_node1(7,tet1)
    n8 = tet_node1(8,tet1)
    n9 = tet_node1(9,tet1)
    nx = tet_node1(10,tet1)

    tet2 = tet2 + 1
    tet_node2(1:4,tet2 ) = (/ n1, n5, n6, n7 /)
    tet2 = tet2 + 1
    tet_node2(1:4,tet2 ) = (/ n2, n5, n8, n9 /)
    tet2 = tet2 + 1
    tet_node2(1:4,tet2 ) = (/ n3, n6, n8, n9 /)
    tet2 = tet2 + 1
    tet_node2(1:4,tet2 ) = (/ n4, n7, n9, nx /)
    tet2 = tet2 + 1
    tet_node2(1:4,tet2)  = (/ n5, n6, n7, n9 /)
    tet2 = tet2 + 1
    tet_node2(1:4,tet2 ) = (/ n5, n6, n8, n9 /)
    tet2 = tet2 + 1
    tet_node2(1:4,tet2 ) = (/ n6, n7, n9, nx /)
    tet2 = tet2 + 1
    tet_node2(1:4,tet2 ) = (/ n6, n8, n9, nx /)

  end do

  return
end
subroutine tet_mesh_order10_to_order4_size ( node_num1, tet_num1, &
  node_num2, tet_num2 )

!*****************************************************************************80
!
!! TET_MESH_ORDER10_TO_ORDER4_SIZE sizes a linear tet mesh from a quadratic one.
!
!  Discussion:
!
!    A linear (4 node) tet mesh can be derived from a quadratic
!    (10 node) tet mesh using the same set of nodes, but reassigning
!    the nodes of each quadratic tet among 8 linear subtets.
!
!    This routine returns the number of nodes and tetrahedrons in the
!    linear mesh.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Anwei Liu, Barry Joe,
!    Quality Local Refinement of Tetrahedral Meshes Based
!    on 8-Subtetrahedron Subdivision,
!    Mathematics of Computation,
!    Volume 65, Number 215, July 1996, pages 1183-1200.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM1, the number of nodes in the 
!    quadratic mesh.
!
!    Input, integer ( kind = 4 ) TET_NUM1, the number of tetrahedrons in the
!    quadratic mesh.
!
!    Output, integer ( kind = 4 ) NODE_NUM2, the number of nodes for the 
!    linear mesh.
!
!    Output, integer ( kind = 4 ) TET_NUM2, the number of tetrahedrons in the
!    linear mesh.
!
  implicit none

  integer ( kind = 4 ) node_num1
  integer ( kind = 4 ) node_num2
  integer ( kind = 4 ) tet_num1
  integer ( kind = 4 ) tet_num2

  node_num2 = node_num1
  tet_num2 = 8 * tet_num1

  return
end
subroutine tet_mesh_quad ( node_num, node_xyz, tet_order, &
  tet_num, tet_node, f, quad_num, quad_xyz, quad_w, quad_value, &
  region_volume )

!*****************************************************************************80
!
!! TET_MESH_QUAD approximates an integral over a tet mesh.
!
!  Discussion:
!
!    The routine will accept tetrahedral meshes of order higher than 4.
!    However, only the first four nodes (the vertices) of each
!    tetrahedron will be used.  This will still produce correct results
!    for higher order tet meshes, as long as the sides of each
!    tetrahedron are flat (linear).
!
!    We assume that the vertices of each tetrahedron are listed first
!    in the description of higher order tetrahedrons.
!
!    The approximation of the integral is made using a quadrature rule 
!    defined on the unit tetrahedron, and supplied by the user.  
!
!    The user also supplies the name of a subroutine, here called "F", which
!    evaluates the integrand at a set of points.  The form of F is:
!
!      subroutine f ( n, xyz_vec, f_vec )
!      integer n
!      real ( kind = 8 ) f_vec(n)
!      real ( kind = 8 ) xyz_vec(3,n)
!
!    and it returns in each entry F_VEC(1:N), the value of the integrand
!    at XYZ_VEC(1:3,1:N).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 December 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes in the tet mesh.
!
!    Input, real ( kind = 8 ) NODE_XYZ(3,NODE_NUM), the coordinates 
!    of the nodes.
!
!    Input, integer ( kind = 4 ) TET_ORDER, the order of tetrahedrons in 
!    the tet mesh.
!
!    Input, integer ( kind = 4 ) TET_NUM, the number of tetrahedrons in the 
!    tet mesh.
!
!    Input, integer ( kind = 4 ) TET_NODE(TET_ORDER,TET_NUM), the indices
!    of the nodes.
!
!    Input, external F, the name of the subroutine that evaluates the integrand.
!
!    Input, integer ( kind = 4 ) QUAD_NUM, the order of the quadrature rule.
!
!    Input, real ( kind = 8 ) QUAD_XYZ(3,QUAD_NUM), the abscissas of the 
!    quadrature rule, in the unit tetrahedron.
!
!    Input, real ( kind = 8 ) QUAD_W(QUAD_NUM), the weights of the 
!    quadrature rule.
!
!    Output, real ( kind = 8 ) QUAD_VALUE, the estimate of the integral
!    of F(X,Y) over the region covered by the tet mesh.
!
!    Output, real ( kind = 8 ) REGION_VOLUME, the volume of the region.
!
  implicit none

  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) quad_num
  integer ( kind = 4 ) tet_num
  integer ( kind = 4 ) tet_order

  external f
  real ( kind = 8 ) node_xyz(3,node_num)
  real ( kind = 8 ) quad_f(quad_num)
  real ( kind = 8 ) quad_value
  real ( kind = 8 ) quad_w(quad_num)
  real ( kind = 8 ) quad_xyz(3,quad_num)
  real ( kind = 8 ) quad2_xyz(3,quad_num)
  real ( kind = 8 ) region_volume
  integer ( kind = 4 ) tet
  real ( kind = 8 ) tet_volume
  integer ( kind = 4 ) tet_node(tet_order,tet_num)
  real ( kind = 8 ) tet_xyz(3,4)

  quad_value = 0.0D+00
  region_volume = 0.0D+00

  do tet = 1, tet_num

    tet_xyz(1:3,1:4) = node_xyz(1:3,tet_node(1:4,tet))

    call tetrahedron_volume ( tet_xyz, tet_volume )

    call tetrahedron_order4_reference_to_physical ( tet_xyz, &
      quad_num, quad_xyz, quad2_xyz )

    call f ( quad_num, quad2_xyz, quad_f )

    quad_value = quad_value + tet_volume &
      * dot_product ( quad_w(1:quad_num), quad_f(1:quad_num) )

    region_volume = region_volume + tet_volume

  end do

  return
end
subroutine tet_mesh_quality1 ( node_num, node_xyz, tet_order, tet_num, &
  tet_node, tet_quality )

!*****************************************************************************80
!
!! TET_MESH_QUALITY1 returns the quality of each tet in a mesh.
!
!  Discussion:
!
!    The overall tet mesh quality measure is the minimum of the corresponding
!    tetrahedron quality measure, over all tetrahedrons in the tet mesh.
!
!    Although tetrahedrons of order 10 are allowed as input,
!    only the first 4 nodes (presumed to be the vertices) are used
!    in the calculation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 November 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, real ( kind = 8 ) NODE_XYZ(3,NODE_NUM), the nodes.
!
!    Input, integer ( kind = 4 ) TET_ORDER, the order of the mesh, either
!    4 or 10.
!
!    Input, integer ( kind = 4 ) TET_NUM, the number of tetrahedrons.
!
!    Input, integer ( kind = 4 ) TET_NODE(TET_ORDER,TET_NUM), the 
!    indices of the nodes.
!
!    Output, real ( kind = 8 ) TET_QUALITY(TET_NUM), the quality
!    measure for each tetrahedron.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3
  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) tet_num
  integer ( kind = 4 ) tet_order

  real ( kind = 8 ) node_xyz(dim_num,node_num)
  integer ( kind = 4 ) tet
  integer ( kind = 4 ) tet_node(tet_order,tet_num)
  real ( kind = 8 ) tet_quality(tet_num)
  real ( kind = 8 ) tet_xyz(dim_num,4)

  do tet = 1, tet_num

    tet_xyz(1:dim_num,1:4) = node_xyz(1:dim_num,tet_node(1:4,tet))

    call tetrahedron_quality1_3d ( tet_xyz, tet_quality(tet) )

  end do

  return
end
subroutine tet_mesh_quality2 ( node_num, node_xyz, tet_order, tet_num, &
  tet_node, tet_quality )

!*****************************************************************************80
!
!! TET_MESH_QUALITY2 returns the quality of each tet in a mesh.
!
!  Discussion:
!
!    The overall tet mesh quality measure is the minimum of the
!    corresponding tetrahedron quality measure, over all tetrahedrons in the
!    tet mesh.
!
!    Although tetrahedrons of order 10 are allowed as input,
!    only the first 4 nodes (presumed to be the vertices) are used
!    in the calculation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 November 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, real ( kind = 8 ) NODE_XYZ(3,NODE_NUM), the nodes.
!
!    Input, integer ( kind = 4 ) TET_ORDER, the order of the mesh, either 
!    4 or 10.
!
!    Input, integer ( kind = 4 ) TET_NUM, the number of tetrahedrons.
!
!    Input, integer ( kind = 4 ) TET_NODE(TET_ORDER,TET_NUM), the 
!    indices of the nodes.
!
!    Output, real ( kind = 8 ) TET_QUALITY(TET_NUM), the quality
!    measure for each tetrahedron.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3
  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) tet_num
  integer ( kind = 4 ) tet_order

  real ( kind = 8 ) node_xyz(dim_num,node_num)
  integer ( kind = 4 ) tet
  integer ( kind = 4 ) tet_node(tet_order,tet_num)
  real ( kind = 8 ) tet_quality(tet_num)
  real ( kind = 8 ) tet_xyz(dim_num,4)

  do tet = 1, tet_num

    tet_xyz(1:dim_num,1:4) = node_xyz(1:dim_num,tet_node(1:4,tet))

    call tetrahedron_quality2_3d ( tet_xyz, tet_quality(tet) )

  end do

  return
end
subroutine tet_mesh_quality3 ( node_num, node_xyz, tet_order, tet_num, &
  tet_node, tet_quality )

!*****************************************************************************80
!
!! TET_MESH_QUALITY3 returns the quality of each tet in a mesh.
!
!  Discussion:
!
!    The overall tet mesh quality measure is the minimum of the
!    corresponding tetrahedron quality measure, over all tetrahedrons in the
!    tet mesh.
!
!    Although tetrahedrons of order 10 are allowed as input,
!    only the first 4 nodes (presumed to be the vertices) are used
!    in the calculation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 November 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, real ( kind = 8 ) NODE_XYZ(3,NODE_NUM), the nodes.
!
!    Input, integer ( kind = 4 ) TET_ORDER, the order of the mesh, either 
!    4 or 10.
!
!    Input, integer ( kind = 4 ) TET_NUM, the number of tetrahedrons.
!
!    Input, integer ( kind = 4 ) TET_NODE(TET_ORDER,TET_NUM), 
!    the indices of the nodes.
!
!    Output, real ( kind = 8 ) TET_QUALITY(TET_NUM), the quality
!    measure for each tetrahedron.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3
  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) tet_num
  integer ( kind = 4 ) tet_order

  real ( kind = 8 ) node_xyz(dim_num,node_num)
  integer ( kind = 4 ) tet
  integer ( kind = 4 ) tet_node(tet_order,tet_num)
  real ( kind = 8 ) tet_quality(tet_num)
  real ( kind = 8 ) tet_xyz(dim_num,4)

  do tet = 1, tet_num

    tet_xyz(1:dim_num,1:4) = node_xyz(1:dim_num,tet_node(1:4,tet))

    call tetrahedron_quality3_3d ( tet_xyz, tet_quality(tet) )

  end do

  return
end
subroutine tet_mesh_quality4 ( node_num, node_xyz, tet_order, tet_num, &
  tet_node, tet_quality )

!*****************************************************************************80
!
!! TET_MESH_QUALITY4 returns the quality of each tet in a mesh.
!
!  Discussion:
!
!    The overall tet mesh quality measure is the minimum of the
!    corresponding tetrahedron quality measure, over all tetrahedrons in the
!    tet mesh.
!
!    Although tetrahedrons of order 10 are allowed as input,
!    only the first 4 nodes (presumed to be the vertices) are used
!    in the calculation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 November 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, real ( kind = 8 ) NODE_XYZ(3,NODE_NUM), the nodes.
!
!    Input, integer ( kind = 4 ) TET_ORDER, the order of the mesh, either 
!    4 or 10.
!
!    Input, integer ( kind = 4 ) TET_NUM, the number of tetrahedrons.
!
!    Input, integer ( kind = 4 ) TET_NODE(TET_ORDER,TET_NUM), the 
!    indices of the nodes.
!
!    Output, real ( kind = 8 ) TET_QUALITY(TET_NUM), the quality
!    measure for each tetrahedron.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3
  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) tet_num
  integer ( kind = 4 ) tet_order

  real ( kind = 8 ) node_xyz(dim_num,node_num)
  integer ( kind = 4 ) tet
  integer ( kind = 4 ) tet_node(tet_order,tet_num)
  real ( kind = 8 ) tet_quality(tet_num)
  real ( kind = 8 ) tet_xyz(dim_num,4)

  do tet = 1, tet_num

    tet_xyz(1:dim_num,1:4) = node_xyz(1:dim_num,tet_node(1:4,tet))

    call tetrahedron_quality4_3d ( tet_xyz, tet_quality(tet) )

  end do

  return
end
subroutine tet_mesh_quality5 ( node_num, node_xyz, tet_order, tet_num, &
  tet_node, tet_quality )

!*****************************************************************************80
!
!! TET_MESH_QUALITY5 returns the quality of each tet in a mesh.
!
!  Discussion:
!
!    The overall tet mesh quality measure is the ratio of the minimum
!    tetrahedron volume to the maximum tetrahedron volume.
!
!    Although tetrahedrons of order 10 are allowed as input,
!    only the first 4 nodes (presumed to be the vertices) are used
!    in the calculation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 November 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, real ( kind = 8 ) NODE_XYZ(3,NODE_NUM), the nodes.
!
!    Input, integer ( kind = 4 ) TET_ORDER, the order of the mesh, 
!    either 4 or 10.
!
!    Input, integer ( kind = 4 ) TET_NUM, the number of tetrahedrons.
!
!    Input, integer ( kind = 4 ) TET_NODE(TET_ORDER,TET_NUM), 
!    the indices of the nodes.
!
!    Output, real ( kind = 8 ) TET_QUALITY(TET_NUM), the quality
!    measure for each tetrahedron.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3
  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) tet_num
  integer ( kind = 4 ) tet_order

  real ( kind = 8 ) node_xyz(dim_num,node_num)
  integer ( kind = 4 ) tet
  integer ( kind = 4 ) tet_node(tet_order,tet_num)
  real ( kind = 8 ) tet_quality(tet_num)
  real ( kind = 8 ) tet_xyz(dim_num,4)
  real ( kind = 8 ) volume_max

  do tet = 1, tet_num

    tet_xyz(1:dim_num,1:4) = node_xyz(1:dim_num,tet_node(1:4,tet))

    call tetrahedron_volume ( tet_xyz, tet_quality(tet) )

  end do

  volume_max = maxval ( tet_quality(1:tet_num) )

  if ( 0.0D+00 < volume_max ) then
    tet_quality(1:tet_num) = tet_quality(1:tet_num) / volume_max
  end if

  return
end
subroutine tet_mesh_search_delaunay ( node_num, node_xyz, tet_order, &
  tet_num, tet_node, tet_neighbor, p, tet_index, face, step_num )

!*****************************************************************************80
!
!! TET_MESH_SEARCH_DELAUNAY searches a Delaunay tet mesh for a point.
!
!  Discussion:
!
!    The algorithm "walks" from one tetrahedron to its neighboring tetrahedron,
!    and so on, until a tetrahedron is found containing point P, or P is found
!    to be outside the convex hull.
!
!    The algorithm computes the barycentric coordinates of the point with
!    respect to the current tetrahedron.  If all 4 quantities are positive,
!    the point is contained in the tetrahedron.  If the I-th coordinate is
!    negative, then P lies on the far side of edge I, which is opposite
!    from vertex I.  This gives a hint as to where to search next.
!
!    For a Delaunay tet mesh, the search is guaranteed to terminate.
!    For other meshes, a cycle may occur.
!
!    Note the surprising fact that, even for a Delaunay tet mesh of
!    a set of nodes, the nearest node to P need not be one of the
!    vertices of the tetrahedron containing P.
!
!    The code can be called for tet meshes of any order, but only
!    the first 4 nodes in each tetrahedron are considered.  Thus, if
!    higher order tetrahedrons are used, and the extra nodes are intended
!    to give the tetrahedron a polygonal shape, these will have no effect,
!    and the results obtained here might be misleading.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 August 2009
!
!  Author:
!
!    John Burkardt.
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
!    Input, real ( kind = 8 ) NODE_XYZ(3,NODE_NUM), the coordinates of 
!    the nodes.
!
!    Input, integer ( kind = 4 ) TET_ORDER, the order of the tetrahedrons.
!
!    Input, integer ( kind = 4 ) TET_NUM, the number of tetrahedrons.
!
!    Input, integer ( kind = 4 ) TET_NODE(TET_ORDER,TET_NUM),
!    the nodes that make up each tetrahedron.
!
!    Input, integer ( kind = 4 ) TET_NEIGHBOR(4,TET_NUM), the 
!    tetrahedron neighbor list.
!
!    Input, real ( kind = 8 ) P(3), the coordinates of a point.
!
!    Output, integer ( kind = 4 ) TET_INDEX, the index of the tetrahedron 
!    where the search ended.  If a cycle occurred, then TET_INDEX = -1.
!
!    Output, integer ( kind = 4 ) FACE, indicates the position of the point P in
!    face TET_INDEX:
!    0, the interior or boundary of the tetrahedron;
!    -1, outside the convex hull of the tet mesh, past face 1;
!    -2, outside the convex hull of the tet mesh, past face 2;
!    -3, outside the convex hull of the tet mesh, past face 3.
!    -4, outside the convex hull of the tet mesh, past face 4.
!
!    Output, integer ( kind = 4 ) STEP_NUM, the number of steps taken.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3
  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) tet_num
  integer ( kind = 4 ) tet_order

  real ( kind = 8 ) alpha(dim_num+1)
  integer ( kind = 4 ) face
  real ( kind = 8 ) node_xyz(dim_num,node_num)
  real ( kind = 8 ) p(dim_num)
  integer ( kind = 4 ) step_num
  integer ( kind = 4 ) tet_node(tet_order,tet_num)
  integer ( kind = 4 ) tet_index
  integer ( kind = 4 ), save :: tet_index_save = -1
  integer ( kind = 4 ) tet_neighbor(dim_num+1,tet_num)
!
!  If possible, start with the previous successful value of TET_INDEX.
!
  if ( tet_index_save < 1 .or. tet_num < tet_index_save ) then
    tet_index = ( tet_num + 1 ) / 2
  else
    tet_index = tet_index_save
  end if

  step_num = -1
  face = 0

  do

    step_num = step_num + 1

    if ( tet_num < step_num ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TET_MESH_SEARCH_DELAUNAY - Fatal error!'
      write ( *, '(a)' ) '  The algorithm seems to be cycling.'
      tet_index = -1
      face = -1
      stop
    end if

    call tetrahedron_barycentric ( node_xyz(1:3,tet_node(1:4,tet_index)), &
      p(1:3), alpha )
!
!  If the barycentric coordinates are all positive, then the point
!  is inside the tetrahedron and we're done.
!
    if ( 0.0D+00 <= alpha(1) .and. &
         0.0D+00 <= alpha(2) .and. &
         0.0D+00 <= alpha(3) .and. &
         0.0D+00 <= alpha(4) ) then
      exit
    end if
!
!  At least one barycentric coordinate is negative.
!
!  If there is a negative barycentric coordinate for which there exists an
!  opposing tetrahedron neighbor closer to the point, move to that tetrahedron.
!
    if ( alpha(1) < 0.0D+00 .and. 0 < tet_neighbor(1,tet_index) ) then
      tet_index = tet_neighbor(1,tet_index)
      cycle
    else if ( alpha(2) < 0.0D+00 .and. &
      0 < tet_neighbor(2,tet_index) ) then
      tet_index = tet_neighbor(2,tet_index)
      cycle
    else if ( alpha(3) < 0.0D+00 .and. &
      0 < tet_neighbor(3,tet_index) ) then
      tet_index = tet_neighbor(3,tet_index)
      cycle
    else if ( alpha(4) < 0.0D+00 .and. &
      0 < tet_neighbor(4,tet_index) ) then
      tet_index = tet_neighbor(4,tet_index)
      cycle
    end if
!
!  All negative barycentric coordinates correspond to vertices opposite
!  faces on the convex hull.
!
!  Note the face and exit.
!
    if ( alpha(1) < 0.0D+00 ) then
      face = -1
      exit
    else if ( alpha(2) < 0.0D+00 ) then
      face = -2
      exit
    else if ( alpha(3) < 0.0D+00 ) then
      face = -3
      exit
    else if ( alpha(4) < 0.0D+00 ) then
      face = -4
      exit
    end if

  end do

  tet_index_save = tet_index

  return
end
subroutine tet_mesh_search_naive ( node_num, node_xyz, &
  tet_order, tet_num, tet_node, p, tet_index, step_num )

!*****************************************************************************80
!
!! TET_MESH_SEARCH_NAIVE naively searches a tet mesh.
!
!  Discussion:
!
!    The algorithm simply checks each tetrahedron to see if point P is
!    contained in it.  
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 August 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, real ( kind = 8 ) NODE_XYZ(3,NODE_NUM), the coordinates 
!    of the nodes.
!
!    Input, integer ( kind = 4 ) TET_ORDER, the order of the tetrahedrons.
!
!    Input, integer ( kind = 4 ) TET_NUM, the number of tetrahedrons in
!    the mesh.
!
!    Input, integer ( kind = 4 ) TET_NODE(TET_ORDER,TET_NUM), 
!    the nodes that make up each tetrahedron.
!
!    Input, real ( kind = 8 ) P(3), the coordinates of a point.
!
!    Output, integer ( kind = 4 ) TET_INDEX, the index of the tetrahedron
!    where the search ended, or -1 if no tetrahedron was found containing
!    the point.
!
!    Output, integer ( kind = 4 ) STEP_NUM, the number of tetrahedrons checked.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3
  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) tet_num
  integer ( kind = 4 ) tet_order

  real ( kind = 8 ) alpha(4)
  real ( kind = 8 ) node_xyz(dim_num,node_num)
  real ( kind = 8 ) p(dim_num)
  integer ( kind = 4 ) step_num
  integer ( kind = 4 ) tet
  integer ( kind = 4 ) tet_node(tet_order,tet_num)
  integer ( kind = 4 ) tet_index

  tet_index = -1
  step_num = 0

  do tet = 1, tet_num

    call tetrahedron_barycentric ( node_xyz(1:3,tet_node(1:4,tet)), &
      p(1:dim_num), alpha )

    if ( all ( 0 <= alpha(1:4) ) ) then
      tet_index = tet
      step_num = tet
      return
    end if

  end do

  return
end
subroutine tetrahedron_barycentric ( tetra, p, c )

!*****************************************************************************80
!
!! TETRAHEDRON_BARYCENTRIC: barycentric coordinates of a point.
!
!  Discussion:
!
!    The barycentric coordinates of a point P with respect to
!    a tetrahedron are a set of four values C(1:4), each associated
!    with a vertex of the tetrahedron.  The values must sum to 1.
!    If all the values are between 0 and 1, the point is contained
!    within the tetrahedron.
!
!    The barycentric coordinate of point P related to vertex A can be
!    interpreted as the ratio of the volume of the tetrahedron with 
!    vertex A replaced by vertex P to the volume of the original 
!    tetrahedron.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) TETRA(3,4) the tetrahedron vertices.
!
!    Input, real ( kind = 8 ) P(3), the point to be checked.
!
!    Output, real ( kind = 8 ) C(4), the barycentric coordinates of P with
!    respect to the tetrahedron.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3
  integer ( kind = 4 ), parameter :: rhs_num = 1

  real ( kind = 8 ) a(dim_num,dim_num+rhs_num)
  real ( kind = 8 ) c(dim_num+1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  real ( kind = 8 ) p(dim_num)
  real ( kind = 8 ) tetra(dim_num,4)
!
!  Set up the linear system
!
!    ( X2-X1  X3-X1  X4-X1 ) C2    X - X1
!    ( Y2-Y1  Y3-Y1  Y4-Y1 ) C3  = Y - Y1
!    ( Z2-Z1  Z3-Z1  Z4-Z1 ) C4    Z - Z1
!
!  which is satisfied by the barycentric coordinates of P.
!
  a(1:dim_num,1:3) = tetra(1:dim_num,2:4)
  a(1:dim_num,4) = p(1:dim_num)

  do i = 1, dim_num
    a(i,1:4) = a(i,1:4) - tetra(i,1)
  end do
!
!  Solve the linear system.
!
  call r8mat_solve ( dim_num, rhs_num, a, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TETRAHEDRON_BARYCENTRIC - Fatal error!'
    write ( *, '(a)' ) '  The linear system is singular.'
    write ( *, '(a)' ) '  The input data does not form a proper tetrahedron.'
    stop
  end if

  c(2:4) = a(1:dim_num,4)

  c(1) = 1.0D+00 - sum ( c(2:4) )

  return
end
subroutine tetrahedron_circumsphere_3d ( tet_xyz, r, pc )

!*****************************************************************************80
!
!! TETRAHEDRON_CIRCUMSPHERE_3D computes the circumsphere of a tetrahedron in 3D.
!
!  Discussion:
!
!    The circumsphere, or circumscribed sphere, of a tetrahedron is the 
!    sphere that passes through the four vertices.  The circumsphere is
!    not necessarily the smallest sphere that contains the tetrahedron.
!
!    Surprisingly, the diameter of the sphere can be found by solving
!    a 3 by 3 linear system.  This is because the vectors P2 - P1,
!    P3 - P1 and P4 - P1 are secants of the sphere, and each forms a
!    right triangle with the diameter through P1.  Hence, the dot product of
!    P2 - P1 with that diameter is equal to the square of the length
!    of P2 - P1, and similarly for P3 - P1 and P4 - P1.  This determines
!    the diameter vector originating at P1, and hence the radius and
!    center.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Adrian Bowyer, John Woodwark,
!    A Programmer's Geometry,
!    Butterworths, 1983,
!    ISBN: 0408012420.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) TET_XYZ(3,4) the coordinates of the vertices.
!
!    Output, real ( kind = 8 ) R, PC(3), the center of the
!    circumscribed sphere, and its radius.  If the linear system is
!    singular, then R = -1, PC(1:3) = 0.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3
  integer ( kind = 4 ), parameter :: rhs_num = 1

  real ( kind = 8 ) a(dim_num,dim_num+rhs_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) j
  real ( kind = 8 ) pc(dim_num)
  real ( kind = 8 ) r
  real ( kind = 8 ) tet_xyz(dim_num,4)
!
!  Set up the linear system.
!
  a(1:dim_num,1:3) = transpose ( tet_xyz(1:dim_num,2:4) )

  do j = 1, dim_num
    a(1:dim_num,j) = a(1:dim_num,j) - tet_xyz(j,1)
  end do

  do i = 1, 3
    a(i,4) = sum ( a(i,1:3)**2 )
  end do
!
!  Solve the linear system.
!
  call r8mat_solve ( dim_num, rhs_num, a, info )
!
!  If the system was singular, return a consolation prize.
!
  if ( info /= 0 ) then
    r = -1.0D+00
    pc(1:dim_num) = 0.0D+00
    return
  end if
!
!  Compute the radius and center.
!
  r = 0.5D+00 * sqrt ( sum ( a(1:dim_num,4)**2 ) )

  pc(1:dim_num) = tet_xyz(1:dim_num,1) + 0.5D+00 * a(1:dim_num,4)

  return
end
subroutine tetrahedron_edge_length_3d ( tet_xyz, edge_length )

!*****************************************************************************80
!
!! TETRAHEDRON_EDGE_LENGTH_3D returns edge lengths of a tetrahedron in 3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) TET_XYZ(3,4), the coordinates of the vertices.
!
!    Output, real ( kind = 8 ) EDGE_LENGTH(6), the length of the edges.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  real ( kind = 8 ) r8vec_length
  real ( kind = 8 ) edge_length(6)
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) k
  real ( kind = 8 ) tet_xyz(dim_num,4)

  k = 0
  do j1 = 1, 3
    do j2 = j1+1, 4
      k = k + 1
      edge_length(k) = r8vec_length ( dim_num, &
        tet_xyz(1:dim_num,j2) - tet_xyz(1:dim_num,j1) )
    end do
  end do

  return
end
subroutine tetrahedron_insphere_3d ( tet_xyz, r, pc )

!*****************************************************************************80
!
!! TETRAHEDRON_INSPHERE_3D finds the insphere of a tetrahedron in 3D.
!
!  Discussion:
!
!    The insphere of a tetrahedron is the inscribed sphere, which touches 
!    each face of the tetrahedron at a single point.
!
!    The points of contact are the centroids of the triangular faces
!    of the tetrahedron.  Therefore, the point of contact for a face
!    can be computed as the average of the vertices of that face.
!
!    The sphere can then be determined as the unique sphere through
!    the four given centroids.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Philip Schneider, David Eberly,
!    Geometric Tools for Computer Graphics,
!    Elsevier, 2002,
!    ISBN: 1558605940,
!    LC: T385.G6974.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) TET_XYZ(3,4), the coordinates of the vertices.
!
!    Output, real ( kind = 8 ) R, PC(3), the radius and the center
!    of the sphere.  
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  real ( kind = 8 ) b(4,4)
  real ( kind = 8 ) gamma
  real ( kind = 8 ) l123
  real ( kind = 8 ) l124
  real ( kind = 8 ) l134
  real ( kind = 8 ) l234
  real ( kind = 8 ) n123(1:dim_num)
  real ( kind = 8 ) n124(1:dim_num)
  real ( kind = 8 ) n134(1:dim_num)
  real ( kind = 8 ) n234(1:dim_num)
  real ( kind = 8 ) pc(1:dim_num)
  real ( kind = 8 ) r
  real ( kind = 8 ) r8mat_det_4d
  real ( kind = 8 ) r8vec_length
  real ( kind = 8 ) tet_xyz( 1:dim_num,4)
  real ( kind = 8 ) v21(1:dim_num)
  real ( kind = 8 ) v31(1:dim_num)
  real ( kind = 8 ) v41(1:dim_num)
  real ( kind = 8 ) v32(1:dim_num)
  real ( kind = 8 ) v42(1:dim_num) 
  real ( kind = 8 ) v43(1:dim_num) 
  
  v21(1:dim_num) = tet_xyz( 1:dim_num,2) - tet_xyz( 1:dim_num,1)
  v31(1:dim_num) = tet_xyz( 1:dim_num,3) - tet_xyz( 1:dim_num,1)
  v41(1:dim_num) = tet_xyz( 1:dim_num,4) - tet_xyz( 1:dim_num,1)
  v32(1:dim_num) = tet_xyz( 1:dim_num,3) - tet_xyz( 1:dim_num,2)
  v42(1:dim_num) = tet_xyz( 1:dim_num,4) - tet_xyz( 1:dim_num,2)
  v43(1:dim_num) = tet_xyz( 1:dim_num,4) - tet_xyz( 1:dim_num,3)

  call r8vec_cross_3d ( v21, v31, n123 )
  call r8vec_cross_3d ( v41, v21, n124 )
  call r8vec_cross_3d ( v31, v41, n134 )
  call r8vec_cross_3d ( v42, v32, n234 )

  l123 = r8vec_length ( dim_num, n123 )
  l124 = r8vec_length ( dim_num, n124 )
  l134 = r8vec_length ( dim_num, n134 )
  l234 = r8vec_length ( dim_num, n234 )

  pc(1:dim_num) = ( l234 * tet_xyz( 1:dim_num,1)   &
                  + l134 * tet_xyz( 1:dim_num,2)   &
                  + l124 * tet_xyz( 1:dim_num,3)   &
                  + l123 * tet_xyz( 1:dim_num,4) ) &
                / ( l234 + l134 + l124 + l123 )

  b(1:dim_num,1:4) = tet_xyz( 1:dim_num,1:4)
  b(4,1:4) = 1.0D+00

  gamma = abs ( r8mat_det_4d ( b ) )

! gamma = abs ( &
!     ( tet_xyz( 1,2) * tet_xyz( 2,3) * tet_xyz( 3,4) &
!     - tet_xyz( 1,3) * tet_xyz( 2,4) * tet_xyz( 3,2) &
!     + tet_xyz( 1,4) * tet_xyz( 2,2) * tet_xyz( 3,3) ) &
!   - ( tet_xyz( 1,1) * tet_xyz( 2,3) * tet_xyz( 3,4) &
!     - tet_xyz( 1,3) * tet_xyz( 2,4) * tet_xyz( 3,1) &
!     + tet_xyz( 1,4) * tet_xyz( 2,1) * tet_xyz( 3,3) ) &
!   + ( tet_xyz( 1,1) * tet_xyz( 2,2) * tet_xyz( 3,4) &
!     - tet_xyz( 1,2) * tet_xyz( 2,4) * tet_xyz( 3,1) & 
!     + tet_xyz( 1,4) * tet_xyz( 2,1) * tet_xyz( 3,2) ) &
!   - ( tet_xyz( 1,1) * tet_xyz( 2,2) * tet_xyz( 3,3) &
!     - tet_xyz( 1,2) * tet_xyz( 2,3) * tet_xyz( 3,1) &
!     + tet_xyz( 1,3) * tet_xyz( 2,1) * tet_xyz( 3,2) ) )
 
  r = gamma / ( l234 + l134 + l124 + l123 )

  return
end
subroutine tetrahedron_order4_physical_to_reference ( tet_xyz, n, phy, ref )

!*****************************************************************************80
!
!! TETRAHEDRON_ORDER4_PHYSICAL_TO_REFERENCE: physical to reference points.
!
!  Discussion:
!
!    Given the vertices of an order 4 physical tetrahedron and a point
!    (X,Y,Z) in the physical tetrahedron, the routine computes the value
!    of the corresponding image point (R,S,T) in reference space.
!
!    This routine may be appropriate for an order 10 tetrahedron,
!    if the mapping between reference and physical space is linear.  
!    This implies, in particular, that the edges of the image tetrahedron 
!    are straight, the faces are flat, and the "midside" nodes in the
!    physical tetrahedron are halfway along the sides of the physical 
!    tetrahedron.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) TET_XYZ(3,4), the coordinates of the vertices.  
!    The vertices are assumed to be the images of
!    (0,0,0), (1,0,0), (0,1,0) and (0,0,1) respectively.
!
!    Input, integer ( kind = 4 ) N, the number of points to transform.
!
!    Input, real ( kind = 8 ) PHY(3,N), the coordinates of physical points
!    to be transformed.
!
!    Output, real ( kind = 8 ) REF(3,N), the coordinates of the corresponding
!    points in the reference space.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(3,3)
  real ( kind = 8 ) det
  real ( kind = 8 ) phy(3,n)
  real ( kind = 8 ) ref(3,n)
  real ( kind = 8 ) tet_xyz( 3,4)
!
!  Set up the matrix.
!
  a(1:3,1) = tet_xyz( 1:3,2) - tet_xyz( 1:3,1)
  a(1:3,2) = tet_xyz( 1:3,3) - tet_xyz( 1:3,1)
  a(1:3,3) = tet_xyz( 1:3,4) - tet_xyz( 1:3,1)
!
!  Compute the determinant.
!
  det =  a(1,1) * ( a(2,2) * a(3,3) - a(2,3) * a(3,2) ) &
       + a(1,2) * ( a(2,3) * a(3,1) - a(2,1) * a(3,3) ) &
       + a(1,3) * ( a(2,1) * a(3,2) - a(2,2) * a(3,1) )
!
!  If the determinant is zero, bail out.
!
  if ( det == 0.0D+00 ) then
    ref(1:3,1:n) = 0.0D+00
    return
  end if
!
!  Compute the solution.
!
  ref(1,1:n) = (   ( a(2,2) * a(3,3) - a(2,3) * a(3,2) ) &
                 * ( phy(1,1:n) - tet_xyz( 1,1) ) &
                 - ( a(1,2) * a(3,3) - a(1,3) * a(3,2) ) &
                 * ( phy(2,1:n) - tet_xyz( 2,1) ) &
                 + ( a(1,2) * a(2,3) - a(1,3) * a(2,2) ) &
                 * ( phy(3,1:n) - tet_xyz( 3,1) ) &
               ) / det

  ref(2,1:n) = ( - ( a(2,1) * a(3,3) - a(2,3) * a(3,1) ) &
                 * ( phy(1,1:n) - tet_xyz( 1,1) ) &
                 + ( a(1,1) * a(3,3) - a(1,3) * a(3,1) ) &
                 * ( phy(2,1:n) - tet_xyz( 2,1) ) &
                 - ( a(1,1) * a(2,3) - a(1,3) * a(2,1) ) &
                 * ( phy(3,1:n) - tet_xyz( 3,1) ) &
               ) / det

  ref(3,1:n) = (   ( a(2,1) * a(3,2) - a(2,2) * a(3,1) ) &
                 * ( phy(1,1:n) - tet_xyz( 1,1) ) &
                 - ( a(1,1) * a(3,2) - a(1,2) * a(3,1) ) &
                 * ( phy(2,1:n) - tet_xyz( 2,1) ) &
                 + ( a(1,1) * a(2,2) - a(1,2) * a(2,1) ) &
                 * ( phy(3,1:n) - tet_xyz( 3,1) ) &
               ) / det

  return
end
subroutine tetrahedron_order4_reference_to_physical ( tet_xyz, n, ref, phy )

!*****************************************************************************80
!
!! TETRAHEDRON_ORDER4_REFERENCE_TO_PHYSICAL: T4 reference to physical points.
!
!  Discussion:
!
!    Given the vertices of an order 4 physical tetrahedron and a point
!    (R,S,T) in the reference tetrahedron, the routine computes the value
!    of the corresponding image point (X,Y,Z) in physical space.
!
!    This routine will also be correct for an order 10 tetrahedron,
!    if the mapping between reference and physical space
!    is linear.  This implies, in particular, that the sides of the
!    image tetrahedron are straight, the faces are flat, and 
!    the "midside" nodes in the physical tetrahedron are
!    halfway along the edges of the physical tetrahedron.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 December 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) TET_XYZ(3,4), the coordinates of the vertices.
!    The vertices are assumed to be the images of (0,0,0), (1,0,0),
!    (0,1,0) and (0,0,1) respectively.
!
!    Input, integer ( kind = 4 ) N, the number of points to transform.
!
!    Input, real ( kind = 8 ) REF(3,N), points in the reference tetrahedron.
!
!    Output, real ( kind = 8 ) PHY(3,N), corresponding points in the
!    physical tetrahedron.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real ( kind = 8 ) phy(3,n)
  real ( kind = 8 ) ref(3,n)
  real ( kind = 8 ) tet_xyz( 3,4)

  do i = 1, 3
    phy(i,1:n) = &
        tet_xyz( i,1) * ( 1.0D+00 - ref(1,1:n) - ref(2,1:n) - ref(3,1:n) ) &
      + tet_xyz( i,2) *             ref(1,1:n)                             &
      + tet_xyz( i,3) *                          ref(2,1:n)                &
      + tet_xyz( i,4) *                                       ref(3,1:n)
  end do

  return
end
subroutine tetrahedron_quality1_3d ( tet_xyz, quality )

!*****************************************************************************80
!
!! TETRAHEDRON_QUALITY1_3D: "quality" of a tetrahedron in 3D.
!
!  Discussion:
!
!    The quality of a tetrahedron is 3.0 times the ratio of the radius of 
!    the inscribed sphere divided by that of the circumscribed sphere.  
!
!    An equilateral tetrahredron achieves the maximum possible quality of 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) TET_XYZ(3,4), the coordinates of the vertices.
!
!    Output, real ( kind = 8 ) QUALITY, the quality of the tetrahedron.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  real ( kind = 8 ) pc(dim_num)
  real ( kind = 8 ) quality
  real ( kind = 8 ) r_in
  real ( kind = 8 ) r_out
  real ( kind = 8 ) tet_xyz(dim_num,4)

  call tetrahedron_circumsphere_3d ( tet_xyz, r_out, pc )

  call tetrahedron_insphere_3d ( tet_xyz, r_in, pc )

  quality = 3.0D+00 * r_in / r_out

  return
end
subroutine tetrahedron_quality2_3d ( tet_xyz, quality2 )

!*****************************************************************************80
!
!! TETRAHEDRON_QUALITY2_3D: "quality" of a tetrahedron in 3D.
!
!  Discussion:
!
!    The quality measure #2 of a tetrahedron is:
!
!      QUALITY2 = 2 * sqrt ( 6 ) * RIN / LMAX
!
!    where 
!
!      RIN = radius of the inscribed sphere;
!      LMAX = length of longest side of the tetrahedron.
!
!    An equilateral tetrahredron achieves the maximum possible quality of 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Qiang Du, Desheng Wang,
!    The Optimal Centroidal Voronoi Tesselations and the Gersho's
!    Conjecture in the Three-Dimensional Space,
!    Computers and Mathematics with Applications,
!    Volume 49, 2005, pages 1355-1373.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) TET_XYZ(3,4), the coordinates of the vertices.
!
!    Output, real ( kind = 8 ) QUALITY2, the quality of the tetrahedron.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  real ( kind = 8 ) edge_length(6)
  real ( kind = 8 ) l_max
  real ( kind = 8 ) pc(dim_num)
  real ( kind = 8 ) quality2
  real ( kind = 8 ) r_in
  real ( kind = 8 ) tet_xyz(dim_num,4)

  call tetrahedron_edge_length_3d ( tet_xyz, edge_length )

  l_max = maxval ( edge_length(1:6) )

  call tetrahedron_insphere_3d ( tet_xyz, r_in, pc )

  quality2 = 2.0D+00 * sqrt ( 6.0D+00 ) * r_in / l_max

  return
end
subroutine tetrahedron_quality3_3d ( tet_xyz, quality3 )

!*****************************************************************************80
!
!! TETRAHEDRON_QUALITY3_3D computes the mean ratio of a tetrahedron.
!
!  Discussion:
!
!    This routine computes QUALITY3, the eigenvalue or mean ratio of 
!    a tetrahedron.
!
!      QUALITY3 = 12 * ( 3 * volume )**(2/3) / (sum of squares of edge lengths).
!
!    This value may be used as a shape quality measure for the tetrahedron.
!
!    For an equilateral tetrahedron, the value of this quality measure
!    will be 1.  For any other tetrahedron, the value will be between
!    0 and 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 August 2005
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
!    Volume 13, Number 5, 1991, pages 325-331.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) TET_XYZ(3,4), the coordinates of the vertices.
!
!    Output, real ( kind = 8 ) QUALITY3, the mean ratio of the tetrahedron.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  real ( kind = 8 ) ab(dim_num)
  real ( kind = 8 ) ac(dim_num)
  real ( kind = 8 ) ad(dim_num)
  real ( kind = 8 ) bc(dim_num)
  real ( kind = 8 ) bd(dim_num)
  real ( kind = 8 ) cd(dim_num)
  real ( kind = 8 ) denom
  real ( kind = 8 ) lab
  real ( kind = 8 ) lac
  real ( kind = 8 ) lad
  real ( kind = 8 ) lbc
  real ( kind = 8 ) lbd
  real ( kind = 8 ) lcd
  real ( kind = 8 ) quality3
  real ( kind = 8 ) tet_xyz(dim_num,4)
  real ( kind = 8 ) volume
!
!  Compute the vectors representing the sides of the tetrahedron.
!
  ab(1:3) = tet_xyz( 1:dim_num,2) - tet_xyz( 1:dim_num,1)
  ac(1:3) = tet_xyz( 1:dim_num,3) - tet_xyz( 1:dim_num,1)
  ad(1:3) = tet_xyz( 1:dim_num,4) - tet_xyz( 1:dim_num,1)
  bc(1:3) = tet_xyz( 1:dim_num,3) - tet_xyz( 1:dim_num,2)
  bd(1:3) = tet_xyz( 1:dim_num,4) - tet_xyz( 1:dim_num,2)
  cd(1:3) = tet_xyz( 1:dim_num,4) - tet_xyz( 1:dim_num,3)
!
!  Compute the squares of the lengths of the sides.
!
  lab = sum ( ab(1:dim_num)**2 )
  lac = sum ( ac(1:dim_num)**2 )
  lad = sum ( ad(1:dim_num)**2 )
  lbc = sum ( bc(1:dim_num)**2 )
  lbd = sum ( bd(1:dim_num)**2 )
  lcd = sum ( cd(1:dim_num)**2 )
!
!  Compute the volume.
!
  volume = abs ( &
      ab(1) * ( ac(2) * ad(3) - ac(3) * ad(2) ) &
    + ab(2) * ( ac(3) * ad(1) - ac(1) * ad(3) ) &
    + ab(3) * ( ac(1) * ad(2) - ac(2) * ad(1) ) ) / 6.0D+00

  denom = lab + lac + lad + lbc + lbd + lcd

  if ( denom == 0.0D+00 ) then
    quality3 = 0.0D+00
  else
    quality3 = 12.0D+00 * ( 3.0D+00 * volume )**( 2.0D+00 / 3.0D+00 ) / denom
  end if

  return
end
subroutine tetrahedron_quality4_3d ( tet_xyz, quality4 )

!*****************************************************************************80
!
!! TETRAHEDRON_QUALITY4_3D computes the minimum solid angle of a tetrahedron.
!
!  Discussion:
!
!    This routine computes a quality measure for a tetrahedron, based
!    on the sine of half the minimum of the four solid angles.
!
!    The quality measure for an equilateral tetrahedron should be 1,
!    since the solid angles of such a tetrahedron are each equal to pi.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 August 2005
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
!    Volume 13, Number 5, 1991, pages 325-331.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) TET_XYZ(3,4), the coordinates of the vertices.
!
!    Output, real ( kind = 8 ) QUALITY4, the value of the quality measure.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  real ( kind = 8 ) a(dim_num)
  real ( kind = 8 ) ab(dim_num)
  real ( kind = 8 ) ac(dim_num)
  real ( kind = 8 ) ad(dim_num)
  real ( kind = 8 ) b(dim_num)
  real ( kind = 8 ) bc(dim_num)
  real ( kind = 8 ) bd(dim_num)
  real ( kind = 8 ) c(dim_num)
  real ( kind = 8 ) cd(dim_num)
  real ( kind = 8 ) d(dim_num)
  real ( kind = 8 ) denom
  real ( kind = 8 ) l1
  real ( kind = 8 ) l2
  real ( kind = 8 ) l3
  real ( kind = 8 ) lab
  real ( kind = 8 ) lac
  real ( kind = 8 ) lad
  real ( kind = 8 ) lbc
  real ( kind = 8 ) lbd
  real ( kind = 8 ) lcd
  real ( kind = 8 ) quality4
  real ( kind = 8 ) tet_xyz(dim_num,4)
  real ( kind = 8 ) volume
!
!  Compute the vectors that represent the sides.
!
  ab(1:dim_num) = tet_xyz( 1:dim_num,2) - tet_xyz( 1:dim_num,1)
  ac(1:dim_num) = tet_xyz( 1:dim_num,3) - tet_xyz( 1:dim_num,1)
  ad(1:dim_num) = tet_xyz( 1:dim_num,4) - tet_xyz( 1:dim_num,1)
  bc(1:dim_num) = tet_xyz( 1:dim_num,3) - tet_xyz( 1:dim_num,2)
  bd(1:dim_num) = tet_xyz( 1:dim_num,4) - tet_xyz( 1:dim_num,2)
  cd(1:dim_num) = tet_xyz( 1:dim_num,4) - tet_xyz( 1:dim_num,3)
!
!  Compute the lengths of the sides.
!
  lab = sqrt ( sum ( ab(1:dim_num)**2 ) )
  lac = sqrt ( sum ( ac(1:dim_num)**2 ) )
  lad = sqrt ( sum ( ad(1:dim_num)**2 ) )
  lbc = sqrt ( sum ( bc(1:dim_num)**2 ) )
  lbd = sqrt ( sum ( bd(1:dim_num)**2 ) )
  lcd = sqrt ( sum ( cd(1:dim_num)**2 ) )
!
!  Compute the volume
!
  volume = abs ( &
      ab(1) * ( ac(2) * ad(3) - ac(3) * ad(2) ) &
    + ab(2) * ( ac(3) * ad(1) - ac(1) * ad(3) ) &
    + ab(3) * ( ac(1) * ad(2) - ac(2) * ad(1) ) ) / 6.0D+00

  quality4 = 1.0D+00

  l1 = lab + lac
  l2 = lab + lad
  l3 = lac + lad

  denom = ( l1 + lbc ) * ( l1 - lbc ) &
        * ( l2 + lbd ) * ( l2 - lbd ) &
        * ( l3 + lcd ) * ( l3 - lcd )

  if ( denom <= 0.0D+00 ) then
    quality4 = 0.0D+00
  else
    quality4 = min ( quality4, 12.0D+00 * volume / sqrt ( denom ) )
  end if

  l1 = lab + lbc
  l2 = lab + lbd
  l3 = lbc + lbd

  denom = ( l1 + lac ) * ( l1 - lac ) &
        * ( l2 + lad ) * ( l2 - lad ) &
        * ( l3 + lcd ) * ( l3 - lcd )

  if ( denom <= 0.0D+00 ) then
    quality4 = 0.0D+00
  else
    quality4 = min ( quality4, 12.0D+00 * volume / sqrt ( denom ) )
  end if

  l1 = lac + lbc
  l2 = lac + lcd
  l3 = lbc + lcd

  denom = ( l1 + lab ) * ( l1 - lab ) &
        * ( l2 + lad ) * ( l2 - lad ) &
        * ( l3 + lbd ) * ( l3 - lbd )

  if ( denom <= 0.0D+00 ) then
    quality4 = 0.0D+00
  else
    quality4 = min ( quality4, 12.0D+00 * volume / sqrt ( denom ) )
  end if

  l1 = lad + lbd
  l2 = lad + lcd
  l3 = lbd + lcd

  denom = ( l1 + lab ) * ( l1 - lab ) &
        * ( l2 + lac ) * ( l2 - lac ) &
        * ( l3 + lbc ) * ( l3 - lbc )

  if ( denom <= 0.0D+00 ) then
    quality4 = 0.0D+00
  else
    quality4 = min ( quality4, 12.0D+00 * volume / sqrt ( denom ) )
  end if

  quality4 = quality4 * 1.5D+00 * sqrt ( 6.0D+00 )

  return
end
subroutine tetrahedron_reference_sample ( n, seed, p )

!*****************************************************************************80
!
!! TETRAHEDRON_REFERENCE_SAMPLE samples points in the reference tetrahedron.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 December 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of points to sample.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) P(3,N), random points in the tetrahedron.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3
  integer ( kind = 4 ) n

  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  real ( kind = 8 ) gamma
  integer ( kind = 4 ) j
  real ( kind = 8 ) p(dim_num,n)
  real ( kind = 8 ) r
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed

  do j = 1, n

    r = r8_uniform_01 ( seed )
!
!  Interpret R as a percentage of the tetrahedron's volume.
!
!  Imagine a plane, parallel to face 1, so that the volume between
!  vertex 1 and the plane is R percent of the full tetrahedron volume.
!
!  The plane will intersect sides 12, 13, and 14 at a fraction
!  ALPHA = R^1/3 of the distance from vertex 1 to vertices 2, 3, and 4.
!
    alpha = r**( 1.0D+00 / 3.0D+00 )
!
!  Determine the coordinates of the points on sides 12, 13 and 14 intersected
!  by the plane, which form a triangle TR.
!
!  Now choose, uniformly at random, a point in this triangle.
!
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
    beta = sqrt ( r )
!
!  Determine the coordinates of the points on sides 2 and 3 intersected
!  by line L.
!
!
!  Now choose, uniformly at random, a point on the line L.
!
    gamma = r8_uniform_01 ( seed )

    p(1:dim_num,j) = (/ &
      alpha * ( 1.0D+00 - beta ) *             gamma, &
      alpha *             beta   * ( 1.0D+00 - gamma ), &
      alpha *             beta   *             gamma /)

  end do

  return
end
subroutine tetrahedron_sample ( tet_xyz, n, seed, p )

!*****************************************************************************80
!
!! TETRAHEDRON_SAMPLE returns random points in a tetrahedron.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 December 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) TET_XYZ(3,4), the coordinates of the vertices.
!
!    Input, integer ( kind = 4 ) N, the  number of points to sample.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) P(3,N), random points in the tetrahedron.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3
  integer ( kind = 4 ) n

  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  real ( kind = 8 ) gamma
  integer ( kind = 4 ) j
  real ( kind = 8 ) p(dim_num,n)
  real ( kind = 8 ) p12(dim_num)
  real ( kind = 8 ) p13(dim_num)
  real ( kind = 8 ) r
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) tet_xyz(dim_num,4)
  real ( kind = 8 ) tr(dim_num,3)

  do j = 1, n

    r = r8_uniform_01 ( seed )
!
!  Interpret R as a percentage of the tetrahedron's volume.
!
!  Imagine a plane, parallel to face 1, so that the volume between
!  vertex 1 and the plane is R percent of the full tetrahedron volume.
!
!  The plane will intersect sides 12, 13, and 14 at a fraction
!  ALPHA = R^1/3 of the distance from vertex 1 to vertices 2, 3, and 4.
!
    alpha = r**( 1.0D+00 / 3.0D+00 )
!
!  Determine the coordinates of the points on sides 12, 13 and 14 intersected
!  by the plane, which form a triangle TR.
!
    tr(1:dim_num,1) = ( 1.0D+00 - alpha ) * tet_xyz( 1:dim_num,1) &
                                + alpha   * tet_xyz( 1:dim_num,2)
    tr(1:dim_num,2) = ( 1.0D+00 - alpha ) * tet_xyz( 1:dim_num,1) &
                                + alpha   * tet_xyz( 1:dim_num,3)
    tr(1:dim_num,3) = ( 1.0D+00 - alpha ) * tet_xyz( 1:dim_num,1) &
                                + alpha   * tet_xyz( 1:dim_num,4)
!
!  Now choose, uniformly at random, a point in this triangle.
!
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
    beta = sqrt ( r )
!
!  Determine the coordinates of the points on sides 2 and 3 intersected
!  by line L.
!
    p12(1:dim_num) = ( 1.0D+00 - beta ) * tr(1:dim_num,1) &
                               + beta   * tr(1:dim_num,2)
    p13(1:dim_num) = ( 1.0D+00 - beta ) * tr(1:dim_num,1) &
                               + beta   * tr(1:dim_num,3)
!
!  Now choose, uniformly at random, a point on the line L.
!
    gamma = r8_uniform_01 ( seed )

    p(1:dim_num,j) = ( 1.0D+00 - gamma ) * p12(1:dim_num) &
                   +             gamma   * p13(1:dim_num)


  end do

  return
end
subroutine tetrahedron_volume ( tet_xyz, volume )

!*****************************************************************************80
!
!! TETRAHEDRON_VOLUME computes the volume of a tetrahedron in 3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) TET_XYZ(3,4), the coordinates of the vertices.
!
!    Output, real ( kind = 8 ) VOLUME, the volume of the tetrahedron.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  real ( kind = 8 ) a(4,4)
  real ( kind = 8 ) r8mat_det_4d
  real ( kind = 8 ) tet_xyz(dim_num,4)
  real ( kind = 8 ) volume

  a(1:dim_num,1:4) = tet_xyz( 1:dim_num,1:4)
  a(4,1:4) = 1.0D+00

  volume = abs ( r8mat_det_4d ( a ) ) / 6.0D+00

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
