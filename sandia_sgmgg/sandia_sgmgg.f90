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
!    An I4 is an integer ( kind = 4 ) value.
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
function i4_wrap ( ival, ilo, ihi )

!*****************************************************************************80
!
!! I4_WRAP forces an I4 to lie between given limits by wrapping.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ) value.
!
!    There appears to be a bug in the GFORTRAN compiler which can lead to
!    erroneous results when the first argument of I4_WRAP is an expression.
!    In particular:
!
!    do i = 1, 3
!      if ( test ) then
!        i4 = i4_wrap ( i + 1, 1, 3 )
!      end if
!    end do
!
!    was, when I = 3, returning I4 = 3.  So I had to replace this with
!
!    do i = 1, 3
!      if ( test ) then
!        i4 = i + 1
!        i4 = i4_wrap ( i4, 1, 3 )
!      end if
!    end do
!
!  Example:
!
!    ILO = 4, IHI = 8
!
!    I  Value
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
!    07 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IVAL, a value.
!
!    Input, integer ( kind = 4 ) ILO, IHI, the desired bounds.
!
!    Output, integer ( kind = 4 ) I4_WRAP, a "wrapped" version of the value.
!
  implicit none

  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) value
  integer ( kind = 4 ) wide

  jlo = min ( ilo, ihi )
  jhi = max ( ilo, ihi )

  wide = jhi - jlo + 1

  if ( wide == 1 ) then
    value = jlo
  else
    value = jlo + i4_modp ( ival - jlo, wide )
  end if

  i4_wrap = value

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
  character ( len = 8 ) ctemp(incx)
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

      write ( *, '(i5,a,10a8)' ) j, ':', ( ctemp(i), i = 1, inc )

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

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i8,a,2x,i12)' ) i, ':', a(i)
  end do

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
!    For now, the input quantity SEED is an integer variable.
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

  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
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
    seed = seed + i4_huge
  end if
!
!  Although SEED can be represented exactly as a 32 bit integer,
!  it generally cannot be represented exactly as a 32 bit real number!
!
  r8_uniform_01 = real ( seed, kind = 8 ) * 4.656612875D-10

  return
end
subroutine r8vec_indexed_heap_d ( n, a, indx )

!*****************************************************************************80
!
!! R8VEC_INDEXED_HEAP_D creates a descending heap from an indexed R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    An indexed R8VEC is an R8VEC of data values, and an R8VEC of N indices, 
!    each referencing an entry of the data vector.
!
!    The function adjusts the index vector INDX so that, for 1 <= J <= N/2,
!    we have:
!      A(INDX(2*J))   <= A(INDX(J))
!    and
!      A(INDX(2*J+1)) <= A(INDX(J))
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Albert Nijenhuis, Herbert Wilf,
!    Combinatorial Algorithms for Computers and Calculators,
!    Academic Press, 1978,
!    ISBN: 0-12-519260-6,
!    LC: QA164.N54.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the index array.
!
!    Input, real ( kind = 8 ) A(*), the data vector.
!
!    Input/output, integer ( kind = 4 ) INDX(N), the index array.
!    Each entry of INDX must be a valid index for the array A.
!    On output, the indices have been reordered into a descending heap.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ifree
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) key
  integer ( kind = 4 ) m
!
!  Only nodes N/2 down to 1 can be "parent" nodes.
!
  do i = n / 2, 1, -1
!
!  Copy the value out of the parent node.
!  Position IFREE is now "open".
!
    key = indx(i)
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
        if ( a(indx(m)) < a(indx(m+1)) ) then
          m = m + 1
        end if

      end if
!
!  If the large descendant is larger than KEY, move it up,
!  and update IFREE, the location of the free position, and
!  consider the descendants of THIS position.
!
      if ( a(indx(m)) <= a(key) ) then
        exit
      end if

      indx(ifree) = indx(m)
      ifree = m

    end do
!
!  Once there is no more shifting to do, KEY moves into the free spot IFREE.
!
    indx(ifree) = key

  end do

  return
end
subroutine r8vec_indexed_heap_d_extract ( n, a, indx, indx_extract )

!*****************************************************************************80
!
!! R8VEC_INDEXED_HEAP_D_EXTRACT: extract from heap descending indexed R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    An indexed R8VEC is an R8VEC of data values, and an R8VEC of N indices, 
!    each referencing an entry of the data vector.
!
!    The routine finds the maximum value in the heap, returns that value to the 
!    user, deletes that value from the heap, and restores the heap to its 
!    proper form.
!
!    Note that the argument N must be a variable, which will be decremented
!    before return, and that INDX will hold one less value on output than it 
!    held on input.
!
!    This is one of three functions needed to model a priority queue.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Thomas Cormen, Charles Leiserson, Ronald Rivest,
!    Introduction to Algorithms,
!    MIT Press, 2001,
!    ISBN: 0262032937,
!    LC: QA76.C662.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N, the number of items in the 
!    index vector.
!
!    Input, real ( kind = 8 ) A(*), the data vector.
!
!    Input/output, integer ( kind = 4 ) INDX(N), the index vector.
!
!    Output, integer ( kind = 4 ) INDX_EXTRACT, the index in A of the item of 
!    maximum value, which has now been removed from the heap.
!
  implicit none

  real ( kind = 8 ) a(*)
  integer ( kind = 4 ) indx(*)
  integer ( kind = 4 ) indx_extract
  integer ( kind = 4 ) n

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_INDEXED_HEAP_D_EXTRACT - Fatal error!'
    write ( *, '(a)' ) '  The heap is empty.'
    stop
  end if
!
!  Get the index of the maximum value.
!
  indx_extract = indx(1)

  if ( n == 1 ) then
    n = 0
    return
  end if
!
!  Shift the last index down.
!
  indx(1) = indx(n)
!
!  Restore the heap structure.
!
  n = n - 1
  call r8vec_indexed_heap_d ( n, a, indx )

  return
end
subroutine r8vec_indexed_heap_d_insert ( n, a, indx, indx_insert )

!*****************************************************************************80
!
!! R8VEC_INDEXED_HEAP_D_INSERT: insert value into heap descending indexed R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    An indexed R8VEC is an R8VEC of data values, and an R8VEC of N indices, 
!    each referencing an entry of the data vector.
!
!    Note that the argument N must be a variable, and will be incremented before
!    return, and that INDX must be able to hold one more entry on output than 
!    it held on input.
!
!    This is one of three functions needed to model a priority queue.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Thomas Cormen, Charles Leiserson, Ronald Rivest,
!    Introduction to Algorithms,
!    MIT Press, 2001,
!    ISBN: 0262032937,
!    LC: QA76.C662.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N, the number of items in the 
!    index vector.
!
!    Input, real ( kind = 8 ) A(*), the data vector.
!
!    Input/output, integer ( kind = 4 ) INDX(N), the index vector.
!
!    Input, integer ( kind = 4 ) INDX_INSERT, the index in A of the value 
!    to be inserted into the heap.
!
  implicit none

  real ( kind = 8 ) a(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(*)
  integer ( kind = 4 ) indx_insert
  integer ( kind = 4 ) n
  integer ( kind = 4 ) parent

  n = n + 1
  i = n

  do while ( 1 < i )

    parent = i / 2

    if ( a(indx_insert) <= a(indx(parent)) ) then
      exit
    end if

    indx(i) = indx(parent)
    i = parent

  end do

  indx(i) = indx_insert

  return
end
subroutine r8vec_indexed_heap_d_max ( n, a, indx, indx_max )

!*****************************************************************************80
!
!! R8VEC_INDEXED_HEAP_D_MAX: maximum value in heap descending indexed R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    An indexed R8VEC is an R8VEC of data values, and an R8VEC of N indices, 
!    each referencing an entry of the data vector.
!
!    This is one of three functions needed to model a priority queue.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Thomas Cormen, Charles Leiserson, Ronald Rivest,
!    Introduction to Algorithms,
!    MIT Press, 2001,
!    ISBN: 0262032937,
!    LC: QA76.C662.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items in the index vector.
!
!    Input, real ( kind = 8 ) A(*), the data vector.
!
!    Input, integer ( kind = 4 ) INDX(N), the index vector.
!
!    Output, integer ( kind = 4 ) INDX_MAX, the index in A of the maximum value 
!    in the heap.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(*)
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) indx_max

  indx_max = indx(1)

  return
end
subroutine sandia_sgmgg_coef_inc2 ( m, n1, s1, c1, s2, c3 )

!*****************************************************************************80
!
!! SANDIA_SGMGG_COEF_INC2 computes tentative coefficient changes.
!
!  Discussion:
!
!    An active set S1 of N1 sparse grid indices is given, each of
!    size M.
!
!    The coefficient C1 of each sparse grid index is also given.
!
!    A candidate sparse grid index S2 is provided.
!
!    This function determines the N+1 coefficients that would be
!    appropriate if the candidate S2 was added to the active set
!    as the (N+1)-st item.
!
!    During the calculation, we may try to update coefficients of inactive 
!    index sets.  By the end of the calculation, all these inactive index
!    sets should have accumulated total coefficients of 0 again.  As a check,
!    we temporarily set aside space for these objects, and check, at the end,
!    that the coefficients are zero.
!
!  Example:
!
!    Input:
!
!      +1 * {0,2}
!      -1 * {0,1}  +1 * {1,1}
!                  -1 * {1,0}  +1 * {2,0}
!
!    Add {3,0}
!
!    Output:
!
!      +1 * {0,2}
!      -1 * {0,1}  +1 * {1,1}
!                  -1 * {1,0}   0 * {2,0}  +1 * {3,0}
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 August 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the dimension of the vector.
!
!    Input, integer ( kind = 4 ) N1, the number of points in the active set.
!
!    Input, integer ( kind = 4 ) S1(M,N1), the indices for the active set.
!
!    Input, integer ( kind = 4 ) C1(N1), the coefficients for the active set.
!
!    Input, integer ( kind = 4 ) S2(M), the indices for the candidate.
!
!    Output, integer ( kind = 4 ) C3(N1+1), the coefficients for the active set
!    plus the candidate.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n1

  integer ( kind = 4 ) c1(n1)
  integer ( kind = 4 ) c3(n1+1)
  integer ( kind = 4 ) c4(n1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n4
  integer ( kind = 4 ) s(m)
  integer ( kind = 4 ) s1(m,n1)
  integer ( kind = 4 ) s2(m)
  integer ( kind = 4 ) s4(m,n1)
!
!  Initialize the inactive data.
!
  n4 = 0
  c4(1:n1) = 0
  s4(1:m,1:n1) = 0
!
!  Copy C1.
!
  c3(1:n1) = c1(1:n1)
  c3(n1+1) = 1
!
!  Consider the effect of the new item S2 on each of the current
!  items in the active set S1.
!
  do j = 1, n1
!
!  Determine S, the element-wise minimum of the J-th item in S1 versus S2.
!
    k = j

    do i = 1, m
 
      if ( s2(i) < s1(i,j) ) then
        s(i) = s2(i)
        k = -1
      else
        s(i) = s1(i,j)
      end if

    end do
!
!  If S = S1(*,J), K is J.
!
    if ( k /= -1 ) then

      c3(k) = c3(k) - c1(j)
!
!  If S is equal to an element of the active set, we set K to that index.
! 
    else

      do j2 = 1, n1

        k = j2

        do i2 = 1, m

          if ( s1(i2,j2) /= s(i2) ) then
            k = -1
            exit
          end if

        end do

        if ( k /= -1 ) then
          c3(k) = c3(k) - c1(j)
          exit
        end if

      end do

    end if
!
!  If S is equal to an element of the inactive set, set K to that index.
!
    if ( k == -1 ) then

      do j2 = 1, n4

        k = j2

        do i2 = 1, m

          if ( s4(i2,j2) /= s(i2) ) then
            k = - 1
            exit
          end if

        end do

        if ( k /= - 1 ) then
          c4(k) = c4(k) - c1(j)
          exit
        end if

      end do

    end if
!
!  S is not equal to S1(*,J), or any element of S1, or any element of S4.
!  Add S to the set of elements S4.
!
    if ( k == -1 ) then
      n4 = n4 + 1
      k = n4
      c4(k) = 0
      s4(1:m,k) = s(1:m)
      c4(k) = c4(k) - c1(j)
    end if

  end do
!
!  At the end, the C4(1:N4) should all be zero.
!
  if ( any ( c4(1:n4) /= 0 ) ) then
    write ( *, '(a)') ' '
    write ( *, '(a)' ) 'SANDIA_SGMGG_COEF_INC2 - Fatal error!'
    write ( *, '(a)' ) '  Some inactive indices were assigned a nonzero coefficient.'
    call i4mat_transpose_print ( m, n4, s4, '  S4:' )
    call i4vec_print ( n4, c4, '  C4:' )
    stop
  end if

  return
end
subroutine sandia_sgmgg_coef_naive ( dim_num, point_num, sparse_index, coef )

!*****************************************************************************80
!
!! SANDIA_SGMGG_COEF_NAIVE returns the combinatorial coefficients.
!
!  Discussion:
!
!    The coefficient of point I is calculated as follows:
!
!    *) point J is a "neighbor" of point I if every entry of the sparse
!       index for point J is either equal to, or 1 greater than, the
!       corresponding entry of the sparse index of point I.
!
!    *) If point J is a neighbor of point I, then it contributes
!       (-1)^D to the coefficient, where D is the sum of the differences
!       between the sparse indices of point I and point J.
!
!    This is a completely naive implementation of the calculation,
!    intended simply as a demonstration for small examples.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2010
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
!    Fabio Nobile, Raul Tempone, Clayton Webster,
!    An Anisotropic Sparse Grid Stochastic Collocation Method for Partial
!    Differential Equations with Random Input Data,
!    SIAM Journal on Numerical Analysis,
!    Volume 46, Number 5, 2008, pages 2411-2442.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the vector.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, integer ( kind = 4 ) SPARSE_INDEX(DIM_NUM,POINT_NUM),
!    the indices that define the points.
!
!    Output, integer ( kind = 4 ) COEF(POINT_NUM), the coefficients.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) coef(point_num)
  integer ( kind = 4 ) dif
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  logical neighbor
  integer ( kind = 4 ) sparse_index(dim_num,point_num)
  integer ( kind = 4 ) term

  coef(1:point_num) = 0

  do j1 = 1, point_num

    do j2 = 1, point_num

      neighbor = .true.
      term = + 1

      do i = 1, dim_num

        dif = sparse_index(i,j2) - sparse_index(i,j1)

        if ( dif == 0 ) then

        else if ( dif == 1 ) then
          term = - term
        else
          neighbor = .false.
          exit
        end if

      end do

      if ( neighbor ) then
        coef(j1) = coef(j1) + term
      end if

    end do

  end do

  return
end
subroutine sandia_sgmgg_neighbor_naive ( dim_num, point_num, sparse_index, &
  neighbor )

!*****************************************************************************80
!
!! SANDIA_SGMGG_NEIGHBOR_NAIVE returns the immediate forward neighbor vector.
!
!  Discussion:
!
!    A sparse grid index vector is a list of DIM_NUM nonnegative indices.
!
!    An immediate forward L-neighbor of a sparse grid index vector I is a
!    sparse grid index vector J with the property that all entries of J 
!    are equal to those of I except for the L-the entry, which is greater by 1.
!
!    A forward neighbor of a sparse grid index vector I is a sparse
!    grid index vector K with the property that every entry of K is
!    equal to or greater by 1 than the corresponding entry of I.
!
!    This function is given a collection of sparse grid index vectors,
!    and returns information defining, for every such vector, the entire
!    set of its immediate forward neighbors.  This is done with a 
!    "NEIGHBOR" array of dimension DIM_NUM.  For sparse grid vector I,
!    entry L of NEIGHBOR is 1 if I has an immediate forward L-neighbor,
!    and 0 otherwise.
!
!    This implementation of the procedure is inefficient, and is provided
!    solely for demonstration on small problems.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 August 2010
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
!    Fabio Nobile, Raul Tempone, Clayton Webster,
!    An Anisotropic Sparse Grid Stochastic Collocation Method for Partial
!    Differential Equations with Random Input Data,
!    SIAM Journal on Numerical Analysis,
!    Volume 46, Number 5, 2008, pages 2411-2442.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the dimension of the vector.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, integer ( kind = 4 ) SPARSE_INDEX(DIM_NUM,POINT_NUM),
!    the indices that define the points.
!
!    Output, integer ( kind = 4 ) NEIGHBOR(DIM_NUM,POINT_NUM), the 
!    immediate forward neighbor array.
!
  implicit none

  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) point_num

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) l
  integer ( kind = 4 ) neighbor(dim_num,point_num)
  integer ( kind = 4 ) sparse_index(dim_num,point_num)

  neighbor(1:dim_num,1:point_num) = 0

  do j1 = 1, point_num

    do j2 = 1, point_num

      l = -1

      do i = 1, dim_num
!
!  If the entries are not equal...
!
        if ( sparse_index(i,j2) /= sparse_index(i,j1) ) then
!
!  ...and we haven't already found a difference...
!
          if ( l /= -1 ) then
            l = - 1
            exit
          end if
!
!  ...and this difference is +1...
!
          if ( sparse_index(i,j2) /= sparse_index(i,j1) + 1 ) then
            exit
          endif
!
!  ...then remember this index.
!
          l = i

        end if

      end do
!
!  If a single unit difference was found, record the direction.
!
      if ( l /= -1 ) then
        neighbor(l,j1) = 1
      end if

    end do

  end do

  return
end
subroutine sgmgg_print ( gg_ma, gg_mi, gg_mo, gg_na, gg_nd, gg_ni, gg_no, &
  gg_a, gg_b, gg_f, gg_g, gg_i, gg_o, gg_s )

!*****************************************************************************80
!
!! SGMGG_PRINT prints out an SGMGG data structure.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 August 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) GG_MA, the maximum dimension for GG_A.
!
!    Input, integer ( kind = 4 ) GG_MI, the maximum dimension for GG_I.
!
!    Input, integer ( kind = 4 ) GG_MO, the maximum dimension for GG_O.
!
!    Input, integer ( kind = 4 ) GG_NA, the current dimension for GG_A.
!
!    Input, integer ( kind = 4 ) GG_ND, the spatial dimension.
!
!    Input, integer ( kind = 4 ) GG_NI, the current dimension for GG_I.
!
!    Input, integer ( kind = 4 ) GG_NO, the current dimension for GG_O.
!
!    Input, integer ( kind = 4 ) GG_A(GG_MA), the active indices.
!
!    Input, integer ( kind = 4 ) GG_B(GG_ND,GG_MI), the forward neighbors.
!
!    Input, integer ( kind = 4 ) GG_F(GG_ND,GG_MI), the backward neighbors.
!
!    Input, real ( kind = 8 ) GG_G(GG_MA), the error estimators.
!
!    Input, integer ( kind = 4 ) GG_I(GG_ND,GG_MI), the index set.
!
!    Input, integer ( kind = 4 ) GG_O(GG_MO), the old indices.
!
!    Input, integer ( kind = 4 ) GG_S(GG_MI), 0 if index I is old, 1
!    if it is active.
!
  implicit none

  integer ( kind = 4 ) gg_ma
  integer ( kind = 4 ) gg_mi
  integer ( kind = 4 ) gg_mo
  integer ( kind = 4 ) gg_nd

  integer ( kind = 4 ) gg_a(gg_ma)
  integer ( kind = 4 ) gg_b(gg_nd,gg_mi)
  integer ( kind = 4 ) gg_f(gg_nd,gg_mi)
  real ( kind = 8 ) gg_g(gg_mi)
  integer ( kind = 4 ) gg_i(gg_nd,gg_mi)
  integer ( kind = 4 ) gg_na
  integer ( kind = 4 ) gg_ni
  integer ( kind = 4 ) gg_no
  integer ( kind = 4 ) gg_o(gg_mo)
  integer ( kind = 4 ) gg_s(gg_mi)
  integer ( kind = 4 ) i
  character s

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GG DATA STRUCTURE:'
  write ( *, '(a,i4)' ) '  ND = ', gg_nd
  write ( *, '(a,i4)' ) '  NI = ', gg_ni
  write ( *, '(a,i4)' ) '  NO = ', gg_no
  write ( *, '(a,i4)' ) '  NA = ', gg_na
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Indices:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I       A/O    G(I)        Index values'
  write ( *, '(a)' ) ' '
  do i = 1, gg_ni
    if ( gg_s(i) == 0 ) then
      s = 'o'
    else if ( gg_s(i) == 1 ) then
      s = 'a'
    else
      s = '?'
    end if
    write ( *, '(2x,i4,6x,2x,a1,2x,g14.6,5(2x,i4))' ) &
      i, s, gg_g(i), gg_i(1:gg_nd,i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Backward neighbors:'
  write ( *, '(a)' ) ' '
  do i = 1, gg_ni
    write ( *, '(2x,i4,6x,5(2x,i4))' ) i, gg_b(1:gg_nd,i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Forward neighbors:'
  write ( *, '(a)' ) ' '
  do i = 1, gg_ni
    write ( *, '(2x,i4,6x,5(2x,i4))' ) i, gg_f(1:gg_nd,i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Active Heap:'
  write ( *, '(a)' ) '     I     A      G'
  write ( *, '(a)' ) ' '

  do i = 1, gg_na
    write ( *, '(2x,i4,2x,i4,2x,g14.6)' ) i, gg_a(i), gg_g(gg_a(i))
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
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9  ), parameter, dimension(12) :: month = (/ &
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
