subroutine b4set_colex_rank ( n, t, rank )

!*****************************************************************************80
!
!! B4SET_COLEX_RANK computes the colexicographic rank of a B4SET.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 August 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the master set.
!
!    Input, integer ( kind = 4 ) T, the set.
!
!    Output, integer ( kind = 4 ) RANK, the rank of the set.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) t

  rank = 0

  do i = 0, n - 1

    if ( btest ( t, i ) ) then
      rank = rank + 2 ** i
    end if

  end do

  return
end
subroutine b4set_colex_successor ( n, t, rank )

!*****************************************************************************80
!
!! B4SET_COLEX_SUCCESSOR computes the colexicographic successor of a B4SET.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 August 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the master set.
!
!    Input/output, integer ( kind = 4 ) T, describes a set.  
!    On input, T describes a set.
!    On output, T describes the next set in the ordering.
!    If the input T was the last in the ordering, then the output T
!    will be the first.
!
!    Input/output, integer ( kind = 4 ) RANK, the rank.
!    If RANK = -1 on input, then the routine understands that this is
!    the first call, and that the user wishes the routine to supply
!    the first element in the ordering, which has RANK = 0.
!    In general, the input value of RANK is increased by 1 for output,
!    unless the very last element of the ordering was input, in which
!    case the output value of RANK is 0.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) t
!
!  Return the first element.
!
  if ( rank == -1 ) then
    t = 0
    rank = 0
    return
  end if

  do i = 0, n - 1

    if ( .not. btest ( t, i ) ) then
      t = ibset ( t, i )
      rank = rank + 1
      return
    else
      t = ibclr ( t, i )
    end if

  end do

  rank = 0

  return
end
subroutine b4set_colex_unrank ( rank, n, t )

!*****************************************************************************80
!
!! B4SET_COLEX_UNRANK computes the B4SET of given colexicographic rank.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 August 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) RANK, the rank of the set.
!
!    Input, integer ( kind = 4 ) N, the order of the master set.
!
!    Output, integer ( kind = 4 ) T, the set of the given rank.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) rank_copy
  integer ( kind = 4 ) sub_num
  integer ( kind = 4 ) t
!
!  Check.
!
  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'B4SET_COLEX_UNRANK - Fatal error!'
    write ( *, '(a)' ) '  Input N is illegal.'
    stop
  end if

  call b4set_enum ( n, sub_num )

  if ( rank < 0 .or. sub_num < rank ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'B4SET_COLEX_UNRANK - Fatal error!'
    write ( *, '(a)' ) '  The input rank is illegal.'
    stop
  end if

  rank_copy = rank

  t = 0

  do i = 0, n - 1
    if ( mod ( rank_copy, 2 ) == 1 ) then
      t = ibset ( t, i )
    else
      t = ibclr ( t, i )
    end if

    rank_copy = rank_copy / 2

  end do

  return
end
subroutine b4set_complement ( n, a, b )

!*****************************************************************************80
!
!! B4SET_COMPLEMENT computes the complement of a B4SET.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 August 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the master set.
!
!    Input, integer ( kind = 4 ) A, the set.
!
!    Output, integer ( kind = 4 ) B, the complement of A.
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n

  do i = 1, n
    if ( .not. btest ( a, i - 1 ) ) then
      b = ibset ( b, i - 1 )
    else
      b = ibclr ( b, i - 1 )
    end if
  end do

  return
end
subroutine b4set_complement_relative ( n, a, b, c )

!*****************************************************************************80
!
!! B4SET_COMPLEMENT_RELATIVE computes the relative complement of a B4SET.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 August 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the master set.
!
!    Input, integer ( kind = 4 ) A, the set.
!
!    Input, integer ( kind = 4 ) B, the set with respect to which 
!    the complement is taken.
!
!    Output, integer ( kind = 4 ) C, the complement of A with respect to B.
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n

  do i = 1, n
    if (  btest ( a, i - 1 ) .and. .not. btest ( b, i - 1 ) ) then
      c = ibset ( c, i - 1 )
    else
      c = ibclr ( c, i - 1 )
    end if
  end do

  return
end
subroutine b4set_delete ( n, a, t )

!*****************************************************************************80
!
!! B4SET_DELETE deletes an element from a B4SET.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 August 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the master set.
!
!    Input, integer ( kind = 4 ) A, an item.
!
!    Input/output, integer ( kind = 4 ) T, a set.
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) n
  integer ( kind = 4 ) t

  if ( a < 1 .or. n < a ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'B4SET_DELETE - Fatal error!'
    write ( *, '(a)' ) '  1 <= A <= N fails.'
    stop
  end if

  t = ibclr ( t, a - 1 )

  return
end
subroutine b4set_distance ( n, t1, t2, dist )

!*****************************************************************************80
!
!! B4SET_DISTANCE computes the Hamming distance between two B4SET's.
!
!  Discussion:
!
!    The sets T1 and T2 are assumed to be subsets of a set of N elements.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 August 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the master set.
!
!    Input, integer ( kind = 4 ) T1, T2, two sets.
!
!    Output, integer ( kind = 4 ) DIST, the Hamming distance between T1 and T2,
!    defined as the number of elements of the master set which are
!    in either T1 or T2 but not both.
!
  implicit none

  integer ( kind = 4 ) dist
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) t1
  integer ( kind = 4 ) t2

  dist = 0

  do i = 1, n

    if ( btest ( t1, i - 1 ) .neqv. btest ( t2, i - 1 ) ) then
      dist = dist + 1
    end if

  end do

  return
end
subroutine b4set_enum ( n, set_num )

!*****************************************************************************80
!
!! B4SET_ENUM enumerates the B4SET's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 August 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the master set.
!
!    Output, integer ( kind = 4 ) SET_NUM, the number of distinct sets.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) set_num

  set_num = 2**n

  return
end
function b4set_index ( n, a, t )

!*****************************************************************************80
!
!! B4SET_INDEX returns the index of an element of a B4SET.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 August 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the master set.
!
!    Input, integer ( kind = 4 ) A, the item.
!
!    Input, integer ( kind = 4 ) T, a set.
!
!    Output, integer B4SET_INDEX, the index of the item in the set,
!    or -1 if the item is not in the set.
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b4set_index
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) t
  integer ( kind = 4 ) value

  if ( a < 1 .or. n < a ) then
    value = -1
  else
    value = 0
    do i = 1, a
      if ( btest ( t, i - 1 ) ) then
        value = value + 1
      end if
    end do
  end if

  b4set_index = value

  return
end
subroutine b4set_insert ( n, a, t )

!*****************************************************************************80
!
!! B4SET_INSERT inserts an item into a B4SET.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 August 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the master set.
!
!    Input, integer ( kind = 4 ) A, the item.
!    1 <= A <= N.
!
!    Input/output, integer ( kind = 4 ) T, a set.
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) n
  integer ( kind = 4 ) t

  if ( a < 1 .or. n < a ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'B4SET_INSERT - Fatal error!'
    write ( *, '(a)' ) '  1 <= A <= N fails.'
    stop
  end if

  t = ibset ( t, a - 1 )

  return
end
subroutine b4set_intersect ( n, a, b, c )

!*****************************************************************************80
!
!! B4SET_INTERSECT computes the intersection of two B4SET's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 September 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the master set.
!
!    Input, integer ( kind = 4 ) A, B, two sets.
!
!    Output, integer ( kind = 4 ) C, the intersection of A and B.
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n

  c = 0

  do i = 1, n
    if ( btest ( a, i - 1 ) .and. btest ( b, i - 1 ) ) then
      c = ibset ( c, i - 1 )
    else
      c = ibclr ( c, i - 1 )
    end if
  end do

  return
end
function b4set_is_empty ( n, t )

!*****************************************************************************80
!
!! B4SET_IS_EMPTY determines if a B4SET is empty.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 August 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the master set.
!
!    Input, integer ( kind = 4 ) T, a set.
!
!    Output, logical B4SET_IS_EMPTY is TRUE if T is empty.
!
  implicit none

  logical b4set_is_empty
  integer ( kind = 4 ) n
  integer ( kind = 4 ) t

  b4set_is_empty = ( t == 0 )

  return
end
function b4set_is_equal ( n, t1, t2 )

!*****************************************************************************80
!
!! B4SET_IS_EQUAL determines if two B4SET's are equal.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 August 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the master set.
!
!    Input, logical T1, T2, two sets.
!
!    Output, logical B4SET_IS_EQUAL, is TRUE if T1 equals T2.
!
  implicit none

  logical b4set_is_equal
  integer ( kind = 4 ) n
  integer ( kind = 4 ) t1
  integer ( kind = 4 ) t2

  b4set_is_equal = ( t1 == t2 )

  return
end
function b4set_is_member ( n, a, t )

!*****************************************************************************80
!
!! B4SET_IS_MEMBER determines if an item is a member of a B4SET.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 August 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the master set.
!
!    Input, integer ( kind = 4 ) A, an item.
!
!    Input, integer ( kind = 4 ) T, a set.
!
!    Output, logical B4SET_IS_MEMBER, is TRUE if A is an element of T.
!
  implicit none

  integer ( kind = 4 ) a
  logical b4set_is_member
  integer ( kind = 4 ) n
  integer ( kind = 4 ) t

  if ( 1 <= a .and. a <= n ) then
    b4set_is_member = btest ( t, a - 1 )
  else
    b4set_is_member = .false.
  end if

  return
end
function b4set_is_subset ( n, t1, t2 )

!*****************************************************************************80
!
!! B4SET_IS_SUBSET determines if one B4SET is a subset of another.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 August 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the master set.
!
!    Input, integer ( kind = 4 ) T1, T2, two sets.
!
!    Output, logical B4SET_IS_SUBSET, is TRUE if T1 is a subset of T2.
!
  implicit none

  logical b4set_is_subset
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) t1
  integer ( kind = 4 ) t2

  b4set_is_subset = .true.

  do i = 1, n
    if ( btest ( t1, i - 1 ) .and. .not. btest ( t2, i - 1 ) ) then
      b4set_is_subset = .false.
      return
    end if
  end do

  return
end
subroutine b4set_lex_rank ( n, t, rank )

!*****************************************************************************80
!
!! B4SET_LEX_RANK computes the lexicographic rank of a B4SET.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 August 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the master set.
!
!    Input, integer ( kind = 4 ) T, the set.
!
!    Output, integer ( kind = 4 ) RANK, the rank of the set.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) t

  rank = 0

  do i = 0, n - 1

    if ( btest ( t, i ) ) then
      rank = rank + 2 ** ( n - i - 1 )
    end if

  end do

  return
end
subroutine b4set_lex_successor ( n, t, rank )

!*****************************************************************************80
!
!! B4SET_LEX_SUCCESSOR computes the lexicographic successor of a B4SET.
!
!  Discussion:
!
!    In the original code, there is a last element with no successor.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 August 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the master set.
!
!    Input/output, integer ( kind = 4 ) T, describes a set.
!    On input, T describes a set.
!    On output, T describes the next set in the ordering.
!    If the input T was the last in the ordering, then the output T
!    will be the first.
!
!    Input/output, integer ( kind = 4 ) RANK, the rank.
!    If RANK = -1 on input, then the routine understands that this is
!    the first call, and that the user wishes the routine to supply
!    the first element in the ordering, which has RANK = 0.
!    In general, the input value of RANK is increased by 1 for output,
!    unless the very last element of the ordering was input, in which
!    case the output value of RANK is 0.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) t
!
!  Return the first element.
!
  if ( rank == -1 ) then
    t = 0
    rank = 0
    return
  end if

  do i = n - 1, 0, -1

    if ( .not. btest ( t, i ) ) then
      t = ibset ( t, i )
      rank = rank + 1
      return
    else
      t = ibclr ( t, i )
    end if

  end do

  rank = 0

  return
end
subroutine b4set_lex_unrank ( rank, n, t )

!*****************************************************************************80
!
!! B4SET_LEX_UNRANK computes the B4SET of given lexicographic rank.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 August 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) RANK, the rank of the set.
!
!    Input, integer ( kind = 4 ) N, the order of the master set.
!
!    Output, integer ( kind = 4 ) T, the set of the given rank.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) nsub
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) rank_copy
  integer ( kind = 4 ) set_num
  integer ( kind = 4 ) t
!
!  Check.
!
  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'B4SET_LEX_UNRANK - Fatal error!'
    write ( *, '(a)' ) '  Input N is illegal.'
    stop
  end if

  call b4set_enum ( n, set_num )

  if ( rank < 0 .or. set_num < rank ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'B4SET_LEX_UNRANK - Fatal error!'
    write ( *, '(a)' ) '  The input rank is illegal.'
    stop
  end if

  rank_copy = rank

  t = 0

  do i = n - 1, 0, -1

    if ( mod ( rank_copy, 2 ) == 1 ) then
      t = ibset ( t, i )
    else
      t = ibclr ( t, i )
    end if

    rank_copy = rank_copy / 2

  end do

  return
end
subroutine b4set_to_lset ( n, t, a )

!*****************************************************************************80
!
!! B4SET_TO_LSET converts a B4SET to an LSET.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 August 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the master set.
!
!    Input, integer ( kind = 4 ) T, the set.
!
!    Input, logical A(N), the LSET version of the set.
!
  implicit none

  integer ( kind = 4 ) n

  logical a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) t

  do i = 1, n
    a(i) = btest ( t, i - 1 )
  end do

  return
end
subroutine b4set_transpose_print ( n, t, title )

!*****************************************************************************80
!
!! B4SET_TRANSPOSE_PRINT prints a B4SET "transposed".
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 August 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the master set.
!
!    Input, integer ( kind = 4 ) T, the set.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) t
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  if ( t == 0 ) then
    write ( *, '(a)' ) '  (Empty set)'
  else
    s = 0
    do i = 1, n
      if ( btest ( t, i - 1 ) ) then
        write ( *, '(2x,i2)', advance = 'no' ) i
        s = s + 4
      end if
      if ( 76 < s .or. ( 0 < s .and. i == n ) ) then
        write ( *, '(1x)', advance = 'yes' )
        s = 0
      end if
    end do
  end if

  return
end
subroutine b4set_random ( n, seed, a )

!*****************************************************************************80
!
!! B4SET_RANDOM sets a rondom B4SET.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 August 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the master set.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, integer ( kind = 4 ) A, the random B4SET.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a
  integer ( kind = 4 ) a_logical(n)
  integer ( kind = 4 ) seed

  call lset_random ( n, seed, a_logical )
  call lset_to_b4set ( n, a_logical, a )

  return
end
subroutine b4set_union ( n, a, b, c )

!*****************************************************************************80
!
!! B4SET_UNION computes the union of two B4SET's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 August 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the master set.
!
!    Input, integer ( kind = 4 ) A, B, two sets.
!
!    Output, integer ( kind = 4 ) C, the union of A and B.
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n

  do i = 1, n
    if ( btest ( a, i - 1 ) .or. btest ( b, i - 1 ) ) then
      c = ibset ( c, i - 1 )
    else
      c = ibclr ( c, i - 1 )
    end if
  end do

  return
end
subroutine b4set_weight ( n, t, weight )

!*****************************************************************************80
!
!! B4SET_WEIGHT computes the Hamming weight of a B4SET.
!
!  Discussion:
!
!    The Hamming weight is simply the number of elements in the set.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 August 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the master set..
!
!    Input, integer ( kind = 4 ) T, the set.
!
!    Output, integer ( kind = 4 ) WEIGHT, the Hamming weight of the set T.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) t
  integer ( kind = 4 ) weight

  weight = 0

  do i = 1, n
    if ( btest ( t, i - 1 ) ) then
      weight = weight + 1
    end if
  end do

  return
end
subroutine b4set_xor ( n, a, b, c )

!*****************************************************************************80
!
!! B4SET_XOR computes the symmetric difference of two B4SET's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 August 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the master set.
!
!    Input, integer ( kind = 4 ) A, B, two sets.
!
!    Output, integer ( kind = 4 ) C, the symmetric difference of A and B.
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n

  do i = 1, n
    if ( btest ( a, i - 1 ) .neqv. btest ( b, i - 1 ) ) then
      c = ibset ( c, i - 1 )
    else
      c = ibclr ( c, i - 1 )
    end if
  end do

  return
end
subroutine digit_to_ch ( digit, ch )

!*****************************************************************************80
!
!! DIGIT_TO_CH returns the character representation of a decimal digit.
!
!  Discussion:
!
!    Instead of CHAR, we now use the ACHAR function, which
!    guarantees the ASCII collating sequence.
!
!  Example:
!
!    DIGIT   CH
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
!    Input, integer ( kind = 4 ) DIGIT, the digit value between 0 and 9.
!
!    Output, character CH, the corresponding character.
!
  implicit none

  character ch
  integer ( kind = 4 ) digit

  if ( 0 <= digit .and. digit <= 9 ) then

    ch = achar ( digit + 48 )

  else

    ch = '*'

  end if

  return
end
subroutine i4_to_s_right ( intval, s )

!*****************************************************************************80
!
!! I4_TO_S_RIGHT converts an I4 to a right justified string.
!
!  Discussion:
!
!    An I4 is an integer ( kind = 4 ).
!
!  Example:
!
!    Assume that S is 6 characters long:
!
!    INTVAL       S
!
!         1       1
!        -1      -1
!         0       0
!      1952    1952
!    123456  123456
!   1234567  ******  <-- Not enough room!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) INTVAL, an integer to be converted.
!
!    Output, character ( len = * ) S, the representation of the integer.
!    The integer will be right-justified.  If there is not enough space,
!    the string will be filled with stars.
!
  implicit none

  character c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) idig
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) intval
  integer ( kind = 4 ) ipos
  integer ( kind = 4 ) ival
  character ( len = * )  s

  s = ' '

  ilo = 1
  ihi = len ( s )

  if ( ihi <= 0 ) then
    return
  end if
!
!  Make a copy of the integer.
!
  ival = intval
!
!  Handle the negative sign.
!
  if ( ival < 0 ) then

    if ( ihi <= 1 ) then
      s(1:1) = '*'
      return
    end if

    ival = -ival
    s(1:1) = '-'
    ilo = 2

  end if
!
!  The absolute value of the integer goes into S(ILO:IHI).
!
  ipos = ihi
!
!  Find the last digit of IVAL, strip it off, and stick it into the string.
!
  do

    idig = mod ( ival, 10 )
    ival = ival / 10

    if ( ipos < ilo ) then
      do i = 1, ihi
        s(i:i) = '*'
      end do
      return
    end if

    call digit_to_ch ( idig, c )
    s(ipos:ipos) = c
    ipos = ipos - 1

    if ( ival == 0 ) then
      exit
    end if

  end do
!
!  Shift the minus sign, if any.
!
  if ( s(1:1) == '-' ) then
    if ( ipos /= 1 ) then
      s(1:1) = ' '
      s(ipos:ipos) = '-'
    end if
  end if

  return
end
subroutine i4vec_to_b4set ( n_num, a_num, n, a )

!*****************************************************************************80
!
!! I4VEC_TO_B4SET converts an I4VEC to a B4SET.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 August 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N_NUM, the number of numeric entries.
!
!    Input, integer ( kind = 4 ) A_NUM(N_NUM), the numeric vector.
!    Entries of A_NUM should be between 1 and 32.
!
!    Input, integer ( kind = 4 ) N, the order of the master set.
!    N <= 32.
!
!    Output, integer ( kind = 4 ) A, the corresponding B4SET.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_num

  integer ( kind = 4 ) a
  integer ( kind = 4 ) a_num(n_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) pos
  integer ( kind = 4 ) pos_max

  a = 0
  pos_max = min ( bit_size ( a ), n )

  do i = 1, n_num
    pos = a_num(i)
    if ( 1 <= pos .and. pos <= pos_max ) then
      a = ibset ( a, pos - 1 )
    end if
  end do

  return
end
subroutine i4vec_to_lset ( n_num, a_num, n, a )

!*****************************************************************************80
!
!! I4VEC_TO_LSET converts an I4VEC to an LSET.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 August 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N_NUM, the number of numeric entries.
!
!    Input, integer ( kind = 4 ) A_NUM(N_NUM), the numeric vector.
!
!    Input, integer ( kind = 4 ) N, the order of the master set.
!
!    Output, logical A(N), the corresponding LSET.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_num

  logical a(n)
  integer ( kind = 4 ) a_num(n_num)
  integer ( kind = 4 ) i

  a(1:n) = .false.

  do i = 1, n_num
    call lset_insert ( n, a_num(i), a )
  end do

  return
end
subroutine lset_colex_rank ( n, t, rank )

!*****************************************************************************80
!
!! LSET_COLEX_RANK computes the colexicographic rank of an LSET.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 August 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the master set.
!
!    Input, logical T(N), the set.
!
!    Output, integer ( kind = 4 ) RANK, the rank of the set.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) rank
  logical t(n)

  rank = 0

  do i = 1, n

    if ( t(i) ) then
      rank = rank + 2 ** ( i - 1 )
    end if

  end do

  return
end
subroutine lset_colex_successor ( n, t, rank )

!*****************************************************************************80
!
!! LSET_COLEX_SUCCESSOR computes the colexicographic successor of an LSET.
!
!  Discussion:
!
!    In the original code, there is a last element with no successor.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 August 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the master set.
!
!    Input/output, logical T(N), describes a set.  
!    On input, T describes a set.
!    On output, T describes the next set in the ordering.
!    If the input T was the last in the ordering, then the output T
!    will be the first.
!
!    Input/output, integer ( kind = 4 ) RANK, the rank.
!    If RANK = -1 on input, then the routine understands that this is
!    the first call, and that the user wishes the routine to supply
!    the first element in the ordering, which has RANK = 0.
!    In general, the input value of RANK is increased by 1 for output,
!    unless the very last element of the ordering was input, in which
!    case the output value of RANK is 0.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) rank
  logical t(n)
!
!  Return the first element.
!
  if ( rank == -1 ) then
    t(1:n) = .false.
    rank = 0
    return
  end if

  do i = 1, n

    if ( .not. t(i) ) then
      t(i) = .true.
      rank = rank + 1
      return
    else
      t(i) = .false.
    end if

  end do

  rank = 0

  return
end
subroutine lset_colex_unrank ( rank, n, t )

!*****************************************************************************80
!
!! LSET_COLEX_UNRANK computes the LSET of given colexicographic rank.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 August 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) RANK, the rank of the set.
!
!    Input, integer ( kind = 4 ) N, the order of the master set.
!
!    Output, logical T(N), the set of the given rank.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) rank_copy
  integer ( kind = 4 ) sub_num
  logical t(n)
!
!  Check.
!
  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LSET_COLEX_UNRANK - Fatal error!'
    write ( *, '(a)' ) '  Input N is illegal.'
    stop
  end if

  call lset_enum ( n, sub_num )

  if ( rank < 0 .or. sub_num < rank ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LSET_COLEX_UNRANK - Fatal error!'
    write ( *, '(a)' ) '  The input rank is illegal.'
    stop
  end if

  rank_copy = rank

  do i = 1, n
    if ( mod ( rank_copy, 2 ) == 1 ) then
      t(i) = .true.
    else
      t(i) = .false.
    end if

    rank_copy = rank_copy / 2

  end do

  return
end
subroutine lset_complement ( n, a, b )

!*****************************************************************************80
!
!! LSET_COMPLEMENT computes the complement of an LSET.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 August 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the master set.
!
!    Input, logical A(N), the set.
!
!    Output, logical B(N), the complement of A.
!
  implicit none

  integer ( kind = 4 ) n

  logical a(n)
  logical b(n)

  b(1:n) = .not. a(1:n)

  return
end
subroutine lset_complement_relative ( n, a, b, c )

!*****************************************************************************80
!
!! LSET_COMPLEMENT_RELATIVE computes the relative complement of an LSET.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 August 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the master set.
!
!    Input, logical A(N), the set.
!
!    Input, logical B(N), the set with respect to which the complement is taken.
!
!    Output, logical C(N), the complement of A with respect to B.
!
  implicit none

  integer ( kind = 4 ) n

  logical a(n)
  logical b(n)
  logical c(n)

  c(1:n) = a(1:n) .and. ( .not. b(1:n) )

  return
end
subroutine lset_delete ( n, a, t )

!*****************************************************************************80
!
!! LSET_DELETE deletes an element from an LSET.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 August 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the master set.
!
!    Input, integer ( kind = 4 ) A, an item.
!
!    Input/output, logical T(N), a set.
!    On output, T(A) = FALSE.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a
  logical t(n)

  if ( a < 1 .or. n < a ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LSET_DELETE - Fatal error!'
    write ( *, '(a)' ) '  1 <= A <= N fails.'
    stop
  end if

  t(a) = .false.

  return
end
subroutine lset_distance ( n, t1, t2, dist )

!*****************************************************************************80
!
!! LSET_DISTANCE computes the Hamming distance between two LSET's.
!
!  Discussion:
!
!    The sets T1 and T2 are assumed to be subsets of a set of N elements.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 August 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the master set.
!
!    Input, logical T1(N), T2(N), two sets.
!
!    Output, integer ( kind = 4 ) DIST, the Hamming distance between T1 and T2,
!    defined as the number of elements of the master set which are
!    in either T1 or T2 but not both.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) dist
  integer ( kind = 4 ) i
  logical t1(n)
  logical t2(n)

  dist = 0

  do i = 1, n

    if ( (         t1(i)   .and. ( .not. t2(i) ) ) .or. &
         ( ( .not. t1(i) ) .and.         t2(i)   )  ) then
      dist = dist + 1
    end if

  end do

  return
end
subroutine lset_enum ( n, set_num )

!*****************************************************************************80
!
!! LSET_ENUM enumerates the LSET's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 August 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the master set.
!
!    Output, integer ( kind = 4 ) SET_NUM, the number of distinct sets.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) set_num

  set_num = 2**n

  return
end
function lset_index ( n, a, t )

!*****************************************************************************80
!
!! LSET_INDEX returns the index of an element of an LSET.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 August 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the master set.
!
!    Input, integer ( kind = 4 ) A, the item.
!
!    Input, logical T(N), a set.
!
!    Output, integer LSET_INDEX, the index of the item in the set,
!    or -1 if the item is not in the set.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a
  integer ( kind = 4 ) i
  integer ( kind = 4 ) lset_index
  logical t(n)
  integer ( kind = 4 ) value

  if ( a < 1 .or. n < a ) then
    value = -1
  else
    value = 0
    do i = 1, a
      if ( t(i) ) then
        value = value + 1
      end if
    end do
  end if

  lset_index = value

  return
end
subroutine lset_insert ( n, a, t )

!*****************************************************************************80
!
!! LSET_INSERT inserts an item into an LSET.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 August 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the master set.
!
!    Input, integer ( kind = 4 ) A, the item.
!    1 <= A <= N.
!
!    Input/output, logical T(N), a set.
!    On output, T(A) = TRUE.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a
  logical t(n)

  if ( a < 1 .or. n < a ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LSET_INSERT - Fatal error!'
    write ( *, '(a)' ) '  1 <= A <= N fails.'
    stop
  end if

  t(a) = .true.

  return
end
subroutine lset_intersect ( n, a, b, c )

!*****************************************************************************80
!
!! LSET_INTERSECT computes the intersection of two LSET's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 August 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the master set.
!
!    Input, logical A(N), B(N), two sets.
!
!    Output, logical C(N), the intersection of A and B.
!
  implicit none

  integer ( kind = 4 ) n

  logical a(n)
  logical b(n)
  logical c(n)

  c(1:n) = a(1:n) .and. b(1:n)

  return
end
function lset_is_empty ( n, t )

!*****************************************************************************80
!
!! LSET_IS_EMPTY determines if an LSET is empty.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 August 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the master set.
!
!    Input, logical T(N), a set.
!
!    Output, logical LSET_IS_EMPTY is TRUE if T is empty.
!
  implicit none

  integer ( kind = 4 ) n

  logical lset_is_empty
  logical t(n)

  lset_is_empty = all ( .not. t(1:n) )

  return
end
function lset_is_equal ( n, t1, t2 )

!*****************************************************************************80
!
!! LSET_IS_EQUAL determines if two LSET's are equal.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 August 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the master set.
!
!    Input, logical T1(N), T2(N), two sets.
!
!    Output, logical LSET_IS_EQUAL, is TRUE if T1 equals T2.
!
  implicit none

  integer ( kind = 4 ) n

  logical lset_is_equal
  logical t1(n)
  logical t2(n)

  lset_is_equal = all ( t1(1:n) .eqv. t2(1:n) )

  return
end
function lset_is_member ( n, a, t )

!*****************************************************************************80
!
!! LSET_IS_MEMBER determines if an item is a member of an LSET.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 August 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the master set.
!
!    Input, integer ( kind = 4 ) A, an item.
!
!    Input, logical T(N), a set.
!
!    Output, logical LSET_IS_MEMBER, is TRUE if A is an element of T.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a
  logical lset_is_member
  logical t(n)

  if ( 1 <= a .and. a <= n ) then
    lset_is_member = t(a)
  else
    lset_is_member = .false.
  end if

  return
end
function lset_is_subset ( n, t1, t2 )

!*****************************************************************************80
!
!! LSET_IS_SUBSET determines if one LSET is a subset of another.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 August 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the master set.
!
!    Input, logical T1(N), T2(N), two sets.
!
!    Output, logical LSET_IS_SUBSET, is TRUE if T1 is a subset of T2.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  logical lset_is_subset
  logical t1(n)
  logical t2(n)

  lset_is_subset = .true.

  do i = 1, n
    if ( t1(i) .and. .not. t2(i) ) then
      lset_is_subset = .false.
      return
    end if
  end do

  return
end
subroutine lset_lex_rank ( n, t, rank )

!*****************************************************************************80
!
!! LSET_LEX_RANK computes the lexicographic rank of an LSET.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 August 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the master set.
!
!    Input, logical T(N), the set.
!
!    Output, integer ( kind = 4 ) RANK, the rank of the set.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) rank
  logical t(n)

  rank = 0

  do i = 1, n

    if ( t(i) ) then
      rank = rank + 2**( n - i )
    end if

  end do

  return
end
subroutine lset_lex_successor ( n, t, rank )

!*****************************************************************************80
!
!! LSET_LEX_SUCCESSOR computes the lexicographic successor of an LSET.
!
!  Discussion:
!
!    In the original code, there is a last element with no successor.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 August 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the master set.
!
!    Input/output, logical T(N), describes a set.
!    On input, T describes a set.
!    On output, T describes the next set in the ordering.
!    If the input T was the last in the ordering, then the output T
!    will be the first.
!
!    Input/output, integer ( kind = 4 ) RANK, the rank.
!    If RANK = -1 on input, then the routine understands that this is
!    the first call, and that the user wishes the routine to supply
!    the first element in the ordering, which has RANK = 0.
!    In general, the input value of RANK is increased by 1 for output,
!    unless the very last element of the ordering was input, in which
!    case the output value of RANK is 0.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) rank
  logical t(n)
!
!  Return the first element.
!
  if ( rank == -1 ) then
    t(1:n) = .false.
    rank = 0
    return
  end if

  do i = n, 1, -1

    if ( .not. t(i) ) then
      t(i) = .true.
      rank = rank + 1
      return
    else
      t(i) = .false.
    end if

  end do

  rank = 0

  return
end
subroutine lset_lex_unrank ( rank, n, t )

!*****************************************************************************80
!
!! LSET_LEX_UNRANK computes the LSET of given lexicographic rank.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 August 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) RANK, the rank of the set.
!
!    Input, integer ( kind = 4 ) N, the order of the master set.
!
!    Output, logical T(N), the set of the given rank.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) nsub
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) rank_copy
  integer ( kind = 4 ) set_num
  logical t(n)
!
!  Check.
!
  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LSET_LEX_UNRANK - Fatal error!'
    write ( *, '(a)' ) '  Input N is illegal.'
    stop
  end if

  call lset_enum ( n, set_num )

  if ( rank < 0 .or. set_num < rank ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LSET_LEX_UNRANK - Fatal error!'
    write ( *, '(a)' ) '  The input rank is illegal.'
    stop
  end if

  rank_copy = rank

  do i = n, 1, -1

    if ( mod ( rank_copy, 2 ) == 1 ) then
      t(i) = .true.
    else
      t(i) = .false.
    end if

    rank_copy = rank_copy / 2

  end do

  return
end
subroutine lset_transpose_print ( n, t, title )

!*****************************************************************************80
!
!! LSET_TRANSPOSE_PRINT prints an LSET "transposed".
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 August 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the master set.
!
!    Input, logical T(N), the set.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  logical lset_is_empty
  integer ( kind = 4 ) s
  logical t(n)
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  if ( lset_is_empty ( n, t ) ) then
    write ( *, '(a)' ) '  (Empty set)'
  else
    s = 0
    do i = 1, n
      if ( t(i) ) then
        write ( *, '(2x,i2)', advance = 'no' ) i
        s = s + 4
      end if
      if ( 76 < s .or. ( 0 < s .and. i == n ) ) then
        write ( *, '(1x)', advance = 'yes' )
        s = 0
      end if
    end do
  end if

  return
end
subroutine lset_random ( n, seed, a )

!*****************************************************************************80
!
!! LSET_RANDOM sets a rondom LSET.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 August 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the master set.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, logical A(N).
!
  implicit none

  integer ( kind = 4 ) n

  logical a(n)
  integer ( kind = 4 ), parameter :: i4_huge      = 2147483647
  integer ( kind = 4 ), parameter :: i4_huge_half = 1073741823
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LSET_RANDOM - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + i4_huge
    end if

    a(i) = ( i4_huge_half < seed )

  end do

  return
end
subroutine lset_to_b4set ( n, a_log, a )

!*****************************************************************************80
!
!! LSET_TO_B4SET converts an I4VEC to a B4SET.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 September 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the master set.
!    N <= 32.
!
!    Input, logical A_LOG(N), the logical representation of the set.
!
!    Output, integer ( kind = 4 ) A, the corresponding B4SET.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a
  logical a_log(n_num)
  integer ( kind = 4 ) i

  a = 0

  do i = 1, n
    if ( a_log(i) ) then
      a = ibset ( a, i - 1 )
    end if
  end do

  return
end
subroutine lset_union ( n, a, b, c )

!*****************************************************************************80
!
!! LSET_UNION computes the union of two LSET's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 August 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the master set.
!
!    Input, logical A(N), B(N), two sets.
!
!    Output, logical C(N), the union of A and B.
!
  implicit none

  integer ( kind = 4 ) n

  logical a(n)
  logical b(n)
  logical c(n)

  c(1:n) = a(1:n) .or. b(1:n)

  return
end
subroutine lset_weight ( n, t, weight )

!*****************************************************************************80
!
!! LSET_WEIGHT computes the Hamming weight of an LSET.
!
!  Discussion:
!
!    The Hamming weight is simply the number of elements in the set.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 August 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the master set.
!
!    Input, logical T(N), the set.
!
!    Output, integer ( kind = 4 ) WEIGHT, the Hamming weight of the set T.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  logical t(n)
  integer ( kind = 4 ) weight

  weight = 0

  do i = 1, n
    if ( t(i) ) then
      weight = weight + 1
    end if
  end do

  return
end
subroutine lset_xor ( n, a, b, c )

!*****************************************************************************80
!
!! LSET_XOR computes the symmetric difference of two LSET's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 August 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the master set.
!
!    Input, logical A(N), B(N), two sets.
!
!    Output, logical C(N), the symmetric difference of A and B.
!
  implicit none

  integer ( kind = 4 ) n

  logical a(n)
  logical b(n)
  logical c(n)

  c(1:n) =  (         a(1:n)   .and. ( .not. b(1:n) ) ) .or. &
            ( ( .not. a(1:n) ) .and.         b(1:n)   )

  return
end
subroutine lvec_transpose_print ( n, t, title )

!*****************************************************************************80
!
!! LVEC_TRANSPOSE_PRINT prints an LVEC "transposed".
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 August 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the master set.
!
!    Input, logical T(N), the set.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  logical t(n)
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do ilo = 1, n, 80
    ihi = min ( ilo + 80 - 1, n )
    write ( *, '(80l1)' ) t(ilo:ihi)
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

  character ( len = 8 )  ampm
  integer   ( kind = 4 ) d
  character ( len = 8 )  date
  integer   ( kind = 4 ) h
  integer   ( kind = 4 ) m
  integer   ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer   ( kind = 4 ) n
  integer   ( kind = 4 ) s
  character ( len = 10 ) time
  integer   ( kind = 4 ) values(8)
  integer   ( kind = 4 ) y
  character ( len = 5 )  zone

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
