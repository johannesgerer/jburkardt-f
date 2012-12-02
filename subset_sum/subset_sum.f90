subroutine i4_to_digits_binary ( i, n, c )

!*****************************************************************************80
!
!! I4_TO_DIGITS_BINARY produces the binary digits of an I4.
!
!  Discussion:
!
!    An I4 is an integer.
!
!  Example:
!
!     I    N     C               Binary
!    --  ---   ---         ------------
!     0    1   0                      0
!     0    2   0, 0                  00
!     1    3   1, 0, 0              100
!     2    3   0, 1, 0              010
!     3    3   1, 1, 0              011
!     4    3   0, 0, 1              100
!     8    3   0, 0, 0           (1)000
!     8    5   0, 0, 0, 1, 0      01000
!    -8    5   0, 0, 0, 1, 0  (-) 01000
!
!     0    3   0, 0, 0
!     1    3   1, 0, 0
!     2    3   0, 1, 0
!     3    3   1, 1, 0
!     4    3   0, 0, 1
!     5    3   1, 0, 1
!     6    3   0, 1, 1
!     7    3   1, 1, 1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 December 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, an integer to be represented.
!
!    Input, integer ( kind = 4 ) N, the number of binary digits to produce.
!
!    Output, integer ( kind = 4 ) C(N), the first N binary digits of I,
!    with C(1) being the units digit.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) c(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_copy
  integer ( kind = 4 ) j

  i_copy = abs ( i )

  do j = 1, n

    c(j) = mod ( i_copy, 2 )
    i_copy = i_copy / 2

  end do

  return
end
subroutine subset_sum_count ( n, w, t, ind_min, ind_max, solution_num )

!*****************************************************************************80
!
!! SUBSET_SUM_COUNT counts solutions to the subset sum problem in a range.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 December 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the set.
!
!    Input, integer ( kind = 4 ) W(N), a set of weights.  The length of this
!    array must be no more than 31.
!
!    Input, integer ( kind = 4 ) T, the target value.
!
!    Input, integer ( kind = 4 ) IND_MIN, IND_MAX, the lower and upper
!    limits to be searched.  0 <= IND_MIN <= IND_MAX <= (2^N)-1.
!
!    Output, integer ( kind = 4 ) SOLUTION_NUM, the number of distinct
!    solutions of the subset sum problem found within the given range.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) c(n)
  integer ( kind = 4 ) ind
  integer ( kind = 4 ) ind_max
  integer ( kind = 4 ) ind_max2
  integer ( kind = 4 ) ind_min
  integer ( kind = 4 ) ind_min2
  integer ( kind = 4 ) solution_num
  integer ( kind = 4 ) t
  integer ( kind = 4 ) w(n)
!
!  Check the data.
!
  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SUBSET_SUM_COUNT - Fatal error!'
    write ( *, '(a)' ) '  N < 1.'
    stop
  end if

  if ( 31 < n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SUBSET_SUM_COUNT - Fatal error!'
    write ( *, '(a)' ) '  31 < N.'
    stop
  end if

  ind_min2 = max ( ind_min, 0 )
  ind_max2 = min ( ind_max, ( 2 ** n ) - 1 )
!
!  Run through the range.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Searching from IND_MIN = ', ind_min2
  write ( *, '(a,i8)' ) '  through IND_MAX = ', ind_max2

  solution_num = 0

  do ind = ind_min2, ind_max2
!
!  Convert INDEX into vector of indices in W.
!
    call i4_to_digits_binary ( ind, n, c )
!
!  If the sum of those weights matches the target, return combination.
!
    if ( dot_product ( c, w ) == t ) then
      solution_num = solution_num + 1
    end if

  end do

  return
end
subroutine subset_sum_find ( n, w, t, ind_min, ind_max, ind, c )

!*****************************************************************************80
!
!! SUBSET_SUM seeks a subset of a set that has a given sum.
!
!  Discussion:
!
!    This function tries to compute a target value as the sum of
!    a selected subset of a given set of weights.
!
!    This function works by brute force, that is, it tries every
!    possible subset to see if it sums to the desired value.
!
!    Given N weights, every possible selection can be described by
!    one of the N-digit binary numbers from 0 to 2^N-1.
!
!    This function includes a range, which allows the user to
!    control which subsets are to be checked.  Thus, if there are
!    N weights, specifying a range of [ 0, 2^N-1] indicates that
!    all subsets should be checked.  On the other hand, this full
!    range could be broken down into smaller subranges, each of
!    which could be checked independently.
!
!    It is possible that, in the given range, there may be multiple
!    solutions of the problem.  This function will only return
!    one such solution, if found.  However, the function may be called
!    again, with an appropriate restriction of the range, to continue
!    the search for other solutions.
!
!  Example:
!
!    w = [ 1, 2, 4, 8, 16, 32 ];
!    t = 22;
!    r = [ 0, 2^6 - 1 ];
!
!    call subset_sum ( w, t, r, c, ind )
!
!    c = [ 2, 3, 5 ]
!    index = 22
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 December 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the set.
!
!    Input, integer ( kind = 4 ) W(N), a set of weights.  The length of this
!    array must be no more than 31.
!
!    Input, integer ( kind = 4 ) T, the target value.
!
!    Input, integer ( kind = 4 ) IND_MIN, IND_MAX, the lower and upper
!    limits to be searched.  0 <= IND_MIN <= IND_MAX <= (2^N)-1.
!
!    Output, integer ( kind = 4 ) IND, the index of the solution.
!    If IND is -1, no solution was found in the range.
!
!    Output, integer ( kind = 4 ) C(N), indicates the solution, assuming
!    that IND is not -1.  In that case, the sum T is made by selecting
!    those weights W(I) for which C(I) is 1.  In fact,
!    T = sum ( 1 <= I <= N ) C(I) * W(I).
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) c(n)
  integer ( kind = 4 ) ind
  integer ( kind = 4 ) ind_max
  integer ( kind = 4 ) ind_max2
  integer ( kind = 4 ) ind_min
  integer ( kind = 4 ) ind_min2
  integer ( kind = 4 ) t
  integer ( kind = 4 ) w(n)
!
!  Check the data.
!
  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SUBSET_SUM - Fatal error!'
    write ( *, '(a)' ) '  N < 1.'
    stop
  end if

  if ( 31 < n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SUBSET_SUM - Fatal error!'
    write ( *, '(a)' ) '  31 < N.'
    stop
  end if

  ind_min2 = max ( ind_min, 0 )
  ind_max2 = min ( ind_max, ( 2 ** n ) - 1 )
!
!  Run through the range.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Searching from IND_MIN = ', ind_min2
  write ( *, '(a,i8)' ) '  through IND_MAX = ', ind_max2

  do ind = ind_min2, ind_max2
!
!  Convert INDEX into vector of indices in W.
!
    call i4_to_digits_binary ( ind, n, c )
!
!  If the sum of those weights matches the target, return combination.
!
    if ( dot_product ( c, w ) == t ) then
      return
    end if

  end do

  ind = - 1

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
