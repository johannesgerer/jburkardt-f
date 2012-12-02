subroutine partition_brute ( n, w, c, discrepancy )

!*****************************************************************************80
!
!! PARTITION_BRUTE approaches the partition problem using brute force.
!
!  Discussion:
!
!    We are given a set of N integers W.
!
!    We seek to partition W into subsets W0 and W1, such that the subsets
!    have equal sums.
!
!    The "discrepancy" is the absolute value of the difference between the
!    two sums, and will be zero if we have solved the problem.
!
!    For a given set of integers, there may be zero, one, or many solutions.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 May 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the set.
!
!    Input, integer ( kind = 4 ) W(N), the integers.
!
!    Output, integer ( kind = 4 ) C(N), indicates the proposed solution.
!    C(I) is 0 for items in set W0 and 1 for items in set W1.
!
!    Output, integer ( kind = 4 ) DISCREPANCY, the discrepancy.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) c(n)
  integer ( kind = 4 ) d(n)
  integer ( kind = 4 ) d_discrepancy
  integer ( kind = 4 ) discrepancy
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) w(n)
  integer ( kind = 4 ) w_sum

  w_sum = sum ( w(1:n) )
  discrepancy = w_sum

  rank = -1

  do

    call subset_next ( n, d, rank )

    if ( rank == -1 ) then
      exit
    end if

    d_discrepancy = abs ( w_sum - 2 * dot_product ( d, w ) )

    if ( d_discrepancy < discrepancy ) then
      discrepancy = d_discrepancy
      c(1:n) = d(1:n)
    end if

    if ( discrepancy == 0 ) then
      exit
    end if

  end do

  return
end
subroutine partition_count ( n, w, count )

!*****************************************************************************80
!
!! PARTITION_COUNT counts the solutions to a partition problem.
!
!  Discussion:
!
!    We are given a set of N integers W.
!
!    We seek to partition W into subsets W0 and W1, such that the subsets
!    have equal sums.
!
!    The "discrepancy" is the absolute value of the difference between the
!    two sums, and will be zero if we have solved the problem.
!
!    For a given set of integers, there may be zero, one, or many solutions.
!
!    In the case where the weights are distinct, the count returned by this
!    function may be regarded as twice as big as it should be, since the
!    partition (W0,W1) is counted a second time as (W1,W0).  A more serious
!    overcount can occur if the set W contains duplicate elements - in the
!    extreme case, W might be entirely 1's, in which case there is really
!    only one (interesting) solution, but this function will count many.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 May 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the set.
!
!    Input, integer ( kind = 4 ) W(N), the integers.
!
!    Output, integer ( kind = 4 ) COUNT, the number of solutions.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) c(n)
  integer ( kind = 4 ) count
  integer ( kind = 4 ) discrepancy
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) w(n)
  integer ( kind = 4 ) w_sum

  w_sum = sum ( w(1:n) )

  rank = -1
  count = 0

  do

    call subset_next ( n, c, rank )

    if ( rank == -1 ) then
      exit
    end if

    discrepancy = abs ( w_sum - 2 * dot_product ( c, w ) )

    if ( discrepancy == 0 ) then
      count = count + 1
    end if

  end do

  return
end
subroutine subset_next ( n, t, rank )

!*****************************************************************************80
!
!! SUBSET_NEXT computes the subset lexicographic successor.
!
!  Discussion:
!
!    This is a lightly modified version of "subset_lex_successor()" from COMBO.
!
!  Example:
!
!    On initial call, N is 5 and the input value of RANK is -1.
!    Then here are the successive outputs from the program:
!
!   Rank   T1   T2   T3   T4   T5
!   ----   --   --   --   --   --
!      0    0    0    0    0    0
!      1    0    0    0    0    1
!      2    0    0    0    1    0
!      3    0    0    0    1    1
!     ..   ..   ..   ..   ..   ..
!     30    1    1    1    1    0
!     31    1    1    1    1    1
!     -1    0    0    0    0    0  <-- Reached end of cycle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 May 2012
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
!    Input, integer ( kind = 4 ) N, the number of elements in the master set.
!    N must be positive.
!
!    Input/output, integer ( kind = 4 ) T(N), describes a subset.  T(I) is 0 if
!    the I-th element of the master set is not in the subset, and is
!    1 if the I-th element is part of the subset.
!    On input, T describes a subset.
!    On output, T describes the next subset in the ordering.
!
!    Input/output, integer ( kind = 4 ) RANK, the rank.
!    If RANK = -1 on input, then the routine understands that this is
!    the first call, and that the user wishes the routine to supply
!    the first element in the ordering, which has RANK = 0.
!    In general, the input value of RANK is increased by 1 for output,
!    unless the very last element of the ordering was input, in which
!    case the output value of RANK is -1.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) t(n)
!
!  Return the first element.
!
  if ( rank == -1 ) then
    t(1:n) = 0
    rank = 0
    return
  end if

  do i = n, 1, -1

    if ( t(i) == 0 ) then
      t(i) = 1
      rank = rank + 1
      return
    else
      t(i) = 0
    end if

  end do

  rank = -1

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
