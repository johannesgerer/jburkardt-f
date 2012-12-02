subroutine find_distances ( l_length, l, x_length, x, y, success )

!*****************************************************************************80
!
!! FIND_DISTANCES determines if the "free" distances include every ||X(I)-Y||.
!
!  Discussion:
!
!    This routine is given a candidate point Y, a set of placed points
!    X(1:X_LENGTH), and a list of unused or "free" distances in
!    L(1:L_LENGTH).  The routine seeks to find in L a copy of the
!    distance from Y to each X.
!
!    If so, then the L array is reordered so that entries
!    L(L_LENGTH-X_LENGTH+1:L_LENGTH) contain theses distances.
!
!    In other words, Y can be added into X, and L_LENGTH reduced to
!    L_LENGTH-X_LENGTH.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 November 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Pavel Pevzner,
!    Computational Molecular Biology,
!    MIT Press, 2000,
!    ISBN: 0-262-16197-4,
!    LC: QH506.P47.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) L_LENGTH, the length of the array.
!
!    Input/output, integer ( kind = 4 ) L(L_LENGTH), the array.  On output,
!    some entries have been shuffled.  In particular, if SUCCESS is TRUE,
!    the entries L(L_LENGTH-X_LENGTH+1:L_LENGTH) contain the distances
!    of X(1:X_LENGTH) to Y.
!
!    Input, integer ( kind = 4 ) X_LENGTH, the number of entries in X.
!
!    Input, integer ( kind = 4 ) X(X_LENGTH), the number of points
!    already accepted.
!
!    Input, integer ( kind = 4 ) Y, a new point that we are considering.
!
!    Output, logical SUCCESS, is TRUE if the entries of L included
!    the values of the distance of Y to each entry of X.
!
  implicit none

  integer ( kind = 4 ) l_length
  integer ( kind = 4 ) x_length

  integer ( kind = 4 ) d
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l(l_length)
  integer ( kind = 4 ) l2_length
  logical success
  integer ( kind = 4 ) x(x_length)
  integer ( kind = 4 ) y

  l2_length = l_length

  do i = 1, x_length

    d = abs ( x(i) - y )

    success = .false.

    do j = 1, l2_length

      if ( l(j) == d ) then
        l(j) = l(l2_length)
        l(l2_length) = d
        l2_length = l2_length - 1
        success = .true.
        exit
      end if

    end do

    if ( .not. success ) then
      return
    end if

  end do

  success = .true.

  return
end
function i4vec_max_delete ( l_length, l )

!*****************************************************************************80
!
!! I4VEC_MAX_DELETE deletes the maximum entry from an I4VEC.
!
!  Discussion:
!
!    This routine finds the largest entry in an array and moves
!    it to the end of the array.
!
!    If we ignore this last array entry, then the effect is the same
!    as "deleting" the maximum entry from the array.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 November 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Pavel Pevzner,
!    Computational Molecular Biology,
!    MIT Press, 2000,
!    ISBN: 0-262-16197-4,
!    LC: QH506.P47.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) L_LENGTH, the length of the array.
!
!    Input/output, integer ( kind = 4 ) L(L_LENGTH), the array.  On output,
!    the maximum entry has been "deleted", that is, the array has
!    been shifted so that this entry occurs at the end.
!
!    Output, integer ( kind = 4 ) I4VEC_MAX_DELETE, the maximum entry in the
!    input array.
!
  implicit none

  integer ( kind = 4 ) l_length

  integer ( kind = 4 ) i4vec_max_delete
  integer ( kind = 4 ) l(l_length)
  integer ( kind = 4 ) max_index(1)
  integer ( kind = 4 ) value

  max_index = maxloc ( l(1:l_length) )
  value = l(max_index(1))

  l(max_index(1):l_length-1) = l(max_index(1)+1:l_length)
  l(l_length) = value

  i4vec_max_delete = value

  return
end
subroutine i4vec_print ( n, a, title )

!*****************************************************************************80
!
!! I4VEC_PRINT prints an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
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
  integer ( kind = 4 ) i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i8,2x,i12)' ) i, a(i)
  end do

  return
end
recursive subroutine place ( l_length, l, x_length, x )

!*****************************************************************************80
!
!! PLACE tries to place the next point for the partial digest problem.
!
!  Discussion:
!
!    Note that this is a recursive subroutine.  A solution to the
!    partial digest problem is sought by calling this routine repeatedly.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 November 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Pavel Pevzner,
!    Computational Molecular Biology,
!    MIT Press, 2000,
!    ISBN: 0-262-16197-4,
!    LC: QH506.P47.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) L_LENGTH, the number of entries in L.
!    On return, L_LENGTH has been reduced by 1.
!
!    Input/output, integer ( kind = 4 ) L(L_LENGTH)...
!
!    Input, integer ( kind = 4 ) X_LENGTH, the number of entries in X.
!
!    Input, integer ( kind = 4 ) X(X_LENGTH), the current partial solution.
!
  implicit none

  integer ( kind = 4 ) l_length
  integer ( kind = 4 ) x_length

  integer ( kind = 4 ) i4vec_max_delete
  integer ( kind = 4 ) l(l_length)
  integer ( kind = 4 ) l_length2
  logical success
  integer ( kind = 4 ) x(x_length)
  integer ( kind = 4 ) y
!
!  Are we done?
!
  if ( l_length <= 0 ) then
    call i4vec_print ( x_length, x, '  Solution:' )
    return
  end if
!
!  Find the maximum remaining distance.
!
  y = i4vec_max_delete ( l_length, l )
!
!  We can add a point at Y if L contains all the distances from Y to
!  the current X's.
!
  call find_distances ( l_length, l, x_length, x, y, success )

  if ( success ) then
    l_length2 = l_length - x_length
    x_length = x_length + 1
    x(x_length) = y
    call place ( l_length2, l, x_length, x )
    x_length = x_length - 1
  end if
!
!  We must also consider the case where Y represents the distance
!  to X(2), not X(1).
!
  y = x(2) - y

  call find_distances ( l_length, l, x_length, x, y, success )

  if ( success ) then
    l_length2 = l_length - x_length
    x_length = x_length + 1
    x(x_length) = y
    call place ( l_length2, l, x_length, x )
    x_length = x_length - 1
  end if

  return
end
subroutine partial_digest_recur ( n, l )

!*****************************************************************************80
!
!! PARTIAL_DIGEST_RECUR uses recursion on the partial digest problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 November 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Pavel Pevzner,
!    Computational Molecular Biology,
!    MIT Press, 2000,
!    ISBN: 0-262-16197-4,
!    LC: QH506.P47.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of nodes.
!
!    Input, integer ( kind = 4 ) L((N*(N-1))/2), the distances between all pairs
!    of distinct nodes.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i4vec_max_delete
  integer ( kind = 4 ) l((n*(n-1))/2)
  integer ( kind = 4 ) l_length
  integer ( kind = 4 ) width
  integer ( kind = 4 ) x(n)
  integer ( kind = 4 ) x_length

  l_length = (n*(n-1))/2
  width = i4vec_max_delete ( l_length, l )
  l_length = l_length - 1

  x(1) = 0
  x(2) = width
  x_length = 2
  call place ( l_length, l, x_length, x )

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
