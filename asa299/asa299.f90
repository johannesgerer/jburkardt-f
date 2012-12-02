subroutine simplex_lattice_point_next ( n, t, more, x )

!*****************************************************************************80
!
!! SIMPLEX_LATTICE_POINT_NEXT generates lattice points in a simplex.
!
!  Discussion:
!
!    The simplex is defined by N-dimensional points X such that:
!
!        0 <= X(1:N)
!
!    and
!
!      sum ( X(1:N) ) <= T
!
!    where T is an integer.
!
!    Lattice points are points X which satisfy the simplex conditions and
!    for which all the components are integers.
!
!    This routine generates all the lattice points in a given simplex, one at 
!    a time, in a reverse lexicographic order.
!
!    To use the routine, initialize by setting N and T to appropriate values, 
!    and MORE to FALSE.  The initial value of X is not important.
!
!    Call the routine. On return, X will contain the first lattice point in 
!    the simplex.  If MORE is TRUE, then the routine may be called again to 
!    get the next point.  In fact, as long as the output value of MORE is 
!    TRUE, there is at least one more lattice point that can be found by 
!    making another call.  When MORE is returned as FALSE, then there are no 
!    more lattice points; the value of X returned at that time is the 
!    "last" such point.
!
!    During the computation of a sequence of lattice points, the user should 
!    not change the values of N, T, MORE or X.  
!
!    The output for N = 3, T = 4 would be:
!
!       1    4  0  0
!       2    3  1  0
!       3    3  0  1
!       4    2  2  0
!       5    2  1  1
!       6    2  0  2
!       7    1  3  0
!       8    1  2  1
!       9    1  1  2
!      10    1  0  3
!      11    0  4  0
!      12    0  3  1
!      13    0  2  2
!      14    0  1  3
!      15    0  0  4
!
!  Modified:
!
!    10 January 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Scott Chasalow, Richard Brand,
!    Algorithm AS 299:
!    Generation of Simplex Lattice Points,
!    Applied Statistics,
!    Volume 44, Number 4, 1995, pages 534-545.
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
!    Input, integer N, the spatial dimension.
!    N must be positive.
!
!    Input, integer T, the characteristic of the simplex.
!    T must be nonnegative.
!
!    Input/output, logical MORE, initialized to FALSE by the user to
!    begin a sequence of calculations, returned by the routine as TRUE,
!    if there are more values of X that can be calculated, or FALSE
!    if the accompanying value of X is the last one for this sequence.
!
!    Input/output, integer X(N), not initialized by the user, but not
!    changed by the user on subsequent calls.  The routine returns
!    a new point on each call, and on subsequent calls uses the input
!    value (old point) to compute the output value (next point).
!
  implicit none

  integer n

  integer i
  integer j
  logical more
  integer t
  integer x(n)

  if ( .not. more ) then

    if ( n < 1 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SIMPLEX_LATTICE_POINT_NEXT - Fatal error!'
      write ( *, '(a)' ) '  N < 1.'
      stop
    end if

    if ( t < 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SIMPLEX_LATTICE_POINT_NEXT - Fatal error!'
      write ( *, '(a)' ) '  T < 0.'
      stop
    end if

    more = .true.
    j = 1

    x(1) = t
    x(2:n) = 0
!
!  The first point can actually also be the last!
!
    if ( n == 1 ) then
      more = .false.
    end if

  else
!
!  Search X(N-1 down to 1) for the first nonzero element.
!  If none, then terminate.  (This should not happen!)
!  Otherwise, set J to this index.
!  Decrement X(J) by 1.
!  Set X(J+1:N) to (T-X(1:J),0,0,...0).
!
    j = n

    do i = n-1, 1, -1

      if ( 0 < x(i) ) then
        j = i
        exit
      end if

    end do

    if ( j == n ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SIMPLEX_LATTICE_POINT_NEXT - Fatal error!'
      write ( *, '(a)' ) '  The input X vector is nonpositive in all entries'
      write ( *, '(a)' ) '  except possibly the last one.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Perhaps the user has miscalled the routine'
      write ( *, '(a)' ) '  or altered data between calls.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ABNORMAL TERMINATION.'
      stop
    end if

    x(j) = x(j) - 1
    x(j+1) = t - sum ( x(1:j) )
    x(j+2:n) = 0
!
!  Is this the last point?
!
    if ( x(n) == t ) then
      more = .false.
    end if

  end if

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
  integer d
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer n
  integer s
  integer values(8)
  integer y

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
