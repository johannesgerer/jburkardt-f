subroutine bisection_integer ( f, a, b, c, fc )

!*****************************************************************************80
!
!! BISECTION_INTEGER seeks an integer root using bisection.
!
!  Discussion:
!
!    A function F(X) confined to integer arguments is given, with an
!    interval [A,B] over which F changes sign.  An integer C is sought
!    such that A <= C <= B and F(C) = 0.
!
!    Because we are restricted to integer arguments, it may the case that
!    there is no such C.
!
!    This routine proceeds by a form of bisection, in which the enclosing
!    interval is restricted to be defined by integer values.
!
!    If the user has given a true change of sign interval [A,B], and if,
!    in the interval, there is a single integer value C for which F(C) = 0,
!    with the additional restrictions that F(C-1) and F(C+1) are of opposite
!    signs, then this procedure should locate and return C.
!
!    In particular, if the function F is monotone, and there is an integer
!    solution C in the interval, then this procedure will find it.
!
!    However, in general, even if there is an integer C in the interval,
!    such that F(C) = 0, this procedure may be unable to find it, particularly
!    if there are also nonintegral solutions within the same interval.
!
!    While any integer function can be used with this program, the bisection
!    approach is most useful if the integer function is monotone, or
!    varies slowly, or can be regarded as the restriction to integer arguments
!    of a continuous (and smoothly varying) function of a real argument.
!    In such cases, knowing that F is negative at A and positive at B
!    suggests that F generally increases from A to B, and might attain 
!    the value 0 at some intermediate argument C.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 May 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, external integer ( kind = 4 ) F, the name of a user-supplied 
!    procedure that evaluates the function, of the form
!      function f ( c )
!      integer ( kind = 4 ) c, f
!
!   Input, integer ( kind = 4 ) A, B, two arguments that define a change of
!   sign interval for F.  In other words, F(A) and F(B) must be of opposite
!   sign.
!
!   Output, integer ( kind = 4 ) C, FC, the candidate for the root, as 
!   determined by the program, and its function value.  If FC is not zero,
!   then the procedure did not find a root in the interval, and C is only
!   an "approximate" root.
!
  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) c
  integer ( kind = 4 ), external :: f
  integer ( kind = 4 ) fa
  integer ( kind = 4 ) fb
  integer ( kind = 4 ) fc
  integer ( kind = 4 ) t
!
!  Ensure that F(A) < 0 < F(B).
!
  fa = f(a)
  fb = f(b)

  if ( fa == 0 ) then
    c = a
    fc = fa
  else if ( fb == 0 ) then
    c = b
    fc = fb
  else if ( fa < 0 .and. 0 < fb ) then

  else if ( fb < 0 .and. 0 < fa ) then
    t = a
    a = b
    b = t
    t = fa
    fa = fb
    fb = t
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BISECTION_INTEGER - Fatal error!'
    write ( *, '(a)' ) '  No change of sign interval supplied.'
    write ( *, '(a,i8,a,i8)' ) '  F(', a, ') = ', fa
    write ( *, '(a,i8,a,i8)' ) '  F(', b, ') = ', fb
    stop
  end if
!
!  Bisection.
!
  do while ( 1 < abs ( b - a ) )

    c = ( a + b ) / 2
    fc = f(c)

    if ( fc == 0 ) then
      return
    else if ( fc < 0 ) then
      a = c
      fa = fc
    else if ( 0 < fc ) then
      b = c
      fb = fc
    end if
    
  end do
!
!  Interval is empty, with FA < 0 and 0 < FB.
!  Bisection did not produce an integer solution.
!  Return the argument with smallest function norm.
!
  if ( - fa < fb ) then
    c = a
    fc = fa
  else
    c = b
    fc = fb
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
