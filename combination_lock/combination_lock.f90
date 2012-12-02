subroutine bicycle_lock ( c, step )

!*****************************************************************************80
!
!! BICYCLE_LOCK finds the combination on a typical bicycle lock.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 May 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) C, the combination, a value between 0 and 999.
!
!    Output, integer ( kind = 4 ) STEP, the step on which the combination 
!    was found.  A value of -1 means the combination was not found.
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) c
  integer ( kind = 4 ) step

  step = -1

  do a = 0, 999

    if ( a == c ) then
      step = a + 1
      exit
    end if
  
  end do

  return
end
subroutine combination_lock ( m, n, c, step )

!*****************************************************************************80
!
!! COMBINATION_LOCK determines the combination of a lock.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 May 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of dials.
!
!    Input, integer ( kind = 4 ) N, the number of symbols on each dial.
!    We assume the symbols are the integers 0 to N-1.
!
!    Input, integer ( kind = 4 ) C(M), the combination.
!
!    Output, integer ( kind = 4 ) STEP, the step on which the combination 
!    was found.  A value of -1 means the combination was not found.
!
  implicit none

  integer ( kind = 4 ) m

  integer ( kind = 4 ) a(m)
  integer ( kind = 4 ) c(m)
  logical i4vec_eq
  logical more
  integer ( kind = 4 ) n
  integer ( kind = 4 ) step
!
!  Starting with the guess (0, 0, ... 0),
!  generate every possible combination, in order, and try it.
!
  more = .false.
  a(1:m) = 0
  step = 0

  do
  
    call combination_next ( m, n, a, more )

    if ( .not. more ) then
      step = -1
      exit
    end if

    step = step + 1

    if ( i4vec_eq ( m, a, c ) ) then
      exit
    end if
  
  end do

  return
end
subroutine combination_next ( m, base, a, more )

!*****************************************************************************80
!
!! COMBINATION_NEXT generates lock combinations in lex order.
!
!  Discussion:
!
!    The vectors are produced in lexical order, starting with
!    (0,0,...,0),
!    (0,0,...,1),
!    ...
!    (BASE-1,BASE-1,...,BASE-1).
!
!  Example:
!
!    M = 2,
!    BASE = 3
!
!    0   0
!    0   1
!    0   2
!    1   0
!    1   1
!    1   2
!    2   0
!    2   1
!    2   2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 May 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Dennis Stanton, Dennis White,
!    Constructive Combinatorics,
!    Springer, 1986,
!    ISBN: 0387963472,
!    LC: QA164.S79.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the size of the vectors to be used.
!
!    Input, integer ( kind = 4 ) BASE, the base to be used.  BASE = 2 will
!    give vectors of 0's and 1's, for instance.
!
!    Input/output, integer ( kind = 4 ) A(M).  The input value of A is
!    not important on the first call.  Thereafter, it should simply be the 
!    output value from the previous call.  The output value is the next vector
!    in the sequence.
!
!    Input/output, logical MORE.  The input value should be FALSE on the first 
!    call, and TRUE on subsequent calls.  The output value will be TRUE as long 
!    as the next vector could be computed, and FALSE once there are no more.
!
  implicit none

  integer ( kind = 4 ) m

  integer ( kind = 4 ) a(m)
  integer ( kind = 4 ) base
  integer ( kind = 4 ) i
  logical more

  if ( .not. more ) then

    a(1:m) = 0
    more = .true.

  else
      
    do i = m, 1, -1

      a(i) = a(i) + 1

      if ( a(i) < base ) then
        return
      end if

      a(i) = 0

    end do

    more = .false.

  end if

  return
end
subroutine get_seed ( seed )

!*****************************************************************************80
!
!! GET_SEED returns a seed for the random number generator.
!
!  Discussion:
!
!    The seed depends on the current time, and ought to be (slightly)
!    different every millisecond.  Once the seed is obtained, a random
!    number generator should be called a few times to further process
!    the seed.
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
!  Parameters:
!
!    Output, integer ( kind = 4 ) SEED, a pseudorandom seed value.
!
  implicit none

  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) seed
  real ( kind = 8 ) temp
  character ( len = 10 ) time
  character ( len = 8 )  today
  integer ( kind = 4 ) values(8)
  character ( len = 5 )  zone

  call date_and_time ( today, time, zone, values )

  temp = 0.0D+00

  temp = temp + real ( values(2) - 1, kind = 8 ) / 11.0D+00
  temp = temp + real ( values(3) - 1, kind = 8 ) / 30.0D+00
  temp = temp + real ( values(5),     kind = 8 ) / 23.0D+00
  temp = temp + real ( values(6),     kind = 8 ) / 59.0D+00
  temp = temp + real ( values(7),     kind = 8 ) / 59.0D+00
  temp = temp + real ( values(8),     kind = 8 ) / 999.0D+00
  temp = temp / 6.0D+00
!
!  Force 0 < TEMP <= 1.
!
  do while ( temp <= 0.0D+00 )
    temp = temp + 1.0D+00
  end do

  do while ( 1.0D+00 < temp )
    temp = temp - 1.0D+00
  end do

  seed = int ( real ( i4_huge, kind = 8 ) * temp )
!
!  Never use a seed of 0 or maximum integer.
!
  if ( seed == 0 ) then
    seed = 1
  end if

  if ( seed == i4_huge ) then
    seed = seed - 1
  end if

  return
end
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
function i4vec_eq ( n, a1, a2 )

!*****************************************************************************80
!
!! I4VEC_EQ is true if two I4VECs are equal.
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
!    13 May 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vectors.
!
!    Input, integer ( kind = 4 ) A1(N), A2(N), two vectors to compare.
!
!    Output, logical I4VEC_EQ, is TRUE if every pair of elements A1(I)
!    and A2(I) are equal, and FALSE otherwise.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a1(n)
  integer ( kind = 4 ) a2(n)
  logical i4vec_eq

  i4vec_eq = ( all ( a1(1:n) == a2(1:n) ) )

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
