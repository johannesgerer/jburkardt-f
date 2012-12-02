subroutine apportion_adams ( state_num, state_pop, rep_num, state_rep )

!*****************************************************************************80
!
!! APPORTION_ADAMS apportions using Adams's method.
!
!  Discussion:
!
!    Because rounding up is involved, no state will ever get 0 representatives,
!    as long as there are at least as many representatives as states!
!
!    However, this procedure will fail if a tie situation occurs, in which,
!    say, the number of representatives assigned jumps by 2 at some point,
!    skipping the desired value, no matter how small we subdivide the interval.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 May 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) STATE_NUM, the number of states entitled 
!    to representation.
!
!    Input, integer ( kind = 4 ) STATE_POP(STATE_NUM), the population of 
!    each state.
!
!    Input, integer ( kind = 4 ) REP_NUM, the number of representatives 
!    in the house.
!
!    Output, integer ( kind = 4 ) STATE_REP(STATE_NUM), the number of 
!    representatives assigned to each state.
!
  implicit none

  integer ( kind = 4 ) state_num

  integer ( kind = 4 ) rep_max
  integer ( kind = 4 ) rep_mid
  integer ( kind = 4 ) rep_min
  integer ( kind = 4 ) rep_num
  real ( kind = 8 ) sd_max
  real ( kind = 8 ) sd_mid
  real ( kind = 8 ) sd_min
  real ( kind = 8 ) sq(state_num)
  integer ( kind = 4 ) state_pop(state_num)
  integer ( kind = 4 ) state_rep(state_num)
  integer ( kind = 4 ) uq(state_num)
  integer ( kind = 4 ) us_pop
!
!  Compute the total population.
!
  us_pop = sum ( state_pop(1:state_num) )
!
!  Compute the standard divisor.
!  If we round up the quotas for each state, we'll probably get 
!  too many representatives.
!
  sd_min = real ( us_pop, kind = 8 ) / real ( rep_num, kind = 8 )
  sq(1:state_num) = real ( state_pop(1:state_num), kind = 8 ) / sd_min
  call r8vec_ceiling ( state_num, sq, uq )
  rep_min = sum ( uq(1:state_num) )
!
!  Try a larger standard divisor.
!  If we round up the quotas for each state, we should get too 
!  few representatives.
!
  sd_max = 2.0D+00 * sd_min
  sq(1:state_num) = real ( state_pop(1:state_num), kind = 8 ) / sd_max
  call r8vec_ceiling ( state_num, sq, uq )
  rep_max = sum ( uq(1:state_num) )
!
!  Now try bisection, and see if you can locate a suitable divisor.
!
  do

    sd_mid = 0.5D+00 * ( sd_min + sd_max )
    sq(1:state_num) = real ( state_pop(1:state_num), kind = 8 ) / sd_mid
    call r8vec_ceiling ( state_num, sq, uq )
    rep_mid = sum ( uq(1:state_num) )

    if ( rep_mid == rep_num ) then
      state_rep(1:state_num) = uq(1:state_num)
      exit
    else if ( rep_mid < rep_num ) then
      sd_max = sd_mid
      rep_max = rep_mid
    else
      sd_min = sd_mid
      rep_min = rep_mid
    end if

  end do

  return
end
subroutine apportion_hamilton ( state_num, state_pop, rep_num, state_rep )

!*****************************************************************************80
!
!! APPORTION_HAMILTON apportions using Hamilton's method.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 April 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) STATE_NUM, the number of states entitled 
!    to representation.
!
!    Input, integer ( kind = 4 ) STATE_POP(STATE_NUM), the population of 
!    each state.
!
!    Input, integer ( kind = 4 ) REP_NUM, the number of representatives 
!    in the house.
!
!    Output, integer ( kind = 4 ) STATE_REP(STATE_NUM), the number of 
!    representatives assigned to each state.
!
  implicit none

  integer ( kind = 4 ) state_num

  real ( kind = 8 ) fraction(state_num)
  integer ( kind = 4 ) indx(state_num)
  integer ( kind = 4 ) rep_num
  integer ( kind = 4 ) rep_remain
  integer ( kind = 4 ) rep_set
  integer ( kind = 4 ) rq(state_num)
  real ( kind = 8 ) sd
  real ( kind = 8 ) sq(state_num)
  integer ( kind = 4 ) state_pop(state_num)
  integer ( kind = 4 ) state_rep(state_num)
  integer ( kind = 4 ) us_pop
!
!  Compute the total population.
!
  us_pop = sum ( state_pop(1:state_num) )
!
!  Compute the standard divisor.
!
  sd = real ( us_pop, kind = 8 ) / real ( rep_num, kind = 8 )
!
!  Compute the standard quota for each state.
!
  sq(1:state_num) = real ( state_pop(1:state_num), kind = 8 ) / sd
!
!  The rounded quota is the lower quota for each state.
!
  call r8vec_floor ( state_num, sq, rq )
!
!  Most representatives have been assigned.
!
  rep_set = sum ( rq(1:state_num) )
  rep_remain = rep_num - rep_set
!
!  If there are still representatives to assign, assign them
!  to states by order of the fractional part of the standard quota.
!
  if ( 0 < rep_remain ) then

    call r8vec_fraction ( state_num, sq, fraction )
 
    call r8vec_sort_heap_index_d ( state_num, fraction, indx )

    rq( indx(1:rep_remain) ) = rq ( indx(1:rep_remain) ) + 1

  end if

  state_rep(1:state_num) = rq(1:state_num)

  return
end
subroutine apportion_huntington_hill ( state_num, state_pop, rep_num, &
  state_rep )

!*****************************************************************************80
!
!! APPORTION_HUNTINGTON_HILL uses the Huntington-Hill apportionment.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 June 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) STATE_NUM, the number of states entitled 
!    to representation.
!
!    Input, integer ( kind = 4 ) STATE_POP(STATE_NUM), the population 
!    of each state.
!
!    Input, integer ( kind = 4 ) REP_NUM, the number of representative 
!    in the house.
!
!    Output, integer ( kind = 4 ) STATE_REP(STATE_NUM), the number of 
!    representative assigned to each state.
!
  implicit none

  integer ( kind = 4 ) state_num

  integer ( kind = 4 ) indx(state_num)
  real ( kind = 8 ) r8_sqrt_i4
  integer ( kind = 4 ) rep
  integer ( kind = 4 ) rep_num
  integer ( kind = 4 ) rep_remain
  integer ( kind = 4 ) state
  real ( kind = 8 ) state_d(state_num)
  real ( kind = 8 ) state_p(state_num)
  integer ( kind = 4 ) state_pop(state_num)
  integer ( kind = 4 ) state_rep(state_num)
  integer ( kind = 4 ) us_pop
!
!  Start by assigning one representative to each state.
!
  state_rep(1:state_num) = 1
!
!  Set each state's D value to sqrt ( state_rep * ( state_rep + 1 ) ).
!
  do state = 1, state_num
    state_d(state) = r8_sqrt_i4 ( state_rep(state) * ( state_rep(state) + 1 ) )
  end do
!
!  Set state's P (priority) value to state_pop / state_d
!
  do state = 1, state_num
    state_p(state) = real ( state_pop(state), kind = 8 ) / state_d(state)
  end do
!
!  Index sort the states by STATE_P.
!
  call r8vec_sort_heap_index_d ( state_num, state_p, indx )
!
!  Assign remaining representatives.
!
  rep_remain = rep_num - state_num

  do rep = 1, rep_remain
!
!  State with index(1) gets next rep.
!
    state = indx(1)
    state_rep(state) = state_rep(state) + 1
    state_d(state) = r8_sqrt_i4 ( state_rep(state) * ( state_rep(state) + 1 ) )
    state_p(state) = real ( state_pop(state), kind = 8 ) / state_d(state)
!
!  Resort array by determining new placement of STATE_P(STATE).
!  For now, be wasteful and just resort the whole array from scratch.
!
    call r8vec_sort_heap_index_d ( state_num, state_p, indx )

  end do
    
  return
end
subroutine apportion_jefferson ( state_num, state_pop, rep_num, state_rep )

!*****************************************************************************80
!
!! APPORTION_JEFFERSON apportions using Jefferson's method.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 May 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) STATE_NUM, the number of states entitled 
!    to representation.
!
!    Input, integer ( kind = 4 ) STATE_POP(STATE_NUM), the population of 
!    each state.
!
!    Input, integer ( kind = 4 ) REP_NUM, the number of representatives 
!    in the house.
!
!    Output, integer ( kind = 4 ) STATE_REP(STATE_NUM), the number of 
!    representatives assigned to each state.
!
  implicit none

  integer ( kind = 4 ) state_num

  integer ( kind = 4 ) rep
  integer ( kind = 4 ) rep_num
  integer ( kind = 4 ) state(1)
  integer ( kind = 4 ) state_pop(state_num)
  integer ( kind = 4 ) state_rep(state_num)
  real ( kind = 8 ) test(state_num)

  if ( rep_num < state_num ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'APPORTIONMENT_JEFFERSON - Fatal error!'
    write ( *, '(a)' ) '  More states than representatives.'
    stop
  end if
!
!  Give every state one representative.
!
  state_rep(1:state_num) = 1
!
!  Give the next representative to whichever state has the highest
!  value of pop / (state_rep+1).
!
  do rep = state_num + 1, rep_num

    test(1:state_num) = real ( state_pop(1:state_num), kind = 8 ) &
      / real ( state_rep(1:state_num) + 1, kind = 8 )

    state = maxloc ( test(1:state_num) )

    state_rep(state(1)) = state_rep(state(1)) + 1

  end do

  return
end
subroutine apportion_webster ( state_num, state_pop, rep_num, state_rep )

!*****************************************************************************80
!
!! APPORTION_WEBSTER apportions using Webster's method.
!
!  Discussion:
!
!    Small states could get 0 representatives under this scheme.
!
!    Also, this procedure will fail if a tie situation occurs, in which,
!    say, the number of representatives assigned jumps by 2 at some point,
!    skipping the desired value, no matter how small we subdivide the interval.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 May 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) STATE_NUM, the number of states entitled 
!    to representation.
!
!    Input, integer ( kind = 4 ) STATE_POP(STATE_NUM), the population of 
!    each state.
!
!    Input, integer ( kind = 4 ) REP_NUM, the number of representatives 
!    in the house.
!
!    Output, integer ( kind = 4 ) STATE_REP(STATE_NUM), the number of 
!    representatives assigned to each state.
!
  implicit none

  integer ( kind = 4 ) state_num

  integer ( kind = 4 ) rep
  integer ( kind = 4 ) rep_max
  integer ( kind = 4 ) rep_mid
  integer ( kind = 4 ) rep_min
  integer ( kind = 4 ) rep_num
  integer ( kind = 4 ) rq(state_num)
  real ( kind = 8 ) sd
  real ( kind = 8 ) sd_max
  real ( kind = 8 ) sd_mid
  real ( kind = 8 ) sd_min
  real ( kind = 8 ) sq(state_num)
  integer ( kind = 4 ) state_pop(state_num)
  integer ( kind = 4 ) state_rep(state_num)
  integer ( kind = 4 ) tries
  integer ( kind = 4 ) us_pop
!
!  Compute the total population.
!
  us_pop = sum ( state_pop(1:state_num) )
!
!  Compute the standard divisor.
!
  sd = real ( us_pop, kind = 8 ) / real ( rep_num, kind = 8 )
  sq(1:state_num) = real ( state_pop(1:state_num), kind = 8 ) / sd
  call r8vec_round ( state_num, sq, rq )
  rep = sum ( rq(1:state_num) )

  if ( rep == rep_num ) then
    state_rep(1:state_num) = rq(1:state_num)
  else if ( rep < rep_num ) then
    sd_min = sd
    rep_min = rep
    sd_max = sd / 2.0D+00
    sq(1:state_num) = real ( state_pop(1:state_num), kind = 8 ) / sd_max
    call r8vec_round ( state_num, sq, rq )
    rep_max = sum ( rq(1:state_num) )
  else if ( rep_num < rep ) then
    sd_max = sd
    rep_max = rep
    sd_min = sd * 2.0D+00
    sq(1:state_num) = real ( state_pop(1:state_num), kind = 8 ) / sd_min
    call r8vec_round ( state_num, sq, rq )
    rep_min = sum ( rq(1:state_num) )
  end if
!
!  Now try bisection, and see if you can locate a suitable divisor.
!
  tries = 0

  do

    sd_mid = 0.5D+00 * ( sd_min + sd_max )
    sq(1:state_num) = real ( state_pop(1:state_num), kind = 8 ) / sd_mid
    call r8vec_round ( state_num, sq, rq )
    rep_mid = sum ( rq(1:state_num) )

    if ( rep_mid == rep_num ) then
      state_rep(1:state_num) = rq(1:state_num)
      exit
    else if ( rep_mid < rep_num ) then
      sd_min = sd_mid
      rep_min = rep_mid
    else
      sd_max = sd_mid
      rep_max = rep_mid
    end if

    tries = tries + 1

    if ( 50 < tries ) then
      write ( *, * ) ' '
      write ( *, * ) 'APPORTIONMENT_WEBSTER - Fatal error!'
      write ( *, * ) '  Too many tries!'
      stop
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
subroutine i4_to_month_name ( m, month_name )

!*****************************************************************************80
!
!! I4_TO_MONTH_NAME returns the name of a given month.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 April 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of the month, which should
!    be between 1 and 12.
!
!    Output, character ( len = * ) MONTH_NAME, a string containing as much 
!    of the month's name as will fit.  To get the typical 3-letter 
!    abbreviations for the months, simply declare
!      character ( len = 3 ) MONTH_NAME
!    or pass in MONTH_NAME(1:3).
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  character ( len = * ) month_name
  character ( len = 9 ), parameter, dimension(12) :: name = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)

  if ( m < 1 .or. 12 < m ) then

    do i = 1, len ( month_name )
      month_name(i:i) = '?'
    end do

  else

    month_name = name(m)

  end if
 
  return
end
subroutine i4_to_s_right_comma ( i4, s )

!*****************************************************************************80
!
!! I4_TO_S_RIGHT_COMMA converts an I4 to a right justified string with commas.
!
!  Example:
!
!    Assume that S is 10 characters long:
!
!          I4           S
!
!           1           1
!          -1          -1
!           0           0
!        1952       1,952
!      123456     123,456
!     1234567   1,234,567
!    12345678  12,345,678
!   123456789  **********  <-- Not enough room!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 April 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I4, an integer to be converted.
!
!    Output, character ( len = * ) S, the representation of the integer.
!    The integer will be right-justified.  Commas will be used to separate
!    sets of three digits.  If there is not enough space, the string will 
!    be filled with stars.
!
  implicit none

  character c
  integer ( kind = 4 ) digit
  integer ( kind = 4 ) digit_num
  integer ( kind = 4 ) hi
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4
  integer ( kind = 4 ) lo
  integer ( kind = 4 ) pos
  character ( len = * ) s
  integer ( kind = 4 ) value

  s = ' '

  lo = 1
  hi = len ( s )

  if ( hi <= 0 ) then
    return
  end if
!
!  Make a copy of the integer.
!
  value = i4
!
!  Handle the negative sign.
!
  if ( value < 0 ) then

    if ( hi <= 1 ) then
      s(1:1) = '*'
      return
    end if

    value = -value
    s(1:1) = '-'
    lo = 2

  end if
!
!  The absolute value of the integer goes into S(LO:HI).
!
  pos = hi
!
!  Find the last digit of VALUE, strip it off, and stick it into the string.
!
  digit_num = 0

  do

    digit = mod ( value, 10 )
    value = value / 10
    digit_num = digit_num + 1

    if ( pos < lo ) then
      do i = 1, hi
        s(i:i) = '*'
      end do
      return
    end if
!
!  Insert a comma?
!
    if ( 1 < digit_num .and. mod ( digit_num, 3 ) == 1 ) then

      if ( pos < lo ) then
        do i = 1, hi
          s(i:i) = '*'
        end do
        return
      end if

      s(pos:pos) = ','
      pos = pos - 1
    end if

    call digit_to_ch ( digit, c )
    s(pos:pos) = c
    pos = pos - 1

    if ( value == 0 ) then
      exit
    end if

  end do
!
!  Shift the minus sign, if any.
!
  if ( s(1:1) == '-' ) then
    if ( pos /= 1 ) then
      s(1:1) = ' '
      s(pos:pos) = '-'
    end if
  end if

  return
end
subroutine i4vec_nonzero_first ( n, x, nz, indx )

!*****************************************************************************80
!
!! I4VEC_NONZERO_FIRST left-shifts all nonzeros in an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of integer values.
!
!    The routine preserves the ordering of the nonzero entries.  It counts
!    the nonzeros, and returns an index vector.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 April 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input/output, integer ( kind = 4 ) X(N), the vector to be shifted.
!
!    Output, integer ( kind = 4 ) NZ, the number of nonzero entries.
!
!    Output, integer ( kind = 4 ) INDX(N), contains the original location 
!    of each entry.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) nz
  integer ( kind = 4 ) x(n)

  nz = 0

  do j = 1, n
    indx(j) = j
  end do

  j = 0

  do while ( j < n )

    j = j + 1

    if ( x(j) /= 0 ) then

      nz = nz + 1

      if ( nz /= j ) then

        x(nz) = x(j)
        x(j) = 0

        k = indx(nz)
        indx(nz) = j
        indx(j) = k

      end if
    end if
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
function r8_div_i4 ( i, j )

!*****************************************************************************80
!
!! R8_DIV_I4 returns an I4 fraction as an R8.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 June 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, J, the numerator and denominator.
!
!    Output, integer ( kind = 4 ) R8_DIV_I4, the value of (I/J).
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) r8_div_i4

  r8_div_i4 = real ( i, kind = 8 ) / real ( j, kind = 8 )

  return
end
function r8_sqrt_i4 ( i )

!*****************************************************************************80
!
!! R8_SQRT_I4 returns the square root of an I4 as an R8.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 June 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the number whose square root is desired.
!
!    Output, integer ( kind = 4 ) R8_SQRT_I4, the value of sqrt(I).
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) r8_sqrt_i4

  r8_sqrt_i4 = sqrt ( real ( i, kind = 8 ) )

  return
end
subroutine r8vec_ceiling ( n, r8vec, ceilingvec )

!*****************************************************************************80
!
!! R8VEC_CEILING rounds "up" (towards +infinity) entries of an R8VEC.
!
!  Example:
!
!    R8    Value
!
!   -1.1  -1
!   -1.0  -1
!   -0.9   0
!    0.0   0
!    5.0   5
!    5.1   6
!    5.9   6
!    6.0   6
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 April 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries.
!
!    Input, real ( kind = 8 ) R8VEC(N), the values to be rounded up.
!
!    Output, integer ( kind = 4 ) CEILINGVEC(N), the rounded values.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) ceilingvec(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) r8vec(n)
  integer ( kind = 4 ) value

  do i = 1, n

    value = int ( r8vec(i) )

    if ( real ( value, kind = 8 ) < r8vec(i) ) then
      value = value + 1
    end if

    ceilingvec(i) = value

  end do

  return
end
subroutine r8vec_floor ( n, r8vec, floorvec )

!*****************************************************************************80
!
!! R8VEC_FLOOR rounds "down" (towards -infinity) entries of an R8VEC.
!
!  Example:
!
!    R8    Value
!
!   -1.1  -2
!   -1.0  -1
!   -0.9  -1
!    0.0   0
!    5.0   5
!    5.1   5
!    5.9   5
!    6.0   6
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 April 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries.
!
!    Input, real ( kind = 8 ) R8VEC(N), the values to be rounded down.
!
!    Output, integer ( kind = 4 ) FLOORVEC(N), the rounded value.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) floorvec(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) r8vec(n)
  integer ( kind = 4 ) value

  do i = 1, n

    value = int ( r8vec(i) )

    if ( r8vec(i) < real ( value, kind = 8 ) ) then
      value = value - 1
    end if

    floorvec(i) = value

  end do

  return
end
subroutine r8vec_fraction ( n, x, fraction )

!*****************************************************************************80
!
!! R8VEC_FRACTION returns the fraction parts of an R8VEC.
!
!  Discussion:
!
!    If we regard a real number as
!
!      R8 = SIGN * ( WHOLE + FRACTION )
!
!    where
!
!      SIGN is +1 or -1,
!      WHOLE is a nonnegative integer
!      FRACTION is a nonnegative real number strictly less than 1,
!
!    then this routine returns the value of FRACTION.
!
!  Example:
!
!     R8    R8_FRACTION
! 
!    0.00      0.00
!    1.01      0.01
!    2.02      0.02
!   19.73      0.73
!   -4.34      0.34
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 April 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of arguments.
!
!    Input, real ( kind = 8 ) X(N), the arguments.
!
!    Output, real ( kind = 8 ) FRACTION(N), the fraction parts.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fraction(n)
  real ( kind = 8 ) x(n)

  fraction(1:n) = abs ( x(1:n) ) - real ( int ( abs ( x(1:n) ) ), kind = 8 )

  return
end
subroutine r8vec_round ( n, a, b )

!*****************************************************************************80
!
!! R8VEC_ROUND rounds entries of an R8VEC.
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
!    10 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, real ( kind = 8 ) A(N), the vector to be rounded.
!
!    Output, integer ( kind = 4 ) B(N), the rounded vector.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) b(n)

  b(1:n) = nint ( a(1:n) )

  return
end
subroutine r8vec_sort_heap_index_a ( n, a, indx )

!*****************************************************************************80
!
!! R8VEC_SORT_HEAP_INDEX_A does an indexed heap ascending sort of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
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
!      call r8vec_permute ( n, a, indx )
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

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) aval
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
subroutine r8vec_sort_heap_index_d ( n, a, indx )

!*****************************************************************************80
!
!! R8VEC_SORT_HEAP_INDEX_D does an indexed heap descending sort of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    The sorting is not actually carried out.  Rather an index array is
!    created which defines the sorting.  This array may be used to sort
!    or index the array, or to sort or index related arrays keyed on the
!    original array.
!
!    Once the index array is computed, the sorting can be carried out
!    "implicitly:
!
!      A(INDX(1:N)) is sorted,
!
!    or explicitly, by the call
!
!      call r8vec_permute ( n, a, index )
!
!    after which A(1:N) is sorted.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 July 2000
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

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) aval
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
        if ( a(indx(j+1)) < a(indx(j)) ) then
          j = j + 1
        end if
      end if

      if ( a(indx(j)) < aval ) then
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
subroutine rep_num_year ( year, rep_num )

!*****************************************************************************80
!
!! REP_NUM_YEAR returns the size of the House of Representatives.
!
!  Discussion:
!
!    Every 10 years, a new census is carried out, and the size of the
!    House of Representatives may be modified.
!
!    We report here only the size authorized at the time of the census.
!    After a particular size was authorized, the number of representatives
!    could be modified, particularly if a new state was admitted.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 April 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) YEAR, the year.
!
!    Output, integer ( kind = 4 ) REP_NUM, the authorized size of the House of
!    Representatives, based on the most recent census preceding or
!    equal to year YEAR.
!
  implicit none

  integer ( kind = 4 ) rep_num
  integer ( kind = 4 ) year

  if ( year < 1789 ) then
    rep_num = 0
  else if ( year < 1790 ) then
    rep_num = 65
  else if ( year < 1800 ) then
    rep_num = 105
  else if ( year < 1810 ) then
    rep_num = 141
  else if ( year < 1820 ) then
    rep_num = 181
  else if ( year < 1830 ) then
    rep_num = 213
  else if ( year < 1840 ) then
    rep_num = 240
  else if ( year < 1850 ) then
    rep_num = 223
  else if ( year < 1860 ) then
    rep_num = 234
  else if ( year < 1870 ) then
    rep_num = 241
  else if ( year < 1880 ) then
    rep_num = 292
  else if ( year < 1890 ) then
    rep_num = 325
  else if ( year < 1900 ) then
    rep_num = 356
  else if ( year < 1910 ) then
    rep_num = 386
  else
    rep_num = 435
  end if

  return
end
function state_id ( state )

!*****************************************************************************80
!
!! STATE_ID returns the 2 letter Postal Code for one of the 50 states and DC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 April 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) STATE, the index of a state.
!
!    Output, character ( len = 2 ) STATE_ID, the 2 letter code.
!
  implicit none

  character ( len = 2 ), parameter, dimension ( 51 ) :: id = (/ &
    'AL', 'AK', 'AZ', 'AR', 'CA', &
    'CO', 'CT', 'DE', 'DC', 'FL', &
    'GA', 'HI', 'ID', 'IL', 'IN', &
    'IA', 'KS', 'KY', 'LA', 'ME', &
    'MD', 'MA', 'MI', 'MN', 'MS', &
    'MO', 'MT', 'NE', 'NV', 'NH', &
    'NJ', 'NM', 'NY', 'NC', 'ND', &
    'OH', 'OK', 'OR', 'PA', 'RI', &
    'SC', 'SD', 'TN', 'TX', 'UT', &
    'VT', 'VA', 'WA', 'WV', 'WI', &
    'WY' /)
  integer ( kind = 4 ) state
  character ( len = 2 ) state_id

  if ( state < 1 ) then
    state_id = '??'
  else if ( state <= 51 ) then
    state_id = id(state)
  else
    state_id = '??'
  end if

  return
end
function state_name ( state )

!*****************************************************************************80
!
!! STATE_NAME returns the name of one of the 50 states plus DC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 April 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) STATE, the index of a state.
!
!    Output, character ( len = 20 ) STATE_NAME, the name of the state.
!
  implicit none

  character ( len = 20 ), parameter, dimension ( 51 ) :: name = (/ &
  'Alabama             ', &
  'Alaska              ', &
  'Arizona             ', &
  'Arkansas            ', &
  'California          ', &
  'Colorado            ', &
  'Connecticut         ', &
  'Delaware            ', &
  'District of Columbia', &
  'Florida             ', &
  'Georgia             ', &
  'Hawaii              ', &
  'Idaho               ', &
  'Illinois            ', &
  'Indiana             ', &
  'Iowa                ', &
  'Kansas              ', &
  'Kentucky            ', &
  'Louisiana           ', &
  'Maine               ', &
  'Maryland            ', &
  'Massachusetts       ', &
  'Michigan            ', &
  'Minnesota           ', &
  'Mississippi         ', &
  'Missouri            ', &
  'Montana             ', &
  'Nebraska            ', &
  'Nevada              ', &
  'New Hampshire       ', &
  'New Jersey          ', &
  'New Mexico          ', &
  'New York            ', &
  'North Carolina      ', &
  'North Dakota        ', &
  'Ohio                ', &
  'Oklahoma            ', &
  'Oregon              ', &
  'Pennsylvania        ', &
  'Rhode Island        ', &
  'South Carolina      ', &
  'South Dakota        ', &
  'Tennessee           ', &
  'Texas               ', &
  'Utah                ', &
  'Vermont             ', &
  'Virginia            ', &
  'Washington          ', &
  'West Virginia       ', &
  'Wisconsin           ', &
  'Wyoming             ' /)
  integer ( kind = 4 ) state
  character ( len = 20 ) state_name

  if ( state < 1 ) then
    state_name = '????????????????????'
  else if ( state <= 51 ) then
    state_name = name(state)
  else
    state_name = '????????????????????'
  end if

  return
end
subroutine state_num_year ( year, state_num )

!*****************************************************************************80
!
!! STATE_NUM_YEAR returns the number of states for a given year.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 April 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) YEAR, the year.
!
!    Output, integer ( kind = 4 ) STATE_NUM, the number of states in the union,
!    at the end of year YEAR.
!
  implicit none

  integer ( kind = 4 ) state_num
  integer ( kind = 4 ) year

  if ( year < 1787 ) then
    state_num = 0
  else if ( year < 1788  ) then
    state_num = 3
  else if ( year < 1789 ) then
    state_num = 11
  else if ( year < 1790 ) then
    state_num = 12
  else if ( year < 1791 ) then
    state_num = 13
  else if ( year < 1792 ) then
    state_num = 14
  else if ( year < 1796 ) then
    state_num = 15
  else if ( year < 1803 ) then
    state_num = 16
  else if ( year < 1812 ) then
    state_num = 17
  else if ( year < 1816 ) then
    state_num = 18
  else if ( year < 1817 ) then
    state_num = 19
  else if ( year < 1818 ) then
    state_num = 20
  else if ( year < 1819 ) then
    state_num = 21
  else if ( year < 1820 ) then
    state_num = 22
  else if ( year < 1821 ) then
    state_num = 23
  else if ( year < 1836 ) then
    state_num = 24
  else if ( year < 1837 ) then
    state_num = 25
  else if ( year < 1845 ) then
    state_num = 26
  else if ( year < 1846 ) then
    state_num = 28
  else if ( year < 1848 ) then
    state_num = 29
  else if ( year < 1850 ) then
    state_num = 30
  else if ( year < 1858 ) then
    state_num = 31
  else if ( year < 1859 ) then
    state_num = 32
  else if ( year < 1861 ) then
    state_num = 33
  else if ( year < 1863 ) then
    state_num = 34
  else if ( year < 1864 ) then
    state_num = 35
  else if ( year < 1867 ) then
    state_num = 36
  else if ( year < 1876 ) then
    state_num = 37
  else if ( year < 1889 ) then
    state_num = 38
  else if ( year < 1890 ) then
    state_num = 42
  else if ( year < 1896 ) then
    state_num = 44
  else if ( year < 1907 ) then
    state_num = 45
  else if ( year < 1912 ) then
    state_num = 46
  else if ( year < 1959 ) then
    state_num = 48
  else
    state_num = 50
  end if

  return
end
subroutine state_pop_year ( year, state_pop )

!*****************************************************************************80
!
!! STATE_POP_YEAR returns the state populations for a given census year.
!
!  Discussion:
!
!    No matter what the year, values are returned for all 50 states
!    and the District of Columbia.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 April 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) YEAR, the census year.
!
!    Output, integer ( kind = 4 ) STATE_POP(51), the state populations.
!
  implicit none

  integer ( kind = 4 ) state
  integer ( kind = 4 ) state_pop(51)
  integer ( kind = 4 ) year

  if ( year < 1790 ) then
    state_pop(1:51) = 0
  else if ( year < 1970 ) then
    state_pop(1:51) = (/ &
    3266740,     226167,    1302161,    1786272,   15717204, &
    1753947,    2535234,     446292,     763956,    4951560, &
    3943116,     632772,     667191,   10081158,    4662498, &
    2757537,    2178611,    3038156,    3257022,     969265, &
    3100689,    5148578,    7823194,    3413864,    2178141, &
    4319813,     674767,    1411330,     285278,     606921, &
    6066782,     951023,   16782304,    4556155,     632446, &
    9706397,    2328284,    1768687,   11319366,     859488, &
    2382594,     680514,    3567089,    9579677,     890627, &
     389881,    3966949,    2853214,    1860421,    3951777, &
     330066 /)
  else
    state_pop(1:51) = 0
  end if

  return
end
subroutine state_rep_year ( year, state_rep )

!*****************************************************************************80
!
!! STATE_REP_YEAR returns the state representatives for a given census year.
!
!  Discussion:
!
!    No matter what the year, values are returned for all 50 states
!    and the District of Columbia.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 April 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) YEAR, the census year.
!
!    Output, integer ( kind = 4 ) STATE_REP(51), the state representatives.
!
  implicit none

  integer ( kind = 4 ) state
  integer ( kind = 4 ) state_rep(51)
  integer ( kind = 4 ) year

  if ( year < 1790 ) then
    state_rep(1:51) = 0
  else if ( year < 1970 ) then
    state_rep(1:51) = (/ &
       8,  1,  3,  4, 38,  4,  6,  1,  0, 12, &
      10,  2,  2, 24, 11,  7,  5,  7,  8,  2, &
       8, 12, 19,  8,  5, 10,  2,  3,  1,  2, &
      15,  2, 41, 11,  2, 24,  6,  4, 27,  2, &
       6,  2,  9, 23,  2,  1, 10,  7,  5, 10, &
       1 /)
  else
    state_rep(1:51) = 0
  end if

  return
end
subroutine state_statehood  ( state, y, m, d )

!*****************************************************************************80
!
!! STATE_STATEHOOD returns the statehood year for each state.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 April 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) STATE, the index of a state.
!
!    Output, integer ( kind = 4 ) Y, M, D, the year, month and day of statehood.
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ), dimension ( 51 ) :: d_state = (/ &
    14, 03, 14, 15, 09, &
    01, 09, 07, 31, 03, &
    02, 21, 03, 03, 11, &
    28, 29, 01, 30, 15, &
    28, 06, 26, 11, 10, &
    10, 08, 01, 31, 21, &
    18, 06, 26, 21, 02, &
    01, 16, 14, 12, 29, &
    23, 02, 01, 29, 04, &
    04, 25, 11, 20, 29, &
    10 /)
  integer ( kind = 4 ) m
  integer ( kind = 4 ), dimension ( 51 ) :: m_state = (/ &
    12, 01, 02, 06, 09, &
    08, 01, 12, 12, 03, &
    01, 08, 07, 12, 12, &
    12, 01, 06, 04, 03, &
    04, 02, 01, 05, 12, &
    08, 11, 03, 10, 06, &
    12, 01, 07, 11, 11, &
    03, 11, 02, 12, 05, &
    05, 11, 06, 12, 01, &
    03, 06, 11, 06, 05, &
    07 /)
  integer ( kind = 4 ) state
  integer ( kind = 4 ) y
  integer ( kind = 4 ), dimension ( 51 ) :: y_state = (/ &
    1819,  1959,  1912,  1836,  1850, &
    1876,  1788,  1787,  3000,  1845, &
    1788,  1959,  1890,  1818,  1816, &
    1846,  1861,  1792,  1812,  1820, &
    1788,  1788,  1837,  1858,  1817, &
    1821,  1889,  1867,  1864,  1788, &
    1787,  1912,  1788,  1789,  1889, &
    1803,  1907,  1859,  1787,  1790, &
    1788,  1889,  1796,  1845,  1896, &
    1791,  1788,  1889,  1863,  1848, &
    1890 /)

  if ( state < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'STATE_STATEHOOD - Fatal error!'
    write ( *, '(a)' ) '  Input STATE < 1.'
    stop
  else if ( state <= 51 ) then
    y = y_state(state)
    m = m_state(state)
    d = d_state(state)
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'STATE_STATEHOOD - Fatal error!'
    write ( *, '(a)' ) '  51 < STATE.'
    stop
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
