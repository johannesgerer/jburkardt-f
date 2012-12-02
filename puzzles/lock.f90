program main

!*****************************************************************************80
!
!! MAIN is the main program for LOCK.
!
!  Discussion:
!
!    A random 4 digit TARGET is given.
!
!    The user starts at another random 4 digit starting point.
!
!    The user can indicate which digit to be advanced, and by how much.
!    That digit is advanced, modulo 10.
!
!    When the target is achieved, the user has won.
!
!    The problem is that the digits are "sticky".  So advancing digit 1
!    will cause digit 2 to move the same amount.  If digit 2 is advanced,
!    it will cause digits 1 and 3 to move the same.  Similarly, digit 4
!    will cause digit 3 to move as well.
!
!    Although the problem is still solvable, and solvable in at most 4 steps,
!    the complications of stickiness make it harder to see what to do.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 May 2011
!
!  Author:
!
!    John Burkardt
!
  integer ( kind = 4 ) dial
  integer ( kind = 4 ) digit(4)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) target(4)
  integer ( kind = 4 ) turns

  call get_seed ( seed )

  call i4vec_uniform ( 4, 0, 9, seed, target )
  call i4vec_uniform ( 4, 0, 9, seed, digit )

  do

    write ( *, '(a)' ) ' '
    write ( *, '(2x,a,i1,a1,i1,a1,i1,a1,i1)' ) &
      'Target:  ', target(1), '/', target(2), '/', target(3), '/', target(4)
    write ( *, '(2x,a,i1,a1,i1,a1,i1,a1,i1)' ) &
      'Current: ', digit(1), '/', digit(2), '/', digit(3), '/', digit(4)

    if ( all ( digit(1:4) == target(1:4) ) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  You have reached the combination!'
      exit
    end if

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter Dial and number of Turns: '

    read ( *, * ) dial, turns

    if ( dial < 1 .or. 4 < dial ) then
      write ( *, '(a)' ) '  Illegal value of dial!'
      cycle
    end if

    if ( 1 < dial ) then
      digit(dial-1) = digit(dial-1) + turns
    end if
    digit(dial) = digit(dial) + turns
    if ( dial < 4 ) then
      digit(dial+1) = digit(dial+1) + turns
    end if

    do i = 1, 4
      digit(i) = i4_modp ( digit(i), 10 )
    end do

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LOCK:'
  write ( *, '(a)' ) '  Normal end of execution.'

  stop
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

  integer   ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer   ( kind = 4 ) seed
  real      ( kind = 8 ) temp
  character ( len = 10 ) time
  character ( len = 8 )  today
  integer   ( kind = 4 ) values(8)
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
subroutine i4vec_uniform ( n, a, b, seed, x )

!*****************************************************************************80
!
!! I4VEC_UNIFORM returns a scaled pseudorandom I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!    The pseudorandom numbers should be scaled to be uniformly distributed
!    between A and B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 November 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the vector.
!
!    Input, integer ( kind = 4 ) A, B, the limits of the interval.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, integer ( kind = 4 ) X(N), a vector of numbers between A and B.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) k
  real    ( kind = 4 ) r
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) value
  integer ( kind = 4 ) x(n)

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4VEC_UNIFORM - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  do i = 1, n

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

    x(i) = value

  end do

  return
end
