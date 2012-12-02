subroutine prime_number ( n, total )

!*****************************************************************************80
!
!! PRIME_NUMBER returns the number of primes between 1 and N.
!
!  Discussion:
!
!    A naive algorithm is used.
!
!    Mathematica can return the number of primes less than or equal to N
!    by the command PrimePi[N].
!
!                N  PRIME_NUMBER
!
!                1           0
!               10           4
!              100          25
!            1,000         168
!           10,000       1,229
!          100,000       9,592
!        1,000,000      78,498
!       10,000,000     664,579
!      100,000,000   5,761,455
!    1,000,000,000  50,847,534
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the maximum number to check.
!
!    Output, integer ( kind = 4 ) TOTAL, the number of prime numbers up to N.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  integer ( kind = 4 ) prime
  integer ( kind = 4 ) total

  total = 0

  do i = 2, n

    prime = 1

    do j = 2, i - 1
      if ( mod ( i, j ) == 0 ) then
        prime = 0
        exit
      end if
    end do

    total = total + prime

  end do

  return
end
subroutine prime_number_sweep ( n_lo, n_hi, n_factor )

!*****************************************************************************80
!
!! PRIME_NUMBER_SWEEP does repeated timed calls to PRIME_NUMBER.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 August 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N_LO, the first value of N.
!
!    Input, integer ( kind = 4 ) N_HI, the last value of N.
!
!    Input, integer ( kind = 4 ) N_FACTOR, the factor by which to increase N
!    after each iteration.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_factor
  integer ( kind = 4 ) n_hi
  integer ( kind = 4 ) n_lo
  integer ( kind = 4 ) primes
  real ( kind = 8 ) time1
  real ( kind = 8 ) time2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PRIME_NUMBER_SWEEP'
  write ( *, '(a)' ) '  Call PRIME_NUMBER to count the primes from 1 to N.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N        Pi          Time'
  write ( *, '(a)' ) ' '

  n = n_lo

  do while ( n <= n_hi )

    call cpu_time ( time1 )

    call prime_number ( n, primes )

    call cpu_time ( time2 )

    write ( *, '(2x,i8,2x,i8,g14.6)' ) n, primes, time2 - time1

    n = n * n_factor

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
