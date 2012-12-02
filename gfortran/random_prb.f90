program main

!*******************************************************************************
!
!! MAIN is the main program for RANDOM_PRB.
!
!  Discussion:
!
!    RANDOM_PRB tests the random number generator.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'RANDOM_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version.'
  write ( *, '(a)' ) '  Tests for random number generation.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'RANDOM_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*******************************************************************************
!
!! TEST01 tests the FORTRAN90 random seed setting.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 February 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer i
  integer iseed
  integer itest
  integer, parameter :: n = 10
  real x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Try to set the FORTRAN90 random number seed.'

  do itest = 1, 4

    if ( itest <= 2 ) then
      iseed = 0
    else
      iseed = 1492
    end if

    write ( *, '(a)' ) ' '
    write ( *, '(a,i12)' ) '  Using a seed of ISEED = ', iseed
  
    call set_random_seed ( iseed )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Here are some random iterates:'
    write ( *, '(a)' ) ' '

    do i = 1, n
      call random_number ( harvest = x )
      write ( *, '(g14.6)' ) x
    end do

  end do

  return
end
subroutine test02 ( )

!*******************************************************************************
!
!! TEST02 tests the UNIFORM_01_SAMPLE random number generator.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 February 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer i
  integer iseed
  integer itest
  integer, parameter :: n = 10
  real uniform_01_sample
  real x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  Set the UNIFORM_01_SAMPLE random number seed.'

  do itest = 1, 4

    if ( itest <= 2 ) then
      iseed = 0
    else
      iseed = 1492
    end if

    write ( *, '(a)' ) ' '
    write ( *, '(a,i12)' ) '  Using a seed of ', iseed
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Here are some random iterates:'
    write ( *, '(a)' ) ' '

    do i = 1, n
      x = uniform_01_sample ( iseed )
      write ( *, '(g14.6)' ) x
    end do

  end do

  return
end
subroutine test03 ( )

!*******************************************************************************
!
!! TEST03 tests the FORTRAN90 random number generator.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 February 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: max_n = 10000

  integer log_n
  real mean
  integer n
  real variance
  real x(max_n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  Test the FORTRAN90 random number generator'
  write ( *, '(a)' ) '  by computing mean and variance.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '               N   Mean        Variance'
  write ( *, '(a)' ) ' '

  do log_n = 2, 4

    n = 10**log_n

    call random_number ( harvest = x(1:n) )

    call rvec_mean ( n, x, mean )

    call rvec_variance ( n, x, variance )

    write ( *, '(a,i6,2g14.6)' ) '  Computed: ', n,  mean, variance

  end do

  mean = 0.5E+00
  variance = 1.0E+00 / 12.0E+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6,2g14.6)' ) '  Expected: ', n,  mean, variance

  return
end
subroutine test04 ( )

!*******************************************************************************
!
!! TEST04 tests the UNIFORM_01_SAMPLE random number generator.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 February 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: max_n = 10000

  integer i
  integer iseed
  integer log_n
  real mean
  integer n
  real uniform_01_sample
  real variance
  real x(max_n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  Test the UNIFORM_01_SAMPLE random number generator'
  write ( *, '(a)' ) '  by computing mean and variance.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '               N   Mean        Variance'
  write ( *, '(a)' ) ' '

  iseed = huge ( iseed ) / 2

  do log_n = 2, 4

    n = 10**log_n

    do i = 1, n
      x(i) = uniform_01_sample ( iseed )
    end do

    call rvec_mean ( n, x, mean )

    call rvec_variance ( n, x, variance )

    write ( *, '(a,i6,2g14.6)' ) '  Computed: ', n,  mean, variance

  end do

  mean = 0.5E+00
  variance = 1.0E+00 / 12.0E+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6,2g14.6)' ) '  Expected: ', n,  mean, variance

  return
end
subroutine set_random_seed ( iseed )

!*******************************************************************************
!
!! SET_RANDOM_SEED initializes the FORTRAN90 random number generator.
!
!  Discussion:
!
!    If ISEED is nonzero, then that value is used to construct a seed.
!
!    If ISEED is zero, then the seed is determined by calling the date 
!    and time routine.  Thus, if the code is run at different times, 
!    different seed values will be set.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 May 2002
!
!  Author:
!
!    Lili Ju
!
!  Parameters:
!
!    Input, integer ISEED, is nonzero for a user seed, or 0 if the
!    seed should be determined by this routine.
!
  implicit none

  integer date_time(8)
  logical, parameter :: debug = .false.
  integer i
  integer iseed
  integer j
  integer k
  integer, parameter :: myrank = 0
  integer, allocatable :: seed(:)
!
!  Initialize the random seed routine.
!
  call random_seed ( )
!
!  Request the size of a typical seed.
!  (On the DEC ALPHA, K is returned as 2.)
!
  call random_seed ( size = k )

  if ( debug ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SET_RANDOM_SEED:'
    write ( *, '(a,i6)' ) '  Random seed size is K = ', k
  end if
!
!  Set up space for a seed vector.
!
  allocate ( seed(k) )

  if ( iseed /= 0 ) then

    seed(1:k) = iseed

  else
!
!  Make up a "random" value based on date and time information.
!
    call date_and_time ( values = date_time )

    do i = 1, k

      seed(i) = 0

      do j = 1, 8
        seed(i) = seed(i) + ( j + i ) * date_time(j) + myrank * 100
        seed(i) = ishftc ( seed(i), 4 * ( j - 1 ) )
      end do

    end do

  end if

  if  ( debug ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SET_RANDOM_SEED:'
    write ( *, '(a)' ) '  The random seed vector:'
    write ( *, '(a)' ) ' '

    do i = 1, k
      write ( *, '(i12)' ) seed(i)
    end do

  end if
!
!  Send this random value back to the RANDOM_SEED routine, to be
!  used as the seed of the random number generator.
!
  call random_seed ( put = seed(1:k) )

  deallocate ( seed )

  return
end
subroutine rvec_mean ( n, a, mean )

!*******************************************************************************
!
!! RVEC_MEAN returns the mean of a real vector.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input, real A(N), the vector whose mean is desired.
!
!    Output, real MEAN, the mean, or average, of the vector entries.
!
  implicit none

  integer n

  real a(n)
  real mean

  mean = sum ( a(1:n) ) / real ( n )

  return
end
subroutine rvec_variance ( n, a, variance )

!*******************************************************************************
!
!! RVEC_VARIANCE returns the variance of a real vector.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input, real A(N), the vector whose variance is desired.
!
!    Output, real VARIANCE, the variance of the vector entries.
!
  implicit none

  integer n

  real a(n)
  real mean
  real variance

  call rvec_mean ( n, a, mean )

  variance = sum ( ( a(1:n) - mean )**2 )

  if ( n > 1 ) then
    variance = variance / real ( n - 1 )
  else
    variance = 0.0E+00
  end if

  return
end
function uniform_01_sample ( iseed )

!*******************************************************************************
!
!! UNIFORM_01_SAMPLE is a portable random number generator.
!
!  Formula:
!
!    ISEED = ISEED * (7**5) mod ( 2**31 - 1 )
!    UNIFORM_01_SAMPLE = ISEED * / ( 2**31 - 1 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ISEED, the integer "seed" used to generate
!    the output random number, and updated in preparation for the
!    next one.  ISEED should not be zero.
!
!    Output, real UNIFORM_01_SAMPLE, a random value between 0 and 1.
!
!  Local Parameters:
!
!    IA = 7**5
!    IB = 2**15
!    IB16 = 2**16
!    IP = 2**31-1
!
  implicit none

  integer, parameter :: ia = 16807
  integer, parameter :: ib15 = 32768
  integer, parameter :: ib16 = 65536
  integer, parameter :: ip = 2147483647
  integer iprhi
  integer iseed
  integer ixhi
  integer k
  integer leftlo
  integer loxa
  real uniform_01_sample
!
!  Don't let ISEED be 0.
!
  if ( iseed == 0 ) then
    iseed = ip / 2
  end if
!
!  Get the 15 high order bits of ISEED.
!
  ixhi = iseed / ib16
!
!  Get the 16 low bits of ISEED and form the low product.
!
  loxa = ( iseed - ixhi * ib16 ) * ia
!
!  Get the 15 high order bits of the low product.
!
  leftlo = loxa / ib16
!
!  Form the 31 highest bits of the full product.
!
  iprhi = ixhi * ia + leftlo
!
!  Get overflow past the 31st bit of full product.
!
  k = iprhi / ib15
!
!  Assemble all the parts and presubtract IP.  The parentheses are
!  essential.
!
  iseed = ( ( ( loxa - leftlo * ib16 ) - ip ) + ( iprhi - k * ib15 ) * ib16 ) &
    + k
!
!  Add IP back in if necessary.
!
  if ( iseed < 0 ) then
    iseed = iseed + ip
  end if
!
!  Multiply by 1 / (2**31-1).
!
  uniform_01_sample = real ( iseed ) * 4.656612875E-10

  return
end
subroutine timestamp ( )

!*******************************************************************************
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
!    26 February 2005
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

  write ( *, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
