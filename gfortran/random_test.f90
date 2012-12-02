program main

!*****************************************************************************80
!
!! MAIN is the main program for RANDOM_TEST.
!
!  Discussion:
!
!    RANDOM_TEST looks at the FORTRAN90 random number generation utilities.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 May 2000
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'RANDOM_TEST'
  write ( *, '(a)' ) '  FORTRAN90 version.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Examine the RANDOM_NUMBER routine'
  write ( *, '(a)' ) '  for the generation or random numbers.'

  call test01

  call test02

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'RANDOM_TEST'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01

!*****************************************************************************80
!
!! TEST01 demonstrates how the RANDOM_NUMBER routine can be called.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 May 2000
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer i
  real x_array(5,2)
  real x_scalar
  real x_vector(2)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01:'
  write ( *, '(a)' ) '  To generate 10 random numbers, you could'
  write ( *, '(a)' ) '  call RANDOM_NUMBER ( X )...'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  10 times, with a scalar X:'
  write ( *, '(a)' ) ' '

  do i = 1, 10
    call random_number ( x_scalar )
    write ( *, '(a)' ) ' '
    write ( *, '(f8.5)' ) x_scalar
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  5 times with a vector X of length 2:'
  write ( *, '(a)' ) ' '

  do i = 1, 5
    call random_number ( x_vector )
    write ( *, '(a)' ) ' '
    write ( *, '(2f8.5)' ) x_vector(1:2)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  or once with a 5 by 2 X array:'
  write ( *, '(a)' ) ' '

  call random_number ( x_array )

  do i = 1, 5
    write ( *, '(2f8.5)' ) x_array(i,1:2)
  end do

  return
end
subroutine test02

!*****************************************************************************80
!
!! TEST02 demonstrates how the RANDOM_SEED routine can be used.
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
!    John Burkardt
!
  implicit none

  integer i
  integer n
  integer, allocatable, dimension ( : ) :: seed
  integer, allocatable, dimension ( : ) :: seed_save
  real x_array(5,2)
  real x_scalar
  real x_vector(2)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02:'
  write ( *, '(a)' ) '  The RANDOM_SEED routine can be used to restart'
  write ( *, '(a)' ) '  the random number generator.  We will repeat'
  write ( *, '(a)' ) '  the previous test, and manipulate the seed so'
  write ( *, '(a)' ) '  that we compute the same 10 numbers each time.'
!
!  This call initializes the seed.
!
  call random_seed 

  call random_seed ( SIZE = n )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  The size of the seed array is N = ', n

  allocate ( seed(1:n) )
  allocate ( seed_save(1:n) )

  call random_seed ( GET = seed )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The current seed is '
  write ( *, '(5i12)' ) seed(1:n)
  write ( *, '(a)' ) '  We''re going to save this seed and reuse it.'

  seed_save = seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  10 times, with a scalar X:'
  write ( *, '(a)' ) ' '

  do i = 1, 10
    call random_number ( x_scalar )
    write ( *, '(a)' ) ' '
    write ( *, '(f8.5)' ) x_scalar
  end do

  call random_seed ( GET = seed )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The current seed is '
  write ( *, '(5i12)' ) seed(1:n)
  write ( *, '(a)' ) '  We reset the seed to '
  write ( *, '(5i12)' ) seed_save(1:n)
     
  call random_seed ( PUT = seed_save )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  five times with a vector X of length 2:'
  write ( *, '(a)' ) ' '

  do i = 1, 5
    call random_number ( x_vector )
    write ( *, '(a)' ) ' '
    write ( *, '(2f8.5)' ) x_vector(1:2)
  end do

  call random_seed ( GET = seed )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The current seed is '
  write ( *, '(5i12)' ) seed(1:n)
  write ( *, '(a)' ) '  We reset the seed to '
  write ( *, '(5i12)' ) seed_save(1:n)
     
  call random_seed ( PUT = seed_save )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  or once with a 5 by 2 X array:'
  write ( *, '(a)' ) ' '

  call random_number ( x_array )

  do i = 1, 5
    write ( *, '(2f8.5)' ) x_array(i,1:2)
  end do

  call random_seed ( GET = seed )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The current seed is '
  write ( *, '(5i12)' ) seed(1:n)

  deallocate ( seed )
  deallocate ( seed_save )

  return
end
subroutine timestamp ( )

!*****************************************************************************80
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
