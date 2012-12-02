program main

!*****************************************************************************80
!
!! MAIN is the main program for RANDOM_SEED.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 April 2008
!
!  Author:
!
!    John Burkardt
!
  use omp_lib

  integer ( kind = 4 ), parameter :: n = 20

  integer ( kind = 4 ) i
  integer ( kind = 4 ) proc_num
  integer ( kind = 4 ) seed_in(n)
  integer ( kind = 4 ) seed_out(n)
  integer ( kind = 4 ) thread_num
  real ( kind = 8 ) x(n)

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'RANDOM_SEED'
  write ( *, '(a)' ) '  FORTRAN90/OpenMP version'

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of processors available = ', omp_get_num_procs ( )
  write ( * ,'(a,i8)' ) '  Number of threads =              ', omp_get_max_threads ( )
!
!  Test 1.
!  Generate 10 random numbers sequentially.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Call a simple random number generator which'
  write ( *, '(a)' ) '  does not have any internal memory.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Do this calculation in sequential mode.'

  thread_num = 1

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a)' ) '  Call OMP_SET_NUM_THREADS, and request ', &
    thread_num, ' threads.'

  call omp_set_num_threads ( thread_num )

  call test01 ( n, seed_in, seed_out, x )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          SEED          SEED       R'
  write ( *, '(a)' ) '          IN            OUT'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i12,2x,i12,2x,g14.6)' ) seed_in(i), seed_out(i), x(i)
  end do
!
!  Test 2.
!  Generate 10 random numbers in parallel.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  Call the same simple random number generator which'
  write ( *, '(a)' ) '  does not have any internal memory.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Do this calculation in parallel mode,'

  proc_num = omp_get_num_procs ( )
  thread_num = proc_num

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a)' ) '  Call OMP_SET_NUM_THREADS, and request ', &
    thread_num, ' threads.'

  call omp_set_num_threads ( thread_num )

  call test02 ( n, seed_in, seed_out, x )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          SEED          SEED       R'
  write ( *, '(a)' ) '          IN            OUT'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i12,2x,i12,2x,g14.6)' ) seed_in(i), seed_out(i), x(i)
  end do
!
!  Test 3.
!  Generate 10 random numbers in parallel, carelessly.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  Repeat the previous test, but "carelessly".'
  write ( *, '(a)' ) '  The results depend, in part, on whether the '
  write ( *, '(a)' ) '  threads intefere with each other or not.'

  call test03 ( n, x )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '                                   R'
  write ( *, '(a)' ) '                           '
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,12x,2x,12x,2x,g14.6)' ) x(i)
  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'RANDOM_SEED'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( n, seed_in, seed_out, x )

!*****************************************************************************80
!
!! TEST01 generates N random values in sequentiol mode.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 April 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) seed_in(n)
  integer ( kind = 4 ) seed_out(n)
  real ( kind = 8 ) x(n)

  seed = 1325

  do i = 1, n
    seed_in(i) = seed
    call random_value ( seed, x(i) )
    seed_out(i) = seed
  end do

  return
end
subroutine test02 ( n, seed_in, seed_out, x )

!*****************************************************************************80
!
!! TEST02 generates N random values in parallel mode.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 April 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) seed_vec(n)
  integer ( kind = 4 ) seed_in(n)
  integer ( kind = 4 ) seed_out(n)
  real ( kind = 8 ) x(n)
!
!  Figure out and save the input seed for EVERY iteration.
!
  seed = 1325
  do i = 1, n
    seed_in(i) = seed
    call random_value ( seed, x(i) )
    seed_out(i) = seed
  end do

  do i = 1, n
    seed_vec(i) = seed_in(i)
  end do
!
!  Call the random number generator with the seed vector.
!
!$omp parallel private ( i ) shared ( n, seed_vec, x )
  !$omp do
  do i = 1, n
    call random_value ( seed_vec(i), x(i) )
  end do
  !$omp end do
!$omp end parallel

  return
end
subroutine test03 ( n, x )

!*****************************************************************************80
!
!! TEST03 generates N random values in parallel mode, carelessly.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 April 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x(n)

  seed = 1325
!
!  Call the random number generator with the seed vector.
!
!$omp parallel private ( i ) shared ( n, seed, x )
  !$omp do
  do i = 1, n
    call random_value ( seed, x(i) )
  end do
  !$omp end do
!$omp end parallel

  return
end
subroutine random_value ( seed, r )

!*****************************************************************************80
!
!! RANDOM_VALUE generates a random value R.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 April 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) R, the random value.
!
  implicit none

  real ( kind = 8 ) r
  integer ( kind = 4 ) seed

  seed = mod ( ( 3125 * seed ), 65536 )
  r = ( seed - 32768.0 ) / 16384.0

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
