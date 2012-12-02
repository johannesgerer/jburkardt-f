program main

!*****************************************************************************80
!
!! MAIN is the main program for MEMORY_TEST.
!
!  Discussion:
!
!    MEMORY_TEST declares larger and larger arrays of various types,
!    to see when an error occurs because of memory limitations or other
!    limits.
!
!    Both "allocated" and "automatic" arrays are examined.
!
!  Usage:
!
!    memory_test n_log_min n_log_max
!
!    where
!
!    * n_log_min is the log base 2 of the minimum N;
!    * n_log_max is the log base 2 of the maximum N.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 January 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) arg
  integer ( kind = 4 ) arg_num
  integer ( kind = 4 ) iargc
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) length
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_log
  integer ( kind = 4 ) n_log_max
  integer ( kind = 4 ) n_log_min
  character ( len = 255 ) string

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MEMORY_TEST'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Try to see how big vectors and matrices can be.'
!
!  Get the number of command line arguments.
!
  arg_num = iargc ( )

  if ( arg_num < 1 ) then

    n_log_min = 0
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Using default value N_LOG_MIN = ', n_log_min

  else

    arg = 1
    call getarg ( arg, string )
    call s_to_i4 ( string, n_log_min, ierror, length )
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  User value N_LOG_MIN = ', n_log_min

  end if

  if ( arg_num < 1 ) then

    n_log_max = 27
    write ( *, '(a,i8)' ) '  Using default value N_LOG_MAX = ', n_log_max

  else if ( arg_num == 1 ) then

    n_log_max = n_log_min
    write ( *, '(a,i8)' ) '  Using N_LOG_MAX = N_LOG_MIN = ', n_log_max

  else

    arg = 2
    call getarg ( arg, string )
    call s_to_i4 ( string, n_log_max, ierror, length )
    write ( *, '(a,i8)' ) '  User value N_LOG_MAX = ', n_log_max

  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    'LogN             N   Type  Mode   Ave     CPU        Real          N1      N2'
  write ( *, '(a)' ) ' '

  do n_log = n_log_min, n_log_max

    n = 2**n_log

    call i4vec_memory_allocated_test ( n_log, n )
    call i4vec_memory_automatic_test ( n_log, n )
    call r4vec_memory_allocated_test ( n_log, n )
    call r4vec_memory_automatic_test ( n_log, n )
    call r8vec_memory_allocated_test ( n_log, n )
    call r8vec_memory_automatic_test ( n_log, n )
    call i4mat_memory_allocated_test ( n_log, n )
    call r4mat_memory_allocated_test ( n_log, n )
    call r8mat_memory_allocated_test ( n_log, n )

  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MEMORY_TEST'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine i4mat_memory_allocated_test ( n_log, n )

!*****************************************************************************80
!
!! I4MAT_ALLOCATED_MEMORY_TEST declares and uses an I4MAT of size N.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 January 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N_LOG, the logarithm base 2 of N.
!
!    Input, integer ( kind = 4 ) N, the value of 2^N_LOG.
!
  implicit none

  real ( kind = 4 ) average
  real ( kind = 4 ) cpu_diff
  real ( kind = 4 ) cpu_time1
  real ( kind = 4 ) cpu_time2
  integer ( kind = 4 ) i
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: i4mat
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_log
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n1_log
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n2_log
  real ( kind = 8 ) real_diff
  real ( kind = 8 ) real_time1
  real ( kind = 8 ) real_time2
  real ( kind = 4 ) x

  n1_log = n_log / 2
  n1 = 2**n1_log
  n2_log = n_log - n1_log
  n2 = 2**n2_log

  call real_time ( real_time1 )
  call cpu_time ( cpu_time1 )

  allocate ( i4mat(1:n1,1:n2) )

  do j = 1, n2
    do i = 1, n1
      call random_number ( harvest = x )
      i4mat(i,j) = int ( 3.0E+00 * x )
    end do
  end do

  average = real ( sum ( i4mat(1:n1,1:n2) ), kind = 4 ) &
    / real ( n1, kind = 4 ) / real ( n2, kind = 4 )

  call cpu_time ( cpu_time2 )
  call real_time ( real_time2 )

  cpu_diff = cpu_time2 - cpu_time1
  real_diff = real_time2 - real_time1

  write ( *, '(2x,i2,2x,i12,2x,a5,2x,a4,2x,f4.2,2x,g10.2,2x,g10.2,2x,i6,2x,i6)' ) &
    n_log, n, 'I4MAT', 'ALLO', average, cpu_diff, real_diff, n1, n2

  deallocate ( i4mat )

  return
end
subroutine i4vec_memory_allocated_test ( n_log, n )

!*****************************************************************************80
!
!! I4VEC_MEMORY_ALLOCATED_TEST uses an allocated I4VEC of size N.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 January 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N_LOG, the logarithm base 2 of N.
!
!    Input, integer ( kind = 4 ) N, the value of 2^N_LOG.
!
  implicit none

  real ( kind = 4 ) average
  real ( kind = 4 ) cpu_diff
  real ( kind = 4 ) cpu_time1
  real ( kind = 4 ) cpu_time2
  integer ( kind = 4 ) i
  integer ( kind = 4 ), allocatable, dimension ( : ) :: i4vec
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_log
  real ( kind = 8 ) real_diff
  real ( kind = 8 ) real_time1
  real ( kind = 8 ) real_time2
  real ( kind = 4 ) x

  call real_time ( real_time1 )
  call cpu_time ( cpu_time1 )

  allocate ( i4vec(1:n) )

  do i = 1, n
    call random_number ( harvest = x )
    i4vec(i) = int ( 3.0E+00 * x )
  end do

  average = real ( sum ( i4vec(1:n) ), kind = 4 ) / real ( n, kind = 4 )

  call cpu_time ( cpu_time2 )
  call real_time ( real_time2 )

  cpu_diff = cpu_time2 - cpu_time1
  real_diff = real_time2 - real_time1

  write ( *, '(2x,i2,2x,i12,2x,a5,2x,a4,2x,f4.2,2x,g10.2,2x,g10.2)' ) &
    n_log, n, 'I4VEC', 'ALLO', average, cpu_diff, real_diff

  deallocate ( i4vec )

  return
end
subroutine i4vec_memory_automatic_test ( n_log, n )

!*****************************************************************************80
!
!! I4VEC_MEMORY_AUTOMATIC_TEST uses an automatic I4VEC of size N.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 January 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N_LOG, the logarithm base 2 of N.
!
!    Input, integer ( kind = 4 ) N, the value of 2^N_LOG.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 4 ) average
  real ( kind = 4 ) cpu_diff
  real ( kind = 4 ) cpu_time1
  real ( kind = 4 ) cpu_time2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4vec(n)
  integer ( kind = 4 ) n_log
  real ( kind = 8 ) real_diff
  real ( kind = 8 ) real_time1
  real ( kind = 8 ) real_time2
  real ( kind = 4 ) x

  call real_time ( real_time1 )
  call cpu_time ( cpu_time1 )

  do i = 1, n
    call random_number ( harvest = x )
    i4vec(i) = int ( 3.0E+00 * x )
  end do

  average = real ( sum ( i4vec(1:n) ), kind = 4 ) / real ( n, kind = 4 )

  call cpu_time ( cpu_time2 )
  call real_time ( real_time2 )

  cpu_diff = cpu_time2 - cpu_time1
  real_diff = real_time2 - real_time1

  write ( *, '(2x,i2,2x,i12,2x,a5,2x,a4,2x,f4.2,2x,g10.2,2x,g10.2)' ) &
    n_log, n, 'I4VEC', 'AUTO', average, cpu_diff, real_diff

  return
end
subroutine r4mat_memory_allocated_test ( n_log, n )

!*****************************************************************************80
!
!! R4MAT_MEMORY_ALLOCATED_TEST declares and uses an R4MAT of size N.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 January 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N_LOG, the logarithm base 2 of N.
!
!    Input, integer ( kind = 4 ) N, the value of 2^N_LOG.
!
  implicit none

  real ( kind = 4 ) average
  real ( kind = 4 ) cpu_diff
  real ( kind = 4 ) cpu_time1
  real ( kind = 4 ) cpu_time2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_log
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n1_log
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n2_log
  real ( kind = 4 ), allocatable, dimension ( :, : ) :: r4mat
  real ( kind = 8 ) real_diff
  real ( kind = 8 ) real_time1
  real ( kind = 8 ) real_time2

  n1_log = n_log / 2
  n1 = 2**n1_log
  n2_log = n_log - n1_log
  n2 = 2**n2_log

  call real_time ( real_time1 )
  call cpu_time ( cpu_time1 )

  allocate ( r4mat(1:n1,1:n2) )

  call random_number ( harvest = r4mat )
  r4mat(1:n1,1:n2) = 2.0E+00 * r4mat(1:n1,1:n2)

  average = sum ( r4mat(1:n1,1:n2) ) &
    / real ( n1, kind = 4 ) / real ( n2, kind = 4 )

  call cpu_time ( cpu_time2 )
  call real_time ( real_time2 )

  cpu_diff = cpu_time2 - cpu_time1
  real_diff = real_time2 - real_time1

  write ( *, '(2x,i2,2x,i12,2x,a5,2x,a4,2x,f4.2,2x,g10.2,2x,g10.2,2x,i6,2x,i6)' ) &
    n_log, n, 'R4MAT', 'ALLO', average, cpu_diff, real_diff, n1, n2

  deallocate ( r4mat )

  return
end
subroutine r4vec_memory_allocated_test ( n_log, n )

!*****************************************************************************80
!
!! R4VEC_MEMORY_ALLOCATED_TEST uses an allocated R4VEC of size N.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 January 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N_LOG, the logarithm base 2 of N.
!
!    Input, integer ( kind = 4 ) N, the value of 2^N_LOG.
!
  implicit none

  real ( kind = 4 ) average
  real ( kind = 4 ) cpu_diff
  real ( kind = 4 ) cpu_time1
  real ( kind = 4 ) cpu_time2
  integer ( kind = 4 ) i
  real ( kind = 4 ), allocatable, dimension ( : ) :: r4vec
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_log
  real ( kind = 8 ) real_diff
  real ( kind = 8 ) real_time1
  real ( kind = 8 ) real_time2

  call real_time ( real_time1 )
  call cpu_time ( cpu_time1 )

  allocate ( r4vec(1:n) )

  call random_number ( harvest = r4vec(1:n) )
  r4vec(1:n) = 2.0E+00 * r4vec(1:n)

  average = sum ( r4vec(1:n) ) / real ( n, kind = 4 )

  call cpu_time ( cpu_time2 )
  call real_time ( real_time2 )

  cpu_diff = cpu_time2 - cpu_time1
  real_diff = real_time2 - real_time1

  write ( *, '(2x,i2,2x,i12,2x,a5,2x,a4,2x,f4.2,2x,g10.2,2x,g10.2)' ) &
    n_log, n, 'R4VEC', 'ALLO', average, cpu_diff, real_diff

  deallocate ( r4vec )

  return
end
subroutine r4vec_memory_automatic_test ( n_log, n )

!*****************************************************************************80
!
!! R4VEC_MEMORY_AUTOMATIC_TEST uses an automatic R4VEC of size N.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 January 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N_LOG, the logarithm base 2 of N.
!
!    Input, integer ( kind = 4 ) N, the value of 2^N_LOG.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 4 ) average
  real ( kind = 4 ) cpu_diff
  real ( kind = 4 ) cpu_time1
  real ( kind = 4 ) cpu_time2
  integer ( kind = 4 ) i
  real ( kind = 4 ) r4vec(n)
  integer ( kind = 4 ) n_log
  real ( kind = 8 ) real_diff
  real ( kind = 8 ) real_time1
  real ( kind = 8 ) real_time2

  call real_time ( real_time1 )
  call cpu_time ( cpu_time1 )

  call random_number ( harvest = r4vec(1:n) )
  r4vec(1:n) = 2.0E+00 * r4vec(1:n)

  average = sum ( r4vec(1:n) ) / real ( n, kind = 4 )

  call cpu_time ( cpu_time2 )
  call real_time ( real_time2 )

  cpu_diff = cpu_time2 - cpu_time1
  real_diff = real_time2 - real_time1

  write ( *, '(2x,i2,2x,i12,2x,a5,2x,a4,2x,f4.2,2x,g10.2,2x,g10.2)' ) &
    n_log, n, 'R4VEC', 'AUTO', average, cpu_diff, real_diff

  return
end
subroutine r8mat_memory_allocated_test ( n_log, n )

!*****************************************************************************80
!
!! R8MAT_MEMORY_ALLOCATED_TEST uses an allocated R8MAT of size N.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 January 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N_LOG, the logarithm base 2 of N.
!
!    Input, integer ( kind = 4 ) N, the value of 2^N_LOG.
!
  implicit none

  real ( kind = 8 ) average
  real ( kind = 4 ) cpu_diff
  real ( kind = 4 ) cpu_time1
  real ( kind = 4 ) cpu_time2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_log
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n1_log
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n2_log
  real ( kind = 8 ), allocatable, dimension ( :, : ) :: r8mat
  real ( kind = 8 ) real_diff
  real ( kind = 8 ) real_time1
  real ( kind = 8 ) real_time2

  n1_log = n_log / 2
  n1 = 2**n1_log
  n2_log = n_log - n1_log
  n2 = 2**n2_log

  call real_time ( real_time1 )
  call cpu_time ( cpu_time1 )

  allocate ( r8mat(1:n1,1:n2) )

  call random_number ( harvest = r8mat(1:n1,1:n2) )
  r8mat(1:n1,1:n2) = 2.0E+00 * r8mat(1:n1,1:n2)

  average = sum ( r8mat(1:n1,1:n2) ) &
    / real ( n1, kind = 4 ) / real ( n2, kind = 4 )

  call cpu_time ( cpu_time2 )
  call real_time ( real_time2 )

  cpu_diff = cpu_time2 - cpu_time1
  real_diff = real_time2 - real_time1

  write ( *, '(2x,i2,2x,i12,2x,a5,2x,a4,2x,f4.2,2x,g10.2,2x,g10.2,2x,i6,2x,i6)' ) &
    n_log, n, 'R8MAT', 'ALLO', average, cpu_diff, real_diff, n1, n2

  deallocate ( r8mat )

  return
end
subroutine r8vec_memory_allocated_test ( n_log, n )

!*****************************************************************************80
!
!! R8VEC_MEMORY_ALLOCATED_TEST uses an allocated R8VEC of size N.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 January 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N_LOG, the logarithm base 2 of N.
!
!    Input, integer ( kind = 4 ) N, the value of 2^N_LOG.
!
  implicit none

  real ( kind = 8 ) average
  real ( kind = 4 ) cpu_diff
  real ( kind = 4 ) cpu_time1
  real ( kind = 4 ) cpu_time2
  integer ( kind = 4 ) i
  real ( kind = 8 ), allocatable, dimension ( : ) :: r8vec
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_log
  real ( kind = 8 ) real_diff
  real ( kind = 8 ) real_time1
  real ( kind = 8 ) real_time2

  call real_time ( real_time1 )
  call cpu_time ( cpu_time1 )

  allocate ( r8vec(1:n) )

  call random_number ( harvest = r8vec(1:n) )
  r8vec(1:n) = 2.0D+00 * r8vec(1:n)

  average = sum ( r8vec(1:n) ) / real ( n, kind = 8 )

  call cpu_time ( cpu_time2 )
  call real_time ( real_time2 )

  cpu_diff = cpu_time2 - cpu_time1
  real_diff = real_time2 - real_time1

  write ( *, '(2x,i2,2x,i12,2x,a5,2x,a4,2x,f4.2,2x,g10.2,2x,g10.2)' ) &
    n_log, n, 'R8VEC', 'ALLO', average, cpu_diff, real_diff

  deallocate ( r8vec )

  return
end
subroutine r8vec_memory_automatic_test ( n_log, n )

!*****************************************************************************80
!
!! R8VEC_MEMORY_AUTOMATIC_TEST uses an automatic R8VEC of size N.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 January 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N_LOG, the logarithm base 2 of N.
!
!    Input, integer ( kind = 4 ) N, the value of 2^N_LOG.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) average
  real ( kind = 4 ) cpu_diff
  real ( kind = 4 ) cpu_time1
  real ( kind = 4 ) cpu_time2
  integer ( kind = 4 ) i
  real ( kind = 8 ) r8vec(n)
  integer ( kind = 4 ) n_log
  real ( kind = 8 ) real_diff
  real ( kind = 8 ) real_time1
  real ( kind = 8 ) real_time2

  call real_time ( real_time1 )
  call cpu_time ( cpu_time1 )

  call random_number ( harvest = r8vec(1:n) )
  r8vec(1:n) = 2.0D+00 * r8vec(1:n)

  average = sum ( r8vec(1:n) ) / real ( n, kind = 8 )

  call cpu_time ( cpu_time2 )
  call real_time ( real_time2 )

  cpu_diff = cpu_time2 - cpu_time1
  real_diff = real_time2 - real_time1

  write ( *, '(2x,i2,2x,i12,2x,a5,2x,a4,2x,f4.2,2x,g10.2,2x,g10.2)' ) &
    n_log, n, 'R8VEC', 'AUTO', average, cpu_diff, real_diff

  return
end
subroutine real_time ( seconds )

!*****************************************************************************80
!
!! REAL_TIME returns the real time in seconds.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 November 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) SECONDS, a reading of the real time clock,
!    in seconds.
!
  implicit none

  integer ( kind = 4 ) clock_count
  integer ( kind = 4 ) clock_max
  integer ( kind = 4 ) clock_rate
  real ( kind = 8 ) seconds

  call system_clock ( clock_count, clock_rate, clock_max )

  seconds = real ( clock_count, kind = 8 ) &
    / real ( clock_rate, kind = 8 )

  return
end
subroutine s_to_i4 ( s, value, ierror, length )

!*****************************************************************************80
!
!! S_TO_I4 reads an integer value from a string.
!
!  Discussion:
!
!    Instead of ICHAR, we now use the IACHAR function, which
!    guarantees the ASCII collating sequence.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, a string to be examined.
!
!    Output, integer ( kind = 4 ) VALUE, the value read from the string.
!    If the string is blank, then VALUE will be returned 0.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!    0, no error.
!    1, an error occurred.
!
!    Output, integer ( kind = 4 ) LENGTH, the number of characters
!    of S used to make the integer.
!
  implicit none

  character c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) length
  character ( len = * ) s
  integer ( kind = 4 ) state
  character :: TAB = achar ( 9 )
  integer ( kind = 4 ) value

  value = 0
  ierror = 0
  length = 0

  state = 0
  isgn = 1

  do i = 1, len_trim ( s )

    c = s(i:i)
!
!  STATE = 0, haven't read anything.
!
    if ( state == 0 ) then

      if ( c == ' ' .or. c == TAB ) then

      else if ( c == '-' ) then
        state = 1
        isgn = -1
      else if ( c == '+' ) then
        state = 1
        isgn = +1
      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        state = 2
        value = iachar ( c ) - iachar ( '0' )
      else
        ierror = 1
        return
      end if
!
!  STATE = 1, have read the sign, expecting digits or spaces.
!
    else if ( state == 1 ) then

      if ( c == ' ' .or. c == TAB ) then

      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        state = 2
        value = iachar ( c ) - iachar ( '0' )
      else
        ierror = 1
        return
      end if
!
!  STATE = 2, have read at least one digit, expecting more.
!
    else if ( state == 2 ) then

      if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then

        value = 10 * value + iachar ( c ) - iachar ( '0' )

      else

        value = isgn * value
        ierror = 0
        length = i - 1
        return

      end if

    end if

  end do
!
!  If we read all the characters in the string, see if we're OK.
!
  if ( state == 2 ) then

    value = isgn * value
    ierror = 0
    length = len_trim ( s )

  else

    value = 0
    ierror = 1
    length = 0

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
