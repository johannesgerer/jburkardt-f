program main

!*****************************************************************************80
!
!! MAIN is the main program for MXV.
!
!  Discussion:
!
!    MXV computes a matrix-vector product in a number of ways, and reports
!    the elapsed CPU time.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 May 2008
!
!  Author:
!
!    John Burkardt
!
!  Usage:
!
!    mxv m n
!
!  Parameters:
!
!    Command line argument, integer M, the number of rows in the matrix.
!
!    Command line argument, integer N, the number of columns in the matrix.
!
  implicit none

  real ( kind = 8 ), allocatable :: a(:,:)
  integer ( kind = 4 ) arg_num
  real ( kind = 8 ) cpu_seconds
  integer ( kind = 4 ) flop_count
  integer ( kind = 4 ) iarg
  integer ( kind = 4 ) iargc
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) last
  integer ( kind = 4 ) m
  real ( kind = 8 ) mflops
  integer ( kind = 4 ) n
  character ( len = 80 ) string
  real ( kind = 8 ), allocatable :: x(:)
  real ( kind = 8 ), allocatable :: y(:)

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MXV:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Compute matrix vector products y = A*x.'
!
!  Get the number of command line arguments.
!
  arg_num = iargc ( )
!
!  Get the number of rows, M.
!
  if ( 1 <= arg_num ) then
  
    iarg = 1
    call getarg ( iarg, string )
    call s_to_i4 ( string, m, ierror, last )
    
  else
  
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter the number of rows, M'
    read ( *, * ) m
    
  end if
!
!  Get the number of columns, N.
!
  if ( 2 <= arg_num ) then
  
    iarg = 2
    call getarg ( iarg, string )
    call s_to_i4 ( string, n, ierror, last )
    
  else
  
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter the number of rows, N'
    read ( *, * ) n
    
  end if
!
!  Record the amount of work.
!  Each of the M entries of Y requires N multiplies and N adds.
!
  flop_count = 2 * m * n

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)'  ) '  Number of matrix rows M =             ', m
  write ( *, '(a,i12)'  ) '  Number of matrix columns N =          ', n
  write ( *, '(a,i12)'  ) '  Number of floating point operations = ', flop_count
!
!  Generate A, X and Y.
!
  allocate ( a(1:m,1:n) )
  allocate ( x(1:n) )
  allocate ( y(1:m) )
!
!  Set A and X.
!
  call matgen ( m, n, a, x )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Method     Cpu Seconds       MegaFlopS'
  write ( *, '(a)' ) '  ------  --------------  --------------'
!
!  DOIDOJ
!
  call mxv_doidoj ( m, n, a, x, y, cpu_seconds )

  if ( 0.0D+00 < cpu_seconds ) then
    mflops = real ( flop_count, kind = 8 ) / cpu_seconds / 1000000.0D+00
  else
    mflops = -1.0D+00
  end if

  write ( *, '(2x,a6,2x,g14.6,2x,g14.6)' ) 'DOIDOJ', cpu_seconds, mflops
!
!  DOJDOI
!
  call mxv_dojdoi ( m, n, a, x, y, cpu_seconds )

  if ( 0.0D+00 < cpu_seconds ) then
    mflops = real ( flop_count, kind = 8 ) / cpu_seconds / 1000000.0D+00
  else
    mflops = -1.0D+00
  end if

  write ( *, '(2x,a6,2x,g14.6,2x,g14.6)' ) 'DOJDOI', cpu_seconds, mflops
!
!  MATMUL
!
  call mxv_matmul ( m, n, a, x, y, cpu_seconds )

  if ( 0.0D+00 < cpu_seconds ) then
    mflops = real ( flop_count, kind = 8 ) / cpu_seconds / 1000000.0D+00
  else
    mflops = -1.0D+00
  end if

  write ( *, '(2x,a6,2x,g14.6,2x,g14.6)' ) 'MATMUL', cpu_seconds, mflops
!
!  Deallocate arrays.
!
  deallocate ( a )
  deallocate ( x )
  deallocate ( y )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MXV:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine matgen ( m, n, a, x )

!*****************************************************************************80
!
!! MATGEN generates a random matrix A and vector X.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 April 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns 
!    of the matrix.
!
!    Output, real ( kind = 8 ) A(M,N), the matrix.
!
!    Output, real ( kind = 8 ) X(N), the vector.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x(n)

  seed = 1325
!
! Set the matrix A.
!
  do j = 1, n
    do i = 1, m
      seed = mod ( ( 3125 * seed ), 65536 )
      a(i,j) = ( seed - 32768.0 ) / 16384.0
    end do
  end do
!
!  Set X.
!
  do i = 1, n
    x(i) = i
  end do

  return
end
subroutine mxv_matmul ( m, n, a, x, y, cpu_seconds )

!*****************************************************************************80
!
!! MXV_MATMUL computes y = A * x, using the FORTRAN90 MATMUL function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 May 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns 
!    of the matrix.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied.
!
!    Output, real ( kind = 8 ) Y(M), the product vector.
!
!    Output, real ( kind = 8 ) CPU_SECONDS, the elapsed CPU time.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) cpu_seconds
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) time1
  real ( kind = 8 ) time2
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(m)

  call cpu_time ( time1 )

  y(1:m) = matmul ( a(1:m,1:n), x(1:n) )

  call cpu_time ( time2 )

  cpu_seconds = time2 - time1

  return
end
subroutine mxv_doidoj ( m, n, a, x, y, cpu_seconds )

!*****************************************************************************80
!
!! MXV_DOIDOJ computes y = A * x, using DO I, DO J loops.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 May 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns
!    of the matrix.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied.
!
!    Output, real ( kind = 8 ) Y(M), the product vector.
!
!    Output, real ( kind = 8 ) CPU_SECONDS, the elapsed CPU time.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) cpu_seconds
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) time1
  real ( kind = 8 ) time2
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(m)

  call cpu_time ( time1 )

  do i = 1, m
    y(i) = 0.0
    do j = 1, n
      y(i) = y(i) + a(i,j) * x(j)
    end do
  end do

  call cpu_time ( time2 )

  cpu_seconds = time2 - time1

  return
end
subroutine mxv_dojdoi ( m, n, a, x, y, cpu_seconds )

!*****************************************************************************80
!
!! MXV_DOJDOI computes y = A * x, using DO J, DO I loops.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 May 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns
!    of the matrix.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied.
!
!    Output, real ( kind = 8 ) Y(M), the product vector.
!
!    Output, real ( kind = 8 ) CPU_SECONDS, the elapsed CPU time.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) cpu_seconds
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) time1
  real ( kind = 8 ) time2
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(m)

  call cpu_time ( time1 )

  do i = 1, m
    y(i) = 0.0
  end do

  do j = 1, n
    do i = 1, m
      y(i) = y(i) + a(i,j) * x(j)
    end do
  end do

  call cpu_time ( time2 )

  cpu_seconds = time2 - time1

  return
end
subroutine s_to_i4 ( s, ival, ierror, length )

!*****************************************************************************80
!
!! S_TO_I4 reads an I4 from a string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, a string to be examined.
!
!    Output, integer ( kind = 4 ) IVAL, the integer value read from the string.
!    If the string is blank, then IVAL will be returned 0.
!
!    Output, integer ( kind = 4 ) IERROR, an error flag.
!    0, no error.
!    1, an error occurred.
!
!    Output, integer ( kind = 4 ) LENGTH, the number of characters of S 
!    used to make IVAL.
!
  implicit none

  character c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) istate
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) length
  character ( len = *  ) s

  ierror = 0
  istate = 0
  isgn = 1
  ival = 0

  do i = 1, len_trim ( s )

    c = s(i:i)
!
!  Haven't read anything.
!
    if ( istate == 0 ) then

      if ( c == ' ' ) then

      else if ( c == '-' ) then
        istate = 1
        isgn = -1
      else if ( c == '+' ) then
        istate = 1
        isgn = + 1
      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        istate = 2
        ival = ichar ( c ) - ichar ( '0' )
      else
        ierror = 1
        return
      end if
!
!  Have read the sign, expecting digits.
!
    else if ( istate == 1 ) then

      if ( c == ' ' ) then

      else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        istate = 2
        ival = ichar ( c ) - ichar ( '0' )
      else
        ierror = 1
        return
      end if
!
!  Have read at least one digit, expecting more.
!
    else if ( istate == 2 ) then

      if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
        ival = 10 * ival + ichar ( c ) - ichar ( '0' )
      else
        ival = isgn * ival
        length = i - 1
        return
      end if

    end if

  end do
!
!  If we read all the characters in the string, see if we're OK.
!
  if ( istate == 2 ) then
    ival = isgn * ival
    length = len_trim ( s )
  else
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
!    May 31 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 May 2001
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
  character ( len = 8 ) date
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  character ( len = 10 )  time
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y
  character ( len = 5 ) zone

  call date_and_time ( date, time, zone, values )

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
