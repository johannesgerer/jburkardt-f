program main

!*****************************************************************************80
!
!! MAIN is the main program for ISING_3D_SIMULATION.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 April 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) arg_num
  integer ( kind = 4 ) iarg
  integer ( kind = 4 ) iargc
  integer ( kind = 4 ) ierror
  logical, allocatable, dimension ( :, :, : ) :: ising
  integer ( kind = 4 ) iterations
  integer ( kind = 4 ) last
  integer ( kind = 4 ) n
  real ( kind = 8) , dimension ( 0:6 ) :: prob
  character ( len = 80 ) string

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ISING_3D_SIMULATION'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Monte Carlo simulation of a 3D Ising model.'

  arg_num = iargc ( )

  if ( 1 <= arg_num ) then
    iarg = 1
    call getarg ( iarg, string )
    call s_to_i4 ( string, n, ierror, last )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter N, the linear dimension of the problem:'
    read ( *, * ) n
  end if

  if ( 2 <= arg_num ) then
    iarg = 2
    call getarg ( iarg, string )
    call s_to_i4 ( string, iterations, ierror, last )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter the number of iterations to take.'
    read ( *, * ) iterations
  end if

! prob = (/ 0.97D+00, 0.95D+00, 0.85D+00, 0.50D+00, 0.15D+00, 0.05D+00, 0.03D+00 /)
  prob = (/ 0.99D+00, 0.99D+00, 0.99D+00, 0.80D+00, 0.30D+00, 0.20D+00, 0.10D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The linear dimension of the system is N = ', n
  write ( *, '(a,i8)' ) '  The number of iterations taken is ITERATIONS = ', iterations
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The transition probability table, based on the number of'
  write ( *, '(a)' ) '  neighbors with the same spin.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      0         1         2         3         4         5         6'
  write ( *, '(a)' ) ' '
  write ( *, '(7f10.4)' ) prob(0:6)
!
!  Initialize the system.
!
  allocate ( ising(n,n,n) )

  call ising_3d_initialize ( n, ising )
!
!  Do the simulation.
!
  call transition ( n, ising, iterations, prob )
!
!  Free memory.
!
  deallocate ( ising )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ISING_3D_SIMULATION'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine ising_3d_initialize ( n, ising )

!*****************************************************************************80
!
!! ISING_3D_INITIALIZE initializes the Ising array.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 December 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of cells in each dimension.
!
!    Output, logical ISING(N,N,N), the initial Ising array.
!
  implicit none

  integer ( kind = 4 ) n

  logical ising(n,n,n)
  real ( kind = 8 ) r(n,n,n)

  call random_number ( harvest = r )

  ising = ( r <= 0.5D+00 )

  return
end
subroutine transition ( n, ising, iterations, prob )

!*****************************************************************************80
!
!! TRANSITION carries out a Monte Carlo simulation of a 3D Ising model.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 November 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    American National Standard for Programming Language: Fortran - Extended,
!    American National Standards Institute, 1992,
!    pages 296-299.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of cells in each spatial 
!    dimension.
!
!    Input/output, logical ISING(N,N,N).  On input, the current state of the
!    system.  On output, the state of the system after the iterations.
!
!    Input, integer ( kind = 4 ) ITERATIONS, the number of iterations to carry out.
!
!    Input, real ( kind = 8 ) PROB(0:6).  PROB(I) represents the probability 
!    that the spin of a given cell will be reversed, given that it has I immediate 
!    neighbors with spin the same as its own.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) count(n,n,n)
  integer ( kind = 4 ) flip(n,n,n)
  integer ( kind = 4 ) flip_count
  integer ( kind = 4 ) i
  logical ising(n,n,n)
  integer ( kind = 4 ) iterations
  integer ( kind = 4 ) ones(n,n,n)
  integer ( kind = 4 ) ones_count
  real ( kind = 8 ) prob(0:6)
  real ( kind = 8 ) r(n,n,n)
  real ( kind = 8 ) threshhold(n,n,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Step     Positives       Flipped'
  write ( *, '(a)' ) ' '

  flip_count = 0

  ones = 0
  where ( ising ) ones = 1

  ones_count = sum ( ones )

  write ( *, '(2x,i4,2x,i12,2x,i12)' ) 0, ones_count

  do i = 1, iterations
!
!  COUNT contains the number of immediate neighbors of (I,J,K) with a value of 1.
!
    count = cshift ( ones, -1,  1 ) &
          + cshift ( ones, +1,  1 ) &
          + cshift ( ones, -1,  2 ) &
          + cshift ( ones, +1,  2 ) &
          + cshift ( ones, -1,  3 ) &
          + cshift ( ones, +1,  3 )
!
!  Now COUNT contains the number of immediate neighbors of (I,J,K) that are the same
!  as the (I,J,K) value.
!
    where ( .not. ising ) count = 6 - count

    where ( count == 0 ) threshhold = prob(0)
    where ( count == 1 ) threshhold = prob(1)
    where ( count == 2 ) threshhold = prob(2)
    where ( count == 3 ) threshhold = prob(3)
    where ( count == 4 ) threshhold = prob(4)
    where ( count == 5 ) threshhold = prob(5)
    where ( count == 6 ) threshhold = prob(6)

    call random_number ( harvest = r )
!
!  "Flip" the value of (I,J,K) with probability PROB(COUNT).
!  
    where ( r < threshhold )
      flip = 1
    elsewhere
      flip = 0
    endwhere

    flip_count = sum ( flip )

    write ( *, '(2x,4x,2x,12x,2x,i12)' ) flip_count

    where ( flip == 1 ) 
      ising = .not. ising
    endwhere

    where ( ising ) 
      ones = 1
    elsewhere
      ones = 0
    endwhere

    ones_count = sum ( ones )

    write ( *, '(2x,i4,2x,i12)' ) i, ones_count

  end do

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
!    Output, integer ( kind = 4 ) VALUE, the integer value read from the string.
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
