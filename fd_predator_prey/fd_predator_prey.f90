program main

!*****************************************************************************80
!
!! FD_PREDATOR_PREY solves a pair of predator-prey ODE's.
!
!  Discussion:
!
!    The physical system under consideration is a pair of animal populations.
!
!    The PREY reproduce rapidly; for each animal alive at the beginning of the
!    year, two more will be born by the end of the year.  The prey do not have
!    a natural death rate; instead, they only die by being eaten by the predator.
!    Every prey animal has 1 chance in 1000 of being eaten in a given year by
!    a given predator.
!
!    The PREDATORS only die of starvation, but this happens very quickly.
!    If unfed, a predator will tend to starve in about 1/10 of a year.
!    On the other hand, the predator reproduction rate is dependent on
!    eating prey, and the chances of this depend on the number of available prey.
!
!    The resulting differential equations can be written:
!
!      PREY(0) = 5000
!      PRED(0) =  100
!
!      d PREY / dT =    2 * PREY(T) - 0.001 * PREY(T) * PRED(T)
!      d PRED / dT = - 10 * PRED(T) + 0.002 * PREY(T) * PRED(T)
!
!    Here, the initial values (5000,100) are a somewhat arbitrary starting point.
!
!    The pair of ordinary differential equations that result have an interesting
!    behavior.  For certain choices of the interaction coefficients (such as
!    those given here), the populations of predator and prey will tend to
!    a periodic oscillation.  The two populations will be out of phase; the number
!    of prey will rise, then after a delay, the predators will rise as the prey
!    begins to fall, causing the predator population to crash again.
!
!    In this program, the pair of ODE's is solved with a simple finite difference
!    approximation using a fixed step size.  In general, this is NOT an efficient
!    or reliable way of solving differential equations.  However, this program is
!    intended to illustrate the ideas of finite difference approximation.
!
!    In particular, if we choose a fixed time step size DT, then a derivative
!    such as dPREY/dT is approximated by:
!
!      d PREY / dT = approximately ( PREY(T+DT) - PREY(T) ) / DT
!
!    which means that the first differential equation can be written as
!
!      PREY(T+DT) = PREY(T) + DT * ( 2 * PREY(T) - 0.001 * PREY(T) * PRED(T) ).
!
!    We can rewrite the equation for PRED as well.  Then, since we know the
!    values of PRED and PREY at time 0, we can use these finite difference
!    equations to estimate the values of PRED and PREY at time DT.  These values
!    can be used to get estimates at time 2*DT, and so on.  To get from time
!    T_START = 0 to time T_STOP = 5, we simply take STEP_NUM steps each of size
!    DT = ( T_STOP - T_START ) / STEP_NUM.
!
!    Because finite differences are only an approximation to derivatives, this
!    process only produces estimates of the solution.  And these estimates tend
!    to become more inaccurate for large values of time.  Usually, we can reduce
!    this error by decreasing DT and taking more, smaller time steps.
!
!    In this example, for instance, taking just 100 steps gives nonsensical
!    answers.  Using STEP_NUM = 1000 gives an approximate solution that seems
!    to have the right kind of oscillatory behavior, except that the amplitude
!    of the waves increases with each repetition.  Using 10000 steps, the
!    approximation begins to become accurate enough that we can see that the
!    waves seem to have a fixed period and amplitude.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 June 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    George Lindfield, John Penny,
!    Numerical Methods Using MATLAB,
!    Second Edition,
!    Prentice Hall, 1999,
!    ISBN: 0-13-012641-1,
!    LC: QA297.P45.
!
!  Parameters:
!
!    Input, integer STEP_NUM, the number of steps.
!
  implicit none

  integer ( kind = 4 ) arg_num
  real ( kind = 8 ) dt
  character ( len = 80 ) filename
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iarg
  integer ( kind = 4 ) iargc
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) length
  real ( kind = 8 ) pred_init
  real ( kind = 8 ) prey_init
  integer ( kind = 4 ) step_num
  character ( len = 80 ) string
  real ( kind = 8 ), allocatable :: trf(:,:)
  real ( kind = 8 ) t_start
  real ( kind = 8 ) t_stop

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FD_PREDATOR_PREY'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A finite difference approximate solution of a pair'
  write ( *, '(a)' ) '  of ordinary differential equations for a population'
  write ( *, '(a)' ) '  of predators and prey.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The exact solution shows wave behavior, with a fixed'
  write ( *, '(a)' ) '  period and amplitude.  The finite difference approximation'
  write ( *, '(a)' ) '  can provide a good estimate for this behavior if the stepsize'
  write ( *, '(a)' ) '  DT is small enough.'
!
!  STEP_NUM is an input argument or else read from the user interactively.
!
  arg_num = iargc ( )

  if ( 1 <= arg_num ) then

    iarg = 1
    call getarg ( iarg, string )
    call s_to_i4 ( string, step_num, ierror, length )

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FD_PREDATOR_PREY:'
    write ( *, '(a)' ) '  Please enter the number of time steps:'

    read ( *, '(a)' ) step_num

  end if

  t_start = 0.0D+00
  t_stop =  5.0D+00
  dt = ( t_stop - t_start ) / real ( step_num, kind = 8 )
  prey_init = 5000.0D+00
  pred_init =  100.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Initial time    = ', t_start
  write ( *, '(a,g14.6)' ) '  Final time      = ', t_stop
  write ( *, '(a,g14.6)' ) '  Initial prey    = ', prey_init
  write ( *, '(a,g14.6)' ) '  Initial pred    = ', pred_init
  write ( *, '(a,i8)' ) '  Number of steps = ', step_num
!
!  TRF(1,1:STEP_NUM+1) contains TIME values for each step.
!  TRF(2,1:STEP_NUM+1) contains PREY values
!  TRF(3,1:STEP_NUM+1) contains PREDATOR values
!
  allocate ( trf(1:3,1:step_num+1) )

  trf(1,1)    = t_start
  trf(2,1) = prey_init
  trf(3,1) =  pred_init

  do i = 1, step_num
    trf(1,i+1) = trf(1,i) + dt
    trf(2,i+1) = trf(2,i) + dt &
      * (    2.0D+00 * trf(2,i) - 0.001D+00 * trf(2,i) * trf(3,i) )
    trf(3,i+1) = trf(3,i) + dt &
      * ( - 10.0D+00 * trf(3,i) + 0.002D+00 * trf(2,i) * trf(3,i) )
  end do
!
!  Write data to files.
!
  write ( filename, '(a,i8,a)' ) 'trf_', step_num, '.txt'
  call s_blank_delete ( filename )

  call r8mat_write ( filename, 3, step_num + 1, trf )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  TRF values written to  "' // trim ( filename ) // '".'
!
!  Free memory.
!
  deallocate ( trf )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FD_PREDATOR_PREY'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  return
end
subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is an integer between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is an integer between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IUNIT, the free unit number.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  logical lopen

  iunit = 0

  do i = 1, 99

    if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then

      inquire ( unit = i, opened = lopen, iostat = ios )

      if ( ios == 0 ) then
        if ( .not. lopen ) then
          iunit = i
          return
        end if
      end if

    end if

  end do

  return
end
subroutine r8mat_write ( output_filename, m, n, table )

!*****************************************************************************80
!
!! R8MAT_WRITE writes an R8MAT file.
!
!  Discussion:
!
!    An R8MAT is an array of R8 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) OUTPUT_FILENAME, the output file name.
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, real ( kind = 8 ) TABLE(M,N), the data.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) j
  character ( len = * )  output_filename
  integer ( kind = 4 ) output_status
  integer ( kind = 4 ) output_unit
  character ( len = 30 ) string
  real ( kind = 8 ) table(m,n)
!
!  Open the file.
!
  call get_unit ( output_unit )

  open ( unit = output_unit, file = output_filename, &
    status = 'replace', iostat = output_status )

  if ( output_status /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_WRITE - Fatal error!'
    write ( *, '(a,i8)' ) '  Could not open the output file "' // &
      trim ( output_filename ) // '" on unit ', output_unit
    output_unit = -1
    stop
  end if
!
!  Create a format string.
!
!  For less precision in the output file, try:
!
!                                            '(', m, 'g', 14, '.', 6, ')'
!
  if ( 0 < m .and. 0 < n ) then

    write ( string, '(a1,i8,a1,i8,a1,i8,a1)' ) '(', m, 'g', 24, '.', 16, ')'
!
!  Write the data.
!
    do j = 1, n
      write ( output_unit, string ) table(1:m,j)
    end do

  end if
!
!  Close the file.
!
  close ( unit = output_unit )

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
subroutine s_blank_delete ( s )

!*****************************************************************************80
!
!! S_BLANK_DELETE removes blanks from a string, left justifying the remainder.
!
!  Discussion:
!
!    All TAB characters are also removed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) S, the string to be transformed.
!
  implicit none

  character ch
  integer ( kind = 4 ) get
  integer ( kind = 4 ) put
  character ( len = * )  s
  integer ( kind = 4 ) s_length
  character, parameter :: tab = achar ( 9 )

  put = 0
  s_length = len_trim ( s )

  do get = 1, s_length

    ch = s(get:get)

    if ( ch /= ' ' .and. ch /= tab ) then
      put = put + 1
      s(put:put) = ch
    end if

  end do

  s(put+1:s_length) = ' '

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
