subroutine fdump ( )

!*****************************************************************************80
!
!! FDUMP produces a symbolic dump.
!
!  Discussion:
!
!    This routine is intended to be replaced by a locally written
!    version which produces a symbolic dump.  Failing this,
!    it should be replaced by a version which prints the
!    subprogram nesting list.
!
!    Normally, the dump information should be printed to all the
!    active error output units.  The number and value of these
!    units can be determined by calling XGETUA.
!
!  Modified:
!
!    05 April 2007
!
!  Author:
!
!    Ron Jones
!
!  Reference:
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Technical Report SAND82-0800,
!    Sandia National Laboratories, 1982.
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Software: Practice and Experience,
!    Volume 13, Number 3, 1983, pages 251-257.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    None
!
  implicit none

  return
end
function i1mach ( i )

!*****************************************************************************80
!
!! I1MACH returns integer machine constants.
!
!  Discussion:
!
!    This routine can be used to obtain machine-dependent parameters
!    for the local machine environment.  It is a function
!    subroutine with one (input) argument, and can be called
!    as follows, for example
!
!      K = I1MACH(I)
!
!    where I=1,...,16.  The output value of K above is
!    determined by the input value of I.  The results for
!    various values of I are discussed below.
!
!    I/O unit numbers.
!
!    I1MACH(1) = the standard input unit.
!    I1MACH(2) = the standard output unit.
!    I1MACH(3) = the standard punch unit.
!    I1MACH(4) = the standard error message unit.
!
!    Words.
!
!    I1MACH(5) = the number of bits per integer storage unit.
!    I1MACH(6) = the number of characters per integer storage unit.
!
!    Integers.
!
!    Assume integers are represented in the S digit base A form:
!
!      Sign * (X(S-1)*A^(S-1) + ... + X(1)*A + X(0))
!
!    where 0 <= X(I) < A for I = 0 to S-1.
!
!    I1MACH(7) = A, the base.
!    I1MACH(8) = S, the number of base A digits.
!    I1MACH(9) = A^S-1, the largest integer.
!
!    Floating point numbers
!
!    Assume floating point numbers are represented in the T digit base B form:
!
!      Sign * (B^E) * ((X(1)/B) + ... + (X(T)/B^T) )
!
!    where 0<= X(I) < B for I = 1 to T, 0 < X(1) and EMIN <= E <= EMAX
!
!    I1MACH(10) = B, the base.
!
!    Single precision
!
!    I1MACH(11) = T, the number of base B digits.
!    I1MACH(12) = EMIN, the smallest exponent E.
!    I1MACH(13) = EMAX, the largest exponent E.
!
!    Double precision
!
!    I1MACH(14) = T, the number of base B digits.
!    I1MACH(15) = EMIN, the smallest exponent E.
!    I1MACH(16) = EMAX, the largest exponent E.
!
!  Modified:
!
!    05 April 2007
!
!  Author:
!
!    Phyllis Fox, Andrew Hall, Norman Schryer
!
!  Reference:
!
!    Phyllis Fox, Andrew Hall, Norman Schryer,
!    Algorithm 528:
!    Framework for a Portable Library,
!    ACM Transactions on Mathematical Software,
!    Volume 4, Number 2, June 1978, page 176-188.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the index of the desired constant.
!
!    Output, integer ( kind = 4 ) I1MACH, the value of the constant.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1mach
  integer ( kind = 4 ) imach(16)
!
!  IEEE arithmetic machines, such as the ATT 3B series, Motorola
!  68000 based machines such as the SUN 3 and ATT PC 7300, and
!  8087 based micros such asthe IBM PC and ATT 6300.
!
   data imach( 1) /    5 /
   data imach( 2) /    6 /
   data imach( 3) /    7 /
   data imach( 4) /    6 /
   data imach( 5) /   32 /
   data imach( 6) /    4 /
   data imach( 7) /    2 /
   data imach( 8) /   31 /
   data imach( 9) / 2147483647 /
   data imach(10) /    2 /
   data imach(11) /   24 /
   data imach(12) / -125 /
   data imach(13) /  128 /
   data imach(14) /   53 /
   data imach(15) / -1021 /
   data imach(16) /  1024 /

  if ( i < 1 .or. 16 < i )then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I1MACH - Fatal error!'
    write ( *, '(a,i8)' ) '  I is out of bounds: ', i
    i1mach = 0
    stop
  else
    i1mach = imach(i)
  end if

  return
end
function j4save ( which, value, set )

!*****************************************************************************80
!
!! J4SAVE saves variables needed by the library error handling routines.
!
!  Discussion:
!
!    The internal parameters are initialized to the following values:
!
!    #1 =  0, NERR, the index of the most recent error;
!    #2 =  0, KONTRL, error control flag (0 means only level 2 errors are fatal,
!             and get a printout, while lower level errors get no printout.)
!    #3 =  0, IUNIT, the main error output unit (0 means use standard output).
!    #4 = 10, MAXMES, the maximum number of times any message is printed.
!    #5 =  1, NUNIT, total number of error output units in use.
!    #6 = -1, second error output unit (-1 means not being used).
!    #7 = -1, third error output unit (-1 means not being used).
!    #8 = -1, fourth error output unit (-1 means not being used).
!    #9 = -1, fifth error output unit (-1 means not being used).
!
!  Modified:
!
!    05 April 2007
!
!  Author:
!
!    Ron Jones
!
!  Reference:
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Technical Report SAND82-0800,
!    Sandia National Laboratories, 1982.
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Software: Practice and Experience,
!    Volume 13, Number 3, 1983, pages 251-257.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) WHICH, the index of the item desired.
!    1, NERR, the current error number.
!    2, KONTRL, the current error control flag.
!    3, IUNIT, the current unit number to which error messages are sent.
!       (0 means use standard.)
!    4, MAXMES, the maximum times any message is printed (as set by xermax).
!    5, NUNIT, the number of units to which each error message is written.
!    6, the 2nd unit for error messages.
!    7, the 3rd unit for error messages.
!    8, the 4th unit for error messages.
!    9, the 5th unit for error messages.
!
!    Input, integer ( kind = 4 ) VALUE, the value to be set for the WHICH-th parameter,
!    if SET is TRUE.
!
!    Input, logical SET.
!    TRUE: the WHICH-th parameter will be given the value, VALUE.
!
!    Output, integer ( kind = 4 ) J4SAVE, the old value of the WHICH-th parameter.
!
  implicit none

  integer ( kind = 4 ) j4save
  integer ( kind = 4 ), save, dimension ( 9 ) :: param = (/  &
    0, 2, 0, 10, 1, -1, -1, -1, -1 /)
  logical set
  integer ( kind = 4 ) value
  integer ( kind = 4 ) which

  if ( which < 1 .or. 9 < which ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'J4SAVE - Fatal error!'
    write ( *, '(a,i10)' ) '  Illegal input value of WHICH = ', which
    stop
  end if

  j4save = param(which)

  if ( set ) then
    param(which) = value
  end if

  return
end
function numxer ( nerr )

!*****************************************************************************80
!
!! NUMXER returns the most recent error number.
!
!  Modified:
!
!    05 April 2007
!
!  Author:
!
!    Ron Jones
!
!  Reference:
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Technical Report SAND82-0800,
!    Sandia National Laboratories, 1982.
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Software: Practice and Experience,
!    Volume 13, Number 3, 1983, pages 251-257.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) NERR, the most recent error number.
!
!    Output, integer ( kind = 4 ) NUMXER, the most recent error number.
!
  implicit none

  integer ( kind = 4 ) j4save
  integer ( kind = 4 ) nerr
  integer ( kind = 4 ) numxer
  logical set
  integer ( kind = 4 ) value
  integer ( kind = 4 ) which

  which = 1
  value = 0
  set = .false.

  nerr = j4save ( which, value, set )

  numxer = nerr

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
subroutine xerabt ( messg, nmessg )

!*****************************************************************************80
!
!! XERABT aborts program execution and prints an error message.
!
!  Discussion:
!
!    This routine is called to abort execution of a running program,
!    indicated by the occurrence of a fatal error.
!
!    The error message associated with the fatal condition is provided
!    in the calling sequence.
!
!    This routine is used when the error message handlers XERROR and
!    XERRWV are employed.  The similar routine XERHLT is to be used
!    when the more modern error message handler XERMSG is used.
!
!  Modified:
!
!    05 April 2007
!
!  Author:
!
!    Ron Jones
!
!  Reference:
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Technical Report SAND82-0800,
!    Sandia National Laboratories, 1982.
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Software: Practice and Experience,
!    Volume 13, Number 3, 1983, pages 251-257.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, character ( len = * ) MESSG, the message to be processed.
!
!    Input, integer ( kind = 4 ) NMESSG, the actual number of characters in MESSG.
!    If NMESSG is 0, no message is being supplied.
!
  implicit none

  character ( len = * ) messg
  integer ( kind = 4 ) nmessg

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'XERABT - Termination after fatal error!'

  if ( 0 < len ( messg ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Associated error message:'
    write ( *, '(a)' ) '"' // trim ( messg ) // '"'
  end if

  stop
end
subroutine xerbla ( subrou, nerr )

!*****************************************************************************80
!
!! XERBLA is an error handler for the Level 2 and Level 3 BLAS routines.
!
!  Discussion:
!
!    This routine is called by Level 2 and 3 BLAS routines if an input
!    parameter is invalid.
!
!  Modified:
!
!    06 April 2007
!
!  Parameters:
!
!    Input, character ( len = * ) SUBROU, the name of the routine which
!    called XERBLA.  The name will not be more than 6 characters.
!
!    Input, integer ( kind = 4 ) NERR, the error number, which here is used to
!    indicate the position of the invalid parameter in the
!    parameter-list of the calling routine.
!
  implicit none

  integer ( kind = 4 ) level
  character ( len = 6 ) librar
  character ( len = 60 ) message
  integer ( kind = 4 ) nerr
  character ( len = * ) subrou

  librar = 'SLATEC'
  write ( message, '(a,a,a,i2,a)' ) 'On entry to ', trim ( subrou ), &
    ', parameter number ', nerr, ' had an illegal value.'
  level = 1

  call xermsg ( librar, subrou, message, nerr, level )

  return
end
subroutine xerclr ( )

!*****************************************************************************80
!
!! XERCLR resets the current error number to zero.
!
!  Discussion:
!
!    This routine simply resets the current error number to zero.
!
!    This may be necessary to do in order to determine that
!    a certain error has occurred again since the last time
!    NUMXER was referenced.
!
!  Modified:
!
!    05 April 2007
!
!  Author:
!
!    Ron Jones
!
!  Reference:
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Technical Report SAND82-0800,
!    Sandia National Laboratories, 1982.
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Software: Practice and Experience,
!    Volume 13, Number 3, 1983, pages 251-257.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    None
!
  implicit none

  integer ( kind = 4 ) j4save
  integer ( kind = 4 ) junk
  logical set
  integer ( kind = 4 ) value
  integer ( kind = 4 ) which

  which = 1
  value = 0
  set = .true.

  junk = j4save ( which, value, set )

  return
end
subroutine xercnt ( librar, subrou, messg, nerr, level, kontrl )

!*****************************************************************************80
!
!! XERCNT allows user control over the handling of errors.
!
!  Description:
!
!    This routine allows user control over handling of individual errors.
!
!    This routine is to be used when the error message routine XERMSG
!    is employed.  The similar routine XERCTL is to be used for the
!    older error message routines XERROR and XERRWV.
!
!    Just after each message is recorded, but before it is
!    processed any further (i.e., before it is printed or
!    a decision to abort is made), a call is made to XERCNT.
!
!    If the user has replaced this default, dummy version of XERCNT
!    with a customized routine, it can then be used to override the
!    value of KONTROL used in processing this message by redefining its value.
!
!    KONTRL may be set to any value from -2 to 2.
!
!    The meanings for KONTRL are the same as in XSETF, except
!    that the value of KONTRL changes only for this message.
!
!    If KONTRL is set to a value outside the range from -2 to 2,
!    it will be moved back into that range.
!
!  Reference:
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Technical Report SAND82-0800,
!    Sandia National Laboratories, 1982.
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Software: Practice and Experience,
!    Volume 13, Number 3, 1983, pages 251-257.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, character ( len = * ) LIBRAR, the library or software package
!    from which the error message is coming.
!
!    Input, character ( len = * ) SUBROU, the subroutine or function within
!    the library, from which the error message is coming.
!
!    Input, character ( len = * ) MESSG, the error message.
!
!    Input, integer ( kind = 4 ) NERR, the error number.
!
!    Input, integer ( kind = 4 ) LEVEL, the error severity level.
!    * 2, this is an unconditionally fatal error.
!    * 1, this is a recoverable error.  It is normally non-fatal, unless
!         KONTRL has been reset by XSETF.
!    * 0, this is a warning message only.
!    *-1, this is a warning message which is to be printed at most once,
!         regardless of how many times this call is executed.
!
!    Input/output, integer ( kind = 4 ) KONTRL.  This routine receives the current
!    value of KONTRL, and may reset it.  The change is effective only
!    for the current error message.  This allows the user to suppress
!    or force printing of certain messages, for instance.
!
  implicit none

  integer ( kind = 4 ) kontrl
  integer ( kind = 4 ) level
  character ( len = * ) librar
  character ( len = * ) messg
  integer ( kind = 4 ) nerr
  character ( len = * ) subrou

  return
end
subroutine xerctl ( messg, nmessg, nerr, level, kontrl )

!*****************************************************************************80
!
!! XERCTL allows user control over handling of individual errors.
!
!  Discussion:
!
!    This routine gives the user control over handling of individual errors.
!
!    This routine is to be used when the error message routines XERROR
!    and XERRWV are used.  The similar routine XERCNT is to be used for
!    the newer error message routine XERMSG.
!
!    This routine is called just after each message has been recorded,
!    but before it is processed any further; that is, before the
!    message is printed or a decision to abort is made.
!
!    If the user wishes to influence the behavior of the error package
!    with respect to certain errors, then this dummy version of the
!    routine should be replaced by a routine that carries out the
!    actions the user desires.
!
!    In particular, the user can override the value of KONTRL used
!    in processing this message by redefining its value.
!
!    KONTRL may be set to any value from -2 to 2.
!    The meanings for KONTRL are the same as in XSETF, except
!    that the value of KONTRL changes only for this message.
!
!    If KONTRL is set to a value outside the range from -2 to 2,
!    it will be moved back into that range.
!
!  Modified:
!
!    05 April 2007
!
!  Author:
!
!    Ron Jones
!
!  Reference:
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Technical Report SAND82-0800,
!    Sandia National Laboratories, 1982.
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Software: Practice and Experience,
!    Volume 13, Number 3, 1983, pages 251-257.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, character ( len = * ) MESSG, the error message.
!
!    Input, integer ( kind = 4 ) NMESSG, same as in the call to XERROR or XERRWV.
!
!    Input, integer ( kind = 4 ) NERR, same as in the call to XERROR or XERRWV.
!
!    Input, integer ( kind = 4 ) LEVEL, same as in the call to XERROR or XERRWV.
!
!    Input/output, integer ( kind = 4 ) KONTRL.  On input, the current value of the control
!    flag as set by a call to XSETF.  On output, the new value of kontrl.
!    If KONTRL is not defined, it will remain at its original value.
!    This changed value affects only the current occurrence of the current
!    message.
!
  implicit none

  integer ( kind = 4 ) kontrl
  integer ( kind = 4 ) level
  character ( len = * ) messg
  integer ( kind = 4 ) nerr
  integer ( kind = 4 ) nmessg

  return
end
subroutine xerdmp ( )

!*****************************************************************************80
!
!! XERDMP prints the error tables and then clears them.
!
!  Modified:
!
!    05 April 2007
!
!  Author:
!
!    Ron Jones
!
!  Reference:
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Technical Report SAND82-0800,
!    Sandia National Laboratories, 1982.
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Software: Practice and Experience,
!    Volume 13, Number 3, 1983, pages 251-257.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    None
!
  implicit none

  integer ( kind = 4 ) count
  integer ( kind = 4 ) level
  character ( len = 1 ) messg
  integer ( kind = 4 ) nerr
  integer ( kind = 4 ) nmessg

  messg = ' '
  nmessg = 0
  nerr = 0
  level = 0
  count = 0

  call xersav ( messg, nmessg, nerr, level, count )

  return
end
subroutine xerhlt ( messg )

!*****************************************************************************80
!
!! XERHLT aborts program execution.
!
!  Discussion:
!
!    This routine aborts the execution of the program.
!
!    The error message causing the abort is given in the calling
!    sequence.
!
!    This routine is used when the error message handler XERMSG is
!    employed.  The similar routine XERABT is to be used when the
!    older error message handlers XERROR and XERRWV are used.
!
!  Modified:
!
!    06 April 2007
!
!  Author:
!
!    Ron Jones
!
!  Reference:
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Technical Report SAND82-0800,
!    Sandia National Laboratories, 1982.
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Software: Practice and Experience,
!    Volume 13, Number 3, 1983, pages 251-257.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, character ( len = * ) MESSG, the error message associated
!    with the halt in execution.
!
  implicit none

  character ( len = * ) messg

  stop
end
subroutine xermax ( maxmes )

!*****************************************************************************80
!
!! XERMAX sets the maximum number of times any error message is to be printed.
!
!  Discussion:
!
!    This routine sets the maximum number of times any error message
!    is to be printed.  That is, a non-fatal message associated with
!    a particular numbered error should not be be printed more than
!    MAXMES times.
!
!    Most error messages won't be printed at all if the error printout
!    suppression mode has been set.  That is the case if the variable
!    KONTRL has been set to zero.
!
!  Modified:
!
!    05 April 2007
!
!  Author:
!
!    Ron Jones
!
!  Reference:
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Technical Report SAND82-0800,
!    Sandia National Laboratories, 1982.
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Software: Practice and Experience,
!    Volume 13, Number 3, 1983, pages 251-257.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) MAXMES, the maximum number of times any one message
!    is to be printed.
!
  implicit none

  integer ( kind = 4 ) j4save
  integer ( kind = 4 ) junk
  integer ( kind = 4 ) maxmes
  logical set
  integer ( kind = 4 ) value
  integer ( kind = 4 ) which

  which = 4
  value = maxmes
  set = .true.

  junk = j4save ( which, value, set )

  return
end
subroutine xermsg ( librar, subrou, messg, nerr, level )

!*****************************************************************************80
!
!! XERMSG processes error messages.
!
!  Description:
!
!    This routine processes a diagnostic message in a manner determined by the
!    value of LEVEL and the current value of the library error control
!    flag, KONTRL.
!
!    See subroutine XSETF for details on KONTRL.
!
!  Modified:
!
!    08 April 2007
!
!  Author:
!
!    Kirby Fong
!
!  Reference:
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Technical Report SAND82-0800,
!    Sandia National Laboratories, 1982.
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Software: Practice and Experience,
!    Volume 13, Number 3, 1983, pages 251-257.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, character ( len = * ) LIBRAR, the name of the library from which the
!    error message was generated.
!
!    Input, character ( len = * ) SUBROU, the name of the subroutine or function
!    from which the error message was generated.
!
!    Input, character ( len = * ) MESSG, the text of the error or warning 
!    message.  In the example below, the message is a character constant that
!    contains a generic message.
!
!      call xermsg ('SLATEC', 'MMPY',
!      'The order of the matrix exceeds the row dimension', 3, 1)
!
!    It is possible (and is sometimes desirable) to generate a
!    specific message--e.g., one that contains actual numeric
!    values.  Specific numeric values can be converted into
!    character strings using formatted WRITE statements into
!    character variables.  This is called standard Fortran
!    internal file I/O and is exemplified in the first three
!    lines of the following example.  You can also catenate
!    substrings of characters to construct the error message.
!    Here is an example showing the use of both writing to
!    an internal file and catenating character strings.
!
!      character ( len = 5 ) charn, charl
!      write (charn,'(i5)') n
!      write (charl,'(i5)') lda
!      call xermsg ('SLATEC', 'MMPY', 'The order'//charn//
!      ' of the matrix exceeds its row dimension of'// charl, 3, 1)
!
!    There are two subtleties worth mentioning.  One is that
!    the // for character catenation is used to construct the
!    error message so that no single character constant is
!    continued to the next line.  This avoids confusion as to
!    whether there are trailing blanks at the end of the line.
!    The second is that by catenating the parts of the message
!    as an actual argument rather than encoding the entire
!     message into one large character variable, we avoid
!    having to know how long the message will be in order to
!    declare an adequate length for that large character
!    variable.  XERMSG calls XERPRN to print the message using
!    multiple lines if necessary.  If the message is very long,
!    XERPRN will break it into pieces of 72 characters (as
!    requested by XERMSG) for printing on multiple lines.
!    Also, XERMSG asks XERPRN to prefix each line with ' *  '
!    so that the total line length could be 76 characters.
!    Note also that XERPRN scans the error message backwards
!    to ignore trailing blanks.  Another feature is that
!    the substring '$$' is treated as a new line sentinel
!    by XERPRN.  If you want to construct a multiline
!    message without having to count out multiples of 72
!    characters, just use '$$' as a separator.  '$$'
!    obviously must occur within 72 characters of the
!    start of each line to have its intended effect since
!    XERPRN is asked to wrap around at 72 characters in
!    addition to looking for '$$'.
!
!    Input, integer ( kind = 4 ) NERR, the error number, chosen by the library routine's
!    author.  It must be in the range -99 to 999 (three printable digits).
!    Each distinct error should have its own error number.  These error
!    numbers should be described in the machine readable documentation
!    for the routine.  The error numbers need be unique only within each
!    routine, so it is reasonable for each routine to start enumerating
!    errors from 1 and proceeding to the next integer.
!
!    Input, integer ( kind = 4 ) LEVEL, a value in the range 0 to 2 that indicates the
!    level (severity) of the error.  Their meanings are
!    * -1: A warning message.  This is used if it is not clear
!    that there really is an error, but the user's attention
!    may be needed.  An attempt is made to only print this
!    message once.
!    * 0: A warning message.  This is used if it is not clear
!    that there really is an error, but the user's attention
!    may be needed.
!    * 1: A recoverable error.  This is used even if the error is
!    so serious that the routine cannot return any useful
!    answer.  If the user has told the error package to
!    return after recoverable errors, then XERMSG will
!    return to the Library routine which can then return to
!    the user's routine.  The user may also permit the error
!    package to terminate the program upon encountering a
!    recoverable error.
!
!    * 2: A fatal error.  XERMSG will not return to its caller
!    after it receives a fatal error.  This level should
!    hardly ever be used; it is much better to allow the
!    user a chance to recover.  An example of one of the few
!    cases in which it is permissible to declare a level 2
!    error is a reverse communication Library routine that
!    is likely to be called repeatedly until it integrates
!    across some interval.  If there is a serious error in
!    the input such that another step cannot be taken and
!    the Library routine is called again without the input
!    error having been corrected by the caller, the Library
!    routine will probably be called forever with improper
!    input.  In this case, it is reasonable to declare the
!    error to be fatal.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j4save
  integer ( kind = 4 ) kdummy
  integer ( kind = 4 ) kount
  integer ( kind = 4 ) lerr
  integer ( kind = 4 ) level
  character ( len = 20 )  lfirst
  character ( len = * ) librar
  integer ( kind = 4 ) lkntrl
  integer ( kind = 4 ) llevel
  integer ( kind = 4 ) ltemp
  integer ( kind = 4 ) maxmes
  character ( len = * ) messg
  integer ( kind = 4 ) mkntrl
  integer ( kind = 4 ) nerr
  logical set
  logical skip
  character ( len = * ) subrou
  character ( len = 72 ) temp
  integer ( kind = 4 ) value
  integer ( kind = 4 ) which
  character ( len = 8 ) xlibr
  character ( len = 8 ) xsubr

  which = 2
  value = 0
  set = .false.

  lkntrl = j4save ( which, value, set )

  which = 4
  value = 0
  set = .false.

  maxmes = j4save ( which, value, set )
!
!  LKNTRL is a local copy of the control flag KONTRL.
!
!  MAXMES is the maximum number of times any particular message
!  should be printed.
!
!  We print a fatal error message and terminate for an error in
!  calling XERMSG.  The error number should be positive,
!  and LEVEL should be between 0 and 2.
!
  if ( nerr .lt. -9999999 .or. &
      99999999 .lt. nerr .or. &
      nerr .eq. 0 .or. &
      level .lt. -1 .or. &
      2 .lt. level ) then

    call xerprn ( ' ***', -1, &
     'Fatal error in...$$ ' // &
     'XERMSG -- Invalid error number or level$$ ' // &
     'Job abort due to fatal error.', 72 )

    call xersve ( ' ', ' ', ' ', 0, 0, 0, kdummy )
    call xerhlt ( ' ***XERMSG -- Invalid input' )
    return

  end if
!
!  Record the message.
!
  which = 1
  value = nerr
  set = .true.

  i = j4save ( which, value, set )

  call xersve ( librar, subrou, messg, 1, nerr, level, kount )
!
!  Handle print-once warning messages.
!
  if ( level .eq. -1 .and. 1 .lt. kount ) then
    return
  end if
!
!  Allow temporary user override of the control flag.
!
  xlibr  = librar
  xsubr  = subrou
  lfirst = messg
  lerr   = nerr
  llevel = level
  call xercnt ( xlibr, xsubr, lfirst, lerr, llevel, lkntrl )

  lkntrl = max ( -2, min ( 2, lkntrl ) )
  mkntrl = abs ( lkntrl )
!
!  Skip printing if the control flag value as reset in xercnt is
!  zero and the error is not fatal.
!
  skip = .false.

  if ( level .lt. 2 .and. lkntrl .eq. 0 ) then
    skip = .true.
  end if

  if ( level .eq. 0 .and. maxmes .lt. kount ) then
    skip = .true.
  end if

  if ( level .eq. 1 .and. maxmes .lt. kount .and. mkntrl .eq. 1 ) then
    skip = .true.
  end if

  if ( level .eq. 2 .and. max ( 1, maxmes ) .lt. kount ) then
    skip = .true.
  end if

  if ( .not. skip ) then
!
!  Announce the names of the library and subroutine by building a
!  message in character variable TEMP (not exceeding 66 characters)
!  and sending it out via XERPRN.  Print only if control flag
!  is not zero.
!
  if ( lkntrl .ne. 0 ) then
    temp(1:21) = 'Message from routine '
    i = min ( len ( subrou ), 16 )
    temp(22:21+i) = subrou(1:i)
    temp(22+i:33+i) = ' in library '
    ltemp = 33 + i
    i = min ( len ( librar ), 16)
    temp(ltemp+1:ltemp+i) = librar (1:i)
    temp(ltemp+i+1:ltemp+i+1) = '.'
    ltemp = ltemp + i + 1
    call xerprn ( ' ***', -1, temp(1:ltemp), 72 )
  end if
!
!  If LKNTRL is positive, print an introductory line before
!  printing the message.  The introductory line tells the choice
!  from each of the following three options.
!
!  1.  Level of the message
!
!  'Informative message'
!  'Potentially recoverable error'
!  'Fatal error'
!
!  2.  Whether control flag will allow program to continue
!
!  'Prog continues'
!  'Prog aborted'
!
!  3.  Whether or not a traceback was requested.  (The traceback
!  may not be implemented at some sites, so this only tells
!  what was requested, not what was delivered.)
!
!  'Traceback requested'
!  'Traceback not requested'
!
!  Notice that the line including four prefix characters will not
!  exceed 74 characters.
!  We skip the next block if the introductory line is not needed.
!
  if ( 0 .lt. lkntrl ) then
!
!  The first part of the message tells about the level.
!
    if ( level .le. 0 ) then
      temp(1:20) = 'Informative message,'
      ltemp = 20
    else if ( level .eq. 1 ) then
      temp(1:30) = 'Potentially recoverable error,'
      ltemp = 30
    else
      temp(1:12) = 'Fatal error,'
      ltemp = 12
    end if
!
!  Then whether the program will continue.
!
    if ( ( mkntrl .eq. 2 .and. 1 .le. level ) .or. &
         ( mkntrl .eq. 1 .and. level .eq. 2 ) ) then
      temp(ltemp+1:ltemp+14) = ' Prog aborted,'
      ltemp = ltemp + 14
    else
      temp(ltemp+1:ltemp+16) = ' Prog continues,'
      ltemp = ltemp + 16
    end if
!
!  Finally tell whether there should be a traceback.
!
    if ( 0 .lt. lkntrl ) then
      temp(ltemp+1:ltemp+20) = ' Traceback requested'
      ltemp = ltemp + 20
    else
      temp(ltemp+1:ltemp+24) = ' Traceback not requested'
      ltemp = ltemp + 24
    end if

    call xerprn ( ' ***', -1, temp(1:ltemp), 72 )

  end if
!
!  Now send out the message.
!
  call xerprn ( ' *  ', -1, messg, 72 )
!
!  IF LKNTRL is positive, write the error number and request a
!  traceback.
!
  if ( 0 .lt. lkntrl ) then

    write ( temp, '(a,i8)' ) '  Error number = ', nerr

    call xerprn ( ' *  ', -1, temp, 72 )
    call fdump ( )

  end if
!
!  IF LKNTRL is not zero, print a blank line and an end of message.
!
  if ( lkntrl .ne. 0 ) then
    call xerprn ( ' *  ', -1, ' ', 72 )
    call xerprn ( ' ***', -1, 'End of message', 72 )
    call xerprn ( '    ',  0, ' ', 72 )
  end if
!
!  If the error is not fatal or the error is recoverable and the
!  control flag is set for recovery, then return.
!
  end if

  if ( level .le. 0 .or. &
    ( level .eq. 1 .and. mkntrl .le. 1 ) ) then
    return
  end if
!
!  The program will be stopped due to an unrecovered error or a
!  fatal error.  Print the reason for the abort and the error
!  summary if the control flag and the maximum error count permit.
!
  if ( 0 .lt. lkntrl .and. kount .lt. max ( 1, maxmes ) ) then

    if ( level .eq. 1 ) then
      call xerprn ( ' ***', -1, 'Job abort due to unrecovered error.', 72 )
    else
      call xerprn ( ' ***', -1, 'Job abort due to fatal error.', 72 )
    end if

    call xersve ( ' ', ' ', ' ', -1, 0, 0, kdummy )
    call xerhlt ( ' ' )

  else

    call xerhlt ( messg )

  end if

  return
end
subroutine xerprn ( prefix, npref, messg, nwrap )

!*****************************************************************************80
!
!! XERPRN prints error messages processed by XERMSG.
!
!  Description:
!
!  Discussion:
!
!    This routine is used by the error handling routine XERMSG.  A related
!    routine, XERPRT, is used by the older error handling routines
!    XERROR and XERRWV.
!
!    This routine sends one or more lines to each of the (up to five)
!    logical units to which error messages are to be sent.  This routine
!    is called several times by XERMSG, sometimes with a single line to
!    print and sometimes with a (potentially very long) message that may
!    wrap around into multiple lines.
!
!  Modified:
!
!    04 April 2007
!
!  Author:
!
!    Kirby Fong
!
!  Reference:
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Technical Report SAND82-0800,
!    Sandia National Laboratories, 1982.
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Software: Practice and Experience,
!    Volume 13, Number 3, 1983, pages 251-257.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, character ( len = * ) PREFIX, a string to be put at the beginning of
!    each line before the body of the message.  No more than 16 characters
!    of PREFIX will be used.
!
!    Input, integer ( kind = 4 ) NPREF, the number of characters to use from PREFIX.
!    If it is negative, the intrinsic function LEN is used to determine
!    its length.  If it is zero, PREFIX is not used.  If it exceeds 16 or if
!    LEN(PREFIX) exceeds 16, only the first 16 characters will be
!    used.  If NPREF is positive and the length of PREFIX is less
!    than NPREF, a copy of PREFIX extended with blanks to length
!    NPREF will be used.
!
!    Input, character ( len = * ) MESSG, the error message.  If it is a long
!    message, it will be broken into pieces for printing on multiple lines.
!    Each line starts with the appropriate prefix and be followed by a piece
!    of the message.  NWRAP is the number of characters per piece; that is,
!    after each NWRAP characters, we break and start a new line.  In addition,
!    the characters '$$' embedded in MESSG are a sentinel for a new line.
!    The counting of characters up to NWRAP starts over for each new line.
!    The value of NWRAP typically used by XERMSG is 72 since many
!    older error messages in the SLATEC Library are laid out to rely on
!    wrap-around every 72 characters.
!
!    Input, integer ( kind = 4 ) NWRAP, the maximum size piece into which to break MESSG
!    for printing on multiple lines.  An embedded '$$' ends a line, and the
!    count restarts at the following character.  If a line break does not occur
!    on a blank (it would split a word) that word is moved to the next line.
!    Values of NWRAP less than 16 will be treated as 16.  Values of NWRAP
!    greater than 132 will be treated as 132.  The actual line length will
!    be NPREF + NWRAP after NPREF has been adjusted to fall between 0 and 16
!    and NWRAP has been adjusted to fall between 16 and 132.
!
  implicit none

  character ( len = 148 ) cbuff
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1mach
  integer ( kind = 4 ) idelta
  integer ( kind = 4 ) iu(5)
  integer ( kind = 4 ) lenmsg
  integer ( kind = 4 ) lpiece
  integer ( kind = 4 ) lpref
  integer ( kind = 4 ) lwrap
  character ( len = * ) messg
  integer ( kind = 4 ) n
  character ( len = 2 ), parameter :: newlin = '$$'
  integer ( kind = 4 ) nextc
  integer ( kind = 4 ) npref
  integer ( kind = 4 ) nunit
  integer ( kind = 4 ) nwrap
  character ( len = * ) prefix

  call xgetua ( iu, nunit )
!
!  A zero value for a logical unit number means to use the standard
!  error message unit instead.  I1MACH(4) retrieves the standard
!  error message unit.
!
  n = i1mach(4)
  do i = 1, nunit
    if ( iu(i) == 0 ) then
      iu(i) = n
    end if
  end do
!
!  LPREF is the length of the prefix.  The prefix is placed at the
!  beginning of CBUFF, the character buffer, and kept there during
!  the rest of this routine.
!
  if ( npref < 0 ) then
    lpref = len ( prefix )
  else
    lpref = npref
  end if

  lpref = min ( 16, lpref )

  if ( lpref /= 0 ) then
    cbuff(1:lpref) = prefix
  end if
!
!  LWRAP is the maximum number of characters we want to take at one
!  time from MESSG to print on one line.
!
  lwrap = max ( 16, min ( 132, nwrap ) )
!
!  Set LENMSG to the length of MESSG, ignore any trailing blanks.
!
  lenmsg = len ( messg )
  n = lenmsg
  do i = 1, n
    if ( messg(lenmsg:lenmsg) /= ' ' ) then
      exit
    end if
    lenmsg = lenmsg - 1
  end do
!
!  If the message is all blanks, then print one blank line.
!
  if ( lenmsg == 0 ) then
    cbuff(lpref+1:lpref+1) = ' '
    do i = 1, nunit
      write ( iu(i), '(a)' ) cbuff(1:lpref+1)
    end do
    return
  end if
!
!  Set NEXTC to the position in MESSG where the next substring
!  starts.  From this position we scan for the new line sentinel.
!  When NEXTC exceeds LENMSG, there is no more to print.
!  We loop back to label 50 until all pieces have been printed.
!
!  We look for the next occurrence of the new line sentinel.  The
!  INDEX intrinsic function returns zero if there is no occurrence
!  or if the length of the first argument is less than the length
!  of the second argument.
!
!  There are several cases which should be checked for in the
!  following order.  We are attempting to set LPIECE to the number
!  of characters that should be taken from MESSG starting at
!  position NEXTC.
!
!  * LPIECE == 0
!  The new line sentinel does not occur in the remainder of the
!  character string.  LPIECE should be set to LWRAP or LENMSG+1-NEXTC,
!  whichever is less.
!
!  * LPIECE == 1
!  The new line sentinel starts at MESSG(NEXTC:NEXTC).  LPIECE is effectively
!  zero, and we print nothing to avoid producing unnecessary blank lines.
!  This takes care of the situation where the library routine has a message of
!  exactly 72 characters followed by a new line sentinel followed by more
!  characters.  NEXTC should be incremented by 2.
!
!  * LWRAP + 1 < LPIECE
!  Reduce LPIECE to LWRAP.
!
!  * Otherwise
!  This last case means 2 <= LPIECE <= LWRAP+1.  Reset LPIECE = LPIECE-1.
!  Note that this properly handles the end case where LPIECE = LWRAP+1.
!  That is, the sentinel falls exactly at the end of a line.
!
  nextc = 1

  do

    lpiece = index ( messg(nextc:lenmsg), newlin )
!
!  There was no new line sentinel found.
!
    if ( lpiece == 0 ) then

      idelta = 0
      lpiece = min ( lwrap, lenmsg + 1 - nextc )

      if ( lpiece < lenmsg + 1 - nextc ) then
        do i = lpiece+1, 2, -1
          if ( messg(nextc+i-1:nextc+i-1) == ' ' ) then
            lpiece = i - 1
            idelta = 1
            exit
          end if
        end do
      end if

      cbuff(lpref+1:lpref+lpiece) = messg(nextc:nextc+lpiece-1)
      nextc = nextc + lpiece + idelta
!
!  We have a new line sentinel at MESSG(NEXTC:NEXTC+1).
!  Don't print a blank line.
!
    else if ( lpiece == 1 ) then

      nextc = nextc + 2
      cycle
!
!  LPIECE should be set down to LWRAP.
!
    else if ( lwrap + 1 < lpiece ) then

      idelta = 0
      lpiece = lwrap

      do i = lpiece + 1, 2, -1
        if ( messg(nextc+i-1:nextc+i-1) == ' ' ) then
          lpiece = i - 1
          idelta = 1
          exit
        end if
      end do

      cbuff(lpref+1:lpref+lpiece) = messg(nextc:nextc+lpiece-1)
      nextc = nextc + lpiece + idelta
!
!  If we arrive here, it means 2 <= LPIECE <= LWRAP+1.
!  We should decrement LPIECE by one.
!
    else

      lpiece = lpiece - 1
      cbuff(lpref+1:lpref+lpiece) = messg(nextc:nextc+lpiece-1)
      nextc = nextc + lpiece + 2

    end if
!
!  Print.
!
    do i = 1, nunit
      write ( iu(i), '(a)' ) cbuff(1:lpref+lpiece)
    end do

    if ( lenmsg < nextc ) then
      exit
    end if

  end do

  return
end
subroutine xerprt ( messg, nmessg )

!*****************************************************************************80
!
!! XERPRT prints a message on each file indicated by xgetua.
!
!  Discussion:
!
!    This routine is used by the error handling routines XERROR and XERRWV.
!    A related routine, XERPRN, is used by the more modern error handling
!    routines XERMSG.
!
!  Modified:
!
!    05 April 2007
!
!  Author:
!
!    Ron Jones
!
!  Reference:
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Technical Report SAND82-0800,
!    Sandia National Laboratories, 1982.
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Software: Practice and Experience,
!    Volume 13, Number 3, 1983, pages 251-257.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, character ( len = * ) MESSG, the message to be printed.
!
!    Input, integer ( kind = 4 ) NMESSG, the actual number of characters in MESSG.
!
  implicit none

  integer ( kind = 4 ) ichar
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) kunit
  integer ( kind = 4 ) last
  integer ( kind = 4 ) lun(5)
  character ( len = * ) messg
  integer ( kind = 4 ) messg_len
  integer ( kind = 4 ) nmessg
  integer ( kind = 4 ) nunit
!
!  Obtain unit numbers and write line to each unit.
!
  call xgetua ( lun, nunit )

  messg_len = len ( messg )

  do kunit = 1, nunit

    iunit = lun(kunit)

    do ichar = 1, messg_len, 72

      last = min ( ichar + 71, messg_len )

      if ( iunit == 0 ) then
        write ( *, '(a)' ) messg(ichar:last)
      else
        write ( iunit, '(a)' ) messg(ichar:last)
      end if

    end do

  end do

  return
end
subroutine xerror ( messg, nerr, level )

!*****************************************************************************80
!
!! XERROR processes a diagnostic error message.
!
!  Discussion:
!
!    This routine processes a diagnostic message, in a manner determined
!    by the value of LEVEL and the current value of the library error
!    control flag KONTRL.
!
!    See XSETF for details about KONTRL.
!
!  Example:
!
!    call xerror ( 'SMOOTH -- NUM was zero.', 1, 2 )
!
!    call xerror ( 'INTEG  -- less than full accuracy achieved.', 2, 1 )
!
!    call xerror ( &
!      'ROOTER -- actual zero of f found before interval fully collapsed.',
!      3, 0 )
!
!    call xerror ( 'EXP    -- underflows being set to zero.', 1, -1 )
!
!  Modified:
!
!    05 April 2007
!
!  Author:
!
!    Ron Jones
!
!  Reference:
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Technical Report SAND82-0800,
!    Sandia National Laboratories, 1982.
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Software: Practice and Experience,
!    Volume 13, Number 3, 1983, pages 251-257.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, character ( len = * ) MESSG, the message to be processed.
!
!    Input, integer ( kind = 4 ) NERR, the error number associated with this message.
!    NERR must not be zero.
!
!    Input, integer ( kind = 4 ) LEVEL, the error category.
!    * 2, this is an unconditionally fatal error.
!    * 1, this is a recoverable error.  It is normally non-fatal, unless
!         KONTRL has been reset by XSETF.
!    * 0, this is a warning message only.
!    *-1, this is a warning message which is to be printed at most once,
!         regardless of how many times this call is executed.
!
  implicit none

  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) level
  character ( len = * ) messg
  integer ( kind = 4 ) nerr
  integer ( kind = 4 ) ni
  integer ( kind = 4 ) nmessg
  integer ( kind = 4 ) nr
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2

  nmessg = len ( messg )
  ni = 0
  i1 = 0
  i2 = 0
  nr = 0
  r1 = 0.0D+00
  r2 = 0.0D+00

  call xerrwv ( messg, nmessg, nerr, level, ni, i1, i2, nr, r1, r2 )

  return
end
subroutine xerrwv ( messg, nmessg, nerr, level, ni, i1, i2, nr, r1, r2 )

!*****************************************************************************80
!
!! XERRWV processes an error message that includes numeric information.
!
!  Discussion:
!
!    This routine processes a diagnostic message, in a manner determined
!    by the value of LEVEL and the current value of the library error
!    control flag KONTRL.
!
!    See XSETF for details about KONTRL.
!
!    In addition, up to two integer values and two real values may be
!    printed along with the message.
!
!  Example:
!
!    call xerrwv ( 'SMOOTH -- NUM (=I1) was zero.', 29, 1, 2, 1, num,
!      0, 0, 0.0, 0.0 )
!
!    call xerrwv ( &
!      'QUADXY -- Requested error (R1) less than minimum(R2).', &
!      54, 77, 1, 0, 0, 0, 2, errreq, errmin )
!
!  Modified:
!
!    05 April 2007
!
!  Author:
!
!    Ron Jones
!
!  Reference:
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Technical Report SAND82-0800,
!    Sandia National Laboratories, 1982.
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Software: Practice and Experience,
!    Volume 13, Number 3, 1983, pages 251-257.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, character ( len = * ) MESSG, the message to be processed.
!
!    Input, integer ( kind = 4 ) NMESSG, the number of characters in MESSG.
!
!    Input, integer ( kind = 4 ) NERR, the error number associated with this message.
!    NERR must not be zero.
!
!    Input, integer ( kind = 4 ) LEVEL, the error category.
!    * 2, this is an unconditionally fatal error.
!    * 1, this is a recoverable error.  It is normally non-fatal, unless
!         KONTRL has been reset by XSETF.
!    * 0, this is a warning message only.
!    *-1, this is a warning message which is to be printed at most once,
!         regardless of how many times this call is executed.
!
!    Input, integer ( kind = 4 ) NI, the number of integer values to be printed. (0 to 2)
!
!    Input, integer ( kind = 4 ) I1, I2, the first and second integer values.
!
!    Input, integer ( kind = 4 ) NR, the number of real values to be printed. (0 to 2)
!
!    Input, real ( kind = 8 ) R1, R2, the first and second real values.
!
  implicit none

  character ( len = 37 ) form
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i1mach
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) ifatal
  integer ( kind = 4 ) isizei
  integer ( kind = 4 ) isizef
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) j4save
  integer ( kind = 4 ) junk
  integer ( kind = 4 ) kdummy
  integer ( kind = 4 ) kount
  integer ( kind = 4 ) kunit
  integer ( kind = 4 ) lerr
  integer ( kind = 4 ) level
  integer ( kind = 4 ) lkntrl
  integer ( kind = 4 ) llevel
  integer ( kind = 4 ) lmessg
  integer ( kind = 4 ) lun(5)
  integer ( kind = 4 ) maxmes
  character ( len = * ) messg
  integer ( kind = 4 ) mkntrl
  integer ( kind = 4 ) nerr
  integer ( kind = 4 ) ni
  integer ( kind = 4 ) nmessg
  integer ( kind = 4 ) nr
  integer ( kind = 4 ) nunit
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  logical set
  integer ( kind = 4 ) value
  integer ( kind = 4 ) which
!
!  Get flags
!
  which = 2
  value = 0
  set = .false.

  lkntrl = j4save ( which, value, set )

  which = 4
  value = 0
  set = .false.

  maxmes = j4save ( which, value, set )
!
!  Check for valid input
!
  if ( nmessg <= 0 .or. &
       nerr == 0 .or.   &
       level < -1 .or. &
       2 < level ) then

    if ( 0 < lkntrl ) then
      call xerprt ( 'Fatal error in...', 17 )
    end if

    call xerprt( 'XERROR -- Invalid input', 23 )

    if ( 0 < lkntrl ) then
      call xerprt ( 'Job abort due to fatal error.', 29 )
    end if

    if ( 0 < lkntrl ) then
      call xersav ( ' ', 0, 0, 0, kdummy )
    end if

    call xerabt ( 'XERROR -- invalid input', 23 )
    return

   end if
!
!  Record the message.
!
  which = 1
  value = nerr
  set = .true.

  junk = j4save ( which, value, set )

  call xersav ( messg, nmessg, nerr, level, kount )
!
!  Let the user override.
!
  lmessg = nmessg
  lerr = nerr
  llevel = level

  call xerctl ( messg, lmessg, lerr, llevel, lkntrl )
!
!  Reset to original values.
!
  lmessg = nmessg
  lerr = nerr
  llevel = level

  lkntrl = max ( -2, min ( 2, lkntrl ) )
  mkntrl = abs ( lkntrl )
!
!  Decide whether to print message
!
  if ( llevel < 2 .and. lkntrl == 0 ) then
    go to 100
  end if

  if ( ( llevel == -1 .and. min ( 1, maxmes ) < kount ) .or. &
       ( llevel == 0 .and. maxmes < kount ) .or. &
       ( llevel == 1 .and. maxmes < kount .and. mkntrl == 1 ) .or. &
       ( llevel == 2  .and. max ( 1, maxmes ) < kount ) ) then
    go to 100
  end if

  if ( 0 < lkntrl ) then

    call xerprt ( ' ', 1 )

    if ( llevel == -1 ) then
      call xerprt &
      ( 'Warning message...this message will only be printed once.',57)
    else if ( llevel == 0 ) then
      call xerprt ( 'Warning in...', 13 )
    else if ( llevel == 1 ) then
      call xerprt ( 'Recoverable error in...', 23 )
    else if ( llevel == 2 ) then
      call xerprt ( 'Fatal error in...', 17 )
    end if

  end if
!
!  Message.
!
  call xerprt ( messg, lmessg )
  call xgetua(lun,nunit)
  isizei = 1 + int ( log10 ( real ( i1mach(9), kind = 8 ) ) )
  isizef = 1 + int ( log10 ( real ( i1mach(10), kind = 8 )**i1mach(14) ) )

  do kunit = 1, nunit

    iunit = lun(kunit)

    do i = 1, min ( ni, 2 )
      write (form,21) i,isizei
   21 format ('(11x,21hin above message, i',i1,'=,i',i2,')   ')
      if ( iunit == 0 ) then
        if ( i == 1 ) write (*,form) i1
        if ( i == 2 ) write (*,form) i2
      else
        if ( i == 1 ) write (iunit,form) i1
        if ( i == 2 ) write (iunit,form) i2
      end if
    end do

    do i = 1, min ( nr, 2 )
      write (form,23) i,isizef+10,isizef
   23 format ('(11x,21hin above message, r',i1,'=,e',i2,'.',i2,')')
      if ( iunit == 0 ) then
        if ( i == 1 ) write (*,form) r1
        if ( i == 2 ) write (*,form) r2
      else
        if ( i == 1 ) write (iunit,form) r1
        if ( i == 2 ) write (iunit,form) r2
      end if
    end do
!
!  Print the error number.
!
    if ( 0 < lkntrl ) then

      if ( iunit == 0 ) then
        write ( *, '(a,i10)' ) '  Error number = ', lerr
      else
        write ( iunit, '(a,i10)' ) '  Error number = ', lerr
      end if

    end if

  end do
!
!  Traceback
!
  if ( 0 < lkntrl ) then
    call fdump ( )
  end if

100 continue

  if ( llevel == 2 .or. ( llevel == 1 .and. mkntrl ==2 ) ) then
    ifatal = 1
  else
    ifatal = 0
  end if
!
!  Quit here if message is not fatal.
!
  if ( ifatal <= 0 ) then
    return
  end if
!
!  Print reason for abort and error summary.
!
  if ( 0 < lkntrl .and. kount <= max ( 1, maxmes ) ) then

    if ( llevel == 1 ) then
      call xerprt ( 'Job abort due to unrecovered error.', 35 )
    end if

    if ( llevel == 2 ) then
      call xerprt ( 'Job abort due to fatal error.', 29 )
    end if

    call xersav ( ' ', -1, 0, 0, kdummy )

  end if
!
!  Abort
!
  if ( llevel == 2 .and. max ( 1, maxmes ) < kount ) then
    lmessg = 0
  end if

  call xerabt ( messg, lmessg )

  return
end
subroutine xersav ( messg, nmessg, nerr, level, count )

!*****************************************************************************80
!
!! XERSAV records that an error occurred.
!
!  Modified:
!
!    05 April 2007
!
!  Author:
!
!    Ron Jones
!
!  Reference:
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Technical Report SAND82-0800,
!    Sandia National Laboratories, 1982.
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Software: Practice and Experience,
!    Volume 13, Number 3, 1983, pages 251-257.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, character ( len = * ) MESSG, as in XERROR.
!
!    Input, integer ( kind = 4 ) NMESSG, as in XERROR, except that, when NMESSG = 0,
!    the tables will be dumped and cleared; and when NMESSG < 0,
!    the tables will be dumped, but not cleared.
!
!    Input, integer ( kind = 4 ) NERR, the error number.  NERR should not be 0.
!
!    Input, integer ( kind = 4 ) LEVEL, the error severity level.
!    * 2, this is an unconditionally fatal error.
!    * 1, this is a recoverable error.  It is normally non-fatal, unless
!         KONTRL has been reset by XSETF.
!    * 0, this is a warning message only.
!    *-1, this is a warning message which is to be printed at most once,
!         regardless of how many times this call is executed.
!
!    Output, integer ( kind = 4 ) COUNT, the number of times this message has
!    been seen, or zero if the table has overflowed and does not contain
!    this message specifically.
!    When NMESSG = 0, COUNT will not be altered.
!
  implicit none

  integer ( kind = 4 ) count
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1mach
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ), save, dimension ( 10 ) :: kount = (/ &
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
  integer ( kind = 4 ), save :: kountx = 0
  integer ( kind = 4 ) kunit
  integer ( kind = 4 ) level
  integer ( kind = 4 ), save, dimension ( 10 ) :: levtab = (/ &
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
  integer ( kind = 4 ) lun(5)
  character ( len = 20 ) mes
  character ( len = * ) messg
  character ( len = 20 ), save, dimension ( 10 ) :: mestab
  integer ( kind = 4 ) nerr
  integer ( kind = 4 ), save, dimension ( 10 ) :: nertab = (/ &
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
  integer ( kind = 4 ) nmessg
  integer ( kind = 4 ) nunit
!
!  Dump the table
!
  if ( nmessg <= 0 ) then

     if ( kount(1) == 0 ) then
       return
     end if
!
!  Print to each unit
!
     call xgetua ( lun, nunit )

     do kunit = 1, nunit

       iunit = lun(kunit)
       if ( iunit == 0 ) then
         iunit = i1mach(4)
       end if
!
!  Print table header
!
       write ( iunit, '(a)' ) ' '
       write ( iunit, '(a)' ) &
       '          Error message summary'
       write ( iunit, '(a)' ) &
       'Message start             NERR     Level     Count'
!
!  Print body of table.
!
        do i = 1, 10

          if ( kount(i) == 0 ) then
            exit
          end if

          write ( iunit, '(a20,3i10)' ) &
            mestab(i), nertab(i), levtab(i), kount(i)

        end do
!
!  Print number of other errors.
!
        if ( kountx /= 0 ) then
          write ( iunit, '(a)' ) ' '
          write ( iunit, '(a,i10)' ) &
            'Other errors not individually tabulated = ', kountx
        end if

        write ( iunit, '(a)' ) ' '

     end do

     if ( nmessg < 0 ) then
       return
     end if
!
!  Clear the error tables.
!
    kount(1:10) = 0
    kountx = 0
!
!  Process a message.
!
!  Search for this message, or else an empty slot for this message,
!  or else determine that the error table is full.
!
  else

    mes(1:20) = messg(1:20)

    do i = 1, 10

      ii = i
!
!  An empty slot was found for the new message.
!
      if ( kount(i) == 0 ) then
        mestab(ii) = mes
        nertab(ii) = nerr
        levtab(ii) = level
        kount(ii)  = 1
        count = 1
        return
      end if
!
!  Message found in table.
!
      if ( mes == mestab(i) .and.  &
           nerr == nertab(i) .and. &
           level == levtab(i) ) then
        kount(ii) = kount(ii) + 1
        count = kount(ii)
        return
      end if

    end do
!
!  The table is full.
!
    kountx = kountx + 1
    count = 1

  end if

  return
end
subroutine xersve ( librar, subrou, messg, kflag, nerr, level, icount )

!*****************************************************************************80
!
!! XERSVE records that an error has occurred.
!
!  Discussion:
!
!    This routine is used by the error handling routines associated
!    with XERMSG.  It is a revised version of the routine XERSAV, which
!    was used with the older pair of error handling routines XERROR
!    and XERRWV.
!
!  Modified:
!
!    06 April 2007
!
!  Author:
!
!    Ron Jones
!
!  Reference:
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Technical Report SAND82-0800,
!    Sandia National Laboratories, 1982.
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Software: Practice and Experience,
!    Volume 13, Number 3, 1983, pages 251-257.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, character (len = * ) LIBRAR, the name of the library or software
!    package from which the error message comes.
!
!    Input, character (len = * ) SUBROU, the name of the subroutine or function
!    from which the error message comes.
!
!    Input, character (len = * ) MESSG, the error message.
!
!    Input, integer ( kind = 4 ) KFLAG, indicates the action to be performed.
!    0 < KFLAG, the message in MESSG is saved.
!    KFLAG=0 the tables will be dumped and cleared.
!    KFLAG < 0, the tables will be dumped and not cleared.
!
!    Input, integer ( kind = 4 ) NERR, the error number.
!
!    Input, integer ( kind = 4 ) LEVEL, the error severity.
!
!    Output, integer ( kind = 4 ) ICOUNT, the number of times this message has been seen,
!    or zero if the table has overflowed and does not contain this message
!    specifically.  When KFLAG=0, ICOUNT will not be altered from its
!    input value.
!
  implicit none

  integer ( kind = 4 ), parameter :: lentab = 10

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1mach
  integer ( kind = 4 ) icount
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) kflag
  integer ( kind = 4 ), save, dimension ( lentab ) :: kount
  integer ( kind = 4 ), save :: kountx = 0
  integer ( kind = 4 ) kunit
  integer ( kind = 4 ) level
  integer ( kind = 4 ), save, dimension ( lentab ) :: levtab
  character ( len = 8 ) lib
  character (len = * ) librar
  character ( len = 8 ), save, dimension ( lentab ) :: libtab
  integer ( kind = 4 ) lun(5)
  character ( len = 20 ) mes
  character (len = * ) messg
  character ( len = 20 ), save, dimension ( lentab ) :: mestab
  integer ( kind = 4 ) nerr
  integer ( kind = 4 ), save, dimension ( lentab ) :: nertab
  integer ( kind = 4 ), save :: nmsg = 0
  integer ( kind = 4 ) nunit
  character ( len = 8 ) sub
  character (len = * ) subrou
  character ( len = 8 ), save, dimension ( lentab ) :: subtab

  if ( kflag <= 0 ) then
!
!  Dump the table.
!
    if ( nmsg == 0 ) then
      return
    end if
!
!  Print to each unit.
!
    call xgetua ( lun, nunit )

    do kunit = 1, nunit

      iunit = lun(kunit)
      if ( iunit == 0 ) then
        iunit = i1mach(4)
      end if
!
!  Print the table header.
!
      write ( iunit, '(a)' ) ' '
      write ( iunit, '(a)' ) '          Error message summary'
      write ( iunit, '(a,a)' ) &
        'Library    Subroutine Message start             NERR', &
        '     Level     Count'

!
!  Print body of table.
!
      do i = 1, nmsg
        write ( iunit, '(a,3x,a,3x,a,3i10)' ) &
          libtab(i), subtab(i), mestab(i), nertab(i), levtab(i), kount(i)
      end do
!
!  Print the number of other errors.
!
      if ( kountx /= 0 ) then
        write ( iunit, '(a)' ) ' '
        write ( iunit, '(a,i10)' ) &
          'Other errors not individually tabulated = ', kountx
      end if

      write ( iunit, '(1x)' )

    end do
!
!  Clear the error tables.
!
    if ( kflag == 0 ) then
      nmsg = 0
      kountx = 0
    end if

  else
!
!  Process a message.
!
!  Search for this message, or else an empty slot for this message,
!  or else determine that the error table is full.
!
    lib = librar
    sub = subrou
    mes = messg

    do i = 1, nmsg

      if ( &
        lib == libtab(i) .and.  &
        sub == subtab(i) .and.  &
        mes == mestab(i) .and.  &
        nerr == nertab(i) .and. &
        level == levtab(i) ) then
        kount(i) = kount(i) + 1
        icount = kount(i)
        return
      end if

    end do
!
!  Empty slot found for new message.
!
    if ( nmsg < lentab ) then

      nmsg = nmsg + 1
      libtab(i) = lib
      subtab(i) = sub
      mestab(i) = mes
      nertab(i) = nerr
      levtab(i) = level
      kount(i)  = 1
      icount    = 1
!
!  Table is full.
!
    else

      kountx = kountx + 1
      icount = 0

    end if

  end if

  return
end
subroutine xgetf ( kontrl )

!*****************************************************************************80
!
!! XGETF returns current value of error control flag.
!
!  Discussion:
!
!    This routine returns the current value of KONTRL, the error
!    control flag.
!
!    The amount of output printed for a given error is determined
!    by LEVEL, the level of severity, and KONTRL, which controls
!    how much output the user wants to see.
!
!    The following table shows how each message is treated,
!    depending on the values of KONTRL and LEVEL.
!
!    If KONTRL is zero or negative, no information other than the
!    message itself (including numeric values, if any) will be
!    printed.  If KONTRL is positive, introductory messages,
!    trace-backs, and so on, will be printed in addition to the message.
!
!            KONTRL   0            -1/+1          -2/+2
!        LEVEL
!          2        fatal          fatal          fatal
!          1     not printed      printed         fatal
!          0     not printed      printed        printed
!         -1     not printed      printed once   printed once
!
!  Modified:
!
!    05 April 2007
!
!  Author:
!
!    Ron Jones
!
!  Reference:
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Technical Report SAND82-0800,
!    Sandia National Laboratories, 1982.
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Software: Practice and Experience,
!    Volume 13, Number 3, 1983, pages 251-257.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) KONTRL, the current value of the error control flag.
!
  implicit none

  integer ( kind = 4 ) j4save
  integer ( kind = 4 ) kontrl
  logical set
  integer ( kind = 4 ) value
  integer ( kind = 4 ) which

  which = 2
  value = 0
  set = .false.

  kontrl = j4save ( which, value, set )

  return
end
subroutine xgetua ( iunit, nunit )

!*****************************************************************************80
!
!! XGETUA returns the unit numbers to which error messages are being sent.
!
!  Discussion:
!
!    These unit numbers may have been set by a call to XSETUN,
!    or a call to XSETUA, or may be default values.
!
!  Modified:
!
!    05 April 2007
!
!  Author:
!
!    Ron Jones
!
!  Reference:
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Technical Report SAND82-0800,
!    Sandia National Laboratories, 1982.
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Software: Practice and Experience,
!    Volume 13, Number 3, 1983, pages 251-257.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IUNIT(5), an array into which the routine will
!    store the values of the NUNIT units to which the error messages
!    are being sent.  The value of NUNIT is never more than 5, so
!    using an array of dimension 5 will be sufficient.
!
!    Output, integer ( kind = 4 ) NUNIT, the number of units to which the
!    error messages are being sent.  NUNIT will be in the
!    range from 1 to 5.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) iunit(5)
  integer ( kind = 4 ) j4save
  integer ( kind = 4 ) nunit
  logical set
  integer ( kind = 4 ) value
  integer ( kind = 4 ) which

  which = 5
  value = 0
  set = .false.

  nunit = j4save ( which, value, set)

  if ( nunit < 1 ) then
    return
  end if

  which = 3
  value = 0
  set = .false.

  iunit(1) = j4save ( which, value, set )

  do i = 2, nunit

    which = i + 4
    value = 0
    set = .false.

    iunit(i) = j4save ( which, value, set )

  end do

  return
end
subroutine xgetun ( iunit )

!*****************************************************************************80
!
!! XGETUN returns the (first) output file to which messages are being sent.
!
!  Discussion:
!
!    This routine returns the unit number associated with the first or
!    main file to which error messages are being sent.
!
!    To find out if more than one file is being used for error output,
!    one must use the XGETUA routine.
!
!  Modified:
!
!    05 April 2007
!
!  Author:
!
!    Ron Jones
!
!  Reference:
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Technical Report SAND82-0800,
!    Sandia National Laboratories, 1982.
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Software: Practice and Experience,
!    Volume 13, Number 3, 1983, pages 251-257.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IUNIT, the logical unit number of the first unit
!    to which error messages are being sent.
!
  implicit none

  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) j4save
  logical set
  integer ( kind = 4 ) value
  integer ( kind = 4 ) which

  which = 3
  value = 0
  set = .false.

  iunit = j4save ( which, value, set )

  return
end
subroutine xsetf ( kontrl )

!*****************************************************************************80
!
!! XSETF sets the error control flag.
!
!  Discussion:
!
!    This routine sets the error control flag value to KONTRL.
!
!    The following table shows how each message is treated,
!    depending on the values of KONTRL and LEVEL.  See XERROR
!    for description of LEVEL.
!
!    If KONTRL is zero or negative, no information other than the
!    message itself (including numeric values, if any) will be
!    printed.  If KONTRL is positive, introductory messages,
!    trace-backs, and so on, will be printed in addition to the message.
!
!            KONTRL   0            -1/+1          -2/+2
!        LEVEL
!          2        fatal          fatal          fatal
!          1     not printed      printed         fatal
!          0     not printed      printed        printed
!         -1     not printed      printed once   printed once
!
!  Modified:
!
!    08 April 2007
!
!  Author:
!
!    Ron Jones
!
!  Reference:
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Technical Report SAND82-0800,
!    Sandia National Laboratories, 1982.
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Software: Practice and Experience,
!    Volume 13, Number 3, 1983, pages 251-257.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) KONTRL, the desired value of the error control flag.
!
  implicit none

  integer ( kind = 4 ) j4save
  integer ( kind = 4 ) junk
  integer ( kind = 4 ) kontrl
  integer ( kind = 4 ) level
  integer ( kind = 4 ) nerr
  logical set
  integer ( kind = 4 ) value
  integer ( kind = 4 ) which
  character ( len = 8 ) xern1

  if ( kontrl < -2 .or. 2 < kontrl ) then

    write ( xern1, '(i8)' ) kontrl
    nerr = 1
    level = 2

    call xermsg ( 'XERROR', 'XSETF', &
      'Invalid value of KONTRL = ' // xern1, nerr, level )

    return

  end if

  which = 2
  value = kontrl
  set = .true.

  junk = j4save ( which, value, set )

  return
end
subroutine xsetua ( iunita, nunit )

!*****************************************************************************80
!
!! XSETUA sets up to 5 unit numbers to which messages are to be sent.
!
!  Discussion:
!
!    This routine may be called to declare a list of up to five
!    logical units, each of which is to receive a copy of
!    each error message processed by this package.
!
!    The existence of multiple error output units makes it possible for
!    the user to ensure simultaneous printing of each error message to,
!    say, a main output file, an interactive terminal, and other files
!    such as graphics communication files.
!
!  Modified:
!
!    08 April 2007
!
!  Author:
!
!    Ron Jones
!
!  Reference:
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Technical Report SAND82-0800,
!    Sandia National Laboratories, 1982.
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Software: Practice and Experience,
!    Volume 13, Number 3, 1983, pages 251-257.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IUNIT(NUNIT), unit numbers to which the error messages
!    should be printed.  Normally these numbers should all be different
!    but duplicates are not prohibited.
!
!    Input, integer ( kind = 4 ) NUNIT, the number of unit numbers provided in IUNIT.
!    1 <= N <= 5.
!
  implicit none

  integer ( kind = 4 ) nunit

  integer ( kind = 4 ) i
  integer ( kind = 4 ) iunita(nunit)
  integer ( kind = 4 ) j4save
  integer ( kind = 4 ) junk
  integer ( kind = 4 ) level
  integer ( kind = 4 ) nerr
  logical set
  integer ( kind = 4 ) value
  integer ( kind = 4 ) which
  character ( len = 8 ) xern1

  if ( nunit < 1 .or. 5 < nunit ) then

    write ( xern1, '(i8)' ) nunit
    nerr = 1
    level = 2

    call xermsg ( 'XERROR', 'XSETUA', &
      'Invalid number of units, NUNIT = ' // xern1, nerr, level )

    return

  end if
!
!  Set the main error output unit.
!
  which = 3
  value = iunita(1)
  set = .true.

  junk = j4save ( which, value, set )
!
!  If 1 < NUNIT, set auxilliary output units.
!
  do i = 2, nunit

    which = i + 4
    value = iunita(i)
    set = .true.

    junk = j4save ( which, value, set )

  end do

  which = 5
  value = nunit
  set = .true.

  junk = j4save ( which, value, set )

  return
end
subroutine xsetun ( iunit )

!*****************************************************************************80
!
!! XSETUN sets the output file to which error messages are to be sent.
!
!  Discussion:
!
!    This routine sets the unit number associated with the main error
!    output file.  If auxilliary error output units were defined,
!    this routine suppresses them, as well.
!
!    Note that this error package initializes this unit number to the standard
!    output file, a reasonable choice.
!
!    Common choices for the unit number to be associated with the main
!    error output file might be 1, 6 or 9.  FORTRAN generally requires
!    that such unit numbers be between 1 and 99.  Depending on the
!    compiler, units -1 or 0 might also be legal choices.  It may
!    be the case that some unit numbers cannot be chosen, because
!    they are reserved by the compiler for other purposes.  In
!    particular, some compilers reserve unit 5 for input.
!
!    Copies of the error output may also be sent to up to four auxilliary
!    units, which can be defined by calling XSETUA.
!
!  Modified:
!
!    05 April 2007
!
!  Author:
!
!    Ron Jones
!
!  Reference:
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Technical Report SAND82-0800,
!    Sandia National Laboratories, 1982.
!
!    Ron Jones, David Kahaner,
!    XERROR, The SLATEC Error Handling Package,
!    Software: Practice and Experience,
!    Volume 13, Number 3, 1983, pages 251-257.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IUNIT, the logical unit number to which error
!    messages are to be sent.
!
  implicit none

  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) j4save
  integer ( kind = 4 ) junk
  integer ( kind = 4 ) nunit
  logical set
  integer ( kind = 4 ) value
  integer ( kind = 4 ) which
!
!  Set the main error output unit.
!
  which = 3
  value = iunit
  set = .true.

  junk = j4save ( which, value, set )
!
!  Suppress all the other error output units.
!
  nunit = 1

  which = 5
  value = nunit
  set = .true.

  junk = j4save ( which, value, set )

  return
end
