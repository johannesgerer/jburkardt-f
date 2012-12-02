program main

!*****************************************************************************80
!
!! MAIN is the main program for TIMESTAMP_PRB.
!
!  Discussion:
!
!    TIMESTAMP_PRB demonstrates the use of TIMESTAMP.
!
!  Modified:
!
!    28 November 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TIMESTAMP_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the TIMESTAMP library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TIMESTAMP_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 demonstrates the use of TIMESTAMP.
!
!  Modified:
!
!    28 November 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  TIMESTAMP prints out the current wallclock time,'
  write ( *, '(a)' ) '  including the year, month, day, hours, minutes,'
  write ( *, '(a)' ) '  seconds, thousandths of a second, and AM/PM.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  This can be useful in keeping track of the date'
  write ( *, '(a)' ) '  of execution of a particular program'
  write ( *, '(a)' ) '  or to give a rough idea of the length of time'
  write ( *, '(a)' ) '  required to run a program.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 demonstrates the use of TIMESTRING.
!
!  Modified:
!
!    28 November 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 40 ) string

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  TIMESTRING returns the current wallclock time,'
  write ( *, '(a)' ) '  including the year, month, day, hours, minutes,'
  write ( *, '(a)' ) '  seconds, thousandths of a second, and AM/PM'
  write ( *, '(a)' ) '  in a string, which the user may print or manipulate.'

  call timestring ( string )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  TIMESTRING returned the value "' &
    // trim ( string ) // '".'

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 demonstrates the use of HMS_CURRENT_PRINT.
!
!  Modified:
!
!    28 November 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer i
  real s
  real t
  real x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  HMS_CURRENT_PRINT prints out the current HMS time'
  write ( *, '(a)' ) '  (hours, minutes, seconds, thousandths of a second,'
  write ( *, '(a)' ) '  and AM/PM) followed by a user-specified string.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  This might be useful in figuring out what routines'
  write ( *, '(a)' ) '  are taking a long time.'

  t = 0.0E+00

  write ( *, '(a)' ) ' '
  call hms_current_print ( '  Call random_number 1000 times.' )

  do i = 1, 10000
    call random_number ( harvest = x )
    t = t + x
  end do

  call hms_current_print ( '  Compute sine of random numbers.' )

  do i = 1, 10000
    call random_number ( harvest = x )
    s = sin ( x )
    t = t + s
  end do

  call hms_current_print ( '  All done!' )

  write ( *, * ) ' '
  write ( *, * ) '  Computed value = ', t

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 demonstrates the use of HMS_DELTA_PRINT.
!
!  Modified:
!
!    20 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer i
  real s
  real t
  real x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  HMS_DELTA_PRINT prints out the delta HMS time'
  write ( *, '(a)' ) '  (hours, minutes, seconds, thousandths of a second,'
  write ( *, '(a)' ) '  and AM/PM) followed by a user-specified string.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  This might be useful in figuring out what routines'
  write ( *, '(a)' ) '  are taking a long time.'

  t = 0.0E+00

  write ( *, '(a)' ) ' '
  call hms_delta_print ( '  Zero out the clock.' )

  do i = 1, 100000
    call random_number ( harvest = x )
    t = t + x
  end do

  call hms_delta_print ( '  Time to compute random numbers.' )

  do i = 1, 100000
    call random_number ( harvest = x )
    s = sin ( x )
    t = t + s
  end do

  call hms_delta_print ( '  Time to compute sine of random numbers.' )

  write ( *, * ) ' '
  write ( *, * ) '  Computed value = ', t

  return
end
