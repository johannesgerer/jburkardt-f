program main

!*****************************************************************************80
!
!! MAIN is the main program for XLF_INTRINSICS_PRB.
!
!  Discussion:
!
!    XLF_INTRINSICS_PRB calls some of the IBM XLF intrinsic routines.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 October 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'XLF_INTRINSICS_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the XLF intrinsic library.'

  call test_erf ( )
  call test_erfc ( )
  call test_etime_ ( )
  call test_fdate_ ( )
  call test_gamma ( )
  call test_hostnm_ ( )
  call test_lgamma ( )
  call test_rand ( )
  call test_sizeof ( )
  call test_sleep_ ( )
  call test_system ( )
  call test_time_ ( )
  call test_timef ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'XLF_INTRINSICS_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test_erf ( )

!*****************************************************************************80
!
!! TEST_ERF checks ERF against ERF_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  real fx2
  integer n_data
  real ( kind = 8 ) x
  real x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_ERF:'
  write ( *, '(a)' ) '  ERF computes the error function.'
  write ( *, '(a)' ) '  ERF_VALUES returns selected values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) &
    '      X             FX', &
    '                        FX2'
  write ( *, '(a,a)' ) &
    '                   (table)', &
    '                   (ERF)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call erf_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    x2 = real ( x )
    fx2 = erf ( x2 )

    write ( *, '(2x,f12.8,2x,g24.16,2x,g24.16)' ) &
      x, fx, fx2

  end do

  return
end
subroutine test_erfc ( )

!*****************************************************************************80
!
!! TEST_ERFC checks ERFC against ERFC_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  real fx2
  integer n_data
  real ( kind = 8 ) x
  real x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_ERFC:'
  write ( *, '(a)' ) '  ERFC computes the complementary error function.'
  write ( *, '(a)' ) '  ERFC_VALUES returns selected values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) &
    '      X             FX', &
    '                        FX2'
  write ( *, '(a,a)' ) &
    '                   (table)', &
    '                   (ERFC)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call erfc_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    x2 = real ( x )
    fx2 = erfc ( x2 )

    write ( *, '(2x,f12.8,2x,g24.16,2x,g24.16)' ) &
      x, fx, fx2

  end do

  return
end
subroutine test_etime_ ( )

!*****************************************************************************80
!
!! TEST_ETIME_ tests ETIME_.
!
!  Discussion:
!
!    call etime_ ( tarray, result ) 
!
!    or
!
!    call etime_ ( tarray )
!
!    or
!
!    result = etime_ ( tarray )
!
!    returns in the real array TARRAY(2) the user and system execution times
!    and a "success" code in the integer value RESULT (0 = success).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 October 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer result
  real tarray(2)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_ETIME_'
  write ( *, '(a)' ) '  ETIME_ returns the user and system execution time.'

  call etime_ ( tarray, result )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  TARRAY(1) (User time)   = ', tarray(1)
  write ( *, '(a,g14.6)' ) '  TARRAY(2) (System time) = ', tarray(2)
  write ( *, '(a,i8)' )    '  RESULT    (0 = success) = ', result

  return
end
subroutine test_fdate_ ( )

!*****************************************************************************80
!
!! TEST_FDATE_ tests FDATE_.
!
!  Discussion:
!
!    call fdate_ ( string ) 
!    or
!    string = fdate_ ( )
!    returns the current date as a string.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 October 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 80 ) string

  call fdate_ ( string )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_FDATEa_'
  write ( *, '(a)' ) '  FDATE_ returns the current date as a string.'
  write ( *, '(a)' ) '  It can be called as a function or subroutine.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  CALL FDATE_ ( STRING ) returns:'
  write ( *, '(a)' ) '  "' // trim ( string ) // '".'

  return
end
subroutine test_gamma ( )

!*****************************************************************************80
!
!! TEST_GAMMA checks GAMMA against GAMMA_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  real fx2
  integer n_data
  real ( kind = 8 ) x
  real x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_GAMMA:'
  write ( *, '(a)' ) &
    '  GAMMA computes the gamma function.'
  write ( *, '(a)' ) '  GAMMA_VALUES returns selected values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) &
    '      X             FX', &
    '                        FX2'
  write ( *, '(a,a)' ) &
    '                   (table)', &
    '                   (GAMMA)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call gamma_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    x2 = real ( x )
    fx2 = gamma ( x2 )

    write ( *, '(2x,f12.8,2x,g24.16,2x,g24.16)' ) &
      x, fx, fx2

  end do

  return
end
subroutine test_hostnm_ ( )

!*****************************************************************************80
!
!! TEST_HOSTNM_ tests HOSTNM_.
!
!  Discussion:
!
!    result = hostnm_ ( string )
!
!    returns a success code in the integer RESULT, and sets STRING
!    to a character string of the host name.
!
!    This function is not an IBM intrinsic.  It can be accessed by a 
!    USE statement invoking the XLFUTILITY module.  
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 October 2007
!
!  Author:
!
!    John Burkardt
!
  use xlfutility

  implicit none

  integer result
  character ( len = 63 ) string

  result = hostnm_ ( string )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_HOSTNM'
  write ( *, '(a)' ) '  HOSTNM_ returns the host name as a string.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  RESULT = HOSTNM_ ( STRING ) returns:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  STRING = "' // trim ( string ) // '".'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  RESULT (0 = success) = ', result

  return
end
subroutine test_lgamma ( )

!*****************************************************************************80
!
!! TEST_LGAMMA checks LGAMMA against GAMMA_LOG_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 October 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  real fx2
  integer n_data
  real ( kind = 8 ) x
  real x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_LGAMMA:'
  write ( *, '(a)' ) &
    '  LGAMMA computes the log of the gamma function.'
  write ( *, '(a)' ) '  GAMMA_LOG_VALUES returns selected values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) &
    '      X             FX', &
    '                        FX2'
  write ( *, '(a,a)' ) &
    '                   (table)', &
    '                   (LGAMMA)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call gamma_log_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    x2 = real ( x )
    fx2 = lgamma ( x2 )

    write ( *, '(2x,f12.8,2x,g24.16,2x,g24.16)' ) &
      x, fx, fx2

  end do

  return
end
subroutine test_rand ( )

!*****************************************************************************80
!
!! TEST_RAND tests RAND.
!
!  Discussion:
!
!    The form of the call is
!
!      X = RAND ( )
!
!    To control the generation of random numbers by RAND requires the
!    use of the SRAND function, for setting the seed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 October 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer i
  real r
  real seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_RAND'
  write ( *, '(a)' ) '  R = RAND() returns a real pseudorandom value.'
  write ( *, '(a)' ) ' '

  do i = 1, 10
    r = rand ( )
    write ( *, '(2x,i8,2x,g14.6)' ) i, r
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  To control the sequence, first call SRAND'
  write ( *, '(a)' ) '  with a REAL number as the seed.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SEED = 1.23456789'
  write ( *, '(a)' ) '  call SRAND ( SEED )'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  do i = 1, 10'
  write ( *, '(a)' ) '    r = RAND ( )'
  write ( *, '(a)' ) '  end do'
  write ( *, '(a)' ) ' '

  seed = 1.23456789
  call srand ( seed )

  do i = 1, 10
    r = rand ( )
    write ( *, '(2x,i8,2x,g14.6)' ) i, r
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Changing the seed changes the sequence:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SEED = 9.87654321'
  write ( *, '(a)' ) '  call SRAND ( SEED )'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  do i = 1, 10'
  write ( *, '(a)' ) '    r = RAND ( )'
  write ( *, '(a)' ) '  end do'
  write ( *, '(a)' ) ' '

  seed = 9.87654321
  call srand ( seed )

  do i = 1, 10
    r = rand ( )
    write ( *, '(2x,i8,2x,g14.6)' ) i, r
  end do

  seed = 123456789
  call srand ( seed )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Restoring the old seed restores the old sequence:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SEED = 1.23456789'
  write ( *, '(a)' ) '  call SRAND ( SEED )'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  do i = 1, 10'
  write ( *, '(a)' ) '    r = RAND ( )'
  write ( *, '(a)' ) '  end do'
  write ( *, '(a)' ) ' '

  seed = 1.23456789
  call srand ( seed )

  do i = 1, 10
    r = rand ( )
    write ( *, '(2x,i8,2x,g14.6)' ) i, r
  end do

  return
end
subroutine test_sizeof ( )

!*****************************************************************************80
!
!! TEST_SIZEOF tests SIZEOF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 October 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  complex c
  character ch
  double complex dc
  double precision dp
  integer size
  integer i
  logical l
  real r
  character ( len = 15 ) s15
  character ( len = 22 ) type

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_SIZEOF'
  write ( *, '(a)' ) '  SIZEOF returns the size, in bytes, of various'
  write ( *, '(a)' ) '  objects.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Type                   Byte Size'
  write ( *, '(a)' ) ' '

  type = 'character'
  size = sizeof ( ch )
  write ( *, '(2x,a22,2x,i8)' ) type, size

  type = 'complex'
  size = sizeof ( c )
  write ( *, '(2x,a22,2x,i8)' ) type, size
  
  type = 'double complex'
  size = sizeof ( dc )
  write ( *, '(2x,a22,2x,i8)' ) type, size

  type = 'double precision'
  size = sizeof ( dp )
  write ( *, '(2x,a22,2x,i8)' ) type, size

  type = 'integer'
  size = sizeof ( i )
  write ( *, '(2x,a22,2x,i8)' ) type, size

  type = 'logical'
  size = sizeof ( l )
  write ( *, '(2x,a22,2x,i8)' ) type, size

  type = 'real'
  size = sizeof ( r )
  write ( *, '(2x,a22,2x,i8)' ) type, size

  type = 'character ( len = 15 )'
  size = sizeof ( s15 )
  write ( *, '(2x,a22,2x,i8)' ) type, size

  return
end
subroutine test_sleep_ ( )

!*****************************************************************************80
!
!! TEST_SLEEP_ tests SLEEP_.
!
!  Discussion:
!
!    The routine
!
!      subroutine sleep_ ( t ) 
!
!    causes the program to "sleep" for T seconds.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 October 2007
!
!  Author:
!
!    John Burkardt
!
  use xlfutility

  implicit none

  integer value

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_SLEEP_'
  write ( *, '(a)' ) '  CALL SLEEP_ ( T ) causes the program to "sleep"'
  write ( *, '(a)' ) '  for T seconds.'

  value = time_ ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  Current TIME in seconds is = ', value

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now CALL SLEEP_ ( 10 )'

  value = 2
  call sleep_ ( value )
  
  value = time_ ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  Current TIME in seconds is = ', value

  return
end
subroutine test_system ( )

!*****************************************************************************80
!
!! TEST_SYSTEM tests SYSTEM.
!
!  Discussion:
!
!    The function
!
!      call system ( string, result ) 
!
!    pauses the execution of the FORTRAN90 program, issues the command
!    (or comamnds) contained in STRING to the operating system, and
!    resumes execution once the commands have completed.
!
!    The output of the command will be included in the output stream
!    of the program.
!
!    RESULT is a success/error return code.
!
!    A sequence of commands can be include in one string by using
!    the ";" to separate them.
!
!    If the command is terminated with an ampersand, then the program
!    will NOT wait for the command to complete.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 October 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer result
  character ( len = 80 ) string

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_SYSTEM'
  write ( *, '(a)' ) '  RESULT = SYSTEM ( STRING )'
  write ( *, '(a)' ) '  issues the command STRING to the operating system,'
  write ( *, '(a)' ) '  waits for the command to complete and sets RESULT'
  write ( *, '(a)' ) '  to the return value of the command.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We issue the command'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    CALL SYSTEM ( ''date'', RESULT )'
  write ( *, '(a)' ) ' '

  string = 'date'
  call system ( string, result )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  RESULT ( 0 = success ) = ', result

  return
end
subroutine test_time_ ( )

!*****************************************************************************80
!
!! TEST_TIME_ tests TIME_.
!
!  Discussion:
!
!    The function
!
!      integer function time_ ( ) 
!
!    returns the current GMT time in seconds.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 October 2007
!
!  Author:
!
!    John Burkardt
!
  use xlfutility

  implicit none

  real t1
  real t2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_TIME_'
  write ( *, '(a)' ) '  T = TIME_ ( ) returns the current GMT time, in seconds.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The code fragment:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    T1 = TIME_ ( )'
  write ( *, '(a)' ) '    stuff happens'
  write ( *, '(a)' ) '    T2 = TIME_ ( )'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  allows us to estimate T2-T1 seconds required for'
  write ( *, '(a)' ) '  "stuff happens".'

  t1 = time_ ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a,f12.4,a)' ) '  Value of TIME_ ( ) = T1 = ', t1, ' seconds.'
  call sleep_ ( 2 )
  t2 = time_ ( )
  write ( *, '(a,f12.4,a)' ) '  Value of TIME_ ( ) = T2 = ', t2, ' seconds.'
  write ( *, '(a,f12.4,a)' ) '  Elapsed time T2 - T1 =    ', t2 - t1, ' seconds.'

  return
end
subroutine test_timef ( )

!*****************************************************************************80
!
!! TEST_TIMEF tests TIMEF.
!
!  Discussion:
!
!    The function
!
!      real ( kind = 8 ) function timef ( ) 
!
!    returns the elapsed time in milliseconds since the previous call to TIMEF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 October 2007
!
!  Author:
!
!    John Burkardt
!
  use xlfutility

  implicit none

  real ( kind = 8 ) value

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_TIMEF'
  write ( *, '(a)' ) '  VALUE = TIMEF ( ) returns the elapsed time'
  write ( *, '(a)' ) '  in milliseconds, since the previous call to TIMEF.'

  value = timef ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Value returned by TIMEF = ', value

  return
end
subroutine erf_values ( n_data, x, fx )

!*****************************************************************************80
!
!! ERF_VALUES returns some values of the ERF or "error" function.
!
!  Discussion:
!
!    The error function is defined by:
!
!      ERF(X) = ( 2 / sqrt ( PI ) * integral ( 0 <= T <= X ) exp ( - T^2 ) dT
!
!    In Mathematica, the function can be evaluated by:
!
!      Erf[x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer, parameter :: n_max = 21

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.0000000000000000D+00, &
    0.1124629160182849D+00, &
    0.2227025892104785D+00, &
    0.3286267594591274D+00, &
    0.4283923550466685D+00, &
    0.5204998778130465D+00, &
    0.6038560908479259D+00, &
    0.6778011938374185D+00, &
    0.7421009647076605D+00, &
    0.7969082124228321D+00, &
    0.8427007929497149D+00, &
    0.8802050695740817D+00, &
    0.9103139782296354D+00, &
    0.9340079449406524D+00, &
    0.9522851197626488D+00, &
    0.9661051464753107D+00, &
    0.9763483833446440D+00, &
    0.9837904585907746D+00, &
    0.9890905016357307D+00, &
    0.9927904292352575D+00, &
    0.9953222650189527D+00 /) 
  integer n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.0D+00, & 
    0.1D+00, & 
    0.2D+00, & 
    0.3D+00, & 
    0.4D+00, & 
    0.5D+00, & 
    0.6D+00, & 
    0.7D+00, & 
    0.8D+00, & 
    0.9D+00, & 
    1.0D+00, & 
    1.1D+00, & 
    1.2D+00, & 
    1.3D+00, & 
    1.4D+00, & 
    1.5D+00, & 
    1.6D+00, & 
    1.7D+00, & 
    1.8D+00, & 
    1.9D+00, & 
    2.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine erfc_values ( n_data, x, fx )

!*****************************************************************************80
!
!! ERFC_VALUES returns some values of the ERFC or "complementary error" function.
!
!  Discussion:
!
!    The complementary error function is defined by:
!
!      ERFC(X) = 1 - ( 2 / sqrt ( PI ) * integral ( 0 <= T <= X ) exp ( - T^2 ) dT
!
!    In Mathematica, the function can be evaluated by:
!
!      Erfc[x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer, parameter :: n_max = 21

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    1.000000000000000D+00, &
    0.7772974107895215D+00, &
    0.5716076449533315D+00, &
    0.3961439091520741D+00, &
    0.2578990352923395D+00, &
    0.1572992070502851D+00, &
    0.08968602177036462D+00, &
    0.04771488023735119D+00, &
    0.02365161665535599D+00, &
    0.01090949836426929D+00, &
    0.004677734981047266D+00, &
    0.001862846297981891D+00, &
    0.0006885138966450786D+00, &
    0.0002360344165293492D+00, &
    0.00007501319466545902D+00, &
    0.00002209049699858544D+00, &
    6.025761151762095D-06, &
    1.521993362862285D-06, &
    3.558629930076853D-07, &
    7.700392745696413D-08, &
    1.541725790028002D-08 /) 
  integer n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.0D+00, & 
    0.2D+00, & 
    0.4D+00, & 
    0.6D+00, & 
    0.8D+00, & 
    1.0D+00, & 
    1.2D+00, & 
    1.4D+00, & 
    1.6D+00, & 
    1.8D+00, & 
    2.0D+00, &
    2.2D+00, & 
    2.4D+00, & 
    2.6D+00, & 
    2.8D+00, & 
    3.0D+00, & 
    3.2D+00, & 
    3.4D+00, & 
    3.6D+00, & 
    3.8D+00, & 
    4.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine gamma_values ( n_data, x, fx )

!*****************************************************************************80
!
!! GAMMA_VALUES returns some values of the Gamma function.
!
!  Discussion:
!
!    The Gamma function is defined as:
!
!      Gamma(Z) = Integral ( 0 <= T < Infinity) T**(Z-1) exp(-T) dT
!
!    It satisfies the recursion:
!
!      Gamma(X+1) = X * Gamma(X)
!
!    Gamma is undefined for nonpositive integral X.
!    Gamma(0.5) = sqrt(PI)
!    For N a positive integer, Gamma(N+1) = N!, the standard factorial.
!
!    In Mathematica, the function can be evaluated by:
!
!      Gamma[x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer, parameter :: n_max = 25

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    -0.3544907701811032D+01, &
    -0.1005871979644108D+03, &
     0.9943258511915060D+02, &
     0.9513507698668732D+01, &
     0.4590843711998803D+01, &
     0.2218159543757688D+01, &
     0.1772453850905516D+01, &
     0.1489192248812817D+01, &
     0.1164229713725303D+01, &
     0.1000000000000000D+01, &
     0.9513507698668732D+00, &
     0.9181687423997606D+00, &
     0.8974706963062772D+00, &
     0.8872638175030753D+00, &
     0.8862269254527580D+00, &
     0.8935153492876903D+00, &
     0.9086387328532904D+00, &
     0.9313837709802427D+00, &
     0.9617658319073874D+00, &
     0.1000000000000000D+01, &
     0.2000000000000000D+01, &
     0.6000000000000000D+01, &
     0.3628800000000000D+06, &
     0.1216451004088320D+18, &
     0.8841761993739702D+31 /)
  integer n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    -0.50D+00, &
    -0.01D+00, &
     0.01D+00, &
     0.10D+00, &
     0.20D+00, &
     0.40D+00, &
     0.50D+00, &
     0.60D+00, &
     0.80D+00, &
     1.00D+00, &
     1.10D+00, &
     1.20D+00, &
     1.30D+00, &
     1.40D+00, &
     1.50D+00, &
     1.60D+00, &
     1.70D+00, &
     1.80D+00, &
     1.90D+00, &
     2.00D+00, &
     3.00D+00, &
     4.00D+00, &
    10.00D+00, &
    20.00D+00, &
    30.00D+00 /) 

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine gamma_log_values ( n_data, x, fx )

!*****************************************************************************80
!
!! GAMMA_LOG_VALUES returns some values of the Log Gamma function.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      Log[Gamma[x]]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer, parameter :: n_max = 20

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
     0.1524063822430784D+01, &
     0.7966778177017837D+00, &
     0.3982338580692348D+00, &
     0.1520596783998375D+00, &
     0.0000000000000000D+00, &
    -0.4987244125983972D-01, &
    -0.8537409000331584D-01, &
    -0.1081748095078604D+00, &
    -0.1196129141723712D+00, &
    -0.1207822376352452D+00, &
    -0.1125917656967557D+00, &
    -0.9580769740706586D-01, &
    -0.7108387291437216D-01, &
    -0.3898427592308333D-01, &
    0.00000000000000000D+00, &
    0.69314718055994530D+00, &
    0.17917594692280550D+01, &
    0.12801827480081469D+02, &
    0.39339884187199494D+02, &
    0.71257038967168009D+02 /) 
  integer n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
     0.20D+00, &
     0.40D+00, &
     0.60D+00, &
     0.80D+00, &
     1.00D+00, &
     1.10D+00, &
     1.20D+00, &
     1.30D+00, &
     1.40D+00, &
     1.50D+00, &
     1.60D+00, &
     1.70D+00, &
     1.80D+00, &
     1.90D+00, &
     2.00D+00, &
     3.00D+00, &
     4.00D+00, &
    10.00D+00, &
    20.00D+00, &
    30.00D+00 /) 

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
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

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
