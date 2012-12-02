program main

!*****************************************************************************80
!
!! MAIN is the main program for GFORTRAN_INTRINSICS.
!
!  Discussion:
!
!    GFORTRAN_INTRINSICS calls some of the GFORTRAN intrinsic routines.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 February 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GFORTRAN_INTRINSICS'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the GFORTRAN intrinsic library.'

  call test_access ( )
  call test_algama ( )
  call test_besj0 ( )
  call test_besj1 ( )
  call test_besjn ( )
  call test_besy0 ( )
  call test_besy1 ( )
  call test_besyn ( )
  call test_dbesj0 ( )
  call test_dbesj1 ( )
  call test_dbesjn ( )
  call test_dbesy0 ( )
  call test_dbesy1 ( )
  call test_dbesyn ( )
  call test_derf ( )
  call test_derfc ( )
  call test_dgamma ( )
  call test_dlgama ( )
  call test_dtime ( )
  call test_erf ( )
  call test_erfc ( )
  call test_etime ( )
  call test_fdate ( )
  call test_gamma ( )
  call test_get_command ( )
  call test_getcwd ( )
  call test_get_environment_variable ( )
  call test_hostnm ( )
  call test_isnan ( )
  call test_rand ( )
  call test_secnds ( )
  call test_sizeof ( )
  call test_sleep ( )
  call test_srand ( )
  call test_system ( )
  call test_time ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GFORTRAN_INTRINSICS'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )
!
!  Demonstrate one more intrinsic, the ABORT call.
!
  call test_abort

  stop
end
subroutine test_abort ( )

!*****************************************************************************80
!
!! TEST_ABORT demonstrates the ABORT routine.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 February 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_ABORT:'
  write ( *, '(a)' ) '  The ABORT routine causes program termination'
  write ( *, '(a)' ) '  possibly with a core dump.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  It does set the error status flag that can be'
  write ( *, '(a)' ) '  detected by checking the $STATUS environment variable.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Presumably, this has to be the LAST demonstration!'

  call abort

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  If this message is printed, then ABORT did not'
  write ( *, '(a)' ) '  do its job.'

  return
end
subroutine test_access ( )

!*****************************************************************************80
!
!! TEST_ACCESS demonstrates the ACCESS routine.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 February 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 27 ) name
  character ( len = 3  ) mode
  integer   ( kind = 4 ) result

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_ACCESS:'
  write ( *, '(a)' ) '  The ACCESS routine checks whether a file exists,'
  write ( *, '(a)' ) '  and whether it may be accessed with given modes.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '        Filename               Mode    Result'
  write ( *, '(a)' ) '  ---------------------------  ----    ------'
  write ( *, '(a)' ) ' '

  name = 'gastropods_are_uneaten.txt'
  mode = ' '
  result = access ( name, mode )
  write ( *, '(2x,a27,2x,a5,2x,i8)' ) name, '''' // mode // '''', result

  name = 'gfortran_intrinsics_prb.f90'
  mode = ' '
  result = access ( name, mode )
  write ( *, '(2x,a27,2x,a5,2x,i8)' ) name, '''' // mode // '''', result

  name = 'gfortran_intrinsics_prb.f90'
  mode = 'x'
  result = access ( name, mode )
  write ( *, '(2x,a27,2x,a5,2x,i8)' ) name, '''' // mode // '''', result

  name = 'gfortran_intrinsics_prb.f90'
  mode = 'r'
  result = access ( name, mode )
  write ( *, '(2x,a27,2x,a5,2x,i8)' ) name, '''' // mode // '''', result

  name = 'gfortran_intrinsics_prb.f90'
  mode = 'w'
  result = access ( name, mode )
  write ( *, '(2x,a27,2x,a5,2x,i8)' ) name, '''' // mode // '''', result

  name = 'gfortran_intrinsics_prb.f90'
  mode = 'rwx'
  result = access ( name, mode )

  write ( *, '(2x,a27,2x,a5,2x,i8)' ) name, '''' // mode // '''', result

  return
end
subroutine test_algama ( )

!*****************************************************************************80
!
!! TEST_ALGAMA checks ALGAMA against GAMMA_LOG_VALUES.
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
  real ( kind = 4 ) fx2
  integer n_data
  real ( kind = 8 ) x
  real ( kind = 4 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_ALGAMA:'
  write ( *, '(a)' ) &
    '  ALGAMA computes the log of the gamma function.'
  write ( *, '(a)' ) '  GAMMA_LOG_VALUES returns selected values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) &
    '      X             FX', &
    '                        FX2'
  write ( *, '(a,a)' ) &
    '                   (table)', &
    '                   (ALGAMA)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call gamma_log_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    x2 = real ( x, kind = 4 )
    fx2 = algama ( x2 )

    write ( *, '(2x,f12.8,2x,g24.16,2x,g24.16)' ) &
      x, fx, fx2

  end do

  return
end
subroutine test_besj0 ( )

!*****************************************************************************80
!
!! TEST_BESJ0 checks BESJ0 against BESSEL_J0_VALUES.
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
  write ( *, '(a)' ) 'TEST_BESJ0:'
  write ( *, '(a)' ) '  BESJ0 computes the Bessel J0 function.'
  write ( *, '(a)' ) '  BESSEL_J0_VALUES returns selected values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) &
    '      X             FX', &
    '                        FX2'
  write ( *, '(a,a)' ) &
    '                   (table)', &
    '                   (BESJ0)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_j0_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    x2 = real ( x )
    fx2 = besj0 ( x2 )

    write ( *, '(2x,f12.8,2x,g24.16,2x,g24.16)' ) &
      x, fx, fx2

  end do

  return
end
subroutine test_besj1 ( )

!*****************************************************************************80
!
!! TEST_BESJ1 checks BESJ1 against BESSEL_J1_VALUES.
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
  write ( *, '(a)' ) 'TEST_BESJ1:'
  write ( *, '(a)' ) '  BESJ1 computes the Bessel J1 function.'
  write ( *, '(a)' ) '  BESSEL_J1_VALUES returns selected values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) &
    '      X             FX', &
    '                        FX2'
  write ( *, '(a,a)' ) &
    '                   (table)', &
    '                   (BESJ1)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_j1_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    x2 = real ( x )
    fx2 = besj1 ( x2 )

    write ( *, '(2x,f12.8,2x,g24.16,2x,g24.16)' ) &
      x, fx, fx2

  end do

  return
end
subroutine test_besjn ( )

!*****************************************************************************80
!
!! TEST_BESJN checks BESJN against BESSEL_JN_VALUES.
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
  integer n
  integer n_data
  real ( kind = 8 ) x
  real x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_BESJN:'
  write ( *, '(a)' ) '  BESJN computes the Bessel Jn function.'
  write ( *, '(a)' ) '  BESSEL_JN_VALUES returns selected values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) &
    '      N     X             FX', &
    '                        FX2'
  write ( *, '(a,a)' ) &
    '                         (table)', &
    '                   (BESJN)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_jn_values ( n_data, n, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    x2 = real ( x )
    fx2 = besjn ( n, x2 )

    write ( *, '(2x,i4,2x,f12.8,2x,g24.16,2x,g24.16)' ) &
      n, x, fx, fx2

  end do

  return
end
subroutine test_besy0 ( )

!*****************************************************************************80
!
!! TEST_BESY0 checks BESY0 against BESSEL_Y0_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 April 2007
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
  write ( *, '(a)' ) 'TEST_BESY0:'
  write ( *, '(a)' ) '  BESY0 computes the Bessel Y0 function.'
  write ( *, '(a)' ) '  BESSEL_Y0_VALUES returns selected values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) &
    '      X             FX', &
    '                        FX2'
  write ( *, '(a,a)' ) &
    '                   (table)', &
    '                   (BESY0)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_y0_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    x2 = real ( x )
    fx2 = besy0 ( x2 )

    write ( *, '(2x,f12.8,2x,g24.16,2x,g24.16)' ) &
      x, fx, fx2

  end do

  return
end
subroutine test_besy1 ( )

!*****************************************************************************80
!
!! TEST_BESY1 checks BESY1 against BESSEL_Y1_VALUES.
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
  write ( *, '(a)' ) 'TEST_BESY1:'
  write ( *, '(a)' ) '  BESY1 computes the Bessel Y1 function.'
  write ( *, '(a)' ) '  BESSEL_Y1_VALUES returns selected values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) &
    '      X             FX', &
    '                        FX2'
  write ( *, '(a,a)' ) &
    '                   (table)', &
    '                   (BESY1)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_y1_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    x2 = real ( x )
    fx2 = besy1 ( x2 )

    write ( *, '(2x,f12.8,2x,g24.16,2x,g24.16)' ) &
      x, fx, fx2

  end do

  return
end
subroutine test_besyn ( )

!*****************************************************************************80
!
!! TEST_BESYN checks BESYN against BESSEL_YN_VALUES.
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
  integer n
  integer n_data
  real ( kind = 8 ) x
  real x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_BESYN:'
  write ( *, '(a)' ) '  BESYN computes the Bessel Yn function.'
  write ( *, '(a)' ) '  BESSEL_YN_VALUES returns selected values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) &
    '      N     X             FX', &
    '                        FX2'
  write ( *, '(a,a)' ) &
    '                         (table)', &
    '                   (BESYN)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_yn_values ( n_data, n, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    x2 = real ( x )
    fx2 = besyn ( n, x2 )

    write ( *, '(2x,i4,2x,f12.8,2x,g24.16,2x,g24.16)' ) &
      n, x, fx, fx2

  end do

  return
end
subroutine test_dbesj0 ( )

!*****************************************************************************80
!
!! TEST_DBESJ0 checks DBESJ0 against BESSEL_J0_VALUES.
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
  real ( kind = 8 ) fx2
  integer n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_DBESJ0:'
  write ( *, '(a)' ) '  DBESJ0 computes the Bessel J0 function.'
  write ( *, '(a)' ) '  BESSEL_J0_VALUES returns selected values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) &
    '      X             FX', &
    '                        FX2'
  write ( *, '(a,a)' ) &
    '                   (table)', &
    '                   (DBESJ0)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_j0_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = dbesj0 ( x )

    write ( *, '(2x,f12.8,2x,g24.16,2x,g24.16)' ) &
      x, fx, fx2

  end do

  return
end
subroutine test_dbesj1 ( )

!*****************************************************************************80
!
!! TEST_DBESJ1 checks DBESJ1 against BESSEL_J1_VALUES.
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
  real ( kind = 8 ) fx2
  integer n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_DBESJ1:'
  write ( *, '(a)' ) '  DBESJ1 computes the Bessel J1 function.'
  write ( *, '(a)' ) '  BESSEL_J1_VALUES returns selected values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) &
    '      X             FX', &
    '                        FX2'
  write ( *, '(a,a)' ) &
    '                   (table)', &
    '                   (DBESJ1)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_j1_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = dbesj1 ( x )

    write ( *, '(2x,f12.8,2x,g24.16,2x,g24.16)' ) &
      x, fx, fx2

  end do

  return
end
subroutine test_dbesjn ( )

!*****************************************************************************80
!
!! TEST_DBESJN checks DBESJN against BESSEL_JN_VALUES.
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
  real ( kind = 8 ) fx2
  integer n
  integer n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_DBESJN:'
  write ( *, '(a)' ) '  DBESJN computes the Bessel Jn function.'
  write ( *, '(a)' ) '  BESSEL_JN_VALUES returns selected values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) &
    '      N     X             FX', &
    '                        FX2'
  write ( *, '(a,a)' ) &
    '                         (table)', &
    '                   (DBESJN)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_jn_values ( n_data, n, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = dbesjn ( n, x )

    write ( *, '(2x,i4,2x,f12.8,2x,g24.16,2x,g24.16)' ) &
      n, x, fx, fx2

  end do

  return
end
subroutine test_dbesy0 ( )

!*****************************************************************************80
!
!! TEST_DBESY0 checks DBESY0 against BESSEL_Y0_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  real ( kind = 8 ) fx2
  integer n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_DBESY0:'
  write ( *, '(a)' ) '  DBESY0 computes the Bessel Y0 function.'
  write ( *, '(a)' ) '  BESSEL_Y0_VALUES returns selected values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) &
    '      X             FX', &
    '                        FX2'
  write ( *, '(a,a)' ) &
    '                   (table)', &
    '                   (DBESY0)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_y0_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = dbesy0 ( x )

    write ( *, '(2x,f12.8,2x,g24.16,2x,g24.16)' ) &
      x, fx, fx2

  end do

  return
end
subroutine test_dbesy1 ( )

!*****************************************************************************80
!
!! TEST_DBESY1 checks DBESY1 against BESSEL_Y1_VALUES.
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
  real ( kind = 8 ) fx2
  integer n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_DBESY1:'
  write ( *, '(a)' ) '  DBESY1 computes the Bessel Y1 function.'
  write ( *, '(a)' ) '  BESSEL_Y1_VALUES returns selected values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) &
    '      X             FX', &
    '                        FX2'
  write ( *, '(a,a)' ) &
    '                   (table)', &
    '                   (DBESY1)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_y1_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = dbesy1 ( x )

    write ( *, '(2x,f12.8,2x,g24.16,2x,g24.16)' ) &
      x, fx, fx2

  end do

  return
end
subroutine test_dbesyn ( )

!*****************************************************************************80
!
!! TEST_DBESYN checks DBESYN against BESSEL_YN_VALUES.
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
  real ( kind = 8 ) fx2
  integer n
  integer n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_DBESYN:'
  write ( *, '(a)' ) '  DBESYN computes the Bessel Yn function.'
  write ( *, '(a)' ) '  BESSEL_YN_VALUES returns selected values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) &
    '      N     X             FX', &
    '                        FX2'
  write ( *, '(a,a)' ) &
    '                         (table)', &
    '                   (DBESYN)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_yn_values ( n_data, n, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = dbesyn ( n, x )

    write ( *, '(2x,i4,2x,f12.8,2x,g24.16,2x,g24.16)' ) &
      n, x, fx, fx2

  end do

  return
end
subroutine test_derf ( )

!*****************************************************************************80
!
!! TEST_DERF checks DERF against ERF_VALUES.
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
  real ( kind = 8 ) fx2
  integer n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_DERF:'
  write ( *, '(a)' ) '  DERF computes the error function.'
  write ( *, '(a)' ) '  ERF_VALUES returns selected values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) &
    '      X             FX', &
    '                        FX2'
  write ( *, '(a,a)' ) &
    '                   (table)', &
    '                   (DERF)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call erf_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = derf ( x )

    write ( *, '(2x,f12.8,2x,g24.16,2x,g24.16)' ) &
      x, fx, fx2

  end do

  return
end
subroutine test_derfc ( )

!*****************************************************************************80
!
!! TEST_DERFC checks DERFC against ERFC_VALUES.
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
  real ( kind = 8 ) fx2
  integer n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_DERFC:'
  write ( *, '(a)' ) '  DERFC computes the complementary error function.'
  write ( *, '(a)' ) '  ERFC_VALUES returns selected values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) &
    '      X             FX', &
    '                        FX2'
  write ( *, '(a,a)' ) &
    '                   (table)', &
    '                   (DERFC)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call erfc_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = derfc ( x )

    write ( *, '(2x,f12.8,2x,g24.16,2x,g24.16)' ) &
      x, fx, fx2

  end do

  return
end
subroutine test_dgamma ( )

!*****************************************************************************80
!
!! TEST_DGAMMA checks DGAMMA against GAMMA_VALUES.
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
  real ( kind = 8 ) fx2
  integer n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_DGAMMA:'
  write ( *, '(a)' ) &
    '  DGAMMA computes the gamma function.'
  write ( *, '(a)' ) '  GAMMA_VALUES returns selected values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) &
    '      X             FX', &
    '                        FX2'
  write ( *, '(a,a)' ) &
    '                   (table)', &
    '                   (DGAMMA)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call gamma_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = dgamma ( x )

    write ( *, '(2x,f12.8,2x,g24.16,2x,g24.16)' ) &
      x, fx, fx2

  end do

  return
end
subroutine test_dlgama ( )

!*****************************************************************************80
!
!! TEST_DLGAMA checks DLGAMA against GAMMA_LOG_VALUES.
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
  real ( kind = 8 ) fx2
  integer n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_DLGAMA:'
  write ( *, '(a)' ) &
    '  DLGAMA computes the log of the gamma function.'
  write ( *, '(a)' ) '  GAMMA_LOG_VALUES returns selected values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) &
    '      X             FX', &
    '                        FX2'
  write ( *, '(a,a)' ) &
    '                   (table)', &
    '                   (DLGAMA)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call gamma_log_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = dlgama ( x )

    write ( *, '(2x,f12.8,2x,g24.16,2x,g24.16)' ) &
      x, fx, fx2

  end do

  return
end
subroutine test_dtime ( )

!*****************************************************************************80
!
!! TEST_DTIME tests DTIME.
!
!  Discussion:
!
!    call dtime ( tarray, result ) 
!
!    returns in the array TARRAY(2) the user and system execution times
!    and in RESULT the run time that has elapsed since the start of the
!    program, or the previous call to DTIME.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 February 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 4 ) result
  real ( kind = 4 ) tarray(2)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_DTIME'
  write ( *, '(a)' ) '  DTIME returns the user and system execution time'
  write ( *, '(a)' ) '  since the program beginning, or the previous'
  write ( *, '(a)' ) '  call to DTIME.'

  write ( *, '(a)' ) '  Well, that''s what the manual says.  But in fact,'
  write ( *, '(a)' ) '  DTIME seems to be identical in behavior to ETIME.'
  write ( *, '(a)' ) '  We always seem to get total elapsed time,'
  write ( *, '(a)' ) '  not the "delta" time or increment in time.'

  call dtime ( tarray, result )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  TARRAY(1) (User time)   = ', tarray(1)
  write ( *, '(a,g14.6)' ) '  TARRAY(2) (System time) = ', tarray(2)
  write ( *, '(a,g14.6)' ) '  RESULT    (Run time)    = ', result

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now sleep 2 seconds and call DTIME again.'
  write ( *, '(a)' ) '  We might expect to get ZERO, at least for CPU time.'

  call sleep ( 2 )

  call dtime ( tarray, result )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  TARRAY(1) (User time)   = ', tarray(1)
  write ( *, '(a,g14.6)' ) '  TARRAY(2) (System time) = ', tarray(2)
  write ( *, '(a,g14.6)' ) '  RESULT    (Run time)    = ', result

  return
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
subroutine test_etime ( )

!*****************************************************************************80
!
!! TEST_ETIME tests ETIME.
!
!  Discussion:
!
!    call etime ( tarray, result ) 
!
!    returns in the real array TARRAY(2) the user and system execution times
!    and in RESULT the run time since starting.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 February 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real result
  real tarray(2)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_ETIME'
  write ( *, '(a)' ) '  ETIME returns the user and system execution time.'

  call etime ( tarray, result )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  TARRAY(1) (User time)   = ', tarray(1)
  write ( *, '(a,g14.6)' ) '  TARRAY(2) (System time) = ', tarray(2)
  write ( *, '(a,g14.6)' ) '  RESULT    (Run time)    = ', result

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now sleep 2 seconds.'

  call sleep ( 2 )

  call etime ( tarray, result )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  TARRAY(1) (User time)   = ', tarray(1)
  write ( *, '(a,g14.6)' ) '  TARRAY(2) (System time) = ', tarray(2)
  write ( *, '(a,g14.6)' ) '  RESULT    (Run time)    = ', result

  return
end
subroutine test_fdate ( )

!*****************************************************************************80
!
!! TEST_FDATE tests FDATE.
!
!  Discussion:
!
!    call fdate ( string ) 
!    or
!    string = fdate ( )
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

  call fdate ( string )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_FDATE'
  write ( *, '(a)' ) '  FDATE returns the current date as a string.'
  write ( *, '(a)' ) '  It can be called as a function or subroutine.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  CALL FDATE ( STRING ) returns:'
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
subroutine test_hostnm ( )

!*****************************************************************************80
!
!! TEST_HOSTNM tests HOSTNM.
!
!  Discussion:
!
!    result = hostnm ( string )
!
!    returns a success code in the integer RESULT, and sets STRING
!    to a character string of the host name.
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

  result = hostnm ( string )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_HOSTNM'
  write ( *, '(a)' ) '  HOSTNM returns the host name as a string.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  RESULT = HOSTNM ( STRING ) returns:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  STRING = "' // trim ( string ) // '".'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  RESULT (0 = success) = ', result

  return
end
subroutine test_get_command ( )

!*****************************************************************************80
!
!! TEST_GET_COMMAND tests GET_COMMAND.
!
!  Discussion:
!
!    call get_command ( string )
!
!    returns the comand line used to invoke the program.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 February 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 80 ) string

  call get_command ( string )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_GET_COMMAND'
  write ( *, '(a)' ) '  GET_COMMAND returns the command line used to'
  write ( *, '(a)' ) '  invoke the current program.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Command = "' // trim ( string ) // '".'

  return
end
subroutine test_getcwd ( )

!*****************************************************************************80
!
!! TEST_GETCWD tests GETCWD.
!
!  Discussion:
!
!    call getcwd ( string )
!
!    returns the current working directory.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 February 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 80 ) string

  call getcwd ( string )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_GETCWD'
  write ( *, '(a)' ) '  GETCWD returns the current working directory.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  CWD = "' // trim ( string ) // '".'

  return
end
subroutine test_get_environment_variable ( )

!*****************************************************************************80
!
!! TEST_GET_ENVIRONMENT_VARIABLE tests GET_ENVIRONMENT_VARIABLE.
!
!  Discussion:
!
!    call get_environment_variable ( string )
!
!    returns an environmental variable.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 February 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 80 ) envvar
  character ( len = 80 ) string

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_GET_ENVIRONMENT_VARIABLE'
  write ( *, '(a)' ) '  GET_ENVIRONMENT_VARIABLE returns an environmental variable.'
  write ( *, '(a)' ) '  Note that this can be case-sensitive!'

  envvar = 'home'
  call get_environment_variable ( envvar, string )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Variable = "' // trim ( envvar ) // '".'
  write ( *, '(a)' ) '  Value =    "' // trim ( string ) // '".'

  envvar = 'HOME'
  call get_environment_variable ( envvar, string )

  return
end
subroutine test_isnan ( )

!*****************************************************************************80
!
!! TEST_ISNAN tests ISNAN.
!
!  Discussion:
!
!    ISNAN ( X ) is a logical function which is TRUE if the number X
!    is a "NaN", or "Not a Number".  This is a special value set to
!    indicate that an illegal argument was input to a function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 September 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  logical l(4)
  real v(4)
  real x(4)

  x(1) = - 1.0E+00
  x(2) =   0.0E+00
  x(3) = + 1.0E+00
  x(4) = + 2.0E+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_ISNAN'
  write ( *, '(a)' ) '  ISNAN(X) is TRUE if X is "Not a Number".'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   Function  -1  0 +1 +2'
  write ( *, '(a)' ) ' '
  
  v(1:4) = acos ( x(1:4) )
  l(1:4) = isnan ( v(1:4) )
  write ( *, '(2x,a10,4(2x,l1))' ) '   ACOS(X)', l(1:4)

  v(1:4) = asin ( x(1:4) )
  l(1:4) = isnan ( v(1:4) )
  write ( *, '(2x,a10,4(2x,l1))' ) '   ASIN(X)', l(1:4)

  v(1:4) = atan ( x(1:4) )
  l(1:4) = isnan ( v(1:4) )
  write ( *, '(2x,a10,4(2x,l1))' ) '   ATAN(X)', l(1:4)

  v(1:4) = log ( x(1:4) )
  l(1:4) = isnan ( v(1:4) )
  write ( *, '(2x,a10,4(2x,l1))' ) '    LOG(X)', l(1:4)

  v(1:4) = sqrt ( x(1:4) )
  l(1:4) = isnan ( v(1:4) )
  write ( *, '(2x,a10,4(2x,l1))' ) '   SQRT(X)', l(1:4)

  v(1:4) = 1.0E+00 / ( x(1:4) )
  l(1:4) = isnan ( v(1:4) )
  write ( *, '(2x,a10,4(2x,l1))' ) '     1 / X', l(1:4)

  return
end
subroutine test_rand ( )

!*****************************************************************************80
!
!! TEST_RAND tests RAND.
!
!  Discussion:
!
!    The seed is optional.
!
!    If it is 0, it is ignored and the next value in the random sequence
!    is returned.
!
!    If it is 1, the generator is restarted by calling SRAND(0).
!
!    Otherwise, it is used as input to SRAND to reseed the sequence.
!
!    Normal usage of RAND would be to call once with the seed nonzero,
!    and then to call repeatedly with the seed 0 or not specified.
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

  integer i
  real r
  integer seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_RAND'
  write ( *, '(a)' ) '  RAND returns a real random value.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The sequence is "seeded" by:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    R = RAND ( 1 ): reseeds by SRAND ( 0 )'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    or'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    R = RAND ( S ): reseeds by SRAND ( S )'
  write ( *, '(a)' ) '                    assuming S not 0, not 1.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The sequence is used by:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    R = RAND ( ): returns next value in current'
  write ( *, '(a)' ) '                  random sequence.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    or'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    R = RAND ( 0 ): same as R = RAND ( );'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    or'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    R = RAND ( S ): reseeds sequence, and returns'
  write ( *, '(a)' ) '                    first value in new sequence,'
  write ( *, '(a)' ) '                    assuming S is not 0.'

  seed = 123456789
  call srand ( seed )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Call SRAND(123456789)'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  do i = 1, 10'
  write ( *, '(a)' ) '    R = RAND ( )'
  write ( *, '(a)' ) '  end do'
  write ( *, '(a)' ) ' '

  do i = 1, 10
    r = rand ( )
    write ( *, '(2x,i8,2x,g14.6)' ) i, r
  end do

  seed = 987654321
  call srand ( seed )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  CHANGING THE SEED CHANGES THE SEQUENCE:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Call SRAND(987654321)'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  do i = 1, 10'
  write ( *, '(a)' ) '    R = RAND ( )'
  write ( *, '(a)' ) '  end do'
  write ( *, '(a)' ) ' '

  do i = 1, 10
    r = rand ( )
    write ( *, '(2x,i8,2x,g14.6)' ) i, r
  end do

  seed = 123456789
  call srand ( seed )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  RESTORING THE OLD SEED RESTARTS THE OLD SEQUENCE:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Call SRAND(123456789)'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  do i = 1, 10'
  write ( *, '(a)' ) '    R = RAND ( )'
  write ( *, '(a)' ) '  end do'
  write ( *, '(a)' ) ' '

  do i = 1, 10
    r = rand ( )
    write ( *, '(2x,i8,2x,g14.6)' ) i, r
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  WE CAN PASS THE SEED IN THROUGH RAND:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  R = RAND(987654321)'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  do i = 2, 10'
  write ( *, '(a)' ) '    R = RAND ( )'
  write ( *, '(a)' ) '  end do'
  write ( *, '(a)' ) ' '

  r = rand ( 987654321 )

  i = 1
  write ( *, '(2x,i8,2x,g14.6)' ) i, r

  do i = 2, 10
    r = rand ( )
    write ( *, '(2x,i8,2x,g14.6)' ) i, r
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  WE CAN GET A "RANDOM" SEED:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  R = RAND(1)'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  do i = 2, 10'
  write ( *, '(a)' ) '    R = RAND ( )'
  write ( *, '(a)' ) '  end do'
  write ( *, '(a)' ) ' '

  r = rand ( 1 )

  i = 1
  write ( *, '(2x,i8,2x,g14.6)' ) i, r

  do i = 2, 10
    r = rand ( )
    write ( *, '(2x,i8,2x,g14.6)' ) i, r
  end do

  seed = 123456789
  call srand ( seed )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  CALLING WITH RAND(0) IS THE SAME AS RAND():'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Call SRAND(123456789)'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  do i = 1, 10'
  write ( *, '(a)' ) '    R = RAND ( )'
  write ( *, '(a)' ) '  end do'
  write ( *, '(a)' ) ' '

  do i = 1, 10
    r = rand ( 0 )
    write ( *, '(2x,i8,2x,g14.6)' ) i, r
  end do

  return
end
subroutine test_secnds ( )

!*****************************************************************************80
!
!! TEST_SECNDS tests SECNDS.
!
!  Discussion:
!
!    The GFORTRAN function
!
!      integer function secnds ( t ) 
!
!    returns the current local time in seconds, since midnight, minus the
!    value of t.
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

  integer s
  real t

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_SECNDS'
  write ( *, '(a)' ) '  I = SECNDS ( T ) returns the local time, in seconds'
  write ( *, '(a)' ) '  since midnight, minus T.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Note that T is real, and I is an integer.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The code fragment:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    T1 = SECNDS ( 0.0 )'
  write ( *, '(a)' ) '    stuff happens'
  write ( *, '(a)' ) '    T2 = SECNDS ( real ( T1 ) )'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  will compute the (wallclock) time elapsed while'
  write ( *, '(a)' ) '  stuff happens.'

  t = 0.0
  s = secnds ( t )
  write ( *, '(a)' ) ' '
  write ( *, '(a,f12.4,a,i8)' ) '  Value returned by SECNDS ( ', t, ' ) = ', s
  t = real ( s )
  s = secnds ( t )
  write ( *, '(a,f12.4,a,i8)' ) '  Value returned by SECNDS ( ', t, ' ) = ', s

  return
end
subroutine test_sleep ( )

!*****************************************************************************80
!
!! TEST_SLEEP tests SLEEP.
!
!  Discussion:
!
!    The GFORTRAN routine
!
!      subroutine sleep ( t ) 
!
!    causes the program to "sleep" for T seconds.
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

  integer value

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_SLEEP'
  write ( *, '(a)' ) '  CALL SLEEP ( T ) causes the program to "sleep"'
  write ( *, '(a)' ) '  for T seconds.'

  value = time ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  Current TIME in seconds is = ', value

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now CALL SLEEP ( 10 )'

  value = 10
  call sleep ( value )
  
  value = time ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  Current TIME in seconds is = ', value

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
subroutine test_srand ( )

!*****************************************************************************80
!
!! TEST_SRAND tests SRAND.
!
!  Discussion:
!
!    RAND is the random number generator included in the GFORTRAN library.
!
!    The values computed by RAND are based on an internal "seed" value.
!
!    RAND has an argument X, whose purpose is as follows:
!
!    * RAND ( 0 ) means update the seed internally and return the next 
!      random value.
!    * RAND ( 1 ) means that the seed is reset by a call to SRAND ( 0 ),
!      following which the next random value is computed.
!    * RAND ( N ) where N is not 0 or 1 means the seed is reset by a
!      call to SRAND ( N ), following which the next random value is
!      computed.
!
!    The routine SRAND allows the user to modify the seed directly.
!    Note that "0" is a bad idea for a seed!
!
!    Typically, SRAND is used in order to specify how a given sequence
!    of random numbers was generated by RAND.  Knowing the seed value
!    that was used, the same stream of random values can be regenerated.
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
  integer seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_SRAND'
  write ( *, '(a)' ) '  SRAND resets the seed used by RAND,'
  write ( *, '(a)' ) '  the random number generator.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  If you do nothing before calling RAND, a seed'
  write ( *, '(a)' ) '  is chosen for you.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  To control the sequence, especially so you can'
  write ( *, '(a)' ) '  regenerate it later, use'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    call SRAND ( SEED )'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  using a NONZERO value of SEED.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  If you use a ZERO value for SEED, the sequence is'
  write ( *, '(a)' ) '  reset, but with a SEED value chosen at random.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  To start, we just call RAND without calling'
  write ( *, '(a)' ) '  SRAND first.  Presumably, we get a "random" seed.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  do i = 1, 10'
  write ( *, '(a)' ) '    R = RAND ( )'
  write ( *, '(a)' ) '  end do'
  write ( *, '(a)' ) ' '

  do i = 1, 10
    r = rand ( )
    write ( *, '(2x,i8,2x,g14.6)' ) i, r
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  To generate a repeatable sequence, we must'
  write ( *, '(a)' ) '  call SRAND.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SEED = 123456789 '
  write ( *, '(a)' ) '  call SRAND ( SEED )'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  do i = 1, 10'
  write ( *, '(a)' ) '    r = RAND ( )'
  write ( *, '(a)' ) '  end do'
  write ( *, '(a)' ) ' '

  seed = 123456789
  call srand ( seed )

  do i = 1, 10
    r = rand ( )
    write ( *, '(2x,i8,2x,g14.6)' ) i, r
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Changing the seed changes the sequence:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SEED = 987654321 '
  write ( *, '(a)' ) '  call SRAND ( SEED )'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  do i = 1, 10'
  write ( *, '(a)' ) '    r = RAND ( )'
  write ( *, '(a)' ) '  end do'
  write ( *, '(a)' ) ' '

  seed = 987654321
  call srand ( seed )

  do i = 1, 10
    r = rand ( )
    write ( *, '(2x,i8,2x,g14.6)' ) i, r
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Restoring the seed restarts the sequence:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SEED = 123456789 '
  write ( *, '(a)' ) '  call SRAND ( SEED )'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  do i = 1, 10'
  write ( *, '(a)' ) '    r = RAND ( )'
  write ( *, '(a)' ) '  end do'
  write ( *, '(a)' ) ' '

  seed = 123456789
  call srand ( seed )

  do i = 1, 10
    r = rand ( )
    write ( *, '(2x,i8,2x,g14.6)' ) i, r
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We can implicitly call SRAND by calling RAND'
  write ( *, '(a)' ) '  with the new seed value in the argument.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SEED = 987654321'
  write ( *, '(a)' ) '  r = RAND ( SEED )'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  do i = 2, 10'
  write ( *, '(a)' ) '    r = RAND ( )'
  write ( *, '(a)' ) '  end do'
  write ( *, '(a)' ) ' '

  seed = 987654321
  r = rand ( seed )

  i = 1
  write ( *, '(2x,i8,2x,g14.6)' ) i, r

  do i = 2, 10
    r = rand ( )
    write ( *, '(2x,i8,2x,g14.6)' ) i, r
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We can get a random seed by using an argument of 0:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SEED = 0'
  write ( *, '(a)' ) '  call SRAND ( SEED )'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  do i = 1, 10'
  write ( *, '(a)' ) '    r = RAND ( )'
  write ( *, '(a)' ) '  end do'
  write ( *, '(a)' ) ' '

  seed = 0
  call srand ( seed )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The output value of SEED is NOT updated.'
  write ( *, '(a,i12)' ) '  SEED on return = ', seed
  write ( *, '(a)' ) ' '

  do i = 1, 10
    r = rand ( )
    write ( *, '(2x,i8,2x,g14.6)' ) i, r
  end do

  return
end
subroutine test_system ( )

!*****************************************************************************80
!
!! TEST_SYSTEM tests SYSTEM.
!
!  Discussion:
!
!    The GFORTRAN function
!
!      result = system ( string ) 
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
  write ( *, '(a)' ) '    RESULT = SYSTEM ( ''date'' )'
  write ( *, '(a)' ) ' '

  string = 'date'
  result = system ( string )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  RESULT ( 0 = success ) = ', result

  return
end
subroutine test_time ( )

!*****************************************************************************80
!
!! TEST_TIME tests TIME.
!
!  Discussion:
!
!    The GFORTRAN function
!
!      integer function time ( ) 
!
!    returns the current time as an integer, similar to the workings of 
!    the UNIX "time" function.  This means the integer counts the number 
!    of seconds since 01 January 1970.
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

  integer value

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_TIME'
  write ( *, '(a)' ) '  I = TIME ( ) returns the time, encoded as the number'
  write ( *, '(a)' ) '  of seconds since 01 January 1970.'

  value = time ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  Value returned by TIME = ', value

  return
end
subroutine bessel_j0_values ( n_data, x, fx )

!*****************************************************************************80
!
!! BESSEL_J0_VALUES returns some values of the J0 Bessel function.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      BesselJ[0,x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 August 2004
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
    -0.1775967713143383D+00, &
    -0.3971498098638474D+00, &
    -0.2600519549019334D+00, &
     0.2238907791412357D+00, &
     0.7651976865579666D+00, &
     0.1000000000000000D+01, &
     0.7651976865579666D+00, &
     0.2238907791412357D+00, &
    -0.2600519549019334D+00, &
    -0.3971498098638474D+00, &
    -0.1775967713143383D+00, &
     0.1506452572509969D+00, &
     0.3000792705195556D+00, &
     0.1716508071375539D+00, &
    -0.9033361118287613D-01, &
    -0.2459357644513483D+00, &
    -0.1711903004071961D+00, &
     0.4768931079683354D-01, &
     0.2069261023770678D+00, &
     0.1710734761104587D+00, &
    -0.1422447282678077D-01 /)
  integer n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    -5.0D+00, &
    -4.0D+00, &
    -3.0D+00, &
    -2.0D+00, &
    -1.0D+00, &
     0.0D+00, &
     1.0D+00, &
     2.0D+00, &
     3.0D+00, &
     4.0D+00, &
     5.0D+00, &
     6.0D+00, &
     7.0D+00, &
     8.0D+00, &
     9.0D+00, &
    10.0D+00, &
    11.0D+00, &
    12.0D+00, &
    13.0D+00, &
    14.0D+00, &
    15.0D+00 /)

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
subroutine bessel_j1_values ( n_data, x, fx )

!*****************************************************************************80
!
!! BESSEL_J1_VALUES returns some values of the J1 Bessel function.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      BesselJ[1,x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 August 2004
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
     0.3275791375914652D+00, & 
     0.6604332802354914D-01, & 
    -0.3390589585259365D+00, & 
    -0.5767248077568734D+00, & 
    -0.4400505857449335D+00, & 
     0.0000000000000000D+00, & 
     0.4400505857449335D+00, & 
     0.5767248077568734D+00, & 
     0.3390589585259365D+00, & 
    -0.6604332802354914D-01, & 
    -0.3275791375914652D+00, & 
    -0.2766838581275656D+00, & 
    -0.4682823482345833D-02, & 
     0.2346363468539146D+00, & 
     0.2453117865733253D+00, & 
     0.4347274616886144D-01, & 
    -0.1767852989567215D+00, & 
    -0.2234471044906276D+00, & 
    -0.7031805212177837D-01, & 
     0.1333751546987933D+00, & 
     0.2051040386135228D+00 /)
  integer n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    -5.0D+00, &
    -4.0D+00, &
    -3.0D+00, &
    -2.0D+00, &
    -1.0D+00, &
     0.0D+00, &
     1.0D+00, &
     2.0D+00, &
     3.0D+00, &
     4.0D+00, &
     5.0D+00, &
     6.0D+00, &
     7.0D+00, &
     8.0D+00, &
     9.0D+00, &
    10.0D+00, &
    11.0D+00, &
    12.0D+00, &
    13.0D+00, &
    14.0D+00, &
    15.0D+00 /)
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
subroutine bessel_jn_values ( n_data, nu, x, fx )

!*****************************************************************************80
!
!! BESSEL_JN_VALUES returns some values of the Jn Bessel function.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      BesselJ[n,x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 April 2001
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
!    Output, integer NU, the order of the function.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer, parameter :: n_max = 20

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
     0.1149034849319005D+00, &
     0.3528340286156377D+00, &
     0.4656511627775222D-01, &
     0.2546303136851206D+00, &
    -0.5971280079425882D-01, &
     0.2497577302112344D-03, &
     0.7039629755871685D-02, &
     0.2611405461201701D+00, &
    -0.2340615281867936D+00, &
    -0.8140024769656964D-01, &
     0.2630615123687453D-09, &
     0.2515386282716737D-06, &
     0.1467802647310474D-02, &
     0.2074861066333589D+00, &
    -0.1138478491494694D+00, &
     0.3873503008524658D-24, &
     0.3918972805090754D-18, &
     0.2770330052128942D-10, &
     0.1151336924781340D-04, &
    -0.1167043527595797D+00 /)
  integer n_data
  integer nu
  integer, save, dimension ( n_max ) :: nu_vec = (/ &
     2,  2,  2,  2, &
     2,  5,  5,  5, &
     5,  5, 10, 10, &
    10, 10, 10, 20, &
    20, 20, 20, 20 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
     1.0D+00, & 
     2.0D+00, &
     5.0D+00, &
    10.0D+00, &
    50.0D+00, &
     1.0D+00, &
     2.0D+00, &
     5.0D+00, &
    10.0D+00, &
    50.0D+00, &
     1.0D+00, &
     2.0D+00, &
     5.0D+00, &
    10.0D+00, &
    50.0D+00, &
     1.0D+00, &
     2.0D+00, &
     5.0D+00, &
    10.0D+00, &
    50.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    nu = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    nu = nu_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine bessel_y0_values ( n_data, x, fx )

!*****************************************************************************80
!
!! BESSEL_Y0_VALUES returns some values of the Y0 Bessel function.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      BesselY[0,x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 August 2004
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

  integer, parameter :: n_max = 16

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    -0.1534238651350367D+01, &
     0.8825696421567696D-01, &
     0.5103756726497451D+00, &
     0.3768500100127904D+00, & 
    -0.1694073932506499D-01, & 
    -0.3085176252490338D+00, & 
    -0.2881946839815792D+00, & 
    -0.2594974396720926D-01, & 
     0.2235214893875662D+00, & 
     0.2499366982850247D+00, & 
     0.5567116728359939D-01, & 
    -0.1688473238920795D+00, & 
    -0.2252373126343614D+00, & 
    -0.7820786452787591D-01, & 
     0.1271925685821837D+00, & 
     0.2054642960389183D+00 /)
  integer n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
     0.1D+00, & 
     1.0D+00, & 
     2.0D+00, & 
     3.0D+00, & 
     4.0D+00, & 
     5.0D+00, & 
     6.0D+00, & 
     7.0D+00, & 
     8.0D+00, & 
     9.0D+00, & 
    10.0D+00, & 
    11.0D+00, & 
    12.0D+00, & 
    13.0D+00, & 
    14.0D+00, & 
    15.0D+00 /)

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
subroutine bessel_y1_values ( n_data, x, fx )

!*****************************************************************************80
!
!! BESSEL_Y1_VALUES returns some values of the Y1 Bessel function.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      BesselY[1,x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 August 2004
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

  integer, parameter :: n_max = 16

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    -0.6458951094702027D+01, &
    -0.7812128213002887D+00, &
    -0.1070324315409375D+00, &
     0.3246744247918000D+00, &
     0.3979257105571000D+00, &
     0.1478631433912268D+00, &
    -0.1750103443003983D+00, &
    -0.3026672370241849D+00, &
    -0.1580604617312475D+00, &
     0.1043145751967159D+00, &
     0.2490154242069539D+00, &
     0.1637055374149429D+00, &
    -0.5709921826089652D-01, &
    -0.2100814084206935D+00, &
    -0.1666448418561723D+00, &
     0.2107362803687351D-01 /)
  integer n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
     0.1D+00, & 
     1.0D+00, & 
     2.0D+00, & 
     3.0D+00, & 
     4.0D+00, & 
     5.0D+00, & 
     6.0D+00, & 
     7.0D+00, & 
     8.0D+00, & 
     9.0D+00, & 
    10.0D+00, & 
    11.0D+00, & 
    12.0D+00, & 
    13.0D+00, & 
    14.0D+00, & 
    15.0D+00 /)

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
subroutine bessel_yn_values ( n_data, nu, x, fx )

!*****************************************************************************80
!
!! BESSEL_YN_VALUES returns some values of the Yn Bessel function.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      BesselY[n,x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 August 2004
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
!    Output, integer NU, the order of the function.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer, parameter :: n_max = 20

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    -0.1650682606816254D+01, &
    -0.6174081041906827D+00, &
     0.3676628826055245D+00, &
    -0.5868082442208615D-02, &
     0.9579316872759649D-01, &
    -0.2604058666258122D+03, &
    -0.9935989128481975D+01, &
    -0.4536948224911019D+00, &
     0.1354030476893623D+00, &
    -0.7854841391308165D-01, &
    -0.1216180142786892D+09, &
    -0.1291845422080393D+06, &
    -0.2512911009561010D+02, &
    -0.3598141521834027D+00, &
     0.5723897182053514D-02, &
    -0.4113970314835505D+23, &
    -0.4081651388998367D+17, &
    -0.5933965296914321D+09, &
    -0.1597483848269626D+04, &
     0.1644263394811578D-01 /)
  integer n_data
  integer nu
  integer, save, dimension ( n_max ) :: nu_vec = (/ &
     2,  2,  2,  2, &
     2,  5,  5,  5, &
     5,  5, 10, 10, &
    10, 10, 10, 20, &
    20, 20, 20, 20 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
     1.0D+00, & 
     2.0D+00, &
     5.0D+00, &
    10.0D+00, &
    50.0D+00, &
     1.0D+00, &
     2.0D+00, &
     5.0D+00, &
    10.0D+00, &
    50.0D+00, &
     1.0D+00, &
     2.0D+00, &
     5.0D+00, &
    10.0D+00, &
    50.0D+00, &
     1.0D+00, &
     2.0D+00, &
     5.0D+00, &
    10.0D+00, &
    50.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    nu = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    nu = nu_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

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
