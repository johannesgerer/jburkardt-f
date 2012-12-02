program main

!*****************************************************************************80
!
!! MAIN is the main program for ASA245_PRB.
!
!  Discussion:
!
!    ASA245_PRB calls the ASA245 routines.
!
!  Modified:
!
!    13 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ASA245_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the ASA245 library.'

  call test01 ( )
  call test02 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ASA245_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 demonstrates the use of ALNGAM.
!
!  Modified:
!
!    13 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) alngam
  real ( kind = 8 ) fx
  real ( kind = 8 ) fx2
  integer ( kind = 4 ) ifault
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01:'
  write ( *, '(a)' ) '  ALNGAM computes the logarithm of the '
  write ( *, '(a)' ) '  Gamma function.  We compare the result'
  write ( *, '(a)' ) '  to tabulated values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) '          X                     ', &
    'FX                        FX2'
  write ( *, '(a,a)' ) '                                ', &
    '(Tabulated)               (ALNGAM)                DIFF'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call gamma_log_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = alngam ( x, ifault )

    write ( *, '(2x,f24.16,2x,g24.16,2x,g24.16,2x,g10.4)' ) &
    x, fx, fx2, abs ( fx - fx2 )

  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 demonstrates the use of LNGAMMA.
!
!  Modified:
!
!    13 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  real ( kind = 8 ) fx2
  integer ( kind = 4 ) ier
  real ( kind = 8 ) lngamma
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02:'
  write ( *, '(a)' ) '  LNGAMMA computes the logarithm of the '
  write ( *, '(a)' ) '  Gamma function.  We compare the result'
  write ( *, '(a)' ) '  to tabulated values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) '          X                     ', &
    'FX                        FX2'
  write ( *, '(a,a)' ) '                                ', &
    '(Tabulated)               (LNGAMMA)                DIFF'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call gamma_log_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = lngamma ( x, ier )

    write ( *, '(2x,f24.16,2x,g24.16,2x,g24.16,2x,g10.4)' ) &
    x, fx, fx2, abs ( fx - fx2 )

  end do

  return
end
