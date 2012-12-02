program main

!*****************************************************************************80
!
!! MAIN is the main program for ASA063_PRB.
!
!  Discussion:
!
!    ASA063_PRB calls the ASA063 routines.
!
!  Modified:
!
!    12 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ASA063_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the ASA063 library.'

  call test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ASA063_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 demonstrates the use of BETAIN.
!
!  Modified:
!
!    12 January 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) alogam
  real ( kind = 8 ) b
  real ( kind = 8 ) beta_log
  real ( kind = 8 ) betain
  real ( kind = 8 ) fx
  real ( kind = 8 ) fx2
  integer ( kind = 4 ) ifault
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01:'
  write ( *, '(a)' ) &
    '  BETAIN evaluates the incomplete Beta function.'
  write ( *, '(a)' ) '  Compare with tabulated values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) '      A       B       X      ', &
    ' BetaInc                   BetaInc                 DIFF'
  write ( *, '(a,a)' ) '                             ', &
    '(tabulated)               (BETAIN)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call beta_inc_values ( n_data, a, b, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    beta_log = alogam ( a, ifault ) &
             + alogam ( b, ifault ) &
             - alogam ( a + b, ifault )

    fx2 = betain ( x, a, b, beta_log, ifault )

    write ( *, &
    '(2x,f6.2,2x,f6.2,2x,f6.2,2x,g24.16,2x,g24.16,2x,g10.4)' ) &
    a, b, x, fx, fx2, abs ( fx - fx2 )

  end do

  return
end
