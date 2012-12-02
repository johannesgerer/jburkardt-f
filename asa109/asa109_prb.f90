program main

!*****************************************************************************80
!
!! MAIN is the main program for ASA109_PRB.
!
!  Discussion:
!
!    ASA109_PRB calls the ASA109 routines.
!
!  Modified:
!
!    15 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ASA109_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the ASA109 library.'

  call test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ASA109_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 demonstrates the use of XINBTA
!
!  Modified:
!
!    27 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) alngam
  real ( kind = 8 ) b
  real ( kind = 8 ) beta_log
  real ( kind = 8 ) fx
  integer ( kind = 4 ) ifault
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ) x2
  real ( kind = 8 ) xinbta

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01:'
  write ( *, '(a)' ) '  XINBTA inverts the incomplete Beta function.'
  write ( *, '(a)' ) '  Compare with tabulated values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) '      A       B           FX      ', &
    '    X                         X                       DIFF'
  write ( *, '(a,a)' ) '                                  ', &
    '   (tabulated)               (XINBTA)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call beta_inc_values ( n_data, a, b, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    beta_log = alngam ( a, ifault ) &
             + alngam ( b, ifault ) &
             - alngam ( a + b, ifault )

    x2 = xinbta ( a, b, beta_log, fx, ifault )

    write ( *, &
    '(2x,f6.2,2x,f6.2,2x,f14.6,2x,g24.16,2x,g24.16,2x,g10.4)' ) &
    a, b, fx, x, x2, abs ( x - x2 )

  end do

  return
end
