program main

!*****************************************************************************80
!
!! MAIN is the main program for ASA032_PRB.
!
!  Discussion:
!
!    ASA032_PRB calls the ASA032 routines.
!
!  Modified:
!
!    17 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  write ( *, '(a)' ) ' '
  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ASA032_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the ASA032 library.'

  call test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ASA032_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 demonstrates the use of GAMAIN.
!
!  Modified:
!
!    17 January 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) fx
  real ( kind = 8 ) fx2
  real ( kind = 8 ) gamain
  integer ( kind = 4 ) ifault
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01:'
  write ( *, '(a)' ) '  GAMAIN computes the incomplete Gamma '
  write ( *, '(a)' ) '  function.  Compare to tabulated values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) &
    '       A            X          ', &
    '   FX                        FX                      DIFF'
  write ( *, '(a,a)' ) &
    '                               ', &
    ' (tabulated)                (GAMAIN)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call gamma_inc_values ( n_data, a, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = gamain ( x, a, ifault )

    write ( *, &
    '(2x,f12.8,2x,f12.8,2x,g24.16,2x,g24.16,2x,g10.4)' ) &
    a, x, fx, fx2, abs ( fx - fx2 )

  end do

  return
end
