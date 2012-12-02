program main

!*****************************************************************************80
!
!! MAIN is the main program for ASA239_PRB.
!
!  Discussion:
!
!    ASA239_PRB calls the ASA239 routines.
!
!  Modified:
!
!    20 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ASA239_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the ASA239 library.'

  call test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ASA239_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 demonstrates the use of GAMMAD.
!
!  Modified:
!
!    20 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) fx
  real ( kind = 8 ) fx2
  real ( kind = 8 ) gammad
  integer ( kind = 4 ) ifault
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01:'
  write ( *, '(a)' ) '  GAMMAD computes the incomplete Gamma '
  write ( *, '(a)' ) '  function.  Compare to tabulated values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) &
  '      A             X          ', &
  '   FX                        FX                      DIFF'
  write ( *, '(a,a)' ) &
  '                               ', &
  ' (tabulated)                (GAMMAD)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call gamma_inc_values ( n_data, a, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = gammad ( x, a, ifault )

    write ( *, '(2x,f12.8,2x,f12.8,2x,g24.16,2x,g24.16,2x,g10.4)' ) &
    a, x, fx, fx2, abs ( fx - fx2 )

  end do

  return
end
