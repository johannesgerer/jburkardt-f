program main

!*****************************************************************************80
!
!! MAIN is the main program for ASA111_PRB.
!
!  Modified:
!
!    21 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ASA111_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the ASA111 library.'

  call test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ASA111_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 compares PPND against tabulated values.
!
!  Modified:
!
!    21 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) ifault
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) ppnd
  real ( kind = 8 ) x
  real ( kind = 8 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01:'
  write ( *, '(a)' ) '  Compare PPND against tabulated values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) &
    '      CDF              X                         X2    ', &
    '                 DIFF'
  write ( *, '(a)' ) &
    '                      (tabulated)               (PPND)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call normal_01_cdf_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    x2 = ppnd ( fx, ifault )

    write ( *, '(2x,g14.6,2x,g24.16,2x,g24.16,2x,g10.4)' ) &
    fx, x, x2, abs ( x - x2 )

  end do

  return
end
