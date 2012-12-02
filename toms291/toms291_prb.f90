program main

!*****************************************************************************80
!
!! MAIN is the main program for TOMS291_PRB.
!
!  Modified:
!
!    28 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  write ( *, '(a)' ) ' '
  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TOMS291_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the TOMS291 library.'

  call test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TOMS291_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 compares results from ALOGAM with tabulated data.
!
!  Modified:
!
!    07 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real    ( kind = 8 ) alogam
  real    ( kind = 8 ) fx1
  real    ( kind = 8 ) fx2
  integer ( kind = 4 ) ifault
  integer ( kind = 4 ) n_data
  real    ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01:'
  write ( *, '(a)' ) &
  '  ALOGAM computes the logarithm of the Gamma function.'
  write ( *, '(a)' ) '  Compare against tabulated data.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) '          X        ', &
    '   FX                        FX'
  write ( *, '(a,a)' ) '                   ', &
    'Tabulated                  ALOGAM'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call gamma_log_values ( n_data, x, fx1 )

    if ( n_data .eq. 0 ) then
       exit
    end if

    fx2 = alogam ( x, ifault )

    write ( *, '(2x,f14.6,2x,g24.16,2x,g24.16)' ) x, fx1, fx2

  end do

  return
end
