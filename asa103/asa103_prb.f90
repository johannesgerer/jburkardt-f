program main

!*****************************************************************************80
!
!! MAIN is the main program for ASA103_PRB.
!
!  Modified:
!
!    18 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ASA103_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Tests the ASA103 library.'

  call test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ASA103_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 compare DIGAMA against tabulated values.
!
!  Modified:
!
!    18 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) digama
  real ( kind = 8 ) fx
  real ( kind = 8 ) fx2
  integer ( kind = 4 ) ifault
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01:'
  write ( *, '(a)' ) '  Compare tabulated values of the digamma'
  write ( *, '(a)' ) '  or Psi function against values computed'
  write ( *, '(a)' ) '  computed by DIGAMA.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) &
    '         X         Psi                     ', &
    '  Psi                    DIFF'
  write ( *, '(a,a)' ) &
    '               (tabulated)                 ', &
    '(DIGAMA)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call psi_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = digama ( x, ifault )

    write ( *, '(2x,f10.4,2x,g24.16,2x,g24.16,2x,g10.4)' ) &
    x, fx, fx2, abs ( fx - fx2 )

  end do

  return
end
