program main

!*****************************************************************************80
!
!! MAIN is the main program for ASA076_PRB.
!
!  Discussion:
!
!    ASA076_PRB calls the ASA076 routines.
!
!  Modified:
!
!    16 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ASA076_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the ASA076 library.'

  call test01 ( )
  call test02 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ASA076_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 demonstrates the use of TFN.
!
!  Modified:
!
!    16 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) h
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) tfn

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01:'
  write ( *, '(a)' ) '  TFN computes the Owen T function.'
  write ( *, '(a)' ) '  Compare with tabulated values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) '          H               A           ', &
  'T                         T                       DIFF'
  write ( *, '(a,a)' ) '                                     ', &
  '(Tabulated)               (TFN)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call owen_values ( n_data, h, a, t1 )

    if ( n_data == 0 ) then
      exit
    end if

    t2 = tfn ( h, a )

    write ( *, '(2x,f14.6,2x,f14.6,2x,g24.16,2x,g24.16,2x,g10.4)' ) &
    h, a, t1, t2, abs ( t1 - t2 )

  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 demonstrates the use of THA.
!
!  Modified:
!
!    16 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) h
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) tha

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02:'
  write ( *, '(a)' ) '  THA computes the Owen T function.'
  write ( *, '(a)' ) '  Compare with tabulated values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) '          H               A           ', &
 'T                         T                       DIFF'
  write ( *, '(a,a)' ) '                                     ', &
 '(Tabulated)               (THA)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call owen_values ( n_data, h, a, t1 )

    if ( n_data == 0 ) then
      exit
    end if

    t2 = tha ( h, 1.0D+00, a, 1.0D+00 )

    write ( *, '(2x,f14.6,2x,f14.6,2x,g24.16,2x,g24.16,2x,g10.4)' ) &
    h, a, t1, t2, abs ( t1 - t2 )

  end do

  return
end
