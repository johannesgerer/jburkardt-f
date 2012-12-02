program main

!*****************************************************************************80
!
!! MAIN is the main program for WEEKDAY_PRB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 March 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'WEEKDAY_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the WEEKDAY library.'

  call test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'WEEKDAY_PRB:'
  write ( *, '(a)' ) '  Noraml end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests YMD_TO_WEEKDAY_COMMON.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 March 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n_data
  character ( len = 9 )  s1
  character ( len = 9 )  s2
  character ( len = 20 ) s3
  integer ( kind = 4 ) w1
  integer ( kind = 4 ) w2
  integer ( kind = 4 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  For dates in the Common calendar:'
  write ( *, '(a)' ) '  YMD_TO_WEEKDAY_COMMON returns the day of the week.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  YMD                   Weekday    Weekday'
  write ( *, '(a)' ) '                        Tabulated  Computed'
  write ( *, '(a)' ) ' '

  do

    call weekday_values ( n_data, y, m, d, w1 )

    if ( n_data == 0 ) then
      exit
    end if
 
    call ymd_to_s_common ( y, m, d, s3 ) 
    call ymd_to_weekday_common ( y, m, d, w2 )
    call weekday_to_name_common ( w1, s1 )
    call weekday_to_name_common ( w2, s2 )

    write ( *, '(2x,a20,2x,a9,2x,a9)' ) s3, s1, s2

  end do

  return
end
