program main

!*****************************************************************************80
!
!! MAIN is the main program for ASA066_PRB.
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
  write ( *, '(a)' ) 'ASA066_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the ASA066 library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ASA066_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 compares ALNORM against tabulated values.
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

  real ( kind = 8 ) alnorm
  real ( kind = 8 ) fx
  real ( kind = 8 ) fx2
  integer ( kind = 4 ) n_data
  logical upper
  real ( kind = 8 ) x

  upper = .false.

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01:'
  write ( *, '(a)' ) '  Compare tabulated values of the normal '
  write ( *, '(a)' ) '  Cumulative Density Function against values'
  write ( *, '(a)' ) '  computed by ALNORM.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) &
    '         X        CDF                       CDF', &
    '                    DIFF'
  write ( *, '(a)' ) &
    '               (tabulated)                 (ALNORM)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do
    call normal_01_cdf_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = alnorm ( x, upper )

    write ( *, '(2x,f10.4,2x,g24.16,2x,g24.16,2x,g10.4)' ) &
    x, fx, fx2, abs ( fx - fx2 )

  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 compares NORMP against tabulated values.
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
  real ( kind = 8 ) fx3
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) pdf
  logical upper
  real ( kind = 8 ) x

  upper = .false.

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02:'
  write ( *, '(a)' ) '  Compare tabulated values of the normal '
  write ( *, '(a)' ) '  Cumulative Density Function against values'
  write ( *, '(a)' ) '  computed by NORMP.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) &
    '         X        CDF                       CDF', &
    '                    DIFF'
  write ( *, '(a)' ) &
    '               (tabulated)                 (NORMP)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call normal_01_cdf_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    call normp ( x, fx2, fx3, pdf )

    write ( *, '(2x,f10.4,2x,g24.16,2x,g24.16,2x,g10.4)' ) &
    x, fx, fx2, abs ( fx - fx2 )

  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 compares NPROB against tabulated values.
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
  real ( kind = 8 ) fx3
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) pdf
  logical upper
  real ( kind = 8 ) x

  upper = .false.

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03:'
  write ( *, '(a)' ) '  Compare tabulated values of the normal '
  write ( *, '(a)' ) '  Cumulative Density Function against values'
  write ( *, '(a)' ) '  computed by NPROB.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) &
    '         X        CDF                       CDF', &
    '                    DIFF'
  write ( *, '(a)' ) &
    '               (tabulated)                 (NPROB)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call normal_01_cdf_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    call nprob ( x, fx2, fx3, pdf )

    write ( *, '(2x,f10.4,2x,g24.16,2x,g24.16,2x,g10.4)' ) &
    x, fx, fx2, abs ( fx - fx2 )

  end do

  return
end
