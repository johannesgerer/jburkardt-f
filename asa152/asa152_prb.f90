program main

!*****************************************************************************80
!
!! MAIN is the main program for ASA152_PRB.
!
!  Discussion:
!
!    ASA152_PRB calls the ASA152 routines.
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

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ASA152_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the ASA152 library.'

  call test01 ( )
  call test02 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ASA152_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 demonstrates CHYPER for cumulative probabilities.
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

  real ( kind = 8 ) chyper
  real ( kind = 8 ) fx
  real ( kind = 8 ) fx2
  integer ( kind = 4 ) ifault
  integer ( kind = 4 ) n_data
  logical point
  integer ( kind = 4 ) pop
  integer ( kind = 4 ) sam
  integer ( kind = 4 ) suc
  integer ( kind = 4 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01:'
  write ( *, '(a)' ) '  CHYPER computes cumulative probabilities'
  write ( *, '(a)' ) '  of the hypergeometric PDF.'
  write ( *, '(a)' ) '  Compare with tabulated values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) '   SAM   SUC   POP     X    ', &
    '  CDF                       CDF                     DIFF'
  write ( *, '(a,a)' ) '                            ', &
    ' (tabulated)               (CHYPER)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call hypergeometric_cdf_values ( n_data, sam, suc, pop, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    point = .false.
    fx2 = chyper ( point, sam, x, pop, suc, ifault )

    write ( *, &
    '(2x,i4,2x,i4,2x,i4,2x,i4,2x,g24.16,2x,g24.16,2x,g10.4)' ) &
    sam, suc, pop, x, fx, fx2, abs ( fx - fx2 )

  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 demonstrates CHYPER for point probabilities.
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

  real ( kind = 8 ) chyper
  real ( kind = 8 ) fx
  real ( kind = 8 ) fx2
  integer ( kind = 4 ) ifault
  integer ( kind = 4 ) n_data
  logical point
  integer ( kind = 4 ) pop
  integer ( kind = 4 ) sam
  integer ( kind = 4 ) suc
  integer ( kind = 4 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02:'
  write ( *, '(a)' ) '  CHYPER computes point probabilities'
  write ( *, '(a)' ) '  of the hypergeometric PDF.'
  write ( *, '(a)' ) '  Compare with tabulated values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) '   SAM   SUC   POP     X', &
    '      PDF                       PDF                    DIFF'
  write ( *, '(a,a)' ) '                        ', &
    '     (tabulated)                (CHYPER)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call hypergeometric_pdf_values ( n_data, sam, suc, pop, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    point = .true.
    fx2 = chyper ( point, sam, x, pop, suc, ifault )

    write ( *, &
    '(2x,i4,2x,i4,2x,i4,2x,i4,2x,g24.16,2x,g24.16,2x,g10.4)' ) &
    sam, suc, pop, x, fx, fx2, abs ( fx - fx2 )

  end do

  return
end
