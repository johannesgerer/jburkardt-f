program main

!*****************************************************************************80
!
!! MAIN is the main program for APPORTIONMENT_PRB.
!
!  Discussion:
!
!    APPORTIONMENT_PRB calls the routines in the APPORTIONMENT library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 June 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'APPORTIONMENT_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the routines in APPORTIONMENT.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
  call test05 ( )
  call test06 ( )
  call test07 ( )
  call test08 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'APPORTIONMENT_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests STATE_NUM_YEAR and rep_NUM_YEAR.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) rep_num
  integer ( kind = 4 ) state_num
  integer ( kind = 4 ) year

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  STATE_NUM_YEAR returns the number of states in'
  write ( *, '(a)' ) '  the union at the end of a given year.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  REP_NUM_YEAR returns the number of reps in'
  write ( *, '(a)' ) '  the House of Representatives (based only on the'
  write ( *, '(a)' ) '  decennial census.)'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Year  States  reps'
  write ( *, '(a)' ) ' '

  do year = 1790, 2010, 10
    call state_num_year ( year, state_num )
    call rep_num_year ( year, rep_num )
    write ( *, '(2x,i4,6x,i2,4x,i3)' ) year, state_num, rep_num
  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests STATE_STATEHOOD.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) m
  character ( len = 9 ) month
  integer ( kind = 4 ) state
  character ( len = 20 ) state_name
  integer ( kind = 4 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  STATE_STATEHOOD returns the statehood date of a state.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   #  Name                     Statehood Date'
  write ( *, '(a)' ) ' '

  do state = 1, 51
    call state_statehood ( state, y, m, d )
    call i4_to_month_name ( m, month )
    write ( *, '(2x,i2,2x,a20,2x,i2,2x,a9,2x,i4)' ) &
      state, state_name(state), d, month, y
  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests APPORTION_HAMILTON.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: state_num = 51

  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) rep_num
  integer ( kind = 4 ) state
  character ( len = 2 ) state_id
  integer ( kind = 4 ) state_pop(state_num)
  integer ( kind = 4 ) state_rep(state_num)
  character ( len = 12 ) string
  integer ( kind = 4 ) us_pop
  integer ( kind = 4 ) year

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  Get the historic representation values.'
!
!  Pick a year.
!
  year = 1960
!
!  What were the state populations in the last decennial census?
!
  call state_pop_year ( year, state_pop )
!
!  What were the state representations based on the last decennial census?
!
  call state_rep_year ( year, state_rep )

  call rep_num_year ( year, rep_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i4)' ) '  Year: ', year
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  ID    Population   Rep     Pop%      Rep%'
  write ( *, '(a)' ) ' '

  us_pop = sum ( state_pop(1:state_num) )

  do state = 1, state_num

    f1 = real ( 100.0 * state_pop(state), kind = 8 ) / real ( us_pop, kind = 8 )
    f2 = real ( 100.0 * state_rep(state), kind = 8 ) / real ( rep_num, kind = 8 )

    call i4_to_s_right_comma ( state_pop(state), string )

    write ( *, '(2x,a2,2x,a12,2x,i3,2x,f8.4,2x,f8.4)' ) &
      state_id (state), string, state_rep(state), f1, f2

  end do

  write ( *, '(a)' ) '  --  ------------  ---  --------  --------'

  call i4_to_s_right_comma ( us_pop, string )

  write ( *, '(2x,a2,2x,a12,2x,i3,2x,f8.4,2x,f8.4)' ) &
    'US', string, rep_num, 100.0, 100.0

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests APPORTION_HAMILTON.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 May 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  real ( kind = 8 ) f3
  integer ( kind = 4 ) indx(51)
  integer ( kind = 4 ) rep_num
  integer ( kind = 4 ) state
  character ( len = 2 ) state_id
  integer ( kind = 4 ) state_num
  integer ( kind = 4 ) state_pop(51)
  integer ( kind = 4 ), allocatable, dimension ( : ) :: state_rep
  character ( len = 12 ) string
  integer ( kind = 4 ) us_pop
  integer ( kind = 4 ) year

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  APPORTION_HAMILTON uses Hamilton''s method to'
  write ( *, '(a)' ) '  apportion representatives.'
!
!  Pick a year.
!
  year = 1960
!
!  What were the state populations in the last decennial census?
!
  call state_pop_year ( year, state_pop )
!
!  Make an index vector.
!
  call i4vec_nonzero_first ( 51, state_pop, state_num, indx )
!
!  "Squeeze" the population vector.
!
  state_pop(1:state_num) = state_pop(indx(1:state_num))

  allocate ( state_rep(1:state_num) )

  call rep_num_year ( year, rep_num )

  call apportion_hamilton ( state_num, state_pop, rep_num, state_rep )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i4)' ) '  Year: ', year
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  ID    Population   Rep     Pop%      Rep%      Pop/Rep'
  write ( *, '(a)' ) ' '
  us_pop = sum ( state_pop(1:state_num) )

  do state = 1, state_num

    f1 = real ( 100.0 * state_pop(state), kind = 8 ) / real ( us_pop, kind = 8 )
    f2 = real ( 100.0 * state_rep(state), kind = 8 ) / real ( rep_num, kind = 8 )
    f3 = real ( state_pop(state), kind = 8 ) / real ( state_rep(state), kind = 8 )

    call i4_to_s_right_comma ( state_pop(state), string )

    write ( *, '(2x,a2,2x,a12,2x,i3,2x,f8.4,2x,f8.4,2x,f12.0)' ) &
      state_id ( indx(state) ), string, state_rep(state), f1, f2, f3

  end do

  write ( *, '(a)' ) '  --  ------------  ---  --------  --------  ------------'

  call i4_to_s_right_comma ( us_pop, string )

  f3 = real ( us_pop, kind = 8 ) / real ( rep_num, kind = 8 )

  write ( *, '(2x,a2,2x,a12,2x,i3,2x,f8.4,2x,f8.4,2x,f12.0)' ) &
    'US', string, rep_num, 100.0, 100.0, f3

  deallocate ( state_rep )

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 tests APPORTION_JEFFERSON.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  real ( kind = 8 ) f3
  integer ( kind = 4 ) indx(51)
  integer ( kind = 4 ) rep_num
  integer ( kind = 4 ) state
  character ( len = 2 ) state_id
  integer ( kind = 4 ) state_num
  integer ( kind = 4 ) state_pop(51)
  integer ( kind = 4 ), allocatable, dimension ( : ) :: state_rep
  character ( len = 12 ) string
  integer ( kind = 4 ) us_pop
  integer ( kind = 4 ) year

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  APPORTION_JEFFERSON uses Jefferson''s method to'
  write ( *, '(a)' ) '  apportion representatives.'
!
!  Pick a year.
!
  year = 1960
!
!  What were the state populations in the last decennial census?
!
  call state_pop_year ( year, state_pop )
!
!  Make an index vector.
!
  call i4vec_nonzero_first ( 51, state_pop, state_num, indx )
!
!  "Squeeze" the population vector.
!
  state_pop(1:state_num) = state_pop(indx(1:state_num))

  allocate ( state_rep(1:state_num) )

  call rep_num_year ( year, rep_num )

  call apportion_jefferson ( state_num, state_pop, rep_num, state_rep )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i4)' ) '  Year: ', year
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  ID    Population   Rep     Pop%      Rep%      Pop/Rep'
  write ( *, '(a)' ) ' '
  us_pop = sum ( state_pop(1:state_num) )

  do state = 1, state_num

    f1 = real ( 100.0 * state_pop(state), kind = 8 ) / real ( us_pop, kind = 8 )
    f2 = real ( 100.0 * state_rep(state), kind = 8 ) / real ( rep_num, kind = 8 )
    f3 = real ( state_pop(state), kind = 8 ) / real ( state_rep(state), kind = 8 )

    call i4_to_s_right_comma ( state_pop(state), string )

    write ( *, '(2x,a2,2x,a12,2x,i3,2x,f8.4,2x,f8.4,2x,f12.0)' ) &
      state_id ( indx(state) ), string, state_rep(state), f1, f2, f3

  end do

  write ( *, '(a)' ) '  --  ------------  ---  --------  --------  ------------'

  call i4_to_s_right_comma ( us_pop, string )

  f3 = real ( us_pop, kind = 8 ) / real ( rep_num, kind = 8 )

  write ( *, '(2x,a2,2x,a12,2x,i3,2x,f8.4,2x,f8.4,2x,f12.0)' ) &
    'US', string, rep_num, 100.0, 100.0, f3

  deallocate ( state_rep )

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 tests APPORTION_ADAMS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 May 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  real ( kind = 8 ) f3
  integer ( kind = 4 ) indx(51)
  integer ( kind = 4 ) rep_num
  integer ( kind = 4 ) state
  character ( len = 2 ) state_id
  integer ( kind = 4 ) state_num
  integer ( kind = 4 ) state_pop(51)
  integer ( kind = 4 ), allocatable, dimension ( : ) :: state_rep
  character ( len = 12 ) string
  integer ( kind = 4 ) us_pop
  integer ( kind = 4 ) year

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  APPORTION_ADAMS uses Adams''s method to'
  write ( *, '(a)' ) '  apportion representatives.'
!
!  Pick a year.
!
  year = 1960
!
!  What were the state populations in the last decennial census?
!
  call state_pop_year ( year, state_pop )
!
!  Make an index vector.
!
  call i4vec_nonzero_first ( 51, state_pop, state_num, indx )
!
!  "Squeeze" the population vector.
!
  state_pop(1:state_num) = state_pop(indx(1:state_num))

  allocate ( state_rep(1:state_num) )

  call rep_num_year ( year, rep_num )

  call apportion_adams ( state_num, state_pop, rep_num, state_rep )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i4)' ) '  Year: ', year
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  ID    Population   Rep     Pop%      Rep%      Pop/Rep'
  write ( *, '(a)' ) ' '
  us_pop = sum ( state_pop(1:state_num) )

  do state = 1, state_num

    f1 = real ( 100.0 * state_pop(state), kind = 8 ) / real ( us_pop, kind = 8 )
    f2 = real ( 100.0 * state_rep(state), kind = 8 ) / real ( rep_num, kind = 8 )
    f3 = real ( state_pop(state), kind = 8 ) / real ( state_rep(state), kind = 8 )

    call i4_to_s_right_comma ( state_pop(state), string )

    write ( *, '(2x,a2,2x,a12,2x,i3,2x,f8.4,2x,f8.4,2x,f12.0)' ) &
      state_id ( indx(state) ), string, state_rep(state), f1, f2, f3

  end do

  write ( *, '(a)' ) '  --  ------------  ---  --------  --------  ------------'

  call i4_to_s_right_comma ( us_pop, string )

  f3 = real ( us_pop, kind = 8 ) / real ( rep_num, kind = 8 )

  write ( *, '(2x,a2,2x,a12,2x,i3,2x,f8.4,2x,f8.4,2x,f12.0)' ) &
    'US', string, rep_num, 100.0, 100.0, f3

  deallocate ( state_rep )

  return
end
subroutine test07 ( )

!*****************************************************************************80
!
!! TEST07 tests APPORTION_WEBSTER.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 May 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  real ( kind = 8 ) f3
  integer ( kind = 4 ) indx(51)
  integer ( kind = 4 ) rep_num
  integer ( kind = 4 ) state
  character ( len = 2 ) state_id
  integer ( kind = 4 ) state_num
  integer ( kind = 4 ) state_pop(51)
  integer ( kind = 4 ), allocatable, dimension ( : ) :: state_rep
  character ( len = 12 ) string
  integer ( kind = 4 ) us_pop
  integer ( kind = 4 ) year

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a)' ) '  APPORTION_WEBSTER uses Webster''s method to'
  write ( *, '(a)' ) '  apportion representatives.'
!
!  Pick a year.
!
  year = 1960
!
!  What were the state populations in the last decennial census?
!
  call state_pop_year ( year, state_pop )
!
!  Make an index vector.
!
  call i4vec_nonzero_first ( 51, state_pop, state_num, indx )
!
!  "Squeeze" the population vector.
!
  state_pop(1:state_num) = state_pop(indx(1:state_num))

  allocate ( state_rep(1:state_num) )

  call rep_num_year ( year, rep_num )

  call apportion_webster ( state_num, state_pop, rep_num, state_rep )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i4)' ) '  Year: ', year
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  ID    Population   Rep     Pop%      Rep%      Pop/Rep'
  write ( *, '(a)' ) ' '
  us_pop = sum ( state_pop(1:state_num) )

  do state = 1, state_num

    f1 = real ( 100.0 * state_pop(state), kind = 8 ) / real ( us_pop, kind = 8 )
    f2 = real ( 100.0 * state_rep(state), kind = 8 ) / real ( rep_num, kind = 8 )
    f3 = real ( state_pop(state), kind = 8 ) / real ( state_rep(state), kind = 8 )

    call i4_to_s_right_comma ( state_pop(state), string )

    write ( *, '(2x,a2,2x,a12,2x,i3,2x,f8.4,2x,f8.4,2x,f12.0)' ) &
      state_id ( indx(state) ), string, state_rep(state), f1, f2, f3

  end do

  write ( *, '(a)' ) '  --  ------------  ---  --------  --------  ------------'

  call i4_to_s_right_comma ( us_pop, string )

  f3 = real ( us_pop, kind = 8 ) / real ( rep_num, kind = 8 )

  write ( *, '(2x,a2,2x,a12,2x,i3,2x,f8.4,2x,f8.4,2x,f12.0)' ) &
    'US', string, rep_num, 100.0, 100.0, f3

  deallocate ( state_rep )

  return
end
subroutine test08 ( )

!*****************************************************************************80
!
!! TEST08 tests APPORTION_HUNTINGTON_HILL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 June 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: state_num = 15

  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  real ( kind = 8 ) f3
  integer ( kind = 4 ), parameter :: rep_num = 105
  integer ( kind = 4 ) state
  character ( len = 2 ), dimension ( state_num ) :: state_id = (/ &
    'CT', 'DE', 'GA', 'KY', 'MD', &
    'MA', 'NH', 'NJ', 'NY', 'NC', &
    'PA', 'RI', 'SC', 'VT', 'VA' /) 
  integer ( kind = 4 ), dimension ( state_num ) :: state_pop = (/ &
    236841,  55540,  70835,  68705, 278514, &
    475327, 141822, 179570, 331589, 353523, &
    432879,  68446, 206236,  85533, 630560 /)	
  integer ( kind = 4 ) state_rep(state_num)
  character ( len = 12 ) string
  integer ( kind = 4 ) us_pop
  integer ( kind = 4 ), parameter :: year = 1790

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08'
  write ( *, '(a)' ) '  APPORTION_HUNTINGON_HILL uses the Huntington-Hill'
  write ( *, '(a)' ) '  apportionment method.'

  call apportion_huntington_hill ( state_num, state_pop, rep_num, state_rep )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i4)' ) '  Year: ', year
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  ID    Population   Rep     Pop%      Rep%      Pop/Rep'
  write ( *, '(a)' ) ' '
  us_pop = sum ( state_pop(1:state_num) )

  do state = 1, state_num

    f1 = real ( 100.0 * state_pop(state), kind = 8 ) / real ( us_pop, kind = 8 )
    f2 = real ( 100.0 * state_rep(state), kind = 8 ) / real ( rep_num, kind = 8 )
    f3 = real ( state_pop(state), kind = 8 ) / real ( state_rep(state), kind = 8 )

    call i4_to_s_right_comma ( state_pop(state), string )

    write ( *, '(2x,a2,2x,a12,2x,i3,2x,f8.4,2x,f8.4,2x,f12.0)' ) &
      state_id(state), string, state_rep(state), f1, f2, f3

  end do

  write ( *, '(a)' ) '  --  ------------  ---  --------  --------  ------------'

  call i4_to_s_right_comma ( us_pop, string )

  f3 = real ( us_pop, kind = 8 ) / real ( rep_num, kind = 8 )

  write ( *, '(2x,a2,2x,a12,2x,i3,2x,f8.4,2x,f8.4,2x,f12.0)' ) &
    'US', string, rep_num, 100.0, 100.0, f3

  return
end
