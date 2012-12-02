program main

!*****************************************************************************80
!
!! CALPAK_PRB calls a series of tests for CALPAK.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 October 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CALPAK_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the CALPAK library.'

  call test0005 ( )
  call test0007 ( )

  call test001 ( )
  call test002 ( )
  call test003 ( )
  call test004 ( )
  call test005 ( )
  call test006 ( )
  call test0065 ( )
  call test007 ( )
  call test0075 ( )
  call test00755 ( )
  call test00756 ( )
  call test0076 ( )
  call test0006 ( )
  call test008 ( )
  call test009 ( )

  call test010 ( )
  call test011 ( )
  call test012 ( )
  call test0125 ( )
  call test013 ( )
  call test014 ( )
  call test015 ( )
  call test016 ( )
  call test017 ( )
  call test0175 ( )
  call test018 ( )
  call test019 ( )
  call test020 ( )

  call test165 ( )
  call test17 ( )
  call test175 ( )
  call test18 ( )
  call test185 ( )
  call test19 ( )
  call test195 ( )

  call test20 ( )
  call test21 ( )
  call test215 ( )
  call test22 ( )
  call test23 ( )
  call test24 ( )
  call test25 ( )
  call test255 ( )
  call test26 ( )
  call test265 ( )
  call test27 ( )
  call test275 ( )
  call test28 ( )
  call test29 ( )

  call test30 ( )
  call test31 ( )
  call test315 ( )
  call test32 ( )
  call test325 ( )
  call test326 ( )
  call test327 ( )
  call test328 ( )
  call test33 ( )
  call test335 ( )
  call test336 ( )
  call test337 ( )
  call test34 ( )
  call test344 ( )
  call test345 ( )
  call test35 ( )
  call test36 ( )
  call test365 ( )
  call test37 ( )
  call test373 ( )
  call test375 ( )
  call test376 ( )
  call test38 ( )
  call test389 ( )
  call test39 ( )
  call test394 ( )
  call test395 ( )

  call test40 ( )
  call test41 ( )
  call test415 ( )
  call test42 ( )
  call test43 ( )
  call test435 ( )
  call test44 ( )
  call test445 ( )
  call test45 ( )
  call test46 ( )
  call test47 ( )
  call test48 ( )
  call test49 ( )
  call test492 ( )
  call test493 ( )
  call test495 ( )

  call test50 ( )
  call test501 ( )
  call test502 ( )
  call test503 ( )
  call test51 ( )
  call test515 ( )
  call test5153 ( )
  call test51535 ( )
  call test5154 ( )
  call test5155 ( )
  call test5156 ( )
  call test52 ( )
  call test525 ( )
  call test53 ( )
  call test535 ( )
  call test54 ( )
  call test555 ( )
  call test56 ( )
  call test565 ( )
  call test566 ( )
  call test57 ( )
  call test58 ( )
  call test585 ( )
  call test59 ( )

  call test60 ( )
  call test605 ( )
  call test61 ( )
  call test615 ( )
  call test616 ( )
  call test62 ( )
  call test621 ( )
  call test622 ( )
  call test623 ( )
  call test624 ( )
  call test63 ( )
  call test635 ( )
  call test636 ( )
  call test64 ( )
  call test65 ( )
  call test66 ( )
  call test67 ( )
  call test675 ( )
  call test68 ( )
  call test685 ( )
  call test686 ( )
  call test687 ( )
  call test688 ( )
  call test689 ( )
  call test69 ( )
  call test695 ( )

  call test70 ( )
  call test71 ( )
  call test72 ( )
  call test73 ( )
  call test74 ( )
  call test75 ( )
  call test76 ( )
  call test77 ( )
  call test775 ( )
  call test78 ( )
  call test79 ( )
  call test795 ( )

  call test80 ( )
  call test805 ( )
  call test81 ( )
  call test82 ( )
  call test83 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CALPAK_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test0005 ( )

!*****************************************************************************80
!
!! TEST0005 tests CWS_TO_JED_GPS and JED_TO_CWS_GPS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 October 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) c2
  integer ( kind = 4 ) i
  real ( kind = 8 ) jed_epoch
  real ( kind = 8 ) jed1
  real ( kind = 8 ) jed3
  character ( len = 25 ) s2
  real ( kind = 8 ) sec2
  integer ( kind = 4 ) w2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0005'
  write ( *, '(a)' ) '  For the GPS calendar:'
  write ( *, '(a)' ) '  JED_TO_CWS_GPS: JED -> CWS.'
  write ( *, '(a)' ) '  CWS_TO_JED_GPS: CWS -> JED.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   JED (in)       CWS                            JED (out)'
  write ( *, '(a)' ) ' '

  call epoch_to_jed_gps ( jed_epoch )

  i = 0

  do

    i = i + 1
    call jed_test ( i, jed1 )

    if ( jed1 < 0.0D+00 ) then
      exit
    end if

    if ( jed_epoch <= jed1 ) then

      call jed_to_cws_gps ( jed1, c2, w2, sec2 )

      call cws_to_s_gps ( c2, w2, sec2, s2 )

      call cws_to_jed_gps ( c2, w2, sec2, jed3 )

      write ( *, '(2x,f11.2,5x,a,5x,f11.2)' ) jed1, s2, jed3

    end if

  end do

  return
end
subroutine test0007 ( )

!*****************************************************************************80
!
!! TEST0007 tests DAY_LIST_COMMON
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 September 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) :: d1 = 25
  integer ( kind = 4 ) :: d2 = 2
  integer ( kind = 4 ) :: m1 = 9
  integer ( kind = 4 ) :: m2 = 10
  character ( len = 20 ) s
  integer ( kind = 4 ) :: y1 = 2006
  integer ( kind = 4 ) :: y2 = 2006

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0007'
  write ( *, '(a)' ) '  DAY_LIST_COMMON prints a list of days between'
  write ( *, '(a)' ) '  two given YMD dates in the common calendar.'
  write ( *, '(a)' ) ' '
  call ymd_to_s_common ( y1, m1, d1, s )
  write ( *, '(a)' ) '  Initial date: ' // trim ( s )
  call ymd_to_s_common ( y2, m2, d2, s )
  write ( *, '(a)' ) '  Final date:   ' // trim ( s )
  write ( *, '(a)' ) ' '

  call day_list_common ( y1, m1, d1, y2, m2, d2 )

  return
end
subroutine test001 ( )

!*****************************************************************************80
!
!! TEST001 tests EASTER_DS, EASTER_EGR, EASTER_EGR2, EASTER_KNUTH, EASTER_STEWART.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 October 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_test = 10

  integer ( kind = 4 ) d
  integer ( kind = 4 ), dimension ( n_test ) :: d_test = &
    (/  30,    12,    4,  23,   15,   31,   20,   11,   27,   16 /)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  integer ( kind = 4 ), dimension ( n_test ) :: m_test = &
    (/   3,     4,    4,    4,    4,    3,    4,    4,    3,    4 /)
  character ( len = 20 ) s
  integer ( kind = 4 ) y
  integer ( kind = 4 ), dimension ( n_test ) :: y_test = &
    (/ 1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST001'
  write ( *, '(a)' ) '  For the Gregorian calendar,'
  write ( *, '(a)' ) '  for a given year, compute the day and month of Easter.'
  write ( *, '(a)' ) '  EASTER_DS uses Duffett-Smith''s algorithm.'
  write ( *, '(a)' ) '  EASTER_EGR uses Richards''s algorithm.'
  write ( *, '(a)' ) '  EASTER_EGR2 uses Richards''s algorithm 2.'
  write ( *, '(a)' ) '  EASTER_KNUTH uses Knuth''s algorithm.'
  write ( *, '(a)' ) '  EASTER_STEWART uses Stewart''s algorithm.'
 
  do i = 1, n_test

    y = y_test(i)
    m = m_test(i)
    d = d_test(i)

    write ( *, '(a)' ) ' '
    call ymd_to_s_gregorian ( y, m, d, s )
    write ( *, '(a)' ) '  CORRECT:        ' // trim ( s )

    call easter_ds ( y, m, d )
    call ymd_to_s_gregorian ( y, m, d, s )
    write ( *, '(a)' ) '  EASTER_DS:      ' // trim ( s )

    call easter_egr ( y, m, d )
    call ymd_to_s_gregorian ( y, m, d, s )
    write ( *, '(a)' ) '  EASTER_EGR:     ' // trim ( s )

    call easter_egr2 ( y, m, d )
    call ymd_to_s_gregorian ( y, m, d, s )
    write ( *, '(a)' ) '  EASTER_EGR2:    ' // trim ( s )

    call easter_knuth ( y, m, d )
    call ymd_to_s_gregorian ( y, m, d, s )
    write ( *, '(a)' ) '  EASTER_KNUTH:   ' // trim ( s )

    call easter_stewart ( y, m, d )
    call ymd_to_s_gregorian ( y, m, d, s )
    write ( *, '(a)' ) '  EASTER_STEWART: ' // trim ( s )

  end do

  return
end
subroutine test002 ( )

!*****************************************************************************80
!
!! TEST002 tests EASTER_JULIAN and EASTER_JULIAN2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 October 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_test = 10

  integer ( kind = 4 ) d
  integer ( kind = 4 ), dimension ( n_test ) :: d_test = &
    (/  27,    19,   11,  30,   15,    5,   27,   11,    1,   23 /)
  real ( kind = 8 ) f
  integer ( kind = 4 ) i
  real ( kind = 8 ) jed
  integer ( kind = 4 ) m
  integer ( kind = 4 ), dimension ( n_test ) :: m_test = &
    (/   4,     4,    4,    4,    4,    5,    4,    4,    5,    4 /)
  character ( len = 20 ) s
  integer ( kind = 4 ) y
  integer ( kind = 4 ), dimension ( n_test ) :: y_test = &
    (/ 1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST002'
  write ( *, '(a)' ) '  For the Julian calendar,'
  write ( *, '(a)' ) '  for a given year, compute the day and month of Easter.'
  write ( *, '(a)' ) '  EASTER_JULIAN uses Richard''s algorithm.'
  write ( *, '(a)' ) '  EASTER_JULIAN2 uses Richards''s algorithm.'
 
  do i = 1, n_test

    y = y_test(i)
    m = m_test(i)
    d = d_test(i)
    f = 0.5D+00

    write ( *, '(a)' ) ' '
    call ymd_to_s_gregorian ( y, m, d, s )
    write ( *, '(a)' ) '  CORRECT (Gregorian): ' // trim ( s )

    call ymdf_to_jed_gregorian ( y, m, d, f, jed )
    call jed_to_ymdf_julian ( jed, y, m, d, f )

    call ymdf_to_s_julian ( y, m, d, f, s )
    write ( *, '(a)' ) '  CORRECT (Julian):    ' // trim ( s )

    call easter_julian ( y, m, d )
    call ymd_to_s_julian ( y, m, d, s )
    write ( *, '(a)' ) '  EASTER_JULIAN:       ' // trim ( s )

    call easter_julian2 ( y, m, d )
    call ymd_to_s_julian ( y, m, d, s )
    write ( *, '(a)' ) '  EASTER_JULIAN2:      ' // trim ( s )

  end do

  return
end
subroutine test003 ( )

!*****************************************************************************80
!
!! TEST003 tests JED_TO_MAYAN_LONG and MAYAN_LONG_TO_JED.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 October 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) baktun
  real ( kind = 8 ) f
  integer ( kind = 4 ) i
  real ( kind = 8 ) jed_epoch
  real ( kind = 8 ) jed1
  real ( kind = 8 ) jed3
  integer ( kind = 4 ) katun
  integer ( kind = 4 ) kin
  integer ( kind = 4 ) pictun
  integer ( kind = 4 ) tun
  integer ( kind = 4 ) uinal

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST003'
  write ( *, '(a)' ) '  For converting between Julian Ephemeris Dates'
  write ( *, '(a)' ) '  and Mayan Long Count dates:'
  write ( *, '(a)' ) '  JED_TO_MAYAN_LONG,'
  write ( *, '(a)' ) '  MAYAN_LONG_TO_JED.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  JED (in)       MAYAN			   JED (out)'
  write ( *, '(a)' ) '                 P   B   K   T   U   D'
  write ( *, '(a)' ) ' '

  call epoch_to_jed_mayan_long ( jed_epoch )

  i = 0

  do

    i = i + 1
    call jed_test ( i, jed1 )
 
    if ( jed1 < 0.0D+00 ) then
      exit
    end if

    if ( jed_epoch <= jed1 ) then

      call jed_to_mayan_long ( jed1, pictun, baktun, katun, tun, uinal, kin, f )

      call mayan_long_to_jed ( pictun, baktun, katun, tun, uinal, kin, f, jed3 )

      write ( *, '(2x,f11.2,5x,6i4,5x,f11.2)' ) jed1, pictun, baktun, katun, &
        tun, uinal, kin, jed3

      end if

  end do

  return
end
subroutine test004 ( )

!*****************************************************************************80
!
!! TEST004 tests JED_TO_MAYAN_ROUND and MAYAN_ROUND_TO_JED.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 October 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) a2
  integer ( kind = 4 ) b2
  integer ( kind = 4 ) c2
  integer ( kind = 4 ) d2
  real ( kind = 8 ) f2
  integer ( kind = 4 ) i
  real ( kind = 8 ) jed_epoch
  real ( kind = 8 ) jed1
  real ( kind = 8 ) jed3
  integer ( kind = 4 ) y2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST004'
  write ( *, '(a)' ) '  For converting between Julian Ephemeris Dates'
  write ( *, '(a)' ) '  and Mayan Round dates:'
  write ( *, '(a)' ) '  JED_TO_MAYAN_ROUND,'
  write ( *, '(a)' ) '  MAYAN_ROUND_TO_JED.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  JED (in)      MAYAN                       JED (out)'
  write ( *, '(a)' ) '                Y   A   B   C  D  F'
  write ( *, '(a)' ) ' '

  call epoch_to_jed_mayan_long ( jed_epoch )

  i = 0

  do

    i = i + 1
    call jed_test ( i, jed1 )
 
    if ( jed1 < 0.0D+00 ) then
      exit
    end if

    if ( jed_epoch <= jed1 ) then

      call jed_to_mayan_round ( jed1, y2, a2, b2, c2, d2, f2 )

      call mayan_round_to_jed ( y2, a2, b2, c2, d2, f2, jed3 )

      write ( *, '(2x,f11.2,5x,5i4,f5.2,5x,f11.2)' ) &
        jed1, y2, a2, b2, c2, d2, f2, jed3

    end if

  end do

  return
end
subroutine test005 ( )

!*****************************************************************************80
!
!! TEST005 tests JED_TO_WEEKDAY.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 October 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) f2
  real ( kind = 8 ) jed1
  real ( kind = 8 ) jed2
  character ( len = 15 ) s2
  integer ( kind = 4 ) w2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST005'
  write ( *, '(a)' ) '  JED_TO_WEEKDAY reports the day of the week'
  write ( *, '(a)' ) '  for a Julian Ephemeris Date.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   JED     W  Name'
  write ( *, '(a)' ) ' '

  i = 0

  do

    i = i + 1
    call jed_test ( i, jed1 )

    if ( jed1 < 0.0D+00 ) then
      exit
    end if

    call jed_to_next_noon ( jed1, jed2 )

    call jed_to_weekday ( jed2, w2, f2 )
 
    call weekday_to_name_common ( w2, s2 )

    write ( *, '(2x,f11.2,2x,i1,2x,a)' ) jed2, w2, s2

  end do

  return
end
subroutine test006 ( )

!*****************************************************************************80
!
!! TEST006 tests JED_TO_YEAR_HEBREW.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 October 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) jed_epoch
  real ( kind = 8 ) jed1
  character ( len = 10 ) s2
  integer ( kind = 4 ) y2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST006'
  write ( *, '(a)' ) '  For the Hebrew calendar,'
  write ( *, '(a)' ) '  JED_TO_YEAR_HEBREW returns the year of a given JED.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     JED      Hebrew Year'
  write ( *, '(a)' ) ' '

  call epoch_to_jed_hebrew ( jed_epoch )

  i = 0

  do

    i = i + 1
    call jed_test ( i, jed1 )

    if ( jed1 < 0.0D+00 ) then
      exit
    end if
 
    if ( jed_epoch <= jed1 ) then

      call jed_to_year_hebrew ( jed1, y2 )

      call y_to_s_hebrew ( y2, s2 )

      write ( *, '(2x,f11.2,5x,a)' ) jed1, s2

    end if

  end do

  return
end
subroutine test0065 ( )

!*****************************************************************************80
!
!! TEST0065 tests JED_TO_YEARCOUNT_BESSEL and JED_TO_YEARCOUNT_JULIAN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 October 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) bessel
  integer ( kind = 4 ) d
  real ( kind = 8 ) f
  integer ( kind = 4 ) i
  real ( kind = 8 ) jed
  real ( kind = 8 ) julian
  integer ( kind = 4 ) m
  character ( len = 25 ) s
  integer ( kind = 4 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0065'
  write ( *, '(a)' ) '  JED_TO_YEARCOUNT_BESSEL'
  write ( *, '(a)' ) '    returns a tropical year count based at 1900.'
  write ( *, '(a)' ) '  JED_TO_YEARCOUNT_JULIAN'
  write ( *, '(a)' ) '    returns a Julian year count based at 2000.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     JED      YMDF Common        Bessel Year  Julian Year'
  write ( *, '(a)' ) ' '

  i = 0

  do

    i = i + 1
    call jed_test ( i, jed )

    if ( jed < 0 ) then
      exit
    end if

    call jed_to_ymdf_common ( jed, y, m, d, f )
    call ymdf_to_s_common ( y, m, d, f, s )
    call jed_to_yearcount_bessel ( jed, bessel )
    call jed_to_yearcount_julian ( jed, julian )

    write ( *, '(2x,f11.2,5x,a20,2x,2f12.4)' ) jed, s, bessel, julian

  end do

  return
end
subroutine test007 ( )

!*****************************************************************************80
!
!! TEST007 tests JED_TO_YJF_COMMON and YJF_TO_JED_COMMON.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 October 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) f2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j2
  real ( kind = 8 ) jed1
  real ( kind = 8 ) jed3
  character ( len = 20 ) s2
  integer ( kind = 4 ) y2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST007'
  write ( *, '(a)' ) '  For the Common calendar:'
  write ( *, '(a)' ) '  JED_TO_YJF_COMMON: JED -> YJF.'
  write ( *, '(a)' ) '  YJF_TO_JED_COMMON: YJF -> JED.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  JED (in)    YJF 	       JED (out)'
  write ( *, '(a)' ) ' '

  i = 0

  do

    i = i + 1
    call jed_test ( i, jed1 )
    
    if ( jed1 < 0.0D+00 ) then
      exit
    end if

    call jed_to_yjf_common ( jed1, y2, j2, f2 )

    call yjf_to_s_common ( y2, j2, f2, s2 )

    call yjf_to_jed_common ( y2, j2, f2, jed3 )

    write ( *, '(2x,f11.2,5x,a,5x,f11.2)' ) jed1, s2, jed3

  end do

  return
end
subroutine test0075 ( )

!*****************************************************************************80
!
!! TEST0075 tests JED_TO_MJD and MJD_TO_JED.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 October 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) jed1
  real ( kind = 8 ) mjd2
  real ( kind = 8 ) jed3

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0075'
  write ( *, '(a)' ) '  For the modified JED:'
  write ( *, '(a)' ) '  JED_TO_MJD: JED -> MJD.'
  write ( *, '(a)' ) '  MJD_TO_JED: MJD -> JED.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  JED (in)    JEDMOD                JED (out)'
  write ( *, '(a)' ) ' '

  i = 0

  do

    i = i + 1
    call jed_test ( i, jed1 )

    if ( jed1 < 0.0D+00 ) then
      exit
    end if

    call jed_to_mjd ( jed1, mjd2 )

    call mjd_to_jed ( mjd2, jed3 )

    write ( *, '(2x,f11.2,5x,f11.2,5x,f11.2)' ) jed1, mjd2, jed3

  end do

  return
end
subroutine test00755 ( )

!*****************************************************************************80
!
!! TEST00755 tests JED_TO_NYT and NYT_TO_JED.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) d
  real ( kind = 8 ) f
  integer ( kind = 4 ) issue2
  real ( kind = 8 ) jed_nyt_epoch
  real ( kind = 8 ) jed1
  real ( kind = 8 ) jed3
  integer ( kind = 4 ) m
  character ( len = 25 ) s
  integer ( kind = 4 ) volume2
  integer ( kind = 4 ) y

  call epoch_to_jed_nyt ( jed_nyt_epoch )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST00755'
  write ( *, '(a)' ) '  For the New York Times issue date:'
  write ( *, '(a)' ) '  JED_TO_NYT: JED -> NYT.'
  write ( *, '(a)' ) '  NYT_TO_JED: NYT -> JED.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '                           JED (in)      Volume  Issue            JED (out)'
  write ( *, '(a)' ) ' '

  i = 0

  do

    i = i + 1
    call jed_test ( i, jed1 )

    if ( jed1 < 0.0D+00 ) then
      exit
    end if

    if ( jed1 < jed_nyt_epoch ) then
      cycle
    end if

    call jed_to_ymdf_common ( jed1, y, m, d, f )

    call ymdf_to_s_common ( y, m, d, f, s )

    call jed_to_nyt ( jed1, volume2, issue2 )

    call nyt_to_jed ( volume2, issue2, jed3 )

    write ( *, '(2x,a25,2x,f11.2,5x,i4,2x,i8,5x,f11.2)' ) &
      s, jed1, volume2, issue2, jed3

  end do

  return
end
subroutine test00756 ( )

!*****************************************************************************80
!
!! TEST00756 tests NYT_TO_JED.
!
!  Discussion:
!
!    Data (some not used):
!
!     1705   7 March     1857
!     3407  25 August    1862
!     3794  20 November  1863
!     3804   3 December  1863
!    16579  24 February  1903
!    16909  15 March     1904
!    17251  18 April     1905
!    17561  22 February  1906
!    25320  22 May       1927
!    26243  30 November  1929
!    27538  17 June      1933
!    29033  21 June      1937
!    29807   3 September 1939
!    31545   6 June      1945
!    31972   7 August    1945
!    32984  15 May       1948
!    36074  30 October   1956
!    38910   5 August    1964
!    39342  11 October   1965
!    50939   8 October   1997
!    51599  11 December  2000
!    51874  12 September 2001
!    53108  28 January   2005
!    53715  27 September 2006
!    53960  30 May       2007
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 34

  integer ( kind = 4 ), dimension ( test_num ) :: d_test = (/ &
    18, 17, 21, 19, 22,  6,  7, 24, 15, 29, &
    22, 18,  9,  3, 22, 23, 14,  8, 15, 20, &
    16, 15, 21, 18,  9,  6, 17, 14,  8, 31, &
     1, 11, 28, 22 /)
  integer ( kind = 4 ) d1
  integer ( kind = 4 ) d2
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ), dimension ( test_num ) :: issue_test = (/ &
        1,  2155,  2990,  4130,  6189, &
    14499, 15000, 16579, 16909, 17292, &
    17561, 18164, 18856, 21619, 24651, &
    29827, 30000, 31881, 31980, 38864, &
    39317, 40076, 40721, 41418, 44027, &
    44028, 48939, 50000, 50939, 51753, &
    51254, 51599, 53108, 54136 /)
  integer ( kind = 4 ) issue1
  integer ( kind = 4 ) issue3
  real ( kind = 8 ) jed1
  real ( kind = 8 ) jed2
  integer ( kind = 4 ), dimension ( test_num ) :: m_test = (/ &
    9,  8,  4, 12,  7,  2,  2,  2,  3,  5,  &
    2, 10,  9,  4,  7,  9,  3,  5,  8,  6,  &
    9, 10,  7,  6,  8, 11,  4,  3, 10, 12,  &
    1, 12,  1, 11 /)
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  character ( len = 25 ) s1
  character ( len = 25 ) s2
  integer ( kind = 4 ) test
  integer ( kind = 4 ), dimension ( test_num ) :: volume_test = (/ &
       1,    7,   10,   14,   20,   47,   47,   52,   53,   54, &
      55,   57,   58,   66,   74,   89,   89,   94,   94,  113, &
     114,  117,  118,  120,  127,  128,  141,  144,  147,  149, &
     149,  150,  154,  157 /)
  integer ( kind = 4 ) volume1
  integer ( kind = 4 ) volume3
  integer ( kind = 4 ), dimension ( test_num ) :: y_test = (/ &
    1851, 1858, 1861, 1864, 1871, 1898, 1898, 1903, 1904, 1905, &
    1906, 1907, 1909, 1917, 1925, 1939, 1940, 1945, 1945, 1964, &
    1965, 1967, 1969, 1971, 1978, 1978, 1992, 1995, 1997, 1999, &
    2000, 2000, 2005, 2007 /)
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST00756'
  write ( *, '(a)' ) '  For the New York Times issue date:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  NYT1 -> JED1 by historical record.'
  write ( *, '(a)' ) '  NYT1 -> JED2 by "NYT_TO_JED"'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Conversions agree between 1905 and 1995.'
  write ( *, '(a)' ) '  but there are problems at 1905 and earlier.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      Volume    Issue   =>  JED1        Date1'
  write ( *, '(a)' ) '      Volume    Issue  <=   JED2        Date2'
  write ( *, '(a)' ) '                            JED diff'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    issue1 = issue_test(test)
    volume1 = volume_test(test)

    y1 = y_test(test)
    m1 = m_test(test)
    d1 = d_test(test)
    f1 = 0.0D+00

    call ymdf_to_jed_common ( y1, m1, d1, f1, jed1 )
    call ymdf_to_s_common ( y1, m1, d1, f1, s1 )
    call nyt_to_jed ( volume1, issue1, jed2 )
    call jed_to_ymdf_common ( jed2, y2, m2, d2, f2 )
    call ymdf_to_s_common ( y2, m2, d2, f2, s2 )
    call jed_to_nyt ( jed2, volume3, issue3 )
    write ( *, '(a)' ) ' '
    write ( *, '(2x,i8,2x,i8,2x,f11.2,2x,a25)' ) volume1, issue1, jed1, s1
    write ( *, '(2x,i8,2x,i8,2x,f11.2,2x,a25)' ) volume3, issue3, jed2, s2
    write ( *, '(2x,8x,2x,8x,2x,f11.2)' )                         jed1 - jed2

  end do

  return
end
subroutine test0076 ( )

!*****************************************************************************80
!
!! TEST0076 tests JED_TO_RD and RD_TO_JED.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 October 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) jed1
  real ( kind = 8 ) rd2
  real ( kind = 8 ) jed3

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0076'
  write ( *, '(a)' ) '  For the RD:'
  write ( *, '(a)' ) '  JED_TO_RD: JED -> RD.'
  write ( *, '(a)' ) '  RD_TO_JED: RD -> JED.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  JED (in)    RD                JED (out)'
  write ( *, '(a)' ) ' '

  i = 0

  do

    i = i + 1
    call jed_test ( i, jed1 )

    if ( jed1 < 0.0D+00 ) then
      exit
    end if

    call jed_to_rd ( jed1, rd2 )

    call rd_to_jed ( rd2, jed3 )

    write ( *, '(2x,f11.2,5x,f11.2,5x,f11.2)' ) jed1, rd2, jed3

  end do

  return
end
subroutine test0006 ( )

!*****************************************************************************80
!
!! TEST0006 tests JED_TO_SS_UNIX and SS_TO_JED_UNIX.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 October 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  real ( kind = 8 ) jed_epoch
  real ( kind = 8 ) jed1
  real ( kind = 8 ) jed3
  character ( len = 20 ) s2
  real ( kind = 8 ) ss2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0006'
  write ( *, '(a)' ) '  For the UNIX SS calendar:'
  write ( *, '(a)' ) '  JED_TO_SS_UNIX: JED -> SS.'
  write ( *, '(a)' ) '  SS_TO_JED_UNIX: SS -> JED.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  JED (in)     SS                JED (out)'
  write ( *, '(a)' ) ' '

  call epoch_to_jed_unix ( jed_epoch )

  i = 0

  do

    i = i + 1
    call jed_test ( i, jed1 )

    if ( jed1 < 0.0D+00 ) then
      exit
    end if

    if ( jed_epoch <= jed1 ) then

      call jed_to_ss_unix ( jed1, ss2 )

      call r8_to_s_left ( ss2, s2 )

      call ss_to_jed_unix ( ss2, jed3 )

      write ( *, '(2x,f11.2,5x,a,5x,f11.2)' ) jed1, s2, jed3

    end if

  end do

  return
end
subroutine test008 ( )

!*****************************************************************************80
!
!! TEST008 tests JED_TO_YJF_ENGLISH and YJF_TO_JED_ENGLISH.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 October 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) f2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j2
  real ( kind = 8 ) jed_epoch
  real ( kind = 8 ) jed1
  real ( kind = 8 ) jed3
  character ( len = 20 ) s2
  integer ( kind = 4 ) y2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST008'
  write ( *, '(a)' ) '  For the English calendar:'
  write ( *, '(a)' ) '  JED_TO_YJF_ENGLISH: JED -> YJF.'
  write ( *, '(a)' ) '  YJF_TO_JED_ENGLISH: YJF -> JED.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  JED (in)    YJF 	       JED (out)'
  write ( *, '(a)' ) ' '

  call epoch_to_jed_english ( jed_epoch )

  i = 0

  do

    i = i + 1
    call jed_test ( i, jed1 )
    
    if ( jed1 < 0.0D+00 ) then
      exit
    end if

    if ( jed_epoch <= jed1 )then

      call jed_to_yjf_english ( jed1, y2, j2, f2 )

      call yjf_to_s_english ( y2, j2, f2, s2 )

      call yjf_to_jed_english ( y2, j2, f2, jed3 )

      write ( *, '(2x,f11.2,5x,a,5x,f11.2)' ) jed1, s2, jed3

    end if

  end do

  return
end
subroutine test009 ( )

!*****************************************************************************80
!
!! TEST009 tests JED_TO_YJF_GREGORIAN and YJF_TO_JED_GREGORIAN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 October 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) f2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j2
  real ( kind = 8 ) jed_epoch
  real ( kind = 8 ) jed1
  real ( kind = 8 ) jed3
  character ( len = 20 ) s2
  integer ( kind = 4 ) y2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST009'
  write ( *, '(a)' ) '  For the Gregorian calendar:'
  write ( *, '(a)' ) '  JED_TO_YJF_GREGORIAN: JED -> YJF.'
  write ( *, '(a)' ) '  YJF_TO_JED_GREGORIAN: YJF -> JED.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  JED (in)    YJF 	       JED (out)'
  write ( *, '(a)' ) ' '

  call epoch_to_jed_gregorian ( jed_epoch )

  i = 0

  do

    i = i + 1
    call jed_test ( i, jed1 )
    
    if ( jed1 < 0.0D+00 ) then
      exit
    end if

    if ( jed_epoch <= jed1 )then

      call jed_to_yjf_gregorian ( jed1, y2, j2, f2 )

      call yjf_to_s_gregorian ( y2, j2, f2, s2 )

      call yjf_to_jed_gregorian ( y2, j2, f2, jed3 )

      write ( *, '(2x,f11.2,5x,a,5x,f11.2)' ) jed1, s2, jed3

    end if

  end do

  return
end
subroutine test010 ( )

!*****************************************************************************80
!
!! TEST010 tests JED_TO_YJF_HEBREW and YJF_TO_JED_HEBREW.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 October 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) f2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j2
  real ( kind = 8 ) jed_epoch
  real ( kind = 8 ) jed1
  real ( kind = 8 ) jed3
  character ( len = 20 ) s2
  integer ( kind = 4 ) y2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST010'
  write ( *, '(a)' ) '  For the Hebrew calendar:'
  write ( *, '(a)' ) '  JED_TO_YJF_HEBREW: JED -> YJF.'
  write ( *, '(a)' ) '  YJF_TO_JED_HEBREW: YJF -> JED.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  JED (in)    YJF                JED (out)'
  write ( *, '(a)' ) ' '

  call epoch_to_jed_hebrew ( jed_epoch )

  i = 0

  do

    i = i + 1
    call jed_test ( i, jed1 )
    
    if ( jed1 < 0.0D+00 ) then
      exit
    end if

    if ( jed_epoch <= jed1 )then

      call jed_to_yjf_hebrew ( jed1, y2, j2, f2 )

      call yjf_to_s_hebrew ( y2, j2, f2, s2 )

      call yjf_to_jed_hebrew ( y2, j2, f2, jed3 )

      write ( *, '(2x,f11.2,5x,a,5x,f11.2)' ) jed1, s2, jed3

    end if

  end do

  return
end
subroutine test011 ( )

!*****************************************************************************80
!
!! TEST011 tests JED_TO_YJF_REPUBLICAN and YJF_TO_JED_REPUBLICAN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) f2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j2
  real ( kind = 8 ) jed_epoch
  real ( kind = 8 ) jed1
  real ( kind = 8 ) jed3
  character ( len = 20 ) s2
  integer ( kind = 4 ) y2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST011'
  write ( *, '(a)' ) '  For the Republican calendar:'
  write ( *, '(a)' ) '  JED_TO_YJF_REPUBLICAN: JED -> YJF.'
  write ( *, '(a)' ) '  YJF_TO_JED_REPUBLICAN: YJF -> JED.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)') '  JED (in)    YJF                JED (out)'
  write ( *, '(a)' ) ' '

  call epoch_to_jed_republican ( jed_epoch )

  i = 0

  do

    i = i + 1
    call jed_test ( i, jed1 )
    
    if ( jed1 < 0.0D+00 ) then
      exit
    end if

    if ( jed_epoch <= jed1 )then

      call jed_to_yjf_republican ( jed1, y2, j2, f2 )

      call yjf_to_s_republican ( y2, j2, f2, s2 )

      call yjf_to_jed_republican ( y2, j2, f2, jed3 )

      write ( *, '(2x,f11.2,5x,a,5x,f11.2)' ) jed1, s2, jed3

    end if

  end do

  return
end
subroutine test012 ( )

!*****************************************************************************80
!
!! TEST012 tests JED_TO_YJF_ROMAN and YJF_TO_JED_ROMAN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) f2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j2
  real ( kind = 8 ) jed_epoch
  real ( kind = 8 ) jed1
  real ( kind = 8 ) jed3
  character ( len = 20 ) s2
  integer ( kind = 4 ) y2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST012'
  write ( *, '(a)' ) '  For the Roman calendar:'
  write ( *, '(a)' ) '  JED_TO_YJF_ROMAN: JED -> YJF.'
  write ( *, '(a)' ) '  YJF_TO_JED_ROMAN: YJF -> JED.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  JED (in)    YJF                JED (out)'
  write ( *, '(a)' ) ' '

  call epoch_to_jed_roman ( jed_epoch )

  i = 0

  do

    i = i + 1
    call jed_test ( i, jed1 )
    
    if ( jed1 < 0.0D+00 ) then
      exit
    end if

    if ( jed_epoch <= jed1 )then

      call jed_to_yjf_roman ( jed1, y2, j2, f2 )

      call yjf_to_s_roman ( y2, j2, f2, s2 )

      call yjf_to_jed_roman ( y2, j2, f2, jed3 )

      write ( *, '(2x,f11.2,5x,a,5x,f11.2)' ) jed1, s2, jed3

    end if

  end do

  return
end
subroutine test0125 ( )

!*****************************************************************************80
!
!! TEST0125 tests JED_TO_YMDF_ALEXANDRIAN and YMDF_TO_JED_ALEXANDRIAN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d2
  real ( kind = 8 ) f2
  integer ( kind = 4 ) i
  real ( kind = 8 ) jed_epoch
  real ( kind = 8 ) jed1
  real ( kind = 8 ) jed3
  integer ( kind = 4 ) m2
  character ( len = 25 ) s2
  integer ( kind = 4 ) y2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0125'
  write ( *, '(a)' ) '  For the Alexandrian calendar:'
  write ( *, '(a)' ) '  JED_TO_YMDF_ALEXANDRIAN: JED -> YMDF.'
  write ( *, '(a)' ) '  YMDF_TO_JED_ALEXANDRIAN: YMDF -> JED.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  JED (in)    YMDF               JED (out)'
  write ( *, '(a)' ) ' '

  call epoch_to_jed_alexandrian ( jed_epoch )

  i = 0

  do

    i = i + 1
    call jed_test ( i, jed1 )

    if ( jed1 < 0.0D+00 ) then
      exit
    end if

    if ( jed_epoch <= jed1 ) then

      call jed_to_ymdf_alexandrian ( jed1, y2, m2, d2, f2 )

      call ymd_to_s_alexandrian ( y2, m2, d2, s2 )

      call ymdf_to_jed_alexandrian ( y2, m2, d2, f2, jed3 )

      write ( *, '(2x,f11.2,5x,a,5x,f11.2)' ) jed1, s2, jed3

    end if

  end do

  return
end
subroutine test013 ( )

!*****************************************************************************80
!
!! TEST013 tests JED_TO_YMDF_ARMENIAN and YMDF_TO_JED_ARMENIAN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d2
  real ( kind = 8 ) f2
  integer ( kind = 4 ) i
  real ( kind = 8 ) jed_epoch
  real ( kind = 8 ) jed1
  real ( kind = 8 ) jed3
  integer ( kind = 4 ) m2
  character ( len = 20 ) s2
  integer ( kind = 4 ) y2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST013'
  write ( *, '(a)' ) '  For the Armenian calendar:'
  write ( *, '(a)' ) '  JED_TO_YMDF_ARMENIAN: JED -> YMDF.'
  write ( *, '(a)' ) '  YMDF_TO_JED_ARMENIAN: YMDF -> JED.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  JED (in)    YMDF               JED (out)'
  write ( *, '(a)' ) ' '

  call epoch_to_jed_armenian ( jed_epoch )

  i = 0

  do

    i = i + 1
    call jed_test ( i, jed1 )
    
    if ( jed1 < 0.0D+00 ) then
      exit
    end if

    if ( jed_epoch <= jed1 ) then

      call jed_to_ymdf_armenian ( jed1, y2, m2, d2, f2 )

      call ymdf_to_s_numeric ( y2, m2, d2, f2, s2 )

      call ymdf_to_jed_armenian ( y2, m2, d2, f2, jed3 )

      write ( *, '(2x,f11.2,5x,a,5x,f11.2)' ) jed1, s2, jed3

    end if

  end do

  return
end
subroutine test014 ( )

!*****************************************************************************80
!
!! TEST014 tests JED_TO_YMDF_BAHAI and YMDF_TO_JED_BAHAI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d2
  real ( kind = 8 ) f2
  integer ( kind = 4 ) i
  real ( kind = 8 ) jed_epoch
  real ( kind = 8 ) jed1
  real ( kind = 8 ) jed3
  integer ( kind = 4 ) m2
  character ( len = 20 ) s2
  integer ( kind = 4 ) y2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST014'
  write ( *, '(a)' ) '  For the Bahai calendar:'
  write ( *, '(a)' ) '  JED_TO_YMDF_BAHAI: JED -> YMDF.'
  write ( *, '(a)' ) '  YMDF_TO_JED_BAHAI: YMDF -> JED.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  JED (in)    YMDF	       JED (out)'
  write ( *, '(a)' ) ' '

  call epoch_to_jed_bahai ( jed_epoch )

  i = 0

  do

    i = i + 1
    call jed_test ( i, jed1 )
    
    if ( jed1 < 0.0D+00 ) then
      exit
    end if

    if ( jed_epoch <= jed1 ) then

      call jed_to_ymdf_bahai ( jed1, y2, m2, d2, f2 )

      call ymdf_to_s_numeric ( y2, m2, d2, f2, s2 )

      call ymdf_to_jed_bahai ( y2, m2, d2, f2, jed3 )

      write ( *, '(2x,f11.2,5x,a,5x,f11.2)' ) jed1, s2, jed3

    end if

  end do

  return
end
subroutine test015 ( )

!*****************************************************************************80
!
!! TEST015 tests JED_TO_YMDF_COMMON and YMDF_TO_JED_COMMON.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d2
  real ( kind = 8 ) f2
  integer ( kind = 4 ) i
  real ( kind = 8 ) jed1
  real ( kind = 8 ) jed3
  integer ( kind = 4 ) m2
  character ( len = 20 ) s2
  integer ( kind = 4 ) y2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST015'
  write ( *, '(a)' ) '  For the Common calendar:'
  write ( *, '(a)' ) '  JED_TO_YMDF_COMMON: JED -> YMDF.'
  write ( *, '(a)' ) '  YMDF_TO_JED_COMMON: YMDF -> JED'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  JED (in)    YMDF               JED (out)'
  write ( *, '(a)' ) ' '

  i = 0

  do

    i = i + 1
    call jed_test ( i, jed1 )

    if ( jed1 < 0.0D+00 ) then
      exit
    end if

    call jed_to_ymdf_common ( jed1, y2, m2, d2, f2 )

    call ymdf_to_s_common ( y2, m2, d2, f2, s2 )

    call ymdf_to_jed_common ( y2, m2, d2, f2, jed3 )

    write ( *, '(2x,f11.2,5x,a,5x,f11.2)' ) jed1, s2, jed3

  end do

  return
end
subroutine test016 ( )

!*****************************************************************************80
!
!! TEST016 tests JED_TO_YMDF_COPTIC and YMDF_TO_JED_COPTIC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d2
  real ( kind = 8 ) f2
  integer ( kind = 4 ) i
  real ( kind = 8 ) jed_epoch
  real ( kind = 8 ) jed1
  real ( kind = 8 ) jed3
  integer ( kind = 4 ) m2
  character ( len = 20 ) s2
  integer ( kind = 4 ) y2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST016'
  write ( *, '(a)' ) '  For the Coptic calendar:'
  write ( *, '(a)' ) '  JED_TO_YMDF_COPTIC: JED -> YMDF.'
  write ( *, '(a)' ) '  YMDF_TO_JED_COPTIC: YMDF -> JED.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  JED (in)    YMDF               JED (out)'
  write ( *, '(a)' ) ' '

  call epoch_to_jed_coptic ( jed_epoch )

  i = 0

  do

    i = i + 1
    call jed_test ( i, jed1 )

    if ( jed1 < 0.0D+00 ) then
      exit
    end if

    if ( jed_epoch <= jed1 ) then

      call jed_to_ymdf_coptic ( jed1, y2, m2, d2, f2 )

      call ymdf_to_s_numeric ( y2, m2, d2, f2, s2 )

      call ymdf_to_jed_coptic ( y2, m2, d2, f2, jed3 )

      write ( *, '(2x,f11.2,5x,a,5x,f11.2)' ) jed1, s2, jed3

    end if

  end do

  return
end
subroutine test017 ( )

!*****************************************************************************80
!
!! TEST017 tests JED_TO_YMDF_EG_CIVIL and YMDF_TO_JED_EG_CIVIL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d2
  real ( kind = 8 ) f2
  integer ( kind = 4 ) i
  real ( kind = 8 ) jed_epoch
  real ( kind = 8 ) jed1
  real ( kind = 8 ) jed3
  integer ( kind = 4 ) m2
  character ( len = 25 ) s2
  integer ( kind = 4 ) y2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST017'
  write ( *, '(a)' ) '  For the Egyptian Civil calendar:'
  write ( *, '(a)' ) '  JED_TO_YMDF_EG_CIVIL: JED -> YMDF.'
  write ( *, '(a)' ) '  YMDF_TO_JED_EG_CIVIL: YMDF -> JED.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  JED (in)    YMDF               JED (out)'
  write ( *, '(a)' ) ' '

  call epoch_to_jed_eg_civil ( jed_epoch )

  i = 0

  do

    i = i + 1
    call jed_test ( i, jed1 )

    if ( jed1 < 0.0D+00 ) then
      exit
    end if
    
    if ( jed_epoch <= jed1 ) then

      call jed_to_ymdf_eg_civil ( jed1, y2, m2, d2, f2 )

      call ymd_to_s_eg_civil ( y2, m2, d2, s2 )

      call ymdf_to_jed_eg_civil ( y2, m2, d2, f2, jed3 )

      write ( *, '(2x,f11.2,5x,a,5x,f11.2)' ) jed1, s2, jed3

    end if

  end do

  return
end
subroutine test0175 ( )

!*****************************************************************************80
!
!! TEST0175 tests JED_TO_YMDF_EG_LUNAR and YMDF_TO_JED_EG_LUNAR.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d2
  real ( kind = 8 ) f2
  integer ( kind = 4 ) i
  real ( kind = 8 ) jed_epoch
  real ( kind = 8 ) jed1
  real ( kind = 8 ) jed3
  integer ( kind = 4 ) m2
  character ( len = 25 ) s2
  integer ( kind = 4 ) y2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0175'
  write ( *, '(a)' ) '  For the Egyptian Lunar calendar:'
  write ( *, '(a)' ) '  JED_TO_YMDF_EG_LUNAR: JED -> YMDF.'
  write ( *, '(a)' ) '  YMDF_TO_JED_EG_LUNAR: YMDF -> JED.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  JED (in)    YMDF               JED (out)'
  write ( *, '(a)' ) ' '

  call epoch_to_jed_eg_lunar ( jed_epoch )

  i = 0

  do

    i = i + 1
    call jed_test ( i, jed1 )

    if ( jed1 < 0.0D+00 ) then
      exit
    end if
    
    if ( jed_epoch <= jed1 ) then

      call jed_to_ymdf_eg_lunar ( jed1, y2, m2, d2, f2 )

      call ymd_to_s_eg_lunar ( y2, m2, d2, s2 )

      call ymdf_to_jed_eg_lunar ( y2, m2, d2, f2, jed3 )

      write ( *, '(2x,f11.2,5x,a,5x,f11.2)' ) jed1, s2, jed3

    end if

  end do

  return
end
subroutine test018 ( )

!*****************************************************************************80
!
!! TEST018 tests JED_TO_YMDF_ENGLISH and YMDF_TO_JED_ENGLISH.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d2
  real ( kind = 8 ) f2
  integer ( kind = 4 ) i
  real ( kind = 8 ) jed1
  real ( kind = 8 ) jed3
  integer ( kind = 4 ) m2
  character ( len = 20 ) s2
  integer ( kind = 4 ) y2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST018'
  write ( *, '(a)' ) '  For the English calendar,'
  write ( *, '(a)' ) '  JED_TO_YMDF_ENGLISH: JED -> YMDF.'
  write ( *, '(a)' ) '  YMDF_TO_JED_ENGLISH: YMDF -> JED.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  JED (in)    YMDF               JED (out)'
  write ( *, '(a)' ) ' '

  i = 0

  do

    i = i + 1
    call jed_test ( i, jed1 )

    if ( jed1 < 0.0D+00 ) then
      exit
    end if

    call jed_to_ymdf_english ( jed1, y2, m2, d2, f2 )

    call ymdf_to_s_english ( y2, m2, d2, f2, s2 )

    call ymdf_to_jed_english ( y2, m2, d2, f2, jed3 )

    write ( *, '(2x,f11.2,5x,a,5x,f11.2)' ) jed1, s2, jed3

  end do

  return
end
subroutine test019 ( )

!*****************************************************************************80
!
!! TEST019 tests JED_TO_YMDF_ETHIOPIAN and YMDF_TO_JED_ETHIOPIAN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d2
  real ( kind = 8 ) f2
  integer ( kind = 4 ) i
  real ( kind = 8 ) jed_epoch
  real ( kind = 8 ) jed1
  real ( kind = 8 ) jed3
  integer ( kind = 4 ) m2
  character ( len = 20 ) s2
  integer ( kind = 4 ) y2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST019'
  write ( *, '(a)' ) '  For the Ethiopian calendar:'
  write ( *, '(a)' ) '  JED_TO_YMDF_ETHIOPIAN: JED -> YMDF.'
  write ( *, '(a)' ) '  YMDF_TO_JED_ETHIOPIAN: YMDF -> JED.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  JED (in)    YMDF		JED (out)'
  write ( *, '(a)' ) ' '

  call epoch_to_jed_ethiopian ( jed_epoch )

  i = 0

  do

    i = i + 1
    call jed_test ( i, jed1 )
    
    if ( jed1 < 0.0D+00 ) then
      exit
    end if

    if ( jed_epoch <= jed1 ) then

      call jed_to_ymdf_ethiopian ( jed1, y2, m2, d2, f2 )

      call ymdf_to_s_numeric ( y2, m2, d2, f2, s2 )

      call ymdf_to_jed_ethiopian ( y2, m2, d2, f2, jed3 )

      write ( *, '(2x,f11.2,5x,a,5x,f11.2)' ) jed1, s2, jed3

    end if

  end do

  return
end
subroutine test020 ( )

!*****************************************************************************80
!
!! TEST020 tests JED_TO_YMDF_GREGORIAN and YMDF_TO_JED_GREGORIAN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d2
  real ( kind = 8 ) f2
  integer ( kind = 4 ) i
  real ( kind = 8 ) jed1
  real ( kind = 8 ) jed3
  integer ( kind = 4 ) m2
  character ( len = 20 ) s2
  integer ( kind = 4 ) y2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST020'
  write ( *, '(a)' ) '  For the Gregorian calendar:'
  write ( *, '(a)' ) '  JED_TO_YMDF_GREGORIAN: JED -> YMDF.'
  write ( *, '(a)' ) '  YMDF_TO_JED_GREGORIAN: YMDF -> JED.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  JED (in)    YMDF                JED (out)'
  write ( *, '(a)' ) ' '

  i = 0

  do

    i = i + 1
    call jed_test ( i, jed1 )
    
    if ( jed1 < 0.0D+00 ) then
      exit
    end if

    call jed_to_ymdf_gregorian ( jed1, y2, m2, d2, f2 )

    call ymdf_to_s_gregorian ( y2, m2, d2, f2, s2 )

    call ymdf_to_jed_gregorian ( y2, m2, d2, f2, jed3 )

    write ( *, '(2x,f11.2,5x,a,5x,f11.2)' ) jed1, s2, jed3

  end do

  return
end
subroutine test165 ( )

!*****************************************************************************80
!
!! TEST165 tests JED_TO_YMDF_GREGORIAN2 and YMDF_TO_JED_GREGORIAN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d2
  real ( kind = 8 ) f2
  integer ( kind = 4 ) i
  real ( kind = 8 ) jed1
  real ( kind = 8 ) jed3
  integer ( kind = 4 ) m2
  character ( len = 25 ) s2
  integer ( kind = 4 ) y2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST165'
  write ( *, '(a)' ) '  For the Gregorian calendar:'
  write ( *, '(a)' ) '  JED_TO_YMDF_GREGORIAN2: JED -> YMDF.'
  write ( *, '(a)' ) '  YMDF_TO_JED_GREGORIAN: YMDF -> JED.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  JED (in)    YMDF                JED (out)'
  write ( *, '(a)' ) ' '

  i = 0

  do

    i = i + 1
    call jed_test ( i, jed1 )
    
    if ( jed1 < 0.0D+00 ) then
      exit
    end if

    call jed_to_ymdf_gregorian2 ( jed1, y2, m2, d2, f2 )

    call ymdf_to_s_gregorian ( y2, m2, d2, f2, s2 )

    call ymdf_to_jed_gregorian ( y2, m2, d2, f2, jed3 )

    write ( *, '(2x,f11.2,5x,a,5x,f11.2)' ) jed1, s2, jed3

  end do

  return
end
subroutine test17 ( )

!*****************************************************************************80
!
!! TEST17 tests JED_TO_YMDF_HEBREW and YMDF_TO_JED_HEBREW.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d2
  real ( kind = 8 ) f2
  integer ( kind = 4 ) i
  real ( kind = 8 ) jed_epoch
  real ( kind = 8 ) jed1
  real ( kind = 8 ) jed3
  integer ( kind = 4 ) m2
  character ( len = 20 ) s2
  integer ( kind = 4 ) y2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST17'
  write ( *, '(a)' ) '  For the Hebrew calendar:'
  write ( *, '(a)' ) '  JED_TO_YMDF_HEBREW: JED -> YMDF.'
  write ( *, '(a)' ) '  YMDF_TO_JED_HEBREW: YMDF -> JED.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  JED (in)    YMDF                       JED (out)'
  write ( *, '(a)' ) ' '

  call epoch_to_jed_hebrew ( jed_epoch )

  i = 0

  do

    i = i + 1
    call jed_test ( i, jed1 )

    if ( jed1 < 0.0D+00 ) then
      exit
    end if

    if ( jed_epoch <= jed1 ) then

      call jed_to_ymdf_hebrew ( jed1, y2, m2, d2, f2 )

      call ymdf_to_s_hebrew ( y2, m2, d2, f2, s2 )

      call ymdf_to_jed_hebrew ( y2, m2, d2, f2, jed3 )

      write ( *, '(2x,f11.2,5x,a,5x,f11.2)' ) jed1, s2, jed3

    end if

  end do

  return
end
subroutine test175 ( )

!*****************************************************************************80
!
!! TEST175 tests JED_TO_YMDF_HINDU_SOLAR and YMDF_TO_JED_HINDU_SOLAR.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d2
  real ( kind = 8 ) f2
  integer ( kind = 4 ) i
  real ( kind = 8 ) jed_epoch
  real ( kind = 8 ) jed1
  real ( kind = 8 ) jed3
  integer ( kind = 4 ) m2
  character ( len = 20 ) s2
  integer ( kind = 4 ) y2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST175'
  write ( *, '(a)' ) '  For the Hindu Solar calendar:'
  write ( *, '(a)' ) '  JED_TO_YMDF_HINDU_SOLAR: JED -> YMDF.'
  write ( *, '(a)' ) '  YMDF_TO_JED_HINDU_SOLAR: YMDF -> JED.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  JED (in)    YMDF		       JED (out)'
  write ( *, '(a)' ) ' '

  call epoch_to_jed_hindu_solar ( jed_epoch )

  i = 0

  do

    i = i + 1
    call jed_test ( i, jed1 )

    if ( jed1 < 0.0D+00 ) then
      exit
    end if

    if ( jed_epoch <= jed1 ) then

      call jed_to_ymdf_hindu_solar ( jed1, y2, m2, d2, f2 )

      call ymdf_to_s_numeric ( y2, m2, d2, f2, s2 )

      call ymdf_to_jed_hindu_solar ( y2, m2, d2, f2, jed3 )

      write ( *, '(2x,f11.2,5x,a,5x,f11.2)' ) jed1, s2, jed3

    end if

  end do

  return
end
subroutine test18 ( )

!*****************************************************************************80
!
!! TEST18 tests JED_TO_YMDF_ISLAMIC_A and YMDF_TO_JED_ISLAMIC_A.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d2
  real ( kind = 8 ) f2
  integer ( kind = 4 ) i
  real ( kind = 8 ) jed_epoch
  real ( kind = 8 ) jed1
  real ( kind = 8 ) jed3
  integer ( kind = 4 ) m2
  character ( len = 20 ) s2
  integer ( kind = 4 ) y2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST18'
  write ( *, '(a)' ) '  For the Islamic A calendar:'
  write ( *, '(a)' ) '  JED_TO_YMDF_ISLAMIC_A: JED -> YMDF.'
  write ( *, '(a)' ) '  YMDF_TO_JED_ISLAMIC_A: YMDF -> JED.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  JED (in)    YMDF                JED (out)'
  write ( *, '(a)' ) ' '

  call epoch_to_jed_islamic_a ( jed_epoch )

  i = 0

  do

    i = i + 1
    call jed_test ( i, jed1 )

    if ( jed1 < 0.0D+00 ) then
      exit
    end if

    if ( jed_epoch <= jed1 ) then

      call jed_to_ymdf_islamic_a ( jed1, y2, m2, d2, f2 )

      call ymdf_to_s_islamic ( y2, m2, d2, f2, s2 )

      call ymdf_to_jed_islamic_a ( y2, m2, d2, f2, jed3 )

      write ( *, '(2x,f11.2,5x,a,5x,f11.2)' ) jed1, s2, jed3

    end if

  end do

  return
end
subroutine test185 ( )

!*****************************************************************************80
!
!! TEST185 tests JED_TO_YMDF_ISLAMIC_A and YMDF_TO_JED_ISLAMIC_A2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d2
  real ( kind = 8 ) f2
  integer ( kind = 4 ) i
  real ( kind = 8 ) jed_epoch
  real ( kind = 8 ) jed1
  real ( kind = 8 ) jed3
  integer ( kind = 4 ) m2
  character ( len = 20 ) s2
  integer ( kind = 4 ) y2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST185'
  write ( *, '(a)' ) '  For the Islamic A calendar:'
  write ( *, '(a)' ) '  JED_TO_YMDF_ISLAMIC_A: JED -> YMDF.'
  write ( *, '(a)' ) '  YMDF_TO_JED_ISLAMIC_A2: YMDF -> JED.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  JED (in)    YMDF		JED (out)'
  write ( *, '(a)' ) ' '

  call epoch_to_jed_islamic_a ( jed_epoch )

  i = 0

  do

    i = i + 1
    call jed_test ( i, jed1 )

    if ( jed1 < 0.0D+00 ) then
      exit
    end if

    if ( jed_epoch <= jed1 ) then

      call jed_to_ymdf_islamic_a ( jed1, y2, m2, d2, f2 )

      call ymdf_to_s_islamic ( y2, m2, d2, f2, s2 )

      call ymdf_to_jed_islamic_a2 ( y2, m2, d2, f2, jed3 )

      write ( *, '(2x,f11.2,5x,a,5x,f11.2)' ) jed1, s2, jed3

    end if

  end do

  return
end
subroutine test19 ( )

!*****************************************************************************80
!
!! TEST19 tests JED_TO_YMDF_ISLAMIC_B and YMDF_TO_JED_ISLAMIC_B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d2
  real ( kind = 8 ) f2
  integer ( kind = 4 ) i
  real ( kind = 8 ) jed_epoch
  real ( kind = 8 ) jed1
  real ( kind = 8 ) jed3
  integer ( kind = 4 ) m2
  character ( len = 20 ) s2
  integer ( kind = 4 ) y2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST19'
  write ( *, '(a)' ) '  For the Islamic B calendar:'
  write ( *, '(a)' ) '  JED_TO_YMDF_ISLAMIC_B: JED -> YMDF.'
  write ( *, '(a)' ) '  YMDF_TO_JED_ISLAMIC_B: YMDF -> JED.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  JED (in)    YMDF                JED (out)'
  write ( *, '(a)' ) ' '

  call epoch_to_jed_islamic_b ( jed_epoch )

  i = 0

  do

    i = i + 1
    call jed_test ( i, jed1 )

    if ( jed1 < 0.0D+00 ) then
      exit
    end if

    if ( jed_epoch <= jed1 ) then

      call jed_to_ymdf_islamic_b ( jed1, y2, m2, d2, f2 )

      call ymdf_to_s_islamic ( y2, m2, d2, f2, s2 )

      call ymdf_to_jed_islamic_b ( y2, m2, d2, f2, jed3 )

      write ( *, '(2x,f11.2,5x,a,5x,f11.2)' ) jed1, s2, jed3

    end if

  end do

  return
end
subroutine test195 ( )

!*****************************************************************************80
!
!! TEST195 tests JED_TO_YMDF_JELALI and YMDF_TO_JED_JELALI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d2
  real ( kind = 8 ) f2
  integer ( kind = 4 ) i
  real ( kind = 8 ) jed_epoch
  real ( kind = 8 ) jed1
  real ( kind = 8 ) jed3
  integer ( kind = 4 ) m2
  character ( len = 20 ) s2
  integer ( kind = 4 ) y2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST195'
  write ( *, '(a)' ) '  For the Jelali calendar:'
  write ( *, '(a)' ) '  JED_TO_YMDF_JELALI: JED -> YMDF.'
  write ( *, '(a)' ) '  YMDF_TO_JED_JELALI: YMDF -> JED.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  JED (in)    YMDF                JED (out)'
  write ( *, '(a)' ) ' '

  call epoch_to_jed_jelali ( jed_epoch )

  i = 0

  do

    i = i + 1
    call jed_test ( i, jed1 )

    if ( jed1 < 0.0D+00 ) then
      exit
    end if

    if ( jed1 < jed_epoch ) then
      cycle
    end if

    call jed_to_ymdf_jelali ( jed1, y2, m2, d2, f2 )

    call ymdf_to_s_numeric ( y2, m2, d2, f2, s2 )

    call ymdf_to_jed_jelali ( y2, m2, d2, f2, jed3 )

    write ( *, '(2x,f11.2,5x,a,5x,f11.2)' ) jed1, s2, jed3

  end do

  return
end
subroutine test20 ( )

!*****************************************************************************80
!
!! TEST20 tests JED_TO_YMDF_JULIAN and YMDF_TO_JED_JULIAN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d2
  real ( kind = 8 ) f2
  integer ( kind = 4 ) i
  real ( kind = 8 ) jed1
  real ( kind = 8 ) jed3
  integer ( kind = 4 ) m2
  character ( len = 20 ) s2
  integer ( kind = 4 ) y2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST20'
  write ( *, '(a)' ) '  For the Julian calendar:'
  write ( *, '(a)' ) '  JED_TO_YMDF_JULIAN: JED -> YMDF.'
  write ( *, '(a)' ) '  YMDF_TO_JED_JULIAN: YMDF -> JED.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  JED (in)    YMDF                JED (out)'
  write ( *, '(a)' ) ' '

  i = 0

  do

    i = i + 1
    call jed_test ( i, jed1 )

    if ( jed1 < 0.0D+00 ) then
      exit
    end if

    call jed_to_ymdf_julian ( jed1, y2, m2, d2, f2 )

    call ymdf_to_s_julian ( y2, m2, d2, f2, s2 )

    call ymdf_to_jed_julian ( y2, m2, d2, f2, jed3 )

    write ( *, '(2x,f11.2,5x,a,5x,f11.2)' ) jed1, s2, jed3

  end do

  return
end
subroutine test21 ( )

!*****************************************************************************80
!
!! TEST21 tests JED_TO_YMDF_JULIAN2 and YMDF_TO_JED_JULIAN2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d2
  real ( kind = 8 ) f2
  integer ( kind = 4 ) i

  real ( kind = 8 ) jed1
  real ( kind = 8 ) jed3
  integer ( kind = 4 ) m2
  character ( len = 20 ) s2
  integer ( kind = 4 ) y2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST21'
  write ( *, '(a)' ) '  For the Julian calendar:'
  write ( *, '(a)' ) '  JED_TO_YMDF_JULIAN2: JED -> YMDF.'
  write ( *, '(a)' ) '  YMDF_TO_JED_JULIAN2: YMDF -> JED.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  JED (in)    YMDF		JED (out)'
  write ( *, '(a)' ) ' '

  i = 0

  do

    i = i + 1
    call jed_test ( i, jed1 )

    if ( jed1 < 0.0D+00 ) then
      exit
    end if

    call jed_to_ymdf_julian2 ( jed1, y2, m2, d2, f2 )

    call ymdf_to_s_julian ( y2, m2, d2, f2, s2 )

    call ymdf_to_jed_julian2 ( y2, m2, d2, f2, jed3 )

    write ( *, '(2x,f11.2,5x,a,5x,f11.2)' ) jed1, s2, jed3

  end do

  return
end
subroutine test215 ( )

!*****************************************************************************80
!
!! TEST215 tests JED_TO_YMDF_JULIAN3 and YMDF_TO_JED_JULIAN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d2
  real ( kind = 8 ) f2
  integer ( kind = 4 ) i
  real ( kind = 8 ) jed1
  real ( kind = 8 ) jed3
  integer ( kind = 4 ) m2
  character ( len = 20 ) s2
  integer ( kind = 4 ) y2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST215'
  write ( *, '(a)' ) '  For the Julian calendar:'
  write ( *, '(a)' ) '  JED_TO_YMDF_JULIAN3: JED -> YMDF.'
  write ( *, '(a)' ) '  YMDF_TO_JED_JULIAN: YMDF -> JED.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  JED (in)    YMDF		JED (out)'
  write ( *, '(a)' ) ' '

  i = 0

  do

    i = i + 1
    call jed_test ( i, jed1 )

    if ( jed1 < 0.0D+00 ) then
      exit
    end if

    call jed_to_ymdf_julian3 ( jed1, y2, m2, d2, f2 )

    call ymdf_to_s_julian ( y2, m2, d2, f2, s2 )

    call ymdf_to_jed_julian ( y2, m2, d2, f2, jed3 )

    write ( *, '(2x,f11.2,5x,a,5x,f11.2)' ) jed1, s2, jed3

  end do

  return
end
subroutine test22 ( )

!*****************************************************************************80
!
!! TEST22 tests JED_TO_YMDF_KHWARIZMIAN and YMDF_TO_JED_KHWARIZMIAN.
!
!  Modified:
!
!    18 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d2
  real ( kind = 8 ) f2
  integer ( kind = 4 ) i
  real ( kind = 8 ) jed_epoch
  real ( kind = 8 ) jed1
  real ( kind = 8 ) jed3
  integer ( kind = 4 ) m2
  character ( len = 20 ) s2
  integer ( kind = 4 ) y2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST22'
  write ( *, '(a)' ) '  For the Khwarizmian calendar:'
  write ( *, '(a)' ) '  JED_TO_YMDF_KHWARIZMIAN: JED -> YMDF.'
  write ( *, '(a)' ) '  YMDF_TO_JED_KHWARIZMIAN: YMDF -> JED.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  JED (in)    YMDF                JED (out)'
  write ( *, '(a)' ) ' '

  call epoch_to_jed_khwarizmian ( jed_epoch )

  i = 0

  do

    i = i + 1
    call jed_test ( i, jed1 )

    if ( jed1 < 0.0D+00 ) then
      exit
    end if

    if ( jed_epoch <= jed1 ) then

      call jed_to_ymdf_khwarizmian ( jed1, y2, m2, d2, f2 )

      call ymdf_to_s_numeric ( y2, m2, d2, f2, s2 )

      call ymdf_to_jed_khwarizmian ( y2, m2, d2, f2, jed3 )

      write ( *, '(2x,f11.2,5x,a,5x,f11.2)' ) jed1, s2, jed3

    end if

  end do

  return
end
subroutine test23 ( )

!*****************************************************************************80
!
!! TEST23 tests JED_TO_YMDF_MACEDONIAN and YMDF_TO_JED_MACEDONIAN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d2
  real ( kind = 8 ) f2
  integer ( kind = 4 ) i
  real ( kind = 8 ) jed_epoch
  real ( kind = 8 ) jed1
  real ( kind = 8 ) jed3
  integer ( kind = 4 ) m2
  character ( len = 20 ) s2
  integer ( kind = 4 ) y2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST23'
  write ( *, '(a)' ) '  For the Macedonian calendar:'
  write ( *, '(a)' ) '  JED_TO_YMDF_MACEDONIAN: JED -> YMDF.'
  write ( *, '(a)' ) '  YMDF_TO_JED_MACEDONIAN: YMDF -> JED.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  JED (in)    YMDF                JED (out)'
  write ( *, '(a)' ) ' '

  call epoch_to_jed_macedonian ( jed_epoch )

  i = 0

  do

    i = i + 1
    call jed_test ( i, jed1 )

    if ( jed1 < 0.0D+00 ) then
      exit
    end if

    if ( jed_epoch <= jed1 ) then

      call jed_to_ymdf_macedonian ( jed1, y2, m2, d2, f2 )
  
      call ymdf_to_s_numeric ( y2, m2, d2, f2, s2 )

      call ymdf_to_jed_macedonian ( y2, m2, d2, f2, jed3 )

      write ( *, '(2x,f11.2,5x,a,5x,f11.2)' ) jed1, s2, jed3

    end if

  end do

  return
end
subroutine test24 ( )

!*****************************************************************************80
!
!! TEST24 tests JED_TO_YMDF_PERSIAN and YMDF_TO_JED_PERSIAN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d2
  real ( kind = 8 ) f2
  integer ( kind = 4 ) i
  real ( kind = 8 ) jed_epoch
  real ( kind = 8 ) jed1
  real ( kind = 8 ) jed3
  integer ( kind = 4 ) m2
  character ( len = 20 ) s2
  integer ( kind = 4 ) y2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST24'
  write ( *, '(a)' ) '  For the Persian calendar:'
  write ( *, '(a)' ) '  JED_TO_YMDF_PERSIAN: JED -> YMDF.'
  write ( *, '(a)' ) '  YMDF_TO_JED_PERSIAN: YMDF -> JED.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  JED (in)    YMDF                JED (out)'
  write ( *, '(a)' ) ' '

  call epoch_to_jed_persian ( jed_epoch )

  i = 0

  do

    i = i + 1
    call jed_test ( i, jed1 )

    if ( jed1 < 0.0D+00 ) then
      exit
    end if

    if ( jed_epoch <= jed1 ) then

      call jed_to_ymdf_persian ( jed1, y2, m2, d2, f2 )

      call ymdf_to_s_numeric ( y2, m2, d2, f2, s2 )

      call ymdf_to_jed_persian ( y2, m2, d2, f2, jed3 )

      write ( *, '(2x,f11.2,5x,a,5x,f11.2)' ) jed1, s2, jed3

    end if

  end do

  return
end
subroutine test25 ( )

!*****************************************************************************80
!
!! TEST25 tests JED_TO_YMDF_REPUBLICAN and YMDF_TO_JED_REPUBLICAN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d2
  real ( kind = 8 ) f2
  integer ( kind = 4 ) i
  real ( kind = 8 ) jed_epoch
  real ( kind = 8 ) jed1
  real ( kind = 8 ) jed3
  integer ( kind = 4 ) m2
  character ( len = 20 ) s2
  integer ( kind = 4 ) y2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST25'
  write ( *, '(a)' ) '  For the Republican calendar:'
  write ( *, '(a)' ) '  JED_TO_YMDF_REPUBLICAN: JED -> YMDF.'
  write ( *, '(a)' ) '  YMDF_TO_JED_REPUBLICAN: YMDF -> JED.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  JED (in)    YMDF                JED (out)'
  write ( *, '(a)' ) ' '

  call epoch_to_jed_republican ( jed_epoch )

  i = 0

  do

    i = i + 1
    call jed_test ( i, jed1 )

    if ( jed1 < 0.0D+00 ) then
      exit
    end if

    if ( jed_epoch <= jed1 ) then

      call jed_to_ymdf_republican ( jed1, y2, m2, d2, f2 )

      call ymdf_to_s_republican ( y2, m2, d2, f2, s2 )

      call ymdf_to_jed_republican ( y2, m2, d2, f2, jed3 )

      write ( *, '(2x,f11.2,5x,a,5x,f11.2)' ) jed1, s2, jed3

    end if

  end do

  return
end
subroutine test255 ( )

!*****************************************************************************80
!
!! TEST255 tests JED_TO_YMDF_ROMAN and YMDF_TO_JED_ROMAN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d2
  real ( kind = 8 ) f2
  integer ( kind = 4 ) i
  real ( kind = 8 ) jed_epoch
  real ( kind = 8 ) jed1
  real ( kind = 8 ) jed3
  integer ( kind = 4 ) m2
  character ( len = 45 ) s2
  integer ( kind = 4 ) y2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST255'
  write ( *, '(a)' ) '  For the Roman calendar:'
  write ( *, '(a)' ) '  JED_TO_YMDF_ROMAN: JED -> YMDF.'
  write ( *, '(a)' ) '  YMDF_TO_JED_ROMAN: YMDF -> JED.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  JED (in)    YMDF                JED (out)'
  write ( *, '(a)' ) ' '

  call epoch_to_jed_roman ( jed_epoch )

  i = 0

  do

    i = i + 1
    call jed_test ( i, jed1 )

    if ( jed1 < 0.0D+00 ) then
      exit
    end if

    if ( jed_epoch <= jed1 ) then

      call jed_to_ymdf_roman ( jed1, y2, m2, d2, f2 )

      call ymdf_to_s_roman ( y2, m2, d2, f2, s2 )
 
      call ymdf_to_jed_roman ( y2, m2, d2, f2, jed3 )

      write ( *, '(2x,f11.2,5x,a,5x,f11.2)' ) jed1, s2, jed3

    end if

  end do

  return
end
subroutine test26 ( )

!*****************************************************************************80
!
!! TEST26 tests JED_TO_YMDF_SAKA and YMDF_TO_JED_SAKA.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d2
  real ( kind = 8 ) f2
  integer ( kind = 4 ) i
  real ( kind = 8 ) jed_epoch
  real ( kind = 8 ) jed1
  real ( kind = 8 ) jed3
  integer ( kind = 4 ) m2
  character ( len = 20 ) s2
  integer ( kind = 4 ) y2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST26'
  write ( *, '(a)' ) '  For the Saka calendar:'
  write ( *, '(a)' ) '  JED_TO_YMDF_SAKA: JED -> YMDF.'
  write ( *, '(a)' ) '  YMDF_TO_JED_SAKA: YMDF -> JED.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  JED (in)    YMDF                JED (out)'
  write ( *, '(a)' ) ' '

  call epoch_to_jed_saka ( jed_epoch )

  i = 0

  do

    i = i + 1
    call jed_test ( i, jed1 )

    if ( jed1 < 0.0D+00 ) then
      exit
    end if

    if ( jed_epoch <= jed1 ) then

      call jed_to_ymdf_saka ( jed1, y2, m2, d2, f2 )

      call ymdf_to_s_numeric ( y2, m2, d2, f2, s2 )

      call ymdf_to_jed_saka ( y2, m2, d2, f2, jed3 )

      write ( *, '(2x,f11.2,5x,a,5x,f11.2)' ) jed1, s2, jed3

    end if

  end do

  return
end
subroutine test265 ( )

!*****************************************************************************80
!
!! TEST265 tests JED_TO_YMDF_SOOR_SAN and YMDF_TO_JED_SOOR_SAN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d2
  real ( kind = 8 ) f2
  integer ( kind = 4 ) i
  real ( kind = 8 ) jed_epoch
  real ( kind = 8 ) jed1
  real ( kind = 8 ) jed3
  integer ( kind = 4 ) m2
  character ( len = 20 ) s2
  integer ( kind = 4 ) y2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST265'
  write ( *, '(a)' ) '  For the Soor San calendar:'
  write ( *, '(a)' ) '  JED_TO_YMDF_SOOR_SAN: JED -> YMDF.'
  write ( *, '(a)' ) '  YMDF_TO_JED_SOOR_SAN: YMDF -> JED.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  JED (in)    YMDF		JED (out)'
  write ( *, '(a)' ) ' '

  call epoch_to_jed_soor_san ( jed_epoch )

  i = 0

  do

    i = i + 1
    call jed_test ( i, jed1 )

    if ( jed1 < 0.0D+00 ) then
      exit
    end if

    if ( jed_epoch <= jed1 ) then

      call jed_to_ymdf_soor_san ( jed1, y2, m2, d2, f2 )

      call ymdf_to_s_numeric ( y2, m2, d2, f2, s2 )

      call ymdf_to_jed_soor_san ( y2, m2, d2, f2, jed3 )

      write ( *, '(2x,f11.2,5x,a,5x,f11.2)' ) jed1, s2, jed3

    end if

  end do

  return
end
subroutine test27 ( )

!*****************************************************************************80
!
!! TEST27 tests JED_TO_YMDF_SYRIAN and YMDF_TO_JED_SYRIAN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d2
  real ( kind = 8 ) f2
  integer ( kind = 4 ) i
  real ( kind = 8 ) jed_epoch
  real ( kind = 8 ) jed1
  real ( kind = 8 ) jed3
  integer ( kind = 4 ) m2
  character ( len = 20 ) s2
  integer ( kind = 4 ) y2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST27'
  write ( *, '(a)' ) '  For the Syrian calendar:'
  write ( *, '(a)' ) '  JED_TO_YMDF_SYRIAN: JED -> YMDF.'
  write ( *, '(a)' ) '  YMDF_TO_JED_SYRIAN: YMDF -> JED.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  JED (in)    YMDF                JED (out)'
  write ( *, '(a)' ) ' '

  call epoch_to_jed_syrian ( jed_epoch )

  i = 0

  do

    i = i + 1
    call jed_test ( i, jed1 )

    if ( jed1 < 0.0D+00 ) then
      exit
    end if

    if ( jed_epoch <= jed1 ) then

      call jed_to_ymdf_syrian ( jed1, y2, m2, d2, f2 )

      call ymdf_to_s_numeric ( y2, m2, d2, f2, s2 )

      call ymdf_to_jed_syrian ( y2, m2, d2, f2, jed3 )

      write ( *, '(2x,f11.2,5x,a,5x,f11.2)' ) jed1, s2, jed3

    end if

  end do

  return
end
subroutine test275 ( )

!*****************************************************************************80
!
!! TEST275 tests JED_TO_YMDF_ZOROASTRIAN and YMDF_TO_JED_ZOROASTRIAN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d2
  real ( kind = 8 ) f2
  integer ( kind = 4 ) i
  real ( kind = 8 ) jed_epoch
  real ( kind = 8 ) jed1
  real ( kind = 8 ) jed3
  integer ( kind = 4 ) m2
  character ( len = 20 ) s2
  integer ( kind = 4 ) y2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST275'
  write ( *, '(a)' ) '  For the Zoroastrian calendar:'
  write ( *, '(a)' ) '  JED_TO_YMDF_ZOROASTRIAN: JED -> YMDF.'
  write ( *, '(a)' ) '  YMDF_TO_JED_ZOROASTRIAN: YMDF -> JED.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)') '  JED (in)    YMDF                JED (out)'
  write ( *, '(a)' ) ' '

  call epoch_to_jed_zoroastrian ( jed_epoch )

  i = 0

  do
    i = i + 1
    call jed_test ( i, jed1 )

    if ( jed1 < 0.0D+00 ) then
      exit
    end if

    if ( jed_epoch <= jed1 ) then

      call jed_to_ymdf_zoroastrian ( jed1, y2, m2, d2, f2 )

      call ymdf_to_s_numeric ( y2, m2, d2, f2, s2 )

      call ymdf_to_jed_zoroastrian ( y2, m2, d2, f2, jed3 )

      write ( *, '(2x,f11.2,5x,a,5x,f11.2)' ) jed1, s2, jed3

    end if

  end do

  return
end
subroutine test28 ( )

!*****************************************************************************80
!
!! TEST28 tests MONTH_CAL_COMMON.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d
  real ( kind = 8 ) f
  real ( kind = 8 ) jed
  integer ( kind = 4 ) m
  integer ( kind = 4 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST28'
  write ( *, '(a)' ) '  For the Common calendar:'
  write ( *, '(a)' ) '  MONTH_CAL_COMMON prints a month calendar.'
  write ( *, '(a)' ) ' '
 
  y = 1582
  m = 10

  call month_cal_common ( y, m )

  y = 1752
  m = 9

  call month_cal_common ( y, m )
 
  call now_to_jed ( jed )
  call jed_to_ymdf_common ( jed, y, m, d, f )
  call month_cal_common ( y, m )

  return
end
subroutine test29 ( )

!*****************************************************************************80
!
!! TEST29 tests MONTH_CAL_ENGLISH.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d
  real ( kind = 8 ) f
  real ( kind = 8 ) jed
  integer ( kind = 4 ) m
  integer ( kind = 4 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST29'
  write ( *, '(a)' ) '  For the English calendar:'
  write ( *, '(a)' ) '  MONTH_CAL_ENGLISH prints a month calendar.'
  write ( *, '(a)' ) ' '
 
  y = 1582
  m = 10

  call month_cal_english ( y, m )

  y = 1752
  m = 9

  call month_cal_english ( y, m )
 
  call now_to_jed ( jed )
  call jed_to_ymdf_english ( jed, y, m, d, f )
  call month_cal_english ( y, m )

  return
end
subroutine test30 ( )

!*****************************************************************************80
!
!! TEST30 tests MONTH_CAL_GREGORIAN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d
  real ( kind = 8 ) f
  real ( kind = 8 ) jed
  integer ( kind = 4 ) m
  integer ( kind = 4 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST30'
  write ( *, '(a)' ) '  For the Gregorian calendar:'
  write ( *, '(a)' ) '  MONTH_CAL_GREGORIAN prints a month calendar.'
  write ( *, '(a)' ) ' '
 
  y = 1582
  m = 10

  call month_cal_gregorian ( y, m )

  y = 1752
  m = 9

  call month_cal_gregorian ( y, m )
 
  call now_to_jed ( jed )
  call jed_to_ymdf_gregorian ( jed, y, m, d, f )
  call month_cal_gregorian ( y, m )

  return
end
subroutine test31 ( )

!*****************************************************************************80
!
!! TEST31 tests MONTH_CAL_HEBREW.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) d2
  real ( kind = 8 ) f
  real ( kind = 8 ) f2
  real ( kind = 8 ) jed
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST31'
  write ( *, '(a)' ) '  For the Hebrew calendar:'
  write ( *, '(a)' ) '  MONTH_CAL_HEBREW prints a month calendar.'
  write ( *, '(a)' ) ' '
 
  y = 1582
  m = 10
  d = 1
  f = 0.5D+00

  call ymdf_to_jed_common ( y, m, d, f, jed )
  call jed_to_ymdf_hebrew ( jed, y2, m2, d2, f2 )
  call month_cal_hebrew ( y2, m2 )

  y = 1752
  m = 9
  d = 1
  f = 0.5D+00

  call ymdf_to_jed_common ( y, m, d, f, jed )
  call jed_to_ymdf_hebrew ( jed, y2, m2, d2, f2 )
  call month_cal_hebrew ( y2, m2 )
 
  call now_to_jed ( jed )
  call jed_to_ymdf_hebrew ( jed, y2, m2, d2, f2 )
  call month_cal_hebrew ( y2, m2 )

  return
end
subroutine test315 ( )

!*****************************************************************************80
!
!! TEST315 tests MONTH_CAL_ISLAMIC_A.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d
  real ( kind = 8 ) f
  real ( kind = 8 ) jed
  integer ( kind = 4 ) m
  integer ( kind = 4 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST315'
  write ( *, '(a)' ) '  For the Islamic A calendar:'
  write ( *, '(a)' ) '  MONTH_CAL_ISLAMIC_A prints a month calendar.'
  write ( *, '(a)' ) ' '
 
  y = 500
  m = 1

  call month_cal_islamic_a ( y, m )

  y = 500
  m = 2

  call month_cal_islamic_a ( y, m )
 
  call now_to_jed ( jed )
  call jed_to_ymdf_islamic_a ( jed, y, m, d, f )
  call month_cal_islamic_a ( y, m )

  return
end
subroutine test32 ( )

!*****************************************************************************80
!
!! TEST32 tests MONTH_CAL_JULIAN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d
  real ( kind = 8 ) f
  real ( kind = 8 ) jed
  integer ( kind = 4 ) m
  integer ( kind = 4 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST32'
  write ( *, '(a)' ) '  For the Julian calendar:'
  write ( *, '(a)' ) '  MONTH_CAL_JULIAN prints a month calendar.'
  write ( *, '(a)' ) ' '
 
  y = 1582
  m = 10

  call month_cal_julian ( y, m )

  y = 1752
  m = 9

  call month_cal_julian ( y, m )
 
  call now_to_jed ( jed )
  call jed_to_ymdf_julian ( jed, y, m, d, f )
  call month_cal_julian ( y, m )

  return
end
subroutine test325 ( )

!*****************************************************************************80
!
!! TEST325 tests MONTH_CAL_REPUBLICAN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d
  real ( kind = 8 ) f
  real ( kind = 8 ) jed
  integer ( kind = 4 ) m
  integer ( kind = 4 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST325'
  write ( *, '(a)' ) '  For the Republican calendar:'
  write ( *, '(a)' ) '  MONTH_CAL_REPUBLICAN prints a month calendar.'
  write ( *, '(a)' ) ' '
 
  y = 3
  m = 12

  call month_cal_republican ( y, m )

  y = 3
  m = 13

  call month_cal_republican ( y, m )
 
  call now_to_jed ( jed )
  call jed_to_ymdf_republican ( jed, y, m, d, f )
  call month_cal_republican ( y, m )

  return
end
subroutine test326 ( )

!*****************************************************************************80
!
!! TEST326 tests MONTH_CAL_ROMAN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST326'
  write ( *, '(a)' ) '  For the Roman calendar:'
  write ( *, '(a)' ) '  MONTH_CAL_ROMAN prints a month calendar.'
 
  y = 100
  m = 12

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Year = ', y
  write ( *, '(a,i6)' ) '  Month = ', m

  call month_cal_roman ( y, m )

  y = 256
  m = 2

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Year = ', y
  write ( *, '(a,i6)' ) '  Month = ', m

  call month_cal_roman ( y, m )
 
  return
end
subroutine test327 ( )

!*****************************************************************************80
!
!! TEST327 tests MONTH_CAL_STORE_COMMON.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  character ( len = 20 ) lines(6)
  integer ( kind = 4 ) m
  integer ( kind = 4 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST327'
  write ( *, '(a)' ) '  For the Common calendar:'
  write ( *, '(a)' ) '  MONTH_CAL_STORE_COMMON writes the day numbers for'
  write ( *, '(a)' ) '  a monthly calendar into a data structure.'

  y = 1984

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Calendar:'
  write ( *, '(a,i6)' ) '  Year = ', y

  do m = 1, 12

    call month_cal_store_common ( y, m, lines )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) '  Month = ', m
    write ( *, '(a)' ) ' '
    do i = 1, 6
      write ( *, '(2x,i1,4x,a)' ) i, lines(i)
    end do

  end do

  return
end
subroutine test328 ( )

!*****************************************************************************80
!
!! TEST328 tests MONTH_LENGTH_BAHAI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_test = 1

  integer ( kind = 4 ) days
  integer ( kind = 4 ) i_test
  integer ( kind = 4 ) m
  integer ( kind = 4 ) month_length_bahai
  character ( len = 15 ) month_name
  integer ( kind = 4 ) months
  character ( len = 15 ) sy
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y_test(n_test)
  integer ( kind = 4 ) year_length_bahai
  integer ( kind = 4 ) year_length_months_bahai

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST328'
  write ( *, '(a)' ) '  For the Bahai calendar:'
  write ( *, '(a)' ) '  MONTH_LENGTH_BAHAI returns month lengths.'

  y_test(1) = 60

  do i_test = 1, n_test

    y = y_test(i_test)
    call y_to_s_bahai ( y, sy )
    months = year_length_months_bahai ( y )
    days = year_length_bahai ( y )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,i6)' ) y
    write ( *, '(2x,a)' ) trim ( sy )
    write ( *, '(a,i6)' ) '    Year length in months = ', months
    write ( *, '(a,i6)' ) '    Year length in days = ', days
    write ( *, '(a)' ) ' '

    do m = 1, months
      call month_to_month_name_bahai ( m, month_name )
      write ( *, '(6x,a,2x,i4)' ) month_name, month_length_bahai ( y, m )
    end do

  end do
 
  return
end
subroutine test33 ( )

!*****************************************************************************80
!
!! TEST33 tests MONTH_LENGTH_COMMON.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_test = 4

  integer ( kind = 4 ) days
  integer ( kind = 4 ) i_test
  integer ( kind = 4 ) m
  integer ( kind = 4 ) month_length_common
  character ( len = 10 ) month_name
  integer ( kind = 4 ) months
  character ( len = 15 ) sy
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y_test(n_test)
  integer ( kind = 4 ) year_length_common
  integer ( kind = 4 ) year_length_months_common

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST33'
  write ( *, '(a)' ) '  For the Common calendar:'
  write ( *, '(a)' ) '  MONTH_LENGTH_COMMON returns month lengths.'

  y_test(1) = 1582
  y_test(2) = 1752
  y_test(3) = 1900
  y_test(4) = 2000

  do i_test = 1, n_test

    y = y_test(i_test)
    call y_to_s_common ( y, sy )
    months = year_length_months_common ( y )
    days = year_length_common ( y )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,i6)' ) y
    write ( *, '(2x,a)' ) trim ( sy )
    write ( *, '(a,i6)' ) '    Year length in months = ', months
    write ( *, '(a,i6)' ) '    Year length in days = ', days
    write ( *, '(a)' ) ' '

    do m = 1, months
      call month_to_month_name_common ( m, month_name )
      write ( *, '(6x,a,2x,i4)' ) month_name, month_length_common ( y, m )
    end do

  end do
 
  return
end
subroutine test335 ( )

!*****************************************************************************80
!
!! TEST335 tests MONTH_LENGTH_COPTIC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_test = 2

  integer ( kind = 4 ) days
  integer ( kind = 4 ) i_test
  integer ( kind = 4 ) m
  integer ( kind = 4 ) month_length_coptic
  character ( len = 15 ) month_name
  integer ( kind = 4 ) months
  character ( len = 15 ) sy
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y_test(n_test)
  integer ( kind = 4 ) year_length_coptic
  integer ( kind = 4 ) year_length_months_coptic

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST335'
  write ( *, '(a)' ) '  For the Coptic calendar,'
  write ( *, '(a)' ) '  MONTH_LENGTH_COPTIC returns month lengths.'

  y_test(1) = 3
  y_test(2) = 4

  do i_test = 1, n_test

    y = y_test(i_test)
    call y_to_s_coptic ( y, sy )
    months = year_length_months_coptic ( y )
    days = year_length_coptic ( y )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,i6)' ) y
    write ( *, '(2x,a)' ) trim ( sy )
    write ( *, '(a,i6)' ) '    Year length in months = ', months
    write ( *, '(a,i6)' ) '    Year length in days = ', days
    write ( *, '(a)' ) ' '

    do m = 1, months
      call month_to_month_name_coptic ( m, month_name )
      write ( *, '(6x,a,2x,i4)' ) month_name, month_length_coptic ( y, m )
    end do

  end do
 
  return
end
subroutine test336 ( )

!*****************************************************************************80
!
!! TEST336 tests MONTH_LENGTH_EG_CIVIL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_test = 2

  integer ( kind = 4 ) days
  integer ( kind = 4 ) i_test
  integer ( kind = 4 ) m
  integer ( kind = 4 ) month_length_eg_civil
  character ( len = 15 ) month_name
  integer ( kind = 4 ) months
  character ( len = 15 ) sy
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y_test(n_test)
  integer ( kind = 4 ) year_length_eg_civil
  integer ( kind = 4 ) year_length_months_eg_civil

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST336'
  write ( *, '(a)' ) '  For the Egyptian Civil calendar,'
  write ( *, '(a)' ) '  MONTH_LENGTH_EG_CIVIL returns month lengths.'

  y_test(1) = 3
  y_test(2) = 4

  do i_test = 1, n_test

    y = y_test(i_test)
    call y_to_s_eg_civil ( y, sy )
    months = year_length_months_eg_civil ( y )
    days = year_length_eg_civil ( y )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,i6)' ) y
    write ( *, '(2x,a)' ) trim ( sy )
    write ( *, '(a,i6)' ) '    Year length in months = ', months
    write ( *, '(a,i6)' ) '    Year length in days = ', days
    write ( *, '(a)' ) ' '

    do m = 1, months
      call month_to_month_name_eg_civil ( m, month_name )
      write ( *, '(6x,a,2x,i4)' ) month_name, month_length_eg_civil ( y, m )
    end do

  end do
 
  return
end
subroutine test337 ( )

!*****************************************************************************80
!
!! TEST337 tests MONTH_EG_LUNAR.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_test = 2

  integer ( kind = 4 ) days
  integer ( kind = 4 ) i_test
  integer ( kind = 4 ) m
  integer ( kind = 4 ) month_length_eg_lunar
  character ( len = 15 ) month_name
  integer ( kind = 4 ) months
  character ( len = 15 ) sy
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y_test(n_test)
  integer ( kind = 4 ) year_length_eg_lunar
  integer ( kind = 4 ) year_length_months_eg_lunar

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST337'
  write ( *, '(a)' ) '  For the Egyptian Lunar calendar,'
  write ( *, '(a)' ) '  MONTH_LENGTH_EG_LUNAR returns month lengths.'

  y_test(1) = 1
  y_test(2) = 2

  do i_test = 1, n_test

    y = y_test(i_test)
    call y_to_s_eg_lunar ( y, sy )
    months = year_length_months_eg_lunar ( y )
    days = year_length_eg_lunar ( y )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,i6)' ) y
    write ( *, '(2x,a)' ) trim ( sy )
    write ( *, '(a,i6)' ) '    Year length in months = ', months
    write ( *, '(a,i6)' ) '    Year length in days = ', days
    write ( *, '(a)' ) ' '

    do m = 1, months
      call month_to_month_name_eg_lunar ( m, month_name )
      write ( *, '(6x,a,2x,i4)' ) month_name, month_length_eg_lunar ( y, m )
    end do

  end do
 
  return
end
subroutine test34 ( )

!*****************************************************************************80
!
!! TEST34 tests MONTH_LENGTH_ENGLISH.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_test = 4

  integer ( kind = 4 ) days
  integer ( kind = 4 ) i_test
  integer ( kind = 4 ) m
  integer ( kind = 4 ) month_length_english
  character ( len = 10 ) month_name
  integer ( kind = 4 ) months
  character ( len = 15 ) sy
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y_test(n_test)
  integer ( kind = 4 ) year_length_english
  integer ( kind = 4 ) year_length_months_english

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST34'
  write ( *, '(a)' ) '  For the English calendar:'
  write ( *, '(a)' ) '  MONTH_LENGTH_ENGLISH returns month lengths.'

  y_test(1) = 1582
  y_test(2) = 1752
  y_test(3) = 1900
  y_test(4) = 2000

  do i_test = 1, n_test

    y = y_test(i_test)
    call y_to_s_english ( y, sy )
    months = year_length_months_english ( y )
    days = year_length_english ( y )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,i6)' ) y
    write ( *, '(2x,a)' ) trim ( sy )
    write ( *, '(a,i6)' ) '    Year length in months = ', months
    write ( *, '(a,i6)' ) '    Year length in days = ', days
    write ( *, '(a)' ) ' '

    do m = 1, months
      call month_to_month_name_common ( m, month_name )
      write ( *, '(6x,a,2x,i4)' ) month_name, month_length_english ( y, m )
    end do

  end do

  return
end
subroutine test344 ( )

!*****************************************************************************80
!
!! TEST344 tests MONTH_LENGTH_ETHIOPIAN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_test = 2

  integer ( kind = 4 ) days
  integer ( kind = 4 ) i_test
  integer ( kind = 4 ) m
  integer ( kind = 4 ) month_length_ethiopian
  character ( len = 15 ) month_name
  integer ( kind = 4 ) months
  character ( len = 15 ) sy
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y_test(n_test)
  integer ( kind = 4 ) year_length_ethiopian
  integer ( kind = 4 ) year_length_months_ethiopian

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST344'
  write ( *, '(a)' ) '  For the Ethiopian calendar,'
  write ( *, '(a)' ) '  MONTH_LENGTH_ETHIOPIAN returns month lengths.'

  y_test(1) = 3
  y_test(2) = 4

  do i_test = 1, n_test

    y = y_test(i_test)
    call y_to_s_ethiopian ( y, sy )
    months = year_length_months_ethiopian ( y )
    days = year_length_ethiopian ( y )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,i6)' ) y
    write ( *, '(2x,a)' ) trim ( sy )
    write ( *, '(a,i6)' ) '    Year length in months = ', months
    write ( *, '(a,i6)' ) '    Year length in days = ', days
    write ( *, '(a)' ) ' '

    do m = 1, months
      call month_to_month_name_ethiopian ( m, month_name )
      write ( *, '(6x,a,2x,i4)' ) month_name, month_length_ethiopian ( y, m )
    end do

  end do
 
  return
end
subroutine test345 ( )

!*****************************************************************************80
!
!! TEST345 tests MONTH_LENGTH_GREEK.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_test = 2

  integer ( kind = 4 ) days
  integer ( kind = 4 ) i_test
  integer ( kind = 4 ) m
  integer ( kind = 4 ) month_length_greek
  character ( len = 15 ) month_name
  integer ( kind = 4 ) months
  character ( len = 10 ) sy
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y_test(n_test)
  integer ( kind = 4 ) year_length_greek
  integer ( kind = 4 ) year_length_months_greek

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST345'
  write ( *, '(a)' ) '  For the Greek calendar,'
  write ( *, '(a)' ) '  MONTH_LENGTH_GREEK returns month lengths.'

  y_test(1) = 3
  y_test(2) = 4

  do i_test = 1, n_test

    y = y_test(i_test)
    call y_to_s_greek ( y, sy )
    months = year_length_months_greek ( y )
    days = year_length_greek ( y )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,i6)' ) y
    write ( *, '(2x,a)' ) trim ( sy )
    write ( *, '(a,i6)' ) '    Year length in months = ', months
    write ( *, '(a,i6)' ) '    Year length in days = ', days
    write ( *, '(a)' ) ' '

    do m = 1, months
      call month_to_month_name_greek ( y, m, month_name )
      write ( *, '(6x,a,2x,i4)' ) month_name, month_length_greek ( y, m )
    end do

  end do
 
  return
end
subroutine test35 ( )

!*****************************************************************************80
!
!! TEST35 tests MONTH_LENGTH_GREGORIAN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_test = 4

  integer ( kind = 4 ) days
  integer ( kind = 4 ) i_test
  integer ( kind = 4 ) m
  integer ( kind = 4 ) month_length_gregorian
  character ( len = 10 ) month_name
  integer ( kind = 4 ) months
  character ( len = 10 ) sy
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y_test(n_test)
  integer ( kind = 4 ) year_length_gregorian
  integer ( kind = 4 ) year_length_months_gregorian

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST35'
  write ( *, '(a)' ) '  For the Gregorian calendar:'
  write ( *, '(a)' ) '  MONTH_LENGTH_GREGORIAN returns month lengths.'

  y_test(1) = 1582
  y_test(2) = 1752
  y_test(3) = 1900
  y_test(4) = 2000

  do i_test = 1, n_test

    y = y_test(i_test)
    call y_to_s_gregorian ( y, sy )
    months = year_length_months_gregorian ( y )
    days = year_length_gregorian ( y )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,i6)' ) y
    write ( *, '(2x,a)' ) trim ( sy )
    write ( *, '(a,i6)' ) '    Year length in months = ', months
    write ( *, '(a,i6)' ) '    Year length in days = ', days
    write ( *, '(a)' ) ' '

    do m = 1, months
      call month_to_month_name_common ( m, month_name )
      write ( *, '(6x,a,2x,i4)' ) month_name, month_length_gregorian ( y, m )
    end do

  end do

  return
end
subroutine test36 ( )

!*****************************************************************************80
!
!! TEST36 tests MONTH_LENGTH_HEBREW.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_test = 3

  integer ( kind = 4 ) days
  integer ( kind = 4 ) i_test
  integer ( kind = 4 ) m
  integer ( kind = 4 ) month_length_hebrew
  character ( len = 10 ) month_name
  integer ( kind = 4 ) months
  character ( len = 10 ) sy
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y_test(n_test)
  integer ( kind = 4 ) year_length_hebrew
  integer ( kind = 4 ) year_length_months_hebrew

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST36'
  write ( *, '(a)' ) '  For the Hebrew calendar,'
  write ( *, '(a)' ) '  MONTH_LENGTH_HEBREW returns month lengths.'

  y_test(1) = 5760
  y_test(2) = 5762
  y_test(3) = 5765

  do i_test = 1, n_test

    y = y_test(i_test)
    call y_to_s_hebrew ( y, sy )
    months = year_length_months_hebrew ( y )
    days = year_length_hebrew ( y )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,i6)' ) y
    write ( *, '(2x,a)' ) trim ( sy )
    write ( *, '(a,i6)' ) '    Year length in months = ', months
    write ( *, '(a,i6)' ) '    Year length in days = ', days
    write ( *, '(a)' ) ' '

    do m = 1, months
      call month_to_month_name_hebrew ( y, m, month_name )
      write ( *, '(6x,a,2x,i4)' ) month_name, month_length_hebrew ( y, m )
    end do

  end do
 
  return
end
subroutine test365 ( )

!*****************************************************************************80
!
!! TEST365 tests MONTH_LENGTH_ISLAMIC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_test = 3

  integer ( kind = 4 ) days
  integer ( kind = 4 ) i_test
  integer ( kind = 4 ) m
  integer ( kind = 4 ) month_length_islamic
  character ( len = 10 ) month_name
  integer ( kind = 4 ) months
  character ( len = 10 ) sy
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y_test(n_test)
  integer ( kind = 4 ) year_length_islamic
  integer ( kind = 4 ) year_length_months_islamic

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST365'
  write ( *, '(a)' ) '  For the Islamic calendar,'
  write ( *, '(a)' ) '  MONTH_LENGTH_ISLAMIC returns month lengths.'

  y_test(1) = 500
  y_test(2) = 501
  y_test(3) = 502

  do i_test = 1, n_test

    y = y_test(i_test)
    call y_to_s_islamic ( y, sy )
    months = year_length_months_islamic ( y )
    days = year_length_islamic ( y )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,i6)' ) y
    write ( *, '(2x,a)' ) trim ( sy )
    write ( *, '(a,i6)' ) '    Year length in months = ', months
    write ( *, '(a,i6)' ) '    Year length in days = ', days
    write ( *, '(a)' ) ' '

    do m = 1, months
      call month_to_month_name_islamic ( m, month_name )
      write ( *, '(6x,a,2x,i4)' ) month_name, month_length_islamic ( y, m )
    end do

  end do
 
  return
end
subroutine test37 ( )

!*****************************************************************************80
!
!! TEST37 tests MONTH_LENGTH_JULIAN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_test = 4

  integer ( kind = 4 ) days
  integer ( kind = 4 ) i_test
  integer ( kind = 4 ) m
  integer ( kind = 4 ) month_length_julian
  character ( len = 10 ) month_name
  integer ( kind = 4 ) months
  character ( len = 10 ) sy
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y_test(n_test)
  integer ( kind = 4 ) year_length_julian
  integer ( kind = 4 ) year_length_months_julian

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST37'
  write ( *, '(a)' ) '  For the Julian calendar,'
  write ( *, '(a)' ) '  MONTH_LENGTH_JULIAN returns month lengths.'

  y_test(1) = 1582
  y_test(2) = 1752
  y_test(3) = 1900
  y_test(4) = 2000

  do i_test = 1, n_test

    y = y_test(i_test)
    call y_to_s_julian ( y, sy )
    months = year_length_months_julian ( y )
    days = year_length_julian ( y )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,i6)' ) y
    write ( *, '(2x,a)' ) trim ( sy )
    write ( *, '(a,i6)' ) '    Year length in months = ', months
    write ( *, '(a,i6)' ) '    Year length in days = ', days
    write ( *, '(a)' ) ' '

    do m = 1, months
      call month_to_month_name_common ( m, month_name )
      write ( *, '(6x,a,2x,i4)' ) month_name, month_length_julian ( y, m )
    end do

  end do

  return
end
subroutine test373 ( )

!*****************************************************************************80
!
!! TEST373 tests MONTH_LENGTH_PERSIAN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_test = 2

  integer ( kind = 4 ) days
  integer ( kind = 4 ) i_test
  integer ( kind = 4 ) m
  integer ( kind = 4 ) month_length_persian
  character ( len = 15 ) month_name
  integer ( kind = 4 ) months
  character ( len = 15 ) sy
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y_test(n_test)
  integer ( kind = 4 ) year_length_persian
  integer ( kind = 4 ) year_length_months_persian

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST373'
  write ( *, '(a)' ) '  For the Persian calendar,'
  write ( *, '(a)' ) '  MONTH_LENGTH_PERSIAN returns month lengths.'

  y_test(1) = 3
  y_test(2) = 4

  do i_test = 1, n_test

    y = y_test(i_test)
    call y_to_s_persian ( y, sy )
    months = year_length_months_persian ( y )
    days = year_length_persian ( y )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,i6)' ) y
    write ( *, '(2x,a)' ) trim ( sy )
    write ( *, '(a,i6)' ) '    Year length in months = ', months
    write ( *, '(a,i6)' ) '    Year length in days = ', days
    write ( *, '(a)' ) ' '

    do m = 1, months
      call month_to_month_name_persian ( m, month_name )
      write ( *, '(6x,a,2x,i4)' ) month_name, month_length_persian ( y, m )
    end do

  end do
 
  return
end
subroutine test375 ( )

!*****************************************************************************80
!
!! TEST375 tests MONTH_LENGTH_REPUBLICAN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_test = 1

  integer ( kind = 4 ) days
  integer ( kind = 4 ) i_test
  integer ( kind = 4 ) m
  integer ( kind = 4 ) month_length_republican
  character ( len = 15 ) month_name
  integer ( kind = 4 ) months
  character ( len = 10 ) sy
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y_test(n_test)
  integer ( kind = 4 ) year_length_months_republican
  integer ( kind = 4 ) year_length_republican

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST375'
  write ( *, '(a)' ) '  For the Republican calendar:'
  write ( *, '(a)' ) '  MONTH_LENGTH_REPUBLICAN returns month lengths.'

  y_test(1) = 4

  do i_test = 1, n_test

    y = y_test(i_test)
    call y_to_s_republican ( y, sy )
    months = year_length_months_republican ( y )
    days = year_length_republican ( y )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,i6)' ) y
    write ( *, '(2x,a)' ) trim ( sy )
    write ( *, '(a,i6)' ) '    Year length in months = ', months
    write ( *, '(a,i6)' ) '    Year length in days = ', days
    write ( *, '(a)' ) ' '

    do m = 1, months
      call month_to_month_name_republican ( m, month_name )
      write ( *, '(6x,a,2x,i4)' ) month_name, month_length_republican ( y, m )
    end do

  end do

  return
end
subroutine test376 ( )

!*****************************************************************************80
!
!! TEST376 tests MONTH_LENGTH_ROMAN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_test = 2

  integer ( kind = 4 ) days
  integer ( kind = 4 ) i_test
  integer ( kind = 4 ) m
  integer ( kind = 4 ) month_length_roman
  character ( len = 15 ) month_name
  integer ( kind = 4 ) months
  character ( len = 10 ) sy
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y_test(n_test)
  integer ( kind = 4 ) year_length_months_roman
  integer ( kind = 4 ) year_length_roman

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST376'
  write ( *, '(a)' ) '  For the Roman calendar,'
  write ( *, '(a)' ) '  MONTH_LENGTH_ROMAN returns month lengths.'

  y_test(1) = 3
  y_test(2) = 4

  do i_test = 1, n_test

    y = y_test(i_test)
    call y_to_s_roman ( y, sy )
    months = year_length_months_roman ( y )
    days = year_length_roman ( y )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,i6)' ) y
    write ( *, '(2x,a)' ) trim ( sy )
    write ( *, '(a,i6)' ) '    Year length in months = ', months
    write ( *, '(a,i6)' ) '    Year length in days = ', days
    write ( *, '(a)' ) ' '

    do m = 1, months
      call month_to_month_name_roman ( m, month_name )
      write ( *, '(6x,a,2x,i4)' ) month_name, month_length_roman ( y, m )
    end do
 
  end do

  return
end
subroutine test38 ( )

!*****************************************************************************80
!
!! TEST38 tests MONTH_NAME_TO_MONTH_COMMON
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: ntest = 9

  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  character ( len = 10 ) month_name
  character ( len = 10 ) test(ntest)

  test(1) = 'J'
  test(2) = 'Febooary'
  test(3) = 'Dec.'
  test(4) = 'April'
  test(5) = 'Aug'
  test(6) = 'Mar'
  test(7) = 'May'
  test(8) = 'o'
  test(9) = 'nO'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST38'
  write ( *, '(a)' ) '  For the Common calendar,'
  write ( *, '(a)' ) '  MONTH_NAME_TO_MONTH_COMMON identifies month names:'
  write ( *, '(a)' ) ' '

  do i = 1, ntest
    call month_name_to_month_common ( test(i), m )
    call month_to_month_name_common ( m, month_name )
    write ( *, '(2x,a3,2x,i2,2x,a9)' ) test(i), m, month_name
  end do
 
  return
end
subroutine test389 ( )

!*****************************************************************************80
!
!! TEST389 tests MONTH_TO_MONTH_NAME_BAHAI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) m
  character ( len = 15 ) month_name
  integer ( kind = 4 ) months
  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length_months_bahai

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST389'
  write ( *, '(a)' ) '  For the Bahai calendar,'
  write ( *, '(a)' ) '  MONTH_TO_MONTH_NAME_BAHAI names the months.'
  write ( *, '(a)' ) ' '

  y = 1
  months = year_length_months_bahai ( y )

  do m = 1, months
    call month_to_month_name_bahai ( m, month_name )
    write ( *, '(2x,i2,2x,a)' ) m, month_name
  end do

  return
end
subroutine test39 ( )

!*****************************************************************************80
!
!! TEST39 tests MONTH_TO_MONTH_NAME_COMMON.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) m
  character ( len = 10 ) month_name
  integer ( kind = 4 ) months
  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length_months_common

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST39'
  write ( *, '(a)' ) '  For the Common calendar,'
  write ( *, '(a)' ) '  MONTH_TO_MONTH_NAME_COMMON names the months:'
  write ( *, '(a)' ) ' '

  y = 1
  months = year_length_months_common ( y )

  do m = 1, months
    call month_to_month_name_common ( m, month_name )
    write ( *, '(2x,i2,2x,a)' ) m, month_name
  end do

  return
end
subroutine test394 ( )

!*****************************************************************************80
!
!! TEST394 tests MONTH_TO_MONTH_NAME_EG_CIVIL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) m
  character ( len = 15 ) month_name
  integer ( kind = 4 ) months
  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length_months_eg_civil

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST394'
  write ( *, '(a)' ) '  For the Egyptian Civil calendar,'
  write ( *, '(a)' ) '  MONTH_TO_MONTH_NAME_EG_CIVIL names the months.'
  write ( *, '(a)' ) ' '

  y = 1
  months = year_length_months_eg_civil ( y )

  do m = 1, months
    call month_to_month_name_eg_civil ( m, month_name )
    write ( *, '(2x,i2,2x,a)' ) m, month_name
  end do

  return
end
subroutine test395 ( )

!*****************************************************************************80
!
!! TEST395 tests MONTH_TO_MONTH_NAME_GREEK.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  character ( len = 15 ) month_name
  integer ( kind = 4 ) months
  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length_months_greek

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST395'
  write ( *, '(a)' ) '  For the Greek calendar,'
  write ( *, '(a)' ) '  MONTH_TO_MONTH_NAME_GREEK names the months.'
  write ( *, '(a)' ) ' '

  y = 1
  months = year_length_months_greek ( y )

  do i = 1, months
    m = i
    call month_to_month_name_greek ( y, m, month_name )
    write ( *, '(2x,i2,2x,a)' ) m, month_name
  end do
 
  return
end
subroutine test40 ( )

!*****************************************************************************80
!
!! TEST40 tests MONTH_TO_MONTH_NAME_HEBREW.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  character ( len = 10 ) month_name
  integer ( kind = 4 ) months
  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length_months_hebrew

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST40'
  write ( *, '(a)' ) '  For the Hebrew calendar,'
  write ( *, '(a)' ) '  MONTH_TO_MONTH_NAME_HEBREW names the months.'
  write ( *, '(a)' ) ' '

  y = 1
  months = year_length_months_hebrew ( y )

  do i = 1, months
    m = i
    call month_to_month_name_hebrew ( y, m, month_name )
    write ( *, '(2x,i2,2x,a)' ) m, month_name
  end do
 
  return
end
subroutine test41 ( )

!*****************************************************************************80
!
!! TEST41 tests MONTH_TO_MONTH_NAME_HINDU_LUNAR.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  character ( len = 10 ) month_name
  integer ( kind = 4 ) months
  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length_months_hindu_lunar

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST41'
  write ( *, '(a)' ) '  For the Hindu lunar calendar,'
  write ( *, '(a)' ) '  MONTH_TO_MONTH_NAME_HINDU_LUNAR names the months.'
  write ( *, '(a)' ) ' '

  y = 1
  months = year_length_months_hindu_lunar ( y )

  do i = 1, months
    m = i
    call month_to_month_name_hindu_lunar ( m, month_name )
    write ( *, '(2x,i2,2x,a)' ) m, month_name
  end do
 
  return
end
subroutine test415 ( )

!*****************************************************************************80
!
!! TEST415 tests MONTH_TO_MONTH_NAME_HINDU_SOLAR.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  character ( len = 10 ) month_name
  integer ( kind = 4 ) months
  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length_months_hindu_solar

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST415'
  write ( *, '(a)' ) '  For the Hindu solar calendar,'
  write ( *, '(a)' ) '  MONTH_TO_MONTH_NAME_HINDU_SOLAR names the months.'
  write ( *, '(a)' ) ' '

  y = 1
  months = year_length_months_hindu_solar ( y )

  do i = 1, months
    m = i
    call month_to_month_name_hindu_solar ( m, month_name )
    write ( *, '(2x,i2,2x,a)' ) m, month_name
  end do
 
  return
end
subroutine test42 ( )

!*****************************************************************************80
!
!! TEST42 tests MONTH_TO_MONTH_NAME_ISLAMIC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  character ( len = 10 ) month_name
  integer ( kind = 4 ) months
  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length_months_islamic

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST42'
  write ( *, '(a)' ) '  For the Islamic calendar,'
  write ( *, '(a)' ) '  MONTH_TO_MONTH_NAME_ISLAMIC names the months:'
  write ( *, '(a)' ) ' '

  y = 1
  months = year_length_months_islamic ( y )

  do i = 1, months
    m = i
    call month_to_month_name_islamic ( m, month_name )
    write ( *, '(i4,2x,a)' ) m, month_name
  end do
 
  return
end
subroutine test43 ( )

!*****************************************************************************80
!
!! TEST43 tests MONTH_TO_MONTH_NAME_REPUBLICAN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) months
  character ( len = 15 ) month_name
  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length_months_republican

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST43'
  write ( *, '(a)' ) '  For the Republican calendar,'
  write ( *, '(a)' ) '  MONTH_TO_MONTH_NAME_REPUBLICAN names the months.'
  write ( *, '(a)' ) ' '

  y = 1
  months = year_length_months_republican ( y )

  do m = 1, months
    call month_to_month_name_republican ( m, month_name )
    write ( *, '(i4,2x,a)' ) m, month_name
  end do

  return
end
subroutine test435 ( )

!*****************************************************************************80
!
!! TEST435 tests MONTH_TO_MONTH_NAME_ROMAN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) months
  character ( len = 10 ) month_name
  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length_months_roman

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST435'
  write ( *, '(a)' ) '  For the Roman calendar,'
  write ( *, '(a)' ) '  MONTH_TO_MONTH_NAME_ROMAN names the months.'
  write ( *, '(a)' ) ' '

  y = 1
  months = year_length_months_roman ( y )

  do m = 1, months
    call month_to_month_name_roman ( m, month_name )
    write ( *, '(i4,2x,a)' ) m, month_name
  end do

  return
end
subroutine test44 ( )

!*****************************************************************************80
!
!! TEST44 tests MOON_PHASE_TO_JED.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d
  real ( kind = 8 ) f
  integer ( kind = 4 ) h
  real ( kind = 8 ) jed
  integer ( kind = 4 ) m
  integer ( kind = 4 ) min
  integer ( kind = 4 ) nphase
  integer ( kind = 4 ) phase
  integer ( kind = 4 ) s
  character ( len = 22 ) string
  integer ( kind = 4 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST44'
  write ( *, '(a)' ) '  MOON_PHASE_TO_JED reports the JED on'
  write ( *, '(a)' ) '  which a phase of the moon occurs.'
  write ( *, '(a)' ) ' '

  phase = 2

  write ( *, '(a)' ) '   N JED        YMDHMS date'
  write ( *, '(a)' ) ' '

  do nphase = 1, 10
    call moon_phase_to_jed ( nphase, phase, jed )
    call jed_to_ymdf_common ( jed, y, m, d, f )
    call frac_to_hms ( f, h, min, s )
    call ymdhms_to_s_common ( y, m, d, h, min, s, string )
    write ( *, '(2x,i3,f11.2,3x,a)' ) nphase, jed, string
  end do

  return
end
subroutine test445 ( )

!*****************************************************************************80
!
!! TEST445 tests MOTHERS_DAY.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) m
  character ( len = 40 ) s
  integer ( kind = 4 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST445'
  write ( *, '(a)' ) '  For a given year of the Common calendar,'
  write ( *, '(a)' ) '  compute the day and month of Mother''s Day (US).'
  write ( *, '(a)' ) ' '

  do y = 1995, 2010

    call mothers_day ( y, m, d )
    call ymd_to_s_common ( y, m, d, s )
    write ( *, '(2x,a)' ) trim ( s )

  end do

  return
end
subroutine test45 ( )

!*****************************************************************************80
!
!! TEST45 tests NEW_YEAR_TO_JED_HEBREW.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d3
  real ( kind = 8 ) f3
  integer ( kind = 4 ) i
  real ( kind = 8 ) jed2
  integer ( kind = 4 ) m3
  character ( len = 10 ) s1
  character ( len = 20 ) s3
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y3

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST45'
  write ( *, '(a)' ) '  For the Hebrew calendar,'
  write ( *, '(a)' ) '  NEW_YEAR_TO_JED_HEBREW determines the JED of'
  write ( *, '(a)' ) '    the first day of a year.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  YEAR       JED         YMDF'
  write ( *, '(a)' ) '  Hebrew                 Common'
  write ( *, '(a)' ) ' '

  do i = 0, 20
    y1 = 5760 + i
    call y_to_s_hebrew ( y1, s1 )
    call new_year_to_jed_hebrew ( y1, jed2 )
    call jed_to_ymdf_common ( jed2, y3, m3, d3, f3 )
    call ymdf_to_s_common ( y3, m3, d3, f3, s3 )
    write ( *, '(2x,a,2x,f11.2,5x,a)' ) trim ( s1 ), jed2, s3
  end do

  return
end
subroutine test46 ( )

!*****************************************************************************80
!
!! TEST46 tests NOW_TO_JED, NOW_TO_YJF_COMMON, NOW_TO_YMDF_COMMON, NOW_TO_YMDHMS_COMMON.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d
  character ( len = 8 ) date
  real ( kind = 8 ) f
  integer ( kind = 4 ) h
  real ( kind = 8 ) jed
  integer ( kind = 4 ) j
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  character ( len = 30 ) s
  integer ( kind = 4 ) second
  character ( len = 10 ) time
  integer ( kind = 4 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST46'
  write ( *, '(a)' ) '  For the current time and date, "NOW", '
  write ( *, '(a)' ) '  NOW_TO_JED returns the JED,'
  write ( *, '(a)' ) '  NOW_TO_YJF_COMMON the YJF date,'
  write ( *, '(a)' ) '  NOW_TO_YMDF_COMMON returns the YMDF date,'
  write ( *, '(a)' ) '  NOW_TO_YMDHMS_COMMON the YMDHMS date.'
  write ( *, '(a)' ) ' '

  call date_and_time ( date, time )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  FORTRAN90 DATE_AND_TIME routine says:'
  write ( *, '(a)' ) '    Now is ' // trim ( date ) // ' ' // trim ( time )

  call now_to_jed ( jed )

  write ( *, '(a)' ) ' '
  write ( *, '(a,f11.2)' ) '  NOW_TO_JED_COMMON:    Now is: ', jed

  call now_to_yjf_common ( y, j, f )
  call yjf_to_s_common ( y, j, f, s )

  write ( *, '(a)' ) '  NOW_TO_YJF_COMMON:     Now is: ' // trim ( s )

  call now_to_ymdf_common ( y, m, d, f )
  call ymdf_to_s_common ( y, m, d, f, s )

  write ( *, '(a)' ) '  NOW_TO_YMDF_COMMON:    Now is: ' // trim ( s )

  call now_to_ymdhms_common ( y, m, d, h, n, second )
  call ymdhms_to_s_common ( y, m, d, h, n, second, s )

  write ( *, '(a)' ) '  NOW_TO_YMDHMS_COMMON: Now is: ' // trim ( s )

  return
end
subroutine test47 ( )

!*****************************************************************************80
!
!! TEST47 tests S_TO_HMS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) h2
  integer ( kind = 4 ) m2
  character ( len = 15 ) p
  character ( len = 15 ) s1
  integer ( kind = 4 ) s2
  character ( len = 8 ) s3

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST47'
  write ( *, '(a)' ) '  S_TO_HMS converts a string to an HMS date.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  ------S--------  ------P--------  HH:MM:SS'
  write ( *, '(a)' ) ' '

  s1 = '10:30:44'
  p = 'hh:mm:ss'

  call s_to_hms ( s1, p, h2, m2, s2 )

  call hms_to_s ( h2, m2, s2, s3 )

  write ( *, '(2x,a,2x,a,2x,a)' ) s1, p, s3

  s1 = '10 past 9'
  p = 'mm xxxx h'

  call s_to_hms ( s1, p, h2, m2, s2 )

  call hms_to_s ( h2, m2, s2, s3 )

  write ( *, '(2x,a,2x,a,2x,a)' ) s1, p, s3

  return
end
subroutine test48 ( )

!*****************************************************************************80
!
!! TEST48 tests S_TO_YMD_COMMON.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_test = 5

  integer ( kind = 4 ) d
  integer ( kind = 4 ) i_test
  integer ( kind = 4 ) m
  character ( len = 35 ) p
  character ( len = 35 ) p_test(n_test)
  character ( len = 35 ) s
  character ( len = 35 ) s_test(n_test)
  integer ( kind = 4 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST48'
  write ( *, '(a)' ) '  S_TO_YMD_COMMON converts a string to a YMD date.'
  write ( *, '(a)' ) ' '
 
  s_test(1) = '1999-10-31'
  p_test(1) = 'yyyy-mm-dd'
  s_test(2) = '01/04/2004'
  p_test(2) = 'dd/mm/yyyy'
  s_test(3) = '8/8/88'
  p_test(3) = 'd/m/yy'
  s_test(4) = '4 7'
  p_test(4) = 'd m'
  s_test(5) = 'On day 1 of month 3 of year 1945'
  p_test(5) = 'xx xxx d xx xxxxx m xx xxxx yyyy'

  do i_test = 1, n_test

    s = s_test(i_test)
    p = p_test(i_test)

    call s_to_ymd_common ( s, p, y, m, d )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,a)' ) s
    write ( *, '(2x,a)' ) p
    write ( *, '(a)' ) ' '
    write ( *, '(2x,3i6)' ) y, m, d

  end do

  return
end
subroutine test49 ( )

!*****************************************************************************80
!
!! TEST49 tests S_TO_YMDHMS_COMMON.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_test = 2

  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) i_test
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  character ( len = 35 ) p
  character ( len = 35 ) p_test(n_test)
  character ( len = 35 ) s
  character ( len = 35 ) s_test(n_test)
  integer ( kind = 4 ) second
  integer ( kind = 4 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST49'
  write ( *, '(a)' ) '  S_TO_YMDHMS_COMMON converts a string to a YMDHMS date.'
  write ( *, '(a)' ) ' '
 
  s_test(1) = '1999-10-31-14-59-47'
  p_test(1) = 'YYYY-MM-DD-hh-mm-ss'
  s_test(2) = '8:30, 01 April 2004'
  p_test(2) = 'h:mm, DD NNNNN YYYY'

  do i_test = 1, n_test

    s = s_test(i_test)
    p = p_test(i_test)

    call s_to_ymdhms_common ( s, p, y, m, d, h, n, second )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,a)' ) s
    write ( *, '(2x,a)' ) p
    write ( *, '(a)' ) ' '
    write ( *, '(2x,6i6)' ) y, m, d, h, n, second

  end do

  return
end
subroutine test492 ( )

!*****************************************************************************80
!
!! TEST492 tests THANKSGIVING_CANADA.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 September 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) m
  character ( len = 40 ) s
  integer ( kind = 4 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST492'
  write ( *, '(a)' ) '  THANKSGIVING_CANADA returns, for a given year of '
  write ( *, '(a)' ) '  the Common calendar, the day and month of '
  write ( *, '(a)' ) '  Thanksgiving in Canada.'
  write ( *, '(a)' ) ' '

  do y = 1995, 2020

    call thanksgiving_canada ( y, m, d )
    call ymd_to_s_common ( y, m, d, s )
    write ( *, '(2x,a)' ) trim ( s )

  end do

  return
end
subroutine test493 ( )

!*****************************************************************************80
!
!! TEST493 tests THANKSGIVING_US.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 September 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) m
  character ( len = 40 ) s
  integer ( kind = 4 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST493'
  write ( *, '(a)' ) '  THANKSGIVING_US returns, for a given year of '
  write ( *, '(a)' ) '  the Common calendar, the day and month of '
  write ( *, '(a)' ) '  Thanksgiving (US).'
  write ( *, '(a)' ) ' '

  do y = 1995, 2020

    call thanksgiving_us ( y, m, d )
    call ymd_to_s_common ( y, m, d, s )
    write ( *, '(2x,a)' ) trim ( s )

  end do

  return
end
subroutine test495 ( )

!*****************************************************************************80
!
!! TEST495 tests WEEKDAY_TO_NAME_BAHAI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  character ( len = 15 ) s

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST495'
  write ( *, '(a)' ) '  For the Bahai calendar:'
  write ( *, '(a)' ) '  WEEKDAY_TO_NAME_BAHAI names the days of the week.'
  write ( *, '(a)' ) ' '

  do i = 1, 7
    call weekday_to_name_bahai ( i, s )
    write ( *, '(2x,i2,2x,a,2x,a)' ) i, s
  end do

  return
end
subroutine test50 ( )

!*****************************************************************************80
!
!! TEST50 tests WEEKDAY_TO_NAME_COMMON, WEEKDAY_TO_NAME_COMMON2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  character ( len = 10 ) s1
  character ( len = 2 ) s2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST50'
  write ( *, '(a)' ) '  For the Common calendar:'
  write ( *, '(a)' ) '  WEEKDAY_TO_NAME_COMMON names the days of the week,'
  write ( *, '(a)' ) '  WEEKDAY_TO_NAME_COMMON2 abbreviates the days of the week.'
  write ( *, '(a)' ) ' '

  do i = 1, 7
    call weekday_to_name_common ( i, s1 )
    call weekday_to_name_common2 ( i, s2 )
    write ( *, '(2x,i4,2x,a,2x,a)' ) i, s1, s2
  end do

  return
end
subroutine test501 ( )

!*****************************************************************************80
!
!! TEST501 tests WEEKDAY_TO_NAME_GERMAN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  character ( len = 15 ) sname

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST501'
  write ( *, '(a)' ) '  For the German calendar,'
  write ( *, '(a)' ) '  WEEKDAY_TO_NAME_GERMAN names the days of the week.'
  write ( *, '(a)' ) ' '

  do i = 1, 7
    call weekday_to_name_german ( i, sname )
    write ( *, '(2x,i2,2x,a)' ) i, sname
  end do

  return
end
subroutine test502 ( )

!*****************************************************************************80
!
!! TEST502 tests WEEKDAY_TO_NAME_HEBREW.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  character ( len = 15 ) sname

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST502'
  write ( *, '(a)' ) '  For the Hebrew calendar,'
  write ( *, '(a)' ) '  WEEKDAY_TO_NAME_HEBREW names the days of the week.'
  write ( *, '(a)' ) ' '

  do i = 1, 7
    call weekday_to_name_hebrew ( i, sname )
    write ( *, '(2x,i2,2x,a)' ) i, sname
  end do

  return
end
subroutine test503 ( )

!*****************************************************************************80
!
!! TEST503 tests WEEKDAY_TO_NAME_ISLAMIC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  character ( len = 15 ) sname

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST503'
  write ( *, '(a)' ) '  For the Islamic calendar,'
  write ( *, '(a)' ) '  WEEKDAY_TO_NAME_ISLAMIC names the days of the week.'
  write ( *, '(a)' ) ' '

  do i = 1, 7
    call weekday_to_name_islamic ( i, sname )
    write ( *, '(2x,i2,2x,a)' ) i, sname
  end do

  return
end
subroutine test51 ( )

!*****************************************************************************80
!
!! TEST51 tests WEEKDAY_TO_NAME_REPUBLICAN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  character ( len = 15 ) sname

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST51'
  write ( *, '(a)' ) '  For the Republican calendar,'
  write ( *, '(a)' ) '  WEEKDAY_TO_NAME_REPUBLICAN names the days of the week.'
  write ( *, '(a)' ) ' '

  do i = 1, 10
    call weekday_to_name_republican ( i, sname )
    write ( *, '(2x,i2,2x,a)' ) i, sname
  end do

  return
end
subroutine test515 ( )

!*****************************************************************************80
!
!! TEST515 tests WEEKDAY_TO_NAME_ROMAN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  character ( len = 15 ) sname

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST515'
  write ( *, '(a)' ) '  WEEKDAY_TO_NAME_ROMAN names the days of '
  write ( *, '(a)' ) '  the week in the Roman calendar.'
  write ( *, '(a)' ) ' '

  do i = 1, 7
    call weekday_to_name_roman ( i, sname )
    write ( *, '(2x,i2,2x,a)' ) i, sname
  end do

  return
end
subroutine test5153 ( )

!*****************************************************************************80
!
!! TEST5153 tests YEAR_CAL_COMMON.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST5153'
  write ( *, '(a)' ) '  For the Common calendar:'
  write ( *, '(a)' ) '  YEAR_CAL_COMMON prints a year calendar.'
  write ( *, '(a)' ) ' '
 
  y = 1968
  call year_cal_common ( y )

  return
end
subroutine test51535 ( )

!*****************************************************************************80
!
!! TEST5154 tests YEAR_IS_EMBOLISMIC_EG_LUNAR.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 10 ) sy
  integer ( kind = 4 ) y
  logical year_is_embolismic_eg_lunar

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST51535'
  write ( *, '(a)' ) '  For the Egyptian Lunar calendar:'
  write ( *, '(a)' ) '  YEAR_IS_EMBOLISMIC_EG_LUNAR determines if a year is'
  write ( *, '(a)' ) '    an embolismic year.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Year  Embolismic?'
  write ( *, '(a)' ) ' '

  do y = 1, 25
    call y_to_s_eg_lunar ( y, sy )
    write ( *, '(4x,a,2x,l1)' ) sy, year_is_embolismic_eg_lunar ( y )
  end do

  return
end
subroutine test5154 ( )

!*****************************************************************************80
!
!! TEST5154 tests YEAR_IS_EMBOLISMIC_GREEK.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 10 ) sy
  integer ( kind = 4 ) y
  logical year_is_embolismic_greek

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST5154'
  write ( *, '(a)' ) '  For the Greek calendar:'
  write ( *, '(a)' ) '  YEAR_IS_EMBOLISMIC_GREEK determines if a year is'
  write ( *, '(a)' ) '    an embolismic year.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Year  Embolismic?'
  write ( *, '(a)' ) ' '

  do y = 1, 20
    call y_to_s_greek ( y, sy )
    write ( *, '(4x,a,2x,l1)' ) sy, year_is_embolismic_greek ( y )
  end do

  return
end
subroutine test5155 ( )

!*****************************************************************************80
!
!! TEST5155 tests YEAR_IS_EMBOLISMIC_HEBREW.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 10 ) sy
  integer ( kind = 4 ) y
  logical year_is_embolismic_hebrew

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST5155'
  write ( *, '(a)' ) '  For the Hebrew calendar:'
  write ( *, '(a)' ) '  YEAR_IS_EMBOLISMIC_HEBREW determines if a Hebrew year is'
  write ( *, '(a)' ) '    an embolismic year.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Year  Embolismic?'
  write ( *, '(a)' ) ' '

  do y = 1, 20
    call y_to_s_hebrew ( y, sy )
    write ( *, '(4x,a,2x,l1)' ) trim ( sy ), year_is_embolismic_hebrew ( y )
  end do

  return
end
subroutine test5156 ( )

!*****************************************************************************80
!
!! TEST5156 tests YEAR_IS_LEAP_BAHAI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 10 ) sy
  integer ( kind = 4 ) y
  logical year_is_leap_bahai

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST5156'
  write ( *, '(a)' ) '  For the Bahai calendar:'
  write ( *, '(a)' ) '  YEAR_IS_LEAP_BAHAI reports leap years.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Year  Leap?'
  write ( *, '(a)' ) ' '

  do y = 1990, 2000
    call y_to_s_bahai ( y, sy )
    write ( *, '(2x,a,2x,l1)' ) sy, year_is_leap_bahai ( y )
  end do

  return
end
subroutine test52 ( )

!*****************************************************************************80
!
!! TEST52 tests YEAR_IS_LEAP_COMMON.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 10 ) sy
  integer ( kind = 4 ) y
  logical year_is_leap_common

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST52'
  write ( *, '(a)' ) '  For the Common calendar:'
  write ( *, '(a)' ) '  YEAR_IS_LEAP_COMMON reports leap years.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Year  Leap?'
  write ( *, '(a)' ) ' '

  do y = 1990, 2000
    call y_to_s_common ( y, sy )
    write ( *, '(2x,a,2x,l1)' ) sy, year_is_leap_common ( y )
  end do

  return
end
subroutine test525 ( )

!*****************************************************************************80
!
!! TEST525 tests YEAR_IS_LEAP_EG_LUNAR.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 10 ) sy
  integer ( kind = 4 ) y
  logical year_is_leap_eg_lunar

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST525'
  write ( *, '(a)' ) '  For the Egyptian Lunar calendar:'
  write ( *, '(a)' ) '  YEAR_IS_LEAP_EG_LUNAR reports leap years.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Year  Leap?'
  write ( *, '(a)' ) ' '

  do y = 1, 10
    call y_to_s_eg_lunar ( y, sy )
    write ( *, '(2x,a,2x,l1)' ) sy, year_is_leap_eg_lunar ( y )
  end do

  return
end
subroutine test53 ( )

!*****************************************************************************80
!
!! TEST53 tests YEAR_IS_LEAP_ENGLISH.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 10 ) sy
  integer ( kind = 4 ) y
  logical year_is_leap_english

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST53'
  write ( *, '(a)' ) '  For the English calendar:'
  write ( *, '(a)' ) '  YEAR_IS_LEAP_ENGLISH reports leap years.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Year  Leap?'
  write ( *, '(a)' ) ' '

  do y = 1990, 2000
    call y_to_s_english ( y, sy )
    write ( *, '(2x,a,2x,l1)' ) sy, year_is_leap_english ( y )
  end do

  return
end
subroutine test535 ( )

!*****************************************************************************80
!
!! TEST535 tests YEAR_IS_LEAP_GREEK.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 10 ) sy
  integer ( kind = 4 ) y
  logical year_is_leap_greek

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST535'
  write ( *, '(a)' ) '  For the Greek calendar:'
  write ( *, '(a)' ) '  YEAR_IS_LEAP_GREEK reports leap years.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Year  Leap?'
  write ( *, '(a)' ) ' '

  do y = 1, 10
    call y_to_s_greek ( y, sy )
    write ( *, '(2x,a,2x,l1)' ) sy, year_is_leap_greek ( y )
  end do

  return
end
subroutine test54 ( )

!*****************************************************************************80
!
!! TEST54 tests YEAR_IS_LEAP_GREGORIAN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 10 ) sy
  integer ( kind = 4 ) y
  logical year_is_leap_gregorian

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST54'
  write ( *, '(a)' ) '  For the Gregorian calendar:'
  write ( *, '(a)' ) '  YEAR_IS_LEAP_GREGORIAN reports leap years.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Year  Leap?'
  write ( *, '(a)' ) ' '

  do y = 1990, 2000
    call y_to_s_gregorian ( y, sy )
    write ( *, '(2x,a,2x,l1)' ) sy, year_is_leap_gregorian ( y )
  end do

  return
end
subroutine test555 ( )

!*****************************************************************************80
!
!! TEST555 tests YEAR_IS_LEAP_ISLAMIC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 10 ) sy
  integer ( kind = 4 ) y
  logical year_is_leap_islamic

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST555'
  write ( *, '(a)' ) '  For the Islamic calendar:'
  write ( *, '(a)' ) '  YEAR_IS_LEAP_ISLAMIC reports leap years.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Year  Leap?'
  write ( *, '(a)' ) ' '

  do y = 500, 510
    call y_to_s_islamic ( y, sy )
    write ( *, '(2x,a,2x,l1)' ) trim ( sy ), year_is_leap_islamic ( y )
  end do

  return
end
subroutine test56 ( )

!*****************************************************************************80
!
!! TEST56 tests YEAR_IS_LEAP_JULIAN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 10 ) sy
  integer ( kind = 4 ) y
  logical year_is_leap_julian

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST56'
  write ( *, '(a)' ) '  For the Julian calendar:'
  write ( *, '(a)' ) '  YEAR_IS_LEAP_JULIAN reports leap years.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Year  Leap?'
  write ( *, '(a)' ) ' '

  do y = 1990, 2000
    call y_to_s_julian ( y, sy )
    write ( *, '(2x,a,2x,l1)' ) sy, year_is_leap_julian ( y )
  end do

  return
end
subroutine test565 ( )

!*****************************************************************************80
!
!! TEST565 tests YEAR_IS_LEAP_REPUBLICAN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 7 ) sy
  integer ( kind = 4 ) y
  logical year_is_leap_republican

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST565'
  write ( *, '(a)' ) '  For the Republican calendar:'
  write ( *, '(a)' ) '  YEAR_IS_LEAP_REPUBLICAN reports leap years.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Year  Leap?'
  write ( *, '(a)' ) ' '

  do y = 1, 6
    call y_to_s_republican ( y, sy )
    write ( *, '(2x,a,2x,l1)' ) sy, year_is_leap_republican ( y )
  end do

  return
end
subroutine test566 ( )

!*****************************************************************************80
!
!! TEST566 tests YEAR_IS_LEAP_ROMAN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 20 ) sy
  integer ( kind = 4 ) y
  logical year_is_leap_roman

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST566'
  write ( *, '(a)' ) '  For the Roman calendar:'
  write ( *, '(a)' ) '  YEAR_IS_LEAP_ROMAN reports leap years.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Year  Leap?'
  write ( *, '(a)' ) ' '

  do y = 96, 100
    call y_to_s_roman ( y, sy )
    write ( *, '(2x,a,2x,l1)' ) sy, year_is_leap_roman ( y )
  end do

  return
end
subroutine test57 ( )

!*****************************************************************************80
!
!! TEST57 tests YEAR_LENGTH_COMMON.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 10 ) sy
  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length_common

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST57'
  write ( *, '(a)' ) '  For the Common calendar:'
  write ( *, '(a)' ) '  YEAR_LENGTH_COMMON determines the length of a year.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Year  Length'
  write ( *, '(a)' ) ' '

  do y = 1580, 1585
    call y_to_s_common ( y, sy )
    write ( *, '(2x,a,2x,i6)' ) sy, year_length_common ( y )
  end do

  do y = 1750, 1755
    call y_to_s_common ( y, sy )
    write ( *, '(2x,a,2x,i6)' ) sy, year_length_common ( y )
  end do

  do y = 1000, 2000, 100
    call y_to_s_common ( y, sy )
    write ( *, '(2x,a,2x,i6)' ) sy, year_length_common ( y )
  end do

  return
end
subroutine test58 ( )

!*****************************************************************************80
!
!! TEST58 tests YEAR_LENGTH_ENGLISH.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 10 ) sy
  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length_english

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST58'
  write ( *, '(a)' ) '  For the English calendar:'
  write ( *, '(a)' ) '  YEAR_LENGTH_ENGLISH determines the length of a year.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Year  Length'
  write ( *, '(a)' ) ' '

  do y = 1580, 1585
    call y_to_s_english ( y, sy )
    write ( *, '(2x,a,2x,i6)' ) sy, year_length_english ( y )
  end do

  do y = 1750, 1755
    call y_to_s_english ( y, sy )
    write ( *, '(2x,a,2x,i6)' ) sy, year_length_english ( y )
  end do

  do y = 1000, 2000, 100
    call y_to_s_english ( y, sy )
    write ( *, '(2x,a,2x,i6)' ) sy, year_length_english ( y )
  end do

  return
end
subroutine test585 ( )

!*****************************************************************************80
!
!! TEST585 tests YEAR_LENGTH_GREEK.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 10 ) sy
  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length_greek

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST585'
  write ( *, '(a)' ) '  For the Greek calendar:'
  write ( *, '(a)' ) '  YEAR_LENGTH_GREEK determines the length of a year.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Year  Length'
  write ( *, '(a)' ) ' '

  do y = 1, 10
    call y_to_s_greek ( y, sy )
    write ( *, '(2x,a,2x,i6)' ) sy, year_length_greek ( y )
  end do

  return
end
subroutine test59 ( )

!*****************************************************************************80
!
!! TEST59 tests YEAR_LENGTH_GREGORIAN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 10 ) sy
  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length_gregorian

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST59'
  write ( *, '(a)' ) '  For the Gregorian calendar:'
  write ( *, '(a)' ) '  YEAR_LENGTH_GREGORIAN determines the length of a year.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Year  Length'
  write ( *, '(a)' ) ' '

  do y = 1580, 1585
    call y_to_s_gregorian ( y, sy )
    write ( *, '(2x,a,2x,i8)' ) sy, year_length_gregorian ( y )
  end do

  do y = 1750, 1755
    call y_to_s_gregorian ( y, sy )
    write ( *, '(2x,a,2x,i8)' ) sy, year_length_gregorian ( y )
  end do

  do y = 1000, 2000, 100
    call y_to_s_gregorian ( y, sy )
    write ( *, '(2x,a,2x,i8)' ) sy, year_length_gregorian ( y )
  end do

  return
end
subroutine test60 ( )

!*****************************************************************************80
!
!! TEST60 tests YEAR_LENGTH_HEBREW.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 10 ) sy
  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length_hebrew

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST60'
  write ( *, '(a)' ) '  For the Hebrew calendar,'
  write ( *, '(a)' ) '  YEAR_LENGTH_HEBREW determines the length of a year.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Year  Length'
  write ( *, '(a)' ) ' '

  do y = 5760, 5780
    call y_to_s_hebrew ( y, sy )
    write ( *, '(2x,a,2x,i6)' ) trim ( sy ), year_length_hebrew ( y )
  end do

  return
end
subroutine test605 ( )

!*****************************************************************************80
!
!! TEST605 tests YEAR_LENGTH_ISLAMIC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 10 ) sy
  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length_islamic

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST605'
  write ( *, '(a)' ) '  For the Islamic calendar:'
  write ( *, '(a)' ) '  YEAR_LENGTH_ISLAMIC determines the length of a year.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Year  Length'
  write ( *, '(a)' ) ' '

  do y = 500, 505
    call y_to_s_islamic ( y, sy )
    write ( *, '(2x,a,2x,i6)' ) trim ( sy ), year_length_islamic ( y )
  end do

  return
end
subroutine test61 ( )

!*****************************************************************************80
!
!! TEST61 tests YEAR_LENGTH_JULIAN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 10 ) sy
  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length_julian

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST61'
  write ( *, '(a)' ) '  For the Julian calendar:'
  write ( *, '(a)' ) '  YEAR_LENGTH_JULIAN determines the length of a year.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Year  Length'
  write ( *, '(a)' ) ' '

  do y = 1580, 1585
    call y_to_s_julian ( y, sy )
    write ( *, '(2x,a,2x,i6)' ) sy, year_length_julian ( y )
  end do

  do y = 1750, 1755
    call y_to_s_julian ( y, sy )
    write ( *, '(2x,a,2x,i6)' ) sy, year_length_julian ( y )
  end do

  do y = 1000, 2000, 100
    call y_to_s_julian ( y, sy )
    write ( *, '(2x,a,2x,i6)' ) sy, year_length_julian ( y )
  end do

  return
end
subroutine test615 ( )

!*****************************************************************************80
!
!! TEST615 tests YEAR_LENGTH_REPUBLICAN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 10 ) sy
  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length_republican

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST615'
  write ( *, '(a)' ) '  For the Republican calendar:'
  write ( *, '(a)' ) '  YEAR_LENGTH_REPUBLICAN determines the length of a year.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Year  Length'
  write ( *, '(a)' ) ' '

  do y = 1, 6
    call y_to_s_republican ( y, sy )
    write ( *, '(2x,a,2x,i6)' ) sy, year_length_republican ( y )
  end do

  return
end
subroutine test616 ( )

!*****************************************************************************80
!
!! TEST616 tests YEAR_LENGTH_ROMAN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 20 ) sy
  integer ( kind = 4 ) y
  integer ( kind = 4 ) year_length_roman

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST616'
  write ( *, '(a)' ) '  For the Roman calendar:'
  write ( *, '(a)' ) '  YEAR_LENGTH_ROMAN determines the length of a year.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Year  Length'
  write ( *, '(a)' ) ' '

  do y = 96, 100
    call y_to_s_roman ( y, sy )
    write ( *, '(2x,a,2x,i6)' ) sy, year_length_roman ( y )
  end do

  return
end
subroutine test62 ( )

!*****************************************************************************80
!
!! TEST62 tests YEAR_TO_DOMINICAL_COMMON.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character d1
  character d2
  character i4_to_a
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  character ( len = 10 ) s
  integer ( kind = 4 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST62'
  write ( *, '(a)' ) '  For the Common calendar,'
  write ( *, '(a)' ) '  YEAR_TO_DOMINICAL_COMMON determines the dominical number.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Year  Dominical Number'
  write ( *, '(a)' ) ' '

  do y = 1577, 1587

    call y_to_s_common ( y, s )
    call year_to_dominical_common ( y, n1, n2 )
    d1 = i4_to_a ( n1 )

    if ( n1 == n2 ) then
      write ( *, '(2x,a,2x,i1,2x,a1)' ) s, n1, d1
    else
      d2 = i4_to_a ( n2 )
      write ( *, '(2x,a,2x,i1,2x,a1,2x,i1,2x,a1)' ) s, n1, d1, n2, d2
    end if

  end do

  return
end
subroutine test621 ( )

!*****************************************************************************80
!
!! TEST621 tests YEAR_TO_DOMINICAL_GREGORIAN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character d1
  character d2
  character i4_to_a
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  character ( len = 10 ) s
  integer ( kind = 4 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST621'
  write ( *, '(a)' ) '  For the Gregorian calendar,'
  write ( *, '(a)' ) '  YEAR_TO_DOMINICAL_GREGORIAN determines the dominical number.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Year  Dominical Number'
  write ( *, '(a)' ) ' '

  do y = 1577, 1587

    call y_to_s_gregorian ( y, s )

    call year_to_dominical_gregorian ( y, n1, n2 )
    d1 = i4_to_a ( n1 )

    if ( n1 == n2 ) then
      write ( *, '(2x,a,2x,i1,2x,a1)' ) s, n1, d1
    else
      d2 = i4_to_a ( n2 )
      write ( *, '(2x,a,2x,i1,2x,a1,2x,i1,2x,a1)' ) s, n1, d1, n2, d2
    end if

  end do

  return
end
subroutine test622 ( )

!*****************************************************************************80
!
!! TEST622 tests YEAR_TO_DOMINICAL_JULIAN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character d1
  character d2
  character i4_to_a
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  character ( len = 10 ) s
  integer ( kind = 4 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST622'
  write ( *, '(a)' ) '  For the Julian calendar,'
  write ( *, '(a)' ) '  YEAR_TO_DOMINICAL_JULIAN determines the dominical number.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Year  Dominical Number'
  write ( *, '(a)' ) ' '

  do y = 1577, 1587

    call y_to_s ( y, s )
    call year_to_dominical_julian ( y, n1, n2 )
    d1 = i4_to_a ( n1 )

    if ( n1 == n2 ) then
      write ( *, '(2x,a,2x,i1,2x,a1)' ) s, n1, d1
    else
      d2 = i4_to_a ( n2 )
      write ( *, '(2x,a,2x,i1,2x,a1,2x,i1,2x,a1)' ) s, n1, d1, n2, d2
    end if

  end do

  return
end
subroutine test623 ( )

!*****************************************************************************80
!
!! TEST623 tests YEAR_TO_EPACT_GREGORIAN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) e
  character ( len = 10 ) s
  integer ( kind = 4 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST623'
  write ( *, '(a)' ) '  For the Gregorian calendar,'
  write ( *, '(a)' ) '  YEAR_TO_EPACT_GREGORIAN determines the epact of a year.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Year  Epact'
  write ( *, '(a)' ) ' '

  do y = -2, 20
    if ( y /= 0 ) then
      call y_to_s_gregorian ( y, s )
      call year_to_epact_gregorian ( y, e )
      write ( *, '(2x,a,2x,i6)' ) s, e
    end if
  end do

  return
end
subroutine test624 ( )

!*****************************************************************************80
!
!! TEST624 tests YEAR_TO_EPACT_JULIAN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) e
  character ( len = 10 ) s
  integer ( kind = 4 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST624'
  write ( *, '(a)' ) '  For the Julian calendar,'
  write ( *, '(a)' ) '  YEAR_TO_EPACT_JULIAN determines the epact of a year.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Year  Epact'
  write ( *, '(a)' ) ' '

  do y = -2, 20
    if ( y /= 0 ) then
      call y_to_s_julian ( y, s )
      call year_to_epact_julian ( y, e )
      write ( *, '(2x,a,2x,i6)' ) s, e
    end if
  end do

  return
end
subroutine test63 ( )

!*****************************************************************************80
!
!! TEST63 tests YEAR_TO_GOLDEN_NUMBER.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) g
  character ( len = 10 ) s
  integer ( kind = 4 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST63'
  write ( *, '(a)' ) '  YEAR_TO_GOLDEN_NUMBER determines the golden'
  write ( *, '(a)' ) '    number of a year.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Year  Golden Number'
  write ( *, '(a)' ) ' '

  do y = -2, 20
    if ( y /= 0 ) then
      call y_to_s_common ( y, s )
      call year_to_golden_number ( y, g )
      write ( *, '(2x,a,2x,i6)' ) s, g
    end if
  end do

  return
end
subroutine test635 ( )

!*****************************************************************************80
!
!! TEST635 tests YEAR_TO_INDICTION_COMMON.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  character ( len = 10 ) sy
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST635'
  write ( *, '(a)' ) '  For a Common year,'
  write ( *, '(a)' ) '  YEAR_TO_INDICTION_COMMON determines the indiction number.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Year  Indiction Number'
  write ( *, '(a)' ) ' '

  do y = -1, 13

    call y_astronomical_to_common ( y, y2 )

    call y_to_s_common ( y2, sy )
    call year_to_indiction_common ( y2, i )
    write ( *, '(4x,a,2x,i2)' ) sy, i

  end do

  return
end
subroutine test636 ( )

!*****************************************************************************80
!
!! TEST636 tests YEAR_TO_SCALIGER_COMMON.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) c1
  integer ( kind = 4 ) c2
  integer ( kind = 4 ) c3
  integer ( kind = 4 ) r1
  integer ( kind = 4 ) r2
  integer ( kind = 4 ) r3
  character ( len = 10 ) sy
  integer ( kind = 4 ) y
  integer ( kind = 4 ) y2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST636'
  write ( *, '(a)' ) '  For a Common year,'
  write ( *, '(a)' ) '  YEAR_TO_SCALIGER_COMMON determines the Scaliger indices.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Year  Julian / Metonic / Indiction'
  write ( *, '(a)' ) ' '

  do y = -4713, -4675

    call y_astronomical_to_common ( y, y2 )

    call y_to_s_common ( y2, sy )
    call year_to_scaliger_common ( y2, c1, c2, c3, r1, r2, r3 )
    write ( *, '(4x,a,2x,2i5,2x,2i5,2x,2i5)' ) sy, c1, r1, c2, r2, c3, r3

  end do

  return
end
subroutine test64 ( )

!*****************************************************************************80
!
!! TEST64 tests YEAR_TO_TYPE_HEBREW.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 10 ) s
  integer ( kind = 4 ) type
  integer ( kind = 4 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST64'
  write ( *, '(a)' ) '  For the Hebrew calendar,'
  write ( *, '(a)' ) '  YEAR_TO_TYPE_HEBREW determines the type of a year.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Year  TYPE'
  write ( *, '(a)' ) ' '

  do y = 5760, 5780
    call y_to_s_hebrew ( y, s )
    call year_to_type_hebrew ( y, type )
    write ( *, '(2x,a,2x,i6)' ) trim ( s ), type
  end do

  return
end
subroutine test65 ( )

!*****************************************************************************80
!
!! TEST65 tests YJF_DIF_COMMON.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) days
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  real ( kind = 8 ), parameter :: fhi = 0.0D+00
  real ( kind = 8 ), parameter :: flo = 0.0D+00
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ), parameter :: jhi = 1
  integer ( kind = 4 ), parameter :: jlo = 1
  character ( len = 20 ) s1
  character ( len = 20 ) s2
  integer ( kind = 4 ) :: seed = 123456789
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2
  integer ( kind = 4 ), parameter :: yhi = 1970
  integer ( kind = 4 ), parameter :: ylo = 1960

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST65'
  write ( *, '(a)' ) '  For Common calendar YJF dates,'
  write ( *, '(a)' ) '  YJF_DIF_COMMON computes the day difference.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  YJF1           YJF2      (YJF2 - YJF1)'
  write ( *, '(a)' ) ' '

  do i = 1, 10
 
    call yjf_uniform_common ( ylo, jlo, flo, yhi, jhi, fhi, seed, y1, j1, f1 )
    call yjf_to_s_common ( y1, j1, f1, s1 )

    call yjf_uniform_common ( ylo, jlo, flo, yhi, jhi, fhi, seed, y2, j2, f2 )
    call yjf_to_s_common ( y2, j2, f2, s2 )

    call yjf_dif_common ( y1, j1, f1, y2, j2, f2, days, ierror )
 
    write ( *, '(2x,a,5x,a,5x,f11.2)' ) s1, s2, days

  end do
 
  return
end
subroutine test66 ( )

!*****************************************************************************80
!
!! TEST66 tests YJF_TO_WEEKDAY_COMMON.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d1
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j2
  real ( kind = 8 ) jed1
  real ( kind = 8 ) jed2
  integer ( kind = 4 ) m1
  character ( len = 20 ) s1
  character ( len = 13 ) s2
  character ( len = 11 ) s3
  integer ( kind = 4 ) w3
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST66'
  write ( *, '(a)' ) '  For the Common calendar,'
  write ( *, '(a)' ) '  YJF_TO_WEEKDAY_COMMON reports day of week for a YJF date.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      JED          YMDF           YJF      W  Name'
  write ( *, '(a)' ) ' '

  i = 0

  do

    i = i + 1
    call jed_test ( i, jed1 )

    if ( jed1 < 0 ) then
      exit
    end if

    call jed_to_next_noon ( jed1, jed2 )

    call jed_to_ymdf_common ( jed2, y1, m1, d1, f1 )

    call ymdf_to_s_common ( y1, m1, d1, f1, s1 )
 
    call ymdf_to_yjf_common ( y1, m1, d1, f1, y2, j2, f2 ) 

    call yjf_to_s_common ( y2, j2, f2, s2 )

    call yjf_to_weekday_common ( y2, j2, f2, w3 )
 
    call weekday_to_name_common ( w3, s3 )

    write ( *, '(f11.2,2x,a,2x,a,2x,i1,2x,a)' ) jed2, s1, s2, w3, s3

  end do

  return
end
subroutine test67 ( )

!*****************************************************************************80
!
!! TEST67 tests YJF_TO_YMDF_COMMON and YMDF_TO_YJF_COMMON.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d1
  integer ( kind = 4 ) d3
  integer ( kind = 4 ), parameter :: dlo = 1
  integer ( kind = 4 ), parameter :: dhi = 1
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  real ( kind = 8 ) f3
  real ( kind = 8 ), parameter :: flo = 0.0
  real ( kind = 8 ), parameter :: fhi = 0.0
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m3
  integer ( kind = 4 ), parameter :: mlo = 1
  integer ( kind = 4 ), parameter :: mhi = 1
  integer ( kind = 4 ) :: seed = 123456789
  character ( len = 20 ) s1
  character ( len = 20 ) s2
  character ( len = 20 ) s3
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2
  integer ( kind = 4 ) y3
  integer ( kind = 4 ), parameter :: ylo = 1960
  integer ( kind = 4 ), parameter :: yhi = 1970

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST67'
  write ( *, '(a)' ) '  For the Common calendar,'
  write ( *, '(a)' ) '  YJF_TO_YMDF_COMMON: YJF => YMDF.'
  write ( *, '(a)' ) '  YMDF_TO_YJF_COMMON: YMDF => YJF.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  YMDF(in)         YJF        YMDF(out)'
  write ( *, '(a)' ) ' '

  do i = 1, 10

    call ymdf_uniform_common ( ylo, mlo, dlo, flo, yhi, mhi, dhi, fhi, &
      seed, y1, m1, d1, f1 )

    call ymdf_to_s_common ( y1, m1, d1, f1, s1 )

    call ymdf_to_yjf_common ( y1, m1, d1, f1, y2, j2, f2 )

    call yjf_to_s_common ( y2, j2, f2, s2 )

    call yjf_to_ymdf_common ( y2, j2, f2, y3, m3, d3, f3 )

    call ymdf_to_s_common ( y3, m3, d3, f3, s3 )

    write ( *, '(2x,3a)' ) s1, s2, s3

  end do

  return
end
subroutine test675 ( )

!*****************************************************************************80
!
!! TEST675 tests YJF_TO_YMDF_ENGLISH and YMDF_TO_YJF_ENGLISH.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d1
  integer ( kind = 4 ) d3
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  real ( kind = 8 ) f3
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j2
  real ( kind = 8 ) jed
  real ( kind = 8 ) jed_epoch
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m3
  character ( len = 20 ) s1
  character ( len = 20 ) s2
  character ( len = 20 ) s3
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2
  integer ( kind = 4 ) y3

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST675'
  write ( *, '(a)' ) '  For the English calendar,'
  write ( *, '(a)' ) '  YJF_TO_YMDF_ENGLISH: YJF => YMDF.'
  write ( *, '(a)' ) '  YMDF_TO_YJF_ENGLISH: YMDF => YJF.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      JED    YMDF(in)         YJF        YMDF(out)'
  write ( *, '(a)' ) ' '

  call epoch_to_jed_english ( jed_epoch )

  i = 0

  do

    i = i + 1
    call jed_test ( i, jed )

    if ( jed < 0 ) then
      exit
    end if

    if ( jed_epoch <= jed ) then

      call jed_to_ymdf_english ( jed, y1, m1, d1, f1 )
      
      call ymdf_to_s_english ( y1, m1, d1, f1, s1 )

      call ymdf_to_yjf_english ( y1, m1, d1, f1, y2, j2, f2 )

      call yjf_to_s_english ( y2, j2, f2, s2 )

      call yjf_to_ymdf_english ( y2, j2, f2, y3, m3, d3, f3 )

      call ymdf_to_s_english ( y3, m3, d3, f3, s3 )

      write ( *, '(f11.2,2x,3a)' ) jed, s1, s2, s3

    end if

  end do

  return
end
subroutine test68 ( )

!*****************************************************************************80
!
!! TEST68 tests YJF_TO_YMDF_HEBREW and YMDF_TO_YJF_HEBREW.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d1
  integer ( kind = 4 ) d3
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  real ( kind = 8 ) f3
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j2
  real ( kind = 8 ) jed_epoch
  real ( kind = 8 ) jed1
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m3
  character ( len = 20 ) s1
  character ( len = 15 ) s2
  character ( len = 20 ) s3
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2
  integer ( kind = 4 ) y3

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST68'
  write ( *, '(a)' ) '  For the Hebrew calendar,'
  write ( *, '(a)' ) '  YJF_TO_YMDF_HEBREW: YJF => YMDF'
  write ( *, '(a)' ) '  YMDF_TO_YJF_HEBREW: YMDF => YJF'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  YMDF(in)         YJF        YMDF(out)'
  write ( *, '(a)' ) ' '

  call epoch_to_jed_hebrew ( jed_epoch )

  i = 0

  do
 
    i = i + 1
    call jed_test ( i, jed1 )

    if ( jed1 < 0.0D+00 ) then
      exit
    end if

    if ( jed_epoch <= jed1 ) then

      call jed_to_ymdf_hebrew ( jed1, y1, m1, d1, f1 )

      call ymdf_to_s_hebrew ( y1, m1, d1, f1, s1 )

      call ymdf_to_yjf_hebrew ( y1, m1, d1, f1, y2, j2, f2 )

      call yjf_to_s_hebrew ( y2, j2, f2, s2 )

      call yjf_to_ymdf_hebrew ( y2, j2, f2, y3, m3, d3, f3 )

      call ymdf_to_s_hebrew ( y3, m3, d3, f3, s3 )

      write ( *, '(2x,3a)' ) s1, s2, s3

    end if

  end do

  return
end
subroutine test685 ( )

!*****************************************************************************80
!
!! TEST685 tests YJF_TO_YMDF_ISLAMIC and YMDF_TO_YJF_ISLAMIC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d1
  integer ( kind = 4 ) d3
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  real ( kind = 8 ) f3
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j2
  real ( kind = 8 ) jed_epoch
  real ( kind = 8 ) jed1
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m3
  character ( len = 20 ) s1
  character ( len = 15 ) s2
  character ( len = 20 ) s3
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2
  integer ( kind = 4 ) y3

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST685'
  write ( *, '(a)' ) '  For the Islamic calendar,'
  write ( *, '(a)' ) '  YJF_TO_YMDF_ISLAMIC: YJF => YMDF'
  write ( *, '(a)' ) '  YMDF_TO_YJF_ISLAMIC: YMDF => YJF'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  YMDF(in)         YJF        YMDF(out)'
  write ( *, '(a)' ) ' '

  call epoch_to_jed_islamic_a ( jed_epoch )

  i = 0

  do
 
    i = i + 1
    call jed_test ( i, jed1 )

    if ( jed1 < 0.0D+00 ) then
      exit
    end if

    if ( jed_epoch <= jed1 ) then

      call jed_to_ymdf_islamic_a ( jed1, y1, m1, d1, f1 )

      call ymdf_to_s_islamic ( y1, m1, d1, f1, s1 )

      call ymdf_to_yjf_islamic ( y1, m1, d1, f1, y2, j2, f2 )

      call yjf_to_s_islamic ( y2, j2, f2, s2 )

      call yjf_to_ymdf_islamic ( y2, j2, f2, y3, m3, d3, f3 )

      call ymdf_to_s_islamic ( y3, m3, d3, f3, s3 )

      write ( *, '(2x,3a)' ) s1, s2, s3

    end if

  end do

  return
end
subroutine test686 ( )

!*****************************************************************************80
!
!! TEST686 tests YJF_TO_YMDF_JULIAN and YMDF_TO_YJF_JULIAN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d1
  integer ( kind = 4 ) d3
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  real ( kind = 8 ) f3
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j2
  real ( kind = 8 ) jed_epoch
  real ( kind = 8 ) jed1
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m3
  character ( len = 20 ) s1
  character ( len = 15 ) s2
  character ( len = 20 ) s3
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2
  integer ( kind = 4 ) y3

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST686'
  write ( *, '(a)' ) '  For the Julian calendar,'
  write ( *, '(a)' ) '  YJF_TO_YMDF_JULIAN: YJF => YMDF'
  write ( *, '(a)' ) '  YMDF_TO_YJF_JULIAN: YMDF => YJF'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  YMDF(in)         YJF        YMDF(out)'
  write ( *, '(a)' ) ' '

  call epoch_to_jed_julian ( jed_epoch )

  i = 0

  do
 
    i = i + 1
    call jed_test ( i, jed1 )

    if ( jed1 < 0.0D+00 ) then
      exit
    end if

    if ( jed_epoch <= jed1 ) then

      call jed_to_ymdf_julian ( jed1, y1, m1, d1, f1 )

      call ymdf_to_s_julian ( y1, m1, d1, f1, s1 )

      call ymdf_to_yjf_julian ( y1, m1, d1, f1, y2, j2, f2 )

      call yjf_to_s_julian ( y2, j2, f2, s2 )

      call yjf_to_ymdf_julian ( y2, j2, f2, y3, m3, d3, f3 )

      call ymdf_to_s_julian ( y3, m3, d3, f3, s3 )

      write ( *, '(2x,3a)' ) s1, s2, s3

    end if

  end do

  return
end
subroutine test687 ( )

!*****************************************************************************80
!
!! TEST687 tests YJF_TO_YMDF_ROMAN and YMDF_TO_YJF_ROMAN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d1
  integer ( kind = 4 ) d3
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  real ( kind = 8 ) f3
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j2
  real ( kind = 8 ) jed_epoch
  real ( kind = 8 ) jed1
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m3
  character ( len = 50 ) s1
  character ( len = 15 ) s2
  character ( len = 50 ) s3
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2
  integer ( kind = 4 ) y3

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST687'
  write ( *, '(a)' ) '  For the Roman calendar,'
  write ( *, '(a)' ) '  YJF_TO_YMDF_ROMAN: YJF => YMDF'
  write ( *, '(a)' ) '  YMDF_TO_YJF_ROMAN: YMDF => YJF'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  YMDF(in)                                YJF'
  write ( *, '(a)' ) '  YMDF(out)'
  write ( *, '(a)' ) ' '

  call epoch_to_jed_roman ( jed_epoch )

  i = 0

  do
 
    i = i + 1
    call jed_test ( i, jed1 )

    if ( jed1 < 0.0D+00 ) then
      exit
    end if

    if ( jed_epoch <= jed1 ) then

      call jed_to_ymdf_roman ( jed1, y1, m1, d1, f1 )

      call ymdf_to_s_roman ( y1, m1, d1, f1, s1 )

      call ymdf_to_yjf_roman ( y1, m1, d1, f1, y2, j2, f2 )

      call yjf_to_s_roman ( y2, j2, f2, s2 )

      call yjf_to_ymdf_roman ( y2, j2, f2, y3, m3, d3, f3 )

      call ymdf_to_s_roman ( y3, m3, d3, f3, s3 )

      write ( *, '(a)' ) ' '
      write ( *, '(2x,3a)' ) s1, s2
      write ( *, '(2x,3a)' ) s3

    end if

  end do

  return
end
subroutine test688 ( )

!*****************************************************************************80
!
!! TEST688 tests YJF_TO_YMDHMS_COMMON and YMDHMS_TO_YJF_COMMON.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d2
  real ( kind = 8 ) f1
  real ( kind = 8 ) f3
  real ( kind = 8 ), parameter :: flo = 0.0D+00
  real ( kind = 8 ), parameter :: fhi = 0.0D+00
  integer ( kind = 4 ) h2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j3
  integer ( kind = 4 ), parameter :: jlo = 1
  integer ( kind = 4 ), parameter :: jhi = 1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) s2
  integer ( kind = 4 ) :: seed = 123456789
  character ( len = 20 ) ss1
  character ( len = 22 ) ss2
  character ( len = 20 ) ss3
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2
  integer ( kind = 4 ) y3
  integer ( kind = 4 ), parameter :: ylo = 1960
  integer ( kind = 4 ), parameter :: yhi = 1970

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST688'
  write ( *, '(a)' ) '  YJF_TO_YMDHMS_COMMON: YJF => YMDHMS'
  write ( *, '(a)' ) '  YMDHMS_TO_YJF_COMMON: YMDHMS => YJF'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  YJF (in)         YMDHMS         YJF(out)'
  write ( *, '(a)' ) ' '

  do i = 1, 5

    call yjf_uniform_common ( ylo, jlo, flo, yhi, jhi, fhi, seed, &
      y1, j1, f1 )

    call yjf_to_s_common ( y1, j1, f1, ss1 )

    call yjf_to_ymdhms_common ( y1, j1, f1, y2, m2, d2, h2, n2, s2 )
    call ymdhms_to_s_common ( y2, m2, d2, h2, n2, s2, ss2 )

    call ymdhms_to_yjf_common ( y2, m2, d2, h2, n2, s2, y3, j3, f3 )
    call yjf_to_s_common ( y3, j3, f3, ss3 )

    write ( *, '(2x,a,2x,a,2x,a)' ) ss1, ss2, ss3

  end do

  return
end
subroutine test689 ( )

!*****************************************************************************80
!
!! TEST689 tests YMD_TO_DECIMAL
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 October 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ), parameter :: dhi = 1
  integer ( kind = 4 ), parameter :: dlo = 1
  real ( kind = 8 ) f
  real ( kind = 8 ), parameter :: fhi = 0.0D+00
  real ( kind = 8 ), parameter :: flo = 0.0D+00
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) m
  integer ( kind = 4 ), parameter :: mhi = 1
  integer ( kind = 4 ), parameter :: mlo = 1
  character ( len = 20 ) s
  integer ( kind = 4 ) :: seed = 123456789
  integer ( kind = 4 ) y
  real ( kind = 8 ) yf
  integer ( kind = 4 ), parameter :: yhi = 1970
  integer ( kind = 4 ), parameter :: ylo = 1960

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST689'
  write ( *, '(a)' ) '  YMD_TO_DECIMAL converts a date to a year and decimal.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  YMD                            Y.F'
  write ( *, '(a)' ) ' '
 
  do i = 1, 10
 
    call ymdf_uniform_common ( ylo, mlo, dlo, flo, yhi, mhi, dhi, fhi, &
      seed, y, m, d, f )

    call ymd_to_s_common ( y, m, d, s )

    call ymd_to_decimal ( y, m, d, yf )

    write ( *, '(2x,a,5x,f14.4)' ) s, yf

  end do

  return
end
subroutine test69 ( )

!*****************************************************************************80
!
!! TEST69 tests YMDF_DIF_COMMON.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d1
  integer ( kind = 4 ) d2
  real ( kind = 8 ) days
  integer ( kind = 4 ), parameter :: dhi = 1
  integer ( kind = 4 ), parameter :: dlo = 1
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  real ( kind = 8 ), parameter :: fhi = 0.0D+00
  real ( kind = 8 ), parameter :: flo = 0.0D+00
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ), parameter :: mhi = 1
  integer ( kind = 4 ), parameter :: mlo = 1
  character ( len = 20 ) s1
  character ( len = 20 ) s2
  integer ( kind = 4 ) :: seed = 123456789
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2
  integer ( kind = 4 ), parameter :: yhi = 1960
  integer ( kind = 4 ), parameter :: ylo = 1970

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST69'
  write ( *, '(a)' ) '  YMDF_DIF_COMMON gets the day difference '
  write ( *, '(a)' ) '  between YMDF dates.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  YMDF1        YMDF2        (YMDF2 - YMDF1)'
  write ( *, '(a)' ) ' '
 
  do i = 1, 10
 
    call ymdf_uniform_common ( ylo, mlo, dlo, flo, yhi, mhi, dhi, fhi, &
      seed, y1, m1, d1, f1 )

    call ymdf_to_s_common ( y1, m1, d1, f1, s1 )

    call ymdf_uniform_common ( ylo, mlo, dlo, flo, yhi, mhi, dhi, fhi, &
      seed, y2, m2, d2, f2 )

    call ymdf_to_s_common ( y2, m2, d2, f2, s2 )

    call ymdf_dif_common ( y1, m1, d1, f1, y2, m2, d2, f2, days, ierror )
 
    write ( *, '(2x,a,5x,a,5x,f11.2)' ) s1, s2, days

  end do

  return
end
subroutine test695 ( )

!*****************************************************************************80
!
!! TEST695 tests YMDF_DIF_ENGLISH.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d1
  integer ( kind = 4 ) d2
  real ( kind = 8 ) days
  integer ( kind = 4 ), parameter :: dhi = 1
  integer ( kind = 4 ), parameter :: dlo = 1
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  real ( kind = 8 ), parameter :: fhi = 0.0D+00
  real ( kind = 8 ), parameter :: flo = 0.0D+00
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ), parameter :: mhi = 1
  integer ( kind = 4 ), parameter :: mlo = 1
  character ( len = 20 ) s1
  character ( len = 20 ) s2
  integer ( kind = 4 ) :: seed = 123456789
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2
  integer ( kind = 4 ), parameter :: yhi = 1960
  integer ( kind = 4 ), parameter :: ylo = 1970

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST695'
  write ( *, '(a)' ) '  YMDF_DIF_ENGLISH gets the day difference'
  write ( *, '(a)' ) '  between YMDF dates.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  YMDF1        YMDF2        (YMDF2 - YMDF1)'
  write ( *, '(a)' ) ' '
 
  do i = 1, 10
 
    call ymdf_uniform_english ( ylo, mlo, dlo, flo, yhi, mhi, dhi, fhi, &
      seed, y1, m1, d1, f1 )

    call ymdf_to_s_english ( y1, m1, d1, f1, s1 )

    call ymdf_uniform_english ( ylo, mlo, dlo, flo, yhi, mhi, dhi, fhi, &
      seed, y2, m2, d2, f2 )

    call ymdf_to_s_english ( y2, m2, d2, f2, s2 )

    call ymdf_dif_english ( y1, m1, d1, f1, y2, m2, d2, f2, days, ierror )
 
    write ( *, '(2x,a,5x,a,5x,f11.2)' ) s1, s2, days

  end do

  return
end
subroutine test70 ( )

!*****************************************************************************80
!
!! TEST70 tests YMDF_INC_COMMON, YMDF_NEXT_COMMON, YMDF_PREV_COMMON.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d1
  integer ( kind = 4 ) d2
  integer ( kind = 4 ) d3
  integer ( kind = 4 ) d4
  real ( kind = 8 ), parameter :: days = 10.25D+00
  integer ( kind = 4 ), parameter :: dhi = 1
  integer ( kind = 4 ), parameter :: dlo = 1
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  real ( kind = 8 ) f3
  real ( kind = 8 ) f4
  real ( kind = 8 ), parameter :: fhi = 0.0D+00
  real ( kind = 8 ), parameter :: flo = 0.0D+00
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) m3
  integer ( kind = 4 ) m4
  integer ( kind = 4 ), parameter :: mhi = 1
  integer ( kind = 4 ), parameter :: mlo = 1
  character ( len = 20 ) s1
  character ( len = 20 ) s2
  character ( len = 20 ) s3
  character ( len = 20 ) s4
  integer ( kind = 4 ) :: seed = 123456789
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2
  integer ( kind = 4 ) y3
  integer ( kind = 4 ) y4
  integer ( kind = 4 ), parameter :: yhi = 1960
  integer ( kind = 4 ), parameter :: ylo = 1970

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST70'
  write ( *, '(a)' ) '  For the Common calendar:'
  write ( *, '(a)' ) '  YMDF_INC_COMMON increments a date by days;'
  write ( *, '(a)' ) '  YMDF_NEXT_COMMON computes the next day,'
  write ( *, '(a)' ) '  YMDF_PREV_COMMON computes the previous day.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  YMDF      Tomorrow       Yesterday      +10.25 days'
  write ( *, '(a)' ) ' '

  do i = 1, 10
 
    call ymdf_uniform_common ( ylo, mlo, dlo, flo, yhi, mhi, dhi, fhi, &
      seed, y1, m1, d1, f1 )

    call ymdf_to_s_common ( y1, m1, d1, f1, s1 )

    call ymdf_next_common ( y1, m1, d1, f1, y2, m2, d2, f2 )
    call ymdf_to_s_common ( y2, m2, d2, f2, s2 )

    call ymdf_prev_common ( y1, m1, d1, f1, y3, m3, d3, f3 )
    call ymdf_to_s_common ( y3, m3, d3, f3, s3 )
 
    call ymdf_inc_common ( y1, m1, d1, f1, days, y4, m4, d4, f4 )
    call ymdf_to_s_common ( y4, m4, d4, f4, s4 )

    write ( *, '(2x,4a)' ) s1, s2, s3, s4

  end do

  return
end
subroutine test71 ( )

!*****************************************************************************80
!
!! TEST71 tests YMDF_INC_ENGLISH, YMDF_NEXT_ENGLISH, YMDF_PREV_ENGLISH.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d1
  integer ( kind = 4 ) d2
  integer ( kind = 4 ) d3
  integer ( kind = 4 ) d4
  real ( kind = 8 ), parameter :: days = 10.25D+00
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  real ( kind = 8 ) f3
  real ( kind = 8 ) f4
  integer ( kind = 4 ) i
  real ( kind = 8 ) jed
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) m3
  integer ( kind = 4 ) m4
  character ( len = 20 ) s1
  character ( len = 20 ) s2
  character ( len = 20 ) s3
  character ( len = 20 ) s4
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2
  integer ( kind = 4 ) y3
  integer ( kind = 4 ) y4

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST71'
  write ( *, '(a)' ) '  For the English calendar:'
  write ( *, '(a)' ) '  YMDF_INC_ENGLISH increments a date by days;'
  write ( *, '(a)' ) '  YMDF_NEXT_ENGLISH computes the next day,'
  write ( *, '(a)' ) '  YMDF_PREV_ENGLISH computes the previous day.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  YMDF	  Tomorrow	 Yesterday	+10.25 days'
  write ( *, '(a)' ) ' '

  i = 0

  do
 
    i = i + 1
    call jed_test ( i, jed )

    if ( jed < 0 ) then
      exit
    end if

    call jed_to_ymdf_english ( jed, y1, m1, d1, f1 )

    call ymdf_to_s_english ( y1, m1, d1, f1, s1 )

    call ymdf_next_english ( y1, m1, d1, f1, y2, m2, d2, f2 )
    call ymdf_to_s_english ( y2, m2, d2, f2, s2 )

    call ymdf_prev_english ( y1, m1, d1, f1, y3, m3, d3, f3 )
    call ymdf_to_s_english ( y3, m3, d3, f3, s3 )
 
    call ymdf_inc_english ( y1, m1, d1, f1, days, y4, m4, d4, f4 )
    call ymdf_to_s_english ( y4, m4, d4, f4, s4 )

    write ( *, '(2x,4a)' ) s1, s2, s3, s4

  end do

  return
end
subroutine test72 ( )

!*****************************************************************************80
!
!! TEST72 tests YMDF_INC_GREGORIAN, YMDF_NEXT_GREGORIAN, YMDF_PREV_GREGORIAN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d1
  integer ( kind = 4 ) d2
  integer ( kind = 4 ) d3
  integer ( kind = 4 ) d4
  real ( kind = 8 ), parameter :: days = 10.25D+00
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  real ( kind = 8 ) f3
  real ( kind = 8 ) f4
  integer ( kind = 4 ) i
  real ( kind = 8 ) jed
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) m3
  integer ( kind = 4 ) m4
  character ( len = 20 ) s1
  character ( len = 20 ) s2
  character ( len = 20 ) s3
  character ( len = 20 ) s4
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2
  integer ( kind = 4 ) y3
  integer ( kind = 4 ) y4

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST72'
  write ( *, '(a)' ) '  For the Gregorian calendar:'
  write ( *, '(a)' ) '  YMDF_INC_GREGORIAN increments a date by days;'
  write ( *, '(a)' ) '  YMDF_NEXT_GREGORIAN computes the next day,'
  write ( *, '(a)' ) '  YMDF_PREV_GREGORIAN computes the previous day.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  YMDF      Tomorrow       Yesterday      +10.25 days'
  write ( *, '(a)' ) ' '

  i = 0

  do

    i = i + 1
    call jed_test ( i, jed )

    if ( jed < 0 ) then
      exit
    end if

    call jed_to_ymdf_gregorian ( jed, y1, m1, d1, f1 )

    call ymdf_to_s_gregorian ( y1, m1, d1, f1, s1 )

    call ymdf_next_gregorian ( y1, m1, d1, f1, y2, m2, d2, f2 )
    call ymdf_to_s_gregorian ( y2, m2, d2, f2, s2 )

    call ymdf_prev_gregorian ( y1, m1, d1, f1, y3, m3, d3, f3 )
    call ymdf_to_s_gregorian ( y3, m3, d3, f3, s3 )
 
    call ymdf_inc_gregorian ( y1, m1, d1, f1, days, y4, m4, d4, f4 )
    call ymdf_to_s_gregorian ( y4, m4, d4, f4, s4 )

    write ( *, '(2x,4a)' ) s1, s2, s3, s4

  end do

  return
end
subroutine test73 ( )

!*****************************************************************************80
!
!! TEST73 tests YMDF_INC_HEBREW, YMDF_NEXT_HEBREW, YMDF_PREV_HEBREW.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d1
  integer ( kind = 4 ) d2
  integer ( kind = 4 ) d3
  integer ( kind = 4 ) d4
  real ( kind = 8 ), parameter :: days = 10.25D+00
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  real ( kind = 8 ) f3
  real ( kind = 8 ) f4
  integer ( kind = 4 ) i
  real ( kind = 8 ) jed
  real ( kind = 8 ) jed_epoch
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) m3
  integer ( kind = 4 ) m4
  character ( len = 20 ) s1
  character ( len = 20 ) s2
  character ( len = 20 ) s3
  character ( len = 20 ) s4
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2
  integer ( kind = 4 ) y3
  integer ( kind = 4 ) y4

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST73'
  write ( *, '(a)' ) '  For the Hebrew calendar:'
  write ( *, '(a)' ) '  YMDF_INC_HEBREW increments a date by days;'
  write ( *, '(a)' ) '  YMDF_NEXT_HEBREW computes the next day,'
  write ( *, '(a)' ) '  YMDF_PREV_HEBREW computes the previous day.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  YMDF        Tomorrow          Yesterday         +10.25 days'
  write ( *, '(a)' ) ' '

  call epoch_to_jed_hebrew ( jed_epoch )

  i = 0

  do

    i = i + 1
    call jed_test ( i, jed )

    if ( jed < 0 ) then
      exit
    end if

    if ( jed_epoch + 1 <= jed ) then

      call jed_to_ymdf_hebrew ( jed, y1, m1, d1, f1 )

      call ymdf_to_s_hebrew ( y1, m1, d1, f1, s1 )

      call ymdf_next_hebrew ( y1, m1, d1, f1, y2, m2, d2, f2 )
      call ymdf_to_s_hebrew ( y2, m2, d2, f2, s2 )

      call ymdf_prev_hebrew ( y1, m1, d1, f1, y3, m3, d3, f3 )
      call ymdf_to_s_hebrew ( y3, m3, d3, f3, s3 )
 
      call ymdf_inc_hebrew ( y1, m1, d1, f1, days, y4, m4, d4, f4 )
      call ymdf_to_s_hebrew ( y4, m4, d4, f4, s4 )

      write ( *, '(2x,4a)' ) s1, s2, s3, trim ( s4 )

    end if

  end do

  return
end
subroutine test74 ( )

!*****************************************************************************80
!
!! TEST74 tests YMDF_INC_JULIAN, YMDF_NEXT_JULIAN, YMDF_PREV_JULIAN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d1
  integer ( kind = 4 ) d2
  integer ( kind = 4 ) d3
  integer ( kind = 4 ) d4
  real ( kind = 8 ), parameter :: days = 10.25D+00
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  real ( kind = 8 ) f3
  real ( kind = 8 ) f4
  integer ( kind = 4 ) i
  real ( kind = 8 ) jed
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) m3
  integer ( kind = 4 ) m4
  character ( len = 20 ) s1
  character ( len = 20 ) s2
  character ( len = 20 ) s3
  character ( len = 20 ) s4
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2
  integer ( kind = 4 ) y3
  integer ( kind = 4 ) y4

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST74'
  write ( *, '(a)' ) '  For the Julian calendar:'
  write ( *, '(a)' ) '  YMDF_INC_JULIAN increments a date by days;'
  write ( *, '(a)' ) '  YMDF_NEXT_JULIAN computes the next day,'
  write ( *, '(a)' ) '  YMDF_PREV_JULIAN computes the previous day.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  YMDF date    Tomorrow	    Yesterday	   +10.25 days'
  write ( *, '(a)' ) ' '

  i = 0

  do
 
    i = i + 1
    call jed_test ( i, jed )

    if ( jed < 0 ) then
      exit
    end if

    call jed_to_ymdf_julian ( jed, y1, m1, d1, f1 )

    call ymdf_to_s_julian ( y1, m1, d1, f1, s1 )

    call ymdf_next_julian ( y1, m1, d1, f1, y2, m2, d2, f2 )
    call ymdf_to_s_julian ( y2, m2, d2, f2, s2 )

    call ymdf_prev_julian ( y1, m1, d1, f1, y3, m3, d3, f3 )
    call ymdf_to_s_julian ( y3, m3, d3, f3, s3 )
 
    call ymdf_inc_julian ( y1, m1, d1, f1, days, y4, m4, d4, f4 )
    call ymdf_to_s_julian ( y4, m4, d4, f4, s4 )

    write ( *, '(2x,4a)' ) s1, s2, s3, s4

  end do

  return
end
subroutine test75 ( )

!*****************************************************************************80
!
!! TEST75 tests YMD_INC_YMD_COMMON and YMDF_DIF_YMDF_COMMON.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d1
  integer ( kind = 4 ) d2
  integer ( kind = 4 ) dn1
  integer ( kind = 4 ) dn2
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  real ( kind = 8 ) fn2
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) mn1
  integer ( kind = 4 ) mn2
  character ( len = 20 ) s1
  character ( len = 20 ) s2
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2
  integer ( kind = 4 ) yn1
  integer ( kind = 4 ) yn2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST75'
  write ( *, '(a)' ) '  For the Common calendar,'
  write ( *, '(a)' ) '  YMD_INC_YMD_COMMON increments a YMDF date by YMDF;'
  write ( *, '(a)' ) '  YMDF_DIF_YMDF_COMMON finds the YMDF difference.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    Date1        increment    Date2          difference'
  write ( *, '(a)' ) ' '

  y1 = 1900
  m1 = 5
  d1 = 27
  f1 = 0.0D+00

  yn1 = 50
  mn1 = 9
  dn1 = 10
  fn2 = 0.0D+00

  call ymd_inc_ymd_common ( y1, m1, d1, yn1, mn1, dn1, y2, m2, d2 )

  call ymdf_dif_ymdf_common ( y1, m1, d1, f1, y2, m2, d2, f2, yn2, mn2, dn2, &
    fn2, ierror )

  call ymdf_to_s_common ( y1, m1, d1, f1, s1 )

  call ymdf_to_s_common ( y2, m2, d2, f2, s2 )

  write ( *, '(2x,a,2x,3i3,2x,a,2x,3i3)' ) s1, yn1, mn1, dn1, s2, yn2, mn2, dn2

  return
end
subroutine test76 ( )

!*****************************************************************************80
!
!! TEST76 tests YMDF_TO_WEEKDAY_COMMON.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d1
  real ( kind = 8 ) f1
  integer ( kind = 4 ) i
  real ( kind = 8 ) jed1
  real ( kind = 8 ) jed2
  integer ( kind = 4 ) m1
  character ( len = 20 ) s1
  character ( len = 9 ) s2
  integer ( kind = 4 ) w2
  integer ( kind = 4 ) y1

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST76'
  write ( *, '(a)' ) '  For the Common calendar:'
  write ( *, '(a)' ) '  YMDF_TO_WEEKDAY_COMMON returns the day of the week.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  JED   YMDF     Day of the week'
  write ( *, '(a)' ) ' '

  i = 0

  do

    i = i + 1
    call jed_test ( i, jed1 )

    if ( jed1 < 0 ) then
      exit
    end if

    call jed_to_next_noon ( jed1, jed2 )

    call jed_to_ymdf_common ( jed2, y1, m1, d1, f1 )
 
    call ymdf_to_s_common ( y1, m1, d1, f1, s1 ) 
 
    call ymdf_to_weekday_common ( y1, m1, d1, f1, w2 )
    call weekday_to_name_common ( w2, s2 )

    write ( *, '(f11.2,2x,a,2x,i1,2x,a)' ) jed2, s1, w2, s2

  end do

  return
end
subroutine test77 ( )

!*****************************************************************************80
!
!! TEST77 tests YMDF_TO_WEEKDAY_ENGLISH.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d1
  real ( kind = 8 ) f1
  integer ( kind = 4 ) i
  real ( kind = 8 ) jed1
  real ( kind = 8 ) jed2
  integer ( kind = 4 ) m1
  character ( len = 20 ) s1
  character ( len = 9 ) s2
  integer ( kind = 4 ) w2
  integer ( kind = 4 ) y1

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST77'
  write ( *, '(a)' ) '  For the English calendar:'
  write ( *, '(a)' ) '  YMDF_TO_WEEKDAY_ENGLISH returns the day of the week.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  JED   YMDF    Day of the week'
  write ( *, '(a)' ) ' '

  i = 0

  do

    i = i + 1
    call jed_test ( i, jed1 )

    if ( jed1 < 0 ) then
      exit
    end if

    call jed_to_next_noon ( jed1, jed2 )

    call jed_to_ymdf_english ( jed2, y1, m1, d1, f1 )
 
    call ymdf_to_s_english ( y1, m1, d1, f1, s1 ) 
 
    call ymdf_to_weekday_english ( y1, m1, d1, f1, w2 )
    call weekday_to_name_common ( w2, s2 )

    write ( *, '(f11.2,2x,a,2x,i2,2x,a)' ) jed2, s1, w2, s2

  end do

  return
end
subroutine test775 ( )

!*****************************************************************************80
!
!! TEST775 tests YMDF_TO_WEEKDAY_ENGLISH2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d1
  real ( kind = 8 ) f1
  integer ( kind = 4 ) i
  real ( kind = 8 ) jed1
  real ( kind = 8 ) jed2
  integer ( kind = 4 ) m1
  character ( len = 20 ) s1
  character ( len = 9 ) s2
  integer ( kind = 4 ) w2
  integer ( kind = 4 ) y1

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST775'
  write ( *, '(a)' ) '  For the English calendar:'
  write ( *, '(a)' ) '  YMDF_TO_WEEKDAY_ENGLISH2 returns the day of the week.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Lewis Carroll''s algorithm is used.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  JED   YMDF    Day of the week'
  write ( *, '(a)' ) ' '

  i = 0

  do

    i = i + 1
    call jed_test ( i, jed1 )

    if ( jed1 < 0 ) then
      exit
    end if

    call jed_to_next_noon ( jed1, jed2 )

    call jed_to_ymdf_english ( jed2, y1, m1, d1, f1 )
 
    call ymdf_to_s_english ( y1, m1, d1, f1, s1 ) 
 
    call ymdf_to_weekday_english2 ( y1, m1, d1, f1, w2 )
    call weekday_to_name_common ( w2, s2 )

    write ( *, '(f11.2,2x,a,2x,i2,2x,a)' ) jed2, s1, w2, s2

  end do

  return
end
subroutine test78 ( )

!*****************************************************************************80
!
!! TEST78 tests YMDF_TO_WEEKDAY_GREGORIAN*.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d1
  real ( kind = 8 ) f1
  integer ( kind = 4 ) i
  real ( kind = 8 ) jed1
  real ( kind = 8 ) jed2
  integer ( kind = 4 ) m1
  character ( len = 20 ) s1
  character ( len = 9 ) s2
  character ( len = 9 ) s3
  character ( len = 9 ) s4
  character ( len = 9 ) s5
  character ( len = 9 ) s6
  integer ( kind = 4 ) w2
  integer ( kind = 4 ) w3
  integer ( kind = 4 ) w4
  integer ( kind = 4 ) w5
  integer ( kind = 4 ) w6
  integer ( kind = 4 ) y1

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST78'
  write ( *, '(a)' ) '  For the Gregorian calendar:'
  write ( *, '(a)' ) '  YMDF_TO_WEEKDAY_GREGORIAN,'
  write ( *, '(a)' ) '  YMDF_TO_WEEKDAY_GREGORIAN2,'
  write ( *, '(a)' ) '  YMDF_TO_WEEKDAY_GREGORIAN3,'
  write ( *, '(a)' ) '  YMDF_TO_WEEKDAY_GREGORIAN4, and'
  write ( *, '(a)' ) '  YMDF_TO_WEEKDAY_GREGORIAN5'
  write ( *, '(a)' ) '  return the day of the week.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  (This is "easy" to do for recent dates,'
  write ( *, '(a)' ) '  but look closely at early dates!)'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  JED   YMDF    Day of the week'
  write ( *, '(a)' ) ' '

  i = 0

  do

    i = i + 1
    call jed_test ( i, jed1 )

    if ( jed1 < 0 ) then
      exit
    end if

    call jed_to_next_noon ( jed1, jed2 )

    call jed_to_ymdf_gregorian ( jed2, y1, m1, d1, f1 )
 
    call ymdf_to_s_gregorian ( y1, m1, d1, f1, s1 ) 
 
    call ymdf_to_weekday_gregorian ( y1, m1, d1, f1, w2 )
    call weekday_to_name_common ( w2, s2 )

    call ymdf_to_weekday_gregorian2 ( y1, m1, d1, f1, w3 )
    call weekday_to_name_common ( w3, s3 )

    call ymdf_to_weekday_gregorian3 ( y1, m1, d1, f1, w4 )
    call weekday_to_name_common ( w4, s4 )

    call ymdf_to_weekday_gregorian4 ( y1, m1, d1, f1, w5 )
    call weekday_to_name_common ( w5, s5 )

    call ymdf_to_weekday_gregorian5 ( y1, m1, d1, f1, w6 )
    call weekday_to_name_common ( w6, s6 )

    write ( *, '(f11.2, 2x,a20, 2(2x,a) )' ) jed2, s1, s2, s3
    write ( *, '(  11x, 2x,20x, 3(2x,a) )' )           s4, s5, s6

  end do

  return
end
subroutine test79 ( )

!*****************************************************************************80
!
!! TEST79 tests YMDF_TO_WEEKDAY_HEBREW.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d1
  real ( kind = 8 ) f1
  integer ( kind = 4 ) i
  real ( kind = 8 ) jed1
  real ( kind = 8 ) jed2
  real ( kind = 8 ) jed_epoch
  integer ( kind = 4 ) m1
  character ( len = 20 ) s1
  character ( len = 15 ) s2
  integer ( kind = 4 ) w2
  integer ( kind = 4 ) y1

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST79'
  write ( *, '(a)' ) '  For the HEBREW calendar:'
  write ( *, '(a)' ) '  YMDF_TO_WEEKDAY_HEBREW'
  write ( *, '(a)' ) '  returns the day of the week.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  JED   YMDF           Day of the week'
  write ( *, '(a)' ) ' '

  call epoch_to_jed_hebrew ( jed_epoch )

  i = 0

  do

    i = i + 1
    call jed_test ( i, jed1 )

    if ( jed1 < 0 ) then
      exit
    end if

    if ( jed_epoch <= jed1 ) then

      call jed_to_next_noon ( jed1, jed2 )

      call jed_to_ymdf_hebrew ( jed2, y1, m1, d1, f1 )
 
      call ymdf_to_s_hebrew ( y1, m1, d1, f1, s1 ) 
 
      call ymdf_to_weekday_hebrew ( y1, m1, d1, f1, w2 )

      call weekday_to_name_hebrew ( w2, s2 )

      write ( *, '(f11.2,2x,a,2x,i1,2x,a)' ) jed2, s1, w2, s2

    end if

  end do

  return
end
subroutine test795 ( )

!*****************************************************************************80
!
!! TEST795 tests YMDF_TO_WEEKDAY_ISLAMIC_A.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d1
  real ( kind = 8 ) f1
  integer ( kind = 4 ) i
  real ( kind = 8 ) jed1
  real ( kind = 8 ) jed2
  real ( kind = 8 ) jed_epoch
  integer ( kind = 4 ) m1
  character ( len = 20 ) s1
  character ( len = 15 ) s2
  integer ( kind = 4 ) w2
  integer ( kind = 4 ) y1

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST795'
  write ( *, '(a)' ) '  For the Islamic-A calendar:'
  write ( *, '(a)' ) '  YMDF_TO_WEEKDAY_ISLAMIC_A'
  write ( *, '(a)' ) '  returns the day of the week.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  JED   YMDF           Day of the week'
  write ( *, '(a)' ) ' '

  call epoch_to_jed_islamic_a ( jed_epoch )

  i = 0

  do

    i = i + 1
    call jed_test ( i, jed1 )

    if ( jed1 < 0 ) then
      exit
    end if

    if ( jed_epoch <= jed1 ) then

      call jed_to_next_noon ( jed1, jed2 )

      call jed_to_ymdf_islamic_a ( jed2, y1, m1, d1, f1 )
 
      call ymdf_to_s_islamic ( y1, m1, d1, f1, s1 ) 
 
      call ymdf_to_weekday_islamic_a ( y1, m1, d1, f1, w2 )

      call weekday_to_name_islamic ( w2, s2 )

      write ( *, '(f11.2,2x,a,2x,i1,2x,a)' ) jed2, s1, w2, s2

    end if

  end do

  return
end
subroutine test80 ( )

!*****************************************************************************80
!
!! TEST80 tests YMDF_TO_WEEKDAY_JULIAN, YMDF_TO_WEEKDAY_JULIAN2, YMDF_TO_WEEKDAY_JULIAN3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d1
  real ( kind = 8 ) f1
  integer ( kind = 4 ) i
  real ( kind = 8 ) jed1
  real ( kind = 8 ) jed2
  integer ( kind = 4 ) m1
  character ( len = 20 ) s1
  character ( len = 9 ) s2
  character ( len = 9 ) s3
  character ( len = 9 ) s4
  integer ( kind = 4 ) w2
  integer ( kind = 4 ) w3
  integer ( kind = 4 ) w4
  integer ( kind = 4 ) y1

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST80'
  write ( *, '(a)' ) '  For the Julian calendar,'
  write ( *, '(a)' ) '  YMDF_TO_WEEKDAY_JULIAN,' 
  write ( *, '(a)' ) '  YMDF_TO_WEEKDAY_JULIAN2, and'
  write ( *, '(a)' ) '  YMDF_TO_WEEKDAY_JULIAN3'
  write ( *, '(a)' ) '    return the day of the week of a given date.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  JED   YMDF    Day of the week'
  write ( *, '(a)' ) ' '

  i = 0

  do

    i = i + 1
    call jed_test ( i, jed1 )

    if ( jed1 < 0 ) then
      exit
    end if

    call jed_to_next_noon ( jed1, jed2 )

    call jed_to_ymdf_julian ( jed2, y1, m1, d1, f1 )
 
    call ymdf_to_s_julian ( y1, m1, d1, f1, s1 ) 
 
    call ymdf_to_weekday_julian ( y1, m1, d1, f1, w2 )
    call weekday_to_name_common ( w2, s2 )

    call ymdf_to_weekday_julian2 ( y1, m1, d1, f1, w3 )
    call weekday_to_name_common ( w3, s3 )

    call ymdf_to_weekday_julian2 ( y1, m1, d1, f1, w4 )
    call weekday_to_name_common ( w4, s4 )

    write ( *, '(f11.2,2x,a,2x,a,2x,a,2x,a)' ) &
      jed2, s1, s2, s3, s4

  end do

  return
end
subroutine test805 ( )

!*****************************************************************************80
!
!! TEST805 tests YMDF_TO_WEEKDAY_REPUBLICAN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d1
  real ( kind = 8 ) f1
  integer ( kind = 4 ) i
  real ( kind = 8 ) jed1
  real ( kind = 8 ) jed2
  real ( kind = 8 ) jed_epoch
  integer ( kind = 4 ) m1
  character ( len = 20 ) s1
  character ( len = 9 ) s2
  integer ( kind = 4 ) w2
  integer ( kind = 4 ) y1

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST805'
  write ( *, '(a)' ) '  For the Republican calendar:'
  write ( *, '(a)' ) '  YMDF_TO_WEEKDAY_REPUBLICAN'
  write ( *, '(a)' ) '  returns the day of the week.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    JED   YMDF           Day of the week'
  write ( *, '(a)' ) ' '

  call epoch_to_jed_republican ( jed_epoch )

  i = 0

  do

    i = i + 1
    call jed_test ( i, jed1 )

    if ( jed1 < 0 ) then
      exit
    end if

    if ( jed_epoch <= jed1 ) then

      call jed_to_next_noon ( jed1, jed2 )

      call jed_to_ymdf_republican ( jed2, y1, m1, d1, f1 )
 
      call ymdf_to_s_republican ( y1, m1, d1, f1, s1 ) 
 
      call ymdf_to_weekday_republican ( y1, m1, d1, f1, w2 )

      call weekday_to_name_republican ( w2, s2 )

      write ( *, '(2x,f11.2,2x,a,1x,i2,2x,a)' ) jed2, s1, w2, s2

    end if

  end do

  return
end
subroutine test81 ( )

!*****************************************************************************80
!
!! TEST81 tests YMDF_TO_WEEK_COMMON.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d1
  real ( kind = 8 ) f1
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iweek
  real ( kind = 8 ) jed
  integer ( kind = 4 ) m1
  character ( len = 20 ) s1
  integer ( kind = 4 ) y1

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST81'
  write ( *, '(a)' ) '  YMDF_TO_WEEK_COMMON reports week number for a YMDF date.' 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    JED     YMDF          Week'
  write ( *, '(a)' ) ' '

  i = 0

  do

    i = i + 1
    call jed_test ( i, jed )

    if ( jed < 0 ) then
      exit
    end if

    call jed_to_ymdf_common ( jed, y1, m1, d1, f1 )
 
    call ymdf_to_s_common ( y1, m1, d1, f1, s1 ) 
 
    call ymdf_to_week_common ( y1, m1, d1, f1, iweek )

    write ( *, '(2x,f11.2,2x,a,2x,i2)' ) jed, s1, iweek

  end do

  return
end
subroutine test82 ( )

!*****************************************************************************80
!
!! TEST82 tests YMDHMS_DIF_DHMS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d1
  integer ( kind = 4 ) d2
  integer ( kind = 4 ) days
  integer ( kind = 4 ) h1
  integer ( kind = 4 ) h2
  integer ( kind = 4 ) hours
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) minutes
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) seconds
  integer ( kind = 4 ) second1
  integer ( kind = 4 ) second2
  character ( len = 22 ) s1
  character ( len = 22 ) s2
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2

  y1 = 1997
  m1 = 02
  d1 = 12
  h1 = 13
  n1 = 12
  second1 = 33

  y2 = 1997
  m2 = 03
  d2 = 14
  h2 = 4
  n2 = 21
  second2 = 33
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST82'
  write ( *, '(a)' ) '  YMDHMS_DIF_DHMS finds the DHMS difference'
  write ( *, '(a)' ) '  between YMDHMS dates.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  YMDHMS 1               YMDHMS 2       ' // &
    '            Difference'
  write ( *, '(a)' ) '                                        ' // &
    '            D   H   M   S'
  write ( *, '(a)' ) ' '

  call ymdhms_to_s_common ( y1, m1, d1, h1, n1, second1, s1 )
 
  call ymdhms_to_s_common ( y2, m2, d2, h2, n2, second2, s2 )

  call ymdhms_dif_dhms ( y1, m1, d1, h1, n1, second1, &
    y2, m2, d2, h2, n2, second2, days, hours, minutes, seconds, ierror )
 
  write ( *, '(2x,a,2x,a,2x,4i4)' ) s1, s2, days, hours, minutes, seconds
 
  return
end
subroutine test83 ( )

!*****************************************************************************80
!
!! TEST83 tests YMDHMS_TO_DECIMAL
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 October 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ), parameter :: dhi = 1
  integer ( kind = 4 ), parameter :: dlo = 1
  real ( kind = 8 ) f
  integer ( kind = 4 ) h
  integer ( kind = 4 ), parameter :: hhi = 0
  integer ( kind = 4 ), parameter :: hlo = 0
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) m
  integer ( kind = 4 ), parameter :: mhi = 1
  integer ( kind = 4 ), parameter :: mlo = 1
  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: nhi = 0
  integer ( kind = 4 ), parameter :: nlo = 0
  character ( len = 22 ) s
  integer ( kind = 4 ) ss
  integer ( kind = 4 ), parameter :: shi = 0
  integer ( kind = 4 ), parameter :: slo = 0
  integer ( kind = 4 ) :: seed = 123456789
  integer ( kind = 4 ) y
  real ( kind = 8 ) yf
  integer ( kind = 4 ), parameter :: yhi = 1970
  integer ( kind = 4 ), parameter :: ylo = 1960

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST83'
  write ( *, '(a)' ) '  YMDHMS_TO_DECIMAL converts a date to a year and decimal.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  YMDHMS                         Y.F'
  write ( *, '(a)' ) ' '
 
  do i = 1, 10
 
    call ymdhms_uniform_common ( ylo, mlo, dlo, hlo, nlo, slo, yhi, mhi, dhi, &
      hhi, nhi, shi, seed, y, m, d, h, n, ss )

    call ymdhms_to_s_common ( y, m, d, h, n, ss, s )

    call ymdhms_to_decimal ( y, m, d, h, n, ss, yf )

    write ( *, '(2x,a,5x,f14.4)' ) s, yf

  end do

  return
end
