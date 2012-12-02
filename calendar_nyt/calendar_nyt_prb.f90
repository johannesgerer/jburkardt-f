program main

!*****************************************************************************80
!
!! MAIN is the main program for CALENDAR_NYT_PRB.
!
!  Discussion:
!
!    CALENDAR_NYT_PRB calls a series of tests for CALENDAR_NYT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 February 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CALENDAR_NYT_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the CALENDAR_NYT library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CALENDAR_NYT_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests JED_TO_NYT and NYT_TO_JED.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 December 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) d
  real ( kind = 8 ) diff
  real ( kind = 8 ) f
  integer   ( kind = 4 ) issue2
  real ( kind = 8 ) jed_now
  real ( kind = 8 ) jed_nyt_epoch
  real ( kind = 8 ) jed1
  real ( kind = 8 ) jed3
  integer   ( kind = 4 ) m
  real ( kind = 8 ) r8_uniform
  character ( len = 25 ) s
  integer   ( kind = 4 ) seed
  integer   ( kind = 4 ), parameter :: test_num = 10
  integer   ( kind = 4 ) volume2
  integer   ( kind = 4 ) y

  call epoch_to_jed_nyt ( jed_nyt_epoch )
  call now_to_jed ( jed_now )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  For the New York Times issue date:'
  write ( *, '(a)' ) '  JED_TO_NYT: JED -> NYT.'
  write ( *, '(a)' ) '  NYT_TO_JED: NYT -> JED.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '                           JED (in)      Volume  Issue            JED (out)   Diff'
  write ( *, '(a)' ) ' '

  seed = 123456789

  do i = 1, test_num

    jed1 = r8_uniform ( jed_nyt_epoch, jed_now, seed )

    jed1 = nint ( jed1 ) + 0.5D+00

    call jed_to_ymdf_common ( jed1, y, m, d, f )

    call ymdf_to_s_common ( y, m, d, f, s )

    call jed_to_nyt ( jed1, volume2, issue2 )

    call nyt_to_jed ( volume2, issue2, jed3 )

    diff = jed1 - jed3

    write ( *, '(2x,a25,2x,f11.2,5x,i4,2x,i8,5x,f11.2,2x,f11.2)' ) &
      s, jed1, volume2, issue2, jed3, diff

  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests NYT_TO_JED.
!
!  Discussion:
!
!  Volume  Issue   D M         Y
!  ------  -----  -- --------- ----
!           1705   7 March     1857
!           3407  25 August    1862
!           3794  20 November  1863
!           3804   3 December  1863
!          16579  24 February  1903
!          16909  15 March     1904
!          17251  18 April     1905
!          17561  22 February  1906
!          25320  22 May       1927
!          26243  30 November  1929
!          27538  17 June      1933
!          29033  21 June      1937
!          29807   3 September 1939
!          31545   6 June      1945
!          31972   7 August    1945
!          32984  15 May       1948
!          36074  30 October   1956
!          38910   5 August    1964
!          39342  11 October   1965
!          50939   8 October   1997
!          51599  11 December  2000
!          51874  12 September 2001
!          53108  28 January   2005
!          53715  27 September 2006
!          53960  30 May       2007
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 December 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 40

  integer ( kind = 4 ), dimension ( test_num ) :: d_test = (/ &
    18, 17, 21, 23, 20, &
    19, 16, 26, 13, 22, &
    22,  6,  7, 24, 15, &
    29, 22, 18,  9,  3, &
    22, 23, 14,  8, 15, &
    20, 16, 15, 21, 18, &
     9,  6, 17, 14,  8, &
    31,  1, 11, 28, 22 /)
  integer ( kind = 4 ) d1
  integer ( kind = 4 ) d2
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ), dimension ( test_num ) :: issue_test = (/ &
        1,  2155,  2990,  3432,  3794, &
     4130,  4230,  4576,  5034,  5250, &
     6189, 14499, 15000, 16579, 16909, &
    17292, 17561, 18164, 18856, 21619, &
    24651, 29827, 30000, 31881, 31980, &
    38864, 39317, 40076, 40721, 41418, &
    44027, 44028, 48939, 50000, 50939, &
    51753, 51254, 51599, 53108, 54136 /)
  integer ( kind = 4 ) issue1
  integer ( kind = 4 ) issue3
  real ( kind = 8 ) jed1
  real ( kind = 8 ) jed2
  integer ( kind = 4 ), dimension ( test_num ) :: m_test = (/ &
    9,  8,  4,  9, 11, &
   12,  4,  5, 11,  7, &
    7,  2,  2,  2,  3, &
    5,  2, 10,  9,  4, &
    7,  9,  3,  5,  8, &
    6,  9, 10,  7,  6, &
    8, 11,  4,  3, 10, &
   12,  1, 12,  1, 11 /)
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  character ( len = 25 ) s1
  character ( len = 25 ) s2
  integer ( kind = 4 ) test
  integer ( kind = 4 ), dimension ( test_num ) :: volume_test = (/ &
       1,    7,   10,   11,   13, &
      14,   14,   15,   17,   17, &
      20,   47,   47,   52,   53, &
      54,   55,   57,   58,   66, &
      74,   89,   89,   94,   94, &
     113,  114,  117,  118,  120, &
     127,  128,  141,  144,  147, &
     149,  149,  150,  154,  157 /)
  integer ( kind = 4 ) volume1
  integer ( kind = 4 ) volume3
  integer ( kind = 4 ), dimension ( test_num ) :: y_test = (/ &
    1851, 1858, 1861, 1862, 1863, &
    1864, 1865, 1866, 1867, 1868, &
    1871, 1898, 1898, 1903, 1904, &
    1905, 1906, 1907, 1909, 1917, &
    1925, 1939, 1940, 1945, 1945, &
    1964, 1965, 1967, 1969, 1971, &
    1978, 1978, 1992, 1995, 1997, &
    1999, 2000, 2000, 2005, 2007 /)
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
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
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests JED_TO_NYT.
!
!  Discussion:
!
!    The calculation JED_TO_NYT has been worked on much more than
!    NYT_TO_JED, and we think it is ALMOST correct.
!
!  Volume  Issue   D M         Y
!  ------  -----  -- --------- ----
!       1      1  18 September 1851
!       1     92   3 January   1852
!       2    404   3 January   1853
!       2    416  17 January   1853
!       3    856  15 June      1854
!
!       4   1210   4 August    1855
!       5   1259   1 October   1855
!       5   1491  28 June      1856
!       6   1706   9 March     1857
!       7   2155  17 August    1858
!
!       8   2421  23 June      1859
!       9   2586   4 January   1860
!      10   2897   3 January   1861
!      10   3000   1 May       1861
!      12   3432  23 September 1862
!
!      13   3794  20 November  1863
!      14   4130  19 December  1864
!      14   4230  16 April     1865
!      15   4576  26 May       1866
!      17   5034  13 November  1867
!
!      17   5250  22 July      1868
!      20   6189  22 July      1871
!      47  14499   5 February  1898
!      47  15000   7 February  1898
!      52  16579  24 February  1903
!
!      53  16909  15 March     1904
!      54  17292  29 May       1905
!      55  17561  22 February  1906
!      57  18164  18 October   1907
!      58  18856   9 September 1909
!
!      66  21619   3 April     1917
!      74  24651  22 July      1925
!      89  29827  23 September 1939
!      89  30000  14 March     1940
!      94  31881   8 May       1945
!
!      94  31980  15 August    1945
!     113  38864  20 June      1964
!     114  39317  16 September 1965
!     117  40076  15 October   1967
!     118  40721  21 July      1969
!
!     120  41418  18 June      1971
!     127  44027   9 August    1978
!     128  44028   6 November  1978
!     141  48939  17 April     1992
!     144  50000  14 March     1995
!
!     147  50939   8 October   1997
!     149  51753  31 December  1999
!     149  51254   1 January   2000
!     150  51599  11 December  2000
!     154  53108  28 January   2005
!
!     157  54136  22 November  2007
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 51

  integer ( kind = 4 ) d
  real ( kind = 8 ) f
  integer ( kind = 4 ), dimension ( test_num ) :: issue_test = (/ &
        1,    92,   404,   416,   856, &
     1210,  1259,  1491,  1706,  2155, &
     2421,  2586,  2897,  3000,  3432, &
     3794,  4130,  4230,  4576,  5034, &
     5250,  6189, 14499, 15000, 16579, &
    16909, 17292, 17561, 18164, 18856, &
    21619, 24651, 29827, 30000, 31881, &
    31980, 38864, 39317, 40076, 40721, &
    41418, 44027, 44028, 48939, 50000, &
    50939, 51753, 51254, 51599, 53108, &
    54136 /)
  integer ( kind = 4 ) issue1
  integer ( kind = 4 ) issue2
  real ( kind = 8 ) jed
  real ( kind = 8 ), dimension ( test_num ) :: jed_test = (/ &
   2397383.50, 2397490.50, 2397856.50, 2397870.50, 2398384.50, &
   2398799.50, 2398857.50, 2399128.50, 2399382.50, 2399908.50, &
   2400218.50, 2400413.50, 2400778.50, 2400896.50, 2401406.50, &
   2401829.50, 2402224.50, 2402342.50, 2402747.50, 2403283.50, &
   2403535.50, 2404630.50, 2414325.50, 2414327.50, 2416169.50, &
   2416554.50, 2416994.50, 2417263.50, 2417866.50, 2418558.50, &
   2421321.50, 2424353.50, 2429529.50, 2429702.50, 2431583.50, &
   2431682.50, 2438566.50, 2439019.50, 2439778.50, 2440423.50, &
   2441120.50, 2443729.50, 2443818.50, 2448729.50, 2449790.50, &
   2450729.50, 2451543.50, 2451544.50, 2451889.50, 2453398.50, &
   2454426.50 /)
  integer ( kind = 4 ) m
  character ( len = 25 ) s
  integer ( kind = 4 ) test
  integer ( kind = 4 ), dimension ( test_num ) :: volume_test = (/ &
       1,    1,    2,    2,    3, &
       4,    5,    5,    6,    7, &
       8,    9,   10,   10,   12, &
      13,   14,   14,   15,   17, &
      17,   20,   47,   47,   52, &
      53,   54,   55,   57,   58, &
      66,   74,   89,   89,   94, &
      94,  113,  114,  117,  118, &
     120,  127,  128,  141,  144, &
     147,  149,  149,  150,  154, &
     157 /)
  integer ( kind = 4 ) volume1
  integer ( kind = 4 ) volume2
  integer ( kind = 4 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  For the New York Times issue date:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  JED -> NYT1 by historical record.'
  write ( *, '(a)' ) '  JED -> NYT2 by "JED_TO_NYT"'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  JED         Date   =>      Volume    Issue (lookup)'
  write ( *, '(a)' ) '                     =>      Volume    Issue (compute)'
  write ( *, '(a)' ) '                             Error     Error'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    jed = jed_test(test)

    call jed_to_ymdf_common ( jed, y, m, d, f )
    call ymdf_to_s_common ( y, m, d, f, s )

    issue1 = issue_test(test)
    volume1 = volume_test(test)

    call jed_to_nyt ( jed, volume2, issue2 )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,f11.2,2x,a25,2x,i8,2x,i8)' ) jed, s, volume1, issue1
    write ( *, '(2x, 11x, 2x,25x,2x,i8,2x,i8)' )         volume2, issue2
    write ( *, '(2x, 11x, 2x,25x,2x,i8,2x,i8)' )         volume2 - volume1, issue2 - issue1

  end do

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!!  TEST04 compares JED_TO_NYT and JED_TO_NYT_ISSUE_IDEAL.
!
!  Discussion:
!
!    JED_TO_NYT returns the issue number printed on the New York Times. 
!    There were many "mistakes" and "accidents" and "corrections" in this system.
!
!    JED_TO_NYT_ISSUE_IDEAL returns an "ideal" issue number which keeps track
!    of every issue, in a sensible, usable way.
!
!  Volume  Issue   D M         Y
!  ------  -----  -- --------- ----
!       1      1  18 September 1851
!       1     92   3 January   1852
!       2    404   3 January   1853
!       2    416  17 January   1853
!       3    856  15 June      1854
!
!       4   1210   4 August    1855
!       5   1259   1 October   1855
!       5   1491  28 June      1856
!       6   1706   9 March     1857
!       7   2155  17 August    1858
!
!       8   2421  23 June      1859
!       9   2586   4 January   1860
!      10   2897   3 January   1861
!      10   3000   1 May       1861
!      12   3432  23 September 1862
!
!      13   3794  20 November  1863
!      14   4130  19 December  1864
!      14   4230  16 April     1865
!      15   4576  26 May       1866
!      17   5034  13 November  1867
!
!      17   5250  22 July      1868
!      20   6189  22 July      1871
!      47  14499   5 February  1898
!      47  15000   7 February  1898
!      52  16579  24 February  1903
!
!      53  16909  15 March     1904
!      54  17292  29 May       1905
!      55  17561  22 February  1906
!      57  18164  18 October   1907
!      58  18856   9 September 1909
!
!      66  21619   3 April     1917
!      74  24651  22 July      1925
!      89  29827  23 September 1939
!      89  30000  14 March     1940
!      94  31881   8 May       1945
!
!      94  31980  15 August    1945
!     113  38864  20 June      1964
!     114  39317  16 September 1965
!     117  40076  15 October   1967
!     118  40721  21 July      1969
!
!     120  41418  18 June      1971
!     127  44027   9 August    1978
!     128  44028   6 November  1978
!     141  48939  17 April     1992
!     144  50000  14 March     1995
!
!     147  50939   8 October   1997
!     149  51753  31 December  1999
!     149  51254   1 January   2000
!     150  51599  11 December  2000
!     154  53108  28 January   2005
!
!     157  54136  22 November  2007
!     157  54267  01 April     2008
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 February 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 52

  integer ( kind = 4 ) d
  real ( kind = 8 ) f
  integer ( kind = 4 ), dimension ( test_num ) :: issue_test = (/ &
        1,    92,   404,   416,   856, &
     1210,  1259,  1491,  1706,  2155, &
     2421,  2586,  2897,  3000,  3432, &
     3794,  4130,  4230,  4576,  5034, &
     5250,  6189, 14499, 15000, 16579, &
    16909, 17292, 17561, 18164, 18856, &
    21619, 24651, 29827, 30000, 31881, &
    31980, 38864, 39317, 40076, 40721, &
    41418, 44027, 44028, 48939, 50000, &
    50939, 51753, 51254, 51599, 53108, &
    54136, 54267 /)
  integer ( kind = 4 ) issue1
  integer ( kind = 4 ) issue2
  integer ( kind = 4 ) issue3
  real ( kind = 8 ) jed
  real ( kind = 8 ), dimension ( test_num ) :: jed_test = (/ &
   2397383.50, 2397490.50, 2397856.50, 2397870.50, 2398384.50, &
   2398799.50, 2398857.50, 2399128.50, 2399382.50, 2399908.50, &
   2400218.50, 2400413.50, 2400778.50, 2400896.50, 2401406.50, &
   2401829.50, 2402224.50, 2402342.50, 2402747.50, 2403283.50, &
   2403535.50, 2404630.50, 2414325.50, 2414327.50, 2416169.50, &
   2416554.50, 2416994.50, 2417263.50, 2417866.50, 2418558.50, &
   2421321.50, 2424353.50, 2429529.50, 2429702.50, 2431583.50, &
   2431682.50, 2438566.50, 2439019.50, 2439778.50, 2440423.50, &
   2441120.50, 2443729.50, 2443818.50, 2448729.50, 2449790.50, &
   2450729.50, 2451543.50, 2451544.50, 2451889.50, 2453398.50, &
   2454426.50, 2454557.50 /)
  integer ( kind = 4 ) m
  character ( len = 25 ) s
  integer ( kind = 4 ) test
  integer ( kind = 4 ), dimension ( test_num ) :: volume_test = (/ &
       1,    1,    2,    2,    3, &
       4,    5,    5,    6,    7, &
       8,    9,   10,   10,   12, &
      13,   14,   14,   15,   17, &
      17,   20,   47,   47,   52, &
      53,   54,   55,   57,   58, &
      66,   74,   89,   89,   94, &
      94,  113,  114,  117,  118, &
     120,  127,  128,  141,  144, &
     147,  149,  149,  150,  154, &
     157,  157 /)
  integer ( kind = 4 ) volume1
  integer ( kind = 4 ) volume2
  integer ( kind = 4 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  For the New York Times issue date:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  JED -> NYT1 by historical record.'
  write ( *, '(a)' ) '  JED -> NYT2 by "JED_TO_NYT"'
  write ( *, '(a)' ) '  JED -> NYT# by "JED_TO_NYT_ISSUE_IDEAL"'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       JED                 Date             Volume     Issue     Issue     Issue'
  write ( *, '(a)' ) '                                                     (lookup)  (compute)  (Ideal)'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    jed = jed_test(test)

    call jed_to_ymdf_common ( jed, y, m, d, f )
    call ymdf_to_s_common ( y, m, d, f, s )

    issue1 = issue_test(test)
    volume1 = volume_test(test)

    call jed_to_nyt ( jed, volume2, issue2 )

    call jed_to_nyt_issue_ideal ( jed, issue3 )

    write ( *, '(2x,f11.2,2x,a25,2x,i8,2x,i8,2x,i8,2x,i8)' ) &
      jed, s, volume1, issue1, issue2, issue3

  end do

  return
end
