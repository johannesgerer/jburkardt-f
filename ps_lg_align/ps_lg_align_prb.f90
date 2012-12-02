program main

!*****************************************************************************80
!
!! MAIN is the main program for PS_LG_ALIGN_PRB.
!
!  Discussion:
!
!    PS_LG_ALIGN_PRB runs the profile/sequence alignment tests.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PS_LG_ALIGN_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the PS_LG_ALIGN library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
  call test05 ( )
  call test06 ( )
  call test07 ( )
  call test08 ( )
  call test09 ( )
  call test10 ( )

  call test11 ( )
  call test12 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PS_LG_ALIGN_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests PROFILE_SCORE_READ, PROFILE_SCORE_PRINT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: acid_num = 23
  integer ( kind = 4 ), parameter :: iunit = 1
  integer ( kind = 4 ), parameter :: position_max = 200

  character acid_code(acid_num)
  character conserved(position_max)
  integer ( kind = 4 ) entropy(position_max)
  character ( len = 80 ) file_name
  integer ( kind = 4 ) gap_extend_percent(0:position_max)
  integer ( kind = 4 ) gap_open_percent(0:position_max)
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) position_num
  integer ( kind = 4 ) score(position_max,acid_num)
  integer ( kind = 4 ) score2(acid_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  PROFILE_SCORE_READ reads profile scoring '
  write ( *, '(a)' ) '    information.'
  write ( *, '(a)' ) '  PROFILE_SCORE_PRINT prints the information.'

  file_name = 'profile.txt'
!
!  Read the data from the file.
!
  open ( unit = iunit, file = file_name, status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST01 - Fatal error!'
    write ( *, '(a)' ) '  Could not open the profile file.'
    return
  end if

  call profile_score_read ( acid_code, acid_num, conserved, entropy, &
    gap_open_percent, gap_extend_percent, iunit, position_max, position_num, &
    score, score2 )

  close ( unit = iunit )
!
!  Print the data.
!
  call profile_score_print ( acid_code, acid_num, conserved, entropy, &
    gap_open_percent, gap_extend_percent, position_max, position_num, &
    score, score2 )

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests PS_LG_FSL, PS_LG_FSQ.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 12
  integer ( kind = 4 ), parameter :: n = 12
  integer ( kind = 4 ), parameter :: lds = m

  character b(n)
  real base
  character ( len = n ) bpack
  real e(0:lds,0:n)
  real ev(0:n)
  real f(0:lds,0:n)
  real fv(0:n)
  real, parameter :: gap_extend = -5.0
  integer ( kind = 4 ) gap_extend_percent(0:m)
  real, parameter :: gap_open = -10.0
  integer ( kind = 4 ) gap_open_percent(0:m)
  real ge(0:m)
  real go(0:m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  real, external :: ps_score_12
  real s(0:lds,0:n)
  real sv(0:n)
  integer ( kind = 4 ) t(0:lds,0:n)
  integer ( kind = 4 ) tv(0:n)

  bpack = 'CATGXCCTGCTT'

  do i = 1, n
    b(i) = bpack(i:i)
  end do

  gap_open_percent( 0) = 100
  gap_open_percent( 1) =  95
  gap_open_percent( 2) =  95
  gap_open_percent( 3) =  80
  gap_open_percent( 4) =  60
  gap_open_percent( 5) =  40
  gap_open_percent( 6) =  40
  gap_open_percent( 7) =  40
  gap_open_percent( 8) =  40
  gap_open_percent( 9) =  40
  gap_open_percent(10) =  75
  gap_open_percent(11) =  90
  gap_open_percent(12) = 100

  gap_extend_percent( 0) = 100
  gap_extend_percent( 1) =  90
  gap_extend_percent( 2) =  80
  gap_extend_percent( 3) =  70
  gap_extend_percent( 4) =  60
  gap_extend_percent( 5) =  50
  gap_extend_percent( 6) =  50
  gap_extend_percent( 7) =  50
  gap_extend_percent( 8) =  80
  gap_extend_percent( 9) =  85
  gap_extend_percent(10) =  90
  gap_extend_percent(11) =  95
  gap_extend_percent(12) = 100

  go(0:m) = gap_open * real ( gap_open_percent(0:m) ) / 100.0
  ge(0:m) = gap_extend * real ( gap_extend_percent(0:m) ) / 100.0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02:'
  write ( *, '(a)' ) '  PS_LG_FSQ - Forward score quadratic;'
  write ( *, '(a)' ) '  PS_LG_FSL - Forward score linear;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  GAP_OPEN penalty =   ', gap_open
  write ( *, '(a,g14.6)' ) '  GAP_EXTEND penalty = ', gap_extend
  write ( *, '(a)' ) '  Matching scores by PS_SCORE_12.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Verify that the FSQ and FSL tables agree.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' I  B'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(i3,2x,a1)' ) i, b(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Matching Scores:'
  write ( *, '(a)' ) ' '
  write ( *, '(3x,13i5)' ) ( j, j = 0, n )
  write ( *, '(3x,5x,12(4x,a1))' ) b(1:n)
  write ( *, '(a)' ) ' '
  do i = 0, m
    if ( i == 0 ) then
      write ( *, '(i3,13f5.1)' ) i, ( 0.0, j = 0, n )
    else
      write ( *, '(i3,13f5.1)' ) i, 0.0, ( ps_score_12 ( i, b(j) ), j = 1, n )
    end if
  end do

  m1 = 0
  m2 = m
  n1 = 0
  n2 = n
  base = gap_open

  call ps_lg_fsq ( b, m, m1, m2, n, n1, n2, ps_score_12, go, &
    ge, base, lds, s, e, f, t )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  PS_LG_FSQ:'
  write ( *, '(a)' ) ' '
  write ( *, '(3x,13i5)' ) ( j, j = n1, n2 )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SF:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13f5.1)' ) i, s(i,n1:n2)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'EF:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13f5.1)' ) i, e(i,n1:n2)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FF:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13f5.1)' ) i, f(i,n1:n2)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TF:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13i5)' ) i, t(i,n1:n2)
  end do

  m1 = 0
  n1 = 0
  n2 = n
  base = gap_open

  do m2 = m1, m

    call ps_lg_fsl ( b, m, m1, m2, n, n1, n2, ps_score_12, go, &
      ge, base, sv, ev, fv, tv )

    s(m2,0:n) = sv(0:n)
    e(m2,0:n) = ev(0:n)
    f(m2,0:n) = fv(0:n)
    t(m2,0:n) = tv(0:n)

  end do

  m2 = m

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  PS_LG_FSL:'
  write ( *, '(a)' ) ' '
  write ( *, '(3x,13i5)' ) ( j, j = n1, n2 )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SF:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13f5.1)' ) i, s(i,n1:n2)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'EF:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13f5.1)' ) i, e(i,n1:n2)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FF:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13f5.1)' ) i, f(i,n1:n2)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TF:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13i5)' ) i, t(i,n1:n2)
  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests PS_LG_BSL, PS_LG_BSQ.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 12
  integer ( kind = 4 ), parameter :: n = 12
  integer ( kind = 4 ), parameter :: lds = m

  character b(n)
  real base
  character ( len = n ) bpack
  real e(0:lds,0:n)
  real ev(0:n)
  real f(0:lds,0:n)
  real fv(0:n)
  real, parameter :: gap_extend = -5.0
  integer ( kind = 4 ) gap_extend_percent(0:m)
  real, parameter :: gap_open = -10.0
  integer ( kind = 4 ) gap_open_percent(0:m)
  real ge(0:m)
  real go(0:m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  real, external :: ps_score_12
  real s(0:lds,0:n)
  real sv(0:n)
  integer ( kind = 4 ) t(0:lds,0:n)
  integer ( kind = 4 ) tv(0:n)

  bpack = 'CATGXCCTGCTT'

  do i = 1, n
    b(i) = bpack(i:i)
  end do

  gap_open_percent( 0) = 100
  gap_open_percent( 1) =  95
  gap_open_percent( 2) =  95
  gap_open_percent( 3) =  80
  gap_open_percent( 4) =  60
  gap_open_percent( 5) =  40
  gap_open_percent( 6) =  40
  gap_open_percent( 7) =  40
  gap_open_percent( 8) =  40
  gap_open_percent( 9) =  40
  gap_open_percent(10) =  75
  gap_open_percent(11) =  90
  gap_open_percent(12) = 100

  gap_extend_percent( 0) = 100
  gap_extend_percent( 1) =  90
  gap_extend_percent( 2) =  80
  gap_extend_percent( 3) =  70
  gap_extend_percent( 4) =  60
  gap_extend_percent( 5) =  50
  gap_extend_percent( 6) =  50
  gap_extend_percent( 7) =  50
  gap_extend_percent( 8) =  80
  gap_extend_percent( 9) =  85
  gap_extend_percent(10) =  90
  gap_extend_percent(11) =  95
  gap_extend_percent(12) = 100

  go(0:m) = gap_open * real ( gap_open_percent(0:m) ) / 100.0
  ge(0:m) = gap_extend * real ( gap_extend_percent(0:m) ) / 100.0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03:'
  write ( *, '(a)' ) '  PS_LG_BSQ - Backward score quadratic;'
  write ( *, '(a)' ) '  PS_LG_BSL - Backward score linear;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  GAP_OPEN penalty =   ', gap_open
  write ( *, '(a,g14.6)' ) '  GAP_EXTEND penalty = ', gap_extend
  write ( *, '(a)' ) '  Matching scores by PS_SCORE_12.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Compare with the FSQ and FSL results.'
  write ( *, '(a)' ) '  (The maxima of the SF and SB tables should agree.)'
  write ( *, '(a)' ) '  Verify that the BSQ and BSL tables agree.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' I  B'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(i3,2x,a1)' ) i, b(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Matching Scores:'
  write ( *, '(a)' ) ' '
  write ( *, '(3x,13i5)' ) ( j, j = 0, n )
  write ( *, '(3x,5x,12(4x,a1))' ) b(1:n)
  write ( *, '(a)' ) ' '
  do i = 0, m
    if ( i == 0 ) then
      write ( *, '(i3,13f5.1)' ) i, ( 0.0, j = 0, n )
    else
      write ( *, '(i3,13f5.1)' ) i, 0.0, &
        ( ps_score_12 ( i, b(j) ), j = 1, n )
    end if
  end do

  m1 = 0
  m2 = m
  n1 = 0
  n2 = n
  base = gap_open

  call ps_lg_bsq ( b, m, m1, m2, n, n1, n2, ps_score_12, go, &
    ge, base, lds, s, e, f, t )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  PS_LG_BSQ:'
  write ( *, '(a)' ) ' '
  write ( *, '(3x,13i5)' ) ( j, j = n1, n2 )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SB:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13f5.1)' ) i, s(i,n1:n2)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'EB:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13f5.1)' ) i, e(i,n1:n2)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FB:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13f5.1)' ) i, f(i,n1:n2)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TB:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13i5)' ) i, t(i,n1:n2)
  end do

  m1 = 0
  n1 = 0
  n2 = n
  base = gap_open

  do m1 = 0, m2

    call ps_lg_bsl ( b, m, m1, m2, n, n1, n2, ps_score_12, go, &
      ge, base, sv, ev, fv, tv )

    s(m2,0:n) = sv(0:n)
    e(m2,0:n) = ev(0:n)
    f(m2,0:n) = fv(0:n)
    t(m2,0:n) = tv(0:n)

  end do

  m1 = 0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  PS_LG_BSL:'
  write ( *, '(a)' ) ' '
  write ( *, '(3x,13i5)' ) ( j, j = n1, n2 )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SB:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13f5.1)' ) i, s(i,n1:n2)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'EB:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13f5.1)' ) i, e(i,n1:n2)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FB:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13f5.1)' ) i, f(i,n1:n2)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TB:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13i5)' ) i, t(i,n1:n2)
  end do

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests PS_LG_BPQ, PS_LG_BSQ, PS_LG_FPQ, PS_LG_FSQ.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 12
  integer ( kind = 4 ), parameter :: n = 12
  integer ( kind = 4 ), parameter :: lds = m

  character b(n)
  real base
  character ( len = n ) bpack
  real e(0:lds,0:n)
  real f(0:lds,0:n)
  real, parameter :: gap_extend = -5.0
  integer ( kind = 4 ) gap_extend_percent(0:m)
  real, parameter :: gap_open = -10.0
  integer ( kind = 4 ) gap_open_percent(0:m)
  real ge(0:m)
  real go(0:m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) npath
  integer ( kind = 4 ) pathi(m+n+1)
  integer ( kind = 4 ) pathj(m+n+1)
  real, external :: ps_score_12
  real s(0:lds,0:n)
  integer ( kind = 4 ) t(0:lds,0:n)

  bpack = 'CATGXCCTGCTT'

  do i = 1, n
    b(i) = bpack(i:i)
  end do

  gap_open_percent( 0) = 100
  gap_open_percent( 1) =  95
  gap_open_percent( 2) =  95
  gap_open_percent( 3) =  80
  gap_open_percent( 4) =  60
  gap_open_percent( 5) =  40
  gap_open_percent( 6) =  40
  gap_open_percent( 7) =  40
  gap_open_percent( 8) =  40
  gap_open_percent( 9) =  40
  gap_open_percent(10) =  75
  gap_open_percent(11) =  90
  gap_open_percent(12) = 100

  gap_extend_percent( 0) = 100
  gap_extend_percent( 1) =  90
  gap_extend_percent( 2) =  80
  gap_extend_percent( 3) =  70
  gap_extend_percent( 4) =  60
  gap_extend_percent( 5) =  50
  gap_extend_percent( 6) =  50
  gap_extend_percent( 7) =  50
  gap_extend_percent( 8) =  80
  gap_extend_percent( 9) =  85
  gap_extend_percent(10) =  90
  gap_extend_percent(11) =  95
  gap_extend_percent(12) = 100

  go(0:m) = gap_open * real ( gap_open_percent(0:m) ) / 100.0
  ge(0:m) = gap_extend * real ( gap_extend_percent(0:m) ) / 100.0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04:'
  write ( *, '(a)' ) '  PS_LG_BPQ - Backward path quadratic;'
  write ( *, '(a)' ) '  PS_LG_BSQ - Backward score quadratic;'
  write ( *, '(a)' ) '  PS_LG_FPQ - Forward path quadratic;'
  write ( *, '(a)' ) '  PS_LG_FSQ - Forward score quadratic;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  GAP_OPEN penalty =   ', gap_open
  write ( *, '(a,g14.6)' ) '  GAP_EXTEND penalty = ', gap_extend
  write ( *, '(a)' ) '  Matching scores by PS_SCORE_12.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Verify that the FSQ and BSQ paths agree.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' I  B'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(i3,2x,a1)' ) i, b(i)
  end do

  m1 = 0
  m2 = m
  n1 = 0
  n2 = n
  base = gap_open

  call ps_lg_fsq ( b, m, m1, m2, n, n1, n2, ps_score_12, go, &
    ge, base, lds, s, e, f, t )

  call r4mat_imax ( lds+1, m+1, n+1, s, i2, j2 )

  i2 = i2 - 1
  j2 = j2 - 1

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PS_LG_FSQ:'
  write ( *, '(a,2i6)' ) '  Optimal score indices:', i2, j2
  write ( *, '(a,g14.6)' ) '  Optimal score: ', s(i2,j2)

  call ps_lg_fpq ( i2, j2, m, m1, m2, n, n1, n2, lds, t, npath, pathi, pathj )

  call i4vec2_print ( npath, pathi, pathj, '  FSQ Matching path:' )

  call ps_lg_match_print ( b, m, n, npath, pathi, pathj, &
    ps_score_12, go, ge )

  m1 = 0
  m2 = m
  n1 = 0
  n2 = n
  base = gap_open

  call ps_lg_bsq ( b, m, m1, m2, n, n1, n2, ps_score_12, go, &
    ge, base, lds, s, e, f, t )

  call r4mat_imax ( lds+1, m+1, n+1, s, i1, j1 )

  i1 = i1 - 1
  j1 = j1 - 1

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PS_LG_BSQ:'
  write ( *, '(a,2i6)' ) '  Optimal score indices:', i1, j1
  write ( *, '(a,g14.6)' ) '  Optimal score: ', s(i1,j1)

  call ps_lg_bpq ( i1, j1, m, m1, m2, n, n1, n2, lds, t, npath, pathi, pathj )

  call i4vec2_print ( npath, pathi, pathj, '  BSQ Matching path:' )

  call ps_lg_match_print ( b, m, n, npath, pathi, pathj, ps_score_12, go, ge )

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 tests PS_LG_CORNERS, PS_LG_RPL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 12
  integer ( kind = 4 ), parameter :: n = 12
  integer ( kind = 4 ), parameter :: lds = m

  character b(n)
  character ( len = n ) bpack
  real, parameter :: gap_extend = -5.0
  integer ( kind = 4 ) gap_extend_percent(0:m)
  real, parameter :: gap_open = -10.0
  integer ( kind = 4 ) gap_open_percent(0:m)
  real ge(0:m)
  real go(0:m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) npath
  integer ( kind = 4 ) pathi(m+n+1)
  integer ( kind = 4 ) pathj(m+n+1)
  real, external :: ps_score_12

  bpack = 'CATGXCCTGCTT'

  do i = 1, n
    b(i) = bpack(i:i)
  end do

  gap_open_percent( 0) = 100
  gap_open_percent( 1) =  95
  gap_open_percent( 2) =  95
  gap_open_percent( 3) =  80
  gap_open_percent( 4) =  60
  gap_open_percent( 5) =  40
  gap_open_percent( 6) =  40
  gap_open_percent( 7) =  40
  gap_open_percent( 8) =  40
  gap_open_percent( 9) =  40
  gap_open_percent(10) =  75
  gap_open_percent(11) =  90
  gap_open_percent(12) = 100

  gap_extend_percent( 0) = 100
  gap_extend_percent( 1) =  90
  gap_extend_percent( 2) =  80
  gap_extend_percent( 3) =  70
  gap_extend_percent( 4) =  60
  gap_extend_percent( 5) =  50
  gap_extend_percent( 6) =  50
  gap_extend_percent( 7) =  50
  gap_extend_percent( 8) =  80
  gap_extend_percent( 9) =  85
  gap_extend_percent(10) =  90
  gap_extend_percent(11) =  95
  gap_extend_percent(12) = 100

  go(0:m) = gap_open * real ( gap_open_percent(0:m) ) / 100.0
  ge(0:m) = gap_extend * real ( gap_extend_percent(0:m) ) / 100.0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05:'
  write ( *, '(a)' ) '  PS_LG_CORNERS - identifies the endpoints of an'
  write ( *, '(a)' ) '	 optimal local alignment.'
  write ( *, '(a)' ) '  PS_LG_RPL - Recursive path linear;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  GAP_OPEN penalty =   ', gap_open
  write ( *, '(a,g14.6)' ) '  GAP_EXTEND penalty = ', gap_extend
  write ( *, '(a)' ) '  Matching scores by PS_SCORE_12.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Verify that the path agrees with the FSQ and BSQ paths.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' I  B'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(i3,2x,a1)' ) i, b(i)
  end do
!
!  Identify the "corners" of the optimal local alignment.
!
  m1 = 0
  m2 = m
  n1 = 0
  n2 = n

  call ps_lg_corners ( b, m, m1, m2, n, n1, n2, ps_score_12, go, ge, &
    i1, j1, i2, j2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PS_LG_CORNERS:'
  write ( *, '(a)' ) '  Optimal alignment "corners" are:'
  write ( *, '(a,2i6)' ) '    I1, J1 = ', i1, j1
  write ( *, '(a,2i6)' ) '    I2, J2 = ', i2, j2

  call ps_lg_rpl ( b, m, i1, i2, n, j1, j2, ps_score_12, go, &
    ge, npath, pathi, pathj )

  call i4vec2_print ( npath, pathi, pathj, '  RPL Matching path:' )

  call ps_lg_match_print ( b, m, n, npath, pathi, pathj, ps_score_12, go, ge )

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 tests PS_LG_FPQ, PS_LG_FSQ.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 4
  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: lds = m

  character b(n)
  real base
  real e(0:lds,0:n)
  real f(0:lds,0:n)
  real, parameter :: gap_extend = -0.5
  integer ( kind = 4 ) gap_extend_percent(0:m)
  real, parameter :: gap_open = -2.0
  integer ( kind = 4 ) gap_open_percent(0:m)
  real ge(0:m)
  real go(0:m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) npath
  integer ( kind = 4 ) pathi(m+n+1)
  integer ( kind = 4 ) pathj(m+n+1)
  real, external :: ps_score_simple
  real s(0:lds,0:n)
  integer ( kind = 4 ) t(0:lds,0:n)

  b(1) = 'A'
  b(2) = 'C'
  b(3) = 'G'
  b(4) = 'C'
  b(5) = 'T'

  gap_open_percent(0) =  95
  gap_open_percent(1) =  90
  gap_open_percent(2) =  75
  gap_open_percent(3) =  50
  gap_open_percent(4) =  25

  gap_extend_percent(0) = 100
  gap_extend_percent(1) =  80
  gap_extend_percent(2) =  70
  gap_extend_percent(3) =  40
  gap_extend_percent(4) =  85

  go(0:m) = gap_open * real ( gap_open_percent(0:m) ) / 100.0
  ge(0:m) = gap_extend * real ( gap_extend_percent(0:m) ) / 100.0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06:'
  write ( *, '(a)' ) '  Forward algorithms using quadratic space:'
  write ( *, '(a)' ) '  PS_LG_FSQ determines an alignment score;'
  write ( *, '(a)' ) '  PS_LG_FPQ determines an alignment path.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' I  B'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(i3,2x,a1)' ) i, b(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Matching Scores:'
  write ( *, '(a)' ) ' '
  write ( *, '(3x,13i5)' ) ( j, j = 0, n )
  write ( *, '(3x,5x,12(4x,a1))' ) b(1:n)
  write ( *, '(a)' ) ' '
  do i = 0, m
    if ( i == 0 ) then
      write ( *, '(i3,13f5.1)' ) i, ( 0.0, j = 0, n )
    else
      write ( *, '(i3,13f5.1)' ) i, 0.0, &
        ( ps_score_simple ( i, b(j) ), j = 1, n )
    end if
  end do

  m1 = 0
  m2 = m
  n1 = 0
  n2 = n
  base = 0.0

  call ps_lg_fsq ( b, m, m1, m2, n, n1, n2, ps_score_simple, go, ge, &
    base, lds, s, e, f, t )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  PS_LG_FSQ:'
  write ( *, '(a)' ) ' '
  write ( *, '(3x,13i5)' ) ( j, j = n1, n2 )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SF:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13f5.1)' ) i, s(i,n1:n2)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'EF:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13f5.1)' ) i, e(i,n1:n2)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FF:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13f5.1)' ) i, f(i,n1:n2)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TF:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13i5)' ) i, t(i,n1:n2)
  end do

  call r4mat_imax ( lds+1, m+1, n+1, s, i2, j2 )

  i2 = i2 - 1
  j2 = j2 - 1

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PS_LG_FSQ:'
  write ( *, '(a,2i6)' ) '  Optimal score indices:', i2, j2
  write ( *, '(a,g14.6)' ) '  Optimal score: ', s(i2,j2)

  call ps_lg_fpq ( i2, j2, m, m1, m2, n, n1, n2, lds, t, npath, pathi, pathj )

  call i4vec2_print ( npath, pathi, pathj, '  FSQ Matching path:' )

  call ps_lg_match_print ( b, m, n, npath, pathi, pathj, &
    ps_score_simple, go, ge )

  return
end
subroutine test07 ( )

!*****************************************************************************80
!
!! TEST07 tests PS_LG_BPQ, PS_LG_BSQ.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 4
  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: lds = m

  character b(n)
  real base
  real e(0:lds,0:n)
  real f(0:lds,0:n)
  real, parameter :: gap_extend = -0.5
  integer ( kind = 4 ) gap_extend_percent(0:m)
  real, parameter :: gap_open = -2.0
  integer ( kind = 4 ) gap_open_percent(0:m)
  real ge(0:m)
  real go(0:m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) npath
  integer ( kind = 4 ) pathi(m+n+1)
  integer ( kind = 4 ) pathj(m+n+1)
  real, external :: ps_score_simple
  real s(0:lds,0:n)
  integer ( kind = 4 ) t(0:lds,0:n)

  b(1) = 'A'
  b(2) = 'C'
  b(3) = 'G'
  b(4) = 'C'
  b(5) = 'T'

  gap_open_percent(0) =  95
  gap_open_percent(1) =  90
  gap_open_percent(2) =  75
  gap_open_percent(3) =  50
  gap_open_percent(4) =  25

  gap_extend_percent(0) = 100
  gap_extend_percent(1) =  80
  gap_extend_percent(2) =  70
  gap_extend_percent(3) =  40
  gap_extend_percent(4) =  85

  go(0:m) = gap_open * real ( gap_open_percent(0:m) ) / 100.0
  ge(0:m) = gap_extend * real ( gap_extend_percent(0:m) ) / 100.0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07:'
  write ( *, '(a)' ) '  Backward algorithms using quadratic space:'
  write ( *, '(a)' ) '  PS_LG_BSQ determines an alignment score;'
  write ( *, '(a)' ) '  PS_LG_BPQ determines an alignment path.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Verify that these results match those for'
  write ( *, '(a)' ) '  PS_LG_FSQ/PS_GG_FPQ on the same data.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' I  B'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(i3,2x,a1)' ) i, b(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Matching Scores:'
  write ( *, '(a)' ) ' '
  write ( *, '(3x,13i5)' ) ( j, j = 0, n )
  write ( *, '(3x,5x,12(4x,a1))' ) b(1:n)
  write ( *, '(a)' ) ' '
  do i = 0, m
    if ( i == 0 ) then
      write ( *, '(i3,13f5.1)' ) i, ( 0.0, j = 0, n )
    else
      write ( *, '(i3,13f5.1)' ) i, 0.0, &
        ( ps_score_simple ( i, b(j) ), j = 1, n )
    end if
  end do

  m1 = 0
  m2 = m
  n1 = 0
  n2 = n
  base = 0.0

  call ps_lg_bsq ( b, m, m1, m2, n, n1, n2, ps_score_simple, go, ge, &
    base, lds, s, e, f, t )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  PS_LG_BSQ:'
  write ( *, '(a)' ) ' '
  write ( *, '(3x,13i5)' ) ( j, j = n1, n2 )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SB:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13f5.1)' ) i, s(i,n1:n2)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'EB:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13f5.1)' ) i, e(i,n1:n2)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FB:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13f5.1)' ) i, f(i,n1:n2)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TB:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13i5)' ) i, t(i,n1:n2)
  end do

  call r4mat_imax ( lds+1, m+1, n+1, s, i1, j1 )

  i1 = i1 - 1
  j1 = j1 - 1

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PS_LG_BSQ:'
  write ( *, '(a,2i6)' ) '  Optimal score indices:', i1, j1
  write ( *, '(a,g14.6)' ) '  Optimal score: ', s(i1,j1)

  call ps_lg_bpq ( i1, j1, m, m1, m2, n, n1, n2, lds, t, npath, pathi, pathj )

  call i4vec2_print ( npath, pathi, pathj, '  BSQ Matching path:' )

  call ps_lg_match_print ( b, m, n, npath, pathi, pathj, ps_score_simple, &
    go, ge )

  return
end
subroutine test08 ( )

!*****************************************************************************80
!
!! TEST08 tests PS_LG_CORNERS, PS_LG_RPL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 4
  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: lds = m

  character b(n)
  real base
  real, parameter :: gap_extend = -0.5
  integer ( kind = 4 ) gap_extend_percent(0:m)
  real, parameter :: gap_open = -2.0
  integer ( kind = 4 ) gap_open_percent(0:m)
  real ge(0:m)
  real go(0:m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) npath
  integer ( kind = 4 ) pathi(m+n+1)
  integer ( kind = 4 ) pathj(m+n+1)
  real, external :: ps_score_simple

  b(1) = 'A'
  b(2) = 'C'
  b(3) = 'G'
  b(4) = 'C'
  b(5) = 'T'

  gap_open_percent(0) =  95
  gap_open_percent(1) =  90
  gap_open_percent(2) =  75
  gap_open_percent(3) =  50
  gap_open_percent(4) =  25

  gap_extend_percent(0) = 100
  gap_extend_percent(1) =  80
  gap_extend_percent(2) =  70
  gap_extend_percent(3) =  40
  gap_extend_percent(4) =  85

  go(0:m) = gap_open * real ( gap_open_percent(0:m) ) / 100.0
  ge(0:m) = gap_extend * real ( gap_extend_percent(0:m) ) / 100.0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08:'
  write ( *, '(a)' ) '  PS_LG_CORNERS determines determines the endpoints'
  write ( *, '(a)' ) '    of an optimal local alignment;'
  write ( *, '(a)' ) '  PS_LG_RPL determines an alignment path.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Verify that these results match those for'
  write ( *, '(a)' ) '  PS_LG_FSQ/PS_GG_FPQ and PS_GG_BSQ/PS_GG_BPQ'
  write ( *, '(a)' ) '  on the same data.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' I  B'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(i3,2x,a1)' ) i, b(i)
  end do
!
!  Identify the "corners" of the optimal local alignment.
!
  m1 = 0
  m2 = m
  n1 = 0
  n2 = n

  call ps_lg_corners ( b, m, m1, m2, n, n1, n2, ps_score_simple, go, ge, &
    i1, j1, i2, j2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PS_LG_CORNERS:'
  write ( *, '(a)' ) '  Optimal alignment "corners" are:'
  write ( *, '(a,2i6)' ) '    I1, J1 = ', i1, j1
  write ( *, '(a,2i6)' ) '    I2, J2 = ', i2, j2

  call ps_lg_rpl ( b, m, i1, i2, n, j1, j2, ps_score_simple, go, ge, &
    npath, pathi, pathj )

  call i4vec2_print ( npath, pathi, pathj, '  RPL Matching path:' )

  call ps_lg_match_print ( b, m, n, npath, pathi, pathj, ps_score_simple, &
    go, ge )

  return
end
subroutine test09 ( )

!*****************************************************************************80
!
!! TEST09 tests PS_LG_BSL, PS_LG_RPL, PS_LG_CORNERS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 4
  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: lds = m

  character b(n)
  real base
  real e(0:lds,0:n)
  real ev(0:n)
  real f(0:lds,0:n)
  real fv(0:n)
  real, parameter :: gap_extend = -0.5
  integer ( kind = 4 ) gap_extend_percent(0:m)
  real, parameter :: gap_open = -2.0
  integer ( kind = 4 ) gap_open_percent(0:m)
  real ge(0:m)
  real go(0:m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) npath
  integer ( kind = 4 ) pathi(m+n+1)
  integer ( kind = 4 ) pathj(m+n+1)
  real, external :: ps_score_simple
  real s(0:lds,0:n)
  real sv(0:n)
  integer ( kind = 4 ) t(0:lds,0:n)
  integer ( kind = 4 ) tv(0:n)

  b(1) = 'A'
  b(2) = 'C'
  b(3) = 'G'
  b(4) = 'C'
  b(5) = 'T'

  gap_open_percent(0) =  95
  gap_open_percent(1) =  90
  gap_open_percent(2) =  75
  gap_open_percent(3) =  50
  gap_open_percent(4) =  25

  gap_extend_percent(0) = 100
  gap_extend_percent(1) =  80
  gap_extend_percent(2) =  70
  gap_extend_percent(3) =  40
  gap_extend_percent(4) =  85

  go(0:m) = gap_open * real ( gap_open_percent(0:m) ) / 100.0
  ge(0:m) = gap_extend * real ( gap_extend_percent(0:m) ) / 100.0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST09:'
  write ( *, '(a)' ) '  Backward algorithm using linear space:'
  write ( *, '(a)' ) '  PS_LG_BSL determines an alignment score;'
  write ( *, '(a)' ) '  Recursive algorithm using linear space:'
  write ( *, '(a)' ) '  PS_LG_RPL determines an alignment path.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Verify that these results match those for'
  write ( *, '(a)' ) '  PS_LG_FSQ/PS_GG_FPQ, PS_GG_BSQ/PS_GG_BPQ'
  write ( *, '(a)' ) '  and PS_LG_FSL/PS_GG_RPL on the same data.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' I  B'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(i3,2x,a1)' ) i, b(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Matching Scores:'
  write ( *, '(a)' ) ' '
  write ( *, '(3x,13i5)' ) ( j, j = 0, n )
  write ( *, '(3x,5x,12(4x,a1))' ) b(1:n)
  write ( *, '(a)' ) ' '
  do i = 0, m
    if ( i == 0 ) then
      write ( *, '(i3,13f5.1)' ) i, ( 0.0, j = 0, n )
    else
      write ( *, '(i3,13f5.1)' ) i, 0.0, &
        ( ps_score_simple ( i, b(j) ), j = 1, n )
    end if
  end do

  m1 = 0
  m2 = m
  n1 = 0
  n2 = n
  base = 0.0

  do m1 = 0, m2

    call ps_lg_bsl ( b, m, m1, m2, n, n1, n2, ps_score_simple, go, ge, &
      base, sv, ev, fv, tv )

    s(m1,0:n) = sv(0:n)
    e(m1,0:n) = ev(0:n)
    f(m1,0:n) = fv(0:n)
    t(m1,0:n) = tv(0:n)

  end do

  m1 = 0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  PS_LG_BSL:'
  write ( *, '(a)' ) ' '
  write ( *, '(3x,13i5)' ) ( j, j = n1, n2 )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SB:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13f5.1)' ) i, s(i,n1:n2)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'EB:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13f5.1)' ) i, e(i,n1:n2)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FB:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13f5.1)' ) i, f(i,n1:n2)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TB:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13i5)' ) i, t(i,n1:n2)
  end do

  call ps_lg_corners ( b, m, m1, m2, n, n1, n2, ps_score_simple, go, &
    ge, i1, j1, i2, j2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PS_LG_CORNERS:'
  write ( *, '(a)' ) '  Optimal alignment "corners" are:'
  write ( *, '(a,2i6)' ) '    I1, J1 = ', i1, j1
  write ( *, '(a,2i6)' ) '    I2, J2 = ', i2, j2

  call ps_lg_rpl ( b, m, i1, i2, n, j1, j2, ps_score_simple, go, ge, &
    npath, pathi, pathj )

  call i4vec2_print ( npath, pathi, pathj, '  RPL Matching path:' )

  call ps_lg_match_print ( b, m, n, npath, pathi, pathj, ps_score_simple, &
    go, ge )

  return
end
subroutine test10 ( )

!*****************************************************************************80
!
!! TEST10 tests PS_LG_FPQ, PS_LG_FSQ.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 200
  integer ( kind = 4 ), parameter :: m_max = 200
  integer ( kind = 4 ), parameter :: lds = m_max

  character b(n_max)
  real base
  real e(0:lds,0:n_max)
  real f(0:lds,0:n_max)
  character ( len = 80 ) file_name
  real, parameter :: gap_extend = -0.5
  integer ( kind = 4 ) gap_extend_percent(0:m_max)
  real, parameter :: gap_open = -2.0
  integer ( kind = 4 ) gap_open_percent(0:m_max)
  real ge(0:m_max)
  real go(0:m_max)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) :: iunit = 1
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) npath
  integer ( kind = 4 ) pathi(m_max+n_max+1)
  integer ( kind = 4 ) pathj(m_max+n_max+1)
  real, external :: profile_score
  real s(0:lds,0:n_max)
  character ( len = 100 ) sequence
  character ( len = 20 ) sequence_name
  integer ( kind = 4 ) t(0:lds,0:n_max)
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST10:'
  write ( *, '(a)' ) '  Forward algorithms using quadratic space:'
  write ( *, '(a)' ) '  PS_LG_FSQ determines an alignment score;'
  write ( *, '(a)' ) '  PS_LG_FPQ determines an alignment path.'

  file_name = 'profile.txt'

  open ( unit = iunit, file = file_name, status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST10 - Fatal error!'
    write ( *, '(a)' ) '  Could not open the profile file.'
    return
  end if

  call profile_score_read2 ( gap_open_percent, gap_extend_percent, iunit, &
    m_max, m )

  close ( unit = iunit )

  go(0:m) = gap_open * real ( gap_open_percent(0:m) ) / 100.0
  ge(0:m) = gap_extend * real ( gap_extend_percent(0:m) ) / 100.0

  do test = 1, 4
!
!  Set the sequence name and values.
!
    if ( test == 1 ) then

      sequence_name = 's00657 (chunk)'

      sequence = 'NGQSYRGTYSTTVTGRTCQAWSSMTPHSHS'

    else if ( test == 2 ) then

      sequence_name = 'c61545 (chunk)'

      sequence = 'DADKSPWCYTTDPRVRWEFCNLKK'

    else if ( test == 3 ) then

      sequence_name = 'plhu-k (chunk)'

      sequence = 'TPHRHQKTPENYPNAGLTMNYCRNPDADKG'

    else if ( test == 4 ) then

      sequence_name = 'Made-Up (chunk)'

      sequence = 'WSSMTPHRHQCATKTPENYPNAGLTM'

    end if
!
!  Convert the sequence string to a vector of characters.
!
    n = len_trim ( sequence )

    call s_to_chvec ( sequence, n, b )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Sequence name: ' // trim ( sequence_name )
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Sequence:'
    write ( *, '(a)' ) ' '
    do i = 1, n
      write ( *, '(i5,5x,a1)' ) i, b(i)
    end do

    m1 = 0
    m2 = m
    n1 = 0
    n2 = n
    base = 0.0

    call ps_lg_fsq ( b, m, m1, m2, n, n1, n2, profile_score, go, ge, &
      base, lds, s, e, f, t )

    call r4mat_imax ( lds+1, m+1, n+1, s, i2, j2 )

    i2 = i2 - 1
    j2 = j2 - 1

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_LG_FSQ:'
    write ( *, '(a,2i6)' ) '  Optimal score indices:', i2, j2
    write ( *, '(a,g14.6)' ) '  Optimal score: ', s(i2,j2)

    call ps_lg_fpq ( i2, j2, m, m1, m2, n, n1, n2, lds, t, npath, pathi, pathj )

    call i4vec2_print ( npath, pathi, pathj, '  FSQ Matching path:' )

    call ps_lg_match_print ( b, m, n, npath, pathi, pathj, &
      profile_score, go, ge )

  end do

  return
end
subroutine test11 ( )

!*****************************************************************************80
!
!! TEST11 tests PS_LG_BPQ, PS_LG_BSQ.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 200
  integer ( kind = 4 ), parameter :: m_max = 200
  integer ( kind = 4 ), parameter :: lds = m_max

  character b(n_max)
  real base
  real e(0:lds,0:n_max)
  real f(0:lds,0:n_max)
  character ( len = 80 ) file_name
  real, parameter :: gap_extend = -0.5
  integer ( kind = 4 ) gap_extend_percent(0:m_max)
  real, parameter :: gap_open = -2.0
  integer ( kind = 4 ) gap_open_percent(0:m_max)
  real ge(0:m_max)
  real go(0:m_max)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) :: iunit = 1
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) npath
  integer ( kind = 4 ) pathi(m_max+n_max+1)
  integer ( kind = 4 ) pathj(m_max+n_max+1)
  real, external :: profile_score
  real s(0:lds,0:n_max)
  character ( len = 100 ) sequence
  character ( len = 20 ) sequence_name
  integer ( kind = 4 ) t(0:lds,0:n_max)
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST11:'
  write ( *, '(a)' ) '  Backward algorithms using quadratic space:'
  write ( *, '(a)' ) '  PS_LG_BSQ determines an alignment score;'
  write ( *, '(a)' ) '  PS_LG_BPQ determines an alignment path.'

  file_name = 'profile.txt'

  open ( unit = iunit, file = file_name, status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST11 - Fatal error!'
    write ( *, '(a)' ) '  Could not open the profile file.'
    return
  end if

  call profile_score_read2 ( gap_open_percent, gap_extend_percent, iunit, &
    m_max, m )

  close ( unit = iunit )

  go(0:m) = gap_open * real ( gap_open_percent(0:m) ) / 100.0
  ge(0:m) = gap_extend * real ( gap_extend_percent(0:m) ) / 100.0

  do test = 1, 4
!
!  Set the sequence name and values.
!
    if ( test == 1 ) then

      sequence_name = 's00657 (chunk)'

      sequence = 'NGQSYRGTYSTTVTGRTCQAWSSMTPHSHS'

    else if ( test == 2 ) then

      sequence_name = 'c61545 (chunk)'

      sequence = 'DADKSPWCYTTDPRVRWEFCNLKK'

    else if ( test == 3 ) then

      sequence_name = 'plhu-k (chunk)'

      sequence = 'TPHRHQKTPENYPNAGLTMNYCRNPDADKG'

    else if ( test == 4 ) then

      sequence_name = 'Made-Up (chunk)'

      sequence = 'WSSMTPHRHQCATKTPENYPNAGLTM'

    end if
!
!  Convert the sequence string to a vector of characters.
!
    n = len_trim ( sequence )

    call s_to_chvec ( sequence, n, b )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Sequence name: ' // trim ( sequence_name )
    write ( *, '(a)' ) ' '

    m1 = 0
    m2 = m
    n1 = 0
    n2 = n
    base = 0.0

    call ps_lg_bsq ( b, m, m1, m2, n, n1, n2, profile_score, go, ge, &
      base, lds, s, e, f, t )

    call r4mat_imax ( lds+1, m+1, n+1, s, i1, j1 )

    i1 = i1 - 1
    j1 = j1 - 1

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_LG_BSQ:'
    write ( *, '(a,2i6)' ) '  Optimal score indices:', i1, j1
    write ( *, '(a,g14.6)' ) '  Optimal score: ', s(i1,j1)

    call ps_lg_bpq ( i1, j1, m, m1, m2, n, n1, n2, lds, t, npath, pathi, pathj )

    call i4vec2_print ( npath, pathi, pathj, '  BSQ Matching path:' )

    call ps_lg_match_print ( b, m, n, npath, pathi, pathj, profile_score, &
      go, ge )

  end do

  return
end
subroutine test12 ( )

!*****************************************************************************80
!
!! TEST12 tests PS_LG_CORNERS, PS_LG_RPL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 200
  integer ( kind = 4 ), parameter :: m_max = 200

  character b(n_max)
  character ( len = 80 ) file_name
  real, parameter :: gap_extend = -0.5
  integer ( kind = 4 ) gap_extend_percent(0:m_max)
  real, parameter :: gap_open = -2.0
  integer ( kind = 4 ) gap_open_percent(0:m_max)
  real ge(0:m_max)
  real go(0:m_max)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) :: iunit = 1
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) npath
  integer ( kind = 4 ) pathi(m_max+n_max+1)
  integer ( kind = 4 ) pathj(m_max+n_max+1)
  real, external :: profile_score
  character ( len = 100 ) sequence
  character ( len = 20 ) sequence_name
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST12:'
  write ( *, '(a)' ) '  Linear space local alignment algorithms.'
  write ( *, '(a)' ) '  PS_LG_CORNERS determines the "corners" of an alignment;'
  write ( *, '(a)' ) '  PS_LG_RPL determines the path of an alignment.'
  write ( *, '(a)' ) ' '

  file_name = 'profile.txt'

  open ( unit = iunit, file = file_name, status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST12 - Fatal error!'
    write ( *, '(a)' ) '  Could not open the profile file.'
    return
  end if

  call profile_score_read2 ( gap_open_percent, gap_extend_percent, iunit, &
    m_max, m )

  close ( unit = iunit )

  go(0:m) = gap_open * real ( gap_open_percent(0:m) ) / 100.0
  ge(0:m) = gap_extend * real ( gap_extend_percent(0:m) ) / 100.0

  do test = 1, 4
!
!  Set the sequence name and values.
!
    if ( test == 1 ) then

      sequence_name = 's00657 (chunk)'

      sequence = 'NGQSYRGTYSTTVTGRTCQAWSSMTPHSHS'

    else if ( test == 2 ) then

      sequence_name = 'c61545 (chunk)'

      sequence = 'DADKSPWCYTTDPRVRWEFCNLKK'

    else if ( test == 3 ) then

      sequence_name = 'plhu-k (chunk)'

      sequence = 'TPHRHQKTPENYPNAGLTMNYCRNPDADKG'

    else if ( test == 4 ) then

      sequence_name = 'Made-Up (chunk)'

      sequence = 'WSSMTPHRHQCATKTPENYPNAGLTM'

    end if
!
!  Convert the sequence string to a vector of characters.
!
    n = len_trim ( sequence )

    call s_to_chvec ( sequence, n, b )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Sequence name: ' // trim ( sequence_name )
    write ( *, '(a)' ) ' '
!
!  Determine the corners of the optimal local alignment.
!
    m1 = 0
    m2 = m
    n1 = 0
    n2 = n

    call ps_lg_corners ( b, m, m1, m2, n, n1, n2, profile_score, go, &
      ge, i1, j1, i2, j2 )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PS_LG_CORNERS:'
    write ( *, '(a)' ) '  Optimal alignment "corners" are:'
    write ( *, '(a,2i6)' ) '    I1, J1 = ', i1, j1
    write ( *, '(a,2i6)' ) '    I2, J2 = ', i2, j2
!
!  Determine the path of the optimal local alignment.
!
    call ps_lg_rpl ( b, m, i1, i2, n, j1, j2, profile_score, go, ge, &
      npath, pathi, pathj )

    call ps_lg_match_print ( b, m, n, npath, pathi, pathj, profile_score, &
      go, ge )

  end do

  return
end
function profile_score ( p1, c2 )

!*****************************************************************************80
!
!! PROFILE_SCORE computes a single entry profile/sequence matching score.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) P1, the position in the profile.
!
!    Input, character C2, the character in the sequence that is
!    to be matched with profile position P1.
!
!    Output, real PROFILE_SCORE, the score for matching the given profile
!    position with the sequence character.
!
  implicit none

  integer ( kind = 4 ), parameter :: acid_num = 23
  integer ( kind = 4 ), parameter :: m_max = 200

  integer ( kind = 4 ) a_to_i4
  character, save, dimension ( acid_num ) :: acid_code
  integer ( kind = 4 ), save, dimension ( 26 ) :: acid_index
  character c2
  character, save, dimension ( m_max ) :: conserved
  logical, parameter :: debug = .false.
  integer ( kind = 4 ), save, dimension ( m_max ) :: entropy
  character ( len = 80 ) file_name
  integer ( kind = 4 ), save, dimension ( m_max ) :: gap_extend_percent
  integer ( kind = 4 ), save, dimension ( m_max ) :: gap_open_percent
  integer ( kind = 4 ) i
  character i4_to_a
  integer ( kind = 4 ) ic2
  integer ( kind = 4 ) index_c2
  logical, save :: initialized = .false.
  integer ( kind = 4 ) ios
  integer ( kind = 4 ), parameter :: iunit = 1
  integer ( kind = 4 ), save :: m
  integer ( kind = 4 ) p1
  real profile_score
  integer ( kind = 4 ), save, dimension ( m_max, acid_num ) :: score1
  integer ( kind = 4 ), save, dimension ( acid_num ) :: score2

  if ( .not. initialized ) then

    file_name = 'profile.txt'
!
!  Read the data from the file.
!
    open ( unit = iunit, file = file_name, status = 'old', iostat = ios )

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PROFILE_SCORE - Fatal error!'
      write ( *, '(a)' ) '  Could not open the profile file.'
      return
    end if

    call profile_score_read ( acid_code, acid_num, conserved, entropy, &
      gap_open_percent, gap_extend_percent, iunit, m_max, m, score1, score2 )

    close ( unit = iunit )

    call a_index ( acid_num, acid_code, acid_index )

    if ( debug ) then

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '    I    Acid_code'
      write ( *, '(a)' ) ' '
      do i = 1, acid_num
        write ( *, '(i5,5x,a1)' ) i, acid_code(i)
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '    I    Acid_index'
      write ( *, '(a)' ) ' '
      do i = 1, 26
        write ( *, '(i5,5x,a1,4x,i5)' ) i, i4_to_a(i), acid_index(i)
      end do

    end if

    initialized = .true.

  end if

  ic2 = a_to_i4 ( c2 )

  index_c2 = acid_index(ic2)

  if ( index_c2 <= 0 ) then
    profile_score = 0.0
  else if ( p1 <= 0 ) then
    profile_score = 0.0
  else
    profile_score = score1(p1,index_c2)
  end if

  return
end
function ps_score_12 ( p1, c2 )

!*****************************************************************************80
!
!! PS_SCORE_12 computes a profile/sequence matching score.
!
!  Discussion:
!
!    Either I have a subtle programming error, or the FORTRAN 90 compiler
!    is misbehaving.  When I try to do the comparisons using character
!    values, when C2 = 'T', the comparison fails, so I converted to
!    integer ( kind = 4 ) comparisons.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) P1, the position in the profile.
!
!    Input, character C2, the character in the sequence that is
!    to be matched with profile position P1.
!
!    Output, real PS_SCORE_12, the score for matching the given profile
!    position with the sequence character.
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) c
  character c2
  integer ( kind = 4 ) g
  integer ( kind = 4 ) i
  integer ( kind = 4 ) p1
  real ps_score_12
  real score
  integer ( kind = 4 ) t

  score = 0.0

  i = ichar ( c2 )

  a = ichar ( 'A' )
  c = ichar ( 'C' )
  g = ichar ( 'G' )
  t = ichar ( 'T' )

  if ( p1 == 1 ) then

         if ( i == a ) then
      score = -1.0
    else if ( i == c ) then
      score = 3.0
    else if ( i == g ) then
      score = 5.0
    else if ( i == t ) then
      score = -5.0
    else
      score = -10.0
    end if

  else if ( p1 == 2 ) then

         if ( i == a ) then
      score = -5.0
    else if ( i == c ) then
      score = 5.0
    else if ( i == g ) then
      score = 3.0
    else if ( i == t ) then
      score = -1.0
    else
      score = -10.0
    end if

  else if ( p1 == 3 ) then

         if ( i == a ) then
      score = 3.0
    else if ( i == c ) then
      score = -1.0
    else if ( i == g ) then
      score = -5.0
    else if ( i == t ) then
      score = 5.0
    else
      score = -10.0
    end if

  else if ( p1 == 4 ) then

         if ( i == a ) then
      score = -1.0
    else if ( i == c ) then
      score = 3.0
    else if ( i == g ) then
      score = 5.0
    else if ( i == t ) then
      score = -5.0
    else
      score = -10.0
    end if

  else if ( p1 == 5 ) then

         if ( i == a ) then
      score = 5.0
    else if ( i == c ) then
      score = -5.0
    else if ( i == g ) then
      score = -1.0
    else if ( i == t ) then
      score = 3.0
    else
      score = -10.0
    end if

  else if ( p1 == 6 ) then

         if ( i == a ) then
      score = 3.0
    else if ( i == c ) then
      score = -1.0
    else if ( i == g ) then
      score = -5.0
    else if ( i == t ) then
      score = 5.0
    else
      score = -10.0
    end if

  else if ( p1 == 7 ) then

         if ( i == a ) then
      score = 5.0
    else if ( i == c ) then
      score = -5.0
    else if ( i == g ) then
      score = -1.0
    else if ( i == t ) then
      score = 3.0
    else
      score = -10.0
    end if

  else if ( p1 == 8 ) then

         if ( i == a ) then
      score = 3.0
    else if ( i == c ) then
      score = -1.0
    else if ( i == g ) then
      score = -5.0
    else if ( i == t ) then
      score = 5.0
    else
      score = -10.0
    end if

  else if ( p1 == 9 ) then

         if ( i == a ) then
      score = 5.0
    else if ( i == c ) then
      score = -5.0
    else if ( i == g ) then
      score = -1.0
    else if ( i == t ) then
      score = 3.0
    else
      score = -10.0
    end if

  else if ( p1 == 10 ) then

         if ( i == a ) then
      score = -1.0
    else if ( i == c ) then
      score = 3.0
    else if ( i == g ) then
      score = 5.0
    else if ( i == t ) then
      score = -5.0
    else
      score = -10.0
    end if

  else if ( p1 == 11 ) then

         if ( i == a ) then
      score = -5.0
    else if ( i == c ) then
      score = 5.0
    else if ( i == g ) then
      score = 3.0
    else if ( i == t ) then
      score = -1.0
    else
      score = -10.0
    end if

  else if ( p1 == 12 ) then

         if ( i == a ) then
      score = 3.0
    else if ( i == c ) then
      score = -1.0
    else if ( i == g ) then
      score = -5.0
    else if ( i == t ) then
      score = 5.0
    else
      score = -10.0
    end if

  else

    score = - 10000.0

  end if

  ps_score_12 = score

  return
end
function ps_score_simple ( p1, c2 )

!*****************************************************************************80
!
!! PS_SCORE_SIMPLE computes a single entry profile/sequence matching score.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) P1, the position in the profile.
!
!    Input, character C2, the character in the sequence that is
!    to be matched with profile position P1.
!
!    Output, real PS_SCORE_SIMPLE, the score for matching the given profile
!    position with the sequence character.
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) c
  character c2
  integer ( kind = 4 ) g
  integer ( kind = 4 ) i
  integer ( kind = 4 ) p1
  real ps_score_simple
  real score
  integer ( kind = 4 ) t

  i = ichar ( c2 )

  a = ichar ( 'A' )
  c = ichar ( 'C' )
  g = ichar ( 'G' )
  t = ichar ( 'T' )

  score = 0.0

       if ( p1 == 0 ) then

         if ( i == a ) then
      score = 5.0
    else if ( i == c ) then
      score = -1.0
    else if ( i == g ) then
      score = 1.0
    else if ( i == t ) then
      score = -5.0
    else
      score = -10.0
    end if

  else if ( p1 == 1 ) then

         if ( i == a ) then
      score = 5.0
    else if ( i == c ) then
      score = -1.0
    else if ( i == g ) then
      score = 1.0
    else if ( i == t ) then
      score = -5.0
    else
      score = -10.0
    end if

  else if ( p1 == 2 ) then

         if ( i == a ) then
      score = -1.0
    else if ( i == c ) then
      score = 4.0
    else if ( i == g ) then
      score = -1.0
    else if ( i == t ) then
      score = -5.0
    else
      score = -10.0
    end if

  else if ( p1 == 3 ) then

         if ( i == a ) then
      score = 4.0
    else if ( i == c ) then
      score = -5.0
    else if ( i == g ) then
      score = 4.0
    else if ( i == t ) then
      score = -5.0
    else
      score = -10.0
    end if

  else if ( p1 == 4 ) then

         if ( i == a ) then
      score = -1.0
    else if ( i == c ) then
      score = -1.0
    else if ( i == g ) then
      score = -1.0
    else if ( i == t ) then
      score = 3.0
    else
      score = -10.0
    end if

  else if ( p1 == 5 ) then

         if ( i == a ) then
      score = 3.0
    else if ( i == c ) then
      score = -2.0
    else if ( i == g ) then
      score = 2.0
    else if ( i == t ) then
      score = -2.0
    else
      score = -10.0
    end if

  end if

  ps_score_simple = score

  return
end
