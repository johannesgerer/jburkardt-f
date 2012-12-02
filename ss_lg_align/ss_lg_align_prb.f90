program main

!*****************************************************************************80
!
!! MAIN is the main program for SS_LG_ALIGN_PRB.
!
!!  Discsussion:
!
!    SS_LG_ALIGN_PRB carries out some alignment tests.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SS_LG_ALIGN_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the SS_LG_ALIGN library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
  call test05 ( )
  call test055 ( )
  call test06 ( )
  call test07 ( )
  call test08 ( )
  call test09 ( )
  call test12 ( )
  call test13 ( )
  call test14 ( )
  call test15 ( )
  call test16 ( )
  call test17 ( )
  call test18 ( )
  call test19 ( )
  call test20 ( )
  call test21 ( )
  call test22 ( )
  call test23 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SS_LG_ALIGN_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests SS_LG_FSL, SS_LG_FSQ.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: m = 12
  integer, parameter :: n = 12
  integer, parameter :: lds = m

  character a(m)
  character ( len = m ) apack
  character b(n)
  real base
  character ( len = n ) bpack
  real e(0:lds,0:n)
  real ev(0:n)
  real f(0:lds,0:n)
  real fv(0:n)
  real, parameter :: gap_extend = -0.5
  real, parameter :: gap_open = -2.0
  integer i
  integer j
  integer m1
  integer m2
  integer n1
  integer n2
  real s(0:lds,0:n)
  real sv(0:n)
  integer t(0:lds,0:n)
  integer tv(0:n)
  real, external :: pam120_score

  apack = 'GCTGATATAGCT'
  bpack = 'CCGTGATTAGCT'

  do i = 1, m
    a(i) = apack(i:i)
  end do

  do i = 1, n
    b(i) = bpack(i:i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01:'
  write ( *, '(a)' ) '  SS_LG_FSQ - Forward score quadratic;'
  write ( *, '(a)' ) '  SS_LG_FSL - Forward score linear;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  GAP_OPEN penalty =   ', gap_open
  write ( *, '(a,g14.6)' ) '  GAP_EXTEND penalty = ', gap_extend
  write ( *, '(a)' ) '  Matching scores by PAM120_SCORE.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Verify that the FSQ and FSL tables agree.'

  call chvec2_print ( m, a, n, b, '  Sequences A and B:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Matching Scores:'
  write ( *, '(a)' ) ' '
  write ( *, '(3x,13i5)' ) ( j, j = 0, n )
  write ( *, '(a)' ) ' '
  do i = 0, m
    if ( i == 0 ) then
      write ( *, '(i3,13f5.1)' ) i, ( 0.0, j = 0, n )
    else
      write ( *, '(i3,13f5.1)' ) i, 0.0, &
        ( pam120_score ( a(i), b(j) ), j = 1, n )
    end if
  end do

  m1 = 0
  m2 = m
  n1 = 0
  n2 = n
  base = gap_open

  call ss_lg_fsq ( a, b, m, m1, m2, n, n1, n2, pam120_score, gap_open, &
    gap_extend, base, lds, s, e, f, t )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SS_LG_FSQ:'
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

    call ss_lg_fsl ( a, b, m, m1, m2, n, n1, n2, pam120_score, gap_open, &
      gap_extend, base, sv, ev, fv, tv )

    s(m2,0:n) = sv(0:n)
    e(m2,0:n) = ev(0:n)
    f(m2,0:n) = fv(0:n)
    t(m2,0:n) = tv(0:n)

  end do

  m2 = m

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SS_LG_FSL:'
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
  do i = m2, m2
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
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests SS_LG_BSQ, SS_LG_BSL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: m = 12
  integer, parameter :: n = 12
  integer, parameter :: lds = m

  character a(m)
  character ( len = m ) apack
  character b(n)
  real base
  character ( len = n ) bpack
  real e(0:lds,0:n)
  real ev(0:n)
  real f(0:lds,0:n)
  real fv(0:n)
  real, parameter :: gap_extend = -0.5
  real, parameter :: gap_open = -2.0
  integer i
  integer j
  integer m1
  integer m2
  integer n1
  integer n2
  real s(0:lds,0:n)
  real sv(0:n)
  integer t(0:lds,0:n)
  integer tv(0:n)
  real, external :: pam120_score

  apack = 'GCTGATATAGCT'
  bpack = 'CCGTGATTAGCT'

  do i = 1, m
    a(i) = apack(i:i)
  end do

  do i = 1, n
    b(i) = bpack(i:i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02:'
  write ( *, '(a)' ) '  SS_LG_BSQ - Backward score quadratic;'
  write ( *, '(a)' ) '  SS_LG_BSL - Backward score linear.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  GAP_OPEN penalty =   ', gap_open
  write ( *, '(a,g14.6)' ) '  GAP_EXTEND penalty = ', gap_extend
  write ( *, '(a)' ) '  Matching scores by PAM120_SCORE.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Verify that the BSQ and BSL tables agree.'
  write ( *, '(a)' ) ' '

  call chvec2_print ( m, a, n, b, '  Sequences A and B:' )

  m1 = 0
  m2 = m
  n1 = 0
  n2 = n
  base = gap_open

  call ss_lg_bsq ( a, b, m, m1, m2, n, n1, n2, pam120_score, gap_open, &
    gap_extend, base, lds, s, e, f, t )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SS_LG_BSQ:'
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

  m2 = m
  n1 = 0
  n2 = n
  base = gap_open

  do m1 = 0, m2

    call ss_lg_bsl ( a, b, m, m1, m2, n, n1, n2, pam120_score, gap_open, &
      gap_extend, base, sv, ev, fv, tv )

    s(m1,0:n) = sv(0:n)
    e(m1,0:n) = ev(0:n)
    f(m1,0:n) = fv(0:n)
    t(m1,0:n) = tv(0:n)

  end do

  m1 = 0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SS_LG_BSL:'
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
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests SS_LG_FSQ, SS_LG_FSL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: m = 3
  integer, parameter :: n = 4
  integer, parameter :: lds = m

  character a(m)
  character ( len = m ) apack
  character b(n)
  real base
  character ( len = n ) bpack
  real e(0:lds,0:n)
  real ev(0:n)
  real f(0:lds,0:n)
  real fv(0:n)
  real, parameter :: gap_extend = -0.5
  real, parameter :: gap_open = -2.0
  integer i
  integer imax
  integer j
  integer m1
  integer m2
  integer n1
  integer n2
  real s(0:lds,0:n)
  real s_max
  real sv(0:n)
  real, external :: pam120_score
  integer t(0:lds,0:n)
  integer tv(0:n)

  apack = 'GCT'
  bpack = 'GGGT'

  do i = 1, m
    a(i) = apack(i:i)
  end do

  do i = 1, n
    b(i) = bpack(i:i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03:'
  write ( *, '(a)' ) '  SS_LG_FSQ - Forward score quadratic;'
  write ( *, '(a)' ) '  SS_LG_FSL - Forward score linear;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  GAP_OPEN penalty =   ', gap_open
  write ( *, '(a,g14.6)' ) '  GAP_EXTEND penalty = ', gap_extend
  write ( *, '(a)' ) '  Matching scores by PAM120_SCORE.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Verify that the FSQ and FSL tables agree.'

  call chvec2_print ( m, a, n, b, '  Sequences A and B:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Matching Scores:'
  write ( *, '(a)' ) ' '
  write ( *, '(3x,13i5)' ) ( j, j = 0, n )
  write ( *, '(a)' ) ' '
  do i = 0, m
    if ( i == 0 ) then
      write ( *, '(i3,13f5.1)' ) i, ( 0.0, j = 0, n )
    else
      write ( *, '(i3,13f5.1)' ) i, 0.0, ( pam120_score ( a(i), b(j) ), j = 1, n )
    end if
  end do

  m1 = 0
  m2 = m
  n1 = 0
  n2 = n
  base = gap_open

  call ss_lg_fsq ( a, b, m, m1, m2, n, n1, n2, pam120_score, gap_open, &
    gap_extend, base, lds, s, e, f, t )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SS_LG_FSQ:'
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

    call ss_lg_fsl ( a, b, m, m1, m2, n, n1, n2, pam120_score, gap_open, &
      gap_extend, base, sv, ev, fv, tv )

    s(m2,0:n) = sv(0:n)
    e(m2,0:n) = ev(0:n)
    f(m2,0:n) = fv(0:n)
    t(m2,0:n) = tv(0:n)

  end do

  m2 = m

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SS_LG_FSL:'
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
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests SS_LG_BSQ, SS_LG_BSL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: m = 3
  integer, parameter :: n = 4
  integer, parameter :: lds = m

  character a(m)
  character ( len = m ) apack
  character b(n)
  real base
  character ( len = n ) bpack
  real e(0:lds,0:n)
  real ev(0:n)
  real f(0:lds,0:n)
  real fv(0:n)
  real, parameter :: gap_extend = -0.5
  real, parameter :: gap_open = -2.0
  integer i
  integer imax
  integer j
  integer m1
  integer m2
  integer n1
  integer n2
  real s(0:lds,0:n)
  real s_max
  real sv(0:n)
  real, external :: pam120_score
  integer t(0:lds,0:n)
  integer tv(0:n)

  apack = 'GCT'
  bpack = 'GGGT'

  do i = 1, m
    a(i) = apack(i:i)
  end do

  do i = 1, n
    b(i) = bpack(i:i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04:'
  write ( *, '(a)' ) '  SS_LG_BSQ - Backward score quadratic;'
  write ( *, '(a)' ) '  SS_LG_BSL - Backward score linear.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  GAP_OPEN penalty =   ', gap_open
  write ( *, '(a,g14.6)' ) '  GAP_EXTEND penalty = ', gap_extend
  write ( *, '(a)' ) '  Matching scores by PAM120_SCORE.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Verify that the BSQ and BSL tables agree.'

  call chvec2_print ( m, a, n, b, '  Sequences A and B:' )

  m1 = 0
  m2 = m
  n1 = 0
  n2 = n
  base = gap_open

  call ss_lg_bsq ( a, b, m, m1, m2, n, n1, n2, pam120_score, gap_open, &
    gap_extend, base, lds, s, e, f, t )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SS_LG_BSQ:'
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

  m2 = m
  n1 = 0
  n2 = n
  base = gap_open

  do m1 = 0, m2

    call ss_lg_bsl ( a, b, m, m1, m2, n, n1, n2, pam120_score, gap_open, &
      gap_extend, base, sv, ev, fv, tv )

    s(m1,0:n) = sv(0:n)
    e(m1,0:n) = ev(0:n)
    f(m1,0:n) = fv(0:n)
    t(m1,0:n) = tv(0:n)

  end do

  m1 = 0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SS_LG_BSL:'
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
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 tests SS_LG_BSQ, SS_LG_FSQ.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: m = 3
  integer, parameter :: n = 4
  integer, parameter :: lds = m

  character a(m)
  character ( len = m ) apack
  character b(n)
  real base
  character ( len = n ) bpack
  real e(0:lds,0:n)
  real f(0:lds,0:n)
  real, parameter :: gap_extend = -0.5
  real, parameter :: gap_open = -2.0
  integer i
  integer imax
  integer j
  integer m1
  integer m2
  integer n1
  integer n2
  real s_max
  real s(0:lds,0:n)
  real, external :: pam120_score
  integer t(0:lds,0:n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05:'
  write ( *, '(a)' ) '  SS_LG_BSQ - Backward score quadratic;'
  write ( *, '(a)' ) '  SS_LG_FSQ - Forward score quadratic.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  GAP_OPEN penalty =   ', gap_open
  write ( *, '(a,g14.6)' ) '  GAP_EXTEND penalty = ', gap_extend
  write ( *, '(a)' ) '  Matching scores by PAM120_SCORE.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Verify that the FSQ and (reversed) BSQ tables agree.'

  apack = 'GCT'
  bpack = 'GGGT'

  do i = 1, m
    a(i) = apack(i:i)
  end do

  do i = 1, n
    b(i) = bpack(i:i)
  end do

  call chvec2_print ( m, a, n, b, '  Sequences A and B:' )

  m1 = 0
  m2 = m
  n1 = 0
  n2 = n
  base = gap_open

  call ss_lg_fsq ( a, b, m, m1, m2, n, n1, n2, pam120_score, gap_open, &
    gap_extend, base, lds, s, e, f, t )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SS_LG_FSQ:'
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
!
!  Now apply BSQ to reversed problem.
!
  apack = 'TCG'
  bpack = 'TGGG'

  do i = 1, m
    a(i) = apack(i:i)
  end do

  do i = 1, n
    b(i) = bpack(i:i)
  end do

  call chvec2_print ( m, a, n, b, '  Reversed sequences A and B:' )

  m1 = 0
  m2 = m
  n1 = 0
  n2 = n
  base = gap_open

  call ss_lg_bsq ( a, b, m, m1, m2, n, n1, n2, pam120_score, gap_open, &
    gap_extend, base, lds, s, e, f, t )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SS_LG_BSQ:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  (The tables are printed in reverse order)'
  write ( *, '(a)' ) ' '
  write ( *, '(3x,13i5)' ) ( j, j = n2, n1, -1 )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SB:'
  write ( *, '(a)' ) ' '
  do i = m2, m1, -1
    write ( *, '(i3,13f5.1)' ) i, s(i,n2:n1:-1)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'EB:'
  write ( *, '(a)' ) ' '
  do i = m2, m1, -1
    write ( *, '(i3,13f5.1)' ) i, e(i,n2:n1:-1)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FB:'
  write ( *, '(a)' ) ' '
  do i = m2, m1, -1
    write ( *, '(i3,13f5.1)' ) i, f(i,n2:n1:-1)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TB:'
  write ( *, '(a)' ) ' '
  do i = m2, m1, -1
    write ( *, '(i3,13i5)' ) i, t(i,n2:n1:-1)
  end do

  return
end
subroutine test055 ( )

!*****************************************************************************80
!
!! TEST055 tests SS_LG_FSQ, SS_LG_BSQ.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: m = 12
  integer, parameter :: n = 12
  integer, parameter :: lds = m

  character a(m)
  character ( len = m ) apack
  character b(n)
  real base
  character ( len = n ) bpack
  real e(0:lds,0:n)
  real f(0:lds,0:n)
  real, parameter :: gap_extend = -0.5
  real, parameter :: gap_open = -2.0
  integer i
  integer i1
  integer i2
  integer j1
  integer j2
  integer imax
  integer j
  integer m1
  integer m2
  integer n1
  integer n2
  integer npath
  integer pathi(m+n+1)
  integer pathj(m+n+1)
  real s_max
  real s(0:lds,0:n)
  real, external :: pam120_score
  integer t(0:lds,0:n)

  apack = 'GCTGATATAGCT'
  bpack = 'CCGTGATTAGCT'

  do i = 1, m
    a(i) = apack(i:i)
  end do

  do i = 1, n
    b(i) = bpack(i:i)
  end do

  s(0:m,0:n) = 0.0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST055:'
  write ( *, '(a)' ) '  SS_LG_FSQ - Forward quadratic;'
  write ( *, '(a)' ) '  SS_LG_BSQ - Backward quadratic.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  GAP_OPEN penalty =   ', gap_open
  write ( *, '(a,g14.6)' ) '  GAP_EXTEND penalty = ', gap_extend
  write ( *, '(a)' ) '  Matching scores by PAM120_SCORE.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Do a score computation on a sub-block.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Check that the highest entry in SF for FSQ agrees with'
  write ( *, '(a)' ) '  the highest entry for SB BSQ.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The paths prescribed from the highest entries in SF and'
  write ( *, '(a)' ) '  SB should agree (although opposite in direction).'

  call chvec2_print ( m, a, n, b, '  Sequences A and B:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Matching Scores:'
  write ( *, '(a)' ) ' '
  write ( *, '(3x,13i5)' ) ( j, j = 0, n )
  write ( *, '(a)' ) ' '
  do i = 0, m
    if ( i == 0 ) then
      write ( *, '(i3,13f5.1)' ) i, ( 0.0, j = 0, n )
    else
      write ( *, '(i3,13f5.1)' ) i, 0.0, ( pam120_score ( a(i), b(j) ), j = 1, n )
    end if
  end do

  m1 = 4
  m2 = 10
  n1 = 3
  n2 = 9
  base = gap_open

  call ss_lg_fsq ( a, b, m, m1, m2, n, n1, n2, pam120_score, gap_open, &
    gap_extend, base, lds, s, e, f, t )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SS_LG_FSQ:'
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
  write ( *, '(a)' ) '  SS_LG_FSQ:'
  write ( *, '(a,2i6)' ) '    Optimal score indices: ', i2, j2
  write ( *, '(a,g14.6)' ) '    Optimal score:         ', s(i2,j2)

  call ss_lg_fpq ( i2, j2, m, m1, m2, n, n1, n2, lds, t, npath, pathi, pathj )

  call i4vec2_print ( npath, pathi, pathj, '  Matching path:' )

  call ss_lg_match_print ( a, b, m, n, npath, pathi, pathj, &
    pam120_score, gap_open, gap_extend )

  base = gap_open

  call ss_lg_bsq ( a, b, m, m1, m2, n, n1, n2, pam120_score, gap_open, &
    gap_extend, base, lds, s, e, f, t )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SS_LG_BSQ:'
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
  write ( *, '(a)' ) '  SS_LG_BSQ:'
  write ( *, '(a,2i6)' ) '    Optimal score indices: ', i1, j1
  write ( *, '(a,g14.6)' ) '    Optimal score:         ', s(i1,j1)

  call ss_lg_bpq ( i1, j1, m, m1, m2, n, n1, n2, lds, t, npath, pathi, pathj )

  call i4vec2_print ( npath, pathi, pathj, '  Matching path:' )

  call ss_lg_match_print ( a, b, m, n, npath, pathi, pathj, &
    pam120_score, gap_open, gap_extend )

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 tests SS_LG_FSL, SS_LG_BSL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: m = 12
  integer, parameter :: n = 12
  integer, parameter :: lds = m

  character a(m)
  character ( len = m ) apack
  character b(n)
  real base
  character ( len = n ) bpack
  real ev(0:n)
  real fv(0:n)
  real, parameter :: gap_extend = -0.5
  real, parameter :: gap_open = -2.0
  integer i
  integer imax
  integer j
  integer m1
  integer m2
  integer n1
  integer n2
  real s_max
  real sv(0:n)
  real, external :: pam120_score
  integer tv(0:n)

  apack = 'GCTGATATAGCT'
  bpack = 'CCGTGATTAGCT'

  do i = 1, m
    a(i) = apack(i:i)
  end do

  do i = 1, n
    b(i) = bpack(i:i)
  end do

  sv(0:n) = 0.0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06:'
  write ( *, '(a)' ) '  SS_LG_FSL - Forward linear;'
  write ( *, '(a)' ) '  SS_LG_BSL - Backward linear.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  GAP_OPEN penalty =   ', gap_open
  write ( *, '(a,g14.6)' ) '  GAP_EXTEND penalty = ', gap_extend
  write ( *, '(a)' ) '  Matching scores by PAM120_SCORE.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Do a score computation on a sub-block.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Verify that the maximum entries for FSL and BSL'
  write ( *, '(a)' ) '  are the same.'

  call chvec2_print ( m, a, n, b, '  Sequences A and B:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Matching Scores:'
  write ( *, '(a)' ) ' '
  write ( *, '(3x,13i5)' ) ( j, j = 0, n )
  write ( *, '(a)' ) ' '
  do i = 0, m
    if ( i == 0 ) then
      write ( *, '(i3,13f5.1)' ) i, ( 0.0, j = 0, n )
    else
      write ( *, '(i3,13f5.1)' ) i, 0.0, ( pam120_score ( a(i), b(j) ), j = 1, n )
    end if
  end do

  m1 = 4
  n1 = 3
  n2 = 9
  base = gap_open

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SS_LG_FSL:'
  write ( *, '(a)' ) ' '
  write ( *, '(3x,13i5)' ) ( j, j = n1, n2 )
  write ( *, '(a)' ) ' '

  do m2 = m1, 10

    call ss_lg_fsl ( a, b, m, m1, m2, n, n1, n2, pam120_score, gap_open, &
      gap_extend, base, sv, ev, fv , tv)

    write ( *, '(i3,13f5.1)' ) m2, sv(n1:n2)

  end do

  m2 = 10
  n1 = 3
  n2 = 9
  base = gap_open

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SS_LG_BSL:'
  write ( *, '(a)' ) ' '
  write ( *, '(3x,13i5)' ) ( j, j = n1, n2 )
  write ( *, '(a)' ) ' '

  do m1 = 4, m2

    call ss_lg_bsl ( a, b, m, m1, m2, n, n1, n2, pam120_score, gap_open, &
      gap_extend, base, sv, ev, fv, tv )

    write ( *, '(i3,13f5.1)' ) m1, sv(n1:n2)

  end do

  return
end
subroutine test07 ( )

!*****************************************************************************80
!
!! TEST07 tests SS_LG_FSQ, SS_LG_FPQ.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: m = 12
  integer, parameter :: n = 12
  integer, parameter :: lds = m

  character a(m)
  character ( len = m ) apack
  character b(n)
  real base
  character ( len = n ) bpack
  real e(0:lds,0:n)
  real f(0:lds,0:n)
  real, parameter :: gap_extend = -0.5
  real, parameter :: gap_open = -2.0
  integer i
  integer i2
  integer j2
  integer m1
  integer m2
  integer n1
  integer n2
  integer npath
  integer pathi(m+n+1)
  integer pathj(m+n+1)
  real s(0:lds,0:n)
  integer t(0:lds,0:n)
  real, external :: pam120_score

  apack = 'GCTGATATAGCT'
  bpack = 'CCGTGATTAGCT'

  do i = 1, m
    a(i) = apack(i:i)
  end do

  do i = 1, n
    b(i) = bpack(i:i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07:'
  write ( *, '(a)' ) '  SS_LG_FSQ - forward score quadratic'
  write ( *, '(a)' ) '  SS_LG_FPQ - forward path quadratic'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  GAP_OPEN penalty =   ', gap_open
  write ( *, '(a,g14.6)' ) '  GAP_EXTEND penalty = ', gap_extend
  write ( *, '(a)' ) '  Matching scores by PAM120_SCORE.'

  call chvec2_print ( m, a, n, b, '  Sequences A and B:' )

  m1 = 0
  m2 = m
  n1 = 0
  n2 = n
  base = gap_open

  call ss_lg_fsq ( a, b, m, m1, m2, n, n1, n2, pam120_score, gap_open, &
    gap_extend, base, lds, s, e, f, t )

  call r4mat_imax ( lds+1, m+1, n+1, s, i2, j2 )
  i2 = i2 - 1
  j2 = j2 - 1
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) &
    '  SS_LG_FSQ reports optimal matching score is ', s(i2,j2)

  call ss_lg_fpq ( i2, j2, m, m1, m2, n, n1, n2, lds, t, npath, pathi, pathj )

  call i4vec2_print ( npath, pathi, pathj, '  Matching path:' )

  call ss_lg_match_print ( a, b, m, n, npath, pathi, pathj, &
    pam120_score, gap_open, gap_extend )

  return
end
subroutine test08 ( )

!*****************************************************************************80
!
!! TEST08 tests SS_LG_BSQ, SS_LG_BPQ.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: m = 12
  integer, parameter :: n = 12
  integer, parameter :: lds = m

  character a(m)
  character ( len = m ) apack
  character b(n)
  real base
  character ( len = n ) bpack
  real e(0:lds,0:n)
  real f(0:lds,0:n)
  real, parameter :: gap_extend = -0.5
  real, parameter :: gap_open = -2.0
  integer i
  integer i1
  integer j1
  integer m1
  integer m2
  integer n1
  integer n2
  integer npath
  real, external :: pam120_score
  integer pathi(m+n+1)
  integer pathj(m+n+1)
  real s(0:lds,0:n)
  integer t(0:lds,0:n)

  apack = 'GCTGATATAGCT'
  bpack = 'CCGTGATTAGCT'

  do i = 1, m
    a(i) = apack(i:i)
  end do

  do i = 1, n
    b(i) = bpack(i:i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08:'
  write ( *, '(a)' ) '  SS_LG_BSQ - backward score quadratic'
  write ( *, '(a)' ) '  SS_LG_BPQ - backward path quadratic'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  GAP_OPEN penalty =   ', gap_open
  write ( *, '(a,g14.6)' ) '  GAP_EXTEND penalty = ', gap_extend
  write ( *, '(a)' ) '  Matching scores by PAM120_SCORE.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Compare path and score with FSQ/FPQ computation.'

  call chvec2_print ( m, a, n, b, '  Sequences A and B:' )

  m1 = 0
  m2 = m
  n1 = 0
  n2 = n
  base = gap_open

  call ss_lg_bsq ( a, b, m, m1, m2, n, n1, n2, pam120_score, gap_open, &
    gap_extend, base, lds, s, e, f, t )

  call r4mat_imax ( lds+1, m+1, n+1, s, i1, j1 )
  i1 = i1 - 1
  j1 = j1 - 1

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) &
    '  SS_LG_BSQ reports optimal matching score is ', s(i1,j1)

  call ss_lg_bpq ( i1, j1, m, m1, m2, n, n1, n2, lds, t, npath, pathi, pathj )

  call i4vec2_print ( npath, pathi, pathj, '  Matching path:' )

  call ss_lg_match_print ( a, b, m, n, npath, pathi, pathj, &
    pam120_score, gap_open, gap_extend )

  return
end
subroutine test09 ( )

!*****************************************************************************80
!
!! TEST09 tests SS_LG_RPL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: m = 12
  integer, parameter :: n = 12

  character a(m)
  character ( len = m ) apack
  character b(n)
  character ( len = n ) bpack
  real, parameter :: gap_extend = -0.5
  real, parameter :: gap_open = -2.0
  integer i
  integer i1
  integer i2
  integer j1
  integer j2
  integer m1
  integer m2
  integer n1
  integer n2
  integer npath
  integer pathi(m+n+1)
  integer pathj(m+n+1)
  real, external :: pam120_score

  apack = 'GCTGATATAGCT'
  bpack = 'CCGTGATTAGCT'

  do i = 1, m
    a(i) = apack(i:i)
  end do

  do i = 1, n
    b(i) = bpack(i:i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST09:'
  write ( *, '(a)' ) '  SS_LG_RPL - path routine;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  GAP_OPEN penalty =   ', gap_open
  write ( *, '(a,g14.6)' ) '  GAP_EXTEND penalty = ', gap_extend
  write ( *, '(a)' ) '  Matching scores by PAM120_SCORE.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Verify that RPL agrees with FSQ/FPQ and BSQ/BPQ.'

  call chvec2_print ( m, a, n, b, '  Sequences A and B:' )

  m1 = 0
  m2 = m
  n1 = 0
  n2 = n

  call ss_lg_corners ( a, b, m, m1, m2, n, n1, n2, pam120_score, gap_open, &
    gap_extend, i1, j1, i2, j2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SS_LG_CORNERS:'
  write ( *, '(a,2i6)' ) '    (I1,J1) = ', i1, j1
  write ( *, '(a,2i6)' ) '    (I2,J2) = ', i2, j2

  call ss_lg_rpl ( a, b, m, i1, i2, n, j1, j2, pam120_score, gap_open, &
    gap_extend, npath, pathi, pathj )

  call i4vec2_print ( npath, pathi, pathj, '  Matching path:' )

  call ss_lg_match_print ( a, b, m, n, npath, pathi, pathj, &
    pam120_score, gap_open, gap_extend )

  return
end
subroutine test12 ( )

!*****************************************************************************80
!
!! TEST12 tests SS_LG_FSQ, SS_LG_FPQ.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: m = 60
  integer, parameter :: n = 59
  integer, parameter :: lds = m

  character a(m)
  character ( len = m ) apack
  character b(n)
  real base
  character ( len = n ) bpack
  real e(0:lds,0:n)
  real f(0:lds,0:n)
  real, parameter :: gap_extend = -0.5
  real, parameter :: gap_open = -2.0
  integer i
  integer i2
  integer j2
  integer m1
  integer m2
  integer n1
  integer n2
  integer npath
  real, external :: pam120_score
  integer pathi(m+n+1)
  integer pathj(m+n+1)
  real s(0:lds,0:n)
  integer t(0:lds,0:n)

  apack = 'MMAAEAGGEEGGPVTAGAAGGGAAAASGAYPAVCRVKIPAALPVAAAAPFPGLAEAGVAA'
  bpack = 'MMAAEAGGPVTAGAAGGGAACCCAASGAYPAVCRVKIPAALPVAAAAPFPGLAEAGVAA'

  do i = 1, m
    a(i) = apack(i:i)
  end do

  do i = 1, n
    b(i) = bpack(i:i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST12:'
  write ( *, '(a)' ) '  SS_LG_FSQ - Forward score quadratic;'
  write ( *, '(a)' ) '  SS_LG_FPQ - Forward path quadratic.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  GAP_OPEN penalty =   ', gap_open
  write ( *, '(a,g14.6)' ) '  GAP_EXTEND penalty = ', gap_extend
  write ( *, '(a)' ) '  Matching scores by PAM120_SCORE.'

  m1 = 0
  m2 = m
  n1 = 0
  n2 = n
  base = gap_open

  call ss_lg_fsq ( a, b, m, m1, m2, n, n1, n2, pam120_score, gap_open, &
    gap_extend, base, lds, s, e, f, t )

  call r4mat_imax ( lds+1, m+1, n+1, s, i2, j2 )
  i2 = i2 - 1
  j2 = j2 - 1
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) &
    '  SS_LG_FSQ reports optimal matching score is ', s(i2,j2)

  call ss_lg_fpq ( i2, j2, m, m1, m2, n, n1, n2, lds, t, npath, pathi, pathj )

  call ss_lg_match_print ( a, b, m, n, npath, pathi, pathj, &
    pam120_score, gap_open, gap_extend )

  return
end
subroutine test13 ( )

!*****************************************************************************80
!
!! TEST13 tests SS_LG_RPL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: m = 60
  integer, parameter :: n = 59

  character a(m)
  character ( len = m ) apack
  character b(n)
  character ( len = n ) bpack
  real, parameter :: gap_extend = -0.5
  real, parameter :: gap_open = -2.0
  integer i
  integer i1
  integer i2
  integer j1
  integer j2
  integer m1
  integer m2
  integer n1
  integer n2
  integer npath
  integer pathi(m+n+1)
  integer pathj(m+n+1)
  real, external :: pam120_score

  apack = 'MMAAEAGGEEGGPVTAGAAGGGAAAASGAYPAVCRVKIPAALPVAAAAPFPGLAEAGVAA'
  bpack = 'MMAAEAGGPVTAGAAGGGAACCCAASGAYPAVCRVKIPAALPVAAAAPFPGLAEAGVAA'

  do i = 1, m
    a(i) = apack(i:i)
  end do

  do i = 1, n
    b(i) = bpack(i:i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST13:'
  write ( *, '(a)' ) '  SS_LG_RPL - path routine;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  GAP_OPEN penalty =   ', gap_open
  write ( *, '(a,g14.6)' ) '  GAP_EXTEND penalty = ', gap_extend
  write ( *, '(a)' ) '  Matching scores by PAM120_SCORE.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Compare with FSQ/FPQ calculation.'

  m1 = 0
  m2 = m
  n1 = 0
  n2 = n

  call ss_lg_corners ( a, b, m, m1, m2, n, n1, n2, pam120_score, gap_open, &
    gap_extend, i1, j1, i2, j2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SS_LG_CORNERS:'
  write ( *, '(a,2i6)' ) '    (I1,J1) = ', i1, j1
  write ( *, '(a,2i6)' ) '    (I2,J2) = ', i2, j2

  call ss_lg_rpl ( a, b, m, i1, i2, n, j1, j2, pam120_score, gap_open, &
    gap_extend, npath, pathi, pathj )

  call ss_lg_match_print ( a, b, m, n, npath, pathi, pathj, &
    pam120_score, gap_open, gap_extend )

  return
end
subroutine test14 ( )

!*****************************************************************************80
!
!! TEST14 tests SS_LG_FSQ, SS_LG_FSL.
!
!  Discussion:
!
!    Attempt to match the table 1 on page 13 (except that the signs should
!    be opposite.)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Eugene Myers and Webb Miller,
!    Optimal Alignments in Linear Space,
!    CABIOS, volume 4, number 1, 1988, pages 11-17.
!
  implicit none

  integer, parameter :: m = 5
  integer, parameter :: n = 3
  integer, parameter :: lds = m

  character a(m)
  character ( len = m ) apack
  character b(n)
  real base
  character ( len = n ) bpack
  real e(0:lds,0:n)
  real ev(0:n)
  real f(0:lds,0:n)
  real fv(0:n)
  real, parameter :: gap_extend = -0.5
  real, parameter :: gap_open = -2.0
  integer i
  integer m1
  integer m2
  integer j
  integer n1
  integer n2
  real s(0:lds,0:n)
  real sv(0:n)
  real, external :: simple_score
  integer t(0:lds,0:n)
  integer tv(0:n)

  apack = 'AGTAC'
  bpack = 'AAG'

  do i = 1, m
    a(i) = apack(i:i)
  end do

  do i = 1, n
    b(i) = bpack(i:i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST14:'
  write ( *, '(a)' ) '  SS_LG_FSQ - Forward score quadratic;'
  write ( *, '(a)' ) '  SS_LG_FSL - Forward score linear;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  GAP_OPEN penalty =   ', gap_open
  write ( *, '(a,g14.6)' ) '  GAP_EXTEND penalty = ', gap_extend
  write ( *, '(a)' ) '  Matching scores by SIMPLE_SCORE.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Verify that FSQ agrees with Myers/Miller calculation.'
  write ( *, '(a)' ) '  Verify that FSQ and FSL agree.'

  call chvec2_print ( m, a, n, b, '  Sequences A and B:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Matching Scores:'
  write ( *, '(a)' ) ' '
  write ( *, '(3x,13i5)' ) ( j, j = 0, n )
  write ( *, '(a)' ) ' '
  do i = 0, m
    if ( i == 0 ) then
      write ( *, '(i3,13f5.1)' ) i, ( 0.0, j = 0, n )
    else
      write ( *, '(i3,13f5.1)' ) i, 0.0, &
        ( simple_score ( a(i), b(j) ), j = 1, n )
    end if
  end do

  m1 = 0
  m2 = m
  n1 = 0
  n2 = n
  base = gap_open

  call ss_lg_fsq ( a, b, m, m1, m2, n, n1, n2, simple_score, gap_open, &
    gap_extend, base, lds, s, e, f, t )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SS_LG_FSQ:'
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

    call ss_lg_fsl ( a, b, m, m1, m2, n, n1, n2, simple_score, gap_open, &
      gap_extend, base, sv, ev, fv, tv )

    s(m2,0:n) = sv(0:n)
    e(m2,0:n) = ev(0:n)
    f(m2,0:n) = fv(0:n)
    t(m2,0:n) = tv(0:n)

  end do

  m2 = m

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SS_LG_FSL:'
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
subroutine test15 ( )

!*****************************************************************************80
!
!! TEST15 tests SS_LG_BSQ, SS_LG_BSL.
!
!  Discussion:
!
!    Attempt to match the table 1 on page 13 (except that the signs should
!    be opposite.)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Eugene Myers and Webb Miller,
!    Optimal Alignments in Linear Space,
!    CABIOS, volume 4, number 1, 1988, pages 11-17.
!
  implicit none

  integer, parameter :: m = 5
  integer, parameter :: n = 3
  integer, parameter :: lds = m

  character a(m)
  character ( len = m ) apack
  character b(n)
  real base
  character ( len = n ) bpack
  real e(0:lds,0:n)
  real ev(0:n)
  real f(0:lds,0:n)
  real fv(0:n)
  real, parameter :: gap_extend = -0.5
  real, parameter :: gap_open = -2.0
  integer i
  integer m1
  integer m2
  integer j
  integer n1
  integer n2
  real s(0:lds,0:n)
  real sv(0:n)
  real, external :: simple_score
  integer t(0:lds,0:n)
  integer tv(0:n)

  apack = 'AGTAC'
  bpack = 'AAG'

  do i = 1, m
    a(i) = apack(i:i)
  end do

  do i = 1, n
    b(i) = bpack(i:i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST15:'
  write ( *, '(a)' ) '  SS_LG_BSQ - Backward score quadratic;'
  write ( *, '(a)' ) '  SS_LG_BSL - Backward score linear;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  GAP_OPEN penalty =   ', gap_open
  write ( *, '(a,g14.6)' ) '  GAP_EXTEND penalty = ', gap_extend
  write ( *, '(a)' ) '  Matching scores by SIMPLE_SCORE.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Verify that BSQ agrees with Myers/Miller calculation.'
  write ( *, '(a)' ) '  Compare BSQ and BSL results.'

  call chvec2_print ( m, a, n, b, '  Sequences A and B:' )

  m1 = 0
  m2 = m
  n1 = 0
  n2 = n
  base = gap_open

  call ss_lg_bsq ( a, b, m, m1, m2, n, n1, n2, simple_score, gap_open, &
    gap_extend, base, lds, s, e, f, t )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SS_LG_BSQ:'
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

    call ss_lg_bsl ( a, b, m, m1, m2, n, n1, n2, simple_score, gap_open, &
      gap_extend, base, sv, ev, fv, tv )

    s(m2,0:n) = sv(0:n)
    e(m2,0:n) = ev(0:n)
    f(m2,0:n) = fv(0:n)
    t(m2,0:n) = tv(0:n)

  end do

  m1 = 0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SS_LG_BSL:'
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
subroutine test16 ( )

!*****************************************************************************80
!
!! TEST16 tests SS_LG_FSQ, SS_LG_FPQ.
!
!  Discussion:
!
!    Produce the alignment corresponding to the example.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Eugene Myers and Webb Miller,
!    Optimal Alignments in Linear Space,
!    CABIOS, volume 4, number 1, 1988, pages 11-17.
!
  implicit none

  integer, parameter :: m = 5
  integer, parameter :: n = 3
  integer, parameter :: lds = m

  character a(m)
  character ( len = m ) apack
  character b(n)
  real base
  character ( len = n ) bpack
  real e(0:lds,0:n)
  real f(0:lds,0:n)
  real, parameter :: gap_extend = -0.5
  real, parameter :: gap_open = -2.0
  integer i
  integer i2
  integer j
  integer j2
  integer m1
  integer m2
  integer n1
  integer n2
  integer npath
  integer pathi(m+n+1)
  integer pathj(m+n+1)
  real s(0:lds,0:n)
  integer t(0:lds,0:n)
  real, external :: simple_score

  apack = 'AGTAC'
  bpack = 'AAG'

  do i = 1, m
    a(i) = apack(i:i)
  end do

  do i = 1, n
    b(i) = bpack(i:i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST16:'
  write ( *, '(a)' ) '  SS_LG_FSQ - forward score quadratic'
  write ( *, '(a)' ) '  SS_LG_FPQ - forward path quadratic'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  GAP_OPEN penalty =   ', gap_open
  write ( *, '(a,g14.6)' ) '  GAP_EXTEND penalty = ', gap_extend
  write ( *, '(a)' ) '  Matching scores by SIMPLE_SCORE.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Compare results with Myers/Miller calculation'

  call chvec2_print ( m, a, n, b, '  Sequences A and B:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Matching Scores:'
  write ( *, '(a)' ) ' '
  write ( *, '(3x,13i5)' ) ( j, j = 0, n )
  write ( *, '(a)' ) ' '
  do i = 0, m
    if ( i == 0 ) then
      write ( *, '(i3,13f5.1)' ) i, ( 0.0, j = 0, n )
    else
      write ( *, '(i3,13f5.1)' ) i, 0.0, &
        ( simple_score ( a(i), b(j) ), j = 1, n )
    end if
  end do

  m1 = 0
  m2 = m
  n1 = 0
  n2 = n
  base = gap_open

  call ss_lg_fsq ( a, b, m, m1, m2, n, n1, n2, simple_score, gap_open, &
    gap_extend, base, lds, s, e, f, t )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SS_LG_FSQ:'
  write ( *, '(a)' ) ' '
  write ( *, '(3x,13i5)' ) ( j, j = n1, n2 )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SF:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13f5.1)' ) i, s(i,n1:n2)
  end do

  call r4mat_imax ( lds+1, m+1, n+1, s, i2, j2 )
  i2 = i2 - 1
  j2 = j2 - 1
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) &
    '  SS_LG_FSQ reports optimal matching score is ', s(i2,j2)

  call ss_lg_fpq ( i2, j2, m, m1, m2, n, n1, n2, lds, t, npath, pathi, pathj )

  call i4vec2_print ( npath, pathi, pathj, '  Matching path:' )

  call ss_lg_match_print ( a, b, m, n, npath, pathi, pathj, &
    simple_score, gap_open, gap_extend )

  return
end
subroutine test17 ( )

!*****************************************************************************80
!
!! TEST17 tests SS_LG_RPL.
!
!  Discussion:
!
!    Produce the alignment corresponding to the example.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Eugene Myers and Webb Miller,
!    Optimal Alignments in Linear Space,
!    CABIOS, volume 4, number 1, 1988, pages 11-17.
!
  implicit none

  integer, parameter :: m = 5
  integer, parameter :: n = 3
  integer, parameter :: lds = m

  character a(m)
  character ( len = m ) apack
  character b(n)
  character ( len = n ) bpack
  real, parameter :: gap_extend = -0.5
  real, parameter :: gap_open = -2.0
  integer i
  integer i1
  integer i2
  integer j1
  integer j2
  integer m1
  integer m2
  integer n1
  integer n2
  integer npath
  integer pathi(m+n+1)
  integer pathj(m+n+1)
  real, external :: simple_score

  apack = 'AGTAC'
  bpack = 'AAG'

  do i = 1, m
    a(i) = apack(i:i)
  end do

  do i = 1, n
    b(i) = bpack(i:i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST17:'
  write ( *, '(a)' ) '  SS_LG_RPL - forward path linear'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  GAP_OPEN penalty =   ', gap_open
  write ( *, '(a,g14.6)' ) '  GAP_EXTEND penalty = ', gap_extend
  write ( *, '(a)' ) '  Matching scores by SIMPLE_SCORE.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Compare with FPQ results.'
  write ( *, '(a)' ) '  Compare results with Myers/Miller calculation'

  call chvec2_print ( m, a, n, b, '  Sequences A and B:' )

  m1 = 0
  m2 = m
  n1 = 0
  n2 = n

  call ss_lg_corners ( a, b, m, m1, m2, n, n1, n2, simple_score, gap_open, &
    gap_extend, i1, j1, i2, j2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SS_LG_CORNERS:'
  write ( *, '(a,2i6)' ) '    (I1,J1) = ', i1, j1
  write ( *, '(a,2i6)' ) '    (I2,J2) = ', i2, j2

  call ss_lg_rpl ( a, b, m, i1, i2, n, j1, j2, simple_score, gap_open, &
    gap_extend, npath, pathi, pathj )

  call i4vec2_print ( npath, pathi, pathj, '  Matching path:' )

  call ss_lg_match_print ( a, b, m, n, npath, pathi, pathj, &
    simple_score, gap_open, gap_extend )

  return
end
subroutine test18 ( )

!*****************************************************************************80
!
!! TEST18 tests SS_LG_RPL.
!
!  Discussion:
!
!    Produce the alignment corresponding to the example.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Eugene Myers and Webb Miller,
!    Optimal Alignments in Linear Space,
!    CABIOS, volume 4, number 1, 1988, pages 11-17.
!
  implicit none

  integer, parameter :: m = 3
  integer, parameter :: n = 5
  integer, parameter :: lds = m

  character a(m)
  character ( len = m ) apack
  character b(n)
  character ( len = n ) bpack
  real, parameter :: gap_extend = -0.5
  real, parameter :: gap_open = -2.0
  integer i
  integer i1
  integer i2
  integer j1
  integer j2
  integer m1
  integer m2
  integer n1
  integer n2
  integer npath
  integer pathi(m+n+1)
  integer pathj(m+n+1)
  real, external :: simple_score

  apack = 'AAG'
  bpack = 'AGTAC'

  do i = 1, m
    a(i) = apack(i:i)
  end do

  do i = 1, n
    b(i) = bpack(i:i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST18:'
  write ( *, '(a)' ) '  SS_LG_RPL - forward path linear'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  GAP_OPEN penalty =   ', gap_open
  write ( *, '(a,g14.6)' ) '  GAP_EXTEND penalty = ', gap_extend
  write ( *, '(a)' ) '  Matching scores by SIMPLE_SCORE.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We simply switched the two sequences.'
  write ( *, '(a)' ) '  Compare with unswitched results.'

  call chvec2_print ( m, a, n, b, '  Sequences A and B:' )

  m1 = 0
  m2 = m
  n1 = 0
  n2 = n

  call ss_lg_corners ( a, b, m, m1, m2, n, n1, n2, simple_score, gap_open, &
    gap_extend, i1, j1, i2, j2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SS_LG_CORNERS:'
  write ( *, '(a,2i6)' ) '    (I1,J1) = ', i1, j1
  write ( *, '(a,2i6)' ) '    (I2,J2) = ', i2, j2

  call ss_lg_rpl ( a, b, m, i1, i2, n, j1, j2, simple_score, gap_open, &
    gap_extend, npath, pathi, pathj )

  call i4vec2_print ( npath, pathi, pathj, '  Matching path:' )

  call ss_lg_match_print ( a, b, m, n, npath, pathi, pathj, &
    simple_score, gap_open, gap_extend )

  return
end
subroutine test19 ( )

!*****************************************************************************80
!
!! TEST19 tests SS_LG_BSQ, SS_LG_BPQ.
!
!  Discussion:
!
!    Produce the alignment corresponding to the example.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Eugene Myers and Webb Miller,
!    Optimal Alignments in Linear Space,
!    CABIOS, volume 4, number 1, 1988, pages 11-17.
!
  implicit none

  integer, parameter :: m = 5
  integer, parameter :: n = 3
  integer, parameter :: lds = m

  character a(m)
  character ( len = m ) apack
  character b(n)
  real base
  character ( len = n ) bpack
  real e(0:lds,0:n)
  real f(0:lds,0:n)
  real, parameter :: gap_extend = -0.5
  real, parameter :: gap_open = -2.0
  integer i
  integer i1
  integer j1
  integer m1
  integer m2
  integer n1
  integer n2
  integer npath
  integer pathi(m+n+1)
  integer pathj(m+n+1)
  real s(0:lds,0:n)
  real, external :: simple_score
  integer t(0:lds,0:n)

  apack = 'AGTAC'
  bpack = 'AAG'

  do i = 1, m
    a(i) = apack(i:i)
  end do

  do i = 1, n
    b(i) = bpack(i:i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST19:'
  write ( *, '(a)' ) '  SS_LG_BSQ - backward score quadratic'
  write ( *, '(a)' ) '  SS_LG_BPQ - backward path quadratic'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  GAP_OPEN penalty =   ', gap_open
  write ( *, '(a,g14.6)' ) '  GAP_EXTEND penalty = ', gap_extend
  write ( *, '(a)' ) '  Matching scores by SIMPLE_SCORE.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Compare path and score with FSQ/FPQ computation.'

  call chvec2_print ( m, a, n, b, '  Sequences A and B:' )

  m1 = 0
  m2 = m
  n1 = 0
  n2 = n
  base = gap_open

  call ss_lg_bsq ( a, b, m, m1, m2, n, n1, n2, simple_score, gap_open, &
    gap_extend, base, lds, s, e, f, t )

  call r4mat_imax ( lds+1, m+1, n+1, s, i1, j1 )
  i1 = i1 - 1
  j1 = j1 - 1

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) &
    '  SS_LG_BSQ reports optimal matching score is ', s(i1,j1)

  call ss_lg_bpq ( i1, j1, m, m1, m2, n, n1, n2, lds, t, npath, pathi, pathj )

  call i4vec2_print ( npath, pathi, pathj, '  Matching path:' )

  call ss_lg_match_print ( a, b, m, n, npath, pathi, pathj, &
    simple_score, gap_open, gap_extend )

  return
end
subroutine test20 ( )

!*****************************************************************************80
!
!! TEST20 compares the quadratic and linear alignment codes.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: m = 12
  integer, parameter :: n_max = 15
  integer, parameter :: lds = m
  integer, parameter :: test_max = 100

  character a(m)
  character ( len = m ) apack
  character b(n_max)
  real base
  character ( len = n_max ) bpack
  real e(0:lds,0:n_max)
  real f(0:lds,0:n_max)
  real, parameter :: gap_extend = -0.5
  real, parameter :: gap_open = -2.0
  integer i
  integer i1
  integer i2
  integer iseed
  integer j1
  integer j2
  integer m1
  integer m2
  integer n
  integer n_bad
  integer n_good
  integer n1
  integer n2
  integer npath_l
  integer npath_q
  real, external :: pam120_score
  integer pathi_l(m+n_max+1)
  integer pathi_q(m+n_max+1)
  integer pathj_l(m+n_max+1)
  integer pathj_q(m+n_max+1)
  real s(0:lds,0:n_max)
  real score_l
  real score_q
  integer t(0:lds,0:n_max)
  integer test_num

  apack = 'GCTGATATAGCT'

  do i = 1, m
    a(i) = apack(i:i)
  end do

  n_bad = 0
  n_good = 0

  call get_seed ( iseed )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST20:'
  write ( *, '(a)' ) '  SS_LG_FSQ - Forward score quadratic;'
  write ( *, '(a)' ) '  SS_LG_FPQ - Forward path quadratic;'
  write ( *, '(a)' ) '  SS_LG_RPL - Recursive path linear;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  GAP_OPEN penalty =   ', gap_open
  write ( *, '(a,g14.6)' ) '  GAP_EXTEND penalty = ', gap_extend
  write ( *, '(a)' ) '  Matching scores by PAM120_SCORE.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Compare the score computed by FSQ/FPQ with'
  write ( *, '(a)' ) '  the score associated with the path determined'
  write ( *, '(a)' ) '  by RPL.  If the scores don''t match, the paths'
  write ( *, '(a)' ) '  differ, and presumably, the RPL algorithm has failed.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The test is carried out by comparing a fixed sequence'
  write ( *, '(a)' ) '  with many "mutated" variations.'

  do test_num = 1, test_max
!
!  Generate a random mutation of the A sequence.
!
    b(1:m) = a(1:m)
    n = m

    call mutate ( n_max, n, b, iseed )
!
!  Align using the quadratic code.
!
    m1 = 0
    m2 = m
    n1 = 0
    n2 = n
    base = gap_open

    call ss_lg_fsq ( a, b, m, m1, m2, n, n1, n2, pam120_score, gap_open, &
      gap_extend, base, lds, s, e, f, t )

    call r4mat_imax ( lds+1, m+1, n+1, s, i2, j2 )
    i2 = i2 - 1
    j2 = j2 - 1

    call ss_lg_fpq ( i2, j2, m, m1, m2, n, n1, n2, lds, t, npath_q, pathi_q, &
      pathj_q )

    score_q = s(i2,j2)
!
!  Align using the linear code.
!
    m1 = 0
    m2 = m
    n1 = 0
    n2 = n

    call ss_lg_corners ( a, b, m, m1, m2, n, n1, n2, pam120_score, gap_open, &
      gap_extend, i1, j1, i2, j2 )

    call ss_lg_rpl ( a, b, m, i1, i2, n, j1, j2, pam120_score, gap_open, &
      gap_extend, npath_l, pathi_l, pathj_l )

    call ss_lg_match_score ( a, b, m, n, npath_l, pathi_l, pathj_l, &
      pam120_score, gap_open, gap_extend, score_l )
!
!  If the scores differ, complain.
!
    if ( score_q == score_l ) then
      n_good = n_good + 1
    else
      n_bad = n_bad + 1
    end if

    if ( score_q /= score_l ) then

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Linear/quadratic score discrepancy!'
      write ( *, '(a)' ) ' '

      call chvec_print ( n, b, '  Mutated sequence B:' )

      call i4vec2_print ( npath_q, pathi_q, pathj_q, &
        '  Quadratic matching path:' )

      call ss_lg_match_print ( a, b, m, n, npath_q, pathi_q, pathj_q, &
        pam120_score, gap_open, gap_extend )

      call i4vec2_print ( npath_l, pathi_l, pathj_l, '  Linear matching path:' )

      call ss_lg_match_print ( a, b, m, n, npath_l, pathi_l, pathj_l, &
        pam120_score, gap_open, gap_extend )

    else if ( mod ( 10 * test_num, test_max ) == 0 ) then

      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  Sample test number ', test_num
      write ( *, '(a)' ) ' '

      call ss_lg_match_print ( a, b, m, n, npath_l, pathi_l, pathj_l, &
        pam120_score, gap_open, gap_extend )

    end if

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of score agreements =    ', n_good
  write ( *, '(a,i6)' ) '  Number of score disagreements = ', n_bad

  return
end
subroutine test21 ( )

!*****************************************************************************80
!
!! TEST21 tests SS_LG_FSQ, SS_LG_FPQ.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: m = 60
  integer, parameter :: n = 17
  integer, parameter :: lds = m

  character a(m)
  character ( len = m ) apack
  character b(n)
  real base
  character ( len = n ) bpack
  real e(0:lds,0:n)
  real f(0:lds,0:n)
  real, parameter :: gap_extend = -0.5
  real, parameter :: gap_open = -2.0
  integer i
  integer i2
  integer j2
  integer m1
  integer m2
  integer n1
  integer n2
  integer npath
  real, external :: pam120_score
  integer pathi(m+n+1)
  integer pathj(m+n+1)
  real s(0:lds,0:n)
  integer t(0:lds,0:n)

  apack = 'MMAAEAGGEEGGPVTAGAAGGGAAAASGAYPAVCRVKIPAALPVAAAAPFPGLAEAGVAA'
  bpack = 'METGPAVLGVKIPAGTM'

  do i = 1, m
    a(i) = apack(i:i)
  end do

  do i = 1, n
    b(i) = bpack(i:i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST21:'
  write ( *, '(a)' ) '  SS_LG_FSQ - Forward score quadratic;'
  write ( *, '(a)' ) '  SS_LG_FPQ - Forward path quadratic.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  GAP_OPEN penalty =   ', gap_open
  write ( *, '(a,g14.6)' ) '  GAP_EXTEND penalty = ', gap_extend
  write ( *, '(a)' ) '  Matching scores by PAM120_SCORE.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The second sequence is much shorter than the first.'

  m1 = 0
  m2 = m
  n1 = 0
  n2 = n
  base = gap_open

  call ss_lg_fsq ( a, b, m, m1, m2, n, n1, n2, pam120_score, gap_open, &
    gap_extend, base, lds, s, e, f, t )

  call r4mat_imax ( lds+1, m+1, n+1, s, i2, j2 )
  i2 = i2 - 1
  j2 = j2 - 1
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) &
    '  SS_LG_FSQ reports optimal matching score is ', s(i2,j2)

  call ss_lg_fpq ( i2, j2, m, m1, m2, n, n1, n2, lds, t, npath, pathi, pathj )

  call ss_lg_match_print ( a, b, m, n, npath, pathi, pathj, &
    pam120_score, gap_open, gap_extend )

  return
end
subroutine test22 ( )

!*****************************************************************************80
!
!! TEST22 tests SS_LG_RPL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: m = 60
  integer, parameter :: n = 17

  character a(m)
  character ( len = m ) apack
  character b(n)
  character ( len = n ) bpack
  real, parameter :: gap_extend = -0.5
  real, parameter :: gap_open = -2.0
  integer i
  integer i1
  integer i2
  integer j1
  integer j2
  integer m1
  integer m2
  integer n1
  integer n2
  integer npath
  integer pathi(m+n+1)
  integer pathj(m+n+1)
  real, external :: pam120_score

  apack = 'MMAAEAGGEEGGPVTAGAAGGGAAAASGAYPAVCRVKIPAALPVAAAAPFPGLAEAGVAA'
  bpack = 'METGPAVLGVKIPAGTM'

  do i = 1, m
    a(i) = apack(i:i)
  end do

  do i = 1, n
    b(i) = bpack(i:i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST22:'
  write ( *, '(a)' ) '  SS_LG_RPL - path routine;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  GAP_OPEN penalty =   ', gap_open
  write ( *, '(a,g14.6)' ) '  GAP_EXTEND penalty = ', gap_extend
  write ( *, '(a)' ) '  Matching scores by PAM120_SCORE.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The second sequence is much shorter than the first.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Compare with FSQ/FPQ calculation.'

  m1 = 0
  m2 = m
  n1 = 0
  n2 = n

  call ss_lg_corners ( a, b, m, m1, m2, n, n1, n2, pam120_score, gap_open, &
    gap_extend, i1, j1, i2, j2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SS_LG_CORNERS:'
  write ( *, '(a,2i6)' ) '    (I1,J1) = ', i1, j1
  write ( *, '(a,2i6)' ) '    (I2,J2) = ', i2, j2

  call ss_lg_rpl ( a, b, m, i1, i2, n, j1, j2, pam120_score, gap_open, &
    gap_extend, npath, pathi, pathj )

  call ss_lg_match_print ( a, b, m, n, npath, pathi, pathj, &
    pam120_score, gap_open, gap_extend )

  return
end
subroutine test23 ( )

!*****************************************************************************80
!
!! TEST23 tests SS_LG_RPL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: m = 17
  integer, parameter :: n = 60

  character a(m)
  character ( len = m ) apack
  character b(n)
  character ( len = n ) bpack
  real, parameter :: gap_extend = -0.5
  real, parameter :: gap_open = -2.0
  integer i
  integer i1
  integer i2
  integer j1
  integer j2
  integer m1
  integer m2
  integer n1
  integer n2
  integer npath
  integer pathi(m+n+1)
  integer pathj(m+n+1)
  real, external :: pam120_score

  apack = 'METGPAVLGVKIPAGTM'
  bpack = 'MMAAEAGGEEGGPVTAGAAGGGAAAASGAYPAVCRVKIPAALPVAAAAPFPGLAEAGVAA'

  do i = 1, m
    a(i) = apack(i:i)
  end do

  do i = 1, n
    b(i) = bpack(i:i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST23:'
  write ( *, '(a)' ) '  SS_LG_RPL - path routine;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  GAP_OPEN penalty =   ', gap_open
  write ( *, '(a,g14.6)' ) '  GAP_EXTEND penalty = ', gap_extend
  write ( *, '(a)' ) '  Matching scores by PAM120_SCORE.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Switch the sequences of the last example.'
  write ( *, '(a)' ) '  The first sequence is much shorter than the second.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Compare with previous and FSQ/FPQ calculations.'

  m1 = 0
  m2 = m
  n1 = 0
  n2 = n

  call ss_lg_corners ( a, b, m, m1, m2, n, n1, n2, pam120_score, gap_open, &
    gap_extend, i1, j1, i2, j2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SS_LG_CORNERS:'
  write ( *, '(a,2i6)' ) '    (I1,J1) = ', i1, j1
  write ( *, '(a,2i6)' ) '    (I2,J2) = ', i2, j2

  call ss_lg_rpl ( a, b, m, i1, i2, n, j1, j2, pam120_score, gap_open, &
    gap_extend, npath, pathi, pathj )

  call ss_lg_match_print ( a, b, m, n, npath, pathi, pathj, &
    pam120_score, gap_open, gap_extend )

  return
end
