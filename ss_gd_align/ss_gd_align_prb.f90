program main

!*****************************************************************************80
!
!! MAIN is the main program for SS_GD_ALIGN_PRB.
!
!  Discussion:
!
!    SS_GD_ALIGN_PRB carries out some alignment tests.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SS_GD_ALIGN_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the SS_GD_ALIGN library.'

  call test01 ( )
  call test02 ( )
  call test04 ( )
  call test05 ( )
  call test06 ( )
  call test07 ( )
  call test08 ( )
  call test09 ( )
  call test10 ( )
  call test11 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SS_GD_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests SS_GD_FSQ, SS_GD_FSL, SS_GD_BSQ, SS_GD_BSL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Michael Waterman,
!    Introduction to Computational Biology,
!    Chapman and Hall, 1995, Table 9.2, page 194.
!
  implicit none

  integer, parameter :: m = 12
  integer, parameter :: n = 12
  integer, parameter :: lds = m

  character a(m)
  character ( len = m ) apack
  character b(n)
  character ( len = m ) bpack
  real e(0:lds,0:n)
  real ev(0:n)
  real f(0:lds,0:n)
  real fv(0:n)
  real, parameter :: gap_extend = -0.5
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
  real, external :: test01_score
!
  apack = 'GCTGATATAGCT'

  do i = 1, m
    a(i) = apack(i:i)
  end do

  bpack = 'GGGTGATTAGCT'

  do i = 1, n
    b(i) = bpack(i:i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01:'
  write ( *, '(a)' ) '  Sequence-Sequence Global Distance'
  write ( *, '(a)' ) '  Forward/backward score Quadratic/Linear'
  write ( *, '(a)' ) '  alignment routines:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SS_GD_FSQ - Forward quadratic;'
  write ( *, '(a)' ) '  SS_GD_FSL - Forward linear;'
  write ( *, '(a)' ) '  SS_GD_BSQ - Backward quadratic;'
  write ( *, '(a)' ) '  SS_GD_BSL - Backward linear.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Gap extend penalty = ', gap_extend

  m1 = 0
  m2 = m
  n1 = 0
  n2 = n

  call chvec2_print ( m, a, n, b, '  Sequences A and B:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  TEST01 Scoring Matrix:'
  write ( *, '(a)' ) ' '
  write ( *, '(3x,13i5)' ) ( j, j = n1, n2 )
  write ( *, '(a)' ) ' '
    write ( *, '(i3,13f5.1)' ) i, test01_score ( 'X', 'X' ), &
      ( test01_score( 'X', b(j) ), j = 1, n )
  do i = 1, m
    write ( *, '(i3,13f5.1)' ) i, test01_score( a(i), 'X' ), &
      ( test01_score( a(i), b(j) ), j = 1, n )
  end do

  call ss_gd_fsq ( a, b, m, m1, m2, n, n1, n2, test01_score, gap_extend, &
    lds, s, e, f, t )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SS_GD_FSQ:'
  write ( *, '(a)' ) ' '
  write ( *, '(3x,13i5)' ) ( j, j = n1, n2 )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SF:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13f5.1)' ) i, s(i,n1:n2)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  EF:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13f5.1)' ) i, e(i,n1:n2)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  FF:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13f5.1)' ) i, f(i,n1:n2)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  TF:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13i5)' ) i, t(i,n1:n2)
  end do

  do m2 = 0, m

    call ss_gd_fsl ( a, b, m, m1, m2, n, n1, n2, test01_score, gap_extend, &
      sv, ev, fv, tv )

    s(m2,n1:n2) = sv(n1:n2)
    e(m2,n1:n2) = ev(n1:n2)
    f(m2,n1:n2) = fv(n1:n2)
    t(m2,n1:n2) = tv(n1:n2)

  end do

  m2 = m

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SS_GD_FSL:'
  write ( *, '(a)' ) ' '
  write ( *, '(3x,13i5)' ) ( j, j = n1, n2 )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SF:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13f5.1)' ) i, s(i,n1:n2)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  EF:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13f5.1)' ) i, e(i,n1:n2)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  FF:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13f5.1)' ) i, f(i,n1:n2)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  TF:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13i5)' ) i, t(i,n1:n2)
  end do

  call ss_gd_bsq ( a, b, m, m1, m2, n, n1, n2, test01_score, gap_extend, &
    lds, s, e, f, t )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SS_GD_BSQ:'
  write ( *, '(a)' ) ' '
  write ( *, '(3x,13i5)' ) ( j, j = n1, n2 )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SB:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13f5.1)' ) i, s(i,n1:n2)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  EB:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13f5.1)' ) i, e(i,n1:n2)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  FB:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13f5.1)' ) i, f(i,n1:n2)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  TB:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13i5)' ) i, t(i,n1:n2)
  end do

  do m1 = 0, m2

    call ss_gd_bsl ( a, b, m, m1, m2, n, n1, n2, test01_score, gap_extend, &
      sv, ev, fv, tv )

    s(m1,n1:n2) = sv(n1:n2)
    e(m1,n1:n2) = ev(n1:n2)
    f(m1,n1:n2) = fv(n1:n2)
    t(m1,n1:n2) = tv(n1:n2)

  end do

  m1 = 0
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SS_GD_BSL:'
  write ( *, '(a)' ) ' '
  write ( *, '(3x,13i5)' ) ( j, j = n1, n2 )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SB:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13f5.1)' ) i, s(i,n1:n2)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  EB:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13f5.1)' ) i, e(i,n1:n2)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  FB:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13f5.1)' ) i, f(i,n1:n2)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  TB:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13i5)' ) i, t(i,n1:n2)
  end do

  return
end
function test01_score ( c1, c2 )

!*****************************************************************************80
!
!! TEST01_SCORE computes a single entry sequence/sequence matching score.
!
!  Discussion:
!
!    For each possible pair of characters C1 from sequence A and
!    C2 from sequence B, this routine returns a matching score,
!    which can be thought of as specifying how similar the two
!    characters are.  It is not necessary for the matching score
!    to be symmetric.  The character '-' is commonly used to signify
!    a gap.
!
!    The scoring matrix is:
!
!      A  C  G  T  X
!  A   0 -1 -1 -1 -1
!  C  -1  0 -1 -1 -1
!  G  -1 -1  0 -1 -1
!  T  -1 -1 -1  0 -1
!  X  -1 -1 -1 -1 -1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C1, C2, two characters to be matched.
!    C1 is from sequence A, C2 from sequence B.
!
!    Output, real TEST01_SCORE, the score for matching the two characters.
!
  implicit none

  character c1
  character c2
  real score
  real test01_score

       if ( c1 == 'A' ) then

         if ( c2 == 'A' ) then
      score = 0.0
    else if ( c2 == 'C' ) then
      score = -1.0
    else if ( c2 == 'G' ) then
      score = -1.0
    else if ( c2 == 'T' ) then
      score = -1.0
    else
      score = -1.0
    end if

  else if ( c1 == 'C' ) then

         if ( c2 == 'A' ) then
      score = -1.0
    else if ( c2 == 'C' ) then
      score = 0.0
    else if ( c2 == 'G' ) then
      score = -1.0
    else if ( c2 == 'T' ) then
      score = -1.0
    else
      score = -1.0
    end if

  else if ( c1 == 'G' ) then

         if ( c2 == 'A' ) then
      score = -1.0
    else if ( c2 == 'C' ) then
      score = -1.0
    else if ( c2 == 'G' ) then
      score = 0.0
    else if ( c2 == 'T' ) then
      score = -1.0
    else
      score = -1.0
    end if

  else if ( c1 == 'T' ) then

         if ( c2 == 'A' ) then
      score = -1.0
    else if ( c2 == 'C' ) then
      score = -1.0
    else if ( c2 == 'G' ) then
      score = -1.0
    else if ( c2 == 'T' ) then
      score = 0.0
    else
      score = -1.0
    end if

  else

    score = -1.0

  end if

  test01_score = score

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests SS_GD_FSQ, SS_GD_FSL, SS_GD_BSQ, SS_GD_BSL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Michael Waterman,
!    Introduction to Computational Biology,
!    Chapman and Hall, 1995, Table 9.3, page 199.
!
  implicit none

  integer, parameter :: m = 12
  integer, parameter :: n = 12
  integer, parameter :: lds = m

  character a(m)
  character ( len = m ) apack
  character b(n)
  character ( len = m ) bpack
  real e(0:lds,0:n)
  real ev(0:n)
  real f(0:lds,0:n)
  real fv(0:n)
  real, parameter :: gap_extend = -0.5
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
  real, external :: test015_score

  apack = 'GCTGATATAGCT'

  do i = 1, m
    a(i) = apack(i:i)
  end do

  bpack = 'GGGTGATTAGCT'

  do i = 1, n
    b(i) = bpack(i:i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02:'
  write ( *, '(a)' ) '  Sequence-Sequence Global Distance'
  write ( *, '(a)' ) '  Forward/backward score Quadratic/Linear'
  write ( *, '(a)' ) '  alignment routines:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SS_GD_FSQ - Forward quadratic;'
  write ( *, '(a)' ) '  SS_GD_FSL - Forward linear;'
  write ( *, '(a)' ) '  SS_GD_BSQ - Backward quadratic;'
  write ( *, '(a)' ) '  SS_GD_BSL - Backward linear.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Gap extend penalty = ', gap_extend

  m1 = 0
  m2 = m
  n1 = 0
  n2 = n

  call chvec2_print ( m, a, n, b, '  Sequences A and B:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  TEST015 Scoring Matrix:'
  write ( *, '(a)' ) ' '
  write ( *, '(3x,13i5)' ) ( j, j = n1, n2 )
  write ( *, '(a)' ) ' '
    write ( *, '(i3,13f5.1)' ) i, test015_score ( 'X', 'X' ), &
      ( test015_score( 'X', b(j) ), j = 1, n )
  do i = 1, m
    write ( *, '(i3,13f5.1)' ) i, test015_score( a(i), 'X' ), &
      ( test015_score( a(i), b(j) ), j = 1, n )
  end do

  call ss_gd_fsq ( a, b, m, m1, m2, n, n1, n2, test015_score, gap_extend, &
    lds, s, e, f, t )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SS_GD_FSQ:'
  write ( *, '(a)' ) ' '
  write ( *, '(3x,13i5)' ) ( j, j = n1, n2 )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SF:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13f5.1)' ) i, s(i,n1:n2)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  EF:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13f5.1)' ) i, e(i,n1:n2)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  FF:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13f5.1)' ) i, f(i,n1:n2)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  TF:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13i5)' ) i, t(i,n1:n2)
  end do

  do m2 = 0, m

    call ss_gd_fsl ( a, b, m, m1, m2, n, n1, n2, test015_score, gap_extend, &
      sv, ev, fv, tv )

    s(m2,n1:n2) = sv(n1:n2)
    e(m2,n1:n2) = ev(n1:n2)
    f(m2,n1:n2) = fv(n1:n2)
    t(m2,n1:n2) = tv(n1:n2)

  end do

  m2 = m

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SS_GD_FSL:'
  write ( *, '(a)' ) ' '
  write ( *, '(3x,13i5)' ) ( j, j = n1, n2 )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SF:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13f5.1)' ) i, s(i,n1:n2)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  EF:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13f5.1)' ) i, e(i,n1:n2)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  FF:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13f5.1)' ) i, f(i,n1:n2)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  TF:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13i5)' ) i, t(i,n1:n2)
  end do

  call ss_gd_bsq ( a, b, m, m1, m2, n, n1, n2, test015_score, gap_extend, &
    lds, s, e, f, t )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SS_GD_BSQ:'
  write ( *, '(a)' ) ' '
  write ( *, '(3x,13i5)' ) ( j, j = n1, n2 )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SB:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13f5.1)' ) i, s(i,n1:n2)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  EB:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13f5.1)' ) i, e(i,n1:n2)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  FB:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13f5.1)' ) i, f(i,n1:n2)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  TB:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13i5)' ) i, t(i,n1:n2)
  end do

  do m1 = 0, m2

    call ss_gd_bsl ( a, b, m, m1, m2, n, n1, n2, test015_score, gap_extend, &
      sv, ev, fv, tv )

    s(m1,n1:n2) = sv(n1:n2)
    e(m1,n1:n2) = ev(n1:n2)
    f(m1,n1:n2) = fv(n1:n2)
    t(m1,n1:n2) = tv(n1:n2)

  end do

  m1 = 0
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SS_GD_BSL:'
  write ( *, '(a)' ) ' '
  write ( *, '(3x,13i5)' ) ( j, j = n1, n2 )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SB:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13f5.1)' ) i, s(i,n1:n2)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  EB:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13f5.1)' ) i, e(i,n1:n2)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  FB:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13f5.1)' ) i, f(i,n1:n2)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  TB:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13i5)' ) i, t(i,n1:n2)
  end do

  return
end
function test015_score ( c1, c2 )

!*****************************************************************************80
!
!! TEST015_SCORE computes a single entry sequence/sequence matching score.
!
!  Discussion:
!
!    For each possible pair of characters C1 from sequence A and
!    C2 from sequence B, this routine returns a matching score,
!    which can be thought of as specifying how similar the two
!    characters are.  It is not necessary for the matching score
!    to be symmetric.  The character '-' is commonly used to signify
!    a gap.
!
!    The scoring matrix is:
!
!      A  C  G  T  X
!  A   1 -1 -1 -1 -2
!  C  -1  1 -1 -1 -2
!  G  -1 -1  1 -1 -2
!  T  -1 -1 -1  1 -2
!  X  -2 -2 -2 -2 -2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C1, C2, two characters to be matched.
!    C1 is from sequence A, C2 from sequence B.
!
!    Output, real TEST015_SCORE, the score for matching the two characters.
!
  implicit none

  character c1
  character c2
  real score
  real test015_score

       if ( c1 == 'A' ) then

         if ( c2 == 'A' ) then
      score = 1.0
    else if ( c2 == 'C' ) then
      score = -1.0
    else if ( c2 == 'G' ) then
      score = -1.0
    else if ( c2 == 'T' ) then
      score = -1.0
    else
      score = -2.0
    end if

  else if ( c1 == 'C' ) then

         if ( c2 == 'A' ) then
      score = -1.0
    else if ( c2 == 'C' ) then
      score = 1.0
    else if ( c2 == 'G' ) then
      score = -1.0
    else if ( c2 == 'T' ) then
      score = -1.0
    else
      score = -2.0
    end if

  else if ( c1 == 'G' ) then

         if ( c2 == 'A' ) then
      score = -1.0
    else if ( c2 == 'C' ) then
      score = -1.0
    else if ( c2 == 'G' ) then
      score = 1.0
    else if ( c2 == 'T' ) then
      score = -1.0
    else
      score = -2.0
    end if

  else if ( c1 == 'T' ) then

         if ( c2 == 'A' ) then
      score = -1.0
    else if ( c2 == 'C' ) then
      score = -1.0
    else if ( c2 == 'G' ) then
      score = -1.0
    else if ( c2 == 'T' ) then
      score = 1.0
    else
      score = -2.0
    end if

  else

    score = -2.0

  end if

  test015_score = score

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests SS_GD_FSQ, SS_GD_FPQ.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 January 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Michael Waterman,
!    Introduction to Computational Biology,
!    Chapman and Hall, 1995, page 199.
!
  implicit none

  integer, parameter :: m = 12
  integer, parameter :: n = 12
  integer, parameter :: lds = m

  character a(m)
  character ( len = m ) apack
  character b(n)
  character ( len = n ) bpack
  real e(0:lds,0:n)
  real f(0:lds,0:n)
  real, parameter :: gap_extend = -0.5
  integer i
  integer m1
  integer m2
  integer n1
  integer n2
  integer npath
  integer pathi(m+n+1)
  integer pathj(m+n+1)
  real s(0:lds,0:n)
  integer t(0:lds,0:n)
  real, external :: test015_score

  apack = 'GCTGATATAGCT'

  do i = 1, m
    a(i) = apack(i:i)
  end do

  bpack = 'GGGTGATTAGCT'

  do i = 1, n
    b(i) = bpack(i:i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04:'
  write ( *, '(a)' ) '  Sequence-Sequence Global Distance'
  write ( *, '(a)' ) '  Forward score and path quadratic alignment routines:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SS_GD_FSQ - score routine;'
  write ( *, '(a)' ) '  SS_GD_FPQ - path routine.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Gap extend penalty = ', gap_extend

  m1 = 0
  m2 = m
  n1 = 0
  n2 = n

  call chvec2_print ( m, a, n, b, '  Sequences A and B:' )

  call ss_gd_fsq ( a, b, m, m1, m2, n, n1, n2, test015_score, gap_extend, &
    lds, s, e, f, t )

  call ss_gd_fpq ( m, m1, m2, n, n1, n2, lds, t, npath, pathi, pathj )

  call i4vec2_print ( npath, pathi, pathj, '  Matching path:' )

  call ss_gd_match_print ( a, b, m, n, npath, pathi, pathj, test015_score, &
    gap_extend )

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 tests SS_GD_RPL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 January 2009
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
  integer i
  integer npath
  integer pathi(m+n+1)
  integer pathj(m+n+1)
  real, external :: test015_score

  apack = 'GCTGATATAGCT'

  do i = 1, m
    a(i) = apack(i:i)
  end do

  bpack = 'GGGTGATTAGCT'

  do i = 1, n
    b(i) = bpack(i:i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05:'
  write ( *, '(a)' ) '  Sequence-Sequence Global Distance'
  write ( *, '(a)' ) '  Recursive Path Linear alignment routine:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SS_GD_RPL - path routine;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Gap extend penalty = ', gap_extend

  call chvec2_print ( m, a, n, b, '  Sequences A and B:' )

  call ss_gd_rpl ( a, b, m, n, test015_score, gap_extend, npath, pathi, pathj )

  call i4vec2_print ( npath, pathi, pathj, '  Matching path:' )

  call ss_gd_match_print ( a, b, m, n, npath, pathi, pathj, test015_score, &
    gap_extend )

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 tests SS_GD_FSQ, SS_GD_FSL, SS_GD_BSQ, SS_GD_BSL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 January 2009
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
  character ( len = n ) bpack
  real e(0:lds,0:n)
  real ev(0:n)
  real f(0:lds,0:n)
  real fv(0:n)
  real, parameter :: gap_extend = -0.5
  integer i
  integer j
  integer m1
  integer m2
  integer n1
  integer n2
  real, external :: pam120_score
  real s(0:lds,0:n)
  real sv(0:n)
  integer t(0:lds,0:n)
  integer tv(0:n)

  apack = 'GCTGATATAGCT'

  do i = 1, m
    a(i) = apack(i:i)
  end do

  bpack = 'GGGTGATTAGCT'

  do i = 1, n
    b(i) = bpack(i:i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06:'
  write ( *, '(a)' ) '  Using PAM120 scoring matrix.'
  write ( *, '(a)' ) '  Sequence-Sequence Global Distance'
  write ( *, '(a)' ) '  Forward Score Quadratic/Linear alignment routines:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SS_GD_FSQ - Forward quadratic;'
  write ( *, '(a)' ) '  SS_GD_FSL - Forward linear;'
  write ( *, '(a)' ) '  SS_GD_BSQ - Backward quadratic;'
  write ( *, '(a)' ) '  SS_GD_BSL - Backward linear.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Gap extend penalty = ', gap_extend

  m1 = 0
  m2 = m
  n1 = 0
  n2 = n

  call chvec2_print ( m, a, n, b, '  Sequences A and B:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  PAM120 Scoring Matrix:'
  write ( *, '(a)' ) ' '
  write ( *, '(3x,13i5)' ) ( j, j = n1, n2 )
  write ( *, '(a)' ) ' '
    write ( *, '(i3,13f5.1)' ) i, pam120_score ( 'X', 'X' ), &
      ( pam120_score( 'X', b(j) ), j = 1, n )
  do i = 1, m
    write ( *, '(i3,13f5.1)' ) i, pam120_score( a(i), 'X' ), &
      ( pam120_score( a(i), b(j) ), j = 1, n )
  end do

  call ss_gd_fsq ( a, b, m, m1, m2, n, n1, n2, pam120_score, gap_extend, &
    lds, s, e, f, t )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SS_GD_FSQ:'
  write ( *, '(a)' ) ' '
  write ( *, '(3x,13i5)' ) ( j, j = n1, n2 )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SF:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13f5.1)' ) i, s(i,n1:n2)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  EF:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13f5.1)' ) i, e(i,n1:n2)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  FF:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13f5.1)' ) i, f(i,n1:n2)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  TF:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13i5)' ) i, t(i,n1:n2)
  end do

  do m2 = 0, m

    call ss_gd_fsl ( a, b, m, m1, m2, n, n1, n2, pam120_score, gap_extend, &
      sv, ev, fv, tv )

    s(m2,n1:n2) = sv(n1:n2)
    e(m2,n1:n2) = ev(n1:n2)
    f(m2,n1:n2) = fv(n1:n2)
    t(m2,n1:n2) = tv(n1:n2)

  end do

  m2 = m

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SS_GD_FSL:'
  write ( *, '(a)' ) ' '
  write ( *, '(3x,13i5)' ) ( j, j = n1, n2 )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SF:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13f5.1)' ) i, s(i,n1:n2)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  EF:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13f5.1)' ) i, e(i,n1:n2)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  FF:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13f5.1)' ) i, f(i,n1:n2)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  TF:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13i5)' ) i, t(i,n1:n2)
  end do

  call ss_gd_bsq ( a, b, m, m1, m2, n, n1, n2, pam120_score, gap_extend, &
    lds, s, e, f, t )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SS_GD_BSQ:'
  write ( *, '(a)' ) ' '
  write ( *, '(3x,13i5)' ) ( j, j = n1, n2 )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SB:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13f5.1)' ) i, s(i,n1:n2)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  EB:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13f5.1)' ) i, e(i,n1:n2)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  FB:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13f5.1)' ) i, f(i,n1:n2)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  TB:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13i5)' ) i, t(i,n1:n2)
  end do

  do m1 = 0, m2

    call ss_gd_bsl ( a, b, m, m1, m2, n, n1, n2, pam120_score, gap_extend, &
      sv, ev, fv, tv )

    s(m1,n1:n2) = sv(n1:n2)
    e(m1,n1:n2) = ev(n1:n2)
    f(m1,n1:n2) = fv(n1:n2)
    t(m1,n1:n2) = tv(n1:n2)

  end do

  m1 = 0
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SS_GD_BSL:'
  write ( *, '(a)' ) ' '
  write ( *, '(3x,13i5)' ) ( j, j = n1, n2 )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SB:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13f5.1)' ) i, s(i,n1:n2)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  EB:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13f5.1)' ) i, e(i,n1:n2)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  FB:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13f5.1)' ) i, f(i,n1:n2)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  TB:'
  write ( *, '(a)' ) ' '
  do i = m1, m2
    write ( *, '(i3,13i5)' ) i, t(i,n1:n2)
  end do

  return
end
subroutine test07 ( )

!*****************************************************************************80
!
!! TEST07 tests SS_GD_FSQ, SS_GD_FPQ.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 January 2009
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
  character ( len = n ) bpack
  real e(0:lds,0:n)
  real f(0:lds,0:n)
  real, parameter :: gap_extend = -0.5
  integer i
  integer m1
  integer m2
  integer n1
  integer n2
  integer npath
  integer pathi(m+n+1)
  integer pathj(m+n+1)
  real, external :: pam120_score
  real s(0:lds,0:n)
  integer t(0:lds,0:n)

  apack = 'GCTGATATAGCT'

  do i = 1, m
    a(i) = apack(i:i)
  end do

  bpack = 'GGGTGATTAGCT'

  do i = 1, n
    b(i) = bpack(i:i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07:'
  write ( *, '(a)' ) '  Using PAM120 scoring matrix.'
  write ( *, '(a)' ) '  Sequence-Sequence Global Distance'
  write ( *, '(a)' ) '  forward score and path quadratic alignment routines:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SS_GD_FSQ - score routine;'
  write ( *, '(a)' ) '  SS_GD_FPQ - path routine.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Gap extend penalty = ', gap_extend

  m1 = 0
  m2 = m
  n1 = 0
  n2 = n

  call chvec2_print ( m, a, n, b, '  Sequences A and B:' )

  call ss_gd_fsq ( a, b, m, m1, m2, n, n1, n2, pam120_score, gap_extend, &
    lds, s, e, f, t )

  call ss_gd_fpq ( m, m1, m2, n, n1, n2, lds, t, npath, pathi, pathj )

  call i4vec2_print ( npath, pathi, pathj, '  Matching path:' )

  call ss_gd_match_print ( a, b, m, n, npath, pathi, pathj, pam120_score, &
    gap_extend )

  return
end
subroutine test08 ( )

!*****************************************************************************80
!
!! TEST08 tests SS_GD_RPL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 January 2009
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
  integer i
  integer npath
  integer pathi(m+n+1)
  integer pathj(m+n+1)
  real, external :: pam120_score

  apack = 'GCTGATATAGCT'

  do i = 1, m
    a(i) = apack(i:i)
  end do

  bpack = 'GGGTGATTAGCT'

  do i = 1, n
    b(i) = bpack(i:i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08:'
  write ( *, '(a)' ) '  Using PAM120 scoring matrix.'
  write ( *, '(a)' ) '  Sequence-sequence Global Distance'
  write ( *, '(a)' ) '  Recursive Path Linear alignment routine:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SS_GG_RPL - path routine;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Gap extend penalty = ', gap_extend

  call chvec2_print ( m, a, n, b, '  Sequences A and B:' )

  call ss_gd_rpl ( a, b, m, n, pam120_score, gap_extend, npath, pathi, pathj )

  call i4vec2_print ( npath, pathi, pathj, '  Matching path:' )

  call ss_gd_match_print ( a, b, m, n, npath, pathi, pathj, pam120_score, &
    gap_extend )

  return
end
subroutine test09 ( )

!*****************************************************************************80
!
!! TEST09 tests SS_GD_RPL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 January 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: m = 60
  integer, parameter :: n = 59

  character a(m)
  character ( len = 60 ) apack
  character b(n)
  character ( len = 59 ) bpack
  real, parameter :: gap_extend = -0.5
  integer i
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
  write ( *, '(a)' ) 'TEST09:'
  write ( *, '(a)' ) '  Using PAM120 scoring matrix.'
  write ( *, '(a)' ) '  Sequence-sequence Global Distance'
  write ( *, '(a)' ) '  Recursive Path Linear alignment routines:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SS_GD_RPL - path routine;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Gap extend penalty = ', gap_extend

  call chvec2_print ( m, a, n, b, '  Sequences A and B:' )

  call ss_gd_rpl ( a, b, m, n, pam120_score, gap_extend, npath, pathi, pathj )

  call i4vec2_print ( npath, pathi, pathj, '  Matching path:' )

  call ss_gd_match_print ( a, b, m, n, npath, pathi, pathj, pam120_score, &
    gap_extend )

  return
end
subroutine test10 ( )

!*****************************************************************************80
!
!! TEST10 tests SS_GD_FSQ, SS_GD_FPQ.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 January 2009
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
  character ( len = n ) bpack
  real e(0:lds,0:n)
  real f(0:lds,0:n)
  real, parameter :: gap_extend = -0.5
  integer i
  integer m1
  integer m2
  integer n1
  integer n2
  integer npath
  integer pathi(m+n+1)
  integer pathj(m+n+1)
  real, external :: pam200_score
  real s(0:lds,0:n)
  integer t(0:lds,0:n)

  apack = 'GCTGATATAGCT'

  do i = 1, m
    a(i) = apack(i:i)
  end do

  bpack = 'GGGTGATTAGCT'

  do i = 1, n
    b(i) = bpack(i:i)
  end do

  m1 = 0
  m2 = m
  n1 = 0
  n2 = n

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST10:'
  write ( *, '(a)' ) '  Using PAM200 scoring matrix.'
  write ( *, '(a)' ) '  Sequence-sequence Global Distance'
  write ( *, '(a)' ) '  Recursive Score and Path Linear alignment routines:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SS_GD_FSQ - score routine;'
  write ( *, '(a)' ) '  SS_GD_FPQ - path routine.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Gap extend penalty = ', gap_extend

  call chvec2_print ( m, a, n, b, '  Sequences A and B:' )

  call ss_gd_fsq ( a, b, m, m1, m2, n, n1, n2, pam200_score, gap_extend, &
    lds, s, e, f, t )

  call ss_gd_fpq ( m, m1, m2, n, n1, n2, lds, t, npath, pathi, pathj )

  call i4vec2_print ( npath, pathi, pathj, '  Matching path:' )

  call ss_gd_match_print ( a, b, m, n, npath, pathi, pathj, pam200_score, &
    gap_extend )

  return
end
subroutine test11 ( )

!*****************************************************************************80
!
!! TEST11 tests SS_GD_RPL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 January 2009
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
  integer i
  integer npath
  integer pathi(m+n+1)
  integer pathj(m+n+1)
  real, external :: pam200_score

  apack = 'GCTGATATAGCT'

  do i = 1, m
    a(i) = apack(i:i)
  end do

  bpack = 'GGGTGATTAGCT'

  do i = 1, n
    b(i) = bpack(i:i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST11:'
  write ( *, '(a)' ) '  Using PAM200 scoring matrix.'
  write ( *, '(a)' ) '  Sequence-sequence Global Distance'
  write ( *, '(a)' ) '  Recursive Path Linear alignment routines:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  SS_GD_RPL - path routine;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Gap extend penalty = ', gap_extend
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Compare with SS_GD_FSQ/SS_GD_FPQ on same data.'

  call chvec2_print ( m, a, n, b, '  Sequences A and B:' )

  call ss_gd_rpl ( a, b, m, n, pam200_score, gap_extend, npath, pathi, pathj )

  call i4vec2_print ( npath, pathi, pathj, '  Matching path:' )

  call ss_gd_match_print ( a, b, m, n, npath, pathi, pathj, pam200_score, &
    gap_extend )

  return
end
