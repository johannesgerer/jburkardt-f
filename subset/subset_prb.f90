program main

!*****************************************************************************80
!
!! MAIN is the main program for SUBSET_PRB.
!
!  Discussion:
!
!    SUBSET_PRB calls the SUBSET test routines.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 August 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SUBSET_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the SUBSET library.'

  call test000 ( )
  call test001 ( )
  call test002 ( )
  call test003 ( )
  call test9001 ( )
  call test9002 ( )
  call test9003 ( )
  call test005 ( )
  call test006 ( )
  call test007 ( )
  call test008 ( )
  call test009 ( )

  call test010 ( )
  call test011 ( )
  call test012 ( )
  call test013 ( )
  call test014 ( )
  call test015 ( )
  call test016 ( )
  call test017 ( )
  call test0174 ( )
  call test0175 ( )
  call test018 ( )
  call test019 ( )

  call test020 ( )
  call test021 ( )
  call test022 ( )
  call test023 ( )
  call test024 ( )
  call test025 ( )
  call test026 ( )
  call test027 ( )
  call test028 ( )
  call test029 ( )

  call test021a ( )
  call test022a ( )
  call test023a ( )
  call test024a ( )
  call test025a ( )
  call test026a ( )
  call test026b ( )
  call test026c ( )
  call test026d ( )
  call test027a ( )
  call test028a ( )
  call test029a ( )
  call test0295 ( )

  call test0304 ( )
  call test0305 ( )
  call test031 ( )
  call test032 ( )
  call test0321 ( )
  call test0322 ( )
  call test03225 ( )
  call test0323 ( )
  call test0324 ( )
  call test0325 ( )
  call test0327 ( )
  call test058 ( )
  call test059 ( )
  call test060 ( )
  call test061 ( )
  call test0615 ( )
  call test062 ( )
  call test06225 ( )
  call test033 ( )
  call test034 ( )
  call test0625 ( )
  call test035 ( )
  call test0364 ( )
  call test0627 ( )
  call test036 ( )
  call test0365 ( )
  call test037 ( )
  call test038 ( )
  call test039 ( )

  call test040 ( )
  call test041 ( )
  call test042 ( )
  call test047 ( )
  call test043 ( )
  call test044 ( )
  call test045 ( )
  call test046 ( )
  call test048 ( )
  call test049 ( )

  call test050 ( )
  call test051 ( )
  call test052 ( )
  call test053 ( )
  call test054 ( )
  call test055 ( )
  call test056 ( )
  call test057 ( )

  call test063 ( )
  call test064 ( )
  call test065 ( )
  call test066 ( )
  call test067 ( )
  call test0675 ( )
  call test068 ( )
  call test0683 ( )
  call test0685 ( )
  call test0686 ( )
  call test0687 ( )
  call test0688 ( )
  call test0689 ( )
  call test06895 ( )
  call test069 ( )

  call test070 ( )
  call test071 ( )
  call test072 ( )
  call test073 ( )
  call test074 ( )
  call test075 ( )
  call test076 ( )
  call test077 ( )
  call test0771 ( )
  call test07715 ( )
  call test0772 ( )
  call test0773 ( )
  call test078 ( )
  call test079 ( )
  call test0795 ( )

  call test080 ( )
  call test081 ( )
  call test0813 ( )
  call test0815 ( )
  call test082 ( )
  call test083 ( )
  call test0835 ( )
  call test084 ( )
  call test085 ( )
  call test093 ( )
  call test086 ( )
  call test087 ( )
  call test088 ( )
  call test089 ( )

  call test090 ( )
  call test091 ( )
  call test092 ( )
  call test094 ( )
  call test095 ( )
  call test0955 ( )
  call test096 ( )
  call test097 ( )
  call test098 ( )
  call test099 ( )

  call test100 ( )
  call test101 ( )
  call test102 ( )
  call test103 ( )
  call test104 ( )
  call test105 ( )
  call test106 ( )
  call test107 ( )
  call test108 ( )
  call test1085 ( )
  call test109 ( )

  call test110 ( )
  call test111 ( )
  call test112 ( )
  call test113 ( )
  call test114 ( )
  call test115 ( )
  call test030 ( )
  call test1245 ( )
  call test116 ( )
  call test1163 ( )
  call test1165 ( )
  call test117 ( )
  call test118 ( )
  call test119 ( )

  call test120 ( )
  call test121 ( )
  call test1215 ( )
  call test122 ( )
  call test123 ( )
  call test124 ( )
  call test125 ( )
  call test126 ( )
  call test127 ( )
  call test128 ( )
  call test129 ( )

  call test130 ( )
  call test131 ( )
  call test132 ( )
  call test133 ( )
  call test134 ( )
  call test135 ( )
  call test136 ( )
  call test137 ( )
  call test138 ( )
  call test139 ( )
  call test1395 ( )

  call test140 ( )
  call test141 ( )
  call test142 ( )
  call test143 ( )
  call test1435 ( )
  call test144 ( )
  call test145 ( )
  call test146 ( )
  call test147 ( )
  call test1475 ( )
  call test1476 ( )
  call test1477 ( )
  call test1478 ( )
  call test148 ( )
  call test149 ( )

  call test150 ( )
  call test151 ( )
  call test152 ( )
  call test153 ( )
  call test1531 ( )
  call test0626 ( )
  call test1535 ( )
  call test1536 ( )
  call test1537 ( )
  call test154 ( )
  call test155 ( )
  call test156 ( )
  call test1565 ( )
  call test1566 ( )
  call test1567 ( )
  call test1568 ( )
  call test1569 ( )
  call test15695 ( )
  call test15696 ( )
  call test15698 ( )
  call test157 ( )
  call test158 ( )
  call test159 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SUBSET_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test000 ( )

!*****************************************************************************80
!
!! TEST000 tests RANDOM_INITIALIZE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST000'
  write ( *, '(a)' ) '  Call RANDOM_INITIALIZE to initialize the'
  write ( *, '(a)' ) '  random number generator.'

  seed = 0 
  call random_initialize ( seed )

  return
end
subroutine test001 ( )

!*****************************************************************************80
!
!! TEST001 tests ASM_ENUM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 7

  integer ( kind = 4 ) asm_num
  integer ( kind = 4 ) n

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST001'
  write ( *, '(a)' ) '  ASM_ENUM returns the number of alternating sign'
  write ( *, '(a)' ) '  matrices of a given order.'

  write ( *, '(a)' ) ' '
  do n = 0, n_max
    call asm_enum ( n, asm_num )
    write ( *, '(2x,i2,2x,i8)' ) n, asm_num
  end do

  return
end
subroutine test002 ( )

!*****************************************************************************80
!
!! TEST002 tests ASM_TRIANGLE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 7

  integer ( kind = 4 ) a(n_max+1)
  integer ( kind = 4 ) n

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST002'
  write ( *, '(a)' ) '  ASM_TRIANGLE returns a row of the alternating sign'
  write ( *, '(a)' ) '  matrix triangle.'
  write ( *, '(a)' ) ' '

  do n = 0, n_max
    call asm_triangle ( n, a )
    write ( *, '(2x,i2,2x,8i8)' ) n, a(1:n+1)
  end do

  return
end
subroutine test003 ( )

!*****************************************************************************80
!
!! TEST003 tests BELL and BELL_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) c
  integer ( kind = 4 ) c2(0:10)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST003'
  write ( *, '(a)' ) '  BELL computes Bell numbers.'
  write ( *, '(a)' ) '  BELL_VALUES returns some exact values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  N  exact C(I)  computed C(I)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bell_values ( n_data, n, c )

    if ( n_data == 0 ) then
      exit
    end if

    call bell ( n, c2 )

    write ( *, '(2x,i4,2i10)' ) n, c, c2(n)

  end do

  return
end
subroutine test9001 ( )

!*****************************************************************************80
!
!! TEST9001 tests BVEC_ADD and BVEC_SUB;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 November 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  integer ( kind = 4 ) bvec1(n)
  integer ( kind = 4 ) bvec2(n)
  integer ( kind = 4 ) bvec3(n)
  integer ( kind = 4 ) bvec4(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) :: seed = 123456789
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 10

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST9001'
  write ( *, '(a)' ) '  BVEC_ADD adds binary vectors representing integers;'
  write ( *, '(a)' ) '  BVEC_SUB subtracts binary vectors representing integers;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '        I        J        K = I + J    L = I - J'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    
    i = i4_uniform ( -100, 100, seed )
    j = i4_uniform ( -100, 100, seed )

    write ( *, '(a)' ) ' '

    write ( *, '(2x,i8,2x,i8)' ) i, j

    k = i + j
    l = i - j

    write ( *, '(a20,2x,i8,2x,i8)' ) '  Directly:         ', k, l

    call i4_to_bvec ( i, n, bvec1 )
    call i4_to_bvec ( j, n, bvec2 )

    call bvec_add ( n, bvec1, bvec2, bvec3 )
    call bvec_to_i4 ( n, bvec3, k )

    call bvec_sub ( n, bvec1, bvec2, bvec4 )
    call bvec_to_i4 ( n, bvec4, l )

    write ( *, '(a20,2x,i8,2x,i8)' ) '  BVEC_ADD, BVEC_SUB', k, l

  end do

  return
end
subroutine test9002 ( )

!*****************************************************************************80
!
!! TEST9002 tests BVEC_COMPLEMENT2;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 November 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  integer ( kind = 4 ) bvec1(n)
  integer ( kind = 4 ) bvec2(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) :: seed = 123456789
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 5

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST9002'
  write ( *, '(a)' ) '  BVEC_COMPLEMENT2 returns the two''s complement'
  write ( *, '(a)' ) '  of a (signed) binary vector;'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    
    i = i4_uniform ( -100, 100, seed )

    call i4_to_bvec ( i, n, bvec1 )

    call bvec_complement2 ( n, bvec1, bvec2 )

    call bvec_to_i4 ( n, bvec2, j )

    write ( *, '(a)' ) ' '
    write ( *, '(a,2x,i8)' ) '  I = ', i
    write ( *, '(a,2x,i8)' ) '  J = ', j
    call bvec_print ( n, bvec1, ' ' )
    call bvec_print ( n, bvec2, ' ' )

  end do

  return
end
subroutine test9003 ( )

!*****************************************************************************80
!
!! TEST9003 tests BVEC_MUL;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 15

  integer ( kind = 4 ) bvec1(n)
  integer ( kind = 4 ) bvec2(n)
  integer ( kind = 4 ) bvec3(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) :: seed = 123456789
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 10

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST9003'
  write ( *, '(a)' ) '  BVEC_MUL multiplies binary vectors '
  write ( *, '(a)' ) '  representing integers;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '        I        J        K = I * J'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    
    i = i4_uniform ( -100, 100, seed )
    j = i4_uniform ( -100, 100, seed )

    write ( *, '(a)' ) ' '

    write ( *, '(2x,i8,2x,i8)' ) i, j

    k = i * j

    write ( *, '(a20,2x,i8)' ) '  Directly:         ', k

    call i4_to_bvec ( i, n, bvec1 )
    call i4_to_bvec ( j, n, bvec2 )
    call bvec_mul ( n, bvec1, bvec2, bvec3 )
    call bvec_to_i4 ( n, bvec3, k )

    write ( *, '(a20,2x,i8)' ) '  BVEC_MUL          ', k

  end do

  return
end
subroutine test005 ( )

!*****************************************************************************80
!
!! TEST005 tests CATALAN and CATALAN_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) c
  integer ( kind = 4 ) c2(0:10)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST005'
  write ( *, '(a)' ) '  CATALAN computes Catalan numbers.'
  write ( *, '(a)' ) '  CATALAN_VALUES returns some exact values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  N  exact C(I)  computed C(I)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call catalan_values ( n_data, n, c )

    if ( n_data == 0 ) then
      exit
    end if

    call catalan ( n, c2 )

    write ( *, '(2x,i4,2i8)' ) n, c, c2(n)

  end do

  return
end
subroutine test006 ( )

!*****************************************************************************80
!
!! TEST006 tests CATALAN_ROW_NEXT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  integer ( kind = 4 ) c(0:n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ido

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST006'
  write ( *, '(a)' ) '  CATALAN_ROW_NEXT computes a row of Catalan''s triangle.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  First, compute row 7:'

  ido = 0
  i = 7
  call catalan_row_next ( ido, i, c )
  write ( *, '(2x,i2,2x,11i6)' ) i, c(0:i)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now compute rows one at a time:'
  write ( *, '(a)' ) ' '

  ido = 0

  do i = 0, n
    call catalan_row_next ( ido, i, c )
    ido = 1
    write ( *, '(2x,i2,2x,11i6)' ) i, c(0:i)
  end do

  return
end
subroutine test007 ( )

!*****************************************************************************80
!
!! TEST007 tests CFRAC_TO_RAT and RAT_TO_CFRAC.
!
!  Discussion:
!
!    Compute the continued fraction form of 4096/15625.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 10

  integer ( kind = 4 ) a(m)
  integer ( kind = 4 ) bot
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) n
  integer ( kind = 4 ) p(m)
  integer ( kind = 4 ) q(m)
  integer ( kind = 4 ) top

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST007'
  write ( *, '(a)' ) '  RAT_TO_CFRAC fraction => continued fraction,'
  write ( *, '(a)' ) '  CFRAC_TO_RAT continued fraction => fraction.'
  write ( *, '(a)' ) ' '
  top = 4096
  bot = 15625
  write ( *, '(a,i8,a,i8)' ) '  Regular fraction is ', top, ' / ', bot
 
  call rat_to_cfrac ( top, bot, m, n, a, ierror )
 
  call i4vec_print ( n, a, '  Continued fraction coefficients:' )

  call cfrac_to_rat ( n, a, p, q )
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The continued fraction convergents.'
  write ( *, '(a)' ) '  The last row contains the value of the continued'
  write ( *, '(a)' ) '  fraction, written as a common fraction.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I, P(I), Q(I), P(I)/Q(I)'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i3,2i8,g14.6)' ) i, p(i), q(i), &
      real ( p(i), kind = 8 ) / real ( q(i), kind = 8 )
  end do
 
  return
end
subroutine test008 ( )

!*****************************************************************************80
!
!! TEST008 tests CFRAC_TO_RFRAC and RFRAC_TO_CFRAC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: maxm = 10

  real ( kind = 8 ) g(2*maxm)
  real ( kind = 8 ) h(2*maxm)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) m
  real ( kind = 8 ) p(maxm)
  real ( kind = 8 ) q(maxm+1)

  m = 3

  p(1:3) = (/ 1.0D+00, 1.0D+00, 2.0D+00 /)
  q(1:4) = (/ 1.0D+00, 3.0D+00, 1.0D+00, 1.0D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST008'
  write ( *, '(a)' ) '  CFRAC_TO_RFRAC: continued fraction to ratio;'
  write ( *, '(a)' ) '  RFRAC_TO_CFRAC: ratio to continued fration.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Rational polynomial fraction coefficients:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,5f12.4)' ) '  P:  ', p(1:m)
  write ( *, '(a,5f12.4)' ) '  Q:  ', q(1:m+1)
 
  call rfrac_to_cfrac ( m, p, q, h, ierror )
 
  call r8vec_print ( 2*m, h, '  Continued fraction coefficients:' )

  g(1:2*m) = 1.0D+00

  call cfrac_to_rfrac ( 2*m, g, h, p, q )
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Recovered rational polynomial:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,5f12.4)' ) '  P:  ', p(1:m)
  write ( *, '(a,5f12.4)' ) '  Q:  ', q(1:m+1)
 
  return
end
subroutine test009 ( )

!*****************************************************************************80
!
!! TEST009 tests CHANGE_GREEDY.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: coin_num = 6

  integer ( kind = 4 ) change(100)
  integer ( kind = 4 ) change_num
  integer ( kind = 4 ), dimension ( coin_num ) :: coin_value = (/ 1, 5, 10, 25, 50, 100 /)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) total
  integer ( kind = 4 ) total2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST009'
  write ( *, '(a)' ) '  CHANGE_GREEDY makes change using the biggest'
  write ( *, '(a)' ) '  coins first.'

  total = 73

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The total for which change is to be made: ', total
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The available coins are:'
  write ( *, '(a)' ) ' '
  do i = 1, coin_num
    write ( *, '(2x,i8)' ) coin_value(i)
  end do

  call change_greedy ( total, coin_num, coin_value, change_num, change )

  write ( *, '(a)' ) ' '
  write ( *, '(2x,i8,4x,(20i3))' ) change_num, change(1:change_num)
  total2 = sum ( coin_value(change(1:change_num) ) )
  write ( *, '(2x,i8,4x,(20i3))' ) total2, coin_value(change(1:change_num) )

  return
end
subroutine test010 ( )

!*****************************************************************************80
!
!! TEST010 tests CHANGE_NEXT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: coin_num = 6

  integer ( kind = 4 ) change(100)
  integer ( kind = 4 ) change_num
  integer ( kind = 4 ), dimension ( coin_num ) :: coin_value = (/ 1, 5, 10, 25, 50, 100 /)
  logical done
  integer ( kind = 4 ) i
  integer ( kind = 4 ) total

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST010'
  write ( *, '(a)' ) '  CHANGE_NEXT displays the next possible way to make'
  write ( *, '(a)' ) '  change for a given total'

  total = 50

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The total for which change is to be made: ', total
  write ( *, '(a)' ) ' '

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The available coins are:'
  write ( *, '(a)' ) ' '
  do i = 1, coin_num
    write ( *, '(2x,i8)' ) coin_value(i)
  end do

  done = .true.
  i = 0

  do

    call change_next ( total, coin_num, coin_value, change_num, change, done )

    if ( done .or. 9 < i ) then
      exit
    end if

    i = i + 1
    write ( *, '(a)' ) ' '
    write ( *, '(i3, ":")' ) i
    write ( *, '(2x,25i3)' ) coin_value(change(1:change_num) )

  end do

  return
end
subroutine test011 ( )

!*****************************************************************************80
!
!! TEST011 tests COMB_NEXT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 April 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), allocatable :: a(:)
  logical              done
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n

  n = 5

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST011'
  write ( *, '(a)' ) '  COMB_NEXT produces combinations.'

  do k = 1, n

    allocate ( a(1:k) )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Combinations of size K = ', k
    write ( *, '(a)' ) ' '

    done = .true.

    do

      call comb_next ( n, k, a, done )
 
      if ( done ) then
        exit
      end if

      write ( *, '(2x,5i3)' ) a(1:k)

    end do

    deallocate ( a )

  end do

  return
end
subroutine test012 ( )

!*****************************************************************************80
!
!! TEST012 tests COMB_ROW.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  integer ( kind = 4 ) c(0:n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ido

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST012'
  write ( *, '(a)' ) '  COMB_ROW computes a row of Pascal''s triangle.'
  write ( *, '(a)' ) ' '
 
  ido = 0
 
  do i = 0, n
    call comb_row ( ido, i, c )
    ido = 1
    write ( *, '(2x,i2,2x,11i5)' ) i, c(0:i)
  end do
 
  return
end
subroutine test013 ( )

!*****************************************************************************80
!
!! TEST013 tests COMB_UNRANK.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) cnk
  integer ( kind = 4 ) i4_choose
  integer ( kind = 4 ), parameter :: m = 10
  integer ( kind = 4 ) rank

  cnk = i4_choose ( m, n )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST013'
  write ( *, '(a)' ) '  COMB_UNRANK returns a combination of N things'
  write ( *, '(a)' ) '  out of M, given the lexicographic rank.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The total set size is M = ', m
  write ( *, '(a,i8)' ) '  The subset size is N =    ', n
  write ( *, '(a,i8)' ) '  The number of combinations of N out of M is ', cnk
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   Rank	  Combination'
  write ( *, '(a)' ) ' '
 
  do rank = 1, 3
    call comb_unrank ( m, n, rank, a )
    write ( *, '(2x,i3,3x,5i4)' ) rank, a(1:n)
  end do
 
  do rank = 6, 8
    call comb_unrank ( m, n, rank, a )
    write ( *, '(2x,i3,3x,5i4)' ) rank, a(1:n)
  end do
 
  do rank = 250, 252
    call comb_unrank ( m, n, rank, a )
    write ( *, '(2x,i3,3x,5i4)' ) rank, a(1:n)
  end do
 
  return
end
subroutine test014 ( )

!*****************************************************************************80
!
!! TEST014 tests R8_CHOOSE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) cnk
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n
  real ( kind = 8 ) r8_choose

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST014'
  write ( *, '(a)' ) '  R8_CHOOSE evaluates C(N,K).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       N         K    CNK'
  write ( *, '(a)' ) ' '

  do n = 0, 4
    do k = 0, n
      cnk = r8_choose ( n, k )
      write ( *, '(2x,i8,2x,i8,g14.6)' ) n, k, cnk
    end do
  end do

  return
end
subroutine test015 ( )

!*****************************************************************************80
!
!! TEST015 tests I4_CHOOSE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) cnk
  integer ( kind = 4 ) i4_choose
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST015'
  write ( *, '(a)' ) '  I4_CHOOSE evaluates C(N,K).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       N       K      CNK'
  write ( *, '(a)' ) ' '

  do n = 0, 4
    do k = 0, n
      cnk = i4_choose ( n, k )
      write ( *, '(2x,i8,i8,i8)' ) n, k, cnk
    end do
  end do

  return
end
subroutine test016 ( )

!*****************************************************************************80
!
!! TEST016 tests COMP_NEXT and COMP_ENUM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 July 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: k = 3

  integer ( kind = 4 ) a(k)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ih
  integer ( kind = 4 ) it
  logical more
  integer ( kind = 4 ) n
  integer ( kind = 4 ) number

  n = 6
  more = .false.

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST016'
  write ( *, '(a)' ) '  COMP_NEXT generates compositions.'
  write ( *, '(a)' ) '  COMP_ENUM enumerates compositions.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Seeking all compositions of N = ', n
  write ( *, '(a,i8,a)' ) '  using ', k, ' parts.'

  call comp_enum ( n, k, number )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of compositions is ', number
  write ( *, '(a)' ) ' '

  i = 0

  do

    call comp_next ( n, k, a, more, ih, it )

    i = i + 1
    write ( *, '(2x,i4,2x,8i4)' ) i, a(1:k)

    if ( .not. more )  then
      exit
    end if

  end do
 
  return
end
subroutine test017 ( )

!*****************************************************************************80
!
!! TEST017 tests COMP_RANDOM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: k = 5

  integer ( kind = 4 ) a(k)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) seed
  integer ( kind = 4 ), parameter :: test_num = 5

  n = 10
  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST017'
  write ( *, '(a)' ) '  COMP_RANDOM generates random compositions.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Seeking random compositions of N = ', n
  write ( *, '(a,i8,a)' ) '  using ', k, ' parts.'
  write ( *, '(a)' ) ' '

  do i = 1, test_num
    call comp_random ( n, k, seed, a )
    write ( *, '(2x,8i4)' ) a(1:k)
  end do
 
  return
end
subroutine test0174 ( )

!*****************************************************************************80
!
!! TEST0174 tests COMPNZ_NEXT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    16 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: k = 3

  integer ( kind = 4 ) a(k)
  integer ( kind = 4 ) i
  logical more
  integer ( kind = 4 ) n
  integer ( kind = 4 ) number

  more = .false.
  n = 6

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0174'
  write ( *, '(a)' ) '  COMPNZ_NEXT generates compositions with nonzero parts.'
  write ( *, '(a)' ) '  COMPNZ_ENUM enumerates them.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Seeking all compositions of N = ', n
  write ( *, '(a,i8,a)' ) '  using ', k, ' nonzero parts.'

  call compnz_enum ( n, k, number )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of these compositions is ', number
  write ( *, '(a)' ) ' '

  i = 0

  do

    call compnz_next ( n, k, a, more )

    i = i + 1

    write ( *, '(2x,i4,2x,8i4)' ) i, a(1:k)

    if ( .not. more )  then
      exit
    end if

  end do
 
  return
end
subroutine test0175 ( )

!*****************************************************************************80
!
!! TEST0175 tests COMPNZ_RANDOM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 December 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: k = 5

  integer ( kind = 4 ) a(k)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) seed
  integer ( kind = 4 ), parameter :: test_num = 5

  n = 10
  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0175'
  write ( *, '(a)' ) '  COMPNZ_RANDOM generates random compositions'
  write ( *, '(a)' ) '  with nonzero parts.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Seeking random compositions of N = ', n
  write ( *, '(a,i8,a)' ) '  using ', k, ' nonzero parts.'
  write ( *, '(a)' ) ' '

  do i = 1, test_num
    call compnz_random ( n, k, seed, a )
    write ( *, '(2x,8i4)' ) a(1:k)
  end do
 
  return
end
subroutine test018 ( )

!*****************************************************************************80
!
!! TEST018 tests CONGRUENCE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 21

  integer ( kind = 4 ) a
  integer ( kind = 4 ), save, dimension ( test_num ) :: a_test = (/ &
     1027,  1027,  1027,  1027, -1027, &
    -1027, -1027, -1027,     6,     0, &
        0,     0,     1,     1,     1, &
     1024,     0,     0,     5,     2, &
        7 /)
  integer ( kind = 4 ) b
  integer ( kind = 4 ), save, dimension ( test_num ) :: b_test = (/ &
      712,   712,  -712,  -712,   712, &
      712,  -712,  -712,     8,     0, &
        1,     1,     0,     0,     1, &
   -15625,     0,     3,     0,     4, &
       19 /)
  integer ( kind = 4 ) c
  integer ( kind = 4 ), save, dimension ( test_num ) :: c_test = (/ &
        7,    -7,     7,    -7,     7, &
       -7,     7,    -7,    50,     0, &
        0,     1,     0,     1,     0, &
    11529,     1,    11,    19,     7, &
        1 /)
  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) result
  integer ( kind = 4 ) test
  integer ( kind = 4 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST018'
  write ( *, '(a)' ) '  CONGRUENCE solves a congruence equation:'
  write ( *, '(a)' ) '    A * X = C mod ( B )'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '   I        A         B         C         X     Mod ( A*X-C,B)'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    a = a_test(test)
    b = b_test(test)
    c = c_test(test)

    call congruence ( a, b, c, ierror, x )

    if ( ierror /= 0 ) then
      write ( *, '(2x,i2,2x,3i10,a,i10)' ) test, a, b, c, &
        ' Error code = ', ierror
    else
      if ( b /= 0 ) then
        result = i4_modp ( a * x - c, b )
      else
        result = 0
      end if
      write ( *, '(2x,i2,2x,5i10)' ) test, a, b, c, x, result
    end if

  end do

  return
end
subroutine test019 ( )

!*****************************************************************************80
!
!! TEST019 tests COUNT_POSE_RANDOM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) blocks(6)
  integer ( kind = 4 ) goal
  integer ( kind = 4 ) i
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST019'
  write ( *, '(a)' ) '  COUNT_POSE_RANDOM poses a random problem for '
  write ( *, '(a)' ) '  the game "The Count is Good".'

  seed = 123456789

  do i = 1, 5

    call count_pose_random ( seed, blocks, goal )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Problem #', i
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '    The goal number:'
    write ( *, '(a)' ) ' '
    write ( *, '(6x,i8)' ) goal
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '    The available numbers are '
    write ( *, '(a)' ) ' '
    write ( *, '(6x,6i4)' ) blocks(1:6)

  end do

  return
end
subroutine test020 ( )

!*****************************************************************************80
!
!! TEST020 tests R8_TO_CFRAC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  integer ( kind = 4 ) a(0:n)
  real ( kind = 8 ) error
  integer ( kind = 4 ) i
  integer ( kind = 4 ) p(-1:n)
  integer ( kind = 4 ) q(-1:n)
  real ( kind = 8 ) r
  real ( kind = 8 ) r8_pi
  real ( kind = 8 ) temp

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST020'
  write ( *, '(a)' ) '  R8_TO_CFRAC converts a real number to a'
  write ( *, '(a)' ) '  a sequence of continued fraction convergents.'

  r = 2.0D+00 * r8_pi ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Use the real number R = ', r

  call r8_to_cfrac ( r, n, a, p, q )

  write ( *, '(a)' ) ' '

  do i = 0, n
    temp = real ( p(i), kind = 8 ) / real ( q(i), kind = 8 )
    error = r - temp
    write ( *, '(2x,3i12,2g14.6)' ) a(i), p(i), q(i), temp, error
  end do

  return
end
subroutine test021 ( )

!*****************************************************************************80
!
!! TEST021 tests DEBRUIJN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 3

  integer ( kind = 4 ) m
  integer ( kind = 4 ), parameter, dimension ( test_num ) :: m_test = (/ &
    2, 3, 2 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter, dimension ( test_num ) :: n_test = (/ &
    3, 3, 4 /)
  integer ( kind = 4 ), dimension (27) :: string
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST021'
  write ( *, '(a)' ) '  DEBRUIJN computes a de Bruijn string.'

  do test = 1, test_num

    m = m_test(test)
    n = n_test(test)

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  The alphabet size is M = ', m
    write ( *, '(a,i8)' ) '  The string length is N = ', n

    call debruijn ( m, n, string )

    write ( *, '(a)' ) ' '
    write ( *, '(4x,80i1)' ) string(1:m**n)

  end do

  return
end
subroutine test022 ( )

!*****************************************************************************80
!
!! TEST022 tests DEC_ADD.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) abot
  integer ( kind = 4 ) atop
  integer ( kind = 4 ) bbot
  integer ( kind = 4 ) btop
  integer ( kind = 4 ) cbot
  integer ( kind = 4 ) ctop
  integer ( kind = 4 ) dec_digit
  integer ( kind = 4 ) ierror
  character ( len = 15 ) string

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST022'
  write ( *, '(a)' ) '  DEC_ADD adds two decimals.'

  dec_digit = 3

  atop = 128
  abot = - 1
  btop = 438
  bbot = - 2

  call dec_add ( atop, abot, btop, bbot, dec_digit, ctop, cbot, ierror )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of decimal places is ', dec_digit
  write ( *, '(a)' ) ' '
  call dec_to_s ( atop, abot, string )
  write ( *, '(a)' ) '  A = ' // trim ( string )
  call dec_to_s ( btop, bbot, string )
  write ( *, '(a)' ) '  B = ' // trim ( string )
  call dec_to_s ( ctop, cbot, string )
  write ( *, '(a)' ) '  C = ' // trim ( string )
 
  return
end
subroutine test023 ( )

!*****************************************************************************80
!
!! TEST023 tests DEC_DIV.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) abot
  integer ( kind = 4 ) atop
  integer ( kind = 4 ) bbot
  integer ( kind = 4 ) btop
  integer ( kind = 4 ) cbot
  integer ( kind = 4 ) ctop
  integer ( kind = 4 ) dec_digit
  integer ( kind = 4 ) ierror
  character ( len = 15 ) string

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST023'
  write ( *, '(a)' ) '  DEC_DIV divides two decimals.'

  dec_digit = 3

  atop = 523
  abot = -1
  btop = 134
  bbot = 2

  call dec_div ( atop, abot, btop, bbot, dec_digit, ctop, cbot, ierror )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of decimal places is ', dec_digit
  write ( *, '(a)' ) ' '
  call dec_to_s ( atop, abot, string )
  write ( *, '(a)' ) '  A = ' // trim ( string )
  call dec_to_s ( btop, bbot, string )
  write ( *, '(a)' ) '  B = ' // trim ( string )
  call dec_to_s ( ctop, cbot, string )
  write ( *, '(a)' ) '  C = ' // trim ( string )
 
  return
end
subroutine test024 ( )

!*****************************************************************************80
!
!! TEST024 tests DEC_MUL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) abot
  integer ( kind = 4 ) atop
  integer ( kind = 4 ) bbot
  integer ( kind = 4 ) btop
  integer ( kind = 4 ) cbot
  integer ( kind = 4 ) ctop
  integer ( kind = 4 ) dec_digit
  character ( len = 15 ) string

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST024'
  write ( *, '(a)' ) '  DEC_MUL multiplies two decimals.'

  dec_digit = 2

  atop = 14
  abot = - 4
  btop = 16
  bbot = 2

  call dec_mul ( atop, abot, btop, bbot, dec_digit, ctop, cbot )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of decimal places is ', dec_digit
  write ( *, '(a)' ) ' '
  call dec_to_s ( atop, abot, string )
  write ( *, '(a)' ) '  A = ' // trim ( string )
  call dec_to_s ( btop, bbot, string )
  write ( *, '(a)' ) '  B = ' // trim ( string )
  call dec_to_s ( ctop, cbot, string )
  write ( *, '(a)' ) '  C = ' // trim ( string )
 
  return
end
subroutine test025 ( )

!*****************************************************************************80
!
!! TEST025 tests DEC_ROUND.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 7

  integer ( kind = 4 ), dimension ( test_num ) :: d_test = (/ & 
      1, 2, 3, 4, 2, 3, 4 /)
  integer ( kind = 4 ) dec_digit
  integer ( kind = 4 ) exponent
  integer ( kind = 4 ), dimension ( test_num ) :: exponent_test = (/ &
     -1,  -1, -1, -1, 2, 2, 2 /)
  integer ( kind = 4 ) mantissa
  integer ( kind = 4 ), dimension ( test_num ) :: mantissa_test = (/ &
    523, 523, 523, 523, 6340, 6340, 6340 /)
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST025'
  write ( *, '(a)' ) '  DEC_ROUND "rounds" a decimal to a number of digits.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '           -----Before-------  -----After--------'
  write ( *, '(a)' ) '  Digits   Mantissa  Exponent  Mantissa  Exponent'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    dec_digit = d_test(test)

    mantissa = mantissa_test(test)
    exponent = exponent_test(test)

    call dec_round ( mantissa, exponent, dec_digit, mantissa, exponent )

    write ( *, '(2x,3i8,4x,2i8)' ) &
      dec_digit, mantissa_test(test), exponent_test(test), mantissa, exponent

  end do

  return
end
subroutine test026 ( )

!*****************************************************************************80
!
!! TEST026 tests DEC_TO_RAT and RAT_TO_DEC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) exponent
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) mantissa
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) r3
  integer ( kind = 4 ) rat_bot
  integer ( kind = 4 ) rat_bot2
  integer ( kind = 4 ) rat_top
  integer ( kind = 4 ) rat_top2
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST026'
  write ( *, '(a)' ) '  RAT_TO_DEC fraction => decimal,'
  write ( *, '(a)' ) '  DEC_TO_RAT decimal => fraction.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this test, choose the top and bottom'
  write ( *, '(a)' ) '  of a rational at random, and compute the'
  write ( *, '(a)' ) '  equivalent real number.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Then convert to decimal, and the equivalent real.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Then convert back to rational and the equivalent real.'
  
  seed = 123456789

  do i = 1, 10

    rat_top = i4_uniform ( -1000, 1000, seed )

    rat_bot = i4_uniform (     1, 1000, seed )

    r1 = real ( rat_top, kind = 8 ) / real ( rat_bot, kind = 8 )

    call rat_to_dec ( rat_top, rat_bot, mantissa, exponent )

    r2 = real ( mantissa, kind = 8 ) * 10.0D+00**( exponent )
 
    call dec_to_rat ( mantissa, exponent, rat_top2, rat_bot2 )
    r3 = real ( rat_top2, kind = 8 ) / real ( rat_bot2, kind = 8 )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,f10.6,a,i12,a,i12)' ) r1, '=', rat_top, '/', rat_bot
    write ( *, '(2x,f10.6,a,i12,a,i12)' ) r2, '=', mantissa, '*10^', exponent
    write ( *, '(2x,f10.6,a,i12,a,i12)' ) r3, '=', rat_top2, '/', rat_bot2

  end do
 
  return
end
subroutine test027 ( )

!*****************************************************************************80
!
!! TEST027 tests DEC_TO_S.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) exponent
  integer ( kind = 4 ) mantissa
  character ( len = 100 ) s

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST027'
  write ( *, '(a)' ) '  DEC_TO_S prints a decimal value.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Mantissa  Exponent  String'
  write ( *, '(a)' ) ' '

  mantissa = 523
  exponent = -1
  call dec_to_s ( mantissa, exponent, s )
  write ( *, '(2x,i8,2x,i8,2x,a)' ) mantissa, exponent, trim ( s )

  mantissa = 134
  exponent = 2
  call dec_to_s ( mantissa, exponent, s )
  write ( *, '(2x,i8,2x,i8,2x,a)' ) mantissa, exponent, trim ( s )

  mantissa = -134
  exponent = 2
  call dec_to_s ( mantissa, exponent, s )
  write ( *, '(2x,i8,2x,i8,2x,a)' ) mantissa, exponent, trim ( s )

  mantissa = 0
  exponent = 10
  call dec_to_s ( mantissa, exponent, s )
  write ( *, '(2x,i8,2x,i8,2x,a)' ) mantissa, exponent, trim ( s )

  do exponent = -8, 3
    mantissa = 123456
    call dec_to_s ( mantissa, exponent, s )
    write ( *, '(2x,i8,2x,i8,2x,a)' ) mantissa, exponent, trim ( s )
  end do

  return
end
subroutine test028 ( )

!*****************************************************************************80
!
!! TEST028 tests DEC_WIDTH.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) dec_width
  integer ( kind = 4 ) exponent
  integer ( kind = 4 ) i
  integer ( kind = 4 ) mantissa

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST028'
  write ( *, '(a)' ) '  DEC_WIDTH determines the "width" of a decimal.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Mantissa  Exponent  Width'
  write ( *, '(a)' ) ' '

  mantissa = 523
  exponent = -1
  i = dec_width ( mantissa, exponent )
  write ( *, '(2x,i8,2x,i8,2x,i8)' ) mantissa, exponent, i

  mantissa = 134
  exponent = 2
  i = dec_width ( mantissa, exponent )
  write ( *, '(2x,i8,2x,i8,2x,i8)' ) mantissa, exponent, i

  mantissa = -134
  exponent = 2
  i = dec_width ( mantissa, exponent )
  write ( *, '(2x,i8,2x,i8,2x,i8)' ) mantissa, exponent, i

  mantissa = 0
  exponent = 10
  i = dec_width ( mantissa, exponent )
  write ( *, '(2x,i8,2x,i8,2x,i8)' ) mantissa, exponent, i

  do exponent = -8, 3
    mantissa = 123456
    i = dec_width ( mantissa, exponent )
    write ( *, '(2x,i8,2x,i8,2x,i8)' ) mantissa, exponent, i
  end do

  return
end
subroutine test029 ( )

!*****************************************************************************80
!
!! TEST029 tests DECMAT_DET.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n3 = 3
  integer ( kind = 4 ), parameter :: n4 = 4

  integer ( kind = 4 ) dec_digit
  integer ( kind = 4 ) i
  integer ( kind = 4 ) a3(n3,n3)
  integer ( kind = 4 ) a4(n4,n4)
  integer ( kind = 4 ) b3(n3,n3)
  integer ( kind = 4 ) b4(n4,n4)
  integer ( kind = 4 ) dbot
  integer ( kind = 4 ) dtop
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) r

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST029'
  write ( *, '(a)' ) '  DECMAT_DET: determinant of a decimal matrix.'
  write ( *, '(a)' ) ' '
 
  dec_digit = 7

  k = 0
  do i = 1, n3
    do j = 1, n3
      k = k + 1
      a3(i,j) = k
    end do
  end do

  b3(1:n3,1:n3) = 0
 
  call decmat_print ( n3, n3, a3, b3, '  The 123/456/789 matrix:' )

  call decmat_det ( n3, a3, b3, dec_digit, dtop, dbot, ierror )
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Determinant of the 123/456/789 matrix'
  write ( *, '(2x,i8,a,i8)' ) dtop, ' * 10** ', dbot
 
  do i = 1, n4
    do j = 1, n4
      r = 1.0D+00 / real ( i + j, kind = 8 )
      call r8_to_dec ( r, dec_digit, a4(i,j), b4(i,j) )
    end do
  end do
 
  call decmat_print ( n4, n4, a4, b4, '  The Hilbert matrix:' )

  call decmat_det ( n4, a4, b4, dec_digit, dtop, dbot, ierror )
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Determinant of the Hilbert matrix:'
  write ( *, '(2x,i8,a,i8)' ) dtop, ' * 10 ** ', dbot 

  do i = 1, n3
    do j = 1, n3
      if ( i == j ) then
        a3(i,j) = 2
      else if ( i == j+1 .or. i == j-1 ) then
        a3(i,j) = -1
      else
        a3(i,j) = 0
      end if
      b3(i,j) = 0
    end do
  end do
 
  call decmat_print ( n3, n3, a3, b3, '  The -1,2,-1 matrix:' )

  call decmat_det ( n3, a3, b3, dec_digit, dtop, dbot, ierror )
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Determinant of the -1,2,-1 matrix:'
  write ( *, '(2x,i8,a,i8)' ) dtop, ' * 10 ** ', dbot
 
  return
end
subroutine test021a ( )

!*****************************************************************************80
!
!! TEST021a tests DERANGE_ENUM, DERANGE_ENUM2 and DERANGE_ENUM3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  integer ( kind = 4 ) d(0:n)
  integer ( kind = 4 ) derange_enum
  integer ( kind = 4 ) derange_enum3
  integer ( kind = 4 ) i

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST021a'
  write ( *, '(a)' ) '  DERANGE_ENUM counts derangements;'
  write ( *, '(a)' ) '  DERANGE_ENUM2 counts derangements.'
  write ( *, '(a)' ) '  DERANGE_ENUM3 counts derangements.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  N    # of derangements'
  write ( *, '(a)' ) ' '

  call derange_enum2 ( n, d )

  do i = 0, n
    write ( *, '(2x,i8,2x,3i10)' ) i, derange_enum(i), d(i), derange_enum3(i)
  end do

  return
end
subroutine test022a ( )

!*****************************************************************************80
!
!! TEST022a tests DERANGE_BACK_NEXT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  integer ( kind = 4 ) a(n)
  logical more
  integer ( kind = 4 ) number

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST022a'
  write ( *, '(a)' ) '  DERANGE_BACK_NEXT generates derangements'
  write ( *, '(a)' ) '  using backtracking.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Here, we seek all derangments of order N = ', n
  write ( *, '(a)' ) ' '

  more = .false.
  number = 0

  do

    call derange_back_next ( n, a, more )

    if ( .not. more ) then
      exit
    end if

    number = number + 1
    write ( *, '(2x,i4,4x,8i4)' ) number, a(1:n)

  end do

  return
end
subroutine test023a ( )

!*****************************************************************************80
!
!! TEST023a tests DERANGE_WEED_NEXT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  integer ( kind = 4 ) a(n)
  logical more
  integer ( kind = 4 ) number

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST023a'
  write ( *, '(a)' ) '  DERANGE_WEED_NEXT generates derangements'
  write ( *, '(a)' ) '  by generating ALL permutations, and "weeding out"'
  write ( *, '(a)' ) '  the ones that are not derangements.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Here, we seek all derangments of order N = ', n
  write ( *, '(a)' ) ' '

  more = .false.
  number = 0
 
  do

    call derange_weed_next ( n, a, more )

    number = number + 1
    write ( *, '(2x,i4,4x,8i4)' ) number, a(1:n)

    if ( .not. more ) then
      exit
    end if
 
  end do

  return
end
subroutine test024a ( )

!*****************************************************************************80
!
!! TEST024a calls DIGRAPH_ARC_EULER.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nedge = 7
  integer ( kind = 4 ), parameter :: nnode = 5

  integer ( kind = 4 ) i
  integer ( kind = 4 ) in
  integer ( kind = 4 ), dimension ( nedge ) :: inode = (/ 2, 1, 2, 1, 3, 5, 4 /)
  integer ( kind = 4 ) j
  integer ( kind = 4 ), dimension ( nedge ) :: jnode = (/ 5, 4, 3, 2, 1, 1, 2 /)
  integer ( kind = 4 ) jp1
  logical success
  integer ( kind = 4 ) trail(nedge)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST024a'
  write ( *, '(a)' ) '  DIGRAPH_ARC_EULER finds an Euler circuit of a digraph.'

  call digraph_arc_print ( nedge, inode, jnode, &
    '  The arc list of the digraph:' )

  call digraph_arc_euler ( nnode, nedge, inode, jnode, success, trail )

  if ( success ) then

    call i4vec_print ( nedge, trail, '  The edge list of the Euler circuit:' )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The node list of the Euler circuit:'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '    I  Edge  Node'
    write ( *, '(a)' ) ' '

    do i = 1, nedge

      j = trail(i)

      if ( i == nedge ) then
        jp1 = trail(1)
      else
        jp1 = trail(i+1)
      end if

      if ( jnode(j) == inode(jp1) ) then
        in = jnode(j)
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'The circuit has failed!'
        exit
      end if

      write ( *, '(2x,3i8)' ) i, j, in

    end do

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The digraph is not eulerian.'
    write ( *, '(a)' ) ' '

  end if

  return
end
subroutine test025a ( )

!*****************************************************************************80
!
!! TEST025a tests DIOPHANTINE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 20

  integer ( kind = 4 ) a
  integer ( kind = 4 ), dimension ( test_num ) :: a_test = (/ &
     1027,  1027,  1027, 1027, -1027, &
    -1027, -1027, -1027,    6,     0, &
        0,     0,     1,    1,     1, &
     1024,     0,     0,    5,     2 /)
  integer ( kind = 4 ) b
  integer ( kind = 4 ), dimension ( test_num) ::  b_test = (/ &
       712,   712, -712, -712, 712, &
       712,  -712, -712,    8,   0, &
         1,     1,    0,    0,   1, &
    -15625,     0,    3,    0,   4 /)
  integer ( kind = 4 ) c
  integer ( kind = 4 ), dimension ( test_num) ::  c_test = (/ &
         7,    -7,    7,   -7,   7, &
        -7,     7,   -7,   50,   0, &
         0,     1,    0,    1,   0, &
     11529,     1,   11,   19,   7 /)
  integer ( kind = 4 ) error
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) test
  integer ( kind = 4 ) x
  integer ( kind = 4 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST025a'
  write ( *, '(a)' ) '  DIOPHANTINE solves a Diophantine equation:'
  write ( *, '(a)' ) '    A * X + B * Y = C'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '        A         B         C         X     Y     Error'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    a = a_test(test)
    b = b_test(test)
    c = c_test(test)

    call diophantine ( a, b, c, ierror, x, y )

    if ( ierror /= 0 ) then
      write ( *, '(2x,3i10,a,i10)' ) a, b, c, ' Error code = ', ierror
    else
      error = a * x + b * y - c
      write ( *, '(2x,6i10)' ) a, b, c, x, y, error
    end if

  end do

  return
end
subroutine test026a ( )

!*****************************************************************************80
!
!! TEST026a tests DIOPHANTINE_SOLUTION_MINIMIZE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) c
  integer ( kind = 4 ) r
  integer ( kind = 4 ) x
  integer ( kind = 4 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST026a'
  write ( *, '(a)' ) '  DIOPHANTINE_SOLUTION_MINIMIZE computes a minimal'
  write ( *, '(a)' ) '  Euclidean norm solution of a Diophantine equation:'
  write ( *, '(a)' ) '    A * X + B * Y = C'

  a = 4096
  b = - 15625
  c = 46116
  x = 665499996
  y = 174456828

  r = a * x + b * y - c

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Coefficients:'
  write ( *, '(a,i12)' ) '    A = ', a
  write ( *, '(a,i12)' ) '    B = ', b
  write ( *, '(a,i12)' ) '    C = ', c
  write ( *, '(a)' ) '  Solution:'
  write ( *, '(a,i12)' ) '    X = ', x
  write ( *, '(a,i12)' ) '    Y = ', y
  write ( *, '(a)' ) '  Residual R = A * X + B * Y - C:'
  write ( *, '(a,i12)' ) '    R = ', r

  call diophantine_solution_minimize ( a, b, x, y )

  r = a * x + b * y - c

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The minimized solution:'
  write ( *, '(a,i12)' ) '    X = ', x
  write ( *, '(a,i12)' ) '    Y = ', y
  write ( *, '(a)' ) '  Residual R = A * X + B * Y - C:'
  write ( *, '(a,i12)' ) '    R = ', r

  x = 15621
  y = 4092

  r = a * x + b * y - c

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The minimal positive solution:'
  write ( *, '(a,i12)' ) '    X = ', x
  write ( *, '(a,i12)' ) '    Y = ', y
  write ( *, '(a)' ) '  Residual R = A * X + B * Y - C:'
  write ( *, '(a,i12)' ) '    R = ', r

  return
end
subroutine test026b ( )

!*****************************************************************************80
!
!! TEST026B tests DVEC_ADD and DVEC_SUB;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 November 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  integer ( kind = 4 ) dvec1(n)
  integer ( kind = 4 ) dvec2(n)
  integer ( kind = 4 ) dvec3(n)
  integer ( kind = 4 ) dvec4(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) :: seed = 123456789
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 10

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST026B'
  write ( *, '(a)' ) '  DVEC_ADD adds decimal vectors representing integers;'
  write ( *, '(a)' ) '  DVEC_SUB subtracts decimal vectors representing integers;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '        I        J        K = I + J    L = I - J'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    
    i = i4_uniform ( -100, 100, seed )
    j = i4_uniform ( -100, 100, seed )

    write ( *, '(a)' ) ' '

    write ( *, '(2x,i8,2x,i8)' ) i, j

    k = i + j
    l = i - j

    write ( *, '(a20,2x,i8,2x,i8)' ) '  Directly:         ', k, l

    call i4_to_dvec ( i, n, dvec1 )
    call i4_to_dvec ( j, n, dvec2 )

    call dvec_add ( n, dvec1, dvec2, dvec3 )
    call dvec_to_i4 ( n, dvec3, k )

    call dvec_sub ( n, dvec1, dvec2, dvec4 )
    call dvec_to_i4 ( n, dvec4, l )

    write ( *, '(a20,2x,i8,2x,i8)' ) '  DVEC_ADD, DVEC_SUB', k, l

  end do

  return
end
subroutine test026c ( )

!*****************************************************************************80
!
!! TEST026C tests DVEC_COMPLEMENTX;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 November 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  integer ( kind = 4 ) dvec1(n)
  integer ( kind = 4 ) dvec2(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) :: seed = 123456789
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 5

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST026C'
  write ( *, '(a)' ) '  DVEC_COMPLEMENTX returns the ten''s complement'
  write ( *, '(a)' ) '  of a (signed) decimal vector;'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    
    i = i4_uniform ( -100, 100, seed )

    call i4_to_dvec ( i, n, dvec1 )

    call dvec_complementx ( n, dvec1, dvec2 )

    call dvec_to_i4 ( n, dvec2, j )

    write ( *, '(a)' ) ' '
    write ( *, '(a,2x,i8)' ) '  I = ', i
    write ( *, '(a,2x,i8)' ) '  J = ', j
    call dvec_print ( n, dvec1, ' ' )
    call dvec_print ( n, dvec2, ' ' )

  end do

  return
end
subroutine test026d ( )

!*****************************************************************************80
!
!! TEST026D tests DVEC_MUL;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 November 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  integer ( kind = 4 ) dvec1(n)
  integer ( kind = 4 ) dvec2(n)
  integer ( kind = 4 ) dvec3(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) :: seed = 123456789
  integer ( kind = 4 ) test
  integer ( kind = 4 ) test2
  integer ( kind = 4 ), parameter :: test_num = 10
  integer ( kind = 4 ), parameter :: test2_num = 2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST026D'
  write ( *, '(a)' ) '  DVEC_MUL multiplies decimal vectors '
  write ( *, '(a)' ) '  representing integers;'

  do test2 = 1, test2_num

    if ( test2 == 1 ) then

      n2 = n

    else if ( test2 == 2 ) then

      n2 = 6

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  NOW REPEAT THE TEST...'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  but use too few digits to represent big products.'
      write ( *, '(a)' ) '  This corresponds to an "overflow".'
      write ( *, '(a)' ) '  The result here should get the final decimal'
      write ( *, '(a)' ) '  digits correctly, though.'
    end if

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '        I        J        K = I * J'
    write ( *, '(a)' ) ' '

    do test = 1, test_num
      
      i = i4_uniform ( -1000, 1000, seed )
      j = i4_uniform ( -1000, 1000, seed )

      write ( *, '(a)' ) ' '

      write ( *, '(2x,i8,2x,i8)' ) i, j

      k = i * j

      write ( *, '(a20,2x,i8)' ) '  Directly:         ', k

      call i4_to_dvec ( i, n2, dvec1 )
      call i4_to_dvec ( j, n2, dvec2 )

      call dvec_mul ( n2, dvec1, dvec2, dvec3 )
      call dvec_to_i4 ( n2, dvec3, k )

      write ( *, '(a20,2x,i8)' ) '  DVEC_MUL          ', k

    end do

  end do

  return
end
subroutine test027a ( )

!*****************************************************************************80
!
!! TEST027a tests EQUIV_NEXT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) jarray(n)
  logical more
  integer ( kind = 4 ) nc
  integer ( kind = 4 ) rank

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST027a'
  write ( *, '(a)' ) '  EQUIV_NEXT generates all partitions of a set.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,6i4)' ) '  Rank/element:', ( i, i = 1, n )
  write ( *, '(a)' ) ' '
 
  rank = 0
  more = .false.
 
  do
 
    call equiv_next ( n, nc, jarray, a, more )
 
    rank = rank + 1
    write ( *, '(2x,15i4)' ) rank, a(1:n)
 
    if ( .not. more ) then
      exit
    end if

  end do

  return
end
subroutine test028a ( )

!*****************************************************************************80
!
!! TEST028a tests EQUIV_NEXT2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4

  integer ( kind = 4 ) a(n)
  logical done
  integer ( kind = 4 ) i
  integer ( kind = 4 ) rank

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST028a'
  write ( *, '(a)' ) '  EQUIV_NEXT2 generates all partitions of a set.'
  write ( *, '(a,i8)' ) '  Here, N = ', n
  rank = 0
  done = .true.
  write ( *, '(a)' ) ' '
  write ( *, '(a,6i4)' ) '  Rank/element:', ( i, i = 1, n )
  write ( *, '(a)' ) ' '
 
  do

    call equiv_next2 ( n, a, done )

    if ( done ) then
      exit
    end if

    rank = rank + 1
    write ( *, '(2x,i4,10x,6i4)' ) rank, a(1:n)

  end do

  return 
end
subroutine test029a ( )

!*****************************************************************************80
!
!! TEST029a tests EQUIV_RANDOM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 May 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4

  integer ( kind = 4 ) a(n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) npart
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST029a'
  write ( *, '(a)' ) '  EQUIV_RANDOM selects a random set partition.'
 
  seed = 123456789

  do i = 1, 5
 
    call equiv_random ( n, seed, npart, a, b )

    call equiv_print ( n, a, '  The partition:' )
 
  end do
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Repeat, but print with EQUIV_PRINT2.'
  write ( *, '(a)' ) ' '
 
  seed = 123456789

  do i = 1, 5
 
    call equiv_random ( n, seed, npart, a, b )

    call equiv_print2 ( n, a, '  The partition:' )
 
  end do
 
  return
end
subroutine test0295 ( )

!*****************************************************************************80
!
!! TEST0295 tests EULER.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: nmax = 9

  integer ( kind = 4 ) ieuler(0:nmax)
  integer ( kind = 4 ) n

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0295'
  write ( *, '(a)' ) '  EULER gets rows of Euler''s triangle.'
  write ( *, '(a)' ) ' '

  do n = 0, nmax
    call euler ( n, ieuler )
    write ( *, '(2x,10i7)' ) ieuler(0:n)
  end do
 
  return
end
subroutine test0304 ( )

!*****************************************************************************80
!
!! TEST0304 tests FROBENIUS_NUMBER_ORDER2 and FROBENIUS_NUMBER_ORDER2_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    20 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) c1
  integer ( kind = 4 ) c2
  integer ( kind = 4 ) f1
  integer ( kind = 4 ) f2
  integer ( kind = 4 ) frobenius_number_order2
  integer ( kind = 4 ) n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0304'
  write ( *, '(a)' ) &
    '  FROBENIUS_NUMBER_ORDER2 computes Frobenius numbers of order 2.'
  write ( *, '(a)' ) &
    '  FROBENIUS_NUMBER_ORDER2_VALUES returns some exact values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '        C1        C1   exact F  comput F'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call frobenius_number_order2_values ( n_data, c1, c2, f1 )

    if ( n_data == 0 ) then
      exit
    end if

    f2 = frobenius_number_order2 ( c1, c2 )

    write ( *, '(2x,i8,2x,i8,2x,i8,2x,i8)' ) c1, c2, f1, f2

  end do

  return
end
subroutine test0305 ( )

!*****************************************************************************80
!
!! TEST0305 tests R8_GAMMA_LOG and GAMMA_LOG_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  real ( kind = 8 ) fx2
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ) r8_gamma_log

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0305:'
  write ( *, '(a)' ) '  R8_GAMMA_LOG evaluates the logarithm of the '
  write ( *, '(a)' ) '  Gamma function.'
  write ( *, '(a)' ) '  GAMMA_LOG_VALUES returns some exact values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     X       Exact F       GAMMA_LOG(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call gamma_log_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = r8_gamma_log ( x )

    write ( *, '(2x,f8.4,2g18.10)' ) x, fx, fx2

  end do

  return
end
subroutine test031 ( )

!*****************************************************************************80
!
!! TEST031 tests GRAY_NEXT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) change
  integer ( kind = 4 ) k

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST031'
  write ( *, '(a)' ) '  GRAY_NEXT returns the index of the single item'
  write ( *, '(a)' ) '  to be changed in order to get the next Gray code.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   K  Change  Gray Code'
  write ( *, '(a)' ) ' '

  change = -n
  k = 0

  do

    call gray_next ( n, change )

    if ( change == -n ) then
      exit
    else if ( change == 0 ) then
      a(1:n) = 0
    else
      a(abs(change)) = 1 - a(abs(change))
    end if

    write ( *, '(2x,i2,2x,i8,2x,4i1)' ) k, change, a(1:n)
    k = k + 1
    
  end do

  return
end
subroutine test032 ( )

!*****************************************************************************80
!
!! TEST032 tests GRAY_RANK and GRAY_UNRANK.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) gray
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) rank2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST032'
  write ( *, '(a)' ) '  GRAY_RANK ranks a Gray code;'
  write ( *, '(a)' ) '  GRAY_UNRANK unranks a Gray code.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    R  =                         RANK'
  write ( *, '(a)' ) '    G  =            GRAY_UNRANK2(RANK)'
  write ( *, '(a)' ) '    R2 = GRAY_RANK2(GRAY_UNRANK2(RANK))'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         R         G         R2'
  write ( *, '(a)' ) ' '
 
  do rank = 0, 24
    call gray_unrank ( rank, gray )
    call gray_rank ( gray, rank2 )
    write ( *, '(2x,i8,2x,i8,2x,i8)' ) rank, gray, rank2
  end do

  return
end
subroutine test0321 ( )

!*****************************************************************************80
!
!! TEST0321 tests GRAY_RANK2 and GRAY_UNRANK2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) gray
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) rank2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0321'
  write ( *, '(a)' ) '  GRAY_RANK2 ranks a Gray code;'
  write ( *, '(a)' ) '  GRAY_UNRANK2 unranks a Gray code.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    R  =                         RANK'
  write ( *, '(a)' ) '    G  =            GRAY_UNRANK2(RANK)'
  write ( *, '(a)' ) '    R2 = GRAY_RANK2(GRAY_UNRANK2(RANK))'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         R         G         R2'
  write ( *, '(a)' ) ' '
 
  do rank = 0, 24
    call gray_unrank2 ( rank, gray )
    call gray_rank2 ( gray, rank2 )
    write ( *, '(2x,i8,2x,i8,2x,i8)' ) rank, gray, rank2
  end do

  return
end
subroutine test0322 ( )

!*****************************************************************************80
!
!! TEST0322 tests I4_BCLR.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 2

  integer ( kind = 4 ) i4
  integer ( kind = 4 ), dimension ( test_num ) :: i4_test = (/ &
    101, -31 /)
  integer ( kind = 4 ) i4_bclr
  integer ( kind = 4 ) ivec(0:31)
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) pos
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0322'
  write ( *, '(a)' ) '  I4_BCLR sets a given bit to 0.'
  write ( *, '(a)' ) '  IBCLR is a FORTRAN90 function which does the same.'

  do test = 1, test_num

    i4 = i4_test(test)

    call i4_to_bvec ( i4, 32, ivec )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Working on I4 = ', i4
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '       Pos     Digit       I4_BCLR         IBCLR'
    write ( *, '(a)' ) ' '

    do pos = 0, 31
  
      j1 = i4_bclr ( i4, pos )
      j2 = ibclr ( i4, pos )

      write ( *, '(2x,i8,2x,i8,2x,i12,2x,i12)' ) pos, ivec(pos), j1, j2

    end do

  end do

  return
end
subroutine test03225 ( )

!*****************************************************************************80
!
!! TEST03225 tests I4_BSET.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 2

  integer ( kind = 4 ) i4
  integer ( kind = 4 ), dimension ( test_num ) :: i4_test = (/ &
    101, -31 /)
  integer ( kind = 4 ) i4_bset
  integer ( kind = 4 ) ivec(0:31)
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) pos
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03225'
  write ( *, '(a)' ) '  I4_BSET sets a given bit to 0.'
  write ( *, '(a)' ) '  IBSET is a FORTRAN90 function which does the same.'

  do test = 1, test_num

    i4 = i4_test(test)

    call i4_to_bvec ( i4, 32, ivec )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Working on I4 = ', i4
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '       Pos     Digit       I4_BSET         IBSET'
    write ( *, '(a)' ) ' '

    do pos = 0, 31
  
      j1 = i4_bset ( i4, pos )
      j2 = ibset ( i4, pos )

      write ( *, '(2x,i8,2x,i8,2x,i12,2x,i12)' ) pos, ivec(pos), j1, j2

    end do

  end do

  return
end
subroutine test0323 ( )

!*****************************************************************************80
!
!! TEST0323 tests I4_BTEST.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 2

  integer ( kind = 4 ) i4
  integer ( kind = 4 ), dimension ( test_num ) :: i4_test = (/ &
    101, -31 /)
  logical i4_btest
  integer ( kind = 4 ) ivec(0:31)
  logical j1
  logical j2
  integer ( kind = 4 ) pos
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0323'
  write ( *, '(a)' ) '  I4_BTEST reports whether a given bit is 0 or 1.'
  write ( *, '(a)' ) '  BTEST is a FORTRAN90 function which does the same.'

  do test = 1, test_num

    i4 = i4_test(test)

    call i4_to_bvec ( i4, 32, ivec )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Analyze the integer I4 = ', i4
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '       Pos     Digit  I4_BTEST     BTEST'
    write ( *, '(a)' ) ' '

    do pos = 0, 31
  
      j1 = i4_btest ( i4, pos )
      j2 = btest ( i4, pos )

      write ( *, '(2x,i8,2x,i8,2x,7x,l1,2x,7x,l1)' ) pos, ivec(pos), j1, j2

    end do

  end do

  return
end
subroutine test0324 ( )

!*****************************************************************************80
!
!! TEST0324 tests I4_FACTOR.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter  :: factor_max = 10

  integer ( kind = 4 ) factor(factor_max)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nfactor
  integer ( kind = 4 ) nleft
  integer ( kind = 4 ) power(factor_max)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0324'
  write ( *, '(a)' ) '  I4_FACTOR factors an integer,'

  n = 2**2 * 17 * 37

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The integer is ', n

  call i4_factor ( n, factor_max, nfactor, factor, power, nleft )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Prime representation:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         I  FACTOR(I)  POWER(I)'
  write ( *, '(a)' ) ' '
  if ( abs ( nleft ) /= 1 ) then
    write ( *, '(2x,i8,i8,a)' ) 0, nleft, ' (UNFACTORED PORTION)'
  end if

  do i = 1, nfactor
    write ( *, '(2x,3i8)' ) i, factor(i), power(i)
  end do

  return
end
subroutine test0325 ( )

!*****************************************************************************80
!
!! TEST0325 tests I4_GCD.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_gcd
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0325'
  write ( *, '(a)' ) '  I4_GCD computes the greatest common divisor'
  write ( *, '(a)' ) '  of two integers.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         I       J    I4_GCD(I,J)'
  write ( *, '(a)' ) ' '

  seed = 123456789

  do k = 1, 15
    i = i4_uniform ( -5, 15, seed )
    j = i4_uniform ( 1, 15, seed )
    write ( *, '(2x,3i8)' ) i, j, i4_gcd ( i, j )
  end do

  return
end
subroutine test0327 ( )

!*****************************************************************************80
!
!! TEST0327 tests I4_LOG_10.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 21

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_log_10
  integer ( kind = 4 ), dimension ( n ) :: x = (/ &
      0,    1,    2,  3,  9, 10, 11,  99, 100, 101, &
    999, 1000, 1001, -1, -2, -3, -9, -10, -11, -99, &
   -101 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0327'
  write ( *, '(a)' ) '  I4_LOG_10: whole part of log base 10,'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         X        I4_LOG_10'
  write ( *, '(a)' ) ' '

  do i = 1, n

    write ( *, '(2x,i8,i12)' ) x(i), i4_log_10 ( x(i) )

  end do

  return
end
subroutine test058 ( )

!*****************************************************************************80
!
!! TEST058 tests I4_PARTITION_CONJ.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 14
  integer ( kind = 4 ), parameter :: npart1 = 4

  integer ( kind = 4 ) a1(npart1)
  integer ( kind = 4 ) a2(n)
  integer ( kind = 4 ) mult1(npart1)
  integer ( kind = 4 ) mult2(n)
  integer ( kind = 4 ) npart2

  a1(1:npart1) = (/ 2, 5, 1, 4 /)
  mult1(1:npart1) = (/ 1, 1, 3, 1 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST058'
  write ( *, '(a)' ) '  I4_PARTITION_CONJ conjugates an integer partition.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Original partition:'
  write ( *, '(a)' ) ' '

  call i4_partition_print ( n, npart1, a1, mult1 )

  call i4_partition_conj ( n, a1, mult1, npart1, a2, mult2, npart2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Conjugate partition:'
  write ( *, '(a)' ) ' '

  call i4_partition_print ( n, npart2, a2, mult2 )

  return
end
subroutine test059 ( )

!*****************************************************************************80
!
!! TEST059 tests I4_PARTITION_NEXT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 7

  integer ( kind = 4 ) a(n)
  logical done
  integer ( kind = 4 ) mult(n)
  integer ( kind = 4 ) npart
  integer ( kind = 4 ) rank

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST059'
  write ( *, '(a)' ) '  I4_PARTITION_NEXT generates partitions of an integer.'
  write ( *, '(a,i8)' ) '  Here N = ', n
  write ( *, '(a)' ) ' '

  rank = 0
  done = .true.
 
  do

    call i4_partition_next ( n, npart, a, mult, done )
 
    if ( done ) then
      exit 
    end if

    rank = rank + 1

    call i4_partition_print ( n, npart, a, mult )

  end do
 
  return
end
subroutine test060 ( )

!*****************************************************************************80
!
!! TEST060 tests I4_PARTITION_NEXT2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 7

  integer ( kind = 4 ) a(n)
  logical more
  integer ( kind = 4 ) mult(n)
  integer ( kind = 4 ) npart

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST060'
  write ( *, '(a)' ) '  I4_PARTITION_NEXT2 produces partitions of an integer.'
  write ( *, '(a)' ) ' '

  more = .false.

  do

    call i4_partition_next2 ( n, npart, a, mult, more )

    call i4_partition_print ( n, npart, a, mult )

    if ( .not. more ) then
      exit
    end if

  end do
  
  return
end
subroutine test061 ( )

!*****************************************************************************80
!
!! TEST061 tests I4_PARTITION_COUNT and I4_PARTITION_COUNT_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ) p
  integer ( kind = 4 ) p2(0:n_max)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST061'
  write ( *, '(a)' ) '  I4_PARTITION_COUNT counts partitions of an integer.'
  write ( *, '(a)' ) '  I4_PARTITION_COUNT_VALUES returns some exact values.'

  n_data = 0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   N     Exact     Count'
  write ( *, '(a)' ) ' '

  do

    call i4_partition_count_values ( n_data, n, p )

    if ( n_data == 0 ) then
      exit
    end if

    call i4_partition_count ( n, p2 )
 
    write ( *, '(2x,i4,2i10)' ) n, p, p2(n)

  end do

  return
end
subroutine test0615 ( )

!*****************************************************************************80
!
!! TEST0615 tests I4_PARTITION_COUNT2 and I4_PARTITION_COUNT_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ) p
  integer ( kind = 4 ) p2(0:n_max)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0615'
  write ( *, '(a)' ) '  I4_PARTITION_COUNT2 counts partitions of an integer.'
  write ( *, '(a)' ) '  I4_PARTITION_COUNT_VALUES returns some exact values.'

  n_data = 0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   N     Exact     Count'
  write ( *, '(a)' ) ' '

  do

    call i4_partition_count_values ( n_data, n, p )

    if ( n_data == 0 ) then
      exit
    end if

    call i4_partition_count2 ( n, p2 )
 
    write ( *, '(2x,i4,2i10)' ) n, p, p2(n)

  end do

  return
end
subroutine test062 ( )

!*****************************************************************************80
!
!! TEST062 tests I4_PARTITION_RANDOM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 8

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) mult(n)
  integer ( kind = 4 ) npart
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) table(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST062'
  write ( *, '(a)' ) '  I4_PARTITION_RANDOM generates a random partition.'
  write ( *, '(a)' ) '  I4_PARTITION_COUNT2 sets up a partition count table'
  write ( *, '(a)' ) '  needed by I4_PARTITION_RANDOM.'
  write ( *, '(a)' ) ' '

  seed = 123456789
!
!  Get the partition table.
!
  call i4_partition_count2 ( n, table )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The number of partitions of N'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       N    Number of partitions'
  write ( *, '(a)' ) ' '

  do j = 1, n
    write ( *, '(2x,i8,4x,i8)' ) j, table(j)
  end do

  write ( *, '(a)' ) ' '

  do i = 1, 5

    call i4_partition_random ( n, table, seed, a, mult, npart )

    call i4_partition_print ( n, npart, a, mult )

  end do
 
  return
end
subroutine test06225 ( )

!*****************************************************************************80
!
!! TEST06225 tests I4_PARTITIONS_NEXT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 August 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: s = 3

  integer ( kind = 4 ) i
  integer ( kind = 4 ) m(s)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06225'
  write ( *, '(a)' ) '  I4_PARTITIONS_NEXT produces the next'
  write ( *, '(a)' ) '  nondecreasing partitions of an integer, and'
  write ( *, '(a)' ) '  if necessary, increments the integer to keep on going.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   I Sum    Partition'
  write ( *, '(a)' ) ' '
  i = 0

  m(1:3) = (/ 0, 0, 0 /)

  write ( *, '(2x,i2,2x,i2,4x,10i2)' ) i, sum ( m ), m(1:s)

  do i = 1, 15

    call i4_partitions_next ( s, m )

    write ( *, '(2x,i2,2x,i2,4x,10i2)' ) i, sum ( m ), m(1:s)

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  You can start from any legal partition.'
  write ( *, '(a)' ) '  Here, we restart at ( 2, 1, 0 ).'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   I Sum    Partition'
  write ( *, '(a)' ) ' '

  i = 0
  m = (/ 2, 1, 0 /)

  write ( *, '(2x,i2,2x,i2,4x,10i2)' ) i, sum ( m ), m(1:s)

  do i = 1, 15

    call i4_partitions_next ( s, m )

    write ( *, '(2x,i2,2x,i2,4x,10i2)' ) i, sum ( m ), m(1:s)

  end do

  return
end
subroutine test033 ( )

!*****************************************************************************80
!
!! TEST033 tests I4_SQRT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) q
  integer ( kind = 4 ) r

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST033'
  write ( *, '(a)' ) '  I4_SQRT computes the square root of an integer.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       N  Sqrt(N) Remainder'
  write ( *, '(a)' ) ' '

  do n = -5, 20

    call i4_sqrt ( n, q, r )
    write ( *, '(2x,3i9)' ) n, q, r

  end do

  return
end
subroutine test034 ( )

!*****************************************************************************80
!
!! TEST034 tests I4_SQRT_CF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: max_term = 100

  integer ( kind = 4 ) b(0:max_term)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_term

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST034'
  write ( *, '(a)' ) '  I4_SQRT_CF computes the continued fraction form'
  write ( *, '(a)' ) '  of the square root of an integer.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   N  Period  Whole  Repeating Part'
  write ( *, '(a)' ) ' '
  do n = 1, 20

    call i4_sqrt_cf ( n, max_term, n_term, b )
    write ( *, '(2x,i5,3x,i5,2x,i5,10i5)' ) n, n_term, b(0), b(1:n_term)

  end do

  return
end
subroutine test0625 ( )

!*****************************************************************************80
!
!! TEST0625 tests I4_TO_BVEC and BVEC_TO_I4;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  integer ( kind = 4 ) bvec(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0625'
  write ( *, '(a)' ) '  I4_TO_BVEC converts an integer to a '
  write ( *, '(a)' ) '  signed binary vector;'
  write ( *, '(a)' ) '  BVEC_TO_I4 converts a signed binary vector'
  write ( *, '(a)' ) '  to an integer;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I --> BVEC  -->  I'
  write ( *, '(a)' ) ' '
  do i = -3, 10
    call i4_to_bvec ( i, n, bvec )
    call bvec_to_i4 ( n, bvec, i2 )
    write ( *, '(2x,i3,2x,10i1,2x,i3)' ) i, bvec(1:n), i2
  end do

  return
end
subroutine test035 ( )

!*****************************************************************************80
!
!! TEST035 tests I4_TO_CHINESE and CHINESE_TO_I4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4

  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ), save, dimension ( n ) :: m = (/ 3, 4, 5, 7 /)
  integer ( kind = 4 ) r(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST035'
  write ( *, '(a)' ) '  I4_TO_CHINESE computes the Chinese Remainder'
  write ( *, '(a)' ) '  representation of an integer.'
  write ( *, '(a)' ) '  CHINESE_TO_I4 computes an integer with the given'
  write ( *, '(a)' ) '  Chinese Remainder representation.'

  call i4vec_print ( n, m, '  The moduli:' )

  j = 37

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number being analyzed is ', j

  call i4_to_chinese ( j, n, m, r )

  call i4vec_print ( n, r, '  The remainders:' )

  call chinese_to_i4 ( n, m, r, j2 )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The reconstructed number is ', j2

  call i4_to_chinese ( j2, n, m, r )

  call i4vec_print ( n, r, '  The remainders of the reconstructed number are:' )

  return
end
subroutine test0364 ( )

!*****************************************************************************80
!
!! TEST0364 tests I4_TO_I4POLY and I4_TO_VAN_DER_CORPUT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: degree_max = 10

  integer ( kind = 4 ) base
  integer ( kind = 4 ) degree
  real ( kind = 8 ) h1
  real ( kind = 8 ) h2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4poly(0:degree_max)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) :: n_test = 10
  integer ( kind = 4 ) prime
  integer ( kind = 4 ) value

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0364'
  write ( *, '(a)' ) '  I4_TO_VAN_DER_CORPUT computes the elements '
  write ( *, '(a)' ) '  of a van der Corput sequence.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I4_TO_I4POLY converts an integer to an integer '
  write ( *, '(a)' ) '  polynomial in some base, and can be used to mimic'
  write ( *, '(a)' ) '  the van der Corput calculation.'
  write ( *, '(a)' ) ' '

  do j = 1, 3

    base = prime(j)
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  BASE = ', base
    write ( *, '(a)' ) ' '

    do i = 1, n_test

      call i4_to_van_der_corput ( i, base, h1 )

      call i4_to_i4poly ( i, base, degree_max, degree, i4poly )

      call i4vec_reverse ( degree+1, i4poly )

      call i4poly_to_i4 ( degree, i4poly, base, value )

      h2 = real ( value, kind = 8 ) / base**( degree + 1 )

      write ( *, '(2x,i8,2x,g14.6,2x,g14.6)' ) i, h1,  h2

    end do

  end do

  return
end
subroutine test0627 ( )

!*****************************************************************************80
!
!! TEST0627 tests I4_TO_I4POLY and I4POLY_TO_I4;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: degree_max = 5
  integer ( kind = 4 ), parameter :: test_num = 9

  integer ( kind = 4 ) a(0:degree_max)
  integer ( kind = 4 ) base
  integer ( kind = 4 ), dimension ( test_num ) :: base_test = (/ &
   2, 2, 2, 3, 4, 5, 6, 23, 24 /)
  integer ( kind = 4 ) degree
  integer ( kind = 4 ) intval
  integer ( kind = 4 ) intval2
  integer ( kind = 4 ), dimension ( test_num ) :: intval_test = (/ &
   1, 6, 23, 23, 23, 23, 23, 23, 23 /)
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0627'
  write ( *, '(a)' ) '  I4_TO_I4POLY converts an integer to a polynomial'
  write ( *, '(a)' ) '  in a given base;'
  write ( *, '(a)' ) '  I4POLY_TO_I4 evaluates an integer polynomial.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       I    BASE  DEGREE  Coefficients'
  write ( *, '(a)' ) ' '
  do test = 1, test_num
    intval = intval_test(test)
    base = base_test(test)
    call i4_to_i4poly ( intval, base, degree_max, degree, a )
    write ( *, '(2x,i4,2x,i4,2x,i4,6x,10i8)' ) &
      intval, base, degree, a(degree:0:-1)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now let I4_TO_I4POLY convert I to a polynomial,'
  write ( *, '(a)' ) '  use I4POLY_TO_I4 to evaluate it, and compare.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       I    I2'
  write ( *, '(a)' ) ' '
  do test = 1, test_num
    intval = intval_test(test)
    base = base_test(test)
    call i4_to_i4poly ( intval, base, degree_max, degree, a )
    call i4poly_to_i4 ( degree, a, base, intval2 )
    write ( *, '(2x,i8,2x,i8)' ) intval, intval2
  end do

  return
end
subroutine test036 ( )

!*****************************************************************************80
!
!! TEST036 tests I4_TO_VAN_DER_CORPUT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_prime = 5

  real ( kind = 8 ) h(n_prime)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) :: n_test = 10
  integer ( kind = 4 ) p
  integer ( kind = 4 ) prime

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST036'
  write ( *, '(a)' ) '  I4_TO_VAN_DER_CORPUT computes the elements '
  write ( *, '(a)' ) '  of a van der Corput sequence.'
  write ( *, '(a)' ) '  The sequence depends on the prime numbers used '
  write ( *, '(a)' ) '  as a base.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Bases:'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,5(8x,i8))' ) ( prime(j), j = 1, n_prime )
  write ( *, '(a)' ) ' '

  do i = 1, n_test
    do j = 1, n_prime
      p = prime(j)
      call i4_to_van_der_corput ( i, p, h(j) )
    end do
    write ( *, '(2x,i3,5g14.6)' ) i, h(1:n_prime)
  end do

  return
end
subroutine test0365 ( )

!*****************************************************************************80
!
!! TEST0365 tests IBSET.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) bit
  integer ( kind = 4 ), parameter :: i = 101
  integer ( kind = 4 ) j

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0365'
  write ( *, '(a)' ) '  IBSET sets a bit to 1.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '           I Bit  IBSET(I,BIT)'
  write ( *, '(a)' ) ' '

  do bit = 0, 31

    j = ibset ( i, bit )

    write ( *, '(2x,i12,2x,i2,2x,i12)' ) i, bit, j

  end do

  return
end
subroutine test037 ( )

!*****************************************************************************80
!
!! TEST037 tests I4MAT_01_ROWCOLSUM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: n = 5

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ), save, dimension ( n ) :: c = (/ 2, 2, 2, 2, 1 /)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ), save, dimension ( m ) :: r = (/ 3, 2, 2, 1, 1 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST037'
  write ( *, '(a)' ) '  I4MAT_01_ROWCOLSUM constructs a 01 matrix with'
  write ( *, '(a)' ) '  given row and column sums.'
  
  call i4vec_print ( m, r, '  The rowsum vector:' )
  call i4vec_print ( n, c, '  The columnsum vector: ' )

  call i4mat_01_rowcolsum ( m, n, r, c, a, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) &
      '  I4MAT_01_ROWCOLSUM returned error flag IERROR = ', ierror
  else
    call i4mat_print ( m, n, a, '  The rowcolsum matrix:' )
  end if

  return
end
subroutine test038 ( )

!*****************************************************************************80
!
!! TEST038 tests I4MAT_01_ROWCOLSUM2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 5
  integer ( kind = 4 ), parameter :: n = 5

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ), save, dimension ( n ) :: c
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ), save, dimension ( m ) :: r

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST038'
  write ( *, '(a)' ) '  I4MAT_01_ROWCOLSUM2 constructs a 01 matrix with'
  write ( *, '(a)' ) '  given row and column sums.'
  
  c = (/ 2, 1, 2, 2, 2 /)
  r = (/ 2, 1, 3, 1, 2 /)

  call i4vec_print ( m, r, '  The rowsum vector:' )
  call i4vec_print ( n, c, '  The columnsum vector: ' )

  call i4mat_01_rowcolsum2 ( m, n, r, c, a, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) &
      '  I4MAT_01_ROWCOLSUM2 returned error flag IERROR = ', ierror
    write ( *, '(a)' ) '  The matrix returned is not an exact solution.'
  end if

  call i4mat_print ( m, n, a, '  The rowcolsum matrix:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now repeat, with data for which there is not'
  write ( *, '(a)' ) '  a solution.  The program will try its best anyway.'

  c = (/ 1, 4, 1, 5, 1 /)
  r = (/ 1, 3, 4, 1, 3 /)

  call i4vec_print ( m, r, '  The rowsum vector:' )
  call i4vec_print ( n, c, '  The columnsum vector: ' )

  call i4mat_01_rowcolsum2 ( m, n, r, c, a, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) &
      '  I4MAT_01_ROWCOLSUM2 returned error flag IERROR = ', ierror
    write ( *, '(a)' ) '  The matrix returned is not an exact solution.'
  end if

  call i4mat_print ( m, n, a, '  The rowcolsum matrix:' )

  return
end
subroutine test039 ( )

!*****************************************************************************80
!
!! TEST039 tests I4MAT_01_ROWCOLSUM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 13
  integer ( kind = 4 ), parameter :: n = 17

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ), save, dimension ( n ) :: c = &
    (/  4,  4, 11, 10, 10,  8,  9, 10,  8,  9,  3, 10,  4,  7,  9,  3,  3 /)
  integer ( kind = 4 ), dimension ( n ) :: c2
  integer ( kind = 4 ), dimension ( n ) :: c3
  integer ( kind = 4 ), dimension ( n ) :: cperm
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  integer ( kind = 4 ), save, dimension ( m ) :: r = &
    (/ 14, 13, 14, 10, 12,  2, 10,  1, 10, 11,  6,  2, 17 /)
  integer ( kind = 4 ), dimension ( m ) :: r2
  integer ( kind = 4 ), dimension ( m ) :: r3
  integer ( kind = 4 ), dimension ( m ) :: rperm

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST039'
  write ( *, '(a)' ) '  I4MAT_01_ROWCOLSUM constructs a 01 matrix with'
  write ( *, '(a)' ) '  given row and column sums.'

  call i4vec_print ( m, r, '  The rowsum vector R:' )
  call i4vec_sort_heap_index_d ( m, r, rperm )

  do i = 1, m
    r2(i) = r(rperm(i))
  end do

  call i4vec_print ( n, c, '  The columnsum vector C: ' )
  call i4vec_sort_heap_index_d ( n, c, cperm )

  do i = 1, n
    c2(i) = c(cperm(i))
  end do

  call i4mat_01_rowcolsum ( m, n, r2, c2, a, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) &
      '  I4MAT_01_ROWCOLSUM returned error flag IERROR = ', ierror
    return
  end if
!
!  Invert the R and C permutations.
!
  call i4mat_perm2 ( m, n, a, rperm, cperm )

  call i4mat_print ( m, n, a, '  The rowcolsum matrix:' )

  do i = 1, m
    r3(i) = sum ( a(i,1:n) )
  end do

  call i4vec_print ( m, r3, '  Computed row sums' )

  do j = 1, n
    c3(j) = sum ( a(1:m,j) )
  end do

  call i4vec_print ( n, c3, '  Computed column sums:' )

  return
end
subroutine test040 ( )

!*****************************************************************************80
!
!! TEST040 tests I4MAT_U1_INVERSE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 6

  integer ( kind = 4 ), dimension (n,n) :: a = reshape ( &
  (/ &
    1, 0, 0, 0, 0, 0, &
    1, 1, 0, 0, 0, 0, &
    0, 0, 1, 0, 0, 0, &
    0, 0, 1, 1, 0, 0, &
    0, 0, 0, 0, 1, 0, &
   75, 0, 0, 0, 1, 1 /), (/ n, n /) )
  integer ( kind = 4 ) b(n,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST040'
  write ( *, '(a)' ) '  I4MAT_U1_INVERSE inverts a unit upper triangular matrix.'

  call i4mat_print ( n, n, a, '  The input matrix:' )
 
  call i4mat_u1_inverse ( n, a, b )
 
  call i4mat_print ( n, n, b, '  The inverse:' )
 
  return
end
subroutine test041 ( )

!*****************************************************************************80
!
!! TEST041 test I4MAT_PERM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 9

  integer ( kind = 4 ) a(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ), save, dimension ( n ) :: p = (/ 2,3,9,6,7,8,5,4,1 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST041'
  write ( *, '(a)' ) '  I4MAT_PERM reorders an integer matrix in place.'
  write ( *, '(a)' ) '  The rows and columns use the same permutation.'

  do i = 1, n
    do j = 1, n
      a(i,j) = i * 10 + j
    end do
  end do

  call i4mat_print ( n, n, a, '  The input matrix:' )
 
  call perm_print ( n, p, '  The row and column permutation:' )
 
  call i4mat_perm ( n, a, p )
 
  call i4mat_print ( n, n, a, '  The permuted matrix:' )
 
  return
end
subroutine test042 ( )

!*****************************************************************************80
!
!! TEST042 test I4MAT_PERM2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 9
  integer ( kind = 4 ), parameter :: n = 7

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ), save, dimension ( m ) :: p = (/ 2,3,9,6,7,8,5,4,1 /)
  integer ( kind = 4 ), save, dimension ( n ) :: q = (/ 3,4,5,6,7,1,2 /)
  integer ( kind = 4 ) j

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST042'
  write ( *, '(a)' ) '  I4MAT_PERM2 reorders an integer matrix in place.'
  write ( *, '(a)' ) '  Rows and columns use different permutations.'

  do i = 1, m
    do j = 1, n
      a(i,j) = i * 10 + j
    end do
  end do
 
  call i4mat_print ( m, n, a, '  The input matrix:' )
 
  call perm_print ( m, p, '  The row permutation:' )

  call perm_print ( n, q, '  The column permutation:' )
 
  call i4mat_perm2 ( m, n, a, p, q )
 
  call i4mat_print ( m, n, a, '  The permuted matrix:' )

  return
end
subroutine test047 ( )

!*****************************************************************************80
!
!! TEST047 tests I4MAT_PERM2 and TRIANG.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  integer ( kind = 4 ), dimension ( n, n ) :: a = reshape ( &
  (/ &
    1,0,1,0,1,0,1,0,0,1, &
    0,1,0,0,1,0,0,0,0,0, &
    0,0,1,0,1,0,1,0,0,1, &
    0,1,1,1,1,1,1,1,0,1, &
    0,0,0,0,1,0,0,0,0,0, &
    0,1,0,0,1,1,1,0,0,0, &
    0,0,0,0,1,0,1,0,0,0, &
    0,1,0,0,1,1,1,1,0,1, &
    0,0,0,0,0,0,0,0,0,0, &
    0,0,0,0,1,0,1,0,0,1 /), (/ n, n /) )
  integer ( kind = 4 ) p(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST047'
  write ( *, '(a)' ) '  TRIANG relabels elements for a partial ordering,'
  write ( *, '(a)' ) '  I4MAT_PERM2 reorders an integer matrix in place.'

  call i4mat_print ( n, n, a, '  The input matrix:' )
 
  call triang ( n, a, p )
 
  call perm_print ( n, p, '  The new ordering:' )

  call i4mat_perm2 ( n, n, a, p, p )
 
  call i4mat_print ( n, n, a, '  The reordered matrix:' )
 
  return
end
subroutine test043 ( )

!*****************************************************************************80
!
!! TEST043 tests INDEX_BOX_NEXT_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  logical more
  integer ( kind = 4 ), parameter :: n1 = 5
  integer ( kind = 4 ), parameter :: n2 = 3
  integer ( kind = 4 ) n

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST043'
  write ( *, '(a)' ) '  INDEX_BOX_NEXT_2D produces IJ indices that'
  write ( *, '(a)' ) '  lie on the surface of a box in 2D.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The box has logical dimensions:'
  write ( *, '(3x,3i3)' ) n1, n2
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   #    I   J'
  write ( *, '(a)' ) ' '

  more = .false.
  n = 0

  do

    call index_box_next_2d ( n1, n2, i, j, more )

    if ( .not. more ) then
      exit
    end if

    n = n + 1
    write ( *, '(2x,4i3)' ) n, i, j

  end do

  return
end
subroutine test044 ( )

!*****************************************************************************80
!
!! TEST044 tests INDEX_BOX_NEXT_3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  logical more
  integer ( kind = 4 ), parameter :: n1 = 5
  integer ( kind = 4 ), parameter :: n2 = 3
  integer ( kind = 4 ), parameter :: n3 = 4
  integer ( kind = 4 ) n

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST044'
  write ( *, '(a)' ) '  INDEX_BOX_NEXT_3D produces IJK indices that'
  write ( *, '(a)' ) '  lie on the surface of a box.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The box has logical dimensions:'
  write ( *, '(3x,3i3)' ) n1, n2, n3
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     #    I   J   K'
  write ( *, '(a)' ) ' '

  more = .false.
  n = 0

  do

    call index_box_next_3d ( n1, n2, n3, i, j, k, more )

    if ( .not. more ) then
      exit
    end if

    n = n + 1
    write ( *, '(2x,4i3)' ) n, i, j, k

  end do

  return
end
subroutine test045 ( )

!*****************************************************************************80
!
!! TEST045 tests INDEX_BOX2_NEXT_2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: ic = 10
  integer ( kind = 4 ) j
  integer ( kind = 4 ), parameter :: jc = 20
  logical more
  integer ( kind = 4 ), parameter :: n1 = 4
  integer ( kind = 4 ), parameter :: n2 = 3
  integer ( kind = 4 ) n

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST045'
  write ( *, '(a)' ) '  INDEX_BOX2_NEXT_2D produces IJ indices that'
  write ( *, '(a)' ) '  lie on the surface of a box2 in 2D.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The box has half-widths:'
  write ( *, '(3x,3i3)' ) n1, n2
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  and has center cell:'
  write ( *, '(3x,2i3)' ) ic, jc
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     #    I   J'
  write ( *, '(a)' ) ' '

  more = .false.
  n = 0

  do

    call index_box2_next_2d ( n1, n2, ic, jc, i, j, more )

    if ( .not. more ) then
      exit
    end if

    n = n + 1
    write ( *, '(2x,4i3)' ) n, i, j

  end do

  return
end
subroutine test046 ( )

!*****************************************************************************80
!
!! TEST046 tests INDEX_BOX2_NEXT_3D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: ic = 10
  integer ( kind = 4 ) j
  integer ( kind = 4 ), parameter :: jc = 20
  integer ( kind = 4 ) k
  integer ( kind = 4 ), parameter :: kc = 30
  logical more
  integer ( kind = 4 ), parameter :: n1 = 5
  integer ( kind = 4 ), parameter :: n2 = 3
  integer ( kind = 4 ), parameter :: n3 = 4
  integer ( kind = 4 ) n

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST046'
  write ( *, '(a)' ) '  INDEX_BOX2_NEXT_3D produces IJK indices that'
  write ( *, '(a)' ) '  lie on the surface of a box.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The box has half widths:'
  write ( *, '(3x,3i3)' ) n1, n2, n3
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  and central cell:'
  write ( *, '(3x,3i3)' ) ic, jc, kc
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  We will only print a PORTION of the data!'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   #    I   J   K'
  write ( *, '(a)' ) ' '

  more = .false.
  n = 0

  do

    call index_box2_next_3d ( n1, n2, n3, ic, jc, kc, i, j, k, more )

    if ( .not. more ) then
      exit
    end if

    n = n + 1

    if ( n <= 10 .or. 370 <= n ) then
      write ( *, '(2x,4i3)' ) n, i, j, k
    end if

  end do

  return
end
subroutine test048 ( )

!*****************************************************************************80
!
!! TEST048 tests INDEX_NEXT0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ), parameter :: hi = 3
  logical more

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST048'
  write ( *, '(a)' ) '  INDEX_NEXT0 generates all indices of an'
  write ( *, '(a)' ) '  array of given shape, with'
  write ( *, '(a)' ) '  lower limit 1 and given upper limit.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of index entries = ', n
  write ( *, '(a,i8)' ) '  Coordinate maximum HI =   ', hi
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Index arrays:'
  write ( *, '(a)' ) ' '

  more = .false.

  do

    call index_next0 ( n, hi, a, more )

    write ( *, '(2x,3i4)' ) a(1:n)

    if ( .not. more ) then
      exit
    end if

  end do

  return
end
subroutine test049 ( )

!*****************************************************************************80
!
!! TEST049 tests INDEX_NEXT1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ), dimension ( n ) :: hi = (/ 4, 2, 3 /)
  logical more

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST049'
  write ( *, '(a)' ) '  INDEX_NEXT1 generates all indices of an'
  write ( *, '(a)' ) '  array of given shape, with'
  write ( *, '(a)' ) '  lower limit 1 and given upper limits.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of index entries = ', n

  call i4vec_print ( n, hi, '  Coordinate maximum indices:' )
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Index arrays:'
  write ( *, '(a)' ) ' '

  more = .false.

  do

    call index_next1 ( n, hi, a, more )

    write ( *, '(2x,3i4)' ) a(1:n)

    if ( .not. more ) then
      exit
    end if

  end do

  return
end
subroutine test050 ( )

!*****************************************************************************80
!
!! TEST050 tests INDEX_NEXT2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ), dimension ( n ) :: hi = (/ 11, -3, 1 /)
  integer ( kind = 4 ) i
  integer ( kind = 4 ), dimension ( n ) :: lo = (/ 10, -5, 0 /)
  logical more

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST050'
  write ( *, '(a)' ) '  INDEX_NEXT2 generates all indices of an'
  write ( *, '(a)' ) '  array of given shape with given'
  write ( *, '(a)' ) '  lower and upper limits.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of index entries = ', n
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Coordinate, Maximum Index'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(3i10)' ) i, lo(i), hi(i)
  end do
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Index arrays:'
  write ( *, '(a)' ) ' '

  more = .false.

  do

    call index_next2 ( n, lo, hi, a, more )

    write ( *, '(2x,3i4)' ) a(1:n)

    if ( .not. more ) then
      exit
    end if

  end do

  return
end
subroutine test051 ( )

!*****************************************************************************80
!
!! TEST051 tests INDEX_RANK0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3

  integer ( kind = 4 ), dimension ( n ) :: a = (/ 3, 1, 2 /)
  integer ( kind = 4 ), parameter :: hi = 3
  integer ( kind = 4 ) rank

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST051'
  write ( *, '(a)' ) '  INDEX_RANK0 ranks an index with'
  write ( *, '(a)' ) '  lower limit 1 and given upper limit.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of index entries = ', n
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Coordinate maximum Index = ', hi
  write ( *, '(a)' ) ' '

  call i4vec_print ( n, a, '  The index array:' )

  call index_rank0 ( n, hi, a, rank )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The rank of this object is ', rank

  return
end
subroutine test052 ( )

!*****************************************************************************80
!
!! TEST052 tests INDEX_RANK1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3

  integer ( kind = 4 ), dimension ( n ) :: a = (/ 4, 1, 2 /)
  integer ( kind = 4 ), dimension ( n ) :: hi = (/ 4, 2, 3 /)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) rank

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST052'
  write ( *, '(a)' ) '  INDEX_RANK1 ranks an index with'
  write ( *, '(a)' ) '  lower limit 1 and given upper limits.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of index entries = ', n
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Coordinate, Maximum Index'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,2i10)' ) i, hi(i)
  end do
 
  call i4vec_print ( n, a, '  The index array:' )

  call index_rank1 ( n, hi, a, rank )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The rank of this object is ', rank

  return
end
subroutine test053 ( )

!*****************************************************************************80
!
!! TEST053 tests INDEX_RANK2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3

  integer ( kind = 4 ), dimension ( n ) :: a = (/ 1, 11, 5 /)
  integer ( kind = 4 ), dimension ( n ) :: hi = (/ 2, 11, 6 /)
  integer ( kind = 4 ) i
  integer ( kind = 4 ), dimension ( n ) :: lo = (/ 1, 10, 4 /)
  integer ( kind = 4 ) rank

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST053'
  write ( *, '(a)' ) '  INDEX_RANK2 ranks an index with given'
  write ( *, '(a)' ) '  lower and upper limits.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of index entries = ', n
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Coordinate, Minimum index, Maximum Index'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,3i10)' ) i, lo, hi(i)
  end do
 
  call i4vec_print ( n, a, '  The index array:' )

  call index_rank2 ( n, lo, hi, a, rank )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The rank of this object is ', rank

  return
end
subroutine test054 ( )

!*****************************************************************************80
!
!! TEST054 tests INDEX_UNRANK0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ), parameter :: hi = 3
  integer ( kind = 4 ) maxrank
  integer ( kind = 4 ) rank

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST054'
  write ( *, '(a)' ) '  INDEX_UNRANK0 unranks a multi-index.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The multi-index has dimension ', n
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The upper limit is HI = ', hi
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Rank, Multi-Index:'
  write ( *, '(a)' ) ' '
 
  maxrank = hi**n

  do rank = 1, maxrank
    call index_unrank0 ( n, hi, rank, a )
    write ( *, '(2x,i3,3i8)' ) rank, a(1:n)
  end do
 
  return
end
subroutine test055 ( )

!*****************************************************************************80
!
!! TEST055 tests INDEX_UNRANK1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ), dimension ( n ) :: hi = (/ 4, 2, 3 /)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) maxrank
  integer ( kind = 4 ) rank

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST055'
  write ( *, '(a)' ) '  INDEX_UNRANK1 unranks a multi-index.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The multi-index has dimension ', n
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The upper limits are:'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,2i10)' ) i, hi(i)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Rank, Multi-Index:'
  write ( *, '(a)' ) ' '
 
  maxrank = product ( hi )

  do rank = 1, maxrank
    call index_unrank1 ( n, hi, rank, a )
    write ( *, '(2x,i3,3i8)' ) rank, a(1:n)
  end do
 
  return
end
subroutine test056 ( )

!*****************************************************************************80
!
!! TEST056 tests INDEX_UNRANK2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ), dimension ( n ) :: hi = (/ 2, 11, 6 /)
  integer ( kind = 4 ) i
  integer ( kind = 4 ), dimension (n) :: lo = (/ 1, 10, 4 /)
  integer ( kind = 4 ) rank

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST056'
  write ( *, '(a)' ) '  INDEX_UNRANK2 unranks a multi-index.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The multi-index has dimension ', n
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The lower and upper limits are:'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,3i10)' ) i, lo(i), hi(i)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Rank, Multi-Index:'
  write ( *, '(a)' ) ' '
 
  rank = 7

  call index_unrank2 ( n, lo, hi, rank, a )
  write ( *, '(2x,i3,3i8)' ) rank, a(1:n)
 
  return
end
subroutine test057 ( )

!*****************************************************************************80
!
!! TEST057 tests INS_PERM and PERM_INS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ins(n)
  integer ( kind = 4 ) perm(n)
  integer ( kind = 4 ) perm2(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST057'
  write ( *, '(a)' ) '  PERM_INS computes the inversion sequence.'
  write ( *, '(a)' ) '  INS_PERM recovers the permutation.'
  write ( *, '(a)' ) ' '

  perm(1:n) = (/ 3, 5, 1, 4, 2 /)

  call perm_ins ( n, perm, ins )

  call ins_perm ( n, ins, perm2 )

  write ( *, '(2x,6i3)' ) ( i, i = 1, n )
  write ( *, '(2x,6i3)' ) perm(1:n)
  write ( *, '(2x,6i3)' ) ins(1:n)
  write ( *, '(2x,6i3)' ) perm2(1:n)
 
  return
end
subroutine test063 ( )

!*****************************************************************************80
!
!! TEST063 tests INVOLUTE_ENUM;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  integer ( kind = 4 ) i
  real ( kind = 8 ) prob
  real ( kind = 8 ) r8_factorial
  integer ( kind = 4 ) s(0:n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST063'
  write ( *, '(a)' ) '  INVOLUTE_ENUM counts involutions;'
  write ( *, '(a)' ) ' '

  call involute_enum ( n, s )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N    Number   Probability'
  write ( *, '(a)' ) ' '

  do i = 0, n
    prob = real ( s(i), kind = 8 ) / r8_factorial ( i )
    write ( *, '(2x,i8,2x,i8,2x,g14.6)' ) i, s(i), prob
  end do

  return
end
subroutine test064 ( )

!*****************************************************************************80
!
!! TEST064 test I4POLY.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 6
  integer ( kind = 4 ), parameter :: test_num = 6

  integer ( kind = 4 ), dimension ( n ) :: a
  integer ( kind = 4 ) iopt
  integer ( kind = 4 ) test
  integer ( kind = 4 ) val
  integer ( kind = 4 ) x0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST064'
  write ( *, '(a)' ) '  I4POLY converts between power sum, factorial'
  write ( *, '(a)' ) '  and Taylor forms, and can evaluate a polynomial'
  write ( *, '(a)' ) ' '
 
  do test = 1, test_num

    if ( test == 1 ) then
      iopt = -3
    else if ( test == 2 ) then
      iopt = -2
    else if ( test == 3 ) then
      iopt = -1
      x0 = 2
    else if ( test == 4 ) then
      iopt = 0
      x0 = 2
    else if ( test == 5 ) then
      iopt = 6
      x0 = 2
    else if ( test == 6 ) then
      iopt = 6
      x0 = -2
    end if

    a = (/ 0, 0, 0, 0, 0, 1 /)

    if ( test == 1 ) then
      write ( *, '(a)' ) '  All calls have input A as follows'
      write ( *, '(10i4)' ) a(1:n)
    end if
 
    call i4poly ( n, a, x0, iopt, val )
 
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Option IOPT = ', iopt
    if ( -1 <= iopt ) then
      write ( *, '(a,i8)' ) '  X0 = ', x0
    end if
    if ( iopt == -3 .or. iopt == -2 .or. 0 < iopt ) then
      write ( *, '(a)' ) '  Output array = '
      write ( *, '(2x,10i4)' ) a(1:n)
    end if
    if ( iopt == -1 .or. iopt == 0 ) then
      write ( *, '(a,i8)' ) '  Value = ', val
    end if
 
  end do

  return
end
subroutine test065 ( )

!*****************************************************************************80
!
!! TEST065 tests I4POLY_CYCLO.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: max_n = 10

  integer ( kind = 4 ) phi(0:max_n)
  integer ( kind = 4 ) n

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST065'
  write ( *, '(a)' ) '  I4POLY_CYCLO computes cyclotomic polynomials.'

  do n = 0, max_n

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  N = ', n
    write ( *, '(a)' ) ' '
    call i4poly_cyclo ( n, phi )

    call i4poly_print ( n, phi, '  The cyclotomic polynomial:' )

  end do

  return
end
subroutine test066 ( )

!*****************************************************************************80
!
!! TEST066 tests I4POLY_DIV.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) a(0:10)
  integer ( kind = 4 ) b(0:10)
  integer ( kind = 4 ) na
  integer ( kind = 4 ) nb
  integer ( kind = 4 ) nq
  integer ( kind = 4 ) nr
  integer ( kind = 4 ), parameter :: test_num = 2
  integer ( kind = 4 ) q(0:10)
  integer ( kind = 4 ) r(0:10)
  integer ( kind = 4 ) test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST066'
  write ( *, '(a)' ) '  I4POLY_DIV computes the quotient and'
  write ( *, '(a)' ) '  remainder for polynomial division.'
  write ( *, '(a)' ) ' '
!
!  1: Divide X**3 + 2*X**2 - 5*X - 6  by X-2.  
!     Quotient is 3+4*X+X**2, remainder is 0.
!
!  2: Divide X**4 + 3*X**3 + 2*X**2 - 2  by  X**2 + X - 3.
!     Quotient is X**2 + 2*X + 3, remainder 8*X + 7.
!
  do test = 1, test_num

    if ( test == 1 ) then
      na = 3
      a(0:na) = (/ -6, -5, 2, 1 /)
      nb = 1
      b(0:nb) = (/ -2, 1 /)
    else if ( test == 2 ) then
      na = 4
      a(0:na) = (/ -2, 5, 2, 3, 1 /)
      nb = 2
      b(0:nb) = (/ -3, 1, 1 /)
    end if

    call i4poly_print ( na, a, '  The polynomial to be divided, A:' )
    call i4poly_print ( nb, b, '  The divisor polynomial, B:' )

    call i4poly_div ( na, a, nb, b, nq, q, nr, r )
 
    call i4poly_print ( nq, q, '  The quotient polynomial, Q:' )
    call i4poly_print ( nr, r, '  The remainder polynomial, R:' )

  end do

  return
end
subroutine test067 ( )

!*****************************************************************************80
!
!! TEST067 tests I4POLY_MUL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: maxn = 5

  integer ( kind = 4 ) a(0:maxn)
  integer ( kind = 4 ) b(0:maxn)
  integer ( kind = 4 ) c(0:maxn)
  integer ( kind = 4 ) test
  integer ( kind = 4 ) na
  integer ( kind = 4 ) nb
  integer ( kind = 4 ), parameter :: test_num = 2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST067'
  write ( *, '(a)' ) '  I4POLY_MUL multiplies two polynomials.'
  write ( *, '(a)' ) ' '
!
!  1: Multiply (1+X) times (1-X).  Answer is 1-X**2.
!  2: Multiply (1+2*X+3*X**2) by (1-2*X). Answer is 1 + 0*X - X**2 - 6*X**3
!
  do test = 1, test_num

    if ( test == 1 ) then
      na = 1
      a(0:na) = (/ 1, 1 /)
      nb = 1
      b(0:nb) = (/ 1, -1 /)
    else if ( test == 2 ) then
      na = 2
      a(0:na) = (/ 1, 2, 3 /)
      nb = 1
      b(0:nb) = (/ 1, -2 /)
    end if

    call i4poly_mul ( na, a, nb, b, c )

    call i4poly_print ( na, a, '  The factor A:' )

    call i4poly_print ( nb, b, '  The factor B:' )

    call i4poly_print ( na+nb, c, '  The product C = A*B:' )

  end do

  return
end
subroutine test0675 ( )

!*****************************************************************************80
!
!! TEST0675 tests I4VEC_DESCENDS;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4

  integer ( kind = 4 ) a(n)
  logical descends
  integer ( kind = 4 ) i
  logical i4vec_descends
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0675'
  write ( *, '(a)' ) '  I4VEC_DESCENDS is true if an integer vector decreases.'
  write ( *, '(a)' ) ' '

  seed = 123456789

  do i = 1, 5

    call i4vec_uniform ( n, 1, n, seed, a )

    call i4vec_print ( n, a, '  The integer array to search:' )
 
    descends = i4vec_descends ( n, a )

    if ( descends ) then
      write ( *, '(a)' ) '  The preceding vector is descending.'
    else
      write ( *, '(a)' ) '  The preceding vector is not descending.'
    end if

  end do

  return
end
subroutine test068 ( )

!*****************************************************************************80
!
!! TEST068 tests I4VEC_FRAC;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) afrac
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST068'
  write ( *, '(a)' ) '  I4VEC_FRAC: K-th smallest integer vector entry.'
  write ( *, '(a)' ) ' '

  seed = 123456789

  call i4vec_uniform ( n, 1, 2*n, seed, a )

  call i4vec_print ( n, a, '  The integer array to search:' )
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       K       K-th smallest'
  write ( *, '(a)' ) ' '

  do k = 1, n

    call i4vec_frac ( n, a, k, afrac )
    write ( *, '(2x,i8,2x,i8)' ) k, afrac

  end do

  return
end
subroutine test0683 ( )

!*****************************************************************************80
!
!! TEST0683 tests I4VEC_INDEX;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 20

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) aval
  integer ( kind = 4 ) first
  integer ( kind = 4 ) i4vec_index
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0683'
  write ( *, '(a)' ) '  I4VEC_INDEX returns the index of the first occurrence'
  write ( *, '(a)' ) '  of a given value in an integer vector.'
  write ( *, '(a)' ) ' '

  seed = 123456789

  call i4vec_uniform ( n, 1, n/2, seed, a )

  aval = a(n/2)

  call i4vec_print ( n, a, '  The integer array to search:' )
 
  first = i4vec_index ( n, a, aval )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The value searched for is ', aval
  write ( *, '(a,i8)' ) '  The index of first occurrence is ', first

  return
end
subroutine test0685 ( )

!*****************************************************************************80
!
!! TEST0685 tests I4VEC_MAXLOC_LAST;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 20

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i4vec_maxloc_last
  integer ( kind = 4 ) last
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0685'
  write ( *, '(a)' ) '  I4VEC_MAXLOC_LAST: index of the last maximal'
  write ( *, '(a)' ) '  entry in an integer vector.'
  write ( *, '(a)' ) ' '

  seed = 123456789

  call i4vec_uniform ( n, 1, n/4, seed, a )

  call i4vec_print ( n, a, '  The integer array to search:' )
 
  last = i4vec_maxloc_last ( n, a )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Index of last maximal entry is ', last

  return
end
subroutine test0686 ( )

!*****************************************************************************80
!
!! TEST0686 tests I4VEC_PAIRWISE_PRIME;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  logical i4vec_pairwise_prime
  logical pairwise_prime
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0686'
  write ( *, '(a)' ) '  I4VEC_PAIRWISE_PRIME is true if an integer vector'
  write ( *, '(a)' ) '  is pairwise prime.'
  write ( *, '(a)' ) ' '

  seed = 123456789

  do i = 1, 5

    call i4vec_uniform ( n, 1, n, seed, a )

    call i4vec_print ( n, a, '  The integer array to check:' )
 
    pairwise_prime = i4vec_pairwise_prime ( n, a )

    if ( pairwise_prime ) then
      write ( *, '(a)' ) '  The preceding vector is pairwise prime.'
    else
      write ( *, '(a)' ) '  The preceding vector is not pairwise prime.'
    end if

  end do

  return
end
subroutine test0687 ( )

!*****************************************************************************80
!
!! TEST0687 tests I4VEC_REVERSE;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0687'
  write ( *, '(a)' ) '  I4VEC_REVERSE reverses an integer vector.'
  write ( *, '(a)' ) ' '

  seed = 123456789

  call i4vec_uniform ( n, 1, n, seed, a )

  call i4vec_print ( n, a, '  The integer array:' )
 
  call i4vec_reverse ( n, a )

  call i4vec_print ( n, a, '  The reversed integer array:' )

  return
end
subroutine test0688 ( )

!*****************************************************************************80
!
!! TEST0688 tests I4VEC_SORT_BUBBLE_A.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 20

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0688'
  write ( *, '(a)' ) '  I4VEC_SORT_BUBBLE_A ascending sorts an integer vector'
  write ( *, '(a)' ) '  using bubble sort.'
  write ( *, '(a)' ) ' '

  seed = 123456789

  call i4vec_uniform ( n, 0, 3*n, seed, a )

  call i4vec_print ( n, a, '  Unsorted array:' )

  call i4vec_sort_bubble_a ( n, a )

  call i4vec_print ( n, a, '  Sorted array:' )

  return
end
subroutine test0689 ( )

!*****************************************************************************80
!
!! TEST0689 tests I4VEC_SORT_HEAP_INDEX_D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 20

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(n)
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0689'
  write ( *, '(a)' ) '  I4VEC_SORT_HEAP_INDEX_D descending index-sorts'
  write ( *, '(a)' ) '  an integer vector using heap sort.'
  write ( *, '(a)' ) ' '

  seed = 123456789

  call i4vec_uniform ( n, 0, 3*n, seed, a )

  call i4vec_print ( n, a, '  Unsorted array:' )

  call i4vec_sort_heap_index_d ( n, a, indx )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       I    INDX    A(INDX)'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,3i8)' ) i, indx(i), a(indx(i))
  end do

  return
end
subroutine test06895 ( )

!*****************************************************************************80
!
!! TEST06895 tests INVERSE_MOD_N.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 November 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) b
  integer ( kind = 4 ) n
  integer ( kind = 4 ) y
  integer ( kind = 4 ) z

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06895'
  write ( *, '(a)' ) '  INVERSE_MOD_N seeks Y, the inverse of B mod N,'
  write ( *, '(a)' ) '  so that mod ( B * Y, N ) = 1, but returns 0'
  write ( *, '(a)' ) '  if the inverse does not exist.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     B     N     Y     Z = mod ( B * Y, N )'
  write ( *, '(a)' ) ' '

  do n = 1, 10
    do b = 1, n - 1
      call inverse_mod_n ( b, n, y )
      z = mod ( b * y, n )
      write ( *, '(2x,i4,2x,i4,2x,i4,2x,i4)' ) b, n, y, z
    end do
  end do

  return
end
subroutine test069 ( )

!*****************************************************************************80
!
!! TEST0695 tests JFRAC_TO_RFRAC and RFRAC_TO_JFRAC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 October 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: maxm = 10

  integer ( kind = 4 ) m
  real ( kind = 8 ) p(maxm)
  real ( kind = 8 ) q(maxm)
  real ( kind = 8 ) r(maxm)
  real ( kind = 8 ) s(maxm)
  integer ( kind = 4 ) seed
!
!  Generate the data, but force Q(M+1) to be 1.  
!  That will make it easier to see that the two operations are inverses
!  of each other.  JFRAC_TO_RFRAC is free to scale its output, and chooses
!  a scaling in which Q(M+1) is 1.
!
  seed = 123456789
  m = 6
  call r8vec_uniform_01 ( m, seed, p )
  call r8vec_uniform_01 ( m + 1, seed, q )

  q(1:m+1) = q(1:m+1) / q(m+1)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0695'
  write ( *, '(a)' ) '  RFRAC_TO_JFRAC converts a rational polynomial'
  write ( *, '(a)' ) '  fraction to a J fraction.'
  write ( *, '(a)' ) '  JFRAC_TO_RFRAC converts a J fraction'
  write ( *, '(a)' ) '  to a rational polynomial fraction.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The original rational polynomial coefficients:'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,5g14.6)' ) p(1:m)
  write ( *, '(2x,5g14.6)' ) q(1:m+1)
 
  call rfrac_to_jfrac ( m, p, q, r, s )
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The J fraction coefficients:'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,5g14.6)' ) r(1:m)
  write ( *, '(2x,5g14.6)' ) s(1:m)
 
  call jfrac_to_rfrac ( m, r, s, p, q )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The recovered rational polynomial:'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,5g14.6)' ) p(1:m)
  write ( *, '(2x,5g14.6)' ) q(1:m+1)

  return
end
subroutine test070 ( )

!*****************************************************************************80
!
!! TEST070 tests JOSEPHUS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST070'
  write ( *, '(a)' ) '  JOSEPHUS solves Josephus problems.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    N    M    K	 X'
  write ( *, '(a)' ) ' '

  m = 3
  n = 41
  k = 41
  call josephus ( n, m, k, x )
  write ( *, '(2x,4i5)' ) n, m, k, x

  m = -38
  n = 41
  k = 41
  call josephus ( n, m, k, x )

  write ( *, '(2x,4i5)' ) n, m, k, x

  m = 3
  n = 41
  k = 40
  call josephus ( n, m, k, x )

  write ( *, '(2x,4i5)' ) n, m, k, x

  m = 2
  n = 64
  k = 64
  call josephus ( n, m, k, x )

  write ( *, '(2x,4i5)' ) n, m, k, x

  m = 2
  n = 1000
  k = 1000
  call josephus ( n, m, k, x )

  write ( *, '(2x,4i5)' ) n, m, k, x

  return
end
subroutine test071 ( )

!*****************************************************************************80
!
!! TEST071 tests KSUB_NEXT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: k = 3

  integer ( kind = 4 ) a(k)
  logical more
  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ) rank

  a(1:k) = 0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST071'
  write ( *, '(a)' ) '  KSUB_NEXT generates all K subsets of an N set'
  write ( *, '(a)' ) '  in lexicographic order.'
  write ( *, '(a)' ) ' '

  more = .false.
  rank = 0
 
  do

    call ksub_next ( n, k, a, more )

    rank = rank + 1
    write ( *, '(2x,i4,4x,8i4)' ) rank, a(1:k)

    if ( .not. more ) then
      exit
    end if

  end do
 
  return
end
subroutine test072 ( )

!*****************************************************************************80
!
!! TEST072 tests KSUB_NEXT2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: k = 3

  integer ( kind = 4 ) a(k)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) in
  integer ( kind = 4 ) iout
  logical more
  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ) rank

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST072'
  write ( *, '(a)' ) '  KSUB_NEXT2 generates the next K subset of an'
  write ( *, '(a)' ) '  N set by the revolving door method.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Rank  Subset  Added  Removed'
  write ( *, '(a)' ) ' '
!
!  KSUB_NEXT2 doesn't have a good way of stopping.  
!  We will save the starting subset, and stop when the
!  new subset is the same as the starting one.
!
  in = 0
  iout = 0
  rank = 0
 
  call i4vec_indicator ( k, a )
 
  do
 
    rank = rank + 1
    write ( *, '(2x,i4,2x,3i2,3x,i2,7x,i2)' ) rank, a(1:k), in, iout
 
    call ksub_next2 ( n, k, a, in, iout )
 
    more = .false.

    do i = 1, k
      if ( a(i) /= i ) then
        more = .true.
      end if
    end do

    if ( .not. more ) then
      exit
    end if

  end do
 
  return
end
subroutine test073 ( )

!*****************************************************************************80
!
!! TEST073 tests KSUB_NEXT3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: k = 3

  integer ( kind = 4 ) a(k)
  integer ( kind = 4 ) in
  integer ( kind = 4 ) iout
  logical more
  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ) rank

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST073'
  write ( *, '(a)' ) '  KSUB_NEXT3 generates all K subsets of an N set'
  write ( *, '(a)' ) '  using the revolving door method.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Rank    Subset  Added Removed'
  write ( *, '(a)' ) ' '

  rank = 0
  more = .false.
 
  do

    call ksub_next3 ( n, k, a, more, in, iout )

    rank = rank + 1
    write ( *, '(2x,i4,2x,3i2,3x,i2,7x,i2)' ) rank, a(1:k), in, iout

    if ( .not. more ) then
      exit
    end if

  end do

  return
end
subroutine test074 ( )

!*****************************************************************************80
!
!! TEST074 tests KSUB_NEXT4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: k = 3

  integer ( kind = 4 ) a(k)
  logical done
  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ) rank

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST074'
  write ( *, '(a)' ) '  KSUB_NEXT4 generates K subsets of an N set.'
  write ( *, '(a,i8)' ) '  N = ', n
  write ( *, '(a,i8)' ) '  K=  ', k
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Rank    Subset'
  write ( *, '(a)' ) ' '

  done = .true.
  rank = 0
 
  do
 
    call ksub_next4 ( n, k, a, done )
 
    if ( done ) then
      exit
    end if

    rank = rank + 1
    write ( *, '(2x,i4,4x,3i4)' ) rank, a(1:k)

  end do
 
  return
end
subroutine test075 ( )

!*****************************************************************************80
!
!! TEST075 tests KSUB_RANDOM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: k = 3

  integer ( kind = 4 ) a(k)
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST075'
  write ( *, '(a)' ) '  KSUB_RANDOM generates a random K subset of an N set.'
  write ( *, '(a,i8)' ) '  Set size is N =    ', n
  write ( *, '(a,i8)' ) '  Subset size is K = ', k
  write ( *, '(a)' ) ' '

  seed = 123456789

  do i = 1, 5

    call ksub_random ( n, k, seed, a )

    write ( *, '(2x,8i3)' ) a(1:k)

  end do
 
  return
end
subroutine test076 ( )

!*****************************************************************************80
!
!! TEST076 tests KSUB_RANDOM2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: k = 3

  integer ( kind = 4 ) a(k)
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST076'
  write ( *, '(a)' ) '  KSUB_RANDOM2 generates a random K subset of an N set.'
  write ( *, '(a,i8)' ) '  Set size is N =    ', n
  write ( *, '(a,i8)' ) '  Subset size is K = ', k
  write ( *, '(a)' ) ' '
 
  seed = 123456789

  do i = 1, 5
    call ksub_random2 ( n, k, seed, a )
    write ( *, '(2x,8i3)' ) a(1:k)
  end do
 
  return
end
subroutine test077 ( )

!*****************************************************************************80
!
!! TEST077 tests KSUB_RANDOM3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: k = 3
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST077'
  write ( *, '(a)' ) '  KSUB_RANDOM3 generates a random K-subset of an N-set.'
  write ( *, '(a,i8)' ) '  Set size is N =    ', n
  write ( *, '(a,i8)' ) '  Subset size is K = ', k
  write ( *, '(a)' ) ' '

  seed = 123456789

  do i = 1, 10
    call ksub_random3 ( n, k, seed, a )
    write ( *, '(2x,15i3)' ) a(1:n)
  end do
 
  return
end
subroutine test0771 ( )

!*****************************************************************************80
!
!! TEST0771 tests KSUB_RANDOM4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: k = 3

  integer ( kind = 4 ) a(k)
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0771'
  write ( *, '(a)' ) '  KSUB_RANDOM4 generates a random K subset of an N set.'
  write ( *, '(a,i8)' ) '  Set size is N =    ', n
  write ( *, '(a,i8)' ) '  Subset size is K = ', k
  write ( *, '(a)' ) ' '

  seed = 123456789

  do i = 1, 5

    call ksub_random4 ( n, k, seed, a )

    write ( *, '(2x,8i3)' ) a(1:k)

  end do
 
  return
end
subroutine test07715 ( )

!*****************************************************************************80
!
!! TEST07715 tests KSUB_RANDOM5.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 June 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: k = 5

  integer ( kind = 4 ) a(k)
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: n = 52
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07715'
  write ( *, '(a)' ) '  KSUB_RANDOM5 generates a random K subset of an N set.'
  write ( *, '(a,i8)' ) '  Set size is N =    ', n
  write ( *, '(a,i8)' ) '  Subset size is K = ', k
  write ( *, '(a)' ) ' '

  seed = 123456789

  do i = 1, 5

    call ksub_random5 ( n, k, seed, a )

    write ( *, '(2x,8i3)' ) a(1:k)

  end do
 
  return
end
subroutine test0772 ( )

!*****************************************************************************80
!
!! TEST0772 tests KSUB_RANK.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ), parameter :: k = 3

  integer ( kind = 4 ), dimension ( k ) :: a = (/ 1, 3, 5 /)
  integer ( kind = 4 ) rank

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0772'
  write ( *, '(a)' ) '  KSUB_RANK: rank of a K subset of an N set.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  For N = ', n
  write ( *, '(a,i8)' ) '  and K = ', k
  write ( *, '(a)' ) '  the subset is:'
  write ( *, '(5i4)' ) a(1:k)
 
  call ksub_rank ( k, a, rank )
 
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The rank is ', rank
 
  return
end
subroutine test0773 ( )

!*****************************************************************************80
!
!! TEST0773 tests KSUB_UNRANK.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: k = 3

  integer ( kind = 4 ) a(k)
  integer ( kind = 4 ), parameter :: n = 5
  integer ( kind = 4 ) rank

  rank = 8
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0773'
  write ( *, '(a)' ) '  KSUB_UNRANK: find the K-subset of an N set'
  write ( *, '(a)' ) '  of a given rank.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  N is ', n
  write ( *, '(a,i8)' ) '  K is ', k
  write ( *, '(a,i8)' ) '  and the desired rank is ', rank
 
  call ksub_unrank ( k, rank, a )
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The subset of the given rank is:'
  write ( *, '(2x,5i4)' ) a(1:k)
 
  return
end
subroutine test078 ( )

!*****************************************************************************80
!
!! TEST078 tests MATRIX_PRODUCT_OPT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 6

  integer ( kind = 4 ) cost
  integer ( kind = 4 ) i
  integer ( kind = 4 ) order(n-1)
  integer ( kind = 4 ), dimension ( n+1) :: rank = (/ 4, 2, 3, 1, 2, 2, 3 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST078'
  write ( *, '(a)' ) '  MATRIX_PRODUCT_OPT seeks the optimal order'
  write ( *, '(a)' ) '  for a chain of matrix products.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Matrix ranks:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         I         R         C'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i8,2x,i8,2x,i8)' ) i, rank(i), rank(i+1)
  end do

  call matrix_product_opt ( n, rank, cost, order )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Optimal cost is ', cost

  call i4vec_print ( n-1, order, '  Ordering:' )

  return
end
subroutine test079 ( )

!*****************************************************************************80
!
!! TEST079 tests MOEBIUS_MATRIX.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 11

  integer ( kind = 4 ), dimension ( n, n ) :: ih = reshape (&
    (/ &
    0,0,1,1,0,0,0,0,0,0,0, &
    0,0,0,0,0,0,0,1,0,0,0, &
    0,1,0,0,0,0,0,0,0,0,0, &
    0,1,0,0,0,0,0,0,0,0,0, &
    0,0,0,1,0,0,0,0,0,0,0, &
    1,0,0,0,1,0,0,0,1,0,0, &
    0,0,0,0,0,1,0,0,0,1,1, &
    0,0,0,0,0,0,0,0,0,0,0, &
    0,0,1,1,0,0,0,0,0,0,0, &
    1,0,0,0,0,0,0,0,1,0,0, &
    0,0,0,0,0,0,0,0,1,0,0 /), (/ n, n /) )
  integer ( kind = 4 ) matrix(n,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST079'
  write ( *, '(a)' ) '  MOEBIUS_MATRIX computes the Moebius matrix.'
 
  call i4mat_print ( n, n, ih, '  The input matrix:' )

  call moebius_matrix ( n, ih, matrix )
 
  call i4mat_print ( n, n, matrix, '  The Moebius matrix:' )
 
  return
end
subroutine test0795 ( )

!*****************************************************************************80
!
!! TEST0795 tests MONOMIAL_COUNT and MONOMIAL_COUNTS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 October 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: degree_max = 9

  integer ( kind = 4 ) counts(0:degree_max)
  integer ( kind = 4 ) degree
  integer ( kind = 4 ) dim
  integer ( kind = 4 ) total

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0795'
  write ( *, '(a)' ) '  MONOMIAL_COUNT counts the number of monomials of'
  write ( *, '(a)' ) '  degrees 0 through DEGREE_MAX in a space of dimension DIM.'
  write ( *, '(a)' ) '  MONOMIAL_COUNTS provides individual counts.'
  do dim = 1, 6

    call monomial_counts ( degree_max, dim, counts )
    call monomial_count ( degree_max, dim, total )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  DIM = ', dim
    write ( *, '(a)' ) ' '
    do degree = 0, degree_max
      write ( *, '(2x,i8,2x,i8)' ) degree, counts(degree)
    end do
    write ( *, '(a)' ) ' '
    write ( *, '(2x,a8,2x,i8,2x,i8)' ) &
      '   Total', sum ( counts(0:degree_max) ), total

  end do

  return
end
subroutine test080 ( )

!*****************************************************************************80
!
!! TEST080 tests MORSE_THUE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 100

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) s(0:n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST080'
  write ( *, '(a)' ) '  MORSE_THUE computes the Morse-Thue numbers.'
  write ( *, '(a)' ) ' '

  do i = 0, n
    call morse_thue ( i, s(i) )
  end do

  do ilo = 0, n, 10
    ihi = min ( ilo + 9, n )
    write ( *, '(4x,40i1)' ) s(ilo:ihi)
  end do

  return
end
subroutine test081 ( )

!*****************************************************************************80
!
!! TEST081 tests MULTINOMIAL_COEF1 and MULTINOMIAL_COEF2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: maxfactor = 5

  integer ( kind = 4 ) factor(maxfactor)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ncomb1
  integer ( kind = 4 ) ncomb2
  integer ( kind = 4 ) nfactor

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST081'
  write ( *, '(a)' ) '  MULTINOMIAL_COEF1 computes multinomial'
  write ( *, '(a)' ) '  coefficients using the Gamma function;'
  write ( *, '(a)' ) '  MULTINOMIAL_COEF2 computes multinomial'
  write ( *, '(a)' ) '  coefficients directly.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Line 10 of the BINOMIAL table:'
  write ( *, '(a)' ) ' '

  n = 10
  nfactor = 2

  do i = 0, n

    factor(1) = i
    factor(2) = n - i

    call multinomial_coef1 ( nfactor, factor, ncomb1 )

    call multinomial_coef2 ( nfactor, factor, ncomb2 )

    write ( *, '(2x,i4,i4,3x,i5,i5)' ) factor(1), factor(2), ncomb1, ncomb2

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Level 5 of the TRINOMIAL coefficients:'

  n = 5
  nfactor = 3

  do i = 0, n

    factor(1) = i

    write ( *, '(a)' ) ' '

    do j = 0, n - factor(1)

      factor(2) = j
      factor(3) = n - factor(1) - factor(2)

      call multinomial_coef1 ( nfactor, factor, ncomb1 )

      call multinomial_coef2 ( nfactor, factor, ncomb2 )

      write ( *, '(2x,i4,i4,i4,3x,i5,i5)' ) factor(1), factor(2), factor(3), &
        ncomb1, ncomb2

    end do

  end do

  return
end
subroutine test0813 ( )

!*****************************************************************************80
!
!! TEST0813 tests MULTIPERM_ENUM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 July 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  integer ( kind = 4 ) counts(n)
  integer ( kind = 4 ), parameter :: i4_1 = 1
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) k
  integer ( kind = 4 ) number
  integer ( kind = 4 ) :: seed = 123456789
  integer ( kind = 4 ) test
  integer ( kind = 4 ) :: test_num = 5
  
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0813'
  write ( *, '(a)' ) '  MULTIPERM_ENUM enumerates multipermutations.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  N is the number of objects to be permuted.'
  write ( *, '(a)' ) '  K is the number of distinct types of objects.'
  write ( *, '(a)' ) '  COUNTS is the number of objects of each type.'
  write ( *, '(a)' ) '  NUMBER is the number of multipermutations.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Number       N       K       Counts(1:K)'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    k = i4_uniform ( i4_1, n, seed )

    call compnz_random ( n, k, seed, counts )

    call multiperm_enum ( n, k, counts, number )

    write ( *, '(2x,i6,2x,i6,2x,i6,5(2x,i4))' ) number, n, k, counts(1:k)

  end do
  
  return
end
subroutine test0815 ( )

!*****************************************************************************80
!
!! TEST0815 tests MULTIPERM_NEXT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 6

  integer ( kind = 4 ), dimension ( n ) :: a
  logical more
  integer ( kind = 4 ) tally
  
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0815'
  write ( *, '(a)' ) '  MULTIPERM_NEXT computes multipermutations in'
  write ( *, '(a)' ) '  lexical order.'
  write ( *, '(a)' ) ' '

  tally = 0
  a(1:n) = (/ 1, 2, 2, 3, 3, 3 /)
  more = .true.
 
  do while ( more )

    tally = tally + 1   

    write ( *, '(2x,i4,2x,6(2x,i2))' ) tally, a(1:n)
 
    call multiperm_next ( n, a, more )
            
  end do
  
  return
end
subroutine test082 ( )

!*****************************************************************************80
!
!! TEST082 tests NETWORK_FLOW_MAX.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: i4_2 = 2
  integer ( kind = 4 ), parameter :: nnode = 6
  integer ( kind = 4 ), parameter :: nedge = 20

  integer ( kind = 4 ) i
  integer ( kind = 4 ) icut(nnode)
  integer ( kind = 4 ), dimension ( 2, nedge ) :: icpflo = reshape ( &
    (/ 3,0,7,0,2,0,5,0,4,0,1,0,4,0,2,0,8,0,3,0, &
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /), (/ i4_2, nedge /) )
  integer ( kind = 4 ), dimension ( 2, nedge ) :: iendpt = reshape ( &
    (/ 1,2, 1,3, 2,3, 2,4, 2,5, 3,4, 3,5, 4,5, 4,6, 5,6, &
    2,1, 3,1, 3,2, 4,2, 5,2, 4,3, 5,3, 5,4, 6,4, 6,5 /), (/ i4_2, nedge /) )
  integer ( kind = 4 ) node_flow(nnode)
  integer ( kind = 4 ) :: sink = 6
  integer ( kind = 4 ) :: source = 1

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST082'
  write ( *, '(a)' ) '  NETWORK_FLOW_MAX finds the maximum flow on a network.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The source is node ', source
  write ( *, '(a,i8)' ) '  The sink is node   ', sink
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Endpoint array:'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,20i3)' ) iendpt(1,1:nedge)
  write ( *, '(2x,20i3)' ) iendpt(2,1:nedge)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Input edge capacity array:'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,20i3)' ) icpflo(1,1:nedge)
 
  call network_flow_max ( nnode, nedge, iendpt, icpflo, source, &
    sink, icut, node_flow )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Reordered endpoint array:'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,20i3)' ) ( iendpt(1,i), i = 1, nedge )
  write ( *, '(2x,20i3)' ) ( iendpt(2,i), i = 1, nedge )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Output edge capacity/flow array:'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,20i3)' ) ( icpflo(1,i), i = 1, nedge )
  write ( *, '(2x,20i3)' ) ( icpflo(2,i), i = 1, nedge )

  call i4vec_print ( nnode, icut, '  Minimal node cut vector:' )

  call i4vec_print ( nnode, node_flow, '  Nodal flow vector:' )

  return
end
subroutine test083 ( )

!*****************************************************************************80
!
!! TEST083 tests NIM_SUM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 November 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 32

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i1vec(n)
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i2vec(n)
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) i3vec(n)
  integer ( kind = 4 ), parameter :: ihi = 1000
  integer ( kind = 4 ), parameter :: ilo = 0
  integer ( kind = 4 ), parameter :: test_num = 5
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST083'
  write ( *, '(a)' ) '  NIM_SUM computes the Nim sum of two integers.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    I    J    Nim(I+J)'
  write ( *, '(a)' ) ' '

  seed = 123456789

  do i = 1, test_num

    i1 = i4_uniform ( ilo, ihi, seed )
    call ui4_to_ubvec ( i1, n, i1vec )

    i2 = i4_uniform ( ilo, ihi, seed )
    call ui4_to_ubvec ( i2, n, i2vec )

    call nim_sum ( i1, i2, i3 )
    call ui4_to_ubvec ( i3, n, i3vec )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  I1, I2, I3 in decimal:'
    write ( *, '(a)' ) ' '
    write ( *, '(i5)' ) i1
    write ( *, '(i5)' ) i2
    write ( *, '(i5)' ) i3
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  I1, I2, I3 in binary:'
    write ( *, '(a)' ) ' '
    call bvec_print ( n, i1vec, ' ' )
    call bvec_print ( n, i2vec, ' ' )
    call bvec_print ( n, i3vec, ' ' )

  end do

  return
end
subroutine test0835 ( )

!*****************************************************************************80
!
!! TEST0835 tests PADOVAN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 15

  integer ( kind = 4 ) i
  integer ( kind = 4 ) p(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0835'
  write ( *, '(a)' ) '  PADOVAN computes the Padovan numbers.'
  write ( *, '(a)' ) ' '

  call padovan ( n, p )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       N        P(N)'
  write ( *, '(a)' ) ' '

  do i = 0, n - 1
    write ( *, '(2x,i8,2x,i10)' ) i, p(i+1)
  end do

  return
end
subroutine test084 ( )

!*****************************************************************************80
!
!! TEST084 tests PELL_BASIC and PELL_NEXT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) q
  integer ( kind = 4 ) r
  integer ( kind = 4 ) x0
  integer ( kind = 4 ) x1
  integer ( kind = 4 ) y0
  integer ( kind = 4 ) y1

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST084'
  write ( *, '(a)' ) '  PELL_BASIC solves the basic Pell equation.'
  write ( *, '(a)' ) '  PELL_NEXT computes the "next" solution.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       D       X        Y         R'
  write ( *, '(a)' ) ' '

  do d = 2, 20

    call i4_sqrt ( d, q, r )

    if ( r /= 0 ) then

      call pell_basic ( d, x0, y0 )

      r = x0**2 - d * y0**2

      write ( *, '(2x,4i9)' ) d, x0, y0, r

      call pell_next ( d, x0, y0, x0, y0, x1, y1 )

      r = x1**2 - d * y1**2

      write ( *, '(2x,9x,3i9)' ) x1, y1, r

    end if

  end do

  return
end
subroutine test085 ( )

!*****************************************************************************80
!
!! TEST085 tests PENT_ENUM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  integer ( kind = 4 ) i
  integer ( kind = 4 ) pi

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST085'
  write ( *, '(a)' ) '  PENT_ENUM counts points in pentagons.'
  write ( *, '(a)' ) ' '

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  N    Pent(N)'
  write ( *, '(a)' ) ' '

  do i = 0, n
    call pent_enum ( i, pi )
    write ( *, '(2x,2i10)' ) i, pi
  end do

  return
end
subroutine test093 ( )

!*****************************************************************************80
!
!! TEST093 tests PERM_ASCEND.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 9

  integer ( kind = 4 ) length
  integer ( kind = 4 ), save, dimension ( n ) :: p = (/ 2,3,9,6,7,8,5,4,1 /)
  integer ( kind = 4 ) subseq(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST093'
  write ( *, '(a)' ) '  PERM_ASCEND determines the length of the longest'
  write ( *, '(a)' ) '  increasing subsequence in a permutation.'

  call perm_print ( n, p, '  The permutation:' )

  call perm_ascend ( n, p, length, subseq )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) &
    '  The length of the longest increasing subsequence is ', length

  call i4vec_print ( length, subseq, '  A longest increasing subsequence:' )

  return
end
subroutine test086 ( )

!*****************************************************************************80
!
!! TEST086 tests PERM_BREAK_COUNT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 6

  integer ( kind = 4 ) break_count
  integer ( kind = 4 ), save, dimension ( n ) :: p = (/ 4, 5, 2, 1, 6, 3 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST086'
  write ( *, '(a)' ) '  PERM_BREAK_COUNT counts the breaks in a permutation.'
 
  call perm_print ( n, p, '  The permutation:' )
 
  call perm_break_count ( n, p, break_count )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of breaks is ', break_count

  return
end
subroutine test087 ( )

!*****************************************************************************80
!
!! TEST087 tests PERM_CANON_TO_CYCLE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 6

  integer ( kind = 4 ), save, dimension ( n ) :: p1 = (/ 4, 5, 2, 1, 6, 3 /)
  integer ( kind = 4 ), dimension ( n ) :: p2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST087'
  write ( *, '(a)' ) '  PERM_CANON_TO_CYCLE converts a permutation from'
  write ( *, '(a)' ) '  canonical to cycle form.'
 
  call perm_print ( n, p1, '  The permutation in canonical form:' )
 
  call perm_canon_to_cycle ( n, p1, p2 )

  call perm_print ( n, p2, '  The permutation in cycle form:' )
 
  return
end
subroutine test088 ( )

!*****************************************************************************80
!
!! TEST088 tests PERM_CYCLE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 9

  integer ( kind = 4 ) iopt
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) ncycle
  integer ( kind = 4 ), save, dimension ( n ) :: p = (/ 2,3,9,6,7,8,5,4,1 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST088'
  write ( *, '(a)' ) '  PERM_CYCLE analyzes a permutation.'
 
  call perm_print ( n, p, '  The permutation:' )
 
  iopt = 1
  call perm_cycle ( n, iopt, p, isgn, ncycle )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  NCYCLE = ', ncycle
  write ( *, '(a,i8)' ) '  ISGN =   ', isgn

  call perm_print ( n, p, '  The permutation in cycle form:' )
 
  return
end
subroutine test089 ( )

!*****************************************************************************80
!
!! TEST089 tests PERM_CYCLE_TO_CANON.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 6

  integer ( kind = 4 ), save, dimension ( n ) :: p1 = (/ -6, 3, 1, -5, 4, -2 /)
  integer ( kind = 4 ), dimension ( n ) :: p2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST089'
  write ( *, '(a)' ) '  PERM_CYCLE_TO_CANON converts a permutation from'
  write ( *, '(a)' ) '  cycle to canonical form.'
 
  call perm_print ( n, p1, '  The permutation in cycle form:' )
 
  call perm_cycle_to_canon ( n, p1, p2 )

  call perm_print ( n, p2, '  The permutation in canonical form:' )
 
  return
end
subroutine test090 ( )

!*****************************************************************************80
!
!! TEST090 tests PERM_CYCLE_TO_INDEX and PERM_INDEX_TO_CYCLE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 9

  integer ( kind = 4 ), save, dimension ( n ) :: p1 = (/ 2,3,9,6,7,8,5,4,1 /)
  integer ( kind = 4 ), dimension ( n ) :: p2
  integer ( kind = 4 ), dimension ( n ) :: p3

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST090'
  write ( *, '(a)' ) '  PERM_CYCLE_TO_INDEX converts a permutation from'
  write ( *, '(a)' ) '  cycle to standard index form.'
  write ( *, '(a)' ) '  PERM_INDEX_TO_CYCLE converts a permutation from'
  write ( *, '(a)' ) '  standard index to cycle form.'
 
  call perm_print ( n, p1, '  The standard index form permutation:' )
 
  call perm_index_to_cycle ( n, p1, p2 )

  call perm_print ( n, p2, '  The permutation in cycle form:' )

  call perm_cycle_to_index ( n, p2, p3 )
 
  call perm_print ( n, p3, '  The standard index form permutation:' )
 
  return
end
subroutine test091 ( )

!*****************************************************************************80
!
!! TEST091 tests PERM_DISTANCE
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  integer ( kind = 4 ) k11
  integer ( kind = 4 ) k12
  integer ( kind = 4 ) k13
  integer ( kind = 4 ) k21
  integer ( kind = 4 ) k23
  integer ( kind = 4 ) p1(n)
  integer ( kind = 4 ) p2(n)
  integer ( kind = 4 ) p3(n)
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST091'
  write ( *, '(a)' ) '  PERM_DISTANCE computes the Ulam metric distance'
  write ( *, '(a)' ) '  between two permutations.'

  seed = 123456789

  call perm_random3 ( n, seed, p1 )
  call perm_print ( n, p1, '  Permutation P1' )
  call perm_random3 ( n, seed, p2 )
  call perm_print ( n, p2, '  Permutation P2' )
  call perm_random3 ( n, seed, p3 )
  call perm_print ( n, p3, '  Permutation P3' )

  call perm_distance ( n, p1, p1, k11 )
  call perm_distance ( n, p1, p2, k12 )
  call perm_distance ( n, p2, p1, k21 )
  call perm_distance ( n, p1, p3, k13 )
  call perm_distance ( n, p2, p3, k23 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  K(P1,P1) should be 0.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  K(P1,P1) = ', k11
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  K(P1,P2) should equal K(P2,P1).'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  K(P1,P2) = ', k12
  write ( *, '(a,i8)' ) '  K(P2,P1) = ', k21
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  K(P1,P3) <= K(P1,P2) + K(P2,P3).'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  K(P1,P3) = ', k13
  write ( *, '(a,i8)' ) '  K(P1,P2) = ', k12
  write ( *, '(a,i8)' ) '  K(P2,P3) = ', k23
  write ( *, '(a,i8)' ) '  K(P1,P2) + K(P2,P3) = ', k12 + k23

  return
end
subroutine test092 ( )

!*****************************************************************************80
!
!! TEST092 tests PERM_FIXED_ENUM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  integer ( kind = 4 ) fnm
  integer ( kind = 4 ) m

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST092'
  write ( *, '(a)' ) '  PERM_FIXED_ENUM enumerates the permutations'
  write ( *, '(a)' ) '  of N objects that leave M unchanged.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  For this test, N = ', n
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  M    F(N,M)'
  write ( *, '(a)' ) ' '

  do m = 0, n

    call perm_fixed_enum ( n, m, fnm )
    write ( *, '(2x,i3,2x,i8)' ) m, fnm

  end do

  return
end
subroutine test094 ( )

!*****************************************************************************80
!
!! TEST094 tests PERM_INVERSE;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 7

  integer ( kind = 4 ), dimension ( n ) :: p = (/ 4, 3, 5, 1, 7, 6, 2 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST094'
  write ( *, '(a)' ) '  PERM_INVERSE inverts a permutation in place;'
  write ( *, '(a)' ) ' '

  call perm_print ( n, p, '  The original permutation:' )
 
  call perm_inverse ( n, p )
 
  call perm_print ( n, p, '  The inverted permutation:' )
 
  return
end
subroutine test095 ( )

!*****************************************************************************80
!
!! TEST095 tests PERM_INVERSE2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 7

  integer ( kind = 4 ), dimension ( n ) :: p = (/ 4, 3, 5, 1, 7, 6, 2 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST095'
  write ( *, '(a)' ) '  PERM_INVERSE2 inverts a permutation in place.'

  call perm_print ( n, p, '  The original permutation:' )
 
  call perm_inverse2 ( n, p )
 
  call perm_print ( n, p, '  The inverted permutation:' )
 
  return
end
subroutine test0955 ( )

!*****************************************************************************80
!
!! TEST0955 tests PERM_INVERSE3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 7

  integer ( kind = 4 ), dimension ( n ) :: p = (/ 4, 3, 5, 1, 7, 6, 2 /)
  integer ( kind = 4 ), dimension ( n ) :: p_inv

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0955'
  write ( *, '(a)' ) '  PERM_INVERSE3 inverts a permutation.'

  call perm_print ( n, p, '  The original permutation:' )
 
  call perm_inverse3 ( n, p, p_inv )
 
  call perm_print ( n, p_inv, '  The inverted permutation:' )
 
  return
end
subroutine test096 ( )

!*****************************************************************************80
!
!! TEST096 tests PERM_LEX_NEXT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4

  logical more
  integer ( kind = 4 ) p(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST096'
  write ( *, '(a)' ) '  PERM_LEX_NEXT generates permutations in order.'
  write ( *, '(a)' ) ' '
  more = .false.
 
  do

    call perm_lex_next ( n, p, more )

    if ( .not. more ) then
      exit
    end if

    call perm_print ( n, p, ' ' )

  end do
 
  return
end
subroutine test097 ( )

!*****************************************************************************80
!
!! TEST097 tests PERM_LEX_NEXT and PERM_SIGN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4

  logical more
  integer ( kind = 4 ) p(n)
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) p_sign

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST097'
  write ( *, '(a)' ) '  PERM_LEX_NEXT generates permutations in order.'
  write ( *, '(a)' ) '  PERM_SIGN computes the sign of a permutation.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  RANK  SIGN  Permutation'
  write ( *, '(a)' ) ' '

  more = .false.
  rank = 0 

  do

    call perm_lex_next ( n, p, more )
    call perm_sign ( n, p, p_sign )

    if ( .not. more ) then
      exit
    end if

    write ( *, '(2x,i4,2x,i4,2x,10i4)' ) rank, p_sign, p(1:n)

    rank = rank + 1

  end do
 
  return
end
subroutine test098 ( )

!*****************************************************************************80
!
!! TEST098 tests PERM_MUL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  integer ( kind = 4 ) p1(n)
  integer ( kind = 4 ) p2(n)
  integer ( kind = 4 ) p3(n)
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST098'
  write ( *, '(a)' ) '  PERM_MUL multiplies two permutations.'
  write ( *, '(a)' ) ' '

  seed = 123456789

  call perm_random ( n, seed, p1 )
  call perm_random ( n, seed, p2 )

  call perm_print ( n, p1, '  Permutation P1:' )

  call perm_print ( n, p2, '  Permutation P2:' )

  call perm_mul ( n, p1, p2, p3 )

  call perm_print ( n, p3, '  Product permutation:' )

  return
end
subroutine test099 ( )

!*****************************************************************************80
!
!! TEST099 tests PERM_NEXT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4

  logical even
  logical more
  integer ( kind = 4 ) p(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST099'
  write ( *, '(a)' ) '  PERM_NEXT generates permutations.'
  write ( *, '(a)' ) ' '
  more = .false.
 
  do

    call perm_next ( n, p, more, even )

    call perm_print ( n, p, ' ' )

    if ( .not. more ) then
      exit
    end if
 
  end do

  return
end
subroutine test100 ( )

!*****************************************************************************80
!
!! TEST100 tests PERM_NEXT2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4

  logical done
  integer ( kind = 4 ) iactiv(n)
  integer ( kind = 4 ) idir(n)
  integer ( kind = 4 ) invers(n)
  integer ( kind = 4 ) p(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST100'
  write ( *, '(a)' ) '  PERM_NEXT2 generates permutations in order.'
  write ( *, '(a)' ) ' '
  done = .true.
 
  do

    call perm_next2 ( n, p, done, iactiv, idir, invers )
 
    if ( done ) then
      exit
    end if

    call perm_print ( n, p, ' ' )

  end do
 
  return
end
subroutine test101 ( )

!*****************************************************************************80
!
!! TEST101 tests PERM_NEXT3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4

  logical more
  integer ( kind = 4 ) p(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST101'
  write ( *, '(a)' ) '  PERM_NEXT3 generates permutations in order.'
  write ( *, '(a)' ) ' '
  more = .false.
 
  do

    call perm_next3 ( n, p, more )

    call perm_print ( n, p, ' ' )

    if ( .not. more ) then
      exit
    end if

  end do
 
  return
end
subroutine test102 ( )

!*****************************************************************************80
!
!! TEST102 tests PERM_RANDOM;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4

  integer ( kind = 4 ) i
  integer ( kind = 4 ) p(n)
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST102'
  write ( *, '(a)' ) '  PERM_RANDOM produces a random permutation;'
  write ( *, '(a,i8)' ) '  For this test, N = ', n
  write ( *, '(a)' ) ' '

  seed = 123456789

  do i = 1, 5
    call perm_random ( n, seed, p )
    call perm_print ( n, p, ' ' )
  end do
 
  return
end
subroutine test103 ( )

!*****************************************************************************80
!
!! TEST103 tests PERM_RANDOM2;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4

  integer ( kind = 4 ) i
  integer ( kind = 4 ) p(n)
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST103'
  write ( *, '(a)' ) '  PERM_RANDOM2 produces a random permutation of labels;'
  write ( *, '(a,i8)' ) '  For this test, N = ', n
  write ( *, '(a)' ) ' '

  seed = 123456789

  do i = 1, 5
    p(1:n) = (/ 101, 202, 303, 404 /)
    call perm_random2 ( n, seed, p )
    call perm_print ( n, p, ' ' )
  end do
 
  return
end
subroutine test104 ( )

!*****************************************************************************80
!
!! TEST104 tests PERM_RANDOM3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4

  integer ( kind = 4 ) i
  integer ( kind = 4 ) p(n)
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST104'
  write ( *, '(a)' ) '  PERM_RANDOM3 produces a random permutation.'
  write ( *, '(a,i8)' ) '  For this test, N = ', n
  write ( *, '(a)' ) ' '
 
  seed = 123456789

  do i = 1, 5
    call perm_random3 ( n, seed, p )
    call perm_print ( n, p, ' ' )
  end do
 
  return
end
subroutine test105 ( )

!*****************************************************************************80
!
!! TEST105 tests PERM_RANK.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4

  integer ( kind = 4 ), save, dimension ( n ) :: p = (/ 1, 4, 2, 3 /)
  integer ( kind = 4 ) rank

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST105'
  write ( *, '(a)' ) '  PERM_RANK ranks a permutation.'

  call perm_print ( n, p, '  The permutation:' )
 
  call perm_rank ( n, p, rank )
 
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The rank is:', rank
 
  return
end
subroutine test106 ( )

!*****************************************************************************80
!
!! TEST106 tests PERM_TO_EQUIV.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 9

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) jarray(n)
  integer ( kind = 4 ) npart
  integer ( kind = 4 ), save, dimension ( n ) :: p = (/ 2,3,9,6,7,8,5,4,1 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST106'
  write ( *, '(a)' ) '  PERM_TO_EQUIV returns the set partition'
  write ( *, '(a)' ) '  or equivalence classes determined by a'
  write ( *, '(a)' ) '  permutation.'

  call perm_print ( n, p, '  The input permutation:' )
 
  call perm_to_equiv ( n, p, npart, jarray, a )

  call equiv_print ( n, a, '  The partition:' )
 
  return
end
subroutine test107 ( )

!*****************************************************************************80
!
!! TEST107 tests PERM_TO_YTB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 7

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) lambda(n)
  integer ( kind = 4 ), dimension ( n ) ::  p = (/ 7, 2, 4, 1, 5, 3, 6 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST107'
  write ( *, '(a)' ) '  PERM_TO_YTB converts a permutation to a'
  write ( *, '(a)' ) '  Young table.'

  call perm_print ( n, p, '  The permutation:' )
 
  call perm_to_ytb ( n, p, lambda, a )

  call ytb_print ( n, a, '  The Young table:' )
 
  return
end
subroutine test108 ( )

!*****************************************************************************80
!
!! TEST108 tests PERM_UNRANK.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4

  integer ( kind = 4 ) p(n)
  integer ( kind = 4 ) rank

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST108'
  write ( *, '(a)' ) '  PERM_UNRANK, given a rank, computes the'
  write ( *, '(a)' ) '  corresponding permutation.'
  write ( *, '(a)' ) ' '
  rank = 6
  write ( *, '(a,i8)' ) '  The requested rank is ', rank
 
  call perm_unrank ( n, rank, p )
 
  call perm_print ( n, p, '  The permutation:' )
 
  return
end
subroutine test1085 ( )

!*****************************************************************************80
!
!! TEST1085 tests PERRIN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 15

  integer ( kind = 4 ) i
  integer ( kind = 4 ) p(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1085'
  write ( *, '(a)' ) '  PERRIN computes the Perrin numbers.'
  write ( *, '(a)' ) ' '

  call perrin ( n, p )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     N        P(N)'
  write ( *, '(a)' ) ' '

  do i = 0, n-1
    write ( *, '(2x,i8,i10)' ) i, p(i+1)
  end do

  return
end
subroutine test109 ( )

!*****************************************************************************80
!
!! TEST109 tests POWER_MOD.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST109'
  write ( *, '(a)' ) '  POWER_MOD computes the remainder of a power'
  write ( *, '(a)' ) '  of an integer modulo another integer.'

  a = 7
  n = 50
  m = 11

  call power_mod ( a, n, m, x )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  A = ', a
  write ( *, '(a,i8)' ) '  N = ', n
  write ( *, '(a,i8)' ) '  M = ', m
  write ( *, '(a,i8)' ) '  mod ( A**N, M ) = ', x

  a = 3
  n = 118
  m = 119

  call power_mod ( a, n, m, x )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  A = ', a
  write ( *, '(a,i8)' ) '  N = ', n
  write ( *, '(a,i8)' ) '  M = ', m
  write ( *, '(a,i8)' ) '  mod ( A**N, M ) = ', x

  return
end
subroutine test110 ( )

!*****************************************************************************80
!
!! TEST110 tests POWER_SERIES1;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) alpha
  real ( kind = 8 ) b(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST110'
  write ( *, '(a)' ) '  POWER_SERIES1 composes a power series;'

  alpha = 7.0D+00
 
  a(1) = 1.0D+00
  a(2:n) = 0.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Power series of G(x) = (1+F(x))**alpha'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  N = ', n
  write ( *, '(a,g14.6)' ) '  ALPHA = ', alpha
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Series for F(x):'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,5g14.6)' ) a(1:n)
 
  call power_series1 ( n, alpha, a, b )
 
  write ( *, '(a)' ) '  Series for G(x):'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,5g14.6)' ) b(1:n)
 
  return
end
subroutine test111 ( )

!*****************************************************************************80
!
!! TEST111 tests POWER_SERIES2;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) b(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST111'
  write ( *, '(a)' ) '  POWER_SERIES2 composes a power series;'
  write ( *, '(a)' ) '  Here we compute the power series of G(x) = exp(F(x))-1'
  write ( *, '(a,i8)' ) '  The number of terms is N = ', n

  a(1) = -4.0D+00
  a(2:n) = 0.0D+00
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Series for F(x):'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,5g14.6)' ) a(1:n)
 
  call power_series2 ( n, a, b )
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Series for G(x):'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,5g14.6)' ) b(1:n)
 
  return
end
subroutine test112 ( )

!*****************************************************************************80
!
!! TEST112 tests POWER_SERIES3;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) c(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST112'
  write ( *, '(a)' ) '  POWER_SERIES3 composes a power series;'
 
  a(1) = 1.0D+00
  a(2) = 1.0D+00
  a(3:n) = 0.0D+00
 
  b(1) = 1.0D+00
  b(2) = 1.0D+00
  b(3:n) = 0.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Power series of H(x) = G(F(x))'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of terms, N = ', n
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Series for F(x):'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,5g14.6)' ) a(1:n)
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Series for G(x):'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,5g14.6)' ) b(1:n)
 
  call power_series3 ( n, a, b, c )
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Series for H(x):'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,5g14.6)' ) c(1:n)
 
  return
end
subroutine test113 ( )

!*****************************************************************************80
!
!! TEST113 tests POWER_SERIES4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) c(n)
  integer ( kind = 4 ) i

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST113'
  write ( *, '(a)' ) '  POWER_SERIES4 composes a power series;'

  do i = 1, n
    a(i) = 1.0D+00 / real ( i, kind = 8 )
  end do

  b(1) = 1.0D+00
  b(2:n) = 0.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Power series of H(x) = G(1/F(x))'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  N = ', n
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Series for F(x):'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,5g14.6)' ) a(1:n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Series for G(x):'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,5g14.6)' ) b(1:n)
 
  call power_series4 ( n, a, b, c )
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Series for H(x):'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,5g14.6)' ) c(1:n)
 
  return
end
subroutine test114 ( )

!*****************************************************************************80
!
!! TEST114 tests PYTHAG_TRIPLE_NEXT
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) c
  integer ( kind = 4 ) d
  integer ( kind = 4 ) e
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST114'
  write ( *, '(a)' ) '  PYTHAG_TRIPLE_NEXT computes the "next"'
  write ( *, '(a)' ) '  Pythagorean triple.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   I   J   A   B   C  A^2+B^2   C^2'
  write ( *, '(a)' ) ' '

  i = 0
  j = 0

  do k = 0, 20
    call pythag_triple_next ( i, j, a, b, c )
    d = a**2 + b**2
    e = c**2
    write ( *, '(2x,5i4,2i8)' ) i, j, a, b, c, d, e
  end do

  return
end
subroutine test115 ( )

!*****************************************************************************80
!
!! TEST115 tests R8_AGM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) j
  real ( kind = 8 ) r8_agm
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST115'
  write ( *, '(a)' ) '  R8_AGM computes the arithmetic-geometric mean (AGM)'
  write ( *, '(a)' ) '  of two nonnegative real numbers.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    X        Y    R8_AGM(X,Y)'
  write ( *, '(a)' ) ' '

  seed = 123456789

  do i = 1, 10
    j = i4_uniform ( 1, 10, seed )
    x = real ( j, kind = 8 )
    j = i4_uniform ( 1, 10, seed )
    y = real ( j, kind = 8 )
    z = r8_agm ( x, y )
    write ( *, '(2x,f8.4,2x,f8.4,2x,f8.4)' ) x, y, z
  end do

  return
end
subroutine test030 ( )

!*****************************************************************************80
!
!! TEST030 tests R8_FALL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n
  real ( kind = 8 ) r8_fall
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST030'
  write ( *, '(a)' ) '  R8_FALL computes the falling factorial function.'
  write ( *, '(a)' ) '  [X]_N = X * (X-1) * (X-2) * ... * ( X-N+1).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    X          N  R8_FALL(X,N)'
  write ( *, '(a)' ) ' '

  x = 4.0D+00

  do n = -2, 5
 
    write ( *, '(2x,f10.6,i8,f10.6)' ) x, n, r8_fall ( x, n )
 
  end do
 
  return
end
subroutine test1245 ( )

!*****************************************************************************80
!
!! TEST1245 tests R8_RISE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n
  real ( kind = 8 ) r8_rise
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1245'
  write ( *, '(a)' ) '  R8_RISE computes the rising factorial function.'
  write ( *, '(a)' ) '  [X]_N = X * (X+1) * (X+2) * ... * ( X+N-1).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    X          N  R8_RISE(X,N)'
  write ( *, '(a)' ) ' '

  x = 4.0D+00

  do n = -2, 5
 
    write ( *, '(2x,f10.6,i8,f12.6)' ) x, n, r8_rise ( x, n )
 
  end do
 
  return
end
subroutine test116 ( )

!*****************************************************************************80
!
!! TEST116 tests R8_TO_CFRAC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 7

  integer ( kind = 4 ) a(0:n)
  real ( kind = 8 ) error
  integer ( kind = 4 ) i
  integer ( kind = 4 ) p(-1:n)
  integer ( kind = 4 ) q(-1:n)
  real ( kind = 8 ) r
  real ( kind = 8 ) r8_pi
  real ( kind = 8 ) temp

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST116'
  write ( *, '(a)' ) '  R8_TO_CFRAC converts a real number to a sequence'
  write ( *, '(a)' ) '  of continued fraction convergents.'

  r = 2.0D+00 * r8_pi ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Use the real number R = ', r

  call r8_to_cfrac ( r, n, a, p, q )

  write ( *, '(a)' ) ' '

  do i = 0, n
    temp = real ( p(i), kind = 8 ) / real ( q(i), kind = 8 )
    error = r - temp
    write ( *, '(2x,3i8,2g14.6)' ) a(i), p(i), q(i), temp, error
  end do

  return
end
subroutine test1163 ( )

!*****************************************************************************80
!
!! TEST1163 tests R8_TO_DEC and DEC_TO_R8.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) dec_digit
  integer ( kind = 4 ) i
  real ( kind = 8 ) r
  real ( kind = 8 ) r2
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1163'
  write ( *, '(a)' ) '  R8_TO_DEC converts a real number to a decimal;'
  write ( *, '(a)' ) '  DEC_TO_R8 converts a decimal to a real number.'

  dec_digit = 5

  write ( *, '(a)' ) ' '
  write ( *, '(a,i3)' ) '  The number of decimal digits is ', dec_digit

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     R   =>  A * 10^B  =>  R2'
  write ( *, '(a)' ) ' '

  do i = 1, 10
    r = r8_uniform_01 ( seed )
    r = 10.0D+00 * ( r - 0.25D+00 )
    call r8_to_dec ( r, dec_digit, a, b )
    call dec_to_r8 ( a, b, r2 )
    write ( *, '(2x,f10.6,2x,i8,2x,i8,2x,f10.6)' ) r, a, b, r2
  end do

  return
end
subroutine test1165 ( )

!*****************************************************************************80
!
!! TEST1165 tests R8_TO_RAT and RAT_TO_R8.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) i
  integer ( kind = 4 ) :: ndig = 4
  real ( kind = 8 ) r
  real ( kind = 8 ) r2
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1165'
  write ( *, '(a)' ) '  R8_TO_RAT converts a real number to a rational;'
  write ( *, '(a)' ) '  RAT_TO_R8 converts a rational to a real number.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i3)' ) '  The maximum number of digits allowed is ', ndig

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     R   =>  A / B  =>  R2'
  write ( *, '(a)' ) ' '

  do i = 1, 10
    r = r8_uniform_01 ( seed )
    r = 10.0D+00 * ( r - 0.25D+00 )
    call r8_to_rat ( r, ndig, a, b )
    call rat_to_r8 ( a, b, r2 )
    write ( *, '(2x,f10.6,i8,2x,i8,f10.6)' ) r, a, b, r2
  end do

  return
end
subroutine test125 ( )

!*****************************************************************************80
!
!! TEST125 tests R8MAT_DET.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n3 = 3
  integer ( kind = 4 ), parameter :: n4 = 4

  real ( kind = 8 ) a3(n3,n3)
  real ( kind = 8 ) a4(n4,n4)
  real ( kind = 8 ) det
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST125'
  write ( *, '(a)' ) '  R8MAT_DET: determinant of a real matrix.'
  write ( *, '(a)' ) ' '
 
  k = 0
  do i = 1, n3
    do j = 1, n3
      k = k+1
      a3(i,j) = real ( k, kind = 8 )
    end do
  end do
 
  call r8mat_print ( n3, n3, a3, '  The 123/456/789 matrix:' )

  call r8mat_det ( n3, a3, det )
 
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Determinant of the 123/456/789 matrix is ', det
 
  do i = 1, n4
    do j = 1, n4
      a4(i,j) = 1.0D+00 / real ( i + j, kind = 8 )
    end do
  end do
 
  call r8mat_print ( n4, n4, a4, '  The Hilbert matrix:' )

  call r8mat_det ( n4, a4, det )
 
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Determinant of the Hilbert matrix is ', det
 
  do i = 1, n3
    do j = 1, n3
      if ( i == j ) then
        a3(i,j) = 2.0D+00
      else if ( i == j+1 .or. i == j-1 ) then
        a3(i,j) = -1.0D+00
      else
        a3(i,j) = 0.0D+00
      end if
    end do
  end do
 
  call r8mat_print ( n3, n3, a3, '  The -1,2,-1 matrix:' )

  call r8mat_det ( n3, a3, det )
 
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Determinant of the -1,2,-1 matrix is ', det
 
  return
end
subroutine test126 ( )

!*****************************************************************************80
!
!! TEST126 tests R8MAT_PERM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 9

  real ( kind = 8 ) a(n,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ), save, dimension ( n ) :: p = (/ 2,3,9,6,7,8,5,4,1 /)
  integer ( kind = 4 ) j

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST126'
  write ( *, '(a)' ) '  R8MAT_PERM reorders a real matrix in place.'
  write ( *, '(a)' ) '  The rows and columns use the same permutation.'
 
  do i = 1, n
    do j = 1, n
      a(i,j) = real ( i * 10 + j, kind = 8 )
    end do
  end do
 
  call r8mat_print ( n, n, a, '  The original matrix' )
 
  call perm_print ( n, p, '  The row and column permutation:' )
 
  call r8mat_perm ( n, a, p )
 
  call r8mat_print ( n, n, a, '  The permuted matrix' )
 
  return
end
subroutine test127 ( )

!*****************************************************************************80
!
!! TEST127 tests R8MAT_PERM2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: m = 9
  integer ( kind = 4 ), parameter :: n = 7

  real ( kind = 8 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ), save, dimension ( m ) :: p = (/ 2, 3, 9, 6, 7, 8, 5, 4, 1 /)
  integer ( kind = 4 ), save, dimension ( n ) :: q = (/ 3, 4, 5, 6, 7, 1, 2 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST127'
  write ( *, '(a)' ) '  R8MAT_PERM2 reorders a real matrix in place.'
  write ( *, '(a)' ) '  Rows and columns use different permutations.'
 
  do i = 1, m
    do j = 1, n
      a(i,j) = real ( i * 10 + j, kind = 8 )
    end do
  end do
 
  call r8mat_print ( m, n, a, '  The original matrix' )
 
  call perm_print ( m, p, '  The row permutation:' )
 
  call perm_print ( n, q, '  The column permutation:' )

  call r8mat_perm2 ( m, n, a, p, q )
 
  call r8mat_print ( m, n, a, '  The permuted matrix' )
 
  return
end
subroutine test128 ( )

!*****************************************************************************80
!
!! TEST128 tests R8MAT_PERMANENT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ), allocatable, dimension ( :, : ) :: a
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  real ( kind = 8 ) perm

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST128'
  write ( *, '(a)' ) '  R8MAT_PERMANENT: the matrix permanent function.'
  write ( *, '(a)' ) '  We will analyze matrices with 0 diagonal and'
  write ( *, '(a)' ) '  1 on all offdiagonals.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Order	    Permanent.'
  write ( *, '(a)' ) ' '
 
  do n = 2, 12
 
    allocate ( a(1:n,1:n) )

    a(1:n,1:n) = 1.0D+00

    do i = 1, n
      a(i,i) = 0.0D+00
    end do
 
    call r8mat_permanent ( n, a, perm )
 
    write ( *, '(7x,i2,8x,g18.10)' ) n, perm

    deallocate ( a )
 
  end do
 
  return
end
subroutine test129 ( )

!*****************************************************************************80
!
!! TEST129 test R8POLY.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 6

  real ( kind = 8 ), dimension ( n ) :: a
  integer ( kind = 4 ) iopt
  integer ( kind = 4 ) test
  real ( kind = 8 ) val
  real ( kind = 8 ) x0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST129'
  write ( *, '(a)' ) '  R8POLY converts between power sum, factorial'
  write ( *, '(a)' ) '  and Taylor forms, and can evaluate a polynomial'
  write ( *, '(a)' ) ' '
 
  do test = 1, 6

    if ( test == 1 ) then
      iopt = -3
    else if ( test == 2 ) then
      iopt = -2
    else if ( test == 3 ) then
      iopt = -1
      x0 = 2.0D+00
    else if ( test == 4 ) then
      iopt = 0
      x0 = 2.0D+00
    else if ( test == 5 ) then
      iopt = 6
      x0 = 2.0D+00
    else if ( test == 6 ) then
      iopt = 6
      x0 = -2.0D+00
    end if

    a(1:n) = (/ 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 1.0D+00 /)

    if ( test == 1 ) then
      write ( *, '(a)' ) '  All calls have input A as follows'
      write ( *, '(2x,6f7.2)' ) a(1:n)
    end if
 
    call r8poly ( n, a, x0, iopt, val )
 
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Option IOPT = ', iopt
    if ( -1 <= iopt ) then
      write ( *, '(a,g14.6)' ) '  X0 = ', x0
    end if

    if ( iopt == -3 .or. iopt == -2 .or. 0 < iopt ) then
      write ( *, '(a)' ) '  Output array = '
      write ( *, '(2x,6f7.2)' ) a(1:n)
    end if

    if ( iopt == -1 .or. iopt == 0 ) then
      write ( *, '(a,g14.6)' ) '  Value = ', val
    end if
 
  end do

  return
end
subroutine test130 ( )

!*****************************************************************************80
!
!! TEST130 tests R8POLY_DIV.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a(0:10)
  real ( kind = 8 ) b(0:10)
  integer ( kind = 4 ) na
  integer ( kind = 4 ) nb
  integer ( kind = 4 ) nq
  integer ( kind = 4 ) nr
  real ( kind = 8 ) q(0:10)
  real ( kind = 8 ) r(0:10)
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST130'
  write ( *, '(a)' ) '  R8POLY_DIV computes the quotient and'
  write ( *, '(a)' ) '  remainder for polynomial division.'
  write ( *, '(a)' ) ' '
!
!  1: Divide X**3 + 2*X**2 - 5*X - 6  by X-2.  
!     Quotient is 3+4*X+X**2, remainder is 0.
!
!  2: Divide X**4 + 3*X**3 + 2*X**2 - 2  by  X**2 + X - 3.
!     Quotient is X**2 + 2*X + 3, remainder 8*X + 7.
!
  do test = 1, test_num

    if ( test == 1 ) then
      na = 3
      a(0:na) = (/ -6.0D+00, -5.0D+00, 2.0D+00, 1.0D+00 /)
      nb = 1
      b(0:nb) = (/ -2.0D+00, 1.0D+00 /)
    else if ( test == 2 ) then
      na = 4
      a(0:na) = (/ -2.0D+00, 5.0D+00, 2.0D+00, 3.0D+00, 1.0D+00 /)
      nb = 2
      b(0:nb) = (/ -3.0D+00, 1.0D+00, 1.0D+00 /)
    end if

    call r8poly_print ( na, a, '  The polynomial to be divided, A:' )
    call r8poly_print ( nb, b, '  The divisor polynomial, B:' )

    call r8poly_div ( na, a, nb, b, nq, q, nr, r )
 
    call r8poly_print ( nq, q, '  The quotient polynomial, Q:' )
    call r8poly_print ( nr, r, '  The remainder polynomial, R:' )

  end do

  return
end
subroutine test131 ( )

!*****************************************************************************80
!
!! TEST131 tests R8POLY_F2P and R8POLY_P2F.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4

  real ( kind = 8 ) a(n)

  call r8vec_indicator ( n, a )
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST131'
  write ( *, '(a)' ) '  R8POLY_P2F: power sum => factorial;'
  write ( *, '(a)' ) '  R8POLY_F2P: factorial => power sum.'

  call r8poly_print ( n-1, a, '  The power sum polynomial:' )
 
  call r8poly_p2f ( n, a )
 
  call r8vec_print ( n, a, '  The factorial polynomial coefficients:' )
 
  call r8poly_f2p ( n, a )
 
  call r8poly_print ( n-1, a, '  The recovered power sum polynomial:' )
 
  return
end
subroutine test132 ( )

!*****************************************************************************80
!
!! TEST132 tests R8POLY_FVAL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) val
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST132'
  write ( *, '(a)' ) '  R8POLY_FVAL evaluates a polynomial in factorial form.'

  call r8vec_indicator ( n, a )
 
  call r8vec_print ( n, a, '  The factorial polynomial coefficients:' )

  x = 2.0D+00

  call r8poly_fval ( n, a, x, val )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6,a,g14.6)' ) '  R8POLY (', x, ' ) = ', val
  write ( *, '(a)' ) '  The correct value is 11.'
 
  return
end
subroutine test133 ( )

!*****************************************************************************80
!
!! TEST133 tests R8POLY_MUL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: maxn = 5

  real ( kind = 8 ) a(0:maxn)
  real ( kind = 8 ) b(0:maxn)
  real ( kind = 8 ) c(0:maxn)
  integer ( kind = 4 ) na
  integer ( kind = 4 ) nb
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST133'
  write ( *, '(a)' ) '  R8POLY_MUL multiplies two polynomials.'
  write ( *, '(a)' ) ' '
!
!  1: Multiply (1+X) times (1-X).  Answer is 1-X**2.
!  2: Multiply (1+2*X+3*X**2) by (1-2*X). Answer is 1 + 0*X - X**2 - 6*X**3
!
  do test = 1, test_num

    if ( test == 1 ) then
      na = 1
      a(0:na) = (/ 1.0D+00, 1.0D+00 /)
      nb = 1
      b(0:nb) = (/ 1.0D+00, -1.0D+00 /)
    else if ( test == 2 ) then
      na = 2
      a(0:na) = (/ 1.0D+00, 2.0D+00, 3.0D+00 /)
      nb = 1
      b(0:nb) = (/ 1.0D+00, -2.0D+00 /)
    end if

    call r8poly_mul ( na, a, nb, b, c )

    call r8poly_print ( na, a, '  The factor A:' )

    call r8poly_print ( nb, b, '  The factor B:' )

    call r8poly_print ( na+nb, c, '  The product C = A*B:' )

  end do

  return
end
subroutine test134 ( )

!*****************************************************************************80
!
!! TEST134 tests R8POLY_N2P and R8POLY_P2N.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) a2(n)

  call r8vec_indicator ( n, a )

  a2(1:n) = 2.0D+00 * a(1:n)
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST134'
  write ( *, '(a)' ) '  R8POLY_N2P: Newton => power sum;'
  write ( *, '(a)' ) '  R8POLY_P2N: Power sum => Newton.'

  call r8poly_print ( n-1, a, '  The power sum polynomial:' )
 
  call r8poly_p2n ( n, a, a2 )
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Newton polynomial coefficients:'
  write ( *, '(a)' ) ' '
  write ( *, '(6f12.4)' ) a(1:n)
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Newton polynomial abscissas:'
  write ( *, '(a)' ) ' '
  write ( *, '(6f12.4)' ) a2(1:n)
 
  call r8poly_n2p ( n, a, a2 )
 
  call r8poly_print ( n-1, a, '  The recovered power sum polynomial:' )

  return
end
subroutine test135 ( )

!*****************************************************************************80
!
!! TEST135 tests R8POLY_NVAL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) a2(n)
  real ( kind = 8 ) val
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST135'
  write ( *, '(a)' ) '  R8POLY_NVAL evaluates a Newton polynomial.'

  call r8vec_indicator ( n, a )

  a2(1:n-1) = a(1:n-1) - 1.0D+00
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Newton polynomial coefficients:'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,6f12.4)' ) a(1:n)
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Newton polynomial abscissas:'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,6f12.4)' ) a2(1:n-1)
 
  x = 2.0D+00
 
  call r8poly_nval ( n, a, a2, x, val )
 
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6,a,g14.6)' ) '  R8POLY ( ', x,' ) = ', val
  write ( *, '(a)' ) '  The correct value is 11.'
 
  return
end
subroutine test136 ( )

!*****************************************************************************80
!
!! TEST136 tests R8POLY_P2T and R8POLY_T2P.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4

  real ( kind = 8 ) a(n)
  real ( kind = 8 ) x

  call r8vec_indicator ( n, a )

  x = 2.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST136'
  write ( *, '(a)' ) '  R8POLY_T2P: Taylor => Power sum;'
  write ( *, '(a)' ) '  R8POLY_P2T: Power sum => Taylor.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Taylor expansion point is X = ', x
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The Taylor coefficients:'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,6f12.4)' ) a(1:n)

  call r8poly_t2p ( n, a, x )

  call r8poly_print ( n-1, a, '  The power sum polynomial:' )

  call r8poly_p2t ( n, a, x )
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The recovered Taylor coefficients:'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,6f12.4)' ) a(1:n)
 
  return
end 
subroutine test137 ( )

!*****************************************************************************80
!
!! TEST137 tests R8POLY_POWER.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: lmax = 10

  real ( kind = 8 ) a(0:lmax)
  real ( kind = 8 ) b(0:10)
  integer ( kind = 4 ) na
  integer ( kind = 4 ) p

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST137'
  write ( *, '(a)' ) '  R8POLY_POWER takes a polynomial to a power.'
!
!  Cube (2-X).  Answer is 8-12*X+6*X**2-X**3.
!
  na = 1
  a(0:na) = (/ 2.0D+00, -1.0D+00 /)
  p = 3

  call r8poly_print ( na, a, '  The polynomial A:' )
 
  call r8poly_power ( na, a, p, b )
 
  call r8poly_print ( p*na, b, '  Raised to the power 3:' )
!
!  Square X+X**2
!
  na = 2
  a(0:na) = (/ 0.0D+00, 1.0D+00, 1.0D+00 /)
  p = 2

  call r8poly_print ( na, a, '  The polynomial A:' )
 
  call r8poly_power ( na, a, p, b )
 
  call r8poly_print ( p*na, b, '  Raised to the power 2:' )
 
  return
end
subroutine test138 ( )

!*****************************************************************************80
!
!! TEST138 tests R8POLY_PVAL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4

  integer ( kind = 4 ) i
  real ( kind = 8 ) a(0:n)
  real ( kind = 8 ) val
  real ( kind = 8 ) x

  do i = 0, n
    a(i) = real ( i + 1, kind = 8 )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST138'
  write ( *, '(a)' ) '  R8POLY_PVAL evaluates a polynomial'
  write ( *, '(a)' ) '  in power sum form.'

  call r8poly_print ( n, a, '  The polynomial to be evaluated:' )

  x = 2.0D+00
 
  call r8poly_pval ( n, a, x, val )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  At X = ', x
  write ( *, '(a,g14.6)' ) '  Computed polynomial value is ', val
  write ( *, '(a)' ) '  Correct value is 129.'
 
  return
end
subroutine test139 ( )

!*****************************************************************************80
!
!! TEST139 tests R8VEC_FRAC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  real ( kind = 8 ) a(n)
  real ( kind = 8 ), parameter :: ahi = 10.0D+00
  real ( kind = 8 ), parameter :: alo = 0.0D+00
  real ( kind = 8 ) afrac
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST139'
  write ( *, '(a)' ) '  R8VEC_FRAC: K-th smallest real vector entry;'

  seed = 123456789

  call r8vec_uniform ( n, alo, ahi, seed, a )

  call r8vec_print ( n, a, '  The real array to search: ' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Frac   R8VEC_FRAC'
  write ( *, '(a)' ) ' '

  do k = 1, n

    call r8vec_frac ( n, a, k, afrac )
    write ( *, '(2x,i4,2x,g14.6)' ) k, afrac

  end do

  return
end
subroutine test1395 ( )

!*****************************************************************************80
!
!! TEST1395 tests R8VEC_MIRROR_NEXT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3

  real ( kind = 8 ) a(n)
  logical done

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1395'
  write ( *, '(a)' ) '  R8VEC_MIRROR_NEXT generates all sign variations'
  write ( *, '(a)' ) '  of a real vector.'

  a(1:n) = (/ 1.0D+00, 2.0D+00, 3.0D+00 /)

  do

    call r8vec_print ( n, a, '  Next vector:' )

    call r8vec_mirror_next ( n, a, done )

    if ( done ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Done.'
      exit
    end if

  end do

  a(1:n) = (/ 1.0D+00, 0.0D+00, 3.0D+00 /)

  do

    call r8vec_print ( n, a, '  Next vector:' )

    call r8vec_mirror_next ( n, a, done )

    if ( done ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Done.'
      exit
    end if

  end do

  return
end
subroutine test117 ( )

!*****************************************************************************80
!
!! TEST117 tests RAT_ADD.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) abot
  integer ( kind = 4 ) atop
  integer ( kind = 4 ) bbot
  integer ( kind = 4 ) btop
  integer ( kind = 4 ) cbot
  integer ( kind = 4 ) ctop
  integer ( kind = 4 ) ierror
  character ( len = 22 ) string

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST117'
  write ( *, '(a)' ) '  RAT_ADD adds two rationals.'

  atop = 3
  abot = 4
  btop = 10
  bbot = 7

  call rat_add ( atop, abot, btop, bbot, ctop, cbot, ierror )

  write ( *, '(a)' ) ' '
  call rat_to_s_left ( atop, abot, string )
  write ( *, '(a)' ) '  A = ' // trim ( string )
  call rat_to_s_left ( btop, bbot, string )
  write ( *, '(a)' ) '  B = ' // trim ( string )
  call rat_to_s_left ( ctop, cbot, string )
  write ( *, '(a)' ) '  C = A + B = ' // trim ( string )
 
  return
end
subroutine test118 ( )

!*****************************************************************************80
!
!! TEST118 tests RAT_DIV.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) abot
  integer ( kind = 4 ) atop
  integer ( kind = 4 ) bbot
  integer ( kind = 4 ) btop
  integer ( kind = 4 ) cbot
  integer ( kind = 4 ) ctop
  integer ( kind = 4 ) ierror
  character ( len = 22 ) string

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST118'
  write ( *, '(a)' ) '  RAT_DIV divides two rationals.'

  atop = 3
  abot = 4
  btop = 10
  bbot = 7

  call rat_div ( atop, abot, btop, bbot, ctop, cbot, ierror )

  write ( *, '(a)' ) ' '
  call rat_to_s_left ( atop, abot, string )
  write ( *, '(a)' ) '  A = ' // trim ( string )
  call rat_to_s_left ( btop, bbot, string )
  write ( *, '(a)' ) '  B = ' // trim ( string )
  call rat_to_s_left ( ctop, cbot, string )
  write ( *, '(a)' ) '  C = A / B = ' // trim ( string )
 
  return
end
subroutine test119 ( )

!*****************************************************************************80
!
!! TEST119 tests RAT_FAREY.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: max_frac = 20

  integer ( kind = 4 ) a(max_frac)
  integer ( kind = 4 ) b(max_frac)
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) n
  integer ( kind = 4 ) num_frac

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST119'
  write ( *, '(a)' ) '  RAT_FAREY computes a row of the Farey fraction table.'

  do n = 1, 7

    call rat_farey ( n, max_frac, num_frac, a, b )
 
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Row ', n
    write ( *, '(a,i8)' ) '  Number of fractions: ', num_frac

    do ilo = 1, num_frac, 20
      ihi = min ( ilo+20-1, num_frac )
      write ( *, '(a)' ) ' '
      write ( *, '(2x,20i3)' ) a(ilo:ihi)
      write ( *, '(2x,20i3)' ) b(ilo:ihi)
    end do

  end do

  return
end
subroutine test120 ( )

!*****************************************************************************80
!
!! TEST120 tests RAT_FAREY2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: max_n = 4

  integer ( kind = 4 ) a(2**max_n+1)
  integer ( kind = 4 ) b(2**max_n+1)
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) n

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST120'
  write ( *, '(a)' ) '  RAT_FAREY2 computes a row of the Farey fraction table.'

  do n = 0, max_n

    call rat_farey2 ( n, a, b )
 
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Row ', n+1

    do ilo = 1, 2**n+1, 20
      ihi = min ( ilo+20-1, 2**n+1 )
      write ( *, '(a)' ) ' '
      write ( *, '(2x,20i3)' ) a(ilo:ihi)
      write ( *, '(2x,20i3)' ) b(ilo:ihi)
    end do

  end do

  return
end
subroutine test121 ( )

!*****************************************************************************80
!
!! TEST121 tests RAT_MUL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) abot
  integer ( kind = 4 ) atop
  integer ( kind = 4 ) bbot
  integer ( kind = 4 ) btop
  integer ( kind = 4 ) cbot
  integer ( kind = 4 ) ctop
  integer ( kind = 4 ) ierror
  character ( len = 22 ) string

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST121'
  write ( *, '(a)' ) '  RAT_MUL multiplies two rationals.'

  atop = 3
  abot = 4
  btop = 10
  bbot = 7

  call rat_mul ( atop, abot, btop, bbot, ctop, cbot, ierror )

  write ( *, '(a)' ) ' '
  call rat_to_s_left ( atop, abot, string )
  write ( *, '(a)' ) '  A = ' // trim ( string )
  call rat_to_s_left ( btop, bbot, string )
  write ( *, '(a)' ) '  B = ' // trim ( string )
  call rat_to_s_left ( ctop, cbot, string )
  write ( *, '(a)' ) '  C = A * B = ' // trim ( string )
 
  return
end
subroutine test1215 ( )

!*****************************************************************************80
!
!! TEST1215 tests RAT_WIDTH.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_test = 17

  integer ( kind = 4 ) a
  integer ( kind = 4 ), dimension ( n_test ) :: a_test = (/ &
    1000, 1000, 1000, 1000, 1000, 1, -1, -10, -100, -1000, &
    1, 10, 100, 1000, 10000, 17, 4000000 /)
  integer ( kind = 4 ) b
  integer ( kind = 4 ), dimension ( n_test ) :: b_test = (/ &
    3, 40, 500, 6000, 70000, 1, 200, 200, 200, 200, &
   -200, -200, -200, -200, -200, 3000, 4000000 /)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) rat_width
  integer ( kind = 4 ) width

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1215'
  write ( *, '(a)' ) '  RAT_WIDTH determines the "width" of a rational.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     Top    Bottom  Width'
  write ( *, '(a)' ) ' '

  do i = 1, n_test
    a = a_test(i)
    b = b_test(i)
    width = rat_width ( a, b )
    write ( *, '(2x,3i8)' ) a, b, width
  end do

  return
end
subroutine test122 ( )

!*****************************************************************************80
!
!! TEST122 tests RATMAT_DET.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 6
  integer ( kind = 4 ) a(0:n,n+1)
  integer ( kind = 4 ) b(0:n,n+1)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST122'
  write ( *, '(a)' ) '  RAT_SUM_FORMULA computes the coefficients for the'
  write ( *, '(a)' ) '  formulas for the sums of powers of integers.'
  
  call rat_sum_formula ( n, a, b )

  call ratmat_print ( n+1, n+1, a, b, '  Power Sum Coefficients:' )

  return
end
subroutine test123 ( )

!*****************************************************************************80
!
!! TEST123 tests RATMAT_DET.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n3 = 3

  integer ( kind = 4 ) a3(n3,n3)
  integer ( kind = 4 ) b3(n3,n3)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) idbot
  integer ( kind = 4 ) idtop
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST123'
  write ( *, '(a)' ) '  RATMAT_DET: determinant of a rational matrix.'
  write ( *, '(a)' ) ' '
 
  k = 0
  do i = 1, n3
    do j = 1, n3
      k = k + 1
      a3(i,j) = k
    end do
  end do

  b3(1:n3,1:n3) = 1
 
  call ratmat_print ( n3, n3, a3, b3, '  The 123/456/789 matrix:' )

  call ratmat_det ( n3, a3, b3, idtop, idbot, ierror )
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Determinant of the 123/456/789 matrix:'
  write ( *, '(2x,i8,a,i8)' ) idtop, ' / ', idbot
 
  do i = 1, n3
    do j = 1, n3
      a3(i,j) = 1
      b3(i,j) = i + j
    end do
  end do
 
  call ratmat_print ( n3, n3, a3, b3, '  The Hilbert matrix:' )

  call ratmat_det ( n3, a3, b3, idtop, idbot, ierror )
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Determinant of the Hilbert matrix:'
  write ( *, '(2x,i8,a,i8)' ) idtop, ' / ', idbot
 
  do i = 1, n3
    do j = 1, n3
      if ( i == j ) then
        a3(i,j) = 2
      else if ( i == j+1 .or. i == j-1 ) then
        a3(i,j) = -1
      else
        a3(i,j) = 0
      end if
      b3(i,j) = 1
    end do
  end do
 
  call ratmat_print ( n3, n3, a3, b3, '  The -1 2 -1 matrix:' )

  call ratmat_det ( n3, a3, b3, idtop, idbot, ierror )
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Determinant of the -1,2,-1 matrix:'
  write ( *, '(2x,i8,a,i8)' ) idtop, ' / ', idbot
 
  return
end
subroutine test124 ( )

!*****************************************************************************80
!
!! TEST124 tests REGRO_NEXT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4

  logical done
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) v(n)
  integer ( kind = 4 ) vmax(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST124'
  write ( *, '(a)' ) '  REGRO_NEXT generates all restricted growth '
  write ( *, '(a)' ) '  functions.'
  write ( *, '(a)' ) ' '

  rank = 0

  done = .true.
 
  do

    call regro_next ( n, v, vmax, done )

    if ( done ) then
      exit
    end if

    rank = rank + 1
    write ( *, '(2x,5i3)' ) rank, v(1:n)

  end do
 
  return
end
subroutine test140 ( )

!*****************************************************************************80
!
!! TEST140 tests SCHROEDER.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  integer ( kind = 4 ) i
  integer ( kind = 4 ) s(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST140'
  write ( *, '(a)' ) '  SCHROEDER computes the Schroeder numbers.'
  write ( *, '(a)' ) ' '

  call schroeder ( n, s )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       N        S(N)'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i8,i10)' ) i, s(i)
  end do

  return
end
subroutine test141 ( )

!*****************************************************************************80
!
!! TEST141 tests SORT_HEAP_EXTERNAL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 20

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST141'
  write ( *, '(a)' ) '  SORT_HEAP_EXTERNAL sorts objects externally.'
  write ( *, '(a)' ) ' '

  indx = 0
  i = 0
  j = 0
  isgn = 0
  seed = 123456789

  call i4vec_uniform ( n, 1, n, seed, a )
 
  call i4vec_print ( n, a, '  Unsorted array:' )
 
  do

    call sort_heap_external ( n, indx, i, j, isgn )
 
    if ( indx < 0 ) then
      isgn = 1
      if ( a(i) <= a(j) ) then
        isgn = -1
      end if
    else if ( 0 < indx ) then
      call i4_swap ( a(i), a(j) )
    else
      exit
    end if

  end do

  call i4vec_print ( n, a, '  Sorted array:' )
 
  return
end
subroutine test142 ( )

!*****************************************************************************80
!
!! TEST142 tests SUBSET_BY_SIZE_NEXT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  integer ( kind = 4 ) a(n)
  logical more
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) size

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST142'
  write ( *, '(a)' ) '  SUBSET_BY_SIZE_NEXT generates all subsets of an N set.'
  write ( *, '(a)' ) ' '

  more = .false.
  rank = 0

  do

    call subset_by_size_next ( n, a, size, more )

    rank = rank + 1

    if ( 0 < size ) then
      write ( *, '(2x,i4,4x,5i2)' ) rank, a(1:size)
    else
      write ( *, '(2x,i4,4x,a)' ) rank, 'The empty set'
    end if

    if ( .not. more ) then
      exit
    end if

  end do

  return
end
subroutine test143 ( )

!*****************************************************************************80
!
!! TEST143 tests SUBSET_LEX_NEXT without size restrictions.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) k
  logical ltest
  integer ( kind = 4 ) ndim
  integer ( kind = 4 ) rank

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST143'
  write ( *, '(a)' ) '  SUBSET_LEX_NEXT generates all subsets of an N set.'
  write ( *, '(a)' ) ' '

  ndim = n
  k = 0
  rank = 0
  ltest = .false.
 
  do
 
    call subset_lex_next ( n, ltest, ndim, k, a )
 
    rank = rank + 1

    if ( 0 < k ) then
      write ( *, '(2x,i4,4x,6i2)' ) rank, a(1:k)
    else
      write ( *, '(2x,i4,4x,a)' ) rank, 'The empty set.'
    end if
 
    if ( k == 0 ) then
      exit
    end if

  end do
 
  return
end
subroutine test1435 ( )

!*****************************************************************************80
!
!! TEST1435 tests SUBSET_LEX_NEXT with size restrictions.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: ndim = 3

  integer ( kind = 4 ) a(ndim)
  integer ( kind = 4 ) k
  logical ltest
  integer ( kind = 4 ) n

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1435'
  write ( *, '(a)' ) '  SUBSET_LEX_NEXT generates all subsets of an N set.'
  write ( *, '(a)' ) '  The user can impose a restriction on the'
  write ( *, '(a)' ) '  maximum size of the subsets.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Here, we require the subsets to be no larger'
  write ( *, '(a,i8)' ) '  than ', ndim

  n = 5
  k = 0
 
  do
 
    ltest = ( k == ndim )

    call subset_lex_next ( n, ltest, ndim, k, a )
 
    if ( 0 < k ) then
      write ( *, '(2x,6i2)' ) a(1:k)
    else
      write ( *, '(a)' ) '  The empty set.'
    end if
 
    if ( k == 0 ) then
      exit
    end if

  end do
 
  return
end
subroutine test144 ( )

!*****************************************************************************80
!
!! TEST144 tests SUBSET_GRAY_NEXT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) iadd
  logical more
  integer ( kind = 4 ) ncard
  integer ( kind = 4 ) rank

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST144'
  write ( *, '(a)' ) '  SUBSET_GRAY_NEXT generates all subsets of an N set.'
  write ( *, '(a)' ) '  using the Gray code ordering:'
  write ( *, '(a)' ) '  0 0 1 0 1 means the subset contains 3 and 5.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Gray code'
  write ( *, '(a)' ) ' '
 
  rank = 0
  more = .false.
 
  do
 
    call subset_gray_next ( n, a, more, ncard, iadd )

    rank = rank + 1 
    write ( *, '(2x,i4,4x,5i2)' ) rank, a(1:n)

    if ( .not. more ) then
      exit
    end if

  end do

  return
end
subroutine test145 ( )

!*****************************************************************************80
!
!! TEST145 tests SUBSET_RANDOM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST145'
  write ( *, '(a)' ) '  SUBSET_RANDOM picks a subset at random.'
  write ( *, '(a,i8)' ) '  The number of elements in the main set is ', n
  write ( *, '(a)' ) ' '

  seed = 123456789

  do i = 1, 5
    call subset_random ( n, seed, a )
    write ( *, '(2x,40i2)' ) a(1:n)
  end do
 
  return
end
subroutine test146 ( )

!*****************************************************************************80
!
!! TEST146 tests SUBSET_GRAY_RANK.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  integer ( kind = 4 ), dimension ( n ) :: a = (/ 1, 0, 1, 1, 0 /)
  integer ( kind = 4 ) rank

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST146'
  write ( *, '(a)' ) '  SUBSET_GRAY_RANK returns rank of a subset of an N set'
  write ( *, '(a)' ) '  using the Gray code ordering.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  For N = ', n
  write ( *, '(a)' ) '  the subset is:'
  write ( *, '(2x,5i2)' ) a(1:n)
 
  call subset_gray_rank ( n, a, rank )
 
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The rank is ', rank
 
  return
end
subroutine test147 ( )

!*****************************************************************************80
!
!! TEST147 tests SUBSET_GRAY_UNRANK.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 5

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) rank

  rank = 8
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST147'
  write ( *, '(a)' ) '  SUBSET_GRAY_UNRANK finds the subset of an N set'
  write ( *, '(a)' ) '  of a given rank under the Gray code ordering.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  N is ', n
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Rank   Subset'
  write ( *, '(a)' ) ' '

  do rank = 1, 10
 
    call subset_gray_unrank ( rank, n, a )

    write ( *, '(2x,i4,4x,5i2)' ) rank, a(1:n)

  end do
 
  return
end
subroutine test1475 ( )

!*****************************************************************************80
!
!! TEST1475 tests SUBCOMP_NEXT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 July 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: k = 3

  integer ( kind = 4 ) a(k)
  integer ( kind = 4 ) count
  integer ( kind = 4 ) ih
  integer ( kind = 4 ) it
  logical more
  integer ( kind = 4 ) :: n = 6

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1475'
  write ( *, '(a)' ) '  SUBCOMP_NEXT generates subcompositions.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Seek all subcompositions of N = ', n
  write ( *, '(a,i8,a)' ) '  using K = ', k, ' parts.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     #   Sum'
  write ( *, '(a)' ) ' '

  more = .false.
  count = 0

  do

    call subcomp_next ( n, k, a, more, ih, it )

    count = count + 1
    write ( *, '(2x,i4,2x,i4,2x,8i4)' ) count, sum ( a(1:k) ), a(1:k)

    if ( .not. more )  then
      exit
    end if

  end do
 
  return
end
subroutine test1476 ( )

!*****************************************************************************80
!
!! TEST1476 tests SUBCOMPNZ_NEXT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 December 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: k = 3

  integer ( kind = 4 ) a(k)
  integer ( kind = 4 ) count
  logical more
  integer ( kind = 4 ) :: n = 6

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1476'
  write ( *, '(a)' ) '  SUBCOMPNZ_NEXT generates subcompositions'
  write ( *, '(a)' ) '  using nonzero parts.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Seek all subcompositions of N = ', n
  write ( *, '(a,i8,a)' ) '  using K = ', k, ' nonzero parts.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     #   Sum'
  write ( *, '(a)' ) ' '

  more = .false.
  count = 0

  do

    call subcompnz_next ( n, k, a, more )

    count = count + 1
    write ( *, '(2x,i4,2x,i4,2x,8i4)' ) count, sum ( a(1:k) ), a(1:k)

    if ( .not. more )  then
      exit
    end if

  end do
 
  return
end
subroutine test1477 ( )

!*****************************************************************************80
!
!! TEST1477 tests SUBCOMPNZ2_NEXT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 December 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: k = 3

  integer ( kind = 4 ) a(k)
  integer ( kind = 4 ) count
  logical more
  integer ( kind = 4 ) n
  integer ( kind = 4 ) :: n_hi = 7
  integer ( kind = 4 ) :: n_lo = 5

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1477'
  write ( *, '(a)' ) '  SUBCOMPNZ2_NEXT generates subcompositions'
  write ( *, '(a)' ) '  using nonzero parts.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Seek all subcompositions of N'
  write ( *, '(a,i8,a)' ) '  using K = ', k, ' nonzero parts.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a,i8)' ) '  N ranges from ', n_lo, ' to ', n_hi
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     #     N'
  write ( *, '(a)' ) ' '

  more = .false.
  count = 0

  do

    call subcompnz2_next ( n_lo, n_hi, k, a, more )

    count = count + 1
    n = sum ( a(1:k) )
    write ( *, '(2x,i4,2x,i4,2x,8i4)' ) count, n, a(1:k)

    if ( .not. more )  then
      exit
    end if

  end do
 
  return
end
subroutine test1478 ( )

!*****************************************************************************80
!
!! TEST1478 tests SUBTRIANGLE_NEXT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j3
  logical more 
  integer ( kind = 4 ) n
  integer ( kind = 4 ) rank

  n = 4
  rank = 0

  more = .false.
  i1 = 0
  j1 = 0
  i2 = 0
  j2 = 0
  i3 = 0
  j3 = 0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1478'
  write ( *, '(a)' ) '  SUBTRIANGLE_NEXT generates the indices of subtriangles'
  write ( *, '(a)' ) '  in a triangle whose edges were divided into N subedges.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  For this test, N = ', n
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Rank    I1  J1    I2  J2    I3  J3'
  write ( *, '(a)' ) ' '

  do

    call subtriangle_next ( n, more, i1, j1, i2, j2, i3, j3 )

    rank = rank + 1

    write ( *, '(2x,i4,4x,i2,2x,i2,4x,i2,2x,i2,4x,i2,2x,i2)' ) &
      rank, i1, j1, i2, j2, i3, j3

    if ( .not. more ) then
      exit
    end if

  end do

  return
end
subroutine test148 ( )

!*****************************************************************************80
!
!! TEST148 tests THUE_BINARY_NEXT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 100

  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) thue(n_max)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST148'
  write ( *, '(a)' ) '  THUE_BINARY_NEXT returns the next Thue binary sequence.'
  write ( *, '(a)' ) ' '

  n = 1
  thue(1) = 0
  write ( *, '(2x,i4,4x,80i1)' ) n, thue(1:n)

  do i = 1, 6
    call thue_binary_next ( n, thue )
    write ( *, '(2x,i4,4x,80i1)' ) n, thue(1:n)
  end do

  return
end
subroutine test149 ( )

!*****************************************************************************80
!
!! TEST149 tests THUE_TERNARY_NEXT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 100

  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
  integer ( kind = 4 ) thue(n_max)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST149'
  write ( *, '(a)' ) '  THUE_TERNARY_NEXT returns the next '
  write ( *, '(a)' ) '  Thue ternary sequence.'
  write ( *, '(a)' ) ' '

  n = 1
  thue(1) = 1
  write ( *, '(2x,i4,4x,80i1)' ) n, thue(1:n)

  do i = 1, 5
    call thue_ternary_next ( n, thue )
    write ( *, '(2x,i4,4x,80i1)' ) n, thue(1:n)
  end do

  return
end
subroutine test150 ( )

!*****************************************************************************80
!
!! TEST150 tests TUPLE_NEXT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 2

  integer ( kind = 4 ), parameter :: m1 = 2
  integer ( kind = 4 ), parameter :: m2 = 4
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST150'
  write ( *, '(a)' ) '  TUPLE_NEXT returns the next "tuple", that is,'
  write ( *, '(a)' ) '  a vector of N integers, each between M1 and M2.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  M1 = ', m1
  write ( *, '(a,i8)' ) '  M2 = ', m2
  write ( *, '(a,i8)' ) '  N =  ', n
  write ( *, '(a)' ) ' '

  rank = 0

  do

    call tuple_next ( m1, m2, n, rank, x )

    if ( rank == 0 ) then
      exit
    end if

    write ( *, '(2x,i4,2x,10i3)' ) rank, x(1:n)

  end do

  return
end
subroutine test151 ( )

!*****************************************************************************80
!
!! TEST151 tests TUPLE_NEXT_FAST.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 2

  integer ( kind = 4 ), parameter :: m = 3
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST151'
  write ( *, '(a)' ) '  TUPLE_NEXT_FAST returns the next "tuple", that is,'
  write ( *, '(a)' ) '  a vector of N integers, each between 1 and M.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  M = ', m
  write ( *, '(a,i8)' ) '  N = ', n
  write ( *, '(a)' ) ' '
!
!  Initialize.
!
  rank = -1
  call tuple_next_fast ( m, n, rank, x )

  do rank = 0, (m**n)-1

    call tuple_next_fast ( m, n, rank, x )

    write ( *, '(2x,i4,2x,10i3)' ) rank, x(1:n)

  end do

  return
end
subroutine test152 ( )

!*****************************************************************************80
!
!! TEST152 tests TUPLE_NEXT_GE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3

  integer ( kind = 4 ), parameter :: m = 3
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST152'
  write ( *, '(a)' ) '  TUPLE_NEXT_GE returns the next "tuple", that is,'
  write ( *, '(a)' ) '  a vector of N integers, each between 1 and M,'
  write ( *, '(a)' ) '  with the constraint that the entries be'
  write ( *, '(a)' ) '  nondecreasing.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  M = ', m
  write ( *, '(a,i8)' ) '  N = ', n
  write ( *, '(a)' ) ' '

  rank = 0

  do

    call tuple_next_ge ( m, n, rank, x )

    if ( rank == 0 ) then
      exit
    end if

    write ( *, '(2x,i4,2x,10i3)' ) rank, x(1:n)

  end do

  return
end
subroutine test153 ( )

!*****************************************************************************80
!
!! TEST153 tests TUPLE_NEXT2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3

  integer ( kind = 4 ), parameter, dimension ( n ) :: xmin = (/ 2, 3, 8 /)
  integer ( kind = 4 ), parameter, dimension ( n ) :: xmax = (/ 4, 3, 5 /)
  integer ( kind = 4 ) rank
  integer ( kind = 4 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST153'
  write ( *, '(a)' ) '  TUPLE_NEXT2 returns the next "tuple", that is,'
  write ( *, '(a)' ) '  a vector of N integers.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  N = ', n
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The minimum tuple is '
  write ( *, '(2x,5i8)' ) xmin(1:n)
  write ( *, '(a)' ) '  The maximum tuple is '
  write ( *, '(2x,5i8)' ) xmax(1:n)
  write ( *, '(a)' ) ' '

  rank = 0

  do

    call tuple_next2 ( n, xmin, xmax, rank, x )

    if ( rank == 0 ) then
      exit
    end if

    write ( *, '(2x,i4,2x,10i3)' ) rank, x(1:n)

  end do

  return
end
subroutine test1531 ( )

!*****************************************************************************80
!
!! TEST1531 tests UBVEC_ADD;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 November 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  integer ( kind = 4 ) bvec1(n)
  integer ( kind = 4 ) bvec2(n)
  integer ( kind = 4 ) bvec3(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_uniform
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) :: seed = 123456789
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 10

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1531'
  write ( *, '(a)' ) '  UBVEC_ADD adds unsigned binary vectors'
  write ( *, '(a)' ) '  representing unsigned integers;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '        I        J        K = I + J'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    
    i = i4_uniform ( 0, 100, seed )
    j = i4_uniform ( 0, 100, seed )

    write ( *, '(a)' ) ' '

    write ( *, '(2x,i8,2x,i8)' ) i, j

    k = i + j

    write ( *, '(a20,2x,i8)' ) '  Directly:         ', k

    call ui4_to_ubvec ( i, n, bvec1 )
    call ui4_to_ubvec ( j, n, bvec2 )

    call ubvec_add ( n, bvec1, bvec2, bvec3 )
    call ubvec_to_ui4 ( n, bvec3, k )

    write ( *, '(a20,2x,i8)' ) '  UBVEC_ADD         ', k

  end do

  return
end
subroutine test0626 ( )

!*****************************************************************************80
!
!! TEST0626 tests UI4_TO_UBVEC and UBVEC_TO_UI4;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 November 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  integer ( kind = 4 ) bvec(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0626'
  write ( *, '(a)' ) '  UI4_TO_UBVEC converts an unsigned integer to an '
  write ( *, '(a)' ) '  unsigned binary vector;'
  write ( *, '(a)' ) '  UBVEC_TO_UI4 converts an unsigned binary vector'
  write ( *, '(a)' ) '  to an unsigned integer;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I --> BVEC  -->  I'
  write ( *, '(a)' ) ' '
  do i = 0, 10
    call ui4_to_ubvec ( i, n, bvec )
    call ubvec_to_ui4 ( n, bvec, i2 )
    write ( *, '(2x,i3,2x,10i1,2x,i3)' ) i, bvec(1:n), i2
  end do

  return
end
subroutine test1535 ( )

!*****************************************************************************80
!
!! TEST1535 tests VEC_COLEX_NEXT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  integer ( kind = 4 ) a(dim_num)
  integer ( kind = 4 ) base
  logical more

  base = 3
  more = .false.

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1535'
  write ( *, '(a)' ) '  VEC_COLEX_NEXT generates all DIM_NUM-vectors'
  write ( *, '(a,i8)' ) '  in colex order in a given base BASE.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The dimension DIM_NUM = ', dim_num
  write ( *, '(a,i8)' ) '  The base BASE =         ', base
  write ( *, '(a)' ) ' '
 
  do

    call vec_colex_next ( dim_num, base, a, more )

    if ( .not. more ) then
      exit
    end if

    write ( *, '(2x,3i4)' ) a(1:dim_num)

  end do

  return
end
subroutine test1536 ( )

!*****************************************************************************80
!
!! TEST1536 tests VEC_COLEX_NEXT2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  integer ( kind = 4 ) a(dim_num)
  integer ( kind = 4 ), dimension(dim_num) :: base = (/ 2, 1, 3 /)
  logical more

  more = .false.

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1536'
  write ( *, '(a)' ) '  VEC_COLEX_NEXT2 generates all DIM_NUM-vectors'
  write ( *, '(a,i8)' ) '  in colex order in given bases BASE(1:DIM_NUM).'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The dimension DIM_NUM = ', dim_num
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The base vector:'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,3i4)' ) base(1:dim_num)
  write ( *, '(a)' ) ' '

  do

    call vec_colex_next2 ( dim_num, base, a, more )

    if ( .not. more ) then
      exit
    end if

    write ( *, '(2x,3i4)' ) a(1:dim_num)

  end do

  return
end
subroutine test1537 ( )

!*****************************************************************************80
!
!! TEST1537 tests VEC_COLEX_NEXT3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 August 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  integer ( kind = 4 ) a(dim_num)
  integer ( kind = 4 ), dimension(dim_num) :: base = (/ 2, 1, 3 /)
  logical more

  more = .false.

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1537'
  write ( *, '(a)' ) '  VEC_COLEX_NEXT3 generates all DIM_NUM-vectors'
  write ( *, '(a,i8)' ) '  in colex order in given bases BASE(1:DIM_NUM).'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The dimension DIM_NUM = ', dim_num
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The base vector:'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,3i4)' ) base(1:dim_num)
  write ( *, '(a)' ) ' '

  do

    call vec_colex_next3 ( dim_num, base, a, more )

    if ( .not. more ) then
      exit
    end if

    write ( *, '(2x,3i4)' ) a(1:dim_num)

  end do

  return
end
subroutine test155 ( )

!*****************************************************************************80
!
!! TEST155 tests VEC_GRAY_NEXT, VEC_GRAY_RANK and VEC_GRAY_UNRANK.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 4

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ), save, dimension ( n ) :: base = (/ 2, 2, 1, 4 /)
  integer ( kind = 4 ) change
  logical done
  integer ( kind = 4 ) prod
  integer ( kind = 4 ) rank

  prod = product ( base )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST155'
  write ( *, '(a)' ) '  VEC_GRAY_NEXT generates product space elements.'
  write ( *, '(a)' ) '  VEC_GRAY_RANK ranks them.'
  write ( *, '(a)' ) '  VEC_GRAY_UNRANK unranks them.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of components is ', n
  write ( *, '(a,i8)' ) '  The number of elements is ', prod
  write ( *, '(a)' ) '  Each component has its own number of degrees of'
  write ( *, '(a)' ) '  freedom, which, for this example, are:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,6i4)' ) '  Rank Change     ', base(1:n)
  write ( *, '(a)' ) ' '
  rank = 0
  done = .true.
 
  do
 
    rank = rank + 1
 
    call vec_gray_next ( n, base, a, done, change )
 
    if ( done ) then
      exit
    end if

    write ( *, '(2x,i4,2x,i4,2x,4x,6i4)' ) rank, change, a(1:n)

  end do
 
  a(1:n) = base(1:n) / 2

  call vec_gray_rank ( n, base, a, rank )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  VEC_GRAY_RANK reports the element '
  write ( *, '(a)' ) ' '
  write ( *, '(4x,3x,6i4)' ) a(1:n)
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  has rank ', rank

  rank = 7
  call vec_gray_unrank ( n, base, rank, a )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  VEC_GRAY_UNRANK reports the element of rank ', rank
  write ( *, '(a)' ) '  is:'
  write ( *, '(a)' ) ' '
  write ( *, '(4x,3x,6i4)' ) a(1:n)
  write ( *, '(a)' ) ' '

  return
end
subroutine test154 ( )

!*****************************************************************************80
!
!! TEST154 tests VEC_LEX_NEXT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) base
  logical more

  base = 3
  more = .false.
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST154'
  write ( *, '(a)' ) '  VEC_LEX_NEXT generates all N-vectors'
  write ( *, '(a,i8)' ) '  in a given base.  Here we use base ', base
  write ( *, '(a)' ) ' '
 
  do

    call vec_lex_next ( n, base, a, more )

    if ( .not. more ) then
      exit
    end if

    write ( *, '(2x,3i4)' ) a(1:n)

  end do

  return
end
subroutine test156 ( )

!*****************************************************************************80
!
!! TEST156 tests VEC_RANDOM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) base
  integer ( kind = 4 ) i
  integer ( kind = 4 ) seed

  base = 3
  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST156'
  write ( *, '(a)' ) '  VEC_RANDOM generates a random N-vector'
  write ( *, '(a)' ) '  in a given base.'
  write ( *, '(a,i8)' ) '  Here, we use base ', base
  write ( *, '(a)' ) ' '

  do i = 1, 5
    call vec_random ( n, base, seed, a )
    write ( *, '(2x,3i4)' ) a(1:n)
  end do
 
  return
end
subroutine test1565 ( )

!*****************************************************************************80
!
!! TEST1565 tests VECTOR_CONSTRAINED_NEXT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 March 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3

  integer ( kind = 4 ) constraint
  integer ( kind = 4 ) i
  logical more
  integer ( kind = 4 ) x(n)
  integer ( kind = 4 ) x_max(n)
  integer ( kind = 4 ) x_min(n)
  integer ( kind = 4 ) x_prod

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1565'
  write ( *, '(a)' ) '  VECTOR_CONSTRAINED_NEXT:'
  write ( *, '(a)' ) '  Consider vectors:'
  write ( *, '(a)' ) '    X_MIN(1:N) <= X(1:N) <= X_MAX(1:N),'
  write ( *, '(a)' ) '  Set'
  write ( *, '(a)' ) '    P = Product X_MAX(1:N)'
  write ( *, '(a)' ) '  Accept only vectors for which:'
  write ( *, '(a)' ) '    sum ( (X(1:N)-1) * P / X_MAX(1:N) ) <= P'

  more = .false.
  x_min(1:n) = (/ 2, 2, 1 /)
  x_max(1:n) = (/ 4, 5, 3 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a,3i4)' ) '  X_MIN:', x_min(1:n)
  write ( *, '(a,3i4)' ) '  X_MAX:', x_max(1:n)

  i = 0
  x_prod = product ( x_max(1:n) )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  Maximum allowed CONSTRAINT = P = ', x_prod
  write ( *, '(a)' ) ' '

  do

    call vector_constrained_next ( n, x_min, x_max, x, constraint, more )

    if ( .not. more ) then
      exit
    end if

    i = i + 1
    write ( *, '(2x,i8,2x,i12,2x,i8,2x,i8,2x,i8)' ) i, constraint, x(1:n)

  end do

  return
end
subroutine test1566 ( )

!*****************************************************************************80
!
!! TEST1566 tests VECTOR_CONSTRAINED_NEXT2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 March 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 3

  integer ( kind = 4 ) constraint
  integer ( kind = 4 ) i
  logical more
  integer ( kind = 4 ) n
  integer ( kind = 4 ) x(n_max)
  integer ( kind = 4 ) x_max(n_max)
  integer ( kind = 4 ) x_min(n_max)
  integer ( kind = 4 ) x_prod

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1566'
  write ( *, '(a)' ) '  VECTOR_CONSTRAINED_NEXT2:'
  write ( *, '(a)' ) '  Consider vectors:'
  write ( *, '(a)' ) '    X_MIN(1:N) <= X(1:N) <= X_MAX(1:N),'
  write ( *, '(a)' ) '  Set'
  write ( *, '(a)' ) '    P = Product X_MAX(1:N)'
  write ( *, '(a)' ) '  Accept only vectors for which:'
  write ( *, '(a)' ) '    sum ( X(1:N) * P / X_MAX(1:N) ) <= P'

  x_min(1:n_max) = (/ 1, 1, 1 /)
  x_max(1:n_max) = (/ 5, 6, 4 /)

  do n = 2, n_max

    more = .false.

    write ( *, '(a)' ) ' '
    write ( *, '(a,3i4)' ) '  X_MIN:', x_min(1:n)
    write ( *, '(a,3i4)' ) '  X_MAX:', x_max(1:n)

    i = 0
    x_prod = product ( x_max(1:n) )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i12)' ) '  Maximum allowed CONSTRAINT = P = ', x_prod
    write ( *, '(a)' ) ' '

    do

      call vector_constrained_next2 ( n, x_min, x_max, x, constraint, more )

      if ( .not. more ) then
        exit
      end if

      i = i + 1
      write ( *, '(2x,i8,2x,i12,2x,i8,2x,i8,2x,i8)' ) i, constraint, x(1:n)

    end do

  end do

  return
end
subroutine test1567 ( )

!*****************************************************************************80
!
!! TEST1567 tests VECTOR_CONSTRAINED_NEXT3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 April 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 3

  real ( kind = 8 ) constraint
  integer ( kind = 4 ) i
  logical more
  integer ( kind = 4 ) n
  integer ( kind = 4 ) x(n_max)
  integer ( kind = 4 ) x_max(n_max)
  integer ( kind = 4 ) x_min(n_max)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1567'
  write ( *, '(a)' ) '  VECTOR_CONSTRAINED_NEXT3:'
  write ( *, '(a)' ) '  Consider vectors:'
  write ( *, '(a)' ) '    X_MIN(1:N) <= X(1:N) <= X_MAX(1:N),'
  write ( *, '(a)' ) '  Set'
  write ( *, '(a)' ) '    CONSTRAINT = sum ( X(1:N) / X_MAX(1:N) )'
  write ( *, '(a)' ) '  Accept only vectors for which:'
  write ( *, '(a)' ) '    CONSTRAINT <= 1'

  x_min(1:n_max) = (/ 1, 1, 1 /)
  x_max(1:n_max) = (/ 5, 6, 4 /)

  do n = 2, n_max

    more = .false.

    write ( *, '(a)' ) ' '
    write ( *, '(a,3i4)' ) '  X_MIN:', x_min(1:n)
    write ( *, '(a,3i4)' ) '  X_MAX:', x_max(1:n)
    write ( *, '(a)' ) ' '

    i = 0

    do

      call vector_constrained_next3 ( n, x_min, x_max, x, constraint, more )

      if ( .not. more ) then
        exit
      end if

      i = i + 1
      write ( *, '(2x,i8,2x,g14.6,2x,i8,2x,i8,2x,i8)' ) i, constraint, x(1:n)

    end do

  end do

  return
end
subroutine test1568 ( )

!*****************************************************************************80
!
!! TEST1568 tests VECTOR_CONSTRAINED_NEXT4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 May 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 3

  real ( kind = 8 ), dimension ( n_max ) :: alpha = (/ &
    4.0D+00, 3.0D+00, 5.0D+00 /)
  integer ( kind = 4 ) i
  logical more
  integer ( kind = 4 ) n
  real ( kind = 8 ) :: q = 20.0D+00
  real ( kind = 8 ) total
  integer ( kind = 4 ) x(n_max)
  integer ( kind = 4 ), dimension ( n_max ) :: x_max = (/ &
    2, 6, 4 /)
  integer ( kind = 4 ), dimension ( n_max ) :: x_min = (/ &
    1, 0, 1 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1568'
  write ( *, '(a)' ) '  VECTOR_CONSTRAINED_NEXT4:'
  write ( *, '(a)' ) '  Consider vectors:'
  write ( *, '(a)' ) '    X_MIN(1:N) <= X(1:N) <= X_MAX(1:N),'
  write ( *, '(a)' ) '  Set'
  write ( *, '(a)' ) '    TOTAL = sum ( ALPHA(1:N) * X(1:N) )'
  write ( *, '(a)' ) '  Accept only vectors for which:'
  write ( *, '(a)' ) '    TOTAL <= Q'

  do n = 2, n_max

    more = .false.

    write ( *, '(a)' ) ' '
    write ( *, '(a,3g14.6)' ) '  ALPHA:', alpha(1:n)
    write ( *, '(a, g14.6)' ) '  Q:    ', q
    write ( *, '(a,3i4)'    ) '  X_MIN:', x_min(1:n)
    write ( *, '(a,3i4)'    ) '  X_MAX:', x_max(1:n)
    write ( *, '(a)' ) ' '

    i = 0

    do

      call vector_constrained_next4 ( n, alpha, x_min, x_max, x, q, more )

      if ( .not. more ) then
        exit
      end if

      total = dot_product ( alpha(1:n), real ( x(1:n), kind = 8 ) )
      i = i + 1
      write ( *, '(2x,i8,2x,g14.6,2x,i8,2x,i8,2x,i8)' ) i, total, x(1:n)

    end do

  end do

  return
end
subroutine test1569 ( )

!*****************************************************************************80
!
!! TEST1569 tests VECTOR_CONSTRAINED_NEXT5.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 December 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3

  integer ( kind = 4 ) i
  logical more
  integer ( kind = 4 ) sum_max
  integer ( kind = 4 ) sum_min
  integer ( kind = 4 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1569'
  write ( *, '(a)' ) '  VECTOR_CONSTRAINED_NEXT5:'
  write ( *, '(a)' ) '  Generate integer vectors X such that:'
  write ( *, '(a)' ) '    SUM_MIN <= sum ( X(1:N) ) <= SUM_MAX,'
  write ( *, '(a)' ) '  We require every X(I) to be at least 1.'

  more = .false.
  sum_min = 5
  sum_max = 7

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  N =       ', n
  write ( *, '(a,i8)' ) '  SUM_MIN = ', sum_min
  write ( *, '(a,i8)' ) '  SUM_MAX = ', sum_max
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         #        X(1)      X(2)      X(3)'
  write ( *, '(a)' ) ' '

  i = 0

  do

    call vector_constrained_next5 ( n, x, sum_min, sum_max, more )

    if ( .not. more ) then
      exit
    end if

    i = i + 1
    write ( *, '(2x,i8,2x,i8,2x,i8,2x,i8)' ) i, x(1:n)

  end do

  return
end
subroutine test15695 ( )

!*****************************************************************************80
!
!! TEST15695 tests VECTOR_CONSTRAINED_NEXT6.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 3

  real ( kind = 8 ), dimension ( n_max ) :: alpha = (/ &
    4.0D+00, 3.0D+00, 5.0D+00 /)
  integer ( kind = 4 ) i
  logical more
  integer ( kind = 4 ) n
  real ( kind = 8 ) :: q_max = 20.0D+00
  real ( kind = 8 ) :: q_min = 16.0D+00
  real ( kind = 8 ) total
  integer ( kind = 4 ) x(n_max)
  integer ( kind = 4 ), dimension ( n_max ) :: x_max = (/ &
    2, 6, 4 /)
  integer ( kind = 4 ), dimension ( n_max ) :: x_min = (/ &
    1, 0, 1 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST15695'
  write ( *, '(a)' ) '  VECTOR_CONSTRAINED_NEXT6:'
  write ( *, '(a)' ) '  Consider vectors:'
  write ( *, '(a)' ) '    X_MIN(1:N) <= X(1:N) <= X_MAX(1:N),'
  write ( *, '(a)' ) '  Set'
  write ( *, '(a)' ) '    TOTAL = sum ( ALPHA(1:N) * X(1:N) )'
  write ( *, '(a)' ) '  Accept only vectors for which:'
  write ( *, '(a)' ) '    Q_MIN <= TOTAL <= Q_MAX'

  do n = 2, n_max

    more = .false.

    write ( *, '(a)' ) ' '
    write ( *, '(a,3g14.6)' ) '  ALPHA:', alpha(1:n)
    write ( *, '(a, g14.6)' ) '  Q_MIN:', q_min
    write ( *, '(a, g14.6)' ) '  Q_MAX:', q_max
    write ( *, '(a,3i4)'    ) '  X_MIN:', x_min(1:n)
    write ( *, '(a,3i4)'    ) '  X_MAX:', x_max(1:n)
    write ( *, '(a)' ) ' '

    i = 0

    do

      call vector_constrained_next6 ( n, alpha, x_min, x_max, x, q_min, &
        q_max, more )

      if ( .not. more ) then
        exit
      end if

      total = dot_product ( alpha(1:n), real ( x(1:n), kind = 8 ) )
      i = i + 1
      write ( *, '(2x,i8,2x,g14.6,2x,i8,2x,i8,2x,i8)' ) i, total, x(1:n)

    end do

  end do

  return
end
subroutine test15696 ( )

!*****************************************************************************80
!
!! TEST15696 tests VECTOR_CONSTRAINED_NEXT7.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 May 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 3

  real ( kind = 8 ), dimension ( n_max ) :: alpha = (/ &
    4.0D+00, 3.0D+00, 5.0D+00 /)
  integer ( kind = 4 ) i
  logical more
  integer ( kind = 4 ) n
  real ( kind = 8 ) :: q_max = 20.0D+00
  real ( kind = 8 ) :: q_min = 16.0D+00
  real ( kind = 8 ) total
  integer ( kind = 4 ) x(n_max)
  integer ( kind = 4 ), dimension ( n_max ) :: x_max = (/ &
    2, 6, 4 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST15696'
  write ( *, '(a)' ) '  VECTOR_CONSTRAINED_NEXT7:'
  write ( *, '(a)' ) '  Consider vectors:'
  write ( *, '(a)' ) '    0 <= X(1:N) <= X_MAX(1:N),'
  write ( *, '(a)' ) '  Set'
  write ( *, '(a)' ) '    TOTAL = sum ( ALPHA(1:N) * X(1:N) )'
  write ( *, '(a)' ) '  Accept only vectors for which:'
  write ( *, '(a)' ) '    Q_MIN <= TOTAL <= Q_MAX'

  do n = 2, n_max

    more = .false.

    write ( *, '(a)' ) ' '
    write ( *, '(a,3g14.6)' ) '  ALPHA:', alpha(1:n)
    write ( *, '(a, g14.6)' ) '  Q_MIN:', q_min
    write ( *, '(a, g14.6)' ) '  Q_MAX:', q_max
    write ( *, '(a,3i4)'    ) '  X_MAX:', x_max(1:n)
    write ( *, '(a)' ) ' '

    i = 0

    do

      call vector_constrained_next7 ( n, alpha, x_max, x, q_min, &
        q_max, more )

      if ( .not. more ) then
        exit
      end if

      total = dot_product ( alpha(1:n), real ( x(1:n), kind = 8 ) )
      i = i + 1
      write ( *, '(2x,i8,2x,g14.6,2x,i8,2x,i8,2x,i8)' ) i, total, x(1:n)

    end do

  end do

  return
end
subroutine test15698 ( )

!*****************************************************************************80
!
!! TEST15698 tests VECTOR_NEXT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 July 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 3

  integer ( kind = 4 ) i
  logical              more
  integer ( kind = 4 ) n
  integer ( kind = 4 ) x(n_max)
  integer ( kind = 4 ), dimension ( n_max ) :: x_max = (/ &
    2, 6, 4 /)
  integer ( kind = 4 ), dimension ( n_max ) :: x_min = (/ &
    1, 4, 3 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST15698'
  write ( *, '(a)' ) '  VECTOR_NEXT:'
  write ( *, '(a)' ) '  Generate all vectors X such that:'
  write ( *, '(a)' ) '    X_MIN(1:N) <= X(1:N) <= X_MAX(1:N),'

  do n = 2, n_max

    more = .false.

    write ( *, '(a)' ) ' '
    write ( *, '(2x,a4,4x,i8,2x,i8,2x,i8)' ) 'XMIN', x_min(1:n)

    i = 0

    do

      call vector_next ( n, x_min, x_max, x, more )

      if ( .not. more ) then
        exit
      end if

      i = i + 1
      write ( *, '(2x,i4,4x,i8,2x,i8,2x,i8)' ) i, x(1:n)

    end do

    write ( *, '(2x,a4,4x,i8,2x,i8,2x,i8)' ) 'XMAX', x_max(1:n)

  end do

  return
end
subroutine test157 ( )

!*****************************************************************************80
!
!! TEST157 tests YTB_ENUM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10

  integer ( kind = 4 ) i
  integer ( kind = 4 ) pi

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST157'
  write ( *, '(a)' ) '  YTB_ENUM counts Young table.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  N    YTB(N)'
  write ( *, '(a)' ) ' '

  do i = 0, n
    call ytb_enum ( i, pi )
    write ( *, '(2x,2i10)' ) i, pi
  end do

  return
end
subroutine test158 ( )

!*****************************************************************************80
!
!! TEST158 tests YTB_NEXT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 6

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ), dimension ( n ) :: lambda = (/ 3, 2, 1, 0, 0, 0 /)
  logical more

  a(1:n) = 0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST158'
  write ( *, '(a)' ) '  YTB_NEXT generates Young tables.'
  write ( *, '(a)' ) ' '
  more = .false.
 
  do
 
    call ytb_next ( n, lambda, a, more )
 
    call ytb_print ( n, a, ' ' )
 
    if ( .not. more ) then
      exit
    end if

  end do

  return
end
subroutine test159 ( )

!*****************************************************************************80
!
!! TEST159 tests YTB_RANDOM.
!
!  Discussion:
!
!    I4_LOG_10 ( I ) + 1 is the number of decimal digits in I.
!
!  Modified:
!
!    05 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 6

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ), dimension ( n ) :: lambda = (/ 3, 2, 1, 0, 0, 0 /)
  integer ( kind = 4 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST159'
  write ( *, '(a)' ) '  YTB_RANDOM generates a random Young table'

  seed = 123456789

  do i = 1, 5
 
    call ytb_random ( n, lambda, seed, a )

    call ytb_print ( n, a, ' ' )
 
  end do
 
  return
end
